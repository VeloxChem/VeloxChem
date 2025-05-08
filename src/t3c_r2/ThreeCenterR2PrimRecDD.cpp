#include "ThreeCenterR2PrimRecDD.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_dd(CSimdArray<double>& pbuffer, 
                const size_t idx_g_dd,
                const size_t idx_sd,
                const size_t idx_pp,
                const size_t idx_pd,
                const size_t idx_ds,
                const size_t idx_dp,
                const size_t idx_dd,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_xy = pbuffer.data(idx_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

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

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    auto ts_xy_x = pbuffer.data(idx_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_dp + 5);

    auto ts_xz_x = pbuffer.data(idx_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    auto ts_yz_x = pbuffer.data(idx_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_dd);

    auto ts_xx_xy = pbuffer.data(idx_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_dd + 5);

    auto ts_xy_xx = pbuffer.data(idx_dd + 6);

    auto ts_xy_xy = pbuffer.data(idx_dd + 7);

    auto ts_xy_xz = pbuffer.data(idx_dd + 8);

    auto ts_xy_yy = pbuffer.data(idx_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_dd + 10);

    auto ts_xy_zz = pbuffer.data(idx_dd + 11);

    auto ts_xz_xx = pbuffer.data(idx_dd + 12);

    auto ts_xz_xy = pbuffer.data(idx_dd + 13);

    auto ts_xz_xz = pbuffer.data(idx_dd + 14);

    auto ts_xz_yy = pbuffer.data(idx_dd + 15);

    auto ts_xz_yz = pbuffer.data(idx_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_dd + 23);

    auto ts_yz_xx = pbuffer.data(idx_dd + 24);

    auto ts_yz_xy = pbuffer.data(idx_dd + 25);

    auto ts_yz_xz = pbuffer.data(idx_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_dd + 35);

    // Set up 0-6 components of targeted buffer : DD

    auto gr_xx_xx = pbuffer.data(idx_g_dd);

    auto gr_xx_xy = pbuffer.data(idx_g_dd + 1);

    auto gr_xx_xz = pbuffer.data(idx_g_dd + 2);

    auto gr_xx_yy = pbuffer.data(idx_g_dd + 3);

    auto gr_xx_yz = pbuffer.data(idx_g_dd + 4);

    auto gr_xx_zz = pbuffer.data(idx_g_dd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_xx, gr_xx_xy, gr_xx_xz, gr_xx_yy, gr_xx_yz, gr_xx_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, ts_xx_0, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xx_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 + 8.0 * ts_x_x[i] * gfe_0 + 4.0 * ts_x_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 + 4.0 * ts_xx_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xx_xx[i] * gfe_0 + ts_xx_xx[i] * rgc2_0;

        gr_xx_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 + 4.0 * ts_x_y[i] * gfe_0 + 4.0 * ts_x_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_xy[i] * gfe_0 + ts_xx_xy[i] * rgc2_0;

        gr_xx_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 + 4.0 * ts_x_z[i] * gfe_0 + 4.0 * ts_x_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_xz[i] * gfe_0 + ts_xx_xz[i] * rgc2_0;

        gr_xx_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 + 4.0 * ts_x_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 + 4.0 * ts_xx_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_yy[i] * gfe_0 + ts_xx_yy[i] * rgc2_0;

        gr_xx_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 + 4.0 * ts_x_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xx_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_yz[i] * gfe_0 + ts_xx_yz[i] * rgc2_0;

        gr_xx_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 + 4.0 * ts_x_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 + 4.0 * ts_xx_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_zz[i] * gfe_0 + ts_xx_zz[i] * rgc2_0;
    }

    // Set up 6-12 components of targeted buffer : DD

    auto gr_xy_xx = pbuffer.data(idx_g_dd + 6);

    auto gr_xy_xy = pbuffer.data(idx_g_dd + 7);

    auto gr_xy_xz = pbuffer.data(idx_g_dd + 8);

    auto gr_xy_yy = pbuffer.data(idx_g_dd + 9);

    auto gr_xy_yz = pbuffer.data(idx_g_dd + 10);

    auto gr_xy_zz = pbuffer.data(idx_g_dd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_xx, gr_xy_xy, gr_xy_xz, gr_xy_yy, gr_xy_yz, gr_xy_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, ts_xy_0, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xy_xx[i] = 4.0 * ts_y_x[i] * gfe_0 + 2.0 * ts_y_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 + 4.0 * ts_xy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xy_xx[i] * gfe_0 + ts_xy_xx[i] * rgc2_0;

        gr_xy_xy[i] = 2.0 * ts_y_y[i] * gfe_0 + 2.0 * ts_y_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 + 2.0 * ts_x_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_xy[i] * gfe_0 + ts_xy_xy[i] * rgc2_0;

        gr_xy_xz[i] = 2.0 * ts_y_z[i] * gfe_0 + 2.0 * ts_y_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_xz[i] * gfe_0 + ts_xy_xz[i] * rgc2_0;

        gr_xy_yy[i] = 2.0 * ts_y_yy[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_y[i] * gfe_0 + 2.0 * ts_x_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 + 4.0 * ts_xy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_yy[i] * gfe_0 + ts_xy_yy[i] * rgc2_0;

        gr_xy_yz[i] = 2.0 * ts_y_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 + 2.0 * ts_x_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_yz[i] * gfe_0 + ts_xy_yz[i] * rgc2_0;

        gr_xy_zz[i] = 2.0 * ts_y_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 + 4.0 * ts_xy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_zz[i] * gfe_0 + ts_xy_zz[i] * rgc2_0;
    }

    // Set up 12-18 components of targeted buffer : DD

    auto gr_xz_xx = pbuffer.data(idx_g_dd + 12);

    auto gr_xz_xy = pbuffer.data(idx_g_dd + 13);

    auto gr_xz_xz = pbuffer.data(idx_g_dd + 14);

    auto gr_xz_yy = pbuffer.data(idx_g_dd + 15);

    auto gr_xz_yz = pbuffer.data(idx_g_dd + 16);

    auto gr_xz_zz = pbuffer.data(idx_g_dd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_xx, gr_xz_xy, gr_xz_xz, gr_xz_yy, gr_xz_yz, gr_xz_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, ts_xz_0, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xz_xx[i] = 4.0 * ts_z_x[i] * gfe_0 + 2.0 * ts_z_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_0[i] * gfe_0 + 4.0 * ts_xz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xz_xx[i] * gfe_0 + ts_xz_xx[i] * rgc2_0;

        gr_xz_xy[i] = 2.0 * ts_z_y[i] * gfe_0 + 2.0 * ts_z_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_xy[i] * gfe_0 + ts_xz_xy[i] * rgc2_0;

        gr_xz_xz[i] = 2.0 * ts_z_z[i] * gfe_0 + 2.0 * ts_z_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 + 2.0 * ts_x_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_xz[i] * gfe_0 + ts_xz_xz[i] * rgc2_0;

        gr_xz_yy[i] = 2.0 * ts_z_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_0[i] * gfe_0 + 4.0 * ts_xz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_yy[i] * gfe_0 + ts_xz_yy[i] * rgc2_0;

        gr_xz_yz[i] = 2.0 * ts_z_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_y[i] * gfe_0 + 2.0 * ts_x_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_yz[i] * gfe_0 + ts_xz_yz[i] * rgc2_0;

        gr_xz_zz[i] = 2.0 * ts_z_zz[i] * gfe_0 * gc_x[i] + 4.0 * ts_x_z[i] * gfe_0 + 2.0 * ts_x_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_0[i] * gfe_0 + 4.0 * ts_xz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_zz[i] * gfe_0 + ts_xz_zz[i] * rgc2_0;
    }

    // Set up 18-24 components of targeted buffer : DD

    auto gr_yy_xx = pbuffer.data(idx_g_dd + 18);

    auto gr_yy_xy = pbuffer.data(idx_g_dd + 19);

    auto gr_yy_xz = pbuffer.data(idx_g_dd + 20);

    auto gr_yy_yy = pbuffer.data(idx_g_dd + 21);

    auto gr_yy_yz = pbuffer.data(idx_g_dd + 22);

    auto gr_yy_zz = pbuffer.data(idx_g_dd + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_xx, gr_yy_xy, gr_yy_xz, gr_yy_yy, gr_yy_yz, gr_yy_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, ts_yy_0, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yy_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 + 4.0 * ts_y_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 + 4.0 * ts_yy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yy_xx[i] * gfe_0 + ts_yy_xx[i] * rgc2_0;

        gr_yy_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 + 4.0 * ts_y_x[i] * gfe_0 + 4.0 * ts_y_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_xy[i] * gfe_0 + ts_yy_xy[i] * rgc2_0;

        gr_yy_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 + 4.0 * ts_y_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_xz[i] * gfe_0 + ts_yy_xz[i] * rgc2_0;

        gr_yy_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 + 8.0 * ts_y_y[i] * gfe_0 + 4.0 * ts_y_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 + 4.0 * ts_yy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_yy[i] * gfe_0 + ts_yy_yy[i] * rgc2_0;

        gr_yy_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 + 4.0 * ts_y_z[i] * gfe_0 + 4.0 * ts_y_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_yz[i] * gfe_0 + ts_yy_yz[i] * rgc2_0;

        gr_yy_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 + 4.0 * ts_y_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 + 4.0 * ts_yy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_zz[i] * gfe_0 + ts_yy_zz[i] * rgc2_0;
    }

    // Set up 24-30 components of targeted buffer : DD

    auto gr_yz_xx = pbuffer.data(idx_g_dd + 24);

    auto gr_yz_xy = pbuffer.data(idx_g_dd + 25);

    auto gr_yz_xz = pbuffer.data(idx_g_dd + 26);

    auto gr_yz_yy = pbuffer.data(idx_g_dd + 27);

    auto gr_yz_yz = pbuffer.data(idx_g_dd + 28);

    auto gr_yz_zz = pbuffer.data(idx_g_dd + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_xx, gr_yz_xy, gr_yz_xz, gr_yz_yy, gr_yz_yz, gr_yz_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, ts_yz_0, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yz_xx[i] = 2.0 * ts_z_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_0[i] * gfe_0 + 4.0 * ts_yz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yz_xx[i] * gfe_0 + ts_yz_xx[i] * rgc2_0;

        gr_yz_xy[i] = 2.0 * ts_z_x[i] * gfe_0 + 2.0 * ts_z_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_xy[i] * gfe_0 + ts_yz_xy[i] * rgc2_0;

        gr_yz_xz[i] = 2.0 * ts_z_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_x[i] * gfe_0 + 2.0 * ts_y_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_xz[i] * gfe_0 + ts_yz_xz[i] * rgc2_0;

        gr_yz_yy[i] = 4.0 * ts_z_y[i] * gfe_0 + 2.0 * ts_z_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_0[i] * gfe_0 + 4.0 * ts_yz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_yy[i] * gfe_0 + ts_yz_yy[i] * rgc2_0;

        gr_yz_yz[i] = 2.0 * ts_z_z[i] * gfe_0 + 2.0 * ts_z_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 + 2.0 * ts_y_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_yz[i] * gfe_0 + ts_yz_yz[i] * rgc2_0;

        gr_yz_zz[i] = 2.0 * ts_z_zz[i] * gfe_0 * gc_y[i] + 4.0 * ts_y_z[i] * gfe_0 + 2.0 * ts_y_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_0[i] * gfe_0 + 4.0 * ts_yz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_zz[i] * gfe_0 + ts_yz_zz[i] * rgc2_0;
    }

    // Set up 30-36 components of targeted buffer : DD

    auto gr_zz_xx = pbuffer.data(idx_g_dd + 30);

    auto gr_zz_xy = pbuffer.data(idx_g_dd + 31);

    auto gr_zz_xz = pbuffer.data(idx_g_dd + 32);

    auto gr_zz_yy = pbuffer.data(idx_g_dd + 33);

    auto gr_zz_yz = pbuffer.data(idx_g_dd + 34);

    auto gr_zz_zz = pbuffer.data(idx_g_dd + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_xx, gr_zz_xy, gr_zz_xz, gr_zz_yy, gr_zz_yz, gr_zz_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, ts_zz_0, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zz_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 + 4.0 * ts_z_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_0[i] * gfe_0 + 4.0 * ts_zz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_zz_xx[i] * gfe_0 + ts_zz_xx[i] * rgc2_0;

        gr_zz_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 + 4.0 * ts_z_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_xy[i] * gfe_0 + ts_zz_xy[i] * rgc2_0;

        gr_zz_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 + 4.0 * ts_z_x[i] * gfe_0 + 4.0 * ts_z_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_zz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_xz[i] * gfe_0 + ts_zz_xz[i] * rgc2_0;

        gr_zz_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 + 4.0 * ts_z_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_0[i] * gfe_0 + 4.0 * ts_zz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_yy[i] * gfe_0 + ts_zz_yy[i] * rgc2_0;

        gr_zz_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 + 4.0 * ts_z_y[i] * gfe_0 + 4.0 * ts_z_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_zz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_yz[i] * gfe_0 + ts_zz_yz[i] * rgc2_0;

        gr_zz_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 + 8.0 * ts_z_z[i] * gfe_0 + 4.0 * ts_z_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_0[i] * gfe_0 + 4.0 * ts_zz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_zz[i] * gfe_0 + ts_zz_zz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

