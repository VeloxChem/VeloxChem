#include "ThreeCenterRR2PrimRecPD.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_pd(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_pd,
                  const size_t idx_sd,
                  const size_t idx_g_sd,
                  const size_t idx_pp,
                  const size_t idx_g_pp,
                  const size_t idx_pd,
                  const size_t idx_g_pd,
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

    // Set up components of auxiliary buffer : SD

    auto gr_0_xx = pbuffer.data(idx_g_sd);

    auto gr_0_xy = pbuffer.data(idx_g_sd + 1);

    auto gr_0_xz = pbuffer.data(idx_g_sd + 2);

    auto gr_0_yy = pbuffer.data(idx_g_sd + 3);

    auto gr_0_yz = pbuffer.data(idx_g_sd + 4);

    auto gr_0_zz = pbuffer.data(idx_g_sd + 5);

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

    // Set up components of auxiliary buffer : PP

    auto gr_x_x = pbuffer.data(idx_g_pp);

    auto gr_x_y = pbuffer.data(idx_g_pp + 1);

    auto gr_x_z = pbuffer.data(idx_g_pp + 2);

    auto gr_y_x = pbuffer.data(idx_g_pp + 3);

    auto gr_y_y = pbuffer.data(idx_g_pp + 4);

    auto gr_y_z = pbuffer.data(idx_g_pp + 5);

    auto gr_z_x = pbuffer.data(idx_g_pp + 6);

    auto gr_z_y = pbuffer.data(idx_g_pp + 7);

    auto gr_z_z = pbuffer.data(idx_g_pp + 8);

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

    // Set up components of auxiliary buffer : PD

    auto gr_x_xx = pbuffer.data(idx_g_pd);

    auto gr_x_xy = pbuffer.data(idx_g_pd + 1);

    auto gr_x_xz = pbuffer.data(idx_g_pd + 2);

    auto gr_x_yy = pbuffer.data(idx_g_pd + 3);

    auto gr_x_yz = pbuffer.data(idx_g_pd + 4);

    auto gr_x_zz = pbuffer.data(idx_g_pd + 5);

    auto gr_y_xx = pbuffer.data(idx_g_pd + 6);

    auto gr_y_xy = pbuffer.data(idx_g_pd + 7);

    auto gr_y_xz = pbuffer.data(idx_g_pd + 8);

    auto gr_y_yy = pbuffer.data(idx_g_pd + 9);

    auto gr_y_yz = pbuffer.data(idx_g_pd + 10);

    auto gr_y_zz = pbuffer.data(idx_g_pd + 11);

    auto gr_z_xx = pbuffer.data(idx_g_pd + 12);

    auto gr_z_xy = pbuffer.data(idx_g_pd + 13);

    auto gr_z_xz = pbuffer.data(idx_g_pd + 14);

    auto gr_z_yy = pbuffer.data(idx_g_pd + 15);

    auto gr_z_yz = pbuffer.data(idx_g_pd + 16);

    auto gr_z_zz = pbuffer.data(idx_g_pd + 17);

    // Set up 0-6 components of targeted buffer : PD

    auto grr_x_x_xx = pbuffer.data(idx_gr_pd);

    auto grr_x_x_xy = pbuffer.data(idx_gr_pd + 1);

    auto grr_x_x_xz = pbuffer.data(idx_gr_pd + 2);

    auto grr_x_x_yy = pbuffer.data(idx_gr_pd + 3);

    auto grr_x_x_yz = pbuffer.data(idx_gr_pd + 4);

    auto grr_x_x_zz = pbuffer.data(idx_gr_pd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xx, gr_0_xy, gr_0_xz, gr_0_yy, gr_0_yz, gr_0_zz, gr_x_x, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_y, gr_x_yy, gr_x_yz, gr_x_z, gr_x_zz, grr_x_x_xx, grr_x_x_xy, grr_x_x_xz, grr_x_x_yy, grr_x_x_yz, grr_x_x_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_x_xx[i] = ts_0_xx[i] * gfe2_0 + gr_0_xx[i] * gfe_0 + 2.0 * ts_x_x[i] * gfe2_0 + 2.0 * gr_x_x[i] * gfe_0 + ts_x_xx[i] * gfe_0 * gc_x[i] + gr_x_xx[i] * gc_x[i];

        grr_x_x_xy[i] = ts_0_xy[i] * gfe2_0 + gr_0_xy[i] * gfe_0 + ts_x_y[i] * gfe2_0 + gr_x_y[i] * gfe_0 + ts_x_xy[i] * gfe_0 * gc_x[i] + gr_x_xy[i] * gc_x[i];

        grr_x_x_xz[i] = ts_0_xz[i] * gfe2_0 + gr_0_xz[i] * gfe_0 + ts_x_z[i] * gfe2_0 + gr_x_z[i] * gfe_0 + ts_x_xz[i] * gfe_0 * gc_x[i] + gr_x_xz[i] * gc_x[i];

        grr_x_x_yy[i] = ts_0_yy[i] * gfe2_0 + gr_0_yy[i] * gfe_0 + ts_x_yy[i] * gfe_0 * gc_x[i] + gr_x_yy[i] * gc_x[i];

        grr_x_x_yz[i] = ts_0_yz[i] * gfe2_0 + gr_0_yz[i] * gfe_0 + ts_x_yz[i] * gfe_0 * gc_x[i] + gr_x_yz[i] * gc_x[i];

        grr_x_x_zz[i] = ts_0_zz[i] * gfe2_0 + gr_0_zz[i] * gfe_0 + ts_x_zz[i] * gfe_0 * gc_x[i] + gr_x_zz[i] * gc_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto grr_x_y_xx = pbuffer.data(idx_gr_pd + 6);

    auto grr_x_y_xy = pbuffer.data(idx_gr_pd + 7);

    auto grr_x_y_xz = pbuffer.data(idx_gr_pd + 8);

    auto grr_x_y_yy = pbuffer.data(idx_gr_pd + 9);

    auto grr_x_y_yz = pbuffer.data(idx_gr_pd + 10);

    auto grr_x_y_zz = pbuffer.data(idx_gr_pd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_x, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_y, gr_y_yy, gr_y_yz, gr_y_z, gr_y_zz, grr_x_y_xx, grr_x_y_xy, grr_x_y_xz, grr_x_y_yy, grr_x_y_yz, grr_x_y_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_y_xx[i] = 2.0 * ts_y_x[i] * gfe2_0 + 2.0 * gr_y_x[i] * gfe_0 + ts_y_xx[i] * gfe_0 * gc_x[i] + gr_y_xx[i] * gc_x[i];

        grr_x_y_xy[i] = ts_y_y[i] * gfe2_0 + gr_y_y[i] * gfe_0 + ts_y_xy[i] * gfe_0 * gc_x[i] + gr_y_xy[i] * gc_x[i];

        grr_x_y_xz[i] = ts_y_z[i] * gfe2_0 + gr_y_z[i] * gfe_0 + ts_y_xz[i] * gfe_0 * gc_x[i] + gr_y_xz[i] * gc_x[i];

        grr_x_y_yy[i] = ts_y_yy[i] * gfe_0 * gc_x[i] + gr_y_yy[i] * gc_x[i];

        grr_x_y_yz[i] = ts_y_yz[i] * gfe_0 * gc_x[i] + gr_y_yz[i] * gc_x[i];

        grr_x_y_zz[i] = ts_y_zz[i] * gfe_0 * gc_x[i] + gr_y_zz[i] * gc_x[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto grr_x_z_xx = pbuffer.data(idx_gr_pd + 12);

    auto grr_x_z_xy = pbuffer.data(idx_gr_pd + 13);

    auto grr_x_z_xz = pbuffer.data(idx_gr_pd + 14);

    auto grr_x_z_yy = pbuffer.data(idx_gr_pd + 15);

    auto grr_x_z_yz = pbuffer.data(idx_gr_pd + 16);

    auto grr_x_z_zz = pbuffer.data(idx_gr_pd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_x, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_y, gr_z_yy, gr_z_yz, gr_z_z, gr_z_zz, grr_x_z_xx, grr_x_z_xy, grr_x_z_xz, grr_x_z_yy, grr_x_z_yz, grr_x_z_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_z_xx[i] = 2.0 * ts_z_x[i] * gfe2_0 + 2.0 * gr_z_x[i] * gfe_0 + ts_z_xx[i] * gfe_0 * gc_x[i] + gr_z_xx[i] * gc_x[i];

        grr_x_z_xy[i] = ts_z_y[i] * gfe2_0 + gr_z_y[i] * gfe_0 + ts_z_xy[i] * gfe_0 * gc_x[i] + gr_z_xy[i] * gc_x[i];

        grr_x_z_xz[i] = ts_z_z[i] * gfe2_0 + gr_z_z[i] * gfe_0 + ts_z_xz[i] * gfe_0 * gc_x[i] + gr_z_xz[i] * gc_x[i];

        grr_x_z_yy[i] = ts_z_yy[i] * gfe_0 * gc_x[i] + gr_z_yy[i] * gc_x[i];

        grr_x_z_yz[i] = ts_z_yz[i] * gfe_0 * gc_x[i] + gr_z_yz[i] * gc_x[i];

        grr_x_z_zz[i] = ts_z_zz[i] * gfe_0 * gc_x[i] + gr_z_zz[i] * gc_x[i];
    }

    // Set up 18-24 components of targeted buffer : PD

    auto grr_y_x_xx = pbuffer.data(idx_gr_pd + 18);

    auto grr_y_x_xy = pbuffer.data(idx_gr_pd + 19);

    auto grr_y_x_xz = pbuffer.data(idx_gr_pd + 20);

    auto grr_y_x_yy = pbuffer.data(idx_gr_pd + 21);

    auto grr_y_x_yz = pbuffer.data(idx_gr_pd + 22);

    auto grr_y_x_zz = pbuffer.data(idx_gr_pd + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_x, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_y, gr_x_yy, gr_x_yz, gr_x_z, gr_x_zz, grr_y_x_xx, grr_y_x_xy, grr_y_x_xz, grr_y_x_yy, grr_y_x_yz, grr_y_x_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_x_xx[i] = ts_x_xx[i] * gfe_0 * gc_y[i] + gr_x_xx[i] * gc_y[i];

        grr_y_x_xy[i] = ts_x_x[i] * gfe2_0 + gr_x_x[i] * gfe_0 + ts_x_xy[i] * gfe_0 * gc_y[i] + gr_x_xy[i] * gc_y[i];

        grr_y_x_xz[i] = ts_x_xz[i] * gfe_0 * gc_y[i] + gr_x_xz[i] * gc_y[i];

        grr_y_x_yy[i] = 2.0 * ts_x_y[i] * gfe2_0 + 2.0 * gr_x_y[i] * gfe_0 + ts_x_yy[i] * gfe_0 * gc_y[i] + gr_x_yy[i] * gc_y[i];

        grr_y_x_yz[i] = ts_x_z[i] * gfe2_0 + gr_x_z[i] * gfe_0 + ts_x_yz[i] * gfe_0 * gc_y[i] + gr_x_yz[i] * gc_y[i];

        grr_y_x_zz[i] = ts_x_zz[i] * gfe_0 * gc_y[i] + gr_x_zz[i] * gc_y[i];
    }

    // Set up 24-30 components of targeted buffer : PD

    auto grr_y_y_xx = pbuffer.data(idx_gr_pd + 24);

    auto grr_y_y_xy = pbuffer.data(idx_gr_pd + 25);

    auto grr_y_y_xz = pbuffer.data(idx_gr_pd + 26);

    auto grr_y_y_yy = pbuffer.data(idx_gr_pd + 27);

    auto grr_y_y_yz = pbuffer.data(idx_gr_pd + 28);

    auto grr_y_y_zz = pbuffer.data(idx_gr_pd + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xx, gr_0_xy, gr_0_xz, gr_0_yy, gr_0_yz, gr_0_zz, gr_y_x, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_y, gr_y_yy, gr_y_yz, gr_y_z, gr_y_zz, grr_y_y_xx, grr_y_y_xy, grr_y_y_xz, grr_y_y_yy, grr_y_y_yz, grr_y_y_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_y_xx[i] = ts_0_xx[i] * gfe2_0 + gr_0_xx[i] * gfe_0 + ts_y_xx[i] * gfe_0 * gc_y[i] + gr_y_xx[i] * gc_y[i];

        grr_y_y_xy[i] = ts_0_xy[i] * gfe2_0 + gr_0_xy[i] * gfe_0 + ts_y_x[i] * gfe2_0 + gr_y_x[i] * gfe_0 + ts_y_xy[i] * gfe_0 * gc_y[i] + gr_y_xy[i] * gc_y[i];

        grr_y_y_xz[i] = ts_0_xz[i] * gfe2_0 + gr_0_xz[i] * gfe_0 + ts_y_xz[i] * gfe_0 * gc_y[i] + gr_y_xz[i] * gc_y[i];

        grr_y_y_yy[i] = ts_0_yy[i] * gfe2_0 + gr_0_yy[i] * gfe_0 + 2.0 * ts_y_y[i] * gfe2_0 + 2.0 * gr_y_y[i] * gfe_0 + ts_y_yy[i] * gfe_0 * gc_y[i] + gr_y_yy[i] * gc_y[i];

        grr_y_y_yz[i] = ts_0_yz[i] * gfe2_0 + gr_0_yz[i] * gfe_0 + ts_y_z[i] * gfe2_0 + gr_y_z[i] * gfe_0 + ts_y_yz[i] * gfe_0 * gc_y[i] + gr_y_yz[i] * gc_y[i];

        grr_y_y_zz[i] = ts_0_zz[i] * gfe2_0 + gr_0_zz[i] * gfe_0 + ts_y_zz[i] * gfe_0 * gc_y[i] + gr_y_zz[i] * gc_y[i];
    }

    // Set up 30-36 components of targeted buffer : PD

    auto grr_y_z_xx = pbuffer.data(idx_gr_pd + 30);

    auto grr_y_z_xy = pbuffer.data(idx_gr_pd + 31);

    auto grr_y_z_xz = pbuffer.data(idx_gr_pd + 32);

    auto grr_y_z_yy = pbuffer.data(idx_gr_pd + 33);

    auto grr_y_z_yz = pbuffer.data(idx_gr_pd + 34);

    auto grr_y_z_zz = pbuffer.data(idx_gr_pd + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_x, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_y, gr_z_yy, gr_z_yz, gr_z_z, gr_z_zz, grr_y_z_xx, grr_y_z_xy, grr_y_z_xz, grr_y_z_yy, grr_y_z_yz, grr_y_z_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_z_xx[i] = ts_z_xx[i] * gfe_0 * gc_y[i] + gr_z_xx[i] * gc_y[i];

        grr_y_z_xy[i] = ts_z_x[i] * gfe2_0 + gr_z_x[i] * gfe_0 + ts_z_xy[i] * gfe_0 * gc_y[i] + gr_z_xy[i] * gc_y[i];

        grr_y_z_xz[i] = ts_z_xz[i] * gfe_0 * gc_y[i] + gr_z_xz[i] * gc_y[i];

        grr_y_z_yy[i] = 2.0 * ts_z_y[i] * gfe2_0 + 2.0 * gr_z_y[i] * gfe_0 + ts_z_yy[i] * gfe_0 * gc_y[i] + gr_z_yy[i] * gc_y[i];

        grr_y_z_yz[i] = ts_z_z[i] * gfe2_0 + gr_z_z[i] * gfe_0 + ts_z_yz[i] * gfe_0 * gc_y[i] + gr_z_yz[i] * gc_y[i];

        grr_y_z_zz[i] = ts_z_zz[i] * gfe_0 * gc_y[i] + gr_z_zz[i] * gc_y[i];
    }

    // Set up 36-42 components of targeted buffer : PD

    auto grr_z_x_xx = pbuffer.data(idx_gr_pd + 36);

    auto grr_z_x_xy = pbuffer.data(idx_gr_pd + 37);

    auto grr_z_x_xz = pbuffer.data(idx_gr_pd + 38);

    auto grr_z_x_yy = pbuffer.data(idx_gr_pd + 39);

    auto grr_z_x_yz = pbuffer.data(idx_gr_pd + 40);

    auto grr_z_x_zz = pbuffer.data(idx_gr_pd + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_x, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_y, gr_x_yy, gr_x_yz, gr_x_z, gr_x_zz, grr_z_x_xx, grr_z_x_xy, grr_z_x_xz, grr_z_x_yy, grr_z_x_yz, grr_z_x_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_x_xx[i] = ts_x_xx[i] * gfe_0 * gc_z[i] + gr_x_xx[i] * gc_z[i];

        grr_z_x_xy[i] = ts_x_xy[i] * gfe_0 * gc_z[i] + gr_x_xy[i] * gc_z[i];

        grr_z_x_xz[i] = ts_x_x[i] * gfe2_0 + gr_x_x[i] * gfe_0 + ts_x_xz[i] * gfe_0 * gc_z[i] + gr_x_xz[i] * gc_z[i];

        grr_z_x_yy[i] = ts_x_yy[i] * gfe_0 * gc_z[i] + gr_x_yy[i] * gc_z[i];

        grr_z_x_yz[i] = ts_x_y[i] * gfe2_0 + gr_x_y[i] * gfe_0 + ts_x_yz[i] * gfe_0 * gc_z[i] + gr_x_yz[i] * gc_z[i];

        grr_z_x_zz[i] = 2.0 * ts_x_z[i] * gfe2_0 + 2.0 * gr_x_z[i] * gfe_0 + ts_x_zz[i] * gfe_0 * gc_z[i] + gr_x_zz[i] * gc_z[i];
    }

    // Set up 42-48 components of targeted buffer : PD

    auto grr_z_y_xx = pbuffer.data(idx_gr_pd + 42);

    auto grr_z_y_xy = pbuffer.data(idx_gr_pd + 43);

    auto grr_z_y_xz = pbuffer.data(idx_gr_pd + 44);

    auto grr_z_y_yy = pbuffer.data(idx_gr_pd + 45);

    auto grr_z_y_yz = pbuffer.data(idx_gr_pd + 46);

    auto grr_z_y_zz = pbuffer.data(idx_gr_pd + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_x, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_y, gr_y_yy, gr_y_yz, gr_y_z, gr_y_zz, grr_z_y_xx, grr_z_y_xy, grr_z_y_xz, grr_z_y_yy, grr_z_y_yz, grr_z_y_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_y_xx[i] = ts_y_xx[i] * gfe_0 * gc_z[i] + gr_y_xx[i] * gc_z[i];

        grr_z_y_xy[i] = ts_y_xy[i] * gfe_0 * gc_z[i] + gr_y_xy[i] * gc_z[i];

        grr_z_y_xz[i] = ts_y_x[i] * gfe2_0 + gr_y_x[i] * gfe_0 + ts_y_xz[i] * gfe_0 * gc_z[i] + gr_y_xz[i] * gc_z[i];

        grr_z_y_yy[i] = ts_y_yy[i] * gfe_0 * gc_z[i] + gr_y_yy[i] * gc_z[i];

        grr_z_y_yz[i] = ts_y_y[i] * gfe2_0 + gr_y_y[i] * gfe_0 + ts_y_yz[i] * gfe_0 * gc_z[i] + gr_y_yz[i] * gc_z[i];

        grr_z_y_zz[i] = 2.0 * ts_y_z[i] * gfe2_0 + 2.0 * gr_y_z[i] * gfe_0 + ts_y_zz[i] * gfe_0 * gc_z[i] + gr_y_zz[i] * gc_z[i];
    }

    // Set up 48-54 components of targeted buffer : PD

    auto grr_z_z_xx = pbuffer.data(idx_gr_pd + 48);

    auto grr_z_z_xy = pbuffer.data(idx_gr_pd + 49);

    auto grr_z_z_xz = pbuffer.data(idx_gr_pd + 50);

    auto grr_z_z_yy = pbuffer.data(idx_gr_pd + 51);

    auto grr_z_z_yz = pbuffer.data(idx_gr_pd + 52);

    auto grr_z_z_zz = pbuffer.data(idx_gr_pd + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xx, gr_0_xy, gr_0_xz, gr_0_yy, gr_0_yz, gr_0_zz, gr_z_x, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_y, gr_z_yy, gr_z_yz, gr_z_z, gr_z_zz, grr_z_z_xx, grr_z_z_xy, grr_z_z_xz, grr_z_z_yy, grr_z_z_yz, grr_z_z_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_z_xx[i] = ts_0_xx[i] * gfe2_0 + gr_0_xx[i] * gfe_0 + ts_z_xx[i] * gfe_0 * gc_z[i] + gr_z_xx[i] * gc_z[i];

        grr_z_z_xy[i] = ts_0_xy[i] * gfe2_0 + gr_0_xy[i] * gfe_0 + ts_z_xy[i] * gfe_0 * gc_z[i] + gr_z_xy[i] * gc_z[i];

        grr_z_z_xz[i] = ts_0_xz[i] * gfe2_0 + gr_0_xz[i] * gfe_0 + ts_z_x[i] * gfe2_0 + gr_z_x[i] * gfe_0 + ts_z_xz[i] * gfe_0 * gc_z[i] + gr_z_xz[i] * gc_z[i];

        grr_z_z_yy[i] = ts_0_yy[i] * gfe2_0 + gr_0_yy[i] * gfe_0 + ts_z_yy[i] * gfe_0 * gc_z[i] + gr_z_yy[i] * gc_z[i];

        grr_z_z_yz[i] = ts_0_yz[i] * gfe2_0 + gr_0_yz[i] * gfe_0 + ts_z_y[i] * gfe2_0 + gr_z_y[i] * gfe_0 + ts_z_yz[i] * gfe_0 * gc_z[i] + gr_z_yz[i] * gc_z[i];

        grr_z_z_zz[i] = ts_0_zz[i] * gfe2_0 + gr_0_zz[i] * gfe_0 + 2.0 * ts_z_z[i] * gfe2_0 + 2.0 * gr_z_z[i] * gfe_0 + ts_z_zz[i] * gfe_0 * gc_z[i] + gr_z_zz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

