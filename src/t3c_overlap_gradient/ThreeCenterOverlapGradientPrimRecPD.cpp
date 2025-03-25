#include "ThreeCenterOverlapGradientPrimRecPD.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_pd(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_pd,
                              const size_t idx_sd,
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

    // Set up 0-6 components of targeted buffer : PD

    auto gs_x_x_xx = pbuffer.data(idx_g_pd);

    auto gs_x_x_xy = pbuffer.data(idx_g_pd + 1);

    auto gs_x_x_xz = pbuffer.data(idx_g_pd + 2);

    auto gs_x_x_yy = pbuffer.data(idx_g_pd + 3);

    auto gs_x_x_yz = pbuffer.data(idx_g_pd + 4);

    auto gs_x_x_zz = pbuffer.data(idx_g_pd + 5);

    #pragma omp simd aligned(gc_x, gs_x_x_xx, gs_x_x_xy, gs_x_x_xz, gs_x_x_yy, gs_x_x_yz, gs_x_x_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_x_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 * tce_0 + 4.0 * ts_x_x[i] * gfe_0 * tce_0 + 2.0 * ts_x_xx[i] * gc_x[i] * tce_0;

        gs_x_x_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 * tce_0 + 2.0 * ts_x_y[i] * gfe_0 * tce_0 + 2.0 * ts_x_xy[i] * gc_x[i] * tce_0;

        gs_x_x_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 * tce_0 + 2.0 * ts_x_z[i] * gfe_0 * tce_0 + 2.0 * ts_x_xz[i] * gc_x[i] * tce_0;

        gs_x_x_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yy[i] * gc_x[i] * tce_0;

        gs_x_x_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yz[i] * gc_x[i] * tce_0;

        gs_x_x_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 6-12 components of targeted buffer : PD

    auto gs_x_y_xx = pbuffer.data(idx_g_pd + 6);

    auto gs_x_y_xy = pbuffer.data(idx_g_pd + 7);

    auto gs_x_y_xz = pbuffer.data(idx_g_pd + 8);

    auto gs_x_y_yy = pbuffer.data(idx_g_pd + 9);

    auto gs_x_y_yz = pbuffer.data(idx_g_pd + 10);

    auto gs_x_y_zz = pbuffer.data(idx_g_pd + 11);

    #pragma omp simd aligned(gc_x, gs_x_y_xx, gs_x_y_xy, gs_x_y_xz, gs_x_y_yy, gs_x_y_yz, gs_x_y_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_y_xx[i] = 4.0 * ts_y_x[i] * gfe_0 * tce_0 + 2.0 * ts_y_xx[i] * gc_x[i] * tce_0;

        gs_x_y_xy[i] = 2.0 * ts_y_y[i] * gfe_0 * tce_0 + 2.0 * ts_y_xy[i] * gc_x[i] * tce_0;

        gs_x_y_xz[i] = 2.0 * ts_y_z[i] * gfe_0 * tce_0 + 2.0 * ts_y_xz[i] * gc_x[i] * tce_0;

        gs_x_y_yy[i] = 2.0 * ts_y_yy[i] * gc_x[i] * tce_0;

        gs_x_y_yz[i] = 2.0 * ts_y_yz[i] * gc_x[i] * tce_0;

        gs_x_y_zz[i] = 2.0 * ts_y_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 12-18 components of targeted buffer : PD

    auto gs_x_z_xx = pbuffer.data(idx_g_pd + 12);

    auto gs_x_z_xy = pbuffer.data(idx_g_pd + 13);

    auto gs_x_z_xz = pbuffer.data(idx_g_pd + 14);

    auto gs_x_z_yy = pbuffer.data(idx_g_pd + 15);

    auto gs_x_z_yz = pbuffer.data(idx_g_pd + 16);

    auto gs_x_z_zz = pbuffer.data(idx_g_pd + 17);

    #pragma omp simd aligned(gc_x, gs_x_z_xx, gs_x_z_xy, gs_x_z_xz, gs_x_z_yy, gs_x_z_yz, gs_x_z_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_z_xx[i] = 4.0 * ts_z_x[i] * gfe_0 * tce_0 + 2.0 * ts_z_xx[i] * gc_x[i] * tce_0;

        gs_x_z_xy[i] = 2.0 * ts_z_y[i] * gfe_0 * tce_0 + 2.0 * ts_z_xy[i] * gc_x[i] * tce_0;

        gs_x_z_xz[i] = 2.0 * ts_z_z[i] * gfe_0 * tce_0 + 2.0 * ts_z_xz[i] * gc_x[i] * tce_0;

        gs_x_z_yy[i] = 2.0 * ts_z_yy[i] * gc_x[i] * tce_0;

        gs_x_z_yz[i] = 2.0 * ts_z_yz[i] * gc_x[i] * tce_0;

        gs_x_z_zz[i] = 2.0 * ts_z_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 18-24 components of targeted buffer : PD

    auto gs_y_x_xx = pbuffer.data(idx_g_pd + 18);

    auto gs_y_x_xy = pbuffer.data(idx_g_pd + 19);

    auto gs_y_x_xz = pbuffer.data(idx_g_pd + 20);

    auto gs_y_x_yy = pbuffer.data(idx_g_pd + 21);

    auto gs_y_x_yz = pbuffer.data(idx_g_pd + 22);

    auto gs_y_x_zz = pbuffer.data(idx_g_pd + 23);

    #pragma omp simd aligned(gc_y, gs_y_x_xx, gs_y_x_xy, gs_y_x_xz, gs_y_x_yy, gs_y_x_yz, gs_y_x_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_x_xx[i] = 2.0 * ts_x_xx[i] * gc_y[i] * tce_0;

        gs_y_x_xy[i] = 2.0 * ts_x_x[i] * gfe_0 * tce_0 + 2.0 * ts_x_xy[i] * gc_y[i] * tce_0;

        gs_y_x_xz[i] = 2.0 * ts_x_xz[i] * gc_y[i] * tce_0;

        gs_y_x_yy[i] = 4.0 * ts_x_y[i] * gfe_0 * tce_0 + 2.0 * ts_x_yy[i] * gc_y[i] * tce_0;

        gs_y_x_yz[i] = 2.0 * ts_x_z[i] * gfe_0 * tce_0 + 2.0 * ts_x_yz[i] * gc_y[i] * tce_0;

        gs_y_x_zz[i] = 2.0 * ts_x_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 24-30 components of targeted buffer : PD

    auto gs_y_y_xx = pbuffer.data(idx_g_pd + 24);

    auto gs_y_y_xy = pbuffer.data(idx_g_pd + 25);

    auto gs_y_y_xz = pbuffer.data(idx_g_pd + 26);

    auto gs_y_y_yy = pbuffer.data(idx_g_pd + 27);

    auto gs_y_y_yz = pbuffer.data(idx_g_pd + 28);

    auto gs_y_y_zz = pbuffer.data(idx_g_pd + 29);

    #pragma omp simd aligned(gc_y, gs_y_y_xx, gs_y_y_xy, gs_y_y_xz, gs_y_y_yy, gs_y_y_yz, gs_y_y_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_y_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xx[i] * gc_y[i] * tce_0;

        gs_y_y_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 * tce_0 + 2.0 * ts_y_x[i] * gfe_0 * tce_0 + 2.0 * ts_y_xy[i] * gc_y[i] * tce_0;

        gs_y_y_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xz[i] * gc_y[i] * tce_0;

        gs_y_y_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 * tce_0 + 4.0 * ts_y_y[i] * gfe_0 * tce_0 + 2.0 * ts_y_yy[i] * gc_y[i] * tce_0;

        gs_y_y_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 * tce_0 + 2.0 * ts_y_z[i] * gfe_0 * tce_0 + 2.0 * ts_y_yz[i] * gc_y[i] * tce_0;

        gs_y_y_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 30-36 components of targeted buffer : PD

    auto gs_y_z_xx = pbuffer.data(idx_g_pd + 30);

    auto gs_y_z_xy = pbuffer.data(idx_g_pd + 31);

    auto gs_y_z_xz = pbuffer.data(idx_g_pd + 32);

    auto gs_y_z_yy = pbuffer.data(idx_g_pd + 33);

    auto gs_y_z_yz = pbuffer.data(idx_g_pd + 34);

    auto gs_y_z_zz = pbuffer.data(idx_g_pd + 35);

    #pragma omp simd aligned(gc_y, gs_y_z_xx, gs_y_z_xy, gs_y_z_xz, gs_y_z_yy, gs_y_z_yz, gs_y_z_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_z_xx[i] = 2.0 * ts_z_xx[i] * gc_y[i] * tce_0;

        gs_y_z_xy[i] = 2.0 * ts_z_x[i] * gfe_0 * tce_0 + 2.0 * ts_z_xy[i] * gc_y[i] * tce_0;

        gs_y_z_xz[i] = 2.0 * ts_z_xz[i] * gc_y[i] * tce_0;

        gs_y_z_yy[i] = 4.0 * ts_z_y[i] * gfe_0 * tce_0 + 2.0 * ts_z_yy[i] * gc_y[i] * tce_0;

        gs_y_z_yz[i] = 2.0 * ts_z_z[i] * gfe_0 * tce_0 + 2.0 * ts_z_yz[i] * gc_y[i] * tce_0;

        gs_y_z_zz[i] = 2.0 * ts_z_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 36-42 components of targeted buffer : PD

    auto gs_z_x_xx = pbuffer.data(idx_g_pd + 36);

    auto gs_z_x_xy = pbuffer.data(idx_g_pd + 37);

    auto gs_z_x_xz = pbuffer.data(idx_g_pd + 38);

    auto gs_z_x_yy = pbuffer.data(idx_g_pd + 39);

    auto gs_z_x_yz = pbuffer.data(idx_g_pd + 40);

    auto gs_z_x_zz = pbuffer.data(idx_g_pd + 41);

    #pragma omp simd aligned(gc_z, gs_z_x_xx, gs_z_x_xy, gs_z_x_xz, gs_z_x_yy, gs_z_x_yz, gs_z_x_zz, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_x_xx[i] = 2.0 * ts_x_xx[i] * gc_z[i] * tce_0;

        gs_z_x_xy[i] = 2.0 * ts_x_xy[i] * gc_z[i] * tce_0;

        gs_z_x_xz[i] = 2.0 * ts_x_x[i] * gfe_0 * tce_0 + 2.0 * ts_x_xz[i] * gc_z[i] * tce_0;

        gs_z_x_yy[i] = 2.0 * ts_x_yy[i] * gc_z[i] * tce_0;

        gs_z_x_yz[i] = 2.0 * ts_x_y[i] * gfe_0 * tce_0 + 2.0 * ts_x_yz[i] * gc_z[i] * tce_0;

        gs_z_x_zz[i] = 4.0 * ts_x_z[i] * gfe_0 * tce_0 + 2.0 * ts_x_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 42-48 components of targeted buffer : PD

    auto gs_z_y_xx = pbuffer.data(idx_g_pd + 42);

    auto gs_z_y_xy = pbuffer.data(idx_g_pd + 43);

    auto gs_z_y_xz = pbuffer.data(idx_g_pd + 44);

    auto gs_z_y_yy = pbuffer.data(idx_g_pd + 45);

    auto gs_z_y_yz = pbuffer.data(idx_g_pd + 46);

    auto gs_z_y_zz = pbuffer.data(idx_g_pd + 47);

    #pragma omp simd aligned(gc_z, gs_z_y_xx, gs_z_y_xy, gs_z_y_xz, gs_z_y_yy, gs_z_y_yz, gs_z_y_zz, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_y_xx[i] = 2.0 * ts_y_xx[i] * gc_z[i] * tce_0;

        gs_z_y_xy[i] = 2.0 * ts_y_xy[i] * gc_z[i] * tce_0;

        gs_z_y_xz[i] = 2.0 * ts_y_x[i] * gfe_0 * tce_0 + 2.0 * ts_y_xz[i] * gc_z[i] * tce_0;

        gs_z_y_yy[i] = 2.0 * ts_y_yy[i] * gc_z[i] * tce_0;

        gs_z_y_yz[i] = 2.0 * ts_y_y[i] * gfe_0 * tce_0 + 2.0 * ts_y_yz[i] * gc_z[i] * tce_0;

        gs_z_y_zz[i] = 4.0 * ts_y_z[i] * gfe_0 * tce_0 + 2.0 * ts_y_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 48-54 components of targeted buffer : PD

    auto gs_z_z_xx = pbuffer.data(idx_g_pd + 48);

    auto gs_z_z_xy = pbuffer.data(idx_g_pd + 49);

    auto gs_z_z_xz = pbuffer.data(idx_g_pd + 50);

    auto gs_z_z_yy = pbuffer.data(idx_g_pd + 51);

    auto gs_z_z_yz = pbuffer.data(idx_g_pd + 52);

    auto gs_z_z_zz = pbuffer.data(idx_g_pd + 53);

    #pragma omp simd aligned(gc_z, gs_z_z_xx, gs_z_z_xy, gs_z_z_xz, gs_z_z_yy, gs_z_z_yz, gs_z_z_zz, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_yy, ts_0_yz, ts_0_zz, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_z_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xx[i] * gc_z[i] * tce_0;

        gs_z_z_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xy[i] * gc_z[i] * tce_0;

        gs_z_z_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 * tce_0 + 2.0 * ts_z_x[i] * gfe_0 * tce_0 + 2.0 * ts_z_xz[i] * gc_z[i] * tce_0;

        gs_z_z_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yy[i] * gc_z[i] * tce_0;

        gs_z_z_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 * tce_0 + 2.0 * ts_z_y[i] * gfe_0 * tce_0 + 2.0 * ts_z_yz[i] * gc_z[i] * tce_0;

        gs_z_z_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 * tce_0 + 4.0 * ts_z_z[i] * gfe_0 * tce_0 + 2.0 * ts_z_zz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

