#include "ThreeCenterOverlapGradientPrimRecDD.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_dd(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_dd,
                              const size_t idx_pd,
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

    auto gs_x_xx_xx = pbuffer.data(idx_g_dd);

    auto gs_x_xx_xy = pbuffer.data(idx_g_dd + 1);

    auto gs_x_xx_xz = pbuffer.data(idx_g_dd + 2);

    auto gs_x_xx_yy = pbuffer.data(idx_g_dd + 3);

    auto gs_x_xx_yz = pbuffer.data(idx_g_dd + 4);

    auto gs_x_xx_zz = pbuffer.data(idx_g_dd + 5);

    #pragma omp simd aligned(gc_x, gs_x_xx_xx, gs_x_xx_xy, gs_x_xx_xz, gs_x_xx_yy, gs_x_xx_yz, gs_x_xx_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xx_xx[i] = 4.0 * ts_x_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xx[i] * gc_x[i] * tce_0;

        gs_x_xx_xy[i] = 4.0 * ts_x_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xy[i] * gc_x[i] * tce_0;

        gs_x_xx_xz[i] = 4.0 * ts_x_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xz[i] * gc_x[i] * tce_0;

        gs_x_xx_yy[i] = 4.0 * ts_x_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yy[i] * gc_x[i] * tce_0;

        gs_x_xx_yz[i] = 4.0 * ts_x_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yz[i] * gc_x[i] * tce_0;

        gs_x_xx_zz[i] = 4.0 * ts_x_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 6-12 components of targeted buffer : DD

    auto gs_x_xy_xx = pbuffer.data(idx_g_dd + 6);

    auto gs_x_xy_xy = pbuffer.data(idx_g_dd + 7);

    auto gs_x_xy_xz = pbuffer.data(idx_g_dd + 8);

    auto gs_x_xy_yy = pbuffer.data(idx_g_dd + 9);

    auto gs_x_xy_yz = pbuffer.data(idx_g_dd + 10);

    auto gs_x_xy_zz = pbuffer.data(idx_g_dd + 11);

    #pragma omp simd aligned(gc_x, gs_x_xy_xx, gs_x_xy_xy, gs_x_xy_xz, gs_x_xy_yy, gs_x_xy_yz, gs_x_xy_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xy_xx[i] = 2.0 * ts_y_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xx[i] * gc_x[i] * tce_0;

        gs_x_xy_xy[i] = 2.0 * ts_y_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xy[i] * gc_x[i] * tce_0;

        gs_x_xy_xz[i] = 2.0 * ts_y_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xz[i] * gc_x[i] * tce_0;

        gs_x_xy_yy[i] = 2.0 * ts_y_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yy[i] * gc_x[i] * tce_0;

        gs_x_xy_yz[i] = 2.0 * ts_y_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yz[i] * gc_x[i] * tce_0;

        gs_x_xy_zz[i] = 2.0 * ts_y_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 12-18 components of targeted buffer : DD

    auto gs_x_xz_xx = pbuffer.data(idx_g_dd + 12);

    auto gs_x_xz_xy = pbuffer.data(idx_g_dd + 13);

    auto gs_x_xz_xz = pbuffer.data(idx_g_dd + 14);

    auto gs_x_xz_yy = pbuffer.data(idx_g_dd + 15);

    auto gs_x_xz_yz = pbuffer.data(idx_g_dd + 16);

    auto gs_x_xz_zz = pbuffer.data(idx_g_dd + 17);

    #pragma omp simd aligned(gc_x, gs_x_xz_xx, gs_x_xz_xy, gs_x_xz_xz, gs_x_xz_yy, gs_x_xz_yz, gs_x_xz_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xz_xx[i] = 2.0 * ts_z_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xx[i] * gc_x[i] * tce_0;

        gs_x_xz_xy[i] = 2.0 * ts_z_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xy[i] * gc_x[i] * tce_0;

        gs_x_xz_xz[i] = 2.0 * ts_z_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xz[i] * gc_x[i] * tce_0;

        gs_x_xz_yy[i] = 2.0 * ts_z_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yy[i] * gc_x[i] * tce_0;

        gs_x_xz_yz[i] = 2.0 * ts_z_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yz[i] * gc_x[i] * tce_0;

        gs_x_xz_zz[i] = 2.0 * ts_z_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 18-24 components of targeted buffer : DD

    auto gs_x_yy_xx = pbuffer.data(idx_g_dd + 18);

    auto gs_x_yy_xy = pbuffer.data(idx_g_dd + 19);

    auto gs_x_yy_xz = pbuffer.data(idx_g_dd + 20);

    auto gs_x_yy_yy = pbuffer.data(idx_g_dd + 21);

    auto gs_x_yy_yz = pbuffer.data(idx_g_dd + 22);

    auto gs_x_yy_zz = pbuffer.data(idx_g_dd + 23);

    #pragma omp simd aligned(gc_x, gs_x_yy_xx, gs_x_yy_xy, gs_x_yy_xz, gs_x_yy_yy, gs_x_yy_yz, gs_x_yy_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yy_xx[i] = 4.0 * ts_yy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xx[i] * gc_x[i] * tce_0;

        gs_x_yy_xy[i] = 2.0 * ts_yy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xy[i] * gc_x[i] * tce_0;

        gs_x_yy_xz[i] = 2.0 * ts_yy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xz[i] * gc_x[i] * tce_0;

        gs_x_yy_yy[i] = 2.0 * ts_yy_yy[i] * gc_x[i] * tce_0;

        gs_x_yy_yz[i] = 2.0 * ts_yy_yz[i] * gc_x[i] * tce_0;

        gs_x_yy_zz[i] = 2.0 * ts_yy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 24-30 components of targeted buffer : DD

    auto gs_x_yz_xx = pbuffer.data(idx_g_dd + 24);

    auto gs_x_yz_xy = pbuffer.data(idx_g_dd + 25);

    auto gs_x_yz_xz = pbuffer.data(idx_g_dd + 26);

    auto gs_x_yz_yy = pbuffer.data(idx_g_dd + 27);

    auto gs_x_yz_yz = pbuffer.data(idx_g_dd + 28);

    auto gs_x_yz_zz = pbuffer.data(idx_g_dd + 29);

    #pragma omp simd aligned(gc_x, gs_x_yz_xx, gs_x_yz_xy, gs_x_yz_xz, gs_x_yz_yy, gs_x_yz_yz, gs_x_yz_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yz_xx[i] = 4.0 * ts_yz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xx[i] * gc_x[i] * tce_0;

        gs_x_yz_xy[i] = 2.0 * ts_yz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xy[i] * gc_x[i] * tce_0;

        gs_x_yz_xz[i] = 2.0 * ts_yz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xz[i] * gc_x[i] * tce_0;

        gs_x_yz_yy[i] = 2.0 * ts_yz_yy[i] * gc_x[i] * tce_0;

        gs_x_yz_yz[i] = 2.0 * ts_yz_yz[i] * gc_x[i] * tce_0;

        gs_x_yz_zz[i] = 2.0 * ts_yz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-36 components of targeted buffer : DD

    auto gs_x_zz_xx = pbuffer.data(idx_g_dd + 30);

    auto gs_x_zz_xy = pbuffer.data(idx_g_dd + 31);

    auto gs_x_zz_xz = pbuffer.data(idx_g_dd + 32);

    auto gs_x_zz_yy = pbuffer.data(idx_g_dd + 33);

    auto gs_x_zz_yz = pbuffer.data(idx_g_dd + 34);

    auto gs_x_zz_zz = pbuffer.data(idx_g_dd + 35);

    #pragma omp simd aligned(gc_x, gs_x_zz_xx, gs_x_zz_xy, gs_x_zz_xz, gs_x_zz_yy, gs_x_zz_yz, gs_x_zz_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zz_xx[i] = 4.0 * ts_zz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xx[i] * gc_x[i] * tce_0;

        gs_x_zz_xy[i] = 2.0 * ts_zz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xy[i] * gc_x[i] * tce_0;

        gs_x_zz_xz[i] = 2.0 * ts_zz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xz[i] * gc_x[i] * tce_0;

        gs_x_zz_yy[i] = 2.0 * ts_zz_yy[i] * gc_x[i] * tce_0;

        gs_x_zz_yz[i] = 2.0 * ts_zz_yz[i] * gc_x[i] * tce_0;

        gs_x_zz_zz[i] = 2.0 * ts_zz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 36-42 components of targeted buffer : DD

    auto gs_y_xx_xx = pbuffer.data(idx_g_dd + 36);

    auto gs_y_xx_xy = pbuffer.data(idx_g_dd + 37);

    auto gs_y_xx_xz = pbuffer.data(idx_g_dd + 38);

    auto gs_y_xx_yy = pbuffer.data(idx_g_dd + 39);

    auto gs_y_xx_yz = pbuffer.data(idx_g_dd + 40);

    auto gs_y_xx_zz = pbuffer.data(idx_g_dd + 41);

    #pragma omp simd aligned(gc_y, gs_y_xx_xx, gs_y_xx_xy, gs_y_xx_xz, gs_y_xx_yy, gs_y_xx_yz, gs_y_xx_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xx_xx[i] = 2.0 * ts_xx_xx[i] * gc_y[i] * tce_0;

        gs_y_xx_xy[i] = 2.0 * ts_xx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xy[i] * gc_y[i] * tce_0;

        gs_y_xx_xz[i] = 2.0 * ts_xx_xz[i] * gc_y[i] * tce_0;

        gs_y_xx_yy[i] = 4.0 * ts_xx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yy[i] * gc_y[i] * tce_0;

        gs_y_xx_yz[i] = 2.0 * ts_xx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yz[i] * gc_y[i] * tce_0;

        gs_y_xx_zz[i] = 2.0 * ts_xx_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 42-48 components of targeted buffer : DD

    auto gs_y_xy_xx = pbuffer.data(idx_g_dd + 42);

    auto gs_y_xy_xy = pbuffer.data(idx_g_dd + 43);

    auto gs_y_xy_xz = pbuffer.data(idx_g_dd + 44);

    auto gs_y_xy_yy = pbuffer.data(idx_g_dd + 45);

    auto gs_y_xy_yz = pbuffer.data(idx_g_dd + 46);

    auto gs_y_xy_zz = pbuffer.data(idx_g_dd + 47);

    #pragma omp simd aligned(gc_y, gs_y_xy_xx, gs_y_xy_xy, gs_y_xy_xz, gs_y_xy_yy, gs_y_xy_yz, gs_y_xy_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xy_xx[i] = 2.0 * ts_x_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xx[i] * gc_y[i] * tce_0;

        gs_y_xy_xy[i] = 2.0 * ts_x_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xy[i] * gc_y[i] * tce_0;

        gs_y_xy_xz[i] = 2.0 * ts_x_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xz[i] * gc_y[i] * tce_0;

        gs_y_xy_yy[i] = 2.0 * ts_x_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yy[i] * gc_y[i] * tce_0;

        gs_y_xy_yz[i] = 2.0 * ts_x_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yz[i] * gc_y[i] * tce_0;

        gs_y_xy_zz[i] = 2.0 * ts_x_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 48-54 components of targeted buffer : DD

    auto gs_y_xz_xx = pbuffer.data(idx_g_dd + 48);

    auto gs_y_xz_xy = pbuffer.data(idx_g_dd + 49);

    auto gs_y_xz_xz = pbuffer.data(idx_g_dd + 50);

    auto gs_y_xz_yy = pbuffer.data(idx_g_dd + 51);

    auto gs_y_xz_yz = pbuffer.data(idx_g_dd + 52);

    auto gs_y_xz_zz = pbuffer.data(idx_g_dd + 53);

    #pragma omp simd aligned(gc_y, gs_y_xz_xx, gs_y_xz_xy, gs_y_xz_xz, gs_y_xz_yy, gs_y_xz_yz, gs_y_xz_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xz_xx[i] = 2.0 * ts_xz_xx[i] * gc_y[i] * tce_0;

        gs_y_xz_xy[i] = 2.0 * ts_xz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xy[i] * gc_y[i] * tce_0;

        gs_y_xz_xz[i] = 2.0 * ts_xz_xz[i] * gc_y[i] * tce_0;

        gs_y_xz_yy[i] = 4.0 * ts_xz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yy[i] * gc_y[i] * tce_0;

        gs_y_xz_yz[i] = 2.0 * ts_xz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yz[i] * gc_y[i] * tce_0;

        gs_y_xz_zz[i] = 2.0 * ts_xz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 54-60 components of targeted buffer : DD

    auto gs_y_yy_xx = pbuffer.data(idx_g_dd + 54);

    auto gs_y_yy_xy = pbuffer.data(idx_g_dd + 55);

    auto gs_y_yy_xz = pbuffer.data(idx_g_dd + 56);

    auto gs_y_yy_yy = pbuffer.data(idx_g_dd + 57);

    auto gs_y_yy_yz = pbuffer.data(idx_g_dd + 58);

    auto gs_y_yy_zz = pbuffer.data(idx_g_dd + 59);

    #pragma omp simd aligned(gc_y, gs_y_yy_xx, gs_y_yy_xy, gs_y_yy_xz, gs_y_yy_yy, gs_y_yy_yz, gs_y_yy_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yy_xx[i] = 4.0 * ts_y_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xx[i] * gc_y[i] * tce_0;

        gs_y_yy_xy[i] = 4.0 * ts_y_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xy[i] * gc_y[i] * tce_0;

        gs_y_yy_xz[i] = 4.0 * ts_y_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xz[i] * gc_y[i] * tce_0;

        gs_y_yy_yy[i] = 4.0 * ts_y_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yy[i] * gc_y[i] * tce_0;

        gs_y_yy_yz[i] = 4.0 * ts_y_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yz[i] * gc_y[i] * tce_0;

        gs_y_yy_zz[i] = 4.0 * ts_y_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 60-66 components of targeted buffer : DD

    auto gs_y_yz_xx = pbuffer.data(idx_g_dd + 60);

    auto gs_y_yz_xy = pbuffer.data(idx_g_dd + 61);

    auto gs_y_yz_xz = pbuffer.data(idx_g_dd + 62);

    auto gs_y_yz_yy = pbuffer.data(idx_g_dd + 63);

    auto gs_y_yz_yz = pbuffer.data(idx_g_dd + 64);

    auto gs_y_yz_zz = pbuffer.data(idx_g_dd + 65);

    #pragma omp simd aligned(gc_y, gs_y_yz_xx, gs_y_yz_xy, gs_y_yz_xz, gs_y_yz_yy, gs_y_yz_yz, gs_y_yz_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yz_xx[i] = 2.0 * ts_z_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xx[i] * gc_y[i] * tce_0;

        gs_y_yz_xy[i] = 2.0 * ts_z_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xy[i] * gc_y[i] * tce_0;

        gs_y_yz_xz[i] = 2.0 * ts_z_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xz[i] * gc_y[i] * tce_0;

        gs_y_yz_yy[i] = 2.0 * ts_z_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yy[i] * gc_y[i] * tce_0;

        gs_y_yz_yz[i] = 2.0 * ts_z_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yz[i] * gc_y[i] * tce_0;

        gs_y_yz_zz[i] = 2.0 * ts_z_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 66-72 components of targeted buffer : DD

    auto gs_y_zz_xx = pbuffer.data(idx_g_dd + 66);

    auto gs_y_zz_xy = pbuffer.data(idx_g_dd + 67);

    auto gs_y_zz_xz = pbuffer.data(idx_g_dd + 68);

    auto gs_y_zz_yy = pbuffer.data(idx_g_dd + 69);

    auto gs_y_zz_yz = pbuffer.data(idx_g_dd + 70);

    auto gs_y_zz_zz = pbuffer.data(idx_g_dd + 71);

    #pragma omp simd aligned(gc_y, gs_y_zz_xx, gs_y_zz_xy, gs_y_zz_xz, gs_y_zz_yy, gs_y_zz_yz, gs_y_zz_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zz_xx[i] = 2.0 * ts_zz_xx[i] * gc_y[i] * tce_0;

        gs_y_zz_xy[i] = 2.0 * ts_zz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xy[i] * gc_y[i] * tce_0;

        gs_y_zz_xz[i] = 2.0 * ts_zz_xz[i] * gc_y[i] * tce_0;

        gs_y_zz_yy[i] = 4.0 * ts_zz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yy[i] * gc_y[i] * tce_0;

        gs_y_zz_yz[i] = 2.0 * ts_zz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yz[i] * gc_y[i] * tce_0;

        gs_y_zz_zz[i] = 2.0 * ts_zz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 72-78 components of targeted buffer : DD

    auto gs_z_xx_xx = pbuffer.data(idx_g_dd + 72);

    auto gs_z_xx_xy = pbuffer.data(idx_g_dd + 73);

    auto gs_z_xx_xz = pbuffer.data(idx_g_dd + 74);

    auto gs_z_xx_yy = pbuffer.data(idx_g_dd + 75);

    auto gs_z_xx_yz = pbuffer.data(idx_g_dd + 76);

    auto gs_z_xx_zz = pbuffer.data(idx_g_dd + 77);

    #pragma omp simd aligned(gc_z, gs_z_xx_xx, gs_z_xx_xy, gs_z_xx_xz, gs_z_xx_yy, gs_z_xx_yz, gs_z_xx_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xx_xx[i] = 2.0 * ts_xx_xx[i] * gc_z[i] * tce_0;

        gs_z_xx_xy[i] = 2.0 * ts_xx_xy[i] * gc_z[i] * tce_0;

        gs_z_xx_xz[i] = 2.0 * ts_xx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xz[i] * gc_z[i] * tce_0;

        gs_z_xx_yy[i] = 2.0 * ts_xx_yy[i] * gc_z[i] * tce_0;

        gs_z_xx_yz[i] = 2.0 * ts_xx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yz[i] * gc_z[i] * tce_0;

        gs_z_xx_zz[i] = 4.0 * ts_xx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 78-84 components of targeted buffer : DD

    auto gs_z_xy_xx = pbuffer.data(idx_g_dd + 78);

    auto gs_z_xy_xy = pbuffer.data(idx_g_dd + 79);

    auto gs_z_xy_xz = pbuffer.data(idx_g_dd + 80);

    auto gs_z_xy_yy = pbuffer.data(idx_g_dd + 81);

    auto gs_z_xy_yz = pbuffer.data(idx_g_dd + 82);

    auto gs_z_xy_zz = pbuffer.data(idx_g_dd + 83);

    #pragma omp simd aligned(gc_z, gs_z_xy_xx, gs_z_xy_xy, gs_z_xy_xz, gs_z_xy_yy, gs_z_xy_yz, gs_z_xy_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xy_xx[i] = 2.0 * ts_xy_xx[i] * gc_z[i] * tce_0;

        gs_z_xy_xy[i] = 2.0 * ts_xy_xy[i] * gc_z[i] * tce_0;

        gs_z_xy_xz[i] = 2.0 * ts_xy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xz[i] * gc_z[i] * tce_0;

        gs_z_xy_yy[i] = 2.0 * ts_xy_yy[i] * gc_z[i] * tce_0;

        gs_z_xy_yz[i] = 2.0 * ts_xy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yz[i] * gc_z[i] * tce_0;

        gs_z_xy_zz[i] = 4.0 * ts_xy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 84-90 components of targeted buffer : DD

    auto gs_z_xz_xx = pbuffer.data(idx_g_dd + 84);

    auto gs_z_xz_xy = pbuffer.data(idx_g_dd + 85);

    auto gs_z_xz_xz = pbuffer.data(idx_g_dd + 86);

    auto gs_z_xz_yy = pbuffer.data(idx_g_dd + 87);

    auto gs_z_xz_yz = pbuffer.data(idx_g_dd + 88);

    auto gs_z_xz_zz = pbuffer.data(idx_g_dd + 89);

    #pragma omp simd aligned(gc_z, gs_z_xz_xx, gs_z_xz_xy, gs_z_xz_xz, gs_z_xz_yy, gs_z_xz_yz, gs_z_xz_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xz_xx[i] = 2.0 * ts_x_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xx[i] * gc_z[i] * tce_0;

        gs_z_xz_xy[i] = 2.0 * ts_x_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xy[i] * gc_z[i] * tce_0;

        gs_z_xz_xz[i] = 2.0 * ts_x_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xz[i] * gc_z[i] * tce_0;

        gs_z_xz_yy[i] = 2.0 * ts_x_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yy[i] * gc_z[i] * tce_0;

        gs_z_xz_yz[i] = 2.0 * ts_x_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yz[i] * gc_z[i] * tce_0;

        gs_z_xz_zz[i] = 2.0 * ts_x_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 90-96 components of targeted buffer : DD

    auto gs_z_yy_xx = pbuffer.data(idx_g_dd + 90);

    auto gs_z_yy_xy = pbuffer.data(idx_g_dd + 91);

    auto gs_z_yy_xz = pbuffer.data(idx_g_dd + 92);

    auto gs_z_yy_yy = pbuffer.data(idx_g_dd + 93);

    auto gs_z_yy_yz = pbuffer.data(idx_g_dd + 94);

    auto gs_z_yy_zz = pbuffer.data(idx_g_dd + 95);

    #pragma omp simd aligned(gc_z, gs_z_yy_xx, gs_z_yy_xy, gs_z_yy_xz, gs_z_yy_yy, gs_z_yy_yz, gs_z_yy_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yy_xx[i] = 2.0 * ts_yy_xx[i] * gc_z[i] * tce_0;

        gs_z_yy_xy[i] = 2.0 * ts_yy_xy[i] * gc_z[i] * tce_0;

        gs_z_yy_xz[i] = 2.0 * ts_yy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xz[i] * gc_z[i] * tce_0;

        gs_z_yy_yy[i] = 2.0 * ts_yy_yy[i] * gc_z[i] * tce_0;

        gs_z_yy_yz[i] = 2.0 * ts_yy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yz[i] * gc_z[i] * tce_0;

        gs_z_yy_zz[i] = 4.0 * ts_yy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 96-102 components of targeted buffer : DD

    auto gs_z_yz_xx = pbuffer.data(idx_g_dd + 96);

    auto gs_z_yz_xy = pbuffer.data(idx_g_dd + 97);

    auto gs_z_yz_xz = pbuffer.data(idx_g_dd + 98);

    auto gs_z_yz_yy = pbuffer.data(idx_g_dd + 99);

    auto gs_z_yz_yz = pbuffer.data(idx_g_dd + 100);

    auto gs_z_yz_zz = pbuffer.data(idx_g_dd + 101);

    #pragma omp simd aligned(gc_z, gs_z_yz_xx, gs_z_yz_xy, gs_z_yz_xz, gs_z_yz_yy, gs_z_yz_yz, gs_z_yz_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yz_xx[i] = 2.0 * ts_y_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xx[i] * gc_z[i] * tce_0;

        gs_z_yz_xy[i] = 2.0 * ts_y_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xy[i] * gc_z[i] * tce_0;

        gs_z_yz_xz[i] = 2.0 * ts_y_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xz[i] * gc_z[i] * tce_0;

        gs_z_yz_yy[i] = 2.0 * ts_y_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yy[i] * gc_z[i] * tce_0;

        gs_z_yz_yz[i] = 2.0 * ts_y_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yz[i] * gc_z[i] * tce_0;

        gs_z_yz_zz[i] = 2.0 * ts_y_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 102-108 components of targeted buffer : DD

    auto gs_z_zz_xx = pbuffer.data(idx_g_dd + 102);

    auto gs_z_zz_xy = pbuffer.data(idx_g_dd + 103);

    auto gs_z_zz_xz = pbuffer.data(idx_g_dd + 104);

    auto gs_z_zz_yy = pbuffer.data(idx_g_dd + 105);

    auto gs_z_zz_yz = pbuffer.data(idx_g_dd + 106);

    auto gs_z_zz_zz = pbuffer.data(idx_g_dd + 107);

    #pragma omp simd aligned(gc_z, gs_z_zz_xx, gs_z_zz_xy, gs_z_zz_xz, gs_z_zz_yy, gs_z_zz_yz, gs_z_zz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zz_xx[i] = 4.0 * ts_z_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xx[i] * gc_z[i] * tce_0;

        gs_z_zz_xy[i] = 4.0 * ts_z_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xy[i] * gc_z[i] * tce_0;

        gs_z_zz_xz[i] = 4.0 * ts_z_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xz[i] * gc_z[i] * tce_0;

        gs_z_zz_yy[i] = 4.0 * ts_z_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yy[i] * gc_z[i] * tce_0;

        gs_z_zz_yz[i] = 4.0 * ts_z_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yz[i] * gc_z[i] * tce_0;

        gs_z_zz_zz[i] = 4.0 * ts_z_zz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zz_zz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

