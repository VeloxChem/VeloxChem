#include "ThreeCenterOverlapGradientPrimRecDP.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_dp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_dp,
                              const size_t idx_pp,
                              const size_t idx_ds,
                              const size_t idx_dp,
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

    // Set up 0-3 components of targeted buffer : DP

    auto gs_x_xx_x = pbuffer.data(idx_g_dp);

    auto gs_x_xx_y = pbuffer.data(idx_g_dp + 1);

    auto gs_x_xx_z = pbuffer.data(idx_g_dp + 2);

    #pragma omp simd aligned(gc_x, gs_x_xx_x, gs_x_xx_y, gs_x_xx_z, ts_x_x, ts_x_y, ts_x_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xx_x[i] = 4.0 * ts_x_x[i] * gfe_0 * tce_0 + 2.0 * ts_xx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xx_x[i] * gc_x[i] * tce_0;

        gs_x_xx_y[i] = 4.0 * ts_x_y[i] * gfe_0 * tce_0 + 2.0 * ts_xx_y[i] * gc_x[i] * tce_0;

        gs_x_xx_z[i] = 4.0 * ts_x_z[i] * gfe_0 * tce_0 + 2.0 * ts_xx_z[i] * gc_x[i] * tce_0;
    }

    // Set up 3-6 components of targeted buffer : DP

    auto gs_x_xy_x = pbuffer.data(idx_g_dp + 3);

    auto gs_x_xy_y = pbuffer.data(idx_g_dp + 4);

    auto gs_x_xy_z = pbuffer.data(idx_g_dp + 5);

    #pragma omp simd aligned(gc_x, gs_x_xy_x, gs_x_xy_y, gs_x_xy_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xy_x[i] = 2.0 * ts_y_x[i] * gfe_0 * tce_0 + 2.0 * ts_xy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xy_x[i] * gc_x[i] * tce_0;

        gs_x_xy_y[i] = 2.0 * ts_y_y[i] * gfe_0 * tce_0 + 2.0 * ts_xy_y[i] * gc_x[i] * tce_0;

        gs_x_xy_z[i] = 2.0 * ts_y_z[i] * gfe_0 * tce_0 + 2.0 * ts_xy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 6-9 components of targeted buffer : DP

    auto gs_x_xz_x = pbuffer.data(idx_g_dp + 6);

    auto gs_x_xz_y = pbuffer.data(idx_g_dp + 7);

    auto gs_x_xz_z = pbuffer.data(idx_g_dp + 8);

    #pragma omp simd aligned(gc_x, gs_x_xz_x, gs_x_xz_y, gs_x_xz_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xz_x[i] = 2.0 * ts_z_x[i] * gfe_0 * tce_0 + 2.0 * ts_xz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xz_x[i] * gc_x[i] * tce_0;

        gs_x_xz_y[i] = 2.0 * ts_z_y[i] * gfe_0 * tce_0 + 2.0 * ts_xz_y[i] * gc_x[i] * tce_0;

        gs_x_xz_z[i] = 2.0 * ts_z_z[i] * gfe_0 * tce_0 + 2.0 * ts_xz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 9-12 components of targeted buffer : DP

    auto gs_x_yy_x = pbuffer.data(idx_g_dp + 9);

    auto gs_x_yy_y = pbuffer.data(idx_g_dp + 10);

    auto gs_x_yy_z = pbuffer.data(idx_g_dp + 11);

    #pragma omp simd aligned(gc_x, gs_x_yy_x, gs_x_yy_y, gs_x_yy_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yy_x[i] = 2.0 * ts_yy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yy_x[i] * gc_x[i] * tce_0;

        gs_x_yy_y[i] = 2.0 * ts_yy_y[i] * gc_x[i] * tce_0;

        gs_x_yy_z[i] = 2.0 * ts_yy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 12-15 components of targeted buffer : DP

    auto gs_x_yz_x = pbuffer.data(idx_g_dp + 12);

    auto gs_x_yz_y = pbuffer.data(idx_g_dp + 13);

    auto gs_x_yz_z = pbuffer.data(idx_g_dp + 14);

    #pragma omp simd aligned(gc_x, gs_x_yz_x, gs_x_yz_y, gs_x_yz_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yz_x[i] = 2.0 * ts_yz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yz_x[i] * gc_x[i] * tce_0;

        gs_x_yz_y[i] = 2.0 * ts_yz_y[i] * gc_x[i] * tce_0;

        gs_x_yz_z[i] = 2.0 * ts_yz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 15-18 components of targeted buffer : DP

    auto gs_x_zz_x = pbuffer.data(idx_g_dp + 15);

    auto gs_x_zz_y = pbuffer.data(idx_g_dp + 16);

    auto gs_x_zz_z = pbuffer.data(idx_g_dp + 17);

    #pragma omp simd aligned(gc_x, gs_x_zz_x, gs_x_zz_y, gs_x_zz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zz_x[i] = 2.0 * ts_zz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zz_x[i] * gc_x[i] * tce_0;

        gs_x_zz_y[i] = 2.0 * ts_zz_y[i] * gc_x[i] * tce_0;

        gs_x_zz_z[i] = 2.0 * ts_zz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 18-21 components of targeted buffer : DP

    auto gs_y_xx_x = pbuffer.data(idx_g_dp + 18);

    auto gs_y_xx_y = pbuffer.data(idx_g_dp + 19);

    auto gs_y_xx_z = pbuffer.data(idx_g_dp + 20);

    #pragma omp simd aligned(gc_y, gs_y_xx_x, gs_y_xx_y, gs_y_xx_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xx_x[i] = 2.0 * ts_xx_x[i] * gc_y[i] * tce_0;

        gs_y_xx_y[i] = 2.0 * ts_xx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xx_y[i] * gc_y[i] * tce_0;

        gs_y_xx_z[i] = 2.0 * ts_xx_z[i] * gc_y[i] * tce_0;
    }

    // Set up 21-24 components of targeted buffer : DP

    auto gs_y_xy_x = pbuffer.data(idx_g_dp + 21);

    auto gs_y_xy_y = pbuffer.data(idx_g_dp + 22);

    auto gs_y_xy_z = pbuffer.data(idx_g_dp + 23);

    #pragma omp simd aligned(gc_y, gs_y_xy_x, gs_y_xy_y, gs_y_xy_z, ts_x_x, ts_x_y, ts_x_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xy_x[i] = 2.0 * ts_x_x[i] * gfe_0 * tce_0 + 2.0 * ts_xy_x[i] * gc_y[i] * tce_0;

        gs_y_xy_y[i] = 2.0 * ts_x_y[i] * gfe_0 * tce_0 + 2.0 * ts_xy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xy_y[i] * gc_y[i] * tce_0;

        gs_y_xy_z[i] = 2.0 * ts_x_z[i] * gfe_0 * tce_0 + 2.0 * ts_xy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 24-27 components of targeted buffer : DP

    auto gs_y_xz_x = pbuffer.data(idx_g_dp + 24);

    auto gs_y_xz_y = pbuffer.data(idx_g_dp + 25);

    auto gs_y_xz_z = pbuffer.data(idx_g_dp + 26);

    #pragma omp simd aligned(gc_y, gs_y_xz_x, gs_y_xz_y, gs_y_xz_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xz_x[i] = 2.0 * ts_xz_x[i] * gc_y[i] * tce_0;

        gs_y_xz_y[i] = 2.0 * ts_xz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xz_y[i] * gc_y[i] * tce_0;

        gs_y_xz_z[i] = 2.0 * ts_xz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 27-30 components of targeted buffer : DP

    auto gs_y_yy_x = pbuffer.data(idx_g_dp + 27);

    auto gs_y_yy_y = pbuffer.data(idx_g_dp + 28);

    auto gs_y_yy_z = pbuffer.data(idx_g_dp + 29);

    #pragma omp simd aligned(gc_y, gs_y_yy_x, gs_y_yy_y, gs_y_yy_z, ts_y_x, ts_y_y, ts_y_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yy_x[i] = 4.0 * ts_y_x[i] * gfe_0 * tce_0 + 2.0 * ts_yy_x[i] * gc_y[i] * tce_0;

        gs_y_yy_y[i] = 4.0 * ts_y_y[i] * gfe_0 * tce_0 + 2.0 * ts_yy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yy_y[i] * gc_y[i] * tce_0;

        gs_y_yy_z[i] = 4.0 * ts_y_z[i] * gfe_0 * tce_0 + 2.0 * ts_yy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 30-33 components of targeted buffer : DP

    auto gs_y_yz_x = pbuffer.data(idx_g_dp + 30);

    auto gs_y_yz_y = pbuffer.data(idx_g_dp + 31);

    auto gs_y_yz_z = pbuffer.data(idx_g_dp + 32);

    #pragma omp simd aligned(gc_y, gs_y_yz_x, gs_y_yz_y, gs_y_yz_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yz_x[i] = 2.0 * ts_z_x[i] * gfe_0 * tce_0 + 2.0 * ts_yz_x[i] * gc_y[i] * tce_0;

        gs_y_yz_y[i] = 2.0 * ts_z_y[i] * gfe_0 * tce_0 + 2.0 * ts_yz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yz_y[i] * gc_y[i] * tce_0;

        gs_y_yz_z[i] = 2.0 * ts_z_z[i] * gfe_0 * tce_0 + 2.0 * ts_yz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 33-36 components of targeted buffer : DP

    auto gs_y_zz_x = pbuffer.data(idx_g_dp + 33);

    auto gs_y_zz_y = pbuffer.data(idx_g_dp + 34);

    auto gs_y_zz_z = pbuffer.data(idx_g_dp + 35);

    #pragma omp simd aligned(gc_y, gs_y_zz_x, gs_y_zz_y, gs_y_zz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zz_x[i] = 2.0 * ts_zz_x[i] * gc_y[i] * tce_0;

        gs_y_zz_y[i] = 2.0 * ts_zz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zz_y[i] * gc_y[i] * tce_0;

        gs_y_zz_z[i] = 2.0 * ts_zz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 36-39 components of targeted buffer : DP

    auto gs_z_xx_x = pbuffer.data(idx_g_dp + 36);

    auto gs_z_xx_y = pbuffer.data(idx_g_dp + 37);

    auto gs_z_xx_z = pbuffer.data(idx_g_dp + 38);

    #pragma omp simd aligned(gc_z, gs_z_xx_x, gs_z_xx_y, gs_z_xx_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xx_x[i] = 2.0 * ts_xx_x[i] * gc_z[i] * tce_0;

        gs_z_xx_y[i] = 2.0 * ts_xx_y[i] * gc_z[i] * tce_0;

        gs_z_xx_z[i] = 2.0 * ts_xx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xx_z[i] * gc_z[i] * tce_0;
    }

    // Set up 39-42 components of targeted buffer : DP

    auto gs_z_xy_x = pbuffer.data(idx_g_dp + 39);

    auto gs_z_xy_y = pbuffer.data(idx_g_dp + 40);

    auto gs_z_xy_z = pbuffer.data(idx_g_dp + 41);

    #pragma omp simd aligned(gc_z, gs_z_xy_x, gs_z_xy_y, gs_z_xy_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xy_x[i] = 2.0 * ts_xy_x[i] * gc_z[i] * tce_0;

        gs_z_xy_y[i] = 2.0 * ts_xy_y[i] * gc_z[i] * tce_0;

        gs_z_xy_z[i] = 2.0 * ts_xy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 42-45 components of targeted buffer : DP

    auto gs_z_xz_x = pbuffer.data(idx_g_dp + 42);

    auto gs_z_xz_y = pbuffer.data(idx_g_dp + 43);

    auto gs_z_xz_z = pbuffer.data(idx_g_dp + 44);

    #pragma omp simd aligned(gc_z, gs_z_xz_x, gs_z_xz_y, gs_z_xz_z, ts_x_x, ts_x_y, ts_x_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xz_x[i] = 2.0 * ts_x_x[i] * gfe_0 * tce_0 + 2.0 * ts_xz_x[i] * gc_z[i] * tce_0;

        gs_z_xz_y[i] = 2.0 * ts_x_y[i] * gfe_0 * tce_0 + 2.0 * ts_xz_y[i] * gc_z[i] * tce_0;

        gs_z_xz_z[i] = 2.0 * ts_x_z[i] * gfe_0 * tce_0 + 2.0 * ts_xz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 45-48 components of targeted buffer : DP

    auto gs_z_yy_x = pbuffer.data(idx_g_dp + 45);

    auto gs_z_yy_y = pbuffer.data(idx_g_dp + 46);

    auto gs_z_yy_z = pbuffer.data(idx_g_dp + 47);

    #pragma omp simd aligned(gc_z, gs_z_yy_x, gs_z_yy_y, gs_z_yy_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yy_x[i] = 2.0 * ts_yy_x[i] * gc_z[i] * tce_0;

        gs_z_yy_y[i] = 2.0 * ts_yy_y[i] * gc_z[i] * tce_0;

        gs_z_yy_z[i] = 2.0 * ts_yy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 48-51 components of targeted buffer : DP

    auto gs_z_yz_x = pbuffer.data(idx_g_dp + 48);

    auto gs_z_yz_y = pbuffer.data(idx_g_dp + 49);

    auto gs_z_yz_z = pbuffer.data(idx_g_dp + 50);

    #pragma omp simd aligned(gc_z, gs_z_yz_x, gs_z_yz_y, gs_z_yz_z, ts_y_x, ts_y_y, ts_y_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yz_x[i] = 2.0 * ts_y_x[i] * gfe_0 * tce_0 + 2.0 * ts_yz_x[i] * gc_z[i] * tce_0;

        gs_z_yz_y[i] = 2.0 * ts_y_y[i] * gfe_0 * tce_0 + 2.0 * ts_yz_y[i] * gc_z[i] * tce_0;

        gs_z_yz_z[i] = 2.0 * ts_y_z[i] * gfe_0 * tce_0 + 2.0 * ts_yz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 51-54 components of targeted buffer : DP

    auto gs_z_zz_x = pbuffer.data(idx_g_dp + 51);

    auto gs_z_zz_y = pbuffer.data(idx_g_dp + 52);

    auto gs_z_zz_z = pbuffer.data(idx_g_dp + 53);

    #pragma omp simd aligned(gc_z, gs_z_zz_x, gs_z_zz_y, gs_z_zz_z, ts_z_x, ts_z_y, ts_z_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zz_x[i] = 4.0 * ts_z_x[i] * gfe_0 * tce_0 + 2.0 * ts_zz_x[i] * gc_z[i] * tce_0;

        gs_z_zz_y[i] = 4.0 * ts_z_y[i] * gfe_0 * tce_0 + 2.0 * ts_zz_y[i] * gc_z[i] * tce_0;

        gs_z_zz_z[i] = 4.0 * ts_z_z[i] * gfe_0 * tce_0 + 2.0 * ts_zz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zz_z[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

