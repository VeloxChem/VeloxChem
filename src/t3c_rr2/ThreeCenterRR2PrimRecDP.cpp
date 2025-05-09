#include "ThreeCenterRR2PrimRecDP.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_dp(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_dp,
                  const size_t idx_pp,
                  const size_t idx_g_pp,
                  const size_t idx_ds,
                  const size_t idx_g_ds,
                  const size_t idx_dp,
                  const size_t idx_g_dp,
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

    // Set up components of auxiliary buffer : DP

    auto gr_xx_x = pbuffer.data(idx_g_dp);

    auto gr_xx_y = pbuffer.data(idx_g_dp + 1);

    auto gr_xx_z = pbuffer.data(idx_g_dp + 2);

    auto gr_xy_x = pbuffer.data(idx_g_dp + 3);

    auto gr_xy_y = pbuffer.data(idx_g_dp + 4);

    auto gr_xy_z = pbuffer.data(idx_g_dp + 5);

    auto gr_xz_x = pbuffer.data(idx_g_dp + 6);

    auto gr_xz_y = pbuffer.data(idx_g_dp + 7);

    auto gr_xz_z = pbuffer.data(idx_g_dp + 8);

    auto gr_yy_x = pbuffer.data(idx_g_dp + 9);

    auto gr_yy_y = pbuffer.data(idx_g_dp + 10);

    auto gr_yy_z = pbuffer.data(idx_g_dp + 11);

    auto gr_yz_x = pbuffer.data(idx_g_dp + 12);

    auto gr_yz_y = pbuffer.data(idx_g_dp + 13);

    auto gr_yz_z = pbuffer.data(idx_g_dp + 14);

    auto gr_zz_x = pbuffer.data(idx_g_dp + 15);

    auto gr_zz_y = pbuffer.data(idx_g_dp + 16);

    auto gr_zz_z = pbuffer.data(idx_g_dp + 17);

    // Set up 0-3 components of targeted buffer : DP

    auto grr_x_xx_x = pbuffer.data(idx_gr_dp);

    auto grr_x_xx_y = pbuffer.data(idx_gr_dp + 1);

    auto grr_x_xx_z = pbuffer.data(idx_gr_dp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_x, gr_x_y, gr_x_z, gr_xx_0, gr_xx_x, gr_xx_y, gr_xx_z, grr_x_xx_x, grr_x_xx_y, grr_x_xx_z, ts_x_x, ts_x_y, ts_x_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xx_x[i] = 2.0 * ts_x_x[i] * gfe2_0 + 2.0 * gr_x_x[i] * gfe_0 + ts_xx_0[i] * gfe2_0 + gr_xx_0[i] * gfe_0 + ts_xx_x[i] * gfe_0 * gc_x[i] + gr_xx_x[i] * gc_x[i];

        grr_x_xx_y[i] = 2.0 * ts_x_y[i] * gfe2_0 + 2.0 * gr_x_y[i] * gfe_0 + ts_xx_y[i] * gfe_0 * gc_x[i] + gr_xx_y[i] * gc_x[i];

        grr_x_xx_z[i] = 2.0 * ts_x_z[i] * gfe2_0 + 2.0 * gr_x_z[i] * gfe_0 + ts_xx_z[i] * gfe_0 * gc_x[i] + gr_xx_z[i] * gc_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto grr_x_xy_x = pbuffer.data(idx_gr_dp + 3);

    auto grr_x_xy_y = pbuffer.data(idx_gr_dp + 4);

    auto grr_x_xy_z = pbuffer.data(idx_gr_dp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_0, gr_xy_x, gr_xy_y, gr_xy_z, gr_y_x, gr_y_y, gr_y_z, grr_x_xy_x, grr_x_xy_y, grr_x_xy_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xy_x[i] = ts_y_x[i] * gfe2_0 + gr_y_x[i] * gfe_0 + ts_xy_0[i] * gfe2_0 + gr_xy_0[i] * gfe_0 + ts_xy_x[i] * gfe_0 * gc_x[i] + gr_xy_x[i] * gc_x[i];

        grr_x_xy_y[i] = ts_y_y[i] * gfe2_0 + gr_y_y[i] * gfe_0 + ts_xy_y[i] * gfe_0 * gc_x[i] + gr_xy_y[i] * gc_x[i];

        grr_x_xy_z[i] = ts_y_z[i] * gfe2_0 + gr_y_z[i] * gfe_0 + ts_xy_z[i] * gfe_0 * gc_x[i] + gr_xy_z[i] * gc_x[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto grr_x_xz_x = pbuffer.data(idx_gr_dp + 6);

    auto grr_x_xz_y = pbuffer.data(idx_gr_dp + 7);

    auto grr_x_xz_z = pbuffer.data(idx_gr_dp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_0, gr_xz_x, gr_xz_y, gr_xz_z, gr_z_x, gr_z_y, gr_z_z, grr_x_xz_x, grr_x_xz_y, grr_x_xz_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xz_x[i] = ts_z_x[i] * gfe2_0 + gr_z_x[i] * gfe_0 + ts_xz_0[i] * gfe2_0 + gr_xz_0[i] * gfe_0 + ts_xz_x[i] * gfe_0 * gc_x[i] + gr_xz_x[i] * gc_x[i];

        grr_x_xz_y[i] = ts_z_y[i] * gfe2_0 + gr_z_y[i] * gfe_0 + ts_xz_y[i] * gfe_0 * gc_x[i] + gr_xz_y[i] * gc_x[i];

        grr_x_xz_z[i] = ts_z_z[i] * gfe2_0 + gr_z_z[i] * gfe_0 + ts_xz_z[i] * gfe_0 * gc_x[i] + gr_xz_z[i] * gc_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto grr_x_yy_x = pbuffer.data(idx_gr_dp + 9);

    auto grr_x_yy_y = pbuffer.data(idx_gr_dp + 10);

    auto grr_x_yy_z = pbuffer.data(idx_gr_dp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_0, gr_yy_x, gr_yy_y, gr_yy_z, grr_x_yy_x, grr_x_yy_y, grr_x_yy_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yy_x[i] = ts_yy_0[i] * gfe2_0 + gr_yy_0[i] * gfe_0 + ts_yy_x[i] * gfe_0 * gc_x[i] + gr_yy_x[i] * gc_x[i];

        grr_x_yy_y[i] = ts_yy_y[i] * gfe_0 * gc_x[i] + gr_yy_y[i] * gc_x[i];

        grr_x_yy_z[i] = ts_yy_z[i] * gfe_0 * gc_x[i] + gr_yy_z[i] * gc_x[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto grr_x_yz_x = pbuffer.data(idx_gr_dp + 12);

    auto grr_x_yz_y = pbuffer.data(idx_gr_dp + 13);

    auto grr_x_yz_z = pbuffer.data(idx_gr_dp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_0, gr_yz_x, gr_yz_y, gr_yz_z, grr_x_yz_x, grr_x_yz_y, grr_x_yz_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yz_x[i] = ts_yz_0[i] * gfe2_0 + gr_yz_0[i] * gfe_0 + ts_yz_x[i] * gfe_0 * gc_x[i] + gr_yz_x[i] * gc_x[i];

        grr_x_yz_y[i] = ts_yz_y[i] * gfe_0 * gc_x[i] + gr_yz_y[i] * gc_x[i];

        grr_x_yz_z[i] = ts_yz_z[i] * gfe_0 * gc_x[i] + gr_yz_z[i] * gc_x[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto grr_x_zz_x = pbuffer.data(idx_gr_dp + 15);

    auto grr_x_zz_y = pbuffer.data(idx_gr_dp + 16);

    auto grr_x_zz_z = pbuffer.data(idx_gr_dp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_0, gr_zz_x, gr_zz_y, gr_zz_z, grr_x_zz_x, grr_x_zz_y, grr_x_zz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zz_x[i] = ts_zz_0[i] * gfe2_0 + gr_zz_0[i] * gfe_0 + ts_zz_x[i] * gfe_0 * gc_x[i] + gr_zz_x[i] * gc_x[i];

        grr_x_zz_y[i] = ts_zz_y[i] * gfe_0 * gc_x[i] + gr_zz_y[i] * gc_x[i];

        grr_x_zz_z[i] = ts_zz_z[i] * gfe_0 * gc_x[i] + gr_zz_z[i] * gc_x[i];
    }

    // Set up 18-21 components of targeted buffer : DP

    auto grr_y_xx_x = pbuffer.data(idx_gr_dp + 18);

    auto grr_y_xx_y = pbuffer.data(idx_gr_dp + 19);

    auto grr_y_xx_z = pbuffer.data(idx_gr_dp + 20);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_0, gr_xx_x, gr_xx_y, gr_xx_z, grr_y_xx_x, grr_y_xx_y, grr_y_xx_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xx_x[i] = ts_xx_x[i] * gfe_0 * gc_y[i] + gr_xx_x[i] * gc_y[i];

        grr_y_xx_y[i] = ts_xx_0[i] * gfe2_0 + gr_xx_0[i] * gfe_0 + ts_xx_y[i] * gfe_0 * gc_y[i] + gr_xx_y[i] * gc_y[i];

        grr_y_xx_z[i] = ts_xx_z[i] * gfe_0 * gc_y[i] + gr_xx_z[i] * gc_y[i];
    }

    // Set up 21-24 components of targeted buffer : DP

    auto grr_y_xy_x = pbuffer.data(idx_gr_dp + 21);

    auto grr_y_xy_y = pbuffer.data(idx_gr_dp + 22);

    auto grr_y_xy_z = pbuffer.data(idx_gr_dp + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_x, gr_x_y, gr_x_z, gr_xy_0, gr_xy_x, gr_xy_y, gr_xy_z, grr_y_xy_x, grr_y_xy_y, grr_y_xy_z, ts_x_x, ts_x_y, ts_x_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xy_x[i] = ts_x_x[i] * gfe2_0 + gr_x_x[i] * gfe_0 + ts_xy_x[i] * gfe_0 * gc_y[i] + gr_xy_x[i] * gc_y[i];

        grr_y_xy_y[i] = ts_x_y[i] * gfe2_0 + gr_x_y[i] * gfe_0 + ts_xy_0[i] * gfe2_0 + gr_xy_0[i] * gfe_0 + ts_xy_y[i] * gfe_0 * gc_y[i] + gr_xy_y[i] * gc_y[i];

        grr_y_xy_z[i] = ts_x_z[i] * gfe2_0 + gr_x_z[i] * gfe_0 + ts_xy_z[i] * gfe_0 * gc_y[i] + gr_xy_z[i] * gc_y[i];
    }

    // Set up 24-27 components of targeted buffer : DP

    auto grr_y_xz_x = pbuffer.data(idx_gr_dp + 24);

    auto grr_y_xz_y = pbuffer.data(idx_gr_dp + 25);

    auto grr_y_xz_z = pbuffer.data(idx_gr_dp + 26);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_0, gr_xz_x, gr_xz_y, gr_xz_z, grr_y_xz_x, grr_y_xz_y, grr_y_xz_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xz_x[i] = ts_xz_x[i] * gfe_0 * gc_y[i] + gr_xz_x[i] * gc_y[i];

        grr_y_xz_y[i] = ts_xz_0[i] * gfe2_0 + gr_xz_0[i] * gfe_0 + ts_xz_y[i] * gfe_0 * gc_y[i] + gr_xz_y[i] * gc_y[i];

        grr_y_xz_z[i] = ts_xz_z[i] * gfe_0 * gc_y[i] + gr_xz_z[i] * gc_y[i];
    }

    // Set up 27-30 components of targeted buffer : DP

    auto grr_y_yy_x = pbuffer.data(idx_gr_dp + 27);

    auto grr_y_yy_y = pbuffer.data(idx_gr_dp + 28);

    auto grr_y_yy_z = pbuffer.data(idx_gr_dp + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_x, gr_y_y, gr_y_z, gr_yy_0, gr_yy_x, gr_yy_y, gr_yy_z, grr_y_yy_x, grr_y_yy_y, grr_y_yy_z, ts_y_x, ts_y_y, ts_y_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yy_x[i] = 2.0 * ts_y_x[i] * gfe2_0 + 2.0 * gr_y_x[i] * gfe_0 + ts_yy_x[i] * gfe_0 * gc_y[i] + gr_yy_x[i] * gc_y[i];

        grr_y_yy_y[i] = 2.0 * ts_y_y[i] * gfe2_0 + 2.0 * gr_y_y[i] * gfe_0 + ts_yy_0[i] * gfe2_0 + gr_yy_0[i] * gfe_0 + ts_yy_y[i] * gfe_0 * gc_y[i] + gr_yy_y[i] * gc_y[i];

        grr_y_yy_z[i] = 2.0 * ts_y_z[i] * gfe2_0 + 2.0 * gr_y_z[i] * gfe_0 + ts_yy_z[i] * gfe_0 * gc_y[i] + gr_yy_z[i] * gc_y[i];
    }

    // Set up 30-33 components of targeted buffer : DP

    auto grr_y_yz_x = pbuffer.data(idx_gr_dp + 30);

    auto grr_y_yz_y = pbuffer.data(idx_gr_dp + 31);

    auto grr_y_yz_z = pbuffer.data(idx_gr_dp + 32);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_0, gr_yz_x, gr_yz_y, gr_yz_z, gr_z_x, gr_z_y, gr_z_z, grr_y_yz_x, grr_y_yz_y, grr_y_yz_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yz_x[i] = ts_z_x[i] * gfe2_0 + gr_z_x[i] * gfe_0 + ts_yz_x[i] * gfe_0 * gc_y[i] + gr_yz_x[i] * gc_y[i];

        grr_y_yz_y[i] = ts_z_y[i] * gfe2_0 + gr_z_y[i] * gfe_0 + ts_yz_0[i] * gfe2_0 + gr_yz_0[i] * gfe_0 + ts_yz_y[i] * gfe_0 * gc_y[i] + gr_yz_y[i] * gc_y[i];

        grr_y_yz_z[i] = ts_z_z[i] * gfe2_0 + gr_z_z[i] * gfe_0 + ts_yz_z[i] * gfe_0 * gc_y[i] + gr_yz_z[i] * gc_y[i];
    }

    // Set up 33-36 components of targeted buffer : DP

    auto grr_y_zz_x = pbuffer.data(idx_gr_dp + 33);

    auto grr_y_zz_y = pbuffer.data(idx_gr_dp + 34);

    auto grr_y_zz_z = pbuffer.data(idx_gr_dp + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_0, gr_zz_x, gr_zz_y, gr_zz_z, grr_y_zz_x, grr_y_zz_y, grr_y_zz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zz_x[i] = ts_zz_x[i] * gfe_0 * gc_y[i] + gr_zz_x[i] * gc_y[i];

        grr_y_zz_y[i] = ts_zz_0[i] * gfe2_0 + gr_zz_0[i] * gfe_0 + ts_zz_y[i] * gfe_0 * gc_y[i] + gr_zz_y[i] * gc_y[i];

        grr_y_zz_z[i] = ts_zz_z[i] * gfe_0 * gc_y[i] + gr_zz_z[i] * gc_y[i];
    }

    // Set up 36-39 components of targeted buffer : DP

    auto grr_z_xx_x = pbuffer.data(idx_gr_dp + 36);

    auto grr_z_xx_y = pbuffer.data(idx_gr_dp + 37);

    auto grr_z_xx_z = pbuffer.data(idx_gr_dp + 38);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_0, gr_xx_x, gr_xx_y, gr_xx_z, grr_z_xx_x, grr_z_xx_y, grr_z_xx_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xx_x[i] = ts_xx_x[i] * gfe_0 * gc_z[i] + gr_xx_x[i] * gc_z[i];

        grr_z_xx_y[i] = ts_xx_y[i] * gfe_0 * gc_z[i] + gr_xx_y[i] * gc_z[i];

        grr_z_xx_z[i] = ts_xx_0[i] * gfe2_0 + gr_xx_0[i] * gfe_0 + ts_xx_z[i] * gfe_0 * gc_z[i] + gr_xx_z[i] * gc_z[i];
    }

    // Set up 39-42 components of targeted buffer : DP

    auto grr_z_xy_x = pbuffer.data(idx_gr_dp + 39);

    auto grr_z_xy_y = pbuffer.data(idx_gr_dp + 40);

    auto grr_z_xy_z = pbuffer.data(idx_gr_dp + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_0, gr_xy_x, gr_xy_y, gr_xy_z, grr_z_xy_x, grr_z_xy_y, grr_z_xy_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xy_x[i] = ts_xy_x[i] * gfe_0 * gc_z[i] + gr_xy_x[i] * gc_z[i];

        grr_z_xy_y[i] = ts_xy_y[i] * gfe_0 * gc_z[i] + gr_xy_y[i] * gc_z[i];

        grr_z_xy_z[i] = ts_xy_0[i] * gfe2_0 + gr_xy_0[i] * gfe_0 + ts_xy_z[i] * gfe_0 * gc_z[i] + gr_xy_z[i] * gc_z[i];
    }

    // Set up 42-45 components of targeted buffer : DP

    auto grr_z_xz_x = pbuffer.data(idx_gr_dp + 42);

    auto grr_z_xz_y = pbuffer.data(idx_gr_dp + 43);

    auto grr_z_xz_z = pbuffer.data(idx_gr_dp + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_x, gr_x_y, gr_x_z, gr_xz_0, gr_xz_x, gr_xz_y, gr_xz_z, grr_z_xz_x, grr_z_xz_y, grr_z_xz_z, ts_x_x, ts_x_y, ts_x_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xz_x[i] = ts_x_x[i] * gfe2_0 + gr_x_x[i] * gfe_0 + ts_xz_x[i] * gfe_0 * gc_z[i] + gr_xz_x[i] * gc_z[i];

        grr_z_xz_y[i] = ts_x_y[i] * gfe2_0 + gr_x_y[i] * gfe_0 + ts_xz_y[i] * gfe_0 * gc_z[i] + gr_xz_y[i] * gc_z[i];

        grr_z_xz_z[i] = ts_x_z[i] * gfe2_0 + gr_x_z[i] * gfe_0 + ts_xz_0[i] * gfe2_0 + gr_xz_0[i] * gfe_0 + ts_xz_z[i] * gfe_0 * gc_z[i] + gr_xz_z[i] * gc_z[i];
    }

    // Set up 45-48 components of targeted buffer : DP

    auto grr_z_yy_x = pbuffer.data(idx_gr_dp + 45);

    auto grr_z_yy_y = pbuffer.data(idx_gr_dp + 46);

    auto grr_z_yy_z = pbuffer.data(idx_gr_dp + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_0, gr_yy_x, gr_yy_y, gr_yy_z, grr_z_yy_x, grr_z_yy_y, grr_z_yy_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yy_x[i] = ts_yy_x[i] * gfe_0 * gc_z[i] + gr_yy_x[i] * gc_z[i];

        grr_z_yy_y[i] = ts_yy_y[i] * gfe_0 * gc_z[i] + gr_yy_y[i] * gc_z[i];

        grr_z_yy_z[i] = ts_yy_0[i] * gfe2_0 + gr_yy_0[i] * gfe_0 + ts_yy_z[i] * gfe_0 * gc_z[i] + gr_yy_z[i] * gc_z[i];
    }

    // Set up 48-51 components of targeted buffer : DP

    auto grr_z_yz_x = pbuffer.data(idx_gr_dp + 48);

    auto grr_z_yz_y = pbuffer.data(idx_gr_dp + 49);

    auto grr_z_yz_z = pbuffer.data(idx_gr_dp + 50);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_x, gr_y_y, gr_y_z, gr_yz_0, gr_yz_x, gr_yz_y, gr_yz_z, grr_z_yz_x, grr_z_yz_y, grr_z_yz_z, ts_y_x, ts_y_y, ts_y_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yz_x[i] = ts_y_x[i] * gfe2_0 + gr_y_x[i] * gfe_0 + ts_yz_x[i] * gfe_0 * gc_z[i] + gr_yz_x[i] * gc_z[i];

        grr_z_yz_y[i] = ts_y_y[i] * gfe2_0 + gr_y_y[i] * gfe_0 + ts_yz_y[i] * gfe_0 * gc_z[i] + gr_yz_y[i] * gc_z[i];

        grr_z_yz_z[i] = ts_y_z[i] * gfe2_0 + gr_y_z[i] * gfe_0 + ts_yz_0[i] * gfe2_0 + gr_yz_0[i] * gfe_0 + ts_yz_z[i] * gfe_0 * gc_z[i] + gr_yz_z[i] * gc_z[i];
    }

    // Set up 51-54 components of targeted buffer : DP

    auto grr_z_zz_x = pbuffer.data(idx_gr_dp + 51);

    auto grr_z_zz_y = pbuffer.data(idx_gr_dp + 52);

    auto grr_z_zz_z = pbuffer.data(idx_gr_dp + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_x, gr_z_y, gr_z_z, gr_zz_0, gr_zz_x, gr_zz_y, gr_zz_z, grr_z_zz_x, grr_z_zz_y, grr_z_zz_z, ts_z_x, ts_z_y, ts_z_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zz_x[i] = 2.0 * ts_z_x[i] * gfe2_0 + 2.0 * gr_z_x[i] * gfe_0 + ts_zz_x[i] * gfe_0 * gc_z[i] + gr_zz_x[i] * gc_z[i];

        grr_z_zz_y[i] = 2.0 * ts_z_y[i] * gfe2_0 + 2.0 * gr_z_y[i] * gfe_0 + ts_zz_y[i] * gfe_0 * gc_z[i] + gr_zz_y[i] * gc_z[i];

        grr_z_zz_z[i] = 2.0 * ts_z_z[i] * gfe2_0 + 2.0 * gr_z_z[i] * gfe_0 + ts_zz_0[i] * gfe2_0 + gr_zz_0[i] * gfe_0 + ts_zz_z[i] * gfe_0 * gc_z[i] + gr_zz_z[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

