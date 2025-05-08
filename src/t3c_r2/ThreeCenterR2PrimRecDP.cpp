#include "ThreeCenterR2PrimRecDP.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_dp(CSimdArray<double>& pbuffer, 
                const size_t idx_g_dp,
                const size_t idx_sp,
                const size_t idx_ps,
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

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

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

    auto gr_xx_x = pbuffer.data(idx_g_dp);

    auto gr_xx_y = pbuffer.data(idx_g_dp + 1);

    auto gr_xx_z = pbuffer.data(idx_g_dp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_x, gr_xx_y, gr_xx_z, ts_0_x, ts_0_y, ts_0_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xx_x[i] = 2.0 * ts_0_x[i] * gfe_0 + 4.0 * ts_x_0[i] * gfe_0 + 4.0 * ts_x_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xx_x[i] * gfe_0 + ts_xx_x[i] * rgc2_0;

        gr_xx_y[i] = 2.0 * ts_0_y[i] * gfe_0 + 4.0 * ts_x_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xx_y[i] * gfe_0 + ts_xx_y[i] * rgc2_0;

        gr_xx_z[i] = 2.0 * ts_0_z[i] * gfe_0 + 4.0 * ts_x_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xx_z[i] * gfe_0 + ts_xx_z[i] * rgc2_0;
    }

    // Set up 3-6 components of targeted buffer : DP

    auto gr_xy_x = pbuffer.data(idx_g_dp + 3);

    auto gr_xy_y = pbuffer.data(idx_g_dp + 4);

    auto gr_xy_z = pbuffer.data(idx_g_dp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_x, gr_xy_y, gr_xy_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xy_x[i] = 2.0 * ts_y_0[i] * gfe_0 + 2.0 * ts_y_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xy_x[i] * gfe_0 + ts_xy_x[i] * rgc2_0;

        gr_xy_y[i] = 2.0 * ts_y_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe_0 + 2.0 * ts_x_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xy_y[i] * gfe_0 + ts_xy_y[i] * rgc2_0;

        gr_xy_z[i] = 2.0 * ts_y_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xy_z[i] * gfe_0 + ts_xy_z[i] * rgc2_0;
    }

    // Set up 6-9 components of targeted buffer : DP

    auto gr_xz_x = pbuffer.data(idx_g_dp + 6);

    auto gr_xz_y = pbuffer.data(idx_g_dp + 7);

    auto gr_xz_z = pbuffer.data(idx_g_dp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_x, gr_xz_y, gr_xz_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xz_x[i] = 2.0 * ts_z_0[i] * gfe_0 + 2.0 * ts_z_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xz_x[i] * gfe_0 + ts_xz_x[i] * rgc2_0;

        gr_xz_y[i] = 2.0 * ts_z_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xz_y[i] * gfe_0 + ts_xz_y[i] * rgc2_0;

        gr_xz_z[i] = 2.0 * ts_z_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe_0 + 2.0 * ts_x_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xz_z[i] * gfe_0 + ts_xz_z[i] * rgc2_0;
    }

    // Set up 9-12 components of targeted buffer : DP

    auto gr_yy_x = pbuffer.data(idx_g_dp + 9);

    auto gr_yy_y = pbuffer.data(idx_g_dp + 10);

    auto gr_yy_z = pbuffer.data(idx_g_dp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_x, gr_yy_y, gr_yy_z, ts_0_x, ts_0_y, ts_0_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yy_x[i] = 2.0 * ts_0_x[i] * gfe_0 + 4.0 * ts_y_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yy_x[i] * gfe_0 + ts_yy_x[i] * rgc2_0;

        gr_yy_y[i] = 2.0 * ts_0_y[i] * gfe_0 + 4.0 * ts_y_0[i] * gfe_0 + 4.0 * ts_y_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yy_y[i] * gfe_0 + ts_yy_y[i] * rgc2_0;

        gr_yy_z[i] = 2.0 * ts_0_z[i] * gfe_0 + 4.0 * ts_y_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yy_z[i] * gfe_0 + ts_yy_z[i] * rgc2_0;
    }

    // Set up 12-15 components of targeted buffer : DP

    auto gr_yz_x = pbuffer.data(idx_g_dp + 12);

    auto gr_yz_y = pbuffer.data(idx_g_dp + 13);

    auto gr_yz_z = pbuffer.data(idx_g_dp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_x, gr_yz_y, gr_yz_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yz_x[i] = 2.0 * ts_z_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yz_x[i] * gfe_0 + ts_yz_x[i] * rgc2_0;

        gr_yz_y[i] = 2.0 * ts_z_0[i] * gfe_0 + 2.0 * ts_z_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yz_y[i] * gfe_0 + ts_yz_y[i] * rgc2_0;

        gr_yz_z[i] = 2.0 * ts_z_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe_0 + 2.0 * ts_y_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_yz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yz_z[i] * gfe_0 + ts_yz_z[i] * rgc2_0;
    }

    // Set up 15-18 components of targeted buffer : DP

    auto gr_zz_x = pbuffer.data(idx_g_dp + 15);

    auto gr_zz_y = pbuffer.data(idx_g_dp + 16);

    auto gr_zz_z = pbuffer.data(idx_g_dp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_x, gr_zz_y, gr_zz_z, ts_0_x, ts_0_y, ts_0_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zz_x[i] = 2.0 * ts_0_x[i] * gfe_0 + 4.0 * ts_z_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_zz_x[i] * gfe_0 + ts_zz_x[i] * rgc2_0;

        gr_zz_y[i] = 2.0 * ts_0_y[i] * gfe_0 + 4.0 * ts_z_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_zz_y[i] * gfe_0 + ts_zz_y[i] * rgc2_0;

        gr_zz_z[i] = 2.0 * ts_0_z[i] * gfe_0 + 4.0 * ts_z_0[i] * gfe_0 + 4.0 * ts_z_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_zz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_zz_z[i] * gfe_0 + ts_zz_z[i] * rgc2_0;
    }

}

} // t3r2rec namespace

