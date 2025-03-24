#include "ThreeCenterOverlapPrimRecDP.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_dp(CSimdArray<double>& pbuffer, 
                     const size_t idx_dp,
                     const size_t idx_sp,
                     const size_t idx_ps,
                     const size_t idx_pp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

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

    // Set up 0-3 components of targeted buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    #pragma omp simd aligned(ga_x, ts_0_x, ts_0_y, ts_0_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, ts_xx_x, ts_xx_y, ts_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xx_x[i] = ts_0_x[i] * gfe_0 + ts_x_0[i] * gfe_0 + ts_x_x[i] * ga_x[i];

        ts_xx_y[i] = ts_0_y[i] * gfe_0 + ts_x_y[i] * ga_x[i];

        ts_xx_z[i] = ts_0_z[i] * gfe_0 + ts_x_z[i] * ga_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto ts_xy_x = pbuffer.data(idx_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_dp + 5);

    #pragma omp simd aligned(ga_x, ga_y, ts_x_x, ts_xy_x, ts_xy_y, ts_xy_z, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xy_x[i] = ts_x_x[i] * ga_y[i];

        ts_xy_y[i] = ts_y_y[i] * ga_x[i];

        ts_xy_z[i] = ts_y_z[i] * ga_x[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto ts_xz_x = pbuffer.data(idx_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    #pragma omp simd aligned(ga_x, ga_z, ts_x_x, ts_xz_x, ts_xz_y, ts_xz_z, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xz_x[i] = ts_x_x[i] * ga_z[i];

        ts_xz_y[i] = ts_z_y[i] * ga_x[i];

        ts_xz_z[i] = ts_z_z[i] * ga_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    #pragma omp simd aligned(ga_y, ts_0_x, ts_0_y, ts_0_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yy_x[i] = ts_0_x[i] * gfe_0 + ts_y_x[i] * ga_y[i];

        ts_yy_y[i] = ts_0_y[i] * gfe_0 + ts_y_0[i] * gfe_0 + ts_y_y[i] * ga_y[i];

        ts_yy_z[i] = ts_0_z[i] * gfe_0 + ts_y_z[i] * ga_y[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto ts_yz_x = pbuffer.data(idx_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    #pragma omp simd aligned(ga_y, ga_z, ts_y_y, ts_yz_x, ts_yz_y, ts_yz_z, ts_z_x, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_yz_x[i] = ts_z_x[i] * ga_y[i];

        ts_yz_y[i] = ts_y_y[i] * ga_z[i];

        ts_yz_z[i] = ts_z_z[i] * ga_y[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

    #pragma omp simd aligned(ga_z, ts_0_x, ts_0_y, ts_0_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zz_x[i] = ts_0_x[i] * gfe_0 + ts_z_x[i] * ga_z[i];

        ts_zz_y[i] = ts_0_y[i] * gfe_0 + ts_z_y[i] * ga_z[i];

        ts_zz_z[i] = ts_0_z[i] * gfe_0 + ts_z_0[i] * gfe_0 + ts_z_z[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

