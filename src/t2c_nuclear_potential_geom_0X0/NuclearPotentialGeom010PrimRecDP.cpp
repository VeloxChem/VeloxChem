#include "NuclearPotentialGeom010PrimRecDP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_dp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_dp,
                                        const size_t              idx_npot_geom_010_0_sp,
                                        const size_t              idx_npot_geom_010_1_sp,
                                        const size_t              idx_npot_geom_010_0_ps,
                                        const size_t              idx_npot_geom_010_1_ps,
                                        const size_t              idx_npot_1_pp,
                                        const size_t              idx_npot_geom_010_0_pp,
                                        const size_t              idx_npot_geom_010_1_pp,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SP

    auto ta1_x_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp);

    auto ta1_x_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 1);

    auto ta1_x_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 2);

    auto ta1_y_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp + 3);

    auto ta1_y_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 4);

    auto ta1_y_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 5);

    auto ta1_z_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp + 6);

    auto ta1_z_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 7);

    auto ta1_z_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 8);

    // Set up components of auxiliary buffer : SP

    auto ta1_x_0_x_1 = pbuffer.data(idx_npot_geom_010_1_sp);

    auto ta1_x_0_y_1 = pbuffer.data(idx_npot_geom_010_1_sp + 1);

    auto ta1_x_0_z_1 = pbuffer.data(idx_npot_geom_010_1_sp + 2);

    auto ta1_y_0_x_1 = pbuffer.data(idx_npot_geom_010_1_sp + 3);

    auto ta1_y_0_y_1 = pbuffer.data(idx_npot_geom_010_1_sp + 4);

    auto ta1_y_0_z_1 = pbuffer.data(idx_npot_geom_010_1_sp + 5);

    auto ta1_z_0_x_1 = pbuffer.data(idx_npot_geom_010_1_sp + 6);

    auto ta1_z_0_y_1 = pbuffer.data(idx_npot_geom_010_1_sp + 7);

    auto ta1_z_0_z_1 = pbuffer.data(idx_npot_geom_010_1_sp + 8);

    // Set up components of auxiliary buffer : PS

    auto ta1_x_x_0_0 = pbuffer.data(idx_npot_geom_010_0_ps);

    auto ta1_x_y_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 1);

    auto ta1_x_z_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 2);

    auto ta1_y_x_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 3);

    auto ta1_y_y_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 4);

    auto ta1_y_z_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 5);

    auto ta1_z_x_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 6);

    auto ta1_z_y_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 7);

    auto ta1_z_z_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 8);

    // Set up components of auxiliary buffer : PS

    auto ta1_x_x_0_1 = pbuffer.data(idx_npot_geom_010_1_ps);

    auto ta1_x_y_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 1);

    auto ta1_x_z_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 2);

    auto ta1_y_x_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 3);

    auto ta1_y_y_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 4);

    auto ta1_y_z_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 5);

    auto ta1_z_x_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 6);

    auto ta1_z_y_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 7);

    auto ta1_z_z_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 8);

    // Set up components of auxiliary buffer : PP

    auto ta_x_x_1 = pbuffer.data(idx_npot_1_pp);

    auto ta_x_y_1 = pbuffer.data(idx_npot_1_pp + 1);

    auto ta_x_z_1 = pbuffer.data(idx_npot_1_pp + 2);

    auto ta_y_x_1 = pbuffer.data(idx_npot_1_pp + 3);

    auto ta_y_y_1 = pbuffer.data(idx_npot_1_pp + 4);

    auto ta_y_z_1 = pbuffer.data(idx_npot_1_pp + 5);

    auto ta_z_x_1 = pbuffer.data(idx_npot_1_pp + 6);

    auto ta_z_y_1 = pbuffer.data(idx_npot_1_pp + 7);

    auto ta_z_z_1 = pbuffer.data(idx_npot_1_pp + 8);

    // Set up components of auxiliary buffer : PP

    auto ta1_x_x_x_0 = pbuffer.data(idx_npot_geom_010_0_pp);

    auto ta1_x_x_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 1);

    auto ta1_x_x_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 2);

    auto ta1_x_y_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 3);

    auto ta1_x_y_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 4);

    auto ta1_x_y_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 5);

    auto ta1_x_z_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 6);

    auto ta1_x_z_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 7);

    auto ta1_x_z_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 8);

    auto ta1_y_x_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 9);

    auto ta1_y_x_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 10);

    auto ta1_y_x_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 11);

    auto ta1_y_y_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 12);

    auto ta1_y_y_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 13);

    auto ta1_y_y_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 14);

    auto ta1_y_z_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 15);

    auto ta1_y_z_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 16);

    auto ta1_y_z_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 17);

    auto ta1_z_x_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 18);

    auto ta1_z_x_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 19);

    auto ta1_z_x_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 20);

    auto ta1_z_y_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 21);

    auto ta1_z_y_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 22);

    auto ta1_z_y_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 23);

    auto ta1_z_z_x_0 = pbuffer.data(idx_npot_geom_010_0_pp + 24);

    auto ta1_z_z_y_0 = pbuffer.data(idx_npot_geom_010_0_pp + 25);

    auto ta1_z_z_z_0 = pbuffer.data(idx_npot_geom_010_0_pp + 26);

    // Set up components of auxiliary buffer : PP

    auto ta1_x_x_x_1 = pbuffer.data(idx_npot_geom_010_1_pp);

    auto ta1_x_x_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 1);

    auto ta1_x_x_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 2);

    auto ta1_x_y_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 3);

    auto ta1_x_y_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 4);

    auto ta1_x_y_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 5);

    auto ta1_x_z_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 6);

    auto ta1_x_z_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 7);

    auto ta1_x_z_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 8);

    auto ta1_y_x_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 9);

    auto ta1_y_x_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 10);

    auto ta1_y_x_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 11);

    auto ta1_y_y_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 12);

    auto ta1_y_y_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 13);

    auto ta1_y_y_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 14);

    auto ta1_y_z_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 15);

    auto ta1_y_z_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 16);

    auto ta1_y_z_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 17);

    auto ta1_z_x_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 18);

    auto ta1_z_x_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 19);

    auto ta1_z_x_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 20);

    auto ta1_z_y_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 21);

    auto ta1_z_y_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 22);

    auto ta1_z_y_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 23);

    auto ta1_z_z_x_1 = pbuffer.data(idx_npot_geom_010_1_pp + 24);

    auto ta1_z_z_y_1 = pbuffer.data(idx_npot_geom_010_1_pp + 25);

    auto ta1_z_z_z_1 = pbuffer.data(idx_npot_geom_010_1_pp + 26);

    // Set up 0-3 components of targeted buffer : DP

    auto ta1_x_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp);

    auto ta1_x_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 1);

    auto ta1_x_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 2);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_x_0_x_0,  \
                             ta1_x_0_x_1,  \
                             ta1_x_0_y_0,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_z_0,  \
                             ta1_x_0_z_1,  \
                             ta1_x_x_0_0,  \
                             ta1_x_x_0_1,  \
                             ta1_x_x_x_0,  \
                             ta1_x_x_x_1,  \
                             ta1_x_x_y_0,  \
                             ta1_x_x_y_1,  \
                             ta1_x_x_z_0,  \
                             ta1_x_x_z_1,  \
                             ta1_x_xx_x_0, \
                             ta1_x_xx_y_0, \
                             ta1_x_xx_z_0, \
                             ta_x_x_1,     \
                             ta_x_y_1,     \
                             ta_x_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xx_x_0[i] = ta1_x_0_x_0[i] * fe_0 - ta1_x_0_x_1[i] * fe_0 + ta1_x_x_0_0[i] * fe_0 - ta1_x_x_0_1[i] * fe_0 + ta_x_x_1[i] +
                          ta1_x_x_x_0[i] * pa_x[i] - ta1_x_x_x_1[i] * pc_x[i];

        ta1_x_xx_y_0[i] = ta1_x_0_y_0[i] * fe_0 - ta1_x_0_y_1[i] * fe_0 + ta_x_y_1[i] + ta1_x_x_y_0[i] * pa_x[i] - ta1_x_x_y_1[i] * pc_x[i];

        ta1_x_xx_z_0[i] = ta1_x_0_z_0[i] * fe_0 - ta1_x_0_z_1[i] * fe_0 + ta_x_z_1[i] + ta1_x_x_z_0[i] * pa_x[i] - ta1_x_x_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto ta1_x_xy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 3);

    auto ta1_x_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 4);

    auto ta1_x_xy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 5);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta1_x_x_x_0,  \
                             ta1_x_x_x_1,  \
                             ta1_x_x_z_0,  \
                             ta1_x_x_z_1,  \
                             ta1_x_xy_x_0, \
                             ta1_x_xy_y_0, \
                             ta1_x_xy_z_0, \
                             ta1_x_y_y_0,  \
                             ta1_x_y_y_1,  \
                             ta_y_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xy_x_0[i] = ta1_x_x_x_0[i] * pa_y[i] - ta1_x_x_x_1[i] * pc_y[i];

        ta1_x_xy_y_0[i] = ta_y_y_1[i] + ta1_x_y_y_0[i] * pa_x[i] - ta1_x_y_y_1[i] * pc_x[i];

        ta1_x_xy_z_0[i] = ta1_x_x_z_0[i] * pa_y[i] - ta1_x_x_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto ta1_x_xz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 6);

    auto ta1_x_xz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 7);

    auto ta1_x_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 8);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta1_x_x_x_0,  \
                             ta1_x_x_x_1,  \
                             ta1_x_x_y_0,  \
                             ta1_x_x_y_1,  \
                             ta1_x_xz_x_0, \
                             ta1_x_xz_y_0, \
                             ta1_x_xz_z_0, \
                             ta1_x_z_z_0,  \
                             ta1_x_z_z_1,  \
                             ta_z_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xz_x_0[i] = ta1_x_x_x_0[i] * pa_z[i] - ta1_x_x_x_1[i] * pc_z[i];

        ta1_x_xz_y_0[i] = ta1_x_x_y_0[i] * pa_z[i] - ta1_x_x_y_1[i] * pc_z[i];

        ta1_x_xz_z_0[i] = ta_z_z_1[i] + ta1_x_z_z_0[i] * pa_x[i] - ta1_x_z_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto ta1_x_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 9);

    auto ta1_x_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 10);

    auto ta1_x_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 11);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_x_0_x_0,  \
                             ta1_x_0_x_1,  \
                             ta1_x_0_y_0,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_z_0,  \
                             ta1_x_0_z_1,  \
                             ta1_x_y_0_0,  \
                             ta1_x_y_0_1,  \
                             ta1_x_y_x_0,  \
                             ta1_x_y_x_1,  \
                             ta1_x_y_y_0,  \
                             ta1_x_y_y_1,  \
                             ta1_x_y_z_0,  \
                             ta1_x_y_z_1,  \
                             ta1_x_yy_x_0, \
                             ta1_x_yy_y_0, \
                             ta1_x_yy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yy_x_0[i] = ta1_x_0_x_0[i] * fe_0 - ta1_x_0_x_1[i] * fe_0 + ta1_x_y_x_0[i] * pa_y[i] - ta1_x_y_x_1[i] * pc_y[i];

        ta1_x_yy_y_0[i] = ta1_x_0_y_0[i] * fe_0 - ta1_x_0_y_1[i] * fe_0 + ta1_x_y_0_0[i] * fe_0 - ta1_x_y_0_1[i] * fe_0 + ta1_x_y_y_0[i] * pa_y[i] -
                          ta1_x_y_y_1[i] * pc_y[i];

        ta1_x_yy_z_0[i] = ta1_x_0_z_0[i] * fe_0 - ta1_x_0_z_1[i] * fe_0 + ta1_x_y_z_0[i] * pa_y[i] - ta1_x_y_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto ta1_x_yz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 12);

    auto ta1_x_yz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 13);

    auto ta1_x_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 14);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta1_x_y_y_0,  \
                             ta1_x_y_y_1,  \
                             ta1_x_yz_x_0, \
                             ta1_x_yz_y_0, \
                             ta1_x_yz_z_0, \
                             ta1_x_z_x_0,  \
                             ta1_x_z_x_1,  \
                             ta1_x_z_z_0,  \
                             ta1_x_z_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_yz_x_0[i] = ta1_x_z_x_0[i] * pa_y[i] - ta1_x_z_x_1[i] * pc_y[i];

        ta1_x_yz_y_0[i] = ta1_x_y_y_0[i] * pa_z[i] - ta1_x_y_y_1[i] * pc_z[i];

        ta1_x_yz_z_0[i] = ta1_x_z_z_0[i] * pa_y[i] - ta1_x_z_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto ta1_x_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 15);

    auto ta1_x_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 16);

    auto ta1_x_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 17);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_x_0_x_0,  \
                             ta1_x_0_x_1,  \
                             ta1_x_0_y_0,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_z_0,  \
                             ta1_x_0_z_1,  \
                             ta1_x_z_0_0,  \
                             ta1_x_z_0_1,  \
                             ta1_x_z_x_0,  \
                             ta1_x_z_x_1,  \
                             ta1_x_z_y_0,  \
                             ta1_x_z_y_1,  \
                             ta1_x_z_z_0,  \
                             ta1_x_z_z_1,  \
                             ta1_x_zz_x_0, \
                             ta1_x_zz_y_0, \
                             ta1_x_zz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zz_x_0[i] = ta1_x_0_x_0[i] * fe_0 - ta1_x_0_x_1[i] * fe_0 + ta1_x_z_x_0[i] * pa_z[i] - ta1_x_z_x_1[i] * pc_z[i];

        ta1_x_zz_y_0[i] = ta1_x_0_y_0[i] * fe_0 - ta1_x_0_y_1[i] * fe_0 + ta1_x_z_y_0[i] * pa_z[i] - ta1_x_z_y_1[i] * pc_z[i];

        ta1_x_zz_z_0[i] = ta1_x_0_z_0[i] * fe_0 - ta1_x_0_z_1[i] * fe_0 + ta1_x_z_0_0[i] * fe_0 - ta1_x_z_0_1[i] * fe_0 + ta1_x_z_z_0[i] * pa_z[i] -
                          ta1_x_z_z_1[i] * pc_z[i];
    }

    // Set up 18-21 components of targeted buffer : DP

    auto ta1_y_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 18);

    auto ta1_y_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 19);

    auto ta1_y_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 20);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_y_0_x_0,  \
                             ta1_y_0_x_1,  \
                             ta1_y_0_y_0,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_z_0,  \
                             ta1_y_0_z_1,  \
                             ta1_y_x_0_0,  \
                             ta1_y_x_0_1,  \
                             ta1_y_x_x_0,  \
                             ta1_y_x_x_1,  \
                             ta1_y_x_y_0,  \
                             ta1_y_x_y_1,  \
                             ta1_y_x_z_0,  \
                             ta1_y_x_z_1,  \
                             ta1_y_xx_x_0, \
                             ta1_y_xx_y_0, \
                             ta1_y_xx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xx_x_0[i] = ta1_y_0_x_0[i] * fe_0 - ta1_y_0_x_1[i] * fe_0 + ta1_y_x_0_0[i] * fe_0 - ta1_y_x_0_1[i] * fe_0 + ta1_y_x_x_0[i] * pa_x[i] -
                          ta1_y_x_x_1[i] * pc_x[i];

        ta1_y_xx_y_0[i] = ta1_y_0_y_0[i] * fe_0 - ta1_y_0_y_1[i] * fe_0 + ta1_y_x_y_0[i] * pa_x[i] - ta1_y_x_y_1[i] * pc_x[i];

        ta1_y_xx_z_0[i] = ta1_y_0_z_0[i] * fe_0 - ta1_y_0_z_1[i] * fe_0 + ta1_y_x_z_0[i] * pa_x[i] - ta1_y_x_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : DP

    auto ta1_y_xy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 21);

    auto ta1_y_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 22);

    auto ta1_y_xy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 23);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta1_y_x_x_0,  \
                             ta1_y_x_x_1,  \
                             ta1_y_xy_x_0, \
                             ta1_y_xy_y_0, \
                             ta1_y_xy_z_0, \
                             ta1_y_y_y_0,  \
                             ta1_y_y_y_1,  \
                             ta1_y_y_z_0,  \
                             ta1_y_y_z_1,  \
                             ta_x_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xy_x_0[i] = ta_x_x_1[i] + ta1_y_x_x_0[i] * pa_y[i] - ta1_y_x_x_1[i] * pc_y[i];

        ta1_y_xy_y_0[i] = ta1_y_y_y_0[i] * pa_x[i] - ta1_y_y_y_1[i] * pc_x[i];

        ta1_y_xy_z_0[i] = ta1_y_y_z_0[i] * pa_x[i] - ta1_y_y_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : DP

    auto ta1_y_xz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 24);

    auto ta1_y_xz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 25);

    auto ta1_y_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 26);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta1_y_x_x_0,  \
                             ta1_y_x_x_1,  \
                             ta1_y_xz_x_0, \
                             ta1_y_xz_y_0, \
                             ta1_y_xz_z_0, \
                             ta1_y_z_y_0,  \
                             ta1_y_z_y_1,  \
                             ta1_y_z_z_0,  \
                             ta1_y_z_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xz_x_0[i] = ta1_y_x_x_0[i] * pa_z[i] - ta1_y_x_x_1[i] * pc_z[i];

        ta1_y_xz_y_0[i] = ta1_y_z_y_0[i] * pa_x[i] - ta1_y_z_y_1[i] * pc_x[i];

        ta1_y_xz_z_0[i] = ta1_y_z_z_0[i] * pa_x[i] - ta1_y_z_z_1[i] * pc_x[i];
    }

    // Set up 27-30 components of targeted buffer : DP

    auto ta1_y_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 27);

    auto ta1_y_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 28);

    auto ta1_y_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 29);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_y_0_x_0,  \
                             ta1_y_0_x_1,  \
                             ta1_y_0_y_0,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_z_0,  \
                             ta1_y_0_z_1,  \
                             ta1_y_y_0_0,  \
                             ta1_y_y_0_1,  \
                             ta1_y_y_x_0,  \
                             ta1_y_y_x_1,  \
                             ta1_y_y_y_0,  \
                             ta1_y_y_y_1,  \
                             ta1_y_y_z_0,  \
                             ta1_y_y_z_1,  \
                             ta1_y_yy_x_0, \
                             ta1_y_yy_y_0, \
                             ta1_y_yy_z_0, \
                             ta_y_x_1,     \
                             ta_y_y_1,     \
                             ta_y_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yy_x_0[i] = ta1_y_0_x_0[i] * fe_0 - ta1_y_0_x_1[i] * fe_0 + ta_y_x_1[i] + ta1_y_y_x_0[i] * pa_y[i] - ta1_y_y_x_1[i] * pc_y[i];

        ta1_y_yy_y_0[i] = ta1_y_0_y_0[i] * fe_0 - ta1_y_0_y_1[i] * fe_0 + ta1_y_y_0_0[i] * fe_0 - ta1_y_y_0_1[i] * fe_0 + ta_y_y_1[i] +
                          ta1_y_y_y_0[i] * pa_y[i] - ta1_y_y_y_1[i] * pc_y[i];

        ta1_y_yy_z_0[i] = ta1_y_0_z_0[i] * fe_0 - ta1_y_0_z_1[i] * fe_0 + ta_y_z_1[i] + ta1_y_y_z_0[i] * pa_y[i] - ta1_y_y_z_1[i] * pc_y[i];
    }

    // Set up 30-33 components of targeted buffer : DP

    auto ta1_y_yz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 30);

    auto ta1_y_yz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 31);

    auto ta1_y_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 32);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta1_y_y_x_0,  \
                             ta1_y_y_x_1,  \
                             ta1_y_y_y_0,  \
                             ta1_y_y_y_1,  \
                             ta1_y_yz_x_0, \
                             ta1_y_yz_y_0, \
                             ta1_y_yz_z_0, \
                             ta1_y_z_z_0,  \
                             ta1_y_z_z_1,  \
                             ta_z_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_yz_x_0[i] = ta1_y_y_x_0[i] * pa_z[i] - ta1_y_y_x_1[i] * pc_z[i];

        ta1_y_yz_y_0[i] = ta1_y_y_y_0[i] * pa_z[i] - ta1_y_y_y_1[i] * pc_z[i];

        ta1_y_yz_z_0[i] = ta_z_z_1[i] + ta1_y_z_z_0[i] * pa_y[i] - ta1_y_z_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : DP

    auto ta1_y_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 33);

    auto ta1_y_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 34);

    auto ta1_y_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 35);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_y_0_x_0,  \
                             ta1_y_0_x_1,  \
                             ta1_y_0_y_0,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_z_0,  \
                             ta1_y_0_z_1,  \
                             ta1_y_z_0_0,  \
                             ta1_y_z_0_1,  \
                             ta1_y_z_x_0,  \
                             ta1_y_z_x_1,  \
                             ta1_y_z_y_0,  \
                             ta1_y_z_y_1,  \
                             ta1_y_z_z_0,  \
                             ta1_y_z_z_1,  \
                             ta1_y_zz_x_0, \
                             ta1_y_zz_y_0, \
                             ta1_y_zz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zz_x_0[i] = ta1_y_0_x_0[i] * fe_0 - ta1_y_0_x_1[i] * fe_0 + ta1_y_z_x_0[i] * pa_z[i] - ta1_y_z_x_1[i] * pc_z[i];

        ta1_y_zz_y_0[i] = ta1_y_0_y_0[i] * fe_0 - ta1_y_0_y_1[i] * fe_0 + ta1_y_z_y_0[i] * pa_z[i] - ta1_y_z_y_1[i] * pc_z[i];

        ta1_y_zz_z_0[i] = ta1_y_0_z_0[i] * fe_0 - ta1_y_0_z_1[i] * fe_0 + ta1_y_z_0_0[i] * fe_0 - ta1_y_z_0_1[i] * fe_0 + ta1_y_z_z_0[i] * pa_z[i] -
                          ta1_y_z_z_1[i] * pc_z[i];
    }

    // Set up 36-39 components of targeted buffer : DP

    auto ta1_z_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 36);

    auto ta1_z_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 37);

    auto ta1_z_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 38);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_z_0_x_0,  \
                             ta1_z_0_x_1,  \
                             ta1_z_0_y_0,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_z_0,  \
                             ta1_z_0_z_1,  \
                             ta1_z_x_0_0,  \
                             ta1_z_x_0_1,  \
                             ta1_z_x_x_0,  \
                             ta1_z_x_x_1,  \
                             ta1_z_x_y_0,  \
                             ta1_z_x_y_1,  \
                             ta1_z_x_z_0,  \
                             ta1_z_x_z_1,  \
                             ta1_z_xx_x_0, \
                             ta1_z_xx_y_0, \
                             ta1_z_xx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xx_x_0[i] = ta1_z_0_x_0[i] * fe_0 - ta1_z_0_x_1[i] * fe_0 + ta1_z_x_0_0[i] * fe_0 - ta1_z_x_0_1[i] * fe_0 + ta1_z_x_x_0[i] * pa_x[i] -
                          ta1_z_x_x_1[i] * pc_x[i];

        ta1_z_xx_y_0[i] = ta1_z_0_y_0[i] * fe_0 - ta1_z_0_y_1[i] * fe_0 + ta1_z_x_y_0[i] * pa_x[i] - ta1_z_x_y_1[i] * pc_x[i];

        ta1_z_xx_z_0[i] = ta1_z_0_z_0[i] * fe_0 - ta1_z_0_z_1[i] * fe_0 + ta1_z_x_z_0[i] * pa_x[i] - ta1_z_x_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : DP

    auto ta1_z_xy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 39);

    auto ta1_z_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 40);

    auto ta1_z_xy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 41);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta1_z_x_x_0,  \
                             ta1_z_x_x_1,  \
                             ta1_z_xy_x_0, \
                             ta1_z_xy_y_0, \
                             ta1_z_xy_z_0, \
                             ta1_z_y_y_0,  \
                             ta1_z_y_y_1,  \
                             ta1_z_y_z_0,  \
                             ta1_z_y_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xy_x_0[i] = ta1_z_x_x_0[i] * pa_y[i] - ta1_z_x_x_1[i] * pc_y[i];

        ta1_z_xy_y_0[i] = ta1_z_y_y_0[i] * pa_x[i] - ta1_z_y_y_1[i] * pc_x[i];

        ta1_z_xy_z_0[i] = ta1_z_y_z_0[i] * pa_x[i] - ta1_z_y_z_1[i] * pc_x[i];
    }

    // Set up 42-45 components of targeted buffer : DP

    auto ta1_z_xz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 42);

    auto ta1_z_xz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 43);

    auto ta1_z_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 44);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta1_z_x_x_0,  \
                             ta1_z_x_x_1,  \
                             ta1_z_xz_x_0, \
                             ta1_z_xz_y_0, \
                             ta1_z_xz_z_0, \
                             ta1_z_z_y_0,  \
                             ta1_z_z_y_1,  \
                             ta1_z_z_z_0,  \
                             ta1_z_z_z_1,  \
                             ta_x_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xz_x_0[i] = ta_x_x_1[i] + ta1_z_x_x_0[i] * pa_z[i] - ta1_z_x_x_1[i] * pc_z[i];

        ta1_z_xz_y_0[i] = ta1_z_z_y_0[i] * pa_x[i] - ta1_z_z_y_1[i] * pc_x[i];

        ta1_z_xz_z_0[i] = ta1_z_z_z_0[i] * pa_x[i] - ta1_z_z_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : DP

    auto ta1_z_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 45);

    auto ta1_z_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 46);

    auto ta1_z_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 47);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_z_0_x_0,  \
                             ta1_z_0_x_1,  \
                             ta1_z_0_y_0,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_z_0,  \
                             ta1_z_0_z_1,  \
                             ta1_z_y_0_0,  \
                             ta1_z_y_0_1,  \
                             ta1_z_y_x_0,  \
                             ta1_z_y_x_1,  \
                             ta1_z_y_y_0,  \
                             ta1_z_y_y_1,  \
                             ta1_z_y_z_0,  \
                             ta1_z_y_z_1,  \
                             ta1_z_yy_x_0, \
                             ta1_z_yy_y_0, \
                             ta1_z_yy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yy_x_0[i] = ta1_z_0_x_0[i] * fe_0 - ta1_z_0_x_1[i] * fe_0 + ta1_z_y_x_0[i] * pa_y[i] - ta1_z_y_x_1[i] * pc_y[i];

        ta1_z_yy_y_0[i] = ta1_z_0_y_0[i] * fe_0 - ta1_z_0_y_1[i] * fe_0 + ta1_z_y_0_0[i] * fe_0 - ta1_z_y_0_1[i] * fe_0 + ta1_z_y_y_0[i] * pa_y[i] -
                          ta1_z_y_y_1[i] * pc_y[i];

        ta1_z_yy_z_0[i] = ta1_z_0_z_0[i] * fe_0 - ta1_z_0_z_1[i] * fe_0 + ta1_z_y_z_0[i] * pa_y[i] - ta1_z_y_z_1[i] * pc_y[i];
    }

    // Set up 48-51 components of targeted buffer : DP

    auto ta1_z_yz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 48);

    auto ta1_z_yz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 49);

    auto ta1_z_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 50);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta1_z_y_y_0,  \
                             ta1_z_y_y_1,  \
                             ta1_z_yz_x_0, \
                             ta1_z_yz_y_0, \
                             ta1_z_yz_z_0, \
                             ta1_z_z_x_0,  \
                             ta1_z_z_x_1,  \
                             ta1_z_z_z_0,  \
                             ta1_z_z_z_1,  \
                             ta_y_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_yz_x_0[i] = ta1_z_z_x_0[i] * pa_y[i] - ta1_z_z_x_1[i] * pc_y[i];

        ta1_z_yz_y_0[i] = ta_y_y_1[i] + ta1_z_y_y_0[i] * pa_z[i] - ta1_z_y_y_1[i] * pc_z[i];

        ta1_z_yz_z_0[i] = ta1_z_z_z_0[i] * pa_y[i] - ta1_z_z_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : DP

    auto ta1_z_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 51);

    auto ta1_z_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 52);

    auto ta1_z_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 53);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_z_0_x_0,  \
                             ta1_z_0_x_1,  \
                             ta1_z_0_y_0,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_z_0,  \
                             ta1_z_0_z_1,  \
                             ta1_z_z_0_0,  \
                             ta1_z_z_0_1,  \
                             ta1_z_z_x_0,  \
                             ta1_z_z_x_1,  \
                             ta1_z_z_y_0,  \
                             ta1_z_z_y_1,  \
                             ta1_z_z_z_0,  \
                             ta1_z_z_z_1,  \
                             ta1_z_zz_x_0, \
                             ta1_z_zz_y_0, \
                             ta1_z_zz_z_0, \
                             ta_z_x_1,     \
                             ta_z_y_1,     \
                             ta_z_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zz_x_0[i] = ta1_z_0_x_0[i] * fe_0 - ta1_z_0_x_1[i] * fe_0 + ta_z_x_1[i] + ta1_z_z_x_0[i] * pa_z[i] - ta1_z_z_x_1[i] * pc_z[i];

        ta1_z_zz_y_0[i] = ta1_z_0_y_0[i] * fe_0 - ta1_z_0_y_1[i] * fe_0 + ta_z_y_1[i] + ta1_z_z_y_0[i] * pa_z[i] - ta1_z_z_y_1[i] * pc_z[i];

        ta1_z_zz_z_0[i] = ta1_z_0_z_0[i] * fe_0 - ta1_z_0_z_1[i] * fe_0 + ta1_z_z_0_0[i] * fe_0 - ta1_z_z_0_1[i] * fe_0 + ta_z_z_1[i] +
                          ta1_z_z_z_0[i] * pa_z[i] - ta1_z_z_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
