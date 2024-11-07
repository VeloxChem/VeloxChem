#include "NuclearPotentialGeom020PrimRecPP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_pp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_pp,
                                        const size_t              idx_npot_geom_020_0_ss,
                                        const size_t              idx_npot_geom_020_1_ss,
                                        const size_t              idx_npot_geom_010_1_sp,
                                        const size_t              idx_npot_geom_020_0_sp,
                                        const size_t              idx_npot_geom_020_1_sp,
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

    // Set up components of auxiliary buffer : SS

    auto ta2_xx_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss);

    auto ta2_xy_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 1);

    auto ta2_xz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 2);

    auto ta2_yy_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 3);

    auto ta2_yz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 4);

    auto ta2_zz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 5);

    // Set up components of auxiliary buffer : SS

    auto ta2_xx_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss);

    auto ta2_xy_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 1);

    auto ta2_xz_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 2);

    auto ta2_yy_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 3);

    auto ta2_yz_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 4);

    auto ta2_zz_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 5);

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

    // Set up components of auxiliary buffer : SP

    auto ta2_xx_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp);

    auto ta2_xx_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 1);

    auto ta2_xx_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 2);

    auto ta2_xy_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 3);

    auto ta2_xy_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 4);

    auto ta2_xy_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 5);

    auto ta2_xz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 6);

    auto ta2_xz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 7);

    auto ta2_xz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 8);

    auto ta2_yy_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 9);

    auto ta2_yy_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 10);

    auto ta2_yy_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 11);

    auto ta2_yz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 12);

    auto ta2_yz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 13);

    auto ta2_yz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 14);

    auto ta2_zz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 15);

    auto ta2_zz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 16);

    auto ta2_zz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 17);

    // Set up components of auxiliary buffer : SP

    auto ta2_xx_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp);

    auto ta2_xx_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 1);

    auto ta2_xx_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 2);

    auto ta2_xy_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 3);

    auto ta2_xy_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 4);

    auto ta2_xy_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 5);

    auto ta2_xz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 6);

    auto ta2_xz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 7);

    auto ta2_xz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 8);

    auto ta2_yy_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 9);

    auto ta2_yy_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 10);

    auto ta2_yy_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 11);

    auto ta2_yz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 12);

    auto ta2_yz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 13);

    auto ta2_yz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 14);

    auto ta2_zz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 15);

    auto ta2_zz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 16);

    auto ta2_zz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 17);

    // Set up 0-3 components of targeted buffer : PP

    auto ta2_xx_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp);

    auto ta2_xx_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 1);

    auto ta2_xx_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 2);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_x_0_x_1,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_z_1,  \
                             ta2_xx_0_0_0, \
                             ta2_xx_0_0_1, \
                             ta2_xx_0_x_0, \
                             ta2_xx_0_x_1, \
                             ta2_xx_0_y_0, \
                             ta2_xx_0_y_1, \
                             ta2_xx_0_z_0, \
                             ta2_xx_0_z_1, \
                             ta2_xx_x_x_0, \
                             ta2_xx_x_y_0, \
                             ta2_xx_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_x_x_0[i] =
            ta2_xx_0_0_0[i] * fe_0 - ta2_xx_0_0_1[i] * fe_0 + 2.0 * ta1_x_0_x_1[i] + ta2_xx_0_x_0[i] * pa_x[i] - ta2_xx_0_x_1[i] * pc_x[i];

        ta2_xx_x_y_0[i] = 2.0 * ta1_x_0_y_1[i] + ta2_xx_0_y_0[i] * pa_x[i] - ta2_xx_0_y_1[i] * pc_x[i];

        ta2_xx_x_z_0[i] = 2.0 * ta1_x_0_z_1[i] + ta2_xx_0_z_0[i] * pa_x[i] - ta2_xx_0_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto ta2_xx_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 3);

    auto ta2_xx_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 4);

    auto ta2_xx_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 5);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta2_xx_0_0_0, \
                             ta2_xx_0_0_1, \
                             ta2_xx_0_x_0, \
                             ta2_xx_0_x_1, \
                             ta2_xx_0_y_0, \
                             ta2_xx_0_y_1, \
                             ta2_xx_0_z_0, \
                             ta2_xx_0_z_1, \
                             ta2_xx_y_x_0, \
                             ta2_xx_y_y_0, \
                             ta2_xx_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_y_x_0[i] = ta2_xx_0_x_0[i] * pa_y[i] - ta2_xx_0_x_1[i] * pc_y[i];

        ta2_xx_y_y_0[i] = ta2_xx_0_0_0[i] * fe_0 - ta2_xx_0_0_1[i] * fe_0 + ta2_xx_0_y_0[i] * pa_y[i] - ta2_xx_0_y_1[i] * pc_y[i];

        ta2_xx_y_z_0[i] = ta2_xx_0_z_0[i] * pa_y[i] - ta2_xx_0_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto ta2_xx_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 6);

    auto ta2_xx_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 7);

    auto ta2_xx_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 8);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta2_xx_0_0_0, \
                             ta2_xx_0_0_1, \
                             ta2_xx_0_x_0, \
                             ta2_xx_0_x_1, \
                             ta2_xx_0_y_0, \
                             ta2_xx_0_y_1, \
                             ta2_xx_0_z_0, \
                             ta2_xx_0_z_1, \
                             ta2_xx_z_x_0, \
                             ta2_xx_z_y_0, \
                             ta2_xx_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_z_x_0[i] = ta2_xx_0_x_0[i] * pa_z[i] - ta2_xx_0_x_1[i] * pc_z[i];

        ta2_xx_z_y_0[i] = ta2_xx_0_y_0[i] * pa_z[i] - ta2_xx_0_y_1[i] * pc_z[i];

        ta2_xx_z_z_0[i] = ta2_xx_0_0_0[i] * fe_0 - ta2_xx_0_0_1[i] * fe_0 + ta2_xx_0_z_0[i] * pa_z[i] - ta2_xx_0_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : PP

    auto ta2_xy_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 9);

    auto ta2_xy_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 10);

    auto ta2_xy_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 11);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_y_0_x_1,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_z_1,  \
                             ta2_xy_0_0_0, \
                             ta2_xy_0_0_1, \
                             ta2_xy_0_x_0, \
                             ta2_xy_0_x_1, \
                             ta2_xy_0_y_0, \
                             ta2_xy_0_y_1, \
                             ta2_xy_0_z_0, \
                             ta2_xy_0_z_1, \
                             ta2_xy_x_x_0, \
                             ta2_xy_x_y_0, \
                             ta2_xy_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_x_x_0[i] = ta2_xy_0_0_0[i] * fe_0 - ta2_xy_0_0_1[i] * fe_0 + ta1_y_0_x_1[i] + ta2_xy_0_x_0[i] * pa_x[i] - ta2_xy_0_x_1[i] * pc_x[i];

        ta2_xy_x_y_0[i] = ta1_y_0_y_1[i] + ta2_xy_0_y_0[i] * pa_x[i] - ta2_xy_0_y_1[i] * pc_x[i];

        ta2_xy_x_z_0[i] = ta1_y_0_z_1[i] + ta2_xy_0_z_0[i] * pa_x[i] - ta2_xy_0_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : PP

    auto ta2_xy_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 12);

    auto ta2_xy_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 13);

    auto ta2_xy_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 14);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_x_0_x_1,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_z_1,  \
                             ta2_xy_0_0_0, \
                             ta2_xy_0_0_1, \
                             ta2_xy_0_x_0, \
                             ta2_xy_0_x_1, \
                             ta2_xy_0_y_0, \
                             ta2_xy_0_y_1, \
                             ta2_xy_0_z_0, \
                             ta2_xy_0_z_1, \
                             ta2_xy_y_x_0, \
                             ta2_xy_y_y_0, \
                             ta2_xy_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_y_x_0[i] = ta1_x_0_x_1[i] + ta2_xy_0_x_0[i] * pa_y[i] - ta2_xy_0_x_1[i] * pc_y[i];

        ta2_xy_y_y_0[i] = ta2_xy_0_0_0[i] * fe_0 - ta2_xy_0_0_1[i] * fe_0 + ta1_x_0_y_1[i] + ta2_xy_0_y_0[i] * pa_y[i] - ta2_xy_0_y_1[i] * pc_y[i];

        ta2_xy_y_z_0[i] = ta1_x_0_z_1[i] + ta2_xy_0_z_0[i] * pa_y[i] - ta2_xy_0_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : PP

    auto ta2_xy_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 15);

    auto ta2_xy_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 16);

    auto ta2_xy_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 17);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta2_xy_0_0_0, \
                             ta2_xy_0_0_1, \
                             ta2_xy_0_x_0, \
                             ta2_xy_0_x_1, \
                             ta2_xy_0_y_0, \
                             ta2_xy_0_y_1, \
                             ta2_xy_0_z_0, \
                             ta2_xy_0_z_1, \
                             ta2_xy_z_x_0, \
                             ta2_xy_z_y_0, \
                             ta2_xy_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_z_x_0[i] = ta2_xy_0_x_0[i] * pa_z[i] - ta2_xy_0_x_1[i] * pc_z[i];

        ta2_xy_z_y_0[i] = ta2_xy_0_y_0[i] * pa_z[i] - ta2_xy_0_y_1[i] * pc_z[i];

        ta2_xy_z_z_0[i] = ta2_xy_0_0_0[i] * fe_0 - ta2_xy_0_0_1[i] * fe_0 + ta2_xy_0_z_0[i] * pa_z[i] - ta2_xy_0_z_1[i] * pc_z[i];
    }

    // Set up 18-21 components of targeted buffer : PP

    auto ta2_xz_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 18);

    auto ta2_xz_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 19);

    auto ta2_xz_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 20);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta1_z_0_x_1,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_z_1,  \
                             ta2_xz_0_0_0, \
                             ta2_xz_0_0_1, \
                             ta2_xz_0_x_0, \
                             ta2_xz_0_x_1, \
                             ta2_xz_0_y_0, \
                             ta2_xz_0_y_1, \
                             ta2_xz_0_z_0, \
                             ta2_xz_0_z_1, \
                             ta2_xz_x_x_0, \
                             ta2_xz_x_y_0, \
                             ta2_xz_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_x_x_0[i] = ta2_xz_0_0_0[i] * fe_0 - ta2_xz_0_0_1[i] * fe_0 + ta1_z_0_x_1[i] + ta2_xz_0_x_0[i] * pa_x[i] - ta2_xz_0_x_1[i] * pc_x[i];

        ta2_xz_x_y_0[i] = ta1_z_0_y_1[i] + ta2_xz_0_y_0[i] * pa_x[i] - ta2_xz_0_y_1[i] * pc_x[i];

        ta2_xz_x_z_0[i] = ta1_z_0_z_1[i] + ta2_xz_0_z_0[i] * pa_x[i] - ta2_xz_0_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : PP

    auto ta2_xz_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 21);

    auto ta2_xz_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 22);

    auto ta2_xz_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 23);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta2_xz_0_0_0, \
                             ta2_xz_0_0_1, \
                             ta2_xz_0_x_0, \
                             ta2_xz_0_x_1, \
                             ta2_xz_0_y_0, \
                             ta2_xz_0_y_1, \
                             ta2_xz_0_z_0, \
                             ta2_xz_0_z_1, \
                             ta2_xz_y_x_0, \
                             ta2_xz_y_y_0, \
                             ta2_xz_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_y_x_0[i] = ta2_xz_0_x_0[i] * pa_y[i] - ta2_xz_0_x_1[i] * pc_y[i];

        ta2_xz_y_y_0[i] = ta2_xz_0_0_0[i] * fe_0 - ta2_xz_0_0_1[i] * fe_0 + ta2_xz_0_y_0[i] * pa_y[i] - ta2_xz_0_y_1[i] * pc_y[i];

        ta2_xz_y_z_0[i] = ta2_xz_0_z_0[i] * pa_y[i] - ta2_xz_0_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : PP

    auto ta2_xz_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 24);

    auto ta2_xz_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 25);

    auto ta2_xz_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 26);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_x_0_x_1,  \
                             ta1_x_0_y_1,  \
                             ta1_x_0_z_1,  \
                             ta2_xz_0_0_0, \
                             ta2_xz_0_0_1, \
                             ta2_xz_0_x_0, \
                             ta2_xz_0_x_1, \
                             ta2_xz_0_y_0, \
                             ta2_xz_0_y_1, \
                             ta2_xz_0_z_0, \
                             ta2_xz_0_z_1, \
                             ta2_xz_z_x_0, \
                             ta2_xz_z_y_0, \
                             ta2_xz_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_z_x_0[i] = ta1_x_0_x_1[i] + ta2_xz_0_x_0[i] * pa_z[i] - ta2_xz_0_x_1[i] * pc_z[i];

        ta2_xz_z_y_0[i] = ta1_x_0_y_1[i] + ta2_xz_0_y_0[i] * pa_z[i] - ta2_xz_0_y_1[i] * pc_z[i];

        ta2_xz_z_z_0[i] = ta2_xz_0_0_0[i] * fe_0 - ta2_xz_0_0_1[i] * fe_0 + ta1_x_0_z_1[i] + ta2_xz_0_z_0[i] * pa_z[i] - ta2_xz_0_z_1[i] * pc_z[i];
    }

    // Set up 27-30 components of targeted buffer : PP

    auto ta2_yy_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 27);

    auto ta2_yy_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 28);

    auto ta2_yy_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 29);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta2_yy_0_0_0, \
                             ta2_yy_0_0_1, \
                             ta2_yy_0_x_0, \
                             ta2_yy_0_x_1, \
                             ta2_yy_0_y_0, \
                             ta2_yy_0_y_1, \
                             ta2_yy_0_z_0, \
                             ta2_yy_0_z_1, \
                             ta2_yy_x_x_0, \
                             ta2_yy_x_y_0, \
                             ta2_yy_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_x_x_0[i] = ta2_yy_0_0_0[i] * fe_0 - ta2_yy_0_0_1[i] * fe_0 + ta2_yy_0_x_0[i] * pa_x[i] - ta2_yy_0_x_1[i] * pc_x[i];

        ta2_yy_x_y_0[i] = ta2_yy_0_y_0[i] * pa_x[i] - ta2_yy_0_y_1[i] * pc_x[i];

        ta2_yy_x_z_0[i] = ta2_yy_0_z_0[i] * pa_x[i] - ta2_yy_0_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : PP

    auto ta2_yy_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 30);

    auto ta2_yy_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 31);

    auto ta2_yy_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 32);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_y_0_x_1,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_z_1,  \
                             ta2_yy_0_0_0, \
                             ta2_yy_0_0_1, \
                             ta2_yy_0_x_0, \
                             ta2_yy_0_x_1, \
                             ta2_yy_0_y_0, \
                             ta2_yy_0_y_1, \
                             ta2_yy_0_z_0, \
                             ta2_yy_0_z_1, \
                             ta2_yy_y_x_0, \
                             ta2_yy_y_y_0, \
                             ta2_yy_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_y_x_0[i] = 2.0 * ta1_y_0_x_1[i] + ta2_yy_0_x_0[i] * pa_y[i] - ta2_yy_0_x_1[i] * pc_y[i];

        ta2_yy_y_y_0[i] =
            ta2_yy_0_0_0[i] * fe_0 - ta2_yy_0_0_1[i] * fe_0 + 2.0 * ta1_y_0_y_1[i] + ta2_yy_0_y_0[i] * pa_y[i] - ta2_yy_0_y_1[i] * pc_y[i];

        ta2_yy_y_z_0[i] = 2.0 * ta1_y_0_z_1[i] + ta2_yy_0_z_0[i] * pa_y[i] - ta2_yy_0_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : PP

    auto ta2_yy_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 33);

    auto ta2_yy_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 34);

    auto ta2_yy_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 35);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta2_yy_0_0_0, \
                             ta2_yy_0_0_1, \
                             ta2_yy_0_x_0, \
                             ta2_yy_0_x_1, \
                             ta2_yy_0_y_0, \
                             ta2_yy_0_y_1, \
                             ta2_yy_0_z_0, \
                             ta2_yy_0_z_1, \
                             ta2_yy_z_x_0, \
                             ta2_yy_z_y_0, \
                             ta2_yy_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_z_x_0[i] = ta2_yy_0_x_0[i] * pa_z[i] - ta2_yy_0_x_1[i] * pc_z[i];

        ta2_yy_z_y_0[i] = ta2_yy_0_y_0[i] * pa_z[i] - ta2_yy_0_y_1[i] * pc_z[i];

        ta2_yy_z_z_0[i] = ta2_yy_0_0_0[i] * fe_0 - ta2_yy_0_0_1[i] * fe_0 + ta2_yy_0_z_0[i] * pa_z[i] - ta2_yy_0_z_1[i] * pc_z[i];
    }

    // Set up 36-39 components of targeted buffer : PP

    auto ta2_yz_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 36);

    auto ta2_yz_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 37);

    auto ta2_yz_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 38);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta2_yz_0_0_0, \
                             ta2_yz_0_0_1, \
                             ta2_yz_0_x_0, \
                             ta2_yz_0_x_1, \
                             ta2_yz_0_y_0, \
                             ta2_yz_0_y_1, \
                             ta2_yz_0_z_0, \
                             ta2_yz_0_z_1, \
                             ta2_yz_x_x_0, \
                             ta2_yz_x_y_0, \
                             ta2_yz_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_x_x_0[i] = ta2_yz_0_0_0[i] * fe_0 - ta2_yz_0_0_1[i] * fe_0 + ta2_yz_0_x_0[i] * pa_x[i] - ta2_yz_0_x_1[i] * pc_x[i];

        ta2_yz_x_y_0[i] = ta2_yz_0_y_0[i] * pa_x[i] - ta2_yz_0_y_1[i] * pc_x[i];

        ta2_yz_x_z_0[i] = ta2_yz_0_z_0[i] * pa_x[i] - ta2_yz_0_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : PP

    auto ta2_yz_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 39);

    auto ta2_yz_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 40);

    auto ta2_yz_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 41);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta1_z_0_x_1,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_z_1,  \
                             ta2_yz_0_0_0, \
                             ta2_yz_0_0_1, \
                             ta2_yz_0_x_0, \
                             ta2_yz_0_x_1, \
                             ta2_yz_0_y_0, \
                             ta2_yz_0_y_1, \
                             ta2_yz_0_z_0, \
                             ta2_yz_0_z_1, \
                             ta2_yz_y_x_0, \
                             ta2_yz_y_y_0, \
                             ta2_yz_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_y_x_0[i] = ta1_z_0_x_1[i] + ta2_yz_0_x_0[i] * pa_y[i] - ta2_yz_0_x_1[i] * pc_y[i];

        ta2_yz_y_y_0[i] = ta2_yz_0_0_0[i] * fe_0 - ta2_yz_0_0_1[i] * fe_0 + ta1_z_0_y_1[i] + ta2_yz_0_y_0[i] * pa_y[i] - ta2_yz_0_y_1[i] * pc_y[i];

        ta2_yz_y_z_0[i] = ta1_z_0_z_1[i] + ta2_yz_0_z_0[i] * pa_y[i] - ta2_yz_0_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : PP

    auto ta2_yz_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 42);

    auto ta2_yz_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 43);

    auto ta2_yz_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 44);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_y_0_x_1,  \
                             ta1_y_0_y_1,  \
                             ta1_y_0_z_1,  \
                             ta2_yz_0_0_0, \
                             ta2_yz_0_0_1, \
                             ta2_yz_0_x_0, \
                             ta2_yz_0_x_1, \
                             ta2_yz_0_y_0, \
                             ta2_yz_0_y_1, \
                             ta2_yz_0_z_0, \
                             ta2_yz_0_z_1, \
                             ta2_yz_z_x_0, \
                             ta2_yz_z_y_0, \
                             ta2_yz_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_z_x_0[i] = ta1_y_0_x_1[i] + ta2_yz_0_x_0[i] * pa_z[i] - ta2_yz_0_x_1[i] * pc_z[i];

        ta2_yz_z_y_0[i] = ta1_y_0_y_1[i] + ta2_yz_0_y_0[i] * pa_z[i] - ta2_yz_0_y_1[i] * pc_z[i];

        ta2_yz_z_z_0[i] = ta2_yz_0_0_0[i] * fe_0 - ta2_yz_0_0_1[i] * fe_0 + ta1_y_0_z_1[i] + ta2_yz_0_z_0[i] * pa_z[i] - ta2_yz_0_z_1[i] * pc_z[i];
    }

    // Set up 45-48 components of targeted buffer : PP

    auto ta2_zz_x_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 45);

    auto ta2_zz_x_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 46);

    auto ta2_zz_x_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 47);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta2_zz_0_0_0, \
                             ta2_zz_0_0_1, \
                             ta2_zz_0_x_0, \
                             ta2_zz_0_x_1, \
                             ta2_zz_0_y_0, \
                             ta2_zz_0_y_1, \
                             ta2_zz_0_z_0, \
                             ta2_zz_0_z_1, \
                             ta2_zz_x_x_0, \
                             ta2_zz_x_y_0, \
                             ta2_zz_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_x_x_0[i] = ta2_zz_0_0_0[i] * fe_0 - ta2_zz_0_0_1[i] * fe_0 + ta2_zz_0_x_0[i] * pa_x[i] - ta2_zz_0_x_1[i] * pc_x[i];

        ta2_zz_x_y_0[i] = ta2_zz_0_y_0[i] * pa_x[i] - ta2_zz_0_y_1[i] * pc_x[i];

        ta2_zz_x_z_0[i] = ta2_zz_0_z_0[i] * pa_x[i] - ta2_zz_0_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : PP

    auto ta2_zz_y_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 48);

    auto ta2_zz_y_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 49);

    auto ta2_zz_y_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 50);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta2_zz_0_0_0, \
                             ta2_zz_0_0_1, \
                             ta2_zz_0_x_0, \
                             ta2_zz_0_x_1, \
                             ta2_zz_0_y_0, \
                             ta2_zz_0_y_1, \
                             ta2_zz_0_z_0, \
                             ta2_zz_0_z_1, \
                             ta2_zz_y_x_0, \
                             ta2_zz_y_y_0, \
                             ta2_zz_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_y_x_0[i] = ta2_zz_0_x_0[i] * pa_y[i] - ta2_zz_0_x_1[i] * pc_y[i];

        ta2_zz_y_y_0[i] = ta2_zz_0_0_0[i] * fe_0 - ta2_zz_0_0_1[i] * fe_0 + ta2_zz_0_y_0[i] * pa_y[i] - ta2_zz_0_y_1[i] * pc_y[i];

        ta2_zz_y_z_0[i] = ta2_zz_0_z_0[i] * pa_y[i] - ta2_zz_0_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : PP

    auto ta2_zz_z_x_0 = pbuffer.data(idx_npot_geom_020_0_pp + 51);

    auto ta2_zz_z_y_0 = pbuffer.data(idx_npot_geom_020_0_pp + 52);

    auto ta2_zz_z_z_0 = pbuffer.data(idx_npot_geom_020_0_pp + 53);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta1_z_0_x_1,  \
                             ta1_z_0_y_1,  \
                             ta1_z_0_z_1,  \
                             ta2_zz_0_0_0, \
                             ta2_zz_0_0_1, \
                             ta2_zz_0_x_0, \
                             ta2_zz_0_x_1, \
                             ta2_zz_0_y_0, \
                             ta2_zz_0_y_1, \
                             ta2_zz_0_z_0, \
                             ta2_zz_0_z_1, \
                             ta2_zz_z_x_0, \
                             ta2_zz_z_y_0, \
                             ta2_zz_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_z_x_0[i] = 2.0 * ta1_z_0_x_1[i] + ta2_zz_0_x_0[i] * pa_z[i] - ta2_zz_0_x_1[i] * pc_z[i];

        ta2_zz_z_y_0[i] = 2.0 * ta1_z_0_y_1[i] + ta2_zz_0_y_0[i] * pa_z[i] - ta2_zz_0_y_1[i] * pc_z[i];

        ta2_zz_z_z_0[i] =
            ta2_zz_0_0_0[i] * fe_0 - ta2_zz_0_0_1[i] * fe_0 + 2.0 * ta1_z_0_z_1[i] + ta2_zz_0_z_0[i] * pa_z[i] - ta2_zz_0_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
