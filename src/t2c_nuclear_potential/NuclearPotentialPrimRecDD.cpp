#include "NuclearPotentialPrimRecDD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_dd(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_dd,
                               const size_t              idx_npot_0_sd,
                               const size_t              idx_npot_1_sd,
                               const size_t              idx_npot_0_pp,
                               const size_t              idx_npot_1_pp,
                               const size_t              idx_npot_0_pd,
                               const size_t              idx_npot_1_pd,
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

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_0 = pbuffer.data(idx_npot_0_sd);

    auto ta_0_xy_0 = pbuffer.data(idx_npot_0_sd + 1);

    auto ta_0_xz_0 = pbuffer.data(idx_npot_0_sd + 2);

    auto ta_0_yy_0 = pbuffer.data(idx_npot_0_sd + 3);

    auto ta_0_yz_0 = pbuffer.data(idx_npot_0_sd + 4);

    auto ta_0_zz_0 = pbuffer.data(idx_npot_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = pbuffer.data(idx_npot_1_sd);

    auto ta_0_xy_1 = pbuffer.data(idx_npot_1_sd + 1);

    auto ta_0_xz_1 = pbuffer.data(idx_npot_1_sd + 2);

    auto ta_0_yy_1 = pbuffer.data(idx_npot_1_sd + 3);

    auto ta_0_yz_1 = pbuffer.data(idx_npot_1_sd + 4);

    auto ta_0_zz_1 = pbuffer.data(idx_npot_1_sd + 5);

    // Set up components of auxiliary buffer : PP

    auto ta_x_x_0 = pbuffer.data(idx_npot_0_pp);

    auto ta_x_y_0 = pbuffer.data(idx_npot_0_pp + 1);

    auto ta_x_z_0 = pbuffer.data(idx_npot_0_pp + 2);

    auto ta_y_x_0 = pbuffer.data(idx_npot_0_pp + 3);

    auto ta_y_y_0 = pbuffer.data(idx_npot_0_pp + 4);

    auto ta_y_z_0 = pbuffer.data(idx_npot_0_pp + 5);

    auto ta_z_x_0 = pbuffer.data(idx_npot_0_pp + 6);

    auto ta_z_y_0 = pbuffer.data(idx_npot_0_pp + 7);

    auto ta_z_z_0 = pbuffer.data(idx_npot_0_pp + 8);

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

    // Set up components of auxiliary buffer : PD

    auto ta_x_xx_0 = pbuffer.data(idx_npot_0_pd);

    auto ta_x_xy_0 = pbuffer.data(idx_npot_0_pd + 1);

    auto ta_x_xz_0 = pbuffer.data(idx_npot_0_pd + 2);

    auto ta_x_yy_0 = pbuffer.data(idx_npot_0_pd + 3);

    auto ta_x_yz_0 = pbuffer.data(idx_npot_0_pd + 4);

    auto ta_x_zz_0 = pbuffer.data(idx_npot_0_pd + 5);

    auto ta_y_xx_0 = pbuffer.data(idx_npot_0_pd + 6);

    auto ta_y_xy_0 = pbuffer.data(idx_npot_0_pd + 7);

    auto ta_y_xz_0 = pbuffer.data(idx_npot_0_pd + 8);

    auto ta_y_yy_0 = pbuffer.data(idx_npot_0_pd + 9);

    auto ta_y_yz_0 = pbuffer.data(idx_npot_0_pd + 10);

    auto ta_y_zz_0 = pbuffer.data(idx_npot_0_pd + 11);

    auto ta_z_xx_0 = pbuffer.data(idx_npot_0_pd + 12);

    auto ta_z_xy_0 = pbuffer.data(idx_npot_0_pd + 13);

    auto ta_z_xz_0 = pbuffer.data(idx_npot_0_pd + 14);

    auto ta_z_yy_0 = pbuffer.data(idx_npot_0_pd + 15);

    auto ta_z_yz_0 = pbuffer.data(idx_npot_0_pd + 16);

    auto ta_z_zz_0 = pbuffer.data(idx_npot_0_pd + 17);

    // Set up components of auxiliary buffer : PD

    auto ta_x_xx_1 = pbuffer.data(idx_npot_1_pd);

    auto ta_x_xy_1 = pbuffer.data(idx_npot_1_pd + 1);

    auto ta_x_xz_1 = pbuffer.data(idx_npot_1_pd + 2);

    auto ta_x_yy_1 = pbuffer.data(idx_npot_1_pd + 3);

    auto ta_x_yz_1 = pbuffer.data(idx_npot_1_pd + 4);

    auto ta_x_zz_1 = pbuffer.data(idx_npot_1_pd + 5);

    auto ta_y_xx_1 = pbuffer.data(idx_npot_1_pd + 6);

    auto ta_y_xy_1 = pbuffer.data(idx_npot_1_pd + 7);

    auto ta_y_xz_1 = pbuffer.data(idx_npot_1_pd + 8);

    auto ta_y_yy_1 = pbuffer.data(idx_npot_1_pd + 9);

    auto ta_y_yz_1 = pbuffer.data(idx_npot_1_pd + 10);

    auto ta_y_zz_1 = pbuffer.data(idx_npot_1_pd + 11);

    auto ta_z_xx_1 = pbuffer.data(idx_npot_1_pd + 12);

    auto ta_z_xy_1 = pbuffer.data(idx_npot_1_pd + 13);

    auto ta_z_xz_1 = pbuffer.data(idx_npot_1_pd + 14);

    auto ta_z_yy_1 = pbuffer.data(idx_npot_1_pd + 15);

    auto ta_z_yz_1 = pbuffer.data(idx_npot_1_pd + 16);

    auto ta_z_zz_1 = pbuffer.data(idx_npot_1_pd + 17);

    // Set up 0-6 components of targeted buffer : DD

    auto ta_xx_xx_0 = pbuffer.data(idx_npot_0_dd);

    auto ta_xx_xy_0 = pbuffer.data(idx_npot_0_dd + 1);

    auto ta_xx_xz_0 = pbuffer.data(idx_npot_0_dd + 2);

    auto ta_xx_yy_0 = pbuffer.data(idx_npot_0_dd + 3);

    auto ta_xx_yz_0 = pbuffer.data(idx_npot_0_dd + 4);

    auto ta_xx_zz_0 = pbuffer.data(idx_npot_0_dd + 5);

#pragma omp simd aligned(pa_x,           \
                             pc_x,       \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xy_0,  \
                             ta_0_xy_1,  \
                             ta_0_xz_0,  \
                             ta_0_xz_1,  \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_x_x_0,   \
                             ta_x_x_1,   \
                             ta_x_xx_0,  \
                             ta_x_xx_1,  \
                             ta_x_xy_0,  \
                             ta_x_xy_1,  \
                             ta_x_xz_0,  \
                             ta_x_xz_1,  \
                             ta_x_y_0,   \
                             ta_x_y_1,   \
                             ta_x_yy_0,  \
                             ta_x_yy_1,  \
                             ta_x_yz_0,  \
                             ta_x_yz_1,  \
                             ta_x_z_0,   \
                             ta_x_z_1,   \
                             ta_x_zz_0,  \
                             ta_x_zz_1,  \
                             ta_xx_xx_0, \
                             ta_xx_xy_0, \
                             ta_xx_xz_0, \
                             ta_xx_yy_0, \
                             ta_xx_yz_0, \
                             ta_xx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xx_xx_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + 2.0 * ta_x_x_0[i] * fe_0 -
                        2.0 * ta_x_x_1[i] * fe_0 + ta_x_xx_0[i] * pa_x[i] - ta_x_xx_1[i] * pc_x[i];

        ta_xx_xy_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_x_y_0[i] * fe_0 - ta_x_y_1[i] * fe_0 +
                        ta_x_xy_0[i] * pa_x[i] - ta_x_xy_1[i] * pc_x[i];

        ta_xx_xz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_x_z_0[i] * fe_0 - ta_x_z_1[i] * fe_0 +
                        ta_x_xz_0[i] * pa_x[i] - ta_x_xz_1[i] * pc_x[i];

        ta_xx_yy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_x_yy_0[i] * pa_x[i] - ta_x_yy_1[i] * pc_x[i];

        ta_xx_yz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_x_yz_0[i] * pa_x[i] - ta_x_yz_1[i] * pc_x[i];

        ta_xx_zz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_x_zz_0[i] * pa_x[i] - ta_x_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto ta_xy_xx_0 = pbuffer.data(idx_npot_0_dd + 6);

    auto ta_xy_xy_0 = pbuffer.data(idx_npot_0_dd + 7);

    auto ta_xy_xz_0 = pbuffer.data(idx_npot_0_dd + 8);

    auto ta_xy_yy_0 = pbuffer.data(idx_npot_0_dd + 9);

    auto ta_xy_yz_0 = pbuffer.data(idx_npot_0_dd + 10);

    auto ta_xy_zz_0 = pbuffer.data(idx_npot_0_dd + 11);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pc_x,       \
                             pc_y,       \
                             ta_x_xx_0,  \
                             ta_x_xx_1,  \
                             ta_x_xz_0,  \
                             ta_x_xz_1,  \
                             ta_xy_xx_0, \
                             ta_xy_xy_0, \
                             ta_xy_xz_0, \
                             ta_xy_yy_0, \
                             ta_xy_yz_0, \
                             ta_xy_zz_0, \
                             ta_y_xy_0,  \
                             ta_y_xy_1,  \
                             ta_y_y_0,   \
                             ta_y_y_1,   \
                             ta_y_yy_0,  \
                             ta_y_yy_1,  \
                             ta_y_yz_0,  \
                             ta_y_yz_1,  \
                             ta_y_zz_0,  \
                             ta_y_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xy_xx_0[i] = ta_x_xx_0[i] * pa_y[i] - ta_x_xx_1[i] * pc_y[i];

        ta_xy_xy_0[i] = ta_y_y_0[i] * fe_0 - ta_y_y_1[i] * fe_0 + ta_y_xy_0[i] * pa_x[i] - ta_y_xy_1[i] * pc_x[i];

        ta_xy_xz_0[i] = ta_x_xz_0[i] * pa_y[i] - ta_x_xz_1[i] * pc_y[i];

        ta_xy_yy_0[i] = ta_y_yy_0[i] * pa_x[i] - ta_y_yy_1[i] * pc_x[i];

        ta_xy_yz_0[i] = ta_y_yz_0[i] * pa_x[i] - ta_y_yz_1[i] * pc_x[i];

        ta_xy_zz_0[i] = ta_y_zz_0[i] * pa_x[i] - ta_y_zz_1[i] * pc_x[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto ta_xz_xx_0 = pbuffer.data(idx_npot_0_dd + 12);

    auto ta_xz_xy_0 = pbuffer.data(idx_npot_0_dd + 13);

    auto ta_xz_xz_0 = pbuffer.data(idx_npot_0_dd + 14);

    auto ta_xz_yy_0 = pbuffer.data(idx_npot_0_dd + 15);

    auto ta_xz_yz_0 = pbuffer.data(idx_npot_0_dd + 16);

    auto ta_xz_zz_0 = pbuffer.data(idx_npot_0_dd + 17);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             pc_x,       \
                             pc_z,       \
                             ta_x_xx_0,  \
                             ta_x_xx_1,  \
                             ta_x_xy_0,  \
                             ta_x_xy_1,  \
                             ta_xz_xx_0, \
                             ta_xz_xy_0, \
                             ta_xz_xz_0, \
                             ta_xz_yy_0, \
                             ta_xz_yz_0, \
                             ta_xz_zz_0, \
                             ta_z_xz_0,  \
                             ta_z_xz_1,  \
                             ta_z_yy_0,  \
                             ta_z_yy_1,  \
                             ta_z_yz_0,  \
                             ta_z_yz_1,  \
                             ta_z_z_0,   \
                             ta_z_z_1,   \
                             ta_z_zz_0,  \
                             ta_z_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xz_xx_0[i] = ta_x_xx_0[i] * pa_z[i] - ta_x_xx_1[i] * pc_z[i];

        ta_xz_xy_0[i] = ta_x_xy_0[i] * pa_z[i] - ta_x_xy_1[i] * pc_z[i];

        ta_xz_xz_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_z_xz_0[i] * pa_x[i] - ta_z_xz_1[i] * pc_x[i];

        ta_xz_yy_0[i] = ta_z_yy_0[i] * pa_x[i] - ta_z_yy_1[i] * pc_x[i];

        ta_xz_yz_0[i] = ta_z_yz_0[i] * pa_x[i] - ta_z_yz_1[i] * pc_x[i];

        ta_xz_zz_0[i] = ta_z_zz_0[i] * pa_x[i] - ta_z_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto ta_yy_xx_0 = pbuffer.data(idx_npot_0_dd + 18);

    auto ta_yy_xy_0 = pbuffer.data(idx_npot_0_dd + 19);

    auto ta_yy_xz_0 = pbuffer.data(idx_npot_0_dd + 20);

    auto ta_yy_yy_0 = pbuffer.data(idx_npot_0_dd + 21);

    auto ta_yy_yz_0 = pbuffer.data(idx_npot_0_dd + 22);

    auto ta_yy_zz_0 = pbuffer.data(idx_npot_0_dd + 23);

#pragma omp simd aligned(pa_y,           \
                             pc_y,       \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xy_0,  \
                             ta_0_xy_1,  \
                             ta_0_xz_0,  \
                             ta_0_xz_1,  \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_y_x_0,   \
                             ta_y_x_1,   \
                             ta_y_xx_0,  \
                             ta_y_xx_1,  \
                             ta_y_xy_0,  \
                             ta_y_xy_1,  \
                             ta_y_xz_0,  \
                             ta_y_xz_1,  \
                             ta_y_y_0,   \
                             ta_y_y_1,   \
                             ta_y_yy_0,  \
                             ta_y_yy_1,  \
                             ta_y_yz_0,  \
                             ta_y_yz_1,  \
                             ta_y_z_0,   \
                             ta_y_z_1,   \
                             ta_y_zz_0,  \
                             ta_y_zz_1,  \
                             ta_yy_xx_0, \
                             ta_yy_xy_0, \
                             ta_yy_xz_0, \
                             ta_yy_yy_0, \
                             ta_yy_yz_0, \
                             ta_yy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yy_xx_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_y_xx_0[i] * pa_y[i] - ta_y_xx_1[i] * pc_y[i];

        ta_yy_xy_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_y_x_0[i] * fe_0 - ta_y_x_1[i] * fe_0 +
                        ta_y_xy_0[i] * pa_y[i] - ta_y_xy_1[i] * pc_y[i];

        ta_yy_xz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_y_xz_0[i] * pa_y[i] - ta_y_xz_1[i] * pc_y[i];

        ta_yy_yy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + 2.0 * ta_y_y_0[i] * fe_0 -
                        2.0 * ta_y_y_1[i] * fe_0 + ta_y_yy_0[i] * pa_y[i] - ta_y_yy_1[i] * pc_y[i];

        ta_yy_yz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_y_z_0[i] * fe_0 - ta_y_z_1[i] * fe_0 +
                        ta_y_yz_0[i] * pa_y[i] - ta_y_yz_1[i] * pc_y[i];

        ta_yy_zz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_y_zz_0[i] * pa_y[i] - ta_y_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto ta_yz_xx_0 = pbuffer.data(idx_npot_0_dd + 24);

    auto ta_yz_xy_0 = pbuffer.data(idx_npot_0_dd + 25);

    auto ta_yz_xz_0 = pbuffer.data(idx_npot_0_dd + 26);

    auto ta_yz_yy_0 = pbuffer.data(idx_npot_0_dd + 27);

    auto ta_yz_yz_0 = pbuffer.data(idx_npot_0_dd + 28);

    auto ta_yz_zz_0 = pbuffer.data(idx_npot_0_dd + 29);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             pc_y,       \
                             pc_z,       \
                             ta_y_xy_0,  \
                             ta_y_xy_1,  \
                             ta_y_yy_0,  \
                             ta_y_yy_1,  \
                             ta_yz_xx_0, \
                             ta_yz_xy_0, \
                             ta_yz_xz_0, \
                             ta_yz_yy_0, \
                             ta_yz_yz_0, \
                             ta_yz_zz_0, \
                             ta_z_xx_0,  \
                             ta_z_xx_1,  \
                             ta_z_xz_0,  \
                             ta_z_xz_1,  \
                             ta_z_yz_0,  \
                             ta_z_yz_1,  \
                             ta_z_z_0,   \
                             ta_z_z_1,   \
                             ta_z_zz_0,  \
                             ta_z_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yz_xx_0[i] = ta_z_xx_0[i] * pa_y[i] - ta_z_xx_1[i] * pc_y[i];

        ta_yz_xy_0[i] = ta_y_xy_0[i] * pa_z[i] - ta_y_xy_1[i] * pc_z[i];

        ta_yz_xz_0[i] = ta_z_xz_0[i] * pa_y[i] - ta_z_xz_1[i] * pc_y[i];

        ta_yz_yy_0[i] = ta_y_yy_0[i] * pa_z[i] - ta_y_yy_1[i] * pc_z[i];

        ta_yz_yz_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_z_yz_0[i] * pa_y[i] - ta_z_yz_1[i] * pc_y[i];

        ta_yz_zz_0[i] = ta_z_zz_0[i] * pa_y[i] - ta_z_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto ta_zz_xx_0 = pbuffer.data(idx_npot_0_dd + 30);

    auto ta_zz_xy_0 = pbuffer.data(idx_npot_0_dd + 31);

    auto ta_zz_xz_0 = pbuffer.data(idx_npot_0_dd + 32);

    auto ta_zz_yy_0 = pbuffer.data(idx_npot_0_dd + 33);

    auto ta_zz_yz_0 = pbuffer.data(idx_npot_0_dd + 34);

    auto ta_zz_zz_0 = pbuffer.data(idx_npot_0_dd + 35);

#pragma omp simd aligned(pa_z,           \
                             pc_z,       \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xy_0,  \
                             ta_0_xy_1,  \
                             ta_0_xz_0,  \
                             ta_0_xz_1,  \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_z_x_0,   \
                             ta_z_x_1,   \
                             ta_z_xx_0,  \
                             ta_z_xx_1,  \
                             ta_z_xy_0,  \
                             ta_z_xy_1,  \
                             ta_z_xz_0,  \
                             ta_z_xz_1,  \
                             ta_z_y_0,   \
                             ta_z_y_1,   \
                             ta_z_yy_0,  \
                             ta_z_yy_1,  \
                             ta_z_yz_0,  \
                             ta_z_yz_1,  \
                             ta_z_z_0,   \
                             ta_z_z_1,   \
                             ta_z_zz_0,  \
                             ta_z_zz_1,  \
                             ta_zz_xx_0, \
                             ta_zz_xy_0, \
                             ta_zz_xz_0, \
                             ta_zz_yy_0, \
                             ta_zz_yz_0, \
                             ta_zz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zz_xx_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_z_xx_0[i] * pa_z[i] - ta_z_xx_1[i] * pc_z[i];

        ta_zz_xy_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_z_xy_0[i] * pa_z[i] - ta_z_xy_1[i] * pc_z[i];

        ta_zz_xz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_z_x_0[i] * fe_0 - ta_z_x_1[i] * fe_0 +
                        ta_z_xz_0[i] * pa_z[i] - ta_z_xz_1[i] * pc_z[i];

        ta_zz_yy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_z_yy_0[i] * pa_z[i] - ta_z_yy_1[i] * pc_z[i];

        ta_zz_yz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_z_y_0[i] * fe_0 - ta_z_y_1[i] * fe_0 +
                        ta_z_yz_0[i] * pa_z[i] - ta_z_yz_1[i] * pc_z[i];

        ta_zz_zz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + 2.0 * ta_z_z_0[i] * fe_0 -
                        2.0 * ta_z_z_1[i] * fe_0 + ta_z_zz_0[i] * pa_z[i] - ta_z_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
