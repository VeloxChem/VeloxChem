#include "NuclearPotentialPrimRecFD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_fd(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_fd,
                               const size_t              idx_npot_0_pd,
                               const size_t              idx_npot_1_pd,
                               const size_t              idx_npot_0_dp,
                               const size_t              idx_npot_1_dp,
                               const size_t              idx_npot_0_dd,
                               const size_t              idx_npot_1_dd,
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

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_0 = pbuffer.data(idx_npot_0_dp);

    auto ta_xx_y_0 = pbuffer.data(idx_npot_0_dp + 1);

    auto ta_xx_z_0 = pbuffer.data(idx_npot_0_dp + 2);

    auto ta_yy_x_0 = pbuffer.data(idx_npot_0_dp + 9);

    auto ta_yy_y_0 = pbuffer.data(idx_npot_0_dp + 10);

    auto ta_yy_z_0 = pbuffer.data(idx_npot_0_dp + 11);

    auto ta_zz_x_0 = pbuffer.data(idx_npot_0_dp + 15);

    auto ta_zz_y_0 = pbuffer.data(idx_npot_0_dp + 16);

    auto ta_zz_z_0 = pbuffer.data(idx_npot_0_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_1 = pbuffer.data(idx_npot_1_dp);

    auto ta_xx_y_1 = pbuffer.data(idx_npot_1_dp + 1);

    auto ta_xx_z_1 = pbuffer.data(idx_npot_1_dp + 2);

    auto ta_yy_x_1 = pbuffer.data(idx_npot_1_dp + 9);

    auto ta_yy_y_1 = pbuffer.data(idx_npot_1_dp + 10);

    auto ta_yy_z_1 = pbuffer.data(idx_npot_1_dp + 11);

    auto ta_zz_x_1 = pbuffer.data(idx_npot_1_dp + 15);

    auto ta_zz_y_1 = pbuffer.data(idx_npot_1_dp + 16);

    auto ta_zz_z_1 = pbuffer.data(idx_npot_1_dp + 17);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_0 = pbuffer.data(idx_npot_0_dd);

    auto ta_xx_xy_0 = pbuffer.data(idx_npot_0_dd + 1);

    auto ta_xx_xz_0 = pbuffer.data(idx_npot_0_dd + 2);

    auto ta_xx_yy_0 = pbuffer.data(idx_npot_0_dd + 3);

    auto ta_xx_yz_0 = pbuffer.data(idx_npot_0_dd + 4);

    auto ta_xx_zz_0 = pbuffer.data(idx_npot_0_dd + 5);

    auto ta_xy_xy_0 = pbuffer.data(idx_npot_0_dd + 7);

    auto ta_xy_yy_0 = pbuffer.data(idx_npot_0_dd + 9);

    auto ta_xy_yz_0 = pbuffer.data(idx_npot_0_dd + 10);

    auto ta_xz_xx_0 = pbuffer.data(idx_npot_0_dd + 12);

    auto ta_xz_xz_0 = pbuffer.data(idx_npot_0_dd + 14);

    auto ta_xz_yz_0 = pbuffer.data(idx_npot_0_dd + 16);

    auto ta_xz_zz_0 = pbuffer.data(idx_npot_0_dd + 17);

    auto ta_yy_xx_0 = pbuffer.data(idx_npot_0_dd + 18);

    auto ta_yy_xy_0 = pbuffer.data(idx_npot_0_dd + 19);

    auto ta_yy_xz_0 = pbuffer.data(idx_npot_0_dd + 20);

    auto ta_yy_yy_0 = pbuffer.data(idx_npot_0_dd + 21);

    auto ta_yy_yz_0 = pbuffer.data(idx_npot_0_dd + 22);

    auto ta_yy_zz_0 = pbuffer.data(idx_npot_0_dd + 23);

    auto ta_yz_xz_0 = pbuffer.data(idx_npot_0_dd + 26);

    auto ta_yz_yy_0 = pbuffer.data(idx_npot_0_dd + 27);

    auto ta_yz_yz_0 = pbuffer.data(idx_npot_0_dd + 28);

    auto ta_yz_zz_0 = pbuffer.data(idx_npot_0_dd + 29);

    auto ta_zz_xx_0 = pbuffer.data(idx_npot_0_dd + 30);

    auto ta_zz_xy_0 = pbuffer.data(idx_npot_0_dd + 31);

    auto ta_zz_xz_0 = pbuffer.data(idx_npot_0_dd + 32);

    auto ta_zz_yy_0 = pbuffer.data(idx_npot_0_dd + 33);

    auto ta_zz_yz_0 = pbuffer.data(idx_npot_0_dd + 34);

    auto ta_zz_zz_0 = pbuffer.data(idx_npot_0_dd + 35);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_1 = pbuffer.data(idx_npot_1_dd);

    auto ta_xx_xy_1 = pbuffer.data(idx_npot_1_dd + 1);

    auto ta_xx_xz_1 = pbuffer.data(idx_npot_1_dd + 2);

    auto ta_xx_yy_1 = pbuffer.data(idx_npot_1_dd + 3);

    auto ta_xx_yz_1 = pbuffer.data(idx_npot_1_dd + 4);

    auto ta_xx_zz_1 = pbuffer.data(idx_npot_1_dd + 5);

    auto ta_xy_xy_1 = pbuffer.data(idx_npot_1_dd + 7);

    auto ta_xy_yy_1 = pbuffer.data(idx_npot_1_dd + 9);

    auto ta_xy_yz_1 = pbuffer.data(idx_npot_1_dd + 10);

    auto ta_xz_xx_1 = pbuffer.data(idx_npot_1_dd + 12);

    auto ta_xz_xz_1 = pbuffer.data(idx_npot_1_dd + 14);

    auto ta_xz_yz_1 = pbuffer.data(idx_npot_1_dd + 16);

    auto ta_xz_zz_1 = pbuffer.data(idx_npot_1_dd + 17);

    auto ta_yy_xx_1 = pbuffer.data(idx_npot_1_dd + 18);

    auto ta_yy_xy_1 = pbuffer.data(idx_npot_1_dd + 19);

    auto ta_yy_xz_1 = pbuffer.data(idx_npot_1_dd + 20);

    auto ta_yy_yy_1 = pbuffer.data(idx_npot_1_dd + 21);

    auto ta_yy_yz_1 = pbuffer.data(idx_npot_1_dd + 22);

    auto ta_yy_zz_1 = pbuffer.data(idx_npot_1_dd + 23);

    auto ta_yz_xz_1 = pbuffer.data(idx_npot_1_dd + 26);

    auto ta_yz_yy_1 = pbuffer.data(idx_npot_1_dd + 27);

    auto ta_yz_yz_1 = pbuffer.data(idx_npot_1_dd + 28);

    auto ta_yz_zz_1 = pbuffer.data(idx_npot_1_dd + 29);

    auto ta_zz_xx_1 = pbuffer.data(idx_npot_1_dd + 30);

    auto ta_zz_xy_1 = pbuffer.data(idx_npot_1_dd + 31);

    auto ta_zz_xz_1 = pbuffer.data(idx_npot_1_dd + 32);

    auto ta_zz_yy_1 = pbuffer.data(idx_npot_1_dd + 33);

    auto ta_zz_yz_1 = pbuffer.data(idx_npot_1_dd + 34);

    auto ta_zz_zz_1 = pbuffer.data(idx_npot_1_dd + 35);

    // Set up 0-6 components of targeted buffer : FD

    auto ta_xxx_xx_0 = pbuffer.data(idx_npot_0_fd);

    auto ta_xxx_xy_0 = pbuffer.data(idx_npot_0_fd + 1);

    auto ta_xxx_xz_0 = pbuffer.data(idx_npot_0_fd + 2);

    auto ta_xxx_yy_0 = pbuffer.data(idx_npot_0_fd + 3);

    auto ta_xxx_yz_0 = pbuffer.data(idx_npot_0_fd + 4);

    auto ta_xxx_zz_0 = pbuffer.data(idx_npot_0_fd + 5);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_x_xx_0,   \
                             ta_x_xx_1,   \
                             ta_x_xy_0,   \
                             ta_x_xy_1,   \
                             ta_x_xz_0,   \
                             ta_x_xz_1,   \
                             ta_x_yy_0,   \
                             ta_x_yy_1,   \
                             ta_x_yz_0,   \
                             ta_x_yz_1,   \
                             ta_x_zz_0,   \
                             ta_x_zz_1,   \
                             ta_xx_x_0,   \
                             ta_xx_x_1,   \
                             ta_xx_xx_0,  \
                             ta_xx_xx_1,  \
                             ta_xx_xy_0,  \
                             ta_xx_xy_1,  \
                             ta_xx_xz_0,  \
                             ta_xx_xz_1,  \
                             ta_xx_y_0,   \
                             ta_xx_y_1,   \
                             ta_xx_yy_0,  \
                             ta_xx_yy_1,  \
                             ta_xx_yz_0,  \
                             ta_xx_yz_1,  \
                             ta_xx_z_0,   \
                             ta_xx_z_1,   \
                             ta_xx_zz_0,  \
                             ta_xx_zz_1,  \
                             ta_xxx_xx_0, \
                             ta_xxx_xy_0, \
                             ta_xxx_xz_0, \
                             ta_xxx_yy_0, \
                             ta_xxx_yz_0, \
                             ta_xxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_xx_0[i] = 2.0 * ta_x_xx_0[i] * fe_0 - 2.0 * ta_x_xx_1[i] * fe_0 + 2.0 * ta_xx_x_0[i] * fe_0 -
                         2.0 * ta_xx_x_1[i] * fe_0 + ta_xx_xx_0[i] * pa_x[i] - ta_xx_xx_1[i] * pc_x[i];

        ta_xxx_xy_0[i] = 2.0 * ta_x_xy_0[i] * fe_0 - 2.0 * ta_x_xy_1[i] * fe_0 + ta_xx_y_0[i] * fe_0 -
                         ta_xx_y_1[i] * fe_0 + ta_xx_xy_0[i] * pa_x[i] - ta_xx_xy_1[i] * pc_x[i];

        ta_xxx_xz_0[i] = 2.0 * ta_x_xz_0[i] * fe_0 - 2.0 * ta_x_xz_1[i] * fe_0 + ta_xx_z_0[i] * fe_0 -
                         ta_xx_z_1[i] * fe_0 + ta_xx_xz_0[i] * pa_x[i] - ta_xx_xz_1[i] * pc_x[i];

        ta_xxx_yy_0[i] =
            2.0 * ta_x_yy_0[i] * fe_0 - 2.0 * ta_x_yy_1[i] * fe_0 + ta_xx_yy_0[i] * pa_x[i] - ta_xx_yy_1[i] * pc_x[i];

        ta_xxx_yz_0[i] =
            2.0 * ta_x_yz_0[i] * fe_0 - 2.0 * ta_x_yz_1[i] * fe_0 + ta_xx_yz_0[i] * pa_x[i] - ta_xx_yz_1[i] * pc_x[i];

        ta_xxx_zz_0[i] =
            2.0 * ta_x_zz_0[i] * fe_0 - 2.0 * ta_x_zz_1[i] * fe_0 + ta_xx_zz_0[i] * pa_x[i] - ta_xx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto ta_xxy_xx_0 = pbuffer.data(idx_npot_0_fd + 6);

    auto ta_xxy_xy_0 = pbuffer.data(idx_npot_0_fd + 7);

    auto ta_xxy_xz_0 = pbuffer.data(idx_npot_0_fd + 8);

    auto ta_xxy_yy_0 = pbuffer.data(idx_npot_0_fd + 9);

    auto ta_xxy_yz_0 = pbuffer.data(idx_npot_0_fd + 10);

    auto ta_xxy_zz_0 = pbuffer.data(idx_npot_0_fd + 11);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pc_x,        \
                             pc_y,        \
                             ta_xx_x_0,   \
                             ta_xx_x_1,   \
                             ta_xx_xx_0,  \
                             ta_xx_xx_1,  \
                             ta_xx_xy_0,  \
                             ta_xx_xy_1,  \
                             ta_xx_xz_0,  \
                             ta_xx_xz_1,  \
                             ta_xx_zz_0,  \
                             ta_xx_zz_1,  \
                             ta_xxy_xx_0, \
                             ta_xxy_xy_0, \
                             ta_xxy_xz_0, \
                             ta_xxy_yy_0, \
                             ta_xxy_yz_0, \
                             ta_xxy_zz_0, \
                             ta_xy_yy_0,  \
                             ta_xy_yy_1,  \
                             ta_xy_yz_0,  \
                             ta_xy_yz_1,  \
                             ta_y_yy_0,   \
                             ta_y_yy_1,   \
                             ta_y_yz_0,   \
                             ta_y_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxy_xx_0[i] = ta_xx_xx_0[i] * pa_y[i] - ta_xx_xx_1[i] * pc_y[i];

        ta_xxy_xy_0[i] = ta_xx_x_0[i] * fe_0 - ta_xx_x_1[i] * fe_0 + ta_xx_xy_0[i] * pa_y[i] - ta_xx_xy_1[i] * pc_y[i];

        ta_xxy_xz_0[i] = ta_xx_xz_0[i] * pa_y[i] - ta_xx_xz_1[i] * pc_y[i];

        ta_xxy_yy_0[i] = ta_y_yy_0[i] * fe_0 - ta_y_yy_1[i] * fe_0 + ta_xy_yy_0[i] * pa_x[i] - ta_xy_yy_1[i] * pc_x[i];

        ta_xxy_yz_0[i] = ta_y_yz_0[i] * fe_0 - ta_y_yz_1[i] * fe_0 + ta_xy_yz_0[i] * pa_x[i] - ta_xy_yz_1[i] * pc_x[i];

        ta_xxy_zz_0[i] = ta_xx_zz_0[i] * pa_y[i] - ta_xx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto ta_xxz_xx_0 = pbuffer.data(idx_npot_0_fd + 12);

    auto ta_xxz_xy_0 = pbuffer.data(idx_npot_0_fd + 13);

    auto ta_xxz_xz_0 = pbuffer.data(idx_npot_0_fd + 14);

    auto ta_xxz_yy_0 = pbuffer.data(idx_npot_0_fd + 15);

    auto ta_xxz_yz_0 = pbuffer.data(idx_npot_0_fd + 16);

    auto ta_xxz_zz_0 = pbuffer.data(idx_npot_0_fd + 17);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             pc_x,        \
                             pc_z,        \
                             ta_xx_x_0,   \
                             ta_xx_x_1,   \
                             ta_xx_xx_0,  \
                             ta_xx_xx_1,  \
                             ta_xx_xy_0,  \
                             ta_xx_xy_1,  \
                             ta_xx_xz_0,  \
                             ta_xx_xz_1,  \
                             ta_xx_yy_0,  \
                             ta_xx_yy_1,  \
                             ta_xxz_xx_0, \
                             ta_xxz_xy_0, \
                             ta_xxz_xz_0, \
                             ta_xxz_yy_0, \
                             ta_xxz_yz_0, \
                             ta_xxz_zz_0, \
                             ta_xz_yz_0,  \
                             ta_xz_yz_1,  \
                             ta_xz_zz_0,  \
                             ta_xz_zz_1,  \
                             ta_z_yz_0,   \
                             ta_z_yz_1,   \
                             ta_z_zz_0,   \
                             ta_z_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxz_xx_0[i] = ta_xx_xx_0[i] * pa_z[i] - ta_xx_xx_1[i] * pc_z[i];

        ta_xxz_xy_0[i] = ta_xx_xy_0[i] * pa_z[i] - ta_xx_xy_1[i] * pc_z[i];

        ta_xxz_xz_0[i] = ta_xx_x_0[i] * fe_0 - ta_xx_x_1[i] * fe_0 + ta_xx_xz_0[i] * pa_z[i] - ta_xx_xz_1[i] * pc_z[i];

        ta_xxz_yy_0[i] = ta_xx_yy_0[i] * pa_z[i] - ta_xx_yy_1[i] * pc_z[i];

        ta_xxz_yz_0[i] = ta_z_yz_0[i] * fe_0 - ta_z_yz_1[i] * fe_0 + ta_xz_yz_0[i] * pa_x[i] - ta_xz_yz_1[i] * pc_x[i];

        ta_xxz_zz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_xz_zz_0[i] * pa_x[i] - ta_xz_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto ta_xyy_xx_0 = pbuffer.data(idx_npot_0_fd + 18);

    auto ta_xyy_xy_0 = pbuffer.data(idx_npot_0_fd + 19);

    auto ta_xyy_xz_0 = pbuffer.data(idx_npot_0_fd + 20);

    auto ta_xyy_yy_0 = pbuffer.data(idx_npot_0_fd + 21);

    auto ta_xyy_yz_0 = pbuffer.data(idx_npot_0_fd + 22);

    auto ta_xyy_zz_0 = pbuffer.data(idx_npot_0_fd + 23);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_xyy_xx_0, \
                             ta_xyy_xy_0, \
                             ta_xyy_xz_0, \
                             ta_xyy_yy_0, \
                             ta_xyy_yz_0, \
                             ta_xyy_zz_0, \
                             ta_yy_x_0,   \
                             ta_yy_x_1,   \
                             ta_yy_xx_0,  \
                             ta_yy_xx_1,  \
                             ta_yy_xy_0,  \
                             ta_yy_xy_1,  \
                             ta_yy_xz_0,  \
                             ta_yy_xz_1,  \
                             ta_yy_y_0,   \
                             ta_yy_y_1,   \
                             ta_yy_yy_0,  \
                             ta_yy_yy_1,  \
                             ta_yy_yz_0,  \
                             ta_yy_yz_1,  \
                             ta_yy_z_0,   \
                             ta_yy_z_1,   \
                             ta_yy_zz_0,  \
                             ta_yy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyy_xx_0[i] =
            2.0 * ta_yy_x_0[i] * fe_0 - 2.0 * ta_yy_x_1[i] * fe_0 + ta_yy_xx_0[i] * pa_x[i] - ta_yy_xx_1[i] * pc_x[i];

        ta_xyy_xy_0[i] = ta_yy_y_0[i] * fe_0 - ta_yy_y_1[i] * fe_0 + ta_yy_xy_0[i] * pa_x[i] - ta_yy_xy_1[i] * pc_x[i];

        ta_xyy_xz_0[i] = ta_yy_z_0[i] * fe_0 - ta_yy_z_1[i] * fe_0 + ta_yy_xz_0[i] * pa_x[i] - ta_yy_xz_1[i] * pc_x[i];

        ta_xyy_yy_0[i] = ta_yy_yy_0[i] * pa_x[i] - ta_yy_yy_1[i] * pc_x[i];

        ta_xyy_yz_0[i] = ta_yy_yz_0[i] * pa_x[i] - ta_yy_yz_1[i] * pc_x[i];

        ta_xyy_zz_0[i] = ta_yy_zz_0[i] * pa_x[i] - ta_yy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto ta_xyz_xx_0 = pbuffer.data(idx_npot_0_fd + 24);

    auto ta_xyz_xy_0 = pbuffer.data(idx_npot_0_fd + 25);

    auto ta_xyz_xz_0 = pbuffer.data(idx_npot_0_fd + 26);

    auto ta_xyz_yy_0 = pbuffer.data(idx_npot_0_fd + 27);

    auto ta_xyz_yz_0 = pbuffer.data(idx_npot_0_fd + 28);

    auto ta_xyz_zz_0 = pbuffer.data(idx_npot_0_fd + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             pc_x,        \
                             pc_y,        \
                             pc_z,        \
                             ta_xy_xy_0,  \
                             ta_xy_xy_1,  \
                             ta_xyz_xx_0, \
                             ta_xyz_xy_0, \
                             ta_xyz_xz_0, \
                             ta_xyz_yy_0, \
                             ta_xyz_yz_0, \
                             ta_xyz_zz_0, \
                             ta_xz_xx_0,  \
                             ta_xz_xx_1,  \
                             ta_xz_xz_0,  \
                             ta_xz_xz_1,  \
                             ta_yz_yy_0,  \
                             ta_yz_yy_1,  \
                             ta_yz_yz_0,  \
                             ta_yz_yz_1,  \
                             ta_yz_zz_0,  \
                             ta_yz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyz_xx_0[i] = ta_xz_xx_0[i] * pa_y[i] - ta_xz_xx_1[i] * pc_y[i];

        ta_xyz_xy_0[i] = ta_xy_xy_0[i] * pa_z[i] - ta_xy_xy_1[i] * pc_z[i];

        ta_xyz_xz_0[i] = ta_xz_xz_0[i] * pa_y[i] - ta_xz_xz_1[i] * pc_y[i];

        ta_xyz_yy_0[i] = ta_yz_yy_0[i] * pa_x[i] - ta_yz_yy_1[i] * pc_x[i];

        ta_xyz_yz_0[i] = ta_yz_yz_0[i] * pa_x[i] - ta_yz_yz_1[i] * pc_x[i];

        ta_xyz_zz_0[i] = ta_yz_zz_0[i] * pa_x[i] - ta_yz_zz_1[i] * pc_x[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto ta_xzz_xx_0 = pbuffer.data(idx_npot_0_fd + 30);

    auto ta_xzz_xy_0 = pbuffer.data(idx_npot_0_fd + 31);

    auto ta_xzz_xz_0 = pbuffer.data(idx_npot_0_fd + 32);

    auto ta_xzz_yy_0 = pbuffer.data(idx_npot_0_fd + 33);

    auto ta_xzz_yz_0 = pbuffer.data(idx_npot_0_fd + 34);

    auto ta_xzz_zz_0 = pbuffer.data(idx_npot_0_fd + 35);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_xzz_xx_0, \
                             ta_xzz_xy_0, \
                             ta_xzz_xz_0, \
                             ta_xzz_yy_0, \
                             ta_xzz_yz_0, \
                             ta_xzz_zz_0, \
                             ta_zz_x_0,   \
                             ta_zz_x_1,   \
                             ta_zz_xx_0,  \
                             ta_zz_xx_1,  \
                             ta_zz_xy_0,  \
                             ta_zz_xy_1,  \
                             ta_zz_xz_0,  \
                             ta_zz_xz_1,  \
                             ta_zz_y_0,   \
                             ta_zz_y_1,   \
                             ta_zz_yy_0,  \
                             ta_zz_yy_1,  \
                             ta_zz_yz_0,  \
                             ta_zz_yz_1,  \
                             ta_zz_z_0,   \
                             ta_zz_z_1,   \
                             ta_zz_zz_0,  \
                             ta_zz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzz_xx_0[i] =
            2.0 * ta_zz_x_0[i] * fe_0 - 2.0 * ta_zz_x_1[i] * fe_0 + ta_zz_xx_0[i] * pa_x[i] - ta_zz_xx_1[i] * pc_x[i];

        ta_xzz_xy_0[i] = ta_zz_y_0[i] * fe_0 - ta_zz_y_1[i] * fe_0 + ta_zz_xy_0[i] * pa_x[i] - ta_zz_xy_1[i] * pc_x[i];

        ta_xzz_xz_0[i] = ta_zz_z_0[i] * fe_0 - ta_zz_z_1[i] * fe_0 + ta_zz_xz_0[i] * pa_x[i] - ta_zz_xz_1[i] * pc_x[i];

        ta_xzz_yy_0[i] = ta_zz_yy_0[i] * pa_x[i] - ta_zz_yy_1[i] * pc_x[i];

        ta_xzz_yz_0[i] = ta_zz_yz_0[i] * pa_x[i] - ta_zz_yz_1[i] * pc_x[i];

        ta_xzz_zz_0[i] = ta_zz_zz_0[i] * pa_x[i] - ta_zz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto ta_yyy_xx_0 = pbuffer.data(idx_npot_0_fd + 36);

    auto ta_yyy_xy_0 = pbuffer.data(idx_npot_0_fd + 37);

    auto ta_yyy_xz_0 = pbuffer.data(idx_npot_0_fd + 38);

    auto ta_yyy_yy_0 = pbuffer.data(idx_npot_0_fd + 39);

    auto ta_yyy_yz_0 = pbuffer.data(idx_npot_0_fd + 40);

    auto ta_yyy_zz_0 = pbuffer.data(idx_npot_0_fd + 41);

#pragma omp simd aligned(pa_y,            \
                             pc_y,        \
                             ta_y_xx_0,   \
                             ta_y_xx_1,   \
                             ta_y_xy_0,   \
                             ta_y_xy_1,   \
                             ta_y_xz_0,   \
                             ta_y_xz_1,   \
                             ta_y_yy_0,   \
                             ta_y_yy_1,   \
                             ta_y_yz_0,   \
                             ta_y_yz_1,   \
                             ta_y_zz_0,   \
                             ta_y_zz_1,   \
                             ta_yy_x_0,   \
                             ta_yy_x_1,   \
                             ta_yy_xx_0,  \
                             ta_yy_xx_1,  \
                             ta_yy_xy_0,  \
                             ta_yy_xy_1,  \
                             ta_yy_xz_0,  \
                             ta_yy_xz_1,  \
                             ta_yy_y_0,   \
                             ta_yy_y_1,   \
                             ta_yy_yy_0,  \
                             ta_yy_yy_1,  \
                             ta_yy_yz_0,  \
                             ta_yy_yz_1,  \
                             ta_yy_z_0,   \
                             ta_yy_z_1,   \
                             ta_yy_zz_0,  \
                             ta_yy_zz_1,  \
                             ta_yyy_xx_0, \
                             ta_yyy_xy_0, \
                             ta_yyy_xz_0, \
                             ta_yyy_yy_0, \
                             ta_yyy_yz_0, \
                             ta_yyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyy_xx_0[i] =
            2.0 * ta_y_xx_0[i] * fe_0 - 2.0 * ta_y_xx_1[i] * fe_0 + ta_yy_xx_0[i] * pa_y[i] - ta_yy_xx_1[i] * pc_y[i];

        ta_yyy_xy_0[i] = 2.0 * ta_y_xy_0[i] * fe_0 - 2.0 * ta_y_xy_1[i] * fe_0 + ta_yy_x_0[i] * fe_0 -
                         ta_yy_x_1[i] * fe_0 + ta_yy_xy_0[i] * pa_y[i] - ta_yy_xy_1[i] * pc_y[i];

        ta_yyy_xz_0[i] =
            2.0 * ta_y_xz_0[i] * fe_0 - 2.0 * ta_y_xz_1[i] * fe_0 + ta_yy_xz_0[i] * pa_y[i] - ta_yy_xz_1[i] * pc_y[i];

        ta_yyy_yy_0[i] = 2.0 * ta_y_yy_0[i] * fe_0 - 2.0 * ta_y_yy_1[i] * fe_0 + 2.0 * ta_yy_y_0[i] * fe_0 -
                         2.0 * ta_yy_y_1[i] * fe_0 + ta_yy_yy_0[i] * pa_y[i] - ta_yy_yy_1[i] * pc_y[i];

        ta_yyy_yz_0[i] = 2.0 * ta_y_yz_0[i] * fe_0 - 2.0 * ta_y_yz_1[i] * fe_0 + ta_yy_z_0[i] * fe_0 -
                         ta_yy_z_1[i] * fe_0 + ta_yy_yz_0[i] * pa_y[i] - ta_yy_yz_1[i] * pc_y[i];

        ta_yyy_zz_0[i] =
            2.0 * ta_y_zz_0[i] * fe_0 - 2.0 * ta_y_zz_1[i] * fe_0 + ta_yy_zz_0[i] * pa_y[i] - ta_yy_zz_1[i] * pc_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto ta_yyz_xx_0 = pbuffer.data(idx_npot_0_fd + 42);

    auto ta_yyz_xy_0 = pbuffer.data(idx_npot_0_fd + 43);

    auto ta_yyz_xz_0 = pbuffer.data(idx_npot_0_fd + 44);

    auto ta_yyz_yy_0 = pbuffer.data(idx_npot_0_fd + 45);

    auto ta_yyz_yz_0 = pbuffer.data(idx_npot_0_fd + 46);

    auto ta_yyz_zz_0 = pbuffer.data(idx_npot_0_fd + 47);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             pc_y,        \
                             pc_z,        \
                             ta_yy_xx_0,  \
                             ta_yy_xx_1,  \
                             ta_yy_xy_0,  \
                             ta_yy_xy_1,  \
                             ta_yy_y_0,   \
                             ta_yy_y_1,   \
                             ta_yy_yy_0,  \
                             ta_yy_yy_1,  \
                             ta_yy_yz_0,  \
                             ta_yy_yz_1,  \
                             ta_yyz_xx_0, \
                             ta_yyz_xy_0, \
                             ta_yyz_xz_0, \
                             ta_yyz_yy_0, \
                             ta_yyz_yz_0, \
                             ta_yyz_zz_0, \
                             ta_yz_xz_0,  \
                             ta_yz_xz_1,  \
                             ta_yz_zz_0,  \
                             ta_yz_zz_1,  \
                             ta_z_xz_0,   \
                             ta_z_xz_1,   \
                             ta_z_zz_0,   \
                             ta_z_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyz_xx_0[i] = ta_yy_xx_0[i] * pa_z[i] - ta_yy_xx_1[i] * pc_z[i];

        ta_yyz_xy_0[i] = ta_yy_xy_0[i] * pa_z[i] - ta_yy_xy_1[i] * pc_z[i];

        ta_yyz_xz_0[i] = ta_z_xz_0[i] * fe_0 - ta_z_xz_1[i] * fe_0 + ta_yz_xz_0[i] * pa_y[i] - ta_yz_xz_1[i] * pc_y[i];

        ta_yyz_yy_0[i] = ta_yy_yy_0[i] * pa_z[i] - ta_yy_yy_1[i] * pc_z[i];

        ta_yyz_yz_0[i] = ta_yy_y_0[i] * fe_0 - ta_yy_y_1[i] * fe_0 + ta_yy_yz_0[i] * pa_z[i] - ta_yy_yz_1[i] * pc_z[i];

        ta_yyz_zz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_yz_zz_0[i] * pa_y[i] - ta_yz_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto ta_yzz_xx_0 = pbuffer.data(idx_npot_0_fd + 48);

    auto ta_yzz_xy_0 = pbuffer.data(idx_npot_0_fd + 49);

    auto ta_yzz_xz_0 = pbuffer.data(idx_npot_0_fd + 50);

    auto ta_yzz_yy_0 = pbuffer.data(idx_npot_0_fd + 51);

    auto ta_yzz_yz_0 = pbuffer.data(idx_npot_0_fd + 52);

    auto ta_yzz_zz_0 = pbuffer.data(idx_npot_0_fd + 53);

#pragma omp simd aligned(pa_y,            \
                             pc_y,        \
                             ta_yzz_xx_0, \
                             ta_yzz_xy_0, \
                             ta_yzz_xz_0, \
                             ta_yzz_yy_0, \
                             ta_yzz_yz_0, \
                             ta_yzz_zz_0, \
                             ta_zz_x_0,   \
                             ta_zz_x_1,   \
                             ta_zz_xx_0,  \
                             ta_zz_xx_1,  \
                             ta_zz_xy_0,  \
                             ta_zz_xy_1,  \
                             ta_zz_xz_0,  \
                             ta_zz_xz_1,  \
                             ta_zz_y_0,   \
                             ta_zz_y_1,   \
                             ta_zz_yy_0,  \
                             ta_zz_yy_1,  \
                             ta_zz_yz_0,  \
                             ta_zz_yz_1,  \
                             ta_zz_z_0,   \
                             ta_zz_z_1,   \
                             ta_zz_zz_0,  \
                             ta_zz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzz_xx_0[i] = ta_zz_xx_0[i] * pa_y[i] - ta_zz_xx_1[i] * pc_y[i];

        ta_yzz_xy_0[i] = ta_zz_x_0[i] * fe_0 - ta_zz_x_1[i] * fe_0 + ta_zz_xy_0[i] * pa_y[i] - ta_zz_xy_1[i] * pc_y[i];

        ta_yzz_xz_0[i] = ta_zz_xz_0[i] * pa_y[i] - ta_zz_xz_1[i] * pc_y[i];

        ta_yzz_yy_0[i] =
            2.0 * ta_zz_y_0[i] * fe_0 - 2.0 * ta_zz_y_1[i] * fe_0 + ta_zz_yy_0[i] * pa_y[i] - ta_zz_yy_1[i] * pc_y[i];

        ta_yzz_yz_0[i] = ta_zz_z_0[i] * fe_0 - ta_zz_z_1[i] * fe_0 + ta_zz_yz_0[i] * pa_y[i] - ta_zz_yz_1[i] * pc_y[i];

        ta_yzz_zz_0[i] = ta_zz_zz_0[i] * pa_y[i] - ta_zz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto ta_zzz_xx_0 = pbuffer.data(idx_npot_0_fd + 54);

    auto ta_zzz_xy_0 = pbuffer.data(idx_npot_0_fd + 55);

    auto ta_zzz_xz_0 = pbuffer.data(idx_npot_0_fd + 56);

    auto ta_zzz_yy_0 = pbuffer.data(idx_npot_0_fd + 57);

    auto ta_zzz_yz_0 = pbuffer.data(idx_npot_0_fd + 58);

    auto ta_zzz_zz_0 = pbuffer.data(idx_npot_0_fd + 59);

#pragma omp simd aligned(pa_z,            \
                             pc_z,        \
                             ta_z_xx_0,   \
                             ta_z_xx_1,   \
                             ta_z_xy_0,   \
                             ta_z_xy_1,   \
                             ta_z_xz_0,   \
                             ta_z_xz_1,   \
                             ta_z_yy_0,   \
                             ta_z_yy_1,   \
                             ta_z_yz_0,   \
                             ta_z_yz_1,   \
                             ta_z_zz_0,   \
                             ta_z_zz_1,   \
                             ta_zz_x_0,   \
                             ta_zz_x_1,   \
                             ta_zz_xx_0,  \
                             ta_zz_xx_1,  \
                             ta_zz_xy_0,  \
                             ta_zz_xy_1,  \
                             ta_zz_xz_0,  \
                             ta_zz_xz_1,  \
                             ta_zz_y_0,   \
                             ta_zz_y_1,   \
                             ta_zz_yy_0,  \
                             ta_zz_yy_1,  \
                             ta_zz_yz_0,  \
                             ta_zz_yz_1,  \
                             ta_zz_z_0,   \
                             ta_zz_z_1,   \
                             ta_zz_zz_0,  \
                             ta_zz_zz_1,  \
                             ta_zzz_xx_0, \
                             ta_zzz_xy_0, \
                             ta_zzz_xz_0, \
                             ta_zzz_yy_0, \
                             ta_zzz_yz_0, \
                             ta_zzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzz_xx_0[i] =
            2.0 * ta_z_xx_0[i] * fe_0 - 2.0 * ta_z_xx_1[i] * fe_0 + ta_zz_xx_0[i] * pa_z[i] - ta_zz_xx_1[i] * pc_z[i];

        ta_zzz_xy_0[i] =
            2.0 * ta_z_xy_0[i] * fe_0 - 2.0 * ta_z_xy_1[i] * fe_0 + ta_zz_xy_0[i] * pa_z[i] - ta_zz_xy_1[i] * pc_z[i];

        ta_zzz_xz_0[i] = 2.0 * ta_z_xz_0[i] * fe_0 - 2.0 * ta_z_xz_1[i] * fe_0 + ta_zz_x_0[i] * fe_0 -
                         ta_zz_x_1[i] * fe_0 + ta_zz_xz_0[i] * pa_z[i] - ta_zz_xz_1[i] * pc_z[i];

        ta_zzz_yy_0[i] =
            2.0 * ta_z_yy_0[i] * fe_0 - 2.0 * ta_z_yy_1[i] * fe_0 + ta_zz_yy_0[i] * pa_z[i] - ta_zz_yy_1[i] * pc_z[i];

        ta_zzz_yz_0[i] = 2.0 * ta_z_yz_0[i] * fe_0 - 2.0 * ta_z_yz_1[i] * fe_0 + ta_zz_y_0[i] * fe_0 -
                         ta_zz_y_1[i] * fe_0 + ta_zz_yz_0[i] * pa_z[i] - ta_zz_yz_1[i] * pc_z[i];

        ta_zzz_zz_0[i] = 2.0 * ta_z_zz_0[i] * fe_0 - 2.0 * ta_z_zz_1[i] * fe_0 + 2.0 * ta_zz_z_0[i] * fe_0 -
                         2.0 * ta_zz_z_1[i] * fe_0 + ta_zz_zz_0[i] * pa_z[i] - ta_zz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
