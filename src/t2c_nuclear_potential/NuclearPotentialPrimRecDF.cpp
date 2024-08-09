#include "NuclearPotentialPrimRecDF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_df(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_df,
                               const size_t              idx_npot_0_sf,
                               const size_t              idx_npot_1_sf,
                               const size_t              idx_npot_0_pd,
                               const size_t              idx_npot_1_pd,
                               const size_t              idx_npot_0_pf,
                               const size_t              idx_npot_1_pf,
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

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_0 = pbuffer.data(idx_npot_0_sf);

    auto ta_0_xxy_0 = pbuffer.data(idx_npot_0_sf + 1);

    auto ta_0_xxz_0 = pbuffer.data(idx_npot_0_sf + 2);

    auto ta_0_xyy_0 = pbuffer.data(idx_npot_0_sf + 3);

    auto ta_0_xyz_0 = pbuffer.data(idx_npot_0_sf + 4);

    auto ta_0_xzz_0 = pbuffer.data(idx_npot_0_sf + 5);

    auto ta_0_yyy_0 = pbuffer.data(idx_npot_0_sf + 6);

    auto ta_0_yyz_0 = pbuffer.data(idx_npot_0_sf + 7);

    auto ta_0_yzz_0 = pbuffer.data(idx_npot_0_sf + 8);

    auto ta_0_zzz_0 = pbuffer.data(idx_npot_0_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_1 = pbuffer.data(idx_npot_1_sf);

    auto ta_0_xxy_1 = pbuffer.data(idx_npot_1_sf + 1);

    auto ta_0_xxz_1 = pbuffer.data(idx_npot_1_sf + 2);

    auto ta_0_xyy_1 = pbuffer.data(idx_npot_1_sf + 3);

    auto ta_0_xyz_1 = pbuffer.data(idx_npot_1_sf + 4);

    auto ta_0_xzz_1 = pbuffer.data(idx_npot_1_sf + 5);

    auto ta_0_yyy_1 = pbuffer.data(idx_npot_1_sf + 6);

    auto ta_0_yyz_1 = pbuffer.data(idx_npot_1_sf + 7);

    auto ta_0_yzz_1 = pbuffer.data(idx_npot_1_sf + 8);

    auto ta_0_zzz_1 = pbuffer.data(idx_npot_1_sf + 9);

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

    // Set up components of auxiliary buffer : PF

    auto ta_x_xxx_0 = pbuffer.data(idx_npot_0_pf);

    auto ta_x_xxy_0 = pbuffer.data(idx_npot_0_pf + 1);

    auto ta_x_xxz_0 = pbuffer.data(idx_npot_0_pf + 2);

    auto ta_x_xyy_0 = pbuffer.data(idx_npot_0_pf + 3);

    auto ta_x_xyz_0 = pbuffer.data(idx_npot_0_pf + 4);

    auto ta_x_xzz_0 = pbuffer.data(idx_npot_0_pf + 5);

    auto ta_x_yyy_0 = pbuffer.data(idx_npot_0_pf + 6);

    auto ta_x_yyz_0 = pbuffer.data(idx_npot_0_pf + 7);

    auto ta_x_yzz_0 = pbuffer.data(idx_npot_0_pf + 8);

    auto ta_x_zzz_0 = pbuffer.data(idx_npot_0_pf + 9);

    auto ta_y_xxx_0 = pbuffer.data(idx_npot_0_pf + 10);

    auto ta_y_xxy_0 = pbuffer.data(idx_npot_0_pf + 11);

    auto ta_y_xxz_0 = pbuffer.data(idx_npot_0_pf + 12);

    auto ta_y_xyy_0 = pbuffer.data(idx_npot_0_pf + 13);

    auto ta_y_xyz_0 = pbuffer.data(idx_npot_0_pf + 14);

    auto ta_y_xzz_0 = pbuffer.data(idx_npot_0_pf + 15);

    auto ta_y_yyy_0 = pbuffer.data(idx_npot_0_pf + 16);

    auto ta_y_yyz_0 = pbuffer.data(idx_npot_0_pf + 17);

    auto ta_y_yzz_0 = pbuffer.data(idx_npot_0_pf + 18);

    auto ta_y_zzz_0 = pbuffer.data(idx_npot_0_pf + 19);

    auto ta_z_xxx_0 = pbuffer.data(idx_npot_0_pf + 20);

    auto ta_z_xxy_0 = pbuffer.data(idx_npot_0_pf + 21);

    auto ta_z_xxz_0 = pbuffer.data(idx_npot_0_pf + 22);

    auto ta_z_xyy_0 = pbuffer.data(idx_npot_0_pf + 23);

    auto ta_z_xyz_0 = pbuffer.data(idx_npot_0_pf + 24);

    auto ta_z_xzz_0 = pbuffer.data(idx_npot_0_pf + 25);

    auto ta_z_yyy_0 = pbuffer.data(idx_npot_0_pf + 26);

    auto ta_z_yyz_0 = pbuffer.data(idx_npot_0_pf + 27);

    auto ta_z_yzz_0 = pbuffer.data(idx_npot_0_pf + 28);

    auto ta_z_zzz_0 = pbuffer.data(idx_npot_0_pf + 29);

    // Set up components of auxiliary buffer : PF

    auto ta_x_xxx_1 = pbuffer.data(idx_npot_1_pf);

    auto ta_x_xxy_1 = pbuffer.data(idx_npot_1_pf + 1);

    auto ta_x_xxz_1 = pbuffer.data(idx_npot_1_pf + 2);

    auto ta_x_xyy_1 = pbuffer.data(idx_npot_1_pf + 3);

    auto ta_x_xyz_1 = pbuffer.data(idx_npot_1_pf + 4);

    auto ta_x_xzz_1 = pbuffer.data(idx_npot_1_pf + 5);

    auto ta_x_yyy_1 = pbuffer.data(idx_npot_1_pf + 6);

    auto ta_x_yyz_1 = pbuffer.data(idx_npot_1_pf + 7);

    auto ta_x_yzz_1 = pbuffer.data(idx_npot_1_pf + 8);

    auto ta_x_zzz_1 = pbuffer.data(idx_npot_1_pf + 9);

    auto ta_y_xxx_1 = pbuffer.data(idx_npot_1_pf + 10);

    auto ta_y_xxy_1 = pbuffer.data(idx_npot_1_pf + 11);

    auto ta_y_xxz_1 = pbuffer.data(idx_npot_1_pf + 12);

    auto ta_y_xyy_1 = pbuffer.data(idx_npot_1_pf + 13);

    auto ta_y_xyz_1 = pbuffer.data(idx_npot_1_pf + 14);

    auto ta_y_xzz_1 = pbuffer.data(idx_npot_1_pf + 15);

    auto ta_y_yyy_1 = pbuffer.data(idx_npot_1_pf + 16);

    auto ta_y_yyz_1 = pbuffer.data(idx_npot_1_pf + 17);

    auto ta_y_yzz_1 = pbuffer.data(idx_npot_1_pf + 18);

    auto ta_y_zzz_1 = pbuffer.data(idx_npot_1_pf + 19);

    auto ta_z_xxx_1 = pbuffer.data(idx_npot_1_pf + 20);

    auto ta_z_xxy_1 = pbuffer.data(idx_npot_1_pf + 21);

    auto ta_z_xxz_1 = pbuffer.data(idx_npot_1_pf + 22);

    auto ta_z_xyy_1 = pbuffer.data(idx_npot_1_pf + 23);

    auto ta_z_xyz_1 = pbuffer.data(idx_npot_1_pf + 24);

    auto ta_z_xzz_1 = pbuffer.data(idx_npot_1_pf + 25);

    auto ta_z_yyy_1 = pbuffer.data(idx_npot_1_pf + 26);

    auto ta_z_yyz_1 = pbuffer.data(idx_npot_1_pf + 27);

    auto ta_z_yzz_1 = pbuffer.data(idx_npot_1_pf + 28);

    auto ta_z_zzz_1 = pbuffer.data(idx_npot_1_pf + 29);

    // Set up 0-10 components of targeted buffer : DF

    auto ta_xx_xxx_0 = pbuffer.data(idx_npot_0_df);

    auto ta_xx_xxy_0 = pbuffer.data(idx_npot_0_df + 1);

    auto ta_xx_xxz_0 = pbuffer.data(idx_npot_0_df + 2);

    auto ta_xx_xyy_0 = pbuffer.data(idx_npot_0_df + 3);

    auto ta_xx_xyz_0 = pbuffer.data(idx_npot_0_df + 4);

    auto ta_xx_xzz_0 = pbuffer.data(idx_npot_0_df + 5);

    auto ta_xx_yyy_0 = pbuffer.data(idx_npot_0_df + 6);

    auto ta_xx_yyz_0 = pbuffer.data(idx_npot_0_df + 7);

    auto ta_xx_yzz_0 = pbuffer.data(idx_npot_0_df + 8);

    auto ta_xx_zzz_0 = pbuffer.data(idx_npot_0_df + 9);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxy_0,  \
                             ta_0_xxy_1,  \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyz_0,  \
                             ta_0_xyz_1,  \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_x_xx_0,   \
                             ta_x_xx_1,   \
                             ta_x_xxx_0,  \
                             ta_x_xxx_1,  \
                             ta_x_xxy_0,  \
                             ta_x_xxy_1,  \
                             ta_x_xxz_0,  \
                             ta_x_xxz_1,  \
                             ta_x_xy_0,   \
                             ta_x_xy_1,   \
                             ta_x_xyy_0,  \
                             ta_x_xyy_1,  \
                             ta_x_xyz_0,  \
                             ta_x_xyz_1,  \
                             ta_x_xz_0,   \
                             ta_x_xz_1,   \
                             ta_x_xzz_0,  \
                             ta_x_xzz_1,  \
                             ta_x_yy_0,   \
                             ta_x_yy_1,   \
                             ta_x_yyy_0,  \
                             ta_x_yyy_1,  \
                             ta_x_yyz_0,  \
                             ta_x_yyz_1,  \
                             ta_x_yz_0,   \
                             ta_x_yz_1,   \
                             ta_x_yzz_0,  \
                             ta_x_yzz_1,  \
                             ta_x_zz_0,   \
                             ta_x_zz_1,   \
                             ta_x_zzz_0,  \
                             ta_x_zzz_1,  \
                             ta_xx_xxx_0, \
                             ta_xx_xxy_0, \
                             ta_xx_xxz_0, \
                             ta_xx_xyy_0, \
                             ta_xx_xyz_0, \
                             ta_xx_xzz_0, \
                             ta_xx_yyy_0, \
                             ta_xx_yyz_0, \
                             ta_xx_yzz_0, \
                             ta_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xx_xxx_0[i] = ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + 3.0 * ta_x_xx_0[i] * fe_0 -
                         3.0 * ta_x_xx_1[i] * fe_0 + ta_x_xxx_0[i] * pa_x[i] - ta_x_xxx_1[i] * pc_x[i];

        ta_xx_xxy_0[i] = ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + 2.0 * ta_x_xy_0[i] * fe_0 -
                         2.0 * ta_x_xy_1[i] * fe_0 + ta_x_xxy_0[i] * pa_x[i] - ta_x_xxy_1[i] * pc_x[i];

        ta_xx_xxz_0[i] = ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + 2.0 * ta_x_xz_0[i] * fe_0 -
                         2.0 * ta_x_xz_1[i] * fe_0 + ta_x_xxz_0[i] * pa_x[i] - ta_x_xxz_1[i] * pc_x[i];

        ta_xx_xyy_0[i] = ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + ta_x_yy_0[i] * fe_0 - ta_x_yy_1[i] * fe_0 +
                         ta_x_xyy_0[i] * pa_x[i] - ta_x_xyy_1[i] * pc_x[i];

        ta_xx_xyz_0[i] = ta_0_xyz_0[i] * fe_0 - ta_0_xyz_1[i] * fe_0 + ta_x_yz_0[i] * fe_0 - ta_x_yz_1[i] * fe_0 +
                         ta_x_xyz_0[i] * pa_x[i] - ta_x_xyz_1[i] * pc_x[i];

        ta_xx_xzz_0[i] = ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + ta_x_zz_0[i] * fe_0 - ta_x_zz_1[i] * fe_0 +
                         ta_x_xzz_0[i] * pa_x[i] - ta_x_xzz_1[i] * pc_x[i];

        ta_xx_yyy_0[i] =
            ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_x_yyy_0[i] * pa_x[i] - ta_x_yyy_1[i] * pc_x[i];

        ta_xx_yyz_0[i] =
            ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + ta_x_yyz_0[i] * pa_x[i] - ta_x_yyz_1[i] * pc_x[i];

        ta_xx_yzz_0[i] =
            ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + ta_x_yzz_0[i] * pa_x[i] - ta_x_yzz_1[i] * pc_x[i];

        ta_xx_zzz_0[i] =
            ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_x_zzz_0[i] * pa_x[i] - ta_x_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto ta_xy_xxx_0 = pbuffer.data(idx_npot_0_df + 10);

    auto ta_xy_xxy_0 = pbuffer.data(idx_npot_0_df + 11);

    auto ta_xy_xxz_0 = pbuffer.data(idx_npot_0_df + 12);

    auto ta_xy_xyy_0 = pbuffer.data(idx_npot_0_df + 13);

    auto ta_xy_xyz_0 = pbuffer.data(idx_npot_0_df + 14);

    auto ta_xy_xzz_0 = pbuffer.data(idx_npot_0_df + 15);

    auto ta_xy_yyy_0 = pbuffer.data(idx_npot_0_df + 16);

    auto ta_xy_yyz_0 = pbuffer.data(idx_npot_0_df + 17);

    auto ta_xy_yzz_0 = pbuffer.data(idx_npot_0_df + 18);

    auto ta_xy_zzz_0 = pbuffer.data(idx_npot_0_df + 19);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pc_x,        \
                             pc_y,        \
                             ta_x_xxx_0,  \
                             ta_x_xxx_1,  \
                             ta_x_xxz_0,  \
                             ta_x_xxz_1,  \
                             ta_x_xzz_0,  \
                             ta_x_xzz_1,  \
                             ta_xy_xxx_0, \
                             ta_xy_xxy_0, \
                             ta_xy_xxz_0, \
                             ta_xy_xyy_0, \
                             ta_xy_xyz_0, \
                             ta_xy_xzz_0, \
                             ta_xy_yyy_0, \
                             ta_xy_yyz_0, \
                             ta_xy_yzz_0, \
                             ta_xy_zzz_0, \
                             ta_y_xxy_0,  \
                             ta_y_xxy_1,  \
                             ta_y_xy_0,   \
                             ta_y_xy_1,   \
                             ta_y_xyy_0,  \
                             ta_y_xyy_1,  \
                             ta_y_xyz_0,  \
                             ta_y_xyz_1,  \
                             ta_y_yy_0,   \
                             ta_y_yy_1,   \
                             ta_y_yyy_0,  \
                             ta_y_yyy_1,  \
                             ta_y_yyz_0,  \
                             ta_y_yyz_1,  \
                             ta_y_yz_0,   \
                             ta_y_yz_1,   \
                             ta_y_yzz_0,  \
                             ta_y_yzz_1,  \
                             ta_y_zzz_0,  \
                             ta_y_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xy_xxx_0[i] = ta_x_xxx_0[i] * pa_y[i] - ta_x_xxx_1[i] * pc_y[i];

        ta_xy_xxy_0[i] =
            2.0 * ta_y_xy_0[i] * fe_0 - 2.0 * ta_y_xy_1[i] * fe_0 + ta_y_xxy_0[i] * pa_x[i] - ta_y_xxy_1[i] * pc_x[i];

        ta_xy_xxz_0[i] = ta_x_xxz_0[i] * pa_y[i] - ta_x_xxz_1[i] * pc_y[i];

        ta_xy_xyy_0[i] = ta_y_yy_0[i] * fe_0 - ta_y_yy_1[i] * fe_0 + ta_y_xyy_0[i] * pa_x[i] - ta_y_xyy_1[i] * pc_x[i];

        ta_xy_xyz_0[i] = ta_y_yz_0[i] * fe_0 - ta_y_yz_1[i] * fe_0 + ta_y_xyz_0[i] * pa_x[i] - ta_y_xyz_1[i] * pc_x[i];

        ta_xy_xzz_0[i] = ta_x_xzz_0[i] * pa_y[i] - ta_x_xzz_1[i] * pc_y[i];

        ta_xy_yyy_0[i] = ta_y_yyy_0[i] * pa_x[i] - ta_y_yyy_1[i] * pc_x[i];

        ta_xy_yyz_0[i] = ta_y_yyz_0[i] * pa_x[i] - ta_y_yyz_1[i] * pc_x[i];

        ta_xy_yzz_0[i] = ta_y_yzz_0[i] * pa_x[i] - ta_y_yzz_1[i] * pc_x[i];

        ta_xy_zzz_0[i] = ta_y_zzz_0[i] * pa_x[i] - ta_y_zzz_1[i] * pc_x[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto ta_xz_xxx_0 = pbuffer.data(idx_npot_0_df + 20);

    auto ta_xz_xxy_0 = pbuffer.data(idx_npot_0_df + 21);

    auto ta_xz_xxz_0 = pbuffer.data(idx_npot_0_df + 22);

    auto ta_xz_xyy_0 = pbuffer.data(idx_npot_0_df + 23);

    auto ta_xz_xyz_0 = pbuffer.data(idx_npot_0_df + 24);

    auto ta_xz_xzz_0 = pbuffer.data(idx_npot_0_df + 25);

    auto ta_xz_yyy_0 = pbuffer.data(idx_npot_0_df + 26);

    auto ta_xz_yyz_0 = pbuffer.data(idx_npot_0_df + 27);

    auto ta_xz_yzz_0 = pbuffer.data(idx_npot_0_df + 28);

    auto ta_xz_zzz_0 = pbuffer.data(idx_npot_0_df + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             pc_x,        \
                             pc_z,        \
                             ta_x_xxx_0,  \
                             ta_x_xxx_1,  \
                             ta_x_xxy_0,  \
                             ta_x_xxy_1,  \
                             ta_x_xyy_0,  \
                             ta_x_xyy_1,  \
                             ta_xz_xxx_0, \
                             ta_xz_xxy_0, \
                             ta_xz_xxz_0, \
                             ta_xz_xyy_0, \
                             ta_xz_xyz_0, \
                             ta_xz_xzz_0, \
                             ta_xz_yyy_0, \
                             ta_xz_yyz_0, \
                             ta_xz_yzz_0, \
                             ta_xz_zzz_0, \
                             ta_z_xxz_0,  \
                             ta_z_xxz_1,  \
                             ta_z_xyz_0,  \
                             ta_z_xyz_1,  \
                             ta_z_xz_0,   \
                             ta_z_xz_1,   \
                             ta_z_xzz_0,  \
                             ta_z_xzz_1,  \
                             ta_z_yyy_0,  \
                             ta_z_yyy_1,  \
                             ta_z_yyz_0,  \
                             ta_z_yyz_1,  \
                             ta_z_yz_0,   \
                             ta_z_yz_1,   \
                             ta_z_yzz_0,  \
                             ta_z_yzz_1,  \
                             ta_z_zz_0,   \
                             ta_z_zz_1,   \
                             ta_z_zzz_0,  \
                             ta_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xz_xxx_0[i] = ta_x_xxx_0[i] * pa_z[i] - ta_x_xxx_1[i] * pc_z[i];

        ta_xz_xxy_0[i] = ta_x_xxy_0[i] * pa_z[i] - ta_x_xxy_1[i] * pc_z[i];

        ta_xz_xxz_0[i] =
            2.0 * ta_z_xz_0[i] * fe_0 - 2.0 * ta_z_xz_1[i] * fe_0 + ta_z_xxz_0[i] * pa_x[i] - ta_z_xxz_1[i] * pc_x[i];

        ta_xz_xyy_0[i] = ta_x_xyy_0[i] * pa_z[i] - ta_x_xyy_1[i] * pc_z[i];

        ta_xz_xyz_0[i] = ta_z_yz_0[i] * fe_0 - ta_z_yz_1[i] * fe_0 + ta_z_xyz_0[i] * pa_x[i] - ta_z_xyz_1[i] * pc_x[i];

        ta_xz_xzz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_z_xzz_0[i] * pa_x[i] - ta_z_xzz_1[i] * pc_x[i];

        ta_xz_yyy_0[i] = ta_z_yyy_0[i] * pa_x[i] - ta_z_yyy_1[i] * pc_x[i];

        ta_xz_yyz_0[i] = ta_z_yyz_0[i] * pa_x[i] - ta_z_yyz_1[i] * pc_x[i];

        ta_xz_yzz_0[i] = ta_z_yzz_0[i] * pa_x[i] - ta_z_yzz_1[i] * pc_x[i];

        ta_xz_zzz_0[i] = ta_z_zzz_0[i] * pa_x[i] - ta_z_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto ta_yy_xxx_0 = pbuffer.data(idx_npot_0_df + 30);

    auto ta_yy_xxy_0 = pbuffer.data(idx_npot_0_df + 31);

    auto ta_yy_xxz_0 = pbuffer.data(idx_npot_0_df + 32);

    auto ta_yy_xyy_0 = pbuffer.data(idx_npot_0_df + 33);

    auto ta_yy_xyz_0 = pbuffer.data(idx_npot_0_df + 34);

    auto ta_yy_xzz_0 = pbuffer.data(idx_npot_0_df + 35);

    auto ta_yy_yyy_0 = pbuffer.data(idx_npot_0_df + 36);

    auto ta_yy_yyz_0 = pbuffer.data(idx_npot_0_df + 37);

    auto ta_yy_yzz_0 = pbuffer.data(idx_npot_0_df + 38);

    auto ta_yy_zzz_0 = pbuffer.data(idx_npot_0_df + 39);

#pragma omp simd aligned(pa_y,            \
                             pc_y,        \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxy_0,  \
                             ta_0_xxy_1,  \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyz_0,  \
                             ta_0_xyz_1,  \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_y_xx_0,   \
                             ta_y_xx_1,   \
                             ta_y_xxx_0,  \
                             ta_y_xxx_1,  \
                             ta_y_xxy_0,  \
                             ta_y_xxy_1,  \
                             ta_y_xxz_0,  \
                             ta_y_xxz_1,  \
                             ta_y_xy_0,   \
                             ta_y_xy_1,   \
                             ta_y_xyy_0,  \
                             ta_y_xyy_1,  \
                             ta_y_xyz_0,  \
                             ta_y_xyz_1,  \
                             ta_y_xz_0,   \
                             ta_y_xz_1,   \
                             ta_y_xzz_0,  \
                             ta_y_xzz_1,  \
                             ta_y_yy_0,   \
                             ta_y_yy_1,   \
                             ta_y_yyy_0,  \
                             ta_y_yyy_1,  \
                             ta_y_yyz_0,  \
                             ta_y_yyz_1,  \
                             ta_y_yz_0,   \
                             ta_y_yz_1,   \
                             ta_y_yzz_0,  \
                             ta_y_yzz_1,  \
                             ta_y_zz_0,   \
                             ta_y_zz_1,   \
                             ta_y_zzz_0,  \
                             ta_y_zzz_1,  \
                             ta_yy_xxx_0, \
                             ta_yy_xxy_0, \
                             ta_yy_xxz_0, \
                             ta_yy_xyy_0, \
                             ta_yy_xyz_0, \
                             ta_yy_xzz_0, \
                             ta_yy_yyy_0, \
                             ta_yy_yyz_0, \
                             ta_yy_yzz_0, \
                             ta_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yy_xxx_0[i] =
            ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + ta_y_xxx_0[i] * pa_y[i] - ta_y_xxx_1[i] * pc_y[i];

        ta_yy_xxy_0[i] = ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + ta_y_xx_0[i] * fe_0 - ta_y_xx_1[i] * fe_0 +
                         ta_y_xxy_0[i] * pa_y[i] - ta_y_xxy_1[i] * pc_y[i];

        ta_yy_xxz_0[i] =
            ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + ta_y_xxz_0[i] * pa_y[i] - ta_y_xxz_1[i] * pc_y[i];

        ta_yy_xyy_0[i] = ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + 2.0 * ta_y_xy_0[i] * fe_0 -
                         2.0 * ta_y_xy_1[i] * fe_0 + ta_y_xyy_0[i] * pa_y[i] - ta_y_xyy_1[i] * pc_y[i];

        ta_yy_xyz_0[i] = ta_0_xyz_0[i] * fe_0 - ta_0_xyz_1[i] * fe_0 + ta_y_xz_0[i] * fe_0 - ta_y_xz_1[i] * fe_0 +
                         ta_y_xyz_0[i] * pa_y[i] - ta_y_xyz_1[i] * pc_y[i];

        ta_yy_xzz_0[i] =
            ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + ta_y_xzz_0[i] * pa_y[i] - ta_y_xzz_1[i] * pc_y[i];

        ta_yy_yyy_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + 3.0 * ta_y_yy_0[i] * fe_0 -
                         3.0 * ta_y_yy_1[i] * fe_0 + ta_y_yyy_0[i] * pa_y[i] - ta_y_yyy_1[i] * pc_y[i];

        ta_yy_yyz_0[i] = ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + 2.0 * ta_y_yz_0[i] * fe_0 -
                         2.0 * ta_y_yz_1[i] * fe_0 + ta_y_yyz_0[i] * pa_y[i] - ta_y_yyz_1[i] * pc_y[i];

        ta_yy_yzz_0[i] = ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + ta_y_zz_0[i] * fe_0 - ta_y_zz_1[i] * fe_0 +
                         ta_y_yzz_0[i] * pa_y[i] - ta_y_yzz_1[i] * pc_y[i];

        ta_yy_zzz_0[i] =
            ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_y_zzz_0[i] * pa_y[i] - ta_y_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto ta_yz_xxx_0 = pbuffer.data(idx_npot_0_df + 40);

    auto ta_yz_xxy_0 = pbuffer.data(idx_npot_0_df + 41);

    auto ta_yz_xxz_0 = pbuffer.data(idx_npot_0_df + 42);

    auto ta_yz_xyy_0 = pbuffer.data(idx_npot_0_df + 43);

    auto ta_yz_xyz_0 = pbuffer.data(idx_npot_0_df + 44);

    auto ta_yz_xzz_0 = pbuffer.data(idx_npot_0_df + 45);

    auto ta_yz_yyy_0 = pbuffer.data(idx_npot_0_df + 46);

    auto ta_yz_yyz_0 = pbuffer.data(idx_npot_0_df + 47);

    auto ta_yz_yzz_0 = pbuffer.data(idx_npot_0_df + 48);

    auto ta_yz_zzz_0 = pbuffer.data(idx_npot_0_df + 49);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             pc_y,        \
                             pc_z,        \
                             ta_y_xxy_0,  \
                             ta_y_xxy_1,  \
                             ta_y_xyy_0,  \
                             ta_y_xyy_1,  \
                             ta_y_yyy_0,  \
                             ta_y_yyy_1,  \
                             ta_yz_xxx_0, \
                             ta_yz_xxy_0, \
                             ta_yz_xxz_0, \
                             ta_yz_xyy_0, \
                             ta_yz_xyz_0, \
                             ta_yz_xzz_0, \
                             ta_yz_yyy_0, \
                             ta_yz_yyz_0, \
                             ta_yz_yzz_0, \
                             ta_yz_zzz_0, \
                             ta_z_xxx_0,  \
                             ta_z_xxx_1,  \
                             ta_z_xxz_0,  \
                             ta_z_xxz_1,  \
                             ta_z_xyz_0,  \
                             ta_z_xyz_1,  \
                             ta_z_xz_0,   \
                             ta_z_xz_1,   \
                             ta_z_xzz_0,  \
                             ta_z_xzz_1,  \
                             ta_z_yyz_0,  \
                             ta_z_yyz_1,  \
                             ta_z_yz_0,   \
                             ta_z_yz_1,   \
                             ta_z_yzz_0,  \
                             ta_z_yzz_1,  \
                             ta_z_zz_0,   \
                             ta_z_zz_1,   \
                             ta_z_zzz_0,  \
                             ta_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yz_xxx_0[i] = ta_z_xxx_0[i] * pa_y[i] - ta_z_xxx_1[i] * pc_y[i];

        ta_yz_xxy_0[i] = ta_y_xxy_0[i] * pa_z[i] - ta_y_xxy_1[i] * pc_z[i];

        ta_yz_xxz_0[i] = ta_z_xxz_0[i] * pa_y[i] - ta_z_xxz_1[i] * pc_y[i];

        ta_yz_xyy_0[i] = ta_y_xyy_0[i] * pa_z[i] - ta_y_xyy_1[i] * pc_z[i];

        ta_yz_xyz_0[i] = ta_z_xz_0[i] * fe_0 - ta_z_xz_1[i] * fe_0 + ta_z_xyz_0[i] * pa_y[i] - ta_z_xyz_1[i] * pc_y[i];

        ta_yz_xzz_0[i] = ta_z_xzz_0[i] * pa_y[i] - ta_z_xzz_1[i] * pc_y[i];

        ta_yz_yyy_0[i] = ta_y_yyy_0[i] * pa_z[i] - ta_y_yyy_1[i] * pc_z[i];

        ta_yz_yyz_0[i] =
            2.0 * ta_z_yz_0[i] * fe_0 - 2.0 * ta_z_yz_1[i] * fe_0 + ta_z_yyz_0[i] * pa_y[i] - ta_z_yyz_1[i] * pc_y[i];

        ta_yz_yzz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_z_yzz_0[i] * pa_y[i] - ta_z_yzz_1[i] * pc_y[i];

        ta_yz_zzz_0[i] = ta_z_zzz_0[i] * pa_y[i] - ta_z_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto ta_zz_xxx_0 = pbuffer.data(idx_npot_0_df + 50);

    auto ta_zz_xxy_0 = pbuffer.data(idx_npot_0_df + 51);

    auto ta_zz_xxz_0 = pbuffer.data(idx_npot_0_df + 52);

    auto ta_zz_xyy_0 = pbuffer.data(idx_npot_0_df + 53);

    auto ta_zz_xyz_0 = pbuffer.data(idx_npot_0_df + 54);

    auto ta_zz_xzz_0 = pbuffer.data(idx_npot_0_df + 55);

    auto ta_zz_yyy_0 = pbuffer.data(idx_npot_0_df + 56);

    auto ta_zz_yyz_0 = pbuffer.data(idx_npot_0_df + 57);

    auto ta_zz_yzz_0 = pbuffer.data(idx_npot_0_df + 58);

    auto ta_zz_zzz_0 = pbuffer.data(idx_npot_0_df + 59);

#pragma omp simd aligned(pa_z,            \
                             pc_z,        \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxy_0,  \
                             ta_0_xxy_1,  \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyz_0,  \
                             ta_0_xyz_1,  \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_z_xx_0,   \
                             ta_z_xx_1,   \
                             ta_z_xxx_0,  \
                             ta_z_xxx_1,  \
                             ta_z_xxy_0,  \
                             ta_z_xxy_1,  \
                             ta_z_xxz_0,  \
                             ta_z_xxz_1,  \
                             ta_z_xy_0,   \
                             ta_z_xy_1,   \
                             ta_z_xyy_0,  \
                             ta_z_xyy_1,  \
                             ta_z_xyz_0,  \
                             ta_z_xyz_1,  \
                             ta_z_xz_0,   \
                             ta_z_xz_1,   \
                             ta_z_xzz_0,  \
                             ta_z_xzz_1,  \
                             ta_z_yy_0,   \
                             ta_z_yy_1,   \
                             ta_z_yyy_0,  \
                             ta_z_yyy_1,  \
                             ta_z_yyz_0,  \
                             ta_z_yyz_1,  \
                             ta_z_yz_0,   \
                             ta_z_yz_1,   \
                             ta_z_yzz_0,  \
                             ta_z_yzz_1,  \
                             ta_z_zz_0,   \
                             ta_z_zz_1,   \
                             ta_z_zzz_0,  \
                             ta_z_zzz_1,  \
                             ta_zz_xxx_0, \
                             ta_zz_xxy_0, \
                             ta_zz_xxz_0, \
                             ta_zz_xyy_0, \
                             ta_zz_xyz_0, \
                             ta_zz_xzz_0, \
                             ta_zz_yyy_0, \
                             ta_zz_yyz_0, \
                             ta_zz_yzz_0, \
                             ta_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zz_xxx_0[i] =
            ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + ta_z_xxx_0[i] * pa_z[i] - ta_z_xxx_1[i] * pc_z[i];

        ta_zz_xxy_0[i] =
            ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + ta_z_xxy_0[i] * pa_z[i] - ta_z_xxy_1[i] * pc_z[i];

        ta_zz_xxz_0[i] = ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + ta_z_xx_0[i] * fe_0 - ta_z_xx_1[i] * fe_0 +
                         ta_z_xxz_0[i] * pa_z[i] - ta_z_xxz_1[i] * pc_z[i];

        ta_zz_xyy_0[i] =
            ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + ta_z_xyy_0[i] * pa_z[i] - ta_z_xyy_1[i] * pc_z[i];

        ta_zz_xyz_0[i] = ta_0_xyz_0[i] * fe_0 - ta_0_xyz_1[i] * fe_0 + ta_z_xy_0[i] * fe_0 - ta_z_xy_1[i] * fe_0 +
                         ta_z_xyz_0[i] * pa_z[i] - ta_z_xyz_1[i] * pc_z[i];

        ta_zz_xzz_0[i] = ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + 2.0 * ta_z_xz_0[i] * fe_0 -
                         2.0 * ta_z_xz_1[i] * fe_0 + ta_z_xzz_0[i] * pa_z[i] - ta_z_xzz_1[i] * pc_z[i];

        ta_zz_yyy_0[i] =
            ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_z_yyy_0[i] * pa_z[i] - ta_z_yyy_1[i] * pc_z[i];

        ta_zz_yyz_0[i] = ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + ta_z_yy_0[i] * fe_0 - ta_z_yy_1[i] * fe_0 +
                         ta_z_yyz_0[i] * pa_z[i] - ta_z_yyz_1[i] * pc_z[i];

        ta_zz_yzz_0[i] = ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + 2.0 * ta_z_yz_0[i] * fe_0 -
                         2.0 * ta_z_yz_1[i] * fe_0 + ta_z_yzz_0[i] * pa_z[i] - ta_z_yzz_1[i] * pc_z[i];

        ta_zz_zzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + 3.0 * ta_z_zz_0[i] * fe_0 -
                         3.0 * ta_z_zz_1[i] * fe_0 + ta_z_zzz_0[i] * pa_z[i] - ta_z_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
