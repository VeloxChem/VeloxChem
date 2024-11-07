#include "NuclearPotentialGeom020PrimRecSD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_sd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_sd,
                                        const size_t              idx_npot_geom_020_0_ss,
                                        const size_t              idx_npot_geom_020_1_ss,
                                        const size_t              idx_npot_geom_010_1_sp,
                                        const size_t              idx_npot_geom_020_0_sp,
                                        const size_t              idx_npot_geom_020_1_sp,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpb,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

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

    // Set up components of targeted buffer : SD

    auto ta2_xx_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd);

    auto ta2_xx_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 1);

    auto ta2_xx_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 2);

    auto ta2_xx_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 3);

    auto ta2_xx_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 4);

    auto ta2_xx_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 5);

    auto ta2_xy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 6);

    auto ta2_xy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 7);

    auto ta2_xy_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 8);

    auto ta2_xy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 9);

    auto ta2_xy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 10);

    auto ta2_xy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 11);

    auto ta2_xz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 12);

    auto ta2_xz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 13);

    auto ta2_xz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 14);

    auto ta2_xz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 15);

    auto ta2_xz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 16);

    auto ta2_xz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 17);

    auto ta2_yy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 18);

    auto ta2_yy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 19);

    auto ta2_yy_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 20);

    auto ta2_yy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 21);

    auto ta2_yy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 22);

    auto ta2_yy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 23);

    auto ta2_yz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 24);

    auto ta2_yz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 25);

    auto ta2_yz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 26);

    auto ta2_yz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 27);

    auto ta2_yz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 28);

    auto ta2_yz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 29);

    auto ta2_zz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 30);

    auto ta2_zz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 31);

    auto ta2_zz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 32);

    auto ta2_zz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 33);

    auto ta2_zz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 34);

    auto ta2_zz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 35);

#pragma omp simd aligned(pb_x,              \
                             pb_y,          \
                             pb_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_x_0_x_1,   \
                             ta1_x_0_y_1,   \
                             ta1_x_0_z_1,   \
                             ta1_y_0_x_1,   \
                             ta1_y_0_y_1,   \
                             ta1_y_0_z_1,   \
                             ta1_z_0_x_1,   \
                             ta1_z_0_y_1,   \
                             ta1_z_0_z_1,   \
                             ta2_xx_0_0_0,  \
                             ta2_xx_0_0_1,  \
                             ta2_xx_0_x_0,  \
                             ta2_xx_0_x_1,  \
                             ta2_xx_0_xx_0, \
                             ta2_xx_0_xy_0, \
                             ta2_xx_0_xz_0, \
                             ta2_xx_0_y_0,  \
                             ta2_xx_0_y_1,  \
                             ta2_xx_0_yy_0, \
                             ta2_xx_0_yz_0, \
                             ta2_xx_0_z_0,  \
                             ta2_xx_0_z_1,  \
                             ta2_xx_0_zz_0, \
                             ta2_xy_0_0_0,  \
                             ta2_xy_0_0_1,  \
                             ta2_xy_0_x_0,  \
                             ta2_xy_0_x_1,  \
                             ta2_xy_0_xx_0, \
                             ta2_xy_0_xy_0, \
                             ta2_xy_0_xz_0, \
                             ta2_xy_0_y_0,  \
                             ta2_xy_0_y_1,  \
                             ta2_xy_0_yy_0, \
                             ta2_xy_0_yz_0, \
                             ta2_xy_0_z_0,  \
                             ta2_xy_0_z_1,  \
                             ta2_xy_0_zz_0, \
                             ta2_xz_0_0_0,  \
                             ta2_xz_0_0_1,  \
                             ta2_xz_0_x_0,  \
                             ta2_xz_0_x_1,  \
                             ta2_xz_0_xx_0, \
                             ta2_xz_0_xy_0, \
                             ta2_xz_0_xz_0, \
                             ta2_xz_0_y_0,  \
                             ta2_xz_0_y_1,  \
                             ta2_xz_0_yy_0, \
                             ta2_xz_0_yz_0, \
                             ta2_xz_0_z_0,  \
                             ta2_xz_0_z_1,  \
                             ta2_xz_0_zz_0, \
                             ta2_yy_0_0_0,  \
                             ta2_yy_0_0_1,  \
                             ta2_yy_0_x_0,  \
                             ta2_yy_0_x_1,  \
                             ta2_yy_0_xx_0, \
                             ta2_yy_0_xy_0, \
                             ta2_yy_0_xz_0, \
                             ta2_yy_0_y_0,  \
                             ta2_yy_0_y_1,  \
                             ta2_yy_0_yy_0, \
                             ta2_yy_0_yz_0, \
                             ta2_yy_0_z_0,  \
                             ta2_yy_0_z_1,  \
                             ta2_yy_0_zz_0, \
                             ta2_yz_0_0_0,  \
                             ta2_yz_0_0_1,  \
                             ta2_yz_0_x_0,  \
                             ta2_yz_0_x_1,  \
                             ta2_yz_0_xx_0, \
                             ta2_yz_0_xy_0, \
                             ta2_yz_0_xz_0, \
                             ta2_yz_0_y_0,  \
                             ta2_yz_0_y_1,  \
                             ta2_yz_0_yy_0, \
                             ta2_yz_0_yz_0, \
                             ta2_yz_0_z_0,  \
                             ta2_yz_0_z_1,  \
                             ta2_yz_0_zz_0, \
                             ta2_zz_0_0_0,  \
                             ta2_zz_0_0_1,  \
                             ta2_zz_0_x_0,  \
                             ta2_zz_0_x_1,  \
                             ta2_zz_0_xx_0, \
                             ta2_zz_0_xy_0, \
                             ta2_zz_0_xz_0, \
                             ta2_zz_0_y_0,  \
                             ta2_zz_0_y_1,  \
                             ta2_zz_0_yy_0, \
                             ta2_zz_0_yz_0, \
                             ta2_zz_0_z_0,  \
                             ta2_zz_0_z_1,  \
                             ta2_zz_0_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_0_xx_0[i] =
            ta2_xx_0_0_0[i] * fe_0 - ta2_xx_0_0_1[i] * fe_0 + 2.0 * ta1_x_0_x_1[i] + ta2_xx_0_x_0[i] * pb_x[i] - ta2_xx_0_x_1[i] * pc_x[i];

        ta2_xx_0_xy_0[i] = ta2_xx_0_x_0[i] * pb_y[i] - ta2_xx_0_x_1[i] * pc_y[i];

        ta2_xx_0_xz_0[i] = ta2_xx_0_x_0[i] * pb_z[i] - ta2_xx_0_x_1[i] * pc_z[i];

        ta2_xx_0_yy_0[i] = ta2_xx_0_0_0[i] * fe_0 - ta2_xx_0_0_1[i] * fe_0 + ta2_xx_0_y_0[i] * pb_y[i] - ta2_xx_0_y_1[i] * pc_y[i];

        ta2_xx_0_yz_0[i] = ta2_xx_0_z_0[i] * pb_y[i] - ta2_xx_0_z_1[i] * pc_y[i];

        ta2_xx_0_zz_0[i] = ta2_xx_0_0_0[i] * fe_0 - ta2_xx_0_0_1[i] * fe_0 + ta2_xx_0_z_0[i] * pb_z[i] - ta2_xx_0_z_1[i] * pc_z[i];

        ta2_xy_0_xx_0[i] = ta2_xy_0_0_0[i] * fe_0 - ta2_xy_0_0_1[i] * fe_0 + ta1_y_0_x_1[i] + ta2_xy_0_x_0[i] * pb_x[i] - ta2_xy_0_x_1[i] * pc_x[i];

        ta2_xy_0_xy_0[i] = ta1_y_0_y_1[i] + ta2_xy_0_y_0[i] * pb_x[i] - ta2_xy_0_y_1[i] * pc_x[i];

        ta2_xy_0_xz_0[i] = ta2_xy_0_x_0[i] * pb_z[i] - ta2_xy_0_x_1[i] * pc_z[i];

        ta2_xy_0_yy_0[i] = ta2_xy_0_0_0[i] * fe_0 - ta2_xy_0_0_1[i] * fe_0 + ta1_x_0_y_1[i] + ta2_xy_0_y_0[i] * pb_y[i] - ta2_xy_0_y_1[i] * pc_y[i];

        ta2_xy_0_yz_0[i] = ta2_xy_0_y_0[i] * pb_z[i] - ta2_xy_0_y_1[i] * pc_z[i];

        ta2_xy_0_zz_0[i] = ta2_xy_0_0_0[i] * fe_0 - ta2_xy_0_0_1[i] * fe_0 + ta2_xy_0_z_0[i] * pb_z[i] - ta2_xy_0_z_1[i] * pc_z[i];

        ta2_xz_0_xx_0[i] = ta2_xz_0_0_0[i] * fe_0 - ta2_xz_0_0_1[i] * fe_0 + ta1_z_0_x_1[i] + ta2_xz_0_x_0[i] * pb_x[i] - ta2_xz_0_x_1[i] * pc_x[i];

        ta2_xz_0_xy_0[i] = ta2_xz_0_x_0[i] * pb_y[i] - ta2_xz_0_x_1[i] * pc_y[i];

        ta2_xz_0_xz_0[i] = ta1_z_0_z_1[i] + ta2_xz_0_z_0[i] * pb_x[i] - ta2_xz_0_z_1[i] * pc_x[i];

        ta2_xz_0_yy_0[i] = ta2_xz_0_0_0[i] * fe_0 - ta2_xz_0_0_1[i] * fe_0 + ta2_xz_0_y_0[i] * pb_y[i] - ta2_xz_0_y_1[i] * pc_y[i];

        ta2_xz_0_yz_0[i] = ta2_xz_0_z_0[i] * pb_y[i] - ta2_xz_0_z_1[i] * pc_y[i];

        ta2_xz_0_zz_0[i] = ta2_xz_0_0_0[i] * fe_0 - ta2_xz_0_0_1[i] * fe_0 + ta1_x_0_z_1[i] + ta2_xz_0_z_0[i] * pb_z[i] - ta2_xz_0_z_1[i] * pc_z[i];

        ta2_yy_0_xx_0[i] = ta2_yy_0_0_0[i] * fe_0 - ta2_yy_0_0_1[i] * fe_0 + ta2_yy_0_x_0[i] * pb_x[i] - ta2_yy_0_x_1[i] * pc_x[i];

        ta2_yy_0_xy_0[i] = ta2_yy_0_y_0[i] * pb_x[i] - ta2_yy_0_y_1[i] * pc_x[i];

        ta2_yy_0_xz_0[i] = ta2_yy_0_z_0[i] * pb_x[i] - ta2_yy_0_z_1[i] * pc_x[i];

        ta2_yy_0_yy_0[i] =
            ta2_yy_0_0_0[i] * fe_0 - ta2_yy_0_0_1[i] * fe_0 + 2.0 * ta1_y_0_y_1[i] + ta2_yy_0_y_0[i] * pb_y[i] - ta2_yy_0_y_1[i] * pc_y[i];

        ta2_yy_0_yz_0[i] = ta2_yy_0_y_0[i] * pb_z[i] - ta2_yy_0_y_1[i] * pc_z[i];

        ta2_yy_0_zz_0[i] = ta2_yy_0_0_0[i] * fe_0 - ta2_yy_0_0_1[i] * fe_0 + ta2_yy_0_z_0[i] * pb_z[i] - ta2_yy_0_z_1[i] * pc_z[i];

        ta2_yz_0_xx_0[i] = ta2_yz_0_0_0[i] * fe_0 - ta2_yz_0_0_1[i] * fe_0 + ta2_yz_0_x_0[i] * pb_x[i] - ta2_yz_0_x_1[i] * pc_x[i];

        ta2_yz_0_xy_0[i] = ta2_yz_0_y_0[i] * pb_x[i] - ta2_yz_0_y_1[i] * pc_x[i];

        ta2_yz_0_xz_0[i] = ta2_yz_0_z_0[i] * pb_x[i] - ta2_yz_0_z_1[i] * pc_x[i];

        ta2_yz_0_yy_0[i] = ta2_yz_0_0_0[i] * fe_0 - ta2_yz_0_0_1[i] * fe_0 + ta1_z_0_y_1[i] + ta2_yz_0_y_0[i] * pb_y[i] - ta2_yz_0_y_1[i] * pc_y[i];

        ta2_yz_0_yz_0[i] = ta1_z_0_z_1[i] + ta2_yz_0_z_0[i] * pb_y[i] - ta2_yz_0_z_1[i] * pc_y[i];

        ta2_yz_0_zz_0[i] = ta2_yz_0_0_0[i] * fe_0 - ta2_yz_0_0_1[i] * fe_0 + ta1_y_0_z_1[i] + ta2_yz_0_z_0[i] * pb_z[i] - ta2_yz_0_z_1[i] * pc_z[i];

        ta2_zz_0_xx_0[i] = ta2_zz_0_0_0[i] * fe_0 - ta2_zz_0_0_1[i] * fe_0 + ta2_zz_0_x_0[i] * pb_x[i] - ta2_zz_0_x_1[i] * pc_x[i];

        ta2_zz_0_xy_0[i] = ta2_zz_0_y_0[i] * pb_x[i] - ta2_zz_0_y_1[i] * pc_x[i];

        ta2_zz_0_xz_0[i] = ta2_zz_0_z_0[i] * pb_x[i] - ta2_zz_0_z_1[i] * pc_x[i];

        ta2_zz_0_yy_0[i] = ta2_zz_0_0_0[i] * fe_0 - ta2_zz_0_0_1[i] * fe_0 + ta2_zz_0_y_0[i] * pb_y[i] - ta2_zz_0_y_1[i] * pc_y[i];

        ta2_zz_0_yz_0[i] = ta2_zz_0_z_0[i] * pb_y[i] - ta2_zz_0_z_1[i] * pc_y[i];

        ta2_zz_0_zz_0[i] =
            ta2_zz_0_0_0[i] * fe_0 - ta2_zz_0_0_1[i] * fe_0 + 2.0 * ta1_z_0_z_1[i] + ta2_zz_0_z_0[i] * pb_z[i] - ta2_zz_0_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
