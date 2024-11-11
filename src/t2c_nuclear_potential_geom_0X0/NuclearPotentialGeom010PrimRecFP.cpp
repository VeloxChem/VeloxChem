#include "NuclearPotentialGeom010PrimRecFP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_fp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_fp,
                                        const size_t              idx_npot_geom_010_0_pp,
                                        const size_t              idx_npot_geom_010_1_pp,
                                        const size_t              idx_npot_geom_010_0_ds,
                                        const size_t              idx_npot_geom_010_1_ds,
                                        const size_t              idx_npot_1_dp,
                                        const size_t              idx_npot_geom_010_0_dp,
                                        const size_t              idx_npot_geom_010_1_dp,
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

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds);

    auto ta1_x_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 3);

    auto ta1_x_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 5);

    auto ta1_y_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 6);

    auto ta1_y_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 9);

    auto ta1_y_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 11);

    auto ta1_z_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 12);

    auto ta1_z_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 15);

    auto ta1_z_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 17);

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds);

    auto ta1_x_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 3);

    auto ta1_x_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 5);

    auto ta1_y_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 6);

    auto ta1_y_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 9);

    auto ta1_y_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 11);

    auto ta1_z_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 12);

    auto ta1_z_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 15);

    auto ta1_z_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 17);

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

    // Set up components of auxiliary buffer : DP

    auto ta1_x_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp);

    auto ta1_x_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 1);

    auto ta1_x_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 2);

    auto ta1_x_xy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 3);

    auto ta1_x_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 4);

    auto ta1_x_xz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 6);

    auto ta1_x_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 8);

    auto ta1_x_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 9);

    auto ta1_x_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 10);

    auto ta1_x_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 11);

    auto ta1_x_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 14);

    auto ta1_x_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 15);

    auto ta1_x_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 16);

    auto ta1_x_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 17);

    auto ta1_y_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 18);

    auto ta1_y_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 19);

    auto ta1_y_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 20);

    auto ta1_y_xy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 21);

    auto ta1_y_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 22);

    auto ta1_y_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 26);

    auto ta1_y_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 27);

    auto ta1_y_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 28);

    auto ta1_y_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 29);

    auto ta1_y_yz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 31);

    auto ta1_y_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 32);

    auto ta1_y_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 33);

    auto ta1_y_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 34);

    auto ta1_y_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 35);

    auto ta1_z_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 36);

    auto ta1_z_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 37);

    auto ta1_z_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 38);

    auto ta1_z_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 40);

    auto ta1_z_xz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 42);

    auto ta1_z_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 44);

    auto ta1_z_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 45);

    auto ta1_z_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 46);

    auto ta1_z_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 47);

    auto ta1_z_yz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 49);

    auto ta1_z_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 50);

    auto ta1_z_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 51);

    auto ta1_z_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 52);

    auto ta1_z_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 53);

    // Set up components of auxiliary buffer : DP

    auto ta1_x_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp);

    auto ta1_x_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 1);

    auto ta1_x_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 2);

    auto ta1_x_xy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 3);

    auto ta1_x_xy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 4);

    auto ta1_x_xz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 6);

    auto ta1_x_xz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 8);

    auto ta1_x_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 9);

    auto ta1_x_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 10);

    auto ta1_x_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 11);

    auto ta1_x_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 14);

    auto ta1_x_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 15);

    auto ta1_x_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 16);

    auto ta1_x_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 17);

    auto ta1_y_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 18);

    auto ta1_y_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 19);

    auto ta1_y_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 20);

    auto ta1_y_xy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 21);

    auto ta1_y_xy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 22);

    auto ta1_y_xz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 26);

    auto ta1_y_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 27);

    auto ta1_y_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 28);

    auto ta1_y_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 29);

    auto ta1_y_yz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 31);

    auto ta1_y_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 32);

    auto ta1_y_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 33);

    auto ta1_y_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 34);

    auto ta1_y_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 35);

    auto ta1_z_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 36);

    auto ta1_z_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 37);

    auto ta1_z_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 38);

    auto ta1_z_xy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 40);

    auto ta1_z_xz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 42);

    auto ta1_z_xz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 44);

    auto ta1_z_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 45);

    auto ta1_z_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 46);

    auto ta1_z_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 47);

    auto ta1_z_yz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 49);

    auto ta1_z_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 50);

    auto ta1_z_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 51);

    auto ta1_z_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 52);

    auto ta1_z_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 53);

    // Set up 0-3 components of targeted buffer : FP

    auto ta1_x_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp);

    auto ta1_x_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 1);

    auto ta1_x_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 2);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_x_x_x_0,   \
                             ta1_x_x_x_1,   \
                             ta1_x_x_y_0,   \
                             ta1_x_x_y_1,   \
                             ta1_x_x_z_0,   \
                             ta1_x_x_z_1,   \
                             ta1_x_xx_0_0,  \
                             ta1_x_xx_0_1,  \
                             ta1_x_xx_x_0,  \
                             ta1_x_xx_x_1,  \
                             ta1_x_xx_y_0,  \
                             ta1_x_xx_y_1,  \
                             ta1_x_xx_z_0,  \
                             ta1_x_xx_z_1,  \
                             ta1_x_xxx_x_0, \
                             ta1_x_xxx_y_0, \
                             ta1_x_xxx_z_0, \
                             ta_xx_x_1,     \
                             ta_xx_y_1,     \
                             ta_xx_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_x_0[i] = 2.0 * ta1_x_x_x_0[i] * fe_0 - 2.0 * ta1_x_x_x_1[i] * fe_0 + ta1_x_xx_0_0[i] * fe_0 - ta1_x_xx_0_1[i] * fe_0 +
                           ta_xx_x_1[i] + ta1_x_xx_x_0[i] * pa_x[i] - ta1_x_xx_x_1[i] * pc_x[i];

        ta1_x_xxx_y_0[i] =
            2.0 * ta1_x_x_y_0[i] * fe_0 - 2.0 * ta1_x_x_y_1[i] * fe_0 + ta_xx_y_1[i] + ta1_x_xx_y_0[i] * pa_x[i] - ta1_x_xx_y_1[i] * pc_x[i];

        ta1_x_xxx_z_0[i] =
            2.0 * ta1_x_x_z_0[i] * fe_0 - 2.0 * ta1_x_x_z_1[i] * fe_0 + ta_xx_z_1[i] + ta1_x_xx_z_0[i] * pa_x[i] - ta1_x_xx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto ta1_x_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 3);

    auto ta1_x_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 4);

    auto ta1_x_xxy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 5);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_x_xx_0_0,  \
                             ta1_x_xx_0_1,  \
                             ta1_x_xx_x_0,  \
                             ta1_x_xx_x_1,  \
                             ta1_x_xx_y_0,  \
                             ta1_x_xx_y_1,  \
                             ta1_x_xx_z_0,  \
                             ta1_x_xx_z_1,  \
                             ta1_x_xxy_x_0, \
                             ta1_x_xxy_y_0, \
                             ta1_x_xxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxy_x_0[i] = ta1_x_xx_x_0[i] * pa_y[i] - ta1_x_xx_x_1[i] * pc_y[i];

        ta1_x_xxy_y_0[i] = ta1_x_xx_0_0[i] * fe_0 - ta1_x_xx_0_1[i] * fe_0 + ta1_x_xx_y_0[i] * pa_y[i] - ta1_x_xx_y_1[i] * pc_y[i];

        ta1_x_xxy_z_0[i] = ta1_x_xx_z_0[i] * pa_y[i] - ta1_x_xx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto ta1_x_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 6);

    auto ta1_x_xxz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 7);

    auto ta1_x_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 8);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_x_xx_0_0,  \
                             ta1_x_xx_0_1,  \
                             ta1_x_xx_x_0,  \
                             ta1_x_xx_x_1,  \
                             ta1_x_xx_y_0,  \
                             ta1_x_xx_y_1,  \
                             ta1_x_xx_z_0,  \
                             ta1_x_xx_z_1,  \
                             ta1_x_xxz_x_0, \
                             ta1_x_xxz_y_0, \
                             ta1_x_xxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxz_x_0[i] = ta1_x_xx_x_0[i] * pa_z[i] - ta1_x_xx_x_1[i] * pc_z[i];

        ta1_x_xxz_y_0[i] = ta1_x_xx_y_0[i] * pa_z[i] - ta1_x_xx_y_1[i] * pc_z[i];

        ta1_x_xxz_z_0[i] = ta1_x_xx_0_0[i] * fe_0 - ta1_x_xx_0_1[i] * fe_0 + ta1_x_xx_z_0[i] * pa_z[i] - ta1_x_xx_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto ta1_x_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 9);

    auto ta1_x_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 10);

    auto ta1_x_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 11);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_x_x_x_0,   \
                             ta1_x_x_x_1,   \
                             ta1_x_xy_x_0,  \
                             ta1_x_xy_x_1,  \
                             ta1_x_xyy_x_0, \
                             ta1_x_xyy_y_0, \
                             ta1_x_xyy_z_0, \
                             ta1_x_yy_y_0,  \
                             ta1_x_yy_y_1,  \
                             ta1_x_yy_z_0,  \
                             ta1_x_yy_z_1,  \
                             ta_yy_y_1,     \
                             ta_yy_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyy_x_0[i] = ta1_x_x_x_0[i] * fe_0 - ta1_x_x_x_1[i] * fe_0 + ta1_x_xy_x_0[i] * pa_y[i] - ta1_x_xy_x_1[i] * pc_y[i];

        ta1_x_xyy_y_0[i] = ta_yy_y_1[i] + ta1_x_yy_y_0[i] * pa_x[i] - ta1_x_yy_y_1[i] * pc_x[i];

        ta1_x_xyy_z_0[i] = ta_yy_z_1[i] + ta1_x_yy_z_0[i] * pa_x[i] - ta1_x_yy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto ta1_x_xyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 12);

    auto ta1_x_xyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 13);

    auto ta1_x_xyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 14);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_x_xy_y_0,  \
                             ta1_x_xy_y_1,  \
                             ta1_x_xyz_x_0, \
                             ta1_x_xyz_y_0, \
                             ta1_x_xyz_z_0, \
                             ta1_x_xz_x_0,  \
                             ta1_x_xz_x_1,  \
                             ta1_x_xz_z_0,  \
                             ta1_x_xz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyz_x_0[i] = ta1_x_xz_x_0[i] * pa_y[i] - ta1_x_xz_x_1[i] * pc_y[i];

        ta1_x_xyz_y_0[i] = ta1_x_xy_y_0[i] * pa_z[i] - ta1_x_xy_y_1[i] * pc_z[i];

        ta1_x_xyz_z_0[i] = ta1_x_xz_z_0[i] * pa_y[i] - ta1_x_xz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto ta1_x_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 15);

    auto ta1_x_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 16);

    auto ta1_x_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 17);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_x_x_x_0,   \
                             ta1_x_x_x_1,   \
                             ta1_x_xz_x_0,  \
                             ta1_x_xz_x_1,  \
                             ta1_x_xzz_x_0, \
                             ta1_x_xzz_y_0, \
                             ta1_x_xzz_z_0, \
                             ta1_x_zz_y_0,  \
                             ta1_x_zz_y_1,  \
                             ta1_x_zz_z_0,  \
                             ta1_x_zz_z_1,  \
                             ta_zz_y_1,     \
                             ta_zz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzz_x_0[i] = ta1_x_x_x_0[i] * fe_0 - ta1_x_x_x_1[i] * fe_0 + ta1_x_xz_x_0[i] * pa_z[i] - ta1_x_xz_x_1[i] * pc_z[i];

        ta1_x_xzz_y_0[i] = ta_zz_y_1[i] + ta1_x_zz_y_0[i] * pa_x[i] - ta1_x_zz_y_1[i] * pc_x[i];

        ta1_x_xzz_z_0[i] = ta_zz_z_1[i] + ta1_x_zz_z_0[i] * pa_x[i] - ta1_x_zz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto ta1_x_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 18);

    auto ta1_x_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 19);

    auto ta1_x_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 20);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_x_y_x_0,   \
                             ta1_x_y_x_1,   \
                             ta1_x_y_y_0,   \
                             ta1_x_y_y_1,   \
                             ta1_x_y_z_0,   \
                             ta1_x_y_z_1,   \
                             ta1_x_yy_0_0,  \
                             ta1_x_yy_0_1,  \
                             ta1_x_yy_x_0,  \
                             ta1_x_yy_x_1,  \
                             ta1_x_yy_y_0,  \
                             ta1_x_yy_y_1,  \
                             ta1_x_yy_z_0,  \
                             ta1_x_yy_z_1,  \
                             ta1_x_yyy_x_0, \
                             ta1_x_yyy_y_0, \
                             ta1_x_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyy_x_0[i] = 2.0 * ta1_x_y_x_0[i] * fe_0 - 2.0 * ta1_x_y_x_1[i] * fe_0 + ta1_x_yy_x_0[i] * pa_y[i] - ta1_x_yy_x_1[i] * pc_y[i];

        ta1_x_yyy_y_0[i] = 2.0 * ta1_x_y_y_0[i] * fe_0 - 2.0 * ta1_x_y_y_1[i] * fe_0 + ta1_x_yy_0_0[i] * fe_0 - ta1_x_yy_0_1[i] * fe_0 +
                           ta1_x_yy_y_0[i] * pa_y[i] - ta1_x_yy_y_1[i] * pc_y[i];

        ta1_x_yyy_z_0[i] = 2.0 * ta1_x_y_z_0[i] * fe_0 - 2.0 * ta1_x_y_z_1[i] * fe_0 + ta1_x_yy_z_0[i] * pa_y[i] - ta1_x_yy_z_1[i] * pc_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto ta1_x_yyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 21);

    auto ta1_x_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 22);

    auto ta1_x_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 23);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_x_yy_x_0,  \
                             ta1_x_yy_x_1,  \
                             ta1_x_yy_y_0,  \
                             ta1_x_yy_y_1,  \
                             ta1_x_yyz_x_0, \
                             ta1_x_yyz_y_0, \
                             ta1_x_yyz_z_0, \
                             ta1_x_yz_z_0,  \
                             ta1_x_yz_z_1,  \
                             ta1_x_z_z_0,   \
                             ta1_x_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyz_x_0[i] = ta1_x_yy_x_0[i] * pa_z[i] - ta1_x_yy_x_1[i] * pc_z[i];

        ta1_x_yyz_y_0[i] = ta1_x_yy_y_0[i] * pa_z[i] - ta1_x_yy_y_1[i] * pc_z[i];

        ta1_x_yyz_z_0[i] = ta1_x_z_z_0[i] * fe_0 - ta1_x_z_z_1[i] * fe_0 + ta1_x_yz_z_0[i] * pa_y[i] - ta1_x_yz_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto ta1_x_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 24);

    auto ta1_x_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 25);

    auto ta1_x_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 26);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_x_yzz_x_0, \
                             ta1_x_yzz_y_0, \
                             ta1_x_yzz_z_0, \
                             ta1_x_zz_0_0,  \
                             ta1_x_zz_0_1,  \
                             ta1_x_zz_x_0,  \
                             ta1_x_zz_x_1,  \
                             ta1_x_zz_y_0,  \
                             ta1_x_zz_y_1,  \
                             ta1_x_zz_z_0,  \
                             ta1_x_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzz_x_0[i] = ta1_x_zz_x_0[i] * pa_y[i] - ta1_x_zz_x_1[i] * pc_y[i];

        ta1_x_yzz_y_0[i] = ta1_x_zz_0_0[i] * fe_0 - ta1_x_zz_0_1[i] * fe_0 + ta1_x_zz_y_0[i] * pa_y[i] - ta1_x_zz_y_1[i] * pc_y[i];

        ta1_x_yzz_z_0[i] = ta1_x_zz_z_0[i] * pa_y[i] - ta1_x_zz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto ta1_x_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 27);

    auto ta1_x_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 28);

    auto ta1_x_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 29);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_x_z_x_0,   \
                             ta1_x_z_x_1,   \
                             ta1_x_z_y_0,   \
                             ta1_x_z_y_1,   \
                             ta1_x_z_z_0,   \
                             ta1_x_z_z_1,   \
                             ta1_x_zz_0_0,  \
                             ta1_x_zz_0_1,  \
                             ta1_x_zz_x_0,  \
                             ta1_x_zz_x_1,  \
                             ta1_x_zz_y_0,  \
                             ta1_x_zz_y_1,  \
                             ta1_x_zz_z_0,  \
                             ta1_x_zz_z_1,  \
                             ta1_x_zzz_x_0, \
                             ta1_x_zzz_y_0, \
                             ta1_x_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzz_x_0[i] = 2.0 * ta1_x_z_x_0[i] * fe_0 - 2.0 * ta1_x_z_x_1[i] * fe_0 + ta1_x_zz_x_0[i] * pa_z[i] - ta1_x_zz_x_1[i] * pc_z[i];

        ta1_x_zzz_y_0[i] = 2.0 * ta1_x_z_y_0[i] * fe_0 - 2.0 * ta1_x_z_y_1[i] * fe_0 + ta1_x_zz_y_0[i] * pa_z[i] - ta1_x_zz_y_1[i] * pc_z[i];

        ta1_x_zzz_z_0[i] = 2.0 * ta1_x_z_z_0[i] * fe_0 - 2.0 * ta1_x_z_z_1[i] * fe_0 + ta1_x_zz_0_0[i] * fe_0 - ta1_x_zz_0_1[i] * fe_0 +
                           ta1_x_zz_z_0[i] * pa_z[i] - ta1_x_zz_z_1[i] * pc_z[i];
    }

    // Set up 30-33 components of targeted buffer : FP

    auto ta1_y_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 30);

    auto ta1_y_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 31);

    auto ta1_y_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 32);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_y_x_x_0,   \
                             ta1_y_x_x_1,   \
                             ta1_y_x_y_0,   \
                             ta1_y_x_y_1,   \
                             ta1_y_x_z_0,   \
                             ta1_y_x_z_1,   \
                             ta1_y_xx_0_0,  \
                             ta1_y_xx_0_1,  \
                             ta1_y_xx_x_0,  \
                             ta1_y_xx_x_1,  \
                             ta1_y_xx_y_0,  \
                             ta1_y_xx_y_1,  \
                             ta1_y_xx_z_0,  \
                             ta1_y_xx_z_1,  \
                             ta1_y_xxx_x_0, \
                             ta1_y_xxx_y_0, \
                             ta1_y_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxx_x_0[i] = 2.0 * ta1_y_x_x_0[i] * fe_0 - 2.0 * ta1_y_x_x_1[i] * fe_0 + ta1_y_xx_0_0[i] * fe_0 - ta1_y_xx_0_1[i] * fe_0 +
                           ta1_y_xx_x_0[i] * pa_x[i] - ta1_y_xx_x_1[i] * pc_x[i];

        ta1_y_xxx_y_0[i] = 2.0 * ta1_y_x_y_0[i] * fe_0 - 2.0 * ta1_y_x_y_1[i] * fe_0 + ta1_y_xx_y_0[i] * pa_x[i] - ta1_y_xx_y_1[i] * pc_x[i];

        ta1_y_xxx_z_0[i] = 2.0 * ta1_y_x_z_0[i] * fe_0 - 2.0 * ta1_y_x_z_1[i] * fe_0 + ta1_y_xx_z_0[i] * pa_x[i] - ta1_y_xx_z_1[i] * pc_x[i];
    }

    // Set up 33-36 components of targeted buffer : FP

    auto ta1_y_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 33);

    auto ta1_y_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 34);

    auto ta1_y_xxy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 35);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_y_xx_x_0,  \
                             ta1_y_xx_x_1,  \
                             ta1_y_xx_z_0,  \
                             ta1_y_xx_z_1,  \
                             ta1_y_xxy_x_0, \
                             ta1_y_xxy_y_0, \
                             ta1_y_xxy_z_0, \
                             ta1_y_xy_y_0,  \
                             ta1_y_xy_y_1,  \
                             ta1_y_y_y_0,   \
                             ta1_y_y_y_1,   \
                             ta_xx_x_1,     \
                             ta_xx_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxy_x_0[i] = ta_xx_x_1[i] + ta1_y_xx_x_0[i] * pa_y[i] - ta1_y_xx_x_1[i] * pc_y[i];

        ta1_y_xxy_y_0[i] = ta1_y_y_y_0[i] * fe_0 - ta1_y_y_y_1[i] * fe_0 + ta1_y_xy_y_0[i] * pa_x[i] - ta1_y_xy_y_1[i] * pc_x[i];

        ta1_y_xxy_z_0[i] = ta_xx_z_1[i] + ta1_y_xx_z_0[i] * pa_y[i] - ta1_y_xx_z_1[i] * pc_y[i];
    }

    // Set up 36-39 components of targeted buffer : FP

    auto ta1_y_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 36);

    auto ta1_y_xxz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 37);

    auto ta1_y_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 38);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_y_xx_x_0,  \
                             ta1_y_xx_x_1,  \
                             ta1_y_xx_y_0,  \
                             ta1_y_xx_y_1,  \
                             ta1_y_xxz_x_0, \
                             ta1_y_xxz_y_0, \
                             ta1_y_xxz_z_0, \
                             ta1_y_xz_z_0,  \
                             ta1_y_xz_z_1,  \
                             ta1_y_z_z_0,   \
                             ta1_y_z_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxz_x_0[i] = ta1_y_xx_x_0[i] * pa_z[i] - ta1_y_xx_x_1[i] * pc_z[i];

        ta1_y_xxz_y_0[i] = ta1_y_xx_y_0[i] * pa_z[i] - ta1_y_xx_y_1[i] * pc_z[i];

        ta1_y_xxz_z_0[i] = ta1_y_z_z_0[i] * fe_0 - ta1_y_z_z_1[i] * fe_0 + ta1_y_xz_z_0[i] * pa_x[i] - ta1_y_xz_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : FP

    auto ta1_y_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 39);

    auto ta1_y_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 40);

    auto ta1_y_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 41);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_y_xyy_x_0, \
                             ta1_y_xyy_y_0, \
                             ta1_y_xyy_z_0, \
                             ta1_y_yy_0_0,  \
                             ta1_y_yy_0_1,  \
                             ta1_y_yy_x_0,  \
                             ta1_y_yy_x_1,  \
                             ta1_y_yy_y_0,  \
                             ta1_y_yy_y_1,  \
                             ta1_y_yy_z_0,  \
                             ta1_y_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyy_x_0[i] = ta1_y_yy_0_0[i] * fe_0 - ta1_y_yy_0_1[i] * fe_0 + ta1_y_yy_x_0[i] * pa_x[i] - ta1_y_yy_x_1[i] * pc_x[i];

        ta1_y_xyy_y_0[i] = ta1_y_yy_y_0[i] * pa_x[i] - ta1_y_yy_y_1[i] * pc_x[i];

        ta1_y_xyy_z_0[i] = ta1_y_yy_z_0[i] * pa_x[i] - ta1_y_yy_z_1[i] * pc_x[i];
    }

    // Set up 42-45 components of targeted buffer : FP

    auto ta1_y_xyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 42);

    auto ta1_y_xyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 43);

    auto ta1_y_xyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 44);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_y_xy_x_0,  \
                             ta1_y_xy_x_1,  \
                             ta1_y_xyz_x_0, \
                             ta1_y_xyz_y_0, \
                             ta1_y_xyz_z_0, \
                             ta1_y_yz_y_0,  \
                             ta1_y_yz_y_1,  \
                             ta1_y_yz_z_0,  \
                             ta1_y_yz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyz_x_0[i] = ta1_y_xy_x_0[i] * pa_z[i] - ta1_y_xy_x_1[i] * pc_z[i];

        ta1_y_xyz_y_0[i] = ta1_y_yz_y_0[i] * pa_x[i] - ta1_y_yz_y_1[i] * pc_x[i];

        ta1_y_xyz_z_0[i] = ta1_y_yz_z_0[i] * pa_x[i] - ta1_y_yz_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : FP

    auto ta1_y_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 45);

    auto ta1_y_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 46);

    auto ta1_y_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 47);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_y_xzz_x_0, \
                             ta1_y_xzz_y_0, \
                             ta1_y_xzz_z_0, \
                             ta1_y_zz_0_0,  \
                             ta1_y_zz_0_1,  \
                             ta1_y_zz_x_0,  \
                             ta1_y_zz_x_1,  \
                             ta1_y_zz_y_0,  \
                             ta1_y_zz_y_1,  \
                             ta1_y_zz_z_0,  \
                             ta1_y_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzz_x_0[i] = ta1_y_zz_0_0[i] * fe_0 - ta1_y_zz_0_1[i] * fe_0 + ta1_y_zz_x_0[i] * pa_x[i] - ta1_y_zz_x_1[i] * pc_x[i];

        ta1_y_xzz_y_0[i] = ta1_y_zz_y_0[i] * pa_x[i] - ta1_y_zz_y_1[i] * pc_x[i];

        ta1_y_xzz_z_0[i] = ta1_y_zz_z_0[i] * pa_x[i] - ta1_y_zz_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : FP

    auto ta1_y_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 48);

    auto ta1_y_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 49);

    auto ta1_y_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 50);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_y_y_x_0,   \
                             ta1_y_y_x_1,   \
                             ta1_y_y_y_0,   \
                             ta1_y_y_y_1,   \
                             ta1_y_y_z_0,   \
                             ta1_y_y_z_1,   \
                             ta1_y_yy_0_0,  \
                             ta1_y_yy_0_1,  \
                             ta1_y_yy_x_0,  \
                             ta1_y_yy_x_1,  \
                             ta1_y_yy_y_0,  \
                             ta1_y_yy_y_1,  \
                             ta1_y_yy_z_0,  \
                             ta1_y_yy_z_1,  \
                             ta1_y_yyy_x_0, \
                             ta1_y_yyy_y_0, \
                             ta1_y_yyy_z_0, \
                             ta_yy_x_1,     \
                             ta_yy_y_1,     \
                             ta_yy_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyy_x_0[i] =
            2.0 * ta1_y_y_x_0[i] * fe_0 - 2.0 * ta1_y_y_x_1[i] * fe_0 + ta_yy_x_1[i] + ta1_y_yy_x_0[i] * pa_y[i] - ta1_y_yy_x_1[i] * pc_y[i];

        ta1_y_yyy_y_0[i] = 2.0 * ta1_y_y_y_0[i] * fe_0 - 2.0 * ta1_y_y_y_1[i] * fe_0 + ta1_y_yy_0_0[i] * fe_0 - ta1_y_yy_0_1[i] * fe_0 +
                           ta_yy_y_1[i] + ta1_y_yy_y_0[i] * pa_y[i] - ta1_y_yy_y_1[i] * pc_y[i];

        ta1_y_yyy_z_0[i] =
            2.0 * ta1_y_y_z_0[i] * fe_0 - 2.0 * ta1_y_y_z_1[i] * fe_0 + ta_yy_z_1[i] + ta1_y_yy_z_0[i] * pa_y[i] - ta1_y_yy_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : FP

    auto ta1_y_yyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 51);

    auto ta1_y_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 52);

    auto ta1_y_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 53);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_y_yy_0_0,  \
                             ta1_y_yy_0_1,  \
                             ta1_y_yy_x_0,  \
                             ta1_y_yy_x_1,  \
                             ta1_y_yy_y_0,  \
                             ta1_y_yy_y_1,  \
                             ta1_y_yy_z_0,  \
                             ta1_y_yy_z_1,  \
                             ta1_y_yyz_x_0, \
                             ta1_y_yyz_y_0, \
                             ta1_y_yyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyz_x_0[i] = ta1_y_yy_x_0[i] * pa_z[i] - ta1_y_yy_x_1[i] * pc_z[i];

        ta1_y_yyz_y_0[i] = ta1_y_yy_y_0[i] * pa_z[i] - ta1_y_yy_y_1[i] * pc_z[i];

        ta1_y_yyz_z_0[i] = ta1_y_yy_0_0[i] * fe_0 - ta1_y_yy_0_1[i] * fe_0 + ta1_y_yy_z_0[i] * pa_z[i] - ta1_y_yy_z_1[i] * pc_z[i];
    }

    // Set up 54-57 components of targeted buffer : FP

    auto ta1_y_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 54);

    auto ta1_y_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 55);

    auto ta1_y_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 56);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_y_y_y_0,   \
                             ta1_y_y_y_1,   \
                             ta1_y_yz_y_0,  \
                             ta1_y_yz_y_1,  \
                             ta1_y_yzz_x_0, \
                             ta1_y_yzz_y_0, \
                             ta1_y_yzz_z_0, \
                             ta1_y_zz_x_0,  \
                             ta1_y_zz_x_1,  \
                             ta1_y_zz_z_0,  \
                             ta1_y_zz_z_1,  \
                             ta_zz_x_1,     \
                             ta_zz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzz_x_0[i] = ta_zz_x_1[i] + ta1_y_zz_x_0[i] * pa_y[i] - ta1_y_zz_x_1[i] * pc_y[i];

        ta1_y_yzz_y_0[i] = ta1_y_y_y_0[i] * fe_0 - ta1_y_y_y_1[i] * fe_0 + ta1_y_yz_y_0[i] * pa_z[i] - ta1_y_yz_y_1[i] * pc_z[i];

        ta1_y_yzz_z_0[i] = ta_zz_z_1[i] + ta1_y_zz_z_0[i] * pa_y[i] - ta1_y_zz_z_1[i] * pc_y[i];
    }

    // Set up 57-60 components of targeted buffer : FP

    auto ta1_y_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 57);

    auto ta1_y_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 58);

    auto ta1_y_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 59);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_y_z_x_0,   \
                             ta1_y_z_x_1,   \
                             ta1_y_z_y_0,   \
                             ta1_y_z_y_1,   \
                             ta1_y_z_z_0,   \
                             ta1_y_z_z_1,   \
                             ta1_y_zz_0_0,  \
                             ta1_y_zz_0_1,  \
                             ta1_y_zz_x_0,  \
                             ta1_y_zz_x_1,  \
                             ta1_y_zz_y_0,  \
                             ta1_y_zz_y_1,  \
                             ta1_y_zz_z_0,  \
                             ta1_y_zz_z_1,  \
                             ta1_y_zzz_x_0, \
                             ta1_y_zzz_y_0, \
                             ta1_y_zzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzz_x_0[i] = 2.0 * ta1_y_z_x_0[i] * fe_0 - 2.0 * ta1_y_z_x_1[i] * fe_0 + ta1_y_zz_x_0[i] * pa_z[i] - ta1_y_zz_x_1[i] * pc_z[i];

        ta1_y_zzz_y_0[i] = 2.0 * ta1_y_z_y_0[i] * fe_0 - 2.0 * ta1_y_z_y_1[i] * fe_0 + ta1_y_zz_y_0[i] * pa_z[i] - ta1_y_zz_y_1[i] * pc_z[i];

        ta1_y_zzz_z_0[i] = 2.0 * ta1_y_z_z_0[i] * fe_0 - 2.0 * ta1_y_z_z_1[i] * fe_0 + ta1_y_zz_0_0[i] * fe_0 - ta1_y_zz_0_1[i] * fe_0 +
                           ta1_y_zz_z_0[i] * pa_z[i] - ta1_y_zz_z_1[i] * pc_z[i];
    }

    // Set up 60-63 components of targeted buffer : FP

    auto ta1_z_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 60);

    auto ta1_z_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 61);

    auto ta1_z_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 62);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_z_x_x_0,   \
                             ta1_z_x_x_1,   \
                             ta1_z_x_y_0,   \
                             ta1_z_x_y_1,   \
                             ta1_z_x_z_0,   \
                             ta1_z_x_z_1,   \
                             ta1_z_xx_0_0,  \
                             ta1_z_xx_0_1,  \
                             ta1_z_xx_x_0,  \
                             ta1_z_xx_x_1,  \
                             ta1_z_xx_y_0,  \
                             ta1_z_xx_y_1,  \
                             ta1_z_xx_z_0,  \
                             ta1_z_xx_z_1,  \
                             ta1_z_xxx_x_0, \
                             ta1_z_xxx_y_0, \
                             ta1_z_xxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxx_x_0[i] = 2.0 * ta1_z_x_x_0[i] * fe_0 - 2.0 * ta1_z_x_x_1[i] * fe_0 + ta1_z_xx_0_0[i] * fe_0 - ta1_z_xx_0_1[i] * fe_0 +
                           ta1_z_xx_x_0[i] * pa_x[i] - ta1_z_xx_x_1[i] * pc_x[i];

        ta1_z_xxx_y_0[i] = 2.0 * ta1_z_x_y_0[i] * fe_0 - 2.0 * ta1_z_x_y_1[i] * fe_0 + ta1_z_xx_y_0[i] * pa_x[i] - ta1_z_xx_y_1[i] * pc_x[i];

        ta1_z_xxx_z_0[i] = 2.0 * ta1_z_x_z_0[i] * fe_0 - 2.0 * ta1_z_x_z_1[i] * fe_0 + ta1_z_xx_z_0[i] * pa_x[i] - ta1_z_xx_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : FP

    auto ta1_z_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 63);

    auto ta1_z_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 64);

    auto ta1_z_xxy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 65);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_z_xx_x_0,  \
                             ta1_z_xx_x_1,  \
                             ta1_z_xx_z_0,  \
                             ta1_z_xx_z_1,  \
                             ta1_z_xxy_x_0, \
                             ta1_z_xxy_y_0, \
                             ta1_z_xxy_z_0, \
                             ta1_z_xy_y_0,  \
                             ta1_z_xy_y_1,  \
                             ta1_z_y_y_0,   \
                             ta1_z_y_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxy_x_0[i] = ta1_z_xx_x_0[i] * pa_y[i] - ta1_z_xx_x_1[i] * pc_y[i];

        ta1_z_xxy_y_0[i] = ta1_z_y_y_0[i] * fe_0 - ta1_z_y_y_1[i] * fe_0 + ta1_z_xy_y_0[i] * pa_x[i] - ta1_z_xy_y_1[i] * pc_x[i];

        ta1_z_xxy_z_0[i] = ta1_z_xx_z_0[i] * pa_y[i] - ta1_z_xx_z_1[i] * pc_y[i];
    }

    // Set up 66-69 components of targeted buffer : FP

    auto ta1_z_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 66);

    auto ta1_z_xxz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 67);

    auto ta1_z_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 68);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta1_z_xx_x_0,  \
                             ta1_z_xx_x_1,  \
                             ta1_z_xx_y_0,  \
                             ta1_z_xx_y_1,  \
                             ta1_z_xxz_x_0, \
                             ta1_z_xxz_y_0, \
                             ta1_z_xxz_z_0, \
                             ta1_z_xz_z_0,  \
                             ta1_z_xz_z_1,  \
                             ta1_z_z_z_0,   \
                             ta1_z_z_z_1,   \
                             ta_xx_x_1,     \
                             ta_xx_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxz_x_0[i] = ta_xx_x_1[i] + ta1_z_xx_x_0[i] * pa_z[i] - ta1_z_xx_x_1[i] * pc_z[i];

        ta1_z_xxz_y_0[i] = ta_xx_y_1[i] + ta1_z_xx_y_0[i] * pa_z[i] - ta1_z_xx_y_1[i] * pc_z[i];

        ta1_z_xxz_z_0[i] = ta1_z_z_z_0[i] * fe_0 - ta1_z_z_z_1[i] * fe_0 + ta1_z_xz_z_0[i] * pa_x[i] - ta1_z_xz_z_1[i] * pc_x[i];
    }

    // Set up 69-72 components of targeted buffer : FP

    auto ta1_z_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 69);

    auto ta1_z_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 70);

    auto ta1_z_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 71);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_z_xyy_x_0, \
                             ta1_z_xyy_y_0, \
                             ta1_z_xyy_z_0, \
                             ta1_z_yy_0_0,  \
                             ta1_z_yy_0_1,  \
                             ta1_z_yy_x_0,  \
                             ta1_z_yy_x_1,  \
                             ta1_z_yy_y_0,  \
                             ta1_z_yy_y_1,  \
                             ta1_z_yy_z_0,  \
                             ta1_z_yy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyy_x_0[i] = ta1_z_yy_0_0[i] * fe_0 - ta1_z_yy_0_1[i] * fe_0 + ta1_z_yy_x_0[i] * pa_x[i] - ta1_z_yy_x_1[i] * pc_x[i];

        ta1_z_xyy_y_0[i] = ta1_z_yy_y_0[i] * pa_x[i] - ta1_z_yy_y_1[i] * pc_x[i];

        ta1_z_xyy_z_0[i] = ta1_z_yy_z_0[i] * pa_x[i] - ta1_z_yy_z_1[i] * pc_x[i];
    }

    // Set up 72-75 components of targeted buffer : FP

    auto ta1_z_xyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 72);

    auto ta1_z_xyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 73);

    auto ta1_z_xyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 74);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta1_z_xyz_x_0, \
                             ta1_z_xyz_y_0, \
                             ta1_z_xyz_z_0, \
                             ta1_z_xz_x_0,  \
                             ta1_z_xz_x_1,  \
                             ta1_z_yz_y_0,  \
                             ta1_z_yz_y_1,  \
                             ta1_z_yz_z_0,  \
                             ta1_z_yz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyz_x_0[i] = ta1_z_xz_x_0[i] * pa_y[i] - ta1_z_xz_x_1[i] * pc_y[i];

        ta1_z_xyz_y_0[i] = ta1_z_yz_y_0[i] * pa_x[i] - ta1_z_yz_y_1[i] * pc_x[i];

        ta1_z_xyz_z_0[i] = ta1_z_yz_z_0[i] * pa_x[i] - ta1_z_yz_z_1[i] * pc_x[i];
    }

    // Set up 75-78 components of targeted buffer : FP

    auto ta1_z_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 75);

    auto ta1_z_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 76);

    auto ta1_z_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 77);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_z_xzz_x_0, \
                             ta1_z_xzz_y_0, \
                             ta1_z_xzz_z_0, \
                             ta1_z_zz_0_0,  \
                             ta1_z_zz_0_1,  \
                             ta1_z_zz_x_0,  \
                             ta1_z_zz_x_1,  \
                             ta1_z_zz_y_0,  \
                             ta1_z_zz_y_1,  \
                             ta1_z_zz_z_0,  \
                             ta1_z_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzz_x_0[i] = ta1_z_zz_0_0[i] * fe_0 - ta1_z_zz_0_1[i] * fe_0 + ta1_z_zz_x_0[i] * pa_x[i] - ta1_z_zz_x_1[i] * pc_x[i];

        ta1_z_xzz_y_0[i] = ta1_z_zz_y_0[i] * pa_x[i] - ta1_z_zz_y_1[i] * pc_x[i];

        ta1_z_xzz_z_0[i] = ta1_z_zz_z_0[i] * pa_x[i] - ta1_z_zz_z_1[i] * pc_x[i];
    }

    // Set up 78-81 components of targeted buffer : FP

    auto ta1_z_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 78);

    auto ta1_z_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 79);

    auto ta1_z_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 80);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_z_y_x_0,   \
                             ta1_z_y_x_1,   \
                             ta1_z_y_y_0,   \
                             ta1_z_y_y_1,   \
                             ta1_z_y_z_0,   \
                             ta1_z_y_z_1,   \
                             ta1_z_yy_0_0,  \
                             ta1_z_yy_0_1,  \
                             ta1_z_yy_x_0,  \
                             ta1_z_yy_x_1,  \
                             ta1_z_yy_y_0,  \
                             ta1_z_yy_y_1,  \
                             ta1_z_yy_z_0,  \
                             ta1_z_yy_z_1,  \
                             ta1_z_yyy_x_0, \
                             ta1_z_yyy_y_0, \
                             ta1_z_yyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyy_x_0[i] = 2.0 * ta1_z_y_x_0[i] * fe_0 - 2.0 * ta1_z_y_x_1[i] * fe_0 + ta1_z_yy_x_0[i] * pa_y[i] - ta1_z_yy_x_1[i] * pc_y[i];

        ta1_z_yyy_y_0[i] = 2.0 * ta1_z_y_y_0[i] * fe_0 - 2.0 * ta1_z_y_y_1[i] * fe_0 + ta1_z_yy_0_0[i] * fe_0 - ta1_z_yy_0_1[i] * fe_0 +
                           ta1_z_yy_y_0[i] * pa_y[i] - ta1_z_yy_y_1[i] * pc_y[i];

        ta1_z_yyy_z_0[i] = 2.0 * ta1_z_y_z_0[i] * fe_0 - 2.0 * ta1_z_y_z_1[i] * fe_0 + ta1_z_yy_z_0[i] * pa_y[i] - ta1_z_yy_z_1[i] * pc_y[i];
    }

    // Set up 81-84 components of targeted buffer : FP

    auto ta1_z_yyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 81);

    auto ta1_z_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 82);

    auto ta1_z_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 83);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_z_yy_x_0,  \
                             ta1_z_yy_x_1,  \
                             ta1_z_yy_y_0,  \
                             ta1_z_yy_y_1,  \
                             ta1_z_yyz_x_0, \
                             ta1_z_yyz_y_0, \
                             ta1_z_yyz_z_0, \
                             ta1_z_yz_z_0,  \
                             ta1_z_yz_z_1,  \
                             ta1_z_z_z_0,   \
                             ta1_z_z_z_1,   \
                             ta_yy_x_1,     \
                             ta_yy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyz_x_0[i] = ta_yy_x_1[i] + ta1_z_yy_x_0[i] * pa_z[i] - ta1_z_yy_x_1[i] * pc_z[i];

        ta1_z_yyz_y_0[i] = ta_yy_y_1[i] + ta1_z_yy_y_0[i] * pa_z[i] - ta1_z_yy_y_1[i] * pc_z[i];

        ta1_z_yyz_z_0[i] = ta1_z_z_z_0[i] * fe_0 - ta1_z_z_z_1[i] * fe_0 + ta1_z_yz_z_0[i] * pa_y[i] - ta1_z_yz_z_1[i] * pc_y[i];
    }

    // Set up 84-87 components of targeted buffer : FP

    auto ta1_z_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 84);

    auto ta1_z_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 85);

    auto ta1_z_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 86);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_z_yzz_x_0, \
                             ta1_z_yzz_y_0, \
                             ta1_z_yzz_z_0, \
                             ta1_z_zz_0_0,  \
                             ta1_z_zz_0_1,  \
                             ta1_z_zz_x_0,  \
                             ta1_z_zz_x_1,  \
                             ta1_z_zz_y_0,  \
                             ta1_z_zz_y_1,  \
                             ta1_z_zz_z_0,  \
                             ta1_z_zz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzz_x_0[i] = ta1_z_zz_x_0[i] * pa_y[i] - ta1_z_zz_x_1[i] * pc_y[i];

        ta1_z_yzz_y_0[i] = ta1_z_zz_0_0[i] * fe_0 - ta1_z_zz_0_1[i] * fe_0 + ta1_z_zz_y_0[i] * pa_y[i] - ta1_z_zz_y_1[i] * pc_y[i];

        ta1_z_yzz_z_0[i] = ta1_z_zz_z_0[i] * pa_y[i] - ta1_z_zz_z_1[i] * pc_y[i];
    }

    // Set up 87-90 components of targeted buffer : FP

    auto ta1_z_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 87);

    auto ta1_z_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 88);

    auto ta1_z_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 89);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_z_z_x_0,   \
                             ta1_z_z_x_1,   \
                             ta1_z_z_y_0,   \
                             ta1_z_z_y_1,   \
                             ta1_z_z_z_0,   \
                             ta1_z_z_z_1,   \
                             ta1_z_zz_0_0,  \
                             ta1_z_zz_0_1,  \
                             ta1_z_zz_x_0,  \
                             ta1_z_zz_x_1,  \
                             ta1_z_zz_y_0,  \
                             ta1_z_zz_y_1,  \
                             ta1_z_zz_z_0,  \
                             ta1_z_zz_z_1,  \
                             ta1_z_zzz_x_0, \
                             ta1_z_zzz_y_0, \
                             ta1_z_zzz_z_0, \
                             ta_zz_x_1,     \
                             ta_zz_y_1,     \
                             ta_zz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzz_x_0[i] =
            2.0 * ta1_z_z_x_0[i] * fe_0 - 2.0 * ta1_z_z_x_1[i] * fe_0 + ta_zz_x_1[i] + ta1_z_zz_x_0[i] * pa_z[i] - ta1_z_zz_x_1[i] * pc_z[i];

        ta1_z_zzz_y_0[i] =
            2.0 * ta1_z_z_y_0[i] * fe_0 - 2.0 * ta1_z_z_y_1[i] * fe_0 + ta_zz_y_1[i] + ta1_z_zz_y_0[i] * pa_z[i] - ta1_z_zz_y_1[i] * pc_z[i];

        ta1_z_zzz_z_0[i] = 2.0 * ta1_z_z_z_0[i] * fe_0 - 2.0 * ta1_z_z_z_1[i] * fe_0 + ta1_z_zz_0_0[i] * fe_0 - ta1_z_zz_0_1[i] * fe_0 +
                           ta_zz_z_1[i] + ta1_z_zz_z_0[i] * pa_z[i] - ta1_z_zz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
