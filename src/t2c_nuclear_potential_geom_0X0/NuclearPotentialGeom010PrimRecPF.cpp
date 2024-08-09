#include "NuclearPotentialGeom010PrimRecPF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_pf(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_pf,
                                        const size_t              idx_npot_geom_010_0_sd,
                                        const size_t              idx_npot_geom_010_1_sd,
                                        const size_t              idx_npot_1_sf,
                                        const size_t              idx_npot_geom_010_0_sf,
                                        const size_t              idx_npot_geom_010_1_sf,
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

    auto ta1_x_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd);

    auto ta1_x_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 1);

    auto ta1_x_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 2);

    auto ta1_x_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 3);

    auto ta1_x_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 4);

    auto ta1_x_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 5);

    auto ta1_y_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 6);

    auto ta1_y_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 7);

    auto ta1_y_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 8);

    auto ta1_y_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 9);

    auto ta1_y_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 10);

    auto ta1_y_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 11);

    auto ta1_z_0_xx_0 = pbuffer.data(idx_npot_geom_010_0_sd + 12);

    auto ta1_z_0_xy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 13);

    auto ta1_z_0_xz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 14);

    auto ta1_z_0_yy_0 = pbuffer.data(idx_npot_geom_010_0_sd + 15);

    auto ta1_z_0_yz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 16);

    auto ta1_z_0_zz_0 = pbuffer.data(idx_npot_geom_010_0_sd + 17);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd);

    auto ta1_x_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 1);

    auto ta1_x_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 2);

    auto ta1_x_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 3);

    auto ta1_x_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 4);

    auto ta1_x_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 5);

    auto ta1_y_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 6);

    auto ta1_y_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 7);

    auto ta1_y_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 8);

    auto ta1_y_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 9);

    auto ta1_y_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 10);

    auto ta1_y_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 11);

    auto ta1_z_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 12);

    auto ta1_z_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 13);

    auto ta1_z_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 14);

    auto ta1_z_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 15);

    auto ta1_z_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 16);

    auto ta1_z_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 17);

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

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf);

    auto ta1_x_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 1);

    auto ta1_x_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 2);

    auto ta1_x_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 3);

    auto ta1_x_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 4);

    auto ta1_x_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 5);

    auto ta1_x_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 6);

    auto ta1_x_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 7);

    auto ta1_x_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 8);

    auto ta1_x_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 9);

    auto ta1_y_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 10);

    auto ta1_y_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 11);

    auto ta1_y_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 12);

    auto ta1_y_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 13);

    auto ta1_y_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 14);

    auto ta1_y_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 15);

    auto ta1_y_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 16);

    auto ta1_y_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 17);

    auto ta1_y_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 18);

    auto ta1_y_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 19);

    auto ta1_z_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 20);

    auto ta1_z_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 21);

    auto ta1_z_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 22);

    auto ta1_z_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 23);

    auto ta1_z_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 24);

    auto ta1_z_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 25);

    auto ta1_z_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 26);

    auto ta1_z_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 27);

    auto ta1_z_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 28);

    auto ta1_z_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 29);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf);

    auto ta1_x_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 1);

    auto ta1_x_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 2);

    auto ta1_x_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 3);

    auto ta1_x_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 4);

    auto ta1_x_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 5);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 7);

    auto ta1_x_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 8);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 11);

    auto ta1_y_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 12);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 14);

    auto ta1_y_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 15);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 17);

    auto ta1_y_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 18);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 21);

    auto ta1_z_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 22);

    auto ta1_z_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 23);

    auto ta1_z_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 24);

    auto ta1_z_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 25);

    auto ta1_z_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 26);

    auto ta1_z_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 27);

    auto ta1_z_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 28);

    auto ta1_z_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 29);

    // Set up 0-10 components of targeted buffer : PF

    auto ta1_x_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf);

    auto ta1_x_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 1);

    auto ta1_x_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 2);

    auto ta1_x_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 3);

    auto ta1_x_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 4);

    auto ta1_x_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 5);

    auto ta1_x_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 6);

    auto ta1_x_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 7);

    auto ta1_x_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 8);

    auto ta1_x_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 9);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_x_0_xx_0,  \
                             ta1_x_0_xx_1,  \
                             ta1_x_0_xxx_0, \
                             ta1_x_0_xxx_1, \
                             ta1_x_0_xxy_0, \
                             ta1_x_0_xxy_1, \
                             ta1_x_0_xxz_0, \
                             ta1_x_0_xxz_1, \
                             ta1_x_0_xy_0,  \
                             ta1_x_0_xy_1,  \
                             ta1_x_0_xyy_0, \
                             ta1_x_0_xyy_1, \
                             ta1_x_0_xyz_0, \
                             ta1_x_0_xyz_1, \
                             ta1_x_0_xz_0,  \
                             ta1_x_0_xz_1,  \
                             ta1_x_0_xzz_0, \
                             ta1_x_0_xzz_1, \
                             ta1_x_0_yy_0,  \
                             ta1_x_0_yy_1,  \
                             ta1_x_0_yyy_0, \
                             ta1_x_0_yyy_1, \
                             ta1_x_0_yyz_0, \
                             ta1_x_0_yyz_1, \
                             ta1_x_0_yz_0,  \
                             ta1_x_0_yz_1,  \
                             ta1_x_0_yzz_0, \
                             ta1_x_0_yzz_1, \
                             ta1_x_0_zz_0,  \
                             ta1_x_0_zz_1,  \
                             ta1_x_0_zzz_0, \
                             ta1_x_0_zzz_1, \
                             ta1_x_x_xxx_0, \
                             ta1_x_x_xxy_0, \
                             ta1_x_x_xxz_0, \
                             ta1_x_x_xyy_0, \
                             ta1_x_x_xyz_0, \
                             ta1_x_x_xzz_0, \
                             ta1_x_x_yyy_0, \
                             ta1_x_x_yyz_0, \
                             ta1_x_x_yzz_0, \
                             ta1_x_x_zzz_0, \
                             ta_0_xxx_1,    \
                             ta_0_xxy_1,    \
                             ta_0_xxz_1,    \
                             ta_0_xyy_1,    \
                             ta_0_xyz_1,    \
                             ta_0_xzz_1,    \
                             ta_0_yyy_1,    \
                             ta_0_yyz_1,    \
                             ta_0_yzz_1,    \
                             ta_0_zzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_x_xxx_0[i] = 3.0 * ta1_x_0_xx_0[i] * fe_0 - 3.0 * ta1_x_0_xx_1[i] * fe_0 + ta_0_xxx_1[i] +
                           ta1_x_0_xxx_0[i] * pa_x[i] - ta1_x_0_xxx_1[i] * pc_x[i];

        ta1_x_x_xxy_0[i] = 2.0 * ta1_x_0_xy_0[i] * fe_0 - 2.0 * ta1_x_0_xy_1[i] * fe_0 + ta_0_xxy_1[i] +
                           ta1_x_0_xxy_0[i] * pa_x[i] - ta1_x_0_xxy_1[i] * pc_x[i];

        ta1_x_x_xxz_0[i] = 2.0 * ta1_x_0_xz_0[i] * fe_0 - 2.0 * ta1_x_0_xz_1[i] * fe_0 + ta_0_xxz_1[i] +
                           ta1_x_0_xxz_0[i] * pa_x[i] - ta1_x_0_xxz_1[i] * pc_x[i];

        ta1_x_x_xyy_0[i] = ta1_x_0_yy_0[i] * fe_0 - ta1_x_0_yy_1[i] * fe_0 + ta_0_xyy_1[i] +
                           ta1_x_0_xyy_0[i] * pa_x[i] - ta1_x_0_xyy_1[i] * pc_x[i];

        ta1_x_x_xyz_0[i] = ta1_x_0_yz_0[i] * fe_0 - ta1_x_0_yz_1[i] * fe_0 + ta_0_xyz_1[i] +
                           ta1_x_0_xyz_0[i] * pa_x[i] - ta1_x_0_xyz_1[i] * pc_x[i];

        ta1_x_x_xzz_0[i] = ta1_x_0_zz_0[i] * fe_0 - ta1_x_0_zz_1[i] * fe_0 + ta_0_xzz_1[i] +
                           ta1_x_0_xzz_0[i] * pa_x[i] - ta1_x_0_xzz_1[i] * pc_x[i];

        ta1_x_x_yyy_0[i] = ta_0_yyy_1[i] + ta1_x_0_yyy_0[i] * pa_x[i] - ta1_x_0_yyy_1[i] * pc_x[i];

        ta1_x_x_yyz_0[i] = ta_0_yyz_1[i] + ta1_x_0_yyz_0[i] * pa_x[i] - ta1_x_0_yyz_1[i] * pc_x[i];

        ta1_x_x_yzz_0[i] = ta_0_yzz_1[i] + ta1_x_0_yzz_0[i] * pa_x[i] - ta1_x_0_yzz_1[i] * pc_x[i];

        ta1_x_x_zzz_0[i] = ta_0_zzz_1[i] + ta1_x_0_zzz_0[i] * pa_x[i] - ta1_x_0_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto ta1_x_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 10);

    auto ta1_x_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 11);

    auto ta1_x_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 12);

    auto ta1_x_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 13);

    auto ta1_x_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 14);

    auto ta1_x_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 15);

    auto ta1_x_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 16);

    auto ta1_x_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 17);

    auto ta1_x_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 18);

    auto ta1_x_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 19);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_x_0_xx_0,  \
                             ta1_x_0_xx_1,  \
                             ta1_x_0_xxx_0, \
                             ta1_x_0_xxx_1, \
                             ta1_x_0_xxy_0, \
                             ta1_x_0_xxy_1, \
                             ta1_x_0_xxz_0, \
                             ta1_x_0_xxz_1, \
                             ta1_x_0_xy_0,  \
                             ta1_x_0_xy_1,  \
                             ta1_x_0_xyy_0, \
                             ta1_x_0_xyy_1, \
                             ta1_x_0_xyz_0, \
                             ta1_x_0_xyz_1, \
                             ta1_x_0_xz_0,  \
                             ta1_x_0_xz_1,  \
                             ta1_x_0_xzz_0, \
                             ta1_x_0_xzz_1, \
                             ta1_x_0_yy_0,  \
                             ta1_x_0_yy_1,  \
                             ta1_x_0_yyy_0, \
                             ta1_x_0_yyy_1, \
                             ta1_x_0_yyz_0, \
                             ta1_x_0_yyz_1, \
                             ta1_x_0_yz_0,  \
                             ta1_x_0_yz_1,  \
                             ta1_x_0_yzz_0, \
                             ta1_x_0_yzz_1, \
                             ta1_x_0_zz_0,  \
                             ta1_x_0_zz_1,  \
                             ta1_x_0_zzz_0, \
                             ta1_x_0_zzz_1, \
                             ta1_x_y_xxx_0, \
                             ta1_x_y_xxy_0, \
                             ta1_x_y_xxz_0, \
                             ta1_x_y_xyy_0, \
                             ta1_x_y_xyz_0, \
                             ta1_x_y_xzz_0, \
                             ta1_x_y_yyy_0, \
                             ta1_x_y_yyz_0, \
                             ta1_x_y_yzz_0, \
                             ta1_x_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_y_xxx_0[i] = ta1_x_0_xxx_0[i] * pa_y[i] - ta1_x_0_xxx_1[i] * pc_y[i];

        ta1_x_y_xxy_0[i] =
            ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + ta1_x_0_xxy_0[i] * pa_y[i] - ta1_x_0_xxy_1[i] * pc_y[i];

        ta1_x_y_xxz_0[i] = ta1_x_0_xxz_0[i] * pa_y[i] - ta1_x_0_xxz_1[i] * pc_y[i];

        ta1_x_y_xyy_0[i] = 2.0 * ta1_x_0_xy_0[i] * fe_0 - 2.0 * ta1_x_0_xy_1[i] * fe_0 + ta1_x_0_xyy_0[i] * pa_y[i] -
                           ta1_x_0_xyy_1[i] * pc_y[i];

        ta1_x_y_xyz_0[i] =
            ta1_x_0_xz_0[i] * fe_0 - ta1_x_0_xz_1[i] * fe_0 + ta1_x_0_xyz_0[i] * pa_y[i] - ta1_x_0_xyz_1[i] * pc_y[i];

        ta1_x_y_xzz_0[i] = ta1_x_0_xzz_0[i] * pa_y[i] - ta1_x_0_xzz_1[i] * pc_y[i];

        ta1_x_y_yyy_0[i] = 3.0 * ta1_x_0_yy_0[i] * fe_0 - 3.0 * ta1_x_0_yy_1[i] * fe_0 + ta1_x_0_yyy_0[i] * pa_y[i] -
                           ta1_x_0_yyy_1[i] * pc_y[i];

        ta1_x_y_yyz_0[i] = 2.0 * ta1_x_0_yz_0[i] * fe_0 - 2.0 * ta1_x_0_yz_1[i] * fe_0 + ta1_x_0_yyz_0[i] * pa_y[i] -
                           ta1_x_0_yyz_1[i] * pc_y[i];

        ta1_x_y_yzz_0[i] =
            ta1_x_0_zz_0[i] * fe_0 - ta1_x_0_zz_1[i] * fe_0 + ta1_x_0_yzz_0[i] * pa_y[i] - ta1_x_0_yzz_1[i] * pc_y[i];

        ta1_x_y_zzz_0[i] = ta1_x_0_zzz_0[i] * pa_y[i] - ta1_x_0_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto ta1_x_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 20);

    auto ta1_x_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 21);

    auto ta1_x_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 22);

    auto ta1_x_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 23);

    auto ta1_x_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 24);

    auto ta1_x_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 25);

    auto ta1_x_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 26);

    auto ta1_x_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 27);

    auto ta1_x_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 28);

    auto ta1_x_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 29);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_x_0_xx_0,  \
                             ta1_x_0_xx_1,  \
                             ta1_x_0_xxx_0, \
                             ta1_x_0_xxx_1, \
                             ta1_x_0_xxy_0, \
                             ta1_x_0_xxy_1, \
                             ta1_x_0_xxz_0, \
                             ta1_x_0_xxz_1, \
                             ta1_x_0_xy_0,  \
                             ta1_x_0_xy_1,  \
                             ta1_x_0_xyy_0, \
                             ta1_x_0_xyy_1, \
                             ta1_x_0_xyz_0, \
                             ta1_x_0_xyz_1, \
                             ta1_x_0_xz_0,  \
                             ta1_x_0_xz_1,  \
                             ta1_x_0_xzz_0, \
                             ta1_x_0_xzz_1, \
                             ta1_x_0_yy_0,  \
                             ta1_x_0_yy_1,  \
                             ta1_x_0_yyy_0, \
                             ta1_x_0_yyy_1, \
                             ta1_x_0_yyz_0, \
                             ta1_x_0_yyz_1, \
                             ta1_x_0_yz_0,  \
                             ta1_x_0_yz_1,  \
                             ta1_x_0_yzz_0, \
                             ta1_x_0_yzz_1, \
                             ta1_x_0_zz_0,  \
                             ta1_x_0_zz_1,  \
                             ta1_x_0_zzz_0, \
                             ta1_x_0_zzz_1, \
                             ta1_x_z_xxx_0, \
                             ta1_x_z_xxy_0, \
                             ta1_x_z_xxz_0, \
                             ta1_x_z_xyy_0, \
                             ta1_x_z_xyz_0, \
                             ta1_x_z_xzz_0, \
                             ta1_x_z_yyy_0, \
                             ta1_x_z_yyz_0, \
                             ta1_x_z_yzz_0, \
                             ta1_x_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_z_xxx_0[i] = ta1_x_0_xxx_0[i] * pa_z[i] - ta1_x_0_xxx_1[i] * pc_z[i];

        ta1_x_z_xxy_0[i] = ta1_x_0_xxy_0[i] * pa_z[i] - ta1_x_0_xxy_1[i] * pc_z[i];

        ta1_x_z_xxz_0[i] =
            ta1_x_0_xx_0[i] * fe_0 - ta1_x_0_xx_1[i] * fe_0 + ta1_x_0_xxz_0[i] * pa_z[i] - ta1_x_0_xxz_1[i] * pc_z[i];

        ta1_x_z_xyy_0[i] = ta1_x_0_xyy_0[i] * pa_z[i] - ta1_x_0_xyy_1[i] * pc_z[i];

        ta1_x_z_xyz_0[i] =
            ta1_x_0_xy_0[i] * fe_0 - ta1_x_0_xy_1[i] * fe_0 + ta1_x_0_xyz_0[i] * pa_z[i] - ta1_x_0_xyz_1[i] * pc_z[i];

        ta1_x_z_xzz_0[i] = 2.0 * ta1_x_0_xz_0[i] * fe_0 - 2.0 * ta1_x_0_xz_1[i] * fe_0 + ta1_x_0_xzz_0[i] * pa_z[i] -
                           ta1_x_0_xzz_1[i] * pc_z[i];

        ta1_x_z_yyy_0[i] = ta1_x_0_yyy_0[i] * pa_z[i] - ta1_x_0_yyy_1[i] * pc_z[i];

        ta1_x_z_yyz_0[i] =
            ta1_x_0_yy_0[i] * fe_0 - ta1_x_0_yy_1[i] * fe_0 + ta1_x_0_yyz_0[i] * pa_z[i] - ta1_x_0_yyz_1[i] * pc_z[i];

        ta1_x_z_yzz_0[i] = 2.0 * ta1_x_0_yz_0[i] * fe_0 - 2.0 * ta1_x_0_yz_1[i] * fe_0 + ta1_x_0_yzz_0[i] * pa_z[i] -
                           ta1_x_0_yzz_1[i] * pc_z[i];

        ta1_x_z_zzz_0[i] = 3.0 * ta1_x_0_zz_0[i] * fe_0 - 3.0 * ta1_x_0_zz_1[i] * fe_0 + ta1_x_0_zzz_0[i] * pa_z[i] -
                           ta1_x_0_zzz_1[i] * pc_z[i];
    }

    // Set up 30-40 components of targeted buffer : PF

    auto ta1_y_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 30);

    auto ta1_y_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 31);

    auto ta1_y_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 32);

    auto ta1_y_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 33);

    auto ta1_y_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 34);

    auto ta1_y_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 35);

    auto ta1_y_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 36);

    auto ta1_y_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 37);

    auto ta1_y_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 38);

    auto ta1_y_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 39);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_y_0_xx_0,  \
                             ta1_y_0_xx_1,  \
                             ta1_y_0_xxx_0, \
                             ta1_y_0_xxx_1, \
                             ta1_y_0_xxy_0, \
                             ta1_y_0_xxy_1, \
                             ta1_y_0_xxz_0, \
                             ta1_y_0_xxz_1, \
                             ta1_y_0_xy_0,  \
                             ta1_y_0_xy_1,  \
                             ta1_y_0_xyy_0, \
                             ta1_y_0_xyy_1, \
                             ta1_y_0_xyz_0, \
                             ta1_y_0_xyz_1, \
                             ta1_y_0_xz_0,  \
                             ta1_y_0_xz_1,  \
                             ta1_y_0_xzz_0, \
                             ta1_y_0_xzz_1, \
                             ta1_y_0_yy_0,  \
                             ta1_y_0_yy_1,  \
                             ta1_y_0_yyy_0, \
                             ta1_y_0_yyy_1, \
                             ta1_y_0_yyz_0, \
                             ta1_y_0_yyz_1, \
                             ta1_y_0_yz_0,  \
                             ta1_y_0_yz_1,  \
                             ta1_y_0_yzz_0, \
                             ta1_y_0_yzz_1, \
                             ta1_y_0_zz_0,  \
                             ta1_y_0_zz_1,  \
                             ta1_y_0_zzz_0, \
                             ta1_y_0_zzz_1, \
                             ta1_y_x_xxx_0, \
                             ta1_y_x_xxy_0, \
                             ta1_y_x_xxz_0, \
                             ta1_y_x_xyy_0, \
                             ta1_y_x_xyz_0, \
                             ta1_y_x_xzz_0, \
                             ta1_y_x_yyy_0, \
                             ta1_y_x_yyz_0, \
                             ta1_y_x_yzz_0, \
                             ta1_y_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_x_xxx_0[i] = 3.0 * ta1_y_0_xx_0[i] * fe_0 - 3.0 * ta1_y_0_xx_1[i] * fe_0 + ta1_y_0_xxx_0[i] * pa_x[i] -
                           ta1_y_0_xxx_1[i] * pc_x[i];

        ta1_y_x_xxy_0[i] = 2.0 * ta1_y_0_xy_0[i] * fe_0 - 2.0 * ta1_y_0_xy_1[i] * fe_0 + ta1_y_0_xxy_0[i] * pa_x[i] -
                           ta1_y_0_xxy_1[i] * pc_x[i];

        ta1_y_x_xxz_0[i] = 2.0 * ta1_y_0_xz_0[i] * fe_0 - 2.0 * ta1_y_0_xz_1[i] * fe_0 + ta1_y_0_xxz_0[i] * pa_x[i] -
                           ta1_y_0_xxz_1[i] * pc_x[i];

        ta1_y_x_xyy_0[i] =
            ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + ta1_y_0_xyy_0[i] * pa_x[i] - ta1_y_0_xyy_1[i] * pc_x[i];

        ta1_y_x_xyz_0[i] =
            ta1_y_0_yz_0[i] * fe_0 - ta1_y_0_yz_1[i] * fe_0 + ta1_y_0_xyz_0[i] * pa_x[i] - ta1_y_0_xyz_1[i] * pc_x[i];

        ta1_y_x_xzz_0[i] =
            ta1_y_0_zz_0[i] * fe_0 - ta1_y_0_zz_1[i] * fe_0 + ta1_y_0_xzz_0[i] * pa_x[i] - ta1_y_0_xzz_1[i] * pc_x[i];

        ta1_y_x_yyy_0[i] = ta1_y_0_yyy_0[i] * pa_x[i] - ta1_y_0_yyy_1[i] * pc_x[i];

        ta1_y_x_yyz_0[i] = ta1_y_0_yyz_0[i] * pa_x[i] - ta1_y_0_yyz_1[i] * pc_x[i];

        ta1_y_x_yzz_0[i] = ta1_y_0_yzz_0[i] * pa_x[i] - ta1_y_0_yzz_1[i] * pc_x[i];

        ta1_y_x_zzz_0[i] = ta1_y_0_zzz_0[i] * pa_x[i] - ta1_y_0_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : PF

    auto ta1_y_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 40);

    auto ta1_y_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 41);

    auto ta1_y_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 42);

    auto ta1_y_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 43);

    auto ta1_y_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 44);

    auto ta1_y_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 45);

    auto ta1_y_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 46);

    auto ta1_y_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 47);

    auto ta1_y_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 48);

    auto ta1_y_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 49);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_y_0_xx_0,  \
                             ta1_y_0_xx_1,  \
                             ta1_y_0_xxx_0, \
                             ta1_y_0_xxx_1, \
                             ta1_y_0_xxy_0, \
                             ta1_y_0_xxy_1, \
                             ta1_y_0_xxz_0, \
                             ta1_y_0_xxz_1, \
                             ta1_y_0_xy_0,  \
                             ta1_y_0_xy_1,  \
                             ta1_y_0_xyy_0, \
                             ta1_y_0_xyy_1, \
                             ta1_y_0_xyz_0, \
                             ta1_y_0_xyz_1, \
                             ta1_y_0_xz_0,  \
                             ta1_y_0_xz_1,  \
                             ta1_y_0_xzz_0, \
                             ta1_y_0_xzz_1, \
                             ta1_y_0_yy_0,  \
                             ta1_y_0_yy_1,  \
                             ta1_y_0_yyy_0, \
                             ta1_y_0_yyy_1, \
                             ta1_y_0_yyz_0, \
                             ta1_y_0_yyz_1, \
                             ta1_y_0_yz_0,  \
                             ta1_y_0_yz_1,  \
                             ta1_y_0_yzz_0, \
                             ta1_y_0_yzz_1, \
                             ta1_y_0_zz_0,  \
                             ta1_y_0_zz_1,  \
                             ta1_y_0_zzz_0, \
                             ta1_y_0_zzz_1, \
                             ta1_y_y_xxx_0, \
                             ta1_y_y_xxy_0, \
                             ta1_y_y_xxz_0, \
                             ta1_y_y_xyy_0, \
                             ta1_y_y_xyz_0, \
                             ta1_y_y_xzz_0, \
                             ta1_y_y_yyy_0, \
                             ta1_y_y_yyz_0, \
                             ta1_y_y_yzz_0, \
                             ta1_y_y_zzz_0, \
                             ta_0_xxx_1,    \
                             ta_0_xxy_1,    \
                             ta_0_xxz_1,    \
                             ta_0_xyy_1,    \
                             ta_0_xyz_1,    \
                             ta_0_xzz_1,    \
                             ta_0_yyy_1,    \
                             ta_0_yyz_1,    \
                             ta_0_yzz_1,    \
                             ta_0_zzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_y_xxx_0[i] = ta_0_xxx_1[i] + ta1_y_0_xxx_0[i] * pa_y[i] - ta1_y_0_xxx_1[i] * pc_y[i];

        ta1_y_y_xxy_0[i] = ta1_y_0_xx_0[i] * fe_0 - ta1_y_0_xx_1[i] * fe_0 + ta_0_xxy_1[i] +
                           ta1_y_0_xxy_0[i] * pa_y[i] - ta1_y_0_xxy_1[i] * pc_y[i];

        ta1_y_y_xxz_0[i] = ta_0_xxz_1[i] + ta1_y_0_xxz_0[i] * pa_y[i] - ta1_y_0_xxz_1[i] * pc_y[i];

        ta1_y_y_xyy_0[i] = 2.0 * ta1_y_0_xy_0[i] * fe_0 - 2.0 * ta1_y_0_xy_1[i] * fe_0 + ta_0_xyy_1[i] +
                           ta1_y_0_xyy_0[i] * pa_y[i] - ta1_y_0_xyy_1[i] * pc_y[i];

        ta1_y_y_xyz_0[i] = ta1_y_0_xz_0[i] * fe_0 - ta1_y_0_xz_1[i] * fe_0 + ta_0_xyz_1[i] +
                           ta1_y_0_xyz_0[i] * pa_y[i] - ta1_y_0_xyz_1[i] * pc_y[i];

        ta1_y_y_xzz_0[i] = ta_0_xzz_1[i] + ta1_y_0_xzz_0[i] * pa_y[i] - ta1_y_0_xzz_1[i] * pc_y[i];

        ta1_y_y_yyy_0[i] = 3.0 * ta1_y_0_yy_0[i] * fe_0 - 3.0 * ta1_y_0_yy_1[i] * fe_0 + ta_0_yyy_1[i] +
                           ta1_y_0_yyy_0[i] * pa_y[i] - ta1_y_0_yyy_1[i] * pc_y[i];

        ta1_y_y_yyz_0[i] = 2.0 * ta1_y_0_yz_0[i] * fe_0 - 2.0 * ta1_y_0_yz_1[i] * fe_0 + ta_0_yyz_1[i] +
                           ta1_y_0_yyz_0[i] * pa_y[i] - ta1_y_0_yyz_1[i] * pc_y[i];

        ta1_y_y_yzz_0[i] = ta1_y_0_zz_0[i] * fe_0 - ta1_y_0_zz_1[i] * fe_0 + ta_0_yzz_1[i] +
                           ta1_y_0_yzz_0[i] * pa_y[i] - ta1_y_0_yzz_1[i] * pc_y[i];

        ta1_y_y_zzz_0[i] = ta_0_zzz_1[i] + ta1_y_0_zzz_0[i] * pa_y[i] - ta1_y_0_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : PF

    auto ta1_y_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 50);

    auto ta1_y_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 51);

    auto ta1_y_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 52);

    auto ta1_y_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 53);

    auto ta1_y_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 54);

    auto ta1_y_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 55);

    auto ta1_y_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 56);

    auto ta1_y_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 57);

    auto ta1_y_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 58);

    auto ta1_y_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 59);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_y_0_xx_0,  \
                             ta1_y_0_xx_1,  \
                             ta1_y_0_xxx_0, \
                             ta1_y_0_xxx_1, \
                             ta1_y_0_xxy_0, \
                             ta1_y_0_xxy_1, \
                             ta1_y_0_xxz_0, \
                             ta1_y_0_xxz_1, \
                             ta1_y_0_xy_0,  \
                             ta1_y_0_xy_1,  \
                             ta1_y_0_xyy_0, \
                             ta1_y_0_xyy_1, \
                             ta1_y_0_xyz_0, \
                             ta1_y_0_xyz_1, \
                             ta1_y_0_xz_0,  \
                             ta1_y_0_xz_1,  \
                             ta1_y_0_xzz_0, \
                             ta1_y_0_xzz_1, \
                             ta1_y_0_yy_0,  \
                             ta1_y_0_yy_1,  \
                             ta1_y_0_yyy_0, \
                             ta1_y_0_yyy_1, \
                             ta1_y_0_yyz_0, \
                             ta1_y_0_yyz_1, \
                             ta1_y_0_yz_0,  \
                             ta1_y_0_yz_1,  \
                             ta1_y_0_yzz_0, \
                             ta1_y_0_yzz_1, \
                             ta1_y_0_zz_0,  \
                             ta1_y_0_zz_1,  \
                             ta1_y_0_zzz_0, \
                             ta1_y_0_zzz_1, \
                             ta1_y_z_xxx_0, \
                             ta1_y_z_xxy_0, \
                             ta1_y_z_xxz_0, \
                             ta1_y_z_xyy_0, \
                             ta1_y_z_xyz_0, \
                             ta1_y_z_xzz_0, \
                             ta1_y_z_yyy_0, \
                             ta1_y_z_yyz_0, \
                             ta1_y_z_yzz_0, \
                             ta1_y_z_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_z_xxx_0[i] = ta1_y_0_xxx_0[i] * pa_z[i] - ta1_y_0_xxx_1[i] * pc_z[i];

        ta1_y_z_xxy_0[i] = ta1_y_0_xxy_0[i] * pa_z[i] - ta1_y_0_xxy_1[i] * pc_z[i];

        ta1_y_z_xxz_0[i] =
            ta1_y_0_xx_0[i] * fe_0 - ta1_y_0_xx_1[i] * fe_0 + ta1_y_0_xxz_0[i] * pa_z[i] - ta1_y_0_xxz_1[i] * pc_z[i];

        ta1_y_z_xyy_0[i] = ta1_y_0_xyy_0[i] * pa_z[i] - ta1_y_0_xyy_1[i] * pc_z[i];

        ta1_y_z_xyz_0[i] =
            ta1_y_0_xy_0[i] * fe_0 - ta1_y_0_xy_1[i] * fe_0 + ta1_y_0_xyz_0[i] * pa_z[i] - ta1_y_0_xyz_1[i] * pc_z[i];

        ta1_y_z_xzz_0[i] = 2.0 * ta1_y_0_xz_0[i] * fe_0 - 2.0 * ta1_y_0_xz_1[i] * fe_0 + ta1_y_0_xzz_0[i] * pa_z[i] -
                           ta1_y_0_xzz_1[i] * pc_z[i];

        ta1_y_z_yyy_0[i] = ta1_y_0_yyy_0[i] * pa_z[i] - ta1_y_0_yyy_1[i] * pc_z[i];

        ta1_y_z_yyz_0[i] =
            ta1_y_0_yy_0[i] * fe_0 - ta1_y_0_yy_1[i] * fe_0 + ta1_y_0_yyz_0[i] * pa_z[i] - ta1_y_0_yyz_1[i] * pc_z[i];

        ta1_y_z_yzz_0[i] = 2.0 * ta1_y_0_yz_0[i] * fe_0 - 2.0 * ta1_y_0_yz_1[i] * fe_0 + ta1_y_0_yzz_0[i] * pa_z[i] -
                           ta1_y_0_yzz_1[i] * pc_z[i];

        ta1_y_z_zzz_0[i] = 3.0 * ta1_y_0_zz_0[i] * fe_0 - 3.0 * ta1_y_0_zz_1[i] * fe_0 + ta1_y_0_zzz_0[i] * pa_z[i] -
                           ta1_y_0_zzz_1[i] * pc_z[i];
    }

    // Set up 60-70 components of targeted buffer : PF

    auto ta1_z_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 60);

    auto ta1_z_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 61);

    auto ta1_z_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 62);

    auto ta1_z_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 63);

    auto ta1_z_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 64);

    auto ta1_z_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 65);

    auto ta1_z_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 66);

    auto ta1_z_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 67);

    auto ta1_z_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 68);

    auto ta1_z_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 69);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta1_z_0_xx_0,  \
                             ta1_z_0_xx_1,  \
                             ta1_z_0_xxx_0, \
                             ta1_z_0_xxx_1, \
                             ta1_z_0_xxy_0, \
                             ta1_z_0_xxy_1, \
                             ta1_z_0_xxz_0, \
                             ta1_z_0_xxz_1, \
                             ta1_z_0_xy_0,  \
                             ta1_z_0_xy_1,  \
                             ta1_z_0_xyy_0, \
                             ta1_z_0_xyy_1, \
                             ta1_z_0_xyz_0, \
                             ta1_z_0_xyz_1, \
                             ta1_z_0_xz_0,  \
                             ta1_z_0_xz_1,  \
                             ta1_z_0_xzz_0, \
                             ta1_z_0_xzz_1, \
                             ta1_z_0_yy_0,  \
                             ta1_z_0_yy_1,  \
                             ta1_z_0_yyy_0, \
                             ta1_z_0_yyy_1, \
                             ta1_z_0_yyz_0, \
                             ta1_z_0_yyz_1, \
                             ta1_z_0_yz_0,  \
                             ta1_z_0_yz_1,  \
                             ta1_z_0_yzz_0, \
                             ta1_z_0_yzz_1, \
                             ta1_z_0_zz_0,  \
                             ta1_z_0_zz_1,  \
                             ta1_z_0_zzz_0, \
                             ta1_z_0_zzz_1, \
                             ta1_z_x_xxx_0, \
                             ta1_z_x_xxy_0, \
                             ta1_z_x_xxz_0, \
                             ta1_z_x_xyy_0, \
                             ta1_z_x_xyz_0, \
                             ta1_z_x_xzz_0, \
                             ta1_z_x_yyy_0, \
                             ta1_z_x_yyz_0, \
                             ta1_z_x_yzz_0, \
                             ta1_z_x_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_x_xxx_0[i] = 3.0 * ta1_z_0_xx_0[i] * fe_0 - 3.0 * ta1_z_0_xx_1[i] * fe_0 + ta1_z_0_xxx_0[i] * pa_x[i] -
                           ta1_z_0_xxx_1[i] * pc_x[i];

        ta1_z_x_xxy_0[i] = 2.0 * ta1_z_0_xy_0[i] * fe_0 - 2.0 * ta1_z_0_xy_1[i] * fe_0 + ta1_z_0_xxy_0[i] * pa_x[i] -
                           ta1_z_0_xxy_1[i] * pc_x[i];

        ta1_z_x_xxz_0[i] = 2.0 * ta1_z_0_xz_0[i] * fe_0 - 2.0 * ta1_z_0_xz_1[i] * fe_0 + ta1_z_0_xxz_0[i] * pa_x[i] -
                           ta1_z_0_xxz_1[i] * pc_x[i];

        ta1_z_x_xyy_0[i] =
            ta1_z_0_yy_0[i] * fe_0 - ta1_z_0_yy_1[i] * fe_0 + ta1_z_0_xyy_0[i] * pa_x[i] - ta1_z_0_xyy_1[i] * pc_x[i];

        ta1_z_x_xyz_0[i] =
            ta1_z_0_yz_0[i] * fe_0 - ta1_z_0_yz_1[i] * fe_0 + ta1_z_0_xyz_0[i] * pa_x[i] - ta1_z_0_xyz_1[i] * pc_x[i];

        ta1_z_x_xzz_0[i] =
            ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + ta1_z_0_xzz_0[i] * pa_x[i] - ta1_z_0_xzz_1[i] * pc_x[i];

        ta1_z_x_yyy_0[i] = ta1_z_0_yyy_0[i] * pa_x[i] - ta1_z_0_yyy_1[i] * pc_x[i];

        ta1_z_x_yyz_0[i] = ta1_z_0_yyz_0[i] * pa_x[i] - ta1_z_0_yyz_1[i] * pc_x[i];

        ta1_z_x_yzz_0[i] = ta1_z_0_yzz_0[i] * pa_x[i] - ta1_z_0_yzz_1[i] * pc_x[i];

        ta1_z_x_zzz_0[i] = ta1_z_0_zzz_0[i] * pa_x[i] - ta1_z_0_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : PF

    auto ta1_z_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 70);

    auto ta1_z_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 71);

    auto ta1_z_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 72);

    auto ta1_z_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 73);

    auto ta1_z_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 74);

    auto ta1_z_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 75);

    auto ta1_z_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 76);

    auto ta1_z_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 77);

    auto ta1_z_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 78);

    auto ta1_z_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 79);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta1_z_0_xx_0,  \
                             ta1_z_0_xx_1,  \
                             ta1_z_0_xxx_0, \
                             ta1_z_0_xxx_1, \
                             ta1_z_0_xxy_0, \
                             ta1_z_0_xxy_1, \
                             ta1_z_0_xxz_0, \
                             ta1_z_0_xxz_1, \
                             ta1_z_0_xy_0,  \
                             ta1_z_0_xy_1,  \
                             ta1_z_0_xyy_0, \
                             ta1_z_0_xyy_1, \
                             ta1_z_0_xyz_0, \
                             ta1_z_0_xyz_1, \
                             ta1_z_0_xz_0,  \
                             ta1_z_0_xz_1,  \
                             ta1_z_0_xzz_0, \
                             ta1_z_0_xzz_1, \
                             ta1_z_0_yy_0,  \
                             ta1_z_0_yy_1,  \
                             ta1_z_0_yyy_0, \
                             ta1_z_0_yyy_1, \
                             ta1_z_0_yyz_0, \
                             ta1_z_0_yyz_1, \
                             ta1_z_0_yz_0,  \
                             ta1_z_0_yz_1,  \
                             ta1_z_0_yzz_0, \
                             ta1_z_0_yzz_1, \
                             ta1_z_0_zz_0,  \
                             ta1_z_0_zz_1,  \
                             ta1_z_0_zzz_0, \
                             ta1_z_0_zzz_1, \
                             ta1_z_y_xxx_0, \
                             ta1_z_y_xxy_0, \
                             ta1_z_y_xxz_0, \
                             ta1_z_y_xyy_0, \
                             ta1_z_y_xyz_0, \
                             ta1_z_y_xzz_0, \
                             ta1_z_y_yyy_0, \
                             ta1_z_y_yyz_0, \
                             ta1_z_y_yzz_0, \
                             ta1_z_y_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_y_xxx_0[i] = ta1_z_0_xxx_0[i] * pa_y[i] - ta1_z_0_xxx_1[i] * pc_y[i];

        ta1_z_y_xxy_0[i] =
            ta1_z_0_xx_0[i] * fe_0 - ta1_z_0_xx_1[i] * fe_0 + ta1_z_0_xxy_0[i] * pa_y[i] - ta1_z_0_xxy_1[i] * pc_y[i];

        ta1_z_y_xxz_0[i] = ta1_z_0_xxz_0[i] * pa_y[i] - ta1_z_0_xxz_1[i] * pc_y[i];

        ta1_z_y_xyy_0[i] = 2.0 * ta1_z_0_xy_0[i] * fe_0 - 2.0 * ta1_z_0_xy_1[i] * fe_0 + ta1_z_0_xyy_0[i] * pa_y[i] -
                           ta1_z_0_xyy_1[i] * pc_y[i];

        ta1_z_y_xyz_0[i] =
            ta1_z_0_xz_0[i] * fe_0 - ta1_z_0_xz_1[i] * fe_0 + ta1_z_0_xyz_0[i] * pa_y[i] - ta1_z_0_xyz_1[i] * pc_y[i];

        ta1_z_y_xzz_0[i] = ta1_z_0_xzz_0[i] * pa_y[i] - ta1_z_0_xzz_1[i] * pc_y[i];

        ta1_z_y_yyy_0[i] = 3.0 * ta1_z_0_yy_0[i] * fe_0 - 3.0 * ta1_z_0_yy_1[i] * fe_0 + ta1_z_0_yyy_0[i] * pa_y[i] -
                           ta1_z_0_yyy_1[i] * pc_y[i];

        ta1_z_y_yyz_0[i] = 2.0 * ta1_z_0_yz_0[i] * fe_0 - 2.0 * ta1_z_0_yz_1[i] * fe_0 + ta1_z_0_yyz_0[i] * pa_y[i] -
                           ta1_z_0_yyz_1[i] * pc_y[i];

        ta1_z_y_yzz_0[i] =
            ta1_z_0_zz_0[i] * fe_0 - ta1_z_0_zz_1[i] * fe_0 + ta1_z_0_yzz_0[i] * pa_y[i] - ta1_z_0_yzz_1[i] * pc_y[i];

        ta1_z_y_zzz_0[i] = ta1_z_0_zzz_0[i] * pa_y[i] - ta1_z_0_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : PF

    auto ta1_z_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 80);

    auto ta1_z_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 81);

    auto ta1_z_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 82);

    auto ta1_z_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 83);

    auto ta1_z_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 84);

    auto ta1_z_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 85);

    auto ta1_z_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 86);

    auto ta1_z_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 87);

    auto ta1_z_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 88);

    auto ta1_z_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 89);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta1_z_0_xx_0,  \
                             ta1_z_0_xx_1,  \
                             ta1_z_0_xxx_0, \
                             ta1_z_0_xxx_1, \
                             ta1_z_0_xxy_0, \
                             ta1_z_0_xxy_1, \
                             ta1_z_0_xxz_0, \
                             ta1_z_0_xxz_1, \
                             ta1_z_0_xy_0,  \
                             ta1_z_0_xy_1,  \
                             ta1_z_0_xyy_0, \
                             ta1_z_0_xyy_1, \
                             ta1_z_0_xyz_0, \
                             ta1_z_0_xyz_1, \
                             ta1_z_0_xz_0,  \
                             ta1_z_0_xz_1,  \
                             ta1_z_0_xzz_0, \
                             ta1_z_0_xzz_1, \
                             ta1_z_0_yy_0,  \
                             ta1_z_0_yy_1,  \
                             ta1_z_0_yyy_0, \
                             ta1_z_0_yyy_1, \
                             ta1_z_0_yyz_0, \
                             ta1_z_0_yyz_1, \
                             ta1_z_0_yz_0,  \
                             ta1_z_0_yz_1,  \
                             ta1_z_0_yzz_0, \
                             ta1_z_0_yzz_1, \
                             ta1_z_0_zz_0,  \
                             ta1_z_0_zz_1,  \
                             ta1_z_0_zzz_0, \
                             ta1_z_0_zzz_1, \
                             ta1_z_z_xxx_0, \
                             ta1_z_z_xxy_0, \
                             ta1_z_z_xxz_0, \
                             ta1_z_z_xyy_0, \
                             ta1_z_z_xyz_0, \
                             ta1_z_z_xzz_0, \
                             ta1_z_z_yyy_0, \
                             ta1_z_z_yyz_0, \
                             ta1_z_z_yzz_0, \
                             ta1_z_z_zzz_0, \
                             ta_0_xxx_1,    \
                             ta_0_xxy_1,    \
                             ta_0_xxz_1,    \
                             ta_0_xyy_1,    \
                             ta_0_xyz_1,    \
                             ta_0_xzz_1,    \
                             ta_0_yyy_1,    \
                             ta_0_yyz_1,    \
                             ta_0_yzz_1,    \
                             ta_0_zzz_1,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_z_xxx_0[i] = ta_0_xxx_1[i] + ta1_z_0_xxx_0[i] * pa_z[i] - ta1_z_0_xxx_1[i] * pc_z[i];

        ta1_z_z_xxy_0[i] = ta_0_xxy_1[i] + ta1_z_0_xxy_0[i] * pa_z[i] - ta1_z_0_xxy_1[i] * pc_z[i];

        ta1_z_z_xxz_0[i] = ta1_z_0_xx_0[i] * fe_0 - ta1_z_0_xx_1[i] * fe_0 + ta_0_xxz_1[i] +
                           ta1_z_0_xxz_0[i] * pa_z[i] - ta1_z_0_xxz_1[i] * pc_z[i];

        ta1_z_z_xyy_0[i] = ta_0_xyy_1[i] + ta1_z_0_xyy_0[i] * pa_z[i] - ta1_z_0_xyy_1[i] * pc_z[i];

        ta1_z_z_xyz_0[i] = ta1_z_0_xy_0[i] * fe_0 - ta1_z_0_xy_1[i] * fe_0 + ta_0_xyz_1[i] +
                           ta1_z_0_xyz_0[i] * pa_z[i] - ta1_z_0_xyz_1[i] * pc_z[i];

        ta1_z_z_xzz_0[i] = 2.0 * ta1_z_0_xz_0[i] * fe_0 - 2.0 * ta1_z_0_xz_1[i] * fe_0 + ta_0_xzz_1[i] +
                           ta1_z_0_xzz_0[i] * pa_z[i] - ta1_z_0_xzz_1[i] * pc_z[i];

        ta1_z_z_yyy_0[i] = ta_0_yyy_1[i] + ta1_z_0_yyy_0[i] * pa_z[i] - ta1_z_0_yyy_1[i] * pc_z[i];

        ta1_z_z_yyz_0[i] = ta1_z_0_yy_0[i] * fe_0 - ta1_z_0_yy_1[i] * fe_0 + ta_0_yyz_1[i] +
                           ta1_z_0_yyz_0[i] * pa_z[i] - ta1_z_0_yyz_1[i] * pc_z[i];

        ta1_z_z_yzz_0[i] = 2.0 * ta1_z_0_yz_0[i] * fe_0 - 2.0 * ta1_z_0_yz_1[i] * fe_0 + ta_0_yzz_1[i] +
                           ta1_z_0_yzz_0[i] * pa_z[i] - ta1_z_0_yzz_1[i] * pc_z[i];

        ta1_z_z_zzz_0[i] = 3.0 * ta1_z_0_zz_0[i] * fe_0 - 3.0 * ta1_z_0_zz_1[i] * fe_0 + ta_0_zzz_1[i] +
                           ta1_z_0_zzz_0[i] * pa_z[i] - ta1_z_0_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
