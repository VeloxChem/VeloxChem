#include "NuclearPotentialPrimRecGD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_gd(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_gd,
                               const size_t              idx_npot_0_dd,
                               const size_t              idx_npot_1_dd,
                               const size_t              idx_npot_0_fp,
                               const size_t              idx_npot_1_fp,
                               const size_t              idx_npot_0_fd,
                               const size_t              idx_npot_1_fd,
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

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_0 = pbuffer.data(idx_npot_0_dd);

    auto ta_xx_xy_0 = pbuffer.data(idx_npot_0_dd + 1);

    auto ta_xx_xz_0 = pbuffer.data(idx_npot_0_dd + 2);

    auto ta_xx_yy_0 = pbuffer.data(idx_npot_0_dd + 3);

    auto ta_xx_yz_0 = pbuffer.data(idx_npot_0_dd + 4);

    auto ta_xx_zz_0 = pbuffer.data(idx_npot_0_dd + 5);

    auto ta_xy_yy_0 = pbuffer.data(idx_npot_0_dd + 9);

    auto ta_xy_yz_0 = pbuffer.data(idx_npot_0_dd + 10);

    auto ta_xz_yz_0 = pbuffer.data(idx_npot_0_dd + 16);

    auto ta_xz_zz_0 = pbuffer.data(idx_npot_0_dd + 17);

    auto ta_yy_xx_0 = pbuffer.data(idx_npot_0_dd + 18);

    auto ta_yy_xy_0 = pbuffer.data(idx_npot_0_dd + 19);

    auto ta_yy_xz_0 = pbuffer.data(idx_npot_0_dd + 20);

    auto ta_yy_yy_0 = pbuffer.data(idx_npot_0_dd + 21);

    auto ta_yy_yz_0 = pbuffer.data(idx_npot_0_dd + 22);

    auto ta_yy_zz_0 = pbuffer.data(idx_npot_0_dd + 23);

    auto ta_yz_xz_0 = pbuffer.data(idx_npot_0_dd + 26);

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

    auto ta_xy_yy_1 = pbuffer.data(idx_npot_1_dd + 9);

    auto ta_xy_yz_1 = pbuffer.data(idx_npot_1_dd + 10);

    auto ta_xz_yz_1 = pbuffer.data(idx_npot_1_dd + 16);

    auto ta_xz_zz_1 = pbuffer.data(idx_npot_1_dd + 17);

    auto ta_yy_xx_1 = pbuffer.data(idx_npot_1_dd + 18);

    auto ta_yy_xy_1 = pbuffer.data(idx_npot_1_dd + 19);

    auto ta_yy_xz_1 = pbuffer.data(idx_npot_1_dd + 20);

    auto ta_yy_yy_1 = pbuffer.data(idx_npot_1_dd + 21);

    auto ta_yy_yz_1 = pbuffer.data(idx_npot_1_dd + 22);

    auto ta_yy_zz_1 = pbuffer.data(idx_npot_1_dd + 23);

    auto ta_yz_xz_1 = pbuffer.data(idx_npot_1_dd + 26);

    auto ta_yz_yz_1 = pbuffer.data(idx_npot_1_dd + 28);

    auto ta_yz_zz_1 = pbuffer.data(idx_npot_1_dd + 29);

    auto ta_zz_xx_1 = pbuffer.data(idx_npot_1_dd + 30);

    auto ta_zz_xy_1 = pbuffer.data(idx_npot_1_dd + 31);

    auto ta_zz_xz_1 = pbuffer.data(idx_npot_1_dd + 32);

    auto ta_zz_yy_1 = pbuffer.data(idx_npot_1_dd + 33);

    auto ta_zz_yz_1 = pbuffer.data(idx_npot_1_dd + 34);

    auto ta_zz_zz_1 = pbuffer.data(idx_npot_1_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto ta_xxx_x_0 = pbuffer.data(idx_npot_0_fp);

    auto ta_xxx_y_0 = pbuffer.data(idx_npot_0_fp + 1);

    auto ta_xxx_z_0 = pbuffer.data(idx_npot_0_fp + 2);

    auto ta_xyy_y_0 = pbuffer.data(idx_npot_0_fp + 10);

    auto ta_xzz_z_0 = pbuffer.data(idx_npot_0_fp + 17);

    auto ta_yyy_x_0 = pbuffer.data(idx_npot_0_fp + 18);

    auto ta_yyy_y_0 = pbuffer.data(idx_npot_0_fp + 19);

    auto ta_yyy_z_0 = pbuffer.data(idx_npot_0_fp + 20);

    auto ta_yyz_z_0 = pbuffer.data(idx_npot_0_fp + 23);

    auto ta_yzz_y_0 = pbuffer.data(idx_npot_0_fp + 25);

    auto ta_yzz_z_0 = pbuffer.data(idx_npot_0_fp + 26);

    auto ta_zzz_x_0 = pbuffer.data(idx_npot_0_fp + 27);

    auto ta_zzz_y_0 = pbuffer.data(idx_npot_0_fp + 28);

    auto ta_zzz_z_0 = pbuffer.data(idx_npot_0_fp + 29);

    // Set up components of auxiliary buffer : FP

    auto ta_xxx_x_1 = pbuffer.data(idx_npot_1_fp);

    auto ta_xxx_y_1 = pbuffer.data(idx_npot_1_fp + 1);

    auto ta_xxx_z_1 = pbuffer.data(idx_npot_1_fp + 2);

    auto ta_xyy_y_1 = pbuffer.data(idx_npot_1_fp + 10);

    auto ta_xzz_z_1 = pbuffer.data(idx_npot_1_fp + 17);

    auto ta_yyy_x_1 = pbuffer.data(idx_npot_1_fp + 18);

    auto ta_yyy_y_1 = pbuffer.data(idx_npot_1_fp + 19);

    auto ta_yyy_z_1 = pbuffer.data(idx_npot_1_fp + 20);

    auto ta_yyz_z_1 = pbuffer.data(idx_npot_1_fp + 23);

    auto ta_yzz_y_1 = pbuffer.data(idx_npot_1_fp + 25);

    auto ta_yzz_z_1 = pbuffer.data(idx_npot_1_fp + 26);

    auto ta_zzz_x_1 = pbuffer.data(idx_npot_1_fp + 27);

    auto ta_zzz_y_1 = pbuffer.data(idx_npot_1_fp + 28);

    auto ta_zzz_z_1 = pbuffer.data(idx_npot_1_fp + 29);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_0 = pbuffer.data(idx_npot_0_fd);

    auto ta_xxx_xy_0 = pbuffer.data(idx_npot_0_fd + 1);

    auto ta_xxx_xz_0 = pbuffer.data(idx_npot_0_fd + 2);

    auto ta_xxx_yy_0 = pbuffer.data(idx_npot_0_fd + 3);

    auto ta_xxx_yz_0 = pbuffer.data(idx_npot_0_fd + 4);

    auto ta_xxx_zz_0 = pbuffer.data(idx_npot_0_fd + 5);

    auto ta_xxy_xx_0 = pbuffer.data(idx_npot_0_fd + 6);

    auto ta_xxy_xy_0 = pbuffer.data(idx_npot_0_fd + 7);

    auto ta_xxy_xz_0 = pbuffer.data(idx_npot_0_fd + 8);

    auto ta_xxy_yy_0 = pbuffer.data(idx_npot_0_fd + 9);

    auto ta_xxy_yz_0 = pbuffer.data(idx_npot_0_fd + 10);

    auto ta_xxz_xx_0 = pbuffer.data(idx_npot_0_fd + 12);

    auto ta_xxz_xy_0 = pbuffer.data(idx_npot_0_fd + 13);

    auto ta_xxz_xz_0 = pbuffer.data(idx_npot_0_fd + 14);

    auto ta_xxz_yz_0 = pbuffer.data(idx_npot_0_fd + 16);

    auto ta_xxz_zz_0 = pbuffer.data(idx_npot_0_fd + 17);

    auto ta_xyy_xx_0 = pbuffer.data(idx_npot_0_fd + 18);

    auto ta_xyy_xy_0 = pbuffer.data(idx_npot_0_fd + 19);

    auto ta_xyy_yy_0 = pbuffer.data(idx_npot_0_fd + 21);

    auto ta_xyy_yz_0 = pbuffer.data(idx_npot_0_fd + 22);

    auto ta_xyy_zz_0 = pbuffer.data(idx_npot_0_fd + 23);

    auto ta_xyz_yz_0 = pbuffer.data(idx_npot_0_fd + 28);

    auto ta_xzz_xx_0 = pbuffer.data(idx_npot_0_fd + 30);

    auto ta_xzz_xz_0 = pbuffer.data(idx_npot_0_fd + 32);

    auto ta_xzz_yy_0 = pbuffer.data(idx_npot_0_fd + 33);

    auto ta_xzz_yz_0 = pbuffer.data(idx_npot_0_fd + 34);

    auto ta_xzz_zz_0 = pbuffer.data(idx_npot_0_fd + 35);

    auto ta_yyy_xx_0 = pbuffer.data(idx_npot_0_fd + 36);

    auto ta_yyy_xy_0 = pbuffer.data(idx_npot_0_fd + 37);

    auto ta_yyy_xz_0 = pbuffer.data(idx_npot_0_fd + 38);

    auto ta_yyy_yy_0 = pbuffer.data(idx_npot_0_fd + 39);

    auto ta_yyy_yz_0 = pbuffer.data(idx_npot_0_fd + 40);

    auto ta_yyy_zz_0 = pbuffer.data(idx_npot_0_fd + 41);

    auto ta_yyz_xy_0 = pbuffer.data(idx_npot_0_fd + 43);

    auto ta_yyz_xz_0 = pbuffer.data(idx_npot_0_fd + 44);

    auto ta_yyz_yy_0 = pbuffer.data(idx_npot_0_fd + 45);

    auto ta_yyz_yz_0 = pbuffer.data(idx_npot_0_fd + 46);

    auto ta_yyz_zz_0 = pbuffer.data(idx_npot_0_fd + 47);

    auto ta_yzz_xx_0 = pbuffer.data(idx_npot_0_fd + 48);

    auto ta_yzz_xy_0 = pbuffer.data(idx_npot_0_fd + 49);

    auto ta_yzz_xz_0 = pbuffer.data(idx_npot_0_fd + 50);

    auto ta_yzz_yy_0 = pbuffer.data(idx_npot_0_fd + 51);

    auto ta_yzz_yz_0 = pbuffer.data(idx_npot_0_fd + 52);

    auto ta_yzz_zz_0 = pbuffer.data(idx_npot_0_fd + 53);

    auto ta_zzz_xx_0 = pbuffer.data(idx_npot_0_fd + 54);

    auto ta_zzz_xy_0 = pbuffer.data(idx_npot_0_fd + 55);

    auto ta_zzz_xz_0 = pbuffer.data(idx_npot_0_fd + 56);

    auto ta_zzz_yy_0 = pbuffer.data(idx_npot_0_fd + 57);

    auto ta_zzz_yz_0 = pbuffer.data(idx_npot_0_fd + 58);

    auto ta_zzz_zz_0 = pbuffer.data(idx_npot_0_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_1 = pbuffer.data(idx_npot_1_fd);

    auto ta_xxx_xy_1 = pbuffer.data(idx_npot_1_fd + 1);

    auto ta_xxx_xz_1 = pbuffer.data(idx_npot_1_fd + 2);

    auto ta_xxx_yy_1 = pbuffer.data(idx_npot_1_fd + 3);

    auto ta_xxx_yz_1 = pbuffer.data(idx_npot_1_fd + 4);

    auto ta_xxx_zz_1 = pbuffer.data(idx_npot_1_fd + 5);

    auto ta_xxy_xx_1 = pbuffer.data(idx_npot_1_fd + 6);

    auto ta_xxy_xy_1 = pbuffer.data(idx_npot_1_fd + 7);

    auto ta_xxy_xz_1 = pbuffer.data(idx_npot_1_fd + 8);

    auto ta_xxy_yy_1 = pbuffer.data(idx_npot_1_fd + 9);

    auto ta_xxy_yz_1 = pbuffer.data(idx_npot_1_fd + 10);

    auto ta_xxz_xx_1 = pbuffer.data(idx_npot_1_fd + 12);

    auto ta_xxz_xy_1 = pbuffer.data(idx_npot_1_fd + 13);

    auto ta_xxz_xz_1 = pbuffer.data(idx_npot_1_fd + 14);

    auto ta_xxz_yz_1 = pbuffer.data(idx_npot_1_fd + 16);

    auto ta_xxz_zz_1 = pbuffer.data(idx_npot_1_fd + 17);

    auto ta_xyy_xx_1 = pbuffer.data(idx_npot_1_fd + 18);

    auto ta_xyy_xy_1 = pbuffer.data(idx_npot_1_fd + 19);

    auto ta_xyy_yy_1 = pbuffer.data(idx_npot_1_fd + 21);

    auto ta_xyy_yz_1 = pbuffer.data(idx_npot_1_fd + 22);

    auto ta_xyy_zz_1 = pbuffer.data(idx_npot_1_fd + 23);

    auto ta_xyz_yz_1 = pbuffer.data(idx_npot_1_fd + 28);

    auto ta_xzz_xx_1 = pbuffer.data(idx_npot_1_fd + 30);

    auto ta_xzz_xz_1 = pbuffer.data(idx_npot_1_fd + 32);

    auto ta_xzz_yy_1 = pbuffer.data(idx_npot_1_fd + 33);

    auto ta_xzz_yz_1 = pbuffer.data(idx_npot_1_fd + 34);

    auto ta_xzz_zz_1 = pbuffer.data(idx_npot_1_fd + 35);

    auto ta_yyy_xx_1 = pbuffer.data(idx_npot_1_fd + 36);

    auto ta_yyy_xy_1 = pbuffer.data(idx_npot_1_fd + 37);

    auto ta_yyy_xz_1 = pbuffer.data(idx_npot_1_fd + 38);

    auto ta_yyy_yy_1 = pbuffer.data(idx_npot_1_fd + 39);

    auto ta_yyy_yz_1 = pbuffer.data(idx_npot_1_fd + 40);

    auto ta_yyy_zz_1 = pbuffer.data(idx_npot_1_fd + 41);

    auto ta_yyz_xy_1 = pbuffer.data(idx_npot_1_fd + 43);

    auto ta_yyz_xz_1 = pbuffer.data(idx_npot_1_fd + 44);

    auto ta_yyz_yy_1 = pbuffer.data(idx_npot_1_fd + 45);

    auto ta_yyz_yz_1 = pbuffer.data(idx_npot_1_fd + 46);

    auto ta_yyz_zz_1 = pbuffer.data(idx_npot_1_fd + 47);

    auto ta_yzz_xx_1 = pbuffer.data(idx_npot_1_fd + 48);

    auto ta_yzz_xy_1 = pbuffer.data(idx_npot_1_fd + 49);

    auto ta_yzz_xz_1 = pbuffer.data(idx_npot_1_fd + 50);

    auto ta_yzz_yy_1 = pbuffer.data(idx_npot_1_fd + 51);

    auto ta_yzz_yz_1 = pbuffer.data(idx_npot_1_fd + 52);

    auto ta_yzz_zz_1 = pbuffer.data(idx_npot_1_fd + 53);

    auto ta_zzz_xx_1 = pbuffer.data(idx_npot_1_fd + 54);

    auto ta_zzz_xy_1 = pbuffer.data(idx_npot_1_fd + 55);

    auto ta_zzz_xz_1 = pbuffer.data(idx_npot_1_fd + 56);

    auto ta_zzz_yy_1 = pbuffer.data(idx_npot_1_fd + 57);

    auto ta_zzz_yz_1 = pbuffer.data(idx_npot_1_fd + 58);

    auto ta_zzz_zz_1 = pbuffer.data(idx_npot_1_fd + 59);

    // Set up 0-6 components of targeted buffer : GD

    auto ta_xxxx_xx_0 = pbuffer.data(idx_npot_0_gd);

    auto ta_xxxx_xy_0 = pbuffer.data(idx_npot_0_gd + 1);

    auto ta_xxxx_xz_0 = pbuffer.data(idx_npot_0_gd + 2);

    auto ta_xxxx_yy_0 = pbuffer.data(idx_npot_0_gd + 3);

    auto ta_xxxx_yz_0 = pbuffer.data(idx_npot_0_gd + 4);

    auto ta_xxxx_zz_0 = pbuffer.data(idx_npot_0_gd + 5);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xx_xx_0,   \
                             ta_xx_xx_1,   \
                             ta_xx_xy_0,   \
                             ta_xx_xy_1,   \
                             ta_xx_xz_0,   \
                             ta_xx_xz_1,   \
                             ta_xx_yy_0,   \
                             ta_xx_yy_1,   \
                             ta_xx_yz_0,   \
                             ta_xx_yz_1,   \
                             ta_xx_zz_0,   \
                             ta_xx_zz_1,   \
                             ta_xxx_x_0,   \
                             ta_xxx_x_1,   \
                             ta_xxx_xx_0,  \
                             ta_xxx_xx_1,  \
                             ta_xxx_xy_0,  \
                             ta_xxx_xy_1,  \
                             ta_xxx_xz_0,  \
                             ta_xxx_xz_1,  \
                             ta_xxx_y_0,   \
                             ta_xxx_y_1,   \
                             ta_xxx_yy_0,  \
                             ta_xxx_yy_1,  \
                             ta_xxx_yz_0,  \
                             ta_xxx_yz_1,  \
                             ta_xxx_z_0,   \
                             ta_xxx_z_1,   \
                             ta_xxx_zz_0,  \
                             ta_xxx_zz_1,  \
                             ta_xxxx_xx_0, \
                             ta_xxxx_xy_0, \
                             ta_xxxx_xz_0, \
                             ta_xxxx_yy_0, \
                             ta_xxxx_yz_0, \
                             ta_xxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_xx_0[i] = 3.0 * ta_xx_xx_0[i] * fe_0 - 3.0 * ta_xx_xx_1[i] * fe_0 + 2.0 * ta_xxx_x_0[i] * fe_0 - 2.0 * ta_xxx_x_1[i] * fe_0 +
                          ta_xxx_xx_0[i] * pa_x[i] - ta_xxx_xx_1[i] * pc_x[i];

        ta_xxxx_xy_0[i] = 3.0 * ta_xx_xy_0[i] * fe_0 - 3.0 * ta_xx_xy_1[i] * fe_0 + ta_xxx_y_0[i] * fe_0 - ta_xxx_y_1[i] * fe_0 +
                          ta_xxx_xy_0[i] * pa_x[i] - ta_xxx_xy_1[i] * pc_x[i];

        ta_xxxx_xz_0[i] = 3.0 * ta_xx_xz_0[i] * fe_0 - 3.0 * ta_xx_xz_1[i] * fe_0 + ta_xxx_z_0[i] * fe_0 - ta_xxx_z_1[i] * fe_0 +
                          ta_xxx_xz_0[i] * pa_x[i] - ta_xxx_xz_1[i] * pc_x[i];

        ta_xxxx_yy_0[i] = 3.0 * ta_xx_yy_0[i] * fe_0 - 3.0 * ta_xx_yy_1[i] * fe_0 + ta_xxx_yy_0[i] * pa_x[i] - ta_xxx_yy_1[i] * pc_x[i];

        ta_xxxx_yz_0[i] = 3.0 * ta_xx_yz_0[i] * fe_0 - 3.0 * ta_xx_yz_1[i] * fe_0 + ta_xxx_yz_0[i] * pa_x[i] - ta_xxx_yz_1[i] * pc_x[i];

        ta_xxxx_zz_0[i] = 3.0 * ta_xx_zz_0[i] * fe_0 - 3.0 * ta_xx_zz_1[i] * fe_0 + ta_xxx_zz_0[i] * pa_x[i] - ta_xxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto ta_xxxy_xx_0 = pbuffer.data(idx_npot_0_gd + 6);

    auto ta_xxxy_xy_0 = pbuffer.data(idx_npot_0_gd + 7);

    auto ta_xxxy_xz_0 = pbuffer.data(idx_npot_0_gd + 8);

    auto ta_xxxy_yy_0 = pbuffer.data(idx_npot_0_gd + 9);

    auto ta_xxxy_yz_0 = pbuffer.data(idx_npot_0_gd + 10);

    auto ta_xxxy_zz_0 = pbuffer.data(idx_npot_0_gd + 11);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xxx_x_0,   \
                             ta_xxx_x_1,   \
                             ta_xxx_xx_0,  \
                             ta_xxx_xx_1,  \
                             ta_xxx_xy_0,  \
                             ta_xxx_xy_1,  \
                             ta_xxx_xz_0,  \
                             ta_xxx_xz_1,  \
                             ta_xxx_zz_0,  \
                             ta_xxx_zz_1,  \
                             ta_xxxy_xx_0, \
                             ta_xxxy_xy_0, \
                             ta_xxxy_xz_0, \
                             ta_xxxy_yy_0, \
                             ta_xxxy_yz_0, \
                             ta_xxxy_zz_0, \
                             ta_xxy_yy_0,  \
                             ta_xxy_yy_1,  \
                             ta_xxy_yz_0,  \
                             ta_xxy_yz_1,  \
                             ta_xy_yy_0,   \
                             ta_xy_yy_1,   \
                             ta_xy_yz_0,   \
                             ta_xy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxy_xx_0[i] = ta_xxx_xx_0[i] * pa_y[i] - ta_xxx_xx_1[i] * pc_y[i];

        ta_xxxy_xy_0[i] = ta_xxx_x_0[i] * fe_0 - ta_xxx_x_1[i] * fe_0 + ta_xxx_xy_0[i] * pa_y[i] - ta_xxx_xy_1[i] * pc_y[i];

        ta_xxxy_xz_0[i] = ta_xxx_xz_0[i] * pa_y[i] - ta_xxx_xz_1[i] * pc_y[i];

        ta_xxxy_yy_0[i] = 2.0 * ta_xy_yy_0[i] * fe_0 - 2.0 * ta_xy_yy_1[i] * fe_0 + ta_xxy_yy_0[i] * pa_x[i] - ta_xxy_yy_1[i] * pc_x[i];

        ta_xxxy_yz_0[i] = 2.0 * ta_xy_yz_0[i] * fe_0 - 2.0 * ta_xy_yz_1[i] * fe_0 + ta_xxy_yz_0[i] * pa_x[i] - ta_xxy_yz_1[i] * pc_x[i];

        ta_xxxy_zz_0[i] = ta_xxx_zz_0[i] * pa_y[i] - ta_xxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto ta_xxxz_xx_0 = pbuffer.data(idx_npot_0_gd + 12);

    auto ta_xxxz_xy_0 = pbuffer.data(idx_npot_0_gd + 13);

    auto ta_xxxz_xz_0 = pbuffer.data(idx_npot_0_gd + 14);

    auto ta_xxxz_yy_0 = pbuffer.data(idx_npot_0_gd + 15);

    auto ta_xxxz_yz_0 = pbuffer.data(idx_npot_0_gd + 16);

    auto ta_xxxz_zz_0 = pbuffer.data(idx_npot_0_gd + 17);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xxx_x_0,   \
                             ta_xxx_x_1,   \
                             ta_xxx_xx_0,  \
                             ta_xxx_xx_1,  \
                             ta_xxx_xy_0,  \
                             ta_xxx_xy_1,  \
                             ta_xxx_xz_0,  \
                             ta_xxx_xz_1,  \
                             ta_xxx_yy_0,  \
                             ta_xxx_yy_1,  \
                             ta_xxxz_xx_0, \
                             ta_xxxz_xy_0, \
                             ta_xxxz_xz_0, \
                             ta_xxxz_yy_0, \
                             ta_xxxz_yz_0, \
                             ta_xxxz_zz_0, \
                             ta_xxz_yz_0,  \
                             ta_xxz_yz_1,  \
                             ta_xxz_zz_0,  \
                             ta_xxz_zz_1,  \
                             ta_xz_yz_0,   \
                             ta_xz_yz_1,   \
                             ta_xz_zz_0,   \
                             ta_xz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxz_xx_0[i] = ta_xxx_xx_0[i] * pa_z[i] - ta_xxx_xx_1[i] * pc_z[i];

        ta_xxxz_xy_0[i] = ta_xxx_xy_0[i] * pa_z[i] - ta_xxx_xy_1[i] * pc_z[i];

        ta_xxxz_xz_0[i] = ta_xxx_x_0[i] * fe_0 - ta_xxx_x_1[i] * fe_0 + ta_xxx_xz_0[i] * pa_z[i] - ta_xxx_xz_1[i] * pc_z[i];

        ta_xxxz_yy_0[i] = ta_xxx_yy_0[i] * pa_z[i] - ta_xxx_yy_1[i] * pc_z[i];

        ta_xxxz_yz_0[i] = 2.0 * ta_xz_yz_0[i] * fe_0 - 2.0 * ta_xz_yz_1[i] * fe_0 + ta_xxz_yz_0[i] * pa_x[i] - ta_xxz_yz_1[i] * pc_x[i];

        ta_xxxz_zz_0[i] = 2.0 * ta_xz_zz_0[i] * fe_0 - 2.0 * ta_xz_zz_1[i] * fe_0 + ta_xxz_zz_0[i] * pa_x[i] - ta_xxz_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto ta_xxyy_xx_0 = pbuffer.data(idx_npot_0_gd + 18);

    auto ta_xxyy_xy_0 = pbuffer.data(idx_npot_0_gd + 19);

    auto ta_xxyy_xz_0 = pbuffer.data(idx_npot_0_gd + 20);

    auto ta_xxyy_yy_0 = pbuffer.data(idx_npot_0_gd + 21);

    auto ta_xxyy_yz_0 = pbuffer.data(idx_npot_0_gd + 22);

    auto ta_xxyy_zz_0 = pbuffer.data(idx_npot_0_gd + 23);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xx_xx_0,   \
                             ta_xx_xx_1,   \
                             ta_xx_xz_0,   \
                             ta_xx_xz_1,   \
                             ta_xxy_xx_0,  \
                             ta_xxy_xx_1,  \
                             ta_xxy_xz_0,  \
                             ta_xxy_xz_1,  \
                             ta_xxyy_xx_0, \
                             ta_xxyy_xy_0, \
                             ta_xxyy_xz_0, \
                             ta_xxyy_yy_0, \
                             ta_xxyy_yz_0, \
                             ta_xxyy_zz_0, \
                             ta_xyy_xy_0,  \
                             ta_xyy_xy_1,  \
                             ta_xyy_y_0,   \
                             ta_xyy_y_1,   \
                             ta_xyy_yy_0,  \
                             ta_xyy_yy_1,  \
                             ta_xyy_yz_0,  \
                             ta_xyy_yz_1,  \
                             ta_xyy_zz_0,  \
                             ta_xyy_zz_1,  \
                             ta_yy_xy_0,   \
                             ta_yy_xy_1,   \
                             ta_yy_yy_0,   \
                             ta_yy_yy_1,   \
                             ta_yy_yz_0,   \
                             ta_yy_yz_1,   \
                             ta_yy_zz_0,   \
                             ta_yy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyy_xx_0[i] = ta_xx_xx_0[i] * fe_0 - ta_xx_xx_1[i] * fe_0 + ta_xxy_xx_0[i] * pa_y[i] - ta_xxy_xx_1[i] * pc_y[i];

        ta_xxyy_xy_0[i] = ta_yy_xy_0[i] * fe_0 - ta_yy_xy_1[i] * fe_0 + ta_xyy_y_0[i] * fe_0 - ta_xyy_y_1[i] * fe_0 + ta_xyy_xy_0[i] * pa_x[i] -
                          ta_xyy_xy_1[i] * pc_x[i];

        ta_xxyy_xz_0[i] = ta_xx_xz_0[i] * fe_0 - ta_xx_xz_1[i] * fe_0 + ta_xxy_xz_0[i] * pa_y[i] - ta_xxy_xz_1[i] * pc_y[i];

        ta_xxyy_yy_0[i] = ta_yy_yy_0[i] * fe_0 - ta_yy_yy_1[i] * fe_0 + ta_xyy_yy_0[i] * pa_x[i] - ta_xyy_yy_1[i] * pc_x[i];

        ta_xxyy_yz_0[i] = ta_yy_yz_0[i] * fe_0 - ta_yy_yz_1[i] * fe_0 + ta_xyy_yz_0[i] * pa_x[i] - ta_xyy_yz_1[i] * pc_x[i];

        ta_xxyy_zz_0[i] = ta_yy_zz_0[i] * fe_0 - ta_yy_zz_1[i] * fe_0 + ta_xyy_zz_0[i] * pa_x[i] - ta_xyy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto ta_xxyz_xx_0 = pbuffer.data(idx_npot_0_gd + 24);

    auto ta_xxyz_xy_0 = pbuffer.data(idx_npot_0_gd + 25);

    auto ta_xxyz_xz_0 = pbuffer.data(idx_npot_0_gd + 26);

    auto ta_xxyz_yy_0 = pbuffer.data(idx_npot_0_gd + 27);

    auto ta_xxyz_yz_0 = pbuffer.data(idx_npot_0_gd + 28);

    auto ta_xxyz_zz_0 = pbuffer.data(idx_npot_0_gd + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             pc_x,         \
                             pc_y,         \
                             pc_z,         \
                             ta_xxy_xy_0,  \
                             ta_xxy_xy_1,  \
                             ta_xxy_yy_0,  \
                             ta_xxy_yy_1,  \
                             ta_xxyz_xx_0, \
                             ta_xxyz_xy_0, \
                             ta_xxyz_xz_0, \
                             ta_xxyz_yy_0, \
                             ta_xxyz_yz_0, \
                             ta_xxyz_zz_0, \
                             ta_xxz_xx_0,  \
                             ta_xxz_xx_1,  \
                             ta_xxz_xz_0,  \
                             ta_xxz_xz_1,  \
                             ta_xxz_zz_0,  \
                             ta_xxz_zz_1,  \
                             ta_xyz_yz_0,  \
                             ta_xyz_yz_1,  \
                             ta_yz_yz_0,   \
                             ta_yz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyz_xx_0[i] = ta_xxz_xx_0[i] * pa_y[i] - ta_xxz_xx_1[i] * pc_y[i];

        ta_xxyz_xy_0[i] = ta_xxy_xy_0[i] * pa_z[i] - ta_xxy_xy_1[i] * pc_z[i];

        ta_xxyz_xz_0[i] = ta_xxz_xz_0[i] * pa_y[i] - ta_xxz_xz_1[i] * pc_y[i];

        ta_xxyz_yy_0[i] = ta_xxy_yy_0[i] * pa_z[i] - ta_xxy_yy_1[i] * pc_z[i];

        ta_xxyz_yz_0[i] = ta_yz_yz_0[i] * fe_0 - ta_yz_yz_1[i] * fe_0 + ta_xyz_yz_0[i] * pa_x[i] - ta_xyz_yz_1[i] * pc_x[i];

        ta_xxyz_zz_0[i] = ta_xxz_zz_0[i] * pa_y[i] - ta_xxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto ta_xxzz_xx_0 = pbuffer.data(idx_npot_0_gd + 30);

    auto ta_xxzz_xy_0 = pbuffer.data(idx_npot_0_gd + 31);

    auto ta_xxzz_xz_0 = pbuffer.data(idx_npot_0_gd + 32);

    auto ta_xxzz_yy_0 = pbuffer.data(idx_npot_0_gd + 33);

    auto ta_xxzz_yz_0 = pbuffer.data(idx_npot_0_gd + 34);

    auto ta_xxzz_zz_0 = pbuffer.data(idx_npot_0_gd + 35);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xx_xx_0,   \
                             ta_xx_xx_1,   \
                             ta_xx_xy_0,   \
                             ta_xx_xy_1,   \
                             ta_xxz_xx_0,  \
                             ta_xxz_xx_1,  \
                             ta_xxz_xy_0,  \
                             ta_xxz_xy_1,  \
                             ta_xxzz_xx_0, \
                             ta_xxzz_xy_0, \
                             ta_xxzz_xz_0, \
                             ta_xxzz_yy_0, \
                             ta_xxzz_yz_0, \
                             ta_xxzz_zz_0, \
                             ta_xzz_xz_0,  \
                             ta_xzz_xz_1,  \
                             ta_xzz_yy_0,  \
                             ta_xzz_yy_1,  \
                             ta_xzz_yz_0,  \
                             ta_xzz_yz_1,  \
                             ta_xzz_z_0,   \
                             ta_xzz_z_1,   \
                             ta_xzz_zz_0,  \
                             ta_xzz_zz_1,  \
                             ta_zz_xz_0,   \
                             ta_zz_xz_1,   \
                             ta_zz_yy_0,   \
                             ta_zz_yy_1,   \
                             ta_zz_yz_0,   \
                             ta_zz_yz_1,   \
                             ta_zz_zz_0,   \
                             ta_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzz_xx_0[i] = ta_xx_xx_0[i] * fe_0 - ta_xx_xx_1[i] * fe_0 + ta_xxz_xx_0[i] * pa_z[i] - ta_xxz_xx_1[i] * pc_z[i];

        ta_xxzz_xy_0[i] = ta_xx_xy_0[i] * fe_0 - ta_xx_xy_1[i] * fe_0 + ta_xxz_xy_0[i] * pa_z[i] - ta_xxz_xy_1[i] * pc_z[i];

        ta_xxzz_xz_0[i] = ta_zz_xz_0[i] * fe_0 - ta_zz_xz_1[i] * fe_0 + ta_xzz_z_0[i] * fe_0 - ta_xzz_z_1[i] * fe_0 + ta_xzz_xz_0[i] * pa_x[i] -
                          ta_xzz_xz_1[i] * pc_x[i];

        ta_xxzz_yy_0[i] = ta_zz_yy_0[i] * fe_0 - ta_zz_yy_1[i] * fe_0 + ta_xzz_yy_0[i] * pa_x[i] - ta_xzz_yy_1[i] * pc_x[i];

        ta_xxzz_yz_0[i] = ta_zz_yz_0[i] * fe_0 - ta_zz_yz_1[i] * fe_0 + ta_xzz_yz_0[i] * pa_x[i] - ta_xzz_yz_1[i] * pc_x[i];

        ta_xxzz_zz_0[i] = ta_zz_zz_0[i] * fe_0 - ta_zz_zz_1[i] * fe_0 + ta_xzz_zz_0[i] * pa_x[i] - ta_xzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto ta_xyyy_xx_0 = pbuffer.data(idx_npot_0_gd + 36);

    auto ta_xyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 37);

    auto ta_xyyy_xz_0 = pbuffer.data(idx_npot_0_gd + 38);

    auto ta_xyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 39);

    auto ta_xyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 40);

    auto ta_xyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 41);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xyyy_xx_0, \
                             ta_xyyy_xy_0, \
                             ta_xyyy_xz_0, \
                             ta_xyyy_yy_0, \
                             ta_xyyy_yz_0, \
                             ta_xyyy_zz_0, \
                             ta_yyy_x_0,   \
                             ta_yyy_x_1,   \
                             ta_yyy_xx_0,  \
                             ta_yyy_xx_1,  \
                             ta_yyy_xy_0,  \
                             ta_yyy_xy_1,  \
                             ta_yyy_xz_0,  \
                             ta_yyy_xz_1,  \
                             ta_yyy_y_0,   \
                             ta_yyy_y_1,   \
                             ta_yyy_yy_0,  \
                             ta_yyy_yy_1,  \
                             ta_yyy_yz_0,  \
                             ta_yyy_yz_1,  \
                             ta_yyy_z_0,   \
                             ta_yyy_z_1,   \
                             ta_yyy_zz_0,  \
                             ta_yyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyy_xx_0[i] = 2.0 * ta_yyy_x_0[i] * fe_0 - 2.0 * ta_yyy_x_1[i] * fe_0 + ta_yyy_xx_0[i] * pa_x[i] - ta_yyy_xx_1[i] * pc_x[i];

        ta_xyyy_xy_0[i] = ta_yyy_y_0[i] * fe_0 - ta_yyy_y_1[i] * fe_0 + ta_yyy_xy_0[i] * pa_x[i] - ta_yyy_xy_1[i] * pc_x[i];

        ta_xyyy_xz_0[i] = ta_yyy_z_0[i] * fe_0 - ta_yyy_z_1[i] * fe_0 + ta_yyy_xz_0[i] * pa_x[i] - ta_yyy_xz_1[i] * pc_x[i];

        ta_xyyy_yy_0[i] = ta_yyy_yy_0[i] * pa_x[i] - ta_yyy_yy_1[i] * pc_x[i];

        ta_xyyy_yz_0[i] = ta_yyy_yz_0[i] * pa_x[i] - ta_yyy_yz_1[i] * pc_x[i];

        ta_xyyy_zz_0[i] = ta_yyy_zz_0[i] * pa_x[i] - ta_yyy_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto ta_xyyz_xx_0 = pbuffer.data(idx_npot_0_gd + 42);

    auto ta_xyyz_xy_0 = pbuffer.data(idx_npot_0_gd + 43);

    auto ta_xyyz_xz_0 = pbuffer.data(idx_npot_0_gd + 44);

    auto ta_xyyz_yy_0 = pbuffer.data(idx_npot_0_gd + 45);

    auto ta_xyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 46);

    auto ta_xyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 47);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xyy_xx_0,  \
                             ta_xyy_xx_1,  \
                             ta_xyy_xy_0,  \
                             ta_xyy_xy_1,  \
                             ta_xyyz_xx_0, \
                             ta_xyyz_xy_0, \
                             ta_xyyz_xz_0, \
                             ta_xyyz_yy_0, \
                             ta_xyyz_yz_0, \
                             ta_xyyz_zz_0, \
                             ta_yyz_xz_0,  \
                             ta_yyz_xz_1,  \
                             ta_yyz_yy_0,  \
                             ta_yyz_yy_1,  \
                             ta_yyz_yz_0,  \
                             ta_yyz_yz_1,  \
                             ta_yyz_z_0,   \
                             ta_yyz_z_1,   \
                             ta_yyz_zz_0,  \
                             ta_yyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyz_xx_0[i] = ta_xyy_xx_0[i] * pa_z[i] - ta_xyy_xx_1[i] * pc_z[i];

        ta_xyyz_xy_0[i] = ta_xyy_xy_0[i] * pa_z[i] - ta_xyy_xy_1[i] * pc_z[i];

        ta_xyyz_xz_0[i] = ta_yyz_z_0[i] * fe_0 - ta_yyz_z_1[i] * fe_0 + ta_yyz_xz_0[i] * pa_x[i] - ta_yyz_xz_1[i] * pc_x[i];

        ta_xyyz_yy_0[i] = ta_yyz_yy_0[i] * pa_x[i] - ta_yyz_yy_1[i] * pc_x[i];

        ta_xyyz_yz_0[i] = ta_yyz_yz_0[i] * pa_x[i] - ta_yyz_yz_1[i] * pc_x[i];

        ta_xyyz_zz_0[i] = ta_yyz_zz_0[i] * pa_x[i] - ta_yyz_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto ta_xyzz_xx_0 = pbuffer.data(idx_npot_0_gd + 48);

    auto ta_xyzz_xy_0 = pbuffer.data(idx_npot_0_gd + 49);

    auto ta_xyzz_xz_0 = pbuffer.data(idx_npot_0_gd + 50);

    auto ta_xyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 51);

    auto ta_xyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 52);

    auto ta_xyzz_zz_0 = pbuffer.data(idx_npot_0_gd + 53);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xyzz_xx_0, \
                             ta_xyzz_xy_0, \
                             ta_xyzz_xz_0, \
                             ta_xyzz_yy_0, \
                             ta_xyzz_yz_0, \
                             ta_xyzz_zz_0, \
                             ta_xzz_xx_0,  \
                             ta_xzz_xx_1,  \
                             ta_xzz_xz_0,  \
                             ta_xzz_xz_1,  \
                             ta_yzz_xy_0,  \
                             ta_yzz_xy_1,  \
                             ta_yzz_y_0,   \
                             ta_yzz_y_1,   \
                             ta_yzz_yy_0,  \
                             ta_yzz_yy_1,  \
                             ta_yzz_yz_0,  \
                             ta_yzz_yz_1,  \
                             ta_yzz_zz_0,  \
                             ta_yzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzz_xx_0[i] = ta_xzz_xx_0[i] * pa_y[i] - ta_xzz_xx_1[i] * pc_y[i];

        ta_xyzz_xy_0[i] = ta_yzz_y_0[i] * fe_0 - ta_yzz_y_1[i] * fe_0 + ta_yzz_xy_0[i] * pa_x[i] - ta_yzz_xy_1[i] * pc_x[i];

        ta_xyzz_xz_0[i] = ta_xzz_xz_0[i] * pa_y[i] - ta_xzz_xz_1[i] * pc_y[i];

        ta_xyzz_yy_0[i] = ta_yzz_yy_0[i] * pa_x[i] - ta_yzz_yy_1[i] * pc_x[i];

        ta_xyzz_yz_0[i] = ta_yzz_yz_0[i] * pa_x[i] - ta_yzz_yz_1[i] * pc_x[i];

        ta_xyzz_zz_0[i] = ta_yzz_zz_0[i] * pa_x[i] - ta_yzz_zz_1[i] * pc_x[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto ta_xzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 54);

    auto ta_xzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 55);

    auto ta_xzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 56);

    auto ta_xzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 57);

    auto ta_xzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 58);

    auto ta_xzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 59);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xzzz_xx_0, \
                             ta_xzzz_xy_0, \
                             ta_xzzz_xz_0, \
                             ta_xzzz_yy_0, \
                             ta_xzzz_yz_0, \
                             ta_xzzz_zz_0, \
                             ta_zzz_x_0,   \
                             ta_zzz_x_1,   \
                             ta_zzz_xx_0,  \
                             ta_zzz_xx_1,  \
                             ta_zzz_xy_0,  \
                             ta_zzz_xy_1,  \
                             ta_zzz_xz_0,  \
                             ta_zzz_xz_1,  \
                             ta_zzz_y_0,   \
                             ta_zzz_y_1,   \
                             ta_zzz_yy_0,  \
                             ta_zzz_yy_1,  \
                             ta_zzz_yz_0,  \
                             ta_zzz_yz_1,  \
                             ta_zzz_z_0,   \
                             ta_zzz_z_1,   \
                             ta_zzz_zz_0,  \
                             ta_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzz_xx_0[i] = 2.0 * ta_zzz_x_0[i] * fe_0 - 2.0 * ta_zzz_x_1[i] * fe_0 + ta_zzz_xx_0[i] * pa_x[i] - ta_zzz_xx_1[i] * pc_x[i];

        ta_xzzz_xy_0[i] = ta_zzz_y_0[i] * fe_0 - ta_zzz_y_1[i] * fe_0 + ta_zzz_xy_0[i] * pa_x[i] - ta_zzz_xy_1[i] * pc_x[i];

        ta_xzzz_xz_0[i] = ta_zzz_z_0[i] * fe_0 - ta_zzz_z_1[i] * fe_0 + ta_zzz_xz_0[i] * pa_x[i] - ta_zzz_xz_1[i] * pc_x[i];

        ta_xzzz_yy_0[i] = ta_zzz_yy_0[i] * pa_x[i] - ta_zzz_yy_1[i] * pc_x[i];

        ta_xzzz_yz_0[i] = ta_zzz_yz_0[i] * pa_x[i] - ta_zzz_yz_1[i] * pc_x[i];

        ta_xzzz_zz_0[i] = ta_zzz_zz_0[i] * pa_x[i] - ta_zzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto ta_yyyy_xx_0 = pbuffer.data(idx_npot_0_gd + 60);

    auto ta_yyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 61);

    auto ta_yyyy_xz_0 = pbuffer.data(idx_npot_0_gd + 62);

    auto ta_yyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 63);

    auto ta_yyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 64);

    auto ta_yyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 65);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_yy_xx_0,   \
                             ta_yy_xx_1,   \
                             ta_yy_xy_0,   \
                             ta_yy_xy_1,   \
                             ta_yy_xz_0,   \
                             ta_yy_xz_1,   \
                             ta_yy_yy_0,   \
                             ta_yy_yy_1,   \
                             ta_yy_yz_0,   \
                             ta_yy_yz_1,   \
                             ta_yy_zz_0,   \
                             ta_yy_zz_1,   \
                             ta_yyy_x_0,   \
                             ta_yyy_x_1,   \
                             ta_yyy_xx_0,  \
                             ta_yyy_xx_1,  \
                             ta_yyy_xy_0,  \
                             ta_yyy_xy_1,  \
                             ta_yyy_xz_0,  \
                             ta_yyy_xz_1,  \
                             ta_yyy_y_0,   \
                             ta_yyy_y_1,   \
                             ta_yyy_yy_0,  \
                             ta_yyy_yy_1,  \
                             ta_yyy_yz_0,  \
                             ta_yyy_yz_1,  \
                             ta_yyy_z_0,   \
                             ta_yyy_z_1,   \
                             ta_yyy_zz_0,  \
                             ta_yyy_zz_1,  \
                             ta_yyyy_xx_0, \
                             ta_yyyy_xy_0, \
                             ta_yyyy_xz_0, \
                             ta_yyyy_yy_0, \
                             ta_yyyy_yz_0, \
                             ta_yyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyy_xx_0[i] = 3.0 * ta_yy_xx_0[i] * fe_0 - 3.0 * ta_yy_xx_1[i] * fe_0 + ta_yyy_xx_0[i] * pa_y[i] - ta_yyy_xx_1[i] * pc_y[i];

        ta_yyyy_xy_0[i] = 3.0 * ta_yy_xy_0[i] * fe_0 - 3.0 * ta_yy_xy_1[i] * fe_0 + ta_yyy_x_0[i] * fe_0 - ta_yyy_x_1[i] * fe_0 +
                          ta_yyy_xy_0[i] * pa_y[i] - ta_yyy_xy_1[i] * pc_y[i];

        ta_yyyy_xz_0[i] = 3.0 * ta_yy_xz_0[i] * fe_0 - 3.0 * ta_yy_xz_1[i] * fe_0 + ta_yyy_xz_0[i] * pa_y[i] - ta_yyy_xz_1[i] * pc_y[i];

        ta_yyyy_yy_0[i] = 3.0 * ta_yy_yy_0[i] * fe_0 - 3.0 * ta_yy_yy_1[i] * fe_0 + 2.0 * ta_yyy_y_0[i] * fe_0 - 2.0 * ta_yyy_y_1[i] * fe_0 +
                          ta_yyy_yy_0[i] * pa_y[i] - ta_yyy_yy_1[i] * pc_y[i];

        ta_yyyy_yz_0[i] = 3.0 * ta_yy_yz_0[i] * fe_0 - 3.0 * ta_yy_yz_1[i] * fe_0 + ta_yyy_z_0[i] * fe_0 - ta_yyy_z_1[i] * fe_0 +
                          ta_yyy_yz_0[i] * pa_y[i] - ta_yyy_yz_1[i] * pc_y[i];

        ta_yyyy_zz_0[i] = 3.0 * ta_yy_zz_0[i] * fe_0 - 3.0 * ta_yy_zz_1[i] * fe_0 + ta_yyy_zz_0[i] * pa_y[i] - ta_yyy_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto ta_yyyz_xx_0 = pbuffer.data(idx_npot_0_gd + 66);

    auto ta_yyyz_xy_0 = pbuffer.data(idx_npot_0_gd + 67);

    auto ta_yyyz_xz_0 = pbuffer.data(idx_npot_0_gd + 68);

    auto ta_yyyz_yy_0 = pbuffer.data(idx_npot_0_gd + 69);

    auto ta_yyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 70);

    auto ta_yyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 71);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_yyy_xx_0,  \
                             ta_yyy_xx_1,  \
                             ta_yyy_xy_0,  \
                             ta_yyy_xy_1,  \
                             ta_yyy_y_0,   \
                             ta_yyy_y_1,   \
                             ta_yyy_yy_0,  \
                             ta_yyy_yy_1,  \
                             ta_yyy_yz_0,  \
                             ta_yyy_yz_1,  \
                             ta_yyyz_xx_0, \
                             ta_yyyz_xy_0, \
                             ta_yyyz_xz_0, \
                             ta_yyyz_yy_0, \
                             ta_yyyz_yz_0, \
                             ta_yyyz_zz_0, \
                             ta_yyz_xz_0,  \
                             ta_yyz_xz_1,  \
                             ta_yyz_zz_0,  \
                             ta_yyz_zz_1,  \
                             ta_yz_xz_0,   \
                             ta_yz_xz_1,   \
                             ta_yz_zz_0,   \
                             ta_yz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyz_xx_0[i] = ta_yyy_xx_0[i] * pa_z[i] - ta_yyy_xx_1[i] * pc_z[i];

        ta_yyyz_xy_0[i] = ta_yyy_xy_0[i] * pa_z[i] - ta_yyy_xy_1[i] * pc_z[i];

        ta_yyyz_xz_0[i] = 2.0 * ta_yz_xz_0[i] * fe_0 - 2.0 * ta_yz_xz_1[i] * fe_0 + ta_yyz_xz_0[i] * pa_y[i] - ta_yyz_xz_1[i] * pc_y[i];

        ta_yyyz_yy_0[i] = ta_yyy_yy_0[i] * pa_z[i] - ta_yyy_yy_1[i] * pc_z[i];

        ta_yyyz_yz_0[i] = ta_yyy_y_0[i] * fe_0 - ta_yyy_y_1[i] * fe_0 + ta_yyy_yz_0[i] * pa_z[i] - ta_yyy_yz_1[i] * pc_z[i];

        ta_yyyz_zz_0[i] = 2.0 * ta_yz_zz_0[i] * fe_0 - 2.0 * ta_yz_zz_1[i] * fe_0 + ta_yyz_zz_0[i] * pa_y[i] - ta_yyz_zz_1[i] * pc_y[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto ta_yyzz_xx_0 = pbuffer.data(idx_npot_0_gd + 72);

    auto ta_yyzz_xy_0 = pbuffer.data(idx_npot_0_gd + 73);

    auto ta_yyzz_xz_0 = pbuffer.data(idx_npot_0_gd + 74);

    auto ta_yyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 75);

    auto ta_yyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 76);

    auto ta_yyzz_zz_0 = pbuffer.data(idx_npot_0_gd + 77);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_yy_xy_0,   \
                             ta_yy_xy_1,   \
                             ta_yy_yy_0,   \
                             ta_yy_yy_1,   \
                             ta_yyz_xy_0,  \
                             ta_yyz_xy_1,  \
                             ta_yyz_yy_0,  \
                             ta_yyz_yy_1,  \
                             ta_yyzz_xx_0, \
                             ta_yyzz_xy_0, \
                             ta_yyzz_xz_0, \
                             ta_yyzz_yy_0, \
                             ta_yyzz_yz_0, \
                             ta_yyzz_zz_0, \
                             ta_yzz_xx_0,  \
                             ta_yzz_xx_1,  \
                             ta_yzz_xz_0,  \
                             ta_yzz_xz_1,  \
                             ta_yzz_yz_0,  \
                             ta_yzz_yz_1,  \
                             ta_yzz_z_0,   \
                             ta_yzz_z_1,   \
                             ta_yzz_zz_0,  \
                             ta_yzz_zz_1,  \
                             ta_zz_xx_0,   \
                             ta_zz_xx_1,   \
                             ta_zz_xz_0,   \
                             ta_zz_xz_1,   \
                             ta_zz_yz_0,   \
                             ta_zz_yz_1,   \
                             ta_zz_zz_0,   \
                             ta_zz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzz_xx_0[i] = ta_zz_xx_0[i] * fe_0 - ta_zz_xx_1[i] * fe_0 + ta_yzz_xx_0[i] * pa_y[i] - ta_yzz_xx_1[i] * pc_y[i];

        ta_yyzz_xy_0[i] = ta_yy_xy_0[i] * fe_0 - ta_yy_xy_1[i] * fe_0 + ta_yyz_xy_0[i] * pa_z[i] - ta_yyz_xy_1[i] * pc_z[i];

        ta_yyzz_xz_0[i] = ta_zz_xz_0[i] * fe_0 - ta_zz_xz_1[i] * fe_0 + ta_yzz_xz_0[i] * pa_y[i] - ta_yzz_xz_1[i] * pc_y[i];

        ta_yyzz_yy_0[i] = ta_yy_yy_0[i] * fe_0 - ta_yy_yy_1[i] * fe_0 + ta_yyz_yy_0[i] * pa_z[i] - ta_yyz_yy_1[i] * pc_z[i];

        ta_yyzz_yz_0[i] = ta_zz_yz_0[i] * fe_0 - ta_zz_yz_1[i] * fe_0 + ta_yzz_z_0[i] * fe_0 - ta_yzz_z_1[i] * fe_0 + ta_yzz_yz_0[i] * pa_y[i] -
                          ta_yzz_yz_1[i] * pc_y[i];

        ta_yyzz_zz_0[i] = ta_zz_zz_0[i] * fe_0 - ta_zz_zz_1[i] * fe_0 + ta_yzz_zz_0[i] * pa_y[i] - ta_yzz_zz_1[i] * pc_y[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto ta_yzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 78);

    auto ta_yzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 79);

    auto ta_yzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 80);

    auto ta_yzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 81);

    auto ta_yzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 82);

    auto ta_yzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 83);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_yzzz_xx_0, \
                             ta_yzzz_xy_0, \
                             ta_yzzz_xz_0, \
                             ta_yzzz_yy_0, \
                             ta_yzzz_yz_0, \
                             ta_yzzz_zz_0, \
                             ta_zzz_x_0,   \
                             ta_zzz_x_1,   \
                             ta_zzz_xx_0,  \
                             ta_zzz_xx_1,  \
                             ta_zzz_xy_0,  \
                             ta_zzz_xy_1,  \
                             ta_zzz_xz_0,  \
                             ta_zzz_xz_1,  \
                             ta_zzz_y_0,   \
                             ta_zzz_y_1,   \
                             ta_zzz_yy_0,  \
                             ta_zzz_yy_1,  \
                             ta_zzz_yz_0,  \
                             ta_zzz_yz_1,  \
                             ta_zzz_z_0,   \
                             ta_zzz_z_1,   \
                             ta_zzz_zz_0,  \
                             ta_zzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzz_xx_0[i] = ta_zzz_xx_0[i] * pa_y[i] - ta_zzz_xx_1[i] * pc_y[i];

        ta_yzzz_xy_0[i] = ta_zzz_x_0[i] * fe_0 - ta_zzz_x_1[i] * fe_0 + ta_zzz_xy_0[i] * pa_y[i] - ta_zzz_xy_1[i] * pc_y[i];

        ta_yzzz_xz_0[i] = ta_zzz_xz_0[i] * pa_y[i] - ta_zzz_xz_1[i] * pc_y[i];

        ta_yzzz_yy_0[i] = 2.0 * ta_zzz_y_0[i] * fe_0 - 2.0 * ta_zzz_y_1[i] * fe_0 + ta_zzz_yy_0[i] * pa_y[i] - ta_zzz_yy_1[i] * pc_y[i];

        ta_yzzz_yz_0[i] = ta_zzz_z_0[i] * fe_0 - ta_zzz_z_1[i] * fe_0 + ta_zzz_yz_0[i] * pa_y[i] - ta_zzz_yz_1[i] * pc_y[i];

        ta_yzzz_zz_0[i] = ta_zzz_zz_0[i] * pa_y[i] - ta_zzz_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto ta_zzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 84);

    auto ta_zzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 85);

    auto ta_zzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 86);

    auto ta_zzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 87);

    auto ta_zzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 88);

    auto ta_zzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 89);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta_zz_xx_0,   \
                             ta_zz_xx_1,   \
                             ta_zz_xy_0,   \
                             ta_zz_xy_1,   \
                             ta_zz_xz_0,   \
                             ta_zz_xz_1,   \
                             ta_zz_yy_0,   \
                             ta_zz_yy_1,   \
                             ta_zz_yz_0,   \
                             ta_zz_yz_1,   \
                             ta_zz_zz_0,   \
                             ta_zz_zz_1,   \
                             ta_zzz_x_0,   \
                             ta_zzz_x_1,   \
                             ta_zzz_xx_0,  \
                             ta_zzz_xx_1,  \
                             ta_zzz_xy_0,  \
                             ta_zzz_xy_1,  \
                             ta_zzz_xz_0,  \
                             ta_zzz_xz_1,  \
                             ta_zzz_y_0,   \
                             ta_zzz_y_1,   \
                             ta_zzz_yy_0,  \
                             ta_zzz_yy_1,  \
                             ta_zzz_yz_0,  \
                             ta_zzz_yz_1,  \
                             ta_zzz_z_0,   \
                             ta_zzz_z_1,   \
                             ta_zzz_zz_0,  \
                             ta_zzz_zz_1,  \
                             ta_zzzz_xx_0, \
                             ta_zzzz_xy_0, \
                             ta_zzzz_xz_0, \
                             ta_zzzz_yy_0, \
                             ta_zzzz_yz_0, \
                             ta_zzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzz_xx_0[i] = 3.0 * ta_zz_xx_0[i] * fe_0 - 3.0 * ta_zz_xx_1[i] * fe_0 + ta_zzz_xx_0[i] * pa_z[i] - ta_zzz_xx_1[i] * pc_z[i];

        ta_zzzz_xy_0[i] = 3.0 * ta_zz_xy_0[i] * fe_0 - 3.0 * ta_zz_xy_1[i] * fe_0 + ta_zzz_xy_0[i] * pa_z[i] - ta_zzz_xy_1[i] * pc_z[i];

        ta_zzzz_xz_0[i] = 3.0 * ta_zz_xz_0[i] * fe_0 - 3.0 * ta_zz_xz_1[i] * fe_0 + ta_zzz_x_0[i] * fe_0 - ta_zzz_x_1[i] * fe_0 +
                          ta_zzz_xz_0[i] * pa_z[i] - ta_zzz_xz_1[i] * pc_z[i];

        ta_zzzz_yy_0[i] = 3.0 * ta_zz_yy_0[i] * fe_0 - 3.0 * ta_zz_yy_1[i] * fe_0 + ta_zzz_yy_0[i] * pa_z[i] - ta_zzz_yy_1[i] * pc_z[i];

        ta_zzzz_yz_0[i] = 3.0 * ta_zz_yz_0[i] * fe_0 - 3.0 * ta_zz_yz_1[i] * fe_0 + ta_zzz_y_0[i] * fe_0 - ta_zzz_y_1[i] * fe_0 +
                          ta_zzz_yz_0[i] * pa_z[i] - ta_zzz_yz_1[i] * pc_z[i];

        ta_zzzz_zz_0[i] = 3.0 * ta_zz_zz_0[i] * fe_0 - 3.0 * ta_zz_zz_1[i] * fe_0 + 2.0 * ta_zzz_z_0[i] * fe_0 - 2.0 * ta_zzz_z_1[i] * fe_0 +
                          ta_zzz_zz_0[i] * pa_z[i] - ta_zzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
