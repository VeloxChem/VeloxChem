#include "NuclearPotentialPrimRecGP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_gp(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_gp,
                               const size_t              idx_npot_0_dp,
                               const size_t              idx_npot_1_dp,
                               const size_t              idx_npot_0_fs,
                               const size_t              idx_npot_1_fs,
                               const size_t              idx_npot_0_fp,
                               const size_t              idx_npot_1_fp,
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

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_0 = pbuffer.data(idx_npot_0_dp);

    auto ta_xx_y_0 = pbuffer.data(idx_npot_0_dp + 1);

    auto ta_xx_z_0 = pbuffer.data(idx_npot_0_dp + 2);

    auto ta_xy_y_0 = pbuffer.data(idx_npot_0_dp + 4);

    auto ta_xz_z_0 = pbuffer.data(idx_npot_0_dp + 8);

    auto ta_yy_x_0 = pbuffer.data(idx_npot_0_dp + 9);

    auto ta_yy_y_0 = pbuffer.data(idx_npot_0_dp + 10);

    auto ta_yy_z_0 = pbuffer.data(idx_npot_0_dp + 11);

    auto ta_yz_z_0 = pbuffer.data(idx_npot_0_dp + 14);

    auto ta_zz_x_0 = pbuffer.data(idx_npot_0_dp + 15);

    auto ta_zz_y_0 = pbuffer.data(idx_npot_0_dp + 16);

    auto ta_zz_z_0 = pbuffer.data(idx_npot_0_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_1 = pbuffer.data(idx_npot_1_dp);

    auto ta_xx_y_1 = pbuffer.data(idx_npot_1_dp + 1);

    auto ta_xx_z_1 = pbuffer.data(idx_npot_1_dp + 2);

    auto ta_xy_y_1 = pbuffer.data(idx_npot_1_dp + 4);

    auto ta_xz_z_1 = pbuffer.data(idx_npot_1_dp + 8);

    auto ta_yy_x_1 = pbuffer.data(idx_npot_1_dp + 9);

    auto ta_yy_y_1 = pbuffer.data(idx_npot_1_dp + 10);

    auto ta_yy_z_1 = pbuffer.data(idx_npot_1_dp + 11);

    auto ta_yz_z_1 = pbuffer.data(idx_npot_1_dp + 14);

    auto ta_zz_x_1 = pbuffer.data(idx_npot_1_dp + 15);

    auto ta_zz_y_1 = pbuffer.data(idx_npot_1_dp + 16);

    auto ta_zz_z_1 = pbuffer.data(idx_npot_1_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_0 = pbuffer.data(idx_npot_0_fs);

    auto ta_yyy_0_0 = pbuffer.data(idx_npot_0_fs + 6);

    auto ta_zzz_0_0 = pbuffer.data(idx_npot_0_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto ta_xxx_0_1 = pbuffer.data(idx_npot_1_fs);

    auto ta_yyy_0_1 = pbuffer.data(idx_npot_1_fs + 6);

    auto ta_zzz_0_1 = pbuffer.data(idx_npot_1_fs + 9);

    // Set up components of auxiliary buffer : FP

    auto ta_xxx_x_0 = pbuffer.data(idx_npot_0_fp);

    auto ta_xxx_y_0 = pbuffer.data(idx_npot_0_fp + 1);

    auto ta_xxx_z_0 = pbuffer.data(idx_npot_0_fp + 2);

    auto ta_xxy_x_0 = pbuffer.data(idx_npot_0_fp + 3);

    auto ta_xxy_y_0 = pbuffer.data(idx_npot_0_fp + 4);

    auto ta_xxz_x_0 = pbuffer.data(idx_npot_0_fp + 6);

    auto ta_xxz_z_0 = pbuffer.data(idx_npot_0_fp + 8);

    auto ta_xyy_x_0 = pbuffer.data(idx_npot_0_fp + 9);

    auto ta_xyy_y_0 = pbuffer.data(idx_npot_0_fp + 10);

    auto ta_xyy_z_0 = pbuffer.data(idx_npot_0_fp + 11);

    auto ta_xzz_x_0 = pbuffer.data(idx_npot_0_fp + 15);

    auto ta_xzz_y_0 = pbuffer.data(idx_npot_0_fp + 16);

    auto ta_xzz_z_0 = pbuffer.data(idx_npot_0_fp + 17);

    auto ta_yyy_x_0 = pbuffer.data(idx_npot_0_fp + 18);

    auto ta_yyy_y_0 = pbuffer.data(idx_npot_0_fp + 19);

    auto ta_yyy_z_0 = pbuffer.data(idx_npot_0_fp + 20);

    auto ta_yyz_y_0 = pbuffer.data(idx_npot_0_fp + 22);

    auto ta_yyz_z_0 = pbuffer.data(idx_npot_0_fp + 23);

    auto ta_yzz_x_0 = pbuffer.data(idx_npot_0_fp + 24);

    auto ta_yzz_y_0 = pbuffer.data(idx_npot_0_fp + 25);

    auto ta_yzz_z_0 = pbuffer.data(idx_npot_0_fp + 26);

    auto ta_zzz_x_0 = pbuffer.data(idx_npot_0_fp + 27);

    auto ta_zzz_y_0 = pbuffer.data(idx_npot_0_fp + 28);

    auto ta_zzz_z_0 = pbuffer.data(idx_npot_0_fp + 29);

    // Set up components of auxiliary buffer : FP

    auto ta_xxx_x_1 = pbuffer.data(idx_npot_1_fp);

    auto ta_xxx_y_1 = pbuffer.data(idx_npot_1_fp + 1);

    auto ta_xxx_z_1 = pbuffer.data(idx_npot_1_fp + 2);

    auto ta_xxy_x_1 = pbuffer.data(idx_npot_1_fp + 3);

    auto ta_xxy_y_1 = pbuffer.data(idx_npot_1_fp + 4);

    auto ta_xxz_x_1 = pbuffer.data(idx_npot_1_fp + 6);

    auto ta_xxz_z_1 = pbuffer.data(idx_npot_1_fp + 8);

    auto ta_xyy_x_1 = pbuffer.data(idx_npot_1_fp + 9);

    auto ta_xyy_y_1 = pbuffer.data(idx_npot_1_fp + 10);

    auto ta_xyy_z_1 = pbuffer.data(idx_npot_1_fp + 11);

    auto ta_xzz_x_1 = pbuffer.data(idx_npot_1_fp + 15);

    auto ta_xzz_y_1 = pbuffer.data(idx_npot_1_fp + 16);

    auto ta_xzz_z_1 = pbuffer.data(idx_npot_1_fp + 17);

    auto ta_yyy_x_1 = pbuffer.data(idx_npot_1_fp + 18);

    auto ta_yyy_y_1 = pbuffer.data(idx_npot_1_fp + 19);

    auto ta_yyy_z_1 = pbuffer.data(idx_npot_1_fp + 20);

    auto ta_yyz_y_1 = pbuffer.data(idx_npot_1_fp + 22);

    auto ta_yyz_z_1 = pbuffer.data(idx_npot_1_fp + 23);

    auto ta_yzz_x_1 = pbuffer.data(idx_npot_1_fp + 24);

    auto ta_yzz_y_1 = pbuffer.data(idx_npot_1_fp + 25);

    auto ta_yzz_z_1 = pbuffer.data(idx_npot_1_fp + 26);

    auto ta_zzz_x_1 = pbuffer.data(idx_npot_1_fp + 27);

    auto ta_zzz_y_1 = pbuffer.data(idx_npot_1_fp + 28);

    auto ta_zzz_z_1 = pbuffer.data(idx_npot_1_fp + 29);

    // Set up 0-3 components of targeted buffer : GP

    auto ta_xxxx_x_0 = pbuffer.data(idx_npot_0_gp);

    auto ta_xxxx_y_0 = pbuffer.data(idx_npot_0_gp + 1);

    auto ta_xxxx_z_0 = pbuffer.data(idx_npot_0_gp + 2);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_xx_x_0,   \
                             ta_xx_x_1,   \
                             ta_xx_y_0,   \
                             ta_xx_y_1,   \
                             ta_xx_z_0,   \
                             ta_xx_z_1,   \
                             ta_xxx_0_0,  \
                             ta_xxx_0_1,  \
                             ta_xxx_x_0,  \
                             ta_xxx_x_1,  \
                             ta_xxx_y_0,  \
                             ta_xxx_y_1,  \
                             ta_xxx_z_0,  \
                             ta_xxx_z_1,  \
                             ta_xxxx_x_0, \
                             ta_xxxx_y_0, \
                             ta_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_x_0[i] = 3.0 * ta_xx_x_0[i] * fe_0 - 3.0 * ta_xx_x_1[i] * fe_0 + ta_xxx_0_0[i] * fe_0 -
                         ta_xxx_0_1[i] * fe_0 + ta_xxx_x_0[i] * pa_x[i] - ta_xxx_x_1[i] * pc_x[i];

        ta_xxxx_y_0[i] =
            3.0 * ta_xx_y_0[i] * fe_0 - 3.0 * ta_xx_y_1[i] * fe_0 + ta_xxx_y_0[i] * pa_x[i] - ta_xxx_y_1[i] * pc_x[i];

        ta_xxxx_z_0[i] =
            3.0 * ta_xx_z_0[i] * fe_0 - 3.0 * ta_xx_z_1[i] * fe_0 + ta_xxx_z_0[i] * pa_x[i] - ta_xxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto ta_xxxy_x_0 = pbuffer.data(idx_npot_0_gp + 3);

    auto ta_xxxy_y_0 = pbuffer.data(idx_npot_0_gp + 4);

    auto ta_xxxy_z_0 = pbuffer.data(idx_npot_0_gp + 5);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pc_x,        \
                             pc_y,        \
                             ta_xxx_x_0,  \
                             ta_xxx_x_1,  \
                             ta_xxx_z_0,  \
                             ta_xxx_z_1,  \
                             ta_xxxy_x_0, \
                             ta_xxxy_y_0, \
                             ta_xxxy_z_0, \
                             ta_xxy_y_0,  \
                             ta_xxy_y_1,  \
                             ta_xy_y_0,   \
                             ta_xy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxy_x_0[i] = ta_xxx_x_0[i] * pa_y[i] - ta_xxx_x_1[i] * pc_y[i];

        ta_xxxy_y_0[i] =
            2.0 * ta_xy_y_0[i] * fe_0 - 2.0 * ta_xy_y_1[i] * fe_0 + ta_xxy_y_0[i] * pa_x[i] - ta_xxy_y_1[i] * pc_x[i];

        ta_xxxy_z_0[i] = ta_xxx_z_0[i] * pa_y[i] - ta_xxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto ta_xxxz_x_0 = pbuffer.data(idx_npot_0_gp + 6);

    auto ta_xxxz_y_0 = pbuffer.data(idx_npot_0_gp + 7);

    auto ta_xxxz_z_0 = pbuffer.data(idx_npot_0_gp + 8);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             pc_x,        \
                             pc_z,        \
                             ta_xxx_x_0,  \
                             ta_xxx_x_1,  \
                             ta_xxx_y_0,  \
                             ta_xxx_y_1,  \
                             ta_xxxz_x_0, \
                             ta_xxxz_y_0, \
                             ta_xxxz_z_0, \
                             ta_xxz_z_0,  \
                             ta_xxz_z_1,  \
                             ta_xz_z_0,   \
                             ta_xz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxz_x_0[i] = ta_xxx_x_0[i] * pa_z[i] - ta_xxx_x_1[i] * pc_z[i];

        ta_xxxz_y_0[i] = ta_xxx_y_0[i] * pa_z[i] - ta_xxx_y_1[i] * pc_z[i];

        ta_xxxz_z_0[i] =
            2.0 * ta_xz_z_0[i] * fe_0 - 2.0 * ta_xz_z_1[i] * fe_0 + ta_xxz_z_0[i] * pa_x[i] - ta_xxz_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto ta_xxyy_x_0 = pbuffer.data(idx_npot_0_gp + 9);

    auto ta_xxyy_y_0 = pbuffer.data(idx_npot_0_gp + 10);

    auto ta_xxyy_z_0 = pbuffer.data(idx_npot_0_gp + 11);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pc_x,        \
                             pc_y,        \
                             ta_xx_x_0,   \
                             ta_xx_x_1,   \
                             ta_xxy_x_0,  \
                             ta_xxy_x_1,  \
                             ta_xxyy_x_0, \
                             ta_xxyy_y_0, \
                             ta_xxyy_z_0, \
                             ta_xyy_y_0,  \
                             ta_xyy_y_1,  \
                             ta_xyy_z_0,  \
                             ta_xyy_z_1,  \
                             ta_yy_y_0,   \
                             ta_yy_y_1,   \
                             ta_yy_z_0,   \
                             ta_yy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyy_x_0[i] = ta_xx_x_0[i] * fe_0 - ta_xx_x_1[i] * fe_0 + ta_xxy_x_0[i] * pa_y[i] - ta_xxy_x_1[i] * pc_y[i];

        ta_xxyy_y_0[i] = ta_yy_y_0[i] * fe_0 - ta_yy_y_1[i] * fe_0 + ta_xyy_y_0[i] * pa_x[i] - ta_xyy_y_1[i] * pc_x[i];

        ta_xxyy_z_0[i] = ta_yy_z_0[i] * fe_0 - ta_yy_z_1[i] * fe_0 + ta_xyy_z_0[i] * pa_x[i] - ta_xyy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto ta_xxyz_x_0 = pbuffer.data(idx_npot_0_gp + 12);

    auto ta_xxyz_y_0 = pbuffer.data(idx_npot_0_gp + 13);

    auto ta_xxyz_z_0 = pbuffer.data(idx_npot_0_gp + 14);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             pc_y,        \
                             pc_z,        \
                             ta_xxy_y_0,  \
                             ta_xxy_y_1,  \
                             ta_xxyz_x_0, \
                             ta_xxyz_y_0, \
                             ta_xxyz_z_0, \
                             ta_xxz_x_0,  \
                             ta_xxz_x_1,  \
                             ta_xxz_z_0,  \
                             ta_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xxyz_x_0[i] = ta_xxz_x_0[i] * pa_y[i] - ta_xxz_x_1[i] * pc_y[i];

        ta_xxyz_y_0[i] = ta_xxy_y_0[i] * pa_z[i] - ta_xxy_y_1[i] * pc_z[i];

        ta_xxyz_z_0[i] = ta_xxz_z_0[i] * pa_y[i] - ta_xxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto ta_xxzz_x_0 = pbuffer.data(idx_npot_0_gp + 15);

    auto ta_xxzz_y_0 = pbuffer.data(idx_npot_0_gp + 16);

    auto ta_xxzz_z_0 = pbuffer.data(idx_npot_0_gp + 17);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             pc_x,        \
                             pc_z,        \
                             ta_xx_x_0,   \
                             ta_xx_x_1,   \
                             ta_xxz_x_0,  \
                             ta_xxz_x_1,  \
                             ta_xxzz_x_0, \
                             ta_xxzz_y_0, \
                             ta_xxzz_z_0, \
                             ta_xzz_y_0,  \
                             ta_xzz_y_1,  \
                             ta_xzz_z_0,  \
                             ta_xzz_z_1,  \
                             ta_zz_y_0,   \
                             ta_zz_y_1,   \
                             ta_zz_z_0,   \
                             ta_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzz_x_0[i] = ta_xx_x_0[i] * fe_0 - ta_xx_x_1[i] * fe_0 + ta_xxz_x_0[i] * pa_z[i] - ta_xxz_x_1[i] * pc_z[i];

        ta_xxzz_y_0[i] = ta_zz_y_0[i] * fe_0 - ta_zz_y_1[i] * fe_0 + ta_xzz_y_0[i] * pa_x[i] - ta_xzz_y_1[i] * pc_x[i];

        ta_xxzz_z_0[i] = ta_zz_z_0[i] * fe_0 - ta_zz_z_1[i] * fe_0 + ta_xzz_z_0[i] * pa_x[i] - ta_xzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto ta_xyyy_x_0 = pbuffer.data(idx_npot_0_gp + 18);

    auto ta_xyyy_y_0 = pbuffer.data(idx_npot_0_gp + 19);

    auto ta_xyyy_z_0 = pbuffer.data(idx_npot_0_gp + 20);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_xyyy_x_0, \
                             ta_xyyy_y_0, \
                             ta_xyyy_z_0, \
                             ta_yyy_0_0,  \
                             ta_yyy_0_1,  \
                             ta_yyy_x_0,  \
                             ta_yyy_x_1,  \
                             ta_yyy_y_0,  \
                             ta_yyy_y_1,  \
                             ta_yyy_z_0,  \
                             ta_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyy_x_0[i] =
            ta_yyy_0_0[i] * fe_0 - ta_yyy_0_1[i] * fe_0 + ta_yyy_x_0[i] * pa_x[i] - ta_yyy_x_1[i] * pc_x[i];

        ta_xyyy_y_0[i] = ta_yyy_y_0[i] * pa_x[i] - ta_yyy_y_1[i] * pc_x[i];

        ta_xyyy_z_0[i] = ta_yyy_z_0[i] * pa_x[i] - ta_yyy_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto ta_xyyz_x_0 = pbuffer.data(idx_npot_0_gp + 21);

    auto ta_xyyz_y_0 = pbuffer.data(idx_npot_0_gp + 22);

    auto ta_xyyz_z_0 = pbuffer.data(idx_npot_0_gp + 23);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             pc_x,        \
                             pc_z,        \
                             ta_xyy_x_0,  \
                             ta_xyy_x_1,  \
                             ta_xyyz_x_0, \
                             ta_xyyz_y_0, \
                             ta_xyyz_z_0, \
                             ta_yyz_y_0,  \
                             ta_yyz_y_1,  \
                             ta_yyz_z_0,  \
                             ta_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyyz_x_0[i] = ta_xyy_x_0[i] * pa_z[i] - ta_xyy_x_1[i] * pc_z[i];

        ta_xyyz_y_0[i] = ta_yyz_y_0[i] * pa_x[i] - ta_yyz_y_1[i] * pc_x[i];

        ta_xyyz_z_0[i] = ta_yyz_z_0[i] * pa_x[i] - ta_yyz_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto ta_xyzz_x_0 = pbuffer.data(idx_npot_0_gp + 24);

    auto ta_xyzz_y_0 = pbuffer.data(idx_npot_0_gp + 25);

    auto ta_xyzz_z_0 = pbuffer.data(idx_npot_0_gp + 26);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pc_x,        \
                             pc_y,        \
                             ta_xyzz_x_0, \
                             ta_xyzz_y_0, \
                             ta_xyzz_z_0, \
                             ta_xzz_x_0,  \
                             ta_xzz_x_1,  \
                             ta_yzz_y_0,  \
                             ta_yzz_y_1,  \
                             ta_yzz_z_0,  \
                             ta_yzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyzz_x_0[i] = ta_xzz_x_0[i] * pa_y[i] - ta_xzz_x_1[i] * pc_y[i];

        ta_xyzz_y_0[i] = ta_yzz_y_0[i] * pa_x[i] - ta_yzz_y_1[i] * pc_x[i];

        ta_xyzz_z_0[i] = ta_yzz_z_0[i] * pa_x[i] - ta_yzz_z_1[i] * pc_x[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto ta_xzzz_x_0 = pbuffer.data(idx_npot_0_gp + 27);

    auto ta_xzzz_y_0 = pbuffer.data(idx_npot_0_gp + 28);

    auto ta_xzzz_z_0 = pbuffer.data(idx_npot_0_gp + 29);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_xzzz_x_0, \
                             ta_xzzz_y_0, \
                             ta_xzzz_z_0, \
                             ta_zzz_0_0,  \
                             ta_zzz_0_1,  \
                             ta_zzz_x_0,  \
                             ta_zzz_x_1,  \
                             ta_zzz_y_0,  \
                             ta_zzz_y_1,  \
                             ta_zzz_z_0,  \
                             ta_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzz_x_0[i] =
            ta_zzz_0_0[i] * fe_0 - ta_zzz_0_1[i] * fe_0 + ta_zzz_x_0[i] * pa_x[i] - ta_zzz_x_1[i] * pc_x[i];

        ta_xzzz_y_0[i] = ta_zzz_y_0[i] * pa_x[i] - ta_zzz_y_1[i] * pc_x[i];

        ta_xzzz_z_0[i] = ta_zzz_z_0[i] * pa_x[i] - ta_zzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto ta_yyyy_x_0 = pbuffer.data(idx_npot_0_gp + 30);

    auto ta_yyyy_y_0 = pbuffer.data(idx_npot_0_gp + 31);

    auto ta_yyyy_z_0 = pbuffer.data(idx_npot_0_gp + 32);

#pragma omp simd aligned(pa_y,            \
                             pc_y,        \
                             ta_yy_x_0,   \
                             ta_yy_x_1,   \
                             ta_yy_y_0,   \
                             ta_yy_y_1,   \
                             ta_yy_z_0,   \
                             ta_yy_z_1,   \
                             ta_yyy_0_0,  \
                             ta_yyy_0_1,  \
                             ta_yyy_x_0,  \
                             ta_yyy_x_1,  \
                             ta_yyy_y_0,  \
                             ta_yyy_y_1,  \
                             ta_yyy_z_0,  \
                             ta_yyy_z_1,  \
                             ta_yyyy_x_0, \
                             ta_yyyy_y_0, \
                             ta_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyy_x_0[i] =
            3.0 * ta_yy_x_0[i] * fe_0 - 3.0 * ta_yy_x_1[i] * fe_0 + ta_yyy_x_0[i] * pa_y[i] - ta_yyy_x_1[i] * pc_y[i];

        ta_yyyy_y_0[i] = 3.0 * ta_yy_y_0[i] * fe_0 - 3.0 * ta_yy_y_1[i] * fe_0 + ta_yyy_0_0[i] * fe_0 -
                         ta_yyy_0_1[i] * fe_0 + ta_yyy_y_0[i] * pa_y[i] - ta_yyy_y_1[i] * pc_y[i];

        ta_yyyy_z_0[i] =
            3.0 * ta_yy_z_0[i] * fe_0 - 3.0 * ta_yy_z_1[i] * fe_0 + ta_yyy_z_0[i] * pa_y[i] - ta_yyy_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto ta_yyyz_x_0 = pbuffer.data(idx_npot_0_gp + 33);

    auto ta_yyyz_y_0 = pbuffer.data(idx_npot_0_gp + 34);

    auto ta_yyyz_z_0 = pbuffer.data(idx_npot_0_gp + 35);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             pc_y,        \
                             pc_z,        \
                             ta_yyy_x_0,  \
                             ta_yyy_x_1,  \
                             ta_yyy_y_0,  \
                             ta_yyy_y_1,  \
                             ta_yyyz_x_0, \
                             ta_yyyz_y_0, \
                             ta_yyyz_z_0, \
                             ta_yyz_z_0,  \
                             ta_yyz_z_1,  \
                             ta_yz_z_0,   \
                             ta_yz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyz_x_0[i] = ta_yyy_x_0[i] * pa_z[i] - ta_yyy_x_1[i] * pc_z[i];

        ta_yyyz_y_0[i] = ta_yyy_y_0[i] * pa_z[i] - ta_yyy_y_1[i] * pc_z[i];

        ta_yyyz_z_0[i] =
            2.0 * ta_yz_z_0[i] * fe_0 - 2.0 * ta_yz_z_1[i] * fe_0 + ta_yyz_z_0[i] * pa_y[i] - ta_yyz_z_1[i] * pc_y[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto ta_yyzz_x_0 = pbuffer.data(idx_npot_0_gp + 36);

    auto ta_yyzz_y_0 = pbuffer.data(idx_npot_0_gp + 37);

    auto ta_yyzz_z_0 = pbuffer.data(idx_npot_0_gp + 38);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             pc_y,        \
                             pc_z,        \
                             ta_yy_y_0,   \
                             ta_yy_y_1,   \
                             ta_yyz_y_0,  \
                             ta_yyz_y_1,  \
                             ta_yyzz_x_0, \
                             ta_yyzz_y_0, \
                             ta_yyzz_z_0, \
                             ta_yzz_x_0,  \
                             ta_yzz_x_1,  \
                             ta_yzz_z_0,  \
                             ta_yzz_z_1,  \
                             ta_zz_x_0,   \
                             ta_zz_x_1,   \
                             ta_zz_z_0,   \
                             ta_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzz_x_0[i] = ta_zz_x_0[i] * fe_0 - ta_zz_x_1[i] * fe_0 + ta_yzz_x_0[i] * pa_y[i] - ta_yzz_x_1[i] * pc_y[i];

        ta_yyzz_y_0[i] = ta_yy_y_0[i] * fe_0 - ta_yy_y_1[i] * fe_0 + ta_yyz_y_0[i] * pa_z[i] - ta_yyz_y_1[i] * pc_z[i];

        ta_yyzz_z_0[i] = ta_zz_z_0[i] * fe_0 - ta_zz_z_1[i] * fe_0 + ta_yzz_z_0[i] * pa_y[i] - ta_yzz_z_1[i] * pc_y[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto ta_yzzz_x_0 = pbuffer.data(idx_npot_0_gp + 39);

    auto ta_yzzz_y_0 = pbuffer.data(idx_npot_0_gp + 40);

    auto ta_yzzz_z_0 = pbuffer.data(idx_npot_0_gp + 41);

#pragma omp simd aligned(pa_y,            \
                             pc_y,        \
                             ta_yzzz_x_0, \
                             ta_yzzz_y_0, \
                             ta_yzzz_z_0, \
                             ta_zzz_0_0,  \
                             ta_zzz_0_1,  \
                             ta_zzz_x_0,  \
                             ta_zzz_x_1,  \
                             ta_zzz_y_0,  \
                             ta_zzz_y_1,  \
                             ta_zzz_z_0,  \
                             ta_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzz_x_0[i] = ta_zzz_x_0[i] * pa_y[i] - ta_zzz_x_1[i] * pc_y[i];

        ta_yzzz_y_0[i] =
            ta_zzz_0_0[i] * fe_0 - ta_zzz_0_1[i] * fe_0 + ta_zzz_y_0[i] * pa_y[i] - ta_zzz_y_1[i] * pc_y[i];

        ta_yzzz_z_0[i] = ta_zzz_z_0[i] * pa_y[i] - ta_zzz_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto ta_zzzz_x_0 = pbuffer.data(idx_npot_0_gp + 42);

    auto ta_zzzz_y_0 = pbuffer.data(idx_npot_0_gp + 43);

    auto ta_zzzz_z_0 = pbuffer.data(idx_npot_0_gp + 44);

#pragma omp simd aligned(pa_z,            \
                             pc_z,        \
                             ta_zz_x_0,   \
                             ta_zz_x_1,   \
                             ta_zz_y_0,   \
                             ta_zz_y_1,   \
                             ta_zz_z_0,   \
                             ta_zz_z_1,   \
                             ta_zzz_0_0,  \
                             ta_zzz_0_1,  \
                             ta_zzz_x_0,  \
                             ta_zzz_x_1,  \
                             ta_zzz_y_0,  \
                             ta_zzz_y_1,  \
                             ta_zzz_z_0,  \
                             ta_zzz_z_1,  \
                             ta_zzzz_x_0, \
                             ta_zzzz_y_0, \
                             ta_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzz_x_0[i] =
            3.0 * ta_zz_x_0[i] * fe_0 - 3.0 * ta_zz_x_1[i] * fe_0 + ta_zzz_x_0[i] * pa_z[i] - ta_zzz_x_1[i] * pc_z[i];

        ta_zzzz_y_0[i] =
            3.0 * ta_zz_y_0[i] * fe_0 - 3.0 * ta_zz_y_1[i] * fe_0 + ta_zzz_y_0[i] * pa_z[i] - ta_zzz_y_1[i] * pc_z[i];

        ta_zzzz_z_0[i] = 3.0 * ta_zz_z_0[i] * fe_0 - 3.0 * ta_zz_z_1[i] * fe_0 + ta_zzz_0_0[i] * fe_0 -
                         ta_zzz_0_1[i] * fe_0 + ta_zzz_z_0[i] * pa_z[i] - ta_zzz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
