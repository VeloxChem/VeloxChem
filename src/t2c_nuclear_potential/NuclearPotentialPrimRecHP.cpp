#include "NuclearPotentialPrimRecHP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hp(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hp,
                               const size_t              idx_npot_0_fp,
                               const size_t              idx_npot_1_fp,
                               const size_t              idx_npot_0_gs,
                               const size_t              idx_npot_1_gs,
                               const size_t              idx_npot_0_gp,
                               const size_t              idx_npot_1_gp,
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

    // Set up components of auxiliary buffer : FP

    auto ta_xxx_x_0 = pbuffer.data(idx_npot_0_fp);

    auto ta_xxx_y_0 = pbuffer.data(idx_npot_0_fp + 1);

    auto ta_xxx_z_0 = pbuffer.data(idx_npot_0_fp + 2);

    auto ta_xxy_x_0 = pbuffer.data(idx_npot_0_fp + 3);

    auto ta_xxy_y_0 = pbuffer.data(idx_npot_0_fp + 4);

    auto ta_xxz_x_0 = pbuffer.data(idx_npot_0_fp + 6);

    auto ta_xxz_z_0 = pbuffer.data(idx_npot_0_fp + 8);

    auto ta_xyy_y_0 = pbuffer.data(idx_npot_0_fp + 10);

    auto ta_xyy_z_0 = pbuffer.data(idx_npot_0_fp + 11);

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

    auto ta_xyy_y_1 = pbuffer.data(idx_npot_1_fp + 10);

    auto ta_xyy_z_1 = pbuffer.data(idx_npot_1_fp + 11);

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

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_0 = pbuffer.data(idx_npot_0_gs);

    auto ta_yyyy_0_0 = pbuffer.data(idx_npot_0_gs + 10);

    auto ta_yyzz_0_0 = pbuffer.data(idx_npot_0_gs + 12);

    auto ta_zzzz_0_0 = pbuffer.data(idx_npot_0_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto ta_xxxx_0_1 = pbuffer.data(idx_npot_1_gs);

    auto ta_yyyy_0_1 = pbuffer.data(idx_npot_1_gs + 10);

    auto ta_yyzz_0_1 = pbuffer.data(idx_npot_1_gs + 12);

    auto ta_zzzz_0_1 = pbuffer.data(idx_npot_1_gs + 14);

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_0 = pbuffer.data(idx_npot_0_gp);

    auto ta_xxxx_y_0 = pbuffer.data(idx_npot_0_gp + 1);

    auto ta_xxxx_z_0 = pbuffer.data(idx_npot_0_gp + 2);

    auto ta_xxxy_x_0 = pbuffer.data(idx_npot_0_gp + 3);

    auto ta_xxxy_y_0 = pbuffer.data(idx_npot_0_gp + 4);

    auto ta_xxxz_x_0 = pbuffer.data(idx_npot_0_gp + 6);

    auto ta_xxxz_z_0 = pbuffer.data(idx_npot_0_gp + 8);

    auto ta_xxyy_x_0 = pbuffer.data(idx_npot_0_gp + 9);

    auto ta_xxyy_y_0 = pbuffer.data(idx_npot_0_gp + 10);

    auto ta_xxyy_z_0 = pbuffer.data(idx_npot_0_gp + 11);

    auto ta_xxzz_x_0 = pbuffer.data(idx_npot_0_gp + 15);

    auto ta_xxzz_y_0 = pbuffer.data(idx_npot_0_gp + 16);

    auto ta_xxzz_z_0 = pbuffer.data(idx_npot_0_gp + 17);

    auto ta_xyyy_x_0 = pbuffer.data(idx_npot_0_gp + 18);

    auto ta_xyyy_y_0 = pbuffer.data(idx_npot_0_gp + 19);

    auto ta_xyyy_z_0 = pbuffer.data(idx_npot_0_gp + 20);

    auto ta_xyyz_z_0 = pbuffer.data(idx_npot_0_gp + 23);

    auto ta_xyzz_y_0 = pbuffer.data(idx_npot_0_gp + 25);

    auto ta_xzzz_x_0 = pbuffer.data(idx_npot_0_gp + 27);

    auto ta_xzzz_y_0 = pbuffer.data(idx_npot_0_gp + 28);

    auto ta_xzzz_z_0 = pbuffer.data(idx_npot_0_gp + 29);

    auto ta_yyyy_x_0 = pbuffer.data(idx_npot_0_gp + 30);

    auto ta_yyyy_y_0 = pbuffer.data(idx_npot_0_gp + 31);

    auto ta_yyyy_z_0 = pbuffer.data(idx_npot_0_gp + 32);

    auto ta_yyyz_y_0 = pbuffer.data(idx_npot_0_gp + 34);

    auto ta_yyyz_z_0 = pbuffer.data(idx_npot_0_gp + 35);

    auto ta_yyzz_x_0 = pbuffer.data(idx_npot_0_gp + 36);

    auto ta_yyzz_y_0 = pbuffer.data(idx_npot_0_gp + 37);

    auto ta_yyzz_z_0 = pbuffer.data(idx_npot_0_gp + 38);

    auto ta_yzzz_x_0 = pbuffer.data(idx_npot_0_gp + 39);

    auto ta_yzzz_y_0 = pbuffer.data(idx_npot_0_gp + 40);

    auto ta_yzzz_z_0 = pbuffer.data(idx_npot_0_gp + 41);

    auto ta_zzzz_x_0 = pbuffer.data(idx_npot_0_gp + 42);

    auto ta_zzzz_y_0 = pbuffer.data(idx_npot_0_gp + 43);

    auto ta_zzzz_z_0 = pbuffer.data(idx_npot_0_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_1 = pbuffer.data(idx_npot_1_gp);

    auto ta_xxxx_y_1 = pbuffer.data(idx_npot_1_gp + 1);

    auto ta_xxxx_z_1 = pbuffer.data(idx_npot_1_gp + 2);

    auto ta_xxxy_x_1 = pbuffer.data(idx_npot_1_gp + 3);

    auto ta_xxxy_y_1 = pbuffer.data(idx_npot_1_gp + 4);

    auto ta_xxxz_x_1 = pbuffer.data(idx_npot_1_gp + 6);

    auto ta_xxxz_z_1 = pbuffer.data(idx_npot_1_gp + 8);

    auto ta_xxyy_x_1 = pbuffer.data(idx_npot_1_gp + 9);

    auto ta_xxyy_y_1 = pbuffer.data(idx_npot_1_gp + 10);

    auto ta_xxyy_z_1 = pbuffer.data(idx_npot_1_gp + 11);

    auto ta_xxzz_x_1 = pbuffer.data(idx_npot_1_gp + 15);

    auto ta_xxzz_y_1 = pbuffer.data(idx_npot_1_gp + 16);

    auto ta_xxzz_z_1 = pbuffer.data(idx_npot_1_gp + 17);

    auto ta_xyyy_x_1 = pbuffer.data(idx_npot_1_gp + 18);

    auto ta_xyyy_y_1 = pbuffer.data(idx_npot_1_gp + 19);

    auto ta_xyyy_z_1 = pbuffer.data(idx_npot_1_gp + 20);

    auto ta_xyyz_z_1 = pbuffer.data(idx_npot_1_gp + 23);

    auto ta_xyzz_y_1 = pbuffer.data(idx_npot_1_gp + 25);

    auto ta_xzzz_x_1 = pbuffer.data(idx_npot_1_gp + 27);

    auto ta_xzzz_y_1 = pbuffer.data(idx_npot_1_gp + 28);

    auto ta_xzzz_z_1 = pbuffer.data(idx_npot_1_gp + 29);

    auto ta_yyyy_x_1 = pbuffer.data(idx_npot_1_gp + 30);

    auto ta_yyyy_y_1 = pbuffer.data(idx_npot_1_gp + 31);

    auto ta_yyyy_z_1 = pbuffer.data(idx_npot_1_gp + 32);

    auto ta_yyyz_y_1 = pbuffer.data(idx_npot_1_gp + 34);

    auto ta_yyyz_z_1 = pbuffer.data(idx_npot_1_gp + 35);

    auto ta_yyzz_x_1 = pbuffer.data(idx_npot_1_gp + 36);

    auto ta_yyzz_y_1 = pbuffer.data(idx_npot_1_gp + 37);

    auto ta_yyzz_z_1 = pbuffer.data(idx_npot_1_gp + 38);

    auto ta_yzzz_x_1 = pbuffer.data(idx_npot_1_gp + 39);

    auto ta_yzzz_y_1 = pbuffer.data(idx_npot_1_gp + 40);

    auto ta_yzzz_z_1 = pbuffer.data(idx_npot_1_gp + 41);

    auto ta_zzzz_x_1 = pbuffer.data(idx_npot_1_gp + 42);

    auto ta_zzzz_y_1 = pbuffer.data(idx_npot_1_gp + 43);

    auto ta_zzzz_z_1 = pbuffer.data(idx_npot_1_gp + 44);

    // Set up 0-3 components of targeted buffer : HP

    auto ta_xxxxx_x_0 = pbuffer.data(idx_npot_0_hp);

    auto ta_xxxxx_y_0 = pbuffer.data(idx_npot_0_hp + 1);

    auto ta_xxxxx_z_0 = pbuffer.data(idx_npot_0_hp + 2);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xxx_x_0,   \
                             ta_xxx_x_1,   \
                             ta_xxx_y_0,   \
                             ta_xxx_y_1,   \
                             ta_xxx_z_0,   \
                             ta_xxx_z_1,   \
                             ta_xxxx_0_0,  \
                             ta_xxxx_0_1,  \
                             ta_xxxx_x_0,  \
                             ta_xxxx_x_1,  \
                             ta_xxxx_y_0,  \
                             ta_xxxx_y_1,  \
                             ta_xxxx_z_0,  \
                             ta_xxxx_z_1,  \
                             ta_xxxxx_x_0, \
                             ta_xxxxx_y_0, \
                             ta_xxxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_x_0[i] = 4.0 * ta_xxx_x_0[i] * fe_0 - 4.0 * ta_xxx_x_1[i] * fe_0 + ta_xxxx_0_0[i] * fe_0 - ta_xxxx_0_1[i] * fe_0 +
                          ta_xxxx_x_0[i] * pa_x[i] - ta_xxxx_x_1[i] * pc_x[i];

        ta_xxxxx_y_0[i] = 4.0 * ta_xxx_y_0[i] * fe_0 - 4.0 * ta_xxx_y_1[i] * fe_0 + ta_xxxx_y_0[i] * pa_x[i] - ta_xxxx_y_1[i] * pc_x[i];

        ta_xxxxx_z_0[i] = 4.0 * ta_xxx_z_0[i] * fe_0 - 4.0 * ta_xxx_z_1[i] * fe_0 + ta_xxxx_z_0[i] * pa_x[i] - ta_xxxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : HP

    auto ta_xxxxy_x_0 = pbuffer.data(idx_npot_0_hp + 3);

    auto ta_xxxxy_y_0 = pbuffer.data(idx_npot_0_hp + 4);

    auto ta_xxxxy_z_0 = pbuffer.data(idx_npot_0_hp + 5);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xxxx_x_0,  \
                             ta_xxxx_x_1,  \
                             ta_xxxx_z_0,  \
                             ta_xxxx_z_1,  \
                             ta_xxxxy_x_0, \
                             ta_xxxxy_y_0, \
                             ta_xxxxy_z_0, \
                             ta_xxxy_y_0,  \
                             ta_xxxy_y_1,  \
                             ta_xxy_y_0,   \
                             ta_xxy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxy_x_0[i] = ta_xxxx_x_0[i] * pa_y[i] - ta_xxxx_x_1[i] * pc_y[i];

        ta_xxxxy_y_0[i] = 3.0 * ta_xxy_y_0[i] * fe_0 - 3.0 * ta_xxy_y_1[i] * fe_0 + ta_xxxy_y_0[i] * pa_x[i] - ta_xxxy_y_1[i] * pc_x[i];

        ta_xxxxy_z_0[i] = ta_xxxx_z_0[i] * pa_y[i] - ta_xxxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : HP

    auto ta_xxxxz_x_0 = pbuffer.data(idx_npot_0_hp + 6);

    auto ta_xxxxz_y_0 = pbuffer.data(idx_npot_0_hp + 7);

    auto ta_xxxxz_z_0 = pbuffer.data(idx_npot_0_hp + 8);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xxxx_x_0,  \
                             ta_xxxx_x_1,  \
                             ta_xxxx_y_0,  \
                             ta_xxxx_y_1,  \
                             ta_xxxxz_x_0, \
                             ta_xxxxz_y_0, \
                             ta_xxxxz_z_0, \
                             ta_xxxz_z_0,  \
                             ta_xxxz_z_1,  \
                             ta_xxz_z_0,   \
                             ta_xxz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxz_x_0[i] = ta_xxxx_x_0[i] * pa_z[i] - ta_xxxx_x_1[i] * pc_z[i];

        ta_xxxxz_y_0[i] = ta_xxxx_y_0[i] * pa_z[i] - ta_xxxx_y_1[i] * pc_z[i];

        ta_xxxxz_z_0[i] = 3.0 * ta_xxz_z_0[i] * fe_0 - 3.0 * ta_xxz_z_1[i] * fe_0 + ta_xxxz_z_0[i] * pa_x[i] - ta_xxxz_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : HP

    auto ta_xxxyy_x_0 = pbuffer.data(idx_npot_0_hp + 9);

    auto ta_xxxyy_y_0 = pbuffer.data(idx_npot_0_hp + 10);

    auto ta_xxxyy_z_0 = pbuffer.data(idx_npot_0_hp + 11);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xxx_x_0,   \
                             ta_xxx_x_1,   \
                             ta_xxxy_x_0,  \
                             ta_xxxy_x_1,  \
                             ta_xxxyy_x_0, \
                             ta_xxxyy_y_0, \
                             ta_xxxyy_z_0, \
                             ta_xxyy_y_0,  \
                             ta_xxyy_y_1,  \
                             ta_xxyy_z_0,  \
                             ta_xxyy_z_1,  \
                             ta_xyy_y_0,   \
                             ta_xyy_y_1,   \
                             ta_xyy_z_0,   \
                             ta_xyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyy_x_0[i] = ta_xxx_x_0[i] * fe_0 - ta_xxx_x_1[i] * fe_0 + ta_xxxy_x_0[i] * pa_y[i] - ta_xxxy_x_1[i] * pc_y[i];

        ta_xxxyy_y_0[i] = 2.0 * ta_xyy_y_0[i] * fe_0 - 2.0 * ta_xyy_y_1[i] * fe_0 + ta_xxyy_y_0[i] * pa_x[i] - ta_xxyy_y_1[i] * pc_x[i];

        ta_xxxyy_z_0[i] = 2.0 * ta_xyy_z_0[i] * fe_0 - 2.0 * ta_xyy_z_1[i] * fe_0 + ta_xxyy_z_0[i] * pa_x[i] - ta_xxyy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : HP

    auto ta_xxxyz_x_0 = pbuffer.data(idx_npot_0_hp + 12);

    auto ta_xxxyz_y_0 = pbuffer.data(idx_npot_0_hp + 13);

    auto ta_xxxyz_z_0 = pbuffer.data(idx_npot_0_hp + 14);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_xxxy_y_0,  \
                             ta_xxxy_y_1,  \
                             ta_xxxyz_x_0, \
                             ta_xxxyz_y_0, \
                             ta_xxxyz_z_0, \
                             ta_xxxz_x_0,  \
                             ta_xxxz_x_1,  \
                             ta_xxxz_z_0,  \
                             ta_xxxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xxxyz_x_0[i] = ta_xxxz_x_0[i] * pa_y[i] - ta_xxxz_x_1[i] * pc_y[i];

        ta_xxxyz_y_0[i] = ta_xxxy_y_0[i] * pa_z[i] - ta_xxxy_y_1[i] * pc_z[i];

        ta_xxxyz_z_0[i] = ta_xxxz_z_0[i] * pa_y[i] - ta_xxxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : HP

    auto ta_xxxzz_x_0 = pbuffer.data(idx_npot_0_hp + 15);

    auto ta_xxxzz_y_0 = pbuffer.data(idx_npot_0_hp + 16);

    auto ta_xxxzz_z_0 = pbuffer.data(idx_npot_0_hp + 17);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xxx_x_0,   \
                             ta_xxx_x_1,   \
                             ta_xxxz_x_0,  \
                             ta_xxxz_x_1,  \
                             ta_xxxzz_x_0, \
                             ta_xxxzz_y_0, \
                             ta_xxxzz_z_0, \
                             ta_xxzz_y_0,  \
                             ta_xxzz_y_1,  \
                             ta_xxzz_z_0,  \
                             ta_xxzz_z_1,  \
                             ta_xzz_y_0,   \
                             ta_xzz_y_1,   \
                             ta_xzz_z_0,   \
                             ta_xzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzz_x_0[i] = ta_xxx_x_0[i] * fe_0 - ta_xxx_x_1[i] * fe_0 + ta_xxxz_x_0[i] * pa_z[i] - ta_xxxz_x_1[i] * pc_z[i];

        ta_xxxzz_y_0[i] = 2.0 * ta_xzz_y_0[i] * fe_0 - 2.0 * ta_xzz_y_1[i] * fe_0 + ta_xxzz_y_0[i] * pa_x[i] - ta_xxzz_y_1[i] * pc_x[i];

        ta_xxxzz_z_0[i] = 2.0 * ta_xzz_z_0[i] * fe_0 - 2.0 * ta_xzz_z_1[i] * fe_0 + ta_xxzz_z_0[i] * pa_x[i] - ta_xxzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : HP

    auto ta_xxyyy_x_0 = pbuffer.data(idx_npot_0_hp + 18);

    auto ta_xxyyy_y_0 = pbuffer.data(idx_npot_0_hp + 19);

    auto ta_xxyyy_z_0 = pbuffer.data(idx_npot_0_hp + 20);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xxy_x_0,   \
                             ta_xxy_x_1,   \
                             ta_xxyy_x_0,  \
                             ta_xxyy_x_1,  \
                             ta_xxyyy_x_0, \
                             ta_xxyyy_y_0, \
                             ta_xxyyy_z_0, \
                             ta_xyyy_y_0,  \
                             ta_xyyy_y_1,  \
                             ta_xyyy_z_0,  \
                             ta_xyyy_z_1,  \
                             ta_yyy_y_0,   \
                             ta_yyy_y_1,   \
                             ta_yyy_z_0,   \
                             ta_yyy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyy_x_0[i] = 2.0 * ta_xxy_x_0[i] * fe_0 - 2.0 * ta_xxy_x_1[i] * fe_0 + ta_xxyy_x_0[i] * pa_y[i] - ta_xxyy_x_1[i] * pc_y[i];

        ta_xxyyy_y_0[i] = ta_yyy_y_0[i] * fe_0 - ta_yyy_y_1[i] * fe_0 + ta_xyyy_y_0[i] * pa_x[i] - ta_xyyy_y_1[i] * pc_x[i];

        ta_xxyyy_z_0[i] = ta_yyy_z_0[i] * fe_0 - ta_yyy_z_1[i] * fe_0 + ta_xyyy_z_0[i] * pa_x[i] - ta_xyyy_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : HP

    auto ta_xxyyz_x_0 = pbuffer.data(idx_npot_0_hp + 21);

    auto ta_xxyyz_y_0 = pbuffer.data(idx_npot_0_hp + 22);

    auto ta_xxyyz_z_0 = pbuffer.data(idx_npot_0_hp + 23);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xxyy_x_0,  \
                             ta_xxyy_x_1,  \
                             ta_xxyy_y_0,  \
                             ta_xxyy_y_1,  \
                             ta_xxyyz_x_0, \
                             ta_xxyyz_y_0, \
                             ta_xxyyz_z_0, \
                             ta_xyyz_z_0,  \
                             ta_xyyz_z_1,  \
                             ta_yyz_z_0,   \
                             ta_yyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyz_x_0[i] = ta_xxyy_x_0[i] * pa_z[i] - ta_xxyy_x_1[i] * pc_z[i];

        ta_xxyyz_y_0[i] = ta_xxyy_y_0[i] * pa_z[i] - ta_xxyy_y_1[i] * pc_z[i];

        ta_xxyyz_z_0[i] = ta_yyz_z_0[i] * fe_0 - ta_yyz_z_1[i] * fe_0 + ta_xyyz_z_0[i] * pa_x[i] - ta_xyyz_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : HP

    auto ta_xxyzz_x_0 = pbuffer.data(idx_npot_0_hp + 24);

    auto ta_xxyzz_y_0 = pbuffer.data(idx_npot_0_hp + 25);

    auto ta_xxyzz_z_0 = pbuffer.data(idx_npot_0_hp + 26);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xxyzz_x_0, \
                             ta_xxyzz_y_0, \
                             ta_xxyzz_z_0, \
                             ta_xxzz_x_0,  \
                             ta_xxzz_x_1,  \
                             ta_xxzz_z_0,  \
                             ta_xxzz_z_1,  \
                             ta_xyzz_y_0,  \
                             ta_xyzz_y_1,  \
                             ta_yzz_y_0,   \
                             ta_yzz_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzz_x_0[i] = ta_xxzz_x_0[i] * pa_y[i] - ta_xxzz_x_1[i] * pc_y[i];

        ta_xxyzz_y_0[i] = ta_yzz_y_0[i] * fe_0 - ta_yzz_y_1[i] * fe_0 + ta_xyzz_y_0[i] * pa_x[i] - ta_xyzz_y_1[i] * pc_x[i];

        ta_xxyzz_z_0[i] = ta_xxzz_z_0[i] * pa_y[i] - ta_xxzz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : HP

    auto ta_xxzzz_x_0 = pbuffer.data(idx_npot_0_hp + 27);

    auto ta_xxzzz_y_0 = pbuffer.data(idx_npot_0_hp + 28);

    auto ta_xxzzz_z_0 = pbuffer.data(idx_npot_0_hp + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xxz_x_0,   \
                             ta_xxz_x_1,   \
                             ta_xxzz_x_0,  \
                             ta_xxzz_x_1,  \
                             ta_xxzzz_x_0, \
                             ta_xxzzz_y_0, \
                             ta_xxzzz_z_0, \
                             ta_xzzz_y_0,  \
                             ta_xzzz_y_1,  \
                             ta_xzzz_z_0,  \
                             ta_xzzz_z_1,  \
                             ta_zzz_y_0,   \
                             ta_zzz_y_1,   \
                             ta_zzz_z_0,   \
                             ta_zzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzz_x_0[i] = 2.0 * ta_xxz_x_0[i] * fe_0 - 2.0 * ta_xxz_x_1[i] * fe_0 + ta_xxzz_x_0[i] * pa_z[i] - ta_xxzz_x_1[i] * pc_z[i];

        ta_xxzzz_y_0[i] = ta_zzz_y_0[i] * fe_0 - ta_zzz_y_1[i] * fe_0 + ta_xzzz_y_0[i] * pa_x[i] - ta_xzzz_y_1[i] * pc_x[i];

        ta_xxzzz_z_0[i] = ta_zzz_z_0[i] * fe_0 - ta_zzz_z_1[i] * fe_0 + ta_xzzz_z_0[i] * pa_x[i] - ta_xzzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : HP

    auto ta_xyyyy_x_0 = pbuffer.data(idx_npot_0_hp + 30);

    auto ta_xyyyy_y_0 = pbuffer.data(idx_npot_0_hp + 31);

    auto ta_xyyyy_z_0 = pbuffer.data(idx_npot_0_hp + 32);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xyyyy_x_0, \
                             ta_xyyyy_y_0, \
                             ta_xyyyy_z_0, \
                             ta_yyyy_0_0,  \
                             ta_yyyy_0_1,  \
                             ta_yyyy_x_0,  \
                             ta_yyyy_x_1,  \
                             ta_yyyy_y_0,  \
                             ta_yyyy_y_1,  \
                             ta_yyyy_z_0,  \
                             ta_yyyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyy_x_0[i] = ta_yyyy_0_0[i] * fe_0 - ta_yyyy_0_1[i] * fe_0 + ta_yyyy_x_0[i] * pa_x[i] - ta_yyyy_x_1[i] * pc_x[i];

        ta_xyyyy_y_0[i] = ta_yyyy_y_0[i] * pa_x[i] - ta_yyyy_y_1[i] * pc_x[i];

        ta_xyyyy_z_0[i] = ta_yyyy_z_0[i] * pa_x[i] - ta_yyyy_z_1[i] * pc_x[i];
    }

    // Set up 33-36 components of targeted buffer : HP

    auto ta_xyyyz_x_0 = pbuffer.data(idx_npot_0_hp + 33);

    auto ta_xyyyz_y_0 = pbuffer.data(idx_npot_0_hp + 34);

    auto ta_xyyyz_z_0 = pbuffer.data(idx_npot_0_hp + 35);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_xyyy_x_0,  \
                             ta_xyyy_x_1,  \
                             ta_xyyyz_x_0, \
                             ta_xyyyz_y_0, \
                             ta_xyyyz_z_0, \
                             ta_yyyz_y_0,  \
                             ta_yyyz_y_1,  \
                             ta_yyyz_z_0,  \
                             ta_yyyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyyyz_x_0[i] = ta_xyyy_x_0[i] * pa_z[i] - ta_xyyy_x_1[i] * pc_z[i];

        ta_xyyyz_y_0[i] = ta_yyyz_y_0[i] * pa_x[i] - ta_yyyz_y_1[i] * pc_x[i];

        ta_xyyyz_z_0[i] = ta_yyyz_z_0[i] * pa_x[i] - ta_yyyz_z_1[i] * pc_x[i];
    }

    // Set up 36-39 components of targeted buffer : HP

    auto ta_xyyzz_x_0 = pbuffer.data(idx_npot_0_hp + 36);

    auto ta_xyyzz_y_0 = pbuffer.data(idx_npot_0_hp + 37);

    auto ta_xyyzz_z_0 = pbuffer.data(idx_npot_0_hp + 38);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xyyzz_x_0, \
                             ta_xyyzz_y_0, \
                             ta_xyyzz_z_0, \
                             ta_yyzz_0_0,  \
                             ta_yyzz_0_1,  \
                             ta_yyzz_x_0,  \
                             ta_yyzz_x_1,  \
                             ta_yyzz_y_0,  \
                             ta_yyzz_y_1,  \
                             ta_yyzz_z_0,  \
                             ta_yyzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzz_x_0[i] = ta_yyzz_0_0[i] * fe_0 - ta_yyzz_0_1[i] * fe_0 + ta_yyzz_x_0[i] * pa_x[i] - ta_yyzz_x_1[i] * pc_x[i];

        ta_xyyzz_y_0[i] = ta_yyzz_y_0[i] * pa_x[i] - ta_yyzz_y_1[i] * pc_x[i];

        ta_xyyzz_z_0[i] = ta_yyzz_z_0[i] * pa_x[i] - ta_yyzz_z_1[i] * pc_x[i];
    }

    // Set up 39-42 components of targeted buffer : HP

    auto ta_xyzzz_x_0 = pbuffer.data(idx_npot_0_hp + 39);

    auto ta_xyzzz_y_0 = pbuffer.data(idx_npot_0_hp + 40);

    auto ta_xyzzz_z_0 = pbuffer.data(idx_npot_0_hp + 41);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_xyzzz_x_0, \
                             ta_xyzzz_y_0, \
                             ta_xyzzz_z_0, \
                             ta_xzzz_x_0,  \
                             ta_xzzz_x_1,  \
                             ta_yzzz_y_0,  \
                             ta_yzzz_y_1,  \
                             ta_yzzz_z_0,  \
                             ta_yzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyzzz_x_0[i] = ta_xzzz_x_0[i] * pa_y[i] - ta_xzzz_x_1[i] * pc_y[i];

        ta_xyzzz_y_0[i] = ta_yzzz_y_0[i] * pa_x[i] - ta_yzzz_y_1[i] * pc_x[i];

        ta_xyzzz_z_0[i] = ta_yzzz_z_0[i] * pa_x[i] - ta_yzzz_z_1[i] * pc_x[i];
    }

    // Set up 42-45 components of targeted buffer : HP

    auto ta_xzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 42);

    auto ta_xzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 43);

    auto ta_xzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 44);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_xzzzz_x_0, \
                             ta_xzzzz_y_0, \
                             ta_xzzzz_z_0, \
                             ta_zzzz_0_0,  \
                             ta_zzzz_0_1,  \
                             ta_zzzz_x_0,  \
                             ta_zzzz_x_1,  \
                             ta_zzzz_y_0,  \
                             ta_zzzz_y_1,  \
                             ta_zzzz_z_0,  \
                             ta_zzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzz_x_0[i] = ta_zzzz_0_0[i] * fe_0 - ta_zzzz_0_1[i] * fe_0 + ta_zzzz_x_0[i] * pa_x[i] - ta_zzzz_x_1[i] * pc_x[i];

        ta_xzzzz_y_0[i] = ta_zzzz_y_0[i] * pa_x[i] - ta_zzzz_y_1[i] * pc_x[i];

        ta_xzzzz_z_0[i] = ta_zzzz_z_0[i] * pa_x[i] - ta_zzzz_z_1[i] * pc_x[i];
    }

    // Set up 45-48 components of targeted buffer : HP

    auto ta_yyyyy_x_0 = pbuffer.data(idx_npot_0_hp + 45);

    auto ta_yyyyy_y_0 = pbuffer.data(idx_npot_0_hp + 46);

    auto ta_yyyyy_z_0 = pbuffer.data(idx_npot_0_hp + 47);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_yyy_x_0,   \
                             ta_yyy_x_1,   \
                             ta_yyy_y_0,   \
                             ta_yyy_y_1,   \
                             ta_yyy_z_0,   \
                             ta_yyy_z_1,   \
                             ta_yyyy_0_0,  \
                             ta_yyyy_0_1,  \
                             ta_yyyy_x_0,  \
                             ta_yyyy_x_1,  \
                             ta_yyyy_y_0,  \
                             ta_yyyy_y_1,  \
                             ta_yyyy_z_0,  \
                             ta_yyyy_z_1,  \
                             ta_yyyyy_x_0, \
                             ta_yyyyy_y_0, \
                             ta_yyyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyy_x_0[i] = 4.0 * ta_yyy_x_0[i] * fe_0 - 4.0 * ta_yyy_x_1[i] * fe_0 + ta_yyyy_x_0[i] * pa_y[i] - ta_yyyy_x_1[i] * pc_y[i];

        ta_yyyyy_y_0[i] = 4.0 * ta_yyy_y_0[i] * fe_0 - 4.0 * ta_yyy_y_1[i] * fe_0 + ta_yyyy_0_0[i] * fe_0 - ta_yyyy_0_1[i] * fe_0 +
                          ta_yyyy_y_0[i] * pa_y[i] - ta_yyyy_y_1[i] * pc_y[i];

        ta_yyyyy_z_0[i] = 4.0 * ta_yyy_z_0[i] * fe_0 - 4.0 * ta_yyy_z_1[i] * fe_0 + ta_yyyy_z_0[i] * pa_y[i] - ta_yyyy_z_1[i] * pc_y[i];
    }

    // Set up 48-51 components of targeted buffer : HP

    auto ta_yyyyz_x_0 = pbuffer.data(idx_npot_0_hp + 48);

    auto ta_yyyyz_y_0 = pbuffer.data(idx_npot_0_hp + 49);

    auto ta_yyyyz_z_0 = pbuffer.data(idx_npot_0_hp + 50);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_yyyy_x_0,  \
                             ta_yyyy_x_1,  \
                             ta_yyyy_y_0,  \
                             ta_yyyy_y_1,  \
                             ta_yyyyz_x_0, \
                             ta_yyyyz_y_0, \
                             ta_yyyyz_z_0, \
                             ta_yyyz_z_0,  \
                             ta_yyyz_z_1,  \
                             ta_yyz_z_0,   \
                             ta_yyz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyz_x_0[i] = ta_yyyy_x_0[i] * pa_z[i] - ta_yyyy_x_1[i] * pc_z[i];

        ta_yyyyz_y_0[i] = ta_yyyy_y_0[i] * pa_z[i] - ta_yyyy_y_1[i] * pc_z[i];

        ta_yyyyz_z_0[i] = 3.0 * ta_yyz_z_0[i] * fe_0 - 3.0 * ta_yyz_z_1[i] * fe_0 + ta_yyyz_z_0[i] * pa_y[i] - ta_yyyz_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : HP

    auto ta_yyyzz_x_0 = pbuffer.data(idx_npot_0_hp + 51);

    auto ta_yyyzz_y_0 = pbuffer.data(idx_npot_0_hp + 52);

    auto ta_yyyzz_z_0 = pbuffer.data(idx_npot_0_hp + 53);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_yyy_y_0,   \
                             ta_yyy_y_1,   \
                             ta_yyyz_y_0,  \
                             ta_yyyz_y_1,  \
                             ta_yyyzz_x_0, \
                             ta_yyyzz_y_0, \
                             ta_yyyzz_z_0, \
                             ta_yyzz_x_0,  \
                             ta_yyzz_x_1,  \
                             ta_yyzz_z_0,  \
                             ta_yyzz_z_1,  \
                             ta_yzz_x_0,   \
                             ta_yzz_x_1,   \
                             ta_yzz_z_0,   \
                             ta_yzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzz_x_0[i] = 2.0 * ta_yzz_x_0[i] * fe_0 - 2.0 * ta_yzz_x_1[i] * fe_0 + ta_yyzz_x_0[i] * pa_y[i] - ta_yyzz_x_1[i] * pc_y[i];

        ta_yyyzz_y_0[i] = ta_yyy_y_0[i] * fe_0 - ta_yyy_y_1[i] * fe_0 + ta_yyyz_y_0[i] * pa_z[i] - ta_yyyz_y_1[i] * pc_z[i];

        ta_yyyzz_z_0[i] = 2.0 * ta_yzz_z_0[i] * fe_0 - 2.0 * ta_yzz_z_1[i] * fe_0 + ta_yyzz_z_0[i] * pa_y[i] - ta_yyzz_z_1[i] * pc_y[i];
    }

    // Set up 54-57 components of targeted buffer : HP

    auto ta_yyzzz_x_0 = pbuffer.data(idx_npot_0_hp + 54);

    auto ta_yyzzz_y_0 = pbuffer.data(idx_npot_0_hp + 55);

    auto ta_yyzzz_z_0 = pbuffer.data(idx_npot_0_hp + 56);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_yyz_y_0,   \
                             ta_yyz_y_1,   \
                             ta_yyzz_y_0,  \
                             ta_yyzz_y_1,  \
                             ta_yyzzz_x_0, \
                             ta_yyzzz_y_0, \
                             ta_yyzzz_z_0, \
                             ta_yzzz_x_0,  \
                             ta_yzzz_x_1,  \
                             ta_yzzz_z_0,  \
                             ta_yzzz_z_1,  \
                             ta_zzz_x_0,   \
                             ta_zzz_x_1,   \
                             ta_zzz_z_0,   \
                             ta_zzz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzz_x_0[i] = ta_zzz_x_0[i] * fe_0 - ta_zzz_x_1[i] * fe_0 + ta_yzzz_x_0[i] * pa_y[i] - ta_yzzz_x_1[i] * pc_y[i];

        ta_yyzzz_y_0[i] = 2.0 * ta_yyz_y_0[i] * fe_0 - 2.0 * ta_yyz_y_1[i] * fe_0 + ta_yyzz_y_0[i] * pa_z[i] - ta_yyzz_y_1[i] * pc_z[i];

        ta_yyzzz_z_0[i] = ta_zzz_z_0[i] * fe_0 - ta_zzz_z_1[i] * fe_0 + ta_yzzz_z_0[i] * pa_y[i] - ta_yzzz_z_1[i] * pc_y[i];
    }

    // Set up 57-60 components of targeted buffer : HP

    auto ta_yzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 57);

    auto ta_yzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 58);

    auto ta_yzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 59);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_yzzzz_x_0, \
                             ta_yzzzz_y_0, \
                             ta_yzzzz_z_0, \
                             ta_zzzz_0_0,  \
                             ta_zzzz_0_1,  \
                             ta_zzzz_x_0,  \
                             ta_zzzz_x_1,  \
                             ta_zzzz_y_0,  \
                             ta_zzzz_y_1,  \
                             ta_zzzz_z_0,  \
                             ta_zzzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzz_x_0[i] = ta_zzzz_x_0[i] * pa_y[i] - ta_zzzz_x_1[i] * pc_y[i];

        ta_yzzzz_y_0[i] = ta_zzzz_0_0[i] * fe_0 - ta_zzzz_0_1[i] * fe_0 + ta_zzzz_y_0[i] * pa_y[i] - ta_zzzz_y_1[i] * pc_y[i];

        ta_yzzzz_z_0[i] = ta_zzzz_z_0[i] * pa_y[i] - ta_zzzz_z_1[i] * pc_y[i];
    }

    // Set up 60-63 components of targeted buffer : HP

    auto ta_zzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 60);

    auto ta_zzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 61);

    auto ta_zzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 62);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta_zzz_x_0,   \
                             ta_zzz_x_1,   \
                             ta_zzz_y_0,   \
                             ta_zzz_y_1,   \
                             ta_zzz_z_0,   \
                             ta_zzz_z_1,   \
                             ta_zzzz_0_0,  \
                             ta_zzzz_0_1,  \
                             ta_zzzz_x_0,  \
                             ta_zzzz_x_1,  \
                             ta_zzzz_y_0,  \
                             ta_zzzz_y_1,  \
                             ta_zzzz_z_0,  \
                             ta_zzzz_z_1,  \
                             ta_zzzzz_x_0, \
                             ta_zzzzz_y_0, \
                             ta_zzzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzz_x_0[i] = 4.0 * ta_zzz_x_0[i] * fe_0 - 4.0 * ta_zzz_x_1[i] * fe_0 + ta_zzzz_x_0[i] * pa_z[i] - ta_zzzz_x_1[i] * pc_z[i];

        ta_zzzzz_y_0[i] = 4.0 * ta_zzz_y_0[i] * fe_0 - 4.0 * ta_zzz_y_1[i] * fe_0 + ta_zzzz_y_0[i] * pa_z[i] - ta_zzzz_y_1[i] * pc_z[i];

        ta_zzzzz_z_0[i] = 4.0 * ta_zzz_z_0[i] * fe_0 - 4.0 * ta_zzz_z_1[i] * fe_0 + ta_zzzz_0_0[i] * fe_0 - ta_zzzz_0_1[i] * fe_0 +
                          ta_zzzz_z_0[i] * pa_z[i] - ta_zzzz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
