#include "NuclearPotentialPrimRecPG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_pg(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_pg,
                               const size_t              idx_npot_0_sf,
                               const size_t              idx_npot_1_sf,
                               const size_t              idx_npot_0_sg,
                               const size_t              idx_npot_1_sg,
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

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_0 = pbuffer.data(idx_npot_0_sg);

    auto ta_0_xxxy_0 = pbuffer.data(idx_npot_0_sg + 1);

    auto ta_0_xxxz_0 = pbuffer.data(idx_npot_0_sg + 2);

    auto ta_0_xxyy_0 = pbuffer.data(idx_npot_0_sg + 3);

    auto ta_0_xxyz_0 = pbuffer.data(idx_npot_0_sg + 4);

    auto ta_0_xxzz_0 = pbuffer.data(idx_npot_0_sg + 5);

    auto ta_0_xyyy_0 = pbuffer.data(idx_npot_0_sg + 6);

    auto ta_0_xyyz_0 = pbuffer.data(idx_npot_0_sg + 7);

    auto ta_0_xyzz_0 = pbuffer.data(idx_npot_0_sg + 8);

    auto ta_0_xzzz_0 = pbuffer.data(idx_npot_0_sg + 9);

    auto ta_0_yyyy_0 = pbuffer.data(idx_npot_0_sg + 10);

    auto ta_0_yyyz_0 = pbuffer.data(idx_npot_0_sg + 11);

    auto ta_0_yyzz_0 = pbuffer.data(idx_npot_0_sg + 12);

    auto ta_0_yzzz_0 = pbuffer.data(idx_npot_0_sg + 13);

    auto ta_0_zzzz_0 = pbuffer.data(idx_npot_0_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto ta_0_xxxx_1 = pbuffer.data(idx_npot_1_sg);

    auto ta_0_xxxy_1 = pbuffer.data(idx_npot_1_sg + 1);

    auto ta_0_xxxz_1 = pbuffer.data(idx_npot_1_sg + 2);

    auto ta_0_xxyy_1 = pbuffer.data(idx_npot_1_sg + 3);

    auto ta_0_xxyz_1 = pbuffer.data(idx_npot_1_sg + 4);

    auto ta_0_xxzz_1 = pbuffer.data(idx_npot_1_sg + 5);

    auto ta_0_xyyy_1 = pbuffer.data(idx_npot_1_sg + 6);

    auto ta_0_xyyz_1 = pbuffer.data(idx_npot_1_sg + 7);

    auto ta_0_xyzz_1 = pbuffer.data(idx_npot_1_sg + 8);

    auto ta_0_xzzz_1 = pbuffer.data(idx_npot_1_sg + 9);

    auto ta_0_yyyy_1 = pbuffer.data(idx_npot_1_sg + 10);

    auto ta_0_yyyz_1 = pbuffer.data(idx_npot_1_sg + 11);

    auto ta_0_yyzz_1 = pbuffer.data(idx_npot_1_sg + 12);

    auto ta_0_yzzz_1 = pbuffer.data(idx_npot_1_sg + 13);

    auto ta_0_zzzz_1 = pbuffer.data(idx_npot_1_sg + 14);

    // Set up 0-15 components of targeted buffer : PG

    auto ta_x_xxxx_0 = pbuffer.data(idx_npot_0_pg);

    auto ta_x_xxxy_0 = pbuffer.data(idx_npot_0_pg + 1);

    auto ta_x_xxxz_0 = pbuffer.data(idx_npot_0_pg + 2);

    auto ta_x_xxyy_0 = pbuffer.data(idx_npot_0_pg + 3);

    auto ta_x_xxyz_0 = pbuffer.data(idx_npot_0_pg + 4);

    auto ta_x_xxzz_0 = pbuffer.data(idx_npot_0_pg + 5);

    auto ta_x_xyyy_0 = pbuffer.data(idx_npot_0_pg + 6);

    auto ta_x_xyyz_0 = pbuffer.data(idx_npot_0_pg + 7);

    auto ta_x_xyzz_0 = pbuffer.data(idx_npot_0_pg + 8);

    auto ta_x_xzzz_0 = pbuffer.data(idx_npot_0_pg + 9);

    auto ta_x_yyyy_0 = pbuffer.data(idx_npot_0_pg + 10);

    auto ta_x_yyyz_0 = pbuffer.data(idx_npot_0_pg + 11);

    auto ta_x_yyzz_0 = pbuffer.data(idx_npot_0_pg + 12);

    auto ta_x_yzzz_0 = pbuffer.data(idx_npot_0_pg + 13);

    auto ta_x_zzzz_0 = pbuffer.data(idx_npot_0_pg + 14);

#pragma omp simd aligned(pa_x,            \
                             pc_x,        \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxxx_0, \
                             ta_0_xxxx_1, \
                             ta_0_xxxy_0, \
                             ta_0_xxxy_1, \
                             ta_0_xxxz_0, \
                             ta_0_xxxz_1, \
                             ta_0_xxy_0,  \
                             ta_0_xxy_1,  \
                             ta_0_xxyy_0, \
                             ta_0_xxyy_1, \
                             ta_0_xxyz_0, \
                             ta_0_xxyz_1, \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xxzz_0, \
                             ta_0_xxzz_1, \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyyy_0, \
                             ta_0_xyyy_1, \
                             ta_0_xyyz_0, \
                             ta_0_xyyz_1, \
                             ta_0_xyz_0,  \
                             ta_0_xyz_1,  \
                             ta_0_xyzz_0, \
                             ta_0_xyzz_1, \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_xzzz_0, \
                             ta_0_xzzz_1, \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyyy_0, \
                             ta_0_yyyy_1, \
                             ta_0_yyyz_0, \
                             ta_0_yyyz_1, \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yyzz_0, \
                             ta_0_yyzz_1, \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_yzzz_0, \
                             ta_0_yzzz_1, \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_0_zzzz_0, \
                             ta_0_zzzz_1, \
                             ta_x_xxxx_0, \
                             ta_x_xxxy_0, \
                             ta_x_xxxz_0, \
                             ta_x_xxyy_0, \
                             ta_x_xxyz_0, \
                             ta_x_xxzz_0, \
                             ta_x_xyyy_0, \
                             ta_x_xyyz_0, \
                             ta_x_xyzz_0, \
                             ta_x_xzzz_0, \
                             ta_x_yyyy_0, \
                             ta_x_yyyz_0, \
                             ta_x_yyzz_0, \
                             ta_x_yzzz_0, \
                             ta_x_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_x_xxxx_0[i] = 4.0 * ta_0_xxx_0[i] * fe_0 - 4.0 * ta_0_xxx_1[i] * fe_0 + ta_0_xxxx_0[i] * pa_x[i] - ta_0_xxxx_1[i] * pc_x[i];

        ta_x_xxxy_0[i] = 3.0 * ta_0_xxy_0[i] * fe_0 - 3.0 * ta_0_xxy_1[i] * fe_0 + ta_0_xxxy_0[i] * pa_x[i] - ta_0_xxxy_1[i] * pc_x[i];

        ta_x_xxxz_0[i] = 3.0 * ta_0_xxz_0[i] * fe_0 - 3.0 * ta_0_xxz_1[i] * fe_0 + ta_0_xxxz_0[i] * pa_x[i] - ta_0_xxxz_1[i] * pc_x[i];

        ta_x_xxyy_0[i] = 2.0 * ta_0_xyy_0[i] * fe_0 - 2.0 * ta_0_xyy_1[i] * fe_0 + ta_0_xxyy_0[i] * pa_x[i] - ta_0_xxyy_1[i] * pc_x[i];

        ta_x_xxyz_0[i] = 2.0 * ta_0_xyz_0[i] * fe_0 - 2.0 * ta_0_xyz_1[i] * fe_0 + ta_0_xxyz_0[i] * pa_x[i] - ta_0_xxyz_1[i] * pc_x[i];

        ta_x_xxzz_0[i] = 2.0 * ta_0_xzz_0[i] * fe_0 - 2.0 * ta_0_xzz_1[i] * fe_0 + ta_0_xxzz_0[i] * pa_x[i] - ta_0_xxzz_1[i] * pc_x[i];

        ta_x_xyyy_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_0_xyyy_0[i] * pa_x[i] - ta_0_xyyy_1[i] * pc_x[i];

        ta_x_xyyz_0[i] = ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + ta_0_xyyz_0[i] * pa_x[i] - ta_0_xyyz_1[i] * pc_x[i];

        ta_x_xyzz_0[i] = ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + ta_0_xyzz_0[i] * pa_x[i] - ta_0_xyzz_1[i] * pc_x[i];

        ta_x_xzzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_0_xzzz_0[i] * pa_x[i] - ta_0_xzzz_1[i] * pc_x[i];

        ta_x_yyyy_0[i] = ta_0_yyyy_0[i] * pa_x[i] - ta_0_yyyy_1[i] * pc_x[i];

        ta_x_yyyz_0[i] = ta_0_yyyz_0[i] * pa_x[i] - ta_0_yyyz_1[i] * pc_x[i];

        ta_x_yyzz_0[i] = ta_0_yyzz_0[i] * pa_x[i] - ta_0_yyzz_1[i] * pc_x[i];

        ta_x_yzzz_0[i] = ta_0_yzzz_0[i] * pa_x[i] - ta_0_yzzz_1[i] * pc_x[i];

        ta_x_zzzz_0[i] = ta_0_zzzz_0[i] * pa_x[i] - ta_0_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

    auto ta_y_xxxx_0 = pbuffer.data(idx_npot_0_pg + 15);

    auto ta_y_xxxy_0 = pbuffer.data(idx_npot_0_pg + 16);

    auto ta_y_xxxz_0 = pbuffer.data(idx_npot_0_pg + 17);

    auto ta_y_xxyy_0 = pbuffer.data(idx_npot_0_pg + 18);

    auto ta_y_xxyz_0 = pbuffer.data(idx_npot_0_pg + 19);

    auto ta_y_xxzz_0 = pbuffer.data(idx_npot_0_pg + 20);

    auto ta_y_xyyy_0 = pbuffer.data(idx_npot_0_pg + 21);

    auto ta_y_xyyz_0 = pbuffer.data(idx_npot_0_pg + 22);

    auto ta_y_xyzz_0 = pbuffer.data(idx_npot_0_pg + 23);

    auto ta_y_xzzz_0 = pbuffer.data(idx_npot_0_pg + 24);

    auto ta_y_yyyy_0 = pbuffer.data(idx_npot_0_pg + 25);

    auto ta_y_yyyz_0 = pbuffer.data(idx_npot_0_pg + 26);

    auto ta_y_yyzz_0 = pbuffer.data(idx_npot_0_pg + 27);

    auto ta_y_yzzz_0 = pbuffer.data(idx_npot_0_pg + 28);

    auto ta_y_zzzz_0 = pbuffer.data(idx_npot_0_pg + 29);

#pragma omp simd aligned(pa_y,            \
                             pc_y,        \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxxx_0, \
                             ta_0_xxxx_1, \
                             ta_0_xxxy_0, \
                             ta_0_xxxy_1, \
                             ta_0_xxxz_0, \
                             ta_0_xxxz_1, \
                             ta_0_xxy_0,  \
                             ta_0_xxy_1,  \
                             ta_0_xxyy_0, \
                             ta_0_xxyy_1, \
                             ta_0_xxyz_0, \
                             ta_0_xxyz_1, \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xxzz_0, \
                             ta_0_xxzz_1, \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyyy_0, \
                             ta_0_xyyy_1, \
                             ta_0_xyyz_0, \
                             ta_0_xyyz_1, \
                             ta_0_xyz_0,  \
                             ta_0_xyz_1,  \
                             ta_0_xyzz_0, \
                             ta_0_xyzz_1, \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_xzzz_0, \
                             ta_0_xzzz_1, \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyyy_0, \
                             ta_0_yyyy_1, \
                             ta_0_yyyz_0, \
                             ta_0_yyyz_1, \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yyzz_0, \
                             ta_0_yyzz_1, \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_yzzz_0, \
                             ta_0_yzzz_1, \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_0_zzzz_0, \
                             ta_0_zzzz_1, \
                             ta_y_xxxx_0, \
                             ta_y_xxxy_0, \
                             ta_y_xxxz_0, \
                             ta_y_xxyy_0, \
                             ta_y_xxyz_0, \
                             ta_y_xxzz_0, \
                             ta_y_xyyy_0, \
                             ta_y_xyyz_0, \
                             ta_y_xyzz_0, \
                             ta_y_xzzz_0, \
                             ta_y_yyyy_0, \
                             ta_y_yyyz_0, \
                             ta_y_yyzz_0, \
                             ta_y_yzzz_0, \
                             ta_y_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_y_xxxx_0[i] = ta_0_xxxx_0[i] * pa_y[i] - ta_0_xxxx_1[i] * pc_y[i];

        ta_y_xxxy_0[i] = ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + ta_0_xxxy_0[i] * pa_y[i] - ta_0_xxxy_1[i] * pc_y[i];

        ta_y_xxxz_0[i] = ta_0_xxxz_0[i] * pa_y[i] - ta_0_xxxz_1[i] * pc_y[i];

        ta_y_xxyy_0[i] = 2.0 * ta_0_xxy_0[i] * fe_0 - 2.0 * ta_0_xxy_1[i] * fe_0 + ta_0_xxyy_0[i] * pa_y[i] - ta_0_xxyy_1[i] * pc_y[i];

        ta_y_xxyz_0[i] = ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + ta_0_xxyz_0[i] * pa_y[i] - ta_0_xxyz_1[i] * pc_y[i];

        ta_y_xxzz_0[i] = ta_0_xxzz_0[i] * pa_y[i] - ta_0_xxzz_1[i] * pc_y[i];

        ta_y_xyyy_0[i] = 3.0 * ta_0_xyy_0[i] * fe_0 - 3.0 * ta_0_xyy_1[i] * fe_0 + ta_0_xyyy_0[i] * pa_y[i] - ta_0_xyyy_1[i] * pc_y[i];

        ta_y_xyyz_0[i] = 2.0 * ta_0_xyz_0[i] * fe_0 - 2.0 * ta_0_xyz_1[i] * fe_0 + ta_0_xyyz_0[i] * pa_y[i] - ta_0_xyyz_1[i] * pc_y[i];

        ta_y_xyzz_0[i] = ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + ta_0_xyzz_0[i] * pa_y[i] - ta_0_xyzz_1[i] * pc_y[i];

        ta_y_xzzz_0[i] = ta_0_xzzz_0[i] * pa_y[i] - ta_0_xzzz_1[i] * pc_y[i];

        ta_y_yyyy_0[i] = 4.0 * ta_0_yyy_0[i] * fe_0 - 4.0 * ta_0_yyy_1[i] * fe_0 + ta_0_yyyy_0[i] * pa_y[i] - ta_0_yyyy_1[i] * pc_y[i];

        ta_y_yyyz_0[i] = 3.0 * ta_0_yyz_0[i] * fe_0 - 3.0 * ta_0_yyz_1[i] * fe_0 + ta_0_yyyz_0[i] * pa_y[i] - ta_0_yyyz_1[i] * pc_y[i];

        ta_y_yyzz_0[i] = 2.0 * ta_0_yzz_0[i] * fe_0 - 2.0 * ta_0_yzz_1[i] * fe_0 + ta_0_yyzz_0[i] * pa_y[i] - ta_0_yyzz_1[i] * pc_y[i];

        ta_y_yzzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_0_yzzz_0[i] * pa_y[i] - ta_0_yzzz_1[i] * pc_y[i];

        ta_y_zzzz_0[i] = ta_0_zzzz_0[i] * pa_y[i] - ta_0_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : PG

    auto ta_z_xxxx_0 = pbuffer.data(idx_npot_0_pg + 30);

    auto ta_z_xxxy_0 = pbuffer.data(idx_npot_0_pg + 31);

    auto ta_z_xxxz_0 = pbuffer.data(idx_npot_0_pg + 32);

    auto ta_z_xxyy_0 = pbuffer.data(idx_npot_0_pg + 33);

    auto ta_z_xxyz_0 = pbuffer.data(idx_npot_0_pg + 34);

    auto ta_z_xxzz_0 = pbuffer.data(idx_npot_0_pg + 35);

    auto ta_z_xyyy_0 = pbuffer.data(idx_npot_0_pg + 36);

    auto ta_z_xyyz_0 = pbuffer.data(idx_npot_0_pg + 37);

    auto ta_z_xyzz_0 = pbuffer.data(idx_npot_0_pg + 38);

    auto ta_z_xzzz_0 = pbuffer.data(idx_npot_0_pg + 39);

    auto ta_z_yyyy_0 = pbuffer.data(idx_npot_0_pg + 40);

    auto ta_z_yyyz_0 = pbuffer.data(idx_npot_0_pg + 41);

    auto ta_z_yyzz_0 = pbuffer.data(idx_npot_0_pg + 42);

    auto ta_z_yzzz_0 = pbuffer.data(idx_npot_0_pg + 43);

    auto ta_z_zzzz_0 = pbuffer.data(idx_npot_0_pg + 44);

#pragma omp simd aligned(pa_z,            \
                             pc_z,        \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxxx_0, \
                             ta_0_xxxx_1, \
                             ta_0_xxxy_0, \
                             ta_0_xxxy_1, \
                             ta_0_xxxz_0, \
                             ta_0_xxxz_1, \
                             ta_0_xxy_0,  \
                             ta_0_xxy_1,  \
                             ta_0_xxyy_0, \
                             ta_0_xxyy_1, \
                             ta_0_xxyz_0, \
                             ta_0_xxyz_1, \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xxzz_0, \
                             ta_0_xxzz_1, \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyyy_0, \
                             ta_0_xyyy_1, \
                             ta_0_xyyz_0, \
                             ta_0_xyyz_1, \
                             ta_0_xyz_0,  \
                             ta_0_xyz_1,  \
                             ta_0_xyzz_0, \
                             ta_0_xyzz_1, \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_xzzz_0, \
                             ta_0_xzzz_1, \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyyy_0, \
                             ta_0_yyyy_1, \
                             ta_0_yyyz_0, \
                             ta_0_yyyz_1, \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yyzz_0, \
                             ta_0_yyzz_1, \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_yzzz_0, \
                             ta_0_yzzz_1, \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_0_zzzz_0, \
                             ta_0_zzzz_1, \
                             ta_z_xxxx_0, \
                             ta_z_xxxy_0, \
                             ta_z_xxxz_0, \
                             ta_z_xxyy_0, \
                             ta_z_xxyz_0, \
                             ta_z_xxzz_0, \
                             ta_z_xyyy_0, \
                             ta_z_xyyz_0, \
                             ta_z_xyzz_0, \
                             ta_z_xzzz_0, \
                             ta_z_yyyy_0, \
                             ta_z_yyyz_0, \
                             ta_z_yyzz_0, \
                             ta_z_yzzz_0, \
                             ta_z_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_z_xxxx_0[i] = ta_0_xxxx_0[i] * pa_z[i] - ta_0_xxxx_1[i] * pc_z[i];

        ta_z_xxxy_0[i] = ta_0_xxxy_0[i] * pa_z[i] - ta_0_xxxy_1[i] * pc_z[i];

        ta_z_xxxz_0[i] = ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + ta_0_xxxz_0[i] * pa_z[i] - ta_0_xxxz_1[i] * pc_z[i];

        ta_z_xxyy_0[i] = ta_0_xxyy_0[i] * pa_z[i] - ta_0_xxyy_1[i] * pc_z[i];

        ta_z_xxyz_0[i] = ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + ta_0_xxyz_0[i] * pa_z[i] - ta_0_xxyz_1[i] * pc_z[i];

        ta_z_xxzz_0[i] = 2.0 * ta_0_xxz_0[i] * fe_0 - 2.0 * ta_0_xxz_1[i] * fe_0 + ta_0_xxzz_0[i] * pa_z[i] - ta_0_xxzz_1[i] * pc_z[i];

        ta_z_xyyy_0[i] = ta_0_xyyy_0[i] * pa_z[i] - ta_0_xyyy_1[i] * pc_z[i];

        ta_z_xyyz_0[i] = ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + ta_0_xyyz_0[i] * pa_z[i] - ta_0_xyyz_1[i] * pc_z[i];

        ta_z_xyzz_0[i] = 2.0 * ta_0_xyz_0[i] * fe_0 - 2.0 * ta_0_xyz_1[i] * fe_0 + ta_0_xyzz_0[i] * pa_z[i] - ta_0_xyzz_1[i] * pc_z[i];

        ta_z_xzzz_0[i] = 3.0 * ta_0_xzz_0[i] * fe_0 - 3.0 * ta_0_xzz_1[i] * fe_0 + ta_0_xzzz_0[i] * pa_z[i] - ta_0_xzzz_1[i] * pc_z[i];

        ta_z_yyyy_0[i] = ta_0_yyyy_0[i] * pa_z[i] - ta_0_yyyy_1[i] * pc_z[i];

        ta_z_yyyz_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_0_yyyz_0[i] * pa_z[i] - ta_0_yyyz_1[i] * pc_z[i];

        ta_z_yyzz_0[i] = 2.0 * ta_0_yyz_0[i] * fe_0 - 2.0 * ta_0_yyz_1[i] * fe_0 + ta_0_yyzz_0[i] * pa_z[i] - ta_0_yyzz_1[i] * pc_z[i];

        ta_z_yzzz_0[i] = 3.0 * ta_0_yzz_0[i] * fe_0 - 3.0 * ta_0_yzz_1[i] * fe_0 + ta_0_yzzz_0[i] * pa_z[i] - ta_0_yzzz_1[i] * pc_z[i];

        ta_z_zzzz_0[i] = 4.0 * ta_0_zzz_0[i] * fe_0 - 4.0 * ta_0_zzz_1[i] * fe_0 + ta_0_zzzz_0[i] * pa_z[i] - ta_0_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
