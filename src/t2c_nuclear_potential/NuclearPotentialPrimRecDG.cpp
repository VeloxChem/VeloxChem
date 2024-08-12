#include "NuclearPotentialPrimRecDG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_dg(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_dg,
                               const size_t              idx_npot_0_sg,
                               const size_t              idx_npot_1_sg,
                               const size_t              idx_npot_0_pf,
                               const size_t              idx_npot_1_pf,
                               const size_t              idx_npot_0_pg,
                               const size_t              idx_npot_1_pg,
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

    // Set up components of auxiliary buffer : PG

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

    // Set up components of auxiliary buffer : PG

    auto ta_x_xxxx_1 = pbuffer.data(idx_npot_1_pg);

    auto ta_x_xxxy_1 = pbuffer.data(idx_npot_1_pg + 1);

    auto ta_x_xxxz_1 = pbuffer.data(idx_npot_1_pg + 2);

    auto ta_x_xxyy_1 = pbuffer.data(idx_npot_1_pg + 3);

    auto ta_x_xxyz_1 = pbuffer.data(idx_npot_1_pg + 4);

    auto ta_x_xxzz_1 = pbuffer.data(idx_npot_1_pg + 5);

    auto ta_x_xyyy_1 = pbuffer.data(idx_npot_1_pg + 6);

    auto ta_x_xyyz_1 = pbuffer.data(idx_npot_1_pg + 7);

    auto ta_x_xyzz_1 = pbuffer.data(idx_npot_1_pg + 8);

    auto ta_x_xzzz_1 = pbuffer.data(idx_npot_1_pg + 9);

    auto ta_x_yyyy_1 = pbuffer.data(idx_npot_1_pg + 10);

    auto ta_x_yyyz_1 = pbuffer.data(idx_npot_1_pg + 11);

    auto ta_x_yyzz_1 = pbuffer.data(idx_npot_1_pg + 12);

    auto ta_x_yzzz_1 = pbuffer.data(idx_npot_1_pg + 13);

    auto ta_x_zzzz_1 = pbuffer.data(idx_npot_1_pg + 14);

    auto ta_y_xxxx_1 = pbuffer.data(idx_npot_1_pg + 15);

    auto ta_y_xxxy_1 = pbuffer.data(idx_npot_1_pg + 16);

    auto ta_y_xxxz_1 = pbuffer.data(idx_npot_1_pg + 17);

    auto ta_y_xxyy_1 = pbuffer.data(idx_npot_1_pg + 18);

    auto ta_y_xxyz_1 = pbuffer.data(idx_npot_1_pg + 19);

    auto ta_y_xxzz_1 = pbuffer.data(idx_npot_1_pg + 20);

    auto ta_y_xyyy_1 = pbuffer.data(idx_npot_1_pg + 21);

    auto ta_y_xyyz_1 = pbuffer.data(idx_npot_1_pg + 22);

    auto ta_y_xyzz_1 = pbuffer.data(idx_npot_1_pg + 23);

    auto ta_y_xzzz_1 = pbuffer.data(idx_npot_1_pg + 24);

    auto ta_y_yyyy_1 = pbuffer.data(idx_npot_1_pg + 25);

    auto ta_y_yyyz_1 = pbuffer.data(idx_npot_1_pg + 26);

    auto ta_y_yyzz_1 = pbuffer.data(idx_npot_1_pg + 27);

    auto ta_y_yzzz_1 = pbuffer.data(idx_npot_1_pg + 28);

    auto ta_y_zzzz_1 = pbuffer.data(idx_npot_1_pg + 29);

    auto ta_z_xxxx_1 = pbuffer.data(idx_npot_1_pg + 30);

    auto ta_z_xxxy_1 = pbuffer.data(idx_npot_1_pg + 31);

    auto ta_z_xxxz_1 = pbuffer.data(idx_npot_1_pg + 32);

    auto ta_z_xxyy_1 = pbuffer.data(idx_npot_1_pg + 33);

    auto ta_z_xxyz_1 = pbuffer.data(idx_npot_1_pg + 34);

    auto ta_z_xxzz_1 = pbuffer.data(idx_npot_1_pg + 35);

    auto ta_z_xyyy_1 = pbuffer.data(idx_npot_1_pg + 36);

    auto ta_z_xyyz_1 = pbuffer.data(idx_npot_1_pg + 37);

    auto ta_z_xyzz_1 = pbuffer.data(idx_npot_1_pg + 38);

    auto ta_z_xzzz_1 = pbuffer.data(idx_npot_1_pg + 39);

    auto ta_z_yyyy_1 = pbuffer.data(idx_npot_1_pg + 40);

    auto ta_z_yyyz_1 = pbuffer.data(idx_npot_1_pg + 41);

    auto ta_z_yyzz_1 = pbuffer.data(idx_npot_1_pg + 42);

    auto ta_z_yzzz_1 = pbuffer.data(idx_npot_1_pg + 43);

    auto ta_z_zzzz_1 = pbuffer.data(idx_npot_1_pg + 44);

    // Set up 0-15 components of targeted buffer : DG

    auto ta_xx_xxxx_0 = pbuffer.data(idx_npot_0_dg);

    auto ta_xx_xxxy_0 = pbuffer.data(idx_npot_0_dg + 1);

    auto ta_xx_xxxz_0 = pbuffer.data(idx_npot_0_dg + 2);

    auto ta_xx_xxyy_0 = pbuffer.data(idx_npot_0_dg + 3);

    auto ta_xx_xxyz_0 = pbuffer.data(idx_npot_0_dg + 4);

    auto ta_xx_xxzz_0 = pbuffer.data(idx_npot_0_dg + 5);

    auto ta_xx_xyyy_0 = pbuffer.data(idx_npot_0_dg + 6);

    auto ta_xx_xyyz_0 = pbuffer.data(idx_npot_0_dg + 7);

    auto ta_xx_xyzz_0 = pbuffer.data(idx_npot_0_dg + 8);

    auto ta_xx_xzzz_0 = pbuffer.data(idx_npot_0_dg + 9);

    auto ta_xx_yyyy_0 = pbuffer.data(idx_npot_0_dg + 10);

    auto ta_xx_yyyz_0 = pbuffer.data(idx_npot_0_dg + 11);

    auto ta_xx_yyzz_0 = pbuffer.data(idx_npot_0_dg + 12);

    auto ta_xx_yzzz_0 = pbuffer.data(idx_npot_0_dg + 13);

    auto ta_xx_zzzz_0 = pbuffer.data(idx_npot_0_dg + 14);

#pragma omp simd aligned(pa_x,             \
                             pc_x,         \
                             ta_0_xxxx_0,  \
                             ta_0_xxxx_1,  \
                             ta_0_xxxy_0,  \
                             ta_0_xxxy_1,  \
                             ta_0_xxxz_0,  \
                             ta_0_xxxz_1,  \
                             ta_0_xxyy_0,  \
                             ta_0_xxyy_1,  \
                             ta_0_xxyz_0,  \
                             ta_0_xxyz_1,  \
                             ta_0_xxzz_0,  \
                             ta_0_xxzz_1,  \
                             ta_0_xyyy_0,  \
                             ta_0_xyyy_1,  \
                             ta_0_xyyz_0,  \
                             ta_0_xyyz_1,  \
                             ta_0_xyzz_0,  \
                             ta_0_xyzz_1,  \
                             ta_0_xzzz_0,  \
                             ta_0_xzzz_1,  \
                             ta_0_yyyy_0,  \
                             ta_0_yyyy_1,  \
                             ta_0_yyyz_0,  \
                             ta_0_yyyz_1,  \
                             ta_0_yyzz_0,  \
                             ta_0_yyzz_1,  \
                             ta_0_yzzz_0,  \
                             ta_0_yzzz_1,  \
                             ta_0_zzzz_0,  \
                             ta_0_zzzz_1,  \
                             ta_x_xxx_0,   \
                             ta_x_xxx_1,   \
                             ta_x_xxxx_0,  \
                             ta_x_xxxx_1,  \
                             ta_x_xxxy_0,  \
                             ta_x_xxxy_1,  \
                             ta_x_xxxz_0,  \
                             ta_x_xxxz_1,  \
                             ta_x_xxy_0,   \
                             ta_x_xxy_1,   \
                             ta_x_xxyy_0,  \
                             ta_x_xxyy_1,  \
                             ta_x_xxyz_0,  \
                             ta_x_xxyz_1,  \
                             ta_x_xxz_0,   \
                             ta_x_xxz_1,   \
                             ta_x_xxzz_0,  \
                             ta_x_xxzz_1,  \
                             ta_x_xyy_0,   \
                             ta_x_xyy_1,   \
                             ta_x_xyyy_0,  \
                             ta_x_xyyy_1,  \
                             ta_x_xyyz_0,  \
                             ta_x_xyyz_1,  \
                             ta_x_xyz_0,   \
                             ta_x_xyz_1,   \
                             ta_x_xyzz_0,  \
                             ta_x_xyzz_1,  \
                             ta_x_xzz_0,   \
                             ta_x_xzz_1,   \
                             ta_x_xzzz_0,  \
                             ta_x_xzzz_1,  \
                             ta_x_yyy_0,   \
                             ta_x_yyy_1,   \
                             ta_x_yyyy_0,  \
                             ta_x_yyyy_1,  \
                             ta_x_yyyz_0,  \
                             ta_x_yyyz_1,  \
                             ta_x_yyz_0,   \
                             ta_x_yyz_1,   \
                             ta_x_yyzz_0,  \
                             ta_x_yyzz_1,  \
                             ta_x_yzz_0,   \
                             ta_x_yzz_1,   \
                             ta_x_yzzz_0,  \
                             ta_x_yzzz_1,  \
                             ta_x_zzz_0,   \
                             ta_x_zzz_1,   \
                             ta_x_zzzz_0,  \
                             ta_x_zzzz_1,  \
                             ta_xx_xxxx_0, \
                             ta_xx_xxxy_0, \
                             ta_xx_xxxz_0, \
                             ta_xx_xxyy_0, \
                             ta_xx_xxyz_0, \
                             ta_xx_xxzz_0, \
                             ta_xx_xyyy_0, \
                             ta_xx_xyyz_0, \
                             ta_xx_xyzz_0, \
                             ta_xx_xzzz_0, \
                             ta_xx_yyyy_0, \
                             ta_xx_yyyz_0, \
                             ta_xx_yyzz_0, \
                             ta_xx_yzzz_0, \
                             ta_xx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xx_xxxx_0[i] = ta_0_xxxx_0[i] * fe_0 - ta_0_xxxx_1[i] * fe_0 + 4.0 * ta_x_xxx_0[i] * fe_0 - 4.0 * ta_x_xxx_1[i] * fe_0 +
                          ta_x_xxxx_0[i] * pa_x[i] - ta_x_xxxx_1[i] * pc_x[i];

        ta_xx_xxxy_0[i] = ta_0_xxxy_0[i] * fe_0 - ta_0_xxxy_1[i] * fe_0 + 3.0 * ta_x_xxy_0[i] * fe_0 - 3.0 * ta_x_xxy_1[i] * fe_0 +
                          ta_x_xxxy_0[i] * pa_x[i] - ta_x_xxxy_1[i] * pc_x[i];

        ta_xx_xxxz_0[i] = ta_0_xxxz_0[i] * fe_0 - ta_0_xxxz_1[i] * fe_0 + 3.0 * ta_x_xxz_0[i] * fe_0 - 3.0 * ta_x_xxz_1[i] * fe_0 +
                          ta_x_xxxz_0[i] * pa_x[i] - ta_x_xxxz_1[i] * pc_x[i];

        ta_xx_xxyy_0[i] = ta_0_xxyy_0[i] * fe_0 - ta_0_xxyy_1[i] * fe_0 + 2.0 * ta_x_xyy_0[i] * fe_0 - 2.0 * ta_x_xyy_1[i] * fe_0 +
                          ta_x_xxyy_0[i] * pa_x[i] - ta_x_xxyy_1[i] * pc_x[i];

        ta_xx_xxyz_0[i] = ta_0_xxyz_0[i] * fe_0 - ta_0_xxyz_1[i] * fe_0 + 2.0 * ta_x_xyz_0[i] * fe_0 - 2.0 * ta_x_xyz_1[i] * fe_0 +
                          ta_x_xxyz_0[i] * pa_x[i] - ta_x_xxyz_1[i] * pc_x[i];

        ta_xx_xxzz_0[i] = ta_0_xxzz_0[i] * fe_0 - ta_0_xxzz_1[i] * fe_0 + 2.0 * ta_x_xzz_0[i] * fe_0 - 2.0 * ta_x_xzz_1[i] * fe_0 +
                          ta_x_xxzz_0[i] * pa_x[i] - ta_x_xxzz_1[i] * pc_x[i];

        ta_xx_xyyy_0[i] = ta_0_xyyy_0[i] * fe_0 - ta_0_xyyy_1[i] * fe_0 + ta_x_yyy_0[i] * fe_0 - ta_x_yyy_1[i] * fe_0 + ta_x_xyyy_0[i] * pa_x[i] -
                          ta_x_xyyy_1[i] * pc_x[i];

        ta_xx_xyyz_0[i] = ta_0_xyyz_0[i] * fe_0 - ta_0_xyyz_1[i] * fe_0 + ta_x_yyz_0[i] * fe_0 - ta_x_yyz_1[i] * fe_0 + ta_x_xyyz_0[i] * pa_x[i] -
                          ta_x_xyyz_1[i] * pc_x[i];

        ta_xx_xyzz_0[i] = ta_0_xyzz_0[i] * fe_0 - ta_0_xyzz_1[i] * fe_0 + ta_x_yzz_0[i] * fe_0 - ta_x_yzz_1[i] * fe_0 + ta_x_xyzz_0[i] * pa_x[i] -
                          ta_x_xyzz_1[i] * pc_x[i];

        ta_xx_xzzz_0[i] = ta_0_xzzz_0[i] * fe_0 - ta_0_xzzz_1[i] * fe_0 + ta_x_zzz_0[i] * fe_0 - ta_x_zzz_1[i] * fe_0 + ta_x_xzzz_0[i] * pa_x[i] -
                          ta_x_xzzz_1[i] * pc_x[i];

        ta_xx_yyyy_0[i] = ta_0_yyyy_0[i] * fe_0 - ta_0_yyyy_1[i] * fe_0 + ta_x_yyyy_0[i] * pa_x[i] - ta_x_yyyy_1[i] * pc_x[i];

        ta_xx_yyyz_0[i] = ta_0_yyyz_0[i] * fe_0 - ta_0_yyyz_1[i] * fe_0 + ta_x_yyyz_0[i] * pa_x[i] - ta_x_yyyz_1[i] * pc_x[i];

        ta_xx_yyzz_0[i] = ta_0_yyzz_0[i] * fe_0 - ta_0_yyzz_1[i] * fe_0 + ta_x_yyzz_0[i] * pa_x[i] - ta_x_yyzz_1[i] * pc_x[i];

        ta_xx_yzzz_0[i] = ta_0_yzzz_0[i] * fe_0 - ta_0_yzzz_1[i] * fe_0 + ta_x_yzzz_0[i] * pa_x[i] - ta_x_yzzz_1[i] * pc_x[i];

        ta_xx_zzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + ta_x_zzzz_0[i] * pa_x[i] - ta_x_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto ta_xy_xxxx_0 = pbuffer.data(idx_npot_0_dg + 15);

    auto ta_xy_xxxy_0 = pbuffer.data(idx_npot_0_dg + 16);

    auto ta_xy_xxxz_0 = pbuffer.data(idx_npot_0_dg + 17);

    auto ta_xy_xxyy_0 = pbuffer.data(idx_npot_0_dg + 18);

    auto ta_xy_xxyz_0 = pbuffer.data(idx_npot_0_dg + 19);

    auto ta_xy_xxzz_0 = pbuffer.data(idx_npot_0_dg + 20);

    auto ta_xy_xyyy_0 = pbuffer.data(idx_npot_0_dg + 21);

    auto ta_xy_xyyz_0 = pbuffer.data(idx_npot_0_dg + 22);

    auto ta_xy_xyzz_0 = pbuffer.data(idx_npot_0_dg + 23);

    auto ta_xy_xzzz_0 = pbuffer.data(idx_npot_0_dg + 24);

    auto ta_xy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 25);

    auto ta_xy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 26);

    auto ta_xy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 27);

    auto ta_xy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 28);

    auto ta_xy_zzzz_0 = pbuffer.data(idx_npot_0_dg + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pc_x,         \
                             pc_y,         \
                             ta_x_xxxx_0,  \
                             ta_x_xxxx_1,  \
                             ta_x_xxxz_0,  \
                             ta_x_xxxz_1,  \
                             ta_x_xxzz_0,  \
                             ta_x_xxzz_1,  \
                             ta_x_xzzz_0,  \
                             ta_x_xzzz_1,  \
                             ta_xy_xxxx_0, \
                             ta_xy_xxxy_0, \
                             ta_xy_xxxz_0, \
                             ta_xy_xxyy_0, \
                             ta_xy_xxyz_0, \
                             ta_xy_xxzz_0, \
                             ta_xy_xyyy_0, \
                             ta_xy_xyyz_0, \
                             ta_xy_xyzz_0, \
                             ta_xy_xzzz_0, \
                             ta_xy_yyyy_0, \
                             ta_xy_yyyz_0, \
                             ta_xy_yyzz_0, \
                             ta_xy_yzzz_0, \
                             ta_xy_zzzz_0, \
                             ta_y_xxxy_0,  \
                             ta_y_xxxy_1,  \
                             ta_y_xxy_0,   \
                             ta_y_xxy_1,   \
                             ta_y_xxyy_0,  \
                             ta_y_xxyy_1,  \
                             ta_y_xxyz_0,  \
                             ta_y_xxyz_1,  \
                             ta_y_xyy_0,   \
                             ta_y_xyy_1,   \
                             ta_y_xyyy_0,  \
                             ta_y_xyyy_1,  \
                             ta_y_xyyz_0,  \
                             ta_y_xyyz_1,  \
                             ta_y_xyz_0,   \
                             ta_y_xyz_1,   \
                             ta_y_xyzz_0,  \
                             ta_y_xyzz_1,  \
                             ta_y_yyy_0,   \
                             ta_y_yyy_1,   \
                             ta_y_yyyy_0,  \
                             ta_y_yyyy_1,  \
                             ta_y_yyyz_0,  \
                             ta_y_yyyz_1,  \
                             ta_y_yyz_0,   \
                             ta_y_yyz_1,   \
                             ta_y_yyzz_0,  \
                             ta_y_yyzz_1,  \
                             ta_y_yzz_0,   \
                             ta_y_yzz_1,   \
                             ta_y_yzzz_0,  \
                             ta_y_yzzz_1,  \
                             ta_y_zzzz_0,  \
                             ta_y_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xy_xxxx_0[i] = ta_x_xxxx_0[i] * pa_y[i] - ta_x_xxxx_1[i] * pc_y[i];

        ta_xy_xxxy_0[i] = 3.0 * ta_y_xxy_0[i] * fe_0 - 3.0 * ta_y_xxy_1[i] * fe_0 + ta_y_xxxy_0[i] * pa_x[i] - ta_y_xxxy_1[i] * pc_x[i];

        ta_xy_xxxz_0[i] = ta_x_xxxz_0[i] * pa_y[i] - ta_x_xxxz_1[i] * pc_y[i];

        ta_xy_xxyy_0[i] = 2.0 * ta_y_xyy_0[i] * fe_0 - 2.0 * ta_y_xyy_1[i] * fe_0 + ta_y_xxyy_0[i] * pa_x[i] - ta_y_xxyy_1[i] * pc_x[i];

        ta_xy_xxyz_0[i] = 2.0 * ta_y_xyz_0[i] * fe_0 - 2.0 * ta_y_xyz_1[i] * fe_0 + ta_y_xxyz_0[i] * pa_x[i] - ta_y_xxyz_1[i] * pc_x[i];

        ta_xy_xxzz_0[i] = ta_x_xxzz_0[i] * pa_y[i] - ta_x_xxzz_1[i] * pc_y[i];

        ta_xy_xyyy_0[i] = ta_y_yyy_0[i] * fe_0 - ta_y_yyy_1[i] * fe_0 + ta_y_xyyy_0[i] * pa_x[i] - ta_y_xyyy_1[i] * pc_x[i];

        ta_xy_xyyz_0[i] = ta_y_yyz_0[i] * fe_0 - ta_y_yyz_1[i] * fe_0 + ta_y_xyyz_0[i] * pa_x[i] - ta_y_xyyz_1[i] * pc_x[i];

        ta_xy_xyzz_0[i] = ta_y_yzz_0[i] * fe_0 - ta_y_yzz_1[i] * fe_0 + ta_y_xyzz_0[i] * pa_x[i] - ta_y_xyzz_1[i] * pc_x[i];

        ta_xy_xzzz_0[i] = ta_x_xzzz_0[i] * pa_y[i] - ta_x_xzzz_1[i] * pc_y[i];

        ta_xy_yyyy_0[i] = ta_y_yyyy_0[i] * pa_x[i] - ta_y_yyyy_1[i] * pc_x[i];

        ta_xy_yyyz_0[i] = ta_y_yyyz_0[i] * pa_x[i] - ta_y_yyyz_1[i] * pc_x[i];

        ta_xy_yyzz_0[i] = ta_y_yyzz_0[i] * pa_x[i] - ta_y_yyzz_1[i] * pc_x[i];

        ta_xy_yzzz_0[i] = ta_y_yzzz_0[i] * pa_x[i] - ta_y_yzzz_1[i] * pc_x[i];

        ta_xy_zzzz_0[i] = ta_y_zzzz_0[i] * pa_x[i] - ta_y_zzzz_1[i] * pc_x[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto ta_xz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 30);

    auto ta_xz_xxxy_0 = pbuffer.data(idx_npot_0_dg + 31);

    auto ta_xz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 32);

    auto ta_xz_xxyy_0 = pbuffer.data(idx_npot_0_dg + 33);

    auto ta_xz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 34);

    auto ta_xz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 35);

    auto ta_xz_xyyy_0 = pbuffer.data(idx_npot_0_dg + 36);

    auto ta_xz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 37);

    auto ta_xz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 38);

    auto ta_xz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 39);

    auto ta_xz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 40);

    auto ta_xz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 41);

    auto ta_xz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 42);

    auto ta_xz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 43);

    auto ta_xz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 44);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             pc_x,         \
                             pc_z,         \
                             ta_x_xxxx_0,  \
                             ta_x_xxxx_1,  \
                             ta_x_xxxy_0,  \
                             ta_x_xxxy_1,  \
                             ta_x_xxyy_0,  \
                             ta_x_xxyy_1,  \
                             ta_x_xyyy_0,  \
                             ta_x_xyyy_1,  \
                             ta_xz_xxxx_0, \
                             ta_xz_xxxy_0, \
                             ta_xz_xxxz_0, \
                             ta_xz_xxyy_0, \
                             ta_xz_xxyz_0, \
                             ta_xz_xxzz_0, \
                             ta_xz_xyyy_0, \
                             ta_xz_xyyz_0, \
                             ta_xz_xyzz_0, \
                             ta_xz_xzzz_0, \
                             ta_xz_yyyy_0, \
                             ta_xz_yyyz_0, \
                             ta_xz_yyzz_0, \
                             ta_xz_yzzz_0, \
                             ta_xz_zzzz_0, \
                             ta_z_xxxz_0,  \
                             ta_z_xxxz_1,  \
                             ta_z_xxyz_0,  \
                             ta_z_xxyz_1,  \
                             ta_z_xxz_0,   \
                             ta_z_xxz_1,   \
                             ta_z_xxzz_0,  \
                             ta_z_xxzz_1,  \
                             ta_z_xyyz_0,  \
                             ta_z_xyyz_1,  \
                             ta_z_xyz_0,   \
                             ta_z_xyz_1,   \
                             ta_z_xyzz_0,  \
                             ta_z_xyzz_1,  \
                             ta_z_xzz_0,   \
                             ta_z_xzz_1,   \
                             ta_z_xzzz_0,  \
                             ta_z_xzzz_1,  \
                             ta_z_yyyy_0,  \
                             ta_z_yyyy_1,  \
                             ta_z_yyyz_0,  \
                             ta_z_yyyz_1,  \
                             ta_z_yyz_0,   \
                             ta_z_yyz_1,   \
                             ta_z_yyzz_0,  \
                             ta_z_yyzz_1,  \
                             ta_z_yzz_0,   \
                             ta_z_yzz_1,   \
                             ta_z_yzzz_0,  \
                             ta_z_yzzz_1,  \
                             ta_z_zzz_0,   \
                             ta_z_zzz_1,   \
                             ta_z_zzzz_0,  \
                             ta_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xz_xxxx_0[i] = ta_x_xxxx_0[i] * pa_z[i] - ta_x_xxxx_1[i] * pc_z[i];

        ta_xz_xxxy_0[i] = ta_x_xxxy_0[i] * pa_z[i] - ta_x_xxxy_1[i] * pc_z[i];

        ta_xz_xxxz_0[i] = 3.0 * ta_z_xxz_0[i] * fe_0 - 3.0 * ta_z_xxz_1[i] * fe_0 + ta_z_xxxz_0[i] * pa_x[i] - ta_z_xxxz_1[i] * pc_x[i];

        ta_xz_xxyy_0[i] = ta_x_xxyy_0[i] * pa_z[i] - ta_x_xxyy_1[i] * pc_z[i];

        ta_xz_xxyz_0[i] = 2.0 * ta_z_xyz_0[i] * fe_0 - 2.0 * ta_z_xyz_1[i] * fe_0 + ta_z_xxyz_0[i] * pa_x[i] - ta_z_xxyz_1[i] * pc_x[i];

        ta_xz_xxzz_0[i] = 2.0 * ta_z_xzz_0[i] * fe_0 - 2.0 * ta_z_xzz_1[i] * fe_0 + ta_z_xxzz_0[i] * pa_x[i] - ta_z_xxzz_1[i] * pc_x[i];

        ta_xz_xyyy_0[i] = ta_x_xyyy_0[i] * pa_z[i] - ta_x_xyyy_1[i] * pc_z[i];

        ta_xz_xyyz_0[i] = ta_z_yyz_0[i] * fe_0 - ta_z_yyz_1[i] * fe_0 + ta_z_xyyz_0[i] * pa_x[i] - ta_z_xyyz_1[i] * pc_x[i];

        ta_xz_xyzz_0[i] = ta_z_yzz_0[i] * fe_0 - ta_z_yzz_1[i] * fe_0 + ta_z_xyzz_0[i] * pa_x[i] - ta_z_xyzz_1[i] * pc_x[i];

        ta_xz_xzzz_0[i] = ta_z_zzz_0[i] * fe_0 - ta_z_zzz_1[i] * fe_0 + ta_z_xzzz_0[i] * pa_x[i] - ta_z_xzzz_1[i] * pc_x[i];

        ta_xz_yyyy_0[i] = ta_z_yyyy_0[i] * pa_x[i] - ta_z_yyyy_1[i] * pc_x[i];

        ta_xz_yyyz_0[i] = ta_z_yyyz_0[i] * pa_x[i] - ta_z_yyyz_1[i] * pc_x[i];

        ta_xz_yyzz_0[i] = ta_z_yyzz_0[i] * pa_x[i] - ta_z_yyzz_1[i] * pc_x[i];

        ta_xz_yzzz_0[i] = ta_z_yzzz_0[i] * pa_x[i] - ta_z_yzzz_1[i] * pc_x[i];

        ta_xz_zzzz_0[i] = ta_z_zzzz_0[i] * pa_x[i] - ta_z_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto ta_yy_xxxx_0 = pbuffer.data(idx_npot_0_dg + 45);

    auto ta_yy_xxxy_0 = pbuffer.data(idx_npot_0_dg + 46);

    auto ta_yy_xxxz_0 = pbuffer.data(idx_npot_0_dg + 47);

    auto ta_yy_xxyy_0 = pbuffer.data(idx_npot_0_dg + 48);

    auto ta_yy_xxyz_0 = pbuffer.data(idx_npot_0_dg + 49);

    auto ta_yy_xxzz_0 = pbuffer.data(idx_npot_0_dg + 50);

    auto ta_yy_xyyy_0 = pbuffer.data(idx_npot_0_dg + 51);

    auto ta_yy_xyyz_0 = pbuffer.data(idx_npot_0_dg + 52);

    auto ta_yy_xyzz_0 = pbuffer.data(idx_npot_0_dg + 53);

    auto ta_yy_xzzz_0 = pbuffer.data(idx_npot_0_dg + 54);

    auto ta_yy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 55);

    auto ta_yy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 56);

    auto ta_yy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 57);

    auto ta_yy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 58);

    auto ta_yy_zzzz_0 = pbuffer.data(idx_npot_0_dg + 59);

#pragma omp simd aligned(pa_y,             \
                             pc_y,         \
                             ta_0_xxxx_0,  \
                             ta_0_xxxx_1,  \
                             ta_0_xxxy_0,  \
                             ta_0_xxxy_1,  \
                             ta_0_xxxz_0,  \
                             ta_0_xxxz_1,  \
                             ta_0_xxyy_0,  \
                             ta_0_xxyy_1,  \
                             ta_0_xxyz_0,  \
                             ta_0_xxyz_1,  \
                             ta_0_xxzz_0,  \
                             ta_0_xxzz_1,  \
                             ta_0_xyyy_0,  \
                             ta_0_xyyy_1,  \
                             ta_0_xyyz_0,  \
                             ta_0_xyyz_1,  \
                             ta_0_xyzz_0,  \
                             ta_0_xyzz_1,  \
                             ta_0_xzzz_0,  \
                             ta_0_xzzz_1,  \
                             ta_0_yyyy_0,  \
                             ta_0_yyyy_1,  \
                             ta_0_yyyz_0,  \
                             ta_0_yyyz_1,  \
                             ta_0_yyzz_0,  \
                             ta_0_yyzz_1,  \
                             ta_0_yzzz_0,  \
                             ta_0_yzzz_1,  \
                             ta_0_zzzz_0,  \
                             ta_0_zzzz_1,  \
                             ta_y_xxx_0,   \
                             ta_y_xxx_1,   \
                             ta_y_xxxx_0,  \
                             ta_y_xxxx_1,  \
                             ta_y_xxxy_0,  \
                             ta_y_xxxy_1,  \
                             ta_y_xxxz_0,  \
                             ta_y_xxxz_1,  \
                             ta_y_xxy_0,   \
                             ta_y_xxy_1,   \
                             ta_y_xxyy_0,  \
                             ta_y_xxyy_1,  \
                             ta_y_xxyz_0,  \
                             ta_y_xxyz_1,  \
                             ta_y_xxz_0,   \
                             ta_y_xxz_1,   \
                             ta_y_xxzz_0,  \
                             ta_y_xxzz_1,  \
                             ta_y_xyy_0,   \
                             ta_y_xyy_1,   \
                             ta_y_xyyy_0,  \
                             ta_y_xyyy_1,  \
                             ta_y_xyyz_0,  \
                             ta_y_xyyz_1,  \
                             ta_y_xyz_0,   \
                             ta_y_xyz_1,   \
                             ta_y_xyzz_0,  \
                             ta_y_xyzz_1,  \
                             ta_y_xzz_0,   \
                             ta_y_xzz_1,   \
                             ta_y_xzzz_0,  \
                             ta_y_xzzz_1,  \
                             ta_y_yyy_0,   \
                             ta_y_yyy_1,   \
                             ta_y_yyyy_0,  \
                             ta_y_yyyy_1,  \
                             ta_y_yyyz_0,  \
                             ta_y_yyyz_1,  \
                             ta_y_yyz_0,   \
                             ta_y_yyz_1,   \
                             ta_y_yyzz_0,  \
                             ta_y_yyzz_1,  \
                             ta_y_yzz_0,   \
                             ta_y_yzz_1,   \
                             ta_y_yzzz_0,  \
                             ta_y_yzzz_1,  \
                             ta_y_zzz_0,   \
                             ta_y_zzz_1,   \
                             ta_y_zzzz_0,  \
                             ta_y_zzzz_1,  \
                             ta_yy_xxxx_0, \
                             ta_yy_xxxy_0, \
                             ta_yy_xxxz_0, \
                             ta_yy_xxyy_0, \
                             ta_yy_xxyz_0, \
                             ta_yy_xxzz_0, \
                             ta_yy_xyyy_0, \
                             ta_yy_xyyz_0, \
                             ta_yy_xyzz_0, \
                             ta_yy_xzzz_0, \
                             ta_yy_yyyy_0, \
                             ta_yy_yyyz_0, \
                             ta_yy_yyzz_0, \
                             ta_yy_yzzz_0, \
                             ta_yy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yy_xxxx_0[i] = ta_0_xxxx_0[i] * fe_0 - ta_0_xxxx_1[i] * fe_0 + ta_y_xxxx_0[i] * pa_y[i] - ta_y_xxxx_1[i] * pc_y[i];

        ta_yy_xxxy_0[i] = ta_0_xxxy_0[i] * fe_0 - ta_0_xxxy_1[i] * fe_0 + ta_y_xxx_0[i] * fe_0 - ta_y_xxx_1[i] * fe_0 + ta_y_xxxy_0[i] * pa_y[i] -
                          ta_y_xxxy_1[i] * pc_y[i];

        ta_yy_xxxz_0[i] = ta_0_xxxz_0[i] * fe_0 - ta_0_xxxz_1[i] * fe_0 + ta_y_xxxz_0[i] * pa_y[i] - ta_y_xxxz_1[i] * pc_y[i];

        ta_yy_xxyy_0[i] = ta_0_xxyy_0[i] * fe_0 - ta_0_xxyy_1[i] * fe_0 + 2.0 * ta_y_xxy_0[i] * fe_0 - 2.0 * ta_y_xxy_1[i] * fe_0 +
                          ta_y_xxyy_0[i] * pa_y[i] - ta_y_xxyy_1[i] * pc_y[i];

        ta_yy_xxyz_0[i] = ta_0_xxyz_0[i] * fe_0 - ta_0_xxyz_1[i] * fe_0 + ta_y_xxz_0[i] * fe_0 - ta_y_xxz_1[i] * fe_0 + ta_y_xxyz_0[i] * pa_y[i] -
                          ta_y_xxyz_1[i] * pc_y[i];

        ta_yy_xxzz_0[i] = ta_0_xxzz_0[i] * fe_0 - ta_0_xxzz_1[i] * fe_0 + ta_y_xxzz_0[i] * pa_y[i] - ta_y_xxzz_1[i] * pc_y[i];

        ta_yy_xyyy_0[i] = ta_0_xyyy_0[i] * fe_0 - ta_0_xyyy_1[i] * fe_0 + 3.0 * ta_y_xyy_0[i] * fe_0 - 3.0 * ta_y_xyy_1[i] * fe_0 +
                          ta_y_xyyy_0[i] * pa_y[i] - ta_y_xyyy_1[i] * pc_y[i];

        ta_yy_xyyz_0[i] = ta_0_xyyz_0[i] * fe_0 - ta_0_xyyz_1[i] * fe_0 + 2.0 * ta_y_xyz_0[i] * fe_0 - 2.0 * ta_y_xyz_1[i] * fe_0 +
                          ta_y_xyyz_0[i] * pa_y[i] - ta_y_xyyz_1[i] * pc_y[i];

        ta_yy_xyzz_0[i] = ta_0_xyzz_0[i] * fe_0 - ta_0_xyzz_1[i] * fe_0 + ta_y_xzz_0[i] * fe_0 - ta_y_xzz_1[i] * fe_0 + ta_y_xyzz_0[i] * pa_y[i] -
                          ta_y_xyzz_1[i] * pc_y[i];

        ta_yy_xzzz_0[i] = ta_0_xzzz_0[i] * fe_0 - ta_0_xzzz_1[i] * fe_0 + ta_y_xzzz_0[i] * pa_y[i] - ta_y_xzzz_1[i] * pc_y[i];

        ta_yy_yyyy_0[i] = ta_0_yyyy_0[i] * fe_0 - ta_0_yyyy_1[i] * fe_0 + 4.0 * ta_y_yyy_0[i] * fe_0 - 4.0 * ta_y_yyy_1[i] * fe_0 +
                          ta_y_yyyy_0[i] * pa_y[i] - ta_y_yyyy_1[i] * pc_y[i];

        ta_yy_yyyz_0[i] = ta_0_yyyz_0[i] * fe_0 - ta_0_yyyz_1[i] * fe_0 + 3.0 * ta_y_yyz_0[i] * fe_0 - 3.0 * ta_y_yyz_1[i] * fe_0 +
                          ta_y_yyyz_0[i] * pa_y[i] - ta_y_yyyz_1[i] * pc_y[i];

        ta_yy_yyzz_0[i] = ta_0_yyzz_0[i] * fe_0 - ta_0_yyzz_1[i] * fe_0 + 2.0 * ta_y_yzz_0[i] * fe_0 - 2.0 * ta_y_yzz_1[i] * fe_0 +
                          ta_y_yyzz_0[i] * pa_y[i] - ta_y_yyzz_1[i] * pc_y[i];

        ta_yy_yzzz_0[i] = ta_0_yzzz_0[i] * fe_0 - ta_0_yzzz_1[i] * fe_0 + ta_y_zzz_0[i] * fe_0 - ta_y_zzz_1[i] * fe_0 + ta_y_yzzz_0[i] * pa_y[i] -
                          ta_y_yzzz_1[i] * pc_y[i];

        ta_yy_zzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + ta_y_zzzz_0[i] * pa_y[i] - ta_y_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : DG

    auto ta_yz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 60);

    auto ta_yz_xxxy_0 = pbuffer.data(idx_npot_0_dg + 61);

    auto ta_yz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 62);

    auto ta_yz_xxyy_0 = pbuffer.data(idx_npot_0_dg + 63);

    auto ta_yz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 64);

    auto ta_yz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 65);

    auto ta_yz_xyyy_0 = pbuffer.data(idx_npot_0_dg + 66);

    auto ta_yz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 67);

    auto ta_yz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 68);

    auto ta_yz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 69);

    auto ta_yz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 70);

    auto ta_yz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 71);

    auto ta_yz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 72);

    auto ta_yz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 73);

    auto ta_yz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 74);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             pc_y,         \
                             pc_z,         \
                             ta_y_xxxy_0,  \
                             ta_y_xxxy_1,  \
                             ta_y_xxyy_0,  \
                             ta_y_xxyy_1,  \
                             ta_y_xyyy_0,  \
                             ta_y_xyyy_1,  \
                             ta_y_yyyy_0,  \
                             ta_y_yyyy_1,  \
                             ta_yz_xxxx_0, \
                             ta_yz_xxxy_0, \
                             ta_yz_xxxz_0, \
                             ta_yz_xxyy_0, \
                             ta_yz_xxyz_0, \
                             ta_yz_xxzz_0, \
                             ta_yz_xyyy_0, \
                             ta_yz_xyyz_0, \
                             ta_yz_xyzz_0, \
                             ta_yz_xzzz_0, \
                             ta_yz_yyyy_0, \
                             ta_yz_yyyz_0, \
                             ta_yz_yyzz_0, \
                             ta_yz_yzzz_0, \
                             ta_yz_zzzz_0, \
                             ta_z_xxxx_0,  \
                             ta_z_xxxx_1,  \
                             ta_z_xxxz_0,  \
                             ta_z_xxxz_1,  \
                             ta_z_xxyz_0,  \
                             ta_z_xxyz_1,  \
                             ta_z_xxz_0,   \
                             ta_z_xxz_1,   \
                             ta_z_xxzz_0,  \
                             ta_z_xxzz_1,  \
                             ta_z_xyyz_0,  \
                             ta_z_xyyz_1,  \
                             ta_z_xyz_0,   \
                             ta_z_xyz_1,   \
                             ta_z_xyzz_0,  \
                             ta_z_xyzz_1,  \
                             ta_z_xzz_0,   \
                             ta_z_xzz_1,   \
                             ta_z_xzzz_0,  \
                             ta_z_xzzz_1,  \
                             ta_z_yyyz_0,  \
                             ta_z_yyyz_1,  \
                             ta_z_yyz_0,   \
                             ta_z_yyz_1,   \
                             ta_z_yyzz_0,  \
                             ta_z_yyzz_1,  \
                             ta_z_yzz_0,   \
                             ta_z_yzz_1,   \
                             ta_z_yzzz_0,  \
                             ta_z_yzzz_1,  \
                             ta_z_zzz_0,   \
                             ta_z_zzz_1,   \
                             ta_z_zzzz_0,  \
                             ta_z_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yz_xxxx_0[i] = ta_z_xxxx_0[i] * pa_y[i] - ta_z_xxxx_1[i] * pc_y[i];

        ta_yz_xxxy_0[i] = ta_y_xxxy_0[i] * pa_z[i] - ta_y_xxxy_1[i] * pc_z[i];

        ta_yz_xxxz_0[i] = ta_z_xxxz_0[i] * pa_y[i] - ta_z_xxxz_1[i] * pc_y[i];

        ta_yz_xxyy_0[i] = ta_y_xxyy_0[i] * pa_z[i] - ta_y_xxyy_1[i] * pc_z[i];

        ta_yz_xxyz_0[i] = ta_z_xxz_0[i] * fe_0 - ta_z_xxz_1[i] * fe_0 + ta_z_xxyz_0[i] * pa_y[i] - ta_z_xxyz_1[i] * pc_y[i];

        ta_yz_xxzz_0[i] = ta_z_xxzz_0[i] * pa_y[i] - ta_z_xxzz_1[i] * pc_y[i];

        ta_yz_xyyy_0[i] = ta_y_xyyy_0[i] * pa_z[i] - ta_y_xyyy_1[i] * pc_z[i];

        ta_yz_xyyz_0[i] = 2.0 * ta_z_xyz_0[i] * fe_0 - 2.0 * ta_z_xyz_1[i] * fe_0 + ta_z_xyyz_0[i] * pa_y[i] - ta_z_xyyz_1[i] * pc_y[i];

        ta_yz_xyzz_0[i] = ta_z_xzz_0[i] * fe_0 - ta_z_xzz_1[i] * fe_0 + ta_z_xyzz_0[i] * pa_y[i] - ta_z_xyzz_1[i] * pc_y[i];

        ta_yz_xzzz_0[i] = ta_z_xzzz_0[i] * pa_y[i] - ta_z_xzzz_1[i] * pc_y[i];

        ta_yz_yyyy_0[i] = ta_y_yyyy_0[i] * pa_z[i] - ta_y_yyyy_1[i] * pc_z[i];

        ta_yz_yyyz_0[i] = 3.0 * ta_z_yyz_0[i] * fe_0 - 3.0 * ta_z_yyz_1[i] * fe_0 + ta_z_yyyz_0[i] * pa_y[i] - ta_z_yyyz_1[i] * pc_y[i];

        ta_yz_yyzz_0[i] = 2.0 * ta_z_yzz_0[i] * fe_0 - 2.0 * ta_z_yzz_1[i] * fe_0 + ta_z_yyzz_0[i] * pa_y[i] - ta_z_yyzz_1[i] * pc_y[i];

        ta_yz_yzzz_0[i] = ta_z_zzz_0[i] * fe_0 - ta_z_zzz_1[i] * fe_0 + ta_z_yzzz_0[i] * pa_y[i] - ta_z_yzzz_1[i] * pc_y[i];

        ta_yz_zzzz_0[i] = ta_z_zzzz_0[i] * pa_y[i] - ta_z_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : DG

    auto ta_zz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 75);

    auto ta_zz_xxxy_0 = pbuffer.data(idx_npot_0_dg + 76);

    auto ta_zz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 77);

    auto ta_zz_xxyy_0 = pbuffer.data(idx_npot_0_dg + 78);

    auto ta_zz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 79);

    auto ta_zz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 80);

    auto ta_zz_xyyy_0 = pbuffer.data(idx_npot_0_dg + 81);

    auto ta_zz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 82);

    auto ta_zz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 83);

    auto ta_zz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 84);

    auto ta_zz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 85);

    auto ta_zz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 86);

    auto ta_zz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 87);

    auto ta_zz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 88);

    auto ta_zz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 89);

#pragma omp simd aligned(pa_z,             \
                             pc_z,         \
                             ta_0_xxxx_0,  \
                             ta_0_xxxx_1,  \
                             ta_0_xxxy_0,  \
                             ta_0_xxxy_1,  \
                             ta_0_xxxz_0,  \
                             ta_0_xxxz_1,  \
                             ta_0_xxyy_0,  \
                             ta_0_xxyy_1,  \
                             ta_0_xxyz_0,  \
                             ta_0_xxyz_1,  \
                             ta_0_xxzz_0,  \
                             ta_0_xxzz_1,  \
                             ta_0_xyyy_0,  \
                             ta_0_xyyy_1,  \
                             ta_0_xyyz_0,  \
                             ta_0_xyyz_1,  \
                             ta_0_xyzz_0,  \
                             ta_0_xyzz_1,  \
                             ta_0_xzzz_0,  \
                             ta_0_xzzz_1,  \
                             ta_0_yyyy_0,  \
                             ta_0_yyyy_1,  \
                             ta_0_yyyz_0,  \
                             ta_0_yyyz_1,  \
                             ta_0_yyzz_0,  \
                             ta_0_yyzz_1,  \
                             ta_0_yzzz_0,  \
                             ta_0_yzzz_1,  \
                             ta_0_zzzz_0,  \
                             ta_0_zzzz_1,  \
                             ta_z_xxx_0,   \
                             ta_z_xxx_1,   \
                             ta_z_xxxx_0,  \
                             ta_z_xxxx_1,  \
                             ta_z_xxxy_0,  \
                             ta_z_xxxy_1,  \
                             ta_z_xxxz_0,  \
                             ta_z_xxxz_1,  \
                             ta_z_xxy_0,   \
                             ta_z_xxy_1,   \
                             ta_z_xxyy_0,  \
                             ta_z_xxyy_1,  \
                             ta_z_xxyz_0,  \
                             ta_z_xxyz_1,  \
                             ta_z_xxz_0,   \
                             ta_z_xxz_1,   \
                             ta_z_xxzz_0,  \
                             ta_z_xxzz_1,  \
                             ta_z_xyy_0,   \
                             ta_z_xyy_1,   \
                             ta_z_xyyy_0,  \
                             ta_z_xyyy_1,  \
                             ta_z_xyyz_0,  \
                             ta_z_xyyz_1,  \
                             ta_z_xyz_0,   \
                             ta_z_xyz_1,   \
                             ta_z_xyzz_0,  \
                             ta_z_xyzz_1,  \
                             ta_z_xzz_0,   \
                             ta_z_xzz_1,   \
                             ta_z_xzzz_0,  \
                             ta_z_xzzz_1,  \
                             ta_z_yyy_0,   \
                             ta_z_yyy_1,   \
                             ta_z_yyyy_0,  \
                             ta_z_yyyy_1,  \
                             ta_z_yyyz_0,  \
                             ta_z_yyyz_1,  \
                             ta_z_yyz_0,   \
                             ta_z_yyz_1,   \
                             ta_z_yyzz_0,  \
                             ta_z_yyzz_1,  \
                             ta_z_yzz_0,   \
                             ta_z_yzz_1,   \
                             ta_z_yzzz_0,  \
                             ta_z_yzzz_1,  \
                             ta_z_zzz_0,   \
                             ta_z_zzz_1,   \
                             ta_z_zzzz_0,  \
                             ta_z_zzzz_1,  \
                             ta_zz_xxxx_0, \
                             ta_zz_xxxy_0, \
                             ta_zz_xxxz_0, \
                             ta_zz_xxyy_0, \
                             ta_zz_xxyz_0, \
                             ta_zz_xxzz_0, \
                             ta_zz_xyyy_0, \
                             ta_zz_xyyz_0, \
                             ta_zz_xyzz_0, \
                             ta_zz_xzzz_0, \
                             ta_zz_yyyy_0, \
                             ta_zz_yyyz_0, \
                             ta_zz_yyzz_0, \
                             ta_zz_yzzz_0, \
                             ta_zz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zz_xxxx_0[i] = ta_0_xxxx_0[i] * fe_0 - ta_0_xxxx_1[i] * fe_0 + ta_z_xxxx_0[i] * pa_z[i] - ta_z_xxxx_1[i] * pc_z[i];

        ta_zz_xxxy_0[i] = ta_0_xxxy_0[i] * fe_0 - ta_0_xxxy_1[i] * fe_0 + ta_z_xxxy_0[i] * pa_z[i] - ta_z_xxxy_1[i] * pc_z[i];

        ta_zz_xxxz_0[i] = ta_0_xxxz_0[i] * fe_0 - ta_0_xxxz_1[i] * fe_0 + ta_z_xxx_0[i] * fe_0 - ta_z_xxx_1[i] * fe_0 + ta_z_xxxz_0[i] * pa_z[i] -
                          ta_z_xxxz_1[i] * pc_z[i];

        ta_zz_xxyy_0[i] = ta_0_xxyy_0[i] * fe_0 - ta_0_xxyy_1[i] * fe_0 + ta_z_xxyy_0[i] * pa_z[i] - ta_z_xxyy_1[i] * pc_z[i];

        ta_zz_xxyz_0[i] = ta_0_xxyz_0[i] * fe_0 - ta_0_xxyz_1[i] * fe_0 + ta_z_xxy_0[i] * fe_0 - ta_z_xxy_1[i] * fe_0 + ta_z_xxyz_0[i] * pa_z[i] -
                          ta_z_xxyz_1[i] * pc_z[i];

        ta_zz_xxzz_0[i] = ta_0_xxzz_0[i] * fe_0 - ta_0_xxzz_1[i] * fe_0 + 2.0 * ta_z_xxz_0[i] * fe_0 - 2.0 * ta_z_xxz_1[i] * fe_0 +
                          ta_z_xxzz_0[i] * pa_z[i] - ta_z_xxzz_1[i] * pc_z[i];

        ta_zz_xyyy_0[i] = ta_0_xyyy_0[i] * fe_0 - ta_0_xyyy_1[i] * fe_0 + ta_z_xyyy_0[i] * pa_z[i] - ta_z_xyyy_1[i] * pc_z[i];

        ta_zz_xyyz_0[i] = ta_0_xyyz_0[i] * fe_0 - ta_0_xyyz_1[i] * fe_0 + ta_z_xyy_0[i] * fe_0 - ta_z_xyy_1[i] * fe_0 + ta_z_xyyz_0[i] * pa_z[i] -
                          ta_z_xyyz_1[i] * pc_z[i];

        ta_zz_xyzz_0[i] = ta_0_xyzz_0[i] * fe_0 - ta_0_xyzz_1[i] * fe_0 + 2.0 * ta_z_xyz_0[i] * fe_0 - 2.0 * ta_z_xyz_1[i] * fe_0 +
                          ta_z_xyzz_0[i] * pa_z[i] - ta_z_xyzz_1[i] * pc_z[i];

        ta_zz_xzzz_0[i] = ta_0_xzzz_0[i] * fe_0 - ta_0_xzzz_1[i] * fe_0 + 3.0 * ta_z_xzz_0[i] * fe_0 - 3.0 * ta_z_xzz_1[i] * fe_0 +
                          ta_z_xzzz_0[i] * pa_z[i] - ta_z_xzzz_1[i] * pc_z[i];

        ta_zz_yyyy_0[i] = ta_0_yyyy_0[i] * fe_0 - ta_0_yyyy_1[i] * fe_0 + ta_z_yyyy_0[i] * pa_z[i] - ta_z_yyyy_1[i] * pc_z[i];

        ta_zz_yyyz_0[i] = ta_0_yyyz_0[i] * fe_0 - ta_0_yyyz_1[i] * fe_0 + ta_z_yyy_0[i] * fe_0 - ta_z_yyy_1[i] * fe_0 + ta_z_yyyz_0[i] * pa_z[i] -
                          ta_z_yyyz_1[i] * pc_z[i];

        ta_zz_yyzz_0[i] = ta_0_yyzz_0[i] * fe_0 - ta_0_yyzz_1[i] * fe_0 + 2.0 * ta_z_yyz_0[i] * fe_0 - 2.0 * ta_z_yyz_1[i] * fe_0 +
                          ta_z_yyzz_0[i] * pa_z[i] - ta_z_yyzz_1[i] * pc_z[i];

        ta_zz_yzzz_0[i] = ta_0_yzzz_0[i] * fe_0 - ta_0_yzzz_1[i] * fe_0 + 3.0 * ta_z_yzz_0[i] * fe_0 - 3.0 * ta_z_yzz_1[i] * fe_0 +
                          ta_z_yzzz_0[i] * pa_z[i] - ta_z_yzzz_1[i] * pc_z[i];

        ta_zz_zzzz_0[i] = ta_0_zzzz_0[i] * fe_0 - ta_0_zzzz_1[i] * fe_0 + 4.0 * ta_z_zzz_0[i] * fe_0 - 4.0 * ta_z_zzz_1[i] * fe_0 +
                          ta_z_zzzz_0[i] * pa_z[i] - ta_z_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
