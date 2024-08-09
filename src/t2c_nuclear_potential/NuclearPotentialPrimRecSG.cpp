#include "NuclearPotentialPrimRecSG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_sg(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_sg,
                               const size_t              idx_npot_0_sd,
                               const size_t              idx_npot_1_sd,
                               const size_t              idx_npot_0_sf,
                               const size_t              idx_npot_1_sf,
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

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_0 = pbuffer.data(idx_npot_0_sd);

    auto ta_0_yy_0 = pbuffer.data(idx_npot_0_sd + 3);

    auto ta_0_zz_0 = pbuffer.data(idx_npot_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = pbuffer.data(idx_npot_1_sd);

    auto ta_0_yy_1 = pbuffer.data(idx_npot_1_sd + 3);

    auto ta_0_zz_1 = pbuffer.data(idx_npot_1_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_0 = pbuffer.data(idx_npot_0_sf);

    auto ta_0_xxz_0 = pbuffer.data(idx_npot_0_sf + 2);

    auto ta_0_xyy_0 = pbuffer.data(idx_npot_0_sf + 3);

    auto ta_0_xzz_0 = pbuffer.data(idx_npot_0_sf + 5);

    auto ta_0_yyy_0 = pbuffer.data(idx_npot_0_sf + 6);

    auto ta_0_yyz_0 = pbuffer.data(idx_npot_0_sf + 7);

    auto ta_0_yzz_0 = pbuffer.data(idx_npot_0_sf + 8);

    auto ta_0_zzz_0 = pbuffer.data(idx_npot_0_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_1 = pbuffer.data(idx_npot_1_sf);

    auto ta_0_xxz_1 = pbuffer.data(idx_npot_1_sf + 2);

    auto ta_0_xyy_1 = pbuffer.data(idx_npot_1_sf + 3);

    auto ta_0_xzz_1 = pbuffer.data(idx_npot_1_sf + 5);

    auto ta_0_yyy_1 = pbuffer.data(idx_npot_1_sf + 6);

    auto ta_0_yyz_1 = pbuffer.data(idx_npot_1_sf + 7);

    auto ta_0_yzz_1 = pbuffer.data(idx_npot_1_sf + 8);

    auto ta_0_zzz_1 = pbuffer.data(idx_npot_1_sf + 9);

    // Set up components of targeted buffer : SG

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

#pragma omp simd aligned(pb_x,            \
                             pb_y,        \
                             pb_z,        \
                             pc_x,        \
                             pc_y,        \
                             pc_z,        \
                             ta_0_xx_0,   \
                             ta_0_xx_1,   \
                             ta_0_xxx_0,  \
                             ta_0_xxx_1,  \
                             ta_0_xxxx_0, \
                             ta_0_xxxy_0, \
                             ta_0_xxxz_0, \
                             ta_0_xxyy_0, \
                             ta_0_xxyz_0, \
                             ta_0_xxz_0,  \
                             ta_0_xxz_1,  \
                             ta_0_xxzz_0, \
                             ta_0_xyy_0,  \
                             ta_0_xyy_1,  \
                             ta_0_xyyy_0, \
                             ta_0_xyyz_0, \
                             ta_0_xyzz_0, \
                             ta_0_xzz_0,  \
                             ta_0_xzz_1,  \
                             ta_0_xzzz_0, \
                             ta_0_yy_0,   \
                             ta_0_yy_1,   \
                             ta_0_yyy_0,  \
                             ta_0_yyy_1,  \
                             ta_0_yyyy_0, \
                             ta_0_yyyz_0, \
                             ta_0_yyz_0,  \
                             ta_0_yyz_1,  \
                             ta_0_yyzz_0, \
                             ta_0_yzz_0,  \
                             ta_0_yzz_1,  \
                             ta_0_yzzz_0, \
                             ta_0_zz_0,   \
                             ta_0_zz_1,   \
                             ta_0_zzz_0,  \
                             ta_0_zzz_1,  \
                             ta_0_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_0_xxxx_0[i] =
            3.0 * ta_0_xx_0[i] * fe_0 - 3.0 * ta_0_xx_1[i] * fe_0 + ta_0_xxx_0[i] * pb_x[i] - ta_0_xxx_1[i] * pc_x[i];

        ta_0_xxxy_0[i] = ta_0_xxx_0[i] * pb_y[i] - ta_0_xxx_1[i] * pc_y[i];

        ta_0_xxxz_0[i] = ta_0_xxx_0[i] * pb_z[i] - ta_0_xxx_1[i] * pc_z[i];

        ta_0_xxyy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_0_xyy_0[i] * pb_x[i] - ta_0_xyy_1[i] * pc_x[i];

        ta_0_xxyz_0[i] = ta_0_xxz_0[i] * pb_y[i] - ta_0_xxz_1[i] * pc_y[i];

        ta_0_xxzz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_0_xzz_0[i] * pb_x[i] - ta_0_xzz_1[i] * pc_x[i];

        ta_0_xyyy_0[i] = ta_0_yyy_0[i] * pb_x[i] - ta_0_yyy_1[i] * pc_x[i];

        ta_0_xyyz_0[i] = ta_0_yyz_0[i] * pb_x[i] - ta_0_yyz_1[i] * pc_x[i];

        ta_0_xyzz_0[i] = ta_0_yzz_0[i] * pb_x[i] - ta_0_yzz_1[i] * pc_x[i];

        ta_0_xzzz_0[i] = ta_0_zzz_0[i] * pb_x[i] - ta_0_zzz_1[i] * pc_x[i];

        ta_0_yyyy_0[i] =
            3.0 * ta_0_yy_0[i] * fe_0 - 3.0 * ta_0_yy_1[i] * fe_0 + ta_0_yyy_0[i] * pb_y[i] - ta_0_yyy_1[i] * pc_y[i];

        ta_0_yyyz_0[i] = ta_0_yyy_0[i] * pb_z[i] - ta_0_yyy_1[i] * pc_z[i];

        ta_0_yyzz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_0_yzz_0[i] * pb_y[i] - ta_0_yzz_1[i] * pc_y[i];

        ta_0_yzzz_0[i] = ta_0_zzz_0[i] * pb_y[i] - ta_0_zzz_1[i] * pc_y[i];

        ta_0_zzzz_0[i] =
            3.0 * ta_0_zz_0[i] * fe_0 - 3.0 * ta_0_zz_1[i] * fe_0 + ta_0_zzz_0[i] * pb_z[i] - ta_0_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
