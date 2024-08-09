#include "NuclearPotentialPrimRecSF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_sf(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_sf,
                               const size_t              idx_npot_0_sp,
                               const size_t              idx_npot_1_sp,
                               const size_t              idx_npot_0_sd,
                               const size_t              idx_npot_1_sd,
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

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_0 = pbuffer.data(idx_npot_0_sp);

    auto ta_0_y_0 = pbuffer.data(idx_npot_0_sp + 1);

    auto ta_0_z_0 = pbuffer.data(idx_npot_0_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_1 = pbuffer.data(idx_npot_1_sp);

    auto ta_0_y_1 = pbuffer.data(idx_npot_1_sp + 1);

    auto ta_0_z_1 = pbuffer.data(idx_npot_1_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_0 = pbuffer.data(idx_npot_0_sd);

    auto ta_0_yy_0 = pbuffer.data(idx_npot_0_sd + 3);

    auto ta_0_yz_0 = pbuffer.data(idx_npot_0_sd + 4);

    auto ta_0_zz_0 = pbuffer.data(idx_npot_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = pbuffer.data(idx_npot_1_sd);

    auto ta_0_yy_1 = pbuffer.data(idx_npot_1_sd + 3);

    auto ta_0_yz_1 = pbuffer.data(idx_npot_1_sd + 4);

    auto ta_0_zz_1 = pbuffer.data(idx_npot_1_sd + 5);

    // Set up components of targeted buffer : SF

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

#pragma omp simd aligned(pb_x,           \
                             pb_y,       \
                             pb_z,       \
                             pc_x,       \
                             pc_y,       \
                             pc_z,       \
                             ta_0_x_0,   \
                             ta_0_x_1,   \
                             ta_0_xx_0,  \
                             ta_0_xx_1,  \
                             ta_0_xxx_0, \
                             ta_0_xxy_0, \
                             ta_0_xxz_0, \
                             ta_0_xyy_0, \
                             ta_0_xyz_0, \
                             ta_0_xzz_0, \
                             ta_0_y_0,   \
                             ta_0_y_1,   \
                             ta_0_yy_0,  \
                             ta_0_yy_1,  \
                             ta_0_yyy_0, \
                             ta_0_yyz_0, \
                             ta_0_yz_0,  \
                             ta_0_yz_1,  \
                             ta_0_yzz_0, \
                             ta_0_z_0,   \
                             ta_0_z_1,   \
                             ta_0_zz_0,  \
                             ta_0_zz_1,  \
                             ta_0_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_0_xxx_0[i] =
            2.0 * ta_0_x_0[i] * fe_0 - 2.0 * ta_0_x_1[i] * fe_0 + ta_0_xx_0[i] * pb_x[i] - ta_0_xx_1[i] * pc_x[i];

        ta_0_xxy_0[i] = ta_0_xx_0[i] * pb_y[i] - ta_0_xx_1[i] * pc_y[i];

        ta_0_xxz_0[i] = ta_0_xx_0[i] * pb_z[i] - ta_0_xx_1[i] * pc_z[i];

        ta_0_xyy_0[i] = ta_0_yy_0[i] * pb_x[i] - ta_0_yy_1[i] * pc_x[i];

        ta_0_xyz_0[i] = ta_0_yz_0[i] * pb_x[i] - ta_0_yz_1[i] * pc_x[i];

        ta_0_xzz_0[i] = ta_0_zz_0[i] * pb_x[i] - ta_0_zz_1[i] * pc_x[i];

        ta_0_yyy_0[i] =
            2.0 * ta_0_y_0[i] * fe_0 - 2.0 * ta_0_y_1[i] * fe_0 + ta_0_yy_0[i] * pb_y[i] - ta_0_yy_1[i] * pc_y[i];

        ta_0_yyz_0[i] = ta_0_yy_0[i] * pb_z[i] - ta_0_yy_1[i] * pc_z[i];

        ta_0_yzz_0[i] = ta_0_zz_0[i] * pb_y[i] - ta_0_zz_1[i] * pc_y[i];

        ta_0_zzz_0[i] =
            2.0 * ta_0_z_0[i] * fe_0 - 2.0 * ta_0_z_1[i] * fe_0 + ta_0_zz_0[i] * pb_z[i] - ta_0_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
