#include "NuclearPotentialPrimRecSD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_sd(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_sd,
                               const size_t              idx_npot_0_ss,
                               const size_t              idx_npot_1_ss,
                               const size_t              idx_npot_0_sp,
                               const size_t              idx_npot_1_sp,
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

    auto ta_0_0_0 = pbuffer.data(idx_npot_0_ss);

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_1 = pbuffer.data(idx_npot_1_ss);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_0 = pbuffer.data(idx_npot_0_sp);

    auto ta_0_y_0 = pbuffer.data(idx_npot_0_sp + 1);

    auto ta_0_z_0 = pbuffer.data(idx_npot_0_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_1 = pbuffer.data(idx_npot_1_sp);

    auto ta_0_y_1 = pbuffer.data(idx_npot_1_sp + 1);

    auto ta_0_z_1 = pbuffer.data(idx_npot_1_sp + 2);

    // Set up components of targeted buffer : SD

    auto ta_0_xx_0 = pbuffer.data(idx_npot_0_sd);

    auto ta_0_xy_0 = pbuffer.data(idx_npot_0_sd + 1);

    auto ta_0_xz_0 = pbuffer.data(idx_npot_0_sd + 2);

    auto ta_0_yy_0 = pbuffer.data(idx_npot_0_sd + 3);

    auto ta_0_yz_0 = pbuffer.data(idx_npot_0_sd + 4);

    auto ta_0_zz_0 = pbuffer.data(idx_npot_0_sd + 5);

#pragma omp simd aligned(pb_x,          \
                             pb_y,      \
                             pb_z,      \
                             pc_x,      \
                             pc_y,      \
                             pc_z,      \
                             ta_0_0_0,  \
                             ta_0_0_1,  \
                             ta_0_x_0,  \
                             ta_0_x_1,  \
                             ta_0_xx_0, \
                             ta_0_xy_0, \
                             ta_0_xz_0, \
                             ta_0_y_0,  \
                             ta_0_y_1,  \
                             ta_0_yy_0, \
                             ta_0_yz_0, \
                             ta_0_z_0,  \
                             ta_0_z_1,  \
                             ta_0_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_0_xx_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_x_0[i] * pb_x[i] - ta_0_x_1[i] * pc_x[i];

        ta_0_xy_0[i] = ta_0_y_0[i] * pb_x[i] - ta_0_y_1[i] * pc_x[i];

        ta_0_xz_0[i] = ta_0_z_0[i] * pb_x[i] - ta_0_z_1[i] * pc_x[i];

        ta_0_yy_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_y_0[i] * pb_y[i] - ta_0_y_1[i] * pc_y[i];

        ta_0_yz_0[i] = ta_0_z_0[i] * pb_y[i] - ta_0_z_1[i] * pc_y[i];

        ta_0_zz_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_z_0[i] * pb_z[i] - ta_0_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec