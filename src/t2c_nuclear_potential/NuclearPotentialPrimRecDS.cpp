#include "NuclearPotentialPrimRecDS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_ds(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_ds,
                               const size_t              idx_npot_0_ss,
                               const size_t              idx_npot_1_ss,
                               const size_t              idx_npot_0_ps,
                               const size_t              idx_npot_1_ps,
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

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_0 = pbuffer.data(idx_npot_0_ss);

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_1 = pbuffer.data(idx_npot_1_ss);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_0 = pbuffer.data(idx_npot_0_ps);

    auto ta_y_0_0 = pbuffer.data(idx_npot_0_ps + 1);

    auto ta_z_0_0 = pbuffer.data(idx_npot_0_ps + 2);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_1 = pbuffer.data(idx_npot_1_ps);

    auto ta_y_0_1 = pbuffer.data(idx_npot_1_ps + 1);

    auto ta_z_0_1 = pbuffer.data(idx_npot_1_ps + 2);

    // Set up components of targeted buffer : DS

    auto ta_xx_0_0 = pbuffer.data(idx_npot_0_ds);

    auto ta_xy_0_0 = pbuffer.data(idx_npot_0_ds + 1);

    auto ta_xz_0_0 = pbuffer.data(idx_npot_0_ds + 2);

    auto ta_yy_0_0 = pbuffer.data(idx_npot_0_ds + 3);

    auto ta_yz_0_0 = pbuffer.data(idx_npot_0_ds + 4);

    auto ta_zz_0_0 = pbuffer.data(idx_npot_0_ds + 5);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             pa_z,      \
                             pc_x,      \
                             pc_y,      \
                             pc_z,      \
                             ta_0_0_0,  \
                             ta_0_0_1,  \
                             ta_x_0_0,  \
                             ta_x_0_1,  \
                             ta_xx_0_0, \
                             ta_xy_0_0, \
                             ta_xz_0_0, \
                             ta_y_0_0,  \
                             ta_y_0_1,  \
                             ta_yy_0_0, \
                             ta_yz_0_0, \
                             ta_z_0_0,  \
                             ta_z_0_1,  \
                             ta_zz_0_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xx_0_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_x_0_0[i] * pa_x[i] - ta_x_0_1[i] * pc_x[i];

        ta_xy_0_0[i] = ta_y_0_0[i] * pa_x[i] - ta_y_0_1[i] * pc_x[i];

        ta_xz_0_0[i] = ta_z_0_0[i] * pa_x[i] - ta_z_0_1[i] * pc_x[i];

        ta_yy_0_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_y_0_0[i] * pa_y[i] - ta_y_0_1[i] * pc_y[i];

        ta_yz_0_0[i] = ta_z_0_0[i] * pa_y[i] - ta_z_0_1[i] * pc_y[i];

        ta_zz_0_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_z_0_0[i] * pa_z[i] - ta_z_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
