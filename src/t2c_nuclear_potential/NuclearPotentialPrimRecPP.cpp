#include "NuclearPotentialPrimRecPP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_pp(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_pp,
                               const size_t              idx_npot_0_ss,
                               const size_t              idx_npot_1_ss,
                               const size_t              idx_npot_0_sp,
                               const size_t              idx_npot_1_sp,
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

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_0 = pbuffer.data(idx_npot_0_sp);

    auto ta_0_y_0 = pbuffer.data(idx_npot_0_sp + 1);

    auto ta_0_z_0 = pbuffer.data(idx_npot_0_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_1 = pbuffer.data(idx_npot_1_sp);

    auto ta_0_y_1 = pbuffer.data(idx_npot_1_sp + 1);

    auto ta_0_z_1 = pbuffer.data(idx_npot_1_sp + 2);

    // Set up 0-3 components of targeted buffer : PP

    auto ta_x_x_0 = pbuffer.data(idx_npot_0_pp);

    auto ta_x_y_0 = pbuffer.data(idx_npot_0_pp + 1);

    auto ta_x_z_0 = pbuffer.data(idx_npot_0_pp + 2);

#pragma omp simd aligned(pa_x,         \
                             pc_x,     \
                             ta_0_0_0, \
                             ta_0_0_1, \
                             ta_0_x_0, \
                             ta_0_x_1, \
                             ta_0_y_0, \
                             ta_0_y_1, \
                             ta_0_z_0, \
                             ta_0_z_1, \
                             ta_x_x_0, \
                             ta_x_y_0, \
                             ta_x_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_x_x_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_x_0[i] * pa_x[i] - ta_0_x_1[i] * pc_x[i];

        ta_x_y_0[i] = ta_0_y_0[i] * pa_x[i] - ta_0_y_1[i] * pc_x[i];

        ta_x_z_0[i] = ta_0_z_0[i] * pa_x[i] - ta_0_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto ta_y_x_0 = pbuffer.data(idx_npot_0_pp + 3);

    auto ta_y_y_0 = pbuffer.data(idx_npot_0_pp + 4);

    auto ta_y_z_0 = pbuffer.data(idx_npot_0_pp + 5);

#pragma omp simd aligned(pa_y,         \
                             pc_y,     \
                             ta_0_0_0, \
                             ta_0_0_1, \
                             ta_0_x_0, \
                             ta_0_x_1, \
                             ta_0_y_0, \
                             ta_0_y_1, \
                             ta_0_z_0, \
                             ta_0_z_1, \
                             ta_y_x_0, \
                             ta_y_y_0, \
                             ta_y_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_y_x_0[i] = ta_0_x_0[i] * pa_y[i] - ta_0_x_1[i] * pc_y[i];

        ta_y_y_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_y_0[i] * pa_y[i] - ta_0_y_1[i] * pc_y[i];

        ta_y_z_0[i] = ta_0_z_0[i] * pa_y[i] - ta_0_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto ta_z_x_0 = pbuffer.data(idx_npot_0_pp + 6);

    auto ta_z_y_0 = pbuffer.data(idx_npot_0_pp + 7);

    auto ta_z_z_0 = pbuffer.data(idx_npot_0_pp + 8);

#pragma omp simd aligned(pa_z,         \
                             pc_z,     \
                             ta_0_0_0, \
                             ta_0_0_1, \
                             ta_0_x_0, \
                             ta_0_x_1, \
                             ta_0_y_0, \
                             ta_0_y_1, \
                             ta_0_z_0, \
                             ta_0_z_1, \
                             ta_z_x_0, \
                             ta_z_y_0, \
                             ta_z_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_z_x_0[i] = ta_0_x_0[i] * pa_z[i] - ta_0_x_1[i] * pc_z[i];

        ta_z_y_0[i] = ta_0_y_0[i] * pa_z[i] - ta_0_y_1[i] * pc_z[i];

        ta_z_z_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_z_0[i] * pa_z[i] - ta_0_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
