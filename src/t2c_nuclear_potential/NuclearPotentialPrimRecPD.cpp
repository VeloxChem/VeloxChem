#include "NuclearPotentialPrimRecPD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_pd(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_pd,
                               const size_t              idx_npot_0_sp,
                               const size_t              idx_npot_1_sp,
                               const size_t              idx_npot_0_sd,
                               const size_t              idx_npot_1_sd,
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

    auto ta_0_xy_0 = pbuffer.data(idx_npot_0_sd + 1);

    auto ta_0_xz_0 = pbuffer.data(idx_npot_0_sd + 2);

    auto ta_0_yy_0 = pbuffer.data(idx_npot_0_sd + 3);

    auto ta_0_yz_0 = pbuffer.data(idx_npot_0_sd + 4);

    auto ta_0_zz_0 = pbuffer.data(idx_npot_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = pbuffer.data(idx_npot_1_sd);

    auto ta_0_xy_1 = pbuffer.data(idx_npot_1_sd + 1);

    auto ta_0_xz_1 = pbuffer.data(idx_npot_1_sd + 2);

    auto ta_0_yy_1 = pbuffer.data(idx_npot_1_sd + 3);

    auto ta_0_yz_1 = pbuffer.data(idx_npot_1_sd + 4);

    auto ta_0_zz_1 = pbuffer.data(idx_npot_1_sd + 5);

    // Set up 0-6 components of targeted buffer : PD

    auto ta_x_xx_0 = pbuffer.data(idx_npot_0_pd);

    auto ta_x_xy_0 = pbuffer.data(idx_npot_0_pd + 1);

    auto ta_x_xz_0 = pbuffer.data(idx_npot_0_pd + 2);

    auto ta_x_yy_0 = pbuffer.data(idx_npot_0_pd + 3);

    auto ta_x_yz_0 = pbuffer.data(idx_npot_0_pd + 4);

    auto ta_x_zz_0 = pbuffer.data(idx_npot_0_pd + 5);

#pragma omp simd aligned(pa_x,          \
                             pc_x,      \
                             ta_0_x_0,  \
                             ta_0_x_1,  \
                             ta_0_xx_0, \
                             ta_0_xx_1, \
                             ta_0_xy_0, \
                             ta_0_xy_1, \
                             ta_0_xz_0, \
                             ta_0_xz_1, \
                             ta_0_y_0,  \
                             ta_0_y_1,  \
                             ta_0_yy_0, \
                             ta_0_yy_1, \
                             ta_0_yz_0, \
                             ta_0_yz_1, \
                             ta_0_z_0,  \
                             ta_0_z_1,  \
                             ta_0_zz_0, \
                             ta_0_zz_1, \
                             ta_x_xx_0, \
                             ta_x_xy_0, \
                             ta_x_xz_0, \
                             ta_x_yy_0, \
                             ta_x_yz_0, \
                             ta_x_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_x_xx_0[i] =
            2.0 * ta_0_x_0[i] * fe_0 - 2.0 * ta_0_x_1[i] * fe_0 + ta_0_xx_0[i] * pa_x[i] - ta_0_xx_1[i] * pc_x[i];

        ta_x_xy_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_0_xy_0[i] * pa_x[i] - ta_0_xy_1[i] * pc_x[i];

        ta_x_xz_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_0_xz_0[i] * pa_x[i] - ta_0_xz_1[i] * pc_x[i];

        ta_x_yy_0[i] = ta_0_yy_0[i] * pa_x[i] - ta_0_yy_1[i] * pc_x[i];

        ta_x_yz_0[i] = ta_0_yz_0[i] * pa_x[i] - ta_0_yz_1[i] * pc_x[i];

        ta_x_zz_0[i] = ta_0_zz_0[i] * pa_x[i] - ta_0_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto ta_y_xx_0 = pbuffer.data(idx_npot_0_pd + 6);

    auto ta_y_xy_0 = pbuffer.data(idx_npot_0_pd + 7);

    auto ta_y_xz_0 = pbuffer.data(idx_npot_0_pd + 8);

    auto ta_y_yy_0 = pbuffer.data(idx_npot_0_pd + 9);

    auto ta_y_yz_0 = pbuffer.data(idx_npot_0_pd + 10);

    auto ta_y_zz_0 = pbuffer.data(idx_npot_0_pd + 11);

#pragma omp simd aligned(pa_y,          \
                             pc_y,      \
                             ta_0_x_0,  \
                             ta_0_x_1,  \
                             ta_0_xx_0, \
                             ta_0_xx_1, \
                             ta_0_xy_0, \
                             ta_0_xy_1, \
                             ta_0_xz_0, \
                             ta_0_xz_1, \
                             ta_0_y_0,  \
                             ta_0_y_1,  \
                             ta_0_yy_0, \
                             ta_0_yy_1, \
                             ta_0_yz_0, \
                             ta_0_yz_1, \
                             ta_0_z_0,  \
                             ta_0_z_1,  \
                             ta_0_zz_0, \
                             ta_0_zz_1, \
                             ta_y_xx_0, \
                             ta_y_xy_0, \
                             ta_y_xz_0, \
                             ta_y_yy_0, \
                             ta_y_yz_0, \
                             ta_y_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_y_xx_0[i] = ta_0_xx_0[i] * pa_y[i] - ta_0_xx_1[i] * pc_y[i];

        ta_y_xy_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_0_xy_0[i] * pa_y[i] - ta_0_xy_1[i] * pc_y[i];

        ta_y_xz_0[i] = ta_0_xz_0[i] * pa_y[i] - ta_0_xz_1[i] * pc_y[i];

        ta_y_yy_0[i] =
            2.0 * ta_0_y_0[i] * fe_0 - 2.0 * ta_0_y_1[i] * fe_0 + ta_0_yy_0[i] * pa_y[i] - ta_0_yy_1[i] * pc_y[i];

        ta_y_yz_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_0_yz_0[i] * pa_y[i] - ta_0_yz_1[i] * pc_y[i];

        ta_y_zz_0[i] = ta_0_zz_0[i] * pa_y[i] - ta_0_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto ta_z_xx_0 = pbuffer.data(idx_npot_0_pd + 12);

    auto ta_z_xy_0 = pbuffer.data(idx_npot_0_pd + 13);

    auto ta_z_xz_0 = pbuffer.data(idx_npot_0_pd + 14);

    auto ta_z_yy_0 = pbuffer.data(idx_npot_0_pd + 15);

    auto ta_z_yz_0 = pbuffer.data(idx_npot_0_pd + 16);

    auto ta_z_zz_0 = pbuffer.data(idx_npot_0_pd + 17);

#pragma omp simd aligned(pa_z,          \
                             pc_z,      \
                             ta_0_x_0,  \
                             ta_0_x_1,  \
                             ta_0_xx_0, \
                             ta_0_xx_1, \
                             ta_0_xy_0, \
                             ta_0_xy_1, \
                             ta_0_xz_0, \
                             ta_0_xz_1, \
                             ta_0_y_0,  \
                             ta_0_y_1,  \
                             ta_0_yy_0, \
                             ta_0_yy_1, \
                             ta_0_yz_0, \
                             ta_0_yz_1, \
                             ta_0_z_0,  \
                             ta_0_z_1,  \
                             ta_0_zz_0, \
                             ta_0_zz_1, \
                             ta_z_xx_0, \
                             ta_z_xy_0, \
                             ta_z_xz_0, \
                             ta_z_yy_0, \
                             ta_z_yz_0, \
                             ta_z_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_z_xx_0[i] = ta_0_xx_0[i] * pa_z[i] - ta_0_xx_1[i] * pc_z[i];

        ta_z_xy_0[i] = ta_0_xy_0[i] * pa_z[i] - ta_0_xy_1[i] * pc_z[i];

        ta_z_xz_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_0_xz_0[i] * pa_z[i] - ta_0_xz_1[i] * pc_z[i];

        ta_z_yy_0[i] = ta_0_yy_0[i] * pa_z[i] - ta_0_yy_1[i] * pc_z[i];

        ta_z_yz_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_0_yz_0[i] * pa_z[i] - ta_0_yz_1[i] * pc_z[i];

        ta_z_zz_0[i] =
            2.0 * ta_0_z_0[i] * fe_0 - 2.0 * ta_0_z_1[i] * fe_0 + ta_0_zz_0[i] * pa_z[i] - ta_0_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
