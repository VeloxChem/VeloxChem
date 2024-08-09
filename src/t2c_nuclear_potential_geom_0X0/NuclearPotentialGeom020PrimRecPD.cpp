#include "NuclearPotentialGeom020PrimRecPD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_pd(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_pd,
                                        const size_t idx_npot_geom_020_0_sp,
                                        const size_t idx_npot_geom_020_1_sp,
                                        const size_t idx_npot_geom_010_1_sd,
                                        const size_t idx_npot_geom_020_0_sd,
                                        const size_t idx_npot_geom_020_1_sd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void
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

    auto ta2_xx_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp);

    auto ta2_xx_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 1);

    auto ta2_xx_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 2);

    auto ta2_xy_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 3);

    auto ta2_xy_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 4);

    auto ta2_xy_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 5);

    auto ta2_xz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 6);

    auto ta2_xz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 7);

    auto ta2_xz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 8);

    auto ta2_yy_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 9);

    auto ta2_yy_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 10);

    auto ta2_yy_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 11);

    auto ta2_yz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 12);

    auto ta2_yz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 13);

    auto ta2_yz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 14);

    auto ta2_zz_0_x_0 = pbuffer.data(idx_npot_geom_020_0_sp + 15);

    auto ta2_zz_0_y_0 = pbuffer.data(idx_npot_geom_020_0_sp + 16);

    auto ta2_zz_0_z_0 = pbuffer.data(idx_npot_geom_020_0_sp + 17);

    // Set up components of auxiliary buffer : SP

    auto ta2_xx_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp);

    auto ta2_xx_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 1);

    auto ta2_xx_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 2);

    auto ta2_xy_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 3);

    auto ta2_xy_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 4);

    auto ta2_xy_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 5);

    auto ta2_xz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 6);

    auto ta2_xz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 7);

    auto ta2_xz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 8);

    auto ta2_yy_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 9);

    auto ta2_yy_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 10);

    auto ta2_yy_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 11);

    auto ta2_yz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 12);

    auto ta2_yz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 13);

    auto ta2_yz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 14);

    auto ta2_zz_0_x_1 = pbuffer.data(idx_npot_geom_020_1_sp + 15);

    auto ta2_zz_0_y_1 = pbuffer.data(idx_npot_geom_020_1_sp + 16);

    auto ta2_zz_0_z_1 = pbuffer.data(idx_npot_geom_020_1_sp + 17);

    // Set up components of auxiliary buffer : SD

    auto ta1_x_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd);

    auto ta1_x_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 1);

    auto ta1_x_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 2);

    auto ta1_x_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 3);

    auto ta1_x_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 4);

    auto ta1_x_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 5);

    auto ta1_y_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 6);

    auto ta1_y_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 7);

    auto ta1_y_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 8);

    auto ta1_y_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 9);

    auto ta1_y_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 10);

    auto ta1_y_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 11);

    auto ta1_z_0_xx_1 = pbuffer.data(idx_npot_geom_010_1_sd + 12);

    auto ta1_z_0_xy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 13);

    auto ta1_z_0_xz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 14);

    auto ta1_z_0_yy_1 = pbuffer.data(idx_npot_geom_010_1_sd + 15);

    auto ta1_z_0_yz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 16);

    auto ta1_z_0_zz_1 = pbuffer.data(idx_npot_geom_010_1_sd + 17);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd);

    auto ta2_xx_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 1);

    auto ta2_xx_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 2);

    auto ta2_xx_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 3);

    auto ta2_xx_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 4);

    auto ta2_xx_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 5);

    auto ta2_xy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 6);

    auto ta2_xy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 7);

    auto ta2_xy_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 8);

    auto ta2_xy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 9);

    auto ta2_xy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 10);

    auto ta2_xy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 11);

    auto ta2_xz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 12);

    auto ta2_xz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 13);

    auto ta2_xz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 14);

    auto ta2_xz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 15);

    auto ta2_xz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 16);

    auto ta2_xz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 17);

    auto ta2_yy_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 18);

    auto ta2_yy_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 19);

    auto ta2_yy_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 20);

    auto ta2_yy_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 21);

    auto ta2_yy_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 22);

    auto ta2_yy_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 23);

    auto ta2_yz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 24);

    auto ta2_yz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 25);

    auto ta2_yz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 26);

    auto ta2_yz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 27);

    auto ta2_yz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 28);

    auto ta2_yz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 29);

    auto ta2_zz_0_xx_0 = pbuffer.data(idx_npot_geom_020_0_sd + 30);

    auto ta2_zz_0_xy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 31);

    auto ta2_zz_0_xz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 32);

    auto ta2_zz_0_yy_0 = pbuffer.data(idx_npot_geom_020_0_sd + 33);

    auto ta2_zz_0_yz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 34);

    auto ta2_zz_0_zz_0 = pbuffer.data(idx_npot_geom_020_0_sd + 35);

    // Set up components of auxiliary buffer : SD

    auto ta2_xx_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd);

    auto ta2_xx_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 1);

    auto ta2_xx_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 2);

    auto ta2_xx_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 3);

    auto ta2_xx_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 4);

    auto ta2_xx_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 5);

    auto ta2_xy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 6);

    auto ta2_xy_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 7);

    auto ta2_xy_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 8);

    auto ta2_xy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 9);

    auto ta2_xy_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 10);

    auto ta2_xy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 11);

    auto ta2_xz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 12);

    auto ta2_xz_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 13);

    auto ta2_xz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 14);

    auto ta2_xz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 15);

    auto ta2_xz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 16);

    auto ta2_xz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 17);

    auto ta2_yy_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 18);

    auto ta2_yy_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 19);

    auto ta2_yy_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 20);

    auto ta2_yy_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 21);

    auto ta2_yy_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 22);

    auto ta2_yy_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 23);

    auto ta2_yz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 24);

    auto ta2_yz_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 25);

    auto ta2_yz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 26);

    auto ta2_yz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 27);

    auto ta2_yz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 28);

    auto ta2_yz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 29);

    auto ta2_zz_0_xx_1 = pbuffer.data(idx_npot_geom_020_1_sd + 30);

    auto ta2_zz_0_xy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 31);

    auto ta2_zz_0_xz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 32);

    auto ta2_zz_0_yy_1 = pbuffer.data(idx_npot_geom_020_1_sd + 33);

    auto ta2_zz_0_yz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 34);

    auto ta2_zz_0_zz_1 = pbuffer.data(idx_npot_geom_020_1_sd + 35);

    // Set up 0-6 components of targeted buffer : PD

    auto ta2_xx_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd);

    auto ta2_xx_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 1);

    auto ta2_xx_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 2);

    auto ta2_xx_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 3);

    auto ta2_xx_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 4);

    auto ta2_xx_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 5);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_0_xx_1, ta1_x_0_xy_1, ta1_x_0_xz_1, ta1_x_0_yy_1, ta1_x_0_yz_1, ta1_x_0_zz_1, ta2_xx_0_x_0, ta2_xx_0_x_1, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xy_0, ta2_xx_0_xy_1, ta2_xx_0_xz_0, ta2_xx_0_xz_1, ta2_xx_0_y_0, ta2_xx_0_y_1, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yz_0, ta2_xx_0_yz_1, ta2_xx_0_z_0, ta2_xx_0_z_1, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_x_xx_0, ta2_xx_x_xy_0, ta2_xx_x_xz_0, ta2_xx_x_yy_0, ta2_xx_x_yz_0, ta2_xx_x_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_x_xx_0[i] = 2.0 * ta2_xx_0_x_0[i] * fe_0 - 2.0 * ta2_xx_0_x_1[i] * fe_0 + 2.0 * ta1_x_0_xx_1[i] + ta2_xx_0_xx_0[i] * pa_x[i] - ta2_xx_0_xx_1[i] * pc_x[i];

        ta2_xx_x_xy_0[i] = ta2_xx_0_y_0[i] * fe_0 - ta2_xx_0_y_1[i] * fe_0 + 2.0 * ta1_x_0_xy_1[i] + ta2_xx_0_xy_0[i] * pa_x[i] - ta2_xx_0_xy_1[i] * pc_x[i];

        ta2_xx_x_xz_0[i] = ta2_xx_0_z_0[i] * fe_0 - ta2_xx_0_z_1[i] * fe_0 + 2.0 * ta1_x_0_xz_1[i] + ta2_xx_0_xz_0[i] * pa_x[i] - ta2_xx_0_xz_1[i] * pc_x[i];

        ta2_xx_x_yy_0[i] = 2.0 * ta1_x_0_yy_1[i] + ta2_xx_0_yy_0[i] * pa_x[i] - ta2_xx_0_yy_1[i] * pc_x[i];

        ta2_xx_x_yz_0[i] = 2.0 * ta1_x_0_yz_1[i] + ta2_xx_0_yz_0[i] * pa_x[i] - ta2_xx_0_yz_1[i] * pc_x[i];

        ta2_xx_x_zz_0[i] = 2.0 * ta1_x_0_zz_1[i] + ta2_xx_0_zz_0[i] * pa_x[i] - ta2_xx_0_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto ta2_xx_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 6);

    auto ta2_xx_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 7);

    auto ta2_xx_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 8);

    auto ta2_xx_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 9);

    auto ta2_xx_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 10);

    auto ta2_xx_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 11);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_0_x_0, ta2_xx_0_x_1, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xy_0, ta2_xx_0_xy_1, ta2_xx_0_xz_0, ta2_xx_0_xz_1, ta2_xx_0_y_0, ta2_xx_0_y_1, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yz_0, ta2_xx_0_yz_1, ta2_xx_0_z_0, ta2_xx_0_z_1, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_y_xx_0, ta2_xx_y_xy_0, ta2_xx_y_xz_0, ta2_xx_y_yy_0, ta2_xx_y_yz_0, ta2_xx_y_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_y_xx_0[i] = ta2_xx_0_xx_0[i] * pa_y[i] - ta2_xx_0_xx_1[i] * pc_y[i];

        ta2_xx_y_xy_0[i] = ta2_xx_0_x_0[i] * fe_0 - ta2_xx_0_x_1[i] * fe_0 + ta2_xx_0_xy_0[i] * pa_y[i] - ta2_xx_0_xy_1[i] * pc_y[i];

        ta2_xx_y_xz_0[i] = ta2_xx_0_xz_0[i] * pa_y[i] - ta2_xx_0_xz_1[i] * pc_y[i];

        ta2_xx_y_yy_0[i] = 2.0 * ta2_xx_0_y_0[i] * fe_0 - 2.0 * ta2_xx_0_y_1[i] * fe_0 + ta2_xx_0_yy_0[i] * pa_y[i] - ta2_xx_0_yy_1[i] * pc_y[i];

        ta2_xx_y_yz_0[i] = ta2_xx_0_z_0[i] * fe_0 - ta2_xx_0_z_1[i] * fe_0 + ta2_xx_0_yz_0[i] * pa_y[i] - ta2_xx_0_yz_1[i] * pc_y[i];

        ta2_xx_y_zz_0[i] = ta2_xx_0_zz_0[i] * pa_y[i] - ta2_xx_0_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto ta2_xx_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 12);

    auto ta2_xx_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 13);

    auto ta2_xx_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 14);

    auto ta2_xx_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 15);

    auto ta2_xx_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 16);

    auto ta2_xx_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 17);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_0_x_0, ta2_xx_0_x_1, ta2_xx_0_xx_0, ta2_xx_0_xx_1, ta2_xx_0_xy_0, ta2_xx_0_xy_1, ta2_xx_0_xz_0, ta2_xx_0_xz_1, ta2_xx_0_y_0, ta2_xx_0_y_1, ta2_xx_0_yy_0, ta2_xx_0_yy_1, ta2_xx_0_yz_0, ta2_xx_0_yz_1, ta2_xx_0_z_0, ta2_xx_0_z_1, ta2_xx_0_zz_0, ta2_xx_0_zz_1, ta2_xx_z_xx_0, ta2_xx_z_xy_0, ta2_xx_z_xz_0, ta2_xx_z_yy_0, ta2_xx_z_yz_0, ta2_xx_z_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_z_xx_0[i] = ta2_xx_0_xx_0[i] * pa_z[i] - ta2_xx_0_xx_1[i] * pc_z[i];

        ta2_xx_z_xy_0[i] = ta2_xx_0_xy_0[i] * pa_z[i] - ta2_xx_0_xy_1[i] * pc_z[i];

        ta2_xx_z_xz_0[i] = ta2_xx_0_x_0[i] * fe_0 - ta2_xx_0_x_1[i] * fe_0 + ta2_xx_0_xz_0[i] * pa_z[i] - ta2_xx_0_xz_1[i] * pc_z[i];

        ta2_xx_z_yy_0[i] = ta2_xx_0_yy_0[i] * pa_z[i] - ta2_xx_0_yy_1[i] * pc_z[i];

        ta2_xx_z_yz_0[i] = ta2_xx_0_y_0[i] * fe_0 - ta2_xx_0_y_1[i] * fe_0 + ta2_xx_0_yz_0[i] * pa_z[i] - ta2_xx_0_yz_1[i] * pc_z[i];

        ta2_xx_z_zz_0[i] = 2.0 * ta2_xx_0_z_0[i] * fe_0 - 2.0 * ta2_xx_0_z_1[i] * fe_0 + ta2_xx_0_zz_0[i] * pa_z[i] - ta2_xx_0_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : PD

    auto ta2_xy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 18);

    auto ta2_xy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 19);

    auto ta2_xy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 20);

    auto ta2_xy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 21);

    auto ta2_xy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 22);

    auto ta2_xy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 23);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_0_xx_1, ta1_y_0_xy_1, ta1_y_0_xz_1, ta1_y_0_yy_1, ta1_y_0_yz_1, ta1_y_0_zz_1, ta2_xy_0_x_0, ta2_xy_0_x_1, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xy_0, ta2_xy_0_xy_1, ta2_xy_0_xz_0, ta2_xy_0_xz_1, ta2_xy_0_y_0, ta2_xy_0_y_1, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yz_0, ta2_xy_0_yz_1, ta2_xy_0_z_0, ta2_xy_0_z_1, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_x_xx_0, ta2_xy_x_xy_0, ta2_xy_x_xz_0, ta2_xy_x_yy_0, ta2_xy_x_yz_0, ta2_xy_x_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_x_xx_0[i] = 2.0 * ta2_xy_0_x_0[i] * fe_0 - 2.0 * ta2_xy_0_x_1[i] * fe_0 + ta1_y_0_xx_1[i] + ta2_xy_0_xx_0[i] * pa_x[i] - ta2_xy_0_xx_1[i] * pc_x[i];

        ta2_xy_x_xy_0[i] = ta2_xy_0_y_0[i] * fe_0 - ta2_xy_0_y_1[i] * fe_0 + ta1_y_0_xy_1[i] + ta2_xy_0_xy_0[i] * pa_x[i] - ta2_xy_0_xy_1[i] * pc_x[i];

        ta2_xy_x_xz_0[i] = ta2_xy_0_z_0[i] * fe_0 - ta2_xy_0_z_1[i] * fe_0 + ta1_y_0_xz_1[i] + ta2_xy_0_xz_0[i] * pa_x[i] - ta2_xy_0_xz_1[i] * pc_x[i];

        ta2_xy_x_yy_0[i] = ta1_y_0_yy_1[i] + ta2_xy_0_yy_0[i] * pa_x[i] - ta2_xy_0_yy_1[i] * pc_x[i];

        ta2_xy_x_yz_0[i] = ta1_y_0_yz_1[i] + ta2_xy_0_yz_0[i] * pa_x[i] - ta2_xy_0_yz_1[i] * pc_x[i];

        ta2_xy_x_zz_0[i] = ta1_y_0_zz_1[i] + ta2_xy_0_zz_0[i] * pa_x[i] - ta2_xy_0_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : PD

    auto ta2_xy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 24);

    auto ta2_xy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 25);

    auto ta2_xy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 26);

    auto ta2_xy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 27);

    auto ta2_xy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 28);

    auto ta2_xy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 29);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_0_xx_1, ta1_x_0_xy_1, ta1_x_0_xz_1, ta1_x_0_yy_1, ta1_x_0_yz_1, ta1_x_0_zz_1, ta2_xy_0_x_0, ta2_xy_0_x_1, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xy_0, ta2_xy_0_xy_1, ta2_xy_0_xz_0, ta2_xy_0_xz_1, ta2_xy_0_y_0, ta2_xy_0_y_1, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yz_0, ta2_xy_0_yz_1, ta2_xy_0_z_0, ta2_xy_0_z_1, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_y_xx_0, ta2_xy_y_xy_0, ta2_xy_y_xz_0, ta2_xy_y_yy_0, ta2_xy_y_yz_0, ta2_xy_y_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_y_xx_0[i] = ta1_x_0_xx_1[i] + ta2_xy_0_xx_0[i] * pa_y[i] - ta2_xy_0_xx_1[i] * pc_y[i];

        ta2_xy_y_xy_0[i] = ta2_xy_0_x_0[i] * fe_0 - ta2_xy_0_x_1[i] * fe_0 + ta1_x_0_xy_1[i] + ta2_xy_0_xy_0[i] * pa_y[i] - ta2_xy_0_xy_1[i] * pc_y[i];

        ta2_xy_y_xz_0[i] = ta1_x_0_xz_1[i] + ta2_xy_0_xz_0[i] * pa_y[i] - ta2_xy_0_xz_1[i] * pc_y[i];

        ta2_xy_y_yy_0[i] = 2.0 * ta2_xy_0_y_0[i] * fe_0 - 2.0 * ta2_xy_0_y_1[i] * fe_0 + ta1_x_0_yy_1[i] + ta2_xy_0_yy_0[i] * pa_y[i] - ta2_xy_0_yy_1[i] * pc_y[i];

        ta2_xy_y_yz_0[i] = ta2_xy_0_z_0[i] * fe_0 - ta2_xy_0_z_1[i] * fe_0 + ta1_x_0_yz_1[i] + ta2_xy_0_yz_0[i] * pa_y[i] - ta2_xy_0_yz_1[i] * pc_y[i];

        ta2_xy_y_zz_0[i] = ta1_x_0_zz_1[i] + ta2_xy_0_zz_0[i] * pa_y[i] - ta2_xy_0_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : PD

    auto ta2_xy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 30);

    auto ta2_xy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 31);

    auto ta2_xy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 32);

    auto ta2_xy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 33);

    auto ta2_xy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 34);

    auto ta2_xy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 35);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_0_x_0, ta2_xy_0_x_1, ta2_xy_0_xx_0, ta2_xy_0_xx_1, ta2_xy_0_xy_0, ta2_xy_0_xy_1, ta2_xy_0_xz_0, ta2_xy_0_xz_1, ta2_xy_0_y_0, ta2_xy_0_y_1, ta2_xy_0_yy_0, ta2_xy_0_yy_1, ta2_xy_0_yz_0, ta2_xy_0_yz_1, ta2_xy_0_z_0, ta2_xy_0_z_1, ta2_xy_0_zz_0, ta2_xy_0_zz_1, ta2_xy_z_xx_0, ta2_xy_z_xy_0, ta2_xy_z_xz_0, ta2_xy_z_yy_0, ta2_xy_z_yz_0, ta2_xy_z_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_z_xx_0[i] = ta2_xy_0_xx_0[i] * pa_z[i] - ta2_xy_0_xx_1[i] * pc_z[i];

        ta2_xy_z_xy_0[i] = ta2_xy_0_xy_0[i] * pa_z[i] - ta2_xy_0_xy_1[i] * pc_z[i];

        ta2_xy_z_xz_0[i] = ta2_xy_0_x_0[i] * fe_0 - ta2_xy_0_x_1[i] * fe_0 + ta2_xy_0_xz_0[i] * pa_z[i] - ta2_xy_0_xz_1[i] * pc_z[i];

        ta2_xy_z_yy_0[i] = ta2_xy_0_yy_0[i] * pa_z[i] - ta2_xy_0_yy_1[i] * pc_z[i];

        ta2_xy_z_yz_0[i] = ta2_xy_0_y_0[i] * fe_0 - ta2_xy_0_y_1[i] * fe_0 + ta2_xy_0_yz_0[i] * pa_z[i] - ta2_xy_0_yz_1[i] * pc_z[i];

        ta2_xy_z_zz_0[i] = 2.0 * ta2_xy_0_z_0[i] * fe_0 - 2.0 * ta2_xy_0_z_1[i] * fe_0 + ta2_xy_0_zz_0[i] * pa_z[i] - ta2_xy_0_zz_1[i] * pc_z[i];
    }

    // Set up 36-42 components of targeted buffer : PD

    auto ta2_xz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 36);

    auto ta2_xz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 37);

    auto ta2_xz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 38);

    auto ta2_xz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 39);

    auto ta2_xz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 40);

    auto ta2_xz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 41);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_0_xx_1, ta1_z_0_xy_1, ta1_z_0_xz_1, ta1_z_0_yy_1, ta1_z_0_yz_1, ta1_z_0_zz_1, ta2_xz_0_x_0, ta2_xz_0_x_1, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xy_0, ta2_xz_0_xy_1, ta2_xz_0_xz_0, ta2_xz_0_xz_1, ta2_xz_0_y_0, ta2_xz_0_y_1, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yz_0, ta2_xz_0_yz_1, ta2_xz_0_z_0, ta2_xz_0_z_1, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_x_xx_0, ta2_xz_x_xy_0, ta2_xz_x_xz_0, ta2_xz_x_yy_0, ta2_xz_x_yz_0, ta2_xz_x_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_x_xx_0[i] = 2.0 * ta2_xz_0_x_0[i] * fe_0 - 2.0 * ta2_xz_0_x_1[i] * fe_0 + ta1_z_0_xx_1[i] + ta2_xz_0_xx_0[i] * pa_x[i] - ta2_xz_0_xx_1[i] * pc_x[i];

        ta2_xz_x_xy_0[i] = ta2_xz_0_y_0[i] * fe_0 - ta2_xz_0_y_1[i] * fe_0 + ta1_z_0_xy_1[i] + ta2_xz_0_xy_0[i] * pa_x[i] - ta2_xz_0_xy_1[i] * pc_x[i];

        ta2_xz_x_xz_0[i] = ta2_xz_0_z_0[i] * fe_0 - ta2_xz_0_z_1[i] * fe_0 + ta1_z_0_xz_1[i] + ta2_xz_0_xz_0[i] * pa_x[i] - ta2_xz_0_xz_1[i] * pc_x[i];

        ta2_xz_x_yy_0[i] = ta1_z_0_yy_1[i] + ta2_xz_0_yy_0[i] * pa_x[i] - ta2_xz_0_yy_1[i] * pc_x[i];

        ta2_xz_x_yz_0[i] = ta1_z_0_yz_1[i] + ta2_xz_0_yz_0[i] * pa_x[i] - ta2_xz_0_yz_1[i] * pc_x[i];

        ta2_xz_x_zz_0[i] = ta1_z_0_zz_1[i] + ta2_xz_0_zz_0[i] * pa_x[i] - ta2_xz_0_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : PD

    auto ta2_xz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 42);

    auto ta2_xz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 43);

    auto ta2_xz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 44);

    auto ta2_xz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 45);

    auto ta2_xz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 46);

    auto ta2_xz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 47);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_0_x_0, ta2_xz_0_x_1, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xy_0, ta2_xz_0_xy_1, ta2_xz_0_xz_0, ta2_xz_0_xz_1, ta2_xz_0_y_0, ta2_xz_0_y_1, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yz_0, ta2_xz_0_yz_1, ta2_xz_0_z_0, ta2_xz_0_z_1, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_y_xx_0, ta2_xz_y_xy_0, ta2_xz_y_xz_0, ta2_xz_y_yy_0, ta2_xz_y_yz_0, ta2_xz_y_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_y_xx_0[i] = ta2_xz_0_xx_0[i] * pa_y[i] - ta2_xz_0_xx_1[i] * pc_y[i];

        ta2_xz_y_xy_0[i] = ta2_xz_0_x_0[i] * fe_0 - ta2_xz_0_x_1[i] * fe_0 + ta2_xz_0_xy_0[i] * pa_y[i] - ta2_xz_0_xy_1[i] * pc_y[i];

        ta2_xz_y_xz_0[i] = ta2_xz_0_xz_0[i] * pa_y[i] - ta2_xz_0_xz_1[i] * pc_y[i];

        ta2_xz_y_yy_0[i] = 2.0 * ta2_xz_0_y_0[i] * fe_0 - 2.0 * ta2_xz_0_y_1[i] * fe_0 + ta2_xz_0_yy_0[i] * pa_y[i] - ta2_xz_0_yy_1[i] * pc_y[i];

        ta2_xz_y_yz_0[i] = ta2_xz_0_z_0[i] * fe_0 - ta2_xz_0_z_1[i] * fe_0 + ta2_xz_0_yz_0[i] * pa_y[i] - ta2_xz_0_yz_1[i] * pc_y[i];

        ta2_xz_y_zz_0[i] = ta2_xz_0_zz_0[i] * pa_y[i] - ta2_xz_0_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : PD

    auto ta2_xz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 48);

    auto ta2_xz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 49);

    auto ta2_xz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 50);

    auto ta2_xz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 51);

    auto ta2_xz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 52);

    auto ta2_xz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 53);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_0_xx_1, ta1_x_0_xy_1, ta1_x_0_xz_1, ta1_x_0_yy_1, ta1_x_0_yz_1, ta1_x_0_zz_1, ta2_xz_0_x_0, ta2_xz_0_x_1, ta2_xz_0_xx_0, ta2_xz_0_xx_1, ta2_xz_0_xy_0, ta2_xz_0_xy_1, ta2_xz_0_xz_0, ta2_xz_0_xz_1, ta2_xz_0_y_0, ta2_xz_0_y_1, ta2_xz_0_yy_0, ta2_xz_0_yy_1, ta2_xz_0_yz_0, ta2_xz_0_yz_1, ta2_xz_0_z_0, ta2_xz_0_z_1, ta2_xz_0_zz_0, ta2_xz_0_zz_1, ta2_xz_z_xx_0, ta2_xz_z_xy_0, ta2_xz_z_xz_0, ta2_xz_z_yy_0, ta2_xz_z_yz_0, ta2_xz_z_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_z_xx_0[i] = ta1_x_0_xx_1[i] + ta2_xz_0_xx_0[i] * pa_z[i] - ta2_xz_0_xx_1[i] * pc_z[i];

        ta2_xz_z_xy_0[i] = ta1_x_0_xy_1[i] + ta2_xz_0_xy_0[i] * pa_z[i] - ta2_xz_0_xy_1[i] * pc_z[i];

        ta2_xz_z_xz_0[i] = ta2_xz_0_x_0[i] * fe_0 - ta2_xz_0_x_1[i] * fe_0 + ta1_x_0_xz_1[i] + ta2_xz_0_xz_0[i] * pa_z[i] - ta2_xz_0_xz_1[i] * pc_z[i];

        ta2_xz_z_yy_0[i] = ta1_x_0_yy_1[i] + ta2_xz_0_yy_0[i] * pa_z[i] - ta2_xz_0_yy_1[i] * pc_z[i];

        ta2_xz_z_yz_0[i] = ta2_xz_0_y_0[i] * fe_0 - ta2_xz_0_y_1[i] * fe_0 + ta1_x_0_yz_1[i] + ta2_xz_0_yz_0[i] * pa_z[i] - ta2_xz_0_yz_1[i] * pc_z[i];

        ta2_xz_z_zz_0[i] = 2.0 * ta2_xz_0_z_0[i] * fe_0 - 2.0 * ta2_xz_0_z_1[i] * fe_0 + ta1_x_0_zz_1[i] + ta2_xz_0_zz_0[i] * pa_z[i] - ta2_xz_0_zz_1[i] * pc_z[i];
    }

    // Set up 54-60 components of targeted buffer : PD

    auto ta2_yy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 54);

    auto ta2_yy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 55);

    auto ta2_yy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 56);

    auto ta2_yy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 57);

    auto ta2_yy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 58);

    auto ta2_yy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 59);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_0_x_0, ta2_yy_0_x_1, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xy_0, ta2_yy_0_xy_1, ta2_yy_0_xz_0, ta2_yy_0_xz_1, ta2_yy_0_y_0, ta2_yy_0_y_1, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yz_0, ta2_yy_0_yz_1, ta2_yy_0_z_0, ta2_yy_0_z_1, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_x_xx_0, ta2_yy_x_xy_0, ta2_yy_x_xz_0, ta2_yy_x_yy_0, ta2_yy_x_yz_0, ta2_yy_x_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_x_xx_0[i] = 2.0 * ta2_yy_0_x_0[i] * fe_0 - 2.0 * ta2_yy_0_x_1[i] * fe_0 + ta2_yy_0_xx_0[i] * pa_x[i] - ta2_yy_0_xx_1[i] * pc_x[i];

        ta2_yy_x_xy_0[i] = ta2_yy_0_y_0[i] * fe_0 - ta2_yy_0_y_1[i] * fe_0 + ta2_yy_0_xy_0[i] * pa_x[i] - ta2_yy_0_xy_1[i] * pc_x[i];

        ta2_yy_x_xz_0[i] = ta2_yy_0_z_0[i] * fe_0 - ta2_yy_0_z_1[i] * fe_0 + ta2_yy_0_xz_0[i] * pa_x[i] - ta2_yy_0_xz_1[i] * pc_x[i];

        ta2_yy_x_yy_0[i] = ta2_yy_0_yy_0[i] * pa_x[i] - ta2_yy_0_yy_1[i] * pc_x[i];

        ta2_yy_x_yz_0[i] = ta2_yy_0_yz_0[i] * pa_x[i] - ta2_yy_0_yz_1[i] * pc_x[i];

        ta2_yy_x_zz_0[i] = ta2_yy_0_zz_0[i] * pa_x[i] - ta2_yy_0_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : PD

    auto ta2_yy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 60);

    auto ta2_yy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 61);

    auto ta2_yy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 62);

    auto ta2_yy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 63);

    auto ta2_yy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 64);

    auto ta2_yy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 65);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_0_xx_1, ta1_y_0_xy_1, ta1_y_0_xz_1, ta1_y_0_yy_1, ta1_y_0_yz_1, ta1_y_0_zz_1, ta2_yy_0_x_0, ta2_yy_0_x_1, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xy_0, ta2_yy_0_xy_1, ta2_yy_0_xz_0, ta2_yy_0_xz_1, ta2_yy_0_y_0, ta2_yy_0_y_1, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yz_0, ta2_yy_0_yz_1, ta2_yy_0_z_0, ta2_yy_0_z_1, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_y_xx_0, ta2_yy_y_xy_0, ta2_yy_y_xz_0, ta2_yy_y_yy_0, ta2_yy_y_yz_0, ta2_yy_y_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_y_xx_0[i] = 2.0 * ta1_y_0_xx_1[i] + ta2_yy_0_xx_0[i] * pa_y[i] - ta2_yy_0_xx_1[i] * pc_y[i];

        ta2_yy_y_xy_0[i] = ta2_yy_0_x_0[i] * fe_0 - ta2_yy_0_x_1[i] * fe_0 + 2.0 * ta1_y_0_xy_1[i] + ta2_yy_0_xy_0[i] * pa_y[i] - ta2_yy_0_xy_1[i] * pc_y[i];

        ta2_yy_y_xz_0[i] = 2.0 * ta1_y_0_xz_1[i] + ta2_yy_0_xz_0[i] * pa_y[i] - ta2_yy_0_xz_1[i] * pc_y[i];

        ta2_yy_y_yy_0[i] = 2.0 * ta2_yy_0_y_0[i] * fe_0 - 2.0 * ta2_yy_0_y_1[i] * fe_0 + 2.0 * ta1_y_0_yy_1[i] + ta2_yy_0_yy_0[i] * pa_y[i] - ta2_yy_0_yy_1[i] * pc_y[i];

        ta2_yy_y_yz_0[i] = ta2_yy_0_z_0[i] * fe_0 - ta2_yy_0_z_1[i] * fe_0 + 2.0 * ta1_y_0_yz_1[i] + ta2_yy_0_yz_0[i] * pa_y[i] - ta2_yy_0_yz_1[i] * pc_y[i];

        ta2_yy_y_zz_0[i] = 2.0 * ta1_y_0_zz_1[i] + ta2_yy_0_zz_0[i] * pa_y[i] - ta2_yy_0_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : PD

    auto ta2_yy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 66);

    auto ta2_yy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 67);

    auto ta2_yy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 68);

    auto ta2_yy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 69);

    auto ta2_yy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 70);

    auto ta2_yy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 71);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_0_x_0, ta2_yy_0_x_1, ta2_yy_0_xx_0, ta2_yy_0_xx_1, ta2_yy_0_xy_0, ta2_yy_0_xy_1, ta2_yy_0_xz_0, ta2_yy_0_xz_1, ta2_yy_0_y_0, ta2_yy_0_y_1, ta2_yy_0_yy_0, ta2_yy_0_yy_1, ta2_yy_0_yz_0, ta2_yy_0_yz_1, ta2_yy_0_z_0, ta2_yy_0_z_1, ta2_yy_0_zz_0, ta2_yy_0_zz_1, ta2_yy_z_xx_0, ta2_yy_z_xy_0, ta2_yy_z_xz_0, ta2_yy_z_yy_0, ta2_yy_z_yz_0, ta2_yy_z_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_z_xx_0[i] = ta2_yy_0_xx_0[i] * pa_z[i] - ta2_yy_0_xx_1[i] * pc_z[i];

        ta2_yy_z_xy_0[i] = ta2_yy_0_xy_0[i] * pa_z[i] - ta2_yy_0_xy_1[i] * pc_z[i];

        ta2_yy_z_xz_0[i] = ta2_yy_0_x_0[i] * fe_0 - ta2_yy_0_x_1[i] * fe_0 + ta2_yy_0_xz_0[i] * pa_z[i] - ta2_yy_0_xz_1[i] * pc_z[i];

        ta2_yy_z_yy_0[i] = ta2_yy_0_yy_0[i] * pa_z[i] - ta2_yy_0_yy_1[i] * pc_z[i];

        ta2_yy_z_yz_0[i] = ta2_yy_0_y_0[i] * fe_0 - ta2_yy_0_y_1[i] * fe_0 + ta2_yy_0_yz_0[i] * pa_z[i] - ta2_yy_0_yz_1[i] * pc_z[i];

        ta2_yy_z_zz_0[i] = 2.0 * ta2_yy_0_z_0[i] * fe_0 - 2.0 * ta2_yy_0_z_1[i] * fe_0 + ta2_yy_0_zz_0[i] * pa_z[i] - ta2_yy_0_zz_1[i] * pc_z[i];
    }

    // Set up 72-78 components of targeted buffer : PD

    auto ta2_yz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 72);

    auto ta2_yz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 73);

    auto ta2_yz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 74);

    auto ta2_yz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 75);

    auto ta2_yz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 76);

    auto ta2_yz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 77);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_0_x_0, ta2_yz_0_x_1, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xy_0, ta2_yz_0_xy_1, ta2_yz_0_xz_0, ta2_yz_0_xz_1, ta2_yz_0_y_0, ta2_yz_0_y_1, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yz_0, ta2_yz_0_yz_1, ta2_yz_0_z_0, ta2_yz_0_z_1, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_x_xx_0, ta2_yz_x_xy_0, ta2_yz_x_xz_0, ta2_yz_x_yy_0, ta2_yz_x_yz_0, ta2_yz_x_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_x_xx_0[i] = 2.0 * ta2_yz_0_x_0[i] * fe_0 - 2.0 * ta2_yz_0_x_1[i] * fe_0 + ta2_yz_0_xx_0[i] * pa_x[i] - ta2_yz_0_xx_1[i] * pc_x[i];

        ta2_yz_x_xy_0[i] = ta2_yz_0_y_0[i] * fe_0 - ta2_yz_0_y_1[i] * fe_0 + ta2_yz_0_xy_0[i] * pa_x[i] - ta2_yz_0_xy_1[i] * pc_x[i];

        ta2_yz_x_xz_0[i] = ta2_yz_0_z_0[i] * fe_0 - ta2_yz_0_z_1[i] * fe_0 + ta2_yz_0_xz_0[i] * pa_x[i] - ta2_yz_0_xz_1[i] * pc_x[i];

        ta2_yz_x_yy_0[i] = ta2_yz_0_yy_0[i] * pa_x[i] - ta2_yz_0_yy_1[i] * pc_x[i];

        ta2_yz_x_yz_0[i] = ta2_yz_0_yz_0[i] * pa_x[i] - ta2_yz_0_yz_1[i] * pc_x[i];

        ta2_yz_x_zz_0[i] = ta2_yz_0_zz_0[i] * pa_x[i] - ta2_yz_0_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : PD

    auto ta2_yz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 78);

    auto ta2_yz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 79);

    auto ta2_yz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 80);

    auto ta2_yz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 81);

    auto ta2_yz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 82);

    auto ta2_yz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 83);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_0_xx_1, ta1_z_0_xy_1, ta1_z_0_xz_1, ta1_z_0_yy_1, ta1_z_0_yz_1, ta1_z_0_zz_1, ta2_yz_0_x_0, ta2_yz_0_x_1, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xy_0, ta2_yz_0_xy_1, ta2_yz_0_xz_0, ta2_yz_0_xz_1, ta2_yz_0_y_0, ta2_yz_0_y_1, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yz_0, ta2_yz_0_yz_1, ta2_yz_0_z_0, ta2_yz_0_z_1, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_y_xx_0, ta2_yz_y_xy_0, ta2_yz_y_xz_0, ta2_yz_y_yy_0, ta2_yz_y_yz_0, ta2_yz_y_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_y_xx_0[i] = ta1_z_0_xx_1[i] + ta2_yz_0_xx_0[i] * pa_y[i] - ta2_yz_0_xx_1[i] * pc_y[i];

        ta2_yz_y_xy_0[i] = ta2_yz_0_x_0[i] * fe_0 - ta2_yz_0_x_1[i] * fe_0 + ta1_z_0_xy_1[i] + ta2_yz_0_xy_0[i] * pa_y[i] - ta2_yz_0_xy_1[i] * pc_y[i];

        ta2_yz_y_xz_0[i] = ta1_z_0_xz_1[i] + ta2_yz_0_xz_0[i] * pa_y[i] - ta2_yz_0_xz_1[i] * pc_y[i];

        ta2_yz_y_yy_0[i] = 2.0 * ta2_yz_0_y_0[i] * fe_0 - 2.0 * ta2_yz_0_y_1[i] * fe_0 + ta1_z_0_yy_1[i] + ta2_yz_0_yy_0[i] * pa_y[i] - ta2_yz_0_yy_1[i] * pc_y[i];

        ta2_yz_y_yz_0[i] = ta2_yz_0_z_0[i] * fe_0 - ta2_yz_0_z_1[i] * fe_0 + ta1_z_0_yz_1[i] + ta2_yz_0_yz_0[i] * pa_y[i] - ta2_yz_0_yz_1[i] * pc_y[i];

        ta2_yz_y_zz_0[i] = ta1_z_0_zz_1[i] + ta2_yz_0_zz_0[i] * pa_y[i] - ta2_yz_0_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : PD

    auto ta2_yz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 84);

    auto ta2_yz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 85);

    auto ta2_yz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 86);

    auto ta2_yz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 87);

    auto ta2_yz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 88);

    auto ta2_yz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 89);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_0_xx_1, ta1_y_0_xy_1, ta1_y_0_xz_1, ta1_y_0_yy_1, ta1_y_0_yz_1, ta1_y_0_zz_1, ta2_yz_0_x_0, ta2_yz_0_x_1, ta2_yz_0_xx_0, ta2_yz_0_xx_1, ta2_yz_0_xy_0, ta2_yz_0_xy_1, ta2_yz_0_xz_0, ta2_yz_0_xz_1, ta2_yz_0_y_0, ta2_yz_0_y_1, ta2_yz_0_yy_0, ta2_yz_0_yy_1, ta2_yz_0_yz_0, ta2_yz_0_yz_1, ta2_yz_0_z_0, ta2_yz_0_z_1, ta2_yz_0_zz_0, ta2_yz_0_zz_1, ta2_yz_z_xx_0, ta2_yz_z_xy_0, ta2_yz_z_xz_0, ta2_yz_z_yy_0, ta2_yz_z_yz_0, ta2_yz_z_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_z_xx_0[i] = ta1_y_0_xx_1[i] + ta2_yz_0_xx_0[i] * pa_z[i] - ta2_yz_0_xx_1[i] * pc_z[i];

        ta2_yz_z_xy_0[i] = ta1_y_0_xy_1[i] + ta2_yz_0_xy_0[i] * pa_z[i] - ta2_yz_0_xy_1[i] * pc_z[i];

        ta2_yz_z_xz_0[i] = ta2_yz_0_x_0[i] * fe_0 - ta2_yz_0_x_1[i] * fe_0 + ta1_y_0_xz_1[i] + ta2_yz_0_xz_0[i] * pa_z[i] - ta2_yz_0_xz_1[i] * pc_z[i];

        ta2_yz_z_yy_0[i] = ta1_y_0_yy_1[i] + ta2_yz_0_yy_0[i] * pa_z[i] - ta2_yz_0_yy_1[i] * pc_z[i];

        ta2_yz_z_yz_0[i] = ta2_yz_0_y_0[i] * fe_0 - ta2_yz_0_y_1[i] * fe_0 + ta1_y_0_yz_1[i] + ta2_yz_0_yz_0[i] * pa_z[i] - ta2_yz_0_yz_1[i] * pc_z[i];

        ta2_yz_z_zz_0[i] = 2.0 * ta2_yz_0_z_0[i] * fe_0 - 2.0 * ta2_yz_0_z_1[i] * fe_0 + ta1_y_0_zz_1[i] + ta2_yz_0_zz_0[i] * pa_z[i] - ta2_yz_0_zz_1[i] * pc_z[i];
    }

    // Set up 90-96 components of targeted buffer : PD

    auto ta2_zz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 90);

    auto ta2_zz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 91);

    auto ta2_zz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 92);

    auto ta2_zz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 93);

    auto ta2_zz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 94);

    auto ta2_zz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 95);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_0_x_0, ta2_zz_0_x_1, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xy_0, ta2_zz_0_xy_1, ta2_zz_0_xz_0, ta2_zz_0_xz_1, ta2_zz_0_y_0, ta2_zz_0_y_1, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yz_0, ta2_zz_0_yz_1, ta2_zz_0_z_0, ta2_zz_0_z_1, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_x_xx_0, ta2_zz_x_xy_0, ta2_zz_x_xz_0, ta2_zz_x_yy_0, ta2_zz_x_yz_0, ta2_zz_x_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_x_xx_0[i] = 2.0 * ta2_zz_0_x_0[i] * fe_0 - 2.0 * ta2_zz_0_x_1[i] * fe_0 + ta2_zz_0_xx_0[i] * pa_x[i] - ta2_zz_0_xx_1[i] * pc_x[i];

        ta2_zz_x_xy_0[i] = ta2_zz_0_y_0[i] * fe_0 - ta2_zz_0_y_1[i] * fe_0 + ta2_zz_0_xy_0[i] * pa_x[i] - ta2_zz_0_xy_1[i] * pc_x[i];

        ta2_zz_x_xz_0[i] = ta2_zz_0_z_0[i] * fe_0 - ta2_zz_0_z_1[i] * fe_0 + ta2_zz_0_xz_0[i] * pa_x[i] - ta2_zz_0_xz_1[i] * pc_x[i];

        ta2_zz_x_yy_0[i] = ta2_zz_0_yy_0[i] * pa_x[i] - ta2_zz_0_yy_1[i] * pc_x[i];

        ta2_zz_x_yz_0[i] = ta2_zz_0_yz_0[i] * pa_x[i] - ta2_zz_0_yz_1[i] * pc_x[i];

        ta2_zz_x_zz_0[i] = ta2_zz_0_zz_0[i] * pa_x[i] - ta2_zz_0_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : PD

    auto ta2_zz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 96);

    auto ta2_zz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 97);

    auto ta2_zz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 98);

    auto ta2_zz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 99);

    auto ta2_zz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 100);

    auto ta2_zz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 101);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_0_x_0, ta2_zz_0_x_1, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xy_0, ta2_zz_0_xy_1, ta2_zz_0_xz_0, ta2_zz_0_xz_1, ta2_zz_0_y_0, ta2_zz_0_y_1, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yz_0, ta2_zz_0_yz_1, ta2_zz_0_z_0, ta2_zz_0_z_1, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_y_xx_0, ta2_zz_y_xy_0, ta2_zz_y_xz_0, ta2_zz_y_yy_0, ta2_zz_y_yz_0, ta2_zz_y_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_y_xx_0[i] = ta2_zz_0_xx_0[i] * pa_y[i] - ta2_zz_0_xx_1[i] * pc_y[i];

        ta2_zz_y_xy_0[i] = ta2_zz_0_x_0[i] * fe_0 - ta2_zz_0_x_1[i] * fe_0 + ta2_zz_0_xy_0[i] * pa_y[i] - ta2_zz_0_xy_1[i] * pc_y[i];

        ta2_zz_y_xz_0[i] = ta2_zz_0_xz_0[i] * pa_y[i] - ta2_zz_0_xz_1[i] * pc_y[i];

        ta2_zz_y_yy_0[i] = 2.0 * ta2_zz_0_y_0[i] * fe_0 - 2.0 * ta2_zz_0_y_1[i] * fe_0 + ta2_zz_0_yy_0[i] * pa_y[i] - ta2_zz_0_yy_1[i] * pc_y[i];

        ta2_zz_y_yz_0[i] = ta2_zz_0_z_0[i] * fe_0 - ta2_zz_0_z_1[i] * fe_0 + ta2_zz_0_yz_0[i] * pa_y[i] - ta2_zz_0_yz_1[i] * pc_y[i];

        ta2_zz_y_zz_0[i] = ta2_zz_0_zz_0[i] * pa_y[i] - ta2_zz_0_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : PD

    auto ta2_zz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 102);

    auto ta2_zz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 103);

    auto ta2_zz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 104);

    auto ta2_zz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 105);

    auto ta2_zz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 106);

    auto ta2_zz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 107);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_0_xx_1, ta1_z_0_xy_1, ta1_z_0_xz_1, ta1_z_0_yy_1, ta1_z_0_yz_1, ta1_z_0_zz_1, ta2_zz_0_x_0, ta2_zz_0_x_1, ta2_zz_0_xx_0, ta2_zz_0_xx_1, ta2_zz_0_xy_0, ta2_zz_0_xy_1, ta2_zz_0_xz_0, ta2_zz_0_xz_1, ta2_zz_0_y_0, ta2_zz_0_y_1, ta2_zz_0_yy_0, ta2_zz_0_yy_1, ta2_zz_0_yz_0, ta2_zz_0_yz_1, ta2_zz_0_z_0, ta2_zz_0_z_1, ta2_zz_0_zz_0, ta2_zz_0_zz_1, ta2_zz_z_xx_0, ta2_zz_z_xy_0, ta2_zz_z_xz_0, ta2_zz_z_yy_0, ta2_zz_z_yz_0, ta2_zz_z_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_z_xx_0[i] = 2.0 * ta1_z_0_xx_1[i] + ta2_zz_0_xx_0[i] * pa_z[i] - ta2_zz_0_xx_1[i] * pc_z[i];

        ta2_zz_z_xy_0[i] = 2.0 * ta1_z_0_xy_1[i] + ta2_zz_0_xy_0[i] * pa_z[i] - ta2_zz_0_xy_1[i] * pc_z[i];

        ta2_zz_z_xz_0[i] = ta2_zz_0_x_0[i] * fe_0 - ta2_zz_0_x_1[i] * fe_0 + 2.0 * ta1_z_0_xz_1[i] + ta2_zz_0_xz_0[i] * pa_z[i] - ta2_zz_0_xz_1[i] * pc_z[i];

        ta2_zz_z_yy_0[i] = 2.0 * ta1_z_0_yy_1[i] + ta2_zz_0_yy_0[i] * pa_z[i] - ta2_zz_0_yy_1[i] * pc_z[i];

        ta2_zz_z_yz_0[i] = ta2_zz_0_y_0[i] * fe_0 - ta2_zz_0_y_1[i] * fe_0 + 2.0 * ta1_z_0_yz_1[i] + ta2_zz_0_yz_0[i] * pa_z[i] - ta2_zz_0_yz_1[i] * pc_z[i];

        ta2_zz_z_zz_0[i] = 2.0 * ta2_zz_0_z_0[i] * fe_0 - 2.0 * ta2_zz_0_z_1[i] * fe_0 + 2.0 * ta1_z_0_zz_1[i] + ta2_zz_0_zz_0[i] * pa_z[i] - ta2_zz_0_zz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

