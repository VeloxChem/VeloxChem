#include "NuclearPotentialGeom020PrimRecPS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_ps(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_ps,
                                        const size_t              idx_npot_geom_010_1_ss,
                                        const size_t              idx_npot_geom_020_0_ss,
                                        const size_t              idx_npot_geom_020_1_ss,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SS

    auto ta1_x_0_0_1 = pbuffer.data(idx_npot_geom_010_1_ss);

    auto ta1_y_0_0_1 = pbuffer.data(idx_npot_geom_010_1_ss + 1);

    auto ta1_z_0_0_1 = pbuffer.data(idx_npot_geom_010_1_ss + 2);

    // Set up components of auxiliary buffer : SS

    auto ta2_xx_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss);

    auto ta2_xy_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 1);

    auto ta2_xz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 2);

    auto ta2_yy_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 3);

    auto ta2_yz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 4);

    auto ta2_zz_0_0_0 = pbuffer.data(idx_npot_geom_020_0_ss + 5);

    // Set up components of auxiliary buffer : SS

    auto ta2_xx_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss);

    auto ta2_xy_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 1);

    auto ta2_xz_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 2);

    auto ta2_yy_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 3);

    auto ta2_yz_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 4);

    auto ta2_zz_0_0_1 = pbuffer.data(idx_npot_geom_020_1_ss + 5);

    // Set up components of targeted buffer : PS

    auto ta2_xx_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps);

    auto ta2_xx_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 1);

    auto ta2_xx_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 2);

    auto ta2_xy_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 3);

    auto ta2_xy_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 4);

    auto ta2_xy_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 5);

    auto ta2_xz_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 6);

    auto ta2_xz_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 7);

    auto ta2_xz_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 8);

    auto ta2_yy_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 9);

    auto ta2_yy_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 10);

    auto ta2_yy_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 11);

    auto ta2_yz_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 12);

    auto ta2_yz_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 13);

    auto ta2_yz_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 14);

    auto ta2_zz_x_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 15);

    auto ta2_zz_y_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 16);

    auto ta2_zz_z_0_0 = pbuffer.data(idx_npot_geom_020_0_ps + 17);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             pc_x,         \
                             pc_y,         \
                             pc_z,         \
                             ta1_x_0_0_1,  \
                             ta1_y_0_0_1,  \
                             ta1_z_0_0_1,  \
                             ta2_xx_0_0_0, \
                             ta2_xx_0_0_1, \
                             ta2_xx_x_0_0, \
                             ta2_xx_y_0_0, \
                             ta2_xx_z_0_0, \
                             ta2_xy_0_0_0, \
                             ta2_xy_0_0_1, \
                             ta2_xy_x_0_0, \
                             ta2_xy_y_0_0, \
                             ta2_xy_z_0_0, \
                             ta2_xz_0_0_0, \
                             ta2_xz_0_0_1, \
                             ta2_xz_x_0_0, \
                             ta2_xz_y_0_0, \
                             ta2_xz_z_0_0, \
                             ta2_yy_0_0_0, \
                             ta2_yy_0_0_1, \
                             ta2_yy_x_0_0, \
                             ta2_yy_y_0_0, \
                             ta2_yy_z_0_0, \
                             ta2_yz_0_0_0, \
                             ta2_yz_0_0_1, \
                             ta2_yz_x_0_0, \
                             ta2_yz_y_0_0, \
                             ta2_yz_z_0_0, \
                             ta2_zz_0_0_0, \
                             ta2_zz_0_0_1, \
                             ta2_zz_x_0_0, \
                             ta2_zz_y_0_0, \
                             ta2_zz_z_0_0 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_x_0_0[i] = 2.0 * ta1_x_0_0_1[i] + ta2_xx_0_0_0[i] * pa_x[i] - ta2_xx_0_0_1[i] * pc_x[i];

        ta2_xx_y_0_0[i] = ta2_xx_0_0_0[i] * pa_y[i] - ta2_xx_0_0_1[i] * pc_y[i];

        ta2_xx_z_0_0[i] = ta2_xx_0_0_0[i] * pa_z[i] - ta2_xx_0_0_1[i] * pc_z[i];

        ta2_xy_x_0_0[i] = ta1_y_0_0_1[i] + ta2_xy_0_0_0[i] * pa_x[i] - ta2_xy_0_0_1[i] * pc_x[i];

        ta2_xy_y_0_0[i] = ta1_x_0_0_1[i] + ta2_xy_0_0_0[i] * pa_y[i] - ta2_xy_0_0_1[i] * pc_y[i];

        ta2_xy_z_0_0[i] = ta2_xy_0_0_0[i] * pa_z[i] - ta2_xy_0_0_1[i] * pc_z[i];

        ta2_xz_x_0_0[i] = ta1_z_0_0_1[i] + ta2_xz_0_0_0[i] * pa_x[i] - ta2_xz_0_0_1[i] * pc_x[i];

        ta2_xz_y_0_0[i] = ta2_xz_0_0_0[i] * pa_y[i] - ta2_xz_0_0_1[i] * pc_y[i];

        ta2_xz_z_0_0[i] = ta1_x_0_0_1[i] + ta2_xz_0_0_0[i] * pa_z[i] - ta2_xz_0_0_1[i] * pc_z[i];

        ta2_yy_x_0_0[i] = ta2_yy_0_0_0[i] * pa_x[i] - ta2_yy_0_0_1[i] * pc_x[i];

        ta2_yy_y_0_0[i] = 2.0 * ta1_y_0_0_1[i] + ta2_yy_0_0_0[i] * pa_y[i] - ta2_yy_0_0_1[i] * pc_y[i];

        ta2_yy_z_0_0[i] = ta2_yy_0_0_0[i] * pa_z[i] - ta2_yy_0_0_1[i] * pc_z[i];

        ta2_yz_x_0_0[i] = ta2_yz_0_0_0[i] * pa_x[i] - ta2_yz_0_0_1[i] * pc_x[i];

        ta2_yz_y_0_0[i] = ta1_z_0_0_1[i] + ta2_yz_0_0_0[i] * pa_y[i] - ta2_yz_0_0_1[i] * pc_y[i];

        ta2_yz_z_0_0[i] = ta1_y_0_0_1[i] + ta2_yz_0_0_0[i] * pa_z[i] - ta2_yz_0_0_1[i] * pc_z[i];

        ta2_zz_x_0_0[i] = ta2_zz_0_0_0[i] * pa_x[i] - ta2_zz_0_0_1[i] * pc_x[i];

        ta2_zz_y_0_0[i] = ta2_zz_0_0_0[i] * pa_y[i] - ta2_zz_0_0_1[i] * pc_y[i];

        ta2_zz_z_0_0[i] = 2.0 * ta1_z_0_0_1[i] + ta2_zz_0_0_0[i] * pa_z[i] - ta2_zz_0_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
