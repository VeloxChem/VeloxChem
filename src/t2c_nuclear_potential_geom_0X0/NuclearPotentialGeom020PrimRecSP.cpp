#include "NuclearPotentialGeom020PrimRecSP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_sp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_sp,
                                        const size_t              idx_npot_geom_010_1_ss,
                                        const size_t              idx_npot_geom_020_0_ss,
                                        const size_t              idx_npot_geom_020_1_ss,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpb,
                                        const size_t              idx_rpc) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

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

    // Set up components of targeted buffer : SP

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

#pragma omp simd aligned(pb_x,             \
                             pb_y,         \
                             pb_z,         \
                             pc_x,         \
                             pc_y,         \
                             pc_z,         \
                             ta1_x_0_0_1,  \
                             ta1_y_0_0_1,  \
                             ta1_z_0_0_1,  \
                             ta2_xx_0_0_0, \
                             ta2_xx_0_0_1, \
                             ta2_xx_0_x_0, \
                             ta2_xx_0_y_0, \
                             ta2_xx_0_z_0, \
                             ta2_xy_0_0_0, \
                             ta2_xy_0_0_1, \
                             ta2_xy_0_x_0, \
                             ta2_xy_0_y_0, \
                             ta2_xy_0_z_0, \
                             ta2_xz_0_0_0, \
                             ta2_xz_0_0_1, \
                             ta2_xz_0_x_0, \
                             ta2_xz_0_y_0, \
                             ta2_xz_0_z_0, \
                             ta2_yy_0_0_0, \
                             ta2_yy_0_0_1, \
                             ta2_yy_0_x_0, \
                             ta2_yy_0_y_0, \
                             ta2_yy_0_z_0, \
                             ta2_yz_0_0_0, \
                             ta2_yz_0_0_1, \
                             ta2_yz_0_x_0, \
                             ta2_yz_0_y_0, \
                             ta2_yz_0_z_0, \
                             ta2_zz_0_0_0, \
                             ta2_zz_0_0_1, \
                             ta2_zz_0_x_0, \
                             ta2_zz_0_y_0, \
                             ta2_zz_0_z_0 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_0_x_0[i] = 2.0 * ta1_x_0_0_1[i] + ta2_xx_0_0_0[i] * pb_x[i] - ta2_xx_0_0_1[i] * pc_x[i];

        ta2_xx_0_y_0[i] = ta2_xx_0_0_0[i] * pb_y[i] - ta2_xx_0_0_1[i] * pc_y[i];

        ta2_xx_0_z_0[i] = ta2_xx_0_0_0[i] * pb_z[i] - ta2_xx_0_0_1[i] * pc_z[i];

        ta2_xy_0_x_0[i] = ta1_y_0_0_1[i] + ta2_xy_0_0_0[i] * pb_x[i] - ta2_xy_0_0_1[i] * pc_x[i];

        ta2_xy_0_y_0[i] = ta1_x_0_0_1[i] + ta2_xy_0_0_0[i] * pb_y[i] - ta2_xy_0_0_1[i] * pc_y[i];

        ta2_xy_0_z_0[i] = ta2_xy_0_0_0[i] * pb_z[i] - ta2_xy_0_0_1[i] * pc_z[i];

        ta2_xz_0_x_0[i] = ta1_z_0_0_1[i] + ta2_xz_0_0_0[i] * pb_x[i] - ta2_xz_0_0_1[i] * pc_x[i];

        ta2_xz_0_y_0[i] = ta2_xz_0_0_0[i] * pb_y[i] - ta2_xz_0_0_1[i] * pc_y[i];

        ta2_xz_0_z_0[i] = ta1_x_0_0_1[i] + ta2_xz_0_0_0[i] * pb_z[i] - ta2_xz_0_0_1[i] * pc_z[i];

        ta2_yy_0_x_0[i] = ta2_yy_0_0_0[i] * pb_x[i] - ta2_yy_0_0_1[i] * pc_x[i];

        ta2_yy_0_y_0[i] = 2.0 * ta1_y_0_0_1[i] + ta2_yy_0_0_0[i] * pb_y[i] - ta2_yy_0_0_1[i] * pc_y[i];

        ta2_yy_0_z_0[i] = ta2_yy_0_0_0[i] * pb_z[i] - ta2_yy_0_0_1[i] * pc_z[i];

        ta2_yz_0_x_0[i] = ta2_yz_0_0_0[i] * pb_x[i] - ta2_yz_0_0_1[i] * pc_x[i];

        ta2_yz_0_y_0[i] = ta1_z_0_0_1[i] + ta2_yz_0_0_0[i] * pb_y[i] - ta2_yz_0_0_1[i] * pc_y[i];

        ta2_yz_0_z_0[i] = ta1_y_0_0_1[i] + ta2_yz_0_0_0[i] * pb_z[i] - ta2_yz_0_0_1[i] * pc_z[i];

        ta2_zz_0_x_0[i] = ta2_zz_0_0_0[i] * pb_x[i] - ta2_zz_0_0_1[i] * pc_x[i];

        ta2_zz_0_y_0[i] = ta2_zz_0_0_0[i] * pb_y[i] - ta2_zz_0_0_1[i] * pc_y[i];

        ta2_zz_0_z_0[i] = 2.0 * ta1_z_0_0_1[i] + ta2_zz_0_0_0[i] * pb_z[i] - ta2_zz_0_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
