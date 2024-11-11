#include "NuclearPotentialGeom010PrimRecSP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_sp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_sp,
                                        const size_t              idx_npot_1_ss,
                                        const size_t              idx_npot_geom_010_0_ss,
                                        const size_t              idx_npot_geom_010_1_ss,
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

    auto ta_0_0_1 = pbuffer.data(idx_npot_1_ss);

    // Set up components of auxiliary buffer : SS

    auto ta1_x_0_0_0 = pbuffer.data(idx_npot_geom_010_0_ss);

    auto ta1_y_0_0_0 = pbuffer.data(idx_npot_geom_010_0_ss + 1);

    auto ta1_z_0_0_0 = pbuffer.data(idx_npot_geom_010_0_ss + 2);

    // Set up components of auxiliary buffer : SS

    auto ta1_x_0_0_1 = pbuffer.data(idx_npot_geom_010_1_ss);

    auto ta1_y_0_0_1 = pbuffer.data(idx_npot_geom_010_1_ss + 1);

    auto ta1_z_0_0_1 = pbuffer.data(idx_npot_geom_010_1_ss + 2);

    // Set up components of targeted buffer : SP

    auto ta1_x_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp);

    auto ta1_x_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 1);

    auto ta1_x_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 2);

    auto ta1_y_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp + 3);

    auto ta1_y_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 4);

    auto ta1_y_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 5);

    auto ta1_z_0_x_0 = pbuffer.data(idx_npot_geom_010_0_sp + 6);

    auto ta1_z_0_y_0 = pbuffer.data(idx_npot_geom_010_0_sp + 7);

    auto ta1_z_0_z_0 = pbuffer.data(idx_npot_geom_010_0_sp + 8);

#pragma omp simd aligned(pb_x,            \
                             pb_y,        \
                             pb_z,        \
                             pc_x,        \
                             pc_y,        \
                             pc_z,        \
                             ta1_x_0_0_0, \
                             ta1_x_0_0_1, \
                             ta1_x_0_x_0, \
                             ta1_x_0_y_0, \
                             ta1_x_0_z_0, \
                             ta1_y_0_0_0, \
                             ta1_y_0_0_1, \
                             ta1_y_0_x_0, \
                             ta1_y_0_y_0, \
                             ta1_y_0_z_0, \
                             ta1_z_0_0_0, \
                             ta1_z_0_0_1, \
                             ta1_z_0_x_0, \
                             ta1_z_0_y_0, \
                             ta1_z_0_z_0, \
                             ta_0_0_1 : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_0_x_0[i] = ta_0_0_1[i] + ta1_x_0_0_0[i] * pb_x[i] - ta1_x_0_0_1[i] * pc_x[i];

        ta1_x_0_y_0[i] = ta1_x_0_0_0[i] * pb_y[i] - ta1_x_0_0_1[i] * pc_y[i];

        ta1_x_0_z_0[i] = ta1_x_0_0_0[i] * pb_z[i] - ta1_x_0_0_1[i] * pc_z[i];

        ta1_y_0_x_0[i] = ta1_y_0_0_0[i] * pb_x[i] - ta1_y_0_0_1[i] * pc_x[i];

        ta1_y_0_y_0[i] = ta_0_0_1[i] + ta1_y_0_0_0[i] * pb_y[i] - ta1_y_0_0_1[i] * pc_y[i];

        ta1_y_0_z_0[i] = ta1_y_0_0_0[i] * pb_z[i] - ta1_y_0_0_1[i] * pc_z[i];

        ta1_z_0_x_0[i] = ta1_z_0_0_0[i] * pb_x[i] - ta1_z_0_0_1[i] * pc_x[i];

        ta1_z_0_y_0[i] = ta1_z_0_0_0[i] * pb_y[i] - ta1_z_0_0_1[i] * pc_y[i];

        ta1_z_0_z_0[i] = ta_0_0_1[i] + ta1_z_0_0_0[i] * pb_z[i] - ta1_z_0_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
