#include "NuclearPotentialGridPrimRecSD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_sd(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_0_ss,
                                       const size_t idx_npot_1_ss,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const double rpb_x,
                                       const double rpb_y,
                                       const double rpb_z,
                                       const double factor) -> void
{
    // set up number of grid points

    const auto nelems = buffer.number_of_columns();

    // set up R(PC) = P - C distances

    auto pc_x = buffer.data();

    auto pc_y = &(buffer.data()[nelems]);

    auto pc_z = &(buffer.data()[2 * nelems]);

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_0 = &(buffer.data()[idx_npot_0_ss * nelems]);

    // Set up components of auxiliary buffer : SS

    auto ta_0_0_1 = &(buffer.data()[idx_npot_1_ss * nelems]);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_0 = &(buffer.data()[idx_npot_0_sp * nelems]);

    auto ta_0_y_0 = &(buffer.data()[(idx_npot_0_sp + 1) * nelems]);

    auto ta_0_z_0 = &(buffer.data()[(idx_npot_0_sp + 2) * nelems]);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_1 = &(buffer.data()[idx_npot_1_sp * nelems]);

    auto ta_0_y_1 = &(buffer.data()[(idx_npot_1_sp + 1) * nelems]);

    auto ta_0_z_1 = &(buffer.data()[(idx_npot_1_sp + 2) * nelems]);

    // Set up components of targeted buffer : SD

    auto ta_0_xx_0 = &(buffer.data()[idx_npot_0_sd * nelems]);

    auto ta_0_xy_0 = &(buffer.data()[(idx_npot_0_sd + 1) * nelems]);

    auto ta_0_xz_0 = &(buffer.data()[(idx_npot_0_sd + 2) * nelems]);

    auto ta_0_yy_0 = &(buffer.data()[(idx_npot_0_sd + 3) * nelems]);

    auto ta_0_yz_0 = &(buffer.data()[(idx_npot_0_sd + 4) * nelems]);

    auto ta_0_zz_0 = &(buffer.data()[(idx_npot_0_sd + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_0_xx_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_x_0[i] * rpb_x - ta_0_x_1[i] * pc_x[i];

        ta_0_xy_0[i] = ta_0_y_0[i] * rpb_x - ta_0_y_1[i] * pc_x[i];

        ta_0_xz_0[i] = ta_0_z_0[i] * rpb_x - ta_0_z_1[i] * pc_x[i];

        ta_0_yy_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_y_0[i] * rpb_y - ta_0_y_1[i] * pc_y[i];

        ta_0_yz_0[i] = ta_0_z_0[i] * rpb_y - ta_0_z_1[i] * pc_y[i];

        ta_0_zz_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_z_0[i] * rpb_z - ta_0_z_1[i] * pc_z[i];
    }
}

} // npotrec namespace

