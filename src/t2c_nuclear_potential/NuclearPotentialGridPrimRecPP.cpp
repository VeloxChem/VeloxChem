#include "NuclearPotentialGridPrimRecPP.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_pp(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_pp,
                                       const size_t idx_npot_0_ss,
                                       const size_t idx_npot_1_ss,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z,
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

    // Set up 0-3 components of targeted buffer : PP

    auto ta_x_x_0 = &(buffer.data()[idx_npot_0_pp * nelems]);

    auto ta_x_y_0 = &(buffer.data()[(idx_npot_0_pp + 1) * nelems]);

    auto ta_x_z_0 = &(buffer.data()[(idx_npot_0_pp + 2) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_x_x_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_x_0[i] * rpa_x - ta_0_x_1[i] * pc_x[i];

        ta_x_y_0[i] = ta_0_y_0[i] * rpa_x - ta_0_y_1[i] * pc_x[i];

        ta_x_z_0[i] = ta_0_z_0[i] * rpa_x - ta_0_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto ta_y_x_0 = &(buffer.data()[(idx_npot_0_pp + 3) * nelems]);

    auto ta_y_y_0 = &(buffer.data()[(idx_npot_0_pp + 4) * nelems]);

    auto ta_y_z_0 = &(buffer.data()[(idx_npot_0_pp + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_y_x_0[i] = ta_0_x_0[i] * rpa_y - ta_0_x_1[i] * pc_y[i];

        ta_y_y_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_y_0[i] * rpa_y - ta_0_y_1[i] * pc_y[i];

        ta_y_z_0[i] = ta_0_z_0[i] * rpa_y - ta_0_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : PP

    auto ta_z_x_0 = &(buffer.data()[(idx_npot_0_pp + 6) * nelems]);

    auto ta_z_y_0 = &(buffer.data()[(idx_npot_0_pp + 7) * nelems]);

    auto ta_z_z_0 = &(buffer.data()[(idx_npot_0_pp + 8) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_z_x_0[i] = ta_0_x_0[i] * rpa_z - ta_0_x_1[i] * pc_z[i];

        ta_z_y_0[i] = ta_0_y_0[i] * rpa_z - ta_0_y_1[i] * pc_z[i];

        ta_z_z_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_0_z_0[i] * rpa_z - ta_0_z_1[i] * pc_z[i];
    }

}

} // npotrec namespace

