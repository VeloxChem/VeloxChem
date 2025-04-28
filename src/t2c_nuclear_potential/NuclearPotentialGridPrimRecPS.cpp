#include "NuclearPotentialGridPrimRecPS.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_ps(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_ps,
                                       const size_t idx_npot_0_ss,
                                       const size_t idx_npot_1_ss,
                                       const double rpa_x,
                                       const double rpa_y,
                                       const double rpa_z) -> void
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

    // Set up components of targeted buffer : PS

    auto ta_x_0_0 = &(buffer.data()[idx_npot_0_ps * nelems]);

    auto ta_y_0_0 = &(buffer.data()[(idx_npot_0_ps + 1) * nelems]);

    auto ta_z_0_0 = &(buffer.data()[(idx_npot_0_ps + 2) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_x_0_0[i] = ta_0_0_0[i] * rpa_x - ta_0_0_1[i] * pc_x[i];

        ta_y_0_0[i] = ta_0_0_0[i] * rpa_y - ta_0_0_1[i] * pc_y[i];

        ta_z_0_0[i] = ta_0_0_0[i] * rpa_z - ta_0_0_1[i] * pc_z[i];
    }
}

} // npotrec namespace

