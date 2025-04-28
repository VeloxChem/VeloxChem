#include "NuclearPotentialGridPrimRecDS.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_ds(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_ds,
                                       const size_t idx_npot_0_ss,
                                       const size_t idx_npot_1_ss,
                                       const size_t idx_npot_0_ps,
                                       const size_t idx_npot_1_ps,
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

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_0 = &(buffer.data()[idx_npot_0_ps * nelems]);

    auto ta_y_0_0 = &(buffer.data()[(idx_npot_0_ps + 1) * nelems]);

    auto ta_z_0_0 = &(buffer.data()[(idx_npot_0_ps + 2) * nelems]);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_1 = &(buffer.data()[idx_npot_1_ps * nelems]);

    auto ta_y_0_1 = &(buffer.data()[(idx_npot_1_ps + 1) * nelems]);

    auto ta_z_0_1 = &(buffer.data()[(idx_npot_1_ps + 2) * nelems]);

    // Set up components of targeted buffer : DS

    auto ta_xx_0_0 = &(buffer.data()[idx_npot_0_ds * nelems]);

    auto ta_xy_0_0 = &(buffer.data()[(idx_npot_0_ds + 1) * nelems]);

    auto ta_xz_0_0 = &(buffer.data()[(idx_npot_0_ds + 2) * nelems]);

    auto ta_yy_0_0 = &(buffer.data()[(idx_npot_0_ds + 3) * nelems]);

    auto ta_yz_0_0 = &(buffer.data()[(idx_npot_0_ds + 4) * nelems]);

    auto ta_zz_0_0 = &(buffer.data()[(idx_npot_0_ds + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xx_0_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_x_0_0[i] * rpa_x - ta_x_0_1[i] * pc_x[i];

        ta_xy_0_0[i] = ta_y_0_0[i] * rpa_x - ta_y_0_1[i] * pc_x[i];

        ta_xz_0_0[i] = ta_z_0_0[i] * rpa_x - ta_z_0_1[i] * pc_x[i];

        ta_yy_0_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_y_0_0[i] * rpa_y - ta_y_0_1[i] * pc_y[i];

        ta_yz_0_0[i] = ta_z_0_0[i] * rpa_y - ta_z_0_1[i] * pc_y[i];

        ta_zz_0_0[i] = ta_0_0_0[i] * fe_0 - ta_0_0_1[i] * fe_0 + ta_z_0_0[i] * rpa_z - ta_z_0_1[i] * pc_z[i];
    }
}

} // npotrec namespace

