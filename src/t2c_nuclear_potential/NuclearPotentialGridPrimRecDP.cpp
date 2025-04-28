#include "NuclearPotentialGridPrimRecDP.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_dp(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_dp,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const size_t idx_npot_0_ps,
                                       const size_t idx_npot_1_ps,
                                       const size_t idx_npot_0_pp,
                                       const size_t idx_npot_1_pp,
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

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_0 = &(buffer.data()[idx_npot_0_sp * nelems]);

    auto ta_0_y_0 = &(buffer.data()[(idx_npot_0_sp + 1) * nelems]);

    auto ta_0_z_0 = &(buffer.data()[(idx_npot_0_sp + 2) * nelems]);

    // Set up components of auxiliary buffer : SP

    auto ta_0_x_1 = &(buffer.data()[idx_npot_1_sp * nelems]);

    auto ta_0_y_1 = &(buffer.data()[(idx_npot_1_sp + 1) * nelems]);

    auto ta_0_z_1 = &(buffer.data()[(idx_npot_1_sp + 2) * nelems]);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_0 = &(buffer.data()[idx_npot_0_ps * nelems]);

    auto ta_y_0_0 = &(buffer.data()[(idx_npot_0_ps + 1) * nelems]);

    auto ta_z_0_0 = &(buffer.data()[(idx_npot_0_ps + 2) * nelems]);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_1 = &(buffer.data()[idx_npot_1_ps * nelems]);

    auto ta_y_0_1 = &(buffer.data()[(idx_npot_1_ps + 1) * nelems]);

    auto ta_z_0_1 = &(buffer.data()[(idx_npot_1_ps + 2) * nelems]);

    // Set up components of auxiliary buffer : PP

    auto ta_x_x_0 = &(buffer.data()[idx_npot_0_pp * nelems]);

    auto ta_x_y_0 = &(buffer.data()[(idx_npot_0_pp + 1) * nelems]);

    auto ta_x_z_0 = &(buffer.data()[(idx_npot_0_pp + 2) * nelems]);

    auto ta_y_x_0 = &(buffer.data()[(idx_npot_0_pp + 3) * nelems]);

    auto ta_y_y_0 = &(buffer.data()[(idx_npot_0_pp + 4) * nelems]);

    auto ta_y_z_0 = &(buffer.data()[(idx_npot_0_pp + 5) * nelems]);

    auto ta_z_x_0 = &(buffer.data()[(idx_npot_0_pp + 6) * nelems]);

    auto ta_z_y_0 = &(buffer.data()[(idx_npot_0_pp + 7) * nelems]);

    auto ta_z_z_0 = &(buffer.data()[(idx_npot_0_pp + 8) * nelems]);

    // Set up components of auxiliary buffer : PP

    auto ta_x_x_1 = &(buffer.data()[idx_npot_1_pp * nelems]);

    auto ta_x_y_1 = &(buffer.data()[(idx_npot_1_pp + 1) * nelems]);

    auto ta_x_z_1 = &(buffer.data()[(idx_npot_1_pp + 2) * nelems]);

    auto ta_y_x_1 = &(buffer.data()[(idx_npot_1_pp + 3) * nelems]);

    auto ta_y_y_1 = &(buffer.data()[(idx_npot_1_pp + 4) * nelems]);

    auto ta_y_z_1 = &(buffer.data()[(idx_npot_1_pp + 5) * nelems]);

    auto ta_z_x_1 = &(buffer.data()[(idx_npot_1_pp + 6) * nelems]);

    auto ta_z_y_1 = &(buffer.data()[(idx_npot_1_pp + 7) * nelems]);

    auto ta_z_z_1 = &(buffer.data()[(idx_npot_1_pp + 8) * nelems]);

    // Set up 0-3 components of targeted buffer : DP

    auto ta_xx_x_0 = &(buffer.data()[idx_npot_0_dp * nelems]);

    auto ta_xx_y_0 = &(buffer.data()[(idx_npot_0_dp + 1) * nelems]);

    auto ta_xx_z_0 = &(buffer.data()[(idx_npot_0_dp + 2) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xx_x_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_x_0_0[i] * fe_0 - ta_x_0_1[i] * fe_0 + ta_x_x_0[i] * rpa_x - ta_x_x_1[i] * pc_x[i];

        ta_xx_y_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_x_y_0[i] * rpa_x - ta_x_y_1[i] * pc_x[i];

        ta_xx_z_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_x_z_0[i] * rpa_x - ta_x_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto ta_xy_x_0 = &(buffer.data()[(idx_npot_0_dp + 3) * nelems]);

    auto ta_xy_y_0 = &(buffer.data()[(idx_npot_0_dp + 4) * nelems]);

    auto ta_xy_z_0 = &(buffer.data()[(idx_npot_0_dp + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xy_x_0[i] = ta_x_x_0[i] * rpa_y - ta_x_x_1[i] * pc_y[i];

        ta_xy_y_0[i] = ta_y_y_0[i] * rpa_x - ta_y_y_1[i] * pc_x[i];

        ta_xy_z_0[i] = ta_y_z_0[i] * rpa_x - ta_y_z_1[i] * pc_x[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto ta_xz_x_0 = &(buffer.data()[(idx_npot_0_dp + 6) * nelems]);

    auto ta_xz_y_0 = &(buffer.data()[(idx_npot_0_dp + 7) * nelems]);

    auto ta_xz_z_0 = &(buffer.data()[(idx_npot_0_dp + 8) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xz_x_0[i] = ta_x_x_0[i] * rpa_z - ta_x_x_1[i] * pc_z[i];

        ta_xz_y_0[i] = ta_z_y_0[i] * rpa_x - ta_z_y_1[i] * pc_x[i];

        ta_xz_z_0[i] = ta_z_z_0[i] * rpa_x - ta_z_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto ta_yy_x_0 = &(buffer.data()[(idx_npot_0_dp + 9) * nelems]);

    auto ta_yy_y_0 = &(buffer.data()[(idx_npot_0_dp + 10) * nelems]);

    auto ta_yy_z_0 = &(buffer.data()[(idx_npot_0_dp + 11) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yy_x_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_y_x_0[i] * rpa_y - ta_y_x_1[i] * pc_y[i];

        ta_yy_y_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_y_0_0[i] * fe_0 - ta_y_0_1[i] * fe_0 + ta_y_y_0[i] * rpa_y - ta_y_y_1[i] * pc_y[i];

        ta_yy_z_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_y_z_0[i] * rpa_y - ta_y_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto ta_yz_x_0 = &(buffer.data()[(idx_npot_0_dp + 12) * nelems]);

    auto ta_yz_y_0 = &(buffer.data()[(idx_npot_0_dp + 13) * nelems]);

    auto ta_yz_z_0 = &(buffer.data()[(idx_npot_0_dp + 14) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_yz_x_0[i] = ta_z_x_0[i] * rpa_y - ta_z_x_1[i] * pc_y[i];

        ta_yz_y_0[i] = ta_y_y_0[i] * rpa_z - ta_y_y_1[i] * pc_z[i];

        ta_yz_z_0[i] = ta_z_z_0[i] * rpa_y - ta_z_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto ta_zz_x_0 = &(buffer.data()[(idx_npot_0_dp + 15) * nelems]);

    auto ta_zz_y_0 = &(buffer.data()[(idx_npot_0_dp + 16) * nelems]);

    auto ta_zz_z_0 = &(buffer.data()[(idx_npot_0_dp + 17) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_zz_x_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_z_x_0[i] * rpa_z - ta_z_x_1[i] * pc_z[i];

        ta_zz_y_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_z_y_0[i] * rpa_z - ta_z_y_1[i] * pc_z[i];

        ta_zz_z_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_z_0_0[i] * fe_0 - ta_z_0_1[i] * fe_0 + ta_z_z_0[i] * rpa_z - ta_z_z_1[i] * pc_z[i];
    }

}

} // npotrec namespace

