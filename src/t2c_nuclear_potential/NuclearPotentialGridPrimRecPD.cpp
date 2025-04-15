#include "NuclearPotentialGridPrimRecPD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_pd(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_1_sd,
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

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_0 = &(buffer.data()[idx_npot_0_sd * nelems]);

    auto ta_0_xy_0 = &(buffer.data()[(idx_npot_0_sd + 1) * nelems]);

    auto ta_0_xz_0 = &(buffer.data()[(idx_npot_0_sd + 2) * nelems]);

    auto ta_0_yy_0 = &(buffer.data()[(idx_npot_0_sd + 3) * nelems]);

    auto ta_0_yz_0 = &(buffer.data()[(idx_npot_0_sd + 4) * nelems]);

    auto ta_0_zz_0 = &(buffer.data()[(idx_npot_0_sd + 5) * nelems]);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = &(buffer.data()[idx_npot_1_sd * nelems]);

    auto ta_0_xy_1 = &(buffer.data()[(idx_npot_1_sd + 1) * nelems]);

    auto ta_0_xz_1 = &(buffer.data()[(idx_npot_1_sd + 2) * nelems]);

    auto ta_0_yy_1 = &(buffer.data()[(idx_npot_1_sd + 3) * nelems]);

    auto ta_0_yz_1 = &(buffer.data()[(idx_npot_1_sd + 4) * nelems]);

    auto ta_0_zz_1 = &(buffer.data()[(idx_npot_1_sd + 5) * nelems]);

    // Set up 0-6 components of targeted buffer : PD

    auto ta_x_xx_0 = &(buffer.data()[idx_npot_0_pd * nelems]);

    auto ta_x_xy_0 = &(buffer.data()[(idx_npot_0_pd + 1) * nelems]);

    auto ta_x_xz_0 = &(buffer.data()[(idx_npot_0_pd + 2) * nelems]);

    auto ta_x_yy_0 = &(buffer.data()[(idx_npot_0_pd + 3) * nelems]);

    auto ta_x_yz_0 = &(buffer.data()[(idx_npot_0_pd + 4) * nelems]);

    auto ta_x_zz_0 = &(buffer.data()[(idx_npot_0_pd + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_x_xx_0[i] = 2.0 * ta_0_x_0[i] * fe_0 - 2.0 * ta_0_x_1[i] * fe_0 + ta_0_xx_0[i] * rpa_x - ta_0_xx_1[i] * pc_x[i];

        ta_x_xy_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_0_xy_0[i] * rpa_x - ta_0_xy_1[i] * pc_x[i];

        ta_x_xz_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_0_xz_0[i] * rpa_x - ta_0_xz_1[i] * pc_x[i];

        ta_x_yy_0[i] = ta_0_yy_0[i] * rpa_x - ta_0_yy_1[i] * pc_x[i];

        ta_x_yz_0[i] = ta_0_yz_0[i] * rpa_x - ta_0_yz_1[i] * pc_x[i];

        ta_x_zz_0[i] = ta_0_zz_0[i] * rpa_x - ta_0_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto ta_y_xx_0 = &(buffer.data()[(idx_npot_0_pd + 6) * nelems]);

    auto ta_y_xy_0 = &(buffer.data()[(idx_npot_0_pd + 7) * nelems]);

    auto ta_y_xz_0 = &(buffer.data()[(idx_npot_0_pd + 8) * nelems]);

    auto ta_y_yy_0 = &(buffer.data()[(idx_npot_0_pd + 9) * nelems]);

    auto ta_y_yz_0 = &(buffer.data()[(idx_npot_0_pd + 10) * nelems]);

    auto ta_y_zz_0 = &(buffer.data()[(idx_npot_0_pd + 11) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_y_xx_0[i] = ta_0_xx_0[i] * rpa_y - ta_0_xx_1[i] * pc_y[i];

        ta_y_xy_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_0_xy_0[i] * rpa_y - ta_0_xy_1[i] * pc_y[i];

        ta_y_xz_0[i] = ta_0_xz_0[i] * rpa_y - ta_0_xz_1[i] * pc_y[i];

        ta_y_yy_0[i] = 2.0 * ta_0_y_0[i] * fe_0 - 2.0 * ta_0_y_1[i] * fe_0 + ta_0_yy_0[i] * rpa_y - ta_0_yy_1[i] * pc_y[i];

        ta_y_yz_0[i] = ta_0_z_0[i] * fe_0 - ta_0_z_1[i] * fe_0 + ta_0_yz_0[i] * rpa_y - ta_0_yz_1[i] * pc_y[i];

        ta_y_zz_0[i] = ta_0_zz_0[i] * rpa_y - ta_0_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto ta_z_xx_0 = &(buffer.data()[(idx_npot_0_pd + 12) * nelems]);

    auto ta_z_xy_0 = &(buffer.data()[(idx_npot_0_pd + 13) * nelems]);

    auto ta_z_xz_0 = &(buffer.data()[(idx_npot_0_pd + 14) * nelems]);

    auto ta_z_yy_0 = &(buffer.data()[(idx_npot_0_pd + 15) * nelems]);

    auto ta_z_yz_0 = &(buffer.data()[(idx_npot_0_pd + 16) * nelems]);

    auto ta_z_zz_0 = &(buffer.data()[(idx_npot_0_pd + 17) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_z_xx_0[i] = ta_0_xx_0[i] * rpa_z - ta_0_xx_1[i] * pc_z[i];

        ta_z_xy_0[i] = ta_0_xy_0[i] * rpa_z - ta_0_xy_1[i] * pc_z[i];

        ta_z_xz_0[i] = ta_0_x_0[i] * fe_0 - ta_0_x_1[i] * fe_0 + ta_0_xz_0[i] * rpa_z - ta_0_xz_1[i] * pc_z[i];

        ta_z_yy_0[i] = ta_0_yy_0[i] * rpa_z - ta_0_yy_1[i] * pc_z[i];

        ta_z_yz_0[i] = ta_0_y_0[i] * fe_0 - ta_0_y_1[i] * fe_0 + ta_0_yz_0[i] * rpa_z - ta_0_yz_1[i] * pc_z[i];

        ta_z_zz_0[i] = 2.0 * ta_0_z_0[i] * fe_0 - 2.0 * ta_0_z_1[i] * fe_0 + ta_0_zz_0[i] * rpa_z - ta_0_zz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

