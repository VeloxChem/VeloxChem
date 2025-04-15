#include "NuclearPotentialGridPrimRecDD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_dd(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_dd,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_1_sd,
                                       const size_t idx_npot_0_pp,
                                       const size_t idx_npot_1_pp,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_1_pd,
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

    // Set up components of auxiliary buffer : PD

    auto ta_x_xx_0 = &(buffer.data()[idx_npot_0_pd * nelems]);

    auto ta_x_xy_0 = &(buffer.data()[(idx_npot_0_pd + 1) * nelems]);

    auto ta_x_xz_0 = &(buffer.data()[(idx_npot_0_pd + 2) * nelems]);

    auto ta_x_yy_0 = &(buffer.data()[(idx_npot_0_pd + 3) * nelems]);

    auto ta_x_yz_0 = &(buffer.data()[(idx_npot_0_pd + 4) * nelems]);

    auto ta_x_zz_0 = &(buffer.data()[(idx_npot_0_pd + 5) * nelems]);

    auto ta_y_xx_0 = &(buffer.data()[(idx_npot_0_pd + 6) * nelems]);

    auto ta_y_xy_0 = &(buffer.data()[(idx_npot_0_pd + 7) * nelems]);

    auto ta_y_xz_0 = &(buffer.data()[(idx_npot_0_pd + 8) * nelems]);

    auto ta_y_yy_0 = &(buffer.data()[(idx_npot_0_pd + 9) * nelems]);

    auto ta_y_yz_0 = &(buffer.data()[(idx_npot_0_pd + 10) * nelems]);

    auto ta_y_zz_0 = &(buffer.data()[(idx_npot_0_pd + 11) * nelems]);

    auto ta_z_xx_0 = &(buffer.data()[(idx_npot_0_pd + 12) * nelems]);

    auto ta_z_xy_0 = &(buffer.data()[(idx_npot_0_pd + 13) * nelems]);

    auto ta_z_xz_0 = &(buffer.data()[(idx_npot_0_pd + 14) * nelems]);

    auto ta_z_yy_0 = &(buffer.data()[(idx_npot_0_pd + 15) * nelems]);

    auto ta_z_yz_0 = &(buffer.data()[(idx_npot_0_pd + 16) * nelems]);

    auto ta_z_zz_0 = &(buffer.data()[(idx_npot_0_pd + 17) * nelems]);

    // Set up components of auxiliary buffer : PD

    auto ta_x_xx_1 = &(buffer.data()[idx_npot_1_pd * nelems]);

    auto ta_x_xy_1 = &(buffer.data()[(idx_npot_1_pd + 1) * nelems]);

    auto ta_x_xz_1 = &(buffer.data()[(idx_npot_1_pd + 2) * nelems]);

    auto ta_x_yy_1 = &(buffer.data()[(idx_npot_1_pd + 3) * nelems]);

    auto ta_x_yz_1 = &(buffer.data()[(idx_npot_1_pd + 4) * nelems]);

    auto ta_x_zz_1 = &(buffer.data()[(idx_npot_1_pd + 5) * nelems]);

    auto ta_y_xx_1 = &(buffer.data()[(idx_npot_1_pd + 6) * nelems]);

    auto ta_y_xy_1 = &(buffer.data()[(idx_npot_1_pd + 7) * nelems]);

    auto ta_y_xz_1 = &(buffer.data()[(idx_npot_1_pd + 8) * nelems]);

    auto ta_y_yy_1 = &(buffer.data()[(idx_npot_1_pd + 9) * nelems]);

    auto ta_y_yz_1 = &(buffer.data()[(idx_npot_1_pd + 10) * nelems]);

    auto ta_y_zz_1 = &(buffer.data()[(idx_npot_1_pd + 11) * nelems]);

    auto ta_z_xx_1 = &(buffer.data()[(idx_npot_1_pd + 12) * nelems]);

    auto ta_z_xy_1 = &(buffer.data()[(idx_npot_1_pd + 13) * nelems]);

    auto ta_z_xz_1 = &(buffer.data()[(idx_npot_1_pd + 14) * nelems]);

    auto ta_z_yy_1 = &(buffer.data()[(idx_npot_1_pd + 15) * nelems]);

    auto ta_z_yz_1 = &(buffer.data()[(idx_npot_1_pd + 16) * nelems]);

    auto ta_z_zz_1 = &(buffer.data()[(idx_npot_1_pd + 17) * nelems]);

    // Set up 0-6 components of targeted buffer : DD

    auto ta_xx_xx_0 = &(buffer.data()[idx_npot_0_dd * nelems]);

    auto ta_xx_xy_0 = &(buffer.data()[(idx_npot_0_dd + 1) * nelems]);

    auto ta_xx_xz_0 = &(buffer.data()[(idx_npot_0_dd + 2) * nelems]);

    auto ta_xx_yy_0 = &(buffer.data()[(idx_npot_0_dd + 3) * nelems]);

    auto ta_xx_yz_0 = &(buffer.data()[(idx_npot_0_dd + 4) * nelems]);

    auto ta_xx_zz_0 = &(buffer.data()[(idx_npot_0_dd + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xx_xx_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + 2.0 * ta_x_x_0[i] * fe_0 - 2.0 * ta_x_x_1[i] * fe_0 + ta_x_xx_0[i] * rpa_x - ta_x_xx_1[i] * pc_x[i];

        ta_xx_xy_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_x_y_0[i] * fe_0 - ta_x_y_1[i] * fe_0 + ta_x_xy_0[i] * rpa_x - ta_x_xy_1[i] * pc_x[i];

        ta_xx_xz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_x_z_0[i] * fe_0 - ta_x_z_1[i] * fe_0 + ta_x_xz_0[i] * rpa_x - ta_x_xz_1[i] * pc_x[i];

        ta_xx_yy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_x_yy_0[i] * rpa_x - ta_x_yy_1[i] * pc_x[i];

        ta_xx_yz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_x_yz_0[i] * rpa_x - ta_x_yz_1[i] * pc_x[i];

        ta_xx_zz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_x_zz_0[i] * rpa_x - ta_x_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto ta_xy_xx_0 = &(buffer.data()[(idx_npot_0_dd + 6) * nelems]);

    auto ta_xy_xy_0 = &(buffer.data()[(idx_npot_0_dd + 7) * nelems]);

    auto ta_xy_xz_0 = &(buffer.data()[(idx_npot_0_dd + 8) * nelems]);

    auto ta_xy_yy_0 = &(buffer.data()[(idx_npot_0_dd + 9) * nelems]);

    auto ta_xy_yz_0 = &(buffer.data()[(idx_npot_0_dd + 10) * nelems]);

    auto ta_xy_zz_0 = &(buffer.data()[(idx_npot_0_dd + 11) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xy_xx_0[i] = ta_x_xx_0[i] * rpa_y - ta_x_xx_1[i] * pc_y[i];

        ta_xy_xy_0[i] = ta_y_y_0[i] * fe_0 - ta_y_y_1[i] * fe_0 + ta_y_xy_0[i] * rpa_x - ta_y_xy_1[i] * pc_x[i];

        ta_xy_xz_0[i] = ta_x_xz_0[i] * rpa_y - ta_x_xz_1[i] * pc_y[i];

        ta_xy_yy_0[i] = ta_y_yy_0[i] * rpa_x - ta_y_yy_1[i] * pc_x[i];

        ta_xy_yz_0[i] = ta_y_yz_0[i] * rpa_x - ta_y_yz_1[i] * pc_x[i];

        ta_xy_zz_0[i] = ta_y_zz_0[i] * rpa_x - ta_y_zz_1[i] * pc_x[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto ta_xz_xx_0 = &(buffer.data()[(idx_npot_0_dd + 12) * nelems]);

    auto ta_xz_xy_0 = &(buffer.data()[(idx_npot_0_dd + 13) * nelems]);

    auto ta_xz_xz_0 = &(buffer.data()[(idx_npot_0_dd + 14) * nelems]);

    auto ta_xz_yy_0 = &(buffer.data()[(idx_npot_0_dd + 15) * nelems]);

    auto ta_xz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 16) * nelems]);

    auto ta_xz_zz_0 = &(buffer.data()[(idx_npot_0_dd + 17) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xz_xx_0[i] = ta_x_xx_0[i] * rpa_z - ta_x_xx_1[i] * pc_z[i];

        ta_xz_xy_0[i] = ta_x_xy_0[i] * rpa_z - ta_x_xy_1[i] * pc_z[i];

        ta_xz_xz_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_z_xz_0[i] * rpa_x - ta_z_xz_1[i] * pc_x[i];

        ta_xz_yy_0[i] = ta_z_yy_0[i] * rpa_x - ta_z_yy_1[i] * pc_x[i];

        ta_xz_yz_0[i] = ta_z_yz_0[i] * rpa_x - ta_z_yz_1[i] * pc_x[i];

        ta_xz_zz_0[i] = ta_z_zz_0[i] * rpa_x - ta_z_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto ta_yy_xx_0 = &(buffer.data()[(idx_npot_0_dd + 18) * nelems]);

    auto ta_yy_xy_0 = &(buffer.data()[(idx_npot_0_dd + 19) * nelems]);

    auto ta_yy_xz_0 = &(buffer.data()[(idx_npot_0_dd + 20) * nelems]);

    auto ta_yy_yy_0 = &(buffer.data()[(idx_npot_0_dd + 21) * nelems]);

    auto ta_yy_yz_0 = &(buffer.data()[(idx_npot_0_dd + 22) * nelems]);

    auto ta_yy_zz_0 = &(buffer.data()[(idx_npot_0_dd + 23) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yy_xx_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_y_xx_0[i] * rpa_y - ta_y_xx_1[i] * pc_y[i];

        ta_yy_xy_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_y_x_0[i] * fe_0 - ta_y_x_1[i] * fe_0 + ta_y_xy_0[i] * rpa_y - ta_y_xy_1[i] * pc_y[i];

        ta_yy_xz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_y_xz_0[i] * rpa_y - ta_y_xz_1[i] * pc_y[i];

        ta_yy_yy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + 2.0 * ta_y_y_0[i] * fe_0 - 2.0 * ta_y_y_1[i] * fe_0 + ta_y_yy_0[i] * rpa_y - ta_y_yy_1[i] * pc_y[i];

        ta_yy_yz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_y_z_0[i] * fe_0 - ta_y_z_1[i] * fe_0 + ta_y_yz_0[i] * rpa_y - ta_y_yz_1[i] * pc_y[i];

        ta_yy_zz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_y_zz_0[i] * rpa_y - ta_y_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto ta_yz_xx_0 = &(buffer.data()[(idx_npot_0_dd + 24) * nelems]);

    auto ta_yz_xy_0 = &(buffer.data()[(idx_npot_0_dd + 25) * nelems]);

    auto ta_yz_xz_0 = &(buffer.data()[(idx_npot_0_dd + 26) * nelems]);

    auto ta_yz_yy_0 = &(buffer.data()[(idx_npot_0_dd + 27) * nelems]);

    auto ta_yz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 28) * nelems]);

    auto ta_yz_zz_0 = &(buffer.data()[(idx_npot_0_dd + 29) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yz_xx_0[i] = ta_z_xx_0[i] * rpa_y - ta_z_xx_1[i] * pc_y[i];

        ta_yz_xy_0[i] = ta_y_xy_0[i] * rpa_z - ta_y_xy_1[i] * pc_z[i];

        ta_yz_xz_0[i] = ta_z_xz_0[i] * rpa_y - ta_z_xz_1[i] * pc_y[i];

        ta_yz_yy_0[i] = ta_y_yy_0[i] * rpa_z - ta_y_yy_1[i] * pc_z[i];

        ta_yz_yz_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_z_yz_0[i] * rpa_y - ta_z_yz_1[i] * pc_y[i];

        ta_yz_zz_0[i] = ta_z_zz_0[i] * rpa_y - ta_z_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto ta_zz_xx_0 = &(buffer.data()[(idx_npot_0_dd + 30) * nelems]);

    auto ta_zz_xy_0 = &(buffer.data()[(idx_npot_0_dd + 31) * nelems]);

    auto ta_zz_xz_0 = &(buffer.data()[(idx_npot_0_dd + 32) * nelems]);

    auto ta_zz_yy_0 = &(buffer.data()[(idx_npot_0_dd + 33) * nelems]);

    auto ta_zz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 34) * nelems]);

    auto ta_zz_zz_0 = &(buffer.data()[(idx_npot_0_dd + 35) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_zz_xx_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_z_xx_0[i] * rpa_z - ta_z_xx_1[i] * pc_z[i];

        ta_zz_xy_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_z_xy_0[i] * rpa_z - ta_z_xy_1[i] * pc_z[i];

        ta_zz_xz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_z_x_0[i] * fe_0 - ta_z_x_1[i] * fe_0 + ta_z_xz_0[i] * rpa_z - ta_z_xz_1[i] * pc_z[i];

        ta_zz_yy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_z_yy_0[i] * rpa_z - ta_z_yy_1[i] * pc_z[i];

        ta_zz_yz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_z_y_0[i] * fe_0 - ta_z_y_1[i] * fe_0 + ta_z_yz_0[i] * rpa_z - ta_z_yz_1[i] * pc_z[i];

        ta_zz_zz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + 2.0 * ta_z_z_0[i] * fe_0 - 2.0 * ta_z_z_1[i] * fe_0 + ta_z_zz_0[i] * rpa_z - ta_z_zz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

