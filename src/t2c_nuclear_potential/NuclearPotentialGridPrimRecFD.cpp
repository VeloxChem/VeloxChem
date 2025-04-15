#include "NuclearPotentialGridPrimRecFD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_fd(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_fd,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_1_pd,
                                       const size_t idx_npot_0_dp,
                                       const size_t idx_npot_1_dp,
                                       const size_t idx_npot_0_dd,
                                       const size_t idx_npot_1_dd,
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

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_0 = &(buffer.data()[idx_npot_0_dp * nelems]);

    auto ta_xx_y_0 = &(buffer.data()[(idx_npot_0_dp + 1) * nelems]);

    auto ta_xx_z_0 = &(buffer.data()[(idx_npot_0_dp + 2) * nelems]);

    auto ta_yy_x_0 = &(buffer.data()[(idx_npot_0_dp + 9) * nelems]);

    auto ta_yy_y_0 = &(buffer.data()[(idx_npot_0_dp + 10) * nelems]);

    auto ta_yy_z_0 = &(buffer.data()[(idx_npot_0_dp + 11) * nelems]);

    auto ta_zz_x_0 = &(buffer.data()[(idx_npot_0_dp + 15) * nelems]);

    auto ta_zz_y_0 = &(buffer.data()[(idx_npot_0_dp + 16) * nelems]);

    auto ta_zz_z_0 = &(buffer.data()[(idx_npot_0_dp + 17) * nelems]);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_1 = &(buffer.data()[idx_npot_1_dp * nelems]);

    auto ta_xx_y_1 = &(buffer.data()[(idx_npot_1_dp + 1) * nelems]);

    auto ta_xx_z_1 = &(buffer.data()[(idx_npot_1_dp + 2) * nelems]);

    auto ta_yy_x_1 = &(buffer.data()[(idx_npot_1_dp + 9) * nelems]);

    auto ta_yy_y_1 = &(buffer.data()[(idx_npot_1_dp + 10) * nelems]);

    auto ta_yy_z_1 = &(buffer.data()[(idx_npot_1_dp + 11) * nelems]);

    auto ta_zz_x_1 = &(buffer.data()[(idx_npot_1_dp + 15) * nelems]);

    auto ta_zz_y_1 = &(buffer.data()[(idx_npot_1_dp + 16) * nelems]);

    auto ta_zz_z_1 = &(buffer.data()[(idx_npot_1_dp + 17) * nelems]);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_0 = &(buffer.data()[idx_npot_0_dd * nelems]);

    auto ta_xx_xy_0 = &(buffer.data()[(idx_npot_0_dd + 1) * nelems]);

    auto ta_xx_xz_0 = &(buffer.data()[(idx_npot_0_dd + 2) * nelems]);

    auto ta_xx_yy_0 = &(buffer.data()[(idx_npot_0_dd + 3) * nelems]);

    auto ta_xx_yz_0 = &(buffer.data()[(idx_npot_0_dd + 4) * nelems]);

    auto ta_xx_zz_0 = &(buffer.data()[(idx_npot_0_dd + 5) * nelems]);

    auto ta_xy_xy_0 = &(buffer.data()[(idx_npot_0_dd + 7) * nelems]);

    auto ta_xy_yy_0 = &(buffer.data()[(idx_npot_0_dd + 9) * nelems]);

    auto ta_xy_yz_0 = &(buffer.data()[(idx_npot_0_dd + 10) * nelems]);

    auto ta_xz_xx_0 = &(buffer.data()[(idx_npot_0_dd + 12) * nelems]);

    auto ta_xz_xz_0 = &(buffer.data()[(idx_npot_0_dd + 14) * nelems]);

    auto ta_xz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 16) * nelems]);

    auto ta_xz_zz_0 = &(buffer.data()[(idx_npot_0_dd + 17) * nelems]);

    auto ta_yy_xx_0 = &(buffer.data()[(idx_npot_0_dd + 18) * nelems]);

    auto ta_yy_xy_0 = &(buffer.data()[(idx_npot_0_dd + 19) * nelems]);

    auto ta_yy_xz_0 = &(buffer.data()[(idx_npot_0_dd + 20) * nelems]);

    auto ta_yy_yy_0 = &(buffer.data()[(idx_npot_0_dd + 21) * nelems]);

    auto ta_yy_yz_0 = &(buffer.data()[(idx_npot_0_dd + 22) * nelems]);

    auto ta_yy_zz_0 = &(buffer.data()[(idx_npot_0_dd + 23) * nelems]);

    auto ta_yz_xz_0 = &(buffer.data()[(idx_npot_0_dd + 26) * nelems]);

    auto ta_yz_yy_0 = &(buffer.data()[(idx_npot_0_dd + 27) * nelems]);

    auto ta_yz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 28) * nelems]);

    auto ta_yz_zz_0 = &(buffer.data()[(idx_npot_0_dd + 29) * nelems]);

    auto ta_zz_xx_0 = &(buffer.data()[(idx_npot_0_dd + 30) * nelems]);

    auto ta_zz_xy_0 = &(buffer.data()[(idx_npot_0_dd + 31) * nelems]);

    auto ta_zz_xz_0 = &(buffer.data()[(idx_npot_0_dd + 32) * nelems]);

    auto ta_zz_yy_0 = &(buffer.data()[(idx_npot_0_dd + 33) * nelems]);

    auto ta_zz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 34) * nelems]);

    auto ta_zz_zz_0 = &(buffer.data()[(idx_npot_0_dd + 35) * nelems]);

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_1 = &(buffer.data()[idx_npot_1_dd * nelems]);

    auto ta_xx_xy_1 = &(buffer.data()[(idx_npot_1_dd + 1) * nelems]);

    auto ta_xx_xz_1 = &(buffer.data()[(idx_npot_1_dd + 2) * nelems]);

    auto ta_xx_yy_1 = &(buffer.data()[(idx_npot_1_dd + 3) * nelems]);

    auto ta_xx_yz_1 = &(buffer.data()[(idx_npot_1_dd + 4) * nelems]);

    auto ta_xx_zz_1 = &(buffer.data()[(idx_npot_1_dd + 5) * nelems]);

    auto ta_xy_xy_1 = &(buffer.data()[(idx_npot_1_dd + 7) * nelems]);

    auto ta_xy_yy_1 = &(buffer.data()[(idx_npot_1_dd + 9) * nelems]);

    auto ta_xy_yz_1 = &(buffer.data()[(idx_npot_1_dd + 10) * nelems]);

    auto ta_xz_xx_1 = &(buffer.data()[(idx_npot_1_dd + 12) * nelems]);

    auto ta_xz_xz_1 = &(buffer.data()[(idx_npot_1_dd + 14) * nelems]);

    auto ta_xz_yz_1 = &(buffer.data()[(idx_npot_1_dd + 16) * nelems]);

    auto ta_xz_zz_1 = &(buffer.data()[(idx_npot_1_dd + 17) * nelems]);

    auto ta_yy_xx_1 = &(buffer.data()[(idx_npot_1_dd + 18) * nelems]);

    auto ta_yy_xy_1 = &(buffer.data()[(idx_npot_1_dd + 19) * nelems]);

    auto ta_yy_xz_1 = &(buffer.data()[(idx_npot_1_dd + 20) * nelems]);

    auto ta_yy_yy_1 = &(buffer.data()[(idx_npot_1_dd + 21) * nelems]);

    auto ta_yy_yz_1 = &(buffer.data()[(idx_npot_1_dd + 22) * nelems]);

    auto ta_yy_zz_1 = &(buffer.data()[(idx_npot_1_dd + 23) * nelems]);

    auto ta_yz_xz_1 = &(buffer.data()[(idx_npot_1_dd + 26) * nelems]);

    auto ta_yz_yy_1 = &(buffer.data()[(idx_npot_1_dd + 27) * nelems]);

    auto ta_yz_yz_1 = &(buffer.data()[(idx_npot_1_dd + 28) * nelems]);

    auto ta_yz_zz_1 = &(buffer.data()[(idx_npot_1_dd + 29) * nelems]);

    auto ta_zz_xx_1 = &(buffer.data()[(idx_npot_1_dd + 30) * nelems]);

    auto ta_zz_xy_1 = &(buffer.data()[(idx_npot_1_dd + 31) * nelems]);

    auto ta_zz_xz_1 = &(buffer.data()[(idx_npot_1_dd + 32) * nelems]);

    auto ta_zz_yy_1 = &(buffer.data()[(idx_npot_1_dd + 33) * nelems]);

    auto ta_zz_yz_1 = &(buffer.data()[(idx_npot_1_dd + 34) * nelems]);

    auto ta_zz_zz_1 = &(buffer.data()[(idx_npot_1_dd + 35) * nelems]);

    // Set up 0-6 components of targeted buffer : FD

    auto ta_xxx_xx_0 = &(buffer.data()[idx_npot_0_fd * nelems]);

    auto ta_xxx_xy_0 = &(buffer.data()[(idx_npot_0_fd + 1) * nelems]);

    auto ta_xxx_xz_0 = &(buffer.data()[(idx_npot_0_fd + 2) * nelems]);

    auto ta_xxx_yy_0 = &(buffer.data()[(idx_npot_0_fd + 3) * nelems]);

    auto ta_xxx_yz_0 = &(buffer.data()[(idx_npot_0_fd + 4) * nelems]);

    auto ta_xxx_zz_0 = &(buffer.data()[(idx_npot_0_fd + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxx_xx_0[i] = 2.0 * ta_x_xx_0[i] * fe_0 - 2.0 * ta_x_xx_1[i] * fe_0 + 2.0 * ta_xx_x_0[i] * fe_0 - 2.0 * ta_xx_x_1[i] * fe_0 + ta_xx_xx_0[i] * rpa_x - ta_xx_xx_1[i] * pc_x[i];

        ta_xxx_xy_0[i] = 2.0 * ta_x_xy_0[i] * fe_0 - 2.0 * ta_x_xy_1[i] * fe_0 + ta_xx_y_0[i] * fe_0 - ta_xx_y_1[i] * fe_0 + ta_xx_xy_0[i] * rpa_x - ta_xx_xy_1[i] * pc_x[i];

        ta_xxx_xz_0[i] = 2.0 * ta_x_xz_0[i] * fe_0 - 2.0 * ta_x_xz_1[i] * fe_0 + ta_xx_z_0[i] * fe_0 - ta_xx_z_1[i] * fe_0 + ta_xx_xz_0[i] * rpa_x - ta_xx_xz_1[i] * pc_x[i];

        ta_xxx_yy_0[i] = 2.0 * ta_x_yy_0[i] * fe_0 - 2.0 * ta_x_yy_1[i] * fe_0 + ta_xx_yy_0[i] * rpa_x - ta_xx_yy_1[i] * pc_x[i];

        ta_xxx_yz_0[i] = 2.0 * ta_x_yz_0[i] * fe_0 - 2.0 * ta_x_yz_1[i] * fe_0 + ta_xx_yz_0[i] * rpa_x - ta_xx_yz_1[i] * pc_x[i];

        ta_xxx_zz_0[i] = 2.0 * ta_x_zz_0[i] * fe_0 - 2.0 * ta_x_zz_1[i] * fe_0 + ta_xx_zz_0[i] * rpa_x - ta_xx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto ta_xxy_xx_0 = &(buffer.data()[(idx_npot_0_fd + 6) * nelems]);

    auto ta_xxy_xy_0 = &(buffer.data()[(idx_npot_0_fd + 7) * nelems]);

    auto ta_xxy_xz_0 = &(buffer.data()[(idx_npot_0_fd + 8) * nelems]);

    auto ta_xxy_yy_0 = &(buffer.data()[(idx_npot_0_fd + 9) * nelems]);

    auto ta_xxy_yz_0 = &(buffer.data()[(idx_npot_0_fd + 10) * nelems]);

    auto ta_xxy_zz_0 = &(buffer.data()[(idx_npot_0_fd + 11) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxy_xx_0[i] = ta_xx_xx_0[i] * rpa_y - ta_xx_xx_1[i] * pc_y[i];

        ta_xxy_xy_0[i] = ta_xx_x_0[i] * fe_0 - ta_xx_x_1[i] * fe_0 + ta_xx_xy_0[i] * rpa_y - ta_xx_xy_1[i] * pc_y[i];

        ta_xxy_xz_0[i] = ta_xx_xz_0[i] * rpa_y - ta_xx_xz_1[i] * pc_y[i];

        ta_xxy_yy_0[i] = ta_y_yy_0[i] * fe_0 - ta_y_yy_1[i] * fe_0 + ta_xy_yy_0[i] * rpa_x - ta_xy_yy_1[i] * pc_x[i];

        ta_xxy_yz_0[i] = ta_y_yz_0[i] * fe_0 - ta_y_yz_1[i] * fe_0 + ta_xy_yz_0[i] * rpa_x - ta_xy_yz_1[i] * pc_x[i];

        ta_xxy_zz_0[i] = ta_xx_zz_0[i] * rpa_y - ta_xx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto ta_xxz_xx_0 = &(buffer.data()[(idx_npot_0_fd + 12) * nelems]);

    auto ta_xxz_xy_0 = &(buffer.data()[(idx_npot_0_fd + 13) * nelems]);

    auto ta_xxz_xz_0 = &(buffer.data()[(idx_npot_0_fd + 14) * nelems]);

    auto ta_xxz_yy_0 = &(buffer.data()[(idx_npot_0_fd + 15) * nelems]);

    auto ta_xxz_yz_0 = &(buffer.data()[(idx_npot_0_fd + 16) * nelems]);

    auto ta_xxz_zz_0 = &(buffer.data()[(idx_npot_0_fd + 17) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxz_xx_0[i] = ta_xx_xx_0[i] * rpa_z - ta_xx_xx_1[i] * pc_z[i];

        ta_xxz_xy_0[i] = ta_xx_xy_0[i] * rpa_z - ta_xx_xy_1[i] * pc_z[i];

        ta_xxz_xz_0[i] = ta_xx_x_0[i] * fe_0 - ta_xx_x_1[i] * fe_0 + ta_xx_xz_0[i] * rpa_z - ta_xx_xz_1[i] * pc_z[i];

        ta_xxz_yy_0[i] = ta_xx_yy_0[i] * rpa_z - ta_xx_yy_1[i] * pc_z[i];

        ta_xxz_yz_0[i] = ta_z_yz_0[i] * fe_0 - ta_z_yz_1[i] * fe_0 + ta_xz_yz_0[i] * rpa_x - ta_xz_yz_1[i] * pc_x[i];

        ta_xxz_zz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_xz_zz_0[i] * rpa_x - ta_xz_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto ta_xyy_xx_0 = &(buffer.data()[(idx_npot_0_fd + 18) * nelems]);

    auto ta_xyy_xy_0 = &(buffer.data()[(idx_npot_0_fd + 19) * nelems]);

    auto ta_xyy_xz_0 = &(buffer.data()[(idx_npot_0_fd + 20) * nelems]);

    auto ta_xyy_yy_0 = &(buffer.data()[(idx_npot_0_fd + 21) * nelems]);

    auto ta_xyy_yz_0 = &(buffer.data()[(idx_npot_0_fd + 22) * nelems]);

    auto ta_xyy_zz_0 = &(buffer.data()[(idx_npot_0_fd + 23) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xyy_xx_0[i] = 2.0 * ta_yy_x_0[i] * fe_0 - 2.0 * ta_yy_x_1[i] * fe_0 + ta_yy_xx_0[i] * rpa_x - ta_yy_xx_1[i] * pc_x[i];

        ta_xyy_xy_0[i] = ta_yy_y_0[i] * fe_0 - ta_yy_y_1[i] * fe_0 + ta_yy_xy_0[i] * rpa_x - ta_yy_xy_1[i] * pc_x[i];

        ta_xyy_xz_0[i] = ta_yy_z_0[i] * fe_0 - ta_yy_z_1[i] * fe_0 + ta_yy_xz_0[i] * rpa_x - ta_yy_xz_1[i] * pc_x[i];

        ta_xyy_yy_0[i] = ta_yy_yy_0[i] * rpa_x - ta_yy_yy_1[i] * pc_x[i];

        ta_xyy_yz_0[i] = ta_yy_yz_0[i] * rpa_x - ta_yy_yz_1[i] * pc_x[i];

        ta_xyy_zz_0[i] = ta_yy_zz_0[i] * rpa_x - ta_yy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto ta_xyz_xx_0 = &(buffer.data()[(idx_npot_0_fd + 24) * nelems]);

    auto ta_xyz_xy_0 = &(buffer.data()[(idx_npot_0_fd + 25) * nelems]);

    auto ta_xyz_xz_0 = &(buffer.data()[(idx_npot_0_fd + 26) * nelems]);

    auto ta_xyz_yy_0 = &(buffer.data()[(idx_npot_0_fd + 27) * nelems]);

    auto ta_xyz_yz_0 = &(buffer.data()[(idx_npot_0_fd + 28) * nelems]);

    auto ta_xyz_zz_0 = &(buffer.data()[(idx_npot_0_fd + 29) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyz_xx_0[i] = ta_xz_xx_0[i] * rpa_y - ta_xz_xx_1[i] * pc_y[i];

        ta_xyz_xy_0[i] = ta_xy_xy_0[i] * rpa_z - ta_xy_xy_1[i] * pc_z[i];

        ta_xyz_xz_0[i] = ta_xz_xz_0[i] * rpa_y - ta_xz_xz_1[i] * pc_y[i];

        ta_xyz_yy_0[i] = ta_yz_yy_0[i] * rpa_x - ta_yz_yy_1[i] * pc_x[i];

        ta_xyz_yz_0[i] = ta_yz_yz_0[i] * rpa_x - ta_yz_yz_1[i] * pc_x[i];

        ta_xyz_zz_0[i] = ta_yz_zz_0[i] * rpa_x - ta_yz_zz_1[i] * pc_x[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto ta_xzz_xx_0 = &(buffer.data()[(idx_npot_0_fd + 30) * nelems]);

    auto ta_xzz_xy_0 = &(buffer.data()[(idx_npot_0_fd + 31) * nelems]);

    auto ta_xzz_xz_0 = &(buffer.data()[(idx_npot_0_fd + 32) * nelems]);

    auto ta_xzz_yy_0 = &(buffer.data()[(idx_npot_0_fd + 33) * nelems]);

    auto ta_xzz_yz_0 = &(buffer.data()[(idx_npot_0_fd + 34) * nelems]);

    auto ta_xzz_zz_0 = &(buffer.data()[(idx_npot_0_fd + 35) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xzz_xx_0[i] = 2.0 * ta_zz_x_0[i] * fe_0 - 2.0 * ta_zz_x_1[i] * fe_0 + ta_zz_xx_0[i] * rpa_x - ta_zz_xx_1[i] * pc_x[i];

        ta_xzz_xy_0[i] = ta_zz_y_0[i] * fe_0 - ta_zz_y_1[i] * fe_0 + ta_zz_xy_0[i] * rpa_x - ta_zz_xy_1[i] * pc_x[i];

        ta_xzz_xz_0[i] = ta_zz_z_0[i] * fe_0 - ta_zz_z_1[i] * fe_0 + ta_zz_xz_0[i] * rpa_x - ta_zz_xz_1[i] * pc_x[i];

        ta_xzz_yy_0[i] = ta_zz_yy_0[i] * rpa_x - ta_zz_yy_1[i] * pc_x[i];

        ta_xzz_yz_0[i] = ta_zz_yz_0[i] * rpa_x - ta_zz_yz_1[i] * pc_x[i];

        ta_xzz_zz_0[i] = ta_zz_zz_0[i] * rpa_x - ta_zz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto ta_yyy_xx_0 = &(buffer.data()[(idx_npot_0_fd + 36) * nelems]);

    auto ta_yyy_xy_0 = &(buffer.data()[(idx_npot_0_fd + 37) * nelems]);

    auto ta_yyy_xz_0 = &(buffer.data()[(idx_npot_0_fd + 38) * nelems]);

    auto ta_yyy_yy_0 = &(buffer.data()[(idx_npot_0_fd + 39) * nelems]);

    auto ta_yyy_yz_0 = &(buffer.data()[(idx_npot_0_fd + 40) * nelems]);

    auto ta_yyy_zz_0 = &(buffer.data()[(idx_npot_0_fd + 41) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yyy_xx_0[i] = 2.0 * ta_y_xx_0[i] * fe_0 - 2.0 * ta_y_xx_1[i] * fe_0 + ta_yy_xx_0[i] * rpa_y - ta_yy_xx_1[i] * pc_y[i];

        ta_yyy_xy_0[i] = 2.0 * ta_y_xy_0[i] * fe_0 - 2.0 * ta_y_xy_1[i] * fe_0 + ta_yy_x_0[i] * fe_0 - ta_yy_x_1[i] * fe_0 + ta_yy_xy_0[i] * rpa_y - ta_yy_xy_1[i] * pc_y[i];

        ta_yyy_xz_0[i] = 2.0 * ta_y_xz_0[i] * fe_0 - 2.0 * ta_y_xz_1[i] * fe_0 + ta_yy_xz_0[i] * rpa_y - ta_yy_xz_1[i] * pc_y[i];

        ta_yyy_yy_0[i] = 2.0 * ta_y_yy_0[i] * fe_0 - 2.0 * ta_y_yy_1[i] * fe_0 + 2.0 * ta_yy_y_0[i] * fe_0 - 2.0 * ta_yy_y_1[i] * fe_0 + ta_yy_yy_0[i] * rpa_y - ta_yy_yy_1[i] * pc_y[i];

        ta_yyy_yz_0[i] = 2.0 * ta_y_yz_0[i] * fe_0 - 2.0 * ta_y_yz_1[i] * fe_0 + ta_yy_z_0[i] * fe_0 - ta_yy_z_1[i] * fe_0 + ta_yy_yz_0[i] * rpa_y - ta_yy_yz_1[i] * pc_y[i];

        ta_yyy_zz_0[i] = 2.0 * ta_y_zz_0[i] * fe_0 - 2.0 * ta_y_zz_1[i] * fe_0 + ta_yy_zz_0[i] * rpa_y - ta_yy_zz_1[i] * pc_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto ta_yyz_xx_0 = &(buffer.data()[(idx_npot_0_fd + 42) * nelems]);

    auto ta_yyz_xy_0 = &(buffer.data()[(idx_npot_0_fd + 43) * nelems]);

    auto ta_yyz_xz_0 = &(buffer.data()[(idx_npot_0_fd + 44) * nelems]);

    auto ta_yyz_yy_0 = &(buffer.data()[(idx_npot_0_fd + 45) * nelems]);

    auto ta_yyz_yz_0 = &(buffer.data()[(idx_npot_0_fd + 46) * nelems]);

    auto ta_yyz_zz_0 = &(buffer.data()[(idx_npot_0_fd + 47) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yyz_xx_0[i] = ta_yy_xx_0[i] * rpa_z - ta_yy_xx_1[i] * pc_z[i];

        ta_yyz_xy_0[i] = ta_yy_xy_0[i] * rpa_z - ta_yy_xy_1[i] * pc_z[i];

        ta_yyz_xz_0[i] = ta_z_xz_0[i] * fe_0 - ta_z_xz_1[i] * fe_0 + ta_yz_xz_0[i] * rpa_y - ta_yz_xz_1[i] * pc_y[i];

        ta_yyz_yy_0[i] = ta_yy_yy_0[i] * rpa_z - ta_yy_yy_1[i] * pc_z[i];

        ta_yyz_yz_0[i] = ta_yy_y_0[i] * fe_0 - ta_yy_y_1[i] * fe_0 + ta_yy_yz_0[i] * rpa_z - ta_yy_yz_1[i] * pc_z[i];

        ta_yyz_zz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_yz_zz_0[i] * rpa_y - ta_yz_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto ta_yzz_xx_0 = &(buffer.data()[(idx_npot_0_fd + 48) * nelems]);

    auto ta_yzz_xy_0 = &(buffer.data()[(idx_npot_0_fd + 49) * nelems]);

    auto ta_yzz_xz_0 = &(buffer.data()[(idx_npot_0_fd + 50) * nelems]);

    auto ta_yzz_yy_0 = &(buffer.data()[(idx_npot_0_fd + 51) * nelems]);

    auto ta_yzz_yz_0 = &(buffer.data()[(idx_npot_0_fd + 52) * nelems]);

    auto ta_yzz_zz_0 = &(buffer.data()[(idx_npot_0_fd + 53) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yzz_xx_0[i] = ta_zz_xx_0[i] * rpa_y - ta_zz_xx_1[i] * pc_y[i];

        ta_yzz_xy_0[i] = ta_zz_x_0[i] * fe_0 - ta_zz_x_1[i] * fe_0 + ta_zz_xy_0[i] * rpa_y - ta_zz_xy_1[i] * pc_y[i];

        ta_yzz_xz_0[i] = ta_zz_xz_0[i] * rpa_y - ta_zz_xz_1[i] * pc_y[i];

        ta_yzz_yy_0[i] = 2.0 * ta_zz_y_0[i] * fe_0 - 2.0 * ta_zz_y_1[i] * fe_0 + ta_zz_yy_0[i] * rpa_y - ta_zz_yy_1[i] * pc_y[i];

        ta_yzz_yz_0[i] = ta_zz_z_0[i] * fe_0 - ta_zz_z_1[i] * fe_0 + ta_zz_yz_0[i] * rpa_y - ta_zz_yz_1[i] * pc_y[i];

        ta_yzz_zz_0[i] = ta_zz_zz_0[i] * rpa_y - ta_zz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto ta_zzz_xx_0 = &(buffer.data()[(idx_npot_0_fd + 54) * nelems]);

    auto ta_zzz_xy_0 = &(buffer.data()[(idx_npot_0_fd + 55) * nelems]);

    auto ta_zzz_xz_0 = &(buffer.data()[(idx_npot_0_fd + 56) * nelems]);

    auto ta_zzz_yy_0 = &(buffer.data()[(idx_npot_0_fd + 57) * nelems]);

    auto ta_zzz_yz_0 = &(buffer.data()[(idx_npot_0_fd + 58) * nelems]);

    auto ta_zzz_zz_0 = &(buffer.data()[(idx_npot_0_fd + 59) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_zzz_xx_0[i] = 2.0 * ta_z_xx_0[i] * fe_0 - 2.0 * ta_z_xx_1[i] * fe_0 + ta_zz_xx_0[i] * rpa_z - ta_zz_xx_1[i] * pc_z[i];

        ta_zzz_xy_0[i] = 2.0 * ta_z_xy_0[i] * fe_0 - 2.0 * ta_z_xy_1[i] * fe_0 + ta_zz_xy_0[i] * rpa_z - ta_zz_xy_1[i] * pc_z[i];

        ta_zzz_xz_0[i] = 2.0 * ta_z_xz_0[i] * fe_0 - 2.0 * ta_z_xz_1[i] * fe_0 + ta_zz_x_0[i] * fe_0 - ta_zz_x_1[i] * fe_0 + ta_zz_xz_0[i] * rpa_z - ta_zz_xz_1[i] * pc_z[i];

        ta_zzz_yy_0[i] = 2.0 * ta_z_yy_0[i] * fe_0 - 2.0 * ta_z_yy_1[i] * fe_0 + ta_zz_yy_0[i] * rpa_z - ta_zz_yy_1[i] * pc_z[i];

        ta_zzz_yz_0[i] = 2.0 * ta_z_yz_0[i] * fe_0 - 2.0 * ta_z_yz_1[i] * fe_0 + ta_zz_y_0[i] * fe_0 - ta_zz_y_1[i] * fe_0 + ta_zz_yz_0[i] * rpa_z - ta_zz_yz_1[i] * pc_z[i];

        ta_zzz_zz_0[i] = 2.0 * ta_z_zz_0[i] * fe_0 - 2.0 * ta_z_zz_1[i] * fe_0 + 2.0 * ta_zz_z_0[i] * fe_0 - 2.0 * ta_zz_z_1[i] * fe_0 + ta_zz_zz_0[i] * rpa_z - ta_zz_zz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

