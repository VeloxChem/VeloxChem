#include "NuclearPotentialGridPrimRecDF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_df(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_df,
                                       const size_t idx_npot_0_sf,
                                       const size_t idx_npot_1_sf,
                                       const size_t idx_npot_0_pd,
                                       const size_t idx_npot_1_pd,
                                       const size_t idx_npot_0_pf,
                                       const size_t idx_npot_1_pf,
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

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_0 = &(buffer.data()[idx_npot_0_sf * nelems]);

    auto ta_0_xxy_0 = &(buffer.data()[(idx_npot_0_sf + 1) * nelems]);

    auto ta_0_xxz_0 = &(buffer.data()[(idx_npot_0_sf + 2) * nelems]);

    auto ta_0_xyy_0 = &(buffer.data()[(idx_npot_0_sf + 3) * nelems]);

    auto ta_0_xyz_0 = &(buffer.data()[(idx_npot_0_sf + 4) * nelems]);

    auto ta_0_xzz_0 = &(buffer.data()[(idx_npot_0_sf + 5) * nelems]);

    auto ta_0_yyy_0 = &(buffer.data()[(idx_npot_0_sf + 6) * nelems]);

    auto ta_0_yyz_0 = &(buffer.data()[(idx_npot_0_sf + 7) * nelems]);

    auto ta_0_yzz_0 = &(buffer.data()[(idx_npot_0_sf + 8) * nelems]);

    auto ta_0_zzz_0 = &(buffer.data()[(idx_npot_0_sf + 9) * nelems]);

    // Set up components of auxiliary buffer : SF

    auto ta_0_xxx_1 = &(buffer.data()[idx_npot_1_sf * nelems]);

    auto ta_0_xxy_1 = &(buffer.data()[(idx_npot_1_sf + 1) * nelems]);

    auto ta_0_xxz_1 = &(buffer.data()[(idx_npot_1_sf + 2) * nelems]);

    auto ta_0_xyy_1 = &(buffer.data()[(idx_npot_1_sf + 3) * nelems]);

    auto ta_0_xyz_1 = &(buffer.data()[(idx_npot_1_sf + 4) * nelems]);

    auto ta_0_xzz_1 = &(buffer.data()[(idx_npot_1_sf + 5) * nelems]);

    auto ta_0_yyy_1 = &(buffer.data()[(idx_npot_1_sf + 6) * nelems]);

    auto ta_0_yyz_1 = &(buffer.data()[(idx_npot_1_sf + 7) * nelems]);

    auto ta_0_yzz_1 = &(buffer.data()[(idx_npot_1_sf + 8) * nelems]);

    auto ta_0_zzz_1 = &(buffer.data()[(idx_npot_1_sf + 9) * nelems]);

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

    // Set up components of auxiliary buffer : PF

    auto ta_x_xxx_0 = &(buffer.data()[idx_npot_0_pf * nelems]);

    auto ta_x_xxy_0 = &(buffer.data()[(idx_npot_0_pf + 1) * nelems]);

    auto ta_x_xxz_0 = &(buffer.data()[(idx_npot_0_pf + 2) * nelems]);

    auto ta_x_xyy_0 = &(buffer.data()[(idx_npot_0_pf + 3) * nelems]);

    auto ta_x_xyz_0 = &(buffer.data()[(idx_npot_0_pf + 4) * nelems]);

    auto ta_x_xzz_0 = &(buffer.data()[(idx_npot_0_pf + 5) * nelems]);

    auto ta_x_yyy_0 = &(buffer.data()[(idx_npot_0_pf + 6) * nelems]);

    auto ta_x_yyz_0 = &(buffer.data()[(idx_npot_0_pf + 7) * nelems]);

    auto ta_x_yzz_0 = &(buffer.data()[(idx_npot_0_pf + 8) * nelems]);

    auto ta_x_zzz_0 = &(buffer.data()[(idx_npot_0_pf + 9) * nelems]);

    auto ta_y_xxx_0 = &(buffer.data()[(idx_npot_0_pf + 10) * nelems]);

    auto ta_y_xxy_0 = &(buffer.data()[(idx_npot_0_pf + 11) * nelems]);

    auto ta_y_xxz_0 = &(buffer.data()[(idx_npot_0_pf + 12) * nelems]);

    auto ta_y_xyy_0 = &(buffer.data()[(idx_npot_0_pf + 13) * nelems]);

    auto ta_y_xyz_0 = &(buffer.data()[(idx_npot_0_pf + 14) * nelems]);

    auto ta_y_xzz_0 = &(buffer.data()[(idx_npot_0_pf + 15) * nelems]);

    auto ta_y_yyy_0 = &(buffer.data()[(idx_npot_0_pf + 16) * nelems]);

    auto ta_y_yyz_0 = &(buffer.data()[(idx_npot_0_pf + 17) * nelems]);

    auto ta_y_yzz_0 = &(buffer.data()[(idx_npot_0_pf + 18) * nelems]);

    auto ta_y_zzz_0 = &(buffer.data()[(idx_npot_0_pf + 19) * nelems]);

    auto ta_z_xxx_0 = &(buffer.data()[(idx_npot_0_pf + 20) * nelems]);

    auto ta_z_xxy_0 = &(buffer.data()[(idx_npot_0_pf + 21) * nelems]);

    auto ta_z_xxz_0 = &(buffer.data()[(idx_npot_0_pf + 22) * nelems]);

    auto ta_z_xyy_0 = &(buffer.data()[(idx_npot_0_pf + 23) * nelems]);

    auto ta_z_xyz_0 = &(buffer.data()[(idx_npot_0_pf + 24) * nelems]);

    auto ta_z_xzz_0 = &(buffer.data()[(idx_npot_0_pf + 25) * nelems]);

    auto ta_z_yyy_0 = &(buffer.data()[(idx_npot_0_pf + 26) * nelems]);

    auto ta_z_yyz_0 = &(buffer.data()[(idx_npot_0_pf + 27) * nelems]);

    auto ta_z_yzz_0 = &(buffer.data()[(idx_npot_0_pf + 28) * nelems]);

    auto ta_z_zzz_0 = &(buffer.data()[(idx_npot_0_pf + 29) * nelems]);

    // Set up components of auxiliary buffer : PF

    auto ta_x_xxx_1 = &(buffer.data()[idx_npot_1_pf * nelems]);

    auto ta_x_xxy_1 = &(buffer.data()[(idx_npot_1_pf + 1) * nelems]);

    auto ta_x_xxz_1 = &(buffer.data()[(idx_npot_1_pf + 2) * nelems]);

    auto ta_x_xyy_1 = &(buffer.data()[(idx_npot_1_pf + 3) * nelems]);

    auto ta_x_xyz_1 = &(buffer.data()[(idx_npot_1_pf + 4) * nelems]);

    auto ta_x_xzz_1 = &(buffer.data()[(idx_npot_1_pf + 5) * nelems]);

    auto ta_x_yyy_1 = &(buffer.data()[(idx_npot_1_pf + 6) * nelems]);

    auto ta_x_yyz_1 = &(buffer.data()[(idx_npot_1_pf + 7) * nelems]);

    auto ta_x_yzz_1 = &(buffer.data()[(idx_npot_1_pf + 8) * nelems]);

    auto ta_x_zzz_1 = &(buffer.data()[(idx_npot_1_pf + 9) * nelems]);

    auto ta_y_xxx_1 = &(buffer.data()[(idx_npot_1_pf + 10) * nelems]);

    auto ta_y_xxy_1 = &(buffer.data()[(idx_npot_1_pf + 11) * nelems]);

    auto ta_y_xxz_1 = &(buffer.data()[(idx_npot_1_pf + 12) * nelems]);

    auto ta_y_xyy_1 = &(buffer.data()[(idx_npot_1_pf + 13) * nelems]);

    auto ta_y_xyz_1 = &(buffer.data()[(idx_npot_1_pf + 14) * nelems]);

    auto ta_y_xzz_1 = &(buffer.data()[(idx_npot_1_pf + 15) * nelems]);

    auto ta_y_yyy_1 = &(buffer.data()[(idx_npot_1_pf + 16) * nelems]);

    auto ta_y_yyz_1 = &(buffer.data()[(idx_npot_1_pf + 17) * nelems]);

    auto ta_y_yzz_1 = &(buffer.data()[(idx_npot_1_pf + 18) * nelems]);

    auto ta_y_zzz_1 = &(buffer.data()[(idx_npot_1_pf + 19) * nelems]);

    auto ta_z_xxx_1 = &(buffer.data()[(idx_npot_1_pf + 20) * nelems]);

    auto ta_z_xxy_1 = &(buffer.data()[(idx_npot_1_pf + 21) * nelems]);

    auto ta_z_xxz_1 = &(buffer.data()[(idx_npot_1_pf + 22) * nelems]);

    auto ta_z_xyy_1 = &(buffer.data()[(idx_npot_1_pf + 23) * nelems]);

    auto ta_z_xyz_1 = &(buffer.data()[(idx_npot_1_pf + 24) * nelems]);

    auto ta_z_xzz_1 = &(buffer.data()[(idx_npot_1_pf + 25) * nelems]);

    auto ta_z_yyy_1 = &(buffer.data()[(idx_npot_1_pf + 26) * nelems]);

    auto ta_z_yyz_1 = &(buffer.data()[(idx_npot_1_pf + 27) * nelems]);

    auto ta_z_yzz_1 = &(buffer.data()[(idx_npot_1_pf + 28) * nelems]);

    auto ta_z_zzz_1 = &(buffer.data()[(idx_npot_1_pf + 29) * nelems]);

    // Set up 0-10 components of targeted buffer : DF

    auto ta_xx_xxx_0 = &(buffer.data()[idx_npot_0_df * nelems]);

    auto ta_xx_xxy_0 = &(buffer.data()[(idx_npot_0_df + 1) * nelems]);

    auto ta_xx_xxz_0 = &(buffer.data()[(idx_npot_0_df + 2) * nelems]);

    auto ta_xx_xyy_0 = &(buffer.data()[(idx_npot_0_df + 3) * nelems]);

    auto ta_xx_xyz_0 = &(buffer.data()[(idx_npot_0_df + 4) * nelems]);

    auto ta_xx_xzz_0 = &(buffer.data()[(idx_npot_0_df + 5) * nelems]);

    auto ta_xx_yyy_0 = &(buffer.data()[(idx_npot_0_df + 6) * nelems]);

    auto ta_xx_yyz_0 = &(buffer.data()[(idx_npot_0_df + 7) * nelems]);

    auto ta_xx_yzz_0 = &(buffer.data()[(idx_npot_0_df + 8) * nelems]);

    auto ta_xx_zzz_0 = &(buffer.data()[(idx_npot_0_df + 9) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xx_xxx_0[i] = ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + 3.0 * ta_x_xx_0[i] * fe_0 - 3.0 * ta_x_xx_1[i] * fe_0 + ta_x_xxx_0[i] * rpa_x - ta_x_xxx_1[i] * pc_x[i];

        ta_xx_xxy_0[i] = ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + 2.0 * ta_x_xy_0[i] * fe_0 - 2.0 * ta_x_xy_1[i] * fe_0 + ta_x_xxy_0[i] * rpa_x - ta_x_xxy_1[i] * pc_x[i];

        ta_xx_xxz_0[i] = ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + 2.0 * ta_x_xz_0[i] * fe_0 - 2.0 * ta_x_xz_1[i] * fe_0 + ta_x_xxz_0[i] * rpa_x - ta_x_xxz_1[i] * pc_x[i];

        ta_xx_xyy_0[i] = ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + ta_x_yy_0[i] * fe_0 - ta_x_yy_1[i] * fe_0 + ta_x_xyy_0[i] * rpa_x - ta_x_xyy_1[i] * pc_x[i];

        ta_xx_xyz_0[i] = ta_0_xyz_0[i] * fe_0 - ta_0_xyz_1[i] * fe_0 + ta_x_yz_0[i] * fe_0 - ta_x_yz_1[i] * fe_0 + ta_x_xyz_0[i] * rpa_x - ta_x_xyz_1[i] * pc_x[i];

        ta_xx_xzz_0[i] = ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + ta_x_zz_0[i] * fe_0 - ta_x_zz_1[i] * fe_0 + ta_x_xzz_0[i] * rpa_x - ta_x_xzz_1[i] * pc_x[i];

        ta_xx_yyy_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_x_yyy_0[i] * rpa_x - ta_x_yyy_1[i] * pc_x[i];

        ta_xx_yyz_0[i] = ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + ta_x_yyz_0[i] * rpa_x - ta_x_yyz_1[i] * pc_x[i];

        ta_xx_yzz_0[i] = ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + ta_x_yzz_0[i] * rpa_x - ta_x_yzz_1[i] * pc_x[i];

        ta_xx_zzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_x_zzz_0[i] * rpa_x - ta_x_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto ta_xy_xxx_0 = &(buffer.data()[(idx_npot_0_df + 10) * nelems]);

    auto ta_xy_xxy_0 = &(buffer.data()[(idx_npot_0_df + 11) * nelems]);

    auto ta_xy_xxz_0 = &(buffer.data()[(idx_npot_0_df + 12) * nelems]);

    auto ta_xy_xyy_0 = &(buffer.data()[(idx_npot_0_df + 13) * nelems]);

    auto ta_xy_xyz_0 = &(buffer.data()[(idx_npot_0_df + 14) * nelems]);

    auto ta_xy_xzz_0 = &(buffer.data()[(idx_npot_0_df + 15) * nelems]);

    auto ta_xy_yyy_0 = &(buffer.data()[(idx_npot_0_df + 16) * nelems]);

    auto ta_xy_yyz_0 = &(buffer.data()[(idx_npot_0_df + 17) * nelems]);

    auto ta_xy_yzz_0 = &(buffer.data()[(idx_npot_0_df + 18) * nelems]);

    auto ta_xy_zzz_0 = &(buffer.data()[(idx_npot_0_df + 19) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xy_xxx_0[i] = ta_x_xxx_0[i] * rpa_y - ta_x_xxx_1[i] * pc_y[i];

        ta_xy_xxy_0[i] = 2.0 * ta_y_xy_0[i] * fe_0 - 2.0 * ta_y_xy_1[i] * fe_0 + ta_y_xxy_0[i] * rpa_x - ta_y_xxy_1[i] * pc_x[i];

        ta_xy_xxz_0[i] = ta_x_xxz_0[i] * rpa_y - ta_x_xxz_1[i] * pc_y[i];

        ta_xy_xyy_0[i] = ta_y_yy_0[i] * fe_0 - ta_y_yy_1[i] * fe_0 + ta_y_xyy_0[i] * rpa_x - ta_y_xyy_1[i] * pc_x[i];

        ta_xy_xyz_0[i] = ta_y_yz_0[i] * fe_0 - ta_y_yz_1[i] * fe_0 + ta_y_xyz_0[i] * rpa_x - ta_y_xyz_1[i] * pc_x[i];

        ta_xy_xzz_0[i] = ta_x_xzz_0[i] * rpa_y - ta_x_xzz_1[i] * pc_y[i];

        ta_xy_yyy_0[i] = ta_y_yyy_0[i] * rpa_x - ta_y_yyy_1[i] * pc_x[i];

        ta_xy_yyz_0[i] = ta_y_yyz_0[i] * rpa_x - ta_y_yyz_1[i] * pc_x[i];

        ta_xy_yzz_0[i] = ta_y_yzz_0[i] * rpa_x - ta_y_yzz_1[i] * pc_x[i];

        ta_xy_zzz_0[i] = ta_y_zzz_0[i] * rpa_x - ta_y_zzz_1[i] * pc_x[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto ta_xz_xxx_0 = &(buffer.data()[(idx_npot_0_df + 20) * nelems]);

    auto ta_xz_xxy_0 = &(buffer.data()[(idx_npot_0_df + 21) * nelems]);

    auto ta_xz_xxz_0 = &(buffer.data()[(idx_npot_0_df + 22) * nelems]);

    auto ta_xz_xyy_0 = &(buffer.data()[(idx_npot_0_df + 23) * nelems]);

    auto ta_xz_xyz_0 = &(buffer.data()[(idx_npot_0_df + 24) * nelems]);

    auto ta_xz_xzz_0 = &(buffer.data()[(idx_npot_0_df + 25) * nelems]);

    auto ta_xz_yyy_0 = &(buffer.data()[(idx_npot_0_df + 26) * nelems]);

    auto ta_xz_yyz_0 = &(buffer.data()[(idx_npot_0_df + 27) * nelems]);

    auto ta_xz_yzz_0 = &(buffer.data()[(idx_npot_0_df + 28) * nelems]);

    auto ta_xz_zzz_0 = &(buffer.data()[(idx_npot_0_df + 29) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xz_xxx_0[i] = ta_x_xxx_0[i] * rpa_z - ta_x_xxx_1[i] * pc_z[i];

        ta_xz_xxy_0[i] = ta_x_xxy_0[i] * rpa_z - ta_x_xxy_1[i] * pc_z[i];

        ta_xz_xxz_0[i] = 2.0 * ta_z_xz_0[i] * fe_0 - 2.0 * ta_z_xz_1[i] * fe_0 + ta_z_xxz_0[i] * rpa_x - ta_z_xxz_1[i] * pc_x[i];

        ta_xz_xyy_0[i] = ta_x_xyy_0[i] * rpa_z - ta_x_xyy_1[i] * pc_z[i];

        ta_xz_xyz_0[i] = ta_z_yz_0[i] * fe_0 - ta_z_yz_1[i] * fe_0 + ta_z_xyz_0[i] * rpa_x - ta_z_xyz_1[i] * pc_x[i];

        ta_xz_xzz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_z_xzz_0[i] * rpa_x - ta_z_xzz_1[i] * pc_x[i];

        ta_xz_yyy_0[i] = ta_z_yyy_0[i] * rpa_x - ta_z_yyy_1[i] * pc_x[i];

        ta_xz_yyz_0[i] = ta_z_yyz_0[i] * rpa_x - ta_z_yyz_1[i] * pc_x[i];

        ta_xz_yzz_0[i] = ta_z_yzz_0[i] * rpa_x - ta_z_yzz_1[i] * pc_x[i];

        ta_xz_zzz_0[i] = ta_z_zzz_0[i] * rpa_x - ta_z_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto ta_yy_xxx_0 = &(buffer.data()[(idx_npot_0_df + 30) * nelems]);

    auto ta_yy_xxy_0 = &(buffer.data()[(idx_npot_0_df + 31) * nelems]);

    auto ta_yy_xxz_0 = &(buffer.data()[(idx_npot_0_df + 32) * nelems]);

    auto ta_yy_xyy_0 = &(buffer.data()[(idx_npot_0_df + 33) * nelems]);

    auto ta_yy_xyz_0 = &(buffer.data()[(idx_npot_0_df + 34) * nelems]);

    auto ta_yy_xzz_0 = &(buffer.data()[(idx_npot_0_df + 35) * nelems]);

    auto ta_yy_yyy_0 = &(buffer.data()[(idx_npot_0_df + 36) * nelems]);

    auto ta_yy_yyz_0 = &(buffer.data()[(idx_npot_0_df + 37) * nelems]);

    auto ta_yy_yzz_0 = &(buffer.data()[(idx_npot_0_df + 38) * nelems]);

    auto ta_yy_zzz_0 = &(buffer.data()[(idx_npot_0_df + 39) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yy_xxx_0[i] = ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + ta_y_xxx_0[i] * rpa_y - ta_y_xxx_1[i] * pc_y[i];

        ta_yy_xxy_0[i] = ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + ta_y_xx_0[i] * fe_0 - ta_y_xx_1[i] * fe_0 + ta_y_xxy_0[i] * rpa_y - ta_y_xxy_1[i] * pc_y[i];

        ta_yy_xxz_0[i] = ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + ta_y_xxz_0[i] * rpa_y - ta_y_xxz_1[i] * pc_y[i];

        ta_yy_xyy_0[i] = ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + 2.0 * ta_y_xy_0[i] * fe_0 - 2.0 * ta_y_xy_1[i] * fe_0 + ta_y_xyy_0[i] * rpa_y - ta_y_xyy_1[i] * pc_y[i];

        ta_yy_xyz_0[i] = ta_0_xyz_0[i] * fe_0 - ta_0_xyz_1[i] * fe_0 + ta_y_xz_0[i] * fe_0 - ta_y_xz_1[i] * fe_0 + ta_y_xyz_0[i] * rpa_y - ta_y_xyz_1[i] * pc_y[i];

        ta_yy_xzz_0[i] = ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + ta_y_xzz_0[i] * rpa_y - ta_y_xzz_1[i] * pc_y[i];

        ta_yy_yyy_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + 3.0 * ta_y_yy_0[i] * fe_0 - 3.0 * ta_y_yy_1[i] * fe_0 + ta_y_yyy_0[i] * rpa_y - ta_y_yyy_1[i] * pc_y[i];

        ta_yy_yyz_0[i] = ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + 2.0 * ta_y_yz_0[i] * fe_0 - 2.0 * ta_y_yz_1[i] * fe_0 + ta_y_yyz_0[i] * rpa_y - ta_y_yyz_1[i] * pc_y[i];

        ta_yy_yzz_0[i] = ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + ta_y_zz_0[i] * fe_0 - ta_y_zz_1[i] * fe_0 + ta_y_yzz_0[i] * rpa_y - ta_y_yzz_1[i] * pc_y[i];

        ta_yy_zzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + ta_y_zzz_0[i] * rpa_y - ta_y_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto ta_yz_xxx_0 = &(buffer.data()[(idx_npot_0_df + 40) * nelems]);

    auto ta_yz_xxy_0 = &(buffer.data()[(idx_npot_0_df + 41) * nelems]);

    auto ta_yz_xxz_0 = &(buffer.data()[(idx_npot_0_df + 42) * nelems]);

    auto ta_yz_xyy_0 = &(buffer.data()[(idx_npot_0_df + 43) * nelems]);

    auto ta_yz_xyz_0 = &(buffer.data()[(idx_npot_0_df + 44) * nelems]);

    auto ta_yz_xzz_0 = &(buffer.data()[(idx_npot_0_df + 45) * nelems]);

    auto ta_yz_yyy_0 = &(buffer.data()[(idx_npot_0_df + 46) * nelems]);

    auto ta_yz_yyz_0 = &(buffer.data()[(idx_npot_0_df + 47) * nelems]);

    auto ta_yz_yzz_0 = &(buffer.data()[(idx_npot_0_df + 48) * nelems]);

    auto ta_yz_zzz_0 = &(buffer.data()[(idx_npot_0_df + 49) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yz_xxx_0[i] = ta_z_xxx_0[i] * rpa_y - ta_z_xxx_1[i] * pc_y[i];

        ta_yz_xxy_0[i] = ta_y_xxy_0[i] * rpa_z - ta_y_xxy_1[i] * pc_z[i];

        ta_yz_xxz_0[i] = ta_z_xxz_0[i] * rpa_y - ta_z_xxz_1[i] * pc_y[i];

        ta_yz_xyy_0[i] = ta_y_xyy_0[i] * rpa_z - ta_y_xyy_1[i] * pc_z[i];

        ta_yz_xyz_0[i] = ta_z_xz_0[i] * fe_0 - ta_z_xz_1[i] * fe_0 + ta_z_xyz_0[i] * rpa_y - ta_z_xyz_1[i] * pc_y[i];

        ta_yz_xzz_0[i] = ta_z_xzz_0[i] * rpa_y - ta_z_xzz_1[i] * pc_y[i];

        ta_yz_yyy_0[i] = ta_y_yyy_0[i] * rpa_z - ta_y_yyy_1[i] * pc_z[i];

        ta_yz_yyz_0[i] = 2.0 * ta_z_yz_0[i] * fe_0 - 2.0 * ta_z_yz_1[i] * fe_0 + ta_z_yyz_0[i] * rpa_y - ta_z_yyz_1[i] * pc_y[i];

        ta_yz_yzz_0[i] = ta_z_zz_0[i] * fe_0 - ta_z_zz_1[i] * fe_0 + ta_z_yzz_0[i] * rpa_y - ta_z_yzz_1[i] * pc_y[i];

        ta_yz_zzz_0[i] = ta_z_zzz_0[i] * rpa_y - ta_z_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto ta_zz_xxx_0 = &(buffer.data()[(idx_npot_0_df + 50) * nelems]);

    auto ta_zz_xxy_0 = &(buffer.data()[(idx_npot_0_df + 51) * nelems]);

    auto ta_zz_xxz_0 = &(buffer.data()[(idx_npot_0_df + 52) * nelems]);

    auto ta_zz_xyy_0 = &(buffer.data()[(idx_npot_0_df + 53) * nelems]);

    auto ta_zz_xyz_0 = &(buffer.data()[(idx_npot_0_df + 54) * nelems]);

    auto ta_zz_xzz_0 = &(buffer.data()[(idx_npot_0_df + 55) * nelems]);

    auto ta_zz_yyy_0 = &(buffer.data()[(idx_npot_0_df + 56) * nelems]);

    auto ta_zz_yyz_0 = &(buffer.data()[(idx_npot_0_df + 57) * nelems]);

    auto ta_zz_yzz_0 = &(buffer.data()[(idx_npot_0_df + 58) * nelems]);

    auto ta_zz_zzz_0 = &(buffer.data()[(idx_npot_0_df + 59) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_zz_xxx_0[i] = ta_0_xxx_0[i] * fe_0 - ta_0_xxx_1[i] * fe_0 + ta_z_xxx_0[i] * rpa_z - ta_z_xxx_1[i] * pc_z[i];

        ta_zz_xxy_0[i] = ta_0_xxy_0[i] * fe_0 - ta_0_xxy_1[i] * fe_0 + ta_z_xxy_0[i] * rpa_z - ta_z_xxy_1[i] * pc_z[i];

        ta_zz_xxz_0[i] = ta_0_xxz_0[i] * fe_0 - ta_0_xxz_1[i] * fe_0 + ta_z_xx_0[i] * fe_0 - ta_z_xx_1[i] * fe_0 + ta_z_xxz_0[i] * rpa_z - ta_z_xxz_1[i] * pc_z[i];

        ta_zz_xyy_0[i] = ta_0_xyy_0[i] * fe_0 - ta_0_xyy_1[i] * fe_0 + ta_z_xyy_0[i] * rpa_z - ta_z_xyy_1[i] * pc_z[i];

        ta_zz_xyz_0[i] = ta_0_xyz_0[i] * fe_0 - ta_0_xyz_1[i] * fe_0 + ta_z_xy_0[i] * fe_0 - ta_z_xy_1[i] * fe_0 + ta_z_xyz_0[i] * rpa_z - ta_z_xyz_1[i] * pc_z[i];

        ta_zz_xzz_0[i] = ta_0_xzz_0[i] * fe_0 - ta_0_xzz_1[i] * fe_0 + 2.0 * ta_z_xz_0[i] * fe_0 - 2.0 * ta_z_xz_1[i] * fe_0 + ta_z_xzz_0[i] * rpa_z - ta_z_xzz_1[i] * pc_z[i];

        ta_zz_yyy_0[i] = ta_0_yyy_0[i] * fe_0 - ta_0_yyy_1[i] * fe_0 + ta_z_yyy_0[i] * rpa_z - ta_z_yyy_1[i] * pc_z[i];

        ta_zz_yyz_0[i] = ta_0_yyz_0[i] * fe_0 - ta_0_yyz_1[i] * fe_0 + ta_z_yy_0[i] * fe_0 - ta_z_yy_1[i] * fe_0 + ta_z_yyz_0[i] * rpa_z - ta_z_yyz_1[i] * pc_z[i];

        ta_zz_yzz_0[i] = ta_0_yzz_0[i] * fe_0 - ta_0_yzz_1[i] * fe_0 + 2.0 * ta_z_yz_0[i] * fe_0 - 2.0 * ta_z_yz_1[i] * fe_0 + ta_z_yzz_0[i] * rpa_z - ta_z_yzz_1[i] * pc_z[i];

        ta_zz_zzz_0[i] = ta_0_zzz_0[i] * fe_0 - ta_0_zzz_1[i] * fe_0 + 3.0 * ta_z_zz_0[i] * fe_0 - 3.0 * ta_z_zz_1[i] * fe_0 + ta_z_zzz_0[i] * rpa_z - ta_z_zzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

