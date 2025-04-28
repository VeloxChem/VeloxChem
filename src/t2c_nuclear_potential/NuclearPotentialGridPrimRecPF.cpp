#include "NuclearPotentialGridPrimRecPF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_pf(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_pf,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_1_sd,
                                       const size_t idx_npot_0_sf,
                                       const size_t idx_npot_1_sf,
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

    // Set up 0-10 components of targeted buffer : PF

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

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_x_xxx_0[i] = 3.0 * ta_0_xx_0[i] * fe_0 - 3.0 * ta_0_xx_1[i] * fe_0 + ta_0_xxx_0[i] * rpa_x - ta_0_xxx_1[i] * pc_x[i];

        ta_x_xxy_0[i] = 2.0 * ta_0_xy_0[i] * fe_0 - 2.0 * ta_0_xy_1[i] * fe_0 + ta_0_xxy_0[i] * rpa_x - ta_0_xxy_1[i] * pc_x[i];

        ta_x_xxz_0[i] = 2.0 * ta_0_xz_0[i] * fe_0 - 2.0 * ta_0_xz_1[i] * fe_0 + ta_0_xxz_0[i] * rpa_x - ta_0_xxz_1[i] * pc_x[i];

        ta_x_xyy_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_0_xyy_0[i] * rpa_x - ta_0_xyy_1[i] * pc_x[i];

        ta_x_xyz_0[i] = ta_0_yz_0[i] * fe_0 - ta_0_yz_1[i] * fe_0 + ta_0_xyz_0[i] * rpa_x - ta_0_xyz_1[i] * pc_x[i];

        ta_x_xzz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_0_xzz_0[i] * rpa_x - ta_0_xzz_1[i] * pc_x[i];

        ta_x_yyy_0[i] = ta_0_yyy_0[i] * rpa_x - ta_0_yyy_1[i] * pc_x[i];

        ta_x_yyz_0[i] = ta_0_yyz_0[i] * rpa_x - ta_0_yyz_1[i] * pc_x[i];

        ta_x_yzz_0[i] = ta_0_yzz_0[i] * rpa_x - ta_0_yzz_1[i] * pc_x[i];

        ta_x_zzz_0[i] = ta_0_zzz_0[i] * rpa_x - ta_0_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

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

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_y_xxx_0[i] = ta_0_xxx_0[i] * rpa_y - ta_0_xxx_1[i] * pc_y[i];

        ta_y_xxy_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_0_xxy_0[i] * rpa_y - ta_0_xxy_1[i] * pc_y[i];

        ta_y_xxz_0[i] = ta_0_xxz_0[i] * rpa_y - ta_0_xxz_1[i] * pc_y[i];

        ta_y_xyy_0[i] = 2.0 * ta_0_xy_0[i] * fe_0 - 2.0 * ta_0_xy_1[i] * fe_0 + ta_0_xyy_0[i] * rpa_y - ta_0_xyy_1[i] * pc_y[i];

        ta_y_xyz_0[i] = ta_0_xz_0[i] * fe_0 - ta_0_xz_1[i] * fe_0 + ta_0_xyz_0[i] * rpa_y - ta_0_xyz_1[i] * pc_y[i];

        ta_y_xzz_0[i] = ta_0_xzz_0[i] * rpa_y - ta_0_xzz_1[i] * pc_y[i];

        ta_y_yyy_0[i] = 3.0 * ta_0_yy_0[i] * fe_0 - 3.0 * ta_0_yy_1[i] * fe_0 + ta_0_yyy_0[i] * rpa_y - ta_0_yyy_1[i] * pc_y[i];

        ta_y_yyz_0[i] = 2.0 * ta_0_yz_0[i] * fe_0 - 2.0 * ta_0_yz_1[i] * fe_0 + ta_0_yyz_0[i] * rpa_y - ta_0_yyz_1[i] * pc_y[i];

        ta_y_yzz_0[i] = ta_0_zz_0[i] * fe_0 - ta_0_zz_1[i] * fe_0 + ta_0_yzz_0[i] * rpa_y - ta_0_yzz_1[i] * pc_y[i];

        ta_y_zzz_0[i] = ta_0_zzz_0[i] * rpa_y - ta_0_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : PF

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

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_z_xxx_0[i] = ta_0_xxx_0[i] * rpa_z - ta_0_xxx_1[i] * pc_z[i];

        ta_z_xxy_0[i] = ta_0_xxy_0[i] * rpa_z - ta_0_xxy_1[i] * pc_z[i];

        ta_z_xxz_0[i] = ta_0_xx_0[i] * fe_0 - ta_0_xx_1[i] * fe_0 + ta_0_xxz_0[i] * rpa_z - ta_0_xxz_1[i] * pc_z[i];

        ta_z_xyy_0[i] = ta_0_xyy_0[i] * rpa_z - ta_0_xyy_1[i] * pc_z[i];

        ta_z_xyz_0[i] = ta_0_xy_0[i] * fe_0 - ta_0_xy_1[i] * fe_0 + ta_0_xyz_0[i] * rpa_z - ta_0_xyz_1[i] * pc_z[i];

        ta_z_xzz_0[i] = 2.0 * ta_0_xz_0[i] * fe_0 - 2.0 * ta_0_xz_1[i] * fe_0 + ta_0_xzz_0[i] * rpa_z - ta_0_xzz_1[i] * pc_z[i];

        ta_z_yyy_0[i] = ta_0_yyy_0[i] * rpa_z - ta_0_yyy_1[i] * pc_z[i];

        ta_z_yyz_0[i] = ta_0_yy_0[i] * fe_0 - ta_0_yy_1[i] * fe_0 + ta_0_yyz_0[i] * rpa_z - ta_0_yyz_1[i] * pc_z[i];

        ta_z_yzz_0[i] = 2.0 * ta_0_yz_0[i] * fe_0 - 2.0 * ta_0_yz_1[i] * fe_0 + ta_0_yzz_0[i] * rpa_z - ta_0_yzz_1[i] * pc_z[i];

        ta_z_zzz_0[i] = 3.0 * ta_0_zz_0[i] * fe_0 - 3.0 * ta_0_zz_1[i] * fe_0 + ta_0_zzz_0[i] * rpa_z - ta_0_zzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

