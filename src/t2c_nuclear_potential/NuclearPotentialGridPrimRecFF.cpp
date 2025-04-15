#include "NuclearPotentialGridPrimRecFF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_ff(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_ff,
                                       const size_t idx_npot_0_pf,
                                       const size_t idx_npot_1_pf,
                                       const size_t idx_npot_0_dd,
                                       const size_t idx_npot_1_dd,
                                       const size_t idx_npot_0_df,
                                       const size_t idx_npot_1_df,
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

    // Set up components of auxiliary buffer : DD

    auto ta_xx_xx_0 = &(buffer.data()[idx_npot_0_dd * nelems]);

    auto ta_xx_xy_0 = &(buffer.data()[(idx_npot_0_dd + 1) * nelems]);

    auto ta_xx_xz_0 = &(buffer.data()[(idx_npot_0_dd + 2) * nelems]);

    auto ta_xx_yy_0 = &(buffer.data()[(idx_npot_0_dd + 3) * nelems]);

    auto ta_xx_yz_0 = &(buffer.data()[(idx_npot_0_dd + 4) * nelems]);

    auto ta_xx_zz_0 = &(buffer.data()[(idx_npot_0_dd + 5) * nelems]);

    auto ta_yy_xx_0 = &(buffer.data()[(idx_npot_0_dd + 18) * nelems]);

    auto ta_yy_xy_0 = &(buffer.data()[(idx_npot_0_dd + 19) * nelems]);

    auto ta_yy_xz_0 = &(buffer.data()[(idx_npot_0_dd + 20) * nelems]);

    auto ta_yy_yy_0 = &(buffer.data()[(idx_npot_0_dd + 21) * nelems]);

    auto ta_yy_yz_0 = &(buffer.data()[(idx_npot_0_dd + 22) * nelems]);

    auto ta_yy_zz_0 = &(buffer.data()[(idx_npot_0_dd + 23) * nelems]);

    auto ta_yz_yz_0 = &(buffer.data()[(idx_npot_0_dd + 28) * nelems]);

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

    auto ta_yy_xx_1 = &(buffer.data()[(idx_npot_1_dd + 18) * nelems]);

    auto ta_yy_xy_1 = &(buffer.data()[(idx_npot_1_dd + 19) * nelems]);

    auto ta_yy_xz_1 = &(buffer.data()[(idx_npot_1_dd + 20) * nelems]);

    auto ta_yy_yy_1 = &(buffer.data()[(idx_npot_1_dd + 21) * nelems]);

    auto ta_yy_yz_1 = &(buffer.data()[(idx_npot_1_dd + 22) * nelems]);

    auto ta_yy_zz_1 = &(buffer.data()[(idx_npot_1_dd + 23) * nelems]);

    auto ta_yz_yz_1 = &(buffer.data()[(idx_npot_1_dd + 28) * nelems]);

    auto ta_zz_xx_1 = &(buffer.data()[(idx_npot_1_dd + 30) * nelems]);

    auto ta_zz_xy_1 = &(buffer.data()[(idx_npot_1_dd + 31) * nelems]);

    auto ta_zz_xz_1 = &(buffer.data()[(idx_npot_1_dd + 32) * nelems]);

    auto ta_zz_yy_1 = &(buffer.data()[(idx_npot_1_dd + 33) * nelems]);

    auto ta_zz_yz_1 = &(buffer.data()[(idx_npot_1_dd + 34) * nelems]);

    auto ta_zz_zz_1 = &(buffer.data()[(idx_npot_1_dd + 35) * nelems]);

    // Set up components of auxiliary buffer : DF

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

    auto ta_xy_xxy_0 = &(buffer.data()[(idx_npot_0_df + 11) * nelems]);

    auto ta_xy_xyy_0 = &(buffer.data()[(idx_npot_0_df + 13) * nelems]);

    auto ta_xy_yyy_0 = &(buffer.data()[(idx_npot_0_df + 16) * nelems]);

    auto ta_xy_yyz_0 = &(buffer.data()[(idx_npot_0_df + 17) * nelems]);

    auto ta_xy_yzz_0 = &(buffer.data()[(idx_npot_0_df + 18) * nelems]);

    auto ta_xz_xxx_0 = &(buffer.data()[(idx_npot_0_df + 20) * nelems]);

    auto ta_xz_xxz_0 = &(buffer.data()[(idx_npot_0_df + 22) * nelems]);

    auto ta_xz_xzz_0 = &(buffer.data()[(idx_npot_0_df + 25) * nelems]);

    auto ta_xz_yyz_0 = &(buffer.data()[(idx_npot_0_df + 27) * nelems]);

    auto ta_xz_yzz_0 = &(buffer.data()[(idx_npot_0_df + 28) * nelems]);

    auto ta_xz_zzz_0 = &(buffer.data()[(idx_npot_0_df + 29) * nelems]);

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

    auto ta_yz_xxz_0 = &(buffer.data()[(idx_npot_0_df + 42) * nelems]);

    auto ta_yz_xyz_0 = &(buffer.data()[(idx_npot_0_df + 44) * nelems]);

    auto ta_yz_xzz_0 = &(buffer.data()[(idx_npot_0_df + 45) * nelems]);

    auto ta_yz_yyy_0 = &(buffer.data()[(idx_npot_0_df + 46) * nelems]);

    auto ta_yz_yyz_0 = &(buffer.data()[(idx_npot_0_df + 47) * nelems]);

    auto ta_yz_yzz_0 = &(buffer.data()[(idx_npot_0_df + 48) * nelems]);

    auto ta_yz_zzz_0 = &(buffer.data()[(idx_npot_0_df + 49) * nelems]);

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

    // Set up components of auxiliary buffer : DF

    auto ta_xx_xxx_1 = &(buffer.data()[idx_npot_1_df * nelems]);

    auto ta_xx_xxy_1 = &(buffer.data()[(idx_npot_1_df + 1) * nelems]);

    auto ta_xx_xxz_1 = &(buffer.data()[(idx_npot_1_df + 2) * nelems]);

    auto ta_xx_xyy_1 = &(buffer.data()[(idx_npot_1_df + 3) * nelems]);

    auto ta_xx_xyz_1 = &(buffer.data()[(idx_npot_1_df + 4) * nelems]);

    auto ta_xx_xzz_1 = &(buffer.data()[(idx_npot_1_df + 5) * nelems]);

    auto ta_xx_yyy_1 = &(buffer.data()[(idx_npot_1_df + 6) * nelems]);

    auto ta_xx_yyz_1 = &(buffer.data()[(idx_npot_1_df + 7) * nelems]);

    auto ta_xx_yzz_1 = &(buffer.data()[(idx_npot_1_df + 8) * nelems]);

    auto ta_xx_zzz_1 = &(buffer.data()[(idx_npot_1_df + 9) * nelems]);

    auto ta_xy_xxy_1 = &(buffer.data()[(idx_npot_1_df + 11) * nelems]);

    auto ta_xy_xyy_1 = &(buffer.data()[(idx_npot_1_df + 13) * nelems]);

    auto ta_xy_yyy_1 = &(buffer.data()[(idx_npot_1_df + 16) * nelems]);

    auto ta_xy_yyz_1 = &(buffer.data()[(idx_npot_1_df + 17) * nelems]);

    auto ta_xy_yzz_1 = &(buffer.data()[(idx_npot_1_df + 18) * nelems]);

    auto ta_xz_xxx_1 = &(buffer.data()[(idx_npot_1_df + 20) * nelems]);

    auto ta_xz_xxz_1 = &(buffer.data()[(idx_npot_1_df + 22) * nelems]);

    auto ta_xz_xzz_1 = &(buffer.data()[(idx_npot_1_df + 25) * nelems]);

    auto ta_xz_yyz_1 = &(buffer.data()[(idx_npot_1_df + 27) * nelems]);

    auto ta_xz_yzz_1 = &(buffer.data()[(idx_npot_1_df + 28) * nelems]);

    auto ta_xz_zzz_1 = &(buffer.data()[(idx_npot_1_df + 29) * nelems]);

    auto ta_yy_xxx_1 = &(buffer.data()[(idx_npot_1_df + 30) * nelems]);

    auto ta_yy_xxy_1 = &(buffer.data()[(idx_npot_1_df + 31) * nelems]);

    auto ta_yy_xxz_1 = &(buffer.data()[(idx_npot_1_df + 32) * nelems]);

    auto ta_yy_xyy_1 = &(buffer.data()[(idx_npot_1_df + 33) * nelems]);

    auto ta_yy_xyz_1 = &(buffer.data()[(idx_npot_1_df + 34) * nelems]);

    auto ta_yy_xzz_1 = &(buffer.data()[(idx_npot_1_df + 35) * nelems]);

    auto ta_yy_yyy_1 = &(buffer.data()[(idx_npot_1_df + 36) * nelems]);

    auto ta_yy_yyz_1 = &(buffer.data()[(idx_npot_1_df + 37) * nelems]);

    auto ta_yy_yzz_1 = &(buffer.data()[(idx_npot_1_df + 38) * nelems]);

    auto ta_yy_zzz_1 = &(buffer.data()[(idx_npot_1_df + 39) * nelems]);

    auto ta_yz_xxz_1 = &(buffer.data()[(idx_npot_1_df + 42) * nelems]);

    auto ta_yz_xyz_1 = &(buffer.data()[(idx_npot_1_df + 44) * nelems]);

    auto ta_yz_xzz_1 = &(buffer.data()[(idx_npot_1_df + 45) * nelems]);

    auto ta_yz_yyy_1 = &(buffer.data()[(idx_npot_1_df + 46) * nelems]);

    auto ta_yz_yyz_1 = &(buffer.data()[(idx_npot_1_df + 47) * nelems]);

    auto ta_yz_yzz_1 = &(buffer.data()[(idx_npot_1_df + 48) * nelems]);

    auto ta_yz_zzz_1 = &(buffer.data()[(idx_npot_1_df + 49) * nelems]);

    auto ta_zz_xxx_1 = &(buffer.data()[(idx_npot_1_df + 50) * nelems]);

    auto ta_zz_xxy_1 = &(buffer.data()[(idx_npot_1_df + 51) * nelems]);

    auto ta_zz_xxz_1 = &(buffer.data()[(idx_npot_1_df + 52) * nelems]);

    auto ta_zz_xyy_1 = &(buffer.data()[(idx_npot_1_df + 53) * nelems]);

    auto ta_zz_xyz_1 = &(buffer.data()[(idx_npot_1_df + 54) * nelems]);

    auto ta_zz_xzz_1 = &(buffer.data()[(idx_npot_1_df + 55) * nelems]);

    auto ta_zz_yyy_1 = &(buffer.data()[(idx_npot_1_df + 56) * nelems]);

    auto ta_zz_yyz_1 = &(buffer.data()[(idx_npot_1_df + 57) * nelems]);

    auto ta_zz_yzz_1 = &(buffer.data()[(idx_npot_1_df + 58) * nelems]);

    auto ta_zz_zzz_1 = &(buffer.data()[(idx_npot_1_df + 59) * nelems]);

    // Set up 0-10 components of targeted buffer : FF

    auto ta_xxx_xxx_0 = &(buffer.data()[idx_npot_0_ff * nelems]);

    auto ta_xxx_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 1) * nelems]);

    auto ta_xxx_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 2) * nelems]);

    auto ta_xxx_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 3) * nelems]);

    auto ta_xxx_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 4) * nelems]);

    auto ta_xxx_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 5) * nelems]);

    auto ta_xxx_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 6) * nelems]);

    auto ta_xxx_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 7) * nelems]);

    auto ta_xxx_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 8) * nelems]);

    auto ta_xxx_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 9) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxx_xxx_0[i] = 2.0 * ta_x_xxx_0[i] * fe_0 - 2.0 * ta_x_xxx_1[i] * fe_0 + 3.0 * ta_xx_xx_0[i] * fe_0 - 3.0 * ta_xx_xx_1[i] * fe_0 + ta_xx_xxx_0[i] * rpa_x - ta_xx_xxx_1[i] * pc_x[i];

        ta_xxx_xxy_0[i] = 2.0 * ta_x_xxy_0[i] * fe_0 - 2.0 * ta_x_xxy_1[i] * fe_0 + 2.0 * ta_xx_xy_0[i] * fe_0 - 2.0 * ta_xx_xy_1[i] * fe_0 + ta_xx_xxy_0[i] * rpa_x - ta_xx_xxy_1[i] * pc_x[i];

        ta_xxx_xxz_0[i] = 2.0 * ta_x_xxz_0[i] * fe_0 - 2.0 * ta_x_xxz_1[i] * fe_0 + 2.0 * ta_xx_xz_0[i] * fe_0 - 2.0 * ta_xx_xz_1[i] * fe_0 + ta_xx_xxz_0[i] * rpa_x - ta_xx_xxz_1[i] * pc_x[i];

        ta_xxx_xyy_0[i] = 2.0 * ta_x_xyy_0[i] * fe_0 - 2.0 * ta_x_xyy_1[i] * fe_0 + ta_xx_yy_0[i] * fe_0 - ta_xx_yy_1[i] * fe_0 + ta_xx_xyy_0[i] * rpa_x - ta_xx_xyy_1[i] * pc_x[i];

        ta_xxx_xyz_0[i] = 2.0 * ta_x_xyz_0[i] * fe_0 - 2.0 * ta_x_xyz_1[i] * fe_0 + ta_xx_yz_0[i] * fe_0 - ta_xx_yz_1[i] * fe_0 + ta_xx_xyz_0[i] * rpa_x - ta_xx_xyz_1[i] * pc_x[i];

        ta_xxx_xzz_0[i] = 2.0 * ta_x_xzz_0[i] * fe_0 - 2.0 * ta_x_xzz_1[i] * fe_0 + ta_xx_zz_0[i] * fe_0 - ta_xx_zz_1[i] * fe_0 + ta_xx_xzz_0[i] * rpa_x - ta_xx_xzz_1[i] * pc_x[i];

        ta_xxx_yyy_0[i] = 2.0 * ta_x_yyy_0[i] * fe_0 - 2.0 * ta_x_yyy_1[i] * fe_0 + ta_xx_yyy_0[i] * rpa_x - ta_xx_yyy_1[i] * pc_x[i];

        ta_xxx_yyz_0[i] = 2.0 * ta_x_yyz_0[i] * fe_0 - 2.0 * ta_x_yyz_1[i] * fe_0 + ta_xx_yyz_0[i] * rpa_x - ta_xx_yyz_1[i] * pc_x[i];

        ta_xxx_yzz_0[i] = 2.0 * ta_x_yzz_0[i] * fe_0 - 2.0 * ta_x_yzz_1[i] * fe_0 + ta_xx_yzz_0[i] * rpa_x - ta_xx_yzz_1[i] * pc_x[i];

        ta_xxx_zzz_0[i] = 2.0 * ta_x_zzz_0[i] * fe_0 - 2.0 * ta_x_zzz_1[i] * fe_0 + ta_xx_zzz_0[i] * rpa_x - ta_xx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto ta_xxy_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 10) * nelems]);

    auto ta_xxy_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 11) * nelems]);

    auto ta_xxy_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 12) * nelems]);

    auto ta_xxy_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 13) * nelems]);

    auto ta_xxy_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 14) * nelems]);

    auto ta_xxy_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 15) * nelems]);

    auto ta_xxy_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 16) * nelems]);

    auto ta_xxy_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 17) * nelems]);

    auto ta_xxy_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 18) * nelems]);

    auto ta_xxy_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 19) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxy_xxx_0[i] = ta_xx_xxx_0[i] * rpa_y - ta_xx_xxx_1[i] * pc_y[i];

        ta_xxy_xxy_0[i] = ta_xx_xx_0[i] * fe_0 - ta_xx_xx_1[i] * fe_0 + ta_xx_xxy_0[i] * rpa_y - ta_xx_xxy_1[i] * pc_y[i];

        ta_xxy_xxz_0[i] = ta_xx_xxz_0[i] * rpa_y - ta_xx_xxz_1[i] * pc_y[i];

        ta_xxy_xyy_0[i] = 2.0 * ta_xx_xy_0[i] * fe_0 - 2.0 * ta_xx_xy_1[i] * fe_0 + ta_xx_xyy_0[i] * rpa_y - ta_xx_xyy_1[i] * pc_y[i];

        ta_xxy_xyz_0[i] = ta_xx_xz_0[i] * fe_0 - ta_xx_xz_1[i] * fe_0 + ta_xx_xyz_0[i] * rpa_y - ta_xx_xyz_1[i] * pc_y[i];

        ta_xxy_xzz_0[i] = ta_xx_xzz_0[i] * rpa_y - ta_xx_xzz_1[i] * pc_y[i];

        ta_xxy_yyy_0[i] = ta_y_yyy_0[i] * fe_0 - ta_y_yyy_1[i] * fe_0 + ta_xy_yyy_0[i] * rpa_x - ta_xy_yyy_1[i] * pc_x[i];

        ta_xxy_yyz_0[i] = ta_y_yyz_0[i] * fe_0 - ta_y_yyz_1[i] * fe_0 + ta_xy_yyz_0[i] * rpa_x - ta_xy_yyz_1[i] * pc_x[i];

        ta_xxy_yzz_0[i] = ta_y_yzz_0[i] * fe_0 - ta_y_yzz_1[i] * fe_0 + ta_xy_yzz_0[i] * rpa_x - ta_xy_yzz_1[i] * pc_x[i];

        ta_xxy_zzz_0[i] = ta_xx_zzz_0[i] * rpa_y - ta_xx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto ta_xxz_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 20) * nelems]);

    auto ta_xxz_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 21) * nelems]);

    auto ta_xxz_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 22) * nelems]);

    auto ta_xxz_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 23) * nelems]);

    auto ta_xxz_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 24) * nelems]);

    auto ta_xxz_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 25) * nelems]);

    auto ta_xxz_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 26) * nelems]);

    auto ta_xxz_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 27) * nelems]);

    auto ta_xxz_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 28) * nelems]);

    auto ta_xxz_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 29) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxz_xxx_0[i] = ta_xx_xxx_0[i] * rpa_z - ta_xx_xxx_1[i] * pc_z[i];

        ta_xxz_xxy_0[i] = ta_xx_xxy_0[i] * rpa_z - ta_xx_xxy_1[i] * pc_z[i];

        ta_xxz_xxz_0[i] = ta_xx_xx_0[i] * fe_0 - ta_xx_xx_1[i] * fe_0 + ta_xx_xxz_0[i] * rpa_z - ta_xx_xxz_1[i] * pc_z[i];

        ta_xxz_xyy_0[i] = ta_xx_xyy_0[i] * rpa_z - ta_xx_xyy_1[i] * pc_z[i];

        ta_xxz_xyz_0[i] = ta_xx_xy_0[i] * fe_0 - ta_xx_xy_1[i] * fe_0 + ta_xx_xyz_0[i] * rpa_z - ta_xx_xyz_1[i] * pc_z[i];

        ta_xxz_xzz_0[i] = 2.0 * ta_xx_xz_0[i] * fe_0 - 2.0 * ta_xx_xz_1[i] * fe_0 + ta_xx_xzz_0[i] * rpa_z - ta_xx_xzz_1[i] * pc_z[i];

        ta_xxz_yyy_0[i] = ta_xx_yyy_0[i] * rpa_z - ta_xx_yyy_1[i] * pc_z[i];

        ta_xxz_yyz_0[i] = ta_z_yyz_0[i] * fe_0 - ta_z_yyz_1[i] * fe_0 + ta_xz_yyz_0[i] * rpa_x - ta_xz_yyz_1[i] * pc_x[i];

        ta_xxz_yzz_0[i] = ta_z_yzz_0[i] * fe_0 - ta_z_yzz_1[i] * fe_0 + ta_xz_yzz_0[i] * rpa_x - ta_xz_yzz_1[i] * pc_x[i];

        ta_xxz_zzz_0[i] = ta_z_zzz_0[i] * fe_0 - ta_z_zzz_1[i] * fe_0 + ta_xz_zzz_0[i] * rpa_x - ta_xz_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto ta_xyy_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 30) * nelems]);

    auto ta_xyy_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 31) * nelems]);

    auto ta_xyy_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 32) * nelems]);

    auto ta_xyy_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 33) * nelems]);

    auto ta_xyy_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 34) * nelems]);

    auto ta_xyy_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 35) * nelems]);

    auto ta_xyy_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 36) * nelems]);

    auto ta_xyy_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 37) * nelems]);

    auto ta_xyy_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 38) * nelems]);

    auto ta_xyy_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 39) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xyy_xxx_0[i] = 3.0 * ta_yy_xx_0[i] * fe_0 - 3.0 * ta_yy_xx_1[i] * fe_0 + ta_yy_xxx_0[i] * rpa_x - ta_yy_xxx_1[i] * pc_x[i];

        ta_xyy_xxy_0[i] = 2.0 * ta_yy_xy_0[i] * fe_0 - 2.0 * ta_yy_xy_1[i] * fe_0 + ta_yy_xxy_0[i] * rpa_x - ta_yy_xxy_1[i] * pc_x[i];

        ta_xyy_xxz_0[i] = 2.0 * ta_yy_xz_0[i] * fe_0 - 2.0 * ta_yy_xz_1[i] * fe_0 + ta_yy_xxz_0[i] * rpa_x - ta_yy_xxz_1[i] * pc_x[i];

        ta_xyy_xyy_0[i] = ta_yy_yy_0[i] * fe_0 - ta_yy_yy_1[i] * fe_0 + ta_yy_xyy_0[i] * rpa_x - ta_yy_xyy_1[i] * pc_x[i];

        ta_xyy_xyz_0[i] = ta_yy_yz_0[i] * fe_0 - ta_yy_yz_1[i] * fe_0 + ta_yy_xyz_0[i] * rpa_x - ta_yy_xyz_1[i] * pc_x[i];

        ta_xyy_xzz_0[i] = ta_yy_zz_0[i] * fe_0 - ta_yy_zz_1[i] * fe_0 + ta_yy_xzz_0[i] * rpa_x - ta_yy_xzz_1[i] * pc_x[i];

        ta_xyy_yyy_0[i] = ta_yy_yyy_0[i] * rpa_x - ta_yy_yyy_1[i] * pc_x[i];

        ta_xyy_yyz_0[i] = ta_yy_yyz_0[i] * rpa_x - ta_yy_yyz_1[i] * pc_x[i];

        ta_xyy_yzz_0[i] = ta_yy_yzz_0[i] * rpa_x - ta_yy_yzz_1[i] * pc_x[i];

        ta_xyy_zzz_0[i] = ta_yy_zzz_0[i] * rpa_x - ta_yy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto ta_xyz_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 40) * nelems]);

    auto ta_xyz_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 41) * nelems]);

    auto ta_xyz_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 42) * nelems]);

    auto ta_xyz_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 43) * nelems]);

    auto ta_xyz_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 44) * nelems]);

    auto ta_xyz_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 45) * nelems]);

    auto ta_xyz_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 46) * nelems]);

    auto ta_xyz_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 47) * nelems]);

    auto ta_xyz_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 48) * nelems]);

    auto ta_xyz_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 49) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xyz_xxx_0[i] = ta_xz_xxx_0[i] * rpa_y - ta_xz_xxx_1[i] * pc_y[i];

        ta_xyz_xxy_0[i] = ta_xy_xxy_0[i] * rpa_z - ta_xy_xxy_1[i] * pc_z[i];

        ta_xyz_xxz_0[i] = ta_xz_xxz_0[i] * rpa_y - ta_xz_xxz_1[i] * pc_y[i];

        ta_xyz_xyy_0[i] = ta_xy_xyy_0[i] * rpa_z - ta_xy_xyy_1[i] * pc_z[i];

        ta_xyz_xyz_0[i] = ta_yz_yz_0[i] * fe_0 - ta_yz_yz_1[i] * fe_0 + ta_yz_xyz_0[i] * rpa_x - ta_yz_xyz_1[i] * pc_x[i];

        ta_xyz_xzz_0[i] = ta_xz_xzz_0[i] * rpa_y - ta_xz_xzz_1[i] * pc_y[i];

        ta_xyz_yyy_0[i] = ta_yz_yyy_0[i] * rpa_x - ta_yz_yyy_1[i] * pc_x[i];

        ta_xyz_yyz_0[i] = ta_yz_yyz_0[i] * rpa_x - ta_yz_yyz_1[i] * pc_x[i];

        ta_xyz_yzz_0[i] = ta_yz_yzz_0[i] * rpa_x - ta_yz_yzz_1[i] * pc_x[i];

        ta_xyz_zzz_0[i] = ta_yz_zzz_0[i] * rpa_x - ta_yz_zzz_1[i] * pc_x[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto ta_xzz_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 50) * nelems]);

    auto ta_xzz_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 51) * nelems]);

    auto ta_xzz_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 52) * nelems]);

    auto ta_xzz_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 53) * nelems]);

    auto ta_xzz_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 54) * nelems]);

    auto ta_xzz_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 55) * nelems]);

    auto ta_xzz_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 56) * nelems]);

    auto ta_xzz_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 57) * nelems]);

    auto ta_xzz_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 58) * nelems]);

    auto ta_xzz_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 59) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xzz_xxx_0[i] = 3.0 * ta_zz_xx_0[i] * fe_0 - 3.0 * ta_zz_xx_1[i] * fe_0 + ta_zz_xxx_0[i] * rpa_x - ta_zz_xxx_1[i] * pc_x[i];

        ta_xzz_xxy_0[i] = 2.0 * ta_zz_xy_0[i] * fe_0 - 2.0 * ta_zz_xy_1[i] * fe_0 + ta_zz_xxy_0[i] * rpa_x - ta_zz_xxy_1[i] * pc_x[i];

        ta_xzz_xxz_0[i] = 2.0 * ta_zz_xz_0[i] * fe_0 - 2.0 * ta_zz_xz_1[i] * fe_0 + ta_zz_xxz_0[i] * rpa_x - ta_zz_xxz_1[i] * pc_x[i];

        ta_xzz_xyy_0[i] = ta_zz_yy_0[i] * fe_0 - ta_zz_yy_1[i] * fe_0 + ta_zz_xyy_0[i] * rpa_x - ta_zz_xyy_1[i] * pc_x[i];

        ta_xzz_xyz_0[i] = ta_zz_yz_0[i] * fe_0 - ta_zz_yz_1[i] * fe_0 + ta_zz_xyz_0[i] * rpa_x - ta_zz_xyz_1[i] * pc_x[i];

        ta_xzz_xzz_0[i] = ta_zz_zz_0[i] * fe_0 - ta_zz_zz_1[i] * fe_0 + ta_zz_xzz_0[i] * rpa_x - ta_zz_xzz_1[i] * pc_x[i];

        ta_xzz_yyy_0[i] = ta_zz_yyy_0[i] * rpa_x - ta_zz_yyy_1[i] * pc_x[i];

        ta_xzz_yyz_0[i] = ta_zz_yyz_0[i] * rpa_x - ta_zz_yyz_1[i] * pc_x[i];

        ta_xzz_yzz_0[i] = ta_zz_yzz_0[i] * rpa_x - ta_zz_yzz_1[i] * pc_x[i];

        ta_xzz_zzz_0[i] = ta_zz_zzz_0[i] * rpa_x - ta_zz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto ta_yyy_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 60) * nelems]);

    auto ta_yyy_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 61) * nelems]);

    auto ta_yyy_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 62) * nelems]);

    auto ta_yyy_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 63) * nelems]);

    auto ta_yyy_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 64) * nelems]);

    auto ta_yyy_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 65) * nelems]);

    auto ta_yyy_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 66) * nelems]);

    auto ta_yyy_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 67) * nelems]);

    auto ta_yyy_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 68) * nelems]);

    auto ta_yyy_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 69) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yyy_xxx_0[i] = 2.0 * ta_y_xxx_0[i] * fe_0 - 2.0 * ta_y_xxx_1[i] * fe_0 + ta_yy_xxx_0[i] * rpa_y - ta_yy_xxx_1[i] * pc_y[i];

        ta_yyy_xxy_0[i] = 2.0 * ta_y_xxy_0[i] * fe_0 - 2.0 * ta_y_xxy_1[i] * fe_0 + ta_yy_xx_0[i] * fe_0 - ta_yy_xx_1[i] * fe_0 + ta_yy_xxy_0[i] * rpa_y - ta_yy_xxy_1[i] * pc_y[i];

        ta_yyy_xxz_0[i] = 2.0 * ta_y_xxz_0[i] * fe_0 - 2.0 * ta_y_xxz_1[i] * fe_0 + ta_yy_xxz_0[i] * rpa_y - ta_yy_xxz_1[i] * pc_y[i];

        ta_yyy_xyy_0[i] = 2.0 * ta_y_xyy_0[i] * fe_0 - 2.0 * ta_y_xyy_1[i] * fe_0 + 2.0 * ta_yy_xy_0[i] * fe_0 - 2.0 * ta_yy_xy_1[i] * fe_0 + ta_yy_xyy_0[i] * rpa_y - ta_yy_xyy_1[i] * pc_y[i];

        ta_yyy_xyz_0[i] = 2.0 * ta_y_xyz_0[i] * fe_0 - 2.0 * ta_y_xyz_1[i] * fe_0 + ta_yy_xz_0[i] * fe_0 - ta_yy_xz_1[i] * fe_0 + ta_yy_xyz_0[i] * rpa_y - ta_yy_xyz_1[i] * pc_y[i];

        ta_yyy_xzz_0[i] = 2.0 * ta_y_xzz_0[i] * fe_0 - 2.0 * ta_y_xzz_1[i] * fe_0 + ta_yy_xzz_0[i] * rpa_y - ta_yy_xzz_1[i] * pc_y[i];

        ta_yyy_yyy_0[i] = 2.0 * ta_y_yyy_0[i] * fe_0 - 2.0 * ta_y_yyy_1[i] * fe_0 + 3.0 * ta_yy_yy_0[i] * fe_0 - 3.0 * ta_yy_yy_1[i] * fe_0 + ta_yy_yyy_0[i] * rpa_y - ta_yy_yyy_1[i] * pc_y[i];

        ta_yyy_yyz_0[i] = 2.0 * ta_y_yyz_0[i] * fe_0 - 2.0 * ta_y_yyz_1[i] * fe_0 + 2.0 * ta_yy_yz_0[i] * fe_0 - 2.0 * ta_yy_yz_1[i] * fe_0 + ta_yy_yyz_0[i] * rpa_y - ta_yy_yyz_1[i] * pc_y[i];

        ta_yyy_yzz_0[i] = 2.0 * ta_y_yzz_0[i] * fe_0 - 2.0 * ta_y_yzz_1[i] * fe_0 + ta_yy_zz_0[i] * fe_0 - ta_yy_zz_1[i] * fe_0 + ta_yy_yzz_0[i] * rpa_y - ta_yy_yzz_1[i] * pc_y[i];

        ta_yyy_zzz_0[i] = 2.0 * ta_y_zzz_0[i] * fe_0 - 2.0 * ta_y_zzz_1[i] * fe_0 + ta_yy_zzz_0[i] * rpa_y - ta_yy_zzz_1[i] * pc_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto ta_yyz_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 70) * nelems]);

    auto ta_yyz_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 71) * nelems]);

    auto ta_yyz_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 72) * nelems]);

    auto ta_yyz_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 73) * nelems]);

    auto ta_yyz_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 74) * nelems]);

    auto ta_yyz_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 75) * nelems]);

    auto ta_yyz_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 76) * nelems]);

    auto ta_yyz_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 77) * nelems]);

    auto ta_yyz_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 78) * nelems]);

    auto ta_yyz_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 79) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yyz_xxx_0[i] = ta_yy_xxx_0[i] * rpa_z - ta_yy_xxx_1[i] * pc_z[i];

        ta_yyz_xxy_0[i] = ta_yy_xxy_0[i] * rpa_z - ta_yy_xxy_1[i] * pc_z[i];

        ta_yyz_xxz_0[i] = ta_z_xxz_0[i] * fe_0 - ta_z_xxz_1[i] * fe_0 + ta_yz_xxz_0[i] * rpa_y - ta_yz_xxz_1[i] * pc_y[i];

        ta_yyz_xyy_0[i] = ta_yy_xyy_0[i] * rpa_z - ta_yy_xyy_1[i] * pc_z[i];

        ta_yyz_xyz_0[i] = ta_yy_xy_0[i] * fe_0 - ta_yy_xy_1[i] * fe_0 + ta_yy_xyz_0[i] * rpa_z - ta_yy_xyz_1[i] * pc_z[i];

        ta_yyz_xzz_0[i] = ta_z_xzz_0[i] * fe_0 - ta_z_xzz_1[i] * fe_0 + ta_yz_xzz_0[i] * rpa_y - ta_yz_xzz_1[i] * pc_y[i];

        ta_yyz_yyy_0[i] = ta_yy_yyy_0[i] * rpa_z - ta_yy_yyy_1[i] * pc_z[i];

        ta_yyz_yyz_0[i] = ta_yy_yy_0[i] * fe_0 - ta_yy_yy_1[i] * fe_0 + ta_yy_yyz_0[i] * rpa_z - ta_yy_yyz_1[i] * pc_z[i];

        ta_yyz_yzz_0[i] = 2.0 * ta_yy_yz_0[i] * fe_0 - 2.0 * ta_yy_yz_1[i] * fe_0 + ta_yy_yzz_0[i] * rpa_z - ta_yy_yzz_1[i] * pc_z[i];

        ta_yyz_zzz_0[i] = ta_z_zzz_0[i] * fe_0 - ta_z_zzz_1[i] * fe_0 + ta_yz_zzz_0[i] * rpa_y - ta_yz_zzz_1[i] * pc_y[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto ta_yzz_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 80) * nelems]);

    auto ta_yzz_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 81) * nelems]);

    auto ta_yzz_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 82) * nelems]);

    auto ta_yzz_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 83) * nelems]);

    auto ta_yzz_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 84) * nelems]);

    auto ta_yzz_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 85) * nelems]);

    auto ta_yzz_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 86) * nelems]);

    auto ta_yzz_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 87) * nelems]);

    auto ta_yzz_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 88) * nelems]);

    auto ta_yzz_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 89) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yzz_xxx_0[i] = ta_zz_xxx_0[i] * rpa_y - ta_zz_xxx_1[i] * pc_y[i];

        ta_yzz_xxy_0[i] = ta_zz_xx_0[i] * fe_0 - ta_zz_xx_1[i] * fe_0 + ta_zz_xxy_0[i] * rpa_y - ta_zz_xxy_1[i] * pc_y[i];

        ta_yzz_xxz_0[i] = ta_zz_xxz_0[i] * rpa_y - ta_zz_xxz_1[i] * pc_y[i];

        ta_yzz_xyy_0[i] = 2.0 * ta_zz_xy_0[i] * fe_0 - 2.0 * ta_zz_xy_1[i] * fe_0 + ta_zz_xyy_0[i] * rpa_y - ta_zz_xyy_1[i] * pc_y[i];

        ta_yzz_xyz_0[i] = ta_zz_xz_0[i] * fe_0 - ta_zz_xz_1[i] * fe_0 + ta_zz_xyz_0[i] * rpa_y - ta_zz_xyz_1[i] * pc_y[i];

        ta_yzz_xzz_0[i] = ta_zz_xzz_0[i] * rpa_y - ta_zz_xzz_1[i] * pc_y[i];

        ta_yzz_yyy_0[i] = 3.0 * ta_zz_yy_0[i] * fe_0 - 3.0 * ta_zz_yy_1[i] * fe_0 + ta_zz_yyy_0[i] * rpa_y - ta_zz_yyy_1[i] * pc_y[i];

        ta_yzz_yyz_0[i] = 2.0 * ta_zz_yz_0[i] * fe_0 - 2.0 * ta_zz_yz_1[i] * fe_0 + ta_zz_yyz_0[i] * rpa_y - ta_zz_yyz_1[i] * pc_y[i];

        ta_yzz_yzz_0[i] = ta_zz_zz_0[i] * fe_0 - ta_zz_zz_1[i] * fe_0 + ta_zz_yzz_0[i] * rpa_y - ta_zz_yzz_1[i] * pc_y[i];

        ta_yzz_zzz_0[i] = ta_zz_zzz_0[i] * rpa_y - ta_zz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto ta_zzz_xxx_0 = &(buffer.data()[(idx_npot_0_ff + 90) * nelems]);

    auto ta_zzz_xxy_0 = &(buffer.data()[(idx_npot_0_ff + 91) * nelems]);

    auto ta_zzz_xxz_0 = &(buffer.data()[(idx_npot_0_ff + 92) * nelems]);

    auto ta_zzz_xyy_0 = &(buffer.data()[(idx_npot_0_ff + 93) * nelems]);

    auto ta_zzz_xyz_0 = &(buffer.data()[(idx_npot_0_ff + 94) * nelems]);

    auto ta_zzz_xzz_0 = &(buffer.data()[(idx_npot_0_ff + 95) * nelems]);

    auto ta_zzz_yyy_0 = &(buffer.data()[(idx_npot_0_ff + 96) * nelems]);

    auto ta_zzz_yyz_0 = &(buffer.data()[(idx_npot_0_ff + 97) * nelems]);

    auto ta_zzz_yzz_0 = &(buffer.data()[(idx_npot_0_ff + 98) * nelems]);

    auto ta_zzz_zzz_0 = &(buffer.data()[(idx_npot_0_ff + 99) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_zzz_xxx_0[i] = 2.0 * ta_z_xxx_0[i] * fe_0 - 2.0 * ta_z_xxx_1[i] * fe_0 + ta_zz_xxx_0[i] * rpa_z - ta_zz_xxx_1[i] * pc_z[i];

        ta_zzz_xxy_0[i] = 2.0 * ta_z_xxy_0[i] * fe_0 - 2.0 * ta_z_xxy_1[i] * fe_0 + ta_zz_xxy_0[i] * rpa_z - ta_zz_xxy_1[i] * pc_z[i];

        ta_zzz_xxz_0[i] = 2.0 * ta_z_xxz_0[i] * fe_0 - 2.0 * ta_z_xxz_1[i] * fe_0 + ta_zz_xx_0[i] * fe_0 - ta_zz_xx_1[i] * fe_0 + ta_zz_xxz_0[i] * rpa_z - ta_zz_xxz_1[i] * pc_z[i];

        ta_zzz_xyy_0[i] = 2.0 * ta_z_xyy_0[i] * fe_0 - 2.0 * ta_z_xyy_1[i] * fe_0 + ta_zz_xyy_0[i] * rpa_z - ta_zz_xyy_1[i] * pc_z[i];

        ta_zzz_xyz_0[i] = 2.0 * ta_z_xyz_0[i] * fe_0 - 2.0 * ta_z_xyz_1[i] * fe_0 + ta_zz_xy_0[i] * fe_0 - ta_zz_xy_1[i] * fe_0 + ta_zz_xyz_0[i] * rpa_z - ta_zz_xyz_1[i] * pc_z[i];

        ta_zzz_xzz_0[i] = 2.0 * ta_z_xzz_0[i] * fe_0 - 2.0 * ta_z_xzz_1[i] * fe_0 + 2.0 * ta_zz_xz_0[i] * fe_0 - 2.0 * ta_zz_xz_1[i] * fe_0 + ta_zz_xzz_0[i] * rpa_z - ta_zz_xzz_1[i] * pc_z[i];

        ta_zzz_yyy_0[i] = 2.0 * ta_z_yyy_0[i] * fe_0 - 2.0 * ta_z_yyy_1[i] * fe_0 + ta_zz_yyy_0[i] * rpa_z - ta_zz_yyy_1[i] * pc_z[i];

        ta_zzz_yyz_0[i] = 2.0 * ta_z_yyz_0[i] * fe_0 - 2.0 * ta_z_yyz_1[i] * fe_0 + ta_zz_yy_0[i] * fe_0 - ta_zz_yy_1[i] * fe_0 + ta_zz_yyz_0[i] * rpa_z - ta_zz_yyz_1[i] * pc_z[i];

        ta_zzz_yzz_0[i] = 2.0 * ta_z_yzz_0[i] * fe_0 - 2.0 * ta_z_yzz_1[i] * fe_0 + 2.0 * ta_zz_yz_0[i] * fe_0 - 2.0 * ta_zz_yz_1[i] * fe_0 + ta_zz_yzz_0[i] * rpa_z - ta_zz_yzz_1[i] * pc_z[i];

        ta_zzz_zzz_0[i] = 2.0 * ta_z_zzz_0[i] * fe_0 - 2.0 * ta_z_zzz_1[i] * fe_0 + 3.0 * ta_zz_zz_0[i] * fe_0 - 3.0 * ta_zz_zz_1[i] * fe_0 + ta_zz_zzz_0[i] * rpa_z - ta_zz_zzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

