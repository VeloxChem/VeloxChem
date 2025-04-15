#include "NuclearPotentialGridPrimRecFP.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_fp(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_fp,
                                       const size_t idx_npot_0_pp,
                                       const size_t idx_npot_1_pp,
                                       const size_t idx_npot_0_ds,
                                       const size_t idx_npot_1_ds,
                                       const size_t idx_npot_0_dp,
                                       const size_t idx_npot_1_dp,
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

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_0 = &(buffer.data()[idx_npot_0_ds * nelems]);

    auto ta_yy_0_0 = &(buffer.data()[(idx_npot_0_ds + 3) * nelems]);

    auto ta_zz_0_0 = &(buffer.data()[(idx_npot_0_ds + 5) * nelems]);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_1 = &(buffer.data()[idx_npot_1_ds * nelems]);

    auto ta_yy_0_1 = &(buffer.data()[(idx_npot_1_ds + 3) * nelems]);

    auto ta_zz_0_1 = &(buffer.data()[(idx_npot_1_ds + 5) * nelems]);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_0 = &(buffer.data()[idx_npot_0_dp * nelems]);

    auto ta_xx_y_0 = &(buffer.data()[(idx_npot_0_dp + 1) * nelems]);

    auto ta_xx_z_0 = &(buffer.data()[(idx_npot_0_dp + 2) * nelems]);

    auto ta_xy_y_0 = &(buffer.data()[(idx_npot_0_dp + 4) * nelems]);

    auto ta_xz_x_0 = &(buffer.data()[(idx_npot_0_dp + 6) * nelems]);

    auto ta_xz_z_0 = &(buffer.data()[(idx_npot_0_dp + 8) * nelems]);

    auto ta_yy_x_0 = &(buffer.data()[(idx_npot_0_dp + 9) * nelems]);

    auto ta_yy_y_0 = &(buffer.data()[(idx_npot_0_dp + 10) * nelems]);

    auto ta_yy_z_0 = &(buffer.data()[(idx_npot_0_dp + 11) * nelems]);

    auto ta_yz_y_0 = &(buffer.data()[(idx_npot_0_dp + 13) * nelems]);

    auto ta_yz_z_0 = &(buffer.data()[(idx_npot_0_dp + 14) * nelems]);

    auto ta_zz_x_0 = &(buffer.data()[(idx_npot_0_dp + 15) * nelems]);

    auto ta_zz_y_0 = &(buffer.data()[(idx_npot_0_dp + 16) * nelems]);

    auto ta_zz_z_0 = &(buffer.data()[(idx_npot_0_dp + 17) * nelems]);

    // Set up components of auxiliary buffer : DP

    auto ta_xx_x_1 = &(buffer.data()[idx_npot_1_dp * nelems]);

    auto ta_xx_y_1 = &(buffer.data()[(idx_npot_1_dp + 1) * nelems]);

    auto ta_xx_z_1 = &(buffer.data()[(idx_npot_1_dp + 2) * nelems]);

    auto ta_xy_y_1 = &(buffer.data()[(idx_npot_1_dp + 4) * nelems]);

    auto ta_xz_x_1 = &(buffer.data()[(idx_npot_1_dp + 6) * nelems]);

    auto ta_xz_z_1 = &(buffer.data()[(idx_npot_1_dp + 8) * nelems]);

    auto ta_yy_x_1 = &(buffer.data()[(idx_npot_1_dp + 9) * nelems]);

    auto ta_yy_y_1 = &(buffer.data()[(idx_npot_1_dp + 10) * nelems]);

    auto ta_yy_z_1 = &(buffer.data()[(idx_npot_1_dp + 11) * nelems]);

    auto ta_yz_y_1 = &(buffer.data()[(idx_npot_1_dp + 13) * nelems]);

    auto ta_yz_z_1 = &(buffer.data()[(idx_npot_1_dp + 14) * nelems]);

    auto ta_zz_x_1 = &(buffer.data()[(idx_npot_1_dp + 15) * nelems]);

    auto ta_zz_y_1 = &(buffer.data()[(idx_npot_1_dp + 16) * nelems]);

    auto ta_zz_z_1 = &(buffer.data()[(idx_npot_1_dp + 17) * nelems]);

    // Set up 0-3 components of targeted buffer : FP

    auto ta_xxx_x_0 = &(buffer.data()[idx_npot_0_fp * nelems]);

    auto ta_xxx_y_0 = &(buffer.data()[(idx_npot_0_fp + 1) * nelems]);

    auto ta_xxx_z_0 = &(buffer.data()[(idx_npot_0_fp + 2) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxx_x_0[i] = 2.0 * ta_x_x_0[i] * fe_0 - 2.0 * ta_x_x_1[i] * fe_0 + ta_xx_0_0[i] * fe_0 - ta_xx_0_1[i] * fe_0 + ta_xx_x_0[i] * rpa_x - ta_xx_x_1[i] * pc_x[i];

        ta_xxx_y_0[i] = 2.0 * ta_x_y_0[i] * fe_0 - 2.0 * ta_x_y_1[i] * fe_0 + ta_xx_y_0[i] * rpa_x - ta_xx_y_1[i] * pc_x[i];

        ta_xxx_z_0[i] = 2.0 * ta_x_z_0[i] * fe_0 - 2.0 * ta_x_z_1[i] * fe_0 + ta_xx_z_0[i] * rpa_x - ta_xx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto ta_xxy_x_0 = &(buffer.data()[(idx_npot_0_fp + 3) * nelems]);

    auto ta_xxy_y_0 = &(buffer.data()[(idx_npot_0_fp + 4) * nelems]);

    auto ta_xxy_z_0 = &(buffer.data()[(idx_npot_0_fp + 5) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxy_x_0[i] = ta_xx_x_0[i] * rpa_y - ta_xx_x_1[i] * pc_y[i];

        ta_xxy_y_0[i] = ta_y_y_0[i] * fe_0 - ta_y_y_1[i] * fe_0 + ta_xy_y_0[i] * rpa_x - ta_xy_y_1[i] * pc_x[i];

        ta_xxy_z_0[i] = ta_xx_z_0[i] * rpa_y - ta_xx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto ta_xxz_x_0 = &(buffer.data()[(idx_npot_0_fp + 6) * nelems]);

    auto ta_xxz_y_0 = &(buffer.data()[(idx_npot_0_fp + 7) * nelems]);

    auto ta_xxz_z_0 = &(buffer.data()[(idx_npot_0_fp + 8) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxz_x_0[i] = ta_xx_x_0[i] * rpa_z - ta_xx_x_1[i] * pc_z[i];

        ta_xxz_y_0[i] = ta_xx_y_0[i] * rpa_z - ta_xx_y_1[i] * pc_z[i];

        ta_xxz_z_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_xz_z_0[i] * rpa_x - ta_xz_z_1[i] * pc_x[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto ta_xyy_x_0 = &(buffer.data()[(idx_npot_0_fp + 9) * nelems]);

    auto ta_xyy_y_0 = &(buffer.data()[(idx_npot_0_fp + 10) * nelems]);

    auto ta_xyy_z_0 = &(buffer.data()[(idx_npot_0_fp + 11) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xyy_x_0[i] = ta_yy_0_0[i] * fe_0 - ta_yy_0_1[i] * fe_0 + ta_yy_x_0[i] * rpa_x - ta_yy_x_1[i] * pc_x[i];

        ta_xyy_y_0[i] = ta_yy_y_0[i] * rpa_x - ta_yy_y_1[i] * pc_x[i];

        ta_xyy_z_0[i] = ta_yy_z_0[i] * rpa_x - ta_yy_z_1[i] * pc_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto ta_xyz_x_0 = &(buffer.data()[(idx_npot_0_fp + 12) * nelems]);

    auto ta_xyz_y_0 = &(buffer.data()[(idx_npot_0_fp + 13) * nelems]);

    auto ta_xyz_z_0 = &(buffer.data()[(idx_npot_0_fp + 14) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        ta_xyz_x_0[i] = ta_xz_x_0[i] * rpa_y - ta_xz_x_1[i] * pc_y[i];

        ta_xyz_y_0[i] = ta_yz_y_0[i] * rpa_x - ta_yz_y_1[i] * pc_x[i];

        ta_xyz_z_0[i] = ta_yz_z_0[i] * rpa_x - ta_yz_z_1[i] * pc_x[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto ta_xzz_x_0 = &(buffer.data()[(idx_npot_0_fp + 15) * nelems]);

    auto ta_xzz_y_0 = &(buffer.data()[(idx_npot_0_fp + 16) * nelems]);

    auto ta_xzz_z_0 = &(buffer.data()[(idx_npot_0_fp + 17) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xzz_x_0[i] = ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_zz_x_0[i] * rpa_x - ta_zz_x_1[i] * pc_x[i];

        ta_xzz_y_0[i] = ta_zz_y_0[i] * rpa_x - ta_zz_y_1[i] * pc_x[i];

        ta_xzz_z_0[i] = ta_zz_z_0[i] * rpa_x - ta_zz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto ta_yyy_x_0 = &(buffer.data()[(idx_npot_0_fp + 18) * nelems]);

    auto ta_yyy_y_0 = &(buffer.data()[(idx_npot_0_fp + 19) * nelems]);

    auto ta_yyy_z_0 = &(buffer.data()[(idx_npot_0_fp + 20) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yyy_x_0[i] = 2.0 * ta_y_x_0[i] * fe_0 - 2.0 * ta_y_x_1[i] * fe_0 + ta_yy_x_0[i] * rpa_y - ta_yy_x_1[i] * pc_y[i];

        ta_yyy_y_0[i] = 2.0 * ta_y_y_0[i] * fe_0 - 2.0 * ta_y_y_1[i] * fe_0 + ta_yy_0_0[i] * fe_0 - ta_yy_0_1[i] * fe_0 + ta_yy_y_0[i] * rpa_y - ta_yy_y_1[i] * pc_y[i];

        ta_yyy_z_0[i] = 2.0 * ta_y_z_0[i] * fe_0 - 2.0 * ta_y_z_1[i] * fe_0 + ta_yy_z_0[i] * rpa_y - ta_yy_z_1[i] * pc_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto ta_yyz_x_0 = &(buffer.data()[(idx_npot_0_fp + 21) * nelems]);

    auto ta_yyz_y_0 = &(buffer.data()[(idx_npot_0_fp + 22) * nelems]);

    auto ta_yyz_z_0 = &(buffer.data()[(idx_npot_0_fp + 23) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yyz_x_0[i] = ta_yy_x_0[i] * rpa_z - ta_yy_x_1[i] * pc_z[i];

        ta_yyz_y_0[i] = ta_yy_y_0[i] * rpa_z - ta_yy_y_1[i] * pc_z[i];

        ta_yyz_z_0[i] = ta_z_z_0[i] * fe_0 - ta_z_z_1[i] * fe_0 + ta_yz_z_0[i] * rpa_y - ta_yz_z_1[i] * pc_y[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto ta_yzz_x_0 = &(buffer.data()[(idx_npot_0_fp + 24) * nelems]);

    auto ta_yzz_y_0 = &(buffer.data()[(idx_npot_0_fp + 25) * nelems]);

    auto ta_yzz_z_0 = &(buffer.data()[(idx_npot_0_fp + 26) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_yzz_x_0[i] = ta_zz_x_0[i] * rpa_y - ta_zz_x_1[i] * pc_y[i];

        ta_yzz_y_0[i] = ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_zz_y_0[i] * rpa_y - ta_zz_y_1[i] * pc_y[i];

        ta_yzz_z_0[i] = ta_zz_z_0[i] * rpa_y - ta_zz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto ta_zzz_x_0 = &(buffer.data()[(idx_npot_0_fp + 27) * nelems]);

    auto ta_zzz_y_0 = &(buffer.data()[(idx_npot_0_fp + 28) * nelems]);

    auto ta_zzz_z_0 = &(buffer.data()[(idx_npot_0_fp + 29) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_zzz_x_0[i] = 2.0 * ta_z_x_0[i] * fe_0 - 2.0 * ta_z_x_1[i] * fe_0 + ta_zz_x_0[i] * rpa_z - ta_zz_x_1[i] * pc_z[i];

        ta_zzz_y_0[i] = 2.0 * ta_z_y_0[i] * fe_0 - 2.0 * ta_z_y_1[i] * fe_0 + ta_zz_y_0[i] * rpa_z - ta_zz_y_1[i] * pc_z[i];

        ta_zzz_z_0[i] = 2.0 * ta_z_z_0[i] * fe_0 - 2.0 * ta_z_z_1[i] * fe_0 + ta_zz_0_0[i] * fe_0 - ta_zz_0_1[i] * fe_0 + ta_zz_z_0[i] * rpa_z - ta_zz_z_1[i] * pc_z[i];
    }

}

} // npotrec namespace

