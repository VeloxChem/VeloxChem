#include "NuclearPotentialGridPrimRecSF.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_sf(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_sf,
                                       const size_t idx_npot_0_sp,
                                       const size_t idx_npot_1_sp,
                                       const size_t idx_npot_0_sd,
                                       const size_t idx_npot_1_sd,
                                       const double rpb_x,
                                       const double rpb_y,
                                       const double rpb_z,
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

    auto ta_0_yy_0 = &(buffer.data()[(idx_npot_0_sd + 3) * nelems]);

    auto ta_0_yz_0 = &(buffer.data()[(idx_npot_0_sd + 4) * nelems]);

    auto ta_0_zz_0 = &(buffer.data()[(idx_npot_0_sd + 5) * nelems]);

    // Set up components of auxiliary buffer : SD

    auto ta_0_xx_1 = &(buffer.data()[idx_npot_1_sd * nelems]);

    auto ta_0_yy_1 = &(buffer.data()[(idx_npot_1_sd + 3) * nelems]);

    auto ta_0_yz_1 = &(buffer.data()[(idx_npot_1_sd + 4) * nelems]);

    auto ta_0_zz_1 = &(buffer.data()[(idx_npot_1_sd + 5) * nelems]);

    // Set up components of targeted buffer : SF

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

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_0_xxx_0[i] = 2.0 * ta_0_x_0[i] * fe_0 - 2.0 * ta_0_x_1[i] * fe_0 + ta_0_xx_0[i] * rpb_x - ta_0_xx_1[i] * pc_x[i];

        ta_0_xxy_0[i] = ta_0_xx_0[i] * rpb_y - ta_0_xx_1[i] * pc_y[i];

        ta_0_xxz_0[i] = ta_0_xx_0[i] * rpb_z - ta_0_xx_1[i] * pc_z[i];

        ta_0_xyy_0[i] = ta_0_yy_0[i] * rpb_x - ta_0_yy_1[i] * pc_x[i];

        ta_0_xyz_0[i] = ta_0_yz_0[i] * rpb_x - ta_0_yz_1[i] * pc_x[i];

        ta_0_xzz_0[i] = ta_0_zz_0[i] * rpb_x - ta_0_zz_1[i] * pc_x[i];

        ta_0_yyy_0[i] = 2.0 * ta_0_y_0[i] * fe_0 - 2.0 * ta_0_y_1[i] * fe_0 + ta_0_yy_0[i] * rpb_y - ta_0_yy_1[i] * pc_y[i];

        ta_0_yyz_0[i] = ta_0_yy_0[i] * rpb_z - ta_0_yy_1[i] * pc_z[i];

        ta_0_yzz_0[i] = ta_0_zz_0[i] * rpb_y - ta_0_zz_1[i] * pc_y[i];

        ta_0_zzz_0[i] = 2.0 * ta_0_z_0[i] * fe_0 - 2.0 * ta_0_z_1[i] * fe_0 + ta_0_zz_0[i] * rpb_z - ta_0_zz_1[i] * pc_z[i];
    }
}

} // npotrec namespace

