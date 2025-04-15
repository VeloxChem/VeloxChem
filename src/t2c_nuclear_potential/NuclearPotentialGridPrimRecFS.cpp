#include "NuclearPotentialGridPrimRecFS.hpp"

namespace npotrec { // npotrec namespace

auto
comp_on_grid_prim_nuclear_potential_fs(CSubMatrix&  buffer,
                                       const size_t idx_npot_0_fs,
                                       const size_t idx_npot_0_ps,
                                       const size_t idx_npot_1_ps,
                                       const size_t idx_npot_0_ds,
                                       const size_t idx_npot_1_ds,
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

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_0 = &(buffer.data()[idx_npot_0_ps * nelems]);

    auto ta_y_0_0 = &(buffer.data()[(idx_npot_0_ps + 1) * nelems]);

    auto ta_z_0_0 = &(buffer.data()[(idx_npot_0_ps + 2) * nelems]);

    // Set up components of auxiliary buffer : PS

    auto ta_x_0_1 = &(buffer.data()[idx_npot_1_ps * nelems]);

    auto ta_y_0_1 = &(buffer.data()[(idx_npot_1_ps + 1) * nelems]);

    auto ta_z_0_1 = &(buffer.data()[(idx_npot_1_ps + 2) * nelems]);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_0 = &(buffer.data()[idx_npot_0_ds * nelems]);

    auto ta_yy_0_0 = &(buffer.data()[(idx_npot_0_ds + 3) * nelems]);

    auto ta_yz_0_0 = &(buffer.data()[(idx_npot_0_ds + 4) * nelems]);

    auto ta_zz_0_0 = &(buffer.data()[(idx_npot_0_ds + 5) * nelems]);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_1 = &(buffer.data()[idx_npot_1_ds * nelems]);

    auto ta_yy_0_1 = &(buffer.data()[(idx_npot_1_ds + 3) * nelems]);

    auto ta_yz_0_1 = &(buffer.data()[(idx_npot_1_ds + 4) * nelems]);

    auto ta_zz_0_1 = &(buffer.data()[(idx_npot_1_ds + 5) * nelems]);

    // Set up components of targeted buffer : FS

    auto ta_xxx_0_0 = &(buffer.data()[idx_npot_0_fs * nelems]);

    auto ta_xxy_0_0 = &(buffer.data()[(idx_npot_0_fs + 1) * nelems]);

    auto ta_xxz_0_0 = &(buffer.data()[(idx_npot_0_fs + 2) * nelems]);

    auto ta_xyy_0_0 = &(buffer.data()[(idx_npot_0_fs + 3) * nelems]);

    auto ta_xyz_0_0 = &(buffer.data()[(idx_npot_0_fs + 4) * nelems]);

    auto ta_xzz_0_0 = &(buffer.data()[(idx_npot_0_fs + 5) * nelems]);

    auto ta_yyy_0_0 = &(buffer.data()[(idx_npot_0_fs + 6) * nelems]);

    auto ta_yyz_0_0 = &(buffer.data()[(idx_npot_0_fs + 7) * nelems]);

    auto ta_yzz_0_0 = &(buffer.data()[(idx_npot_0_fs + 8) * nelems]);

    auto ta_zzz_0_0 = &(buffer.data()[(idx_npot_0_fs + 9) * nelems]);

    #pragma omp simd
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / factor;

        ta_xxx_0_0[i] = 2.0 * ta_x_0_0[i] * fe_0 - 2.0 * ta_x_0_1[i] * fe_0 + ta_xx_0_0[i] * rpa_x - ta_xx_0_1[i] * pc_x[i];

        ta_xxy_0_0[i] = ta_xx_0_0[i] * rpa_y - ta_xx_0_1[i] * pc_y[i];

        ta_xxz_0_0[i] = ta_xx_0_0[i] * rpa_z - ta_xx_0_1[i] * pc_z[i];

        ta_xyy_0_0[i] = ta_yy_0_0[i] * rpa_x - ta_yy_0_1[i] * pc_x[i];

        ta_xyz_0_0[i] = ta_yz_0_0[i] * rpa_x - ta_yz_0_1[i] * pc_x[i];

        ta_xzz_0_0[i] = ta_zz_0_0[i] * rpa_x - ta_zz_0_1[i] * pc_x[i];

        ta_yyy_0_0[i] = 2.0 * ta_y_0_0[i] * fe_0 - 2.0 * ta_y_0_1[i] * fe_0 + ta_yy_0_0[i] * rpa_y - ta_yy_0_1[i] * pc_y[i];

        ta_yyz_0_0[i] = ta_yy_0_0[i] * rpa_z - ta_yy_0_1[i] * pc_z[i];

        ta_yzz_0_0[i] = ta_zz_0_0[i] * rpa_y - ta_zz_0_1[i] * pc_y[i];

        ta_zzz_0_0[i] = 2.0 * ta_z_0_0[i] * fe_0 - 2.0 * ta_z_0_1[i] * fe_0 + ta_zz_0_0[i] * rpa_z - ta_zz_0_1[i] * pc_z[i];
    }
}

} // npotrec namespace

