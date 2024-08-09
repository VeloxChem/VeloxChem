#include "NuclearPotentialGeom010PrimRecFS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_fs(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_fs,
                                        const size_t              idx_npot_geom_010_0_ps,
                                        const size_t              idx_npot_geom_010_1_ps,
                                        const size_t              idx_npot_1_ds,
                                        const size_t              idx_npot_geom_010_0_ds,
                                        const size_t              idx_npot_geom_010_1_ds,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : PS

    auto ta1_x_x_0_0 = pbuffer.data(idx_npot_geom_010_0_ps);

    auto ta1_x_y_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 1);

    auto ta1_x_z_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 2);

    auto ta1_y_x_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 3);

    auto ta1_y_y_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 4);

    auto ta1_y_z_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 5);

    auto ta1_z_x_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 6);

    auto ta1_z_y_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 7);

    auto ta1_z_z_0_0 = pbuffer.data(idx_npot_geom_010_0_ps + 8);

    // Set up components of auxiliary buffer : PS

    auto ta1_x_x_0_1 = pbuffer.data(idx_npot_geom_010_1_ps);

    auto ta1_x_y_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 1);

    auto ta1_x_z_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 2);

    auto ta1_y_x_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 3);

    auto ta1_y_y_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 4);

    auto ta1_y_z_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 5);

    auto ta1_z_x_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 6);

    auto ta1_z_y_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 7);

    auto ta1_z_z_0_1 = pbuffer.data(idx_npot_geom_010_1_ps + 8);

    // Set up components of auxiliary buffer : DS

    auto ta_xx_0_1 = pbuffer.data(idx_npot_1_ds);

    auto ta_yy_0_1 = pbuffer.data(idx_npot_1_ds + 3);

    auto ta_zz_0_1 = pbuffer.data(idx_npot_1_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds);

    auto ta1_x_xz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 2);

    auto ta1_x_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 3);

    auto ta1_x_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 5);

    auto ta1_y_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 6);

    auto ta1_y_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 9);

    auto ta1_y_yz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 10);

    auto ta1_y_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 11);

    auto ta1_z_xx_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 12);

    auto ta1_z_yy_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 15);

    auto ta1_z_yz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 16);

    auto ta1_z_zz_0_0 = pbuffer.data(idx_npot_geom_010_0_ds + 17);

    // Set up components of auxiliary buffer : DS

    auto ta1_x_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds);

    auto ta1_x_xz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 2);

    auto ta1_x_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 3);

    auto ta1_x_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 5);

    auto ta1_y_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 6);

    auto ta1_y_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 9);

    auto ta1_y_yz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 10);

    auto ta1_y_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 11);

    auto ta1_z_xx_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 12);

    auto ta1_z_yy_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 15);

    auto ta1_z_yz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 16);

    auto ta1_z_zz_0_1 = pbuffer.data(idx_npot_geom_010_1_ds + 17);

    // Set up components of targeted buffer : FS

    auto ta1_x_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs);

    auto ta1_x_xxy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 1);

    auto ta1_x_xxz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 2);

    auto ta1_x_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 3);

    auto ta1_x_xyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 4);

    auto ta1_x_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 5);

    auto ta1_x_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 6);

    auto ta1_x_yyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 7);

    auto ta1_x_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 8);

    auto ta1_x_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 9);

    auto ta1_y_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 10);

    auto ta1_y_xxy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 11);

    auto ta1_y_xxz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 12);

    auto ta1_y_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 13);

    auto ta1_y_xyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 14);

    auto ta1_y_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 15);

    auto ta1_y_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 16);

    auto ta1_y_yyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 17);

    auto ta1_y_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 18);

    auto ta1_y_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 19);

    auto ta1_z_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 20);

    auto ta1_z_xxy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 21);

    auto ta1_z_xxz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 22);

    auto ta1_z_xyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 23);

    auto ta1_z_xyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 24);

    auto ta1_z_xzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 25);

    auto ta1_z_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 26);

    auto ta1_z_yyz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 27);

    auto ta1_z_yzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 28);

    auto ta1_z_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 29);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta1_x_x_0_0,   \
                             ta1_x_x_0_1,   \
                             ta1_x_xx_0_0,  \
                             ta1_x_xx_0_1,  \
                             ta1_x_xxx_0_0, \
                             ta1_x_xxy_0_0, \
                             ta1_x_xxz_0_0, \
                             ta1_x_xyy_0_0, \
                             ta1_x_xyz_0_0, \
                             ta1_x_xz_0_0,  \
                             ta1_x_xz_0_1,  \
                             ta1_x_xzz_0_0, \
                             ta1_x_y_0_0,   \
                             ta1_x_y_0_1,   \
                             ta1_x_yy_0_0,  \
                             ta1_x_yy_0_1,  \
                             ta1_x_yyy_0_0, \
                             ta1_x_yyz_0_0, \
                             ta1_x_yzz_0_0, \
                             ta1_x_z_0_0,   \
                             ta1_x_z_0_1,   \
                             ta1_x_zz_0_0,  \
                             ta1_x_zz_0_1,  \
                             ta1_x_zzz_0_0, \
                             ta1_y_x_0_0,   \
                             ta1_y_x_0_1,   \
                             ta1_y_xx_0_0,  \
                             ta1_y_xx_0_1,  \
                             ta1_y_xxx_0_0, \
                             ta1_y_xxy_0_0, \
                             ta1_y_xxz_0_0, \
                             ta1_y_xyy_0_0, \
                             ta1_y_xyz_0_0, \
                             ta1_y_xzz_0_0, \
                             ta1_y_y_0_0,   \
                             ta1_y_y_0_1,   \
                             ta1_y_yy_0_0,  \
                             ta1_y_yy_0_1,  \
                             ta1_y_yyy_0_0, \
                             ta1_y_yyz_0_0, \
                             ta1_y_yz_0_0,  \
                             ta1_y_yz_0_1,  \
                             ta1_y_yzz_0_0, \
                             ta1_y_z_0_0,   \
                             ta1_y_z_0_1,   \
                             ta1_y_zz_0_0,  \
                             ta1_y_zz_0_1,  \
                             ta1_y_zzz_0_0, \
                             ta1_z_x_0_0,   \
                             ta1_z_x_0_1,   \
                             ta1_z_xx_0_0,  \
                             ta1_z_xx_0_1,  \
                             ta1_z_xxx_0_0, \
                             ta1_z_xxy_0_0, \
                             ta1_z_xxz_0_0, \
                             ta1_z_xyy_0_0, \
                             ta1_z_xyz_0_0, \
                             ta1_z_xzz_0_0, \
                             ta1_z_y_0_0,   \
                             ta1_z_y_0_1,   \
                             ta1_z_yy_0_0,  \
                             ta1_z_yy_0_1,  \
                             ta1_z_yyy_0_0, \
                             ta1_z_yyz_0_0, \
                             ta1_z_yz_0_0,  \
                             ta1_z_yz_0_1,  \
                             ta1_z_yzz_0_0, \
                             ta1_z_z_0_0,   \
                             ta1_z_z_0_1,   \
                             ta1_z_zz_0_0,  \
                             ta1_z_zz_0_1,  \
                             ta1_z_zzz_0_0, \
                             ta_xx_0_1,     \
                             ta_yy_0_1,     \
                             ta_zz_0_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_0_0[i] = 2.0 * ta1_x_x_0_0[i] * fe_0 - 2.0 * ta1_x_x_0_1[i] * fe_0 + ta_xx_0_1[i] +
                           ta1_x_xx_0_0[i] * pa_x[i] - ta1_x_xx_0_1[i] * pc_x[i];

        ta1_x_xxy_0_0[i] = ta1_x_xx_0_0[i] * pa_y[i] - ta1_x_xx_0_1[i] * pc_y[i];

        ta1_x_xxz_0_0[i] = ta1_x_xx_0_0[i] * pa_z[i] - ta1_x_xx_0_1[i] * pc_z[i];

        ta1_x_xyy_0_0[i] = ta_yy_0_1[i] + ta1_x_yy_0_0[i] * pa_x[i] - ta1_x_yy_0_1[i] * pc_x[i];

        ta1_x_xyz_0_0[i] = ta1_x_xz_0_0[i] * pa_y[i] - ta1_x_xz_0_1[i] * pc_y[i];

        ta1_x_xzz_0_0[i] = ta_zz_0_1[i] + ta1_x_zz_0_0[i] * pa_x[i] - ta1_x_zz_0_1[i] * pc_x[i];

        ta1_x_yyy_0_0[i] = 2.0 * ta1_x_y_0_0[i] * fe_0 - 2.0 * ta1_x_y_0_1[i] * fe_0 + ta1_x_yy_0_0[i] * pa_y[i] -
                           ta1_x_yy_0_1[i] * pc_y[i];

        ta1_x_yyz_0_0[i] = ta1_x_yy_0_0[i] * pa_z[i] - ta1_x_yy_0_1[i] * pc_z[i];

        ta1_x_yzz_0_0[i] = ta1_x_zz_0_0[i] * pa_y[i] - ta1_x_zz_0_1[i] * pc_y[i];

        ta1_x_zzz_0_0[i] = 2.0 * ta1_x_z_0_0[i] * fe_0 - 2.0 * ta1_x_z_0_1[i] * fe_0 + ta1_x_zz_0_0[i] * pa_z[i] -
                           ta1_x_zz_0_1[i] * pc_z[i];

        ta1_y_xxx_0_0[i] = 2.0 * ta1_y_x_0_0[i] * fe_0 - 2.0 * ta1_y_x_0_1[i] * fe_0 + ta1_y_xx_0_0[i] * pa_x[i] -
                           ta1_y_xx_0_1[i] * pc_x[i];

        ta1_y_xxy_0_0[i] = ta_xx_0_1[i] + ta1_y_xx_0_0[i] * pa_y[i] - ta1_y_xx_0_1[i] * pc_y[i];

        ta1_y_xxz_0_0[i] = ta1_y_xx_0_0[i] * pa_z[i] - ta1_y_xx_0_1[i] * pc_z[i];

        ta1_y_xyy_0_0[i] = ta1_y_yy_0_0[i] * pa_x[i] - ta1_y_yy_0_1[i] * pc_x[i];

        ta1_y_xyz_0_0[i] = ta1_y_yz_0_0[i] * pa_x[i] - ta1_y_yz_0_1[i] * pc_x[i];

        ta1_y_xzz_0_0[i] = ta1_y_zz_0_0[i] * pa_x[i] - ta1_y_zz_0_1[i] * pc_x[i];

        ta1_y_yyy_0_0[i] = 2.0 * ta1_y_y_0_0[i] * fe_0 - 2.0 * ta1_y_y_0_1[i] * fe_0 + ta_yy_0_1[i] +
                           ta1_y_yy_0_0[i] * pa_y[i] - ta1_y_yy_0_1[i] * pc_y[i];

        ta1_y_yyz_0_0[i] = ta1_y_yy_0_0[i] * pa_z[i] - ta1_y_yy_0_1[i] * pc_z[i];

        ta1_y_yzz_0_0[i] = ta_zz_0_1[i] + ta1_y_zz_0_0[i] * pa_y[i] - ta1_y_zz_0_1[i] * pc_y[i];

        ta1_y_zzz_0_0[i] = 2.0 * ta1_y_z_0_0[i] * fe_0 - 2.0 * ta1_y_z_0_1[i] * fe_0 + ta1_y_zz_0_0[i] * pa_z[i] -
                           ta1_y_zz_0_1[i] * pc_z[i];

        ta1_z_xxx_0_0[i] = 2.0 * ta1_z_x_0_0[i] * fe_0 - 2.0 * ta1_z_x_0_1[i] * fe_0 + ta1_z_xx_0_0[i] * pa_x[i] -
                           ta1_z_xx_0_1[i] * pc_x[i];

        ta1_z_xxy_0_0[i] = ta1_z_xx_0_0[i] * pa_y[i] - ta1_z_xx_0_1[i] * pc_y[i];

        ta1_z_xxz_0_0[i] = ta_xx_0_1[i] + ta1_z_xx_0_0[i] * pa_z[i] - ta1_z_xx_0_1[i] * pc_z[i];

        ta1_z_xyy_0_0[i] = ta1_z_yy_0_0[i] * pa_x[i] - ta1_z_yy_0_1[i] * pc_x[i];

        ta1_z_xyz_0_0[i] = ta1_z_yz_0_0[i] * pa_x[i] - ta1_z_yz_0_1[i] * pc_x[i];

        ta1_z_xzz_0_0[i] = ta1_z_zz_0_0[i] * pa_x[i] - ta1_z_zz_0_1[i] * pc_x[i];

        ta1_z_yyy_0_0[i] = 2.0 * ta1_z_y_0_0[i] * fe_0 - 2.0 * ta1_z_y_0_1[i] * fe_0 + ta1_z_yy_0_0[i] * pa_y[i] -
                           ta1_z_yy_0_1[i] * pc_y[i];

        ta1_z_yyz_0_0[i] = ta_yy_0_1[i] + ta1_z_yy_0_0[i] * pa_z[i] - ta1_z_yy_0_1[i] * pc_z[i];

        ta1_z_yzz_0_0[i] = ta1_z_zz_0_0[i] * pa_y[i] - ta1_z_zz_0_1[i] * pc_y[i];

        ta1_z_zzz_0_0[i] = 2.0 * ta1_z_z_0_0[i] * fe_0 - 2.0 * ta1_z_z_0_1[i] * fe_0 + ta_zz_0_1[i] +
                           ta1_z_zz_0_0[i] * pa_z[i] - ta1_z_zz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
