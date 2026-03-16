#include "GeometricalDerivatives010ForDS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_ds(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_ds,
                         const int idx_op_ps,
                         const int idx_op_dp,
                         const int idx_op_fs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : PS

    auto tr_x_0 = pbuffer.data(idx_op_ps);

    auto tr_y_0 = pbuffer.data(idx_op_ps + 1);

    auto tr_z_0 = pbuffer.data(idx_op_ps + 2);

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto tr_xxx_0 = pbuffer.data(idx_op_fs);

    auto tr_xxy_0 = pbuffer.data(idx_op_fs + 1);

    auto tr_xxz_0 = pbuffer.data(idx_op_fs + 2);

    auto tr_xyy_0 = pbuffer.data(idx_op_fs + 3);

    auto tr_xyz_0 = pbuffer.data(idx_op_fs + 4);

    auto tr_xzz_0 = pbuffer.data(idx_op_fs + 5);

    auto tr_yyy_0 = pbuffer.data(idx_op_fs + 6);

    auto tr_yyz_0 = pbuffer.data(idx_op_fs + 7);

    auto tr_yzz_0 = pbuffer.data(idx_op_fs + 8);

    auto tr_zzz_0 = pbuffer.data(idx_op_fs + 9);

    // Set up components of targeted buffer : DS

    auto tr_0_0_x_xx_0 = pbuffer.data(idx_op_geom_010_ds);

    auto tr_0_0_x_xy_0 = pbuffer.data(idx_op_geom_010_ds + 1);

    auto tr_0_0_x_xz_0 = pbuffer.data(idx_op_geom_010_ds + 2);

    auto tr_0_0_x_yy_0 = pbuffer.data(idx_op_geom_010_ds + 3);

    auto tr_0_0_x_yz_0 = pbuffer.data(idx_op_geom_010_ds + 4);

    auto tr_0_0_x_zz_0 = pbuffer.data(idx_op_geom_010_ds + 5);

    auto tr_0_0_y_xx_0 = pbuffer.data(idx_op_geom_010_ds + 6);

    auto tr_0_0_y_xy_0 = pbuffer.data(idx_op_geom_010_ds + 7);

    auto tr_0_0_y_xz_0 = pbuffer.data(idx_op_geom_010_ds + 8);

    auto tr_0_0_y_yy_0 = pbuffer.data(idx_op_geom_010_ds + 9);

    auto tr_0_0_y_yz_0 = pbuffer.data(idx_op_geom_010_ds + 10);

    auto tr_0_0_y_zz_0 = pbuffer.data(idx_op_geom_010_ds + 11);

    auto tr_0_0_z_xx_0 = pbuffer.data(idx_op_geom_010_ds + 12);

    auto tr_0_0_z_xy_0 = pbuffer.data(idx_op_geom_010_ds + 13);

    auto tr_0_0_z_xz_0 = pbuffer.data(idx_op_geom_010_ds + 14);

    auto tr_0_0_z_yy_0 = pbuffer.data(idx_op_geom_010_ds + 15);

    auto tr_0_0_z_yz_0 = pbuffer.data(idx_op_geom_010_ds + 16);

    auto tr_0_0_z_zz_0 = pbuffer.data(idx_op_geom_010_ds + 17);

    #pragma omp simd aligned(tr_0_0_x_xx_0, tr_0_0_x_xy_0, tr_0_0_x_xz_0, tr_0_0_x_yy_0, tr_0_0_x_yz_0, tr_0_0_x_zz_0, tr_0_0_y_xx_0, tr_0_0_y_xy_0, tr_0_0_y_xz_0, tr_0_0_y_yy_0, tr_0_0_y_yz_0, tr_0_0_y_zz_0, tr_0_0_z_xx_0, tr_0_0_z_xy_0, tr_0_0_z_xz_0, tr_0_0_z_yy_0, tr_0_0_z_yz_0, tr_0_0_z_zz_0, tr_x_0, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxy_0, tr_xxz_0, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyz_0, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_y_0, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyz_0, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_z_0, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xx_0[i] = 2.0 * tr_xxx_0[i] * tbe_0 + 2.0 * tr_xx_x[i] * tke_0 - 2.0 * tr_x_0[i];

        tr_0_0_x_xy_0[i] = 2.0 * tr_xxy_0[i] * tbe_0 + 2.0 * tr_xy_x[i] * tke_0 - tr_y_0[i];

        tr_0_0_x_xz_0[i] = 2.0 * tr_xxz_0[i] * tbe_0 + 2.0 * tr_xz_x[i] * tke_0 - tr_z_0[i];

        tr_0_0_x_yy_0[i] = 2.0 * tr_xyy_0[i] * tbe_0 + 2.0 * tr_yy_x[i] * tke_0;

        tr_0_0_x_yz_0[i] = 2.0 * tr_xyz_0[i] * tbe_0 + 2.0 * tr_yz_x[i] * tke_0;

        tr_0_0_x_zz_0[i] = 2.0 * tr_xzz_0[i] * tbe_0 + 2.0 * tr_zz_x[i] * tke_0;

        tr_0_0_y_xx_0[i] = 2.0 * tr_xxy_0[i] * tbe_0 + 2.0 * tr_xx_y[i] * tke_0;

        tr_0_0_y_xy_0[i] = 2.0 * tr_xyy_0[i] * tbe_0 + 2.0 * tr_xy_y[i] * tke_0 - tr_x_0[i];

        tr_0_0_y_xz_0[i] = 2.0 * tr_xyz_0[i] * tbe_0 + 2.0 * tr_xz_y[i] * tke_0;

        tr_0_0_y_yy_0[i] = 2.0 * tr_yyy_0[i] * tbe_0 + 2.0 * tr_yy_y[i] * tke_0 - 2.0 * tr_y_0[i];

        tr_0_0_y_yz_0[i] = 2.0 * tr_yyz_0[i] * tbe_0 + 2.0 * tr_yz_y[i] * tke_0 - tr_z_0[i];

        tr_0_0_y_zz_0[i] = 2.0 * tr_yzz_0[i] * tbe_0 + 2.0 * tr_zz_y[i] * tke_0;

        tr_0_0_z_xx_0[i] = 2.0 * tr_xxz_0[i] * tbe_0 + 2.0 * tr_xx_z[i] * tke_0;

        tr_0_0_z_xy_0[i] = 2.0 * tr_xyz_0[i] * tbe_0 + 2.0 * tr_xy_z[i] * tke_0;

        tr_0_0_z_xz_0[i] = 2.0 * tr_xzz_0[i] * tbe_0 + 2.0 * tr_xz_z[i] * tke_0 - tr_x_0[i];

        tr_0_0_z_yy_0[i] = 2.0 * tr_yyz_0[i] * tbe_0 + 2.0 * tr_yy_z[i] * tke_0;

        tr_0_0_z_yz_0[i] = 2.0 * tr_yzz_0[i] * tbe_0 + 2.0 * tr_yz_z[i] * tke_0 - tr_y_0[i];

        tr_0_0_z_zz_0[i] = 2.0 * tr_zzz_0[i] * tbe_0 + 2.0 * tr_zz_z[i] * tke_0 - 2.0 * tr_z_0[i];
    }
}

} // t2cgeom namespace

