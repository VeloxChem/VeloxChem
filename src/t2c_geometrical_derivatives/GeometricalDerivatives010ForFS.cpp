#include "GeometricalDerivatives010ForFS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_fs(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_fs,
                         const int idx_op_ds,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : DS

    auto tr_xx_0 = pbuffer.data(idx_op_ds);

    auto tr_xy_0 = pbuffer.data(idx_op_ds + 1);

    auto tr_xz_0 = pbuffer.data(idx_op_ds + 2);

    auto tr_yy_0 = pbuffer.data(idx_op_ds + 3);

    auto tr_yz_0 = pbuffer.data(idx_op_ds + 4);

    auto tr_zz_0 = pbuffer.data(idx_op_ds + 5);

    // Set up components of auxiliary buffer : FP

    auto tr_xxx_x = pbuffer.data(idx_op_fp);

    auto tr_xxx_y = pbuffer.data(idx_op_fp + 1);

    auto tr_xxx_z = pbuffer.data(idx_op_fp + 2);

    auto tr_xxy_x = pbuffer.data(idx_op_fp + 3);

    auto tr_xxy_y = pbuffer.data(idx_op_fp + 4);

    auto tr_xxy_z = pbuffer.data(idx_op_fp + 5);

    auto tr_xxz_x = pbuffer.data(idx_op_fp + 6);

    auto tr_xxz_y = pbuffer.data(idx_op_fp + 7);

    auto tr_xxz_z = pbuffer.data(idx_op_fp + 8);

    auto tr_xyy_x = pbuffer.data(idx_op_fp + 9);

    auto tr_xyy_y = pbuffer.data(idx_op_fp + 10);

    auto tr_xyy_z = pbuffer.data(idx_op_fp + 11);

    auto tr_xyz_x = pbuffer.data(idx_op_fp + 12);

    auto tr_xyz_y = pbuffer.data(idx_op_fp + 13);

    auto tr_xyz_z = pbuffer.data(idx_op_fp + 14);

    auto tr_xzz_x = pbuffer.data(idx_op_fp + 15);

    auto tr_xzz_y = pbuffer.data(idx_op_fp + 16);

    auto tr_xzz_z = pbuffer.data(idx_op_fp + 17);

    auto tr_yyy_x = pbuffer.data(idx_op_fp + 18);

    auto tr_yyy_y = pbuffer.data(idx_op_fp + 19);

    auto tr_yyy_z = pbuffer.data(idx_op_fp + 20);

    auto tr_yyz_x = pbuffer.data(idx_op_fp + 21);

    auto tr_yyz_y = pbuffer.data(idx_op_fp + 22);

    auto tr_yyz_z = pbuffer.data(idx_op_fp + 23);

    auto tr_yzz_x = pbuffer.data(idx_op_fp + 24);

    auto tr_yzz_y = pbuffer.data(idx_op_fp + 25);

    auto tr_yzz_z = pbuffer.data(idx_op_fp + 26);

    auto tr_zzz_x = pbuffer.data(idx_op_fp + 27);

    auto tr_zzz_y = pbuffer.data(idx_op_fp + 28);

    auto tr_zzz_z = pbuffer.data(idx_op_fp + 29);

    // Set up components of auxiliary buffer : GS

    auto tr_xxxx_0 = pbuffer.data(idx_op_gs);

    auto tr_xxxy_0 = pbuffer.data(idx_op_gs + 1);

    auto tr_xxxz_0 = pbuffer.data(idx_op_gs + 2);

    auto tr_xxyy_0 = pbuffer.data(idx_op_gs + 3);

    auto tr_xxyz_0 = pbuffer.data(idx_op_gs + 4);

    auto tr_xxzz_0 = pbuffer.data(idx_op_gs + 5);

    auto tr_xyyy_0 = pbuffer.data(idx_op_gs + 6);

    auto tr_xyyz_0 = pbuffer.data(idx_op_gs + 7);

    auto tr_xyzz_0 = pbuffer.data(idx_op_gs + 8);

    auto tr_xzzz_0 = pbuffer.data(idx_op_gs + 9);

    auto tr_yyyy_0 = pbuffer.data(idx_op_gs + 10);

    auto tr_yyyz_0 = pbuffer.data(idx_op_gs + 11);

    auto tr_yyzz_0 = pbuffer.data(idx_op_gs + 12);

    auto tr_yzzz_0 = pbuffer.data(idx_op_gs + 13);

    auto tr_zzzz_0 = pbuffer.data(idx_op_gs + 14);

    // Set up components of targeted buffer : FS

    auto tr_0_0_x_xxx_0 = pbuffer.data(idx_op_geom_010_fs);

    auto tr_0_0_x_xxy_0 = pbuffer.data(idx_op_geom_010_fs + 1);

    auto tr_0_0_x_xxz_0 = pbuffer.data(idx_op_geom_010_fs + 2);

    auto tr_0_0_x_xyy_0 = pbuffer.data(idx_op_geom_010_fs + 3);

    auto tr_0_0_x_xyz_0 = pbuffer.data(idx_op_geom_010_fs + 4);

    auto tr_0_0_x_xzz_0 = pbuffer.data(idx_op_geom_010_fs + 5);

    auto tr_0_0_x_yyy_0 = pbuffer.data(idx_op_geom_010_fs + 6);

    auto tr_0_0_x_yyz_0 = pbuffer.data(idx_op_geom_010_fs + 7);

    auto tr_0_0_x_yzz_0 = pbuffer.data(idx_op_geom_010_fs + 8);

    auto tr_0_0_x_zzz_0 = pbuffer.data(idx_op_geom_010_fs + 9);

    auto tr_0_0_y_xxx_0 = pbuffer.data(idx_op_geom_010_fs + 10);

    auto tr_0_0_y_xxy_0 = pbuffer.data(idx_op_geom_010_fs + 11);

    auto tr_0_0_y_xxz_0 = pbuffer.data(idx_op_geom_010_fs + 12);

    auto tr_0_0_y_xyy_0 = pbuffer.data(idx_op_geom_010_fs + 13);

    auto tr_0_0_y_xyz_0 = pbuffer.data(idx_op_geom_010_fs + 14);

    auto tr_0_0_y_xzz_0 = pbuffer.data(idx_op_geom_010_fs + 15);

    auto tr_0_0_y_yyy_0 = pbuffer.data(idx_op_geom_010_fs + 16);

    auto tr_0_0_y_yyz_0 = pbuffer.data(idx_op_geom_010_fs + 17);

    auto tr_0_0_y_yzz_0 = pbuffer.data(idx_op_geom_010_fs + 18);

    auto tr_0_0_y_zzz_0 = pbuffer.data(idx_op_geom_010_fs + 19);

    auto tr_0_0_z_xxx_0 = pbuffer.data(idx_op_geom_010_fs + 20);

    auto tr_0_0_z_xxy_0 = pbuffer.data(idx_op_geom_010_fs + 21);

    auto tr_0_0_z_xxz_0 = pbuffer.data(idx_op_geom_010_fs + 22);

    auto tr_0_0_z_xyy_0 = pbuffer.data(idx_op_geom_010_fs + 23);

    auto tr_0_0_z_xyz_0 = pbuffer.data(idx_op_geom_010_fs + 24);

    auto tr_0_0_z_xzz_0 = pbuffer.data(idx_op_geom_010_fs + 25);

    auto tr_0_0_z_yyy_0 = pbuffer.data(idx_op_geom_010_fs + 26);

    auto tr_0_0_z_yyz_0 = pbuffer.data(idx_op_geom_010_fs + 27);

    auto tr_0_0_z_yzz_0 = pbuffer.data(idx_op_geom_010_fs + 28);

    auto tr_0_0_z_zzz_0 = pbuffer.data(idx_op_geom_010_fs + 29);

    #pragma omp simd aligned(tr_0_0_x_xxx_0, tr_0_0_x_xxy_0, tr_0_0_x_xxz_0, tr_0_0_x_xyy_0, tr_0_0_x_xyz_0, tr_0_0_x_xzz_0, tr_0_0_x_yyy_0, tr_0_0_x_yyz_0, tr_0_0_x_yzz_0, tr_0_0_x_zzz_0, tr_0_0_y_xxx_0, tr_0_0_y_xxy_0, tr_0_0_y_xxz_0, tr_0_0_y_xyy_0, tr_0_0_y_xyz_0, tr_0_0_y_xzz_0, tr_0_0_y_yyy_0, tr_0_0_y_yyz_0, tr_0_0_y_yzz_0, tr_0_0_y_zzz_0, tr_0_0_z_xxx_0, tr_0_0_z_xxy_0, tr_0_0_z_xxz_0, tr_0_0_z_xyy_0, tr_0_0_z_xyz_0, tr_0_0_z_xzz_0, tr_0_0_z_yyy_0, tr_0_0_z_yyz_0, tr_0_0_z_yzz_0, tr_0_0_z_zzz_0, tr_xx_0, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxy_0, tr_xxxz_0, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyz_0, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xy_0, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyz_0, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xz_0, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_yy_0, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyy_0, tr_yyyz_0, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yz_0, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_zz_0, tr_zzz_x, tr_zzz_y, tr_zzz_z, tr_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxx_0[i] = 2.0 * tr_xxxx_0[i] * tbe_0 + 2.0 * tr_xxx_x[i] * tke_0 - 3.0 * tr_xx_0[i];

        tr_0_0_x_xxy_0[i] = 2.0 * tr_xxxy_0[i] * tbe_0 + 2.0 * tr_xxy_x[i] * tke_0 - 2.0 * tr_xy_0[i];

        tr_0_0_x_xxz_0[i] = 2.0 * tr_xxxz_0[i] * tbe_0 + 2.0 * tr_xxz_x[i] * tke_0 - 2.0 * tr_xz_0[i];

        tr_0_0_x_xyy_0[i] = 2.0 * tr_xxyy_0[i] * tbe_0 + 2.0 * tr_xyy_x[i] * tke_0 - tr_yy_0[i];

        tr_0_0_x_xyz_0[i] = 2.0 * tr_xxyz_0[i] * tbe_0 + 2.0 * tr_xyz_x[i] * tke_0 - tr_yz_0[i];

        tr_0_0_x_xzz_0[i] = 2.0 * tr_xxzz_0[i] * tbe_0 + 2.0 * tr_xzz_x[i] * tke_0 - tr_zz_0[i];

        tr_0_0_x_yyy_0[i] = 2.0 * tr_xyyy_0[i] * tbe_0 + 2.0 * tr_yyy_x[i] * tke_0;

        tr_0_0_x_yyz_0[i] = 2.0 * tr_xyyz_0[i] * tbe_0 + 2.0 * tr_yyz_x[i] * tke_0;

        tr_0_0_x_yzz_0[i] = 2.0 * tr_xyzz_0[i] * tbe_0 + 2.0 * tr_yzz_x[i] * tke_0;

        tr_0_0_x_zzz_0[i] = 2.0 * tr_xzzz_0[i] * tbe_0 + 2.0 * tr_zzz_x[i] * tke_0;

        tr_0_0_y_xxx_0[i] = 2.0 * tr_xxxy_0[i] * tbe_0 + 2.0 * tr_xxx_y[i] * tke_0;

        tr_0_0_y_xxy_0[i] = 2.0 * tr_xxyy_0[i] * tbe_0 + 2.0 * tr_xxy_y[i] * tke_0 - tr_xx_0[i];

        tr_0_0_y_xxz_0[i] = 2.0 * tr_xxyz_0[i] * tbe_0 + 2.0 * tr_xxz_y[i] * tke_0;

        tr_0_0_y_xyy_0[i] = 2.0 * tr_xyyy_0[i] * tbe_0 + 2.0 * tr_xyy_y[i] * tke_0 - 2.0 * tr_xy_0[i];

        tr_0_0_y_xyz_0[i] = 2.0 * tr_xyyz_0[i] * tbe_0 + 2.0 * tr_xyz_y[i] * tke_0 - tr_xz_0[i];

        tr_0_0_y_xzz_0[i] = 2.0 * tr_xyzz_0[i] * tbe_0 + 2.0 * tr_xzz_y[i] * tke_0;

        tr_0_0_y_yyy_0[i] = 2.0 * tr_yyyy_0[i] * tbe_0 + 2.0 * tr_yyy_y[i] * tke_0 - 3.0 * tr_yy_0[i];

        tr_0_0_y_yyz_0[i] = 2.0 * tr_yyyz_0[i] * tbe_0 + 2.0 * tr_yyz_y[i] * tke_0 - 2.0 * tr_yz_0[i];

        tr_0_0_y_yzz_0[i] = 2.0 * tr_yyzz_0[i] * tbe_0 + 2.0 * tr_yzz_y[i] * tke_0 - tr_zz_0[i];

        tr_0_0_y_zzz_0[i] = 2.0 * tr_yzzz_0[i] * tbe_0 + 2.0 * tr_zzz_y[i] * tke_0;

        tr_0_0_z_xxx_0[i] = 2.0 * tr_xxxz_0[i] * tbe_0 + 2.0 * tr_xxx_z[i] * tke_0;

        tr_0_0_z_xxy_0[i] = 2.0 * tr_xxyz_0[i] * tbe_0 + 2.0 * tr_xxy_z[i] * tke_0;

        tr_0_0_z_xxz_0[i] = 2.0 * tr_xxzz_0[i] * tbe_0 + 2.0 * tr_xxz_z[i] * tke_0 - tr_xx_0[i];

        tr_0_0_z_xyy_0[i] = 2.0 * tr_xyyz_0[i] * tbe_0 + 2.0 * tr_xyy_z[i] * tke_0;

        tr_0_0_z_xyz_0[i] = 2.0 * tr_xyzz_0[i] * tbe_0 + 2.0 * tr_xyz_z[i] * tke_0 - tr_xy_0[i];

        tr_0_0_z_xzz_0[i] = 2.0 * tr_xzzz_0[i] * tbe_0 + 2.0 * tr_xzz_z[i] * tke_0 - 2.0 * tr_xz_0[i];

        tr_0_0_z_yyy_0[i] = 2.0 * tr_yyyz_0[i] * tbe_0 + 2.0 * tr_yyy_z[i] * tke_0;

        tr_0_0_z_yyz_0[i] = 2.0 * tr_yyzz_0[i] * tbe_0 + 2.0 * tr_yyz_z[i] * tke_0 - tr_yy_0[i];

        tr_0_0_z_yzz_0[i] = 2.0 * tr_yzzz_0[i] * tbe_0 + 2.0 * tr_yzz_z[i] * tke_0 - 2.0 * tr_yz_0[i];

        tr_0_0_z_zzz_0[i] = 2.0 * tr_zzzz_0[i] * tbe_0 + 2.0 * tr_zzz_z[i] * tke_0 - 3.0 * tr_zz_0[i];
    }
}

} // t2cgeom namespace

