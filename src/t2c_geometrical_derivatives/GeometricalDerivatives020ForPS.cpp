#include "GeometricalDerivatives020ForPS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_ps(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_ps,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const int idx_op_fs,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SP

    auto tr_0_x = pbuffer.data(idx_op_sp);

    auto tr_0_y = pbuffer.data(idx_op_sp + 1);

    auto tr_0_z = pbuffer.data(idx_op_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto tr_x_0 = pbuffer.data(idx_op_ps);

    auto tr_y_0 = pbuffer.data(idx_op_ps + 1);

    auto tr_z_0 = pbuffer.data(idx_op_ps + 2);

    // Set up components of auxiliary buffer : PD

    auto tr_x_xx = pbuffer.data(idx_op_pd);

    auto tr_x_xy = pbuffer.data(idx_op_pd + 1);

    auto tr_x_xz = pbuffer.data(idx_op_pd + 2);

    auto tr_x_yy = pbuffer.data(idx_op_pd + 3);

    auto tr_x_yz = pbuffer.data(idx_op_pd + 4);

    auto tr_x_zz = pbuffer.data(idx_op_pd + 5);

    auto tr_y_xx = pbuffer.data(idx_op_pd + 6);

    auto tr_y_xy = pbuffer.data(idx_op_pd + 7);

    auto tr_y_xz = pbuffer.data(idx_op_pd + 8);

    auto tr_y_yy = pbuffer.data(idx_op_pd + 9);

    auto tr_y_yz = pbuffer.data(idx_op_pd + 10);

    auto tr_y_zz = pbuffer.data(idx_op_pd + 11);

    auto tr_z_xx = pbuffer.data(idx_op_pd + 12);

    auto tr_z_xy = pbuffer.data(idx_op_pd + 13);

    auto tr_z_xz = pbuffer.data(idx_op_pd + 14);

    auto tr_z_yy = pbuffer.data(idx_op_pd + 15);

    auto tr_z_yz = pbuffer.data(idx_op_pd + 16);

    auto tr_z_zz = pbuffer.data(idx_op_pd + 17);

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

    // Set up components of targeted buffer : PS

    auto tr_0_0_xx_x_0 = pbuffer.data(idx_op_geom_020_ps);

    auto tr_0_0_xx_y_0 = pbuffer.data(idx_op_geom_020_ps + 1);

    auto tr_0_0_xx_z_0 = pbuffer.data(idx_op_geom_020_ps + 2);

    auto tr_0_0_xy_x_0 = pbuffer.data(idx_op_geom_020_ps + 3);

    auto tr_0_0_xy_y_0 = pbuffer.data(idx_op_geom_020_ps + 4);

    auto tr_0_0_xy_z_0 = pbuffer.data(idx_op_geom_020_ps + 5);

    auto tr_0_0_xz_x_0 = pbuffer.data(idx_op_geom_020_ps + 6);

    auto tr_0_0_xz_y_0 = pbuffer.data(idx_op_geom_020_ps + 7);

    auto tr_0_0_xz_z_0 = pbuffer.data(idx_op_geom_020_ps + 8);

    auto tr_0_0_yy_x_0 = pbuffer.data(idx_op_geom_020_ps + 9);

    auto tr_0_0_yy_y_0 = pbuffer.data(idx_op_geom_020_ps + 10);

    auto tr_0_0_yy_z_0 = pbuffer.data(idx_op_geom_020_ps + 11);

    auto tr_0_0_yz_x_0 = pbuffer.data(idx_op_geom_020_ps + 12);

    auto tr_0_0_yz_y_0 = pbuffer.data(idx_op_geom_020_ps + 13);

    auto tr_0_0_yz_z_0 = pbuffer.data(idx_op_geom_020_ps + 14);

    auto tr_0_0_zz_x_0 = pbuffer.data(idx_op_geom_020_ps + 15);

    auto tr_0_0_zz_y_0 = pbuffer.data(idx_op_geom_020_ps + 16);

    auto tr_0_0_zz_z_0 = pbuffer.data(idx_op_geom_020_ps + 17);

    #pragma omp simd aligned(tr_0_0_xx_x_0, tr_0_0_xx_y_0, tr_0_0_xx_z_0, tr_0_0_xy_x_0, tr_0_0_xy_y_0, tr_0_0_xy_z_0, tr_0_0_xz_x_0, tr_0_0_xz_y_0, tr_0_0_xz_z_0, tr_0_0_yy_x_0, tr_0_0_yy_y_0, tr_0_0_yy_z_0, tr_0_0_yz_x_0, tr_0_0_yz_y_0, tr_0_0_yz_z_0, tr_0_0_zz_x_0, tr_0_0_zz_y_0, tr_0_0_zz_z_0, tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_y, tr_xx_z, tr_xxx_0, tr_xxy_0, tr_xxz_0, tr_xy_x, tr_xy_y, tr_xy_z, tr_xyy_0, tr_xyz_0, tr_xz_x, tr_xz_y, tr_xz_z, tr_xzz_0, tr_y_0, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_y, tr_yy_z, tr_yyy_0, tr_yyz_0, tr_yz_x, tr_yz_y, tr_yz_z, tr_yzz_0, tr_z_0, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_y, tr_zz_z, tr_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_x_0[i] = -4.0 * tr_0_x[i] * tke_0 - 6.0 * tr_x_0[i] * tbe_0 - 2.0 * tr_x_0[i] * tke_0 + 4.0 * tr_x_xx[i] * tke_0 * tke_0 + 8.0 * tr_xx_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_0[i] = -2.0 * tr_y_0[i] * tbe_0 - 2.0 * tr_y_0[i] * tke_0 + 4.0 * tr_y_xx[i] * tke_0 * tke_0 + 8.0 * tr_xy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_0[i] = -2.0 * tr_z_0[i] * tbe_0 - 2.0 * tr_z_0[i] * tke_0 + 4.0 * tr_z_xx[i] * tke_0 * tke_0 + 8.0 * tr_xz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_0[i] = -2.0 * tr_0_y[i] * tke_0 - 2.0 * tr_y_0[i] * tbe_0 + 4.0 * tr_x_xy[i] * tke_0 * tke_0 + 4.0 * tr_xy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xx_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_0[i] = -2.0 * tr_0_x[i] * tke_0 + 4.0 * tr_y_xy[i] * tke_0 * tke_0 + 4.0 * tr_yy_x[i] * tbe_0 * tke_0 - 2.0 * tr_x_0[i] * tbe_0 + 4.0 * tr_xy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_0[i] = 4.0 * tr_z_xy[i] * tke_0 * tke_0 + 4.0 * tr_yz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_0[i] = -2.0 * tr_0_z[i] * tke_0 - 2.0 * tr_z_0[i] * tbe_0 + 4.0 * tr_x_xz[i] * tke_0 * tke_0 + 4.0 * tr_xz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xx_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_0[i] = 4.0 * tr_y_xz[i] * tke_0 * tke_0 + 4.0 * tr_yz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_0[i] = -2.0 * tr_0_x[i] * tke_0 + 4.0 * tr_z_xz[i] * tke_0 * tke_0 + 4.0 * tr_zz_x[i] * tbe_0 * tke_0 - 2.0 * tr_x_0[i] * tbe_0 + 4.0 * tr_xz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_0[i] = -2.0 * tr_x_0[i] * tbe_0 - 2.0 * tr_x_0[i] * tke_0 + 4.0 * tr_x_yy[i] * tke_0 * tke_0 + 8.0 * tr_xy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_0[i] = -4.0 * tr_0_y[i] * tke_0 - 6.0 * tr_y_0[i] * tbe_0 - 2.0 * tr_y_0[i] * tke_0 + 4.0 * tr_y_yy[i] * tke_0 * tke_0 + 8.0 * tr_yy_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_0[i] = -2.0 * tr_z_0[i] * tbe_0 - 2.0 * tr_z_0[i] * tke_0 + 4.0 * tr_z_yy[i] * tke_0 * tke_0 + 8.0 * tr_yz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_0[i] = 4.0 * tr_x_yz[i] * tke_0 * tke_0 + 4.0 * tr_xz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_0[i] = -2.0 * tr_0_z[i] * tke_0 - 2.0 * tr_z_0[i] * tbe_0 + 4.0 * tr_y_yz[i] * tke_0 * tke_0 + 4.0 * tr_yz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yy_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_0[i] = -2.0 * tr_0_y[i] * tke_0 + 4.0 * tr_z_yz[i] * tke_0 * tke_0 + 4.0 * tr_zz_y[i] * tbe_0 * tke_0 - 2.0 * tr_y_0[i] * tbe_0 + 4.0 * tr_yz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_0[i] = -2.0 * tr_x_0[i] * tbe_0 - 2.0 * tr_x_0[i] * tke_0 + 4.0 * tr_x_zz[i] * tke_0 * tke_0 + 8.0 * tr_xz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_0[i] = -2.0 * tr_y_0[i] * tbe_0 - 2.0 * tr_y_0[i] * tke_0 + 4.0 * tr_y_zz[i] * tke_0 * tke_0 + 8.0 * tr_yz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_0[i] = -4.0 * tr_0_z[i] * tke_0 - 6.0 * tr_z_0[i] * tbe_0 - 2.0 * tr_z_0[i] * tke_0 + 4.0 * tr_z_zz[i] * tke_0 * tke_0 + 8.0 * tr_zz_z[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_0[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

