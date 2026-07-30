#include "GeometricalDerivatives010ForSD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_sd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_sd,
                         const int idx_op_sp,
                         const int idx_op_sf,
                         const int idx_op_pd,
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

    // Set up components of auxiliary buffer : SF

    auto tr_0_xxx = pbuffer.data(idx_op_sf);

    auto tr_0_xxy = pbuffer.data(idx_op_sf + 1);

    auto tr_0_xxz = pbuffer.data(idx_op_sf + 2);

    auto tr_0_xyy = pbuffer.data(idx_op_sf + 3);

    auto tr_0_xyz = pbuffer.data(idx_op_sf + 4);

    auto tr_0_xzz = pbuffer.data(idx_op_sf + 5);

    auto tr_0_yyy = pbuffer.data(idx_op_sf + 6);

    auto tr_0_yyz = pbuffer.data(idx_op_sf + 7);

    auto tr_0_yzz = pbuffer.data(idx_op_sf + 8);

    auto tr_0_zzz = pbuffer.data(idx_op_sf + 9);

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

    // Set up components of targeted buffer : SD

    auto tr_0_0_x_0_xx = pbuffer.data(idx_op_geom_010_sd);

    auto tr_0_0_x_0_xy = pbuffer.data(idx_op_geom_010_sd + 1);

    auto tr_0_0_x_0_xz = pbuffer.data(idx_op_geom_010_sd + 2);

    auto tr_0_0_x_0_yy = pbuffer.data(idx_op_geom_010_sd + 3);

    auto tr_0_0_x_0_yz = pbuffer.data(idx_op_geom_010_sd + 4);

    auto tr_0_0_x_0_zz = pbuffer.data(idx_op_geom_010_sd + 5);

    auto tr_0_0_y_0_xx = pbuffer.data(idx_op_geom_010_sd + 6);

    auto tr_0_0_y_0_xy = pbuffer.data(idx_op_geom_010_sd + 7);

    auto tr_0_0_y_0_xz = pbuffer.data(idx_op_geom_010_sd + 8);

    auto tr_0_0_y_0_yy = pbuffer.data(idx_op_geom_010_sd + 9);

    auto tr_0_0_y_0_yz = pbuffer.data(idx_op_geom_010_sd + 10);

    auto tr_0_0_y_0_zz = pbuffer.data(idx_op_geom_010_sd + 11);

    auto tr_0_0_z_0_xx = pbuffer.data(idx_op_geom_010_sd + 12);

    auto tr_0_0_z_0_xy = pbuffer.data(idx_op_geom_010_sd + 13);

    auto tr_0_0_z_0_xz = pbuffer.data(idx_op_geom_010_sd + 14);

    auto tr_0_0_z_0_yy = pbuffer.data(idx_op_geom_010_sd + 15);

    auto tr_0_0_z_0_yz = pbuffer.data(idx_op_geom_010_sd + 16);

    auto tr_0_0_z_0_zz = pbuffer.data(idx_op_geom_010_sd + 17);

    #pragma omp simd aligned(tr_0_0_x_0_xx, tr_0_0_x_0_xy, tr_0_0_x_0_xz, tr_0_0_x_0_yy, tr_0_0_x_0_yz, tr_0_0_x_0_zz, tr_0_0_y_0_xx, tr_0_0_y_0_xy, tr_0_0_y_0_xz, tr_0_0_y_0_yy, tr_0_0_y_0_yz, tr_0_0_y_0_zz, tr_0_0_z_0_xx, tr_0_0_z_0_xy, tr_0_0_z_0_xz, tr_0_0_z_0_yy, tr_0_0_z_0_yz, tr_0_0_z_0_zz, tr_0_x, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_y, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_z, tr_0_zzz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_0_xx[i] = 2.0 * tr_x_xx[i] * tbe_0 + 2.0 * tr_0_xxx[i] * tke_0 - 2.0 * tr_0_x[i];

        tr_0_0_x_0_xy[i] = 2.0 * tr_x_xy[i] * tbe_0 + 2.0 * tr_0_xxy[i] * tke_0 - tr_0_y[i];

        tr_0_0_x_0_xz[i] = 2.0 * tr_x_xz[i] * tbe_0 + 2.0 * tr_0_xxz[i] * tke_0 - tr_0_z[i];

        tr_0_0_x_0_yy[i] = 2.0 * tr_x_yy[i] * tbe_0 + 2.0 * tr_0_xyy[i] * tke_0;

        tr_0_0_x_0_yz[i] = 2.0 * tr_x_yz[i] * tbe_0 + 2.0 * tr_0_xyz[i] * tke_0;

        tr_0_0_x_0_zz[i] = 2.0 * tr_x_zz[i] * tbe_0 + 2.0 * tr_0_xzz[i] * tke_0;

        tr_0_0_y_0_xx[i] = 2.0 * tr_y_xx[i] * tbe_0 + 2.0 * tr_0_xxy[i] * tke_0;

        tr_0_0_y_0_xy[i] = 2.0 * tr_y_xy[i] * tbe_0 + 2.0 * tr_0_xyy[i] * tke_0 - tr_0_x[i];

        tr_0_0_y_0_xz[i] = 2.0 * tr_y_xz[i] * tbe_0 + 2.0 * tr_0_xyz[i] * tke_0;

        tr_0_0_y_0_yy[i] = 2.0 * tr_y_yy[i] * tbe_0 + 2.0 * tr_0_yyy[i] * tke_0 - 2.0 * tr_0_y[i];

        tr_0_0_y_0_yz[i] = 2.0 * tr_y_yz[i] * tbe_0 + 2.0 * tr_0_yyz[i] * tke_0 - tr_0_z[i];

        tr_0_0_y_0_zz[i] = 2.0 * tr_y_zz[i] * tbe_0 + 2.0 * tr_0_yzz[i] * tke_0;

        tr_0_0_z_0_xx[i] = 2.0 * tr_z_xx[i] * tbe_0 + 2.0 * tr_0_xxz[i] * tke_0;

        tr_0_0_z_0_xy[i] = 2.0 * tr_z_xy[i] * tbe_0 + 2.0 * tr_0_xyz[i] * tke_0;

        tr_0_0_z_0_xz[i] = 2.0 * tr_z_xz[i] * tbe_0 + 2.0 * tr_0_xzz[i] * tke_0 - tr_0_x[i];

        tr_0_0_z_0_yy[i] = 2.0 * tr_z_yy[i] * tbe_0 + 2.0 * tr_0_yyz[i] * tke_0;

        tr_0_0_z_0_yz[i] = 2.0 * tr_z_yz[i] * tbe_0 + 2.0 * tr_0_yzz[i] * tke_0 - tr_0_y[i];

        tr_0_0_z_0_zz[i] = 2.0 * tr_z_zz[i] * tbe_0 + 2.0 * tr_0_zzz[i] * tke_0 - 2.0 * tr_0_z[i];
    }
}

} // t2cgeom namespace

