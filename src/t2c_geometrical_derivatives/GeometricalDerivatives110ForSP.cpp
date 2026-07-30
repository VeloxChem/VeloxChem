#include "GeometricalDerivatives110ForSP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_sp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_sp,
                         const int idx_op_sp,
                         const int idx_op_ps,
                         const int idx_op_pd,
                         const int idx_op_dp,
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

    // Set up components of targeted buffer : SP

    auto tr_x_0_x_0_x = pbuffer.data(idx_op_geom_110_sp);

    auto tr_x_0_x_0_y = pbuffer.data(idx_op_geom_110_sp + 1);

    auto tr_x_0_x_0_z = pbuffer.data(idx_op_geom_110_sp + 2);

    auto tr_x_0_y_0_x = pbuffer.data(idx_op_geom_110_sp + 3);

    auto tr_x_0_y_0_y = pbuffer.data(idx_op_geom_110_sp + 4);

    auto tr_x_0_y_0_z = pbuffer.data(idx_op_geom_110_sp + 5);

    auto tr_x_0_z_0_x = pbuffer.data(idx_op_geom_110_sp + 6);

    auto tr_x_0_z_0_y = pbuffer.data(idx_op_geom_110_sp + 7);

    auto tr_x_0_z_0_z = pbuffer.data(idx_op_geom_110_sp + 8);

    auto tr_y_0_x_0_x = pbuffer.data(idx_op_geom_110_sp + 9);

    auto tr_y_0_x_0_y = pbuffer.data(idx_op_geom_110_sp + 10);

    auto tr_y_0_x_0_z = pbuffer.data(idx_op_geom_110_sp + 11);

    auto tr_y_0_y_0_x = pbuffer.data(idx_op_geom_110_sp + 12);

    auto tr_y_0_y_0_y = pbuffer.data(idx_op_geom_110_sp + 13);

    auto tr_y_0_y_0_z = pbuffer.data(idx_op_geom_110_sp + 14);

    auto tr_y_0_z_0_x = pbuffer.data(idx_op_geom_110_sp + 15);

    auto tr_y_0_z_0_y = pbuffer.data(idx_op_geom_110_sp + 16);

    auto tr_y_0_z_0_z = pbuffer.data(idx_op_geom_110_sp + 17);

    auto tr_z_0_x_0_x = pbuffer.data(idx_op_geom_110_sp + 18);

    auto tr_z_0_x_0_y = pbuffer.data(idx_op_geom_110_sp + 19);

    auto tr_z_0_x_0_z = pbuffer.data(idx_op_geom_110_sp + 20);

    auto tr_z_0_y_0_x = pbuffer.data(idx_op_geom_110_sp + 21);

    auto tr_z_0_y_0_y = pbuffer.data(idx_op_geom_110_sp + 22);

    auto tr_z_0_y_0_z = pbuffer.data(idx_op_geom_110_sp + 23);

    auto tr_z_0_z_0_x = pbuffer.data(idx_op_geom_110_sp + 24);

    auto tr_z_0_z_0_y = pbuffer.data(idx_op_geom_110_sp + 25);

    auto tr_z_0_z_0_z = pbuffer.data(idx_op_geom_110_sp + 26);

    #pragma omp simd aligned(tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_x_0_x_0_x, tr_x_0_x_0_y, tr_x_0_x_0_z, tr_x_0_y_0_x, tr_x_0_y_0_y, tr_x_0_y_0_z, tr_x_0_z_0_x, tr_x_0_z_0_y, tr_x_0_z_0_z, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_y, tr_xx_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_y_0, tr_y_0_x_0_x, tr_y_0_x_0_y, tr_y_0_x_0_z, tr_y_0_y_0_x, tr_y_0_y_0_y, tr_y_0_y_0_z, tr_y_0_z_0_x, tr_y_0_z_0_y, tr_y_0_z_0_z, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_y, tr_yy_z, tr_yz_x, tr_yz_y, tr_yz_z, tr_z_0, tr_z_0_x_0_x, tr_z_0_x_0_y, tr_z_0_x_0_z, tr_z_0_y_0_x, tr_z_0_y_0_y, tr_z_0_y_0_z, tr_z_0_z_0_x, tr_z_0_z_0_y, tr_z_0_z_0_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_0_x[i] = -2.0 * tr_0_x[i] * tbe_0 - 2.0 * tr_x_0[i] * tbe_0 + 4.0 * tr_x_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_x[i] * tbe_0 * tbe_0;

        tr_x_0_x_0_y[i] = -2.0 * tr_0_y[i] * tbe_0 + 4.0 * tr_x_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_y[i] * tbe_0 * tbe_0;

        tr_x_0_x_0_z[i] = -2.0 * tr_0_z[i] * tbe_0 + 4.0 * tr_x_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_z[i] * tbe_0 * tbe_0;

        tr_x_0_y_0_x[i] = 4.0 * tr_x_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_x[i] * tbe_0 * tbe_0;

        tr_x_0_y_0_y[i] = -2.0 * tr_x_0[i] * tbe_0 + 4.0 * tr_x_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_y[i] * tbe_0 * tbe_0;

        tr_x_0_y_0_z[i] = 4.0 * tr_x_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_z[i] * tbe_0 * tbe_0;

        tr_x_0_z_0_x[i] = 4.0 * tr_x_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_x[i] * tbe_0 * tbe_0;

        tr_x_0_z_0_y[i] = 4.0 * tr_x_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_y[i] * tbe_0 * tbe_0;

        tr_x_0_z_0_z[i] = -2.0 * tr_x_0[i] * tbe_0 + 4.0 * tr_x_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_z[i] * tbe_0 * tbe_0;

        tr_y_0_x_0_x[i] = -2.0 * tr_y_0[i] * tbe_0 + 4.0 * tr_y_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xy_x[i] * tbe_0 * tbe_0;

        tr_y_0_x_0_y[i] = 4.0 * tr_y_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_y[i] * tbe_0 * tbe_0;

        tr_y_0_x_0_z[i] = 4.0 * tr_y_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_z[i] * tbe_0 * tbe_0;

        tr_y_0_y_0_x[i] = -2.0 * tr_0_x[i] * tbe_0 + 4.0 * tr_y_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_x[i] * tbe_0 * tbe_0;

        tr_y_0_y_0_y[i] = -2.0 * tr_0_y[i] * tbe_0 - 2.0 * tr_y_0[i] * tbe_0 + 4.0 * tr_y_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_y[i] * tbe_0 * tbe_0;

        tr_y_0_y_0_z[i] = -2.0 * tr_0_z[i] * tbe_0 + 4.0 * tr_y_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_z[i] * tbe_0 * tbe_0;

        tr_y_0_z_0_x[i] = 4.0 * tr_y_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_x[i] * tbe_0 * tbe_0;

        tr_y_0_z_0_y[i] = 4.0 * tr_y_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_y[i] * tbe_0 * tbe_0;

        tr_y_0_z_0_z[i] = -2.0 * tr_y_0[i] * tbe_0 + 4.0 * tr_y_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_z[i] * tbe_0 * tbe_0;

        tr_z_0_x_0_x[i] = -2.0 * tr_z_0[i] * tbe_0 + 4.0 * tr_z_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xz_x[i] * tbe_0 * tbe_0;

        tr_z_0_x_0_y[i] = 4.0 * tr_z_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xz_y[i] * tbe_0 * tbe_0;

        tr_z_0_x_0_z[i] = 4.0 * tr_z_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_z[i] * tbe_0 * tbe_0;

        tr_z_0_y_0_x[i] = 4.0 * tr_z_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yz_x[i] * tbe_0 * tbe_0;

        tr_z_0_y_0_y[i] = -2.0 * tr_z_0[i] * tbe_0 + 4.0 * tr_z_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yz_y[i] * tbe_0 * tbe_0;

        tr_z_0_y_0_z[i] = 4.0 * tr_z_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_z[i] * tbe_0 * tbe_0;

        tr_z_0_z_0_x[i] = -2.0 * tr_0_x[i] * tbe_0 + 4.0 * tr_z_xz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_x[i] * tbe_0 * tbe_0;

        tr_z_0_z_0_y[i] = -2.0 * tr_0_y[i] * tbe_0 + 4.0 * tr_z_yz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_y[i] * tbe_0 * tbe_0;

        tr_z_0_z_0_z[i] = -2.0 * tr_0_z[i] * tbe_0 - 2.0 * tr_z_0[i] * tbe_0 + 4.0 * tr_z_zz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_z[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

