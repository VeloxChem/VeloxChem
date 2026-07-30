#include "GeometricalDerivatives010ForPP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_pp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_pp,
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

    // Set up 0-3 components of targeted buffer : PP

    auto tr_0_0_x_x_x = pbuffer.data(idx_op_geom_010_pp);

    auto tr_0_0_x_x_y = pbuffer.data(idx_op_geom_010_pp + 1);

    auto tr_0_0_x_x_z = pbuffer.data(idx_op_geom_010_pp + 2);

    #pragma omp simd aligned(tr_0_0_x_x_x, tr_0_0_x_x_y, tr_0_0_x_x_z, tr_0_x, tr_0_y, tr_0_z, tr_x_0, tr_x_xx, tr_x_xy, tr_x_xz, tr_xx_x, tr_xx_y, tr_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_x_x[i] = 2.0 * tr_xx_x[i] * tbe_0 + 2.0 * tr_x_xx[i] * tke_0 - tr_0_x[i] - tr_x_0[i];

        tr_0_0_x_x_y[i] = 2.0 * tr_xx_y[i] * tbe_0 + 2.0 * tr_x_xy[i] * tke_0 - tr_0_y[i];

        tr_0_0_x_x_z[i] = 2.0 * tr_xx_z[i] * tbe_0 + 2.0 * tr_x_xz[i] * tke_0 - tr_0_z[i];
    }

    // Set up 3-6 components of targeted buffer : PP

    auto tr_0_0_x_y_x = pbuffer.data(idx_op_geom_010_pp + 3);

    auto tr_0_0_x_y_y = pbuffer.data(idx_op_geom_010_pp + 4);

    auto tr_0_0_x_y_z = pbuffer.data(idx_op_geom_010_pp + 5);

    #pragma omp simd aligned(tr_0_0_x_y_x, tr_0_0_x_y_y, tr_0_0_x_y_z, tr_xy_x, tr_xy_y, tr_xy_z, tr_y_0, tr_y_xx, tr_y_xy, tr_y_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_y_x[i] = 2.0 * tr_xy_x[i] * tbe_0 + 2.0 * tr_y_xx[i] * tke_0 - tr_y_0[i];

        tr_0_0_x_y_y[i] = 2.0 * tr_xy_y[i] * tbe_0 + 2.0 * tr_y_xy[i] * tke_0;

        tr_0_0_x_y_z[i] = 2.0 * tr_xy_z[i] * tbe_0 + 2.0 * tr_y_xz[i] * tke_0;
    }

    // Set up 6-9 components of targeted buffer : PP

    auto tr_0_0_x_z_x = pbuffer.data(idx_op_geom_010_pp + 6);

    auto tr_0_0_x_z_y = pbuffer.data(idx_op_geom_010_pp + 7);

    auto tr_0_0_x_z_z = pbuffer.data(idx_op_geom_010_pp + 8);

    #pragma omp simd aligned(tr_0_0_x_z_x, tr_0_0_x_z_y, tr_0_0_x_z_z, tr_xz_x, tr_xz_y, tr_xz_z, tr_z_0, tr_z_xx, tr_z_xy, tr_z_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_z_x[i] = 2.0 * tr_xz_x[i] * tbe_0 + 2.0 * tr_z_xx[i] * tke_0 - tr_z_0[i];

        tr_0_0_x_z_y[i] = 2.0 * tr_xz_y[i] * tbe_0 + 2.0 * tr_z_xy[i] * tke_0;

        tr_0_0_x_z_z[i] = 2.0 * tr_xz_z[i] * tbe_0 + 2.0 * tr_z_xz[i] * tke_0;
    }

    // Set up 9-12 components of targeted buffer : PP

    auto tr_0_0_y_x_x = pbuffer.data(idx_op_geom_010_pp + 9);

    auto tr_0_0_y_x_y = pbuffer.data(idx_op_geom_010_pp + 10);

    auto tr_0_0_y_x_z = pbuffer.data(idx_op_geom_010_pp + 11);

    #pragma omp simd aligned(tr_0_0_y_x_x, tr_0_0_y_x_y, tr_0_0_y_x_z, tr_x_0, tr_x_xy, tr_x_yy, tr_x_yz, tr_xy_x, tr_xy_y, tr_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_x_x[i] = 2.0 * tr_xy_x[i] * tbe_0 + 2.0 * tr_x_xy[i] * tke_0;

        tr_0_0_y_x_y[i] = 2.0 * tr_xy_y[i] * tbe_0 + 2.0 * tr_x_yy[i] * tke_0 - tr_x_0[i];

        tr_0_0_y_x_z[i] = 2.0 * tr_xy_z[i] * tbe_0 + 2.0 * tr_x_yz[i] * tke_0;
    }

    // Set up 12-15 components of targeted buffer : PP

    auto tr_0_0_y_y_x = pbuffer.data(idx_op_geom_010_pp + 12);

    auto tr_0_0_y_y_y = pbuffer.data(idx_op_geom_010_pp + 13);

    auto tr_0_0_y_y_z = pbuffer.data(idx_op_geom_010_pp + 14);

    #pragma omp simd aligned(tr_0_0_y_y_x, tr_0_0_y_y_y, tr_0_0_y_y_z, tr_0_x, tr_0_y, tr_0_z, tr_y_0, tr_y_xy, tr_y_yy, tr_y_yz, tr_yy_x, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_y_x[i] = 2.0 * tr_yy_x[i] * tbe_0 + 2.0 * tr_y_xy[i] * tke_0 - tr_0_x[i];

        tr_0_0_y_y_y[i] = 2.0 * tr_yy_y[i] * tbe_0 + 2.0 * tr_y_yy[i] * tke_0 - tr_0_y[i] - tr_y_0[i];

        tr_0_0_y_y_z[i] = 2.0 * tr_yy_z[i] * tbe_0 + 2.0 * tr_y_yz[i] * tke_0 - tr_0_z[i];
    }

    // Set up 15-18 components of targeted buffer : PP

    auto tr_0_0_y_z_x = pbuffer.data(idx_op_geom_010_pp + 15);

    auto tr_0_0_y_z_y = pbuffer.data(idx_op_geom_010_pp + 16);

    auto tr_0_0_y_z_z = pbuffer.data(idx_op_geom_010_pp + 17);

    #pragma omp simd aligned(tr_0_0_y_z_x, tr_0_0_y_z_y, tr_0_0_y_z_z, tr_yz_x, tr_yz_y, tr_yz_z, tr_z_0, tr_z_xy, tr_z_yy, tr_z_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_z_x[i] = 2.0 * tr_yz_x[i] * tbe_0 + 2.0 * tr_z_xy[i] * tke_0;

        tr_0_0_y_z_y[i] = 2.0 * tr_yz_y[i] * tbe_0 + 2.0 * tr_z_yy[i] * tke_0 - tr_z_0[i];

        tr_0_0_y_z_z[i] = 2.0 * tr_yz_z[i] * tbe_0 + 2.0 * tr_z_yz[i] * tke_0;
    }

    // Set up 18-21 components of targeted buffer : PP

    auto tr_0_0_z_x_x = pbuffer.data(idx_op_geom_010_pp + 18);

    auto tr_0_0_z_x_y = pbuffer.data(idx_op_geom_010_pp + 19);

    auto tr_0_0_z_x_z = pbuffer.data(idx_op_geom_010_pp + 20);

    #pragma omp simd aligned(tr_0_0_z_x_x, tr_0_0_z_x_y, tr_0_0_z_x_z, tr_x_0, tr_x_xz, tr_x_yz, tr_x_zz, tr_xz_x, tr_xz_y, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_x_x[i] = 2.0 * tr_xz_x[i] * tbe_0 + 2.0 * tr_x_xz[i] * tke_0;

        tr_0_0_z_x_y[i] = 2.0 * tr_xz_y[i] * tbe_0 + 2.0 * tr_x_yz[i] * tke_0;

        tr_0_0_z_x_z[i] = 2.0 * tr_xz_z[i] * tbe_0 + 2.0 * tr_x_zz[i] * tke_0 - tr_x_0[i];
    }

    // Set up 21-24 components of targeted buffer : PP

    auto tr_0_0_z_y_x = pbuffer.data(idx_op_geom_010_pp + 21);

    auto tr_0_0_z_y_y = pbuffer.data(idx_op_geom_010_pp + 22);

    auto tr_0_0_z_y_z = pbuffer.data(idx_op_geom_010_pp + 23);

    #pragma omp simd aligned(tr_0_0_z_y_x, tr_0_0_z_y_y, tr_0_0_z_y_z, tr_y_0, tr_y_xz, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_y_x[i] = 2.0 * tr_yz_x[i] * tbe_0 + 2.0 * tr_y_xz[i] * tke_0;

        tr_0_0_z_y_y[i] = 2.0 * tr_yz_y[i] * tbe_0 + 2.0 * tr_y_yz[i] * tke_0;

        tr_0_0_z_y_z[i] = 2.0 * tr_yz_z[i] * tbe_0 + 2.0 * tr_y_zz[i] * tke_0 - tr_y_0[i];
    }

    // Set up 24-27 components of targeted buffer : PP

    auto tr_0_0_z_z_x = pbuffer.data(idx_op_geom_010_pp + 24);

    auto tr_0_0_z_z_y = pbuffer.data(idx_op_geom_010_pp + 25);

    auto tr_0_0_z_z_z = pbuffer.data(idx_op_geom_010_pp + 26);

    #pragma omp simd aligned(tr_0_0_z_z_x, tr_0_0_z_z_y, tr_0_0_z_z_z, tr_0_x, tr_0_y, tr_0_z, tr_z_0, tr_z_xz, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_z_x[i] = 2.0 * tr_zz_x[i] * tbe_0 + 2.0 * tr_z_xz[i] * tke_0 - tr_0_x[i];

        tr_0_0_z_z_y[i] = 2.0 * tr_zz_y[i] * tbe_0 + 2.0 * tr_z_yz[i] * tke_0 - tr_0_y[i];

        tr_0_0_z_z_z[i] = 2.0 * tr_zz_z[i] * tbe_0 + 2.0 * tr_z_zz[i] * tke_0 - tr_0_z[i] - tr_z_0[i];
    }

}

} // t2cgeom namespace

