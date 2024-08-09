#include "GeometricalDerivatives1X1ForPP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_pp(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_pp,
                        const size_t idx_op_ss,
                        const size_t idx_op_sd,
                        const size_t idx_op_ds,
                        const size_t idx_op_dd,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SS

        auto to_0_0 = pbuffer.data(idx_op_ss + i * 1 + 0);

        // Set up components of auxiliary buffer : SD

        auto to_0_xx = pbuffer.data(idx_op_sd + i * 6 + 0);

        auto to_0_xy = pbuffer.data(idx_op_sd + i * 6 + 1);

        auto to_0_xz = pbuffer.data(idx_op_sd + i * 6 + 2);

        auto to_0_yy = pbuffer.data(idx_op_sd + i * 6 + 3);

        auto to_0_yz = pbuffer.data(idx_op_sd + i * 6 + 4);

        auto to_0_zz = pbuffer.data(idx_op_sd + i * 6 + 5);

        // Set up components of auxiliary buffer : DS

        auto to_xx_0 = pbuffer.data(idx_op_ds + i * 6 + 0);

        auto to_xy_0 = pbuffer.data(idx_op_ds + i * 6 + 1);

        auto to_xz_0 = pbuffer.data(idx_op_ds + i * 6 + 2);

        auto to_yy_0 = pbuffer.data(idx_op_ds + i * 6 + 3);

        auto to_yz_0 = pbuffer.data(idx_op_ds + i * 6 + 4);

        auto to_zz_0 = pbuffer.data(idx_op_ds + i * 6 + 5);

        // Set up components of auxiliary buffer : DD

        auto to_xx_xx = pbuffer.data(idx_op_dd + i * 36 + 0);

        auto to_xx_xy = pbuffer.data(idx_op_dd + i * 36 + 1);

        auto to_xx_xz = pbuffer.data(idx_op_dd + i * 36 + 2);

        auto to_xx_yy = pbuffer.data(idx_op_dd + i * 36 + 3);

        auto to_xx_yz = pbuffer.data(idx_op_dd + i * 36 + 4);

        auto to_xx_zz = pbuffer.data(idx_op_dd + i * 36 + 5);

        auto to_xy_xx = pbuffer.data(idx_op_dd + i * 36 + 6);

        auto to_xy_xy = pbuffer.data(idx_op_dd + i * 36 + 7);

        auto to_xy_xz = pbuffer.data(idx_op_dd + i * 36 + 8);

        auto to_xy_yy = pbuffer.data(idx_op_dd + i * 36 + 9);

        auto to_xy_yz = pbuffer.data(idx_op_dd + i * 36 + 10);

        auto to_xy_zz = pbuffer.data(idx_op_dd + i * 36 + 11);

        auto to_xz_xx = pbuffer.data(idx_op_dd + i * 36 + 12);

        auto to_xz_xy = pbuffer.data(idx_op_dd + i * 36 + 13);

        auto to_xz_xz = pbuffer.data(idx_op_dd + i * 36 + 14);

        auto to_xz_yy = pbuffer.data(idx_op_dd + i * 36 + 15);

        auto to_xz_yz = pbuffer.data(idx_op_dd + i * 36 + 16);

        auto to_xz_zz = pbuffer.data(idx_op_dd + i * 36 + 17);

        auto to_yy_xx = pbuffer.data(idx_op_dd + i * 36 + 18);

        auto to_yy_xy = pbuffer.data(idx_op_dd + i * 36 + 19);

        auto to_yy_xz = pbuffer.data(idx_op_dd + i * 36 + 20);

        auto to_yy_yy = pbuffer.data(idx_op_dd + i * 36 + 21);

        auto to_yy_yz = pbuffer.data(idx_op_dd + i * 36 + 22);

        auto to_yy_zz = pbuffer.data(idx_op_dd + i * 36 + 23);

        auto to_yz_xx = pbuffer.data(idx_op_dd + i * 36 + 24);

        auto to_yz_xy = pbuffer.data(idx_op_dd + i * 36 + 25);

        auto to_yz_xz = pbuffer.data(idx_op_dd + i * 36 + 26);

        auto to_yz_yy = pbuffer.data(idx_op_dd + i * 36 + 27);

        auto to_yz_yz = pbuffer.data(idx_op_dd + i * 36 + 28);

        auto to_yz_zz = pbuffer.data(idx_op_dd + i * 36 + 29);

        auto to_zz_xx = pbuffer.data(idx_op_dd + i * 36 + 30);

        auto to_zz_xy = pbuffer.data(idx_op_dd + i * 36 + 31);

        auto to_zz_xz = pbuffer.data(idx_op_dd + i * 36 + 32);

        auto to_zz_yy = pbuffer.data(idx_op_dd + i * 36 + 33);

        auto to_zz_yz = pbuffer.data(idx_op_dd + i * 36 + 34);

        auto to_zz_zz = pbuffer.data(idx_op_dd + i * 36 + 35);

        // Set up 0-3 components of targeted buffer : PP

        auto to_x_x_x_x = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 0);

        auto to_x_x_x_y = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 1);

        auto to_x_x_x_z = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_0_0, to_0_xx, to_0_xy, to_0_xz, to_x_x_x_x, to_x_x_x_y, to_x_x_x_z, to_xx_0, to_xx_xx, to_xx_xy, to_xx_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_x_x[k] = to_0_0[k] - 2.0 * to_0_xx[k] * tke_0 - 2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xx_xx[k] * tbe_0 * tke_0;

            to_x_x_x_y[k] = -2.0 * to_0_xy[k] * tke_0 + 4.0 * to_xx_xy[k] * tbe_0 * tke_0;

            to_x_x_x_z[k] = -2.0 * to_0_xz[k] * tke_0 + 4.0 * to_xx_xz[k] * tbe_0 * tke_0;
        }

        // Set up 3-6 components of targeted buffer : PP

        auto to_x_x_y_x = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 3);

        auto to_x_x_y_y = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 4);

        auto to_x_x_y_z = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_x_x_y_x, to_x_x_y_y, to_x_x_y_z, to_xy_0, to_xy_xx, to_xy_xy, to_xy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_y_x[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xy_xx[k] * tbe_0 * tke_0;

            to_x_x_y_y[k] = 4.0 * to_xy_xy[k] * tbe_0 * tke_0;

            to_x_x_y_z[k] = 4.0 * to_xy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 6-9 components of targeted buffer : PP

        auto to_x_x_z_x = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 6);

        auto to_x_x_z_y = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 7);

        auto to_x_x_z_z = pbuffer.data(idx_op_geom_101_pp + 0 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_x_x_z_x, to_x_x_z_y, to_x_x_z_z, to_xz_0, to_xz_xx, to_xz_xy, to_xz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_z_x[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xz_xx[k] * tbe_0 * tke_0;

            to_x_x_z_y[k] = 4.0 * to_xz_xy[k] * tbe_0 * tke_0;

            to_x_x_z_z[k] = 4.0 * to_xz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 9-12 components of targeted buffer : PP

        auto to_x_y_x_x = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 0);

        auto to_x_y_x_y = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 1);

        auto to_x_y_x_z = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_0_0, to_0_xy, to_0_yy, to_0_yz, to_x_y_x_x, to_x_y_x_y, to_x_y_x_z, to_xx_0, to_xx_xy, to_xx_yy, to_xx_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_x_x[k] = -2.0 * to_0_xy[k] * tke_0 + 4.0 * to_xx_xy[k] * tbe_0 * tke_0;

            to_x_y_x_y[k] = to_0_0[k] - 2.0 * to_0_yy[k] * tke_0 - 2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xx_yy[k] * tbe_0 * tke_0;

            to_x_y_x_z[k] = -2.0 * to_0_yz[k] * tke_0 + 4.0 * to_xx_yz[k] * tbe_0 * tke_0;
        }

        // Set up 12-15 components of targeted buffer : PP

        auto to_x_y_y_x = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 3);

        auto to_x_y_y_y = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 4);

        auto to_x_y_y_z = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_x_y_y_x, to_x_y_y_y, to_x_y_y_z, to_xy_0, to_xy_xy, to_xy_yy, to_xy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_y_x[k] = 4.0 * to_xy_xy[k] * tbe_0 * tke_0;

            to_x_y_y_y[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xy_yy[k] * tbe_0 * tke_0;

            to_x_y_y_z[k] = 4.0 * to_xy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 15-18 components of targeted buffer : PP

        auto to_x_y_z_x = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 6);

        auto to_x_y_z_y = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 7);

        auto to_x_y_z_z = pbuffer.data(idx_op_geom_101_pp + 1 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_x_y_z_x, to_x_y_z_y, to_x_y_z_z, to_xz_0, to_xz_xy, to_xz_yy, to_xz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_z_x[k] = 4.0 * to_xz_xy[k] * tbe_0 * tke_0;

            to_x_y_z_y[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xz_yy[k] * tbe_0 * tke_0;

            to_x_y_z_z[k] = 4.0 * to_xz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 18-21 components of targeted buffer : PP

        auto to_x_z_x_x = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 0);

        auto to_x_z_x_y = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 1);

        auto to_x_z_x_z = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_0_0, to_0_xz, to_0_yz, to_0_zz, to_x_z_x_x, to_x_z_x_y, to_x_z_x_z, to_xx_0, to_xx_xz, to_xx_yz, to_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_x_x[k] = -2.0 * to_0_xz[k] * tke_0 + 4.0 * to_xx_xz[k] * tbe_0 * tke_0;

            to_x_z_x_y[k] = -2.0 * to_0_yz[k] * tke_0 + 4.0 * to_xx_yz[k] * tbe_0 * tke_0;

            to_x_z_x_z[k] = to_0_0[k] - 2.0 * to_0_zz[k] * tke_0 - 2.0 * to_xx_0[k] * tbe_0 + 4.0 * to_xx_zz[k] * tbe_0 * tke_0;
        }

        // Set up 21-24 components of targeted buffer : PP

        auto to_x_z_y_x = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 3);

        auto to_x_z_y_y = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 4);

        auto to_x_z_y_z = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_x_z_y_x, to_x_z_y_y, to_x_z_y_z, to_xy_0, to_xy_xz, to_xy_yz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_y_x[k] = 4.0 * to_xy_xz[k] * tbe_0 * tke_0;

            to_x_z_y_y[k] = 4.0 * to_xy_yz[k] * tbe_0 * tke_0;

            to_x_z_y_z[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 24-27 components of targeted buffer : PP

        auto to_x_z_z_x = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 6);

        auto to_x_z_z_y = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 7);

        auto to_x_z_z_z = pbuffer.data(idx_op_geom_101_pp + 2 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_x_z_z_x, to_x_z_z_y, to_x_z_z_z, to_xz_0, to_xz_xz, to_xz_yz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_z_x[k] = 4.0 * to_xz_xz[k] * tbe_0 * tke_0;

            to_x_z_z_y[k] = 4.0 * to_xz_yz[k] * tbe_0 * tke_0;

            to_x_z_z_z[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 27-30 components of targeted buffer : PP

        auto to_y_x_x_x = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 0);

        auto to_y_x_x_y = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 1);

        auto to_y_x_x_z = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_xy_0, to_xy_xx, to_xy_xy, to_xy_xz, to_y_x_x_x, to_y_x_x_y, to_y_x_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_x_x[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xy_xx[k] * tbe_0 * tke_0;

            to_y_x_x_y[k] = 4.0 * to_xy_xy[k] * tbe_0 * tke_0;

            to_y_x_x_z[k] = 4.0 * to_xy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 30-33 components of targeted buffer : PP

        auto to_y_x_y_x = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 3);

        auto to_y_x_y_y = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 4);

        auto to_y_x_y_z = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_0_0, to_0_xx, to_0_xy, to_0_xz, to_y_x_y_x, to_y_x_y_y, to_y_x_y_z, to_yy_0, to_yy_xx, to_yy_xy, to_yy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_y_x[k] = to_0_0[k] - 2.0 * to_0_xx[k] * tke_0 - 2.0 * to_yy_0[k] * tbe_0 + 4.0 * to_yy_xx[k] * tbe_0 * tke_0;

            to_y_x_y_y[k] = -2.0 * to_0_xy[k] * tke_0 + 4.0 * to_yy_xy[k] * tbe_0 * tke_0;

            to_y_x_y_z[k] = -2.0 * to_0_xz[k] * tke_0 + 4.0 * to_yy_xz[k] * tbe_0 * tke_0;
        }

        // Set up 33-36 components of targeted buffer : PP

        auto to_y_x_z_x = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 6);

        auto to_y_x_z_y = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 7);

        auto to_y_x_z_z = pbuffer.data(idx_op_geom_101_pp + 3 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_y_x_z_x, to_y_x_z_y, to_y_x_z_z, to_yz_0, to_yz_xx, to_yz_xy, to_yz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_z_x[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yz_xx[k] * tbe_0 * tke_0;

            to_y_x_z_y[k] = 4.0 * to_yz_xy[k] * tbe_0 * tke_0;

            to_y_x_z_z[k] = 4.0 * to_yz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 36-39 components of targeted buffer : PP

        auto to_y_y_x_x = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 0);

        auto to_y_y_x_y = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 1);

        auto to_y_y_x_z = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_xy_0, to_xy_xy, to_xy_yy, to_xy_yz, to_y_y_x_x, to_y_y_x_y, to_y_y_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_x_x[k] = 4.0 * to_xy_xy[k] * tbe_0 * tke_0;

            to_y_y_x_y[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xy_yy[k] * tbe_0 * tke_0;

            to_y_y_x_z[k] = 4.0 * to_xy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 39-42 components of targeted buffer : PP

        auto to_y_y_y_x = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 3);

        auto to_y_y_y_y = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 4);

        auto to_y_y_y_z = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_0_0, to_0_xy, to_0_yy, to_0_yz, to_y_y_y_x, to_y_y_y_y, to_y_y_y_z, to_yy_0, to_yy_xy, to_yy_yy, to_yy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_y_x[k] = -2.0 * to_0_xy[k] * tke_0 + 4.0 * to_yy_xy[k] * tbe_0 * tke_0;

            to_y_y_y_y[k] = to_0_0[k] - 2.0 * to_0_yy[k] * tke_0 - 2.0 * to_yy_0[k] * tbe_0 + 4.0 * to_yy_yy[k] * tbe_0 * tke_0;

            to_y_y_y_z[k] = -2.0 * to_0_yz[k] * tke_0 + 4.0 * to_yy_yz[k] * tbe_0 * tke_0;
        }

        // Set up 42-45 components of targeted buffer : PP

        auto to_y_y_z_x = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 6);

        auto to_y_y_z_y = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 7);

        auto to_y_y_z_z = pbuffer.data(idx_op_geom_101_pp + 4 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_y_y_z_x, to_y_y_z_y, to_y_y_z_z, to_yz_0, to_yz_xy, to_yz_yy, to_yz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_z_x[k] = 4.0 * to_yz_xy[k] * tbe_0 * tke_0;

            to_y_y_z_y[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yz_yy[k] * tbe_0 * tke_0;

            to_y_y_z_z[k] = 4.0 * to_yz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 45-48 components of targeted buffer : PP

        auto to_y_z_x_x = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 0);

        auto to_y_z_x_y = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 1);

        auto to_y_z_x_z = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_xy_0, to_xy_xz, to_xy_yz, to_xy_zz, to_y_z_x_x, to_y_z_x_y, to_y_z_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_x_x[k] = 4.0 * to_xy_xz[k] * tbe_0 * tke_0;

            to_y_z_x_y[k] = 4.0 * to_xy_yz[k] * tbe_0 * tke_0;

            to_y_z_x_z[k] = -2.0 * to_xy_0[k] * tbe_0 + 4.0 * to_xy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 48-51 components of targeted buffer : PP

        auto to_y_z_y_x = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 3);

        auto to_y_z_y_y = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 4);

        auto to_y_z_y_z = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_0_0, to_0_xz, to_0_yz, to_0_zz, to_y_z_y_x, to_y_z_y_y, to_y_z_y_z, to_yy_0, to_yy_xz, to_yy_yz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_y_x[k] = -2.0 * to_0_xz[k] * tke_0 + 4.0 * to_yy_xz[k] * tbe_0 * tke_0;

            to_y_z_y_y[k] = -2.0 * to_0_yz[k] * tke_0 + 4.0 * to_yy_yz[k] * tbe_0 * tke_0;

            to_y_z_y_z[k] = to_0_0[k] - 2.0 * to_0_zz[k] * tke_0 - 2.0 * to_yy_0[k] * tbe_0 + 4.0 * to_yy_zz[k] * tbe_0 * tke_0;
        }

        // Set up 51-54 components of targeted buffer : PP

        auto to_y_z_z_x = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 6);

        auto to_y_z_z_y = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 7);

        auto to_y_z_z_z = pbuffer.data(idx_op_geom_101_pp + 5 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_y_z_z_x, to_y_z_z_y, to_y_z_z_z, to_yz_0, to_yz_xz, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_z_x[k] = 4.0 * to_yz_xz[k] * tbe_0 * tke_0;

            to_y_z_z_y[k] = 4.0 * to_yz_yz[k] * tbe_0 * tke_0;

            to_y_z_z_z[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 54-57 components of targeted buffer : PP

        auto to_z_x_x_x = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 0);

        auto to_z_x_x_y = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 1);

        auto to_z_x_x_z = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_xz_0, to_xz_xx, to_xz_xy, to_xz_xz, to_z_x_x_x, to_z_x_x_y, to_z_x_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_x_x[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xz_xx[k] * tbe_0 * tke_0;

            to_z_x_x_y[k] = 4.0 * to_xz_xy[k] * tbe_0 * tke_0;

            to_z_x_x_z[k] = 4.0 * to_xz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 57-60 components of targeted buffer : PP

        auto to_z_x_y_x = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 3);

        auto to_z_x_y_y = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 4);

        auto to_z_x_y_z = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_yz_0, to_yz_xx, to_yz_xy, to_yz_xz, to_z_x_y_x, to_z_x_y_y, to_z_x_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_y_x[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yz_xx[k] * tbe_0 * tke_0;

            to_z_x_y_y[k] = 4.0 * to_yz_xy[k] * tbe_0 * tke_0;

            to_z_x_y_z[k] = 4.0 * to_yz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 60-63 components of targeted buffer : PP

        auto to_z_x_z_x = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 6);

        auto to_z_x_z_y = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 7);

        auto to_z_x_z_z = pbuffer.data(idx_op_geom_101_pp + 6 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_0_0, to_0_xx, to_0_xy, to_0_xz, to_z_x_z_x, to_z_x_z_y, to_z_x_z_z, to_zz_0, to_zz_xx, to_zz_xy, to_zz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_z_x[k] = to_0_0[k] - 2.0 * to_0_xx[k] * tke_0 - 2.0 * to_zz_0[k] * tbe_0 + 4.0 * to_zz_xx[k] * tbe_0 * tke_0;

            to_z_x_z_y[k] = -2.0 * to_0_xy[k] * tke_0 + 4.0 * to_zz_xy[k] * tbe_0 * tke_0;

            to_z_x_z_z[k] = -2.0 * to_0_xz[k] * tke_0 + 4.0 * to_zz_xz[k] * tbe_0 * tke_0;
        }

        // Set up 63-66 components of targeted buffer : PP

        auto to_z_y_x_x = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 0);

        auto to_z_y_x_y = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 1);

        auto to_z_y_x_z = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_xz_0, to_xz_xy, to_xz_yy, to_xz_yz, to_z_y_x_x, to_z_y_x_y, to_z_y_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_x_x[k] = 4.0 * to_xz_xy[k] * tbe_0 * tke_0;

            to_z_y_x_y[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xz_yy[k] * tbe_0 * tke_0;

            to_z_y_x_z[k] = 4.0 * to_xz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 66-69 components of targeted buffer : PP

        auto to_z_y_y_x = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 3);

        auto to_z_y_y_y = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 4);

        auto to_z_y_y_z = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_yz_0, to_yz_xy, to_yz_yy, to_yz_yz, to_z_y_y_x, to_z_y_y_y, to_z_y_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_y_x[k] = 4.0 * to_yz_xy[k] * tbe_0 * tke_0;

            to_z_y_y_y[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yz_yy[k] * tbe_0 * tke_0;

            to_z_y_y_z[k] = 4.0 * to_yz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 69-72 components of targeted buffer : PP

        auto to_z_y_z_x = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 6);

        auto to_z_y_z_y = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 7);

        auto to_z_y_z_z = pbuffer.data(idx_op_geom_101_pp + 7 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_0_0, to_0_xy, to_0_yy, to_0_yz, to_z_y_z_x, to_z_y_z_y, to_z_y_z_z, to_zz_0, to_zz_xy, to_zz_yy, to_zz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_z_x[k] = -2.0 * to_0_xy[k] * tke_0 + 4.0 * to_zz_xy[k] * tbe_0 * tke_0;

            to_z_y_z_y[k] = to_0_0[k] - 2.0 * to_0_yy[k] * tke_0 - 2.0 * to_zz_0[k] * tbe_0 + 4.0 * to_zz_yy[k] * tbe_0 * tke_0;

            to_z_y_z_z[k] = -2.0 * to_0_yz[k] * tke_0 + 4.0 * to_zz_yz[k] * tbe_0 * tke_0;
        }

        // Set up 72-75 components of targeted buffer : PP

        auto to_z_z_x_x = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 0);

        auto to_z_z_x_y = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 1);

        auto to_z_z_x_z = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_xz_0, to_xz_xz, to_xz_yz, to_xz_zz, to_z_z_x_x, to_z_z_x_y, to_z_z_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_x_x[k] = 4.0 * to_xz_xz[k] * tbe_0 * tke_0;

            to_z_z_x_y[k] = 4.0 * to_xz_yz[k] * tbe_0 * tke_0;

            to_z_z_x_z[k] = -2.0 * to_xz_0[k] * tbe_0 + 4.0 * to_xz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 75-78 components of targeted buffer : PP

        auto to_z_z_y_x = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 3);

        auto to_z_z_y_y = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 4);

        auto to_z_z_y_z = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_yz_0, to_yz_xz, to_yz_yz, to_yz_zz, to_z_z_y_x, to_z_z_y_y, to_z_z_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_y_x[k] = 4.0 * to_yz_xz[k] * tbe_0 * tke_0;

            to_z_z_y_y[k] = 4.0 * to_yz_yz[k] * tbe_0 * tke_0;

            to_z_z_y_z[k] = -2.0 * to_yz_0[k] * tbe_0 + 4.0 * to_yz_zz[k] * tbe_0 * tke_0;
        }

        // Set up 78-81 components of targeted buffer : PP

        auto to_z_z_z_x = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 6);

        auto to_z_z_z_y = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 7);

        auto to_z_z_z_z = pbuffer.data(idx_op_geom_101_pp + 8 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_0_0, to_0_xz, to_0_yz, to_0_zz, to_z_z_z_x, to_z_z_z_y, to_z_z_z_z, to_zz_0, to_zz_xz, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_z_x[k] = -2.0 * to_0_xz[k] * tke_0 + 4.0 * to_zz_xz[k] * tbe_0 * tke_0;

            to_z_z_z_y[k] = -2.0 * to_0_yz[k] * tke_0 + 4.0 * to_zz_yz[k] * tbe_0 * tke_0;

            to_z_z_z_z[k] = to_0_0[k] - 2.0 * to_0_zz[k] * tke_0 - 2.0 * to_zz_0[k] * tbe_0 + 4.0 * to_zz_zz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

