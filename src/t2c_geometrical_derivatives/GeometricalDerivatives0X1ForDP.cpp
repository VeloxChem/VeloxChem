#include "GeometricalDerivatives0X1ForDP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_dp(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_dp,
                       const int idx_op_ds,
                       const int idx_op_dd,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up 0-3 components of targeted buffer : DP

        auto to_0_x_xx_x = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 0);

        auto to_0_x_xx_y = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 1);

        auto to_0_x_xx_z = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_0_x_xx_x, to_0_x_xx_y, to_0_x_xx_z, to_xx_0, to_xx_xx, to_xx_xy, to_xx_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xx_x[k] = -to_xx_0[k] + 2.0 * to_xx_xx[k] * tke_0;

            to_0_x_xx_y[k] = 2.0 * to_xx_xy[k] * tke_0;

            to_0_x_xx_z[k] = 2.0 * to_xx_xz[k] * tke_0;
        }

        // Set up 3-6 components of targeted buffer : DP

        auto to_0_x_xy_x = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 3);

        auto to_0_x_xy_y = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 4);

        auto to_0_x_xy_z = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_0_x_xy_x, to_0_x_xy_y, to_0_x_xy_z, to_xy_0, to_xy_xx, to_xy_xy, to_xy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xy_x[k] = -to_xy_0[k] + 2.0 * to_xy_xx[k] * tke_0;

            to_0_x_xy_y[k] = 2.0 * to_xy_xy[k] * tke_0;

            to_0_x_xy_z[k] = 2.0 * to_xy_xz[k] * tke_0;
        }

        // Set up 6-9 components of targeted buffer : DP

        auto to_0_x_xz_x = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 6);

        auto to_0_x_xz_y = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 7);

        auto to_0_x_xz_z = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_0_x_xz_x, to_0_x_xz_y, to_0_x_xz_z, to_xz_0, to_xz_xx, to_xz_xy, to_xz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_xz_x[k] = -to_xz_0[k] + 2.0 * to_xz_xx[k] * tke_0;

            to_0_x_xz_y[k] = 2.0 * to_xz_xy[k] * tke_0;

            to_0_x_xz_z[k] = 2.0 * to_xz_xz[k] * tke_0;
        }

        // Set up 9-12 components of targeted buffer : DP

        auto to_0_x_yy_x = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 9);

        auto to_0_x_yy_y = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 10);

        auto to_0_x_yy_z = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_0_x_yy_x, to_0_x_yy_y, to_0_x_yy_z, to_yy_0, to_yy_xx, to_yy_xy, to_yy_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yy_x[k] = -to_yy_0[k] + 2.0 * to_yy_xx[k] * tke_0;

            to_0_x_yy_y[k] = 2.0 * to_yy_xy[k] * tke_0;

            to_0_x_yy_z[k] = 2.0 * to_yy_xz[k] * tke_0;
        }

        // Set up 12-15 components of targeted buffer : DP

        auto to_0_x_yz_x = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 12);

        auto to_0_x_yz_y = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 13);

        auto to_0_x_yz_z = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_0_x_yz_x, to_0_x_yz_y, to_0_x_yz_z, to_yz_0, to_yz_xx, to_yz_xy, to_yz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_yz_x[k] = -to_yz_0[k] + 2.0 * to_yz_xx[k] * tke_0;

            to_0_x_yz_y[k] = 2.0 * to_yz_xy[k] * tke_0;

            to_0_x_yz_z[k] = 2.0 * to_yz_xz[k] * tke_0;
        }

        // Set up 15-18 components of targeted buffer : DP

        auto to_0_x_zz_x = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 15);

        auto to_0_x_zz_y = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 16);

        auto to_0_x_zz_z = pbuffer.data(idx_op_geom_001_dp + 0 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_0_x_zz_x, to_0_x_zz_y, to_0_x_zz_z, to_zz_0, to_zz_xx, to_zz_xy, to_zz_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_zz_x[k] = -to_zz_0[k] + 2.0 * to_zz_xx[k] * tke_0;

            to_0_x_zz_y[k] = 2.0 * to_zz_xy[k] * tke_0;

            to_0_x_zz_z[k] = 2.0 * to_zz_xz[k] * tke_0;
        }

        // Set up 18-21 components of targeted buffer : DP

        auto to_0_y_xx_x = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 0);

        auto to_0_y_xx_y = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 1);

        auto to_0_y_xx_z = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_0_y_xx_x, to_0_y_xx_y, to_0_y_xx_z, to_xx_0, to_xx_xy, to_xx_yy, to_xx_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xx_x[k] = 2.0 * to_xx_xy[k] * tke_0;

            to_0_y_xx_y[k] = -to_xx_0[k] + 2.0 * to_xx_yy[k] * tke_0;

            to_0_y_xx_z[k] = 2.0 * to_xx_yz[k] * tke_0;
        }

        // Set up 21-24 components of targeted buffer : DP

        auto to_0_y_xy_x = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 3);

        auto to_0_y_xy_y = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 4);

        auto to_0_y_xy_z = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_0_y_xy_x, to_0_y_xy_y, to_0_y_xy_z, to_xy_0, to_xy_xy, to_xy_yy, to_xy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xy_x[k] = 2.0 * to_xy_xy[k] * tke_0;

            to_0_y_xy_y[k] = -to_xy_0[k] + 2.0 * to_xy_yy[k] * tke_0;

            to_0_y_xy_z[k] = 2.0 * to_xy_yz[k] * tke_0;
        }

        // Set up 24-27 components of targeted buffer : DP

        auto to_0_y_xz_x = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 6);

        auto to_0_y_xz_y = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 7);

        auto to_0_y_xz_z = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_0_y_xz_x, to_0_y_xz_y, to_0_y_xz_z, to_xz_0, to_xz_xy, to_xz_yy, to_xz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_xz_x[k] = 2.0 * to_xz_xy[k] * tke_0;

            to_0_y_xz_y[k] = -to_xz_0[k] + 2.0 * to_xz_yy[k] * tke_0;

            to_0_y_xz_z[k] = 2.0 * to_xz_yz[k] * tke_0;
        }

        // Set up 27-30 components of targeted buffer : DP

        auto to_0_y_yy_x = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 9);

        auto to_0_y_yy_y = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 10);

        auto to_0_y_yy_z = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_0_y_yy_x, to_0_y_yy_y, to_0_y_yy_z, to_yy_0, to_yy_xy, to_yy_yy, to_yy_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yy_x[k] = 2.0 * to_yy_xy[k] * tke_0;

            to_0_y_yy_y[k] = -to_yy_0[k] + 2.0 * to_yy_yy[k] * tke_0;

            to_0_y_yy_z[k] = 2.0 * to_yy_yz[k] * tke_0;
        }

        // Set up 30-33 components of targeted buffer : DP

        auto to_0_y_yz_x = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 12);

        auto to_0_y_yz_y = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 13);

        auto to_0_y_yz_z = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_0_y_yz_x, to_0_y_yz_y, to_0_y_yz_z, to_yz_0, to_yz_xy, to_yz_yy, to_yz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_yz_x[k] = 2.0 * to_yz_xy[k] * tke_0;

            to_0_y_yz_y[k] = -to_yz_0[k] + 2.0 * to_yz_yy[k] * tke_0;

            to_0_y_yz_z[k] = 2.0 * to_yz_yz[k] * tke_0;
        }

        // Set up 33-36 components of targeted buffer : DP

        auto to_0_y_zz_x = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 15);

        auto to_0_y_zz_y = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 16);

        auto to_0_y_zz_z = pbuffer.data(idx_op_geom_001_dp + 1 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_0_y_zz_x, to_0_y_zz_y, to_0_y_zz_z, to_zz_0, to_zz_xy, to_zz_yy, to_zz_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_zz_x[k] = 2.0 * to_zz_xy[k] * tke_0;

            to_0_y_zz_y[k] = -to_zz_0[k] + 2.0 * to_zz_yy[k] * tke_0;

            to_0_y_zz_z[k] = 2.0 * to_zz_yz[k] * tke_0;
        }

        // Set up 36-39 components of targeted buffer : DP

        auto to_0_z_xx_x = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 0);

        auto to_0_z_xx_y = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 1);

        auto to_0_z_xx_z = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 2);

        #pragma omp simd aligned(to_0_z_xx_x, to_0_z_xx_y, to_0_z_xx_z, to_xx_0, to_xx_xz, to_xx_yz, to_xx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xx_x[k] = 2.0 * to_xx_xz[k] * tke_0;

            to_0_z_xx_y[k] = 2.0 * to_xx_yz[k] * tke_0;

            to_0_z_xx_z[k] = -to_xx_0[k] + 2.0 * to_xx_zz[k] * tke_0;
        }

        // Set up 39-42 components of targeted buffer : DP

        auto to_0_z_xy_x = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 3);

        auto to_0_z_xy_y = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 4);

        auto to_0_z_xy_z = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_0_z_xy_x, to_0_z_xy_y, to_0_z_xy_z, to_xy_0, to_xy_xz, to_xy_yz, to_xy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xy_x[k] = 2.0 * to_xy_xz[k] * tke_0;

            to_0_z_xy_y[k] = 2.0 * to_xy_yz[k] * tke_0;

            to_0_z_xy_z[k] = -to_xy_0[k] + 2.0 * to_xy_zz[k] * tke_0;
        }

        // Set up 42-45 components of targeted buffer : DP

        auto to_0_z_xz_x = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 6);

        auto to_0_z_xz_y = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 7);

        auto to_0_z_xz_z = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 8);

        #pragma omp simd aligned(to_0_z_xz_x, to_0_z_xz_y, to_0_z_xz_z, to_xz_0, to_xz_xz, to_xz_yz, to_xz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_xz_x[k] = 2.0 * to_xz_xz[k] * tke_0;

            to_0_z_xz_y[k] = 2.0 * to_xz_yz[k] * tke_0;

            to_0_z_xz_z[k] = -to_xz_0[k] + 2.0 * to_xz_zz[k] * tke_0;
        }

        // Set up 45-48 components of targeted buffer : DP

        auto to_0_z_yy_x = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 9);

        auto to_0_z_yy_y = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 10);

        auto to_0_z_yy_z = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_0_z_yy_x, to_0_z_yy_y, to_0_z_yy_z, to_yy_0, to_yy_xz, to_yy_yz, to_yy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yy_x[k] = 2.0 * to_yy_xz[k] * tke_0;

            to_0_z_yy_y[k] = 2.0 * to_yy_yz[k] * tke_0;

            to_0_z_yy_z[k] = -to_yy_0[k] + 2.0 * to_yy_zz[k] * tke_0;
        }

        // Set up 48-51 components of targeted buffer : DP

        auto to_0_z_yz_x = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 12);

        auto to_0_z_yz_y = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 13);

        auto to_0_z_yz_z = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 14);

        #pragma omp simd aligned(to_0_z_yz_x, to_0_z_yz_y, to_0_z_yz_z, to_yz_0, to_yz_xz, to_yz_yz, to_yz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_yz_x[k] = 2.0 * to_yz_xz[k] * tke_0;

            to_0_z_yz_y[k] = 2.0 * to_yz_yz[k] * tke_0;

            to_0_z_yz_z[k] = -to_yz_0[k] + 2.0 * to_yz_zz[k] * tke_0;
        }

        // Set up 51-54 components of targeted buffer : DP

        auto to_0_z_zz_x = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 15);

        auto to_0_z_zz_y = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 16);

        auto to_0_z_zz_z = pbuffer.data(idx_op_geom_001_dp + 2 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_0_z_zz_x, to_0_z_zz_y, to_0_z_zz_z, to_zz_0, to_zz_xz, to_zz_yz, to_zz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_zz_x[k] = 2.0 * to_zz_xz[k] * tke_0;

            to_0_z_zz_y[k] = 2.0 * to_zz_yz[k] * tke_0;

            to_0_z_zz_z[k] = -to_zz_0[k] + 2.0 * to_zz_zz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

