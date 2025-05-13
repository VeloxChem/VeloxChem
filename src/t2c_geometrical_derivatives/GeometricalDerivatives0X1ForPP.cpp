#include "GeometricalDerivatives0X1ForPP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_pp(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_pp,
                       const int idx_op_ps,
                       const int idx_op_pd,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PS

        auto to_x_0 = pbuffer.data(idx_op_ps + i * 3 + 0);

        auto to_y_0 = pbuffer.data(idx_op_ps + i * 3 + 1);

        auto to_z_0 = pbuffer.data(idx_op_ps + i * 3 + 2);

        // Set up components of auxiliary buffer : PD

        auto to_x_xx = pbuffer.data(idx_op_pd + i * 18 + 0);

        auto to_x_xy = pbuffer.data(idx_op_pd + i * 18 + 1);

        auto to_x_xz = pbuffer.data(idx_op_pd + i * 18 + 2);

        auto to_x_yy = pbuffer.data(idx_op_pd + i * 18 + 3);

        auto to_x_yz = pbuffer.data(idx_op_pd + i * 18 + 4);

        auto to_x_zz = pbuffer.data(idx_op_pd + i * 18 + 5);

        auto to_y_xx = pbuffer.data(idx_op_pd + i * 18 + 6);

        auto to_y_xy = pbuffer.data(idx_op_pd + i * 18 + 7);

        auto to_y_xz = pbuffer.data(idx_op_pd + i * 18 + 8);

        auto to_y_yy = pbuffer.data(idx_op_pd + i * 18 + 9);

        auto to_y_yz = pbuffer.data(idx_op_pd + i * 18 + 10);

        auto to_y_zz = pbuffer.data(idx_op_pd + i * 18 + 11);

        auto to_z_xx = pbuffer.data(idx_op_pd + i * 18 + 12);

        auto to_z_xy = pbuffer.data(idx_op_pd + i * 18 + 13);

        auto to_z_xz = pbuffer.data(idx_op_pd + i * 18 + 14);

        auto to_z_yy = pbuffer.data(idx_op_pd + i * 18 + 15);

        auto to_z_yz = pbuffer.data(idx_op_pd + i * 18 + 16);

        auto to_z_zz = pbuffer.data(idx_op_pd + i * 18 + 17);

        // Set up 0-3 components of targeted buffer : PP

        auto to_0_x_x_x = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 0);

        auto to_0_x_x_y = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 1);

        auto to_0_x_x_z = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_0_x_x_x, to_0_x_x_y, to_0_x_x_z, to_x_0, to_x_xx, to_x_xy, to_x_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_x_x[k] = -to_x_0[k] + 2.0 * to_x_xx[k] * tke_0;

            to_0_x_x_y[k] = 2.0 * to_x_xy[k] * tke_0;

            to_0_x_x_z[k] = 2.0 * to_x_xz[k] * tke_0;
        }

        // Set up 3-6 components of targeted buffer : PP

        auto to_0_x_y_x = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 3);

        auto to_0_x_y_y = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 4);

        auto to_0_x_y_z = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_0_x_y_x, to_0_x_y_y, to_0_x_y_z, to_y_0, to_y_xx, to_y_xy, to_y_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_y_x[k] = -to_y_0[k] + 2.0 * to_y_xx[k] * tke_0;

            to_0_x_y_y[k] = 2.0 * to_y_xy[k] * tke_0;

            to_0_x_y_z[k] = 2.0 * to_y_xz[k] * tke_0;
        }

        // Set up 6-9 components of targeted buffer : PP

        auto to_0_x_z_x = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 6);

        auto to_0_x_z_y = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 7);

        auto to_0_x_z_z = pbuffer.data(idx_op_geom_001_pp + 0 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_0_x_z_x, to_0_x_z_y, to_0_x_z_z, to_z_0, to_z_xx, to_z_xy, to_z_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_z_x[k] = -to_z_0[k] + 2.0 * to_z_xx[k] * tke_0;

            to_0_x_z_y[k] = 2.0 * to_z_xy[k] * tke_0;

            to_0_x_z_z[k] = 2.0 * to_z_xz[k] * tke_0;
        }

        // Set up 9-12 components of targeted buffer : PP

        auto to_0_y_x_x = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 0);

        auto to_0_y_x_y = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 1);

        auto to_0_y_x_z = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_0_y_x_x, to_0_y_x_y, to_0_y_x_z, to_x_0, to_x_xy, to_x_yy, to_x_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_x_x[k] = 2.0 * to_x_xy[k] * tke_0;

            to_0_y_x_y[k] = -to_x_0[k] + 2.0 * to_x_yy[k] * tke_0;

            to_0_y_x_z[k] = 2.0 * to_x_yz[k] * tke_0;
        }

        // Set up 12-15 components of targeted buffer : PP

        auto to_0_y_y_x = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 3);

        auto to_0_y_y_y = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 4);

        auto to_0_y_y_z = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_0_y_y_x, to_0_y_y_y, to_0_y_y_z, to_y_0, to_y_xy, to_y_yy, to_y_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_y_x[k] = 2.0 * to_y_xy[k] * tke_0;

            to_0_y_y_y[k] = -to_y_0[k] + 2.0 * to_y_yy[k] * tke_0;

            to_0_y_y_z[k] = 2.0 * to_y_yz[k] * tke_0;
        }

        // Set up 15-18 components of targeted buffer : PP

        auto to_0_y_z_x = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 6);

        auto to_0_y_z_y = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 7);

        auto to_0_y_z_z = pbuffer.data(idx_op_geom_001_pp + 1 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_0_y_z_x, to_0_y_z_y, to_0_y_z_z, to_z_0, to_z_xy, to_z_yy, to_z_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_z_x[k] = 2.0 * to_z_xy[k] * tke_0;

            to_0_y_z_y[k] = -to_z_0[k] + 2.0 * to_z_yy[k] * tke_0;

            to_0_y_z_z[k] = 2.0 * to_z_yz[k] * tke_0;
        }

        // Set up 18-21 components of targeted buffer : PP

        auto to_0_z_x_x = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 0);

        auto to_0_z_x_y = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 1);

        auto to_0_z_x_z = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 2);

        #pragma omp simd aligned(to_0_z_x_x, to_0_z_x_y, to_0_z_x_z, to_x_0, to_x_xz, to_x_yz, to_x_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_x_x[k] = 2.0 * to_x_xz[k] * tke_0;

            to_0_z_x_y[k] = 2.0 * to_x_yz[k] * tke_0;

            to_0_z_x_z[k] = -to_x_0[k] + 2.0 * to_x_zz[k] * tke_0;
        }

        // Set up 21-24 components of targeted buffer : PP

        auto to_0_z_y_x = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 3);

        auto to_0_z_y_y = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 4);

        auto to_0_z_y_z = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 5);

        #pragma omp simd aligned(to_0_z_y_x, to_0_z_y_y, to_0_z_y_z, to_y_0, to_y_xz, to_y_yz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_y_x[k] = 2.0 * to_y_xz[k] * tke_0;

            to_0_z_y_y[k] = 2.0 * to_y_yz[k] * tke_0;

            to_0_z_y_z[k] = -to_y_0[k] + 2.0 * to_y_zz[k] * tke_0;
        }

        // Set up 24-27 components of targeted buffer : PP

        auto to_0_z_z_x = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 6);

        auto to_0_z_z_y = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 7);

        auto to_0_z_z_z = pbuffer.data(idx_op_geom_001_pp + 2 * op_comps * 9 + i * 9 + 8);

        #pragma omp simd aligned(to_0_z_z_x, to_0_z_z_y, to_0_z_z_z, to_z_0, to_z_xz, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_z_x[k] = 2.0 * to_z_xz[k] * tke_0;

            to_0_z_z_y[k] = 2.0 * to_z_yz[k] * tke_0;

            to_0_z_z_z[k] = -to_z_0[k] + 2.0 * to_z_zz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

