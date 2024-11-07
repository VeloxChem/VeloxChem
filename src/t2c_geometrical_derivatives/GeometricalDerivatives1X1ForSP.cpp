#include "GeometricalDerivatives1X1ForSP.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_sp(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_sp,
                        const size_t              idx_op_ps,
                        const size_t              idx_op_pd,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
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

        // Set up 0-3 components of targeted buffer : SP

        auto to_x_x_0_x = pbuffer.data(idx_op_geom_101_sp + 0 * op_comps * 3 + i * 3 + 0);

        auto to_x_x_0_y = pbuffer.data(idx_op_geom_101_sp + 0 * op_comps * 3 + i * 3 + 1);

        auto to_x_x_0_z = pbuffer.data(idx_op_geom_101_sp + 0 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_x_0, to_x_x_0_x, to_x_x_0_y, to_x_x_0_z, to_x_xx, to_x_xy, to_x_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_0_x[k] = -2.0 * to_x_0[k] * tbe_0 + 4.0 * to_x_xx[k] * tbe_0 * tke_0;

            to_x_x_0_y[k] = 4.0 * to_x_xy[k] * tbe_0 * tke_0;

            to_x_x_0_z[k] = 4.0 * to_x_xz[k] * tbe_0 * tke_0;
        }

        // Set up 3-6 components of targeted buffer : SP

        auto to_x_y_0_x = pbuffer.data(idx_op_geom_101_sp + 1 * op_comps * 3 + i * 3 + 0);

        auto to_x_y_0_y = pbuffer.data(idx_op_geom_101_sp + 1 * op_comps * 3 + i * 3 + 1);

        auto to_x_y_0_z = pbuffer.data(idx_op_geom_101_sp + 1 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_x_0, to_x_xy, to_x_y_0_x, to_x_y_0_y, to_x_y_0_z, to_x_yy, to_x_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_0_x[k] = 4.0 * to_x_xy[k] * tbe_0 * tke_0;

            to_x_y_0_y[k] = -2.0 * to_x_0[k] * tbe_0 + 4.0 * to_x_yy[k] * tbe_0 * tke_0;

            to_x_y_0_z[k] = 4.0 * to_x_yz[k] * tbe_0 * tke_0;
        }

        // Set up 6-9 components of targeted buffer : SP

        auto to_x_z_0_x = pbuffer.data(idx_op_geom_101_sp + 2 * op_comps * 3 + i * 3 + 0);

        auto to_x_z_0_y = pbuffer.data(idx_op_geom_101_sp + 2 * op_comps * 3 + i * 3 + 1);

        auto to_x_z_0_z = pbuffer.data(idx_op_geom_101_sp + 2 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_x_0, to_x_xz, to_x_yz, to_x_z_0_x, to_x_z_0_y, to_x_z_0_z, to_x_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_0_x[k] = 4.0 * to_x_xz[k] * tbe_0 * tke_0;

            to_x_z_0_y[k] = 4.0 * to_x_yz[k] * tbe_0 * tke_0;

            to_x_z_0_z[k] = -2.0 * to_x_0[k] * tbe_0 + 4.0 * to_x_zz[k] * tbe_0 * tke_0;
        }

        // Set up 9-12 components of targeted buffer : SP

        auto to_y_x_0_x = pbuffer.data(idx_op_geom_101_sp + 3 * op_comps * 3 + i * 3 + 0);

        auto to_y_x_0_y = pbuffer.data(idx_op_geom_101_sp + 3 * op_comps * 3 + i * 3 + 1);

        auto to_y_x_0_z = pbuffer.data(idx_op_geom_101_sp + 3 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_y_0, to_y_x_0_x, to_y_x_0_y, to_y_x_0_z, to_y_xx, to_y_xy, to_y_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_0_x[k] = -2.0 * to_y_0[k] * tbe_0 + 4.0 * to_y_xx[k] * tbe_0 * tke_0;

            to_y_x_0_y[k] = 4.0 * to_y_xy[k] * tbe_0 * tke_0;

            to_y_x_0_z[k] = 4.0 * to_y_xz[k] * tbe_0 * tke_0;
        }

        // Set up 12-15 components of targeted buffer : SP

        auto to_y_y_0_x = pbuffer.data(idx_op_geom_101_sp + 4 * op_comps * 3 + i * 3 + 0);

        auto to_y_y_0_y = pbuffer.data(idx_op_geom_101_sp + 4 * op_comps * 3 + i * 3 + 1);

        auto to_y_y_0_z = pbuffer.data(idx_op_geom_101_sp + 4 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_y_0, to_y_xy, to_y_y_0_x, to_y_y_0_y, to_y_y_0_z, to_y_yy, to_y_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_0_x[k] = 4.0 * to_y_xy[k] * tbe_0 * tke_0;

            to_y_y_0_y[k] = -2.0 * to_y_0[k] * tbe_0 + 4.0 * to_y_yy[k] * tbe_0 * tke_0;

            to_y_y_0_z[k] = 4.0 * to_y_yz[k] * tbe_0 * tke_0;
        }

        // Set up 15-18 components of targeted buffer : SP

        auto to_y_z_0_x = pbuffer.data(idx_op_geom_101_sp + 5 * op_comps * 3 + i * 3 + 0);

        auto to_y_z_0_y = pbuffer.data(idx_op_geom_101_sp + 5 * op_comps * 3 + i * 3 + 1);

        auto to_y_z_0_z = pbuffer.data(idx_op_geom_101_sp + 5 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_y_0, to_y_xz, to_y_yz, to_y_z_0_x, to_y_z_0_y, to_y_z_0_z, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_0_x[k] = 4.0 * to_y_xz[k] * tbe_0 * tke_0;

            to_y_z_0_y[k] = 4.0 * to_y_yz[k] * tbe_0 * tke_0;

            to_y_z_0_z[k] = -2.0 * to_y_0[k] * tbe_0 + 4.0 * to_y_zz[k] * tbe_0 * tke_0;
        }

        // Set up 18-21 components of targeted buffer : SP

        auto to_z_x_0_x = pbuffer.data(idx_op_geom_101_sp + 6 * op_comps * 3 + i * 3 + 0);

        auto to_z_x_0_y = pbuffer.data(idx_op_geom_101_sp + 6 * op_comps * 3 + i * 3 + 1);

        auto to_z_x_0_z = pbuffer.data(idx_op_geom_101_sp + 6 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_z_0, to_z_x_0_x, to_z_x_0_y, to_z_x_0_z, to_z_xx, to_z_xy, to_z_xz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_0_x[k] = -2.0 * to_z_0[k] * tbe_0 + 4.0 * to_z_xx[k] * tbe_0 * tke_0;

            to_z_x_0_y[k] = 4.0 * to_z_xy[k] * tbe_0 * tke_0;

            to_z_x_0_z[k] = 4.0 * to_z_xz[k] * tbe_0 * tke_0;
        }

        // Set up 21-24 components of targeted buffer : SP

        auto to_z_y_0_x = pbuffer.data(idx_op_geom_101_sp + 7 * op_comps * 3 + i * 3 + 0);

        auto to_z_y_0_y = pbuffer.data(idx_op_geom_101_sp + 7 * op_comps * 3 + i * 3 + 1);

        auto to_z_y_0_z = pbuffer.data(idx_op_geom_101_sp + 7 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_z_0, to_z_xy, to_z_y_0_x, to_z_y_0_y, to_z_y_0_z, to_z_yy, to_z_yz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_0_x[k] = 4.0 * to_z_xy[k] * tbe_0 * tke_0;

            to_z_y_0_y[k] = -2.0 * to_z_0[k] * tbe_0 + 4.0 * to_z_yy[k] * tbe_0 * tke_0;

            to_z_y_0_z[k] = 4.0 * to_z_yz[k] * tbe_0 * tke_0;
        }

        // Set up 24-27 components of targeted buffer : SP

        auto to_z_z_0_x = pbuffer.data(idx_op_geom_101_sp + 8 * op_comps * 3 + i * 3 + 0);

        auto to_z_z_0_y = pbuffer.data(idx_op_geom_101_sp + 8 * op_comps * 3 + i * 3 + 1);

        auto to_z_z_0_z = pbuffer.data(idx_op_geom_101_sp + 8 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_z_0, to_z_xz, to_z_yz, to_z_z_0_x, to_z_z_0_y, to_z_z_0_z, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_0_x[k] = 4.0 * to_z_xz[k] * tbe_0 * tke_0;

            to_z_z_0_y[k] = 4.0 * to_z_yz[k] * tbe_0 * tke_0;

            to_z_z_0_z[k] = -2.0 * to_z_0[k] * tbe_0 + 4.0 * to_z_zz[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
