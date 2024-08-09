#include "GeometricalDerivatives1X1ForSD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_sd(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_sd,
                        const size_t idx_op_pp,
                        const size_t idx_op_pf,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PP

        auto to_x_x = pbuffer.data(idx_op_pp + i * 9 + 0);

        auto to_x_y = pbuffer.data(idx_op_pp + i * 9 + 1);

        auto to_x_z = pbuffer.data(idx_op_pp + i * 9 + 2);

        auto to_y_x = pbuffer.data(idx_op_pp + i * 9 + 3);

        auto to_y_y = pbuffer.data(idx_op_pp + i * 9 + 4);

        auto to_y_z = pbuffer.data(idx_op_pp + i * 9 + 5);

        auto to_z_x = pbuffer.data(idx_op_pp + i * 9 + 6);

        auto to_z_y = pbuffer.data(idx_op_pp + i * 9 + 7);

        auto to_z_z = pbuffer.data(idx_op_pp + i * 9 + 8);

        // Set up components of auxiliary buffer : PF

        auto to_x_xxx = pbuffer.data(idx_op_pf + i * 30 + 0);

        auto to_x_xxy = pbuffer.data(idx_op_pf + i * 30 + 1);

        auto to_x_xxz = pbuffer.data(idx_op_pf + i * 30 + 2);

        auto to_x_xyy = pbuffer.data(idx_op_pf + i * 30 + 3);

        auto to_x_xyz = pbuffer.data(idx_op_pf + i * 30 + 4);

        auto to_x_xzz = pbuffer.data(idx_op_pf + i * 30 + 5);

        auto to_x_yyy = pbuffer.data(idx_op_pf + i * 30 + 6);

        auto to_x_yyz = pbuffer.data(idx_op_pf + i * 30 + 7);

        auto to_x_yzz = pbuffer.data(idx_op_pf + i * 30 + 8);

        auto to_x_zzz = pbuffer.data(idx_op_pf + i * 30 + 9);

        auto to_y_xxx = pbuffer.data(idx_op_pf + i * 30 + 10);

        auto to_y_xxy = pbuffer.data(idx_op_pf + i * 30 + 11);

        auto to_y_xxz = pbuffer.data(idx_op_pf + i * 30 + 12);

        auto to_y_xyy = pbuffer.data(idx_op_pf + i * 30 + 13);

        auto to_y_xyz = pbuffer.data(idx_op_pf + i * 30 + 14);

        auto to_y_xzz = pbuffer.data(idx_op_pf + i * 30 + 15);

        auto to_y_yyy = pbuffer.data(idx_op_pf + i * 30 + 16);

        auto to_y_yyz = pbuffer.data(idx_op_pf + i * 30 + 17);

        auto to_y_yzz = pbuffer.data(idx_op_pf + i * 30 + 18);

        auto to_y_zzz = pbuffer.data(idx_op_pf + i * 30 + 19);

        auto to_z_xxx = pbuffer.data(idx_op_pf + i * 30 + 20);

        auto to_z_xxy = pbuffer.data(idx_op_pf + i * 30 + 21);

        auto to_z_xxz = pbuffer.data(idx_op_pf + i * 30 + 22);

        auto to_z_xyy = pbuffer.data(idx_op_pf + i * 30 + 23);

        auto to_z_xyz = pbuffer.data(idx_op_pf + i * 30 + 24);

        auto to_z_xzz = pbuffer.data(idx_op_pf + i * 30 + 25);

        auto to_z_yyy = pbuffer.data(idx_op_pf + i * 30 + 26);

        auto to_z_yyz = pbuffer.data(idx_op_pf + i * 30 + 27);

        auto to_z_yzz = pbuffer.data(idx_op_pf + i * 30 + 28);

        auto to_z_zzz = pbuffer.data(idx_op_pf + i * 30 + 29);

        // Set up 0-6 components of targeted buffer : SD

        auto to_x_x_0_xx = pbuffer.data(idx_op_geom_101_sd + 0 * op_comps * 6 + i * 6 + 0);

        auto to_x_x_0_xy = pbuffer.data(idx_op_geom_101_sd + 0 * op_comps * 6 + i * 6 + 1);

        auto to_x_x_0_xz = pbuffer.data(idx_op_geom_101_sd + 0 * op_comps * 6 + i * 6 + 2);

        auto to_x_x_0_yy = pbuffer.data(idx_op_geom_101_sd + 0 * op_comps * 6 + i * 6 + 3);

        auto to_x_x_0_yz = pbuffer.data(idx_op_geom_101_sd + 0 * op_comps * 6 + i * 6 + 4);

        auto to_x_x_0_zz = pbuffer.data(idx_op_geom_101_sd + 0 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_x_x, to_x_x_0_xx, to_x_x_0_xy, to_x_x_0_xz, to_x_x_0_yy, to_x_x_0_yz, to_x_x_0_zz, to_x_xxx, to_x_xxy, to_x_xxz, to_x_xyy, to_x_xyz, to_x_xzz, to_x_y, to_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_0_xx[k] = -4.0 * to_x_x[k] * tbe_0 + 4.0 * to_x_xxx[k] * tbe_0 * tke_0;

            to_x_x_0_xy[k] = -2.0 * to_x_y[k] * tbe_0 + 4.0 * to_x_xxy[k] * tbe_0 * tke_0;

            to_x_x_0_xz[k] = -2.0 * to_x_z[k] * tbe_0 + 4.0 * to_x_xxz[k] * tbe_0 * tke_0;

            to_x_x_0_yy[k] = 4.0 * to_x_xyy[k] * tbe_0 * tke_0;

            to_x_x_0_yz[k] = 4.0 * to_x_xyz[k] * tbe_0 * tke_0;

            to_x_x_0_zz[k] = 4.0 * to_x_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 6-12 components of targeted buffer : SD

        auto to_x_y_0_xx = pbuffer.data(idx_op_geom_101_sd + 1 * op_comps * 6 + i * 6 + 0);

        auto to_x_y_0_xy = pbuffer.data(idx_op_geom_101_sd + 1 * op_comps * 6 + i * 6 + 1);

        auto to_x_y_0_xz = pbuffer.data(idx_op_geom_101_sd + 1 * op_comps * 6 + i * 6 + 2);

        auto to_x_y_0_yy = pbuffer.data(idx_op_geom_101_sd + 1 * op_comps * 6 + i * 6 + 3);

        auto to_x_y_0_yz = pbuffer.data(idx_op_geom_101_sd + 1 * op_comps * 6 + i * 6 + 4);

        auto to_x_y_0_zz = pbuffer.data(idx_op_geom_101_sd + 1 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_x_x, to_x_xxy, to_x_xyy, to_x_xyz, to_x_y, to_x_y_0_xx, to_x_y_0_xy, to_x_y_0_xz, to_x_y_0_yy, to_x_y_0_yz, to_x_y_0_zz, to_x_yyy, to_x_yyz, to_x_yzz, to_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_0_xx[k] = 4.0 * to_x_xxy[k] * tbe_0 * tke_0;

            to_x_y_0_xy[k] = -2.0 * to_x_x[k] * tbe_0 + 4.0 * to_x_xyy[k] * tbe_0 * tke_0;

            to_x_y_0_xz[k] = 4.0 * to_x_xyz[k] * tbe_0 * tke_0;

            to_x_y_0_yy[k] = -4.0 * to_x_y[k] * tbe_0 + 4.0 * to_x_yyy[k] * tbe_0 * tke_0;

            to_x_y_0_yz[k] = -2.0 * to_x_z[k] * tbe_0 + 4.0 * to_x_yyz[k] * tbe_0 * tke_0;

            to_x_y_0_zz[k] = 4.0 * to_x_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 12-18 components of targeted buffer : SD

        auto to_x_z_0_xx = pbuffer.data(idx_op_geom_101_sd + 2 * op_comps * 6 + i * 6 + 0);

        auto to_x_z_0_xy = pbuffer.data(idx_op_geom_101_sd + 2 * op_comps * 6 + i * 6 + 1);

        auto to_x_z_0_xz = pbuffer.data(idx_op_geom_101_sd + 2 * op_comps * 6 + i * 6 + 2);

        auto to_x_z_0_yy = pbuffer.data(idx_op_geom_101_sd + 2 * op_comps * 6 + i * 6 + 3);

        auto to_x_z_0_yz = pbuffer.data(idx_op_geom_101_sd + 2 * op_comps * 6 + i * 6 + 4);

        auto to_x_z_0_zz = pbuffer.data(idx_op_geom_101_sd + 2 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_x_x, to_x_xxz, to_x_xyz, to_x_xzz, to_x_y, to_x_yyz, to_x_yzz, to_x_z, to_x_z_0_xx, to_x_z_0_xy, to_x_z_0_xz, to_x_z_0_yy, to_x_z_0_yz, to_x_z_0_zz, to_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_0_xx[k] = 4.0 * to_x_xxz[k] * tbe_0 * tke_0;

            to_x_z_0_xy[k] = 4.0 * to_x_xyz[k] * tbe_0 * tke_0;

            to_x_z_0_xz[k] = -2.0 * to_x_x[k] * tbe_0 + 4.0 * to_x_xzz[k] * tbe_0 * tke_0;

            to_x_z_0_yy[k] = 4.0 * to_x_yyz[k] * tbe_0 * tke_0;

            to_x_z_0_yz[k] = -2.0 * to_x_y[k] * tbe_0 + 4.0 * to_x_yzz[k] * tbe_0 * tke_0;

            to_x_z_0_zz[k] = -4.0 * to_x_z[k] * tbe_0 + 4.0 * to_x_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 18-24 components of targeted buffer : SD

        auto to_y_x_0_xx = pbuffer.data(idx_op_geom_101_sd + 3 * op_comps * 6 + i * 6 + 0);

        auto to_y_x_0_xy = pbuffer.data(idx_op_geom_101_sd + 3 * op_comps * 6 + i * 6 + 1);

        auto to_y_x_0_xz = pbuffer.data(idx_op_geom_101_sd + 3 * op_comps * 6 + i * 6 + 2);

        auto to_y_x_0_yy = pbuffer.data(idx_op_geom_101_sd + 3 * op_comps * 6 + i * 6 + 3);

        auto to_y_x_0_yz = pbuffer.data(idx_op_geom_101_sd + 3 * op_comps * 6 + i * 6 + 4);

        auto to_y_x_0_zz = pbuffer.data(idx_op_geom_101_sd + 3 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_y_x, to_y_x_0_xx, to_y_x_0_xy, to_y_x_0_xz, to_y_x_0_yy, to_y_x_0_yz, to_y_x_0_zz, to_y_xxx, to_y_xxy, to_y_xxz, to_y_xyy, to_y_xyz, to_y_xzz, to_y_y, to_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_0_xx[k] = -4.0 * to_y_x[k] * tbe_0 + 4.0 * to_y_xxx[k] * tbe_0 * tke_0;

            to_y_x_0_xy[k] = -2.0 * to_y_y[k] * tbe_0 + 4.0 * to_y_xxy[k] * tbe_0 * tke_0;

            to_y_x_0_xz[k] = -2.0 * to_y_z[k] * tbe_0 + 4.0 * to_y_xxz[k] * tbe_0 * tke_0;

            to_y_x_0_yy[k] = 4.0 * to_y_xyy[k] * tbe_0 * tke_0;

            to_y_x_0_yz[k] = 4.0 * to_y_xyz[k] * tbe_0 * tke_0;

            to_y_x_0_zz[k] = 4.0 * to_y_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 24-30 components of targeted buffer : SD

        auto to_y_y_0_xx = pbuffer.data(idx_op_geom_101_sd + 4 * op_comps * 6 + i * 6 + 0);

        auto to_y_y_0_xy = pbuffer.data(idx_op_geom_101_sd + 4 * op_comps * 6 + i * 6 + 1);

        auto to_y_y_0_xz = pbuffer.data(idx_op_geom_101_sd + 4 * op_comps * 6 + i * 6 + 2);

        auto to_y_y_0_yy = pbuffer.data(idx_op_geom_101_sd + 4 * op_comps * 6 + i * 6 + 3);

        auto to_y_y_0_yz = pbuffer.data(idx_op_geom_101_sd + 4 * op_comps * 6 + i * 6 + 4);

        auto to_y_y_0_zz = pbuffer.data(idx_op_geom_101_sd + 4 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_y_x, to_y_xxy, to_y_xyy, to_y_xyz, to_y_y, to_y_y_0_xx, to_y_y_0_xy, to_y_y_0_xz, to_y_y_0_yy, to_y_y_0_yz, to_y_y_0_zz, to_y_yyy, to_y_yyz, to_y_yzz, to_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_0_xx[k] = 4.0 * to_y_xxy[k] * tbe_0 * tke_0;

            to_y_y_0_xy[k] = -2.0 * to_y_x[k] * tbe_0 + 4.0 * to_y_xyy[k] * tbe_0 * tke_0;

            to_y_y_0_xz[k] = 4.0 * to_y_xyz[k] * tbe_0 * tke_0;

            to_y_y_0_yy[k] = -4.0 * to_y_y[k] * tbe_0 + 4.0 * to_y_yyy[k] * tbe_0 * tke_0;

            to_y_y_0_yz[k] = -2.0 * to_y_z[k] * tbe_0 + 4.0 * to_y_yyz[k] * tbe_0 * tke_0;

            to_y_y_0_zz[k] = 4.0 * to_y_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-36 components of targeted buffer : SD

        auto to_y_z_0_xx = pbuffer.data(idx_op_geom_101_sd + 5 * op_comps * 6 + i * 6 + 0);

        auto to_y_z_0_xy = pbuffer.data(idx_op_geom_101_sd + 5 * op_comps * 6 + i * 6 + 1);

        auto to_y_z_0_xz = pbuffer.data(idx_op_geom_101_sd + 5 * op_comps * 6 + i * 6 + 2);

        auto to_y_z_0_yy = pbuffer.data(idx_op_geom_101_sd + 5 * op_comps * 6 + i * 6 + 3);

        auto to_y_z_0_yz = pbuffer.data(idx_op_geom_101_sd + 5 * op_comps * 6 + i * 6 + 4);

        auto to_y_z_0_zz = pbuffer.data(idx_op_geom_101_sd + 5 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_y_x, to_y_xxz, to_y_xyz, to_y_xzz, to_y_y, to_y_yyz, to_y_yzz, to_y_z, to_y_z_0_xx, to_y_z_0_xy, to_y_z_0_xz, to_y_z_0_yy, to_y_z_0_yz, to_y_z_0_zz, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_0_xx[k] = 4.0 * to_y_xxz[k] * tbe_0 * tke_0;

            to_y_z_0_xy[k] = 4.0 * to_y_xyz[k] * tbe_0 * tke_0;

            to_y_z_0_xz[k] = -2.0 * to_y_x[k] * tbe_0 + 4.0 * to_y_xzz[k] * tbe_0 * tke_0;

            to_y_z_0_yy[k] = 4.0 * to_y_yyz[k] * tbe_0 * tke_0;

            to_y_z_0_yz[k] = -2.0 * to_y_y[k] * tbe_0 + 4.0 * to_y_yzz[k] * tbe_0 * tke_0;

            to_y_z_0_zz[k] = -4.0 * to_y_z[k] * tbe_0 + 4.0 * to_y_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 36-42 components of targeted buffer : SD

        auto to_z_x_0_xx = pbuffer.data(idx_op_geom_101_sd + 6 * op_comps * 6 + i * 6 + 0);

        auto to_z_x_0_xy = pbuffer.data(idx_op_geom_101_sd + 6 * op_comps * 6 + i * 6 + 1);

        auto to_z_x_0_xz = pbuffer.data(idx_op_geom_101_sd + 6 * op_comps * 6 + i * 6 + 2);

        auto to_z_x_0_yy = pbuffer.data(idx_op_geom_101_sd + 6 * op_comps * 6 + i * 6 + 3);

        auto to_z_x_0_yz = pbuffer.data(idx_op_geom_101_sd + 6 * op_comps * 6 + i * 6 + 4);

        auto to_z_x_0_zz = pbuffer.data(idx_op_geom_101_sd + 6 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_z_x, to_z_x_0_xx, to_z_x_0_xy, to_z_x_0_xz, to_z_x_0_yy, to_z_x_0_yz, to_z_x_0_zz, to_z_xxx, to_z_xxy, to_z_xxz, to_z_xyy, to_z_xyz, to_z_xzz, to_z_y, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_0_xx[k] = -4.0 * to_z_x[k] * tbe_0 + 4.0 * to_z_xxx[k] * tbe_0 * tke_0;

            to_z_x_0_xy[k] = -2.0 * to_z_y[k] * tbe_0 + 4.0 * to_z_xxy[k] * tbe_0 * tke_0;

            to_z_x_0_xz[k] = -2.0 * to_z_z[k] * tbe_0 + 4.0 * to_z_xxz[k] * tbe_0 * tke_0;

            to_z_x_0_yy[k] = 4.0 * to_z_xyy[k] * tbe_0 * tke_0;

            to_z_x_0_yz[k] = 4.0 * to_z_xyz[k] * tbe_0 * tke_0;

            to_z_x_0_zz[k] = 4.0 * to_z_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 42-48 components of targeted buffer : SD

        auto to_z_y_0_xx = pbuffer.data(idx_op_geom_101_sd + 7 * op_comps * 6 + i * 6 + 0);

        auto to_z_y_0_xy = pbuffer.data(idx_op_geom_101_sd + 7 * op_comps * 6 + i * 6 + 1);

        auto to_z_y_0_xz = pbuffer.data(idx_op_geom_101_sd + 7 * op_comps * 6 + i * 6 + 2);

        auto to_z_y_0_yy = pbuffer.data(idx_op_geom_101_sd + 7 * op_comps * 6 + i * 6 + 3);

        auto to_z_y_0_yz = pbuffer.data(idx_op_geom_101_sd + 7 * op_comps * 6 + i * 6 + 4);

        auto to_z_y_0_zz = pbuffer.data(idx_op_geom_101_sd + 7 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_z_x, to_z_xxy, to_z_xyy, to_z_xyz, to_z_y, to_z_y_0_xx, to_z_y_0_xy, to_z_y_0_xz, to_z_y_0_yy, to_z_y_0_yz, to_z_y_0_zz, to_z_yyy, to_z_yyz, to_z_yzz, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_0_xx[k] = 4.0 * to_z_xxy[k] * tbe_0 * tke_0;

            to_z_y_0_xy[k] = -2.0 * to_z_x[k] * tbe_0 + 4.0 * to_z_xyy[k] * tbe_0 * tke_0;

            to_z_y_0_xz[k] = 4.0 * to_z_xyz[k] * tbe_0 * tke_0;

            to_z_y_0_yy[k] = -4.0 * to_z_y[k] * tbe_0 + 4.0 * to_z_yyy[k] * tbe_0 * tke_0;

            to_z_y_0_yz[k] = -2.0 * to_z_z[k] * tbe_0 + 4.0 * to_z_yyz[k] * tbe_0 * tke_0;

            to_z_y_0_zz[k] = 4.0 * to_z_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 48-54 components of targeted buffer : SD

        auto to_z_z_0_xx = pbuffer.data(idx_op_geom_101_sd + 8 * op_comps * 6 + i * 6 + 0);

        auto to_z_z_0_xy = pbuffer.data(idx_op_geom_101_sd + 8 * op_comps * 6 + i * 6 + 1);

        auto to_z_z_0_xz = pbuffer.data(idx_op_geom_101_sd + 8 * op_comps * 6 + i * 6 + 2);

        auto to_z_z_0_yy = pbuffer.data(idx_op_geom_101_sd + 8 * op_comps * 6 + i * 6 + 3);

        auto to_z_z_0_yz = pbuffer.data(idx_op_geom_101_sd + 8 * op_comps * 6 + i * 6 + 4);

        auto to_z_z_0_zz = pbuffer.data(idx_op_geom_101_sd + 8 * op_comps * 6 + i * 6 + 5);

        #pragma omp simd aligned(to_z_x, to_z_xxz, to_z_xyz, to_z_xzz, to_z_y, to_z_yyz, to_z_yzz, to_z_z, to_z_z_0_xx, to_z_z_0_xy, to_z_z_0_xz, to_z_z_0_yy, to_z_z_0_yz, to_z_z_0_zz, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_0_xx[k] = 4.0 * to_z_xxz[k] * tbe_0 * tke_0;

            to_z_z_0_xy[k] = 4.0 * to_z_xyz[k] * tbe_0 * tke_0;

            to_z_z_0_xz[k] = -2.0 * to_z_x[k] * tbe_0 + 4.0 * to_z_xzz[k] * tbe_0 * tke_0;

            to_z_z_0_yy[k] = 4.0 * to_z_yyz[k] * tbe_0 * tke_0;

            to_z_z_0_yz[k] = -2.0 * to_z_y[k] * tbe_0 + 4.0 * to_z_yzz[k] * tbe_0 * tke_0;

            to_z_z_0_zz[k] = -4.0 * to_z_z[k] * tbe_0 + 4.0 * to_z_zzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

