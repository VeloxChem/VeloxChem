#include "GeometricalDerivatives0X1ForPD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_pd(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_pd,
                       const int idx_op_pp,
                       const int idx_op_pf,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
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

        // Set up 0-6 components of targeted buffer : PD

        auto to_0_x_x_xx = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 0);

        auto to_0_x_x_xy = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 1);

        auto to_0_x_x_xz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 2);

        auto to_0_x_x_yy = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 3);

        auto to_0_x_x_yz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 4);

        auto to_0_x_x_zz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_0_x_x_xx, to_0_x_x_xy, to_0_x_x_xz, to_0_x_x_yy, to_0_x_x_yz, to_0_x_x_zz, to_x_x, to_x_xxx, to_x_xxy, to_x_xxz, to_x_xyy, to_x_xyz, to_x_xzz, to_x_y, to_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_x_xx[k] = -2.0 * to_x_x[k] + 2.0 * to_x_xxx[k] * tke_0;

            to_0_x_x_xy[k] = -to_x_y[k] + 2.0 * to_x_xxy[k] * tke_0;

            to_0_x_x_xz[k] = -to_x_z[k] + 2.0 * to_x_xxz[k] * tke_0;

            to_0_x_x_yy[k] = 2.0 * to_x_xyy[k] * tke_0;

            to_0_x_x_yz[k] = 2.0 * to_x_xyz[k] * tke_0;

            to_0_x_x_zz[k] = 2.0 * to_x_xzz[k] * tke_0;
        }

        // Set up 6-12 components of targeted buffer : PD

        auto to_0_x_y_xx = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 6);

        auto to_0_x_y_xy = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 7);

        auto to_0_x_y_xz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 8);

        auto to_0_x_y_yy = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 9);

        auto to_0_x_y_yz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 10);

        auto to_0_x_y_zz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_0_x_y_xx, to_0_x_y_xy, to_0_x_y_xz, to_0_x_y_yy, to_0_x_y_yz, to_0_x_y_zz, to_y_x, to_y_xxx, to_y_xxy, to_y_xxz, to_y_xyy, to_y_xyz, to_y_xzz, to_y_y, to_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_y_xx[k] = -2.0 * to_y_x[k] + 2.0 * to_y_xxx[k] * tke_0;

            to_0_x_y_xy[k] = -to_y_y[k] + 2.0 * to_y_xxy[k] * tke_0;

            to_0_x_y_xz[k] = -to_y_z[k] + 2.0 * to_y_xxz[k] * tke_0;

            to_0_x_y_yy[k] = 2.0 * to_y_xyy[k] * tke_0;

            to_0_x_y_yz[k] = 2.0 * to_y_xyz[k] * tke_0;

            to_0_x_y_zz[k] = 2.0 * to_y_xzz[k] * tke_0;
        }

        // Set up 12-18 components of targeted buffer : PD

        auto to_0_x_z_xx = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 12);

        auto to_0_x_z_xy = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 13);

        auto to_0_x_z_xz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 14);

        auto to_0_x_z_yy = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 15);

        auto to_0_x_z_yz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 16);

        auto to_0_x_z_zz = pbuffer.data(idx_op_geom_001_pd + 0 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_0_x_z_xx, to_0_x_z_xy, to_0_x_z_xz, to_0_x_z_yy, to_0_x_z_yz, to_0_x_z_zz, to_z_x, to_z_xxx, to_z_xxy, to_z_xxz, to_z_xyy, to_z_xyz, to_z_xzz, to_z_y, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_z_xx[k] = -2.0 * to_z_x[k] + 2.0 * to_z_xxx[k] * tke_0;

            to_0_x_z_xy[k] = -to_z_y[k] + 2.0 * to_z_xxy[k] * tke_0;

            to_0_x_z_xz[k] = -to_z_z[k] + 2.0 * to_z_xxz[k] * tke_0;

            to_0_x_z_yy[k] = 2.0 * to_z_xyy[k] * tke_0;

            to_0_x_z_yz[k] = 2.0 * to_z_xyz[k] * tke_0;

            to_0_x_z_zz[k] = 2.0 * to_z_xzz[k] * tke_0;
        }

        // Set up 18-24 components of targeted buffer : PD

        auto to_0_y_x_xx = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 0);

        auto to_0_y_x_xy = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 1);

        auto to_0_y_x_xz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 2);

        auto to_0_y_x_yy = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 3);

        auto to_0_y_x_yz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 4);

        auto to_0_y_x_zz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_0_y_x_xx, to_0_y_x_xy, to_0_y_x_xz, to_0_y_x_yy, to_0_y_x_yz, to_0_y_x_zz, to_x_x, to_x_xxy, to_x_xyy, to_x_xyz, to_x_y, to_x_yyy, to_x_yyz, to_x_yzz, to_x_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_x_xx[k] = 2.0 * to_x_xxy[k] * tke_0;

            to_0_y_x_xy[k] = -to_x_x[k] + 2.0 * to_x_xyy[k] * tke_0;

            to_0_y_x_xz[k] = 2.0 * to_x_xyz[k] * tke_0;

            to_0_y_x_yy[k] = -2.0 * to_x_y[k] + 2.0 * to_x_yyy[k] * tke_0;

            to_0_y_x_yz[k] = -to_x_z[k] + 2.0 * to_x_yyz[k] * tke_0;

            to_0_y_x_zz[k] = 2.0 * to_x_yzz[k] * tke_0;
        }

        // Set up 24-30 components of targeted buffer : PD

        auto to_0_y_y_xx = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 6);

        auto to_0_y_y_xy = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 7);

        auto to_0_y_y_xz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 8);

        auto to_0_y_y_yy = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 9);

        auto to_0_y_y_yz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 10);

        auto to_0_y_y_zz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_0_y_y_xx, to_0_y_y_xy, to_0_y_y_xz, to_0_y_y_yy, to_0_y_y_yz, to_0_y_y_zz, to_y_x, to_y_xxy, to_y_xyy, to_y_xyz, to_y_y, to_y_yyy, to_y_yyz, to_y_yzz, to_y_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_y_xx[k] = 2.0 * to_y_xxy[k] * tke_0;

            to_0_y_y_xy[k] = -to_y_x[k] + 2.0 * to_y_xyy[k] * tke_0;

            to_0_y_y_xz[k] = 2.0 * to_y_xyz[k] * tke_0;

            to_0_y_y_yy[k] = -2.0 * to_y_y[k] + 2.0 * to_y_yyy[k] * tke_0;

            to_0_y_y_yz[k] = -to_y_z[k] + 2.0 * to_y_yyz[k] * tke_0;

            to_0_y_y_zz[k] = 2.0 * to_y_yzz[k] * tke_0;
        }

        // Set up 30-36 components of targeted buffer : PD

        auto to_0_y_z_xx = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 12);

        auto to_0_y_z_xy = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 13);

        auto to_0_y_z_xz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 14);

        auto to_0_y_z_yy = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 15);

        auto to_0_y_z_yz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 16);

        auto to_0_y_z_zz = pbuffer.data(idx_op_geom_001_pd + 1 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_0_y_z_xx, to_0_y_z_xy, to_0_y_z_xz, to_0_y_z_yy, to_0_y_z_yz, to_0_y_z_zz, to_z_x, to_z_xxy, to_z_xyy, to_z_xyz, to_z_y, to_z_yyy, to_z_yyz, to_z_yzz, to_z_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_z_xx[k] = 2.0 * to_z_xxy[k] * tke_0;

            to_0_y_z_xy[k] = -to_z_x[k] + 2.0 * to_z_xyy[k] * tke_0;

            to_0_y_z_xz[k] = 2.0 * to_z_xyz[k] * tke_0;

            to_0_y_z_yy[k] = -2.0 * to_z_y[k] + 2.0 * to_z_yyy[k] * tke_0;

            to_0_y_z_yz[k] = -to_z_z[k] + 2.0 * to_z_yyz[k] * tke_0;

            to_0_y_z_zz[k] = 2.0 * to_z_yzz[k] * tke_0;
        }

        // Set up 36-42 components of targeted buffer : PD

        auto to_0_z_x_xx = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 0);

        auto to_0_z_x_xy = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 1);

        auto to_0_z_x_xz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 2);

        auto to_0_z_x_yy = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 3);

        auto to_0_z_x_yz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 4);

        auto to_0_z_x_zz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 5);

        #pragma omp simd aligned(to_0_z_x_xx, to_0_z_x_xy, to_0_z_x_xz, to_0_z_x_yy, to_0_z_x_yz, to_0_z_x_zz, to_x_x, to_x_xxz, to_x_xyz, to_x_xzz, to_x_y, to_x_yyz, to_x_yzz, to_x_z, to_x_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_x_xx[k] = 2.0 * to_x_xxz[k] * tke_0;

            to_0_z_x_xy[k] = 2.0 * to_x_xyz[k] * tke_0;

            to_0_z_x_xz[k] = -to_x_x[k] + 2.0 * to_x_xzz[k] * tke_0;

            to_0_z_x_yy[k] = 2.0 * to_x_yyz[k] * tke_0;

            to_0_z_x_yz[k] = -to_x_y[k] + 2.0 * to_x_yzz[k] * tke_0;

            to_0_z_x_zz[k] = -2.0 * to_x_z[k] + 2.0 * to_x_zzz[k] * tke_0;
        }

        // Set up 42-48 components of targeted buffer : PD

        auto to_0_z_y_xx = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 6);

        auto to_0_z_y_xy = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 7);

        auto to_0_z_y_xz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 8);

        auto to_0_z_y_yy = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 9);

        auto to_0_z_y_yz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 10);

        auto to_0_z_y_zz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 11);

        #pragma omp simd aligned(to_0_z_y_xx, to_0_z_y_xy, to_0_z_y_xz, to_0_z_y_yy, to_0_z_y_yz, to_0_z_y_zz, to_y_x, to_y_xxz, to_y_xyz, to_y_xzz, to_y_y, to_y_yyz, to_y_yzz, to_y_z, to_y_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_y_xx[k] = 2.0 * to_y_xxz[k] * tke_0;

            to_0_z_y_xy[k] = 2.0 * to_y_xyz[k] * tke_0;

            to_0_z_y_xz[k] = -to_y_x[k] + 2.0 * to_y_xzz[k] * tke_0;

            to_0_z_y_yy[k] = 2.0 * to_y_yyz[k] * tke_0;

            to_0_z_y_yz[k] = -to_y_y[k] + 2.0 * to_y_yzz[k] * tke_0;

            to_0_z_y_zz[k] = -2.0 * to_y_z[k] + 2.0 * to_y_zzz[k] * tke_0;
        }

        // Set up 48-54 components of targeted buffer : PD

        auto to_0_z_z_xx = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 12);

        auto to_0_z_z_xy = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 13);

        auto to_0_z_z_xz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 14);

        auto to_0_z_z_yy = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 15);

        auto to_0_z_z_yz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 16);

        auto to_0_z_z_zz = pbuffer.data(idx_op_geom_001_pd + 2 * op_comps * 18 + i * 18 + 17);

        #pragma omp simd aligned(to_0_z_z_xx, to_0_z_z_xy, to_0_z_z_xz, to_0_z_z_yy, to_0_z_z_yz, to_0_z_z_zz, to_z_x, to_z_xxz, to_z_xyz, to_z_xzz, to_z_y, to_z_yyz, to_z_yzz, to_z_z, to_z_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_z_xx[k] = 2.0 * to_z_xxz[k] * tke_0;

            to_0_z_z_xy[k] = 2.0 * to_z_xyz[k] * tke_0;

            to_0_z_z_xz[k] = -to_z_x[k] + 2.0 * to_z_xzz[k] * tke_0;

            to_0_z_z_yy[k] = 2.0 * to_z_yyz[k] * tke_0;

            to_0_z_z_yz[k] = -to_z_y[k] + 2.0 * to_z_yzz[k] * tke_0;

            to_0_z_z_zz[k] = -2.0 * to_z_z[k] + 2.0 * to_z_zzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

