#include "GeometricalDerivatives1X1ForSF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_sf(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_sf,
                        const size_t idx_op_pd,
                        const size_t idx_op_pg,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
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

        // Set up components of auxiliary buffer : PG

        auto to_x_xxxx = pbuffer.data(idx_op_pg + i * 45 + 0);

        auto to_x_xxxy = pbuffer.data(idx_op_pg + i * 45 + 1);

        auto to_x_xxxz = pbuffer.data(idx_op_pg + i * 45 + 2);

        auto to_x_xxyy = pbuffer.data(idx_op_pg + i * 45 + 3);

        auto to_x_xxyz = pbuffer.data(idx_op_pg + i * 45 + 4);

        auto to_x_xxzz = pbuffer.data(idx_op_pg + i * 45 + 5);

        auto to_x_xyyy = pbuffer.data(idx_op_pg + i * 45 + 6);

        auto to_x_xyyz = pbuffer.data(idx_op_pg + i * 45 + 7);

        auto to_x_xyzz = pbuffer.data(idx_op_pg + i * 45 + 8);

        auto to_x_xzzz = pbuffer.data(idx_op_pg + i * 45 + 9);

        auto to_x_yyyy = pbuffer.data(idx_op_pg + i * 45 + 10);

        auto to_x_yyyz = pbuffer.data(idx_op_pg + i * 45 + 11);

        auto to_x_yyzz = pbuffer.data(idx_op_pg + i * 45 + 12);

        auto to_x_yzzz = pbuffer.data(idx_op_pg + i * 45 + 13);

        auto to_x_zzzz = pbuffer.data(idx_op_pg + i * 45 + 14);

        auto to_y_xxxx = pbuffer.data(idx_op_pg + i * 45 + 15);

        auto to_y_xxxy = pbuffer.data(idx_op_pg + i * 45 + 16);

        auto to_y_xxxz = pbuffer.data(idx_op_pg + i * 45 + 17);

        auto to_y_xxyy = pbuffer.data(idx_op_pg + i * 45 + 18);

        auto to_y_xxyz = pbuffer.data(idx_op_pg + i * 45 + 19);

        auto to_y_xxzz = pbuffer.data(idx_op_pg + i * 45 + 20);

        auto to_y_xyyy = pbuffer.data(idx_op_pg + i * 45 + 21);

        auto to_y_xyyz = pbuffer.data(idx_op_pg + i * 45 + 22);

        auto to_y_xyzz = pbuffer.data(idx_op_pg + i * 45 + 23);

        auto to_y_xzzz = pbuffer.data(idx_op_pg + i * 45 + 24);

        auto to_y_yyyy = pbuffer.data(idx_op_pg + i * 45 + 25);

        auto to_y_yyyz = pbuffer.data(idx_op_pg + i * 45 + 26);

        auto to_y_yyzz = pbuffer.data(idx_op_pg + i * 45 + 27);

        auto to_y_yzzz = pbuffer.data(idx_op_pg + i * 45 + 28);

        auto to_y_zzzz = pbuffer.data(idx_op_pg + i * 45 + 29);

        auto to_z_xxxx = pbuffer.data(idx_op_pg + i * 45 + 30);

        auto to_z_xxxy = pbuffer.data(idx_op_pg + i * 45 + 31);

        auto to_z_xxxz = pbuffer.data(idx_op_pg + i * 45 + 32);

        auto to_z_xxyy = pbuffer.data(idx_op_pg + i * 45 + 33);

        auto to_z_xxyz = pbuffer.data(idx_op_pg + i * 45 + 34);

        auto to_z_xxzz = pbuffer.data(idx_op_pg + i * 45 + 35);

        auto to_z_xyyy = pbuffer.data(idx_op_pg + i * 45 + 36);

        auto to_z_xyyz = pbuffer.data(idx_op_pg + i * 45 + 37);

        auto to_z_xyzz = pbuffer.data(idx_op_pg + i * 45 + 38);

        auto to_z_xzzz = pbuffer.data(idx_op_pg + i * 45 + 39);

        auto to_z_yyyy = pbuffer.data(idx_op_pg + i * 45 + 40);

        auto to_z_yyyz = pbuffer.data(idx_op_pg + i * 45 + 41);

        auto to_z_yyzz = pbuffer.data(idx_op_pg + i * 45 + 42);

        auto to_z_yzzz = pbuffer.data(idx_op_pg + i * 45 + 43);

        auto to_z_zzzz = pbuffer.data(idx_op_pg + i * 45 + 44);

        // Set up 0-10 components of targeted buffer : SF

        auto to_x_x_0_xxx = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 0);

        auto to_x_x_0_xxy = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 1);

        auto to_x_x_0_xxz = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 2);

        auto to_x_x_0_xyy = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 3);

        auto to_x_x_0_xyz = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 4);

        auto to_x_x_0_xzz = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 5);

        auto to_x_x_0_yyy = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 6);

        auto to_x_x_0_yyz = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 7);

        auto to_x_x_0_yzz = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 8);

        auto to_x_x_0_zzz = pbuffer.data(idx_op_geom_101_sf + 0 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_x_x_0_xxx, to_x_x_0_xxy, to_x_x_0_xxz, to_x_x_0_xyy, to_x_x_0_xyz, to_x_x_0_xzz, to_x_x_0_yyy, to_x_x_0_yyz, to_x_x_0_yzz, to_x_x_0_zzz, to_x_xx, to_x_xxxx, to_x_xxxy, to_x_xxxz, to_x_xxyy, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yz, to_x_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_0_xxx[k] = -6.0 * to_x_xx[k] * tbe_0 + 4.0 * to_x_xxxx[k] * tbe_0 * tke_0;

            to_x_x_0_xxy[k] = -4.0 * to_x_xy[k] * tbe_0 + 4.0 * to_x_xxxy[k] * tbe_0 * tke_0;

            to_x_x_0_xxz[k] = -4.0 * to_x_xz[k] * tbe_0 + 4.0 * to_x_xxxz[k] * tbe_0 * tke_0;

            to_x_x_0_xyy[k] = -2.0 * to_x_yy[k] * tbe_0 + 4.0 * to_x_xxyy[k] * tbe_0 * tke_0;

            to_x_x_0_xyz[k] = -2.0 * to_x_yz[k] * tbe_0 + 4.0 * to_x_xxyz[k] * tbe_0 * tke_0;

            to_x_x_0_xzz[k] = -2.0 * to_x_zz[k] * tbe_0 + 4.0 * to_x_xxzz[k] * tbe_0 * tke_0;

            to_x_x_0_yyy[k] = 4.0 * to_x_xyyy[k] * tbe_0 * tke_0;

            to_x_x_0_yyz[k] = 4.0 * to_x_xyyz[k] * tbe_0 * tke_0;

            to_x_x_0_yzz[k] = 4.0 * to_x_xyzz[k] * tbe_0 * tke_0;

            to_x_x_0_zzz[k] = 4.0 * to_x_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 10-20 components of targeted buffer : SF

        auto to_x_y_0_xxx = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 0);

        auto to_x_y_0_xxy = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 1);

        auto to_x_y_0_xxz = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 2);

        auto to_x_y_0_xyy = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 3);

        auto to_x_y_0_xyz = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 4);

        auto to_x_y_0_xzz = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 5);

        auto to_x_y_0_yyy = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 6);

        auto to_x_y_0_yyz = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 7);

        auto to_x_y_0_yzz = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 8);

        auto to_x_y_0_zzz = pbuffer.data(idx_op_geom_101_sf + 1 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_x_xx, to_x_xxxy, to_x_xxyy, to_x_xxyz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_y_0_xxx, to_x_y_0_xxy, to_x_y_0_xxz, to_x_y_0_xyy, to_x_y_0_xyz, to_x_y_0_xzz, to_x_y_0_yyy, to_x_y_0_yyz, to_x_y_0_yzz, to_x_y_0_zzz, to_x_yy, to_x_yyyy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_0_xxx[k] = 4.0 * to_x_xxxy[k] * tbe_0 * tke_0;

            to_x_y_0_xxy[k] = -2.0 * to_x_xx[k] * tbe_0 + 4.0 * to_x_xxyy[k] * tbe_0 * tke_0;

            to_x_y_0_xxz[k] = 4.0 * to_x_xxyz[k] * tbe_0 * tke_0;

            to_x_y_0_xyy[k] = -4.0 * to_x_xy[k] * tbe_0 + 4.0 * to_x_xyyy[k] * tbe_0 * tke_0;

            to_x_y_0_xyz[k] = -2.0 * to_x_xz[k] * tbe_0 + 4.0 * to_x_xyyz[k] * tbe_0 * tke_0;

            to_x_y_0_xzz[k] = 4.0 * to_x_xyzz[k] * tbe_0 * tke_0;

            to_x_y_0_yyy[k] = -6.0 * to_x_yy[k] * tbe_0 + 4.0 * to_x_yyyy[k] * tbe_0 * tke_0;

            to_x_y_0_yyz[k] = -4.0 * to_x_yz[k] * tbe_0 + 4.0 * to_x_yyyz[k] * tbe_0 * tke_0;

            to_x_y_0_yzz[k] = -2.0 * to_x_zz[k] * tbe_0 + 4.0 * to_x_yyzz[k] * tbe_0 * tke_0;

            to_x_y_0_zzz[k] = 4.0 * to_x_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 20-30 components of targeted buffer : SF

        auto to_x_z_0_xxx = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 0);

        auto to_x_z_0_xxy = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 1);

        auto to_x_z_0_xxz = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 2);

        auto to_x_z_0_xyy = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 3);

        auto to_x_z_0_xyz = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 4);

        auto to_x_z_0_xzz = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 5);

        auto to_x_z_0_yyy = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 6);

        auto to_x_z_0_yyz = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 7);

        auto to_x_z_0_yzz = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 8);

        auto to_x_z_0_zzz = pbuffer.data(idx_op_geom_101_sf + 2 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_x_xx, to_x_xxxz, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_z_0_xxx, to_x_z_0_xxy, to_x_z_0_xxz, to_x_z_0_xyy, to_x_z_0_xyz, to_x_z_0_xzz, to_x_z_0_yyy, to_x_z_0_yyz, to_x_z_0_yzz, to_x_z_0_zzz, to_x_zz, to_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_0_xxx[k] = 4.0 * to_x_xxxz[k] * tbe_0 * tke_0;

            to_x_z_0_xxy[k] = 4.0 * to_x_xxyz[k] * tbe_0 * tke_0;

            to_x_z_0_xxz[k] = -2.0 * to_x_xx[k] * tbe_0 + 4.0 * to_x_xxzz[k] * tbe_0 * tke_0;

            to_x_z_0_xyy[k] = 4.0 * to_x_xyyz[k] * tbe_0 * tke_0;

            to_x_z_0_xyz[k] = -2.0 * to_x_xy[k] * tbe_0 + 4.0 * to_x_xyzz[k] * tbe_0 * tke_0;

            to_x_z_0_xzz[k] = -4.0 * to_x_xz[k] * tbe_0 + 4.0 * to_x_xzzz[k] * tbe_0 * tke_0;

            to_x_z_0_yyy[k] = 4.0 * to_x_yyyz[k] * tbe_0 * tke_0;

            to_x_z_0_yyz[k] = -2.0 * to_x_yy[k] * tbe_0 + 4.0 * to_x_yyzz[k] * tbe_0 * tke_0;

            to_x_z_0_yzz[k] = -4.0 * to_x_yz[k] * tbe_0 + 4.0 * to_x_yzzz[k] * tbe_0 * tke_0;

            to_x_z_0_zzz[k] = -6.0 * to_x_zz[k] * tbe_0 + 4.0 * to_x_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-40 components of targeted buffer : SF

        auto to_y_x_0_xxx = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 0);

        auto to_y_x_0_xxy = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 1);

        auto to_y_x_0_xxz = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 2);

        auto to_y_x_0_xyy = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 3);

        auto to_y_x_0_xyz = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 4);

        auto to_y_x_0_xzz = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 5);

        auto to_y_x_0_yyy = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 6);

        auto to_y_x_0_yyz = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 7);

        auto to_y_x_0_yzz = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 8);

        auto to_y_x_0_zzz = pbuffer.data(idx_op_geom_101_sf + 3 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_y_x_0_xxx, to_y_x_0_xxy, to_y_x_0_xxz, to_y_x_0_xyy, to_y_x_0_xyz, to_y_x_0_xzz, to_y_x_0_yyy, to_y_x_0_yyz, to_y_x_0_yzz, to_y_x_0_zzz, to_y_xx, to_y_xxxx, to_y_xxxy, to_y_xxxz, to_y_xxyy, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_0_xxx[k] = -6.0 * to_y_xx[k] * tbe_0 + 4.0 * to_y_xxxx[k] * tbe_0 * tke_0;

            to_y_x_0_xxy[k] = -4.0 * to_y_xy[k] * tbe_0 + 4.0 * to_y_xxxy[k] * tbe_0 * tke_0;

            to_y_x_0_xxz[k] = -4.0 * to_y_xz[k] * tbe_0 + 4.0 * to_y_xxxz[k] * tbe_0 * tke_0;

            to_y_x_0_xyy[k] = -2.0 * to_y_yy[k] * tbe_0 + 4.0 * to_y_xxyy[k] * tbe_0 * tke_0;

            to_y_x_0_xyz[k] = -2.0 * to_y_yz[k] * tbe_0 + 4.0 * to_y_xxyz[k] * tbe_0 * tke_0;

            to_y_x_0_xzz[k] = -2.0 * to_y_zz[k] * tbe_0 + 4.0 * to_y_xxzz[k] * tbe_0 * tke_0;

            to_y_x_0_yyy[k] = 4.0 * to_y_xyyy[k] * tbe_0 * tke_0;

            to_y_x_0_yyz[k] = 4.0 * to_y_xyyz[k] * tbe_0 * tke_0;

            to_y_x_0_yzz[k] = 4.0 * to_y_xyzz[k] * tbe_0 * tke_0;

            to_y_x_0_zzz[k] = 4.0 * to_y_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 40-50 components of targeted buffer : SF

        auto to_y_y_0_xxx = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 0);

        auto to_y_y_0_xxy = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 1);

        auto to_y_y_0_xxz = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 2);

        auto to_y_y_0_xyy = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 3);

        auto to_y_y_0_xyz = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 4);

        auto to_y_y_0_xzz = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 5);

        auto to_y_y_0_yyy = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 6);

        auto to_y_y_0_yyz = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 7);

        auto to_y_y_0_yzz = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 8);

        auto to_y_y_0_zzz = pbuffer.data(idx_op_geom_101_sf + 4 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_y_xx, to_y_xxxy, to_y_xxyy, to_y_xxyz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_y_0_xxx, to_y_y_0_xxy, to_y_y_0_xxz, to_y_y_0_xyy, to_y_y_0_xyz, to_y_y_0_xzz, to_y_y_0_yyy, to_y_y_0_yyz, to_y_y_0_yzz, to_y_y_0_zzz, to_y_yy, to_y_yyyy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_0_xxx[k] = 4.0 * to_y_xxxy[k] * tbe_0 * tke_0;

            to_y_y_0_xxy[k] = -2.0 * to_y_xx[k] * tbe_0 + 4.0 * to_y_xxyy[k] * tbe_0 * tke_0;

            to_y_y_0_xxz[k] = 4.0 * to_y_xxyz[k] * tbe_0 * tke_0;

            to_y_y_0_xyy[k] = -4.0 * to_y_xy[k] * tbe_0 + 4.0 * to_y_xyyy[k] * tbe_0 * tke_0;

            to_y_y_0_xyz[k] = -2.0 * to_y_xz[k] * tbe_0 + 4.0 * to_y_xyyz[k] * tbe_0 * tke_0;

            to_y_y_0_xzz[k] = 4.0 * to_y_xyzz[k] * tbe_0 * tke_0;

            to_y_y_0_yyy[k] = -6.0 * to_y_yy[k] * tbe_0 + 4.0 * to_y_yyyy[k] * tbe_0 * tke_0;

            to_y_y_0_yyz[k] = -4.0 * to_y_yz[k] * tbe_0 + 4.0 * to_y_yyyz[k] * tbe_0 * tke_0;

            to_y_y_0_yzz[k] = -2.0 * to_y_zz[k] * tbe_0 + 4.0 * to_y_yyzz[k] * tbe_0 * tke_0;

            to_y_y_0_zzz[k] = 4.0 * to_y_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 50-60 components of targeted buffer : SF

        auto to_y_z_0_xxx = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 0);

        auto to_y_z_0_xxy = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 1);

        auto to_y_z_0_xxz = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 2);

        auto to_y_z_0_xyy = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 3);

        auto to_y_z_0_xyz = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 4);

        auto to_y_z_0_xzz = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 5);

        auto to_y_z_0_yyy = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 6);

        auto to_y_z_0_yyz = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 7);

        auto to_y_z_0_yzz = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 8);

        auto to_y_z_0_zzz = pbuffer.data(idx_op_geom_101_sf + 5 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_y_xx, to_y_xxxz, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_z_0_xxx, to_y_z_0_xxy, to_y_z_0_xxz, to_y_z_0_xyy, to_y_z_0_xyz, to_y_z_0_xzz, to_y_z_0_yyy, to_y_z_0_yyz, to_y_z_0_yzz, to_y_z_0_zzz, to_y_zz, to_y_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_0_xxx[k] = 4.0 * to_y_xxxz[k] * tbe_0 * tke_0;

            to_y_z_0_xxy[k] = 4.0 * to_y_xxyz[k] * tbe_0 * tke_0;

            to_y_z_0_xxz[k] = -2.0 * to_y_xx[k] * tbe_0 + 4.0 * to_y_xxzz[k] * tbe_0 * tke_0;

            to_y_z_0_xyy[k] = 4.0 * to_y_xyyz[k] * tbe_0 * tke_0;

            to_y_z_0_xyz[k] = -2.0 * to_y_xy[k] * tbe_0 + 4.0 * to_y_xyzz[k] * tbe_0 * tke_0;

            to_y_z_0_xzz[k] = -4.0 * to_y_xz[k] * tbe_0 + 4.0 * to_y_xzzz[k] * tbe_0 * tke_0;

            to_y_z_0_yyy[k] = 4.0 * to_y_yyyz[k] * tbe_0 * tke_0;

            to_y_z_0_yyz[k] = -2.0 * to_y_yy[k] * tbe_0 + 4.0 * to_y_yyzz[k] * tbe_0 * tke_0;

            to_y_z_0_yzz[k] = -4.0 * to_y_yz[k] * tbe_0 + 4.0 * to_y_yzzz[k] * tbe_0 * tke_0;

            to_y_z_0_zzz[k] = -6.0 * to_y_zz[k] * tbe_0 + 4.0 * to_y_zzzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-70 components of targeted buffer : SF

        auto to_z_x_0_xxx = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 0);

        auto to_z_x_0_xxy = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 1);

        auto to_z_x_0_xxz = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 2);

        auto to_z_x_0_xyy = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 3);

        auto to_z_x_0_xyz = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 4);

        auto to_z_x_0_xzz = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 5);

        auto to_z_x_0_yyy = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 6);

        auto to_z_x_0_yyz = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 7);

        auto to_z_x_0_yzz = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 8);

        auto to_z_x_0_zzz = pbuffer.data(idx_op_geom_101_sf + 6 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_z_x_0_xxx, to_z_x_0_xxy, to_z_x_0_xxz, to_z_x_0_xyy, to_z_x_0_xyz, to_z_x_0_xzz, to_z_x_0_yyy, to_z_x_0_yyz, to_z_x_0_yzz, to_z_x_0_zzz, to_z_xx, to_z_xxxx, to_z_xxxy, to_z_xxxz, to_z_xxyy, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_0_xxx[k] = -6.0 * to_z_xx[k] * tbe_0 + 4.0 * to_z_xxxx[k] * tbe_0 * tke_0;

            to_z_x_0_xxy[k] = -4.0 * to_z_xy[k] * tbe_0 + 4.0 * to_z_xxxy[k] * tbe_0 * tke_0;

            to_z_x_0_xxz[k] = -4.0 * to_z_xz[k] * tbe_0 + 4.0 * to_z_xxxz[k] * tbe_0 * tke_0;

            to_z_x_0_xyy[k] = -2.0 * to_z_yy[k] * tbe_0 + 4.0 * to_z_xxyy[k] * tbe_0 * tke_0;

            to_z_x_0_xyz[k] = -2.0 * to_z_yz[k] * tbe_0 + 4.0 * to_z_xxyz[k] * tbe_0 * tke_0;

            to_z_x_0_xzz[k] = -2.0 * to_z_zz[k] * tbe_0 + 4.0 * to_z_xxzz[k] * tbe_0 * tke_0;

            to_z_x_0_yyy[k] = 4.0 * to_z_xyyy[k] * tbe_0 * tke_0;

            to_z_x_0_yyz[k] = 4.0 * to_z_xyyz[k] * tbe_0 * tke_0;

            to_z_x_0_yzz[k] = 4.0 * to_z_xyzz[k] * tbe_0 * tke_0;

            to_z_x_0_zzz[k] = 4.0 * to_z_xzzz[k] * tbe_0 * tke_0;
        }

        // Set up 70-80 components of targeted buffer : SF

        auto to_z_y_0_xxx = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 0);

        auto to_z_y_0_xxy = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 1);

        auto to_z_y_0_xxz = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 2);

        auto to_z_y_0_xyy = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 3);

        auto to_z_y_0_xyz = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 4);

        auto to_z_y_0_xzz = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 5);

        auto to_z_y_0_yyy = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 6);

        auto to_z_y_0_yyz = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 7);

        auto to_z_y_0_yzz = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 8);

        auto to_z_y_0_zzz = pbuffer.data(idx_op_geom_101_sf + 7 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_z_xx, to_z_xxxy, to_z_xxyy, to_z_xxyz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_y_0_xxx, to_z_y_0_xxy, to_z_y_0_xxz, to_z_y_0_xyy, to_z_y_0_xyz, to_z_y_0_xzz, to_z_y_0_yyy, to_z_y_0_yyz, to_z_y_0_yzz, to_z_y_0_zzz, to_z_yy, to_z_yyyy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_0_xxx[k] = 4.0 * to_z_xxxy[k] * tbe_0 * tke_0;

            to_z_y_0_xxy[k] = -2.0 * to_z_xx[k] * tbe_0 + 4.0 * to_z_xxyy[k] * tbe_0 * tke_0;

            to_z_y_0_xxz[k] = 4.0 * to_z_xxyz[k] * tbe_0 * tke_0;

            to_z_y_0_xyy[k] = -4.0 * to_z_xy[k] * tbe_0 + 4.0 * to_z_xyyy[k] * tbe_0 * tke_0;

            to_z_y_0_xyz[k] = -2.0 * to_z_xz[k] * tbe_0 + 4.0 * to_z_xyyz[k] * tbe_0 * tke_0;

            to_z_y_0_xzz[k] = 4.0 * to_z_xyzz[k] * tbe_0 * tke_0;

            to_z_y_0_yyy[k] = -6.0 * to_z_yy[k] * tbe_0 + 4.0 * to_z_yyyy[k] * tbe_0 * tke_0;

            to_z_y_0_yyz[k] = -4.0 * to_z_yz[k] * tbe_0 + 4.0 * to_z_yyyz[k] * tbe_0 * tke_0;

            to_z_y_0_yzz[k] = -2.0 * to_z_zz[k] * tbe_0 + 4.0 * to_z_yyzz[k] * tbe_0 * tke_0;

            to_z_y_0_zzz[k] = 4.0 * to_z_yzzz[k] * tbe_0 * tke_0;
        }

        // Set up 80-90 components of targeted buffer : SF

        auto to_z_z_0_xxx = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 0);

        auto to_z_z_0_xxy = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 1);

        auto to_z_z_0_xxz = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 2);

        auto to_z_z_0_xyy = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 3);

        auto to_z_z_0_xyz = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 4);

        auto to_z_z_0_xzz = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 5);

        auto to_z_z_0_yyy = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 6);

        auto to_z_z_0_yyz = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 7);

        auto to_z_z_0_yzz = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 8);

        auto to_z_z_0_zzz = pbuffer.data(idx_op_geom_101_sf + 8 * op_comps * 10 + i * 10 + 9);

        #pragma omp simd aligned(to_z_xx, to_z_xxxz, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_z_0_xxx, to_z_z_0_xxy, to_z_z_0_xxz, to_z_z_0_xyy, to_z_z_0_xyz, to_z_z_0_xzz, to_z_z_0_yyy, to_z_z_0_yyz, to_z_z_0_yzz, to_z_z_0_zzz, to_z_zz, to_z_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_0_xxx[k] = 4.0 * to_z_xxxz[k] * tbe_0 * tke_0;

            to_z_z_0_xxy[k] = 4.0 * to_z_xxyz[k] * tbe_0 * tke_0;

            to_z_z_0_xxz[k] = -2.0 * to_z_xx[k] * tbe_0 + 4.0 * to_z_xxzz[k] * tbe_0 * tke_0;

            to_z_z_0_xyy[k] = 4.0 * to_z_xyyz[k] * tbe_0 * tke_0;

            to_z_z_0_xyz[k] = -2.0 * to_z_xy[k] * tbe_0 + 4.0 * to_z_xyzz[k] * tbe_0 * tke_0;

            to_z_z_0_xzz[k] = -4.0 * to_z_xz[k] * tbe_0 + 4.0 * to_z_xzzz[k] * tbe_0 * tke_0;

            to_z_z_0_yyy[k] = 4.0 * to_z_yyyz[k] * tbe_0 * tke_0;

            to_z_z_0_yyz[k] = -2.0 * to_z_yy[k] * tbe_0 + 4.0 * to_z_yyzz[k] * tbe_0 * tke_0;

            to_z_z_0_yzz[k] = -4.0 * to_z_yz[k] * tbe_0 + 4.0 * to_z_yzzz[k] * tbe_0 * tke_0;

            to_z_z_0_zzz[k] = -6.0 * to_z_zz[k] * tbe_0 + 4.0 * to_z_zzzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

