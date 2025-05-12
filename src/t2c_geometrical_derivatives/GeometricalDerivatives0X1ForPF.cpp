#include "GeometricalDerivatives0X1ForPF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_geom_deriv_0x1_pf(CSimdArray<double>& pbuffer,
                       const int idx_op_geom_001_pf,
                       const int idx_op_pd,
                       const int idx_op_pg,
                       const CSimdArray<double>& factors,
                       const int op_comps) -> void
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

        // Set up 0-10 components of targeted buffer : PF

        auto to_0_x_x_xxx = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 0);

        auto to_0_x_x_xxy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 1);

        auto to_0_x_x_xxz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 2);

        auto to_0_x_x_xyy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 3);

        auto to_0_x_x_xyz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 4);

        auto to_0_x_x_xzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 5);

        auto to_0_x_x_yyy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 6);

        auto to_0_x_x_yyz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 7);

        auto to_0_x_x_yzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 8);

        auto to_0_x_x_zzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_0_x_x_xxx, to_0_x_x_xxy, to_0_x_x_xxz, to_0_x_x_xyy, to_0_x_x_xyz, to_0_x_x_xzz, to_0_x_x_yyy, to_0_x_x_yyz, to_0_x_x_yzz, to_0_x_x_zzz, to_x_xx, to_x_xxxx, to_x_xxxy, to_x_xxxz, to_x_xxyy, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yz, to_x_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_x_xxx[k] = -3.0 * to_x_xx[k] + 2.0 * to_x_xxxx[k] * tke_0;

            to_0_x_x_xxy[k] = -2.0 * to_x_xy[k] + 2.0 * to_x_xxxy[k] * tke_0;

            to_0_x_x_xxz[k] = -2.0 * to_x_xz[k] + 2.0 * to_x_xxxz[k] * tke_0;

            to_0_x_x_xyy[k] = -to_x_yy[k] + 2.0 * to_x_xxyy[k] * tke_0;

            to_0_x_x_xyz[k] = -to_x_yz[k] + 2.0 * to_x_xxyz[k] * tke_0;

            to_0_x_x_xzz[k] = -to_x_zz[k] + 2.0 * to_x_xxzz[k] * tke_0;

            to_0_x_x_yyy[k] = 2.0 * to_x_xyyy[k] * tke_0;

            to_0_x_x_yyz[k] = 2.0 * to_x_xyyz[k] * tke_0;

            to_0_x_x_yzz[k] = 2.0 * to_x_xyzz[k] * tke_0;

            to_0_x_x_zzz[k] = 2.0 * to_x_xzzz[k] * tke_0;
        }

        // Set up 10-20 components of targeted buffer : PF

        auto to_0_x_y_xxx = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 10);

        auto to_0_x_y_xxy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 11);

        auto to_0_x_y_xxz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 12);

        auto to_0_x_y_xyy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 13);

        auto to_0_x_y_xyz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 14);

        auto to_0_x_y_xzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 15);

        auto to_0_x_y_yyy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 16);

        auto to_0_x_y_yyz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 17);

        auto to_0_x_y_yzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 18);

        auto to_0_x_y_zzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_0_x_y_xxx, to_0_x_y_xxy, to_0_x_y_xxz, to_0_x_y_xyy, to_0_x_y_xyz, to_0_x_y_xzz, to_0_x_y_yyy, to_0_x_y_yyz, to_0_x_y_yzz, to_0_x_y_zzz, to_y_xx, to_y_xxxx, to_y_xxxy, to_y_xxxz, to_y_xxyy, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_y_xxx[k] = -3.0 * to_y_xx[k] + 2.0 * to_y_xxxx[k] * tke_0;

            to_0_x_y_xxy[k] = -2.0 * to_y_xy[k] + 2.0 * to_y_xxxy[k] * tke_0;

            to_0_x_y_xxz[k] = -2.0 * to_y_xz[k] + 2.0 * to_y_xxxz[k] * tke_0;

            to_0_x_y_xyy[k] = -to_y_yy[k] + 2.0 * to_y_xxyy[k] * tke_0;

            to_0_x_y_xyz[k] = -to_y_yz[k] + 2.0 * to_y_xxyz[k] * tke_0;

            to_0_x_y_xzz[k] = -to_y_zz[k] + 2.0 * to_y_xxzz[k] * tke_0;

            to_0_x_y_yyy[k] = 2.0 * to_y_xyyy[k] * tke_0;

            to_0_x_y_yyz[k] = 2.0 * to_y_xyyz[k] * tke_0;

            to_0_x_y_yzz[k] = 2.0 * to_y_xyzz[k] * tke_0;

            to_0_x_y_zzz[k] = 2.0 * to_y_xzzz[k] * tke_0;
        }

        // Set up 20-30 components of targeted buffer : PF

        auto to_0_x_z_xxx = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 20);

        auto to_0_x_z_xxy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 21);

        auto to_0_x_z_xxz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 22);

        auto to_0_x_z_xyy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 23);

        auto to_0_x_z_xyz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 24);

        auto to_0_x_z_xzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 25);

        auto to_0_x_z_yyy = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 26);

        auto to_0_x_z_yyz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 27);

        auto to_0_x_z_yzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 28);

        auto to_0_x_z_zzz = pbuffer.data(idx_op_geom_001_pf + 0 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_x_z_xxx, to_0_x_z_xxy, to_0_x_z_xxz, to_0_x_z_xyy, to_0_x_z_xyz, to_0_x_z_xzz, to_0_x_z_yyy, to_0_x_z_yyz, to_0_x_z_yzz, to_0_x_z_zzz, to_z_xx, to_z_xxxx, to_z_xxxy, to_z_xxxz, to_z_xxyy, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_x_z_xxx[k] = -3.0 * to_z_xx[k] + 2.0 * to_z_xxxx[k] * tke_0;

            to_0_x_z_xxy[k] = -2.0 * to_z_xy[k] + 2.0 * to_z_xxxy[k] * tke_0;

            to_0_x_z_xxz[k] = -2.0 * to_z_xz[k] + 2.0 * to_z_xxxz[k] * tke_0;

            to_0_x_z_xyy[k] = -to_z_yy[k] + 2.0 * to_z_xxyy[k] * tke_0;

            to_0_x_z_xyz[k] = -to_z_yz[k] + 2.0 * to_z_xxyz[k] * tke_0;

            to_0_x_z_xzz[k] = -to_z_zz[k] + 2.0 * to_z_xxzz[k] * tke_0;

            to_0_x_z_yyy[k] = 2.0 * to_z_xyyy[k] * tke_0;

            to_0_x_z_yyz[k] = 2.0 * to_z_xyyz[k] * tke_0;

            to_0_x_z_yzz[k] = 2.0 * to_z_xyzz[k] * tke_0;

            to_0_x_z_zzz[k] = 2.0 * to_z_xzzz[k] * tke_0;
        }

        // Set up 30-40 components of targeted buffer : PF

        auto to_0_y_x_xxx = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 0);

        auto to_0_y_x_xxy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 1);

        auto to_0_y_x_xxz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 2);

        auto to_0_y_x_xyy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 3);

        auto to_0_y_x_xyz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 4);

        auto to_0_y_x_xzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 5);

        auto to_0_y_x_yyy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 6);

        auto to_0_y_x_yyz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 7);

        auto to_0_y_x_yzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 8);

        auto to_0_y_x_zzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_0_y_x_xxx, to_0_y_x_xxy, to_0_y_x_xxz, to_0_y_x_xyy, to_0_y_x_xyz, to_0_y_x_xzz, to_0_y_x_yyy, to_0_y_x_yyz, to_0_y_x_yzz, to_0_y_x_zzz, to_x_xx, to_x_xxxy, to_x_xxyy, to_x_xxyz, to_x_xy, to_x_xyyy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_yy, to_x_yyyy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_x_xxx[k] = 2.0 * to_x_xxxy[k] * tke_0;

            to_0_y_x_xxy[k] = -to_x_xx[k] + 2.0 * to_x_xxyy[k] * tke_0;

            to_0_y_x_xxz[k] = 2.0 * to_x_xxyz[k] * tke_0;

            to_0_y_x_xyy[k] = -2.0 * to_x_xy[k] + 2.0 * to_x_xyyy[k] * tke_0;

            to_0_y_x_xyz[k] = -to_x_xz[k] + 2.0 * to_x_xyyz[k] * tke_0;

            to_0_y_x_xzz[k] = 2.0 * to_x_xyzz[k] * tke_0;

            to_0_y_x_yyy[k] = -3.0 * to_x_yy[k] + 2.0 * to_x_yyyy[k] * tke_0;

            to_0_y_x_yyz[k] = -2.0 * to_x_yz[k] + 2.0 * to_x_yyyz[k] * tke_0;

            to_0_y_x_yzz[k] = -to_x_zz[k] + 2.0 * to_x_yyzz[k] * tke_0;

            to_0_y_x_zzz[k] = 2.0 * to_x_yzzz[k] * tke_0;
        }

        // Set up 40-50 components of targeted buffer : PF

        auto to_0_y_y_xxx = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 10);

        auto to_0_y_y_xxy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 11);

        auto to_0_y_y_xxz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 12);

        auto to_0_y_y_xyy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 13);

        auto to_0_y_y_xyz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 14);

        auto to_0_y_y_xzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 15);

        auto to_0_y_y_yyy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 16);

        auto to_0_y_y_yyz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 17);

        auto to_0_y_y_yzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 18);

        auto to_0_y_y_zzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_0_y_y_xxx, to_0_y_y_xxy, to_0_y_y_xxz, to_0_y_y_xyy, to_0_y_y_xyz, to_0_y_y_xzz, to_0_y_y_yyy, to_0_y_y_yyz, to_0_y_y_yzz, to_0_y_y_zzz, to_y_xx, to_y_xxxy, to_y_xxyy, to_y_xxyz, to_y_xy, to_y_xyyy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_yy, to_y_yyyy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_y_xxx[k] = 2.0 * to_y_xxxy[k] * tke_0;

            to_0_y_y_xxy[k] = -to_y_xx[k] + 2.0 * to_y_xxyy[k] * tke_0;

            to_0_y_y_xxz[k] = 2.0 * to_y_xxyz[k] * tke_0;

            to_0_y_y_xyy[k] = -2.0 * to_y_xy[k] + 2.0 * to_y_xyyy[k] * tke_0;

            to_0_y_y_xyz[k] = -to_y_xz[k] + 2.0 * to_y_xyyz[k] * tke_0;

            to_0_y_y_xzz[k] = 2.0 * to_y_xyzz[k] * tke_0;

            to_0_y_y_yyy[k] = -3.0 * to_y_yy[k] + 2.0 * to_y_yyyy[k] * tke_0;

            to_0_y_y_yyz[k] = -2.0 * to_y_yz[k] + 2.0 * to_y_yyyz[k] * tke_0;

            to_0_y_y_yzz[k] = -to_y_zz[k] + 2.0 * to_y_yyzz[k] * tke_0;

            to_0_y_y_zzz[k] = 2.0 * to_y_yzzz[k] * tke_0;
        }

        // Set up 50-60 components of targeted buffer : PF

        auto to_0_y_z_xxx = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 20);

        auto to_0_y_z_xxy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 21);

        auto to_0_y_z_xxz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 22);

        auto to_0_y_z_xyy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 23);

        auto to_0_y_z_xyz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 24);

        auto to_0_y_z_xzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 25);

        auto to_0_y_z_yyy = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 26);

        auto to_0_y_z_yyz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 27);

        auto to_0_y_z_yzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 28);

        auto to_0_y_z_zzz = pbuffer.data(idx_op_geom_001_pf + 1 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_y_z_xxx, to_0_y_z_xxy, to_0_y_z_xxz, to_0_y_z_xyy, to_0_y_z_xyz, to_0_y_z_xzz, to_0_y_z_yyy, to_0_y_z_yyz, to_0_y_z_yzz, to_0_y_z_zzz, to_z_xx, to_z_xxxy, to_z_xxyy, to_z_xxyz, to_z_xy, to_z_xyyy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_yy, to_z_yyyy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_y_z_xxx[k] = 2.0 * to_z_xxxy[k] * tke_0;

            to_0_y_z_xxy[k] = -to_z_xx[k] + 2.0 * to_z_xxyy[k] * tke_0;

            to_0_y_z_xxz[k] = 2.0 * to_z_xxyz[k] * tke_0;

            to_0_y_z_xyy[k] = -2.0 * to_z_xy[k] + 2.0 * to_z_xyyy[k] * tke_0;

            to_0_y_z_xyz[k] = -to_z_xz[k] + 2.0 * to_z_xyyz[k] * tke_0;

            to_0_y_z_xzz[k] = 2.0 * to_z_xyzz[k] * tke_0;

            to_0_y_z_yyy[k] = -3.0 * to_z_yy[k] + 2.0 * to_z_yyyy[k] * tke_0;

            to_0_y_z_yyz[k] = -2.0 * to_z_yz[k] + 2.0 * to_z_yyyz[k] * tke_0;

            to_0_y_z_yzz[k] = -to_z_zz[k] + 2.0 * to_z_yyzz[k] * tke_0;

            to_0_y_z_zzz[k] = 2.0 * to_z_yzzz[k] * tke_0;
        }

        // Set up 60-70 components of targeted buffer : PF

        auto to_0_z_x_xxx = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 0);

        auto to_0_z_x_xxy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 1);

        auto to_0_z_x_xxz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 2);

        auto to_0_z_x_xyy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 3);

        auto to_0_z_x_xyz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 4);

        auto to_0_z_x_xzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 5);

        auto to_0_z_x_yyy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 6);

        auto to_0_z_x_yyz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 7);

        auto to_0_z_x_yzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 8);

        auto to_0_z_x_zzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 9);

        #pragma omp simd aligned(to_0_z_x_xxx, to_0_z_x_xxy, to_0_z_x_xxz, to_0_z_x_xyy, to_0_z_x_xyz, to_0_z_x_xzz, to_0_z_x_yyy, to_0_z_x_yyz, to_0_z_x_yzz, to_0_z_x_zzz, to_x_xx, to_x_xxxz, to_x_xxyz, to_x_xxzz, to_x_xy, to_x_xyyz, to_x_xyzz, to_x_xz, to_x_xzzz, to_x_yy, to_x_yyyz, to_x_yyzz, to_x_yz, to_x_yzzz, to_x_zz, to_x_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_x_xxx[k] = 2.0 * to_x_xxxz[k] * tke_0;

            to_0_z_x_xxy[k] = 2.0 * to_x_xxyz[k] * tke_0;

            to_0_z_x_xxz[k] = -to_x_xx[k] + 2.0 * to_x_xxzz[k] * tke_0;

            to_0_z_x_xyy[k] = 2.0 * to_x_xyyz[k] * tke_0;

            to_0_z_x_xyz[k] = -to_x_xy[k] + 2.0 * to_x_xyzz[k] * tke_0;

            to_0_z_x_xzz[k] = -2.0 * to_x_xz[k] + 2.0 * to_x_xzzz[k] * tke_0;

            to_0_z_x_yyy[k] = 2.0 * to_x_yyyz[k] * tke_0;

            to_0_z_x_yyz[k] = -to_x_yy[k] + 2.0 * to_x_yyzz[k] * tke_0;

            to_0_z_x_yzz[k] = -2.0 * to_x_yz[k] + 2.0 * to_x_yzzz[k] * tke_0;

            to_0_z_x_zzz[k] = -3.0 * to_x_zz[k] + 2.0 * to_x_zzzz[k] * tke_0;
        }

        // Set up 70-80 components of targeted buffer : PF

        auto to_0_z_y_xxx = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 10);

        auto to_0_z_y_xxy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 11);

        auto to_0_z_y_xxz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 12);

        auto to_0_z_y_xyy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 13);

        auto to_0_z_y_xyz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 14);

        auto to_0_z_y_xzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 15);

        auto to_0_z_y_yyy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 16);

        auto to_0_z_y_yyz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 17);

        auto to_0_z_y_yzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 18);

        auto to_0_z_y_zzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 19);

        #pragma omp simd aligned(to_0_z_y_xxx, to_0_z_y_xxy, to_0_z_y_xxz, to_0_z_y_xyy, to_0_z_y_xyz, to_0_z_y_xzz, to_0_z_y_yyy, to_0_z_y_yyz, to_0_z_y_yzz, to_0_z_y_zzz, to_y_xx, to_y_xxxz, to_y_xxyz, to_y_xxzz, to_y_xy, to_y_xyyz, to_y_xyzz, to_y_xz, to_y_xzzz, to_y_yy, to_y_yyyz, to_y_yyzz, to_y_yz, to_y_yzzz, to_y_zz, to_y_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_y_xxx[k] = 2.0 * to_y_xxxz[k] * tke_0;

            to_0_z_y_xxy[k] = 2.0 * to_y_xxyz[k] * tke_0;

            to_0_z_y_xxz[k] = -to_y_xx[k] + 2.0 * to_y_xxzz[k] * tke_0;

            to_0_z_y_xyy[k] = 2.0 * to_y_xyyz[k] * tke_0;

            to_0_z_y_xyz[k] = -to_y_xy[k] + 2.0 * to_y_xyzz[k] * tke_0;

            to_0_z_y_xzz[k] = -2.0 * to_y_xz[k] + 2.0 * to_y_xzzz[k] * tke_0;

            to_0_z_y_yyy[k] = 2.0 * to_y_yyyz[k] * tke_0;

            to_0_z_y_yyz[k] = -to_y_yy[k] + 2.0 * to_y_yyzz[k] * tke_0;

            to_0_z_y_yzz[k] = -2.0 * to_y_yz[k] + 2.0 * to_y_yzzz[k] * tke_0;

            to_0_z_y_zzz[k] = -3.0 * to_y_zz[k] + 2.0 * to_y_zzzz[k] * tke_0;
        }

        // Set up 80-90 components of targeted buffer : PF

        auto to_0_z_z_xxx = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 20);

        auto to_0_z_z_xxy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 21);

        auto to_0_z_z_xxz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 22);

        auto to_0_z_z_xyy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 23);

        auto to_0_z_z_xyz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 24);

        auto to_0_z_z_xzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 25);

        auto to_0_z_z_yyy = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 26);

        auto to_0_z_z_yyz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 27);

        auto to_0_z_z_yzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 28);

        auto to_0_z_z_zzz = pbuffer.data(idx_op_geom_001_pf + 2 * op_comps * 30 + i * 30 + 29);

        #pragma omp simd aligned(to_0_z_z_xxx, to_0_z_z_xxy, to_0_z_z_xxz, to_0_z_z_xyy, to_0_z_z_xyz, to_0_z_z_xzz, to_0_z_z_yyy, to_0_z_z_yyz, to_0_z_z_yzz, to_0_z_z_zzz, to_z_xx, to_z_xxxz, to_z_xxyz, to_z_xxzz, to_z_xy, to_z_xyyz, to_z_xyzz, to_z_xz, to_z_xzzz, to_z_yy, to_z_yyyz, to_z_yyzz, to_z_yz, to_z_yzzz, to_z_zz, to_z_zzzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tke_0 = b_exps[k];

            to_0_z_z_xxx[k] = 2.0 * to_z_xxxz[k] * tke_0;

            to_0_z_z_xxy[k] = 2.0 * to_z_xxyz[k] * tke_0;

            to_0_z_z_xxz[k] = -to_z_xx[k] + 2.0 * to_z_xxzz[k] * tke_0;

            to_0_z_z_xyy[k] = 2.0 * to_z_xyyz[k] * tke_0;

            to_0_z_z_xyz[k] = -to_z_xy[k] + 2.0 * to_z_xyzz[k] * tke_0;

            to_0_z_z_xzz[k] = -2.0 * to_z_xz[k] + 2.0 * to_z_xzzz[k] * tke_0;

            to_0_z_z_yyy[k] = 2.0 * to_z_yyyz[k] * tke_0;

            to_0_z_z_yyz[k] = -to_z_yy[k] + 2.0 * to_z_yyzz[k] * tke_0;

            to_0_z_z_yzz[k] = -2.0 * to_z_yz[k] + 2.0 * to_z_yzzz[k] * tke_0;

            to_0_z_z_zzz[k] = -3.0 * to_z_zz[k] + 2.0 * to_z_zzzz[k] * tke_0;
        }

    }

}

} // t2cgeom namespace

