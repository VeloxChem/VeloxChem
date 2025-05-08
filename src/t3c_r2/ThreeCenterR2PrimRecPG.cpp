#include "ThreeCenterR2PrimRecPG.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_pg(CSimdArray<double>& pbuffer, 
                const size_t idx_g_pg,
                const size_t idx_sf,
                const size_t idx_sg,
                const size_t idx_pd,
                const size_t idx_pf,
                const size_t idx_pg,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_sf);

    auto ts_0_xxy = pbuffer.data(idx_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_sg);

    auto ts_0_xxxy = pbuffer.data(idx_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_pd);

    auto ts_x_xy = pbuffer.data(idx_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_pd + 17);

    // Set up components of auxiliary buffer : PF

    auto ts_x_xxx = pbuffer.data(idx_pf);

    auto ts_x_xxy = pbuffer.data(idx_pf + 1);

    auto ts_x_xxz = pbuffer.data(idx_pf + 2);

    auto ts_x_xyy = pbuffer.data(idx_pf + 3);

    auto ts_x_xyz = pbuffer.data(idx_pf + 4);

    auto ts_x_xzz = pbuffer.data(idx_pf + 5);

    auto ts_x_yyy = pbuffer.data(idx_pf + 6);

    auto ts_x_yyz = pbuffer.data(idx_pf + 7);

    auto ts_x_yzz = pbuffer.data(idx_pf + 8);

    auto ts_x_zzz = pbuffer.data(idx_pf + 9);

    auto ts_y_xxx = pbuffer.data(idx_pf + 10);

    auto ts_y_xxy = pbuffer.data(idx_pf + 11);

    auto ts_y_xxz = pbuffer.data(idx_pf + 12);

    auto ts_y_xyy = pbuffer.data(idx_pf + 13);

    auto ts_y_xyz = pbuffer.data(idx_pf + 14);

    auto ts_y_xzz = pbuffer.data(idx_pf + 15);

    auto ts_y_yyy = pbuffer.data(idx_pf + 16);

    auto ts_y_yyz = pbuffer.data(idx_pf + 17);

    auto ts_y_yzz = pbuffer.data(idx_pf + 18);

    auto ts_y_zzz = pbuffer.data(idx_pf + 19);

    auto ts_z_xxx = pbuffer.data(idx_pf + 20);

    auto ts_z_xxy = pbuffer.data(idx_pf + 21);

    auto ts_z_xxz = pbuffer.data(idx_pf + 22);

    auto ts_z_xyy = pbuffer.data(idx_pf + 23);

    auto ts_z_xyz = pbuffer.data(idx_pf + 24);

    auto ts_z_xzz = pbuffer.data(idx_pf + 25);

    auto ts_z_yyy = pbuffer.data(idx_pf + 26);

    auto ts_z_yyz = pbuffer.data(idx_pf + 27);

    auto ts_z_yzz = pbuffer.data(idx_pf + 28);

    auto ts_z_zzz = pbuffer.data(idx_pf + 29);

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_pg);

    auto ts_x_xxxy = pbuffer.data(idx_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_pg + 44);

    // Set up 0-15 components of targeted buffer : PG

    auto gr_x_xxxx = pbuffer.data(idx_g_pg);

    auto gr_x_xxxy = pbuffer.data(idx_g_pg + 1);

    auto gr_x_xxxz = pbuffer.data(idx_g_pg + 2);

    auto gr_x_xxyy = pbuffer.data(idx_g_pg + 3);

    auto gr_x_xxyz = pbuffer.data(idx_g_pg + 4);

    auto gr_x_xxzz = pbuffer.data(idx_g_pg + 5);

    auto gr_x_xyyy = pbuffer.data(idx_g_pg + 6);

    auto gr_x_xyyz = pbuffer.data(idx_g_pg + 7);

    auto gr_x_xyzz = pbuffer.data(idx_g_pg + 8);

    auto gr_x_xzzz = pbuffer.data(idx_g_pg + 9);

    auto gr_x_yyyy = pbuffer.data(idx_g_pg + 10);

    auto gr_x_yyyz = pbuffer.data(idx_g_pg + 11);

    auto gr_x_yyzz = pbuffer.data(idx_g_pg + 12);

    auto gr_x_yzzz = pbuffer.data(idx_g_pg + 13);

    auto gr_x_zzzz = pbuffer.data(idx_g_pg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxyy, gr_x_xxyz, gr_x_xxzz, gr_x_xyyy, gr_x_xyyz, gr_x_xyzz, gr_x_xzzz, gr_x_yyyy, gr_x_yyyz, gr_x_yyzz, gr_x_yzzz, gr_x_zzzz, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, ts_x_xx, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xy, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xz, ts_x_xzz, ts_x_xzzz, ts_x_yy, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yz, ts_x_yzz, ts_x_yzzz, ts_x_zz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_x_xxxx[i] = 8.0 * ts_0_xxx[i] * gfe_0 + 2.0 * ts_0_xxxx[i] * gfe_0 * gc_x[i] + 12.0 * ts_x_xx[i] * gfe_0 + 8.0 * ts_x_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_x_xxxx[i] * gfe_0 + ts_x_xxxx[i] * rgc2_0;

        gr_x_xxxy[i] = 6.0 * ts_0_xxy[i] * gfe_0 + 2.0 * ts_0_xxxy[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_xy[i] * gfe_0 + 6.0 * ts_x_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_xxxy[i] * gfe_0 + ts_x_xxxy[i] * rgc2_0;

        gr_x_xxxz[i] = 6.0 * ts_0_xxz[i] * gfe_0 + 2.0 * ts_0_xxxz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_xz[i] * gfe_0 + 6.0 * ts_x_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xxxz[i] * gfe_0 + ts_x_xxxz[i] * rgc2_0;

        gr_x_xxyy[i] = 4.0 * ts_0_xyy[i] * gfe_0 + 2.0 * ts_0_xxyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yy[i] * gfe_0 + 4.0 * ts_x_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 + 4.0 * ts_x_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_xxyy[i] * gfe_0 + ts_x_xxyy[i] * rgc2_0;

        gr_x_xxyz[i] = 4.0 * ts_0_xyz[i] * gfe_0 + 2.0 * ts_0_xxyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yz[i] * gfe_0 + 4.0 * ts_x_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xxyz[i] * gfe_0 + ts_x_xxyz[i] * rgc2_0;

        gr_x_xxzz[i] = 4.0 * ts_0_xzz[i] * gfe_0 + 2.0 * ts_0_xxzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe_0 + 4.0 * ts_x_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 + 4.0 * ts_x_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xxzz[i] * gfe_0 + ts_x_xxzz[i] * rgc2_0;

        gr_x_xyyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 + 2.0 * ts_0_xyyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_xy[i] * gfe_0 + 6.0 * ts_x_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_xyyy[i] * gfe_0 + ts_x_xyyy[i] * rgc2_0;

        gr_x_xyyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 + 2.0 * ts_0_xyyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xz[i] * gfe_0 + 4.0 * ts_x_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xyyz[i] * gfe_0 + ts_x_xyyz[i] * rgc2_0;

        gr_x_xyzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 + 2.0 * ts_0_xyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_xy[i] * gfe_0 + 4.0 * ts_x_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xyzz[i] * gfe_0 + ts_x_xyzz[i] * rgc2_0;

        gr_x_xzzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 + 2.0 * ts_0_xzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_xz[i] * gfe_0 + 6.0 * ts_x_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xzzz[i] * gfe_0 + ts_x_xzzz[i] * rgc2_0;

        gr_x_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * gc_x[i] + 12.0 * ts_x_yy[i] * gfe_0 + 8.0 * ts_x_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_yyyy[i] * gfe_0 + ts_x_yyyy[i] * rgc2_0;

        gr_x_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_yz[i] * gfe_0 + 6.0 * ts_x_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_yyyz[i] * gfe_0 + ts_x_yyyz[i] * rgc2_0;

        gr_x_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe_0 + 4.0 * ts_x_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_yy[i] * gfe_0 + 4.0 * ts_x_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_yyzz[i] * gfe_0 + ts_x_yyzz[i] * rgc2_0;

        gr_x_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_x_yz[i] * gfe_0 + 6.0 * ts_x_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_yzzz[i] * gfe_0 + ts_x_yzzz[i] * rgc2_0;

        gr_x_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * gc_x[i] + 12.0 * ts_x_zz[i] * gfe_0 + 8.0 * ts_x_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_zzzz[i] * gfe_0 + ts_x_zzzz[i] * rgc2_0;
    }

    // Set up 15-30 components of targeted buffer : PG

    auto gr_y_xxxx = pbuffer.data(idx_g_pg + 15);

    auto gr_y_xxxy = pbuffer.data(idx_g_pg + 16);

    auto gr_y_xxxz = pbuffer.data(idx_g_pg + 17);

    auto gr_y_xxyy = pbuffer.data(idx_g_pg + 18);

    auto gr_y_xxyz = pbuffer.data(idx_g_pg + 19);

    auto gr_y_xxzz = pbuffer.data(idx_g_pg + 20);

    auto gr_y_xyyy = pbuffer.data(idx_g_pg + 21);

    auto gr_y_xyyz = pbuffer.data(idx_g_pg + 22);

    auto gr_y_xyzz = pbuffer.data(idx_g_pg + 23);

    auto gr_y_xzzz = pbuffer.data(idx_g_pg + 24);

    auto gr_y_yyyy = pbuffer.data(idx_g_pg + 25);

    auto gr_y_yyyz = pbuffer.data(idx_g_pg + 26);

    auto gr_y_yyzz = pbuffer.data(idx_g_pg + 27);

    auto gr_y_yzzz = pbuffer.data(idx_g_pg + 28);

    auto gr_y_zzzz = pbuffer.data(idx_g_pg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxyy, gr_y_xxyz, gr_y_xxzz, gr_y_xyyy, gr_y_xyyz, gr_y_xyzz, gr_y_xzzz, gr_y_yyyy, gr_y_yyyz, gr_y_yyzz, gr_y_yzzz, gr_y_zzzz, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, ts_y_xx, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xy, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xz, ts_y_xzz, ts_y_xzzz, ts_y_yy, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yz, ts_y_yzz, ts_y_yzzz, ts_y_zz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_y_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * gc_y[i] + 12.0 * ts_y_xx[i] * gfe_0 + 8.0 * ts_y_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_y_xxxx[i] * gfe_0 + ts_y_xxxx[i] * rgc2_0;

        gr_y_xxxy[i] = 2.0 * ts_0_xxx[i] * gfe_0 + 2.0 * ts_0_xxxy[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_xy[i] * gfe_0 + 6.0 * ts_y_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_xxxy[i] * gfe_0 + ts_y_xxxy[i] * rgc2_0;

        gr_y_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_xz[i] * gfe_0 + 6.0 * ts_y_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xxxz[i] * gfe_0 + ts_y_xxxz[i] * rgc2_0;

        gr_y_xxyy[i] = 4.0 * ts_0_xxy[i] * gfe_0 + 2.0 * ts_0_xxyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe_0 + 4.0 * ts_y_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xx[i] * gfe_0 + 4.0 * ts_y_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_xxyy[i] * gfe_0 + ts_y_xxyy[i] * rgc2_0;

        gr_y_xxyz[i] = 2.0 * ts_0_xxz[i] * gfe_0 + 2.0 * ts_0_xxyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yz[i] * gfe_0 + 4.0 * ts_y_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xxyz[i] * gfe_0 + ts_y_xxyz[i] * rgc2_0;

        gr_y_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zz[i] * gfe_0 + 4.0 * ts_y_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xx[i] * gfe_0 + 4.0 * ts_y_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xxzz[i] * gfe_0 + ts_y_xxzz[i] * rgc2_0;

        gr_y_xyyy[i] = 6.0 * ts_0_xyy[i] * gfe_0 + 2.0 * ts_0_xyyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_y_xy[i] * gfe_0 + 6.0 * ts_y_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_xyyy[i] * gfe_0 + ts_y_xyyy[i] * rgc2_0;

        gr_y_xyyz[i] = 4.0 * ts_0_xyz[i] * gfe_0 + 2.0 * ts_0_xyyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xz[i] * gfe_0 + 4.0 * ts_y_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xyyz[i] * gfe_0 + ts_y_xyyz[i] * rgc2_0;

        gr_y_xyzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 + 2.0 * ts_0_xyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xy[i] * gfe_0 + 4.0 * ts_y_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xyzz[i] * gfe_0 + ts_y_xyzz[i] * rgc2_0;

        gr_y_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_y_xz[i] * gfe_0 + 6.0 * ts_y_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xzzz[i] * gfe_0 + ts_y_xzzz[i] * rgc2_0;

        gr_y_yyyy[i] = 8.0 * ts_0_yyy[i] * gfe_0 + 2.0 * ts_0_yyyy[i] * gfe_0 * gc_y[i] + 12.0 * ts_y_yy[i] * gfe_0 + 8.0 * ts_y_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_yyyy[i] * gfe_0 + ts_y_yyyy[i] * rgc2_0;

        gr_y_yyyz[i] = 6.0 * ts_0_yyz[i] * gfe_0 + 2.0 * ts_0_yyyz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_yz[i] * gfe_0 + 6.0 * ts_y_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_yyyz[i] * gfe_0 + ts_y_yyyz[i] * rgc2_0;

        gr_y_yyzz[i] = 4.0 * ts_0_yzz[i] * gfe_0 + 2.0 * ts_0_yyzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zz[i] * gfe_0 + 4.0 * ts_y_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe_0 + 4.0 * ts_y_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_yyzz[i] * gfe_0 + ts_y_yyzz[i] * rgc2_0;

        gr_y_yzzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 + 2.0 * ts_0_yzzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_yz[i] * gfe_0 + 6.0 * ts_y_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_yzzz[i] * gfe_0 + ts_y_yzzz[i] * rgc2_0;

        gr_y_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * gc_y[i] + 12.0 * ts_y_zz[i] * gfe_0 + 8.0 * ts_y_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_zzzz[i] * gfe_0 + ts_y_zzzz[i] * rgc2_0;
    }

    // Set up 30-45 components of targeted buffer : PG

    auto gr_z_xxxx = pbuffer.data(idx_g_pg + 30);

    auto gr_z_xxxy = pbuffer.data(idx_g_pg + 31);

    auto gr_z_xxxz = pbuffer.data(idx_g_pg + 32);

    auto gr_z_xxyy = pbuffer.data(idx_g_pg + 33);

    auto gr_z_xxyz = pbuffer.data(idx_g_pg + 34);

    auto gr_z_xxzz = pbuffer.data(idx_g_pg + 35);

    auto gr_z_xyyy = pbuffer.data(idx_g_pg + 36);

    auto gr_z_xyyz = pbuffer.data(idx_g_pg + 37);

    auto gr_z_xyzz = pbuffer.data(idx_g_pg + 38);

    auto gr_z_xzzz = pbuffer.data(idx_g_pg + 39);

    auto gr_z_yyyy = pbuffer.data(idx_g_pg + 40);

    auto gr_z_yyyz = pbuffer.data(idx_g_pg + 41);

    auto gr_z_yyzz = pbuffer.data(idx_g_pg + 42);

    auto gr_z_yzzz = pbuffer.data(idx_g_pg + 43);

    auto gr_z_zzzz = pbuffer.data(idx_g_pg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxyy, gr_z_xxyz, gr_z_xxzz, gr_z_xyyy, gr_z_xyyz, gr_z_xyzz, gr_z_xzzz, gr_z_yyyy, gr_z_yyyz, gr_z_yyzz, gr_z_yzzz, gr_z_zzzz, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, ts_z_xx, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xy, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xz, ts_z_xzz, ts_z_xzzz, ts_z_yy, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yz, ts_z_yzz, ts_z_yzzz, ts_z_zz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_z_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * gc_z[i] + 12.0 * ts_z_xx[i] * gfe_0 + 8.0 * ts_z_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_z_xxxx[i] * gfe_0 + ts_z_xxxx[i] * rgc2_0;

        gr_z_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 * gc_z[i] + 6.0 * ts_z_xy[i] * gfe_0 + 6.0 * ts_z_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_xxxy[i] * gfe_0 + ts_z_xxxy[i] * rgc2_0;

        gr_z_xxxz[i] = 2.0 * ts_0_xxx[i] * gfe_0 + 2.0 * ts_0_xxxz[i] * gfe_0 * gc_z[i] + 6.0 * ts_z_xz[i] * gfe_0 + 6.0 * ts_z_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xxxz[i] * gfe_0 + ts_z_xxxz[i] * rgc2_0;

        gr_z_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yy[i] * gfe_0 + 4.0 * ts_z_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xx[i] * gfe_0 + 4.0 * ts_z_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_xxyy[i] * gfe_0 + ts_z_xxyy[i] * rgc2_0;

        gr_z_xxyz[i] = 2.0 * ts_0_xxy[i] * gfe_0 + 2.0 * ts_0_xxyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yz[i] * gfe_0 + 4.0 * ts_z_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xxyz[i] * gfe_0 + ts_z_xxyz[i] * rgc2_0;

        gr_z_xxzz[i] = 4.0 * ts_0_xxz[i] * gfe_0 + 2.0 * ts_0_xxzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_zz[i] * gfe_0 + 4.0 * ts_z_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xx[i] * gfe_0 + 4.0 * ts_z_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xxzz[i] * gfe_0 + ts_z_xxzz[i] * rgc2_0;

        gr_z_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_z_xy[i] * gfe_0 + 6.0 * ts_z_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_xyyy[i] * gfe_0 + ts_z_xyyy[i] * rgc2_0;

        gr_z_xyyz[i] = 2.0 * ts_0_xyy[i] * gfe_0 + 2.0 * ts_0_xyyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xz[i] * gfe_0 + 4.0 * ts_z_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xyyz[i] * gfe_0 + ts_z_xyyz[i] * rgc2_0;

        gr_z_xyzz[i] = 4.0 * ts_0_xyz[i] * gfe_0 + 2.0 * ts_0_xyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_xy[i] * gfe_0 + 4.0 * ts_z_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xyzz[i] * gfe_0 + ts_z_xyzz[i] * rgc2_0;

        gr_z_xzzz[i] = 6.0 * ts_0_xzz[i] * gfe_0 + 2.0 * ts_0_xzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_z_xz[i] * gfe_0 + 6.0 * ts_z_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xzzz[i] * gfe_0 + ts_z_xzzz[i] * rgc2_0;

        gr_z_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * gc_z[i] + 12.0 * ts_z_yy[i] * gfe_0 + 8.0 * ts_z_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_yyyy[i] * gfe_0 + ts_z_yyyy[i] * rgc2_0;

        gr_z_yyyz[i] = 2.0 * ts_0_yyy[i] * gfe_0 + 2.0 * ts_0_yyyz[i] * gfe_0 * gc_z[i] + 6.0 * ts_z_yz[i] * gfe_0 + 6.0 * ts_z_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_yyyz[i] * gfe_0 + ts_z_yyyz[i] * rgc2_0;

        gr_z_yyzz[i] = 4.0 * ts_0_yyz[i] * gfe_0 + 2.0 * ts_0_yyzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_zz[i] * gfe_0 + 4.0 * ts_z_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_yy[i] * gfe_0 + 4.0 * ts_z_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_yyzz[i] * gfe_0 + ts_z_yyzz[i] * rgc2_0;

        gr_z_yzzz[i] = 6.0 * ts_0_yzz[i] * gfe_0 + 2.0 * ts_0_yzzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_z_yz[i] * gfe_0 + 6.0 * ts_z_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_yzzz[i] * gfe_0 + ts_z_yzzz[i] * rgc2_0;

        gr_z_zzzz[i] = 8.0 * ts_0_zzz[i] * gfe_0 + 2.0 * ts_0_zzzz[i] * gfe_0 * gc_z[i] + 12.0 * ts_z_zz[i] * gfe_0 + 8.0 * ts_z_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_zzzz[i] * gfe_0 + ts_z_zzzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

