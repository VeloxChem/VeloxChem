#include "ThreeCenterOverlapGradientPrimRecPF.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_pf(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_pf,
                              const size_t idx_sf,
                              const size_t idx_pd,
                              const size_t idx_pf,
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

    // Set up 0-10 components of targeted buffer : PF

    auto gs_x_x_xxx = pbuffer.data(idx_g_pf);

    auto gs_x_x_xxy = pbuffer.data(idx_g_pf + 1);

    auto gs_x_x_xxz = pbuffer.data(idx_g_pf + 2);

    auto gs_x_x_xyy = pbuffer.data(idx_g_pf + 3);

    auto gs_x_x_xyz = pbuffer.data(idx_g_pf + 4);

    auto gs_x_x_xzz = pbuffer.data(idx_g_pf + 5);

    auto gs_x_x_yyy = pbuffer.data(idx_g_pf + 6);

    auto gs_x_x_yyz = pbuffer.data(idx_g_pf + 7);

    auto gs_x_x_yzz = pbuffer.data(idx_g_pf + 8);

    auto gs_x_x_zzz = pbuffer.data(idx_g_pf + 9);

    #pragma omp simd aligned(gc_x, gs_x_x_xxx, gs_x_x_xxy, gs_x_x_xxz, gs_x_x_xyy, gs_x_x_xyz, gs_x_x_xzz, gs_x_x_yyy, gs_x_x_yyz, gs_x_x_yzz, gs_x_x_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_x_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * tce_0 + 6.0 * ts_x_xx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxx[i] * gc_x[i] * tce_0;

        gs_x_x_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 * tce_0 + 4.0 * ts_x_xy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxy[i] * gc_x[i] * tce_0;

        gs_x_x_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 * tce_0 + 4.0 * ts_x_xz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxz[i] * gc_x[i] * tce_0;

        gs_x_x_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyy[i] * gc_x[i] * tce_0;

        gs_x_x_xyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyz[i] * gc_x[i] * tce_0;

        gs_x_x_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xzz[i] * gc_x[i] * tce_0;

        gs_x_x_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyy[i] * gc_x[i] * tce_0;

        gs_x_x_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyz[i] * gc_x[i] * tce_0;

        gs_x_x_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzz[i] * gc_x[i] * tce_0;

        gs_x_x_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 10-20 components of targeted buffer : PF

    auto gs_x_y_xxx = pbuffer.data(idx_g_pf + 10);

    auto gs_x_y_xxy = pbuffer.data(idx_g_pf + 11);

    auto gs_x_y_xxz = pbuffer.data(idx_g_pf + 12);

    auto gs_x_y_xyy = pbuffer.data(idx_g_pf + 13);

    auto gs_x_y_xyz = pbuffer.data(idx_g_pf + 14);

    auto gs_x_y_xzz = pbuffer.data(idx_g_pf + 15);

    auto gs_x_y_yyy = pbuffer.data(idx_g_pf + 16);

    auto gs_x_y_yyz = pbuffer.data(idx_g_pf + 17);

    auto gs_x_y_yzz = pbuffer.data(idx_g_pf + 18);

    auto gs_x_y_zzz = pbuffer.data(idx_g_pf + 19);

    #pragma omp simd aligned(gc_x, gs_x_y_xxx, gs_x_y_xxy, gs_x_y_xxz, gs_x_y_xyy, gs_x_y_xyz, gs_x_y_xzz, gs_x_y_yyy, gs_x_y_yyz, gs_x_y_yzz, gs_x_y_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_y_xxx[i] = 6.0 * ts_y_xx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxx[i] * gc_x[i] * tce_0;

        gs_x_y_xxy[i] = 4.0 * ts_y_xy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxy[i] * gc_x[i] * tce_0;

        gs_x_y_xxz[i] = 4.0 * ts_y_xz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxz[i] * gc_x[i] * tce_0;

        gs_x_y_xyy[i] = 2.0 * ts_y_yy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyy[i] * gc_x[i] * tce_0;

        gs_x_y_xyz[i] = 2.0 * ts_y_yz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyz[i] * gc_x[i] * tce_0;

        gs_x_y_xzz[i] = 2.0 * ts_y_zz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzz[i] * gc_x[i] * tce_0;

        gs_x_y_yyy[i] = 2.0 * ts_y_yyy[i] * gc_x[i] * tce_0;

        gs_x_y_yyz[i] = 2.0 * ts_y_yyz[i] * gc_x[i] * tce_0;

        gs_x_y_yzz[i] = 2.0 * ts_y_yzz[i] * gc_x[i] * tce_0;

        gs_x_y_zzz[i] = 2.0 * ts_y_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 20-30 components of targeted buffer : PF

    auto gs_x_z_xxx = pbuffer.data(idx_g_pf + 20);

    auto gs_x_z_xxy = pbuffer.data(idx_g_pf + 21);

    auto gs_x_z_xxz = pbuffer.data(idx_g_pf + 22);

    auto gs_x_z_xyy = pbuffer.data(idx_g_pf + 23);

    auto gs_x_z_xyz = pbuffer.data(idx_g_pf + 24);

    auto gs_x_z_xzz = pbuffer.data(idx_g_pf + 25);

    auto gs_x_z_yyy = pbuffer.data(idx_g_pf + 26);

    auto gs_x_z_yyz = pbuffer.data(idx_g_pf + 27);

    auto gs_x_z_yzz = pbuffer.data(idx_g_pf + 28);

    auto gs_x_z_zzz = pbuffer.data(idx_g_pf + 29);

    #pragma omp simd aligned(gc_x, gs_x_z_xxx, gs_x_z_xxy, gs_x_z_xxz, gs_x_z_xyy, gs_x_z_xyz, gs_x_z_xzz, gs_x_z_yyy, gs_x_z_yyz, gs_x_z_yzz, gs_x_z_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_z_xxx[i] = 6.0 * ts_z_xx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxx[i] * gc_x[i] * tce_0;

        gs_x_z_xxy[i] = 4.0 * ts_z_xy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxy[i] * gc_x[i] * tce_0;

        gs_x_z_xxz[i] = 4.0 * ts_z_xz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxz[i] * gc_x[i] * tce_0;

        gs_x_z_xyy[i] = 2.0 * ts_z_yy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyy[i] * gc_x[i] * tce_0;

        gs_x_z_xyz[i] = 2.0 * ts_z_yz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyz[i] * gc_x[i] * tce_0;

        gs_x_z_xzz[i] = 2.0 * ts_z_zz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xzz[i] * gc_x[i] * tce_0;

        gs_x_z_yyy[i] = 2.0 * ts_z_yyy[i] * gc_x[i] * tce_0;

        gs_x_z_yyz[i] = 2.0 * ts_z_yyz[i] * gc_x[i] * tce_0;

        gs_x_z_yzz[i] = 2.0 * ts_z_yzz[i] * gc_x[i] * tce_0;

        gs_x_z_zzz[i] = 2.0 * ts_z_zzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-40 components of targeted buffer : PF

    auto gs_y_x_xxx = pbuffer.data(idx_g_pf + 30);

    auto gs_y_x_xxy = pbuffer.data(idx_g_pf + 31);

    auto gs_y_x_xxz = pbuffer.data(idx_g_pf + 32);

    auto gs_y_x_xyy = pbuffer.data(idx_g_pf + 33);

    auto gs_y_x_xyz = pbuffer.data(idx_g_pf + 34);

    auto gs_y_x_xzz = pbuffer.data(idx_g_pf + 35);

    auto gs_y_x_yyy = pbuffer.data(idx_g_pf + 36);

    auto gs_y_x_yyz = pbuffer.data(idx_g_pf + 37);

    auto gs_y_x_yzz = pbuffer.data(idx_g_pf + 38);

    auto gs_y_x_zzz = pbuffer.data(idx_g_pf + 39);

    #pragma omp simd aligned(gc_y, gs_y_x_xxx, gs_y_x_xxy, gs_y_x_xxz, gs_y_x_xyy, gs_y_x_xyz, gs_y_x_xzz, gs_y_x_yyy, gs_y_x_yyz, gs_y_x_yzz, gs_y_x_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_x_xxx[i] = 2.0 * ts_x_xxx[i] * gc_y[i] * tce_0;

        gs_y_x_xxy[i] = 2.0 * ts_x_xx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxy[i] * gc_y[i] * tce_0;

        gs_y_x_xxz[i] = 2.0 * ts_x_xxz[i] * gc_y[i] * tce_0;

        gs_y_x_xyy[i] = 4.0 * ts_x_xy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyy[i] * gc_y[i] * tce_0;

        gs_y_x_xyz[i] = 2.0 * ts_x_xz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyz[i] * gc_y[i] * tce_0;

        gs_y_x_xzz[i] = 2.0 * ts_x_xzz[i] * gc_y[i] * tce_0;

        gs_y_x_yyy[i] = 6.0 * ts_x_yy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyy[i] * gc_y[i] * tce_0;

        gs_y_x_yyz[i] = 4.0 * ts_x_yz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyz[i] * gc_y[i] * tce_0;

        gs_y_x_yzz[i] = 2.0 * ts_x_zz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzz[i] * gc_y[i] * tce_0;

        gs_y_x_zzz[i] = 2.0 * ts_x_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 40-50 components of targeted buffer : PF

    auto gs_y_y_xxx = pbuffer.data(idx_g_pf + 40);

    auto gs_y_y_xxy = pbuffer.data(idx_g_pf + 41);

    auto gs_y_y_xxz = pbuffer.data(idx_g_pf + 42);

    auto gs_y_y_xyy = pbuffer.data(idx_g_pf + 43);

    auto gs_y_y_xyz = pbuffer.data(idx_g_pf + 44);

    auto gs_y_y_xzz = pbuffer.data(idx_g_pf + 45);

    auto gs_y_y_yyy = pbuffer.data(idx_g_pf + 46);

    auto gs_y_y_yyz = pbuffer.data(idx_g_pf + 47);

    auto gs_y_y_yzz = pbuffer.data(idx_g_pf + 48);

    auto gs_y_y_zzz = pbuffer.data(idx_g_pf + 49);

    #pragma omp simd aligned(gc_y, gs_y_y_xxx, gs_y_y_xxy, gs_y_y_xxz, gs_y_y_xyy, gs_y_y_xyz, gs_y_y_xzz, gs_y_y_yyy, gs_y_y_yyz, gs_y_y_yzz, gs_y_y_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_y_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxx[i] * gc_y[i] * tce_0;

        gs_y_y_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxy[i] * gc_y[i] * tce_0;

        gs_y_y_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxz[i] * gc_y[i] * tce_0;

        gs_y_y_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 * tce_0 + 4.0 * ts_y_xy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyy[i] * gc_y[i] * tce_0;

        gs_y_y_xyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyz[i] * gc_y[i] * tce_0;

        gs_y_y_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzz[i] * gc_y[i] * tce_0;

        gs_y_y_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * tce_0 + 6.0 * ts_y_yy[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyy[i] * gc_y[i] * tce_0;

        gs_y_y_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * tce_0 + 4.0 * ts_y_yz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyz[i] * gc_y[i] * tce_0;

        gs_y_y_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yzz[i] * gc_y[i] * tce_0;

        gs_y_y_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 50-60 components of targeted buffer : PF

    auto gs_y_z_xxx = pbuffer.data(idx_g_pf + 50);

    auto gs_y_z_xxy = pbuffer.data(idx_g_pf + 51);

    auto gs_y_z_xxz = pbuffer.data(idx_g_pf + 52);

    auto gs_y_z_xyy = pbuffer.data(idx_g_pf + 53);

    auto gs_y_z_xyz = pbuffer.data(idx_g_pf + 54);

    auto gs_y_z_xzz = pbuffer.data(idx_g_pf + 55);

    auto gs_y_z_yyy = pbuffer.data(idx_g_pf + 56);

    auto gs_y_z_yyz = pbuffer.data(idx_g_pf + 57);

    auto gs_y_z_yzz = pbuffer.data(idx_g_pf + 58);

    auto gs_y_z_zzz = pbuffer.data(idx_g_pf + 59);

    #pragma omp simd aligned(gc_y, gs_y_z_xxx, gs_y_z_xxy, gs_y_z_xxz, gs_y_z_xyy, gs_y_z_xyz, gs_y_z_xzz, gs_y_z_yyy, gs_y_z_yyz, gs_y_z_yzz, gs_y_z_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_z_xxx[i] = 2.0 * ts_z_xxx[i] * gc_y[i] * tce_0;

        gs_y_z_xxy[i] = 2.0 * ts_z_xx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxy[i] * gc_y[i] * tce_0;

        gs_y_z_xxz[i] = 2.0 * ts_z_xxz[i] * gc_y[i] * tce_0;

        gs_y_z_xyy[i] = 4.0 * ts_z_xy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyy[i] * gc_y[i] * tce_0;

        gs_y_z_xyz[i] = 2.0 * ts_z_xz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyz[i] * gc_y[i] * tce_0;

        gs_y_z_xzz[i] = 2.0 * ts_z_xzz[i] * gc_y[i] * tce_0;

        gs_y_z_yyy[i] = 6.0 * ts_z_yy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyy[i] * gc_y[i] * tce_0;

        gs_y_z_yyz[i] = 4.0 * ts_z_yz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyz[i] * gc_y[i] * tce_0;

        gs_y_z_yzz[i] = 2.0 * ts_z_zz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yzz[i] * gc_y[i] * tce_0;

        gs_y_z_zzz[i] = 2.0 * ts_z_zzz[i] * gc_y[i] * tce_0;
    }

    // Set up 60-70 components of targeted buffer : PF

    auto gs_z_x_xxx = pbuffer.data(idx_g_pf + 60);

    auto gs_z_x_xxy = pbuffer.data(idx_g_pf + 61);

    auto gs_z_x_xxz = pbuffer.data(idx_g_pf + 62);

    auto gs_z_x_xyy = pbuffer.data(idx_g_pf + 63);

    auto gs_z_x_xyz = pbuffer.data(idx_g_pf + 64);

    auto gs_z_x_xzz = pbuffer.data(idx_g_pf + 65);

    auto gs_z_x_yyy = pbuffer.data(idx_g_pf + 66);

    auto gs_z_x_yyz = pbuffer.data(idx_g_pf + 67);

    auto gs_z_x_yzz = pbuffer.data(idx_g_pf + 68);

    auto gs_z_x_zzz = pbuffer.data(idx_g_pf + 69);

    #pragma omp simd aligned(gc_z, gs_z_x_xxx, gs_z_x_xxy, gs_z_x_xxz, gs_z_x_xyy, gs_z_x_xyz, gs_z_x_xzz, gs_z_x_yyy, gs_z_x_yyz, gs_z_x_yzz, gs_z_x_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_x_xxx[i] = 2.0 * ts_x_xxx[i] * gc_z[i] * tce_0;

        gs_z_x_xxy[i] = 2.0 * ts_x_xxy[i] * gc_z[i] * tce_0;

        gs_z_x_xxz[i] = 2.0 * ts_x_xx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxz[i] * gc_z[i] * tce_0;

        gs_z_x_xyy[i] = 2.0 * ts_x_xyy[i] * gc_z[i] * tce_0;

        gs_z_x_xyz[i] = 2.0 * ts_x_xy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyz[i] * gc_z[i] * tce_0;

        gs_z_x_xzz[i] = 4.0 * ts_x_xz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xzz[i] * gc_z[i] * tce_0;

        gs_z_x_yyy[i] = 2.0 * ts_x_yyy[i] * gc_z[i] * tce_0;

        gs_z_x_yyz[i] = 2.0 * ts_x_yy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyz[i] * gc_z[i] * tce_0;

        gs_z_x_yzz[i] = 4.0 * ts_x_yz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzz[i] * gc_z[i] * tce_0;

        gs_z_x_zzz[i] = 6.0 * ts_x_zz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 70-80 components of targeted buffer : PF

    auto gs_z_y_xxx = pbuffer.data(idx_g_pf + 70);

    auto gs_z_y_xxy = pbuffer.data(idx_g_pf + 71);

    auto gs_z_y_xxz = pbuffer.data(idx_g_pf + 72);

    auto gs_z_y_xyy = pbuffer.data(idx_g_pf + 73);

    auto gs_z_y_xyz = pbuffer.data(idx_g_pf + 74);

    auto gs_z_y_xzz = pbuffer.data(idx_g_pf + 75);

    auto gs_z_y_yyy = pbuffer.data(idx_g_pf + 76);

    auto gs_z_y_yyz = pbuffer.data(idx_g_pf + 77);

    auto gs_z_y_yzz = pbuffer.data(idx_g_pf + 78);

    auto gs_z_y_zzz = pbuffer.data(idx_g_pf + 79);

    #pragma omp simd aligned(gc_z, gs_z_y_xxx, gs_z_y_xxy, gs_z_y_xxz, gs_z_y_xyy, gs_z_y_xyz, gs_z_y_xzz, gs_z_y_yyy, gs_z_y_yyz, gs_z_y_yzz, gs_z_y_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_y_xxx[i] = 2.0 * ts_y_xxx[i] * gc_z[i] * tce_0;

        gs_z_y_xxy[i] = 2.0 * ts_y_xxy[i] * gc_z[i] * tce_0;

        gs_z_y_xxz[i] = 2.0 * ts_y_xx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxz[i] * gc_z[i] * tce_0;

        gs_z_y_xyy[i] = 2.0 * ts_y_xyy[i] * gc_z[i] * tce_0;

        gs_z_y_xyz[i] = 2.0 * ts_y_xy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyz[i] * gc_z[i] * tce_0;

        gs_z_y_xzz[i] = 4.0 * ts_y_xz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzz[i] * gc_z[i] * tce_0;

        gs_z_y_yyy[i] = 2.0 * ts_y_yyy[i] * gc_z[i] * tce_0;

        gs_z_y_yyz[i] = 2.0 * ts_y_yy[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyz[i] * gc_z[i] * tce_0;

        gs_z_y_yzz[i] = 4.0 * ts_y_yz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yzz[i] * gc_z[i] * tce_0;

        gs_z_y_zzz[i] = 6.0 * ts_y_zz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzz[i] * gc_z[i] * tce_0;
    }

    // Set up 80-90 components of targeted buffer : PF

    auto gs_z_z_xxx = pbuffer.data(idx_g_pf + 80);

    auto gs_z_z_xxy = pbuffer.data(idx_g_pf + 81);

    auto gs_z_z_xxz = pbuffer.data(idx_g_pf + 82);

    auto gs_z_z_xyy = pbuffer.data(idx_g_pf + 83);

    auto gs_z_z_xyz = pbuffer.data(idx_g_pf + 84);

    auto gs_z_z_xzz = pbuffer.data(idx_g_pf + 85);

    auto gs_z_z_yyy = pbuffer.data(idx_g_pf + 86);

    auto gs_z_z_yyz = pbuffer.data(idx_g_pf + 87);

    auto gs_z_z_yzz = pbuffer.data(idx_g_pf + 88);

    auto gs_z_z_zzz = pbuffer.data(idx_g_pf + 89);

    #pragma omp simd aligned(gc_z, gs_z_z_xxx, gs_z_z_xxy, gs_z_z_xxz, gs_z_z_xyy, gs_z_z_xyz, gs_z_z_xzz, gs_z_z_yyy, gs_z_z_yyz, gs_z_z_yzz, gs_z_z_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_z_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxx[i] * gc_z[i] * tce_0;

        gs_z_z_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxy[i] * gc_z[i] * tce_0;

        gs_z_z_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxz[i] * gc_z[i] * tce_0;

        gs_z_z_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyy[i] * gc_z[i] * tce_0;

        gs_z_z_xyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyz[i] * gc_z[i] * tce_0;

        gs_z_z_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_xz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xzz[i] * gc_z[i] * tce_0;

        gs_z_z_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyy[i] * gc_z[i] * tce_0;

        gs_z_z_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyz[i] * gc_z[i] * tce_0;

        gs_z_z_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_yz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yzz[i] * gc_z[i] * tce_0;

        gs_z_z_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * tce_0 + 6.0 * ts_z_zz[i] * gfe_0 * tce_0 + 2.0 * ts_z_zzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

