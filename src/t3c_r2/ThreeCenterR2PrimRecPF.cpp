#include "ThreeCenterR2PrimRecPF.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_pf(CSimdArray<double>& pbuffer, 
                const size_t idx_g_pf,
                const size_t idx_sd,
                const size_t idx_sf,
                const size_t idx_pp,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_xy = pbuffer.data(idx_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

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

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_pp);

    auto ts_x_y = pbuffer.data(idx_pp + 1);

    auto ts_x_z = pbuffer.data(idx_pp + 2);

    auto ts_y_x = pbuffer.data(idx_pp + 3);

    auto ts_y_y = pbuffer.data(idx_pp + 4);

    auto ts_y_z = pbuffer.data(idx_pp + 5);

    auto ts_z_x = pbuffer.data(idx_pp + 6);

    auto ts_z_y = pbuffer.data(idx_pp + 7);

    auto ts_z_z = pbuffer.data(idx_pp + 8);

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

    auto gr_x_xxx = pbuffer.data(idx_g_pf);

    auto gr_x_xxy = pbuffer.data(idx_g_pf + 1);

    auto gr_x_xxz = pbuffer.data(idx_g_pf + 2);

    auto gr_x_xyy = pbuffer.data(idx_g_pf + 3);

    auto gr_x_xyz = pbuffer.data(idx_g_pf + 4);

    auto gr_x_xzz = pbuffer.data(idx_g_pf + 5);

    auto gr_x_yyy = pbuffer.data(idx_g_pf + 6);

    auto gr_x_yyz = pbuffer.data(idx_g_pf + 7);

    auto gr_x_yzz = pbuffer.data(idx_g_pf + 8);

    auto gr_x_zzz = pbuffer.data(idx_g_pf + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xyy, gr_x_xyz, gr_x_xzz, gr_x_yyy, gr_x_yyz, gr_x_yzz, gr_x_zzz, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, ts_x_x, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_y, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_z, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_x_xxx[i] = 6.0 * ts_0_xx[i] * gfe_0 + 2.0 * ts_0_xxx[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_x[i] * gfe_0 + 6.0 * ts_x_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_x_xxx[i] * gfe_0 + ts_x_xxx[i] * rgc2_0;

        gr_x_xxy[i] = 4.0 * ts_0_xy[i] * gfe_0 + 2.0 * ts_0_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_y[i] * gfe_0 + 4.0 * ts_x_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_xxy[i] * gfe_0 + ts_x_xxy[i] * rgc2_0;

        gr_x_xxz[i] = 4.0 * ts_0_xz[i] * gfe_0 + 2.0 * ts_0_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 + 4.0 * ts_x_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xxz[i] * gfe_0 + ts_x_xxz[i] * rgc2_0;

        gr_x_xyy[i] = 2.0 * ts_0_yy[i] * gfe_0 + 2.0 * ts_0_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 + 4.0 * ts_x_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_xyy[i] * gfe_0 + ts_x_xyy[i] * rgc2_0;

        gr_x_xyz[i] = 2.0 * ts_0_yz[i] * gfe_0 + 2.0 * ts_0_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xyz[i] * gfe_0 + ts_x_xyz[i] * rgc2_0;

        gr_x_xzz[i] = 2.0 * ts_0_zz[i] * gfe_0 + 2.0 * ts_0_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 + 4.0 * ts_x_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xzz[i] * gfe_0 + ts_x_xzz[i] * rgc2_0;

        gr_x_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_y[i] * gfe_0 + 6.0 * ts_x_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_yyy[i] * gfe_0 + ts_x_yyy[i] * rgc2_0;

        gr_x_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 + 4.0 * ts_x_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_yyz[i] * gfe_0 + ts_x_yyz[i] * rgc2_0;

        gr_x_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_y[i] * gfe_0 + 4.0 * ts_x_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_yzz[i] * gfe_0 + ts_x_yzz[i] * rgc2_0;

        gr_x_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_x_z[i] * gfe_0 + 6.0 * ts_x_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_zzz[i] * gfe_0 + ts_x_zzz[i] * rgc2_0;
    }

    // Set up 10-20 components of targeted buffer : PF

    auto gr_y_xxx = pbuffer.data(idx_g_pf + 10);

    auto gr_y_xxy = pbuffer.data(idx_g_pf + 11);

    auto gr_y_xxz = pbuffer.data(idx_g_pf + 12);

    auto gr_y_xyy = pbuffer.data(idx_g_pf + 13);

    auto gr_y_xyz = pbuffer.data(idx_g_pf + 14);

    auto gr_y_xzz = pbuffer.data(idx_g_pf + 15);

    auto gr_y_yyy = pbuffer.data(idx_g_pf + 16);

    auto gr_y_yyz = pbuffer.data(idx_g_pf + 17);

    auto gr_y_yzz = pbuffer.data(idx_g_pf + 18);

    auto gr_y_zzz = pbuffer.data(idx_g_pf + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xyy, gr_y_xyz, gr_y_xzz, gr_y_yyy, gr_y_yyz, gr_y_yzz, gr_y_zzz, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, ts_y_x, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_y, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_z, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_y_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_x[i] * gfe_0 + 6.0 * ts_y_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_y_xxx[i] * gfe_0 + ts_y_xxx[i] * rgc2_0;

        gr_y_xxy[i] = 2.0 * ts_0_xx[i] * gfe_0 + 2.0 * ts_0_xxy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 + 4.0 * ts_y_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_xxy[i] * gfe_0 + ts_y_xxy[i] * rgc2_0;

        gr_y_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_z[i] * gfe_0 + 4.0 * ts_y_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xxz[i] * gfe_0 + ts_y_xxz[i] * rgc2_0;

        gr_y_xyy[i] = 4.0 * ts_0_xy[i] * gfe_0 + 2.0 * ts_0_xyy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_x[i] * gfe_0 + 4.0 * ts_y_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_xyy[i] * gfe_0 + ts_y_xyy[i] * rgc2_0;

        gr_y_xyz[i] = 2.0 * ts_0_xz[i] * gfe_0 + 2.0 * ts_0_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xyz[i] * gfe_0 + ts_y_xyz[i] * rgc2_0;

        gr_y_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_x[i] * gfe_0 + 4.0 * ts_y_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xzz[i] * gfe_0 + ts_y_xzz[i] * rgc2_0;

        gr_y_yyy[i] = 6.0 * ts_0_yy[i] * gfe_0 + 2.0 * ts_0_yyy[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_y[i] * gfe_0 + 6.0 * ts_y_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_yyy[i] * gfe_0 + ts_y_yyy[i] * rgc2_0;

        gr_y_yyz[i] = 4.0 * ts_0_yz[i] * gfe_0 + 2.0 * ts_0_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_z[i] * gfe_0 + 4.0 * ts_y_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_yyz[i] * gfe_0 + ts_y_yyz[i] * rgc2_0;

        gr_y_yzz[i] = 2.0 * ts_0_zz[i] * gfe_0 + 2.0 * ts_0_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 + 4.0 * ts_y_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_yzz[i] * gfe_0 + ts_y_yzz[i] * rgc2_0;

        gr_y_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_y_z[i] * gfe_0 + 6.0 * ts_y_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_zzz[i] * gfe_0 + ts_y_zzz[i] * rgc2_0;
    }

    // Set up 20-30 components of targeted buffer : PF

    auto gr_z_xxx = pbuffer.data(idx_g_pf + 20);

    auto gr_z_xxy = pbuffer.data(idx_g_pf + 21);

    auto gr_z_xxz = pbuffer.data(idx_g_pf + 22);

    auto gr_z_xyy = pbuffer.data(idx_g_pf + 23);

    auto gr_z_xyz = pbuffer.data(idx_g_pf + 24);

    auto gr_z_xzz = pbuffer.data(idx_g_pf + 25);

    auto gr_z_yyy = pbuffer.data(idx_g_pf + 26);

    auto gr_z_yyz = pbuffer.data(idx_g_pf + 27);

    auto gr_z_yzz = pbuffer.data(idx_g_pf + 28);

    auto gr_z_zzz = pbuffer.data(idx_g_pf + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xyy, gr_z_xyz, gr_z_xzz, gr_z_yyy, gr_z_yyz, gr_z_yzz, gr_z_zzz, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, ts_z_x, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_y, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_z, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_z_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * gc_z[i] + 6.0 * ts_z_x[i] * gfe_0 + 6.0 * ts_z_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_z_xxx[i] * gfe_0 + ts_z_xxx[i] * rgc2_0;

        gr_z_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_y[i] * gfe_0 + 4.0 * ts_z_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_xxy[i] * gfe_0 + ts_z_xxy[i] * rgc2_0;

        gr_z_xxz[i] = 2.0 * ts_0_xx[i] * gfe_0 + 2.0 * ts_0_xxz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_z[i] * gfe_0 + 4.0 * ts_z_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xxz[i] * gfe_0 + ts_z_xxz[i] * rgc2_0;

        gr_z_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_x[i] * gfe_0 + 4.0 * ts_z_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_xyy[i] * gfe_0 + ts_z_xyy[i] * rgc2_0;

        gr_z_xyz[i] = 2.0 * ts_0_xy[i] * gfe_0 + 2.0 * ts_0_xyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xyz[i] * gfe_0 + ts_z_xyz[i] * rgc2_0;

        gr_z_xzz[i] = 4.0 * ts_0_xz[i] * gfe_0 + 2.0 * ts_0_xzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_x[i] * gfe_0 + 4.0 * ts_z_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xzz[i] * gfe_0 + ts_z_xzz[i] * rgc2_0;

        gr_z_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * gc_z[i] + 6.0 * ts_z_y[i] * gfe_0 + 6.0 * ts_z_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_yyy[i] * gfe_0 + ts_z_yyy[i] * rgc2_0;

        gr_z_yyz[i] = 2.0 * ts_0_yy[i] * gfe_0 + 2.0 * ts_0_yyz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_z[i] * gfe_0 + 4.0 * ts_z_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_yyz[i] * gfe_0 + ts_z_yyz[i] * rgc2_0;

        gr_z_yzz[i] = 4.0 * ts_0_yz[i] * gfe_0 + 2.0 * ts_0_yzz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_y[i] * gfe_0 + 4.0 * ts_z_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_yzz[i] * gfe_0 + ts_z_yzz[i] * rgc2_0;

        gr_z_zzz[i] = 6.0 * ts_0_zz[i] * gfe_0 + 2.0 * ts_0_zzz[i] * gfe_0 * gc_z[i] + 6.0 * ts_z_z[i] * gfe_0 + 6.0 * ts_z_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_zzz[i] * gfe_0 + ts_z_zzz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

