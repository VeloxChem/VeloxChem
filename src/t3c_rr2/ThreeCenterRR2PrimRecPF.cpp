#include "ThreeCenterRR2PrimRecPF.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_pf(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_pf,
                  const size_t idx_sf,
                  const size_t idx_g_sf,
                  const size_t idx_pd,
                  const size_t idx_g_pd,
                  const size_t idx_pf,
                  const size_t idx_g_pf,
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

    // Set up components of auxiliary buffer : SF

    auto gr_0_xxx = pbuffer.data(idx_g_sf);

    auto gr_0_xxy = pbuffer.data(idx_g_sf + 1);

    auto gr_0_xxz = pbuffer.data(idx_g_sf + 2);

    auto gr_0_xyy = pbuffer.data(idx_g_sf + 3);

    auto gr_0_xyz = pbuffer.data(idx_g_sf + 4);

    auto gr_0_xzz = pbuffer.data(idx_g_sf + 5);

    auto gr_0_yyy = pbuffer.data(idx_g_sf + 6);

    auto gr_0_yyz = pbuffer.data(idx_g_sf + 7);

    auto gr_0_yzz = pbuffer.data(idx_g_sf + 8);

    auto gr_0_zzz = pbuffer.data(idx_g_sf + 9);

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

    // Set up components of auxiliary buffer : PD

    auto gr_x_xx = pbuffer.data(idx_g_pd);

    auto gr_x_xy = pbuffer.data(idx_g_pd + 1);

    auto gr_x_xz = pbuffer.data(idx_g_pd + 2);

    auto gr_x_yy = pbuffer.data(idx_g_pd + 3);

    auto gr_x_yz = pbuffer.data(idx_g_pd + 4);

    auto gr_x_zz = pbuffer.data(idx_g_pd + 5);

    auto gr_y_xx = pbuffer.data(idx_g_pd + 6);

    auto gr_y_xy = pbuffer.data(idx_g_pd + 7);

    auto gr_y_xz = pbuffer.data(idx_g_pd + 8);

    auto gr_y_yy = pbuffer.data(idx_g_pd + 9);

    auto gr_y_yz = pbuffer.data(idx_g_pd + 10);

    auto gr_y_zz = pbuffer.data(idx_g_pd + 11);

    auto gr_z_xx = pbuffer.data(idx_g_pd + 12);

    auto gr_z_xy = pbuffer.data(idx_g_pd + 13);

    auto gr_z_xz = pbuffer.data(idx_g_pd + 14);

    auto gr_z_yy = pbuffer.data(idx_g_pd + 15);

    auto gr_z_yz = pbuffer.data(idx_g_pd + 16);

    auto gr_z_zz = pbuffer.data(idx_g_pd + 17);

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

    // Set up components of auxiliary buffer : PF

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

    // Set up 0-10 components of targeted buffer : PF

    auto grr_x_x_xxx = pbuffer.data(idx_gr_pf);

    auto grr_x_x_xxy = pbuffer.data(idx_gr_pf + 1);

    auto grr_x_x_xxz = pbuffer.data(idx_gr_pf + 2);

    auto grr_x_x_xyy = pbuffer.data(idx_gr_pf + 3);

    auto grr_x_x_xyz = pbuffer.data(idx_gr_pf + 4);

    auto grr_x_x_xzz = pbuffer.data(idx_gr_pf + 5);

    auto grr_x_x_yyy = pbuffer.data(idx_gr_pf + 6);

    auto grr_x_x_yyz = pbuffer.data(idx_gr_pf + 7);

    auto grr_x_x_yzz = pbuffer.data(idx_gr_pf + 8);

    auto grr_x_x_zzz = pbuffer.data(idx_gr_pf + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxx, gr_0_xxy, gr_0_xxz, gr_0_xyy, gr_0_xyz, gr_0_xzz, gr_0_yyy, gr_0_yyz, gr_0_yzz, gr_0_zzz, gr_x_xx, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xy, gr_x_xyy, gr_x_xyz, gr_x_xz, gr_x_xzz, gr_x_yy, gr_x_yyy, gr_x_yyz, gr_x_yz, gr_x_yzz, gr_x_zz, gr_x_zzz, grr_x_x_xxx, grr_x_x_xxy, grr_x_x_xxz, grr_x_x_xyy, grr_x_x_xyz, grr_x_x_xzz, grr_x_x_yyy, grr_x_x_yyz, grr_x_x_yzz, grr_x_x_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_x_xxx[i] = 2.0 * ts_0_xxx[i] * gfe2_0 + gr_0_xxx[i] * gfe_0 + 6.0 * ts_x_xx[i] * gfe2_0 + 3.0 * gr_x_xx[i] * gfe_0 + 2.0 * ts_x_xxx[i] * gfe_0 * gc_x[i] + gr_x_xxx[i] * gc_x[i];

        grr_x_x_xxy[i] = 2.0 * ts_0_xxy[i] * gfe2_0 + gr_0_xxy[i] * gfe_0 + 4.0 * ts_x_xy[i] * gfe2_0 + 2.0 * gr_x_xy[i] * gfe_0 + 2.0 * ts_x_xxy[i] * gfe_0 * gc_x[i] + gr_x_xxy[i] * gc_x[i];

        grr_x_x_xxz[i] = 2.0 * ts_0_xxz[i] * gfe2_0 + gr_0_xxz[i] * gfe_0 + 4.0 * ts_x_xz[i] * gfe2_0 + 2.0 * gr_x_xz[i] * gfe_0 + 2.0 * ts_x_xxz[i] * gfe_0 * gc_x[i] + gr_x_xxz[i] * gc_x[i];

        grr_x_x_xyy[i] = 2.0 * ts_0_xyy[i] * gfe2_0 + gr_0_xyy[i] * gfe_0 + 2.0 * ts_x_yy[i] * gfe2_0 + gr_x_yy[i] * gfe_0 + 2.0 * ts_x_xyy[i] * gfe_0 * gc_x[i] + gr_x_xyy[i] * gc_x[i];

        grr_x_x_xyz[i] = 2.0 * ts_0_xyz[i] * gfe2_0 + gr_0_xyz[i] * gfe_0 + 2.0 * ts_x_yz[i] * gfe2_0 + gr_x_yz[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe_0 * gc_x[i] + gr_x_xyz[i] * gc_x[i];

        grr_x_x_xzz[i] = 2.0 * ts_0_xzz[i] * gfe2_0 + gr_0_xzz[i] * gfe_0 + 2.0 * ts_x_zz[i] * gfe2_0 + gr_x_zz[i] * gfe_0 + 2.0 * ts_x_xzz[i] * gfe_0 * gc_x[i] + gr_x_xzz[i] * gc_x[i];

        grr_x_x_yyy[i] = 2.0 * ts_0_yyy[i] * gfe2_0 + gr_0_yyy[i] * gfe_0 + 2.0 * ts_x_yyy[i] * gfe_0 * gc_x[i] + gr_x_yyy[i] * gc_x[i];

        grr_x_x_yyz[i] = 2.0 * ts_0_yyz[i] * gfe2_0 + gr_0_yyz[i] * gfe_0 + 2.0 * ts_x_yyz[i] * gfe_0 * gc_x[i] + gr_x_yyz[i] * gc_x[i];

        grr_x_x_yzz[i] = 2.0 * ts_0_yzz[i] * gfe2_0 + gr_0_yzz[i] * gfe_0 + 2.0 * ts_x_yzz[i] * gfe_0 * gc_x[i] + gr_x_yzz[i] * gc_x[i];

        grr_x_x_zzz[i] = 2.0 * ts_0_zzz[i] * gfe2_0 + gr_0_zzz[i] * gfe_0 + 2.0 * ts_x_zzz[i] * gfe_0 * gc_x[i] + gr_x_zzz[i] * gc_x[i];
    }

    // Set up 10-20 components of targeted buffer : PF

    auto grr_x_y_xxx = pbuffer.data(idx_gr_pf + 10);

    auto grr_x_y_xxy = pbuffer.data(idx_gr_pf + 11);

    auto grr_x_y_xxz = pbuffer.data(idx_gr_pf + 12);

    auto grr_x_y_xyy = pbuffer.data(idx_gr_pf + 13);

    auto grr_x_y_xyz = pbuffer.data(idx_gr_pf + 14);

    auto grr_x_y_xzz = pbuffer.data(idx_gr_pf + 15);

    auto grr_x_y_yyy = pbuffer.data(idx_gr_pf + 16);

    auto grr_x_y_yyz = pbuffer.data(idx_gr_pf + 17);

    auto grr_x_y_yzz = pbuffer.data(idx_gr_pf + 18);

    auto grr_x_y_zzz = pbuffer.data(idx_gr_pf + 19);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xx, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xy, gr_y_xyy, gr_y_xyz, gr_y_xz, gr_y_xzz, gr_y_yy, gr_y_yyy, gr_y_yyz, gr_y_yz, gr_y_yzz, gr_y_zz, gr_y_zzz, grr_x_y_xxx, grr_x_y_xxy, grr_x_y_xxz, grr_x_y_xyy, grr_x_y_xyz, grr_x_y_xzz, grr_x_y_yyy, grr_x_y_yyz, grr_x_y_yzz, grr_x_y_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_y_xxx[i] = 6.0 * ts_y_xx[i] * gfe2_0 + 3.0 * gr_y_xx[i] * gfe_0 + 2.0 * ts_y_xxx[i] * gfe_0 * gc_x[i] + gr_y_xxx[i] * gc_x[i];

        grr_x_y_xxy[i] = 4.0 * ts_y_xy[i] * gfe2_0 + 2.0 * gr_y_xy[i] * gfe_0 + 2.0 * ts_y_xxy[i] * gfe_0 * gc_x[i] + gr_y_xxy[i] * gc_x[i];

        grr_x_y_xxz[i] = 4.0 * ts_y_xz[i] * gfe2_0 + 2.0 * gr_y_xz[i] * gfe_0 + 2.0 * ts_y_xxz[i] * gfe_0 * gc_x[i] + gr_y_xxz[i] * gc_x[i];

        grr_x_y_xyy[i] = 2.0 * ts_y_yy[i] * gfe2_0 + gr_y_yy[i] * gfe_0 + 2.0 * ts_y_xyy[i] * gfe_0 * gc_x[i] + gr_y_xyy[i] * gc_x[i];

        grr_x_y_xyz[i] = 2.0 * ts_y_yz[i] * gfe2_0 + gr_y_yz[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe_0 * gc_x[i] + gr_y_xyz[i] * gc_x[i];

        grr_x_y_xzz[i] = 2.0 * ts_y_zz[i] * gfe2_0 + gr_y_zz[i] * gfe_0 + 2.0 * ts_y_xzz[i] * gfe_0 * gc_x[i] + gr_y_xzz[i] * gc_x[i];

        grr_x_y_yyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 * gc_x[i] + gr_y_yyy[i] * gc_x[i];

        grr_x_y_yyz[i] = 2.0 * ts_y_yyz[i] * gfe_0 * gc_x[i] + gr_y_yyz[i] * gc_x[i];

        grr_x_y_yzz[i] = 2.0 * ts_y_yzz[i] * gfe_0 * gc_x[i] + gr_y_yzz[i] * gc_x[i];

        grr_x_y_zzz[i] = 2.0 * ts_y_zzz[i] * gfe_0 * gc_x[i] + gr_y_zzz[i] * gc_x[i];
    }

    // Set up 20-30 components of targeted buffer : PF

    auto grr_x_z_xxx = pbuffer.data(idx_gr_pf + 20);

    auto grr_x_z_xxy = pbuffer.data(idx_gr_pf + 21);

    auto grr_x_z_xxz = pbuffer.data(idx_gr_pf + 22);

    auto grr_x_z_xyy = pbuffer.data(idx_gr_pf + 23);

    auto grr_x_z_xyz = pbuffer.data(idx_gr_pf + 24);

    auto grr_x_z_xzz = pbuffer.data(idx_gr_pf + 25);

    auto grr_x_z_yyy = pbuffer.data(idx_gr_pf + 26);

    auto grr_x_z_yyz = pbuffer.data(idx_gr_pf + 27);

    auto grr_x_z_yzz = pbuffer.data(idx_gr_pf + 28);

    auto grr_x_z_zzz = pbuffer.data(idx_gr_pf + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xx, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xy, gr_z_xyy, gr_z_xyz, gr_z_xz, gr_z_xzz, gr_z_yy, gr_z_yyy, gr_z_yyz, gr_z_yz, gr_z_yzz, gr_z_zz, gr_z_zzz, grr_x_z_xxx, grr_x_z_xxy, grr_x_z_xxz, grr_x_z_xyy, grr_x_z_xyz, grr_x_z_xzz, grr_x_z_yyy, grr_x_z_yyz, grr_x_z_yzz, grr_x_z_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_z_xxx[i] = 6.0 * ts_z_xx[i] * gfe2_0 + 3.0 * gr_z_xx[i] * gfe_0 + 2.0 * ts_z_xxx[i] * gfe_0 * gc_x[i] + gr_z_xxx[i] * gc_x[i];

        grr_x_z_xxy[i] = 4.0 * ts_z_xy[i] * gfe2_0 + 2.0 * gr_z_xy[i] * gfe_0 + 2.0 * ts_z_xxy[i] * gfe_0 * gc_x[i] + gr_z_xxy[i] * gc_x[i];

        grr_x_z_xxz[i] = 4.0 * ts_z_xz[i] * gfe2_0 + 2.0 * gr_z_xz[i] * gfe_0 + 2.0 * ts_z_xxz[i] * gfe_0 * gc_x[i] + gr_z_xxz[i] * gc_x[i];

        grr_x_z_xyy[i] = 2.0 * ts_z_yy[i] * gfe2_0 + gr_z_yy[i] * gfe_0 + 2.0 * ts_z_xyy[i] * gfe_0 * gc_x[i] + gr_z_xyy[i] * gc_x[i];

        grr_x_z_xyz[i] = 2.0 * ts_z_yz[i] * gfe2_0 + gr_z_yz[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe_0 * gc_x[i] + gr_z_xyz[i] * gc_x[i];

        grr_x_z_xzz[i] = 2.0 * ts_z_zz[i] * gfe2_0 + gr_z_zz[i] * gfe_0 + 2.0 * ts_z_xzz[i] * gfe_0 * gc_x[i] + gr_z_xzz[i] * gc_x[i];

        grr_x_z_yyy[i] = 2.0 * ts_z_yyy[i] * gfe_0 * gc_x[i] + gr_z_yyy[i] * gc_x[i];

        grr_x_z_yyz[i] = 2.0 * ts_z_yyz[i] * gfe_0 * gc_x[i] + gr_z_yyz[i] * gc_x[i];

        grr_x_z_yzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 * gc_x[i] + gr_z_yzz[i] * gc_x[i];

        grr_x_z_zzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * gc_x[i] + gr_z_zzz[i] * gc_x[i];
    }

    // Set up 30-40 components of targeted buffer : PF

    auto grr_y_x_xxx = pbuffer.data(idx_gr_pf + 30);

    auto grr_y_x_xxy = pbuffer.data(idx_gr_pf + 31);

    auto grr_y_x_xxz = pbuffer.data(idx_gr_pf + 32);

    auto grr_y_x_xyy = pbuffer.data(idx_gr_pf + 33);

    auto grr_y_x_xyz = pbuffer.data(idx_gr_pf + 34);

    auto grr_y_x_xzz = pbuffer.data(idx_gr_pf + 35);

    auto grr_y_x_yyy = pbuffer.data(idx_gr_pf + 36);

    auto grr_y_x_yyz = pbuffer.data(idx_gr_pf + 37);

    auto grr_y_x_yzz = pbuffer.data(idx_gr_pf + 38);

    auto grr_y_x_zzz = pbuffer.data(idx_gr_pf + 39);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xx, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xy, gr_x_xyy, gr_x_xyz, gr_x_xz, gr_x_xzz, gr_x_yy, gr_x_yyy, gr_x_yyz, gr_x_yz, gr_x_yzz, gr_x_zz, gr_x_zzz, grr_y_x_xxx, grr_y_x_xxy, grr_y_x_xxz, grr_y_x_xyy, grr_y_x_xyz, grr_y_x_xzz, grr_y_x_yyy, grr_y_x_yyz, grr_y_x_yzz, grr_y_x_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_x_xxx[i] = 2.0 * ts_x_xxx[i] * gfe_0 * gc_y[i] + gr_x_xxx[i] * gc_y[i];

        grr_y_x_xxy[i] = 2.0 * ts_x_xx[i] * gfe2_0 + gr_x_xx[i] * gfe_0 + 2.0 * ts_x_xxy[i] * gfe_0 * gc_y[i] + gr_x_xxy[i] * gc_y[i];

        grr_y_x_xxz[i] = 2.0 * ts_x_xxz[i] * gfe_0 * gc_y[i] + gr_x_xxz[i] * gc_y[i];

        grr_y_x_xyy[i] = 4.0 * ts_x_xy[i] * gfe2_0 + 2.0 * gr_x_xy[i] * gfe_0 + 2.0 * ts_x_xyy[i] * gfe_0 * gc_y[i] + gr_x_xyy[i] * gc_y[i];

        grr_y_x_xyz[i] = 2.0 * ts_x_xz[i] * gfe2_0 + gr_x_xz[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe_0 * gc_y[i] + gr_x_xyz[i] * gc_y[i];

        grr_y_x_xzz[i] = 2.0 * ts_x_xzz[i] * gfe_0 * gc_y[i] + gr_x_xzz[i] * gc_y[i];

        grr_y_x_yyy[i] = 6.0 * ts_x_yy[i] * gfe2_0 + 3.0 * gr_x_yy[i] * gfe_0 + 2.0 * ts_x_yyy[i] * gfe_0 * gc_y[i] + gr_x_yyy[i] * gc_y[i];

        grr_y_x_yyz[i] = 4.0 * ts_x_yz[i] * gfe2_0 + 2.0 * gr_x_yz[i] * gfe_0 + 2.0 * ts_x_yyz[i] * gfe_0 * gc_y[i] + gr_x_yyz[i] * gc_y[i];

        grr_y_x_yzz[i] = 2.0 * ts_x_zz[i] * gfe2_0 + gr_x_zz[i] * gfe_0 + 2.0 * ts_x_yzz[i] * gfe_0 * gc_y[i] + gr_x_yzz[i] * gc_y[i];

        grr_y_x_zzz[i] = 2.0 * ts_x_zzz[i] * gfe_0 * gc_y[i] + gr_x_zzz[i] * gc_y[i];
    }

    // Set up 40-50 components of targeted buffer : PF

    auto grr_y_y_xxx = pbuffer.data(idx_gr_pf + 40);

    auto grr_y_y_xxy = pbuffer.data(idx_gr_pf + 41);

    auto grr_y_y_xxz = pbuffer.data(idx_gr_pf + 42);

    auto grr_y_y_xyy = pbuffer.data(idx_gr_pf + 43);

    auto grr_y_y_xyz = pbuffer.data(idx_gr_pf + 44);

    auto grr_y_y_xzz = pbuffer.data(idx_gr_pf + 45);

    auto grr_y_y_yyy = pbuffer.data(idx_gr_pf + 46);

    auto grr_y_y_yyz = pbuffer.data(idx_gr_pf + 47);

    auto grr_y_y_yzz = pbuffer.data(idx_gr_pf + 48);

    auto grr_y_y_zzz = pbuffer.data(idx_gr_pf + 49);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxx, gr_0_xxy, gr_0_xxz, gr_0_xyy, gr_0_xyz, gr_0_xzz, gr_0_yyy, gr_0_yyz, gr_0_yzz, gr_0_zzz, gr_y_xx, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xy, gr_y_xyy, gr_y_xyz, gr_y_xz, gr_y_xzz, gr_y_yy, gr_y_yyy, gr_y_yyz, gr_y_yz, gr_y_yzz, gr_y_zz, gr_y_zzz, grr_y_y_xxx, grr_y_y_xxy, grr_y_y_xxz, grr_y_y_xyy, grr_y_y_xyz, grr_y_y_xzz, grr_y_y_yyy, grr_y_y_yyz, grr_y_y_yzz, grr_y_y_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_y_xxx[i] = 2.0 * ts_0_xxx[i] * gfe2_0 + gr_0_xxx[i] * gfe_0 + 2.0 * ts_y_xxx[i] * gfe_0 * gc_y[i] + gr_y_xxx[i] * gc_y[i];

        grr_y_y_xxy[i] = 2.0 * ts_0_xxy[i] * gfe2_0 + gr_0_xxy[i] * gfe_0 + 2.0 * ts_y_xx[i] * gfe2_0 + gr_y_xx[i] * gfe_0 + 2.0 * ts_y_xxy[i] * gfe_0 * gc_y[i] + gr_y_xxy[i] * gc_y[i];

        grr_y_y_xxz[i] = 2.0 * ts_0_xxz[i] * gfe2_0 + gr_0_xxz[i] * gfe_0 + 2.0 * ts_y_xxz[i] * gfe_0 * gc_y[i] + gr_y_xxz[i] * gc_y[i];

        grr_y_y_xyy[i] = 2.0 * ts_0_xyy[i] * gfe2_0 + gr_0_xyy[i] * gfe_0 + 4.0 * ts_y_xy[i] * gfe2_0 + 2.0 * gr_y_xy[i] * gfe_0 + 2.0 * ts_y_xyy[i] * gfe_0 * gc_y[i] + gr_y_xyy[i] * gc_y[i];

        grr_y_y_xyz[i] = 2.0 * ts_0_xyz[i] * gfe2_0 + gr_0_xyz[i] * gfe_0 + 2.0 * ts_y_xz[i] * gfe2_0 + gr_y_xz[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe_0 * gc_y[i] + gr_y_xyz[i] * gc_y[i];

        grr_y_y_xzz[i] = 2.0 * ts_0_xzz[i] * gfe2_0 + gr_0_xzz[i] * gfe_0 + 2.0 * ts_y_xzz[i] * gfe_0 * gc_y[i] + gr_y_xzz[i] * gc_y[i];

        grr_y_y_yyy[i] = 2.0 * ts_0_yyy[i] * gfe2_0 + gr_0_yyy[i] * gfe_0 + 6.0 * ts_y_yy[i] * gfe2_0 + 3.0 * gr_y_yy[i] * gfe_0 + 2.0 * ts_y_yyy[i] * gfe_0 * gc_y[i] + gr_y_yyy[i] * gc_y[i];

        grr_y_y_yyz[i] = 2.0 * ts_0_yyz[i] * gfe2_0 + gr_0_yyz[i] * gfe_0 + 4.0 * ts_y_yz[i] * gfe2_0 + 2.0 * gr_y_yz[i] * gfe_0 + 2.0 * ts_y_yyz[i] * gfe_0 * gc_y[i] + gr_y_yyz[i] * gc_y[i];

        grr_y_y_yzz[i] = 2.0 * ts_0_yzz[i] * gfe2_0 + gr_0_yzz[i] * gfe_0 + 2.0 * ts_y_zz[i] * gfe2_0 + gr_y_zz[i] * gfe_0 + 2.0 * ts_y_yzz[i] * gfe_0 * gc_y[i] + gr_y_yzz[i] * gc_y[i];

        grr_y_y_zzz[i] = 2.0 * ts_0_zzz[i] * gfe2_0 + gr_0_zzz[i] * gfe_0 + 2.0 * ts_y_zzz[i] * gfe_0 * gc_y[i] + gr_y_zzz[i] * gc_y[i];
    }

    // Set up 50-60 components of targeted buffer : PF

    auto grr_y_z_xxx = pbuffer.data(idx_gr_pf + 50);

    auto grr_y_z_xxy = pbuffer.data(idx_gr_pf + 51);

    auto grr_y_z_xxz = pbuffer.data(idx_gr_pf + 52);

    auto grr_y_z_xyy = pbuffer.data(idx_gr_pf + 53);

    auto grr_y_z_xyz = pbuffer.data(idx_gr_pf + 54);

    auto grr_y_z_xzz = pbuffer.data(idx_gr_pf + 55);

    auto grr_y_z_yyy = pbuffer.data(idx_gr_pf + 56);

    auto grr_y_z_yyz = pbuffer.data(idx_gr_pf + 57);

    auto grr_y_z_yzz = pbuffer.data(idx_gr_pf + 58);

    auto grr_y_z_zzz = pbuffer.data(idx_gr_pf + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xx, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xy, gr_z_xyy, gr_z_xyz, gr_z_xz, gr_z_xzz, gr_z_yy, gr_z_yyy, gr_z_yyz, gr_z_yz, gr_z_yzz, gr_z_zz, gr_z_zzz, grr_y_z_xxx, grr_y_z_xxy, grr_y_z_xxz, grr_y_z_xyy, grr_y_z_xyz, grr_y_z_xzz, grr_y_z_yyy, grr_y_z_yyz, grr_y_z_yzz, grr_y_z_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_z_xxx[i] = 2.0 * ts_z_xxx[i] * gfe_0 * gc_y[i] + gr_z_xxx[i] * gc_y[i];

        grr_y_z_xxy[i] = 2.0 * ts_z_xx[i] * gfe2_0 + gr_z_xx[i] * gfe_0 + 2.0 * ts_z_xxy[i] * gfe_0 * gc_y[i] + gr_z_xxy[i] * gc_y[i];

        grr_y_z_xxz[i] = 2.0 * ts_z_xxz[i] * gfe_0 * gc_y[i] + gr_z_xxz[i] * gc_y[i];

        grr_y_z_xyy[i] = 4.0 * ts_z_xy[i] * gfe2_0 + 2.0 * gr_z_xy[i] * gfe_0 + 2.0 * ts_z_xyy[i] * gfe_0 * gc_y[i] + gr_z_xyy[i] * gc_y[i];

        grr_y_z_xyz[i] = 2.0 * ts_z_xz[i] * gfe2_0 + gr_z_xz[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe_0 * gc_y[i] + gr_z_xyz[i] * gc_y[i];

        grr_y_z_xzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 * gc_y[i] + gr_z_xzz[i] * gc_y[i];

        grr_y_z_yyy[i] = 6.0 * ts_z_yy[i] * gfe2_0 + 3.0 * gr_z_yy[i] * gfe_0 + 2.0 * ts_z_yyy[i] * gfe_0 * gc_y[i] + gr_z_yyy[i] * gc_y[i];

        grr_y_z_yyz[i] = 4.0 * ts_z_yz[i] * gfe2_0 + 2.0 * gr_z_yz[i] * gfe_0 + 2.0 * ts_z_yyz[i] * gfe_0 * gc_y[i] + gr_z_yyz[i] * gc_y[i];

        grr_y_z_yzz[i] = 2.0 * ts_z_zz[i] * gfe2_0 + gr_z_zz[i] * gfe_0 + 2.0 * ts_z_yzz[i] * gfe_0 * gc_y[i] + gr_z_yzz[i] * gc_y[i];

        grr_y_z_zzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * gc_y[i] + gr_z_zzz[i] * gc_y[i];
    }

    // Set up 60-70 components of targeted buffer : PF

    auto grr_z_x_xxx = pbuffer.data(idx_gr_pf + 60);

    auto grr_z_x_xxy = pbuffer.data(idx_gr_pf + 61);

    auto grr_z_x_xxz = pbuffer.data(idx_gr_pf + 62);

    auto grr_z_x_xyy = pbuffer.data(idx_gr_pf + 63);

    auto grr_z_x_xyz = pbuffer.data(idx_gr_pf + 64);

    auto grr_z_x_xzz = pbuffer.data(idx_gr_pf + 65);

    auto grr_z_x_yyy = pbuffer.data(idx_gr_pf + 66);

    auto grr_z_x_yyz = pbuffer.data(idx_gr_pf + 67);

    auto grr_z_x_yzz = pbuffer.data(idx_gr_pf + 68);

    auto grr_z_x_zzz = pbuffer.data(idx_gr_pf + 69);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xx, gr_x_xxx, gr_x_xxy, gr_x_xxz, gr_x_xy, gr_x_xyy, gr_x_xyz, gr_x_xz, gr_x_xzz, gr_x_yy, gr_x_yyy, gr_x_yyz, gr_x_yz, gr_x_yzz, gr_x_zz, gr_x_zzz, grr_z_x_xxx, grr_z_x_xxy, grr_z_x_xxz, grr_z_x_xyy, grr_z_x_xyz, grr_z_x_xzz, grr_z_x_yyy, grr_z_x_yyz, grr_z_x_yzz, grr_z_x_zzz, ts_x_xx, ts_x_xxx, ts_x_xxy, ts_x_xxz, ts_x_xy, ts_x_xyy, ts_x_xyz, ts_x_xz, ts_x_xzz, ts_x_yy, ts_x_yyy, ts_x_yyz, ts_x_yz, ts_x_yzz, ts_x_zz, ts_x_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_x_xxx[i] = 2.0 * ts_x_xxx[i] * gfe_0 * gc_z[i] + gr_x_xxx[i] * gc_z[i];

        grr_z_x_xxy[i] = 2.0 * ts_x_xxy[i] * gfe_0 * gc_z[i] + gr_x_xxy[i] * gc_z[i];

        grr_z_x_xxz[i] = 2.0 * ts_x_xx[i] * gfe2_0 + gr_x_xx[i] * gfe_0 + 2.0 * ts_x_xxz[i] * gfe_0 * gc_z[i] + gr_x_xxz[i] * gc_z[i];

        grr_z_x_xyy[i] = 2.0 * ts_x_xyy[i] * gfe_0 * gc_z[i] + gr_x_xyy[i] * gc_z[i];

        grr_z_x_xyz[i] = 2.0 * ts_x_xy[i] * gfe2_0 + gr_x_xy[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe_0 * gc_z[i] + gr_x_xyz[i] * gc_z[i];

        grr_z_x_xzz[i] = 4.0 * ts_x_xz[i] * gfe2_0 + 2.0 * gr_x_xz[i] * gfe_0 + 2.0 * ts_x_xzz[i] * gfe_0 * gc_z[i] + gr_x_xzz[i] * gc_z[i];

        grr_z_x_yyy[i] = 2.0 * ts_x_yyy[i] * gfe_0 * gc_z[i] + gr_x_yyy[i] * gc_z[i];

        grr_z_x_yyz[i] = 2.0 * ts_x_yy[i] * gfe2_0 + gr_x_yy[i] * gfe_0 + 2.0 * ts_x_yyz[i] * gfe_0 * gc_z[i] + gr_x_yyz[i] * gc_z[i];

        grr_z_x_yzz[i] = 4.0 * ts_x_yz[i] * gfe2_0 + 2.0 * gr_x_yz[i] * gfe_0 + 2.0 * ts_x_yzz[i] * gfe_0 * gc_z[i] + gr_x_yzz[i] * gc_z[i];

        grr_z_x_zzz[i] = 6.0 * ts_x_zz[i] * gfe2_0 + 3.0 * gr_x_zz[i] * gfe_0 + 2.0 * ts_x_zzz[i] * gfe_0 * gc_z[i] + gr_x_zzz[i] * gc_z[i];
    }

    // Set up 70-80 components of targeted buffer : PF

    auto grr_z_y_xxx = pbuffer.data(idx_gr_pf + 70);

    auto grr_z_y_xxy = pbuffer.data(idx_gr_pf + 71);

    auto grr_z_y_xxz = pbuffer.data(idx_gr_pf + 72);

    auto grr_z_y_xyy = pbuffer.data(idx_gr_pf + 73);

    auto grr_z_y_xyz = pbuffer.data(idx_gr_pf + 74);

    auto grr_z_y_xzz = pbuffer.data(idx_gr_pf + 75);

    auto grr_z_y_yyy = pbuffer.data(idx_gr_pf + 76);

    auto grr_z_y_yyz = pbuffer.data(idx_gr_pf + 77);

    auto grr_z_y_yzz = pbuffer.data(idx_gr_pf + 78);

    auto grr_z_y_zzz = pbuffer.data(idx_gr_pf + 79);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xx, gr_y_xxx, gr_y_xxy, gr_y_xxz, gr_y_xy, gr_y_xyy, gr_y_xyz, gr_y_xz, gr_y_xzz, gr_y_yy, gr_y_yyy, gr_y_yyz, gr_y_yz, gr_y_yzz, gr_y_zz, gr_y_zzz, grr_z_y_xxx, grr_z_y_xxy, grr_z_y_xxz, grr_z_y_xyy, grr_z_y_xyz, grr_z_y_xzz, grr_z_y_yyy, grr_z_y_yyz, grr_z_y_yzz, grr_z_y_zzz, ts_y_xx, ts_y_xxx, ts_y_xxy, ts_y_xxz, ts_y_xy, ts_y_xyy, ts_y_xyz, ts_y_xz, ts_y_xzz, ts_y_yy, ts_y_yyy, ts_y_yyz, ts_y_yz, ts_y_yzz, ts_y_zz, ts_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_y_xxx[i] = 2.0 * ts_y_xxx[i] * gfe_0 * gc_z[i] + gr_y_xxx[i] * gc_z[i];

        grr_z_y_xxy[i] = 2.0 * ts_y_xxy[i] * gfe_0 * gc_z[i] + gr_y_xxy[i] * gc_z[i];

        grr_z_y_xxz[i] = 2.0 * ts_y_xx[i] * gfe2_0 + gr_y_xx[i] * gfe_0 + 2.0 * ts_y_xxz[i] * gfe_0 * gc_z[i] + gr_y_xxz[i] * gc_z[i];

        grr_z_y_xyy[i] = 2.0 * ts_y_xyy[i] * gfe_0 * gc_z[i] + gr_y_xyy[i] * gc_z[i];

        grr_z_y_xyz[i] = 2.0 * ts_y_xy[i] * gfe2_0 + gr_y_xy[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe_0 * gc_z[i] + gr_y_xyz[i] * gc_z[i];

        grr_z_y_xzz[i] = 4.0 * ts_y_xz[i] * gfe2_0 + 2.0 * gr_y_xz[i] * gfe_0 + 2.0 * ts_y_xzz[i] * gfe_0 * gc_z[i] + gr_y_xzz[i] * gc_z[i];

        grr_z_y_yyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 * gc_z[i] + gr_y_yyy[i] * gc_z[i];

        grr_z_y_yyz[i] = 2.0 * ts_y_yy[i] * gfe2_0 + gr_y_yy[i] * gfe_0 + 2.0 * ts_y_yyz[i] * gfe_0 * gc_z[i] + gr_y_yyz[i] * gc_z[i];

        grr_z_y_yzz[i] = 4.0 * ts_y_yz[i] * gfe2_0 + 2.0 * gr_y_yz[i] * gfe_0 + 2.0 * ts_y_yzz[i] * gfe_0 * gc_z[i] + gr_y_yzz[i] * gc_z[i];

        grr_z_y_zzz[i] = 6.0 * ts_y_zz[i] * gfe2_0 + 3.0 * gr_y_zz[i] * gfe_0 + 2.0 * ts_y_zzz[i] * gfe_0 * gc_z[i] + gr_y_zzz[i] * gc_z[i];
    }

    // Set up 80-90 components of targeted buffer : PF

    auto grr_z_z_xxx = pbuffer.data(idx_gr_pf + 80);

    auto grr_z_z_xxy = pbuffer.data(idx_gr_pf + 81);

    auto grr_z_z_xxz = pbuffer.data(idx_gr_pf + 82);

    auto grr_z_z_xyy = pbuffer.data(idx_gr_pf + 83);

    auto grr_z_z_xyz = pbuffer.data(idx_gr_pf + 84);

    auto grr_z_z_xzz = pbuffer.data(idx_gr_pf + 85);

    auto grr_z_z_yyy = pbuffer.data(idx_gr_pf + 86);

    auto grr_z_z_yyz = pbuffer.data(idx_gr_pf + 87);

    auto grr_z_z_yzz = pbuffer.data(idx_gr_pf + 88);

    auto grr_z_z_zzz = pbuffer.data(idx_gr_pf + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxx, gr_0_xxy, gr_0_xxz, gr_0_xyy, gr_0_xyz, gr_0_xzz, gr_0_yyy, gr_0_yyz, gr_0_yzz, gr_0_zzz, gr_z_xx, gr_z_xxx, gr_z_xxy, gr_z_xxz, gr_z_xy, gr_z_xyy, gr_z_xyz, gr_z_xz, gr_z_xzz, gr_z_yy, gr_z_yyy, gr_z_yyz, gr_z_yz, gr_z_yzz, gr_z_zz, gr_z_zzz, grr_z_z_xxx, grr_z_z_xxy, grr_z_z_xxz, grr_z_z_xyy, grr_z_z_xyz, grr_z_z_xzz, grr_z_z_yyy, grr_z_z_yyz, grr_z_z_yzz, grr_z_z_zzz, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xyy, ts_0_xyz, ts_0_xzz, ts_0_yyy, ts_0_yyz, ts_0_yzz, ts_0_zzz, ts_z_xx, ts_z_xxx, ts_z_xxy, ts_z_xxz, ts_z_xy, ts_z_xyy, ts_z_xyz, ts_z_xz, ts_z_xzz, ts_z_yy, ts_z_yyy, ts_z_yyz, ts_z_yz, ts_z_yzz, ts_z_zz, ts_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_z_xxx[i] = 2.0 * ts_0_xxx[i] * gfe2_0 + gr_0_xxx[i] * gfe_0 + 2.0 * ts_z_xxx[i] * gfe_0 * gc_z[i] + gr_z_xxx[i] * gc_z[i];

        grr_z_z_xxy[i] = 2.0 * ts_0_xxy[i] * gfe2_0 + gr_0_xxy[i] * gfe_0 + 2.0 * ts_z_xxy[i] * gfe_0 * gc_z[i] + gr_z_xxy[i] * gc_z[i];

        grr_z_z_xxz[i] = 2.0 * ts_0_xxz[i] * gfe2_0 + gr_0_xxz[i] * gfe_0 + 2.0 * ts_z_xx[i] * gfe2_0 + gr_z_xx[i] * gfe_0 + 2.0 * ts_z_xxz[i] * gfe_0 * gc_z[i] + gr_z_xxz[i] * gc_z[i];

        grr_z_z_xyy[i] = 2.0 * ts_0_xyy[i] * gfe2_0 + gr_0_xyy[i] * gfe_0 + 2.0 * ts_z_xyy[i] * gfe_0 * gc_z[i] + gr_z_xyy[i] * gc_z[i];

        grr_z_z_xyz[i] = 2.0 * ts_0_xyz[i] * gfe2_0 + gr_0_xyz[i] * gfe_0 + 2.0 * ts_z_xy[i] * gfe2_0 + gr_z_xy[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe_0 * gc_z[i] + gr_z_xyz[i] * gc_z[i];

        grr_z_z_xzz[i] = 2.0 * ts_0_xzz[i] * gfe2_0 + gr_0_xzz[i] * gfe_0 + 4.0 * ts_z_xz[i] * gfe2_0 + 2.0 * gr_z_xz[i] * gfe_0 + 2.0 * ts_z_xzz[i] * gfe_0 * gc_z[i] + gr_z_xzz[i] * gc_z[i];

        grr_z_z_yyy[i] = 2.0 * ts_0_yyy[i] * gfe2_0 + gr_0_yyy[i] * gfe_0 + 2.0 * ts_z_yyy[i] * gfe_0 * gc_z[i] + gr_z_yyy[i] * gc_z[i];

        grr_z_z_yyz[i] = 2.0 * ts_0_yyz[i] * gfe2_0 + gr_0_yyz[i] * gfe_0 + 2.0 * ts_z_yy[i] * gfe2_0 + gr_z_yy[i] * gfe_0 + 2.0 * ts_z_yyz[i] * gfe_0 * gc_z[i] + gr_z_yyz[i] * gc_z[i];

        grr_z_z_yzz[i] = 2.0 * ts_0_yzz[i] * gfe2_0 + gr_0_yzz[i] * gfe_0 + 4.0 * ts_z_yz[i] * gfe2_0 + 2.0 * gr_z_yz[i] * gfe_0 + 2.0 * ts_z_yzz[i] * gfe_0 * gc_z[i] + gr_z_yzz[i] * gc_z[i];

        grr_z_z_zzz[i] = 2.0 * ts_0_zzz[i] * gfe2_0 + gr_0_zzz[i] * gfe_0 + 6.0 * ts_z_zz[i] * gfe2_0 + 3.0 * gr_z_zz[i] * gfe_0 + 2.0 * ts_z_zzz[i] * gfe_0 * gc_z[i] + gr_z_zzz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

