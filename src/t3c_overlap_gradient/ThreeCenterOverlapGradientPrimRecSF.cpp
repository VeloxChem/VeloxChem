#include "ThreeCenterOverlapGradientPrimRecSF.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_sf(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_sf,
                              const size_t idx_sd,
                              const size_t idx_sf,
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

    // Set up components of targeted buffer : SF

    auto gs_x_0_xxx = pbuffer.data(idx_g_sf);

    auto gs_x_0_xxy = pbuffer.data(idx_g_sf + 1);

    auto gs_x_0_xxz = pbuffer.data(idx_g_sf + 2);

    auto gs_x_0_xyy = pbuffer.data(idx_g_sf + 3);

    auto gs_x_0_xyz = pbuffer.data(idx_g_sf + 4);

    auto gs_x_0_xzz = pbuffer.data(idx_g_sf + 5);

    auto gs_x_0_yyy = pbuffer.data(idx_g_sf + 6);

    auto gs_x_0_yyz = pbuffer.data(idx_g_sf + 7);

    auto gs_x_0_yzz = pbuffer.data(idx_g_sf + 8);

    auto gs_x_0_zzz = pbuffer.data(idx_g_sf + 9);

    auto gs_y_0_xxx = pbuffer.data(idx_g_sf + 10);

    auto gs_y_0_xxy = pbuffer.data(idx_g_sf + 11);

    auto gs_y_0_xxz = pbuffer.data(idx_g_sf + 12);

    auto gs_y_0_xyy = pbuffer.data(idx_g_sf + 13);

    auto gs_y_0_xyz = pbuffer.data(idx_g_sf + 14);

    auto gs_y_0_xzz = pbuffer.data(idx_g_sf + 15);

    auto gs_y_0_yyy = pbuffer.data(idx_g_sf + 16);

    auto gs_y_0_yyz = pbuffer.data(idx_g_sf + 17);

    auto gs_y_0_yzz = pbuffer.data(idx_g_sf + 18);

    auto gs_y_0_zzz = pbuffer.data(idx_g_sf + 19);

    auto gs_z_0_xxx = pbuffer.data(idx_g_sf + 20);

    auto gs_z_0_xxy = pbuffer.data(idx_g_sf + 21);

    auto gs_z_0_xxz = pbuffer.data(idx_g_sf + 22);

    auto gs_z_0_xyy = pbuffer.data(idx_g_sf + 23);

    auto gs_z_0_xyz = pbuffer.data(idx_g_sf + 24);

    auto gs_z_0_xzz = pbuffer.data(idx_g_sf + 25);

    auto gs_z_0_yyy = pbuffer.data(idx_g_sf + 26);

    auto gs_z_0_yyz = pbuffer.data(idx_g_sf + 27);

    auto gs_z_0_yzz = pbuffer.data(idx_g_sf + 28);

    auto gs_z_0_zzz = pbuffer.data(idx_g_sf + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_0_xxx, gs_x_0_xxy, gs_x_0_xxz, gs_x_0_xyy, gs_x_0_xyz, gs_x_0_xzz, gs_x_0_yyy, gs_x_0_yyz, gs_x_0_yzz, gs_x_0_zzz, gs_y_0_xxx, gs_y_0_xxy, gs_y_0_xxz, gs_y_0_xyy, gs_y_0_xyz, gs_y_0_xzz, gs_y_0_yyy, gs_y_0_yyz, gs_y_0_yzz, gs_y_0_zzz, gs_z_0_xxx, gs_z_0_xxy, gs_z_0_xxz, gs_z_0_xyy, gs_z_0_xyz, gs_z_0_xzz, gs_z_0_yyy, gs_z_0_yyz, gs_z_0_yzz, gs_z_0_zzz, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_0_xxx[i] = 6.0 * ts_0_xx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxx[i] * gc_x[i] * tce_0;

        gs_x_0_xxy[i] = 4.0 * ts_0_xy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxy[i] * gc_x[i] * tce_0;

        gs_x_0_xxz[i] = 4.0 * ts_0_xz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxz[i] * gc_x[i] * tce_0;

        gs_x_0_xyy[i] = 2.0 * ts_0_yy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyy[i] * gc_x[i] * tce_0;

        gs_x_0_xyz[i] = 2.0 * ts_0_yz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyz[i] * gc_x[i] * tce_0;

        gs_x_0_xzz[i] = 2.0 * ts_0_zz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xzz[i] * gc_x[i] * tce_0;

        gs_x_0_yyy[i] = 2.0 * ts_0_yyy[i] * gc_x[i] * tce_0;

        gs_x_0_yyz[i] = 2.0 * ts_0_yyz[i] * gc_x[i] * tce_0;

        gs_x_0_yzz[i] = 2.0 * ts_0_yzz[i] * gc_x[i] * tce_0;

        gs_x_0_zzz[i] = 2.0 * ts_0_zzz[i] * gc_x[i] * tce_0;

        gs_y_0_xxx[i] = 2.0 * ts_0_xxx[i] * gc_y[i] * tce_0;

        gs_y_0_xxy[i] = 2.0 * ts_0_xx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxy[i] * gc_y[i] * tce_0;

        gs_y_0_xxz[i] = 2.0 * ts_0_xxz[i] * gc_y[i] * tce_0;

        gs_y_0_xyy[i] = 4.0 * ts_0_xy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyy[i] * gc_y[i] * tce_0;

        gs_y_0_xyz[i] = 2.0 * ts_0_xz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyz[i] * gc_y[i] * tce_0;

        gs_y_0_xzz[i] = 2.0 * ts_0_xzz[i] * gc_y[i] * tce_0;

        gs_y_0_yyy[i] = 6.0 * ts_0_yy[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyy[i] * gc_y[i] * tce_0;

        gs_y_0_yyz[i] = 4.0 * ts_0_yz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyz[i] * gc_y[i] * tce_0;

        gs_y_0_yzz[i] = 2.0 * ts_0_zz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yzz[i] * gc_y[i] * tce_0;

        gs_y_0_zzz[i] = 2.0 * ts_0_zzz[i] * gc_y[i] * tce_0;

        gs_z_0_xxx[i] = 2.0 * ts_0_xxx[i] * gc_z[i] * tce_0;

        gs_z_0_xxy[i] = 2.0 * ts_0_xxy[i] * gc_z[i] * tce_0;

        gs_z_0_xxz[i] = 2.0 * ts_0_xx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxz[i] * gc_z[i] * tce_0;

        gs_z_0_xyy[i] = 2.0 * ts_0_xyy[i] * gc_z[i] * tce_0;

        gs_z_0_xyz[i] = 2.0 * ts_0_xy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyz[i] * gc_z[i] * tce_0;

        gs_z_0_xzz[i] = 4.0 * ts_0_xz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xzz[i] * gc_z[i] * tce_0;

        gs_z_0_yyy[i] = 2.0 * ts_0_yyy[i] * gc_z[i] * tce_0;

        gs_z_0_yyz[i] = 2.0 * ts_0_yy[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyz[i] * gc_z[i] * tce_0;

        gs_z_0_yzz[i] = 4.0 * ts_0_yz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yzz[i] * gc_z[i] * tce_0;

        gs_z_0_zzz[i] = 6.0 * ts_0_zz[i] * gfe_0 * tce_0 + 2.0 * ts_0_zzz[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

