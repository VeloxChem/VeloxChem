#include "ThreeCenterRR2PrimRecSF.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_sf(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_sf,
                  const size_t idx_sd,
                  const size_t idx_g_sd,
                  const size_t idx_sf,
                  const size_t idx_g_sf,
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

    // Set up components of auxiliary buffer : SD

    auto gr_0_xx = pbuffer.data(idx_g_sd);

    auto gr_0_xy = pbuffer.data(idx_g_sd + 1);

    auto gr_0_xz = pbuffer.data(idx_g_sd + 2);

    auto gr_0_yy = pbuffer.data(idx_g_sd + 3);

    auto gr_0_yz = pbuffer.data(idx_g_sd + 4);

    auto gr_0_zz = pbuffer.data(idx_g_sd + 5);

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

    // Set up components of targeted buffer : SF

    auto grr_x_0_xxx = pbuffer.data(idx_gr_sf);

    auto grr_x_0_xxy = pbuffer.data(idx_gr_sf + 1);

    auto grr_x_0_xxz = pbuffer.data(idx_gr_sf + 2);

    auto grr_x_0_xyy = pbuffer.data(idx_gr_sf + 3);

    auto grr_x_0_xyz = pbuffer.data(idx_gr_sf + 4);

    auto grr_x_0_xzz = pbuffer.data(idx_gr_sf + 5);

    auto grr_x_0_yyy = pbuffer.data(idx_gr_sf + 6);

    auto grr_x_0_yyz = pbuffer.data(idx_gr_sf + 7);

    auto grr_x_0_yzz = pbuffer.data(idx_gr_sf + 8);

    auto grr_x_0_zzz = pbuffer.data(idx_gr_sf + 9);

    auto grr_y_0_xxx = pbuffer.data(idx_gr_sf + 10);

    auto grr_y_0_xxy = pbuffer.data(idx_gr_sf + 11);

    auto grr_y_0_xxz = pbuffer.data(idx_gr_sf + 12);

    auto grr_y_0_xyy = pbuffer.data(idx_gr_sf + 13);

    auto grr_y_0_xyz = pbuffer.data(idx_gr_sf + 14);

    auto grr_y_0_xzz = pbuffer.data(idx_gr_sf + 15);

    auto grr_y_0_yyy = pbuffer.data(idx_gr_sf + 16);

    auto grr_y_0_yyz = pbuffer.data(idx_gr_sf + 17);

    auto grr_y_0_yzz = pbuffer.data(idx_gr_sf + 18);

    auto grr_y_0_zzz = pbuffer.data(idx_gr_sf + 19);

    auto grr_z_0_xxx = pbuffer.data(idx_gr_sf + 20);

    auto grr_z_0_xxy = pbuffer.data(idx_gr_sf + 21);

    auto grr_z_0_xxz = pbuffer.data(idx_gr_sf + 22);

    auto grr_z_0_xyy = pbuffer.data(idx_gr_sf + 23);

    auto grr_z_0_xyz = pbuffer.data(idx_gr_sf + 24);

    auto grr_z_0_xzz = pbuffer.data(idx_gr_sf + 25);

    auto grr_z_0_yyy = pbuffer.data(idx_gr_sf + 26);

    auto grr_z_0_yyz = pbuffer.data(idx_gr_sf + 27);

    auto grr_z_0_yzz = pbuffer.data(idx_gr_sf + 28);

    auto grr_z_0_zzz = pbuffer.data(idx_gr_sf + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xx, gr_0_xxx, gr_0_xxy, gr_0_xxz, gr_0_xy, gr_0_xyy, gr_0_xyz, gr_0_xz, gr_0_xzz, gr_0_yy, gr_0_yyy, gr_0_yyz, gr_0_yz, gr_0_yzz, gr_0_zz, gr_0_zzz, grr_x_0_xxx, grr_x_0_xxy, grr_x_0_xxz, grr_x_0_xyy, grr_x_0_xyz, grr_x_0_xzz, grr_x_0_yyy, grr_x_0_yyz, grr_x_0_yzz, grr_x_0_zzz, grr_y_0_xxx, grr_y_0_xxy, grr_y_0_xxz, grr_y_0_xyy, grr_y_0_xyz, grr_y_0_xzz, grr_y_0_yyy, grr_y_0_yyz, grr_y_0_yzz, grr_y_0_zzz, grr_z_0_xxx, grr_z_0_xxy, grr_z_0_xxz, grr_z_0_xyy, grr_z_0_xyz, grr_z_0_xzz, grr_z_0_yyy, grr_z_0_yyz, grr_z_0_yzz, grr_z_0_zzz, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_zz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_0_xxx[i] = 6.0 * ts_0_xx[i] * gfe2_0 + 3.0 * gr_0_xx[i] * gfe_0 + 2.0 * ts_0_xxx[i] * gfe_0 * gc_x[i] + gr_0_xxx[i] * gc_x[i];

        grr_x_0_xxy[i] = 4.0 * ts_0_xy[i] * gfe2_0 + 2.0 * gr_0_xy[i] * gfe_0 + 2.0 * ts_0_xxy[i] * gfe_0 * gc_x[i] + gr_0_xxy[i] * gc_x[i];

        grr_x_0_xxz[i] = 4.0 * ts_0_xz[i] * gfe2_0 + 2.0 * gr_0_xz[i] * gfe_0 + 2.0 * ts_0_xxz[i] * gfe_0 * gc_x[i] + gr_0_xxz[i] * gc_x[i];

        grr_x_0_xyy[i] = 2.0 * ts_0_yy[i] * gfe2_0 + gr_0_yy[i] * gfe_0 + 2.0 * ts_0_xyy[i] * gfe_0 * gc_x[i] + gr_0_xyy[i] * gc_x[i];

        grr_x_0_xyz[i] = 2.0 * ts_0_yz[i] * gfe2_0 + gr_0_yz[i] * gfe_0 + 2.0 * ts_0_xyz[i] * gfe_0 * gc_x[i] + gr_0_xyz[i] * gc_x[i];

        grr_x_0_xzz[i] = 2.0 * ts_0_zz[i] * gfe2_0 + gr_0_zz[i] * gfe_0 + 2.0 * ts_0_xzz[i] * gfe_0 * gc_x[i] + gr_0_xzz[i] * gc_x[i];

        grr_x_0_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * gc_x[i] + gr_0_yyy[i] * gc_x[i];

        grr_x_0_yyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * gc_x[i] + gr_0_yyz[i] * gc_x[i];

        grr_x_0_yzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * gc_x[i] + gr_0_yzz[i] * gc_x[i];

        grr_x_0_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * gc_x[i] + gr_0_zzz[i] * gc_x[i];

        grr_y_0_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * gc_y[i] + gr_0_xxx[i] * gc_y[i];

        grr_y_0_xxy[i] = 2.0 * ts_0_xx[i] * gfe2_0 + gr_0_xx[i] * gfe_0 + 2.0 * ts_0_xxy[i] * gfe_0 * gc_y[i] + gr_0_xxy[i] * gc_y[i];

        grr_y_0_xxz[i] = 2.0 * ts_0_xxz[i] * gfe_0 * gc_y[i] + gr_0_xxz[i] * gc_y[i];

        grr_y_0_xyy[i] = 4.0 * ts_0_xy[i] * gfe2_0 + 2.0 * gr_0_xy[i] * gfe_0 + 2.0 * ts_0_xyy[i] * gfe_0 * gc_y[i] + gr_0_xyy[i] * gc_y[i];

        grr_y_0_xyz[i] = 2.0 * ts_0_xz[i] * gfe2_0 + gr_0_xz[i] * gfe_0 + 2.0 * ts_0_xyz[i] * gfe_0 * gc_y[i] + gr_0_xyz[i] * gc_y[i];

        grr_y_0_xzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 * gc_y[i] + gr_0_xzz[i] * gc_y[i];

        grr_y_0_yyy[i] = 6.0 * ts_0_yy[i] * gfe2_0 + 3.0 * gr_0_yy[i] * gfe_0 + 2.0 * ts_0_yyy[i] * gfe_0 * gc_y[i] + gr_0_yyy[i] * gc_y[i];

        grr_y_0_yyz[i] = 4.0 * ts_0_yz[i] * gfe2_0 + 2.0 * gr_0_yz[i] * gfe_0 + 2.0 * ts_0_yyz[i] * gfe_0 * gc_y[i] + gr_0_yyz[i] * gc_y[i];

        grr_y_0_yzz[i] = 2.0 * ts_0_zz[i] * gfe2_0 + gr_0_zz[i] * gfe_0 + 2.0 * ts_0_yzz[i] * gfe_0 * gc_y[i] + gr_0_yzz[i] * gc_y[i];

        grr_y_0_zzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * gc_y[i] + gr_0_zzz[i] * gc_y[i];

        grr_z_0_xxx[i] = 2.0 * ts_0_xxx[i] * gfe_0 * gc_z[i] + gr_0_xxx[i] * gc_z[i];

        grr_z_0_xxy[i] = 2.0 * ts_0_xxy[i] * gfe_0 * gc_z[i] + gr_0_xxy[i] * gc_z[i];

        grr_z_0_xxz[i] = 2.0 * ts_0_xx[i] * gfe2_0 + gr_0_xx[i] * gfe_0 + 2.0 * ts_0_xxz[i] * gfe_0 * gc_z[i] + gr_0_xxz[i] * gc_z[i];

        grr_z_0_xyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 * gc_z[i] + gr_0_xyy[i] * gc_z[i];

        grr_z_0_xyz[i] = 2.0 * ts_0_xy[i] * gfe2_0 + gr_0_xy[i] * gfe_0 + 2.0 * ts_0_xyz[i] * gfe_0 * gc_z[i] + gr_0_xyz[i] * gc_z[i];

        grr_z_0_xzz[i] = 4.0 * ts_0_xz[i] * gfe2_0 + 2.0 * gr_0_xz[i] * gfe_0 + 2.0 * ts_0_xzz[i] * gfe_0 * gc_z[i] + gr_0_xzz[i] * gc_z[i];

        grr_z_0_yyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * gc_z[i] + gr_0_yyy[i] * gc_z[i];

        grr_z_0_yyz[i] = 2.0 * ts_0_yy[i] * gfe2_0 + gr_0_yy[i] * gfe_0 + 2.0 * ts_0_yyz[i] * gfe_0 * gc_z[i] + gr_0_yyz[i] * gc_z[i];

        grr_z_0_yzz[i] = 4.0 * ts_0_yz[i] * gfe2_0 + 2.0 * gr_0_yz[i] * gfe_0 + 2.0 * ts_0_yzz[i] * gfe_0 * gc_z[i] + gr_0_yzz[i] * gc_z[i];

        grr_z_0_zzz[i] = 6.0 * ts_0_zz[i] * gfe2_0 + 3.0 * gr_0_zz[i] * gfe_0 + 2.0 * ts_0_zzz[i] * gfe_0 * gc_z[i] + gr_0_zzz[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

