#include "ThreeCenterR2PrimRecSF.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_sf(CSimdArray<double>& pbuffer, 
                const size_t idx_g_sf,
                const size_t idx_sp,
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

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

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

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxx, gr_0_xxy, gr_0_xxz, gr_0_xyy, gr_0_xyz, gr_0_xzz, gr_0_yyy, gr_0_yyz, gr_0_yzz, gr_0_zzz, ts_0_x, ts_0_xx, ts_0_xxx, ts_0_xxy, ts_0_xxz, ts_0_xy, ts_0_xyy, ts_0_xyz, ts_0_xz, ts_0_xzz, ts_0_y, ts_0_yy, ts_0_yyy, ts_0_yyz, ts_0_yz, ts_0_yzz, ts_0_z, ts_0_zz, ts_0_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_0_xxx[i] = 6.0 * ts_0_x[i] * gfe_0 + 6.0 * ts_0_xx[i] * gfe_0 * gc_x[i] + 3.0 * ts_0_xxx[i] * gfe_0 + ts_0_xxx[i] * rgc2_0;

        gr_0_xxy[i] = 2.0 * ts_0_y[i] * gfe_0 + 4.0 * ts_0_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xx[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_xxy[i] * gfe_0 + ts_0_xxy[i] * rgc2_0;

        gr_0_xxz[i] = 2.0 * ts_0_z[i] * gfe_0 + 4.0 * ts_0_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xx[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xxz[i] * gfe_0 + ts_0_xxz[i] * rgc2_0;

        gr_0_xyy[i] = 2.0 * ts_0_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_x[i] * gfe_0 + 4.0 * ts_0_xy[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_xyy[i] * gfe_0 + ts_0_xyy[i] * rgc2_0;

        gr_0_xyz[i] = 2.0 * ts_0_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_xy[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xyz[i] * gfe_0 + ts_0_xyz[i] * rgc2_0;

        gr_0_xzz[i] = 2.0 * ts_0_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_x[i] * gfe_0 + 4.0 * ts_0_xz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xzz[i] * gfe_0 + ts_0_xzz[i] * rgc2_0;

        gr_0_yyy[i] = 6.0 * ts_0_y[i] * gfe_0 + 6.0 * ts_0_yy[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_yyy[i] * gfe_0 + ts_0_yyy[i] * rgc2_0;

        gr_0_yyz[i] = 2.0 * ts_0_z[i] * gfe_0 + 4.0 * ts_0_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_yy[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_yyz[i] * gfe_0 + ts_0_yyz[i] * rgc2_0;

        gr_0_yzz[i] = 2.0 * ts_0_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_y[i] * gfe_0 + 4.0 * ts_0_yz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_yzz[i] * gfe_0 + ts_0_yzz[i] * rgc2_0;

        gr_0_zzz[i] = 6.0 * ts_0_z[i] * gfe_0 + 6.0 * ts_0_zz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_zzz[i] * gfe_0 + ts_0_zzz[i] * rgc2_0;
    }
}

} // t3r2rec namespace

