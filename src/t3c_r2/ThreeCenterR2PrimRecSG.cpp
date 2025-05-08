#include "ThreeCenterR2PrimRecSG.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_sg(CSimdArray<double>& pbuffer, 
                const size_t idx_g_sg,
                const size_t idx_sd,
                const size_t idx_sf,
                const size_t idx_sg,
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

    // Set up components of targeted buffer : SG

    auto gr_0_xxxx = pbuffer.data(idx_g_sg);

    auto gr_0_xxxy = pbuffer.data(idx_g_sg + 1);

    auto gr_0_xxxz = pbuffer.data(idx_g_sg + 2);

    auto gr_0_xxyy = pbuffer.data(idx_g_sg + 3);

    auto gr_0_xxyz = pbuffer.data(idx_g_sg + 4);

    auto gr_0_xxzz = pbuffer.data(idx_g_sg + 5);

    auto gr_0_xyyy = pbuffer.data(idx_g_sg + 6);

    auto gr_0_xyyz = pbuffer.data(idx_g_sg + 7);

    auto gr_0_xyzz = pbuffer.data(idx_g_sg + 8);

    auto gr_0_xzzz = pbuffer.data(idx_g_sg + 9);

    auto gr_0_yyyy = pbuffer.data(idx_g_sg + 10);

    auto gr_0_yyyz = pbuffer.data(idx_g_sg + 11);

    auto gr_0_yyzz = pbuffer.data(idx_g_sg + 12);

    auto gr_0_yzzz = pbuffer.data(idx_g_sg + 13);

    auto gr_0_zzzz = pbuffer.data(idx_g_sg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxxx, gr_0_xxxy, gr_0_xxxz, gr_0_xxyy, gr_0_xxyz, gr_0_xxzz, gr_0_xyyy, gr_0_xyyz, gr_0_xyzz, gr_0_xzzz, gr_0_yyyy, gr_0_yyyz, gr_0_yyzz, gr_0_yzzz, gr_0_zzzz, ts_0_xx, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xy, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xz, ts_0_xzz, ts_0_xzzz, ts_0_yy, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yz, ts_0_yzz, ts_0_yzzz, ts_0_zz, ts_0_zzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_0_xxxx[i] = 12.0 * ts_0_xx[i] * gfe2_0 + 8.0 * ts_0_xxx[i] * gfe_0 * gc_x[i] + 3.0 * ts_0_xxxx[i] * gfe_0 + ts_0_xxxx[i] * rgc2_0;

        gr_0_xxxy[i] = 6.0 * ts_0_xy[i] * gfe2_0 + 6.0 * ts_0_xxy[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xxx[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_xxxy[i] * gfe_0 + ts_0_xxxy[i] * rgc2_0;

        gr_0_xxxz[i] = 6.0 * ts_0_xz[i] * gfe2_0 + 6.0 * ts_0_xxz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xxx[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xxxz[i] * gfe_0 + ts_0_xxxz[i] * rgc2_0;

        gr_0_xxyy[i] = 2.0 * ts_0_yy[i] * gfe2_0 + 4.0 * ts_0_xyy[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xx[i] * gfe2_0 + 4.0 * ts_0_xxy[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_xxyy[i] * gfe_0 + ts_0_xxyy[i] * rgc2_0;

        gr_0_xxyz[i] = 2.0 * ts_0_yz[i] * gfe2_0 + 4.0 * ts_0_xyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xxz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_xxy[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xxyz[i] * gfe_0 + ts_0_xxyz[i] * rgc2_0;

        gr_0_xxzz[i] = 2.0 * ts_0_zz[i] * gfe2_0 + 4.0 * ts_0_xzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xx[i] * gfe2_0 + 4.0 * ts_0_xxz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xxzz[i] * gfe_0 + ts_0_xxzz[i] * rgc2_0;

        gr_0_xyyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * gc_x[i] + 6.0 * ts_0_xy[i] * gfe2_0 + 6.0 * ts_0_xyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_xyyy[i] * gfe_0 + ts_0_xyyy[i] * rgc2_0;

        gr_0_xyyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xz[i] * gfe2_0 + 4.0 * ts_0_xyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_xyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xyyz[i] * gfe_0 + ts_0_xyyz[i] * rgc2_0;

        gr_0_xyzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_xzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_xy[i] * gfe2_0 + 4.0 * ts_0_xyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xyzz[i] * gfe_0 + ts_0_xyzz[i] * rgc2_0;

        gr_0_xzzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * gc_x[i] + 6.0 * ts_0_xz[i] * gfe2_0 + 6.0 * ts_0_xzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xzzz[i] * gfe_0 + ts_0_xzzz[i] * rgc2_0;

        gr_0_yyyy[i] = 12.0 * ts_0_yy[i] * gfe2_0 + 8.0 * ts_0_yyy[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_yyyy[i] * gfe_0 + ts_0_yyyy[i] * rgc2_0;

        gr_0_yyyz[i] = 6.0 * ts_0_yz[i] * gfe2_0 + 6.0 * ts_0_yyz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_yyy[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_yyyz[i] * gfe_0 + ts_0_yyyz[i] * rgc2_0;

        gr_0_yyzz[i] = 2.0 * ts_0_zz[i] * gfe2_0 + 4.0 * ts_0_yzz[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_yy[i] * gfe2_0 + 4.0 * ts_0_yyz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_yyzz[i] * gfe_0 + ts_0_yyzz[i] * rgc2_0;

        gr_0_yzzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * gc_y[i] + 6.0 * ts_0_yz[i] * gfe2_0 + 6.0 * ts_0_yzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_yzzz[i] * gfe_0 + ts_0_yzzz[i] * rgc2_0;

        gr_0_zzzz[i] = 12.0 * ts_0_zz[i] * gfe2_0 + 8.0 * ts_0_zzz[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_zzzz[i] * gfe_0 + ts_0_zzzz[i] * rgc2_0;
    }
}

} // t3r2rec namespace

