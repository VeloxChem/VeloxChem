#include "ThreeCenterOverlapGradientPrimRecSG.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_sg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_sg,
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

    auto gs_x_0_xxxx = pbuffer.data(idx_g_sg);

    auto gs_x_0_xxxy = pbuffer.data(idx_g_sg + 1);

    auto gs_x_0_xxxz = pbuffer.data(idx_g_sg + 2);

    auto gs_x_0_xxyy = pbuffer.data(idx_g_sg + 3);

    auto gs_x_0_xxyz = pbuffer.data(idx_g_sg + 4);

    auto gs_x_0_xxzz = pbuffer.data(idx_g_sg + 5);

    auto gs_x_0_xyyy = pbuffer.data(idx_g_sg + 6);

    auto gs_x_0_xyyz = pbuffer.data(idx_g_sg + 7);

    auto gs_x_0_xyzz = pbuffer.data(idx_g_sg + 8);

    auto gs_x_0_xzzz = pbuffer.data(idx_g_sg + 9);

    auto gs_x_0_yyyy = pbuffer.data(idx_g_sg + 10);

    auto gs_x_0_yyyz = pbuffer.data(idx_g_sg + 11);

    auto gs_x_0_yyzz = pbuffer.data(idx_g_sg + 12);

    auto gs_x_0_yzzz = pbuffer.data(idx_g_sg + 13);

    auto gs_x_0_zzzz = pbuffer.data(idx_g_sg + 14);

    auto gs_y_0_xxxx = pbuffer.data(idx_g_sg + 15);

    auto gs_y_0_xxxy = pbuffer.data(idx_g_sg + 16);

    auto gs_y_0_xxxz = pbuffer.data(idx_g_sg + 17);

    auto gs_y_0_xxyy = pbuffer.data(idx_g_sg + 18);

    auto gs_y_0_xxyz = pbuffer.data(idx_g_sg + 19);

    auto gs_y_0_xxzz = pbuffer.data(idx_g_sg + 20);

    auto gs_y_0_xyyy = pbuffer.data(idx_g_sg + 21);

    auto gs_y_0_xyyz = pbuffer.data(idx_g_sg + 22);

    auto gs_y_0_xyzz = pbuffer.data(idx_g_sg + 23);

    auto gs_y_0_xzzz = pbuffer.data(idx_g_sg + 24);

    auto gs_y_0_yyyy = pbuffer.data(idx_g_sg + 25);

    auto gs_y_0_yyyz = pbuffer.data(idx_g_sg + 26);

    auto gs_y_0_yyzz = pbuffer.data(idx_g_sg + 27);

    auto gs_y_0_yzzz = pbuffer.data(idx_g_sg + 28);

    auto gs_y_0_zzzz = pbuffer.data(idx_g_sg + 29);

    auto gs_z_0_xxxx = pbuffer.data(idx_g_sg + 30);

    auto gs_z_0_xxxy = pbuffer.data(idx_g_sg + 31);

    auto gs_z_0_xxxz = pbuffer.data(idx_g_sg + 32);

    auto gs_z_0_xxyy = pbuffer.data(idx_g_sg + 33);

    auto gs_z_0_xxyz = pbuffer.data(idx_g_sg + 34);

    auto gs_z_0_xxzz = pbuffer.data(idx_g_sg + 35);

    auto gs_z_0_xyyy = pbuffer.data(idx_g_sg + 36);

    auto gs_z_0_xyyz = pbuffer.data(idx_g_sg + 37);

    auto gs_z_0_xyzz = pbuffer.data(idx_g_sg + 38);

    auto gs_z_0_xzzz = pbuffer.data(idx_g_sg + 39);

    auto gs_z_0_yyyy = pbuffer.data(idx_g_sg + 40);

    auto gs_z_0_yyyz = pbuffer.data(idx_g_sg + 41);

    auto gs_z_0_yyzz = pbuffer.data(idx_g_sg + 42);

    auto gs_z_0_yzzz = pbuffer.data(idx_g_sg + 43);

    auto gs_z_0_zzzz = pbuffer.data(idx_g_sg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_0_xxxx, gs_x_0_xxxy, gs_x_0_xxxz, gs_x_0_xxyy, gs_x_0_xxyz, gs_x_0_xxzz, gs_x_0_xyyy, gs_x_0_xyyz, gs_x_0_xyzz, gs_x_0_xzzz, gs_x_0_yyyy, gs_x_0_yyyz, gs_x_0_yyzz, gs_x_0_yzzz, gs_x_0_zzzz, gs_y_0_xxxx, gs_y_0_xxxy, gs_y_0_xxxz, gs_y_0_xxyy, gs_y_0_xxyz, gs_y_0_xxzz, gs_y_0_xyyy, gs_y_0_xyyz, gs_y_0_xyzz, gs_y_0_xzzz, gs_y_0_yyyy, gs_y_0_yyyz, gs_y_0_yyzz, gs_y_0_yzzz, gs_y_0_zzzz, gs_z_0_xxxx, gs_z_0_xxxy, gs_z_0_xxxz, gs_z_0_xxyy, gs_z_0_xxyz, gs_z_0_xxzz, gs_z_0_xyyy, gs_z_0_xyyz, gs_z_0_xyzz, gs_z_0_xzzz, gs_z_0_yyyy, gs_z_0_yyyz, gs_z_0_yyzz, gs_z_0_yzzz, gs_z_0_zzzz, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_0_xxxx[i] = 8.0 * ts_0_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxx[i] * gc_x[i] * tce_0;

        gs_x_0_xxxy[i] = 6.0 * ts_0_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxy[i] * gc_x[i] * tce_0;

        gs_x_0_xxxz[i] = 6.0 * ts_0_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxz[i] * gc_x[i] * tce_0;

        gs_x_0_xxyy[i] = 4.0 * ts_0_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyy[i] * gc_x[i] * tce_0;

        gs_x_0_xxyz[i] = 4.0 * ts_0_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyz[i] * gc_x[i] * tce_0;

        gs_x_0_xxzz[i] = 4.0 * ts_0_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxzz[i] * gc_x[i] * tce_0;

        gs_x_0_xyyy[i] = 2.0 * ts_0_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyy[i] * gc_x[i] * tce_0;

        gs_x_0_xyyz[i] = 2.0 * ts_0_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyz[i] * gc_x[i] * tce_0;

        gs_x_0_xyzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyzz[i] * gc_x[i] * tce_0;

        gs_x_0_xzzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xzzz[i] * gc_x[i] * tce_0;

        gs_x_0_yyyy[i] = 2.0 * ts_0_yyyy[i] * gc_x[i] * tce_0;

        gs_x_0_yyyz[i] = 2.0 * ts_0_yyyz[i] * gc_x[i] * tce_0;

        gs_x_0_yyzz[i] = 2.0 * ts_0_yyzz[i] * gc_x[i] * tce_0;

        gs_x_0_yzzz[i] = 2.0 * ts_0_yzzz[i] * gc_x[i] * tce_0;

        gs_x_0_zzzz[i] = 2.0 * ts_0_zzzz[i] * gc_x[i] * tce_0;

        gs_y_0_xxxx[i] = 2.0 * ts_0_xxxx[i] * gc_y[i] * tce_0;

        gs_y_0_xxxy[i] = 2.0 * ts_0_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxy[i] * gc_y[i] * tce_0;

        gs_y_0_xxxz[i] = 2.0 * ts_0_xxxz[i] * gc_y[i] * tce_0;

        gs_y_0_xxyy[i] = 4.0 * ts_0_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyy[i] * gc_y[i] * tce_0;

        gs_y_0_xxyz[i] = 2.0 * ts_0_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyz[i] * gc_y[i] * tce_0;

        gs_y_0_xxzz[i] = 2.0 * ts_0_xxzz[i] * gc_y[i] * tce_0;

        gs_y_0_xyyy[i] = 6.0 * ts_0_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyy[i] * gc_y[i] * tce_0;

        gs_y_0_xyyz[i] = 4.0 * ts_0_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyz[i] * gc_y[i] * tce_0;

        gs_y_0_xyzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyzz[i] * gc_y[i] * tce_0;

        gs_y_0_xzzz[i] = 2.0 * ts_0_xzzz[i] * gc_y[i] * tce_0;

        gs_y_0_yyyy[i] = 8.0 * ts_0_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyy[i] * gc_y[i] * tce_0;

        gs_y_0_yyyz[i] = 6.0 * ts_0_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyz[i] * gc_y[i] * tce_0;

        gs_y_0_yyzz[i] = 4.0 * ts_0_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyzz[i] * gc_y[i] * tce_0;

        gs_y_0_yzzz[i] = 2.0 * ts_0_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yzzz[i] * gc_y[i] * tce_0;

        gs_y_0_zzzz[i] = 2.0 * ts_0_zzzz[i] * gc_y[i] * tce_0;

        gs_z_0_xxxx[i] = 2.0 * ts_0_xxxx[i] * gc_z[i] * tce_0;

        gs_z_0_xxxy[i] = 2.0 * ts_0_xxxy[i] * gc_z[i] * tce_0;

        gs_z_0_xxxz[i] = 2.0 * ts_0_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxz[i] * gc_z[i] * tce_0;

        gs_z_0_xxyy[i] = 2.0 * ts_0_xxyy[i] * gc_z[i] * tce_0;

        gs_z_0_xxyz[i] = 2.0 * ts_0_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyz[i] * gc_z[i] * tce_0;

        gs_z_0_xxzz[i] = 4.0 * ts_0_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxzz[i] * gc_z[i] * tce_0;

        gs_z_0_xyyy[i] = 2.0 * ts_0_xyyy[i] * gc_z[i] * tce_0;

        gs_z_0_xyyz[i] = 2.0 * ts_0_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyz[i] * gc_z[i] * tce_0;

        gs_z_0_xyzz[i] = 4.0 * ts_0_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyzz[i] * gc_z[i] * tce_0;

        gs_z_0_xzzz[i] = 6.0 * ts_0_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xzzz[i] * gc_z[i] * tce_0;

        gs_z_0_yyyy[i] = 2.0 * ts_0_yyyy[i] * gc_z[i] * tce_0;

        gs_z_0_yyyz[i] = 2.0 * ts_0_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyz[i] * gc_z[i] * tce_0;

        gs_z_0_yyzz[i] = 4.0 * ts_0_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyzz[i] * gc_z[i] * tce_0;

        gs_z_0_yzzz[i] = 6.0 * ts_0_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yzzz[i] * gc_z[i] * tce_0;

        gs_z_0_zzzz[i] = 8.0 * ts_0_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_zzzz[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

