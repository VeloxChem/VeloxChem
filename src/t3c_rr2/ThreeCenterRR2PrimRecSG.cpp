#include "ThreeCenterRR2PrimRecSG.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_sg(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_sg,
                  const size_t idx_sf,
                  const size_t idx_g_sf,
                  const size_t idx_sg,
                  const size_t idx_g_sg,
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

    // Set up components of auxiliary buffer : SG

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

    // Set up components of targeted buffer : SG

    auto grr_x_0_xxxx = pbuffer.data(idx_gr_sg);

    auto grr_x_0_xxxy = pbuffer.data(idx_gr_sg + 1);

    auto grr_x_0_xxxz = pbuffer.data(idx_gr_sg + 2);

    auto grr_x_0_xxyy = pbuffer.data(idx_gr_sg + 3);

    auto grr_x_0_xxyz = pbuffer.data(idx_gr_sg + 4);

    auto grr_x_0_xxzz = pbuffer.data(idx_gr_sg + 5);

    auto grr_x_0_xyyy = pbuffer.data(idx_gr_sg + 6);

    auto grr_x_0_xyyz = pbuffer.data(idx_gr_sg + 7);

    auto grr_x_0_xyzz = pbuffer.data(idx_gr_sg + 8);

    auto grr_x_0_xzzz = pbuffer.data(idx_gr_sg + 9);

    auto grr_x_0_yyyy = pbuffer.data(idx_gr_sg + 10);

    auto grr_x_0_yyyz = pbuffer.data(idx_gr_sg + 11);

    auto grr_x_0_yyzz = pbuffer.data(idx_gr_sg + 12);

    auto grr_x_0_yzzz = pbuffer.data(idx_gr_sg + 13);

    auto grr_x_0_zzzz = pbuffer.data(idx_gr_sg + 14);

    auto grr_y_0_xxxx = pbuffer.data(idx_gr_sg + 15);

    auto grr_y_0_xxxy = pbuffer.data(idx_gr_sg + 16);

    auto grr_y_0_xxxz = pbuffer.data(idx_gr_sg + 17);

    auto grr_y_0_xxyy = pbuffer.data(idx_gr_sg + 18);

    auto grr_y_0_xxyz = pbuffer.data(idx_gr_sg + 19);

    auto grr_y_0_xxzz = pbuffer.data(idx_gr_sg + 20);

    auto grr_y_0_xyyy = pbuffer.data(idx_gr_sg + 21);

    auto grr_y_0_xyyz = pbuffer.data(idx_gr_sg + 22);

    auto grr_y_0_xyzz = pbuffer.data(idx_gr_sg + 23);

    auto grr_y_0_xzzz = pbuffer.data(idx_gr_sg + 24);

    auto grr_y_0_yyyy = pbuffer.data(idx_gr_sg + 25);

    auto grr_y_0_yyyz = pbuffer.data(idx_gr_sg + 26);

    auto grr_y_0_yyzz = pbuffer.data(idx_gr_sg + 27);

    auto grr_y_0_yzzz = pbuffer.data(idx_gr_sg + 28);

    auto grr_y_0_zzzz = pbuffer.data(idx_gr_sg + 29);

    auto grr_z_0_xxxx = pbuffer.data(idx_gr_sg + 30);

    auto grr_z_0_xxxy = pbuffer.data(idx_gr_sg + 31);

    auto grr_z_0_xxxz = pbuffer.data(idx_gr_sg + 32);

    auto grr_z_0_xxyy = pbuffer.data(idx_gr_sg + 33);

    auto grr_z_0_xxyz = pbuffer.data(idx_gr_sg + 34);

    auto grr_z_0_xxzz = pbuffer.data(idx_gr_sg + 35);

    auto grr_z_0_xyyy = pbuffer.data(idx_gr_sg + 36);

    auto grr_z_0_xyyz = pbuffer.data(idx_gr_sg + 37);

    auto grr_z_0_xyzz = pbuffer.data(idx_gr_sg + 38);

    auto grr_z_0_xzzz = pbuffer.data(idx_gr_sg + 39);

    auto grr_z_0_yyyy = pbuffer.data(idx_gr_sg + 40);

    auto grr_z_0_yyyz = pbuffer.data(idx_gr_sg + 41);

    auto grr_z_0_yyzz = pbuffer.data(idx_gr_sg + 42);

    auto grr_z_0_yzzz = pbuffer.data(idx_gr_sg + 43);

    auto grr_z_0_zzzz = pbuffer.data(idx_gr_sg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxx, gr_0_xxxx, gr_0_xxxy, gr_0_xxxz, gr_0_xxy, gr_0_xxyy, gr_0_xxyz, gr_0_xxz, gr_0_xxzz, gr_0_xyy, gr_0_xyyy, gr_0_xyyz, gr_0_xyz, gr_0_xyzz, gr_0_xzz, gr_0_xzzz, gr_0_yyy, gr_0_yyyy, gr_0_yyyz, gr_0_yyz, gr_0_yyzz, gr_0_yzz, gr_0_yzzz, gr_0_zzz, gr_0_zzzz, grr_x_0_xxxx, grr_x_0_xxxy, grr_x_0_xxxz, grr_x_0_xxyy, grr_x_0_xxyz, grr_x_0_xxzz, grr_x_0_xyyy, grr_x_0_xyyz, grr_x_0_xyzz, grr_x_0_xzzz, grr_x_0_yyyy, grr_x_0_yyyz, grr_x_0_yyzz, grr_x_0_yzzz, grr_x_0_zzzz, grr_y_0_xxxx, grr_y_0_xxxy, grr_y_0_xxxz, grr_y_0_xxyy, grr_y_0_xxyz, grr_y_0_xxzz, grr_y_0_xyyy, grr_y_0_xyyz, grr_y_0_xyzz, grr_y_0_xzzz, grr_y_0_yyyy, grr_y_0_yyyz, grr_y_0_yyzz, grr_y_0_yzzz, grr_y_0_zzzz, grr_z_0_xxxx, grr_z_0_xxxy, grr_z_0_xxxz, grr_z_0_xxyy, grr_z_0_xxyz, grr_z_0_xxzz, grr_z_0_xyyy, grr_z_0_xyyz, grr_z_0_xyzz, grr_z_0_xzzz, grr_z_0_yyyy, grr_z_0_yyyz, grr_z_0_yyzz, grr_z_0_yzzz, grr_z_0_zzzz, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_0_xxxx[i] = 4.0 * ts_0_xxx[i] * gfe2_0 + 4.0 * gr_0_xxx[i] * gfe_0 + ts_0_xxxx[i] * gfe_0 * gc_x[i] + gr_0_xxxx[i] * gc_x[i];

        grr_x_0_xxxy[i] = 3.0 * ts_0_xxy[i] * gfe2_0 + 3.0 * gr_0_xxy[i] * gfe_0 + ts_0_xxxy[i] * gfe_0 * gc_x[i] + gr_0_xxxy[i] * gc_x[i];

        grr_x_0_xxxz[i] = 3.0 * ts_0_xxz[i] * gfe2_0 + 3.0 * gr_0_xxz[i] * gfe_0 + ts_0_xxxz[i] * gfe_0 * gc_x[i] + gr_0_xxxz[i] * gc_x[i];

        grr_x_0_xxyy[i] = 2.0 * ts_0_xyy[i] * gfe2_0 + 2.0 * gr_0_xyy[i] * gfe_0 + ts_0_xxyy[i] * gfe_0 * gc_x[i] + gr_0_xxyy[i] * gc_x[i];

        grr_x_0_xxyz[i] = 2.0 * ts_0_xyz[i] * gfe2_0 + 2.0 * gr_0_xyz[i] * gfe_0 + ts_0_xxyz[i] * gfe_0 * gc_x[i] + gr_0_xxyz[i] * gc_x[i];

        grr_x_0_xxzz[i] = 2.0 * ts_0_xzz[i] * gfe2_0 + 2.0 * gr_0_xzz[i] * gfe_0 + ts_0_xxzz[i] * gfe_0 * gc_x[i] + gr_0_xxzz[i] * gc_x[i];

        grr_x_0_xyyy[i] = ts_0_yyy[i] * gfe2_0 + gr_0_yyy[i] * gfe_0 + ts_0_xyyy[i] * gfe_0 * gc_x[i] + gr_0_xyyy[i] * gc_x[i];

        grr_x_0_xyyz[i] = ts_0_yyz[i] * gfe2_0 + gr_0_yyz[i] * gfe_0 + ts_0_xyyz[i] * gfe_0 * gc_x[i] + gr_0_xyyz[i] * gc_x[i];

        grr_x_0_xyzz[i] = ts_0_yzz[i] * gfe2_0 + gr_0_yzz[i] * gfe_0 + ts_0_xyzz[i] * gfe_0 * gc_x[i] + gr_0_xyzz[i] * gc_x[i];

        grr_x_0_xzzz[i] = ts_0_zzz[i] * gfe2_0 + gr_0_zzz[i] * gfe_0 + ts_0_xzzz[i] * gfe_0 * gc_x[i] + gr_0_xzzz[i] * gc_x[i];

        grr_x_0_yyyy[i] = ts_0_yyyy[i] * gfe_0 * gc_x[i] + gr_0_yyyy[i] * gc_x[i];

        grr_x_0_yyyz[i] = ts_0_yyyz[i] * gfe_0 * gc_x[i] + gr_0_yyyz[i] * gc_x[i];

        grr_x_0_yyzz[i] = ts_0_yyzz[i] * gfe_0 * gc_x[i] + gr_0_yyzz[i] * gc_x[i];

        grr_x_0_yzzz[i] = ts_0_yzzz[i] * gfe_0 * gc_x[i] + gr_0_yzzz[i] * gc_x[i];

        grr_x_0_zzzz[i] = ts_0_zzzz[i] * gfe_0 * gc_x[i] + gr_0_zzzz[i] * gc_x[i];

        grr_y_0_xxxx[i] = ts_0_xxxx[i] * gfe_0 * gc_y[i] + gr_0_xxxx[i] * gc_y[i];

        grr_y_0_xxxy[i] = ts_0_xxx[i] * gfe2_0 + gr_0_xxx[i] * gfe_0 + ts_0_xxxy[i] * gfe_0 * gc_y[i] + gr_0_xxxy[i] * gc_y[i];

        grr_y_0_xxxz[i] = ts_0_xxxz[i] * gfe_0 * gc_y[i] + gr_0_xxxz[i] * gc_y[i];

        grr_y_0_xxyy[i] = 2.0 * ts_0_xxy[i] * gfe2_0 + 2.0 * gr_0_xxy[i] * gfe_0 + ts_0_xxyy[i] * gfe_0 * gc_y[i] + gr_0_xxyy[i] * gc_y[i];

        grr_y_0_xxyz[i] = ts_0_xxz[i] * gfe2_0 + gr_0_xxz[i] * gfe_0 + ts_0_xxyz[i] * gfe_0 * gc_y[i] + gr_0_xxyz[i] * gc_y[i];

        grr_y_0_xxzz[i] = ts_0_xxzz[i] * gfe_0 * gc_y[i] + gr_0_xxzz[i] * gc_y[i];

        grr_y_0_xyyy[i] = 3.0 * ts_0_xyy[i] * gfe2_0 + 3.0 * gr_0_xyy[i] * gfe_0 + ts_0_xyyy[i] * gfe_0 * gc_y[i] + gr_0_xyyy[i] * gc_y[i];

        grr_y_0_xyyz[i] = 2.0 * ts_0_xyz[i] * gfe2_0 + 2.0 * gr_0_xyz[i] * gfe_0 + ts_0_xyyz[i] * gfe_0 * gc_y[i] + gr_0_xyyz[i] * gc_y[i];

        grr_y_0_xyzz[i] = ts_0_xzz[i] * gfe2_0 + gr_0_xzz[i] * gfe_0 + ts_0_xyzz[i] * gfe_0 * gc_y[i] + gr_0_xyzz[i] * gc_y[i];

        grr_y_0_xzzz[i] = ts_0_xzzz[i] * gfe_0 * gc_y[i] + gr_0_xzzz[i] * gc_y[i];

        grr_y_0_yyyy[i] = 4.0 * ts_0_yyy[i] * gfe2_0 + 4.0 * gr_0_yyy[i] * gfe_0 + ts_0_yyyy[i] * gfe_0 * gc_y[i] + gr_0_yyyy[i] * gc_y[i];

        grr_y_0_yyyz[i] = 3.0 * ts_0_yyz[i] * gfe2_0 + 3.0 * gr_0_yyz[i] * gfe_0 + ts_0_yyyz[i] * gfe_0 * gc_y[i] + gr_0_yyyz[i] * gc_y[i];

        grr_y_0_yyzz[i] = 2.0 * ts_0_yzz[i] * gfe2_0 + 2.0 * gr_0_yzz[i] * gfe_0 + ts_0_yyzz[i] * gfe_0 * gc_y[i] + gr_0_yyzz[i] * gc_y[i];

        grr_y_0_yzzz[i] = ts_0_zzz[i] * gfe2_0 + gr_0_zzz[i] * gfe_0 + ts_0_yzzz[i] * gfe_0 * gc_y[i] + gr_0_yzzz[i] * gc_y[i];

        grr_y_0_zzzz[i] = ts_0_zzzz[i] * gfe_0 * gc_y[i] + gr_0_zzzz[i] * gc_y[i];

        grr_z_0_xxxx[i] = ts_0_xxxx[i] * gfe_0 * gc_z[i] + gr_0_xxxx[i] * gc_z[i];

        grr_z_0_xxxy[i] = ts_0_xxxy[i] * gfe_0 * gc_z[i] + gr_0_xxxy[i] * gc_z[i];

        grr_z_0_xxxz[i] = ts_0_xxx[i] * gfe2_0 + gr_0_xxx[i] * gfe_0 + ts_0_xxxz[i] * gfe_0 * gc_z[i] + gr_0_xxxz[i] * gc_z[i];

        grr_z_0_xxyy[i] = ts_0_xxyy[i] * gfe_0 * gc_z[i] + gr_0_xxyy[i] * gc_z[i];

        grr_z_0_xxyz[i] = ts_0_xxy[i] * gfe2_0 + gr_0_xxy[i] * gfe_0 + ts_0_xxyz[i] * gfe_0 * gc_z[i] + gr_0_xxyz[i] * gc_z[i];

        grr_z_0_xxzz[i] = 2.0 * ts_0_xxz[i] * gfe2_0 + 2.0 * gr_0_xxz[i] * gfe_0 + ts_0_xxzz[i] * gfe_0 * gc_z[i] + gr_0_xxzz[i] * gc_z[i];

        grr_z_0_xyyy[i] = ts_0_xyyy[i] * gfe_0 * gc_z[i] + gr_0_xyyy[i] * gc_z[i];

        grr_z_0_xyyz[i] = ts_0_xyy[i] * gfe2_0 + gr_0_xyy[i] * gfe_0 + ts_0_xyyz[i] * gfe_0 * gc_z[i] + gr_0_xyyz[i] * gc_z[i];

        grr_z_0_xyzz[i] = 2.0 * ts_0_xyz[i] * gfe2_0 + 2.0 * gr_0_xyz[i] * gfe_0 + ts_0_xyzz[i] * gfe_0 * gc_z[i] + gr_0_xyzz[i] * gc_z[i];

        grr_z_0_xzzz[i] = 3.0 * ts_0_xzz[i] * gfe2_0 + 3.0 * gr_0_xzz[i] * gfe_0 + ts_0_xzzz[i] * gfe_0 * gc_z[i] + gr_0_xzzz[i] * gc_z[i];

        grr_z_0_yyyy[i] = ts_0_yyyy[i] * gfe_0 * gc_z[i] + gr_0_yyyy[i] * gc_z[i];

        grr_z_0_yyyz[i] = ts_0_yyy[i] * gfe2_0 + gr_0_yyy[i] * gfe_0 + ts_0_yyyz[i] * gfe_0 * gc_z[i] + gr_0_yyyz[i] * gc_z[i];

        grr_z_0_yyzz[i] = 2.0 * ts_0_yyz[i] * gfe2_0 + 2.0 * gr_0_yyz[i] * gfe_0 + ts_0_yyzz[i] * gfe_0 * gc_z[i] + gr_0_yyzz[i] * gc_z[i];

        grr_z_0_yzzz[i] = 3.0 * ts_0_yzz[i] * gfe2_0 + 3.0 * gr_0_yzz[i] * gfe_0 + ts_0_yzzz[i] * gfe_0 * gc_z[i] + gr_0_yzzz[i] * gc_z[i];

        grr_z_0_zzzz[i] = 4.0 * ts_0_zzz[i] * gfe2_0 + 4.0 * gr_0_zzz[i] * gfe_0 + ts_0_zzzz[i] * gfe_0 * gc_z[i] + gr_0_zzzz[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

