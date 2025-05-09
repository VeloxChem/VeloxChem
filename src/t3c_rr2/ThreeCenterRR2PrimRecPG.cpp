#include "ThreeCenterRR2PrimRecPG.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_pg(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_pg,
                  const size_t idx_sg,
                  const size_t idx_g_sg,
                  const size_t idx_pf,
                  const size_t idx_g_pf,
                  const size_t idx_pg,
                  const size_t idx_g_pg,
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

    // Set up components of auxiliary buffer : PG

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

    // Set up 0-15 components of targeted buffer : PG

    auto grr_x_x_xxxx = pbuffer.data(idx_gr_pg);

    auto grr_x_x_xxxy = pbuffer.data(idx_gr_pg + 1);

    auto grr_x_x_xxxz = pbuffer.data(idx_gr_pg + 2);

    auto grr_x_x_xxyy = pbuffer.data(idx_gr_pg + 3);

    auto grr_x_x_xxyz = pbuffer.data(idx_gr_pg + 4);

    auto grr_x_x_xxzz = pbuffer.data(idx_gr_pg + 5);

    auto grr_x_x_xyyy = pbuffer.data(idx_gr_pg + 6);

    auto grr_x_x_xyyz = pbuffer.data(idx_gr_pg + 7);

    auto grr_x_x_xyzz = pbuffer.data(idx_gr_pg + 8);

    auto grr_x_x_xzzz = pbuffer.data(idx_gr_pg + 9);

    auto grr_x_x_yyyy = pbuffer.data(idx_gr_pg + 10);

    auto grr_x_x_yyyz = pbuffer.data(idx_gr_pg + 11);

    auto grr_x_x_yyzz = pbuffer.data(idx_gr_pg + 12);

    auto grr_x_x_yzzz = pbuffer.data(idx_gr_pg + 13);

    auto grr_x_x_zzzz = pbuffer.data(idx_gr_pg + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxxx, gr_0_xxxy, gr_0_xxxz, gr_0_xxyy, gr_0_xxyz, gr_0_xxzz, gr_0_xyyy, gr_0_xyyz, gr_0_xyzz, gr_0_xzzz, gr_0_yyyy, gr_0_yyyz, gr_0_yyzz, gr_0_yzzz, gr_0_zzzz, gr_x_xxx, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxy, gr_x_xxyy, gr_x_xxyz, gr_x_xxz, gr_x_xxzz, gr_x_xyy, gr_x_xyyy, gr_x_xyyz, gr_x_xyz, gr_x_xyzz, gr_x_xzz, gr_x_xzzz, gr_x_yyy, gr_x_yyyy, gr_x_yyyz, gr_x_yyz, gr_x_yyzz, gr_x_yzz, gr_x_yzzz, gr_x_zzz, gr_x_zzzz, grr_x_x_xxxx, grr_x_x_xxxy, grr_x_x_xxxz, grr_x_x_xxyy, grr_x_x_xxyz, grr_x_x_xxzz, grr_x_x_xyyy, grr_x_x_xyyz, grr_x_x_xyzz, grr_x_x_xzzz, grr_x_x_yyyy, grr_x_x_yyyz, grr_x_x_yyzz, grr_x_x_yzzz, grr_x_x_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_x_xxxx[i] = ts_0_xxxx[i] * gfe2_0 + gr_0_xxxx[i] * gfe_0 + 4.0 * ts_x_xxx[i] * gfe2_0 + 4.0 * gr_x_xxx[i] * gfe_0 + ts_x_xxxx[i] * gfe_0 * gc_x[i] + gr_x_xxxx[i] * gc_x[i];

        grr_x_x_xxxy[i] = ts_0_xxxy[i] * gfe2_0 + gr_0_xxxy[i] * gfe_0 + 3.0 * ts_x_xxy[i] * gfe2_0 + 3.0 * gr_x_xxy[i] * gfe_0 + ts_x_xxxy[i] * gfe_0 * gc_x[i] + gr_x_xxxy[i] * gc_x[i];

        grr_x_x_xxxz[i] = ts_0_xxxz[i] * gfe2_0 + gr_0_xxxz[i] * gfe_0 + 3.0 * ts_x_xxz[i] * gfe2_0 + 3.0 * gr_x_xxz[i] * gfe_0 + ts_x_xxxz[i] * gfe_0 * gc_x[i] + gr_x_xxxz[i] * gc_x[i];

        grr_x_x_xxyy[i] = ts_0_xxyy[i] * gfe2_0 + gr_0_xxyy[i] * gfe_0 + 2.0 * ts_x_xyy[i] * gfe2_0 + 2.0 * gr_x_xyy[i] * gfe_0 + ts_x_xxyy[i] * gfe_0 * gc_x[i] + gr_x_xxyy[i] * gc_x[i];

        grr_x_x_xxyz[i] = ts_0_xxyz[i] * gfe2_0 + gr_0_xxyz[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe2_0 + 2.0 * gr_x_xyz[i] * gfe_0 + ts_x_xxyz[i] * gfe_0 * gc_x[i] + gr_x_xxyz[i] * gc_x[i];

        grr_x_x_xxzz[i] = ts_0_xxzz[i] * gfe2_0 + gr_0_xxzz[i] * gfe_0 + 2.0 * ts_x_xzz[i] * gfe2_0 + 2.0 * gr_x_xzz[i] * gfe_0 + ts_x_xxzz[i] * gfe_0 * gc_x[i] + gr_x_xxzz[i] * gc_x[i];

        grr_x_x_xyyy[i] = ts_0_xyyy[i] * gfe2_0 + gr_0_xyyy[i] * gfe_0 + ts_x_yyy[i] * gfe2_0 + gr_x_yyy[i] * gfe_0 + ts_x_xyyy[i] * gfe_0 * gc_x[i] + gr_x_xyyy[i] * gc_x[i];

        grr_x_x_xyyz[i] = ts_0_xyyz[i] * gfe2_0 + gr_0_xyyz[i] * gfe_0 + ts_x_yyz[i] * gfe2_0 + gr_x_yyz[i] * gfe_0 + ts_x_xyyz[i] * gfe_0 * gc_x[i] + gr_x_xyyz[i] * gc_x[i];

        grr_x_x_xyzz[i] = ts_0_xyzz[i] * gfe2_0 + gr_0_xyzz[i] * gfe_0 + ts_x_yzz[i] * gfe2_0 + gr_x_yzz[i] * gfe_0 + ts_x_xyzz[i] * gfe_0 * gc_x[i] + gr_x_xyzz[i] * gc_x[i];

        grr_x_x_xzzz[i] = ts_0_xzzz[i] * gfe2_0 + gr_0_xzzz[i] * gfe_0 + ts_x_zzz[i] * gfe2_0 + gr_x_zzz[i] * gfe_0 + ts_x_xzzz[i] * gfe_0 * gc_x[i] + gr_x_xzzz[i] * gc_x[i];

        grr_x_x_yyyy[i] = ts_0_yyyy[i] * gfe2_0 + gr_0_yyyy[i] * gfe_0 + ts_x_yyyy[i] * gfe_0 * gc_x[i] + gr_x_yyyy[i] * gc_x[i];

        grr_x_x_yyyz[i] = ts_0_yyyz[i] * gfe2_0 + gr_0_yyyz[i] * gfe_0 + ts_x_yyyz[i] * gfe_0 * gc_x[i] + gr_x_yyyz[i] * gc_x[i];

        grr_x_x_yyzz[i] = ts_0_yyzz[i] * gfe2_0 + gr_0_yyzz[i] * gfe_0 + ts_x_yyzz[i] * gfe_0 * gc_x[i] + gr_x_yyzz[i] * gc_x[i];

        grr_x_x_yzzz[i] = ts_0_yzzz[i] * gfe2_0 + gr_0_yzzz[i] * gfe_0 + ts_x_yzzz[i] * gfe_0 * gc_x[i] + gr_x_yzzz[i] * gc_x[i];

        grr_x_x_zzzz[i] = ts_0_zzzz[i] * gfe2_0 + gr_0_zzzz[i] * gfe_0 + ts_x_zzzz[i] * gfe_0 * gc_x[i] + gr_x_zzzz[i] * gc_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

    auto grr_x_y_xxxx = pbuffer.data(idx_gr_pg + 15);

    auto grr_x_y_xxxy = pbuffer.data(idx_gr_pg + 16);

    auto grr_x_y_xxxz = pbuffer.data(idx_gr_pg + 17);

    auto grr_x_y_xxyy = pbuffer.data(idx_gr_pg + 18);

    auto grr_x_y_xxyz = pbuffer.data(idx_gr_pg + 19);

    auto grr_x_y_xxzz = pbuffer.data(idx_gr_pg + 20);

    auto grr_x_y_xyyy = pbuffer.data(idx_gr_pg + 21);

    auto grr_x_y_xyyz = pbuffer.data(idx_gr_pg + 22);

    auto grr_x_y_xyzz = pbuffer.data(idx_gr_pg + 23);

    auto grr_x_y_xzzz = pbuffer.data(idx_gr_pg + 24);

    auto grr_x_y_yyyy = pbuffer.data(idx_gr_pg + 25);

    auto grr_x_y_yyyz = pbuffer.data(idx_gr_pg + 26);

    auto grr_x_y_yyzz = pbuffer.data(idx_gr_pg + 27);

    auto grr_x_y_yzzz = pbuffer.data(idx_gr_pg + 28);

    auto grr_x_y_zzzz = pbuffer.data(idx_gr_pg + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxx, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxy, gr_y_xxyy, gr_y_xxyz, gr_y_xxz, gr_y_xxzz, gr_y_xyy, gr_y_xyyy, gr_y_xyyz, gr_y_xyz, gr_y_xyzz, gr_y_xzz, gr_y_xzzz, gr_y_yyy, gr_y_yyyy, gr_y_yyyz, gr_y_yyz, gr_y_yyzz, gr_y_yzz, gr_y_yzzz, gr_y_zzz, gr_y_zzzz, grr_x_y_xxxx, grr_x_y_xxxy, grr_x_y_xxxz, grr_x_y_xxyy, grr_x_y_xxyz, grr_x_y_xxzz, grr_x_y_xyyy, grr_x_y_xyyz, grr_x_y_xyzz, grr_x_y_xzzz, grr_x_y_yyyy, grr_x_y_yyyz, grr_x_y_yyzz, grr_x_y_yzzz, grr_x_y_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_y_xxxx[i] = 4.0 * ts_y_xxx[i] * gfe2_0 + 4.0 * gr_y_xxx[i] * gfe_0 + ts_y_xxxx[i] * gfe_0 * gc_x[i] + gr_y_xxxx[i] * gc_x[i];

        grr_x_y_xxxy[i] = 3.0 * ts_y_xxy[i] * gfe2_0 + 3.0 * gr_y_xxy[i] * gfe_0 + ts_y_xxxy[i] * gfe_0 * gc_x[i] + gr_y_xxxy[i] * gc_x[i];

        grr_x_y_xxxz[i] = 3.0 * ts_y_xxz[i] * gfe2_0 + 3.0 * gr_y_xxz[i] * gfe_0 + ts_y_xxxz[i] * gfe_0 * gc_x[i] + gr_y_xxxz[i] * gc_x[i];

        grr_x_y_xxyy[i] = 2.0 * ts_y_xyy[i] * gfe2_0 + 2.0 * gr_y_xyy[i] * gfe_0 + ts_y_xxyy[i] * gfe_0 * gc_x[i] + gr_y_xxyy[i] * gc_x[i];

        grr_x_y_xxyz[i] = 2.0 * ts_y_xyz[i] * gfe2_0 + 2.0 * gr_y_xyz[i] * gfe_0 + ts_y_xxyz[i] * gfe_0 * gc_x[i] + gr_y_xxyz[i] * gc_x[i];

        grr_x_y_xxzz[i] = 2.0 * ts_y_xzz[i] * gfe2_0 + 2.0 * gr_y_xzz[i] * gfe_0 + ts_y_xxzz[i] * gfe_0 * gc_x[i] + gr_y_xxzz[i] * gc_x[i];

        grr_x_y_xyyy[i] = ts_y_yyy[i] * gfe2_0 + gr_y_yyy[i] * gfe_0 + ts_y_xyyy[i] * gfe_0 * gc_x[i] + gr_y_xyyy[i] * gc_x[i];

        grr_x_y_xyyz[i] = ts_y_yyz[i] * gfe2_0 + gr_y_yyz[i] * gfe_0 + ts_y_xyyz[i] * gfe_0 * gc_x[i] + gr_y_xyyz[i] * gc_x[i];

        grr_x_y_xyzz[i] = ts_y_yzz[i] * gfe2_0 + gr_y_yzz[i] * gfe_0 + ts_y_xyzz[i] * gfe_0 * gc_x[i] + gr_y_xyzz[i] * gc_x[i];

        grr_x_y_xzzz[i] = ts_y_zzz[i] * gfe2_0 + gr_y_zzz[i] * gfe_0 + ts_y_xzzz[i] * gfe_0 * gc_x[i] + gr_y_xzzz[i] * gc_x[i];

        grr_x_y_yyyy[i] = ts_y_yyyy[i] * gfe_0 * gc_x[i] + gr_y_yyyy[i] * gc_x[i];

        grr_x_y_yyyz[i] = ts_y_yyyz[i] * gfe_0 * gc_x[i] + gr_y_yyyz[i] * gc_x[i];

        grr_x_y_yyzz[i] = ts_y_yyzz[i] * gfe_0 * gc_x[i] + gr_y_yyzz[i] * gc_x[i];

        grr_x_y_yzzz[i] = ts_y_yzzz[i] * gfe_0 * gc_x[i] + gr_y_yzzz[i] * gc_x[i];

        grr_x_y_zzzz[i] = ts_y_zzzz[i] * gfe_0 * gc_x[i] + gr_y_zzzz[i] * gc_x[i];
    }

    // Set up 30-45 components of targeted buffer : PG

    auto grr_x_z_xxxx = pbuffer.data(idx_gr_pg + 30);

    auto grr_x_z_xxxy = pbuffer.data(idx_gr_pg + 31);

    auto grr_x_z_xxxz = pbuffer.data(idx_gr_pg + 32);

    auto grr_x_z_xxyy = pbuffer.data(idx_gr_pg + 33);

    auto grr_x_z_xxyz = pbuffer.data(idx_gr_pg + 34);

    auto grr_x_z_xxzz = pbuffer.data(idx_gr_pg + 35);

    auto grr_x_z_xyyy = pbuffer.data(idx_gr_pg + 36);

    auto grr_x_z_xyyz = pbuffer.data(idx_gr_pg + 37);

    auto grr_x_z_xyzz = pbuffer.data(idx_gr_pg + 38);

    auto grr_x_z_xzzz = pbuffer.data(idx_gr_pg + 39);

    auto grr_x_z_yyyy = pbuffer.data(idx_gr_pg + 40);

    auto grr_x_z_yyyz = pbuffer.data(idx_gr_pg + 41);

    auto grr_x_z_yyzz = pbuffer.data(idx_gr_pg + 42);

    auto grr_x_z_yzzz = pbuffer.data(idx_gr_pg + 43);

    auto grr_x_z_zzzz = pbuffer.data(idx_gr_pg + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xxx, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxy, gr_z_xxyy, gr_z_xxyz, gr_z_xxz, gr_z_xxzz, gr_z_xyy, gr_z_xyyy, gr_z_xyyz, gr_z_xyz, gr_z_xyzz, gr_z_xzz, gr_z_xzzz, gr_z_yyy, gr_z_yyyy, gr_z_yyyz, gr_z_yyz, gr_z_yyzz, gr_z_yzz, gr_z_yzzz, gr_z_zzz, gr_z_zzzz, grr_x_z_xxxx, grr_x_z_xxxy, grr_x_z_xxxz, grr_x_z_xxyy, grr_x_z_xxyz, grr_x_z_xxzz, grr_x_z_xyyy, grr_x_z_xyyz, grr_x_z_xyzz, grr_x_z_xzzz, grr_x_z_yyyy, grr_x_z_yyyz, grr_x_z_yyzz, grr_x_z_yzzz, grr_x_z_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_z_xxxx[i] = 4.0 * ts_z_xxx[i] * gfe2_0 + 4.0 * gr_z_xxx[i] * gfe_0 + ts_z_xxxx[i] * gfe_0 * gc_x[i] + gr_z_xxxx[i] * gc_x[i];

        grr_x_z_xxxy[i] = 3.0 * ts_z_xxy[i] * gfe2_0 + 3.0 * gr_z_xxy[i] * gfe_0 + ts_z_xxxy[i] * gfe_0 * gc_x[i] + gr_z_xxxy[i] * gc_x[i];

        grr_x_z_xxxz[i] = 3.0 * ts_z_xxz[i] * gfe2_0 + 3.0 * gr_z_xxz[i] * gfe_0 + ts_z_xxxz[i] * gfe_0 * gc_x[i] + gr_z_xxxz[i] * gc_x[i];

        grr_x_z_xxyy[i] = 2.0 * ts_z_xyy[i] * gfe2_0 + 2.0 * gr_z_xyy[i] * gfe_0 + ts_z_xxyy[i] * gfe_0 * gc_x[i] + gr_z_xxyy[i] * gc_x[i];

        grr_x_z_xxyz[i] = 2.0 * ts_z_xyz[i] * gfe2_0 + 2.0 * gr_z_xyz[i] * gfe_0 + ts_z_xxyz[i] * gfe_0 * gc_x[i] + gr_z_xxyz[i] * gc_x[i];

        grr_x_z_xxzz[i] = 2.0 * ts_z_xzz[i] * gfe2_0 + 2.0 * gr_z_xzz[i] * gfe_0 + ts_z_xxzz[i] * gfe_0 * gc_x[i] + gr_z_xxzz[i] * gc_x[i];

        grr_x_z_xyyy[i] = ts_z_yyy[i] * gfe2_0 + gr_z_yyy[i] * gfe_0 + ts_z_xyyy[i] * gfe_0 * gc_x[i] + gr_z_xyyy[i] * gc_x[i];

        grr_x_z_xyyz[i] = ts_z_yyz[i] * gfe2_0 + gr_z_yyz[i] * gfe_0 + ts_z_xyyz[i] * gfe_0 * gc_x[i] + gr_z_xyyz[i] * gc_x[i];

        grr_x_z_xyzz[i] = ts_z_yzz[i] * gfe2_0 + gr_z_yzz[i] * gfe_0 + ts_z_xyzz[i] * gfe_0 * gc_x[i] + gr_z_xyzz[i] * gc_x[i];

        grr_x_z_xzzz[i] = ts_z_zzz[i] * gfe2_0 + gr_z_zzz[i] * gfe_0 + ts_z_xzzz[i] * gfe_0 * gc_x[i] + gr_z_xzzz[i] * gc_x[i];

        grr_x_z_yyyy[i] = ts_z_yyyy[i] * gfe_0 * gc_x[i] + gr_z_yyyy[i] * gc_x[i];

        grr_x_z_yyyz[i] = ts_z_yyyz[i] * gfe_0 * gc_x[i] + gr_z_yyyz[i] * gc_x[i];

        grr_x_z_yyzz[i] = ts_z_yyzz[i] * gfe_0 * gc_x[i] + gr_z_yyzz[i] * gc_x[i];

        grr_x_z_yzzz[i] = ts_z_yzzz[i] * gfe_0 * gc_x[i] + gr_z_yzzz[i] * gc_x[i];

        grr_x_z_zzzz[i] = ts_z_zzzz[i] * gfe_0 * gc_x[i] + gr_z_zzzz[i] * gc_x[i];
    }

    // Set up 45-60 components of targeted buffer : PG

    auto grr_y_x_xxxx = pbuffer.data(idx_gr_pg + 45);

    auto grr_y_x_xxxy = pbuffer.data(idx_gr_pg + 46);

    auto grr_y_x_xxxz = pbuffer.data(idx_gr_pg + 47);

    auto grr_y_x_xxyy = pbuffer.data(idx_gr_pg + 48);

    auto grr_y_x_xxyz = pbuffer.data(idx_gr_pg + 49);

    auto grr_y_x_xxzz = pbuffer.data(idx_gr_pg + 50);

    auto grr_y_x_xyyy = pbuffer.data(idx_gr_pg + 51);

    auto grr_y_x_xyyz = pbuffer.data(idx_gr_pg + 52);

    auto grr_y_x_xyzz = pbuffer.data(idx_gr_pg + 53);

    auto grr_y_x_xzzz = pbuffer.data(idx_gr_pg + 54);

    auto grr_y_x_yyyy = pbuffer.data(idx_gr_pg + 55);

    auto grr_y_x_yyyz = pbuffer.data(idx_gr_pg + 56);

    auto grr_y_x_yyzz = pbuffer.data(idx_gr_pg + 57);

    auto grr_y_x_yzzz = pbuffer.data(idx_gr_pg + 58);

    auto grr_y_x_zzzz = pbuffer.data(idx_gr_pg + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxx, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxy, gr_x_xxyy, gr_x_xxyz, gr_x_xxz, gr_x_xxzz, gr_x_xyy, gr_x_xyyy, gr_x_xyyz, gr_x_xyz, gr_x_xyzz, gr_x_xzz, gr_x_xzzz, gr_x_yyy, gr_x_yyyy, gr_x_yyyz, gr_x_yyz, gr_x_yyzz, gr_x_yzz, gr_x_yzzz, gr_x_zzz, gr_x_zzzz, grr_y_x_xxxx, grr_y_x_xxxy, grr_y_x_xxxz, grr_y_x_xxyy, grr_y_x_xxyz, grr_y_x_xxzz, grr_y_x_xyyy, grr_y_x_xyyz, grr_y_x_xyzz, grr_y_x_xzzz, grr_y_x_yyyy, grr_y_x_yyyz, grr_y_x_yyzz, grr_y_x_yzzz, grr_y_x_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_x_xxxx[i] = ts_x_xxxx[i] * gfe_0 * gc_y[i] + gr_x_xxxx[i] * gc_y[i];

        grr_y_x_xxxy[i] = ts_x_xxx[i] * gfe2_0 + gr_x_xxx[i] * gfe_0 + ts_x_xxxy[i] * gfe_0 * gc_y[i] + gr_x_xxxy[i] * gc_y[i];

        grr_y_x_xxxz[i] = ts_x_xxxz[i] * gfe_0 * gc_y[i] + gr_x_xxxz[i] * gc_y[i];

        grr_y_x_xxyy[i] = 2.0 * ts_x_xxy[i] * gfe2_0 + 2.0 * gr_x_xxy[i] * gfe_0 + ts_x_xxyy[i] * gfe_0 * gc_y[i] + gr_x_xxyy[i] * gc_y[i];

        grr_y_x_xxyz[i] = ts_x_xxz[i] * gfe2_0 + gr_x_xxz[i] * gfe_0 + ts_x_xxyz[i] * gfe_0 * gc_y[i] + gr_x_xxyz[i] * gc_y[i];

        grr_y_x_xxzz[i] = ts_x_xxzz[i] * gfe_0 * gc_y[i] + gr_x_xxzz[i] * gc_y[i];

        grr_y_x_xyyy[i] = 3.0 * ts_x_xyy[i] * gfe2_0 + 3.0 * gr_x_xyy[i] * gfe_0 + ts_x_xyyy[i] * gfe_0 * gc_y[i] + gr_x_xyyy[i] * gc_y[i];

        grr_y_x_xyyz[i] = 2.0 * ts_x_xyz[i] * gfe2_0 + 2.0 * gr_x_xyz[i] * gfe_0 + ts_x_xyyz[i] * gfe_0 * gc_y[i] + gr_x_xyyz[i] * gc_y[i];

        grr_y_x_xyzz[i] = ts_x_xzz[i] * gfe2_0 + gr_x_xzz[i] * gfe_0 + ts_x_xyzz[i] * gfe_0 * gc_y[i] + gr_x_xyzz[i] * gc_y[i];

        grr_y_x_xzzz[i] = ts_x_xzzz[i] * gfe_0 * gc_y[i] + gr_x_xzzz[i] * gc_y[i];

        grr_y_x_yyyy[i] = 4.0 * ts_x_yyy[i] * gfe2_0 + 4.0 * gr_x_yyy[i] * gfe_0 + ts_x_yyyy[i] * gfe_0 * gc_y[i] + gr_x_yyyy[i] * gc_y[i];

        grr_y_x_yyyz[i] = 3.0 * ts_x_yyz[i] * gfe2_0 + 3.0 * gr_x_yyz[i] * gfe_0 + ts_x_yyyz[i] * gfe_0 * gc_y[i] + gr_x_yyyz[i] * gc_y[i];

        grr_y_x_yyzz[i] = 2.0 * ts_x_yzz[i] * gfe2_0 + 2.0 * gr_x_yzz[i] * gfe_0 + ts_x_yyzz[i] * gfe_0 * gc_y[i] + gr_x_yyzz[i] * gc_y[i];

        grr_y_x_yzzz[i] = ts_x_zzz[i] * gfe2_0 + gr_x_zzz[i] * gfe_0 + ts_x_yzzz[i] * gfe_0 * gc_y[i] + gr_x_yzzz[i] * gc_y[i];

        grr_y_x_zzzz[i] = ts_x_zzzz[i] * gfe_0 * gc_y[i] + gr_x_zzzz[i] * gc_y[i];
    }

    // Set up 60-75 components of targeted buffer : PG

    auto grr_y_y_xxxx = pbuffer.data(idx_gr_pg + 60);

    auto grr_y_y_xxxy = pbuffer.data(idx_gr_pg + 61);

    auto grr_y_y_xxxz = pbuffer.data(idx_gr_pg + 62);

    auto grr_y_y_xxyy = pbuffer.data(idx_gr_pg + 63);

    auto grr_y_y_xxyz = pbuffer.data(idx_gr_pg + 64);

    auto grr_y_y_xxzz = pbuffer.data(idx_gr_pg + 65);

    auto grr_y_y_xyyy = pbuffer.data(idx_gr_pg + 66);

    auto grr_y_y_xyyz = pbuffer.data(idx_gr_pg + 67);

    auto grr_y_y_xyzz = pbuffer.data(idx_gr_pg + 68);

    auto grr_y_y_xzzz = pbuffer.data(idx_gr_pg + 69);

    auto grr_y_y_yyyy = pbuffer.data(idx_gr_pg + 70);

    auto grr_y_y_yyyz = pbuffer.data(idx_gr_pg + 71);

    auto grr_y_y_yyzz = pbuffer.data(idx_gr_pg + 72);

    auto grr_y_y_yzzz = pbuffer.data(idx_gr_pg + 73);

    auto grr_y_y_zzzz = pbuffer.data(idx_gr_pg + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxxx, gr_0_xxxy, gr_0_xxxz, gr_0_xxyy, gr_0_xxyz, gr_0_xxzz, gr_0_xyyy, gr_0_xyyz, gr_0_xyzz, gr_0_xzzz, gr_0_yyyy, gr_0_yyyz, gr_0_yyzz, gr_0_yzzz, gr_0_zzzz, gr_y_xxx, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxy, gr_y_xxyy, gr_y_xxyz, gr_y_xxz, gr_y_xxzz, gr_y_xyy, gr_y_xyyy, gr_y_xyyz, gr_y_xyz, gr_y_xyzz, gr_y_xzz, gr_y_xzzz, gr_y_yyy, gr_y_yyyy, gr_y_yyyz, gr_y_yyz, gr_y_yyzz, gr_y_yzz, gr_y_yzzz, gr_y_zzz, gr_y_zzzz, grr_y_y_xxxx, grr_y_y_xxxy, grr_y_y_xxxz, grr_y_y_xxyy, grr_y_y_xxyz, grr_y_y_xxzz, grr_y_y_xyyy, grr_y_y_xyyz, grr_y_y_xyzz, grr_y_y_xzzz, grr_y_y_yyyy, grr_y_y_yyyz, grr_y_y_yyzz, grr_y_y_yzzz, grr_y_y_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_y_xxxx[i] = ts_0_xxxx[i] * gfe2_0 + gr_0_xxxx[i] * gfe_0 + ts_y_xxxx[i] * gfe_0 * gc_y[i] + gr_y_xxxx[i] * gc_y[i];

        grr_y_y_xxxy[i] = ts_0_xxxy[i] * gfe2_0 + gr_0_xxxy[i] * gfe_0 + ts_y_xxx[i] * gfe2_0 + gr_y_xxx[i] * gfe_0 + ts_y_xxxy[i] * gfe_0 * gc_y[i] + gr_y_xxxy[i] * gc_y[i];

        grr_y_y_xxxz[i] = ts_0_xxxz[i] * gfe2_0 + gr_0_xxxz[i] * gfe_0 + ts_y_xxxz[i] * gfe_0 * gc_y[i] + gr_y_xxxz[i] * gc_y[i];

        grr_y_y_xxyy[i] = ts_0_xxyy[i] * gfe2_0 + gr_0_xxyy[i] * gfe_0 + 2.0 * ts_y_xxy[i] * gfe2_0 + 2.0 * gr_y_xxy[i] * gfe_0 + ts_y_xxyy[i] * gfe_0 * gc_y[i] + gr_y_xxyy[i] * gc_y[i];

        grr_y_y_xxyz[i] = ts_0_xxyz[i] * gfe2_0 + gr_0_xxyz[i] * gfe_0 + ts_y_xxz[i] * gfe2_0 + gr_y_xxz[i] * gfe_0 + ts_y_xxyz[i] * gfe_0 * gc_y[i] + gr_y_xxyz[i] * gc_y[i];

        grr_y_y_xxzz[i] = ts_0_xxzz[i] * gfe2_0 + gr_0_xxzz[i] * gfe_0 + ts_y_xxzz[i] * gfe_0 * gc_y[i] + gr_y_xxzz[i] * gc_y[i];

        grr_y_y_xyyy[i] = ts_0_xyyy[i] * gfe2_0 + gr_0_xyyy[i] * gfe_0 + 3.0 * ts_y_xyy[i] * gfe2_0 + 3.0 * gr_y_xyy[i] * gfe_0 + ts_y_xyyy[i] * gfe_0 * gc_y[i] + gr_y_xyyy[i] * gc_y[i];

        grr_y_y_xyyz[i] = ts_0_xyyz[i] * gfe2_0 + gr_0_xyyz[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe2_0 + 2.0 * gr_y_xyz[i] * gfe_0 + ts_y_xyyz[i] * gfe_0 * gc_y[i] + gr_y_xyyz[i] * gc_y[i];

        grr_y_y_xyzz[i] = ts_0_xyzz[i] * gfe2_0 + gr_0_xyzz[i] * gfe_0 + ts_y_xzz[i] * gfe2_0 + gr_y_xzz[i] * gfe_0 + ts_y_xyzz[i] * gfe_0 * gc_y[i] + gr_y_xyzz[i] * gc_y[i];

        grr_y_y_xzzz[i] = ts_0_xzzz[i] * gfe2_0 + gr_0_xzzz[i] * gfe_0 + ts_y_xzzz[i] * gfe_0 * gc_y[i] + gr_y_xzzz[i] * gc_y[i];

        grr_y_y_yyyy[i] = ts_0_yyyy[i] * gfe2_0 + gr_0_yyyy[i] * gfe_0 + 4.0 * ts_y_yyy[i] * gfe2_0 + 4.0 * gr_y_yyy[i] * gfe_0 + ts_y_yyyy[i] * gfe_0 * gc_y[i] + gr_y_yyyy[i] * gc_y[i];

        grr_y_y_yyyz[i] = ts_0_yyyz[i] * gfe2_0 + gr_0_yyyz[i] * gfe_0 + 3.0 * ts_y_yyz[i] * gfe2_0 + 3.0 * gr_y_yyz[i] * gfe_0 + ts_y_yyyz[i] * gfe_0 * gc_y[i] + gr_y_yyyz[i] * gc_y[i];

        grr_y_y_yyzz[i] = ts_0_yyzz[i] * gfe2_0 + gr_0_yyzz[i] * gfe_0 + 2.0 * ts_y_yzz[i] * gfe2_0 + 2.0 * gr_y_yzz[i] * gfe_0 + ts_y_yyzz[i] * gfe_0 * gc_y[i] + gr_y_yyzz[i] * gc_y[i];

        grr_y_y_yzzz[i] = ts_0_yzzz[i] * gfe2_0 + gr_0_yzzz[i] * gfe_0 + ts_y_zzz[i] * gfe2_0 + gr_y_zzz[i] * gfe_0 + ts_y_yzzz[i] * gfe_0 * gc_y[i] + gr_y_yzzz[i] * gc_y[i];

        grr_y_y_zzzz[i] = ts_0_zzzz[i] * gfe2_0 + gr_0_zzzz[i] * gfe_0 + ts_y_zzzz[i] * gfe_0 * gc_y[i] + gr_y_zzzz[i] * gc_y[i];
    }

    // Set up 75-90 components of targeted buffer : PG

    auto grr_y_z_xxxx = pbuffer.data(idx_gr_pg + 75);

    auto grr_y_z_xxxy = pbuffer.data(idx_gr_pg + 76);

    auto grr_y_z_xxxz = pbuffer.data(idx_gr_pg + 77);

    auto grr_y_z_xxyy = pbuffer.data(idx_gr_pg + 78);

    auto grr_y_z_xxyz = pbuffer.data(idx_gr_pg + 79);

    auto grr_y_z_xxzz = pbuffer.data(idx_gr_pg + 80);

    auto grr_y_z_xyyy = pbuffer.data(idx_gr_pg + 81);

    auto grr_y_z_xyyz = pbuffer.data(idx_gr_pg + 82);

    auto grr_y_z_xyzz = pbuffer.data(idx_gr_pg + 83);

    auto grr_y_z_xzzz = pbuffer.data(idx_gr_pg + 84);

    auto grr_y_z_yyyy = pbuffer.data(idx_gr_pg + 85);

    auto grr_y_z_yyyz = pbuffer.data(idx_gr_pg + 86);

    auto grr_y_z_yyzz = pbuffer.data(idx_gr_pg + 87);

    auto grr_y_z_yzzz = pbuffer.data(idx_gr_pg + 88);

    auto grr_y_z_zzzz = pbuffer.data(idx_gr_pg + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xxx, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxy, gr_z_xxyy, gr_z_xxyz, gr_z_xxz, gr_z_xxzz, gr_z_xyy, gr_z_xyyy, gr_z_xyyz, gr_z_xyz, gr_z_xyzz, gr_z_xzz, gr_z_xzzz, gr_z_yyy, gr_z_yyyy, gr_z_yyyz, gr_z_yyz, gr_z_yyzz, gr_z_yzz, gr_z_yzzz, gr_z_zzz, gr_z_zzzz, grr_y_z_xxxx, grr_y_z_xxxy, grr_y_z_xxxz, grr_y_z_xxyy, grr_y_z_xxyz, grr_y_z_xxzz, grr_y_z_xyyy, grr_y_z_xyyz, grr_y_z_xyzz, grr_y_z_xzzz, grr_y_z_yyyy, grr_y_z_yyyz, grr_y_z_yyzz, grr_y_z_yzzz, grr_y_z_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_z_xxxx[i] = ts_z_xxxx[i] * gfe_0 * gc_y[i] + gr_z_xxxx[i] * gc_y[i];

        grr_y_z_xxxy[i] = ts_z_xxx[i] * gfe2_0 + gr_z_xxx[i] * gfe_0 + ts_z_xxxy[i] * gfe_0 * gc_y[i] + gr_z_xxxy[i] * gc_y[i];

        grr_y_z_xxxz[i] = ts_z_xxxz[i] * gfe_0 * gc_y[i] + gr_z_xxxz[i] * gc_y[i];

        grr_y_z_xxyy[i] = 2.0 * ts_z_xxy[i] * gfe2_0 + 2.0 * gr_z_xxy[i] * gfe_0 + ts_z_xxyy[i] * gfe_0 * gc_y[i] + gr_z_xxyy[i] * gc_y[i];

        grr_y_z_xxyz[i] = ts_z_xxz[i] * gfe2_0 + gr_z_xxz[i] * gfe_0 + ts_z_xxyz[i] * gfe_0 * gc_y[i] + gr_z_xxyz[i] * gc_y[i];

        grr_y_z_xxzz[i] = ts_z_xxzz[i] * gfe_0 * gc_y[i] + gr_z_xxzz[i] * gc_y[i];

        grr_y_z_xyyy[i] = 3.0 * ts_z_xyy[i] * gfe2_0 + 3.0 * gr_z_xyy[i] * gfe_0 + ts_z_xyyy[i] * gfe_0 * gc_y[i] + gr_z_xyyy[i] * gc_y[i];

        grr_y_z_xyyz[i] = 2.0 * ts_z_xyz[i] * gfe2_0 + 2.0 * gr_z_xyz[i] * gfe_0 + ts_z_xyyz[i] * gfe_0 * gc_y[i] + gr_z_xyyz[i] * gc_y[i];

        grr_y_z_xyzz[i] = ts_z_xzz[i] * gfe2_0 + gr_z_xzz[i] * gfe_0 + ts_z_xyzz[i] * gfe_0 * gc_y[i] + gr_z_xyzz[i] * gc_y[i];

        grr_y_z_xzzz[i] = ts_z_xzzz[i] * gfe_0 * gc_y[i] + gr_z_xzzz[i] * gc_y[i];

        grr_y_z_yyyy[i] = 4.0 * ts_z_yyy[i] * gfe2_0 + 4.0 * gr_z_yyy[i] * gfe_0 + ts_z_yyyy[i] * gfe_0 * gc_y[i] + gr_z_yyyy[i] * gc_y[i];

        grr_y_z_yyyz[i] = 3.0 * ts_z_yyz[i] * gfe2_0 + 3.0 * gr_z_yyz[i] * gfe_0 + ts_z_yyyz[i] * gfe_0 * gc_y[i] + gr_z_yyyz[i] * gc_y[i];

        grr_y_z_yyzz[i] = 2.0 * ts_z_yzz[i] * gfe2_0 + 2.0 * gr_z_yzz[i] * gfe_0 + ts_z_yyzz[i] * gfe_0 * gc_y[i] + gr_z_yyzz[i] * gc_y[i];

        grr_y_z_yzzz[i] = ts_z_zzz[i] * gfe2_0 + gr_z_zzz[i] * gfe_0 + ts_z_yzzz[i] * gfe_0 * gc_y[i] + gr_z_yzzz[i] * gc_y[i];

        grr_y_z_zzzz[i] = ts_z_zzzz[i] * gfe_0 * gc_y[i] + gr_z_zzzz[i] * gc_y[i];
    }

    // Set up 90-105 components of targeted buffer : PG

    auto grr_z_x_xxxx = pbuffer.data(idx_gr_pg + 90);

    auto grr_z_x_xxxy = pbuffer.data(idx_gr_pg + 91);

    auto grr_z_x_xxxz = pbuffer.data(idx_gr_pg + 92);

    auto grr_z_x_xxyy = pbuffer.data(idx_gr_pg + 93);

    auto grr_z_x_xxyz = pbuffer.data(idx_gr_pg + 94);

    auto grr_z_x_xxzz = pbuffer.data(idx_gr_pg + 95);

    auto grr_z_x_xyyy = pbuffer.data(idx_gr_pg + 96);

    auto grr_z_x_xyyz = pbuffer.data(idx_gr_pg + 97);

    auto grr_z_x_xyzz = pbuffer.data(idx_gr_pg + 98);

    auto grr_z_x_xzzz = pbuffer.data(idx_gr_pg + 99);

    auto grr_z_x_yyyy = pbuffer.data(idx_gr_pg + 100);

    auto grr_z_x_yyyz = pbuffer.data(idx_gr_pg + 101);

    auto grr_z_x_yyzz = pbuffer.data(idx_gr_pg + 102);

    auto grr_z_x_yzzz = pbuffer.data(idx_gr_pg + 103);

    auto grr_z_x_zzzz = pbuffer.data(idx_gr_pg + 104);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xxx, gr_x_xxxx, gr_x_xxxy, gr_x_xxxz, gr_x_xxy, gr_x_xxyy, gr_x_xxyz, gr_x_xxz, gr_x_xxzz, gr_x_xyy, gr_x_xyyy, gr_x_xyyz, gr_x_xyz, gr_x_xyzz, gr_x_xzz, gr_x_xzzz, gr_x_yyy, gr_x_yyyy, gr_x_yyyz, gr_x_yyz, gr_x_yyzz, gr_x_yzz, gr_x_yzzz, gr_x_zzz, gr_x_zzzz, grr_z_x_xxxx, grr_z_x_xxxy, grr_z_x_xxxz, grr_z_x_xxyy, grr_z_x_xxyz, grr_z_x_xxzz, grr_z_x_xyyy, grr_z_x_xyyz, grr_z_x_xyzz, grr_z_x_xzzz, grr_z_x_yyyy, grr_z_x_yyyz, grr_z_x_yyzz, grr_z_x_yzzz, grr_z_x_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_x_xxxx[i] = ts_x_xxxx[i] * gfe_0 * gc_z[i] + gr_x_xxxx[i] * gc_z[i];

        grr_z_x_xxxy[i] = ts_x_xxxy[i] * gfe_0 * gc_z[i] + gr_x_xxxy[i] * gc_z[i];

        grr_z_x_xxxz[i] = ts_x_xxx[i] * gfe2_0 + gr_x_xxx[i] * gfe_0 + ts_x_xxxz[i] * gfe_0 * gc_z[i] + gr_x_xxxz[i] * gc_z[i];

        grr_z_x_xxyy[i] = ts_x_xxyy[i] * gfe_0 * gc_z[i] + gr_x_xxyy[i] * gc_z[i];

        grr_z_x_xxyz[i] = ts_x_xxy[i] * gfe2_0 + gr_x_xxy[i] * gfe_0 + ts_x_xxyz[i] * gfe_0 * gc_z[i] + gr_x_xxyz[i] * gc_z[i];

        grr_z_x_xxzz[i] = 2.0 * ts_x_xxz[i] * gfe2_0 + 2.0 * gr_x_xxz[i] * gfe_0 + ts_x_xxzz[i] * gfe_0 * gc_z[i] + gr_x_xxzz[i] * gc_z[i];

        grr_z_x_xyyy[i] = ts_x_xyyy[i] * gfe_0 * gc_z[i] + gr_x_xyyy[i] * gc_z[i];

        grr_z_x_xyyz[i] = ts_x_xyy[i] * gfe2_0 + gr_x_xyy[i] * gfe_0 + ts_x_xyyz[i] * gfe_0 * gc_z[i] + gr_x_xyyz[i] * gc_z[i];

        grr_z_x_xyzz[i] = 2.0 * ts_x_xyz[i] * gfe2_0 + 2.0 * gr_x_xyz[i] * gfe_0 + ts_x_xyzz[i] * gfe_0 * gc_z[i] + gr_x_xyzz[i] * gc_z[i];

        grr_z_x_xzzz[i] = 3.0 * ts_x_xzz[i] * gfe2_0 + 3.0 * gr_x_xzz[i] * gfe_0 + ts_x_xzzz[i] * gfe_0 * gc_z[i] + gr_x_xzzz[i] * gc_z[i];

        grr_z_x_yyyy[i] = ts_x_yyyy[i] * gfe_0 * gc_z[i] + gr_x_yyyy[i] * gc_z[i];

        grr_z_x_yyyz[i] = ts_x_yyy[i] * gfe2_0 + gr_x_yyy[i] * gfe_0 + ts_x_yyyz[i] * gfe_0 * gc_z[i] + gr_x_yyyz[i] * gc_z[i];

        grr_z_x_yyzz[i] = 2.0 * ts_x_yyz[i] * gfe2_0 + 2.0 * gr_x_yyz[i] * gfe_0 + ts_x_yyzz[i] * gfe_0 * gc_z[i] + gr_x_yyzz[i] * gc_z[i];

        grr_z_x_yzzz[i] = 3.0 * ts_x_yzz[i] * gfe2_0 + 3.0 * gr_x_yzz[i] * gfe_0 + ts_x_yzzz[i] * gfe_0 * gc_z[i] + gr_x_yzzz[i] * gc_z[i];

        grr_z_x_zzzz[i] = 4.0 * ts_x_zzz[i] * gfe2_0 + 4.0 * gr_x_zzz[i] * gfe_0 + ts_x_zzzz[i] * gfe_0 * gc_z[i] + gr_x_zzzz[i] * gc_z[i];
    }

    // Set up 105-120 components of targeted buffer : PG

    auto grr_z_y_xxxx = pbuffer.data(idx_gr_pg + 105);

    auto grr_z_y_xxxy = pbuffer.data(idx_gr_pg + 106);

    auto grr_z_y_xxxz = pbuffer.data(idx_gr_pg + 107);

    auto grr_z_y_xxyy = pbuffer.data(idx_gr_pg + 108);

    auto grr_z_y_xxyz = pbuffer.data(idx_gr_pg + 109);

    auto grr_z_y_xxzz = pbuffer.data(idx_gr_pg + 110);

    auto grr_z_y_xyyy = pbuffer.data(idx_gr_pg + 111);

    auto grr_z_y_xyyz = pbuffer.data(idx_gr_pg + 112);

    auto grr_z_y_xyzz = pbuffer.data(idx_gr_pg + 113);

    auto grr_z_y_xzzz = pbuffer.data(idx_gr_pg + 114);

    auto grr_z_y_yyyy = pbuffer.data(idx_gr_pg + 115);

    auto grr_z_y_yyyz = pbuffer.data(idx_gr_pg + 116);

    auto grr_z_y_yyzz = pbuffer.data(idx_gr_pg + 117);

    auto grr_z_y_yzzz = pbuffer.data(idx_gr_pg + 118);

    auto grr_z_y_zzzz = pbuffer.data(idx_gr_pg + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xxx, gr_y_xxxx, gr_y_xxxy, gr_y_xxxz, gr_y_xxy, gr_y_xxyy, gr_y_xxyz, gr_y_xxz, gr_y_xxzz, gr_y_xyy, gr_y_xyyy, gr_y_xyyz, gr_y_xyz, gr_y_xyzz, gr_y_xzz, gr_y_xzzz, gr_y_yyy, gr_y_yyyy, gr_y_yyyz, gr_y_yyz, gr_y_yyzz, gr_y_yzz, gr_y_yzzz, gr_y_zzz, gr_y_zzzz, grr_z_y_xxxx, grr_z_y_xxxy, grr_z_y_xxxz, grr_z_y_xxyy, grr_z_y_xxyz, grr_z_y_xxzz, grr_z_y_xyyy, grr_z_y_xyyz, grr_z_y_xyzz, grr_z_y_xzzz, grr_z_y_yyyy, grr_z_y_yyyz, grr_z_y_yyzz, grr_z_y_yzzz, grr_z_y_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_y_xxxx[i] = ts_y_xxxx[i] * gfe_0 * gc_z[i] + gr_y_xxxx[i] * gc_z[i];

        grr_z_y_xxxy[i] = ts_y_xxxy[i] * gfe_0 * gc_z[i] + gr_y_xxxy[i] * gc_z[i];

        grr_z_y_xxxz[i] = ts_y_xxx[i] * gfe2_0 + gr_y_xxx[i] * gfe_0 + ts_y_xxxz[i] * gfe_0 * gc_z[i] + gr_y_xxxz[i] * gc_z[i];

        grr_z_y_xxyy[i] = ts_y_xxyy[i] * gfe_0 * gc_z[i] + gr_y_xxyy[i] * gc_z[i];

        grr_z_y_xxyz[i] = ts_y_xxy[i] * gfe2_0 + gr_y_xxy[i] * gfe_0 + ts_y_xxyz[i] * gfe_0 * gc_z[i] + gr_y_xxyz[i] * gc_z[i];

        grr_z_y_xxzz[i] = 2.0 * ts_y_xxz[i] * gfe2_0 + 2.0 * gr_y_xxz[i] * gfe_0 + ts_y_xxzz[i] * gfe_0 * gc_z[i] + gr_y_xxzz[i] * gc_z[i];

        grr_z_y_xyyy[i] = ts_y_xyyy[i] * gfe_0 * gc_z[i] + gr_y_xyyy[i] * gc_z[i];

        grr_z_y_xyyz[i] = ts_y_xyy[i] * gfe2_0 + gr_y_xyy[i] * gfe_0 + ts_y_xyyz[i] * gfe_0 * gc_z[i] + gr_y_xyyz[i] * gc_z[i];

        grr_z_y_xyzz[i] = 2.0 * ts_y_xyz[i] * gfe2_0 + 2.0 * gr_y_xyz[i] * gfe_0 + ts_y_xyzz[i] * gfe_0 * gc_z[i] + gr_y_xyzz[i] * gc_z[i];

        grr_z_y_xzzz[i] = 3.0 * ts_y_xzz[i] * gfe2_0 + 3.0 * gr_y_xzz[i] * gfe_0 + ts_y_xzzz[i] * gfe_0 * gc_z[i] + gr_y_xzzz[i] * gc_z[i];

        grr_z_y_yyyy[i] = ts_y_yyyy[i] * gfe_0 * gc_z[i] + gr_y_yyyy[i] * gc_z[i];

        grr_z_y_yyyz[i] = ts_y_yyy[i] * gfe2_0 + gr_y_yyy[i] * gfe_0 + ts_y_yyyz[i] * gfe_0 * gc_z[i] + gr_y_yyyz[i] * gc_z[i];

        grr_z_y_yyzz[i] = 2.0 * ts_y_yyz[i] * gfe2_0 + 2.0 * gr_y_yyz[i] * gfe_0 + ts_y_yyzz[i] * gfe_0 * gc_z[i] + gr_y_yyzz[i] * gc_z[i];

        grr_z_y_yzzz[i] = 3.0 * ts_y_yzz[i] * gfe2_0 + 3.0 * gr_y_yzz[i] * gfe_0 + ts_y_yzzz[i] * gfe_0 * gc_z[i] + gr_y_yzzz[i] * gc_z[i];

        grr_z_y_zzzz[i] = 4.0 * ts_y_zzz[i] * gfe2_0 + 4.0 * gr_y_zzz[i] * gfe_0 + ts_y_zzzz[i] * gfe_0 * gc_z[i] + gr_y_zzzz[i] * gc_z[i];
    }

    // Set up 120-135 components of targeted buffer : PG

    auto grr_z_z_xxxx = pbuffer.data(idx_gr_pg + 120);

    auto grr_z_z_xxxy = pbuffer.data(idx_gr_pg + 121);

    auto grr_z_z_xxxz = pbuffer.data(idx_gr_pg + 122);

    auto grr_z_z_xxyy = pbuffer.data(idx_gr_pg + 123);

    auto grr_z_z_xxyz = pbuffer.data(idx_gr_pg + 124);

    auto grr_z_z_xxzz = pbuffer.data(idx_gr_pg + 125);

    auto grr_z_z_xyyy = pbuffer.data(idx_gr_pg + 126);

    auto grr_z_z_xyyz = pbuffer.data(idx_gr_pg + 127);

    auto grr_z_z_xyzz = pbuffer.data(idx_gr_pg + 128);

    auto grr_z_z_xzzz = pbuffer.data(idx_gr_pg + 129);

    auto grr_z_z_yyyy = pbuffer.data(idx_gr_pg + 130);

    auto grr_z_z_yyyz = pbuffer.data(idx_gr_pg + 131);

    auto grr_z_z_yyzz = pbuffer.data(idx_gr_pg + 132);

    auto grr_z_z_yzzz = pbuffer.data(idx_gr_pg + 133);

    auto grr_z_z_zzzz = pbuffer.data(idx_gr_pg + 134);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xxxx, gr_0_xxxy, gr_0_xxxz, gr_0_xxyy, gr_0_xxyz, gr_0_xxzz, gr_0_xyyy, gr_0_xyyz, gr_0_xyzz, gr_0_xzzz, gr_0_yyyy, gr_0_yyyz, gr_0_yyzz, gr_0_yzzz, gr_0_zzzz, gr_z_xxx, gr_z_xxxx, gr_z_xxxy, gr_z_xxxz, gr_z_xxy, gr_z_xxyy, gr_z_xxyz, gr_z_xxz, gr_z_xxzz, gr_z_xyy, gr_z_xyyy, gr_z_xyyz, gr_z_xyz, gr_z_xyzz, gr_z_xzz, gr_z_xzzz, gr_z_yyy, gr_z_yyyy, gr_z_yyyz, gr_z_yyz, gr_z_yyzz, gr_z_yzz, gr_z_yzzz, gr_z_zzz, gr_z_zzzz, grr_z_z_xxxx, grr_z_z_xxxy, grr_z_z_xxxz, grr_z_z_xxyy, grr_z_z_xxyz, grr_z_z_xxzz, grr_z_z_xyyy, grr_z_z_xyyz, grr_z_z_xyzz, grr_z_z_xzzz, grr_z_z_yyyy, grr_z_z_yyyz, grr_z_z_yyzz, grr_z_z_yzzz, grr_z_z_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_z_xxxx[i] = ts_0_xxxx[i] * gfe2_0 + gr_0_xxxx[i] * gfe_0 + ts_z_xxxx[i] * gfe_0 * gc_z[i] + gr_z_xxxx[i] * gc_z[i];

        grr_z_z_xxxy[i] = ts_0_xxxy[i] * gfe2_0 + gr_0_xxxy[i] * gfe_0 + ts_z_xxxy[i] * gfe_0 * gc_z[i] + gr_z_xxxy[i] * gc_z[i];

        grr_z_z_xxxz[i] = ts_0_xxxz[i] * gfe2_0 + gr_0_xxxz[i] * gfe_0 + ts_z_xxx[i] * gfe2_0 + gr_z_xxx[i] * gfe_0 + ts_z_xxxz[i] * gfe_0 * gc_z[i] + gr_z_xxxz[i] * gc_z[i];

        grr_z_z_xxyy[i] = ts_0_xxyy[i] * gfe2_0 + gr_0_xxyy[i] * gfe_0 + ts_z_xxyy[i] * gfe_0 * gc_z[i] + gr_z_xxyy[i] * gc_z[i];

        grr_z_z_xxyz[i] = ts_0_xxyz[i] * gfe2_0 + gr_0_xxyz[i] * gfe_0 + ts_z_xxy[i] * gfe2_0 + gr_z_xxy[i] * gfe_0 + ts_z_xxyz[i] * gfe_0 * gc_z[i] + gr_z_xxyz[i] * gc_z[i];

        grr_z_z_xxzz[i] = ts_0_xxzz[i] * gfe2_0 + gr_0_xxzz[i] * gfe_0 + 2.0 * ts_z_xxz[i] * gfe2_0 + 2.0 * gr_z_xxz[i] * gfe_0 + ts_z_xxzz[i] * gfe_0 * gc_z[i] + gr_z_xxzz[i] * gc_z[i];

        grr_z_z_xyyy[i] = ts_0_xyyy[i] * gfe2_0 + gr_0_xyyy[i] * gfe_0 + ts_z_xyyy[i] * gfe_0 * gc_z[i] + gr_z_xyyy[i] * gc_z[i];

        grr_z_z_xyyz[i] = ts_0_xyyz[i] * gfe2_0 + gr_0_xyyz[i] * gfe_0 + ts_z_xyy[i] * gfe2_0 + gr_z_xyy[i] * gfe_0 + ts_z_xyyz[i] * gfe_0 * gc_z[i] + gr_z_xyyz[i] * gc_z[i];

        grr_z_z_xyzz[i] = ts_0_xyzz[i] * gfe2_0 + gr_0_xyzz[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe2_0 + 2.0 * gr_z_xyz[i] * gfe_0 + ts_z_xyzz[i] * gfe_0 * gc_z[i] + gr_z_xyzz[i] * gc_z[i];

        grr_z_z_xzzz[i] = ts_0_xzzz[i] * gfe2_0 + gr_0_xzzz[i] * gfe_0 + 3.0 * ts_z_xzz[i] * gfe2_0 + 3.0 * gr_z_xzz[i] * gfe_0 + ts_z_xzzz[i] * gfe_0 * gc_z[i] + gr_z_xzzz[i] * gc_z[i];

        grr_z_z_yyyy[i] = ts_0_yyyy[i] * gfe2_0 + gr_0_yyyy[i] * gfe_0 + ts_z_yyyy[i] * gfe_0 * gc_z[i] + gr_z_yyyy[i] * gc_z[i];

        grr_z_z_yyyz[i] = ts_0_yyyz[i] * gfe2_0 + gr_0_yyyz[i] * gfe_0 + ts_z_yyy[i] * gfe2_0 + gr_z_yyy[i] * gfe_0 + ts_z_yyyz[i] * gfe_0 * gc_z[i] + gr_z_yyyz[i] * gc_z[i];

        grr_z_z_yyzz[i] = ts_0_yyzz[i] * gfe2_0 + gr_0_yyzz[i] * gfe_0 + 2.0 * ts_z_yyz[i] * gfe2_0 + 2.0 * gr_z_yyz[i] * gfe_0 + ts_z_yyzz[i] * gfe_0 * gc_z[i] + gr_z_yyzz[i] * gc_z[i];

        grr_z_z_yzzz[i] = ts_0_yzzz[i] * gfe2_0 + gr_0_yzzz[i] * gfe_0 + 3.0 * ts_z_yzz[i] * gfe2_0 + 3.0 * gr_z_yzz[i] * gfe_0 + ts_z_yzzz[i] * gfe_0 * gc_z[i] + gr_z_yzzz[i] * gc_z[i];

        grr_z_z_zzzz[i] = ts_0_zzzz[i] * gfe2_0 + gr_0_zzzz[i] * gfe_0 + 4.0 * ts_z_zzz[i] * gfe2_0 + 4.0 * gr_z_zzz[i] * gfe_0 + ts_z_zzzz[i] * gfe_0 * gc_z[i] + gr_z_zzzz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

