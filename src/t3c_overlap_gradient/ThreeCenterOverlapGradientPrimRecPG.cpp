#include "ThreeCenterOverlapGradientPrimRecPG.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_pg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_pg,
                              const size_t idx_sg,
                              const size_t idx_pf,
                              const size_t idx_pg,
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

    // Set up 0-15 components of targeted buffer : PG

    auto gs_x_x_xxxx = pbuffer.data(idx_g_pg);

    auto gs_x_x_xxxy = pbuffer.data(idx_g_pg + 1);

    auto gs_x_x_xxxz = pbuffer.data(idx_g_pg + 2);

    auto gs_x_x_xxyy = pbuffer.data(idx_g_pg + 3);

    auto gs_x_x_xxyz = pbuffer.data(idx_g_pg + 4);

    auto gs_x_x_xxzz = pbuffer.data(idx_g_pg + 5);

    auto gs_x_x_xyyy = pbuffer.data(idx_g_pg + 6);

    auto gs_x_x_xyyz = pbuffer.data(idx_g_pg + 7);

    auto gs_x_x_xyzz = pbuffer.data(idx_g_pg + 8);

    auto gs_x_x_xzzz = pbuffer.data(idx_g_pg + 9);

    auto gs_x_x_yyyy = pbuffer.data(idx_g_pg + 10);

    auto gs_x_x_yyyz = pbuffer.data(idx_g_pg + 11);

    auto gs_x_x_yyzz = pbuffer.data(idx_g_pg + 12);

    auto gs_x_x_yzzz = pbuffer.data(idx_g_pg + 13);

    auto gs_x_x_zzzz = pbuffer.data(idx_g_pg + 14);

    #pragma omp simd aligned(gc_x, gs_x_x_xxxx, gs_x_x_xxxy, gs_x_x_xxxz, gs_x_x_xxyy, gs_x_x_xxyz, gs_x_x_xxzz, gs_x_x_xyyy, gs_x_x_xyyz, gs_x_x_xyzz, gs_x_x_xzzz, gs_x_x_yyyy, gs_x_x_yyyz, gs_x_x_yyzz, gs_x_x_yzzz, gs_x_x_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_x_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_x_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxx[i] * gc_x[i] * tce_0;

        gs_x_x_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_x_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxy[i] * gc_x[i] * tce_0;

        gs_x_x_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_x_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxz[i] * gc_x[i] * tce_0;

        gs_x_x_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_x_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyy[i] * gc_x[i] * tce_0;

        gs_x_x_xxyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_x_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyz[i] * gc_x[i] * tce_0;

        gs_x_x_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_x_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxzz[i] * gc_x[i] * tce_0;

        gs_x_x_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyy[i] * gc_x[i] * tce_0;

        gs_x_x_xyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyz[i] * gc_x[i] * tce_0;

        gs_x_x_xyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyzz[i] * gc_x[i] * tce_0;

        gs_x_x_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xzzz[i] * gc_x[i] * tce_0;

        gs_x_x_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyy[i] * gc_x[i] * tce_0;

        gs_x_x_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyz[i] * gc_x[i] * tce_0;

        gs_x_x_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzz[i] * gc_x[i] * tce_0;

        gs_x_x_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzz[i] * gc_x[i] * tce_0;

        gs_x_x_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 15-30 components of targeted buffer : PG

    auto gs_x_y_xxxx = pbuffer.data(idx_g_pg + 15);

    auto gs_x_y_xxxy = pbuffer.data(idx_g_pg + 16);

    auto gs_x_y_xxxz = pbuffer.data(idx_g_pg + 17);

    auto gs_x_y_xxyy = pbuffer.data(idx_g_pg + 18);

    auto gs_x_y_xxyz = pbuffer.data(idx_g_pg + 19);

    auto gs_x_y_xxzz = pbuffer.data(idx_g_pg + 20);

    auto gs_x_y_xyyy = pbuffer.data(idx_g_pg + 21);

    auto gs_x_y_xyyz = pbuffer.data(idx_g_pg + 22);

    auto gs_x_y_xyzz = pbuffer.data(idx_g_pg + 23);

    auto gs_x_y_xzzz = pbuffer.data(idx_g_pg + 24);

    auto gs_x_y_yyyy = pbuffer.data(idx_g_pg + 25);

    auto gs_x_y_yyyz = pbuffer.data(idx_g_pg + 26);

    auto gs_x_y_yyzz = pbuffer.data(idx_g_pg + 27);

    auto gs_x_y_yzzz = pbuffer.data(idx_g_pg + 28);

    auto gs_x_y_zzzz = pbuffer.data(idx_g_pg + 29);

    #pragma omp simd aligned(gc_x, gs_x_y_xxxx, gs_x_y_xxxy, gs_x_y_xxxz, gs_x_y_xxyy, gs_x_y_xxyz, gs_x_y_xxzz, gs_x_y_xyyy, gs_x_y_xyyz, gs_x_y_xyzz, gs_x_y_xzzz, gs_x_y_yyyy, gs_x_y_yyyz, gs_x_y_yyzz, gs_x_y_yzzz, gs_x_y_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_y_xxxx[i] = 8.0 * ts_y_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxx[i] * gc_x[i] * tce_0;

        gs_x_y_xxxy[i] = 6.0 * ts_y_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxy[i] * gc_x[i] * tce_0;

        gs_x_y_xxxz[i] = 6.0 * ts_y_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxz[i] * gc_x[i] * tce_0;

        gs_x_y_xxyy[i] = 4.0 * ts_y_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyy[i] * gc_x[i] * tce_0;

        gs_x_y_xxyz[i] = 4.0 * ts_y_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyz[i] * gc_x[i] * tce_0;

        gs_x_y_xxzz[i] = 4.0 * ts_y_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzz[i] * gc_x[i] * tce_0;

        gs_x_y_xyyy[i] = 2.0 * ts_y_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyy[i] * gc_x[i] * tce_0;

        gs_x_y_xyyz[i] = 2.0 * ts_y_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyz[i] * gc_x[i] * tce_0;

        gs_x_y_xyzz[i] = 2.0 * ts_y_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyzz[i] * gc_x[i] * tce_0;

        gs_x_y_xzzz[i] = 2.0 * ts_y_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzz[i] * gc_x[i] * tce_0;

        gs_x_y_yyyy[i] = 2.0 * ts_y_yyyy[i] * gc_x[i] * tce_0;

        gs_x_y_yyyz[i] = 2.0 * ts_y_yyyz[i] * gc_x[i] * tce_0;

        gs_x_y_yyzz[i] = 2.0 * ts_y_yyzz[i] * gc_x[i] * tce_0;

        gs_x_y_yzzz[i] = 2.0 * ts_y_yzzz[i] * gc_x[i] * tce_0;

        gs_x_y_zzzz[i] = 2.0 * ts_y_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-45 components of targeted buffer : PG

    auto gs_x_z_xxxx = pbuffer.data(idx_g_pg + 30);

    auto gs_x_z_xxxy = pbuffer.data(idx_g_pg + 31);

    auto gs_x_z_xxxz = pbuffer.data(idx_g_pg + 32);

    auto gs_x_z_xxyy = pbuffer.data(idx_g_pg + 33);

    auto gs_x_z_xxyz = pbuffer.data(idx_g_pg + 34);

    auto gs_x_z_xxzz = pbuffer.data(idx_g_pg + 35);

    auto gs_x_z_xyyy = pbuffer.data(idx_g_pg + 36);

    auto gs_x_z_xyyz = pbuffer.data(idx_g_pg + 37);

    auto gs_x_z_xyzz = pbuffer.data(idx_g_pg + 38);

    auto gs_x_z_xzzz = pbuffer.data(idx_g_pg + 39);

    auto gs_x_z_yyyy = pbuffer.data(idx_g_pg + 40);

    auto gs_x_z_yyyz = pbuffer.data(idx_g_pg + 41);

    auto gs_x_z_yyzz = pbuffer.data(idx_g_pg + 42);

    auto gs_x_z_yzzz = pbuffer.data(idx_g_pg + 43);

    auto gs_x_z_zzzz = pbuffer.data(idx_g_pg + 44);

    #pragma omp simd aligned(gc_x, gs_x_z_xxxx, gs_x_z_xxxy, gs_x_z_xxxz, gs_x_z_xxyy, gs_x_z_xxyz, gs_x_z_xxzz, gs_x_z_xyyy, gs_x_z_xyyz, gs_x_z_xyzz, gs_x_z_xzzz, gs_x_z_yyyy, gs_x_z_yyyz, gs_x_z_yyzz, gs_x_z_yzzz, gs_x_z_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_z_xxxx[i] = 8.0 * ts_z_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxx[i] * gc_x[i] * tce_0;

        gs_x_z_xxxy[i] = 6.0 * ts_z_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxy[i] * gc_x[i] * tce_0;

        gs_x_z_xxxz[i] = 6.0 * ts_z_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxz[i] * gc_x[i] * tce_0;

        gs_x_z_xxyy[i] = 4.0 * ts_z_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyy[i] * gc_x[i] * tce_0;

        gs_x_z_xxyz[i] = 4.0 * ts_z_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyz[i] * gc_x[i] * tce_0;

        gs_x_z_xxzz[i] = 4.0 * ts_z_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxzz[i] * gc_x[i] * tce_0;

        gs_x_z_xyyy[i] = 2.0 * ts_z_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyy[i] * gc_x[i] * tce_0;

        gs_x_z_xyyz[i] = 2.0 * ts_z_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyz[i] * gc_x[i] * tce_0;

        gs_x_z_xyzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyzz[i] * gc_x[i] * tce_0;

        gs_x_z_xzzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xzzz[i] * gc_x[i] * tce_0;

        gs_x_z_yyyy[i] = 2.0 * ts_z_yyyy[i] * gc_x[i] * tce_0;

        gs_x_z_yyyz[i] = 2.0 * ts_z_yyyz[i] * gc_x[i] * tce_0;

        gs_x_z_yyzz[i] = 2.0 * ts_z_yyzz[i] * gc_x[i] * tce_0;

        gs_x_z_yzzz[i] = 2.0 * ts_z_yzzz[i] * gc_x[i] * tce_0;

        gs_x_z_zzzz[i] = 2.0 * ts_z_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 45-60 components of targeted buffer : PG

    auto gs_y_x_xxxx = pbuffer.data(idx_g_pg + 45);

    auto gs_y_x_xxxy = pbuffer.data(idx_g_pg + 46);

    auto gs_y_x_xxxz = pbuffer.data(idx_g_pg + 47);

    auto gs_y_x_xxyy = pbuffer.data(idx_g_pg + 48);

    auto gs_y_x_xxyz = pbuffer.data(idx_g_pg + 49);

    auto gs_y_x_xxzz = pbuffer.data(idx_g_pg + 50);

    auto gs_y_x_xyyy = pbuffer.data(idx_g_pg + 51);

    auto gs_y_x_xyyz = pbuffer.data(idx_g_pg + 52);

    auto gs_y_x_xyzz = pbuffer.data(idx_g_pg + 53);

    auto gs_y_x_xzzz = pbuffer.data(idx_g_pg + 54);

    auto gs_y_x_yyyy = pbuffer.data(idx_g_pg + 55);

    auto gs_y_x_yyyz = pbuffer.data(idx_g_pg + 56);

    auto gs_y_x_yyzz = pbuffer.data(idx_g_pg + 57);

    auto gs_y_x_yzzz = pbuffer.data(idx_g_pg + 58);

    auto gs_y_x_zzzz = pbuffer.data(idx_g_pg + 59);

    #pragma omp simd aligned(gc_y, gs_y_x_xxxx, gs_y_x_xxxy, gs_y_x_xxxz, gs_y_x_xxyy, gs_y_x_xxyz, gs_y_x_xxzz, gs_y_x_xyyy, gs_y_x_xyyz, gs_y_x_xyzz, gs_y_x_xzzz, gs_y_x_yyyy, gs_y_x_yyyz, gs_y_x_yyzz, gs_y_x_yzzz, gs_y_x_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_x_xxxx[i] = 2.0 * ts_x_xxxx[i] * gc_y[i] * tce_0;

        gs_y_x_xxxy[i] = 2.0 * ts_x_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxy[i] * gc_y[i] * tce_0;

        gs_y_x_xxxz[i] = 2.0 * ts_x_xxxz[i] * gc_y[i] * tce_0;

        gs_y_x_xxyy[i] = 4.0 * ts_x_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyy[i] * gc_y[i] * tce_0;

        gs_y_x_xxyz[i] = 2.0 * ts_x_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyz[i] * gc_y[i] * tce_0;

        gs_y_x_xxzz[i] = 2.0 * ts_x_xxzz[i] * gc_y[i] * tce_0;

        gs_y_x_xyyy[i] = 6.0 * ts_x_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyy[i] * gc_y[i] * tce_0;

        gs_y_x_xyyz[i] = 4.0 * ts_x_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyz[i] * gc_y[i] * tce_0;

        gs_y_x_xyzz[i] = 2.0 * ts_x_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyzz[i] * gc_y[i] * tce_0;

        gs_y_x_xzzz[i] = 2.0 * ts_x_xzzz[i] * gc_y[i] * tce_0;

        gs_y_x_yyyy[i] = 8.0 * ts_x_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyy[i] * gc_y[i] * tce_0;

        gs_y_x_yyyz[i] = 6.0 * ts_x_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyz[i] * gc_y[i] * tce_0;

        gs_y_x_yyzz[i] = 4.0 * ts_x_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzz[i] * gc_y[i] * tce_0;

        gs_y_x_yzzz[i] = 2.0 * ts_x_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzz[i] * gc_y[i] * tce_0;

        gs_y_x_zzzz[i] = 2.0 * ts_x_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 60-75 components of targeted buffer : PG

    auto gs_y_y_xxxx = pbuffer.data(idx_g_pg + 60);

    auto gs_y_y_xxxy = pbuffer.data(idx_g_pg + 61);

    auto gs_y_y_xxxz = pbuffer.data(idx_g_pg + 62);

    auto gs_y_y_xxyy = pbuffer.data(idx_g_pg + 63);

    auto gs_y_y_xxyz = pbuffer.data(idx_g_pg + 64);

    auto gs_y_y_xxzz = pbuffer.data(idx_g_pg + 65);

    auto gs_y_y_xyyy = pbuffer.data(idx_g_pg + 66);

    auto gs_y_y_xyyz = pbuffer.data(idx_g_pg + 67);

    auto gs_y_y_xyzz = pbuffer.data(idx_g_pg + 68);

    auto gs_y_y_xzzz = pbuffer.data(idx_g_pg + 69);

    auto gs_y_y_yyyy = pbuffer.data(idx_g_pg + 70);

    auto gs_y_y_yyyz = pbuffer.data(idx_g_pg + 71);

    auto gs_y_y_yyzz = pbuffer.data(idx_g_pg + 72);

    auto gs_y_y_yzzz = pbuffer.data(idx_g_pg + 73);

    auto gs_y_y_zzzz = pbuffer.data(idx_g_pg + 74);

    #pragma omp simd aligned(gc_y, gs_y_y_xxxx, gs_y_y_xxxy, gs_y_y_xxxz, gs_y_y_xxyy, gs_y_y_xxyz, gs_y_y_xxzz, gs_y_y_xyyy, gs_y_y_xyyz, gs_y_y_xyzz, gs_y_y_xzzz, gs_y_y_yyyy, gs_y_y_yyyz, gs_y_y_yyzz, gs_y_y_yzzz, gs_y_y_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_y_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxx[i] * gc_y[i] * tce_0;

        gs_y_y_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxy[i] * gc_y[i] * tce_0;

        gs_y_y_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxz[i] * gc_y[i] * tce_0;

        gs_y_y_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_y_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyy[i] * gc_y[i] * tce_0;

        gs_y_y_xxyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyz[i] * gc_y[i] * tce_0;

        gs_y_y_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzz[i] * gc_y[i] * tce_0;

        gs_y_y_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_y_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyy[i] * gc_y[i] * tce_0;

        gs_y_y_xyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_y_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyz[i] * gc_y[i] * tce_0;

        gs_y_y_xyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyzz[i] * gc_y[i] * tce_0;

        gs_y_y_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzz[i] * gc_y[i] * tce_0;

        gs_y_y_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_y_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyy[i] * gc_y[i] * tce_0;

        gs_y_y_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_y_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyz[i] * gc_y[i] * tce_0;

        gs_y_y_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_y_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyzz[i] * gc_y[i] * tce_0;

        gs_y_y_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yzzz[i] * gc_y[i] * tce_0;

        gs_y_y_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 75-90 components of targeted buffer : PG

    auto gs_y_z_xxxx = pbuffer.data(idx_g_pg + 75);

    auto gs_y_z_xxxy = pbuffer.data(idx_g_pg + 76);

    auto gs_y_z_xxxz = pbuffer.data(idx_g_pg + 77);

    auto gs_y_z_xxyy = pbuffer.data(idx_g_pg + 78);

    auto gs_y_z_xxyz = pbuffer.data(idx_g_pg + 79);

    auto gs_y_z_xxzz = pbuffer.data(idx_g_pg + 80);

    auto gs_y_z_xyyy = pbuffer.data(idx_g_pg + 81);

    auto gs_y_z_xyyz = pbuffer.data(idx_g_pg + 82);

    auto gs_y_z_xyzz = pbuffer.data(idx_g_pg + 83);

    auto gs_y_z_xzzz = pbuffer.data(idx_g_pg + 84);

    auto gs_y_z_yyyy = pbuffer.data(idx_g_pg + 85);

    auto gs_y_z_yyyz = pbuffer.data(idx_g_pg + 86);

    auto gs_y_z_yyzz = pbuffer.data(idx_g_pg + 87);

    auto gs_y_z_yzzz = pbuffer.data(idx_g_pg + 88);

    auto gs_y_z_zzzz = pbuffer.data(idx_g_pg + 89);

    #pragma omp simd aligned(gc_y, gs_y_z_xxxx, gs_y_z_xxxy, gs_y_z_xxxz, gs_y_z_xxyy, gs_y_z_xxyz, gs_y_z_xxzz, gs_y_z_xyyy, gs_y_z_xyyz, gs_y_z_xyzz, gs_y_z_xzzz, gs_y_z_yyyy, gs_y_z_yyyz, gs_y_z_yyzz, gs_y_z_yzzz, gs_y_z_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_z_xxxx[i] = 2.0 * ts_z_xxxx[i] * gc_y[i] * tce_0;

        gs_y_z_xxxy[i] = 2.0 * ts_z_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxy[i] * gc_y[i] * tce_0;

        gs_y_z_xxxz[i] = 2.0 * ts_z_xxxz[i] * gc_y[i] * tce_0;

        gs_y_z_xxyy[i] = 4.0 * ts_z_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyy[i] * gc_y[i] * tce_0;

        gs_y_z_xxyz[i] = 2.0 * ts_z_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyz[i] * gc_y[i] * tce_0;

        gs_y_z_xxzz[i] = 2.0 * ts_z_xxzz[i] * gc_y[i] * tce_0;

        gs_y_z_xyyy[i] = 6.0 * ts_z_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyy[i] * gc_y[i] * tce_0;

        gs_y_z_xyyz[i] = 4.0 * ts_z_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyz[i] * gc_y[i] * tce_0;

        gs_y_z_xyzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyzz[i] * gc_y[i] * tce_0;

        gs_y_z_xzzz[i] = 2.0 * ts_z_xzzz[i] * gc_y[i] * tce_0;

        gs_y_z_yyyy[i] = 8.0 * ts_z_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyy[i] * gc_y[i] * tce_0;

        gs_y_z_yyyz[i] = 6.0 * ts_z_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyz[i] * gc_y[i] * tce_0;

        gs_y_z_yyzz[i] = 4.0 * ts_z_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyzz[i] * gc_y[i] * tce_0;

        gs_y_z_yzzz[i] = 2.0 * ts_z_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yzzz[i] * gc_y[i] * tce_0;

        gs_y_z_zzzz[i] = 2.0 * ts_z_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 90-105 components of targeted buffer : PG

    auto gs_z_x_xxxx = pbuffer.data(idx_g_pg + 90);

    auto gs_z_x_xxxy = pbuffer.data(idx_g_pg + 91);

    auto gs_z_x_xxxz = pbuffer.data(idx_g_pg + 92);

    auto gs_z_x_xxyy = pbuffer.data(idx_g_pg + 93);

    auto gs_z_x_xxyz = pbuffer.data(idx_g_pg + 94);

    auto gs_z_x_xxzz = pbuffer.data(idx_g_pg + 95);

    auto gs_z_x_xyyy = pbuffer.data(idx_g_pg + 96);

    auto gs_z_x_xyyz = pbuffer.data(idx_g_pg + 97);

    auto gs_z_x_xyzz = pbuffer.data(idx_g_pg + 98);

    auto gs_z_x_xzzz = pbuffer.data(idx_g_pg + 99);

    auto gs_z_x_yyyy = pbuffer.data(idx_g_pg + 100);

    auto gs_z_x_yyyz = pbuffer.data(idx_g_pg + 101);

    auto gs_z_x_yyzz = pbuffer.data(idx_g_pg + 102);

    auto gs_z_x_yzzz = pbuffer.data(idx_g_pg + 103);

    auto gs_z_x_zzzz = pbuffer.data(idx_g_pg + 104);

    #pragma omp simd aligned(gc_z, gs_z_x_xxxx, gs_z_x_xxxy, gs_z_x_xxxz, gs_z_x_xxyy, gs_z_x_xxyz, gs_z_x_xxzz, gs_z_x_xyyy, gs_z_x_xyyz, gs_z_x_xyzz, gs_z_x_xzzz, gs_z_x_yyyy, gs_z_x_yyyz, gs_z_x_yyzz, gs_z_x_yzzz, gs_z_x_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_x_xxxx[i] = 2.0 * ts_x_xxxx[i] * gc_z[i] * tce_0;

        gs_z_x_xxxy[i] = 2.0 * ts_x_xxxy[i] * gc_z[i] * tce_0;

        gs_z_x_xxxz[i] = 2.0 * ts_x_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxxz[i] * gc_z[i] * tce_0;

        gs_z_x_xxyy[i] = 2.0 * ts_x_xxyy[i] * gc_z[i] * tce_0;

        gs_z_x_xxyz[i] = 2.0 * ts_x_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxyz[i] * gc_z[i] * tce_0;

        gs_z_x_xxzz[i] = 4.0 * ts_x_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xxzz[i] * gc_z[i] * tce_0;

        gs_z_x_xyyy[i] = 2.0 * ts_x_xyyy[i] * gc_z[i] * tce_0;

        gs_z_x_xyyz[i] = 2.0 * ts_x_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyyz[i] * gc_z[i] * tce_0;

        gs_z_x_xyzz[i] = 4.0 * ts_x_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xyzz[i] * gc_z[i] * tce_0;

        gs_z_x_xzzz[i] = 6.0 * ts_x_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_xzzz[i] * gc_z[i] * tce_0;

        gs_z_x_yyyy[i] = 2.0 * ts_x_yyyy[i] * gc_z[i] * tce_0;

        gs_z_x_yyyz[i] = 2.0 * ts_x_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyyz[i] * gc_z[i] * tce_0;

        gs_z_x_yyzz[i] = 4.0 * ts_x_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yyzz[i] * gc_z[i] * tce_0;

        gs_z_x_yzzz[i] = 6.0 * ts_x_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_yzzz[i] * gc_z[i] * tce_0;

        gs_z_x_zzzz[i] = 8.0 * ts_x_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_x_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 105-120 components of targeted buffer : PG

    auto gs_z_y_xxxx = pbuffer.data(idx_g_pg + 105);

    auto gs_z_y_xxxy = pbuffer.data(idx_g_pg + 106);

    auto gs_z_y_xxxz = pbuffer.data(idx_g_pg + 107);

    auto gs_z_y_xxyy = pbuffer.data(idx_g_pg + 108);

    auto gs_z_y_xxyz = pbuffer.data(idx_g_pg + 109);

    auto gs_z_y_xxzz = pbuffer.data(idx_g_pg + 110);

    auto gs_z_y_xyyy = pbuffer.data(idx_g_pg + 111);

    auto gs_z_y_xyyz = pbuffer.data(idx_g_pg + 112);

    auto gs_z_y_xyzz = pbuffer.data(idx_g_pg + 113);

    auto gs_z_y_xzzz = pbuffer.data(idx_g_pg + 114);

    auto gs_z_y_yyyy = pbuffer.data(idx_g_pg + 115);

    auto gs_z_y_yyyz = pbuffer.data(idx_g_pg + 116);

    auto gs_z_y_yyzz = pbuffer.data(idx_g_pg + 117);

    auto gs_z_y_yzzz = pbuffer.data(idx_g_pg + 118);

    auto gs_z_y_zzzz = pbuffer.data(idx_g_pg + 119);

    #pragma omp simd aligned(gc_z, gs_z_y_xxxx, gs_z_y_xxxy, gs_z_y_xxxz, gs_z_y_xxyy, gs_z_y_xxyz, gs_z_y_xxzz, gs_z_y_xyyy, gs_z_y_xyyz, gs_z_y_xyzz, gs_z_y_xzzz, gs_z_y_yyyy, gs_z_y_yyyz, gs_z_y_yyzz, gs_z_y_yzzz, gs_z_y_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_y_xxxx[i] = 2.0 * ts_y_xxxx[i] * gc_z[i] * tce_0;

        gs_z_y_xxxy[i] = 2.0 * ts_y_xxxy[i] * gc_z[i] * tce_0;

        gs_z_y_xxxz[i] = 2.0 * ts_y_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxxz[i] * gc_z[i] * tce_0;

        gs_z_y_xxyy[i] = 2.0 * ts_y_xxyy[i] * gc_z[i] * tce_0;

        gs_z_y_xxyz[i] = 2.0 * ts_y_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxyz[i] * gc_z[i] * tce_0;

        gs_z_y_xxzz[i] = 4.0 * ts_y_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xxzz[i] * gc_z[i] * tce_0;

        gs_z_y_xyyy[i] = 2.0 * ts_y_xyyy[i] * gc_z[i] * tce_0;

        gs_z_y_xyyz[i] = 2.0 * ts_y_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyyz[i] * gc_z[i] * tce_0;

        gs_z_y_xyzz[i] = 4.0 * ts_y_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xyzz[i] * gc_z[i] * tce_0;

        gs_z_y_xzzz[i] = 6.0 * ts_y_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_xzzz[i] * gc_z[i] * tce_0;

        gs_z_y_yyyy[i] = 2.0 * ts_y_yyyy[i] * gc_z[i] * tce_0;

        gs_z_y_yyyz[i] = 2.0 * ts_y_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyyz[i] * gc_z[i] * tce_0;

        gs_z_y_yyzz[i] = 4.0 * ts_y_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yyzz[i] * gc_z[i] * tce_0;

        gs_z_y_yzzz[i] = 6.0 * ts_y_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_yzzz[i] * gc_z[i] * tce_0;

        gs_z_y_zzzz[i] = 8.0 * ts_y_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_y_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 120-135 components of targeted buffer : PG

    auto gs_z_z_xxxx = pbuffer.data(idx_g_pg + 120);

    auto gs_z_z_xxxy = pbuffer.data(idx_g_pg + 121);

    auto gs_z_z_xxxz = pbuffer.data(idx_g_pg + 122);

    auto gs_z_z_xxyy = pbuffer.data(idx_g_pg + 123);

    auto gs_z_z_xxyz = pbuffer.data(idx_g_pg + 124);

    auto gs_z_z_xxzz = pbuffer.data(idx_g_pg + 125);

    auto gs_z_z_xyyy = pbuffer.data(idx_g_pg + 126);

    auto gs_z_z_xyyz = pbuffer.data(idx_g_pg + 127);

    auto gs_z_z_xyzz = pbuffer.data(idx_g_pg + 128);

    auto gs_z_z_xzzz = pbuffer.data(idx_g_pg + 129);

    auto gs_z_z_yyyy = pbuffer.data(idx_g_pg + 130);

    auto gs_z_z_yyyz = pbuffer.data(idx_g_pg + 131);

    auto gs_z_z_yyzz = pbuffer.data(idx_g_pg + 132);

    auto gs_z_z_yzzz = pbuffer.data(idx_g_pg + 133);

    auto gs_z_z_zzzz = pbuffer.data(idx_g_pg + 134);

    #pragma omp simd aligned(gc_z, gs_z_z_xxxx, gs_z_z_xxxy, gs_z_z_xxxz, gs_z_z_xxyy, gs_z_z_xxyz, gs_z_z_xxzz, gs_z_z_xyyy, gs_z_z_xyyz, gs_z_z_xyzz, gs_z_z_xzzz, gs_z_z_yyyy, gs_z_z_yyyz, gs_z_z_yyzz, gs_z_z_yzzz, gs_z_z_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_z_xxxx[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxx[i] * gc_z[i] * tce_0;

        gs_z_z_xxxy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxy[i] * gc_z[i] * tce_0;

        gs_z_z_xxxz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxxz[i] * gc_z[i] * tce_0;

        gs_z_z_xxyy[i] = 2.0 * ts_0_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyy[i] * gc_z[i] * tce_0;

        gs_z_z_xxyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxyz[i] * gc_z[i] * tce_0;

        gs_z_z_xxzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xxzz[i] * gc_z[i] * tce_0;

        gs_z_z_xyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyy[i] * gc_z[i] * tce_0;

        gs_z_z_xyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyyz[i] * gc_z[i] * tce_0;

        gs_z_z_xyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xyzz[i] * gc_z[i] * tce_0;

        gs_z_z_xzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_z_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_xzzz[i] * gc_z[i] * tce_0;

        gs_z_z_yyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyy[i] * gc_z[i] * tce_0;

        gs_z_z_yyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyyz[i] * gc_z[i] * tce_0;

        gs_z_z_yyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_z_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yyzz[i] * gc_z[i] * tce_0;

        gs_z_z_yzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_z_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_yzzz[i] * gc_z[i] * tce_0;

        gs_z_z_zzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_z_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_z_zzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

