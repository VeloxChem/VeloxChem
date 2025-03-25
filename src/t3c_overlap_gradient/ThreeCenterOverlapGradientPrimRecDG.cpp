#include "ThreeCenterOverlapGradientPrimRecDG.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_dg(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_dg,
                              const size_t idx_pg,
                              const size_t idx_df,
                              const size_t idx_dg,
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

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_df);

    auto ts_xx_xxy = pbuffer.data(idx_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_df + 9);

    auto ts_xy_xxx = pbuffer.data(idx_df + 10);

    auto ts_xy_xxy = pbuffer.data(idx_df + 11);

    auto ts_xy_xxz = pbuffer.data(idx_df + 12);

    auto ts_xy_xyy = pbuffer.data(idx_df + 13);

    auto ts_xy_xyz = pbuffer.data(idx_df + 14);

    auto ts_xy_xzz = pbuffer.data(idx_df + 15);

    auto ts_xy_yyy = pbuffer.data(idx_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_df + 18);

    auto ts_xy_zzz = pbuffer.data(idx_df + 19);

    auto ts_xz_xxx = pbuffer.data(idx_df + 20);

    auto ts_xz_xxy = pbuffer.data(idx_df + 21);

    auto ts_xz_xxz = pbuffer.data(idx_df + 22);

    auto ts_xz_xyy = pbuffer.data(idx_df + 23);

    auto ts_xz_xyz = pbuffer.data(idx_df + 24);

    auto ts_xz_xzz = pbuffer.data(idx_df + 25);

    auto ts_xz_yyy = pbuffer.data(idx_df + 26);

    auto ts_xz_yyz = pbuffer.data(idx_df + 27);

    auto ts_xz_yzz = pbuffer.data(idx_df + 28);

    auto ts_xz_zzz = pbuffer.data(idx_df + 29);

    auto ts_yy_xxx = pbuffer.data(idx_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_df + 39);

    auto ts_yz_xxx = pbuffer.data(idx_df + 40);

    auto ts_yz_xxy = pbuffer.data(idx_df + 41);

    auto ts_yz_xxz = pbuffer.data(idx_df + 42);

    auto ts_yz_xyy = pbuffer.data(idx_df + 43);

    auto ts_yz_xyz = pbuffer.data(idx_df + 44);

    auto ts_yz_xzz = pbuffer.data(idx_df + 45);

    auto ts_yz_yyy = pbuffer.data(idx_df + 46);

    auto ts_yz_yyz = pbuffer.data(idx_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_df + 48);

    auto ts_yz_zzz = pbuffer.data(idx_df + 49);

    auto ts_zz_xxx = pbuffer.data(idx_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_df + 59);

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_dg + 14);

    auto ts_xy_xxxx = pbuffer.data(idx_dg + 15);

    auto ts_xy_xxxy = pbuffer.data(idx_dg + 16);

    auto ts_xy_xxxz = pbuffer.data(idx_dg + 17);

    auto ts_xy_xxyy = pbuffer.data(idx_dg + 18);

    auto ts_xy_xxyz = pbuffer.data(idx_dg + 19);

    auto ts_xy_xxzz = pbuffer.data(idx_dg + 20);

    auto ts_xy_xyyy = pbuffer.data(idx_dg + 21);

    auto ts_xy_xyyz = pbuffer.data(idx_dg + 22);

    auto ts_xy_xyzz = pbuffer.data(idx_dg + 23);

    auto ts_xy_xzzz = pbuffer.data(idx_dg + 24);

    auto ts_xy_yyyy = pbuffer.data(idx_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_dg + 28);

    auto ts_xy_zzzz = pbuffer.data(idx_dg + 29);

    auto ts_xz_xxxx = pbuffer.data(idx_dg + 30);

    auto ts_xz_xxxy = pbuffer.data(idx_dg + 31);

    auto ts_xz_xxxz = pbuffer.data(idx_dg + 32);

    auto ts_xz_xxyy = pbuffer.data(idx_dg + 33);

    auto ts_xz_xxyz = pbuffer.data(idx_dg + 34);

    auto ts_xz_xxzz = pbuffer.data(idx_dg + 35);

    auto ts_xz_xyyy = pbuffer.data(idx_dg + 36);

    auto ts_xz_xyyz = pbuffer.data(idx_dg + 37);

    auto ts_xz_xyzz = pbuffer.data(idx_dg + 38);

    auto ts_xz_xzzz = pbuffer.data(idx_dg + 39);

    auto ts_xz_yyyy = pbuffer.data(idx_dg + 40);

    auto ts_xz_yyyz = pbuffer.data(idx_dg + 41);

    auto ts_xz_yyzz = pbuffer.data(idx_dg + 42);

    auto ts_xz_yzzz = pbuffer.data(idx_dg + 43);

    auto ts_xz_zzzz = pbuffer.data(idx_dg + 44);

    auto ts_yy_xxxx = pbuffer.data(idx_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_dg + 59);

    auto ts_yz_xxxx = pbuffer.data(idx_dg + 60);

    auto ts_yz_xxxy = pbuffer.data(idx_dg + 61);

    auto ts_yz_xxxz = pbuffer.data(idx_dg + 62);

    auto ts_yz_xxyy = pbuffer.data(idx_dg + 63);

    auto ts_yz_xxyz = pbuffer.data(idx_dg + 64);

    auto ts_yz_xxzz = pbuffer.data(idx_dg + 65);

    auto ts_yz_xyyy = pbuffer.data(idx_dg + 66);

    auto ts_yz_xyyz = pbuffer.data(idx_dg + 67);

    auto ts_yz_xyzz = pbuffer.data(idx_dg + 68);

    auto ts_yz_xzzz = pbuffer.data(idx_dg + 69);

    auto ts_yz_yyyy = pbuffer.data(idx_dg + 70);

    auto ts_yz_yyyz = pbuffer.data(idx_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_dg + 73);

    auto ts_yz_zzzz = pbuffer.data(idx_dg + 74);

    auto ts_zz_xxxx = pbuffer.data(idx_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_dg + 89);

    // Set up 0-15 components of targeted buffer : DG

    auto gs_x_xx_xxxx = pbuffer.data(idx_g_dg);

    auto gs_x_xx_xxxy = pbuffer.data(idx_g_dg + 1);

    auto gs_x_xx_xxxz = pbuffer.data(idx_g_dg + 2);

    auto gs_x_xx_xxyy = pbuffer.data(idx_g_dg + 3);

    auto gs_x_xx_xxyz = pbuffer.data(idx_g_dg + 4);

    auto gs_x_xx_xxzz = pbuffer.data(idx_g_dg + 5);

    auto gs_x_xx_xyyy = pbuffer.data(idx_g_dg + 6);

    auto gs_x_xx_xyyz = pbuffer.data(idx_g_dg + 7);

    auto gs_x_xx_xyzz = pbuffer.data(idx_g_dg + 8);

    auto gs_x_xx_xzzz = pbuffer.data(idx_g_dg + 9);

    auto gs_x_xx_yyyy = pbuffer.data(idx_g_dg + 10);

    auto gs_x_xx_yyyz = pbuffer.data(idx_g_dg + 11);

    auto gs_x_xx_yyzz = pbuffer.data(idx_g_dg + 12);

    auto gs_x_xx_yzzz = pbuffer.data(idx_g_dg + 13);

    auto gs_x_xx_zzzz = pbuffer.data(idx_g_dg + 14);

    #pragma omp simd aligned(gc_x, gs_x_xx_xxxx, gs_x_xx_xxxy, gs_x_xx_xxxz, gs_x_xx_xxyy, gs_x_xx_xxyz, gs_x_xx_xxzz, gs_x_xx_xyyy, gs_x_xx_xyyz, gs_x_xx_xyzz, gs_x_xx_xzzz, gs_x_xx_yyyy, gs_x_xx_yyyz, gs_x_xx_yyzz, gs_x_xx_yzzz, gs_x_xx_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xx_xxxx[i] = 4.0 * ts_x_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxy[i] = 4.0 * ts_x_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xx_xxxz[i] = 4.0 * ts_x_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxyy[i] = 4.0 * ts_x_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xx_xxyz[i] = 4.0 * ts_x_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xx_xxzz[i] = 4.0 * ts_x_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xyyy[i] = 4.0 * ts_x_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xx_xyyz[i] = 4.0 * ts_x_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xx_xyzz[i] = 4.0 * ts_x_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xx_xzzz[i] = 4.0 * ts_x_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_yyyy[i] = 4.0 * ts_x_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xx_yyyz[i] = 4.0 * ts_x_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xx_yyzz[i] = 4.0 * ts_x_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xx_yzzz[i] = 4.0 * ts_x_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xx_zzzz[i] = 4.0 * ts_x_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 15-30 components of targeted buffer : DG

    auto gs_x_xy_xxxx = pbuffer.data(idx_g_dg + 15);

    auto gs_x_xy_xxxy = pbuffer.data(idx_g_dg + 16);

    auto gs_x_xy_xxxz = pbuffer.data(idx_g_dg + 17);

    auto gs_x_xy_xxyy = pbuffer.data(idx_g_dg + 18);

    auto gs_x_xy_xxyz = pbuffer.data(idx_g_dg + 19);

    auto gs_x_xy_xxzz = pbuffer.data(idx_g_dg + 20);

    auto gs_x_xy_xyyy = pbuffer.data(idx_g_dg + 21);

    auto gs_x_xy_xyyz = pbuffer.data(idx_g_dg + 22);

    auto gs_x_xy_xyzz = pbuffer.data(idx_g_dg + 23);

    auto gs_x_xy_xzzz = pbuffer.data(idx_g_dg + 24);

    auto gs_x_xy_yyyy = pbuffer.data(idx_g_dg + 25);

    auto gs_x_xy_yyyz = pbuffer.data(idx_g_dg + 26);

    auto gs_x_xy_yyzz = pbuffer.data(idx_g_dg + 27);

    auto gs_x_xy_yzzz = pbuffer.data(idx_g_dg + 28);

    auto gs_x_xy_zzzz = pbuffer.data(idx_g_dg + 29);

    #pragma omp simd aligned(gc_x, gs_x_xy_xxxx, gs_x_xy_xxxy, gs_x_xy_xxxz, gs_x_xy_xxyy, gs_x_xy_xxyz, gs_x_xy_xxzz, gs_x_xy_xyyy, gs_x_xy_xyyz, gs_x_xy_xyzz, gs_x_xy_xzzz, gs_x_xy_yyyy, gs_x_xy_yyyz, gs_x_xy_yyzz, gs_x_xy_yzzz, gs_x_xy_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xy_xxxx[i] = 2.0 * ts_y_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxy[i] = 2.0 * ts_y_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xy_xxxz[i] = 2.0 * ts_y_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxyy[i] = 2.0 * ts_y_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xy_xxyz[i] = 2.0 * ts_y_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xy_xxzz[i] = 2.0 * ts_y_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xyyy[i] = 2.0 * ts_y_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xy_xyyz[i] = 2.0 * ts_y_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xy_xyzz[i] = 2.0 * ts_y_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xy_xzzz[i] = 2.0 * ts_y_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_yyyy[i] = 2.0 * ts_y_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xy_yyyz[i] = 2.0 * ts_y_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xy_yyzz[i] = 2.0 * ts_y_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xy_yzzz[i] = 2.0 * ts_y_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xy_zzzz[i] = 2.0 * ts_y_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-45 components of targeted buffer : DG

    auto gs_x_xz_xxxx = pbuffer.data(idx_g_dg + 30);

    auto gs_x_xz_xxxy = pbuffer.data(idx_g_dg + 31);

    auto gs_x_xz_xxxz = pbuffer.data(idx_g_dg + 32);

    auto gs_x_xz_xxyy = pbuffer.data(idx_g_dg + 33);

    auto gs_x_xz_xxyz = pbuffer.data(idx_g_dg + 34);

    auto gs_x_xz_xxzz = pbuffer.data(idx_g_dg + 35);

    auto gs_x_xz_xyyy = pbuffer.data(idx_g_dg + 36);

    auto gs_x_xz_xyyz = pbuffer.data(idx_g_dg + 37);

    auto gs_x_xz_xyzz = pbuffer.data(idx_g_dg + 38);

    auto gs_x_xz_xzzz = pbuffer.data(idx_g_dg + 39);

    auto gs_x_xz_yyyy = pbuffer.data(idx_g_dg + 40);

    auto gs_x_xz_yyyz = pbuffer.data(idx_g_dg + 41);

    auto gs_x_xz_yyzz = pbuffer.data(idx_g_dg + 42);

    auto gs_x_xz_yzzz = pbuffer.data(idx_g_dg + 43);

    auto gs_x_xz_zzzz = pbuffer.data(idx_g_dg + 44);

    #pragma omp simd aligned(gc_x, gs_x_xz_xxxx, gs_x_xz_xxxy, gs_x_xz_xxxz, gs_x_xz_xxyy, gs_x_xz_xxyz, gs_x_xz_xxzz, gs_x_xz_xyyy, gs_x_xz_xyyz, gs_x_xz_xyzz, gs_x_xz_xzzz, gs_x_xz_yyyy, gs_x_xz_yyyz, gs_x_xz_yyzz, gs_x_xz_yzzz, gs_x_xz_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe_0 * tce_0 + 8.0 * ts_xz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxy[i] = 2.0 * ts_z_xxxy[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_xz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxyy[i] = 2.0 * ts_z_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_xz_xxyz[i] = 2.0 * ts_z_xxyz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_xz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xyyy[i] = 2.0 * ts_z_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_xz_xyyz[i] = 2.0 * ts_z_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_xz_xyzz[i] = 2.0 * ts_z_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_xz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_xz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_xz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_xz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_xz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 45-60 components of targeted buffer : DG

    auto gs_x_yy_xxxx = pbuffer.data(idx_g_dg + 45);

    auto gs_x_yy_xxxy = pbuffer.data(idx_g_dg + 46);

    auto gs_x_yy_xxxz = pbuffer.data(idx_g_dg + 47);

    auto gs_x_yy_xxyy = pbuffer.data(idx_g_dg + 48);

    auto gs_x_yy_xxyz = pbuffer.data(idx_g_dg + 49);

    auto gs_x_yy_xxzz = pbuffer.data(idx_g_dg + 50);

    auto gs_x_yy_xyyy = pbuffer.data(idx_g_dg + 51);

    auto gs_x_yy_xyyz = pbuffer.data(idx_g_dg + 52);

    auto gs_x_yy_xyzz = pbuffer.data(idx_g_dg + 53);

    auto gs_x_yy_xzzz = pbuffer.data(idx_g_dg + 54);

    auto gs_x_yy_yyyy = pbuffer.data(idx_g_dg + 55);

    auto gs_x_yy_yyyz = pbuffer.data(idx_g_dg + 56);

    auto gs_x_yy_yyzz = pbuffer.data(idx_g_dg + 57);

    auto gs_x_yy_yzzz = pbuffer.data(idx_g_dg + 58);

    auto gs_x_yy_zzzz = pbuffer.data(idx_g_dg + 59);

    #pragma omp simd aligned(gc_x, gs_x_yy_xxxx, gs_x_yy_xxxy, gs_x_yy_xxxz, gs_x_yy_xxyy, gs_x_yy_xxyz, gs_x_yy_xxzz, gs_x_yy_xyyy, gs_x_yy_xyyz, gs_x_yy_xyzz, gs_x_yy_xzzz, gs_x_yy_yyyy, gs_x_yy_yyyz, gs_x_yy_yyzz, gs_x_yy_yzzz, gs_x_yy_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yy_xxxx[i] = 8.0 * ts_yy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxy[i] = 6.0 * ts_yy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yy_xxxz[i] = 6.0 * ts_yy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxyy[i] = 4.0 * ts_yy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yy_xxyz[i] = 4.0 * ts_yy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yy_xxzz[i] = 4.0 * ts_yy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xyyy[i] = 2.0 * ts_yy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yy_xyyz[i] = 2.0 * ts_yy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yy_xyzz[i] = 2.0 * ts_yy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yy_xzzz[i] = 2.0 * ts_yy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yy_yyyz[i] = 2.0 * ts_yy_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yy_yyzz[i] = 2.0 * ts_yy_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yy_yzzz[i] = 2.0 * ts_yy_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yy_zzzz[i] = 2.0 * ts_yy_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-75 components of targeted buffer : DG

    auto gs_x_yz_xxxx = pbuffer.data(idx_g_dg + 60);

    auto gs_x_yz_xxxy = pbuffer.data(idx_g_dg + 61);

    auto gs_x_yz_xxxz = pbuffer.data(idx_g_dg + 62);

    auto gs_x_yz_xxyy = pbuffer.data(idx_g_dg + 63);

    auto gs_x_yz_xxyz = pbuffer.data(idx_g_dg + 64);

    auto gs_x_yz_xxzz = pbuffer.data(idx_g_dg + 65);

    auto gs_x_yz_xyyy = pbuffer.data(idx_g_dg + 66);

    auto gs_x_yz_xyyz = pbuffer.data(idx_g_dg + 67);

    auto gs_x_yz_xyzz = pbuffer.data(idx_g_dg + 68);

    auto gs_x_yz_xzzz = pbuffer.data(idx_g_dg + 69);

    auto gs_x_yz_yyyy = pbuffer.data(idx_g_dg + 70);

    auto gs_x_yz_yyyz = pbuffer.data(idx_g_dg + 71);

    auto gs_x_yz_yyzz = pbuffer.data(idx_g_dg + 72);

    auto gs_x_yz_yzzz = pbuffer.data(idx_g_dg + 73);

    auto gs_x_yz_zzzz = pbuffer.data(idx_g_dg + 74);

    #pragma omp simd aligned(gc_x, gs_x_yz_xxxx, gs_x_yz_xxxy, gs_x_yz_xxxz, gs_x_yz_xxyy, gs_x_yz_xxyz, gs_x_yz_xxzz, gs_x_yz_xyyy, gs_x_yz_xyyz, gs_x_yz_xyzz, gs_x_yz_xzzz, gs_x_yz_yyyy, gs_x_yz_yyyz, gs_x_yz_yyzz, gs_x_yz_yzzz, gs_x_yz_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yz_xxxx[i] = 8.0 * ts_yz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxy[i] = 6.0 * ts_yz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_yz_xxxz[i] = 6.0 * ts_yz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxyy[i] = 4.0 * ts_yz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_yz_xxyz[i] = 4.0 * ts_yz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_yz_xxzz[i] = 4.0 * ts_yz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xyyy[i] = 2.0 * ts_yz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_yz_xyyz[i] = 2.0 * ts_yz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_yz_xyzz[i] = 2.0 * ts_yz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_yz_xzzz[i] = 2.0 * ts_yz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_yyyy[i] = 2.0 * ts_yz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_yz_yyyz[i] = 2.0 * ts_yz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_yz_yyzz[i] = 2.0 * ts_yz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_yz_yzzz[i] = 2.0 * ts_yz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_yz_zzzz[i] = 2.0 * ts_yz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 75-90 components of targeted buffer : DG

    auto gs_x_zz_xxxx = pbuffer.data(idx_g_dg + 75);

    auto gs_x_zz_xxxy = pbuffer.data(idx_g_dg + 76);

    auto gs_x_zz_xxxz = pbuffer.data(idx_g_dg + 77);

    auto gs_x_zz_xxyy = pbuffer.data(idx_g_dg + 78);

    auto gs_x_zz_xxyz = pbuffer.data(idx_g_dg + 79);

    auto gs_x_zz_xxzz = pbuffer.data(idx_g_dg + 80);

    auto gs_x_zz_xyyy = pbuffer.data(idx_g_dg + 81);

    auto gs_x_zz_xyyz = pbuffer.data(idx_g_dg + 82);

    auto gs_x_zz_xyzz = pbuffer.data(idx_g_dg + 83);

    auto gs_x_zz_xzzz = pbuffer.data(idx_g_dg + 84);

    auto gs_x_zz_yyyy = pbuffer.data(idx_g_dg + 85);

    auto gs_x_zz_yyyz = pbuffer.data(idx_g_dg + 86);

    auto gs_x_zz_yyzz = pbuffer.data(idx_g_dg + 87);

    auto gs_x_zz_yzzz = pbuffer.data(idx_g_dg + 88);

    auto gs_x_zz_zzzz = pbuffer.data(idx_g_dg + 89);

    #pragma omp simd aligned(gc_x, gs_x_zz_xxxx, gs_x_zz_xxxy, gs_x_zz_xxxz, gs_x_zz_xxyy, gs_x_zz_xxyz, gs_x_zz_xxzz, gs_x_zz_xyyy, gs_x_zz_xyyz, gs_x_zz_xyzz, gs_x_zz_xzzz, gs_x_zz_yyyy, gs_x_zz_yyyz, gs_x_zz_yyzz, gs_x_zz_yzzz, gs_x_zz_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zz_xxxx[i] = 8.0 * ts_zz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxx[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxy[i] = 6.0 * ts_zz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxy[i] * gc_x[i] * tce_0;

        gs_x_zz_xxxz[i] = 6.0 * ts_zz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxyy[i] = 4.0 * ts_zz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyy[i] * gc_x[i] * tce_0;

        gs_x_zz_xxyz[i] = 4.0 * ts_zz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyz[i] * gc_x[i] * tce_0;

        gs_x_zz_xxzz[i] = 4.0 * ts_zz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xyyy[i] = 2.0 * ts_zz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyy[i] * gc_x[i] * tce_0;

        gs_x_zz_xyyz[i] = 2.0 * ts_zz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyz[i] * gc_x[i] * tce_0;

        gs_x_zz_xyzz[i] = 2.0 * ts_zz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyzz[i] * gc_x[i] * tce_0;

        gs_x_zz_xzzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_yyyy[i] = 2.0 * ts_zz_yyyy[i] * gc_x[i] * tce_0;

        gs_x_zz_yyyz[i] = 2.0 * ts_zz_yyyz[i] * gc_x[i] * tce_0;

        gs_x_zz_yyzz[i] = 2.0 * ts_zz_yyzz[i] * gc_x[i] * tce_0;

        gs_x_zz_yzzz[i] = 2.0 * ts_zz_yzzz[i] * gc_x[i] * tce_0;

        gs_x_zz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-105 components of targeted buffer : DG

    auto gs_y_xx_xxxx = pbuffer.data(idx_g_dg + 90);

    auto gs_y_xx_xxxy = pbuffer.data(idx_g_dg + 91);

    auto gs_y_xx_xxxz = pbuffer.data(idx_g_dg + 92);

    auto gs_y_xx_xxyy = pbuffer.data(idx_g_dg + 93);

    auto gs_y_xx_xxyz = pbuffer.data(idx_g_dg + 94);

    auto gs_y_xx_xxzz = pbuffer.data(idx_g_dg + 95);

    auto gs_y_xx_xyyy = pbuffer.data(idx_g_dg + 96);

    auto gs_y_xx_xyyz = pbuffer.data(idx_g_dg + 97);

    auto gs_y_xx_xyzz = pbuffer.data(idx_g_dg + 98);

    auto gs_y_xx_xzzz = pbuffer.data(idx_g_dg + 99);

    auto gs_y_xx_yyyy = pbuffer.data(idx_g_dg + 100);

    auto gs_y_xx_yyyz = pbuffer.data(idx_g_dg + 101);

    auto gs_y_xx_yyzz = pbuffer.data(idx_g_dg + 102);

    auto gs_y_xx_yzzz = pbuffer.data(idx_g_dg + 103);

    auto gs_y_xx_zzzz = pbuffer.data(idx_g_dg + 104);

    #pragma omp simd aligned(gc_y, gs_y_xx_xxxx, gs_y_xx_xxxy, gs_y_xx_xxxz, gs_y_xx_xxyy, gs_y_xx_xxyz, gs_y_xx_xxzz, gs_y_xx_xyyy, gs_y_xx_xyyz, gs_y_xx_xyzz, gs_y_xx_xzzz, gs_y_xx_yyyy, gs_y_xx_yyyz, gs_y_xx_yyzz, gs_y_xx_yzzz, gs_y_xx_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xx_xxxx[i] = 2.0 * ts_xx_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxy[i] = 2.0 * ts_xx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xx_xxxz[i] = 2.0 * ts_xx_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxyy[i] = 4.0 * ts_xx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xx_xxyz[i] = 2.0 * ts_xx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xx_xxzz[i] = 2.0 * ts_xx_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xyyy[i] = 6.0 * ts_xx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xx_xyyz[i] = 4.0 * ts_xx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xx_xyzz[i] = 2.0 * ts_xx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xx_xzzz[i] = 2.0 * ts_xx_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_yyyy[i] = 8.0 * ts_xx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xx_yyyz[i] = 6.0 * ts_xx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xx_yyzz[i] = 4.0 * ts_xx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xx_yzzz[i] = 2.0 * ts_xx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xx_zzzz[i] = 2.0 * ts_xx_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 105-120 components of targeted buffer : DG

    auto gs_y_xy_xxxx = pbuffer.data(idx_g_dg + 105);

    auto gs_y_xy_xxxy = pbuffer.data(idx_g_dg + 106);

    auto gs_y_xy_xxxz = pbuffer.data(idx_g_dg + 107);

    auto gs_y_xy_xxyy = pbuffer.data(idx_g_dg + 108);

    auto gs_y_xy_xxyz = pbuffer.data(idx_g_dg + 109);

    auto gs_y_xy_xxzz = pbuffer.data(idx_g_dg + 110);

    auto gs_y_xy_xyyy = pbuffer.data(idx_g_dg + 111);

    auto gs_y_xy_xyyz = pbuffer.data(idx_g_dg + 112);

    auto gs_y_xy_xyzz = pbuffer.data(idx_g_dg + 113);

    auto gs_y_xy_xzzz = pbuffer.data(idx_g_dg + 114);

    auto gs_y_xy_yyyy = pbuffer.data(idx_g_dg + 115);

    auto gs_y_xy_yyyz = pbuffer.data(idx_g_dg + 116);

    auto gs_y_xy_yyzz = pbuffer.data(idx_g_dg + 117);

    auto gs_y_xy_yzzz = pbuffer.data(idx_g_dg + 118);

    auto gs_y_xy_zzzz = pbuffer.data(idx_g_dg + 119);

    #pragma omp simd aligned(gc_y, gs_y_xy_xxxx, gs_y_xy_xxxy, gs_y_xy_xxxz, gs_y_xy_xxyy, gs_y_xy_xxyz, gs_y_xy_xxzz, gs_y_xy_xyyy, gs_y_xy_xyyz, gs_y_xy_xyzz, gs_y_xy_xzzz, gs_y_xy_yyyy, gs_y_xy_yyyz, gs_y_xy_yyzz, gs_y_xy_yzzz, gs_y_xy_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xy_xxxx[i] = 2.0 * ts_x_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxy[i] = 2.0 * ts_x_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xy_xxxz[i] = 2.0 * ts_x_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxyy[i] = 2.0 * ts_x_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xy_xxyz[i] = 2.0 * ts_x_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xy_xxzz[i] = 2.0 * ts_x_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xyyy[i] = 2.0 * ts_x_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_xy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xy_xyyz[i] = 2.0 * ts_x_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xy_xyzz[i] = 2.0 * ts_x_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xy_xzzz[i] = 2.0 * ts_x_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_yyyy[i] = 2.0 * ts_x_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_xy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xy_yyyz[i] = 2.0 * ts_x_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_xy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xy_yyzz[i] = 2.0 * ts_x_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xy_yzzz[i] = 2.0 * ts_x_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xy_zzzz[i] = 2.0 * ts_x_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 120-135 components of targeted buffer : DG

    auto gs_y_xz_xxxx = pbuffer.data(idx_g_dg + 120);

    auto gs_y_xz_xxxy = pbuffer.data(idx_g_dg + 121);

    auto gs_y_xz_xxxz = pbuffer.data(idx_g_dg + 122);

    auto gs_y_xz_xxyy = pbuffer.data(idx_g_dg + 123);

    auto gs_y_xz_xxyz = pbuffer.data(idx_g_dg + 124);

    auto gs_y_xz_xxzz = pbuffer.data(idx_g_dg + 125);

    auto gs_y_xz_xyyy = pbuffer.data(idx_g_dg + 126);

    auto gs_y_xz_xyyz = pbuffer.data(idx_g_dg + 127);

    auto gs_y_xz_xyzz = pbuffer.data(idx_g_dg + 128);

    auto gs_y_xz_xzzz = pbuffer.data(idx_g_dg + 129);

    auto gs_y_xz_yyyy = pbuffer.data(idx_g_dg + 130);

    auto gs_y_xz_yyyz = pbuffer.data(idx_g_dg + 131);

    auto gs_y_xz_yyzz = pbuffer.data(idx_g_dg + 132);

    auto gs_y_xz_yzzz = pbuffer.data(idx_g_dg + 133);

    auto gs_y_xz_zzzz = pbuffer.data(idx_g_dg + 134);

    #pragma omp simd aligned(gc_y, gs_y_xz_xxxx, gs_y_xz_xxxy, gs_y_xz_xxxz, gs_y_xz_xxyy, gs_y_xz_xxyz, gs_y_xz_xxzz, gs_y_xz_xyyy, gs_y_xz_xyyz, gs_y_xz_xyzz, gs_y_xz_xzzz, gs_y_xz_yyyy, gs_y_xz_yyyz, gs_y_xz_yyzz, gs_y_xz_yzzz, gs_y_xz_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xz_xxxx[i] = 2.0 * ts_xz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxy[i] = 2.0 * ts_xz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_xz_xxxz[i] = 2.0 * ts_xz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxyy[i] = 4.0 * ts_xz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_xz_xxyz[i] = 2.0 * ts_xz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_xz_xxzz[i] = 2.0 * ts_xz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xyyy[i] = 6.0 * ts_xz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_xz_xyyz[i] = 4.0 * ts_xz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_xz_xyzz[i] = 2.0 * ts_xz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_xz_xzzz[i] = 2.0 * ts_xz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_yyyy[i] = 8.0 * ts_xz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_xz_yyyz[i] = 6.0 * ts_xz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_xz_yyzz[i] = 4.0 * ts_xz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_xz_yzzz[i] = 2.0 * ts_xz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_xz_zzzz[i] = 2.0 * ts_xz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 135-150 components of targeted buffer : DG

    auto gs_y_yy_xxxx = pbuffer.data(idx_g_dg + 135);

    auto gs_y_yy_xxxy = pbuffer.data(idx_g_dg + 136);

    auto gs_y_yy_xxxz = pbuffer.data(idx_g_dg + 137);

    auto gs_y_yy_xxyy = pbuffer.data(idx_g_dg + 138);

    auto gs_y_yy_xxyz = pbuffer.data(idx_g_dg + 139);

    auto gs_y_yy_xxzz = pbuffer.data(idx_g_dg + 140);

    auto gs_y_yy_xyyy = pbuffer.data(idx_g_dg + 141);

    auto gs_y_yy_xyyz = pbuffer.data(idx_g_dg + 142);

    auto gs_y_yy_xyzz = pbuffer.data(idx_g_dg + 143);

    auto gs_y_yy_xzzz = pbuffer.data(idx_g_dg + 144);

    auto gs_y_yy_yyyy = pbuffer.data(idx_g_dg + 145);

    auto gs_y_yy_yyyz = pbuffer.data(idx_g_dg + 146);

    auto gs_y_yy_yyzz = pbuffer.data(idx_g_dg + 147);

    auto gs_y_yy_yzzz = pbuffer.data(idx_g_dg + 148);

    auto gs_y_yy_zzzz = pbuffer.data(idx_g_dg + 149);

    #pragma omp simd aligned(gc_y, gs_y_yy_xxxx, gs_y_yy_xxxy, gs_y_yy_xxxz, gs_y_yy_xxyy, gs_y_yy_xxyz, gs_y_yy_xxzz, gs_y_yy_xyyy, gs_y_yy_xyyz, gs_y_yy_xyzz, gs_y_yy_xzzz, gs_y_yy_yyyy, gs_y_yy_yyyz, gs_y_yy_yyzz, gs_y_yy_yzzz, gs_y_yy_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yy_xxxx[i] = 4.0 * ts_y_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxy[i] = 4.0 * ts_y_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yy_xxxz[i] = 4.0 * ts_y_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxyy[i] = 4.0 * ts_y_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yy_xxyz[i] = 4.0 * ts_y_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yy_xxzz[i] = 4.0 * ts_y_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xyyy[i] = 4.0 * ts_y_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yy_xyyz[i] = 4.0 * ts_y_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yy_xyzz[i] = 4.0 * ts_y_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yy_xzzz[i] = 4.0 * ts_y_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_yyyy[i] = 4.0 * ts_y_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yy_yyyz[i] = 4.0 * ts_y_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yy_yyzz[i] = 4.0 * ts_y_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yy_yzzz[i] = 4.0 * ts_y_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yy_zzzz[i] = 4.0 * ts_y_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 150-165 components of targeted buffer : DG

    auto gs_y_yz_xxxx = pbuffer.data(idx_g_dg + 150);

    auto gs_y_yz_xxxy = pbuffer.data(idx_g_dg + 151);

    auto gs_y_yz_xxxz = pbuffer.data(idx_g_dg + 152);

    auto gs_y_yz_xxyy = pbuffer.data(idx_g_dg + 153);

    auto gs_y_yz_xxyz = pbuffer.data(idx_g_dg + 154);

    auto gs_y_yz_xxzz = pbuffer.data(idx_g_dg + 155);

    auto gs_y_yz_xyyy = pbuffer.data(idx_g_dg + 156);

    auto gs_y_yz_xyyz = pbuffer.data(idx_g_dg + 157);

    auto gs_y_yz_xyzz = pbuffer.data(idx_g_dg + 158);

    auto gs_y_yz_xzzz = pbuffer.data(idx_g_dg + 159);

    auto gs_y_yz_yyyy = pbuffer.data(idx_g_dg + 160);

    auto gs_y_yz_yyyz = pbuffer.data(idx_g_dg + 161);

    auto gs_y_yz_yyzz = pbuffer.data(idx_g_dg + 162);

    auto gs_y_yz_yzzz = pbuffer.data(idx_g_dg + 163);

    auto gs_y_yz_zzzz = pbuffer.data(idx_g_dg + 164);

    #pragma omp simd aligned(gc_y, gs_y_yz_xxxx, gs_y_yz_xxxy, gs_y_yz_xxxz, gs_y_yz_xxyy, gs_y_yz_xxyz, gs_y_yz_xxzz, gs_y_yz_xyyy, gs_y_yz_xyyz, gs_y_yz_xyzz, gs_y_yz_xzzz, gs_y_yz_yyyy, gs_y_yz_yyyz, gs_y_yz_yyzz, gs_y_yz_yzzz, gs_y_yz_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yz_xxxx[i] = 2.0 * ts_z_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxy[i] = 2.0 * ts_z_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_yz_xxxz[i] = 2.0 * ts_z_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxyy[i] = 2.0 * ts_z_xxyy[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_yz_xxyz[i] = 2.0 * ts_z_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_yz_xxzz[i] = 2.0 * ts_z_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xyyy[i] = 2.0 * ts_z_xyyy[i] * gfe_0 * tce_0 + 6.0 * ts_yz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_yz_xyyz[i] = 2.0 * ts_z_xyyz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_yz_xyzz[i] = 2.0 * ts_z_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_yz_xzzz[i] = 2.0 * ts_z_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_yyyy[i] = 2.0 * ts_z_yyyy[i] * gfe_0 * tce_0 + 8.0 * ts_yz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_yz_yyyz[i] = 2.0 * ts_z_yyyz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_yz_yyzz[i] = 2.0 * ts_z_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_yz_yzzz[i] = 2.0 * ts_z_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_yz_zzzz[i] = 2.0 * ts_z_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 165-180 components of targeted buffer : DG

    auto gs_y_zz_xxxx = pbuffer.data(idx_g_dg + 165);

    auto gs_y_zz_xxxy = pbuffer.data(idx_g_dg + 166);

    auto gs_y_zz_xxxz = pbuffer.data(idx_g_dg + 167);

    auto gs_y_zz_xxyy = pbuffer.data(idx_g_dg + 168);

    auto gs_y_zz_xxyz = pbuffer.data(idx_g_dg + 169);

    auto gs_y_zz_xxzz = pbuffer.data(idx_g_dg + 170);

    auto gs_y_zz_xyyy = pbuffer.data(idx_g_dg + 171);

    auto gs_y_zz_xyyz = pbuffer.data(idx_g_dg + 172);

    auto gs_y_zz_xyzz = pbuffer.data(idx_g_dg + 173);

    auto gs_y_zz_xzzz = pbuffer.data(idx_g_dg + 174);

    auto gs_y_zz_yyyy = pbuffer.data(idx_g_dg + 175);

    auto gs_y_zz_yyyz = pbuffer.data(idx_g_dg + 176);

    auto gs_y_zz_yyzz = pbuffer.data(idx_g_dg + 177);

    auto gs_y_zz_yzzz = pbuffer.data(idx_g_dg + 178);

    auto gs_y_zz_zzzz = pbuffer.data(idx_g_dg + 179);

    #pragma omp simd aligned(gc_y, gs_y_zz_xxxx, gs_y_zz_xxxy, gs_y_zz_xxxz, gs_y_zz_xxyy, gs_y_zz_xxyz, gs_y_zz_xxzz, gs_y_zz_xyyy, gs_y_zz_xyyz, gs_y_zz_xyzz, gs_y_zz_xzzz, gs_y_zz_yyyy, gs_y_zz_yyyz, gs_y_zz_yyzz, gs_y_zz_yzzz, gs_y_zz_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zz_xxxx[i] = 2.0 * ts_zz_xxxx[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxy[i] = 2.0 * ts_zz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxy[i] * gc_y[i] * tce_0;

        gs_y_zz_xxxz[i] = 2.0 * ts_zz_xxxz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxyy[i] = 4.0 * ts_zz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyy[i] * gc_y[i] * tce_0;

        gs_y_zz_xxyz[i] = 2.0 * ts_zz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyz[i] * gc_y[i] * tce_0;

        gs_y_zz_xxzz[i] = 2.0 * ts_zz_xxzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xyyy[i] = 6.0 * ts_zz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyy[i] * gc_y[i] * tce_0;

        gs_y_zz_xyyz[i] = 4.0 * ts_zz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyz[i] * gc_y[i] * tce_0;

        gs_y_zz_xyzz[i] = 2.0 * ts_zz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyzz[i] * gc_y[i] * tce_0;

        gs_y_zz_xzzz[i] = 2.0 * ts_zz_xzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_yyyy[i] = 8.0 * ts_zz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyy[i] * gc_y[i] * tce_0;

        gs_y_zz_yyyz[i] = 6.0 * ts_zz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyz[i] * gc_y[i] * tce_0;

        gs_y_zz_yyzz[i] = 4.0 * ts_zz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyzz[i] * gc_y[i] * tce_0;

        gs_y_zz_yzzz[i] = 2.0 * ts_zz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yzzz[i] * gc_y[i] * tce_0;

        gs_y_zz_zzzz[i] = 2.0 * ts_zz_zzzz[i] * gc_y[i] * tce_0;
    }

    // Set up 180-195 components of targeted buffer : DG

    auto gs_z_xx_xxxx = pbuffer.data(idx_g_dg + 180);

    auto gs_z_xx_xxxy = pbuffer.data(idx_g_dg + 181);

    auto gs_z_xx_xxxz = pbuffer.data(idx_g_dg + 182);

    auto gs_z_xx_xxyy = pbuffer.data(idx_g_dg + 183);

    auto gs_z_xx_xxyz = pbuffer.data(idx_g_dg + 184);

    auto gs_z_xx_xxzz = pbuffer.data(idx_g_dg + 185);

    auto gs_z_xx_xyyy = pbuffer.data(idx_g_dg + 186);

    auto gs_z_xx_xyyz = pbuffer.data(idx_g_dg + 187);

    auto gs_z_xx_xyzz = pbuffer.data(idx_g_dg + 188);

    auto gs_z_xx_xzzz = pbuffer.data(idx_g_dg + 189);

    auto gs_z_xx_yyyy = pbuffer.data(idx_g_dg + 190);

    auto gs_z_xx_yyyz = pbuffer.data(idx_g_dg + 191);

    auto gs_z_xx_yyzz = pbuffer.data(idx_g_dg + 192);

    auto gs_z_xx_yzzz = pbuffer.data(idx_g_dg + 193);

    auto gs_z_xx_zzzz = pbuffer.data(idx_g_dg + 194);

    #pragma omp simd aligned(gc_z, gs_z_xx_xxxx, gs_z_xx_xxxy, gs_z_xx_xxxz, gs_z_xx_xxyy, gs_z_xx_xxyz, gs_z_xx_xxzz, gs_z_xx_xyyy, gs_z_xx_xyyz, gs_z_xx_xyzz, gs_z_xx_xzzz, gs_z_xx_yyyy, gs_z_xx_yyyz, gs_z_xx_yyzz, gs_z_xx_yzzz, gs_z_xx_zzzz, ts_xx_xxx, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxy, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxz, ts_xx_xxzz, ts_xx_xyy, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyz, ts_xx_xyzz, ts_xx_xzz, ts_xx_xzzz, ts_xx_yyy, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyz, ts_xx_yyzz, ts_xx_yzz, ts_xx_yzzz, ts_xx_zzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xx_xxxx[i] = 2.0 * ts_xx_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxy[i] = 2.0 * ts_xx_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xx_xxxz[i] = 2.0 * ts_xx_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxyy[i] = 2.0 * ts_xx_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xx_xxyz[i] = 2.0 * ts_xx_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xx_xxzz[i] = 4.0 * ts_xx_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xyyy[i] = 2.0 * ts_xx_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xx_xyyz[i] = 2.0 * ts_xx_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xx_xyzz[i] = 4.0 * ts_xx_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xx_xzzz[i] = 6.0 * ts_xx_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_yyyy[i] = 2.0 * ts_xx_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xx_yyyz[i] = 2.0 * ts_xx_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xx_yyzz[i] = 4.0 * ts_xx_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xx_yzzz[i] = 6.0 * ts_xx_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xx_zzzz[i] = 8.0 * ts_xx_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xx_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 195-210 components of targeted buffer : DG

    auto gs_z_xy_xxxx = pbuffer.data(idx_g_dg + 195);

    auto gs_z_xy_xxxy = pbuffer.data(idx_g_dg + 196);

    auto gs_z_xy_xxxz = pbuffer.data(idx_g_dg + 197);

    auto gs_z_xy_xxyy = pbuffer.data(idx_g_dg + 198);

    auto gs_z_xy_xxyz = pbuffer.data(idx_g_dg + 199);

    auto gs_z_xy_xxzz = pbuffer.data(idx_g_dg + 200);

    auto gs_z_xy_xyyy = pbuffer.data(idx_g_dg + 201);

    auto gs_z_xy_xyyz = pbuffer.data(idx_g_dg + 202);

    auto gs_z_xy_xyzz = pbuffer.data(idx_g_dg + 203);

    auto gs_z_xy_xzzz = pbuffer.data(idx_g_dg + 204);

    auto gs_z_xy_yyyy = pbuffer.data(idx_g_dg + 205);

    auto gs_z_xy_yyyz = pbuffer.data(idx_g_dg + 206);

    auto gs_z_xy_yyzz = pbuffer.data(idx_g_dg + 207);

    auto gs_z_xy_yzzz = pbuffer.data(idx_g_dg + 208);

    auto gs_z_xy_zzzz = pbuffer.data(idx_g_dg + 209);

    #pragma omp simd aligned(gc_z, gs_z_xy_xxxx, gs_z_xy_xxxy, gs_z_xy_xxxz, gs_z_xy_xxyy, gs_z_xy_xxyz, gs_z_xy_xxzz, gs_z_xy_xyyy, gs_z_xy_xyyz, gs_z_xy_xyzz, gs_z_xy_xzzz, gs_z_xy_yyyy, gs_z_xy_yyyz, gs_z_xy_yyzz, gs_z_xy_yzzz, gs_z_xy_zzzz, ts_xy_xxx, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxy, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxz, ts_xy_xxzz, ts_xy_xyy, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyz, ts_xy_xyzz, ts_xy_xzz, ts_xy_xzzz, ts_xy_yyy, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyz, ts_xy_yyzz, ts_xy_yzz, ts_xy_yzzz, ts_xy_zzz, ts_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xy_xxxx[i] = 2.0 * ts_xy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxy[i] = 2.0 * ts_xy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xy_xxxz[i] = 2.0 * ts_xy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxyy[i] = 2.0 * ts_xy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xy_xxyz[i] = 2.0 * ts_xy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xy_xxzz[i] = 4.0 * ts_xy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xyyy[i] = 2.0 * ts_xy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xy_xyyz[i] = 2.0 * ts_xy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xy_xyzz[i] = 4.0 * ts_xy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xy_xzzz[i] = 6.0 * ts_xy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_yyyy[i] = 2.0 * ts_xy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xy_yyyz[i] = 2.0 * ts_xy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xy_yyzz[i] = 4.0 * ts_xy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xy_yzzz[i] = 6.0 * ts_xy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xy_zzzz[i] = 8.0 * ts_xy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 210-225 components of targeted buffer : DG

    auto gs_z_xz_xxxx = pbuffer.data(idx_g_dg + 210);

    auto gs_z_xz_xxxy = pbuffer.data(idx_g_dg + 211);

    auto gs_z_xz_xxxz = pbuffer.data(idx_g_dg + 212);

    auto gs_z_xz_xxyy = pbuffer.data(idx_g_dg + 213);

    auto gs_z_xz_xxyz = pbuffer.data(idx_g_dg + 214);

    auto gs_z_xz_xxzz = pbuffer.data(idx_g_dg + 215);

    auto gs_z_xz_xyyy = pbuffer.data(idx_g_dg + 216);

    auto gs_z_xz_xyyz = pbuffer.data(idx_g_dg + 217);

    auto gs_z_xz_xyzz = pbuffer.data(idx_g_dg + 218);

    auto gs_z_xz_xzzz = pbuffer.data(idx_g_dg + 219);

    auto gs_z_xz_yyyy = pbuffer.data(idx_g_dg + 220);

    auto gs_z_xz_yyyz = pbuffer.data(idx_g_dg + 221);

    auto gs_z_xz_yyzz = pbuffer.data(idx_g_dg + 222);

    auto gs_z_xz_yzzz = pbuffer.data(idx_g_dg + 223);

    auto gs_z_xz_zzzz = pbuffer.data(idx_g_dg + 224);

    #pragma omp simd aligned(gc_z, gs_z_xz_xxxx, gs_z_xz_xxxy, gs_z_xz_xxxz, gs_z_xz_xxyy, gs_z_xz_xxyz, gs_z_xz_xxzz, gs_z_xz_xyyy, gs_z_xz_xyyz, gs_z_xz_xyzz, gs_z_xz_xzzz, gs_z_xz_yyyy, gs_z_xz_yyyz, gs_z_xz_yyzz, gs_z_xz_yzzz, gs_z_xz_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, ts_xz_xxx, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxy, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxz, ts_xz_xxzz, ts_xz_xyy, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyz, ts_xz_xyzz, ts_xz_xzz, ts_xz_xzzz, ts_xz_yyy, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyz, ts_xz_yyzz, ts_xz_yzz, ts_xz_yzzz, ts_xz_zzz, ts_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xz_xxxx[i] = 2.0 * ts_x_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxy[i] = 2.0 * ts_x_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_xz_xxxz[i] = 2.0 * ts_x_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxyy[i] = 2.0 * ts_x_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_xz_xxyz[i] = 2.0 * ts_x_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_xz_xxzz[i] = 2.0 * ts_x_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xyyy[i] = 2.0 * ts_x_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_xz_xyyz[i] = 2.0 * ts_x_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_xz_xyzz[i] = 2.0 * ts_x_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_xz_xzzz[i] = 2.0 * ts_x_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_yyyy[i] = 2.0 * ts_x_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_xz_yyyz[i] = 2.0 * ts_x_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_xz_yyzz[i] = 2.0 * ts_x_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_xz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_xz_yzzz[i] = 2.0 * ts_x_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_xz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_xz_zzzz[i] = 2.0 * ts_x_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_xz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_xz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 225-240 components of targeted buffer : DG

    auto gs_z_yy_xxxx = pbuffer.data(idx_g_dg + 225);

    auto gs_z_yy_xxxy = pbuffer.data(idx_g_dg + 226);

    auto gs_z_yy_xxxz = pbuffer.data(idx_g_dg + 227);

    auto gs_z_yy_xxyy = pbuffer.data(idx_g_dg + 228);

    auto gs_z_yy_xxyz = pbuffer.data(idx_g_dg + 229);

    auto gs_z_yy_xxzz = pbuffer.data(idx_g_dg + 230);

    auto gs_z_yy_xyyy = pbuffer.data(idx_g_dg + 231);

    auto gs_z_yy_xyyz = pbuffer.data(idx_g_dg + 232);

    auto gs_z_yy_xyzz = pbuffer.data(idx_g_dg + 233);

    auto gs_z_yy_xzzz = pbuffer.data(idx_g_dg + 234);

    auto gs_z_yy_yyyy = pbuffer.data(idx_g_dg + 235);

    auto gs_z_yy_yyyz = pbuffer.data(idx_g_dg + 236);

    auto gs_z_yy_yyzz = pbuffer.data(idx_g_dg + 237);

    auto gs_z_yy_yzzz = pbuffer.data(idx_g_dg + 238);

    auto gs_z_yy_zzzz = pbuffer.data(idx_g_dg + 239);

    #pragma omp simd aligned(gc_z, gs_z_yy_xxxx, gs_z_yy_xxxy, gs_z_yy_xxxz, gs_z_yy_xxyy, gs_z_yy_xxyz, gs_z_yy_xxzz, gs_z_yy_xyyy, gs_z_yy_xyyz, gs_z_yy_xyzz, gs_z_yy_xzzz, gs_z_yy_yyyy, gs_z_yy_yyyz, gs_z_yy_yyzz, gs_z_yy_yzzz, gs_z_yy_zzzz, ts_yy_xxx, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxy, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxz, ts_yy_xxzz, ts_yy_xyy, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyz, ts_yy_xyzz, ts_yy_xzz, ts_yy_xzzz, ts_yy_yyy, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyz, ts_yy_yyzz, ts_yy_yzz, ts_yy_yzzz, ts_yy_zzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yy_xxxx[i] = 2.0 * ts_yy_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxy[i] = 2.0 * ts_yy_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yy_xxxz[i] = 2.0 * ts_yy_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxyy[i] = 2.0 * ts_yy_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yy_xxyz[i] = 2.0 * ts_yy_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yy_xxzz[i] = 4.0 * ts_yy_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xyyy[i] = 2.0 * ts_yy_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yy_xyyz[i] = 2.0 * ts_yy_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yy_xyzz[i] = 4.0 * ts_yy_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yy_xzzz[i] = 6.0 * ts_yy_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_yyyy[i] = 2.0 * ts_yy_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yy_yyyz[i] = 2.0 * ts_yy_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yy_yyzz[i] = 4.0 * ts_yy_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yy_yzzz[i] = 6.0 * ts_yy_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yy_zzzz[i] = 8.0 * ts_yy_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yy_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 240-255 components of targeted buffer : DG

    auto gs_z_yz_xxxx = pbuffer.data(idx_g_dg + 240);

    auto gs_z_yz_xxxy = pbuffer.data(idx_g_dg + 241);

    auto gs_z_yz_xxxz = pbuffer.data(idx_g_dg + 242);

    auto gs_z_yz_xxyy = pbuffer.data(idx_g_dg + 243);

    auto gs_z_yz_xxyz = pbuffer.data(idx_g_dg + 244);

    auto gs_z_yz_xxzz = pbuffer.data(idx_g_dg + 245);

    auto gs_z_yz_xyyy = pbuffer.data(idx_g_dg + 246);

    auto gs_z_yz_xyyz = pbuffer.data(idx_g_dg + 247);

    auto gs_z_yz_xyzz = pbuffer.data(idx_g_dg + 248);

    auto gs_z_yz_xzzz = pbuffer.data(idx_g_dg + 249);

    auto gs_z_yz_yyyy = pbuffer.data(idx_g_dg + 250);

    auto gs_z_yz_yyyz = pbuffer.data(idx_g_dg + 251);

    auto gs_z_yz_yyzz = pbuffer.data(idx_g_dg + 252);

    auto gs_z_yz_yzzz = pbuffer.data(idx_g_dg + 253);

    auto gs_z_yz_zzzz = pbuffer.data(idx_g_dg + 254);

    #pragma omp simd aligned(gc_z, gs_z_yz_xxxx, gs_z_yz_xxxy, gs_z_yz_xxxz, gs_z_yz_xxyy, gs_z_yz_xxyz, gs_z_yz_xxzz, gs_z_yz_xyyy, gs_z_yz_xyyz, gs_z_yz_xyzz, gs_z_yz_xzzz, gs_z_yz_yyyy, gs_z_yz_yyyz, gs_z_yz_yyzz, gs_z_yz_yzzz, gs_z_yz_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, ts_yz_xxx, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxy, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxz, ts_yz_xxzz, ts_yz_xyy, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyz, ts_yz_xyzz, ts_yz_xzz, ts_yz_xzzz, ts_yz_yyy, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyz, ts_yz_yyzz, ts_yz_yzz, ts_yz_yzzz, ts_yz_zzz, ts_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yz_xxxx[i] = 2.0 * ts_y_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxy[i] = 2.0 * ts_y_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_yz_xxxz[i] = 2.0 * ts_y_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxyy[i] = 2.0 * ts_y_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_yz_xxyz[i] = 2.0 * ts_y_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_yz_xxzz[i] = 2.0 * ts_y_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xyyy[i] = 2.0 * ts_y_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_yz_xyyz[i] = 2.0 * ts_y_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_yz_xyzz[i] = 2.0 * ts_y_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_yz_xzzz[i] = 2.0 * ts_y_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_yyyy[i] = 2.0 * ts_y_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_yz_yyyz[i] = 2.0 * ts_y_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_yz_yyzz[i] = 2.0 * ts_y_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_yz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_yz_yzzz[i] = 2.0 * ts_y_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_yz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_yz_zzzz[i] = 2.0 * ts_y_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_yz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_yz_zzzz[i] * gc_z[i] * tce_0;
    }

    // Set up 255-270 components of targeted buffer : DG

    auto gs_z_zz_xxxx = pbuffer.data(idx_g_dg + 255);

    auto gs_z_zz_xxxy = pbuffer.data(idx_g_dg + 256);

    auto gs_z_zz_xxxz = pbuffer.data(idx_g_dg + 257);

    auto gs_z_zz_xxyy = pbuffer.data(idx_g_dg + 258);

    auto gs_z_zz_xxyz = pbuffer.data(idx_g_dg + 259);

    auto gs_z_zz_xxzz = pbuffer.data(idx_g_dg + 260);

    auto gs_z_zz_xyyy = pbuffer.data(idx_g_dg + 261);

    auto gs_z_zz_xyyz = pbuffer.data(idx_g_dg + 262);

    auto gs_z_zz_xyzz = pbuffer.data(idx_g_dg + 263);

    auto gs_z_zz_xzzz = pbuffer.data(idx_g_dg + 264);

    auto gs_z_zz_yyyy = pbuffer.data(idx_g_dg + 265);

    auto gs_z_zz_yyyz = pbuffer.data(idx_g_dg + 266);

    auto gs_z_zz_yyzz = pbuffer.data(idx_g_dg + 267);

    auto gs_z_zz_yzzz = pbuffer.data(idx_g_dg + 268);

    auto gs_z_zz_zzzz = pbuffer.data(idx_g_dg + 269);

    #pragma omp simd aligned(gc_z, gs_z_zz_xxxx, gs_z_zz_xxxy, gs_z_zz_xxxz, gs_z_zz_xxyy, gs_z_zz_xxyz, gs_z_zz_xxzz, gs_z_zz_xyyy, gs_z_zz_xyyz, gs_z_zz_xyzz, gs_z_zz_xzzz, gs_z_zz_yyyy, gs_z_zz_yyyz, gs_z_zz_yyzz, gs_z_zz_yzzz, gs_z_zz_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, ts_zz_xxx, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxy, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxz, ts_zz_xxzz, ts_zz_xyy, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyz, ts_zz_xyzz, ts_zz_xzz, ts_zz_xzzz, ts_zz_yyy, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyz, ts_zz_yyzz, ts_zz_yzz, ts_zz_yzzz, ts_zz_zzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zz_xxxx[i] = 4.0 * ts_z_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxx[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxy[i] = 4.0 * ts_z_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxy[i] * gc_z[i] * tce_0;

        gs_z_zz_xxxz[i] = 4.0 * ts_z_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxx[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxxz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxyy[i] = 4.0 * ts_z_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyy[i] * gc_z[i] * tce_0;

        gs_z_zz_xxyz[i] = 4.0 * ts_z_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxyz[i] * gc_z[i] * tce_0;

        gs_z_zz_xxzz[i] = 4.0 * ts_z_xxzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_xxz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xxzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xyyy[i] = 4.0 * ts_z_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyy[i] * gc_z[i] * tce_0;

        gs_z_zz_xyyz[i] = 4.0 * ts_z_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyyz[i] * gc_z[i] * tce_0;

        gs_z_zz_xyzz[i] = 4.0 * ts_z_xyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_xyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xyzz[i] * gc_z[i] * tce_0;

        gs_z_zz_xzzz[i] = 4.0 * ts_z_xzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zz_xzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_xzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_yyyy[i] = 4.0 * ts_z_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyy[i] * gc_z[i] * tce_0;

        gs_z_zz_yyyz[i] = 4.0 * ts_z_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyy[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyyz[i] * gc_z[i] * tce_0;

        gs_z_zz_yyzz[i] = 4.0 * ts_z_yyzz[i] * gfe_0 * tce_0 + 4.0 * ts_zz_yyz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yyzz[i] * gc_z[i] * tce_0;

        gs_z_zz_yzzz[i] = 4.0 * ts_z_yzzz[i] * gfe_0 * tce_0 + 6.0 * ts_zz_yzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_yzzz[i] * gc_z[i] * tce_0;

        gs_z_zz_zzzz[i] = 4.0 * ts_z_zzzz[i] * gfe_0 * tce_0 + 8.0 * ts_zz_zzz[i] * gfe_0 * tce_0 + 2.0 * ts_zz_zzzz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

