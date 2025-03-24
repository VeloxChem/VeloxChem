#include "ThreeCenterOverlapPrimRecDG.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_dg(CSimdArray<double>& pbuffer, 
                     const size_t idx_dg,
                     const size_t idx_sg,
                     const size_t idx_pf,
                     const size_t idx_pg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

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

    // Set up 0-15 components of targeted buffer : DG

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

    #pragma omp simd aligned(ga_x, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_x_xxx, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxy, ts_x_xxyy, ts_x_xxyz, ts_x_xxz, ts_x_xxzz, ts_x_xyy, ts_x_xyyy, ts_x_xyyz, ts_x_xyz, ts_x_xyzz, ts_x_xzz, ts_x_xzzz, ts_x_yyy, ts_x_yyyy, ts_x_yyyz, ts_x_yyz, ts_x_yyzz, ts_x_yzz, ts_x_yzzz, ts_x_zzz, ts_x_zzzz, ts_xx_xxxx, ts_xx_xxxy, ts_xx_xxxz, ts_xx_xxyy, ts_xx_xxyz, ts_xx_xxzz, ts_xx_xyyy, ts_xx_xyyz, ts_xx_xyzz, ts_xx_xzzz, ts_xx_yyyy, ts_xx_yyyz, ts_xx_yyzz, ts_xx_yzzz, ts_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xx_xxxx[i] = ts_0_xxxx[i] * gfe_0 + 4.0 * ts_x_xxx[i] * gfe_0 + ts_x_xxxx[i] * ga_x[i];

        ts_xx_xxxy[i] = ts_0_xxxy[i] * gfe_0 + 3.0 * ts_x_xxy[i] * gfe_0 + ts_x_xxxy[i] * ga_x[i];

        ts_xx_xxxz[i] = ts_0_xxxz[i] * gfe_0 + 3.0 * ts_x_xxz[i] * gfe_0 + ts_x_xxxz[i] * ga_x[i];

        ts_xx_xxyy[i] = ts_0_xxyy[i] * gfe_0 + 2.0 * ts_x_xyy[i] * gfe_0 + ts_x_xxyy[i] * ga_x[i];

        ts_xx_xxyz[i] = ts_0_xxyz[i] * gfe_0 + 2.0 * ts_x_xyz[i] * gfe_0 + ts_x_xxyz[i] * ga_x[i];

        ts_xx_xxzz[i] = ts_0_xxzz[i] * gfe_0 + 2.0 * ts_x_xzz[i] * gfe_0 + ts_x_xxzz[i] * ga_x[i];

        ts_xx_xyyy[i] = ts_0_xyyy[i] * gfe_0 + ts_x_yyy[i] * gfe_0 + ts_x_xyyy[i] * ga_x[i];

        ts_xx_xyyz[i] = ts_0_xyyz[i] * gfe_0 + ts_x_yyz[i] * gfe_0 + ts_x_xyyz[i] * ga_x[i];

        ts_xx_xyzz[i] = ts_0_xyzz[i] * gfe_0 + ts_x_yzz[i] * gfe_0 + ts_x_xyzz[i] * ga_x[i];

        ts_xx_xzzz[i] = ts_0_xzzz[i] * gfe_0 + ts_x_zzz[i] * gfe_0 + ts_x_xzzz[i] * ga_x[i];

        ts_xx_yyyy[i] = ts_0_yyyy[i] * gfe_0 + ts_x_yyyy[i] * ga_x[i];

        ts_xx_yyyz[i] = ts_0_yyyz[i] * gfe_0 + ts_x_yyyz[i] * ga_x[i];

        ts_xx_yyzz[i] = ts_0_yyzz[i] * gfe_0 + ts_x_yyzz[i] * ga_x[i];

        ts_xx_yzzz[i] = ts_0_yzzz[i] * gfe_0 + ts_x_yzzz[i] * ga_x[i];

        ts_xx_zzzz[i] = ts_0_zzzz[i] * gfe_0 + ts_x_zzzz[i] * ga_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

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

    #pragma omp simd aligned(ga_x, ga_y, ts_x_xxxx, ts_x_xxxz, ts_x_xxzz, ts_x_xzzz, ts_xy_xxxx, ts_xy_xxxy, ts_xy_xxxz, ts_xy_xxyy, ts_xy_xxyz, ts_xy_xxzz, ts_xy_xyyy, ts_xy_xyyz, ts_xy_xyzz, ts_xy_xzzz, ts_xy_yyyy, ts_xy_yyyz, ts_xy_yyzz, ts_xy_yzzz, ts_xy_zzzz, ts_y_xxxy, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xy_xxxx[i] = ts_x_xxxx[i] * ga_y[i];

        ts_xy_xxxy[i] = 3.0 * ts_y_xxy[i] * gfe_0 + ts_y_xxxy[i] * ga_x[i];

        ts_xy_xxxz[i] = ts_x_xxxz[i] * ga_y[i];

        ts_xy_xxyy[i] = 2.0 * ts_y_xyy[i] * gfe_0 + ts_y_xxyy[i] * ga_x[i];

        ts_xy_xxyz[i] = 2.0 * ts_y_xyz[i] * gfe_0 + ts_y_xxyz[i] * ga_x[i];

        ts_xy_xxzz[i] = ts_x_xxzz[i] * ga_y[i];

        ts_xy_xyyy[i] = ts_y_yyy[i] * gfe_0 + ts_y_xyyy[i] * ga_x[i];

        ts_xy_xyyz[i] = ts_y_yyz[i] * gfe_0 + ts_y_xyyz[i] * ga_x[i];

        ts_xy_xyzz[i] = ts_y_yzz[i] * gfe_0 + ts_y_xyzz[i] * ga_x[i];

        ts_xy_xzzz[i] = ts_x_xzzz[i] * ga_y[i];

        ts_xy_yyyy[i] = ts_y_yyyy[i] * ga_x[i];

        ts_xy_yyyz[i] = ts_y_yyyz[i] * ga_x[i];

        ts_xy_yyzz[i] = ts_y_yyzz[i] * ga_x[i];

        ts_xy_yzzz[i] = ts_y_yzzz[i] * ga_x[i];

        ts_xy_zzzz[i] = ts_y_zzzz[i] * ga_x[i];
    }

    // Set up 30-45 components of targeted buffer : DG

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

    #pragma omp simd aligned(ga_x, ga_z, ts_x_xxxx, ts_x_xxxy, ts_x_xxyy, ts_x_xyyy, ts_xz_xxxx, ts_xz_xxxy, ts_xz_xxxz, ts_xz_xxyy, ts_xz_xxyz, ts_xz_xxzz, ts_xz_xyyy, ts_xz_xyyz, ts_xz_xyzz, ts_xz_xzzz, ts_xz_yyyy, ts_xz_yyyz, ts_xz_yyzz, ts_xz_yzzz, ts_xz_zzzz, ts_z_xxxz, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xz_xxxx[i] = ts_x_xxxx[i] * ga_z[i];

        ts_xz_xxxy[i] = ts_x_xxxy[i] * ga_z[i];

        ts_xz_xxxz[i] = 3.0 * ts_z_xxz[i] * gfe_0 + ts_z_xxxz[i] * ga_x[i];

        ts_xz_xxyy[i] = ts_x_xxyy[i] * ga_z[i];

        ts_xz_xxyz[i] = 2.0 * ts_z_xyz[i] * gfe_0 + ts_z_xxyz[i] * ga_x[i];

        ts_xz_xxzz[i] = 2.0 * ts_z_xzz[i] * gfe_0 + ts_z_xxzz[i] * ga_x[i];

        ts_xz_xyyy[i] = ts_x_xyyy[i] * ga_z[i];

        ts_xz_xyyz[i] = ts_z_yyz[i] * gfe_0 + ts_z_xyyz[i] * ga_x[i];

        ts_xz_xyzz[i] = ts_z_yzz[i] * gfe_0 + ts_z_xyzz[i] * ga_x[i];

        ts_xz_xzzz[i] = ts_z_zzz[i] * gfe_0 + ts_z_xzzz[i] * ga_x[i];

        ts_xz_yyyy[i] = ts_z_yyyy[i] * ga_x[i];

        ts_xz_yyyz[i] = ts_z_yyyz[i] * ga_x[i];

        ts_xz_yyzz[i] = ts_z_yyzz[i] * ga_x[i];

        ts_xz_yzzz[i] = ts_z_yzzz[i] * ga_x[i];

        ts_xz_zzzz[i] = ts_z_zzzz[i] * ga_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

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

    #pragma omp simd aligned(ga_y, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_y_xxx, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxy, ts_y_xxyy, ts_y_xxyz, ts_y_xxz, ts_y_xxzz, ts_y_xyy, ts_y_xyyy, ts_y_xyyz, ts_y_xyz, ts_y_xyzz, ts_y_xzz, ts_y_xzzz, ts_y_yyy, ts_y_yyyy, ts_y_yyyz, ts_y_yyz, ts_y_yyzz, ts_y_yzz, ts_y_yzzz, ts_y_zzz, ts_y_zzzz, ts_yy_xxxx, ts_yy_xxxy, ts_yy_xxxz, ts_yy_xxyy, ts_yy_xxyz, ts_yy_xxzz, ts_yy_xyyy, ts_yy_xyyz, ts_yy_xyzz, ts_yy_xzzz, ts_yy_yyyy, ts_yy_yyyz, ts_yy_yyzz, ts_yy_yzzz, ts_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yy_xxxx[i] = ts_0_xxxx[i] * gfe_0 + ts_y_xxxx[i] * ga_y[i];

        ts_yy_xxxy[i] = ts_0_xxxy[i] * gfe_0 + ts_y_xxx[i] * gfe_0 + ts_y_xxxy[i] * ga_y[i];

        ts_yy_xxxz[i] = ts_0_xxxz[i] * gfe_0 + ts_y_xxxz[i] * ga_y[i];

        ts_yy_xxyy[i] = ts_0_xxyy[i] * gfe_0 + 2.0 * ts_y_xxy[i] * gfe_0 + ts_y_xxyy[i] * ga_y[i];

        ts_yy_xxyz[i] = ts_0_xxyz[i] * gfe_0 + ts_y_xxz[i] * gfe_0 + ts_y_xxyz[i] * ga_y[i];

        ts_yy_xxzz[i] = ts_0_xxzz[i] * gfe_0 + ts_y_xxzz[i] * ga_y[i];

        ts_yy_xyyy[i] = ts_0_xyyy[i] * gfe_0 + 3.0 * ts_y_xyy[i] * gfe_0 + ts_y_xyyy[i] * ga_y[i];

        ts_yy_xyyz[i] = ts_0_xyyz[i] * gfe_0 + 2.0 * ts_y_xyz[i] * gfe_0 + ts_y_xyyz[i] * ga_y[i];

        ts_yy_xyzz[i] = ts_0_xyzz[i] * gfe_0 + ts_y_xzz[i] * gfe_0 + ts_y_xyzz[i] * ga_y[i];

        ts_yy_xzzz[i] = ts_0_xzzz[i] * gfe_0 + ts_y_xzzz[i] * ga_y[i];

        ts_yy_yyyy[i] = ts_0_yyyy[i] * gfe_0 + 4.0 * ts_y_yyy[i] * gfe_0 + ts_y_yyyy[i] * ga_y[i];

        ts_yy_yyyz[i] = ts_0_yyyz[i] * gfe_0 + 3.0 * ts_y_yyz[i] * gfe_0 + ts_y_yyyz[i] * ga_y[i];

        ts_yy_yyzz[i] = ts_0_yyzz[i] * gfe_0 + 2.0 * ts_y_yzz[i] * gfe_0 + ts_y_yyzz[i] * ga_y[i];

        ts_yy_yzzz[i] = ts_0_yzzz[i] * gfe_0 + ts_y_zzz[i] * gfe_0 + ts_y_yzzz[i] * ga_y[i];

        ts_yy_zzzz[i] = ts_0_zzzz[i] * gfe_0 + ts_y_zzzz[i] * ga_y[i];
    }

    // Set up 60-75 components of targeted buffer : DG

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

    #pragma omp simd aligned(ga_y, ga_z, ts_y_xxxy, ts_y_xxyy, ts_y_xyyy, ts_y_yyyy, ts_yz_xxxx, ts_yz_xxxy, ts_yz_xxxz, ts_yz_xxyy, ts_yz_xxyz, ts_yz_xxzz, ts_yz_xyyy, ts_yz_xyyz, ts_yz_xyzz, ts_yz_xzzz, ts_yz_yyyy, ts_yz_yyyz, ts_yz_yyzz, ts_yz_yzzz, ts_yz_zzzz, ts_z_xxxx, ts_z_xxxz, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yz_xxxx[i] = ts_z_xxxx[i] * ga_y[i];

        ts_yz_xxxy[i] = ts_y_xxxy[i] * ga_z[i];

        ts_yz_xxxz[i] = ts_z_xxxz[i] * ga_y[i];

        ts_yz_xxyy[i] = ts_y_xxyy[i] * ga_z[i];

        ts_yz_xxyz[i] = ts_z_xxz[i] * gfe_0 + ts_z_xxyz[i] * ga_y[i];

        ts_yz_xxzz[i] = ts_z_xxzz[i] * ga_y[i];

        ts_yz_xyyy[i] = ts_y_xyyy[i] * ga_z[i];

        ts_yz_xyyz[i] = 2.0 * ts_z_xyz[i] * gfe_0 + ts_z_xyyz[i] * ga_y[i];

        ts_yz_xyzz[i] = ts_z_xzz[i] * gfe_0 + ts_z_xyzz[i] * ga_y[i];

        ts_yz_xzzz[i] = ts_z_xzzz[i] * ga_y[i];

        ts_yz_yyyy[i] = ts_y_yyyy[i] * ga_z[i];

        ts_yz_yyyz[i] = 3.0 * ts_z_yyz[i] * gfe_0 + ts_z_yyyz[i] * ga_y[i];

        ts_yz_yyzz[i] = 2.0 * ts_z_yzz[i] * gfe_0 + ts_z_yyzz[i] * ga_y[i];

        ts_yz_yzzz[i] = ts_z_zzz[i] * gfe_0 + ts_z_yzzz[i] * ga_y[i];

        ts_yz_zzzz[i] = ts_z_zzzz[i] * ga_y[i];
    }

    // Set up 75-90 components of targeted buffer : DG

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

    #pragma omp simd aligned(ga_z, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, ts_z_xxx, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxy, ts_z_xxyy, ts_z_xxyz, ts_z_xxz, ts_z_xxzz, ts_z_xyy, ts_z_xyyy, ts_z_xyyz, ts_z_xyz, ts_z_xyzz, ts_z_xzz, ts_z_xzzz, ts_z_yyy, ts_z_yyyy, ts_z_yyyz, ts_z_yyz, ts_z_yyzz, ts_z_yzz, ts_z_yzzz, ts_z_zzz, ts_z_zzzz, ts_zz_xxxx, ts_zz_xxxy, ts_zz_xxxz, ts_zz_xxyy, ts_zz_xxyz, ts_zz_xxzz, ts_zz_xyyy, ts_zz_xyyz, ts_zz_xyzz, ts_zz_xzzz, ts_zz_yyyy, ts_zz_yyyz, ts_zz_yyzz, ts_zz_yzzz, ts_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zz_xxxx[i] = ts_0_xxxx[i] * gfe_0 + ts_z_xxxx[i] * ga_z[i];

        ts_zz_xxxy[i] = ts_0_xxxy[i] * gfe_0 + ts_z_xxxy[i] * ga_z[i];

        ts_zz_xxxz[i] = ts_0_xxxz[i] * gfe_0 + ts_z_xxx[i] * gfe_0 + ts_z_xxxz[i] * ga_z[i];

        ts_zz_xxyy[i] = ts_0_xxyy[i] * gfe_0 + ts_z_xxyy[i] * ga_z[i];

        ts_zz_xxyz[i] = ts_0_xxyz[i] * gfe_0 + ts_z_xxy[i] * gfe_0 + ts_z_xxyz[i] * ga_z[i];

        ts_zz_xxzz[i] = ts_0_xxzz[i] * gfe_0 + 2.0 * ts_z_xxz[i] * gfe_0 + ts_z_xxzz[i] * ga_z[i];

        ts_zz_xyyy[i] = ts_0_xyyy[i] * gfe_0 + ts_z_xyyy[i] * ga_z[i];

        ts_zz_xyyz[i] = ts_0_xyyz[i] * gfe_0 + ts_z_xyy[i] * gfe_0 + ts_z_xyyz[i] * ga_z[i];

        ts_zz_xyzz[i] = ts_0_xyzz[i] * gfe_0 + 2.0 * ts_z_xyz[i] * gfe_0 + ts_z_xyzz[i] * ga_z[i];

        ts_zz_xzzz[i] = ts_0_xzzz[i] * gfe_0 + 3.0 * ts_z_xzz[i] * gfe_0 + ts_z_xzzz[i] * ga_z[i];

        ts_zz_yyyy[i] = ts_0_yyyy[i] * gfe_0 + ts_z_yyyy[i] * ga_z[i];

        ts_zz_yyyz[i] = ts_0_yyyz[i] * gfe_0 + ts_z_yyy[i] * gfe_0 + ts_z_yyyz[i] * ga_z[i];

        ts_zz_yyzz[i] = ts_0_yyzz[i] * gfe_0 + 2.0 * ts_z_yyz[i] * gfe_0 + ts_z_yyzz[i] * ga_z[i];

        ts_zz_yzzz[i] = ts_0_yzzz[i] * gfe_0 + 3.0 * ts_z_yzz[i] * gfe_0 + ts_z_yzzz[i] * ga_z[i];

        ts_zz_zzzz[i] = ts_0_zzzz[i] * gfe_0 + 4.0 * ts_z_zzz[i] * gfe_0 + ts_z_zzzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

