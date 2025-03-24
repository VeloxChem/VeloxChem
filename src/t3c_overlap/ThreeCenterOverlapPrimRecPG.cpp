#include "ThreeCenterOverlapPrimRecPG.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_pg(CSimdArray<double>& pbuffer, 
                     const size_t idx_pg,
                     const size_t idx_sf,
                     const size_t idx_sg,
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

    // Set up 0-15 components of targeted buffer : PG

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

    #pragma omp simd aligned(ga_x, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, ts_x_xxxx, ts_x_xxxy, ts_x_xxxz, ts_x_xxyy, ts_x_xxyz, ts_x_xxzz, ts_x_xyyy, ts_x_xyyz, ts_x_xyzz, ts_x_xzzz, ts_x_yyyy, ts_x_yyyz, ts_x_yyzz, ts_x_yzzz, ts_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_x_xxxx[i] = 4.0 * ts_0_xxx[i] * gfe_0 + ts_0_xxxx[i] * ga_x[i];

        ts_x_xxxy[i] = 3.0 * ts_0_xxy[i] * gfe_0 + ts_0_xxxy[i] * ga_x[i];

        ts_x_xxxz[i] = 3.0 * ts_0_xxz[i] * gfe_0 + ts_0_xxxz[i] * ga_x[i];

        ts_x_xxyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 + ts_0_xxyy[i] * ga_x[i];

        ts_x_xxyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 + ts_0_xxyz[i] * ga_x[i];

        ts_x_xxzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 + ts_0_xxzz[i] * ga_x[i];

        ts_x_xyyy[i] = ts_0_yyy[i] * gfe_0 + ts_0_xyyy[i] * ga_x[i];

        ts_x_xyyz[i] = ts_0_yyz[i] * gfe_0 + ts_0_xyyz[i] * ga_x[i];

        ts_x_xyzz[i] = ts_0_yzz[i] * gfe_0 + ts_0_xyzz[i] * ga_x[i];

        ts_x_xzzz[i] = ts_0_zzz[i] * gfe_0 + ts_0_xzzz[i] * ga_x[i];

        ts_x_yyyy[i] = ts_0_yyyy[i] * ga_x[i];

        ts_x_yyyz[i] = ts_0_yyyz[i] * ga_x[i];

        ts_x_yyzz[i] = ts_0_yyzz[i] * ga_x[i];

        ts_x_yzzz[i] = ts_0_yzzz[i] * ga_x[i];

        ts_x_zzzz[i] = ts_0_zzzz[i] * ga_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

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

    #pragma omp simd aligned(ga_y, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, ts_y_xxxx, ts_y_xxxy, ts_y_xxxz, ts_y_xxyy, ts_y_xxyz, ts_y_xxzz, ts_y_xyyy, ts_y_xyyz, ts_y_xyzz, ts_y_xzzz, ts_y_yyyy, ts_y_yyyz, ts_y_yyzz, ts_y_yzzz, ts_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_y_xxxx[i] = ts_0_xxxx[i] * ga_y[i];

        ts_y_xxxy[i] = ts_0_xxx[i] * gfe_0 + ts_0_xxxy[i] * ga_y[i];

        ts_y_xxxz[i] = ts_0_xxxz[i] * ga_y[i];

        ts_y_xxyy[i] = 2.0 * ts_0_xxy[i] * gfe_0 + ts_0_xxyy[i] * ga_y[i];

        ts_y_xxyz[i] = ts_0_xxz[i] * gfe_0 + ts_0_xxyz[i] * ga_y[i];

        ts_y_xxzz[i] = ts_0_xxzz[i] * ga_y[i];

        ts_y_xyyy[i] = 3.0 * ts_0_xyy[i] * gfe_0 + ts_0_xyyy[i] * ga_y[i];

        ts_y_xyyz[i] = 2.0 * ts_0_xyz[i] * gfe_0 + ts_0_xyyz[i] * ga_y[i];

        ts_y_xyzz[i] = ts_0_xzz[i] * gfe_0 + ts_0_xyzz[i] * ga_y[i];

        ts_y_xzzz[i] = ts_0_xzzz[i] * ga_y[i];

        ts_y_yyyy[i] = 4.0 * ts_0_yyy[i] * gfe_0 + ts_0_yyyy[i] * ga_y[i];

        ts_y_yyyz[i] = 3.0 * ts_0_yyz[i] * gfe_0 + ts_0_yyyz[i] * ga_y[i];

        ts_y_yyzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 + ts_0_yyzz[i] * ga_y[i];

        ts_y_yzzz[i] = ts_0_zzz[i] * gfe_0 + ts_0_yzzz[i] * ga_y[i];

        ts_y_zzzz[i] = ts_0_zzzz[i] * ga_y[i];
    }

    // Set up 30-45 components of targeted buffer : PG

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

    #pragma omp simd aligned(ga_z, ts_0_xxx, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxy, ts_0_xxyy, ts_0_xxyz, ts_0_xxz, ts_0_xxzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyz, ts_0_xyz, ts_0_xyzz, ts_0_xzz, ts_0_xzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyz, ts_0_yyz, ts_0_yyzz, ts_0_yzz, ts_0_yzzz, ts_0_zzz, ts_0_zzzz, ts_z_xxxx, ts_z_xxxy, ts_z_xxxz, ts_z_xxyy, ts_z_xxyz, ts_z_xxzz, ts_z_xyyy, ts_z_xyyz, ts_z_xyzz, ts_z_xzzz, ts_z_yyyy, ts_z_yyyz, ts_z_yyzz, ts_z_yzzz, ts_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_z_xxxx[i] = ts_0_xxxx[i] * ga_z[i];

        ts_z_xxxy[i] = ts_0_xxxy[i] * ga_z[i];

        ts_z_xxxz[i] = ts_0_xxx[i] * gfe_0 + ts_0_xxxz[i] * ga_z[i];

        ts_z_xxyy[i] = ts_0_xxyy[i] * ga_z[i];

        ts_z_xxyz[i] = ts_0_xxy[i] * gfe_0 + ts_0_xxyz[i] * ga_z[i];

        ts_z_xxzz[i] = 2.0 * ts_0_xxz[i] * gfe_0 + ts_0_xxzz[i] * ga_z[i];

        ts_z_xyyy[i] = ts_0_xyyy[i] * ga_z[i];

        ts_z_xyyz[i] = ts_0_xyy[i] * gfe_0 + ts_0_xyyz[i] * ga_z[i];

        ts_z_xyzz[i] = 2.0 * ts_0_xyz[i] * gfe_0 + ts_0_xyzz[i] * ga_z[i];

        ts_z_xzzz[i] = 3.0 * ts_0_xzz[i] * gfe_0 + ts_0_xzzz[i] * ga_z[i];

        ts_z_yyyy[i] = ts_0_yyyy[i] * ga_z[i];

        ts_z_yyyz[i] = ts_0_yyy[i] * gfe_0 + ts_0_yyyz[i] * ga_z[i];

        ts_z_yyzz[i] = 2.0 * ts_0_yyz[i] * gfe_0 + ts_0_yyzz[i] * ga_z[i];

        ts_z_yzzz[i] = 3.0 * ts_0_yzz[i] * gfe_0 + ts_0_yzzz[i] * ga_z[i];

        ts_z_zzzz[i] = 4.0 * ts_0_zzz[i] * gfe_0 + ts_0_zzzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

