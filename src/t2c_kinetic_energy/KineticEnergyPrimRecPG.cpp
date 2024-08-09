#include "KineticEnergyPrimRecPG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_pg(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_pg,
                            const size_t              idx_kin_sf,
                            const size_t              idx_kin_sg,
                            const size_t              idx_ovl_pg,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SF

    auto tk_0_xxx = pbuffer.data(idx_kin_sf);

    auto tk_0_xxy = pbuffer.data(idx_kin_sf + 1);

    auto tk_0_xxz = pbuffer.data(idx_kin_sf + 2);

    auto tk_0_xyy = pbuffer.data(idx_kin_sf + 3);

    auto tk_0_xyz = pbuffer.data(idx_kin_sf + 4);

    auto tk_0_xzz = pbuffer.data(idx_kin_sf + 5);

    auto tk_0_yyy = pbuffer.data(idx_kin_sf + 6);

    auto tk_0_yyz = pbuffer.data(idx_kin_sf + 7);

    auto tk_0_yzz = pbuffer.data(idx_kin_sf + 8);

    auto tk_0_zzz = pbuffer.data(idx_kin_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto tk_0_xxxx = pbuffer.data(idx_kin_sg);

    auto tk_0_xxxy = pbuffer.data(idx_kin_sg + 1);

    auto tk_0_xxxz = pbuffer.data(idx_kin_sg + 2);

    auto tk_0_xxyy = pbuffer.data(idx_kin_sg + 3);

    auto tk_0_xxyz = pbuffer.data(idx_kin_sg + 4);

    auto tk_0_xxzz = pbuffer.data(idx_kin_sg + 5);

    auto tk_0_xyyy = pbuffer.data(idx_kin_sg + 6);

    auto tk_0_xyyz = pbuffer.data(idx_kin_sg + 7);

    auto tk_0_xyzz = pbuffer.data(idx_kin_sg + 8);

    auto tk_0_xzzz = pbuffer.data(idx_kin_sg + 9);

    auto tk_0_yyyy = pbuffer.data(idx_kin_sg + 10);

    auto tk_0_yyyz = pbuffer.data(idx_kin_sg + 11);

    auto tk_0_yyzz = pbuffer.data(idx_kin_sg + 12);

    auto tk_0_yzzz = pbuffer.data(idx_kin_sg + 13);

    auto tk_0_zzzz = pbuffer.data(idx_kin_sg + 14);

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_ovl_pg);

    auto ts_x_xxxy = pbuffer.data(idx_ovl_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_ovl_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_ovl_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_ovl_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_ovl_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_ovl_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_ovl_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_ovl_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_ovl_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_ovl_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_ovl_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_ovl_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_ovl_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_ovl_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_ovl_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_ovl_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_ovl_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_ovl_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_ovl_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_ovl_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_ovl_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_ovl_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_ovl_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_ovl_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_ovl_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_ovl_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_ovl_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_ovl_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_ovl_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_ovl_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_ovl_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_ovl_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_ovl_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_ovl_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_ovl_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_ovl_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_ovl_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_ovl_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_ovl_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_ovl_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_ovl_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_ovl_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_ovl_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_ovl_pg + 44);

    // Set up 0-15 components of targeted buffer : PG

    auto tk_x_xxxx = pbuffer.data(idx_kin_pg);

    auto tk_x_xxxy = pbuffer.data(idx_kin_pg + 1);

    auto tk_x_xxxz = pbuffer.data(idx_kin_pg + 2);

    auto tk_x_xxyy = pbuffer.data(idx_kin_pg + 3);

    auto tk_x_xxyz = pbuffer.data(idx_kin_pg + 4);

    auto tk_x_xxzz = pbuffer.data(idx_kin_pg + 5);

    auto tk_x_xyyy = pbuffer.data(idx_kin_pg + 6);

    auto tk_x_xyyz = pbuffer.data(idx_kin_pg + 7);

    auto tk_x_xyzz = pbuffer.data(idx_kin_pg + 8);

    auto tk_x_xzzz = pbuffer.data(idx_kin_pg + 9);

    auto tk_x_yyyy = pbuffer.data(idx_kin_pg + 10);

    auto tk_x_yyyz = pbuffer.data(idx_kin_pg + 11);

    auto tk_x_yyzz = pbuffer.data(idx_kin_pg + 12);

    auto tk_x_yzzz = pbuffer.data(idx_kin_pg + 13);

    auto tk_x_zzzz = pbuffer.data(idx_kin_pg + 14);

#pragma omp simd aligned(pa_x,          \
                             tk_0_xxx,  \
                             tk_0_xxxx, \
                             tk_0_xxxy, \
                             tk_0_xxxz, \
                             tk_0_xxy,  \
                             tk_0_xxyy, \
                             tk_0_xxyz, \
                             tk_0_xxz,  \
                             tk_0_xxzz, \
                             tk_0_xyy,  \
                             tk_0_xyyy, \
                             tk_0_xyyz, \
                             tk_0_xyz,  \
                             tk_0_xyzz, \
                             tk_0_xzz,  \
                             tk_0_xzzz, \
                             tk_0_yyy,  \
                             tk_0_yyyy, \
                             tk_0_yyyz, \
                             tk_0_yyz,  \
                             tk_0_yyzz, \
                             tk_0_yzz,  \
                             tk_0_yzzz, \
                             tk_0_zzz,  \
                             tk_0_zzzz, \
                             tk_x_xxxx, \
                             tk_x_xxxy, \
                             tk_x_xxxz, \
                             tk_x_xxyy, \
                             tk_x_xxyz, \
                             tk_x_xxzz, \
                             tk_x_xyyy, \
                             tk_x_xyyz, \
                             tk_x_xyzz, \
                             tk_x_xzzz, \
                             tk_x_yyyy, \
                             tk_x_yyyz, \
                             tk_x_yyzz, \
                             tk_x_yzzz, \
                             tk_x_zzzz, \
                             ts_x_xxxx, \
                             ts_x_xxxy, \
                             ts_x_xxxz, \
                             ts_x_xxyy, \
                             ts_x_xxyz, \
                             ts_x_xxzz, \
                             ts_x_xyyy, \
                             ts_x_xyyz, \
                             ts_x_xyzz, \
                             ts_x_xzzz, \
                             ts_x_yyyy, \
                             ts_x_yyyz, \
                             ts_x_yyzz, \
                             ts_x_yzzz, \
                             ts_x_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_x_xxxx[i] = 4.0 * tk_0_xxx[i] * fe_0 + tk_0_xxxx[i] * pa_x[i] + 2.0 * ts_x_xxxx[i] * fz_0;

        tk_x_xxxy[i] = 3.0 * tk_0_xxy[i] * fe_0 + tk_0_xxxy[i] * pa_x[i] + 2.0 * ts_x_xxxy[i] * fz_0;

        tk_x_xxxz[i] = 3.0 * tk_0_xxz[i] * fe_0 + tk_0_xxxz[i] * pa_x[i] + 2.0 * ts_x_xxxz[i] * fz_0;

        tk_x_xxyy[i] = 2.0 * tk_0_xyy[i] * fe_0 + tk_0_xxyy[i] * pa_x[i] + 2.0 * ts_x_xxyy[i] * fz_0;

        tk_x_xxyz[i] = 2.0 * tk_0_xyz[i] * fe_0 + tk_0_xxyz[i] * pa_x[i] + 2.0 * ts_x_xxyz[i] * fz_0;

        tk_x_xxzz[i] = 2.0 * tk_0_xzz[i] * fe_0 + tk_0_xxzz[i] * pa_x[i] + 2.0 * ts_x_xxzz[i] * fz_0;

        tk_x_xyyy[i] = tk_0_yyy[i] * fe_0 + tk_0_xyyy[i] * pa_x[i] + 2.0 * ts_x_xyyy[i] * fz_0;

        tk_x_xyyz[i] = tk_0_yyz[i] * fe_0 + tk_0_xyyz[i] * pa_x[i] + 2.0 * ts_x_xyyz[i] * fz_0;

        tk_x_xyzz[i] = tk_0_yzz[i] * fe_0 + tk_0_xyzz[i] * pa_x[i] + 2.0 * ts_x_xyzz[i] * fz_0;

        tk_x_xzzz[i] = tk_0_zzz[i] * fe_0 + tk_0_xzzz[i] * pa_x[i] + 2.0 * ts_x_xzzz[i] * fz_0;

        tk_x_yyyy[i] = tk_0_yyyy[i] * pa_x[i] + 2.0 * ts_x_yyyy[i] * fz_0;

        tk_x_yyyz[i] = tk_0_yyyz[i] * pa_x[i] + 2.0 * ts_x_yyyz[i] * fz_0;

        tk_x_yyzz[i] = tk_0_yyzz[i] * pa_x[i] + 2.0 * ts_x_yyzz[i] * fz_0;

        tk_x_yzzz[i] = tk_0_yzzz[i] * pa_x[i] + 2.0 * ts_x_yzzz[i] * fz_0;

        tk_x_zzzz[i] = tk_0_zzzz[i] * pa_x[i] + 2.0 * ts_x_zzzz[i] * fz_0;
    }

    // Set up 15-30 components of targeted buffer : PG

    auto tk_y_xxxx = pbuffer.data(idx_kin_pg + 15);

    auto tk_y_xxxy = pbuffer.data(idx_kin_pg + 16);

    auto tk_y_xxxz = pbuffer.data(idx_kin_pg + 17);

    auto tk_y_xxyy = pbuffer.data(idx_kin_pg + 18);

    auto tk_y_xxyz = pbuffer.data(idx_kin_pg + 19);

    auto tk_y_xxzz = pbuffer.data(idx_kin_pg + 20);

    auto tk_y_xyyy = pbuffer.data(idx_kin_pg + 21);

    auto tk_y_xyyz = pbuffer.data(idx_kin_pg + 22);

    auto tk_y_xyzz = pbuffer.data(idx_kin_pg + 23);

    auto tk_y_xzzz = pbuffer.data(idx_kin_pg + 24);

    auto tk_y_yyyy = pbuffer.data(idx_kin_pg + 25);

    auto tk_y_yyyz = pbuffer.data(idx_kin_pg + 26);

    auto tk_y_yyzz = pbuffer.data(idx_kin_pg + 27);

    auto tk_y_yzzz = pbuffer.data(idx_kin_pg + 28);

    auto tk_y_zzzz = pbuffer.data(idx_kin_pg + 29);

#pragma omp simd aligned(pa_y,          \
                             tk_0_xxx,  \
                             tk_0_xxxx, \
                             tk_0_xxxy, \
                             tk_0_xxxz, \
                             tk_0_xxy,  \
                             tk_0_xxyy, \
                             tk_0_xxyz, \
                             tk_0_xxz,  \
                             tk_0_xxzz, \
                             tk_0_xyy,  \
                             tk_0_xyyy, \
                             tk_0_xyyz, \
                             tk_0_xyz,  \
                             tk_0_xyzz, \
                             tk_0_xzz,  \
                             tk_0_xzzz, \
                             tk_0_yyy,  \
                             tk_0_yyyy, \
                             tk_0_yyyz, \
                             tk_0_yyz,  \
                             tk_0_yyzz, \
                             tk_0_yzz,  \
                             tk_0_yzzz, \
                             tk_0_zzz,  \
                             tk_0_zzzz, \
                             tk_y_xxxx, \
                             tk_y_xxxy, \
                             tk_y_xxxz, \
                             tk_y_xxyy, \
                             tk_y_xxyz, \
                             tk_y_xxzz, \
                             tk_y_xyyy, \
                             tk_y_xyyz, \
                             tk_y_xyzz, \
                             tk_y_xzzz, \
                             tk_y_yyyy, \
                             tk_y_yyyz, \
                             tk_y_yyzz, \
                             tk_y_yzzz, \
                             tk_y_zzzz, \
                             ts_y_xxxx, \
                             ts_y_xxxy, \
                             ts_y_xxxz, \
                             ts_y_xxyy, \
                             ts_y_xxyz, \
                             ts_y_xxzz, \
                             ts_y_xyyy, \
                             ts_y_xyyz, \
                             ts_y_xyzz, \
                             ts_y_xzzz, \
                             ts_y_yyyy, \
                             ts_y_yyyz, \
                             ts_y_yyzz, \
                             ts_y_yzzz, \
                             ts_y_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_y_xxxx[i] = tk_0_xxxx[i] * pa_y[i] + 2.0 * ts_y_xxxx[i] * fz_0;

        tk_y_xxxy[i] = tk_0_xxx[i] * fe_0 + tk_0_xxxy[i] * pa_y[i] + 2.0 * ts_y_xxxy[i] * fz_0;

        tk_y_xxxz[i] = tk_0_xxxz[i] * pa_y[i] + 2.0 * ts_y_xxxz[i] * fz_0;

        tk_y_xxyy[i] = 2.0 * tk_0_xxy[i] * fe_0 + tk_0_xxyy[i] * pa_y[i] + 2.0 * ts_y_xxyy[i] * fz_0;

        tk_y_xxyz[i] = tk_0_xxz[i] * fe_0 + tk_0_xxyz[i] * pa_y[i] + 2.0 * ts_y_xxyz[i] * fz_0;

        tk_y_xxzz[i] = tk_0_xxzz[i] * pa_y[i] + 2.0 * ts_y_xxzz[i] * fz_0;

        tk_y_xyyy[i] = 3.0 * tk_0_xyy[i] * fe_0 + tk_0_xyyy[i] * pa_y[i] + 2.0 * ts_y_xyyy[i] * fz_0;

        tk_y_xyyz[i] = 2.0 * tk_0_xyz[i] * fe_0 + tk_0_xyyz[i] * pa_y[i] + 2.0 * ts_y_xyyz[i] * fz_0;

        tk_y_xyzz[i] = tk_0_xzz[i] * fe_0 + tk_0_xyzz[i] * pa_y[i] + 2.0 * ts_y_xyzz[i] * fz_0;

        tk_y_xzzz[i] = tk_0_xzzz[i] * pa_y[i] + 2.0 * ts_y_xzzz[i] * fz_0;

        tk_y_yyyy[i] = 4.0 * tk_0_yyy[i] * fe_0 + tk_0_yyyy[i] * pa_y[i] + 2.0 * ts_y_yyyy[i] * fz_0;

        tk_y_yyyz[i] = 3.0 * tk_0_yyz[i] * fe_0 + tk_0_yyyz[i] * pa_y[i] + 2.0 * ts_y_yyyz[i] * fz_0;

        tk_y_yyzz[i] = 2.0 * tk_0_yzz[i] * fe_0 + tk_0_yyzz[i] * pa_y[i] + 2.0 * ts_y_yyzz[i] * fz_0;

        tk_y_yzzz[i] = tk_0_zzz[i] * fe_0 + tk_0_yzzz[i] * pa_y[i] + 2.0 * ts_y_yzzz[i] * fz_0;

        tk_y_zzzz[i] = tk_0_zzzz[i] * pa_y[i] + 2.0 * ts_y_zzzz[i] * fz_0;
    }

    // Set up 30-45 components of targeted buffer : PG

    auto tk_z_xxxx = pbuffer.data(idx_kin_pg + 30);

    auto tk_z_xxxy = pbuffer.data(idx_kin_pg + 31);

    auto tk_z_xxxz = pbuffer.data(idx_kin_pg + 32);

    auto tk_z_xxyy = pbuffer.data(idx_kin_pg + 33);

    auto tk_z_xxyz = pbuffer.data(idx_kin_pg + 34);

    auto tk_z_xxzz = pbuffer.data(idx_kin_pg + 35);

    auto tk_z_xyyy = pbuffer.data(idx_kin_pg + 36);

    auto tk_z_xyyz = pbuffer.data(idx_kin_pg + 37);

    auto tk_z_xyzz = pbuffer.data(idx_kin_pg + 38);

    auto tk_z_xzzz = pbuffer.data(idx_kin_pg + 39);

    auto tk_z_yyyy = pbuffer.data(idx_kin_pg + 40);

    auto tk_z_yyyz = pbuffer.data(idx_kin_pg + 41);

    auto tk_z_yyzz = pbuffer.data(idx_kin_pg + 42);

    auto tk_z_yzzz = pbuffer.data(idx_kin_pg + 43);

    auto tk_z_zzzz = pbuffer.data(idx_kin_pg + 44);

#pragma omp simd aligned(pa_z,          \
                             tk_0_xxx,  \
                             tk_0_xxxx, \
                             tk_0_xxxy, \
                             tk_0_xxxz, \
                             tk_0_xxy,  \
                             tk_0_xxyy, \
                             tk_0_xxyz, \
                             tk_0_xxz,  \
                             tk_0_xxzz, \
                             tk_0_xyy,  \
                             tk_0_xyyy, \
                             tk_0_xyyz, \
                             tk_0_xyz,  \
                             tk_0_xyzz, \
                             tk_0_xzz,  \
                             tk_0_xzzz, \
                             tk_0_yyy,  \
                             tk_0_yyyy, \
                             tk_0_yyyz, \
                             tk_0_yyz,  \
                             tk_0_yyzz, \
                             tk_0_yzz,  \
                             tk_0_yzzz, \
                             tk_0_zzz,  \
                             tk_0_zzzz, \
                             tk_z_xxxx, \
                             tk_z_xxxy, \
                             tk_z_xxxz, \
                             tk_z_xxyy, \
                             tk_z_xxyz, \
                             tk_z_xxzz, \
                             tk_z_xyyy, \
                             tk_z_xyyz, \
                             tk_z_xyzz, \
                             tk_z_xzzz, \
                             tk_z_yyyy, \
                             tk_z_yyyz, \
                             tk_z_yyzz, \
                             tk_z_yzzz, \
                             tk_z_zzzz, \
                             ts_z_xxxx, \
                             ts_z_xxxy, \
                             ts_z_xxxz, \
                             ts_z_xxyy, \
                             ts_z_xxyz, \
                             ts_z_xxzz, \
                             ts_z_xyyy, \
                             ts_z_xyyz, \
                             ts_z_xyzz, \
                             ts_z_xzzz, \
                             ts_z_yyyy, \
                             ts_z_yyyz, \
                             ts_z_yyzz, \
                             ts_z_yzzz, \
                             ts_z_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_z_xxxx[i] = tk_0_xxxx[i] * pa_z[i] + 2.0 * ts_z_xxxx[i] * fz_0;

        tk_z_xxxy[i] = tk_0_xxxy[i] * pa_z[i] + 2.0 * ts_z_xxxy[i] * fz_0;

        tk_z_xxxz[i] = tk_0_xxx[i] * fe_0 + tk_0_xxxz[i] * pa_z[i] + 2.0 * ts_z_xxxz[i] * fz_0;

        tk_z_xxyy[i] = tk_0_xxyy[i] * pa_z[i] + 2.0 * ts_z_xxyy[i] * fz_0;

        tk_z_xxyz[i] = tk_0_xxy[i] * fe_0 + tk_0_xxyz[i] * pa_z[i] + 2.0 * ts_z_xxyz[i] * fz_0;

        tk_z_xxzz[i] = 2.0 * tk_0_xxz[i] * fe_0 + tk_0_xxzz[i] * pa_z[i] + 2.0 * ts_z_xxzz[i] * fz_0;

        tk_z_xyyy[i] = tk_0_xyyy[i] * pa_z[i] + 2.0 * ts_z_xyyy[i] * fz_0;

        tk_z_xyyz[i] = tk_0_xyy[i] * fe_0 + tk_0_xyyz[i] * pa_z[i] + 2.0 * ts_z_xyyz[i] * fz_0;

        tk_z_xyzz[i] = 2.0 * tk_0_xyz[i] * fe_0 + tk_0_xyzz[i] * pa_z[i] + 2.0 * ts_z_xyzz[i] * fz_0;

        tk_z_xzzz[i] = 3.0 * tk_0_xzz[i] * fe_0 + tk_0_xzzz[i] * pa_z[i] + 2.0 * ts_z_xzzz[i] * fz_0;

        tk_z_yyyy[i] = tk_0_yyyy[i] * pa_z[i] + 2.0 * ts_z_yyyy[i] * fz_0;

        tk_z_yyyz[i] = tk_0_yyy[i] * fe_0 + tk_0_yyyz[i] * pa_z[i] + 2.0 * ts_z_yyyz[i] * fz_0;

        tk_z_yyzz[i] = 2.0 * tk_0_yyz[i] * fe_0 + tk_0_yyzz[i] * pa_z[i] + 2.0 * ts_z_yyzz[i] * fz_0;

        tk_z_yzzz[i] = 3.0 * tk_0_yzz[i] * fe_0 + tk_0_yzzz[i] * pa_z[i] + 2.0 * ts_z_yzzz[i] * fz_0;

        tk_z_zzzz[i] = 4.0 * tk_0_zzz[i] * fe_0 + tk_0_zzzz[i] * pa_z[i] + 2.0 * ts_z_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
