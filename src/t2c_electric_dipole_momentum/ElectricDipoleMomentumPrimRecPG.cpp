#include "ElectricDipoleMomentumPrimRecPG.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_pg(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_pg,
                                      const size_t idx_dip_sf,
                                      const size_t idx_ovl_sg,
                                      const size_t idx_dip_sg,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SF

    auto tr_x_0_xxx = pbuffer.data(idx_dip_sf);

    auto tr_x_0_xxy = pbuffer.data(idx_dip_sf + 1);

    auto tr_x_0_xxz = pbuffer.data(idx_dip_sf + 2);

    auto tr_x_0_xyy = pbuffer.data(idx_dip_sf + 3);

    auto tr_x_0_xyz = pbuffer.data(idx_dip_sf + 4);

    auto tr_x_0_xzz = pbuffer.data(idx_dip_sf + 5);

    auto tr_x_0_yyy = pbuffer.data(idx_dip_sf + 6);

    auto tr_x_0_yyz = pbuffer.data(idx_dip_sf + 7);

    auto tr_x_0_yzz = pbuffer.data(idx_dip_sf + 8);

    auto tr_x_0_zzz = pbuffer.data(idx_dip_sf + 9);

    auto tr_y_0_xxx = pbuffer.data(idx_dip_sf + 10);

    auto tr_y_0_xxy = pbuffer.data(idx_dip_sf + 11);

    auto tr_y_0_xxz = pbuffer.data(idx_dip_sf + 12);

    auto tr_y_0_xyy = pbuffer.data(idx_dip_sf + 13);

    auto tr_y_0_xyz = pbuffer.data(idx_dip_sf + 14);

    auto tr_y_0_xzz = pbuffer.data(idx_dip_sf + 15);

    auto tr_y_0_yyy = pbuffer.data(idx_dip_sf + 16);

    auto tr_y_0_yyz = pbuffer.data(idx_dip_sf + 17);

    auto tr_y_0_yzz = pbuffer.data(idx_dip_sf + 18);

    auto tr_y_0_zzz = pbuffer.data(idx_dip_sf + 19);

    auto tr_z_0_xxx = pbuffer.data(idx_dip_sf + 20);

    auto tr_z_0_xxy = pbuffer.data(idx_dip_sf + 21);

    auto tr_z_0_xxz = pbuffer.data(idx_dip_sf + 22);

    auto tr_z_0_xyy = pbuffer.data(idx_dip_sf + 23);

    auto tr_z_0_xyz = pbuffer.data(idx_dip_sf + 24);

    auto tr_z_0_xzz = pbuffer.data(idx_dip_sf + 25);

    auto tr_z_0_yyy = pbuffer.data(idx_dip_sf + 26);

    auto tr_z_0_yyz = pbuffer.data(idx_dip_sf + 27);

    auto tr_z_0_yzz = pbuffer.data(idx_dip_sf + 28);

    auto tr_z_0_zzz = pbuffer.data(idx_dip_sf + 29);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxxy = pbuffer.data(idx_ovl_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_ovl_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_ovl_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_ovl_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_ovl_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_ovl_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto tr_x_0_xxxx = pbuffer.data(idx_dip_sg);

    auto tr_x_0_xxxy = pbuffer.data(idx_dip_sg + 1);

    auto tr_x_0_xxxz = pbuffer.data(idx_dip_sg + 2);

    auto tr_x_0_xxyy = pbuffer.data(idx_dip_sg + 3);

    auto tr_x_0_xxyz = pbuffer.data(idx_dip_sg + 4);

    auto tr_x_0_xxzz = pbuffer.data(idx_dip_sg + 5);

    auto tr_x_0_xyyy = pbuffer.data(idx_dip_sg + 6);

    auto tr_x_0_xyyz = pbuffer.data(idx_dip_sg + 7);

    auto tr_x_0_xyzz = pbuffer.data(idx_dip_sg + 8);

    auto tr_x_0_xzzz = pbuffer.data(idx_dip_sg + 9);

    auto tr_x_0_yyyy = pbuffer.data(idx_dip_sg + 10);

    auto tr_x_0_yyyz = pbuffer.data(idx_dip_sg + 11);

    auto tr_x_0_yyzz = pbuffer.data(idx_dip_sg + 12);

    auto tr_x_0_yzzz = pbuffer.data(idx_dip_sg + 13);

    auto tr_x_0_zzzz = pbuffer.data(idx_dip_sg + 14);

    auto tr_y_0_xxxx = pbuffer.data(idx_dip_sg + 15);

    auto tr_y_0_xxxy = pbuffer.data(idx_dip_sg + 16);

    auto tr_y_0_xxxz = pbuffer.data(idx_dip_sg + 17);

    auto tr_y_0_xxyy = pbuffer.data(idx_dip_sg + 18);

    auto tr_y_0_xxyz = pbuffer.data(idx_dip_sg + 19);

    auto tr_y_0_xxzz = pbuffer.data(idx_dip_sg + 20);

    auto tr_y_0_xyyy = pbuffer.data(idx_dip_sg + 21);

    auto tr_y_0_xyyz = pbuffer.data(idx_dip_sg + 22);

    auto tr_y_0_xyzz = pbuffer.data(idx_dip_sg + 23);

    auto tr_y_0_xzzz = pbuffer.data(idx_dip_sg + 24);

    auto tr_y_0_yyyy = pbuffer.data(idx_dip_sg + 25);

    auto tr_y_0_yyyz = pbuffer.data(idx_dip_sg + 26);

    auto tr_y_0_yyzz = pbuffer.data(idx_dip_sg + 27);

    auto tr_y_0_yzzz = pbuffer.data(idx_dip_sg + 28);

    auto tr_y_0_zzzz = pbuffer.data(idx_dip_sg + 29);

    auto tr_z_0_xxxx = pbuffer.data(idx_dip_sg + 30);

    auto tr_z_0_xxxy = pbuffer.data(idx_dip_sg + 31);

    auto tr_z_0_xxxz = pbuffer.data(idx_dip_sg + 32);

    auto tr_z_0_xxyy = pbuffer.data(idx_dip_sg + 33);

    auto tr_z_0_xxyz = pbuffer.data(idx_dip_sg + 34);

    auto tr_z_0_xxzz = pbuffer.data(idx_dip_sg + 35);

    auto tr_z_0_xyyy = pbuffer.data(idx_dip_sg + 36);

    auto tr_z_0_xyyz = pbuffer.data(idx_dip_sg + 37);

    auto tr_z_0_xyzz = pbuffer.data(idx_dip_sg + 38);

    auto tr_z_0_xzzz = pbuffer.data(idx_dip_sg + 39);

    auto tr_z_0_yyyy = pbuffer.data(idx_dip_sg + 40);

    auto tr_z_0_yyyz = pbuffer.data(idx_dip_sg + 41);

    auto tr_z_0_yyzz = pbuffer.data(idx_dip_sg + 42);

    auto tr_z_0_yzzz = pbuffer.data(idx_dip_sg + 43);

    auto tr_z_0_zzzz = pbuffer.data(idx_dip_sg + 44);

    // Set up 0-15 components of targeted buffer : PG

    auto tr_x_x_xxxx = pbuffer.data(idx_dip_pg);

    auto tr_x_x_xxxy = pbuffer.data(idx_dip_pg + 1);

    auto tr_x_x_xxxz = pbuffer.data(idx_dip_pg + 2);

    auto tr_x_x_xxyy = pbuffer.data(idx_dip_pg + 3);

    auto tr_x_x_xxyz = pbuffer.data(idx_dip_pg + 4);

    auto tr_x_x_xxzz = pbuffer.data(idx_dip_pg + 5);

    auto tr_x_x_xyyy = pbuffer.data(idx_dip_pg + 6);

    auto tr_x_x_xyyz = pbuffer.data(idx_dip_pg + 7);

    auto tr_x_x_xyzz = pbuffer.data(idx_dip_pg + 8);

    auto tr_x_x_xzzz = pbuffer.data(idx_dip_pg + 9);

    auto tr_x_x_yyyy = pbuffer.data(idx_dip_pg + 10);

    auto tr_x_x_yyyz = pbuffer.data(idx_dip_pg + 11);

    auto tr_x_x_yyzz = pbuffer.data(idx_dip_pg + 12);

    auto tr_x_x_yzzz = pbuffer.data(idx_dip_pg + 13);

    auto tr_x_x_zzzz = pbuffer.data(idx_dip_pg + 14);

    #pragma omp simd aligned(pa_x, tr_x_0_xxx, tr_x_0_xxxx, tr_x_0_xxxy, tr_x_0_xxxz, tr_x_0_xxy, tr_x_0_xxyy, tr_x_0_xxyz, tr_x_0_xxz, tr_x_0_xxzz, tr_x_0_xyy, tr_x_0_xyyy, tr_x_0_xyyz, tr_x_0_xyz, tr_x_0_xyzz, tr_x_0_xzz, tr_x_0_xzzz, tr_x_0_yyy, tr_x_0_yyyy, tr_x_0_yyyz, tr_x_0_yyz, tr_x_0_yyzz, tr_x_0_yzz, tr_x_0_yzzz, tr_x_0_zzz, tr_x_0_zzzz, tr_x_x_xxxx, tr_x_x_xxxy, tr_x_x_xxxz, tr_x_x_xxyy, tr_x_x_xxyz, tr_x_x_xxzz, tr_x_x_xyyy, tr_x_x_xyyz, tr_x_x_xyzz, tr_x_x_xzzz, tr_x_x_yyyy, tr_x_x_yyyz, tr_x_x_yyzz, tr_x_x_yzzz, tr_x_x_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_xxxx[i] = 4.0 * tr_x_0_xxx[i] * fe_0 + ts_0_xxxx[i] * fe_0 + tr_x_0_xxxx[i] * pa_x[i];

        tr_x_x_xxxy[i] = 3.0 * tr_x_0_xxy[i] * fe_0 + ts_0_xxxy[i] * fe_0 + tr_x_0_xxxy[i] * pa_x[i];

        tr_x_x_xxxz[i] = 3.0 * tr_x_0_xxz[i] * fe_0 + ts_0_xxxz[i] * fe_0 + tr_x_0_xxxz[i] * pa_x[i];

        tr_x_x_xxyy[i] = 2.0 * tr_x_0_xyy[i] * fe_0 + ts_0_xxyy[i] * fe_0 + tr_x_0_xxyy[i] * pa_x[i];

        tr_x_x_xxyz[i] = 2.0 * tr_x_0_xyz[i] * fe_0 + ts_0_xxyz[i] * fe_0 + tr_x_0_xxyz[i] * pa_x[i];

        tr_x_x_xxzz[i] = 2.0 * tr_x_0_xzz[i] * fe_0 + ts_0_xxzz[i] * fe_0 + tr_x_0_xxzz[i] * pa_x[i];

        tr_x_x_xyyy[i] = tr_x_0_yyy[i] * fe_0 + ts_0_xyyy[i] * fe_0 + tr_x_0_xyyy[i] * pa_x[i];

        tr_x_x_xyyz[i] = tr_x_0_yyz[i] * fe_0 + ts_0_xyyz[i] * fe_0 + tr_x_0_xyyz[i] * pa_x[i];

        tr_x_x_xyzz[i] = tr_x_0_yzz[i] * fe_0 + ts_0_xyzz[i] * fe_0 + tr_x_0_xyzz[i] * pa_x[i];

        tr_x_x_xzzz[i] = tr_x_0_zzz[i] * fe_0 + ts_0_xzzz[i] * fe_0 + tr_x_0_xzzz[i] * pa_x[i];

        tr_x_x_yyyy[i] = ts_0_yyyy[i] * fe_0 + tr_x_0_yyyy[i] * pa_x[i];

        tr_x_x_yyyz[i] = ts_0_yyyz[i] * fe_0 + tr_x_0_yyyz[i] * pa_x[i];

        tr_x_x_yyzz[i] = ts_0_yyzz[i] * fe_0 + tr_x_0_yyzz[i] * pa_x[i];

        tr_x_x_yzzz[i] = ts_0_yzzz[i] * fe_0 + tr_x_0_yzzz[i] * pa_x[i];

        tr_x_x_zzzz[i] = ts_0_zzzz[i] * fe_0 + tr_x_0_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : PG

    auto tr_x_y_xxxx = pbuffer.data(idx_dip_pg + 15);

    auto tr_x_y_xxxy = pbuffer.data(idx_dip_pg + 16);

    auto tr_x_y_xxxz = pbuffer.data(idx_dip_pg + 17);

    auto tr_x_y_xxyy = pbuffer.data(idx_dip_pg + 18);

    auto tr_x_y_xxyz = pbuffer.data(idx_dip_pg + 19);

    auto tr_x_y_xxzz = pbuffer.data(idx_dip_pg + 20);

    auto tr_x_y_xyyy = pbuffer.data(idx_dip_pg + 21);

    auto tr_x_y_xyyz = pbuffer.data(idx_dip_pg + 22);

    auto tr_x_y_xyzz = pbuffer.data(idx_dip_pg + 23);

    auto tr_x_y_xzzz = pbuffer.data(idx_dip_pg + 24);

    auto tr_x_y_yyyy = pbuffer.data(idx_dip_pg + 25);

    auto tr_x_y_yyyz = pbuffer.data(idx_dip_pg + 26);

    auto tr_x_y_yyzz = pbuffer.data(idx_dip_pg + 27);

    auto tr_x_y_yzzz = pbuffer.data(idx_dip_pg + 28);

    auto tr_x_y_zzzz = pbuffer.data(idx_dip_pg + 29);

    #pragma omp simd aligned(pa_y, tr_x_0_xxx, tr_x_0_xxxx, tr_x_0_xxxy, tr_x_0_xxxz, tr_x_0_xxy, tr_x_0_xxyy, tr_x_0_xxyz, tr_x_0_xxz, tr_x_0_xxzz, tr_x_0_xyy, tr_x_0_xyyy, tr_x_0_xyyz, tr_x_0_xyz, tr_x_0_xyzz, tr_x_0_xzz, tr_x_0_xzzz, tr_x_0_yyy, tr_x_0_yyyy, tr_x_0_yyyz, tr_x_0_yyz, tr_x_0_yyzz, tr_x_0_yzz, tr_x_0_yzzz, tr_x_0_zzz, tr_x_0_zzzz, tr_x_y_xxxx, tr_x_y_xxxy, tr_x_y_xxxz, tr_x_y_xxyy, tr_x_y_xxyz, tr_x_y_xxzz, tr_x_y_xyyy, tr_x_y_xyyz, tr_x_y_xyzz, tr_x_y_xzzz, tr_x_y_yyyy, tr_x_y_yyyz, tr_x_y_yyzz, tr_x_y_yzzz, tr_x_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_y_xxxx[i] = tr_x_0_xxxx[i] * pa_y[i];

        tr_x_y_xxxy[i] = tr_x_0_xxx[i] * fe_0 + tr_x_0_xxxy[i] * pa_y[i];

        tr_x_y_xxxz[i] = tr_x_0_xxxz[i] * pa_y[i];

        tr_x_y_xxyy[i] = 2.0 * tr_x_0_xxy[i] * fe_0 + tr_x_0_xxyy[i] * pa_y[i];

        tr_x_y_xxyz[i] = tr_x_0_xxz[i] * fe_0 + tr_x_0_xxyz[i] * pa_y[i];

        tr_x_y_xxzz[i] = tr_x_0_xxzz[i] * pa_y[i];

        tr_x_y_xyyy[i] = 3.0 * tr_x_0_xyy[i] * fe_0 + tr_x_0_xyyy[i] * pa_y[i];

        tr_x_y_xyyz[i] = 2.0 * tr_x_0_xyz[i] * fe_0 + tr_x_0_xyyz[i] * pa_y[i];

        tr_x_y_xyzz[i] = tr_x_0_xzz[i] * fe_0 + tr_x_0_xyzz[i] * pa_y[i];

        tr_x_y_xzzz[i] = tr_x_0_xzzz[i] * pa_y[i];

        tr_x_y_yyyy[i] = 4.0 * tr_x_0_yyy[i] * fe_0 + tr_x_0_yyyy[i] * pa_y[i];

        tr_x_y_yyyz[i] = 3.0 * tr_x_0_yyz[i] * fe_0 + tr_x_0_yyyz[i] * pa_y[i];

        tr_x_y_yyzz[i] = 2.0 * tr_x_0_yzz[i] * fe_0 + tr_x_0_yyzz[i] * pa_y[i];

        tr_x_y_yzzz[i] = tr_x_0_zzz[i] * fe_0 + tr_x_0_yzzz[i] * pa_y[i];

        tr_x_y_zzzz[i] = tr_x_0_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : PG

    auto tr_x_z_xxxx = pbuffer.data(idx_dip_pg + 30);

    auto tr_x_z_xxxy = pbuffer.data(idx_dip_pg + 31);

    auto tr_x_z_xxxz = pbuffer.data(idx_dip_pg + 32);

    auto tr_x_z_xxyy = pbuffer.data(idx_dip_pg + 33);

    auto tr_x_z_xxyz = pbuffer.data(idx_dip_pg + 34);

    auto tr_x_z_xxzz = pbuffer.data(idx_dip_pg + 35);

    auto tr_x_z_xyyy = pbuffer.data(idx_dip_pg + 36);

    auto tr_x_z_xyyz = pbuffer.data(idx_dip_pg + 37);

    auto tr_x_z_xyzz = pbuffer.data(idx_dip_pg + 38);

    auto tr_x_z_xzzz = pbuffer.data(idx_dip_pg + 39);

    auto tr_x_z_yyyy = pbuffer.data(idx_dip_pg + 40);

    auto tr_x_z_yyyz = pbuffer.data(idx_dip_pg + 41);

    auto tr_x_z_yyzz = pbuffer.data(idx_dip_pg + 42);

    auto tr_x_z_yzzz = pbuffer.data(idx_dip_pg + 43);

    auto tr_x_z_zzzz = pbuffer.data(idx_dip_pg + 44);

    #pragma omp simd aligned(pa_z, tr_x_0_xxx, tr_x_0_xxxx, tr_x_0_xxxy, tr_x_0_xxxz, tr_x_0_xxy, tr_x_0_xxyy, tr_x_0_xxyz, tr_x_0_xxz, tr_x_0_xxzz, tr_x_0_xyy, tr_x_0_xyyy, tr_x_0_xyyz, tr_x_0_xyz, tr_x_0_xyzz, tr_x_0_xzz, tr_x_0_xzzz, tr_x_0_yyy, tr_x_0_yyyy, tr_x_0_yyyz, tr_x_0_yyz, tr_x_0_yyzz, tr_x_0_yzz, tr_x_0_yzzz, tr_x_0_zzz, tr_x_0_zzzz, tr_x_z_xxxx, tr_x_z_xxxy, tr_x_z_xxxz, tr_x_z_xxyy, tr_x_z_xxyz, tr_x_z_xxzz, tr_x_z_xyyy, tr_x_z_xyyz, tr_x_z_xyzz, tr_x_z_xzzz, tr_x_z_yyyy, tr_x_z_yyyz, tr_x_z_yyzz, tr_x_z_yzzz, tr_x_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_z_xxxx[i] = tr_x_0_xxxx[i] * pa_z[i];

        tr_x_z_xxxy[i] = tr_x_0_xxxy[i] * pa_z[i];

        tr_x_z_xxxz[i] = tr_x_0_xxx[i] * fe_0 + tr_x_0_xxxz[i] * pa_z[i];

        tr_x_z_xxyy[i] = tr_x_0_xxyy[i] * pa_z[i];

        tr_x_z_xxyz[i] = tr_x_0_xxy[i] * fe_0 + tr_x_0_xxyz[i] * pa_z[i];

        tr_x_z_xxzz[i] = 2.0 * tr_x_0_xxz[i] * fe_0 + tr_x_0_xxzz[i] * pa_z[i];

        tr_x_z_xyyy[i] = tr_x_0_xyyy[i] * pa_z[i];

        tr_x_z_xyyz[i] = tr_x_0_xyy[i] * fe_0 + tr_x_0_xyyz[i] * pa_z[i];

        tr_x_z_xyzz[i] = 2.0 * tr_x_0_xyz[i] * fe_0 + tr_x_0_xyzz[i] * pa_z[i];

        tr_x_z_xzzz[i] = 3.0 * tr_x_0_xzz[i] * fe_0 + tr_x_0_xzzz[i] * pa_z[i];

        tr_x_z_yyyy[i] = tr_x_0_yyyy[i] * pa_z[i];

        tr_x_z_yyyz[i] = tr_x_0_yyy[i] * fe_0 + tr_x_0_yyyz[i] * pa_z[i];

        tr_x_z_yyzz[i] = 2.0 * tr_x_0_yyz[i] * fe_0 + tr_x_0_yyzz[i] * pa_z[i];

        tr_x_z_yzzz[i] = 3.0 * tr_x_0_yzz[i] * fe_0 + tr_x_0_yzzz[i] * pa_z[i];

        tr_x_z_zzzz[i] = 4.0 * tr_x_0_zzz[i] * fe_0 + tr_x_0_zzzz[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : PG

    auto tr_y_x_xxxx = pbuffer.data(idx_dip_pg + 45);

    auto tr_y_x_xxxy = pbuffer.data(idx_dip_pg + 46);

    auto tr_y_x_xxxz = pbuffer.data(idx_dip_pg + 47);

    auto tr_y_x_xxyy = pbuffer.data(idx_dip_pg + 48);

    auto tr_y_x_xxyz = pbuffer.data(idx_dip_pg + 49);

    auto tr_y_x_xxzz = pbuffer.data(idx_dip_pg + 50);

    auto tr_y_x_xyyy = pbuffer.data(idx_dip_pg + 51);

    auto tr_y_x_xyyz = pbuffer.data(idx_dip_pg + 52);

    auto tr_y_x_xyzz = pbuffer.data(idx_dip_pg + 53);

    auto tr_y_x_xzzz = pbuffer.data(idx_dip_pg + 54);

    auto tr_y_x_yyyy = pbuffer.data(idx_dip_pg + 55);

    auto tr_y_x_yyyz = pbuffer.data(idx_dip_pg + 56);

    auto tr_y_x_yyzz = pbuffer.data(idx_dip_pg + 57);

    auto tr_y_x_yzzz = pbuffer.data(idx_dip_pg + 58);

    auto tr_y_x_zzzz = pbuffer.data(idx_dip_pg + 59);

    #pragma omp simd aligned(pa_x, tr_y_0_xxx, tr_y_0_xxxx, tr_y_0_xxxy, tr_y_0_xxxz, tr_y_0_xxy, tr_y_0_xxyy, tr_y_0_xxyz, tr_y_0_xxz, tr_y_0_xxzz, tr_y_0_xyy, tr_y_0_xyyy, tr_y_0_xyyz, tr_y_0_xyz, tr_y_0_xyzz, tr_y_0_xzz, tr_y_0_xzzz, tr_y_0_yyy, tr_y_0_yyyy, tr_y_0_yyyz, tr_y_0_yyz, tr_y_0_yyzz, tr_y_0_yzz, tr_y_0_yzzz, tr_y_0_zzz, tr_y_0_zzzz, tr_y_x_xxxx, tr_y_x_xxxy, tr_y_x_xxxz, tr_y_x_xxyy, tr_y_x_xxyz, tr_y_x_xxzz, tr_y_x_xyyy, tr_y_x_xyyz, tr_y_x_xyzz, tr_y_x_xzzz, tr_y_x_yyyy, tr_y_x_yyyz, tr_y_x_yyzz, tr_y_x_yzzz, tr_y_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_x_xxxx[i] = 4.0 * tr_y_0_xxx[i] * fe_0 + tr_y_0_xxxx[i] * pa_x[i];

        tr_y_x_xxxy[i] = 3.0 * tr_y_0_xxy[i] * fe_0 + tr_y_0_xxxy[i] * pa_x[i];

        tr_y_x_xxxz[i] = 3.0 * tr_y_0_xxz[i] * fe_0 + tr_y_0_xxxz[i] * pa_x[i];

        tr_y_x_xxyy[i] = 2.0 * tr_y_0_xyy[i] * fe_0 + tr_y_0_xxyy[i] * pa_x[i];

        tr_y_x_xxyz[i] = 2.0 * tr_y_0_xyz[i] * fe_0 + tr_y_0_xxyz[i] * pa_x[i];

        tr_y_x_xxzz[i] = 2.0 * tr_y_0_xzz[i] * fe_0 + tr_y_0_xxzz[i] * pa_x[i];

        tr_y_x_xyyy[i] = tr_y_0_yyy[i] * fe_0 + tr_y_0_xyyy[i] * pa_x[i];

        tr_y_x_xyyz[i] = tr_y_0_yyz[i] * fe_0 + tr_y_0_xyyz[i] * pa_x[i];

        tr_y_x_xyzz[i] = tr_y_0_yzz[i] * fe_0 + tr_y_0_xyzz[i] * pa_x[i];

        tr_y_x_xzzz[i] = tr_y_0_zzz[i] * fe_0 + tr_y_0_xzzz[i] * pa_x[i];

        tr_y_x_yyyy[i] = tr_y_0_yyyy[i] * pa_x[i];

        tr_y_x_yyyz[i] = tr_y_0_yyyz[i] * pa_x[i];

        tr_y_x_yyzz[i] = tr_y_0_yyzz[i] * pa_x[i];

        tr_y_x_yzzz[i] = tr_y_0_yzzz[i] * pa_x[i];

        tr_y_x_zzzz[i] = tr_y_0_zzzz[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : PG

    auto tr_y_y_xxxx = pbuffer.data(idx_dip_pg + 60);

    auto tr_y_y_xxxy = pbuffer.data(idx_dip_pg + 61);

    auto tr_y_y_xxxz = pbuffer.data(idx_dip_pg + 62);

    auto tr_y_y_xxyy = pbuffer.data(idx_dip_pg + 63);

    auto tr_y_y_xxyz = pbuffer.data(idx_dip_pg + 64);

    auto tr_y_y_xxzz = pbuffer.data(idx_dip_pg + 65);

    auto tr_y_y_xyyy = pbuffer.data(idx_dip_pg + 66);

    auto tr_y_y_xyyz = pbuffer.data(idx_dip_pg + 67);

    auto tr_y_y_xyzz = pbuffer.data(idx_dip_pg + 68);

    auto tr_y_y_xzzz = pbuffer.data(idx_dip_pg + 69);

    auto tr_y_y_yyyy = pbuffer.data(idx_dip_pg + 70);

    auto tr_y_y_yyyz = pbuffer.data(idx_dip_pg + 71);

    auto tr_y_y_yyzz = pbuffer.data(idx_dip_pg + 72);

    auto tr_y_y_yzzz = pbuffer.data(idx_dip_pg + 73);

    auto tr_y_y_zzzz = pbuffer.data(idx_dip_pg + 74);

    #pragma omp simd aligned(pa_y, tr_y_0_xxx, tr_y_0_xxxx, tr_y_0_xxxy, tr_y_0_xxxz, tr_y_0_xxy, tr_y_0_xxyy, tr_y_0_xxyz, tr_y_0_xxz, tr_y_0_xxzz, tr_y_0_xyy, tr_y_0_xyyy, tr_y_0_xyyz, tr_y_0_xyz, tr_y_0_xyzz, tr_y_0_xzz, tr_y_0_xzzz, tr_y_0_yyy, tr_y_0_yyyy, tr_y_0_yyyz, tr_y_0_yyz, tr_y_0_yyzz, tr_y_0_yzz, tr_y_0_yzzz, tr_y_0_zzz, tr_y_0_zzzz, tr_y_y_xxxx, tr_y_y_xxxy, tr_y_y_xxxz, tr_y_y_xxyy, tr_y_y_xxyz, tr_y_y_xxzz, tr_y_y_xyyy, tr_y_y_xyyz, tr_y_y_xyzz, tr_y_y_xzzz, tr_y_y_yyyy, tr_y_y_yyyz, tr_y_y_yyzz, tr_y_y_yzzz, tr_y_y_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_y_xxxx[i] = ts_0_xxxx[i] * fe_0 + tr_y_0_xxxx[i] * pa_y[i];

        tr_y_y_xxxy[i] = tr_y_0_xxx[i] * fe_0 + ts_0_xxxy[i] * fe_0 + tr_y_0_xxxy[i] * pa_y[i];

        tr_y_y_xxxz[i] = ts_0_xxxz[i] * fe_0 + tr_y_0_xxxz[i] * pa_y[i];

        tr_y_y_xxyy[i] = 2.0 * tr_y_0_xxy[i] * fe_0 + ts_0_xxyy[i] * fe_0 + tr_y_0_xxyy[i] * pa_y[i];

        tr_y_y_xxyz[i] = tr_y_0_xxz[i] * fe_0 + ts_0_xxyz[i] * fe_0 + tr_y_0_xxyz[i] * pa_y[i];

        tr_y_y_xxzz[i] = ts_0_xxzz[i] * fe_0 + tr_y_0_xxzz[i] * pa_y[i];

        tr_y_y_xyyy[i] = 3.0 * tr_y_0_xyy[i] * fe_0 + ts_0_xyyy[i] * fe_0 + tr_y_0_xyyy[i] * pa_y[i];

        tr_y_y_xyyz[i] = 2.0 * tr_y_0_xyz[i] * fe_0 + ts_0_xyyz[i] * fe_0 + tr_y_0_xyyz[i] * pa_y[i];

        tr_y_y_xyzz[i] = tr_y_0_xzz[i] * fe_0 + ts_0_xyzz[i] * fe_0 + tr_y_0_xyzz[i] * pa_y[i];

        tr_y_y_xzzz[i] = ts_0_xzzz[i] * fe_0 + tr_y_0_xzzz[i] * pa_y[i];

        tr_y_y_yyyy[i] = 4.0 * tr_y_0_yyy[i] * fe_0 + ts_0_yyyy[i] * fe_0 + tr_y_0_yyyy[i] * pa_y[i];

        tr_y_y_yyyz[i] = 3.0 * tr_y_0_yyz[i] * fe_0 + ts_0_yyyz[i] * fe_0 + tr_y_0_yyyz[i] * pa_y[i];

        tr_y_y_yyzz[i] = 2.0 * tr_y_0_yzz[i] * fe_0 + ts_0_yyzz[i] * fe_0 + tr_y_0_yyzz[i] * pa_y[i];

        tr_y_y_yzzz[i] = tr_y_0_zzz[i] * fe_0 + ts_0_yzzz[i] * fe_0 + tr_y_0_yzzz[i] * pa_y[i];

        tr_y_y_zzzz[i] = ts_0_zzzz[i] * fe_0 + tr_y_0_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : PG

    auto tr_y_z_xxxx = pbuffer.data(idx_dip_pg + 75);

    auto tr_y_z_xxxy = pbuffer.data(idx_dip_pg + 76);

    auto tr_y_z_xxxz = pbuffer.data(idx_dip_pg + 77);

    auto tr_y_z_xxyy = pbuffer.data(idx_dip_pg + 78);

    auto tr_y_z_xxyz = pbuffer.data(idx_dip_pg + 79);

    auto tr_y_z_xxzz = pbuffer.data(idx_dip_pg + 80);

    auto tr_y_z_xyyy = pbuffer.data(idx_dip_pg + 81);

    auto tr_y_z_xyyz = pbuffer.data(idx_dip_pg + 82);

    auto tr_y_z_xyzz = pbuffer.data(idx_dip_pg + 83);

    auto tr_y_z_xzzz = pbuffer.data(idx_dip_pg + 84);

    auto tr_y_z_yyyy = pbuffer.data(idx_dip_pg + 85);

    auto tr_y_z_yyyz = pbuffer.data(idx_dip_pg + 86);

    auto tr_y_z_yyzz = pbuffer.data(idx_dip_pg + 87);

    auto tr_y_z_yzzz = pbuffer.data(idx_dip_pg + 88);

    auto tr_y_z_zzzz = pbuffer.data(idx_dip_pg + 89);

    #pragma omp simd aligned(pa_z, tr_y_0_xxx, tr_y_0_xxxx, tr_y_0_xxxy, tr_y_0_xxxz, tr_y_0_xxy, tr_y_0_xxyy, tr_y_0_xxyz, tr_y_0_xxz, tr_y_0_xxzz, tr_y_0_xyy, tr_y_0_xyyy, tr_y_0_xyyz, tr_y_0_xyz, tr_y_0_xyzz, tr_y_0_xzz, tr_y_0_xzzz, tr_y_0_yyy, tr_y_0_yyyy, tr_y_0_yyyz, tr_y_0_yyz, tr_y_0_yyzz, tr_y_0_yzz, tr_y_0_yzzz, tr_y_0_zzz, tr_y_0_zzzz, tr_y_z_xxxx, tr_y_z_xxxy, tr_y_z_xxxz, tr_y_z_xxyy, tr_y_z_xxyz, tr_y_z_xxzz, tr_y_z_xyyy, tr_y_z_xyyz, tr_y_z_xyzz, tr_y_z_xzzz, tr_y_z_yyyy, tr_y_z_yyyz, tr_y_z_yyzz, tr_y_z_yzzz, tr_y_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_z_xxxx[i] = tr_y_0_xxxx[i] * pa_z[i];

        tr_y_z_xxxy[i] = tr_y_0_xxxy[i] * pa_z[i];

        tr_y_z_xxxz[i] = tr_y_0_xxx[i] * fe_0 + tr_y_0_xxxz[i] * pa_z[i];

        tr_y_z_xxyy[i] = tr_y_0_xxyy[i] * pa_z[i];

        tr_y_z_xxyz[i] = tr_y_0_xxy[i] * fe_0 + tr_y_0_xxyz[i] * pa_z[i];

        tr_y_z_xxzz[i] = 2.0 * tr_y_0_xxz[i] * fe_0 + tr_y_0_xxzz[i] * pa_z[i];

        tr_y_z_xyyy[i] = tr_y_0_xyyy[i] * pa_z[i];

        tr_y_z_xyyz[i] = tr_y_0_xyy[i] * fe_0 + tr_y_0_xyyz[i] * pa_z[i];

        tr_y_z_xyzz[i] = 2.0 * tr_y_0_xyz[i] * fe_0 + tr_y_0_xyzz[i] * pa_z[i];

        tr_y_z_xzzz[i] = 3.0 * tr_y_0_xzz[i] * fe_0 + tr_y_0_xzzz[i] * pa_z[i];

        tr_y_z_yyyy[i] = tr_y_0_yyyy[i] * pa_z[i];

        tr_y_z_yyyz[i] = tr_y_0_yyy[i] * fe_0 + tr_y_0_yyyz[i] * pa_z[i];

        tr_y_z_yyzz[i] = 2.0 * tr_y_0_yyz[i] * fe_0 + tr_y_0_yyzz[i] * pa_z[i];

        tr_y_z_yzzz[i] = 3.0 * tr_y_0_yzz[i] * fe_0 + tr_y_0_yzzz[i] * pa_z[i];

        tr_y_z_zzzz[i] = 4.0 * tr_y_0_zzz[i] * fe_0 + tr_y_0_zzzz[i] * pa_z[i];
    }

    // Set up 90-105 components of targeted buffer : PG

    auto tr_z_x_xxxx = pbuffer.data(idx_dip_pg + 90);

    auto tr_z_x_xxxy = pbuffer.data(idx_dip_pg + 91);

    auto tr_z_x_xxxz = pbuffer.data(idx_dip_pg + 92);

    auto tr_z_x_xxyy = pbuffer.data(idx_dip_pg + 93);

    auto tr_z_x_xxyz = pbuffer.data(idx_dip_pg + 94);

    auto tr_z_x_xxzz = pbuffer.data(idx_dip_pg + 95);

    auto tr_z_x_xyyy = pbuffer.data(idx_dip_pg + 96);

    auto tr_z_x_xyyz = pbuffer.data(idx_dip_pg + 97);

    auto tr_z_x_xyzz = pbuffer.data(idx_dip_pg + 98);

    auto tr_z_x_xzzz = pbuffer.data(idx_dip_pg + 99);

    auto tr_z_x_yyyy = pbuffer.data(idx_dip_pg + 100);

    auto tr_z_x_yyyz = pbuffer.data(idx_dip_pg + 101);

    auto tr_z_x_yyzz = pbuffer.data(idx_dip_pg + 102);

    auto tr_z_x_yzzz = pbuffer.data(idx_dip_pg + 103);

    auto tr_z_x_zzzz = pbuffer.data(idx_dip_pg + 104);

    #pragma omp simd aligned(pa_x, tr_z_0_xxx, tr_z_0_xxxx, tr_z_0_xxxy, tr_z_0_xxxz, tr_z_0_xxy, tr_z_0_xxyy, tr_z_0_xxyz, tr_z_0_xxz, tr_z_0_xxzz, tr_z_0_xyy, tr_z_0_xyyy, tr_z_0_xyyz, tr_z_0_xyz, tr_z_0_xyzz, tr_z_0_xzz, tr_z_0_xzzz, tr_z_0_yyy, tr_z_0_yyyy, tr_z_0_yyyz, tr_z_0_yyz, tr_z_0_yyzz, tr_z_0_yzz, tr_z_0_yzzz, tr_z_0_zzz, tr_z_0_zzzz, tr_z_x_xxxx, tr_z_x_xxxy, tr_z_x_xxxz, tr_z_x_xxyy, tr_z_x_xxyz, tr_z_x_xxzz, tr_z_x_xyyy, tr_z_x_xyyz, tr_z_x_xyzz, tr_z_x_xzzz, tr_z_x_yyyy, tr_z_x_yyyz, tr_z_x_yyzz, tr_z_x_yzzz, tr_z_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_x_xxxx[i] = 4.0 * tr_z_0_xxx[i] * fe_0 + tr_z_0_xxxx[i] * pa_x[i];

        tr_z_x_xxxy[i] = 3.0 * tr_z_0_xxy[i] * fe_0 + tr_z_0_xxxy[i] * pa_x[i];

        tr_z_x_xxxz[i] = 3.0 * tr_z_0_xxz[i] * fe_0 + tr_z_0_xxxz[i] * pa_x[i];

        tr_z_x_xxyy[i] = 2.0 * tr_z_0_xyy[i] * fe_0 + tr_z_0_xxyy[i] * pa_x[i];

        tr_z_x_xxyz[i] = 2.0 * tr_z_0_xyz[i] * fe_0 + tr_z_0_xxyz[i] * pa_x[i];

        tr_z_x_xxzz[i] = 2.0 * tr_z_0_xzz[i] * fe_0 + tr_z_0_xxzz[i] * pa_x[i];

        tr_z_x_xyyy[i] = tr_z_0_yyy[i] * fe_0 + tr_z_0_xyyy[i] * pa_x[i];

        tr_z_x_xyyz[i] = tr_z_0_yyz[i] * fe_0 + tr_z_0_xyyz[i] * pa_x[i];

        tr_z_x_xyzz[i] = tr_z_0_yzz[i] * fe_0 + tr_z_0_xyzz[i] * pa_x[i];

        tr_z_x_xzzz[i] = tr_z_0_zzz[i] * fe_0 + tr_z_0_xzzz[i] * pa_x[i];

        tr_z_x_yyyy[i] = tr_z_0_yyyy[i] * pa_x[i];

        tr_z_x_yyyz[i] = tr_z_0_yyyz[i] * pa_x[i];

        tr_z_x_yyzz[i] = tr_z_0_yyzz[i] * pa_x[i];

        tr_z_x_yzzz[i] = tr_z_0_yzzz[i] * pa_x[i];

        tr_z_x_zzzz[i] = tr_z_0_zzzz[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : PG

    auto tr_z_y_xxxx = pbuffer.data(idx_dip_pg + 105);

    auto tr_z_y_xxxy = pbuffer.data(idx_dip_pg + 106);

    auto tr_z_y_xxxz = pbuffer.data(idx_dip_pg + 107);

    auto tr_z_y_xxyy = pbuffer.data(idx_dip_pg + 108);

    auto tr_z_y_xxyz = pbuffer.data(idx_dip_pg + 109);

    auto tr_z_y_xxzz = pbuffer.data(idx_dip_pg + 110);

    auto tr_z_y_xyyy = pbuffer.data(idx_dip_pg + 111);

    auto tr_z_y_xyyz = pbuffer.data(idx_dip_pg + 112);

    auto tr_z_y_xyzz = pbuffer.data(idx_dip_pg + 113);

    auto tr_z_y_xzzz = pbuffer.data(idx_dip_pg + 114);

    auto tr_z_y_yyyy = pbuffer.data(idx_dip_pg + 115);

    auto tr_z_y_yyyz = pbuffer.data(idx_dip_pg + 116);

    auto tr_z_y_yyzz = pbuffer.data(idx_dip_pg + 117);

    auto tr_z_y_yzzz = pbuffer.data(idx_dip_pg + 118);

    auto tr_z_y_zzzz = pbuffer.data(idx_dip_pg + 119);

    #pragma omp simd aligned(pa_y, tr_z_0_xxx, tr_z_0_xxxx, tr_z_0_xxxy, tr_z_0_xxxz, tr_z_0_xxy, tr_z_0_xxyy, tr_z_0_xxyz, tr_z_0_xxz, tr_z_0_xxzz, tr_z_0_xyy, tr_z_0_xyyy, tr_z_0_xyyz, tr_z_0_xyz, tr_z_0_xyzz, tr_z_0_xzz, tr_z_0_xzzz, tr_z_0_yyy, tr_z_0_yyyy, tr_z_0_yyyz, tr_z_0_yyz, tr_z_0_yyzz, tr_z_0_yzz, tr_z_0_yzzz, tr_z_0_zzz, tr_z_0_zzzz, tr_z_y_xxxx, tr_z_y_xxxy, tr_z_y_xxxz, tr_z_y_xxyy, tr_z_y_xxyz, tr_z_y_xxzz, tr_z_y_xyyy, tr_z_y_xyyz, tr_z_y_xyzz, tr_z_y_xzzz, tr_z_y_yyyy, tr_z_y_yyyz, tr_z_y_yyzz, tr_z_y_yzzz, tr_z_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_y_xxxx[i] = tr_z_0_xxxx[i] * pa_y[i];

        tr_z_y_xxxy[i] = tr_z_0_xxx[i] * fe_0 + tr_z_0_xxxy[i] * pa_y[i];

        tr_z_y_xxxz[i] = tr_z_0_xxxz[i] * pa_y[i];

        tr_z_y_xxyy[i] = 2.0 * tr_z_0_xxy[i] * fe_0 + tr_z_0_xxyy[i] * pa_y[i];

        tr_z_y_xxyz[i] = tr_z_0_xxz[i] * fe_0 + tr_z_0_xxyz[i] * pa_y[i];

        tr_z_y_xxzz[i] = tr_z_0_xxzz[i] * pa_y[i];

        tr_z_y_xyyy[i] = 3.0 * tr_z_0_xyy[i] * fe_0 + tr_z_0_xyyy[i] * pa_y[i];

        tr_z_y_xyyz[i] = 2.0 * tr_z_0_xyz[i] * fe_0 + tr_z_0_xyyz[i] * pa_y[i];

        tr_z_y_xyzz[i] = tr_z_0_xzz[i] * fe_0 + tr_z_0_xyzz[i] * pa_y[i];

        tr_z_y_xzzz[i] = tr_z_0_xzzz[i] * pa_y[i];

        tr_z_y_yyyy[i] = 4.0 * tr_z_0_yyy[i] * fe_0 + tr_z_0_yyyy[i] * pa_y[i];

        tr_z_y_yyyz[i] = 3.0 * tr_z_0_yyz[i] * fe_0 + tr_z_0_yyyz[i] * pa_y[i];

        tr_z_y_yyzz[i] = 2.0 * tr_z_0_yzz[i] * fe_0 + tr_z_0_yyzz[i] * pa_y[i];

        tr_z_y_yzzz[i] = tr_z_0_zzz[i] * fe_0 + tr_z_0_yzzz[i] * pa_y[i];

        tr_z_y_zzzz[i] = tr_z_0_zzzz[i] * pa_y[i];
    }

    // Set up 120-135 components of targeted buffer : PG

    auto tr_z_z_xxxx = pbuffer.data(idx_dip_pg + 120);

    auto tr_z_z_xxxy = pbuffer.data(idx_dip_pg + 121);

    auto tr_z_z_xxxz = pbuffer.data(idx_dip_pg + 122);

    auto tr_z_z_xxyy = pbuffer.data(idx_dip_pg + 123);

    auto tr_z_z_xxyz = pbuffer.data(idx_dip_pg + 124);

    auto tr_z_z_xxzz = pbuffer.data(idx_dip_pg + 125);

    auto tr_z_z_xyyy = pbuffer.data(idx_dip_pg + 126);

    auto tr_z_z_xyyz = pbuffer.data(idx_dip_pg + 127);

    auto tr_z_z_xyzz = pbuffer.data(idx_dip_pg + 128);

    auto tr_z_z_xzzz = pbuffer.data(idx_dip_pg + 129);

    auto tr_z_z_yyyy = pbuffer.data(idx_dip_pg + 130);

    auto tr_z_z_yyyz = pbuffer.data(idx_dip_pg + 131);

    auto tr_z_z_yyzz = pbuffer.data(idx_dip_pg + 132);

    auto tr_z_z_yzzz = pbuffer.data(idx_dip_pg + 133);

    auto tr_z_z_zzzz = pbuffer.data(idx_dip_pg + 134);

    #pragma omp simd aligned(pa_z, tr_z_0_xxx, tr_z_0_xxxx, tr_z_0_xxxy, tr_z_0_xxxz, tr_z_0_xxy, tr_z_0_xxyy, tr_z_0_xxyz, tr_z_0_xxz, tr_z_0_xxzz, tr_z_0_xyy, tr_z_0_xyyy, tr_z_0_xyyz, tr_z_0_xyz, tr_z_0_xyzz, tr_z_0_xzz, tr_z_0_xzzz, tr_z_0_yyy, tr_z_0_yyyy, tr_z_0_yyyz, tr_z_0_yyz, tr_z_0_yyzz, tr_z_0_yzz, tr_z_0_yzzz, tr_z_0_zzz, tr_z_0_zzzz, tr_z_z_xxxx, tr_z_z_xxxy, tr_z_z_xxxz, tr_z_z_xxyy, tr_z_z_xxyz, tr_z_z_xxzz, tr_z_z_xyyy, tr_z_z_xyyz, tr_z_z_xyzz, tr_z_z_xzzz, tr_z_z_yyyy, tr_z_z_yyyz, tr_z_z_yyzz, tr_z_z_yzzz, tr_z_z_zzzz, ts_0_xxxx, ts_0_xxxy, ts_0_xxxz, ts_0_xxyy, ts_0_xxyz, ts_0_xxzz, ts_0_xyyy, ts_0_xyyz, ts_0_xyzz, ts_0_xzzz, ts_0_yyyy, ts_0_yyyz, ts_0_yyzz, ts_0_yzzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_z_xxxx[i] = ts_0_xxxx[i] * fe_0 + tr_z_0_xxxx[i] * pa_z[i];

        tr_z_z_xxxy[i] = ts_0_xxxy[i] * fe_0 + tr_z_0_xxxy[i] * pa_z[i];

        tr_z_z_xxxz[i] = tr_z_0_xxx[i] * fe_0 + ts_0_xxxz[i] * fe_0 + tr_z_0_xxxz[i] * pa_z[i];

        tr_z_z_xxyy[i] = ts_0_xxyy[i] * fe_0 + tr_z_0_xxyy[i] * pa_z[i];

        tr_z_z_xxyz[i] = tr_z_0_xxy[i] * fe_0 + ts_0_xxyz[i] * fe_0 + tr_z_0_xxyz[i] * pa_z[i];

        tr_z_z_xxzz[i] = 2.0 * tr_z_0_xxz[i] * fe_0 + ts_0_xxzz[i] * fe_0 + tr_z_0_xxzz[i] * pa_z[i];

        tr_z_z_xyyy[i] = ts_0_xyyy[i] * fe_0 + tr_z_0_xyyy[i] * pa_z[i];

        tr_z_z_xyyz[i] = tr_z_0_xyy[i] * fe_0 + ts_0_xyyz[i] * fe_0 + tr_z_0_xyyz[i] * pa_z[i];

        tr_z_z_xyzz[i] = 2.0 * tr_z_0_xyz[i] * fe_0 + ts_0_xyzz[i] * fe_0 + tr_z_0_xyzz[i] * pa_z[i];

        tr_z_z_xzzz[i] = 3.0 * tr_z_0_xzz[i] * fe_0 + ts_0_xzzz[i] * fe_0 + tr_z_0_xzzz[i] * pa_z[i];

        tr_z_z_yyyy[i] = ts_0_yyyy[i] * fe_0 + tr_z_0_yyyy[i] * pa_z[i];

        tr_z_z_yyyz[i] = tr_z_0_yyy[i] * fe_0 + ts_0_yyyz[i] * fe_0 + tr_z_0_yyyz[i] * pa_z[i];

        tr_z_z_yyzz[i] = 2.0 * tr_z_0_yyz[i] * fe_0 + ts_0_yyzz[i] * fe_0 + tr_z_0_yyzz[i] * pa_z[i];

        tr_z_z_yzzz[i] = 3.0 * tr_z_0_yzz[i] * fe_0 + ts_0_yzzz[i] * fe_0 + tr_z_0_yzzz[i] * pa_z[i];

        tr_z_z_zzzz[i] = 4.0 * tr_z_0_zzz[i] * fe_0 + ts_0_zzzz[i] * fe_0 + tr_z_0_zzzz[i] * pa_z[i];
    }

}

} // diprec namespace

