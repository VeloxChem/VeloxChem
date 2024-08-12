#include "KineticEnergyPrimRecDD.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_dd(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_dd,
                            const size_t              idx_ovl_sd,
                            const size_t              idx_kin_sd,
                            const size_t              idx_kin_pp,
                            const size_t              idx_kin_pd,
                            const size_t              idx_ovl_dd,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto tk_0_xx = pbuffer.data(idx_kin_sd);

    auto tk_0_xy = pbuffer.data(idx_kin_sd + 1);

    auto tk_0_xz = pbuffer.data(idx_kin_sd + 2);

    auto tk_0_yy = pbuffer.data(idx_kin_sd + 3);

    auto tk_0_yz = pbuffer.data(idx_kin_sd + 4);

    auto tk_0_zz = pbuffer.data(idx_kin_sd + 5);

    // Set up components of auxiliary buffer : PP

    auto tk_x_x = pbuffer.data(idx_kin_pp);

    auto tk_x_y = pbuffer.data(idx_kin_pp + 1);

    auto tk_x_z = pbuffer.data(idx_kin_pp + 2);

    auto tk_y_x = pbuffer.data(idx_kin_pp + 3);

    auto tk_y_y = pbuffer.data(idx_kin_pp + 4);

    auto tk_y_z = pbuffer.data(idx_kin_pp + 5);

    auto tk_z_x = pbuffer.data(idx_kin_pp + 6);

    auto tk_z_y = pbuffer.data(idx_kin_pp + 7);

    auto tk_z_z = pbuffer.data(idx_kin_pp + 8);

    // Set up components of auxiliary buffer : PD

    auto tk_x_xx = pbuffer.data(idx_kin_pd);

    auto tk_x_xy = pbuffer.data(idx_kin_pd + 1);

    auto tk_x_xz = pbuffer.data(idx_kin_pd + 2);

    auto tk_x_yy = pbuffer.data(idx_kin_pd + 3);

    auto tk_x_yz = pbuffer.data(idx_kin_pd + 4);

    auto tk_x_zz = pbuffer.data(idx_kin_pd + 5);

    auto tk_y_xx = pbuffer.data(idx_kin_pd + 6);

    auto tk_y_xy = pbuffer.data(idx_kin_pd + 7);

    auto tk_y_xz = pbuffer.data(idx_kin_pd + 8);

    auto tk_y_yy = pbuffer.data(idx_kin_pd + 9);

    auto tk_y_yz = pbuffer.data(idx_kin_pd + 10);

    auto tk_y_zz = pbuffer.data(idx_kin_pd + 11);

    auto tk_z_xx = pbuffer.data(idx_kin_pd + 12);

    auto tk_z_xy = pbuffer.data(idx_kin_pd + 13);

    auto tk_z_xz = pbuffer.data(idx_kin_pd + 14);

    auto tk_z_yy = pbuffer.data(idx_kin_pd + 15);

    auto tk_z_yz = pbuffer.data(idx_kin_pd + 16);

    auto tk_z_zz = pbuffer.data(idx_kin_pd + 17);

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_ovl_dd);

    auto ts_xx_xy = pbuffer.data(idx_ovl_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_ovl_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_ovl_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_ovl_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_ovl_dd + 5);

    auto ts_xy_xx = pbuffer.data(idx_ovl_dd + 6);

    auto ts_xy_xy = pbuffer.data(idx_ovl_dd + 7);

    auto ts_xy_xz = pbuffer.data(idx_ovl_dd + 8);

    auto ts_xy_yy = pbuffer.data(idx_ovl_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_ovl_dd + 10);

    auto ts_xy_zz = pbuffer.data(idx_ovl_dd + 11);

    auto ts_xz_xx = pbuffer.data(idx_ovl_dd + 12);

    auto ts_xz_xy = pbuffer.data(idx_ovl_dd + 13);

    auto ts_xz_xz = pbuffer.data(idx_ovl_dd + 14);

    auto ts_xz_yy = pbuffer.data(idx_ovl_dd + 15);

    auto ts_xz_yz = pbuffer.data(idx_ovl_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_ovl_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_ovl_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_ovl_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_ovl_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_ovl_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_ovl_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_ovl_dd + 23);

    auto ts_yz_xx = pbuffer.data(idx_ovl_dd + 24);

    auto ts_yz_xy = pbuffer.data(idx_ovl_dd + 25);

    auto ts_yz_xz = pbuffer.data(idx_ovl_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_ovl_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_ovl_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_ovl_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_ovl_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_ovl_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_ovl_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_ovl_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_ovl_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_ovl_dd + 35);

    // Set up 0-6 components of targeted buffer : DD

    auto tk_xx_xx = pbuffer.data(idx_kin_dd);

    auto tk_xx_xy = pbuffer.data(idx_kin_dd + 1);

    auto tk_xx_xz = pbuffer.data(idx_kin_dd + 2);

    auto tk_xx_yy = pbuffer.data(idx_kin_dd + 3);

    auto tk_xx_yz = pbuffer.data(idx_kin_dd + 4);

    auto tk_xx_zz = pbuffer.data(idx_kin_dd + 5);

#pragma omp simd aligned(pa_x,         \
                             tk_0_xx,  \
                             tk_0_xy,  \
                             tk_0_xz,  \
                             tk_0_yy,  \
                             tk_0_yz,  \
                             tk_0_zz,  \
                             tk_x_x,   \
                             tk_x_xx,  \
                             tk_x_xy,  \
                             tk_x_xz,  \
                             tk_x_y,   \
                             tk_x_yy,  \
                             tk_x_yz,  \
                             tk_x_z,   \
                             tk_x_zz,  \
                             tk_xx_xx, \
                             tk_xx_xy, \
                             tk_xx_xz, \
                             tk_xx_yy, \
                             tk_xx_yz, \
                             tk_xx_zz, \
                             ts_0_xx,  \
                             ts_0_xy,  \
                             ts_0_xz,  \
                             ts_0_yy,  \
                             ts_0_yz,  \
                             ts_0_zz,  \
                             ts_xx_xx, \
                             ts_xx_xy, \
                             ts_xx_xz, \
                             ts_xx_yy, \
                             ts_xx_yz, \
                             ts_xx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_xx[i] = -2.0 * ts_0_xx[i] * fbe_0 * fz_0 + tk_0_xx[i] * fe_0 + 2.0 * tk_x_x[i] * fe_0 + tk_x_xx[i] * pa_x[i] + 2.0 * ts_xx_xx[i] * fz_0;

        tk_xx_xy[i] = -2.0 * ts_0_xy[i] * fbe_0 * fz_0 + tk_0_xy[i] * fe_0 + tk_x_y[i] * fe_0 + tk_x_xy[i] * pa_x[i] + 2.0 * ts_xx_xy[i] * fz_0;

        tk_xx_xz[i] = -2.0 * ts_0_xz[i] * fbe_0 * fz_0 + tk_0_xz[i] * fe_0 + tk_x_z[i] * fe_0 + tk_x_xz[i] * pa_x[i] + 2.0 * ts_xx_xz[i] * fz_0;

        tk_xx_yy[i] = -2.0 * ts_0_yy[i] * fbe_0 * fz_0 + tk_0_yy[i] * fe_0 + tk_x_yy[i] * pa_x[i] + 2.0 * ts_xx_yy[i] * fz_0;

        tk_xx_yz[i] = -2.0 * ts_0_yz[i] * fbe_0 * fz_0 + tk_0_yz[i] * fe_0 + tk_x_yz[i] * pa_x[i] + 2.0 * ts_xx_yz[i] * fz_0;

        tk_xx_zz[i] = -2.0 * ts_0_zz[i] * fbe_0 * fz_0 + tk_0_zz[i] * fe_0 + tk_x_zz[i] * pa_x[i] + 2.0 * ts_xx_zz[i] * fz_0;
    }

    // Set up 6-12 components of targeted buffer : DD

    auto tk_xy_xx = pbuffer.data(idx_kin_dd + 6);

    auto tk_xy_xy = pbuffer.data(idx_kin_dd + 7);

    auto tk_xy_xz = pbuffer.data(idx_kin_dd + 8);

    auto tk_xy_yy = pbuffer.data(idx_kin_dd + 9);

    auto tk_xy_yz = pbuffer.data(idx_kin_dd + 10);

    auto tk_xy_zz = pbuffer.data(idx_kin_dd + 11);

#pragma omp simd aligned(pa_x,         \
                             pa_y,     \
                             tk_x_xx,  \
                             tk_x_xz,  \
                             tk_xy_xx, \
                             tk_xy_xy, \
                             tk_xy_xz, \
                             tk_xy_yy, \
                             tk_xy_yz, \
                             tk_xy_zz, \
                             tk_y_xy,  \
                             tk_y_y,   \
                             tk_y_yy,  \
                             tk_y_yz,  \
                             tk_y_zz,  \
                             ts_xy_xx, \
                             ts_xy_xy, \
                             ts_xy_xz, \
                             ts_xy_yy, \
                             ts_xy_yz, \
                             ts_xy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xy_xx[i] = tk_x_xx[i] * pa_y[i] + 2.0 * ts_xy_xx[i] * fz_0;

        tk_xy_xy[i] = tk_y_y[i] * fe_0 + tk_y_xy[i] * pa_x[i] + 2.0 * ts_xy_xy[i] * fz_0;

        tk_xy_xz[i] = tk_x_xz[i] * pa_y[i] + 2.0 * ts_xy_xz[i] * fz_0;

        tk_xy_yy[i] = tk_y_yy[i] * pa_x[i] + 2.0 * ts_xy_yy[i] * fz_0;

        tk_xy_yz[i] = tk_y_yz[i] * pa_x[i] + 2.0 * ts_xy_yz[i] * fz_0;

        tk_xy_zz[i] = tk_y_zz[i] * pa_x[i] + 2.0 * ts_xy_zz[i] * fz_0;
    }

    // Set up 12-18 components of targeted buffer : DD

    auto tk_xz_xx = pbuffer.data(idx_kin_dd + 12);

    auto tk_xz_xy = pbuffer.data(idx_kin_dd + 13);

    auto tk_xz_xz = pbuffer.data(idx_kin_dd + 14);

    auto tk_xz_yy = pbuffer.data(idx_kin_dd + 15);

    auto tk_xz_yz = pbuffer.data(idx_kin_dd + 16);

    auto tk_xz_zz = pbuffer.data(idx_kin_dd + 17);

#pragma omp simd aligned(pa_x,         \
                             pa_z,     \
                             tk_x_xx,  \
                             tk_x_xy,  \
                             tk_xz_xx, \
                             tk_xz_xy, \
                             tk_xz_xz, \
                             tk_xz_yy, \
                             tk_xz_yz, \
                             tk_xz_zz, \
                             tk_z_xz,  \
                             tk_z_yy,  \
                             tk_z_yz,  \
                             tk_z_z,   \
                             tk_z_zz,  \
                             ts_xz_xx, \
                             ts_xz_xy, \
                             ts_xz_xz, \
                             ts_xz_yy, \
                             ts_xz_yz, \
                             ts_xz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xz_xx[i] = tk_x_xx[i] * pa_z[i] + 2.0 * ts_xz_xx[i] * fz_0;

        tk_xz_xy[i] = tk_x_xy[i] * pa_z[i] + 2.0 * ts_xz_xy[i] * fz_0;

        tk_xz_xz[i] = tk_z_z[i] * fe_0 + tk_z_xz[i] * pa_x[i] + 2.0 * ts_xz_xz[i] * fz_0;

        tk_xz_yy[i] = tk_z_yy[i] * pa_x[i] + 2.0 * ts_xz_yy[i] * fz_0;

        tk_xz_yz[i] = tk_z_yz[i] * pa_x[i] + 2.0 * ts_xz_yz[i] * fz_0;

        tk_xz_zz[i] = tk_z_zz[i] * pa_x[i] + 2.0 * ts_xz_zz[i] * fz_0;
    }

    // Set up 18-24 components of targeted buffer : DD

    auto tk_yy_xx = pbuffer.data(idx_kin_dd + 18);

    auto tk_yy_xy = pbuffer.data(idx_kin_dd + 19);

    auto tk_yy_xz = pbuffer.data(idx_kin_dd + 20);

    auto tk_yy_yy = pbuffer.data(idx_kin_dd + 21);

    auto tk_yy_yz = pbuffer.data(idx_kin_dd + 22);

    auto tk_yy_zz = pbuffer.data(idx_kin_dd + 23);

#pragma omp simd aligned(pa_y,         \
                             tk_0_xx,  \
                             tk_0_xy,  \
                             tk_0_xz,  \
                             tk_0_yy,  \
                             tk_0_yz,  \
                             tk_0_zz,  \
                             tk_y_x,   \
                             tk_y_xx,  \
                             tk_y_xy,  \
                             tk_y_xz,  \
                             tk_y_y,   \
                             tk_y_yy,  \
                             tk_y_yz,  \
                             tk_y_z,   \
                             tk_y_zz,  \
                             tk_yy_xx, \
                             tk_yy_xy, \
                             tk_yy_xz, \
                             tk_yy_yy, \
                             tk_yy_yz, \
                             tk_yy_zz, \
                             ts_0_xx,  \
                             ts_0_xy,  \
                             ts_0_xz,  \
                             ts_0_yy,  \
                             ts_0_yz,  \
                             ts_0_zz,  \
                             ts_yy_xx, \
                             ts_yy_xy, \
                             ts_yy_xz, \
                             ts_yy_yy, \
                             ts_yy_yz, \
                             ts_yy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yy_xx[i] = -2.0 * ts_0_xx[i] * fbe_0 * fz_0 + tk_0_xx[i] * fe_0 + tk_y_xx[i] * pa_y[i] + 2.0 * ts_yy_xx[i] * fz_0;

        tk_yy_xy[i] = -2.0 * ts_0_xy[i] * fbe_0 * fz_0 + tk_0_xy[i] * fe_0 + tk_y_x[i] * fe_0 + tk_y_xy[i] * pa_y[i] + 2.0 * ts_yy_xy[i] * fz_0;

        tk_yy_xz[i] = -2.0 * ts_0_xz[i] * fbe_0 * fz_0 + tk_0_xz[i] * fe_0 + tk_y_xz[i] * pa_y[i] + 2.0 * ts_yy_xz[i] * fz_0;

        tk_yy_yy[i] = -2.0 * ts_0_yy[i] * fbe_0 * fz_0 + tk_0_yy[i] * fe_0 + 2.0 * tk_y_y[i] * fe_0 + tk_y_yy[i] * pa_y[i] + 2.0 * ts_yy_yy[i] * fz_0;

        tk_yy_yz[i] = -2.0 * ts_0_yz[i] * fbe_0 * fz_0 + tk_0_yz[i] * fe_0 + tk_y_z[i] * fe_0 + tk_y_yz[i] * pa_y[i] + 2.0 * ts_yy_yz[i] * fz_0;

        tk_yy_zz[i] = -2.0 * ts_0_zz[i] * fbe_0 * fz_0 + tk_0_zz[i] * fe_0 + tk_y_zz[i] * pa_y[i] + 2.0 * ts_yy_zz[i] * fz_0;
    }

    // Set up 24-30 components of targeted buffer : DD

    auto tk_yz_xx = pbuffer.data(idx_kin_dd + 24);

    auto tk_yz_xy = pbuffer.data(idx_kin_dd + 25);

    auto tk_yz_xz = pbuffer.data(idx_kin_dd + 26);

    auto tk_yz_yy = pbuffer.data(idx_kin_dd + 27);

    auto tk_yz_yz = pbuffer.data(idx_kin_dd + 28);

    auto tk_yz_zz = pbuffer.data(idx_kin_dd + 29);

#pragma omp simd aligned(pa_y,         \
                             pa_z,     \
                             tk_y_xy,  \
                             tk_y_yy,  \
                             tk_yz_xx, \
                             tk_yz_xy, \
                             tk_yz_xz, \
                             tk_yz_yy, \
                             tk_yz_yz, \
                             tk_yz_zz, \
                             tk_z_xx,  \
                             tk_z_xz,  \
                             tk_z_yz,  \
                             tk_z_z,   \
                             tk_z_zz,  \
                             ts_yz_xx, \
                             ts_yz_xy, \
                             ts_yz_xz, \
                             ts_yz_yy, \
                             ts_yz_yz, \
                             ts_yz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yz_xx[i] = tk_z_xx[i] * pa_y[i] + 2.0 * ts_yz_xx[i] * fz_0;

        tk_yz_xy[i] = tk_y_xy[i] * pa_z[i] + 2.0 * ts_yz_xy[i] * fz_0;

        tk_yz_xz[i] = tk_z_xz[i] * pa_y[i] + 2.0 * ts_yz_xz[i] * fz_0;

        tk_yz_yy[i] = tk_y_yy[i] * pa_z[i] + 2.0 * ts_yz_yy[i] * fz_0;

        tk_yz_yz[i] = tk_z_z[i] * fe_0 + tk_z_yz[i] * pa_y[i] + 2.0 * ts_yz_yz[i] * fz_0;

        tk_yz_zz[i] = tk_z_zz[i] * pa_y[i] + 2.0 * ts_yz_zz[i] * fz_0;
    }

    // Set up 30-36 components of targeted buffer : DD

    auto tk_zz_xx = pbuffer.data(idx_kin_dd + 30);

    auto tk_zz_xy = pbuffer.data(idx_kin_dd + 31);

    auto tk_zz_xz = pbuffer.data(idx_kin_dd + 32);

    auto tk_zz_yy = pbuffer.data(idx_kin_dd + 33);

    auto tk_zz_yz = pbuffer.data(idx_kin_dd + 34);

    auto tk_zz_zz = pbuffer.data(idx_kin_dd + 35);

#pragma omp simd aligned(pa_z,         \
                             tk_0_xx,  \
                             tk_0_xy,  \
                             tk_0_xz,  \
                             tk_0_yy,  \
                             tk_0_yz,  \
                             tk_0_zz,  \
                             tk_z_x,   \
                             tk_z_xx,  \
                             tk_z_xy,  \
                             tk_z_xz,  \
                             tk_z_y,   \
                             tk_z_yy,  \
                             tk_z_yz,  \
                             tk_z_z,   \
                             tk_z_zz,  \
                             tk_zz_xx, \
                             tk_zz_xy, \
                             tk_zz_xz, \
                             tk_zz_yy, \
                             tk_zz_yz, \
                             tk_zz_zz, \
                             ts_0_xx,  \
                             ts_0_xy,  \
                             ts_0_xz,  \
                             ts_0_yy,  \
                             ts_0_yz,  \
                             ts_0_zz,  \
                             ts_zz_xx, \
                             ts_zz_xy, \
                             ts_zz_xz, \
                             ts_zz_yy, \
                             ts_zz_yz, \
                             ts_zz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zz_xx[i] = -2.0 * ts_0_xx[i] * fbe_0 * fz_0 + tk_0_xx[i] * fe_0 + tk_z_xx[i] * pa_z[i] + 2.0 * ts_zz_xx[i] * fz_0;

        tk_zz_xy[i] = -2.0 * ts_0_xy[i] * fbe_0 * fz_0 + tk_0_xy[i] * fe_0 + tk_z_xy[i] * pa_z[i] + 2.0 * ts_zz_xy[i] * fz_0;

        tk_zz_xz[i] = -2.0 * ts_0_xz[i] * fbe_0 * fz_0 + tk_0_xz[i] * fe_0 + tk_z_x[i] * fe_0 + tk_z_xz[i] * pa_z[i] + 2.0 * ts_zz_xz[i] * fz_0;

        tk_zz_yy[i] = -2.0 * ts_0_yy[i] * fbe_0 * fz_0 + tk_0_yy[i] * fe_0 + tk_z_yy[i] * pa_z[i] + 2.0 * ts_zz_yy[i] * fz_0;

        tk_zz_yz[i] = -2.0 * ts_0_yz[i] * fbe_0 * fz_0 + tk_0_yz[i] * fe_0 + tk_z_y[i] * fe_0 + tk_z_yz[i] * pa_z[i] + 2.0 * ts_zz_yz[i] * fz_0;

        tk_zz_zz[i] = -2.0 * ts_0_zz[i] * fbe_0 * fz_0 + tk_0_zz[i] * fe_0 + 2.0 * tk_z_z[i] * fe_0 + tk_z_zz[i] * pa_z[i] + 2.0 * ts_zz_zz[i] * fz_0;
    }
}

}  // namespace kinrec
