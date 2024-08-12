#include "ElectricDipoleMomentumPrimRecDD.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_dd(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_dd,
                                      const size_t              idx_dip_sd,
                                      const size_t              idx_dip_pp,
                                      const size_t              idx_ovl_pd,
                                      const size_t              idx_dip_pd,
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

    auto tr_x_0_xx = pbuffer.data(idx_dip_sd);

    auto tr_x_0_xy = pbuffer.data(idx_dip_sd + 1);

    auto tr_x_0_xz = pbuffer.data(idx_dip_sd + 2);

    auto tr_x_0_yy = pbuffer.data(idx_dip_sd + 3);

    auto tr_x_0_yz = pbuffer.data(idx_dip_sd + 4);

    auto tr_x_0_zz = pbuffer.data(idx_dip_sd + 5);

    auto tr_y_0_xx = pbuffer.data(idx_dip_sd + 6);

    auto tr_y_0_xy = pbuffer.data(idx_dip_sd + 7);

    auto tr_y_0_xz = pbuffer.data(idx_dip_sd + 8);

    auto tr_y_0_yy = pbuffer.data(idx_dip_sd + 9);

    auto tr_y_0_yz = pbuffer.data(idx_dip_sd + 10);

    auto tr_y_0_zz = pbuffer.data(idx_dip_sd + 11);

    auto tr_z_0_xx = pbuffer.data(idx_dip_sd + 12);

    auto tr_z_0_xy = pbuffer.data(idx_dip_sd + 13);

    auto tr_z_0_xz = pbuffer.data(idx_dip_sd + 14);

    auto tr_z_0_yy = pbuffer.data(idx_dip_sd + 15);

    auto tr_z_0_yz = pbuffer.data(idx_dip_sd + 16);

    auto tr_z_0_zz = pbuffer.data(idx_dip_sd + 17);

    // Set up components of auxiliary buffer : PP

    auto tr_x_x_x = pbuffer.data(idx_dip_pp);

    auto tr_x_x_y = pbuffer.data(idx_dip_pp + 1);

    auto tr_x_x_z = pbuffer.data(idx_dip_pp + 2);

    auto tr_x_y_x = pbuffer.data(idx_dip_pp + 3);

    auto tr_x_y_y = pbuffer.data(idx_dip_pp + 4);

    auto tr_x_y_z = pbuffer.data(idx_dip_pp + 5);

    auto tr_x_z_x = pbuffer.data(idx_dip_pp + 6);

    auto tr_x_z_y = pbuffer.data(idx_dip_pp + 7);

    auto tr_x_z_z = pbuffer.data(idx_dip_pp + 8);

    auto tr_y_x_x = pbuffer.data(idx_dip_pp + 9);

    auto tr_y_x_y = pbuffer.data(idx_dip_pp + 10);

    auto tr_y_x_z = pbuffer.data(idx_dip_pp + 11);

    auto tr_y_y_x = pbuffer.data(idx_dip_pp + 12);

    auto tr_y_y_y = pbuffer.data(idx_dip_pp + 13);

    auto tr_y_y_z = pbuffer.data(idx_dip_pp + 14);

    auto tr_y_z_x = pbuffer.data(idx_dip_pp + 15);

    auto tr_y_z_y = pbuffer.data(idx_dip_pp + 16);

    auto tr_y_z_z = pbuffer.data(idx_dip_pp + 17);

    auto tr_z_x_x = pbuffer.data(idx_dip_pp + 18);

    auto tr_z_x_y = pbuffer.data(idx_dip_pp + 19);

    auto tr_z_x_z = pbuffer.data(idx_dip_pp + 20);

    auto tr_z_y_x = pbuffer.data(idx_dip_pp + 21);

    auto tr_z_y_y = pbuffer.data(idx_dip_pp + 22);

    auto tr_z_y_z = pbuffer.data(idx_dip_pp + 23);

    auto tr_z_z_x = pbuffer.data(idx_dip_pp + 24);

    auto tr_z_z_y = pbuffer.data(idx_dip_pp + 25);

    auto tr_z_z_z = pbuffer.data(idx_dip_pp + 26);

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_ovl_pd);

    auto ts_x_xy = pbuffer.data(idx_ovl_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_ovl_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_ovl_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_ovl_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_ovl_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_ovl_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_ovl_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_ovl_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_ovl_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_ovl_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_ovl_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_ovl_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_ovl_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_ovl_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_ovl_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_ovl_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_ovl_pd + 17);

    // Set up components of auxiliary buffer : PD

    auto tr_x_x_xx = pbuffer.data(idx_dip_pd);

    auto tr_x_x_xy = pbuffer.data(idx_dip_pd + 1);

    auto tr_x_x_xz = pbuffer.data(idx_dip_pd + 2);

    auto tr_x_x_yy = pbuffer.data(idx_dip_pd + 3);

    auto tr_x_x_yz = pbuffer.data(idx_dip_pd + 4);

    auto tr_x_x_zz = pbuffer.data(idx_dip_pd + 5);

    auto tr_x_y_xx = pbuffer.data(idx_dip_pd + 6);

    auto tr_x_y_xy = pbuffer.data(idx_dip_pd + 7);

    auto tr_x_y_xz = pbuffer.data(idx_dip_pd + 8);

    auto tr_x_y_yy = pbuffer.data(idx_dip_pd + 9);

    auto tr_x_y_yz = pbuffer.data(idx_dip_pd + 10);

    auto tr_x_y_zz = pbuffer.data(idx_dip_pd + 11);

    auto tr_x_z_xx = pbuffer.data(idx_dip_pd + 12);

    auto tr_x_z_xy = pbuffer.data(idx_dip_pd + 13);

    auto tr_x_z_xz = pbuffer.data(idx_dip_pd + 14);

    auto tr_x_z_yy = pbuffer.data(idx_dip_pd + 15);

    auto tr_x_z_yz = pbuffer.data(idx_dip_pd + 16);

    auto tr_x_z_zz = pbuffer.data(idx_dip_pd + 17);

    auto tr_y_x_xx = pbuffer.data(idx_dip_pd + 18);

    auto tr_y_x_xy = pbuffer.data(idx_dip_pd + 19);

    auto tr_y_x_xz = pbuffer.data(idx_dip_pd + 20);

    auto tr_y_x_yy = pbuffer.data(idx_dip_pd + 21);

    auto tr_y_x_yz = pbuffer.data(idx_dip_pd + 22);

    auto tr_y_x_zz = pbuffer.data(idx_dip_pd + 23);

    auto tr_y_y_xx = pbuffer.data(idx_dip_pd + 24);

    auto tr_y_y_xy = pbuffer.data(idx_dip_pd + 25);

    auto tr_y_y_xz = pbuffer.data(idx_dip_pd + 26);

    auto tr_y_y_yy = pbuffer.data(idx_dip_pd + 27);

    auto tr_y_y_yz = pbuffer.data(idx_dip_pd + 28);

    auto tr_y_y_zz = pbuffer.data(idx_dip_pd + 29);

    auto tr_y_z_xx = pbuffer.data(idx_dip_pd + 30);

    auto tr_y_z_xy = pbuffer.data(idx_dip_pd + 31);

    auto tr_y_z_xz = pbuffer.data(idx_dip_pd + 32);

    auto tr_y_z_yy = pbuffer.data(idx_dip_pd + 33);

    auto tr_y_z_yz = pbuffer.data(idx_dip_pd + 34);

    auto tr_y_z_zz = pbuffer.data(idx_dip_pd + 35);

    auto tr_z_x_xx = pbuffer.data(idx_dip_pd + 36);

    auto tr_z_x_xy = pbuffer.data(idx_dip_pd + 37);

    auto tr_z_x_xz = pbuffer.data(idx_dip_pd + 38);

    auto tr_z_x_yy = pbuffer.data(idx_dip_pd + 39);

    auto tr_z_x_yz = pbuffer.data(idx_dip_pd + 40);

    auto tr_z_x_zz = pbuffer.data(idx_dip_pd + 41);

    auto tr_z_y_xx = pbuffer.data(idx_dip_pd + 42);

    auto tr_z_y_xy = pbuffer.data(idx_dip_pd + 43);

    auto tr_z_y_xz = pbuffer.data(idx_dip_pd + 44);

    auto tr_z_y_yy = pbuffer.data(idx_dip_pd + 45);

    auto tr_z_y_yz = pbuffer.data(idx_dip_pd + 46);

    auto tr_z_y_zz = pbuffer.data(idx_dip_pd + 47);

    auto tr_z_z_xx = pbuffer.data(idx_dip_pd + 48);

    auto tr_z_z_xy = pbuffer.data(idx_dip_pd + 49);

    auto tr_z_z_xz = pbuffer.data(idx_dip_pd + 50);

    auto tr_z_z_yy = pbuffer.data(idx_dip_pd + 51);

    auto tr_z_z_yz = pbuffer.data(idx_dip_pd + 52);

    auto tr_z_z_zz = pbuffer.data(idx_dip_pd + 53);

    // Set up 0-6 components of targeted buffer : DD

    auto tr_x_xx_xx = pbuffer.data(idx_dip_dd);

    auto tr_x_xx_xy = pbuffer.data(idx_dip_dd + 1);

    auto tr_x_xx_xz = pbuffer.data(idx_dip_dd + 2);

    auto tr_x_xx_yy = pbuffer.data(idx_dip_dd + 3);

    auto tr_x_xx_yz = pbuffer.data(idx_dip_dd + 4);

    auto tr_x_xx_zz = pbuffer.data(idx_dip_dd + 5);

#pragma omp simd aligned(pa_x,           \
                             tr_x_0_xx,  \
                             tr_x_0_xy,  \
                             tr_x_0_xz,  \
                             tr_x_0_yy,  \
                             tr_x_0_yz,  \
                             tr_x_0_zz,  \
                             tr_x_x_x,   \
                             tr_x_x_xx,  \
                             tr_x_x_xy,  \
                             tr_x_x_xz,  \
                             tr_x_x_y,   \
                             tr_x_x_yy,  \
                             tr_x_x_yz,  \
                             tr_x_x_z,   \
                             tr_x_x_zz,  \
                             tr_x_xx_xx, \
                             tr_x_xx_xy, \
                             tr_x_xx_xz, \
                             tr_x_xx_yy, \
                             tr_x_xx_yz, \
                             tr_x_xx_zz, \
                             ts_x_xx,    \
                             ts_x_xy,    \
                             ts_x_xz,    \
                             ts_x_yy,    \
                             ts_x_yz,    \
                             ts_x_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_xx[i] = tr_x_0_xx[i] * fe_0 + 2.0 * tr_x_x_x[i] * fe_0 + ts_x_xx[i] * fe_0 + tr_x_x_xx[i] * pa_x[i];

        tr_x_xx_xy[i] = tr_x_0_xy[i] * fe_0 + tr_x_x_y[i] * fe_0 + ts_x_xy[i] * fe_0 + tr_x_x_xy[i] * pa_x[i];

        tr_x_xx_xz[i] = tr_x_0_xz[i] * fe_0 + tr_x_x_z[i] * fe_0 + ts_x_xz[i] * fe_0 + tr_x_x_xz[i] * pa_x[i];

        tr_x_xx_yy[i] = tr_x_0_yy[i] * fe_0 + ts_x_yy[i] * fe_0 + tr_x_x_yy[i] * pa_x[i];

        tr_x_xx_yz[i] = tr_x_0_yz[i] * fe_0 + ts_x_yz[i] * fe_0 + tr_x_x_yz[i] * pa_x[i];

        tr_x_xx_zz[i] = tr_x_0_zz[i] * fe_0 + ts_x_zz[i] * fe_0 + tr_x_x_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto tr_x_xy_xx = pbuffer.data(idx_dip_dd + 6);

    auto tr_x_xy_xy = pbuffer.data(idx_dip_dd + 7);

    auto tr_x_xy_xz = pbuffer.data(idx_dip_dd + 8);

    auto tr_x_xy_yy = pbuffer.data(idx_dip_dd + 9);

    auto tr_x_xy_yz = pbuffer.data(idx_dip_dd + 10);

    auto tr_x_xy_zz = pbuffer.data(idx_dip_dd + 11);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tr_x_x_x,   \
                             tr_x_x_xx,  \
                             tr_x_x_xy,  \
                             tr_x_x_xz,  \
                             tr_x_x_zz,  \
                             tr_x_xy_xx, \
                             tr_x_xy_xy, \
                             tr_x_xy_xz, \
                             tr_x_xy_yy, \
                             tr_x_xy_yz, \
                             tr_x_xy_zz, \
                             tr_x_y_yy,  \
                             tr_x_y_yz,  \
                             ts_y_yy,    \
                             ts_y_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xy_xx[i] = tr_x_x_xx[i] * pa_y[i];

        tr_x_xy_xy[i] = tr_x_x_x[i] * fe_0 + tr_x_x_xy[i] * pa_y[i];

        tr_x_xy_xz[i] = tr_x_x_xz[i] * pa_y[i];

        tr_x_xy_yy[i] = ts_y_yy[i] * fe_0 + tr_x_y_yy[i] * pa_x[i];

        tr_x_xy_yz[i] = ts_y_yz[i] * fe_0 + tr_x_y_yz[i] * pa_x[i];

        tr_x_xy_zz[i] = tr_x_x_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto tr_x_xz_xx = pbuffer.data(idx_dip_dd + 12);

    auto tr_x_xz_xy = pbuffer.data(idx_dip_dd + 13);

    auto tr_x_xz_xz = pbuffer.data(idx_dip_dd + 14);

    auto tr_x_xz_yy = pbuffer.data(idx_dip_dd + 15);

    auto tr_x_xz_yz = pbuffer.data(idx_dip_dd + 16);

    auto tr_x_xz_zz = pbuffer.data(idx_dip_dd + 17);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tr_x_x_x,   \
                             tr_x_x_xx,  \
                             tr_x_x_xy,  \
                             tr_x_x_xz,  \
                             tr_x_x_yy,  \
                             tr_x_xz_xx, \
                             tr_x_xz_xy, \
                             tr_x_xz_xz, \
                             tr_x_xz_yy, \
                             tr_x_xz_yz, \
                             tr_x_xz_zz, \
                             tr_x_z_yz,  \
                             tr_x_z_zz,  \
                             ts_z_yz,    \
                             ts_z_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xz_xx[i] = tr_x_x_xx[i] * pa_z[i];

        tr_x_xz_xy[i] = tr_x_x_xy[i] * pa_z[i];

        tr_x_xz_xz[i] = tr_x_x_x[i] * fe_0 + tr_x_x_xz[i] * pa_z[i];

        tr_x_xz_yy[i] = tr_x_x_yy[i] * pa_z[i];

        tr_x_xz_yz[i] = ts_z_yz[i] * fe_0 + tr_x_z_yz[i] * pa_x[i];

        tr_x_xz_zz[i] = ts_z_zz[i] * fe_0 + tr_x_z_zz[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto tr_x_yy_xx = pbuffer.data(idx_dip_dd + 18);

    auto tr_x_yy_xy = pbuffer.data(idx_dip_dd + 19);

    auto tr_x_yy_xz = pbuffer.data(idx_dip_dd + 20);

    auto tr_x_yy_yy = pbuffer.data(idx_dip_dd + 21);

    auto tr_x_yy_yz = pbuffer.data(idx_dip_dd + 22);

    auto tr_x_yy_zz = pbuffer.data(idx_dip_dd + 23);

#pragma omp simd aligned(pa_y,           \
                             tr_x_0_xx,  \
                             tr_x_0_xy,  \
                             tr_x_0_xz,  \
                             tr_x_0_yy,  \
                             tr_x_0_yz,  \
                             tr_x_0_zz,  \
                             tr_x_y_x,   \
                             tr_x_y_xx,  \
                             tr_x_y_xy,  \
                             tr_x_y_xz,  \
                             tr_x_y_y,   \
                             tr_x_y_yy,  \
                             tr_x_y_yz,  \
                             tr_x_y_z,   \
                             tr_x_y_zz,  \
                             tr_x_yy_xx, \
                             tr_x_yy_xy, \
                             tr_x_yy_xz, \
                             tr_x_yy_yy, \
                             tr_x_yy_yz, \
                             tr_x_yy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yy_xx[i] = tr_x_0_xx[i] * fe_0 + tr_x_y_xx[i] * pa_y[i];

        tr_x_yy_xy[i] = tr_x_0_xy[i] * fe_0 + tr_x_y_x[i] * fe_0 + tr_x_y_xy[i] * pa_y[i];

        tr_x_yy_xz[i] = tr_x_0_xz[i] * fe_0 + tr_x_y_xz[i] * pa_y[i];

        tr_x_yy_yy[i] = tr_x_0_yy[i] * fe_0 + 2.0 * tr_x_y_y[i] * fe_0 + tr_x_y_yy[i] * pa_y[i];

        tr_x_yy_yz[i] = tr_x_0_yz[i] * fe_0 + tr_x_y_z[i] * fe_0 + tr_x_y_yz[i] * pa_y[i];

        tr_x_yy_zz[i] = tr_x_0_zz[i] * fe_0 + tr_x_y_zz[i] * pa_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto tr_x_yz_xx = pbuffer.data(idx_dip_dd + 24);

    auto tr_x_yz_xy = pbuffer.data(idx_dip_dd + 25);

    auto tr_x_yz_xz = pbuffer.data(idx_dip_dd + 26);

    auto tr_x_yz_yy = pbuffer.data(idx_dip_dd + 27);

    auto tr_x_yz_yz = pbuffer.data(idx_dip_dd + 28);

    auto tr_x_yz_zz = pbuffer.data(idx_dip_dd + 29);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tr_x_y_xy,  \
                             tr_x_y_yy,  \
                             tr_x_yz_xx, \
                             tr_x_yz_xy, \
                             tr_x_yz_xz, \
                             tr_x_yz_yy, \
                             tr_x_yz_yz, \
                             tr_x_yz_zz, \
                             tr_x_z_xx,  \
                             tr_x_z_xz,  \
                             tr_x_z_yz,  \
                             tr_x_z_z,   \
                             tr_x_z_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yz_xx[i] = tr_x_z_xx[i] * pa_y[i];

        tr_x_yz_xy[i] = tr_x_y_xy[i] * pa_z[i];

        tr_x_yz_xz[i] = tr_x_z_xz[i] * pa_y[i];

        tr_x_yz_yy[i] = tr_x_y_yy[i] * pa_z[i];

        tr_x_yz_yz[i] = tr_x_z_z[i] * fe_0 + tr_x_z_yz[i] * pa_y[i];

        tr_x_yz_zz[i] = tr_x_z_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto tr_x_zz_xx = pbuffer.data(idx_dip_dd + 30);

    auto tr_x_zz_xy = pbuffer.data(idx_dip_dd + 31);

    auto tr_x_zz_xz = pbuffer.data(idx_dip_dd + 32);

    auto tr_x_zz_yy = pbuffer.data(idx_dip_dd + 33);

    auto tr_x_zz_yz = pbuffer.data(idx_dip_dd + 34);

    auto tr_x_zz_zz = pbuffer.data(idx_dip_dd + 35);

#pragma omp simd aligned(pa_z,           \
                             tr_x_0_xx,  \
                             tr_x_0_xy,  \
                             tr_x_0_xz,  \
                             tr_x_0_yy,  \
                             tr_x_0_yz,  \
                             tr_x_0_zz,  \
                             tr_x_z_x,   \
                             tr_x_z_xx,  \
                             tr_x_z_xy,  \
                             tr_x_z_xz,  \
                             tr_x_z_y,   \
                             tr_x_z_yy,  \
                             tr_x_z_yz,  \
                             tr_x_z_z,   \
                             tr_x_z_zz,  \
                             tr_x_zz_xx, \
                             tr_x_zz_xy, \
                             tr_x_zz_xz, \
                             tr_x_zz_yy, \
                             tr_x_zz_yz, \
                             tr_x_zz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zz_xx[i] = tr_x_0_xx[i] * fe_0 + tr_x_z_xx[i] * pa_z[i];

        tr_x_zz_xy[i] = tr_x_0_xy[i] * fe_0 + tr_x_z_xy[i] * pa_z[i];

        tr_x_zz_xz[i] = tr_x_0_xz[i] * fe_0 + tr_x_z_x[i] * fe_0 + tr_x_z_xz[i] * pa_z[i];

        tr_x_zz_yy[i] = tr_x_0_yy[i] * fe_0 + tr_x_z_yy[i] * pa_z[i];

        tr_x_zz_yz[i] = tr_x_0_yz[i] * fe_0 + tr_x_z_y[i] * fe_0 + tr_x_z_yz[i] * pa_z[i];

        tr_x_zz_zz[i] = tr_x_0_zz[i] * fe_0 + 2.0 * tr_x_z_z[i] * fe_0 + tr_x_z_zz[i] * pa_z[i];
    }

    // Set up 36-42 components of targeted buffer : DD

    auto tr_y_xx_xx = pbuffer.data(idx_dip_dd + 36);

    auto tr_y_xx_xy = pbuffer.data(idx_dip_dd + 37);

    auto tr_y_xx_xz = pbuffer.data(idx_dip_dd + 38);

    auto tr_y_xx_yy = pbuffer.data(idx_dip_dd + 39);

    auto tr_y_xx_yz = pbuffer.data(idx_dip_dd + 40);

    auto tr_y_xx_zz = pbuffer.data(idx_dip_dd + 41);

#pragma omp simd aligned(pa_x,           \
                             tr_y_0_xx,  \
                             tr_y_0_xy,  \
                             tr_y_0_xz,  \
                             tr_y_0_yy,  \
                             tr_y_0_yz,  \
                             tr_y_0_zz,  \
                             tr_y_x_x,   \
                             tr_y_x_xx,  \
                             tr_y_x_xy,  \
                             tr_y_x_xz,  \
                             tr_y_x_y,   \
                             tr_y_x_yy,  \
                             tr_y_x_yz,  \
                             tr_y_x_z,   \
                             tr_y_x_zz,  \
                             tr_y_xx_xx, \
                             tr_y_xx_xy, \
                             tr_y_xx_xz, \
                             tr_y_xx_yy, \
                             tr_y_xx_yz, \
                             tr_y_xx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xx_xx[i] = tr_y_0_xx[i] * fe_0 + 2.0 * tr_y_x_x[i] * fe_0 + tr_y_x_xx[i] * pa_x[i];

        tr_y_xx_xy[i] = tr_y_0_xy[i] * fe_0 + tr_y_x_y[i] * fe_0 + tr_y_x_xy[i] * pa_x[i];

        tr_y_xx_xz[i] = tr_y_0_xz[i] * fe_0 + tr_y_x_z[i] * fe_0 + tr_y_x_xz[i] * pa_x[i];

        tr_y_xx_yy[i] = tr_y_0_yy[i] * fe_0 + tr_y_x_yy[i] * pa_x[i];

        tr_y_xx_yz[i] = tr_y_0_yz[i] * fe_0 + tr_y_x_yz[i] * pa_x[i];

        tr_y_xx_zz[i] = tr_y_0_zz[i] * fe_0 + tr_y_x_zz[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : DD

    auto tr_y_xy_xx = pbuffer.data(idx_dip_dd + 42);

    auto tr_y_xy_xy = pbuffer.data(idx_dip_dd + 43);

    auto tr_y_xy_xz = pbuffer.data(idx_dip_dd + 44);

    auto tr_y_xy_yy = pbuffer.data(idx_dip_dd + 45);

    auto tr_y_xy_yz = pbuffer.data(idx_dip_dd + 46);

    auto tr_y_xy_zz = pbuffer.data(idx_dip_dd + 47);

#pragma omp simd aligned(pa_x,           \
                             tr_y_xy_xx, \
                             tr_y_xy_xy, \
                             tr_y_xy_xz, \
                             tr_y_xy_yy, \
                             tr_y_xy_yz, \
                             tr_y_xy_zz, \
                             tr_y_y_x,   \
                             tr_y_y_xx,  \
                             tr_y_y_xy,  \
                             tr_y_y_xz,  \
                             tr_y_y_y,   \
                             tr_y_y_yy,  \
                             tr_y_y_yz,  \
                             tr_y_y_z,   \
                             tr_y_y_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xy_xx[i] = 2.0 * tr_y_y_x[i] * fe_0 + tr_y_y_xx[i] * pa_x[i];

        tr_y_xy_xy[i] = tr_y_y_y[i] * fe_0 + tr_y_y_xy[i] * pa_x[i];

        tr_y_xy_xz[i] = tr_y_y_z[i] * fe_0 + tr_y_y_xz[i] * pa_x[i];

        tr_y_xy_yy[i] = tr_y_y_yy[i] * pa_x[i];

        tr_y_xy_yz[i] = tr_y_y_yz[i] * pa_x[i];

        tr_y_xy_zz[i] = tr_y_y_zz[i] * pa_x[i];
    }

    // Set up 48-54 components of targeted buffer : DD

    auto tr_y_xz_xx = pbuffer.data(idx_dip_dd + 48);

    auto tr_y_xz_xy = pbuffer.data(idx_dip_dd + 49);

    auto tr_y_xz_xz = pbuffer.data(idx_dip_dd + 50);

    auto tr_y_xz_yy = pbuffer.data(idx_dip_dd + 51);

    auto tr_y_xz_yz = pbuffer.data(idx_dip_dd + 52);

    auto tr_y_xz_zz = pbuffer.data(idx_dip_dd + 53);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tr_y_x_xx,  \
                             tr_y_x_xy,  \
                             tr_y_xz_xx, \
                             tr_y_xz_xy, \
                             tr_y_xz_xz, \
                             tr_y_xz_yy, \
                             tr_y_xz_yz, \
                             tr_y_xz_zz, \
                             tr_y_z_xz,  \
                             tr_y_z_yy,  \
                             tr_y_z_yz,  \
                             tr_y_z_z,   \
                             tr_y_z_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xz_xx[i] = tr_y_x_xx[i] * pa_z[i];

        tr_y_xz_xy[i] = tr_y_x_xy[i] * pa_z[i];

        tr_y_xz_xz[i] = tr_y_z_z[i] * fe_0 + tr_y_z_xz[i] * pa_x[i];

        tr_y_xz_yy[i] = tr_y_z_yy[i] * pa_x[i];

        tr_y_xz_yz[i] = tr_y_z_yz[i] * pa_x[i];

        tr_y_xz_zz[i] = tr_y_z_zz[i] * pa_x[i];
    }

    // Set up 54-60 components of targeted buffer : DD

    auto tr_y_yy_xx = pbuffer.data(idx_dip_dd + 54);

    auto tr_y_yy_xy = pbuffer.data(idx_dip_dd + 55);

    auto tr_y_yy_xz = pbuffer.data(idx_dip_dd + 56);

    auto tr_y_yy_yy = pbuffer.data(idx_dip_dd + 57);

    auto tr_y_yy_yz = pbuffer.data(idx_dip_dd + 58);

    auto tr_y_yy_zz = pbuffer.data(idx_dip_dd + 59);

#pragma omp simd aligned(pa_y,           \
                             tr_y_0_xx,  \
                             tr_y_0_xy,  \
                             tr_y_0_xz,  \
                             tr_y_0_yy,  \
                             tr_y_0_yz,  \
                             tr_y_0_zz,  \
                             tr_y_y_x,   \
                             tr_y_y_xx,  \
                             tr_y_y_xy,  \
                             tr_y_y_xz,  \
                             tr_y_y_y,   \
                             tr_y_y_yy,  \
                             tr_y_y_yz,  \
                             tr_y_y_z,   \
                             tr_y_y_zz,  \
                             tr_y_yy_xx, \
                             tr_y_yy_xy, \
                             tr_y_yy_xz, \
                             tr_y_yy_yy, \
                             tr_y_yy_yz, \
                             tr_y_yy_zz, \
                             ts_y_xx,    \
                             ts_y_xy,    \
                             ts_y_xz,    \
                             ts_y_yy,    \
                             ts_y_yz,    \
                             ts_y_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yy_xx[i] = tr_y_0_xx[i] * fe_0 + ts_y_xx[i] * fe_0 + tr_y_y_xx[i] * pa_y[i];

        tr_y_yy_xy[i] = tr_y_0_xy[i] * fe_0 + tr_y_y_x[i] * fe_0 + ts_y_xy[i] * fe_0 + tr_y_y_xy[i] * pa_y[i];

        tr_y_yy_xz[i] = tr_y_0_xz[i] * fe_0 + ts_y_xz[i] * fe_0 + tr_y_y_xz[i] * pa_y[i];

        tr_y_yy_yy[i] = tr_y_0_yy[i] * fe_0 + 2.0 * tr_y_y_y[i] * fe_0 + ts_y_yy[i] * fe_0 + tr_y_y_yy[i] * pa_y[i];

        tr_y_yy_yz[i] = tr_y_0_yz[i] * fe_0 + tr_y_y_z[i] * fe_0 + ts_y_yz[i] * fe_0 + tr_y_y_yz[i] * pa_y[i];

        tr_y_yy_zz[i] = tr_y_0_zz[i] * fe_0 + ts_y_zz[i] * fe_0 + tr_y_y_zz[i] * pa_y[i];
    }

    // Set up 60-66 components of targeted buffer : DD

    auto tr_y_yz_xx = pbuffer.data(idx_dip_dd + 60);

    auto tr_y_yz_xy = pbuffer.data(idx_dip_dd + 61);

    auto tr_y_yz_xz = pbuffer.data(idx_dip_dd + 62);

    auto tr_y_yz_yy = pbuffer.data(idx_dip_dd + 63);

    auto tr_y_yz_yz = pbuffer.data(idx_dip_dd + 64);

    auto tr_y_yz_zz = pbuffer.data(idx_dip_dd + 65);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tr_y_y_xx,  \
                             tr_y_y_xy,  \
                             tr_y_y_y,   \
                             tr_y_y_yy,  \
                             tr_y_y_yz,  \
                             tr_y_yz_xx, \
                             tr_y_yz_xy, \
                             tr_y_yz_xz, \
                             tr_y_yz_yy, \
                             tr_y_yz_yz, \
                             tr_y_yz_zz, \
                             tr_y_z_xz,  \
                             tr_y_z_zz,  \
                             ts_z_xz,    \
                             ts_z_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yz_xx[i] = tr_y_y_xx[i] * pa_z[i];

        tr_y_yz_xy[i] = tr_y_y_xy[i] * pa_z[i];

        tr_y_yz_xz[i] = ts_z_xz[i] * fe_0 + tr_y_z_xz[i] * pa_y[i];

        tr_y_yz_yy[i] = tr_y_y_yy[i] * pa_z[i];

        tr_y_yz_yz[i] = tr_y_y_y[i] * fe_0 + tr_y_y_yz[i] * pa_z[i];

        tr_y_yz_zz[i] = ts_z_zz[i] * fe_0 + tr_y_z_zz[i] * pa_y[i];
    }

    // Set up 66-72 components of targeted buffer : DD

    auto tr_y_zz_xx = pbuffer.data(idx_dip_dd + 66);

    auto tr_y_zz_xy = pbuffer.data(idx_dip_dd + 67);

    auto tr_y_zz_xz = pbuffer.data(idx_dip_dd + 68);

    auto tr_y_zz_yy = pbuffer.data(idx_dip_dd + 69);

    auto tr_y_zz_yz = pbuffer.data(idx_dip_dd + 70);

    auto tr_y_zz_zz = pbuffer.data(idx_dip_dd + 71);

#pragma omp simd aligned(pa_z,           \
                             tr_y_0_xx,  \
                             tr_y_0_xy,  \
                             tr_y_0_xz,  \
                             tr_y_0_yy,  \
                             tr_y_0_yz,  \
                             tr_y_0_zz,  \
                             tr_y_z_x,   \
                             tr_y_z_xx,  \
                             tr_y_z_xy,  \
                             tr_y_z_xz,  \
                             tr_y_z_y,   \
                             tr_y_z_yy,  \
                             tr_y_z_yz,  \
                             tr_y_z_z,   \
                             tr_y_z_zz,  \
                             tr_y_zz_xx, \
                             tr_y_zz_xy, \
                             tr_y_zz_xz, \
                             tr_y_zz_yy, \
                             tr_y_zz_yz, \
                             tr_y_zz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zz_xx[i] = tr_y_0_xx[i] * fe_0 + tr_y_z_xx[i] * pa_z[i];

        tr_y_zz_xy[i] = tr_y_0_xy[i] * fe_0 + tr_y_z_xy[i] * pa_z[i];

        tr_y_zz_xz[i] = tr_y_0_xz[i] * fe_0 + tr_y_z_x[i] * fe_0 + tr_y_z_xz[i] * pa_z[i];

        tr_y_zz_yy[i] = tr_y_0_yy[i] * fe_0 + tr_y_z_yy[i] * pa_z[i];

        tr_y_zz_yz[i] = tr_y_0_yz[i] * fe_0 + tr_y_z_y[i] * fe_0 + tr_y_z_yz[i] * pa_z[i];

        tr_y_zz_zz[i] = tr_y_0_zz[i] * fe_0 + 2.0 * tr_y_z_z[i] * fe_0 + tr_y_z_zz[i] * pa_z[i];
    }

    // Set up 72-78 components of targeted buffer : DD

    auto tr_z_xx_xx = pbuffer.data(idx_dip_dd + 72);

    auto tr_z_xx_xy = pbuffer.data(idx_dip_dd + 73);

    auto tr_z_xx_xz = pbuffer.data(idx_dip_dd + 74);

    auto tr_z_xx_yy = pbuffer.data(idx_dip_dd + 75);

    auto tr_z_xx_yz = pbuffer.data(idx_dip_dd + 76);

    auto tr_z_xx_zz = pbuffer.data(idx_dip_dd + 77);

#pragma omp simd aligned(pa_x,           \
                             tr_z_0_xx,  \
                             tr_z_0_xy,  \
                             tr_z_0_xz,  \
                             tr_z_0_yy,  \
                             tr_z_0_yz,  \
                             tr_z_0_zz,  \
                             tr_z_x_x,   \
                             tr_z_x_xx,  \
                             tr_z_x_xy,  \
                             tr_z_x_xz,  \
                             tr_z_x_y,   \
                             tr_z_x_yy,  \
                             tr_z_x_yz,  \
                             tr_z_x_z,   \
                             tr_z_x_zz,  \
                             tr_z_xx_xx, \
                             tr_z_xx_xy, \
                             tr_z_xx_xz, \
                             tr_z_xx_yy, \
                             tr_z_xx_yz, \
                             tr_z_xx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xx_xx[i] = tr_z_0_xx[i] * fe_0 + 2.0 * tr_z_x_x[i] * fe_0 + tr_z_x_xx[i] * pa_x[i];

        tr_z_xx_xy[i] = tr_z_0_xy[i] * fe_0 + tr_z_x_y[i] * fe_0 + tr_z_x_xy[i] * pa_x[i];

        tr_z_xx_xz[i] = tr_z_0_xz[i] * fe_0 + tr_z_x_z[i] * fe_0 + tr_z_x_xz[i] * pa_x[i];

        tr_z_xx_yy[i] = tr_z_0_yy[i] * fe_0 + tr_z_x_yy[i] * pa_x[i];

        tr_z_xx_yz[i] = tr_z_0_yz[i] * fe_0 + tr_z_x_yz[i] * pa_x[i];

        tr_z_xx_zz[i] = tr_z_0_zz[i] * fe_0 + tr_z_x_zz[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : DD

    auto tr_z_xy_xx = pbuffer.data(idx_dip_dd + 78);

    auto tr_z_xy_xy = pbuffer.data(idx_dip_dd + 79);

    auto tr_z_xy_xz = pbuffer.data(idx_dip_dd + 80);

    auto tr_z_xy_yy = pbuffer.data(idx_dip_dd + 81);

    auto tr_z_xy_yz = pbuffer.data(idx_dip_dd + 82);

    auto tr_z_xy_zz = pbuffer.data(idx_dip_dd + 83);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tr_z_x_xx,  \
                             tr_z_x_xz,  \
                             tr_z_xy_xx, \
                             tr_z_xy_xy, \
                             tr_z_xy_xz, \
                             tr_z_xy_yy, \
                             tr_z_xy_yz, \
                             tr_z_xy_zz, \
                             tr_z_y_xy,  \
                             tr_z_y_y,   \
                             tr_z_y_yy,  \
                             tr_z_y_yz,  \
                             tr_z_y_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xy_xx[i] = tr_z_x_xx[i] * pa_y[i];

        tr_z_xy_xy[i] = tr_z_y_y[i] * fe_0 + tr_z_y_xy[i] * pa_x[i];

        tr_z_xy_xz[i] = tr_z_x_xz[i] * pa_y[i];

        tr_z_xy_yy[i] = tr_z_y_yy[i] * pa_x[i];

        tr_z_xy_yz[i] = tr_z_y_yz[i] * pa_x[i];

        tr_z_xy_zz[i] = tr_z_y_zz[i] * pa_x[i];
    }

    // Set up 84-90 components of targeted buffer : DD

    auto tr_z_xz_xx = pbuffer.data(idx_dip_dd + 84);

    auto tr_z_xz_xy = pbuffer.data(idx_dip_dd + 85);

    auto tr_z_xz_xz = pbuffer.data(idx_dip_dd + 86);

    auto tr_z_xz_yy = pbuffer.data(idx_dip_dd + 87);

    auto tr_z_xz_yz = pbuffer.data(idx_dip_dd + 88);

    auto tr_z_xz_zz = pbuffer.data(idx_dip_dd + 89);

#pragma omp simd aligned(pa_x,           \
                             tr_z_xz_xx, \
                             tr_z_xz_xy, \
                             tr_z_xz_xz, \
                             tr_z_xz_yy, \
                             tr_z_xz_yz, \
                             tr_z_xz_zz, \
                             tr_z_z_x,   \
                             tr_z_z_xx,  \
                             tr_z_z_xy,  \
                             tr_z_z_xz,  \
                             tr_z_z_y,   \
                             tr_z_z_yy,  \
                             tr_z_z_yz,  \
                             tr_z_z_z,   \
                             tr_z_z_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xz_xx[i] = 2.0 * tr_z_z_x[i] * fe_0 + tr_z_z_xx[i] * pa_x[i];

        tr_z_xz_xy[i] = tr_z_z_y[i] * fe_0 + tr_z_z_xy[i] * pa_x[i];

        tr_z_xz_xz[i] = tr_z_z_z[i] * fe_0 + tr_z_z_xz[i] * pa_x[i];

        tr_z_xz_yy[i] = tr_z_z_yy[i] * pa_x[i];

        tr_z_xz_yz[i] = tr_z_z_yz[i] * pa_x[i];

        tr_z_xz_zz[i] = tr_z_z_zz[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : DD

    auto tr_z_yy_xx = pbuffer.data(idx_dip_dd + 90);

    auto tr_z_yy_xy = pbuffer.data(idx_dip_dd + 91);

    auto tr_z_yy_xz = pbuffer.data(idx_dip_dd + 92);

    auto tr_z_yy_yy = pbuffer.data(idx_dip_dd + 93);

    auto tr_z_yy_yz = pbuffer.data(idx_dip_dd + 94);

    auto tr_z_yy_zz = pbuffer.data(idx_dip_dd + 95);

#pragma omp simd aligned(pa_y,           \
                             tr_z_0_xx,  \
                             tr_z_0_xy,  \
                             tr_z_0_xz,  \
                             tr_z_0_yy,  \
                             tr_z_0_yz,  \
                             tr_z_0_zz,  \
                             tr_z_y_x,   \
                             tr_z_y_xx,  \
                             tr_z_y_xy,  \
                             tr_z_y_xz,  \
                             tr_z_y_y,   \
                             tr_z_y_yy,  \
                             tr_z_y_yz,  \
                             tr_z_y_z,   \
                             tr_z_y_zz,  \
                             tr_z_yy_xx, \
                             tr_z_yy_xy, \
                             tr_z_yy_xz, \
                             tr_z_yy_yy, \
                             tr_z_yy_yz, \
                             tr_z_yy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yy_xx[i] = tr_z_0_xx[i] * fe_0 + tr_z_y_xx[i] * pa_y[i];

        tr_z_yy_xy[i] = tr_z_0_xy[i] * fe_0 + tr_z_y_x[i] * fe_0 + tr_z_y_xy[i] * pa_y[i];

        tr_z_yy_xz[i] = tr_z_0_xz[i] * fe_0 + tr_z_y_xz[i] * pa_y[i];

        tr_z_yy_yy[i] = tr_z_0_yy[i] * fe_0 + 2.0 * tr_z_y_y[i] * fe_0 + tr_z_y_yy[i] * pa_y[i];

        tr_z_yy_yz[i] = tr_z_0_yz[i] * fe_0 + tr_z_y_z[i] * fe_0 + tr_z_y_yz[i] * pa_y[i];

        tr_z_yy_zz[i] = tr_z_0_zz[i] * fe_0 + tr_z_y_zz[i] * pa_y[i];
    }

    // Set up 96-102 components of targeted buffer : DD

    auto tr_z_yz_xx = pbuffer.data(idx_dip_dd + 96);

    auto tr_z_yz_xy = pbuffer.data(idx_dip_dd + 97);

    auto tr_z_yz_xz = pbuffer.data(idx_dip_dd + 98);

    auto tr_z_yz_yy = pbuffer.data(idx_dip_dd + 99);

    auto tr_z_yz_yz = pbuffer.data(idx_dip_dd + 100);

    auto tr_z_yz_zz = pbuffer.data(idx_dip_dd + 101);

#pragma omp simd aligned(pa_y,           \
                             tr_z_yz_xx, \
                             tr_z_yz_xy, \
                             tr_z_yz_xz, \
                             tr_z_yz_yy, \
                             tr_z_yz_yz, \
                             tr_z_yz_zz, \
                             tr_z_z_x,   \
                             tr_z_z_xx,  \
                             tr_z_z_xy,  \
                             tr_z_z_xz,  \
                             tr_z_z_y,   \
                             tr_z_z_yy,  \
                             tr_z_z_yz,  \
                             tr_z_z_z,   \
                             tr_z_z_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yz_xx[i] = tr_z_z_xx[i] * pa_y[i];

        tr_z_yz_xy[i] = tr_z_z_x[i] * fe_0 + tr_z_z_xy[i] * pa_y[i];

        tr_z_yz_xz[i] = tr_z_z_xz[i] * pa_y[i];

        tr_z_yz_yy[i] = 2.0 * tr_z_z_y[i] * fe_0 + tr_z_z_yy[i] * pa_y[i];

        tr_z_yz_yz[i] = tr_z_z_z[i] * fe_0 + tr_z_z_yz[i] * pa_y[i];

        tr_z_yz_zz[i] = tr_z_z_zz[i] * pa_y[i];
    }

    // Set up 102-108 components of targeted buffer : DD

    auto tr_z_zz_xx = pbuffer.data(idx_dip_dd + 102);

    auto tr_z_zz_xy = pbuffer.data(idx_dip_dd + 103);

    auto tr_z_zz_xz = pbuffer.data(idx_dip_dd + 104);

    auto tr_z_zz_yy = pbuffer.data(idx_dip_dd + 105);

    auto tr_z_zz_yz = pbuffer.data(idx_dip_dd + 106);

    auto tr_z_zz_zz = pbuffer.data(idx_dip_dd + 107);

#pragma omp simd aligned(pa_z,           \
                             tr_z_0_xx,  \
                             tr_z_0_xy,  \
                             tr_z_0_xz,  \
                             tr_z_0_yy,  \
                             tr_z_0_yz,  \
                             tr_z_0_zz,  \
                             tr_z_z_x,   \
                             tr_z_z_xx,  \
                             tr_z_z_xy,  \
                             tr_z_z_xz,  \
                             tr_z_z_y,   \
                             tr_z_z_yy,  \
                             tr_z_z_yz,  \
                             tr_z_z_z,   \
                             tr_z_z_zz,  \
                             tr_z_zz_xx, \
                             tr_z_zz_xy, \
                             tr_z_zz_xz, \
                             tr_z_zz_yy, \
                             tr_z_zz_yz, \
                             tr_z_zz_zz, \
                             ts_z_xx,    \
                             ts_z_xy,    \
                             ts_z_xz,    \
                             ts_z_yy,    \
                             ts_z_yz,    \
                             ts_z_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zz_xx[i] = tr_z_0_xx[i] * fe_0 + ts_z_xx[i] * fe_0 + tr_z_z_xx[i] * pa_z[i];

        tr_z_zz_xy[i] = tr_z_0_xy[i] * fe_0 + ts_z_xy[i] * fe_0 + tr_z_z_xy[i] * pa_z[i];

        tr_z_zz_xz[i] = tr_z_0_xz[i] * fe_0 + tr_z_z_x[i] * fe_0 + ts_z_xz[i] * fe_0 + tr_z_z_xz[i] * pa_z[i];

        tr_z_zz_yy[i] = tr_z_0_yy[i] * fe_0 + ts_z_yy[i] * fe_0 + tr_z_z_yy[i] * pa_z[i];

        tr_z_zz_yz[i] = tr_z_0_yz[i] * fe_0 + tr_z_z_y[i] * fe_0 + ts_z_yz[i] * fe_0 + tr_z_z_yz[i] * pa_z[i];

        tr_z_zz_zz[i] = tr_z_0_zz[i] * fe_0 + 2.0 * tr_z_z_z[i] * fe_0 + ts_z_zz[i] * fe_0 + tr_z_z_zz[i] * pa_z[i];
    }
}

}  // namespace diprec
