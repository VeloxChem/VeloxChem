#include "ElectricDipoleMomentumPrimRecPD.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_pd(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_pd,
                                      const size_t              idx_dip_sp,
                                      const size_t              idx_ovl_sd,
                                      const size_t              idx_dip_sd,
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

    // Set up components of auxiliary buffer : SP

    auto tr_x_0_x = pbuffer.data(idx_dip_sp);

    auto tr_x_0_y = pbuffer.data(idx_dip_sp + 1);

    auto tr_x_0_z = pbuffer.data(idx_dip_sp + 2);

    auto tr_y_0_x = pbuffer.data(idx_dip_sp + 3);

    auto tr_y_0_y = pbuffer.data(idx_dip_sp + 4);

    auto tr_y_0_z = pbuffer.data(idx_dip_sp + 5);

    auto tr_z_0_x = pbuffer.data(idx_dip_sp + 6);

    auto tr_z_0_y = pbuffer.data(idx_dip_sp + 7);

    auto tr_z_0_z = pbuffer.data(idx_dip_sp + 8);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_xy = pbuffer.data(idx_ovl_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_ovl_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_ovl_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

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

    // Set up 0-6 components of targeted buffer : PD

    auto tr_x_x_xx = pbuffer.data(idx_dip_pd);

    auto tr_x_x_xy = pbuffer.data(idx_dip_pd + 1);

    auto tr_x_x_xz = pbuffer.data(idx_dip_pd + 2);

    auto tr_x_x_yy = pbuffer.data(idx_dip_pd + 3);

    auto tr_x_x_yz = pbuffer.data(idx_dip_pd + 4);

    auto tr_x_x_zz = pbuffer.data(idx_dip_pd + 5);

#pragma omp simd aligned(pa_x,          \
                             tr_x_0_x,  \
                             tr_x_0_xx, \
                             tr_x_0_xy, \
                             tr_x_0_xz, \
                             tr_x_0_y,  \
                             tr_x_0_yy, \
                             tr_x_0_yz, \
                             tr_x_0_z,  \
                             tr_x_0_zz, \
                             tr_x_x_xx, \
                             tr_x_x_xy, \
                             tr_x_x_xz, \
                             tr_x_x_yy, \
                             tr_x_x_yz, \
                             tr_x_x_zz, \
                             ts_0_xx,   \
                             ts_0_xy,   \
                             ts_0_xz,   \
                             ts_0_yy,   \
                             ts_0_yz,   \
                             ts_0_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_x_xx[i] = 2.0 * tr_x_0_x[i] * fe_0 + ts_0_xx[i] * fe_0 + tr_x_0_xx[i] * pa_x[i];

        tr_x_x_xy[i] = tr_x_0_y[i] * fe_0 + ts_0_xy[i] * fe_0 + tr_x_0_xy[i] * pa_x[i];

        tr_x_x_xz[i] = tr_x_0_z[i] * fe_0 + ts_0_xz[i] * fe_0 + tr_x_0_xz[i] * pa_x[i];

        tr_x_x_yy[i] = ts_0_yy[i] * fe_0 + tr_x_0_yy[i] * pa_x[i];

        tr_x_x_yz[i] = ts_0_yz[i] * fe_0 + tr_x_0_yz[i] * pa_x[i];

        tr_x_x_zz[i] = ts_0_zz[i] * fe_0 + tr_x_0_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : PD

    auto tr_x_y_xx = pbuffer.data(idx_dip_pd + 6);

    auto tr_x_y_xy = pbuffer.data(idx_dip_pd + 7);

    auto tr_x_y_xz = pbuffer.data(idx_dip_pd + 8);

    auto tr_x_y_yy = pbuffer.data(idx_dip_pd + 9);

    auto tr_x_y_yz = pbuffer.data(idx_dip_pd + 10);

    auto tr_x_y_zz = pbuffer.data(idx_dip_pd + 11);

#pragma omp simd aligned(pa_y,          \
                             tr_x_0_x,  \
                             tr_x_0_xx, \
                             tr_x_0_xy, \
                             tr_x_0_xz, \
                             tr_x_0_y,  \
                             tr_x_0_yy, \
                             tr_x_0_yz, \
                             tr_x_0_z,  \
                             tr_x_0_zz, \
                             tr_x_y_xx, \
                             tr_x_y_xy, \
                             tr_x_y_xz, \
                             tr_x_y_yy, \
                             tr_x_y_yz, \
                             tr_x_y_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_y_xx[i] = tr_x_0_xx[i] * pa_y[i];

        tr_x_y_xy[i] = tr_x_0_x[i] * fe_0 + tr_x_0_xy[i] * pa_y[i];

        tr_x_y_xz[i] = tr_x_0_xz[i] * pa_y[i];

        tr_x_y_yy[i] = 2.0 * tr_x_0_y[i] * fe_0 + tr_x_0_yy[i] * pa_y[i];

        tr_x_y_yz[i] = tr_x_0_z[i] * fe_0 + tr_x_0_yz[i] * pa_y[i];

        tr_x_y_zz[i] = tr_x_0_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : PD

    auto tr_x_z_xx = pbuffer.data(idx_dip_pd + 12);

    auto tr_x_z_xy = pbuffer.data(idx_dip_pd + 13);

    auto tr_x_z_xz = pbuffer.data(idx_dip_pd + 14);

    auto tr_x_z_yy = pbuffer.data(idx_dip_pd + 15);

    auto tr_x_z_yz = pbuffer.data(idx_dip_pd + 16);

    auto tr_x_z_zz = pbuffer.data(idx_dip_pd + 17);

#pragma omp simd aligned(pa_z,          \
                             tr_x_0_x,  \
                             tr_x_0_xx, \
                             tr_x_0_xy, \
                             tr_x_0_xz, \
                             tr_x_0_y,  \
                             tr_x_0_yy, \
                             tr_x_0_yz, \
                             tr_x_0_z,  \
                             tr_x_0_zz, \
                             tr_x_z_xx, \
                             tr_x_z_xy, \
                             tr_x_z_xz, \
                             tr_x_z_yy, \
                             tr_x_z_yz, \
                             tr_x_z_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_z_xx[i] = tr_x_0_xx[i] * pa_z[i];

        tr_x_z_xy[i] = tr_x_0_xy[i] * pa_z[i];

        tr_x_z_xz[i] = tr_x_0_x[i] * fe_0 + tr_x_0_xz[i] * pa_z[i];

        tr_x_z_yy[i] = tr_x_0_yy[i] * pa_z[i];

        tr_x_z_yz[i] = tr_x_0_y[i] * fe_0 + tr_x_0_yz[i] * pa_z[i];

        tr_x_z_zz[i] = 2.0 * tr_x_0_z[i] * fe_0 + tr_x_0_zz[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : PD

    auto tr_y_x_xx = pbuffer.data(idx_dip_pd + 18);

    auto tr_y_x_xy = pbuffer.data(idx_dip_pd + 19);

    auto tr_y_x_xz = pbuffer.data(idx_dip_pd + 20);

    auto tr_y_x_yy = pbuffer.data(idx_dip_pd + 21);

    auto tr_y_x_yz = pbuffer.data(idx_dip_pd + 22);

    auto tr_y_x_zz = pbuffer.data(idx_dip_pd + 23);

#pragma omp simd aligned(pa_x,          \
                             tr_y_0_x,  \
                             tr_y_0_xx, \
                             tr_y_0_xy, \
                             tr_y_0_xz, \
                             tr_y_0_y,  \
                             tr_y_0_yy, \
                             tr_y_0_yz, \
                             tr_y_0_z,  \
                             tr_y_0_zz, \
                             tr_y_x_xx, \
                             tr_y_x_xy, \
                             tr_y_x_xz, \
                             tr_y_x_yy, \
                             tr_y_x_yz, \
                             tr_y_x_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_x_xx[i] = 2.0 * tr_y_0_x[i] * fe_0 + tr_y_0_xx[i] * pa_x[i];

        tr_y_x_xy[i] = tr_y_0_y[i] * fe_0 + tr_y_0_xy[i] * pa_x[i];

        tr_y_x_xz[i] = tr_y_0_z[i] * fe_0 + tr_y_0_xz[i] * pa_x[i];

        tr_y_x_yy[i] = tr_y_0_yy[i] * pa_x[i];

        tr_y_x_yz[i] = tr_y_0_yz[i] * pa_x[i];

        tr_y_x_zz[i] = tr_y_0_zz[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : PD

    auto tr_y_y_xx = pbuffer.data(idx_dip_pd + 24);

    auto tr_y_y_xy = pbuffer.data(idx_dip_pd + 25);

    auto tr_y_y_xz = pbuffer.data(idx_dip_pd + 26);

    auto tr_y_y_yy = pbuffer.data(idx_dip_pd + 27);

    auto tr_y_y_yz = pbuffer.data(idx_dip_pd + 28);

    auto tr_y_y_zz = pbuffer.data(idx_dip_pd + 29);

#pragma omp simd aligned(pa_y,          \
                             tr_y_0_x,  \
                             tr_y_0_xx, \
                             tr_y_0_xy, \
                             tr_y_0_xz, \
                             tr_y_0_y,  \
                             tr_y_0_yy, \
                             tr_y_0_yz, \
                             tr_y_0_z,  \
                             tr_y_0_zz, \
                             tr_y_y_xx, \
                             tr_y_y_xy, \
                             tr_y_y_xz, \
                             tr_y_y_yy, \
                             tr_y_y_yz, \
                             tr_y_y_zz, \
                             ts_0_xx,   \
                             ts_0_xy,   \
                             ts_0_xz,   \
                             ts_0_yy,   \
                             ts_0_yz,   \
                             ts_0_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_y_xx[i] = ts_0_xx[i] * fe_0 + tr_y_0_xx[i] * pa_y[i];

        tr_y_y_xy[i] = tr_y_0_x[i] * fe_0 + ts_0_xy[i] * fe_0 + tr_y_0_xy[i] * pa_y[i];

        tr_y_y_xz[i] = ts_0_xz[i] * fe_0 + tr_y_0_xz[i] * pa_y[i];

        tr_y_y_yy[i] = 2.0 * tr_y_0_y[i] * fe_0 + ts_0_yy[i] * fe_0 + tr_y_0_yy[i] * pa_y[i];

        tr_y_y_yz[i] = tr_y_0_z[i] * fe_0 + ts_0_yz[i] * fe_0 + tr_y_0_yz[i] * pa_y[i];

        tr_y_y_zz[i] = ts_0_zz[i] * fe_0 + tr_y_0_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : PD

    auto tr_y_z_xx = pbuffer.data(idx_dip_pd + 30);

    auto tr_y_z_xy = pbuffer.data(idx_dip_pd + 31);

    auto tr_y_z_xz = pbuffer.data(idx_dip_pd + 32);

    auto tr_y_z_yy = pbuffer.data(idx_dip_pd + 33);

    auto tr_y_z_yz = pbuffer.data(idx_dip_pd + 34);

    auto tr_y_z_zz = pbuffer.data(idx_dip_pd + 35);

#pragma omp simd aligned(pa_z,          \
                             tr_y_0_x,  \
                             tr_y_0_xx, \
                             tr_y_0_xy, \
                             tr_y_0_xz, \
                             tr_y_0_y,  \
                             tr_y_0_yy, \
                             tr_y_0_yz, \
                             tr_y_0_z,  \
                             tr_y_0_zz, \
                             tr_y_z_xx, \
                             tr_y_z_xy, \
                             tr_y_z_xz, \
                             tr_y_z_yy, \
                             tr_y_z_yz, \
                             tr_y_z_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_z_xx[i] = tr_y_0_xx[i] * pa_z[i];

        tr_y_z_xy[i] = tr_y_0_xy[i] * pa_z[i];

        tr_y_z_xz[i] = tr_y_0_x[i] * fe_0 + tr_y_0_xz[i] * pa_z[i];

        tr_y_z_yy[i] = tr_y_0_yy[i] * pa_z[i];

        tr_y_z_yz[i] = tr_y_0_y[i] * fe_0 + tr_y_0_yz[i] * pa_z[i];

        tr_y_z_zz[i] = 2.0 * tr_y_0_z[i] * fe_0 + tr_y_0_zz[i] * pa_z[i];
    }

    // Set up 36-42 components of targeted buffer : PD

    auto tr_z_x_xx = pbuffer.data(idx_dip_pd + 36);

    auto tr_z_x_xy = pbuffer.data(idx_dip_pd + 37);

    auto tr_z_x_xz = pbuffer.data(idx_dip_pd + 38);

    auto tr_z_x_yy = pbuffer.data(idx_dip_pd + 39);

    auto tr_z_x_yz = pbuffer.data(idx_dip_pd + 40);

    auto tr_z_x_zz = pbuffer.data(idx_dip_pd + 41);

#pragma omp simd aligned(pa_x,          \
                             tr_z_0_x,  \
                             tr_z_0_xx, \
                             tr_z_0_xy, \
                             tr_z_0_xz, \
                             tr_z_0_y,  \
                             tr_z_0_yy, \
                             tr_z_0_yz, \
                             tr_z_0_z,  \
                             tr_z_0_zz, \
                             tr_z_x_xx, \
                             tr_z_x_xy, \
                             tr_z_x_xz, \
                             tr_z_x_yy, \
                             tr_z_x_yz, \
                             tr_z_x_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_x_xx[i] = 2.0 * tr_z_0_x[i] * fe_0 + tr_z_0_xx[i] * pa_x[i];

        tr_z_x_xy[i] = tr_z_0_y[i] * fe_0 + tr_z_0_xy[i] * pa_x[i];

        tr_z_x_xz[i] = tr_z_0_z[i] * fe_0 + tr_z_0_xz[i] * pa_x[i];

        tr_z_x_yy[i] = tr_z_0_yy[i] * pa_x[i];

        tr_z_x_yz[i] = tr_z_0_yz[i] * pa_x[i];

        tr_z_x_zz[i] = tr_z_0_zz[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : PD

    auto tr_z_y_xx = pbuffer.data(idx_dip_pd + 42);

    auto tr_z_y_xy = pbuffer.data(idx_dip_pd + 43);

    auto tr_z_y_xz = pbuffer.data(idx_dip_pd + 44);

    auto tr_z_y_yy = pbuffer.data(idx_dip_pd + 45);

    auto tr_z_y_yz = pbuffer.data(idx_dip_pd + 46);

    auto tr_z_y_zz = pbuffer.data(idx_dip_pd + 47);

#pragma omp simd aligned(pa_y,          \
                             tr_z_0_x,  \
                             tr_z_0_xx, \
                             tr_z_0_xy, \
                             tr_z_0_xz, \
                             tr_z_0_y,  \
                             tr_z_0_yy, \
                             tr_z_0_yz, \
                             tr_z_0_z,  \
                             tr_z_0_zz, \
                             tr_z_y_xx, \
                             tr_z_y_xy, \
                             tr_z_y_xz, \
                             tr_z_y_yy, \
                             tr_z_y_yz, \
                             tr_z_y_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_y_xx[i] = tr_z_0_xx[i] * pa_y[i];

        tr_z_y_xy[i] = tr_z_0_x[i] * fe_0 + tr_z_0_xy[i] * pa_y[i];

        tr_z_y_xz[i] = tr_z_0_xz[i] * pa_y[i];

        tr_z_y_yy[i] = 2.0 * tr_z_0_y[i] * fe_0 + tr_z_0_yy[i] * pa_y[i];

        tr_z_y_yz[i] = tr_z_0_z[i] * fe_0 + tr_z_0_yz[i] * pa_y[i];

        tr_z_y_zz[i] = tr_z_0_zz[i] * pa_y[i];
    }

    // Set up 48-54 components of targeted buffer : PD

    auto tr_z_z_xx = pbuffer.data(idx_dip_pd + 48);

    auto tr_z_z_xy = pbuffer.data(idx_dip_pd + 49);

    auto tr_z_z_xz = pbuffer.data(idx_dip_pd + 50);

    auto tr_z_z_yy = pbuffer.data(idx_dip_pd + 51);

    auto tr_z_z_yz = pbuffer.data(idx_dip_pd + 52);

    auto tr_z_z_zz = pbuffer.data(idx_dip_pd + 53);

#pragma omp simd aligned(pa_z,          \
                             tr_z_0_x,  \
                             tr_z_0_xx, \
                             tr_z_0_xy, \
                             tr_z_0_xz, \
                             tr_z_0_y,  \
                             tr_z_0_yy, \
                             tr_z_0_yz, \
                             tr_z_0_z,  \
                             tr_z_0_zz, \
                             tr_z_z_xx, \
                             tr_z_z_xy, \
                             tr_z_z_xz, \
                             tr_z_z_yy, \
                             tr_z_z_yz, \
                             tr_z_z_zz, \
                             ts_0_xx,   \
                             ts_0_xy,   \
                             ts_0_xz,   \
                             ts_0_yy,   \
                             ts_0_yz,   \
                             ts_0_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_z_xx[i] = ts_0_xx[i] * fe_0 + tr_z_0_xx[i] * pa_z[i];

        tr_z_z_xy[i] = ts_0_xy[i] * fe_0 + tr_z_0_xy[i] * pa_z[i];

        tr_z_z_xz[i] = tr_z_0_x[i] * fe_0 + ts_0_xz[i] * fe_0 + tr_z_0_xz[i] * pa_z[i];

        tr_z_z_yy[i] = ts_0_yy[i] * fe_0 + tr_z_0_yy[i] * pa_z[i];

        tr_z_z_yz[i] = tr_z_0_y[i] * fe_0 + ts_0_yz[i] * fe_0 + tr_z_0_yz[i] * pa_z[i];

        tr_z_z_zz[i] = 2.0 * tr_z_0_z[i] * fe_0 + ts_0_zz[i] * fe_0 + tr_z_0_zz[i] * pa_z[i];
    }
}

}  // namespace diprec
