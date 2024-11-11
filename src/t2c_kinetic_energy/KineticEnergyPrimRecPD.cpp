#include "KineticEnergyPrimRecPD.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_pd(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_pd,
                            const size_t              idx_kin_sp,
                            const size_t              idx_kin_sd,
                            const size_t              idx_ovl_pd,
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

    auto tk_0_x = pbuffer.data(idx_kin_sp);

    auto tk_0_y = pbuffer.data(idx_kin_sp + 1);

    auto tk_0_z = pbuffer.data(idx_kin_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto tk_0_xx = pbuffer.data(idx_kin_sd);

    auto tk_0_xy = pbuffer.data(idx_kin_sd + 1);

    auto tk_0_xz = pbuffer.data(idx_kin_sd + 2);

    auto tk_0_yy = pbuffer.data(idx_kin_sd + 3);

    auto tk_0_yz = pbuffer.data(idx_kin_sd + 4);

    auto tk_0_zz = pbuffer.data(idx_kin_sd + 5);

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

    // Set up 0-6 components of targeted buffer : PD

    auto tk_x_xx = pbuffer.data(idx_kin_pd);

    auto tk_x_xy = pbuffer.data(idx_kin_pd + 1);

    auto tk_x_xz = pbuffer.data(idx_kin_pd + 2);

    auto tk_x_yy = pbuffer.data(idx_kin_pd + 3);

    auto tk_x_yz = pbuffer.data(idx_kin_pd + 4);

    auto tk_x_zz = pbuffer.data(idx_kin_pd + 5);

#pragma omp simd aligned(pa_x,        \
                             tk_0_x,  \
                             tk_0_xx, \
                             tk_0_xy, \
                             tk_0_xz, \
                             tk_0_y,  \
                             tk_0_yy, \
                             tk_0_yz, \
                             tk_0_z,  \
                             tk_0_zz, \
                             tk_x_xx, \
                             tk_x_xy, \
                             tk_x_xz, \
                             tk_x_yy, \
                             tk_x_yz, \
                             tk_x_zz, \
                             ts_x_xx, \
                             ts_x_xy, \
                             ts_x_xz, \
                             ts_x_yy, \
                             ts_x_yz, \
                             ts_x_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_x_xx[i] = 2.0 * tk_0_x[i] * fe_0 + tk_0_xx[i] * pa_x[i] + 2.0 * ts_x_xx[i] * fz_0;

        tk_x_xy[i] = tk_0_y[i] * fe_0 + tk_0_xy[i] * pa_x[i] + 2.0 * ts_x_xy[i] * fz_0;

        tk_x_xz[i] = tk_0_z[i] * fe_0 + tk_0_xz[i] * pa_x[i] + 2.0 * ts_x_xz[i] * fz_0;

        tk_x_yy[i] = tk_0_yy[i] * pa_x[i] + 2.0 * ts_x_yy[i] * fz_0;

        tk_x_yz[i] = tk_0_yz[i] * pa_x[i] + 2.0 * ts_x_yz[i] * fz_0;

        tk_x_zz[i] = tk_0_zz[i] * pa_x[i] + 2.0 * ts_x_zz[i] * fz_0;
    }

    // Set up 6-12 components of targeted buffer : PD

    auto tk_y_xx = pbuffer.data(idx_kin_pd + 6);

    auto tk_y_xy = pbuffer.data(idx_kin_pd + 7);

    auto tk_y_xz = pbuffer.data(idx_kin_pd + 8);

    auto tk_y_yy = pbuffer.data(idx_kin_pd + 9);

    auto tk_y_yz = pbuffer.data(idx_kin_pd + 10);

    auto tk_y_zz = pbuffer.data(idx_kin_pd + 11);

#pragma omp simd aligned(pa_y,        \
                             tk_0_x,  \
                             tk_0_xx, \
                             tk_0_xy, \
                             tk_0_xz, \
                             tk_0_y,  \
                             tk_0_yy, \
                             tk_0_yz, \
                             tk_0_z,  \
                             tk_0_zz, \
                             tk_y_xx, \
                             tk_y_xy, \
                             tk_y_xz, \
                             tk_y_yy, \
                             tk_y_yz, \
                             tk_y_zz, \
                             ts_y_xx, \
                             ts_y_xy, \
                             ts_y_xz, \
                             ts_y_yy, \
                             ts_y_yz, \
                             ts_y_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_y_xx[i] = tk_0_xx[i] * pa_y[i] + 2.0 * ts_y_xx[i] * fz_0;

        tk_y_xy[i] = tk_0_x[i] * fe_0 + tk_0_xy[i] * pa_y[i] + 2.0 * ts_y_xy[i] * fz_0;

        tk_y_xz[i] = tk_0_xz[i] * pa_y[i] + 2.0 * ts_y_xz[i] * fz_0;

        tk_y_yy[i] = 2.0 * tk_0_y[i] * fe_0 + tk_0_yy[i] * pa_y[i] + 2.0 * ts_y_yy[i] * fz_0;

        tk_y_yz[i] = tk_0_z[i] * fe_0 + tk_0_yz[i] * pa_y[i] + 2.0 * ts_y_yz[i] * fz_0;

        tk_y_zz[i] = tk_0_zz[i] * pa_y[i] + 2.0 * ts_y_zz[i] * fz_0;
    }

    // Set up 12-18 components of targeted buffer : PD

    auto tk_z_xx = pbuffer.data(idx_kin_pd + 12);

    auto tk_z_xy = pbuffer.data(idx_kin_pd + 13);

    auto tk_z_xz = pbuffer.data(idx_kin_pd + 14);

    auto tk_z_yy = pbuffer.data(idx_kin_pd + 15);

    auto tk_z_yz = pbuffer.data(idx_kin_pd + 16);

    auto tk_z_zz = pbuffer.data(idx_kin_pd + 17);

#pragma omp simd aligned(pa_z,        \
                             tk_0_x,  \
                             tk_0_xx, \
                             tk_0_xy, \
                             tk_0_xz, \
                             tk_0_y,  \
                             tk_0_yy, \
                             tk_0_yz, \
                             tk_0_z,  \
                             tk_0_zz, \
                             tk_z_xx, \
                             tk_z_xy, \
                             tk_z_xz, \
                             tk_z_yy, \
                             tk_z_yz, \
                             tk_z_zz, \
                             ts_z_xx, \
                             ts_z_xy, \
                             ts_z_xz, \
                             ts_z_yy, \
                             ts_z_yz, \
                             ts_z_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_z_xx[i] = tk_0_xx[i] * pa_z[i] + 2.0 * ts_z_xx[i] * fz_0;

        tk_z_xy[i] = tk_0_xy[i] * pa_z[i] + 2.0 * ts_z_xy[i] * fz_0;

        tk_z_xz[i] = tk_0_x[i] * fe_0 + tk_0_xz[i] * pa_z[i] + 2.0 * ts_z_xz[i] * fz_0;

        tk_z_yy[i] = tk_0_yy[i] * pa_z[i] + 2.0 * ts_z_yy[i] * fz_0;

        tk_z_yz[i] = tk_0_y[i] * fe_0 + tk_0_yz[i] * pa_z[i] + 2.0 * ts_z_yz[i] * fz_0;

        tk_z_zz[i] = 2.0 * tk_0_z[i] * fe_0 + tk_0_zz[i] * pa_z[i] + 2.0 * ts_z_zz[i] * fz_0;
    }
}

}  // namespace kinrec
