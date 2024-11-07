#include "KineticEnergyPrimRecPF.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_pf(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_pf,
                            const size_t              idx_kin_sd,
                            const size_t              idx_kin_sf,
                            const size_t              idx_ovl_pf,
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

    auto tk_0_xx = pbuffer.data(idx_kin_sd);

    auto tk_0_xy = pbuffer.data(idx_kin_sd + 1);

    auto tk_0_xz = pbuffer.data(idx_kin_sd + 2);

    auto tk_0_yy = pbuffer.data(idx_kin_sd + 3);

    auto tk_0_yz = pbuffer.data(idx_kin_sd + 4);

    auto tk_0_zz = pbuffer.data(idx_kin_sd + 5);

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

    // Set up components of auxiliary buffer : PF

    auto ts_x_xxx = pbuffer.data(idx_ovl_pf);

    auto ts_x_xxy = pbuffer.data(idx_ovl_pf + 1);

    auto ts_x_xxz = pbuffer.data(idx_ovl_pf + 2);

    auto ts_x_xyy = pbuffer.data(idx_ovl_pf + 3);

    auto ts_x_xyz = pbuffer.data(idx_ovl_pf + 4);

    auto ts_x_xzz = pbuffer.data(idx_ovl_pf + 5);

    auto ts_x_yyy = pbuffer.data(idx_ovl_pf + 6);

    auto ts_x_yyz = pbuffer.data(idx_ovl_pf + 7);

    auto ts_x_yzz = pbuffer.data(idx_ovl_pf + 8);

    auto ts_x_zzz = pbuffer.data(idx_ovl_pf + 9);

    auto ts_y_xxx = pbuffer.data(idx_ovl_pf + 10);

    auto ts_y_xxy = pbuffer.data(idx_ovl_pf + 11);

    auto ts_y_xxz = pbuffer.data(idx_ovl_pf + 12);

    auto ts_y_xyy = pbuffer.data(idx_ovl_pf + 13);

    auto ts_y_xyz = pbuffer.data(idx_ovl_pf + 14);

    auto ts_y_xzz = pbuffer.data(idx_ovl_pf + 15);

    auto ts_y_yyy = pbuffer.data(idx_ovl_pf + 16);

    auto ts_y_yyz = pbuffer.data(idx_ovl_pf + 17);

    auto ts_y_yzz = pbuffer.data(idx_ovl_pf + 18);

    auto ts_y_zzz = pbuffer.data(idx_ovl_pf + 19);

    auto ts_z_xxx = pbuffer.data(idx_ovl_pf + 20);

    auto ts_z_xxy = pbuffer.data(idx_ovl_pf + 21);

    auto ts_z_xxz = pbuffer.data(idx_ovl_pf + 22);

    auto ts_z_xyy = pbuffer.data(idx_ovl_pf + 23);

    auto ts_z_xyz = pbuffer.data(idx_ovl_pf + 24);

    auto ts_z_xzz = pbuffer.data(idx_ovl_pf + 25);

    auto ts_z_yyy = pbuffer.data(idx_ovl_pf + 26);

    auto ts_z_yyz = pbuffer.data(idx_ovl_pf + 27);

    auto ts_z_yzz = pbuffer.data(idx_ovl_pf + 28);

    auto ts_z_zzz = pbuffer.data(idx_ovl_pf + 29);

    // Set up 0-10 components of targeted buffer : PF

    auto tk_x_xxx = pbuffer.data(idx_kin_pf);

    auto tk_x_xxy = pbuffer.data(idx_kin_pf + 1);

    auto tk_x_xxz = pbuffer.data(idx_kin_pf + 2);

    auto tk_x_xyy = pbuffer.data(idx_kin_pf + 3);

    auto tk_x_xyz = pbuffer.data(idx_kin_pf + 4);

    auto tk_x_xzz = pbuffer.data(idx_kin_pf + 5);

    auto tk_x_yyy = pbuffer.data(idx_kin_pf + 6);

    auto tk_x_yyz = pbuffer.data(idx_kin_pf + 7);

    auto tk_x_yzz = pbuffer.data(idx_kin_pf + 8);

    auto tk_x_zzz = pbuffer.data(idx_kin_pf + 9);

#pragma omp simd aligned(pa_x,         \
                             tk_0_xx,  \
                             tk_0_xxx, \
                             tk_0_xxy, \
                             tk_0_xxz, \
                             tk_0_xy,  \
                             tk_0_xyy, \
                             tk_0_xyz, \
                             tk_0_xz,  \
                             tk_0_xzz, \
                             tk_0_yy,  \
                             tk_0_yyy, \
                             tk_0_yyz, \
                             tk_0_yz,  \
                             tk_0_yzz, \
                             tk_0_zz,  \
                             tk_0_zzz, \
                             tk_x_xxx, \
                             tk_x_xxy, \
                             tk_x_xxz, \
                             tk_x_xyy, \
                             tk_x_xyz, \
                             tk_x_xzz, \
                             tk_x_yyy, \
                             tk_x_yyz, \
                             tk_x_yzz, \
                             tk_x_zzz, \
                             ts_x_xxx, \
                             ts_x_xxy, \
                             ts_x_xxz, \
                             ts_x_xyy, \
                             ts_x_xyz, \
                             ts_x_xzz, \
                             ts_x_yyy, \
                             ts_x_yyz, \
                             ts_x_yzz, \
                             ts_x_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_x_xxx[i] = 3.0 * tk_0_xx[i] * fe_0 + tk_0_xxx[i] * pa_x[i] + 2.0 * ts_x_xxx[i] * fz_0;

        tk_x_xxy[i] = 2.0 * tk_0_xy[i] * fe_0 + tk_0_xxy[i] * pa_x[i] + 2.0 * ts_x_xxy[i] * fz_0;

        tk_x_xxz[i] = 2.0 * tk_0_xz[i] * fe_0 + tk_0_xxz[i] * pa_x[i] + 2.0 * ts_x_xxz[i] * fz_0;

        tk_x_xyy[i] = tk_0_yy[i] * fe_0 + tk_0_xyy[i] * pa_x[i] + 2.0 * ts_x_xyy[i] * fz_0;

        tk_x_xyz[i] = tk_0_yz[i] * fe_0 + tk_0_xyz[i] * pa_x[i] + 2.0 * ts_x_xyz[i] * fz_0;

        tk_x_xzz[i] = tk_0_zz[i] * fe_0 + tk_0_xzz[i] * pa_x[i] + 2.0 * ts_x_xzz[i] * fz_0;

        tk_x_yyy[i] = tk_0_yyy[i] * pa_x[i] + 2.0 * ts_x_yyy[i] * fz_0;

        tk_x_yyz[i] = tk_0_yyz[i] * pa_x[i] + 2.0 * ts_x_yyz[i] * fz_0;

        tk_x_yzz[i] = tk_0_yzz[i] * pa_x[i] + 2.0 * ts_x_yzz[i] * fz_0;

        tk_x_zzz[i] = tk_0_zzz[i] * pa_x[i] + 2.0 * ts_x_zzz[i] * fz_0;
    }

    // Set up 10-20 components of targeted buffer : PF

    auto tk_y_xxx = pbuffer.data(idx_kin_pf + 10);

    auto tk_y_xxy = pbuffer.data(idx_kin_pf + 11);

    auto tk_y_xxz = pbuffer.data(idx_kin_pf + 12);

    auto tk_y_xyy = pbuffer.data(idx_kin_pf + 13);

    auto tk_y_xyz = pbuffer.data(idx_kin_pf + 14);

    auto tk_y_xzz = pbuffer.data(idx_kin_pf + 15);

    auto tk_y_yyy = pbuffer.data(idx_kin_pf + 16);

    auto tk_y_yyz = pbuffer.data(idx_kin_pf + 17);

    auto tk_y_yzz = pbuffer.data(idx_kin_pf + 18);

    auto tk_y_zzz = pbuffer.data(idx_kin_pf + 19);

#pragma omp simd aligned(pa_y,         \
                             tk_0_xx,  \
                             tk_0_xxx, \
                             tk_0_xxy, \
                             tk_0_xxz, \
                             tk_0_xy,  \
                             tk_0_xyy, \
                             tk_0_xyz, \
                             tk_0_xz,  \
                             tk_0_xzz, \
                             tk_0_yy,  \
                             tk_0_yyy, \
                             tk_0_yyz, \
                             tk_0_yz,  \
                             tk_0_yzz, \
                             tk_0_zz,  \
                             tk_0_zzz, \
                             tk_y_xxx, \
                             tk_y_xxy, \
                             tk_y_xxz, \
                             tk_y_xyy, \
                             tk_y_xyz, \
                             tk_y_xzz, \
                             tk_y_yyy, \
                             tk_y_yyz, \
                             tk_y_yzz, \
                             tk_y_zzz, \
                             ts_y_xxx, \
                             ts_y_xxy, \
                             ts_y_xxz, \
                             ts_y_xyy, \
                             ts_y_xyz, \
                             ts_y_xzz, \
                             ts_y_yyy, \
                             ts_y_yyz, \
                             ts_y_yzz, \
                             ts_y_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_y_xxx[i] = tk_0_xxx[i] * pa_y[i] + 2.0 * ts_y_xxx[i] * fz_0;

        tk_y_xxy[i] = tk_0_xx[i] * fe_0 + tk_0_xxy[i] * pa_y[i] + 2.0 * ts_y_xxy[i] * fz_0;

        tk_y_xxz[i] = tk_0_xxz[i] * pa_y[i] + 2.0 * ts_y_xxz[i] * fz_0;

        tk_y_xyy[i] = 2.0 * tk_0_xy[i] * fe_0 + tk_0_xyy[i] * pa_y[i] + 2.0 * ts_y_xyy[i] * fz_0;

        tk_y_xyz[i] = tk_0_xz[i] * fe_0 + tk_0_xyz[i] * pa_y[i] + 2.0 * ts_y_xyz[i] * fz_0;

        tk_y_xzz[i] = tk_0_xzz[i] * pa_y[i] + 2.0 * ts_y_xzz[i] * fz_0;

        tk_y_yyy[i] = 3.0 * tk_0_yy[i] * fe_0 + tk_0_yyy[i] * pa_y[i] + 2.0 * ts_y_yyy[i] * fz_0;

        tk_y_yyz[i] = 2.0 * tk_0_yz[i] * fe_0 + tk_0_yyz[i] * pa_y[i] + 2.0 * ts_y_yyz[i] * fz_0;

        tk_y_yzz[i] = tk_0_zz[i] * fe_0 + tk_0_yzz[i] * pa_y[i] + 2.0 * ts_y_yzz[i] * fz_0;

        tk_y_zzz[i] = tk_0_zzz[i] * pa_y[i] + 2.0 * ts_y_zzz[i] * fz_0;
    }

    // Set up 20-30 components of targeted buffer : PF

    auto tk_z_xxx = pbuffer.data(idx_kin_pf + 20);

    auto tk_z_xxy = pbuffer.data(idx_kin_pf + 21);

    auto tk_z_xxz = pbuffer.data(idx_kin_pf + 22);

    auto tk_z_xyy = pbuffer.data(idx_kin_pf + 23);

    auto tk_z_xyz = pbuffer.data(idx_kin_pf + 24);

    auto tk_z_xzz = pbuffer.data(idx_kin_pf + 25);

    auto tk_z_yyy = pbuffer.data(idx_kin_pf + 26);

    auto tk_z_yyz = pbuffer.data(idx_kin_pf + 27);

    auto tk_z_yzz = pbuffer.data(idx_kin_pf + 28);

    auto tk_z_zzz = pbuffer.data(idx_kin_pf + 29);

#pragma omp simd aligned(pa_z,         \
                             tk_0_xx,  \
                             tk_0_xxx, \
                             tk_0_xxy, \
                             tk_0_xxz, \
                             tk_0_xy,  \
                             tk_0_xyy, \
                             tk_0_xyz, \
                             tk_0_xz,  \
                             tk_0_xzz, \
                             tk_0_yy,  \
                             tk_0_yyy, \
                             tk_0_yyz, \
                             tk_0_yz,  \
                             tk_0_yzz, \
                             tk_0_zz,  \
                             tk_0_zzz, \
                             tk_z_xxx, \
                             tk_z_xxy, \
                             tk_z_xxz, \
                             tk_z_xyy, \
                             tk_z_xyz, \
                             tk_z_xzz, \
                             tk_z_yyy, \
                             tk_z_yyz, \
                             tk_z_yzz, \
                             tk_z_zzz, \
                             ts_z_xxx, \
                             ts_z_xxy, \
                             ts_z_xxz, \
                             ts_z_xyy, \
                             ts_z_xyz, \
                             ts_z_xzz, \
                             ts_z_yyy, \
                             ts_z_yyz, \
                             ts_z_yzz, \
                             ts_z_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_z_xxx[i] = tk_0_xxx[i] * pa_z[i] + 2.0 * ts_z_xxx[i] * fz_0;

        tk_z_xxy[i] = tk_0_xxy[i] * pa_z[i] + 2.0 * ts_z_xxy[i] * fz_0;

        tk_z_xxz[i] = tk_0_xx[i] * fe_0 + tk_0_xxz[i] * pa_z[i] + 2.0 * ts_z_xxz[i] * fz_0;

        tk_z_xyy[i] = tk_0_xyy[i] * pa_z[i] + 2.0 * ts_z_xyy[i] * fz_0;

        tk_z_xyz[i] = tk_0_xy[i] * fe_0 + tk_0_xyz[i] * pa_z[i] + 2.0 * ts_z_xyz[i] * fz_0;

        tk_z_xzz[i] = 2.0 * tk_0_xz[i] * fe_0 + tk_0_xzz[i] * pa_z[i] + 2.0 * ts_z_xzz[i] * fz_0;

        tk_z_yyy[i] = tk_0_yyy[i] * pa_z[i] + 2.0 * ts_z_yyy[i] * fz_0;

        tk_z_yyz[i] = tk_0_yy[i] * fe_0 + tk_0_yyz[i] * pa_z[i] + 2.0 * ts_z_yyz[i] * fz_0;

        tk_z_yzz[i] = 2.0 * tk_0_yz[i] * fe_0 + tk_0_yzz[i] * pa_z[i] + 2.0 * ts_z_yzz[i] * fz_0;

        tk_z_zzz[i] = 3.0 * tk_0_zz[i] * fe_0 + tk_0_zzz[i] * pa_z[i] + 2.0 * ts_z_zzz[i] * fz_0;
    }
}

}  // namespace kinrec
