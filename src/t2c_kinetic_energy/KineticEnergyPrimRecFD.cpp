#include "KineticEnergyPrimRecFD.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_fd(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_fd,
                            const size_t              idx_ovl_pd,
                            const size_t              idx_kin_pd,
                            const size_t              idx_kin_dp,
                            const size_t              idx_kin_dd,
                            const size_t              idx_ovl_fd,
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

    // Set up components of auxiliary buffer : DP

    auto tk_xx_x = pbuffer.data(idx_kin_dp);

    auto tk_xx_y = pbuffer.data(idx_kin_dp + 1);

    auto tk_xx_z = pbuffer.data(idx_kin_dp + 2);

    auto tk_yy_x = pbuffer.data(idx_kin_dp + 9);

    auto tk_yy_y = pbuffer.data(idx_kin_dp + 10);

    auto tk_yy_z = pbuffer.data(idx_kin_dp + 11);

    auto tk_zz_x = pbuffer.data(idx_kin_dp + 15);

    auto tk_zz_y = pbuffer.data(idx_kin_dp + 16);

    auto tk_zz_z = pbuffer.data(idx_kin_dp + 17);

    // Set up components of auxiliary buffer : DD

    auto tk_xx_xx = pbuffer.data(idx_kin_dd);

    auto tk_xx_xy = pbuffer.data(idx_kin_dd + 1);

    auto tk_xx_xz = pbuffer.data(idx_kin_dd + 2);

    auto tk_xx_yy = pbuffer.data(idx_kin_dd + 3);

    auto tk_xx_yz = pbuffer.data(idx_kin_dd + 4);

    auto tk_xx_zz = pbuffer.data(idx_kin_dd + 5);

    auto tk_xy_xy = pbuffer.data(idx_kin_dd + 7);

    auto tk_xz_xx = pbuffer.data(idx_kin_dd + 12);

    auto tk_xz_xz = pbuffer.data(idx_kin_dd + 14);

    auto tk_yy_xx = pbuffer.data(idx_kin_dd + 18);

    auto tk_yy_xy = pbuffer.data(idx_kin_dd + 19);

    auto tk_yy_xz = pbuffer.data(idx_kin_dd + 20);

    auto tk_yy_yy = pbuffer.data(idx_kin_dd + 21);

    auto tk_yy_yz = pbuffer.data(idx_kin_dd + 22);

    auto tk_yy_zz = pbuffer.data(idx_kin_dd + 23);

    auto tk_yz_yy = pbuffer.data(idx_kin_dd + 27);

    auto tk_yz_yz = pbuffer.data(idx_kin_dd + 28);

    auto tk_yz_zz = pbuffer.data(idx_kin_dd + 29);

    auto tk_zz_xx = pbuffer.data(idx_kin_dd + 30);

    auto tk_zz_xy = pbuffer.data(idx_kin_dd + 31);

    auto tk_zz_xz = pbuffer.data(idx_kin_dd + 32);

    auto tk_zz_yy = pbuffer.data(idx_kin_dd + 33);

    auto tk_zz_yz = pbuffer.data(idx_kin_dd + 34);

    auto tk_zz_zz = pbuffer.data(idx_kin_dd + 35);

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_ovl_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_ovl_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_ovl_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_ovl_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_ovl_fd + 10);

    auto ts_xxy_zz = pbuffer.data(idx_ovl_fd + 11);

    auto ts_xxz_xx = pbuffer.data(idx_ovl_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_ovl_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_ovl_fd + 14);

    auto ts_xxz_yy = pbuffer.data(idx_ovl_fd + 15);

    auto ts_xxz_yz = pbuffer.data(idx_ovl_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_ovl_fd + 17);

    auto ts_xyy_xx = pbuffer.data(idx_ovl_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_ovl_fd + 19);

    auto ts_xyy_xz = pbuffer.data(idx_ovl_fd + 20);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_ovl_fd + 23);

    auto ts_xyz_xx = pbuffer.data(idx_ovl_fd + 24);

    auto ts_xyz_xy = pbuffer.data(idx_ovl_fd + 25);

    auto ts_xyz_xz = pbuffer.data(idx_ovl_fd + 26);

    auto ts_xyz_yy = pbuffer.data(idx_ovl_fd + 27);

    auto ts_xyz_yz = pbuffer.data(idx_ovl_fd + 28);

    auto ts_xyz_zz = pbuffer.data(idx_ovl_fd + 29);

    auto ts_xzz_xx = pbuffer.data(idx_ovl_fd + 30);

    auto ts_xzz_xy = pbuffer.data(idx_ovl_fd + 31);

    auto ts_xzz_xz = pbuffer.data(idx_ovl_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_ovl_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_ovl_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_ovl_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_ovl_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_ovl_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_ovl_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_ovl_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_ovl_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_ovl_fd + 41);

    auto ts_yyz_xx = pbuffer.data(idx_ovl_fd + 42);

    auto ts_yyz_xy = pbuffer.data(idx_ovl_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_ovl_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_ovl_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_ovl_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_ovl_fd + 47);

    auto ts_yzz_xx = pbuffer.data(idx_ovl_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_ovl_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_ovl_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_ovl_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_ovl_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_ovl_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_ovl_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_ovl_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_ovl_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_ovl_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_ovl_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_ovl_fd + 59);

    // Set up 0-6 components of targeted buffer : FD

    auto tk_xxx_xx = pbuffer.data(idx_kin_fd);

    auto tk_xxx_xy = pbuffer.data(idx_kin_fd + 1);

    auto tk_xxx_xz = pbuffer.data(idx_kin_fd + 2);

    auto tk_xxx_yy = pbuffer.data(idx_kin_fd + 3);

    auto tk_xxx_yz = pbuffer.data(idx_kin_fd + 4);

    auto tk_xxx_zz = pbuffer.data(idx_kin_fd + 5);

#pragma omp simd aligned(pa_x,          \
                             tk_x_xx,   \
                             tk_x_xy,   \
                             tk_x_xz,   \
                             tk_x_yy,   \
                             tk_x_yz,   \
                             tk_x_zz,   \
                             tk_xx_x,   \
                             tk_xx_xx,  \
                             tk_xx_xy,  \
                             tk_xx_xz,  \
                             tk_xx_y,   \
                             tk_xx_yy,  \
                             tk_xx_yz,  \
                             tk_xx_z,   \
                             tk_xx_zz,  \
                             tk_xxx_xx, \
                             tk_xxx_xy, \
                             tk_xxx_xz, \
                             tk_xxx_yy, \
                             tk_xxx_yz, \
                             tk_xxx_zz, \
                             ts_x_xx,   \
                             ts_x_xy,   \
                             ts_x_xz,   \
                             ts_x_yy,   \
                             ts_x_yz,   \
                             ts_x_zz,   \
                             ts_xxx_xx, \
                             ts_xxx_xy, \
                             ts_xxx_xz, \
                             ts_xxx_yy, \
                             ts_xxx_yz, \
                             ts_xxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_xx[i] = -4.0 * ts_x_xx[i] * fbe_0 * fz_0 + 2.0 * tk_x_xx[i] * fe_0 + 2.0 * tk_xx_x[i] * fe_0 +
                       tk_xx_xx[i] * pa_x[i] + 2.0 * ts_xxx_xx[i] * fz_0;

        tk_xxx_xy[i] = -4.0 * ts_x_xy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xy[i] * fe_0 + tk_xx_y[i] * fe_0 +
                       tk_xx_xy[i] * pa_x[i] + 2.0 * ts_xxx_xy[i] * fz_0;

        tk_xxx_xz[i] = -4.0 * ts_x_xz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xz[i] * fe_0 + tk_xx_z[i] * fe_0 +
                       tk_xx_xz[i] * pa_x[i] + 2.0 * ts_xxx_xz[i] * fz_0;

        tk_xxx_yy[i] = -4.0 * ts_x_yy[i] * fbe_0 * fz_0 + 2.0 * tk_x_yy[i] * fe_0 + tk_xx_yy[i] * pa_x[i] +
                       2.0 * ts_xxx_yy[i] * fz_0;

        tk_xxx_yz[i] = -4.0 * ts_x_yz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yz[i] * fe_0 + tk_xx_yz[i] * pa_x[i] +
                       2.0 * ts_xxx_yz[i] * fz_0;

        tk_xxx_zz[i] = -4.0 * ts_x_zz[i] * fbe_0 * fz_0 + 2.0 * tk_x_zz[i] * fe_0 + tk_xx_zz[i] * pa_x[i] +
                       2.0 * ts_xxx_zz[i] * fz_0;
    }

    // Set up 6-12 components of targeted buffer : FD

    auto tk_xxy_xx = pbuffer.data(idx_kin_fd + 6);

    auto tk_xxy_xy = pbuffer.data(idx_kin_fd + 7);

    auto tk_xxy_xz = pbuffer.data(idx_kin_fd + 8);

    auto tk_xxy_yy = pbuffer.data(idx_kin_fd + 9);

    auto tk_xxy_yz = pbuffer.data(idx_kin_fd + 10);

    auto tk_xxy_zz = pbuffer.data(idx_kin_fd + 11);

#pragma omp simd aligned(pa_y,          \
                             tk_xx_x,   \
                             tk_xx_xx,  \
                             tk_xx_xy,  \
                             tk_xx_xz,  \
                             tk_xx_y,   \
                             tk_xx_yy,  \
                             tk_xx_yz,  \
                             tk_xx_z,   \
                             tk_xx_zz,  \
                             tk_xxy_xx, \
                             tk_xxy_xy, \
                             tk_xxy_xz, \
                             tk_xxy_yy, \
                             tk_xxy_yz, \
                             tk_xxy_zz, \
                             ts_xxy_xx, \
                             ts_xxy_xy, \
                             ts_xxy_xz, \
                             ts_xxy_yy, \
                             ts_xxy_yz, \
                             ts_xxy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxy_xx[i] = tk_xx_xx[i] * pa_y[i] + 2.0 * ts_xxy_xx[i] * fz_0;

        tk_xxy_xy[i] = tk_xx_x[i] * fe_0 + tk_xx_xy[i] * pa_y[i] + 2.0 * ts_xxy_xy[i] * fz_0;

        tk_xxy_xz[i] = tk_xx_xz[i] * pa_y[i] + 2.0 * ts_xxy_xz[i] * fz_0;

        tk_xxy_yy[i] = 2.0 * tk_xx_y[i] * fe_0 + tk_xx_yy[i] * pa_y[i] + 2.0 * ts_xxy_yy[i] * fz_0;

        tk_xxy_yz[i] = tk_xx_z[i] * fe_0 + tk_xx_yz[i] * pa_y[i] + 2.0 * ts_xxy_yz[i] * fz_0;

        tk_xxy_zz[i] = tk_xx_zz[i] * pa_y[i] + 2.0 * ts_xxy_zz[i] * fz_0;
    }

    // Set up 12-18 components of targeted buffer : FD

    auto tk_xxz_xx = pbuffer.data(idx_kin_fd + 12);

    auto tk_xxz_xy = pbuffer.data(idx_kin_fd + 13);

    auto tk_xxz_xz = pbuffer.data(idx_kin_fd + 14);

    auto tk_xxz_yy = pbuffer.data(idx_kin_fd + 15);

    auto tk_xxz_yz = pbuffer.data(idx_kin_fd + 16);

    auto tk_xxz_zz = pbuffer.data(idx_kin_fd + 17);

#pragma omp simd aligned(pa_z,          \
                             tk_xx_x,   \
                             tk_xx_xx,  \
                             tk_xx_xy,  \
                             tk_xx_xz,  \
                             tk_xx_y,   \
                             tk_xx_yy,  \
                             tk_xx_yz,  \
                             tk_xx_z,   \
                             tk_xx_zz,  \
                             tk_xxz_xx, \
                             tk_xxz_xy, \
                             tk_xxz_xz, \
                             tk_xxz_yy, \
                             tk_xxz_yz, \
                             tk_xxz_zz, \
                             ts_xxz_xx, \
                             ts_xxz_xy, \
                             ts_xxz_xz, \
                             ts_xxz_yy, \
                             ts_xxz_yz, \
                             ts_xxz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxz_xx[i] = tk_xx_xx[i] * pa_z[i] + 2.0 * ts_xxz_xx[i] * fz_0;

        tk_xxz_xy[i] = tk_xx_xy[i] * pa_z[i] + 2.0 * ts_xxz_xy[i] * fz_0;

        tk_xxz_xz[i] = tk_xx_x[i] * fe_0 + tk_xx_xz[i] * pa_z[i] + 2.0 * ts_xxz_xz[i] * fz_0;

        tk_xxz_yy[i] = tk_xx_yy[i] * pa_z[i] + 2.0 * ts_xxz_yy[i] * fz_0;

        tk_xxz_yz[i] = tk_xx_y[i] * fe_0 + tk_xx_yz[i] * pa_z[i] + 2.0 * ts_xxz_yz[i] * fz_0;

        tk_xxz_zz[i] = 2.0 * tk_xx_z[i] * fe_0 + tk_xx_zz[i] * pa_z[i] + 2.0 * ts_xxz_zz[i] * fz_0;
    }

    // Set up 18-24 components of targeted buffer : FD

    auto tk_xyy_xx = pbuffer.data(idx_kin_fd + 18);

    auto tk_xyy_xy = pbuffer.data(idx_kin_fd + 19);

    auto tk_xyy_xz = pbuffer.data(idx_kin_fd + 20);

    auto tk_xyy_yy = pbuffer.data(idx_kin_fd + 21);

    auto tk_xyy_yz = pbuffer.data(idx_kin_fd + 22);

    auto tk_xyy_zz = pbuffer.data(idx_kin_fd + 23);

#pragma omp simd aligned(pa_x,          \
                             tk_xyy_xx, \
                             tk_xyy_xy, \
                             tk_xyy_xz, \
                             tk_xyy_yy, \
                             tk_xyy_yz, \
                             tk_xyy_zz, \
                             tk_yy_x,   \
                             tk_yy_xx,  \
                             tk_yy_xy,  \
                             tk_yy_xz,  \
                             tk_yy_y,   \
                             tk_yy_yy,  \
                             tk_yy_yz,  \
                             tk_yy_z,   \
                             tk_yy_zz,  \
                             ts_xyy_xx, \
                             ts_xyy_xy, \
                             ts_xyy_xz, \
                             ts_xyy_yy, \
                             ts_xyy_yz, \
                             ts_xyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyy_xx[i] = 2.0 * tk_yy_x[i] * fe_0 + tk_yy_xx[i] * pa_x[i] + 2.0 * ts_xyy_xx[i] * fz_0;

        tk_xyy_xy[i] = tk_yy_y[i] * fe_0 + tk_yy_xy[i] * pa_x[i] + 2.0 * ts_xyy_xy[i] * fz_0;

        tk_xyy_xz[i] = tk_yy_z[i] * fe_0 + tk_yy_xz[i] * pa_x[i] + 2.0 * ts_xyy_xz[i] * fz_0;

        tk_xyy_yy[i] = tk_yy_yy[i] * pa_x[i] + 2.0 * ts_xyy_yy[i] * fz_0;

        tk_xyy_yz[i] = tk_yy_yz[i] * pa_x[i] + 2.0 * ts_xyy_yz[i] * fz_0;

        tk_xyy_zz[i] = tk_yy_zz[i] * pa_x[i] + 2.0 * ts_xyy_zz[i] * fz_0;
    }

    // Set up 24-30 components of targeted buffer : FD

    auto tk_xyz_xx = pbuffer.data(idx_kin_fd + 24);

    auto tk_xyz_xy = pbuffer.data(idx_kin_fd + 25);

    auto tk_xyz_xz = pbuffer.data(idx_kin_fd + 26);

    auto tk_xyz_yy = pbuffer.data(idx_kin_fd + 27);

    auto tk_xyz_yz = pbuffer.data(idx_kin_fd + 28);

    auto tk_xyz_zz = pbuffer.data(idx_kin_fd + 29);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             pa_z,      \
                             tk_xy_xy,  \
                             tk_xyz_xx, \
                             tk_xyz_xy, \
                             tk_xyz_xz, \
                             tk_xyz_yy, \
                             tk_xyz_yz, \
                             tk_xyz_zz, \
                             tk_xz_xx,  \
                             tk_xz_xz,  \
                             tk_yz_yy,  \
                             tk_yz_yz,  \
                             tk_yz_zz,  \
                             ts_xyz_xx, \
                             ts_xyz_xy, \
                             ts_xyz_xz, \
                             ts_xyz_yy, \
                             ts_xyz_yz, \
                             ts_xyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xyz_xx[i] = tk_xz_xx[i] * pa_y[i] + 2.0 * ts_xyz_xx[i] * fz_0;

        tk_xyz_xy[i] = tk_xy_xy[i] * pa_z[i] + 2.0 * ts_xyz_xy[i] * fz_0;

        tk_xyz_xz[i] = tk_xz_xz[i] * pa_y[i] + 2.0 * ts_xyz_xz[i] * fz_0;

        tk_xyz_yy[i] = tk_yz_yy[i] * pa_x[i] + 2.0 * ts_xyz_yy[i] * fz_0;

        tk_xyz_yz[i] = tk_yz_yz[i] * pa_x[i] + 2.0 * ts_xyz_yz[i] * fz_0;

        tk_xyz_zz[i] = tk_yz_zz[i] * pa_x[i] + 2.0 * ts_xyz_zz[i] * fz_0;
    }

    // Set up 30-36 components of targeted buffer : FD

    auto tk_xzz_xx = pbuffer.data(idx_kin_fd + 30);

    auto tk_xzz_xy = pbuffer.data(idx_kin_fd + 31);

    auto tk_xzz_xz = pbuffer.data(idx_kin_fd + 32);

    auto tk_xzz_yy = pbuffer.data(idx_kin_fd + 33);

    auto tk_xzz_yz = pbuffer.data(idx_kin_fd + 34);

    auto tk_xzz_zz = pbuffer.data(idx_kin_fd + 35);

#pragma omp simd aligned(pa_x,          \
                             tk_xzz_xx, \
                             tk_xzz_xy, \
                             tk_xzz_xz, \
                             tk_xzz_yy, \
                             tk_xzz_yz, \
                             tk_xzz_zz, \
                             tk_zz_x,   \
                             tk_zz_xx,  \
                             tk_zz_xy,  \
                             tk_zz_xz,  \
                             tk_zz_y,   \
                             tk_zz_yy,  \
                             tk_zz_yz,  \
                             tk_zz_z,   \
                             tk_zz_zz,  \
                             ts_xzz_xx, \
                             ts_xzz_xy, \
                             ts_xzz_xz, \
                             ts_xzz_yy, \
                             ts_xzz_yz, \
                             ts_xzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzz_xx[i] = 2.0 * tk_zz_x[i] * fe_0 + tk_zz_xx[i] * pa_x[i] + 2.0 * ts_xzz_xx[i] * fz_0;

        tk_xzz_xy[i] = tk_zz_y[i] * fe_0 + tk_zz_xy[i] * pa_x[i] + 2.0 * ts_xzz_xy[i] * fz_0;

        tk_xzz_xz[i] = tk_zz_z[i] * fe_0 + tk_zz_xz[i] * pa_x[i] + 2.0 * ts_xzz_xz[i] * fz_0;

        tk_xzz_yy[i] = tk_zz_yy[i] * pa_x[i] + 2.0 * ts_xzz_yy[i] * fz_0;

        tk_xzz_yz[i] = tk_zz_yz[i] * pa_x[i] + 2.0 * ts_xzz_yz[i] * fz_0;

        tk_xzz_zz[i] = tk_zz_zz[i] * pa_x[i] + 2.0 * ts_xzz_zz[i] * fz_0;
    }

    // Set up 36-42 components of targeted buffer : FD

    auto tk_yyy_xx = pbuffer.data(idx_kin_fd + 36);

    auto tk_yyy_xy = pbuffer.data(idx_kin_fd + 37);

    auto tk_yyy_xz = pbuffer.data(idx_kin_fd + 38);

    auto tk_yyy_yy = pbuffer.data(idx_kin_fd + 39);

    auto tk_yyy_yz = pbuffer.data(idx_kin_fd + 40);

    auto tk_yyy_zz = pbuffer.data(idx_kin_fd + 41);

#pragma omp simd aligned(pa_y,          \
                             tk_y_xx,   \
                             tk_y_xy,   \
                             tk_y_xz,   \
                             tk_y_yy,   \
                             tk_y_yz,   \
                             tk_y_zz,   \
                             tk_yy_x,   \
                             tk_yy_xx,  \
                             tk_yy_xy,  \
                             tk_yy_xz,  \
                             tk_yy_y,   \
                             tk_yy_yy,  \
                             tk_yy_yz,  \
                             tk_yy_z,   \
                             tk_yy_zz,  \
                             tk_yyy_xx, \
                             tk_yyy_xy, \
                             tk_yyy_xz, \
                             tk_yyy_yy, \
                             tk_yyy_yz, \
                             tk_yyy_zz, \
                             ts_y_xx,   \
                             ts_y_xy,   \
                             ts_y_xz,   \
                             ts_y_yy,   \
                             ts_y_yz,   \
                             ts_y_zz,   \
                             ts_yyy_xx, \
                             ts_yyy_xy, \
                             ts_yyy_xz, \
                             ts_yyy_yy, \
                             ts_yyy_yz, \
                             ts_yyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyy_xx[i] = -4.0 * ts_y_xx[i] * fbe_0 * fz_0 + 2.0 * tk_y_xx[i] * fe_0 + tk_yy_xx[i] * pa_y[i] +
                       2.0 * ts_yyy_xx[i] * fz_0;

        tk_yyy_xy[i] = -4.0 * ts_y_xy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xy[i] * fe_0 + tk_yy_x[i] * fe_0 +
                       tk_yy_xy[i] * pa_y[i] + 2.0 * ts_yyy_xy[i] * fz_0;

        tk_yyy_xz[i] = -4.0 * ts_y_xz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xz[i] * fe_0 + tk_yy_xz[i] * pa_y[i] +
                       2.0 * ts_yyy_xz[i] * fz_0;

        tk_yyy_yy[i] = -4.0 * ts_y_yy[i] * fbe_0 * fz_0 + 2.0 * tk_y_yy[i] * fe_0 + 2.0 * tk_yy_y[i] * fe_0 +
                       tk_yy_yy[i] * pa_y[i] + 2.0 * ts_yyy_yy[i] * fz_0;

        tk_yyy_yz[i] = -4.0 * ts_y_yz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yz[i] * fe_0 + tk_yy_z[i] * fe_0 +
                       tk_yy_yz[i] * pa_y[i] + 2.0 * ts_yyy_yz[i] * fz_0;

        tk_yyy_zz[i] = -4.0 * ts_y_zz[i] * fbe_0 * fz_0 + 2.0 * tk_y_zz[i] * fe_0 + tk_yy_zz[i] * pa_y[i] +
                       2.0 * ts_yyy_zz[i] * fz_0;
    }

    // Set up 42-48 components of targeted buffer : FD

    auto tk_yyz_xx = pbuffer.data(idx_kin_fd + 42);

    auto tk_yyz_xy = pbuffer.data(idx_kin_fd + 43);

    auto tk_yyz_xz = pbuffer.data(idx_kin_fd + 44);

    auto tk_yyz_yy = pbuffer.data(idx_kin_fd + 45);

    auto tk_yyz_yz = pbuffer.data(idx_kin_fd + 46);

    auto tk_yyz_zz = pbuffer.data(idx_kin_fd + 47);

#pragma omp simd aligned(pa_z,          \
                             tk_yy_x,   \
                             tk_yy_xx,  \
                             tk_yy_xy,  \
                             tk_yy_xz,  \
                             tk_yy_y,   \
                             tk_yy_yy,  \
                             tk_yy_yz,  \
                             tk_yy_z,   \
                             tk_yy_zz,  \
                             tk_yyz_xx, \
                             tk_yyz_xy, \
                             tk_yyz_xz, \
                             tk_yyz_yy, \
                             tk_yyz_yz, \
                             tk_yyz_zz, \
                             ts_yyz_xx, \
                             ts_yyz_xy, \
                             ts_yyz_xz, \
                             ts_yyz_yy, \
                             ts_yyz_yz, \
                             ts_yyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyz_xx[i] = tk_yy_xx[i] * pa_z[i] + 2.0 * ts_yyz_xx[i] * fz_0;

        tk_yyz_xy[i] = tk_yy_xy[i] * pa_z[i] + 2.0 * ts_yyz_xy[i] * fz_0;

        tk_yyz_xz[i] = tk_yy_x[i] * fe_0 + tk_yy_xz[i] * pa_z[i] + 2.0 * ts_yyz_xz[i] * fz_0;

        tk_yyz_yy[i] = tk_yy_yy[i] * pa_z[i] + 2.0 * ts_yyz_yy[i] * fz_0;

        tk_yyz_yz[i] = tk_yy_y[i] * fe_0 + tk_yy_yz[i] * pa_z[i] + 2.0 * ts_yyz_yz[i] * fz_0;

        tk_yyz_zz[i] = 2.0 * tk_yy_z[i] * fe_0 + tk_yy_zz[i] * pa_z[i] + 2.0 * ts_yyz_zz[i] * fz_0;
    }

    // Set up 48-54 components of targeted buffer : FD

    auto tk_yzz_xx = pbuffer.data(idx_kin_fd + 48);

    auto tk_yzz_xy = pbuffer.data(idx_kin_fd + 49);

    auto tk_yzz_xz = pbuffer.data(idx_kin_fd + 50);

    auto tk_yzz_yy = pbuffer.data(idx_kin_fd + 51);

    auto tk_yzz_yz = pbuffer.data(idx_kin_fd + 52);

    auto tk_yzz_zz = pbuffer.data(idx_kin_fd + 53);

#pragma omp simd aligned(pa_y,          \
                             tk_yzz_xx, \
                             tk_yzz_xy, \
                             tk_yzz_xz, \
                             tk_yzz_yy, \
                             tk_yzz_yz, \
                             tk_yzz_zz, \
                             tk_zz_x,   \
                             tk_zz_xx,  \
                             tk_zz_xy,  \
                             tk_zz_xz,  \
                             tk_zz_y,   \
                             tk_zz_yy,  \
                             tk_zz_yz,  \
                             tk_zz_z,   \
                             tk_zz_zz,  \
                             ts_yzz_xx, \
                             ts_yzz_xy, \
                             ts_yzz_xz, \
                             ts_yzz_yy, \
                             ts_yzz_yz, \
                             ts_yzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzz_xx[i] = tk_zz_xx[i] * pa_y[i] + 2.0 * ts_yzz_xx[i] * fz_0;

        tk_yzz_xy[i] = tk_zz_x[i] * fe_0 + tk_zz_xy[i] * pa_y[i] + 2.0 * ts_yzz_xy[i] * fz_0;

        tk_yzz_xz[i] = tk_zz_xz[i] * pa_y[i] + 2.0 * ts_yzz_xz[i] * fz_0;

        tk_yzz_yy[i] = 2.0 * tk_zz_y[i] * fe_0 + tk_zz_yy[i] * pa_y[i] + 2.0 * ts_yzz_yy[i] * fz_0;

        tk_yzz_yz[i] = tk_zz_z[i] * fe_0 + tk_zz_yz[i] * pa_y[i] + 2.0 * ts_yzz_yz[i] * fz_0;

        tk_yzz_zz[i] = tk_zz_zz[i] * pa_y[i] + 2.0 * ts_yzz_zz[i] * fz_0;
    }

    // Set up 54-60 components of targeted buffer : FD

    auto tk_zzz_xx = pbuffer.data(idx_kin_fd + 54);

    auto tk_zzz_xy = pbuffer.data(idx_kin_fd + 55);

    auto tk_zzz_xz = pbuffer.data(idx_kin_fd + 56);

    auto tk_zzz_yy = pbuffer.data(idx_kin_fd + 57);

    auto tk_zzz_yz = pbuffer.data(idx_kin_fd + 58);

    auto tk_zzz_zz = pbuffer.data(idx_kin_fd + 59);

#pragma omp simd aligned(pa_z,          \
                             tk_z_xx,   \
                             tk_z_xy,   \
                             tk_z_xz,   \
                             tk_z_yy,   \
                             tk_z_yz,   \
                             tk_z_zz,   \
                             tk_zz_x,   \
                             tk_zz_xx,  \
                             tk_zz_xy,  \
                             tk_zz_xz,  \
                             tk_zz_y,   \
                             tk_zz_yy,  \
                             tk_zz_yz,  \
                             tk_zz_z,   \
                             tk_zz_zz,  \
                             tk_zzz_xx, \
                             tk_zzz_xy, \
                             tk_zzz_xz, \
                             tk_zzz_yy, \
                             tk_zzz_yz, \
                             tk_zzz_zz, \
                             ts_z_xx,   \
                             ts_z_xy,   \
                             ts_z_xz,   \
                             ts_z_yy,   \
                             ts_z_yz,   \
                             ts_z_zz,   \
                             ts_zzz_xx, \
                             ts_zzz_xy, \
                             ts_zzz_xz, \
                             ts_zzz_yy, \
                             ts_zzz_yz, \
                             ts_zzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzz_xx[i] = -4.0 * ts_z_xx[i] * fbe_0 * fz_0 + 2.0 * tk_z_xx[i] * fe_0 + tk_zz_xx[i] * pa_z[i] +
                       2.0 * ts_zzz_xx[i] * fz_0;

        tk_zzz_xy[i] = -4.0 * ts_z_xy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xy[i] * fe_0 + tk_zz_xy[i] * pa_z[i] +
                       2.0 * ts_zzz_xy[i] * fz_0;

        tk_zzz_xz[i] = -4.0 * ts_z_xz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xz[i] * fe_0 + tk_zz_x[i] * fe_0 +
                       tk_zz_xz[i] * pa_z[i] + 2.0 * ts_zzz_xz[i] * fz_0;

        tk_zzz_yy[i] = -4.0 * ts_z_yy[i] * fbe_0 * fz_0 + 2.0 * tk_z_yy[i] * fe_0 + tk_zz_yy[i] * pa_z[i] +
                       2.0 * ts_zzz_yy[i] * fz_0;

        tk_zzz_yz[i] = -4.0 * ts_z_yz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yz[i] * fe_0 + tk_zz_y[i] * fe_0 +
                       tk_zz_yz[i] * pa_z[i] + 2.0 * ts_zzz_yz[i] * fz_0;

        tk_zzz_zz[i] = -4.0 * ts_z_zz[i] * fbe_0 * fz_0 + 2.0 * tk_z_zz[i] * fe_0 + 2.0 * tk_zz_z[i] * fe_0 +
                       tk_zz_zz[i] * pa_z[i] + 2.0 * ts_zzz_zz[i] * fz_0;
    }
}

}  // namespace kinrec
