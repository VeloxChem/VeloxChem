#include "KineticEnergyPrimRecGD.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_gd(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_gd,
                            const size_t              idx_ovl_dd,
                            const size_t              idx_kin_dd,
                            const size_t              idx_kin_fp,
                            const size_t              idx_kin_fd,
                            const size_t              idx_ovl_gd,
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

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_ovl_dd);

    auto ts_xx_xy = pbuffer.data(idx_ovl_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_ovl_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_ovl_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_ovl_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_ovl_dd + 5);

    auto ts_yy_xx = pbuffer.data(idx_ovl_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_ovl_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_ovl_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_ovl_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_ovl_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_ovl_dd + 23);

    auto ts_zz_xx = pbuffer.data(idx_ovl_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_ovl_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_ovl_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_ovl_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_ovl_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_ovl_dd + 35);

    // Set up components of auxiliary buffer : DD

    auto tk_xx_xx = pbuffer.data(idx_kin_dd);

    auto tk_xx_xy = pbuffer.data(idx_kin_dd + 1);

    auto tk_xx_xz = pbuffer.data(idx_kin_dd + 2);

    auto tk_xx_yy = pbuffer.data(idx_kin_dd + 3);

    auto tk_xx_yz = pbuffer.data(idx_kin_dd + 4);

    auto tk_xx_zz = pbuffer.data(idx_kin_dd + 5);

    auto tk_yy_xx = pbuffer.data(idx_kin_dd + 18);

    auto tk_yy_xy = pbuffer.data(idx_kin_dd + 19);

    auto tk_yy_xz = pbuffer.data(idx_kin_dd + 20);

    auto tk_yy_yy = pbuffer.data(idx_kin_dd + 21);

    auto tk_yy_yz = pbuffer.data(idx_kin_dd + 22);

    auto tk_yy_zz = pbuffer.data(idx_kin_dd + 23);

    auto tk_zz_xx = pbuffer.data(idx_kin_dd + 30);

    auto tk_zz_xy = pbuffer.data(idx_kin_dd + 31);

    auto tk_zz_xz = pbuffer.data(idx_kin_dd + 32);

    auto tk_zz_yy = pbuffer.data(idx_kin_dd + 33);

    auto tk_zz_yz = pbuffer.data(idx_kin_dd + 34);

    auto tk_zz_zz = pbuffer.data(idx_kin_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto tk_xxx_x = pbuffer.data(idx_kin_fp);

    auto tk_xxx_y = pbuffer.data(idx_kin_fp + 1);

    auto tk_xxx_z = pbuffer.data(idx_kin_fp + 2);

    auto tk_xxz_z = pbuffer.data(idx_kin_fp + 8);

    auto tk_xyy_y = pbuffer.data(idx_kin_fp + 10);

    auto tk_xzz_z = pbuffer.data(idx_kin_fp + 17);

    auto tk_yyy_x = pbuffer.data(idx_kin_fp + 18);

    auto tk_yyy_y = pbuffer.data(idx_kin_fp + 19);

    auto tk_yyy_z = pbuffer.data(idx_kin_fp + 20);

    auto tk_yyz_z = pbuffer.data(idx_kin_fp + 23);

    auto tk_yzz_y = pbuffer.data(idx_kin_fp + 25);

    auto tk_yzz_z = pbuffer.data(idx_kin_fp + 26);

    auto tk_zzz_x = pbuffer.data(idx_kin_fp + 27);

    auto tk_zzz_y = pbuffer.data(idx_kin_fp + 28);

    auto tk_zzz_z = pbuffer.data(idx_kin_fp + 29);

    // Set up components of auxiliary buffer : FD

    auto tk_xxx_xx = pbuffer.data(idx_kin_fd);

    auto tk_xxx_xy = pbuffer.data(idx_kin_fd + 1);

    auto tk_xxx_xz = pbuffer.data(idx_kin_fd + 2);

    auto tk_xxx_yy = pbuffer.data(idx_kin_fd + 3);

    auto tk_xxx_yz = pbuffer.data(idx_kin_fd + 4);

    auto tk_xxx_zz = pbuffer.data(idx_kin_fd + 5);

    auto tk_xxy_xx = pbuffer.data(idx_kin_fd + 6);

    auto tk_xxy_xy = pbuffer.data(idx_kin_fd + 7);

    auto tk_xxy_xz = pbuffer.data(idx_kin_fd + 8);

    auto tk_xxy_yy = pbuffer.data(idx_kin_fd + 9);

    auto tk_xxz_xx = pbuffer.data(idx_kin_fd + 12);

    auto tk_xxz_xy = pbuffer.data(idx_kin_fd + 13);

    auto tk_xxz_xz = pbuffer.data(idx_kin_fd + 14);

    auto tk_xxz_yz = pbuffer.data(idx_kin_fd + 16);

    auto tk_xxz_zz = pbuffer.data(idx_kin_fd + 17);

    auto tk_xyy_xx = pbuffer.data(idx_kin_fd + 18);

    auto tk_xyy_xy = pbuffer.data(idx_kin_fd + 19);

    auto tk_xyy_yy = pbuffer.data(idx_kin_fd + 21);

    auto tk_xyy_yz = pbuffer.data(idx_kin_fd + 22);

    auto tk_xyy_zz = pbuffer.data(idx_kin_fd + 23);

    auto tk_xzz_xx = pbuffer.data(idx_kin_fd + 30);

    auto tk_xzz_xz = pbuffer.data(idx_kin_fd + 32);

    auto tk_xzz_yy = pbuffer.data(idx_kin_fd + 33);

    auto tk_xzz_yz = pbuffer.data(idx_kin_fd + 34);

    auto tk_xzz_zz = pbuffer.data(idx_kin_fd + 35);

    auto tk_yyy_xx = pbuffer.data(idx_kin_fd + 36);

    auto tk_yyy_xy = pbuffer.data(idx_kin_fd + 37);

    auto tk_yyy_xz = pbuffer.data(idx_kin_fd + 38);

    auto tk_yyy_yy = pbuffer.data(idx_kin_fd + 39);

    auto tk_yyy_yz = pbuffer.data(idx_kin_fd + 40);

    auto tk_yyy_zz = pbuffer.data(idx_kin_fd + 41);

    auto tk_yyz_xy = pbuffer.data(idx_kin_fd + 43);

    auto tk_yyz_xz = pbuffer.data(idx_kin_fd + 44);

    auto tk_yyz_yy = pbuffer.data(idx_kin_fd + 45);

    auto tk_yyz_yz = pbuffer.data(idx_kin_fd + 46);

    auto tk_yyz_zz = pbuffer.data(idx_kin_fd + 47);

    auto tk_yzz_xx = pbuffer.data(idx_kin_fd + 48);

    auto tk_yzz_xy = pbuffer.data(idx_kin_fd + 49);

    auto tk_yzz_xz = pbuffer.data(idx_kin_fd + 50);

    auto tk_yzz_yy = pbuffer.data(idx_kin_fd + 51);

    auto tk_yzz_yz = pbuffer.data(idx_kin_fd + 52);

    auto tk_yzz_zz = pbuffer.data(idx_kin_fd + 53);

    auto tk_zzz_xx = pbuffer.data(idx_kin_fd + 54);

    auto tk_zzz_xy = pbuffer.data(idx_kin_fd + 55);

    auto tk_zzz_xz = pbuffer.data(idx_kin_fd + 56);

    auto tk_zzz_yy = pbuffer.data(idx_kin_fd + 57);

    auto tk_zzz_yz = pbuffer.data(idx_kin_fd + 58);

    auto tk_zzz_zz = pbuffer.data(idx_kin_fd + 59);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_ovl_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_ovl_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_ovl_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_ovl_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_ovl_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_ovl_gd + 11);

    auto ts_xxxz_xx = pbuffer.data(idx_ovl_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_ovl_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_ovl_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_ovl_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_ovl_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_ovl_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_ovl_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_ovl_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_ovl_gd + 23);

    auto ts_xxyz_xx = pbuffer.data(idx_ovl_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_ovl_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_ovl_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_ovl_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_ovl_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_ovl_gd + 29);

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_ovl_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_ovl_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_ovl_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_ovl_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_ovl_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_ovl_gd + 41);

    auto ts_xyyz_xx = pbuffer.data(idx_ovl_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_ovl_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_ovl_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_ovl_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_ovl_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_ovl_gd + 47);

    auto ts_xyzz_xx = pbuffer.data(idx_ovl_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_ovl_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_ovl_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_ovl_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_ovl_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_ovl_gd + 53);

    auto ts_xzzz_xx = pbuffer.data(idx_ovl_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_ovl_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_ovl_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_ovl_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

    auto ts_yyyz_xx = pbuffer.data(idx_ovl_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_ovl_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_ovl_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_ovl_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_ovl_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_ovl_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_ovl_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_ovl_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_ovl_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_ovl_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_ovl_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_ovl_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_ovl_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_ovl_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_ovl_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_ovl_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_ovl_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_ovl_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_ovl_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_ovl_gd + 89);

    // Set up 0-6 components of targeted buffer : GD

    auto tk_xxxx_xx = pbuffer.data(idx_kin_gd);

    auto tk_xxxx_xy = pbuffer.data(idx_kin_gd + 1);

    auto tk_xxxx_xz = pbuffer.data(idx_kin_gd + 2);

    auto tk_xxxx_yy = pbuffer.data(idx_kin_gd + 3);

    auto tk_xxxx_yz = pbuffer.data(idx_kin_gd + 4);

    auto tk_xxxx_zz = pbuffer.data(idx_kin_gd + 5);

#pragma omp simd aligned(pa_x,           \
                             tk_xx_xx,   \
                             tk_xx_xy,   \
                             tk_xx_xz,   \
                             tk_xx_yy,   \
                             tk_xx_yz,   \
                             tk_xx_zz,   \
                             tk_xxx_x,   \
                             tk_xxx_xx,  \
                             tk_xxx_xy,  \
                             tk_xxx_xz,  \
                             tk_xxx_y,   \
                             tk_xxx_yy,  \
                             tk_xxx_yz,  \
                             tk_xxx_z,   \
                             tk_xxx_zz,  \
                             tk_xxxx_xx, \
                             tk_xxxx_xy, \
                             tk_xxxx_xz, \
                             tk_xxxx_yy, \
                             tk_xxxx_yz, \
                             tk_xxxx_zz, \
                             ts_xx_xx,   \
                             ts_xx_xy,   \
                             ts_xx_xz,   \
                             ts_xx_yy,   \
                             ts_xx_yz,   \
                             ts_xx_zz,   \
                             ts_xxxx_xx, \
                             ts_xxxx_xy, \
                             ts_xxxx_xz, \
                             ts_xxxx_yy, \
                             ts_xxxx_yz, \
                             ts_xxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxx_xx[i] = -6.0 * ts_xx_xx[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xx[i] * fe_0 + 2.0 * tk_xxx_x[i] * fe_0 + tk_xxx_xx[i] * pa_x[i] +
                        2.0 * ts_xxxx_xx[i] * fz_0;

        tk_xxxx_xy[i] =
            -6.0 * ts_xx_xy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xy[i] * fe_0 + tk_xxx_y[i] * fe_0 + tk_xxx_xy[i] * pa_x[i] + 2.0 * ts_xxxx_xy[i] * fz_0;

        tk_xxxx_xz[i] =
            -6.0 * ts_xx_xz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_xz[i] * fe_0 + tk_xxx_z[i] * fe_0 + tk_xxx_xz[i] * pa_x[i] + 2.0 * ts_xxxx_xz[i] * fz_0;

        tk_xxxx_yy[i] = -6.0 * ts_xx_yy[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yy[i] * fe_0 + tk_xxx_yy[i] * pa_x[i] + 2.0 * ts_xxxx_yy[i] * fz_0;

        tk_xxxx_yz[i] = -6.0 * ts_xx_yz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_yz[i] * fe_0 + tk_xxx_yz[i] * pa_x[i] + 2.0 * ts_xxxx_yz[i] * fz_0;

        tk_xxxx_zz[i] = -6.0 * ts_xx_zz[i] * fbe_0 * fz_0 + 3.0 * tk_xx_zz[i] * fe_0 + tk_xxx_zz[i] * pa_x[i] + 2.0 * ts_xxxx_zz[i] * fz_0;
    }

    // Set up 6-12 components of targeted buffer : GD

    auto tk_xxxy_xx = pbuffer.data(idx_kin_gd + 6);

    auto tk_xxxy_xy = pbuffer.data(idx_kin_gd + 7);

    auto tk_xxxy_xz = pbuffer.data(idx_kin_gd + 8);

    auto tk_xxxy_yy = pbuffer.data(idx_kin_gd + 9);

    auto tk_xxxy_yz = pbuffer.data(idx_kin_gd + 10);

    auto tk_xxxy_zz = pbuffer.data(idx_kin_gd + 11);

#pragma omp simd aligned(pa_y,           \
                             tk_xxx_x,   \
                             tk_xxx_xx,  \
                             tk_xxx_xy,  \
                             tk_xxx_xz,  \
                             tk_xxx_y,   \
                             tk_xxx_yy,  \
                             tk_xxx_yz,  \
                             tk_xxx_z,   \
                             tk_xxx_zz,  \
                             tk_xxxy_xx, \
                             tk_xxxy_xy, \
                             tk_xxxy_xz, \
                             tk_xxxy_yy, \
                             tk_xxxy_yz, \
                             tk_xxxy_zz, \
                             ts_xxxy_xx, \
                             ts_xxxy_xy, \
                             ts_xxxy_xz, \
                             ts_xxxy_yy, \
                             ts_xxxy_yz, \
                             ts_xxxy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxy_xx[i] = tk_xxx_xx[i] * pa_y[i] + 2.0 * ts_xxxy_xx[i] * fz_0;

        tk_xxxy_xy[i] = tk_xxx_x[i] * fe_0 + tk_xxx_xy[i] * pa_y[i] + 2.0 * ts_xxxy_xy[i] * fz_0;

        tk_xxxy_xz[i] = tk_xxx_xz[i] * pa_y[i] + 2.0 * ts_xxxy_xz[i] * fz_0;

        tk_xxxy_yy[i] = 2.0 * tk_xxx_y[i] * fe_0 + tk_xxx_yy[i] * pa_y[i] + 2.0 * ts_xxxy_yy[i] * fz_0;

        tk_xxxy_yz[i] = tk_xxx_z[i] * fe_0 + tk_xxx_yz[i] * pa_y[i] + 2.0 * ts_xxxy_yz[i] * fz_0;

        tk_xxxy_zz[i] = tk_xxx_zz[i] * pa_y[i] + 2.0 * ts_xxxy_zz[i] * fz_0;
    }

    // Set up 12-18 components of targeted buffer : GD

    auto tk_xxxz_xx = pbuffer.data(idx_kin_gd + 12);

    auto tk_xxxz_xy = pbuffer.data(idx_kin_gd + 13);

    auto tk_xxxz_xz = pbuffer.data(idx_kin_gd + 14);

    auto tk_xxxz_yy = pbuffer.data(idx_kin_gd + 15);

    auto tk_xxxz_yz = pbuffer.data(idx_kin_gd + 16);

    auto tk_xxxz_zz = pbuffer.data(idx_kin_gd + 17);

#pragma omp simd aligned(pa_z,           \
                             tk_xxx_x,   \
                             tk_xxx_xx,  \
                             tk_xxx_xy,  \
                             tk_xxx_xz,  \
                             tk_xxx_y,   \
                             tk_xxx_yy,  \
                             tk_xxx_yz,  \
                             tk_xxx_z,   \
                             tk_xxx_zz,  \
                             tk_xxxz_xx, \
                             tk_xxxz_xy, \
                             tk_xxxz_xz, \
                             tk_xxxz_yy, \
                             tk_xxxz_yz, \
                             tk_xxxz_zz, \
                             ts_xxxz_xx, \
                             ts_xxxz_xy, \
                             ts_xxxz_xz, \
                             ts_xxxz_yy, \
                             ts_xxxz_yz, \
                             ts_xxxz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxz_xx[i] = tk_xxx_xx[i] * pa_z[i] + 2.0 * ts_xxxz_xx[i] * fz_0;

        tk_xxxz_xy[i] = tk_xxx_xy[i] * pa_z[i] + 2.0 * ts_xxxz_xy[i] * fz_0;

        tk_xxxz_xz[i] = tk_xxx_x[i] * fe_0 + tk_xxx_xz[i] * pa_z[i] + 2.0 * ts_xxxz_xz[i] * fz_0;

        tk_xxxz_yy[i] = tk_xxx_yy[i] * pa_z[i] + 2.0 * ts_xxxz_yy[i] * fz_0;

        tk_xxxz_yz[i] = tk_xxx_y[i] * fe_0 + tk_xxx_yz[i] * pa_z[i] + 2.0 * ts_xxxz_yz[i] * fz_0;

        tk_xxxz_zz[i] = 2.0 * tk_xxx_z[i] * fe_0 + tk_xxx_zz[i] * pa_z[i] + 2.0 * ts_xxxz_zz[i] * fz_0;
    }

    // Set up 18-24 components of targeted buffer : GD

    auto tk_xxyy_xx = pbuffer.data(idx_kin_gd + 18);

    auto tk_xxyy_xy = pbuffer.data(idx_kin_gd + 19);

    auto tk_xxyy_xz = pbuffer.data(idx_kin_gd + 20);

    auto tk_xxyy_yy = pbuffer.data(idx_kin_gd + 21);

    auto tk_xxyy_yz = pbuffer.data(idx_kin_gd + 22);

    auto tk_xxyy_zz = pbuffer.data(idx_kin_gd + 23);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tk_xx_xx,   \
                             tk_xx_xz,   \
                             tk_xxy_xx,  \
                             tk_xxy_xz,  \
                             tk_xxyy_xx, \
                             tk_xxyy_xy, \
                             tk_xxyy_xz, \
                             tk_xxyy_yy, \
                             tk_xxyy_yz, \
                             tk_xxyy_zz, \
                             tk_xyy_xy,  \
                             tk_xyy_y,   \
                             tk_xyy_yy,  \
                             tk_xyy_yz,  \
                             tk_xyy_zz,  \
                             tk_yy_xy,   \
                             tk_yy_yy,   \
                             tk_yy_yz,   \
                             tk_yy_zz,   \
                             ts_xx_xx,   \
                             ts_xx_xz,   \
                             ts_xxyy_xx, \
                             ts_xxyy_xy, \
                             ts_xxyy_xz, \
                             ts_xxyy_yy, \
                             ts_xxyy_yz, \
                             ts_xxyy_zz, \
                             ts_yy_xy,   \
                             ts_yy_yy,   \
                             ts_yy_yz,   \
                             ts_yy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyy_xx[i] = -2.0 * ts_xx_xx[i] * fbe_0 * fz_0 + tk_xx_xx[i] * fe_0 + tk_xxy_xx[i] * pa_y[i] + 2.0 * ts_xxyy_xx[i] * fz_0;

        tk_xxyy_xy[i] =
            -2.0 * ts_yy_xy[i] * fbe_0 * fz_0 + tk_yy_xy[i] * fe_0 + tk_xyy_y[i] * fe_0 + tk_xyy_xy[i] * pa_x[i] + 2.0 * ts_xxyy_xy[i] * fz_0;

        tk_xxyy_xz[i] = -2.0 * ts_xx_xz[i] * fbe_0 * fz_0 + tk_xx_xz[i] * fe_0 + tk_xxy_xz[i] * pa_y[i] + 2.0 * ts_xxyy_xz[i] * fz_0;

        tk_xxyy_yy[i] = -2.0 * ts_yy_yy[i] * fbe_0 * fz_0 + tk_yy_yy[i] * fe_0 + tk_xyy_yy[i] * pa_x[i] + 2.0 * ts_xxyy_yy[i] * fz_0;

        tk_xxyy_yz[i] = -2.0 * ts_yy_yz[i] * fbe_0 * fz_0 + tk_yy_yz[i] * fe_0 + tk_xyy_yz[i] * pa_x[i] + 2.0 * ts_xxyy_yz[i] * fz_0;

        tk_xxyy_zz[i] = -2.0 * ts_yy_zz[i] * fbe_0 * fz_0 + tk_yy_zz[i] * fe_0 + tk_xyy_zz[i] * pa_x[i] + 2.0 * ts_xxyy_zz[i] * fz_0;
    }

    // Set up 24-30 components of targeted buffer : GD

    auto tk_xxyz_xx = pbuffer.data(idx_kin_gd + 24);

    auto tk_xxyz_xy = pbuffer.data(idx_kin_gd + 25);

    auto tk_xxyz_xz = pbuffer.data(idx_kin_gd + 26);

    auto tk_xxyz_yy = pbuffer.data(idx_kin_gd + 27);

    auto tk_xxyz_yz = pbuffer.data(idx_kin_gd + 28);

    auto tk_xxyz_zz = pbuffer.data(idx_kin_gd + 29);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tk_xxy_xy,  \
                             tk_xxy_yy,  \
                             tk_xxyz_xx, \
                             tk_xxyz_xy, \
                             tk_xxyz_xz, \
                             tk_xxyz_yy, \
                             tk_xxyz_yz, \
                             tk_xxyz_zz, \
                             tk_xxz_xx,  \
                             tk_xxz_xz,  \
                             tk_xxz_yz,  \
                             tk_xxz_z,   \
                             tk_xxz_zz,  \
                             ts_xxyz_xx, \
                             ts_xxyz_xy, \
                             ts_xxyz_xz, \
                             ts_xxyz_yy, \
                             ts_xxyz_yz, \
                             ts_xxyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyz_xx[i] = tk_xxz_xx[i] * pa_y[i] + 2.0 * ts_xxyz_xx[i] * fz_0;

        tk_xxyz_xy[i] = tk_xxy_xy[i] * pa_z[i] + 2.0 * ts_xxyz_xy[i] * fz_0;

        tk_xxyz_xz[i] = tk_xxz_xz[i] * pa_y[i] + 2.0 * ts_xxyz_xz[i] * fz_0;

        tk_xxyz_yy[i] = tk_xxy_yy[i] * pa_z[i] + 2.0 * ts_xxyz_yy[i] * fz_0;

        tk_xxyz_yz[i] = tk_xxz_z[i] * fe_0 + tk_xxz_yz[i] * pa_y[i] + 2.0 * ts_xxyz_yz[i] * fz_0;

        tk_xxyz_zz[i] = tk_xxz_zz[i] * pa_y[i] + 2.0 * ts_xxyz_zz[i] * fz_0;
    }

    // Set up 30-36 components of targeted buffer : GD

    auto tk_xxzz_xx = pbuffer.data(idx_kin_gd + 30);

    auto tk_xxzz_xy = pbuffer.data(idx_kin_gd + 31);

    auto tk_xxzz_xz = pbuffer.data(idx_kin_gd + 32);

    auto tk_xxzz_yy = pbuffer.data(idx_kin_gd + 33);

    auto tk_xxzz_yz = pbuffer.data(idx_kin_gd + 34);

    auto tk_xxzz_zz = pbuffer.data(idx_kin_gd + 35);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tk_xx_xx,   \
                             tk_xx_xy,   \
                             tk_xxz_xx,  \
                             tk_xxz_xy,  \
                             tk_xxzz_xx, \
                             tk_xxzz_xy, \
                             tk_xxzz_xz, \
                             tk_xxzz_yy, \
                             tk_xxzz_yz, \
                             tk_xxzz_zz, \
                             tk_xzz_xz,  \
                             tk_xzz_yy,  \
                             tk_xzz_yz,  \
                             tk_xzz_z,   \
                             tk_xzz_zz,  \
                             tk_zz_xz,   \
                             tk_zz_yy,   \
                             tk_zz_yz,   \
                             tk_zz_zz,   \
                             ts_xx_xx,   \
                             ts_xx_xy,   \
                             ts_xxzz_xx, \
                             ts_xxzz_xy, \
                             ts_xxzz_xz, \
                             ts_xxzz_yy, \
                             ts_xxzz_yz, \
                             ts_xxzz_zz, \
                             ts_zz_xz,   \
                             ts_zz_yy,   \
                             ts_zz_yz,   \
                             ts_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzz_xx[i] = -2.0 * ts_xx_xx[i] * fbe_0 * fz_0 + tk_xx_xx[i] * fe_0 + tk_xxz_xx[i] * pa_z[i] + 2.0 * ts_xxzz_xx[i] * fz_0;

        tk_xxzz_xy[i] = -2.0 * ts_xx_xy[i] * fbe_0 * fz_0 + tk_xx_xy[i] * fe_0 + tk_xxz_xy[i] * pa_z[i] + 2.0 * ts_xxzz_xy[i] * fz_0;

        tk_xxzz_xz[i] =
            -2.0 * ts_zz_xz[i] * fbe_0 * fz_0 + tk_zz_xz[i] * fe_0 + tk_xzz_z[i] * fe_0 + tk_xzz_xz[i] * pa_x[i] + 2.0 * ts_xxzz_xz[i] * fz_0;

        tk_xxzz_yy[i] = -2.0 * ts_zz_yy[i] * fbe_0 * fz_0 + tk_zz_yy[i] * fe_0 + tk_xzz_yy[i] * pa_x[i] + 2.0 * ts_xxzz_yy[i] * fz_0;

        tk_xxzz_yz[i] = -2.0 * ts_zz_yz[i] * fbe_0 * fz_0 + tk_zz_yz[i] * fe_0 + tk_xzz_yz[i] * pa_x[i] + 2.0 * ts_xxzz_yz[i] * fz_0;

        tk_xxzz_zz[i] = -2.0 * ts_zz_zz[i] * fbe_0 * fz_0 + tk_zz_zz[i] * fe_0 + tk_xzz_zz[i] * pa_x[i] + 2.0 * ts_xxzz_zz[i] * fz_0;
    }

    // Set up 36-42 components of targeted buffer : GD

    auto tk_xyyy_xx = pbuffer.data(idx_kin_gd + 36);

    auto tk_xyyy_xy = pbuffer.data(idx_kin_gd + 37);

    auto tk_xyyy_xz = pbuffer.data(idx_kin_gd + 38);

    auto tk_xyyy_yy = pbuffer.data(idx_kin_gd + 39);

    auto tk_xyyy_yz = pbuffer.data(idx_kin_gd + 40);

    auto tk_xyyy_zz = pbuffer.data(idx_kin_gd + 41);

#pragma omp simd aligned(pa_x,           \
                             tk_xyyy_xx, \
                             tk_xyyy_xy, \
                             tk_xyyy_xz, \
                             tk_xyyy_yy, \
                             tk_xyyy_yz, \
                             tk_xyyy_zz, \
                             tk_yyy_x,   \
                             tk_yyy_xx,  \
                             tk_yyy_xy,  \
                             tk_yyy_xz,  \
                             tk_yyy_y,   \
                             tk_yyy_yy,  \
                             tk_yyy_yz,  \
                             tk_yyy_z,   \
                             tk_yyy_zz,  \
                             ts_xyyy_xx, \
                             ts_xyyy_xy, \
                             ts_xyyy_xz, \
                             ts_xyyy_yy, \
                             ts_xyyy_yz, \
                             ts_xyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyy_xx[i] = 2.0 * tk_yyy_x[i] * fe_0 + tk_yyy_xx[i] * pa_x[i] + 2.0 * ts_xyyy_xx[i] * fz_0;

        tk_xyyy_xy[i] = tk_yyy_y[i] * fe_0 + tk_yyy_xy[i] * pa_x[i] + 2.0 * ts_xyyy_xy[i] * fz_0;

        tk_xyyy_xz[i] = tk_yyy_z[i] * fe_0 + tk_yyy_xz[i] * pa_x[i] + 2.0 * ts_xyyy_xz[i] * fz_0;

        tk_xyyy_yy[i] = tk_yyy_yy[i] * pa_x[i] + 2.0 * ts_xyyy_yy[i] * fz_0;

        tk_xyyy_yz[i] = tk_yyy_yz[i] * pa_x[i] + 2.0 * ts_xyyy_yz[i] * fz_0;

        tk_xyyy_zz[i] = tk_yyy_zz[i] * pa_x[i] + 2.0 * ts_xyyy_zz[i] * fz_0;
    }

    // Set up 42-48 components of targeted buffer : GD

    auto tk_xyyz_xx = pbuffer.data(idx_kin_gd + 42);

    auto tk_xyyz_xy = pbuffer.data(idx_kin_gd + 43);

    auto tk_xyyz_xz = pbuffer.data(idx_kin_gd + 44);

    auto tk_xyyz_yy = pbuffer.data(idx_kin_gd + 45);

    auto tk_xyyz_yz = pbuffer.data(idx_kin_gd + 46);

    auto tk_xyyz_zz = pbuffer.data(idx_kin_gd + 47);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tk_xyy_xx,  \
                             tk_xyy_xy,  \
                             tk_xyyz_xx, \
                             tk_xyyz_xy, \
                             tk_xyyz_xz, \
                             tk_xyyz_yy, \
                             tk_xyyz_yz, \
                             tk_xyyz_zz, \
                             tk_yyz_xz,  \
                             tk_yyz_yy,  \
                             tk_yyz_yz,  \
                             tk_yyz_z,   \
                             tk_yyz_zz,  \
                             ts_xyyz_xx, \
                             ts_xyyz_xy, \
                             ts_xyyz_xz, \
                             ts_xyyz_yy, \
                             ts_xyyz_yz, \
                             ts_xyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyz_xx[i] = tk_xyy_xx[i] * pa_z[i] + 2.0 * ts_xyyz_xx[i] * fz_0;

        tk_xyyz_xy[i] = tk_xyy_xy[i] * pa_z[i] + 2.0 * ts_xyyz_xy[i] * fz_0;

        tk_xyyz_xz[i] = tk_yyz_z[i] * fe_0 + tk_yyz_xz[i] * pa_x[i] + 2.0 * ts_xyyz_xz[i] * fz_0;

        tk_xyyz_yy[i] = tk_yyz_yy[i] * pa_x[i] + 2.0 * ts_xyyz_yy[i] * fz_0;

        tk_xyyz_yz[i] = tk_yyz_yz[i] * pa_x[i] + 2.0 * ts_xyyz_yz[i] * fz_0;

        tk_xyyz_zz[i] = tk_yyz_zz[i] * pa_x[i] + 2.0 * ts_xyyz_zz[i] * fz_0;
    }

    // Set up 48-54 components of targeted buffer : GD

    auto tk_xyzz_xx = pbuffer.data(idx_kin_gd + 48);

    auto tk_xyzz_xy = pbuffer.data(idx_kin_gd + 49);

    auto tk_xyzz_xz = pbuffer.data(idx_kin_gd + 50);

    auto tk_xyzz_yy = pbuffer.data(idx_kin_gd + 51);

    auto tk_xyzz_yz = pbuffer.data(idx_kin_gd + 52);

    auto tk_xyzz_zz = pbuffer.data(idx_kin_gd + 53);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tk_xyzz_xx, \
                             tk_xyzz_xy, \
                             tk_xyzz_xz, \
                             tk_xyzz_yy, \
                             tk_xyzz_yz, \
                             tk_xyzz_zz, \
                             tk_xzz_xx,  \
                             tk_xzz_xz,  \
                             tk_yzz_xy,  \
                             tk_yzz_y,   \
                             tk_yzz_yy,  \
                             tk_yzz_yz,  \
                             tk_yzz_zz,  \
                             ts_xyzz_xx, \
                             ts_xyzz_xy, \
                             ts_xyzz_xz, \
                             ts_xyzz_yy, \
                             ts_xyzz_yz, \
                             ts_xyzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzz_xx[i] = tk_xzz_xx[i] * pa_y[i] + 2.0 * ts_xyzz_xx[i] * fz_0;

        tk_xyzz_xy[i] = tk_yzz_y[i] * fe_0 + tk_yzz_xy[i] * pa_x[i] + 2.0 * ts_xyzz_xy[i] * fz_0;

        tk_xyzz_xz[i] = tk_xzz_xz[i] * pa_y[i] + 2.0 * ts_xyzz_xz[i] * fz_0;

        tk_xyzz_yy[i] = tk_yzz_yy[i] * pa_x[i] + 2.0 * ts_xyzz_yy[i] * fz_0;

        tk_xyzz_yz[i] = tk_yzz_yz[i] * pa_x[i] + 2.0 * ts_xyzz_yz[i] * fz_0;

        tk_xyzz_zz[i] = tk_yzz_zz[i] * pa_x[i] + 2.0 * ts_xyzz_zz[i] * fz_0;
    }

    // Set up 54-60 components of targeted buffer : GD

    auto tk_xzzz_xx = pbuffer.data(idx_kin_gd + 54);

    auto tk_xzzz_xy = pbuffer.data(idx_kin_gd + 55);

    auto tk_xzzz_xz = pbuffer.data(idx_kin_gd + 56);

    auto tk_xzzz_yy = pbuffer.data(idx_kin_gd + 57);

    auto tk_xzzz_yz = pbuffer.data(idx_kin_gd + 58);

    auto tk_xzzz_zz = pbuffer.data(idx_kin_gd + 59);

#pragma omp simd aligned(pa_x,           \
                             tk_xzzz_xx, \
                             tk_xzzz_xy, \
                             tk_xzzz_xz, \
                             tk_xzzz_yy, \
                             tk_xzzz_yz, \
                             tk_xzzz_zz, \
                             tk_zzz_x,   \
                             tk_zzz_xx,  \
                             tk_zzz_xy,  \
                             tk_zzz_xz,  \
                             tk_zzz_y,   \
                             tk_zzz_yy,  \
                             tk_zzz_yz,  \
                             tk_zzz_z,   \
                             tk_zzz_zz,  \
                             ts_xzzz_xx, \
                             ts_xzzz_xy, \
                             ts_xzzz_xz, \
                             ts_xzzz_yy, \
                             ts_xzzz_yz, \
                             ts_xzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzz_xx[i] = 2.0 * tk_zzz_x[i] * fe_0 + tk_zzz_xx[i] * pa_x[i] + 2.0 * ts_xzzz_xx[i] * fz_0;

        tk_xzzz_xy[i] = tk_zzz_y[i] * fe_0 + tk_zzz_xy[i] * pa_x[i] + 2.0 * ts_xzzz_xy[i] * fz_0;

        tk_xzzz_xz[i] = tk_zzz_z[i] * fe_0 + tk_zzz_xz[i] * pa_x[i] + 2.0 * ts_xzzz_xz[i] * fz_0;

        tk_xzzz_yy[i] = tk_zzz_yy[i] * pa_x[i] + 2.0 * ts_xzzz_yy[i] * fz_0;

        tk_xzzz_yz[i] = tk_zzz_yz[i] * pa_x[i] + 2.0 * ts_xzzz_yz[i] * fz_0;

        tk_xzzz_zz[i] = tk_zzz_zz[i] * pa_x[i] + 2.0 * ts_xzzz_zz[i] * fz_0;
    }

    // Set up 60-66 components of targeted buffer : GD

    auto tk_yyyy_xx = pbuffer.data(idx_kin_gd + 60);

    auto tk_yyyy_xy = pbuffer.data(idx_kin_gd + 61);

    auto tk_yyyy_xz = pbuffer.data(idx_kin_gd + 62);

    auto tk_yyyy_yy = pbuffer.data(idx_kin_gd + 63);

    auto tk_yyyy_yz = pbuffer.data(idx_kin_gd + 64);

    auto tk_yyyy_zz = pbuffer.data(idx_kin_gd + 65);

#pragma omp simd aligned(pa_y,           \
                             tk_yy_xx,   \
                             tk_yy_xy,   \
                             tk_yy_xz,   \
                             tk_yy_yy,   \
                             tk_yy_yz,   \
                             tk_yy_zz,   \
                             tk_yyy_x,   \
                             tk_yyy_xx,  \
                             tk_yyy_xy,  \
                             tk_yyy_xz,  \
                             tk_yyy_y,   \
                             tk_yyy_yy,  \
                             tk_yyy_yz,  \
                             tk_yyy_z,   \
                             tk_yyy_zz,  \
                             tk_yyyy_xx, \
                             tk_yyyy_xy, \
                             tk_yyyy_xz, \
                             tk_yyyy_yy, \
                             tk_yyyy_yz, \
                             tk_yyyy_zz, \
                             ts_yy_xx,   \
                             ts_yy_xy,   \
                             ts_yy_xz,   \
                             ts_yy_yy,   \
                             ts_yy_yz,   \
                             ts_yy_zz,   \
                             ts_yyyy_xx, \
                             ts_yyyy_xy, \
                             ts_yyyy_xz, \
                             ts_yyyy_yy, \
                             ts_yyyy_yz, \
                             ts_yyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyy_xx[i] = -6.0 * ts_yy_xx[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xx[i] * fe_0 + tk_yyy_xx[i] * pa_y[i] + 2.0 * ts_yyyy_xx[i] * fz_0;

        tk_yyyy_xy[i] =
            -6.0 * ts_yy_xy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xy[i] * fe_0 + tk_yyy_x[i] * fe_0 + tk_yyy_xy[i] * pa_y[i] + 2.0 * ts_yyyy_xy[i] * fz_0;

        tk_yyyy_xz[i] = -6.0 * ts_yy_xz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_xz[i] * fe_0 + tk_yyy_xz[i] * pa_y[i] + 2.0 * ts_yyyy_xz[i] * fz_0;

        tk_yyyy_yy[i] = -6.0 * ts_yy_yy[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yy[i] * fe_0 + 2.0 * tk_yyy_y[i] * fe_0 + tk_yyy_yy[i] * pa_y[i] +
                        2.0 * ts_yyyy_yy[i] * fz_0;

        tk_yyyy_yz[i] =
            -6.0 * ts_yy_yz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_yz[i] * fe_0 + tk_yyy_z[i] * fe_0 + tk_yyy_yz[i] * pa_y[i] + 2.0 * ts_yyyy_yz[i] * fz_0;

        tk_yyyy_zz[i] = -6.0 * ts_yy_zz[i] * fbe_0 * fz_0 + 3.0 * tk_yy_zz[i] * fe_0 + tk_yyy_zz[i] * pa_y[i] + 2.0 * ts_yyyy_zz[i] * fz_0;
    }

    // Set up 66-72 components of targeted buffer : GD

    auto tk_yyyz_xx = pbuffer.data(idx_kin_gd + 66);

    auto tk_yyyz_xy = pbuffer.data(idx_kin_gd + 67);

    auto tk_yyyz_xz = pbuffer.data(idx_kin_gd + 68);

    auto tk_yyyz_yy = pbuffer.data(idx_kin_gd + 69);

    auto tk_yyyz_yz = pbuffer.data(idx_kin_gd + 70);

    auto tk_yyyz_zz = pbuffer.data(idx_kin_gd + 71);

#pragma omp simd aligned(pa_z,           \
                             tk_yyy_x,   \
                             tk_yyy_xx,  \
                             tk_yyy_xy,  \
                             tk_yyy_xz,  \
                             tk_yyy_y,   \
                             tk_yyy_yy,  \
                             tk_yyy_yz,  \
                             tk_yyy_z,   \
                             tk_yyy_zz,  \
                             tk_yyyz_xx, \
                             tk_yyyz_xy, \
                             tk_yyyz_xz, \
                             tk_yyyz_yy, \
                             tk_yyyz_yz, \
                             tk_yyyz_zz, \
                             ts_yyyz_xx, \
                             ts_yyyz_xy, \
                             ts_yyyz_xz, \
                             ts_yyyz_yy, \
                             ts_yyyz_yz, \
                             ts_yyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyz_xx[i] = tk_yyy_xx[i] * pa_z[i] + 2.0 * ts_yyyz_xx[i] * fz_0;

        tk_yyyz_xy[i] = tk_yyy_xy[i] * pa_z[i] + 2.0 * ts_yyyz_xy[i] * fz_0;

        tk_yyyz_xz[i] = tk_yyy_x[i] * fe_0 + tk_yyy_xz[i] * pa_z[i] + 2.0 * ts_yyyz_xz[i] * fz_0;

        tk_yyyz_yy[i] = tk_yyy_yy[i] * pa_z[i] + 2.0 * ts_yyyz_yy[i] * fz_0;

        tk_yyyz_yz[i] = tk_yyy_y[i] * fe_0 + tk_yyy_yz[i] * pa_z[i] + 2.0 * ts_yyyz_yz[i] * fz_0;

        tk_yyyz_zz[i] = 2.0 * tk_yyy_z[i] * fe_0 + tk_yyy_zz[i] * pa_z[i] + 2.0 * ts_yyyz_zz[i] * fz_0;
    }

    // Set up 72-78 components of targeted buffer : GD

    auto tk_yyzz_xx = pbuffer.data(idx_kin_gd + 72);

    auto tk_yyzz_xy = pbuffer.data(idx_kin_gd + 73);

    auto tk_yyzz_xz = pbuffer.data(idx_kin_gd + 74);

    auto tk_yyzz_yy = pbuffer.data(idx_kin_gd + 75);

    auto tk_yyzz_yz = pbuffer.data(idx_kin_gd + 76);

    auto tk_yyzz_zz = pbuffer.data(idx_kin_gd + 77);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tk_yy_xy,   \
                             tk_yy_yy,   \
                             tk_yyz_xy,  \
                             tk_yyz_yy,  \
                             tk_yyzz_xx, \
                             tk_yyzz_xy, \
                             tk_yyzz_xz, \
                             tk_yyzz_yy, \
                             tk_yyzz_yz, \
                             tk_yyzz_zz, \
                             tk_yzz_xx,  \
                             tk_yzz_xz,  \
                             tk_yzz_yz,  \
                             tk_yzz_z,   \
                             tk_yzz_zz,  \
                             tk_zz_xx,   \
                             tk_zz_xz,   \
                             tk_zz_yz,   \
                             tk_zz_zz,   \
                             ts_yy_xy,   \
                             ts_yy_yy,   \
                             ts_yyzz_xx, \
                             ts_yyzz_xy, \
                             ts_yyzz_xz, \
                             ts_yyzz_yy, \
                             ts_yyzz_yz, \
                             ts_yyzz_zz, \
                             ts_zz_xx,   \
                             ts_zz_xz,   \
                             ts_zz_yz,   \
                             ts_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzz_xx[i] = -2.0 * ts_zz_xx[i] * fbe_0 * fz_0 + tk_zz_xx[i] * fe_0 + tk_yzz_xx[i] * pa_y[i] + 2.0 * ts_yyzz_xx[i] * fz_0;

        tk_yyzz_xy[i] = -2.0 * ts_yy_xy[i] * fbe_0 * fz_0 + tk_yy_xy[i] * fe_0 + tk_yyz_xy[i] * pa_z[i] + 2.0 * ts_yyzz_xy[i] * fz_0;

        tk_yyzz_xz[i] = -2.0 * ts_zz_xz[i] * fbe_0 * fz_0 + tk_zz_xz[i] * fe_0 + tk_yzz_xz[i] * pa_y[i] + 2.0 * ts_yyzz_xz[i] * fz_0;

        tk_yyzz_yy[i] = -2.0 * ts_yy_yy[i] * fbe_0 * fz_0 + tk_yy_yy[i] * fe_0 + tk_yyz_yy[i] * pa_z[i] + 2.0 * ts_yyzz_yy[i] * fz_0;

        tk_yyzz_yz[i] =
            -2.0 * ts_zz_yz[i] * fbe_0 * fz_0 + tk_zz_yz[i] * fe_0 + tk_yzz_z[i] * fe_0 + tk_yzz_yz[i] * pa_y[i] + 2.0 * ts_yyzz_yz[i] * fz_0;

        tk_yyzz_zz[i] = -2.0 * ts_zz_zz[i] * fbe_0 * fz_0 + tk_zz_zz[i] * fe_0 + tk_yzz_zz[i] * pa_y[i] + 2.0 * ts_yyzz_zz[i] * fz_0;
    }

    // Set up 78-84 components of targeted buffer : GD

    auto tk_yzzz_xx = pbuffer.data(idx_kin_gd + 78);

    auto tk_yzzz_xy = pbuffer.data(idx_kin_gd + 79);

    auto tk_yzzz_xz = pbuffer.data(idx_kin_gd + 80);

    auto tk_yzzz_yy = pbuffer.data(idx_kin_gd + 81);

    auto tk_yzzz_yz = pbuffer.data(idx_kin_gd + 82);

    auto tk_yzzz_zz = pbuffer.data(idx_kin_gd + 83);

#pragma omp simd aligned(pa_y,           \
                             tk_yzzz_xx, \
                             tk_yzzz_xy, \
                             tk_yzzz_xz, \
                             tk_yzzz_yy, \
                             tk_yzzz_yz, \
                             tk_yzzz_zz, \
                             tk_zzz_x,   \
                             tk_zzz_xx,  \
                             tk_zzz_xy,  \
                             tk_zzz_xz,  \
                             tk_zzz_y,   \
                             tk_zzz_yy,  \
                             tk_zzz_yz,  \
                             tk_zzz_z,   \
                             tk_zzz_zz,  \
                             ts_yzzz_xx, \
                             ts_yzzz_xy, \
                             ts_yzzz_xz, \
                             ts_yzzz_yy, \
                             ts_yzzz_yz, \
                             ts_yzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzz_xx[i] = tk_zzz_xx[i] * pa_y[i] + 2.0 * ts_yzzz_xx[i] * fz_0;

        tk_yzzz_xy[i] = tk_zzz_x[i] * fe_0 + tk_zzz_xy[i] * pa_y[i] + 2.0 * ts_yzzz_xy[i] * fz_0;

        tk_yzzz_xz[i] = tk_zzz_xz[i] * pa_y[i] + 2.0 * ts_yzzz_xz[i] * fz_0;

        tk_yzzz_yy[i] = 2.0 * tk_zzz_y[i] * fe_0 + tk_zzz_yy[i] * pa_y[i] + 2.0 * ts_yzzz_yy[i] * fz_0;

        tk_yzzz_yz[i] = tk_zzz_z[i] * fe_0 + tk_zzz_yz[i] * pa_y[i] + 2.0 * ts_yzzz_yz[i] * fz_0;

        tk_yzzz_zz[i] = tk_zzz_zz[i] * pa_y[i] + 2.0 * ts_yzzz_zz[i] * fz_0;
    }

    // Set up 84-90 components of targeted buffer : GD

    auto tk_zzzz_xx = pbuffer.data(idx_kin_gd + 84);

    auto tk_zzzz_xy = pbuffer.data(idx_kin_gd + 85);

    auto tk_zzzz_xz = pbuffer.data(idx_kin_gd + 86);

    auto tk_zzzz_yy = pbuffer.data(idx_kin_gd + 87);

    auto tk_zzzz_yz = pbuffer.data(idx_kin_gd + 88);

    auto tk_zzzz_zz = pbuffer.data(idx_kin_gd + 89);

#pragma omp simd aligned(pa_z,           \
                             tk_zz_xx,   \
                             tk_zz_xy,   \
                             tk_zz_xz,   \
                             tk_zz_yy,   \
                             tk_zz_yz,   \
                             tk_zz_zz,   \
                             tk_zzz_x,   \
                             tk_zzz_xx,  \
                             tk_zzz_xy,  \
                             tk_zzz_xz,  \
                             tk_zzz_y,   \
                             tk_zzz_yy,  \
                             tk_zzz_yz,  \
                             tk_zzz_z,   \
                             tk_zzz_zz,  \
                             tk_zzzz_xx, \
                             tk_zzzz_xy, \
                             tk_zzzz_xz, \
                             tk_zzzz_yy, \
                             tk_zzzz_yz, \
                             tk_zzzz_zz, \
                             ts_zz_xx,   \
                             ts_zz_xy,   \
                             ts_zz_xz,   \
                             ts_zz_yy,   \
                             ts_zz_yz,   \
                             ts_zz_zz,   \
                             ts_zzzz_xx, \
                             ts_zzzz_xy, \
                             ts_zzzz_xz, \
                             ts_zzzz_yy, \
                             ts_zzzz_yz, \
                             ts_zzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzz_xx[i] = -6.0 * ts_zz_xx[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xx[i] * fe_0 + tk_zzz_xx[i] * pa_z[i] + 2.0 * ts_zzzz_xx[i] * fz_0;

        tk_zzzz_xy[i] = -6.0 * ts_zz_xy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xy[i] * fe_0 + tk_zzz_xy[i] * pa_z[i] + 2.0 * ts_zzzz_xy[i] * fz_0;

        tk_zzzz_xz[i] =
            -6.0 * ts_zz_xz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_xz[i] * fe_0 + tk_zzz_x[i] * fe_0 + tk_zzz_xz[i] * pa_z[i] + 2.0 * ts_zzzz_xz[i] * fz_0;

        tk_zzzz_yy[i] = -6.0 * ts_zz_yy[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yy[i] * fe_0 + tk_zzz_yy[i] * pa_z[i] + 2.0 * ts_zzzz_yy[i] * fz_0;

        tk_zzzz_yz[i] =
            -6.0 * ts_zz_yz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_yz[i] * fe_0 + tk_zzz_y[i] * fe_0 + tk_zzz_yz[i] * pa_z[i] + 2.0 * ts_zzzz_yz[i] * fz_0;

        tk_zzzz_zz[i] = -6.0 * ts_zz_zz[i] * fbe_0 * fz_0 + 3.0 * tk_zz_zz[i] * fe_0 + 2.0 * tk_zzz_z[i] * fe_0 + tk_zzz_zz[i] * pa_z[i] +
                        2.0 * ts_zzzz_zz[i] * fz_0;
    }
}

}  // namespace kinrec
