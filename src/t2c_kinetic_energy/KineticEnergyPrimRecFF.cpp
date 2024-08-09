#include "KineticEnergyPrimRecFF.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_ff(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_ff,
                            const size_t              idx_ovl_pf,
                            const size_t              idx_kin_pf,
                            const size_t              idx_kin_dd,
                            const size_t              idx_kin_df,
                            const size_t              idx_ovl_ff,
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

    // Set up components of auxiliary buffer : PF

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

    auto tk_yz_yz = pbuffer.data(idx_kin_dd + 28);

    auto tk_zz_xx = pbuffer.data(idx_kin_dd + 30);

    auto tk_zz_xy = pbuffer.data(idx_kin_dd + 31);

    auto tk_zz_xz = pbuffer.data(idx_kin_dd + 32);

    auto tk_zz_yy = pbuffer.data(idx_kin_dd + 33);

    auto tk_zz_yz = pbuffer.data(idx_kin_dd + 34);

    auto tk_zz_zz = pbuffer.data(idx_kin_dd + 35);

    // Set up components of auxiliary buffer : DF

    auto tk_xx_xxx = pbuffer.data(idx_kin_df);

    auto tk_xx_xxy = pbuffer.data(idx_kin_df + 1);

    auto tk_xx_xxz = pbuffer.data(idx_kin_df + 2);

    auto tk_xx_xyy = pbuffer.data(idx_kin_df + 3);

    auto tk_xx_xyz = pbuffer.data(idx_kin_df + 4);

    auto tk_xx_xzz = pbuffer.data(idx_kin_df + 5);

    auto tk_xx_yyy = pbuffer.data(idx_kin_df + 6);

    auto tk_xx_yyz = pbuffer.data(idx_kin_df + 7);

    auto tk_xx_yzz = pbuffer.data(idx_kin_df + 8);

    auto tk_xx_zzz = pbuffer.data(idx_kin_df + 9);

    auto tk_xy_xxy = pbuffer.data(idx_kin_df + 11);

    auto tk_xy_xyy = pbuffer.data(idx_kin_df + 13);

    auto tk_xz_xxx = pbuffer.data(idx_kin_df + 20);

    auto tk_xz_xxz = pbuffer.data(idx_kin_df + 22);

    auto tk_xz_xzz = pbuffer.data(idx_kin_df + 25);

    auto tk_yy_xxx = pbuffer.data(idx_kin_df + 30);

    auto tk_yy_xxy = pbuffer.data(idx_kin_df + 31);

    auto tk_yy_xxz = pbuffer.data(idx_kin_df + 32);

    auto tk_yy_xyy = pbuffer.data(idx_kin_df + 33);

    auto tk_yy_xyz = pbuffer.data(idx_kin_df + 34);

    auto tk_yy_xzz = pbuffer.data(idx_kin_df + 35);

    auto tk_yy_yyy = pbuffer.data(idx_kin_df + 36);

    auto tk_yy_yyz = pbuffer.data(idx_kin_df + 37);

    auto tk_yy_yzz = pbuffer.data(idx_kin_df + 38);

    auto tk_yy_zzz = pbuffer.data(idx_kin_df + 39);

    auto tk_yz_xyz = pbuffer.data(idx_kin_df + 44);

    auto tk_yz_yyy = pbuffer.data(idx_kin_df + 46);

    auto tk_yz_yyz = pbuffer.data(idx_kin_df + 47);

    auto tk_yz_yzz = pbuffer.data(idx_kin_df + 48);

    auto tk_yz_zzz = pbuffer.data(idx_kin_df + 49);

    auto tk_zz_xxx = pbuffer.data(idx_kin_df + 50);

    auto tk_zz_xxy = pbuffer.data(idx_kin_df + 51);

    auto tk_zz_xxz = pbuffer.data(idx_kin_df + 52);

    auto tk_zz_xyy = pbuffer.data(idx_kin_df + 53);

    auto tk_zz_xyz = pbuffer.data(idx_kin_df + 54);

    auto tk_zz_xzz = pbuffer.data(idx_kin_df + 55);

    auto tk_zz_yyy = pbuffer.data(idx_kin_df + 56);

    auto tk_zz_yyz = pbuffer.data(idx_kin_df + 57);

    auto tk_zz_yzz = pbuffer.data(idx_kin_df + 58);

    auto tk_zz_zzz = pbuffer.data(idx_kin_df + 59);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ovl_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ovl_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ovl_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ovl_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ovl_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ovl_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ovl_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ovl_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ovl_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ovl_ff + 9);

    auto ts_xxy_xxx = pbuffer.data(idx_ovl_ff + 10);

    auto ts_xxy_xxy = pbuffer.data(idx_ovl_ff + 11);

    auto ts_xxy_xxz = pbuffer.data(idx_ovl_ff + 12);

    auto ts_xxy_xyy = pbuffer.data(idx_ovl_ff + 13);

    auto ts_xxy_xyz = pbuffer.data(idx_ovl_ff + 14);

    auto ts_xxy_xzz = pbuffer.data(idx_ovl_ff + 15);

    auto ts_xxy_yyy = pbuffer.data(idx_ovl_ff + 16);

    auto ts_xxy_yyz = pbuffer.data(idx_ovl_ff + 17);

    auto ts_xxy_yzz = pbuffer.data(idx_ovl_ff + 18);

    auto ts_xxy_zzz = pbuffer.data(idx_ovl_ff + 19);

    auto ts_xxz_xxx = pbuffer.data(idx_ovl_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ovl_ff + 21);

    auto ts_xxz_xxz = pbuffer.data(idx_ovl_ff + 22);

    auto ts_xxz_xyy = pbuffer.data(idx_ovl_ff + 23);

    auto ts_xxz_xyz = pbuffer.data(idx_ovl_ff + 24);

    auto ts_xxz_xzz = pbuffer.data(idx_ovl_ff + 25);

    auto ts_xxz_yyy = pbuffer.data(idx_ovl_ff + 26);

    auto ts_xxz_yyz = pbuffer.data(idx_ovl_ff + 27);

    auto ts_xxz_yzz = pbuffer.data(idx_ovl_ff + 28);

    auto ts_xxz_zzz = pbuffer.data(idx_ovl_ff + 29);

    auto ts_xyy_xxx = pbuffer.data(idx_ovl_ff + 30);

    auto ts_xyy_xxy = pbuffer.data(idx_ovl_ff + 31);

    auto ts_xyy_xxz = pbuffer.data(idx_ovl_ff + 32);

    auto ts_xyy_xyy = pbuffer.data(idx_ovl_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ovl_ff + 34);

    auto ts_xyy_xzz = pbuffer.data(idx_ovl_ff + 35);

    auto ts_xyy_yyy = pbuffer.data(idx_ovl_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ovl_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ovl_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ovl_ff + 39);

    auto ts_xyz_xxx = pbuffer.data(idx_ovl_ff + 40);

    auto ts_xyz_xxy = pbuffer.data(idx_ovl_ff + 41);

    auto ts_xyz_xxz = pbuffer.data(idx_ovl_ff + 42);

    auto ts_xyz_xyy = pbuffer.data(idx_ovl_ff + 43);

    auto ts_xyz_xyz = pbuffer.data(idx_ovl_ff + 44);

    auto ts_xyz_xzz = pbuffer.data(idx_ovl_ff + 45);

    auto ts_xyz_yyy = pbuffer.data(idx_ovl_ff + 46);

    auto ts_xyz_yyz = pbuffer.data(idx_ovl_ff + 47);

    auto ts_xyz_yzz = pbuffer.data(idx_ovl_ff + 48);

    auto ts_xyz_zzz = pbuffer.data(idx_ovl_ff + 49);

    auto ts_xzz_xxx = pbuffer.data(idx_ovl_ff + 50);

    auto ts_xzz_xxy = pbuffer.data(idx_ovl_ff + 51);

    auto ts_xzz_xxz = pbuffer.data(idx_ovl_ff + 52);

    auto ts_xzz_xyy = pbuffer.data(idx_ovl_ff + 53);

    auto ts_xzz_xyz = pbuffer.data(idx_ovl_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ovl_ff + 55);

    auto ts_xzz_yyy = pbuffer.data(idx_ovl_ff + 56);

    auto ts_xzz_yyz = pbuffer.data(idx_ovl_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ovl_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ovl_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ovl_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ovl_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ovl_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ovl_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ovl_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ovl_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ovl_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ovl_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ovl_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ovl_ff + 69);

    auto ts_yyz_xxx = pbuffer.data(idx_ovl_ff + 70);

    auto ts_yyz_xxy = pbuffer.data(idx_ovl_ff + 71);

    auto ts_yyz_xxz = pbuffer.data(idx_ovl_ff + 72);

    auto ts_yyz_xyy = pbuffer.data(idx_ovl_ff + 73);

    auto ts_yyz_xyz = pbuffer.data(idx_ovl_ff + 74);

    auto ts_yyz_xzz = pbuffer.data(idx_ovl_ff + 75);

    auto ts_yyz_yyy = pbuffer.data(idx_ovl_ff + 76);

    auto ts_yyz_yyz = pbuffer.data(idx_ovl_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ovl_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ovl_ff + 79);

    auto ts_yzz_xxx = pbuffer.data(idx_ovl_ff + 80);

    auto ts_yzz_xxy = pbuffer.data(idx_ovl_ff + 81);

    auto ts_yzz_xxz = pbuffer.data(idx_ovl_ff + 82);

    auto ts_yzz_xyy = pbuffer.data(idx_ovl_ff + 83);

    auto ts_yzz_xyz = pbuffer.data(idx_ovl_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ovl_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ovl_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ovl_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ovl_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ovl_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ovl_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ovl_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ovl_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ovl_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ovl_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ovl_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ovl_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ovl_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ovl_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ovl_ff + 99);

    // Set up 0-10 components of targeted buffer : FF

    auto tk_xxx_xxx = pbuffer.data(idx_kin_ff);

    auto tk_xxx_xxy = pbuffer.data(idx_kin_ff + 1);

    auto tk_xxx_xxz = pbuffer.data(idx_kin_ff + 2);

    auto tk_xxx_xyy = pbuffer.data(idx_kin_ff + 3);

    auto tk_xxx_xyz = pbuffer.data(idx_kin_ff + 4);

    auto tk_xxx_xzz = pbuffer.data(idx_kin_ff + 5);

    auto tk_xxx_yyy = pbuffer.data(idx_kin_ff + 6);

    auto tk_xxx_yyz = pbuffer.data(idx_kin_ff + 7);

    auto tk_xxx_yzz = pbuffer.data(idx_kin_ff + 8);

    auto tk_xxx_zzz = pbuffer.data(idx_kin_ff + 9);

#pragma omp simd aligned(pa_x,           \
                             tk_x_xxx,   \
                             tk_x_xxy,   \
                             tk_x_xxz,   \
                             tk_x_xyy,   \
                             tk_x_xyz,   \
                             tk_x_xzz,   \
                             tk_x_yyy,   \
                             tk_x_yyz,   \
                             tk_x_yzz,   \
                             tk_x_zzz,   \
                             tk_xx_xx,   \
                             tk_xx_xxx,  \
                             tk_xx_xxy,  \
                             tk_xx_xxz,  \
                             tk_xx_xy,   \
                             tk_xx_xyy,  \
                             tk_xx_xyz,  \
                             tk_xx_xz,   \
                             tk_xx_xzz,  \
                             tk_xx_yy,   \
                             tk_xx_yyy,  \
                             tk_xx_yyz,  \
                             tk_xx_yz,   \
                             tk_xx_yzz,  \
                             tk_xx_zz,   \
                             tk_xx_zzz,  \
                             tk_xxx_xxx, \
                             tk_xxx_xxy, \
                             tk_xxx_xxz, \
                             tk_xxx_xyy, \
                             tk_xxx_xyz, \
                             tk_xxx_xzz, \
                             tk_xxx_yyy, \
                             tk_xxx_yyz, \
                             tk_xxx_yzz, \
                             tk_xxx_zzz, \
                             ts_x_xxx,   \
                             ts_x_xxy,   \
                             ts_x_xxz,   \
                             ts_x_xyy,   \
                             ts_x_xyz,   \
                             ts_x_xzz,   \
                             ts_x_yyy,   \
                             ts_x_yyz,   \
                             ts_x_yzz,   \
                             ts_x_zzz,   \
                             ts_xxx_xxx, \
                             ts_xxx_xxy, \
                             ts_xxx_xxz, \
                             ts_xxx_xyy, \
                             ts_xxx_xyz, \
                             ts_xxx_xzz, \
                             ts_xxx_yyy, \
                             ts_xxx_yyz, \
                             ts_xxx_yzz, \
                             ts_xxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_xxx[i] = -4.0 * ts_x_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxx[i] * fe_0 + 3.0 * tk_xx_xx[i] * fe_0 +
                        tk_xx_xxx[i] * pa_x[i] + 2.0 * ts_xxx_xxx[i] * fz_0;

        tk_xxx_xxy[i] = -4.0 * ts_x_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxy[i] * fe_0 + 2.0 * tk_xx_xy[i] * fe_0 +
                        tk_xx_xxy[i] * pa_x[i] + 2.0 * ts_xxx_xxy[i] * fz_0;

        tk_xxx_xxz[i] = -4.0 * ts_x_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxz[i] * fe_0 + 2.0 * tk_xx_xz[i] * fe_0 +
                        tk_xx_xxz[i] * pa_x[i] + 2.0 * ts_xxx_xxz[i] * fz_0;

        tk_xxx_xyy[i] = -4.0 * ts_x_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyy[i] * fe_0 + tk_xx_yy[i] * fe_0 +
                        tk_xx_xyy[i] * pa_x[i] + 2.0 * ts_xxx_xyy[i] * fz_0;

        tk_xxx_xyz[i] = -4.0 * ts_x_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyz[i] * fe_0 + tk_xx_yz[i] * fe_0 +
                        tk_xx_xyz[i] * pa_x[i] + 2.0 * ts_xxx_xyz[i] * fz_0;

        tk_xxx_xzz[i] = -4.0 * ts_x_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xzz[i] * fe_0 + tk_xx_zz[i] * fe_0 +
                        tk_xx_xzz[i] * pa_x[i] + 2.0 * ts_xxx_xzz[i] * fz_0;

        tk_xxx_yyy[i] = -4.0 * ts_x_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyy[i] * fe_0 + tk_xx_yyy[i] * pa_x[i] +
                        2.0 * ts_xxx_yyy[i] * fz_0;

        tk_xxx_yyz[i] = -4.0 * ts_x_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyz[i] * fe_0 + tk_xx_yyz[i] * pa_x[i] +
                        2.0 * ts_xxx_yyz[i] * fz_0;

        tk_xxx_yzz[i] = -4.0 * ts_x_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yzz[i] * fe_0 + tk_xx_yzz[i] * pa_x[i] +
                        2.0 * ts_xxx_yzz[i] * fz_0;

        tk_xxx_zzz[i] = -4.0 * ts_x_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_zzz[i] * fe_0 + tk_xx_zzz[i] * pa_x[i] +
                        2.0 * ts_xxx_zzz[i] * fz_0;
    }

    // Set up 10-20 components of targeted buffer : FF

    auto tk_xxy_xxx = pbuffer.data(idx_kin_ff + 10);

    auto tk_xxy_xxy = pbuffer.data(idx_kin_ff + 11);

    auto tk_xxy_xxz = pbuffer.data(idx_kin_ff + 12);

    auto tk_xxy_xyy = pbuffer.data(idx_kin_ff + 13);

    auto tk_xxy_xyz = pbuffer.data(idx_kin_ff + 14);

    auto tk_xxy_xzz = pbuffer.data(idx_kin_ff + 15);

    auto tk_xxy_yyy = pbuffer.data(idx_kin_ff + 16);

    auto tk_xxy_yyz = pbuffer.data(idx_kin_ff + 17);

    auto tk_xxy_yzz = pbuffer.data(idx_kin_ff + 18);

    auto tk_xxy_zzz = pbuffer.data(idx_kin_ff + 19);

#pragma omp simd aligned(pa_y,           \
                             tk_xx_xx,   \
                             tk_xx_xxx,  \
                             tk_xx_xxy,  \
                             tk_xx_xxz,  \
                             tk_xx_xy,   \
                             tk_xx_xyy,  \
                             tk_xx_xyz,  \
                             tk_xx_xz,   \
                             tk_xx_xzz,  \
                             tk_xx_yy,   \
                             tk_xx_yyy,  \
                             tk_xx_yyz,  \
                             tk_xx_yz,   \
                             tk_xx_yzz,  \
                             tk_xx_zz,   \
                             tk_xx_zzz,  \
                             tk_xxy_xxx, \
                             tk_xxy_xxy, \
                             tk_xxy_xxz, \
                             tk_xxy_xyy, \
                             tk_xxy_xyz, \
                             tk_xxy_xzz, \
                             tk_xxy_yyy, \
                             tk_xxy_yyz, \
                             tk_xxy_yzz, \
                             tk_xxy_zzz, \
                             ts_xxy_xxx, \
                             ts_xxy_xxy, \
                             ts_xxy_xxz, \
                             ts_xxy_xyy, \
                             ts_xxy_xyz, \
                             ts_xxy_xzz, \
                             ts_xxy_yyy, \
                             ts_xxy_yyz, \
                             ts_xxy_yzz, \
                             ts_xxy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxy_xxx[i] = tk_xx_xxx[i] * pa_y[i] + 2.0 * ts_xxy_xxx[i] * fz_0;

        tk_xxy_xxy[i] = tk_xx_xx[i] * fe_0 + tk_xx_xxy[i] * pa_y[i] + 2.0 * ts_xxy_xxy[i] * fz_0;

        tk_xxy_xxz[i] = tk_xx_xxz[i] * pa_y[i] + 2.0 * ts_xxy_xxz[i] * fz_0;

        tk_xxy_xyy[i] = 2.0 * tk_xx_xy[i] * fe_0 + tk_xx_xyy[i] * pa_y[i] + 2.0 * ts_xxy_xyy[i] * fz_0;

        tk_xxy_xyz[i] = tk_xx_xz[i] * fe_0 + tk_xx_xyz[i] * pa_y[i] + 2.0 * ts_xxy_xyz[i] * fz_0;

        tk_xxy_xzz[i] = tk_xx_xzz[i] * pa_y[i] + 2.0 * ts_xxy_xzz[i] * fz_0;

        tk_xxy_yyy[i] = 3.0 * tk_xx_yy[i] * fe_0 + tk_xx_yyy[i] * pa_y[i] + 2.0 * ts_xxy_yyy[i] * fz_0;

        tk_xxy_yyz[i] = 2.0 * tk_xx_yz[i] * fe_0 + tk_xx_yyz[i] * pa_y[i] + 2.0 * ts_xxy_yyz[i] * fz_0;

        tk_xxy_yzz[i] = tk_xx_zz[i] * fe_0 + tk_xx_yzz[i] * pa_y[i] + 2.0 * ts_xxy_yzz[i] * fz_0;

        tk_xxy_zzz[i] = tk_xx_zzz[i] * pa_y[i] + 2.0 * ts_xxy_zzz[i] * fz_0;
    }

    // Set up 20-30 components of targeted buffer : FF

    auto tk_xxz_xxx = pbuffer.data(idx_kin_ff + 20);

    auto tk_xxz_xxy = pbuffer.data(idx_kin_ff + 21);

    auto tk_xxz_xxz = pbuffer.data(idx_kin_ff + 22);

    auto tk_xxz_xyy = pbuffer.data(idx_kin_ff + 23);

    auto tk_xxz_xyz = pbuffer.data(idx_kin_ff + 24);

    auto tk_xxz_xzz = pbuffer.data(idx_kin_ff + 25);

    auto tk_xxz_yyy = pbuffer.data(idx_kin_ff + 26);

    auto tk_xxz_yyz = pbuffer.data(idx_kin_ff + 27);

    auto tk_xxz_yzz = pbuffer.data(idx_kin_ff + 28);

    auto tk_xxz_zzz = pbuffer.data(idx_kin_ff + 29);

#pragma omp simd aligned(pa_z,           \
                             tk_xx_xx,   \
                             tk_xx_xxx,  \
                             tk_xx_xxy,  \
                             tk_xx_xxz,  \
                             tk_xx_xy,   \
                             tk_xx_xyy,  \
                             tk_xx_xyz,  \
                             tk_xx_xz,   \
                             tk_xx_xzz,  \
                             tk_xx_yy,   \
                             tk_xx_yyy,  \
                             tk_xx_yyz,  \
                             tk_xx_yz,   \
                             tk_xx_yzz,  \
                             tk_xx_zz,   \
                             tk_xx_zzz,  \
                             tk_xxz_xxx, \
                             tk_xxz_xxy, \
                             tk_xxz_xxz, \
                             tk_xxz_xyy, \
                             tk_xxz_xyz, \
                             tk_xxz_xzz, \
                             tk_xxz_yyy, \
                             tk_xxz_yyz, \
                             tk_xxz_yzz, \
                             tk_xxz_zzz, \
                             ts_xxz_xxx, \
                             ts_xxz_xxy, \
                             ts_xxz_xxz, \
                             ts_xxz_xyy, \
                             ts_xxz_xyz, \
                             ts_xxz_xzz, \
                             ts_xxz_yyy, \
                             ts_xxz_yyz, \
                             ts_xxz_yzz, \
                             ts_xxz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxz_xxx[i] = tk_xx_xxx[i] * pa_z[i] + 2.0 * ts_xxz_xxx[i] * fz_0;

        tk_xxz_xxy[i] = tk_xx_xxy[i] * pa_z[i] + 2.0 * ts_xxz_xxy[i] * fz_0;

        tk_xxz_xxz[i] = tk_xx_xx[i] * fe_0 + tk_xx_xxz[i] * pa_z[i] + 2.0 * ts_xxz_xxz[i] * fz_0;

        tk_xxz_xyy[i] = tk_xx_xyy[i] * pa_z[i] + 2.0 * ts_xxz_xyy[i] * fz_0;

        tk_xxz_xyz[i] = tk_xx_xy[i] * fe_0 + tk_xx_xyz[i] * pa_z[i] + 2.0 * ts_xxz_xyz[i] * fz_0;

        tk_xxz_xzz[i] = 2.0 * tk_xx_xz[i] * fe_0 + tk_xx_xzz[i] * pa_z[i] + 2.0 * ts_xxz_xzz[i] * fz_0;

        tk_xxz_yyy[i] = tk_xx_yyy[i] * pa_z[i] + 2.0 * ts_xxz_yyy[i] * fz_0;

        tk_xxz_yyz[i] = tk_xx_yy[i] * fe_0 + tk_xx_yyz[i] * pa_z[i] + 2.0 * ts_xxz_yyz[i] * fz_0;

        tk_xxz_yzz[i] = 2.0 * tk_xx_yz[i] * fe_0 + tk_xx_yzz[i] * pa_z[i] + 2.0 * ts_xxz_yzz[i] * fz_0;

        tk_xxz_zzz[i] = 3.0 * tk_xx_zz[i] * fe_0 + tk_xx_zzz[i] * pa_z[i] + 2.0 * ts_xxz_zzz[i] * fz_0;
    }

    // Set up 30-40 components of targeted buffer : FF

    auto tk_xyy_xxx = pbuffer.data(idx_kin_ff + 30);

    auto tk_xyy_xxy = pbuffer.data(idx_kin_ff + 31);

    auto tk_xyy_xxz = pbuffer.data(idx_kin_ff + 32);

    auto tk_xyy_xyy = pbuffer.data(idx_kin_ff + 33);

    auto tk_xyy_xyz = pbuffer.data(idx_kin_ff + 34);

    auto tk_xyy_xzz = pbuffer.data(idx_kin_ff + 35);

    auto tk_xyy_yyy = pbuffer.data(idx_kin_ff + 36);

    auto tk_xyy_yyz = pbuffer.data(idx_kin_ff + 37);

    auto tk_xyy_yzz = pbuffer.data(idx_kin_ff + 38);

    auto tk_xyy_zzz = pbuffer.data(idx_kin_ff + 39);

#pragma omp simd aligned(pa_x,           \
                             tk_xyy_xxx, \
                             tk_xyy_xxy, \
                             tk_xyy_xxz, \
                             tk_xyy_xyy, \
                             tk_xyy_xyz, \
                             tk_xyy_xzz, \
                             tk_xyy_yyy, \
                             tk_xyy_yyz, \
                             tk_xyy_yzz, \
                             tk_xyy_zzz, \
                             tk_yy_xx,   \
                             tk_yy_xxx,  \
                             tk_yy_xxy,  \
                             tk_yy_xxz,  \
                             tk_yy_xy,   \
                             tk_yy_xyy,  \
                             tk_yy_xyz,  \
                             tk_yy_xz,   \
                             tk_yy_xzz,  \
                             tk_yy_yy,   \
                             tk_yy_yyy,  \
                             tk_yy_yyz,  \
                             tk_yy_yz,   \
                             tk_yy_yzz,  \
                             tk_yy_zz,   \
                             tk_yy_zzz,  \
                             ts_xyy_xxx, \
                             ts_xyy_xxy, \
                             ts_xyy_xxz, \
                             ts_xyy_xyy, \
                             ts_xyy_xyz, \
                             ts_xyy_xzz, \
                             ts_xyy_yyy, \
                             ts_xyy_yyz, \
                             ts_xyy_yzz, \
                             ts_xyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyy_xxx[i] = 3.0 * tk_yy_xx[i] * fe_0 + tk_yy_xxx[i] * pa_x[i] + 2.0 * ts_xyy_xxx[i] * fz_0;

        tk_xyy_xxy[i] = 2.0 * tk_yy_xy[i] * fe_0 + tk_yy_xxy[i] * pa_x[i] + 2.0 * ts_xyy_xxy[i] * fz_0;

        tk_xyy_xxz[i] = 2.0 * tk_yy_xz[i] * fe_0 + tk_yy_xxz[i] * pa_x[i] + 2.0 * ts_xyy_xxz[i] * fz_0;

        tk_xyy_xyy[i] = tk_yy_yy[i] * fe_0 + tk_yy_xyy[i] * pa_x[i] + 2.0 * ts_xyy_xyy[i] * fz_0;

        tk_xyy_xyz[i] = tk_yy_yz[i] * fe_0 + tk_yy_xyz[i] * pa_x[i] + 2.0 * ts_xyy_xyz[i] * fz_0;

        tk_xyy_xzz[i] = tk_yy_zz[i] * fe_0 + tk_yy_xzz[i] * pa_x[i] + 2.0 * ts_xyy_xzz[i] * fz_0;

        tk_xyy_yyy[i] = tk_yy_yyy[i] * pa_x[i] + 2.0 * ts_xyy_yyy[i] * fz_0;

        tk_xyy_yyz[i] = tk_yy_yyz[i] * pa_x[i] + 2.0 * ts_xyy_yyz[i] * fz_0;

        tk_xyy_yzz[i] = tk_yy_yzz[i] * pa_x[i] + 2.0 * ts_xyy_yzz[i] * fz_0;

        tk_xyy_zzz[i] = tk_yy_zzz[i] * pa_x[i] + 2.0 * ts_xyy_zzz[i] * fz_0;
    }

    // Set up 40-50 components of targeted buffer : FF

    auto tk_xyz_xxx = pbuffer.data(idx_kin_ff + 40);

    auto tk_xyz_xxy = pbuffer.data(idx_kin_ff + 41);

    auto tk_xyz_xxz = pbuffer.data(idx_kin_ff + 42);

    auto tk_xyz_xyy = pbuffer.data(idx_kin_ff + 43);

    auto tk_xyz_xyz = pbuffer.data(idx_kin_ff + 44);

    auto tk_xyz_xzz = pbuffer.data(idx_kin_ff + 45);

    auto tk_xyz_yyy = pbuffer.data(idx_kin_ff + 46);

    auto tk_xyz_yyz = pbuffer.data(idx_kin_ff + 47);

    auto tk_xyz_yzz = pbuffer.data(idx_kin_ff + 48);

    auto tk_xyz_zzz = pbuffer.data(idx_kin_ff + 49);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pa_z,       \
                             tk_xy_xxy,  \
                             tk_xy_xyy,  \
                             tk_xyz_xxx, \
                             tk_xyz_xxy, \
                             tk_xyz_xxz, \
                             tk_xyz_xyy, \
                             tk_xyz_xyz, \
                             tk_xyz_xzz, \
                             tk_xyz_yyy, \
                             tk_xyz_yyz, \
                             tk_xyz_yzz, \
                             tk_xyz_zzz, \
                             tk_xz_xxx,  \
                             tk_xz_xxz,  \
                             tk_xz_xzz,  \
                             tk_yz_xyz,  \
                             tk_yz_yyy,  \
                             tk_yz_yyz,  \
                             tk_yz_yz,   \
                             tk_yz_yzz,  \
                             tk_yz_zzz,  \
                             ts_xyz_xxx, \
                             ts_xyz_xxy, \
                             ts_xyz_xxz, \
                             ts_xyz_xyy, \
                             ts_xyz_xyz, \
                             ts_xyz_xzz, \
                             ts_xyz_yyy, \
                             ts_xyz_yyz, \
                             ts_xyz_yzz, \
                             ts_xyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyz_xxx[i] = tk_xz_xxx[i] * pa_y[i] + 2.0 * ts_xyz_xxx[i] * fz_0;

        tk_xyz_xxy[i] = tk_xy_xxy[i] * pa_z[i] + 2.0 * ts_xyz_xxy[i] * fz_0;

        tk_xyz_xxz[i] = tk_xz_xxz[i] * pa_y[i] + 2.0 * ts_xyz_xxz[i] * fz_0;

        tk_xyz_xyy[i] = tk_xy_xyy[i] * pa_z[i] + 2.0 * ts_xyz_xyy[i] * fz_0;

        tk_xyz_xyz[i] = tk_yz_yz[i] * fe_0 + tk_yz_xyz[i] * pa_x[i] + 2.0 * ts_xyz_xyz[i] * fz_0;

        tk_xyz_xzz[i] = tk_xz_xzz[i] * pa_y[i] + 2.0 * ts_xyz_xzz[i] * fz_0;

        tk_xyz_yyy[i] = tk_yz_yyy[i] * pa_x[i] + 2.0 * ts_xyz_yyy[i] * fz_0;

        tk_xyz_yyz[i] = tk_yz_yyz[i] * pa_x[i] + 2.0 * ts_xyz_yyz[i] * fz_0;

        tk_xyz_yzz[i] = tk_yz_yzz[i] * pa_x[i] + 2.0 * ts_xyz_yzz[i] * fz_0;

        tk_xyz_zzz[i] = tk_yz_zzz[i] * pa_x[i] + 2.0 * ts_xyz_zzz[i] * fz_0;
    }

    // Set up 50-60 components of targeted buffer : FF

    auto tk_xzz_xxx = pbuffer.data(idx_kin_ff + 50);

    auto tk_xzz_xxy = pbuffer.data(idx_kin_ff + 51);

    auto tk_xzz_xxz = pbuffer.data(idx_kin_ff + 52);

    auto tk_xzz_xyy = pbuffer.data(idx_kin_ff + 53);

    auto tk_xzz_xyz = pbuffer.data(idx_kin_ff + 54);

    auto tk_xzz_xzz = pbuffer.data(idx_kin_ff + 55);

    auto tk_xzz_yyy = pbuffer.data(idx_kin_ff + 56);

    auto tk_xzz_yyz = pbuffer.data(idx_kin_ff + 57);

    auto tk_xzz_yzz = pbuffer.data(idx_kin_ff + 58);

    auto tk_xzz_zzz = pbuffer.data(idx_kin_ff + 59);

#pragma omp simd aligned(pa_x,           \
                             tk_xzz_xxx, \
                             tk_xzz_xxy, \
                             tk_xzz_xxz, \
                             tk_xzz_xyy, \
                             tk_xzz_xyz, \
                             tk_xzz_xzz, \
                             tk_xzz_yyy, \
                             tk_xzz_yyz, \
                             tk_xzz_yzz, \
                             tk_xzz_zzz, \
                             tk_zz_xx,   \
                             tk_zz_xxx,  \
                             tk_zz_xxy,  \
                             tk_zz_xxz,  \
                             tk_zz_xy,   \
                             tk_zz_xyy,  \
                             tk_zz_xyz,  \
                             tk_zz_xz,   \
                             tk_zz_xzz,  \
                             tk_zz_yy,   \
                             tk_zz_yyy,  \
                             tk_zz_yyz,  \
                             tk_zz_yz,   \
                             tk_zz_yzz,  \
                             tk_zz_zz,   \
                             tk_zz_zzz,  \
                             ts_xzz_xxx, \
                             ts_xzz_xxy, \
                             ts_xzz_xxz, \
                             ts_xzz_xyy, \
                             ts_xzz_xyz, \
                             ts_xzz_xzz, \
                             ts_xzz_yyy, \
                             ts_xzz_yyz, \
                             ts_xzz_yzz, \
                             ts_xzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzz_xxx[i] = 3.0 * tk_zz_xx[i] * fe_0 + tk_zz_xxx[i] * pa_x[i] + 2.0 * ts_xzz_xxx[i] * fz_0;

        tk_xzz_xxy[i] = 2.0 * tk_zz_xy[i] * fe_0 + tk_zz_xxy[i] * pa_x[i] + 2.0 * ts_xzz_xxy[i] * fz_0;

        tk_xzz_xxz[i] = 2.0 * tk_zz_xz[i] * fe_0 + tk_zz_xxz[i] * pa_x[i] + 2.0 * ts_xzz_xxz[i] * fz_0;

        tk_xzz_xyy[i] = tk_zz_yy[i] * fe_0 + tk_zz_xyy[i] * pa_x[i] + 2.0 * ts_xzz_xyy[i] * fz_0;

        tk_xzz_xyz[i] = tk_zz_yz[i] * fe_0 + tk_zz_xyz[i] * pa_x[i] + 2.0 * ts_xzz_xyz[i] * fz_0;

        tk_xzz_xzz[i] = tk_zz_zz[i] * fe_0 + tk_zz_xzz[i] * pa_x[i] + 2.0 * ts_xzz_xzz[i] * fz_0;

        tk_xzz_yyy[i] = tk_zz_yyy[i] * pa_x[i] + 2.0 * ts_xzz_yyy[i] * fz_0;

        tk_xzz_yyz[i] = tk_zz_yyz[i] * pa_x[i] + 2.0 * ts_xzz_yyz[i] * fz_0;

        tk_xzz_yzz[i] = tk_zz_yzz[i] * pa_x[i] + 2.0 * ts_xzz_yzz[i] * fz_0;

        tk_xzz_zzz[i] = tk_zz_zzz[i] * pa_x[i] + 2.0 * ts_xzz_zzz[i] * fz_0;
    }

    // Set up 60-70 components of targeted buffer : FF

    auto tk_yyy_xxx = pbuffer.data(idx_kin_ff + 60);

    auto tk_yyy_xxy = pbuffer.data(idx_kin_ff + 61);

    auto tk_yyy_xxz = pbuffer.data(idx_kin_ff + 62);

    auto tk_yyy_xyy = pbuffer.data(idx_kin_ff + 63);

    auto tk_yyy_xyz = pbuffer.data(idx_kin_ff + 64);

    auto tk_yyy_xzz = pbuffer.data(idx_kin_ff + 65);

    auto tk_yyy_yyy = pbuffer.data(idx_kin_ff + 66);

    auto tk_yyy_yyz = pbuffer.data(idx_kin_ff + 67);

    auto tk_yyy_yzz = pbuffer.data(idx_kin_ff + 68);

    auto tk_yyy_zzz = pbuffer.data(idx_kin_ff + 69);

#pragma omp simd aligned(pa_y,           \
                             tk_y_xxx,   \
                             tk_y_xxy,   \
                             tk_y_xxz,   \
                             tk_y_xyy,   \
                             tk_y_xyz,   \
                             tk_y_xzz,   \
                             tk_y_yyy,   \
                             tk_y_yyz,   \
                             tk_y_yzz,   \
                             tk_y_zzz,   \
                             tk_yy_xx,   \
                             tk_yy_xxx,  \
                             tk_yy_xxy,  \
                             tk_yy_xxz,  \
                             tk_yy_xy,   \
                             tk_yy_xyy,  \
                             tk_yy_xyz,  \
                             tk_yy_xz,   \
                             tk_yy_xzz,  \
                             tk_yy_yy,   \
                             tk_yy_yyy,  \
                             tk_yy_yyz,  \
                             tk_yy_yz,   \
                             tk_yy_yzz,  \
                             tk_yy_zz,   \
                             tk_yy_zzz,  \
                             tk_yyy_xxx, \
                             tk_yyy_xxy, \
                             tk_yyy_xxz, \
                             tk_yyy_xyy, \
                             tk_yyy_xyz, \
                             tk_yyy_xzz, \
                             tk_yyy_yyy, \
                             tk_yyy_yyz, \
                             tk_yyy_yzz, \
                             tk_yyy_zzz, \
                             ts_y_xxx,   \
                             ts_y_xxy,   \
                             ts_y_xxz,   \
                             ts_y_xyy,   \
                             ts_y_xyz,   \
                             ts_y_xzz,   \
                             ts_y_yyy,   \
                             ts_y_yyz,   \
                             ts_y_yzz,   \
                             ts_y_zzz,   \
                             ts_yyy_xxx, \
                             ts_yyy_xxy, \
                             ts_yyy_xxz, \
                             ts_yyy_xyy, \
                             ts_yyy_xyz, \
                             ts_yyy_xzz, \
                             ts_yyy_yyy, \
                             ts_yyy_yyz, \
                             ts_yyy_yzz, \
                             ts_yyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyy_xxx[i] = -4.0 * ts_y_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxx[i] * fe_0 + tk_yy_xxx[i] * pa_y[i] +
                        2.0 * ts_yyy_xxx[i] * fz_0;

        tk_yyy_xxy[i] = -4.0 * ts_y_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxy[i] * fe_0 + tk_yy_xx[i] * fe_0 +
                        tk_yy_xxy[i] * pa_y[i] + 2.0 * ts_yyy_xxy[i] * fz_0;

        tk_yyy_xxz[i] = -4.0 * ts_y_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxz[i] * fe_0 + tk_yy_xxz[i] * pa_y[i] +
                        2.0 * ts_yyy_xxz[i] * fz_0;

        tk_yyy_xyy[i] = -4.0 * ts_y_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyy[i] * fe_0 + 2.0 * tk_yy_xy[i] * fe_0 +
                        tk_yy_xyy[i] * pa_y[i] + 2.0 * ts_yyy_xyy[i] * fz_0;

        tk_yyy_xyz[i] = -4.0 * ts_y_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyz[i] * fe_0 + tk_yy_xz[i] * fe_0 +
                        tk_yy_xyz[i] * pa_y[i] + 2.0 * ts_yyy_xyz[i] * fz_0;

        tk_yyy_xzz[i] = -4.0 * ts_y_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xzz[i] * fe_0 + tk_yy_xzz[i] * pa_y[i] +
                        2.0 * ts_yyy_xzz[i] * fz_0;

        tk_yyy_yyy[i] = -4.0 * ts_y_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyy[i] * fe_0 + 3.0 * tk_yy_yy[i] * fe_0 +
                        tk_yy_yyy[i] * pa_y[i] + 2.0 * ts_yyy_yyy[i] * fz_0;

        tk_yyy_yyz[i] = -4.0 * ts_y_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyz[i] * fe_0 + 2.0 * tk_yy_yz[i] * fe_0 +
                        tk_yy_yyz[i] * pa_y[i] + 2.0 * ts_yyy_yyz[i] * fz_0;

        tk_yyy_yzz[i] = -4.0 * ts_y_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yzz[i] * fe_0 + tk_yy_zz[i] * fe_0 +
                        tk_yy_yzz[i] * pa_y[i] + 2.0 * ts_yyy_yzz[i] * fz_0;

        tk_yyy_zzz[i] = -4.0 * ts_y_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_zzz[i] * fe_0 + tk_yy_zzz[i] * pa_y[i] +
                        2.0 * ts_yyy_zzz[i] * fz_0;
    }

    // Set up 70-80 components of targeted buffer : FF

    auto tk_yyz_xxx = pbuffer.data(idx_kin_ff + 70);

    auto tk_yyz_xxy = pbuffer.data(idx_kin_ff + 71);

    auto tk_yyz_xxz = pbuffer.data(idx_kin_ff + 72);

    auto tk_yyz_xyy = pbuffer.data(idx_kin_ff + 73);

    auto tk_yyz_xyz = pbuffer.data(idx_kin_ff + 74);

    auto tk_yyz_xzz = pbuffer.data(idx_kin_ff + 75);

    auto tk_yyz_yyy = pbuffer.data(idx_kin_ff + 76);

    auto tk_yyz_yyz = pbuffer.data(idx_kin_ff + 77);

    auto tk_yyz_yzz = pbuffer.data(idx_kin_ff + 78);

    auto tk_yyz_zzz = pbuffer.data(idx_kin_ff + 79);

#pragma omp simd aligned(pa_z,           \
                             tk_yy_xx,   \
                             tk_yy_xxx,  \
                             tk_yy_xxy,  \
                             tk_yy_xxz,  \
                             tk_yy_xy,   \
                             tk_yy_xyy,  \
                             tk_yy_xyz,  \
                             tk_yy_xz,   \
                             tk_yy_xzz,  \
                             tk_yy_yy,   \
                             tk_yy_yyy,  \
                             tk_yy_yyz,  \
                             tk_yy_yz,   \
                             tk_yy_yzz,  \
                             tk_yy_zz,   \
                             tk_yy_zzz,  \
                             tk_yyz_xxx, \
                             tk_yyz_xxy, \
                             tk_yyz_xxz, \
                             tk_yyz_xyy, \
                             tk_yyz_xyz, \
                             tk_yyz_xzz, \
                             tk_yyz_yyy, \
                             tk_yyz_yyz, \
                             tk_yyz_yzz, \
                             tk_yyz_zzz, \
                             ts_yyz_xxx, \
                             ts_yyz_xxy, \
                             ts_yyz_xxz, \
                             ts_yyz_xyy, \
                             ts_yyz_xyz, \
                             ts_yyz_xzz, \
                             ts_yyz_yyy, \
                             ts_yyz_yyz, \
                             ts_yyz_yzz, \
                             ts_yyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyz_xxx[i] = tk_yy_xxx[i] * pa_z[i] + 2.0 * ts_yyz_xxx[i] * fz_0;

        tk_yyz_xxy[i] = tk_yy_xxy[i] * pa_z[i] + 2.0 * ts_yyz_xxy[i] * fz_0;

        tk_yyz_xxz[i] = tk_yy_xx[i] * fe_0 + tk_yy_xxz[i] * pa_z[i] + 2.0 * ts_yyz_xxz[i] * fz_0;

        tk_yyz_xyy[i] = tk_yy_xyy[i] * pa_z[i] + 2.0 * ts_yyz_xyy[i] * fz_0;

        tk_yyz_xyz[i] = tk_yy_xy[i] * fe_0 + tk_yy_xyz[i] * pa_z[i] + 2.0 * ts_yyz_xyz[i] * fz_0;

        tk_yyz_xzz[i] = 2.0 * tk_yy_xz[i] * fe_0 + tk_yy_xzz[i] * pa_z[i] + 2.0 * ts_yyz_xzz[i] * fz_0;

        tk_yyz_yyy[i] = tk_yy_yyy[i] * pa_z[i] + 2.0 * ts_yyz_yyy[i] * fz_0;

        tk_yyz_yyz[i] = tk_yy_yy[i] * fe_0 + tk_yy_yyz[i] * pa_z[i] + 2.0 * ts_yyz_yyz[i] * fz_0;

        tk_yyz_yzz[i] = 2.0 * tk_yy_yz[i] * fe_0 + tk_yy_yzz[i] * pa_z[i] + 2.0 * ts_yyz_yzz[i] * fz_0;

        tk_yyz_zzz[i] = 3.0 * tk_yy_zz[i] * fe_0 + tk_yy_zzz[i] * pa_z[i] + 2.0 * ts_yyz_zzz[i] * fz_0;
    }

    // Set up 80-90 components of targeted buffer : FF

    auto tk_yzz_xxx = pbuffer.data(idx_kin_ff + 80);

    auto tk_yzz_xxy = pbuffer.data(idx_kin_ff + 81);

    auto tk_yzz_xxz = pbuffer.data(idx_kin_ff + 82);

    auto tk_yzz_xyy = pbuffer.data(idx_kin_ff + 83);

    auto tk_yzz_xyz = pbuffer.data(idx_kin_ff + 84);

    auto tk_yzz_xzz = pbuffer.data(idx_kin_ff + 85);

    auto tk_yzz_yyy = pbuffer.data(idx_kin_ff + 86);

    auto tk_yzz_yyz = pbuffer.data(idx_kin_ff + 87);

    auto tk_yzz_yzz = pbuffer.data(idx_kin_ff + 88);

    auto tk_yzz_zzz = pbuffer.data(idx_kin_ff + 89);

#pragma omp simd aligned(pa_y,           \
                             tk_yzz_xxx, \
                             tk_yzz_xxy, \
                             tk_yzz_xxz, \
                             tk_yzz_xyy, \
                             tk_yzz_xyz, \
                             tk_yzz_xzz, \
                             tk_yzz_yyy, \
                             tk_yzz_yyz, \
                             tk_yzz_yzz, \
                             tk_yzz_zzz, \
                             tk_zz_xx,   \
                             tk_zz_xxx,  \
                             tk_zz_xxy,  \
                             tk_zz_xxz,  \
                             tk_zz_xy,   \
                             tk_zz_xyy,  \
                             tk_zz_xyz,  \
                             tk_zz_xz,   \
                             tk_zz_xzz,  \
                             tk_zz_yy,   \
                             tk_zz_yyy,  \
                             tk_zz_yyz,  \
                             tk_zz_yz,   \
                             tk_zz_yzz,  \
                             tk_zz_zz,   \
                             tk_zz_zzz,  \
                             ts_yzz_xxx, \
                             ts_yzz_xxy, \
                             ts_yzz_xxz, \
                             ts_yzz_xyy, \
                             ts_yzz_xyz, \
                             ts_yzz_xzz, \
                             ts_yzz_yyy, \
                             ts_yzz_yyz, \
                             ts_yzz_yzz, \
                             ts_yzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzz_xxx[i] = tk_zz_xxx[i] * pa_y[i] + 2.0 * ts_yzz_xxx[i] * fz_0;

        tk_yzz_xxy[i] = tk_zz_xx[i] * fe_0 + tk_zz_xxy[i] * pa_y[i] + 2.0 * ts_yzz_xxy[i] * fz_0;

        tk_yzz_xxz[i] = tk_zz_xxz[i] * pa_y[i] + 2.0 * ts_yzz_xxz[i] * fz_0;

        tk_yzz_xyy[i] = 2.0 * tk_zz_xy[i] * fe_0 + tk_zz_xyy[i] * pa_y[i] + 2.0 * ts_yzz_xyy[i] * fz_0;

        tk_yzz_xyz[i] = tk_zz_xz[i] * fe_0 + tk_zz_xyz[i] * pa_y[i] + 2.0 * ts_yzz_xyz[i] * fz_0;

        tk_yzz_xzz[i] = tk_zz_xzz[i] * pa_y[i] + 2.0 * ts_yzz_xzz[i] * fz_0;

        tk_yzz_yyy[i] = 3.0 * tk_zz_yy[i] * fe_0 + tk_zz_yyy[i] * pa_y[i] + 2.0 * ts_yzz_yyy[i] * fz_0;

        tk_yzz_yyz[i] = 2.0 * tk_zz_yz[i] * fe_0 + tk_zz_yyz[i] * pa_y[i] + 2.0 * ts_yzz_yyz[i] * fz_0;

        tk_yzz_yzz[i] = tk_zz_zz[i] * fe_0 + tk_zz_yzz[i] * pa_y[i] + 2.0 * ts_yzz_yzz[i] * fz_0;

        tk_yzz_zzz[i] = tk_zz_zzz[i] * pa_y[i] + 2.0 * ts_yzz_zzz[i] * fz_0;
    }

    // Set up 90-100 components of targeted buffer : FF

    auto tk_zzz_xxx = pbuffer.data(idx_kin_ff + 90);

    auto tk_zzz_xxy = pbuffer.data(idx_kin_ff + 91);

    auto tk_zzz_xxz = pbuffer.data(idx_kin_ff + 92);

    auto tk_zzz_xyy = pbuffer.data(idx_kin_ff + 93);

    auto tk_zzz_xyz = pbuffer.data(idx_kin_ff + 94);

    auto tk_zzz_xzz = pbuffer.data(idx_kin_ff + 95);

    auto tk_zzz_yyy = pbuffer.data(idx_kin_ff + 96);

    auto tk_zzz_yyz = pbuffer.data(idx_kin_ff + 97);

    auto tk_zzz_yzz = pbuffer.data(idx_kin_ff + 98);

    auto tk_zzz_zzz = pbuffer.data(idx_kin_ff + 99);

#pragma omp simd aligned(pa_z,           \
                             tk_z_xxx,   \
                             tk_z_xxy,   \
                             tk_z_xxz,   \
                             tk_z_xyy,   \
                             tk_z_xyz,   \
                             tk_z_xzz,   \
                             tk_z_yyy,   \
                             tk_z_yyz,   \
                             tk_z_yzz,   \
                             tk_z_zzz,   \
                             tk_zz_xx,   \
                             tk_zz_xxx,  \
                             tk_zz_xxy,  \
                             tk_zz_xxz,  \
                             tk_zz_xy,   \
                             tk_zz_xyy,  \
                             tk_zz_xyz,  \
                             tk_zz_xz,   \
                             tk_zz_xzz,  \
                             tk_zz_yy,   \
                             tk_zz_yyy,  \
                             tk_zz_yyz,  \
                             tk_zz_yz,   \
                             tk_zz_yzz,  \
                             tk_zz_zz,   \
                             tk_zz_zzz,  \
                             tk_zzz_xxx, \
                             tk_zzz_xxy, \
                             tk_zzz_xxz, \
                             tk_zzz_xyy, \
                             tk_zzz_xyz, \
                             tk_zzz_xzz, \
                             tk_zzz_yyy, \
                             tk_zzz_yyz, \
                             tk_zzz_yzz, \
                             tk_zzz_zzz, \
                             ts_z_xxx,   \
                             ts_z_xxy,   \
                             ts_z_xxz,   \
                             ts_z_xyy,   \
                             ts_z_xyz,   \
                             ts_z_xzz,   \
                             ts_z_yyy,   \
                             ts_z_yyz,   \
                             ts_z_yzz,   \
                             ts_z_zzz,   \
                             ts_zzz_xxx, \
                             ts_zzz_xxy, \
                             ts_zzz_xxz, \
                             ts_zzz_xyy, \
                             ts_zzz_xyz, \
                             ts_zzz_xzz, \
                             ts_zzz_yyy, \
                             ts_zzz_yyz, \
                             ts_zzz_yzz, \
                             ts_zzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzz_xxx[i] = -4.0 * ts_z_xxx[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxx[i] * fe_0 + tk_zz_xxx[i] * pa_z[i] +
                        2.0 * ts_zzz_xxx[i] * fz_0;

        tk_zzz_xxy[i] = -4.0 * ts_z_xxy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxy[i] * fe_0 + tk_zz_xxy[i] * pa_z[i] +
                        2.0 * ts_zzz_xxy[i] * fz_0;

        tk_zzz_xxz[i] = -4.0 * ts_z_xxz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxz[i] * fe_0 + tk_zz_xx[i] * fe_0 +
                        tk_zz_xxz[i] * pa_z[i] + 2.0 * ts_zzz_xxz[i] * fz_0;

        tk_zzz_xyy[i] = -4.0 * ts_z_xyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyy[i] * fe_0 + tk_zz_xyy[i] * pa_z[i] +
                        2.0 * ts_zzz_xyy[i] * fz_0;

        tk_zzz_xyz[i] = -4.0 * ts_z_xyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyz[i] * fe_0 + tk_zz_xy[i] * fe_0 +
                        tk_zz_xyz[i] * pa_z[i] + 2.0 * ts_zzz_xyz[i] * fz_0;

        tk_zzz_xzz[i] = -4.0 * ts_z_xzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xzz[i] * fe_0 + 2.0 * tk_zz_xz[i] * fe_0 +
                        tk_zz_xzz[i] * pa_z[i] + 2.0 * ts_zzz_xzz[i] * fz_0;

        tk_zzz_yyy[i] = -4.0 * ts_z_yyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyy[i] * fe_0 + tk_zz_yyy[i] * pa_z[i] +
                        2.0 * ts_zzz_yyy[i] * fz_0;

        tk_zzz_yyz[i] = -4.0 * ts_z_yyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyz[i] * fe_0 + tk_zz_yy[i] * fe_0 +
                        tk_zz_yyz[i] * pa_z[i] + 2.0 * ts_zzz_yyz[i] * fz_0;

        tk_zzz_yzz[i] = -4.0 * ts_z_yzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yzz[i] * fe_0 + 2.0 * tk_zz_yz[i] * fe_0 +
                        tk_zz_yzz[i] * pa_z[i] + 2.0 * ts_zzz_yzz[i] * fz_0;

        tk_zzz_zzz[i] = -4.0 * ts_z_zzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_zzz[i] * fe_0 + 3.0 * tk_zz_zz[i] * fe_0 +
                        tk_zz_zzz[i] * pa_z[i] + 2.0 * ts_zzz_zzz[i] * fz_0;
    }
}

}  // namespace kinrec
