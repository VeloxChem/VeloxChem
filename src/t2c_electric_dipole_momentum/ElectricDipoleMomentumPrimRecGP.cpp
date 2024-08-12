#include "ElectricDipoleMomentumPrimRecGP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_gp(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_gp,
                                      const size_t              idx_dip_dp,
                                      const size_t              idx_dip_fs,
                                      const size_t              idx_ovl_fp,
                                      const size_t              idx_dip_fp,
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

    // Set up components of auxiliary buffer : DP

    auto tr_x_xx_x = pbuffer.data(idx_dip_dp);

    auto tr_x_xx_y = pbuffer.data(idx_dip_dp + 1);

    auto tr_x_xx_z = pbuffer.data(idx_dip_dp + 2);

    auto tr_x_xy_x = pbuffer.data(idx_dip_dp + 3);

    auto tr_x_xz_x = pbuffer.data(idx_dip_dp + 6);

    auto tr_x_yy_x = pbuffer.data(idx_dip_dp + 9);

    auto tr_x_yy_y = pbuffer.data(idx_dip_dp + 10);

    auto tr_x_yy_z = pbuffer.data(idx_dip_dp + 11);

    auto tr_x_yz_z = pbuffer.data(idx_dip_dp + 14);

    auto tr_x_zz_x = pbuffer.data(idx_dip_dp + 15);

    auto tr_x_zz_y = pbuffer.data(idx_dip_dp + 16);

    auto tr_x_zz_z = pbuffer.data(idx_dip_dp + 17);

    auto tr_y_xx_x = pbuffer.data(idx_dip_dp + 18);

    auto tr_y_xx_y = pbuffer.data(idx_dip_dp + 19);

    auto tr_y_xx_z = pbuffer.data(idx_dip_dp + 20);

    auto tr_y_xy_y = pbuffer.data(idx_dip_dp + 22);

    auto tr_y_xy_z = pbuffer.data(idx_dip_dp + 23);

    auto tr_y_xz_z = pbuffer.data(idx_dip_dp + 26);

    auto tr_y_yy_x = pbuffer.data(idx_dip_dp + 27);

    auto tr_y_yy_y = pbuffer.data(idx_dip_dp + 28);

    auto tr_y_yy_z = pbuffer.data(idx_dip_dp + 29);

    auto tr_y_yz_y = pbuffer.data(idx_dip_dp + 31);

    auto tr_y_yz_z = pbuffer.data(idx_dip_dp + 32);

    auto tr_y_zz_x = pbuffer.data(idx_dip_dp + 33);

    auto tr_y_zz_y = pbuffer.data(idx_dip_dp + 34);

    auto tr_y_zz_z = pbuffer.data(idx_dip_dp + 35);

    auto tr_z_xx_x = pbuffer.data(idx_dip_dp + 36);

    auto tr_z_xx_y = pbuffer.data(idx_dip_dp + 37);

    auto tr_z_xx_z = pbuffer.data(idx_dip_dp + 38);

    auto tr_z_xy_y = pbuffer.data(idx_dip_dp + 40);

    auto tr_z_xz_y = pbuffer.data(idx_dip_dp + 43);

    auto tr_z_xz_z = pbuffer.data(idx_dip_dp + 44);

    auto tr_z_yy_x = pbuffer.data(idx_dip_dp + 45);

    auto tr_z_yy_y = pbuffer.data(idx_dip_dp + 46);

    auto tr_z_yy_z = pbuffer.data(idx_dip_dp + 47);

    auto tr_z_yz_x = pbuffer.data(idx_dip_dp + 48);

    auto tr_z_yz_y = pbuffer.data(idx_dip_dp + 49);

    auto tr_z_yz_z = pbuffer.data(idx_dip_dp + 50);

    auto tr_z_zz_x = pbuffer.data(idx_dip_dp + 51);

    auto tr_z_zz_y = pbuffer.data(idx_dip_dp + 52);

    auto tr_z_zz_z = pbuffer.data(idx_dip_dp + 53);

    // Set up components of auxiliary buffer : FS

    auto tr_x_xxx_0 = pbuffer.data(idx_dip_fs);

    auto tr_x_yyy_0 = pbuffer.data(idx_dip_fs + 6);

    auto tr_x_zzz_0 = pbuffer.data(idx_dip_fs + 9);

    auto tr_y_xxx_0 = pbuffer.data(idx_dip_fs + 10);

    auto tr_y_xyy_0 = pbuffer.data(idx_dip_fs + 13);

    auto tr_y_yyy_0 = pbuffer.data(idx_dip_fs + 16);

    auto tr_y_yzz_0 = pbuffer.data(idx_dip_fs + 18);

    auto tr_y_zzz_0 = pbuffer.data(idx_dip_fs + 19);

    auto tr_z_xxx_0 = pbuffer.data(idx_dip_fs + 20);

    auto tr_z_xzz_0 = pbuffer.data(idx_dip_fs + 25);

    auto tr_z_yyy_0 = pbuffer.data(idx_dip_fs + 26);

    auto tr_z_yyz_0 = pbuffer.data(idx_dip_fs + 27);

    auto tr_z_yzz_0 = pbuffer.data(idx_dip_fs + 28);

    auto tr_z_zzz_0 = pbuffer.data(idx_dip_fs + 29);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_ovl_fp);

    auto ts_xxx_y = pbuffer.data(idx_ovl_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_ovl_fp + 2);

    auto ts_xyy_y = pbuffer.data(idx_ovl_fp + 10);

    auto ts_xzz_z = pbuffer.data(idx_ovl_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_ovl_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_ovl_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_ovl_fp + 20);

    auto ts_yyz_z = pbuffer.data(idx_ovl_fp + 23);

    auto ts_yzz_y = pbuffer.data(idx_ovl_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_ovl_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_ovl_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_ovl_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_ovl_fp + 29);

    // Set up components of auxiliary buffer : FP

    auto tr_x_xxx_x = pbuffer.data(idx_dip_fp);

    auto tr_x_xxx_y = pbuffer.data(idx_dip_fp + 1);

    auto tr_x_xxx_z = pbuffer.data(idx_dip_fp + 2);

    auto tr_x_xxy_x = pbuffer.data(idx_dip_fp + 3);

    auto tr_x_xxy_y = pbuffer.data(idx_dip_fp + 4);

    auto tr_x_xxy_z = pbuffer.data(idx_dip_fp + 5);

    auto tr_x_xxz_x = pbuffer.data(idx_dip_fp + 6);

    auto tr_x_xxz_y = pbuffer.data(idx_dip_fp + 7);

    auto tr_x_xxz_z = pbuffer.data(idx_dip_fp + 8);

    auto tr_x_xyy_x = pbuffer.data(idx_dip_fp + 9);

    auto tr_x_xyy_y = pbuffer.data(idx_dip_fp + 10);

    auto tr_x_xzz_x = pbuffer.data(idx_dip_fp + 15);

    auto tr_x_xzz_z = pbuffer.data(idx_dip_fp + 17);

    auto tr_x_yyy_x = pbuffer.data(idx_dip_fp + 18);

    auto tr_x_yyy_y = pbuffer.data(idx_dip_fp + 19);

    auto tr_x_yyy_z = pbuffer.data(idx_dip_fp + 20);

    auto tr_x_yyz_y = pbuffer.data(idx_dip_fp + 22);

    auto tr_x_yyz_z = pbuffer.data(idx_dip_fp + 23);

    auto tr_x_yzz_x = pbuffer.data(idx_dip_fp + 24);

    auto tr_x_yzz_y = pbuffer.data(idx_dip_fp + 25);

    auto tr_x_yzz_z = pbuffer.data(idx_dip_fp + 26);

    auto tr_x_zzz_x = pbuffer.data(idx_dip_fp + 27);

    auto tr_x_zzz_y = pbuffer.data(idx_dip_fp + 28);

    auto tr_x_zzz_z = pbuffer.data(idx_dip_fp + 29);

    auto tr_y_xxx_x = pbuffer.data(idx_dip_fp + 30);

    auto tr_y_xxx_y = pbuffer.data(idx_dip_fp + 31);

    auto tr_y_xxx_z = pbuffer.data(idx_dip_fp + 32);

    auto tr_y_xxy_x = pbuffer.data(idx_dip_fp + 33);

    auto tr_y_xxy_y = pbuffer.data(idx_dip_fp + 34);

    auto tr_y_xxy_z = pbuffer.data(idx_dip_fp + 35);

    auto tr_y_xxz_x = pbuffer.data(idx_dip_fp + 36);

    auto tr_y_xxz_z = pbuffer.data(idx_dip_fp + 38);

    auto tr_y_xyy_x = pbuffer.data(idx_dip_fp + 39);

    auto tr_y_xyy_y = pbuffer.data(idx_dip_fp + 40);

    auto tr_y_xyy_z = pbuffer.data(idx_dip_fp + 41);

    auto tr_y_xyz_z = pbuffer.data(idx_dip_fp + 44);

    auto tr_y_xzz_y = pbuffer.data(idx_dip_fp + 46);

    auto tr_y_xzz_z = pbuffer.data(idx_dip_fp + 47);

    auto tr_y_yyy_x = pbuffer.data(idx_dip_fp + 48);

    auto tr_y_yyy_y = pbuffer.data(idx_dip_fp + 49);

    auto tr_y_yyy_z = pbuffer.data(idx_dip_fp + 50);

    auto tr_y_yyz_x = pbuffer.data(idx_dip_fp + 51);

    auto tr_y_yyz_y = pbuffer.data(idx_dip_fp + 52);

    auto tr_y_yyz_z = pbuffer.data(idx_dip_fp + 53);

    auto tr_y_yzz_x = pbuffer.data(idx_dip_fp + 54);

    auto tr_y_yzz_y = pbuffer.data(idx_dip_fp + 55);

    auto tr_y_yzz_z = pbuffer.data(idx_dip_fp + 56);

    auto tr_y_zzz_x = pbuffer.data(idx_dip_fp + 57);

    auto tr_y_zzz_y = pbuffer.data(idx_dip_fp + 58);

    auto tr_y_zzz_z = pbuffer.data(idx_dip_fp + 59);

    auto tr_z_xxx_x = pbuffer.data(idx_dip_fp + 60);

    auto tr_z_xxx_y = pbuffer.data(idx_dip_fp + 61);

    auto tr_z_xxx_z = pbuffer.data(idx_dip_fp + 62);

    auto tr_z_xxy_x = pbuffer.data(idx_dip_fp + 63);

    auto tr_z_xxy_y = pbuffer.data(idx_dip_fp + 64);

    auto tr_z_xxz_x = pbuffer.data(idx_dip_fp + 66);

    auto tr_z_xxz_y = pbuffer.data(idx_dip_fp + 67);

    auto tr_z_xxz_z = pbuffer.data(idx_dip_fp + 68);

    auto tr_z_xyy_y = pbuffer.data(idx_dip_fp + 70);

    auto tr_z_xyy_z = pbuffer.data(idx_dip_fp + 71);

    auto tr_z_xyz_y = pbuffer.data(idx_dip_fp + 73);

    auto tr_z_xzz_x = pbuffer.data(idx_dip_fp + 75);

    auto tr_z_xzz_y = pbuffer.data(idx_dip_fp + 76);

    auto tr_z_xzz_z = pbuffer.data(idx_dip_fp + 77);

    auto tr_z_yyy_x = pbuffer.data(idx_dip_fp + 78);

    auto tr_z_yyy_y = pbuffer.data(idx_dip_fp + 79);

    auto tr_z_yyy_z = pbuffer.data(idx_dip_fp + 80);

    auto tr_z_yyz_x = pbuffer.data(idx_dip_fp + 81);

    auto tr_z_yyz_y = pbuffer.data(idx_dip_fp + 82);

    auto tr_z_yyz_z = pbuffer.data(idx_dip_fp + 83);

    auto tr_z_yzz_x = pbuffer.data(idx_dip_fp + 84);

    auto tr_z_yzz_y = pbuffer.data(idx_dip_fp + 85);

    auto tr_z_yzz_z = pbuffer.data(idx_dip_fp + 86);

    auto tr_z_zzz_x = pbuffer.data(idx_dip_fp + 87);

    auto tr_z_zzz_y = pbuffer.data(idx_dip_fp + 88);

    auto tr_z_zzz_z = pbuffer.data(idx_dip_fp + 89);

    // Set up 0-3 components of targeted buffer : GP

    auto tr_x_xxxx_x = pbuffer.data(idx_dip_gp);

    auto tr_x_xxxx_y = pbuffer.data(idx_dip_gp + 1);

    auto tr_x_xxxx_z = pbuffer.data(idx_dip_gp + 2);

#pragma omp simd aligned(pa_x,            \
                             tr_x_xx_x,   \
                             tr_x_xx_y,   \
                             tr_x_xx_z,   \
                             tr_x_xxx_0,  \
                             tr_x_xxx_x,  \
                             tr_x_xxx_y,  \
                             tr_x_xxx_z,  \
                             tr_x_xxxx_x, \
                             tr_x_xxxx_y, \
                             tr_x_xxxx_z, \
                             ts_xxx_x,    \
                             ts_xxx_y,    \
                             ts_xxx_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_x[i] = 3.0 * tr_x_xx_x[i] * fe_0 + tr_x_xxx_0[i] * fe_0 + ts_xxx_x[i] * fe_0 + tr_x_xxx_x[i] * pa_x[i];

        tr_x_xxxx_y[i] = 3.0 * tr_x_xx_y[i] * fe_0 + ts_xxx_y[i] * fe_0 + tr_x_xxx_y[i] * pa_x[i];

        tr_x_xxxx_z[i] = 3.0 * tr_x_xx_z[i] * fe_0 + ts_xxx_z[i] * fe_0 + tr_x_xxx_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto tr_x_xxxy_x = pbuffer.data(idx_dip_gp + 3);

    auto tr_x_xxxy_y = pbuffer.data(idx_dip_gp + 4);

    auto tr_x_xxxy_z = pbuffer.data(idx_dip_gp + 5);

#pragma omp simd aligned(pa_y, tr_x_xxx_0, tr_x_xxx_x, tr_x_xxx_y, tr_x_xxx_z, tr_x_xxxy_x, tr_x_xxxy_y, tr_x_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxy_x[i] = tr_x_xxx_x[i] * pa_y[i];

        tr_x_xxxy_y[i] = tr_x_xxx_0[i] * fe_0 + tr_x_xxx_y[i] * pa_y[i];

        tr_x_xxxy_z[i] = tr_x_xxx_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto tr_x_xxxz_x = pbuffer.data(idx_dip_gp + 6);

    auto tr_x_xxxz_y = pbuffer.data(idx_dip_gp + 7);

    auto tr_x_xxxz_z = pbuffer.data(idx_dip_gp + 8);

#pragma omp simd aligned(pa_z, tr_x_xxx_0, tr_x_xxx_x, tr_x_xxx_y, tr_x_xxx_z, tr_x_xxxz_x, tr_x_xxxz_y, tr_x_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxz_x[i] = tr_x_xxx_x[i] * pa_z[i];

        tr_x_xxxz_y[i] = tr_x_xxx_y[i] * pa_z[i];

        tr_x_xxxz_z[i] = tr_x_xxx_0[i] * fe_0 + tr_x_xxx_z[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto tr_x_xxyy_x = pbuffer.data(idx_dip_gp + 9);

    auto tr_x_xxyy_y = pbuffer.data(idx_dip_gp + 10);

    auto tr_x_xxyy_z = pbuffer.data(idx_dip_gp + 11);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_x_xx_x,   \
                             tr_x_xx_z,   \
                             tr_x_xxy_x,  \
                             tr_x_xxy_z,  \
                             tr_x_xxyy_x, \
                             tr_x_xxyy_y, \
                             tr_x_xxyy_z, \
                             tr_x_xyy_y,  \
                             tr_x_yy_y,   \
                             ts_xyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyy_x[i] = tr_x_xx_x[i] * fe_0 + tr_x_xxy_x[i] * pa_y[i];

        tr_x_xxyy_y[i] = tr_x_yy_y[i] * fe_0 + ts_xyy_y[i] * fe_0 + tr_x_xyy_y[i] * pa_x[i];

        tr_x_xxyy_z[i] = tr_x_xx_z[i] * fe_0 + tr_x_xxy_z[i] * pa_y[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto tr_x_xxyz_x = pbuffer.data(idx_dip_gp + 12);

    auto tr_x_xxyz_y = pbuffer.data(idx_dip_gp + 13);

    auto tr_x_xxyz_z = pbuffer.data(idx_dip_gp + 14);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xxy_y, tr_x_xxyz_x, tr_x_xxyz_y, tr_x_xxyz_z, tr_x_xxz_x, tr_x_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_x_xxyz_x[i] = tr_x_xxz_x[i] * pa_y[i];

        tr_x_xxyz_y[i] = tr_x_xxy_y[i] * pa_z[i];

        tr_x_xxyz_z[i] = tr_x_xxz_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto tr_x_xxzz_x = pbuffer.data(idx_dip_gp + 15);

    auto tr_x_xxzz_y = pbuffer.data(idx_dip_gp + 16);

    auto tr_x_xxzz_z = pbuffer.data(idx_dip_gp + 17);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_x_xx_x,   \
                             tr_x_xx_y,   \
                             tr_x_xxz_x,  \
                             tr_x_xxz_y,  \
                             tr_x_xxzz_x, \
                             tr_x_xxzz_y, \
                             tr_x_xxzz_z, \
                             tr_x_xzz_z,  \
                             tr_x_zz_z,   \
                             ts_xzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzz_x[i] = tr_x_xx_x[i] * fe_0 + tr_x_xxz_x[i] * pa_z[i];

        tr_x_xxzz_y[i] = tr_x_xx_y[i] * fe_0 + tr_x_xxz_y[i] * pa_z[i];

        tr_x_xxzz_z[i] = tr_x_zz_z[i] * fe_0 + ts_xzz_z[i] * fe_0 + tr_x_xzz_z[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto tr_x_xyyy_x = pbuffer.data(idx_dip_gp + 18);

    auto tr_x_xyyy_y = pbuffer.data(idx_dip_gp + 19);

    auto tr_x_xyyy_z = pbuffer.data(idx_dip_gp + 20);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_x_xy_x, tr_x_xyy_x, tr_x_xyyy_x, tr_x_xyyy_y, tr_x_xyyy_z, tr_x_yyy_y, tr_x_yyy_z, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyy_x[i] = 2.0 * tr_x_xy_x[i] * fe_0 + tr_x_xyy_x[i] * pa_y[i];

        tr_x_xyyy_y[i] = ts_yyy_y[i] * fe_0 + tr_x_yyy_y[i] * pa_x[i];

        tr_x_xyyy_z[i] = ts_yyy_z[i] * fe_0 + tr_x_yyy_z[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto tr_x_xyyz_x = pbuffer.data(idx_dip_gp + 21);

    auto tr_x_xyyz_y = pbuffer.data(idx_dip_gp + 22);

    auto tr_x_xyyz_z = pbuffer.data(idx_dip_gp + 23);

#pragma omp simd aligned(pa_x, pa_z, tr_x_xyy_x, tr_x_xyy_y, tr_x_xyyz_x, tr_x_xyyz_y, tr_x_xyyz_z, tr_x_yyz_z, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyz_x[i] = tr_x_xyy_x[i] * pa_z[i];

        tr_x_xyyz_y[i] = tr_x_xyy_y[i] * pa_z[i];

        tr_x_xyyz_z[i] = ts_yyz_z[i] * fe_0 + tr_x_yyz_z[i] * pa_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto tr_x_xyzz_x = pbuffer.data(idx_dip_gp + 24);

    auto tr_x_xyzz_y = pbuffer.data(idx_dip_gp + 25);

    auto tr_x_xyzz_z = pbuffer.data(idx_dip_gp + 26);

#pragma omp simd aligned(pa_x, pa_y, tr_x_xyzz_x, tr_x_xyzz_y, tr_x_xyzz_z, tr_x_xzz_x, tr_x_xzz_z, tr_x_yzz_y, ts_yzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzz_x[i] = tr_x_xzz_x[i] * pa_y[i];

        tr_x_xyzz_y[i] = ts_yzz_y[i] * fe_0 + tr_x_yzz_y[i] * pa_x[i];

        tr_x_xyzz_z[i] = tr_x_xzz_z[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto tr_x_xzzz_x = pbuffer.data(idx_dip_gp + 27);

    auto tr_x_xzzz_y = pbuffer.data(idx_dip_gp + 28);

    auto tr_x_xzzz_z = pbuffer.data(idx_dip_gp + 29);

#pragma omp simd aligned( \
        pa_x, pa_z, tr_x_xz_x, tr_x_xzz_x, tr_x_xzzz_x, tr_x_xzzz_y, tr_x_xzzz_z, tr_x_zzz_y, tr_x_zzz_z, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzz_x[i] = 2.0 * tr_x_xz_x[i] * fe_0 + tr_x_xzz_x[i] * pa_z[i];

        tr_x_xzzz_y[i] = ts_zzz_y[i] * fe_0 + tr_x_zzz_y[i] * pa_x[i];

        tr_x_xzzz_z[i] = ts_zzz_z[i] * fe_0 + tr_x_zzz_z[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto tr_x_yyyy_x = pbuffer.data(idx_dip_gp + 30);

    auto tr_x_yyyy_y = pbuffer.data(idx_dip_gp + 31);

    auto tr_x_yyyy_z = pbuffer.data(idx_dip_gp + 32);

#pragma omp simd aligned( \
        pa_y, tr_x_yy_x, tr_x_yy_y, tr_x_yy_z, tr_x_yyy_0, tr_x_yyy_x, tr_x_yyy_y, tr_x_yyy_z, tr_x_yyyy_x, tr_x_yyyy_y, tr_x_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyy_x[i] = 3.0 * tr_x_yy_x[i] * fe_0 + tr_x_yyy_x[i] * pa_y[i];

        tr_x_yyyy_y[i] = 3.0 * tr_x_yy_y[i] * fe_0 + tr_x_yyy_0[i] * fe_0 + tr_x_yyy_y[i] * pa_y[i];

        tr_x_yyyy_z[i] = 3.0 * tr_x_yy_z[i] * fe_0 + tr_x_yyy_z[i] * pa_y[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto tr_x_yyyz_x = pbuffer.data(idx_dip_gp + 33);

    auto tr_x_yyyz_y = pbuffer.data(idx_dip_gp + 34);

    auto tr_x_yyyz_z = pbuffer.data(idx_dip_gp + 35);

#pragma omp simd aligned(pa_y, pa_z, tr_x_yyy_x, tr_x_yyy_y, tr_x_yyyz_x, tr_x_yyyz_y, tr_x_yyyz_z, tr_x_yyz_z, tr_x_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyz_x[i] = tr_x_yyy_x[i] * pa_z[i];

        tr_x_yyyz_y[i] = tr_x_yyy_y[i] * pa_z[i];

        tr_x_yyyz_z[i] = 2.0 * tr_x_yz_z[i] * fe_0 + tr_x_yyz_z[i] * pa_y[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto tr_x_yyzz_x = pbuffer.data(idx_dip_gp + 36);

    auto tr_x_yyzz_y = pbuffer.data(idx_dip_gp + 37);

    auto tr_x_yyzz_z = pbuffer.data(idx_dip_gp + 38);

#pragma omp simd aligned( \
        pa_y, pa_z, tr_x_yy_y, tr_x_yyz_y, tr_x_yyzz_x, tr_x_yyzz_y, tr_x_yyzz_z, tr_x_yzz_x, tr_x_yzz_z, tr_x_zz_x, tr_x_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzz_x[i] = tr_x_zz_x[i] * fe_0 + tr_x_yzz_x[i] * pa_y[i];

        tr_x_yyzz_y[i] = tr_x_yy_y[i] * fe_0 + tr_x_yyz_y[i] * pa_z[i];

        tr_x_yyzz_z[i] = tr_x_zz_z[i] * fe_0 + tr_x_yzz_z[i] * pa_y[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto tr_x_yzzz_x = pbuffer.data(idx_dip_gp + 39);

    auto tr_x_yzzz_y = pbuffer.data(idx_dip_gp + 40);

    auto tr_x_yzzz_z = pbuffer.data(idx_dip_gp + 41);

#pragma omp simd aligned(pa_y, tr_x_yzzz_x, tr_x_yzzz_y, tr_x_yzzz_z, tr_x_zzz_0, tr_x_zzz_x, tr_x_zzz_y, tr_x_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzz_x[i] = tr_x_zzz_x[i] * pa_y[i];

        tr_x_yzzz_y[i] = tr_x_zzz_0[i] * fe_0 + tr_x_zzz_y[i] * pa_y[i];

        tr_x_yzzz_z[i] = tr_x_zzz_z[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto tr_x_zzzz_x = pbuffer.data(idx_dip_gp + 42);

    auto tr_x_zzzz_y = pbuffer.data(idx_dip_gp + 43);

    auto tr_x_zzzz_z = pbuffer.data(idx_dip_gp + 44);

#pragma omp simd aligned( \
        pa_z, tr_x_zz_x, tr_x_zz_y, tr_x_zz_z, tr_x_zzz_0, tr_x_zzz_x, tr_x_zzz_y, tr_x_zzz_z, tr_x_zzzz_x, tr_x_zzzz_y, tr_x_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzz_x[i] = 3.0 * tr_x_zz_x[i] * fe_0 + tr_x_zzz_x[i] * pa_z[i];

        tr_x_zzzz_y[i] = 3.0 * tr_x_zz_y[i] * fe_0 + tr_x_zzz_y[i] * pa_z[i];

        tr_x_zzzz_z[i] = 3.0 * tr_x_zz_z[i] * fe_0 + tr_x_zzz_0[i] * fe_0 + tr_x_zzz_z[i] * pa_z[i];
    }

    // Set up 45-48 components of targeted buffer : GP

    auto tr_y_xxxx_x = pbuffer.data(idx_dip_gp + 45);

    auto tr_y_xxxx_y = pbuffer.data(idx_dip_gp + 46);

    auto tr_y_xxxx_z = pbuffer.data(idx_dip_gp + 47);

#pragma omp simd aligned( \
        pa_x, tr_y_xx_x, tr_y_xx_y, tr_y_xx_z, tr_y_xxx_0, tr_y_xxx_x, tr_y_xxx_y, tr_y_xxx_z, tr_y_xxxx_x, tr_y_xxxx_y, tr_y_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxx_x[i] = 3.0 * tr_y_xx_x[i] * fe_0 + tr_y_xxx_0[i] * fe_0 + tr_y_xxx_x[i] * pa_x[i];

        tr_y_xxxx_y[i] = 3.0 * tr_y_xx_y[i] * fe_0 + tr_y_xxx_y[i] * pa_x[i];

        tr_y_xxxx_z[i] = 3.0 * tr_y_xx_z[i] * fe_0 + tr_y_xxx_z[i] * pa_x[i];
    }

    // Set up 48-51 components of targeted buffer : GP

    auto tr_y_xxxy_x = pbuffer.data(idx_dip_gp + 48);

    auto tr_y_xxxy_y = pbuffer.data(idx_dip_gp + 49);

    auto tr_y_xxxy_z = pbuffer.data(idx_dip_gp + 50);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_y_xxx_x, tr_y_xxxy_x, tr_y_xxxy_y, tr_y_xxxy_z, tr_y_xxy_y, tr_y_xxy_z, tr_y_xy_y, tr_y_xy_z, ts_xxx_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxy_x[i] = ts_xxx_x[i] * fe_0 + tr_y_xxx_x[i] * pa_y[i];

        tr_y_xxxy_y[i] = 2.0 * tr_y_xy_y[i] * fe_0 + tr_y_xxy_y[i] * pa_x[i];

        tr_y_xxxy_z[i] = 2.0 * tr_y_xy_z[i] * fe_0 + tr_y_xxy_z[i] * pa_x[i];
    }

    // Set up 51-54 components of targeted buffer : GP

    auto tr_y_xxxz_x = pbuffer.data(idx_dip_gp + 51);

    auto tr_y_xxxz_y = pbuffer.data(idx_dip_gp + 52);

    auto tr_y_xxxz_z = pbuffer.data(idx_dip_gp + 53);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxx_x, tr_y_xxx_y, tr_y_xxxz_x, tr_y_xxxz_y, tr_y_xxxz_z, tr_y_xxz_z, tr_y_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxz_x[i] = tr_y_xxx_x[i] * pa_z[i];

        tr_y_xxxz_y[i] = tr_y_xxx_y[i] * pa_z[i];

        tr_y_xxxz_z[i] = 2.0 * tr_y_xz_z[i] * fe_0 + tr_y_xxz_z[i] * pa_x[i];
    }

    // Set up 54-57 components of targeted buffer : GP

    auto tr_y_xxyy_x = pbuffer.data(idx_dip_gp + 54);

    auto tr_y_xxyy_y = pbuffer.data(idx_dip_gp + 55);

    auto tr_y_xxyy_z = pbuffer.data(idx_dip_gp + 56);

#pragma omp simd aligned( \
        pa_x, tr_y_xxyy_x, tr_y_xxyy_y, tr_y_xxyy_z, tr_y_xyy_0, tr_y_xyy_x, tr_y_xyy_y, tr_y_xyy_z, tr_y_yy_x, tr_y_yy_y, tr_y_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyy_x[i] = tr_y_yy_x[i] * fe_0 + tr_y_xyy_0[i] * fe_0 + tr_y_xyy_x[i] * pa_x[i];

        tr_y_xxyy_y[i] = tr_y_yy_y[i] * fe_0 + tr_y_xyy_y[i] * pa_x[i];

        tr_y_xxyy_z[i] = tr_y_yy_z[i] * fe_0 + tr_y_xyy_z[i] * pa_x[i];
    }

    // Set up 57-60 components of targeted buffer : GP

    auto tr_y_xxyz_x = pbuffer.data(idx_dip_gp + 57);

    auto tr_y_xxyz_y = pbuffer.data(idx_dip_gp + 58);

    auto tr_y_xxyz_z = pbuffer.data(idx_dip_gp + 59);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxy_x, tr_y_xxy_y, tr_y_xxyz_x, tr_y_xxyz_y, tr_y_xxyz_z, tr_y_xyz_z, tr_y_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyz_x[i] = tr_y_xxy_x[i] * pa_z[i];

        tr_y_xxyz_y[i] = tr_y_xxy_y[i] * pa_z[i];

        tr_y_xxyz_z[i] = tr_y_yz_z[i] * fe_0 + tr_y_xyz_z[i] * pa_x[i];
    }

    // Set up 60-63 components of targeted buffer : GP

    auto tr_y_xxzz_x = pbuffer.data(idx_dip_gp + 60);

    auto tr_y_xxzz_y = pbuffer.data(idx_dip_gp + 61);

    auto tr_y_xxzz_z = pbuffer.data(idx_dip_gp + 62);

#pragma omp simd aligned( \
        pa_x, pa_z, tr_y_xx_x, tr_y_xxz_x, tr_y_xxzz_x, tr_y_xxzz_y, tr_y_xxzz_z, tr_y_xzz_y, tr_y_xzz_z, tr_y_zz_y, tr_y_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzz_x[i] = tr_y_xx_x[i] * fe_0 + tr_y_xxz_x[i] * pa_z[i];

        tr_y_xxzz_y[i] = tr_y_zz_y[i] * fe_0 + tr_y_xzz_y[i] * pa_x[i];

        tr_y_xxzz_z[i] = tr_y_zz_z[i] * fe_0 + tr_y_xzz_z[i] * pa_x[i];
    }

    // Set up 63-66 components of targeted buffer : GP

    auto tr_y_xyyy_x = pbuffer.data(idx_dip_gp + 63);

    auto tr_y_xyyy_y = pbuffer.data(idx_dip_gp + 64);

    auto tr_y_xyyy_z = pbuffer.data(idx_dip_gp + 65);

#pragma omp simd aligned(pa_x, tr_y_xyyy_x, tr_y_xyyy_y, tr_y_xyyy_z, tr_y_yyy_0, tr_y_yyy_x, tr_y_yyy_y, tr_y_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyy_x[i] = tr_y_yyy_0[i] * fe_0 + tr_y_yyy_x[i] * pa_x[i];

        tr_y_xyyy_y[i] = tr_y_yyy_y[i] * pa_x[i];

        tr_y_xyyy_z[i] = tr_y_yyy_z[i] * pa_x[i];
    }

    // Set up 66-69 components of targeted buffer : GP

    auto tr_y_xyyz_x = pbuffer.data(idx_dip_gp + 66);

    auto tr_y_xyyz_y = pbuffer.data(idx_dip_gp + 67);

    auto tr_y_xyyz_z = pbuffer.data(idx_dip_gp + 68);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xyy_x, tr_y_xyyz_x, tr_y_xyyz_y, tr_y_xyyz_z, tr_y_yyz_y, tr_y_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_y_xyyz_x[i] = tr_y_xyy_x[i] * pa_z[i];

        tr_y_xyyz_y[i] = tr_y_yyz_y[i] * pa_x[i];

        tr_y_xyyz_z[i] = tr_y_yyz_z[i] * pa_x[i];
    }

    // Set up 69-72 components of targeted buffer : GP

    auto tr_y_xyzz_x = pbuffer.data(idx_dip_gp + 69);

    auto tr_y_xyzz_y = pbuffer.data(idx_dip_gp + 70);

    auto tr_y_xyzz_z = pbuffer.data(idx_dip_gp + 71);

#pragma omp simd aligned(pa_x, tr_y_xyzz_x, tr_y_xyzz_y, tr_y_xyzz_z, tr_y_yzz_0, tr_y_yzz_x, tr_y_yzz_y, tr_y_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzz_x[i] = tr_y_yzz_0[i] * fe_0 + tr_y_yzz_x[i] * pa_x[i];

        tr_y_xyzz_y[i] = tr_y_yzz_y[i] * pa_x[i];

        tr_y_xyzz_z[i] = tr_y_yzz_z[i] * pa_x[i];
    }

    // Set up 72-75 components of targeted buffer : GP

    auto tr_y_xzzz_x = pbuffer.data(idx_dip_gp + 72);

    auto tr_y_xzzz_y = pbuffer.data(idx_dip_gp + 73);

    auto tr_y_xzzz_z = pbuffer.data(idx_dip_gp + 74);

#pragma omp simd aligned(pa_x, tr_y_xzzz_x, tr_y_xzzz_y, tr_y_xzzz_z, tr_y_zzz_0, tr_y_zzz_x, tr_y_zzz_y, tr_y_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzz_x[i] = tr_y_zzz_0[i] * fe_0 + tr_y_zzz_x[i] * pa_x[i];

        tr_y_xzzz_y[i] = tr_y_zzz_y[i] * pa_x[i];

        tr_y_xzzz_z[i] = tr_y_zzz_z[i] * pa_x[i];
    }

    // Set up 75-78 components of targeted buffer : GP

    auto tr_y_yyyy_x = pbuffer.data(idx_dip_gp + 75);

    auto tr_y_yyyy_y = pbuffer.data(idx_dip_gp + 76);

    auto tr_y_yyyy_z = pbuffer.data(idx_dip_gp + 77);

#pragma omp simd aligned(pa_y,            \
                             tr_y_yy_x,   \
                             tr_y_yy_y,   \
                             tr_y_yy_z,   \
                             tr_y_yyy_0,  \
                             tr_y_yyy_x,  \
                             tr_y_yyy_y,  \
                             tr_y_yyy_z,  \
                             tr_y_yyyy_x, \
                             tr_y_yyyy_y, \
                             tr_y_yyyy_z, \
                             ts_yyy_x,    \
                             ts_yyy_y,    \
                             ts_yyy_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyy_x[i] = 3.0 * tr_y_yy_x[i] * fe_0 + ts_yyy_x[i] * fe_0 + tr_y_yyy_x[i] * pa_y[i];

        tr_y_yyyy_y[i] = 3.0 * tr_y_yy_y[i] * fe_0 + tr_y_yyy_0[i] * fe_0 + ts_yyy_y[i] * fe_0 + tr_y_yyy_y[i] * pa_y[i];

        tr_y_yyyy_z[i] = 3.0 * tr_y_yy_z[i] * fe_0 + ts_yyy_z[i] * fe_0 + tr_y_yyy_z[i] * pa_y[i];
    }

    // Set up 78-81 components of targeted buffer : GP

    auto tr_y_yyyz_x = pbuffer.data(idx_dip_gp + 78);

    auto tr_y_yyyz_y = pbuffer.data(idx_dip_gp + 79);

    auto tr_y_yyyz_z = pbuffer.data(idx_dip_gp + 80);

#pragma omp simd aligned(pa_z, tr_y_yyy_0, tr_y_yyy_x, tr_y_yyy_y, tr_y_yyy_z, tr_y_yyyz_x, tr_y_yyyz_y, tr_y_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyz_x[i] = tr_y_yyy_x[i] * pa_z[i];

        tr_y_yyyz_y[i] = tr_y_yyy_y[i] * pa_z[i];

        tr_y_yyyz_z[i] = tr_y_yyy_0[i] * fe_0 + tr_y_yyy_z[i] * pa_z[i];
    }

    // Set up 81-84 components of targeted buffer : GP

    auto tr_y_yyzz_x = pbuffer.data(idx_dip_gp + 81);

    auto tr_y_yyzz_y = pbuffer.data(idx_dip_gp + 82);

    auto tr_y_yyzz_z = pbuffer.data(idx_dip_gp + 83);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tr_y_yy_x,   \
                             tr_y_yy_y,   \
                             tr_y_yyz_x,  \
                             tr_y_yyz_y,  \
                             tr_y_yyzz_x, \
                             tr_y_yyzz_y, \
                             tr_y_yyzz_z, \
                             tr_y_yzz_z,  \
                             tr_y_zz_z,   \
                             ts_yzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzz_x[i] = tr_y_yy_x[i] * fe_0 + tr_y_yyz_x[i] * pa_z[i];

        tr_y_yyzz_y[i] = tr_y_yy_y[i] * fe_0 + tr_y_yyz_y[i] * pa_z[i];

        tr_y_yyzz_z[i] = tr_y_zz_z[i] * fe_0 + ts_yzz_z[i] * fe_0 + tr_y_yzz_z[i] * pa_y[i];
    }

    // Set up 84-87 components of targeted buffer : GP

    auto tr_y_yzzz_x = pbuffer.data(idx_dip_gp + 84);

    auto tr_y_yzzz_y = pbuffer.data(idx_dip_gp + 85);

    auto tr_y_yzzz_z = pbuffer.data(idx_dip_gp + 86);

#pragma omp simd aligned( \
        pa_y, pa_z, tr_y_yz_y, tr_y_yzz_y, tr_y_yzzz_x, tr_y_yzzz_y, tr_y_yzzz_z, tr_y_zzz_x, tr_y_zzz_z, ts_zzz_x, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzz_x[i] = ts_zzz_x[i] * fe_0 + tr_y_zzz_x[i] * pa_y[i];

        tr_y_yzzz_y[i] = 2.0 * tr_y_yz_y[i] * fe_0 + tr_y_yzz_y[i] * pa_z[i];

        tr_y_yzzz_z[i] = ts_zzz_z[i] * fe_0 + tr_y_zzz_z[i] * pa_y[i];
    }

    // Set up 87-90 components of targeted buffer : GP

    auto tr_y_zzzz_x = pbuffer.data(idx_dip_gp + 87);

    auto tr_y_zzzz_y = pbuffer.data(idx_dip_gp + 88);

    auto tr_y_zzzz_z = pbuffer.data(idx_dip_gp + 89);

#pragma omp simd aligned( \
        pa_z, tr_y_zz_x, tr_y_zz_y, tr_y_zz_z, tr_y_zzz_0, tr_y_zzz_x, tr_y_zzz_y, tr_y_zzz_z, tr_y_zzzz_x, tr_y_zzzz_y, tr_y_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzz_x[i] = 3.0 * tr_y_zz_x[i] * fe_0 + tr_y_zzz_x[i] * pa_z[i];

        tr_y_zzzz_y[i] = 3.0 * tr_y_zz_y[i] * fe_0 + tr_y_zzz_y[i] * pa_z[i];

        tr_y_zzzz_z[i] = 3.0 * tr_y_zz_z[i] * fe_0 + tr_y_zzz_0[i] * fe_0 + tr_y_zzz_z[i] * pa_z[i];
    }

    // Set up 90-93 components of targeted buffer : GP

    auto tr_z_xxxx_x = pbuffer.data(idx_dip_gp + 90);

    auto tr_z_xxxx_y = pbuffer.data(idx_dip_gp + 91);

    auto tr_z_xxxx_z = pbuffer.data(idx_dip_gp + 92);

#pragma omp simd aligned( \
        pa_x, tr_z_xx_x, tr_z_xx_y, tr_z_xx_z, tr_z_xxx_0, tr_z_xxx_x, tr_z_xxx_y, tr_z_xxx_z, tr_z_xxxx_x, tr_z_xxxx_y, tr_z_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxx_x[i] = 3.0 * tr_z_xx_x[i] * fe_0 + tr_z_xxx_0[i] * fe_0 + tr_z_xxx_x[i] * pa_x[i];

        tr_z_xxxx_y[i] = 3.0 * tr_z_xx_y[i] * fe_0 + tr_z_xxx_y[i] * pa_x[i];

        tr_z_xxxx_z[i] = 3.0 * tr_z_xx_z[i] * fe_0 + tr_z_xxx_z[i] * pa_x[i];
    }

    // Set up 93-96 components of targeted buffer : GP

    auto tr_z_xxxy_x = pbuffer.data(idx_dip_gp + 93);

    auto tr_z_xxxy_y = pbuffer.data(idx_dip_gp + 94);

    auto tr_z_xxxy_z = pbuffer.data(idx_dip_gp + 95);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxx_x, tr_z_xxx_z, tr_z_xxxy_x, tr_z_xxxy_y, tr_z_xxxy_z, tr_z_xxy_y, tr_z_xy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxy_x[i] = tr_z_xxx_x[i] * pa_y[i];

        tr_z_xxxy_y[i] = 2.0 * tr_z_xy_y[i] * fe_0 + tr_z_xxy_y[i] * pa_x[i];

        tr_z_xxxy_z[i] = tr_z_xxx_z[i] * pa_y[i];
    }

    // Set up 96-99 components of targeted buffer : GP

    auto tr_z_xxxz_x = pbuffer.data(idx_dip_gp + 96);

    auto tr_z_xxxz_y = pbuffer.data(idx_dip_gp + 97);

    auto tr_z_xxxz_z = pbuffer.data(idx_dip_gp + 98);

#pragma omp simd aligned( \
        pa_x, pa_z, tr_z_xxx_x, tr_z_xxxz_x, tr_z_xxxz_y, tr_z_xxxz_z, tr_z_xxz_y, tr_z_xxz_z, tr_z_xz_y, tr_z_xz_z, ts_xxx_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxz_x[i] = ts_xxx_x[i] * fe_0 + tr_z_xxx_x[i] * pa_z[i];

        tr_z_xxxz_y[i] = 2.0 * tr_z_xz_y[i] * fe_0 + tr_z_xxz_y[i] * pa_x[i];

        tr_z_xxxz_z[i] = 2.0 * tr_z_xz_z[i] * fe_0 + tr_z_xxz_z[i] * pa_x[i];
    }

    // Set up 99-102 components of targeted buffer : GP

    auto tr_z_xxyy_x = pbuffer.data(idx_dip_gp + 99);

    auto tr_z_xxyy_y = pbuffer.data(idx_dip_gp + 100);

    auto tr_z_xxyy_z = pbuffer.data(idx_dip_gp + 101);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_z_xx_x, tr_z_xxy_x, tr_z_xxyy_x, tr_z_xxyy_y, tr_z_xxyy_z, tr_z_xyy_y, tr_z_xyy_z, tr_z_yy_y, tr_z_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyy_x[i] = tr_z_xx_x[i] * fe_0 + tr_z_xxy_x[i] * pa_y[i];

        tr_z_xxyy_y[i] = tr_z_yy_y[i] * fe_0 + tr_z_xyy_y[i] * pa_x[i];

        tr_z_xxyy_z[i] = tr_z_yy_z[i] * fe_0 + tr_z_xyy_z[i] * pa_x[i];
    }

    // Set up 102-105 components of targeted buffer : GP

    auto tr_z_xxyz_x = pbuffer.data(idx_dip_gp + 102);

    auto tr_z_xxyz_y = pbuffer.data(idx_dip_gp + 103);

    auto tr_z_xxyz_z = pbuffer.data(idx_dip_gp + 104);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxyz_x, tr_z_xxyz_y, tr_z_xxyz_z, tr_z_xxz_x, tr_z_xxz_z, tr_z_xyz_y, tr_z_yz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyz_x[i] = tr_z_xxz_x[i] * pa_y[i];

        tr_z_xxyz_y[i] = tr_z_yz_y[i] * fe_0 + tr_z_xyz_y[i] * pa_x[i];

        tr_z_xxyz_z[i] = tr_z_xxz_z[i] * pa_y[i];
    }

    // Set up 105-108 components of targeted buffer : GP

    auto tr_z_xxzz_x = pbuffer.data(idx_dip_gp + 105);

    auto tr_z_xxzz_y = pbuffer.data(idx_dip_gp + 106);

    auto tr_z_xxzz_z = pbuffer.data(idx_dip_gp + 107);

#pragma omp simd aligned( \
        pa_x, tr_z_xxzz_x, tr_z_xxzz_y, tr_z_xxzz_z, tr_z_xzz_0, tr_z_xzz_x, tr_z_xzz_y, tr_z_xzz_z, tr_z_zz_x, tr_z_zz_y, tr_z_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzz_x[i] = tr_z_zz_x[i] * fe_0 + tr_z_xzz_0[i] * fe_0 + tr_z_xzz_x[i] * pa_x[i];

        tr_z_xxzz_y[i] = tr_z_zz_y[i] * fe_0 + tr_z_xzz_y[i] * pa_x[i];

        tr_z_xxzz_z[i] = tr_z_zz_z[i] * fe_0 + tr_z_xzz_z[i] * pa_x[i];
    }

    // Set up 108-111 components of targeted buffer : GP

    auto tr_z_xyyy_x = pbuffer.data(idx_dip_gp + 108);

    auto tr_z_xyyy_y = pbuffer.data(idx_dip_gp + 109);

    auto tr_z_xyyy_z = pbuffer.data(idx_dip_gp + 110);

#pragma omp simd aligned(pa_x, tr_z_xyyy_x, tr_z_xyyy_y, tr_z_xyyy_z, tr_z_yyy_0, tr_z_yyy_x, tr_z_yyy_y, tr_z_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyy_x[i] = tr_z_yyy_0[i] * fe_0 + tr_z_yyy_x[i] * pa_x[i];

        tr_z_xyyy_y[i] = tr_z_yyy_y[i] * pa_x[i];

        tr_z_xyyy_z[i] = tr_z_yyy_z[i] * pa_x[i];
    }

    // Set up 111-114 components of targeted buffer : GP

    auto tr_z_xyyz_x = pbuffer.data(idx_dip_gp + 111);

    auto tr_z_xyyz_y = pbuffer.data(idx_dip_gp + 112);

    auto tr_z_xyyz_z = pbuffer.data(idx_dip_gp + 113);

#pragma omp simd aligned(pa_x, tr_z_xyyz_x, tr_z_xyyz_y, tr_z_xyyz_z, tr_z_yyz_0, tr_z_yyz_x, tr_z_yyz_y, tr_z_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyz_x[i] = tr_z_yyz_0[i] * fe_0 + tr_z_yyz_x[i] * pa_x[i];

        tr_z_xyyz_y[i] = tr_z_yyz_y[i] * pa_x[i];

        tr_z_xyyz_z[i] = tr_z_yyz_z[i] * pa_x[i];
    }

    // Set up 114-117 components of targeted buffer : GP

    auto tr_z_xyzz_x = pbuffer.data(idx_dip_gp + 114);

    auto tr_z_xyzz_y = pbuffer.data(idx_dip_gp + 115);

    auto tr_z_xyzz_z = pbuffer.data(idx_dip_gp + 116);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xyzz_x, tr_z_xyzz_y, tr_z_xyzz_z, tr_z_xzz_x, tr_z_yzz_y, tr_z_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_z_xyzz_x[i] = tr_z_xzz_x[i] * pa_y[i];

        tr_z_xyzz_y[i] = tr_z_yzz_y[i] * pa_x[i];

        tr_z_xyzz_z[i] = tr_z_yzz_z[i] * pa_x[i];
    }

    // Set up 117-120 components of targeted buffer : GP

    auto tr_z_xzzz_x = pbuffer.data(idx_dip_gp + 117);

    auto tr_z_xzzz_y = pbuffer.data(idx_dip_gp + 118);

    auto tr_z_xzzz_z = pbuffer.data(idx_dip_gp + 119);

#pragma omp simd aligned(pa_x, tr_z_xzzz_x, tr_z_xzzz_y, tr_z_xzzz_z, tr_z_zzz_0, tr_z_zzz_x, tr_z_zzz_y, tr_z_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzz_x[i] = tr_z_zzz_0[i] * fe_0 + tr_z_zzz_x[i] * pa_x[i];

        tr_z_xzzz_y[i] = tr_z_zzz_y[i] * pa_x[i];

        tr_z_xzzz_z[i] = tr_z_zzz_z[i] * pa_x[i];
    }

    // Set up 120-123 components of targeted buffer : GP

    auto tr_z_yyyy_x = pbuffer.data(idx_dip_gp + 120);

    auto tr_z_yyyy_y = pbuffer.data(idx_dip_gp + 121);

    auto tr_z_yyyy_z = pbuffer.data(idx_dip_gp + 122);

#pragma omp simd aligned( \
        pa_y, tr_z_yy_x, tr_z_yy_y, tr_z_yy_z, tr_z_yyy_0, tr_z_yyy_x, tr_z_yyy_y, tr_z_yyy_z, tr_z_yyyy_x, tr_z_yyyy_y, tr_z_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyy_x[i] = 3.0 * tr_z_yy_x[i] * fe_0 + tr_z_yyy_x[i] * pa_y[i];

        tr_z_yyyy_y[i] = 3.0 * tr_z_yy_y[i] * fe_0 + tr_z_yyy_0[i] * fe_0 + tr_z_yyy_y[i] * pa_y[i];

        tr_z_yyyy_z[i] = 3.0 * tr_z_yy_z[i] * fe_0 + tr_z_yyy_z[i] * pa_y[i];
    }

    // Set up 123-126 components of targeted buffer : GP

    auto tr_z_yyyz_x = pbuffer.data(idx_dip_gp + 123);

    auto tr_z_yyyz_y = pbuffer.data(idx_dip_gp + 124);

    auto tr_z_yyyz_z = pbuffer.data(idx_dip_gp + 125);

#pragma omp simd aligned( \
        pa_y, pa_z, tr_z_yyy_y, tr_z_yyyz_x, tr_z_yyyz_y, tr_z_yyyz_z, tr_z_yyz_x, tr_z_yyz_z, tr_z_yz_x, tr_z_yz_z, ts_yyy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyz_x[i] = 2.0 * tr_z_yz_x[i] * fe_0 + tr_z_yyz_x[i] * pa_y[i];

        tr_z_yyyz_y[i] = ts_yyy_y[i] * fe_0 + tr_z_yyy_y[i] * pa_z[i];

        tr_z_yyyz_z[i] = 2.0 * tr_z_yz_z[i] * fe_0 + tr_z_yyz_z[i] * pa_y[i];
    }

    // Set up 126-129 components of targeted buffer : GP

    auto tr_z_yyzz_x = pbuffer.data(idx_dip_gp + 126);

    auto tr_z_yyzz_y = pbuffer.data(idx_dip_gp + 127);

    auto tr_z_yyzz_z = pbuffer.data(idx_dip_gp + 128);

#pragma omp simd aligned( \
        pa_y, tr_z_yyzz_x, tr_z_yyzz_y, tr_z_yyzz_z, tr_z_yzz_0, tr_z_yzz_x, tr_z_yzz_y, tr_z_yzz_z, tr_z_zz_x, tr_z_zz_y, tr_z_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzz_x[i] = tr_z_zz_x[i] * fe_0 + tr_z_yzz_x[i] * pa_y[i];

        tr_z_yyzz_y[i] = tr_z_zz_y[i] * fe_0 + tr_z_yzz_0[i] * fe_0 + tr_z_yzz_y[i] * pa_y[i];

        tr_z_yyzz_z[i] = tr_z_zz_z[i] * fe_0 + tr_z_yzz_z[i] * pa_y[i];
    }

    // Set up 129-132 components of targeted buffer : GP

    auto tr_z_yzzz_x = pbuffer.data(idx_dip_gp + 129);

    auto tr_z_yzzz_y = pbuffer.data(idx_dip_gp + 130);

    auto tr_z_yzzz_z = pbuffer.data(idx_dip_gp + 131);

#pragma omp simd aligned(pa_y, tr_z_yzzz_x, tr_z_yzzz_y, tr_z_yzzz_z, tr_z_zzz_0, tr_z_zzz_x, tr_z_zzz_y, tr_z_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzz_x[i] = tr_z_zzz_x[i] * pa_y[i];

        tr_z_yzzz_y[i] = tr_z_zzz_0[i] * fe_0 + tr_z_zzz_y[i] * pa_y[i];

        tr_z_yzzz_z[i] = tr_z_zzz_z[i] * pa_y[i];
    }

    // Set up 132-135 components of targeted buffer : GP

    auto tr_z_zzzz_x = pbuffer.data(idx_dip_gp + 132);

    auto tr_z_zzzz_y = pbuffer.data(idx_dip_gp + 133);

    auto tr_z_zzzz_z = pbuffer.data(idx_dip_gp + 134);

#pragma omp simd aligned(pa_z,            \
                             tr_z_zz_x,   \
                             tr_z_zz_y,   \
                             tr_z_zz_z,   \
                             tr_z_zzz_0,  \
                             tr_z_zzz_x,  \
                             tr_z_zzz_y,  \
                             tr_z_zzz_z,  \
                             tr_z_zzzz_x, \
                             tr_z_zzzz_y, \
                             tr_z_zzzz_z, \
                             ts_zzz_x,    \
                             ts_zzz_y,    \
                             ts_zzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzz_x[i] = 3.0 * tr_z_zz_x[i] * fe_0 + ts_zzz_x[i] * fe_0 + tr_z_zzz_x[i] * pa_z[i];

        tr_z_zzzz_y[i] = 3.0 * tr_z_zz_y[i] * fe_0 + ts_zzz_y[i] * fe_0 + tr_z_zzz_y[i] * pa_z[i];

        tr_z_zzzz_z[i] = 3.0 * tr_z_zz_z[i] * fe_0 + tr_z_zzz_0[i] * fe_0 + ts_zzz_z[i] * fe_0 + tr_z_zzz_z[i] * pa_z[i];
    }
}

}  // namespace diprec
