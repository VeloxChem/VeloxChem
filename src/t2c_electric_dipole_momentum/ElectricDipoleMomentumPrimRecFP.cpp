#include "ElectricDipoleMomentumPrimRecFP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_fp(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_fp,
                                      const size_t              idx_dip_pp,
                                      const size_t              idx_dip_ds,
                                      const size_t              idx_ovl_dp,
                                      const size_t              idx_dip_dp,
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

    // Set up components of auxiliary buffer : DS

    auto tr_x_xx_0 = pbuffer.data(idx_dip_ds);

    auto tr_x_yy_0 = pbuffer.data(idx_dip_ds + 3);

    auto tr_x_zz_0 = pbuffer.data(idx_dip_ds + 5);

    auto tr_y_xx_0 = pbuffer.data(idx_dip_ds + 6);

    auto tr_y_yy_0 = pbuffer.data(idx_dip_ds + 9);

    auto tr_y_zz_0 = pbuffer.data(idx_dip_ds + 11);

    auto tr_z_xx_0 = pbuffer.data(idx_dip_ds + 12);

    auto tr_z_yy_0 = pbuffer.data(idx_dip_ds + 15);

    auto tr_z_zz_0 = pbuffer.data(idx_dip_ds + 17);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_ovl_dp);

    auto ts_xx_y = pbuffer.data(idx_ovl_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_ovl_dp + 2);

    auto ts_yy_x = pbuffer.data(idx_ovl_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_ovl_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_ovl_dp + 11);

    auto ts_zz_x = pbuffer.data(idx_ovl_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_ovl_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_ovl_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto tr_x_xx_x = pbuffer.data(idx_dip_dp);

    auto tr_x_xx_y = pbuffer.data(idx_dip_dp + 1);

    auto tr_x_xx_z = pbuffer.data(idx_dip_dp + 2);

    auto tr_x_xy_x = pbuffer.data(idx_dip_dp + 3);

    auto tr_x_xy_y = pbuffer.data(idx_dip_dp + 4);

    auto tr_x_xz_x = pbuffer.data(idx_dip_dp + 6);

    auto tr_x_xz_z = pbuffer.data(idx_dip_dp + 8);

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

    auto tr_y_xy_x = pbuffer.data(idx_dip_dp + 21);

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

    auto tr_z_xz_x = pbuffer.data(idx_dip_dp + 42);

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

    // Set up 0-3 components of targeted buffer : FP

    auto tr_x_xxx_x = pbuffer.data(idx_dip_fp);

    auto tr_x_xxx_y = pbuffer.data(idx_dip_fp + 1);

    auto tr_x_xxx_z = pbuffer.data(idx_dip_fp + 2);

#pragma omp simd aligned(pa_x,           \
                             tr_x_x_x,   \
                             tr_x_x_y,   \
                             tr_x_x_z,   \
                             tr_x_xx_0,  \
                             tr_x_xx_x,  \
                             tr_x_xx_y,  \
                             tr_x_xx_z,  \
                             tr_x_xxx_x, \
                             tr_x_xxx_y, \
                             tr_x_xxx_z, \
                             ts_xx_x,    \
                             ts_xx_y,    \
                             ts_xx_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_x[i] = 2.0 * tr_x_x_x[i] * fe_0 + tr_x_xx_0[i] * fe_0 + ts_xx_x[i] * fe_0 + tr_x_xx_x[i] * pa_x[i];

        tr_x_xxx_y[i] = 2.0 * tr_x_x_y[i] * fe_0 + ts_xx_y[i] * fe_0 + tr_x_xx_y[i] * pa_x[i];

        tr_x_xxx_z[i] = 2.0 * tr_x_x_z[i] * fe_0 + ts_xx_z[i] * fe_0 + tr_x_xx_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto tr_x_xxy_x = pbuffer.data(idx_dip_fp + 3);

    auto tr_x_xxy_y = pbuffer.data(idx_dip_fp + 4);

    auto tr_x_xxy_z = pbuffer.data(idx_dip_fp + 5);

#pragma omp simd aligned(pa_y, tr_x_xx_0, tr_x_xx_x, tr_x_xx_y, tr_x_xx_z, tr_x_xxy_x, tr_x_xxy_y, tr_x_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxy_x[i] = tr_x_xx_x[i] * pa_y[i];

        tr_x_xxy_y[i] = tr_x_xx_0[i] * fe_0 + tr_x_xx_y[i] * pa_y[i];

        tr_x_xxy_z[i] = tr_x_xx_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto tr_x_xxz_x = pbuffer.data(idx_dip_fp + 6);

    auto tr_x_xxz_y = pbuffer.data(idx_dip_fp + 7);

    auto tr_x_xxz_z = pbuffer.data(idx_dip_fp + 8);

#pragma omp simd aligned(pa_z, tr_x_xx_0, tr_x_xx_x, tr_x_xx_y, tr_x_xx_z, tr_x_xxz_x, tr_x_xxz_y, tr_x_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxz_x[i] = tr_x_xx_x[i] * pa_z[i];

        tr_x_xxz_y[i] = tr_x_xx_y[i] * pa_z[i];

        tr_x_xxz_z[i] = tr_x_xx_0[i] * fe_0 + tr_x_xx_z[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto tr_x_xyy_x = pbuffer.data(idx_dip_fp + 9);

    auto tr_x_xyy_y = pbuffer.data(idx_dip_fp + 10);

    auto tr_x_xyy_z = pbuffer.data(idx_dip_fp + 11);

#pragma omp simd aligned(pa_x, pa_y, tr_x_x_x, tr_x_xy_x, tr_x_xyy_x, tr_x_xyy_y, tr_x_xyy_z, tr_x_yy_y, tr_x_yy_z, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyy_x[i] = tr_x_x_x[i] * fe_0 + tr_x_xy_x[i] * pa_y[i];

        tr_x_xyy_y[i] = ts_yy_y[i] * fe_0 + tr_x_yy_y[i] * pa_x[i];

        tr_x_xyy_z[i] = ts_yy_z[i] * fe_0 + tr_x_yy_z[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto tr_x_xyz_x = pbuffer.data(idx_dip_fp + 12);

    auto tr_x_xyz_y = pbuffer.data(idx_dip_fp + 13);

    auto tr_x_xyz_z = pbuffer.data(idx_dip_fp + 14);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xy_y, tr_x_xyz_x, tr_x_xyz_y, tr_x_xyz_z, tr_x_xz_x, tr_x_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_x_xyz_x[i] = tr_x_xz_x[i] * pa_y[i];

        tr_x_xyz_y[i] = tr_x_xy_y[i] * pa_z[i];

        tr_x_xyz_z[i] = tr_x_xz_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto tr_x_xzz_x = pbuffer.data(idx_dip_fp + 15);

    auto tr_x_xzz_y = pbuffer.data(idx_dip_fp + 16);

    auto tr_x_xzz_z = pbuffer.data(idx_dip_fp + 17);

#pragma omp simd aligned(pa_x, pa_z, tr_x_x_x, tr_x_xz_x, tr_x_xzz_x, tr_x_xzz_y, tr_x_xzz_z, tr_x_zz_y, tr_x_zz_z, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzz_x[i] = tr_x_x_x[i] * fe_0 + tr_x_xz_x[i] * pa_z[i];

        tr_x_xzz_y[i] = ts_zz_y[i] * fe_0 + tr_x_zz_y[i] * pa_x[i];

        tr_x_xzz_z[i] = ts_zz_z[i] * fe_0 + tr_x_zz_z[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto tr_x_yyy_x = pbuffer.data(idx_dip_fp + 18);

    auto tr_x_yyy_y = pbuffer.data(idx_dip_fp + 19);

    auto tr_x_yyy_z = pbuffer.data(idx_dip_fp + 20);

#pragma omp simd aligned( \
        pa_y, tr_x_y_x, tr_x_y_y, tr_x_y_z, tr_x_yy_0, tr_x_yy_x, tr_x_yy_y, tr_x_yy_z, tr_x_yyy_x, tr_x_yyy_y, tr_x_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyy_x[i] = 2.0 * tr_x_y_x[i] * fe_0 + tr_x_yy_x[i] * pa_y[i];

        tr_x_yyy_y[i] = 2.0 * tr_x_y_y[i] * fe_0 + tr_x_yy_0[i] * fe_0 + tr_x_yy_y[i] * pa_y[i];

        tr_x_yyy_z[i] = 2.0 * tr_x_y_z[i] * fe_0 + tr_x_yy_z[i] * pa_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto tr_x_yyz_x = pbuffer.data(idx_dip_fp + 21);

    auto tr_x_yyz_y = pbuffer.data(idx_dip_fp + 22);

    auto tr_x_yyz_z = pbuffer.data(idx_dip_fp + 23);

#pragma omp simd aligned(pa_y, pa_z, tr_x_yy_x, tr_x_yy_y, tr_x_yyz_x, tr_x_yyz_y, tr_x_yyz_z, tr_x_yz_z, tr_x_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyz_x[i] = tr_x_yy_x[i] * pa_z[i];

        tr_x_yyz_y[i] = tr_x_yy_y[i] * pa_z[i];

        tr_x_yyz_z[i] = tr_x_z_z[i] * fe_0 + tr_x_yz_z[i] * pa_y[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto tr_x_yzz_x = pbuffer.data(idx_dip_fp + 24);

    auto tr_x_yzz_y = pbuffer.data(idx_dip_fp + 25);

    auto tr_x_yzz_z = pbuffer.data(idx_dip_fp + 26);

#pragma omp simd aligned(pa_y, tr_x_yzz_x, tr_x_yzz_y, tr_x_yzz_z, tr_x_zz_0, tr_x_zz_x, tr_x_zz_y, tr_x_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzz_x[i] = tr_x_zz_x[i] * pa_y[i];

        tr_x_yzz_y[i] = tr_x_zz_0[i] * fe_0 + tr_x_zz_y[i] * pa_y[i];

        tr_x_yzz_z[i] = tr_x_zz_z[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto tr_x_zzz_x = pbuffer.data(idx_dip_fp + 27);

    auto tr_x_zzz_y = pbuffer.data(idx_dip_fp + 28);

    auto tr_x_zzz_z = pbuffer.data(idx_dip_fp + 29);

#pragma omp simd aligned( \
        pa_z, tr_x_z_x, tr_x_z_y, tr_x_z_z, tr_x_zz_0, tr_x_zz_x, tr_x_zz_y, tr_x_zz_z, tr_x_zzz_x, tr_x_zzz_y, tr_x_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzz_x[i] = 2.0 * tr_x_z_x[i] * fe_0 + tr_x_zz_x[i] * pa_z[i];

        tr_x_zzz_y[i] = 2.0 * tr_x_z_y[i] * fe_0 + tr_x_zz_y[i] * pa_z[i];

        tr_x_zzz_z[i] = 2.0 * tr_x_z_z[i] * fe_0 + tr_x_zz_0[i] * fe_0 + tr_x_zz_z[i] * pa_z[i];
    }

    // Set up 30-33 components of targeted buffer : FP

    auto tr_y_xxx_x = pbuffer.data(idx_dip_fp + 30);

    auto tr_y_xxx_y = pbuffer.data(idx_dip_fp + 31);

    auto tr_y_xxx_z = pbuffer.data(idx_dip_fp + 32);

#pragma omp simd aligned( \
        pa_x, tr_y_x_x, tr_y_x_y, tr_y_x_z, tr_y_xx_0, tr_y_xx_x, tr_y_xx_y, tr_y_xx_z, tr_y_xxx_x, tr_y_xxx_y, tr_y_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxx_x[i] = 2.0 * tr_y_x_x[i] * fe_0 + tr_y_xx_0[i] * fe_0 + tr_y_xx_x[i] * pa_x[i];

        tr_y_xxx_y[i] = 2.0 * tr_y_x_y[i] * fe_0 + tr_y_xx_y[i] * pa_x[i];

        tr_y_xxx_z[i] = 2.0 * tr_y_x_z[i] * fe_0 + tr_y_xx_z[i] * pa_x[i];
    }

    // Set up 33-36 components of targeted buffer : FP

    auto tr_y_xxy_x = pbuffer.data(idx_dip_fp + 33);

    auto tr_y_xxy_y = pbuffer.data(idx_dip_fp + 34);

    auto tr_y_xxy_z = pbuffer.data(idx_dip_fp + 35);

#pragma omp simd aligned(pa_x, pa_y, tr_y_xx_x, tr_y_xxy_x, tr_y_xxy_y, tr_y_xxy_z, tr_y_xy_y, tr_y_xy_z, tr_y_y_y, tr_y_y_z, ts_xx_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxy_x[i] = ts_xx_x[i] * fe_0 + tr_y_xx_x[i] * pa_y[i];

        tr_y_xxy_y[i] = tr_y_y_y[i] * fe_0 + tr_y_xy_y[i] * pa_x[i];

        tr_y_xxy_z[i] = tr_y_y_z[i] * fe_0 + tr_y_xy_z[i] * pa_x[i];
    }

    // Set up 36-39 components of targeted buffer : FP

    auto tr_y_xxz_x = pbuffer.data(idx_dip_fp + 36);

    auto tr_y_xxz_y = pbuffer.data(idx_dip_fp + 37);

    auto tr_y_xxz_z = pbuffer.data(idx_dip_fp + 38);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xx_x, tr_y_xx_y, tr_y_xxz_x, tr_y_xxz_y, tr_y_xxz_z, tr_y_xz_z, tr_y_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxz_x[i] = tr_y_xx_x[i] * pa_z[i];

        tr_y_xxz_y[i] = tr_y_xx_y[i] * pa_z[i];

        tr_y_xxz_z[i] = tr_y_z_z[i] * fe_0 + tr_y_xz_z[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : FP

    auto tr_y_xyy_x = pbuffer.data(idx_dip_fp + 39);

    auto tr_y_xyy_y = pbuffer.data(idx_dip_fp + 40);

    auto tr_y_xyy_z = pbuffer.data(idx_dip_fp + 41);

#pragma omp simd aligned(pa_x, tr_y_xyy_x, tr_y_xyy_y, tr_y_xyy_z, tr_y_yy_0, tr_y_yy_x, tr_y_yy_y, tr_y_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyy_x[i] = tr_y_yy_0[i] * fe_0 + tr_y_yy_x[i] * pa_x[i];

        tr_y_xyy_y[i] = tr_y_yy_y[i] * pa_x[i];

        tr_y_xyy_z[i] = tr_y_yy_z[i] * pa_x[i];
    }

    // Set up 42-45 components of targeted buffer : FP

    auto tr_y_xyz_x = pbuffer.data(idx_dip_fp + 42);

    auto tr_y_xyz_y = pbuffer.data(idx_dip_fp + 43);

    auto tr_y_xyz_z = pbuffer.data(idx_dip_fp + 44);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xy_x, tr_y_xyz_x, tr_y_xyz_y, tr_y_xyz_z, tr_y_yz_y, tr_y_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_y_xyz_x[i] = tr_y_xy_x[i] * pa_z[i];

        tr_y_xyz_y[i] = tr_y_yz_y[i] * pa_x[i];

        tr_y_xyz_z[i] = tr_y_yz_z[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : FP

    auto tr_y_xzz_x = pbuffer.data(idx_dip_fp + 45);

    auto tr_y_xzz_y = pbuffer.data(idx_dip_fp + 46);

    auto tr_y_xzz_z = pbuffer.data(idx_dip_fp + 47);

#pragma omp simd aligned(pa_x, tr_y_xzz_x, tr_y_xzz_y, tr_y_xzz_z, tr_y_zz_0, tr_y_zz_x, tr_y_zz_y, tr_y_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzz_x[i] = tr_y_zz_0[i] * fe_0 + tr_y_zz_x[i] * pa_x[i];

        tr_y_xzz_y[i] = tr_y_zz_y[i] * pa_x[i];

        tr_y_xzz_z[i] = tr_y_zz_z[i] * pa_x[i];
    }

    // Set up 48-51 components of targeted buffer : FP

    auto tr_y_yyy_x = pbuffer.data(idx_dip_fp + 48);

    auto tr_y_yyy_y = pbuffer.data(idx_dip_fp + 49);

    auto tr_y_yyy_z = pbuffer.data(idx_dip_fp + 50);

#pragma omp simd aligned(pa_y,           \
                             tr_y_y_x,   \
                             tr_y_y_y,   \
                             tr_y_y_z,   \
                             tr_y_yy_0,  \
                             tr_y_yy_x,  \
                             tr_y_yy_y,  \
                             tr_y_yy_z,  \
                             tr_y_yyy_x, \
                             tr_y_yyy_y, \
                             tr_y_yyy_z, \
                             ts_yy_x,    \
                             ts_yy_y,    \
                             ts_yy_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyy_x[i] = 2.0 * tr_y_y_x[i] * fe_0 + ts_yy_x[i] * fe_0 + tr_y_yy_x[i] * pa_y[i];

        tr_y_yyy_y[i] = 2.0 * tr_y_y_y[i] * fe_0 + tr_y_yy_0[i] * fe_0 + ts_yy_y[i] * fe_0 + tr_y_yy_y[i] * pa_y[i];

        tr_y_yyy_z[i] = 2.0 * tr_y_y_z[i] * fe_0 + ts_yy_z[i] * fe_0 + tr_y_yy_z[i] * pa_y[i];
    }

    // Set up 51-54 components of targeted buffer : FP

    auto tr_y_yyz_x = pbuffer.data(idx_dip_fp + 51);

    auto tr_y_yyz_y = pbuffer.data(idx_dip_fp + 52);

    auto tr_y_yyz_z = pbuffer.data(idx_dip_fp + 53);

#pragma omp simd aligned(pa_z, tr_y_yy_0, tr_y_yy_x, tr_y_yy_y, tr_y_yy_z, tr_y_yyz_x, tr_y_yyz_y, tr_y_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyz_x[i] = tr_y_yy_x[i] * pa_z[i];

        tr_y_yyz_y[i] = tr_y_yy_y[i] * pa_z[i];

        tr_y_yyz_z[i] = tr_y_yy_0[i] * fe_0 + tr_y_yy_z[i] * pa_z[i];
    }

    // Set up 54-57 components of targeted buffer : FP

    auto tr_y_yzz_x = pbuffer.data(idx_dip_fp + 54);

    auto tr_y_yzz_y = pbuffer.data(idx_dip_fp + 55);

    auto tr_y_yzz_z = pbuffer.data(idx_dip_fp + 56);

#pragma omp simd aligned(pa_y, pa_z, tr_y_y_y, tr_y_yz_y, tr_y_yzz_x, tr_y_yzz_y, tr_y_yzz_z, tr_y_zz_x, tr_y_zz_z, ts_zz_x, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzz_x[i] = ts_zz_x[i] * fe_0 + tr_y_zz_x[i] * pa_y[i];

        tr_y_yzz_y[i] = tr_y_y_y[i] * fe_0 + tr_y_yz_y[i] * pa_z[i];

        tr_y_yzz_z[i] = ts_zz_z[i] * fe_0 + tr_y_zz_z[i] * pa_y[i];
    }

    // Set up 57-60 components of targeted buffer : FP

    auto tr_y_zzz_x = pbuffer.data(idx_dip_fp + 57);

    auto tr_y_zzz_y = pbuffer.data(idx_dip_fp + 58);

    auto tr_y_zzz_z = pbuffer.data(idx_dip_fp + 59);

#pragma omp simd aligned( \
        pa_z, tr_y_z_x, tr_y_z_y, tr_y_z_z, tr_y_zz_0, tr_y_zz_x, tr_y_zz_y, tr_y_zz_z, tr_y_zzz_x, tr_y_zzz_y, tr_y_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzz_x[i] = 2.0 * tr_y_z_x[i] * fe_0 + tr_y_zz_x[i] * pa_z[i];

        tr_y_zzz_y[i] = 2.0 * tr_y_z_y[i] * fe_0 + tr_y_zz_y[i] * pa_z[i];

        tr_y_zzz_z[i] = 2.0 * tr_y_z_z[i] * fe_0 + tr_y_zz_0[i] * fe_0 + tr_y_zz_z[i] * pa_z[i];
    }

    // Set up 60-63 components of targeted buffer : FP

    auto tr_z_xxx_x = pbuffer.data(idx_dip_fp + 60);

    auto tr_z_xxx_y = pbuffer.data(idx_dip_fp + 61);

    auto tr_z_xxx_z = pbuffer.data(idx_dip_fp + 62);

#pragma omp simd aligned( \
        pa_x, tr_z_x_x, tr_z_x_y, tr_z_x_z, tr_z_xx_0, tr_z_xx_x, tr_z_xx_y, tr_z_xx_z, tr_z_xxx_x, tr_z_xxx_y, tr_z_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxx_x[i] = 2.0 * tr_z_x_x[i] * fe_0 + tr_z_xx_0[i] * fe_0 + tr_z_xx_x[i] * pa_x[i];

        tr_z_xxx_y[i] = 2.0 * tr_z_x_y[i] * fe_0 + tr_z_xx_y[i] * pa_x[i];

        tr_z_xxx_z[i] = 2.0 * tr_z_x_z[i] * fe_0 + tr_z_xx_z[i] * pa_x[i];
    }

    // Set up 63-66 components of targeted buffer : FP

    auto tr_z_xxy_x = pbuffer.data(idx_dip_fp + 63);

    auto tr_z_xxy_y = pbuffer.data(idx_dip_fp + 64);

    auto tr_z_xxy_z = pbuffer.data(idx_dip_fp + 65);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xx_x, tr_z_xx_z, tr_z_xxy_x, tr_z_xxy_y, tr_z_xxy_z, tr_z_xy_y, tr_z_y_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxy_x[i] = tr_z_xx_x[i] * pa_y[i];

        tr_z_xxy_y[i] = tr_z_y_y[i] * fe_0 + tr_z_xy_y[i] * pa_x[i];

        tr_z_xxy_z[i] = tr_z_xx_z[i] * pa_y[i];
    }

    // Set up 66-69 components of targeted buffer : FP

    auto tr_z_xxz_x = pbuffer.data(idx_dip_fp + 66);

    auto tr_z_xxz_y = pbuffer.data(idx_dip_fp + 67);

    auto tr_z_xxz_z = pbuffer.data(idx_dip_fp + 68);

#pragma omp simd aligned(pa_x, pa_z, tr_z_xx_x, tr_z_xxz_x, tr_z_xxz_y, tr_z_xxz_z, tr_z_xz_y, tr_z_xz_z, tr_z_z_y, tr_z_z_z, ts_xx_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxz_x[i] = ts_xx_x[i] * fe_0 + tr_z_xx_x[i] * pa_z[i];

        tr_z_xxz_y[i] = tr_z_z_y[i] * fe_0 + tr_z_xz_y[i] * pa_x[i];

        tr_z_xxz_z[i] = tr_z_z_z[i] * fe_0 + tr_z_xz_z[i] * pa_x[i];
    }

    // Set up 69-72 components of targeted buffer : FP

    auto tr_z_xyy_x = pbuffer.data(idx_dip_fp + 69);

    auto tr_z_xyy_y = pbuffer.data(idx_dip_fp + 70);

    auto tr_z_xyy_z = pbuffer.data(idx_dip_fp + 71);

#pragma omp simd aligned(pa_x, tr_z_xyy_x, tr_z_xyy_y, tr_z_xyy_z, tr_z_yy_0, tr_z_yy_x, tr_z_yy_y, tr_z_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyy_x[i] = tr_z_yy_0[i] * fe_0 + tr_z_yy_x[i] * pa_x[i];

        tr_z_xyy_y[i] = tr_z_yy_y[i] * pa_x[i];

        tr_z_xyy_z[i] = tr_z_yy_z[i] * pa_x[i];
    }

    // Set up 72-75 components of targeted buffer : FP

    auto tr_z_xyz_x = pbuffer.data(idx_dip_fp + 72);

    auto tr_z_xyz_y = pbuffer.data(idx_dip_fp + 73);

    auto tr_z_xyz_z = pbuffer.data(idx_dip_fp + 74);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xyz_x, tr_z_xyz_y, tr_z_xyz_z, tr_z_xz_x, tr_z_yz_y, tr_z_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_z_xyz_x[i] = tr_z_xz_x[i] * pa_y[i];

        tr_z_xyz_y[i] = tr_z_yz_y[i] * pa_x[i];

        tr_z_xyz_z[i] = tr_z_yz_z[i] * pa_x[i];
    }

    // Set up 75-78 components of targeted buffer : FP

    auto tr_z_xzz_x = pbuffer.data(idx_dip_fp + 75);

    auto tr_z_xzz_y = pbuffer.data(idx_dip_fp + 76);

    auto tr_z_xzz_z = pbuffer.data(idx_dip_fp + 77);

#pragma omp simd aligned(pa_x, tr_z_xzz_x, tr_z_xzz_y, tr_z_xzz_z, tr_z_zz_0, tr_z_zz_x, tr_z_zz_y, tr_z_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzz_x[i] = tr_z_zz_0[i] * fe_0 + tr_z_zz_x[i] * pa_x[i];

        tr_z_xzz_y[i] = tr_z_zz_y[i] * pa_x[i];

        tr_z_xzz_z[i] = tr_z_zz_z[i] * pa_x[i];
    }

    // Set up 78-81 components of targeted buffer : FP

    auto tr_z_yyy_x = pbuffer.data(idx_dip_fp + 78);

    auto tr_z_yyy_y = pbuffer.data(idx_dip_fp + 79);

    auto tr_z_yyy_z = pbuffer.data(idx_dip_fp + 80);

#pragma omp simd aligned( \
        pa_y, tr_z_y_x, tr_z_y_y, tr_z_y_z, tr_z_yy_0, tr_z_yy_x, tr_z_yy_y, tr_z_yy_z, tr_z_yyy_x, tr_z_yyy_y, tr_z_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyy_x[i] = 2.0 * tr_z_y_x[i] * fe_0 + tr_z_yy_x[i] * pa_y[i];

        tr_z_yyy_y[i] = 2.0 * tr_z_y_y[i] * fe_0 + tr_z_yy_0[i] * fe_0 + tr_z_yy_y[i] * pa_y[i];

        tr_z_yyy_z[i] = 2.0 * tr_z_y_z[i] * fe_0 + tr_z_yy_z[i] * pa_y[i];
    }

    // Set up 81-84 components of targeted buffer : FP

    auto tr_z_yyz_x = pbuffer.data(idx_dip_fp + 81);

    auto tr_z_yyz_y = pbuffer.data(idx_dip_fp + 82);

    auto tr_z_yyz_z = pbuffer.data(idx_dip_fp + 83);

#pragma omp simd aligned(pa_y, pa_z, tr_z_yy_y, tr_z_yyz_x, tr_z_yyz_y, tr_z_yyz_z, tr_z_yz_x, tr_z_yz_z, tr_z_z_x, tr_z_z_z, ts_yy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyz_x[i] = tr_z_z_x[i] * fe_0 + tr_z_yz_x[i] * pa_y[i];

        tr_z_yyz_y[i] = ts_yy_y[i] * fe_0 + tr_z_yy_y[i] * pa_z[i];

        tr_z_yyz_z[i] = tr_z_z_z[i] * fe_0 + tr_z_yz_z[i] * pa_y[i];
    }

    // Set up 84-87 components of targeted buffer : FP

    auto tr_z_yzz_x = pbuffer.data(idx_dip_fp + 84);

    auto tr_z_yzz_y = pbuffer.data(idx_dip_fp + 85);

    auto tr_z_yzz_z = pbuffer.data(idx_dip_fp + 86);

#pragma omp simd aligned(pa_y, tr_z_yzz_x, tr_z_yzz_y, tr_z_yzz_z, tr_z_zz_0, tr_z_zz_x, tr_z_zz_y, tr_z_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzz_x[i] = tr_z_zz_x[i] * pa_y[i];

        tr_z_yzz_y[i] = tr_z_zz_0[i] * fe_0 + tr_z_zz_y[i] * pa_y[i];

        tr_z_yzz_z[i] = tr_z_zz_z[i] * pa_y[i];
    }

    // Set up 87-90 components of targeted buffer : FP

    auto tr_z_zzz_x = pbuffer.data(idx_dip_fp + 87);

    auto tr_z_zzz_y = pbuffer.data(idx_dip_fp + 88);

    auto tr_z_zzz_z = pbuffer.data(idx_dip_fp + 89);

#pragma omp simd aligned(pa_z,           \
                             tr_z_z_x,   \
                             tr_z_z_y,   \
                             tr_z_z_z,   \
                             tr_z_zz_0,  \
                             tr_z_zz_x,  \
                             tr_z_zz_y,  \
                             tr_z_zz_z,  \
                             tr_z_zzz_x, \
                             tr_z_zzz_y, \
                             tr_z_zzz_z, \
                             ts_zz_x,    \
                             ts_zz_y,    \
                             ts_zz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzz_x[i] = 2.0 * tr_z_z_x[i] * fe_0 + ts_zz_x[i] * fe_0 + tr_z_zz_x[i] * pa_z[i];

        tr_z_zzz_y[i] = 2.0 * tr_z_z_y[i] * fe_0 + ts_zz_y[i] * fe_0 + tr_z_zz_y[i] * pa_z[i];

        tr_z_zzz_z[i] = 2.0 * tr_z_z_z[i] * fe_0 + tr_z_zz_0[i] * fe_0 + ts_zz_z[i] * fe_0 + tr_z_zz_z[i] * pa_z[i];
    }
}

}  // namespace diprec
