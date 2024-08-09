#include "OverlapPrimRecFP.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_fp(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_fp,
                     const size_t              idx_ovl_pp,
                     const size_t              idx_ovl_ds,
                     const size_t              idx_ovl_dp,
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

    auto ts_x_x = pbuffer.data(idx_ovl_pp);

    auto ts_x_y = pbuffer.data(idx_ovl_pp + 1);

    auto ts_x_z = pbuffer.data(idx_ovl_pp + 2);

    auto ts_y_x = pbuffer.data(idx_ovl_pp + 3);

    auto ts_y_y = pbuffer.data(idx_ovl_pp + 4);

    auto ts_y_z = pbuffer.data(idx_ovl_pp + 5);

    auto ts_z_x = pbuffer.data(idx_ovl_pp + 6);

    auto ts_z_y = pbuffer.data(idx_ovl_pp + 7);

    auto ts_z_z = pbuffer.data(idx_ovl_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ovl_ds);

    auto ts_yy_0 = pbuffer.data(idx_ovl_ds + 3);

    auto ts_zz_0 = pbuffer.data(idx_ovl_ds + 5);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_ovl_dp);

    auto ts_xx_y = pbuffer.data(idx_ovl_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_ovl_dp + 2);

    auto ts_xy_y = pbuffer.data(idx_ovl_dp + 4);

    auto ts_xz_x = pbuffer.data(idx_ovl_dp + 6);

    auto ts_xz_z = pbuffer.data(idx_ovl_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_ovl_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_ovl_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_ovl_dp + 11);

    auto ts_yz_y = pbuffer.data(idx_ovl_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_ovl_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_ovl_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_ovl_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_ovl_dp + 17);

    // Set up 0-3 components of targeted buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_ovl_fp);

    auto ts_xxx_y = pbuffer.data(idx_ovl_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_ovl_fp + 2);

#pragma omp simd aligned( \
        pa_x, ts_x_x, ts_x_y, ts_x_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxx_x[i] = 2.0 * ts_x_x[i] * fe_0 + ts_xx_0[i] * fe_0 + ts_xx_x[i] * pa_x[i];

        ts_xxx_y[i] = 2.0 * ts_x_y[i] * fe_0 + ts_xx_y[i] * pa_x[i];

        ts_xxx_z[i] = 2.0 * ts_x_z[i] * fe_0 + ts_xx_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto ts_xxy_x = pbuffer.data(idx_ovl_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_ovl_fp + 4);

    auto ts_xxy_z = pbuffer.data(idx_ovl_fp + 5);

#pragma omp simd aligned(pa_x, pa_y, ts_xx_x, ts_xx_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xy_y, ts_y_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxy_x[i] = ts_xx_x[i] * pa_y[i];

        ts_xxy_y[i] = ts_y_y[i] * fe_0 + ts_xy_y[i] * pa_x[i];

        ts_xxy_z[i] = ts_xx_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto ts_xxz_x = pbuffer.data(idx_ovl_fp + 6);

    auto ts_xxz_y = pbuffer.data(idx_ovl_fp + 7);

    auto ts_xxz_z = pbuffer.data(idx_ovl_fp + 8);

#pragma omp simd aligned(pa_x, pa_z, ts_xx_x, ts_xx_y, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xz_z, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxz_x[i] = ts_xx_x[i] * pa_z[i];

        ts_xxz_y[i] = ts_xx_y[i] * pa_z[i];

        ts_xxz_z[i] = ts_z_z[i] * fe_0 + ts_xz_z[i] * pa_x[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto ts_xyy_x = pbuffer.data(idx_ovl_fp + 9);

    auto ts_xyy_y = pbuffer.data(idx_ovl_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_ovl_fp + 11);

#pragma omp simd aligned(pa_x, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyy_x[i] = ts_yy_0[i] * fe_0 + ts_yy_x[i] * pa_x[i];

        ts_xyy_y[i] = ts_yy_y[i] * pa_x[i];

        ts_xyy_z[i] = ts_yy_z[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto ts_xyz_x = pbuffer.data(idx_ovl_fp + 12);

    auto ts_xyz_y = pbuffer.data(idx_ovl_fp + 13);

    auto ts_xyz_z = pbuffer.data(idx_ovl_fp + 14);

#pragma omp simd aligned(pa_x, pa_y, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyz_x[i] = ts_xz_x[i] * pa_y[i];

        ts_xyz_y[i] = ts_yz_y[i] * pa_x[i];

        ts_xyz_z[i] = ts_yz_z[i] * pa_x[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto ts_xzz_x = pbuffer.data(idx_ovl_fp + 15);

    auto ts_xzz_y = pbuffer.data(idx_ovl_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_ovl_fp + 17);

#pragma omp simd aligned(pa_x, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzz_x[i] = ts_zz_0[i] * fe_0 + ts_zz_x[i] * pa_x[i];

        ts_xzz_y[i] = ts_zz_y[i] * pa_x[i];

        ts_xzz_z[i] = ts_zz_z[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto ts_yyy_x = pbuffer.data(idx_ovl_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_ovl_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_ovl_fp + 20);

#pragma omp simd aligned( \
        pa_y, ts_y_x, ts_y_y, ts_y_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyy_x[i] = 2.0 * ts_y_x[i] * fe_0 + ts_yy_x[i] * pa_y[i];

        ts_yyy_y[i] = 2.0 * ts_y_y[i] * fe_0 + ts_yy_0[i] * fe_0 + ts_yy_y[i] * pa_y[i];

        ts_yyy_z[i] = 2.0 * ts_y_z[i] * fe_0 + ts_yy_z[i] * pa_y[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto ts_yyz_x = pbuffer.data(idx_ovl_fp + 21);

    auto ts_yyz_y = pbuffer.data(idx_ovl_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_ovl_fp + 23);

#pragma omp simd aligned(pa_y, pa_z, ts_yy_x, ts_yy_y, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yz_z, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyz_x[i] = ts_yy_x[i] * pa_z[i];

        ts_yyz_y[i] = ts_yy_y[i] * pa_z[i];

        ts_yyz_z[i] = ts_z_z[i] * fe_0 + ts_yz_z[i] * pa_y[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto ts_yzz_x = pbuffer.data(idx_ovl_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_ovl_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_ovl_fp + 26);

#pragma omp simd aligned(pa_y, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzz_x[i] = ts_zz_x[i] * pa_y[i];

        ts_yzz_y[i] = ts_zz_0[i] * fe_0 + ts_zz_y[i] * pa_y[i];

        ts_yzz_z[i] = ts_zz_z[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto ts_zzz_x = pbuffer.data(idx_ovl_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_ovl_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_ovl_fp + 29);

#pragma omp simd aligned( \
        pa_z, ts_z_x, ts_z_y, ts_z_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzz_x[i] = 2.0 * ts_z_x[i] * fe_0 + ts_zz_x[i] * pa_z[i];

        ts_zzz_y[i] = 2.0 * ts_z_y[i] * fe_0 + ts_zz_y[i] * pa_z[i];

        ts_zzz_z[i] = 2.0 * ts_z_z[i] * fe_0 + ts_zz_0[i] * fe_0 + ts_zz_z[i] * pa_z[i];
    }
}

}  // namespace ovlrec
