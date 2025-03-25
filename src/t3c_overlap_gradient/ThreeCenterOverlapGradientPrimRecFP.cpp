#include "ThreeCenterOverlapGradientPrimRecFP.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_fp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_fp,
                              const size_t idx_dp,
                              const size_t idx_fs,
                              const size_t idx_fp,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    auto ts_xy_x = pbuffer.data(idx_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_dp + 5);

    auto ts_xz_x = pbuffer.data(idx_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    auto ts_yz_x = pbuffer.data(idx_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_fs);

    auto ts_xxy_0 = pbuffer.data(idx_fs + 1);

    auto ts_xxz_0 = pbuffer.data(idx_fs + 2);

    auto ts_xyy_0 = pbuffer.data(idx_fs + 3);

    auto ts_xyz_0 = pbuffer.data(idx_fs + 4);

    auto ts_xzz_0 = pbuffer.data(idx_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_fs + 6);

    auto ts_yyz_0 = pbuffer.data(idx_fs + 7);

    auto ts_yzz_0 = pbuffer.data(idx_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_fp);

    auto ts_xxx_y = pbuffer.data(idx_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_fp + 4);

    auto ts_xxy_z = pbuffer.data(idx_fp + 5);

    auto ts_xxz_x = pbuffer.data(idx_fp + 6);

    auto ts_xxz_y = pbuffer.data(idx_fp + 7);

    auto ts_xxz_z = pbuffer.data(idx_fp + 8);

    auto ts_xyy_x = pbuffer.data(idx_fp + 9);

    auto ts_xyy_y = pbuffer.data(idx_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_fp + 11);

    auto ts_xyz_x = pbuffer.data(idx_fp + 12);

    auto ts_xyz_y = pbuffer.data(idx_fp + 13);

    auto ts_xyz_z = pbuffer.data(idx_fp + 14);

    auto ts_xzz_x = pbuffer.data(idx_fp + 15);

    auto ts_xzz_y = pbuffer.data(idx_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_fp + 20);

    auto ts_yyz_x = pbuffer.data(idx_fp + 21);

    auto ts_yyz_y = pbuffer.data(idx_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_fp + 23);

    auto ts_yzz_x = pbuffer.data(idx_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_fp + 29);

    // Set up 0-3 components of targeted buffer : FP

    auto gs_x_xxx_x = pbuffer.data(idx_g_fp);

    auto gs_x_xxx_y = pbuffer.data(idx_g_fp + 1);

    auto gs_x_xxx_z = pbuffer.data(idx_g_fp + 2);

    #pragma omp simd aligned(gc_x, gs_x_xxx_x, gs_x_xxx_y, gs_x_xxx_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxx_x[i] = 6.0 * ts_xx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_x[i] * gc_x[i] * tce_0;

        gs_x_xxx_y[i] = 6.0 * ts_xx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_y[i] * gc_x[i] * tce_0;

        gs_x_xxx_z[i] = 6.0 * ts_xx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_z[i] * gc_x[i] * tce_0;
    }

    // Set up 3-6 components of targeted buffer : FP

    auto gs_x_xxy_x = pbuffer.data(idx_g_fp + 3);

    auto gs_x_xxy_y = pbuffer.data(idx_g_fp + 4);

    auto gs_x_xxy_z = pbuffer.data(idx_g_fp + 5);

    #pragma omp simd aligned(gc_x, gs_x_xxy_x, gs_x_xxy_y, gs_x_xxy_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxy_x[i] = 4.0 * ts_xy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_x[i] * gc_x[i] * tce_0;

        gs_x_xxy_y[i] = 4.0 * ts_xy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_y[i] * gc_x[i] * tce_0;

        gs_x_xxy_z[i] = 4.0 * ts_xy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 6-9 components of targeted buffer : FP

    auto gs_x_xxz_x = pbuffer.data(idx_g_fp + 6);

    auto gs_x_xxz_y = pbuffer.data(idx_g_fp + 7);

    auto gs_x_xxz_z = pbuffer.data(idx_g_fp + 8);

    #pragma omp simd aligned(gc_x, gs_x_xxz_x, gs_x_xxz_y, gs_x_xxz_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxz_x[i] = 4.0 * ts_xz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_x[i] * gc_x[i] * tce_0;

        gs_x_xxz_y[i] = 4.0 * ts_xz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_y[i] * gc_x[i] * tce_0;

        gs_x_xxz_z[i] = 4.0 * ts_xz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 9-12 components of targeted buffer : FP

    auto gs_x_xyy_x = pbuffer.data(idx_g_fp + 9);

    auto gs_x_xyy_y = pbuffer.data(idx_g_fp + 10);

    auto gs_x_xyy_z = pbuffer.data(idx_g_fp + 11);

    #pragma omp simd aligned(gc_x, gs_x_xyy_x, gs_x_xyy_y, gs_x_xyy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyy_x[i] = 2.0 * ts_yy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_x[i] * gc_x[i] * tce_0;

        gs_x_xyy_y[i] = 2.0 * ts_yy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_y[i] * gc_x[i] * tce_0;

        gs_x_xyy_z[i] = 2.0 * ts_yy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 12-15 components of targeted buffer : FP

    auto gs_x_xyz_x = pbuffer.data(idx_g_fp + 12);

    auto gs_x_xyz_y = pbuffer.data(idx_g_fp + 13);

    auto gs_x_xyz_z = pbuffer.data(idx_g_fp + 14);

    #pragma omp simd aligned(gc_x, gs_x_xyz_x, gs_x_xyz_y, gs_x_xyz_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyz_x[i] = 2.0 * ts_yz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_x[i] * gc_x[i] * tce_0;

        gs_x_xyz_y[i] = 2.0 * ts_yz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_y[i] * gc_x[i] * tce_0;

        gs_x_xyz_z[i] = 2.0 * ts_yz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 15-18 components of targeted buffer : FP

    auto gs_x_xzz_x = pbuffer.data(idx_g_fp + 15);

    auto gs_x_xzz_y = pbuffer.data(idx_g_fp + 16);

    auto gs_x_xzz_z = pbuffer.data(idx_g_fp + 17);

    #pragma omp simd aligned(gc_x, gs_x_xzz_x, gs_x_xzz_y, gs_x_xzz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzz_x[i] = 2.0 * ts_zz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_x[i] * gc_x[i] * tce_0;

        gs_x_xzz_y[i] = 2.0 * ts_zz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_y[i] * gc_x[i] * tce_0;

        gs_x_xzz_z[i] = 2.0 * ts_zz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 18-21 components of targeted buffer : FP

    auto gs_x_yyy_x = pbuffer.data(idx_g_fp + 18);

    auto gs_x_yyy_y = pbuffer.data(idx_g_fp + 19);

    auto gs_x_yyy_z = pbuffer.data(idx_g_fp + 20);

    #pragma omp simd aligned(gc_x, gs_x_yyy_x, gs_x_yyy_y, gs_x_yyy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyy_x[i] = 2.0 * ts_yyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_x[i] * gc_x[i] * tce_0;

        gs_x_yyy_y[i] = 2.0 * ts_yyy_y[i] * gc_x[i] * tce_0;

        gs_x_yyy_z[i] = 2.0 * ts_yyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 21-24 components of targeted buffer : FP

    auto gs_x_yyz_x = pbuffer.data(idx_g_fp + 21);

    auto gs_x_yyz_y = pbuffer.data(idx_g_fp + 22);

    auto gs_x_yyz_z = pbuffer.data(idx_g_fp + 23);

    #pragma omp simd aligned(gc_x, gs_x_yyz_x, gs_x_yyz_y, gs_x_yyz_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyz_x[i] = 2.0 * ts_yyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_x[i] * gc_x[i] * tce_0;

        gs_x_yyz_y[i] = 2.0 * ts_yyz_y[i] * gc_x[i] * tce_0;

        gs_x_yyz_z[i] = 2.0 * ts_yyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 24-27 components of targeted buffer : FP

    auto gs_x_yzz_x = pbuffer.data(idx_g_fp + 24);

    auto gs_x_yzz_y = pbuffer.data(idx_g_fp + 25);

    auto gs_x_yzz_z = pbuffer.data(idx_g_fp + 26);

    #pragma omp simd aligned(gc_x, gs_x_yzz_x, gs_x_yzz_y, gs_x_yzz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzz_x[i] = 2.0 * ts_yzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_x[i] * gc_x[i] * tce_0;

        gs_x_yzz_y[i] = 2.0 * ts_yzz_y[i] * gc_x[i] * tce_0;

        gs_x_yzz_z[i] = 2.0 * ts_yzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 27-30 components of targeted buffer : FP

    auto gs_x_zzz_x = pbuffer.data(idx_g_fp + 27);

    auto gs_x_zzz_y = pbuffer.data(idx_g_fp + 28);

    auto gs_x_zzz_z = pbuffer.data(idx_g_fp + 29);

    #pragma omp simd aligned(gc_x, gs_x_zzz_x, gs_x_zzz_y, gs_x_zzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzz_x[i] = 2.0 * ts_zzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_x[i] * gc_x[i] * tce_0;

        gs_x_zzz_y[i] = 2.0 * ts_zzz_y[i] * gc_x[i] * tce_0;

        gs_x_zzz_z[i] = 2.0 * ts_zzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 30-33 components of targeted buffer : FP

    auto gs_y_xxx_x = pbuffer.data(idx_g_fp + 30);

    auto gs_y_xxx_y = pbuffer.data(idx_g_fp + 31);

    auto gs_y_xxx_z = pbuffer.data(idx_g_fp + 32);

    #pragma omp simd aligned(gc_y, gs_y_xxx_x, gs_y_xxx_y, gs_y_xxx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxx_x[i] = 2.0 * ts_xxx_x[i] * gc_y[i] * tce_0;

        gs_y_xxx_y[i] = 2.0 * ts_xxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_y[i] * gc_y[i] * tce_0;

        gs_y_xxx_z[i] = 2.0 * ts_xxx_z[i] * gc_y[i] * tce_0;
    }

    // Set up 33-36 components of targeted buffer : FP

    auto gs_y_xxy_x = pbuffer.data(idx_g_fp + 33);

    auto gs_y_xxy_y = pbuffer.data(idx_g_fp + 34);

    auto gs_y_xxy_z = pbuffer.data(idx_g_fp + 35);

    #pragma omp simd aligned(gc_y, gs_y_xxy_x, gs_y_xxy_y, gs_y_xxy_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxy_x[i] = 2.0 * ts_xx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_x[i] * gc_y[i] * tce_0;

        gs_y_xxy_y[i] = 2.0 * ts_xx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_y[i] * gc_y[i] * tce_0;

        gs_y_xxy_z[i] = 2.0 * ts_xx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 36-39 components of targeted buffer : FP

    auto gs_y_xxz_x = pbuffer.data(idx_g_fp + 36);

    auto gs_y_xxz_y = pbuffer.data(idx_g_fp + 37);

    auto gs_y_xxz_z = pbuffer.data(idx_g_fp + 38);

    #pragma omp simd aligned(gc_y, gs_y_xxz_x, gs_y_xxz_y, gs_y_xxz_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxz_x[i] = 2.0 * ts_xxz_x[i] * gc_y[i] * tce_0;

        gs_y_xxz_y[i] = 2.0 * ts_xxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_y[i] * gc_y[i] * tce_0;

        gs_y_xxz_z[i] = 2.0 * ts_xxz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 39-42 components of targeted buffer : FP

    auto gs_y_xyy_x = pbuffer.data(idx_g_fp + 39);

    auto gs_y_xyy_y = pbuffer.data(idx_g_fp + 40);

    auto gs_y_xyy_z = pbuffer.data(idx_g_fp + 41);

    #pragma omp simd aligned(gc_y, gs_y_xyy_x, gs_y_xyy_y, gs_y_xyy_z, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyy_x[i] = 4.0 * ts_xy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_x[i] * gc_y[i] * tce_0;

        gs_y_xyy_y[i] = 4.0 * ts_xy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_y[i] * gc_y[i] * tce_0;

        gs_y_xyy_z[i] = 4.0 * ts_xy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 42-45 components of targeted buffer : FP

    auto gs_y_xyz_x = pbuffer.data(idx_g_fp + 42);

    auto gs_y_xyz_y = pbuffer.data(idx_g_fp + 43);

    auto gs_y_xyz_z = pbuffer.data(idx_g_fp + 44);

    #pragma omp simd aligned(gc_y, gs_y_xyz_x, gs_y_xyz_y, gs_y_xyz_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyz_x[i] = 2.0 * ts_xz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_x[i] * gc_y[i] * tce_0;

        gs_y_xyz_y[i] = 2.0 * ts_xz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_y[i] * gc_y[i] * tce_0;

        gs_y_xyz_z[i] = 2.0 * ts_xz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 45-48 components of targeted buffer : FP

    auto gs_y_xzz_x = pbuffer.data(idx_g_fp + 45);

    auto gs_y_xzz_y = pbuffer.data(idx_g_fp + 46);

    auto gs_y_xzz_z = pbuffer.data(idx_g_fp + 47);

    #pragma omp simd aligned(gc_y, gs_y_xzz_x, gs_y_xzz_y, gs_y_xzz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzz_x[i] = 2.0 * ts_xzz_x[i] * gc_y[i] * tce_0;

        gs_y_xzz_y[i] = 2.0 * ts_xzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_y[i] * gc_y[i] * tce_0;

        gs_y_xzz_z[i] = 2.0 * ts_xzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 48-51 components of targeted buffer : FP

    auto gs_y_yyy_x = pbuffer.data(idx_g_fp + 48);

    auto gs_y_yyy_y = pbuffer.data(idx_g_fp + 49);

    auto gs_y_yyy_z = pbuffer.data(idx_g_fp + 50);

    #pragma omp simd aligned(gc_y, gs_y_yyy_x, gs_y_yyy_y, gs_y_yyy_z, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyy_x[i] = 6.0 * ts_yy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_x[i] * gc_y[i] * tce_0;

        gs_y_yyy_y[i] = 6.0 * ts_yy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_y[i] * gc_y[i] * tce_0;

        gs_y_yyy_z[i] = 6.0 * ts_yy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 51-54 components of targeted buffer : FP

    auto gs_y_yyz_x = pbuffer.data(idx_g_fp + 51);

    auto gs_y_yyz_y = pbuffer.data(idx_g_fp + 52);

    auto gs_y_yyz_z = pbuffer.data(idx_g_fp + 53);

    #pragma omp simd aligned(gc_y, gs_y_yyz_x, gs_y_yyz_y, gs_y_yyz_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyz_x[i] = 4.0 * ts_yz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_x[i] * gc_y[i] * tce_0;

        gs_y_yyz_y[i] = 4.0 * ts_yz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_y[i] * gc_y[i] * tce_0;

        gs_y_yyz_z[i] = 4.0 * ts_yz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 54-57 components of targeted buffer : FP

    auto gs_y_yzz_x = pbuffer.data(idx_g_fp + 54);

    auto gs_y_yzz_y = pbuffer.data(idx_g_fp + 55);

    auto gs_y_yzz_z = pbuffer.data(idx_g_fp + 56);

    #pragma omp simd aligned(gc_y, gs_y_yzz_x, gs_y_yzz_y, gs_y_yzz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzz_x[i] = 2.0 * ts_zz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_x[i] * gc_y[i] * tce_0;

        gs_y_yzz_y[i] = 2.0 * ts_zz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_y[i] * gc_y[i] * tce_0;

        gs_y_yzz_z[i] = 2.0 * ts_zz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 57-60 components of targeted buffer : FP

    auto gs_y_zzz_x = pbuffer.data(idx_g_fp + 57);

    auto gs_y_zzz_y = pbuffer.data(idx_g_fp + 58);

    auto gs_y_zzz_z = pbuffer.data(idx_g_fp + 59);

    #pragma omp simd aligned(gc_y, gs_y_zzz_x, gs_y_zzz_y, gs_y_zzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzz_x[i] = 2.0 * ts_zzz_x[i] * gc_y[i] * tce_0;

        gs_y_zzz_y[i] = 2.0 * ts_zzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_y[i] * gc_y[i] * tce_0;

        gs_y_zzz_z[i] = 2.0 * ts_zzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 60-63 components of targeted buffer : FP

    auto gs_z_xxx_x = pbuffer.data(idx_g_fp + 60);

    auto gs_z_xxx_y = pbuffer.data(idx_g_fp + 61);

    auto gs_z_xxx_z = pbuffer.data(idx_g_fp + 62);

    #pragma omp simd aligned(gc_z, gs_z_xxx_x, gs_z_xxx_y, gs_z_xxx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxx_x[i] = 2.0 * ts_xxx_x[i] * gc_z[i] * tce_0;

        gs_z_xxx_y[i] = 2.0 * ts_xxx_y[i] * gc_z[i] * tce_0;

        gs_z_xxx_z[i] = 2.0 * ts_xxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_z[i] * gc_z[i] * tce_0;
    }

    // Set up 63-66 components of targeted buffer : FP

    auto gs_z_xxy_x = pbuffer.data(idx_g_fp + 63);

    auto gs_z_xxy_y = pbuffer.data(idx_g_fp + 64);

    auto gs_z_xxy_z = pbuffer.data(idx_g_fp + 65);

    #pragma omp simd aligned(gc_z, gs_z_xxy_x, gs_z_xxy_y, gs_z_xxy_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxy_x[i] = 2.0 * ts_xxy_x[i] * gc_z[i] * tce_0;

        gs_z_xxy_y[i] = 2.0 * ts_xxy_y[i] * gc_z[i] * tce_0;

        gs_z_xxy_z[i] = 2.0 * ts_xxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 66-69 components of targeted buffer : FP

    auto gs_z_xxz_x = pbuffer.data(idx_g_fp + 66);

    auto gs_z_xxz_y = pbuffer.data(idx_g_fp + 67);

    auto gs_z_xxz_z = pbuffer.data(idx_g_fp + 68);

    #pragma omp simd aligned(gc_z, gs_z_xxz_x, gs_z_xxz_y, gs_z_xxz_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxz_x[i] = 2.0 * ts_xx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_x[i] * gc_z[i] * tce_0;

        gs_z_xxz_y[i] = 2.0 * ts_xx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_y[i] * gc_z[i] * tce_0;

        gs_z_xxz_z[i] = 2.0 * ts_xx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 69-72 components of targeted buffer : FP

    auto gs_z_xyy_x = pbuffer.data(idx_g_fp + 69);

    auto gs_z_xyy_y = pbuffer.data(idx_g_fp + 70);

    auto gs_z_xyy_z = pbuffer.data(idx_g_fp + 71);

    #pragma omp simd aligned(gc_z, gs_z_xyy_x, gs_z_xyy_y, gs_z_xyy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyy_x[i] = 2.0 * ts_xyy_x[i] * gc_z[i] * tce_0;

        gs_z_xyy_y[i] = 2.0 * ts_xyy_y[i] * gc_z[i] * tce_0;

        gs_z_xyy_z[i] = 2.0 * ts_xyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 72-75 components of targeted buffer : FP

    auto gs_z_xyz_x = pbuffer.data(idx_g_fp + 72);

    auto gs_z_xyz_y = pbuffer.data(idx_g_fp + 73);

    auto gs_z_xyz_z = pbuffer.data(idx_g_fp + 74);

    #pragma omp simd aligned(gc_z, gs_z_xyz_x, gs_z_xyz_y, gs_z_xyz_z, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyz_x[i] = 2.0 * ts_xy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_x[i] * gc_z[i] * tce_0;

        gs_z_xyz_y[i] = 2.0 * ts_xy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_y[i] * gc_z[i] * tce_0;

        gs_z_xyz_z[i] = 2.0 * ts_xy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 75-78 components of targeted buffer : FP

    auto gs_z_xzz_x = pbuffer.data(idx_g_fp + 75);

    auto gs_z_xzz_y = pbuffer.data(idx_g_fp + 76);

    auto gs_z_xzz_z = pbuffer.data(idx_g_fp + 77);

    #pragma omp simd aligned(gc_z, gs_z_xzz_x, gs_z_xzz_y, gs_z_xzz_z, ts_xz_x, ts_xz_y, ts_xz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzz_x[i] = 4.0 * ts_xz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_x[i] * gc_z[i] * tce_0;

        gs_z_xzz_y[i] = 4.0 * ts_xz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_y[i] * gc_z[i] * tce_0;

        gs_z_xzz_z[i] = 4.0 * ts_xz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 78-81 components of targeted buffer : FP

    auto gs_z_yyy_x = pbuffer.data(idx_g_fp + 78);

    auto gs_z_yyy_y = pbuffer.data(idx_g_fp + 79);

    auto gs_z_yyy_z = pbuffer.data(idx_g_fp + 80);

    #pragma omp simd aligned(gc_z, gs_z_yyy_x, gs_z_yyy_y, gs_z_yyy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyy_x[i] = 2.0 * ts_yyy_x[i] * gc_z[i] * tce_0;

        gs_z_yyy_y[i] = 2.0 * ts_yyy_y[i] * gc_z[i] * tce_0;

        gs_z_yyy_z[i] = 2.0 * ts_yyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 81-84 components of targeted buffer : FP

    auto gs_z_yyz_x = pbuffer.data(idx_g_fp + 81);

    auto gs_z_yyz_y = pbuffer.data(idx_g_fp + 82);

    auto gs_z_yyz_z = pbuffer.data(idx_g_fp + 83);

    #pragma omp simd aligned(gc_z, gs_z_yyz_x, gs_z_yyz_y, gs_z_yyz_z, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyz_x[i] = 2.0 * ts_yy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_x[i] * gc_z[i] * tce_0;

        gs_z_yyz_y[i] = 2.0 * ts_yy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_y[i] * gc_z[i] * tce_0;

        gs_z_yyz_z[i] = 2.0 * ts_yy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 84-87 components of targeted buffer : FP

    auto gs_z_yzz_x = pbuffer.data(idx_g_fp + 84);

    auto gs_z_yzz_y = pbuffer.data(idx_g_fp + 85);

    auto gs_z_yzz_z = pbuffer.data(idx_g_fp + 86);

    #pragma omp simd aligned(gc_z, gs_z_yzz_x, gs_z_yzz_y, gs_z_yzz_z, ts_yz_x, ts_yz_y, ts_yz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzz_x[i] = 4.0 * ts_yz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_x[i] * gc_z[i] * tce_0;

        gs_z_yzz_y[i] = 4.0 * ts_yz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_y[i] * gc_z[i] * tce_0;

        gs_z_yzz_z[i] = 4.0 * ts_yz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 87-90 components of targeted buffer : FP

    auto gs_z_zzz_x = pbuffer.data(idx_g_fp + 87);

    auto gs_z_zzz_y = pbuffer.data(idx_g_fp + 88);

    auto gs_z_zzz_z = pbuffer.data(idx_g_fp + 89);

    #pragma omp simd aligned(gc_z, gs_z_zzz_x, gs_z_zzz_y, gs_z_zzz_z, ts_zz_x, ts_zz_y, ts_zz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzz_x[i] = 6.0 * ts_zz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_x[i] * gc_z[i] * tce_0;

        gs_z_zzz_y[i] = 6.0 * ts_zz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_y[i] * gc_z[i] * tce_0;

        gs_z_zzz_z[i] = 6.0 * ts_zz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_z[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

