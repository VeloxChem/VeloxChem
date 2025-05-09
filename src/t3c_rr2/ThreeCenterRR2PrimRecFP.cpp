#include "ThreeCenterRR2PrimRecFP.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_fp(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_fp,
                  const size_t idx_dp,
                  const size_t idx_g_dp,
                  const size_t idx_fs,
                  const size_t idx_g_fs,
                  const size_t idx_fp,
                  const size_t idx_g_fp,
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

    // Set up components of auxiliary buffer : DP

    auto gr_xx_x = pbuffer.data(idx_g_dp);

    auto gr_xx_y = pbuffer.data(idx_g_dp + 1);

    auto gr_xx_z = pbuffer.data(idx_g_dp + 2);

    auto gr_xy_x = pbuffer.data(idx_g_dp + 3);

    auto gr_xy_y = pbuffer.data(idx_g_dp + 4);

    auto gr_xy_z = pbuffer.data(idx_g_dp + 5);

    auto gr_xz_x = pbuffer.data(idx_g_dp + 6);

    auto gr_xz_y = pbuffer.data(idx_g_dp + 7);

    auto gr_xz_z = pbuffer.data(idx_g_dp + 8);

    auto gr_yy_x = pbuffer.data(idx_g_dp + 9);

    auto gr_yy_y = pbuffer.data(idx_g_dp + 10);

    auto gr_yy_z = pbuffer.data(idx_g_dp + 11);

    auto gr_yz_x = pbuffer.data(idx_g_dp + 12);

    auto gr_yz_y = pbuffer.data(idx_g_dp + 13);

    auto gr_yz_z = pbuffer.data(idx_g_dp + 14);

    auto gr_zz_x = pbuffer.data(idx_g_dp + 15);

    auto gr_zz_y = pbuffer.data(idx_g_dp + 16);

    auto gr_zz_z = pbuffer.data(idx_g_dp + 17);

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

    // Set up components of auxiliary buffer : FS

    auto gr_xxx_0 = pbuffer.data(idx_g_fs);

    auto gr_xxy_0 = pbuffer.data(idx_g_fs + 1);

    auto gr_xxz_0 = pbuffer.data(idx_g_fs + 2);

    auto gr_xyy_0 = pbuffer.data(idx_g_fs + 3);

    auto gr_xyz_0 = pbuffer.data(idx_g_fs + 4);

    auto gr_xzz_0 = pbuffer.data(idx_g_fs + 5);

    auto gr_yyy_0 = pbuffer.data(idx_g_fs + 6);

    auto gr_yyz_0 = pbuffer.data(idx_g_fs + 7);

    auto gr_yzz_0 = pbuffer.data(idx_g_fs + 8);

    auto gr_zzz_0 = pbuffer.data(idx_g_fs + 9);

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

    // Set up components of auxiliary buffer : FP

    auto gr_xxx_x = pbuffer.data(idx_g_fp);

    auto gr_xxx_y = pbuffer.data(idx_g_fp + 1);

    auto gr_xxx_z = pbuffer.data(idx_g_fp + 2);

    auto gr_xxy_x = pbuffer.data(idx_g_fp + 3);

    auto gr_xxy_y = pbuffer.data(idx_g_fp + 4);

    auto gr_xxy_z = pbuffer.data(idx_g_fp + 5);

    auto gr_xxz_x = pbuffer.data(idx_g_fp + 6);

    auto gr_xxz_y = pbuffer.data(idx_g_fp + 7);

    auto gr_xxz_z = pbuffer.data(idx_g_fp + 8);

    auto gr_xyy_x = pbuffer.data(idx_g_fp + 9);

    auto gr_xyy_y = pbuffer.data(idx_g_fp + 10);

    auto gr_xyy_z = pbuffer.data(idx_g_fp + 11);

    auto gr_xyz_x = pbuffer.data(idx_g_fp + 12);

    auto gr_xyz_y = pbuffer.data(idx_g_fp + 13);

    auto gr_xyz_z = pbuffer.data(idx_g_fp + 14);

    auto gr_xzz_x = pbuffer.data(idx_g_fp + 15);

    auto gr_xzz_y = pbuffer.data(idx_g_fp + 16);

    auto gr_xzz_z = pbuffer.data(idx_g_fp + 17);

    auto gr_yyy_x = pbuffer.data(idx_g_fp + 18);

    auto gr_yyy_y = pbuffer.data(idx_g_fp + 19);

    auto gr_yyy_z = pbuffer.data(idx_g_fp + 20);

    auto gr_yyz_x = pbuffer.data(idx_g_fp + 21);

    auto gr_yyz_y = pbuffer.data(idx_g_fp + 22);

    auto gr_yyz_z = pbuffer.data(idx_g_fp + 23);

    auto gr_yzz_x = pbuffer.data(idx_g_fp + 24);

    auto gr_yzz_y = pbuffer.data(idx_g_fp + 25);

    auto gr_yzz_z = pbuffer.data(idx_g_fp + 26);

    auto gr_zzz_x = pbuffer.data(idx_g_fp + 27);

    auto gr_zzz_y = pbuffer.data(idx_g_fp + 28);

    auto gr_zzz_z = pbuffer.data(idx_g_fp + 29);

    // Set up 0-3 components of targeted buffer : FP

    auto grr_x_xxx_x = pbuffer.data(idx_gr_fp);

    auto grr_x_xxx_y = pbuffer.data(idx_gr_fp + 1);

    auto grr_x_xxx_z = pbuffer.data(idx_gr_fp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_x, gr_xx_y, gr_xx_z, gr_xxx_0, gr_xxx_x, gr_xxx_y, gr_xxx_z, grr_x_xxx_x, grr_x_xxx_y, grr_x_xxx_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxx_x[i] = 3.0 * ts_xx_x[i] * gfe2_0 + 3.0 * gr_xx_x[i] * gfe_0 + ts_xxx_0[i] * gfe2_0 + gr_xxx_0[i] * gfe_0 + ts_xxx_x[i] * gfe_0 * gc_x[i] + gr_xxx_x[i] * gc_x[i];

        grr_x_xxx_y[i] = 3.0 * ts_xx_y[i] * gfe2_0 + 3.0 * gr_xx_y[i] * gfe_0 + ts_xxx_y[i] * gfe_0 * gc_x[i] + gr_xxx_y[i] * gc_x[i];

        grr_x_xxx_z[i] = 3.0 * ts_xx_z[i] * gfe2_0 + 3.0 * gr_xx_z[i] * gfe_0 + ts_xxx_z[i] * gfe_0 * gc_x[i] + gr_xxx_z[i] * gc_x[i];
    }

    // Set up 3-6 components of targeted buffer : FP

    auto grr_x_xxy_x = pbuffer.data(idx_gr_fp + 3);

    auto grr_x_xxy_y = pbuffer.data(idx_gr_fp + 4);

    auto grr_x_xxy_z = pbuffer.data(idx_gr_fp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_0, gr_xxy_x, gr_xxy_y, gr_xxy_z, gr_xy_x, gr_xy_y, gr_xy_z, grr_x_xxy_x, grr_x_xxy_y, grr_x_xxy_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxy_x[i] = 2.0 * ts_xy_x[i] * gfe2_0 + 2.0 * gr_xy_x[i] * gfe_0 + ts_xxy_0[i] * gfe2_0 + gr_xxy_0[i] * gfe_0 + ts_xxy_x[i] * gfe_0 * gc_x[i] + gr_xxy_x[i] * gc_x[i];

        grr_x_xxy_y[i] = 2.0 * ts_xy_y[i] * gfe2_0 + 2.0 * gr_xy_y[i] * gfe_0 + ts_xxy_y[i] * gfe_0 * gc_x[i] + gr_xxy_y[i] * gc_x[i];

        grr_x_xxy_z[i] = 2.0 * ts_xy_z[i] * gfe2_0 + 2.0 * gr_xy_z[i] * gfe_0 + ts_xxy_z[i] * gfe_0 * gc_x[i] + gr_xxy_z[i] * gc_x[i];
    }

    // Set up 6-9 components of targeted buffer : FP

    auto grr_x_xxz_x = pbuffer.data(idx_gr_fp + 6);

    auto grr_x_xxz_y = pbuffer.data(idx_gr_fp + 7);

    auto grr_x_xxz_z = pbuffer.data(idx_gr_fp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_0, gr_xxz_x, gr_xxz_y, gr_xxz_z, gr_xz_x, gr_xz_y, gr_xz_z, grr_x_xxz_x, grr_x_xxz_y, grr_x_xxz_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxz_x[i] = 2.0 * ts_xz_x[i] * gfe2_0 + 2.0 * gr_xz_x[i] * gfe_0 + ts_xxz_0[i] * gfe2_0 + gr_xxz_0[i] * gfe_0 + ts_xxz_x[i] * gfe_0 * gc_x[i] + gr_xxz_x[i] * gc_x[i];

        grr_x_xxz_y[i] = 2.0 * ts_xz_y[i] * gfe2_0 + 2.0 * gr_xz_y[i] * gfe_0 + ts_xxz_y[i] * gfe_0 * gc_x[i] + gr_xxz_y[i] * gc_x[i];

        grr_x_xxz_z[i] = 2.0 * ts_xz_z[i] * gfe2_0 + 2.0 * gr_xz_z[i] * gfe_0 + ts_xxz_z[i] * gfe_0 * gc_x[i] + gr_xxz_z[i] * gc_x[i];
    }

    // Set up 9-12 components of targeted buffer : FP

    auto grr_x_xyy_x = pbuffer.data(idx_gr_fp + 9);

    auto grr_x_xyy_y = pbuffer.data(idx_gr_fp + 10);

    auto grr_x_xyy_z = pbuffer.data(idx_gr_fp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_0, gr_xyy_x, gr_xyy_y, gr_xyy_z, gr_yy_x, gr_yy_y, gr_yy_z, grr_x_xyy_x, grr_x_xyy_y, grr_x_xyy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyy_x[i] = ts_yy_x[i] * gfe2_0 + gr_yy_x[i] * gfe_0 + ts_xyy_0[i] * gfe2_0 + gr_xyy_0[i] * gfe_0 + ts_xyy_x[i] * gfe_0 * gc_x[i] + gr_xyy_x[i] * gc_x[i];

        grr_x_xyy_y[i] = ts_yy_y[i] * gfe2_0 + gr_yy_y[i] * gfe_0 + ts_xyy_y[i] * gfe_0 * gc_x[i] + gr_xyy_y[i] * gc_x[i];

        grr_x_xyy_z[i] = ts_yy_z[i] * gfe2_0 + gr_yy_z[i] * gfe_0 + ts_xyy_z[i] * gfe_0 * gc_x[i] + gr_xyy_z[i] * gc_x[i];
    }

    // Set up 12-15 components of targeted buffer : FP

    auto grr_x_xyz_x = pbuffer.data(idx_gr_fp + 12);

    auto grr_x_xyz_y = pbuffer.data(idx_gr_fp + 13);

    auto grr_x_xyz_z = pbuffer.data(idx_gr_fp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_0, gr_xyz_x, gr_xyz_y, gr_xyz_z, gr_yz_x, gr_yz_y, gr_yz_z, grr_x_xyz_x, grr_x_xyz_y, grr_x_xyz_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyz_x[i] = ts_yz_x[i] * gfe2_0 + gr_yz_x[i] * gfe_0 + ts_xyz_0[i] * gfe2_0 + gr_xyz_0[i] * gfe_0 + ts_xyz_x[i] * gfe_0 * gc_x[i] + gr_xyz_x[i] * gc_x[i];

        grr_x_xyz_y[i] = ts_yz_y[i] * gfe2_0 + gr_yz_y[i] * gfe_0 + ts_xyz_y[i] * gfe_0 * gc_x[i] + gr_xyz_y[i] * gc_x[i];

        grr_x_xyz_z[i] = ts_yz_z[i] * gfe2_0 + gr_yz_z[i] * gfe_0 + ts_xyz_z[i] * gfe_0 * gc_x[i] + gr_xyz_z[i] * gc_x[i];
    }

    // Set up 15-18 components of targeted buffer : FP

    auto grr_x_xzz_x = pbuffer.data(idx_gr_fp + 15);

    auto grr_x_xzz_y = pbuffer.data(idx_gr_fp + 16);

    auto grr_x_xzz_z = pbuffer.data(idx_gr_fp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_0, gr_xzz_x, gr_xzz_y, gr_xzz_z, gr_zz_x, gr_zz_y, gr_zz_z, grr_x_xzz_x, grr_x_xzz_y, grr_x_xzz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xzz_x[i] = ts_zz_x[i] * gfe2_0 + gr_zz_x[i] * gfe_0 + ts_xzz_0[i] * gfe2_0 + gr_xzz_0[i] * gfe_0 + ts_xzz_x[i] * gfe_0 * gc_x[i] + gr_xzz_x[i] * gc_x[i];

        grr_x_xzz_y[i] = ts_zz_y[i] * gfe2_0 + gr_zz_y[i] * gfe_0 + ts_xzz_y[i] * gfe_0 * gc_x[i] + gr_xzz_y[i] * gc_x[i];

        grr_x_xzz_z[i] = ts_zz_z[i] * gfe2_0 + gr_zz_z[i] * gfe_0 + ts_xzz_z[i] * gfe_0 * gc_x[i] + gr_xzz_z[i] * gc_x[i];
    }

    // Set up 18-21 components of targeted buffer : FP

    auto grr_x_yyy_x = pbuffer.data(idx_gr_fp + 18);

    auto grr_x_yyy_y = pbuffer.data(idx_gr_fp + 19);

    auto grr_x_yyy_z = pbuffer.data(idx_gr_fp + 20);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_0, gr_yyy_x, gr_yyy_y, gr_yyy_z, grr_x_yyy_x, grr_x_yyy_y, grr_x_yyy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyy_x[i] = ts_yyy_0[i] * gfe2_0 + gr_yyy_0[i] * gfe_0 + ts_yyy_x[i] * gfe_0 * gc_x[i] + gr_yyy_x[i] * gc_x[i];

        grr_x_yyy_y[i] = ts_yyy_y[i] * gfe_0 * gc_x[i] + gr_yyy_y[i] * gc_x[i];

        grr_x_yyy_z[i] = ts_yyy_z[i] * gfe_0 * gc_x[i] + gr_yyy_z[i] * gc_x[i];
    }

    // Set up 21-24 components of targeted buffer : FP

    auto grr_x_yyz_x = pbuffer.data(idx_gr_fp + 21);

    auto grr_x_yyz_y = pbuffer.data(idx_gr_fp + 22);

    auto grr_x_yyz_z = pbuffer.data(idx_gr_fp + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_0, gr_yyz_x, gr_yyz_y, gr_yyz_z, grr_x_yyz_x, grr_x_yyz_y, grr_x_yyz_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyz_x[i] = ts_yyz_0[i] * gfe2_0 + gr_yyz_0[i] * gfe_0 + ts_yyz_x[i] * gfe_0 * gc_x[i] + gr_yyz_x[i] * gc_x[i];

        grr_x_yyz_y[i] = ts_yyz_y[i] * gfe_0 * gc_x[i] + gr_yyz_y[i] * gc_x[i];

        grr_x_yyz_z[i] = ts_yyz_z[i] * gfe_0 * gc_x[i] + gr_yyz_z[i] * gc_x[i];
    }

    // Set up 24-27 components of targeted buffer : FP

    auto grr_x_yzz_x = pbuffer.data(idx_gr_fp + 24);

    auto grr_x_yzz_y = pbuffer.data(idx_gr_fp + 25);

    auto grr_x_yzz_z = pbuffer.data(idx_gr_fp + 26);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_0, gr_yzz_x, gr_yzz_y, gr_yzz_z, grr_x_yzz_x, grr_x_yzz_y, grr_x_yzz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yzz_x[i] = ts_yzz_0[i] * gfe2_0 + gr_yzz_0[i] * gfe_0 + ts_yzz_x[i] * gfe_0 * gc_x[i] + gr_yzz_x[i] * gc_x[i];

        grr_x_yzz_y[i] = ts_yzz_y[i] * gfe_0 * gc_x[i] + gr_yzz_y[i] * gc_x[i];

        grr_x_yzz_z[i] = ts_yzz_z[i] * gfe_0 * gc_x[i] + gr_yzz_z[i] * gc_x[i];
    }

    // Set up 27-30 components of targeted buffer : FP

    auto grr_x_zzz_x = pbuffer.data(idx_gr_fp + 27);

    auto grr_x_zzz_y = pbuffer.data(idx_gr_fp + 28);

    auto grr_x_zzz_z = pbuffer.data(idx_gr_fp + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_0, gr_zzz_x, gr_zzz_y, gr_zzz_z, grr_x_zzz_x, grr_x_zzz_y, grr_x_zzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zzz_x[i] = ts_zzz_0[i] * gfe2_0 + gr_zzz_0[i] * gfe_0 + ts_zzz_x[i] * gfe_0 * gc_x[i] + gr_zzz_x[i] * gc_x[i];

        grr_x_zzz_y[i] = ts_zzz_y[i] * gfe_0 * gc_x[i] + gr_zzz_y[i] * gc_x[i];

        grr_x_zzz_z[i] = ts_zzz_z[i] * gfe_0 * gc_x[i] + gr_zzz_z[i] * gc_x[i];
    }

    // Set up 30-33 components of targeted buffer : FP

    auto grr_y_xxx_x = pbuffer.data(idx_gr_fp + 30);

    auto grr_y_xxx_y = pbuffer.data(idx_gr_fp + 31);

    auto grr_y_xxx_z = pbuffer.data(idx_gr_fp + 32);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_0, gr_xxx_x, gr_xxx_y, gr_xxx_z, grr_y_xxx_x, grr_y_xxx_y, grr_y_xxx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxx_x[i] = ts_xxx_x[i] * gfe_0 * gc_y[i] + gr_xxx_x[i] * gc_y[i];

        grr_y_xxx_y[i] = ts_xxx_0[i] * gfe2_0 + gr_xxx_0[i] * gfe_0 + ts_xxx_y[i] * gfe_0 * gc_y[i] + gr_xxx_y[i] * gc_y[i];

        grr_y_xxx_z[i] = ts_xxx_z[i] * gfe_0 * gc_y[i] + gr_xxx_z[i] * gc_y[i];
    }

    // Set up 33-36 components of targeted buffer : FP

    auto grr_y_xxy_x = pbuffer.data(idx_gr_fp + 33);

    auto grr_y_xxy_y = pbuffer.data(idx_gr_fp + 34);

    auto grr_y_xxy_z = pbuffer.data(idx_gr_fp + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_x, gr_xx_y, gr_xx_z, gr_xxy_0, gr_xxy_x, gr_xxy_y, gr_xxy_z, grr_y_xxy_x, grr_y_xxy_y, grr_y_xxy_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxy_x[i] = ts_xx_x[i] * gfe2_0 + gr_xx_x[i] * gfe_0 + ts_xxy_x[i] * gfe_0 * gc_y[i] + gr_xxy_x[i] * gc_y[i];

        grr_y_xxy_y[i] = ts_xx_y[i] * gfe2_0 + gr_xx_y[i] * gfe_0 + ts_xxy_0[i] * gfe2_0 + gr_xxy_0[i] * gfe_0 + ts_xxy_y[i] * gfe_0 * gc_y[i] + gr_xxy_y[i] * gc_y[i];

        grr_y_xxy_z[i] = ts_xx_z[i] * gfe2_0 + gr_xx_z[i] * gfe_0 + ts_xxy_z[i] * gfe_0 * gc_y[i] + gr_xxy_z[i] * gc_y[i];
    }

    // Set up 36-39 components of targeted buffer : FP

    auto grr_y_xxz_x = pbuffer.data(idx_gr_fp + 36);

    auto grr_y_xxz_y = pbuffer.data(idx_gr_fp + 37);

    auto grr_y_xxz_z = pbuffer.data(idx_gr_fp + 38);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_0, gr_xxz_x, gr_xxz_y, gr_xxz_z, grr_y_xxz_x, grr_y_xxz_y, grr_y_xxz_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxz_x[i] = ts_xxz_x[i] * gfe_0 * gc_y[i] + gr_xxz_x[i] * gc_y[i];

        grr_y_xxz_y[i] = ts_xxz_0[i] * gfe2_0 + gr_xxz_0[i] * gfe_0 + ts_xxz_y[i] * gfe_0 * gc_y[i] + gr_xxz_y[i] * gc_y[i];

        grr_y_xxz_z[i] = ts_xxz_z[i] * gfe_0 * gc_y[i] + gr_xxz_z[i] * gc_y[i];
    }

    // Set up 39-42 components of targeted buffer : FP

    auto grr_y_xyy_x = pbuffer.data(idx_gr_fp + 39);

    auto grr_y_xyy_y = pbuffer.data(idx_gr_fp + 40);

    auto grr_y_xyy_z = pbuffer.data(idx_gr_fp + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_x, gr_xy_y, gr_xy_z, gr_xyy_0, gr_xyy_x, gr_xyy_y, gr_xyy_z, grr_y_xyy_x, grr_y_xyy_y, grr_y_xyy_z, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyy_x[i] = 2.0 * ts_xy_x[i] * gfe2_0 + 2.0 * gr_xy_x[i] * gfe_0 + ts_xyy_x[i] * gfe_0 * gc_y[i] + gr_xyy_x[i] * gc_y[i];

        grr_y_xyy_y[i] = 2.0 * ts_xy_y[i] * gfe2_0 + 2.0 * gr_xy_y[i] * gfe_0 + ts_xyy_0[i] * gfe2_0 + gr_xyy_0[i] * gfe_0 + ts_xyy_y[i] * gfe_0 * gc_y[i] + gr_xyy_y[i] * gc_y[i];

        grr_y_xyy_z[i] = 2.0 * ts_xy_z[i] * gfe2_0 + 2.0 * gr_xy_z[i] * gfe_0 + ts_xyy_z[i] * gfe_0 * gc_y[i] + gr_xyy_z[i] * gc_y[i];
    }

    // Set up 42-45 components of targeted buffer : FP

    auto grr_y_xyz_x = pbuffer.data(idx_gr_fp + 42);

    auto grr_y_xyz_y = pbuffer.data(idx_gr_fp + 43);

    auto grr_y_xyz_z = pbuffer.data(idx_gr_fp + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_0, gr_xyz_x, gr_xyz_y, gr_xyz_z, gr_xz_x, gr_xz_y, gr_xz_z, grr_y_xyz_x, grr_y_xyz_y, grr_y_xyz_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyz_x[i] = ts_xz_x[i] * gfe2_0 + gr_xz_x[i] * gfe_0 + ts_xyz_x[i] * gfe_0 * gc_y[i] + gr_xyz_x[i] * gc_y[i];

        grr_y_xyz_y[i] = ts_xz_y[i] * gfe2_0 + gr_xz_y[i] * gfe_0 + ts_xyz_0[i] * gfe2_0 + gr_xyz_0[i] * gfe_0 + ts_xyz_y[i] * gfe_0 * gc_y[i] + gr_xyz_y[i] * gc_y[i];

        grr_y_xyz_z[i] = ts_xz_z[i] * gfe2_0 + gr_xz_z[i] * gfe_0 + ts_xyz_z[i] * gfe_0 * gc_y[i] + gr_xyz_z[i] * gc_y[i];
    }

    // Set up 45-48 components of targeted buffer : FP

    auto grr_y_xzz_x = pbuffer.data(idx_gr_fp + 45);

    auto grr_y_xzz_y = pbuffer.data(idx_gr_fp + 46);

    auto grr_y_xzz_z = pbuffer.data(idx_gr_fp + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_0, gr_xzz_x, gr_xzz_y, gr_xzz_z, grr_y_xzz_x, grr_y_xzz_y, grr_y_xzz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xzz_x[i] = ts_xzz_x[i] * gfe_0 * gc_y[i] + gr_xzz_x[i] * gc_y[i];

        grr_y_xzz_y[i] = ts_xzz_0[i] * gfe2_0 + gr_xzz_0[i] * gfe_0 + ts_xzz_y[i] * gfe_0 * gc_y[i] + gr_xzz_y[i] * gc_y[i];

        grr_y_xzz_z[i] = ts_xzz_z[i] * gfe_0 * gc_y[i] + gr_xzz_z[i] * gc_y[i];
    }

    // Set up 48-51 components of targeted buffer : FP

    auto grr_y_yyy_x = pbuffer.data(idx_gr_fp + 48);

    auto grr_y_yyy_y = pbuffer.data(idx_gr_fp + 49);

    auto grr_y_yyy_z = pbuffer.data(idx_gr_fp + 50);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_x, gr_yy_y, gr_yy_z, gr_yyy_0, gr_yyy_x, gr_yyy_y, gr_yyy_z, grr_y_yyy_x, grr_y_yyy_y, grr_y_yyy_z, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyy_x[i] = 3.0 * ts_yy_x[i] * gfe2_0 + 3.0 * gr_yy_x[i] * gfe_0 + ts_yyy_x[i] * gfe_0 * gc_y[i] + gr_yyy_x[i] * gc_y[i];

        grr_y_yyy_y[i] = 3.0 * ts_yy_y[i] * gfe2_0 + 3.0 * gr_yy_y[i] * gfe_0 + ts_yyy_0[i] * gfe2_0 + gr_yyy_0[i] * gfe_0 + ts_yyy_y[i] * gfe_0 * gc_y[i] + gr_yyy_y[i] * gc_y[i];

        grr_y_yyy_z[i] = 3.0 * ts_yy_z[i] * gfe2_0 + 3.0 * gr_yy_z[i] * gfe_0 + ts_yyy_z[i] * gfe_0 * gc_y[i] + gr_yyy_z[i] * gc_y[i];
    }

    // Set up 51-54 components of targeted buffer : FP

    auto grr_y_yyz_x = pbuffer.data(idx_gr_fp + 51);

    auto grr_y_yyz_y = pbuffer.data(idx_gr_fp + 52);

    auto grr_y_yyz_z = pbuffer.data(idx_gr_fp + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_0, gr_yyz_x, gr_yyz_y, gr_yyz_z, gr_yz_x, gr_yz_y, gr_yz_z, grr_y_yyz_x, grr_y_yyz_y, grr_y_yyz_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyz_x[i] = 2.0 * ts_yz_x[i] * gfe2_0 + 2.0 * gr_yz_x[i] * gfe_0 + ts_yyz_x[i] * gfe_0 * gc_y[i] + gr_yyz_x[i] * gc_y[i];

        grr_y_yyz_y[i] = 2.0 * ts_yz_y[i] * gfe2_0 + 2.0 * gr_yz_y[i] * gfe_0 + ts_yyz_0[i] * gfe2_0 + gr_yyz_0[i] * gfe_0 + ts_yyz_y[i] * gfe_0 * gc_y[i] + gr_yyz_y[i] * gc_y[i];

        grr_y_yyz_z[i] = 2.0 * ts_yz_z[i] * gfe2_0 + 2.0 * gr_yz_z[i] * gfe_0 + ts_yyz_z[i] * gfe_0 * gc_y[i] + gr_yyz_z[i] * gc_y[i];
    }

    // Set up 54-57 components of targeted buffer : FP

    auto grr_y_yzz_x = pbuffer.data(idx_gr_fp + 54);

    auto grr_y_yzz_y = pbuffer.data(idx_gr_fp + 55);

    auto grr_y_yzz_z = pbuffer.data(idx_gr_fp + 56);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_0, gr_yzz_x, gr_yzz_y, gr_yzz_z, gr_zz_x, gr_zz_y, gr_zz_z, grr_y_yzz_x, grr_y_yzz_y, grr_y_yzz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yzz_x[i] = ts_zz_x[i] * gfe2_0 + gr_zz_x[i] * gfe_0 + ts_yzz_x[i] * gfe_0 * gc_y[i] + gr_yzz_x[i] * gc_y[i];

        grr_y_yzz_y[i] = ts_zz_y[i] * gfe2_0 + gr_zz_y[i] * gfe_0 + ts_yzz_0[i] * gfe2_0 + gr_yzz_0[i] * gfe_0 + ts_yzz_y[i] * gfe_0 * gc_y[i] + gr_yzz_y[i] * gc_y[i];

        grr_y_yzz_z[i] = ts_zz_z[i] * gfe2_0 + gr_zz_z[i] * gfe_0 + ts_yzz_z[i] * gfe_0 * gc_y[i] + gr_yzz_z[i] * gc_y[i];
    }

    // Set up 57-60 components of targeted buffer : FP

    auto grr_y_zzz_x = pbuffer.data(idx_gr_fp + 57);

    auto grr_y_zzz_y = pbuffer.data(idx_gr_fp + 58);

    auto grr_y_zzz_z = pbuffer.data(idx_gr_fp + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_0, gr_zzz_x, gr_zzz_y, gr_zzz_z, grr_y_zzz_x, grr_y_zzz_y, grr_y_zzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zzz_x[i] = ts_zzz_x[i] * gfe_0 * gc_y[i] + gr_zzz_x[i] * gc_y[i];

        grr_y_zzz_y[i] = ts_zzz_0[i] * gfe2_0 + gr_zzz_0[i] * gfe_0 + ts_zzz_y[i] * gfe_0 * gc_y[i] + gr_zzz_y[i] * gc_y[i];

        grr_y_zzz_z[i] = ts_zzz_z[i] * gfe_0 * gc_y[i] + gr_zzz_z[i] * gc_y[i];
    }

    // Set up 60-63 components of targeted buffer : FP

    auto grr_z_xxx_x = pbuffer.data(idx_gr_fp + 60);

    auto grr_z_xxx_y = pbuffer.data(idx_gr_fp + 61);

    auto grr_z_xxx_z = pbuffer.data(idx_gr_fp + 62);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_0, gr_xxx_x, gr_xxx_y, gr_xxx_z, grr_z_xxx_x, grr_z_xxx_y, grr_z_xxx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxx_x[i] = ts_xxx_x[i] * gfe_0 * gc_z[i] + gr_xxx_x[i] * gc_z[i];

        grr_z_xxx_y[i] = ts_xxx_y[i] * gfe_0 * gc_z[i] + gr_xxx_y[i] * gc_z[i];

        grr_z_xxx_z[i] = ts_xxx_0[i] * gfe2_0 + gr_xxx_0[i] * gfe_0 + ts_xxx_z[i] * gfe_0 * gc_z[i] + gr_xxx_z[i] * gc_z[i];
    }

    // Set up 63-66 components of targeted buffer : FP

    auto grr_z_xxy_x = pbuffer.data(idx_gr_fp + 63);

    auto grr_z_xxy_y = pbuffer.data(idx_gr_fp + 64);

    auto grr_z_xxy_z = pbuffer.data(idx_gr_fp + 65);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_0, gr_xxy_x, gr_xxy_y, gr_xxy_z, grr_z_xxy_x, grr_z_xxy_y, grr_z_xxy_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxy_x[i] = ts_xxy_x[i] * gfe_0 * gc_z[i] + gr_xxy_x[i] * gc_z[i];

        grr_z_xxy_y[i] = ts_xxy_y[i] * gfe_0 * gc_z[i] + gr_xxy_y[i] * gc_z[i];

        grr_z_xxy_z[i] = ts_xxy_0[i] * gfe2_0 + gr_xxy_0[i] * gfe_0 + ts_xxy_z[i] * gfe_0 * gc_z[i] + gr_xxy_z[i] * gc_z[i];
    }

    // Set up 66-69 components of targeted buffer : FP

    auto grr_z_xxz_x = pbuffer.data(idx_gr_fp + 66);

    auto grr_z_xxz_y = pbuffer.data(idx_gr_fp + 67);

    auto grr_z_xxz_z = pbuffer.data(idx_gr_fp + 68);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_x, gr_xx_y, gr_xx_z, gr_xxz_0, gr_xxz_x, gr_xxz_y, gr_xxz_z, grr_z_xxz_x, grr_z_xxz_y, grr_z_xxz_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxz_x[i] = ts_xx_x[i] * gfe2_0 + gr_xx_x[i] * gfe_0 + ts_xxz_x[i] * gfe_0 * gc_z[i] + gr_xxz_x[i] * gc_z[i];

        grr_z_xxz_y[i] = ts_xx_y[i] * gfe2_0 + gr_xx_y[i] * gfe_0 + ts_xxz_y[i] * gfe_0 * gc_z[i] + gr_xxz_y[i] * gc_z[i];

        grr_z_xxz_z[i] = ts_xx_z[i] * gfe2_0 + gr_xx_z[i] * gfe_0 + ts_xxz_0[i] * gfe2_0 + gr_xxz_0[i] * gfe_0 + ts_xxz_z[i] * gfe_0 * gc_z[i] + gr_xxz_z[i] * gc_z[i];
    }

    // Set up 69-72 components of targeted buffer : FP

    auto grr_z_xyy_x = pbuffer.data(idx_gr_fp + 69);

    auto grr_z_xyy_y = pbuffer.data(idx_gr_fp + 70);

    auto grr_z_xyy_z = pbuffer.data(idx_gr_fp + 71);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_0, gr_xyy_x, gr_xyy_y, gr_xyy_z, grr_z_xyy_x, grr_z_xyy_y, grr_z_xyy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyy_x[i] = ts_xyy_x[i] * gfe_0 * gc_z[i] + gr_xyy_x[i] * gc_z[i];

        grr_z_xyy_y[i] = ts_xyy_y[i] * gfe_0 * gc_z[i] + gr_xyy_y[i] * gc_z[i];

        grr_z_xyy_z[i] = ts_xyy_0[i] * gfe2_0 + gr_xyy_0[i] * gfe_0 + ts_xyy_z[i] * gfe_0 * gc_z[i] + gr_xyy_z[i] * gc_z[i];
    }

    // Set up 72-75 components of targeted buffer : FP

    auto grr_z_xyz_x = pbuffer.data(idx_gr_fp + 72);

    auto grr_z_xyz_y = pbuffer.data(idx_gr_fp + 73);

    auto grr_z_xyz_z = pbuffer.data(idx_gr_fp + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xy_x, gr_xy_y, gr_xy_z, gr_xyz_0, gr_xyz_x, gr_xyz_y, gr_xyz_z, grr_z_xyz_x, grr_z_xyz_y, grr_z_xyz_z, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyz_x[i] = ts_xy_x[i] * gfe2_0 + gr_xy_x[i] * gfe_0 + ts_xyz_x[i] * gfe_0 * gc_z[i] + gr_xyz_x[i] * gc_z[i];

        grr_z_xyz_y[i] = ts_xy_y[i] * gfe2_0 + gr_xy_y[i] * gfe_0 + ts_xyz_y[i] * gfe_0 * gc_z[i] + gr_xyz_y[i] * gc_z[i];

        grr_z_xyz_z[i] = ts_xy_z[i] * gfe2_0 + gr_xy_z[i] * gfe_0 + ts_xyz_0[i] * gfe2_0 + gr_xyz_0[i] * gfe_0 + ts_xyz_z[i] * gfe_0 * gc_z[i] + gr_xyz_z[i] * gc_z[i];
    }

    // Set up 75-78 components of targeted buffer : FP

    auto grr_z_xzz_x = pbuffer.data(idx_gr_fp + 75);

    auto grr_z_xzz_y = pbuffer.data(idx_gr_fp + 76);

    auto grr_z_xzz_z = pbuffer.data(idx_gr_fp + 77);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xz_x, gr_xz_y, gr_xz_z, gr_xzz_0, gr_xzz_x, gr_xzz_y, gr_xzz_z, grr_z_xzz_x, grr_z_xzz_y, grr_z_xzz_z, ts_xz_x, ts_xz_y, ts_xz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xzz_x[i] = 2.0 * ts_xz_x[i] * gfe2_0 + 2.0 * gr_xz_x[i] * gfe_0 + ts_xzz_x[i] * gfe_0 * gc_z[i] + gr_xzz_x[i] * gc_z[i];

        grr_z_xzz_y[i] = 2.0 * ts_xz_y[i] * gfe2_0 + 2.0 * gr_xz_y[i] * gfe_0 + ts_xzz_y[i] * gfe_0 * gc_z[i] + gr_xzz_y[i] * gc_z[i];

        grr_z_xzz_z[i] = 2.0 * ts_xz_z[i] * gfe2_0 + 2.0 * gr_xz_z[i] * gfe_0 + ts_xzz_0[i] * gfe2_0 + gr_xzz_0[i] * gfe_0 + ts_xzz_z[i] * gfe_0 * gc_z[i] + gr_xzz_z[i] * gc_z[i];
    }

    // Set up 78-81 components of targeted buffer : FP

    auto grr_z_yyy_x = pbuffer.data(idx_gr_fp + 78);

    auto grr_z_yyy_y = pbuffer.data(idx_gr_fp + 79);

    auto grr_z_yyy_z = pbuffer.data(idx_gr_fp + 80);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_0, gr_yyy_x, gr_yyy_y, gr_yyy_z, grr_z_yyy_x, grr_z_yyy_y, grr_z_yyy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyy_x[i] = ts_yyy_x[i] * gfe_0 * gc_z[i] + gr_yyy_x[i] * gc_z[i];

        grr_z_yyy_y[i] = ts_yyy_y[i] * gfe_0 * gc_z[i] + gr_yyy_y[i] * gc_z[i];

        grr_z_yyy_z[i] = ts_yyy_0[i] * gfe2_0 + gr_yyy_0[i] * gfe_0 + ts_yyy_z[i] * gfe_0 * gc_z[i] + gr_yyy_z[i] * gc_z[i];
    }

    // Set up 81-84 components of targeted buffer : FP

    auto grr_z_yyz_x = pbuffer.data(idx_gr_fp + 81);

    auto grr_z_yyz_y = pbuffer.data(idx_gr_fp + 82);

    auto grr_z_yyz_z = pbuffer.data(idx_gr_fp + 83);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yy_x, gr_yy_y, gr_yy_z, gr_yyz_0, gr_yyz_x, gr_yyz_y, gr_yyz_z, grr_z_yyz_x, grr_z_yyz_y, grr_z_yyz_z, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyz_x[i] = ts_yy_x[i] * gfe2_0 + gr_yy_x[i] * gfe_0 + ts_yyz_x[i] * gfe_0 * gc_z[i] + gr_yyz_x[i] * gc_z[i];

        grr_z_yyz_y[i] = ts_yy_y[i] * gfe2_0 + gr_yy_y[i] * gfe_0 + ts_yyz_y[i] * gfe_0 * gc_z[i] + gr_yyz_y[i] * gc_z[i];

        grr_z_yyz_z[i] = ts_yy_z[i] * gfe2_0 + gr_yy_z[i] * gfe_0 + ts_yyz_0[i] * gfe2_0 + gr_yyz_0[i] * gfe_0 + ts_yyz_z[i] * gfe_0 * gc_z[i] + gr_yyz_z[i] * gc_z[i];
    }

    // Set up 84-87 components of targeted buffer : FP

    auto grr_z_yzz_x = pbuffer.data(idx_gr_fp + 84);

    auto grr_z_yzz_y = pbuffer.data(idx_gr_fp + 85);

    auto grr_z_yzz_z = pbuffer.data(idx_gr_fp + 86);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yz_x, gr_yz_y, gr_yz_z, gr_yzz_0, gr_yzz_x, gr_yzz_y, gr_yzz_z, grr_z_yzz_x, grr_z_yzz_y, grr_z_yzz_z, ts_yz_x, ts_yz_y, ts_yz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yzz_x[i] = 2.0 * ts_yz_x[i] * gfe2_0 + 2.0 * gr_yz_x[i] * gfe_0 + ts_yzz_x[i] * gfe_0 * gc_z[i] + gr_yzz_x[i] * gc_z[i];

        grr_z_yzz_y[i] = 2.0 * ts_yz_y[i] * gfe2_0 + 2.0 * gr_yz_y[i] * gfe_0 + ts_yzz_y[i] * gfe_0 * gc_z[i] + gr_yzz_y[i] * gc_z[i];

        grr_z_yzz_z[i] = 2.0 * ts_yz_z[i] * gfe2_0 + 2.0 * gr_yz_z[i] * gfe_0 + ts_yzz_0[i] * gfe2_0 + gr_yzz_0[i] * gfe_0 + ts_yzz_z[i] * gfe_0 * gc_z[i] + gr_yzz_z[i] * gc_z[i];
    }

    // Set up 87-90 components of targeted buffer : FP

    auto grr_z_zzz_x = pbuffer.data(idx_gr_fp + 87);

    auto grr_z_zzz_y = pbuffer.data(idx_gr_fp + 88);

    auto grr_z_zzz_z = pbuffer.data(idx_gr_fp + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zz_x, gr_zz_y, gr_zz_z, gr_zzz_0, gr_zzz_x, gr_zzz_y, gr_zzz_z, grr_z_zzz_x, grr_z_zzz_y, grr_z_zzz_z, ts_zz_x, ts_zz_y, ts_zz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zzz_x[i] = 3.0 * ts_zz_x[i] * gfe2_0 + 3.0 * gr_zz_x[i] * gfe_0 + ts_zzz_x[i] * gfe_0 * gc_z[i] + gr_zzz_x[i] * gc_z[i];

        grr_z_zzz_y[i] = 3.0 * ts_zz_y[i] * gfe2_0 + 3.0 * gr_zz_y[i] * gfe_0 + ts_zzz_y[i] * gfe_0 * gc_z[i] + gr_zzz_y[i] * gc_z[i];

        grr_z_zzz_z[i] = 3.0 * ts_zz_z[i] * gfe2_0 + 3.0 * gr_zz_z[i] * gfe_0 + ts_zzz_0[i] * gfe2_0 + gr_zzz_0[i] * gfe_0 + ts_zzz_z[i] * gfe_0 * gc_z[i] + gr_zzz_z[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

