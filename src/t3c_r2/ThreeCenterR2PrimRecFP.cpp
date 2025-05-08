#include "ThreeCenterR2PrimRecFP.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_fp(CSimdArray<double>& pbuffer, 
                const size_t idx_g_fp,
                const size_t idx_pp,
                const size_t idx_ds,
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

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_pp);

    auto ts_x_y = pbuffer.data(idx_pp + 1);

    auto ts_x_z = pbuffer.data(idx_pp + 2);

    auto ts_y_x = pbuffer.data(idx_pp + 3);

    auto ts_y_y = pbuffer.data(idx_pp + 4);

    auto ts_y_z = pbuffer.data(idx_pp + 5);

    auto ts_z_x = pbuffer.data(idx_pp + 6);

    auto ts_z_y = pbuffer.data(idx_pp + 7);

    auto ts_z_z = pbuffer.data(idx_pp + 8);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

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

    auto gr_xxx_x = pbuffer.data(idx_g_fp);

    auto gr_xxx_y = pbuffer.data(idx_g_fp + 1);

    auto gr_xxx_z = pbuffer.data(idx_g_fp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_x, gr_xxx_y, gr_xxx_z, ts_x_x, ts_x_y, ts_x_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxx_x[i] = 6.0 * ts_x_x[i] * gfe_0 + 6.0 * ts_xx_0[i] * gfe_0 + 6.0 * ts_xx_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxx_x[i] * gfe_0 + ts_xxx_x[i] * rgc2_0;

        gr_xxx_y[i] = 6.0 * ts_x_y[i] * gfe_0 + 6.0 * ts_xx_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_y[i] * gfe_0 + ts_xxx_y[i] * rgc2_0;

        gr_xxx_z[i] = 6.0 * ts_x_z[i] * gfe_0 + 6.0 * ts_xx_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_z[i] * gfe_0 + ts_xxx_z[i] * rgc2_0;
    }

    // Set up 3-6 components of targeted buffer : FP

    auto gr_xxy_x = pbuffer.data(idx_g_fp + 3);

    auto gr_xxy_y = pbuffer.data(idx_g_fp + 4);

    auto gr_xxy_z = pbuffer.data(idx_g_fp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_x, gr_xxy_y, gr_xxy_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxy_x[i] = 2.0 * ts_y_x[i] * gfe_0 + 4.0 * ts_xy_0[i] * gfe_0 + 4.0 * ts_xy_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxy_x[i] * gfe_0 + ts_xxy_x[i] * rgc2_0;

        gr_xxy_y[i] = 2.0 * ts_y_y[i] * gfe_0 + 4.0 * ts_xy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 + 2.0 * ts_xx_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_y[i] * gfe_0 + ts_xxy_y[i] * rgc2_0;

        gr_xxy_z[i] = 2.0 * ts_y_z[i] * gfe_0 + 4.0 * ts_xy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_z[i] * gfe_0 + ts_xxy_z[i] * rgc2_0;
    }

    // Set up 6-9 components of targeted buffer : FP

    auto gr_xxz_x = pbuffer.data(idx_g_fp + 6);

    auto gr_xxz_y = pbuffer.data(idx_g_fp + 7);

    auto gr_xxz_z = pbuffer.data(idx_g_fp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_x, gr_xxz_y, gr_xxz_z, ts_xx_0, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxz_x[i] = 2.0 * ts_z_x[i] * gfe_0 + 4.0 * ts_xz_0[i] * gfe_0 + 4.0 * ts_xz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxz_x[i] * gfe_0 + ts_xxz_x[i] * rgc2_0;

        gr_xxz_y[i] = 2.0 * ts_z_y[i] * gfe_0 + 4.0 * ts_xz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_y[i] * gfe_0 + ts_xxz_y[i] * rgc2_0;

        gr_xxz_z[i] = 2.0 * ts_z_z[i] * gfe_0 + 4.0 * ts_xz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 + 2.0 * ts_xx_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_z[i] * gfe_0 + ts_xxz_z[i] * rgc2_0;
    }

    // Set up 9-12 components of targeted buffer : FP

    auto gr_xyy_x = pbuffer.data(idx_g_fp + 9);

    auto gr_xyy_y = pbuffer.data(idx_g_fp + 10);

    auto gr_xyy_z = pbuffer.data(idx_g_fp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_x, gr_xyy_y, gr_xyy_z, ts_x_x, ts_x_y, ts_x_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyy_x[i] = 2.0 * ts_yy_0[i] * gfe_0 + 2.0 * ts_yy_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 + 4.0 * ts_xy_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyy_x[i] * gfe_0 + ts_xyy_x[i] * rgc2_0;

        gr_xyy_y[i] = 2.0 * ts_yy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_y[i] * gfe_0 + 4.0 * ts_xy_0[i] * gfe_0 + 4.0 * ts_xy_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_y[i] * gfe_0 + ts_xyy_y[i] * rgc2_0;

        gr_xyy_z[i] = 2.0 * ts_yy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 + 4.0 * ts_xy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_z[i] * gfe_0 + ts_xyy_z[i] * rgc2_0;
    }

    // Set up 12-15 components of targeted buffer : FP

    auto gr_xyz_x = pbuffer.data(idx_g_fp + 12);

    auto gr_xyz_y = pbuffer.data(idx_g_fp + 13);

    auto gr_xyz_z = pbuffer.data(idx_g_fp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_x, gr_xyz_y, gr_xyz_z, ts_xy_0, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyz_x[i] = 2.0 * ts_yz_0[i] * gfe_0 + 2.0 * ts_yz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyz_x[i] * gfe_0 + ts_xyz_x[i] * rgc2_0;

        gr_xyz_y[i] = 2.0 * ts_yz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_0[i] * gfe_0 + 2.0 * ts_xz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_y[i] * gfe_0 + ts_xyz_y[i] * rgc2_0;

        gr_xyz_z[i] = 2.0 * ts_yz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 + 2.0 * ts_xy_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_z[i] * gfe_0 + ts_xyz_z[i] * rgc2_0;
    }

    // Set up 15-18 components of targeted buffer : FP

    auto gr_xzz_x = pbuffer.data(idx_g_fp + 15);

    auto gr_xzz_y = pbuffer.data(idx_g_fp + 16);

    auto gr_xzz_z = pbuffer.data(idx_g_fp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_x, gr_xzz_y, gr_xzz_z, ts_x_x, ts_x_y, ts_x_z, ts_xz_0, ts_xz_x, ts_xz_y, ts_xz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xzz_x[i] = 2.0 * ts_zz_0[i] * gfe_0 + 2.0 * ts_zz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 + 4.0 * ts_xz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzz_x[i] * gfe_0 + ts_xzz_x[i] * rgc2_0;

        gr_xzz_y[i] = 2.0 * ts_zz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_y[i] * gfe_0 + 4.0 * ts_xz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_y[i] * gfe_0 + ts_xzz_y[i] * rgc2_0;

        gr_xzz_z[i] = 2.0 * ts_zz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 + 4.0 * ts_xz_0[i] * gfe_0 + 4.0 * ts_xz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_z[i] * gfe_0 + ts_xzz_z[i] * rgc2_0;
    }

    // Set up 18-21 components of targeted buffer : FP

    auto gr_yyy_x = pbuffer.data(idx_g_fp + 18);

    auto gr_yyy_y = pbuffer.data(idx_g_fp + 19);

    auto gr_yyy_z = pbuffer.data(idx_g_fp + 20);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_x, gr_yyy_y, gr_yyy_z, ts_y_x, ts_y_y, ts_y_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyy_x[i] = 6.0 * ts_y_x[i] * gfe_0 + 6.0 * ts_yy_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyy_x[i] * gfe_0 + ts_yyy_x[i] * rgc2_0;

        gr_yyy_y[i] = 6.0 * ts_y_y[i] * gfe_0 + 6.0 * ts_yy_0[i] * gfe_0 + 6.0 * ts_yy_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_y[i] * gfe_0 + ts_yyy_y[i] * rgc2_0;

        gr_yyy_z[i] = 6.0 * ts_y_z[i] * gfe_0 + 6.0 * ts_yy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_z[i] * gfe_0 + ts_yyy_z[i] * rgc2_0;
    }

    // Set up 21-24 components of targeted buffer : FP

    auto gr_yyz_x = pbuffer.data(idx_g_fp + 21);

    auto gr_yyz_y = pbuffer.data(idx_g_fp + 22);

    auto gr_yyz_z = pbuffer.data(idx_g_fp + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_x, gr_yyz_y, gr_yyz_z, ts_yy_0, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyz_x[i] = 2.0 * ts_z_x[i] * gfe_0 + 4.0 * ts_yz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyz_x[i] * gfe_0 + ts_yyz_x[i] * rgc2_0;

        gr_yyz_y[i] = 2.0 * ts_z_y[i] * gfe_0 + 4.0 * ts_yz_0[i] * gfe_0 + 4.0 * ts_yz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_y[i] * gfe_0 + ts_yyz_y[i] * rgc2_0;

        gr_yyz_z[i] = 2.0 * ts_z_z[i] * gfe_0 + 4.0 * ts_yz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 + 2.0 * ts_yy_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_z[i] * gfe_0 + ts_yyz_z[i] * rgc2_0;
    }

    // Set up 24-27 components of targeted buffer : FP

    auto gr_yzz_x = pbuffer.data(idx_g_fp + 24);

    auto gr_yzz_y = pbuffer.data(idx_g_fp + 25);

    auto gr_yzz_z = pbuffer.data(idx_g_fp + 26);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_x, gr_yzz_y, gr_yzz_z, ts_y_x, ts_y_y, ts_y_z, ts_yz_0, ts_yz_x, ts_yz_y, ts_yz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yzz_x[i] = 2.0 * ts_zz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_x[i] * gfe_0 + 4.0 * ts_yz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzz_x[i] * gfe_0 + ts_yzz_x[i] * rgc2_0;

        gr_yzz_y[i] = 2.0 * ts_zz_0[i] * gfe_0 + 2.0 * ts_zz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 + 4.0 * ts_yz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_y[i] * gfe_0 + ts_yzz_y[i] * rgc2_0;

        gr_yzz_z[i] = 2.0 * ts_zz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_z[i] * gfe_0 + 4.0 * ts_yz_0[i] * gfe_0 + 4.0 * ts_yz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_z[i] * gfe_0 + ts_yzz_z[i] * rgc2_0;
    }

    // Set up 27-30 components of targeted buffer : FP

    auto gr_zzz_x = pbuffer.data(idx_g_fp + 27);

    auto gr_zzz_y = pbuffer.data(idx_g_fp + 28);

    auto gr_zzz_z = pbuffer.data(idx_g_fp + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_x, gr_zzz_y, gr_zzz_z, ts_z_x, ts_z_y, ts_z_z, ts_zz_0, ts_zz_x, ts_zz_y, ts_zz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zzz_x[i] = 6.0 * ts_z_x[i] * gfe_0 + 6.0 * ts_zz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzz_x[i] * gfe_0 + ts_zzz_x[i] * rgc2_0;

        gr_zzz_y[i] = 6.0 * ts_z_y[i] * gfe_0 + 6.0 * ts_zz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_y[i] * gfe_0 + ts_zzz_y[i] * rgc2_0;

        gr_zzz_z[i] = 6.0 * ts_z_z[i] * gfe_0 + 6.0 * ts_zz_0[i] * gfe_0 + 6.0 * ts_zz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_z[i] * gfe_0 + ts_zzz_z[i] * rgc2_0;
    }

}

} // t3r2rec namespace

