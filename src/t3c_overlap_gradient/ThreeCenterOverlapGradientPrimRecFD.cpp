#include "ThreeCenterOverlapGradientPrimRecFD.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_fd(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_fd,
                              const size_t idx_dd,
                              const size_t idx_fp,
                              const size_t idx_fd,
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

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_dd);

    auto ts_xx_xy = pbuffer.data(idx_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_dd + 5);

    auto ts_xy_xx = pbuffer.data(idx_dd + 6);

    auto ts_xy_xy = pbuffer.data(idx_dd + 7);

    auto ts_xy_xz = pbuffer.data(idx_dd + 8);

    auto ts_xy_yy = pbuffer.data(idx_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_dd + 10);

    auto ts_xy_zz = pbuffer.data(idx_dd + 11);

    auto ts_xz_xx = pbuffer.data(idx_dd + 12);

    auto ts_xz_xy = pbuffer.data(idx_dd + 13);

    auto ts_xz_xz = pbuffer.data(idx_dd + 14);

    auto ts_xz_yy = pbuffer.data(idx_dd + 15);

    auto ts_xz_yz = pbuffer.data(idx_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_dd + 23);

    auto ts_yz_xx = pbuffer.data(idx_dd + 24);

    auto ts_yz_xy = pbuffer.data(idx_dd + 25);

    auto ts_yz_xz = pbuffer.data(idx_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_dd + 35);

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

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_fd);

    auto ts_xxx_xy = pbuffer.data(idx_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_fd + 10);

    auto ts_xxy_zz = pbuffer.data(idx_fd + 11);

    auto ts_xxz_xx = pbuffer.data(idx_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_fd + 14);

    auto ts_xxz_yy = pbuffer.data(idx_fd + 15);

    auto ts_xxz_yz = pbuffer.data(idx_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_fd + 17);

    auto ts_xyy_xx = pbuffer.data(idx_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_fd + 19);

    auto ts_xyy_xz = pbuffer.data(idx_fd + 20);

    auto ts_xyy_yy = pbuffer.data(idx_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_fd + 23);

    auto ts_xyz_xx = pbuffer.data(idx_fd + 24);

    auto ts_xyz_xy = pbuffer.data(idx_fd + 25);

    auto ts_xyz_xz = pbuffer.data(idx_fd + 26);

    auto ts_xyz_yy = pbuffer.data(idx_fd + 27);

    auto ts_xyz_yz = pbuffer.data(idx_fd + 28);

    auto ts_xyz_zz = pbuffer.data(idx_fd + 29);

    auto ts_xzz_xx = pbuffer.data(idx_fd + 30);

    auto ts_xzz_xy = pbuffer.data(idx_fd + 31);

    auto ts_xzz_xz = pbuffer.data(idx_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_fd + 41);

    auto ts_yyz_xx = pbuffer.data(idx_fd + 42);

    auto ts_yyz_xy = pbuffer.data(idx_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_fd + 47);

    auto ts_yzz_xx = pbuffer.data(idx_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up 0-6 components of targeted buffer : FD

    auto gs_x_xxx_xx = pbuffer.data(idx_g_fd);

    auto gs_x_xxx_xy = pbuffer.data(idx_g_fd + 1);

    auto gs_x_xxx_xz = pbuffer.data(idx_g_fd + 2);

    auto gs_x_xxx_yy = pbuffer.data(idx_g_fd + 3);

    auto gs_x_xxx_yz = pbuffer.data(idx_g_fd + 4);

    auto gs_x_xxx_zz = pbuffer.data(idx_g_fd + 5);

    #pragma omp simd aligned(gc_x, gs_x_xxx_xx, gs_x_xxx_xy, gs_x_xxx_xz, gs_x_xxx_yy, gs_x_xxx_yz, gs_x_xxx_zz, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_yy, ts_xx_yz, ts_xx_zz, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxx_xx[i] = 6.0 * ts_xx_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xx[i] * gc_x[i] * tce_0;

        gs_x_xxx_xy[i] = 6.0 * ts_xx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xy[i] * gc_x[i] * tce_0;

        gs_x_xxx_xz[i] = 6.0 * ts_xx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xz[i] * gc_x[i] * tce_0;

        gs_x_xxx_yy[i] = 6.0 * ts_xx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yy[i] * gc_x[i] * tce_0;

        gs_x_xxx_yz[i] = 6.0 * ts_xx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yz[i] * gc_x[i] * tce_0;

        gs_x_xxx_zz[i] = 6.0 * ts_xx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 6-12 components of targeted buffer : FD

    auto gs_x_xxy_xx = pbuffer.data(idx_g_fd + 6);

    auto gs_x_xxy_xy = pbuffer.data(idx_g_fd + 7);

    auto gs_x_xxy_xz = pbuffer.data(idx_g_fd + 8);

    auto gs_x_xxy_yy = pbuffer.data(idx_g_fd + 9);

    auto gs_x_xxy_yz = pbuffer.data(idx_g_fd + 10);

    auto gs_x_xxy_zz = pbuffer.data(idx_g_fd + 11);

    #pragma omp simd aligned(gc_x, gs_x_xxy_xx, gs_x_xxy_xy, gs_x_xxy_xz, gs_x_xxy_yy, gs_x_xxy_yz, gs_x_xxy_zz, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_yy, ts_xy_yz, ts_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxy_xx[i] = 4.0 * ts_xy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xx[i] * gc_x[i] * tce_0;

        gs_x_xxy_xy[i] = 4.0 * ts_xy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xy[i] * gc_x[i] * tce_0;

        gs_x_xxy_xz[i] = 4.0 * ts_xy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xz[i] * gc_x[i] * tce_0;

        gs_x_xxy_yy[i] = 4.0 * ts_xy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yy[i] * gc_x[i] * tce_0;

        gs_x_xxy_yz[i] = 4.0 * ts_xy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yz[i] * gc_x[i] * tce_0;

        gs_x_xxy_zz[i] = 4.0 * ts_xy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 12-18 components of targeted buffer : FD

    auto gs_x_xxz_xx = pbuffer.data(idx_g_fd + 12);

    auto gs_x_xxz_xy = pbuffer.data(idx_g_fd + 13);

    auto gs_x_xxz_xz = pbuffer.data(idx_g_fd + 14);

    auto gs_x_xxz_yy = pbuffer.data(idx_g_fd + 15);

    auto gs_x_xxz_yz = pbuffer.data(idx_g_fd + 16);

    auto gs_x_xxz_zz = pbuffer.data(idx_g_fd + 17);

    #pragma omp simd aligned(gc_x, gs_x_xxz_xx, gs_x_xxz_xy, gs_x_xxz_xz, gs_x_xxz_yy, gs_x_xxz_yz, gs_x_xxz_zz, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_yy, ts_xz_yz, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxz_xx[i] = 4.0 * ts_xz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxz_xy[i] = 4.0 * ts_xz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxz_xz[i] = 4.0 * ts_xz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxz_yy[i] = 4.0 * ts_xz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxz_yz[i] = 4.0 * ts_xz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxz_zz[i] = 4.0 * ts_xz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 18-24 components of targeted buffer : FD

    auto gs_x_xyy_xx = pbuffer.data(idx_g_fd + 18);

    auto gs_x_xyy_xy = pbuffer.data(idx_g_fd + 19);

    auto gs_x_xyy_xz = pbuffer.data(idx_g_fd + 20);

    auto gs_x_xyy_yy = pbuffer.data(idx_g_fd + 21);

    auto gs_x_xyy_yz = pbuffer.data(idx_g_fd + 22);

    auto gs_x_xyy_zz = pbuffer.data(idx_g_fd + 23);

    #pragma omp simd aligned(gc_x, gs_x_xyy_xx, gs_x_xyy_xy, gs_x_xyy_xz, gs_x_xyy_yy, gs_x_xyy_yz, gs_x_xyy_zz, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_yy, ts_yy_yz, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyy_xx[i] = 2.0 * ts_yy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xx[i] * gc_x[i] * tce_0;

        gs_x_xyy_xy[i] = 2.0 * ts_yy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xy[i] * gc_x[i] * tce_0;

        gs_x_xyy_xz[i] = 2.0 * ts_yy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xz[i] * gc_x[i] * tce_0;

        gs_x_xyy_yy[i] = 2.0 * ts_yy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yy[i] * gc_x[i] * tce_0;

        gs_x_xyy_yz[i] = 2.0 * ts_yy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yz[i] * gc_x[i] * tce_0;

        gs_x_xyy_zz[i] = 2.0 * ts_yy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 24-30 components of targeted buffer : FD

    auto gs_x_xyz_xx = pbuffer.data(idx_g_fd + 24);

    auto gs_x_xyz_xy = pbuffer.data(idx_g_fd + 25);

    auto gs_x_xyz_xz = pbuffer.data(idx_g_fd + 26);

    auto gs_x_xyz_yy = pbuffer.data(idx_g_fd + 27);

    auto gs_x_xyz_yz = pbuffer.data(idx_g_fd + 28);

    auto gs_x_xyz_zz = pbuffer.data(idx_g_fd + 29);

    #pragma omp simd aligned(gc_x, gs_x_xyz_xx, gs_x_xyz_xy, gs_x_xyz_xz, gs_x_xyz_yy, gs_x_xyz_yz, gs_x_xyz_zz, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_yy, ts_yz_yz, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyz_xx[i] = 2.0 * ts_yz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xx[i] * gc_x[i] * tce_0;

        gs_x_xyz_xy[i] = 2.0 * ts_yz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xy[i] * gc_x[i] * tce_0;

        gs_x_xyz_xz[i] = 2.0 * ts_yz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xz[i] * gc_x[i] * tce_0;

        gs_x_xyz_yy[i] = 2.0 * ts_yz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yy[i] * gc_x[i] * tce_0;

        gs_x_xyz_yz[i] = 2.0 * ts_yz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yz[i] * gc_x[i] * tce_0;

        gs_x_xyz_zz[i] = 2.0 * ts_yz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-36 components of targeted buffer : FD

    auto gs_x_xzz_xx = pbuffer.data(idx_g_fd + 30);

    auto gs_x_xzz_xy = pbuffer.data(idx_g_fd + 31);

    auto gs_x_xzz_xz = pbuffer.data(idx_g_fd + 32);

    auto gs_x_xzz_yy = pbuffer.data(idx_g_fd + 33);

    auto gs_x_xzz_yz = pbuffer.data(idx_g_fd + 34);

    auto gs_x_xzz_zz = pbuffer.data(idx_g_fd + 35);

    #pragma omp simd aligned(gc_x, gs_x_xzz_xx, gs_x_xzz_xy, gs_x_xzz_xz, gs_x_xzz_yy, gs_x_xzz_yz, gs_x_xzz_zz, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_yy, ts_zz_yz, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzz_xx[i] = 2.0 * ts_zz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xzz_xy[i] = 2.0 * ts_zz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xzz_xz[i] = 2.0 * ts_zz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xzz_yy[i] = 2.0 * ts_zz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xzz_yz[i] = 2.0 * ts_zz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xzz_zz[i] = 2.0 * ts_zz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 36-42 components of targeted buffer : FD

    auto gs_x_yyy_xx = pbuffer.data(idx_g_fd + 36);

    auto gs_x_yyy_xy = pbuffer.data(idx_g_fd + 37);

    auto gs_x_yyy_xz = pbuffer.data(idx_g_fd + 38);

    auto gs_x_yyy_yy = pbuffer.data(idx_g_fd + 39);

    auto gs_x_yyy_yz = pbuffer.data(idx_g_fd + 40);

    auto gs_x_yyy_zz = pbuffer.data(idx_g_fd + 41);

    #pragma omp simd aligned(gc_x, gs_x_yyy_xx, gs_x_yyy_xy, gs_x_yyy_xz, gs_x_yyy_yy, gs_x_yyy_yz, gs_x_yyy_zz, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyy_xx[i] = 4.0 * ts_yyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xx[i] * gc_x[i] * tce_0;

        gs_x_yyy_xy[i] = 2.0 * ts_yyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xy[i] * gc_x[i] * tce_0;

        gs_x_yyy_xz[i] = 2.0 * ts_yyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xz[i] * gc_x[i] * tce_0;

        gs_x_yyy_yy[i] = 2.0 * ts_yyy_yy[i] * gc_x[i] * tce_0;

        gs_x_yyy_yz[i] = 2.0 * ts_yyy_yz[i] * gc_x[i] * tce_0;

        gs_x_yyy_zz[i] = 2.0 * ts_yyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-48 components of targeted buffer : FD

    auto gs_x_yyz_xx = pbuffer.data(idx_g_fd + 42);

    auto gs_x_yyz_xy = pbuffer.data(idx_g_fd + 43);

    auto gs_x_yyz_xz = pbuffer.data(idx_g_fd + 44);

    auto gs_x_yyz_yy = pbuffer.data(idx_g_fd + 45);

    auto gs_x_yyz_yz = pbuffer.data(idx_g_fd + 46);

    auto gs_x_yyz_zz = pbuffer.data(idx_g_fd + 47);

    #pragma omp simd aligned(gc_x, gs_x_yyz_xx, gs_x_yyz_xy, gs_x_yyz_xz, gs_x_yyz_yy, gs_x_yyz_yz, gs_x_yyz_zz, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyz_xx[i] = 4.0 * ts_yyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xx[i] * gc_x[i] * tce_0;

        gs_x_yyz_xy[i] = 2.0 * ts_yyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xy[i] * gc_x[i] * tce_0;

        gs_x_yyz_xz[i] = 2.0 * ts_yyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xz[i] * gc_x[i] * tce_0;

        gs_x_yyz_yy[i] = 2.0 * ts_yyz_yy[i] * gc_x[i] * tce_0;

        gs_x_yyz_yz[i] = 2.0 * ts_yyz_yz[i] * gc_x[i] * tce_0;

        gs_x_yyz_zz[i] = 2.0 * ts_yyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 48-54 components of targeted buffer : FD

    auto gs_x_yzz_xx = pbuffer.data(idx_g_fd + 48);

    auto gs_x_yzz_xy = pbuffer.data(idx_g_fd + 49);

    auto gs_x_yzz_xz = pbuffer.data(idx_g_fd + 50);

    auto gs_x_yzz_yy = pbuffer.data(idx_g_fd + 51);

    auto gs_x_yzz_yz = pbuffer.data(idx_g_fd + 52);

    auto gs_x_yzz_zz = pbuffer.data(idx_g_fd + 53);

    #pragma omp simd aligned(gc_x, gs_x_yzz_xx, gs_x_yzz_xy, gs_x_yzz_xz, gs_x_yzz_yy, gs_x_yzz_yz, gs_x_yzz_zz, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzz_xx[i] = 4.0 * ts_yzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xx[i] * gc_x[i] * tce_0;

        gs_x_yzz_xy[i] = 2.0 * ts_yzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xy[i] * gc_x[i] * tce_0;

        gs_x_yzz_xz[i] = 2.0 * ts_yzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xz[i] * gc_x[i] * tce_0;

        gs_x_yzz_yy[i] = 2.0 * ts_yzz_yy[i] * gc_x[i] * tce_0;

        gs_x_yzz_yz[i] = 2.0 * ts_yzz_yz[i] * gc_x[i] * tce_0;

        gs_x_yzz_zz[i] = 2.0 * ts_yzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 54-60 components of targeted buffer : FD

    auto gs_x_zzz_xx = pbuffer.data(idx_g_fd + 54);

    auto gs_x_zzz_xy = pbuffer.data(idx_g_fd + 55);

    auto gs_x_zzz_xz = pbuffer.data(idx_g_fd + 56);

    auto gs_x_zzz_yy = pbuffer.data(idx_g_fd + 57);

    auto gs_x_zzz_yz = pbuffer.data(idx_g_fd + 58);

    auto gs_x_zzz_zz = pbuffer.data(idx_g_fd + 59);

    #pragma omp simd aligned(gc_x, gs_x_zzz_xx, gs_x_zzz_xy, gs_x_zzz_xz, gs_x_zzz_yy, gs_x_zzz_yz, gs_x_zzz_zz, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzz_xx[i] = 4.0 * ts_zzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xx[i] * gc_x[i] * tce_0;

        gs_x_zzz_xy[i] = 2.0 * ts_zzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xy[i] * gc_x[i] * tce_0;

        gs_x_zzz_xz[i] = 2.0 * ts_zzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xz[i] * gc_x[i] * tce_0;

        gs_x_zzz_yy[i] = 2.0 * ts_zzz_yy[i] * gc_x[i] * tce_0;

        gs_x_zzz_yz[i] = 2.0 * ts_zzz_yz[i] * gc_x[i] * tce_0;

        gs_x_zzz_zz[i] = 2.0 * ts_zzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-66 components of targeted buffer : FD

    auto gs_y_xxx_xx = pbuffer.data(idx_g_fd + 60);

    auto gs_y_xxx_xy = pbuffer.data(idx_g_fd + 61);

    auto gs_y_xxx_xz = pbuffer.data(idx_g_fd + 62);

    auto gs_y_xxx_yy = pbuffer.data(idx_g_fd + 63);

    auto gs_y_xxx_yz = pbuffer.data(idx_g_fd + 64);

    auto gs_y_xxx_zz = pbuffer.data(idx_g_fd + 65);

    #pragma omp simd aligned(gc_y, gs_y_xxx_xx, gs_y_xxx_xy, gs_y_xxx_xz, gs_y_xxx_yy, gs_y_xxx_yz, gs_y_xxx_zz, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxx_xx[i] = 2.0 * ts_xxx_xx[i] * gc_y[i] * tce_0;

        gs_y_xxx_xy[i] = 2.0 * ts_xxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xy[i] * gc_y[i] * tce_0;

        gs_y_xxx_xz[i] = 2.0 * ts_xxx_xz[i] * gc_y[i] * tce_0;

        gs_y_xxx_yy[i] = 4.0 * ts_xxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yy[i] * gc_y[i] * tce_0;

        gs_y_xxx_yz[i] = 2.0 * ts_xxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yz[i] * gc_y[i] * tce_0;

        gs_y_xxx_zz[i] = 2.0 * ts_xxx_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 66-72 components of targeted buffer : FD

    auto gs_y_xxy_xx = pbuffer.data(idx_g_fd + 66);

    auto gs_y_xxy_xy = pbuffer.data(idx_g_fd + 67);

    auto gs_y_xxy_xz = pbuffer.data(idx_g_fd + 68);

    auto gs_y_xxy_yy = pbuffer.data(idx_g_fd + 69);

    auto gs_y_xxy_yz = pbuffer.data(idx_g_fd + 70);

    auto gs_y_xxy_zz = pbuffer.data(idx_g_fd + 71);

    #pragma omp simd aligned(gc_y, gs_y_xxy_xx, gs_y_xxy_xy, gs_y_xxy_xz, gs_y_xxy_yy, gs_y_xxy_yz, gs_y_xxy_zz, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_yy, ts_xx_yz, ts_xx_zz, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxy_xx[i] = 2.0 * ts_xx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xx[i] * gc_y[i] * tce_0;

        gs_y_xxy_xy[i] = 2.0 * ts_xx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xy[i] * gc_y[i] * tce_0;

        gs_y_xxy_xz[i] = 2.0 * ts_xx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xz[i] * gc_y[i] * tce_0;

        gs_y_xxy_yy[i] = 2.0 * ts_xx_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yy[i] * gc_y[i] * tce_0;

        gs_y_xxy_yz[i] = 2.0 * ts_xx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yz[i] * gc_y[i] * tce_0;

        gs_y_xxy_zz[i] = 2.0 * ts_xx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 72-78 components of targeted buffer : FD

    auto gs_y_xxz_xx = pbuffer.data(idx_g_fd + 72);

    auto gs_y_xxz_xy = pbuffer.data(idx_g_fd + 73);

    auto gs_y_xxz_xz = pbuffer.data(idx_g_fd + 74);

    auto gs_y_xxz_yy = pbuffer.data(idx_g_fd + 75);

    auto gs_y_xxz_yz = pbuffer.data(idx_g_fd + 76);

    auto gs_y_xxz_zz = pbuffer.data(idx_g_fd + 77);

    #pragma omp simd aligned(gc_y, gs_y_xxz_xx, gs_y_xxz_xy, gs_y_xxz_xz, gs_y_xxz_yy, gs_y_xxz_yz, gs_y_xxz_zz, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxz_xx[i] = 2.0 * ts_xxz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxz_xy[i] = 2.0 * ts_xxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxz_xz[i] = 2.0 * ts_xxz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxz_yy[i] = 4.0 * ts_xxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxz_yz[i] = 2.0 * ts_xxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxz_zz[i] = 2.0 * ts_xxz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 78-84 components of targeted buffer : FD

    auto gs_y_xyy_xx = pbuffer.data(idx_g_fd + 78);

    auto gs_y_xyy_xy = pbuffer.data(idx_g_fd + 79);

    auto gs_y_xyy_xz = pbuffer.data(idx_g_fd + 80);

    auto gs_y_xyy_yy = pbuffer.data(idx_g_fd + 81);

    auto gs_y_xyy_yz = pbuffer.data(idx_g_fd + 82);

    auto gs_y_xyy_zz = pbuffer.data(idx_g_fd + 83);

    #pragma omp simd aligned(gc_y, gs_y_xyy_xx, gs_y_xyy_xy, gs_y_xyy_xz, gs_y_xyy_yy, gs_y_xyy_yz, gs_y_xyy_zz, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_yy, ts_xy_yz, ts_xy_zz, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyy_xx[i] = 4.0 * ts_xy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xx[i] * gc_y[i] * tce_0;

        gs_y_xyy_xy[i] = 4.0 * ts_xy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xy[i] * gc_y[i] * tce_0;

        gs_y_xyy_xz[i] = 4.0 * ts_xy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xz[i] * gc_y[i] * tce_0;

        gs_y_xyy_yy[i] = 4.0 * ts_xy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yy[i] * gc_y[i] * tce_0;

        gs_y_xyy_yz[i] = 4.0 * ts_xy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yz[i] * gc_y[i] * tce_0;

        gs_y_xyy_zz[i] = 4.0 * ts_xy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 84-90 components of targeted buffer : FD

    auto gs_y_xyz_xx = pbuffer.data(idx_g_fd + 84);

    auto gs_y_xyz_xy = pbuffer.data(idx_g_fd + 85);

    auto gs_y_xyz_xz = pbuffer.data(idx_g_fd + 86);

    auto gs_y_xyz_yy = pbuffer.data(idx_g_fd + 87);

    auto gs_y_xyz_yz = pbuffer.data(idx_g_fd + 88);

    auto gs_y_xyz_zz = pbuffer.data(idx_g_fd + 89);

    #pragma omp simd aligned(gc_y, gs_y_xyz_xx, gs_y_xyz_xy, gs_y_xyz_xz, gs_y_xyz_yy, gs_y_xyz_yz, gs_y_xyz_zz, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_yy, ts_xz_yz, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyz_xx[i] = 2.0 * ts_xz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xx[i] * gc_y[i] * tce_0;

        gs_y_xyz_xy[i] = 2.0 * ts_xz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xy[i] * gc_y[i] * tce_0;

        gs_y_xyz_xz[i] = 2.0 * ts_xz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xz[i] * gc_y[i] * tce_0;

        gs_y_xyz_yy[i] = 2.0 * ts_xz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yy[i] * gc_y[i] * tce_0;

        gs_y_xyz_yz[i] = 2.0 * ts_xz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yz[i] * gc_y[i] * tce_0;

        gs_y_xyz_zz[i] = 2.0 * ts_xz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 90-96 components of targeted buffer : FD

    auto gs_y_xzz_xx = pbuffer.data(idx_g_fd + 90);

    auto gs_y_xzz_xy = pbuffer.data(idx_g_fd + 91);

    auto gs_y_xzz_xz = pbuffer.data(idx_g_fd + 92);

    auto gs_y_xzz_yy = pbuffer.data(idx_g_fd + 93);

    auto gs_y_xzz_yz = pbuffer.data(idx_g_fd + 94);

    auto gs_y_xzz_zz = pbuffer.data(idx_g_fd + 95);

    #pragma omp simd aligned(gc_y, gs_y_xzz_xx, gs_y_xzz_xy, gs_y_xzz_xz, gs_y_xzz_yy, gs_y_xzz_yz, gs_y_xzz_zz, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzz_xx[i] = 2.0 * ts_xzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xzz_xy[i] = 2.0 * ts_xzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xzz_xz[i] = 2.0 * ts_xzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xzz_yy[i] = 4.0 * ts_xzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xzz_yz[i] = 2.0 * ts_xzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xzz_zz[i] = 2.0 * ts_xzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 96-102 components of targeted buffer : FD

    auto gs_y_yyy_xx = pbuffer.data(idx_g_fd + 96);

    auto gs_y_yyy_xy = pbuffer.data(idx_g_fd + 97);

    auto gs_y_yyy_xz = pbuffer.data(idx_g_fd + 98);

    auto gs_y_yyy_yy = pbuffer.data(idx_g_fd + 99);

    auto gs_y_yyy_yz = pbuffer.data(idx_g_fd + 100);

    auto gs_y_yyy_zz = pbuffer.data(idx_g_fd + 101);

    #pragma omp simd aligned(gc_y, gs_y_yyy_xx, gs_y_yyy_xy, gs_y_yyy_xz, gs_y_yyy_yy, gs_y_yyy_yz, gs_y_yyy_zz, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_yy, ts_yy_yz, ts_yy_zz, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyy_xx[i] = 6.0 * ts_yy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xx[i] * gc_y[i] * tce_0;

        gs_y_yyy_xy[i] = 6.0 * ts_yy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xy[i] * gc_y[i] * tce_0;

        gs_y_yyy_xz[i] = 6.0 * ts_yy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xz[i] * gc_y[i] * tce_0;

        gs_y_yyy_yy[i] = 6.0 * ts_yy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yy[i] * gc_y[i] * tce_0;

        gs_y_yyy_yz[i] = 6.0 * ts_yy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yz[i] * gc_y[i] * tce_0;

        gs_y_yyy_zz[i] = 6.0 * ts_yy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 102-108 components of targeted buffer : FD

    auto gs_y_yyz_xx = pbuffer.data(idx_g_fd + 102);

    auto gs_y_yyz_xy = pbuffer.data(idx_g_fd + 103);

    auto gs_y_yyz_xz = pbuffer.data(idx_g_fd + 104);

    auto gs_y_yyz_yy = pbuffer.data(idx_g_fd + 105);

    auto gs_y_yyz_yz = pbuffer.data(idx_g_fd + 106);

    auto gs_y_yyz_zz = pbuffer.data(idx_g_fd + 107);

    #pragma omp simd aligned(gc_y, gs_y_yyz_xx, gs_y_yyz_xy, gs_y_yyz_xz, gs_y_yyz_yy, gs_y_yyz_yz, gs_y_yyz_zz, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_yy, ts_yz_yz, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyz_xx[i] = 4.0 * ts_yz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xx[i] * gc_y[i] * tce_0;

        gs_y_yyz_xy[i] = 4.0 * ts_yz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xy[i] * gc_y[i] * tce_0;

        gs_y_yyz_xz[i] = 4.0 * ts_yz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xz[i] * gc_y[i] * tce_0;

        gs_y_yyz_yy[i] = 4.0 * ts_yz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yy[i] * gc_y[i] * tce_0;

        gs_y_yyz_yz[i] = 4.0 * ts_yz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yz[i] * gc_y[i] * tce_0;

        gs_y_yyz_zz[i] = 4.0 * ts_yz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 108-114 components of targeted buffer : FD

    auto gs_y_yzz_xx = pbuffer.data(idx_g_fd + 108);

    auto gs_y_yzz_xy = pbuffer.data(idx_g_fd + 109);

    auto gs_y_yzz_xz = pbuffer.data(idx_g_fd + 110);

    auto gs_y_yzz_yy = pbuffer.data(idx_g_fd + 111);

    auto gs_y_yzz_yz = pbuffer.data(idx_g_fd + 112);

    auto gs_y_yzz_zz = pbuffer.data(idx_g_fd + 113);

    #pragma omp simd aligned(gc_y, gs_y_yzz_xx, gs_y_yzz_xy, gs_y_yzz_xz, gs_y_yzz_yy, gs_y_yzz_yz, gs_y_yzz_zz, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_yy, ts_zz_yz, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzz_xx[i] = 2.0 * ts_zz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xx[i] * gc_y[i] * tce_0;

        gs_y_yzz_xy[i] = 2.0 * ts_zz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xy[i] * gc_y[i] * tce_0;

        gs_y_yzz_xz[i] = 2.0 * ts_zz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xz[i] * gc_y[i] * tce_0;

        gs_y_yzz_yy[i] = 2.0 * ts_zz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yy[i] * gc_y[i] * tce_0;

        gs_y_yzz_yz[i] = 2.0 * ts_zz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yz[i] * gc_y[i] * tce_0;

        gs_y_yzz_zz[i] = 2.0 * ts_zz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 114-120 components of targeted buffer : FD

    auto gs_y_zzz_xx = pbuffer.data(idx_g_fd + 114);

    auto gs_y_zzz_xy = pbuffer.data(idx_g_fd + 115);

    auto gs_y_zzz_xz = pbuffer.data(idx_g_fd + 116);

    auto gs_y_zzz_yy = pbuffer.data(idx_g_fd + 117);

    auto gs_y_zzz_yz = pbuffer.data(idx_g_fd + 118);

    auto gs_y_zzz_zz = pbuffer.data(idx_g_fd + 119);

    #pragma omp simd aligned(gc_y, gs_y_zzz_xx, gs_y_zzz_xy, gs_y_zzz_xz, gs_y_zzz_yy, gs_y_zzz_yz, gs_y_zzz_zz, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzz_xx[i] = 2.0 * ts_zzz_xx[i] * gc_y[i] * tce_0;

        gs_y_zzz_xy[i] = 2.0 * ts_zzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xy[i] * gc_y[i] * tce_0;

        gs_y_zzz_xz[i] = 2.0 * ts_zzz_xz[i] * gc_y[i] * tce_0;

        gs_y_zzz_yy[i] = 4.0 * ts_zzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yy[i] * gc_y[i] * tce_0;

        gs_y_zzz_yz[i] = 2.0 * ts_zzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yz[i] * gc_y[i] * tce_0;

        gs_y_zzz_zz[i] = 2.0 * ts_zzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 120-126 components of targeted buffer : FD

    auto gs_z_xxx_xx = pbuffer.data(idx_g_fd + 120);

    auto gs_z_xxx_xy = pbuffer.data(idx_g_fd + 121);

    auto gs_z_xxx_xz = pbuffer.data(idx_g_fd + 122);

    auto gs_z_xxx_yy = pbuffer.data(idx_g_fd + 123);

    auto gs_z_xxx_yz = pbuffer.data(idx_g_fd + 124);

    auto gs_z_xxx_zz = pbuffer.data(idx_g_fd + 125);

    #pragma omp simd aligned(gc_z, gs_z_xxx_xx, gs_z_xxx_xy, gs_z_xxx_xz, gs_z_xxx_yy, gs_z_xxx_yz, gs_z_xxx_zz, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxx_xx[i] = 2.0 * ts_xxx_xx[i] * gc_z[i] * tce_0;

        gs_z_xxx_xy[i] = 2.0 * ts_xxx_xy[i] * gc_z[i] * tce_0;

        gs_z_xxx_xz[i] = 2.0 * ts_xxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_xz[i] * gc_z[i] * tce_0;

        gs_z_xxx_yy[i] = 2.0 * ts_xxx_yy[i] * gc_z[i] * tce_0;

        gs_z_xxx_yz[i] = 2.0 * ts_xxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_yz[i] * gc_z[i] * tce_0;

        gs_z_xxx_zz[i] = 4.0 * ts_xxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 126-132 components of targeted buffer : FD

    auto gs_z_xxy_xx = pbuffer.data(idx_g_fd + 126);

    auto gs_z_xxy_xy = pbuffer.data(idx_g_fd + 127);

    auto gs_z_xxy_xz = pbuffer.data(idx_g_fd + 128);

    auto gs_z_xxy_yy = pbuffer.data(idx_g_fd + 129);

    auto gs_z_xxy_yz = pbuffer.data(idx_g_fd + 130);

    auto gs_z_xxy_zz = pbuffer.data(idx_g_fd + 131);

    #pragma omp simd aligned(gc_z, gs_z_xxy_xx, gs_z_xxy_xy, gs_z_xxy_xz, gs_z_xxy_yy, gs_z_xxy_yz, gs_z_xxy_zz, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxy_xx[i] = 2.0 * ts_xxy_xx[i] * gc_z[i] * tce_0;

        gs_z_xxy_xy[i] = 2.0 * ts_xxy_xy[i] * gc_z[i] * tce_0;

        gs_z_xxy_xz[i] = 2.0 * ts_xxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_xz[i] * gc_z[i] * tce_0;

        gs_z_xxy_yy[i] = 2.0 * ts_xxy_yy[i] * gc_z[i] * tce_0;

        gs_z_xxy_yz[i] = 2.0 * ts_xxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_yz[i] * gc_z[i] * tce_0;

        gs_z_xxy_zz[i] = 4.0 * ts_xxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 132-138 components of targeted buffer : FD

    auto gs_z_xxz_xx = pbuffer.data(idx_g_fd + 132);

    auto gs_z_xxz_xy = pbuffer.data(idx_g_fd + 133);

    auto gs_z_xxz_xz = pbuffer.data(idx_g_fd + 134);

    auto gs_z_xxz_yy = pbuffer.data(idx_g_fd + 135);

    auto gs_z_xxz_yz = pbuffer.data(idx_g_fd + 136);

    auto gs_z_xxz_zz = pbuffer.data(idx_g_fd + 137);

    #pragma omp simd aligned(gc_z, gs_z_xxz_xx, gs_z_xxz_xy, gs_z_xxz_xz, gs_z_xxz_yy, gs_z_xxz_yz, gs_z_xxz_zz, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_yy, ts_xx_yz, ts_xx_zz, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxz_xx[i] = 2.0 * ts_xx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxz_xy[i] = 2.0 * ts_xx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxz_xz[i] = 2.0 * ts_xx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxz_yy[i] = 2.0 * ts_xx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxz_yz[i] = 2.0 * ts_xx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxz_zz[i] = 2.0 * ts_xx_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 138-144 components of targeted buffer : FD

    auto gs_z_xyy_xx = pbuffer.data(idx_g_fd + 138);

    auto gs_z_xyy_xy = pbuffer.data(idx_g_fd + 139);

    auto gs_z_xyy_xz = pbuffer.data(idx_g_fd + 140);

    auto gs_z_xyy_yy = pbuffer.data(idx_g_fd + 141);

    auto gs_z_xyy_yz = pbuffer.data(idx_g_fd + 142);

    auto gs_z_xyy_zz = pbuffer.data(idx_g_fd + 143);

    #pragma omp simd aligned(gc_z, gs_z_xyy_xx, gs_z_xyy_xy, gs_z_xyy_xz, gs_z_xyy_yy, gs_z_xyy_yz, gs_z_xyy_zz, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyy_xx[i] = 2.0 * ts_xyy_xx[i] * gc_z[i] * tce_0;

        gs_z_xyy_xy[i] = 2.0 * ts_xyy_xy[i] * gc_z[i] * tce_0;

        gs_z_xyy_xz[i] = 2.0 * ts_xyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_xz[i] * gc_z[i] * tce_0;

        gs_z_xyy_yy[i] = 2.0 * ts_xyy_yy[i] * gc_z[i] * tce_0;

        gs_z_xyy_yz[i] = 2.0 * ts_xyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_yz[i] * gc_z[i] * tce_0;

        gs_z_xyy_zz[i] = 4.0 * ts_xyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 144-150 components of targeted buffer : FD

    auto gs_z_xyz_xx = pbuffer.data(idx_g_fd + 144);

    auto gs_z_xyz_xy = pbuffer.data(idx_g_fd + 145);

    auto gs_z_xyz_xz = pbuffer.data(idx_g_fd + 146);

    auto gs_z_xyz_yy = pbuffer.data(idx_g_fd + 147);

    auto gs_z_xyz_yz = pbuffer.data(idx_g_fd + 148);

    auto gs_z_xyz_zz = pbuffer.data(idx_g_fd + 149);

    #pragma omp simd aligned(gc_z, gs_z_xyz_xx, gs_z_xyz_xy, gs_z_xyz_xz, gs_z_xyz_yy, gs_z_xyz_yz, gs_z_xyz_zz, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_yy, ts_xy_yz, ts_xy_zz, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyz_xx[i] = 2.0 * ts_xy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xx[i] * gc_z[i] * tce_0;

        gs_z_xyz_xy[i] = 2.0 * ts_xy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xy[i] * gc_z[i] * tce_0;

        gs_z_xyz_xz[i] = 2.0 * ts_xy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_xz[i] * gc_z[i] * tce_0;

        gs_z_xyz_yy[i] = 2.0 * ts_xy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yy[i] * gc_z[i] * tce_0;

        gs_z_xyz_yz[i] = 2.0 * ts_xy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_yz[i] * gc_z[i] * tce_0;

        gs_z_xyz_zz[i] = 2.0 * ts_xy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 150-156 components of targeted buffer : FD

    auto gs_z_xzz_xx = pbuffer.data(idx_g_fd + 150);

    auto gs_z_xzz_xy = pbuffer.data(idx_g_fd + 151);

    auto gs_z_xzz_xz = pbuffer.data(idx_g_fd + 152);

    auto gs_z_xzz_yy = pbuffer.data(idx_g_fd + 153);

    auto gs_z_xzz_yz = pbuffer.data(idx_g_fd + 154);

    auto gs_z_xzz_zz = pbuffer.data(idx_g_fd + 155);

    #pragma omp simd aligned(gc_z, gs_z_xzz_xx, gs_z_xzz_xy, gs_z_xzz_xz, gs_z_xzz_yy, gs_z_xzz_yz, gs_z_xzz_zz, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_yy, ts_xz_yz, ts_xz_zz, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzz_xx[i] = 4.0 * ts_xz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xzz_xy[i] = 4.0 * ts_xz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xzz_xz[i] = 4.0 * ts_xz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xzz_yy[i] = 4.0 * ts_xz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xzz_yz[i] = 4.0 * ts_xz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xzz_zz[i] = 4.0 * ts_xz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 156-162 components of targeted buffer : FD

    auto gs_z_yyy_xx = pbuffer.data(idx_g_fd + 156);

    auto gs_z_yyy_xy = pbuffer.data(idx_g_fd + 157);

    auto gs_z_yyy_xz = pbuffer.data(idx_g_fd + 158);

    auto gs_z_yyy_yy = pbuffer.data(idx_g_fd + 159);

    auto gs_z_yyy_yz = pbuffer.data(idx_g_fd + 160);

    auto gs_z_yyy_zz = pbuffer.data(idx_g_fd + 161);

    #pragma omp simd aligned(gc_z, gs_z_yyy_xx, gs_z_yyy_xy, gs_z_yyy_xz, gs_z_yyy_yy, gs_z_yyy_yz, gs_z_yyy_zz, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyy_xx[i] = 2.0 * ts_yyy_xx[i] * gc_z[i] * tce_0;

        gs_z_yyy_xy[i] = 2.0 * ts_yyy_xy[i] * gc_z[i] * tce_0;

        gs_z_yyy_xz[i] = 2.0 * ts_yyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_xz[i] * gc_z[i] * tce_0;

        gs_z_yyy_yy[i] = 2.0 * ts_yyy_yy[i] * gc_z[i] * tce_0;

        gs_z_yyy_yz[i] = 2.0 * ts_yyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_yz[i] * gc_z[i] * tce_0;

        gs_z_yyy_zz[i] = 4.0 * ts_yyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 162-168 components of targeted buffer : FD

    auto gs_z_yyz_xx = pbuffer.data(idx_g_fd + 162);

    auto gs_z_yyz_xy = pbuffer.data(idx_g_fd + 163);

    auto gs_z_yyz_xz = pbuffer.data(idx_g_fd + 164);

    auto gs_z_yyz_yy = pbuffer.data(idx_g_fd + 165);

    auto gs_z_yyz_yz = pbuffer.data(idx_g_fd + 166);

    auto gs_z_yyz_zz = pbuffer.data(idx_g_fd + 167);

    #pragma omp simd aligned(gc_z, gs_z_yyz_xx, gs_z_yyz_xy, gs_z_yyz_xz, gs_z_yyz_yy, gs_z_yyz_yz, gs_z_yyz_zz, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_yy, ts_yy_yz, ts_yy_zz, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyz_xx[i] = 2.0 * ts_yy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xx[i] * gc_z[i] * tce_0;

        gs_z_yyz_xy[i] = 2.0 * ts_yy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xy[i] * gc_z[i] * tce_0;

        gs_z_yyz_xz[i] = 2.0 * ts_yy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_xz[i] * gc_z[i] * tce_0;

        gs_z_yyz_yy[i] = 2.0 * ts_yy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yy[i] * gc_z[i] * tce_0;

        gs_z_yyz_yz[i] = 2.0 * ts_yy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_yz[i] * gc_z[i] * tce_0;

        gs_z_yyz_zz[i] = 2.0 * ts_yy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 168-174 components of targeted buffer : FD

    auto gs_z_yzz_xx = pbuffer.data(idx_g_fd + 168);

    auto gs_z_yzz_xy = pbuffer.data(idx_g_fd + 169);

    auto gs_z_yzz_xz = pbuffer.data(idx_g_fd + 170);

    auto gs_z_yzz_yy = pbuffer.data(idx_g_fd + 171);

    auto gs_z_yzz_yz = pbuffer.data(idx_g_fd + 172);

    auto gs_z_yzz_zz = pbuffer.data(idx_g_fd + 173);

    #pragma omp simd aligned(gc_z, gs_z_yzz_xx, gs_z_yzz_xy, gs_z_yzz_xz, gs_z_yzz_yy, gs_z_yzz_yz, gs_z_yzz_zz, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_yy, ts_yz_yz, ts_yz_zz, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzz_xx[i] = 4.0 * ts_yz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xx[i] * gc_z[i] * tce_0;

        gs_z_yzz_xy[i] = 4.0 * ts_yz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xy[i] * gc_z[i] * tce_0;

        gs_z_yzz_xz[i] = 4.0 * ts_yz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_xz[i] * gc_z[i] * tce_0;

        gs_z_yzz_yy[i] = 4.0 * ts_yz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yy[i] * gc_z[i] * tce_0;

        gs_z_yzz_yz[i] = 4.0 * ts_yz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_yz[i] * gc_z[i] * tce_0;

        gs_z_yzz_zz[i] = 4.0 * ts_yz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 174-180 components of targeted buffer : FD

    auto gs_z_zzz_xx = pbuffer.data(idx_g_fd + 174);

    auto gs_z_zzz_xy = pbuffer.data(idx_g_fd + 175);

    auto gs_z_zzz_xz = pbuffer.data(idx_g_fd + 176);

    auto gs_z_zzz_yy = pbuffer.data(idx_g_fd + 177);

    auto gs_z_zzz_yz = pbuffer.data(idx_g_fd + 178);

    auto gs_z_zzz_zz = pbuffer.data(idx_g_fd + 179);

    #pragma omp simd aligned(gc_z, gs_z_zzz_xx, gs_z_zzz_xy, gs_z_zzz_xz, gs_z_zzz_yy, gs_z_zzz_yz, gs_z_zzz_zz, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_yy, ts_zz_yz, ts_zz_zz, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzz_xx[i] = 6.0 * ts_zz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xx[i] * gc_z[i] * tce_0;

        gs_z_zzz_xy[i] = 6.0 * ts_zz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xy[i] * gc_z[i] * tce_0;

        gs_z_zzz_xz[i] = 6.0 * ts_zz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_xz[i] * gc_z[i] * tce_0;

        gs_z_zzz_yy[i] = 6.0 * ts_zz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yy[i] * gc_z[i] * tce_0;

        gs_z_zzz_yz[i] = 6.0 * ts_zz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_yz[i] * gc_z[i] * tce_0;

        gs_z_zzz_zz[i] = 6.0 * ts_zz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_zzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_zz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

