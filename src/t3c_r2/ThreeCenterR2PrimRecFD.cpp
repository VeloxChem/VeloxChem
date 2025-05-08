#include "ThreeCenterR2PrimRecFD.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_fd(CSimdArray<double>& pbuffer, 
                const size_t idx_g_fd,
                const size_t idx_pd,
                const size_t idx_dp,
                const size_t idx_dd,
                const size_t idx_fs,
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

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_pd);

    auto ts_x_xy = pbuffer.data(idx_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_pd + 17);

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

    auto gr_xxx_xx = pbuffer.data(idx_g_fd);

    auto gr_xxx_xy = pbuffer.data(idx_g_fd + 1);

    auto gr_xxx_xz = pbuffer.data(idx_g_fd + 2);

    auto gr_xxx_yy = pbuffer.data(idx_g_fd + 3);

    auto gr_xxx_yz = pbuffer.data(idx_g_fd + 4);

    auto gr_xxx_zz = pbuffer.data(idx_g_fd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xx, gr_xxx_xy, gr_xxx_xz, gr_xxx_yy, gr_xxx_yz, gr_xxx_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, ts_xxx_0, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxx_xx[i] = 6.0 * ts_x_xx[i] * gfe2_0 + 12.0 * ts_xx_x[i] * gfe2_0 + 6.0 * ts_xx_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe2_0 + 4.0 * ts_xxx_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxx_xx[i] * gfe_0 + ts_xxx_xx[i] * rgc2_0;

        gr_xxx_xy[i] = 6.0 * ts_x_xy[i] * gfe2_0 + 6.0 * ts_xx_y[i] * gfe2_0 + 6.0 * ts_xx_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_xy[i] * gfe_0 + ts_xxx_xy[i] * rgc2_0;

        gr_xxx_xz[i] = 6.0 * ts_x_xz[i] * gfe2_0 + 6.0 * ts_xx_z[i] * gfe2_0 + 6.0 * ts_xx_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_xz[i] * gfe_0 + ts_xxx_xz[i] * rgc2_0;

        gr_xxx_yy[i] = 6.0 * ts_x_yy[i] * gfe2_0 + 6.0 * ts_xx_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe2_0 + 4.0 * ts_xxx_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxx_yy[i] * gfe_0 + ts_xxx_yy[i] * rgc2_0;

        gr_xxx_yz[i] = 6.0 * ts_x_yz[i] * gfe2_0 + 6.0 * ts_xx_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxx_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_yz[i] * gfe_0 + ts_xxx_yz[i] * rgc2_0;

        gr_xxx_zz[i] = 6.0 * ts_x_zz[i] * gfe2_0 + 6.0 * ts_xx_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe2_0 + 4.0 * ts_xxx_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxx_zz[i] * gfe_0 + ts_xxx_zz[i] * rgc2_0;
    }

    // Set up 6-12 components of targeted buffer : FD

    auto gr_xxy_xx = pbuffer.data(idx_g_fd + 6);

    auto gr_xxy_xy = pbuffer.data(idx_g_fd + 7);

    auto gr_xxy_xz = pbuffer.data(idx_g_fd + 8);

    auto gr_xxy_yy = pbuffer.data(idx_g_fd + 9);

    auto gr_xxy_yz = pbuffer.data(idx_g_fd + 10);

    auto gr_xxy_zz = pbuffer.data(idx_g_fd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xx, gr_xxy_xy, gr_xxy_xz, gr_xxy_yy, gr_xxy_yz, gr_xxy_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, ts_xxy_0, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxy_xx[i] = 2.0 * ts_y_xx[i] * gfe2_0 + 8.0 * ts_xy_x[i] * gfe2_0 + 4.0 * ts_xy_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe2_0 + 4.0 * ts_xxy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxy_xx[i] * gfe_0 + ts_xxy_xx[i] * rgc2_0;

        gr_xxy_xy[i] = 2.0 * ts_y_xy[i] * gfe2_0 + 4.0 * ts_xy_y[i] * gfe2_0 + 4.0 * ts_xy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe2_0 + 2.0 * ts_xx_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_xy[i] * gfe_0 + ts_xxy_xy[i] * rgc2_0;

        gr_xxy_xz[i] = 2.0 * ts_y_xz[i] * gfe2_0 + 4.0 * ts_xy_z[i] * gfe2_0 + 4.0 * ts_xy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_xz[i] * gfe_0 + ts_xxy_xz[i] * rgc2_0;

        gr_xxy_yy[i] = 2.0 * ts_y_yy[i] * gfe2_0 + 4.0 * ts_xy_yy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_y[i] * gfe2_0 + 2.0 * ts_xx_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe2_0 + 4.0 * ts_xxy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_yy[i] * gfe_0 + ts_xxy_yy[i] * rgc2_0;

        gr_xxy_yz[i] = 2.0 * ts_y_yz[i] * gfe2_0 + 4.0 * ts_xy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe2_0 + 2.0 * ts_xx_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_yz[i] * gfe_0 + ts_xxy_yz[i] * rgc2_0;

        gr_xxy_zz[i] = 2.0 * ts_y_zz[i] * gfe2_0 + 4.0 * ts_xy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe2_0 + 4.0 * ts_xxy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxy_zz[i] * gfe_0 + ts_xxy_zz[i] * rgc2_0;
    }

    // Set up 12-18 components of targeted buffer : FD

    auto gr_xxz_xx = pbuffer.data(idx_g_fd + 12);

    auto gr_xxz_xy = pbuffer.data(idx_g_fd + 13);

    auto gr_xxz_xz = pbuffer.data(idx_g_fd + 14);

    auto gr_xxz_yy = pbuffer.data(idx_g_fd + 15);

    auto gr_xxz_yz = pbuffer.data(idx_g_fd + 16);

    auto gr_xxz_zz = pbuffer.data(idx_g_fd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xx, gr_xxz_xy, gr_xxz_xz, gr_xxz_yy, gr_xxz_yz, gr_xxz_zz, ts_xx_x, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_y, ts_xx_yy, ts_xx_yz, ts_xx_z, ts_xx_zz, ts_xxz_0, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxz_xx[i] = 2.0 * ts_z_xx[i] * gfe2_0 + 8.0 * ts_xz_x[i] * gfe2_0 + 4.0 * ts_xz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_0[i] * gfe2_0 + 4.0 * ts_xxz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxz_xx[i] * gfe_0 + ts_xxz_xx[i] * rgc2_0;

        gr_xxz_xy[i] = 2.0 * ts_z_xy[i] * gfe2_0 + 4.0 * ts_xz_y[i] * gfe2_0 + 4.0 * ts_xz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_xy[i] * gfe_0 + ts_xxz_xy[i] * rgc2_0;

        gr_xxz_xz[i] = 2.0 * ts_z_xz[i] * gfe2_0 + 4.0 * ts_xz_z[i] * gfe2_0 + 4.0 * ts_xz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe2_0 + 2.0 * ts_xx_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_xz[i] * gfe_0 + ts_xxz_xz[i] * rgc2_0;

        gr_xxz_yy[i] = 2.0 * ts_z_yy[i] * gfe2_0 + 4.0 * ts_xz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_0[i] * gfe2_0 + 4.0 * ts_xxz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxz_yy[i] * gfe_0 + ts_xxz_yy[i] * rgc2_0;

        gr_xxz_yz[i] = 2.0 * ts_z_yz[i] * gfe2_0 + 4.0 * ts_xz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_y[i] * gfe2_0 + 2.0 * ts_xx_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_yz[i] * gfe_0 + ts_xxz_yz[i] * rgc2_0;

        gr_xxz_zz[i] = 2.0 * ts_z_zz[i] * gfe2_0 + 4.0 * ts_xz_zz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xx_z[i] * gfe2_0 + 2.0 * ts_xx_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxz_0[i] * gfe2_0 + 4.0 * ts_xxz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_zz[i] * gfe_0 + ts_xxz_zz[i] * rgc2_0;
    }

    // Set up 18-24 components of targeted buffer : FD

    auto gr_xyy_xx = pbuffer.data(idx_g_fd + 18);

    auto gr_xyy_xy = pbuffer.data(idx_g_fd + 19);

    auto gr_xyy_xz = pbuffer.data(idx_g_fd + 20);

    auto gr_xyy_yy = pbuffer.data(idx_g_fd + 21);

    auto gr_xyy_yz = pbuffer.data(idx_g_fd + 22);

    auto gr_xyy_zz = pbuffer.data(idx_g_fd + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xx, gr_xyy_xy, gr_xyy_xz, gr_xyy_yy, gr_xyy_yz, gr_xyy_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, ts_xyy_0, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyy_xx[i] = 4.0 * ts_yy_x[i] * gfe2_0 + 2.0 * ts_yy_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe2_0 + 4.0 * ts_xy_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe2_0 + 4.0 * ts_xyy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyy_xx[i] * gfe_0 + ts_xyy_xx[i] * rgc2_0;

        gr_xyy_xy[i] = 2.0 * ts_yy_y[i] * gfe2_0 + 2.0 * ts_yy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xy[i] * gfe2_0 + 4.0 * ts_xy_x[i] * gfe2_0 + 4.0 * ts_xy_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_xy[i] * gfe_0 + ts_xyy_xy[i] * rgc2_0;

        gr_xyy_xz[i] = 2.0 * ts_yy_z[i] * gfe2_0 + 2.0 * ts_yy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xz[i] * gfe2_0 + 4.0 * ts_xy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_xz[i] * gfe_0 + ts_xyy_xz[i] * rgc2_0;

        gr_xyy_yy[i] = 2.0 * ts_yy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yy[i] * gfe2_0 + 8.0 * ts_xy_y[i] * gfe2_0 + 4.0 * ts_xy_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe2_0 + 4.0 * ts_xyy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_yy[i] * gfe_0 + ts_xyy_yy[i] * rgc2_0;

        gr_xyy_yz[i] = 2.0 * ts_yy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yz[i] * gfe2_0 + 4.0 * ts_xy_z[i] * gfe2_0 + 4.0 * ts_xy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_yz[i] * gfe_0 + ts_xyy_yz[i] * rgc2_0;

        gr_xyy_zz[i] = 2.0 * ts_yy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe2_0 + 4.0 * ts_xy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe2_0 + 4.0 * ts_xyy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyy_zz[i] * gfe_0 + ts_xyy_zz[i] * rgc2_0;
    }

    // Set up 24-30 components of targeted buffer : FD

    auto gr_xyz_xx = pbuffer.data(idx_g_fd + 24);

    auto gr_xyz_xy = pbuffer.data(idx_g_fd + 25);

    auto gr_xyz_xz = pbuffer.data(idx_g_fd + 26);

    auto gr_xyz_yy = pbuffer.data(idx_g_fd + 27);

    auto gr_xyz_yz = pbuffer.data(idx_g_fd + 28);

    auto gr_xyz_zz = pbuffer.data(idx_g_fd + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xx, gr_xyz_xy, gr_xyz_xz, gr_xyz_yy, gr_xyz_yz, gr_xyz_zz, ts_xy_x, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_y, ts_xy_yy, ts_xy_yz, ts_xy_z, ts_xy_zz, ts_xyz_0, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xyz_xx[i] = 4.0 * ts_yz_x[i] * gfe2_0 + 2.0 * ts_yz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_0[i] * gfe2_0 + 4.0 * ts_xyz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyz_xx[i] * gfe_0 + ts_xyz_xx[i] * rgc2_0;

        gr_xyz_xy[i] = 2.0 * ts_yz_y[i] * gfe2_0 + 2.0 * ts_yz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe2_0 + 2.0 * ts_xz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_xy[i] * gfe_0 + ts_xyz_xy[i] * rgc2_0;

        gr_xyz_xz[i] = 2.0 * ts_yz_z[i] * gfe2_0 + 2.0 * ts_yz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_x[i] * gfe2_0 + 2.0 * ts_xy_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_xz[i] * gfe_0 + ts_xyz_xz[i] * rgc2_0;

        gr_xyz_yy[i] = 2.0 * ts_yz_yy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xz_y[i] * gfe2_0 + 2.0 * ts_xz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_0[i] * gfe2_0 + 4.0 * ts_xyz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyz_yy[i] * gfe_0 + ts_xyz_yy[i] * rgc2_0;

        gr_xyz_yz[i] = 2.0 * ts_yz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_z[i] * gfe2_0 + 2.0 * ts_xz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe2_0 + 2.0 * ts_xy_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_yz[i] * gfe_0 + ts_xyz_yz[i] * rgc2_0;

        gr_xyz_zz[i] = 2.0 * ts_yz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xy_z[i] * gfe2_0 + 2.0 * ts_xy_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyz_0[i] * gfe2_0 + 4.0 * ts_xyz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_zz[i] * gfe_0 + ts_xyz_zz[i] * rgc2_0;
    }

    // Set up 30-36 components of targeted buffer : FD

    auto gr_xzz_xx = pbuffer.data(idx_g_fd + 30);

    auto gr_xzz_xy = pbuffer.data(idx_g_fd + 31);

    auto gr_xzz_xz = pbuffer.data(idx_g_fd + 32);

    auto gr_xzz_yy = pbuffer.data(idx_g_fd + 33);

    auto gr_xzz_yz = pbuffer.data(idx_g_fd + 34);

    auto gr_xzz_zz = pbuffer.data(idx_g_fd + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xx, gr_xzz_xy, gr_xzz_xz, gr_xzz_yy, gr_xzz_yz, gr_xzz_zz, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_yy, ts_x_yz, ts_x_zz, ts_xz_x, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_y, ts_xz_yy, ts_xz_yz, ts_xz_z, ts_xz_zz, ts_xzz_0, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xzz_xx[i] = 4.0 * ts_zz_x[i] * gfe2_0 + 2.0 * ts_zz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xx[i] * gfe2_0 + 4.0 * ts_xz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_0[i] * gfe2_0 + 4.0 * ts_xzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzz_xx[i] * gfe_0 + ts_xzz_xx[i] * rgc2_0;

        gr_xzz_xy[i] = 2.0 * ts_zz_y[i] * gfe2_0 + 2.0 * ts_zz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xy[i] * gfe2_0 + 4.0 * ts_xz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_xy[i] * gfe_0 + ts_xzz_xy[i] * rgc2_0;

        gr_xzz_xz[i] = 2.0 * ts_zz_z[i] * gfe2_0 + 2.0 * ts_zz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_xz[i] * gfe2_0 + 4.0 * ts_xz_x[i] * gfe2_0 + 4.0 * ts_xz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_xz[i] * gfe_0 + ts_xzz_xz[i] * rgc2_0;

        gr_xzz_yy[i] = 2.0 * ts_zz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yy[i] * gfe2_0 + 4.0 * ts_xz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_0[i] * gfe2_0 + 4.0 * ts_xzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzz_yy[i] * gfe_0 + ts_xzz_yy[i] * rgc2_0;

        gr_xzz_yz[i] = 2.0 * ts_zz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_yz[i] * gfe2_0 + 4.0 * ts_xz_y[i] * gfe2_0 + 4.0 * ts_xz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_yz[i] * gfe_0 + ts_xzz_yz[i] * rgc2_0;

        gr_xzz_zz[i] = 2.0 * ts_zz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_zz[i] * gfe2_0 + 8.0 * ts_xz_z[i] * gfe2_0 + 4.0 * ts_xz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzz_0[i] * gfe2_0 + 4.0 * ts_xzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_zz[i] * gfe_0 + ts_xzz_zz[i] * rgc2_0;
    }

    // Set up 36-42 components of targeted buffer : FD

    auto gr_yyy_xx = pbuffer.data(idx_g_fd + 36);

    auto gr_yyy_xy = pbuffer.data(idx_g_fd + 37);

    auto gr_yyy_xz = pbuffer.data(idx_g_fd + 38);

    auto gr_yyy_yy = pbuffer.data(idx_g_fd + 39);

    auto gr_yyy_yz = pbuffer.data(idx_g_fd + 40);

    auto gr_yyy_zz = pbuffer.data(idx_g_fd + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xx, gr_yyy_xy, gr_yyy_xz, gr_yyy_yy, gr_yyy_yz, gr_yyy_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, ts_yyy_0, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyy_xx[i] = 6.0 * ts_y_xx[i] * gfe2_0 + 6.0 * ts_yy_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe2_0 + 4.0 * ts_yyy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyy_xx[i] * gfe_0 + ts_yyy_xx[i] * rgc2_0;

        gr_yyy_xy[i] = 6.0 * ts_y_xy[i] * gfe2_0 + 6.0 * ts_yy_x[i] * gfe2_0 + 6.0 * ts_yy_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_xy[i] * gfe_0 + ts_yyy_xy[i] * rgc2_0;

        gr_yyy_xz[i] = 6.0 * ts_y_xz[i] * gfe2_0 + 6.0 * ts_yy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_xz[i] * gfe_0 + ts_yyy_xz[i] * rgc2_0;

        gr_yyy_yy[i] = 6.0 * ts_y_yy[i] * gfe2_0 + 12.0 * ts_yy_y[i] * gfe2_0 + 6.0 * ts_yy_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe2_0 + 4.0 * ts_yyy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_yy[i] * gfe_0 + ts_yyy_yy[i] * rgc2_0;

        gr_yyy_yz[i] = 6.0 * ts_y_yz[i] * gfe2_0 + 6.0 * ts_yy_z[i] * gfe2_0 + 6.0 * ts_yy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_yz[i] * gfe_0 + ts_yyy_yz[i] * rgc2_0;

        gr_yyy_zz[i] = 6.0 * ts_y_zz[i] * gfe2_0 + 6.0 * ts_yy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe2_0 + 4.0 * ts_yyy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyy_zz[i] * gfe_0 + ts_yyy_zz[i] * rgc2_0;
    }

    // Set up 42-48 components of targeted buffer : FD

    auto gr_yyz_xx = pbuffer.data(idx_g_fd + 42);

    auto gr_yyz_xy = pbuffer.data(idx_g_fd + 43);

    auto gr_yyz_xz = pbuffer.data(idx_g_fd + 44);

    auto gr_yyz_yy = pbuffer.data(idx_g_fd + 45);

    auto gr_yyz_yz = pbuffer.data(idx_g_fd + 46);

    auto gr_yyz_zz = pbuffer.data(idx_g_fd + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xx, gr_yyz_xy, gr_yyz_xz, gr_yyz_yy, gr_yyz_yz, gr_yyz_zz, ts_yy_x, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_y, ts_yy_yy, ts_yy_yz, ts_yy_z, ts_yy_zz, ts_yyz_0, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yyz_xx[i] = 2.0 * ts_z_xx[i] * gfe2_0 + 4.0 * ts_yz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_0[i] * gfe2_0 + 4.0 * ts_yyz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyz_xx[i] * gfe_0 + ts_yyz_xx[i] * rgc2_0;

        gr_yyz_xy[i] = 2.0 * ts_z_xy[i] * gfe2_0 + 4.0 * ts_yz_x[i] * gfe2_0 + 4.0 * ts_yz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_xy[i] * gfe_0 + ts_yyz_xy[i] * rgc2_0;

        gr_yyz_xz[i] = 2.0 * ts_z_xz[i] * gfe2_0 + 4.0 * ts_yz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_x[i] * gfe2_0 + 2.0 * ts_yy_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_xz[i] * gfe_0 + ts_yyz_xz[i] * rgc2_0;

        gr_yyz_yy[i] = 2.0 * ts_z_yy[i] * gfe2_0 + 8.0 * ts_yz_y[i] * gfe2_0 + 4.0 * ts_yz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_0[i] * gfe2_0 + 4.0 * ts_yyz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyz_yy[i] * gfe_0 + ts_yyz_yy[i] * rgc2_0;

        gr_yyz_yz[i] = 2.0 * ts_z_yz[i] * gfe2_0 + 4.0 * ts_yz_z[i] * gfe2_0 + 4.0 * ts_yz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe2_0 + 2.0 * ts_yy_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_yz[i] * gfe_0 + ts_yyz_yz[i] * rgc2_0;

        gr_yyz_zz[i] = 2.0 * ts_z_zz[i] * gfe2_0 + 4.0 * ts_yz_zz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yy_z[i] * gfe2_0 + 2.0 * ts_yy_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyz_0[i] * gfe2_0 + 4.0 * ts_yyz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_zz[i] * gfe_0 + ts_yyz_zz[i] * rgc2_0;
    }

    // Set up 48-54 components of targeted buffer : FD

    auto gr_yzz_xx = pbuffer.data(idx_g_fd + 48);

    auto gr_yzz_xy = pbuffer.data(idx_g_fd + 49);

    auto gr_yzz_xz = pbuffer.data(idx_g_fd + 50);

    auto gr_yzz_yy = pbuffer.data(idx_g_fd + 51);

    auto gr_yzz_yz = pbuffer.data(idx_g_fd + 52);

    auto gr_yzz_zz = pbuffer.data(idx_g_fd + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xx, gr_yzz_xy, gr_yzz_xz, gr_yzz_yy, gr_yzz_yz, gr_yzz_zz, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_yy, ts_y_yz, ts_y_zz, ts_yz_x, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_y, ts_yz_yy, ts_yz_yz, ts_yz_z, ts_yz_zz, ts_yzz_0, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_yzz_xx[i] = 2.0 * ts_zz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xx[i] * gfe2_0 + 4.0 * ts_yz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_0[i] * gfe2_0 + 4.0 * ts_yzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzz_xx[i] * gfe_0 + ts_yzz_xx[i] * rgc2_0;

        gr_yzz_xy[i] = 2.0 * ts_zz_x[i] * gfe2_0 + 2.0 * ts_zz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xy[i] * gfe2_0 + 4.0 * ts_yz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_xy[i] * gfe_0 + ts_yzz_xy[i] * rgc2_0;

        gr_yzz_xz[i] = 2.0 * ts_zz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_xz[i] * gfe2_0 + 4.0 * ts_yz_x[i] * gfe2_0 + 4.0 * ts_yz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_xz[i] * gfe_0 + ts_yzz_xz[i] * rgc2_0;

        gr_yzz_yy[i] = 4.0 * ts_zz_y[i] * gfe2_0 + 2.0 * ts_zz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yy[i] * gfe2_0 + 4.0 * ts_yz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_0[i] * gfe2_0 + 4.0 * ts_yzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzz_yy[i] * gfe_0 + ts_yzz_yy[i] * rgc2_0;

        gr_yzz_yz[i] = 2.0 * ts_zz_z[i] * gfe2_0 + 2.0 * ts_zz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_yz[i] * gfe2_0 + 4.0 * ts_yz_y[i] * gfe2_0 + 4.0 * ts_yz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_yz[i] * gfe_0 + ts_yzz_yz[i] * rgc2_0;

        gr_yzz_zz[i] = 2.0 * ts_zz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_zz[i] * gfe2_0 + 8.0 * ts_yz_z[i] * gfe2_0 + 4.0 * ts_yz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzz_0[i] * gfe2_0 + 4.0 * ts_yzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_zz[i] * gfe_0 + ts_yzz_zz[i] * rgc2_0;
    }

    // Set up 54-60 components of targeted buffer : FD

    auto gr_zzz_xx = pbuffer.data(idx_g_fd + 54);

    auto gr_zzz_xy = pbuffer.data(idx_g_fd + 55);

    auto gr_zzz_xz = pbuffer.data(idx_g_fd + 56);

    auto gr_zzz_yy = pbuffer.data(idx_g_fd + 57);

    auto gr_zzz_yz = pbuffer.data(idx_g_fd + 58);

    auto gr_zzz_zz = pbuffer.data(idx_g_fd + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xx, gr_zzz_xy, gr_zzz_xz, gr_zzz_yy, gr_zzz_yz, gr_zzz_zz, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_yy, ts_z_yz, ts_z_zz, ts_zz_x, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_y, ts_zz_yy, ts_zz_yz, ts_zz_z, ts_zz_zz, ts_zzz_0, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_zzz_xx[i] = 6.0 * ts_z_xx[i] * gfe2_0 + 6.0 * ts_zz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_0[i] * gfe2_0 + 4.0 * ts_zzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzz_xx[i] * gfe_0 + ts_zzz_xx[i] * rgc2_0;

        gr_zzz_xy[i] = 6.0 * ts_z_xy[i] * gfe2_0 + 6.0 * ts_zz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_xy[i] * gfe_0 + ts_zzz_xy[i] * rgc2_0;

        gr_zzz_xz[i] = 6.0 * ts_z_xz[i] * gfe2_0 + 6.0 * ts_zz_x[i] * gfe2_0 + 6.0 * ts_zz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_xz[i] * gfe_0 + ts_zzz_xz[i] * rgc2_0;

        gr_zzz_yy[i] = 6.0 * ts_z_yy[i] * gfe2_0 + 6.0 * ts_zz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_0[i] * gfe2_0 + 4.0 * ts_zzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzz_yy[i] * gfe_0 + ts_zzz_yy[i] * rgc2_0;

        gr_zzz_yz[i] = 6.0 * ts_z_yz[i] * gfe2_0 + 6.0 * ts_zz_y[i] * gfe2_0 + 6.0 * ts_zz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_yz[i] * gfe_0 + ts_zzz_yz[i] * rgc2_0;

        gr_zzz_zz[i] = 6.0 * ts_z_zz[i] * gfe2_0 + 12.0 * ts_zz_z[i] * gfe2_0 + 6.0 * ts_zz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzz_0[i] * gfe2_0 + 4.0 * ts_zzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_zz[i] * gfe_0 + ts_zzz_zz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

