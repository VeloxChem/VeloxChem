#include "ThreeCenterR2PrimRecGD.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_gd(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gd,
                const size_t idx_dd,
                const size_t idx_fp,
                const size_t idx_fd,
                const size_t idx_gs,
                const size_t idx_gp,
                const size_t idx_gd,
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

    // Set up components of auxiliary buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_gs);

    auto ts_xxxy_0 = pbuffer.data(idx_gs + 1);

    auto ts_xxxz_0 = pbuffer.data(idx_gs + 2);

    auto ts_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto ts_xxyz_0 = pbuffer.data(idx_gs + 4);

    auto ts_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto ts_xyyz_0 = pbuffer.data(idx_gs + 7);

    auto ts_xyzz_0 = pbuffer.data(idx_gs + 8);

    auto ts_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto ts_yyyz_0 = pbuffer.data(idx_gs + 11);

    auto ts_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_gp + 5);

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    auto ts_xxyz_x = pbuffer.data(idx_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_gp + 14);

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    auto ts_xyyz_x = pbuffer.data(idx_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    auto ts_xyzz_x = pbuffer.data(idx_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_gp + 26);

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

    auto ts_yyyz_x = pbuffer.data(idx_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_gd + 11);

    auto ts_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto ts_xxyz_xx = pbuffer.data(idx_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_gd + 29);

    auto ts_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto ts_xyyz_xx = pbuffer.data(idx_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto ts_xyzz_xx = pbuffer.data(idx_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_gd + 53);

    auto ts_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto ts_yyyz_xx = pbuffer.data(idx_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up 0-6 components of targeted buffer : GD

    auto gr_xxxx_xx = pbuffer.data(idx_g_gd);

    auto gr_xxxx_xy = pbuffer.data(idx_g_gd + 1);

    auto gr_xxxx_xz = pbuffer.data(idx_g_gd + 2);

    auto gr_xxxx_yy = pbuffer.data(idx_g_gd + 3);

    auto gr_xxxx_yz = pbuffer.data(idx_g_gd + 4);

    auto gr_xxxx_zz = pbuffer.data(idx_g_gd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_xx, gr_xxxx_xy, gr_xxxx_xz, gr_xxxx_yy, gr_xxxx_yz, gr_xxxx_zz, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_yy, ts_xx_yz, ts_xx_zz, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, ts_xxxx_0, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxx_xx[i] = 12.0 * ts_xx_xx[i] * gfe_0 + 16.0 * ts_xxx_x[i] * gfe_0 + 8.0 * ts_xxx_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_0[i] * gfe_0 + 4.0 * ts_xxxx_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxx_xx[i] * gfe_0 + ts_xxxx_xx[i] * rgc2_0;

        gr_xxxx_xy[i] = 12.0 * ts_xx_xy[i] * gfe_0 + 8.0 * ts_xxx_y[i] * gfe_0 + 8.0 * ts_xxx_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_xy[i] * gfe_0 + ts_xxxx_xy[i] * rgc2_0;

        gr_xxxx_xz[i] = 12.0 * ts_xx_xz[i] * gfe_0 + 8.0 * ts_xxx_z[i] * gfe_0 + 8.0 * ts_xxx_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_xz[i] * gfe_0 + ts_xxxx_xz[i] * rgc2_0;

        gr_xxxx_yy[i] = 12.0 * ts_xx_yy[i] * gfe_0 + 8.0 * ts_xxx_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_0[i] * gfe_0 + 4.0 * ts_xxxx_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_yy[i] * gfe_0 + ts_xxxx_yy[i] * rgc2_0;

        gr_xxxx_yz[i] = 12.0 * ts_xx_yz[i] * gfe_0 + 8.0 * ts_xxx_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxx_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_yz[i] * gfe_0 + ts_xxxx_yz[i] * rgc2_0;

        gr_xxxx_zz[i] = 12.0 * ts_xx_zz[i] * gfe_0 + 8.0 * ts_xxx_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_0[i] * gfe_0 + 4.0 * ts_xxxx_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_zz[i] * gfe_0 + ts_xxxx_zz[i] * rgc2_0;
    }

    // Set up 6-12 components of targeted buffer : GD

    auto gr_xxxy_xx = pbuffer.data(idx_g_gd + 6);

    auto gr_xxxy_xy = pbuffer.data(idx_g_gd + 7);

    auto gr_xxxy_xz = pbuffer.data(idx_g_gd + 8);

    auto gr_xxxy_yy = pbuffer.data(idx_g_gd + 9);

    auto gr_xxxy_yz = pbuffer.data(idx_g_gd + 10);

    auto gr_xxxy_zz = pbuffer.data(idx_g_gd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_xx, gr_xxxy_xy, gr_xxxy_xz, gr_xxxy_yy, gr_xxxy_yz, gr_xxxy_zz, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, ts_xxxy_0, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_yy, ts_xy_yz, ts_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxy_xx[i] = 6.0 * ts_xy_xx[i] * gfe_0 + 12.0 * ts_xxy_x[i] * gfe_0 + 6.0 * ts_xxy_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_0[i] * gfe_0 + 4.0 * ts_xxxy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxy_xx[i] * gfe_0 + ts_xxxy_xx[i] * rgc2_0;

        gr_xxxy_xy[i] = 6.0 * ts_xy_xy[i] * gfe_0 + 6.0 * ts_xxy_y[i] * gfe_0 + 6.0 * ts_xxy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe_0 + 2.0 * ts_xxx_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_xy[i] * gfe_0 + ts_xxxy_xy[i] * rgc2_0;

        gr_xxxy_xz[i] = 6.0 * ts_xy_xz[i] * gfe_0 + 6.0 * ts_xxy_z[i] * gfe_0 + 6.0 * ts_xxy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_xz[i] * gfe_0 + ts_xxxy_xz[i] * rgc2_0;

        gr_xxxy_yy[i] = 6.0 * ts_xy_yy[i] * gfe_0 + 6.0 * ts_xxy_yy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_y[i] * gfe_0 + 2.0 * ts_xxx_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_0[i] * gfe_0 + 4.0 * ts_xxxy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_yy[i] * gfe_0 + ts_xxxy_yy[i] * rgc2_0;

        gr_xxxy_yz[i] = 6.0 * ts_xy_yz[i] * gfe_0 + 6.0 * ts_xxy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_z[i] * gfe_0 + 2.0 * ts_xxx_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_yz[i] * gfe_0 + ts_xxxy_yz[i] * rgc2_0;

        gr_xxxy_zz[i] = 6.0 * ts_xy_zz[i] * gfe_0 + 6.0 * ts_xxy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_0[i] * gfe_0 + 4.0 * ts_xxxy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_zz[i] * gfe_0 + ts_xxxy_zz[i] * rgc2_0;
    }

    // Set up 12-18 components of targeted buffer : GD

    auto gr_xxxz_xx = pbuffer.data(idx_g_gd + 12);

    auto gr_xxxz_xy = pbuffer.data(idx_g_gd + 13);

    auto gr_xxxz_xz = pbuffer.data(idx_g_gd + 14);

    auto gr_xxxz_yy = pbuffer.data(idx_g_gd + 15);

    auto gr_xxxz_yz = pbuffer.data(idx_g_gd + 16);

    auto gr_xxxz_zz = pbuffer.data(idx_g_gd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_xx, gr_xxxz_xy, gr_xxxz_xz, gr_xxxz_yy, gr_xxxz_yz, gr_xxxz_zz, ts_xxx_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_y, ts_xxx_yy, ts_xxx_yz, ts_xxx_z, ts_xxx_zz, ts_xxxz_0, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_yy, ts_xz_yz, ts_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxz_xx[i] = 6.0 * ts_xz_xx[i] * gfe_0 + 12.0 * ts_xxz_x[i] * gfe_0 + 6.0 * ts_xxz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_0[i] * gfe_0 + 4.0 * ts_xxxz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxz_xx[i] * gfe_0 + ts_xxxz_xx[i] * rgc2_0;

        gr_xxxz_xy[i] = 6.0 * ts_xz_xy[i] * gfe_0 + 6.0 * ts_xxz_y[i] * gfe_0 + 6.0 * ts_xxz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_xy[i] * gfe_0 + ts_xxxz_xy[i] * rgc2_0;

        gr_xxxz_xz[i] = 6.0 * ts_xz_xz[i] * gfe_0 + 6.0 * ts_xxz_z[i] * gfe_0 + 6.0 * ts_xxz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe_0 + 2.0 * ts_xxx_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_xz[i] * gfe_0 + ts_xxxz_xz[i] * rgc2_0;

        gr_xxxz_yy[i] = 6.0 * ts_xz_yy[i] * gfe_0 + 6.0 * ts_xxz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_0[i] * gfe_0 + 4.0 * ts_xxxz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_yy[i] * gfe_0 + ts_xxxz_yy[i] * rgc2_0;

        gr_xxxz_yz[i] = 6.0 * ts_xz_yz[i] * gfe_0 + 6.0 * ts_xxz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_y[i] * gfe_0 + 2.0 * ts_xxx_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_yz[i] * gfe_0 + ts_xxxz_yz[i] * rgc2_0;

        gr_xxxz_zz[i] = 6.0 * ts_xz_zz[i] * gfe_0 + 6.0 * ts_xxz_zz[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxx_z[i] * gfe_0 + 2.0 * ts_xxx_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_0[i] * gfe_0 + 4.0 * ts_xxxz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_zz[i] * gfe_0 + ts_xxxz_zz[i] * rgc2_0;
    }

    // Set up 18-24 components of targeted buffer : GD

    auto gr_xxyy_xx = pbuffer.data(idx_g_gd + 18);

    auto gr_xxyy_xy = pbuffer.data(idx_g_gd + 19);

    auto gr_xxyy_xz = pbuffer.data(idx_g_gd + 20);

    auto gr_xxyy_yy = pbuffer.data(idx_g_gd + 21);

    auto gr_xxyy_yz = pbuffer.data(idx_g_gd + 22);

    auto gr_xxyy_zz = pbuffer.data(idx_g_gd + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_xx, gr_xxyy_xy, gr_xxyy_xz, gr_xxyy_yy, gr_xxyy_yz, gr_xxyy_zz, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_yy, ts_xx_yz, ts_xx_zz, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, ts_xxyy_0, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_yy, ts_yy_yz, ts_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxyy_xx[i] = 2.0 * ts_yy_xx[i] * gfe_0 + 8.0 * ts_xyy_x[i] * gfe_0 + 4.0 * ts_xyy_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 + 4.0 * ts_xxy_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_0[i] * gfe_0 + 4.0 * ts_xxyy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyy_xx[i] * gfe_0 + ts_xxyy_xx[i] * rgc2_0;

        gr_xxyy_xy[i] = 2.0 * ts_yy_xy[i] * gfe_0 + 4.0 * ts_xyy_y[i] * gfe_0 + 4.0 * ts_xyy_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xy[i] * gfe_0 + 4.0 * ts_xxy_x[i] * gfe_0 + 4.0 * ts_xxy_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_xy[i] * gfe_0 + ts_xxyy_xy[i] * rgc2_0;

        gr_xxyy_xz[i] = 2.0 * ts_yy_xz[i] * gfe_0 + 4.0 * ts_xyy_z[i] * gfe_0 + 4.0 * ts_xyy_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xz[i] * gfe_0 + 4.0 * ts_xxy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_xz[i] * gfe_0 + ts_xxyy_xz[i] * rgc2_0;

        gr_xxyy_yy[i] = 2.0 * ts_yy_yy[i] * gfe_0 + 4.0 * ts_xyy_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yy[i] * gfe_0 + 8.0 * ts_xxy_y[i] * gfe_0 + 4.0 * ts_xxy_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_0[i] * gfe_0 + 4.0 * ts_xxyy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_yy[i] * gfe_0 + ts_xxyy_yy[i] * rgc2_0;

        gr_xxyy_yz[i] = 2.0 * ts_yy_yz[i] * gfe_0 + 4.0 * ts_xyy_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yz[i] * gfe_0 + 4.0 * ts_xxy_z[i] * gfe_0 + 4.0 * ts_xxy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_yz[i] * gfe_0 + ts_xxyy_yz[i] * rgc2_0;

        gr_xxyy_zz[i] = 2.0 * ts_yy_zz[i] * gfe_0 + 4.0 * ts_xyy_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 + 4.0 * ts_xxy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_0[i] * gfe_0 + 4.0 * ts_xxyy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_zz[i] * gfe_0 + ts_xxyy_zz[i] * rgc2_0;
    }

    // Set up 24-30 components of targeted buffer : GD

    auto gr_xxyz_xx = pbuffer.data(idx_g_gd + 24);

    auto gr_xxyz_xy = pbuffer.data(idx_g_gd + 25);

    auto gr_xxyz_xz = pbuffer.data(idx_g_gd + 26);

    auto gr_xxyz_yy = pbuffer.data(idx_g_gd + 27);

    auto gr_xxyz_yz = pbuffer.data(idx_g_gd + 28);

    auto gr_xxyz_zz = pbuffer.data(idx_g_gd + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_xx, gr_xxyz_xy, gr_xxyz_xz, gr_xxyz_yy, gr_xxyz_yz, gr_xxyz_zz, ts_xxy_x, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_y, ts_xxy_yy, ts_xxy_yz, ts_xxy_z, ts_xxy_zz, ts_xxyz_0, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_yy, ts_yz_yz, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxyz_xx[i] = 2.0 * ts_yz_xx[i] * gfe_0 + 8.0 * ts_xyz_x[i] * gfe_0 + 4.0 * ts_xyz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_0[i] * gfe_0 + 4.0 * ts_xxyz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyz_xx[i] * gfe_0 + ts_xxyz_xx[i] * rgc2_0;

        gr_xxyz_xy[i] = 2.0 * ts_yz_xy[i] * gfe_0 + 4.0 * ts_xyz_y[i] * gfe_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_x[i] * gfe_0 + 2.0 * ts_xxz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_xy[i] * gfe_0 + ts_xxyz_xy[i] * rgc2_0;

        gr_xxyz_xz[i] = 2.0 * ts_yz_xz[i] * gfe_0 + 4.0 * ts_xyz_z[i] * gfe_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_x[i] * gfe_0 + 2.0 * ts_xxy_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxyz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_xz[i] * gfe_0 + ts_xxyz_xz[i] * rgc2_0;

        gr_xxyz_yy[i] = 2.0 * ts_yz_yy[i] * gfe_0 + 4.0 * ts_xyz_yy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xxz_y[i] * gfe_0 + 2.0 * ts_xxz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_0[i] * gfe_0 + 4.0 * ts_xxyz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_yy[i] * gfe_0 + ts_xxyz_yy[i] * rgc2_0;

        gr_xxyz_yz[i] = 2.0 * ts_yz_yz[i] * gfe_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_z[i] * gfe_0 + 2.0 * ts_xxz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_y[i] * gfe_0 + 2.0 * ts_xxy_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_yz[i] * gfe_0 + ts_xxyz_yz[i] * rgc2_0;

        gr_xxyz_zz[i] = 2.0 * ts_yz_zz[i] * gfe_0 + 4.0 * ts_xyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_zz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xxy_z[i] * gfe_0 + 2.0 * ts_xxy_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_0[i] * gfe_0 + 4.0 * ts_xxyz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_zz[i] * gfe_0 + ts_xxyz_zz[i] * rgc2_0;
    }

    // Set up 30-36 components of targeted buffer : GD

    auto gr_xxzz_xx = pbuffer.data(idx_g_gd + 30);

    auto gr_xxzz_xy = pbuffer.data(idx_g_gd + 31);

    auto gr_xxzz_xz = pbuffer.data(idx_g_gd + 32);

    auto gr_xxzz_yy = pbuffer.data(idx_g_gd + 33);

    auto gr_xxzz_yz = pbuffer.data(idx_g_gd + 34);

    auto gr_xxzz_zz = pbuffer.data(idx_g_gd + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_xx, gr_xxzz_xy, gr_xxzz_xz, gr_xxzz_yy, gr_xxzz_yz, gr_xxzz_zz, ts_xx_xx, ts_xx_xy, ts_xx_xz, ts_xx_yy, ts_xx_yz, ts_xx_zz, ts_xxz_x, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_y, ts_xxz_yy, ts_xxz_yz, ts_xxz_z, ts_xxz_zz, ts_xxzz_0, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_yy, ts_zz_yz, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxzz_xx[i] = 2.0 * ts_zz_xx[i] * gfe_0 + 8.0 * ts_xzz_x[i] * gfe_0 + 4.0 * ts_xzz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xx[i] * gfe_0 + 4.0 * ts_xxz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_0[i] * gfe_0 + 4.0 * ts_xxzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxzz_xx[i] * gfe_0 + ts_xxzz_xx[i] * rgc2_0;

        gr_xxzz_xy[i] = 2.0 * ts_zz_xy[i] * gfe_0 + 4.0 * ts_xzz_y[i] * gfe_0 + 4.0 * ts_xzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xy[i] * gfe_0 + 4.0 * ts_xxz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_xy[i] * gfe_0 + ts_xxzz_xy[i] * rgc2_0;

        gr_xxzz_xz[i] = 2.0 * ts_zz_xz[i] * gfe_0 + 4.0 * ts_xzz_z[i] * gfe_0 + 4.0 * ts_xzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_xz[i] * gfe_0 + 4.0 * ts_xxz_x[i] * gfe_0 + 4.0 * ts_xxz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_xz[i] * gfe_0 + ts_xxzz_xz[i] * rgc2_0;

        gr_xxzz_yy[i] = 2.0 * ts_zz_yy[i] * gfe_0 + 4.0 * ts_xzz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yy[i] * gfe_0 + 4.0 * ts_xxz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_0[i] * gfe_0 + 4.0 * ts_xxzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_yy[i] * gfe_0 + ts_xxzz_yy[i] * rgc2_0;

        gr_xxzz_yz[i] = 2.0 * ts_zz_yz[i] * gfe_0 + 4.0 * ts_xzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_yz[i] * gfe_0 + 4.0 * ts_xxz_y[i] * gfe_0 + 4.0 * ts_xxz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_yz[i] * gfe_0 + ts_xxzz_yz[i] * rgc2_0;

        gr_xxzz_zz[i] = 2.0 * ts_zz_zz[i] * gfe_0 + 4.0 * ts_xzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_zz[i] * gfe_0 + 8.0 * ts_xxz_z[i] * gfe_0 + 4.0 * ts_xxz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_0[i] * gfe_0 + 4.0 * ts_xxzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_zz[i] * gfe_0 + ts_xxzz_zz[i] * rgc2_0;
    }

    // Set up 36-42 components of targeted buffer : GD

    auto gr_xyyy_xx = pbuffer.data(idx_g_gd + 36);

    auto gr_xyyy_xy = pbuffer.data(idx_g_gd + 37);

    auto gr_xyyy_xz = pbuffer.data(idx_g_gd + 38);

    auto gr_xyyy_yy = pbuffer.data(idx_g_gd + 39);

    auto gr_xyyy_yz = pbuffer.data(idx_g_gd + 40);

    auto gr_xyyy_zz = pbuffer.data(idx_g_gd + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_xx, gr_xyyy_xy, gr_xyyy_xz, gr_xyyy_yy, gr_xyyy_yz, gr_xyyy_zz, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_yy, ts_xy_yz, ts_xy_zz, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, ts_xyyy_0, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyyy_xx[i] = 4.0 * ts_yyy_x[i] * gfe_0 + 2.0 * ts_yyy_xx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xx[i] * gfe_0 + 6.0 * ts_xyy_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_0[i] * gfe_0 + 4.0 * ts_xyyy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyy_xx[i] * gfe_0 + ts_xyyy_xx[i] * rgc2_0;

        gr_xyyy_xy[i] = 2.0 * ts_yyy_y[i] * gfe_0 + 2.0 * ts_yyy_xy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xy[i] * gfe_0 + 6.0 * ts_xyy_x[i] * gfe_0 + 6.0 * ts_xyy_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_xy[i] * gfe_0 + ts_xyyy_xy[i] * rgc2_0;

        gr_xyyy_xz[i] = 2.0 * ts_yyy_z[i] * gfe_0 + 2.0 * ts_yyy_xz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_xz[i] * gfe_0 + 6.0 * ts_xyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_xz[i] * gfe_0 + ts_xyyy_xz[i] * rgc2_0;

        gr_xyyy_yy[i] = 2.0 * ts_yyy_yy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yy[i] * gfe_0 + 12.0 * ts_xyy_y[i] * gfe_0 + 6.0 * ts_xyy_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_0[i] * gfe_0 + 4.0 * ts_xyyy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_yy[i] * gfe_0 + ts_xyyy_yy[i] * rgc2_0;

        gr_xyyy_yz[i] = 2.0 * ts_yyy_yz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_yz[i] * gfe_0 + 6.0 * ts_xyy_z[i] * gfe_0 + 6.0 * ts_xyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_yz[i] * gfe_0 + ts_xyyy_yz[i] * rgc2_0;

        gr_xyyy_zz[i] = 2.0 * ts_yyy_zz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_zz[i] * gfe_0 + 6.0 * ts_xyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_0[i] * gfe_0 + 4.0 * ts_xyyy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_zz[i] * gfe_0 + ts_xyyy_zz[i] * rgc2_0;
    }

    // Set up 42-48 components of targeted buffer : GD

    auto gr_xyyz_xx = pbuffer.data(idx_g_gd + 42);

    auto gr_xyyz_xy = pbuffer.data(idx_g_gd + 43);

    auto gr_xyyz_xz = pbuffer.data(idx_g_gd + 44);

    auto gr_xyyz_yy = pbuffer.data(idx_g_gd + 45);

    auto gr_xyyz_yz = pbuffer.data(idx_g_gd + 46);

    auto gr_xyyz_zz = pbuffer.data(idx_g_gd + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_xx, gr_xyyz_xy, gr_xyyz_xz, gr_xyyz_yy, gr_xyyz_yz, gr_xyyz_zz, ts_xyy_x, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_y, ts_xyy_yy, ts_xyy_yz, ts_xyy_z, ts_xyy_zz, ts_xyyz_0, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_yy, ts_xz_yz, ts_xz_zz, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyyz_xx[i] = 4.0 * ts_yyz_x[i] * gfe_0 + 2.0 * ts_yyz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xx[i] * gfe_0 + 4.0 * ts_xyz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_0[i] * gfe_0 + 4.0 * ts_xyyz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyz_xx[i] * gfe_0 + ts_xyyz_xx[i] * rgc2_0;

        gr_xyyz_xy[i] = 2.0 * ts_yyz_y[i] * gfe_0 + 2.0 * ts_yyz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xy[i] * gfe_0 + 4.0 * ts_xyz_x[i] * gfe_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_xy[i] * gfe_0 + ts_xyyz_xy[i] * rgc2_0;

        gr_xyyz_xz[i] = 2.0 * ts_yyz_z[i] * gfe_0 + 2.0 * ts_yyz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_xz[i] * gfe_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_x[i] * gfe_0 + 2.0 * ts_xyy_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyyz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_xz[i] * gfe_0 + ts_xyyz_xz[i] * rgc2_0;

        gr_xyyz_yy[i] = 2.0 * ts_yyz_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yy[i] * gfe_0 + 8.0 * ts_xyz_y[i] * gfe_0 + 4.0 * ts_xyz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_0[i] * gfe_0 + 4.0 * ts_xyyz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_yy[i] * gfe_0 + ts_xyyz_yy[i] * rgc2_0;

        gr_xyyz_yz[i] = 2.0 * ts_yyz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_yz[i] * gfe_0 + 4.0 * ts_xyz_z[i] * gfe_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_y[i] * gfe_0 + 2.0 * ts_xyy_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_yz[i] * gfe_0 + ts_xyyz_yz[i] * rgc2_0;

        gr_xyyz_zz[i] = 2.0 * ts_yyz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_zz[i] * gfe_0 + 4.0 * ts_xyz_zz[i] * gfe_0 * gc_y[i] + 4.0 * ts_xyy_z[i] * gfe_0 + 2.0 * ts_xyy_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_0[i] * gfe_0 + 4.0 * ts_xyyz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_zz[i] * gfe_0 + ts_xyyz_zz[i] * rgc2_0;
    }

    // Set up 48-54 components of targeted buffer : GD

    auto gr_xyzz_xx = pbuffer.data(idx_g_gd + 48);

    auto gr_xyzz_xy = pbuffer.data(idx_g_gd + 49);

    auto gr_xyzz_xz = pbuffer.data(idx_g_gd + 50);

    auto gr_xyzz_yy = pbuffer.data(idx_g_gd + 51);

    auto gr_xyzz_yz = pbuffer.data(idx_g_gd + 52);

    auto gr_xyzz_zz = pbuffer.data(idx_g_gd + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_xx, gr_xyzz_xy, gr_xyzz_xz, gr_xyzz_yy, gr_xyzz_yz, gr_xyzz_zz, ts_xy_xx, ts_xy_xy, ts_xy_xz, ts_xy_yy, ts_xy_yz, ts_xy_zz, ts_xyz_x, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_y, ts_xyz_yy, ts_xyz_yz, ts_xyz_z, ts_xyz_zz, ts_xyzz_0, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyzz_xx[i] = 4.0 * ts_yzz_x[i] * gfe_0 + 2.0 * ts_yzz_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xx[i] * gfe_0 + 4.0 * ts_xyz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_0[i] * gfe_0 + 4.0 * ts_xyzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyzz_xx[i] * gfe_0 + ts_xyzz_xx[i] * rgc2_0;

        gr_xyzz_xy[i] = 2.0 * ts_yzz_y[i] * gfe_0 + 2.0 * ts_yzz_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_x[i] * gfe_0 + 2.0 * ts_xzz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xy[i] * gfe_0 + 4.0 * ts_xyz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_xy[i] * gfe_0 + ts_xyzz_xy[i] * rgc2_0;

        gr_xyzz_xz[i] = 2.0 * ts_yzz_z[i] * gfe_0 + 2.0 * ts_yzz_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_xz[i] * gfe_0 + 4.0 * ts_xyz_x[i] * gfe_0 + 4.0 * ts_xyz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xyzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_xz[i] * gfe_0 + ts_xyzz_xz[i] * rgc2_0;

        gr_xyzz_yy[i] = 2.0 * ts_yzz_yy[i] * gfe_0 * gc_x[i] + 4.0 * ts_xzz_y[i] * gfe_0 + 2.0 * ts_xzz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yy[i] * gfe_0 + 4.0 * ts_xyz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_0[i] * gfe_0 + 4.0 * ts_xyzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_yy[i] * gfe_0 + ts_xyzz_yy[i] * rgc2_0;

        gr_xyzz_yz[i] = 2.0 * ts_yzz_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_z[i] * gfe_0 + 2.0 * ts_xzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_yz[i] * gfe_0 + 4.0 * ts_xyz_y[i] * gfe_0 + 4.0 * ts_xyz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_yz[i] * gfe_0 + ts_xyzz_yz[i] * rgc2_0;

        gr_xyzz_zz[i] = 2.0 * ts_yzz_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_zz[i] * gfe_0 + 8.0 * ts_xyz_z[i] * gfe_0 + 4.0 * ts_xyz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_0[i] * gfe_0 + 4.0 * ts_xyzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_zz[i] * gfe_0 + ts_xyzz_zz[i] * rgc2_0;
    }

    // Set up 54-60 components of targeted buffer : GD

    auto gr_xzzz_xx = pbuffer.data(idx_g_gd + 54);

    auto gr_xzzz_xy = pbuffer.data(idx_g_gd + 55);

    auto gr_xzzz_xz = pbuffer.data(idx_g_gd + 56);

    auto gr_xzzz_yy = pbuffer.data(idx_g_gd + 57);

    auto gr_xzzz_yz = pbuffer.data(idx_g_gd + 58);

    auto gr_xzzz_zz = pbuffer.data(idx_g_gd + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_xx, gr_xzzz_xy, gr_xzzz_xz, gr_xzzz_yy, gr_xzzz_yz, gr_xzzz_zz, ts_xz_xx, ts_xz_xy, ts_xz_xz, ts_xz_yy, ts_xz_yz, ts_xz_zz, ts_xzz_x, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_y, ts_xzz_yy, ts_xzz_yz, ts_xzz_z, ts_xzz_zz, ts_xzzz_0, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xzzz_xx[i] = 4.0 * ts_zzz_x[i] * gfe_0 + 2.0 * ts_zzz_xx[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xx[i] * gfe_0 + 6.0 * ts_xzz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_0[i] * gfe_0 + 4.0 * ts_xzzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzzz_xx[i] * gfe_0 + ts_xzzz_xx[i] * rgc2_0;

        gr_xzzz_xy[i] = 2.0 * ts_zzz_y[i] * gfe_0 + 2.0 * ts_zzz_xy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xy[i] * gfe_0 + 6.0 * ts_xzz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_xy[i] * gfe_0 + ts_xzzz_xy[i] * rgc2_0;

        gr_xzzz_xz[i] = 2.0 * ts_zzz_z[i] * gfe_0 + 2.0 * ts_zzz_xz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_xz[i] * gfe_0 + 6.0 * ts_xzz_x[i] * gfe_0 + 6.0 * ts_xzz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_xz[i] * gfe_0 + ts_xzzz_xz[i] * rgc2_0;

        gr_xzzz_yy[i] = 2.0 * ts_zzz_yy[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yy[i] * gfe_0 + 6.0 * ts_xzz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_0[i] * gfe_0 + 4.0 * ts_xzzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_yy[i] * gfe_0 + ts_xzzz_yy[i] * rgc2_0;

        gr_xzzz_yz[i] = 2.0 * ts_zzz_yz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_yz[i] * gfe_0 + 6.0 * ts_xzz_y[i] * gfe_0 + 6.0 * ts_xzz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xzzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_yz[i] * gfe_0 + ts_xzzz_yz[i] * rgc2_0;

        gr_xzzz_zz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_zz[i] * gfe_0 + 12.0 * ts_xzz_z[i] * gfe_0 + 6.0 * ts_xzz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_0[i] * gfe_0 + 4.0 * ts_xzzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_zz[i] * gfe_0 + ts_xzzz_zz[i] * rgc2_0;
    }

    // Set up 60-66 components of targeted buffer : GD

    auto gr_yyyy_xx = pbuffer.data(idx_g_gd + 60);

    auto gr_yyyy_xy = pbuffer.data(idx_g_gd + 61);

    auto gr_yyyy_xz = pbuffer.data(idx_g_gd + 62);

    auto gr_yyyy_yy = pbuffer.data(idx_g_gd + 63);

    auto gr_yyyy_yz = pbuffer.data(idx_g_gd + 64);

    auto gr_yyyy_zz = pbuffer.data(idx_g_gd + 65);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_xx, gr_yyyy_xy, gr_yyyy_xz, gr_yyyy_yy, gr_yyyy_yz, gr_yyyy_zz, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_yy, ts_yy_yz, ts_yy_zz, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, ts_yyyy_0, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyyy_xx[i] = 12.0 * ts_yy_xx[i] * gfe_0 + 8.0 * ts_yyy_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_0[i] * gfe_0 + 4.0 * ts_yyyy_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyy_xx[i] * gfe_0 + ts_yyyy_xx[i] * rgc2_0;

        gr_yyyy_xy[i] = 12.0 * ts_yy_xy[i] * gfe_0 + 8.0 * ts_yyy_x[i] * gfe_0 + 8.0 * ts_yyy_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_xy[i] * gfe_0 + ts_yyyy_xy[i] * rgc2_0;

        gr_yyyy_xz[i] = 12.0 * ts_yy_xz[i] * gfe_0 + 8.0 * ts_yyy_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyy_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_xz[i] * gfe_0 + ts_yyyy_xz[i] * rgc2_0;

        gr_yyyy_yy[i] = 12.0 * ts_yy_yy[i] * gfe_0 + 16.0 * ts_yyy_y[i] * gfe_0 + 8.0 * ts_yyy_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_0[i] * gfe_0 + 4.0 * ts_yyyy_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_yy[i] * gfe_0 + ts_yyyy_yy[i] * rgc2_0;

        gr_yyyy_yz[i] = 12.0 * ts_yy_yz[i] * gfe_0 + 8.0 * ts_yyy_z[i] * gfe_0 + 8.0 * ts_yyy_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_yz[i] * gfe_0 + ts_yyyy_yz[i] * rgc2_0;

        gr_yyyy_zz[i] = 12.0 * ts_yy_zz[i] * gfe_0 + 8.0 * ts_yyy_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_0[i] * gfe_0 + 4.0 * ts_yyyy_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_zz[i] * gfe_0 + ts_yyyy_zz[i] * rgc2_0;
    }

    // Set up 66-72 components of targeted buffer : GD

    auto gr_yyyz_xx = pbuffer.data(idx_g_gd + 66);

    auto gr_yyyz_xy = pbuffer.data(idx_g_gd + 67);

    auto gr_yyyz_xz = pbuffer.data(idx_g_gd + 68);

    auto gr_yyyz_yy = pbuffer.data(idx_g_gd + 69);

    auto gr_yyyz_yz = pbuffer.data(idx_g_gd + 70);

    auto gr_yyyz_zz = pbuffer.data(idx_g_gd + 71);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_xx, gr_yyyz_xy, gr_yyyz_xz, gr_yyyz_yy, gr_yyyz_yz, gr_yyyz_zz, ts_yyy_x, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_y, ts_yyy_yy, ts_yyy_yz, ts_yyy_z, ts_yyy_zz, ts_yyyz_0, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_yy, ts_yz_yz, ts_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyyz_xx[i] = 6.0 * ts_yz_xx[i] * gfe_0 + 6.0 * ts_yyz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_0[i] * gfe_0 + 4.0 * ts_yyyz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyz_xx[i] * gfe_0 + ts_yyyz_xx[i] * rgc2_0;

        gr_yyyz_xy[i] = 6.0 * ts_yz_xy[i] * gfe_0 + 6.0 * ts_yyz_x[i] * gfe_0 + 6.0 * ts_yyz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_xy[i] * gfe_0 + ts_yyyz_xy[i] * rgc2_0;

        gr_yyyz_xz[i] = 6.0 * ts_yz_xz[i] * gfe_0 + 6.0 * ts_yyz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_x[i] * gfe_0 + 2.0 * ts_yyy_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyyz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_xz[i] * gfe_0 + ts_yyyz_xz[i] * rgc2_0;

        gr_yyyz_yy[i] = 6.0 * ts_yz_yy[i] * gfe_0 + 12.0 * ts_yyz_y[i] * gfe_0 + 6.0 * ts_yyz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_0[i] * gfe_0 + 4.0 * ts_yyyz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_yy[i] * gfe_0 + ts_yyyz_yy[i] * rgc2_0;

        gr_yyyz_yz[i] = 6.0 * ts_yz_yz[i] * gfe_0 + 6.0 * ts_yyz_z[i] * gfe_0 + 6.0 * ts_yyz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_y[i] * gfe_0 + 2.0 * ts_yyy_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_yz[i] * gfe_0 + ts_yyyz_yz[i] * rgc2_0;

        gr_yyyz_zz[i] = 6.0 * ts_yz_zz[i] * gfe_0 + 6.0 * ts_yyz_zz[i] * gfe_0 * gc_y[i] + 4.0 * ts_yyy_z[i] * gfe_0 + 2.0 * ts_yyy_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_0[i] * gfe_0 + 4.0 * ts_yyyz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_zz[i] * gfe_0 + ts_yyyz_zz[i] * rgc2_0;
    }

    // Set up 72-78 components of targeted buffer : GD

    auto gr_yyzz_xx = pbuffer.data(idx_g_gd + 72);

    auto gr_yyzz_xy = pbuffer.data(idx_g_gd + 73);

    auto gr_yyzz_xz = pbuffer.data(idx_g_gd + 74);

    auto gr_yyzz_yy = pbuffer.data(idx_g_gd + 75);

    auto gr_yyzz_yz = pbuffer.data(idx_g_gd + 76);

    auto gr_yyzz_zz = pbuffer.data(idx_g_gd + 77);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_xx, gr_yyzz_xy, gr_yyzz_xz, gr_yyzz_yy, gr_yyzz_yz, gr_yyzz_zz, ts_yy_xx, ts_yy_xy, ts_yy_xz, ts_yy_yy, ts_yy_yz, ts_yy_zz, ts_yyz_x, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_y, ts_yyz_yy, ts_yyz_yz, ts_yyz_z, ts_yyz_zz, ts_yyzz_0, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_yy, ts_zz_yz, ts_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyzz_xx[i] = 2.0 * ts_zz_xx[i] * gfe_0 + 4.0 * ts_yzz_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xx[i] * gfe_0 + 4.0 * ts_yyz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_0[i] * gfe_0 + 4.0 * ts_yyzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyzz_xx[i] * gfe_0 + ts_yyzz_xx[i] * rgc2_0;

        gr_yyzz_xy[i] = 2.0 * ts_zz_xy[i] * gfe_0 + 4.0 * ts_yzz_x[i] * gfe_0 + 4.0 * ts_yzz_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xy[i] * gfe_0 + 4.0 * ts_yyz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_xy[i] * gfe_0 + ts_yyzz_xy[i] * rgc2_0;

        gr_yyzz_xz[i] = 2.0 * ts_zz_xz[i] * gfe_0 + 4.0 * ts_yzz_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_xz[i] * gfe_0 + 4.0 * ts_yyz_x[i] * gfe_0 + 4.0 * ts_yyz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yyzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_xz[i] * gfe_0 + ts_yyzz_xz[i] * rgc2_0;

        gr_yyzz_yy[i] = 2.0 * ts_zz_yy[i] * gfe_0 + 8.0 * ts_yzz_y[i] * gfe_0 + 4.0 * ts_yzz_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yy[i] * gfe_0 + 4.0 * ts_yyz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_0[i] * gfe_0 + 4.0 * ts_yyzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_yy[i] * gfe_0 + ts_yyzz_yy[i] * rgc2_0;

        gr_yyzz_yz[i] = 2.0 * ts_zz_yz[i] * gfe_0 + 4.0 * ts_yzz_z[i] * gfe_0 + 4.0 * ts_yzz_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_yz[i] * gfe_0 + 4.0 * ts_yyz_y[i] * gfe_0 + 4.0 * ts_yyz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_yz[i] * gfe_0 + ts_yyzz_yz[i] * rgc2_0;

        gr_yyzz_zz[i] = 2.0 * ts_zz_zz[i] * gfe_0 + 4.0 * ts_yzz_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_zz[i] * gfe_0 + 8.0 * ts_yyz_z[i] * gfe_0 + 4.0 * ts_yyz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_0[i] * gfe_0 + 4.0 * ts_yyzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_zz[i] * gfe_0 + ts_yyzz_zz[i] * rgc2_0;
    }

    // Set up 78-84 components of targeted buffer : GD

    auto gr_yzzz_xx = pbuffer.data(idx_g_gd + 78);

    auto gr_yzzz_xy = pbuffer.data(idx_g_gd + 79);

    auto gr_yzzz_xz = pbuffer.data(idx_g_gd + 80);

    auto gr_yzzz_yy = pbuffer.data(idx_g_gd + 81);

    auto gr_yzzz_yz = pbuffer.data(idx_g_gd + 82);

    auto gr_yzzz_zz = pbuffer.data(idx_g_gd + 83);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_xx, gr_yzzz_xy, gr_yzzz_xz, gr_yzzz_yy, gr_yzzz_yz, gr_yzzz_zz, ts_yz_xx, ts_yz_xy, ts_yz_xz, ts_yz_yy, ts_yz_yz, ts_yz_zz, ts_yzz_x, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_y, ts_yzz_yy, ts_yzz_yz, ts_yzz_z, ts_yzz_zz, ts_yzzz_0, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yzzz_xx[i] = 2.0 * ts_zzz_xx[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xx[i] * gfe_0 + 6.0 * ts_yzz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_0[i] * gfe_0 + 4.0 * ts_yzzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzzz_xx[i] * gfe_0 + ts_yzzz_xx[i] * rgc2_0;

        gr_yzzz_xy[i] = 2.0 * ts_zzz_x[i] * gfe_0 + 2.0 * ts_zzz_xy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xy[i] * gfe_0 + 6.0 * ts_yzz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_xy[i] * gfe_0 + ts_yzzz_xy[i] * rgc2_0;

        gr_yzzz_xz[i] = 2.0 * ts_zzz_xz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_xz[i] * gfe_0 + 6.0 * ts_yzz_x[i] * gfe_0 + 6.0 * ts_yzz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_yzzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_xz[i] * gfe_0 + ts_yzzz_xz[i] * rgc2_0;

        gr_yzzz_yy[i] = 4.0 * ts_zzz_y[i] * gfe_0 + 2.0 * ts_zzz_yy[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yy[i] * gfe_0 + 6.0 * ts_yzz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_0[i] * gfe_0 + 4.0 * ts_yzzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_yy[i] * gfe_0 + ts_yzzz_yy[i] * rgc2_0;

        gr_yzzz_yz[i] = 2.0 * ts_zzz_z[i] * gfe_0 + 2.0 * ts_zzz_yz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_yz[i] * gfe_0 + 6.0 * ts_yzz_y[i] * gfe_0 + 6.0 * ts_yzz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yzzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_yz[i] * gfe_0 + ts_yzzz_yz[i] * rgc2_0;

        gr_yzzz_zz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_zz[i] * gfe_0 + 12.0 * ts_yzz_z[i] * gfe_0 + 6.0 * ts_yzz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_0[i] * gfe_0 + 4.0 * ts_yzzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_zz[i] * gfe_0 + ts_yzzz_zz[i] * rgc2_0;
    }

    // Set up 84-90 components of targeted buffer : GD

    auto gr_zzzz_xx = pbuffer.data(idx_g_gd + 84);

    auto gr_zzzz_xy = pbuffer.data(idx_g_gd + 85);

    auto gr_zzzz_xz = pbuffer.data(idx_g_gd + 86);

    auto gr_zzzz_yy = pbuffer.data(idx_g_gd + 87);

    auto gr_zzzz_yz = pbuffer.data(idx_g_gd + 88);

    auto gr_zzzz_zz = pbuffer.data(idx_g_gd + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_xx, gr_zzzz_xy, gr_zzzz_xz, gr_zzzz_yy, gr_zzzz_yz, gr_zzzz_zz, ts_zz_xx, ts_zz_xy, ts_zz_xz, ts_zz_yy, ts_zz_yz, ts_zz_zz, ts_zzz_x, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_y, ts_zzz_yy, ts_zzz_yz, ts_zzz_z, ts_zzz_zz, ts_zzzz_0, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zzzz_xx[i] = 12.0 * ts_zz_xx[i] * gfe_0 + 8.0 * ts_zzz_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_0[i] * gfe_0 + 4.0 * ts_zzzz_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzzz_xx[i] * gfe_0 + ts_zzzz_xx[i] * rgc2_0;

        gr_zzzz_xy[i] = 12.0 * ts_zz_xy[i] * gfe_0 + 8.0 * ts_zzz_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_xy[i] * gfe_0 + ts_zzzz_xy[i] * rgc2_0;

        gr_zzzz_xz[i] = 12.0 * ts_zz_xz[i] * gfe_0 + 8.0 * ts_zzz_x[i] * gfe_0 + 8.0 * ts_zzz_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_zzzz_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_xz[i] * gfe_0 + ts_zzzz_xz[i] * rgc2_0;

        gr_zzzz_yy[i] = 12.0 * ts_zz_yy[i] * gfe_0 + 8.0 * ts_zzz_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_0[i] * gfe_0 + 4.0 * ts_zzzz_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_yy[i] * gfe_0 + ts_zzzz_yy[i] * rgc2_0;

        gr_zzzz_yz[i] = 12.0 * ts_zz_yz[i] * gfe_0 + 8.0 * ts_zzz_y[i] * gfe_0 + 8.0 * ts_zzz_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_zzzz_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_yz[i] * gfe_0 + ts_zzzz_yz[i] * rgc2_0;

        gr_zzzz_zz[i] = 12.0 * ts_zz_zz[i] * gfe_0 + 16.0 * ts_zzz_z[i] * gfe_0 + 8.0 * ts_zzz_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_0[i] * gfe_0 + 4.0 * ts_zzzz_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_zz[i] * gfe_0 + ts_zzzz_zz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

