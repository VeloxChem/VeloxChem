#include "ThreeCenterOverlapGradientPrimRecGD.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_gd(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_gd,
                              const size_t idx_fd,
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

    auto gs_x_xxxx_xx = pbuffer.data(idx_g_gd);

    auto gs_x_xxxx_xy = pbuffer.data(idx_g_gd + 1);

    auto gs_x_xxxx_xz = pbuffer.data(idx_g_gd + 2);

    auto gs_x_xxxx_yy = pbuffer.data(idx_g_gd + 3);

    auto gs_x_xxxx_yz = pbuffer.data(idx_g_gd + 4);

    auto gs_x_xxxx_zz = pbuffer.data(idx_g_gd + 5);

    #pragma omp simd aligned(gc_x, gs_x_xxxx_xx, gs_x_xxxx_xy, gs_x_xxxx_xz, gs_x_xxxx_yy, gs_x_xxxx_yz, gs_x_xxxx_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxx_xx[i] = 8.0 * ts_xxx_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xy[i] = 8.0 * ts_xxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_xz[i] = 8.0 * ts_xxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yy[i] = 8.0 * ts_xxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxx_yz[i] = 8.0 * ts_xxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxx_zz[i] = 8.0 * ts_xxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 6-12 components of targeted buffer : GD

    auto gs_x_xxxy_xx = pbuffer.data(idx_g_gd + 6);

    auto gs_x_xxxy_xy = pbuffer.data(idx_g_gd + 7);

    auto gs_x_xxxy_xz = pbuffer.data(idx_g_gd + 8);

    auto gs_x_xxxy_yy = pbuffer.data(idx_g_gd + 9);

    auto gs_x_xxxy_yz = pbuffer.data(idx_g_gd + 10);

    auto gs_x_xxxy_zz = pbuffer.data(idx_g_gd + 11);

    #pragma omp simd aligned(gc_x, gs_x_xxxy_xx, gs_x_xxxy_xy, gs_x_xxxy_xz, gs_x_xxxy_yy, gs_x_xxxy_yz, gs_x_xxxy_zz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_yy, ts_xxy_yz, ts_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxy_xx[i] = 6.0 * ts_xxy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xy[i] = 6.0 * ts_xxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_xz[i] = 6.0 * ts_xxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yy[i] = 6.0 * ts_xxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxy_yz[i] = 6.0 * ts_xxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxy_zz[i] = 6.0 * ts_xxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 12-18 components of targeted buffer : GD

    auto gs_x_xxxz_xx = pbuffer.data(idx_g_gd + 12);

    auto gs_x_xxxz_xy = pbuffer.data(idx_g_gd + 13);

    auto gs_x_xxxz_xz = pbuffer.data(idx_g_gd + 14);

    auto gs_x_xxxz_yy = pbuffer.data(idx_g_gd + 15);

    auto gs_x_xxxz_yz = pbuffer.data(idx_g_gd + 16);

    auto gs_x_xxxz_zz = pbuffer.data(idx_g_gd + 17);

    #pragma omp simd aligned(gc_x, gs_x_xxxz_xx, gs_x_xxxz_xy, gs_x_xxxz_xz, gs_x_xxxz_yy, gs_x_xxxz_yz, gs_x_xxxz_zz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_yy, ts_xxz_yz, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxz_xx[i] = 6.0 * ts_xxz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xy[i] = 6.0 * ts_xxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_xz[i] = 6.0 * ts_xxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yy[i] = 6.0 * ts_xxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxxz_yz[i] = 6.0 * ts_xxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxxz_zz[i] = 6.0 * ts_xxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 18-24 components of targeted buffer : GD

    auto gs_x_xxyy_xx = pbuffer.data(idx_g_gd + 18);

    auto gs_x_xxyy_xy = pbuffer.data(idx_g_gd + 19);

    auto gs_x_xxyy_xz = pbuffer.data(idx_g_gd + 20);

    auto gs_x_xxyy_yy = pbuffer.data(idx_g_gd + 21);

    auto gs_x_xxyy_yz = pbuffer.data(idx_g_gd + 22);

    auto gs_x_xxyy_zz = pbuffer.data(idx_g_gd + 23);

    #pragma omp simd aligned(gc_x, gs_x_xxyy_xx, gs_x_xxyy_xy, gs_x_xxyy_xz, gs_x_xxyy_yy, gs_x_xxyy_yz, gs_x_xxyy_zz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyy_xx[i] = 4.0 * ts_xyy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xx[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xy[i] = 4.0 * ts_xyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_xz[i] = 4.0 * ts_xyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yy[i] = 4.0 * ts_xyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yy[i] * gc_x[i] * tce_0;

        gs_x_xxyy_yz[i] = 4.0 * ts_xyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yz[i] * gc_x[i] * tce_0;

        gs_x_xxyy_zz[i] = 4.0 * ts_xyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 24-30 components of targeted buffer : GD

    auto gs_x_xxyz_xx = pbuffer.data(idx_g_gd + 24);

    auto gs_x_xxyz_xy = pbuffer.data(idx_g_gd + 25);

    auto gs_x_xxyz_xz = pbuffer.data(idx_g_gd + 26);

    auto gs_x_xxyz_yy = pbuffer.data(idx_g_gd + 27);

    auto gs_x_xxyz_yz = pbuffer.data(idx_g_gd + 28);

    auto gs_x_xxyz_zz = pbuffer.data(idx_g_gd + 29);

    #pragma omp simd aligned(gc_x, gs_x_xxyz_xx, gs_x_xxyz_xy, gs_x_xxyz_xz, gs_x_xxyz_yy, gs_x_xxyz_yz, gs_x_xxyz_zz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_yy, ts_xyz_yz, ts_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyz_xx[i] = 4.0 * ts_xyz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xy[i] = 4.0 * ts_xyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_xz[i] = 4.0 * ts_xyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yy[i] = 4.0 * ts_xyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxyz_yz[i] = 4.0 * ts_xyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxyz_zz[i] = 4.0 * ts_xyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 30-36 components of targeted buffer : GD

    auto gs_x_xxzz_xx = pbuffer.data(idx_g_gd + 30);

    auto gs_x_xxzz_xy = pbuffer.data(idx_g_gd + 31);

    auto gs_x_xxzz_xz = pbuffer.data(idx_g_gd + 32);

    auto gs_x_xxzz_yy = pbuffer.data(idx_g_gd + 33);

    auto gs_x_xxzz_yz = pbuffer.data(idx_g_gd + 34);

    auto gs_x_xxzz_zz = pbuffer.data(idx_g_gd + 35);

    #pragma omp simd aligned(gc_x, gs_x_xxzz_xx, gs_x_xxzz_xy, gs_x_xxzz_xz, gs_x_xxzz_yy, gs_x_xxzz_yz, gs_x_xxzz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzz_xx[i] = 4.0 * ts_xzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xy[i] = 4.0 * ts_xzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_xz[i] = 4.0 * ts_xzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yy[i] = 4.0 * ts_xzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xxzz_yz[i] = 4.0 * ts_xzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xxzz_zz[i] = 4.0 * ts_xzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 36-42 components of targeted buffer : GD

    auto gs_x_xyyy_xx = pbuffer.data(idx_g_gd + 36);

    auto gs_x_xyyy_xy = pbuffer.data(idx_g_gd + 37);

    auto gs_x_xyyy_xz = pbuffer.data(idx_g_gd + 38);

    auto gs_x_xyyy_yy = pbuffer.data(idx_g_gd + 39);

    auto gs_x_xyyy_yz = pbuffer.data(idx_g_gd + 40);

    auto gs_x_xyyy_zz = pbuffer.data(idx_g_gd + 41);

    #pragma omp simd aligned(gc_x, gs_x_xyyy_xx, gs_x_xyyy_xy, gs_x_xyyy_xz, gs_x_xyyy_yy, gs_x_xyyy_yz, gs_x_xyyy_zz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyy_xx[i] = 2.0 * ts_yyy_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xx[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xy[i] = 2.0 * ts_yyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_xz[i] = 2.0 * ts_yyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yy[i] = 2.0 * ts_yyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yy[i] * gc_x[i] * tce_0;

        gs_x_xyyy_yz[i] = 2.0 * ts_yyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yz[i] * gc_x[i] * tce_0;

        gs_x_xyyy_zz[i] = 2.0 * ts_yyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 42-48 components of targeted buffer : GD

    auto gs_x_xyyz_xx = pbuffer.data(idx_g_gd + 42);

    auto gs_x_xyyz_xy = pbuffer.data(idx_g_gd + 43);

    auto gs_x_xyyz_xz = pbuffer.data(idx_g_gd + 44);

    auto gs_x_xyyz_yy = pbuffer.data(idx_g_gd + 45);

    auto gs_x_xyyz_yz = pbuffer.data(idx_g_gd + 46);

    auto gs_x_xyyz_zz = pbuffer.data(idx_g_gd + 47);

    #pragma omp simd aligned(gc_x, gs_x_xyyz_xx, gs_x_xyyz_xy, gs_x_xyyz_xz, gs_x_xyyz_yy, gs_x_xyyz_yz, gs_x_xyyz_zz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_yy, ts_yyz_yz, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyz_xx[i] = 2.0 * ts_yyz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xx[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xy[i] = 2.0 * ts_yyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_xz[i] = 2.0 * ts_yyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yy[i] = 2.0 * ts_yyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yy[i] * gc_x[i] * tce_0;

        gs_x_xyyz_yz[i] = 2.0 * ts_yyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yz[i] * gc_x[i] * tce_0;

        gs_x_xyyz_zz[i] = 2.0 * ts_yyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 48-54 components of targeted buffer : GD

    auto gs_x_xyzz_xx = pbuffer.data(idx_g_gd + 48);

    auto gs_x_xyzz_xy = pbuffer.data(idx_g_gd + 49);

    auto gs_x_xyzz_xz = pbuffer.data(idx_g_gd + 50);

    auto gs_x_xyzz_yy = pbuffer.data(idx_g_gd + 51);

    auto gs_x_xyzz_yz = pbuffer.data(idx_g_gd + 52);

    auto gs_x_xyzz_zz = pbuffer.data(idx_g_gd + 53);

    #pragma omp simd aligned(gc_x, gs_x_xyzz_xx, gs_x_xyzz_xy, gs_x_xyzz_xz, gs_x_xyzz_yy, gs_x_xyzz_yz, gs_x_xyzz_zz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_yy, ts_yzz_yz, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzz_xx[i] = 2.0 * ts_yzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xy[i] = 2.0 * ts_yzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_xz[i] = 2.0 * ts_yzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yy[i] = 2.0 * ts_yzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xyzz_yz[i] = 2.0 * ts_yzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xyzz_zz[i] = 2.0 * ts_yzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 54-60 components of targeted buffer : GD

    auto gs_x_xzzz_xx = pbuffer.data(idx_g_gd + 54);

    auto gs_x_xzzz_xy = pbuffer.data(idx_g_gd + 55);

    auto gs_x_xzzz_xz = pbuffer.data(idx_g_gd + 56);

    auto gs_x_xzzz_yy = pbuffer.data(idx_g_gd + 57);

    auto gs_x_xzzz_yz = pbuffer.data(idx_g_gd + 58);

    auto gs_x_xzzz_zz = pbuffer.data(idx_g_gd + 59);

    #pragma omp simd aligned(gc_x, gs_x_xzzz_xx, gs_x_xzzz_xy, gs_x_xzzz_xz, gs_x_xzzz_yy, gs_x_xzzz_yz, gs_x_xzzz_zz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzz_xx[i] = 2.0 * ts_zzz_xx[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xy[i] = 2.0 * ts_zzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_xz[i] = 2.0 * ts_zzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yy[i] = 2.0 * ts_zzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_xzzz_yz[i] = 2.0 * ts_zzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_xzzz_zz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 60-66 components of targeted buffer : GD

    auto gs_x_yyyy_xx = pbuffer.data(idx_g_gd + 60);

    auto gs_x_yyyy_xy = pbuffer.data(idx_g_gd + 61);

    auto gs_x_yyyy_xz = pbuffer.data(idx_g_gd + 62);

    auto gs_x_yyyy_yy = pbuffer.data(idx_g_gd + 63);

    auto gs_x_yyyy_yz = pbuffer.data(idx_g_gd + 64);

    auto gs_x_yyyy_zz = pbuffer.data(idx_g_gd + 65);

    #pragma omp simd aligned(gc_x, gs_x_yyyy_xx, gs_x_yyyy_xy, gs_x_yyyy_xz, gs_x_yyyy_yy, gs_x_yyyy_yz, gs_x_yyyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyy_xx[i] = 4.0 * ts_yyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xx[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xy[i] = 2.0 * ts_yyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_xz[i] = 2.0 * ts_yyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yy[i] = 2.0 * ts_yyyy_yy[i] * gc_x[i] * tce_0;

        gs_x_yyyy_yz[i] = 2.0 * ts_yyyy_yz[i] * gc_x[i] * tce_0;

        gs_x_yyyy_zz[i] = 2.0 * ts_yyyy_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 66-72 components of targeted buffer : GD

    auto gs_x_yyyz_xx = pbuffer.data(idx_g_gd + 66);

    auto gs_x_yyyz_xy = pbuffer.data(idx_g_gd + 67);

    auto gs_x_yyyz_xz = pbuffer.data(idx_g_gd + 68);

    auto gs_x_yyyz_yy = pbuffer.data(idx_g_gd + 69);

    auto gs_x_yyyz_yz = pbuffer.data(idx_g_gd + 70);

    auto gs_x_yyyz_zz = pbuffer.data(idx_g_gd + 71);

    #pragma omp simd aligned(gc_x, gs_x_yyyz_xx, gs_x_yyyz_xy, gs_x_yyyz_xz, gs_x_yyyz_yy, gs_x_yyyz_yz, gs_x_yyyz_zz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyz_xx[i] = 4.0 * ts_yyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xx[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xy[i] = 2.0 * ts_yyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_xz[i] = 2.0 * ts_yyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yy[i] = 2.0 * ts_yyyz_yy[i] * gc_x[i] * tce_0;

        gs_x_yyyz_yz[i] = 2.0 * ts_yyyz_yz[i] * gc_x[i] * tce_0;

        gs_x_yyyz_zz[i] = 2.0 * ts_yyyz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 72-78 components of targeted buffer : GD

    auto gs_x_yyzz_xx = pbuffer.data(idx_g_gd + 72);

    auto gs_x_yyzz_xy = pbuffer.data(idx_g_gd + 73);

    auto gs_x_yyzz_xz = pbuffer.data(idx_g_gd + 74);

    auto gs_x_yyzz_yy = pbuffer.data(idx_g_gd + 75);

    auto gs_x_yyzz_yz = pbuffer.data(idx_g_gd + 76);

    auto gs_x_yyzz_zz = pbuffer.data(idx_g_gd + 77);

    #pragma omp simd aligned(gc_x, gs_x_yyzz_xx, gs_x_yyzz_xy, gs_x_yyzz_xz, gs_x_yyzz_yy, gs_x_yyzz_yz, gs_x_yyzz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzz_xx[i] = 4.0 * ts_yyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xx[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xy[i] = 2.0 * ts_yyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_xz[i] = 2.0 * ts_yyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yy[i] = 2.0 * ts_yyzz_yy[i] * gc_x[i] * tce_0;

        gs_x_yyzz_yz[i] = 2.0 * ts_yyzz_yz[i] * gc_x[i] * tce_0;

        gs_x_yyzz_zz[i] = 2.0 * ts_yyzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 78-84 components of targeted buffer : GD

    auto gs_x_yzzz_xx = pbuffer.data(idx_g_gd + 78);

    auto gs_x_yzzz_xy = pbuffer.data(idx_g_gd + 79);

    auto gs_x_yzzz_xz = pbuffer.data(idx_g_gd + 80);

    auto gs_x_yzzz_yy = pbuffer.data(idx_g_gd + 81);

    auto gs_x_yzzz_yz = pbuffer.data(idx_g_gd + 82);

    auto gs_x_yzzz_zz = pbuffer.data(idx_g_gd + 83);

    #pragma omp simd aligned(gc_x, gs_x_yzzz_xx, gs_x_yzzz_xy, gs_x_yzzz_xz, gs_x_yzzz_yy, gs_x_yzzz_yz, gs_x_yzzz_zz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzz_xx[i] = 4.0 * ts_yzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xy[i] = 2.0 * ts_yzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_xz[i] = 2.0 * ts_yzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yy[i] = 2.0 * ts_yzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_yzzz_yz[i] = 2.0 * ts_yzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_yzzz_zz[i] = 2.0 * ts_yzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 84-90 components of targeted buffer : GD

    auto gs_x_zzzz_xx = pbuffer.data(idx_g_gd + 84);

    auto gs_x_zzzz_xy = pbuffer.data(idx_g_gd + 85);

    auto gs_x_zzzz_xz = pbuffer.data(idx_g_gd + 86);

    auto gs_x_zzzz_yy = pbuffer.data(idx_g_gd + 87);

    auto gs_x_zzzz_yz = pbuffer.data(idx_g_gd + 88);

    auto gs_x_zzzz_zz = pbuffer.data(idx_g_gd + 89);

    #pragma omp simd aligned(gc_x, gs_x_zzzz_xx, gs_x_zzzz_xy, gs_x_zzzz_xz, gs_x_zzzz_yy, gs_x_zzzz_yz, gs_x_zzzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzz_xx[i] = 4.0 * ts_zzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xx[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xy[i] = 2.0 * ts_zzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_xz[i] = 2.0 * ts_zzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yy[i] = 2.0 * ts_zzzz_yy[i] * gc_x[i] * tce_0;

        gs_x_zzzz_yz[i] = 2.0 * ts_zzzz_yz[i] * gc_x[i] * tce_0;

        gs_x_zzzz_zz[i] = 2.0 * ts_zzzz_zz[i] * gc_x[i] * tce_0;
    }

    // Set up 90-96 components of targeted buffer : GD

    auto gs_y_xxxx_xx = pbuffer.data(idx_g_gd + 90);

    auto gs_y_xxxx_xy = pbuffer.data(idx_g_gd + 91);

    auto gs_y_xxxx_xz = pbuffer.data(idx_g_gd + 92);

    auto gs_y_xxxx_yy = pbuffer.data(idx_g_gd + 93);

    auto gs_y_xxxx_yz = pbuffer.data(idx_g_gd + 94);

    auto gs_y_xxxx_zz = pbuffer.data(idx_g_gd + 95);

    #pragma omp simd aligned(gc_y, gs_y_xxxx_xx, gs_y_xxxx_xy, gs_y_xxxx_xz, gs_y_xxxx_yy, gs_y_xxxx_yz, gs_y_xxxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxx_xx[i] = 2.0 * ts_xxxx_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xy[i] = 2.0 * ts_xxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_xz[i] = 2.0 * ts_xxxx_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yy[i] = 4.0 * ts_xxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxx_yz[i] = 2.0 * ts_xxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxx_zz[i] = 2.0 * ts_xxxx_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 96-102 components of targeted buffer : GD

    auto gs_y_xxxy_xx = pbuffer.data(idx_g_gd + 96);

    auto gs_y_xxxy_xy = pbuffer.data(idx_g_gd + 97);

    auto gs_y_xxxy_xz = pbuffer.data(idx_g_gd + 98);

    auto gs_y_xxxy_yy = pbuffer.data(idx_g_gd + 99);

    auto gs_y_xxxy_yz = pbuffer.data(idx_g_gd + 100);

    auto gs_y_xxxy_zz = pbuffer.data(idx_g_gd + 101);

    #pragma omp simd aligned(gc_y, gs_y_xxxy_xx, gs_y_xxxy_xy, gs_y_xxxy_xz, gs_y_xxxy_yy, gs_y_xxxy_yz, gs_y_xxxy_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxy_xx[i] = 2.0 * ts_xxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xy[i] = 2.0 * ts_xxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_xz[i] = 2.0 * ts_xxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yy[i] = 2.0 * ts_xxx_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxy_yz[i] = 2.0 * ts_xxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxy_zz[i] = 2.0 * ts_xxx_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 102-108 components of targeted buffer : GD

    auto gs_y_xxxz_xx = pbuffer.data(idx_g_gd + 102);

    auto gs_y_xxxz_xy = pbuffer.data(idx_g_gd + 103);

    auto gs_y_xxxz_xz = pbuffer.data(idx_g_gd + 104);

    auto gs_y_xxxz_yy = pbuffer.data(idx_g_gd + 105);

    auto gs_y_xxxz_yz = pbuffer.data(idx_g_gd + 106);

    auto gs_y_xxxz_zz = pbuffer.data(idx_g_gd + 107);

    #pragma omp simd aligned(gc_y, gs_y_xxxz_xx, gs_y_xxxz_xy, gs_y_xxxz_xz, gs_y_xxxz_yy, gs_y_xxxz_yz, gs_y_xxxz_zz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxz_xx[i] = 2.0 * ts_xxxz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xy[i] = 2.0 * ts_xxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_xz[i] = 2.0 * ts_xxxz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yy[i] = 4.0 * ts_xxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxxz_yz[i] = 2.0 * ts_xxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxxz_zz[i] = 2.0 * ts_xxxz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 108-114 components of targeted buffer : GD

    auto gs_y_xxyy_xx = pbuffer.data(idx_g_gd + 108);

    auto gs_y_xxyy_xy = pbuffer.data(idx_g_gd + 109);

    auto gs_y_xxyy_xz = pbuffer.data(idx_g_gd + 110);

    auto gs_y_xxyy_yy = pbuffer.data(idx_g_gd + 111);

    auto gs_y_xxyy_yz = pbuffer.data(idx_g_gd + 112);

    auto gs_y_xxyy_zz = pbuffer.data(idx_g_gd + 113);

    #pragma omp simd aligned(gc_y, gs_y_xxyy_xx, gs_y_xxyy_xy, gs_y_xxyy_xz, gs_y_xxyy_yy, gs_y_xxyy_yz, gs_y_xxyy_zz, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_yy, ts_xxy_yz, ts_xxy_zz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyy_xx[i] = 4.0 * ts_xxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xx[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xy[i] = 4.0 * ts_xxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_xz[i] = 4.0 * ts_xxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yy[i] = 4.0 * ts_xxy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yy[i] * gc_y[i] * tce_0;

        gs_y_xxyy_yz[i] = 4.0 * ts_xxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yz[i] * gc_y[i] * tce_0;

        gs_y_xxyy_zz[i] = 4.0 * ts_xxy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 114-120 components of targeted buffer : GD

    auto gs_y_xxyz_xx = pbuffer.data(idx_g_gd + 114);

    auto gs_y_xxyz_xy = pbuffer.data(idx_g_gd + 115);

    auto gs_y_xxyz_xz = pbuffer.data(idx_g_gd + 116);

    auto gs_y_xxyz_yy = pbuffer.data(idx_g_gd + 117);

    auto gs_y_xxyz_yz = pbuffer.data(idx_g_gd + 118);

    auto gs_y_xxyz_zz = pbuffer.data(idx_g_gd + 119);

    #pragma omp simd aligned(gc_y, gs_y_xxyz_xx, gs_y_xxyz_xy, gs_y_xxyz_xz, gs_y_xxyz_yy, gs_y_xxyz_yz, gs_y_xxyz_zz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_yy, ts_xxz_yz, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyz_xx[i] = 2.0 * ts_xxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xy[i] = 2.0 * ts_xxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_xz[i] = 2.0 * ts_xxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yy[i] = 2.0 * ts_xxz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxyz_yz[i] = 2.0 * ts_xxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxyz_zz[i] = 2.0 * ts_xxz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 120-126 components of targeted buffer : GD

    auto gs_y_xxzz_xx = pbuffer.data(idx_g_gd + 120);

    auto gs_y_xxzz_xy = pbuffer.data(idx_g_gd + 121);

    auto gs_y_xxzz_xz = pbuffer.data(idx_g_gd + 122);

    auto gs_y_xxzz_yy = pbuffer.data(idx_g_gd + 123);

    auto gs_y_xxzz_yz = pbuffer.data(idx_g_gd + 124);

    auto gs_y_xxzz_zz = pbuffer.data(idx_g_gd + 125);

    #pragma omp simd aligned(gc_y, gs_y_xxzz_xx, gs_y_xxzz_xy, gs_y_xxzz_xz, gs_y_xxzz_yy, gs_y_xxzz_yz, gs_y_xxzz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzz_xx[i] = 2.0 * ts_xxzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xy[i] = 2.0 * ts_xxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_xz[i] = 2.0 * ts_xxzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yy[i] = 4.0 * ts_xxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xxzz_yz[i] = 2.0 * ts_xxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xxzz_zz[i] = 2.0 * ts_xxzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 126-132 components of targeted buffer : GD

    auto gs_y_xyyy_xx = pbuffer.data(idx_g_gd + 126);

    auto gs_y_xyyy_xy = pbuffer.data(idx_g_gd + 127);

    auto gs_y_xyyy_xz = pbuffer.data(idx_g_gd + 128);

    auto gs_y_xyyy_yy = pbuffer.data(idx_g_gd + 129);

    auto gs_y_xyyy_yz = pbuffer.data(idx_g_gd + 130);

    auto gs_y_xyyy_zz = pbuffer.data(idx_g_gd + 131);

    #pragma omp simd aligned(gc_y, gs_y_xyyy_xx, gs_y_xyyy_xy, gs_y_xyyy_xz, gs_y_xyyy_yy, gs_y_xyyy_yz, gs_y_xyyy_zz, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyy_xx[i] = 6.0 * ts_xyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xx[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xy[i] = 6.0 * ts_xyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_xz[i] = 6.0 * ts_xyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yy[i] = 6.0 * ts_xyy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yy[i] * gc_y[i] * tce_0;

        gs_y_xyyy_yz[i] = 6.0 * ts_xyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yz[i] * gc_y[i] * tce_0;

        gs_y_xyyy_zz[i] = 6.0 * ts_xyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 132-138 components of targeted buffer : GD

    auto gs_y_xyyz_xx = pbuffer.data(idx_g_gd + 132);

    auto gs_y_xyyz_xy = pbuffer.data(idx_g_gd + 133);

    auto gs_y_xyyz_xz = pbuffer.data(idx_g_gd + 134);

    auto gs_y_xyyz_yy = pbuffer.data(idx_g_gd + 135);

    auto gs_y_xyyz_yz = pbuffer.data(idx_g_gd + 136);

    auto gs_y_xyyz_zz = pbuffer.data(idx_g_gd + 137);

    #pragma omp simd aligned(gc_y, gs_y_xyyz_xx, gs_y_xyyz_xy, gs_y_xyyz_xz, gs_y_xyyz_yy, gs_y_xyyz_yz, gs_y_xyyz_zz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_yy, ts_xyz_yz, ts_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyz_xx[i] = 4.0 * ts_xyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xx[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xy[i] = 4.0 * ts_xyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_xz[i] = 4.0 * ts_xyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yy[i] = 4.0 * ts_xyz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yy[i] * gc_y[i] * tce_0;

        gs_y_xyyz_yz[i] = 4.0 * ts_xyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yz[i] * gc_y[i] * tce_0;

        gs_y_xyyz_zz[i] = 4.0 * ts_xyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 138-144 components of targeted buffer : GD

    auto gs_y_xyzz_xx = pbuffer.data(idx_g_gd + 138);

    auto gs_y_xyzz_xy = pbuffer.data(idx_g_gd + 139);

    auto gs_y_xyzz_xz = pbuffer.data(idx_g_gd + 140);

    auto gs_y_xyzz_yy = pbuffer.data(idx_g_gd + 141);

    auto gs_y_xyzz_yz = pbuffer.data(idx_g_gd + 142);

    auto gs_y_xyzz_zz = pbuffer.data(idx_g_gd + 143);

    #pragma omp simd aligned(gc_y, gs_y_xyzz_xx, gs_y_xyzz_xy, gs_y_xyzz_xz, gs_y_xyzz_yy, gs_y_xyzz_yz, gs_y_xyzz_zz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzz_xx[i] = 2.0 * ts_xzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xy[i] = 2.0 * ts_xzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_xz[i] = 2.0 * ts_xzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yy[i] = 2.0 * ts_xzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xyzz_yz[i] = 2.0 * ts_xzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xyzz_zz[i] = 2.0 * ts_xzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 144-150 components of targeted buffer : GD

    auto gs_y_xzzz_xx = pbuffer.data(idx_g_gd + 144);

    auto gs_y_xzzz_xy = pbuffer.data(idx_g_gd + 145);

    auto gs_y_xzzz_xz = pbuffer.data(idx_g_gd + 146);

    auto gs_y_xzzz_yy = pbuffer.data(idx_g_gd + 147);

    auto gs_y_xzzz_yz = pbuffer.data(idx_g_gd + 148);

    auto gs_y_xzzz_zz = pbuffer.data(idx_g_gd + 149);

    #pragma omp simd aligned(gc_y, gs_y_xzzz_xx, gs_y_xzzz_xy, gs_y_xzzz_xz, gs_y_xzzz_yy, gs_y_xzzz_yz, gs_y_xzzz_zz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzz_xx[i] = 2.0 * ts_xzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xy[i] = 2.0 * ts_xzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_xz[i] = 2.0 * ts_xzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yy[i] = 4.0 * ts_xzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_xzzz_yz[i] = 2.0 * ts_xzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_xzzz_zz[i] = 2.0 * ts_xzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 150-156 components of targeted buffer : GD

    auto gs_y_yyyy_xx = pbuffer.data(idx_g_gd + 150);

    auto gs_y_yyyy_xy = pbuffer.data(idx_g_gd + 151);

    auto gs_y_yyyy_xz = pbuffer.data(idx_g_gd + 152);

    auto gs_y_yyyy_yy = pbuffer.data(idx_g_gd + 153);

    auto gs_y_yyyy_yz = pbuffer.data(idx_g_gd + 154);

    auto gs_y_yyyy_zz = pbuffer.data(idx_g_gd + 155);

    #pragma omp simd aligned(gc_y, gs_y_yyyy_xx, gs_y_yyyy_xy, gs_y_yyyy_xz, gs_y_yyyy_yy, gs_y_yyyy_yz, gs_y_yyyy_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyy_xx[i] = 8.0 * ts_yyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xx[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xy[i] = 8.0 * ts_yyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_xz[i] = 8.0 * ts_yyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yy[i] = 8.0 * ts_yyy_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yy[i] * gc_y[i] * tce_0;

        gs_y_yyyy_yz[i] = 8.0 * ts_yyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yz[i] * gc_y[i] * tce_0;

        gs_y_yyyy_zz[i] = 8.0 * ts_yyy_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 156-162 components of targeted buffer : GD

    auto gs_y_yyyz_xx = pbuffer.data(idx_g_gd + 156);

    auto gs_y_yyyz_xy = pbuffer.data(idx_g_gd + 157);

    auto gs_y_yyyz_xz = pbuffer.data(idx_g_gd + 158);

    auto gs_y_yyyz_yy = pbuffer.data(idx_g_gd + 159);

    auto gs_y_yyyz_yz = pbuffer.data(idx_g_gd + 160);

    auto gs_y_yyyz_zz = pbuffer.data(idx_g_gd + 161);

    #pragma omp simd aligned(gc_y, gs_y_yyyz_xx, gs_y_yyyz_xy, gs_y_yyyz_xz, gs_y_yyyz_yy, gs_y_yyyz_yz, gs_y_yyyz_zz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_yy, ts_yyz_yz, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyz_xx[i] = 6.0 * ts_yyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xx[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xy[i] = 6.0 * ts_yyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_xz[i] = 6.0 * ts_yyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yy[i] = 6.0 * ts_yyz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yy[i] * gc_y[i] * tce_0;

        gs_y_yyyz_yz[i] = 6.0 * ts_yyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yz[i] * gc_y[i] * tce_0;

        gs_y_yyyz_zz[i] = 6.0 * ts_yyz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 162-168 components of targeted buffer : GD

    auto gs_y_yyzz_xx = pbuffer.data(idx_g_gd + 162);

    auto gs_y_yyzz_xy = pbuffer.data(idx_g_gd + 163);

    auto gs_y_yyzz_xz = pbuffer.data(idx_g_gd + 164);

    auto gs_y_yyzz_yy = pbuffer.data(idx_g_gd + 165);

    auto gs_y_yyzz_yz = pbuffer.data(idx_g_gd + 166);

    auto gs_y_yyzz_zz = pbuffer.data(idx_g_gd + 167);

    #pragma omp simd aligned(gc_y, gs_y_yyzz_xx, gs_y_yyzz_xy, gs_y_yyzz_xz, gs_y_yyzz_yy, gs_y_yyzz_yz, gs_y_yyzz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_yy, ts_yzz_yz, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzz_xx[i] = 4.0 * ts_yzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xx[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xy[i] = 4.0 * ts_yzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_xz[i] = 4.0 * ts_yzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yy[i] = 4.0 * ts_yzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yy[i] * gc_y[i] * tce_0;

        gs_y_yyzz_yz[i] = 4.0 * ts_yzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yz[i] * gc_y[i] * tce_0;

        gs_y_yyzz_zz[i] = 4.0 * ts_yzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 168-174 components of targeted buffer : GD

    auto gs_y_yzzz_xx = pbuffer.data(idx_g_gd + 168);

    auto gs_y_yzzz_xy = pbuffer.data(idx_g_gd + 169);

    auto gs_y_yzzz_xz = pbuffer.data(idx_g_gd + 170);

    auto gs_y_yzzz_yy = pbuffer.data(idx_g_gd + 171);

    auto gs_y_yzzz_yz = pbuffer.data(idx_g_gd + 172);

    auto gs_y_yzzz_zz = pbuffer.data(idx_g_gd + 173);

    #pragma omp simd aligned(gc_y, gs_y_yzzz_xx, gs_y_yzzz_xy, gs_y_yzzz_xz, gs_y_yzzz_yy, gs_y_yzzz_yz, gs_y_yzzz_zz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzz_xx[i] = 2.0 * ts_zzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xy[i] = 2.0 * ts_zzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_xz[i] = 2.0 * ts_zzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yy[i] = 2.0 * ts_zzz_yy[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_yzzz_yz[i] = 2.0 * ts_zzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_yzzz_zz[i] = 2.0 * ts_zzz_zz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 174-180 components of targeted buffer : GD

    auto gs_y_zzzz_xx = pbuffer.data(idx_g_gd + 174);

    auto gs_y_zzzz_xy = pbuffer.data(idx_g_gd + 175);

    auto gs_y_zzzz_xz = pbuffer.data(idx_g_gd + 176);

    auto gs_y_zzzz_yy = pbuffer.data(idx_g_gd + 177);

    auto gs_y_zzzz_yz = pbuffer.data(idx_g_gd + 178);

    auto gs_y_zzzz_zz = pbuffer.data(idx_g_gd + 179);

    #pragma omp simd aligned(gc_y, gs_y_zzzz_xx, gs_y_zzzz_xy, gs_y_zzzz_xz, gs_y_zzzz_yy, gs_y_zzzz_yz, gs_y_zzzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzz_xx[i] = 2.0 * ts_zzzz_xx[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xy[i] = 2.0 * ts_zzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_xz[i] = 2.0 * ts_zzzz_xz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yy[i] = 4.0 * ts_zzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yy[i] * gc_y[i] * tce_0;

        gs_y_zzzz_yz[i] = 2.0 * ts_zzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yz[i] * gc_y[i] * tce_0;

        gs_y_zzzz_zz[i] = 2.0 * ts_zzzz_zz[i] * gc_y[i] * tce_0;
    }

    // Set up 180-186 components of targeted buffer : GD

    auto gs_z_xxxx_xx = pbuffer.data(idx_g_gd + 180);

    auto gs_z_xxxx_xy = pbuffer.data(idx_g_gd + 181);

    auto gs_z_xxxx_xz = pbuffer.data(idx_g_gd + 182);

    auto gs_z_xxxx_yy = pbuffer.data(idx_g_gd + 183);

    auto gs_z_xxxx_yz = pbuffer.data(idx_g_gd + 184);

    auto gs_z_xxxx_zz = pbuffer.data(idx_g_gd + 185);

    #pragma omp simd aligned(gc_z, gs_z_xxxx_xx, gs_z_xxxx_xy, gs_z_xxxx_xz, gs_z_xxxx_yy, gs_z_xxxx_yz, gs_z_xxxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxx_xx[i] = 2.0 * ts_xxxx_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xy[i] = 2.0 * ts_xxxx_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_xz[i] = 2.0 * ts_xxxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yy[i] = 2.0 * ts_xxxx_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxx_yz[i] = 2.0 * ts_xxxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxx_zz[i] = 4.0 * ts_xxxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 186-192 components of targeted buffer : GD

    auto gs_z_xxxy_xx = pbuffer.data(idx_g_gd + 186);

    auto gs_z_xxxy_xy = pbuffer.data(idx_g_gd + 187);

    auto gs_z_xxxy_xz = pbuffer.data(idx_g_gd + 188);

    auto gs_z_xxxy_yy = pbuffer.data(idx_g_gd + 189);

    auto gs_z_xxxy_yz = pbuffer.data(idx_g_gd + 190);

    auto gs_z_xxxy_zz = pbuffer.data(idx_g_gd + 191);

    #pragma omp simd aligned(gc_z, gs_z_xxxy_xx, gs_z_xxxy_xy, gs_z_xxxy_xz, gs_z_xxxy_yy, gs_z_xxxy_yz, gs_z_xxxy_zz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxy_xx[i] = 2.0 * ts_xxxy_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xy[i] = 2.0 * ts_xxxy_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_xz[i] = 2.0 * ts_xxxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yy[i] = 2.0 * ts_xxxy_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxy_yz[i] = 2.0 * ts_xxxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxy_zz[i] = 4.0 * ts_xxxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 192-198 components of targeted buffer : GD

    auto gs_z_xxxz_xx = pbuffer.data(idx_g_gd + 192);

    auto gs_z_xxxz_xy = pbuffer.data(idx_g_gd + 193);

    auto gs_z_xxxz_xz = pbuffer.data(idx_g_gd + 194);

    auto gs_z_xxxz_yy = pbuffer.data(idx_g_gd + 195);

    auto gs_z_xxxz_yz = pbuffer.data(idx_g_gd + 196);

    auto gs_z_xxxz_zz = pbuffer.data(idx_g_gd + 197);

    #pragma omp simd aligned(gc_z, gs_z_xxxz_xx, gs_z_xxxz_xy, gs_z_xxxz_xz, gs_z_xxxz_yy, gs_z_xxxz_yz, gs_z_xxxz_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxz_xx[i] = 2.0 * ts_xxx_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xy[i] = 2.0 * ts_xxx_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_xz[i] = 2.0 * ts_xxx_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yy[i] = 2.0 * ts_xxx_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxxz_yz[i] = 2.0 * ts_xxx_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxxz_zz[i] = 2.0 * ts_xxx_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 198-204 components of targeted buffer : GD

    auto gs_z_xxyy_xx = pbuffer.data(idx_g_gd + 198);

    auto gs_z_xxyy_xy = pbuffer.data(idx_g_gd + 199);

    auto gs_z_xxyy_xz = pbuffer.data(idx_g_gd + 200);

    auto gs_z_xxyy_yy = pbuffer.data(idx_g_gd + 201);

    auto gs_z_xxyy_yz = pbuffer.data(idx_g_gd + 202);

    auto gs_z_xxyy_zz = pbuffer.data(idx_g_gd + 203);

    #pragma omp simd aligned(gc_z, gs_z_xxyy_xx, gs_z_xxyy_xy, gs_z_xxyy_xz, gs_z_xxyy_yy, gs_z_xxyy_yz, gs_z_xxyy_zz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyy_xx[i] = 2.0 * ts_xxyy_xx[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xy[i] = 2.0 * ts_xxyy_xy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_xz[i] = 2.0 * ts_xxyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_xz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yy[i] = 2.0 * ts_xxyy_yy[i] * gc_z[i] * tce_0;

        gs_z_xxyy_yz[i] = 2.0 * ts_xxyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_yz[i] * gc_z[i] * tce_0;

        gs_z_xxyy_zz[i] = 4.0 * ts_xxyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 204-210 components of targeted buffer : GD

    auto gs_z_xxyz_xx = pbuffer.data(idx_g_gd + 204);

    auto gs_z_xxyz_xy = pbuffer.data(idx_g_gd + 205);

    auto gs_z_xxyz_xz = pbuffer.data(idx_g_gd + 206);

    auto gs_z_xxyz_yy = pbuffer.data(idx_g_gd + 207);

    auto gs_z_xxyz_yz = pbuffer.data(idx_g_gd + 208);

    auto gs_z_xxyz_zz = pbuffer.data(idx_g_gd + 209);

    #pragma omp simd aligned(gc_z, gs_z_xxyz_xx, gs_z_xxyz_xy, gs_z_xxyz_xz, gs_z_xxyz_yy, gs_z_xxyz_yz, gs_z_xxyz_zz, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_yy, ts_xxy_yz, ts_xxy_zz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyz_xx[i] = 2.0 * ts_xxy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xy[i] = 2.0 * ts_xxy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_xz[i] = 2.0 * ts_xxy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yy[i] = 2.0 * ts_xxy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxyz_yz[i] = 2.0 * ts_xxy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxyz_zz[i] = 2.0 * ts_xxy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 210-216 components of targeted buffer : GD

    auto gs_z_xxzz_xx = pbuffer.data(idx_g_gd + 210);

    auto gs_z_xxzz_xy = pbuffer.data(idx_g_gd + 211);

    auto gs_z_xxzz_xz = pbuffer.data(idx_g_gd + 212);

    auto gs_z_xxzz_yy = pbuffer.data(idx_g_gd + 213);

    auto gs_z_xxzz_yz = pbuffer.data(idx_g_gd + 214);

    auto gs_z_xxzz_zz = pbuffer.data(idx_g_gd + 215);

    #pragma omp simd aligned(gc_z, gs_z_xxzz_xx, gs_z_xxzz_xy, gs_z_xxzz_xz, gs_z_xxzz_yy, gs_z_xxzz_yz, gs_z_xxzz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_yy, ts_xxz_yz, ts_xxz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzz_xx[i] = 4.0 * ts_xxz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xy[i] = 4.0 * ts_xxz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_xz[i] = 4.0 * ts_xxz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yy[i] = 4.0 * ts_xxz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xxzz_yz[i] = 4.0 * ts_xxz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xxzz_zz[i] = 4.0 * ts_xxz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xxzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 216-222 components of targeted buffer : GD

    auto gs_z_xyyy_xx = pbuffer.data(idx_g_gd + 216);

    auto gs_z_xyyy_xy = pbuffer.data(idx_g_gd + 217);

    auto gs_z_xyyy_xz = pbuffer.data(idx_g_gd + 218);

    auto gs_z_xyyy_yy = pbuffer.data(idx_g_gd + 219);

    auto gs_z_xyyy_yz = pbuffer.data(idx_g_gd + 220);

    auto gs_z_xyyy_zz = pbuffer.data(idx_g_gd + 221);

    #pragma omp simd aligned(gc_z, gs_z_xyyy_xx, gs_z_xyyy_xy, gs_z_xyyy_xz, gs_z_xyyy_yy, gs_z_xyyy_yz, gs_z_xyyy_zz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyy_xx[i] = 2.0 * ts_xyyy_xx[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xy[i] = 2.0 * ts_xyyy_xy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_xz[i] = 2.0 * ts_xyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_xz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yy[i] = 2.0 * ts_xyyy_yy[i] * gc_z[i] * tce_0;

        gs_z_xyyy_yz[i] = 2.0 * ts_xyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_yz[i] * gc_z[i] * tce_0;

        gs_z_xyyy_zz[i] = 4.0 * ts_xyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 222-228 components of targeted buffer : GD

    auto gs_z_xyyz_xx = pbuffer.data(idx_g_gd + 222);

    auto gs_z_xyyz_xy = pbuffer.data(idx_g_gd + 223);

    auto gs_z_xyyz_xz = pbuffer.data(idx_g_gd + 224);

    auto gs_z_xyyz_yy = pbuffer.data(idx_g_gd + 225);

    auto gs_z_xyyz_yz = pbuffer.data(idx_g_gd + 226);

    auto gs_z_xyyz_zz = pbuffer.data(idx_g_gd + 227);

    #pragma omp simd aligned(gc_z, gs_z_xyyz_xx, gs_z_xyyz_xy, gs_z_xyyz_xz, gs_z_xyyz_yy, gs_z_xyyz_yz, gs_z_xyyz_zz, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyz_xx[i] = 2.0 * ts_xyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xx[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xy[i] = 2.0 * ts_xyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_xz[i] = 2.0 * ts_xyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_xz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yy[i] = 2.0 * ts_xyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yy[i] * gc_z[i] * tce_0;

        gs_z_xyyz_yz[i] = 2.0 * ts_xyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_yz[i] * gc_z[i] * tce_0;

        gs_z_xyyz_zz[i] = 2.0 * ts_xyy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 228-234 components of targeted buffer : GD

    auto gs_z_xyzz_xx = pbuffer.data(idx_g_gd + 228);

    auto gs_z_xyzz_xy = pbuffer.data(idx_g_gd + 229);

    auto gs_z_xyzz_xz = pbuffer.data(idx_g_gd + 230);

    auto gs_z_xyzz_yy = pbuffer.data(idx_g_gd + 231);

    auto gs_z_xyzz_yz = pbuffer.data(idx_g_gd + 232);

    auto gs_z_xyzz_zz = pbuffer.data(idx_g_gd + 233);

    #pragma omp simd aligned(gc_z, gs_z_xyzz_xx, gs_z_xyzz_xy, gs_z_xyzz_xz, gs_z_xyzz_yy, gs_z_xyzz_yz, gs_z_xyzz_zz, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_yy, ts_xyz_yz, ts_xyz_zz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzz_xx[i] = 4.0 * ts_xyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xy[i] = 4.0 * ts_xyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_xz[i] = 4.0 * ts_xyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yy[i] = 4.0 * ts_xyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xyzz_yz[i] = 4.0 * ts_xyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xyzz_zz[i] = 4.0 * ts_xyz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 234-240 components of targeted buffer : GD

    auto gs_z_xzzz_xx = pbuffer.data(idx_g_gd + 234);

    auto gs_z_xzzz_xy = pbuffer.data(idx_g_gd + 235);

    auto gs_z_xzzz_xz = pbuffer.data(idx_g_gd + 236);

    auto gs_z_xzzz_yy = pbuffer.data(idx_g_gd + 237);

    auto gs_z_xzzz_yz = pbuffer.data(idx_g_gd + 238);

    auto gs_z_xzzz_zz = pbuffer.data(idx_g_gd + 239);

    #pragma omp simd aligned(gc_z, gs_z_xzzz_xx, gs_z_xzzz_xy, gs_z_xzzz_xz, gs_z_xzzz_yy, gs_z_xzzz_yz, gs_z_xzzz_zz, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzz_xx[i] = 6.0 * ts_xzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xy[i] = 6.0 * ts_xzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_xz[i] = 6.0 * ts_xzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yy[i] = 6.0 * ts_xzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_xzzz_yz[i] = 6.0 * ts_xzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_xzzz_zz[i] = 6.0 * ts_xzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_xzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 240-246 components of targeted buffer : GD

    auto gs_z_yyyy_xx = pbuffer.data(idx_g_gd + 240);

    auto gs_z_yyyy_xy = pbuffer.data(idx_g_gd + 241);

    auto gs_z_yyyy_xz = pbuffer.data(idx_g_gd + 242);

    auto gs_z_yyyy_yy = pbuffer.data(idx_g_gd + 243);

    auto gs_z_yyyy_yz = pbuffer.data(idx_g_gd + 244);

    auto gs_z_yyyy_zz = pbuffer.data(idx_g_gd + 245);

    #pragma omp simd aligned(gc_z, gs_z_yyyy_xx, gs_z_yyyy_xy, gs_z_yyyy_xz, gs_z_yyyy_yy, gs_z_yyyy_yz, gs_z_yyyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyy_xx[i] = 2.0 * ts_yyyy_xx[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xy[i] = 2.0 * ts_yyyy_xy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_xz[i] = 2.0 * ts_yyyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_xz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yy[i] = 2.0 * ts_yyyy_yy[i] * gc_z[i] * tce_0;

        gs_z_yyyy_yz[i] = 2.0 * ts_yyyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_yz[i] * gc_z[i] * tce_0;

        gs_z_yyyy_zz[i] = 4.0 * ts_yyyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 246-252 components of targeted buffer : GD

    auto gs_z_yyyz_xx = pbuffer.data(idx_g_gd + 246);

    auto gs_z_yyyz_xy = pbuffer.data(idx_g_gd + 247);

    auto gs_z_yyyz_xz = pbuffer.data(idx_g_gd + 248);

    auto gs_z_yyyz_yy = pbuffer.data(idx_g_gd + 249);

    auto gs_z_yyyz_yz = pbuffer.data(idx_g_gd + 250);

    auto gs_z_yyyz_zz = pbuffer.data(idx_g_gd + 251);

    #pragma omp simd aligned(gc_z, gs_z_yyyz_xx, gs_z_yyyz_xy, gs_z_yyyz_xz, gs_z_yyyz_yy, gs_z_yyyz_yz, gs_z_yyyz_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyz_xx[i] = 2.0 * ts_yyy_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xx[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xy[i] = 2.0 * ts_yyy_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_xz[i] = 2.0 * ts_yyy_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_xz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yy[i] = 2.0 * ts_yyy_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yy[i] * gc_z[i] * tce_0;

        gs_z_yyyz_yz[i] = 2.0 * ts_yyy_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_yz[i] * gc_z[i] * tce_0;

        gs_z_yyyz_zz[i] = 2.0 * ts_yyy_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yyyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 252-258 components of targeted buffer : GD

    auto gs_z_yyzz_xx = pbuffer.data(idx_g_gd + 252);

    auto gs_z_yyzz_xy = pbuffer.data(idx_g_gd + 253);

    auto gs_z_yyzz_xz = pbuffer.data(idx_g_gd + 254);

    auto gs_z_yyzz_yy = pbuffer.data(idx_g_gd + 255);

    auto gs_z_yyzz_yz = pbuffer.data(idx_g_gd + 256);

    auto gs_z_yyzz_zz = pbuffer.data(idx_g_gd + 257);

    #pragma omp simd aligned(gc_z, gs_z_yyzz_xx, gs_z_yyzz_xy, gs_z_yyzz_xz, gs_z_yyzz_yy, gs_z_yyzz_yz, gs_z_yyzz_zz, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_yy, ts_yyz_yz, ts_yyz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzz_xx[i] = 4.0 * ts_yyz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xx[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xy[i] = 4.0 * ts_yyz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_xz[i] = 4.0 * ts_yyz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_xz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yy[i] = 4.0 * ts_yyz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yy[i] * gc_z[i] * tce_0;

        gs_z_yyzz_yz[i] = 4.0 * ts_yyz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_yz[i] * gc_z[i] * tce_0;

        gs_z_yyzz_zz[i] = 4.0 * ts_yyz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yyzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 258-264 components of targeted buffer : GD

    auto gs_z_yzzz_xx = pbuffer.data(idx_g_gd + 258);

    auto gs_z_yzzz_xy = pbuffer.data(idx_g_gd + 259);

    auto gs_z_yzzz_xz = pbuffer.data(idx_g_gd + 260);

    auto gs_z_yzzz_yy = pbuffer.data(idx_g_gd + 261);

    auto gs_z_yzzz_yz = pbuffer.data(idx_g_gd + 262);

    auto gs_z_yzzz_zz = pbuffer.data(idx_g_gd + 263);

    #pragma omp simd aligned(gc_z, gs_z_yzzz_xx, gs_z_yzzz_xy, gs_z_yzzz_xz, gs_z_yzzz_yy, gs_z_yzzz_yz, gs_z_yzzz_zz, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_yy, ts_yzz_yz, ts_yzz_zz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzz_xx[i] = 6.0 * ts_yzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xy[i] = 6.0 * ts_yzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_xz[i] = 6.0 * ts_yzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yy[i] = 6.0 * ts_yzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_yzzz_yz[i] = 6.0 * ts_yzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_yzzz_zz[i] = 6.0 * ts_yzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_yzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_zz[i] * gc_z[i] * tce_0;
    }

    // Set up 264-270 components of targeted buffer : GD

    auto gs_z_zzzz_xx = pbuffer.data(idx_g_gd + 264);

    auto gs_z_zzzz_xy = pbuffer.data(idx_g_gd + 265);

    auto gs_z_zzzz_xz = pbuffer.data(idx_g_gd + 266);

    auto gs_z_zzzz_yy = pbuffer.data(idx_g_gd + 267);

    auto gs_z_zzzz_yz = pbuffer.data(idx_g_gd + 268);

    auto gs_z_zzzz_zz = pbuffer.data(idx_g_gd + 269);

    #pragma omp simd aligned(gc_z, gs_z_zzzz_xx, gs_z_zzzz_xy, gs_z_zzzz_xz, gs_z_zzzz_yy, gs_z_zzzz_yz, gs_z_zzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzz_xx[i] = 8.0 * ts_zzz_xx[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xx[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xy[i] = 8.0 * ts_zzz_xy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_xz[i] = 8.0 * ts_zzz_xz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_xz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yy[i] = 8.0 * ts_zzz_yy[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yy[i] * gc_z[i] * tce_0;

        gs_z_zzzz_yz[i] = 8.0 * ts_zzz_yz[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_yz[i] * gc_z[i] * tce_0;

        gs_z_zzzz_zz[i] = 8.0 * ts_zzz_zz[i] * gfe_0 * tce_0 + 4.0 * ts_zzzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_zz[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

