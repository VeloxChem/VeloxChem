#include "GeometricalDerivatives010ForFD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_fd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_fd,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_gd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : DD

    auto tr_xx_xx = pbuffer.data(idx_op_dd);

    auto tr_xx_xy = pbuffer.data(idx_op_dd + 1);

    auto tr_xx_xz = pbuffer.data(idx_op_dd + 2);

    auto tr_xx_yy = pbuffer.data(idx_op_dd + 3);

    auto tr_xx_yz = pbuffer.data(idx_op_dd + 4);

    auto tr_xx_zz = pbuffer.data(idx_op_dd + 5);

    auto tr_xy_xx = pbuffer.data(idx_op_dd + 6);

    auto tr_xy_xy = pbuffer.data(idx_op_dd + 7);

    auto tr_xy_xz = pbuffer.data(idx_op_dd + 8);

    auto tr_xy_yy = pbuffer.data(idx_op_dd + 9);

    auto tr_xy_yz = pbuffer.data(idx_op_dd + 10);

    auto tr_xy_zz = pbuffer.data(idx_op_dd + 11);

    auto tr_xz_xx = pbuffer.data(idx_op_dd + 12);

    auto tr_xz_xy = pbuffer.data(idx_op_dd + 13);

    auto tr_xz_xz = pbuffer.data(idx_op_dd + 14);

    auto tr_xz_yy = pbuffer.data(idx_op_dd + 15);

    auto tr_xz_yz = pbuffer.data(idx_op_dd + 16);

    auto tr_xz_zz = pbuffer.data(idx_op_dd + 17);

    auto tr_yy_xx = pbuffer.data(idx_op_dd + 18);

    auto tr_yy_xy = pbuffer.data(idx_op_dd + 19);

    auto tr_yy_xz = pbuffer.data(idx_op_dd + 20);

    auto tr_yy_yy = pbuffer.data(idx_op_dd + 21);

    auto tr_yy_yz = pbuffer.data(idx_op_dd + 22);

    auto tr_yy_zz = pbuffer.data(idx_op_dd + 23);

    auto tr_yz_xx = pbuffer.data(idx_op_dd + 24);

    auto tr_yz_xy = pbuffer.data(idx_op_dd + 25);

    auto tr_yz_xz = pbuffer.data(idx_op_dd + 26);

    auto tr_yz_yy = pbuffer.data(idx_op_dd + 27);

    auto tr_yz_yz = pbuffer.data(idx_op_dd + 28);

    auto tr_yz_zz = pbuffer.data(idx_op_dd + 29);

    auto tr_zz_xx = pbuffer.data(idx_op_dd + 30);

    auto tr_zz_xy = pbuffer.data(idx_op_dd + 31);

    auto tr_zz_xz = pbuffer.data(idx_op_dd + 32);

    auto tr_zz_yy = pbuffer.data(idx_op_dd + 33);

    auto tr_zz_yz = pbuffer.data(idx_op_dd + 34);

    auto tr_zz_zz = pbuffer.data(idx_op_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto tr_xxx_x = pbuffer.data(idx_op_fp);

    auto tr_xxx_y = pbuffer.data(idx_op_fp + 1);

    auto tr_xxx_z = pbuffer.data(idx_op_fp + 2);

    auto tr_xxy_x = pbuffer.data(idx_op_fp + 3);

    auto tr_xxy_y = pbuffer.data(idx_op_fp + 4);

    auto tr_xxy_z = pbuffer.data(idx_op_fp + 5);

    auto tr_xxz_x = pbuffer.data(idx_op_fp + 6);

    auto tr_xxz_y = pbuffer.data(idx_op_fp + 7);

    auto tr_xxz_z = pbuffer.data(idx_op_fp + 8);

    auto tr_xyy_x = pbuffer.data(idx_op_fp + 9);

    auto tr_xyy_y = pbuffer.data(idx_op_fp + 10);

    auto tr_xyy_z = pbuffer.data(idx_op_fp + 11);

    auto tr_xyz_x = pbuffer.data(idx_op_fp + 12);

    auto tr_xyz_y = pbuffer.data(idx_op_fp + 13);

    auto tr_xyz_z = pbuffer.data(idx_op_fp + 14);

    auto tr_xzz_x = pbuffer.data(idx_op_fp + 15);

    auto tr_xzz_y = pbuffer.data(idx_op_fp + 16);

    auto tr_xzz_z = pbuffer.data(idx_op_fp + 17);

    auto tr_yyy_x = pbuffer.data(idx_op_fp + 18);

    auto tr_yyy_y = pbuffer.data(idx_op_fp + 19);

    auto tr_yyy_z = pbuffer.data(idx_op_fp + 20);

    auto tr_yyz_x = pbuffer.data(idx_op_fp + 21);

    auto tr_yyz_y = pbuffer.data(idx_op_fp + 22);

    auto tr_yyz_z = pbuffer.data(idx_op_fp + 23);

    auto tr_yzz_x = pbuffer.data(idx_op_fp + 24);

    auto tr_yzz_y = pbuffer.data(idx_op_fp + 25);

    auto tr_yzz_z = pbuffer.data(idx_op_fp + 26);

    auto tr_zzz_x = pbuffer.data(idx_op_fp + 27);

    auto tr_zzz_y = pbuffer.data(idx_op_fp + 28);

    auto tr_zzz_z = pbuffer.data(idx_op_fp + 29);

    // Set up components of auxiliary buffer : FF

    auto tr_xxx_xxx = pbuffer.data(idx_op_ff);

    auto tr_xxx_xxy = pbuffer.data(idx_op_ff + 1);

    auto tr_xxx_xxz = pbuffer.data(idx_op_ff + 2);

    auto tr_xxx_xyy = pbuffer.data(idx_op_ff + 3);

    auto tr_xxx_xyz = pbuffer.data(idx_op_ff + 4);

    auto tr_xxx_xzz = pbuffer.data(idx_op_ff + 5);

    auto tr_xxx_yyy = pbuffer.data(idx_op_ff + 6);

    auto tr_xxx_yyz = pbuffer.data(idx_op_ff + 7);

    auto tr_xxx_yzz = pbuffer.data(idx_op_ff + 8);

    auto tr_xxx_zzz = pbuffer.data(idx_op_ff + 9);

    auto tr_xxy_xxx = pbuffer.data(idx_op_ff + 10);

    auto tr_xxy_xxy = pbuffer.data(idx_op_ff + 11);

    auto tr_xxy_xxz = pbuffer.data(idx_op_ff + 12);

    auto tr_xxy_xyy = pbuffer.data(idx_op_ff + 13);

    auto tr_xxy_xyz = pbuffer.data(idx_op_ff + 14);

    auto tr_xxy_xzz = pbuffer.data(idx_op_ff + 15);

    auto tr_xxy_yyy = pbuffer.data(idx_op_ff + 16);

    auto tr_xxy_yyz = pbuffer.data(idx_op_ff + 17);

    auto tr_xxy_yzz = pbuffer.data(idx_op_ff + 18);

    auto tr_xxy_zzz = pbuffer.data(idx_op_ff + 19);

    auto tr_xxz_xxx = pbuffer.data(idx_op_ff + 20);

    auto tr_xxz_xxy = pbuffer.data(idx_op_ff + 21);

    auto tr_xxz_xxz = pbuffer.data(idx_op_ff + 22);

    auto tr_xxz_xyy = pbuffer.data(idx_op_ff + 23);

    auto tr_xxz_xyz = pbuffer.data(idx_op_ff + 24);

    auto tr_xxz_xzz = pbuffer.data(idx_op_ff + 25);

    auto tr_xxz_yyy = pbuffer.data(idx_op_ff + 26);

    auto tr_xxz_yyz = pbuffer.data(idx_op_ff + 27);

    auto tr_xxz_yzz = pbuffer.data(idx_op_ff + 28);

    auto tr_xxz_zzz = pbuffer.data(idx_op_ff + 29);

    auto tr_xyy_xxx = pbuffer.data(idx_op_ff + 30);

    auto tr_xyy_xxy = pbuffer.data(idx_op_ff + 31);

    auto tr_xyy_xxz = pbuffer.data(idx_op_ff + 32);

    auto tr_xyy_xyy = pbuffer.data(idx_op_ff + 33);

    auto tr_xyy_xyz = pbuffer.data(idx_op_ff + 34);

    auto tr_xyy_xzz = pbuffer.data(idx_op_ff + 35);

    auto tr_xyy_yyy = pbuffer.data(idx_op_ff + 36);

    auto tr_xyy_yyz = pbuffer.data(idx_op_ff + 37);

    auto tr_xyy_yzz = pbuffer.data(idx_op_ff + 38);

    auto tr_xyy_zzz = pbuffer.data(idx_op_ff + 39);

    auto tr_xyz_xxx = pbuffer.data(idx_op_ff + 40);

    auto tr_xyz_xxy = pbuffer.data(idx_op_ff + 41);

    auto tr_xyz_xxz = pbuffer.data(idx_op_ff + 42);

    auto tr_xyz_xyy = pbuffer.data(idx_op_ff + 43);

    auto tr_xyz_xyz = pbuffer.data(idx_op_ff + 44);

    auto tr_xyz_xzz = pbuffer.data(idx_op_ff + 45);

    auto tr_xyz_yyy = pbuffer.data(idx_op_ff + 46);

    auto tr_xyz_yyz = pbuffer.data(idx_op_ff + 47);

    auto tr_xyz_yzz = pbuffer.data(idx_op_ff + 48);

    auto tr_xyz_zzz = pbuffer.data(idx_op_ff + 49);

    auto tr_xzz_xxx = pbuffer.data(idx_op_ff + 50);

    auto tr_xzz_xxy = pbuffer.data(idx_op_ff + 51);

    auto tr_xzz_xxz = pbuffer.data(idx_op_ff + 52);

    auto tr_xzz_xyy = pbuffer.data(idx_op_ff + 53);

    auto tr_xzz_xyz = pbuffer.data(idx_op_ff + 54);

    auto tr_xzz_xzz = pbuffer.data(idx_op_ff + 55);

    auto tr_xzz_yyy = pbuffer.data(idx_op_ff + 56);

    auto tr_xzz_yyz = pbuffer.data(idx_op_ff + 57);

    auto tr_xzz_yzz = pbuffer.data(idx_op_ff + 58);

    auto tr_xzz_zzz = pbuffer.data(idx_op_ff + 59);

    auto tr_yyy_xxx = pbuffer.data(idx_op_ff + 60);

    auto tr_yyy_xxy = pbuffer.data(idx_op_ff + 61);

    auto tr_yyy_xxz = pbuffer.data(idx_op_ff + 62);

    auto tr_yyy_xyy = pbuffer.data(idx_op_ff + 63);

    auto tr_yyy_xyz = pbuffer.data(idx_op_ff + 64);

    auto tr_yyy_xzz = pbuffer.data(idx_op_ff + 65);

    auto tr_yyy_yyy = pbuffer.data(idx_op_ff + 66);

    auto tr_yyy_yyz = pbuffer.data(idx_op_ff + 67);

    auto tr_yyy_yzz = pbuffer.data(idx_op_ff + 68);

    auto tr_yyy_zzz = pbuffer.data(idx_op_ff + 69);

    auto tr_yyz_xxx = pbuffer.data(idx_op_ff + 70);

    auto tr_yyz_xxy = pbuffer.data(idx_op_ff + 71);

    auto tr_yyz_xxz = pbuffer.data(idx_op_ff + 72);

    auto tr_yyz_xyy = pbuffer.data(idx_op_ff + 73);

    auto tr_yyz_xyz = pbuffer.data(idx_op_ff + 74);

    auto tr_yyz_xzz = pbuffer.data(idx_op_ff + 75);

    auto tr_yyz_yyy = pbuffer.data(idx_op_ff + 76);

    auto tr_yyz_yyz = pbuffer.data(idx_op_ff + 77);

    auto tr_yyz_yzz = pbuffer.data(idx_op_ff + 78);

    auto tr_yyz_zzz = pbuffer.data(idx_op_ff + 79);

    auto tr_yzz_xxx = pbuffer.data(idx_op_ff + 80);

    auto tr_yzz_xxy = pbuffer.data(idx_op_ff + 81);

    auto tr_yzz_xxz = pbuffer.data(idx_op_ff + 82);

    auto tr_yzz_xyy = pbuffer.data(idx_op_ff + 83);

    auto tr_yzz_xyz = pbuffer.data(idx_op_ff + 84);

    auto tr_yzz_xzz = pbuffer.data(idx_op_ff + 85);

    auto tr_yzz_yyy = pbuffer.data(idx_op_ff + 86);

    auto tr_yzz_yyz = pbuffer.data(idx_op_ff + 87);

    auto tr_yzz_yzz = pbuffer.data(idx_op_ff + 88);

    auto tr_yzz_zzz = pbuffer.data(idx_op_ff + 89);

    auto tr_zzz_xxx = pbuffer.data(idx_op_ff + 90);

    auto tr_zzz_xxy = pbuffer.data(idx_op_ff + 91);

    auto tr_zzz_xxz = pbuffer.data(idx_op_ff + 92);

    auto tr_zzz_xyy = pbuffer.data(idx_op_ff + 93);

    auto tr_zzz_xyz = pbuffer.data(idx_op_ff + 94);

    auto tr_zzz_xzz = pbuffer.data(idx_op_ff + 95);

    auto tr_zzz_yyy = pbuffer.data(idx_op_ff + 96);

    auto tr_zzz_yyz = pbuffer.data(idx_op_ff + 97);

    auto tr_zzz_yzz = pbuffer.data(idx_op_ff + 98);

    auto tr_zzz_zzz = pbuffer.data(idx_op_ff + 99);

    // Set up components of auxiliary buffer : GD

    auto tr_xxxx_xx = pbuffer.data(idx_op_gd);

    auto tr_xxxx_xy = pbuffer.data(idx_op_gd + 1);

    auto tr_xxxx_xz = pbuffer.data(idx_op_gd + 2);

    auto tr_xxxx_yy = pbuffer.data(idx_op_gd + 3);

    auto tr_xxxx_yz = pbuffer.data(idx_op_gd + 4);

    auto tr_xxxx_zz = pbuffer.data(idx_op_gd + 5);

    auto tr_xxxy_xx = pbuffer.data(idx_op_gd + 6);

    auto tr_xxxy_xy = pbuffer.data(idx_op_gd + 7);

    auto tr_xxxy_xz = pbuffer.data(idx_op_gd + 8);

    auto tr_xxxy_yy = pbuffer.data(idx_op_gd + 9);

    auto tr_xxxy_yz = pbuffer.data(idx_op_gd + 10);

    auto tr_xxxy_zz = pbuffer.data(idx_op_gd + 11);

    auto tr_xxxz_xx = pbuffer.data(idx_op_gd + 12);

    auto tr_xxxz_xy = pbuffer.data(idx_op_gd + 13);

    auto tr_xxxz_xz = pbuffer.data(idx_op_gd + 14);

    auto tr_xxxz_yy = pbuffer.data(idx_op_gd + 15);

    auto tr_xxxz_yz = pbuffer.data(idx_op_gd + 16);

    auto tr_xxxz_zz = pbuffer.data(idx_op_gd + 17);

    auto tr_xxyy_xx = pbuffer.data(idx_op_gd + 18);

    auto tr_xxyy_xy = pbuffer.data(idx_op_gd + 19);

    auto tr_xxyy_xz = pbuffer.data(idx_op_gd + 20);

    auto tr_xxyy_yy = pbuffer.data(idx_op_gd + 21);

    auto tr_xxyy_yz = pbuffer.data(idx_op_gd + 22);

    auto tr_xxyy_zz = pbuffer.data(idx_op_gd + 23);

    auto tr_xxyz_xx = pbuffer.data(idx_op_gd + 24);

    auto tr_xxyz_xy = pbuffer.data(idx_op_gd + 25);

    auto tr_xxyz_xz = pbuffer.data(idx_op_gd + 26);

    auto tr_xxyz_yy = pbuffer.data(idx_op_gd + 27);

    auto tr_xxyz_yz = pbuffer.data(idx_op_gd + 28);

    auto tr_xxyz_zz = pbuffer.data(idx_op_gd + 29);

    auto tr_xxzz_xx = pbuffer.data(idx_op_gd + 30);

    auto tr_xxzz_xy = pbuffer.data(idx_op_gd + 31);

    auto tr_xxzz_xz = pbuffer.data(idx_op_gd + 32);

    auto tr_xxzz_yy = pbuffer.data(idx_op_gd + 33);

    auto tr_xxzz_yz = pbuffer.data(idx_op_gd + 34);

    auto tr_xxzz_zz = pbuffer.data(idx_op_gd + 35);

    auto tr_xyyy_xx = pbuffer.data(idx_op_gd + 36);

    auto tr_xyyy_xy = pbuffer.data(idx_op_gd + 37);

    auto tr_xyyy_xz = pbuffer.data(idx_op_gd + 38);

    auto tr_xyyy_yy = pbuffer.data(idx_op_gd + 39);

    auto tr_xyyy_yz = pbuffer.data(idx_op_gd + 40);

    auto tr_xyyy_zz = pbuffer.data(idx_op_gd + 41);

    auto tr_xyyz_xx = pbuffer.data(idx_op_gd + 42);

    auto tr_xyyz_xy = pbuffer.data(idx_op_gd + 43);

    auto tr_xyyz_xz = pbuffer.data(idx_op_gd + 44);

    auto tr_xyyz_yy = pbuffer.data(idx_op_gd + 45);

    auto tr_xyyz_yz = pbuffer.data(idx_op_gd + 46);

    auto tr_xyyz_zz = pbuffer.data(idx_op_gd + 47);

    auto tr_xyzz_xx = pbuffer.data(idx_op_gd + 48);

    auto tr_xyzz_xy = pbuffer.data(idx_op_gd + 49);

    auto tr_xyzz_xz = pbuffer.data(idx_op_gd + 50);

    auto tr_xyzz_yy = pbuffer.data(idx_op_gd + 51);

    auto tr_xyzz_yz = pbuffer.data(idx_op_gd + 52);

    auto tr_xyzz_zz = pbuffer.data(idx_op_gd + 53);

    auto tr_xzzz_xx = pbuffer.data(idx_op_gd + 54);

    auto tr_xzzz_xy = pbuffer.data(idx_op_gd + 55);

    auto tr_xzzz_xz = pbuffer.data(idx_op_gd + 56);

    auto tr_xzzz_yy = pbuffer.data(idx_op_gd + 57);

    auto tr_xzzz_yz = pbuffer.data(idx_op_gd + 58);

    auto tr_xzzz_zz = pbuffer.data(idx_op_gd + 59);

    auto tr_yyyy_xx = pbuffer.data(idx_op_gd + 60);

    auto tr_yyyy_xy = pbuffer.data(idx_op_gd + 61);

    auto tr_yyyy_xz = pbuffer.data(idx_op_gd + 62);

    auto tr_yyyy_yy = pbuffer.data(idx_op_gd + 63);

    auto tr_yyyy_yz = pbuffer.data(idx_op_gd + 64);

    auto tr_yyyy_zz = pbuffer.data(idx_op_gd + 65);

    auto tr_yyyz_xx = pbuffer.data(idx_op_gd + 66);

    auto tr_yyyz_xy = pbuffer.data(idx_op_gd + 67);

    auto tr_yyyz_xz = pbuffer.data(idx_op_gd + 68);

    auto tr_yyyz_yy = pbuffer.data(idx_op_gd + 69);

    auto tr_yyyz_yz = pbuffer.data(idx_op_gd + 70);

    auto tr_yyyz_zz = pbuffer.data(idx_op_gd + 71);

    auto tr_yyzz_xx = pbuffer.data(idx_op_gd + 72);

    auto tr_yyzz_xy = pbuffer.data(idx_op_gd + 73);

    auto tr_yyzz_xz = pbuffer.data(idx_op_gd + 74);

    auto tr_yyzz_yy = pbuffer.data(idx_op_gd + 75);

    auto tr_yyzz_yz = pbuffer.data(idx_op_gd + 76);

    auto tr_yyzz_zz = pbuffer.data(idx_op_gd + 77);

    auto tr_yzzz_xx = pbuffer.data(idx_op_gd + 78);

    auto tr_yzzz_xy = pbuffer.data(idx_op_gd + 79);

    auto tr_yzzz_xz = pbuffer.data(idx_op_gd + 80);

    auto tr_yzzz_yy = pbuffer.data(idx_op_gd + 81);

    auto tr_yzzz_yz = pbuffer.data(idx_op_gd + 82);

    auto tr_yzzz_zz = pbuffer.data(idx_op_gd + 83);

    auto tr_zzzz_xx = pbuffer.data(idx_op_gd + 84);

    auto tr_zzzz_xy = pbuffer.data(idx_op_gd + 85);

    auto tr_zzzz_xz = pbuffer.data(idx_op_gd + 86);

    auto tr_zzzz_yy = pbuffer.data(idx_op_gd + 87);

    auto tr_zzzz_yz = pbuffer.data(idx_op_gd + 88);

    auto tr_zzzz_zz = pbuffer.data(idx_op_gd + 89);

    // Set up 0-6 components of targeted buffer : FD

    auto tr_0_0_x_xxx_xx = pbuffer.data(idx_op_geom_010_fd);

    auto tr_0_0_x_xxx_xy = pbuffer.data(idx_op_geom_010_fd + 1);

    auto tr_0_0_x_xxx_xz = pbuffer.data(idx_op_geom_010_fd + 2);

    auto tr_0_0_x_xxx_yy = pbuffer.data(idx_op_geom_010_fd + 3);

    auto tr_0_0_x_xxx_yz = pbuffer.data(idx_op_geom_010_fd + 4);

    auto tr_0_0_x_xxx_zz = pbuffer.data(idx_op_geom_010_fd + 5);

    #pragma omp simd aligned(tr_0_0_x_xxx_xx, tr_0_0_x_xxx_xy, tr_0_0_x_xxx_xz, tr_0_0_x_xxx_yy, tr_0_0_x_xxx_yz, tr_0_0_x_xxx_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_z, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxx_xx[i] = 2.0 * tr_xxxx_xx[i] * tbe_0 + 2.0 * tr_xxx_xxx[i] * tke_0 - 3.0 * tr_xx_xx[i] - 2.0 * tr_xxx_x[i];

        tr_0_0_x_xxx_xy[i] = 2.0 * tr_xxxx_xy[i] * tbe_0 + 2.0 * tr_xxx_xxy[i] * tke_0 - 3.0 * tr_xx_xy[i] - tr_xxx_y[i];

        tr_0_0_x_xxx_xz[i] = 2.0 * tr_xxxx_xz[i] * tbe_0 + 2.0 * tr_xxx_xxz[i] * tke_0 - 3.0 * tr_xx_xz[i] - tr_xxx_z[i];

        tr_0_0_x_xxx_yy[i] = 2.0 * tr_xxxx_yy[i] * tbe_0 + 2.0 * tr_xxx_xyy[i] * tke_0 - 3.0 * tr_xx_yy[i];

        tr_0_0_x_xxx_yz[i] = 2.0 * tr_xxxx_yz[i] * tbe_0 + 2.0 * tr_xxx_xyz[i] * tke_0 - 3.0 * tr_xx_yz[i];

        tr_0_0_x_xxx_zz[i] = 2.0 * tr_xxxx_zz[i] * tbe_0 + 2.0 * tr_xxx_xzz[i] * tke_0 - 3.0 * tr_xx_zz[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto tr_0_0_x_xxy_xx = pbuffer.data(idx_op_geom_010_fd + 6);

    auto tr_0_0_x_xxy_xy = pbuffer.data(idx_op_geom_010_fd + 7);

    auto tr_0_0_x_xxy_xz = pbuffer.data(idx_op_geom_010_fd + 8);

    auto tr_0_0_x_xxy_yy = pbuffer.data(idx_op_geom_010_fd + 9);

    auto tr_0_0_x_xxy_yz = pbuffer.data(idx_op_geom_010_fd + 10);

    auto tr_0_0_x_xxy_zz = pbuffer.data(idx_op_geom_010_fd + 11);

    #pragma omp simd aligned(tr_0_0_x_xxy_xx, tr_0_0_x_xxy_xy, tr_0_0_x_xxy_xz, tr_0_0_x_xxy_yy, tr_0_0_x_xxy_yz, tr_0_0_x_xxy_zz, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxy_xx[i] = 2.0 * tr_xxxy_xx[i] * tbe_0 + 2.0 * tr_xxy_xxx[i] * tke_0 - 2.0 * tr_xy_xx[i] - 2.0 * tr_xxy_x[i];

        tr_0_0_x_xxy_xy[i] = 2.0 * tr_xxxy_xy[i] * tbe_0 + 2.0 * tr_xxy_xxy[i] * tke_0 - 2.0 * tr_xy_xy[i] - tr_xxy_y[i];

        tr_0_0_x_xxy_xz[i] = 2.0 * tr_xxxy_xz[i] * tbe_0 + 2.0 * tr_xxy_xxz[i] * tke_0 - 2.0 * tr_xy_xz[i] - tr_xxy_z[i];

        tr_0_0_x_xxy_yy[i] = 2.0 * tr_xxxy_yy[i] * tbe_0 + 2.0 * tr_xxy_xyy[i] * tke_0 - 2.0 * tr_xy_yy[i];

        tr_0_0_x_xxy_yz[i] = 2.0 * tr_xxxy_yz[i] * tbe_0 + 2.0 * tr_xxy_xyz[i] * tke_0 - 2.0 * tr_xy_yz[i];

        tr_0_0_x_xxy_zz[i] = 2.0 * tr_xxxy_zz[i] * tbe_0 + 2.0 * tr_xxy_xzz[i] * tke_0 - 2.0 * tr_xy_zz[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto tr_0_0_x_xxz_xx = pbuffer.data(idx_op_geom_010_fd + 12);

    auto tr_0_0_x_xxz_xy = pbuffer.data(idx_op_geom_010_fd + 13);

    auto tr_0_0_x_xxz_xz = pbuffer.data(idx_op_geom_010_fd + 14);

    auto tr_0_0_x_xxz_yy = pbuffer.data(idx_op_geom_010_fd + 15);

    auto tr_0_0_x_xxz_yz = pbuffer.data(idx_op_geom_010_fd + 16);

    auto tr_0_0_x_xxz_zz = pbuffer.data(idx_op_geom_010_fd + 17);

    #pragma omp simd aligned(tr_0_0_x_xxz_xx, tr_0_0_x_xxz_xy, tr_0_0_x_xxz_xz, tr_0_0_x_xxz_yy, tr_0_0_x_xxz_yz, tr_0_0_x_xxz_zz, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxz_xx[i] = 2.0 * tr_xxxz_xx[i] * tbe_0 + 2.0 * tr_xxz_xxx[i] * tke_0 - 2.0 * tr_xz_xx[i] - 2.0 * tr_xxz_x[i];

        tr_0_0_x_xxz_xy[i] = 2.0 * tr_xxxz_xy[i] * tbe_0 + 2.0 * tr_xxz_xxy[i] * tke_0 - 2.0 * tr_xz_xy[i] - tr_xxz_y[i];

        tr_0_0_x_xxz_xz[i] = 2.0 * tr_xxxz_xz[i] * tbe_0 + 2.0 * tr_xxz_xxz[i] * tke_0 - 2.0 * tr_xz_xz[i] - tr_xxz_z[i];

        tr_0_0_x_xxz_yy[i] = 2.0 * tr_xxxz_yy[i] * tbe_0 + 2.0 * tr_xxz_xyy[i] * tke_0 - 2.0 * tr_xz_yy[i];

        tr_0_0_x_xxz_yz[i] = 2.0 * tr_xxxz_yz[i] * tbe_0 + 2.0 * tr_xxz_xyz[i] * tke_0 - 2.0 * tr_xz_yz[i];

        tr_0_0_x_xxz_zz[i] = 2.0 * tr_xxxz_zz[i] * tbe_0 + 2.0 * tr_xxz_xzz[i] * tke_0 - 2.0 * tr_xz_zz[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto tr_0_0_x_xyy_xx = pbuffer.data(idx_op_geom_010_fd + 18);

    auto tr_0_0_x_xyy_xy = pbuffer.data(idx_op_geom_010_fd + 19);

    auto tr_0_0_x_xyy_xz = pbuffer.data(idx_op_geom_010_fd + 20);

    auto tr_0_0_x_xyy_yy = pbuffer.data(idx_op_geom_010_fd + 21);

    auto tr_0_0_x_xyy_yz = pbuffer.data(idx_op_geom_010_fd + 22);

    auto tr_0_0_x_xyy_zz = pbuffer.data(idx_op_geom_010_fd + 23);

    #pragma omp simd aligned(tr_0_0_x_xyy_xx, tr_0_0_x_xyy_xy, tr_0_0_x_xyy_xz, tr_0_0_x_xyy_yy, tr_0_0_x_xyy_yz, tr_0_0_x_xyy_zz, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyy_xx[i] = 2.0 * tr_xxyy_xx[i] * tbe_0 + 2.0 * tr_xyy_xxx[i] * tke_0 - tr_yy_xx[i] - 2.0 * tr_xyy_x[i];

        tr_0_0_x_xyy_xy[i] = 2.0 * tr_xxyy_xy[i] * tbe_0 + 2.0 * tr_xyy_xxy[i] * tke_0 - tr_yy_xy[i] - tr_xyy_y[i];

        tr_0_0_x_xyy_xz[i] = 2.0 * tr_xxyy_xz[i] * tbe_0 + 2.0 * tr_xyy_xxz[i] * tke_0 - tr_yy_xz[i] - tr_xyy_z[i];

        tr_0_0_x_xyy_yy[i] = 2.0 * tr_xxyy_yy[i] * tbe_0 + 2.0 * tr_xyy_xyy[i] * tke_0 - tr_yy_yy[i];

        tr_0_0_x_xyy_yz[i] = 2.0 * tr_xxyy_yz[i] * tbe_0 + 2.0 * tr_xyy_xyz[i] * tke_0 - tr_yy_yz[i];

        tr_0_0_x_xyy_zz[i] = 2.0 * tr_xxyy_zz[i] * tbe_0 + 2.0 * tr_xyy_xzz[i] * tke_0 - tr_yy_zz[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto tr_0_0_x_xyz_xx = pbuffer.data(idx_op_geom_010_fd + 24);

    auto tr_0_0_x_xyz_xy = pbuffer.data(idx_op_geom_010_fd + 25);

    auto tr_0_0_x_xyz_xz = pbuffer.data(idx_op_geom_010_fd + 26);

    auto tr_0_0_x_xyz_yy = pbuffer.data(idx_op_geom_010_fd + 27);

    auto tr_0_0_x_xyz_yz = pbuffer.data(idx_op_geom_010_fd + 28);

    auto tr_0_0_x_xyz_zz = pbuffer.data(idx_op_geom_010_fd + 29);

    #pragma omp simd aligned(tr_0_0_x_xyz_xx, tr_0_0_x_xyz_xy, tr_0_0_x_xyz_xz, tr_0_0_x_xyz_yy, tr_0_0_x_xyz_yz, tr_0_0_x_xyz_zz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyz_xx[i] = 2.0 * tr_xxyz_xx[i] * tbe_0 + 2.0 * tr_xyz_xxx[i] * tke_0 - tr_yz_xx[i] - 2.0 * tr_xyz_x[i];

        tr_0_0_x_xyz_xy[i] = 2.0 * tr_xxyz_xy[i] * tbe_0 + 2.0 * tr_xyz_xxy[i] * tke_0 - tr_yz_xy[i] - tr_xyz_y[i];

        tr_0_0_x_xyz_xz[i] = 2.0 * tr_xxyz_xz[i] * tbe_0 + 2.0 * tr_xyz_xxz[i] * tke_0 - tr_yz_xz[i] - tr_xyz_z[i];

        tr_0_0_x_xyz_yy[i] = 2.0 * tr_xxyz_yy[i] * tbe_0 + 2.0 * tr_xyz_xyy[i] * tke_0 - tr_yz_yy[i];

        tr_0_0_x_xyz_yz[i] = 2.0 * tr_xxyz_yz[i] * tbe_0 + 2.0 * tr_xyz_xyz[i] * tke_0 - tr_yz_yz[i];

        tr_0_0_x_xyz_zz[i] = 2.0 * tr_xxyz_zz[i] * tbe_0 + 2.0 * tr_xyz_xzz[i] * tke_0 - tr_yz_zz[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto tr_0_0_x_xzz_xx = pbuffer.data(idx_op_geom_010_fd + 30);

    auto tr_0_0_x_xzz_xy = pbuffer.data(idx_op_geom_010_fd + 31);

    auto tr_0_0_x_xzz_xz = pbuffer.data(idx_op_geom_010_fd + 32);

    auto tr_0_0_x_xzz_yy = pbuffer.data(idx_op_geom_010_fd + 33);

    auto tr_0_0_x_xzz_yz = pbuffer.data(idx_op_geom_010_fd + 34);

    auto tr_0_0_x_xzz_zz = pbuffer.data(idx_op_geom_010_fd + 35);

    #pragma omp simd aligned(tr_0_0_x_xzz_xx, tr_0_0_x_xzz_xy, tr_0_0_x_xzz_xz, tr_0_0_x_xzz_yy, tr_0_0_x_xzz_yz, tr_0_0_x_xzz_zz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzz_xx[i] = 2.0 * tr_xxzz_xx[i] * tbe_0 + 2.0 * tr_xzz_xxx[i] * tke_0 - tr_zz_xx[i] - 2.0 * tr_xzz_x[i];

        tr_0_0_x_xzz_xy[i] = 2.0 * tr_xxzz_xy[i] * tbe_0 + 2.0 * tr_xzz_xxy[i] * tke_0 - tr_zz_xy[i] - tr_xzz_y[i];

        tr_0_0_x_xzz_xz[i] = 2.0 * tr_xxzz_xz[i] * tbe_0 + 2.0 * tr_xzz_xxz[i] * tke_0 - tr_zz_xz[i] - tr_xzz_z[i];

        tr_0_0_x_xzz_yy[i] = 2.0 * tr_xxzz_yy[i] * tbe_0 + 2.0 * tr_xzz_xyy[i] * tke_0 - tr_zz_yy[i];

        tr_0_0_x_xzz_yz[i] = 2.0 * tr_xxzz_yz[i] * tbe_0 + 2.0 * tr_xzz_xyz[i] * tke_0 - tr_zz_yz[i];

        tr_0_0_x_xzz_zz[i] = 2.0 * tr_xxzz_zz[i] * tbe_0 + 2.0 * tr_xzz_xzz[i] * tke_0 - tr_zz_zz[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto tr_0_0_x_yyy_xx = pbuffer.data(idx_op_geom_010_fd + 36);

    auto tr_0_0_x_yyy_xy = pbuffer.data(idx_op_geom_010_fd + 37);

    auto tr_0_0_x_yyy_xz = pbuffer.data(idx_op_geom_010_fd + 38);

    auto tr_0_0_x_yyy_yy = pbuffer.data(idx_op_geom_010_fd + 39);

    auto tr_0_0_x_yyy_yz = pbuffer.data(idx_op_geom_010_fd + 40);

    auto tr_0_0_x_yyy_zz = pbuffer.data(idx_op_geom_010_fd + 41);

    #pragma omp simd aligned(tr_0_0_x_yyy_xx, tr_0_0_x_yyy_xy, tr_0_0_x_yyy_xz, tr_0_0_x_yyy_yy, tr_0_0_x_yyy_yz, tr_0_0_x_yyy_zz, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyy_xx[i] = 2.0 * tr_xyyy_xx[i] * tbe_0 + 2.0 * tr_yyy_xxx[i] * tke_0 - 2.0 * tr_yyy_x[i];

        tr_0_0_x_yyy_xy[i] = 2.0 * tr_xyyy_xy[i] * tbe_0 + 2.0 * tr_yyy_xxy[i] * tke_0 - tr_yyy_y[i];

        tr_0_0_x_yyy_xz[i] = 2.0 * tr_xyyy_xz[i] * tbe_0 + 2.0 * tr_yyy_xxz[i] * tke_0 - tr_yyy_z[i];

        tr_0_0_x_yyy_yy[i] = 2.0 * tr_xyyy_yy[i] * tbe_0 + 2.0 * tr_yyy_xyy[i] * tke_0;

        tr_0_0_x_yyy_yz[i] = 2.0 * tr_xyyy_yz[i] * tbe_0 + 2.0 * tr_yyy_xyz[i] * tke_0;

        tr_0_0_x_yyy_zz[i] = 2.0 * tr_xyyy_zz[i] * tbe_0 + 2.0 * tr_yyy_xzz[i] * tke_0;
    }

    // Set up 42-48 components of targeted buffer : FD

    auto tr_0_0_x_yyz_xx = pbuffer.data(idx_op_geom_010_fd + 42);

    auto tr_0_0_x_yyz_xy = pbuffer.data(idx_op_geom_010_fd + 43);

    auto tr_0_0_x_yyz_xz = pbuffer.data(idx_op_geom_010_fd + 44);

    auto tr_0_0_x_yyz_yy = pbuffer.data(idx_op_geom_010_fd + 45);

    auto tr_0_0_x_yyz_yz = pbuffer.data(idx_op_geom_010_fd + 46);

    auto tr_0_0_x_yyz_zz = pbuffer.data(idx_op_geom_010_fd + 47);

    #pragma omp simd aligned(tr_0_0_x_yyz_xx, tr_0_0_x_yyz_xy, tr_0_0_x_yyz_xz, tr_0_0_x_yyz_yy, tr_0_0_x_yyz_yz, tr_0_0_x_yyz_zz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyz_xx[i] = 2.0 * tr_xyyz_xx[i] * tbe_0 + 2.0 * tr_yyz_xxx[i] * tke_0 - 2.0 * tr_yyz_x[i];

        tr_0_0_x_yyz_xy[i] = 2.0 * tr_xyyz_xy[i] * tbe_0 + 2.0 * tr_yyz_xxy[i] * tke_0 - tr_yyz_y[i];

        tr_0_0_x_yyz_xz[i] = 2.0 * tr_xyyz_xz[i] * tbe_0 + 2.0 * tr_yyz_xxz[i] * tke_0 - tr_yyz_z[i];

        tr_0_0_x_yyz_yy[i] = 2.0 * tr_xyyz_yy[i] * tbe_0 + 2.0 * tr_yyz_xyy[i] * tke_0;

        tr_0_0_x_yyz_yz[i] = 2.0 * tr_xyyz_yz[i] * tbe_0 + 2.0 * tr_yyz_xyz[i] * tke_0;

        tr_0_0_x_yyz_zz[i] = 2.0 * tr_xyyz_zz[i] * tbe_0 + 2.0 * tr_yyz_xzz[i] * tke_0;
    }

    // Set up 48-54 components of targeted buffer : FD

    auto tr_0_0_x_yzz_xx = pbuffer.data(idx_op_geom_010_fd + 48);

    auto tr_0_0_x_yzz_xy = pbuffer.data(idx_op_geom_010_fd + 49);

    auto tr_0_0_x_yzz_xz = pbuffer.data(idx_op_geom_010_fd + 50);

    auto tr_0_0_x_yzz_yy = pbuffer.data(idx_op_geom_010_fd + 51);

    auto tr_0_0_x_yzz_yz = pbuffer.data(idx_op_geom_010_fd + 52);

    auto tr_0_0_x_yzz_zz = pbuffer.data(idx_op_geom_010_fd + 53);

    #pragma omp simd aligned(tr_0_0_x_yzz_xx, tr_0_0_x_yzz_xy, tr_0_0_x_yzz_xz, tr_0_0_x_yzz_yy, tr_0_0_x_yzz_yz, tr_0_0_x_yzz_zz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzz_xx[i] = 2.0 * tr_xyzz_xx[i] * tbe_0 + 2.0 * tr_yzz_xxx[i] * tke_0 - 2.0 * tr_yzz_x[i];

        tr_0_0_x_yzz_xy[i] = 2.0 * tr_xyzz_xy[i] * tbe_0 + 2.0 * tr_yzz_xxy[i] * tke_0 - tr_yzz_y[i];

        tr_0_0_x_yzz_xz[i] = 2.0 * tr_xyzz_xz[i] * tbe_0 + 2.0 * tr_yzz_xxz[i] * tke_0 - tr_yzz_z[i];

        tr_0_0_x_yzz_yy[i] = 2.0 * tr_xyzz_yy[i] * tbe_0 + 2.0 * tr_yzz_xyy[i] * tke_0;

        tr_0_0_x_yzz_yz[i] = 2.0 * tr_xyzz_yz[i] * tbe_0 + 2.0 * tr_yzz_xyz[i] * tke_0;

        tr_0_0_x_yzz_zz[i] = 2.0 * tr_xyzz_zz[i] * tbe_0 + 2.0 * tr_yzz_xzz[i] * tke_0;
    }

    // Set up 54-60 components of targeted buffer : FD

    auto tr_0_0_x_zzz_xx = pbuffer.data(idx_op_geom_010_fd + 54);

    auto tr_0_0_x_zzz_xy = pbuffer.data(idx_op_geom_010_fd + 55);

    auto tr_0_0_x_zzz_xz = pbuffer.data(idx_op_geom_010_fd + 56);

    auto tr_0_0_x_zzz_yy = pbuffer.data(idx_op_geom_010_fd + 57);

    auto tr_0_0_x_zzz_yz = pbuffer.data(idx_op_geom_010_fd + 58);

    auto tr_0_0_x_zzz_zz = pbuffer.data(idx_op_geom_010_fd + 59);

    #pragma omp simd aligned(tr_0_0_x_zzz_xx, tr_0_0_x_zzz_xy, tr_0_0_x_zzz_xz, tr_0_0_x_zzz_yy, tr_0_0_x_zzz_yz, tr_0_0_x_zzz_zz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzz_xx[i] = 2.0 * tr_xzzz_xx[i] * tbe_0 + 2.0 * tr_zzz_xxx[i] * tke_0 - 2.0 * tr_zzz_x[i];

        tr_0_0_x_zzz_xy[i] = 2.0 * tr_xzzz_xy[i] * tbe_0 + 2.0 * tr_zzz_xxy[i] * tke_0 - tr_zzz_y[i];

        tr_0_0_x_zzz_xz[i] = 2.0 * tr_xzzz_xz[i] * tbe_0 + 2.0 * tr_zzz_xxz[i] * tke_0 - tr_zzz_z[i];

        tr_0_0_x_zzz_yy[i] = 2.0 * tr_xzzz_yy[i] * tbe_0 + 2.0 * tr_zzz_xyy[i] * tke_0;

        tr_0_0_x_zzz_yz[i] = 2.0 * tr_xzzz_yz[i] * tbe_0 + 2.0 * tr_zzz_xyz[i] * tke_0;

        tr_0_0_x_zzz_zz[i] = 2.0 * tr_xzzz_zz[i] * tbe_0 + 2.0 * tr_zzz_xzz[i] * tke_0;
    }

    // Set up 60-66 components of targeted buffer : FD

    auto tr_0_0_y_xxx_xx = pbuffer.data(idx_op_geom_010_fd + 60);

    auto tr_0_0_y_xxx_xy = pbuffer.data(idx_op_geom_010_fd + 61);

    auto tr_0_0_y_xxx_xz = pbuffer.data(idx_op_geom_010_fd + 62);

    auto tr_0_0_y_xxx_yy = pbuffer.data(idx_op_geom_010_fd + 63);

    auto tr_0_0_y_xxx_yz = pbuffer.data(idx_op_geom_010_fd + 64);

    auto tr_0_0_y_xxx_zz = pbuffer.data(idx_op_geom_010_fd + 65);

    #pragma omp simd aligned(tr_0_0_y_xxx_xx, tr_0_0_y_xxx_xy, tr_0_0_y_xxx_xz, tr_0_0_y_xxx_yy, tr_0_0_y_xxx_yz, tr_0_0_y_xxx_zz, tr_xxx_x, tr_xxx_xxy, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxx_xx[i] = 2.0 * tr_xxxy_xx[i] * tbe_0 + 2.0 * tr_xxx_xxy[i] * tke_0;

        tr_0_0_y_xxx_xy[i] = 2.0 * tr_xxxy_xy[i] * tbe_0 + 2.0 * tr_xxx_xyy[i] * tke_0 - tr_xxx_x[i];

        tr_0_0_y_xxx_xz[i] = 2.0 * tr_xxxy_xz[i] * tbe_0 + 2.0 * tr_xxx_xyz[i] * tke_0;

        tr_0_0_y_xxx_yy[i] = 2.0 * tr_xxxy_yy[i] * tbe_0 + 2.0 * tr_xxx_yyy[i] * tke_0 - 2.0 * tr_xxx_y[i];

        tr_0_0_y_xxx_yz[i] = 2.0 * tr_xxxy_yz[i] * tbe_0 + 2.0 * tr_xxx_yyz[i] * tke_0 - tr_xxx_z[i];

        tr_0_0_y_xxx_zz[i] = 2.0 * tr_xxxy_zz[i] * tbe_0 + 2.0 * tr_xxx_yzz[i] * tke_0;
    }

    // Set up 66-72 components of targeted buffer : FD

    auto tr_0_0_y_xxy_xx = pbuffer.data(idx_op_geom_010_fd + 66);

    auto tr_0_0_y_xxy_xy = pbuffer.data(idx_op_geom_010_fd + 67);

    auto tr_0_0_y_xxy_xz = pbuffer.data(idx_op_geom_010_fd + 68);

    auto tr_0_0_y_xxy_yy = pbuffer.data(idx_op_geom_010_fd + 69);

    auto tr_0_0_y_xxy_yz = pbuffer.data(idx_op_geom_010_fd + 70);

    auto tr_0_0_y_xxy_zz = pbuffer.data(idx_op_geom_010_fd + 71);

    #pragma omp simd aligned(tr_0_0_y_xxy_xx, tr_0_0_y_xxy_xy, tr_0_0_y_xxy_xz, tr_0_0_y_xxy_yy, tr_0_0_y_xxy_yz, tr_0_0_y_xxy_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxy_xx[i] = 2.0 * tr_xxyy_xx[i] * tbe_0 + 2.0 * tr_xxy_xxy[i] * tke_0 - tr_xx_xx[i];

        tr_0_0_y_xxy_xy[i] = 2.0 * tr_xxyy_xy[i] * tbe_0 + 2.0 * tr_xxy_xyy[i] * tke_0 - tr_xx_xy[i] - tr_xxy_x[i];

        tr_0_0_y_xxy_xz[i] = 2.0 * tr_xxyy_xz[i] * tbe_0 + 2.0 * tr_xxy_xyz[i] * tke_0 - tr_xx_xz[i];

        tr_0_0_y_xxy_yy[i] = 2.0 * tr_xxyy_yy[i] * tbe_0 + 2.0 * tr_xxy_yyy[i] * tke_0 - tr_xx_yy[i] - 2.0 * tr_xxy_y[i];

        tr_0_0_y_xxy_yz[i] = 2.0 * tr_xxyy_yz[i] * tbe_0 + 2.0 * tr_xxy_yyz[i] * tke_0 - tr_xx_yz[i] - tr_xxy_z[i];

        tr_0_0_y_xxy_zz[i] = 2.0 * tr_xxyy_zz[i] * tbe_0 + 2.0 * tr_xxy_yzz[i] * tke_0 - tr_xx_zz[i];
    }

    // Set up 72-78 components of targeted buffer : FD

    auto tr_0_0_y_xxz_xx = pbuffer.data(idx_op_geom_010_fd + 72);

    auto tr_0_0_y_xxz_xy = pbuffer.data(idx_op_geom_010_fd + 73);

    auto tr_0_0_y_xxz_xz = pbuffer.data(idx_op_geom_010_fd + 74);

    auto tr_0_0_y_xxz_yy = pbuffer.data(idx_op_geom_010_fd + 75);

    auto tr_0_0_y_xxz_yz = pbuffer.data(idx_op_geom_010_fd + 76);

    auto tr_0_0_y_xxz_zz = pbuffer.data(idx_op_geom_010_fd + 77);

    #pragma omp simd aligned(tr_0_0_y_xxz_xx, tr_0_0_y_xxz_xy, tr_0_0_y_xxz_xz, tr_0_0_y_xxz_yy, tr_0_0_y_xxz_yz, tr_0_0_y_xxz_zz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxz_xx[i] = 2.0 * tr_xxyz_xx[i] * tbe_0 + 2.0 * tr_xxz_xxy[i] * tke_0;

        tr_0_0_y_xxz_xy[i] = 2.0 * tr_xxyz_xy[i] * tbe_0 + 2.0 * tr_xxz_xyy[i] * tke_0 - tr_xxz_x[i];

        tr_0_0_y_xxz_xz[i] = 2.0 * tr_xxyz_xz[i] * tbe_0 + 2.0 * tr_xxz_xyz[i] * tke_0;

        tr_0_0_y_xxz_yy[i] = 2.0 * tr_xxyz_yy[i] * tbe_0 + 2.0 * tr_xxz_yyy[i] * tke_0 - 2.0 * tr_xxz_y[i];

        tr_0_0_y_xxz_yz[i] = 2.0 * tr_xxyz_yz[i] * tbe_0 + 2.0 * tr_xxz_yyz[i] * tke_0 - tr_xxz_z[i];

        tr_0_0_y_xxz_zz[i] = 2.0 * tr_xxyz_zz[i] * tbe_0 + 2.0 * tr_xxz_yzz[i] * tke_0;
    }

    // Set up 78-84 components of targeted buffer : FD

    auto tr_0_0_y_xyy_xx = pbuffer.data(idx_op_geom_010_fd + 78);

    auto tr_0_0_y_xyy_xy = pbuffer.data(idx_op_geom_010_fd + 79);

    auto tr_0_0_y_xyy_xz = pbuffer.data(idx_op_geom_010_fd + 80);

    auto tr_0_0_y_xyy_yy = pbuffer.data(idx_op_geom_010_fd + 81);

    auto tr_0_0_y_xyy_yz = pbuffer.data(idx_op_geom_010_fd + 82);

    auto tr_0_0_y_xyy_zz = pbuffer.data(idx_op_geom_010_fd + 83);

    #pragma omp simd aligned(tr_0_0_y_xyy_xx, tr_0_0_y_xyy_xy, tr_0_0_y_xyy_xz, tr_0_0_y_xyy_yy, tr_0_0_y_xyy_yz, tr_0_0_y_xyy_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyy_xx[i] = 2.0 * tr_xyyy_xx[i] * tbe_0 + 2.0 * tr_xyy_xxy[i] * tke_0 - 2.0 * tr_xy_xx[i];

        tr_0_0_y_xyy_xy[i] = 2.0 * tr_xyyy_xy[i] * tbe_0 + 2.0 * tr_xyy_xyy[i] * tke_0 - 2.0 * tr_xy_xy[i] - tr_xyy_x[i];

        tr_0_0_y_xyy_xz[i] = 2.0 * tr_xyyy_xz[i] * tbe_0 + 2.0 * tr_xyy_xyz[i] * tke_0 - 2.0 * tr_xy_xz[i];

        tr_0_0_y_xyy_yy[i] = 2.0 * tr_xyyy_yy[i] * tbe_0 + 2.0 * tr_xyy_yyy[i] * tke_0 - 2.0 * tr_xy_yy[i] - 2.0 * tr_xyy_y[i];

        tr_0_0_y_xyy_yz[i] = 2.0 * tr_xyyy_yz[i] * tbe_0 + 2.0 * tr_xyy_yyz[i] * tke_0 - 2.0 * tr_xy_yz[i] - tr_xyy_z[i];

        tr_0_0_y_xyy_zz[i] = 2.0 * tr_xyyy_zz[i] * tbe_0 + 2.0 * tr_xyy_yzz[i] * tke_0 - 2.0 * tr_xy_zz[i];
    }

    // Set up 84-90 components of targeted buffer : FD

    auto tr_0_0_y_xyz_xx = pbuffer.data(idx_op_geom_010_fd + 84);

    auto tr_0_0_y_xyz_xy = pbuffer.data(idx_op_geom_010_fd + 85);

    auto tr_0_0_y_xyz_xz = pbuffer.data(idx_op_geom_010_fd + 86);

    auto tr_0_0_y_xyz_yy = pbuffer.data(idx_op_geom_010_fd + 87);

    auto tr_0_0_y_xyz_yz = pbuffer.data(idx_op_geom_010_fd + 88);

    auto tr_0_0_y_xyz_zz = pbuffer.data(idx_op_geom_010_fd + 89);

    #pragma omp simd aligned(tr_0_0_y_xyz_xx, tr_0_0_y_xyz_xy, tr_0_0_y_xyz_xz, tr_0_0_y_xyz_yy, tr_0_0_y_xyz_yz, tr_0_0_y_xyz_zz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyz_xx[i] = 2.0 * tr_xyyz_xx[i] * tbe_0 + 2.0 * tr_xyz_xxy[i] * tke_0 - tr_xz_xx[i];

        tr_0_0_y_xyz_xy[i] = 2.0 * tr_xyyz_xy[i] * tbe_0 + 2.0 * tr_xyz_xyy[i] * tke_0 - tr_xz_xy[i] - tr_xyz_x[i];

        tr_0_0_y_xyz_xz[i] = 2.0 * tr_xyyz_xz[i] * tbe_0 + 2.0 * tr_xyz_xyz[i] * tke_0 - tr_xz_xz[i];

        tr_0_0_y_xyz_yy[i] = 2.0 * tr_xyyz_yy[i] * tbe_0 + 2.0 * tr_xyz_yyy[i] * tke_0 - tr_xz_yy[i] - 2.0 * tr_xyz_y[i];

        tr_0_0_y_xyz_yz[i] = 2.0 * tr_xyyz_yz[i] * tbe_0 + 2.0 * tr_xyz_yyz[i] * tke_0 - tr_xz_yz[i] - tr_xyz_z[i];

        tr_0_0_y_xyz_zz[i] = 2.0 * tr_xyyz_zz[i] * tbe_0 + 2.0 * tr_xyz_yzz[i] * tke_0 - tr_xz_zz[i];
    }

    // Set up 90-96 components of targeted buffer : FD

    auto tr_0_0_y_xzz_xx = pbuffer.data(idx_op_geom_010_fd + 90);

    auto tr_0_0_y_xzz_xy = pbuffer.data(idx_op_geom_010_fd + 91);

    auto tr_0_0_y_xzz_xz = pbuffer.data(idx_op_geom_010_fd + 92);

    auto tr_0_0_y_xzz_yy = pbuffer.data(idx_op_geom_010_fd + 93);

    auto tr_0_0_y_xzz_yz = pbuffer.data(idx_op_geom_010_fd + 94);

    auto tr_0_0_y_xzz_zz = pbuffer.data(idx_op_geom_010_fd + 95);

    #pragma omp simd aligned(tr_0_0_y_xzz_xx, tr_0_0_y_xzz_xy, tr_0_0_y_xzz_xz, tr_0_0_y_xzz_yy, tr_0_0_y_xzz_yz, tr_0_0_y_xzz_zz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzz_xx[i] = 2.0 * tr_xyzz_xx[i] * tbe_0 + 2.0 * tr_xzz_xxy[i] * tke_0;

        tr_0_0_y_xzz_xy[i] = 2.0 * tr_xyzz_xy[i] * tbe_0 + 2.0 * tr_xzz_xyy[i] * tke_0 - tr_xzz_x[i];

        tr_0_0_y_xzz_xz[i] = 2.0 * tr_xyzz_xz[i] * tbe_0 + 2.0 * tr_xzz_xyz[i] * tke_0;

        tr_0_0_y_xzz_yy[i] = 2.0 * tr_xyzz_yy[i] * tbe_0 + 2.0 * tr_xzz_yyy[i] * tke_0 - 2.0 * tr_xzz_y[i];

        tr_0_0_y_xzz_yz[i] = 2.0 * tr_xyzz_yz[i] * tbe_0 + 2.0 * tr_xzz_yyz[i] * tke_0 - tr_xzz_z[i];

        tr_0_0_y_xzz_zz[i] = 2.0 * tr_xyzz_zz[i] * tbe_0 + 2.0 * tr_xzz_yzz[i] * tke_0;
    }

    // Set up 96-102 components of targeted buffer : FD

    auto tr_0_0_y_yyy_xx = pbuffer.data(idx_op_geom_010_fd + 96);

    auto tr_0_0_y_yyy_xy = pbuffer.data(idx_op_geom_010_fd + 97);

    auto tr_0_0_y_yyy_xz = pbuffer.data(idx_op_geom_010_fd + 98);

    auto tr_0_0_y_yyy_yy = pbuffer.data(idx_op_geom_010_fd + 99);

    auto tr_0_0_y_yyy_yz = pbuffer.data(idx_op_geom_010_fd + 100);

    auto tr_0_0_y_yyy_zz = pbuffer.data(idx_op_geom_010_fd + 101);

    #pragma omp simd aligned(tr_0_0_y_yyy_xx, tr_0_0_y_yyy_xy, tr_0_0_y_yyy_xz, tr_0_0_y_yyy_yy, tr_0_0_y_yyy_yz, tr_0_0_y_yyy_zz, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyy_x, tr_yyy_xxy, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyy_xx[i] = 2.0 * tr_yyyy_xx[i] * tbe_0 + 2.0 * tr_yyy_xxy[i] * tke_0 - 3.0 * tr_yy_xx[i];

        tr_0_0_y_yyy_xy[i] = 2.0 * tr_yyyy_xy[i] * tbe_0 + 2.0 * tr_yyy_xyy[i] * tke_0 - 3.0 * tr_yy_xy[i] - tr_yyy_x[i];

        tr_0_0_y_yyy_xz[i] = 2.0 * tr_yyyy_xz[i] * tbe_0 + 2.0 * tr_yyy_xyz[i] * tke_0 - 3.0 * tr_yy_xz[i];

        tr_0_0_y_yyy_yy[i] = 2.0 * tr_yyyy_yy[i] * tbe_0 + 2.0 * tr_yyy_yyy[i] * tke_0 - 3.0 * tr_yy_yy[i] - 2.0 * tr_yyy_y[i];

        tr_0_0_y_yyy_yz[i] = 2.0 * tr_yyyy_yz[i] * tbe_0 + 2.0 * tr_yyy_yyz[i] * tke_0 - 3.0 * tr_yy_yz[i] - tr_yyy_z[i];

        tr_0_0_y_yyy_zz[i] = 2.0 * tr_yyyy_zz[i] * tbe_0 + 2.0 * tr_yyy_yzz[i] * tke_0 - 3.0 * tr_yy_zz[i];
    }

    // Set up 102-108 components of targeted buffer : FD

    auto tr_0_0_y_yyz_xx = pbuffer.data(idx_op_geom_010_fd + 102);

    auto tr_0_0_y_yyz_xy = pbuffer.data(idx_op_geom_010_fd + 103);

    auto tr_0_0_y_yyz_xz = pbuffer.data(idx_op_geom_010_fd + 104);

    auto tr_0_0_y_yyz_yy = pbuffer.data(idx_op_geom_010_fd + 105);

    auto tr_0_0_y_yyz_yz = pbuffer.data(idx_op_geom_010_fd + 106);

    auto tr_0_0_y_yyz_zz = pbuffer.data(idx_op_geom_010_fd + 107);

    #pragma omp simd aligned(tr_0_0_y_yyz_xx, tr_0_0_y_yyz_xy, tr_0_0_y_yyz_xz, tr_0_0_y_yyz_yy, tr_0_0_y_yyz_yz, tr_0_0_y_yyz_zz, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyz_xx[i] = 2.0 * tr_yyyz_xx[i] * tbe_0 + 2.0 * tr_yyz_xxy[i] * tke_0 - 2.0 * tr_yz_xx[i];

        tr_0_0_y_yyz_xy[i] = 2.0 * tr_yyyz_xy[i] * tbe_0 + 2.0 * tr_yyz_xyy[i] * tke_0 - 2.0 * tr_yz_xy[i] - tr_yyz_x[i];

        tr_0_0_y_yyz_xz[i] = 2.0 * tr_yyyz_xz[i] * tbe_0 + 2.0 * tr_yyz_xyz[i] * tke_0 - 2.0 * tr_yz_xz[i];

        tr_0_0_y_yyz_yy[i] = 2.0 * tr_yyyz_yy[i] * tbe_0 + 2.0 * tr_yyz_yyy[i] * tke_0 - 2.0 * tr_yz_yy[i] - 2.0 * tr_yyz_y[i];

        tr_0_0_y_yyz_yz[i] = 2.0 * tr_yyyz_yz[i] * tbe_0 + 2.0 * tr_yyz_yyz[i] * tke_0 - 2.0 * tr_yz_yz[i] - tr_yyz_z[i];

        tr_0_0_y_yyz_zz[i] = 2.0 * tr_yyyz_zz[i] * tbe_0 + 2.0 * tr_yyz_yzz[i] * tke_0 - 2.0 * tr_yz_zz[i];
    }

    // Set up 108-114 components of targeted buffer : FD

    auto tr_0_0_y_yzz_xx = pbuffer.data(idx_op_geom_010_fd + 108);

    auto tr_0_0_y_yzz_xy = pbuffer.data(idx_op_geom_010_fd + 109);

    auto tr_0_0_y_yzz_xz = pbuffer.data(idx_op_geom_010_fd + 110);

    auto tr_0_0_y_yzz_yy = pbuffer.data(idx_op_geom_010_fd + 111);

    auto tr_0_0_y_yzz_yz = pbuffer.data(idx_op_geom_010_fd + 112);

    auto tr_0_0_y_yzz_zz = pbuffer.data(idx_op_geom_010_fd + 113);

    #pragma omp simd aligned(tr_0_0_y_yzz_xx, tr_0_0_y_yzz_xy, tr_0_0_y_yzz_xz, tr_0_0_y_yzz_yy, tr_0_0_y_yzz_yz, tr_0_0_y_yzz_zz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzz_xx[i] = 2.0 * tr_yyzz_xx[i] * tbe_0 + 2.0 * tr_yzz_xxy[i] * tke_0 - tr_zz_xx[i];

        tr_0_0_y_yzz_xy[i] = 2.0 * tr_yyzz_xy[i] * tbe_0 + 2.0 * tr_yzz_xyy[i] * tke_0 - tr_zz_xy[i] - tr_yzz_x[i];

        tr_0_0_y_yzz_xz[i] = 2.0 * tr_yyzz_xz[i] * tbe_0 + 2.0 * tr_yzz_xyz[i] * tke_0 - tr_zz_xz[i];

        tr_0_0_y_yzz_yy[i] = 2.0 * tr_yyzz_yy[i] * tbe_0 + 2.0 * tr_yzz_yyy[i] * tke_0 - tr_zz_yy[i] - 2.0 * tr_yzz_y[i];

        tr_0_0_y_yzz_yz[i] = 2.0 * tr_yyzz_yz[i] * tbe_0 + 2.0 * tr_yzz_yyz[i] * tke_0 - tr_zz_yz[i] - tr_yzz_z[i];

        tr_0_0_y_yzz_zz[i] = 2.0 * tr_yyzz_zz[i] * tbe_0 + 2.0 * tr_yzz_yzz[i] * tke_0 - tr_zz_zz[i];
    }

    // Set up 114-120 components of targeted buffer : FD

    auto tr_0_0_y_zzz_xx = pbuffer.data(idx_op_geom_010_fd + 114);

    auto tr_0_0_y_zzz_xy = pbuffer.data(idx_op_geom_010_fd + 115);

    auto tr_0_0_y_zzz_xz = pbuffer.data(idx_op_geom_010_fd + 116);

    auto tr_0_0_y_zzz_yy = pbuffer.data(idx_op_geom_010_fd + 117);

    auto tr_0_0_y_zzz_yz = pbuffer.data(idx_op_geom_010_fd + 118);

    auto tr_0_0_y_zzz_zz = pbuffer.data(idx_op_geom_010_fd + 119);

    #pragma omp simd aligned(tr_0_0_y_zzz_xx, tr_0_0_y_zzz_xy, tr_0_0_y_zzz_xz, tr_0_0_y_zzz_yy, tr_0_0_y_zzz_yz, tr_0_0_y_zzz_zz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, tr_zzz_x, tr_zzz_xxy, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzz_xx[i] = 2.0 * tr_yzzz_xx[i] * tbe_0 + 2.0 * tr_zzz_xxy[i] * tke_0;

        tr_0_0_y_zzz_xy[i] = 2.0 * tr_yzzz_xy[i] * tbe_0 + 2.0 * tr_zzz_xyy[i] * tke_0 - tr_zzz_x[i];

        tr_0_0_y_zzz_xz[i] = 2.0 * tr_yzzz_xz[i] * tbe_0 + 2.0 * tr_zzz_xyz[i] * tke_0;

        tr_0_0_y_zzz_yy[i] = 2.0 * tr_yzzz_yy[i] * tbe_0 + 2.0 * tr_zzz_yyy[i] * tke_0 - 2.0 * tr_zzz_y[i];

        tr_0_0_y_zzz_yz[i] = 2.0 * tr_yzzz_yz[i] * tbe_0 + 2.0 * tr_zzz_yyz[i] * tke_0 - tr_zzz_z[i];

        tr_0_0_y_zzz_zz[i] = 2.0 * tr_yzzz_zz[i] * tbe_0 + 2.0 * tr_zzz_yzz[i] * tke_0;
    }

    // Set up 120-126 components of targeted buffer : FD

    auto tr_0_0_z_xxx_xx = pbuffer.data(idx_op_geom_010_fd + 120);

    auto tr_0_0_z_xxx_xy = pbuffer.data(idx_op_geom_010_fd + 121);

    auto tr_0_0_z_xxx_xz = pbuffer.data(idx_op_geom_010_fd + 122);

    auto tr_0_0_z_xxx_yy = pbuffer.data(idx_op_geom_010_fd + 123);

    auto tr_0_0_z_xxx_yz = pbuffer.data(idx_op_geom_010_fd + 124);

    auto tr_0_0_z_xxx_zz = pbuffer.data(idx_op_geom_010_fd + 125);

    #pragma omp simd aligned(tr_0_0_z_xxx_xx, tr_0_0_z_xxx_xy, tr_0_0_z_xxx_xz, tr_0_0_z_xxx_yy, tr_0_0_z_xxx_yz, tr_0_0_z_xxx_zz, tr_xxx_x, tr_xxx_xxz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxx_xx[i] = 2.0 * tr_xxxz_xx[i] * tbe_0 + 2.0 * tr_xxx_xxz[i] * tke_0;

        tr_0_0_z_xxx_xy[i] = 2.0 * tr_xxxz_xy[i] * tbe_0 + 2.0 * tr_xxx_xyz[i] * tke_0;

        tr_0_0_z_xxx_xz[i] = 2.0 * tr_xxxz_xz[i] * tbe_0 + 2.0 * tr_xxx_xzz[i] * tke_0 - tr_xxx_x[i];

        tr_0_0_z_xxx_yy[i] = 2.0 * tr_xxxz_yy[i] * tbe_0 + 2.0 * tr_xxx_yyz[i] * tke_0;

        tr_0_0_z_xxx_yz[i] = 2.0 * tr_xxxz_yz[i] * tbe_0 + 2.0 * tr_xxx_yzz[i] * tke_0 - tr_xxx_y[i];

        tr_0_0_z_xxx_zz[i] = 2.0 * tr_xxxz_zz[i] * tbe_0 + 2.0 * tr_xxx_zzz[i] * tke_0 - 2.0 * tr_xxx_z[i];
    }

    // Set up 126-132 components of targeted buffer : FD

    auto tr_0_0_z_xxy_xx = pbuffer.data(idx_op_geom_010_fd + 126);

    auto tr_0_0_z_xxy_xy = pbuffer.data(idx_op_geom_010_fd + 127);

    auto tr_0_0_z_xxy_xz = pbuffer.data(idx_op_geom_010_fd + 128);

    auto tr_0_0_z_xxy_yy = pbuffer.data(idx_op_geom_010_fd + 129);

    auto tr_0_0_z_xxy_yz = pbuffer.data(idx_op_geom_010_fd + 130);

    auto tr_0_0_z_xxy_zz = pbuffer.data(idx_op_geom_010_fd + 131);

    #pragma omp simd aligned(tr_0_0_z_xxy_xx, tr_0_0_z_xxy_xy, tr_0_0_z_xxy_xz, tr_0_0_z_xxy_yy, tr_0_0_z_xxy_yz, tr_0_0_z_xxy_zz, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxy_xx[i] = 2.0 * tr_xxyz_xx[i] * tbe_0 + 2.0 * tr_xxy_xxz[i] * tke_0;

        tr_0_0_z_xxy_xy[i] = 2.0 * tr_xxyz_xy[i] * tbe_0 + 2.0 * tr_xxy_xyz[i] * tke_0;

        tr_0_0_z_xxy_xz[i] = 2.0 * tr_xxyz_xz[i] * tbe_0 + 2.0 * tr_xxy_xzz[i] * tke_0 - tr_xxy_x[i];

        tr_0_0_z_xxy_yy[i] = 2.0 * tr_xxyz_yy[i] * tbe_0 + 2.0 * tr_xxy_yyz[i] * tke_0;

        tr_0_0_z_xxy_yz[i] = 2.0 * tr_xxyz_yz[i] * tbe_0 + 2.0 * tr_xxy_yzz[i] * tke_0 - tr_xxy_y[i];

        tr_0_0_z_xxy_zz[i] = 2.0 * tr_xxyz_zz[i] * tbe_0 + 2.0 * tr_xxy_zzz[i] * tke_0 - 2.0 * tr_xxy_z[i];
    }

    // Set up 132-138 components of targeted buffer : FD

    auto tr_0_0_z_xxz_xx = pbuffer.data(idx_op_geom_010_fd + 132);

    auto tr_0_0_z_xxz_xy = pbuffer.data(idx_op_geom_010_fd + 133);

    auto tr_0_0_z_xxz_xz = pbuffer.data(idx_op_geom_010_fd + 134);

    auto tr_0_0_z_xxz_yy = pbuffer.data(idx_op_geom_010_fd + 135);

    auto tr_0_0_z_xxz_yz = pbuffer.data(idx_op_geom_010_fd + 136);

    auto tr_0_0_z_xxz_zz = pbuffer.data(idx_op_geom_010_fd + 137);

    #pragma omp simd aligned(tr_0_0_z_xxz_xx, tr_0_0_z_xxz_xy, tr_0_0_z_xxz_xz, tr_0_0_z_xxz_yy, tr_0_0_z_xxz_yz, tr_0_0_z_xxz_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxz_xx[i] = 2.0 * tr_xxzz_xx[i] * tbe_0 + 2.0 * tr_xxz_xxz[i] * tke_0 - tr_xx_xx[i];

        tr_0_0_z_xxz_xy[i] = 2.0 * tr_xxzz_xy[i] * tbe_0 + 2.0 * tr_xxz_xyz[i] * tke_0 - tr_xx_xy[i];

        tr_0_0_z_xxz_xz[i] = 2.0 * tr_xxzz_xz[i] * tbe_0 + 2.0 * tr_xxz_xzz[i] * tke_0 - tr_xx_xz[i] - tr_xxz_x[i];

        tr_0_0_z_xxz_yy[i] = 2.0 * tr_xxzz_yy[i] * tbe_0 + 2.0 * tr_xxz_yyz[i] * tke_0 - tr_xx_yy[i];

        tr_0_0_z_xxz_yz[i] = 2.0 * tr_xxzz_yz[i] * tbe_0 + 2.0 * tr_xxz_yzz[i] * tke_0 - tr_xx_yz[i] - tr_xxz_y[i];

        tr_0_0_z_xxz_zz[i] = 2.0 * tr_xxzz_zz[i] * tbe_0 + 2.0 * tr_xxz_zzz[i] * tke_0 - tr_xx_zz[i] - 2.0 * tr_xxz_z[i];
    }

    // Set up 138-144 components of targeted buffer : FD

    auto tr_0_0_z_xyy_xx = pbuffer.data(idx_op_geom_010_fd + 138);

    auto tr_0_0_z_xyy_xy = pbuffer.data(idx_op_geom_010_fd + 139);

    auto tr_0_0_z_xyy_xz = pbuffer.data(idx_op_geom_010_fd + 140);

    auto tr_0_0_z_xyy_yy = pbuffer.data(idx_op_geom_010_fd + 141);

    auto tr_0_0_z_xyy_yz = pbuffer.data(idx_op_geom_010_fd + 142);

    auto tr_0_0_z_xyy_zz = pbuffer.data(idx_op_geom_010_fd + 143);

    #pragma omp simd aligned(tr_0_0_z_xyy_xx, tr_0_0_z_xyy_xy, tr_0_0_z_xyy_xz, tr_0_0_z_xyy_yy, tr_0_0_z_xyy_yz, tr_0_0_z_xyy_zz, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyy_xx[i] = 2.0 * tr_xyyz_xx[i] * tbe_0 + 2.0 * tr_xyy_xxz[i] * tke_0;

        tr_0_0_z_xyy_xy[i] = 2.0 * tr_xyyz_xy[i] * tbe_0 + 2.0 * tr_xyy_xyz[i] * tke_0;

        tr_0_0_z_xyy_xz[i] = 2.0 * tr_xyyz_xz[i] * tbe_0 + 2.0 * tr_xyy_xzz[i] * tke_0 - tr_xyy_x[i];

        tr_0_0_z_xyy_yy[i] = 2.0 * tr_xyyz_yy[i] * tbe_0 + 2.0 * tr_xyy_yyz[i] * tke_0;

        tr_0_0_z_xyy_yz[i] = 2.0 * tr_xyyz_yz[i] * tbe_0 + 2.0 * tr_xyy_yzz[i] * tke_0 - tr_xyy_y[i];

        tr_0_0_z_xyy_zz[i] = 2.0 * tr_xyyz_zz[i] * tbe_0 + 2.0 * tr_xyy_zzz[i] * tke_0 - 2.0 * tr_xyy_z[i];
    }

    // Set up 144-150 components of targeted buffer : FD

    auto tr_0_0_z_xyz_xx = pbuffer.data(idx_op_geom_010_fd + 144);

    auto tr_0_0_z_xyz_xy = pbuffer.data(idx_op_geom_010_fd + 145);

    auto tr_0_0_z_xyz_xz = pbuffer.data(idx_op_geom_010_fd + 146);

    auto tr_0_0_z_xyz_yy = pbuffer.data(idx_op_geom_010_fd + 147);

    auto tr_0_0_z_xyz_yz = pbuffer.data(idx_op_geom_010_fd + 148);

    auto tr_0_0_z_xyz_zz = pbuffer.data(idx_op_geom_010_fd + 149);

    #pragma omp simd aligned(tr_0_0_z_xyz_xx, tr_0_0_z_xyz_xy, tr_0_0_z_xyz_xz, tr_0_0_z_xyz_yy, tr_0_0_z_xyz_yz, tr_0_0_z_xyz_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyz_xx[i] = 2.0 * tr_xyzz_xx[i] * tbe_0 + 2.0 * tr_xyz_xxz[i] * tke_0 - tr_xy_xx[i];

        tr_0_0_z_xyz_xy[i] = 2.0 * tr_xyzz_xy[i] * tbe_0 + 2.0 * tr_xyz_xyz[i] * tke_0 - tr_xy_xy[i];

        tr_0_0_z_xyz_xz[i] = 2.0 * tr_xyzz_xz[i] * tbe_0 + 2.0 * tr_xyz_xzz[i] * tke_0 - tr_xy_xz[i] - tr_xyz_x[i];

        tr_0_0_z_xyz_yy[i] = 2.0 * tr_xyzz_yy[i] * tbe_0 + 2.0 * tr_xyz_yyz[i] * tke_0 - tr_xy_yy[i];

        tr_0_0_z_xyz_yz[i] = 2.0 * tr_xyzz_yz[i] * tbe_0 + 2.0 * tr_xyz_yzz[i] * tke_0 - tr_xy_yz[i] - tr_xyz_y[i];

        tr_0_0_z_xyz_zz[i] = 2.0 * tr_xyzz_zz[i] * tbe_0 + 2.0 * tr_xyz_zzz[i] * tke_0 - tr_xy_zz[i] - 2.0 * tr_xyz_z[i];
    }

    // Set up 150-156 components of targeted buffer : FD

    auto tr_0_0_z_xzz_xx = pbuffer.data(idx_op_geom_010_fd + 150);

    auto tr_0_0_z_xzz_xy = pbuffer.data(idx_op_geom_010_fd + 151);

    auto tr_0_0_z_xzz_xz = pbuffer.data(idx_op_geom_010_fd + 152);

    auto tr_0_0_z_xzz_yy = pbuffer.data(idx_op_geom_010_fd + 153);

    auto tr_0_0_z_xzz_yz = pbuffer.data(idx_op_geom_010_fd + 154);

    auto tr_0_0_z_xzz_zz = pbuffer.data(idx_op_geom_010_fd + 155);

    #pragma omp simd aligned(tr_0_0_z_xzz_xx, tr_0_0_z_xzz_xy, tr_0_0_z_xzz_xz, tr_0_0_z_xzz_yy, tr_0_0_z_xzz_yz, tr_0_0_z_xzz_zz, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzz_xx[i] = 2.0 * tr_xzzz_xx[i] * tbe_0 + 2.0 * tr_xzz_xxz[i] * tke_0 - 2.0 * tr_xz_xx[i];

        tr_0_0_z_xzz_xy[i] = 2.0 * tr_xzzz_xy[i] * tbe_0 + 2.0 * tr_xzz_xyz[i] * tke_0 - 2.0 * tr_xz_xy[i];

        tr_0_0_z_xzz_xz[i] = 2.0 * tr_xzzz_xz[i] * tbe_0 + 2.0 * tr_xzz_xzz[i] * tke_0 - 2.0 * tr_xz_xz[i] - tr_xzz_x[i];

        tr_0_0_z_xzz_yy[i] = 2.0 * tr_xzzz_yy[i] * tbe_0 + 2.0 * tr_xzz_yyz[i] * tke_0 - 2.0 * tr_xz_yy[i];

        tr_0_0_z_xzz_yz[i] = 2.0 * tr_xzzz_yz[i] * tbe_0 + 2.0 * tr_xzz_yzz[i] * tke_0 - 2.0 * tr_xz_yz[i] - tr_xzz_y[i];

        tr_0_0_z_xzz_zz[i] = 2.0 * tr_xzzz_zz[i] * tbe_0 + 2.0 * tr_xzz_zzz[i] * tke_0 - 2.0 * tr_xz_zz[i] - 2.0 * tr_xzz_z[i];
    }

    // Set up 156-162 components of targeted buffer : FD

    auto tr_0_0_z_yyy_xx = pbuffer.data(idx_op_geom_010_fd + 156);

    auto tr_0_0_z_yyy_xy = pbuffer.data(idx_op_geom_010_fd + 157);

    auto tr_0_0_z_yyy_xz = pbuffer.data(idx_op_geom_010_fd + 158);

    auto tr_0_0_z_yyy_yy = pbuffer.data(idx_op_geom_010_fd + 159);

    auto tr_0_0_z_yyy_yz = pbuffer.data(idx_op_geom_010_fd + 160);

    auto tr_0_0_z_yyy_zz = pbuffer.data(idx_op_geom_010_fd + 161);

    #pragma omp simd aligned(tr_0_0_z_yyy_xx, tr_0_0_z_yyy_xy, tr_0_0_z_yyy_xz, tr_0_0_z_yyy_yy, tr_0_0_z_yyy_yz, tr_0_0_z_yyy_zz, tr_yyy_x, tr_yyy_xxz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyy_xx[i] = 2.0 * tr_yyyz_xx[i] * tbe_0 + 2.0 * tr_yyy_xxz[i] * tke_0;

        tr_0_0_z_yyy_xy[i] = 2.0 * tr_yyyz_xy[i] * tbe_0 + 2.0 * tr_yyy_xyz[i] * tke_0;

        tr_0_0_z_yyy_xz[i] = 2.0 * tr_yyyz_xz[i] * tbe_0 + 2.0 * tr_yyy_xzz[i] * tke_0 - tr_yyy_x[i];

        tr_0_0_z_yyy_yy[i] = 2.0 * tr_yyyz_yy[i] * tbe_0 + 2.0 * tr_yyy_yyz[i] * tke_0;

        tr_0_0_z_yyy_yz[i] = 2.0 * tr_yyyz_yz[i] * tbe_0 + 2.0 * tr_yyy_yzz[i] * tke_0 - tr_yyy_y[i];

        tr_0_0_z_yyy_zz[i] = 2.0 * tr_yyyz_zz[i] * tbe_0 + 2.0 * tr_yyy_zzz[i] * tke_0 - 2.0 * tr_yyy_z[i];
    }

    // Set up 162-168 components of targeted buffer : FD

    auto tr_0_0_z_yyz_xx = pbuffer.data(idx_op_geom_010_fd + 162);

    auto tr_0_0_z_yyz_xy = pbuffer.data(idx_op_geom_010_fd + 163);

    auto tr_0_0_z_yyz_xz = pbuffer.data(idx_op_geom_010_fd + 164);

    auto tr_0_0_z_yyz_yy = pbuffer.data(idx_op_geom_010_fd + 165);

    auto tr_0_0_z_yyz_yz = pbuffer.data(idx_op_geom_010_fd + 166);

    auto tr_0_0_z_yyz_zz = pbuffer.data(idx_op_geom_010_fd + 167);

    #pragma omp simd aligned(tr_0_0_z_yyz_xx, tr_0_0_z_yyz_xy, tr_0_0_z_yyz_xz, tr_0_0_z_yyz_yy, tr_0_0_z_yyz_yz, tr_0_0_z_yyz_zz, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyz_xx[i] = 2.0 * tr_yyzz_xx[i] * tbe_0 + 2.0 * tr_yyz_xxz[i] * tke_0 - tr_yy_xx[i];

        tr_0_0_z_yyz_xy[i] = 2.0 * tr_yyzz_xy[i] * tbe_0 + 2.0 * tr_yyz_xyz[i] * tke_0 - tr_yy_xy[i];

        tr_0_0_z_yyz_xz[i] = 2.0 * tr_yyzz_xz[i] * tbe_0 + 2.0 * tr_yyz_xzz[i] * tke_0 - tr_yy_xz[i] - tr_yyz_x[i];

        tr_0_0_z_yyz_yy[i] = 2.0 * tr_yyzz_yy[i] * tbe_0 + 2.0 * tr_yyz_yyz[i] * tke_0 - tr_yy_yy[i];

        tr_0_0_z_yyz_yz[i] = 2.0 * tr_yyzz_yz[i] * tbe_0 + 2.0 * tr_yyz_yzz[i] * tke_0 - tr_yy_yz[i] - tr_yyz_y[i];

        tr_0_0_z_yyz_zz[i] = 2.0 * tr_yyzz_zz[i] * tbe_0 + 2.0 * tr_yyz_zzz[i] * tke_0 - tr_yy_zz[i] - 2.0 * tr_yyz_z[i];
    }

    // Set up 168-174 components of targeted buffer : FD

    auto tr_0_0_z_yzz_xx = pbuffer.data(idx_op_geom_010_fd + 168);

    auto tr_0_0_z_yzz_xy = pbuffer.data(idx_op_geom_010_fd + 169);

    auto tr_0_0_z_yzz_xz = pbuffer.data(idx_op_geom_010_fd + 170);

    auto tr_0_0_z_yzz_yy = pbuffer.data(idx_op_geom_010_fd + 171);

    auto tr_0_0_z_yzz_yz = pbuffer.data(idx_op_geom_010_fd + 172);

    auto tr_0_0_z_yzz_zz = pbuffer.data(idx_op_geom_010_fd + 173);

    #pragma omp simd aligned(tr_0_0_z_yzz_xx, tr_0_0_z_yzz_xy, tr_0_0_z_yzz_xz, tr_0_0_z_yzz_yy, tr_0_0_z_yzz_yz, tr_0_0_z_yzz_zz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzz_xx[i] = 2.0 * tr_yzzz_xx[i] * tbe_0 + 2.0 * tr_yzz_xxz[i] * tke_0 - 2.0 * tr_yz_xx[i];

        tr_0_0_z_yzz_xy[i] = 2.0 * tr_yzzz_xy[i] * tbe_0 + 2.0 * tr_yzz_xyz[i] * tke_0 - 2.0 * tr_yz_xy[i];

        tr_0_0_z_yzz_xz[i] = 2.0 * tr_yzzz_xz[i] * tbe_0 + 2.0 * tr_yzz_xzz[i] * tke_0 - 2.0 * tr_yz_xz[i] - tr_yzz_x[i];

        tr_0_0_z_yzz_yy[i] = 2.0 * tr_yzzz_yy[i] * tbe_0 + 2.0 * tr_yzz_yyz[i] * tke_0 - 2.0 * tr_yz_yy[i];

        tr_0_0_z_yzz_yz[i] = 2.0 * tr_yzzz_yz[i] * tbe_0 + 2.0 * tr_yzz_yzz[i] * tke_0 - 2.0 * tr_yz_yz[i] - tr_yzz_y[i];

        tr_0_0_z_yzz_zz[i] = 2.0 * tr_yzzz_zz[i] * tbe_0 + 2.0 * tr_yzz_zzz[i] * tke_0 - 2.0 * tr_yz_zz[i] - 2.0 * tr_yzz_z[i];
    }

    // Set up 174-180 components of targeted buffer : FD

    auto tr_0_0_z_zzz_xx = pbuffer.data(idx_op_geom_010_fd + 174);

    auto tr_0_0_z_zzz_xy = pbuffer.data(idx_op_geom_010_fd + 175);

    auto tr_0_0_z_zzz_xz = pbuffer.data(idx_op_geom_010_fd + 176);

    auto tr_0_0_z_zzz_yy = pbuffer.data(idx_op_geom_010_fd + 177);

    auto tr_0_0_z_zzz_yz = pbuffer.data(idx_op_geom_010_fd + 178);

    auto tr_0_0_z_zzz_zz = pbuffer.data(idx_op_geom_010_fd + 179);

    #pragma omp simd aligned(tr_0_0_z_zzz_xx, tr_0_0_z_zzz_xy, tr_0_0_z_zzz_xz, tr_0_0_z_zzz_yy, tr_0_0_z_zzz_yz, tr_0_0_z_zzz_zz, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xxz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzz_xx[i] = 2.0 * tr_zzzz_xx[i] * tbe_0 + 2.0 * tr_zzz_xxz[i] * tke_0 - 3.0 * tr_zz_xx[i];

        tr_0_0_z_zzz_xy[i] = 2.0 * tr_zzzz_xy[i] * tbe_0 + 2.0 * tr_zzz_xyz[i] * tke_0 - 3.0 * tr_zz_xy[i];

        tr_0_0_z_zzz_xz[i] = 2.0 * tr_zzzz_xz[i] * tbe_0 + 2.0 * tr_zzz_xzz[i] * tke_0 - 3.0 * tr_zz_xz[i] - tr_zzz_x[i];

        tr_0_0_z_zzz_yy[i] = 2.0 * tr_zzzz_yy[i] * tbe_0 + 2.0 * tr_zzz_yyz[i] * tke_0 - 3.0 * tr_zz_yy[i];

        tr_0_0_z_zzz_yz[i] = 2.0 * tr_zzzz_yz[i] * tbe_0 + 2.0 * tr_zzz_yzz[i] * tke_0 - 3.0 * tr_zz_yz[i] - tr_zzz_y[i];

        tr_0_0_z_zzz_zz[i] = 2.0 * tr_zzzz_zz[i] * tbe_0 + 2.0 * tr_zzz_zzz[i] * tke_0 - 3.0 * tr_zz_zz[i] - 2.0 * tr_zzz_z[i];
    }

}

} // t2cgeom namespace

