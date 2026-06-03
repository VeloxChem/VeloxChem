#include "GeometricalDerivatives020ForDD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_dd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_dd,
                         const int idx_op_sd,
                         const int idx_op_pp,
                         const int idx_op_pf,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_dg,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_gd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SD

    auto tr_0_xx = pbuffer.data(idx_op_sd);

    auto tr_0_xy = pbuffer.data(idx_op_sd + 1);

    auto tr_0_xz = pbuffer.data(idx_op_sd + 2);

    auto tr_0_yy = pbuffer.data(idx_op_sd + 3);

    auto tr_0_yz = pbuffer.data(idx_op_sd + 4);

    auto tr_0_zz = pbuffer.data(idx_op_sd + 5);

    // Set up components of auxiliary buffer : PP

    auto tr_x_x = pbuffer.data(idx_op_pp);

    auto tr_x_y = pbuffer.data(idx_op_pp + 1);

    auto tr_x_z = pbuffer.data(idx_op_pp + 2);

    auto tr_y_x = pbuffer.data(idx_op_pp + 3);

    auto tr_y_y = pbuffer.data(idx_op_pp + 4);

    auto tr_y_z = pbuffer.data(idx_op_pp + 5);

    auto tr_z_x = pbuffer.data(idx_op_pp + 6);

    auto tr_z_y = pbuffer.data(idx_op_pp + 7);

    auto tr_z_z = pbuffer.data(idx_op_pp + 8);

    // Set up components of auxiliary buffer : PF

    auto tr_x_xxx = pbuffer.data(idx_op_pf);

    auto tr_x_xxy = pbuffer.data(idx_op_pf + 1);

    auto tr_x_xxz = pbuffer.data(idx_op_pf + 2);

    auto tr_x_xyy = pbuffer.data(idx_op_pf + 3);

    auto tr_x_xyz = pbuffer.data(idx_op_pf + 4);

    auto tr_x_xzz = pbuffer.data(idx_op_pf + 5);

    auto tr_x_yyy = pbuffer.data(idx_op_pf + 6);

    auto tr_x_yyz = pbuffer.data(idx_op_pf + 7);

    auto tr_x_yzz = pbuffer.data(idx_op_pf + 8);

    auto tr_x_zzz = pbuffer.data(idx_op_pf + 9);

    auto tr_y_xxx = pbuffer.data(idx_op_pf + 10);

    auto tr_y_xxy = pbuffer.data(idx_op_pf + 11);

    auto tr_y_xxz = pbuffer.data(idx_op_pf + 12);

    auto tr_y_xyy = pbuffer.data(idx_op_pf + 13);

    auto tr_y_xyz = pbuffer.data(idx_op_pf + 14);

    auto tr_y_xzz = pbuffer.data(idx_op_pf + 15);

    auto tr_y_yyy = pbuffer.data(idx_op_pf + 16);

    auto tr_y_yyz = pbuffer.data(idx_op_pf + 17);

    auto tr_y_yzz = pbuffer.data(idx_op_pf + 18);

    auto tr_y_zzz = pbuffer.data(idx_op_pf + 19);

    auto tr_z_xxx = pbuffer.data(idx_op_pf + 20);

    auto tr_z_xxy = pbuffer.data(idx_op_pf + 21);

    auto tr_z_xxz = pbuffer.data(idx_op_pf + 22);

    auto tr_z_xyy = pbuffer.data(idx_op_pf + 23);

    auto tr_z_xyz = pbuffer.data(idx_op_pf + 24);

    auto tr_z_xzz = pbuffer.data(idx_op_pf + 25);

    auto tr_z_yyy = pbuffer.data(idx_op_pf + 26);

    auto tr_z_yyz = pbuffer.data(idx_op_pf + 27);

    auto tr_z_yzz = pbuffer.data(idx_op_pf + 28);

    auto tr_z_zzz = pbuffer.data(idx_op_pf + 29);

    // Set up components of auxiliary buffer : DS

    auto tr_xx_0 = pbuffer.data(idx_op_ds);

    auto tr_xy_0 = pbuffer.data(idx_op_ds + 1);

    auto tr_xz_0 = pbuffer.data(idx_op_ds + 2);

    auto tr_yy_0 = pbuffer.data(idx_op_ds + 3);

    auto tr_yz_0 = pbuffer.data(idx_op_ds + 4);

    auto tr_zz_0 = pbuffer.data(idx_op_ds + 5);

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

    // Set up components of auxiliary buffer : DG

    auto tr_xx_xxxx = pbuffer.data(idx_op_dg);

    auto tr_xx_xxxy = pbuffer.data(idx_op_dg + 1);

    auto tr_xx_xxxz = pbuffer.data(idx_op_dg + 2);

    auto tr_xx_xxyy = pbuffer.data(idx_op_dg + 3);

    auto tr_xx_xxyz = pbuffer.data(idx_op_dg + 4);

    auto tr_xx_xxzz = pbuffer.data(idx_op_dg + 5);

    auto tr_xx_xyyy = pbuffer.data(idx_op_dg + 6);

    auto tr_xx_xyyz = pbuffer.data(idx_op_dg + 7);

    auto tr_xx_xyzz = pbuffer.data(idx_op_dg + 8);

    auto tr_xx_xzzz = pbuffer.data(idx_op_dg + 9);

    auto tr_xx_yyyy = pbuffer.data(idx_op_dg + 10);

    auto tr_xx_yyyz = pbuffer.data(idx_op_dg + 11);

    auto tr_xx_yyzz = pbuffer.data(idx_op_dg + 12);

    auto tr_xx_yzzz = pbuffer.data(idx_op_dg + 13);

    auto tr_xx_zzzz = pbuffer.data(idx_op_dg + 14);

    auto tr_xy_xxxx = pbuffer.data(idx_op_dg + 15);

    auto tr_xy_xxxy = pbuffer.data(idx_op_dg + 16);

    auto tr_xy_xxxz = pbuffer.data(idx_op_dg + 17);

    auto tr_xy_xxyy = pbuffer.data(idx_op_dg + 18);

    auto tr_xy_xxyz = pbuffer.data(idx_op_dg + 19);

    auto tr_xy_xxzz = pbuffer.data(idx_op_dg + 20);

    auto tr_xy_xyyy = pbuffer.data(idx_op_dg + 21);

    auto tr_xy_xyyz = pbuffer.data(idx_op_dg + 22);

    auto tr_xy_xyzz = pbuffer.data(idx_op_dg + 23);

    auto tr_xy_xzzz = pbuffer.data(idx_op_dg + 24);

    auto tr_xy_yyyy = pbuffer.data(idx_op_dg + 25);

    auto tr_xy_yyyz = pbuffer.data(idx_op_dg + 26);

    auto tr_xy_yyzz = pbuffer.data(idx_op_dg + 27);

    auto tr_xy_yzzz = pbuffer.data(idx_op_dg + 28);

    auto tr_xy_zzzz = pbuffer.data(idx_op_dg + 29);

    auto tr_xz_xxxx = pbuffer.data(idx_op_dg + 30);

    auto tr_xz_xxxy = pbuffer.data(idx_op_dg + 31);

    auto tr_xz_xxxz = pbuffer.data(idx_op_dg + 32);

    auto tr_xz_xxyy = pbuffer.data(idx_op_dg + 33);

    auto tr_xz_xxyz = pbuffer.data(idx_op_dg + 34);

    auto tr_xz_xxzz = pbuffer.data(idx_op_dg + 35);

    auto tr_xz_xyyy = pbuffer.data(idx_op_dg + 36);

    auto tr_xz_xyyz = pbuffer.data(idx_op_dg + 37);

    auto tr_xz_xyzz = pbuffer.data(idx_op_dg + 38);

    auto tr_xz_xzzz = pbuffer.data(idx_op_dg + 39);

    auto tr_xz_yyyy = pbuffer.data(idx_op_dg + 40);

    auto tr_xz_yyyz = pbuffer.data(idx_op_dg + 41);

    auto tr_xz_yyzz = pbuffer.data(idx_op_dg + 42);

    auto tr_xz_yzzz = pbuffer.data(idx_op_dg + 43);

    auto tr_xz_zzzz = pbuffer.data(idx_op_dg + 44);

    auto tr_yy_xxxx = pbuffer.data(idx_op_dg + 45);

    auto tr_yy_xxxy = pbuffer.data(idx_op_dg + 46);

    auto tr_yy_xxxz = pbuffer.data(idx_op_dg + 47);

    auto tr_yy_xxyy = pbuffer.data(idx_op_dg + 48);

    auto tr_yy_xxyz = pbuffer.data(idx_op_dg + 49);

    auto tr_yy_xxzz = pbuffer.data(idx_op_dg + 50);

    auto tr_yy_xyyy = pbuffer.data(idx_op_dg + 51);

    auto tr_yy_xyyz = pbuffer.data(idx_op_dg + 52);

    auto tr_yy_xyzz = pbuffer.data(idx_op_dg + 53);

    auto tr_yy_xzzz = pbuffer.data(idx_op_dg + 54);

    auto tr_yy_yyyy = pbuffer.data(idx_op_dg + 55);

    auto tr_yy_yyyz = pbuffer.data(idx_op_dg + 56);

    auto tr_yy_yyzz = pbuffer.data(idx_op_dg + 57);

    auto tr_yy_yzzz = pbuffer.data(idx_op_dg + 58);

    auto tr_yy_zzzz = pbuffer.data(idx_op_dg + 59);

    auto tr_yz_xxxx = pbuffer.data(idx_op_dg + 60);

    auto tr_yz_xxxy = pbuffer.data(idx_op_dg + 61);

    auto tr_yz_xxxz = pbuffer.data(idx_op_dg + 62);

    auto tr_yz_xxyy = pbuffer.data(idx_op_dg + 63);

    auto tr_yz_xxyz = pbuffer.data(idx_op_dg + 64);

    auto tr_yz_xxzz = pbuffer.data(idx_op_dg + 65);

    auto tr_yz_xyyy = pbuffer.data(idx_op_dg + 66);

    auto tr_yz_xyyz = pbuffer.data(idx_op_dg + 67);

    auto tr_yz_xyzz = pbuffer.data(idx_op_dg + 68);

    auto tr_yz_xzzz = pbuffer.data(idx_op_dg + 69);

    auto tr_yz_yyyy = pbuffer.data(idx_op_dg + 70);

    auto tr_yz_yyyz = pbuffer.data(idx_op_dg + 71);

    auto tr_yz_yyzz = pbuffer.data(idx_op_dg + 72);

    auto tr_yz_yzzz = pbuffer.data(idx_op_dg + 73);

    auto tr_yz_zzzz = pbuffer.data(idx_op_dg + 74);

    auto tr_zz_xxxx = pbuffer.data(idx_op_dg + 75);

    auto tr_zz_xxxy = pbuffer.data(idx_op_dg + 76);

    auto tr_zz_xxxz = pbuffer.data(idx_op_dg + 77);

    auto tr_zz_xxyy = pbuffer.data(idx_op_dg + 78);

    auto tr_zz_xxyz = pbuffer.data(idx_op_dg + 79);

    auto tr_zz_xxzz = pbuffer.data(idx_op_dg + 80);

    auto tr_zz_xyyy = pbuffer.data(idx_op_dg + 81);

    auto tr_zz_xyyz = pbuffer.data(idx_op_dg + 82);

    auto tr_zz_xyzz = pbuffer.data(idx_op_dg + 83);

    auto tr_zz_xzzz = pbuffer.data(idx_op_dg + 84);

    auto tr_zz_yyyy = pbuffer.data(idx_op_dg + 85);

    auto tr_zz_yyyz = pbuffer.data(idx_op_dg + 86);

    auto tr_zz_yyzz = pbuffer.data(idx_op_dg + 87);

    auto tr_zz_yzzz = pbuffer.data(idx_op_dg + 88);

    auto tr_zz_zzzz = pbuffer.data(idx_op_dg + 89);

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

    // Set up 0-6 components of targeted buffer : DD

    auto tr_0_0_xx_xx_xx = pbuffer.data(idx_op_geom_020_dd);

    auto tr_0_0_xx_xx_xy = pbuffer.data(idx_op_geom_020_dd + 1);

    auto tr_0_0_xx_xx_xz = pbuffer.data(idx_op_geom_020_dd + 2);

    auto tr_0_0_xx_xx_yy = pbuffer.data(idx_op_geom_020_dd + 3);

    auto tr_0_0_xx_xx_yz = pbuffer.data(idx_op_geom_020_dd + 4);

    auto tr_0_0_xx_xx_zz = pbuffer.data(idx_op_geom_020_dd + 5);

    #pragma omp simd aligned(tr_0_0_xx_xx_xx, tr_0_0_xx_xx_xy, tr_0_0_xx_xx_xz, tr_0_0_xx_xx_yy, tr_0_0_xx_xx_yz, tr_0_0_xx_xx_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_x, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_z, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xx_xx[i] = 2.0 * tr_0_xx[i] + 8.0 * tr_x_x[i] - 8.0 * tr_x_xxx[i] * tke_0 + 2.0 * tr_xx_0[i] - 10.0 * tr_xx_xx[i] * tbe_0 - 10.0 * tr_xx_xx[i] * tke_0 + 4.0 * tr_xx_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxx_x[i] * tbe_0 + 8.0 * tr_xxx_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xy[i] = 2.0 * tr_0_xy[i] + 4.0 * tr_x_y[i] - 8.0 * tr_x_xxy[i] * tke_0 - 10.0 * tr_xx_xy[i] * tbe_0 - 6.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxx_y[i] * tbe_0 + 8.0 * tr_xxx_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xz[i] = 2.0 * tr_0_xz[i] + 4.0 * tr_x_z[i] - 8.0 * tr_x_xxz[i] * tke_0 - 10.0 * tr_xx_xz[i] * tbe_0 - 6.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxx_z[i] * tbe_0 + 8.0 * tr_xxx_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_yy[i] = 2.0 * tr_0_yy[i] - 8.0 * tr_x_xyy[i] * tke_0 - 10.0 * tr_xx_yy[i] * tbe_0 - 2.0 * tr_xx_yy[i] * tke_0 + 4.0 * tr_xx_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_yz[i] = 2.0 * tr_0_yz[i] - 8.0 * tr_x_xyz[i] * tke_0 - 10.0 * tr_xx_yz[i] * tbe_0 - 2.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_zz[i] = 2.0 * tr_0_zz[i] - 8.0 * tr_x_xzz[i] * tke_0 - 10.0 * tr_xx_zz[i] * tbe_0 - 2.0 * tr_xx_zz[i] * tke_0 + 4.0 * tr_xx_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 6-12 components of targeted buffer : DD

    auto tr_0_0_xx_xy_xx = pbuffer.data(idx_op_geom_020_dd + 6);

    auto tr_0_0_xx_xy_xy = pbuffer.data(idx_op_geom_020_dd + 7);

    auto tr_0_0_xx_xy_xz = pbuffer.data(idx_op_geom_020_dd + 8);

    auto tr_0_0_xx_xy_yy = pbuffer.data(idx_op_geom_020_dd + 9);

    auto tr_0_0_xx_xy_yz = pbuffer.data(idx_op_geom_020_dd + 10);

    auto tr_0_0_xx_xy_zz = pbuffer.data(idx_op_geom_020_dd + 11);

    #pragma omp simd aligned(tr_0_0_xx_xy_xx, tr_0_0_xx_xy_xy, tr_0_0_xx_xy_xz, tr_0_0_xx_xy_yy, tr_0_0_xx_xy_yz, tr_0_0_xx_xy_zz, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_z, tr_xy_0, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_x, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xy_xx[i] = 4.0 * tr_y_x[i] - 4.0 * tr_y_xxx[i] * tke_0 + 2.0 * tr_xy_0[i] - 6.0 * tr_xy_xx[i] * tbe_0 - 10.0 * tr_xy_xx[i] * tke_0 + 4.0 * tr_xy_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxy_x[i] * tbe_0 + 8.0 * tr_xxy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xy[i] = 2.0 * tr_y_y[i] - 4.0 * tr_y_xxy[i] * tke_0 - 6.0 * tr_xy_xy[i] * tbe_0 - 6.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxy_y[i] * tbe_0 + 8.0 * tr_xxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xz[i] = 2.0 * tr_y_z[i] - 4.0 * tr_y_xxz[i] * tke_0 - 6.0 * tr_xy_xz[i] * tbe_0 - 6.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 8.0 * tr_xxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_yy[i] = -4.0 * tr_y_xyy[i] * tke_0 - 6.0 * tr_xy_yy[i] * tbe_0 - 2.0 * tr_xy_yy[i] * tke_0 + 4.0 * tr_xy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_yz[i] = -4.0 * tr_y_xyz[i] * tke_0 - 6.0 * tr_xy_yz[i] * tbe_0 - 2.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_zz[i] = -4.0 * tr_y_xzz[i] * tke_0 - 6.0 * tr_xy_zz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tke_0 + 4.0 * tr_xy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 12-18 components of targeted buffer : DD

    auto tr_0_0_xx_xz_xx = pbuffer.data(idx_op_geom_020_dd + 12);

    auto tr_0_0_xx_xz_xy = pbuffer.data(idx_op_geom_020_dd + 13);

    auto tr_0_0_xx_xz_xz = pbuffer.data(idx_op_geom_020_dd + 14);

    auto tr_0_0_xx_xz_yy = pbuffer.data(idx_op_geom_020_dd + 15);

    auto tr_0_0_xx_xz_yz = pbuffer.data(idx_op_geom_020_dd + 16);

    auto tr_0_0_xx_xz_zz = pbuffer.data(idx_op_geom_020_dd + 17);

    #pragma omp simd aligned(tr_0_0_xx_xz_xx, tr_0_0_xx_xz_xy, tr_0_0_xx_xz_xz, tr_0_0_xx_xz_yy, tr_0_0_xx_xz_yz, tr_0_0_xx_xz_zz, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_z, tr_xz_0, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_x, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xz_xx[i] = 4.0 * tr_z_x[i] - 4.0 * tr_z_xxx[i] * tke_0 + 2.0 * tr_xz_0[i] - 6.0 * tr_xz_xx[i] * tbe_0 - 10.0 * tr_xz_xx[i] * tke_0 + 4.0 * tr_xz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxz_x[i] * tbe_0 + 8.0 * tr_xxz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xy[i] = 2.0 * tr_z_y[i] - 4.0 * tr_z_xxy[i] * tke_0 - 6.0 * tr_xz_xy[i] * tbe_0 - 6.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 8.0 * tr_xxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xz[i] = 2.0 * tr_z_z[i] - 4.0 * tr_z_xxz[i] * tke_0 - 6.0 * tr_xz_xz[i] * tbe_0 - 6.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_z[i] * tbe_0 + 8.0 * tr_xxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_yy[i] = -4.0 * tr_z_xyy[i] * tke_0 - 6.0 * tr_xz_yy[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tke_0 + 4.0 * tr_xz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_yz[i] = -4.0 * tr_z_xyz[i] * tke_0 - 6.0 * tr_xz_yz[i] * tbe_0 - 2.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_zz[i] = -4.0 * tr_z_xzz[i] * tke_0 - 6.0 * tr_xz_zz[i] * tbe_0 - 2.0 * tr_xz_zz[i] * tke_0 + 4.0 * tr_xz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 18-24 components of targeted buffer : DD

    auto tr_0_0_xx_yy_xx = pbuffer.data(idx_op_geom_020_dd + 18);

    auto tr_0_0_xx_yy_xy = pbuffer.data(idx_op_geom_020_dd + 19);

    auto tr_0_0_xx_yy_xz = pbuffer.data(idx_op_geom_020_dd + 20);

    auto tr_0_0_xx_yy_yy = pbuffer.data(idx_op_geom_020_dd + 21);

    auto tr_0_0_xx_yy_yz = pbuffer.data(idx_op_geom_020_dd + 22);

    auto tr_0_0_xx_yy_zz = pbuffer.data(idx_op_geom_020_dd + 23);

    #pragma omp simd aligned(tr_0_0_xx_yy_xx, tr_0_0_xx_yy_xy, tr_0_0_xx_yy_xz, tr_0_0_xx_yy_yy, tr_0_0_xx_yy_yz, tr_0_0_xx_yy_zz, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_z, tr_yy_0, tr_yy_xx, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yy_xx[i] = 2.0 * tr_yy_0[i] - 2.0 * tr_yy_xx[i] * tbe_0 - 10.0 * tr_yy_xx[i] * tke_0 + 4.0 * tr_yy_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xyy_x[i] * tbe_0 + 8.0 * tr_xyy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xy[i] = -2.0 * tr_yy_xy[i] * tbe_0 - 6.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyy_y[i] * tbe_0 + 8.0 * tr_xyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xz[i] = -2.0 * tr_yy_xz[i] * tbe_0 - 6.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_z[i] * tbe_0 + 8.0 * tr_xyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_yy[i] = -2.0 * tr_yy_yy[i] * tbe_0 - 2.0 * tr_yy_yy[i] * tke_0 + 4.0 * tr_yy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_yz[i] = -2.0 * tr_yy_yz[i] * tbe_0 - 2.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_zz[i] = -2.0 * tr_yy_zz[i] * tbe_0 - 2.0 * tr_yy_zz[i] * tke_0 + 4.0 * tr_yy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 24-30 components of targeted buffer : DD

    auto tr_0_0_xx_yz_xx = pbuffer.data(idx_op_geom_020_dd + 24);

    auto tr_0_0_xx_yz_xy = pbuffer.data(idx_op_geom_020_dd + 25);

    auto tr_0_0_xx_yz_xz = pbuffer.data(idx_op_geom_020_dd + 26);

    auto tr_0_0_xx_yz_yy = pbuffer.data(idx_op_geom_020_dd + 27);

    auto tr_0_0_xx_yz_yz = pbuffer.data(idx_op_geom_020_dd + 28);

    auto tr_0_0_xx_yz_zz = pbuffer.data(idx_op_geom_020_dd + 29);

    #pragma omp simd aligned(tr_0_0_xx_yz_xx, tr_0_0_xx_yz_xy, tr_0_0_xx_yz_xz, tr_0_0_xx_yz_yy, tr_0_0_xx_yz_yz, tr_0_0_xx_yz_zz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_yz_0, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yz_xx[i] = 2.0 * tr_yz_0[i] - 2.0 * tr_yz_xx[i] * tbe_0 - 10.0 * tr_yz_xx[i] * tke_0 + 4.0 * tr_yz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xyz_x[i] * tbe_0 + 8.0 * tr_xyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xy[i] = -2.0 * tr_yz_xy[i] * tbe_0 - 6.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 8.0 * tr_xyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xz[i] = -2.0 * tr_yz_xz[i] * tbe_0 - 6.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 8.0 * tr_xyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_yy[i] = -2.0 * tr_yz_yy[i] * tbe_0 - 2.0 * tr_yz_yy[i] * tke_0 + 4.0 * tr_yz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_yz[i] = -2.0 * tr_yz_yz[i] * tbe_0 - 2.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_zz[i] = -2.0 * tr_yz_zz[i] * tbe_0 - 2.0 * tr_yz_zz[i] * tke_0 + 4.0 * tr_yz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-36 components of targeted buffer : DD

    auto tr_0_0_xx_zz_xx = pbuffer.data(idx_op_geom_020_dd + 30);

    auto tr_0_0_xx_zz_xy = pbuffer.data(idx_op_geom_020_dd + 31);

    auto tr_0_0_xx_zz_xz = pbuffer.data(idx_op_geom_020_dd + 32);

    auto tr_0_0_xx_zz_yy = pbuffer.data(idx_op_geom_020_dd + 33);

    auto tr_0_0_xx_zz_yz = pbuffer.data(idx_op_geom_020_dd + 34);

    auto tr_0_0_xx_zz_zz = pbuffer.data(idx_op_geom_020_dd + 35);

    #pragma omp simd aligned(tr_0_0_xx_zz_xx, tr_0_0_xx_zz_xy, tr_0_0_xx_zz_xz, tr_0_0_xx_zz_yy, tr_0_0_xx_zz_yz, tr_0_0_xx_zz_zz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_z, tr_zz_0, tr_zz_xx, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_zz_xx[i] = 2.0 * tr_zz_0[i] - 2.0 * tr_zz_xx[i] * tbe_0 - 10.0 * tr_zz_xx[i] * tke_0 + 4.0 * tr_zz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xzz_x[i] * tbe_0 + 8.0 * tr_xzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xy[i] = -2.0 * tr_zz_xy[i] * tbe_0 - 6.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xzz_y[i] * tbe_0 + 8.0 * tr_xzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xz[i] = -2.0 * tr_zz_xz[i] * tbe_0 - 6.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_z[i] * tbe_0 + 8.0 * tr_xzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_yy[i] = -2.0 * tr_zz_yy[i] * tbe_0 - 2.0 * tr_zz_yy[i] * tke_0 + 4.0 * tr_zz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_yz[i] = -2.0 * tr_zz_yz[i] * tbe_0 - 2.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_zz[i] = -2.0 * tr_zz_zz[i] * tbe_0 - 2.0 * tr_zz_zz[i] * tke_0 + 4.0 * tr_zz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 36-42 components of targeted buffer : DD

    auto tr_0_0_xy_xx_xx = pbuffer.data(idx_op_geom_020_dd + 36);

    auto tr_0_0_xy_xx_xy = pbuffer.data(idx_op_geom_020_dd + 37);

    auto tr_0_0_xy_xx_xz = pbuffer.data(idx_op_geom_020_dd + 38);

    auto tr_0_0_xy_xx_yy = pbuffer.data(idx_op_geom_020_dd + 39);

    auto tr_0_0_xy_xx_yz = pbuffer.data(idx_op_geom_020_dd + 40);

    auto tr_0_0_xy_xx_zz = pbuffer.data(idx_op_geom_020_dd + 41);

    #pragma omp simd aligned(tr_0_0_xy_xx_xx, tr_0_0_xy_xx_xy, tr_0_0_xy_xx_xz, tr_0_0_xy_xx_yy, tr_0_0_xy_xx_yz, tr_0_0_xy_xx_zz, tr_x_x, tr_x_xxy, tr_x_xyy, tr_x_xyz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xxxy, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xxx_x, tr_xxx_xxy, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xx_xx[i] = -4.0 * tr_x_xxy[i] * tke_0 - 4.0 * tr_xy_xx[i] * tbe_0 - 4.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xy[i] = 2.0 * tr_x_x[i] - 4.0 * tr_x_xyy[i] * tke_0 - 4.0 * tr_xy_xy[i] * tbe_0 + tr_xx_0[i] - 2.0 * tr_xx_yy[i] * tke_0 - 2.0 * tr_xx_xx[i] * tke_0 + 4.0 * tr_xx_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxx_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xz[i] = -4.0 * tr_x_xyz[i] * tke_0 - 4.0 * tr_xy_xz[i] * tbe_0 - 2.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_yy[i] = 4.0 * tr_x_y[i] - 4.0 * tr_x_yyy[i] * tke_0 - 4.0 * tr_xy_yy[i] * tbe_0 - 4.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxx_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_yz[i] = 2.0 * tr_x_z[i] - 4.0 * tr_x_yyz[i] * tke_0 - 4.0 * tr_xy_yz[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxx_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_zz[i] = -4.0 * tr_x_yzz[i] * tke_0 - 4.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xx_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 42-48 components of targeted buffer : DD

    auto tr_0_0_xy_xy_xx = pbuffer.data(idx_op_geom_020_dd + 42);

    auto tr_0_0_xy_xy_xy = pbuffer.data(idx_op_geom_020_dd + 43);

    auto tr_0_0_xy_xy_xz = pbuffer.data(idx_op_geom_020_dd + 44);

    auto tr_0_0_xy_xy_yy = pbuffer.data(idx_op_geom_020_dd + 45);

    auto tr_0_0_xy_xy_yz = pbuffer.data(idx_op_geom_020_dd + 46);

    auto tr_0_0_xy_xy_zz = pbuffer.data(idx_op_geom_020_dd + 47);

    #pragma omp simd aligned(tr_0_0_xy_xy_xx, tr_0_0_xy_xy_xy, tr_0_0_xy_xy_xz, tr_0_0_xy_xy_yy, tr_0_0_xy_xy_yz, tr_0_0_xy_xy_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_x, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_z, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xy_0, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_z, tr_y_x, tr_y_xxy, tr_y_xyy, tr_y_xyz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xy_xx[i] = tr_0_xx[i] - 2.0 * tr_y_xxy[i] * tke_0 - 2.0 * tr_yy_xx[i] * tbe_0 + 2.0 * tr_x_x[i] - 2.0 * tr_x_xxx[i] * tke_0 - 4.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyy_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xy[i] = tr_0_xy[i] + tr_y_x[i] - 2.0 * tr_y_xyy[i] * tke_0 - 2.0 * tr_yy_xy[i] * tbe_0 + tr_x_y[i] - 2.0 * tr_x_xxy[i] * tke_0 + tr_xy_0[i] - 2.0 * tr_xy_yy[i] * tke_0 - 2.0 * tr_xy_xx[i] * tke_0 + 4.0 * tr_xy_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xy[i] * tbe_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xz[i] = tr_0_xz[i] - 2.0 * tr_y_xyz[i] * tke_0 - 2.0 * tr_yy_xz[i] * tbe_0 + tr_x_z[i] - 2.0 * tr_x_xxz[i] * tke_0 - 2.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xz[i] * tbe_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_yy[i] = tr_0_yy[i] + 2.0 * tr_y_y[i] - 2.0 * tr_y_yyy[i] * tke_0 - 2.0 * tr_yy_yy[i] * tbe_0 - 2.0 * tr_x_xyy[i] * tke_0 - 4.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yy[i] * tbe_0 - 4.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_yz[i] = tr_0_yz[i] + tr_y_z[i] - 2.0 * tr_y_yyz[i] * tke_0 - 2.0 * tr_yy_yz[i] * tbe_0 - 2.0 * tr_x_xyz[i] * tke_0 - 2.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yz[i] * tbe_0 - 2.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_zz[i] = tr_0_zz[i] - 2.0 * tr_y_yzz[i] * tke_0 - 2.0 * tr_yy_zz[i] * tbe_0 - 2.0 * tr_x_xzz[i] * tke_0 + 4.0 * tr_xy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zz[i] * tbe_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 48-54 components of targeted buffer : DD

    auto tr_0_0_xy_xz_xx = pbuffer.data(idx_op_geom_020_dd + 48);

    auto tr_0_0_xy_xz_xy = pbuffer.data(idx_op_geom_020_dd + 49);

    auto tr_0_0_xy_xz_xz = pbuffer.data(idx_op_geom_020_dd + 50);

    auto tr_0_0_xy_xz_yy = pbuffer.data(idx_op_geom_020_dd + 51);

    auto tr_0_0_xy_xz_yz = pbuffer.data(idx_op_geom_020_dd + 52);

    auto tr_0_0_xy_xz_zz = pbuffer.data(idx_op_geom_020_dd + 53);

    #pragma omp simd aligned(tr_0_0_xy_xz_xx, tr_0_0_xy_xz_xy, tr_0_0_xy_xz_xz, tr_0_0_xy_xz_yy, tr_0_0_xy_xz_yz, tr_0_0_xy_xz_zz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_z_x, tr_z_xxy, tr_z_xyy, tr_z_xyz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xz_xx[i] = -2.0 * tr_z_xxy[i] * tke_0 - 2.0 * tr_yz_xx[i] * tbe_0 - 4.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xy[i] = tr_z_x[i] - 2.0 * tr_z_xyy[i] * tke_0 - 2.0 * tr_yz_xy[i] * tbe_0 + tr_xz_0[i] - 2.0 * tr_xz_yy[i] * tke_0 - 2.0 * tr_xz_xx[i] * tke_0 + 4.0 * tr_xz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xz[i] = -2.0 * tr_z_xyz[i] * tke_0 - 2.0 * tr_yz_xz[i] * tbe_0 - 2.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_yy[i] = 2.0 * tr_z_y[i] - 2.0 * tr_z_yyy[i] * tke_0 - 2.0 * tr_yz_yy[i] * tbe_0 - 4.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_yz[i] = tr_z_z[i] - 2.0 * tr_z_yyz[i] * tke_0 - 2.0 * tr_yz_yz[i] * tbe_0 - 2.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_zz[i] = -2.0 * tr_z_yzz[i] * tke_0 - 2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_xz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 54-60 components of targeted buffer : DD

    auto tr_0_0_xy_yy_xx = pbuffer.data(idx_op_geom_020_dd + 54);

    auto tr_0_0_xy_yy_xy = pbuffer.data(idx_op_geom_020_dd + 55);

    auto tr_0_0_xy_yy_xz = pbuffer.data(idx_op_geom_020_dd + 56);

    auto tr_0_0_xy_yy_yy = pbuffer.data(idx_op_geom_020_dd + 57);

    auto tr_0_0_xy_yy_yz = pbuffer.data(idx_op_geom_020_dd + 58);

    auto tr_0_0_xy_yy_zz = pbuffer.data(idx_op_geom_020_dd + 59);

    #pragma omp simd aligned(tr_0_0_xy_yy_xx, tr_0_0_xy_yy_xy, tr_0_0_xy_yy_xz, tr_0_0_xy_yy_yy, tr_0_0_xy_yy_yz, tr_0_0_xy_yy_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_y_x, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xx, tr_yy_xxxy, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yy_xx[i] = 4.0 * tr_y_x[i] - 4.0 * tr_y_xxx[i] * tke_0 - 4.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyy_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xy[i] = 2.0 * tr_y_y[i] - 4.0 * tr_y_xxy[i] * tke_0 + tr_yy_0[i] - 2.0 * tr_yy_yy[i] * tke_0 - 2.0 * tr_yy_xx[i] * tke_0 + 4.0 * tr_yy_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_yyy_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xy[i] * tbe_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xz[i] = 2.0 * tr_y_z[i] - 4.0 * tr_y_xxz[i] * tke_0 - 2.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyy_z[i] * tbe_0 + 4.0 * tr_yyy_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_yy[i] = -4.0 * tr_y_xyy[i] * tke_0 - 4.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yy[i] * tbe_0 - 4.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_yz[i] = -4.0 * tr_y_xyz[i] * tke_0 - 2.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yz[i] * tbe_0 - 2.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_zz[i] = -4.0 * tr_y_xzz[i] * tke_0 + 4.0 * tr_yy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-66 components of targeted buffer : DD

    auto tr_0_0_xy_yz_xx = pbuffer.data(idx_op_geom_020_dd + 60);

    auto tr_0_0_xy_yz_xy = pbuffer.data(idx_op_geom_020_dd + 61);

    auto tr_0_0_xy_yz_xz = pbuffer.data(idx_op_geom_020_dd + 62);

    auto tr_0_0_xy_yz_yy = pbuffer.data(idx_op_geom_020_dd + 63);

    auto tr_0_0_xy_yz_yz = pbuffer.data(idx_op_geom_020_dd + 64);

    auto tr_0_0_xy_yz_zz = pbuffer.data(idx_op_geom_020_dd + 65);

    #pragma omp simd aligned(tr_0_0_xy_yz_xx, tr_0_0_xy_yz_xy, tr_0_0_xy_yz_xz, tr_0_0_xy_yz_yy, tr_0_0_xy_yz_yz, tr_0_0_xy_yz_zz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yz_xx[i] = 2.0 * tr_z_x[i] - 2.0 * tr_z_xxx[i] * tke_0 - 4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xy[i] = tr_z_y[i] - 2.0 * tr_z_xxy[i] * tke_0 + tr_yz_0[i] - 2.0 * tr_yz_yy[i] * tke_0 - 2.0 * tr_yz_xx[i] * tke_0 + 4.0 * tr_yz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xy[i] * tbe_0 - 2.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xz[i] = tr_z_z[i] - 2.0 * tr_z_xxz[i] * tke_0 - 2.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_yy[i] = -2.0 * tr_z_xyy[i] * tke_0 - 4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yy[i] * tbe_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_yz[i] = -2.0 * tr_z_xyz[i] * tke_0 - 2.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yz[i] * tbe_0 - 2.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_zz[i] = -2.0 * tr_z_xzz[i] * tke_0 + 4.0 * tr_yz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 66-72 components of targeted buffer : DD

    auto tr_0_0_xy_zz_xx = pbuffer.data(idx_op_geom_020_dd + 66);

    auto tr_0_0_xy_zz_xy = pbuffer.data(idx_op_geom_020_dd + 67);

    auto tr_0_0_xy_zz_xz = pbuffer.data(idx_op_geom_020_dd + 68);

    auto tr_0_0_xy_zz_yy = pbuffer.data(idx_op_geom_020_dd + 69);

    auto tr_0_0_xy_zz_yz = pbuffer.data(idx_op_geom_020_dd + 70);

    auto tr_0_0_xy_zz_zz = pbuffer.data(idx_op_geom_020_dd + 71);

    #pragma omp simd aligned(tr_0_0_xy_zz_xx, tr_0_0_xy_zz_xy, tr_0_0_xy_zz_xz, tr_0_0_xy_zz_yy, tr_0_0_xy_zz_yz, tr_0_0_xy_zz_zz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_z, tr_zz_0, tr_zz_xx, tr_zz_xxxy, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_zz_xx[i] = -4.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xy[i] = tr_zz_0[i] - 2.0 * tr_zz_yy[i] * tke_0 - 2.0 * tr_zz_xx[i] * tke_0 + 4.0 * tr_zz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xz[i] = -2.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_yy[i] = -4.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_yz[i] = -2.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_zz[i] = 4.0 * tr_zz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 72-78 components of targeted buffer : DD

    auto tr_0_0_xz_xx_xx = pbuffer.data(idx_op_geom_020_dd + 72);

    auto tr_0_0_xz_xx_xy = pbuffer.data(idx_op_geom_020_dd + 73);

    auto tr_0_0_xz_xx_xz = pbuffer.data(idx_op_geom_020_dd + 74);

    auto tr_0_0_xz_xx_yy = pbuffer.data(idx_op_geom_020_dd + 75);

    auto tr_0_0_xz_xx_yz = pbuffer.data(idx_op_geom_020_dd + 76);

    auto tr_0_0_xz_xx_zz = pbuffer.data(idx_op_geom_020_dd + 77);

    #pragma omp simd aligned(tr_0_0_xz_xx_xx, tr_0_0_xz_xx_xy, tr_0_0_xz_xx_xz, tr_0_0_xz_xx_yy, tr_0_0_xz_xx_yz, tr_0_0_xz_xx_zz, tr_x_x, tr_x_xxz, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xx_0, tr_xx_xx, tr_xx_xxxz, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xx_xx[i] = -4.0 * tr_x_xxz[i] * tke_0 - 4.0 * tr_xz_xx[i] * tbe_0 - 4.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xy[i] = -4.0 * tr_x_xyz[i] * tke_0 - 4.0 * tr_xz_xy[i] * tbe_0 - 2.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xz[i] = 2.0 * tr_x_x[i] - 4.0 * tr_x_xzz[i] * tke_0 - 4.0 * tr_xz_xz[i] * tbe_0 + tr_xx_0[i] - 2.0 * tr_xx_zz[i] * tke_0 - 2.0 * tr_xx_xx[i] * tke_0 + 4.0 * tr_xx_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxx_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_yy[i] = -4.0 * tr_x_yyz[i] * tke_0 - 4.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xx_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_yz[i] = 2.0 * tr_x_y[i] - 4.0 * tr_x_yzz[i] * tke_0 - 4.0 * tr_xz_yz[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxx_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_zz[i] = 4.0 * tr_x_z[i] - 4.0 * tr_x_zzz[i] * tke_0 - 4.0 * tr_xz_zz[i] * tbe_0 - 4.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxx_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 78-84 components of targeted buffer : DD

    auto tr_0_0_xz_xy_xx = pbuffer.data(idx_op_geom_020_dd + 78);

    auto tr_0_0_xz_xy_xy = pbuffer.data(idx_op_geom_020_dd + 79);

    auto tr_0_0_xz_xy_xz = pbuffer.data(idx_op_geom_020_dd + 80);

    auto tr_0_0_xz_xy_yy = pbuffer.data(idx_op_geom_020_dd + 81);

    auto tr_0_0_xz_xy_yz = pbuffer.data(idx_op_geom_020_dd + 82);

    auto tr_0_0_xz_xy_zz = pbuffer.data(idx_op_geom_020_dd + 83);

    #pragma omp simd aligned(tr_0_0_xz_xy_xx, tr_0_0_xz_xy_xy, tr_0_0_xz_xy_xz, tr_0_0_xz_xy_yy, tr_0_0_xz_xy_yz, tr_0_0_xz_xy_zz, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xy_0, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_y_x, tr_y_xxz, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xy_xx[i] = -2.0 * tr_y_xxz[i] * tke_0 - 2.0 * tr_yz_xx[i] * tbe_0 - 4.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xy[i] = -2.0 * tr_y_xyz[i] * tke_0 - 2.0 * tr_yz_xy[i] * tbe_0 - 2.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xz[i] = tr_y_x[i] - 2.0 * tr_y_xzz[i] * tke_0 - 2.0 * tr_yz_xz[i] * tbe_0 + tr_xy_0[i] - 2.0 * tr_xy_zz[i] * tke_0 - 2.0 * tr_xy_xx[i] * tke_0 + 4.0 * tr_xy_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_yy[i] = -2.0 * tr_y_yyz[i] * tke_0 - 2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_xy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_yz[i] = tr_y_y[i] - 2.0 * tr_y_yzz[i] * tke_0 - 2.0 * tr_yz_yz[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_zz[i] = 2.0 * tr_y_z[i] - 2.0 * tr_y_zzz[i] * tke_0 - 2.0 * tr_yz_zz[i] * tbe_0 - 4.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 84-90 components of targeted buffer : DD

    auto tr_0_0_xz_xz_xx = pbuffer.data(idx_op_geom_020_dd + 84);

    auto tr_0_0_xz_xz_xy = pbuffer.data(idx_op_geom_020_dd + 85);

    auto tr_0_0_xz_xz_xz = pbuffer.data(idx_op_geom_020_dd + 86);

    auto tr_0_0_xz_xz_yy = pbuffer.data(idx_op_geom_020_dd + 87);

    auto tr_0_0_xz_xz_yz = pbuffer.data(idx_op_geom_020_dd + 88);

    auto tr_0_0_xz_xz_zz = pbuffer.data(idx_op_geom_020_dd + 89);

    #pragma omp simd aligned(tr_0_0_xz_xz_xx, tr_0_0_xz_xz_xy, tr_0_0_xz_xz_xz, tr_0_0_xz_xz_yy, tr_0_0_xz_xz_yz, tr_0_0_xz_xz_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_x, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_z, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xz_0, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_z, tr_z_x, tr_z_xxz, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xz_xx[i] = tr_0_xx[i] - 2.0 * tr_z_xxz[i] * tke_0 - 2.0 * tr_zz_xx[i] * tbe_0 + 2.0 * tr_x_x[i] - 2.0 * tr_x_xxx[i] * tke_0 - 4.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xy[i] = tr_0_xy[i] - 2.0 * tr_z_xyz[i] * tke_0 - 2.0 * tr_zz_xy[i] * tbe_0 + tr_x_y[i] - 2.0 * tr_x_xxy[i] * tke_0 - 2.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xy[i] * tbe_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xz[i] = tr_0_xz[i] + tr_z_x[i] - 2.0 * tr_z_xzz[i] * tke_0 - 2.0 * tr_zz_xz[i] * tbe_0 + tr_x_z[i] - 2.0 * tr_x_xxz[i] * tke_0 + tr_xz_0[i] - 2.0 * tr_xz_zz[i] * tke_0 - 2.0 * tr_xz_xx[i] * tke_0 + 4.0 * tr_xz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xz[i] * tbe_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_yy[i] = tr_0_yy[i] - 2.0 * tr_z_yyz[i] * tke_0 - 2.0 * tr_zz_yy[i] * tbe_0 - 2.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_xz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yy[i] * tbe_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_yz[i] = tr_0_yz[i] + tr_z_y[i] - 2.0 * tr_z_yzz[i] * tke_0 - 2.0 * tr_zz_yz[i] * tbe_0 - 2.0 * tr_x_xyz[i] * tke_0 - 2.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yz[i] * tbe_0 - 2.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_zz[i] = tr_0_zz[i] + 2.0 * tr_z_z[i] - 2.0 * tr_z_zzz[i] * tke_0 - 2.0 * tr_zz_zz[i] * tbe_0 - 2.0 * tr_x_xzz[i] * tke_0 - 4.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zz[i] * tbe_0 - 4.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-96 components of targeted buffer : DD

    auto tr_0_0_xz_yy_xx = pbuffer.data(idx_op_geom_020_dd + 90);

    auto tr_0_0_xz_yy_xy = pbuffer.data(idx_op_geom_020_dd + 91);

    auto tr_0_0_xz_yy_xz = pbuffer.data(idx_op_geom_020_dd + 92);

    auto tr_0_0_xz_yy_yy = pbuffer.data(idx_op_geom_020_dd + 93);

    auto tr_0_0_xz_yy_yz = pbuffer.data(idx_op_geom_020_dd + 94);

    auto tr_0_0_xz_yy_zz = pbuffer.data(idx_op_geom_020_dd + 95);

    #pragma omp simd aligned(tr_0_0_xz_yy_xx, tr_0_0_xz_yy_xy, tr_0_0_xz_yy_xz, tr_0_0_xz_yy_yy, tr_0_0_xz_yy_yz, tr_0_0_xz_yy_zz, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yy_0, tr_yy_xx, tr_yy_xxxz, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yy_xx[i] = -4.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xy[i] = -2.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xz[i] = tr_yy_0[i] - 2.0 * tr_yy_zz[i] * tke_0 - 2.0 * tr_yy_xx[i] * tke_0 + 4.0 * tr_yy_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_yy[i] = 4.0 * tr_yy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_yz[i] = -2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_zz[i] = -4.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 96-102 components of targeted buffer : DD

    auto tr_0_0_xz_yz_xx = pbuffer.data(idx_op_geom_020_dd + 96);

    auto tr_0_0_xz_yz_xy = pbuffer.data(idx_op_geom_020_dd + 97);

    auto tr_0_0_xz_yz_xz = pbuffer.data(idx_op_geom_020_dd + 98);

    auto tr_0_0_xz_yz_yy = pbuffer.data(idx_op_geom_020_dd + 99);

    auto tr_0_0_xz_yz_yz = pbuffer.data(idx_op_geom_020_dd + 100);

    auto tr_0_0_xz_yz_zz = pbuffer.data(idx_op_geom_020_dd + 101);

    #pragma omp simd aligned(tr_0_0_xz_yz_xx, tr_0_0_xz_yz_xy, tr_0_0_xz_yz_xz, tr_0_0_xz_yz_yy, tr_0_0_xz_yz_yz, tr_0_0_xz_yz_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_y_x, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yz_xx[i] = 2.0 * tr_y_x[i] - 2.0 * tr_y_xxx[i] * tke_0 - 4.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xy[i] = tr_y_y[i] - 2.0 * tr_y_xxy[i] * tke_0 - 2.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xz[i] = tr_y_z[i] - 2.0 * tr_y_xxz[i] * tke_0 + tr_yz_0[i] - 2.0 * tr_yz_zz[i] * tke_0 - 2.0 * tr_yz_xx[i] * tke_0 + 4.0 * tr_yz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xz[i] * tbe_0 - 2.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_yy[i] = -2.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_yz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_yz[i] = -2.0 * tr_y_xyz[i] * tke_0 - 2.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yz[i] * tbe_0 - 2.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_zz[i] = -2.0 * tr_y_xzz[i] * tke_0 - 4.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_zz[i] * tbe_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 102-108 components of targeted buffer : DD

    auto tr_0_0_xz_zz_xx = pbuffer.data(idx_op_geom_020_dd + 102);

    auto tr_0_0_xz_zz_xy = pbuffer.data(idx_op_geom_020_dd + 103);

    auto tr_0_0_xz_zz_xz = pbuffer.data(idx_op_geom_020_dd + 104);

    auto tr_0_0_xz_zz_yy = pbuffer.data(idx_op_geom_020_dd + 105);

    auto tr_0_0_xz_zz_yz = pbuffer.data(idx_op_geom_020_dd + 106);

    auto tr_0_0_xz_zz_zz = pbuffer.data(idx_op_geom_020_dd + 107);

    #pragma omp simd aligned(tr_0_0_xz_zz_xx, tr_0_0_xz_zz_xy, tr_0_0_xz_zz_xz, tr_0_0_xz_zz_yy, tr_0_0_xz_zz_yz, tr_0_0_xz_zz_zz, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_z_x, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xx, tr_zz_xxxz, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_zz_xx[i] = 4.0 * tr_z_x[i] - 4.0 * tr_z_xxx[i] * tke_0 - 4.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_x[i] * tbe_0 + 4.0 * tr_zzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xy[i] = 2.0 * tr_z_y[i] - 4.0 * tr_z_xxy[i] * tke_0 - 2.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_y[i] * tbe_0 + 4.0 * tr_zzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xz[i] = 2.0 * tr_z_z[i] - 4.0 * tr_z_xxz[i] * tke_0 + tr_zz_0[i] - 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_zz_xx[i] * tke_0 + 4.0 * tr_zz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_z[i] * tbe_0 + 4.0 * tr_zzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xz[i] * tbe_0 - 2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_yy[i] = -4.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_zz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_yz[i] = -4.0 * tr_z_xyz[i] * tke_0 - 2.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yz[i] * tbe_0 - 2.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_zz[i] = -4.0 * tr_z_xzz[i] * tke_0 - 4.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_zz[i] * tbe_0 - 4.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 108-114 components of targeted buffer : DD

    auto tr_0_0_yy_xx_xx = pbuffer.data(idx_op_geom_020_dd + 108);

    auto tr_0_0_yy_xx_xy = pbuffer.data(idx_op_geom_020_dd + 109);

    auto tr_0_0_yy_xx_xz = pbuffer.data(idx_op_geom_020_dd + 110);

    auto tr_0_0_yy_xx_yy = pbuffer.data(idx_op_geom_020_dd + 111);

    auto tr_0_0_yy_xx_yz = pbuffer.data(idx_op_geom_020_dd + 112);

    auto tr_0_0_yy_xx_zz = pbuffer.data(idx_op_geom_020_dd + 113);

    #pragma omp simd aligned(tr_0_0_yy_xx_xx, tr_0_0_yy_xx_xy, tr_0_0_yy_xx_xz, tr_0_0_yy_xx_yy, tr_0_0_yy_xx_yz, tr_0_0_yy_xx_zz, tr_xx_0, tr_xx_xx, tr_xx_xxyy, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xz, tr_xx_yy, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_zz, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xx_xx[i] = -2.0 * tr_xx_xx[i] * tbe_0 - 2.0 * tr_xx_xx[i] * tke_0 + 4.0 * tr_xx_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xy[i] = -2.0 * tr_xx_xy[i] * tbe_0 - 6.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxy_x[i] * tbe_0 + 8.0 * tr_xxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xz[i] = -2.0 * tr_xx_xz[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_yy[i] = 2.0 * tr_xx_0[i] - 2.0 * tr_xx_yy[i] * tbe_0 - 10.0 * tr_xx_yy[i] * tke_0 + 4.0 * tr_xx_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxy_y[i] * tbe_0 + 8.0 * tr_xxy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_yz[i] = -2.0 * tr_xx_yz[i] * tbe_0 - 6.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 8.0 * tr_xxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_zz[i] = -2.0 * tr_xx_zz[i] * tbe_0 - 2.0 * tr_xx_zz[i] * tke_0 + 4.0 * tr_xx_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 114-120 components of targeted buffer : DD

    auto tr_0_0_yy_xy_xx = pbuffer.data(idx_op_geom_020_dd + 114);

    auto tr_0_0_yy_xy_xy = pbuffer.data(idx_op_geom_020_dd + 115);

    auto tr_0_0_yy_xy_xz = pbuffer.data(idx_op_geom_020_dd + 116);

    auto tr_0_0_yy_xy_yy = pbuffer.data(idx_op_geom_020_dd + 117);

    auto tr_0_0_yy_xy_yz = pbuffer.data(idx_op_geom_020_dd + 118);

    auto tr_0_0_yy_xy_zz = pbuffer.data(idx_op_geom_020_dd + 119);

    #pragma omp simd aligned(tr_0_0_yy_xy_xx, tr_0_0_yy_xy_xy, tr_0_0_yy_xy_xz, tr_0_0_yy_xy_yy, tr_0_0_yy_xy_yz, tr_0_0_yy_xy_zz, tr_x_x, tr_x_xxy, tr_x_xyy, tr_x_xyz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_xy_0, tr_xy_xx, tr_xy_xxyy, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xy_xx[i] = -4.0 * tr_x_xxy[i] * tke_0 - 6.0 * tr_xy_xx[i] * tbe_0 - 2.0 * tr_xy_xx[i] * tke_0 + 4.0 * tr_xy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xy[i] = 2.0 * tr_x_x[i] - 4.0 * tr_x_xyy[i] * tke_0 - 6.0 * tr_xy_xy[i] * tbe_0 - 6.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyy_x[i] * tbe_0 + 8.0 * tr_xyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xz[i] = -4.0 * tr_x_xyz[i] * tke_0 - 6.0 * tr_xy_xz[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_yy[i] = 4.0 * tr_x_y[i] - 4.0 * tr_x_yyy[i] * tke_0 + 2.0 * tr_xy_0[i] - 6.0 * tr_xy_yy[i] * tbe_0 - 10.0 * tr_xy_yy[i] * tke_0 + 4.0 * tr_xy_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyy_y[i] * tbe_0 + 8.0 * tr_xyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_yz[i] = 2.0 * tr_x_z[i] - 4.0 * tr_x_yyz[i] * tke_0 - 6.0 * tr_xy_yz[i] * tbe_0 - 6.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_z[i] * tbe_0 + 8.0 * tr_xyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_zz[i] = -4.0 * tr_x_yzz[i] * tke_0 - 6.0 * tr_xy_zz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tke_0 + 4.0 * tr_xy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-126 components of targeted buffer : DD

    auto tr_0_0_yy_xz_xx = pbuffer.data(idx_op_geom_020_dd + 120);

    auto tr_0_0_yy_xz_xy = pbuffer.data(idx_op_geom_020_dd + 121);

    auto tr_0_0_yy_xz_xz = pbuffer.data(idx_op_geom_020_dd + 122);

    auto tr_0_0_yy_xz_yy = pbuffer.data(idx_op_geom_020_dd + 123);

    auto tr_0_0_yy_xz_yz = pbuffer.data(idx_op_geom_020_dd + 124);

    auto tr_0_0_yy_xz_zz = pbuffer.data(idx_op_geom_020_dd + 125);

    #pragma omp simd aligned(tr_0_0_yy_xz_xx, tr_0_0_yy_xz_xy, tr_0_0_yy_xz_xz, tr_0_0_yy_xz_yy, tr_0_0_yy_xz_yz, tr_0_0_yy_xz_zz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_0, tr_xz_xx, tr_xz_xxyy, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xz_xx[i] = -2.0 * tr_xz_xx[i] * tbe_0 - 2.0 * tr_xz_xx[i] * tke_0 + 4.0 * tr_xz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xy[i] = -2.0 * tr_xz_xy[i] * tbe_0 - 6.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 8.0 * tr_xyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xz[i] = -2.0 * tr_xz_xz[i] * tbe_0 - 2.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_yy[i] = 2.0 * tr_xz_0[i] - 2.0 * tr_xz_yy[i] * tbe_0 - 10.0 * tr_xz_yy[i] * tke_0 + 4.0 * tr_xz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyz_y[i] * tbe_0 + 8.0 * tr_xyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_yz[i] = -2.0 * tr_xz_yz[i] * tbe_0 - 6.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 8.0 * tr_xyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_zz[i] = -2.0 * tr_xz_zz[i] * tbe_0 - 2.0 * tr_xz_zz[i] * tke_0 + 4.0 * tr_xz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 126-132 components of targeted buffer : DD

    auto tr_0_0_yy_yy_xx = pbuffer.data(idx_op_geom_020_dd + 126);

    auto tr_0_0_yy_yy_xy = pbuffer.data(idx_op_geom_020_dd + 127);

    auto tr_0_0_yy_yy_xz = pbuffer.data(idx_op_geom_020_dd + 128);

    auto tr_0_0_yy_yy_yy = pbuffer.data(idx_op_geom_020_dd + 129);

    auto tr_0_0_yy_yy_yz = pbuffer.data(idx_op_geom_020_dd + 130);

    auto tr_0_0_yy_yy_zz = pbuffer.data(idx_op_geom_020_dd + 131);

    #pragma omp simd aligned(tr_0_0_yy_yy_xx, tr_0_0_yy_yy_xy, tr_0_0_yy_yy_xz, tr_0_0_yy_yy_yy, tr_0_0_yy_yy_yz, tr_0_0_yy_yy_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_y_x, tr_y_xxy, tr_y_xyy, tr_y_xyz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_yy_0, tr_yy_xx, tr_yy_xxyy, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xz, tr_yy_yy, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_zz, tr_yyy_x, tr_yyy_xxy, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yy_xx[i] = 2.0 * tr_0_xx[i] - 8.0 * tr_y_xxy[i] * tke_0 - 10.0 * tr_yy_xx[i] * tbe_0 - 2.0 * tr_yy_xx[i] * tke_0 + 4.0 * tr_yy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xy[i] = 2.0 * tr_0_xy[i] + 4.0 * tr_y_x[i] - 8.0 * tr_y_xyy[i] * tke_0 - 10.0 * tr_yy_xy[i] * tbe_0 - 6.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyy_x[i] * tbe_0 + 8.0 * tr_yyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xz[i] = 2.0 * tr_0_xz[i] - 8.0 * tr_y_xyz[i] * tke_0 - 10.0 * tr_yy_xz[i] * tbe_0 - 2.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_yy[i] = 2.0 * tr_0_yy[i] + 8.0 * tr_y_y[i] - 8.0 * tr_y_yyy[i] * tke_0 + 2.0 * tr_yy_0[i] - 10.0 * tr_yy_yy[i] * tbe_0 - 10.0 * tr_yy_yy[i] * tke_0 + 4.0 * tr_yy_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyy_y[i] * tbe_0 + 8.0 * tr_yyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_yz[i] = 2.0 * tr_0_yz[i] + 4.0 * tr_y_z[i] - 8.0 * tr_y_yyz[i] * tke_0 - 10.0 * tr_yy_yz[i] * tbe_0 - 6.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyy_z[i] * tbe_0 + 8.0 * tr_yyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_zz[i] = 2.0 * tr_0_zz[i] - 8.0 * tr_y_yzz[i] * tke_0 - 10.0 * tr_yy_zz[i] * tbe_0 - 2.0 * tr_yy_zz[i] * tke_0 + 4.0 * tr_yy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 132-138 components of targeted buffer : DD

    auto tr_0_0_yy_yz_xx = pbuffer.data(idx_op_geom_020_dd + 132);

    auto tr_0_0_yy_yz_xy = pbuffer.data(idx_op_geom_020_dd + 133);

    auto tr_0_0_yy_yz_xz = pbuffer.data(idx_op_geom_020_dd + 134);

    auto tr_0_0_yy_yz_yy = pbuffer.data(idx_op_geom_020_dd + 135);

    auto tr_0_0_yy_yz_yz = pbuffer.data(idx_op_geom_020_dd + 136);

    auto tr_0_0_yy_yz_zz = pbuffer.data(idx_op_geom_020_dd + 137);

    #pragma omp simd aligned(tr_0_0_yy_yz_xx, tr_0_0_yy_yz_xy, tr_0_0_yy_yz_xz, tr_0_0_yy_yz_yy, tr_0_0_yy_yz_yz, tr_0_0_yy_yz_zz, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yz_0, tr_yz_xx, tr_yz_xxyy, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_zz, tr_z_x, tr_z_xxy, tr_z_xyy, tr_z_xyz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yz_xx[i] = -4.0 * tr_z_xxy[i] * tke_0 - 6.0 * tr_yz_xx[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tke_0 + 4.0 * tr_yz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xy[i] = 2.0 * tr_z_x[i] - 4.0 * tr_z_xyy[i] * tke_0 - 6.0 * tr_yz_xy[i] * tbe_0 - 6.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 8.0 * tr_yyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xz[i] = -4.0 * tr_z_xyz[i] * tke_0 - 6.0 * tr_yz_xz[i] * tbe_0 - 2.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_yy[i] = 4.0 * tr_z_y[i] - 4.0 * tr_z_yyy[i] * tke_0 + 2.0 * tr_yz_0[i] - 6.0 * tr_yz_yy[i] * tbe_0 - 10.0 * tr_yz_yy[i] * tke_0 + 4.0 * tr_yz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyz_y[i] * tbe_0 + 8.0 * tr_yyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_yz[i] = 2.0 * tr_z_z[i] - 4.0 * tr_z_yyz[i] * tke_0 - 6.0 * tr_yz_yz[i] * tbe_0 - 6.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_z[i] * tbe_0 + 8.0 * tr_yyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_zz[i] = -4.0 * tr_z_yzz[i] * tke_0 - 6.0 * tr_yz_zz[i] * tbe_0 - 2.0 * tr_yz_zz[i] * tke_0 + 4.0 * tr_yz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 138-144 components of targeted buffer : DD

    auto tr_0_0_yy_zz_xx = pbuffer.data(idx_op_geom_020_dd + 138);

    auto tr_0_0_yy_zz_xy = pbuffer.data(idx_op_geom_020_dd + 139);

    auto tr_0_0_yy_zz_xz = pbuffer.data(idx_op_geom_020_dd + 140);

    auto tr_0_0_yy_zz_yy = pbuffer.data(idx_op_geom_020_dd + 141);

    auto tr_0_0_yy_zz_yz = pbuffer.data(idx_op_geom_020_dd + 142);

    auto tr_0_0_yy_zz_zz = pbuffer.data(idx_op_geom_020_dd + 143);

    #pragma omp simd aligned(tr_0_0_yy_zz_xx, tr_0_0_yy_zz_xy, tr_0_0_yy_zz_xz, tr_0_0_yy_zz_yy, tr_0_0_yy_zz_yz, tr_0_0_yy_zz_zz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_zz_0, tr_zz_xx, tr_zz_xxyy, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xz, tr_zz_yy, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_zz_xx[i] = -2.0 * tr_zz_xx[i] * tbe_0 - 2.0 * tr_zz_xx[i] * tke_0 + 4.0 * tr_zz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xy[i] = -2.0 * tr_zz_xy[i] * tbe_0 - 6.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yzz_x[i] * tbe_0 + 8.0 * tr_yzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xz[i] = -2.0 * tr_zz_xz[i] * tbe_0 - 2.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_yy[i] = 2.0 * tr_zz_0[i] - 2.0 * tr_zz_yy[i] * tbe_0 - 10.0 * tr_zz_yy[i] * tke_0 + 4.0 * tr_zz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yzz_y[i] * tbe_0 + 8.0 * tr_yzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_yz[i] = -2.0 * tr_zz_yz[i] * tbe_0 - 6.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_z[i] * tbe_0 + 8.0 * tr_yzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_zz[i] = -2.0 * tr_zz_zz[i] * tbe_0 - 2.0 * tr_zz_zz[i] * tke_0 + 4.0 * tr_zz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 144-150 components of targeted buffer : DD

    auto tr_0_0_yz_xx_xx = pbuffer.data(idx_op_geom_020_dd + 144);

    auto tr_0_0_yz_xx_xy = pbuffer.data(idx_op_geom_020_dd + 145);

    auto tr_0_0_yz_xx_xz = pbuffer.data(idx_op_geom_020_dd + 146);

    auto tr_0_0_yz_xx_yy = pbuffer.data(idx_op_geom_020_dd + 147);

    auto tr_0_0_yz_xx_yz = pbuffer.data(idx_op_geom_020_dd + 148);

    auto tr_0_0_yz_xx_zz = pbuffer.data(idx_op_geom_020_dd + 149);

    #pragma omp simd aligned(tr_0_0_yz_xx_xx, tr_0_0_yz_xx_xy, tr_0_0_yz_xx_xz, tr_0_0_yz_xx_yy, tr_0_0_yz_xx_yz, tr_0_0_yz_xx_zz, tr_xx_0, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xx_xx[i] = 4.0 * tr_xx_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xy[i] = -2.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xz[i] = -2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_yy[i] = -4.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_yz[i] = tr_xx_0[i] - 2.0 * tr_xx_zz[i] * tke_0 - 2.0 * tr_xx_yy[i] * tke_0 + 4.0 * tr_xx_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_zz[i] = -4.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-156 components of targeted buffer : DD

    auto tr_0_0_yz_xy_xx = pbuffer.data(idx_op_geom_020_dd + 150);

    auto tr_0_0_yz_xy_xy = pbuffer.data(idx_op_geom_020_dd + 151);

    auto tr_0_0_yz_xy_xz = pbuffer.data(idx_op_geom_020_dd + 152);

    auto tr_0_0_yz_xy_yy = pbuffer.data(idx_op_geom_020_dd + 153);

    auto tr_0_0_yz_xy_yz = pbuffer.data(idx_op_geom_020_dd + 154);

    auto tr_0_0_yz_xy_zz = pbuffer.data(idx_op_geom_020_dd + 155);

    #pragma omp simd aligned(tr_0_0_yz_xy_xx, tr_0_0_yz_xy_xy, tr_0_0_yz_xy_xz, tr_0_0_yz_xy_yy, tr_0_0_yz_xy_yz, tr_0_0_yz_xy_zz, tr_x_x, tr_x_xxz, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xy_0, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xy_xx[i] = -2.0 * tr_x_xxz[i] * tke_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xy[i] = -2.0 * tr_x_xyz[i] * tke_0 - 2.0 * tr_xz_xy[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xz[i] = tr_x_x[i] - 2.0 * tr_x_xzz[i] * tke_0 - 2.0 * tr_xz_xz[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_yy[i] = -2.0 * tr_x_yyz[i] * tke_0 - 2.0 * tr_xz_yy[i] * tbe_0 - 4.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_yz[i] = tr_x_y[i] - 2.0 * tr_x_yzz[i] * tke_0 - 2.0 * tr_xz_yz[i] * tbe_0 + tr_xy_0[i] - 2.0 * tr_xy_zz[i] * tke_0 - 2.0 * tr_xy_yy[i] * tke_0 + 4.0 * tr_xy_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_zz[i] = 2.0 * tr_x_z[i] - 2.0 * tr_x_zzz[i] * tke_0 - 2.0 * tr_xz_zz[i] * tbe_0 - 4.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 156-162 components of targeted buffer : DD

    auto tr_0_0_yz_xz_xx = pbuffer.data(idx_op_geom_020_dd + 156);

    auto tr_0_0_yz_xz_xy = pbuffer.data(idx_op_geom_020_dd + 157);

    auto tr_0_0_yz_xz_xz = pbuffer.data(idx_op_geom_020_dd + 158);

    auto tr_0_0_yz_xz_yy = pbuffer.data(idx_op_geom_020_dd + 159);

    auto tr_0_0_yz_xz_yz = pbuffer.data(idx_op_geom_020_dd + 160);

    auto tr_0_0_yz_xz_zz = pbuffer.data(idx_op_geom_020_dd + 161);

    #pragma omp simd aligned(tr_0_0_yz_xz_xx, tr_0_0_yz_xz_xy, tr_0_0_yz_xz_xz, tr_0_0_yz_xz_yy, tr_0_0_yz_xz_yz, tr_0_0_yz_xz_zz, tr_x_x, tr_x_xxy, tr_x_xyy, tr_x_xyz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xz_0, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xz_xx[i] = -2.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_xz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xy[i] = tr_x_x[i] - 2.0 * tr_x_xyy[i] * tke_0 - 2.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xz[i] = -2.0 * tr_x_xyz[i] * tke_0 - 2.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xz[i] * tbe_0 - 2.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_yy[i] = 2.0 * tr_x_y[i] - 2.0 * tr_x_yyy[i] * tke_0 - 4.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_yz[i] = tr_x_z[i] - 2.0 * tr_x_yyz[i] * tke_0 + tr_xz_0[i] - 2.0 * tr_xz_zz[i] * tke_0 - 2.0 * tr_xz_yy[i] * tke_0 + 4.0 * tr_xz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yz[i] * tbe_0 - 2.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_zz[i] = -2.0 * tr_x_yzz[i] * tke_0 - 4.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_zz[i] * tbe_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 162-168 components of targeted buffer : DD

    auto tr_0_0_yz_yy_xx = pbuffer.data(idx_op_geom_020_dd + 162);

    auto tr_0_0_yz_yy_xy = pbuffer.data(idx_op_geom_020_dd + 163);

    auto tr_0_0_yz_yy_xz = pbuffer.data(idx_op_geom_020_dd + 164);

    auto tr_0_0_yz_yy_yy = pbuffer.data(idx_op_geom_020_dd + 165);

    auto tr_0_0_yz_yy_yz = pbuffer.data(idx_op_geom_020_dd + 166);

    auto tr_0_0_yz_yy_zz = pbuffer.data(idx_op_geom_020_dd + 167);

    #pragma omp simd aligned(tr_0_0_yz_yy_xx, tr_0_0_yz_yy_xy, tr_0_0_yz_yy_xz, tr_0_0_yz_yy_yy, tr_0_0_yz_yy_yz, tr_0_0_yz_yy_zz, tr_y_x, tr_y_xxz, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yy_0, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yyy_x, tr_yyy_xxz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yy_xx[i] = -4.0 * tr_y_xxz[i] * tke_0 - 4.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yy_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xy[i] = -4.0 * tr_y_xyz[i] * tke_0 - 4.0 * tr_yz_xy[i] * tbe_0 - 2.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xz[i] = 2.0 * tr_y_x[i] - 4.0 * tr_y_xzz[i] * tke_0 - 4.0 * tr_yz_xz[i] * tbe_0 - 2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_yy[i] = -4.0 * tr_y_yyz[i] * tke_0 - 4.0 * tr_yz_yy[i] * tbe_0 - 4.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_yz[i] = 2.0 * tr_y_y[i] - 4.0 * tr_y_yzz[i] * tke_0 - 4.0 * tr_yz_yz[i] * tbe_0 + tr_yy_0[i] - 2.0 * tr_yy_zz[i] * tke_0 - 2.0 * tr_yy_yy[i] * tke_0 + 4.0 * tr_yy_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_yyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_zz[i] = 4.0 * tr_y_z[i] - 4.0 * tr_y_zzz[i] * tke_0 - 4.0 * tr_yz_zz[i] * tbe_0 - 4.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyy_z[i] * tbe_0 + 4.0 * tr_yyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 168-174 components of targeted buffer : DD

    auto tr_0_0_yz_yz_xx = pbuffer.data(idx_op_geom_020_dd + 168);

    auto tr_0_0_yz_yz_xy = pbuffer.data(idx_op_geom_020_dd + 169);

    auto tr_0_0_yz_yz_xz = pbuffer.data(idx_op_geom_020_dd + 170);

    auto tr_0_0_yz_yz_yy = pbuffer.data(idx_op_geom_020_dd + 171);

    auto tr_0_0_yz_yz_yz = pbuffer.data(idx_op_geom_020_dd + 172);

    auto tr_0_0_yz_yz_zz = pbuffer.data(idx_op_geom_020_dd + 173);

    #pragma omp simd aligned(tr_0_0_yz_yz_xx, tr_0_0_yz_yz_xy, tr_0_0_yz_yz_xz, tr_0_0_yz_yz_yy, tr_0_0_yz_yz_yz, tr_0_0_yz_yz_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_y_x, tr_y_xxy, tr_y_xyy, tr_y_xyz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yz_0, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_z_x, tr_z_xxz, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yz_xx[i] = tr_0_xx[i] - 2.0 * tr_z_xxz[i] * tke_0 - 2.0 * tr_zz_xx[i] * tbe_0 - 2.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_yz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xx[i] * tbe_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xy[i] = tr_0_xy[i] - 2.0 * tr_z_xyz[i] * tke_0 - 2.0 * tr_zz_xy[i] * tbe_0 + tr_y_x[i] - 2.0 * tr_y_xyy[i] * tke_0 - 2.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xy[i] * tbe_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xz[i] = tr_0_xz[i] + tr_z_x[i] - 2.0 * tr_z_xzz[i] * tke_0 - 2.0 * tr_zz_xz[i] * tbe_0 - 2.0 * tr_y_xyz[i] * tke_0 - 2.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xz[i] * tbe_0 - 2.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_yy[i] = tr_0_yy[i] - 2.0 * tr_z_yyz[i] * tke_0 - 2.0 * tr_zz_yy[i] * tbe_0 + 2.0 * tr_y_y[i] - 2.0 * tr_y_yyy[i] * tke_0 - 4.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yy[i] * tbe_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_yz[i] = tr_0_yz[i] + tr_z_y[i] - 2.0 * tr_z_yzz[i] * tke_0 - 2.0 * tr_zz_yz[i] * tbe_0 + tr_y_z[i] - 2.0 * tr_y_yyz[i] * tke_0 + tr_yz_0[i] - 2.0 * tr_yz_zz[i] * tke_0 - 2.0 * tr_yz_yy[i] * tke_0 + 4.0 * tr_yz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yz[i] * tbe_0 - 2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_zz[i] = tr_0_zz[i] + 2.0 * tr_z_z[i] - 2.0 * tr_z_zzz[i] * tke_0 - 2.0 * tr_zz_zz[i] * tbe_0 - 2.0 * tr_y_yzz[i] * tke_0 - 4.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_zz[i] * tbe_0 - 4.0 * tr_yyz_z[i] * tbe_0 + 4.0 * tr_yyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 174-180 components of targeted buffer : DD

    auto tr_0_0_yz_zz_xx = pbuffer.data(idx_op_geom_020_dd + 174);

    auto tr_0_0_yz_zz_xy = pbuffer.data(idx_op_geom_020_dd + 175);

    auto tr_0_0_yz_zz_xz = pbuffer.data(idx_op_geom_020_dd + 176);

    auto tr_0_0_yz_zz_yy = pbuffer.data(idx_op_geom_020_dd + 177);

    auto tr_0_0_yz_zz_yz = pbuffer.data(idx_op_geom_020_dd + 178);

    auto tr_0_0_yz_zz_zz = pbuffer.data(idx_op_geom_020_dd + 179);

    #pragma omp simd aligned(tr_0_0_yz_zz_xx, tr_0_0_yz_zz_xy, tr_0_0_yz_zz_xz, tr_0_0_yz_zz_yy, tr_0_0_yz_zz_yz, tr_0_0_yz_zz_zz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, tr_z_x, tr_z_xxy, tr_z_xyy, tr_z_xyz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, tr_zz_0, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zzz_x, tr_zzz_xxy, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_zz_xx[i] = -4.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_zz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xy[i] = 2.0 * tr_z_x[i] - 4.0 * tr_z_xyy[i] * tke_0 - 2.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_x[i] * tbe_0 + 4.0 * tr_zzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xz[i] = -4.0 * tr_z_xyz[i] * tke_0 - 2.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xz[i] * tbe_0 - 2.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_yy[i] = 4.0 * tr_z_y[i] - 4.0 * tr_z_yyy[i] * tke_0 - 4.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_y[i] * tbe_0 + 4.0 * tr_zzz_yyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_yz[i] = 2.0 * tr_z_z[i] - 4.0 * tr_z_yyz[i] * tke_0 + tr_zz_0[i] - 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_zz_yy[i] * tke_0 + 4.0 * tr_zz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_z[i] * tbe_0 + 4.0 * tr_zzz_yyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yz[i] * tbe_0 - 2.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_zz[i] = -4.0 * tr_z_yzz[i] * tke_0 - 4.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_zz[i] * tbe_0 - 4.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_yzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-186 components of targeted buffer : DD

    auto tr_0_0_zz_xx_xx = pbuffer.data(idx_op_geom_020_dd + 180);

    auto tr_0_0_zz_xx_xy = pbuffer.data(idx_op_geom_020_dd + 181);

    auto tr_0_0_zz_xx_xz = pbuffer.data(idx_op_geom_020_dd + 182);

    auto tr_0_0_zz_xx_yy = pbuffer.data(idx_op_geom_020_dd + 183);

    auto tr_0_0_zz_xx_yz = pbuffer.data(idx_op_geom_020_dd + 184);

    auto tr_0_0_zz_xx_zz = pbuffer.data(idx_op_geom_020_dd + 185);

    #pragma omp simd aligned(tr_0_0_zz_xx_xx, tr_0_0_zz_xx_xy, tr_0_0_zz_xx_xz, tr_0_0_zz_xx_yy, tr_0_0_zz_xx_yz, tr_0_0_zz_xx_zz, tr_xx_0, tr_xx_xx, tr_xx_xxzz, tr_xx_xy, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xx_zzzz, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xx_xx[i] = -2.0 * tr_xx_xx[i] * tbe_0 - 2.0 * tr_xx_xx[i] * tke_0 + 4.0 * tr_xx_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xy[i] = -2.0 * tr_xx_xy[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xx_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xz[i] = -2.0 * tr_xx_xz[i] * tbe_0 - 6.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xx_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_x[i] * tbe_0 + 8.0 * tr_xxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_yy[i] = -2.0 * tr_xx_yy[i] * tbe_0 - 2.0 * tr_xx_yy[i] * tke_0 + 4.0 * tr_xx_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_yz[i] = -2.0 * tr_xx_yz[i] * tbe_0 - 6.0 * tr_xx_yz[i] * tke_0 + 4.0 * tr_xx_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 8.0 * tr_xxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_zz[i] = 2.0 * tr_xx_0[i] - 2.0 * tr_xx_zz[i] * tbe_0 - 10.0 * tr_xx_zz[i] * tke_0 + 4.0 * tr_xx_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxz_z[i] * tbe_0 + 8.0 * tr_xxz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 186-192 components of targeted buffer : DD

    auto tr_0_0_zz_xy_xx = pbuffer.data(idx_op_geom_020_dd + 186);

    auto tr_0_0_zz_xy_xy = pbuffer.data(idx_op_geom_020_dd + 187);

    auto tr_0_0_zz_xy_xz = pbuffer.data(idx_op_geom_020_dd + 188);

    auto tr_0_0_zz_xy_yy = pbuffer.data(idx_op_geom_020_dd + 189);

    auto tr_0_0_zz_xy_yz = pbuffer.data(idx_op_geom_020_dd + 190);

    auto tr_0_0_zz_xy_zz = pbuffer.data(idx_op_geom_020_dd + 191);

    #pragma omp simd aligned(tr_0_0_zz_xy_xx, tr_0_0_zz_xy_xy, tr_0_0_zz_xy_xz, tr_0_0_zz_xy_yy, tr_0_0_zz_xy_yz, tr_0_0_zz_xy_zz, tr_xy_0, tr_xy_xx, tr_xy_xxzz, tr_xy_xy, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xy_xx[i] = -2.0 * tr_xy_xx[i] * tbe_0 - 2.0 * tr_xy_xx[i] * tke_0 + 4.0 * tr_xy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xy[i] = -2.0 * tr_xy_xy[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xy_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xz[i] = -2.0 * tr_xy_xz[i] * tbe_0 - 6.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xy_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 8.0 * tr_xyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_yy[i] = -2.0 * tr_xy_yy[i] * tbe_0 - 2.0 * tr_xy_yy[i] * tke_0 + 4.0 * tr_xy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_yz[i] = -2.0 * tr_xy_yz[i] * tbe_0 - 6.0 * tr_xy_yz[i] * tke_0 + 4.0 * tr_xy_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 8.0 * tr_xyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_zz[i] = 2.0 * tr_xy_0[i] - 2.0 * tr_xy_zz[i] * tbe_0 - 10.0 * tr_xy_zz[i] * tke_0 + 4.0 * tr_xy_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyz_z[i] * tbe_0 + 8.0 * tr_xyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 192-198 components of targeted buffer : DD

    auto tr_0_0_zz_xz_xx = pbuffer.data(idx_op_geom_020_dd + 192);

    auto tr_0_0_zz_xz_xy = pbuffer.data(idx_op_geom_020_dd + 193);

    auto tr_0_0_zz_xz_xz = pbuffer.data(idx_op_geom_020_dd + 194);

    auto tr_0_0_zz_xz_yy = pbuffer.data(idx_op_geom_020_dd + 195);

    auto tr_0_0_zz_xz_yz = pbuffer.data(idx_op_geom_020_dd + 196);

    auto tr_0_0_zz_xz_zz = pbuffer.data(idx_op_geom_020_dd + 197);

    #pragma omp simd aligned(tr_0_0_zz_xz_xx, tr_0_0_zz_xz_xy, tr_0_0_zz_xz_xz, tr_0_0_zz_xz_yy, tr_0_0_zz_xz_yz, tr_0_0_zz_xz_zz, tr_x_x, tr_x_xxz, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xz_0, tr_xz_xx, tr_xz_xxzz, tr_xz_xy, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xz_xx[i] = -4.0 * tr_x_xxz[i] * tke_0 - 6.0 * tr_xz_xx[i] * tbe_0 - 2.0 * tr_xz_xx[i] * tke_0 + 4.0 * tr_xz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xy[i] = -4.0 * tr_x_xyz[i] * tke_0 - 6.0 * tr_xz_xy[i] * tbe_0 - 2.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xz[i] = 2.0 * tr_x_x[i] - 4.0 * tr_x_xzz[i] * tke_0 - 6.0 * tr_xz_xz[i] * tbe_0 - 6.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_x[i] * tbe_0 + 8.0 * tr_xzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_yy[i] = -4.0 * tr_x_yyz[i] * tke_0 - 6.0 * tr_xz_yy[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tke_0 + 4.0 * tr_xz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_yz[i] = 2.0 * tr_x_y[i] - 4.0 * tr_x_yzz[i] * tke_0 - 6.0 * tr_xz_yz[i] * tbe_0 - 6.0 * tr_xz_yz[i] * tke_0 + 4.0 * tr_xz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_y[i] * tbe_0 + 8.0 * tr_xzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_zz[i] = 4.0 * tr_x_z[i] - 4.0 * tr_x_zzz[i] * tke_0 + 2.0 * tr_xz_0[i] - 6.0 * tr_xz_zz[i] * tbe_0 - 10.0 * tr_xz_zz[i] * tke_0 + 4.0 * tr_xz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzz_z[i] * tbe_0 + 8.0 * tr_xzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 198-204 components of targeted buffer : DD

    auto tr_0_0_zz_yy_xx = pbuffer.data(idx_op_geom_020_dd + 198);

    auto tr_0_0_zz_yy_xy = pbuffer.data(idx_op_geom_020_dd + 199);

    auto tr_0_0_zz_yy_xz = pbuffer.data(idx_op_geom_020_dd + 200);

    auto tr_0_0_zz_yy_yy = pbuffer.data(idx_op_geom_020_dd + 201);

    auto tr_0_0_zz_yy_yz = pbuffer.data(idx_op_geom_020_dd + 202);

    auto tr_0_0_zz_yy_zz = pbuffer.data(idx_op_geom_020_dd + 203);

    #pragma omp simd aligned(tr_0_0_zz_yy_xx, tr_0_0_zz_yy_xy, tr_0_0_zz_yy_xz, tr_0_0_zz_yy_yy, tr_0_0_zz_yy_yz, tr_0_0_zz_yy_zz, tr_yy_0, tr_yy_xx, tr_yy_xxzz, tr_yy_xy, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yy_zzzz, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yy_xx[i] = -2.0 * tr_yy_xx[i] * tbe_0 - 2.0 * tr_yy_xx[i] * tke_0 + 4.0 * tr_yy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xy[i] = -2.0 * tr_yy_xy[i] * tbe_0 - 2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yy_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xz[i] = -2.0 * tr_yy_xz[i] * tbe_0 - 6.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yy_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 8.0 * tr_yyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_yy[i] = -2.0 * tr_yy_yy[i] * tbe_0 - 2.0 * tr_yy_yy[i] * tke_0 + 4.0 * tr_yy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_yz[i] = -2.0 * tr_yy_yz[i] * tbe_0 - 6.0 * tr_yy_yz[i] * tke_0 + 4.0 * tr_yy_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_y[i] * tbe_0 + 8.0 * tr_yyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_zz[i] = 2.0 * tr_yy_0[i] - 2.0 * tr_yy_zz[i] * tbe_0 - 10.0 * tr_yy_zz[i] * tke_0 + 4.0 * tr_yy_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyz_z[i] * tbe_0 + 8.0 * tr_yyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 204-210 components of targeted buffer : DD

    auto tr_0_0_zz_yz_xx = pbuffer.data(idx_op_geom_020_dd + 204);

    auto tr_0_0_zz_yz_xy = pbuffer.data(idx_op_geom_020_dd + 205);

    auto tr_0_0_zz_yz_xz = pbuffer.data(idx_op_geom_020_dd + 206);

    auto tr_0_0_zz_yz_yy = pbuffer.data(idx_op_geom_020_dd + 207);

    auto tr_0_0_zz_yz_yz = pbuffer.data(idx_op_geom_020_dd + 208);

    auto tr_0_0_zz_yz_zz = pbuffer.data(idx_op_geom_020_dd + 209);

    #pragma omp simd aligned(tr_0_0_zz_yz_xx, tr_0_0_zz_yz_xy, tr_0_0_zz_yz_xz, tr_0_0_zz_yz_yy, tr_0_0_zz_yz_yz, tr_0_0_zz_yz_zz, tr_y_x, tr_y_xxz, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yz_0, tr_yz_xx, tr_yz_xxzz, tr_yz_xy, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yz_xx[i] = -4.0 * tr_y_xxz[i] * tke_0 - 6.0 * tr_yz_xx[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tke_0 + 4.0 * tr_yz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xy[i] = -4.0 * tr_y_xyz[i] * tke_0 - 6.0 * tr_yz_xy[i] * tbe_0 - 2.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xz[i] = 2.0 * tr_y_x[i] - 4.0 * tr_y_xzz[i] * tke_0 - 6.0 * tr_yz_xz[i] * tbe_0 - 6.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_x[i] * tbe_0 + 8.0 * tr_yzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_yy[i] = -4.0 * tr_y_yyz[i] * tke_0 - 6.0 * tr_yz_yy[i] * tbe_0 - 2.0 * tr_yz_yy[i] * tke_0 + 4.0 * tr_yz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_yz[i] = 2.0 * tr_y_y[i] - 4.0 * tr_y_yzz[i] * tke_0 - 6.0 * tr_yz_yz[i] * tbe_0 - 6.0 * tr_yz_yz[i] * tke_0 + 4.0 * tr_yz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_y[i] * tbe_0 + 8.0 * tr_yzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_zz[i] = 4.0 * tr_y_z[i] - 4.0 * tr_y_zzz[i] * tke_0 + 2.0 * tr_yz_0[i] - 6.0 * tr_yz_zz[i] * tbe_0 - 10.0 * tr_yz_zz[i] * tke_0 + 4.0 * tr_yz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzz_z[i] * tbe_0 + 8.0 * tr_yzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-216 components of targeted buffer : DD

    auto tr_0_0_zz_zz_xx = pbuffer.data(idx_op_geom_020_dd + 210);

    auto tr_0_0_zz_zz_xy = pbuffer.data(idx_op_geom_020_dd + 211);

    auto tr_0_0_zz_zz_xz = pbuffer.data(idx_op_geom_020_dd + 212);

    auto tr_0_0_zz_zz_yy = pbuffer.data(idx_op_geom_020_dd + 213);

    auto tr_0_0_zz_zz_yz = pbuffer.data(idx_op_geom_020_dd + 214);

    auto tr_0_0_zz_zz_zz = pbuffer.data(idx_op_geom_020_dd + 215);

    #pragma omp simd aligned(tr_0_0_zz_zz_xx, tr_0_0_zz_zz_xy, tr_0_0_zz_zz_xz, tr_0_0_zz_zz_yy, tr_0_0_zz_zz_yz, tr_0_0_zz_zz_zz, tr_0_xx, tr_0_xy, tr_0_xz, tr_0_yy, tr_0_yz, tr_0_zz, tr_z_x, tr_z_xxz, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, tr_zz_0, tr_zz_xx, tr_zz_xxzz, tr_zz_xy, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zz_zzzz, tr_zzz_x, tr_zzz_xxz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_zz_xx[i] = 2.0 * tr_0_xx[i] - 8.0 * tr_z_xxz[i] * tke_0 - 10.0 * tr_zz_xx[i] * tbe_0 - 2.0 * tr_zz_xx[i] * tke_0 + 4.0 * tr_zz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xy[i] = 2.0 * tr_0_xy[i] - 8.0 * tr_z_xyz[i] * tke_0 - 10.0 * tr_zz_xy[i] * tbe_0 - 2.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xz[i] = 2.0 * tr_0_xz[i] + 4.0 * tr_z_x[i] - 8.0 * tr_z_xzz[i] * tke_0 - 10.0 * tr_zz_xz[i] * tbe_0 - 6.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_zz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_x[i] * tbe_0 + 8.0 * tr_zzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_yy[i] = 2.0 * tr_0_yy[i] - 8.0 * tr_z_yyz[i] * tke_0 - 10.0 * tr_zz_yy[i] * tbe_0 - 2.0 * tr_zz_yy[i] * tke_0 + 4.0 * tr_zz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_yz[i] = 2.0 * tr_0_yz[i] + 4.0 * tr_z_y[i] - 8.0 * tr_z_yzz[i] * tke_0 - 10.0 * tr_zz_yz[i] * tbe_0 - 6.0 * tr_zz_yz[i] * tke_0 + 4.0 * tr_zz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_y[i] * tbe_0 + 8.0 * tr_zzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_zz[i] = 2.0 * tr_0_zz[i] + 8.0 * tr_z_z[i] - 8.0 * tr_z_zzz[i] * tke_0 + 2.0 * tr_zz_0[i] - 10.0 * tr_zz_zz[i] * tbe_0 - 10.0 * tr_zz_zz[i] * tke_0 + 4.0 * tr_zz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzz_z[i] * tbe_0 + 8.0 * tr_zzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_zz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

