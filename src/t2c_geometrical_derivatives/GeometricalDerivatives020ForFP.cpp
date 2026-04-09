#include "GeometricalDerivatives020ForFP.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_fp(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_fp,
                         const int idx_op_pp,
                         const int idx_op_ds,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_hp,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up components of auxiliary buffer : GS

    auto tr_xxxx_0 = pbuffer.data(idx_op_gs);

    auto tr_xxxy_0 = pbuffer.data(idx_op_gs + 1);

    auto tr_xxxz_0 = pbuffer.data(idx_op_gs + 2);

    auto tr_xxyy_0 = pbuffer.data(idx_op_gs + 3);

    auto tr_xxyz_0 = pbuffer.data(idx_op_gs + 4);

    auto tr_xxzz_0 = pbuffer.data(idx_op_gs + 5);

    auto tr_xyyy_0 = pbuffer.data(idx_op_gs + 6);

    auto tr_xyyz_0 = pbuffer.data(idx_op_gs + 7);

    auto tr_xyzz_0 = pbuffer.data(idx_op_gs + 8);

    auto tr_xzzz_0 = pbuffer.data(idx_op_gs + 9);

    auto tr_yyyy_0 = pbuffer.data(idx_op_gs + 10);

    auto tr_yyyz_0 = pbuffer.data(idx_op_gs + 11);

    auto tr_yyzz_0 = pbuffer.data(idx_op_gs + 12);

    auto tr_yzzz_0 = pbuffer.data(idx_op_gs + 13);

    auto tr_zzzz_0 = pbuffer.data(idx_op_gs + 14);

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

    // Set up components of auxiliary buffer : HP

    auto tr_xxxxx_x = pbuffer.data(idx_op_hp);

    auto tr_xxxxx_y = pbuffer.data(idx_op_hp + 1);

    auto tr_xxxxx_z = pbuffer.data(idx_op_hp + 2);

    auto tr_xxxxy_x = pbuffer.data(idx_op_hp + 3);

    auto tr_xxxxy_y = pbuffer.data(idx_op_hp + 4);

    auto tr_xxxxy_z = pbuffer.data(idx_op_hp + 5);

    auto tr_xxxxz_x = pbuffer.data(idx_op_hp + 6);

    auto tr_xxxxz_y = pbuffer.data(idx_op_hp + 7);

    auto tr_xxxxz_z = pbuffer.data(idx_op_hp + 8);

    auto tr_xxxyy_x = pbuffer.data(idx_op_hp + 9);

    auto tr_xxxyy_y = pbuffer.data(idx_op_hp + 10);

    auto tr_xxxyy_z = pbuffer.data(idx_op_hp + 11);

    auto tr_xxxyz_x = pbuffer.data(idx_op_hp + 12);

    auto tr_xxxyz_y = pbuffer.data(idx_op_hp + 13);

    auto tr_xxxyz_z = pbuffer.data(idx_op_hp + 14);

    auto tr_xxxzz_x = pbuffer.data(idx_op_hp + 15);

    auto tr_xxxzz_y = pbuffer.data(idx_op_hp + 16);

    auto tr_xxxzz_z = pbuffer.data(idx_op_hp + 17);

    auto tr_xxyyy_x = pbuffer.data(idx_op_hp + 18);

    auto tr_xxyyy_y = pbuffer.data(idx_op_hp + 19);

    auto tr_xxyyy_z = pbuffer.data(idx_op_hp + 20);

    auto tr_xxyyz_x = pbuffer.data(idx_op_hp + 21);

    auto tr_xxyyz_y = pbuffer.data(idx_op_hp + 22);

    auto tr_xxyyz_z = pbuffer.data(idx_op_hp + 23);

    auto tr_xxyzz_x = pbuffer.data(idx_op_hp + 24);

    auto tr_xxyzz_y = pbuffer.data(idx_op_hp + 25);

    auto tr_xxyzz_z = pbuffer.data(idx_op_hp + 26);

    auto tr_xxzzz_x = pbuffer.data(idx_op_hp + 27);

    auto tr_xxzzz_y = pbuffer.data(idx_op_hp + 28);

    auto tr_xxzzz_z = pbuffer.data(idx_op_hp + 29);

    auto tr_xyyyy_x = pbuffer.data(idx_op_hp + 30);

    auto tr_xyyyy_y = pbuffer.data(idx_op_hp + 31);

    auto tr_xyyyy_z = pbuffer.data(idx_op_hp + 32);

    auto tr_xyyyz_x = pbuffer.data(idx_op_hp + 33);

    auto tr_xyyyz_y = pbuffer.data(idx_op_hp + 34);

    auto tr_xyyyz_z = pbuffer.data(idx_op_hp + 35);

    auto tr_xyyzz_x = pbuffer.data(idx_op_hp + 36);

    auto tr_xyyzz_y = pbuffer.data(idx_op_hp + 37);

    auto tr_xyyzz_z = pbuffer.data(idx_op_hp + 38);

    auto tr_xyzzz_x = pbuffer.data(idx_op_hp + 39);

    auto tr_xyzzz_y = pbuffer.data(idx_op_hp + 40);

    auto tr_xyzzz_z = pbuffer.data(idx_op_hp + 41);

    auto tr_xzzzz_x = pbuffer.data(idx_op_hp + 42);

    auto tr_xzzzz_y = pbuffer.data(idx_op_hp + 43);

    auto tr_xzzzz_z = pbuffer.data(idx_op_hp + 44);

    auto tr_yyyyy_x = pbuffer.data(idx_op_hp + 45);

    auto tr_yyyyy_y = pbuffer.data(idx_op_hp + 46);

    auto tr_yyyyy_z = pbuffer.data(idx_op_hp + 47);

    auto tr_yyyyz_x = pbuffer.data(idx_op_hp + 48);

    auto tr_yyyyz_y = pbuffer.data(idx_op_hp + 49);

    auto tr_yyyyz_z = pbuffer.data(idx_op_hp + 50);

    auto tr_yyyzz_x = pbuffer.data(idx_op_hp + 51);

    auto tr_yyyzz_y = pbuffer.data(idx_op_hp + 52);

    auto tr_yyyzz_z = pbuffer.data(idx_op_hp + 53);

    auto tr_yyzzz_x = pbuffer.data(idx_op_hp + 54);

    auto tr_yyzzz_y = pbuffer.data(idx_op_hp + 55);

    auto tr_yyzzz_z = pbuffer.data(idx_op_hp + 56);

    auto tr_yzzzz_x = pbuffer.data(idx_op_hp + 57);

    auto tr_yzzzz_y = pbuffer.data(idx_op_hp + 58);

    auto tr_yzzzz_z = pbuffer.data(idx_op_hp + 59);

    auto tr_zzzzz_x = pbuffer.data(idx_op_hp + 60);

    auto tr_zzzzz_y = pbuffer.data(idx_op_hp + 61);

    auto tr_zzzzz_z = pbuffer.data(idx_op_hp + 62);

    // Set up 0-3 components of targeted buffer : FP

    auto tr_0_0_xx_xxx_x = pbuffer.data(idx_op_geom_020_fp);

    auto tr_0_0_xx_xxx_y = pbuffer.data(idx_op_geom_020_fp + 1);

    auto tr_0_0_xx_xxx_z = pbuffer.data(idx_op_geom_020_fp + 2);

    #pragma omp simd aligned(tr_0_0_xx_xxx_x, tr_0_0_xx_xxx_y, tr_0_0_xx_xxx_z, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxxx_x, tr_xxxxx_y, tr_xxxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxx_x[i] = 6.0 * tr_x_x[i] + 6.0 * tr_xx_0[i] - 12.0 * tr_xx_xx[i] * tke_0 - 14.0 * tr_xxx_x[i] * tbe_0 - 6.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxx_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xxxx_0[i] * tbe_0 + 8.0 * tr_xxxx_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_y[i] = 6.0 * tr_x_y[i] - 12.0 * tr_xx_xy[i] * tke_0 - 14.0 * tr_xxx_y[i] * tbe_0 - 2.0 * tr_xxx_y[i] * tke_0 + 4.0 * tr_xxx_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_z[i] = 6.0 * tr_x_z[i] - 12.0 * tr_xx_xz[i] * tke_0 - 14.0 * tr_xxx_z[i] * tbe_0 - 2.0 * tr_xxx_z[i] * tke_0 + 4.0 * tr_xxx_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_z[i] * tbe_0 * tbe_0;
    }

    // Set up 3-6 components of targeted buffer : FP

    auto tr_0_0_xx_xxy_x = pbuffer.data(idx_op_geom_020_fp + 3);

    auto tr_0_0_xx_xxy_y = pbuffer.data(idx_op_geom_020_fp + 4);

    auto tr_0_0_xx_xxy_z = pbuffer.data(idx_op_geom_020_fp + 5);

    #pragma omp simd aligned(tr_0_0_xx_xxy_x, tr_0_0_xx_xxy_y, tr_0_0_xx_xxy_z, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_y, tr_xxy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_y_x, tr_y_y, tr_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxy_x[i] = 2.0 * tr_y_x[i] + 4.0 * tr_xy_0[i] - 8.0 * tr_xy_xx[i] * tke_0 - 10.0 * tr_xxy_x[i] * tbe_0 - 6.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxy_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_0[i] * tbe_0 + 8.0 * tr_xxxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_y[i] = 2.0 * tr_y_y[i] - 8.0 * tr_xy_xy[i] * tke_0 - 10.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxy_y[i] * tke_0 + 4.0 * tr_xxy_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_z[i] = 2.0 * tr_y_z[i] - 8.0 * tr_xy_xz[i] * tke_0 - 10.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxy_z[i] * tke_0 + 4.0 * tr_xxy_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 6-9 components of targeted buffer : FP

    auto tr_0_0_xx_xxz_x = pbuffer.data(idx_op_geom_020_fp + 6);

    auto tr_0_0_xx_xxz_y = pbuffer.data(idx_op_geom_020_fp + 7);

    auto tr_0_0_xx_xxz_z = pbuffer.data(idx_op_geom_020_fp + 8);

    #pragma omp simd aligned(tr_0_0_xx_xxz_x, tr_0_0_xx_xxz_y, tr_0_0_xx_xxz_z, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_y, tr_xxz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxz_x[i] = 2.0 * tr_z_x[i] + 4.0 * tr_xz_0[i] - 8.0 * tr_xz_xx[i] * tke_0 - 10.0 * tr_xxz_x[i] * tbe_0 - 6.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxz_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_0[i] * tbe_0 + 8.0 * tr_xxxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_y[i] = 2.0 * tr_z_y[i] - 8.0 * tr_xz_xy[i] * tke_0 - 10.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxz_y[i] * tke_0 + 4.0 * tr_xxz_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_z[i] = 2.0 * tr_z_z[i] - 8.0 * tr_xz_xz[i] * tke_0 - 10.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxz_z[i] * tke_0 + 4.0 * tr_xxz_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 9-12 components of targeted buffer : FP

    auto tr_0_0_xx_xyy_x = pbuffer.data(idx_op_geom_020_fp + 9);

    auto tr_0_0_xx_xyy_y = pbuffer.data(idx_op_geom_020_fp + 10);

    auto tr_0_0_xx_xyy_z = pbuffer.data(idx_op_geom_020_fp + 11);

    #pragma omp simd aligned(tr_0_0_xx_xyy_x, tr_0_0_xx_xyy_y, tr_0_0_xx_xyy_z, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_y, tr_xyy_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyy_x[i] = 2.0 * tr_yy_0[i] - 4.0 * tr_yy_xx[i] * tke_0 - 6.0 * tr_xyy_x[i] * tbe_0 - 6.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyy_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_0[i] * tbe_0 + 8.0 * tr_xxyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_y[i] = -4.0 * tr_yy_xy[i] * tke_0 - 6.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xyy_y[i] * tke_0 + 4.0 * tr_xyy_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_z[i] = -4.0 * tr_yy_xz[i] * tke_0 - 6.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyy_z[i] * tke_0 + 4.0 * tr_xyy_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 12-15 components of targeted buffer : FP

    auto tr_0_0_xx_xyz_x = pbuffer.data(idx_op_geom_020_fp + 12);

    auto tr_0_0_xx_xyz_y = pbuffer.data(idx_op_geom_020_fp + 13);

    auto tr_0_0_xx_xyz_z = pbuffer.data(idx_op_geom_020_fp + 14);

    #pragma omp simd aligned(tr_0_0_xx_xyz_x, tr_0_0_xx_xyz_y, tr_0_0_xx_xyz_z, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_y, tr_xyz_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyz_x[i] = 2.0 * tr_yz_0[i] - 4.0 * tr_yz_xx[i] * tke_0 - 6.0 * tr_xyz_x[i] * tbe_0 - 6.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyz_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 8.0 * tr_xxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_y[i] = -4.0 * tr_yz_xy[i] * tke_0 - 6.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyz_y[i] * tke_0 + 4.0 * tr_xyz_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_z[i] = -4.0 * tr_yz_xz[i] * tke_0 - 6.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyz_z[i] * tke_0 + 4.0 * tr_xyz_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 15-18 components of targeted buffer : FP

    auto tr_0_0_xx_xzz_x = pbuffer.data(idx_op_geom_020_fp + 15);

    auto tr_0_0_xx_xzz_y = pbuffer.data(idx_op_geom_020_fp + 16);

    auto tr_0_0_xx_xzz_z = pbuffer.data(idx_op_geom_020_fp + 17);

    #pragma omp simd aligned(tr_0_0_xx_xzz_x, tr_0_0_xx_xzz_y, tr_0_0_xx_xzz_z, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_y, tr_xzz_z, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xzz_x[i] = 2.0 * tr_zz_0[i] - 4.0 * tr_zz_xx[i] * tke_0 - 6.0 * tr_xzz_x[i] * tbe_0 - 6.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xzz_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_0[i] * tbe_0 + 8.0 * tr_xxzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_y[i] = -4.0 * tr_zz_xy[i] * tke_0 - 6.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xzz_y[i] * tke_0 + 4.0 * tr_xzz_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_z[i] = -4.0 * tr_zz_xz[i] * tke_0 - 6.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xzz_z[i] * tke_0 + 4.0 * tr_xzz_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 18-21 components of targeted buffer : FP

    auto tr_0_0_xx_yyy_x = pbuffer.data(idx_op_geom_020_fp + 18);

    auto tr_0_0_xx_yyy_y = pbuffer.data(idx_op_geom_020_fp + 19);

    auto tr_0_0_xx_yyy_z = pbuffer.data(idx_op_geom_020_fp + 20);

    #pragma omp simd aligned(tr_0_0_xx_yyy_x, tr_0_0_xx_yyy_y, tr_0_0_xx_yyy_z, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyy_x[i] = -2.0 * tr_yyy_x[i] * tbe_0 - 6.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyy_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_0[i] * tbe_0 + 8.0 * tr_xyyy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_y[i] = -2.0 * tr_yyy_y[i] * tbe_0 - 2.0 * tr_yyy_y[i] * tke_0 + 4.0 * tr_yyy_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_z[i] = -2.0 * tr_yyy_z[i] * tbe_0 - 2.0 * tr_yyy_z[i] * tke_0 + 4.0 * tr_yyy_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 21-24 components of targeted buffer : FP

    auto tr_0_0_xx_yyz_x = pbuffer.data(idx_op_geom_020_fp + 21);

    auto tr_0_0_xx_yyz_y = pbuffer.data(idx_op_geom_020_fp + 22);

    auto tr_0_0_xx_yyz_z = pbuffer.data(idx_op_geom_020_fp + 23);

    #pragma omp simd aligned(tr_0_0_xx_yyz_x, tr_0_0_xx_yyz_y, tr_0_0_xx_yyz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyz_x[i] = -2.0 * tr_yyz_x[i] * tbe_0 - 6.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyz_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 + 8.0 * tr_xyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_y[i] = -2.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_yyz_y[i] * tke_0 + 4.0 * tr_yyz_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_z[i] = -2.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyz_z[i] * tke_0 + 4.0 * tr_yyz_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 24-27 components of targeted buffer : FP

    auto tr_0_0_xx_yzz_x = pbuffer.data(idx_op_geom_020_fp + 24);

    auto tr_0_0_xx_yzz_y = pbuffer.data(idx_op_geom_020_fp + 25);

    auto tr_0_0_xx_yzz_z = pbuffer.data(idx_op_geom_020_fp + 26);

    #pragma omp simd aligned(tr_0_0_xx_yzz_x, tr_0_0_xx_yzz_y, tr_0_0_xx_yzz_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yzz_x[i] = -2.0 * tr_yzz_x[i] * tbe_0 - 6.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yzz_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 + 8.0 * tr_xyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_y[i] = -2.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_yzz_y[i] * tke_0 + 4.0 * tr_yzz_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_z[i] = -2.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yzz_z[i] * tke_0 + 4.0 * tr_yzz_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 27-30 components of targeted buffer : FP

    auto tr_0_0_xx_zzz_x = pbuffer.data(idx_op_geom_020_fp + 27);

    auto tr_0_0_xx_zzz_y = pbuffer.data(idx_op_geom_020_fp + 28);

    auto tr_0_0_xx_zzz_z = pbuffer.data(idx_op_geom_020_fp + 29);

    #pragma omp simd aligned(tr_0_0_xx_zzz_x, tr_0_0_xx_zzz_y, tr_0_0_xx_zzz_z, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_zzz_x[i] = -2.0 * tr_zzz_x[i] * tbe_0 - 6.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_zzz_xxx[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_0[i] * tbe_0 + 8.0 * tr_xzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_y[i] = -2.0 * tr_zzz_y[i] * tbe_0 - 2.0 * tr_zzz_y[i] * tke_0 + 4.0 * tr_zzz_xxy[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_z[i] = -2.0 * tr_zzz_z[i] * tbe_0 - 2.0 * tr_zzz_z[i] * tke_0 + 4.0 * tr_zzz_xxz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 30-33 components of targeted buffer : FP

    auto tr_0_0_xy_xxx_x = pbuffer.data(idx_op_geom_020_fp + 30);

    auto tr_0_0_xy_xxx_y = pbuffer.data(idx_op_geom_020_fp + 31);

    auto tr_0_0_xy_xxx_z = pbuffer.data(idx_op_geom_020_fp + 32);

    #pragma omp simd aligned(tr_0_0_xy_xxx_x, tr_0_0_xy_xxx_y, tr_0_0_xy_xxx_z, tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxx_x, tr_xxx_xxy, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_y, tr_xxxx_0, tr_xxxx_xy, tr_xxxx_yy, tr_xxxx_yz, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxy_x, tr_xxy_y, tr_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxx_x[i] = -6.0 * tr_xx_xy[i] * tke_0 - 6.0 * tr_xxy_x[i] * tbe_0 - 2.0 * tr_xxx_y[i] * tke_0 + 4.0 * tr_xxx_xxy[i] * tke_0 * tke_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_y[i] = 3.0 * tr_xx_0[i] - 6.0 * tr_xx_yy[i] * tke_0 - 6.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxx_xyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxx_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_z[i] = -6.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxx_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 33-36 components of targeted buffer : FP

    auto tr_0_0_xy_xxy_x = pbuffer.data(idx_op_geom_020_fp + 33);

    auto tr_0_0_xy_xxy_y = pbuffer.data(idx_op_geom_020_fp + 34);

    auto tr_0_0_xy_xxy_z = pbuffer.data(idx_op_geom_020_fp + 35);

    #pragma omp simd aligned(tr_0_0_xy_xxy_x, tr_0_0_xy_xxy_y, tr_0_0_xy_xxy_z, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxy_0, tr_xxxy_xy, tr_xxxy_yy, tr_xxxy_yz, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxy_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xy_xy[i] * tke_0 - 4.0 * tr_xyy_x[i] * tbe_0 + tr_xx_0[i] - 2.0 * tr_xx_xx[i] * tke_0 - 2.0 * tr_xxy_y[i] * tke_0 + 4.0 * tr_xxy_xxy[i] * tke_0 * tke_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_y[i] = 2.0 * tr_x_y[i] + 2.0 * tr_xy_0[i] - 4.0 * tr_xy_yy[i] * tke_0 - 4.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tke_0 - 2.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxy_xyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_y[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_z[i] = 2.0 * tr_x_z[i] - 4.0 * tr_xy_yz[i] * tke_0 - 4.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tke_0 + 4.0 * tr_xxy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_z[i] * tbe_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 36-39 components of targeted buffer : FP

    auto tr_0_0_xy_xxz_x = pbuffer.data(idx_op_geom_020_fp + 36);

    auto tr_0_0_xy_xxz_y = pbuffer.data(idx_op_geom_020_fp + 37);

    auto tr_0_0_xy_xxz_z = pbuffer.data(idx_op_geom_020_fp + 38);

    #pragma omp simd aligned(tr_0_0_xy_xxz_x, tr_0_0_xy_xxz_y, tr_0_0_xy_xxz_z, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xy, tr_xxxz_yy, tr_xxxz_yz, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxz_x[i] = -4.0 * tr_xz_xy[i] * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xxz_y[i] * tke_0 + 4.0 * tr_xxz_xxy[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_y[i] = 2.0 * tr_xz_0[i] - 4.0 * tr_xz_yy[i] * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxz_xyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_z[i] = -4.0 * tr_xz_yz[i] * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xxz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 39-42 components of targeted buffer : FP

    auto tr_0_0_xy_xyy_x = pbuffer.data(idx_op_geom_020_fp + 39);

    auto tr_0_0_xy_xyy_y = pbuffer.data(idx_op_geom_020_fp + 40);

    auto tr_0_0_xy_xyy_z = pbuffer.data(idx_op_geom_020_fp + 41);

    #pragma omp simd aligned(tr_0_0_xy_xyy_x, tr_0_0_xy_xyy_y, tr_0_0_xy_xyy_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyy_xy, tr_xxyy_yy, tr_xxyy_yz, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyy_x[i] = 2.0 * tr_y_x[i] - 2.0 * tr_yy_xy[i] * tke_0 - 2.0 * tr_yyy_x[i] * tbe_0 + 2.0 * tr_xy_0[i] - 4.0 * tr_xy_xx[i] * tke_0 - 2.0 * tr_xyy_y[i] * tke_0 + 4.0 * tr_xyy_xxy[i] * tke_0 * tke_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_y[i] = 2.0 * tr_y_y[i] + tr_yy_0[i] - 2.0 * tr_yy_yy[i] * tke_0 - 2.0 * tr_yyy_y[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tke_0 - 2.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyy_xyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_z[i] = 2.0 * tr_y_z[i] - 2.0 * tr_yy_yz[i] * tke_0 - 2.0 * tr_yyy_z[i] * tbe_0 - 4.0 * tr_xy_xz[i] * tke_0 + 4.0 * tr_xyy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_z[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 42-45 components of targeted buffer : FP

    auto tr_0_0_xy_xyz_x = pbuffer.data(idx_op_geom_020_fp + 42);

    auto tr_0_0_xy_xyz_y = pbuffer.data(idx_op_geom_020_fp + 43);

    auto tr_0_0_xy_xyz_z = pbuffer.data(idx_op_geom_020_fp + 44);

    #pragma omp simd aligned(tr_0_0_xy_xyz_x, tr_0_0_xy_xyz_y, tr_0_0_xy_xyz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyz_x[i] = tr_z_x[i] - 2.0 * tr_yz_xy[i] * tke_0 - 2.0 * tr_yyz_x[i] * tbe_0 + tr_xz_0[i] - 2.0 * tr_xz_xx[i] * tke_0 - 2.0 * tr_xyz_y[i] * tke_0 + 4.0 * tr_xyz_xxy[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_y[i] = tr_z_y[i] + tr_yz_0[i] - 2.0 * tr_yz_yy[i] * tke_0 - 2.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_xz_xy[i] * tke_0 - 2.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyz_xyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_z[i] = tr_z_z[i] - 2.0 * tr_yz_yz[i] * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_xz_xz[i] * tke_0 + 4.0 * tr_xyz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 45-48 components of targeted buffer : FP

    auto tr_0_0_xy_xzz_x = pbuffer.data(idx_op_geom_020_fp + 45);

    auto tr_0_0_xy_xzz_y = pbuffer.data(idx_op_geom_020_fp + 46);

    auto tr_0_0_xy_xzz_z = pbuffer.data(idx_op_geom_020_fp + 47);

    #pragma omp simd aligned(tr_0_0_xy_xzz_x, tr_0_0_xy_xzz_y, tr_0_0_xy_xzz_z, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xxzz_0, tr_xxzz_xy, tr_xxzz_yy, tr_xxzz_yz, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xzz_x[i] = -2.0 * tr_zz_xy[i] * tke_0 - 2.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_xzz_y[i] * tke_0 + 4.0 * tr_xzz_xxy[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_y[i] = tr_zz_0[i] - 2.0 * tr_zz_yy[i] * tke_0 - 2.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xzz_xyy[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_z[i] = -2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 + 4.0 * tr_xzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 48-51 components of targeted buffer : FP

    auto tr_0_0_xy_yyy_x = pbuffer.data(idx_op_geom_020_fp + 48);

    auto tr_0_0_xy_yyy_y = pbuffer.data(idx_op_geom_020_fp + 49);

    auto tr_0_0_xy_yyy_z = pbuffer.data(idx_op_geom_020_fp + 50);

    #pragma omp simd aligned(tr_0_0_xy_yyy_x, tr_0_0_xy_yyy_y, tr_0_0_xy_yyy_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyy_xy, tr_xyyy_yy, tr_xyyy_yz, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yyy_x, tr_yyy_xxy, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_y, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyy_x[i] = 3.0 * tr_yy_0[i] - 6.0 * tr_yy_xx[i] * tke_0 - 2.0 * tr_yyy_y[i] * tke_0 + 4.0 * tr_yyy_xxy[i] * tke_0 * tke_0 - 2.0 * tr_yyyy_0[i] * tbe_0 + 4.0 * tr_yyyy_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_y[i] = -6.0 * tr_yy_xy[i] * tke_0 - 2.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyy_xyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_z[i] = -6.0 * tr_yy_xz[i] * tke_0 + 4.0 * tr_yyy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_z[i] * tbe_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 51-54 components of targeted buffer : FP

    auto tr_0_0_xy_yyz_x = pbuffer.data(idx_op_geom_020_fp + 51);

    auto tr_0_0_xy_yyz_y = pbuffer.data(idx_op_geom_020_fp + 52);

    auto tr_0_0_xy_yyz_z = pbuffer.data(idx_op_geom_020_fp + 53);

    #pragma omp simd aligned(tr_0_0_xy_yyz_x, tr_0_0_xy_yyz_y, tr_0_0_xy_yyz_z, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyz_x[i] = 2.0 * tr_yz_0[i] - 4.0 * tr_yz_xx[i] * tke_0 - 2.0 * tr_yyz_y[i] * tke_0 + 4.0 * tr_yyz_xxy[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_y[i] = -4.0 * tr_yz_xy[i] * tke_0 - 2.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyz_xyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_z[i] = -4.0 * tr_yz_xz[i] * tke_0 + 4.0 * tr_yyz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 54-57 components of targeted buffer : FP

    auto tr_0_0_xy_yzz_x = pbuffer.data(idx_op_geom_020_fp + 54);

    auto tr_0_0_xy_yzz_y = pbuffer.data(idx_op_geom_020_fp + 55);

    auto tr_0_0_xy_yzz_z = pbuffer.data(idx_op_geom_020_fp + 56);

    #pragma omp simd aligned(tr_0_0_xy_yzz_x, tr_0_0_xy_yzz_y, tr_0_0_xy_yzz_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yzz_x[i] = tr_zz_0[i] - 2.0 * tr_zz_xx[i] * tke_0 - 2.0 * tr_yzz_y[i] * tke_0 + 4.0 * tr_yzz_xxy[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_y[i] = -2.0 * tr_zz_xy[i] * tke_0 - 2.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yzz_xyy[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_z[i] = -2.0 * tr_zz_xz[i] * tke_0 + 4.0 * tr_yzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 57-60 components of targeted buffer : FP

    auto tr_0_0_xy_zzz_x = pbuffer.data(idx_op_geom_020_fp + 57);

    auto tr_0_0_xy_zzz_y = pbuffer.data(idx_op_geom_020_fp + 58);

    auto tr_0_0_xy_zzz_z = pbuffer.data(idx_op_geom_020_fp + 59);

    #pragma omp simd aligned(tr_0_0_xy_zzz_x, tr_0_0_xy_zzz_y, tr_0_0_xy_zzz_z, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_xzzz_0, tr_xzzz_xy, tr_xzzz_yy, tr_xzzz_yz, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_zzz_x, tr_zzz_xxy, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_zzz_x[i] = -2.0 * tr_zzz_y[i] * tke_0 + 4.0 * tr_zzz_xxy[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_y[i] = -2.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_zzz_xyy[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_z[i] = 4.0 * tr_zzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 60-63 components of targeted buffer : FP

    auto tr_0_0_xz_xxx_x = pbuffer.data(idx_op_geom_020_fp + 60);

    auto tr_0_0_xz_xxx_y = pbuffer.data(idx_op_geom_020_fp + 61);

    auto tr_0_0_xz_xxx_z = pbuffer.data(idx_op_geom_020_fp + 62);

    #pragma omp simd aligned(tr_0_0_xz_xxx_x, tr_0_0_xz_xxx_y, tr_0_0_xz_xxx_z, tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_z, tr_xxxx_0, tr_xxxx_xz, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxz_x, tr_xxz_y, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxx_x[i] = -6.0 * tr_xx_xz[i] * tke_0 - 6.0 * tr_xxz_x[i] * tbe_0 - 2.0 * tr_xxx_z[i] * tke_0 + 4.0 * tr_xxx_xxz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_y[i] = -6.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxx_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_z[i] = 3.0 * tr_xx_0[i] - 6.0 * tr_xx_zz[i] * tke_0 - 6.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxx_xzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxx_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 63-66 components of targeted buffer : FP

    auto tr_0_0_xz_xxy_x = pbuffer.data(idx_op_geom_020_fp + 63);

    auto tr_0_0_xz_xxy_y = pbuffer.data(idx_op_geom_020_fp + 64);

    auto tr_0_0_xz_xxy_z = pbuffer.data(idx_op_geom_020_fp + 65);

    #pragma omp simd aligned(tr_0_0_xz_xxy_x, tr_0_0_xz_xxy_y, tr_0_0_xz_xxy_z, tr_xxxy_0, tr_xxxy_xz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxy_x[i] = -4.0 * tr_xy_xz[i] * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xxy_z[i] * tke_0 + 4.0 * tr_xxy_xxz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_y[i] = -4.0 * tr_xy_yz[i] * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xxy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_z[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xy_zz[i] * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxy_xzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 66-69 components of targeted buffer : FP

    auto tr_0_0_xz_xxz_x = pbuffer.data(idx_op_geom_020_fp + 66);

    auto tr_0_0_xz_xxz_y = pbuffer.data(idx_op_geom_020_fp + 67);

    auto tr_0_0_xz_xxz_z = pbuffer.data(idx_op_geom_020_fp + 68);

    #pragma omp simd aligned(tr_0_0_xz_xxz_x, tr_0_0_xz_xxz_y, tr_0_0_xz_xxz_z, tr_x_x, tr_x_y, tr_x_z, tr_xx_0, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxz_0, tr_xxxz_xz, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxz_x[i] = 2.0 * tr_x_x[i] - 4.0 * tr_xz_xz[i] * tke_0 - 4.0 * tr_xzz_x[i] * tbe_0 + tr_xx_0[i] - 2.0 * tr_xx_xx[i] * tke_0 - 2.0 * tr_xxz_z[i] * tke_0 + 4.0 * tr_xxz_xxz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_x[i] * tbe_0 + 4.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_y[i] = 2.0 * tr_x_y[i] - 4.0 * tr_xz_yz[i] * tke_0 - 4.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xxz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_y[i] * tbe_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_z[i] = 2.0 * tr_x_z[i] + 2.0 * tr_xz_0[i] - 4.0 * tr_xz_zz[i] * tke_0 - 4.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xx_xz[i] * tke_0 - 2.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxz_xzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_z[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 69-72 components of targeted buffer : FP

    auto tr_0_0_xz_xyy_x = pbuffer.data(idx_op_geom_020_fp + 69);

    auto tr_0_0_xz_xyy_y = pbuffer.data(idx_op_geom_020_fp + 70);

    auto tr_0_0_xz_xyy_z = pbuffer.data(idx_op_geom_020_fp + 71);

    #pragma omp simd aligned(tr_0_0_xz_xyy_x, tr_0_0_xz_xyy_y, tr_0_0_xz_xyy_z, tr_xxyy_0, tr_xxyy_xz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyy_x[i] = -2.0 * tr_yy_xz[i] * tke_0 - 2.0 * tr_yyz_x[i] * tbe_0 - 2.0 * tr_xyy_z[i] * tke_0 + 4.0 * tr_xyy_xxz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_y[i] = -2.0 * tr_yy_yz[i] * tke_0 - 2.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_xyy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_z[i] = tr_yy_0[i] - 2.0 * tr_yy_zz[i] * tke_0 - 2.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyy_xzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 72-75 components of targeted buffer : FP

    auto tr_0_0_xz_xyz_x = pbuffer.data(idx_op_geom_020_fp + 72);

    auto tr_0_0_xz_xyz_y = pbuffer.data(idx_op_geom_020_fp + 73);

    auto tr_0_0_xz_xyz_z = pbuffer.data(idx_op_geom_020_fp + 74);

    #pragma omp simd aligned(tr_0_0_xz_xyz_x, tr_0_0_xz_xyz_y, tr_0_0_xz_xyz_z, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xy_0, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyz_x[i] = tr_y_x[i] - 2.0 * tr_yz_xz[i] * tke_0 - 2.0 * tr_yzz_x[i] * tbe_0 + tr_xy_0[i] - 2.0 * tr_xy_xx[i] * tke_0 - 2.0 * tr_xyz_z[i] * tke_0 + 4.0 * tr_xyz_xxz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_y[i] = tr_y_y[i] - 2.0 * tr_yz_yz[i] * tke_0 - 2.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xyz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_z[i] = tr_y_z[i] + tr_yz_0[i] - 2.0 * tr_yz_zz[i] * tke_0 - 2.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tke_0 - 2.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyz_xzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 75-78 components of targeted buffer : FP

    auto tr_0_0_xz_xzz_x = pbuffer.data(idx_op_geom_020_fp + 75);

    auto tr_0_0_xz_xzz_y = pbuffer.data(idx_op_geom_020_fp + 76);

    auto tr_0_0_xz_xzz_z = pbuffer.data(idx_op_geom_020_fp + 77);

    #pragma omp simd aligned(tr_0_0_xz_xzz_x, tr_0_0_xz_xzz_y, tr_0_0_xz_xzz_z, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xxzz_xz, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_xz_0, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xzz_x[i] = 2.0 * tr_z_x[i] - 2.0 * tr_zz_xz[i] * tke_0 - 2.0 * tr_zzz_x[i] * tbe_0 + 2.0 * tr_xz_0[i] - 4.0 * tr_xz_xx[i] * tke_0 - 2.0 * tr_xzz_z[i] * tke_0 + 4.0 * tr_xzz_xxz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_y[i] = 2.0 * tr_z_y[i] - 2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_zzz_y[i] * tbe_0 - 4.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_y[i] * tbe_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_z[i] = 2.0 * tr_z_z[i] + tr_zz_0[i] - 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_zzz_z[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tke_0 - 2.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xzz_xzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 78-81 components of targeted buffer : FP

    auto tr_0_0_xz_yyy_x = pbuffer.data(idx_op_geom_020_fp + 78);

    auto tr_0_0_xz_yyy_y = pbuffer.data(idx_op_geom_020_fp + 79);

    auto tr_0_0_xz_yyy_z = pbuffer.data(idx_op_geom_020_fp + 80);

    #pragma omp simd aligned(tr_0_0_xz_yyy_x, tr_0_0_xz_yyy_y, tr_0_0_xz_yyy_z, tr_xyyy_0, tr_xyyy_xz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_yyy_x, tr_yyy_xxz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyy_x[i] = -2.0 * tr_yyy_z[i] * tke_0 + 4.0 * tr_yyy_xxz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_xx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_y[i] = 4.0 * tr_yyy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_z[i] = -2.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyy_xzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 81-84 components of targeted buffer : FP

    auto tr_0_0_xz_yyz_x = pbuffer.data(idx_op_geom_020_fp + 81);

    auto tr_0_0_xz_yyz_y = pbuffer.data(idx_op_geom_020_fp + 82);

    auto tr_0_0_xz_yyz_z = pbuffer.data(idx_op_geom_020_fp + 83);

    #pragma omp simd aligned(tr_0_0_xz_yyz_x, tr_0_0_xz_yyz_y, tr_0_0_xz_yyz_z, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_yy_0, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyz_x[i] = tr_yy_0[i] - 2.0 * tr_yy_xx[i] * tke_0 - 2.0 * tr_yyz_z[i] * tke_0 + 4.0 * tr_yyz_xxz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_xx[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_y[i] = -2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yyz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_z[i] = -2.0 * tr_yy_xz[i] * tke_0 - 2.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyz_xzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 84-87 components of targeted buffer : FP

    auto tr_0_0_xz_yzz_x = pbuffer.data(idx_op_geom_020_fp + 84);

    auto tr_0_0_xz_yzz_y = pbuffer.data(idx_op_geom_020_fp + 85);

    auto tr_0_0_xz_yzz_z = pbuffer.data(idx_op_geom_020_fp + 86);

    #pragma omp simd aligned(tr_0_0_xz_yzz_x, tr_0_0_xz_yzz_y, tr_0_0_xz_yzz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_yz_0, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yzz_x[i] = 2.0 * tr_yz_0[i] - 4.0 * tr_yz_xx[i] * tke_0 - 2.0 * tr_yzz_z[i] * tke_0 + 4.0 * tr_yzz_xxz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_xx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_y[i] = -4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_z[i] = -4.0 * tr_yz_xz[i] * tke_0 - 2.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yzz_xzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 87-90 components of targeted buffer : FP

    auto tr_0_0_xz_zzz_x = pbuffer.data(idx_op_geom_020_fp + 87);

    auto tr_0_0_xz_zzz_y = pbuffer.data(idx_op_geom_020_fp + 88);

    auto tr_0_0_xz_zzz_z = pbuffer.data(idx_op_geom_020_fp + 89);

    #pragma omp simd aligned(tr_0_0_xz_zzz_x, tr_0_0_xz_zzz_y, tr_0_0_xz_zzz_z, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_xzzz_xz, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, tr_zz_0, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zzz_x, tr_zzz_xxz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_zzz_x[i] = 3.0 * tr_zz_0[i] - 6.0 * tr_zz_xx[i] * tke_0 - 2.0 * tr_zzz_z[i] * tke_0 + 4.0 * tr_zzz_xxz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_0[i] * tbe_0 + 4.0 * tr_zzzz_xx[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_y[i] = -6.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_z[i] = -6.0 * tr_zz_xz[i] * tke_0 - 2.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_zzz_xzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 90-93 components of targeted buffer : FP

    auto tr_0_0_yy_xxx_x = pbuffer.data(idx_op_geom_020_fp + 90);

    auto tr_0_0_yy_xxx_y = pbuffer.data(idx_op_geom_020_fp + 91);

    auto tr_0_0_yy_xxx_z = pbuffer.data(idx_op_geom_020_fp + 92);

    #pragma omp simd aligned(tr_0_0_yy_xxx_x, tr_0_0_yy_xxx_y, tr_0_0_yy_xxx_z, tr_xxx_x, tr_xxx_xyy, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_z, tr_xxxy_0, tr_xxxy_xy, tr_xxxy_yy, tr_xxxy_yz, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxx_x[i] = -2.0 * tr_xxx_x[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxx_xyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_y[i] = -2.0 * tr_xxx_y[i] * tbe_0 - 6.0 * tr_xxx_y[i] * tke_0 + 4.0 * tr_xxx_yyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_0[i] * tbe_0 + 8.0 * tr_xxxy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_z[i] = -2.0 * tr_xxx_z[i] * tbe_0 - 2.0 * tr_xxx_z[i] * tke_0 + 4.0 * tr_xxx_yyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 93-96 components of targeted buffer : FP

    auto tr_0_0_yy_xxy_x = pbuffer.data(idx_op_geom_020_fp + 93);

    auto tr_0_0_yy_xxy_y = pbuffer.data(idx_op_geom_020_fp + 94);

    auto tr_0_0_yy_xxy_z = pbuffer.data(idx_op_geom_020_fp + 95);

    #pragma omp simd aligned(tr_0_0_yy_xxy_x, tr_0_0_yy_xxy_y, tr_0_0_yy_xxy_z, tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxy_x, tr_xxy_xyy, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_z, tr_xxyy_0, tr_xxyy_xy, tr_xxyy_yy, tr_xxyy_yz, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxy_x[i] = -4.0 * tr_xx_xy[i] * tke_0 - 6.0 * tr_xxy_x[i] * tbe_0 - 2.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxy_xyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_y[i] = 2.0 * tr_xx_0[i] - 4.0 * tr_xx_yy[i] * tke_0 - 6.0 * tr_xxy_y[i] * tbe_0 - 6.0 * tr_xxy_y[i] * tke_0 + 4.0 * tr_xxy_yyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_0[i] * tbe_0 + 8.0 * tr_xxyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_z[i] = -4.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxy_z[i] * tke_0 + 4.0 * tr_xxy_yyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 96-99 components of targeted buffer : FP

    auto tr_0_0_yy_xxz_x = pbuffer.data(idx_op_geom_020_fp + 96);

    auto tr_0_0_yy_xxz_y = pbuffer.data(idx_op_geom_020_fp + 97);

    auto tr_0_0_yy_xxz_z = pbuffer.data(idx_op_geom_020_fp + 98);

    #pragma omp simd aligned(tr_0_0_yy_xxz_x, tr_0_0_yy_xxz_y, tr_0_0_yy_xxz_z, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_xyy, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxz_x[i] = -2.0 * tr_xxz_x[i] * tbe_0 - 2.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxz_xyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_y[i] = -2.0 * tr_xxz_y[i] * tbe_0 - 6.0 * tr_xxz_y[i] * tke_0 + 4.0 * tr_xxz_yyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 8.0 * tr_xxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_z[i] = -2.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxz_z[i] * tke_0 + 4.0 * tr_xxz_yyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 99-102 components of targeted buffer : FP

    auto tr_0_0_yy_xyy_x = pbuffer.data(idx_op_geom_020_fp + 99);

    auto tr_0_0_yy_xyy_y = pbuffer.data(idx_op_geom_020_fp + 100);

    auto tr_0_0_yy_xyy_z = pbuffer.data(idx_op_geom_020_fp + 101);

    #pragma omp simd aligned(tr_0_0_yy_xyy_x, tr_0_0_yy_xyy_y, tr_0_0_yy_xyy_z, tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_xyy, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_z, tr_xyyy_0, tr_xyyy_xy, tr_xyyy_yy, tr_xyyy_yz, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyy_x[i] = 2.0 * tr_x_x[i] - 8.0 * tr_xy_xy[i] * tke_0 - 10.0 * tr_xyy_x[i] * tbe_0 - 2.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyy_xyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_y[i] = 2.0 * tr_x_y[i] + 4.0 * tr_xy_0[i] - 8.0 * tr_xy_yy[i] * tke_0 - 10.0 * tr_xyy_y[i] * tbe_0 - 6.0 * tr_xyy_y[i] * tke_0 + 4.0 * tr_xyy_yyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_0[i] * tbe_0 + 8.0 * tr_xyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_z[i] = 2.0 * tr_x_z[i] - 8.0 * tr_xy_yz[i] * tke_0 - 10.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyy_z[i] * tke_0 + 4.0 * tr_xyy_yyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 102-105 components of targeted buffer : FP

    auto tr_0_0_yy_xyz_x = pbuffer.data(idx_op_geom_020_fp + 102);

    auto tr_0_0_yy_xyz_y = pbuffer.data(idx_op_geom_020_fp + 103);

    auto tr_0_0_yy_xyz_z = pbuffer.data(idx_op_geom_020_fp + 104);

    #pragma omp simd aligned(tr_0_0_yy_xyz_x, tr_0_0_yy_xyz_y, tr_0_0_yy_xyz_z, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_xyy, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyz_x[i] = -4.0 * tr_xz_xy[i] * tke_0 - 6.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyz_xyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_y[i] = 2.0 * tr_xz_0[i] - 4.0 * tr_xz_yy[i] * tke_0 - 6.0 * tr_xyz_y[i] * tbe_0 - 6.0 * tr_xyz_y[i] * tke_0 + 4.0 * tr_xyz_yyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 + 8.0 * tr_xyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_z[i] = -4.0 * tr_xz_yz[i] * tke_0 - 6.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyz_z[i] * tke_0 + 4.0 * tr_xyz_yyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 105-108 components of targeted buffer : FP

    auto tr_0_0_yy_xzz_x = pbuffer.data(idx_op_geom_020_fp + 105);

    auto tr_0_0_yy_xzz_y = pbuffer.data(idx_op_geom_020_fp + 106);

    auto tr_0_0_yy_xzz_z = pbuffer.data(idx_op_geom_020_fp + 107);

    #pragma omp simd aligned(tr_0_0_yy_xzz_x, tr_0_0_yy_xzz_y, tr_0_0_yy_xzz_z, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_xyy, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xzz_x[i] = -2.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xzz_xyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_y[i] = -2.0 * tr_xzz_y[i] * tbe_0 - 6.0 * tr_xzz_y[i] * tke_0 + 4.0 * tr_xzz_yyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 + 8.0 * tr_xyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_z[i] = -2.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xzz_z[i] * tke_0 + 4.0 * tr_xzz_yyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 108-111 components of targeted buffer : FP

    auto tr_0_0_yy_yyy_x = pbuffer.data(idx_op_geom_020_fp + 108);

    auto tr_0_0_yy_yyy_y = pbuffer.data(idx_op_geom_020_fp + 109);

    auto tr_0_0_yy_yyy_z = pbuffer.data(idx_op_geom_020_fp + 110);

    #pragma omp simd aligned(tr_0_0_yy_yyy_x, tr_0_0_yy_yyy_y, tr_0_0_yy_yyy_z, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_xyy, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_z, tr_yyyy_0, tr_yyyy_xy, tr_yyyy_yy, tr_yyyy_yz, tr_yyyyy_x, tr_yyyyy_y, tr_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyy_x[i] = 6.0 * tr_y_x[i] - 12.0 * tr_yy_xy[i] * tke_0 - 14.0 * tr_yyy_x[i] * tbe_0 - 2.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyy_xyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_y[i] = 6.0 * tr_y_y[i] + 6.0 * tr_yy_0[i] - 12.0 * tr_yy_yy[i] * tke_0 - 14.0 * tr_yyy_y[i] * tbe_0 - 6.0 * tr_yyy_y[i] * tke_0 + 4.0 * tr_yyy_yyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyy_0[i] * tbe_0 + 8.0 * tr_yyyy_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_z[i] = 6.0 * tr_y_z[i] - 12.0 * tr_yy_yz[i] * tke_0 - 14.0 * tr_yyy_z[i] * tbe_0 - 2.0 * tr_yyy_z[i] * tke_0 + 4.0 * tr_yyy_yyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_z[i] * tbe_0 * tbe_0;
    }

    // Set up 111-114 components of targeted buffer : FP

    auto tr_0_0_yy_yyz_x = pbuffer.data(idx_op_geom_020_fp + 111);

    auto tr_0_0_yy_yyz_y = pbuffer.data(idx_op_geom_020_fp + 112);

    auto tr_0_0_yy_yyz_z = pbuffer.data(idx_op_geom_020_fp + 113);

    #pragma omp simd aligned(tr_0_0_yy_yyz_x, tr_0_0_yy_yyz_y, tr_0_0_yy_yyz_z, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xy, tr_yyyz_yy, tr_yyyz_yz, tr_yyz_x, tr_yyz_xyy, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_z_x, tr_z_y, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyz_x[i] = 2.0 * tr_z_x[i] - 8.0 * tr_yz_xy[i] * tke_0 - 10.0 * tr_yyz_x[i] * tbe_0 - 2.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyz_xyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_y[i] = 2.0 * tr_z_y[i] + 4.0 * tr_yz_0[i] - 8.0 * tr_yz_yy[i] * tke_0 - 10.0 * tr_yyz_y[i] * tbe_0 - 6.0 * tr_yyz_y[i] * tke_0 + 4.0 * tr_yyz_yyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_0[i] * tbe_0 + 8.0 * tr_yyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_z[i] = 2.0 * tr_z_z[i] - 8.0 * tr_yz_yz[i] * tke_0 - 10.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyz_z[i] * tke_0 + 4.0 * tr_yyz_yyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 114-117 components of targeted buffer : FP

    auto tr_0_0_yy_yzz_x = pbuffer.data(idx_op_geom_020_fp + 114);

    auto tr_0_0_yy_yzz_y = pbuffer.data(idx_op_geom_020_fp + 115);

    auto tr_0_0_yy_yzz_z = pbuffer.data(idx_op_geom_020_fp + 116);

    #pragma omp simd aligned(tr_0_0_yy_yzz_x, tr_0_0_yy_yzz_y, tr_0_0_yy_yzz_z, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yyzz_0, tr_yyzz_xy, tr_yyzz_yy, tr_yyzz_yz, tr_yzz_x, tr_yzz_xyy, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yzz_x[i] = -4.0 * tr_zz_xy[i] * tke_0 - 6.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yzz_xyy[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_y[i] = 2.0 * tr_zz_0[i] - 4.0 * tr_zz_yy[i] * tke_0 - 6.0 * tr_yzz_y[i] * tbe_0 - 6.0 * tr_yzz_y[i] * tke_0 + 4.0 * tr_yzz_yyy[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_0[i] * tbe_0 + 8.0 * tr_yyzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_z[i] = -4.0 * tr_zz_yz[i] * tke_0 - 6.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yzz_z[i] * tke_0 + 4.0 * tr_yzz_yyz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 117-120 components of targeted buffer : FP

    auto tr_0_0_yy_zzz_x = pbuffer.data(idx_op_geom_020_fp + 117);

    auto tr_0_0_yy_zzz_y = pbuffer.data(idx_op_geom_020_fp + 118);

    auto tr_0_0_yy_zzz_z = pbuffer.data(idx_op_geom_020_fp + 119);

    #pragma omp simd aligned(tr_0_0_yy_zzz_x, tr_0_0_yy_zzz_y, tr_0_0_yy_zzz_z, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_yzzz_0, tr_yzzz_xy, tr_yzzz_yy, tr_yzzz_yz, tr_zzz_x, tr_zzz_xyy, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_zzz_x[i] = -2.0 * tr_zzz_x[i] * tbe_0 - 2.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_zzz_xyy[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_y[i] = -2.0 * tr_zzz_y[i] * tbe_0 - 6.0 * tr_zzz_y[i] * tke_0 + 4.0 * tr_zzz_yyy[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_0[i] * tbe_0 + 8.0 * tr_yzzz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_z[i] = -2.0 * tr_zzz_z[i] * tbe_0 - 2.0 * tr_zzz_z[i] * tke_0 + 4.0 * tr_zzz_yyz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 120-123 components of targeted buffer : FP

    auto tr_0_0_yz_xxx_x = pbuffer.data(idx_op_geom_020_fp + 120);

    auto tr_0_0_yz_xxx_y = pbuffer.data(idx_op_geom_020_fp + 121);

    auto tr_0_0_yz_xxx_z = pbuffer.data(idx_op_geom_020_fp + 122);

    #pragma omp simd aligned(tr_0_0_yz_xxx_x, tr_0_0_yz_xxx_y, tr_0_0_yz_xxx_z, tr_xxx_xyz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxy_0, tr_xxxy_xz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xy, tr_xxxz_yy, tr_xxxz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxx_x[i] = 4.0 * tr_xxx_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_y[i] = -2.0 * tr_xxx_z[i] * tke_0 + 4.0 * tr_xxx_yyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_z[i] = -2.0 * tr_xxx_y[i] * tke_0 + 4.0 * tr_xxx_yzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 123-126 components of targeted buffer : FP

    auto tr_0_0_yz_xxy_x = pbuffer.data(idx_op_geom_020_fp + 123);

    auto tr_0_0_yz_xxy_y = pbuffer.data(idx_op_geom_020_fp + 124);

    auto tr_0_0_yz_xxy_z = pbuffer.data(idx_op_geom_020_fp + 125);

    #pragma omp simd aligned(tr_0_0_yz_xxy_x, tr_0_0_yz_xxy_y, tr_0_0_yz_xxy_z, tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_0, tr_xxyy_xz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xy, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_y, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxy_x[i] = -2.0 * tr_xx_xz[i] * tke_0 - 2.0 * tr_xxz_x[i] * tbe_0 + 4.0 * tr_xxy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_y[i] = -2.0 * tr_xx_yz[i] * tke_0 - 2.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxy_z[i] * tke_0 + 4.0 * tr_xxy_yyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_z[i] = tr_xx_0[i] - 2.0 * tr_xx_zz[i] * tke_0 - 2.0 * tr_xxz_z[i] * tbe_0 - 2.0 * tr_xxy_y[i] * tke_0 + 4.0 * tr_xxy_yzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 126-129 components of targeted buffer : FP

    auto tr_0_0_yz_xxz_x = pbuffer.data(idx_op_geom_020_fp + 126);

    auto tr_0_0_yz_xxz_y = pbuffer.data(idx_op_geom_020_fp + 127);

    auto tr_0_0_yz_xxz_z = pbuffer.data(idx_op_geom_020_fp + 128);

    #pragma omp simd aligned(tr_0_0_yz_xxz_x, tr_0_0_yz_xxz_y, tr_0_0_yz_xxz_z, tr_xx_0, tr_xx_xy, tr_xx_yy, tr_xx_yz, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxzz_0, tr_xxzz_xy, tr_xxzz_yy, tr_xxzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxz_x[i] = -2.0 * tr_xx_xy[i] * tke_0 + 4.0 * tr_xxz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_x[i] * tbe_0 + 4.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_y[i] = tr_xx_0[i] - 2.0 * tr_xx_yy[i] * tke_0 - 2.0 * tr_xxz_z[i] * tke_0 + 4.0 * tr_xxz_yyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_y[i] * tbe_0 + 4.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_z[i] = -2.0 * tr_xx_yz[i] * tke_0 - 2.0 * tr_xxz_y[i] * tke_0 + 4.0 * tr_xxz_yzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_z[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 129-132 components of targeted buffer : FP

    auto tr_0_0_yz_xyy_x = pbuffer.data(idx_op_geom_020_fp + 129);

    auto tr_0_0_yz_xyy_y = pbuffer.data(idx_op_geom_020_fp + 130);

    auto tr_0_0_yz_xyy_z = pbuffer.data(idx_op_geom_020_fp + 131);

    #pragma omp simd aligned(tr_0_0_yz_xyy_x, tr_0_0_yz_xyy_y, tr_0_0_yz_xyy_z, tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_0, tr_xyyy_xz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xy, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_y, tr_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyy_x[i] = -4.0 * tr_xy_xz[i] * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_y[i] = -4.0 * tr_xy_yz[i] * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyy_z[i] * tke_0 + 4.0 * tr_xyy_yyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_z[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xy_zz[i] * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyy_y[i] * tke_0 + 4.0 * tr_xyy_yzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 132-135 components of targeted buffer : FP

    auto tr_0_0_yz_xyz_x = pbuffer.data(idx_op_geom_020_fp + 132);

    auto tr_0_0_yz_xyz_y = pbuffer.data(idx_op_geom_020_fp + 133);

    auto tr_0_0_yz_xyz_z = pbuffer.data(idx_op_geom_020_fp + 134);

    #pragma omp simd aligned(tr_0_0_yz_xyz_x, tr_0_0_yz_xyz_y, tr_0_0_yz_xyz_z, tr_x_x, tr_x_y, tr_x_z, tr_xy_0, tr_xy_xy, tr_xy_yy, tr_xy_yz, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyzz_0, tr_xyzz_xy, tr_xyzz_yy, tr_xyzz_yz, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_y, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyz_x[i] = tr_x_x[i] - 2.0 * tr_xz_xz[i] * tke_0 - 2.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xy_xy[i] * tke_0 + 4.0 * tr_xyz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_x[i] * tbe_0 + 4.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_y[i] = tr_x_y[i] - 2.0 * tr_xz_yz[i] * tke_0 - 2.0 * tr_xzz_y[i] * tbe_0 + tr_xy_0[i] - 2.0 * tr_xy_yy[i] * tke_0 - 2.0 * tr_xyz_z[i] * tke_0 + 4.0 * tr_xyz_yyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_y[i] * tbe_0 + 4.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_z[i] = tr_x_z[i] + tr_xz_0[i] - 2.0 * tr_xz_zz[i] * tke_0 - 2.0 * tr_xzz_z[i] * tbe_0 - 2.0 * tr_xy_yz[i] * tke_0 - 2.0 * tr_xyz_y[i] * tke_0 + 4.0 * tr_xyz_yzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_z[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 135-138 components of targeted buffer : FP

    auto tr_0_0_yz_xzz_x = pbuffer.data(idx_op_geom_020_fp + 135);

    auto tr_0_0_yz_xzz_y = pbuffer.data(idx_op_geom_020_fp + 136);

    auto tr_0_0_yz_xzz_z = pbuffer.data(idx_op_geom_020_fp + 137);

    #pragma omp simd aligned(tr_0_0_yz_xzz_x, tr_0_0_yz_xzz_y, tr_0_0_yz_xzz_z, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_xz_0, tr_xz_xy, tr_xz_yy, tr_xz_yz, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzzz_0, tr_xzzz_xy, tr_xzzz_yy, tr_xzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xzz_x[i] = -4.0 * tr_xz_xy[i] * tke_0 + 4.0 * tr_xzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_x[i] * tbe_0 + 4.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_y[i] = 2.0 * tr_xz_0[i] - 4.0 * tr_xz_yy[i] * tke_0 - 2.0 * tr_xzz_z[i] * tke_0 + 4.0 * tr_xzz_yyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_y[i] * tbe_0 + 4.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_z[i] = -4.0 * tr_xz_yz[i] * tke_0 - 2.0 * tr_xzz_y[i] * tke_0 + 4.0 * tr_xzz_yzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_z[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 138-141 components of targeted buffer : FP

    auto tr_0_0_yz_yyy_x = pbuffer.data(idx_op_geom_020_fp + 138);

    auto tr_0_0_yz_yyy_y = pbuffer.data(idx_op_geom_020_fp + 139);

    auto tr_0_0_yz_yyy_z = pbuffer.data(idx_op_geom_020_fp + 140);

    #pragma omp simd aligned(tr_0_0_yz_yyy_x, tr_0_0_yz_yyy_y, tr_0_0_yz_yyy_z, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyy_xyz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_0, tr_yyyy_xz, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xy, tr_yyyz_yy, tr_yyyz_yz, tr_yyz_x, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyy_x[i] = -6.0 * tr_yy_xz[i] * tke_0 - 6.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyy_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_y[i] = -6.0 * tr_yy_yz[i] * tke_0 - 6.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_yyy_z[i] * tke_0 + 4.0 * tr_yyy_yyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_yy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_z[i] = 3.0 * tr_yy_0[i] - 6.0 * tr_yy_zz[i] * tke_0 - 6.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyy_y[i] * tke_0 + 4.0 * tr_yyy_yzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_0[i] * tbe_0 + 4.0 * tr_yyyy_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 141-144 components of targeted buffer : FP

    auto tr_0_0_yz_yyz_x = pbuffer.data(idx_op_geom_020_fp + 141);

    auto tr_0_0_yz_yyz_y = pbuffer.data(idx_op_geom_020_fp + 142);

    auto tr_0_0_yz_yyz_z = pbuffer.data(idx_op_geom_020_fp + 143);

    #pragma omp simd aligned(tr_0_0_yz_yyz_x, tr_0_0_yz_yyz_y, tr_0_0_yz_yyz_z, tr_y_x, tr_y_y, tr_y_z, tr_yy_0, tr_yy_xy, tr_yy_yy, tr_yy_yz, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyz_0, tr_yyyz_xz, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyzz_0, tr_yyzz_xy, tr_yyzz_yy, tr_yyzz_yz, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyz_x[i] = 2.0 * tr_y_x[i] - 4.0 * tr_yz_xz[i] * tke_0 - 4.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_yy_xy[i] * tke_0 + 4.0 * tr_yyz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_x[i] * tbe_0 + 4.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_y[i] = 2.0 * tr_y_y[i] - 4.0 * tr_yz_yz[i] * tke_0 - 4.0 * tr_yzz_y[i] * tbe_0 + tr_yy_0[i] - 2.0 * tr_yy_yy[i] * tke_0 - 2.0 * tr_yyz_z[i] * tke_0 + 4.0 * tr_yyz_yyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_yy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_y[i] * tbe_0 + 4.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_z[i] = 2.0 * tr_y_z[i] + 2.0 * tr_yz_0[i] - 4.0 * tr_yz_zz[i] * tke_0 - 4.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yy_yz[i] * tke_0 - 2.0 * tr_yyz_y[i] * tke_0 + 4.0 * tr_yyz_yzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_z[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 144-147 components of targeted buffer : FP

    auto tr_0_0_yz_yzz_x = pbuffer.data(idx_op_geom_020_fp + 144);

    auto tr_0_0_yz_yzz_y = pbuffer.data(idx_op_geom_020_fp + 145);

    auto tr_0_0_yz_yzz_z = pbuffer.data(idx_op_geom_020_fp + 146);

    #pragma omp simd aligned(tr_0_0_yz_yzz_x, tr_0_0_yz_yzz_y, tr_0_0_yz_yzz_z, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yyzz_xz, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_yz_0, tr_yz_xy, tr_yz_yy, tr_yz_yz, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzzz_0, tr_yzzz_xy, tr_yzzz_yy, tr_yzzz_yz, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yzz_x[i] = 2.0 * tr_z_x[i] - 2.0 * tr_zz_xz[i] * tke_0 - 2.0 * tr_zzz_x[i] * tbe_0 - 4.0 * tr_yz_xy[i] * tke_0 + 4.0 * tr_yzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_x[i] * tbe_0 + 4.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_y[i] = 2.0 * tr_z_y[i] - 2.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_zzz_y[i] * tbe_0 + 2.0 * tr_yz_0[i] - 4.0 * tr_yz_yy[i] * tke_0 - 2.0 * tr_yzz_z[i] * tke_0 + 4.0 * tr_yzz_yyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_yy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_y[i] * tbe_0 + 4.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_z[i] = 2.0 * tr_z_z[i] + tr_zz_0[i] - 2.0 * tr_zz_zz[i] * tke_0 - 2.0 * tr_zzz_z[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tke_0 - 2.0 * tr_yzz_y[i] * tke_0 + 4.0 * tr_yzz_yzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_z[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 147-150 components of targeted buffer : FP

    auto tr_0_0_yz_zzz_x = pbuffer.data(idx_op_geom_020_fp + 147);

    auto tr_0_0_yz_zzz_y = pbuffer.data(idx_op_geom_020_fp + 148);

    auto tr_0_0_yz_zzz_z = pbuffer.data(idx_op_geom_020_fp + 149);

    #pragma omp simd aligned(tr_0_0_yz_zzz_x, tr_0_0_yz_zzz_y, tr_0_0_yz_zzz_z, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_yzzz_xz, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, tr_zz_0, tr_zz_xy, tr_zz_yy, tr_zz_yz, tr_zzz_xyz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzzz_0, tr_zzzz_xy, tr_zzzz_yy, tr_zzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_zzz_x[i] = -6.0 * tr_zz_xy[i] * tke_0 + 4.0 * tr_zzz_xyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_y[i] = 3.0 * tr_zz_0[i] - 6.0 * tr_zz_yy[i] * tke_0 - 2.0 * tr_zzz_z[i] * tke_0 + 4.0 * tr_zzz_yyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_0[i] * tbe_0 + 4.0 * tr_zzzz_yy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_z[i] = -6.0 * tr_zz_yz[i] * tke_0 - 2.0 * tr_zzz_y[i] * tke_0 + 4.0 * tr_zzz_yzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_yz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_z[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 150-153 components of targeted buffer : FP

    auto tr_0_0_zz_xxx_x = pbuffer.data(idx_op_geom_020_fp + 150);

    auto tr_0_0_zz_xxx_y = pbuffer.data(idx_op_geom_020_fp + 151);

    auto tr_0_0_zz_xxx_z = pbuffer.data(idx_op_geom_020_fp + 152);

    #pragma omp simd aligned(tr_0_0_zz_xxx_x, tr_0_0_zz_xxx_y, tr_0_0_zz_xxx_z, tr_xxx_x, tr_xxx_xzz, tr_xxx_y, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxz_0, tr_xxxz_xz, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxx_x[i] = -2.0 * tr_xxx_x[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxx_xzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_y[i] = -2.0 * tr_xxx_y[i] * tbe_0 - 2.0 * tr_xxx_y[i] * tke_0 + 4.0 * tr_xxx_yzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_z[i] = -2.0 * tr_xxx_z[i] * tbe_0 - 6.0 * tr_xxx_z[i] * tke_0 + 4.0 * tr_xxx_zzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_0[i] * tbe_0 + 8.0 * tr_xxxz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 153-156 components of targeted buffer : FP

    auto tr_0_0_zz_xxy_x = pbuffer.data(idx_op_geom_020_fp + 153);

    auto tr_0_0_zz_xxy_y = pbuffer.data(idx_op_geom_020_fp + 154);

    auto tr_0_0_zz_xxy_z = pbuffer.data(idx_op_geom_020_fp + 155);

    #pragma omp simd aligned(tr_0_0_zz_xxy_x, tr_0_0_zz_xxy_y, tr_0_0_zz_xxy_z, tr_xxy_x, tr_xxy_xzz, tr_xxy_y, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_0, tr_xxyz_xz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxy_x[i] = -2.0 * tr_xxy_x[i] * tbe_0 - 2.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxy_xzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_y[i] = -2.0 * tr_xxy_y[i] * tbe_0 - 2.0 * tr_xxy_y[i] * tke_0 + 4.0 * tr_xxy_yzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_z[i] = -2.0 * tr_xxy_z[i] * tbe_0 - 6.0 * tr_xxy_z[i] * tke_0 + 4.0 * tr_xxy_zzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 8.0 * tr_xxyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 156-159 components of targeted buffer : FP

    auto tr_0_0_zz_xxz_x = pbuffer.data(idx_op_geom_020_fp + 156);

    auto tr_0_0_zz_xxz_y = pbuffer.data(idx_op_geom_020_fp + 157);

    auto tr_0_0_zz_xxz_z = pbuffer.data(idx_op_geom_020_fp + 158);

    #pragma omp simd aligned(tr_0_0_zz_xxz_x, tr_0_0_zz_xxz_y, tr_0_0_zz_xxz_z, tr_xx_0, tr_xx_xz, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_xzz, tr_xxz_y, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_0, tr_xxzz_xz, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxz_x[i] = -4.0 * tr_xx_xz[i] * tke_0 - 6.0 * tr_xxz_x[i] * tbe_0 - 2.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxz_xzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_y[i] = -4.0 * tr_xx_yz[i] * tke_0 - 6.0 * tr_xxz_y[i] * tbe_0 - 2.0 * tr_xxz_y[i] * tke_0 + 4.0 * tr_xxz_yzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_z[i] = 2.0 * tr_xx_0[i] - 4.0 * tr_xx_zz[i] * tke_0 - 6.0 * tr_xxz_z[i] * tbe_0 - 6.0 * tr_xxz_z[i] * tke_0 + 4.0 * tr_xxz_zzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_0[i] * tbe_0 + 8.0 * tr_xxzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 159-162 components of targeted buffer : FP

    auto tr_0_0_zz_xyy_x = pbuffer.data(idx_op_geom_020_fp + 159);

    auto tr_0_0_zz_xyy_y = pbuffer.data(idx_op_geom_020_fp + 160);

    auto tr_0_0_zz_xyy_z = pbuffer.data(idx_op_geom_020_fp + 161);

    #pragma omp simd aligned(tr_0_0_zz_xyy_x, tr_0_0_zz_xyy_y, tr_0_0_zz_xyy_z, tr_xyy_x, tr_xyy_xzz, tr_xyy_y, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_0, tr_xyyz_xz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyy_x[i] = -2.0 * tr_xyy_x[i] * tbe_0 - 2.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyy_xzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_y[i] = -2.0 * tr_xyy_y[i] * tbe_0 - 2.0 * tr_xyy_y[i] * tke_0 + 4.0 * tr_xyy_yzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_z[i] = -2.0 * tr_xyy_z[i] * tbe_0 - 6.0 * tr_xyy_z[i] * tke_0 + 4.0 * tr_xyy_zzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 + 8.0 * tr_xyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 162-165 components of targeted buffer : FP

    auto tr_0_0_zz_xyz_x = pbuffer.data(idx_op_geom_020_fp + 162);

    auto tr_0_0_zz_xyz_y = pbuffer.data(idx_op_geom_020_fp + 163);

    auto tr_0_0_zz_xyz_z = pbuffer.data(idx_op_geom_020_fp + 164);

    #pragma omp simd aligned(tr_0_0_zz_xyz_x, tr_0_0_zz_xyz_y, tr_0_0_zz_xyz_z, tr_xy_0, tr_xy_xz, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xzz, tr_xyz_y, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_0, tr_xyzz_xz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyz_x[i] = -4.0 * tr_xy_xz[i] * tke_0 - 6.0 * tr_xyz_x[i] * tbe_0 - 2.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyz_xzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_y[i] = -4.0 * tr_xy_yz[i] * tke_0 - 6.0 * tr_xyz_y[i] * tbe_0 - 2.0 * tr_xyz_y[i] * tke_0 + 4.0 * tr_xyz_yzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_z[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xy_zz[i] * tke_0 - 6.0 * tr_xyz_z[i] * tbe_0 - 6.0 * tr_xyz_z[i] * tke_0 + 4.0 * tr_xyz_zzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 + 8.0 * tr_xyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 165-168 components of targeted buffer : FP

    auto tr_0_0_zz_xzz_x = pbuffer.data(idx_op_geom_020_fp + 165);

    auto tr_0_0_zz_xzz_y = pbuffer.data(idx_op_geom_020_fp + 166);

    auto tr_0_0_zz_xzz_z = pbuffer.data(idx_op_geom_020_fp + 167);

    #pragma omp simd aligned(tr_0_0_zz_xzz_x, tr_0_0_zz_xzz_y, tr_0_0_zz_xzz_z, tr_x_x, tr_x_y, tr_x_z, tr_xz_0, tr_xz_xz, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xzz, tr_xzz_y, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_0, tr_xzzz_xz, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xzz_x[i] = 2.0 * tr_x_x[i] - 8.0 * tr_xz_xz[i] * tke_0 - 10.0 * tr_xzz_x[i] * tbe_0 - 2.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xzz_xzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_y[i] = 2.0 * tr_x_y[i] - 8.0 * tr_xz_yz[i] * tke_0 - 10.0 * tr_xzz_y[i] * tbe_0 - 2.0 * tr_xzz_y[i] * tke_0 + 4.0 * tr_xzz_yzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_z[i] = 2.0 * tr_x_z[i] + 4.0 * tr_xz_0[i] - 8.0 * tr_xz_zz[i] * tke_0 - 10.0 * tr_xzz_z[i] * tbe_0 - 6.0 * tr_xzz_z[i] * tke_0 + 4.0 * tr_xzz_zzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_0[i] * tbe_0 + 8.0 * tr_xzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 168-171 components of targeted buffer : FP

    auto tr_0_0_zz_yyy_x = pbuffer.data(idx_op_geom_020_fp + 168);

    auto tr_0_0_zz_yyy_y = pbuffer.data(idx_op_geom_020_fp + 169);

    auto tr_0_0_zz_yyy_z = pbuffer.data(idx_op_geom_020_fp + 170);

    #pragma omp simd aligned(tr_0_0_zz_yyy_x, tr_0_0_zz_yyy_y, tr_0_0_zz_yyy_z, tr_yyy_x, tr_yyy_xzz, tr_yyy_y, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyz_0, tr_yyyz_xz, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyy_x[i] = -2.0 * tr_yyy_x[i] * tbe_0 - 2.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyy_xzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_y[i] = -2.0 * tr_yyy_y[i] * tbe_0 - 2.0 * tr_yyy_y[i] * tke_0 + 4.0 * tr_yyy_yzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_z[i] = -2.0 * tr_yyy_z[i] * tbe_0 - 6.0 * tr_yyy_z[i] * tke_0 + 4.0 * tr_yyy_zzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_0[i] * tbe_0 + 8.0 * tr_yyyz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 171-174 components of targeted buffer : FP

    auto tr_0_0_zz_yyz_x = pbuffer.data(idx_op_geom_020_fp + 171);

    auto tr_0_0_zz_yyz_y = pbuffer.data(idx_op_geom_020_fp + 172);

    auto tr_0_0_zz_yyz_z = pbuffer.data(idx_op_geom_020_fp + 173);

    #pragma omp simd aligned(tr_0_0_zz_yyz_x, tr_0_0_zz_yyz_y, tr_0_0_zz_yyz_z, tr_yy_0, tr_yy_xz, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xzz, tr_yyz_y, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_0, tr_yyzz_xz, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyz_x[i] = -4.0 * tr_yy_xz[i] * tke_0 - 6.0 * tr_yyz_x[i] * tbe_0 - 2.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyz_xzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_y[i] = -4.0 * tr_yy_yz[i] * tke_0 - 6.0 * tr_yyz_y[i] * tbe_0 - 2.0 * tr_yyz_y[i] * tke_0 + 4.0 * tr_yyz_yzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_z[i] = 2.0 * tr_yy_0[i] - 4.0 * tr_yy_zz[i] * tke_0 - 6.0 * tr_yyz_z[i] * tbe_0 - 6.0 * tr_yyz_z[i] * tke_0 + 4.0 * tr_yyz_zzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_0[i] * tbe_0 + 8.0 * tr_yyzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 174-177 components of targeted buffer : FP

    auto tr_0_0_zz_yzz_x = pbuffer.data(idx_op_geom_020_fp + 174);

    auto tr_0_0_zz_yzz_y = pbuffer.data(idx_op_geom_020_fp + 175);

    auto tr_0_0_zz_yzz_z = pbuffer.data(idx_op_geom_020_fp + 176);

    #pragma omp simd aligned(tr_0_0_zz_yzz_x, tr_0_0_zz_yzz_y, tr_0_0_zz_yzz_z, tr_y_x, tr_y_y, tr_y_z, tr_yz_0, tr_yz_xz, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xzz, tr_yzz_y, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_0, tr_yzzz_xz, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yzz_x[i] = 2.0 * tr_y_x[i] - 8.0 * tr_yz_xz[i] * tke_0 - 10.0 * tr_yzz_x[i] * tbe_0 - 2.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yzz_xzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_y[i] = 2.0 * tr_y_y[i] - 8.0 * tr_yz_yz[i] * tke_0 - 10.0 * tr_yzz_y[i] * tbe_0 - 2.0 * tr_yzz_y[i] * tke_0 + 4.0 * tr_yzz_yzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_z[i] = 2.0 * tr_y_z[i] + 4.0 * tr_yz_0[i] - 8.0 * tr_yz_zz[i] * tke_0 - 10.0 * tr_yzz_z[i] * tbe_0 - 6.0 * tr_yzz_z[i] * tke_0 + 4.0 * tr_yzz_zzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_0[i] * tbe_0 + 8.0 * tr_yzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_z[i] * tbe_0 * tbe_0;
    }

    // Set up 177-180 components of targeted buffer : FP

    auto tr_0_0_zz_zzz_x = pbuffer.data(idx_op_geom_020_fp + 177);

    auto tr_0_0_zz_zzz_y = pbuffer.data(idx_op_geom_020_fp + 178);

    auto tr_0_0_zz_zzz_z = pbuffer.data(idx_op_geom_020_fp + 179);

    #pragma omp simd aligned(tr_0_0_zz_zzz_x, tr_0_0_zz_zzz_y, tr_0_0_zz_zzz_z, tr_z_x, tr_z_y, tr_z_z, tr_zz_0, tr_zz_xz, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xzz, tr_zzz_y, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_0, tr_zzzz_xz, tr_zzzz_yz, tr_zzzz_zz, tr_zzzzz_x, tr_zzzzz_y, tr_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_zzz_x[i] = 6.0 * tr_z_x[i] - 12.0 * tr_zz_xz[i] * tke_0 - 14.0 * tr_zzz_x[i] * tbe_0 - 2.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_zzz_xzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_xz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_x[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_y[i] = 6.0 * tr_z_y[i] - 12.0 * tr_zz_yz[i] * tke_0 - 14.0 * tr_zzz_y[i] * tbe_0 - 2.0 * tr_zzz_y[i] * tke_0 + 4.0 * tr_zzz_yzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_yz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_y[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_z[i] = 6.0 * tr_z_z[i] + 6.0 * tr_zz_0[i] - 12.0 * tr_zz_zz[i] * tke_0 - 14.0 * tr_zzz_z[i] * tbe_0 - 6.0 * tr_zzz_z[i] * tke_0 + 4.0 * tr_zzz_zzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_0[i] * tbe_0 + 8.0 * tr_zzzz_zz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_z[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

