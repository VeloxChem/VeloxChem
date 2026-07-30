#include "GeometricalDerivatives020ForFF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_ff(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_ff,
                         const int idx_op_pf,
                         const int idx_op_dd,
                         const int idx_op_dg,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_fh,
                         const int idx_op_gd,
                         const int idx_op_gg,
                         const int idx_op_hf,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up components of auxiliary buffer : FH

    auto tr_xxx_xxxxx = pbuffer.data(idx_op_fh);

    auto tr_xxx_xxxxy = pbuffer.data(idx_op_fh + 1);

    auto tr_xxx_xxxxz = pbuffer.data(idx_op_fh + 2);

    auto tr_xxx_xxxyy = pbuffer.data(idx_op_fh + 3);

    auto tr_xxx_xxxyz = pbuffer.data(idx_op_fh + 4);

    auto tr_xxx_xxxzz = pbuffer.data(idx_op_fh + 5);

    auto tr_xxx_xxyyy = pbuffer.data(idx_op_fh + 6);

    auto tr_xxx_xxyyz = pbuffer.data(idx_op_fh + 7);

    auto tr_xxx_xxyzz = pbuffer.data(idx_op_fh + 8);

    auto tr_xxx_xxzzz = pbuffer.data(idx_op_fh + 9);

    auto tr_xxx_xyyyy = pbuffer.data(idx_op_fh + 10);

    auto tr_xxx_xyyyz = pbuffer.data(idx_op_fh + 11);

    auto tr_xxx_xyyzz = pbuffer.data(idx_op_fh + 12);

    auto tr_xxx_xyzzz = pbuffer.data(idx_op_fh + 13);

    auto tr_xxx_xzzzz = pbuffer.data(idx_op_fh + 14);

    auto tr_xxx_yyyyy = pbuffer.data(idx_op_fh + 15);

    auto tr_xxx_yyyyz = pbuffer.data(idx_op_fh + 16);

    auto tr_xxx_yyyzz = pbuffer.data(idx_op_fh + 17);

    auto tr_xxx_yyzzz = pbuffer.data(idx_op_fh + 18);

    auto tr_xxx_yzzzz = pbuffer.data(idx_op_fh + 19);

    auto tr_xxx_zzzzz = pbuffer.data(idx_op_fh + 20);

    auto tr_xxy_xxxxx = pbuffer.data(idx_op_fh + 21);

    auto tr_xxy_xxxxy = pbuffer.data(idx_op_fh + 22);

    auto tr_xxy_xxxxz = pbuffer.data(idx_op_fh + 23);

    auto tr_xxy_xxxyy = pbuffer.data(idx_op_fh + 24);

    auto tr_xxy_xxxyz = pbuffer.data(idx_op_fh + 25);

    auto tr_xxy_xxxzz = pbuffer.data(idx_op_fh + 26);

    auto tr_xxy_xxyyy = pbuffer.data(idx_op_fh + 27);

    auto tr_xxy_xxyyz = pbuffer.data(idx_op_fh + 28);

    auto tr_xxy_xxyzz = pbuffer.data(idx_op_fh + 29);

    auto tr_xxy_xxzzz = pbuffer.data(idx_op_fh + 30);

    auto tr_xxy_xyyyy = pbuffer.data(idx_op_fh + 31);

    auto tr_xxy_xyyyz = pbuffer.data(idx_op_fh + 32);

    auto tr_xxy_xyyzz = pbuffer.data(idx_op_fh + 33);

    auto tr_xxy_xyzzz = pbuffer.data(idx_op_fh + 34);

    auto tr_xxy_xzzzz = pbuffer.data(idx_op_fh + 35);

    auto tr_xxy_yyyyy = pbuffer.data(idx_op_fh + 36);

    auto tr_xxy_yyyyz = pbuffer.data(idx_op_fh + 37);

    auto tr_xxy_yyyzz = pbuffer.data(idx_op_fh + 38);

    auto tr_xxy_yyzzz = pbuffer.data(idx_op_fh + 39);

    auto tr_xxy_yzzzz = pbuffer.data(idx_op_fh + 40);

    auto tr_xxy_zzzzz = pbuffer.data(idx_op_fh + 41);

    auto tr_xxz_xxxxx = pbuffer.data(idx_op_fh + 42);

    auto tr_xxz_xxxxy = pbuffer.data(idx_op_fh + 43);

    auto tr_xxz_xxxxz = pbuffer.data(idx_op_fh + 44);

    auto tr_xxz_xxxyy = pbuffer.data(idx_op_fh + 45);

    auto tr_xxz_xxxyz = pbuffer.data(idx_op_fh + 46);

    auto tr_xxz_xxxzz = pbuffer.data(idx_op_fh + 47);

    auto tr_xxz_xxyyy = pbuffer.data(idx_op_fh + 48);

    auto tr_xxz_xxyyz = pbuffer.data(idx_op_fh + 49);

    auto tr_xxz_xxyzz = pbuffer.data(idx_op_fh + 50);

    auto tr_xxz_xxzzz = pbuffer.data(idx_op_fh + 51);

    auto tr_xxz_xyyyy = pbuffer.data(idx_op_fh + 52);

    auto tr_xxz_xyyyz = pbuffer.data(idx_op_fh + 53);

    auto tr_xxz_xyyzz = pbuffer.data(idx_op_fh + 54);

    auto tr_xxz_xyzzz = pbuffer.data(idx_op_fh + 55);

    auto tr_xxz_xzzzz = pbuffer.data(idx_op_fh + 56);

    auto tr_xxz_yyyyy = pbuffer.data(idx_op_fh + 57);

    auto tr_xxz_yyyyz = pbuffer.data(idx_op_fh + 58);

    auto tr_xxz_yyyzz = pbuffer.data(idx_op_fh + 59);

    auto tr_xxz_yyzzz = pbuffer.data(idx_op_fh + 60);

    auto tr_xxz_yzzzz = pbuffer.data(idx_op_fh + 61);

    auto tr_xxz_zzzzz = pbuffer.data(idx_op_fh + 62);

    auto tr_xyy_xxxxx = pbuffer.data(idx_op_fh + 63);

    auto tr_xyy_xxxxy = pbuffer.data(idx_op_fh + 64);

    auto tr_xyy_xxxxz = pbuffer.data(idx_op_fh + 65);

    auto tr_xyy_xxxyy = pbuffer.data(idx_op_fh + 66);

    auto tr_xyy_xxxyz = pbuffer.data(idx_op_fh + 67);

    auto tr_xyy_xxxzz = pbuffer.data(idx_op_fh + 68);

    auto tr_xyy_xxyyy = pbuffer.data(idx_op_fh + 69);

    auto tr_xyy_xxyyz = pbuffer.data(idx_op_fh + 70);

    auto tr_xyy_xxyzz = pbuffer.data(idx_op_fh + 71);

    auto tr_xyy_xxzzz = pbuffer.data(idx_op_fh + 72);

    auto tr_xyy_xyyyy = pbuffer.data(idx_op_fh + 73);

    auto tr_xyy_xyyyz = pbuffer.data(idx_op_fh + 74);

    auto tr_xyy_xyyzz = pbuffer.data(idx_op_fh + 75);

    auto tr_xyy_xyzzz = pbuffer.data(idx_op_fh + 76);

    auto tr_xyy_xzzzz = pbuffer.data(idx_op_fh + 77);

    auto tr_xyy_yyyyy = pbuffer.data(idx_op_fh + 78);

    auto tr_xyy_yyyyz = pbuffer.data(idx_op_fh + 79);

    auto tr_xyy_yyyzz = pbuffer.data(idx_op_fh + 80);

    auto tr_xyy_yyzzz = pbuffer.data(idx_op_fh + 81);

    auto tr_xyy_yzzzz = pbuffer.data(idx_op_fh + 82);

    auto tr_xyy_zzzzz = pbuffer.data(idx_op_fh + 83);

    auto tr_xyz_xxxxx = pbuffer.data(idx_op_fh + 84);

    auto tr_xyz_xxxxy = pbuffer.data(idx_op_fh + 85);

    auto tr_xyz_xxxxz = pbuffer.data(idx_op_fh + 86);

    auto tr_xyz_xxxyy = pbuffer.data(idx_op_fh + 87);

    auto tr_xyz_xxxyz = pbuffer.data(idx_op_fh + 88);

    auto tr_xyz_xxxzz = pbuffer.data(idx_op_fh + 89);

    auto tr_xyz_xxyyy = pbuffer.data(idx_op_fh + 90);

    auto tr_xyz_xxyyz = pbuffer.data(idx_op_fh + 91);

    auto tr_xyz_xxyzz = pbuffer.data(idx_op_fh + 92);

    auto tr_xyz_xxzzz = pbuffer.data(idx_op_fh + 93);

    auto tr_xyz_xyyyy = pbuffer.data(idx_op_fh + 94);

    auto tr_xyz_xyyyz = pbuffer.data(idx_op_fh + 95);

    auto tr_xyz_xyyzz = pbuffer.data(idx_op_fh + 96);

    auto tr_xyz_xyzzz = pbuffer.data(idx_op_fh + 97);

    auto tr_xyz_xzzzz = pbuffer.data(idx_op_fh + 98);

    auto tr_xyz_yyyyy = pbuffer.data(idx_op_fh + 99);

    auto tr_xyz_yyyyz = pbuffer.data(idx_op_fh + 100);

    auto tr_xyz_yyyzz = pbuffer.data(idx_op_fh + 101);

    auto tr_xyz_yyzzz = pbuffer.data(idx_op_fh + 102);

    auto tr_xyz_yzzzz = pbuffer.data(idx_op_fh + 103);

    auto tr_xyz_zzzzz = pbuffer.data(idx_op_fh + 104);

    auto tr_xzz_xxxxx = pbuffer.data(idx_op_fh + 105);

    auto tr_xzz_xxxxy = pbuffer.data(idx_op_fh + 106);

    auto tr_xzz_xxxxz = pbuffer.data(idx_op_fh + 107);

    auto tr_xzz_xxxyy = pbuffer.data(idx_op_fh + 108);

    auto tr_xzz_xxxyz = pbuffer.data(idx_op_fh + 109);

    auto tr_xzz_xxxzz = pbuffer.data(idx_op_fh + 110);

    auto tr_xzz_xxyyy = pbuffer.data(idx_op_fh + 111);

    auto tr_xzz_xxyyz = pbuffer.data(idx_op_fh + 112);

    auto tr_xzz_xxyzz = pbuffer.data(idx_op_fh + 113);

    auto tr_xzz_xxzzz = pbuffer.data(idx_op_fh + 114);

    auto tr_xzz_xyyyy = pbuffer.data(idx_op_fh + 115);

    auto tr_xzz_xyyyz = pbuffer.data(idx_op_fh + 116);

    auto tr_xzz_xyyzz = pbuffer.data(idx_op_fh + 117);

    auto tr_xzz_xyzzz = pbuffer.data(idx_op_fh + 118);

    auto tr_xzz_xzzzz = pbuffer.data(idx_op_fh + 119);

    auto tr_xzz_yyyyy = pbuffer.data(idx_op_fh + 120);

    auto tr_xzz_yyyyz = pbuffer.data(idx_op_fh + 121);

    auto tr_xzz_yyyzz = pbuffer.data(idx_op_fh + 122);

    auto tr_xzz_yyzzz = pbuffer.data(idx_op_fh + 123);

    auto tr_xzz_yzzzz = pbuffer.data(idx_op_fh + 124);

    auto tr_xzz_zzzzz = pbuffer.data(idx_op_fh + 125);

    auto tr_yyy_xxxxx = pbuffer.data(idx_op_fh + 126);

    auto tr_yyy_xxxxy = pbuffer.data(idx_op_fh + 127);

    auto tr_yyy_xxxxz = pbuffer.data(idx_op_fh + 128);

    auto tr_yyy_xxxyy = pbuffer.data(idx_op_fh + 129);

    auto tr_yyy_xxxyz = pbuffer.data(idx_op_fh + 130);

    auto tr_yyy_xxxzz = pbuffer.data(idx_op_fh + 131);

    auto tr_yyy_xxyyy = pbuffer.data(idx_op_fh + 132);

    auto tr_yyy_xxyyz = pbuffer.data(idx_op_fh + 133);

    auto tr_yyy_xxyzz = pbuffer.data(idx_op_fh + 134);

    auto tr_yyy_xxzzz = pbuffer.data(idx_op_fh + 135);

    auto tr_yyy_xyyyy = pbuffer.data(idx_op_fh + 136);

    auto tr_yyy_xyyyz = pbuffer.data(idx_op_fh + 137);

    auto tr_yyy_xyyzz = pbuffer.data(idx_op_fh + 138);

    auto tr_yyy_xyzzz = pbuffer.data(idx_op_fh + 139);

    auto tr_yyy_xzzzz = pbuffer.data(idx_op_fh + 140);

    auto tr_yyy_yyyyy = pbuffer.data(idx_op_fh + 141);

    auto tr_yyy_yyyyz = pbuffer.data(idx_op_fh + 142);

    auto tr_yyy_yyyzz = pbuffer.data(idx_op_fh + 143);

    auto tr_yyy_yyzzz = pbuffer.data(idx_op_fh + 144);

    auto tr_yyy_yzzzz = pbuffer.data(idx_op_fh + 145);

    auto tr_yyy_zzzzz = pbuffer.data(idx_op_fh + 146);

    auto tr_yyz_xxxxx = pbuffer.data(idx_op_fh + 147);

    auto tr_yyz_xxxxy = pbuffer.data(idx_op_fh + 148);

    auto tr_yyz_xxxxz = pbuffer.data(idx_op_fh + 149);

    auto tr_yyz_xxxyy = pbuffer.data(idx_op_fh + 150);

    auto tr_yyz_xxxyz = pbuffer.data(idx_op_fh + 151);

    auto tr_yyz_xxxzz = pbuffer.data(idx_op_fh + 152);

    auto tr_yyz_xxyyy = pbuffer.data(idx_op_fh + 153);

    auto tr_yyz_xxyyz = pbuffer.data(idx_op_fh + 154);

    auto tr_yyz_xxyzz = pbuffer.data(idx_op_fh + 155);

    auto tr_yyz_xxzzz = pbuffer.data(idx_op_fh + 156);

    auto tr_yyz_xyyyy = pbuffer.data(idx_op_fh + 157);

    auto tr_yyz_xyyyz = pbuffer.data(idx_op_fh + 158);

    auto tr_yyz_xyyzz = pbuffer.data(idx_op_fh + 159);

    auto tr_yyz_xyzzz = pbuffer.data(idx_op_fh + 160);

    auto tr_yyz_xzzzz = pbuffer.data(idx_op_fh + 161);

    auto tr_yyz_yyyyy = pbuffer.data(idx_op_fh + 162);

    auto tr_yyz_yyyyz = pbuffer.data(idx_op_fh + 163);

    auto tr_yyz_yyyzz = pbuffer.data(idx_op_fh + 164);

    auto tr_yyz_yyzzz = pbuffer.data(idx_op_fh + 165);

    auto tr_yyz_yzzzz = pbuffer.data(idx_op_fh + 166);

    auto tr_yyz_zzzzz = pbuffer.data(idx_op_fh + 167);

    auto tr_yzz_xxxxx = pbuffer.data(idx_op_fh + 168);

    auto tr_yzz_xxxxy = pbuffer.data(idx_op_fh + 169);

    auto tr_yzz_xxxxz = pbuffer.data(idx_op_fh + 170);

    auto tr_yzz_xxxyy = pbuffer.data(idx_op_fh + 171);

    auto tr_yzz_xxxyz = pbuffer.data(idx_op_fh + 172);

    auto tr_yzz_xxxzz = pbuffer.data(idx_op_fh + 173);

    auto tr_yzz_xxyyy = pbuffer.data(idx_op_fh + 174);

    auto tr_yzz_xxyyz = pbuffer.data(idx_op_fh + 175);

    auto tr_yzz_xxyzz = pbuffer.data(idx_op_fh + 176);

    auto tr_yzz_xxzzz = pbuffer.data(idx_op_fh + 177);

    auto tr_yzz_xyyyy = pbuffer.data(idx_op_fh + 178);

    auto tr_yzz_xyyyz = pbuffer.data(idx_op_fh + 179);

    auto tr_yzz_xyyzz = pbuffer.data(idx_op_fh + 180);

    auto tr_yzz_xyzzz = pbuffer.data(idx_op_fh + 181);

    auto tr_yzz_xzzzz = pbuffer.data(idx_op_fh + 182);

    auto tr_yzz_yyyyy = pbuffer.data(idx_op_fh + 183);

    auto tr_yzz_yyyyz = pbuffer.data(idx_op_fh + 184);

    auto tr_yzz_yyyzz = pbuffer.data(idx_op_fh + 185);

    auto tr_yzz_yyzzz = pbuffer.data(idx_op_fh + 186);

    auto tr_yzz_yzzzz = pbuffer.data(idx_op_fh + 187);

    auto tr_yzz_zzzzz = pbuffer.data(idx_op_fh + 188);

    auto tr_zzz_xxxxx = pbuffer.data(idx_op_fh + 189);

    auto tr_zzz_xxxxy = pbuffer.data(idx_op_fh + 190);

    auto tr_zzz_xxxxz = pbuffer.data(idx_op_fh + 191);

    auto tr_zzz_xxxyy = pbuffer.data(idx_op_fh + 192);

    auto tr_zzz_xxxyz = pbuffer.data(idx_op_fh + 193);

    auto tr_zzz_xxxzz = pbuffer.data(idx_op_fh + 194);

    auto tr_zzz_xxyyy = pbuffer.data(idx_op_fh + 195);

    auto tr_zzz_xxyyz = pbuffer.data(idx_op_fh + 196);

    auto tr_zzz_xxyzz = pbuffer.data(idx_op_fh + 197);

    auto tr_zzz_xxzzz = pbuffer.data(idx_op_fh + 198);

    auto tr_zzz_xyyyy = pbuffer.data(idx_op_fh + 199);

    auto tr_zzz_xyyyz = pbuffer.data(idx_op_fh + 200);

    auto tr_zzz_xyyzz = pbuffer.data(idx_op_fh + 201);

    auto tr_zzz_xyzzz = pbuffer.data(idx_op_fh + 202);

    auto tr_zzz_xzzzz = pbuffer.data(idx_op_fh + 203);

    auto tr_zzz_yyyyy = pbuffer.data(idx_op_fh + 204);

    auto tr_zzz_yyyyz = pbuffer.data(idx_op_fh + 205);

    auto tr_zzz_yyyzz = pbuffer.data(idx_op_fh + 206);

    auto tr_zzz_yyzzz = pbuffer.data(idx_op_fh + 207);

    auto tr_zzz_yzzzz = pbuffer.data(idx_op_fh + 208);

    auto tr_zzz_zzzzz = pbuffer.data(idx_op_fh + 209);

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

    // Set up components of auxiliary buffer : GG

    auto tr_xxxx_xxxx = pbuffer.data(idx_op_gg);

    auto tr_xxxx_xxxy = pbuffer.data(idx_op_gg + 1);

    auto tr_xxxx_xxxz = pbuffer.data(idx_op_gg + 2);

    auto tr_xxxx_xxyy = pbuffer.data(idx_op_gg + 3);

    auto tr_xxxx_xxyz = pbuffer.data(idx_op_gg + 4);

    auto tr_xxxx_xxzz = pbuffer.data(idx_op_gg + 5);

    auto tr_xxxx_xyyy = pbuffer.data(idx_op_gg + 6);

    auto tr_xxxx_xyyz = pbuffer.data(idx_op_gg + 7);

    auto tr_xxxx_xyzz = pbuffer.data(idx_op_gg + 8);

    auto tr_xxxx_xzzz = pbuffer.data(idx_op_gg + 9);

    auto tr_xxxx_yyyy = pbuffer.data(idx_op_gg + 10);

    auto tr_xxxx_yyyz = pbuffer.data(idx_op_gg + 11);

    auto tr_xxxx_yyzz = pbuffer.data(idx_op_gg + 12);

    auto tr_xxxx_yzzz = pbuffer.data(idx_op_gg + 13);

    auto tr_xxxx_zzzz = pbuffer.data(idx_op_gg + 14);

    auto tr_xxxy_xxxx = pbuffer.data(idx_op_gg + 15);

    auto tr_xxxy_xxxy = pbuffer.data(idx_op_gg + 16);

    auto tr_xxxy_xxxz = pbuffer.data(idx_op_gg + 17);

    auto tr_xxxy_xxyy = pbuffer.data(idx_op_gg + 18);

    auto tr_xxxy_xxyz = pbuffer.data(idx_op_gg + 19);

    auto tr_xxxy_xxzz = pbuffer.data(idx_op_gg + 20);

    auto tr_xxxy_xyyy = pbuffer.data(idx_op_gg + 21);

    auto tr_xxxy_xyyz = pbuffer.data(idx_op_gg + 22);

    auto tr_xxxy_xyzz = pbuffer.data(idx_op_gg + 23);

    auto tr_xxxy_xzzz = pbuffer.data(idx_op_gg + 24);

    auto tr_xxxy_yyyy = pbuffer.data(idx_op_gg + 25);

    auto tr_xxxy_yyyz = pbuffer.data(idx_op_gg + 26);

    auto tr_xxxy_yyzz = pbuffer.data(idx_op_gg + 27);

    auto tr_xxxy_yzzz = pbuffer.data(idx_op_gg + 28);

    auto tr_xxxy_zzzz = pbuffer.data(idx_op_gg + 29);

    auto tr_xxxz_xxxx = pbuffer.data(idx_op_gg + 30);

    auto tr_xxxz_xxxy = pbuffer.data(idx_op_gg + 31);

    auto tr_xxxz_xxxz = pbuffer.data(idx_op_gg + 32);

    auto tr_xxxz_xxyy = pbuffer.data(idx_op_gg + 33);

    auto tr_xxxz_xxyz = pbuffer.data(idx_op_gg + 34);

    auto tr_xxxz_xxzz = pbuffer.data(idx_op_gg + 35);

    auto tr_xxxz_xyyy = pbuffer.data(idx_op_gg + 36);

    auto tr_xxxz_xyyz = pbuffer.data(idx_op_gg + 37);

    auto tr_xxxz_xyzz = pbuffer.data(idx_op_gg + 38);

    auto tr_xxxz_xzzz = pbuffer.data(idx_op_gg + 39);

    auto tr_xxxz_yyyy = pbuffer.data(idx_op_gg + 40);

    auto tr_xxxz_yyyz = pbuffer.data(idx_op_gg + 41);

    auto tr_xxxz_yyzz = pbuffer.data(idx_op_gg + 42);

    auto tr_xxxz_yzzz = pbuffer.data(idx_op_gg + 43);

    auto tr_xxxz_zzzz = pbuffer.data(idx_op_gg + 44);

    auto tr_xxyy_xxxx = pbuffer.data(idx_op_gg + 45);

    auto tr_xxyy_xxxy = pbuffer.data(idx_op_gg + 46);

    auto tr_xxyy_xxxz = pbuffer.data(idx_op_gg + 47);

    auto tr_xxyy_xxyy = pbuffer.data(idx_op_gg + 48);

    auto tr_xxyy_xxyz = pbuffer.data(idx_op_gg + 49);

    auto tr_xxyy_xxzz = pbuffer.data(idx_op_gg + 50);

    auto tr_xxyy_xyyy = pbuffer.data(idx_op_gg + 51);

    auto tr_xxyy_xyyz = pbuffer.data(idx_op_gg + 52);

    auto tr_xxyy_xyzz = pbuffer.data(idx_op_gg + 53);

    auto tr_xxyy_xzzz = pbuffer.data(idx_op_gg + 54);

    auto tr_xxyy_yyyy = pbuffer.data(idx_op_gg + 55);

    auto tr_xxyy_yyyz = pbuffer.data(idx_op_gg + 56);

    auto tr_xxyy_yyzz = pbuffer.data(idx_op_gg + 57);

    auto tr_xxyy_yzzz = pbuffer.data(idx_op_gg + 58);

    auto tr_xxyy_zzzz = pbuffer.data(idx_op_gg + 59);

    auto tr_xxyz_xxxx = pbuffer.data(idx_op_gg + 60);

    auto tr_xxyz_xxxy = pbuffer.data(idx_op_gg + 61);

    auto tr_xxyz_xxxz = pbuffer.data(idx_op_gg + 62);

    auto tr_xxyz_xxyy = pbuffer.data(idx_op_gg + 63);

    auto tr_xxyz_xxyz = pbuffer.data(idx_op_gg + 64);

    auto tr_xxyz_xxzz = pbuffer.data(idx_op_gg + 65);

    auto tr_xxyz_xyyy = pbuffer.data(idx_op_gg + 66);

    auto tr_xxyz_xyyz = pbuffer.data(idx_op_gg + 67);

    auto tr_xxyz_xyzz = pbuffer.data(idx_op_gg + 68);

    auto tr_xxyz_xzzz = pbuffer.data(idx_op_gg + 69);

    auto tr_xxyz_yyyy = pbuffer.data(idx_op_gg + 70);

    auto tr_xxyz_yyyz = pbuffer.data(idx_op_gg + 71);

    auto tr_xxyz_yyzz = pbuffer.data(idx_op_gg + 72);

    auto tr_xxyz_yzzz = pbuffer.data(idx_op_gg + 73);

    auto tr_xxyz_zzzz = pbuffer.data(idx_op_gg + 74);

    auto tr_xxzz_xxxx = pbuffer.data(idx_op_gg + 75);

    auto tr_xxzz_xxxy = pbuffer.data(idx_op_gg + 76);

    auto tr_xxzz_xxxz = pbuffer.data(idx_op_gg + 77);

    auto tr_xxzz_xxyy = pbuffer.data(idx_op_gg + 78);

    auto tr_xxzz_xxyz = pbuffer.data(idx_op_gg + 79);

    auto tr_xxzz_xxzz = pbuffer.data(idx_op_gg + 80);

    auto tr_xxzz_xyyy = pbuffer.data(idx_op_gg + 81);

    auto tr_xxzz_xyyz = pbuffer.data(idx_op_gg + 82);

    auto tr_xxzz_xyzz = pbuffer.data(idx_op_gg + 83);

    auto tr_xxzz_xzzz = pbuffer.data(idx_op_gg + 84);

    auto tr_xxzz_yyyy = pbuffer.data(idx_op_gg + 85);

    auto tr_xxzz_yyyz = pbuffer.data(idx_op_gg + 86);

    auto tr_xxzz_yyzz = pbuffer.data(idx_op_gg + 87);

    auto tr_xxzz_yzzz = pbuffer.data(idx_op_gg + 88);

    auto tr_xxzz_zzzz = pbuffer.data(idx_op_gg + 89);

    auto tr_xyyy_xxxx = pbuffer.data(idx_op_gg + 90);

    auto tr_xyyy_xxxy = pbuffer.data(idx_op_gg + 91);

    auto tr_xyyy_xxxz = pbuffer.data(idx_op_gg + 92);

    auto tr_xyyy_xxyy = pbuffer.data(idx_op_gg + 93);

    auto tr_xyyy_xxyz = pbuffer.data(idx_op_gg + 94);

    auto tr_xyyy_xxzz = pbuffer.data(idx_op_gg + 95);

    auto tr_xyyy_xyyy = pbuffer.data(idx_op_gg + 96);

    auto tr_xyyy_xyyz = pbuffer.data(idx_op_gg + 97);

    auto tr_xyyy_xyzz = pbuffer.data(idx_op_gg + 98);

    auto tr_xyyy_xzzz = pbuffer.data(idx_op_gg + 99);

    auto tr_xyyy_yyyy = pbuffer.data(idx_op_gg + 100);

    auto tr_xyyy_yyyz = pbuffer.data(idx_op_gg + 101);

    auto tr_xyyy_yyzz = pbuffer.data(idx_op_gg + 102);

    auto tr_xyyy_yzzz = pbuffer.data(idx_op_gg + 103);

    auto tr_xyyy_zzzz = pbuffer.data(idx_op_gg + 104);

    auto tr_xyyz_xxxx = pbuffer.data(idx_op_gg + 105);

    auto tr_xyyz_xxxy = pbuffer.data(idx_op_gg + 106);

    auto tr_xyyz_xxxz = pbuffer.data(idx_op_gg + 107);

    auto tr_xyyz_xxyy = pbuffer.data(idx_op_gg + 108);

    auto tr_xyyz_xxyz = pbuffer.data(idx_op_gg + 109);

    auto tr_xyyz_xxzz = pbuffer.data(idx_op_gg + 110);

    auto tr_xyyz_xyyy = pbuffer.data(idx_op_gg + 111);

    auto tr_xyyz_xyyz = pbuffer.data(idx_op_gg + 112);

    auto tr_xyyz_xyzz = pbuffer.data(idx_op_gg + 113);

    auto tr_xyyz_xzzz = pbuffer.data(idx_op_gg + 114);

    auto tr_xyyz_yyyy = pbuffer.data(idx_op_gg + 115);

    auto tr_xyyz_yyyz = pbuffer.data(idx_op_gg + 116);

    auto tr_xyyz_yyzz = pbuffer.data(idx_op_gg + 117);

    auto tr_xyyz_yzzz = pbuffer.data(idx_op_gg + 118);

    auto tr_xyyz_zzzz = pbuffer.data(idx_op_gg + 119);

    auto tr_xyzz_xxxx = pbuffer.data(idx_op_gg + 120);

    auto tr_xyzz_xxxy = pbuffer.data(idx_op_gg + 121);

    auto tr_xyzz_xxxz = pbuffer.data(idx_op_gg + 122);

    auto tr_xyzz_xxyy = pbuffer.data(idx_op_gg + 123);

    auto tr_xyzz_xxyz = pbuffer.data(idx_op_gg + 124);

    auto tr_xyzz_xxzz = pbuffer.data(idx_op_gg + 125);

    auto tr_xyzz_xyyy = pbuffer.data(idx_op_gg + 126);

    auto tr_xyzz_xyyz = pbuffer.data(idx_op_gg + 127);

    auto tr_xyzz_xyzz = pbuffer.data(idx_op_gg + 128);

    auto tr_xyzz_xzzz = pbuffer.data(idx_op_gg + 129);

    auto tr_xyzz_yyyy = pbuffer.data(idx_op_gg + 130);

    auto tr_xyzz_yyyz = pbuffer.data(idx_op_gg + 131);

    auto tr_xyzz_yyzz = pbuffer.data(idx_op_gg + 132);

    auto tr_xyzz_yzzz = pbuffer.data(idx_op_gg + 133);

    auto tr_xyzz_zzzz = pbuffer.data(idx_op_gg + 134);

    auto tr_xzzz_xxxx = pbuffer.data(idx_op_gg + 135);

    auto tr_xzzz_xxxy = pbuffer.data(idx_op_gg + 136);

    auto tr_xzzz_xxxz = pbuffer.data(idx_op_gg + 137);

    auto tr_xzzz_xxyy = pbuffer.data(idx_op_gg + 138);

    auto tr_xzzz_xxyz = pbuffer.data(idx_op_gg + 139);

    auto tr_xzzz_xxzz = pbuffer.data(idx_op_gg + 140);

    auto tr_xzzz_xyyy = pbuffer.data(idx_op_gg + 141);

    auto tr_xzzz_xyyz = pbuffer.data(idx_op_gg + 142);

    auto tr_xzzz_xyzz = pbuffer.data(idx_op_gg + 143);

    auto tr_xzzz_xzzz = pbuffer.data(idx_op_gg + 144);

    auto tr_xzzz_yyyy = pbuffer.data(idx_op_gg + 145);

    auto tr_xzzz_yyyz = pbuffer.data(idx_op_gg + 146);

    auto tr_xzzz_yyzz = pbuffer.data(idx_op_gg + 147);

    auto tr_xzzz_yzzz = pbuffer.data(idx_op_gg + 148);

    auto tr_xzzz_zzzz = pbuffer.data(idx_op_gg + 149);

    auto tr_yyyy_xxxx = pbuffer.data(idx_op_gg + 150);

    auto tr_yyyy_xxxy = pbuffer.data(idx_op_gg + 151);

    auto tr_yyyy_xxxz = pbuffer.data(idx_op_gg + 152);

    auto tr_yyyy_xxyy = pbuffer.data(idx_op_gg + 153);

    auto tr_yyyy_xxyz = pbuffer.data(idx_op_gg + 154);

    auto tr_yyyy_xxzz = pbuffer.data(idx_op_gg + 155);

    auto tr_yyyy_xyyy = pbuffer.data(idx_op_gg + 156);

    auto tr_yyyy_xyyz = pbuffer.data(idx_op_gg + 157);

    auto tr_yyyy_xyzz = pbuffer.data(idx_op_gg + 158);

    auto tr_yyyy_xzzz = pbuffer.data(idx_op_gg + 159);

    auto tr_yyyy_yyyy = pbuffer.data(idx_op_gg + 160);

    auto tr_yyyy_yyyz = pbuffer.data(idx_op_gg + 161);

    auto tr_yyyy_yyzz = pbuffer.data(idx_op_gg + 162);

    auto tr_yyyy_yzzz = pbuffer.data(idx_op_gg + 163);

    auto tr_yyyy_zzzz = pbuffer.data(idx_op_gg + 164);

    auto tr_yyyz_xxxx = pbuffer.data(idx_op_gg + 165);

    auto tr_yyyz_xxxy = pbuffer.data(idx_op_gg + 166);

    auto tr_yyyz_xxxz = pbuffer.data(idx_op_gg + 167);

    auto tr_yyyz_xxyy = pbuffer.data(idx_op_gg + 168);

    auto tr_yyyz_xxyz = pbuffer.data(idx_op_gg + 169);

    auto tr_yyyz_xxzz = pbuffer.data(idx_op_gg + 170);

    auto tr_yyyz_xyyy = pbuffer.data(idx_op_gg + 171);

    auto tr_yyyz_xyyz = pbuffer.data(idx_op_gg + 172);

    auto tr_yyyz_xyzz = pbuffer.data(idx_op_gg + 173);

    auto tr_yyyz_xzzz = pbuffer.data(idx_op_gg + 174);

    auto tr_yyyz_yyyy = pbuffer.data(idx_op_gg + 175);

    auto tr_yyyz_yyyz = pbuffer.data(idx_op_gg + 176);

    auto tr_yyyz_yyzz = pbuffer.data(idx_op_gg + 177);

    auto tr_yyyz_yzzz = pbuffer.data(idx_op_gg + 178);

    auto tr_yyyz_zzzz = pbuffer.data(idx_op_gg + 179);

    auto tr_yyzz_xxxx = pbuffer.data(idx_op_gg + 180);

    auto tr_yyzz_xxxy = pbuffer.data(idx_op_gg + 181);

    auto tr_yyzz_xxxz = pbuffer.data(idx_op_gg + 182);

    auto tr_yyzz_xxyy = pbuffer.data(idx_op_gg + 183);

    auto tr_yyzz_xxyz = pbuffer.data(idx_op_gg + 184);

    auto tr_yyzz_xxzz = pbuffer.data(idx_op_gg + 185);

    auto tr_yyzz_xyyy = pbuffer.data(idx_op_gg + 186);

    auto tr_yyzz_xyyz = pbuffer.data(idx_op_gg + 187);

    auto tr_yyzz_xyzz = pbuffer.data(idx_op_gg + 188);

    auto tr_yyzz_xzzz = pbuffer.data(idx_op_gg + 189);

    auto tr_yyzz_yyyy = pbuffer.data(idx_op_gg + 190);

    auto tr_yyzz_yyyz = pbuffer.data(idx_op_gg + 191);

    auto tr_yyzz_yyzz = pbuffer.data(idx_op_gg + 192);

    auto tr_yyzz_yzzz = pbuffer.data(idx_op_gg + 193);

    auto tr_yyzz_zzzz = pbuffer.data(idx_op_gg + 194);

    auto tr_yzzz_xxxx = pbuffer.data(idx_op_gg + 195);

    auto tr_yzzz_xxxy = pbuffer.data(idx_op_gg + 196);

    auto tr_yzzz_xxxz = pbuffer.data(idx_op_gg + 197);

    auto tr_yzzz_xxyy = pbuffer.data(idx_op_gg + 198);

    auto tr_yzzz_xxyz = pbuffer.data(idx_op_gg + 199);

    auto tr_yzzz_xxzz = pbuffer.data(idx_op_gg + 200);

    auto tr_yzzz_xyyy = pbuffer.data(idx_op_gg + 201);

    auto tr_yzzz_xyyz = pbuffer.data(idx_op_gg + 202);

    auto tr_yzzz_xyzz = pbuffer.data(idx_op_gg + 203);

    auto tr_yzzz_xzzz = pbuffer.data(idx_op_gg + 204);

    auto tr_yzzz_yyyy = pbuffer.data(idx_op_gg + 205);

    auto tr_yzzz_yyyz = pbuffer.data(idx_op_gg + 206);

    auto tr_yzzz_yyzz = pbuffer.data(idx_op_gg + 207);

    auto tr_yzzz_yzzz = pbuffer.data(idx_op_gg + 208);

    auto tr_yzzz_zzzz = pbuffer.data(idx_op_gg + 209);

    auto tr_zzzz_xxxx = pbuffer.data(idx_op_gg + 210);

    auto tr_zzzz_xxxy = pbuffer.data(idx_op_gg + 211);

    auto tr_zzzz_xxxz = pbuffer.data(idx_op_gg + 212);

    auto tr_zzzz_xxyy = pbuffer.data(idx_op_gg + 213);

    auto tr_zzzz_xxyz = pbuffer.data(idx_op_gg + 214);

    auto tr_zzzz_xxzz = pbuffer.data(idx_op_gg + 215);

    auto tr_zzzz_xyyy = pbuffer.data(idx_op_gg + 216);

    auto tr_zzzz_xyyz = pbuffer.data(idx_op_gg + 217);

    auto tr_zzzz_xyzz = pbuffer.data(idx_op_gg + 218);

    auto tr_zzzz_xzzz = pbuffer.data(idx_op_gg + 219);

    auto tr_zzzz_yyyy = pbuffer.data(idx_op_gg + 220);

    auto tr_zzzz_yyyz = pbuffer.data(idx_op_gg + 221);

    auto tr_zzzz_yyzz = pbuffer.data(idx_op_gg + 222);

    auto tr_zzzz_yzzz = pbuffer.data(idx_op_gg + 223);

    auto tr_zzzz_zzzz = pbuffer.data(idx_op_gg + 224);

    // Set up components of auxiliary buffer : HF

    auto tr_xxxxx_xxx = pbuffer.data(idx_op_hf);

    auto tr_xxxxx_xxy = pbuffer.data(idx_op_hf + 1);

    auto tr_xxxxx_xxz = pbuffer.data(idx_op_hf + 2);

    auto tr_xxxxx_xyy = pbuffer.data(idx_op_hf + 3);

    auto tr_xxxxx_xyz = pbuffer.data(idx_op_hf + 4);

    auto tr_xxxxx_xzz = pbuffer.data(idx_op_hf + 5);

    auto tr_xxxxx_yyy = pbuffer.data(idx_op_hf + 6);

    auto tr_xxxxx_yyz = pbuffer.data(idx_op_hf + 7);

    auto tr_xxxxx_yzz = pbuffer.data(idx_op_hf + 8);

    auto tr_xxxxx_zzz = pbuffer.data(idx_op_hf + 9);

    auto tr_xxxxy_xxx = pbuffer.data(idx_op_hf + 10);

    auto tr_xxxxy_xxy = pbuffer.data(idx_op_hf + 11);

    auto tr_xxxxy_xxz = pbuffer.data(idx_op_hf + 12);

    auto tr_xxxxy_xyy = pbuffer.data(idx_op_hf + 13);

    auto tr_xxxxy_xyz = pbuffer.data(idx_op_hf + 14);

    auto tr_xxxxy_xzz = pbuffer.data(idx_op_hf + 15);

    auto tr_xxxxy_yyy = pbuffer.data(idx_op_hf + 16);

    auto tr_xxxxy_yyz = pbuffer.data(idx_op_hf + 17);

    auto tr_xxxxy_yzz = pbuffer.data(idx_op_hf + 18);

    auto tr_xxxxy_zzz = pbuffer.data(idx_op_hf + 19);

    auto tr_xxxxz_xxx = pbuffer.data(idx_op_hf + 20);

    auto tr_xxxxz_xxy = pbuffer.data(idx_op_hf + 21);

    auto tr_xxxxz_xxz = pbuffer.data(idx_op_hf + 22);

    auto tr_xxxxz_xyy = pbuffer.data(idx_op_hf + 23);

    auto tr_xxxxz_xyz = pbuffer.data(idx_op_hf + 24);

    auto tr_xxxxz_xzz = pbuffer.data(idx_op_hf + 25);

    auto tr_xxxxz_yyy = pbuffer.data(idx_op_hf + 26);

    auto tr_xxxxz_yyz = pbuffer.data(idx_op_hf + 27);

    auto tr_xxxxz_yzz = pbuffer.data(idx_op_hf + 28);

    auto tr_xxxxz_zzz = pbuffer.data(idx_op_hf + 29);

    auto tr_xxxyy_xxx = pbuffer.data(idx_op_hf + 30);

    auto tr_xxxyy_xxy = pbuffer.data(idx_op_hf + 31);

    auto tr_xxxyy_xxz = pbuffer.data(idx_op_hf + 32);

    auto tr_xxxyy_xyy = pbuffer.data(idx_op_hf + 33);

    auto tr_xxxyy_xyz = pbuffer.data(idx_op_hf + 34);

    auto tr_xxxyy_xzz = pbuffer.data(idx_op_hf + 35);

    auto tr_xxxyy_yyy = pbuffer.data(idx_op_hf + 36);

    auto tr_xxxyy_yyz = pbuffer.data(idx_op_hf + 37);

    auto tr_xxxyy_yzz = pbuffer.data(idx_op_hf + 38);

    auto tr_xxxyy_zzz = pbuffer.data(idx_op_hf + 39);

    auto tr_xxxyz_xxx = pbuffer.data(idx_op_hf + 40);

    auto tr_xxxyz_xxy = pbuffer.data(idx_op_hf + 41);

    auto tr_xxxyz_xxz = pbuffer.data(idx_op_hf + 42);

    auto tr_xxxyz_xyy = pbuffer.data(idx_op_hf + 43);

    auto tr_xxxyz_xyz = pbuffer.data(idx_op_hf + 44);

    auto tr_xxxyz_xzz = pbuffer.data(idx_op_hf + 45);

    auto tr_xxxyz_yyy = pbuffer.data(idx_op_hf + 46);

    auto tr_xxxyz_yyz = pbuffer.data(idx_op_hf + 47);

    auto tr_xxxyz_yzz = pbuffer.data(idx_op_hf + 48);

    auto tr_xxxyz_zzz = pbuffer.data(idx_op_hf + 49);

    auto tr_xxxzz_xxx = pbuffer.data(idx_op_hf + 50);

    auto tr_xxxzz_xxy = pbuffer.data(idx_op_hf + 51);

    auto tr_xxxzz_xxz = pbuffer.data(idx_op_hf + 52);

    auto tr_xxxzz_xyy = pbuffer.data(idx_op_hf + 53);

    auto tr_xxxzz_xyz = pbuffer.data(idx_op_hf + 54);

    auto tr_xxxzz_xzz = pbuffer.data(idx_op_hf + 55);

    auto tr_xxxzz_yyy = pbuffer.data(idx_op_hf + 56);

    auto tr_xxxzz_yyz = pbuffer.data(idx_op_hf + 57);

    auto tr_xxxzz_yzz = pbuffer.data(idx_op_hf + 58);

    auto tr_xxxzz_zzz = pbuffer.data(idx_op_hf + 59);

    auto tr_xxyyy_xxx = pbuffer.data(idx_op_hf + 60);

    auto tr_xxyyy_xxy = pbuffer.data(idx_op_hf + 61);

    auto tr_xxyyy_xxz = pbuffer.data(idx_op_hf + 62);

    auto tr_xxyyy_xyy = pbuffer.data(idx_op_hf + 63);

    auto tr_xxyyy_xyz = pbuffer.data(idx_op_hf + 64);

    auto tr_xxyyy_xzz = pbuffer.data(idx_op_hf + 65);

    auto tr_xxyyy_yyy = pbuffer.data(idx_op_hf + 66);

    auto tr_xxyyy_yyz = pbuffer.data(idx_op_hf + 67);

    auto tr_xxyyy_yzz = pbuffer.data(idx_op_hf + 68);

    auto tr_xxyyy_zzz = pbuffer.data(idx_op_hf + 69);

    auto tr_xxyyz_xxx = pbuffer.data(idx_op_hf + 70);

    auto tr_xxyyz_xxy = pbuffer.data(idx_op_hf + 71);

    auto tr_xxyyz_xxz = pbuffer.data(idx_op_hf + 72);

    auto tr_xxyyz_xyy = pbuffer.data(idx_op_hf + 73);

    auto tr_xxyyz_xyz = pbuffer.data(idx_op_hf + 74);

    auto tr_xxyyz_xzz = pbuffer.data(idx_op_hf + 75);

    auto tr_xxyyz_yyy = pbuffer.data(idx_op_hf + 76);

    auto tr_xxyyz_yyz = pbuffer.data(idx_op_hf + 77);

    auto tr_xxyyz_yzz = pbuffer.data(idx_op_hf + 78);

    auto tr_xxyyz_zzz = pbuffer.data(idx_op_hf + 79);

    auto tr_xxyzz_xxx = pbuffer.data(idx_op_hf + 80);

    auto tr_xxyzz_xxy = pbuffer.data(idx_op_hf + 81);

    auto tr_xxyzz_xxz = pbuffer.data(idx_op_hf + 82);

    auto tr_xxyzz_xyy = pbuffer.data(idx_op_hf + 83);

    auto tr_xxyzz_xyz = pbuffer.data(idx_op_hf + 84);

    auto tr_xxyzz_xzz = pbuffer.data(idx_op_hf + 85);

    auto tr_xxyzz_yyy = pbuffer.data(idx_op_hf + 86);

    auto tr_xxyzz_yyz = pbuffer.data(idx_op_hf + 87);

    auto tr_xxyzz_yzz = pbuffer.data(idx_op_hf + 88);

    auto tr_xxyzz_zzz = pbuffer.data(idx_op_hf + 89);

    auto tr_xxzzz_xxx = pbuffer.data(idx_op_hf + 90);

    auto tr_xxzzz_xxy = pbuffer.data(idx_op_hf + 91);

    auto tr_xxzzz_xxz = pbuffer.data(idx_op_hf + 92);

    auto tr_xxzzz_xyy = pbuffer.data(idx_op_hf + 93);

    auto tr_xxzzz_xyz = pbuffer.data(idx_op_hf + 94);

    auto tr_xxzzz_xzz = pbuffer.data(idx_op_hf + 95);

    auto tr_xxzzz_yyy = pbuffer.data(idx_op_hf + 96);

    auto tr_xxzzz_yyz = pbuffer.data(idx_op_hf + 97);

    auto tr_xxzzz_yzz = pbuffer.data(idx_op_hf + 98);

    auto tr_xxzzz_zzz = pbuffer.data(idx_op_hf + 99);

    auto tr_xyyyy_xxx = pbuffer.data(idx_op_hf + 100);

    auto tr_xyyyy_xxy = pbuffer.data(idx_op_hf + 101);

    auto tr_xyyyy_xxz = pbuffer.data(idx_op_hf + 102);

    auto tr_xyyyy_xyy = pbuffer.data(idx_op_hf + 103);

    auto tr_xyyyy_xyz = pbuffer.data(idx_op_hf + 104);

    auto tr_xyyyy_xzz = pbuffer.data(idx_op_hf + 105);

    auto tr_xyyyy_yyy = pbuffer.data(idx_op_hf + 106);

    auto tr_xyyyy_yyz = pbuffer.data(idx_op_hf + 107);

    auto tr_xyyyy_yzz = pbuffer.data(idx_op_hf + 108);

    auto tr_xyyyy_zzz = pbuffer.data(idx_op_hf + 109);

    auto tr_xyyyz_xxx = pbuffer.data(idx_op_hf + 110);

    auto tr_xyyyz_xxy = pbuffer.data(idx_op_hf + 111);

    auto tr_xyyyz_xxz = pbuffer.data(idx_op_hf + 112);

    auto tr_xyyyz_xyy = pbuffer.data(idx_op_hf + 113);

    auto tr_xyyyz_xyz = pbuffer.data(idx_op_hf + 114);

    auto tr_xyyyz_xzz = pbuffer.data(idx_op_hf + 115);

    auto tr_xyyyz_yyy = pbuffer.data(idx_op_hf + 116);

    auto tr_xyyyz_yyz = pbuffer.data(idx_op_hf + 117);

    auto tr_xyyyz_yzz = pbuffer.data(idx_op_hf + 118);

    auto tr_xyyyz_zzz = pbuffer.data(idx_op_hf + 119);

    auto tr_xyyzz_xxx = pbuffer.data(idx_op_hf + 120);

    auto tr_xyyzz_xxy = pbuffer.data(idx_op_hf + 121);

    auto tr_xyyzz_xxz = pbuffer.data(idx_op_hf + 122);

    auto tr_xyyzz_xyy = pbuffer.data(idx_op_hf + 123);

    auto tr_xyyzz_xyz = pbuffer.data(idx_op_hf + 124);

    auto tr_xyyzz_xzz = pbuffer.data(idx_op_hf + 125);

    auto tr_xyyzz_yyy = pbuffer.data(idx_op_hf + 126);

    auto tr_xyyzz_yyz = pbuffer.data(idx_op_hf + 127);

    auto tr_xyyzz_yzz = pbuffer.data(idx_op_hf + 128);

    auto tr_xyyzz_zzz = pbuffer.data(idx_op_hf + 129);

    auto tr_xyzzz_xxx = pbuffer.data(idx_op_hf + 130);

    auto tr_xyzzz_xxy = pbuffer.data(idx_op_hf + 131);

    auto tr_xyzzz_xxz = pbuffer.data(idx_op_hf + 132);

    auto tr_xyzzz_xyy = pbuffer.data(idx_op_hf + 133);

    auto tr_xyzzz_xyz = pbuffer.data(idx_op_hf + 134);

    auto tr_xyzzz_xzz = pbuffer.data(idx_op_hf + 135);

    auto tr_xyzzz_yyy = pbuffer.data(idx_op_hf + 136);

    auto tr_xyzzz_yyz = pbuffer.data(idx_op_hf + 137);

    auto tr_xyzzz_yzz = pbuffer.data(idx_op_hf + 138);

    auto tr_xyzzz_zzz = pbuffer.data(idx_op_hf + 139);

    auto tr_xzzzz_xxx = pbuffer.data(idx_op_hf + 140);

    auto tr_xzzzz_xxy = pbuffer.data(idx_op_hf + 141);

    auto tr_xzzzz_xxz = pbuffer.data(idx_op_hf + 142);

    auto tr_xzzzz_xyy = pbuffer.data(idx_op_hf + 143);

    auto tr_xzzzz_xyz = pbuffer.data(idx_op_hf + 144);

    auto tr_xzzzz_xzz = pbuffer.data(idx_op_hf + 145);

    auto tr_xzzzz_yyy = pbuffer.data(idx_op_hf + 146);

    auto tr_xzzzz_yyz = pbuffer.data(idx_op_hf + 147);

    auto tr_xzzzz_yzz = pbuffer.data(idx_op_hf + 148);

    auto tr_xzzzz_zzz = pbuffer.data(idx_op_hf + 149);

    auto tr_yyyyy_xxx = pbuffer.data(idx_op_hf + 150);

    auto tr_yyyyy_xxy = pbuffer.data(idx_op_hf + 151);

    auto tr_yyyyy_xxz = pbuffer.data(idx_op_hf + 152);

    auto tr_yyyyy_xyy = pbuffer.data(idx_op_hf + 153);

    auto tr_yyyyy_xyz = pbuffer.data(idx_op_hf + 154);

    auto tr_yyyyy_xzz = pbuffer.data(idx_op_hf + 155);

    auto tr_yyyyy_yyy = pbuffer.data(idx_op_hf + 156);

    auto tr_yyyyy_yyz = pbuffer.data(idx_op_hf + 157);

    auto tr_yyyyy_yzz = pbuffer.data(idx_op_hf + 158);

    auto tr_yyyyy_zzz = pbuffer.data(idx_op_hf + 159);

    auto tr_yyyyz_xxx = pbuffer.data(idx_op_hf + 160);

    auto tr_yyyyz_xxy = pbuffer.data(idx_op_hf + 161);

    auto tr_yyyyz_xxz = pbuffer.data(idx_op_hf + 162);

    auto tr_yyyyz_xyy = pbuffer.data(idx_op_hf + 163);

    auto tr_yyyyz_xyz = pbuffer.data(idx_op_hf + 164);

    auto tr_yyyyz_xzz = pbuffer.data(idx_op_hf + 165);

    auto tr_yyyyz_yyy = pbuffer.data(idx_op_hf + 166);

    auto tr_yyyyz_yyz = pbuffer.data(idx_op_hf + 167);

    auto tr_yyyyz_yzz = pbuffer.data(idx_op_hf + 168);

    auto tr_yyyyz_zzz = pbuffer.data(idx_op_hf + 169);

    auto tr_yyyzz_xxx = pbuffer.data(idx_op_hf + 170);

    auto tr_yyyzz_xxy = pbuffer.data(idx_op_hf + 171);

    auto tr_yyyzz_xxz = pbuffer.data(idx_op_hf + 172);

    auto tr_yyyzz_xyy = pbuffer.data(idx_op_hf + 173);

    auto tr_yyyzz_xyz = pbuffer.data(idx_op_hf + 174);

    auto tr_yyyzz_xzz = pbuffer.data(idx_op_hf + 175);

    auto tr_yyyzz_yyy = pbuffer.data(idx_op_hf + 176);

    auto tr_yyyzz_yyz = pbuffer.data(idx_op_hf + 177);

    auto tr_yyyzz_yzz = pbuffer.data(idx_op_hf + 178);

    auto tr_yyyzz_zzz = pbuffer.data(idx_op_hf + 179);

    auto tr_yyzzz_xxx = pbuffer.data(idx_op_hf + 180);

    auto tr_yyzzz_xxy = pbuffer.data(idx_op_hf + 181);

    auto tr_yyzzz_xxz = pbuffer.data(idx_op_hf + 182);

    auto tr_yyzzz_xyy = pbuffer.data(idx_op_hf + 183);

    auto tr_yyzzz_xyz = pbuffer.data(idx_op_hf + 184);

    auto tr_yyzzz_xzz = pbuffer.data(idx_op_hf + 185);

    auto tr_yyzzz_yyy = pbuffer.data(idx_op_hf + 186);

    auto tr_yyzzz_yyz = pbuffer.data(idx_op_hf + 187);

    auto tr_yyzzz_yzz = pbuffer.data(idx_op_hf + 188);

    auto tr_yyzzz_zzz = pbuffer.data(idx_op_hf + 189);

    auto tr_yzzzz_xxx = pbuffer.data(idx_op_hf + 190);

    auto tr_yzzzz_xxy = pbuffer.data(idx_op_hf + 191);

    auto tr_yzzzz_xxz = pbuffer.data(idx_op_hf + 192);

    auto tr_yzzzz_xyy = pbuffer.data(idx_op_hf + 193);

    auto tr_yzzzz_xyz = pbuffer.data(idx_op_hf + 194);

    auto tr_yzzzz_xzz = pbuffer.data(idx_op_hf + 195);

    auto tr_yzzzz_yyy = pbuffer.data(idx_op_hf + 196);

    auto tr_yzzzz_yyz = pbuffer.data(idx_op_hf + 197);

    auto tr_yzzzz_yzz = pbuffer.data(idx_op_hf + 198);

    auto tr_yzzzz_zzz = pbuffer.data(idx_op_hf + 199);

    auto tr_zzzzz_xxx = pbuffer.data(idx_op_hf + 200);

    auto tr_zzzzz_xxy = pbuffer.data(idx_op_hf + 201);

    auto tr_zzzzz_xxz = pbuffer.data(idx_op_hf + 202);

    auto tr_zzzzz_xyy = pbuffer.data(idx_op_hf + 203);

    auto tr_zzzzz_xyz = pbuffer.data(idx_op_hf + 204);

    auto tr_zzzzz_xzz = pbuffer.data(idx_op_hf + 205);

    auto tr_zzzzz_yyy = pbuffer.data(idx_op_hf + 206);

    auto tr_zzzzz_yyz = pbuffer.data(idx_op_hf + 207);

    auto tr_zzzzz_yzz = pbuffer.data(idx_op_hf + 208);

    auto tr_zzzzz_zzz = pbuffer.data(idx_op_hf + 209);

    // Set up 0-10 components of targeted buffer : FF

    auto tr_0_0_xx_xxx_xxx = pbuffer.data(idx_op_geom_020_ff);

    auto tr_0_0_xx_xxx_xxy = pbuffer.data(idx_op_geom_020_ff + 1);

    auto tr_0_0_xx_xxx_xxz = pbuffer.data(idx_op_geom_020_ff + 2);

    auto tr_0_0_xx_xxx_xyy = pbuffer.data(idx_op_geom_020_ff + 3);

    auto tr_0_0_xx_xxx_xyz = pbuffer.data(idx_op_geom_020_ff + 4);

    auto tr_0_0_xx_xxx_xzz = pbuffer.data(idx_op_geom_020_ff + 5);

    auto tr_0_0_xx_xxx_yyy = pbuffer.data(idx_op_geom_020_ff + 6);

    auto tr_0_0_xx_xxx_yyz = pbuffer.data(idx_op_geom_020_ff + 7);

    auto tr_0_0_xx_xxx_yzz = pbuffer.data(idx_op_geom_020_ff + 8);

    auto tr_0_0_xx_xxx_zzz = pbuffer.data(idx_op_geom_020_ff + 9);

    #pragma omp simd aligned(tr_0_0_xx_xxx_xxx, tr_0_0_xx_xxx_xxy, tr_0_0_xx_xxx_xxz, tr_0_0_xx_xxx_xyy, tr_0_0_xx_xxx_xyz, tr_0_0_xx_xxx_xzz, tr_0_0_xx_xxx_yyy, tr_0_0_xx_xxx_yyz, tr_0_0_xx_xxx_yzz, tr_0_0_xx_xxx_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xx, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxxxx, tr_xxx_xxxxy, tr_xxx_xxxxz, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxx_xxx, tr_xxxxx_xxy, tr_xxxxx_xxz, tr_xxxxx_xyy, tr_xxxxx_xyz, tr_xxxxx_xzz, tr_xxxxx_yyy, tr_xxxxx_yyz, tr_xxxxx_yzz, tr_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxx_xxx[i] = 6.0 * tr_x_xxx[i] + 18.0 * tr_xx_xx[i] - 12.0 * tr_xx_xxxx[i] * tke_0 + 6.0 * tr_xxx_x[i] - 14.0 * tr_xxx_xxx[i] * tbe_0 - 14.0 * tr_xxx_xxx[i] * tke_0 + 4.0 * tr_xxx_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxxx_xx[i] * tbe_0 + 8.0 * tr_xxxx_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_xxy[i] = 6.0 * tr_x_xxy[i] + 12.0 * tr_xx_xy[i] - 12.0 * tr_xx_xxxy[i] * tke_0 + 2.0 * tr_xxx_y[i] - 14.0 * tr_xxx_xxy[i] * tbe_0 - 10.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxxx_xy[i] * tbe_0 + 8.0 * tr_xxxx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_xxz[i] = 6.0 * tr_x_xxz[i] + 12.0 * tr_xx_xz[i] - 12.0 * tr_xx_xxxz[i] * tke_0 + 2.0 * tr_xxx_z[i] - 14.0 * tr_xxx_xxz[i] * tbe_0 - 10.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxxx_xz[i] * tbe_0 + 8.0 * tr_xxxx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_xyy[i] = 6.0 * tr_x_xyy[i] + 6.0 * tr_xx_yy[i] - 12.0 * tr_xx_xxyy[i] * tke_0 - 14.0 * tr_xxx_xyy[i] * tbe_0 - 6.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxx_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxx_yy[i] * tbe_0 + 8.0 * tr_xxxx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_xyz[i] = 6.0 * tr_x_xyz[i] + 6.0 * tr_xx_yz[i] - 12.0 * tr_xx_xxyz[i] * tke_0 - 14.0 * tr_xxx_xyz[i] * tbe_0 - 6.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxx_yz[i] * tbe_0 + 8.0 * tr_xxxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_xzz[i] = 6.0 * tr_x_xzz[i] + 6.0 * tr_xx_zz[i] - 12.0 * tr_xx_xxzz[i] * tke_0 - 14.0 * tr_xxx_xzz[i] * tbe_0 - 6.0 * tr_xxx_xzz[i] * tke_0 + 4.0 * tr_xxx_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxx_zz[i] * tbe_0 + 8.0 * tr_xxxx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_yyy[i] = 6.0 * tr_x_yyy[i] - 12.0 * tr_xx_xyyy[i] * tke_0 - 14.0 * tr_xxx_yyy[i] * tbe_0 - 2.0 * tr_xxx_yyy[i] * tke_0 + 4.0 * tr_xxx_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_yyz[i] = 6.0 * tr_x_yyz[i] - 12.0 * tr_xx_xyyz[i] * tke_0 - 14.0 * tr_xxx_yyz[i] * tbe_0 - 2.0 * tr_xxx_yyz[i] * tke_0 + 4.0 * tr_xxx_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_yzz[i] = 6.0 * tr_x_yzz[i] - 12.0 * tr_xx_xyzz[i] * tke_0 - 14.0 * tr_xxx_yzz[i] * tbe_0 - 2.0 * tr_xxx_yzz[i] * tke_0 + 4.0 * tr_xxx_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxx_zzz[i] = 6.0 * tr_x_zzz[i] - 12.0 * tr_xx_xzzz[i] * tke_0 - 14.0 * tr_xxx_zzz[i] * tbe_0 - 2.0 * tr_xxx_zzz[i] * tke_0 + 4.0 * tr_xxx_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 10-20 components of targeted buffer : FF

    auto tr_0_0_xx_xxy_xxx = pbuffer.data(idx_op_geom_020_ff + 10);

    auto tr_0_0_xx_xxy_xxy = pbuffer.data(idx_op_geom_020_ff + 11);

    auto tr_0_0_xx_xxy_xxz = pbuffer.data(idx_op_geom_020_ff + 12);

    auto tr_0_0_xx_xxy_xyy = pbuffer.data(idx_op_geom_020_ff + 13);

    auto tr_0_0_xx_xxy_xyz = pbuffer.data(idx_op_geom_020_ff + 14);

    auto tr_0_0_xx_xxy_xzz = pbuffer.data(idx_op_geom_020_ff + 15);

    auto tr_0_0_xx_xxy_yyy = pbuffer.data(idx_op_geom_020_ff + 16);

    auto tr_0_0_xx_xxy_yyz = pbuffer.data(idx_op_geom_020_ff + 17);

    auto tr_0_0_xx_xxy_yzz = pbuffer.data(idx_op_geom_020_ff + 18);

    auto tr_0_0_xx_xxy_zzz = pbuffer.data(idx_op_geom_020_ff + 19);

    #pragma omp simd aligned(tr_0_0_xx_xxy_xxx, tr_0_0_xx_xxy_xxy, tr_0_0_xx_xxy_xxz, tr_0_0_xx_xxy_xyy, tr_0_0_xx_xxy_xyz, tr_0_0_xx_xxy_xzz, tr_0_0_xx_xxy_yyy, tr_0_0_xx_xxy_yyz, tr_0_0_xx_xxy_yzz, tr_0_0_xx_xxy_zzz, tr_xxxxy_xxx, tr_xxxxy_xxy, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_zzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxy_xxx[i] = 2.0 * tr_y_xxx[i] + 12.0 * tr_xy_xx[i] - 8.0 * tr_xy_xxxx[i] * tke_0 + 6.0 * tr_xxy_x[i] - 10.0 * tr_xxy_xxx[i] * tbe_0 - 14.0 * tr_xxy_xxx[i] * tke_0 + 4.0 * tr_xxy_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxxy_xx[i] * tbe_0 + 8.0 * tr_xxxy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_xxy[i] = 2.0 * tr_y_xxy[i] + 8.0 * tr_xy_xy[i] - 8.0 * tr_xy_xxxy[i] * tke_0 + 2.0 * tr_xxy_y[i] - 10.0 * tr_xxy_xxy[i] * tbe_0 - 10.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxxy_xy[i] * tbe_0 + 8.0 * tr_xxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_xxz[i] = 2.0 * tr_y_xxz[i] + 8.0 * tr_xy_xz[i] - 8.0 * tr_xy_xxxz[i] * tke_0 + 2.0 * tr_xxy_z[i] - 10.0 * tr_xxy_xxz[i] * tbe_0 - 10.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxxy_xz[i] * tbe_0 + 8.0 * tr_xxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_xyy[i] = 2.0 * tr_y_xyy[i] + 4.0 * tr_xy_yy[i] - 8.0 * tr_xy_xxyy[i] * tke_0 - 10.0 * tr_xxy_xyy[i] * tbe_0 - 6.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_yy[i] * tbe_0 + 8.0 * tr_xxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_xyz[i] = 2.0 * tr_y_xyz[i] + 4.0 * tr_xy_yz[i] - 8.0 * tr_xy_xxyz[i] * tke_0 - 10.0 * tr_xxy_xyz[i] * tbe_0 - 6.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_yz[i] * tbe_0 + 8.0 * tr_xxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_xzz[i] = 2.0 * tr_y_xzz[i] + 4.0 * tr_xy_zz[i] - 8.0 * tr_xy_xxzz[i] * tke_0 - 10.0 * tr_xxy_xzz[i] * tbe_0 - 6.0 * tr_xxy_xzz[i] * tke_0 + 4.0 * tr_xxy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_zz[i] * tbe_0 + 8.0 * tr_xxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_yyy[i] = 2.0 * tr_y_yyy[i] - 8.0 * tr_xy_xyyy[i] * tke_0 - 10.0 * tr_xxy_yyy[i] * tbe_0 - 2.0 * tr_xxy_yyy[i] * tke_0 + 4.0 * tr_xxy_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_yyz[i] = 2.0 * tr_y_yyz[i] - 8.0 * tr_xy_xyyz[i] * tke_0 - 10.0 * tr_xxy_yyz[i] * tbe_0 - 2.0 * tr_xxy_yyz[i] * tke_0 + 4.0 * tr_xxy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_yzz[i] = 2.0 * tr_y_yzz[i] - 8.0 * tr_xy_xyzz[i] * tke_0 - 10.0 * tr_xxy_yzz[i] * tbe_0 - 2.0 * tr_xxy_yzz[i] * tke_0 + 4.0 * tr_xxy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxy_zzz[i] = 2.0 * tr_y_zzz[i] - 8.0 * tr_xy_xzzz[i] * tke_0 - 10.0 * tr_xxy_zzz[i] * tbe_0 - 2.0 * tr_xxy_zzz[i] * tke_0 + 4.0 * tr_xxy_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 20-30 components of targeted buffer : FF

    auto tr_0_0_xx_xxz_xxx = pbuffer.data(idx_op_geom_020_ff + 20);

    auto tr_0_0_xx_xxz_xxy = pbuffer.data(idx_op_geom_020_ff + 21);

    auto tr_0_0_xx_xxz_xxz = pbuffer.data(idx_op_geom_020_ff + 22);

    auto tr_0_0_xx_xxz_xyy = pbuffer.data(idx_op_geom_020_ff + 23);

    auto tr_0_0_xx_xxz_xyz = pbuffer.data(idx_op_geom_020_ff + 24);

    auto tr_0_0_xx_xxz_xzz = pbuffer.data(idx_op_geom_020_ff + 25);

    auto tr_0_0_xx_xxz_yyy = pbuffer.data(idx_op_geom_020_ff + 26);

    auto tr_0_0_xx_xxz_yyz = pbuffer.data(idx_op_geom_020_ff + 27);

    auto tr_0_0_xx_xxz_yzz = pbuffer.data(idx_op_geom_020_ff + 28);

    auto tr_0_0_xx_xxz_zzz = pbuffer.data(idx_op_geom_020_ff + 29);

    #pragma omp simd aligned(tr_0_0_xx_xxz_xxx, tr_0_0_xx_xxz_xxy, tr_0_0_xx_xxz_xxz, tr_0_0_xx_xxz_xyy, tr_0_0_xx_xxz_xyz, tr_0_0_xx_xxz_xzz, tr_0_0_xx_xxz_yyy, tr_0_0_xx_xxz_yyz, tr_0_0_xx_xxz_yzz, tr_0_0_xx_xxz_zzz, tr_xxxxz_xxx, tr_xxxxz_xxy, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_zzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxz_xxx[i] = 2.0 * tr_z_xxx[i] + 12.0 * tr_xz_xx[i] - 8.0 * tr_xz_xxxx[i] * tke_0 + 6.0 * tr_xxz_x[i] - 10.0 * tr_xxz_xxx[i] * tbe_0 - 14.0 * tr_xxz_xxx[i] * tke_0 + 4.0 * tr_xxz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxxz_xx[i] * tbe_0 + 8.0 * tr_xxxz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_xxy[i] = 2.0 * tr_z_xxy[i] + 8.0 * tr_xz_xy[i] - 8.0 * tr_xz_xxxy[i] * tke_0 + 2.0 * tr_xxz_y[i] - 10.0 * tr_xxz_xxy[i] * tbe_0 - 10.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxxz_xy[i] * tbe_0 + 8.0 * tr_xxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_xxz[i] = 2.0 * tr_z_xxz[i] + 8.0 * tr_xz_xz[i] - 8.0 * tr_xz_xxxz[i] * tke_0 + 2.0 * tr_xxz_z[i] - 10.0 * tr_xxz_xxz[i] * tbe_0 - 10.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxxz_xz[i] * tbe_0 + 8.0 * tr_xxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_xyy[i] = 2.0 * tr_z_xyy[i] + 4.0 * tr_xz_yy[i] - 8.0 * tr_xz_xxyy[i] * tke_0 - 10.0 * tr_xxz_xyy[i] * tbe_0 - 6.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_yy[i] * tbe_0 + 8.0 * tr_xxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_xyz[i] = 2.0 * tr_z_xyz[i] + 4.0 * tr_xz_yz[i] - 8.0 * tr_xz_xxyz[i] * tke_0 - 10.0 * tr_xxz_xyz[i] * tbe_0 - 6.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_yz[i] * tbe_0 + 8.0 * tr_xxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_xzz[i] = 2.0 * tr_z_xzz[i] + 4.0 * tr_xz_zz[i] - 8.0 * tr_xz_xxzz[i] * tke_0 - 10.0 * tr_xxz_xzz[i] * tbe_0 - 6.0 * tr_xxz_xzz[i] * tke_0 + 4.0 * tr_xxz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_zz[i] * tbe_0 + 8.0 * tr_xxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_yyy[i] = 2.0 * tr_z_yyy[i] - 8.0 * tr_xz_xyyy[i] * tke_0 - 10.0 * tr_xxz_yyy[i] * tbe_0 - 2.0 * tr_xxz_yyy[i] * tke_0 + 4.0 * tr_xxz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_yyz[i] = 2.0 * tr_z_yyz[i] - 8.0 * tr_xz_xyyz[i] * tke_0 - 10.0 * tr_xxz_yyz[i] * tbe_0 - 2.0 * tr_xxz_yyz[i] * tke_0 + 4.0 * tr_xxz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_yzz[i] = 2.0 * tr_z_yzz[i] - 8.0 * tr_xz_xyzz[i] * tke_0 - 10.0 * tr_xxz_yzz[i] * tbe_0 - 2.0 * tr_xxz_yzz[i] * tke_0 + 4.0 * tr_xxz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxz_zzz[i] = 2.0 * tr_z_zzz[i] - 8.0 * tr_xz_xzzz[i] * tke_0 - 10.0 * tr_xxz_zzz[i] * tbe_0 - 2.0 * tr_xxz_zzz[i] * tke_0 + 4.0 * tr_xxz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-40 components of targeted buffer : FF

    auto tr_0_0_xx_xyy_xxx = pbuffer.data(idx_op_geom_020_ff + 30);

    auto tr_0_0_xx_xyy_xxy = pbuffer.data(idx_op_geom_020_ff + 31);

    auto tr_0_0_xx_xyy_xxz = pbuffer.data(idx_op_geom_020_ff + 32);

    auto tr_0_0_xx_xyy_xyy = pbuffer.data(idx_op_geom_020_ff + 33);

    auto tr_0_0_xx_xyy_xyz = pbuffer.data(idx_op_geom_020_ff + 34);

    auto tr_0_0_xx_xyy_xzz = pbuffer.data(idx_op_geom_020_ff + 35);

    auto tr_0_0_xx_xyy_yyy = pbuffer.data(idx_op_geom_020_ff + 36);

    auto tr_0_0_xx_xyy_yyz = pbuffer.data(idx_op_geom_020_ff + 37);

    auto tr_0_0_xx_xyy_yzz = pbuffer.data(idx_op_geom_020_ff + 38);

    auto tr_0_0_xx_xyy_zzz = pbuffer.data(idx_op_geom_020_ff + 39);

    #pragma omp simd aligned(tr_0_0_xx_xyy_xxx, tr_0_0_xx_xyy_xxy, tr_0_0_xx_xyy_xxz, tr_0_0_xx_xyy_xyy, tr_0_0_xx_xyy_xyz, tr_0_0_xx_xyy_xzz, tr_0_0_xx_xyy_yyy, tr_0_0_xx_xyy_yyz, tr_0_0_xx_xyy_yzz, tr_0_0_xx_xyy_zzz, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_yy_xx, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyy_xxx[i] = 6.0 * tr_yy_xx[i] - 4.0 * tr_yy_xxxx[i] * tke_0 + 6.0 * tr_xyy_x[i] - 6.0 * tr_xyy_xxx[i] * tbe_0 - 14.0 * tr_xyy_xxx[i] * tke_0 + 4.0 * tr_xyy_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxyy_xx[i] * tbe_0 + 8.0 * tr_xxyy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_xxy[i] = 4.0 * tr_yy_xy[i] - 4.0 * tr_yy_xxxy[i] * tke_0 + 2.0 * tr_xyy_y[i] - 6.0 * tr_xyy_xxy[i] * tbe_0 - 10.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxyy_xy[i] * tbe_0 + 8.0 * tr_xxyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_xxz[i] = 4.0 * tr_yy_xz[i] - 4.0 * tr_yy_xxxz[i] * tke_0 + 2.0 * tr_xyy_z[i] - 6.0 * tr_xyy_xxz[i] * tbe_0 - 10.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxyy_xz[i] * tbe_0 + 8.0 * tr_xxyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_xyy[i] = 2.0 * tr_yy_yy[i] - 4.0 * tr_yy_xxyy[i] * tke_0 - 6.0 * tr_xyy_xyy[i] * tbe_0 - 6.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_yy[i] * tbe_0 + 8.0 * tr_xxyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_xyz[i] = 2.0 * tr_yy_yz[i] - 4.0 * tr_yy_xxyz[i] * tke_0 - 6.0 * tr_xyy_xyz[i] * tbe_0 - 6.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_yz[i] * tbe_0 + 8.0 * tr_xxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_xzz[i] = 2.0 * tr_yy_zz[i] - 4.0 * tr_yy_xxzz[i] * tke_0 - 6.0 * tr_xyy_xzz[i] * tbe_0 - 6.0 * tr_xyy_xzz[i] * tke_0 + 4.0 * tr_xyy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_zz[i] * tbe_0 + 8.0 * tr_xxyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_yyy[i] = -4.0 * tr_yy_xyyy[i] * tke_0 - 6.0 * tr_xyy_yyy[i] * tbe_0 - 2.0 * tr_xyy_yyy[i] * tke_0 + 4.0 * tr_xyy_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_yyz[i] = -4.0 * tr_yy_xyyz[i] * tke_0 - 6.0 * tr_xyy_yyz[i] * tbe_0 - 2.0 * tr_xyy_yyz[i] * tke_0 + 4.0 * tr_xyy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_yzz[i] = -4.0 * tr_yy_xyzz[i] * tke_0 - 6.0 * tr_xyy_yzz[i] * tbe_0 - 2.0 * tr_xyy_yzz[i] * tke_0 + 4.0 * tr_xyy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyy_zzz[i] = -4.0 * tr_yy_xzzz[i] * tke_0 - 6.0 * tr_xyy_zzz[i] * tbe_0 - 2.0 * tr_xyy_zzz[i] * tke_0 + 4.0 * tr_xyy_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 40-50 components of targeted buffer : FF

    auto tr_0_0_xx_xyz_xxx = pbuffer.data(idx_op_geom_020_ff + 40);

    auto tr_0_0_xx_xyz_xxy = pbuffer.data(idx_op_geom_020_ff + 41);

    auto tr_0_0_xx_xyz_xxz = pbuffer.data(idx_op_geom_020_ff + 42);

    auto tr_0_0_xx_xyz_xyy = pbuffer.data(idx_op_geom_020_ff + 43);

    auto tr_0_0_xx_xyz_xyz = pbuffer.data(idx_op_geom_020_ff + 44);

    auto tr_0_0_xx_xyz_xzz = pbuffer.data(idx_op_geom_020_ff + 45);

    auto tr_0_0_xx_xyz_yyy = pbuffer.data(idx_op_geom_020_ff + 46);

    auto tr_0_0_xx_xyz_yyz = pbuffer.data(idx_op_geom_020_ff + 47);

    auto tr_0_0_xx_xyz_yzz = pbuffer.data(idx_op_geom_020_ff + 48);

    auto tr_0_0_xx_xyz_zzz = pbuffer.data(idx_op_geom_020_ff + 49);

    #pragma omp simd aligned(tr_0_0_xx_xyz_xxx, tr_0_0_xx_xyz_xxy, tr_0_0_xx_xyz_xxz, tr_0_0_xx_xyz_xyy, tr_0_0_xx_xyz_xyz, tr_0_0_xx_xyz_xzz, tr_0_0_xx_xyz_yyy, tr_0_0_xx_xyz_yyz, tr_0_0_xx_xyz_yzz, tr_0_0_xx_xyz_zzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyz_xxx[i] = 6.0 * tr_yz_xx[i] - 4.0 * tr_yz_xxxx[i] * tke_0 + 6.0 * tr_xyz_x[i] - 6.0 * tr_xyz_xxx[i] * tbe_0 - 14.0 * tr_xyz_xxx[i] * tke_0 + 4.0 * tr_xyz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxyz_xx[i] * tbe_0 + 8.0 * tr_xxyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_xxy[i] = 4.0 * tr_yz_xy[i] - 4.0 * tr_yz_xxxy[i] * tke_0 + 2.0 * tr_xyz_y[i] - 6.0 * tr_xyz_xxy[i] * tbe_0 - 10.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxyz_xy[i] * tbe_0 + 8.0 * tr_xxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_xxz[i] = 4.0 * tr_yz_xz[i] - 4.0 * tr_yz_xxxz[i] * tke_0 + 2.0 * tr_xyz_z[i] - 6.0 * tr_xyz_xxz[i] * tbe_0 - 10.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxyz_xz[i] * tbe_0 + 8.0 * tr_xxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_xyy[i] = 2.0 * tr_yz_yy[i] - 4.0 * tr_yz_xxyy[i] * tke_0 - 6.0 * tr_xyz_xyy[i] * tbe_0 - 6.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_yy[i] * tbe_0 + 8.0 * tr_xxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_xyz[i] = 2.0 * tr_yz_yz[i] - 4.0 * tr_yz_xxyz[i] * tke_0 - 6.0 * tr_xyz_xyz[i] * tbe_0 - 6.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 + 8.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_xzz[i] = 2.0 * tr_yz_zz[i] - 4.0 * tr_yz_xxzz[i] * tke_0 - 6.0 * tr_xyz_xzz[i] * tbe_0 - 6.0 * tr_xyz_xzz[i] * tke_0 + 4.0 * tr_xyz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_zz[i] * tbe_0 + 8.0 * tr_xxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_yyy[i] = -4.0 * tr_yz_xyyy[i] * tke_0 - 6.0 * tr_xyz_yyy[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tke_0 + 4.0 * tr_xyz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_yyz[i] = -4.0 * tr_yz_xyyz[i] * tke_0 - 6.0 * tr_xyz_yyz[i] * tbe_0 - 2.0 * tr_xyz_yyz[i] * tke_0 + 4.0 * tr_xyz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_yzz[i] = -4.0 * tr_yz_xyzz[i] * tke_0 - 6.0 * tr_xyz_yzz[i] * tbe_0 - 2.0 * tr_xyz_yzz[i] * tke_0 + 4.0 * tr_xyz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyz_zzz[i] = -4.0 * tr_yz_xzzz[i] * tke_0 - 6.0 * tr_xyz_zzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tke_0 + 4.0 * tr_xyz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 50-60 components of targeted buffer : FF

    auto tr_0_0_xx_xzz_xxx = pbuffer.data(idx_op_geom_020_ff + 50);

    auto tr_0_0_xx_xzz_xxy = pbuffer.data(idx_op_geom_020_ff + 51);

    auto tr_0_0_xx_xzz_xxz = pbuffer.data(idx_op_geom_020_ff + 52);

    auto tr_0_0_xx_xzz_xyy = pbuffer.data(idx_op_geom_020_ff + 53);

    auto tr_0_0_xx_xzz_xyz = pbuffer.data(idx_op_geom_020_ff + 54);

    auto tr_0_0_xx_xzz_xzz = pbuffer.data(idx_op_geom_020_ff + 55);

    auto tr_0_0_xx_xzz_yyy = pbuffer.data(idx_op_geom_020_ff + 56);

    auto tr_0_0_xx_xzz_yyz = pbuffer.data(idx_op_geom_020_ff + 57);

    auto tr_0_0_xx_xzz_yzz = pbuffer.data(idx_op_geom_020_ff + 58);

    auto tr_0_0_xx_xzz_zzz = pbuffer.data(idx_op_geom_020_ff + 59);

    #pragma omp simd aligned(tr_0_0_xx_xzz_xxx, tr_0_0_xx_xzz_xxy, tr_0_0_xx_xzz_xxz, tr_0_0_xx_xzz_xyy, tr_0_0_xx_xzz_xyz, tr_0_0_xx_xzz_xzz, tr_0_0_xx_xzz_yyy, tr_0_0_xx_xzz_yyz, tr_0_0_xx_xzz_yzz, tr_0_0_xx_xzz_zzz, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_zz_xx, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xzz_xxx[i] = 6.0 * tr_zz_xx[i] - 4.0 * tr_zz_xxxx[i] * tke_0 + 6.0 * tr_xzz_x[i] - 6.0 * tr_xzz_xxx[i] * tbe_0 - 14.0 * tr_xzz_xxx[i] * tke_0 + 4.0 * tr_xzz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxzz_xx[i] * tbe_0 + 8.0 * tr_xxzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_xxy[i] = 4.0 * tr_zz_xy[i] - 4.0 * tr_zz_xxxy[i] * tke_0 + 2.0 * tr_xzz_y[i] - 6.0 * tr_xzz_xxy[i] * tbe_0 - 10.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxzz_xy[i] * tbe_0 + 8.0 * tr_xxzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_xxz[i] = 4.0 * tr_zz_xz[i] - 4.0 * tr_zz_xxxz[i] * tke_0 + 2.0 * tr_xzz_z[i] - 6.0 * tr_xzz_xxz[i] * tbe_0 - 10.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxzz_xz[i] * tbe_0 + 8.0 * tr_xxzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_xyy[i] = 2.0 * tr_zz_yy[i] - 4.0 * tr_zz_xxyy[i] * tke_0 - 6.0 * tr_xzz_xyy[i] * tbe_0 - 6.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_yy[i] * tbe_0 + 8.0 * tr_xxzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_xyz[i] = 2.0 * tr_zz_yz[i] - 4.0 * tr_zz_xxyz[i] * tke_0 - 6.0 * tr_xzz_xyz[i] * tbe_0 - 6.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_yz[i] * tbe_0 + 8.0 * tr_xxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_xzz[i] = 2.0 * tr_zz_zz[i] - 4.0 * tr_zz_xxzz[i] * tke_0 - 6.0 * tr_xzz_xzz[i] * tbe_0 - 6.0 * tr_xzz_xzz[i] * tke_0 + 4.0 * tr_xzz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_zz[i] * tbe_0 + 8.0 * tr_xxzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_yyy[i] = -4.0 * tr_zz_xyyy[i] * tke_0 - 6.0 * tr_xzz_yyy[i] * tbe_0 - 2.0 * tr_xzz_yyy[i] * tke_0 + 4.0 * tr_xzz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_yyz[i] = -4.0 * tr_zz_xyyz[i] * tke_0 - 6.0 * tr_xzz_yyz[i] * tbe_0 - 2.0 * tr_xzz_yyz[i] * tke_0 + 4.0 * tr_xzz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_yzz[i] = -4.0 * tr_zz_xyzz[i] * tke_0 - 6.0 * tr_xzz_yzz[i] * tbe_0 - 2.0 * tr_xzz_yzz[i] * tke_0 + 4.0 * tr_xzz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzz_zzz[i] = -4.0 * tr_zz_xzzz[i] * tke_0 - 6.0 * tr_xzz_zzz[i] * tbe_0 - 2.0 * tr_xzz_zzz[i] * tke_0 + 4.0 * tr_xzz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-70 components of targeted buffer : FF

    auto tr_0_0_xx_yyy_xxx = pbuffer.data(idx_op_geom_020_ff + 60);

    auto tr_0_0_xx_yyy_xxy = pbuffer.data(idx_op_geom_020_ff + 61);

    auto tr_0_0_xx_yyy_xxz = pbuffer.data(idx_op_geom_020_ff + 62);

    auto tr_0_0_xx_yyy_xyy = pbuffer.data(idx_op_geom_020_ff + 63);

    auto tr_0_0_xx_yyy_xyz = pbuffer.data(idx_op_geom_020_ff + 64);

    auto tr_0_0_xx_yyy_xzz = pbuffer.data(idx_op_geom_020_ff + 65);

    auto tr_0_0_xx_yyy_yyy = pbuffer.data(idx_op_geom_020_ff + 66);

    auto tr_0_0_xx_yyy_yyz = pbuffer.data(idx_op_geom_020_ff + 67);

    auto tr_0_0_xx_yyy_yzz = pbuffer.data(idx_op_geom_020_ff + 68);

    auto tr_0_0_xx_yyy_zzz = pbuffer.data(idx_op_geom_020_ff + 69);

    #pragma omp simd aligned(tr_0_0_xx_yyy_xxx, tr_0_0_xx_yyy_xxy, tr_0_0_xx_yyy_xxz, tr_0_0_xx_yyy_xyy, tr_0_0_xx_yyy_xyz, tr_0_0_xx_yyy_xzz, tr_0_0_xx_yyy_yyy, tr_0_0_xx_yyy_yyz, tr_0_0_xx_yyy_yzz, tr_0_0_xx_yyy_zzz, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxxxx, tr_yyy_xxxxy, tr_yyy_xxxxz, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyy_xxx[i] = 6.0 * tr_yyy_x[i] - 2.0 * tr_yyy_xxx[i] * tbe_0 - 14.0 * tr_yyy_xxx[i] * tke_0 + 4.0 * tr_yyy_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xyyy_xx[i] * tbe_0 + 8.0 * tr_xyyy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_xxy[i] = 2.0 * tr_yyy_y[i] - 2.0 * tr_yyy_xxy[i] * tbe_0 - 10.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyyy_xy[i] * tbe_0 + 8.0 * tr_xyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_xxz[i] = 2.0 * tr_yyy_z[i] - 2.0 * tr_yyy_xxz[i] * tbe_0 - 10.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyyy_xz[i] * tbe_0 + 8.0 * tr_xyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_xyy[i] = -2.0 * tr_yyy_xyy[i] * tbe_0 - 6.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_yy[i] * tbe_0 + 8.0 * tr_xyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_xyz[i] = -2.0 * tr_yyy_xyz[i] * tbe_0 - 6.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_yz[i] * tbe_0 + 8.0 * tr_xyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_xzz[i] = -2.0 * tr_yyy_xzz[i] * tbe_0 - 6.0 * tr_yyy_xzz[i] * tke_0 + 4.0 * tr_yyy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_zz[i] * tbe_0 + 8.0 * tr_xyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_yyy[i] = -2.0 * tr_yyy_yyy[i] * tbe_0 - 2.0 * tr_yyy_yyy[i] * tke_0 + 4.0 * tr_yyy_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_yyz[i] = -2.0 * tr_yyy_yyz[i] * tbe_0 - 2.0 * tr_yyy_yyz[i] * tke_0 + 4.0 * tr_yyy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_yzz[i] = -2.0 * tr_yyy_yzz[i] * tbe_0 - 2.0 * tr_yyy_yzz[i] * tke_0 + 4.0 * tr_yyy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyy_zzz[i] = -2.0 * tr_yyy_zzz[i] * tbe_0 - 2.0 * tr_yyy_zzz[i] * tke_0 + 4.0 * tr_yyy_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 70-80 components of targeted buffer : FF

    auto tr_0_0_xx_yyz_xxx = pbuffer.data(idx_op_geom_020_ff + 70);

    auto tr_0_0_xx_yyz_xxy = pbuffer.data(idx_op_geom_020_ff + 71);

    auto tr_0_0_xx_yyz_xxz = pbuffer.data(idx_op_geom_020_ff + 72);

    auto tr_0_0_xx_yyz_xyy = pbuffer.data(idx_op_geom_020_ff + 73);

    auto tr_0_0_xx_yyz_xyz = pbuffer.data(idx_op_geom_020_ff + 74);

    auto tr_0_0_xx_yyz_xzz = pbuffer.data(idx_op_geom_020_ff + 75);

    auto tr_0_0_xx_yyz_yyy = pbuffer.data(idx_op_geom_020_ff + 76);

    auto tr_0_0_xx_yyz_yyz = pbuffer.data(idx_op_geom_020_ff + 77);

    auto tr_0_0_xx_yyz_yzz = pbuffer.data(idx_op_geom_020_ff + 78);

    auto tr_0_0_xx_yyz_zzz = pbuffer.data(idx_op_geom_020_ff + 79);

    #pragma omp simd aligned(tr_0_0_xx_yyz_xxx, tr_0_0_xx_yyz_xxy, tr_0_0_xx_yyz_xxz, tr_0_0_xx_yyz_xyy, tr_0_0_xx_yyz_xyz, tr_0_0_xx_yyz_xzz, tr_0_0_xx_yyz_yyy, tr_0_0_xx_yyz_yyz, tr_0_0_xx_yyz_yzz, tr_0_0_xx_yyz_zzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyz_xxx[i] = 6.0 * tr_yyz_x[i] - 2.0 * tr_yyz_xxx[i] * tbe_0 - 14.0 * tr_yyz_xxx[i] * tke_0 + 4.0 * tr_yyz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xyyz_xx[i] * tbe_0 + 8.0 * tr_xyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_xxy[i] = 2.0 * tr_yyz_y[i] - 2.0 * tr_yyz_xxy[i] * tbe_0 - 10.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyyz_xy[i] * tbe_0 + 8.0 * tr_xyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_xxz[i] = 2.0 * tr_yyz_z[i] - 2.0 * tr_yyz_xxz[i] * tbe_0 - 10.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyyz_xz[i] * tbe_0 + 8.0 * tr_xyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_xyy[i] = -2.0 * tr_yyz_xyy[i] * tbe_0 - 6.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_yy[i] * tbe_0 + 8.0 * tr_xyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_xyz[i] = -2.0 * tr_yyz_xyz[i] * tbe_0 - 6.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 + 8.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_xzz[i] = -2.0 * tr_yyz_xzz[i] * tbe_0 - 6.0 * tr_yyz_xzz[i] * tke_0 + 4.0 * tr_yyz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_zz[i] * tbe_0 + 8.0 * tr_xyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_yyy[i] = -2.0 * tr_yyz_yyy[i] * tbe_0 - 2.0 * tr_yyz_yyy[i] * tke_0 + 4.0 * tr_yyz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_yyz[i] = -2.0 * tr_yyz_yyz[i] * tbe_0 - 2.0 * tr_yyz_yyz[i] * tke_0 + 4.0 * tr_yyz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_yzz[i] = -2.0 * tr_yyz_yzz[i] * tbe_0 - 2.0 * tr_yyz_yzz[i] * tke_0 + 4.0 * tr_yyz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyz_zzz[i] = -2.0 * tr_yyz_zzz[i] * tbe_0 - 2.0 * tr_yyz_zzz[i] * tke_0 + 4.0 * tr_yyz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 80-90 components of targeted buffer : FF

    auto tr_0_0_xx_yzz_xxx = pbuffer.data(idx_op_geom_020_ff + 80);

    auto tr_0_0_xx_yzz_xxy = pbuffer.data(idx_op_geom_020_ff + 81);

    auto tr_0_0_xx_yzz_xxz = pbuffer.data(idx_op_geom_020_ff + 82);

    auto tr_0_0_xx_yzz_xyy = pbuffer.data(idx_op_geom_020_ff + 83);

    auto tr_0_0_xx_yzz_xyz = pbuffer.data(idx_op_geom_020_ff + 84);

    auto tr_0_0_xx_yzz_xzz = pbuffer.data(idx_op_geom_020_ff + 85);

    auto tr_0_0_xx_yzz_yyy = pbuffer.data(idx_op_geom_020_ff + 86);

    auto tr_0_0_xx_yzz_yyz = pbuffer.data(idx_op_geom_020_ff + 87);

    auto tr_0_0_xx_yzz_yzz = pbuffer.data(idx_op_geom_020_ff + 88);

    auto tr_0_0_xx_yzz_zzz = pbuffer.data(idx_op_geom_020_ff + 89);

    #pragma omp simd aligned(tr_0_0_xx_yzz_xxx, tr_0_0_xx_yzz_xxy, tr_0_0_xx_yzz_xxz, tr_0_0_xx_yzz_xyy, tr_0_0_xx_yzz_xyz, tr_0_0_xx_yzz_xzz, tr_0_0_xx_yzz_yyy, tr_0_0_xx_yzz_yyz, tr_0_0_xx_yzz_yzz, tr_0_0_xx_yzz_zzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yzz_xxx[i] = 6.0 * tr_yzz_x[i] - 2.0 * tr_yzz_xxx[i] * tbe_0 - 14.0 * tr_yzz_xxx[i] * tke_0 + 4.0 * tr_yzz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xyzz_xx[i] * tbe_0 + 8.0 * tr_xyzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_xxy[i] = 2.0 * tr_yzz_y[i] - 2.0 * tr_yzz_xxy[i] * tbe_0 - 10.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyzz_xy[i] * tbe_0 + 8.0 * tr_xyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_xxz[i] = 2.0 * tr_yzz_z[i] - 2.0 * tr_yzz_xxz[i] * tbe_0 - 10.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyzz_xz[i] * tbe_0 + 8.0 * tr_xyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_xyy[i] = -2.0 * tr_yzz_xyy[i] * tbe_0 - 6.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_yy[i] * tbe_0 + 8.0 * tr_xyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_xyz[i] = -2.0 * tr_yzz_xyz[i] * tbe_0 - 6.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 + 8.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_xzz[i] = -2.0 * tr_yzz_xzz[i] * tbe_0 - 6.0 * tr_yzz_xzz[i] * tke_0 + 4.0 * tr_yzz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_zz[i] * tbe_0 + 8.0 * tr_xyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_yyy[i] = -2.0 * tr_yzz_yyy[i] * tbe_0 - 2.0 * tr_yzz_yyy[i] * tke_0 + 4.0 * tr_yzz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_yyz[i] = -2.0 * tr_yzz_yyz[i] * tbe_0 - 2.0 * tr_yzz_yyz[i] * tke_0 + 4.0 * tr_yzz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_yzz[i] = -2.0 * tr_yzz_yzz[i] * tbe_0 - 2.0 * tr_yzz_yzz[i] * tke_0 + 4.0 * tr_yzz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzz_zzz[i] = -2.0 * tr_yzz_zzz[i] * tbe_0 - 2.0 * tr_yzz_zzz[i] * tke_0 + 4.0 * tr_yzz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-100 components of targeted buffer : FF

    auto tr_0_0_xx_zzz_xxx = pbuffer.data(idx_op_geom_020_ff + 90);

    auto tr_0_0_xx_zzz_xxy = pbuffer.data(idx_op_geom_020_ff + 91);

    auto tr_0_0_xx_zzz_xxz = pbuffer.data(idx_op_geom_020_ff + 92);

    auto tr_0_0_xx_zzz_xyy = pbuffer.data(idx_op_geom_020_ff + 93);

    auto tr_0_0_xx_zzz_xyz = pbuffer.data(idx_op_geom_020_ff + 94);

    auto tr_0_0_xx_zzz_xzz = pbuffer.data(idx_op_geom_020_ff + 95);

    auto tr_0_0_xx_zzz_yyy = pbuffer.data(idx_op_geom_020_ff + 96);

    auto tr_0_0_xx_zzz_yyz = pbuffer.data(idx_op_geom_020_ff + 97);

    auto tr_0_0_xx_zzz_yzz = pbuffer.data(idx_op_geom_020_ff + 98);

    auto tr_0_0_xx_zzz_zzz = pbuffer.data(idx_op_geom_020_ff + 99);

    #pragma omp simd aligned(tr_0_0_xx_zzz_xxx, tr_0_0_xx_zzz_xxy, tr_0_0_xx_zzz_xxz, tr_0_0_xx_zzz_xyy, tr_0_0_xx_zzz_xyz, tr_0_0_xx_zzz_xzz, tr_0_0_xx_zzz_yyy, tr_0_0_xx_zzz_yyz, tr_0_0_xx_zzz_yzz, tr_0_0_xx_zzz_zzz, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxxxx, tr_zzz_xxxxy, tr_zzz_xxxxz, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_zzz_xxx[i] = 6.0 * tr_zzz_x[i] - 2.0 * tr_zzz_xxx[i] * tbe_0 - 14.0 * tr_zzz_xxx[i] * tke_0 + 4.0 * tr_zzz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xzzz_xx[i] * tbe_0 + 8.0 * tr_xzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_xxy[i] = 2.0 * tr_zzz_y[i] - 2.0 * tr_zzz_xxy[i] * tbe_0 - 10.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xzzz_xy[i] * tbe_0 + 8.0 * tr_xzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_xxz[i] = 2.0 * tr_zzz_z[i] - 2.0 * tr_zzz_xxz[i] * tbe_0 - 10.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xzzz_xz[i] * tbe_0 + 8.0 * tr_xzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_xyy[i] = -2.0 * tr_zzz_xyy[i] * tbe_0 - 6.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_yy[i] * tbe_0 + 8.0 * tr_xzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_xyz[i] = -2.0 * tr_zzz_xyz[i] * tbe_0 - 6.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_yz[i] * tbe_0 + 8.0 * tr_xzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_xzz[i] = -2.0 * tr_zzz_xzz[i] * tbe_0 - 6.0 * tr_zzz_xzz[i] * tke_0 + 4.0 * tr_zzz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_zz[i] * tbe_0 + 8.0 * tr_xzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_yyy[i] = -2.0 * tr_zzz_yyy[i] * tbe_0 - 2.0 * tr_zzz_yyy[i] * tke_0 + 4.0 * tr_zzz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_yyz[i] = -2.0 * tr_zzz_yyz[i] * tbe_0 - 2.0 * tr_zzz_yyz[i] * tke_0 + 4.0 * tr_zzz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_yzz[i] = -2.0 * tr_zzz_yzz[i] * tbe_0 - 2.0 * tr_zzz_yzz[i] * tke_0 + 4.0 * tr_zzz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzz_zzz[i] = -2.0 * tr_zzz_zzz[i] * tbe_0 - 2.0 * tr_zzz_zzz[i] * tke_0 + 4.0 * tr_zzz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 100-110 components of targeted buffer : FF

    auto tr_0_0_xy_xxx_xxx = pbuffer.data(idx_op_geom_020_ff + 100);

    auto tr_0_0_xy_xxx_xxy = pbuffer.data(idx_op_geom_020_ff + 101);

    auto tr_0_0_xy_xxx_xxz = pbuffer.data(idx_op_geom_020_ff + 102);

    auto tr_0_0_xy_xxx_xyy = pbuffer.data(idx_op_geom_020_ff + 103);

    auto tr_0_0_xy_xxx_xyz = pbuffer.data(idx_op_geom_020_ff + 104);

    auto tr_0_0_xy_xxx_xzz = pbuffer.data(idx_op_geom_020_ff + 105);

    auto tr_0_0_xy_xxx_yyy = pbuffer.data(idx_op_geom_020_ff + 106);

    auto tr_0_0_xy_xxx_yyz = pbuffer.data(idx_op_geom_020_ff + 107);

    auto tr_0_0_xy_xxx_yzz = pbuffer.data(idx_op_geom_020_ff + 108);

    auto tr_0_0_xy_xxx_zzz = pbuffer.data(idx_op_geom_020_ff + 109);

    #pragma omp simd aligned(tr_0_0_xy_xxx_xxx, tr_0_0_xy_xxx_xxy, tr_0_0_xy_xxx_xxz, tr_0_0_xy_xxx_xyy, tr_0_0_xy_xxx_xyz, tr_0_0_xy_xxx_xzz, tr_0_0_xy_xxx_yyy, tr_0_0_xy_xxx_yyz, tr_0_0_xy_xxx_yzz, tr_0_0_xy_xxx_zzz, tr_xx_xx, tr_xx_xxxy, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxxxy, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxx_xx, tr_xxxx_xxxy, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxxy_xxx, tr_xxxxy_xxy, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_zzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxx_xxx[i] = -6.0 * tr_xx_xxxy[i] * tke_0 - 6.0 * tr_xxy_xxx[i] * tbe_0 - 6.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_xxy[i] = 3.0 * tr_xx_xx[i] - 6.0 * tr_xx_xxyy[i] * tke_0 - 6.0 * tr_xxy_xxy[i] * tbe_0 + 2.0 * tr_xxx_x[i] - 4.0 * tr_xxx_xyy[i] * tke_0 - 2.0 * tr_xxx_xxx[i] * tke_0 + 4.0 * tr_xxx_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_xy[i] * tbe_0 + 4.0 * tr_xxxy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xx[i] * tbe_0 + 4.0 * tr_xxxx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_xxz[i] = -6.0 * tr_xx_xxyz[i] * tke_0 - 6.0 * tr_xxy_xxz[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_xz[i] * tbe_0 + 4.0 * tr_xxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_xyy[i] = 6.0 * tr_xx_xy[i] - 6.0 * tr_xx_xyyy[i] * tke_0 - 6.0 * tr_xxy_xyy[i] * tbe_0 + 2.0 * tr_xxx_y[i] - 2.0 * tr_xxx_yyy[i] * tke_0 - 4.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxy_yy[i] * tbe_0 + 4.0 * tr_xxxy_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxx_xy[i] * tbe_0 + 4.0 * tr_xxxx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_xyz[i] = 3.0 * tr_xx_xz[i] - 6.0 * tr_xx_xyyz[i] * tke_0 - 6.0 * tr_xxy_xyz[i] * tbe_0 + tr_xxx_z[i] - 2.0 * tr_xxx_yyz[i] * tke_0 - 2.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxy_yz[i] * tbe_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xz[i] * tbe_0 + 4.0 * tr_xxxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_xzz[i] = -6.0 * tr_xx_xyzz[i] * tke_0 - 6.0 * tr_xxy_xzz[i] * tbe_0 - 2.0 * tr_xxx_yzz[i] * tke_0 + 4.0 * tr_xxx_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxy_zz[i] * tbe_0 + 4.0 * tr_xxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_yyy[i] = 9.0 * tr_xx_yy[i] - 6.0 * tr_xx_yyyy[i] * tke_0 - 6.0 * tr_xxy_yyy[i] * tbe_0 - 6.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxx_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxxx_yy[i] * tbe_0 + 4.0 * tr_xxxx_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_yyz[i] = 6.0 * tr_xx_yz[i] - 6.0 * tr_xx_yyyz[i] * tke_0 - 6.0 * tr_xxy_yyz[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxx_yz[i] * tbe_0 + 4.0 * tr_xxxx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_yzz[i] = 3.0 * tr_xx_zz[i] - 6.0 * tr_xx_yyzz[i] * tke_0 - 6.0 * tr_xxy_yzz[i] * tbe_0 - 2.0 * tr_xxx_xzz[i] * tke_0 + 4.0 * tr_xxx_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zz[i] * tbe_0 + 4.0 * tr_xxxx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxx_zzz[i] = -6.0 * tr_xx_yzzz[i] * tke_0 - 6.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxx_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 110-120 components of targeted buffer : FF

    auto tr_0_0_xy_xxy_xxx = pbuffer.data(idx_op_geom_020_ff + 110);

    auto tr_0_0_xy_xxy_xxy = pbuffer.data(idx_op_geom_020_ff + 111);

    auto tr_0_0_xy_xxy_xxz = pbuffer.data(idx_op_geom_020_ff + 112);

    auto tr_0_0_xy_xxy_xyy = pbuffer.data(idx_op_geom_020_ff + 113);

    auto tr_0_0_xy_xxy_xyz = pbuffer.data(idx_op_geom_020_ff + 114);

    auto tr_0_0_xy_xxy_xzz = pbuffer.data(idx_op_geom_020_ff + 115);

    auto tr_0_0_xy_xxy_yyy = pbuffer.data(idx_op_geom_020_ff + 116);

    auto tr_0_0_xy_xxy_yyz = pbuffer.data(idx_op_geom_020_ff + 117);

    auto tr_0_0_xy_xxy_yzz = pbuffer.data(idx_op_geom_020_ff + 118);

    auto tr_0_0_xy_xxy_zzz = pbuffer.data(idx_op_geom_020_ff + 119);

    #pragma omp simd aligned(tr_0_0_xy_xxy_xxx, tr_0_0_xy_xxy_xxy, tr_0_0_xy_xxy_xxz, tr_0_0_xy_xxy_xyy, tr_0_0_xy_xxy_xyz, tr_0_0_xy_xxy_xzz, tr_0_0_xy_xxy_yyy, tr_0_0_xy_xxy_yyz, tr_0_0_xy_xxy_yzz, tr_0_0_xy_xxy_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xx, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxy_xx, tr_xxxy_xxxy, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxy_xxx[i] = 2.0 * tr_x_xxx[i] - 4.0 * tr_xy_xxxy[i] * tke_0 - 4.0 * tr_xyy_xxx[i] * tbe_0 + 3.0 * tr_xx_xx[i] - 2.0 * tr_xx_xxxx[i] * tke_0 - 6.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_xxy[i] = 2.0 * tr_x_xxy[i] + 2.0 * tr_xy_xx[i] - 4.0 * tr_xy_xxyy[i] * tke_0 - 4.0 * tr_xyy_xxy[i] * tbe_0 + 2.0 * tr_xx_xy[i] - 2.0 * tr_xx_xxxy[i] * tke_0 + 2.0 * tr_xxy_x[i] - 4.0 * tr_xxy_xyy[i] * tke_0 - 2.0 * tr_xxy_xxx[i] * tke_0 + 4.0 * tr_xxy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_xy[i] * tbe_0 + 4.0 * tr_xxyy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxy[i] * tbe_0 - 2.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_xxz[i] = 2.0 * tr_x_xxz[i] - 4.0 * tr_xy_xxyz[i] * tke_0 - 4.0 * tr_xyy_xxz[i] * tbe_0 + 2.0 * tr_xx_xz[i] - 2.0 * tr_xx_xxxz[i] * tke_0 - 4.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_xz[i] * tbe_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_xyy[i] = 2.0 * tr_x_xyy[i] + 4.0 * tr_xy_xy[i] - 4.0 * tr_xy_xyyy[i] * tke_0 - 4.0 * tr_xyy_xyy[i] * tbe_0 + tr_xx_yy[i] - 2.0 * tr_xx_xxyy[i] * tke_0 + 2.0 * tr_xxy_y[i] - 2.0 * tr_xxy_yyy[i] * tke_0 - 4.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyy_yy[i] * tbe_0 + 4.0 * tr_xxyy_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyy[i] * tbe_0 - 4.0 * tr_xxxy_xy[i] * tbe_0 + 4.0 * tr_xxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_xyz[i] = 2.0 * tr_x_xyz[i] + 2.0 * tr_xy_xz[i] - 4.0 * tr_xy_xyyz[i] * tke_0 - 4.0 * tr_xyy_xyz[i] * tbe_0 + tr_xx_yz[i] - 2.0 * tr_xx_xxyz[i] * tke_0 + tr_xxy_z[i] - 2.0 * tr_xxy_yyz[i] * tke_0 - 2.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyy_yz[i] * tbe_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyz[i] * tbe_0 - 2.0 * tr_xxxy_xz[i] * tbe_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_xzz[i] = 2.0 * tr_x_xzz[i] - 4.0 * tr_xy_xyzz[i] * tke_0 - 4.0 * tr_xyy_xzz[i] * tbe_0 + tr_xx_zz[i] - 2.0 * tr_xx_xxzz[i] * tke_0 - 2.0 * tr_xxy_yzz[i] * tke_0 + 4.0 * tr_xxy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyy_zz[i] * tbe_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_yyy[i] = 2.0 * tr_x_yyy[i] + 6.0 * tr_xy_yy[i] - 4.0 * tr_xy_yyyy[i] * tke_0 - 4.0 * tr_xyy_yyy[i] * tbe_0 - 2.0 * tr_xx_xyyy[i] * tke_0 - 6.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxy_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyy[i] * tbe_0 - 6.0 * tr_xxxy_yy[i] * tbe_0 + 4.0 * tr_xxxy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_yyz[i] = 2.0 * tr_x_yyz[i] + 4.0 * tr_xy_yz[i] - 4.0 * tr_xy_yyyz[i] * tke_0 - 4.0 * tr_xyy_yyz[i] * tbe_0 - 2.0 * tr_xx_xyyz[i] * tke_0 - 4.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyz[i] * tbe_0 - 4.0 * tr_xxxy_yz[i] * tbe_0 + 4.0 * tr_xxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_yzz[i] = 2.0 * tr_x_yzz[i] + 2.0 * tr_xy_zz[i] - 4.0 * tr_xy_yyzz[i] * tke_0 - 4.0 * tr_xyy_yzz[i] * tbe_0 - 2.0 * tr_xx_xyzz[i] * tke_0 - 2.0 * tr_xxy_xzz[i] * tke_0 + 4.0 * tr_xxy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yzz[i] * tbe_0 - 2.0 * tr_xxxy_zz[i] * tbe_0 + 4.0 * tr_xxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxy_zzz[i] = 2.0 * tr_x_zzz[i] - 4.0 * tr_xy_yzzz[i] * tke_0 - 4.0 * tr_xyy_zzz[i] * tbe_0 - 2.0 * tr_xx_xzzz[i] * tke_0 + 4.0 * tr_xxy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zzz[i] * tbe_0 + 4.0 * tr_xxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-130 components of targeted buffer : FF

    auto tr_0_0_xy_xxz_xxx = pbuffer.data(idx_op_geom_020_ff + 120);

    auto tr_0_0_xy_xxz_xxy = pbuffer.data(idx_op_geom_020_ff + 121);

    auto tr_0_0_xy_xxz_xxz = pbuffer.data(idx_op_geom_020_ff + 122);

    auto tr_0_0_xy_xxz_xyy = pbuffer.data(idx_op_geom_020_ff + 123);

    auto tr_0_0_xy_xxz_xyz = pbuffer.data(idx_op_geom_020_ff + 124);

    auto tr_0_0_xy_xxz_xzz = pbuffer.data(idx_op_geom_020_ff + 125);

    auto tr_0_0_xy_xxz_yyy = pbuffer.data(idx_op_geom_020_ff + 126);

    auto tr_0_0_xy_xxz_yyz = pbuffer.data(idx_op_geom_020_ff + 127);

    auto tr_0_0_xy_xxz_yzz = pbuffer.data(idx_op_geom_020_ff + 128);

    auto tr_0_0_xy_xxz_zzz = pbuffer.data(idx_op_geom_020_ff + 129);

    #pragma omp simd aligned(tr_0_0_xy_xxz_xxx, tr_0_0_xy_xxz_xxy, tr_0_0_xy_xxz_xxz, tr_0_0_xy_xxz_xyy, tr_0_0_xy_xxz_xyz, tr_0_0_xy_xxz_xzz, tr_0_0_xy_xxz_yyy, tr_0_0_xy_xxz_yyz, tr_0_0_xy_xxz_yzz, tr_0_0_xy_xxz_zzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxxz_xx, tr_xxxz_xxxy, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxz_xxx[i] = -4.0 * tr_xz_xxxy[i] * tke_0 - 4.0 * tr_xyz_xxx[i] * tbe_0 - 6.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_xxy[i] = 2.0 * tr_xz_xx[i] - 4.0 * tr_xz_xxyy[i] * tke_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 + 2.0 * tr_xxz_x[i] - 4.0 * tr_xxz_xyy[i] * tke_0 - 2.0 * tr_xxz_xxx[i] * tke_0 + 4.0 * tr_xxz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_xxz[i] = -4.0 * tr_xz_xxyz[i] * tke_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_xyy[i] = 4.0 * tr_xz_xy[i] - 4.0 * tr_xz_xyyy[i] * tke_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 + 2.0 * tr_xxz_y[i] - 2.0 * tr_xxz_yyy[i] * tke_0 - 4.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xy[i] * tbe_0 + 4.0 * tr_xxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_xyz[i] = 2.0 * tr_xz_xz[i] - 4.0 * tr_xz_xyyz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + tr_xxz_z[i] - 2.0 * tr_xxz_yyz[i] * tke_0 - 2.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xz[i] * tbe_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_xzz[i] = -4.0 * tr_xz_xyzz[i] * tke_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 - 2.0 * tr_xxz_yzz[i] * tke_0 + 4.0 * tr_xxz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_yyy[i] = 6.0 * tr_xz_yy[i] - 4.0 * tr_xz_yyyy[i] * tke_0 - 4.0 * tr_xyz_yyy[i] * tbe_0 - 6.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxxz_yy[i] * tbe_0 + 4.0 * tr_xxxz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_yyz[i] = 4.0 * tr_xz_yz[i] - 4.0 * tr_xz_yyyz[i] * tke_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yz[i] * tbe_0 + 4.0 * tr_xxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_yzz[i] = 2.0 * tr_xz_zz[i] - 4.0 * tr_xz_yyzz[i] * tke_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 - 2.0 * tr_xxz_xzz[i] * tke_0 + 4.0 * tr_xxz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_zz[i] * tbe_0 + 4.0 * tr_xxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxz_zzz[i] = -4.0 * tr_xz_yzzz[i] * tke_0 - 4.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xxz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 130-140 components of targeted buffer : FF

    auto tr_0_0_xy_xyy_xxx = pbuffer.data(idx_op_geom_020_ff + 130);

    auto tr_0_0_xy_xyy_xxy = pbuffer.data(idx_op_geom_020_ff + 131);

    auto tr_0_0_xy_xyy_xxz = pbuffer.data(idx_op_geom_020_ff + 132);

    auto tr_0_0_xy_xyy_xyy = pbuffer.data(idx_op_geom_020_ff + 133);

    auto tr_0_0_xy_xyy_xyz = pbuffer.data(idx_op_geom_020_ff + 134);

    auto tr_0_0_xy_xyy_xzz = pbuffer.data(idx_op_geom_020_ff + 135);

    auto tr_0_0_xy_xyy_yyy = pbuffer.data(idx_op_geom_020_ff + 136);

    auto tr_0_0_xy_xyy_yyz = pbuffer.data(idx_op_geom_020_ff + 137);

    auto tr_0_0_xy_xyy_yzz = pbuffer.data(idx_op_geom_020_ff + 138);

    auto tr_0_0_xy_xyy_zzz = pbuffer.data(idx_op_geom_020_ff + 139);

    #pragma omp simd aligned(tr_0_0_xy_xyy_xxx, tr_0_0_xy_xyy_xxy, tr_0_0_xy_xyy_xxz, tr_0_0_xy_xyy_xyy, tr_0_0_xy_xyy_xyz, tr_0_0_xy_xyy_xzz, tr_0_0_xy_xyy_yyy, tr_0_0_xy_xyy_yyz, tr_0_0_xy_xyy_yzz, tr_0_0_xy_xyy_zzz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyy_xx, tr_xxyy_xxxy, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yy_xx, tr_yy_xxxy, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyy_xxx[i] = 2.0 * tr_y_xxx[i] - 2.0 * tr_yy_xxxy[i] * tke_0 - 2.0 * tr_yyy_xxx[i] * tbe_0 + 6.0 * tr_xy_xx[i] - 4.0 * tr_xy_xxxx[i] * tke_0 - 6.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_xxy[i] = 2.0 * tr_y_xxy[i] + tr_yy_xx[i] - 2.0 * tr_yy_xxyy[i] * tke_0 - 2.0 * tr_yyy_xxy[i] * tbe_0 + 4.0 * tr_xy_xy[i] - 4.0 * tr_xy_xxxy[i] * tke_0 + 2.0 * tr_xyy_x[i] - 4.0 * tr_xyy_xyy[i] * tke_0 - 2.0 * tr_xyy_xxx[i] * tke_0 + 4.0 * tr_xyy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_xy[i] * tbe_0 + 4.0 * tr_xyyy_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxy[i] * tbe_0 - 2.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_xxz[i] = 2.0 * tr_y_xxz[i] - 2.0 * tr_yy_xxyz[i] * tke_0 - 2.0 * tr_yyy_xxz[i] * tbe_0 + 4.0 * tr_xy_xz[i] - 4.0 * tr_xy_xxxz[i] * tke_0 - 4.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_xz[i] * tbe_0 + 4.0 * tr_xyyy_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_xyy[i] = 2.0 * tr_y_xyy[i] + 2.0 * tr_yy_xy[i] - 2.0 * tr_yy_xyyy[i] * tke_0 - 2.0 * tr_yyy_xyy[i] * tbe_0 + 2.0 * tr_xy_yy[i] - 4.0 * tr_xy_xxyy[i] * tke_0 + 2.0 * tr_xyy_y[i] - 2.0 * tr_xyy_yyy[i] * tke_0 - 4.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyy_yy[i] * tbe_0 + 4.0 * tr_xyyy_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xyy[i] * tbe_0 - 4.0 * tr_xxyy_xy[i] * tbe_0 + 4.0 * tr_xxyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_xyz[i] = 2.0 * tr_y_xyz[i] + tr_yy_xz[i] - 2.0 * tr_yy_xyyz[i] * tke_0 - 2.0 * tr_yyy_xyz[i] * tbe_0 + 2.0 * tr_xy_yz[i] - 4.0 * tr_xy_xxyz[i] * tke_0 + tr_xyy_z[i] - 2.0 * tr_xyy_yyz[i] * tke_0 - 2.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyy_yz[i] * tbe_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xyz[i] * tbe_0 - 2.0 * tr_xxyy_xz[i] * tbe_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_xzz[i] = 2.0 * tr_y_xzz[i] - 2.0 * tr_yy_xyzz[i] * tke_0 - 2.0 * tr_yyy_xzz[i] * tbe_0 + 2.0 * tr_xy_zz[i] - 4.0 * tr_xy_xxzz[i] * tke_0 - 2.0 * tr_xyy_yzz[i] * tke_0 + 4.0 * tr_xyy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_yyy[i] = 2.0 * tr_y_yyy[i] + 3.0 * tr_yy_yy[i] - 2.0 * tr_yy_yyyy[i] * tke_0 - 2.0 * tr_yyy_yyy[i] * tbe_0 - 4.0 * tr_xy_xyyy[i] * tke_0 - 6.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyy_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yyy[i] * tbe_0 - 6.0 * tr_xxyy_yy[i] * tbe_0 + 4.0 * tr_xxyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_yyz[i] = 2.0 * tr_y_yyz[i] + 2.0 * tr_yy_yz[i] - 2.0 * tr_yy_yyyz[i] * tke_0 - 2.0 * tr_yyy_yyz[i] * tbe_0 - 4.0 * tr_xy_xyyz[i] * tke_0 - 4.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yyz[i] * tbe_0 - 4.0 * tr_xxyy_yz[i] * tbe_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_yzz[i] = 2.0 * tr_y_yzz[i] + tr_yy_zz[i] - 2.0 * tr_yy_yyzz[i] * tke_0 - 2.0 * tr_yyy_yzz[i] * tbe_0 - 4.0 * tr_xy_xyzz[i] * tke_0 - 2.0 * tr_xyy_xzz[i] * tke_0 + 4.0 * tr_xyy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yzz[i] * tbe_0 - 2.0 * tr_xxyy_zz[i] * tbe_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyy_zzz[i] = 2.0 * tr_y_zzz[i] - 2.0 * tr_yy_yzzz[i] * tke_0 - 2.0 * tr_yyy_zzz[i] * tbe_0 - 4.0 * tr_xy_xzzz[i] * tke_0 + 4.0 * tr_xyy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyy_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 140-150 components of targeted buffer : FF

    auto tr_0_0_xy_xyz_xxx = pbuffer.data(idx_op_geom_020_ff + 140);

    auto tr_0_0_xy_xyz_xxy = pbuffer.data(idx_op_geom_020_ff + 141);

    auto tr_0_0_xy_xyz_xxz = pbuffer.data(idx_op_geom_020_ff + 142);

    auto tr_0_0_xy_xyz_xyy = pbuffer.data(idx_op_geom_020_ff + 143);

    auto tr_0_0_xy_xyz_xyz = pbuffer.data(idx_op_geom_020_ff + 144);

    auto tr_0_0_xy_xyz_xzz = pbuffer.data(idx_op_geom_020_ff + 145);

    auto tr_0_0_xy_xyz_yyy = pbuffer.data(idx_op_geom_020_ff + 146);

    auto tr_0_0_xy_xyz_yyz = pbuffer.data(idx_op_geom_020_ff + 147);

    auto tr_0_0_xy_xyz_yzz = pbuffer.data(idx_op_geom_020_ff + 148);

    auto tr_0_0_xy_xyz_zzz = pbuffer.data(idx_op_geom_020_ff + 149);

    #pragma omp simd aligned(tr_0_0_xy_xyz_xxx, tr_0_0_xy_xyz_xxy, tr_0_0_xy_xyz_xxz, tr_0_0_xy_xyz_xyy, tr_0_0_xy_xyz_xyz, tr_0_0_xy_xyz_xzz, tr_0_0_xy_xyz_yyy, tr_0_0_xy_xyz_yyz, tr_0_0_xy_xyz_yzz, tr_0_0_xy_xyz_zzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xxyz_xx, tr_xxyz_xxxy, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyz_xxx[i] = tr_z_xxx[i] - 2.0 * tr_yz_xxxy[i] * tke_0 - 2.0 * tr_yyz_xxx[i] * tbe_0 + 3.0 * tr_xz_xx[i] - 2.0 * tr_xz_xxxx[i] * tke_0 - 6.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_xxy[i] = tr_z_xxy[i] + tr_yz_xx[i] - 2.0 * tr_yz_xxyy[i] * tke_0 - 2.0 * tr_yyz_xxy[i] * tbe_0 + 2.0 * tr_xz_xy[i] - 2.0 * tr_xz_xxxy[i] * tke_0 + 2.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xyy[i] * tke_0 - 2.0 * tr_xyz_xxx[i] * tke_0 + 4.0 * tr_xyz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxy[i] * tbe_0 - 2.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_xxz[i] = tr_z_xxz[i] - 2.0 * tr_yz_xxyz[i] * tke_0 - 2.0 * tr_yyz_xxz[i] * tbe_0 + 2.0 * tr_xz_xz[i] - 2.0 * tr_xz_xxxz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_xyy[i] = tr_z_xyy[i] + 2.0 * tr_yz_xy[i] - 2.0 * tr_yz_xyyy[i] * tke_0 - 2.0 * tr_yyz_xyy[i] * tbe_0 + tr_xz_yy[i] - 2.0 * tr_xz_xxyy[i] * tke_0 + 2.0 * tr_xyz_y[i] - 2.0 * tr_xyz_yyy[i] * tke_0 - 4.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xyy[i] * tbe_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_xyz[i] = tr_z_xyz[i] + tr_yz_xz[i] - 2.0 * tr_yz_xyyz[i] * tke_0 - 2.0 * tr_yyz_xyz[i] * tbe_0 + tr_xz_yz[i] - 2.0 * tr_xz_xxyz[i] * tke_0 + tr_xyz_z[i] - 2.0 * tr_xyz_yyz[i] * tke_0 - 2.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xyz[i] * tbe_0 - 2.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_xzz[i] = tr_z_xzz[i] - 2.0 * tr_yz_xyzz[i] * tke_0 - 2.0 * tr_yyz_xzz[i] * tbe_0 + tr_xz_zz[i] - 2.0 * tr_xz_xxzz[i] * tke_0 - 2.0 * tr_xyz_yzz[i] * tke_0 + 4.0 * tr_xyz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_yyy[i] = tr_z_yyy[i] + 3.0 * tr_yz_yy[i] - 2.0 * tr_yz_yyyy[i] * tke_0 - 2.0 * tr_yyz_yyy[i] * tbe_0 - 2.0 * tr_xz_xyyy[i] * tke_0 - 6.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yyy[i] * tbe_0 - 6.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_yyz[i] = tr_z_yyz[i] + 2.0 * tr_yz_yz[i] - 2.0 * tr_yz_yyyz[i] * tke_0 - 2.0 * tr_yyz_yyz[i] * tbe_0 - 2.0 * tr_xz_xyyz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yyz[i] * tbe_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_yzz[i] = tr_z_yzz[i] + tr_yz_zz[i] - 2.0 * tr_yz_yyzz[i] * tke_0 - 2.0 * tr_yyz_yzz[i] * tbe_0 - 2.0 * tr_xz_xyzz[i] * tke_0 - 2.0 * tr_xyz_xzz[i] * tke_0 + 4.0 * tr_xyz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yzz[i] * tbe_0 - 2.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyz_zzz[i] = tr_z_zzz[i] - 2.0 * tr_yz_yzzz[i] * tke_0 - 2.0 * tr_yyz_zzz[i] * tbe_0 - 2.0 * tr_xz_xzzz[i] * tke_0 + 4.0 * tr_xyz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-160 components of targeted buffer : FF

    auto tr_0_0_xy_xzz_xxx = pbuffer.data(idx_op_geom_020_ff + 150);

    auto tr_0_0_xy_xzz_xxy = pbuffer.data(idx_op_geom_020_ff + 151);

    auto tr_0_0_xy_xzz_xxz = pbuffer.data(idx_op_geom_020_ff + 152);

    auto tr_0_0_xy_xzz_xyy = pbuffer.data(idx_op_geom_020_ff + 153);

    auto tr_0_0_xy_xzz_xyz = pbuffer.data(idx_op_geom_020_ff + 154);

    auto tr_0_0_xy_xzz_xzz = pbuffer.data(idx_op_geom_020_ff + 155);

    auto tr_0_0_xy_xzz_yyy = pbuffer.data(idx_op_geom_020_ff + 156);

    auto tr_0_0_xy_xzz_yyz = pbuffer.data(idx_op_geom_020_ff + 157);

    auto tr_0_0_xy_xzz_yzz = pbuffer.data(idx_op_geom_020_ff + 158);

    auto tr_0_0_xy_xzz_zzz = pbuffer.data(idx_op_geom_020_ff + 159);

    #pragma omp simd aligned(tr_0_0_xy_xzz_xxx, tr_0_0_xy_xzz_xxy, tr_0_0_xy_xzz_xxz, tr_0_0_xy_xzz_xyy, tr_0_0_xy_xzz_xyz, tr_0_0_xy_xzz_xzz, tr_0_0_xy_xzz_yyy, tr_0_0_xy_xzz_yyz, tr_0_0_xy_xzz_yzz, tr_0_0_xy_xzz_zzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xxzz_xx, tr_xxzz_xxxy, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_zz_xx, tr_zz_xxxy, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xzz_xxx[i] = -2.0 * tr_zz_xxxy[i] * tke_0 - 2.0 * tr_yzz_xxx[i] * tbe_0 - 6.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_xxy[i] = tr_zz_xx[i] - 2.0 * tr_zz_xxyy[i] * tke_0 - 2.0 * tr_yzz_xxy[i] * tbe_0 + 2.0 * tr_xzz_x[i] - 4.0 * tr_xzz_xyy[i] * tke_0 - 2.0 * tr_xzz_xxx[i] * tke_0 + 4.0 * tr_xzz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_xxz[i] = -2.0 * tr_zz_xxyz[i] * tke_0 - 2.0 * tr_yzz_xxz[i] * tbe_0 - 4.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_xyy[i] = 2.0 * tr_zz_xy[i] - 2.0 * tr_zz_xyyy[i] * tke_0 - 2.0 * tr_yzz_xyy[i] * tbe_0 + 2.0 * tr_xzz_y[i] - 2.0 * tr_xzz_yyy[i] * tke_0 - 4.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxzz_xy[i] * tbe_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_xyz[i] = tr_zz_xz[i] - 2.0 * tr_zz_xyyz[i] * tke_0 - 2.0 * tr_yzz_xyz[i] * tbe_0 + tr_xzz_z[i] - 2.0 * tr_xzz_yyz[i] * tke_0 - 2.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xz[i] * tbe_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_xzz[i] = -2.0 * tr_zz_xyzz[i] * tke_0 - 2.0 * tr_yzz_xzz[i] * tbe_0 - 2.0 * tr_xzz_yzz[i] * tke_0 + 4.0 * tr_xzz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_yyy[i] = 3.0 * tr_zz_yy[i] - 2.0 * tr_zz_yyyy[i] * tke_0 - 2.0 * tr_yzz_yyy[i] * tbe_0 - 6.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yy[i] * tbe_0 + 4.0 * tr_xxzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_yyz[i] = 2.0 * tr_zz_yz[i] - 2.0 * tr_zz_yyyz[i] * tke_0 - 2.0 * tr_yzz_yyz[i] * tbe_0 - 4.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxzz_yz[i] * tbe_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_yzz[i] = tr_zz_zz[i] - 2.0 * tr_zz_yyzz[i] * tke_0 - 2.0 * tr_yzz_yzz[i] * tbe_0 - 2.0 * tr_xzz_xzz[i] * tke_0 + 4.0 * tr_xzz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzz_zzz[i] = -2.0 * tr_zz_yzzz[i] * tke_0 - 2.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 160-170 components of targeted buffer : FF

    auto tr_0_0_xy_yyy_xxx = pbuffer.data(idx_op_geom_020_ff + 160);

    auto tr_0_0_xy_yyy_xxy = pbuffer.data(idx_op_geom_020_ff + 161);

    auto tr_0_0_xy_yyy_xxz = pbuffer.data(idx_op_geom_020_ff + 162);

    auto tr_0_0_xy_yyy_xyy = pbuffer.data(idx_op_geom_020_ff + 163);

    auto tr_0_0_xy_yyy_xyz = pbuffer.data(idx_op_geom_020_ff + 164);

    auto tr_0_0_xy_yyy_xzz = pbuffer.data(idx_op_geom_020_ff + 165);

    auto tr_0_0_xy_yyy_yyy = pbuffer.data(idx_op_geom_020_ff + 166);

    auto tr_0_0_xy_yyy_yyz = pbuffer.data(idx_op_geom_020_ff + 167);

    auto tr_0_0_xy_yyy_yzz = pbuffer.data(idx_op_geom_020_ff + 168);

    auto tr_0_0_xy_yyy_zzz = pbuffer.data(idx_op_geom_020_ff + 169);

    #pragma omp simd aligned(tr_0_0_xy_yyy_xxx, tr_0_0_xy_yyy_xxy, tr_0_0_xy_yyy_xxz, tr_0_0_xy_yyy_xyy, tr_0_0_xy_yyy_xyz, tr_0_0_xy_yyy_xzz, tr_0_0_xy_yyy_yyy, tr_0_0_xy_yyy_yyz, tr_0_0_xy_yyy_yzz, tr_0_0_xy_yyy_zzz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyy_xx, tr_xyyy_xxxy, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyyy_xxx, tr_xyyyy_xxy, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_zzz, tr_yy_xx, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxxxy, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyy_xxx[i] = 9.0 * tr_yy_xx[i] - 6.0 * tr_yy_xxxx[i] * tke_0 - 6.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yyyy_xx[i] * tbe_0 + 4.0 * tr_yyyy_xxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_xxy[i] = 6.0 * tr_yy_xy[i] - 6.0 * tr_yy_xxxy[i] * tke_0 + 2.0 * tr_yyy_x[i] - 4.0 * tr_yyy_xyy[i] * tke_0 - 2.0 * tr_yyy_xxx[i] * tke_0 + 4.0 * tr_yyy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyy_xy[i] * tbe_0 + 4.0 * tr_yyyy_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxy[i] * tbe_0 - 2.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_xxz[i] = 6.0 * tr_yy_xz[i] - 6.0 * tr_yy_xxxz[i] * tke_0 - 4.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyy_xz[i] * tbe_0 + 4.0 * tr_yyyy_xxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_xyy[i] = 3.0 * tr_yy_yy[i] - 6.0 * tr_yy_xxyy[i] * tke_0 + 2.0 * tr_yyy_y[i] - 2.0 * tr_yyy_yyy[i] * tke_0 - 4.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyy_yy[i] * tbe_0 + 4.0 * tr_yyyy_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xyy[i] * tbe_0 - 4.0 * tr_xyyy_xy[i] * tbe_0 + 4.0 * tr_xyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_xyz[i] = 3.0 * tr_yy_yz[i] - 6.0 * tr_yy_xxyz[i] * tke_0 + tr_yyy_z[i] - 2.0 * tr_yyy_yyz[i] * tke_0 - 2.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyy_yz[i] * tbe_0 + 4.0 * tr_yyyy_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xyz[i] * tbe_0 - 2.0 * tr_xyyy_xz[i] * tbe_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_xzz[i] = 3.0 * tr_yy_zz[i] - 6.0 * tr_yy_xxzz[i] * tke_0 - 2.0 * tr_yyy_yzz[i] * tke_0 + 4.0 * tr_yyy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyy_zz[i] * tbe_0 + 4.0 * tr_yyyy_xxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_yyy[i] = -6.0 * tr_yy_xyyy[i] * tke_0 - 6.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyy_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yyy[i] * tbe_0 - 6.0 * tr_xyyy_yy[i] * tbe_0 + 4.0 * tr_xyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_yyz[i] = -6.0 * tr_yy_xyyz[i] * tke_0 - 4.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yyz[i] * tbe_0 - 4.0 * tr_xyyy_yz[i] * tbe_0 + 4.0 * tr_xyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_yzz[i] = -6.0 * tr_yy_xyzz[i] * tke_0 - 2.0 * tr_yyy_xzz[i] * tke_0 + 4.0 * tr_yyy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yzz[i] * tbe_0 - 2.0 * tr_xyyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyy_zzz[i] = -6.0 * tr_yy_xzzz[i] * tke_0 + 4.0 * tr_yyy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyy_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 170-180 components of targeted buffer : FF

    auto tr_0_0_xy_yyz_xxx = pbuffer.data(idx_op_geom_020_ff + 170);

    auto tr_0_0_xy_yyz_xxy = pbuffer.data(idx_op_geom_020_ff + 171);

    auto tr_0_0_xy_yyz_xxz = pbuffer.data(idx_op_geom_020_ff + 172);

    auto tr_0_0_xy_yyz_xyy = pbuffer.data(idx_op_geom_020_ff + 173);

    auto tr_0_0_xy_yyz_xyz = pbuffer.data(idx_op_geom_020_ff + 174);

    auto tr_0_0_xy_yyz_xzz = pbuffer.data(idx_op_geom_020_ff + 175);

    auto tr_0_0_xy_yyz_yyy = pbuffer.data(idx_op_geom_020_ff + 176);

    auto tr_0_0_xy_yyz_yyz = pbuffer.data(idx_op_geom_020_ff + 177);

    auto tr_0_0_xy_yyz_yzz = pbuffer.data(idx_op_geom_020_ff + 178);

    auto tr_0_0_xy_yyz_zzz = pbuffer.data(idx_op_geom_020_ff + 179);

    #pragma omp simd aligned(tr_0_0_xy_yyz_xxx, tr_0_0_xy_yyz_xxy, tr_0_0_xy_yyz_xxz, tr_0_0_xy_yyz_xyy, tr_0_0_xy_yyz_xyz, tr_0_0_xy_yyz_xzz, tr_0_0_xy_yyz_yyy, tr_0_0_xy_yyz_yyz, tr_0_0_xy_yyz_yzz, tr_0_0_xy_yyz_zzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxy, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyz_xxx[i] = 6.0 * tr_yz_xx[i] - 4.0 * tr_yz_xxxx[i] * tke_0 - 6.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_xxy[i] = 4.0 * tr_yz_xy[i] - 4.0 * tr_yz_xxxy[i] * tke_0 + 2.0 * tr_yyz_x[i] - 4.0 * tr_yyz_xyy[i] * tke_0 - 2.0 * tr_yyz_xxx[i] * tke_0 + 4.0 * tr_yyz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xy[i] * tbe_0 + 4.0 * tr_yyyz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 - 2.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_xxz[i] = 4.0 * tr_yz_xz[i] - 4.0 * tr_yz_xxxz[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xz[i] * tbe_0 + 4.0 * tr_yyyz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_xyy[i] = 2.0 * tr_yz_yy[i] - 4.0 * tr_yz_xxyy[i] * tke_0 + 2.0 * tr_yyz_y[i] - 2.0 * tr_yyz_yyy[i] * tke_0 - 4.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_yy[i] * tbe_0 + 4.0 * tr_yyyz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_xyz[i] = 2.0 * tr_yz_yz[i] - 4.0 * tr_yz_xxyz[i] * tke_0 + tr_yyz_z[i] - 2.0 * tr_yyz_yyz[i] * tke_0 - 2.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_yz[i] * tbe_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 - 2.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_xzz[i] = 2.0 * tr_yz_zz[i] - 4.0 * tr_yz_xxzz[i] * tke_0 - 2.0 * tr_yyz_yzz[i] * tke_0 + 4.0 * tr_yyz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_zz[i] * tbe_0 + 4.0 * tr_yyyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_yyy[i] = -4.0 * tr_yz_xyyy[i] * tke_0 - 6.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyy[i] * tbe_0 - 6.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_yyz[i] = -4.0 * tr_yz_xyyz[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_yzz[i] = -4.0 * tr_yz_xyzz[i] * tke_0 - 2.0 * tr_yyz_xzz[i] * tke_0 + 4.0 * tr_yyz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 - 2.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyz_zzz[i] = -4.0 * tr_yz_xzzz[i] * tke_0 + 4.0 * tr_yyz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-190 components of targeted buffer : FF

    auto tr_0_0_xy_yzz_xxx = pbuffer.data(idx_op_geom_020_ff + 180);

    auto tr_0_0_xy_yzz_xxy = pbuffer.data(idx_op_geom_020_ff + 181);

    auto tr_0_0_xy_yzz_xxz = pbuffer.data(idx_op_geom_020_ff + 182);

    auto tr_0_0_xy_yzz_xyy = pbuffer.data(idx_op_geom_020_ff + 183);

    auto tr_0_0_xy_yzz_xyz = pbuffer.data(idx_op_geom_020_ff + 184);

    auto tr_0_0_xy_yzz_xzz = pbuffer.data(idx_op_geom_020_ff + 185);

    auto tr_0_0_xy_yzz_yyy = pbuffer.data(idx_op_geom_020_ff + 186);

    auto tr_0_0_xy_yzz_yyz = pbuffer.data(idx_op_geom_020_ff + 187);

    auto tr_0_0_xy_yzz_yzz = pbuffer.data(idx_op_geom_020_ff + 188);

    auto tr_0_0_xy_yzz_zzz = pbuffer.data(idx_op_geom_020_ff + 189);

    #pragma omp simd aligned(tr_0_0_xy_yzz_xxx, tr_0_0_xy_yzz_xxy, tr_0_0_xy_yzz_xxz, tr_0_0_xy_yzz_xyy, tr_0_0_xy_yzz_xyz, tr_0_0_xy_yzz_xzz, tr_0_0_xy_yzz_yyy, tr_0_0_xy_yzz_yyz, tr_0_0_xy_yzz_yzz, tr_0_0_xy_yzz_zzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxy, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_zz_xx, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yzz_xxx[i] = 3.0 * tr_zz_xx[i] - 2.0 * tr_zz_xxxx[i] * tke_0 - 6.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yyzz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_xxy[i] = 2.0 * tr_zz_xy[i] - 2.0 * tr_zz_xxxy[i] * tke_0 + 2.0 * tr_yzz_x[i] - 4.0 * tr_yzz_xyy[i] * tke_0 - 2.0 * tr_yzz_xxx[i] * tke_0 + 4.0 * tr_yzz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xy[i] * tbe_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxy[i] * tbe_0 - 2.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_xxz[i] = 2.0 * tr_zz_xz[i] - 2.0 * tr_zz_xxxz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xz[i] * tbe_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_xyy[i] = tr_zz_yy[i] - 2.0 * tr_zz_xxyy[i] * tke_0 + 2.0 * tr_yzz_y[i] - 2.0 * tr_yzz_yyy[i] * tke_0 - 4.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_yy[i] * tbe_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xyy[i] * tbe_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_xyz[i] = tr_zz_yz[i] - 2.0 * tr_zz_xxyz[i] * tke_0 + tr_yzz_z[i] - 2.0 * tr_yzz_yyz[i] * tke_0 - 2.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_yz[i] * tbe_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xyz[i] * tbe_0 - 2.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_xzz[i] = tr_zz_zz[i] - 2.0 * tr_zz_xxzz[i] * tke_0 - 2.0 * tr_yzz_yzz[i] * tke_0 + 4.0 * tr_yzz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_zz[i] * tbe_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_yyy[i] = -2.0 * tr_zz_xyyy[i] * tke_0 - 6.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yyy[i] * tbe_0 - 6.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_yyz[i] = -2.0 * tr_zz_xyyz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yyz[i] * tbe_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_yzz[i] = -2.0 * tr_zz_xyzz[i] * tke_0 - 2.0 * tr_yzz_xzz[i] * tke_0 + 4.0 * tr_yzz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yzz[i] * tbe_0 - 2.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzz_zzz[i] = -2.0 * tr_zz_xzzz[i] * tke_0 + 4.0 * tr_yzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 190-200 components of targeted buffer : FF

    auto tr_0_0_xy_zzz_xxx = pbuffer.data(idx_op_geom_020_ff + 190);

    auto tr_0_0_xy_zzz_xxy = pbuffer.data(idx_op_geom_020_ff + 191);

    auto tr_0_0_xy_zzz_xxz = pbuffer.data(idx_op_geom_020_ff + 192);

    auto tr_0_0_xy_zzz_xyy = pbuffer.data(idx_op_geom_020_ff + 193);

    auto tr_0_0_xy_zzz_xyz = pbuffer.data(idx_op_geom_020_ff + 194);

    auto tr_0_0_xy_zzz_xzz = pbuffer.data(idx_op_geom_020_ff + 195);

    auto tr_0_0_xy_zzz_yyy = pbuffer.data(idx_op_geom_020_ff + 196);

    auto tr_0_0_xy_zzz_yyz = pbuffer.data(idx_op_geom_020_ff + 197);

    auto tr_0_0_xy_zzz_yzz = pbuffer.data(idx_op_geom_020_ff + 198);

    auto tr_0_0_xy_zzz_zzz = pbuffer.data(idx_op_geom_020_ff + 199);

    #pragma omp simd aligned(tr_0_0_xy_zzz_xxx, tr_0_0_xy_zzz_xxy, tr_0_0_xy_zzz_xxz, tr_0_0_xy_zzz_xyy, tr_0_0_xy_zzz_xyz, tr_0_0_xy_zzz_xzz, tr_0_0_xy_zzz_yyy, tr_0_0_xy_zzz_yyz, tr_0_0_xy_zzz_yzz, tr_0_0_xy_zzz_zzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxy, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxxxy, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_zzz_xxx[i] = -6.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yzzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_xxy[i] = 2.0 * tr_zzz_x[i] - 4.0 * tr_zzz_xyy[i] * tke_0 - 2.0 * tr_zzz_xxx[i] * tke_0 + 4.0 * tr_zzz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xy[i] * tbe_0 + 4.0 * tr_yzzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_xxz[i] = -4.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xz[i] * tbe_0 + 4.0 * tr_yzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_xyy[i] = 2.0 * tr_zzz_y[i] - 2.0 * tr_zzz_yyy[i] * tke_0 - 4.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_yy[i] * tbe_0 + 4.0 * tr_yzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xzzz_xy[i] * tbe_0 + 4.0 * tr_xzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_xyz[i] = tr_zzz_z[i] - 2.0 * tr_zzz_yyz[i] * tke_0 - 2.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_yz[i] * tbe_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xz[i] * tbe_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_xzz[i] = -2.0 * tr_zzz_yzz[i] * tke_0 + 4.0 * tr_zzz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_zz[i] * tbe_0 + 4.0 * tr_yzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_yyy[i] = -6.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_yyz[i] = -4.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xzzz_yz[i] * tbe_0 + 4.0 * tr_xzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_yzz[i] = -2.0 * tr_zzz_xzz[i] * tke_0 + 4.0 * tr_zzz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_zz[i] * tbe_0 + 4.0 * tr_xzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzz_zzz[i] = 4.0 * tr_zzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 200-210 components of targeted buffer : FF

    auto tr_0_0_xz_xxx_xxx = pbuffer.data(idx_op_geom_020_ff + 200);

    auto tr_0_0_xz_xxx_xxy = pbuffer.data(idx_op_geom_020_ff + 201);

    auto tr_0_0_xz_xxx_xxz = pbuffer.data(idx_op_geom_020_ff + 202);

    auto tr_0_0_xz_xxx_xyy = pbuffer.data(idx_op_geom_020_ff + 203);

    auto tr_0_0_xz_xxx_xyz = pbuffer.data(idx_op_geom_020_ff + 204);

    auto tr_0_0_xz_xxx_xzz = pbuffer.data(idx_op_geom_020_ff + 205);

    auto tr_0_0_xz_xxx_yyy = pbuffer.data(idx_op_geom_020_ff + 206);

    auto tr_0_0_xz_xxx_yyz = pbuffer.data(idx_op_geom_020_ff + 207);

    auto tr_0_0_xz_xxx_yzz = pbuffer.data(idx_op_geom_020_ff + 208);

    auto tr_0_0_xz_xxx_zzz = pbuffer.data(idx_op_geom_020_ff + 209);

    #pragma omp simd aligned(tr_0_0_xz_xxx_xxx, tr_0_0_xz_xxx_xxy, tr_0_0_xz_xxx_xxz, tr_0_0_xz_xxx_xyy, tr_0_0_xz_xxx_xyz, tr_0_0_xz_xxx_xzz, tr_0_0_xz_xxx_yyy, tr_0_0_xz_xxx_yyz, tr_0_0_xz_xxx_yzz, tr_0_0_xz_xxx_zzz, tr_xx_xx, tr_xx_xxxz, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xx_zzzz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxxxz, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxx_xx, tr_xxxx_xxxz, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxy, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_zzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxx_xxx[i] = -6.0 * tr_xx_xxxz[i] * tke_0 - 6.0 * tr_xxz_xxx[i] * tbe_0 - 6.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_xxy[i] = -6.0 * tr_xx_xxyz[i] * tke_0 - 6.0 * tr_xxz_xxy[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_xy[i] * tbe_0 + 4.0 * tr_xxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_xxz[i] = 3.0 * tr_xx_xx[i] - 6.0 * tr_xx_xxzz[i] * tke_0 - 6.0 * tr_xxz_xxz[i] * tbe_0 + 2.0 * tr_xxx_x[i] - 4.0 * tr_xxx_xzz[i] * tke_0 - 2.0 * tr_xxx_xxx[i] * tke_0 + 4.0 * tr_xxx_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_xz[i] * tbe_0 + 4.0 * tr_xxxz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xx[i] * tbe_0 + 4.0 * tr_xxxx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_xyy[i] = -6.0 * tr_xx_xyyz[i] * tke_0 - 6.0 * tr_xxz_xyy[i] * tbe_0 - 2.0 * tr_xxx_yyz[i] * tke_0 + 4.0 * tr_xxx_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_yy[i] * tbe_0 + 4.0 * tr_xxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_xyz[i] = 3.0 * tr_xx_xy[i] - 6.0 * tr_xx_xyzz[i] * tke_0 - 6.0 * tr_xxz_xyz[i] * tbe_0 + tr_xxx_y[i] - 2.0 * tr_xxx_yzz[i] * tke_0 - 2.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_yz[i] * tbe_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xy[i] * tbe_0 + 4.0 * tr_xxxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_xzz[i] = 6.0 * tr_xx_xz[i] - 6.0 * tr_xx_xzzz[i] * tke_0 - 6.0 * tr_xxz_xzz[i] * tbe_0 + 2.0 * tr_xxx_z[i] - 2.0 * tr_xxx_zzz[i] * tke_0 - 4.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_zz[i] * tbe_0 + 4.0 * tr_xxxz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxx_xz[i] * tbe_0 + 4.0 * tr_xxxx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_yyy[i] = -6.0 * tr_xx_yyyz[i] * tke_0 - 6.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxx_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_yyz[i] = 3.0 * tr_xx_yy[i] - 6.0 * tr_xx_yyzz[i] * tke_0 - 6.0 * tr_xxz_yyz[i] * tbe_0 - 2.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxx_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yy[i] * tbe_0 + 4.0 * tr_xxxx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_yzz[i] = 6.0 * tr_xx_yz[i] - 6.0 * tr_xx_yzzz[i] * tke_0 - 6.0 * tr_xxz_yzz[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxx_yz[i] * tbe_0 + 4.0 * tr_xxxx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxx_zzz[i] = 9.0 * tr_xx_zz[i] - 6.0 * tr_xx_zzzz[i] * tke_0 - 6.0 * tr_xxz_zzz[i] * tbe_0 - 6.0 * tr_xxx_xzz[i] * tke_0 + 4.0 * tr_xxx_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxx_zz[i] * tbe_0 + 4.0 * tr_xxxx_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-220 components of targeted buffer : FF

    auto tr_0_0_xz_xxy_xxx = pbuffer.data(idx_op_geom_020_ff + 210);

    auto tr_0_0_xz_xxy_xxy = pbuffer.data(idx_op_geom_020_ff + 211);

    auto tr_0_0_xz_xxy_xxz = pbuffer.data(idx_op_geom_020_ff + 212);

    auto tr_0_0_xz_xxy_xyy = pbuffer.data(idx_op_geom_020_ff + 213);

    auto tr_0_0_xz_xxy_xyz = pbuffer.data(idx_op_geom_020_ff + 214);

    auto tr_0_0_xz_xxy_xzz = pbuffer.data(idx_op_geom_020_ff + 215);

    auto tr_0_0_xz_xxy_yyy = pbuffer.data(idx_op_geom_020_ff + 216);

    auto tr_0_0_xz_xxy_yyz = pbuffer.data(idx_op_geom_020_ff + 217);

    auto tr_0_0_xz_xxy_yzz = pbuffer.data(idx_op_geom_020_ff + 218);

    auto tr_0_0_xz_xxy_zzz = pbuffer.data(idx_op_geom_020_ff + 219);

    #pragma omp simd aligned(tr_0_0_xz_xxy_xxx, tr_0_0_xz_xxy_xxy, tr_0_0_xz_xxy_xxz, tr_0_0_xz_xxy_xyy, tr_0_0_xz_xxy_xyz, tr_0_0_xz_xxy_xzz, tr_0_0_xz_xxy_yyy, tr_0_0_xz_xxy_yyz, tr_0_0_xz_xxy_yzz, tr_0_0_xz_xxy_zzz, tr_xxxy_xx, tr_xxxy_xxxz, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxy_xxx[i] = -4.0 * tr_xy_xxxz[i] * tke_0 - 4.0 * tr_xyz_xxx[i] * tbe_0 - 6.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_xxy[i] = -4.0 * tr_xy_xxyz[i] * tke_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_xxz[i] = 2.0 * tr_xy_xx[i] - 4.0 * tr_xy_xxzz[i] * tke_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 + 2.0 * tr_xxy_x[i] - 4.0 * tr_xxy_xzz[i] * tke_0 - 2.0 * tr_xxy_xxx[i] * tke_0 + 4.0 * tr_xxy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_xyy[i] = -4.0 * tr_xy_xyyz[i] * tke_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 - 2.0 * tr_xxy_yyz[i] * tke_0 + 4.0 * tr_xxy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_xyz[i] = 2.0 * tr_xy_xy[i] - 4.0 * tr_xy_xyzz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + tr_xxy_y[i] - 2.0 * tr_xxy_yzz[i] * tke_0 - 2.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xy[i] * tbe_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_xzz[i] = 4.0 * tr_xy_xz[i] - 4.0 * tr_xy_xzzz[i] * tke_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 + 2.0 * tr_xxy_z[i] - 2.0 * tr_xxy_zzz[i] * tke_0 - 4.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xz[i] * tbe_0 + 4.0 * tr_xxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_yyy[i] = -4.0 * tr_xy_yyyz[i] * tke_0 - 4.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xxy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_yyz[i] = 2.0 * tr_xy_yy[i] - 4.0 * tr_xy_yyzz[i] * tke_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 - 2.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yy[i] * tbe_0 + 4.0 * tr_xxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_yzz[i] = 4.0 * tr_xy_yz[i] - 4.0 * tr_xy_yzzz[i] * tke_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yz[i] * tbe_0 + 4.0 * tr_xxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxy_zzz[i] = 6.0 * tr_xy_zz[i] - 4.0 * tr_xy_zzzz[i] * tke_0 - 4.0 * tr_xyz_zzz[i] * tbe_0 - 6.0 * tr_xxy_xzz[i] * tke_0 + 4.0 * tr_xxy_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxy_zz[i] * tbe_0 + 4.0 * tr_xxxy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 220-230 components of targeted buffer : FF

    auto tr_0_0_xz_xxz_xxx = pbuffer.data(idx_op_geom_020_ff + 220);

    auto tr_0_0_xz_xxz_xxy = pbuffer.data(idx_op_geom_020_ff + 221);

    auto tr_0_0_xz_xxz_xxz = pbuffer.data(idx_op_geom_020_ff + 222);

    auto tr_0_0_xz_xxz_xyy = pbuffer.data(idx_op_geom_020_ff + 223);

    auto tr_0_0_xz_xxz_xyz = pbuffer.data(idx_op_geom_020_ff + 224);

    auto tr_0_0_xz_xxz_xzz = pbuffer.data(idx_op_geom_020_ff + 225);

    auto tr_0_0_xz_xxz_yyy = pbuffer.data(idx_op_geom_020_ff + 226);

    auto tr_0_0_xz_xxz_yyz = pbuffer.data(idx_op_geom_020_ff + 227);

    auto tr_0_0_xz_xxz_yzz = pbuffer.data(idx_op_geom_020_ff + 228);

    auto tr_0_0_xz_xxz_zzz = pbuffer.data(idx_op_geom_020_ff + 229);

    #pragma omp simd aligned(tr_0_0_xz_xxz_xxx, tr_0_0_xz_xxz_xxy, tr_0_0_xz_xxz_xxz, tr_0_0_xz_xxz_xyy, tr_0_0_xz_xxz_xyz, tr_0_0_xz_xxz_xzz, tr_0_0_xz_xxz_yyy, tr_0_0_xz_xxz_yyz, tr_0_0_xz_xxz_yzz, tr_0_0_xz_xxz_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xx, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxz_xx, tr_xxxz_xxxz, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxz_xxx[i] = 2.0 * tr_x_xxx[i] - 4.0 * tr_xz_xxxz[i] * tke_0 - 4.0 * tr_xzz_xxx[i] * tbe_0 + 3.0 * tr_xx_xx[i] - 2.0 * tr_xx_xxxx[i] * tke_0 - 6.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xxzz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_xxy[i] = 2.0 * tr_x_xxy[i] - 4.0 * tr_xz_xxyz[i] * tke_0 - 4.0 * tr_xzz_xxy[i] * tbe_0 + 2.0 * tr_xx_xy[i] - 2.0 * tr_xx_xxxy[i] * tke_0 - 4.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_xy[i] * tbe_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_xxz[i] = 2.0 * tr_x_xxz[i] + 2.0 * tr_xz_xx[i] - 4.0 * tr_xz_xxzz[i] * tke_0 - 4.0 * tr_xzz_xxz[i] * tbe_0 + 2.0 * tr_xx_xz[i] - 2.0 * tr_xx_xxxz[i] * tke_0 + 2.0 * tr_xxz_x[i] - 4.0 * tr_xxz_xzz[i] * tke_0 - 2.0 * tr_xxz_xxx[i] * tke_0 + 4.0 * tr_xxz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_xz[i] * tbe_0 + 4.0 * tr_xxzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxz[i] * tbe_0 - 2.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_xyy[i] = 2.0 * tr_x_xyy[i] - 4.0 * tr_xz_xyyz[i] * tke_0 - 4.0 * tr_xzz_xyy[i] * tbe_0 + tr_xx_yy[i] - 2.0 * tr_xx_xxyy[i] * tke_0 - 2.0 * tr_xxz_yyz[i] * tke_0 + 4.0 * tr_xxz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_yy[i] * tbe_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_xyz[i] = 2.0 * tr_x_xyz[i] + 2.0 * tr_xz_xy[i] - 4.0 * tr_xz_xyzz[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tbe_0 + tr_xx_yz[i] - 2.0 * tr_xx_xxyz[i] * tke_0 + tr_xxz_y[i] - 2.0 * tr_xxz_yzz[i] * tke_0 - 2.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_yz[i] * tbe_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyz[i] * tbe_0 - 2.0 * tr_xxxz_xy[i] * tbe_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_xzz[i] = 2.0 * tr_x_xzz[i] + 4.0 * tr_xz_xz[i] - 4.0 * tr_xz_xzzz[i] * tke_0 - 4.0 * tr_xzz_xzz[i] * tbe_0 + tr_xx_zz[i] - 2.0 * tr_xx_xxzz[i] * tke_0 + 2.0 * tr_xxz_z[i] - 2.0 * tr_xxz_zzz[i] * tke_0 - 4.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xzz[i] * tbe_0 - 4.0 * tr_xxxz_xz[i] * tbe_0 + 4.0 * tr_xxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_yyy[i] = 2.0 * tr_x_yyy[i] - 4.0 * tr_xz_yyyz[i] * tke_0 - 4.0 * tr_xzz_yyy[i] * tbe_0 - 2.0 * tr_xx_xyyy[i] * tke_0 + 4.0 * tr_xxz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_yyz[i] = 2.0 * tr_x_yyz[i] + 2.0 * tr_xz_yy[i] - 4.0 * tr_xz_yyzz[i] * tke_0 - 4.0 * tr_xzz_yyz[i] * tbe_0 - 2.0 * tr_xx_xyyz[i] * tke_0 - 2.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyz[i] * tbe_0 - 2.0 * tr_xxxz_yy[i] * tbe_0 + 4.0 * tr_xxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_yzz[i] = 2.0 * tr_x_yzz[i] + 4.0 * tr_xz_yz[i] - 4.0 * tr_xz_yzzz[i] * tke_0 - 4.0 * tr_xzz_yzz[i] * tbe_0 - 2.0 * tr_xx_xyzz[i] * tke_0 - 4.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yzz[i] * tbe_0 - 4.0 * tr_xxxz_yz[i] * tbe_0 + 4.0 * tr_xxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxz_zzz[i] = 2.0 * tr_x_zzz[i] + 6.0 * tr_xz_zz[i] - 4.0 * tr_xz_zzzz[i] * tke_0 - 4.0 * tr_xzz_zzz[i] * tbe_0 - 2.0 * tr_xx_xzzz[i] * tke_0 - 6.0 * tr_xxz_xzz[i] * tke_0 + 4.0 * tr_xxz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zzz[i] * tbe_0 - 6.0 * tr_xxxz_zz[i] * tbe_0 + 4.0 * tr_xxxz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 230-240 components of targeted buffer : FF

    auto tr_0_0_xz_xyy_xxx = pbuffer.data(idx_op_geom_020_ff + 230);

    auto tr_0_0_xz_xyy_xxy = pbuffer.data(idx_op_geom_020_ff + 231);

    auto tr_0_0_xz_xyy_xxz = pbuffer.data(idx_op_geom_020_ff + 232);

    auto tr_0_0_xz_xyy_xyy = pbuffer.data(idx_op_geom_020_ff + 233);

    auto tr_0_0_xz_xyy_xyz = pbuffer.data(idx_op_geom_020_ff + 234);

    auto tr_0_0_xz_xyy_xzz = pbuffer.data(idx_op_geom_020_ff + 235);

    auto tr_0_0_xz_xyy_yyy = pbuffer.data(idx_op_geom_020_ff + 236);

    auto tr_0_0_xz_xyy_yyz = pbuffer.data(idx_op_geom_020_ff + 237);

    auto tr_0_0_xz_xyy_yzz = pbuffer.data(idx_op_geom_020_ff + 238);

    auto tr_0_0_xz_xyy_zzz = pbuffer.data(idx_op_geom_020_ff + 239);

    #pragma omp simd aligned(tr_0_0_xz_xyy_xxx, tr_0_0_xz_xyy_xxy, tr_0_0_xz_xyy_xxz, tr_0_0_xz_xyy_xyy, tr_0_0_xz_xyy_xyz, tr_0_0_xz_xyy_xzz, tr_0_0_xz_xyy_yyy, tr_0_0_xz_xyy_yyz, tr_0_0_xz_xyy_yzz, tr_0_0_xz_xyy_zzz, tr_xxyy_xx, tr_xxyy_xxxz, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yy_xx, tr_yy_xxxz, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyy_xxx[i] = -2.0 * tr_yy_xxxz[i] * tke_0 - 2.0 * tr_yyz_xxx[i] * tbe_0 - 6.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_xxy[i] = -2.0 * tr_yy_xxyz[i] * tke_0 - 2.0 * tr_yyz_xxy[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_xxz[i] = tr_yy_xx[i] - 2.0 * tr_yy_xxzz[i] * tke_0 - 2.0 * tr_yyz_xxz[i] * tbe_0 + 2.0 * tr_xyy_x[i] - 4.0 * tr_xyy_xzz[i] * tke_0 - 2.0 * tr_xyy_xxx[i] * tke_0 + 4.0 * tr_xyy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_xyy[i] = -2.0 * tr_yy_xyyz[i] * tke_0 - 2.0 * tr_yyz_xyy[i] * tbe_0 - 2.0 * tr_xyy_yyz[i] * tke_0 + 4.0 * tr_xyy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_xyz[i] = tr_yy_xy[i] - 2.0 * tr_yy_xyzz[i] * tke_0 - 2.0 * tr_yyz_xyz[i] * tbe_0 + tr_xyy_y[i] - 2.0 * tr_xyy_yzz[i] * tke_0 - 2.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xy[i] * tbe_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_xzz[i] = 2.0 * tr_yy_xz[i] - 2.0 * tr_yy_xzzz[i] * tke_0 - 2.0 * tr_yyz_xzz[i] * tbe_0 + 2.0 * tr_xyy_z[i] - 2.0 * tr_xyy_zzz[i] * tke_0 - 4.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyy_xz[i] * tbe_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_yyy[i] = -2.0 * tr_yy_yyyz[i] * tke_0 - 2.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_xyy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_yyz[i] = tr_yy_yy[i] - 2.0 * tr_yy_yyzz[i] * tke_0 - 2.0 * tr_yyz_yyz[i] * tbe_0 - 2.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yy[i] * tbe_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_yzz[i] = 2.0 * tr_yy_yz[i] - 2.0 * tr_yy_yzzz[i] * tke_0 - 2.0 * tr_yyz_yzz[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyy_yz[i] * tbe_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyy_zzz[i] = 3.0 * tr_yy_zz[i] - 2.0 * tr_yy_zzzz[i] * tke_0 - 2.0 * tr_yyz_zzz[i] * tbe_0 - 6.0 * tr_xyy_xzz[i] * tke_0 + 4.0 * tr_xyy_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_zz[i] * tbe_0 + 4.0 * tr_xxyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-250 components of targeted buffer : FF

    auto tr_0_0_xz_xyz_xxx = pbuffer.data(idx_op_geom_020_ff + 240);

    auto tr_0_0_xz_xyz_xxy = pbuffer.data(idx_op_geom_020_ff + 241);

    auto tr_0_0_xz_xyz_xxz = pbuffer.data(idx_op_geom_020_ff + 242);

    auto tr_0_0_xz_xyz_xyy = pbuffer.data(idx_op_geom_020_ff + 243);

    auto tr_0_0_xz_xyz_xyz = pbuffer.data(idx_op_geom_020_ff + 244);

    auto tr_0_0_xz_xyz_xzz = pbuffer.data(idx_op_geom_020_ff + 245);

    auto tr_0_0_xz_xyz_yyy = pbuffer.data(idx_op_geom_020_ff + 246);

    auto tr_0_0_xz_xyz_yyz = pbuffer.data(idx_op_geom_020_ff + 247);

    auto tr_0_0_xz_xyz_yzz = pbuffer.data(idx_op_geom_020_ff + 248);

    auto tr_0_0_xz_xyz_zzz = pbuffer.data(idx_op_geom_020_ff + 249);

    #pragma omp simd aligned(tr_0_0_xz_xyz_xxx, tr_0_0_xz_xyz_xxy, tr_0_0_xz_xyz_xxz, tr_0_0_xz_xyz_xyy, tr_0_0_xz_xyz_xyz, tr_0_0_xz_xyz_xzz, tr_0_0_xz_xyz_yyy, tr_0_0_xz_xyz_yyz, tr_0_0_xz_xyz_yzz, tr_0_0_xz_xyz_zzz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xxxz, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyz_xxx[i] = tr_y_xxx[i] - 2.0 * tr_yz_xxxz[i] * tke_0 - 2.0 * tr_yzz_xxx[i] * tbe_0 + 3.0 * tr_xy_xx[i] - 2.0 * tr_xy_xxxx[i] * tke_0 - 6.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_xxy[i] = tr_y_xxy[i] - 2.0 * tr_yz_xxyz[i] * tke_0 - 2.0 * tr_yzz_xxy[i] * tbe_0 + 2.0 * tr_xy_xy[i] - 2.0 * tr_xy_xxxy[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_xxz[i] = tr_y_xxz[i] + tr_yz_xx[i] - 2.0 * tr_yz_xxzz[i] * tke_0 - 2.0 * tr_yzz_xxz[i] * tbe_0 + 2.0 * tr_xy_xz[i] - 2.0 * tr_xy_xxxz[i] * tke_0 + 2.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xzz[i] * tke_0 - 2.0 * tr_xyz_xxx[i] * tke_0 + 4.0 * tr_xyz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxz[i] * tbe_0 - 2.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_xyy[i] = tr_y_xyy[i] - 2.0 * tr_yz_xyyz[i] * tke_0 - 2.0 * tr_yzz_xyy[i] * tbe_0 + tr_xy_yy[i] - 2.0 * tr_xy_xxyy[i] * tke_0 - 2.0 * tr_xyz_yyz[i] * tke_0 + 4.0 * tr_xyz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_xyz[i] = tr_y_xyz[i] + tr_yz_xy[i] - 2.0 * tr_yz_xyzz[i] * tke_0 - 2.0 * tr_yzz_xyz[i] * tbe_0 + tr_xy_yz[i] - 2.0 * tr_xy_xxyz[i] * tke_0 + tr_xyz_y[i] - 2.0 * tr_xyz_yzz[i] * tke_0 - 2.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyz[i] * tbe_0 - 2.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_xzz[i] = tr_y_xzz[i] + 2.0 * tr_yz_xz[i] - 2.0 * tr_yz_xzzz[i] * tke_0 - 2.0 * tr_yzz_xzz[i] * tbe_0 + tr_xy_zz[i] - 2.0 * tr_xy_xxzz[i] * tke_0 + 2.0 * tr_xyz_z[i] - 2.0 * tr_xyz_zzz[i] * tke_0 - 4.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xzz[i] * tbe_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_yyy[i] = tr_y_yyy[i] - 2.0 * tr_yz_yyyz[i] * tke_0 - 2.0 * tr_yzz_yyy[i] * tbe_0 - 2.0 * tr_xy_xyyy[i] * tke_0 + 4.0 * tr_xyz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_yyz[i] = tr_y_yyz[i] + tr_yz_yy[i] - 2.0 * tr_yz_yyzz[i] * tke_0 - 2.0 * tr_yzz_yyz[i] * tbe_0 - 2.0 * tr_xy_xyyz[i] * tke_0 - 2.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyz[i] * tbe_0 - 2.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_yzz[i] = tr_y_yzz[i] + 2.0 * tr_yz_yz[i] - 2.0 * tr_yz_yzzz[i] * tke_0 - 2.0 * tr_yzz_yzz[i] * tbe_0 - 2.0 * tr_xy_xyzz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yzz[i] * tbe_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyz_zzz[i] = tr_y_zzz[i] + 3.0 * tr_yz_zz[i] - 2.0 * tr_yz_zzzz[i] * tke_0 - 2.0 * tr_yzz_zzz[i] * tbe_0 - 2.0 * tr_xy_xzzz[i] * tke_0 - 6.0 * tr_xyz_xzz[i] * tke_0 + 4.0 * tr_xyz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_zzz[i] * tbe_0 - 6.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 250-260 components of targeted buffer : FF

    auto tr_0_0_xz_xzz_xxx = pbuffer.data(idx_op_geom_020_ff + 250);

    auto tr_0_0_xz_xzz_xxy = pbuffer.data(idx_op_geom_020_ff + 251);

    auto tr_0_0_xz_xzz_xxz = pbuffer.data(idx_op_geom_020_ff + 252);

    auto tr_0_0_xz_xzz_xyy = pbuffer.data(idx_op_geom_020_ff + 253);

    auto tr_0_0_xz_xzz_xyz = pbuffer.data(idx_op_geom_020_ff + 254);

    auto tr_0_0_xz_xzz_xzz = pbuffer.data(idx_op_geom_020_ff + 255);

    auto tr_0_0_xz_xzz_yyy = pbuffer.data(idx_op_geom_020_ff + 256);

    auto tr_0_0_xz_xzz_yyz = pbuffer.data(idx_op_geom_020_ff + 257);

    auto tr_0_0_xz_xzz_yzz = pbuffer.data(idx_op_geom_020_ff + 258);

    auto tr_0_0_xz_xzz_zzz = pbuffer.data(idx_op_geom_020_ff + 259);

    #pragma omp simd aligned(tr_0_0_xz_xzz_xxx, tr_0_0_xz_xzz_xxy, tr_0_0_xz_xzz_xxz, tr_0_0_xz_xzz_xyy, tr_0_0_xz_xzz_xyz, tr_0_0_xz_xzz_xzz, tr_0_0_xz_xzz_yyy, tr_0_0_xz_xzz_yyz, tr_0_0_xz_xzz_yzz, tr_0_0_xz_xzz_zzz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xxxz, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, tr_zz_xx, tr_zz_xxxz, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xzz_xxx[i] = 2.0 * tr_z_xxx[i] - 2.0 * tr_zz_xxxz[i] * tke_0 - 2.0 * tr_zzz_xxx[i] * tbe_0 + 6.0 * tr_xz_xx[i] - 4.0 * tr_xz_xxxx[i] * tke_0 - 6.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xzzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_xxy[i] = 2.0 * tr_z_xxy[i] - 2.0 * tr_zz_xxyz[i] * tke_0 - 2.0 * tr_zzz_xxy[i] * tbe_0 + 4.0 * tr_xz_xy[i] - 4.0 * tr_xz_xxxy[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_xy[i] * tbe_0 + 4.0 * tr_xzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_xxz[i] = 2.0 * tr_z_xxz[i] + tr_zz_xx[i] - 2.0 * tr_zz_xxzz[i] * tke_0 - 2.0 * tr_zzz_xxz[i] * tbe_0 + 4.0 * tr_xz_xz[i] - 4.0 * tr_xz_xxxz[i] * tke_0 + 2.0 * tr_xzz_x[i] - 4.0 * tr_xzz_xzz[i] * tke_0 - 2.0 * tr_xzz_xxx[i] * tke_0 + 4.0 * tr_xzz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_xz[i] * tbe_0 + 4.0 * tr_xzzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxz[i] * tbe_0 - 2.0 * tr_xxzz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_xyy[i] = 2.0 * tr_z_xyy[i] - 2.0 * tr_zz_xyyz[i] * tke_0 - 2.0 * tr_zzz_xyy[i] * tbe_0 + 2.0 * tr_xz_yy[i] - 4.0 * tr_xz_xxyy[i] * tke_0 - 2.0 * tr_xzz_yyz[i] * tke_0 + 4.0 * tr_xzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_xyz[i] = 2.0 * tr_z_xyz[i] + tr_zz_xy[i] - 2.0 * tr_zz_xyzz[i] * tke_0 - 2.0 * tr_zzz_xyz[i] * tbe_0 + 2.0 * tr_xz_yz[i] - 4.0 * tr_xz_xxyz[i] * tke_0 + tr_xzz_y[i] - 2.0 * tr_xzz_yzz[i] * tke_0 - 2.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_yz[i] * tbe_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xyz[i] * tbe_0 - 2.0 * tr_xxzz_xy[i] * tbe_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_xzz[i] = 2.0 * tr_z_xzz[i] + 2.0 * tr_zz_xz[i] - 2.0 * tr_zz_xzzz[i] * tke_0 - 2.0 * tr_zzz_xzz[i] * tbe_0 + 2.0 * tr_xz_zz[i] - 4.0 * tr_xz_xxzz[i] * tke_0 + 2.0 * tr_xzz_z[i] - 2.0 * tr_xzz_zzz[i] * tke_0 - 4.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_zz[i] * tbe_0 + 4.0 * tr_xzzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xzz[i] * tbe_0 - 4.0 * tr_xxzz_xz[i] * tbe_0 + 4.0 * tr_xxzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_yyy[i] = 2.0 * tr_z_yyy[i] - 2.0 * tr_zz_yyyz[i] * tke_0 - 2.0 * tr_zzz_yyy[i] * tbe_0 - 4.0 * tr_xz_xyyy[i] * tke_0 + 4.0 * tr_xzz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_yyz[i] = 2.0 * tr_z_yyz[i] + tr_zz_yy[i] - 2.0 * tr_zz_yyzz[i] * tke_0 - 2.0 * tr_zzz_yyz[i] * tbe_0 - 4.0 * tr_xz_xyyz[i] * tke_0 - 2.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yyz[i] * tbe_0 - 2.0 * tr_xxzz_yy[i] * tbe_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_yzz[i] = 2.0 * tr_z_yzz[i] + 2.0 * tr_zz_yz[i] - 2.0 * tr_zz_yzzz[i] * tke_0 - 2.0 * tr_zzz_yzz[i] * tbe_0 - 4.0 * tr_xz_xyzz[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yzz[i] * tbe_0 - 4.0 * tr_xxzz_yz[i] * tbe_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzz_zzz[i] = 2.0 * tr_z_zzz[i] + 3.0 * tr_zz_zz[i] - 2.0 * tr_zz_zzzz[i] * tke_0 - 2.0 * tr_zzz_zzz[i] * tbe_0 - 4.0 * tr_xz_xzzz[i] * tke_0 - 6.0 * tr_xzz_xzz[i] * tke_0 + 4.0 * tr_xzz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_zzz[i] * tbe_0 - 6.0 * tr_xxzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 260-270 components of targeted buffer : FF

    auto tr_0_0_xz_yyy_xxx = pbuffer.data(idx_op_geom_020_ff + 260);

    auto tr_0_0_xz_yyy_xxy = pbuffer.data(idx_op_geom_020_ff + 261);

    auto tr_0_0_xz_yyy_xxz = pbuffer.data(idx_op_geom_020_ff + 262);

    auto tr_0_0_xz_yyy_xyy = pbuffer.data(idx_op_geom_020_ff + 263);

    auto tr_0_0_xz_yyy_xyz = pbuffer.data(idx_op_geom_020_ff + 264);

    auto tr_0_0_xz_yyy_xzz = pbuffer.data(idx_op_geom_020_ff + 265);

    auto tr_0_0_xz_yyy_yyy = pbuffer.data(idx_op_geom_020_ff + 266);

    auto tr_0_0_xz_yyy_yyz = pbuffer.data(idx_op_geom_020_ff + 267);

    auto tr_0_0_xz_yyy_yzz = pbuffer.data(idx_op_geom_020_ff + 268);

    auto tr_0_0_xz_yyy_zzz = pbuffer.data(idx_op_geom_020_ff + 269);

    #pragma omp simd aligned(tr_0_0_xz_yyy_xxx, tr_0_0_xz_yyy_xxy, tr_0_0_xz_yyy_xxz, tr_0_0_xz_yyy_xyy, tr_0_0_xz_yyy_xyz, tr_0_0_xz_yyy_xzz, tr_0_0_xz_yyy_yyy, tr_0_0_xz_yyy_yyz, tr_0_0_xz_yyy_yzz, tr_0_0_xz_yyy_zzz, tr_xyyy_xx, tr_xyyy_xxxz, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxxxz, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyy_xxx[i] = -6.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_xxy[i] = -4.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xy[i] * tbe_0 + 4.0 * tr_yyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_xxz[i] = 2.0 * tr_yyy_x[i] - 4.0 * tr_yyy_xzz[i] * tke_0 - 2.0 * tr_yyy_xxx[i] * tke_0 + 4.0 * tr_yyy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xz[i] * tbe_0 + 4.0 * tr_yyyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_xyy[i] = -2.0 * tr_yyy_yyz[i] * tke_0 + 4.0 * tr_yyy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_yy[i] * tbe_0 + 4.0 * tr_yyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_xyz[i] = tr_yyy_y[i] - 2.0 * tr_yyy_yzz[i] * tke_0 - 2.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_yz[i] * tbe_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xy[i] * tbe_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_xzz[i] = 2.0 * tr_yyy_z[i] - 2.0 * tr_yyy_zzz[i] * tke_0 - 4.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_zz[i] * tbe_0 + 4.0 * tr_yyyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyy_xz[i] * tbe_0 + 4.0 * tr_xyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_yyy[i] = 4.0 * tr_yyy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_yyz[i] = -2.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yy[i] * tbe_0 + 4.0 * tr_xyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_yzz[i] = -4.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyy_yz[i] * tbe_0 + 4.0 * tr_xyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyy_zzz[i] = -6.0 * tr_yyy_xzz[i] * tke_0 + 4.0 * tr_yyy_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-280 components of targeted buffer : FF

    auto tr_0_0_xz_yyz_xxx = pbuffer.data(idx_op_geom_020_ff + 270);

    auto tr_0_0_xz_yyz_xxy = pbuffer.data(idx_op_geom_020_ff + 271);

    auto tr_0_0_xz_yyz_xxz = pbuffer.data(idx_op_geom_020_ff + 272);

    auto tr_0_0_xz_yyz_xyy = pbuffer.data(idx_op_geom_020_ff + 273);

    auto tr_0_0_xz_yyz_xyz = pbuffer.data(idx_op_geom_020_ff + 274);

    auto tr_0_0_xz_yyz_xzz = pbuffer.data(idx_op_geom_020_ff + 275);

    auto tr_0_0_xz_yyz_yyy = pbuffer.data(idx_op_geom_020_ff + 276);

    auto tr_0_0_xz_yyz_yyz = pbuffer.data(idx_op_geom_020_ff + 277);

    auto tr_0_0_xz_yyz_yzz = pbuffer.data(idx_op_geom_020_ff + 278);

    auto tr_0_0_xz_yyz_zzz = pbuffer.data(idx_op_geom_020_ff + 279);

    #pragma omp simd aligned(tr_0_0_xz_yyz_xxx, tr_0_0_xz_yyz_xxy, tr_0_0_xz_yyz_xxz, tr_0_0_xz_yyz_xyy, tr_0_0_xz_yyz_xyz, tr_0_0_xz_yyz_xzz, tr_0_0_xz_yyz_yyy, tr_0_0_xz_yyz_yyz, tr_0_0_xz_yyz_yzz, tr_0_0_xz_yyz_zzz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xxxz, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_yy_xx, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyz_xxx[i] = 3.0 * tr_yy_xx[i] - 2.0 * tr_yy_xxxx[i] * tke_0 - 6.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_yyzz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_xxy[i] = 2.0 * tr_yy_xy[i] - 2.0 * tr_yy_xxxy[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xy[i] * tbe_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_xxz[i] = 2.0 * tr_yy_xz[i] - 2.0 * tr_yy_xxxz[i] * tke_0 + 2.0 * tr_yyz_x[i] - 4.0 * tr_yyz_xzz[i] * tke_0 - 2.0 * tr_yyz_xxx[i] * tke_0 + 4.0 * tr_yyz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xz[i] * tbe_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxz[i] * tbe_0 - 2.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_xyy[i] = tr_yy_yy[i] - 2.0 * tr_yy_xxyy[i] * tke_0 - 2.0 * tr_yyz_yyz[i] * tke_0 + 4.0 * tr_yyz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_yy[i] * tbe_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_xyz[i] = tr_yy_yz[i] - 2.0 * tr_yy_xxyz[i] * tke_0 + tr_yyz_y[i] - 2.0 * tr_yyz_yzz[i] * tke_0 - 2.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_yz[i] * tbe_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyz[i] * tbe_0 - 2.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_xzz[i] = tr_yy_zz[i] - 2.0 * tr_yy_xxzz[i] * tke_0 + 2.0 * tr_yyz_z[i] - 2.0 * tr_yyz_zzz[i] * tke_0 - 4.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_zz[i] * tbe_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xzz[i] * tbe_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_yyy[i] = -2.0 * tr_yy_xyyy[i] * tke_0 + 4.0 * tr_yyz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_yyz[i] = -2.0 * tr_yy_xyyz[i] * tke_0 - 2.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyz[i] * tbe_0 - 2.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_yzz[i] = -2.0 * tr_yy_xyzz[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yzz[i] * tbe_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyz_zzz[i] = -2.0 * tr_yy_xzzz[i] * tke_0 - 6.0 * tr_yyz_xzz[i] * tke_0 + 4.0 * tr_yyz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_zzz[i] * tbe_0 - 6.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 280-290 components of targeted buffer : FF

    auto tr_0_0_xz_yzz_xxx = pbuffer.data(idx_op_geom_020_ff + 280);

    auto tr_0_0_xz_yzz_xxy = pbuffer.data(idx_op_geom_020_ff + 281);

    auto tr_0_0_xz_yzz_xxz = pbuffer.data(idx_op_geom_020_ff + 282);

    auto tr_0_0_xz_yzz_xyy = pbuffer.data(idx_op_geom_020_ff + 283);

    auto tr_0_0_xz_yzz_xyz = pbuffer.data(idx_op_geom_020_ff + 284);

    auto tr_0_0_xz_yzz_xzz = pbuffer.data(idx_op_geom_020_ff + 285);

    auto tr_0_0_xz_yzz_yyy = pbuffer.data(idx_op_geom_020_ff + 286);

    auto tr_0_0_xz_yzz_yyz = pbuffer.data(idx_op_geom_020_ff + 287);

    auto tr_0_0_xz_yzz_yzz = pbuffer.data(idx_op_geom_020_ff + 288);

    auto tr_0_0_xz_yzz_zzz = pbuffer.data(idx_op_geom_020_ff + 289);

    #pragma omp simd aligned(tr_0_0_xz_yzz_xxx, tr_0_0_xz_yzz_xxy, tr_0_0_xz_yzz_xxz, tr_0_0_xz_yzz_xyy, tr_0_0_xz_yzz_xyz, tr_0_0_xz_yzz_xzz, tr_0_0_xz_yzz_yyy, tr_0_0_xz_yzz_yyz, tr_0_0_xz_yzz_yzz, tr_0_0_xz_yzz_zzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxz, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yzz_xxx[i] = 6.0 * tr_yz_xx[i] - 4.0 * tr_yz_xxxx[i] * tke_0 - 6.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_yzzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_xxy[i] = 4.0 * tr_yz_xy[i] - 4.0 * tr_yz_xxxy[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xy[i] * tbe_0 + 4.0 * tr_yzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_xxz[i] = 4.0 * tr_yz_xz[i] - 4.0 * tr_yz_xxxz[i] * tke_0 + 2.0 * tr_yzz_x[i] - 4.0 * tr_yzz_xzz[i] * tke_0 - 2.0 * tr_yzz_xxx[i] * tke_0 + 4.0 * tr_yzz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xz[i] * tbe_0 + 4.0 * tr_yzzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 - 2.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_xyy[i] = 2.0 * tr_yz_yy[i] - 4.0 * tr_yz_xxyy[i] * tke_0 - 2.0 * tr_yzz_yyz[i] * tke_0 + 4.0 * tr_yzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_yy[i] * tbe_0 + 4.0 * tr_yzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_xyz[i] = 2.0 * tr_yz_yz[i] - 4.0 * tr_yz_xxyz[i] * tke_0 + tr_yzz_y[i] - 2.0 * tr_yzz_yzz[i] * tke_0 - 2.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_yz[i] * tbe_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 - 2.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_xzz[i] = 2.0 * tr_yz_zz[i] - 4.0 * tr_yz_xxzz[i] * tke_0 + 2.0 * tr_yzz_z[i] - 2.0 * tr_yzz_zzz[i] * tke_0 - 4.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_zz[i] * tbe_0 + 4.0 * tr_yzzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_yyy[i] = -4.0 * tr_yz_xyyy[i] * tke_0 + 4.0 * tr_yzz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_yyz[i] = -4.0 * tr_yz_xyyz[i] * tke_0 - 2.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 - 2.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_yzz[i] = -4.0 * tr_yz_xyzz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzz_zzz[i] = -4.0 * tr_yz_xzzz[i] * tke_0 - 6.0 * tr_yzz_xzz[i] * tke_0 + 4.0 * tr_yzz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zzz[i] * tbe_0 - 6.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 290-300 components of targeted buffer : FF

    auto tr_0_0_xz_zzz_xxx = pbuffer.data(idx_op_geom_020_ff + 290);

    auto tr_0_0_xz_zzz_xxy = pbuffer.data(idx_op_geom_020_ff + 291);

    auto tr_0_0_xz_zzz_xxz = pbuffer.data(idx_op_geom_020_ff + 292);

    auto tr_0_0_xz_zzz_xyy = pbuffer.data(idx_op_geom_020_ff + 293);

    auto tr_0_0_xz_zzz_xyz = pbuffer.data(idx_op_geom_020_ff + 294);

    auto tr_0_0_xz_zzz_xzz = pbuffer.data(idx_op_geom_020_ff + 295);

    auto tr_0_0_xz_zzz_yyy = pbuffer.data(idx_op_geom_020_ff + 296);

    auto tr_0_0_xz_zzz_yyz = pbuffer.data(idx_op_geom_020_ff + 297);

    auto tr_0_0_xz_zzz_yzz = pbuffer.data(idx_op_geom_020_ff + 298);

    auto tr_0_0_xz_zzz_zzz = pbuffer.data(idx_op_geom_020_ff + 299);

    #pragma omp simd aligned(tr_0_0_xz_zzz_xxx, tr_0_0_xz_zzz_xxy, tr_0_0_xz_zzz_xxz, tr_0_0_xz_zzz_xyy, tr_0_0_xz_zzz_xyz, tr_0_0_xz_zzz_xzz, tr_0_0_xz_zzz_yyy, tr_0_0_xz_zzz_yyz, tr_0_0_xz_zzz_yzz, tr_0_0_xz_zzz_zzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xxxz, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxy, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_zzz, tr_zz_xx, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxxxz, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_zzz_xxx[i] = 9.0 * tr_zz_xx[i] - 6.0 * tr_zz_xxxx[i] * tke_0 - 6.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_zzzz_xx[i] * tbe_0 + 4.0 * tr_zzzz_xxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_xxy[i] = 6.0 * tr_zz_xy[i] - 6.0 * tr_zz_xxxy[i] * tke_0 - 4.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_xy[i] * tbe_0 + 4.0 * tr_zzzz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_xxz[i] = 6.0 * tr_zz_xz[i] - 6.0 * tr_zz_xxxz[i] * tke_0 + 2.0 * tr_zzz_x[i] - 4.0 * tr_zzz_xzz[i] * tke_0 - 2.0 * tr_zzz_xxx[i] * tke_0 + 4.0 * tr_zzz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_xz[i] * tbe_0 + 4.0 * tr_zzzz_xxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxz[i] * tbe_0 - 2.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_xyy[i] = 3.0 * tr_zz_yy[i] - 6.0 * tr_zz_xxyy[i] * tke_0 - 2.0 * tr_zzz_yyz[i] * tke_0 + 4.0 * tr_zzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_yy[i] * tbe_0 + 4.0 * tr_zzzz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_xyz[i] = 3.0 * tr_zz_yz[i] - 6.0 * tr_zz_xxyz[i] * tke_0 + tr_zzz_y[i] - 2.0 * tr_zzz_yzz[i] * tke_0 - 2.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_yz[i] * tbe_0 + 4.0 * tr_zzzz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xyz[i] * tbe_0 - 2.0 * tr_xzzz_xy[i] * tbe_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_xzz[i] = 3.0 * tr_zz_zz[i] - 6.0 * tr_zz_xxzz[i] * tke_0 + 2.0 * tr_zzz_z[i] - 2.0 * tr_zzz_zzz[i] * tke_0 - 4.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_zz[i] * tbe_0 + 4.0 * tr_zzzz_xxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xzz[i] * tbe_0 - 4.0 * tr_xzzz_xz[i] * tbe_0 + 4.0 * tr_xzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_yyy[i] = -6.0 * tr_zz_xyyy[i] * tke_0 + 4.0 * tr_zzz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_yyz[i] = -6.0 * tr_zz_xyyz[i] * tke_0 - 2.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yyz[i] * tbe_0 - 2.0 * tr_xzzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_yzz[i] = -6.0 * tr_zz_xyzz[i] * tke_0 - 4.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yzz[i] * tbe_0 - 4.0 * tr_xzzz_yz[i] * tbe_0 + 4.0 * tr_xzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzz_zzz[i] = -6.0 * tr_zz_xzzz[i] * tke_0 - 6.0 * tr_zzz_xzz[i] * tke_0 + 4.0 * tr_zzz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_zzz[i] * tbe_0 - 6.0 * tr_xzzz_zz[i] * tbe_0 + 4.0 * tr_xzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-310 components of targeted buffer : FF

    auto tr_0_0_yy_xxx_xxx = pbuffer.data(idx_op_geom_020_ff + 300);

    auto tr_0_0_yy_xxx_xxy = pbuffer.data(idx_op_geom_020_ff + 301);

    auto tr_0_0_yy_xxx_xxz = pbuffer.data(idx_op_geom_020_ff + 302);

    auto tr_0_0_yy_xxx_xyy = pbuffer.data(idx_op_geom_020_ff + 303);

    auto tr_0_0_yy_xxx_xyz = pbuffer.data(idx_op_geom_020_ff + 304);

    auto tr_0_0_yy_xxx_xzz = pbuffer.data(idx_op_geom_020_ff + 305);

    auto tr_0_0_yy_xxx_yyy = pbuffer.data(idx_op_geom_020_ff + 306);

    auto tr_0_0_yy_xxx_yyz = pbuffer.data(idx_op_geom_020_ff + 307);

    auto tr_0_0_yy_xxx_yzz = pbuffer.data(idx_op_geom_020_ff + 308);

    auto tr_0_0_yy_xxx_zzz = pbuffer.data(idx_op_geom_020_ff + 309);

    #pragma omp simd aligned(tr_0_0_yy_xxx_xxx, tr_0_0_yy_xxx_xxy, tr_0_0_yy_xxx_xxz, tr_0_0_yy_xxx_xyy, tr_0_0_yy_xxx_xyz, tr_0_0_yy_xxx_xzz, tr_0_0_yy_xxx_yyy, tr_0_0_yy_xxx_yyz, tr_0_0_yy_xxx_yzz, tr_0_0_yy_xxx_zzz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxxyy, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyyyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxy_xx, tr_xxxy_xxxy, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxx_xxx[i] = -2.0 * tr_xxx_xxx[i] * tbe_0 - 2.0 * tr_xxx_xxx[i] * tke_0 + 4.0 * tr_xxx_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_xxy[i] = -2.0 * tr_xxx_xxy[i] * tbe_0 - 6.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_xx[i] * tbe_0 + 8.0 * tr_xxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_xxz[i] = -2.0 * tr_xxx_xxz[i] * tbe_0 - 2.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_xyy[i] = 2.0 * tr_xxx_x[i] - 2.0 * tr_xxx_xyy[i] * tbe_0 - 10.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxx_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxy_xy[i] * tbe_0 + 8.0 * tr_xxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_xyz[i] = -2.0 * tr_xxx_xyz[i] * tbe_0 - 6.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_xz[i] * tbe_0 + 8.0 * tr_xxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_xzz[i] = -2.0 * tr_xxx_xzz[i] * tbe_0 - 2.0 * tr_xxx_xzz[i] * tke_0 + 4.0 * tr_xxx_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_yyy[i] = 6.0 * tr_xxx_y[i] - 2.0 * tr_xxx_yyy[i] * tbe_0 - 14.0 * tr_xxx_yyy[i] * tke_0 + 4.0 * tr_xxx_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxxy_yy[i] * tbe_0 + 8.0 * tr_xxxy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_yyz[i] = 2.0 * tr_xxx_z[i] - 2.0 * tr_xxx_yyz[i] * tbe_0 - 10.0 * tr_xxx_yyz[i] * tke_0 + 4.0 * tr_xxx_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxy_yz[i] * tbe_0 + 8.0 * tr_xxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_yzz[i] = -2.0 * tr_xxx_yzz[i] * tbe_0 - 6.0 * tr_xxx_yzz[i] * tke_0 + 4.0 * tr_xxx_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxy_zz[i] * tbe_0 + 8.0 * tr_xxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxx_zzz[i] = -2.0 * tr_xxx_zzz[i] * tbe_0 - 2.0 * tr_xxx_zzz[i] * tke_0 + 4.0 * tr_xxx_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 310-320 components of targeted buffer : FF

    auto tr_0_0_yy_xxy_xxx = pbuffer.data(idx_op_geom_020_ff + 310);

    auto tr_0_0_yy_xxy_xxy = pbuffer.data(idx_op_geom_020_ff + 311);

    auto tr_0_0_yy_xxy_xxz = pbuffer.data(idx_op_geom_020_ff + 312);

    auto tr_0_0_yy_xxy_xyy = pbuffer.data(idx_op_geom_020_ff + 313);

    auto tr_0_0_yy_xxy_xyz = pbuffer.data(idx_op_geom_020_ff + 314);

    auto tr_0_0_yy_xxy_xzz = pbuffer.data(idx_op_geom_020_ff + 315);

    auto tr_0_0_yy_xxy_yyy = pbuffer.data(idx_op_geom_020_ff + 316);

    auto tr_0_0_yy_xxy_yyz = pbuffer.data(idx_op_geom_020_ff + 317);

    auto tr_0_0_yy_xxy_yzz = pbuffer.data(idx_op_geom_020_ff + 318);

    auto tr_0_0_yy_xxy_zzz = pbuffer.data(idx_op_geom_020_ff + 319);

    #pragma omp simd aligned(tr_0_0_yy_xxy_xxx, tr_0_0_yy_xxy_xxy, tr_0_0_yy_xxy_xxz, tr_0_0_yy_xxy_xyy, tr_0_0_yy_xxy_xyz, tr_0_0_yy_xxy_xzz, tr_0_0_yy_xxy_yyy, tr_0_0_yy_xxy_yyz, tr_0_0_yy_xxy_yzz, tr_0_0_yy_xxy_zzz, tr_xx_xx, tr_xx_xxxy, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxxyy, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyy_xx, tr_xxyy_xxxy, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxy_xxx[i] = -4.0 * tr_xx_xxxy[i] * tke_0 - 6.0 * tr_xxy_xxx[i] * tbe_0 - 2.0 * tr_xxy_xxx[i] * tke_0 + 4.0 * tr_xxy_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_xxy[i] = 2.0 * tr_xx_xx[i] - 4.0 * tr_xx_xxyy[i] * tke_0 - 6.0 * tr_xxy_xxy[i] * tbe_0 - 6.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_xx[i] * tbe_0 + 8.0 * tr_xxyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_xxz[i] = -4.0 * tr_xx_xxyz[i] * tke_0 - 6.0 * tr_xxy_xxz[i] * tbe_0 - 2.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_xyy[i] = 4.0 * tr_xx_xy[i] - 4.0 * tr_xx_xyyy[i] * tke_0 + 2.0 * tr_xxy_x[i] - 6.0 * tr_xxy_xyy[i] * tbe_0 - 10.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxy_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyy_xy[i] * tbe_0 + 8.0 * tr_xxyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_xyz[i] = 2.0 * tr_xx_xz[i] - 4.0 * tr_xx_xyyz[i] * tke_0 - 6.0 * tr_xxy_xyz[i] * tbe_0 - 6.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_xz[i] * tbe_0 + 8.0 * tr_xxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_xzz[i] = -4.0 * tr_xx_xyzz[i] * tke_0 - 6.0 * tr_xxy_xzz[i] * tbe_0 - 2.0 * tr_xxy_xzz[i] * tke_0 + 4.0 * tr_xxy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_yyy[i] = 6.0 * tr_xx_yy[i] - 4.0 * tr_xx_yyyy[i] * tke_0 + 6.0 * tr_xxy_y[i] - 6.0 * tr_xxy_yyy[i] * tbe_0 - 14.0 * tr_xxy_yyy[i] * tke_0 + 4.0 * tr_xxy_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxyy_yy[i] * tbe_0 + 8.0 * tr_xxyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_yyz[i] = 4.0 * tr_xx_yz[i] - 4.0 * tr_xx_yyyz[i] * tke_0 + 2.0 * tr_xxy_z[i] - 6.0 * tr_xxy_yyz[i] * tbe_0 - 10.0 * tr_xxy_yyz[i] * tke_0 + 4.0 * tr_xxy_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyy_yz[i] * tbe_0 + 8.0 * tr_xxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_yzz[i] = 2.0 * tr_xx_zz[i] - 4.0 * tr_xx_yyzz[i] * tke_0 - 6.0 * tr_xxy_yzz[i] * tbe_0 - 6.0 * tr_xxy_yzz[i] * tke_0 + 4.0 * tr_xxy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyy_zz[i] * tbe_0 + 8.0 * tr_xxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxy_zzz[i] = -4.0 * tr_xx_yzzz[i] * tke_0 - 6.0 * tr_xxy_zzz[i] * tbe_0 - 2.0 * tr_xxy_zzz[i] * tke_0 + 4.0 * tr_xxy_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 320-330 components of targeted buffer : FF

    auto tr_0_0_yy_xxz_xxx = pbuffer.data(idx_op_geom_020_ff + 320);

    auto tr_0_0_yy_xxz_xxy = pbuffer.data(idx_op_geom_020_ff + 321);

    auto tr_0_0_yy_xxz_xxz = pbuffer.data(idx_op_geom_020_ff + 322);

    auto tr_0_0_yy_xxz_xyy = pbuffer.data(idx_op_geom_020_ff + 323);

    auto tr_0_0_yy_xxz_xyz = pbuffer.data(idx_op_geom_020_ff + 324);

    auto tr_0_0_yy_xxz_xzz = pbuffer.data(idx_op_geom_020_ff + 325);

    auto tr_0_0_yy_xxz_yyy = pbuffer.data(idx_op_geom_020_ff + 326);

    auto tr_0_0_yy_xxz_yyz = pbuffer.data(idx_op_geom_020_ff + 327);

    auto tr_0_0_yy_xxz_yzz = pbuffer.data(idx_op_geom_020_ff + 328);

    auto tr_0_0_yy_xxz_zzz = pbuffer.data(idx_op_geom_020_ff + 329);

    #pragma omp simd aligned(tr_0_0_yy_xxz_xxx, tr_0_0_yy_xxz_xxy, tr_0_0_yy_xxz_xxz, tr_0_0_yy_xxz_xyy, tr_0_0_yy_xxz_xyz, tr_0_0_yy_xxz_xzz, tr_0_0_yy_xxz_yyy, tr_0_0_yy_xxz_yyz, tr_0_0_yy_xxz_yzz, tr_0_0_yy_xxz_zzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xxyz_xx, tr_xxyz_xxxy, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxxyy, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxz_xxx[i] = -2.0 * tr_xxz_xxx[i] * tbe_0 - 2.0 * tr_xxz_xxx[i] * tke_0 + 4.0 * tr_xxz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_xxy[i] = -2.0 * tr_xxz_xxy[i] * tbe_0 - 6.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xx[i] * tbe_0 + 8.0 * tr_xxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_xxz[i] = -2.0 * tr_xxz_xxz[i] * tbe_0 - 2.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_xyy[i] = 2.0 * tr_xxz_x[i] - 2.0 * tr_xxz_xyy[i] * tbe_0 - 10.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyz_xy[i] * tbe_0 + 8.0 * tr_xxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_xyz[i] = -2.0 * tr_xxz_xyz[i] * tbe_0 - 6.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 + 8.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_xzz[i] = -2.0 * tr_xxz_xzz[i] * tbe_0 - 2.0 * tr_xxz_xzz[i] * tke_0 + 4.0 * tr_xxz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_yyy[i] = 6.0 * tr_xxz_y[i] - 2.0 * tr_xxz_yyy[i] * tbe_0 - 14.0 * tr_xxz_yyy[i] * tke_0 + 4.0 * tr_xxz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxyz_yy[i] * tbe_0 + 8.0 * tr_xxyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_yyz[i] = 2.0 * tr_xxz_z[i] - 2.0 * tr_xxz_yyz[i] * tbe_0 - 10.0 * tr_xxz_yyz[i] * tke_0 + 4.0 * tr_xxz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyz_yz[i] * tbe_0 + 8.0 * tr_xxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_yzz[i] = -2.0 * tr_xxz_yzz[i] * tbe_0 - 6.0 * tr_xxz_yzz[i] * tke_0 + 4.0 * tr_xxz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_zz[i] * tbe_0 + 8.0 * tr_xxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxz_zzz[i] = -2.0 * tr_xxz_zzz[i] * tbe_0 - 2.0 * tr_xxz_zzz[i] * tke_0 + 4.0 * tr_xxz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-340 components of targeted buffer : FF

    auto tr_0_0_yy_xyy_xxx = pbuffer.data(idx_op_geom_020_ff + 330);

    auto tr_0_0_yy_xyy_xxy = pbuffer.data(idx_op_geom_020_ff + 331);

    auto tr_0_0_yy_xyy_xxz = pbuffer.data(idx_op_geom_020_ff + 332);

    auto tr_0_0_yy_xyy_xyy = pbuffer.data(idx_op_geom_020_ff + 333);

    auto tr_0_0_yy_xyy_xyz = pbuffer.data(idx_op_geom_020_ff + 334);

    auto tr_0_0_yy_xyy_xzz = pbuffer.data(idx_op_geom_020_ff + 335);

    auto tr_0_0_yy_xyy_yyy = pbuffer.data(idx_op_geom_020_ff + 336);

    auto tr_0_0_yy_xyy_yyz = pbuffer.data(idx_op_geom_020_ff + 337);

    auto tr_0_0_yy_xyy_yzz = pbuffer.data(idx_op_geom_020_ff + 338);

    auto tr_0_0_yy_xyy_zzz = pbuffer.data(idx_op_geom_020_ff + 339);

    #pragma omp simd aligned(tr_0_0_yy_xyy_xxx, tr_0_0_yy_xyy_xxy, tr_0_0_yy_xyy_xxz, tr_0_0_yy_xyy_xyy, tr_0_0_yy_xyy_xyz, tr_0_0_yy_xyy_xzz, tr_0_0_yy_xyy_yyy, tr_0_0_yy_xyy_yyz, tr_0_0_yy_xyy_yzz, tr_0_0_yy_xyy_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxxyy, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyy_xx, tr_xyyy_xxxy, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyyy_xxx, tr_xyyyy_xxy, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyy_xxx[i] = 2.0 * tr_x_xxx[i] - 8.0 * tr_xy_xxxy[i] * tke_0 - 10.0 * tr_xyy_xxx[i] * tbe_0 - 2.0 * tr_xyy_xxx[i] * tke_0 + 4.0 * tr_xyy_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_xxy[i] = 2.0 * tr_x_xxy[i] + 4.0 * tr_xy_xx[i] - 8.0 * tr_xy_xxyy[i] * tke_0 - 10.0 * tr_xyy_xxy[i] * tbe_0 - 6.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_xx[i] * tbe_0 + 8.0 * tr_xyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_xxz[i] = 2.0 * tr_x_xxz[i] - 8.0 * tr_xy_xxyz[i] * tke_0 - 10.0 * tr_xyy_xxz[i] * tbe_0 - 2.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_xyy[i] = 2.0 * tr_x_xyy[i] + 8.0 * tr_xy_xy[i] - 8.0 * tr_xy_xyyy[i] * tke_0 + 2.0 * tr_xyy_x[i] - 10.0 * tr_xyy_xyy[i] * tbe_0 - 10.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyy_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyy_xy[i] * tbe_0 + 8.0 * tr_xyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_xyz[i] = 2.0 * tr_x_xyz[i] + 4.0 * tr_xy_xz[i] - 8.0 * tr_xy_xyyz[i] * tke_0 - 10.0 * tr_xyy_xyz[i] * tbe_0 - 6.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_xz[i] * tbe_0 + 8.0 * tr_xyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_xzz[i] = 2.0 * tr_x_xzz[i] - 8.0 * tr_xy_xyzz[i] * tke_0 - 10.0 * tr_xyy_xzz[i] * tbe_0 - 2.0 * tr_xyy_xzz[i] * tke_0 + 4.0 * tr_xyy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_yyy[i] = 2.0 * tr_x_yyy[i] + 12.0 * tr_xy_yy[i] - 8.0 * tr_xy_yyyy[i] * tke_0 + 6.0 * tr_xyy_y[i] - 10.0 * tr_xyy_yyy[i] * tbe_0 - 14.0 * tr_xyy_yyy[i] * tke_0 + 4.0 * tr_xyy_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyyy_yy[i] * tbe_0 + 8.0 * tr_xyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_yyz[i] = 2.0 * tr_x_yyz[i] + 8.0 * tr_xy_yz[i] - 8.0 * tr_xy_yyyz[i] * tke_0 + 2.0 * tr_xyy_z[i] - 10.0 * tr_xyy_yyz[i] * tbe_0 - 10.0 * tr_xyy_yyz[i] * tke_0 + 4.0 * tr_xyy_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyy_yz[i] * tbe_0 + 8.0 * tr_xyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_yzz[i] = 2.0 * tr_x_yzz[i] + 4.0 * tr_xy_zz[i] - 8.0 * tr_xy_yyzz[i] * tke_0 - 10.0 * tr_xyy_yzz[i] * tbe_0 - 6.0 * tr_xyy_yzz[i] * tke_0 + 4.0 * tr_xyy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyy_zz[i] * tbe_0 + 8.0 * tr_xyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyy_zzz[i] = 2.0 * tr_x_zzz[i] - 8.0 * tr_xy_yzzz[i] * tke_0 - 10.0 * tr_xyy_zzz[i] * tbe_0 - 2.0 * tr_xyy_zzz[i] * tke_0 + 4.0 * tr_xyy_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 340-350 components of targeted buffer : FF

    auto tr_0_0_yy_xyz_xxx = pbuffer.data(idx_op_geom_020_ff + 340);

    auto tr_0_0_yy_xyz_xxy = pbuffer.data(idx_op_geom_020_ff + 341);

    auto tr_0_0_yy_xyz_xxz = pbuffer.data(idx_op_geom_020_ff + 342);

    auto tr_0_0_yy_xyz_xyy = pbuffer.data(idx_op_geom_020_ff + 343);

    auto tr_0_0_yy_xyz_xyz = pbuffer.data(idx_op_geom_020_ff + 344);

    auto tr_0_0_yy_xyz_xzz = pbuffer.data(idx_op_geom_020_ff + 345);

    auto tr_0_0_yy_xyz_yyy = pbuffer.data(idx_op_geom_020_ff + 346);

    auto tr_0_0_yy_xyz_yyz = pbuffer.data(idx_op_geom_020_ff + 347);

    auto tr_0_0_yy_xyz_yzz = pbuffer.data(idx_op_geom_020_ff + 348);

    auto tr_0_0_yy_xyz_zzz = pbuffer.data(idx_op_geom_020_ff + 349);

    #pragma omp simd aligned(tr_0_0_yy_xyz_xxx, tr_0_0_yy_xyz_xxy, tr_0_0_yy_xyz_xxz, tr_0_0_yy_xyz_xyy, tr_0_0_yy_xyz_xyz, tr_0_0_yy_xyz_xzz, tr_0_0_yy_xyz_yyy, tr_0_0_yy_xyz_yyz, tr_0_0_yy_xyz_yzz, tr_0_0_yy_xyz_zzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxy, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxxyy, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyz_xxx[i] = -4.0 * tr_xz_xxxy[i] * tke_0 - 6.0 * tr_xyz_xxx[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tke_0 + 4.0 * tr_xyz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_xxy[i] = 2.0 * tr_xz_xx[i] - 4.0 * tr_xz_xxyy[i] * tke_0 - 6.0 * tr_xyz_xxy[i] * tbe_0 - 6.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xx[i] * tbe_0 + 8.0 * tr_xyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_xxz[i] = -4.0 * tr_xz_xxyz[i] * tke_0 - 6.0 * tr_xyz_xxz[i] * tbe_0 - 2.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_xyy[i] = 4.0 * tr_xz_xy[i] - 4.0 * tr_xz_xyyy[i] * tke_0 + 2.0 * tr_xyz_x[i] - 6.0 * tr_xyz_xyy[i] * tbe_0 - 10.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyz_xy[i] * tbe_0 + 8.0 * tr_xyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_xyz[i] = 2.0 * tr_xz_xz[i] - 4.0 * tr_xz_xyyz[i] * tke_0 - 6.0 * tr_xyz_xyz[i] * tbe_0 - 6.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + 8.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_xzz[i] = -4.0 * tr_xz_xyzz[i] * tke_0 - 6.0 * tr_xyz_xzz[i] * tbe_0 - 2.0 * tr_xyz_xzz[i] * tke_0 + 4.0 * tr_xyz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_yyy[i] = 6.0 * tr_xz_yy[i] - 4.0 * tr_xz_yyyy[i] * tke_0 + 6.0 * tr_xyz_y[i] - 6.0 * tr_xyz_yyy[i] * tbe_0 - 14.0 * tr_xyz_yyy[i] * tke_0 + 4.0 * tr_xyz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyyz_yy[i] * tbe_0 + 8.0 * tr_xyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_yyz[i] = 4.0 * tr_xz_yz[i] - 4.0 * tr_xz_yyyz[i] * tke_0 + 2.0 * tr_xyz_z[i] - 6.0 * tr_xyz_yyz[i] * tbe_0 - 10.0 * tr_xyz_yyz[i] * tke_0 + 4.0 * tr_xyz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyz_yz[i] * tbe_0 + 8.0 * tr_xyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_yzz[i] = 2.0 * tr_xz_zz[i] - 4.0 * tr_xz_yyzz[i] * tke_0 - 6.0 * tr_xyz_yzz[i] * tbe_0 - 6.0 * tr_xyz_yzz[i] * tke_0 + 4.0 * tr_xyz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_zz[i] * tbe_0 + 8.0 * tr_xyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyz_zzz[i] = -4.0 * tr_xz_yzzz[i] * tke_0 - 6.0 * tr_xyz_zzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tke_0 + 4.0 * tr_xyz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 350-360 components of targeted buffer : FF

    auto tr_0_0_yy_xzz_xxx = pbuffer.data(idx_op_geom_020_ff + 350);

    auto tr_0_0_yy_xzz_xxy = pbuffer.data(idx_op_geom_020_ff + 351);

    auto tr_0_0_yy_xzz_xxz = pbuffer.data(idx_op_geom_020_ff + 352);

    auto tr_0_0_yy_xzz_xyy = pbuffer.data(idx_op_geom_020_ff + 353);

    auto tr_0_0_yy_xzz_xyz = pbuffer.data(idx_op_geom_020_ff + 354);

    auto tr_0_0_yy_xzz_xzz = pbuffer.data(idx_op_geom_020_ff + 355);

    auto tr_0_0_yy_xzz_yyy = pbuffer.data(idx_op_geom_020_ff + 356);

    auto tr_0_0_yy_xzz_yyz = pbuffer.data(idx_op_geom_020_ff + 357);

    auto tr_0_0_yy_xzz_yzz = pbuffer.data(idx_op_geom_020_ff + 358);

    auto tr_0_0_yy_xzz_zzz = pbuffer.data(idx_op_geom_020_ff + 359);

    #pragma omp simd aligned(tr_0_0_yy_xzz_xxx, tr_0_0_yy_xzz_xxy, tr_0_0_yy_xzz_xxz, tr_0_0_yy_xzz_xyy, tr_0_0_yy_xzz_xyz, tr_0_0_yy_xzz_xzz, tr_0_0_yy_xzz_yyy, tr_0_0_yy_xzz_yyz, tr_0_0_yy_xzz_yzz, tr_0_0_yy_xzz_zzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxy, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxxyy, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xzz_xxx[i] = -2.0 * tr_xzz_xxx[i] * tbe_0 - 2.0 * tr_xzz_xxx[i] * tke_0 + 4.0 * tr_xzz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_xxy[i] = -2.0 * tr_xzz_xxy[i] * tbe_0 - 6.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xx[i] * tbe_0 + 8.0 * tr_xyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_xxz[i] = -2.0 * tr_xzz_xxz[i] * tbe_0 - 2.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_xyy[i] = 2.0 * tr_xzz_x[i] - 2.0 * tr_xzz_xyy[i] * tbe_0 - 10.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyzz_xy[i] * tbe_0 + 8.0 * tr_xyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_xyz[i] = -2.0 * tr_xzz_xyz[i] * tbe_0 - 6.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + 8.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_xzz[i] = -2.0 * tr_xzz_xzz[i] * tbe_0 - 2.0 * tr_xzz_xzz[i] * tke_0 + 4.0 * tr_xzz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_yyy[i] = 6.0 * tr_xzz_y[i] - 2.0 * tr_xzz_yyy[i] * tbe_0 - 14.0 * tr_xzz_yyy[i] * tke_0 + 4.0 * tr_xzz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyzz_yy[i] * tbe_0 + 8.0 * tr_xyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_yyz[i] = 2.0 * tr_xzz_z[i] - 2.0 * tr_xzz_yyz[i] * tbe_0 - 10.0 * tr_xzz_yyz[i] * tke_0 + 4.0 * tr_xzz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyzz_yz[i] * tbe_0 + 8.0 * tr_xyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_yzz[i] = -2.0 * tr_xzz_yzz[i] * tbe_0 - 6.0 * tr_xzz_yzz[i] * tke_0 + 4.0 * tr_xzz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_zz[i] * tbe_0 + 8.0 * tr_xyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzz_zzz[i] = -2.0 * tr_xzz_zzz[i] * tbe_0 - 2.0 * tr_xzz_zzz[i] * tke_0 + 4.0 * tr_xzz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-370 components of targeted buffer : FF

    auto tr_0_0_yy_yyy_xxx = pbuffer.data(idx_op_geom_020_ff + 360);

    auto tr_0_0_yy_yyy_xxy = pbuffer.data(idx_op_geom_020_ff + 361);

    auto tr_0_0_yy_yyy_xxz = pbuffer.data(idx_op_geom_020_ff + 362);

    auto tr_0_0_yy_yyy_xyy = pbuffer.data(idx_op_geom_020_ff + 363);

    auto tr_0_0_yy_yyy_xyz = pbuffer.data(idx_op_geom_020_ff + 364);

    auto tr_0_0_yy_yyy_xzz = pbuffer.data(idx_op_geom_020_ff + 365);

    auto tr_0_0_yy_yyy_yyy = pbuffer.data(idx_op_geom_020_ff + 366);

    auto tr_0_0_yy_yyy_yyz = pbuffer.data(idx_op_geom_020_ff + 367);

    auto tr_0_0_yy_yyy_yzz = pbuffer.data(idx_op_geom_020_ff + 368);

    auto tr_0_0_yy_yyy_zzz = pbuffer.data(idx_op_geom_020_ff + 369);

    #pragma omp simd aligned(tr_0_0_yy_yyy_xxx, tr_0_0_yy_yyy_xxy, tr_0_0_yy_yyy_xxz, tr_0_0_yy_yyy_xyy, tr_0_0_yy_yyy_xyz, tr_0_0_yy_yyy_xzz, tr_0_0_yy_yyy_yyy, tr_0_0_yy_yyy_yyz, tr_0_0_yy_yyy_yzz, tr_0_0_yy_yyy_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yy_xx, tr_yy_xxxy, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxxyy, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyyyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyy_xx, tr_yyyy_xxxy, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyyy_xxx, tr_yyyyy_xxy, tr_yyyyy_xxz, tr_yyyyy_xyy, tr_yyyyy_xyz, tr_yyyyy_xzz, tr_yyyyy_yyy, tr_yyyyy_yyz, tr_yyyyy_yzz, tr_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyy_xxx[i] = 6.0 * tr_y_xxx[i] - 12.0 * tr_yy_xxxy[i] * tke_0 - 14.0 * tr_yyy_xxx[i] * tbe_0 - 2.0 * tr_yyy_xxx[i] * tke_0 + 4.0 * tr_yyy_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_xxy[i] = 6.0 * tr_y_xxy[i] + 6.0 * tr_yy_xx[i] - 12.0 * tr_yy_xxyy[i] * tke_0 - 14.0 * tr_yyy_xxy[i] * tbe_0 - 6.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyy_xx[i] * tbe_0 + 8.0 * tr_yyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_xxz[i] = 6.0 * tr_y_xxz[i] - 12.0 * tr_yy_xxyz[i] * tke_0 - 14.0 * tr_yyy_xxz[i] * tbe_0 - 2.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_xyy[i] = 6.0 * tr_y_xyy[i] + 12.0 * tr_yy_xy[i] - 12.0 * tr_yy_xyyy[i] * tke_0 + 2.0 * tr_yyy_x[i] - 14.0 * tr_yyy_xyy[i] * tbe_0 - 10.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyy_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyy_xy[i] * tbe_0 + 8.0 * tr_yyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_xyz[i] = 6.0 * tr_y_xyz[i] + 6.0 * tr_yy_xz[i] - 12.0 * tr_yy_xyyz[i] * tke_0 - 14.0 * tr_yyy_xyz[i] * tbe_0 - 6.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyy_xz[i] * tbe_0 + 8.0 * tr_yyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_xzz[i] = 6.0 * tr_y_xzz[i] - 12.0 * tr_yy_xyzz[i] * tke_0 - 14.0 * tr_yyy_xzz[i] * tbe_0 - 2.0 * tr_yyy_xzz[i] * tke_0 + 4.0 * tr_yyy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_yyy[i] = 6.0 * tr_y_yyy[i] + 18.0 * tr_yy_yy[i] - 12.0 * tr_yy_yyyy[i] * tke_0 + 6.0 * tr_yyy_y[i] - 14.0 * tr_yyy_yyy[i] * tbe_0 - 14.0 * tr_yyy_yyy[i] * tke_0 + 4.0 * tr_yyy_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyyy_yy[i] * tbe_0 + 8.0 * tr_yyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_yyz[i] = 6.0 * tr_y_yyz[i] + 12.0 * tr_yy_yz[i] - 12.0 * tr_yy_yyyz[i] * tke_0 + 2.0 * tr_yyy_z[i] - 14.0 * tr_yyy_yyz[i] * tbe_0 - 10.0 * tr_yyy_yyz[i] * tke_0 + 4.0 * tr_yyy_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyy_yz[i] * tbe_0 + 8.0 * tr_yyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_yzz[i] = 6.0 * tr_y_yzz[i] + 6.0 * tr_yy_zz[i] - 12.0 * tr_yy_yyzz[i] * tke_0 - 14.0 * tr_yyy_yzz[i] * tbe_0 - 6.0 * tr_yyy_yzz[i] * tke_0 + 4.0 * tr_yyy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyy_zz[i] * tbe_0 + 8.0 * tr_yyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyy_zzz[i] = 6.0 * tr_y_zzz[i] - 12.0 * tr_yy_yzzz[i] * tke_0 - 14.0 * tr_yyy_zzz[i] * tbe_0 - 2.0 * tr_yyy_zzz[i] * tke_0 + 4.0 * tr_yyy_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 370-380 components of targeted buffer : FF

    auto tr_0_0_yy_yyz_xxx = pbuffer.data(idx_op_geom_020_ff + 370);

    auto tr_0_0_yy_yyz_xxy = pbuffer.data(idx_op_geom_020_ff + 371);

    auto tr_0_0_yy_yyz_xxz = pbuffer.data(idx_op_geom_020_ff + 372);

    auto tr_0_0_yy_yyz_xyy = pbuffer.data(idx_op_geom_020_ff + 373);

    auto tr_0_0_yy_yyz_xyz = pbuffer.data(idx_op_geom_020_ff + 374);

    auto tr_0_0_yy_yyz_xzz = pbuffer.data(idx_op_geom_020_ff + 375);

    auto tr_0_0_yy_yyz_yyy = pbuffer.data(idx_op_geom_020_ff + 376);

    auto tr_0_0_yy_yyz_yyz = pbuffer.data(idx_op_geom_020_ff + 377);

    auto tr_0_0_yy_yyz_yzz = pbuffer.data(idx_op_geom_020_ff + 378);

    auto tr_0_0_yy_yyz_zzz = pbuffer.data(idx_op_geom_020_ff + 379);

    #pragma omp simd aligned(tr_0_0_yy_yyz_xxx, tr_0_0_yy_yyz_xxy, tr_0_0_yy_yyz_xxz, tr_0_0_yy_yyz_xyy, tr_0_0_yy_yyz_xyz, tr_0_0_yy_yyz_xzz, tr_0_0_yy_yyz_yyy, tr_0_0_yy_yyz_yyz, tr_0_0_yy_yyz_yzz, tr_0_0_yy_yyz_zzz, tr_yyyyz_xxx, tr_yyyyz_xxy, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxy, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxxyy, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyz_xxx[i] = 2.0 * tr_z_xxx[i] - 8.0 * tr_yz_xxxy[i] * tke_0 - 10.0 * tr_yyz_xxx[i] * tbe_0 - 2.0 * tr_yyz_xxx[i] * tke_0 + 4.0 * tr_yyz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_xxy[i] = 2.0 * tr_z_xxy[i] + 4.0 * tr_yz_xx[i] - 8.0 * tr_yz_xxyy[i] * tke_0 - 10.0 * tr_yyz_xxy[i] * tbe_0 - 6.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xx[i] * tbe_0 + 8.0 * tr_yyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_xxz[i] = 2.0 * tr_z_xxz[i] - 8.0 * tr_yz_xxyz[i] * tke_0 - 10.0 * tr_yyz_xxz[i] * tbe_0 - 2.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_xyy[i] = 2.0 * tr_z_xyy[i] + 8.0 * tr_yz_xy[i] - 8.0 * tr_yz_xyyy[i] * tke_0 + 2.0 * tr_yyz_x[i] - 10.0 * tr_yyz_xyy[i] * tbe_0 - 10.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyz_xy[i] * tbe_0 + 8.0 * tr_yyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_xyz[i] = 2.0 * tr_z_xyz[i] + 4.0 * tr_yz_xz[i] - 8.0 * tr_yz_xyyz[i] * tke_0 - 10.0 * tr_yyz_xyz[i] * tbe_0 - 6.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xz[i] * tbe_0 + 8.0 * tr_yyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_xzz[i] = 2.0 * tr_z_xzz[i] - 8.0 * tr_yz_xyzz[i] * tke_0 - 10.0 * tr_yyz_xzz[i] * tbe_0 - 2.0 * tr_yyz_xzz[i] * tke_0 + 4.0 * tr_yyz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_yyy[i] = 2.0 * tr_z_yyy[i] + 12.0 * tr_yz_yy[i] - 8.0 * tr_yz_yyyy[i] * tke_0 + 6.0 * tr_yyz_y[i] - 10.0 * tr_yyz_yyy[i] * tbe_0 - 14.0 * tr_yyz_yyy[i] * tke_0 + 4.0 * tr_yyz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyyz_yy[i] * tbe_0 + 8.0 * tr_yyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_yyz[i] = 2.0 * tr_z_yyz[i] + 8.0 * tr_yz_yz[i] - 8.0 * tr_yz_yyyz[i] * tke_0 + 2.0 * tr_yyz_z[i] - 10.0 * tr_yyz_yyz[i] * tbe_0 - 10.0 * tr_yyz_yyz[i] * tke_0 + 4.0 * tr_yyz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyz_yz[i] * tbe_0 + 8.0 * tr_yyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_yzz[i] = 2.0 * tr_z_yzz[i] + 4.0 * tr_yz_zz[i] - 8.0 * tr_yz_yyzz[i] * tke_0 - 10.0 * tr_yyz_yzz[i] * tbe_0 - 6.0 * tr_yyz_yzz[i] * tke_0 + 4.0 * tr_yyz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_zz[i] * tbe_0 + 8.0 * tr_yyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyz_zzz[i] = 2.0 * tr_z_zzz[i] - 8.0 * tr_yz_yzzz[i] * tke_0 - 10.0 * tr_yyz_zzz[i] * tbe_0 - 2.0 * tr_yyz_zzz[i] * tke_0 + 4.0 * tr_yyz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 380-390 components of targeted buffer : FF

    auto tr_0_0_yy_yzz_xxx = pbuffer.data(idx_op_geom_020_ff + 380);

    auto tr_0_0_yy_yzz_xxy = pbuffer.data(idx_op_geom_020_ff + 381);

    auto tr_0_0_yy_yzz_xxz = pbuffer.data(idx_op_geom_020_ff + 382);

    auto tr_0_0_yy_yzz_xyy = pbuffer.data(idx_op_geom_020_ff + 383);

    auto tr_0_0_yy_yzz_xyz = pbuffer.data(idx_op_geom_020_ff + 384);

    auto tr_0_0_yy_yzz_xzz = pbuffer.data(idx_op_geom_020_ff + 385);

    auto tr_0_0_yy_yzz_yyy = pbuffer.data(idx_op_geom_020_ff + 386);

    auto tr_0_0_yy_yzz_yyz = pbuffer.data(idx_op_geom_020_ff + 387);

    auto tr_0_0_yy_yzz_yzz = pbuffer.data(idx_op_geom_020_ff + 388);

    auto tr_0_0_yy_yzz_zzz = pbuffer.data(idx_op_geom_020_ff + 389);

    #pragma omp simd aligned(tr_0_0_yy_yzz_xxx, tr_0_0_yy_yzz_xxy, tr_0_0_yy_yzz_xxz, tr_0_0_yy_yzz_xyy, tr_0_0_yy_yzz_xyz, tr_0_0_yy_yzz_xzz, tr_0_0_yy_yzz_yyy, tr_0_0_yy_yzz_yyz, tr_0_0_yy_yzz_yzz, tr_0_0_yy_yzz_zzz, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, tr_yyzz_xx, tr_yyzz_xxxy, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxxyy, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_zz_xx, tr_zz_xxxy, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yzz_xxx[i] = -4.0 * tr_zz_xxxy[i] * tke_0 - 6.0 * tr_yzz_xxx[i] * tbe_0 - 2.0 * tr_yzz_xxx[i] * tke_0 + 4.0 * tr_yzz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_xxy[i] = 2.0 * tr_zz_xx[i] - 4.0 * tr_zz_xxyy[i] * tke_0 - 6.0 * tr_yzz_xxy[i] * tbe_0 - 6.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xx[i] * tbe_0 + 8.0 * tr_yyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_xxz[i] = -4.0 * tr_zz_xxyz[i] * tke_0 - 6.0 * tr_yzz_xxz[i] * tbe_0 - 2.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_xyy[i] = 4.0 * tr_zz_xy[i] - 4.0 * tr_zz_xyyy[i] * tke_0 + 2.0 * tr_yzz_x[i] - 6.0 * tr_yzz_xyy[i] * tbe_0 - 10.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyzz_xy[i] * tbe_0 + 8.0 * tr_yyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_xyz[i] = 2.0 * tr_zz_xz[i] - 4.0 * tr_zz_xyyz[i] * tke_0 - 6.0 * tr_yzz_xyz[i] * tbe_0 - 6.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xz[i] * tbe_0 + 8.0 * tr_yyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_xzz[i] = -4.0 * tr_zz_xyzz[i] * tke_0 - 6.0 * tr_yzz_xzz[i] * tbe_0 - 2.0 * tr_yzz_xzz[i] * tke_0 + 4.0 * tr_yzz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_yyy[i] = 6.0 * tr_zz_yy[i] - 4.0 * tr_zz_yyyy[i] * tke_0 + 6.0 * tr_yzz_y[i] - 6.0 * tr_yzz_yyy[i] * tbe_0 - 14.0 * tr_yzz_yyy[i] * tke_0 + 4.0 * tr_yzz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyzz_yy[i] * tbe_0 + 8.0 * tr_yyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_yyz[i] = 4.0 * tr_zz_yz[i] - 4.0 * tr_zz_yyyz[i] * tke_0 + 2.0 * tr_yzz_z[i] - 6.0 * tr_yzz_yyz[i] * tbe_0 - 10.0 * tr_yzz_yyz[i] * tke_0 + 4.0 * tr_yzz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyzz_yz[i] * tbe_0 + 8.0 * tr_yyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_yzz[i] = 2.0 * tr_zz_zz[i] - 4.0 * tr_zz_yyzz[i] * tke_0 - 6.0 * tr_yzz_yzz[i] * tbe_0 - 6.0 * tr_yzz_yzz[i] * tke_0 + 4.0 * tr_yzz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_zz[i] * tbe_0 + 8.0 * tr_yyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzz_zzz[i] = -4.0 * tr_zz_yzzz[i] * tke_0 - 6.0 * tr_yzz_zzz[i] * tbe_0 - 2.0 * tr_yzz_zzz[i] * tke_0 + 4.0 * tr_yzz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-400 components of targeted buffer : FF

    auto tr_0_0_yy_zzz_xxx = pbuffer.data(idx_op_geom_020_ff + 390);

    auto tr_0_0_yy_zzz_xxy = pbuffer.data(idx_op_geom_020_ff + 391);

    auto tr_0_0_yy_zzz_xxz = pbuffer.data(idx_op_geom_020_ff + 392);

    auto tr_0_0_yy_zzz_xyy = pbuffer.data(idx_op_geom_020_ff + 393);

    auto tr_0_0_yy_zzz_xyz = pbuffer.data(idx_op_geom_020_ff + 394);

    auto tr_0_0_yy_zzz_xzz = pbuffer.data(idx_op_geom_020_ff + 395);

    auto tr_0_0_yy_zzz_yyy = pbuffer.data(idx_op_geom_020_ff + 396);

    auto tr_0_0_yy_zzz_yyz = pbuffer.data(idx_op_geom_020_ff + 397);

    auto tr_0_0_yy_zzz_yzz = pbuffer.data(idx_op_geom_020_ff + 398);

    auto tr_0_0_yy_zzz_zzz = pbuffer.data(idx_op_geom_020_ff + 399);

    #pragma omp simd aligned(tr_0_0_yy_zzz_xxx, tr_0_0_yy_zzz_xxy, tr_0_0_yy_zzz_xxz, tr_0_0_yy_zzz_xyy, tr_0_0_yy_zzz_xyz, tr_0_0_yy_zzz_xzz, tr_0_0_yy_zzz_yyy, tr_0_0_yy_zzz_yyz, tr_0_0_yy_zzz_yzz, tr_0_0_yy_zzz_zzz, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, tr_yzzz_xx, tr_yzzz_xxxy, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxxyy, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyyyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_zzz_xxx[i] = -2.0 * tr_zzz_xxx[i] * tbe_0 - 2.0 * tr_zzz_xxx[i] * tke_0 + 4.0 * tr_zzz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_xxy[i] = -2.0 * tr_zzz_xxy[i] * tbe_0 - 6.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xx[i] * tbe_0 + 8.0 * tr_yzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_xxz[i] = -2.0 * tr_zzz_xxz[i] * tbe_0 - 2.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_xyy[i] = 2.0 * tr_zzz_x[i] - 2.0 * tr_zzz_xyy[i] * tbe_0 - 10.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yzzz_xy[i] * tbe_0 + 8.0 * tr_yzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_xyz[i] = -2.0 * tr_zzz_xyz[i] * tbe_0 - 6.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xz[i] * tbe_0 + 8.0 * tr_yzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_xzz[i] = -2.0 * tr_zzz_xzz[i] * tbe_0 - 2.0 * tr_zzz_xzz[i] * tke_0 + 4.0 * tr_zzz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_yyy[i] = 6.0 * tr_zzz_y[i] - 2.0 * tr_zzz_yyy[i] * tbe_0 - 14.0 * tr_zzz_yyy[i] * tke_0 + 4.0 * tr_zzz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yzzz_yy[i] * tbe_0 + 8.0 * tr_yzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_yyz[i] = 2.0 * tr_zzz_z[i] - 2.0 * tr_zzz_yyz[i] * tbe_0 - 10.0 * tr_zzz_yyz[i] * tke_0 + 4.0 * tr_zzz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yzzz_yz[i] * tbe_0 + 8.0 * tr_yzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_yzz[i] = -2.0 * tr_zzz_yzz[i] * tbe_0 - 6.0 * tr_zzz_yzz[i] * tke_0 + 4.0 * tr_zzz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_zz[i] * tbe_0 + 8.0 * tr_yzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzz_zzz[i] = -2.0 * tr_zzz_zzz[i] * tbe_0 - 2.0 * tr_zzz_zzz[i] * tke_0 + 4.0 * tr_zzz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 400-410 components of targeted buffer : FF

    auto tr_0_0_yz_xxx_xxx = pbuffer.data(idx_op_geom_020_ff + 400);

    auto tr_0_0_yz_xxx_xxy = pbuffer.data(idx_op_geom_020_ff + 401);

    auto tr_0_0_yz_xxx_xxz = pbuffer.data(idx_op_geom_020_ff + 402);

    auto tr_0_0_yz_xxx_xyy = pbuffer.data(idx_op_geom_020_ff + 403);

    auto tr_0_0_yz_xxx_xyz = pbuffer.data(idx_op_geom_020_ff + 404);

    auto tr_0_0_yz_xxx_xzz = pbuffer.data(idx_op_geom_020_ff + 405);

    auto tr_0_0_yz_xxx_yyy = pbuffer.data(idx_op_geom_020_ff + 406);

    auto tr_0_0_yz_xxx_yyz = pbuffer.data(idx_op_geom_020_ff + 407);

    auto tr_0_0_yz_xxx_yzz = pbuffer.data(idx_op_geom_020_ff + 408);

    auto tr_0_0_yz_xxx_zzz = pbuffer.data(idx_op_geom_020_ff + 409);

    #pragma omp simd aligned(tr_0_0_yz_xxx_xxx, tr_0_0_yz_xxx_xxy, tr_0_0_yz_xxx_xxz, tr_0_0_yz_xxx_xyy, tr_0_0_yz_xxx_xyz, tr_0_0_yz_xxx_xzz, tr_0_0_yz_xxx_yyy, tr_0_0_yz_xxx_yyz, tr_0_0_yz_xxx_yzz, tr_0_0_yz_xxx_zzz, tr_xxx_x, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_z, tr_xxx_zzz, tr_xxxy_xx, tr_xxxy_xxxz, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxxz_xx, tr_xxxz_xxxy, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxx_xxx[i] = 4.0 * tr_xxx_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_xxy[i] = -2.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_xxz[i] = -2.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_xyy[i] = -4.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_xy[i] * tbe_0 + 4.0 * tr_xxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_xyz[i] = tr_xxx_x[i] - 2.0 * tr_xxx_xzz[i] * tke_0 - 2.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxx_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_xz[i] * tbe_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xy[i] * tbe_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_xzz[i] = -4.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xz[i] * tbe_0 + 4.0 * tr_xxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_yyy[i] = -6.0 * tr_xxx_yyz[i] * tke_0 + 4.0 * tr_xxx_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxz_yy[i] * tbe_0 + 4.0 * tr_xxxz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_yyz[i] = 2.0 * tr_xxx_y[i] - 4.0 * tr_xxx_yzz[i] * tke_0 - 2.0 * tr_xxx_yyy[i] * tke_0 + 4.0 * tr_xxx_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_yz[i] * tbe_0 + 4.0 * tr_xxxz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yy[i] * tbe_0 + 4.0 * tr_xxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_yzz[i] = 2.0 * tr_xxx_z[i] - 2.0 * tr_xxx_zzz[i] * tke_0 - 4.0 * tr_xxx_yyz[i] * tke_0 + 4.0 * tr_xxx_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxz_zz[i] * tbe_0 + 4.0 * tr_xxxz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yz[i] * tbe_0 + 4.0 * tr_xxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxx_zzz[i] = -6.0 * tr_xxx_yzz[i] * tke_0 + 4.0 * tr_xxx_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxy_zz[i] * tbe_0 + 4.0 * tr_xxxy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 410-420 components of targeted buffer : FF

    auto tr_0_0_yz_xxy_xxx = pbuffer.data(idx_op_geom_020_ff + 410);

    auto tr_0_0_yz_xxy_xxy = pbuffer.data(idx_op_geom_020_ff + 411);

    auto tr_0_0_yz_xxy_xxz = pbuffer.data(idx_op_geom_020_ff + 412);

    auto tr_0_0_yz_xxy_xyy = pbuffer.data(idx_op_geom_020_ff + 413);

    auto tr_0_0_yz_xxy_xyz = pbuffer.data(idx_op_geom_020_ff + 414);

    auto tr_0_0_yz_xxy_xzz = pbuffer.data(idx_op_geom_020_ff + 415);

    auto tr_0_0_yz_xxy_yyy = pbuffer.data(idx_op_geom_020_ff + 416);

    auto tr_0_0_yz_xxy_yyz = pbuffer.data(idx_op_geom_020_ff + 417);

    auto tr_0_0_yz_xxy_yzz = pbuffer.data(idx_op_geom_020_ff + 418);

    auto tr_0_0_yz_xxy_zzz = pbuffer.data(idx_op_geom_020_ff + 419);

    #pragma omp simd aligned(tr_0_0_yz_xxy_xxx, tr_0_0_yz_xxy_xxy, tr_0_0_yz_xxy_xxz, tr_0_0_yz_xxy_xyy, tr_0_0_yz_xxy_xyz, tr_0_0_yz_xxy_xzz, tr_0_0_yz_xxy_yyy, tr_0_0_yz_xxy_yyz, tr_0_0_yz_xxy_yzz, tr_0_0_yz_xxy_zzz, tr_xx_xx, tr_xx_xxxz, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xx_zzzz, tr_xxy_x, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_z, tr_xxy_zzz, tr_xxyy_xx, tr_xxyy_xxxz, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xxyz_xx, tr_xxyz_xxxy, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxy_xxx[i] = -2.0 * tr_xx_xxxz[i] * tke_0 - 2.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_xxy[i] = -2.0 * tr_xx_xxyz[i] * tke_0 - 2.0 * tr_xxz_xxy[i] * tbe_0 - 2.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_xxz[i] = tr_xx_xx[i] - 2.0 * tr_xx_xxzz[i] * tke_0 - 2.0 * tr_xxz_xxz[i] * tbe_0 - 2.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_xyy[i] = -2.0 * tr_xx_xyyz[i] * tke_0 - 2.0 * tr_xxz_xyy[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_xyz[i] = tr_xx_xy[i] - 2.0 * tr_xx_xyzz[i] * tke_0 - 2.0 * tr_xxz_xyz[i] * tbe_0 + tr_xxy_x[i] - 2.0 * tr_xxy_xzz[i] * tke_0 - 2.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxy_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xy[i] * tbe_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_xzz[i] = 2.0 * tr_xx_xz[i] - 2.0 * tr_xx_xzzz[i] * tke_0 - 2.0 * tr_xxz_xzz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyy_xz[i] * tbe_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_yyy[i] = -2.0 * tr_xx_yyyz[i] * tke_0 - 2.0 * tr_xxz_yyy[i] * tbe_0 - 6.0 * tr_xxy_yyz[i] * tke_0 + 4.0 * tr_xxy_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_yyz[i] = tr_xx_yy[i] - 2.0 * tr_xx_yyzz[i] * tke_0 - 2.0 * tr_xxz_yyz[i] * tbe_0 + 2.0 * tr_xxy_y[i] - 4.0 * tr_xxy_yzz[i] * tke_0 - 2.0 * tr_xxy_yyy[i] * tke_0 + 4.0 * tr_xxy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yy[i] * tbe_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_yzz[i] = 2.0 * tr_xx_yz[i] - 2.0 * tr_xx_yzzz[i] * tke_0 - 2.0 * tr_xxz_yzz[i] * tbe_0 + 2.0 * tr_xxy_z[i] - 2.0 * tr_xxy_zzz[i] * tke_0 - 4.0 * tr_xxy_yyz[i] * tke_0 + 4.0 * tr_xxy_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyy_yz[i] * tbe_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxy_zzz[i] = 3.0 * tr_xx_zz[i] - 2.0 * tr_xx_zzzz[i] * tke_0 - 2.0 * tr_xxz_zzz[i] * tbe_0 - 6.0 * tr_xxy_yzz[i] * tke_0 + 4.0 * tr_xxy_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_zz[i] * tbe_0 + 4.0 * tr_xxyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-430 components of targeted buffer : FF

    auto tr_0_0_yz_xxz_xxx = pbuffer.data(idx_op_geom_020_ff + 420);

    auto tr_0_0_yz_xxz_xxy = pbuffer.data(idx_op_geom_020_ff + 421);

    auto tr_0_0_yz_xxz_xxz = pbuffer.data(idx_op_geom_020_ff + 422);

    auto tr_0_0_yz_xxz_xyy = pbuffer.data(idx_op_geom_020_ff + 423);

    auto tr_0_0_yz_xxz_xyz = pbuffer.data(idx_op_geom_020_ff + 424);

    auto tr_0_0_yz_xxz_xzz = pbuffer.data(idx_op_geom_020_ff + 425);

    auto tr_0_0_yz_xxz_yyy = pbuffer.data(idx_op_geom_020_ff + 426);

    auto tr_0_0_yz_xxz_yyz = pbuffer.data(idx_op_geom_020_ff + 427);

    auto tr_0_0_yz_xxz_yzz = pbuffer.data(idx_op_geom_020_ff + 428);

    auto tr_0_0_yz_xxz_zzz = pbuffer.data(idx_op_geom_020_ff + 429);

    #pragma omp simd aligned(tr_0_0_yz_xxz_xxx, tr_0_0_yz_xxz_xxy, tr_0_0_yz_xxz_xxz, tr_0_0_yz_xxz_xyy, tr_0_0_yz_xxz_xyz, tr_0_0_yz_xxz_xzz, tr_0_0_yz_xxz_yyy, tr_0_0_yz_xxz_yyz, tr_0_0_yz_xxz_yzz, tr_0_0_yz_xxz_zzz, tr_xx_xx, tr_xx_xxxy, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xxxz, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xxz_x, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xxxy, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxz_xxx[i] = -2.0 * tr_xx_xxxy[i] * tke_0 + 4.0 * tr_xxz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_xxy[i] = tr_xx_xx[i] - 2.0 * tr_xx_xxyy[i] * tke_0 - 2.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_xxz[i] = -2.0 * tr_xx_xxyz[i] * tke_0 - 2.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxz[i] * tbe_0 - 2.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_xyy[i] = 2.0 * tr_xx_xy[i] - 2.0 * tr_xx_xyyy[i] * tke_0 - 4.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_xy[i] * tbe_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_xyz[i] = tr_xx_xz[i] - 2.0 * tr_xx_xyyz[i] * tke_0 + tr_xxz_x[i] - 2.0 * tr_xxz_xzz[i] * tke_0 - 2.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_xz[i] * tbe_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyz[i] * tbe_0 - 2.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_xzz[i] = -2.0 * tr_xx_xyzz[i] * tke_0 - 4.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xzz[i] * tbe_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_yyy[i] = 3.0 * tr_xx_yy[i] - 2.0 * tr_xx_yyyy[i] * tke_0 - 6.0 * tr_xxz_yyz[i] * tke_0 + 4.0 * tr_xxz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxzz_yy[i] * tbe_0 + 4.0 * tr_xxzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_yyz[i] = 2.0 * tr_xx_yz[i] - 2.0 * tr_xx_yyyz[i] * tke_0 + 2.0 * tr_xxz_y[i] - 4.0 * tr_xxz_yzz[i] * tke_0 - 2.0 * tr_xxz_yyy[i] * tke_0 + 4.0 * tr_xxz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_yz[i] * tbe_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyz[i] * tbe_0 - 2.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_yzz[i] = tr_xx_zz[i] - 2.0 * tr_xx_yyzz[i] * tke_0 + 2.0 * tr_xxz_z[i] - 2.0 * tr_xxz_zzz[i] * tke_0 - 4.0 * tr_xxz_yyz[i] * tke_0 + 4.0 * tr_xxz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yzz[i] * tbe_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxz_zzz[i] = -2.0 * tr_xx_yzzz[i] * tke_0 - 6.0 * tr_xxz_yzz[i] * tke_0 + 4.0 * tr_xxz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_zzz[i] * tbe_0 - 6.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 430-440 components of targeted buffer : FF

    auto tr_0_0_yz_xyy_xxx = pbuffer.data(idx_op_geom_020_ff + 430);

    auto tr_0_0_yz_xyy_xxy = pbuffer.data(idx_op_geom_020_ff + 431);

    auto tr_0_0_yz_xyy_xxz = pbuffer.data(idx_op_geom_020_ff + 432);

    auto tr_0_0_yz_xyy_xyy = pbuffer.data(idx_op_geom_020_ff + 433);

    auto tr_0_0_yz_xyy_xyz = pbuffer.data(idx_op_geom_020_ff + 434);

    auto tr_0_0_yz_xyy_xzz = pbuffer.data(idx_op_geom_020_ff + 435);

    auto tr_0_0_yz_xyy_yyy = pbuffer.data(idx_op_geom_020_ff + 436);

    auto tr_0_0_yz_xyy_yyz = pbuffer.data(idx_op_geom_020_ff + 437);

    auto tr_0_0_yz_xyy_yzz = pbuffer.data(idx_op_geom_020_ff + 438);

    auto tr_0_0_yz_xyy_zzz = pbuffer.data(idx_op_geom_020_ff + 439);

    #pragma omp simd aligned(tr_0_0_yz_xyy_xxx, tr_0_0_yz_xyy_xxy, tr_0_0_yz_xyy_xxz, tr_0_0_yz_xyy_xyy, tr_0_0_yz_xyy_xyz, tr_0_0_yz_xyy_xzz, tr_0_0_yz_xyy_yyy, tr_0_0_yz_xyy_yyz, tr_0_0_yz_xyy_yzz, tr_0_0_yz_xyy_zzz, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyy_x, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_z, tr_xyy_zzz, tr_xyyy_xx, tr_xyyy_xxxz, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxy, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyy_xxx[i] = -4.0 * tr_xy_xxxz[i] * tke_0 - 4.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_xxy[i] = -4.0 * tr_xy_xxyz[i] * tke_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 - 2.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_xxz[i] = 2.0 * tr_xy_xx[i] - 4.0 * tr_xy_xxzz[i] * tke_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 - 2.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_xyy[i] = -4.0 * tr_xy_xyyz[i] * tke_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_xyz[i] = 2.0 * tr_xy_xy[i] - 4.0 * tr_xy_xyzz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + tr_xyy_x[i] - 2.0 * tr_xyy_xzz[i] * tke_0 - 2.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyy_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xy[i] * tbe_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_xzz[i] = 4.0 * tr_xy_xz[i] - 4.0 * tr_xy_xzzz[i] * tke_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyy_xz[i] * tbe_0 + 4.0 * tr_xyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_yyy[i] = -4.0 * tr_xy_yyyz[i] * tke_0 - 4.0 * tr_xyz_yyy[i] * tbe_0 - 6.0 * tr_xyy_yyz[i] * tke_0 + 4.0 * tr_xyy_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_yyz[i] = 2.0 * tr_xy_yy[i] - 4.0 * tr_xy_yyzz[i] * tke_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 + 2.0 * tr_xyy_y[i] - 4.0 * tr_xyy_yzz[i] * tke_0 - 2.0 * tr_xyy_yyy[i] * tke_0 + 4.0 * tr_xyy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yy[i] * tbe_0 + 4.0 * tr_xyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_yzz[i] = 4.0 * tr_xy_yz[i] - 4.0 * tr_xy_yzzz[i] * tke_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 + 2.0 * tr_xyy_z[i] - 2.0 * tr_xyy_zzz[i] * tke_0 - 4.0 * tr_xyy_yyz[i] * tke_0 + 4.0 * tr_xyy_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyy_yz[i] * tbe_0 + 4.0 * tr_xyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyy_zzz[i] = 6.0 * tr_xy_zz[i] - 4.0 * tr_xy_zzzz[i] * tke_0 - 4.0 * tr_xyz_zzz[i] * tbe_0 - 6.0 * tr_xyy_yzz[i] * tke_0 + 4.0 * tr_xyy_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 440-450 components of targeted buffer : FF

    auto tr_0_0_yz_xyz_xxx = pbuffer.data(idx_op_geom_020_ff + 440);

    auto tr_0_0_yz_xyz_xxy = pbuffer.data(idx_op_geom_020_ff + 441);

    auto tr_0_0_yz_xyz_xxz = pbuffer.data(idx_op_geom_020_ff + 442);

    auto tr_0_0_yz_xyz_xyy = pbuffer.data(idx_op_geom_020_ff + 443);

    auto tr_0_0_yz_xyz_xyz = pbuffer.data(idx_op_geom_020_ff + 444);

    auto tr_0_0_yz_xyz_xzz = pbuffer.data(idx_op_geom_020_ff + 445);

    auto tr_0_0_yz_xyz_yyy = pbuffer.data(idx_op_geom_020_ff + 446);

    auto tr_0_0_yz_xyz_yyz = pbuffer.data(idx_op_geom_020_ff + 447);

    auto tr_0_0_yz_xyz_yzz = pbuffer.data(idx_op_geom_020_ff + 448);

    auto tr_0_0_yz_xyz_zzz = pbuffer.data(idx_op_geom_020_ff + 449);

    #pragma omp simd aligned(tr_0_0_yz_xyz_xxx, tr_0_0_yz_xyz_xxy, tr_0_0_yz_xyz_xxz, tr_0_0_yz_xyz_xyy, tr_0_0_yz_xyz_xyz, tr_0_0_yz_xyz_xzz, tr_0_0_yz_xyz_yyy, tr_0_0_yz_xyz_yyz, tr_0_0_yz_xyz_yzz, tr_0_0_yz_xyz_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xxxz, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_xyz_x, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxy, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyz_xxx[i] = tr_x_xxx[i] - 2.0 * tr_xz_xxxz[i] * tke_0 - 2.0 * tr_xzz_xxx[i] * tbe_0 - 2.0 * tr_xy_xxxy[i] * tke_0 + 4.0 * tr_xyz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_xxy[i] = tr_x_xxy[i] - 2.0 * tr_xz_xxyz[i] * tke_0 - 2.0 * tr_xzz_xxy[i] * tbe_0 + tr_xy_xx[i] - 2.0 * tr_xy_xxyy[i] * tke_0 - 2.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_xxz[i] = tr_x_xxz[i] + tr_xz_xx[i] - 2.0 * tr_xz_xxzz[i] * tke_0 - 2.0 * tr_xzz_xxz[i] * tbe_0 - 2.0 * tr_xy_xxyz[i] * tke_0 - 2.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxz[i] * tbe_0 - 2.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_xyy[i] = tr_x_xyy[i] - 2.0 * tr_xz_xyyz[i] * tke_0 - 2.0 * tr_xzz_xyy[i] * tbe_0 + 2.0 * tr_xy_xy[i] - 2.0 * tr_xy_xyyy[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_xyz[i] = tr_x_xyz[i] + tr_xz_xy[i] - 2.0 * tr_xz_xyzz[i] * tke_0 - 2.0 * tr_xzz_xyz[i] * tbe_0 + tr_xy_xz[i] - 2.0 * tr_xy_xyyz[i] * tke_0 + tr_xyz_x[i] - 2.0 * tr_xyz_xzz[i] * tke_0 - 2.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyz[i] * tbe_0 - 2.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_xzz[i] = tr_x_xzz[i] + 2.0 * tr_xz_xz[i] - 2.0 * tr_xz_xzzz[i] * tke_0 - 2.0 * tr_xzz_xzz[i] * tbe_0 - 2.0 * tr_xy_xyzz[i] * tke_0 - 4.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xzz[i] * tbe_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_yyy[i] = tr_x_yyy[i] - 2.0 * tr_xz_yyyz[i] * tke_0 - 2.0 * tr_xzz_yyy[i] * tbe_0 + 3.0 * tr_xy_yy[i] - 2.0 * tr_xy_yyyy[i] * tke_0 - 6.0 * tr_xyz_yyz[i] * tke_0 + 4.0 * tr_xyz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_yyz[i] = tr_x_yyz[i] + tr_xz_yy[i] - 2.0 * tr_xz_yyzz[i] * tke_0 - 2.0 * tr_xzz_yyz[i] * tbe_0 + 2.0 * tr_xy_yz[i] - 2.0 * tr_xy_yyyz[i] * tke_0 + 2.0 * tr_xyz_y[i] - 4.0 * tr_xyz_yzz[i] * tke_0 - 2.0 * tr_xyz_yyy[i] * tke_0 + 4.0 * tr_xyz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyz[i] * tbe_0 - 2.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_yzz[i] = tr_x_yzz[i] + 2.0 * tr_xz_yz[i] - 2.0 * tr_xz_yzzz[i] * tke_0 - 2.0 * tr_xzz_yzz[i] * tbe_0 + tr_xy_zz[i] - 2.0 * tr_xy_yyzz[i] * tke_0 + 2.0 * tr_xyz_z[i] - 2.0 * tr_xyz_zzz[i] * tke_0 - 4.0 * tr_xyz_yyz[i] * tke_0 + 4.0 * tr_xyz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yzz[i] * tbe_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyz_zzz[i] = tr_x_zzz[i] + 3.0 * tr_xz_zz[i] - 2.0 * tr_xz_zzzz[i] * tke_0 - 2.0 * tr_xzz_zzz[i] * tbe_0 - 2.0 * tr_xy_yzzz[i] * tke_0 - 6.0 * tr_xyz_yzz[i] * tke_0 + 4.0 * tr_xyz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_zzz[i] * tbe_0 - 6.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-460 components of targeted buffer : FF

    auto tr_0_0_yz_xzz_xxx = pbuffer.data(idx_op_geom_020_ff + 450);

    auto tr_0_0_yz_xzz_xxy = pbuffer.data(idx_op_geom_020_ff + 451);

    auto tr_0_0_yz_xzz_xxz = pbuffer.data(idx_op_geom_020_ff + 452);

    auto tr_0_0_yz_xzz_xyy = pbuffer.data(idx_op_geom_020_ff + 453);

    auto tr_0_0_yz_xzz_xyz = pbuffer.data(idx_op_geom_020_ff + 454);

    auto tr_0_0_yz_xzz_xzz = pbuffer.data(idx_op_geom_020_ff + 455);

    auto tr_0_0_yz_xzz_yyy = pbuffer.data(idx_op_geom_020_ff + 456);

    auto tr_0_0_yz_xzz_yyz = pbuffer.data(idx_op_geom_020_ff + 457);

    auto tr_0_0_yz_xzz_yzz = pbuffer.data(idx_op_geom_020_ff + 458);

    auto tr_0_0_yz_xzz_zzz = pbuffer.data(idx_op_geom_020_ff + 459);

    #pragma omp simd aligned(tr_0_0_yz_xzz_xxx, tr_0_0_yz_xzz_xxy, tr_0_0_yz_xzz_xxz, tr_0_0_yz_xzz_xyy, tr_0_0_yz_xzz_xyz, tr_0_0_yz_xzz_xzz, tr_0_0_yz_xzz_yyy, tr_0_0_yz_xzz_yyz, tr_0_0_yz_xzz_yzz, tr_0_0_yz_xzz_zzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxz, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xzz_x, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xxxy, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xzz_xxx[i] = -4.0 * tr_xz_xxxy[i] * tke_0 + 4.0 * tr_xzz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_xxy[i] = 2.0 * tr_xz_xx[i] - 4.0 * tr_xz_xxyy[i] * tke_0 - 2.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_xxz[i] = -4.0 * tr_xz_xxyz[i] * tke_0 - 2.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 - 2.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_xyy[i] = 4.0 * tr_xz_xy[i] - 4.0 * tr_xz_xyyy[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_xy[i] * tbe_0 + 4.0 * tr_xzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_xyz[i] = 2.0 * tr_xz_xz[i] - 4.0 * tr_xz_xyyz[i] * tke_0 + tr_xzz_x[i] - 2.0 * tr_xzz_xzz[i] * tke_0 - 2.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_xz[i] * tbe_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 - 2.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_xzz[i] = -4.0 * tr_xz_xyzz[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_yyy[i] = 6.0 * tr_xz_yy[i] - 4.0 * tr_xz_yyyy[i] * tke_0 - 6.0 * tr_xzz_yyz[i] * tke_0 + 4.0 * tr_xzz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xzzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_yyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_yyz[i] = 4.0 * tr_xz_yz[i] - 4.0 * tr_xz_yyyz[i] * tke_0 + 2.0 * tr_xzz_y[i] - 4.0 * tr_xzz_yzz[i] * tke_0 - 2.0 * tr_xzz_yyy[i] * tke_0 + 4.0 * tr_xzz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_yz[i] * tbe_0 + 4.0 * tr_xzzz_yyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 - 2.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_yzz[i] = 2.0 * tr_xz_zz[i] - 4.0 * tr_xz_yyzz[i] * tke_0 + 2.0 * tr_xzz_z[i] - 2.0 * tr_xzz_zzz[i] * tke_0 - 4.0 * tr_xzz_yyz[i] * tke_0 + 4.0 * tr_xzz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzz_zz[i] * tbe_0 + 4.0 * tr_xzzz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzz_zzz[i] = -4.0 * tr_xz_yzzz[i] * tke_0 - 6.0 * tr_xzz_yzz[i] * tke_0 + 4.0 * tr_xzz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzz_yzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zzz[i] * tbe_0 - 6.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 460-470 components of targeted buffer : FF

    auto tr_0_0_yz_yyy_xxx = pbuffer.data(idx_op_geom_020_ff + 460);

    auto tr_0_0_yz_yyy_xxy = pbuffer.data(idx_op_geom_020_ff + 461);

    auto tr_0_0_yz_yyy_xxz = pbuffer.data(idx_op_geom_020_ff + 462);

    auto tr_0_0_yz_yyy_xyy = pbuffer.data(idx_op_geom_020_ff + 463);

    auto tr_0_0_yz_yyy_xyz = pbuffer.data(idx_op_geom_020_ff + 464);

    auto tr_0_0_yz_yyy_xzz = pbuffer.data(idx_op_geom_020_ff + 465);

    auto tr_0_0_yz_yyy_yyy = pbuffer.data(idx_op_geom_020_ff + 466);

    auto tr_0_0_yz_yyy_yyz = pbuffer.data(idx_op_geom_020_ff + 467);

    auto tr_0_0_yz_yyy_yzz = pbuffer.data(idx_op_geom_020_ff + 468);

    auto tr_0_0_yz_yyy_zzz = pbuffer.data(idx_op_geom_020_ff + 469);

    #pragma omp simd aligned(tr_0_0_yz_yyy_xxx, tr_0_0_yz_yyy_xxy, tr_0_0_yz_yyy_xxz, tr_0_0_yz_yyy_xyy, tr_0_0_yz_yyy_xyz, tr_0_0_yz_yyy_xzz, tr_0_0_yz_yyy_yyy, tr_0_0_yz_yyy_yyz, tr_0_0_yz_yyy_yzz, tr_0_0_yz_yyy_zzz, tr_yy_xx, tr_yy_xxxz, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yy_zzzz, tr_yyy_x, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_z, tr_yyy_zzz, tr_yyyy_xx, tr_yyyy_xxxz, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxy, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxy, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyy_xxx[i] = -6.0 * tr_yy_xxxz[i] * tke_0 - 6.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyy_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_xxy[i] = -6.0 * tr_yy_xxyz[i] * tke_0 - 6.0 * tr_yyz_xxy[i] * tbe_0 - 2.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_xxz[i] = 3.0 * tr_yy_xx[i] - 6.0 * tr_yy_xxzz[i] * tke_0 - 6.0 * tr_yyz_xxz[i] * tbe_0 - 2.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xx[i] * tbe_0 + 4.0 * tr_yyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_xyy[i] = -6.0 * tr_yy_xyyz[i] * tke_0 - 6.0 * tr_yyz_xyy[i] * tbe_0 - 4.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xy[i] * tbe_0 + 4.0 * tr_yyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_xyz[i] = 3.0 * tr_yy_xy[i] - 6.0 * tr_yy_xyzz[i] * tke_0 - 6.0 * tr_yyz_xyz[i] * tbe_0 + tr_yyy_x[i] - 2.0 * tr_yyy_xzz[i] * tke_0 - 2.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyy_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_xz[i] * tbe_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xy[i] * tbe_0 + 4.0 * tr_yyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_xzz[i] = 6.0 * tr_yy_xz[i] - 6.0 * tr_yy_xzzz[i] * tke_0 - 6.0 * tr_yyz_xzz[i] * tbe_0 - 4.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyy_xz[i] * tbe_0 + 4.0 * tr_yyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_yyy[i] = -6.0 * tr_yy_yyyz[i] * tke_0 - 6.0 * tr_yyz_yyy[i] * tbe_0 - 6.0 * tr_yyy_yyz[i] * tke_0 + 4.0 * tr_yyy_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyz_yy[i] * tbe_0 + 4.0 * tr_yyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_yyz[i] = 3.0 * tr_yy_yy[i] - 6.0 * tr_yy_yyzz[i] * tke_0 - 6.0 * tr_yyz_yyz[i] * tbe_0 + 2.0 * tr_yyy_y[i] - 4.0 * tr_yyy_yzz[i] * tke_0 - 2.0 * tr_yyy_yyy[i] * tke_0 + 4.0 * tr_yyy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_yz[i] * tbe_0 + 4.0 * tr_yyyz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yy[i] * tbe_0 + 4.0 * tr_yyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_yzz[i] = 6.0 * tr_yy_yz[i] - 6.0 * tr_yy_yzzz[i] * tke_0 - 6.0 * tr_yyz_yzz[i] * tbe_0 + 2.0 * tr_yyy_z[i] - 2.0 * tr_yyy_zzz[i] * tke_0 - 4.0 * tr_yyy_yyz[i] * tke_0 + 4.0 * tr_yyy_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyz_zz[i] * tbe_0 + 4.0 * tr_yyyz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyy_yz[i] * tbe_0 + 4.0 * tr_yyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyy_zzz[i] = 9.0 * tr_yy_zz[i] - 6.0 * tr_yy_zzzz[i] * tke_0 - 6.0 * tr_yyz_zzz[i] * tbe_0 - 6.0 * tr_yyy_yzz[i] * tke_0 + 4.0 * tr_yyy_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyyy_zz[i] * tbe_0 + 4.0 * tr_yyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 470-480 components of targeted buffer : FF

    auto tr_0_0_yz_yyz_xxx = pbuffer.data(idx_op_geom_020_ff + 470);

    auto tr_0_0_yz_yyz_xxy = pbuffer.data(idx_op_geom_020_ff + 471);

    auto tr_0_0_yz_yyz_xxz = pbuffer.data(idx_op_geom_020_ff + 472);

    auto tr_0_0_yz_yyz_xyy = pbuffer.data(idx_op_geom_020_ff + 473);

    auto tr_0_0_yz_yyz_xyz = pbuffer.data(idx_op_geom_020_ff + 474);

    auto tr_0_0_yz_yyz_xzz = pbuffer.data(idx_op_geom_020_ff + 475);

    auto tr_0_0_yz_yyz_yyy = pbuffer.data(idx_op_geom_020_ff + 476);

    auto tr_0_0_yz_yyz_yyz = pbuffer.data(idx_op_geom_020_ff + 477);

    auto tr_0_0_yz_yyz_yzz = pbuffer.data(idx_op_geom_020_ff + 478);

    auto tr_0_0_yz_yyz_zzz = pbuffer.data(idx_op_geom_020_ff + 479);

    #pragma omp simd aligned(tr_0_0_yz_yyz_xxx, tr_0_0_yz_yyz_xxy, tr_0_0_yz_yyz_xxz, tr_0_0_yz_yyz_xyy, tr_0_0_yz_yyz_xyz, tr_0_0_yz_yyz_xzz, tr_0_0_yz_yyz_yyy, tr_0_0_yz_yyz_yyz, tr_0_0_yz_yyz_yzz, tr_0_0_yz_yyz_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yy_xx, tr_yy_xxxy, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xxxz, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, tr_yyz_x, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xxxy, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyz_xxx[i] = 2.0 * tr_y_xxx[i] - 4.0 * tr_yz_xxxz[i] * tke_0 - 4.0 * tr_yzz_xxx[i] * tbe_0 - 2.0 * tr_yy_xxxy[i] * tke_0 + 4.0 * tr_yyz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_xxy[i] = 2.0 * tr_y_xxy[i] - 4.0 * tr_yz_xxyz[i] * tke_0 - 4.0 * tr_yzz_xxy[i] * tbe_0 + tr_yy_xx[i] - 2.0 * tr_yy_xxyy[i] * tke_0 - 2.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_xxz[i] = 2.0 * tr_y_xxz[i] + 2.0 * tr_yz_xx[i] - 4.0 * tr_yz_xxzz[i] * tke_0 - 4.0 * tr_yzz_xxz[i] * tbe_0 - 2.0 * tr_yy_xxyz[i] * tke_0 - 2.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxz[i] * tbe_0 - 2.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_xyy[i] = 2.0 * tr_y_xyy[i] - 4.0 * tr_yz_xyyz[i] * tke_0 - 4.0 * tr_yzz_xyy[i] * tbe_0 + 2.0 * tr_yy_xy[i] - 2.0 * tr_yy_xyyy[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xy[i] * tbe_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_xyz[i] = 2.0 * tr_y_xyz[i] + 2.0 * tr_yz_xy[i] - 4.0 * tr_yz_xyzz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tbe_0 + tr_yy_xz[i] - 2.0 * tr_yy_xyyz[i] * tke_0 + tr_yyz_x[i] - 2.0 * tr_yyz_xzz[i] * tke_0 - 2.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_xz[i] * tbe_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xyz[i] * tbe_0 - 2.0 * tr_yyyz_xy[i] * tbe_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_xzz[i] = 2.0 * tr_y_xzz[i] + 4.0 * tr_yz_xz[i] - 4.0 * tr_yz_xzzz[i] * tke_0 - 4.0 * tr_yzz_xzz[i] * tbe_0 - 2.0 * tr_yy_xyzz[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xzz[i] * tbe_0 - 4.0 * tr_yyyz_xz[i] * tbe_0 + 4.0 * tr_yyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_yyy[i] = 2.0 * tr_y_yyy[i] - 4.0 * tr_yz_yyyz[i] * tke_0 - 4.0 * tr_yzz_yyy[i] * tbe_0 + 3.0 * tr_yy_yy[i] - 2.0 * tr_yy_yyyy[i] * tke_0 - 6.0 * tr_yyz_yyz[i] * tke_0 + 4.0 * tr_yyz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yyzz_yy[i] * tbe_0 + 4.0 * tr_yyzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_yyz[i] = 2.0 * tr_y_yyz[i] + 2.0 * tr_yz_yy[i] - 4.0 * tr_yz_yyzz[i] * tke_0 - 4.0 * tr_yzz_yyz[i] * tbe_0 + 2.0 * tr_yy_yz[i] - 2.0 * tr_yy_yyyz[i] * tke_0 + 2.0 * tr_yyz_y[i] - 4.0 * tr_yyz_yzz[i] * tke_0 - 2.0 * tr_yyz_yyy[i] * tke_0 + 4.0 * tr_yyz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_yz[i] * tbe_0 + 4.0 * tr_yyzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yyz[i] * tbe_0 - 2.0 * tr_yyyz_yy[i] * tbe_0 + 4.0 * tr_yyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_yzz[i] = 2.0 * tr_y_yzz[i] + 4.0 * tr_yz_yz[i] - 4.0 * tr_yz_yzzz[i] * tke_0 - 4.0 * tr_yzz_yzz[i] * tbe_0 + tr_yy_zz[i] - 2.0 * tr_yy_yyzz[i] * tke_0 + 2.0 * tr_yyz_z[i] - 2.0 * tr_yyz_zzz[i] * tke_0 - 4.0 * tr_yyz_yyz[i] * tke_0 + 4.0 * tr_yyz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzz_zz[i] * tbe_0 + 4.0 * tr_yyzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yzz[i] * tbe_0 - 4.0 * tr_yyyz_yz[i] * tbe_0 + 4.0 * tr_yyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyz_zzz[i] = 2.0 * tr_y_zzz[i] + 6.0 * tr_yz_zz[i] - 4.0 * tr_yz_zzzz[i] * tke_0 - 4.0 * tr_yzz_zzz[i] * tbe_0 - 2.0 * tr_yy_yzzz[i] * tke_0 - 6.0 * tr_yyz_yzz[i] * tke_0 + 4.0 * tr_yyz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_zzz[i] * tbe_0 - 6.0 * tr_yyyz_zz[i] * tbe_0 + 4.0 * tr_yyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-490 components of targeted buffer : FF

    auto tr_0_0_yz_yzz_xxx = pbuffer.data(idx_op_geom_020_ff + 480);

    auto tr_0_0_yz_yzz_xxy = pbuffer.data(idx_op_geom_020_ff + 481);

    auto tr_0_0_yz_yzz_xxz = pbuffer.data(idx_op_geom_020_ff + 482);

    auto tr_0_0_yz_yzz_xyy = pbuffer.data(idx_op_geom_020_ff + 483);

    auto tr_0_0_yz_yzz_xyz = pbuffer.data(idx_op_geom_020_ff + 484);

    auto tr_0_0_yz_yzz_xzz = pbuffer.data(idx_op_geom_020_ff + 485);

    auto tr_0_0_yz_yzz_yyy = pbuffer.data(idx_op_geom_020_ff + 486);

    auto tr_0_0_yz_yzz_yyz = pbuffer.data(idx_op_geom_020_ff + 487);

    auto tr_0_0_yz_yzz_yzz = pbuffer.data(idx_op_geom_020_ff + 488);

    auto tr_0_0_yz_yzz_zzz = pbuffer.data(idx_op_geom_020_ff + 489);

    #pragma omp simd aligned(tr_0_0_yz_yzz_xxx, tr_0_0_yz_yzz_xxy, tr_0_0_yz_yzz_xxz, tr_0_0_yz_yzz_xyy, tr_0_0_yz_yzz_xyz, tr_0_0_yz_yzz_xzz, tr_0_0_yz_yzz_yyy, tr_0_0_yz_yzz_yyz, tr_0_0_yz_yzz_yzz, tr_0_0_yz_yzz_zzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xxxz, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yzz_x, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xxxy, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, tr_zz_xx, tr_zz_xxxz, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yzz_xxx[i] = 2.0 * tr_z_xxx[i] - 2.0 * tr_zz_xxxz[i] * tke_0 - 2.0 * tr_zzz_xxx[i] * tbe_0 - 4.0 * tr_yz_xxxy[i] * tke_0 + 4.0 * tr_yzz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_xxy[i] = 2.0 * tr_z_xxy[i] - 2.0 * tr_zz_xxyz[i] * tke_0 - 2.0 * tr_zzz_xxy[i] * tbe_0 + 2.0 * tr_yz_xx[i] - 4.0 * tr_yz_xxyy[i] * tke_0 - 2.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_xxz[i] = 2.0 * tr_z_xxz[i] + tr_zz_xx[i] - 2.0 * tr_zz_xxzz[i] * tke_0 - 2.0 * tr_zzz_xxz[i] * tbe_0 - 4.0 * tr_yz_xxyz[i] * tke_0 - 2.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxz[i] * tbe_0 - 2.0 * tr_yyzz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_xyy[i] = 2.0 * tr_z_xyy[i] - 2.0 * tr_zz_xyyz[i] * tke_0 - 2.0 * tr_zzz_xyy[i] * tbe_0 + 4.0 * tr_yz_xy[i] - 4.0 * tr_yz_xyyy[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xy[i] * tbe_0 + 4.0 * tr_yzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_xyz[i] = 2.0 * tr_z_xyz[i] + tr_zz_xy[i] - 2.0 * tr_zz_xyzz[i] * tke_0 - 2.0 * tr_zzz_xyz[i] * tbe_0 + 2.0 * tr_yz_xz[i] - 4.0 * tr_yz_xyyz[i] * tke_0 + tr_yzz_x[i] - 2.0 * tr_yzz_xzz[i] * tke_0 - 2.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_xz[i] * tbe_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xyz[i] * tbe_0 - 2.0 * tr_yyzz_xy[i] * tbe_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_xzz[i] = 2.0 * tr_z_xzz[i] + 2.0 * tr_zz_xz[i] - 2.0 * tr_zz_xzzz[i] * tke_0 - 2.0 * tr_zzz_xzz[i] * tbe_0 - 4.0 * tr_yz_xyzz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xzz[i] * tbe_0 - 4.0 * tr_yyzz_xz[i] * tbe_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_yyy[i] = 2.0 * tr_z_yyy[i] - 2.0 * tr_zz_yyyz[i] * tke_0 - 2.0 * tr_zzz_yyy[i] * tbe_0 + 6.0 * tr_yz_yy[i] - 4.0 * tr_yz_yyyy[i] * tke_0 - 6.0 * tr_yzz_yyz[i] * tke_0 + 4.0 * tr_yzz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yzzz_yy[i] * tbe_0 + 4.0 * tr_yzzz_yyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_yyz[i] = 2.0 * tr_z_yyz[i] + tr_zz_yy[i] - 2.0 * tr_zz_yyzz[i] * tke_0 - 2.0 * tr_zzz_yyz[i] * tbe_0 + 4.0 * tr_yz_yz[i] - 4.0 * tr_yz_yyyz[i] * tke_0 + 2.0 * tr_yzz_y[i] - 4.0 * tr_yzz_yzz[i] * tke_0 - 2.0 * tr_yzz_yyy[i] * tke_0 + 4.0 * tr_yzz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_yz[i] * tbe_0 + 4.0 * tr_yzzz_yyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yyz[i] * tbe_0 - 2.0 * tr_yyzz_yy[i] * tbe_0 + 4.0 * tr_yyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_yzz[i] = 2.0 * tr_z_yzz[i] + 2.0 * tr_zz_yz[i] - 2.0 * tr_zz_yzzz[i] * tke_0 - 2.0 * tr_zzz_yzz[i] * tbe_0 + 2.0 * tr_yz_zz[i] - 4.0 * tr_yz_yyzz[i] * tke_0 + 2.0 * tr_yzz_z[i] - 2.0 * tr_yzz_zzz[i] * tke_0 - 4.0 * tr_yzz_yyz[i] * tke_0 + 4.0 * tr_yzz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzz_zz[i] * tbe_0 + 4.0 * tr_yzzz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yzz[i] * tbe_0 - 4.0 * tr_yyzz_yz[i] * tbe_0 + 4.0 * tr_yyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzz_zzz[i] = 2.0 * tr_z_zzz[i] + 3.0 * tr_zz_zz[i] - 2.0 * tr_zz_zzzz[i] * tke_0 - 2.0 * tr_zzz_zzz[i] * tbe_0 - 4.0 * tr_yz_yzzz[i] * tke_0 - 6.0 * tr_yzz_yzz[i] * tke_0 + 4.0 * tr_yzz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzz_yzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_zzz[i] * tbe_0 - 6.0 * tr_yyzz_zz[i] * tbe_0 + 4.0 * tr_yyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 490-500 components of targeted buffer : FF

    auto tr_0_0_yz_zzz_xxx = pbuffer.data(idx_op_geom_020_ff + 490);

    auto tr_0_0_yz_zzz_xxy = pbuffer.data(idx_op_geom_020_ff + 491);

    auto tr_0_0_yz_zzz_xxz = pbuffer.data(idx_op_geom_020_ff + 492);

    auto tr_0_0_yz_zzz_xyy = pbuffer.data(idx_op_geom_020_ff + 493);

    auto tr_0_0_yz_zzz_xyz = pbuffer.data(idx_op_geom_020_ff + 494);

    auto tr_0_0_yz_zzz_xzz = pbuffer.data(idx_op_geom_020_ff + 495);

    auto tr_0_0_yz_zzz_yyy = pbuffer.data(idx_op_geom_020_ff + 496);

    auto tr_0_0_yz_zzz_yyz = pbuffer.data(idx_op_geom_020_ff + 497);

    auto tr_0_0_yz_zzz_yzz = pbuffer.data(idx_op_geom_020_ff + 498);

    auto tr_0_0_yz_zzz_zzz = pbuffer.data(idx_op_geom_020_ff + 499);

    #pragma omp simd aligned(tr_0_0_yz_zzz_xxx, tr_0_0_yz_zzz_xxy, tr_0_0_yz_zzz_xxz, tr_0_0_yz_zzz_xyy, tr_0_0_yz_zzz_xyz, tr_0_0_yz_zzz_xzz, tr_0_0_yz_zzz_yyy, tr_0_0_yz_zzz_yyz, tr_0_0_yz_zzz_yzz, tr_0_0_yz_zzz_zzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xxxz, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxy, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_zzz, tr_zz_xx, tr_zz_xxxy, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zzz_x, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xxxy, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_zzz_xxx[i] = -6.0 * tr_zz_xxxy[i] * tke_0 + 4.0 * tr_zzz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_xxy[i] = 3.0 * tr_zz_xx[i] - 6.0 * tr_zz_xxyy[i] * tke_0 - 2.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_xx[i] * tbe_0 + 4.0 * tr_zzzz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_xxz[i] = -6.0 * tr_zz_xxyz[i] * tke_0 - 2.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxz[i] * tbe_0 - 2.0 * tr_yzzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_xyy[i] = 6.0 * tr_zz_xy[i] - 6.0 * tr_zz_xyyy[i] * tke_0 - 4.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_xy[i] * tbe_0 + 4.0 * tr_zzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_xyz[i] = 3.0 * tr_zz_xz[i] - 6.0 * tr_zz_xyyz[i] * tke_0 + tr_zzz_x[i] - 2.0 * tr_zzz_xzz[i] * tke_0 - 2.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_xz[i] * tbe_0 + 4.0 * tr_zzzz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xyz[i] * tbe_0 - 2.0 * tr_yzzz_xy[i] * tbe_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_xzz[i] = -6.0 * tr_zz_xyzz[i] * tke_0 - 4.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xzz[i] * tbe_0 - 4.0 * tr_yzzz_xz[i] * tbe_0 + 4.0 * tr_yzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_yyy[i] = 9.0 * tr_zz_yy[i] - 6.0 * tr_zz_yyyy[i] * tke_0 - 6.0 * tr_zzz_yyz[i] * tke_0 + 4.0 * tr_zzz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_zzzz_yy[i] * tbe_0 + 4.0 * tr_zzzz_yyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_yyz[i] = 6.0 * tr_zz_yz[i] - 6.0 * tr_zz_yyyz[i] * tke_0 + 2.0 * tr_zzz_y[i] - 4.0 * tr_zzz_yzz[i] * tke_0 - 2.0 * tr_zzz_yyy[i] * tke_0 + 4.0 * tr_zzz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_yz[i] * tbe_0 + 4.0 * tr_zzzz_yyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yyz[i] * tbe_0 - 2.0 * tr_yzzz_yy[i] * tbe_0 + 4.0 * tr_yzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_yzz[i] = 3.0 * tr_zz_zz[i] - 6.0 * tr_zz_yyzz[i] * tke_0 + 2.0 * tr_zzz_z[i] - 2.0 * tr_zzz_zzz[i] * tke_0 - 4.0 * tr_zzz_yyz[i] * tke_0 + 4.0 * tr_zzz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzz_zz[i] * tbe_0 + 4.0 * tr_zzzz_yyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yzz[i] * tbe_0 - 4.0 * tr_yzzz_yz[i] * tbe_0 + 4.0 * tr_yzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzz_zzz[i] = -6.0 * tr_zz_yzzz[i] * tke_0 - 6.0 * tr_zzz_yzz[i] * tke_0 + 4.0 * tr_zzz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_zzz[i] * tbe_0 - 6.0 * tr_yzzz_zz[i] * tbe_0 + 4.0 * tr_yzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 500-510 components of targeted buffer : FF

    auto tr_0_0_zz_xxx_xxx = pbuffer.data(idx_op_geom_020_ff + 500);

    auto tr_0_0_zz_xxx_xxy = pbuffer.data(idx_op_geom_020_ff + 501);

    auto tr_0_0_zz_xxx_xxz = pbuffer.data(idx_op_geom_020_ff + 502);

    auto tr_0_0_zz_xxx_xyy = pbuffer.data(idx_op_geom_020_ff + 503);

    auto tr_0_0_zz_xxx_xyz = pbuffer.data(idx_op_geom_020_ff + 504);

    auto tr_0_0_zz_xxx_xzz = pbuffer.data(idx_op_geom_020_ff + 505);

    auto tr_0_0_zz_xxx_yyy = pbuffer.data(idx_op_geom_020_ff + 506);

    auto tr_0_0_zz_xxx_yyz = pbuffer.data(idx_op_geom_020_ff + 507);

    auto tr_0_0_zz_xxx_yzz = pbuffer.data(idx_op_geom_020_ff + 508);

    auto tr_0_0_zz_xxx_zzz = pbuffer.data(idx_op_geom_020_ff + 509);

    #pragma omp simd aligned(tr_0_0_zz_xxx_xxx, tr_0_0_zz_xxx_xxy, tr_0_0_zz_xxx_xxz, tr_0_0_zz_xxx_xyy, tr_0_0_zz_xxx_xyz, tr_0_0_zz_xxx_xzz, tr_0_0_zz_xxx_yyy, tr_0_0_zz_xxx_yyz, tr_0_0_zz_xxx_yzz, tr_0_0_zz_xxx_zzz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_z, tr_xxx_zzz, tr_xxx_zzzzz, tr_xxxz_xx, tr_xxxz_xxxz, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxx_xxx[i] = -2.0 * tr_xxx_xxx[i] * tbe_0 - 2.0 * tr_xxx_xxx[i] * tke_0 + 4.0 * tr_xxx_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_xxy[i] = -2.0 * tr_xxx_xxy[i] * tbe_0 - 2.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxx_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_xxz[i] = -2.0 * tr_xxx_xxz[i] * tbe_0 - 6.0 * tr_xxx_xxz[i] * tke_0 + 4.0 * tr_xxx_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_xx[i] * tbe_0 + 8.0 * tr_xxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_xyy[i] = -2.0 * tr_xxx_xyy[i] * tbe_0 - 2.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxx_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_xyz[i] = -2.0 * tr_xxx_xyz[i] * tbe_0 - 6.0 * tr_xxx_xyz[i] * tke_0 + 4.0 * tr_xxx_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_xy[i] * tbe_0 + 8.0 * tr_xxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_xzz[i] = 2.0 * tr_xxx_x[i] - 2.0 * tr_xxx_xzz[i] * tbe_0 - 10.0 * tr_xxx_xzz[i] * tke_0 + 4.0 * tr_xxx_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxz_xz[i] * tbe_0 + 8.0 * tr_xxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_yyy[i] = -2.0 * tr_xxx_yyy[i] * tbe_0 - 2.0 * tr_xxx_yyy[i] * tke_0 + 4.0 * tr_xxx_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_yyz[i] = -2.0 * tr_xxx_yyz[i] * tbe_0 - 6.0 * tr_xxx_yyz[i] * tke_0 + 4.0 * tr_xxx_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxz_yy[i] * tbe_0 + 8.0 * tr_xxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_yzz[i] = 2.0 * tr_xxx_y[i] - 2.0 * tr_xxx_yzz[i] * tbe_0 - 10.0 * tr_xxx_yzz[i] * tke_0 + 4.0 * tr_xxx_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxz_yz[i] * tbe_0 + 8.0 * tr_xxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxx_zzz[i] = 6.0 * tr_xxx_z[i] - 2.0 * tr_xxx_zzz[i] * tbe_0 - 14.0 * tr_xxx_zzz[i] * tke_0 + 4.0 * tr_xxx_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxz_zz[i] * tbe_0 + 8.0 * tr_xxxz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-520 components of targeted buffer : FF

    auto tr_0_0_zz_xxy_xxx = pbuffer.data(idx_op_geom_020_ff + 510);

    auto tr_0_0_zz_xxy_xxy = pbuffer.data(idx_op_geom_020_ff + 511);

    auto tr_0_0_zz_xxy_xxz = pbuffer.data(idx_op_geom_020_ff + 512);

    auto tr_0_0_zz_xxy_xyy = pbuffer.data(idx_op_geom_020_ff + 513);

    auto tr_0_0_zz_xxy_xyz = pbuffer.data(idx_op_geom_020_ff + 514);

    auto tr_0_0_zz_xxy_xzz = pbuffer.data(idx_op_geom_020_ff + 515);

    auto tr_0_0_zz_xxy_yyy = pbuffer.data(idx_op_geom_020_ff + 516);

    auto tr_0_0_zz_xxy_yyz = pbuffer.data(idx_op_geom_020_ff + 517);

    auto tr_0_0_zz_xxy_yzz = pbuffer.data(idx_op_geom_020_ff + 518);

    auto tr_0_0_zz_xxy_zzz = pbuffer.data(idx_op_geom_020_ff + 519);

    #pragma omp simd aligned(tr_0_0_zz_xxy_xxx, tr_0_0_zz_xxy_xxy, tr_0_0_zz_xxy_xxz, tr_0_0_zz_xxy_xyy, tr_0_0_zz_xxy_xyz, tr_0_0_zz_xxy_xzz, tr_0_0_zz_xxy_yyy, tr_0_0_zz_xxy_yyz, tr_0_0_zz_xxy_yzz, tr_0_0_zz_xxy_zzz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_z, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyz_xx, tr_xxyz_xxxz, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxy_xxx[i] = -2.0 * tr_xxy_xxx[i] * tbe_0 - 2.0 * tr_xxy_xxx[i] * tke_0 + 4.0 * tr_xxy_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_xxy[i] = -2.0 * tr_xxy_xxy[i] * tbe_0 - 2.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_xxz[i] = -2.0 * tr_xxy_xxz[i] * tbe_0 - 6.0 * tr_xxy_xxz[i] * tke_0 + 4.0 * tr_xxy_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xx[i] * tbe_0 + 8.0 * tr_xxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_xyy[i] = -2.0 * tr_xxy_xyy[i] * tbe_0 - 2.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_xyz[i] = -2.0 * tr_xxy_xyz[i] * tbe_0 - 6.0 * tr_xxy_xyz[i] * tke_0 + 4.0 * tr_xxy_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 8.0 * tr_xxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_xzz[i] = 2.0 * tr_xxy_x[i] - 2.0 * tr_xxy_xzz[i] * tbe_0 - 10.0 * tr_xxy_xzz[i] * tke_0 + 4.0 * tr_xxy_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyz_xz[i] * tbe_0 + 8.0 * tr_xxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_yyy[i] = -2.0 * tr_xxy_yyy[i] * tbe_0 - 2.0 * tr_xxy_yyy[i] * tke_0 + 4.0 * tr_xxy_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_yyz[i] = -2.0 * tr_xxy_yyz[i] * tbe_0 - 6.0 * tr_xxy_yyz[i] * tke_0 + 4.0 * tr_xxy_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyz_yy[i] * tbe_0 + 8.0 * tr_xxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_yzz[i] = 2.0 * tr_xxy_y[i] - 2.0 * tr_xxy_yzz[i] * tbe_0 - 10.0 * tr_xxy_yzz[i] * tke_0 + 4.0 * tr_xxy_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyz_yz[i] * tbe_0 + 8.0 * tr_xxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxy_zzz[i] = 6.0 * tr_xxy_z[i] - 2.0 * tr_xxy_zzz[i] * tbe_0 - 14.0 * tr_xxy_zzz[i] * tke_0 + 4.0 * tr_xxy_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxyz_zz[i] * tbe_0 + 8.0 * tr_xxyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 520-530 components of targeted buffer : FF

    auto tr_0_0_zz_xxz_xxx = pbuffer.data(idx_op_geom_020_ff + 520);

    auto tr_0_0_zz_xxz_xxy = pbuffer.data(idx_op_geom_020_ff + 521);

    auto tr_0_0_zz_xxz_xxz = pbuffer.data(idx_op_geom_020_ff + 522);

    auto tr_0_0_zz_xxz_xyy = pbuffer.data(idx_op_geom_020_ff + 523);

    auto tr_0_0_zz_xxz_xyz = pbuffer.data(idx_op_geom_020_ff + 524);

    auto tr_0_0_zz_xxz_xzz = pbuffer.data(idx_op_geom_020_ff + 525);

    auto tr_0_0_zz_xxz_yyy = pbuffer.data(idx_op_geom_020_ff + 526);

    auto tr_0_0_zz_xxz_yyz = pbuffer.data(idx_op_geom_020_ff + 527);

    auto tr_0_0_zz_xxz_yzz = pbuffer.data(idx_op_geom_020_ff + 528);

    auto tr_0_0_zz_xxz_zzz = pbuffer.data(idx_op_geom_020_ff + 529);

    #pragma omp simd aligned(tr_0_0_zz_xxz_xxx, tr_0_0_zz_xxz_xxy, tr_0_0_zz_xxz_xxz, tr_0_0_zz_xxz_xyy, tr_0_0_zz_xxz_xyz, tr_0_0_zz_xxz_xzz, tr_0_0_zz_xxz_yyy, tr_0_0_zz_xxz_yyz, tr_0_0_zz_xxz_yzz, tr_0_0_zz_xxz_zzz, tr_xx_xx, tr_xx_xxxz, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xx_zzzz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_z, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xx, tr_xxzz_xxxz, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxz_xxx[i] = -4.0 * tr_xx_xxxz[i] * tke_0 - 6.0 * tr_xxz_xxx[i] * tbe_0 - 2.0 * tr_xxz_xxx[i] * tke_0 + 4.0 * tr_xxz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_xxy[i] = -4.0 * tr_xx_xxyz[i] * tke_0 - 6.0 * tr_xxz_xxy[i] * tbe_0 - 2.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_xxz[i] = 2.0 * tr_xx_xx[i] - 4.0 * tr_xx_xxzz[i] * tke_0 - 6.0 * tr_xxz_xxz[i] * tbe_0 - 6.0 * tr_xxz_xxz[i] * tke_0 + 4.0 * tr_xxz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_xx[i] * tbe_0 + 8.0 * tr_xxzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_xyy[i] = -4.0 * tr_xx_xyyz[i] * tke_0 - 6.0 * tr_xxz_xyy[i] * tbe_0 - 2.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_xyz[i] = 2.0 * tr_xx_xy[i] - 4.0 * tr_xx_xyzz[i] * tke_0 - 6.0 * tr_xxz_xyz[i] * tbe_0 - 6.0 * tr_xxz_xyz[i] * tke_0 + 4.0 * tr_xxz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_xy[i] * tbe_0 + 8.0 * tr_xxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_xzz[i] = 4.0 * tr_xx_xz[i] - 4.0 * tr_xx_xzzz[i] * tke_0 + 2.0 * tr_xxz_x[i] - 6.0 * tr_xxz_xzz[i] * tbe_0 - 10.0 * tr_xxz_xzz[i] * tke_0 + 4.0 * tr_xxz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzz_xz[i] * tbe_0 + 8.0 * tr_xxzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_yyy[i] = -4.0 * tr_xx_yyyz[i] * tke_0 - 6.0 * tr_xxz_yyy[i] * tbe_0 - 2.0 * tr_xxz_yyy[i] * tke_0 + 4.0 * tr_xxz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_yyz[i] = 2.0 * tr_xx_yy[i] - 4.0 * tr_xx_yyzz[i] * tke_0 - 6.0 * tr_xxz_yyz[i] * tbe_0 - 6.0 * tr_xxz_yyz[i] * tke_0 + 4.0 * tr_xxz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzz_yy[i] * tbe_0 + 8.0 * tr_xxzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_yzz[i] = 4.0 * tr_xx_yz[i] - 4.0 * tr_xx_yzzz[i] * tke_0 + 2.0 * tr_xxz_y[i] - 6.0 * tr_xxz_yzz[i] * tbe_0 - 10.0 * tr_xxz_yzz[i] * tke_0 + 4.0 * tr_xxz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzz_yz[i] * tbe_0 + 8.0 * tr_xxzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxz_zzz[i] = 6.0 * tr_xx_zz[i] - 4.0 * tr_xx_zzzz[i] * tke_0 + 6.0 * tr_xxz_z[i] - 6.0 * tr_xxz_zzz[i] * tbe_0 - 14.0 * tr_xxz_zzz[i] * tke_0 + 4.0 * tr_xxz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxzz_zz[i] * tbe_0 + 8.0 * tr_xxzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 530-540 components of targeted buffer : FF

    auto tr_0_0_zz_xyy_xxx = pbuffer.data(idx_op_geom_020_ff + 530);

    auto tr_0_0_zz_xyy_xxy = pbuffer.data(idx_op_geom_020_ff + 531);

    auto tr_0_0_zz_xyy_xxz = pbuffer.data(idx_op_geom_020_ff + 532);

    auto tr_0_0_zz_xyy_xyy = pbuffer.data(idx_op_geom_020_ff + 533);

    auto tr_0_0_zz_xyy_xyz = pbuffer.data(idx_op_geom_020_ff + 534);

    auto tr_0_0_zz_xyy_xzz = pbuffer.data(idx_op_geom_020_ff + 535);

    auto tr_0_0_zz_xyy_yyy = pbuffer.data(idx_op_geom_020_ff + 536);

    auto tr_0_0_zz_xyy_yyz = pbuffer.data(idx_op_geom_020_ff + 537);

    auto tr_0_0_zz_xyy_yzz = pbuffer.data(idx_op_geom_020_ff + 538);

    auto tr_0_0_zz_xyy_zzz = pbuffer.data(idx_op_geom_020_ff + 539);

    #pragma omp simd aligned(tr_0_0_zz_xyy_xxx, tr_0_0_zz_xyy_xxy, tr_0_0_zz_xyy_xxz, tr_0_0_zz_xyy_xyy, tr_0_0_zz_xyy_xyz, tr_0_0_zz_xyy_xzz, tr_0_0_zz_xyy_yyy, tr_0_0_zz_xyy_yyz, tr_0_0_zz_xyy_yzz, tr_0_0_zz_xyy_zzz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_z, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyz_xx, tr_xyyz_xxxz, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyy_xxx[i] = -2.0 * tr_xyy_xxx[i] * tbe_0 - 2.0 * tr_xyy_xxx[i] * tke_0 + 4.0 * tr_xyy_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_xxy[i] = -2.0 * tr_xyy_xxy[i] * tbe_0 - 2.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_xxz[i] = -2.0 * tr_xyy_xxz[i] * tbe_0 - 6.0 * tr_xyy_xxz[i] * tke_0 + 4.0 * tr_xyy_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xx[i] * tbe_0 + 8.0 * tr_xyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_xyy[i] = -2.0 * tr_xyy_xyy[i] * tbe_0 - 2.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_xyz[i] = -2.0 * tr_xyy_xyz[i] * tbe_0 - 6.0 * tr_xyy_xyz[i] * tke_0 + 4.0 * tr_xyy_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 8.0 * tr_xyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_xzz[i] = 2.0 * tr_xyy_x[i] - 2.0 * tr_xyy_xzz[i] * tbe_0 - 10.0 * tr_xyy_xzz[i] * tke_0 + 4.0 * tr_xyy_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyz_xz[i] * tbe_0 + 8.0 * tr_xyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_yyy[i] = -2.0 * tr_xyy_yyy[i] * tbe_0 - 2.0 * tr_xyy_yyy[i] * tke_0 + 4.0 * tr_xyy_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_yyz[i] = -2.0 * tr_xyy_yyz[i] * tbe_0 - 6.0 * tr_xyy_yyz[i] * tke_0 + 4.0 * tr_xyy_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyz_yy[i] * tbe_0 + 8.0 * tr_xyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_yzz[i] = 2.0 * tr_xyy_y[i] - 2.0 * tr_xyy_yzz[i] * tbe_0 - 10.0 * tr_xyy_yzz[i] * tke_0 + 4.0 * tr_xyy_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyz_yz[i] * tbe_0 + 8.0 * tr_xyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyy_zzz[i] = 6.0 * tr_xyy_z[i] - 2.0 * tr_xyy_zzz[i] * tbe_0 - 14.0 * tr_xyy_zzz[i] * tke_0 + 4.0 * tr_xyy_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyyz_zz[i] * tbe_0 + 8.0 * tr_xyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 540-550 components of targeted buffer : FF

    auto tr_0_0_zz_xyz_xxx = pbuffer.data(idx_op_geom_020_ff + 540);

    auto tr_0_0_zz_xyz_xxy = pbuffer.data(idx_op_geom_020_ff + 541);

    auto tr_0_0_zz_xyz_xxz = pbuffer.data(idx_op_geom_020_ff + 542);

    auto tr_0_0_zz_xyz_xyy = pbuffer.data(idx_op_geom_020_ff + 543);

    auto tr_0_0_zz_xyz_xyz = pbuffer.data(idx_op_geom_020_ff + 544);

    auto tr_0_0_zz_xyz_xzz = pbuffer.data(idx_op_geom_020_ff + 545);

    auto tr_0_0_zz_xyz_yyy = pbuffer.data(idx_op_geom_020_ff + 546);

    auto tr_0_0_zz_xyz_yyz = pbuffer.data(idx_op_geom_020_ff + 547);

    auto tr_0_0_zz_xyz_yzz = pbuffer.data(idx_op_geom_020_ff + 548);

    auto tr_0_0_zz_xyz_zzz = pbuffer.data(idx_op_geom_020_ff + 549);

    #pragma omp simd aligned(tr_0_0_zz_xyz_xxx, tr_0_0_zz_xyz_xxy, tr_0_0_zz_xyz_xxz, tr_0_0_zz_xyz_xyy, tr_0_0_zz_xyz_xyz, tr_0_0_zz_xyz_xzz, tr_0_0_zz_xyz_yyy, tr_0_0_zz_xyz_yyz, tr_0_0_zz_xyz_yzz, tr_0_0_zz_xyz_zzz, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_z, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xx, tr_xyzz_xxxz, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyz_xxx[i] = -4.0 * tr_xy_xxxz[i] * tke_0 - 6.0 * tr_xyz_xxx[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tke_0 + 4.0 * tr_xyz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_xxy[i] = -4.0 * tr_xy_xxyz[i] * tke_0 - 6.0 * tr_xyz_xxy[i] * tbe_0 - 2.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_xxz[i] = 2.0 * tr_xy_xx[i] - 4.0 * tr_xy_xxzz[i] * tke_0 - 6.0 * tr_xyz_xxz[i] * tbe_0 - 6.0 * tr_xyz_xxz[i] * tke_0 + 4.0 * tr_xyz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xx[i] * tbe_0 + 8.0 * tr_xyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_xyy[i] = -4.0 * tr_xy_xyyz[i] * tke_0 - 6.0 * tr_xyz_xyy[i] * tbe_0 - 2.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_xyz[i] = 2.0 * tr_xy_xy[i] - 4.0 * tr_xy_xyzz[i] * tke_0 - 6.0 * tr_xyz_xyz[i] * tbe_0 - 6.0 * tr_xyz_xyz[i] * tke_0 + 4.0 * tr_xyz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + 8.0 * tr_xyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_xzz[i] = 4.0 * tr_xy_xz[i] - 4.0 * tr_xy_xzzz[i] * tke_0 + 2.0 * tr_xyz_x[i] - 6.0 * tr_xyz_xzz[i] * tbe_0 - 10.0 * tr_xyz_xzz[i] * tke_0 + 4.0 * tr_xyz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzz_xz[i] * tbe_0 + 8.0 * tr_xyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_yyy[i] = -4.0 * tr_xy_yyyz[i] * tke_0 - 6.0 * tr_xyz_yyy[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tke_0 + 4.0 * tr_xyz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_yyz[i] = 2.0 * tr_xy_yy[i] - 4.0 * tr_xy_yyzz[i] * tke_0 - 6.0 * tr_xyz_yyz[i] * tbe_0 - 6.0 * tr_xyz_yyz[i] * tke_0 + 4.0 * tr_xyz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzz_yy[i] * tbe_0 + 8.0 * tr_xyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_yzz[i] = 4.0 * tr_xy_yz[i] - 4.0 * tr_xy_yzzz[i] * tke_0 + 2.0 * tr_xyz_y[i] - 6.0 * tr_xyz_yzz[i] * tbe_0 - 10.0 * tr_xyz_yzz[i] * tke_0 + 4.0 * tr_xyz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzz_yz[i] * tbe_0 + 8.0 * tr_xyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyz_zzz[i] = 6.0 * tr_xy_zz[i] - 4.0 * tr_xy_zzzz[i] * tke_0 + 6.0 * tr_xyz_z[i] - 6.0 * tr_xyz_zzz[i] * tbe_0 - 14.0 * tr_xyz_zzz[i] * tke_0 + 4.0 * tr_xyz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyzz_zz[i] * tbe_0 + 8.0 * tr_xyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 550-560 components of targeted buffer : FF

    auto tr_0_0_zz_xzz_xxx = pbuffer.data(idx_op_geom_020_ff + 550);

    auto tr_0_0_zz_xzz_xxy = pbuffer.data(idx_op_geom_020_ff + 551);

    auto tr_0_0_zz_xzz_xxz = pbuffer.data(idx_op_geom_020_ff + 552);

    auto tr_0_0_zz_xzz_xyy = pbuffer.data(idx_op_geom_020_ff + 553);

    auto tr_0_0_zz_xzz_xyz = pbuffer.data(idx_op_geom_020_ff + 554);

    auto tr_0_0_zz_xzz_xzz = pbuffer.data(idx_op_geom_020_ff + 555);

    auto tr_0_0_zz_xzz_yyy = pbuffer.data(idx_op_geom_020_ff + 556);

    auto tr_0_0_zz_xzz_yyz = pbuffer.data(idx_op_geom_020_ff + 557);

    auto tr_0_0_zz_xzz_yzz = pbuffer.data(idx_op_geom_020_ff + 558);

    auto tr_0_0_zz_xzz_zzz = pbuffer.data(idx_op_geom_020_ff + 559);

    #pragma omp simd aligned(tr_0_0_zz_xzz_xxx, tr_0_0_zz_xzz_xxy, tr_0_0_zz_xzz_xxz, tr_0_0_zz_xzz_xyy, tr_0_0_zz_xzz_xyz, tr_0_0_zz_xzz_xzz, tr_0_0_zz_xzz_yyy, tr_0_0_zz_xzz_yyz, tr_0_0_zz_xzz_yzz, tr_0_0_zz_xzz_zzz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_z, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xx, tr_xzzz_xxxz, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxy, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xzz_xxx[i] = 2.0 * tr_x_xxx[i] - 8.0 * tr_xz_xxxz[i] * tke_0 - 10.0 * tr_xzz_xxx[i] * tbe_0 - 2.0 * tr_xzz_xxx[i] * tke_0 + 4.0 * tr_xzz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_xxy[i] = 2.0 * tr_x_xxy[i] - 8.0 * tr_xz_xxyz[i] * tke_0 - 10.0 * tr_xzz_xxy[i] * tbe_0 - 2.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_xxz[i] = 2.0 * tr_x_xxz[i] + 4.0 * tr_xz_xx[i] - 8.0 * tr_xz_xxzz[i] * tke_0 - 10.0 * tr_xzz_xxz[i] * tbe_0 - 6.0 * tr_xzz_xxz[i] * tke_0 + 4.0 * tr_xzz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_xx[i] * tbe_0 + 8.0 * tr_xzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_xyy[i] = 2.0 * tr_x_xyy[i] - 8.0 * tr_xz_xyyz[i] * tke_0 - 10.0 * tr_xzz_xyy[i] * tbe_0 - 2.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_xyz[i] = 2.0 * tr_x_xyz[i] + 4.0 * tr_xz_xy[i] - 8.0 * tr_xz_xyzz[i] * tke_0 - 10.0 * tr_xzz_xyz[i] * tbe_0 - 6.0 * tr_xzz_xyz[i] * tke_0 + 4.0 * tr_xzz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_xy[i] * tbe_0 + 8.0 * tr_xzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_xzz[i] = 2.0 * tr_x_xzz[i] + 8.0 * tr_xz_xz[i] - 8.0 * tr_xz_xzzz[i] * tke_0 + 2.0 * tr_xzz_x[i] - 10.0 * tr_xzz_xzz[i] * tbe_0 - 10.0 * tr_xzz_xzz[i] * tke_0 + 4.0 * tr_xzz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzz_xz[i] * tbe_0 + 8.0 * tr_xzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_yyy[i] = 2.0 * tr_x_yyy[i] - 8.0 * tr_xz_yyyz[i] * tke_0 - 10.0 * tr_xzz_yyy[i] * tbe_0 - 2.0 * tr_xzz_yyy[i] * tke_0 + 4.0 * tr_xzz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_yyz[i] = 2.0 * tr_x_yyz[i] + 4.0 * tr_xz_yy[i] - 8.0 * tr_xz_yyzz[i] * tke_0 - 10.0 * tr_xzz_yyz[i] * tbe_0 - 6.0 * tr_xzz_yyz[i] * tke_0 + 4.0 * tr_xzz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzz_yy[i] * tbe_0 + 8.0 * tr_xzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_yzz[i] = 2.0 * tr_x_yzz[i] + 8.0 * tr_xz_yz[i] - 8.0 * tr_xz_yzzz[i] * tke_0 + 2.0 * tr_xzz_y[i] - 10.0 * tr_xzz_yzz[i] * tbe_0 - 10.0 * tr_xzz_yzz[i] * tke_0 + 4.0 * tr_xzz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzz_yz[i] * tbe_0 + 8.0 * tr_xzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzz_zzz[i] = 2.0 * tr_x_zzz[i] + 12.0 * tr_xz_zz[i] - 8.0 * tr_xz_zzzz[i] * tke_0 + 6.0 * tr_xzz_z[i] - 10.0 * tr_xzz_zzz[i] * tbe_0 - 14.0 * tr_xzz_zzz[i] * tke_0 + 4.0 * tr_xzz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xzzz_zz[i] * tbe_0 + 8.0 * tr_xzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 560-570 components of targeted buffer : FF

    auto tr_0_0_zz_yyy_xxx = pbuffer.data(idx_op_geom_020_ff + 560);

    auto tr_0_0_zz_yyy_xxy = pbuffer.data(idx_op_geom_020_ff + 561);

    auto tr_0_0_zz_yyy_xxz = pbuffer.data(idx_op_geom_020_ff + 562);

    auto tr_0_0_zz_yyy_xyy = pbuffer.data(idx_op_geom_020_ff + 563);

    auto tr_0_0_zz_yyy_xyz = pbuffer.data(idx_op_geom_020_ff + 564);

    auto tr_0_0_zz_yyy_xzz = pbuffer.data(idx_op_geom_020_ff + 565);

    auto tr_0_0_zz_yyy_yyy = pbuffer.data(idx_op_geom_020_ff + 566);

    auto tr_0_0_zz_yyy_yyz = pbuffer.data(idx_op_geom_020_ff + 567);

    auto tr_0_0_zz_yyy_yzz = pbuffer.data(idx_op_geom_020_ff + 568);

    auto tr_0_0_zz_yyy_zzz = pbuffer.data(idx_op_geom_020_ff + 569);

    #pragma omp simd aligned(tr_0_0_zz_yyy_xxx, tr_0_0_zz_yyy_xxy, tr_0_0_zz_yyy_xxz, tr_0_0_zz_yyy_xyy, tr_0_0_zz_yyy_xyz, tr_0_0_zz_yyy_xzz, tr_0_0_zz_yyy_yyy, tr_0_0_zz_yyy_yyz, tr_0_0_zz_yyy_yzz, tr_0_0_zz_yyy_zzz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_z, tr_yyy_zzz, tr_yyy_zzzzz, tr_yyyz_xx, tr_yyyz_xxxz, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyy_xxx[i] = -2.0 * tr_yyy_xxx[i] * tbe_0 - 2.0 * tr_yyy_xxx[i] * tke_0 + 4.0 * tr_yyy_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_xxy[i] = -2.0 * tr_yyy_xxy[i] * tbe_0 - 2.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_xxz[i] = -2.0 * tr_yyy_xxz[i] * tbe_0 - 6.0 * tr_yyy_xxz[i] * tke_0 + 4.0 * tr_yyy_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xx[i] * tbe_0 + 8.0 * tr_yyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_xyy[i] = -2.0 * tr_yyy_xyy[i] * tbe_0 - 2.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_xyz[i] = -2.0 * tr_yyy_xyz[i] * tbe_0 - 6.0 * tr_yyy_xyz[i] * tke_0 + 4.0 * tr_yyy_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_xy[i] * tbe_0 + 8.0 * tr_yyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_xzz[i] = 2.0 * tr_yyy_x[i] - 2.0 * tr_yyy_xzz[i] * tbe_0 - 10.0 * tr_yyy_xzz[i] * tke_0 + 4.0 * tr_yyy_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyz_xz[i] * tbe_0 + 8.0 * tr_yyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_yyy[i] = -2.0 * tr_yyy_yyy[i] * tbe_0 - 2.0 * tr_yyy_yyy[i] * tke_0 + 4.0 * tr_yyy_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_yyz[i] = -2.0 * tr_yyy_yyz[i] * tbe_0 - 6.0 * tr_yyy_yyz[i] * tke_0 + 4.0 * tr_yyy_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyz_yy[i] * tbe_0 + 8.0 * tr_yyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_yzz[i] = 2.0 * tr_yyy_y[i] - 2.0 * tr_yyy_yzz[i] * tbe_0 - 10.0 * tr_yyy_yzz[i] * tke_0 + 4.0 * tr_yyy_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyz_yz[i] * tbe_0 + 8.0 * tr_yyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyy_zzz[i] = 6.0 * tr_yyy_z[i] - 2.0 * tr_yyy_zzz[i] * tbe_0 - 14.0 * tr_yyy_zzz[i] * tke_0 + 4.0 * tr_yyy_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyyz_zz[i] * tbe_0 + 8.0 * tr_yyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 570-580 components of targeted buffer : FF

    auto tr_0_0_zz_yyz_xxx = pbuffer.data(idx_op_geom_020_ff + 570);

    auto tr_0_0_zz_yyz_xxy = pbuffer.data(idx_op_geom_020_ff + 571);

    auto tr_0_0_zz_yyz_xxz = pbuffer.data(idx_op_geom_020_ff + 572);

    auto tr_0_0_zz_yyz_xyy = pbuffer.data(idx_op_geom_020_ff + 573);

    auto tr_0_0_zz_yyz_xyz = pbuffer.data(idx_op_geom_020_ff + 574);

    auto tr_0_0_zz_yyz_xzz = pbuffer.data(idx_op_geom_020_ff + 575);

    auto tr_0_0_zz_yyz_yyy = pbuffer.data(idx_op_geom_020_ff + 576);

    auto tr_0_0_zz_yyz_yyz = pbuffer.data(idx_op_geom_020_ff + 577);

    auto tr_0_0_zz_yyz_yzz = pbuffer.data(idx_op_geom_020_ff + 578);

    auto tr_0_0_zz_yyz_zzz = pbuffer.data(idx_op_geom_020_ff + 579);

    #pragma omp simd aligned(tr_0_0_zz_yyz_xxx, tr_0_0_zz_yyz_xxy, tr_0_0_zz_yyz_xxz, tr_0_0_zz_yyz_xyy, tr_0_0_zz_yyz_xyz, tr_0_0_zz_yyz_xzz, tr_0_0_zz_yyz_yyy, tr_0_0_zz_yyz_yyz, tr_0_0_zz_yyz_yzz, tr_0_0_zz_yyz_zzz, tr_yy_xx, tr_yy_xxxz, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yy_zzzz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_z, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xx, tr_yyzz_xxxz, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyz_xxx[i] = -4.0 * tr_yy_xxxz[i] * tke_0 - 6.0 * tr_yyz_xxx[i] * tbe_0 - 2.0 * tr_yyz_xxx[i] * tke_0 + 4.0 * tr_yyz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_xxy[i] = -4.0 * tr_yy_xxyz[i] * tke_0 - 6.0 * tr_yyz_xxy[i] * tbe_0 - 2.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_xxz[i] = 2.0 * tr_yy_xx[i] - 4.0 * tr_yy_xxzz[i] * tke_0 - 6.0 * tr_yyz_xxz[i] * tbe_0 - 6.0 * tr_yyz_xxz[i] * tke_0 + 4.0 * tr_yyz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xx[i] * tbe_0 + 8.0 * tr_yyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_xyy[i] = -4.0 * tr_yy_xyyz[i] * tke_0 - 6.0 * tr_yyz_xyy[i] * tbe_0 - 2.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_xyz[i] = 2.0 * tr_yy_xy[i] - 4.0 * tr_yy_xyzz[i] * tke_0 - 6.0 * tr_yyz_xyz[i] * tbe_0 - 6.0 * tr_yyz_xyz[i] * tke_0 + 4.0 * tr_yyz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_xy[i] * tbe_0 + 8.0 * tr_yyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_xzz[i] = 4.0 * tr_yy_xz[i] - 4.0 * tr_yy_xzzz[i] * tke_0 + 2.0 * tr_yyz_x[i] - 6.0 * tr_yyz_xzz[i] * tbe_0 - 10.0 * tr_yyz_xzz[i] * tke_0 + 4.0 * tr_yyz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzz_xz[i] * tbe_0 + 8.0 * tr_yyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_yyy[i] = -4.0 * tr_yy_yyyz[i] * tke_0 - 6.0 * tr_yyz_yyy[i] * tbe_0 - 2.0 * tr_yyz_yyy[i] * tke_0 + 4.0 * tr_yyz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_yyz[i] = 2.0 * tr_yy_yy[i] - 4.0 * tr_yy_yyzz[i] * tke_0 - 6.0 * tr_yyz_yyz[i] * tbe_0 - 6.0 * tr_yyz_yyz[i] * tke_0 + 4.0 * tr_yyz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzz_yy[i] * tbe_0 + 8.0 * tr_yyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_yzz[i] = 4.0 * tr_yy_yz[i] - 4.0 * tr_yy_yzzz[i] * tke_0 + 2.0 * tr_yyz_y[i] - 6.0 * tr_yyz_yzz[i] * tbe_0 - 10.0 * tr_yyz_yzz[i] * tke_0 + 4.0 * tr_yyz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzz_yz[i] * tbe_0 + 8.0 * tr_yyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyz_zzz[i] = 6.0 * tr_yy_zz[i] - 4.0 * tr_yy_zzzz[i] * tke_0 + 6.0 * tr_yyz_z[i] - 6.0 * tr_yyz_zzz[i] * tbe_0 - 14.0 * tr_yyz_zzz[i] * tke_0 + 4.0 * tr_yyz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyzz_zz[i] * tbe_0 + 8.0 * tr_yyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 580-590 components of targeted buffer : FF

    auto tr_0_0_zz_yzz_xxx = pbuffer.data(idx_op_geom_020_ff + 580);

    auto tr_0_0_zz_yzz_xxy = pbuffer.data(idx_op_geom_020_ff + 581);

    auto tr_0_0_zz_yzz_xxz = pbuffer.data(idx_op_geom_020_ff + 582);

    auto tr_0_0_zz_yzz_xyy = pbuffer.data(idx_op_geom_020_ff + 583);

    auto tr_0_0_zz_yzz_xyz = pbuffer.data(idx_op_geom_020_ff + 584);

    auto tr_0_0_zz_yzz_xzz = pbuffer.data(idx_op_geom_020_ff + 585);

    auto tr_0_0_zz_yzz_yyy = pbuffer.data(idx_op_geom_020_ff + 586);

    auto tr_0_0_zz_yzz_yyz = pbuffer.data(idx_op_geom_020_ff + 587);

    auto tr_0_0_zz_yzz_yzz = pbuffer.data(idx_op_geom_020_ff + 588);

    auto tr_0_0_zz_yzz_zzz = pbuffer.data(idx_op_geom_020_ff + 589);

    #pragma omp simd aligned(tr_0_0_zz_yzz_xxx, tr_0_0_zz_yzz_xxy, tr_0_0_zz_yzz_xxz, tr_0_0_zz_yzz_xyy, tr_0_0_zz_yzz_xyz, tr_0_0_zz_yzz_xzz, tr_0_0_zz_yzz_yyy, tr_0_0_zz_yzz_yyz, tr_0_0_zz_yzz_yzz, tr_0_0_zz_yzz_zzz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_z, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xx, tr_yzzz_xxxz, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxy, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yzz_xxx[i] = 2.0 * tr_y_xxx[i] - 8.0 * tr_yz_xxxz[i] * tke_0 - 10.0 * tr_yzz_xxx[i] * tbe_0 - 2.0 * tr_yzz_xxx[i] * tke_0 + 4.0 * tr_yzz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_xxy[i] = 2.0 * tr_y_xxy[i] - 8.0 * tr_yz_xxyz[i] * tke_0 - 10.0 * tr_yzz_xxy[i] * tbe_0 - 2.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_xxz[i] = 2.0 * tr_y_xxz[i] + 4.0 * tr_yz_xx[i] - 8.0 * tr_yz_xxzz[i] * tke_0 - 10.0 * tr_yzz_xxz[i] * tbe_0 - 6.0 * tr_yzz_xxz[i] * tke_0 + 4.0 * tr_yzz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xx[i] * tbe_0 + 8.0 * tr_yzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_xyy[i] = 2.0 * tr_y_xyy[i] - 8.0 * tr_yz_xyyz[i] * tke_0 - 10.0 * tr_yzz_xyy[i] * tbe_0 - 2.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_xyz[i] = 2.0 * tr_y_xyz[i] + 4.0 * tr_yz_xy[i] - 8.0 * tr_yz_xyzz[i] * tke_0 - 10.0 * tr_yzz_xyz[i] * tbe_0 - 6.0 * tr_yzz_xyz[i] * tke_0 + 4.0 * tr_yzz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_xy[i] * tbe_0 + 8.0 * tr_yzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_xzz[i] = 2.0 * tr_y_xzz[i] + 8.0 * tr_yz_xz[i] - 8.0 * tr_yz_xzzz[i] * tke_0 + 2.0 * tr_yzz_x[i] - 10.0 * tr_yzz_xzz[i] * tbe_0 - 10.0 * tr_yzz_xzz[i] * tke_0 + 4.0 * tr_yzz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzz_xz[i] * tbe_0 + 8.0 * tr_yzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_yyy[i] = 2.0 * tr_y_yyy[i] - 8.0 * tr_yz_yyyz[i] * tke_0 - 10.0 * tr_yzz_yyy[i] * tbe_0 - 2.0 * tr_yzz_yyy[i] * tke_0 + 4.0 * tr_yzz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_yyz[i] = 2.0 * tr_y_yyz[i] + 4.0 * tr_yz_yy[i] - 8.0 * tr_yz_yyzz[i] * tke_0 - 10.0 * tr_yzz_yyz[i] * tbe_0 - 6.0 * tr_yzz_yyz[i] * tke_0 + 4.0 * tr_yzz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzz_yy[i] * tbe_0 + 8.0 * tr_yzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_yzz[i] = 2.0 * tr_y_yzz[i] + 8.0 * tr_yz_yz[i] - 8.0 * tr_yz_yzzz[i] * tke_0 + 2.0 * tr_yzz_y[i] - 10.0 * tr_yzz_yzz[i] * tbe_0 - 10.0 * tr_yzz_yzz[i] * tke_0 + 4.0 * tr_yzz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzz_yz[i] * tbe_0 + 8.0 * tr_yzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzz_zzz[i] = 2.0 * tr_y_zzz[i] + 12.0 * tr_yz_zz[i] - 8.0 * tr_yz_zzzz[i] * tke_0 + 6.0 * tr_yzz_z[i] - 10.0 * tr_yzz_zzz[i] * tbe_0 - 14.0 * tr_yzz_zzz[i] * tke_0 + 4.0 * tr_yzz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yzzz_zz[i] * tbe_0 + 8.0 * tr_yzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 590-600 components of targeted buffer : FF

    auto tr_0_0_zz_zzz_xxx = pbuffer.data(idx_op_geom_020_ff + 590);

    auto tr_0_0_zz_zzz_xxy = pbuffer.data(idx_op_geom_020_ff + 591);

    auto tr_0_0_zz_zzz_xxz = pbuffer.data(idx_op_geom_020_ff + 592);

    auto tr_0_0_zz_zzz_xyy = pbuffer.data(idx_op_geom_020_ff + 593);

    auto tr_0_0_zz_zzz_xyz = pbuffer.data(idx_op_geom_020_ff + 594);

    auto tr_0_0_zz_zzz_xzz = pbuffer.data(idx_op_geom_020_ff + 595);

    auto tr_0_0_zz_zzz_yyy = pbuffer.data(idx_op_geom_020_ff + 596);

    auto tr_0_0_zz_zzz_yyz = pbuffer.data(idx_op_geom_020_ff + 597);

    auto tr_0_0_zz_zzz_yzz = pbuffer.data(idx_op_geom_020_ff + 598);

    auto tr_0_0_zz_zzz_zzz = pbuffer.data(idx_op_geom_020_ff + 599);

    #pragma omp simd aligned(tr_0_0_zz_zzz_xxx, tr_0_0_zz_zzz_xxy, tr_0_0_zz_zzz_xxz, tr_0_0_zz_zzz_xyy, tr_0_0_zz_zzz_xyz, tr_0_0_zz_zzz_xzz, tr_0_0_zz_zzz_yyy, tr_0_0_zz_zzz_yyz, tr_0_0_zz_zzz_yzz, tr_0_0_zz_zzz_zzz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, tr_zz_xx, tr_zz_xxxz, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zz_zzzz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_z, tr_zzz_zzz, tr_zzz_zzzzz, tr_zzzz_xx, tr_zzzz_xxxz, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzz_zzzz, tr_zzzzz_xxx, tr_zzzzz_xxy, tr_zzzzz_xxz, tr_zzzzz_xyy, tr_zzzzz_xyz, tr_zzzzz_xzz, tr_zzzzz_yyy, tr_zzzzz_yyz, tr_zzzzz_yzz, tr_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_zzz_xxx[i] = 6.0 * tr_z_xxx[i] - 12.0 * tr_zz_xxxz[i] * tke_0 - 14.0 * tr_zzz_xxx[i] * tbe_0 - 2.0 * tr_zzz_xxx[i] * tke_0 + 4.0 * tr_zzz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_xxy[i] = 6.0 * tr_z_xxy[i] - 12.0 * tr_zz_xxyz[i] * tke_0 - 14.0 * tr_zzz_xxy[i] * tbe_0 - 2.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_xxz[i] = 6.0 * tr_z_xxz[i] + 6.0 * tr_zz_xx[i] - 12.0 * tr_zz_xxzz[i] * tke_0 - 14.0 * tr_zzz_xxz[i] * tbe_0 - 6.0 * tr_zzz_xxz[i] * tke_0 + 4.0 * tr_zzz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_xx[i] * tbe_0 + 8.0 * tr_zzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_xyy[i] = 6.0 * tr_z_xyy[i] - 12.0 * tr_zz_xyyz[i] * tke_0 - 14.0 * tr_zzz_xyy[i] * tbe_0 - 2.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_xyz[i] = 6.0 * tr_z_xyz[i] + 6.0 * tr_zz_xy[i] - 12.0 * tr_zz_xyzz[i] * tke_0 - 14.0 * tr_zzz_xyz[i] * tbe_0 - 6.0 * tr_zzz_xyz[i] * tke_0 + 4.0 * tr_zzz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_xy[i] * tbe_0 + 8.0 * tr_zzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_xzz[i] = 6.0 * tr_z_xzz[i] + 12.0 * tr_zz_xz[i] - 12.0 * tr_zz_xzzz[i] * tke_0 + 2.0 * tr_zzz_x[i] - 14.0 * tr_zzz_xzz[i] * tbe_0 - 10.0 * tr_zzz_xzz[i] * tke_0 + 4.0 * tr_zzz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzzz_xz[i] * tbe_0 + 8.0 * tr_zzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_yyy[i] = 6.0 * tr_z_yyy[i] - 12.0 * tr_zz_yyyz[i] * tke_0 - 14.0 * tr_zzz_yyy[i] * tbe_0 - 2.0 * tr_zzz_yyy[i] * tke_0 + 4.0 * tr_zzz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_yyz[i] = 6.0 * tr_z_yyz[i] + 6.0 * tr_zz_yy[i] - 12.0 * tr_zz_yyzz[i] * tke_0 - 14.0 * tr_zzz_yyz[i] * tbe_0 - 6.0 * tr_zzz_yyz[i] * tke_0 + 4.0 * tr_zzz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzz_yy[i] * tbe_0 + 8.0 * tr_zzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_yzz[i] = 6.0 * tr_z_yzz[i] + 12.0 * tr_zz_yz[i] - 12.0 * tr_zz_yzzz[i] * tke_0 + 2.0 * tr_zzz_y[i] - 14.0 * tr_zzz_yzz[i] * tbe_0 - 10.0 * tr_zzz_yzz[i] * tke_0 + 4.0 * tr_zzz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzzz_yz[i] * tbe_0 + 8.0 * tr_zzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzz_zzz[i] = 6.0 * tr_z_zzz[i] + 18.0 * tr_zz_zz[i] - 12.0 * tr_zz_zzzz[i] * tke_0 + 6.0 * tr_zzz_z[i] - 14.0 * tr_zzz_zzz[i] * tbe_0 - 14.0 * tr_zzz_zzz[i] * tke_0 + 4.0 * tr_zzz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zzzz_zz[i] * tbe_0 + 8.0 * tr_zzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_zzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

