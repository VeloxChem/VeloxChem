#include "GeometricalDerivatives110ForDG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_dg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_dg,
                         const int idx_op_sg,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dg,
                         const int idx_op_ff,
                         const int idx_op_fh,
                         const int idx_op_gg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : SG

    auto tr_0_xxxx = pbuffer.data(idx_op_sg);

    auto tr_0_xxxy = pbuffer.data(idx_op_sg + 1);

    auto tr_0_xxxz = pbuffer.data(idx_op_sg + 2);

    auto tr_0_xxyy = pbuffer.data(idx_op_sg + 3);

    auto tr_0_xxyz = pbuffer.data(idx_op_sg + 4);

    auto tr_0_xxzz = pbuffer.data(idx_op_sg + 5);

    auto tr_0_xyyy = pbuffer.data(idx_op_sg + 6);

    auto tr_0_xyyz = pbuffer.data(idx_op_sg + 7);

    auto tr_0_xyzz = pbuffer.data(idx_op_sg + 8);

    auto tr_0_xzzz = pbuffer.data(idx_op_sg + 9);

    auto tr_0_yyyy = pbuffer.data(idx_op_sg + 10);

    auto tr_0_yyyz = pbuffer.data(idx_op_sg + 11);

    auto tr_0_yyzz = pbuffer.data(idx_op_sg + 12);

    auto tr_0_yzzz = pbuffer.data(idx_op_sg + 13);

    auto tr_0_zzzz = pbuffer.data(idx_op_sg + 14);

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

    // Set up components of auxiliary buffer : PH

    auto tr_x_xxxxx = pbuffer.data(idx_op_ph);

    auto tr_x_xxxxy = pbuffer.data(idx_op_ph + 1);

    auto tr_x_xxxxz = pbuffer.data(idx_op_ph + 2);

    auto tr_x_xxxyy = pbuffer.data(idx_op_ph + 3);

    auto tr_x_xxxyz = pbuffer.data(idx_op_ph + 4);

    auto tr_x_xxxzz = pbuffer.data(idx_op_ph + 5);

    auto tr_x_xxyyy = pbuffer.data(idx_op_ph + 6);

    auto tr_x_xxyyz = pbuffer.data(idx_op_ph + 7);

    auto tr_x_xxyzz = pbuffer.data(idx_op_ph + 8);

    auto tr_x_xxzzz = pbuffer.data(idx_op_ph + 9);

    auto tr_x_xyyyy = pbuffer.data(idx_op_ph + 10);

    auto tr_x_xyyyz = pbuffer.data(idx_op_ph + 11);

    auto tr_x_xyyzz = pbuffer.data(idx_op_ph + 12);

    auto tr_x_xyzzz = pbuffer.data(idx_op_ph + 13);

    auto tr_x_xzzzz = pbuffer.data(idx_op_ph + 14);

    auto tr_x_yyyyy = pbuffer.data(idx_op_ph + 15);

    auto tr_x_yyyyz = pbuffer.data(idx_op_ph + 16);

    auto tr_x_yyyzz = pbuffer.data(idx_op_ph + 17);

    auto tr_x_yyzzz = pbuffer.data(idx_op_ph + 18);

    auto tr_x_yzzzz = pbuffer.data(idx_op_ph + 19);

    auto tr_x_zzzzz = pbuffer.data(idx_op_ph + 20);

    auto tr_y_xxxxx = pbuffer.data(idx_op_ph + 21);

    auto tr_y_xxxxy = pbuffer.data(idx_op_ph + 22);

    auto tr_y_xxxxz = pbuffer.data(idx_op_ph + 23);

    auto tr_y_xxxyy = pbuffer.data(idx_op_ph + 24);

    auto tr_y_xxxyz = pbuffer.data(idx_op_ph + 25);

    auto tr_y_xxxzz = pbuffer.data(idx_op_ph + 26);

    auto tr_y_xxyyy = pbuffer.data(idx_op_ph + 27);

    auto tr_y_xxyyz = pbuffer.data(idx_op_ph + 28);

    auto tr_y_xxyzz = pbuffer.data(idx_op_ph + 29);

    auto tr_y_xxzzz = pbuffer.data(idx_op_ph + 30);

    auto tr_y_xyyyy = pbuffer.data(idx_op_ph + 31);

    auto tr_y_xyyyz = pbuffer.data(idx_op_ph + 32);

    auto tr_y_xyyzz = pbuffer.data(idx_op_ph + 33);

    auto tr_y_xyzzz = pbuffer.data(idx_op_ph + 34);

    auto tr_y_xzzzz = pbuffer.data(idx_op_ph + 35);

    auto tr_y_yyyyy = pbuffer.data(idx_op_ph + 36);

    auto tr_y_yyyyz = pbuffer.data(idx_op_ph + 37);

    auto tr_y_yyyzz = pbuffer.data(idx_op_ph + 38);

    auto tr_y_yyzzz = pbuffer.data(idx_op_ph + 39);

    auto tr_y_yzzzz = pbuffer.data(idx_op_ph + 40);

    auto tr_y_zzzzz = pbuffer.data(idx_op_ph + 41);

    auto tr_z_xxxxx = pbuffer.data(idx_op_ph + 42);

    auto tr_z_xxxxy = pbuffer.data(idx_op_ph + 43);

    auto tr_z_xxxxz = pbuffer.data(idx_op_ph + 44);

    auto tr_z_xxxyy = pbuffer.data(idx_op_ph + 45);

    auto tr_z_xxxyz = pbuffer.data(idx_op_ph + 46);

    auto tr_z_xxxzz = pbuffer.data(idx_op_ph + 47);

    auto tr_z_xxyyy = pbuffer.data(idx_op_ph + 48);

    auto tr_z_xxyyz = pbuffer.data(idx_op_ph + 49);

    auto tr_z_xxyzz = pbuffer.data(idx_op_ph + 50);

    auto tr_z_xxzzz = pbuffer.data(idx_op_ph + 51);

    auto tr_z_xyyyy = pbuffer.data(idx_op_ph + 52);

    auto tr_z_xyyyz = pbuffer.data(idx_op_ph + 53);

    auto tr_z_xyyzz = pbuffer.data(idx_op_ph + 54);

    auto tr_z_xyzzz = pbuffer.data(idx_op_ph + 55);

    auto tr_z_xzzzz = pbuffer.data(idx_op_ph + 56);

    auto tr_z_yyyyy = pbuffer.data(idx_op_ph + 57);

    auto tr_z_yyyyz = pbuffer.data(idx_op_ph + 58);

    auto tr_z_yyyzz = pbuffer.data(idx_op_ph + 59);

    auto tr_z_yyzzz = pbuffer.data(idx_op_ph + 60);

    auto tr_z_yzzzz = pbuffer.data(idx_op_ph + 61);

    auto tr_z_zzzzz = pbuffer.data(idx_op_ph + 62);

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

    // Set up 0-15 components of targeted buffer : DG

    auto tr_x_0_x_xx_xxxx = pbuffer.data(idx_op_geom_110_dg);

    auto tr_x_0_x_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 1);

    auto tr_x_0_x_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 2);

    auto tr_x_0_x_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 3);

    auto tr_x_0_x_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 4);

    auto tr_x_0_x_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 5);

    auto tr_x_0_x_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 6);

    auto tr_x_0_x_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 7);

    auto tr_x_0_x_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 8);

    auto tr_x_0_x_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 9);

    auto tr_x_0_x_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 10);

    auto tr_x_0_x_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 11);

    auto tr_x_0_x_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 12);

    auto tr_x_0_x_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 13);

    auto tr_x_0_x_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 14);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_x_0_x_xx_xxxx, tr_x_0_x_xx_xxxy, tr_x_0_x_xx_xxxz, tr_x_0_x_xx_xxyy, tr_x_0_x_xx_xxyz, tr_x_0_x_xx_xxzz, tr_x_0_x_xx_xyyy, tr_x_0_x_xx_xyyz, tr_x_0_x_xx_xyzz, tr_x_0_x_xx_xzzz, tr_x_0_x_xx_yyyy, tr_x_0_x_xx_yyyz, tr_x_0_x_xx_yyzz, tr_x_0_x_xx_yzzz, tr_x_0_x_xx_zzzz, tr_x_xxx, tr_x_xxxxx, tr_x_xxxxy, tr_x_xxxxz, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxx_xxx, tr_xxx_xxxxx, tr_xxx_xxxxy, tr_xxx_xxxxz, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xzzz, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yzzz, tr_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xx_xxxx[i] = 2.0 * tr_0_xxxx[i] + 8.0 * tr_x_xxx[i] - 4.0 * tr_x_xxxxx[i] * tke_0 - 10.0 * tr_xx_xxxx[i] * tbe_0 - 8.0 * tr_xxx_xxx[i] * tbe_0 + 4.0 * tr_xxx_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xxxy[i] = 2.0 * tr_0_xxxy[i] + 6.0 * tr_x_xxy[i] - 4.0 * tr_x_xxxxy[i] * tke_0 - 10.0 * tr_xx_xxxy[i] * tbe_0 - 6.0 * tr_xxx_xxy[i] * tbe_0 + 4.0 * tr_xxx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xxxz[i] = 2.0 * tr_0_xxxz[i] + 6.0 * tr_x_xxz[i] - 4.0 * tr_x_xxxxz[i] * tke_0 - 10.0 * tr_xx_xxxz[i] * tbe_0 - 6.0 * tr_xxx_xxz[i] * tbe_0 + 4.0 * tr_xxx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xxyy[i] = 2.0 * tr_0_xxyy[i] + 4.0 * tr_x_xyy[i] - 4.0 * tr_x_xxxyy[i] * tke_0 - 10.0 * tr_xx_xxyy[i] * tbe_0 - 4.0 * tr_xxx_xyy[i] * tbe_0 + 4.0 * tr_xxx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xxyz[i] = 2.0 * tr_0_xxyz[i] + 4.0 * tr_x_xyz[i] - 4.0 * tr_x_xxxyz[i] * tke_0 - 10.0 * tr_xx_xxyz[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tbe_0 + 4.0 * tr_xxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xxzz[i] = 2.0 * tr_0_xxzz[i] + 4.0 * tr_x_xzz[i] - 4.0 * tr_x_xxxzz[i] * tke_0 - 10.0 * tr_xx_xxzz[i] * tbe_0 - 4.0 * tr_xxx_xzz[i] * tbe_0 + 4.0 * tr_xxx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xyyy[i] = 2.0 * tr_0_xyyy[i] + 2.0 * tr_x_yyy[i] - 4.0 * tr_x_xxyyy[i] * tke_0 - 10.0 * tr_xx_xyyy[i] * tbe_0 - 2.0 * tr_xxx_yyy[i] * tbe_0 + 4.0 * tr_xxx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xyyz[i] = 2.0 * tr_0_xyyz[i] + 2.0 * tr_x_yyz[i] - 4.0 * tr_x_xxyyz[i] * tke_0 - 10.0 * tr_xx_xyyz[i] * tbe_0 - 2.0 * tr_xxx_yyz[i] * tbe_0 + 4.0 * tr_xxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xyzz[i] = 2.0 * tr_0_xyzz[i] + 2.0 * tr_x_yzz[i] - 4.0 * tr_x_xxyzz[i] * tke_0 - 10.0 * tr_xx_xyzz[i] * tbe_0 - 2.0 * tr_xxx_yzz[i] * tbe_0 + 4.0 * tr_xxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_xzzz[i] = 2.0 * tr_0_xzzz[i] + 2.0 * tr_x_zzz[i] - 4.0 * tr_x_xxzzz[i] * tke_0 - 10.0 * tr_xx_xzzz[i] * tbe_0 - 2.0 * tr_xxx_zzz[i] * tbe_0 + 4.0 * tr_xxx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_yyyy[i] = 2.0 * tr_0_yyyy[i] - 4.0 * tr_x_xyyyy[i] * tke_0 - 10.0 * tr_xx_yyyy[i] * tbe_0 + 4.0 * tr_xxx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_yyyz[i] = 2.0 * tr_0_yyyz[i] - 4.0 * tr_x_xyyyz[i] * tke_0 - 10.0 * tr_xx_yyyz[i] * tbe_0 + 4.0 * tr_xxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_yyzz[i] = 2.0 * tr_0_yyzz[i] - 4.0 * tr_x_xyyzz[i] * tke_0 - 10.0 * tr_xx_yyzz[i] * tbe_0 + 4.0 * tr_xxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_yzzz[i] = 2.0 * tr_0_yzzz[i] - 4.0 * tr_x_xyzzz[i] * tke_0 - 10.0 * tr_xx_yzzz[i] * tbe_0 + 4.0 * tr_xxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xx_zzzz[i] = 2.0 * tr_0_zzzz[i] - 4.0 * tr_x_xzzzz[i] * tke_0 - 10.0 * tr_xx_zzzz[i] * tbe_0 + 4.0 * tr_xxx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 15-30 components of targeted buffer : DG

    auto tr_x_0_x_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 15);

    auto tr_x_0_x_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 16);

    auto tr_x_0_x_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 17);

    auto tr_x_0_x_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 18);

    auto tr_x_0_x_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 19);

    auto tr_x_0_x_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 20);

    auto tr_x_0_x_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 21);

    auto tr_x_0_x_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 22);

    auto tr_x_0_x_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 23);

    auto tr_x_0_x_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 24);

    auto tr_x_0_x_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 25);

    auto tr_x_0_x_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 26);

    auto tr_x_0_x_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 27);

    auto tr_x_0_x_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 28);

    auto tr_x_0_x_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 29);

    #pragma omp simd aligned(tr_x_0_x_xy_xxxx, tr_x_0_x_xy_xxxy, tr_x_0_x_xy_xxxz, tr_x_0_x_xy_xxyy, tr_x_0_x_xy_xxyz, tr_x_0_x_xy_xxzz, tr_x_0_x_xy_xyyy, tr_x_0_x_xy_xyyz, tr_x_0_x_xy_xyzz, tr_x_0_x_xy_xzzz, tr_x_0_x_xy_yyyy, tr_x_0_x_xy_yyyz, tr_x_0_x_xy_yyzz, tr_x_0_x_xy_yzzz, tr_x_0_x_xy_zzzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_y_xxx, tr_y_xxxxx, tr_y_xxxxy, tr_y_xxxxz, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xy_xxxx[i] = 4.0 * tr_y_xxx[i] - 2.0 * tr_y_xxxxx[i] * tke_0 - 6.0 * tr_xy_xxxx[i] * tbe_0 - 8.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xxxy[i] = 3.0 * tr_y_xxy[i] - 2.0 * tr_y_xxxxy[i] * tke_0 - 6.0 * tr_xy_xxxy[i] * tbe_0 - 6.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xxxz[i] = 3.0 * tr_y_xxz[i] - 2.0 * tr_y_xxxxz[i] * tke_0 - 6.0 * tr_xy_xxxz[i] * tbe_0 - 6.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xxyy[i] = 2.0 * tr_y_xyy[i] - 2.0 * tr_y_xxxyy[i] * tke_0 - 6.0 * tr_xy_xxyy[i] * tbe_0 - 4.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xxyz[i] = 2.0 * tr_y_xyz[i] - 2.0 * tr_y_xxxyz[i] * tke_0 - 6.0 * tr_xy_xxyz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tbe_0 + 4.0 * tr_xxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xxzz[i] = 2.0 * tr_y_xzz[i] - 2.0 * tr_y_xxxzz[i] * tke_0 - 6.0 * tr_xy_xxzz[i] * tbe_0 - 4.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xyyy[i] = tr_y_yyy[i] - 2.0 * tr_y_xxyyy[i] * tke_0 - 6.0 * tr_xy_xyyy[i] * tbe_0 - 2.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xyyz[i] = tr_y_yyz[i] - 2.0 * tr_y_xxyyz[i] * tke_0 - 6.0 * tr_xy_xyyz[i] * tbe_0 - 2.0 * tr_xxy_yyz[i] * tbe_0 + 4.0 * tr_xxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xyzz[i] = tr_y_yzz[i] - 2.0 * tr_y_xxyzz[i] * tke_0 - 6.0 * tr_xy_xyzz[i] * tbe_0 - 2.0 * tr_xxy_yzz[i] * tbe_0 + 4.0 * tr_xxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_xzzz[i] = tr_y_zzz[i] - 2.0 * tr_y_xxzzz[i] * tke_0 - 6.0 * tr_xy_xzzz[i] * tbe_0 - 2.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_yyyy[i] = -2.0 * tr_y_xyyyy[i] * tke_0 - 6.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_yyyz[i] = -2.0 * tr_y_xyyyz[i] * tke_0 - 6.0 * tr_xy_yyyz[i] * tbe_0 + 4.0 * tr_xxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_yyzz[i] = -2.0 * tr_y_xyyzz[i] * tke_0 - 6.0 * tr_xy_yyzz[i] * tbe_0 + 4.0 * tr_xxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_yzzz[i] = -2.0 * tr_y_xyzzz[i] * tke_0 - 6.0 * tr_xy_yzzz[i] * tbe_0 + 4.0 * tr_xxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xy_zzzz[i] = -2.0 * tr_y_xzzzz[i] * tke_0 - 6.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-45 components of targeted buffer : DG

    auto tr_x_0_x_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 30);

    auto tr_x_0_x_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 31);

    auto tr_x_0_x_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 32);

    auto tr_x_0_x_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 33);

    auto tr_x_0_x_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 34);

    auto tr_x_0_x_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 35);

    auto tr_x_0_x_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 36);

    auto tr_x_0_x_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 37);

    auto tr_x_0_x_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 38);

    auto tr_x_0_x_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 39);

    auto tr_x_0_x_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 40);

    auto tr_x_0_x_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 41);

    auto tr_x_0_x_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 42);

    auto tr_x_0_x_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 43);

    auto tr_x_0_x_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 44);

    #pragma omp simd aligned(tr_x_0_x_xz_xxxx, tr_x_0_x_xz_xxxy, tr_x_0_x_xz_xxxz, tr_x_0_x_xz_xxyy, tr_x_0_x_xz_xxyz, tr_x_0_x_xz_xxzz, tr_x_0_x_xz_xyyy, tr_x_0_x_xz_xyyz, tr_x_0_x_xz_xyzz, tr_x_0_x_xz_xzzz, tr_x_0_x_xz_yyyy, tr_x_0_x_xz_yyyz, tr_x_0_x_xz_yyzz, tr_x_0_x_xz_yzzz, tr_x_0_x_xz_zzzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_z_xxx, tr_z_xxxxx, tr_z_xxxxy, tr_z_xxxxz, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xz_xxxx[i] = 4.0 * tr_z_xxx[i] - 2.0 * tr_z_xxxxx[i] * tke_0 - 6.0 * tr_xz_xxxx[i] * tbe_0 - 8.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xxxy[i] = 3.0 * tr_z_xxy[i] - 2.0 * tr_z_xxxxy[i] * tke_0 - 6.0 * tr_xz_xxxy[i] * tbe_0 - 6.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xxxz[i] = 3.0 * tr_z_xxz[i] - 2.0 * tr_z_xxxxz[i] * tke_0 - 6.0 * tr_xz_xxxz[i] * tbe_0 - 6.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xxyy[i] = 2.0 * tr_z_xyy[i] - 2.0 * tr_z_xxxyy[i] * tke_0 - 6.0 * tr_xz_xxyy[i] * tbe_0 - 4.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xxyz[i] = 2.0 * tr_z_xyz[i] - 2.0 * tr_z_xxxyz[i] * tke_0 - 6.0 * tr_xz_xxyz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tbe_0 + 4.0 * tr_xxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xxzz[i] = 2.0 * tr_z_xzz[i] - 2.0 * tr_z_xxxzz[i] * tke_0 - 6.0 * tr_xz_xxzz[i] * tbe_0 - 4.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xyyy[i] = tr_z_yyy[i] - 2.0 * tr_z_xxyyy[i] * tke_0 - 6.0 * tr_xz_xyyy[i] * tbe_0 - 2.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xyyz[i] = tr_z_yyz[i] - 2.0 * tr_z_xxyyz[i] * tke_0 - 6.0 * tr_xz_xyyz[i] * tbe_0 - 2.0 * tr_xxz_yyz[i] * tbe_0 + 4.0 * tr_xxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xyzz[i] = tr_z_yzz[i] - 2.0 * tr_z_xxyzz[i] * tke_0 - 6.0 * tr_xz_xyzz[i] * tbe_0 - 2.0 * tr_xxz_yzz[i] * tbe_0 + 4.0 * tr_xxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_xzzz[i] = tr_z_zzz[i] - 2.0 * tr_z_xxzzz[i] * tke_0 - 6.0 * tr_xz_xzzz[i] * tbe_0 - 2.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_yyyy[i] = -2.0 * tr_z_xyyyy[i] * tke_0 - 6.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_yyyz[i] = -2.0 * tr_z_xyyyz[i] * tke_0 - 6.0 * tr_xz_yyyz[i] * tbe_0 + 4.0 * tr_xxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_yyzz[i] = -2.0 * tr_z_xyyzz[i] * tke_0 - 6.0 * tr_xz_yyzz[i] * tbe_0 + 4.0 * tr_xxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_yzzz[i] = -2.0 * tr_z_xyzzz[i] * tke_0 - 6.0 * tr_xz_yzzz[i] * tbe_0 + 4.0 * tr_xxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xz_zzzz[i] = -2.0 * tr_z_xzzzz[i] * tke_0 - 6.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 45-60 components of targeted buffer : DG

    auto tr_x_0_x_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 45);

    auto tr_x_0_x_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 46);

    auto tr_x_0_x_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 47);

    auto tr_x_0_x_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 48);

    auto tr_x_0_x_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 49);

    auto tr_x_0_x_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 50);

    auto tr_x_0_x_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 51);

    auto tr_x_0_x_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 52);

    auto tr_x_0_x_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 53);

    auto tr_x_0_x_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 54);

    auto tr_x_0_x_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 55);

    auto tr_x_0_x_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 56);

    auto tr_x_0_x_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 57);

    auto tr_x_0_x_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 58);

    auto tr_x_0_x_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 59);

    #pragma omp simd aligned(tr_x_0_x_yy_xxxx, tr_x_0_x_yy_xxxy, tr_x_0_x_yy_xxxz, tr_x_0_x_yy_xxyy, tr_x_0_x_yy_xxyz, tr_x_0_x_yy_xxzz, tr_x_0_x_yy_xyyy, tr_x_0_x_yy_xyyz, tr_x_0_x_yy_xyzz, tr_x_0_x_yy_xzzz, tr_x_0_x_yy_yyyy, tr_x_0_x_yy_yyyz, tr_x_0_x_yy_yyzz, tr_x_0_x_yy_yzzz, tr_x_0_x_yy_zzzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yy_xxxx[i] = -2.0 * tr_yy_xxxx[i] * tbe_0 - 8.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xxxy[i] = -2.0 * tr_yy_xxxy[i] * tbe_0 - 6.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xxxz[i] = -2.0 * tr_yy_xxxz[i] * tbe_0 - 6.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xxyy[i] = -2.0 * tr_yy_xxyy[i] * tbe_0 - 4.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xxyz[i] = -2.0 * tr_yy_xxyz[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tbe_0 + 4.0 * tr_xyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xxzz[i] = -2.0 * tr_yy_xxzz[i] * tbe_0 - 4.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xyyy[i] = -2.0 * tr_yy_xyyy[i] * tbe_0 - 2.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xyyz[i] = -2.0 * tr_yy_xyyz[i] * tbe_0 - 2.0 * tr_xyy_yyz[i] * tbe_0 + 4.0 * tr_xyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xyzz[i] = -2.0 * tr_yy_xyzz[i] * tbe_0 - 2.0 * tr_xyy_yzz[i] * tbe_0 + 4.0 * tr_xyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_xzzz[i] = -2.0 * tr_yy_xzzz[i] * tbe_0 - 2.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_yyyy[i] = -2.0 * tr_yy_yyyy[i] * tbe_0 + 4.0 * tr_xyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_yyyz[i] = -2.0 * tr_yy_yyyz[i] * tbe_0 + 4.0 * tr_xyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_yyzz[i] = -2.0 * tr_yy_yyzz[i] * tbe_0 + 4.0 * tr_xyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_yzzz[i] = -2.0 * tr_yy_yzzz[i] * tbe_0 + 4.0 * tr_xyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yy_zzzz[i] = -2.0 * tr_yy_zzzz[i] * tbe_0 + 4.0 * tr_xyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-75 components of targeted buffer : DG

    auto tr_x_0_x_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 60);

    auto tr_x_0_x_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 61);

    auto tr_x_0_x_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 62);

    auto tr_x_0_x_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 63);

    auto tr_x_0_x_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 64);

    auto tr_x_0_x_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 65);

    auto tr_x_0_x_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 66);

    auto tr_x_0_x_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 67);

    auto tr_x_0_x_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 68);

    auto tr_x_0_x_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 69);

    auto tr_x_0_x_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 70);

    auto tr_x_0_x_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 71);

    auto tr_x_0_x_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 72);

    auto tr_x_0_x_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 73);

    auto tr_x_0_x_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 74);

    #pragma omp simd aligned(tr_x_0_x_yz_xxxx, tr_x_0_x_yz_xxxy, tr_x_0_x_yz_xxxz, tr_x_0_x_yz_xxyy, tr_x_0_x_yz_xxyz, tr_x_0_x_yz_xxzz, tr_x_0_x_yz_xyyy, tr_x_0_x_yz_xyyz, tr_x_0_x_yz_xyzz, tr_x_0_x_yz_xzzz, tr_x_0_x_yz_yyyy, tr_x_0_x_yz_yyyz, tr_x_0_x_yz_yyzz, tr_x_0_x_yz_yzzz, tr_x_0_x_yz_zzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yz_xxxx[i] = -2.0 * tr_yz_xxxx[i] * tbe_0 - 8.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xxxy[i] = -2.0 * tr_yz_xxxy[i] * tbe_0 - 6.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xxxz[i] = -2.0 * tr_yz_xxxz[i] * tbe_0 - 6.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xxyy[i] = -2.0 * tr_yz_xxyy[i] * tbe_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xxyz[i] = -2.0 * tr_yz_xxyz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xxzz[i] = -2.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xyyy[i] = -2.0 * tr_yz_xyyy[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xyyz[i] = -2.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xyzz[i] = -2.0 * tr_yz_xyzz[i] * tbe_0 - 2.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_xzzz[i] = -2.0 * tr_yz_xzzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_yyyy[i] = -2.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_yyyz[i] = -2.0 * tr_yz_yyyz[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_yyzz[i] = -2.0 * tr_yz_yyzz[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_yzzz[i] = -2.0 * tr_yz_yzzz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yz_zzzz[i] = -2.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_xyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 75-90 components of targeted buffer : DG

    auto tr_x_0_x_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 75);

    auto tr_x_0_x_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 76);

    auto tr_x_0_x_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 77);

    auto tr_x_0_x_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 78);

    auto tr_x_0_x_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 79);

    auto tr_x_0_x_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 80);

    auto tr_x_0_x_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 81);

    auto tr_x_0_x_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 82);

    auto tr_x_0_x_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 83);

    auto tr_x_0_x_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 84);

    auto tr_x_0_x_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 85);

    auto tr_x_0_x_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 86);

    auto tr_x_0_x_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 87);

    auto tr_x_0_x_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 88);

    auto tr_x_0_x_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 89);

    #pragma omp simd aligned(tr_x_0_x_zz_xxxx, tr_x_0_x_zz_xxxy, tr_x_0_x_zz_xxxz, tr_x_0_x_zz_xxyy, tr_x_0_x_zz_xxyz, tr_x_0_x_zz_xxzz, tr_x_0_x_zz_xyyy, tr_x_0_x_zz_xyyz, tr_x_0_x_zz_xyzz, tr_x_0_x_zz_xzzz, tr_x_0_x_zz_yyyy, tr_x_0_x_zz_yyyz, tr_x_0_x_zz_yyzz, tr_x_0_x_zz_yzzz, tr_x_0_x_zz_zzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zz_xxxx[i] = -2.0 * tr_zz_xxxx[i] * tbe_0 - 8.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xxxy[i] = -2.0 * tr_zz_xxxy[i] * tbe_0 - 6.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xxxz[i] = -2.0 * tr_zz_xxxz[i] * tbe_0 - 6.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xxyy[i] = -2.0 * tr_zz_xxyy[i] * tbe_0 - 4.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xxyz[i] = -2.0 * tr_zz_xxyz[i] * tbe_0 - 4.0 * tr_xzz_xyz[i] * tbe_0 + 4.0 * tr_xzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xxzz[i] = -2.0 * tr_zz_xxzz[i] * tbe_0 - 4.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xyyy[i] = -2.0 * tr_zz_xyyy[i] * tbe_0 - 2.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xyyz[i] = -2.0 * tr_zz_xyyz[i] * tbe_0 - 2.0 * tr_xzz_yyz[i] * tbe_0 + 4.0 * tr_xzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xyzz[i] = -2.0 * tr_zz_xyzz[i] * tbe_0 - 2.0 * tr_xzz_yzz[i] * tbe_0 + 4.0 * tr_xzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_xzzz[i] = -2.0 * tr_zz_xzzz[i] * tbe_0 - 2.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_yyyy[i] = -2.0 * tr_zz_yyyy[i] * tbe_0 + 4.0 * tr_xzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_yyyz[i] = -2.0 * tr_zz_yyyz[i] * tbe_0 + 4.0 * tr_xzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_yyzz[i] = -2.0 * tr_zz_yyzz[i] * tbe_0 + 4.0 * tr_xzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_yzzz[i] = -2.0 * tr_zz_yzzz[i] * tbe_0 + 4.0 * tr_xzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zz_zzzz[i] = -2.0 * tr_zz_zzzz[i] * tbe_0 + 4.0 * tr_xzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-105 components of targeted buffer : DG

    auto tr_x_0_y_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 90);

    auto tr_x_0_y_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 91);

    auto tr_x_0_y_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 92);

    auto tr_x_0_y_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 93);

    auto tr_x_0_y_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 94);

    auto tr_x_0_y_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 95);

    auto tr_x_0_y_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 96);

    auto tr_x_0_y_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 97);

    auto tr_x_0_y_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 98);

    auto tr_x_0_y_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 99);

    auto tr_x_0_y_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 100);

    auto tr_x_0_y_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 101);

    auto tr_x_0_y_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 102);

    auto tr_x_0_y_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 103);

    auto tr_x_0_y_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 104);

    #pragma omp simd aligned(tr_x_0_y_xx_xxxx, tr_x_0_y_xx_xxxy, tr_x_0_y_xx_xxxz, tr_x_0_y_xx_xxyy, tr_x_0_y_xx_xxyz, tr_x_0_y_xx_xxzz, tr_x_0_y_xx_xyyy, tr_x_0_y_xx_xyyz, tr_x_0_y_xx_xyzz, tr_x_0_y_xx_xzzz, tr_x_0_y_xx_yyyy, tr_x_0_y_xx_yyyz, tr_x_0_y_xx_yyzz, tr_x_0_y_xx_yzzz, tr_x_0_y_xx_zzzz, tr_x_xxx, tr_x_xxxxy, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_yyy, tr_x_yyyyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_xxx_xxx, tr_xxx_xxxxy, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyyyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xx_xxxx[i] = -4.0 * tr_x_xxxxy[i] * tke_0 - 4.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xxx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xxxy[i] = 2.0 * tr_x_xxx[i] - 4.0 * tr_x_xxxyy[i] * tke_0 - 4.0 * tr_xy_xxxy[i] * tbe_0 - 2.0 * tr_xxx_xxx[i] * tbe_0 + 4.0 * tr_xxx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xxxz[i] = -4.0 * tr_x_xxxyz[i] * tke_0 - 4.0 * tr_xy_xxxz[i] * tbe_0 + 4.0 * tr_xxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xxyy[i] = 4.0 * tr_x_xxy[i] - 4.0 * tr_x_xxyyy[i] * tke_0 - 4.0 * tr_xy_xxyy[i] * tbe_0 - 4.0 * tr_xxx_xxy[i] * tbe_0 + 4.0 * tr_xxx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xxyz[i] = 2.0 * tr_x_xxz[i] - 4.0 * tr_x_xxyyz[i] * tke_0 - 4.0 * tr_xy_xxyz[i] * tbe_0 - 2.0 * tr_xxx_xxz[i] * tbe_0 + 4.0 * tr_xxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xxzz[i] = -4.0 * tr_x_xxyzz[i] * tke_0 - 4.0 * tr_xy_xxzz[i] * tbe_0 + 4.0 * tr_xxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xyyy[i] = 6.0 * tr_x_xyy[i] - 4.0 * tr_x_xyyyy[i] * tke_0 - 4.0 * tr_xy_xyyy[i] * tbe_0 - 6.0 * tr_xxx_xyy[i] * tbe_0 + 4.0 * tr_xxx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xyyz[i] = 4.0 * tr_x_xyz[i] - 4.0 * tr_x_xyyyz[i] * tke_0 - 4.0 * tr_xy_xyyz[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tbe_0 + 4.0 * tr_xxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xyzz[i] = 2.0 * tr_x_xzz[i] - 4.0 * tr_x_xyyzz[i] * tke_0 - 4.0 * tr_xy_xyzz[i] * tbe_0 - 2.0 * tr_xxx_xzz[i] * tbe_0 + 4.0 * tr_xxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_xzzz[i] = -4.0 * tr_x_xyzzz[i] * tke_0 - 4.0 * tr_xy_xzzz[i] * tbe_0 + 4.0 * tr_xxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_yyyy[i] = 8.0 * tr_x_yyy[i] - 4.0 * tr_x_yyyyy[i] * tke_0 - 4.0 * tr_xy_yyyy[i] * tbe_0 - 8.0 * tr_xxx_yyy[i] * tbe_0 + 4.0 * tr_xxx_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_yyyz[i] = 6.0 * tr_x_yyz[i] - 4.0 * tr_x_yyyyz[i] * tke_0 - 4.0 * tr_xy_yyyz[i] * tbe_0 - 6.0 * tr_xxx_yyz[i] * tbe_0 + 4.0 * tr_xxx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_yyzz[i] = 4.0 * tr_x_yzz[i] - 4.0 * tr_x_yyyzz[i] * tke_0 - 4.0 * tr_xy_yyzz[i] * tbe_0 - 4.0 * tr_xxx_yzz[i] * tbe_0 + 4.0 * tr_xxx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_yzzz[i] = 2.0 * tr_x_zzz[i] - 4.0 * tr_x_yyzzz[i] * tke_0 - 4.0 * tr_xy_yzzz[i] * tbe_0 - 2.0 * tr_xxx_zzz[i] * tbe_0 + 4.0 * tr_xxx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xx_zzzz[i] = -4.0 * tr_x_yzzzz[i] * tke_0 - 4.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xxx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 105-120 components of targeted buffer : DG

    auto tr_x_0_y_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 105);

    auto tr_x_0_y_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 106);

    auto tr_x_0_y_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 107);

    auto tr_x_0_y_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 108);

    auto tr_x_0_y_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 109);

    auto tr_x_0_y_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 110);

    auto tr_x_0_y_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 111);

    auto tr_x_0_y_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 112);

    auto tr_x_0_y_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 113);

    auto tr_x_0_y_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 114);

    auto tr_x_0_y_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 115);

    auto tr_x_0_y_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 116);

    auto tr_x_0_y_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 117);

    auto tr_x_0_y_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 118);

    auto tr_x_0_y_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 119);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_x_0_y_xy_xxxx, tr_x_0_y_xy_xxxy, tr_x_0_y_xy_xxxz, tr_x_0_y_xy_xxyy, tr_x_0_y_xy_xxyz, tr_x_0_y_xy_xxzz, tr_x_0_y_xy_xyyy, tr_x_0_y_xy_xyyz, tr_x_0_y_xy_xyzz, tr_x_0_y_xy_xzzz, tr_x_0_y_xy_yyyy, tr_x_0_y_xy_yyyz, tr_x_0_y_xy_yyzz, tr_x_0_y_xy_yzzz, tr_x_0_y_xy_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_y_xxx, tr_y_xxxxy, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_yyy, tr_y_yyyyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xy_xxxx[i] = tr_0_xxxx[i] - 2.0 * tr_y_xxxxy[i] * tke_0 - 2.0 * tr_yy_xxxx[i] * tbe_0 - 2.0 * tr_xx_xxxx[i] * tbe_0 + 4.0 * tr_xxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xxxy[i] = tr_0_xxxy[i] + tr_y_xxx[i] - 2.0 * tr_y_xxxyy[i] * tke_0 - 2.0 * tr_yy_xxxy[i] * tbe_0 - 2.0 * tr_xx_xxxy[i] * tbe_0 - 2.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xxxz[i] = tr_0_xxxz[i] - 2.0 * tr_y_xxxyz[i] * tke_0 - 2.0 * tr_yy_xxxz[i] * tbe_0 - 2.0 * tr_xx_xxxz[i] * tbe_0 + 4.0 * tr_xxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xxyy[i] = tr_0_xxyy[i] + 2.0 * tr_y_xxy[i] - 2.0 * tr_y_xxyyy[i] * tke_0 - 2.0 * tr_yy_xxyy[i] * tbe_0 - 2.0 * tr_xx_xxyy[i] * tbe_0 - 4.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xxyz[i] = tr_0_xxyz[i] + tr_y_xxz[i] - 2.0 * tr_y_xxyyz[i] * tke_0 - 2.0 * tr_yy_xxyz[i] * tbe_0 - 2.0 * tr_xx_xxyz[i] * tbe_0 - 2.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xxzz[i] = tr_0_xxzz[i] - 2.0 * tr_y_xxyzz[i] * tke_0 - 2.0 * tr_yy_xxzz[i] * tbe_0 - 2.0 * tr_xx_xxzz[i] * tbe_0 + 4.0 * tr_xxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xyyy[i] = tr_0_xyyy[i] + 3.0 * tr_y_xyy[i] - 2.0 * tr_y_xyyyy[i] * tke_0 - 2.0 * tr_yy_xyyy[i] * tbe_0 - 2.0 * tr_xx_xyyy[i] * tbe_0 - 6.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xyyz[i] = tr_0_xyyz[i] + 2.0 * tr_y_xyz[i] - 2.0 * tr_y_xyyyz[i] * tke_0 - 2.0 * tr_yy_xyyz[i] * tbe_0 - 2.0 * tr_xx_xyyz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tbe_0 + 4.0 * tr_xxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xyzz[i] = tr_0_xyzz[i] + tr_y_xzz[i] - 2.0 * tr_y_xyyzz[i] * tke_0 - 2.0 * tr_yy_xyzz[i] * tbe_0 - 2.0 * tr_xx_xyzz[i] * tbe_0 - 2.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_xzzz[i] = tr_0_xzzz[i] - 2.0 * tr_y_xyzzz[i] * tke_0 - 2.0 * tr_yy_xzzz[i] * tbe_0 - 2.0 * tr_xx_xzzz[i] * tbe_0 + 4.0 * tr_xxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_yyyy[i] = tr_0_yyyy[i] + 4.0 * tr_y_yyy[i] - 2.0 * tr_y_yyyyy[i] * tke_0 - 2.0 * tr_yy_yyyy[i] * tbe_0 - 2.0 * tr_xx_yyyy[i] * tbe_0 - 8.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_yyyz[i] = tr_0_yyyz[i] + 3.0 * tr_y_yyz[i] - 2.0 * tr_y_yyyyz[i] * tke_0 - 2.0 * tr_yy_yyyz[i] * tbe_0 - 2.0 * tr_xx_yyyz[i] * tbe_0 - 6.0 * tr_xxy_yyz[i] * tbe_0 + 4.0 * tr_xxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_yyzz[i] = tr_0_yyzz[i] + 2.0 * tr_y_yzz[i] - 2.0 * tr_y_yyyzz[i] * tke_0 - 2.0 * tr_yy_yyzz[i] * tbe_0 - 2.0 * tr_xx_yyzz[i] * tbe_0 - 4.0 * tr_xxy_yzz[i] * tbe_0 + 4.0 * tr_xxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_yzzz[i] = tr_0_yzzz[i] + tr_y_zzz[i] - 2.0 * tr_y_yyzzz[i] * tke_0 - 2.0 * tr_yy_yzzz[i] * tbe_0 - 2.0 * tr_xx_yzzz[i] * tbe_0 - 2.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xy_zzzz[i] = tr_0_zzzz[i] - 2.0 * tr_y_yzzzz[i] * tke_0 - 2.0 * tr_yy_zzzz[i] * tbe_0 - 2.0 * tr_xx_zzzz[i] * tbe_0 + 4.0 * tr_xxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-135 components of targeted buffer : DG

    auto tr_x_0_y_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 120);

    auto tr_x_0_y_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 121);

    auto tr_x_0_y_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 122);

    auto tr_x_0_y_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 123);

    auto tr_x_0_y_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 124);

    auto tr_x_0_y_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 125);

    auto tr_x_0_y_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 126);

    auto tr_x_0_y_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 127);

    auto tr_x_0_y_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 128);

    auto tr_x_0_y_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 129);

    auto tr_x_0_y_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 130);

    auto tr_x_0_y_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 131);

    auto tr_x_0_y_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 132);

    auto tr_x_0_y_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 133);

    auto tr_x_0_y_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 134);

    #pragma omp simd aligned(tr_x_0_y_xz_xxxx, tr_x_0_y_xz_xxxy, tr_x_0_y_xz_xxxz, tr_x_0_y_xz_xxyy, tr_x_0_y_xz_xxyz, tr_x_0_y_xz_xxzz, tr_x_0_y_xz_xyyy, tr_x_0_y_xz_xyyz, tr_x_0_y_xz_xyzz, tr_x_0_y_xz_xzzz, tr_x_0_y_xz_yyyy, tr_x_0_y_xz_yyyz, tr_x_0_y_xz_yyzz, tr_x_0_y_xz_yzzz, tr_x_0_y_xz_zzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_z_xxx, tr_z_xxxxy, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_yyy, tr_z_yyyyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xz_xxxx[i] = -2.0 * tr_z_xxxxy[i] * tke_0 - 2.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_xxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xxxy[i] = tr_z_xxx[i] - 2.0 * tr_z_xxxyy[i] * tke_0 - 2.0 * tr_yz_xxxy[i] * tbe_0 - 2.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xxxz[i] = -2.0 * tr_z_xxxyz[i] * tke_0 - 2.0 * tr_yz_xxxz[i] * tbe_0 + 4.0 * tr_xxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xxyy[i] = 2.0 * tr_z_xxy[i] - 2.0 * tr_z_xxyyy[i] * tke_0 - 2.0 * tr_yz_xxyy[i] * tbe_0 - 4.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xxyz[i] = tr_z_xxz[i] - 2.0 * tr_z_xxyyz[i] * tke_0 - 2.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xxzz[i] = -2.0 * tr_z_xxyzz[i] * tke_0 - 2.0 * tr_yz_xxzz[i] * tbe_0 + 4.0 * tr_xxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xyyy[i] = 3.0 * tr_z_xyy[i] - 2.0 * tr_z_xyyyy[i] * tke_0 - 2.0 * tr_yz_xyyy[i] * tbe_0 - 6.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xyyz[i] = 2.0 * tr_z_xyz[i] - 2.0 * tr_z_xyyyz[i] * tke_0 - 2.0 * tr_yz_xyyz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tbe_0 + 4.0 * tr_xxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xyzz[i] = tr_z_xzz[i] - 2.0 * tr_z_xyyzz[i] * tke_0 - 2.0 * tr_yz_xyzz[i] * tbe_0 - 2.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_xzzz[i] = -2.0 * tr_z_xyzzz[i] * tke_0 - 2.0 * tr_yz_xzzz[i] * tbe_0 + 4.0 * tr_xxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_yyyy[i] = 4.0 * tr_z_yyy[i] - 2.0 * tr_z_yyyyy[i] * tke_0 - 2.0 * tr_yz_yyyy[i] * tbe_0 - 8.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_yyyz[i] = 3.0 * tr_z_yyz[i] - 2.0 * tr_z_yyyyz[i] * tke_0 - 2.0 * tr_yz_yyyz[i] * tbe_0 - 6.0 * tr_xxz_yyz[i] * tbe_0 + 4.0 * tr_xxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_yyzz[i] = 2.0 * tr_z_yzz[i] - 2.0 * tr_z_yyyzz[i] * tke_0 - 2.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_xxz_yzz[i] * tbe_0 + 4.0 * tr_xxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_yzzz[i] = tr_z_zzz[i] - 2.0 * tr_z_yyzzz[i] * tke_0 - 2.0 * tr_yz_yzzz[i] * tbe_0 - 2.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xz_zzzz[i] = -2.0 * tr_z_yzzzz[i] * tke_0 - 2.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_xxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 135-150 components of targeted buffer : DG

    auto tr_x_0_y_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 135);

    auto tr_x_0_y_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 136);

    auto tr_x_0_y_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 137);

    auto tr_x_0_y_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 138);

    auto tr_x_0_y_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 139);

    auto tr_x_0_y_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 140);

    auto tr_x_0_y_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 141);

    auto tr_x_0_y_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 142);

    auto tr_x_0_y_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 143);

    auto tr_x_0_y_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 144);

    auto tr_x_0_y_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 145);

    auto tr_x_0_y_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 146);

    auto tr_x_0_y_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 147);

    auto tr_x_0_y_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 148);

    auto tr_x_0_y_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 149);

    #pragma omp simd aligned(tr_x_0_y_yy_xxxx, tr_x_0_y_yy_xxxy, tr_x_0_y_yy_xxxz, tr_x_0_y_yy_xxyy, tr_x_0_y_yy_xxyz, tr_x_0_y_yy_xxzz, tr_x_0_y_yy_xyyy, tr_x_0_y_yy_xyyz, tr_x_0_y_yy_xyzz, tr_x_0_y_yy_xzzz, tr_x_0_y_yy_yyyy, tr_x_0_y_yy_yyyz, tr_x_0_y_yy_yyzz, tr_x_0_y_yy_yzzz, tr_x_0_y_yy_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yy_xxxx[i] = -4.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xxxy[i] = -4.0 * tr_xy_xxxy[i] * tbe_0 - 2.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xxxz[i] = -4.0 * tr_xy_xxxz[i] * tbe_0 + 4.0 * tr_xyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xxyy[i] = -4.0 * tr_xy_xxyy[i] * tbe_0 - 4.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xxyz[i] = -4.0 * tr_xy_xxyz[i] * tbe_0 - 2.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xxzz[i] = -4.0 * tr_xy_xxzz[i] * tbe_0 + 4.0 * tr_xyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xyyy[i] = -4.0 * tr_xy_xyyy[i] * tbe_0 - 6.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xyyz[i] = -4.0 * tr_xy_xyyz[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tbe_0 + 4.0 * tr_xyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xyzz[i] = -4.0 * tr_xy_xyzz[i] * tbe_0 - 2.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_xzzz[i] = -4.0 * tr_xy_xzzz[i] * tbe_0 + 4.0 * tr_xyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_yyyy[i] = -4.0 * tr_xy_yyyy[i] * tbe_0 - 8.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_yyyz[i] = -4.0 * tr_xy_yyyz[i] * tbe_0 - 6.0 * tr_xyy_yyz[i] * tbe_0 + 4.0 * tr_xyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_yyzz[i] = -4.0 * tr_xy_yyzz[i] * tbe_0 - 4.0 * tr_xyy_yzz[i] * tbe_0 + 4.0 * tr_xyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_yzzz[i] = -4.0 * tr_xy_yzzz[i] * tbe_0 - 2.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yy_zzzz[i] = -4.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-165 components of targeted buffer : DG

    auto tr_x_0_y_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 150);

    auto tr_x_0_y_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 151);

    auto tr_x_0_y_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 152);

    auto tr_x_0_y_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 153);

    auto tr_x_0_y_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 154);

    auto tr_x_0_y_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 155);

    auto tr_x_0_y_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 156);

    auto tr_x_0_y_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 157);

    auto tr_x_0_y_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 158);

    auto tr_x_0_y_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 159);

    auto tr_x_0_y_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 160);

    auto tr_x_0_y_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 161);

    auto tr_x_0_y_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 162);

    auto tr_x_0_y_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 163);

    auto tr_x_0_y_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 164);

    #pragma omp simd aligned(tr_x_0_y_yz_xxxx, tr_x_0_y_yz_xxxy, tr_x_0_y_yz_xxxz, tr_x_0_y_yz_xxyy, tr_x_0_y_yz_xxyz, tr_x_0_y_yz_xxzz, tr_x_0_y_yz_xyyy, tr_x_0_y_yz_xyyz, tr_x_0_y_yz_xyzz, tr_x_0_y_yz_xzzz, tr_x_0_y_yz_yyyy, tr_x_0_y_yz_yyyz, tr_x_0_y_yz_yyzz, tr_x_0_y_yz_yzzz, tr_x_0_y_yz_zzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yz_xxxx[i] = -2.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xxxy[i] = -2.0 * tr_xz_xxxy[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xxxz[i] = -2.0 * tr_xz_xxxz[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xxyy[i] = -2.0 * tr_xz_xxyy[i] * tbe_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xxyz[i] = -2.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xxzz[i] = -2.0 * tr_xz_xxzz[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xyyy[i] = -2.0 * tr_xz_xyyy[i] * tbe_0 - 6.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xyyz[i] = -2.0 * tr_xz_xyyz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xyzz[i] = -2.0 * tr_xz_xyzz[i] * tbe_0 - 2.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_xzzz[i] = -2.0 * tr_xz_xzzz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_yyyy[i] = -2.0 * tr_xz_yyyy[i] * tbe_0 - 8.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_yyyz[i] = -2.0 * tr_xz_yyyz[i] * tbe_0 - 6.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_yyzz[i] = -2.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_yzzz[i] = -2.0 * tr_xz_yzzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yz_zzzz[i] = -2.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 165-180 components of targeted buffer : DG

    auto tr_x_0_y_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 165);

    auto tr_x_0_y_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 166);

    auto tr_x_0_y_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 167);

    auto tr_x_0_y_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 168);

    auto tr_x_0_y_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 169);

    auto tr_x_0_y_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 170);

    auto tr_x_0_y_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 171);

    auto tr_x_0_y_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 172);

    auto tr_x_0_y_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 173);

    auto tr_x_0_y_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 174);

    auto tr_x_0_y_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 175);

    auto tr_x_0_y_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 176);

    auto tr_x_0_y_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 177);

    auto tr_x_0_y_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 178);

    auto tr_x_0_y_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 179);

    #pragma omp simd aligned(tr_x_0_y_zz_xxxx, tr_x_0_y_zz_xxxy, tr_x_0_y_zz_xxxz, tr_x_0_y_zz_xxyy, tr_x_0_y_zz_xxyz, tr_x_0_y_zz_xxzz, tr_x_0_y_zz_xyyy, tr_x_0_y_zz_xyyz, tr_x_0_y_zz_xyzz, tr_x_0_y_zz_xzzz, tr_x_0_y_zz_yyyy, tr_x_0_y_zz_yyyz, tr_x_0_y_zz_yyzz, tr_x_0_y_zz_yzzz, tr_x_0_y_zz_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zz_xxxx[i] = 4.0 * tr_xzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xxxy[i] = -2.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xxxz[i] = 4.0 * tr_xzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xxyy[i] = -4.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xxyz[i] = -2.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xxzz[i] = 4.0 * tr_xzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xyyy[i] = -6.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xyyz[i] = -4.0 * tr_xzz_xyz[i] * tbe_0 + 4.0 * tr_xzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xyzz[i] = -2.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_xzzz[i] = 4.0 * tr_xzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_yyyy[i] = -8.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_yyyz[i] = -6.0 * tr_xzz_yyz[i] * tbe_0 + 4.0 * tr_xzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_yyzz[i] = -4.0 * tr_xzz_yzz[i] * tbe_0 + 4.0 * tr_xzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_yzzz[i] = -2.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zz_zzzz[i] = 4.0 * tr_xzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-195 components of targeted buffer : DG

    auto tr_x_0_z_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 180);

    auto tr_x_0_z_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 181);

    auto tr_x_0_z_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 182);

    auto tr_x_0_z_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 183);

    auto tr_x_0_z_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 184);

    auto tr_x_0_z_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 185);

    auto tr_x_0_z_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 186);

    auto tr_x_0_z_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 187);

    auto tr_x_0_z_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 188);

    auto tr_x_0_z_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 189);

    auto tr_x_0_z_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 190);

    auto tr_x_0_z_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 191);

    auto tr_x_0_z_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 192);

    auto tr_x_0_z_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 193);

    auto tr_x_0_z_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 194);

    #pragma omp simd aligned(tr_x_0_z_xx_xxxx, tr_x_0_z_xx_xxxy, tr_x_0_z_xx_xxxz, tr_x_0_z_xx_xxyy, tr_x_0_z_xx_xxyz, tr_x_0_z_xx_xxzz, tr_x_0_z_xx_xyyy, tr_x_0_z_xx_xyyz, tr_x_0_z_xx_xyzz, tr_x_0_z_xx_xzzz, tr_x_0_z_xx_yyyy, tr_x_0_z_xx_yyyz, tr_x_0_z_xx_yyzz, tr_x_0_z_xx_yzzz, tr_x_0_z_xx_zzzz, tr_x_xxx, tr_x_xxxxz, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_x_zzzzz, tr_xxx_xxx, tr_xxx_xxxxz, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxx_zzzzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xx_xxxx[i] = -4.0 * tr_x_xxxxz[i] * tke_0 - 4.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xxx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xxxy[i] = -4.0 * tr_x_xxxyz[i] * tke_0 - 4.0 * tr_xz_xxxy[i] * tbe_0 + 4.0 * tr_xxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xxxz[i] = 2.0 * tr_x_xxx[i] - 4.0 * tr_x_xxxzz[i] * tke_0 - 4.0 * tr_xz_xxxz[i] * tbe_0 - 2.0 * tr_xxx_xxx[i] * tbe_0 + 4.0 * tr_xxx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xxyy[i] = -4.0 * tr_x_xxyyz[i] * tke_0 - 4.0 * tr_xz_xxyy[i] * tbe_0 + 4.0 * tr_xxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xxyz[i] = 2.0 * tr_x_xxy[i] - 4.0 * tr_x_xxyzz[i] * tke_0 - 4.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xxx_xxy[i] * tbe_0 + 4.0 * tr_xxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xxzz[i] = 4.0 * tr_x_xxz[i] - 4.0 * tr_x_xxzzz[i] * tke_0 - 4.0 * tr_xz_xxzz[i] * tbe_0 - 4.0 * tr_xxx_xxz[i] * tbe_0 + 4.0 * tr_xxx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xyyy[i] = -4.0 * tr_x_xyyyz[i] * tke_0 - 4.0 * tr_xz_xyyy[i] * tbe_0 + 4.0 * tr_xxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xyyz[i] = 2.0 * tr_x_xyy[i] - 4.0 * tr_x_xyyzz[i] * tke_0 - 4.0 * tr_xz_xyyz[i] * tbe_0 - 2.0 * tr_xxx_xyy[i] * tbe_0 + 4.0 * tr_xxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xyzz[i] = 4.0 * tr_x_xyz[i] - 4.0 * tr_x_xyzzz[i] * tke_0 - 4.0 * tr_xz_xyzz[i] * tbe_0 - 4.0 * tr_xxx_xyz[i] * tbe_0 + 4.0 * tr_xxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_xzzz[i] = 6.0 * tr_x_xzz[i] - 4.0 * tr_x_xzzzz[i] * tke_0 - 4.0 * tr_xz_xzzz[i] * tbe_0 - 6.0 * tr_xxx_xzz[i] * tbe_0 + 4.0 * tr_xxx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_yyyy[i] = -4.0 * tr_x_yyyyz[i] * tke_0 - 4.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xxx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_yyyz[i] = 2.0 * tr_x_yyy[i] - 4.0 * tr_x_yyyzz[i] * tke_0 - 4.0 * tr_xz_yyyz[i] * tbe_0 - 2.0 * tr_xxx_yyy[i] * tbe_0 + 4.0 * tr_xxx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_yyzz[i] = 4.0 * tr_x_yyz[i] - 4.0 * tr_x_yyzzz[i] * tke_0 - 4.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xxx_yyz[i] * tbe_0 + 4.0 * tr_xxx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_yzzz[i] = 6.0 * tr_x_yzz[i] - 4.0 * tr_x_yzzzz[i] * tke_0 - 4.0 * tr_xz_yzzz[i] * tbe_0 - 6.0 * tr_xxx_yzz[i] * tbe_0 + 4.0 * tr_xxx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xx_zzzz[i] = 8.0 * tr_x_zzz[i] - 4.0 * tr_x_zzzzz[i] * tke_0 - 4.0 * tr_xz_zzzz[i] * tbe_0 - 8.0 * tr_xxx_zzz[i] * tbe_0 + 4.0 * tr_xxx_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 195-210 components of targeted buffer : DG

    auto tr_x_0_z_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 195);

    auto tr_x_0_z_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 196);

    auto tr_x_0_z_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 197);

    auto tr_x_0_z_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 198);

    auto tr_x_0_z_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 199);

    auto tr_x_0_z_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 200);

    auto tr_x_0_z_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 201);

    auto tr_x_0_z_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 202);

    auto tr_x_0_z_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 203);

    auto tr_x_0_z_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 204);

    auto tr_x_0_z_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 205);

    auto tr_x_0_z_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 206);

    auto tr_x_0_z_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 207);

    auto tr_x_0_z_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 208);

    auto tr_x_0_z_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 209);

    #pragma omp simd aligned(tr_x_0_z_xy_xxxx, tr_x_0_z_xy_xxxy, tr_x_0_z_xy_xxxz, tr_x_0_z_xy_xxyy, tr_x_0_z_xy_xxyz, tr_x_0_z_xy_xxzz, tr_x_0_z_xy_xyyy, tr_x_0_z_xy_xyyz, tr_x_0_z_xy_xyzz, tr_x_0_z_xy_xzzz, tr_x_0_z_xy_yyyy, tr_x_0_z_xy_yyyz, tr_x_0_z_xy_yyzz, tr_x_0_z_xy_yzzz, tr_x_0_z_xy_zzzz, tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_y_xxx, tr_y_xxxxz, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_y_zzzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xy_xxxx[i] = -2.0 * tr_y_xxxxz[i] * tke_0 - 2.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_xxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xxxy[i] = -2.0 * tr_y_xxxyz[i] * tke_0 - 2.0 * tr_yz_xxxy[i] * tbe_0 + 4.0 * tr_xxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xxxz[i] = tr_y_xxx[i] - 2.0 * tr_y_xxxzz[i] * tke_0 - 2.0 * tr_yz_xxxz[i] * tbe_0 - 2.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xxyy[i] = -2.0 * tr_y_xxyyz[i] * tke_0 - 2.0 * tr_yz_xxyy[i] * tbe_0 + 4.0 * tr_xxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xxyz[i] = tr_y_xxy[i] - 2.0 * tr_y_xxyzz[i] * tke_0 - 2.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xxzz[i] = 2.0 * tr_y_xxz[i] - 2.0 * tr_y_xxzzz[i] * tke_0 - 2.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xyyy[i] = -2.0 * tr_y_xyyyz[i] * tke_0 - 2.0 * tr_yz_xyyy[i] * tbe_0 + 4.0 * tr_xxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xyyz[i] = tr_y_xyy[i] - 2.0 * tr_y_xyyzz[i] * tke_0 - 2.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xyzz[i] = 2.0 * tr_y_xyz[i] - 2.0 * tr_y_xyzzz[i] * tke_0 - 2.0 * tr_yz_xyzz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tbe_0 + 4.0 * tr_xxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_xzzz[i] = 3.0 * tr_y_xzz[i] - 2.0 * tr_y_xzzzz[i] * tke_0 - 2.0 * tr_yz_xzzz[i] * tbe_0 - 6.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_yyyy[i] = -2.0 * tr_y_yyyyz[i] * tke_0 - 2.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_xxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_yyyz[i] = tr_y_yyy[i] - 2.0 * tr_y_yyyzz[i] * tke_0 - 2.0 * tr_yz_yyyz[i] * tbe_0 - 2.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_yyzz[i] = 2.0 * tr_y_yyz[i] - 2.0 * tr_y_yyzzz[i] * tke_0 - 2.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_xxy_yyz[i] * tbe_0 + 4.0 * tr_xxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_yzzz[i] = 3.0 * tr_y_yzz[i] - 2.0 * tr_y_yzzzz[i] * tke_0 - 2.0 * tr_yz_yzzz[i] * tbe_0 - 6.0 * tr_xxy_yzz[i] * tbe_0 + 4.0 * tr_xxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xy_zzzz[i] = 4.0 * tr_y_zzz[i] - 2.0 * tr_y_zzzzz[i] * tke_0 - 2.0 * tr_yz_zzzz[i] * tbe_0 - 8.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-225 components of targeted buffer : DG

    auto tr_x_0_z_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 210);

    auto tr_x_0_z_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 211);

    auto tr_x_0_z_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 212);

    auto tr_x_0_z_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 213);

    auto tr_x_0_z_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 214);

    auto tr_x_0_z_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 215);

    auto tr_x_0_z_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 216);

    auto tr_x_0_z_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 217);

    auto tr_x_0_z_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 218);

    auto tr_x_0_z_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 219);

    auto tr_x_0_z_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 220);

    auto tr_x_0_z_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 221);

    auto tr_x_0_z_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 222);

    auto tr_x_0_z_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 223);

    auto tr_x_0_z_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 224);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_x_0_z_xz_xxxx, tr_x_0_z_xz_xxxy, tr_x_0_z_xz_xxxz, tr_x_0_z_xz_xxyy, tr_x_0_z_xz_xxyz, tr_x_0_z_xz_xxzz, tr_x_0_z_xz_xyyy, tr_x_0_z_xz_xyyz, tr_x_0_z_xz_xyzz, tr_x_0_z_xz_xzzz, tr_x_0_z_xz_yyyy, tr_x_0_z_xz_yyyz, tr_x_0_z_xz_yyzz, tr_x_0_z_xz_yzzz, tr_x_0_z_xz_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_z_xxx, tr_z_xxxxz, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, tr_z_zzzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xz_xxxx[i] = tr_0_xxxx[i] - 2.0 * tr_z_xxxxz[i] * tke_0 - 2.0 * tr_zz_xxxx[i] * tbe_0 - 2.0 * tr_xx_xxxx[i] * tbe_0 + 4.0 * tr_xxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xxxy[i] = tr_0_xxxy[i] - 2.0 * tr_z_xxxyz[i] * tke_0 - 2.0 * tr_zz_xxxy[i] * tbe_0 - 2.0 * tr_xx_xxxy[i] * tbe_0 + 4.0 * tr_xxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xxxz[i] = tr_0_xxxz[i] + tr_z_xxx[i] - 2.0 * tr_z_xxxzz[i] * tke_0 - 2.0 * tr_zz_xxxz[i] * tbe_0 - 2.0 * tr_xx_xxxz[i] * tbe_0 - 2.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xxyy[i] = tr_0_xxyy[i] - 2.0 * tr_z_xxyyz[i] * tke_0 - 2.0 * tr_zz_xxyy[i] * tbe_0 - 2.0 * tr_xx_xxyy[i] * tbe_0 + 4.0 * tr_xxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xxyz[i] = tr_0_xxyz[i] + tr_z_xxy[i] - 2.0 * tr_z_xxyzz[i] * tke_0 - 2.0 * tr_zz_xxyz[i] * tbe_0 - 2.0 * tr_xx_xxyz[i] * tbe_0 - 2.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xxzz[i] = tr_0_xxzz[i] + 2.0 * tr_z_xxz[i] - 2.0 * tr_z_xxzzz[i] * tke_0 - 2.0 * tr_zz_xxzz[i] * tbe_0 - 2.0 * tr_xx_xxzz[i] * tbe_0 - 4.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xyyy[i] = tr_0_xyyy[i] - 2.0 * tr_z_xyyyz[i] * tke_0 - 2.0 * tr_zz_xyyy[i] * tbe_0 - 2.0 * tr_xx_xyyy[i] * tbe_0 + 4.0 * tr_xxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xyyz[i] = tr_0_xyyz[i] + tr_z_xyy[i] - 2.0 * tr_z_xyyzz[i] * tke_0 - 2.0 * tr_zz_xyyz[i] * tbe_0 - 2.0 * tr_xx_xyyz[i] * tbe_0 - 2.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xyzz[i] = tr_0_xyzz[i] + 2.0 * tr_z_xyz[i] - 2.0 * tr_z_xyzzz[i] * tke_0 - 2.0 * tr_zz_xyzz[i] * tbe_0 - 2.0 * tr_xx_xyzz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tbe_0 + 4.0 * tr_xxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_xzzz[i] = tr_0_xzzz[i] + 3.0 * tr_z_xzz[i] - 2.0 * tr_z_xzzzz[i] * tke_0 - 2.0 * tr_zz_xzzz[i] * tbe_0 - 2.0 * tr_xx_xzzz[i] * tbe_0 - 6.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_yyyy[i] = tr_0_yyyy[i] - 2.0 * tr_z_yyyyz[i] * tke_0 - 2.0 * tr_zz_yyyy[i] * tbe_0 - 2.0 * tr_xx_yyyy[i] * tbe_0 + 4.0 * tr_xxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_yyyz[i] = tr_0_yyyz[i] + tr_z_yyy[i] - 2.0 * tr_z_yyyzz[i] * tke_0 - 2.0 * tr_zz_yyyz[i] * tbe_0 - 2.0 * tr_xx_yyyz[i] * tbe_0 - 2.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_yyzz[i] = tr_0_yyzz[i] + 2.0 * tr_z_yyz[i] - 2.0 * tr_z_yyzzz[i] * tke_0 - 2.0 * tr_zz_yyzz[i] * tbe_0 - 2.0 * tr_xx_yyzz[i] * tbe_0 - 4.0 * tr_xxz_yyz[i] * tbe_0 + 4.0 * tr_xxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_yzzz[i] = tr_0_yzzz[i] + 3.0 * tr_z_yzz[i] - 2.0 * tr_z_yzzzz[i] * tke_0 - 2.0 * tr_zz_yzzz[i] * tbe_0 - 2.0 * tr_xx_yzzz[i] * tbe_0 - 6.0 * tr_xxz_yzz[i] * tbe_0 + 4.0 * tr_xxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xz_zzzz[i] = tr_0_zzzz[i] + 4.0 * tr_z_zzz[i] - 2.0 * tr_z_zzzzz[i] * tke_0 - 2.0 * tr_zz_zzzz[i] * tbe_0 - 2.0 * tr_xx_zzzz[i] * tbe_0 - 8.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 225-240 components of targeted buffer : DG

    auto tr_x_0_z_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 225);

    auto tr_x_0_z_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 226);

    auto tr_x_0_z_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 227);

    auto tr_x_0_z_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 228);

    auto tr_x_0_z_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 229);

    auto tr_x_0_z_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 230);

    auto tr_x_0_z_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 231);

    auto tr_x_0_z_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 232);

    auto tr_x_0_z_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 233);

    auto tr_x_0_z_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 234);

    auto tr_x_0_z_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 235);

    auto tr_x_0_z_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 236);

    auto tr_x_0_z_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 237);

    auto tr_x_0_z_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 238);

    auto tr_x_0_z_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 239);

    #pragma omp simd aligned(tr_x_0_z_yy_xxxx, tr_x_0_z_yy_xxxy, tr_x_0_z_yy_xxxz, tr_x_0_z_yy_xxyy, tr_x_0_z_yy_xxyz, tr_x_0_z_yy_xxzz, tr_x_0_z_yy_xyyy, tr_x_0_z_yy_xyyz, tr_x_0_z_yy_xyzz, tr_x_0_z_yy_xzzz, tr_x_0_z_yy_yyyy, tr_x_0_z_yy_yyyz, tr_x_0_z_yy_yyzz, tr_x_0_z_yy_yzzz, tr_x_0_z_yy_zzzz, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yy_xxxx[i] = 4.0 * tr_xyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xxxy[i] = 4.0 * tr_xyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xxxz[i] = -2.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xxyy[i] = 4.0 * tr_xyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xxyz[i] = -2.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xxzz[i] = -4.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xyyy[i] = 4.0 * tr_xyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xyyz[i] = -2.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xyzz[i] = -4.0 * tr_xyy_xyz[i] * tbe_0 + 4.0 * tr_xyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_xzzz[i] = -6.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_yyyy[i] = 4.0 * tr_xyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_yyyz[i] = -2.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_yyzz[i] = -4.0 * tr_xyy_yyz[i] * tbe_0 + 4.0 * tr_xyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_yzzz[i] = -6.0 * tr_xyy_yzz[i] * tbe_0 + 4.0 * tr_xyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yy_zzzz[i] = -8.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-255 components of targeted buffer : DG

    auto tr_x_0_z_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 240);

    auto tr_x_0_z_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 241);

    auto tr_x_0_z_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 242);

    auto tr_x_0_z_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 243);

    auto tr_x_0_z_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 244);

    auto tr_x_0_z_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 245);

    auto tr_x_0_z_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 246);

    auto tr_x_0_z_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 247);

    auto tr_x_0_z_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 248);

    auto tr_x_0_z_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 249);

    auto tr_x_0_z_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 250);

    auto tr_x_0_z_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 251);

    auto tr_x_0_z_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 252);

    auto tr_x_0_z_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 253);

    auto tr_x_0_z_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 254);

    #pragma omp simd aligned(tr_x_0_z_yz_xxxx, tr_x_0_z_yz_xxxy, tr_x_0_z_yz_xxxz, tr_x_0_z_yz_xxyy, tr_x_0_z_yz_xxyz, tr_x_0_z_yz_xxzz, tr_x_0_z_yz_xyyy, tr_x_0_z_yz_xyyz, tr_x_0_z_yz_xyzz, tr_x_0_z_yz_xzzz, tr_x_0_z_yz_yyyy, tr_x_0_z_yz_yyyz, tr_x_0_z_yz_yyzz, tr_x_0_z_yz_yzzz, tr_x_0_z_yz_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yz_xxxx[i] = -2.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xxxy[i] = -2.0 * tr_xy_xxxy[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xxxz[i] = -2.0 * tr_xy_xxxz[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xxyy[i] = -2.0 * tr_xy_xxyy[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xxyz[i] = -2.0 * tr_xy_xxyz[i] * tbe_0 - 2.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xxzz[i] = -2.0 * tr_xy_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xyyy[i] = -2.0 * tr_xy_xyyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xyyz[i] = -2.0 * tr_xy_xyyz[i] * tbe_0 - 2.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xyzz[i] = -2.0 * tr_xy_xyzz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_xzzz[i] = -2.0 * tr_xy_xzzz[i] * tbe_0 - 6.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_yyyy[i] = -2.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_yyyz[i] = -2.0 * tr_xy_yyyz[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_yyzz[i] = -2.0 * tr_xy_yyzz[i] * tbe_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_yzzz[i] = -2.0 * tr_xy_yzzz[i] * tbe_0 - 6.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yz_zzzz[i] = -2.0 * tr_xy_zzzz[i] * tbe_0 - 8.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 255-270 components of targeted buffer : DG

    auto tr_x_0_z_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 255);

    auto tr_x_0_z_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 256);

    auto tr_x_0_z_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 257);

    auto tr_x_0_z_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 258);

    auto tr_x_0_z_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 259);

    auto tr_x_0_z_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 260);

    auto tr_x_0_z_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 261);

    auto tr_x_0_z_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 262);

    auto tr_x_0_z_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 263);

    auto tr_x_0_z_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 264);

    auto tr_x_0_z_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 265);

    auto tr_x_0_z_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 266);

    auto tr_x_0_z_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 267);

    auto tr_x_0_z_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 268);

    auto tr_x_0_z_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 269);

    #pragma omp simd aligned(tr_x_0_z_zz_xxxx, tr_x_0_z_zz_xxxy, tr_x_0_z_zz_xxxz, tr_x_0_z_zz_xxyy, tr_x_0_z_zz_xxyz, tr_x_0_z_zz_xxzz, tr_x_0_z_zz_xyyy, tr_x_0_z_zz_xyyz, tr_x_0_z_zz_xyzz, tr_x_0_z_zz_xzzz, tr_x_0_z_zz_yyyy, tr_x_0_z_zz_yyyz, tr_x_0_z_zz_yyzz, tr_x_0_z_zz_yzzz, tr_x_0_z_zz_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zz_xxxx[i] = -4.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xxxy[i] = -4.0 * tr_xz_xxxy[i] * tbe_0 + 4.0 * tr_xzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xxxz[i] = -4.0 * tr_xz_xxxz[i] * tbe_0 - 2.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xxyy[i] = -4.0 * tr_xz_xxyy[i] * tbe_0 + 4.0 * tr_xzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xxyz[i] = -4.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xxzz[i] = -4.0 * tr_xz_xxzz[i] * tbe_0 - 4.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xyyy[i] = -4.0 * tr_xz_xyyy[i] * tbe_0 + 4.0 * tr_xzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xyyz[i] = -4.0 * tr_xz_xyyz[i] * tbe_0 - 2.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xyzz[i] = -4.0 * tr_xz_xyzz[i] * tbe_0 - 4.0 * tr_xzz_xyz[i] * tbe_0 + 4.0 * tr_xzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_xzzz[i] = -4.0 * tr_xz_xzzz[i] * tbe_0 - 6.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_yyyy[i] = -4.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_yyyz[i] = -4.0 * tr_xz_yyyz[i] * tbe_0 - 2.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_yyzz[i] = -4.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xzz_yyz[i] * tbe_0 + 4.0 * tr_xzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_yzzz[i] = -4.0 * tr_xz_yzzz[i] * tbe_0 - 6.0 * tr_xzz_yzz[i] * tbe_0 + 4.0 * tr_xzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zz_zzzz[i] = -4.0 * tr_xz_zzzz[i] * tbe_0 - 8.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-285 components of targeted buffer : DG

    auto tr_y_0_x_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 270);

    auto tr_y_0_x_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 271);

    auto tr_y_0_x_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 272);

    auto tr_y_0_x_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 273);

    auto tr_y_0_x_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 274);

    auto tr_y_0_x_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 275);

    auto tr_y_0_x_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 276);

    auto tr_y_0_x_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 277);

    auto tr_y_0_x_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 278);

    auto tr_y_0_x_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 279);

    auto tr_y_0_x_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 280);

    auto tr_y_0_x_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 281);

    auto tr_y_0_x_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 282);

    auto tr_y_0_x_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 283);

    auto tr_y_0_x_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 284);

    #pragma omp simd aligned(tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_y_0_x_xx_xxxx, tr_y_0_x_xx_xxxy, tr_y_0_x_xx_xxxz, tr_y_0_x_xx_xxyy, tr_y_0_x_xx_xxyz, tr_y_0_x_xx_xxzz, tr_y_0_x_xx_xyyy, tr_y_0_x_xx_xyyz, tr_y_0_x_xx_xyzz, tr_y_0_x_xx_xzzz, tr_y_0_x_xx_yyyy, tr_y_0_x_xx_yyyz, tr_y_0_x_xx_yyzz, tr_y_0_x_xx_yzzz, tr_y_0_x_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xx_xxxx[i] = -4.0 * tr_xy_xxxx[i] * tbe_0 - 8.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xxxy[i] = -4.0 * tr_xy_xxxy[i] * tbe_0 - 6.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xxxz[i] = -4.0 * tr_xy_xxxz[i] * tbe_0 - 6.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xxyy[i] = -4.0 * tr_xy_xxyy[i] * tbe_0 - 4.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xxyz[i] = -4.0 * tr_xy_xxyz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tbe_0 + 4.0 * tr_xxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xxzz[i] = -4.0 * tr_xy_xxzz[i] * tbe_0 - 4.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xyyy[i] = -4.0 * tr_xy_xyyy[i] * tbe_0 - 2.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xyyz[i] = -4.0 * tr_xy_xyyz[i] * tbe_0 - 2.0 * tr_xxy_yyz[i] * tbe_0 + 4.0 * tr_xxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xyzz[i] = -4.0 * tr_xy_xyzz[i] * tbe_0 - 2.0 * tr_xxy_yzz[i] * tbe_0 + 4.0 * tr_xxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_xzzz[i] = -4.0 * tr_xy_xzzz[i] * tbe_0 - 2.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_yyyy[i] = -4.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_yyyz[i] = -4.0 * tr_xy_yyyz[i] * tbe_0 + 4.0 * tr_xxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_yyzz[i] = -4.0 * tr_xy_yyzz[i] * tbe_0 + 4.0 * tr_xxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_yzzz[i] = -4.0 * tr_xy_yzzz[i] * tbe_0 + 4.0 * tr_xxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xx_zzzz[i] = -4.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 285-300 components of targeted buffer : DG

    auto tr_y_0_x_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 285);

    auto tr_y_0_x_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 286);

    auto tr_y_0_x_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 287);

    auto tr_y_0_x_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 288);

    auto tr_y_0_x_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 289);

    auto tr_y_0_x_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 290);

    auto tr_y_0_x_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 291);

    auto tr_y_0_x_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 292);

    auto tr_y_0_x_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 293);

    auto tr_y_0_x_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 294);

    auto tr_y_0_x_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 295);

    auto tr_y_0_x_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 296);

    auto tr_y_0_x_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 297);

    auto tr_y_0_x_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 298);

    auto tr_y_0_x_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 299);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_x_xxx, tr_x_xxxxx, tr_x_xxxxy, tr_x_xxxxz, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_y_0_x_xy_xxxx, tr_y_0_x_xy_xxxy, tr_y_0_x_xy_xxxz, tr_y_0_x_xy_xxyy, tr_y_0_x_xy_xxyz, tr_y_0_x_xy_xxzz, tr_y_0_x_xy_xyyy, tr_y_0_x_xy_xyyz, tr_y_0_x_xy_xyzz, tr_y_0_x_xy_xzzz, tr_y_0_x_xy_yyyy, tr_y_0_x_xy_yyyz, tr_y_0_x_xy_yyzz, tr_y_0_x_xy_yzzz, tr_y_0_x_xy_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xy_xxxx[i] = tr_0_xxxx[i] - 2.0 * tr_yy_xxxx[i] * tbe_0 + 4.0 * tr_x_xxx[i] - 2.0 * tr_x_xxxxx[i] * tke_0 - 8.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xxxy[i] = tr_0_xxxy[i] - 2.0 * tr_yy_xxxy[i] * tbe_0 + 3.0 * tr_x_xxy[i] - 2.0 * tr_x_xxxxy[i] * tke_0 - 6.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxxy[i] * tbe_0 + 4.0 * tr_xxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xxxz[i] = tr_0_xxxz[i] - 2.0 * tr_yy_xxxz[i] * tbe_0 + 3.0 * tr_x_xxz[i] - 2.0 * tr_x_xxxxz[i] * tke_0 - 6.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxxz[i] * tbe_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xxyy[i] = tr_0_xxyy[i] - 2.0 * tr_yy_xxyy[i] * tbe_0 + 2.0 * tr_x_xyy[i] - 2.0 * tr_x_xxxyy[i] * tke_0 - 4.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyy_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxyy[i] * tbe_0 + 4.0 * tr_xxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xxyz[i] = tr_0_xxyz[i] - 2.0 * tr_yy_xxyz[i] * tbe_0 + 2.0 * tr_x_xyz[i] - 2.0 * tr_x_xxxyz[i] * tke_0 - 4.0 * tr_xyy_xyz[i] * tbe_0 + 4.0 * tr_xyy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxyz[i] * tbe_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xxzz[i] = tr_0_xxzz[i] - 2.0 * tr_yy_xxzz[i] * tbe_0 + 2.0 * tr_x_xzz[i] - 2.0 * tr_x_xxxzz[i] * tke_0 - 4.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxzz[i] * tbe_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xyyy[i] = tr_0_xyyy[i] - 2.0 * tr_yy_xyyy[i] * tbe_0 + tr_x_yyy[i] - 2.0 * tr_x_xxyyy[i] * tke_0 - 2.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyy_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xyyz[i] = tr_0_xyyz[i] - 2.0 * tr_yy_xyyz[i] * tbe_0 + tr_x_yyz[i] - 2.0 * tr_x_xxyyz[i] * tke_0 - 2.0 * tr_xyy_yyz[i] * tbe_0 + 4.0 * tr_xyy_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyyz[i] * tbe_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xyzz[i] = tr_0_xyzz[i] - 2.0 * tr_yy_xyzz[i] * tbe_0 + tr_x_yzz[i] - 2.0 * tr_x_xxyzz[i] * tke_0 - 2.0 * tr_xyy_yzz[i] * tbe_0 + 4.0 * tr_xyy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyzz[i] * tbe_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_xzzz[i] = tr_0_xzzz[i] - 2.0 * tr_yy_xzzz[i] * tbe_0 + tr_x_zzz[i] - 2.0 * tr_x_xxzzz[i] * tke_0 - 2.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyy_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xzzz[i] * tbe_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_yyyy[i] = tr_0_yyyy[i] - 2.0 * tr_yy_yyyy[i] * tbe_0 - 2.0 * tr_x_xyyyy[i] * tke_0 + 4.0 * tr_xyy_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_yyyz[i] = tr_0_yyyz[i] - 2.0 * tr_yy_yyyz[i] * tbe_0 - 2.0 * tr_x_xyyyz[i] * tke_0 + 4.0 * tr_xyy_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyyz[i] * tbe_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_yyzz[i] = tr_0_yyzz[i] - 2.0 * tr_yy_yyzz[i] * tbe_0 - 2.0 * tr_x_xyyzz[i] * tke_0 + 4.0 * tr_xyy_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyzz[i] * tbe_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_yzzz[i] = tr_0_yzzz[i] - 2.0 * tr_yy_yzzz[i] * tbe_0 - 2.0 * tr_x_xyzzz[i] * tke_0 + 4.0 * tr_xyy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yzzz[i] * tbe_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xy_zzzz[i] = tr_0_zzzz[i] - 2.0 * tr_yy_zzzz[i] * tbe_0 - 2.0 * tr_x_xzzzz[i] * tke_0 + 4.0 * tr_xyy_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zzzz[i] * tbe_0 + 4.0 * tr_xxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-315 components of targeted buffer : DG

    auto tr_y_0_x_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 300);

    auto tr_y_0_x_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 301);

    auto tr_y_0_x_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 302);

    auto tr_y_0_x_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 303);

    auto tr_y_0_x_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 304);

    auto tr_y_0_x_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 305);

    auto tr_y_0_x_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 306);

    auto tr_y_0_x_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 307);

    auto tr_y_0_x_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 308);

    auto tr_y_0_x_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 309);

    auto tr_y_0_x_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 310);

    auto tr_y_0_x_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 311);

    auto tr_y_0_x_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 312);

    auto tr_y_0_x_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 313);

    auto tr_y_0_x_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 314);

    #pragma omp simd aligned(tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_y_0_x_xz_xxxx, tr_y_0_x_xz_xxxy, tr_y_0_x_xz_xxxz, tr_y_0_x_xz_xxyy, tr_y_0_x_xz_xxyz, tr_y_0_x_xz_xxzz, tr_y_0_x_xz_xyyy, tr_y_0_x_xz_xyyz, tr_y_0_x_xz_xyzz, tr_y_0_x_xz_xzzz, tr_y_0_x_xz_yyyy, tr_y_0_x_xz_yyyz, tr_y_0_x_xz_yyzz, tr_y_0_x_xz_yzzz, tr_y_0_x_xz_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xz_xxxx[i] = -2.0 * tr_yz_xxxx[i] * tbe_0 - 8.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xxxy[i] = -2.0 * tr_yz_xxxy[i] * tbe_0 - 6.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xxxz[i] = -2.0 * tr_yz_xxxz[i] * tbe_0 - 6.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xxyy[i] = -2.0 * tr_yz_xxyy[i] * tbe_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xxyz[i] = -2.0 * tr_yz_xxyz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xxzz[i] = -2.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xyyy[i] = -2.0 * tr_yz_xyyy[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xyyz[i] = -2.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xyzz[i] = -2.0 * tr_yz_xyzz[i] * tbe_0 - 2.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_xzzz[i] = -2.0 * tr_yz_xzzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_yyyy[i] = -2.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_yyyz[i] = -2.0 * tr_yz_yyyz[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_yyzz[i] = -2.0 * tr_yz_yyzz[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_yzzz[i] = -2.0 * tr_yz_yzzz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xz_zzzz[i] = -2.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_xyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 315-330 components of targeted buffer : DG

    auto tr_y_0_x_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 315);

    auto tr_y_0_x_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 316);

    auto tr_y_0_x_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 317);

    auto tr_y_0_x_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 318);

    auto tr_y_0_x_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 319);

    auto tr_y_0_x_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 320);

    auto tr_y_0_x_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 321);

    auto tr_y_0_x_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 322);

    auto tr_y_0_x_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 323);

    auto tr_y_0_x_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 324);

    auto tr_y_0_x_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 325);

    auto tr_y_0_x_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 326);

    auto tr_y_0_x_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 327);

    auto tr_y_0_x_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 328);

    auto tr_y_0_x_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 329);

    #pragma omp simd aligned(tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, tr_y_0_x_yy_xxxx, tr_y_0_x_yy_xxxy, tr_y_0_x_yy_xxxz, tr_y_0_x_yy_xxyy, tr_y_0_x_yy_xxyz, tr_y_0_x_yy_xxzz, tr_y_0_x_yy_xyyy, tr_y_0_x_yy_xyyz, tr_y_0_x_yy_xyzz, tr_y_0_x_yy_xzzz, tr_y_0_x_yy_yyyy, tr_y_0_x_yy_yyyz, tr_y_0_x_yy_yyzz, tr_y_0_x_yy_yzzz, tr_y_0_x_yy_zzzz, tr_y_xxx, tr_y_xxxxx, tr_y_xxxxy, tr_y_xxxxz, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yyy_xxx, tr_yyy_xxxxx, tr_yyy_xxxxy, tr_yyy_xxxxz, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yy_xxxx[i] = 8.0 * tr_y_xxx[i] - 4.0 * tr_y_xxxxx[i] * tke_0 - 8.0 * tr_yyy_xxx[i] * tbe_0 + 4.0 * tr_yyy_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xxxy[i] = 6.0 * tr_y_xxy[i] - 4.0 * tr_y_xxxxy[i] * tke_0 - 6.0 * tr_yyy_xxy[i] * tbe_0 + 4.0 * tr_yyy_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxxy[i] * tbe_0 + 4.0 * tr_xyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xxxz[i] = 6.0 * tr_y_xxz[i] - 4.0 * tr_y_xxxxz[i] * tke_0 - 6.0 * tr_yyy_xxz[i] * tbe_0 + 4.0 * tr_yyy_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxxz[i] * tbe_0 + 4.0 * tr_xyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xxyy[i] = 4.0 * tr_y_xyy[i] - 4.0 * tr_y_xxxyy[i] * tke_0 - 4.0 * tr_yyy_xyy[i] * tbe_0 + 4.0 * tr_yyy_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxyy[i] * tbe_0 + 4.0 * tr_xyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xxyz[i] = 4.0 * tr_y_xyz[i] - 4.0 * tr_y_xxxyz[i] * tke_0 - 4.0 * tr_yyy_xyz[i] * tbe_0 + 4.0 * tr_yyy_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxyz[i] * tbe_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xxzz[i] = 4.0 * tr_y_xzz[i] - 4.0 * tr_y_xxxzz[i] * tke_0 - 4.0 * tr_yyy_xzz[i] * tbe_0 + 4.0 * tr_yyy_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxzz[i] * tbe_0 + 4.0 * tr_xyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xyyy[i] = 2.0 * tr_y_yyy[i] - 4.0 * tr_y_xxyyy[i] * tke_0 - 2.0 * tr_yyy_yyy[i] * tbe_0 + 4.0 * tr_yyy_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xyyz[i] = 2.0 * tr_y_yyz[i] - 4.0 * tr_y_xxyyz[i] * tke_0 - 2.0 * tr_yyy_yyz[i] * tbe_0 + 4.0 * tr_yyy_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyyz[i] * tbe_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xyzz[i] = 2.0 * tr_y_yzz[i] - 4.0 * tr_y_xxyzz[i] * tke_0 - 2.0 * tr_yyy_yzz[i] * tbe_0 + 4.0 * tr_yyy_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyzz[i] * tbe_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_xzzz[i] = 2.0 * tr_y_zzz[i] - 4.0 * tr_y_xxzzz[i] * tke_0 - 2.0 * tr_yyy_zzz[i] * tbe_0 + 4.0 * tr_yyy_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xzzz[i] * tbe_0 + 4.0 * tr_xyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_yyyy[i] = -4.0 * tr_y_xyyyy[i] * tke_0 + 4.0 * tr_yyy_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_yyyz[i] = -4.0 * tr_y_xyyyz[i] * tke_0 + 4.0 * tr_yyy_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyyz[i] * tbe_0 + 4.0 * tr_xyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_yyzz[i] = -4.0 * tr_y_xyyzz[i] * tke_0 + 4.0 * tr_yyy_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyzz[i] * tbe_0 + 4.0 * tr_xyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_yzzz[i] = -4.0 * tr_y_xyzzz[i] * tke_0 + 4.0 * tr_yyy_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yzzz[i] * tbe_0 + 4.0 * tr_xyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yy_zzzz[i] = -4.0 * tr_y_xzzzz[i] * tke_0 + 4.0 * tr_yyy_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-345 components of targeted buffer : DG

    auto tr_y_0_x_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 330);

    auto tr_y_0_x_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 331);

    auto tr_y_0_x_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 332);

    auto tr_y_0_x_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 333);

    auto tr_y_0_x_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 334);

    auto tr_y_0_x_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 335);

    auto tr_y_0_x_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 336);

    auto tr_y_0_x_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 337);

    auto tr_y_0_x_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 338);

    auto tr_y_0_x_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 339);

    auto tr_y_0_x_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 340);

    auto tr_y_0_x_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 341);

    auto tr_y_0_x_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 342);

    auto tr_y_0_x_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 343);

    auto tr_y_0_x_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 344);

    #pragma omp simd aligned(tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_y_0_x_yz_xxxx, tr_y_0_x_yz_xxxy, tr_y_0_x_yz_xxxz, tr_y_0_x_yz_xxyy, tr_y_0_x_yz_xxyz, tr_y_0_x_yz_xxzz, tr_y_0_x_yz_xyyy, tr_y_0_x_yz_xyyz, tr_y_0_x_yz_xyzz, tr_y_0_x_yz_xzzz, tr_y_0_x_yz_yyyy, tr_y_0_x_yz_yyyz, tr_y_0_x_yz_yyzz, tr_y_0_x_yz_yzzz, tr_y_0_x_yz_zzzz, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_z_xxx, tr_z_xxxxx, tr_z_xxxxy, tr_z_xxxxz, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yz_xxxx[i] = 4.0 * tr_z_xxx[i] - 2.0 * tr_z_xxxxx[i] * tke_0 - 8.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xxxy[i] = 3.0 * tr_z_xxy[i] - 2.0 * tr_z_xxxxy[i] * tke_0 - 6.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xxxz[i] = 3.0 * tr_z_xxz[i] - 2.0 * tr_z_xxxxz[i] * tke_0 - 6.0 * tr_yyz_xxz[i] * tbe_0 + 4.0 * tr_yyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xxyy[i] = 2.0 * tr_z_xyy[i] - 2.0 * tr_z_xxxyy[i] * tke_0 - 4.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xxyz[i] = 2.0 * tr_z_xyz[i] - 2.0 * tr_z_xxxyz[i] * tke_0 - 4.0 * tr_yyz_xyz[i] * tbe_0 + 4.0 * tr_yyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxyz[i] * tbe_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xxzz[i] = 2.0 * tr_z_xzz[i] - 2.0 * tr_z_xxxzz[i] * tke_0 - 4.0 * tr_yyz_xzz[i] * tbe_0 + 4.0 * tr_yyz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxzz[i] * tbe_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xyyy[i] = tr_z_yyy[i] - 2.0 * tr_z_xxyyy[i] * tke_0 - 2.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xyyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xyyz[i] = tr_z_yyz[i] - 2.0 * tr_z_xxyyz[i] * tke_0 - 2.0 * tr_yyz_yyz[i] * tbe_0 + 4.0 * tr_yyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xyyz[i] * tbe_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xyzz[i] = tr_z_yzz[i] - 2.0 * tr_z_xxyzz[i] * tke_0 - 2.0 * tr_yyz_yzz[i] * tbe_0 + 4.0 * tr_yyz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xyzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_xzzz[i] = tr_z_zzz[i] - 2.0 * tr_z_xxzzz[i] * tke_0 - 2.0 * tr_yyz_zzz[i] * tbe_0 + 4.0 * tr_yyz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xzzz[i] * tbe_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_yyyy[i] = -2.0 * tr_z_xyyyy[i] * tke_0 + 4.0 * tr_yyz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_yyyz[i] = -2.0 * tr_z_xyyyz[i] * tke_0 + 4.0 * tr_yyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yyyz[i] * tbe_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_yyzz[i] = -2.0 * tr_z_xyyzz[i] * tke_0 + 4.0 * tr_yyz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yyzz[i] * tbe_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_yzzz[i] = -2.0 * tr_z_xyzzz[i] * tke_0 + 4.0 * tr_yyz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yzzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yz_zzzz[i] = -2.0 * tr_z_xzzzz[i] * tke_0 + 4.0 * tr_yyz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 345-360 components of targeted buffer : DG

    auto tr_y_0_x_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 345);

    auto tr_y_0_x_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 346);

    auto tr_y_0_x_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 347);

    auto tr_y_0_x_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 348);

    auto tr_y_0_x_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 349);

    auto tr_y_0_x_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 350);

    auto tr_y_0_x_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 351);

    auto tr_y_0_x_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 352);

    auto tr_y_0_x_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 353);

    auto tr_y_0_x_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 354);

    auto tr_y_0_x_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 355);

    auto tr_y_0_x_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 356);

    auto tr_y_0_x_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 357);

    auto tr_y_0_x_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 358);

    auto tr_y_0_x_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 359);

    #pragma omp simd aligned(tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_y_0_x_zz_xxxx, tr_y_0_x_zz_xxxy, tr_y_0_x_zz_xxxz, tr_y_0_x_zz_xxyy, tr_y_0_x_zz_xxyz, tr_y_0_x_zz_xxzz, tr_y_0_x_zz_xyyy, tr_y_0_x_zz_xyyz, tr_y_0_x_zz_xyzz, tr_y_0_x_zz_xzzz, tr_y_0_x_zz_yyyy, tr_y_0_x_zz_yyyz, tr_y_0_x_zz_yyzz, tr_y_0_x_zz_yzzz, tr_y_0_x_zz_zzzz, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zz_xxxx[i] = -8.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xxxy[i] = -6.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xxxz[i] = -6.0 * tr_yzz_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xxyy[i] = -4.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xxyz[i] = -4.0 * tr_yzz_xyz[i] * tbe_0 + 4.0 * tr_yzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xxzz[i] = -4.0 * tr_yzz_xzz[i] * tbe_0 + 4.0 * tr_yzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xyyy[i] = -2.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xyyz[i] = -2.0 * tr_yzz_yyz[i] * tbe_0 + 4.0 * tr_yzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xyzz[i] = -2.0 * tr_yzz_yzz[i] * tbe_0 + 4.0 * tr_yzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_xzzz[i] = -2.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_yzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_yyyy[i] = 4.0 * tr_yzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_yyyz[i] = 4.0 * tr_yzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_yyzz[i] = 4.0 * tr_yzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_yzzz[i] = 4.0 * tr_yzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zz_zzzz[i] = 4.0 * tr_yzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-375 components of targeted buffer : DG

    auto tr_y_0_y_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 360);

    auto tr_y_0_y_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 361);

    auto tr_y_0_y_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 362);

    auto tr_y_0_y_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 363);

    auto tr_y_0_y_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 364);

    auto tr_y_0_y_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 365);

    auto tr_y_0_y_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 366);

    auto tr_y_0_y_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 367);

    auto tr_y_0_y_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 368);

    auto tr_y_0_y_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 369);

    auto tr_y_0_y_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 370);

    auto tr_y_0_y_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 371);

    auto tr_y_0_y_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 372);

    auto tr_y_0_y_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 373);

    auto tr_y_0_y_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 374);

    #pragma omp simd aligned(tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_y_0_y_xx_xxxx, tr_y_0_y_xx_xxxy, tr_y_0_y_xx_xxxz, tr_y_0_y_xx_xxyy, tr_y_0_y_xx_xxyz, tr_y_0_y_xx_xxzz, tr_y_0_y_xx_xyyy, tr_y_0_y_xx_xyyz, tr_y_0_y_xx_xyzz, tr_y_0_y_xx_xzzz, tr_y_0_y_xx_yyyy, tr_y_0_y_xx_yyyz, tr_y_0_y_xx_yyzz, tr_y_0_y_xx_yzzz, tr_y_0_y_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xx_xxxx[i] = -2.0 * tr_xx_xxxx[i] * tbe_0 + 4.0 * tr_xxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xxxy[i] = -2.0 * tr_xx_xxxy[i] * tbe_0 - 2.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xxxz[i] = -2.0 * tr_xx_xxxz[i] * tbe_0 + 4.0 * tr_xxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xxyy[i] = -2.0 * tr_xx_xxyy[i] * tbe_0 - 4.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xxyz[i] = -2.0 * tr_xx_xxyz[i] * tbe_0 - 2.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xxzz[i] = -2.0 * tr_xx_xxzz[i] * tbe_0 + 4.0 * tr_xxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xyyy[i] = -2.0 * tr_xx_xyyy[i] * tbe_0 - 6.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xyyz[i] = -2.0 * tr_xx_xyyz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tbe_0 + 4.0 * tr_xxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xyzz[i] = -2.0 * tr_xx_xyzz[i] * tbe_0 - 2.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_xzzz[i] = -2.0 * tr_xx_xzzz[i] * tbe_0 + 4.0 * tr_xxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_yyyy[i] = -2.0 * tr_xx_yyyy[i] * tbe_0 - 8.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_yyyz[i] = -2.0 * tr_xx_yyyz[i] * tbe_0 - 6.0 * tr_xxy_yyz[i] * tbe_0 + 4.0 * tr_xxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_yyzz[i] = -2.0 * tr_xx_yyzz[i] * tbe_0 - 4.0 * tr_xxy_yzz[i] * tbe_0 + 4.0 * tr_xxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_yzzz[i] = -2.0 * tr_xx_yzzz[i] * tbe_0 - 2.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xx_zzzz[i] = -2.0 * tr_xx_zzzz[i] * tbe_0 + 4.0 * tr_xxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 375-390 components of targeted buffer : DG

    auto tr_y_0_y_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 375);

    auto tr_y_0_y_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 376);

    auto tr_y_0_y_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 377);

    auto tr_y_0_y_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 378);

    auto tr_y_0_y_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 379);

    auto tr_y_0_y_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 380);

    auto tr_y_0_y_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 381);

    auto tr_y_0_y_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 382);

    auto tr_y_0_y_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 383);

    auto tr_y_0_y_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 384);

    auto tr_y_0_y_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 385);

    auto tr_y_0_y_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 386);

    auto tr_y_0_y_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 387);

    auto tr_y_0_y_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 388);

    auto tr_y_0_y_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 389);

    #pragma omp simd aligned(tr_x_xxx, tr_x_xxxxy, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_yyy, tr_x_yyyyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, tr_y_0_y_xy_xxxx, tr_y_0_y_xy_xxxy, tr_y_0_y_xy_xxxz, tr_y_0_y_xy_xxyy, tr_y_0_y_xy_xxyz, tr_y_0_y_xy_xxzz, tr_y_0_y_xy_xyyy, tr_y_0_y_xy_xyyz, tr_y_0_y_xy_xyzz, tr_y_0_y_xy_xzzz, tr_y_0_y_xy_yyyy, tr_y_0_y_xy_yyyz, tr_y_0_y_xy_yyzz, tr_y_0_y_xy_yzzz, tr_y_0_y_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xy_xxxx[i] = -2.0 * tr_x_xxxxy[i] * tke_0 - 6.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xxxy[i] = tr_x_xxx[i] - 2.0 * tr_x_xxxyy[i] * tke_0 - 6.0 * tr_xy_xxxy[i] * tbe_0 - 2.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xxxz[i] = -2.0 * tr_x_xxxyz[i] * tke_0 - 6.0 * tr_xy_xxxz[i] * tbe_0 + 4.0 * tr_xyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xxyy[i] = 2.0 * tr_x_xxy[i] - 2.0 * tr_x_xxyyy[i] * tke_0 - 6.0 * tr_xy_xxyy[i] * tbe_0 - 4.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xxyz[i] = tr_x_xxz[i] - 2.0 * tr_x_xxyyz[i] * tke_0 - 6.0 * tr_xy_xxyz[i] * tbe_0 - 2.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xxzz[i] = -2.0 * tr_x_xxyzz[i] * tke_0 - 6.0 * tr_xy_xxzz[i] * tbe_0 + 4.0 * tr_xyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xyyy[i] = 3.0 * tr_x_xyy[i] - 2.0 * tr_x_xyyyy[i] * tke_0 - 6.0 * tr_xy_xyyy[i] * tbe_0 - 6.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xyyz[i] = 2.0 * tr_x_xyz[i] - 2.0 * tr_x_xyyyz[i] * tke_0 - 6.0 * tr_xy_xyyz[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tbe_0 + 4.0 * tr_xyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xyzz[i] = tr_x_xzz[i] - 2.0 * tr_x_xyyzz[i] * tke_0 - 6.0 * tr_xy_xyzz[i] * tbe_0 - 2.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_xzzz[i] = -2.0 * tr_x_xyzzz[i] * tke_0 - 6.0 * tr_xy_xzzz[i] * tbe_0 + 4.0 * tr_xyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_yyyy[i] = 4.0 * tr_x_yyy[i] - 2.0 * tr_x_yyyyy[i] * tke_0 - 6.0 * tr_xy_yyyy[i] * tbe_0 - 8.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_yyyz[i] = 3.0 * tr_x_yyz[i] - 2.0 * tr_x_yyyyz[i] * tke_0 - 6.0 * tr_xy_yyyz[i] * tbe_0 - 6.0 * tr_xyy_yyz[i] * tbe_0 + 4.0 * tr_xyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_yyzz[i] = 2.0 * tr_x_yzz[i] - 2.0 * tr_x_yyyzz[i] * tke_0 - 6.0 * tr_xy_yyzz[i] * tbe_0 - 4.0 * tr_xyy_yzz[i] * tbe_0 + 4.0 * tr_xyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_yzzz[i] = tr_x_zzz[i] - 2.0 * tr_x_yyzzz[i] * tke_0 - 6.0 * tr_xy_yzzz[i] * tbe_0 - 2.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xy_zzzz[i] = -2.0 * tr_x_yzzzz[i] * tke_0 - 6.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-405 components of targeted buffer : DG

    auto tr_y_0_y_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 390);

    auto tr_y_0_y_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 391);

    auto tr_y_0_y_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 392);

    auto tr_y_0_y_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 393);

    auto tr_y_0_y_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 394);

    auto tr_y_0_y_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 395);

    auto tr_y_0_y_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 396);

    auto tr_y_0_y_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 397);

    auto tr_y_0_y_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 398);

    auto tr_y_0_y_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 399);

    auto tr_y_0_y_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 400);

    auto tr_y_0_y_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 401);

    auto tr_y_0_y_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 402);

    auto tr_y_0_y_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 403);

    auto tr_y_0_y_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 404);

    #pragma omp simd aligned(tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_y_0_y_xz_xxxx, tr_y_0_y_xz_xxxy, tr_y_0_y_xz_xxxz, tr_y_0_y_xz_xxyy, tr_y_0_y_xz_xxyz, tr_y_0_y_xz_xxzz, tr_y_0_y_xz_xyyy, tr_y_0_y_xz_xyyz, tr_y_0_y_xz_xyzz, tr_y_0_y_xz_xzzz, tr_y_0_y_xz_yyyy, tr_y_0_y_xz_yyyz, tr_y_0_y_xz_yyzz, tr_y_0_y_xz_yzzz, tr_y_0_y_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xz_xxxx[i] = -2.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xxxy[i] = -2.0 * tr_xz_xxxy[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xxxz[i] = -2.0 * tr_xz_xxxz[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xxyy[i] = -2.0 * tr_xz_xxyy[i] * tbe_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xxyz[i] = -2.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xxzz[i] = -2.0 * tr_xz_xxzz[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xyyy[i] = -2.0 * tr_xz_xyyy[i] * tbe_0 - 6.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xyyz[i] = -2.0 * tr_xz_xyyz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xyzz[i] = -2.0 * tr_xz_xyzz[i] * tbe_0 - 2.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_xzzz[i] = -2.0 * tr_xz_xzzz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_yyyy[i] = -2.0 * tr_xz_yyyy[i] * tbe_0 - 8.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_yyyz[i] = -2.0 * tr_xz_yyyz[i] * tbe_0 - 6.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_yyzz[i] = -2.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_yzzz[i] = -2.0 * tr_xz_yzzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xz_zzzz[i] = -2.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 405-420 components of targeted buffer : DG

    auto tr_y_0_y_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 405);

    auto tr_y_0_y_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 406);

    auto tr_y_0_y_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 407);

    auto tr_y_0_y_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 408);

    auto tr_y_0_y_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 409);

    auto tr_y_0_y_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 410);

    auto tr_y_0_y_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 411);

    auto tr_y_0_y_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 412);

    auto tr_y_0_y_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 413);

    auto tr_y_0_y_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 414);

    auto tr_y_0_y_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 415);

    auto tr_y_0_y_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 416);

    auto tr_y_0_y_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 417);

    auto tr_y_0_y_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 418);

    auto tr_y_0_y_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 419);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_y_0_y_yy_xxxx, tr_y_0_y_yy_xxxy, tr_y_0_y_yy_xxxz, tr_y_0_y_yy_xxyy, tr_y_0_y_yy_xxyz, tr_y_0_y_yy_xxzz, tr_y_0_y_yy_xyyy, tr_y_0_y_yy_xyyz, tr_y_0_y_yy_xyzz, tr_y_0_y_yy_xzzz, tr_y_0_y_yy_yyyy, tr_y_0_y_yy_yyyz, tr_y_0_y_yy_yyzz, tr_y_0_y_yy_yzzz, tr_y_0_y_yy_zzzz, tr_y_xxx, tr_y_xxxxy, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_yyy, tr_y_yyyyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyy_xxx, tr_yyy_xxxxy, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyyyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xzzz, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yzzz, tr_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yy_xxxx[i] = 2.0 * tr_0_xxxx[i] - 4.0 * tr_y_xxxxy[i] * tke_0 - 10.0 * tr_yy_xxxx[i] * tbe_0 + 4.0 * tr_yyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xxxy[i] = 2.0 * tr_0_xxxy[i] + 2.0 * tr_y_xxx[i] - 4.0 * tr_y_xxxyy[i] * tke_0 - 10.0 * tr_yy_xxxy[i] * tbe_0 - 2.0 * tr_yyy_xxx[i] * tbe_0 + 4.0 * tr_yyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xxxz[i] = 2.0 * tr_0_xxxz[i] - 4.0 * tr_y_xxxyz[i] * tke_0 - 10.0 * tr_yy_xxxz[i] * tbe_0 + 4.0 * tr_yyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xxyy[i] = 2.0 * tr_0_xxyy[i] + 4.0 * tr_y_xxy[i] - 4.0 * tr_y_xxyyy[i] * tke_0 - 10.0 * tr_yy_xxyy[i] * tbe_0 - 4.0 * tr_yyy_xxy[i] * tbe_0 + 4.0 * tr_yyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xxyz[i] = 2.0 * tr_0_xxyz[i] + 2.0 * tr_y_xxz[i] - 4.0 * tr_y_xxyyz[i] * tke_0 - 10.0 * tr_yy_xxyz[i] * tbe_0 - 2.0 * tr_yyy_xxz[i] * tbe_0 + 4.0 * tr_yyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xxzz[i] = 2.0 * tr_0_xxzz[i] - 4.0 * tr_y_xxyzz[i] * tke_0 - 10.0 * tr_yy_xxzz[i] * tbe_0 + 4.0 * tr_yyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xyyy[i] = 2.0 * tr_0_xyyy[i] + 6.0 * tr_y_xyy[i] - 4.0 * tr_y_xyyyy[i] * tke_0 - 10.0 * tr_yy_xyyy[i] * tbe_0 - 6.0 * tr_yyy_xyy[i] * tbe_0 + 4.0 * tr_yyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xyyz[i] = 2.0 * tr_0_xyyz[i] + 4.0 * tr_y_xyz[i] - 4.0 * tr_y_xyyyz[i] * tke_0 - 10.0 * tr_yy_xyyz[i] * tbe_0 - 4.0 * tr_yyy_xyz[i] * tbe_0 + 4.0 * tr_yyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xyzz[i] = 2.0 * tr_0_xyzz[i] + 2.0 * tr_y_xzz[i] - 4.0 * tr_y_xyyzz[i] * tke_0 - 10.0 * tr_yy_xyzz[i] * tbe_0 - 2.0 * tr_yyy_xzz[i] * tbe_0 + 4.0 * tr_yyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_xzzz[i] = 2.0 * tr_0_xzzz[i] - 4.0 * tr_y_xyzzz[i] * tke_0 - 10.0 * tr_yy_xzzz[i] * tbe_0 + 4.0 * tr_yyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_yyyy[i] = 2.0 * tr_0_yyyy[i] + 8.0 * tr_y_yyy[i] - 4.0 * tr_y_yyyyy[i] * tke_0 - 10.0 * tr_yy_yyyy[i] * tbe_0 - 8.0 * tr_yyy_yyy[i] * tbe_0 + 4.0 * tr_yyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_yyyz[i] = 2.0 * tr_0_yyyz[i] + 6.0 * tr_y_yyz[i] - 4.0 * tr_y_yyyyz[i] * tke_0 - 10.0 * tr_yy_yyyz[i] * tbe_0 - 6.0 * tr_yyy_yyz[i] * tbe_0 + 4.0 * tr_yyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_yyzz[i] = 2.0 * tr_0_yyzz[i] + 4.0 * tr_y_yzz[i] - 4.0 * tr_y_yyyzz[i] * tke_0 - 10.0 * tr_yy_yyzz[i] * tbe_0 - 4.0 * tr_yyy_yzz[i] * tbe_0 + 4.0 * tr_yyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_yzzz[i] = 2.0 * tr_0_yzzz[i] + 2.0 * tr_y_zzz[i] - 4.0 * tr_y_yyzzz[i] * tke_0 - 10.0 * tr_yy_yzzz[i] * tbe_0 - 2.0 * tr_yyy_zzz[i] * tbe_0 + 4.0 * tr_yyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yy_zzzz[i] = 2.0 * tr_0_zzzz[i] - 4.0 * tr_y_yzzzz[i] * tke_0 - 10.0 * tr_yy_zzzz[i] * tbe_0 + 4.0 * tr_yyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-435 components of targeted buffer : DG

    auto tr_y_0_y_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 420);

    auto tr_y_0_y_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 421);

    auto tr_y_0_y_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 422);

    auto tr_y_0_y_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 423);

    auto tr_y_0_y_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 424);

    auto tr_y_0_y_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 425);

    auto tr_y_0_y_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 426);

    auto tr_y_0_y_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 427);

    auto tr_y_0_y_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 428);

    auto tr_y_0_y_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 429);

    auto tr_y_0_y_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 430);

    auto tr_y_0_y_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 431);

    auto tr_y_0_y_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 432);

    auto tr_y_0_y_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 433);

    auto tr_y_0_y_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 434);

    #pragma omp simd aligned(tr_y_0_y_yz_xxxx, tr_y_0_y_yz_xxxy, tr_y_0_y_yz_xxxz, tr_y_0_y_yz_xxyy, tr_y_0_y_yz_xxyz, tr_y_0_y_yz_xxzz, tr_y_0_y_yz_xyyy, tr_y_0_y_yz_xyyz, tr_y_0_y_yz_xyzz, tr_y_0_y_yz_xzzz, tr_y_0_y_yz_yyyy, tr_y_0_y_yz_yyyz, tr_y_0_y_yz_yyzz, tr_y_0_y_yz_yzzz, tr_y_0_y_yz_zzzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_z_xxx, tr_z_xxxxy, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_yyy, tr_z_yyyyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yz_xxxx[i] = -2.0 * tr_z_xxxxy[i] * tke_0 - 6.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_yyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xxxy[i] = tr_z_xxx[i] - 2.0 * tr_z_xxxyy[i] * tke_0 - 6.0 * tr_yz_xxxy[i] * tbe_0 - 2.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xxxz[i] = -2.0 * tr_z_xxxyz[i] * tke_0 - 6.0 * tr_yz_xxxz[i] * tbe_0 + 4.0 * tr_yyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xxyy[i] = 2.0 * tr_z_xxy[i] - 2.0 * tr_z_xxyyy[i] * tke_0 - 6.0 * tr_yz_xxyy[i] * tbe_0 - 4.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xxyz[i] = tr_z_xxz[i] - 2.0 * tr_z_xxyyz[i] * tke_0 - 6.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_yyz_xxz[i] * tbe_0 + 4.0 * tr_yyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xxzz[i] = -2.0 * tr_z_xxyzz[i] * tke_0 - 6.0 * tr_yz_xxzz[i] * tbe_0 + 4.0 * tr_yyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xyyy[i] = 3.0 * tr_z_xyy[i] - 2.0 * tr_z_xyyyy[i] * tke_0 - 6.0 * tr_yz_xyyy[i] * tbe_0 - 6.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xyyz[i] = 2.0 * tr_z_xyz[i] - 2.0 * tr_z_xyyyz[i] * tke_0 - 6.0 * tr_yz_xyyz[i] * tbe_0 - 4.0 * tr_yyz_xyz[i] * tbe_0 + 4.0 * tr_yyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xyzz[i] = tr_z_xzz[i] - 2.0 * tr_z_xyyzz[i] * tke_0 - 6.0 * tr_yz_xyzz[i] * tbe_0 - 2.0 * tr_yyz_xzz[i] * tbe_0 + 4.0 * tr_yyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_xzzz[i] = -2.0 * tr_z_xyzzz[i] * tke_0 - 6.0 * tr_yz_xzzz[i] * tbe_0 + 4.0 * tr_yyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_yyyy[i] = 4.0 * tr_z_yyy[i] - 2.0 * tr_z_yyyyy[i] * tke_0 - 6.0 * tr_yz_yyyy[i] * tbe_0 - 8.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_yyyz[i] = 3.0 * tr_z_yyz[i] - 2.0 * tr_z_yyyyz[i] * tke_0 - 6.0 * tr_yz_yyyz[i] * tbe_0 - 6.0 * tr_yyz_yyz[i] * tbe_0 + 4.0 * tr_yyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_yyzz[i] = 2.0 * tr_z_yzz[i] - 2.0 * tr_z_yyyzz[i] * tke_0 - 6.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_yyz_yzz[i] * tbe_0 + 4.0 * tr_yyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_yzzz[i] = tr_z_zzz[i] - 2.0 * tr_z_yyzzz[i] * tke_0 - 6.0 * tr_yz_yzzz[i] * tbe_0 - 2.0 * tr_yyz_zzz[i] * tbe_0 + 4.0 * tr_yyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yz_zzzz[i] = -2.0 * tr_z_yzzzz[i] * tke_0 - 6.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_yyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 435-450 components of targeted buffer : DG

    auto tr_y_0_y_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 435);

    auto tr_y_0_y_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 436);

    auto tr_y_0_y_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 437);

    auto tr_y_0_y_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 438);

    auto tr_y_0_y_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 439);

    auto tr_y_0_y_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 440);

    auto tr_y_0_y_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 441);

    auto tr_y_0_y_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 442);

    auto tr_y_0_y_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 443);

    auto tr_y_0_y_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 444);

    auto tr_y_0_y_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 445);

    auto tr_y_0_y_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 446);

    auto tr_y_0_y_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 447);

    auto tr_y_0_y_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 448);

    auto tr_y_0_y_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 449);

    #pragma omp simd aligned(tr_y_0_y_zz_xxxx, tr_y_0_y_zz_xxxy, tr_y_0_y_zz_xxxz, tr_y_0_y_zz_xxyy, tr_y_0_y_zz_xxyz, tr_y_0_y_zz_xxzz, tr_y_0_y_zz_xyyy, tr_y_0_y_zz_xyyz, tr_y_0_y_zz_xyzz, tr_y_0_y_zz_xzzz, tr_y_0_y_zz_yyyy, tr_y_0_y_zz_yyyz, tr_y_0_y_zz_yyzz, tr_y_0_y_zz_yzzz, tr_y_0_y_zz_zzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zz_xxxx[i] = -2.0 * tr_zz_xxxx[i] * tbe_0 + 4.0 * tr_yzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xxxy[i] = -2.0 * tr_zz_xxxy[i] * tbe_0 - 2.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xxxz[i] = -2.0 * tr_zz_xxxz[i] * tbe_0 + 4.0 * tr_yzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xxyy[i] = -2.0 * tr_zz_xxyy[i] * tbe_0 - 4.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xxyz[i] = -2.0 * tr_zz_xxyz[i] * tbe_0 - 2.0 * tr_yzz_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xxzz[i] = -2.0 * tr_zz_xxzz[i] * tbe_0 + 4.0 * tr_yzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xyyy[i] = -2.0 * tr_zz_xyyy[i] * tbe_0 - 6.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xyyz[i] = -2.0 * tr_zz_xyyz[i] * tbe_0 - 4.0 * tr_yzz_xyz[i] * tbe_0 + 4.0 * tr_yzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xyzz[i] = -2.0 * tr_zz_xyzz[i] * tbe_0 - 2.0 * tr_yzz_xzz[i] * tbe_0 + 4.0 * tr_yzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_xzzz[i] = -2.0 * tr_zz_xzzz[i] * tbe_0 + 4.0 * tr_yzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_yyyy[i] = -2.0 * tr_zz_yyyy[i] * tbe_0 - 8.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_yyyz[i] = -2.0 * tr_zz_yyyz[i] * tbe_0 - 6.0 * tr_yzz_yyz[i] * tbe_0 + 4.0 * tr_yzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_yyzz[i] = -2.0 * tr_zz_yyzz[i] * tbe_0 - 4.0 * tr_yzz_yzz[i] * tbe_0 + 4.0 * tr_yzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_yzzz[i] = -2.0 * tr_zz_yzzz[i] * tbe_0 - 2.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_yzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zz_zzzz[i] = -2.0 * tr_zz_zzzz[i] * tbe_0 + 4.0 * tr_yzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-465 components of targeted buffer : DG

    auto tr_y_0_z_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 450);

    auto tr_y_0_z_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 451);

    auto tr_y_0_z_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 452);

    auto tr_y_0_z_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 453);

    auto tr_y_0_z_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 454);

    auto tr_y_0_z_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 455);

    auto tr_y_0_z_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 456);

    auto tr_y_0_z_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 457);

    auto tr_y_0_z_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 458);

    auto tr_y_0_z_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 459);

    auto tr_y_0_z_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 460);

    auto tr_y_0_z_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 461);

    auto tr_y_0_z_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 462);

    auto tr_y_0_z_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 463);

    auto tr_y_0_z_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 464);

    #pragma omp simd aligned(tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_y_0_z_xx_xxxx, tr_y_0_z_xx_xxxy, tr_y_0_z_xx_xxxz, tr_y_0_z_xx_xxyy, tr_y_0_z_xx_xxyz, tr_y_0_z_xx_xxzz, tr_y_0_z_xx_xyyy, tr_y_0_z_xx_xyyz, tr_y_0_z_xx_xyzz, tr_y_0_z_xx_xzzz, tr_y_0_z_xx_yyyy, tr_y_0_z_xx_yyyz, tr_y_0_z_xx_yyzz, tr_y_0_z_xx_yzzz, tr_y_0_z_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xx_xxxx[i] = 4.0 * tr_xxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xxxy[i] = 4.0 * tr_xxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xxxz[i] = -2.0 * tr_xxy_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xxyy[i] = 4.0 * tr_xxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xxyz[i] = -2.0 * tr_xxy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xxzz[i] = -4.0 * tr_xxy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xyyy[i] = 4.0 * tr_xxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xyyz[i] = -2.0 * tr_xxy_xyy[i] * tbe_0 + 4.0 * tr_xxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xyzz[i] = -4.0 * tr_xxy_xyz[i] * tbe_0 + 4.0 * tr_xxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_xzzz[i] = -6.0 * tr_xxy_xzz[i] * tbe_0 + 4.0 * tr_xxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_yyyy[i] = 4.0 * tr_xxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_yyyz[i] = -2.0 * tr_xxy_yyy[i] * tbe_0 + 4.0 * tr_xxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_yyzz[i] = -4.0 * tr_xxy_yyz[i] * tbe_0 + 4.0 * tr_xxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_yzzz[i] = -6.0 * tr_xxy_yzz[i] * tbe_0 + 4.0 * tr_xxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xx_zzzz[i] = -8.0 * tr_xxy_zzz[i] * tbe_0 + 4.0 * tr_xxy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 465-480 components of targeted buffer : DG

    auto tr_y_0_z_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 465);

    auto tr_y_0_z_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 466);

    auto tr_y_0_z_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 467);

    auto tr_y_0_z_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 468);

    auto tr_y_0_z_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 469);

    auto tr_y_0_z_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 470);

    auto tr_y_0_z_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 471);

    auto tr_y_0_z_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 472);

    auto tr_y_0_z_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 473);

    auto tr_y_0_z_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 474);

    auto tr_y_0_z_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 475);

    auto tr_y_0_z_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 476);

    auto tr_y_0_z_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 477);

    auto tr_y_0_z_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 478);

    auto tr_y_0_z_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 479);

    #pragma omp simd aligned(tr_x_xxx, tr_x_xxxxz, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_x_zzzzz, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_y_0_z_xy_xxxx, tr_y_0_z_xy_xxxy, tr_y_0_z_xy_xxxz, tr_y_0_z_xy_xxyy, tr_y_0_z_xy_xxyz, tr_y_0_z_xy_xxzz, tr_y_0_z_xy_xyyy, tr_y_0_z_xy_xyyz, tr_y_0_z_xy_xyzz, tr_y_0_z_xy_xzzz, tr_y_0_z_xy_yyyy, tr_y_0_z_xy_yyyz, tr_y_0_z_xy_yyzz, tr_y_0_z_xy_yzzz, tr_y_0_z_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xy_xxxx[i] = -2.0 * tr_x_xxxxz[i] * tke_0 - 2.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xxxy[i] = -2.0 * tr_x_xxxyz[i] * tke_0 - 2.0 * tr_xz_xxxy[i] * tbe_0 + 4.0 * tr_xyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xxxz[i] = tr_x_xxx[i] - 2.0 * tr_x_xxxzz[i] * tke_0 - 2.0 * tr_xz_xxxz[i] * tbe_0 - 2.0 * tr_xyy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xxyy[i] = -2.0 * tr_x_xxyyz[i] * tke_0 - 2.0 * tr_xz_xxyy[i] * tbe_0 + 4.0 * tr_xyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xxyz[i] = tr_x_xxy[i] - 2.0 * tr_x_xxyzz[i] * tke_0 - 2.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xyy_xxy[i] * tbe_0 + 4.0 * tr_xyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xxzz[i] = 2.0 * tr_x_xxz[i] - 2.0 * tr_x_xxzzz[i] * tke_0 - 2.0 * tr_xz_xxzz[i] * tbe_0 - 4.0 * tr_xyy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xyyy[i] = -2.0 * tr_x_xyyyz[i] * tke_0 - 2.0 * tr_xz_xyyy[i] * tbe_0 + 4.0 * tr_xyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xyyz[i] = tr_x_xyy[i] - 2.0 * tr_x_xyyzz[i] * tke_0 - 2.0 * tr_xz_xyyz[i] * tbe_0 - 2.0 * tr_xyy_xyy[i] * tbe_0 + 4.0 * tr_xyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xyzz[i] = 2.0 * tr_x_xyz[i] - 2.0 * tr_x_xyzzz[i] * tke_0 - 2.0 * tr_xz_xyzz[i] * tbe_0 - 4.0 * tr_xyy_xyz[i] * tbe_0 + 4.0 * tr_xyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_xzzz[i] = 3.0 * tr_x_xzz[i] - 2.0 * tr_x_xzzzz[i] * tke_0 - 2.0 * tr_xz_xzzz[i] * tbe_0 - 6.0 * tr_xyy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_yyyy[i] = -2.0 * tr_x_yyyyz[i] * tke_0 - 2.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_yyyz[i] = tr_x_yyy[i] - 2.0 * tr_x_yyyzz[i] * tke_0 - 2.0 * tr_xz_yyyz[i] * tbe_0 - 2.0 * tr_xyy_yyy[i] * tbe_0 + 4.0 * tr_xyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_yyzz[i] = 2.0 * tr_x_yyz[i] - 2.0 * tr_x_yyzzz[i] * tke_0 - 2.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xyy_yyz[i] * tbe_0 + 4.0 * tr_xyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_yzzz[i] = 3.0 * tr_x_yzz[i] - 2.0 * tr_x_yzzzz[i] * tke_0 - 2.0 * tr_xz_yzzz[i] * tbe_0 - 6.0 * tr_xyy_yzz[i] * tbe_0 + 4.0 * tr_xyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xy_zzzz[i] = 4.0 * tr_x_zzz[i] - 2.0 * tr_x_zzzzz[i] * tke_0 - 2.0 * tr_xz_zzzz[i] * tbe_0 - 8.0 * tr_xyy_zzz[i] * tbe_0 + 4.0 * tr_xyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-495 components of targeted buffer : DG

    auto tr_y_0_z_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 480);

    auto tr_y_0_z_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 481);

    auto tr_y_0_z_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 482);

    auto tr_y_0_z_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 483);

    auto tr_y_0_z_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 484);

    auto tr_y_0_z_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 485);

    auto tr_y_0_z_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 486);

    auto tr_y_0_z_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 487);

    auto tr_y_0_z_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 488);

    auto tr_y_0_z_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 489);

    auto tr_y_0_z_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 490);

    auto tr_y_0_z_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 491);

    auto tr_y_0_z_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 492);

    auto tr_y_0_z_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 493);

    auto tr_y_0_z_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 494);

    #pragma omp simd aligned(tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_y_0_z_xz_xxxx, tr_y_0_z_xz_xxxy, tr_y_0_z_xz_xxxz, tr_y_0_z_xz_xxyy, tr_y_0_z_xz_xxyz, tr_y_0_z_xz_xxzz, tr_y_0_z_xz_xyyy, tr_y_0_z_xz_xyyz, tr_y_0_z_xz_xyzz, tr_y_0_z_xz_xzzz, tr_y_0_z_xz_yyyy, tr_y_0_z_xz_yyyz, tr_y_0_z_xz_yyzz, tr_y_0_z_xz_yzzz, tr_y_0_z_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xz_xxxx[i] = -2.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xxxy[i] = -2.0 * tr_xy_xxxy[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xxxz[i] = -2.0 * tr_xy_xxxz[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xxyy[i] = -2.0 * tr_xy_xxyy[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xxyz[i] = -2.0 * tr_xy_xxyz[i] * tbe_0 - 2.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xxzz[i] = -2.0 * tr_xy_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xyyy[i] = -2.0 * tr_xy_xyyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xyyz[i] = -2.0 * tr_xy_xyyz[i] * tbe_0 - 2.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xyzz[i] = -2.0 * tr_xy_xyzz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_xzzz[i] = -2.0 * tr_xy_xzzz[i] * tbe_0 - 6.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_yyyy[i] = -2.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_yyyz[i] = -2.0 * tr_xy_yyyz[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_yyzz[i] = -2.0 * tr_xy_yyzz[i] * tbe_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_yzzz[i] = -2.0 * tr_xy_yzzz[i] * tbe_0 - 6.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xz_zzzz[i] = -2.0 * tr_xy_zzzz[i] * tbe_0 - 8.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 495-510 components of targeted buffer : DG

    auto tr_y_0_z_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 495);

    auto tr_y_0_z_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 496);

    auto tr_y_0_z_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 497);

    auto tr_y_0_z_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 498);

    auto tr_y_0_z_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 499);

    auto tr_y_0_z_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 500);

    auto tr_y_0_z_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 501);

    auto tr_y_0_z_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 502);

    auto tr_y_0_z_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 503);

    auto tr_y_0_z_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 504);

    auto tr_y_0_z_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 505);

    auto tr_y_0_z_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 506);

    auto tr_y_0_z_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 507);

    auto tr_y_0_z_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 508);

    auto tr_y_0_z_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 509);

    #pragma omp simd aligned(tr_y_0_z_yy_xxxx, tr_y_0_z_yy_xxxy, tr_y_0_z_yy_xxxz, tr_y_0_z_yy_xxyy, tr_y_0_z_yy_xxyz, tr_y_0_z_yy_xxzz, tr_y_0_z_yy_xyyy, tr_y_0_z_yy_xyyz, tr_y_0_z_yy_xyzz, tr_y_0_z_yy_xzzz, tr_y_0_z_yy_yyyy, tr_y_0_z_yy_yyyz, tr_y_0_z_yy_yyzz, tr_y_0_z_yy_yzzz, tr_y_0_z_yy_zzzz, tr_y_xxx, tr_y_xxxxz, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_y_zzzzz, tr_yyy_xxx, tr_yyy_xxxxz, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyy_zzzzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yy_xxxx[i] = -4.0 * tr_y_xxxxz[i] * tke_0 - 4.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_yyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xxxy[i] = -4.0 * tr_y_xxxyz[i] * tke_0 - 4.0 * tr_yz_xxxy[i] * tbe_0 + 4.0 * tr_yyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xxxz[i] = 2.0 * tr_y_xxx[i] - 4.0 * tr_y_xxxzz[i] * tke_0 - 4.0 * tr_yz_xxxz[i] * tbe_0 - 2.0 * tr_yyy_xxx[i] * tbe_0 + 4.0 * tr_yyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xxyy[i] = -4.0 * tr_y_xxyyz[i] * tke_0 - 4.0 * tr_yz_xxyy[i] * tbe_0 + 4.0 * tr_yyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xxyz[i] = 2.0 * tr_y_xxy[i] - 4.0 * tr_y_xxyzz[i] * tke_0 - 4.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_yyy_xxy[i] * tbe_0 + 4.0 * tr_yyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xxzz[i] = 4.0 * tr_y_xxz[i] - 4.0 * tr_y_xxzzz[i] * tke_0 - 4.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_yyy_xxz[i] * tbe_0 + 4.0 * tr_yyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xyyy[i] = -4.0 * tr_y_xyyyz[i] * tke_0 - 4.0 * tr_yz_xyyy[i] * tbe_0 + 4.0 * tr_yyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xyyz[i] = 2.0 * tr_y_xyy[i] - 4.0 * tr_y_xyyzz[i] * tke_0 - 4.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_yyy_xyy[i] * tbe_0 + 4.0 * tr_yyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xyzz[i] = 4.0 * tr_y_xyz[i] - 4.0 * tr_y_xyzzz[i] * tke_0 - 4.0 * tr_yz_xyzz[i] * tbe_0 - 4.0 * tr_yyy_xyz[i] * tbe_0 + 4.0 * tr_yyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_xzzz[i] = 6.0 * tr_y_xzz[i] - 4.0 * tr_y_xzzzz[i] * tke_0 - 4.0 * tr_yz_xzzz[i] * tbe_0 - 6.0 * tr_yyy_xzz[i] * tbe_0 + 4.0 * tr_yyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_yyyy[i] = -4.0 * tr_y_yyyyz[i] * tke_0 - 4.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_yyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_yyyz[i] = 2.0 * tr_y_yyy[i] - 4.0 * tr_y_yyyzz[i] * tke_0 - 4.0 * tr_yz_yyyz[i] * tbe_0 - 2.0 * tr_yyy_yyy[i] * tbe_0 + 4.0 * tr_yyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_yyzz[i] = 4.0 * tr_y_yyz[i] - 4.0 * tr_y_yyzzz[i] * tke_0 - 4.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_yyy_yyz[i] * tbe_0 + 4.0 * tr_yyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_yzzz[i] = 6.0 * tr_y_yzz[i] - 4.0 * tr_y_yzzzz[i] * tke_0 - 4.0 * tr_yz_yzzz[i] * tbe_0 - 6.0 * tr_yyy_yzz[i] * tbe_0 + 4.0 * tr_yyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yy_zzzz[i] = 8.0 * tr_y_zzz[i] - 4.0 * tr_y_zzzzz[i] * tke_0 - 4.0 * tr_yz_zzzz[i] * tbe_0 - 8.0 * tr_yyy_zzz[i] * tbe_0 + 4.0 * tr_yyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-525 components of targeted buffer : DG

    auto tr_y_0_z_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 510);

    auto tr_y_0_z_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 511);

    auto tr_y_0_z_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 512);

    auto tr_y_0_z_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 513);

    auto tr_y_0_z_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 514);

    auto tr_y_0_z_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 515);

    auto tr_y_0_z_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 516);

    auto tr_y_0_z_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 517);

    auto tr_y_0_z_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 518);

    auto tr_y_0_z_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 519);

    auto tr_y_0_z_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 520);

    auto tr_y_0_z_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 521);

    auto tr_y_0_z_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 522);

    auto tr_y_0_z_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 523);

    auto tr_y_0_z_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 524);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_y_0_z_yz_xxxx, tr_y_0_z_yz_xxxy, tr_y_0_z_yz_xxxz, tr_y_0_z_yz_xxyy, tr_y_0_z_yz_xxyz, tr_y_0_z_yz_xxzz, tr_y_0_z_yz_xyyy, tr_y_0_z_yz_xyyz, tr_y_0_z_yz_xyzz, tr_y_0_z_yz_xzzz, tr_y_0_z_yz_yyyy, tr_y_0_z_yz_yyyz, tr_y_0_z_yz_yyzz, tr_y_0_z_yz_yzzz, tr_y_0_z_yz_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_z_xxx, tr_z_xxxxz, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, tr_z_zzzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yz_xxxx[i] = tr_0_xxxx[i] - 2.0 * tr_z_xxxxz[i] * tke_0 - 2.0 * tr_zz_xxxx[i] * tbe_0 - 2.0 * tr_yy_xxxx[i] * tbe_0 + 4.0 * tr_yyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xxxy[i] = tr_0_xxxy[i] - 2.0 * tr_z_xxxyz[i] * tke_0 - 2.0 * tr_zz_xxxy[i] * tbe_0 - 2.0 * tr_yy_xxxy[i] * tbe_0 + 4.0 * tr_yyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xxxz[i] = tr_0_xxxz[i] + tr_z_xxx[i] - 2.0 * tr_z_xxxzz[i] * tke_0 - 2.0 * tr_zz_xxxz[i] * tbe_0 - 2.0 * tr_yy_xxxz[i] * tbe_0 - 2.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xxyy[i] = tr_0_xxyy[i] - 2.0 * tr_z_xxyyz[i] * tke_0 - 2.0 * tr_zz_xxyy[i] * tbe_0 - 2.0 * tr_yy_xxyy[i] * tbe_0 + 4.0 * tr_yyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xxyz[i] = tr_0_xxyz[i] + tr_z_xxy[i] - 2.0 * tr_z_xxyzz[i] * tke_0 - 2.0 * tr_zz_xxyz[i] * tbe_0 - 2.0 * tr_yy_xxyz[i] * tbe_0 - 2.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xxzz[i] = tr_0_xxzz[i] + 2.0 * tr_z_xxz[i] - 2.0 * tr_z_xxzzz[i] * tke_0 - 2.0 * tr_zz_xxzz[i] * tbe_0 - 2.0 * tr_yy_xxzz[i] * tbe_0 - 4.0 * tr_yyz_xxz[i] * tbe_0 + 4.0 * tr_yyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xyyy[i] = tr_0_xyyy[i] - 2.0 * tr_z_xyyyz[i] * tke_0 - 2.0 * tr_zz_xyyy[i] * tbe_0 - 2.0 * tr_yy_xyyy[i] * tbe_0 + 4.0 * tr_yyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xyyz[i] = tr_0_xyyz[i] + tr_z_xyy[i] - 2.0 * tr_z_xyyzz[i] * tke_0 - 2.0 * tr_zz_xyyz[i] * tbe_0 - 2.0 * tr_yy_xyyz[i] * tbe_0 - 2.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xyzz[i] = tr_0_xyzz[i] + 2.0 * tr_z_xyz[i] - 2.0 * tr_z_xyzzz[i] * tke_0 - 2.0 * tr_zz_xyzz[i] * tbe_0 - 2.0 * tr_yy_xyzz[i] * tbe_0 - 4.0 * tr_yyz_xyz[i] * tbe_0 + 4.0 * tr_yyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_xzzz[i] = tr_0_xzzz[i] + 3.0 * tr_z_xzz[i] - 2.0 * tr_z_xzzzz[i] * tke_0 - 2.0 * tr_zz_xzzz[i] * tbe_0 - 2.0 * tr_yy_xzzz[i] * tbe_0 - 6.0 * tr_yyz_xzz[i] * tbe_0 + 4.0 * tr_yyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_yyyy[i] = tr_0_yyyy[i] - 2.0 * tr_z_yyyyz[i] * tke_0 - 2.0 * tr_zz_yyyy[i] * tbe_0 - 2.0 * tr_yy_yyyy[i] * tbe_0 + 4.0 * tr_yyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_yyyz[i] = tr_0_yyyz[i] + tr_z_yyy[i] - 2.0 * tr_z_yyyzz[i] * tke_0 - 2.0 * tr_zz_yyyz[i] * tbe_0 - 2.0 * tr_yy_yyyz[i] * tbe_0 - 2.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_yyzz[i] = tr_0_yyzz[i] + 2.0 * tr_z_yyz[i] - 2.0 * tr_z_yyzzz[i] * tke_0 - 2.0 * tr_zz_yyzz[i] * tbe_0 - 2.0 * tr_yy_yyzz[i] * tbe_0 - 4.0 * tr_yyz_yyz[i] * tbe_0 + 4.0 * tr_yyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_yzzz[i] = tr_0_yzzz[i] + 3.0 * tr_z_yzz[i] - 2.0 * tr_z_yzzzz[i] * tke_0 - 2.0 * tr_zz_yzzz[i] * tbe_0 - 2.0 * tr_yy_yzzz[i] * tbe_0 - 6.0 * tr_yyz_yzz[i] * tbe_0 + 4.0 * tr_yyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yz_zzzz[i] = tr_0_zzzz[i] + 4.0 * tr_z_zzz[i] - 2.0 * tr_z_zzzzz[i] * tke_0 - 2.0 * tr_zz_zzzz[i] * tbe_0 - 2.0 * tr_yy_zzzz[i] * tbe_0 - 8.0 * tr_yyz_zzz[i] * tbe_0 + 4.0 * tr_yyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 525-540 components of targeted buffer : DG

    auto tr_y_0_z_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 525);

    auto tr_y_0_z_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 526);

    auto tr_y_0_z_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 527);

    auto tr_y_0_z_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 528);

    auto tr_y_0_z_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 529);

    auto tr_y_0_z_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 530);

    auto tr_y_0_z_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 531);

    auto tr_y_0_z_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 532);

    auto tr_y_0_z_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 533);

    auto tr_y_0_z_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 534);

    auto tr_y_0_z_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 535);

    auto tr_y_0_z_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 536);

    auto tr_y_0_z_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 537);

    auto tr_y_0_z_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 538);

    auto tr_y_0_z_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 539);

    #pragma omp simd aligned(tr_y_0_z_zz_xxxx, tr_y_0_z_zz_xxxy, tr_y_0_z_zz_xxxz, tr_y_0_z_zz_xxyy, tr_y_0_z_zz_xxyz, tr_y_0_z_zz_xxzz, tr_y_0_z_zz_xyyy, tr_y_0_z_zz_xyyz, tr_y_0_z_zz_xyzz, tr_y_0_z_zz_xzzz, tr_y_0_z_zz_yyyy, tr_y_0_z_zz_yyyz, tr_y_0_z_zz_yyzz, tr_y_0_z_zz_yzzz, tr_y_0_z_zz_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zz_xxxx[i] = -4.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_yzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xxxy[i] = -4.0 * tr_yz_xxxy[i] * tbe_0 + 4.0 * tr_yzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xxxz[i] = -4.0 * tr_yz_xxxz[i] * tbe_0 - 2.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xxyy[i] = -4.0 * tr_yz_xxyy[i] * tbe_0 + 4.0 * tr_yzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xxyz[i] = -4.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xxzz[i] = -4.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_yzz_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xyyy[i] = -4.0 * tr_yz_xyyy[i] * tbe_0 + 4.0 * tr_yzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xyyz[i] = -4.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xyzz[i] = -4.0 * tr_yz_xyzz[i] * tbe_0 - 4.0 * tr_yzz_xyz[i] * tbe_0 + 4.0 * tr_yzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_xzzz[i] = -4.0 * tr_yz_xzzz[i] * tbe_0 - 6.0 * tr_yzz_xzz[i] * tbe_0 + 4.0 * tr_yzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_yyyy[i] = -4.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_yzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_yyyz[i] = -4.0 * tr_yz_yyyz[i] * tbe_0 - 2.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_yyzz[i] = -4.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_yzz_yyz[i] * tbe_0 + 4.0 * tr_yzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_yzzz[i] = -4.0 * tr_yz_yzzz[i] * tbe_0 - 6.0 * tr_yzz_yzz[i] * tbe_0 + 4.0 * tr_yzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zz_zzzz[i] = -4.0 * tr_yz_zzzz[i] * tbe_0 - 8.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_yzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 540-555 components of targeted buffer : DG

    auto tr_z_0_x_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 540);

    auto tr_z_0_x_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 541);

    auto tr_z_0_x_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 542);

    auto tr_z_0_x_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 543);

    auto tr_z_0_x_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 544);

    auto tr_z_0_x_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 545);

    auto tr_z_0_x_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 546);

    auto tr_z_0_x_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 547);

    auto tr_z_0_x_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 548);

    auto tr_z_0_x_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 549);

    auto tr_z_0_x_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 550);

    auto tr_z_0_x_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 551);

    auto tr_z_0_x_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 552);

    auto tr_z_0_x_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 553);

    auto tr_z_0_x_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 554);

    #pragma omp simd aligned(tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_z_0_x_xx_xxxx, tr_z_0_x_xx_xxxy, tr_z_0_x_xx_xxxz, tr_z_0_x_xx_xxyy, tr_z_0_x_xx_xxyz, tr_z_0_x_xx_xxzz, tr_z_0_x_xx_xyyy, tr_z_0_x_xx_xyyz, tr_z_0_x_xx_xyzz, tr_z_0_x_xx_xzzz, tr_z_0_x_xx_yyyy, tr_z_0_x_xx_yyyz, tr_z_0_x_xx_yyzz, tr_z_0_x_xx_yzzz, tr_z_0_x_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xx_xxxx[i] = -4.0 * tr_xz_xxxx[i] * tbe_0 - 8.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xxxy[i] = -4.0 * tr_xz_xxxy[i] * tbe_0 - 6.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xxxz[i] = -4.0 * tr_xz_xxxz[i] * tbe_0 - 6.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xxyy[i] = -4.0 * tr_xz_xxyy[i] * tbe_0 - 4.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xxyz[i] = -4.0 * tr_xz_xxyz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tbe_0 + 4.0 * tr_xxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xxzz[i] = -4.0 * tr_xz_xxzz[i] * tbe_0 - 4.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xyyy[i] = -4.0 * tr_xz_xyyy[i] * tbe_0 - 2.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xyyz[i] = -4.0 * tr_xz_xyyz[i] * tbe_0 - 2.0 * tr_xxz_yyz[i] * tbe_0 + 4.0 * tr_xxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xyzz[i] = -4.0 * tr_xz_xyzz[i] * tbe_0 - 2.0 * tr_xxz_yzz[i] * tbe_0 + 4.0 * tr_xxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_xzzz[i] = -4.0 * tr_xz_xzzz[i] * tbe_0 - 2.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_yyyy[i] = -4.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_yyyz[i] = -4.0 * tr_xz_yyyz[i] * tbe_0 + 4.0 * tr_xxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_yyzz[i] = -4.0 * tr_xz_yyzz[i] * tbe_0 + 4.0 * tr_xxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_yzzz[i] = -4.0 * tr_xz_yzzz[i] * tbe_0 + 4.0 * tr_xxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xx_zzzz[i] = -4.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 555-570 components of targeted buffer : DG

    auto tr_z_0_x_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 555);

    auto tr_z_0_x_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 556);

    auto tr_z_0_x_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 557);

    auto tr_z_0_x_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 558);

    auto tr_z_0_x_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 559);

    auto tr_z_0_x_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 560);

    auto tr_z_0_x_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 561);

    auto tr_z_0_x_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 562);

    auto tr_z_0_x_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 563);

    auto tr_z_0_x_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 564);

    auto tr_z_0_x_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 565);

    auto tr_z_0_x_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 566);

    auto tr_z_0_x_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 567);

    auto tr_z_0_x_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 568);

    auto tr_z_0_x_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 569);

    #pragma omp simd aligned(tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_z_0_x_xy_xxxx, tr_z_0_x_xy_xxxy, tr_z_0_x_xy_xxxz, tr_z_0_x_xy_xxyy, tr_z_0_x_xy_xxyz, tr_z_0_x_xy_xxzz, tr_z_0_x_xy_xyyy, tr_z_0_x_xy_xyyz, tr_z_0_x_xy_xyzz, tr_z_0_x_xy_xzzz, tr_z_0_x_xy_yyyy, tr_z_0_x_xy_yyyz, tr_z_0_x_xy_yyzz, tr_z_0_x_xy_yzzz, tr_z_0_x_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xy_xxxx[i] = -2.0 * tr_yz_xxxx[i] * tbe_0 - 8.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xxxy[i] = -2.0 * tr_yz_xxxy[i] * tbe_0 - 6.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xxxz[i] = -2.0 * tr_yz_xxxz[i] * tbe_0 - 6.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xxyy[i] = -2.0 * tr_yz_xxyy[i] * tbe_0 - 4.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xxyz[i] = -2.0 * tr_yz_xxyz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xxzz[i] = -2.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xyyy[i] = -2.0 * tr_yz_xyyy[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xyyz[i] = -2.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xyzz[i] = -2.0 * tr_yz_xyzz[i] * tbe_0 - 2.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_xzzz[i] = -2.0 * tr_yz_xzzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_yyyy[i] = -2.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_yyyz[i] = -2.0 * tr_yz_yyyz[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_yyzz[i] = -2.0 * tr_yz_yyzz[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_yzzz[i] = -2.0 * tr_yz_yzzz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xy_zzzz[i] = -2.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_xyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 570-585 components of targeted buffer : DG

    auto tr_z_0_x_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 570);

    auto tr_z_0_x_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 571);

    auto tr_z_0_x_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 572);

    auto tr_z_0_x_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 573);

    auto tr_z_0_x_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 574);

    auto tr_z_0_x_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 575);

    auto tr_z_0_x_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 576);

    auto tr_z_0_x_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 577);

    auto tr_z_0_x_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 578);

    auto tr_z_0_x_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 579);

    auto tr_z_0_x_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 580);

    auto tr_z_0_x_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 581);

    auto tr_z_0_x_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 582);

    auto tr_z_0_x_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 583);

    auto tr_z_0_x_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 584);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_x_xxx, tr_x_xxxxx, tr_x_xxxxy, tr_x_xxxxz, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_z_0_x_xz_xxxx, tr_z_0_x_xz_xxxy, tr_z_0_x_xz_xxxz, tr_z_0_x_xz_xxyy, tr_z_0_x_xz_xxyz, tr_z_0_x_xz_xxzz, tr_z_0_x_xz_xyyy, tr_z_0_x_xz_xyyz, tr_z_0_x_xz_xyzz, tr_z_0_x_xz_xzzz, tr_z_0_x_xz_yyyy, tr_z_0_x_xz_yyyz, tr_z_0_x_xz_yyzz, tr_z_0_x_xz_yzzz, tr_z_0_x_xz_zzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xz_xxxx[i] = tr_0_xxxx[i] - 2.0 * tr_zz_xxxx[i] * tbe_0 + 4.0 * tr_x_xxx[i] - 2.0 * tr_x_xxxxx[i] * tke_0 - 8.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xxxy[i] = tr_0_xxxy[i] - 2.0 * tr_zz_xxxy[i] * tbe_0 + 3.0 * tr_x_xxy[i] - 2.0 * tr_x_xxxxy[i] * tke_0 - 6.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxxy[i] * tbe_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xxxz[i] = tr_0_xxxz[i] - 2.0 * tr_zz_xxxz[i] * tbe_0 + 3.0 * tr_x_xxz[i] - 2.0 * tr_x_xxxxz[i] * tke_0 - 6.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxxz[i] * tbe_0 + 4.0 * tr_xxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xxyy[i] = tr_0_xxyy[i] - 2.0 * tr_zz_xxyy[i] * tbe_0 + 2.0 * tr_x_xyy[i] - 2.0 * tr_x_xxxyy[i] * tke_0 - 4.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxyy[i] * tbe_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xxyz[i] = tr_0_xxyz[i] - 2.0 * tr_zz_xxyz[i] * tbe_0 + 2.0 * tr_x_xyz[i] - 2.0 * tr_x_xxxyz[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tbe_0 + 4.0 * tr_xzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxyz[i] * tbe_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xxzz[i] = tr_0_xxzz[i] - 2.0 * tr_zz_xxzz[i] * tbe_0 + 2.0 * tr_x_xzz[i] - 2.0 * tr_x_xxxzz[i] * tke_0 - 4.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxzz[i] * tbe_0 + 4.0 * tr_xxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xyyy[i] = tr_0_xyyy[i] - 2.0 * tr_zz_xyyy[i] * tbe_0 + tr_x_yyy[i] - 2.0 * tr_x_xxyyy[i] * tke_0 - 2.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xyyz[i] = tr_0_xyyz[i] - 2.0 * tr_zz_xyyz[i] * tbe_0 + tr_x_yyz[i] - 2.0 * tr_x_xxyyz[i] * tke_0 - 2.0 * tr_xzz_yyz[i] * tbe_0 + 4.0 * tr_xzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyyz[i] * tbe_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xyzz[i] = tr_0_xyzz[i] - 2.0 * tr_zz_xyzz[i] * tbe_0 + tr_x_yzz[i] - 2.0 * tr_x_xxyzz[i] * tke_0 - 2.0 * tr_xzz_yzz[i] * tbe_0 + 4.0 * tr_xzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyzz[i] * tbe_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_xzzz[i] = tr_0_xzzz[i] - 2.0 * tr_zz_xzzz[i] * tbe_0 + tr_x_zzz[i] - 2.0 * tr_x_xxzzz[i] * tke_0 - 2.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xzzz[i] * tbe_0 + 4.0 * tr_xxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_yyyy[i] = tr_0_yyyy[i] - 2.0 * tr_zz_yyyy[i] * tbe_0 - 2.0 * tr_x_xyyyy[i] * tke_0 + 4.0 * tr_xzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_yyyz[i] = tr_0_yyyz[i] - 2.0 * tr_zz_yyyz[i] * tbe_0 - 2.0 * tr_x_xyyyz[i] * tke_0 + 4.0 * tr_xzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyyz[i] * tbe_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_yyzz[i] = tr_0_yyzz[i] - 2.0 * tr_zz_yyzz[i] * tbe_0 - 2.0 * tr_x_xyyzz[i] * tke_0 + 4.0 * tr_xzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyzz[i] * tbe_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_yzzz[i] = tr_0_yzzz[i] - 2.0 * tr_zz_yzzz[i] * tbe_0 - 2.0 * tr_x_xyzzz[i] * tke_0 + 4.0 * tr_xzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yzzz[i] * tbe_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xz_zzzz[i] = tr_0_zzzz[i] - 2.0 * tr_zz_zzzz[i] * tbe_0 - 2.0 * tr_x_xzzzz[i] * tke_0 + 4.0 * tr_xzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zzzz[i] * tbe_0 + 4.0 * tr_xxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 585-600 components of targeted buffer : DG

    auto tr_z_0_x_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 585);

    auto tr_z_0_x_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 586);

    auto tr_z_0_x_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 587);

    auto tr_z_0_x_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 588);

    auto tr_z_0_x_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 589);

    auto tr_z_0_x_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 590);

    auto tr_z_0_x_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 591);

    auto tr_z_0_x_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 592);

    auto tr_z_0_x_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 593);

    auto tr_z_0_x_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 594);

    auto tr_z_0_x_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 595);

    auto tr_z_0_x_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 596);

    auto tr_z_0_x_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 597);

    auto tr_z_0_x_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 598);

    auto tr_z_0_x_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 599);

    #pragma omp simd aligned(tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_z_0_x_yy_xxxx, tr_z_0_x_yy_xxxy, tr_z_0_x_yy_xxxz, tr_z_0_x_yy_xxyy, tr_z_0_x_yy_xxyz, tr_z_0_x_yy_xxzz, tr_z_0_x_yy_xyyy, tr_z_0_x_yy_xyyz, tr_z_0_x_yy_xyzz, tr_z_0_x_yy_xzzz, tr_z_0_x_yy_yyyy, tr_z_0_x_yy_yyyz, tr_z_0_x_yy_yyzz, tr_z_0_x_yy_yzzz, tr_z_0_x_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yy_xxxx[i] = -8.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xxxy[i] = -6.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xxxz[i] = -6.0 * tr_yyz_xxz[i] * tbe_0 + 4.0 * tr_yyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xxyy[i] = -4.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xxyz[i] = -4.0 * tr_yyz_xyz[i] * tbe_0 + 4.0 * tr_yyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xxzz[i] = -4.0 * tr_yyz_xzz[i] * tbe_0 + 4.0 * tr_yyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xyyy[i] = -2.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xyyz[i] = -2.0 * tr_yyz_yyz[i] * tbe_0 + 4.0 * tr_yyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xyzz[i] = -2.0 * tr_yyz_yzz[i] * tbe_0 + 4.0 * tr_yyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_xzzz[i] = -2.0 * tr_yyz_zzz[i] * tbe_0 + 4.0 * tr_yyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_yyyy[i] = 4.0 * tr_yyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_yyyz[i] = 4.0 * tr_yyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_yyzz[i] = 4.0 * tr_yyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_yzzz[i] = 4.0 * tr_yyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yy_zzzz[i] = 4.0 * tr_yyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 600-615 components of targeted buffer : DG

    auto tr_z_0_x_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 600);

    auto tr_z_0_x_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 601);

    auto tr_z_0_x_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 602);

    auto tr_z_0_x_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 603);

    auto tr_z_0_x_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 604);

    auto tr_z_0_x_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 605);

    auto tr_z_0_x_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 606);

    auto tr_z_0_x_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 607);

    auto tr_z_0_x_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 608);

    auto tr_z_0_x_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 609);

    auto tr_z_0_x_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 610);

    auto tr_z_0_x_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 611);

    auto tr_z_0_x_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 612);

    auto tr_z_0_x_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 613);

    auto tr_z_0_x_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 614);

    #pragma omp simd aligned(tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_y_xxx, tr_y_xxxxx, tr_y_xxxxy, tr_y_xxxxz, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_z_0_x_yz_xxxx, tr_z_0_x_yz_xxxy, tr_z_0_x_yz_xxxz, tr_z_0_x_yz_xxyy, tr_z_0_x_yz_xxyz, tr_z_0_x_yz_xxzz, tr_z_0_x_yz_xyyy, tr_z_0_x_yz_xyyz, tr_z_0_x_yz_xyzz, tr_z_0_x_yz_xzzz, tr_z_0_x_yz_yyyy, tr_z_0_x_yz_yyyz, tr_z_0_x_yz_yyzz, tr_z_0_x_yz_yzzz, tr_z_0_x_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yz_xxxx[i] = 4.0 * tr_y_xxx[i] - 2.0 * tr_y_xxxxx[i] * tke_0 - 8.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xxxy[i] = 3.0 * tr_y_xxy[i] - 2.0 * tr_y_xxxxy[i] * tke_0 - 6.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xxxz[i] = 3.0 * tr_y_xxz[i] - 2.0 * tr_y_xxxxz[i] * tke_0 - 6.0 * tr_yzz_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xxyy[i] = 2.0 * tr_y_xyy[i] - 2.0 * tr_y_xxxyy[i] * tke_0 - 4.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xxyz[i] = 2.0 * tr_y_xyz[i] - 2.0 * tr_y_xxxyz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tbe_0 + 4.0 * tr_yzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxyz[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xxzz[i] = 2.0 * tr_y_xzz[i] - 2.0 * tr_y_xxxzz[i] * tke_0 - 4.0 * tr_yzz_xzz[i] * tbe_0 + 4.0 * tr_yzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxzz[i] * tbe_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xyyy[i] = tr_y_yyy[i] - 2.0 * tr_y_xxyyy[i] * tke_0 - 2.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xyyz[i] = tr_y_yyz[i] - 2.0 * tr_y_xxyyz[i] * tke_0 - 2.0 * tr_yzz_yyz[i] * tbe_0 + 4.0 * tr_yzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xyzz[i] = tr_y_yzz[i] - 2.0 * tr_y_xxyzz[i] * tke_0 - 2.0 * tr_yzz_yzz[i] * tbe_0 + 4.0 * tr_yzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_xzzz[i] = tr_y_zzz[i] - 2.0 * tr_y_xxzzz[i] * tke_0 - 2.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_yzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xzzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_yyyy[i] = -2.0 * tr_y_xyyyy[i] * tke_0 + 4.0 * tr_yzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_yyyz[i] = -2.0 * tr_y_xyyyz[i] * tke_0 + 4.0 * tr_yzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyyz[i] * tbe_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_yyzz[i] = -2.0 * tr_y_xyyzz[i] * tke_0 + 4.0 * tr_yzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyzz[i] * tbe_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_yzzz[i] = -2.0 * tr_y_xyzzz[i] * tke_0 + 4.0 * tr_yzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yzzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yz_zzzz[i] = -2.0 * tr_y_xzzzz[i] * tke_0 + 4.0 * tr_yzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 615-630 components of targeted buffer : DG

    auto tr_z_0_x_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 615);

    auto tr_z_0_x_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 616);

    auto tr_z_0_x_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 617);

    auto tr_z_0_x_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 618);

    auto tr_z_0_x_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 619);

    auto tr_z_0_x_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 620);

    auto tr_z_0_x_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 621);

    auto tr_z_0_x_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 622);

    auto tr_z_0_x_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 623);

    auto tr_z_0_x_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 624);

    auto tr_z_0_x_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 625);

    auto tr_z_0_x_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 626);

    auto tr_z_0_x_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 627);

    auto tr_z_0_x_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 628);

    auto tr_z_0_x_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 629);

    #pragma omp simd aligned(tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, tr_z_0_x_zz_xxxx, tr_z_0_x_zz_xxxy, tr_z_0_x_zz_xxxz, tr_z_0_x_zz_xxyy, tr_z_0_x_zz_xxyz, tr_z_0_x_zz_xxzz, tr_z_0_x_zz_xyyy, tr_z_0_x_zz_xyyz, tr_z_0_x_zz_xyzz, tr_z_0_x_zz_xzzz, tr_z_0_x_zz_yyyy, tr_z_0_x_zz_yyyz, tr_z_0_x_zz_yyzz, tr_z_0_x_zz_yzzz, tr_z_0_x_zz_zzzz, tr_z_xxx, tr_z_xxxxx, tr_z_xxxxy, tr_z_xxxxz, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, tr_zzz_xxx, tr_zzz_xxxxx, tr_zzz_xxxxy, tr_zzz_xxxxz, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zz_xxxx[i] = 8.0 * tr_z_xxx[i] - 4.0 * tr_z_xxxxx[i] * tke_0 - 8.0 * tr_zzz_xxx[i] * tbe_0 + 4.0 * tr_zzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xxxy[i] = 6.0 * tr_z_xxy[i] - 4.0 * tr_z_xxxxy[i] * tke_0 - 6.0 * tr_zzz_xxy[i] * tbe_0 + 4.0 * tr_zzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxxy[i] * tbe_0 + 4.0 * tr_xzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xxxz[i] = 6.0 * tr_z_xxz[i] - 4.0 * tr_z_xxxxz[i] * tke_0 - 6.0 * tr_zzz_xxz[i] * tbe_0 + 4.0 * tr_zzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxxz[i] * tbe_0 + 4.0 * tr_xzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xxyy[i] = 4.0 * tr_z_xyy[i] - 4.0 * tr_z_xxxyy[i] * tke_0 - 4.0 * tr_zzz_xyy[i] * tbe_0 + 4.0 * tr_zzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxyy[i] * tbe_0 + 4.0 * tr_xzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xxyz[i] = 4.0 * tr_z_xyz[i] - 4.0 * tr_z_xxxyz[i] * tke_0 - 4.0 * tr_zzz_xyz[i] * tbe_0 + 4.0 * tr_zzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxyz[i] * tbe_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xxzz[i] = 4.0 * tr_z_xzz[i] - 4.0 * tr_z_xxxzz[i] * tke_0 - 4.0 * tr_zzz_xzz[i] * tbe_0 + 4.0 * tr_zzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxzz[i] * tbe_0 + 4.0 * tr_xzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xyyy[i] = 2.0 * tr_z_yyy[i] - 4.0 * tr_z_xxyyy[i] * tke_0 - 2.0 * tr_zzz_yyy[i] * tbe_0 + 4.0 * tr_zzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xyyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xyyz[i] = 2.0 * tr_z_yyz[i] - 4.0 * tr_z_xxyyz[i] * tke_0 - 2.0 * tr_zzz_yyz[i] * tbe_0 + 4.0 * tr_zzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xyyz[i] * tbe_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xyzz[i] = 2.0 * tr_z_yzz[i] - 4.0 * tr_z_xxyzz[i] * tke_0 - 2.0 * tr_zzz_yzz[i] * tbe_0 + 4.0 * tr_zzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xyzz[i] * tbe_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_xzzz[i] = 2.0 * tr_z_zzz[i] - 4.0 * tr_z_xxzzz[i] * tke_0 - 2.0 * tr_zzz_zzz[i] * tbe_0 + 4.0 * tr_zzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xzzz[i] * tbe_0 + 4.0 * tr_xzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_yyyy[i] = -4.0 * tr_z_xyyyy[i] * tke_0 + 4.0 * tr_zzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_yyyz[i] = -4.0 * tr_z_xyyyz[i] * tke_0 + 4.0 * tr_zzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yyyz[i] * tbe_0 + 4.0 * tr_xzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_yyzz[i] = -4.0 * tr_z_xyyzz[i] * tke_0 + 4.0 * tr_zzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yyzz[i] * tbe_0 + 4.0 * tr_xzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_yzzz[i] = -4.0 * tr_z_xyzzz[i] * tke_0 + 4.0 * tr_zzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yzzz[i] * tbe_0 + 4.0 * tr_xzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zz_zzzz[i] = -4.0 * tr_z_xzzzz[i] * tke_0 + 4.0 * tr_zzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 630-645 components of targeted buffer : DG

    auto tr_z_0_y_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 630);

    auto tr_z_0_y_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 631);

    auto tr_z_0_y_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 632);

    auto tr_z_0_y_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 633);

    auto tr_z_0_y_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 634);

    auto tr_z_0_y_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 635);

    auto tr_z_0_y_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 636);

    auto tr_z_0_y_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 637);

    auto tr_z_0_y_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 638);

    auto tr_z_0_y_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 639);

    auto tr_z_0_y_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 640);

    auto tr_z_0_y_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 641);

    auto tr_z_0_y_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 642);

    auto tr_z_0_y_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 643);

    auto tr_z_0_y_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 644);

    #pragma omp simd aligned(tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_z_0_y_xx_xxxx, tr_z_0_y_xx_xxxy, tr_z_0_y_xx_xxxz, tr_z_0_y_xx_xxyy, tr_z_0_y_xx_xxyz, tr_z_0_y_xx_xxzz, tr_z_0_y_xx_xyyy, tr_z_0_y_xx_xyyz, tr_z_0_y_xx_xyzz, tr_z_0_y_xx_xzzz, tr_z_0_y_xx_yyyy, tr_z_0_y_xx_yyyz, tr_z_0_y_xx_yyzz, tr_z_0_y_xx_yzzz, tr_z_0_y_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xx_xxxx[i] = 4.0 * tr_xxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xxxy[i] = -2.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xxxz[i] = 4.0 * tr_xxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xxyy[i] = -4.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xxyz[i] = -2.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xxzz[i] = 4.0 * tr_xxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xyyy[i] = -6.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xyyz[i] = -4.0 * tr_xxz_xyz[i] * tbe_0 + 4.0 * tr_xxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xyzz[i] = -2.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_xzzz[i] = 4.0 * tr_xxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_yyyy[i] = -8.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_yyyz[i] = -6.0 * tr_xxz_yyz[i] * tbe_0 + 4.0 * tr_xxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_yyzz[i] = -4.0 * tr_xxz_yzz[i] * tbe_0 + 4.0 * tr_xxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_yzzz[i] = -2.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xx_zzzz[i] = 4.0 * tr_xxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 645-660 components of targeted buffer : DG

    auto tr_z_0_y_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 645);

    auto tr_z_0_y_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 646);

    auto tr_z_0_y_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 647);

    auto tr_z_0_y_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 648);

    auto tr_z_0_y_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 649);

    auto tr_z_0_y_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 650);

    auto tr_z_0_y_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 651);

    auto tr_z_0_y_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 652);

    auto tr_z_0_y_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 653);

    auto tr_z_0_y_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 654);

    auto tr_z_0_y_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 655);

    auto tr_z_0_y_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 656);

    auto tr_z_0_y_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 657);

    auto tr_z_0_y_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 658);

    auto tr_z_0_y_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 659);

    #pragma omp simd aligned(tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_z_0_y_xy_xxxx, tr_z_0_y_xy_xxxy, tr_z_0_y_xy_xxxz, tr_z_0_y_xy_xxyy, tr_z_0_y_xy_xxyz, tr_z_0_y_xy_xxzz, tr_z_0_y_xy_xyyy, tr_z_0_y_xy_xyyz, tr_z_0_y_xy_xyzz, tr_z_0_y_xy_xzzz, tr_z_0_y_xy_yyyy, tr_z_0_y_xy_yyyz, tr_z_0_y_xy_yyzz, tr_z_0_y_xy_yzzz, tr_z_0_y_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xy_xxxx[i] = -2.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xxxy[i] = -2.0 * tr_xz_xxxy[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xxxz[i] = -2.0 * tr_xz_xxxz[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xxyy[i] = -2.0 * tr_xz_xxyy[i] * tbe_0 - 4.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xxyz[i] = -2.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xxzz[i] = -2.0 * tr_xz_xxzz[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xyyy[i] = -2.0 * tr_xz_xyyy[i] * tbe_0 - 6.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xyyz[i] = -2.0 * tr_xz_xyyz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xyzz[i] = -2.0 * tr_xz_xyzz[i] * tbe_0 - 2.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_xzzz[i] = -2.0 * tr_xz_xzzz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_yyyy[i] = -2.0 * tr_xz_yyyy[i] * tbe_0 - 8.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_yyyz[i] = -2.0 * tr_xz_yyyz[i] * tbe_0 - 6.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_yyzz[i] = -2.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_yzzz[i] = -2.0 * tr_xz_yzzz[i] * tbe_0 - 2.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xy_zzzz[i] = -2.0 * tr_xz_zzzz[i] * tbe_0 + 4.0 * tr_xyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 660-675 components of targeted buffer : DG

    auto tr_z_0_y_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 660);

    auto tr_z_0_y_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 661);

    auto tr_z_0_y_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 662);

    auto tr_z_0_y_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 663);

    auto tr_z_0_y_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 664);

    auto tr_z_0_y_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 665);

    auto tr_z_0_y_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 666);

    auto tr_z_0_y_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 667);

    auto tr_z_0_y_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 668);

    auto tr_z_0_y_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 669);

    auto tr_z_0_y_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 670);

    auto tr_z_0_y_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 671);

    auto tr_z_0_y_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 672);

    auto tr_z_0_y_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 673);

    auto tr_z_0_y_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 674);

    #pragma omp simd aligned(tr_x_xxx, tr_x_xxxxy, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_yyy, tr_x_yyyyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_z_0_y_xz_xxxx, tr_z_0_y_xz_xxxy, tr_z_0_y_xz_xxxz, tr_z_0_y_xz_xxyy, tr_z_0_y_xz_xxyz, tr_z_0_y_xz_xxzz, tr_z_0_y_xz_xyyy, tr_z_0_y_xz_xyyz, tr_z_0_y_xz_xyzz, tr_z_0_y_xz_xzzz, tr_z_0_y_xz_yyyy, tr_z_0_y_xz_yyyz, tr_z_0_y_xz_yyzz, tr_z_0_y_xz_yzzz, tr_z_0_y_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xz_xxxx[i] = -2.0 * tr_x_xxxxy[i] * tke_0 + 4.0 * tr_xzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xxxy[i] = tr_x_xxx[i] - 2.0 * tr_x_xxxyy[i] * tke_0 - 2.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xxxz[i] = -2.0 * tr_x_xxxyz[i] * tke_0 + 4.0 * tr_xzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xxyy[i] = 2.0 * tr_x_xxy[i] - 2.0 * tr_x_xxyyy[i] * tke_0 - 4.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xxyz[i] = tr_x_xxz[i] - 2.0 * tr_x_xxyyz[i] * tke_0 - 2.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxyz[i] * tbe_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xxzz[i] = -2.0 * tr_x_xxyzz[i] * tke_0 + 4.0 * tr_xzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxzz[i] * tbe_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xyyy[i] = 3.0 * tr_x_xyy[i] - 2.0 * tr_x_xyyyy[i] * tke_0 - 6.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xyyz[i] = 2.0 * tr_x_xyz[i] - 2.0 * tr_x_xyyyz[i] * tke_0 - 4.0 * tr_xzz_xyz[i] * tbe_0 + 4.0 * tr_xzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xyzz[i] = tr_x_xzz[i] - 2.0 * tr_x_xyyzz[i] * tke_0 - 2.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_xzzz[i] = -2.0 * tr_x_xyzzz[i] * tke_0 + 4.0 * tr_xzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xzzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_yyyy[i] = 4.0 * tr_x_yyy[i] - 2.0 * tr_x_yyyyy[i] * tke_0 - 8.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_yyyz[i] = 3.0 * tr_x_yyz[i] - 2.0 * tr_x_yyyyz[i] * tke_0 - 6.0 * tr_xzz_yyz[i] * tbe_0 + 4.0 * tr_xzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyyz[i] * tbe_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_yyzz[i] = 2.0 * tr_x_yzz[i] - 2.0 * tr_x_yyyzz[i] * tke_0 - 4.0 * tr_xzz_yzz[i] * tbe_0 + 4.0 * tr_xzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyzz[i] * tbe_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_yzzz[i] = tr_x_zzz[i] - 2.0 * tr_x_yyzzz[i] * tke_0 - 2.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yzzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xz_zzzz[i] = -2.0 * tr_x_yzzzz[i] * tke_0 + 4.0 * tr_xzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_zzzz[i] * tbe_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 675-690 components of targeted buffer : DG

    auto tr_z_0_y_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 675);

    auto tr_z_0_y_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 676);

    auto tr_z_0_y_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 677);

    auto tr_z_0_y_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 678);

    auto tr_z_0_y_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 679);

    auto tr_z_0_y_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 680);

    auto tr_z_0_y_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 681);

    auto tr_z_0_y_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 682);

    auto tr_z_0_y_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 683);

    auto tr_z_0_y_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 684);

    auto tr_z_0_y_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 685);

    auto tr_z_0_y_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 686);

    auto tr_z_0_y_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 687);

    auto tr_z_0_y_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 688);

    auto tr_z_0_y_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 689);

    #pragma omp simd aligned(tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_z_0_y_yy_xxxx, tr_z_0_y_yy_xxxy, tr_z_0_y_yy_xxxz, tr_z_0_y_yy_xxyy, tr_z_0_y_yy_xxyz, tr_z_0_y_yy_xxzz, tr_z_0_y_yy_xyyy, tr_z_0_y_yy_xyyz, tr_z_0_y_yy_xyzz, tr_z_0_y_yy_xzzz, tr_z_0_y_yy_yyyy, tr_z_0_y_yy_yyyz, tr_z_0_y_yy_yyzz, tr_z_0_y_yy_yzzz, tr_z_0_y_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yy_xxxx[i] = -4.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_yyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xxxy[i] = -4.0 * tr_yz_xxxy[i] * tbe_0 - 2.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xxxz[i] = -4.0 * tr_yz_xxxz[i] * tbe_0 + 4.0 * tr_yyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xxyy[i] = -4.0 * tr_yz_xxyy[i] * tbe_0 - 4.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xxyz[i] = -4.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_yyz_xxz[i] * tbe_0 + 4.0 * tr_yyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xxzz[i] = -4.0 * tr_yz_xxzz[i] * tbe_0 + 4.0 * tr_yyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xyyy[i] = -4.0 * tr_yz_xyyy[i] * tbe_0 - 6.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xyyz[i] = -4.0 * tr_yz_xyyz[i] * tbe_0 - 4.0 * tr_yyz_xyz[i] * tbe_0 + 4.0 * tr_yyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xyzz[i] = -4.0 * tr_yz_xyzz[i] * tbe_0 - 2.0 * tr_yyz_xzz[i] * tbe_0 + 4.0 * tr_yyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_xzzz[i] = -4.0 * tr_yz_xzzz[i] * tbe_0 + 4.0 * tr_yyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_yyyy[i] = -4.0 * tr_yz_yyyy[i] * tbe_0 - 8.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_yyyz[i] = -4.0 * tr_yz_yyyz[i] * tbe_0 - 6.0 * tr_yyz_yyz[i] * tbe_0 + 4.0 * tr_yyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_yyzz[i] = -4.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_yyz_yzz[i] * tbe_0 + 4.0 * tr_yyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_yzzz[i] = -4.0 * tr_yz_yzzz[i] * tbe_0 - 2.0 * tr_yyz_zzz[i] * tbe_0 + 4.0 * tr_yyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yy_zzzz[i] = -4.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_yyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 690-705 components of targeted buffer : DG

    auto tr_z_0_y_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 690);

    auto tr_z_0_y_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 691);

    auto tr_z_0_y_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 692);

    auto tr_z_0_y_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 693);

    auto tr_z_0_y_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 694);

    auto tr_z_0_y_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 695);

    auto tr_z_0_y_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 696);

    auto tr_z_0_y_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 697);

    auto tr_z_0_y_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 698);

    auto tr_z_0_y_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 699);

    auto tr_z_0_y_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 700);

    auto tr_z_0_y_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 701);

    auto tr_z_0_y_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 702);

    auto tr_z_0_y_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 703);

    auto tr_z_0_y_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 704);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_y_xxx, tr_y_xxxxy, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_yyy, tr_y_yyyyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_z_0_y_yz_xxxx, tr_z_0_y_yz_xxxy, tr_z_0_y_yz_xxxz, tr_z_0_y_yz_xxyy, tr_z_0_y_yz_xxyz, tr_z_0_y_yz_xxzz, tr_z_0_y_yz_xyyy, tr_z_0_y_yz_xyyz, tr_z_0_y_yz_xyzz, tr_z_0_y_yz_xzzz, tr_z_0_y_yz_yyyy, tr_z_0_y_yz_yyyz, tr_z_0_y_yz_yyzz, tr_z_0_y_yz_yzzz, tr_z_0_y_yz_zzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yz_xxxx[i] = tr_0_xxxx[i] - 2.0 * tr_zz_xxxx[i] * tbe_0 - 2.0 * tr_y_xxxxy[i] * tke_0 + 4.0 * tr_yzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xxxy[i] = tr_0_xxxy[i] - 2.0 * tr_zz_xxxy[i] * tbe_0 + tr_y_xxx[i] - 2.0 * tr_y_xxxyy[i] * tke_0 - 2.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxxy[i] * tbe_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xxxz[i] = tr_0_xxxz[i] - 2.0 * tr_zz_xxxz[i] * tbe_0 - 2.0 * tr_y_xxxyz[i] * tke_0 + 4.0 * tr_yzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxxz[i] * tbe_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xxyy[i] = tr_0_xxyy[i] - 2.0 * tr_zz_xxyy[i] * tbe_0 + 2.0 * tr_y_xxy[i] - 2.0 * tr_y_xxyyy[i] * tke_0 - 4.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxyy[i] * tbe_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xxyz[i] = tr_0_xxyz[i] - 2.0 * tr_zz_xxyz[i] * tbe_0 + tr_y_xxz[i] - 2.0 * tr_y_xxyyz[i] * tke_0 - 2.0 * tr_yzz_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxyz[i] * tbe_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xxzz[i] = tr_0_xxzz[i] - 2.0 * tr_zz_xxzz[i] * tbe_0 - 2.0 * tr_y_xxyzz[i] * tke_0 + 4.0 * tr_yzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxzz[i] * tbe_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xyyy[i] = tr_0_xyyy[i] - 2.0 * tr_zz_xyyy[i] * tbe_0 + 3.0 * tr_y_xyy[i] - 2.0 * tr_y_xyyyy[i] * tke_0 - 6.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xyyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xyyz[i] = tr_0_xyyz[i] - 2.0 * tr_zz_xyyz[i] * tbe_0 + 2.0 * tr_y_xyz[i] - 2.0 * tr_y_xyyyz[i] * tke_0 - 4.0 * tr_yzz_xyz[i] * tbe_0 + 4.0 * tr_yzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xyyz[i] * tbe_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xyzz[i] = tr_0_xyzz[i] - 2.0 * tr_zz_xyzz[i] * tbe_0 + tr_y_xzz[i] - 2.0 * tr_y_xyyzz[i] * tke_0 - 2.0 * tr_yzz_xzz[i] * tbe_0 + 4.0 * tr_yzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xyzz[i] * tbe_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_xzzz[i] = tr_0_xzzz[i] - 2.0 * tr_zz_xzzz[i] * tbe_0 - 2.0 * tr_y_xyzzz[i] * tke_0 + 4.0 * tr_yzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xzzz[i] * tbe_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_yyyy[i] = tr_0_yyyy[i] - 2.0 * tr_zz_yyyy[i] * tbe_0 + 4.0 * tr_y_yyy[i] - 2.0 * tr_y_yyyyy[i] * tke_0 - 8.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yyyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_yyyz[i] = tr_0_yyyz[i] - 2.0 * tr_zz_yyyz[i] * tbe_0 + 3.0 * tr_y_yyz[i] - 2.0 * tr_y_yyyyz[i] * tke_0 - 6.0 * tr_yzz_yyz[i] * tbe_0 + 4.0 * tr_yzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yyyz[i] * tbe_0 + 4.0 * tr_yyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_yyzz[i] = tr_0_yyzz[i] - 2.0 * tr_zz_yyzz[i] * tbe_0 + 2.0 * tr_y_yzz[i] - 2.0 * tr_y_yyyzz[i] * tke_0 - 4.0 * tr_yzz_yzz[i] * tbe_0 + 4.0 * tr_yzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yyzz[i] * tbe_0 + 4.0 * tr_yyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_yzzz[i] = tr_0_yzzz[i] - 2.0 * tr_zz_yzzz[i] * tbe_0 + tr_y_zzz[i] - 2.0 * tr_y_yyzzz[i] * tke_0 - 2.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_yzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yzzz[i] * tbe_0 + 4.0 * tr_yyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yz_zzzz[i] = tr_0_zzzz[i] - 2.0 * tr_zz_zzzz[i] * tbe_0 - 2.0 * tr_y_yzzzz[i] * tke_0 + 4.0 * tr_yzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_zzzz[i] * tbe_0 + 4.0 * tr_yyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 705-720 components of targeted buffer : DG

    auto tr_z_0_y_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 705);

    auto tr_z_0_y_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 706);

    auto tr_z_0_y_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 707);

    auto tr_z_0_y_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 708);

    auto tr_z_0_y_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 709);

    auto tr_z_0_y_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 710);

    auto tr_z_0_y_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 711);

    auto tr_z_0_y_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 712);

    auto tr_z_0_y_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 713);

    auto tr_z_0_y_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 714);

    auto tr_z_0_y_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 715);

    auto tr_z_0_y_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 716);

    auto tr_z_0_y_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 717);

    auto tr_z_0_y_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 718);

    auto tr_z_0_y_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 719);

    #pragma omp simd aligned(tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, tr_z_0_y_zz_xxxx, tr_z_0_y_zz_xxxy, tr_z_0_y_zz_xxxz, tr_z_0_y_zz_xxyy, tr_z_0_y_zz_xxyz, tr_z_0_y_zz_xxzz, tr_z_0_y_zz_xyyy, tr_z_0_y_zz_xyyz, tr_z_0_y_zz_xyzz, tr_z_0_y_zz_xzzz, tr_z_0_y_zz_yyyy, tr_z_0_y_zz_yyyz, tr_z_0_y_zz_yyzz, tr_z_0_y_zz_yzzz, tr_z_0_y_zz_zzzz, tr_z_xxx, tr_z_xxxxy, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_yyy, tr_z_yyyyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, tr_zzz_xxx, tr_zzz_xxxxy, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyyyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zz_xxxx[i] = -4.0 * tr_z_xxxxy[i] * tke_0 + 4.0 * tr_zzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xxxy[i] = 2.0 * tr_z_xxx[i] - 4.0 * tr_z_xxxyy[i] * tke_0 - 2.0 * tr_zzz_xxx[i] * tbe_0 + 4.0 * tr_zzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxxy[i] * tbe_0 + 4.0 * tr_yzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xxxz[i] = -4.0 * tr_z_xxxyz[i] * tke_0 + 4.0 * tr_zzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxxz[i] * tbe_0 + 4.0 * tr_yzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xxyy[i] = 4.0 * tr_z_xxy[i] - 4.0 * tr_z_xxyyy[i] * tke_0 - 4.0 * tr_zzz_xxy[i] * tbe_0 + 4.0 * tr_zzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxyy[i] * tbe_0 + 4.0 * tr_yzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xxyz[i] = 2.0 * tr_z_xxz[i] - 4.0 * tr_z_xxyyz[i] * tke_0 - 2.0 * tr_zzz_xxz[i] * tbe_0 + 4.0 * tr_zzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxyz[i] * tbe_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xxzz[i] = -4.0 * tr_z_xxyzz[i] * tke_0 + 4.0 * tr_zzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxzz[i] * tbe_0 + 4.0 * tr_yzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xyyy[i] = 6.0 * tr_z_xyy[i] - 4.0 * tr_z_xyyyy[i] * tke_0 - 6.0 * tr_zzz_xyy[i] * tbe_0 + 4.0 * tr_zzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xyyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xyyz[i] = 4.0 * tr_z_xyz[i] - 4.0 * tr_z_xyyyz[i] * tke_0 - 4.0 * tr_zzz_xyz[i] * tbe_0 + 4.0 * tr_zzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xyyz[i] * tbe_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xyzz[i] = 2.0 * tr_z_xzz[i] - 4.0 * tr_z_xyyzz[i] * tke_0 - 2.0 * tr_zzz_xzz[i] * tbe_0 + 4.0 * tr_zzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xyzz[i] * tbe_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_xzzz[i] = -4.0 * tr_z_xyzzz[i] * tke_0 + 4.0 * tr_zzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xzzz[i] * tbe_0 + 4.0 * tr_yzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_yyyy[i] = 8.0 * tr_z_yyy[i] - 4.0 * tr_z_yyyyy[i] * tke_0 - 8.0 * tr_zzz_yyy[i] * tbe_0 + 4.0 * tr_zzz_yyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_yyyz[i] = 6.0 * tr_z_yyz[i] - 4.0 * tr_z_yyyyz[i] * tke_0 - 6.0 * tr_zzz_yyz[i] * tbe_0 + 4.0 * tr_zzz_yyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yyyz[i] * tbe_0 + 4.0 * tr_yzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_yyzz[i] = 4.0 * tr_z_yzz[i] - 4.0 * tr_z_yyyzz[i] * tke_0 - 4.0 * tr_zzz_yzz[i] * tbe_0 + 4.0 * tr_zzz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yyzz[i] * tbe_0 + 4.0 * tr_yzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_yzzz[i] = 2.0 * tr_z_zzz[i] - 4.0 * tr_z_yyzzz[i] * tke_0 - 2.0 * tr_zzz_zzz[i] * tbe_0 + 4.0 * tr_zzz_yyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yzzz[i] * tbe_0 + 4.0 * tr_yzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zz_zzzz[i] = -4.0 * tr_z_yzzzz[i] * tke_0 + 4.0 * tr_zzz_yzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_zzzz[i] * tbe_0 + 4.0 * tr_yzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 720-735 components of targeted buffer : DG

    auto tr_z_0_z_xx_xxxx = pbuffer.data(idx_op_geom_110_dg + 720);

    auto tr_z_0_z_xx_xxxy = pbuffer.data(idx_op_geom_110_dg + 721);

    auto tr_z_0_z_xx_xxxz = pbuffer.data(idx_op_geom_110_dg + 722);

    auto tr_z_0_z_xx_xxyy = pbuffer.data(idx_op_geom_110_dg + 723);

    auto tr_z_0_z_xx_xxyz = pbuffer.data(idx_op_geom_110_dg + 724);

    auto tr_z_0_z_xx_xxzz = pbuffer.data(idx_op_geom_110_dg + 725);

    auto tr_z_0_z_xx_xyyy = pbuffer.data(idx_op_geom_110_dg + 726);

    auto tr_z_0_z_xx_xyyz = pbuffer.data(idx_op_geom_110_dg + 727);

    auto tr_z_0_z_xx_xyzz = pbuffer.data(idx_op_geom_110_dg + 728);

    auto tr_z_0_z_xx_xzzz = pbuffer.data(idx_op_geom_110_dg + 729);

    auto tr_z_0_z_xx_yyyy = pbuffer.data(idx_op_geom_110_dg + 730);

    auto tr_z_0_z_xx_yyyz = pbuffer.data(idx_op_geom_110_dg + 731);

    auto tr_z_0_z_xx_yyzz = pbuffer.data(idx_op_geom_110_dg + 732);

    auto tr_z_0_z_xx_yzzz = pbuffer.data(idx_op_geom_110_dg + 733);

    auto tr_z_0_z_xx_zzzz = pbuffer.data(idx_op_geom_110_dg + 734);

    #pragma omp simd aligned(tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_z_0_z_xx_xxxx, tr_z_0_z_xx_xxxy, tr_z_0_z_xx_xxxz, tr_z_0_z_xx_xxyy, tr_z_0_z_xx_xxyz, tr_z_0_z_xx_xxzz, tr_z_0_z_xx_xyyy, tr_z_0_z_xx_xyyz, tr_z_0_z_xx_xyzz, tr_z_0_z_xx_xzzz, tr_z_0_z_xx_yyyy, tr_z_0_z_xx_yyyz, tr_z_0_z_xx_yyzz, tr_z_0_z_xx_yzzz, tr_z_0_z_xx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xx_xxxx[i] = -2.0 * tr_xx_xxxx[i] * tbe_0 + 4.0 * tr_xxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xxxy[i] = -2.0 * tr_xx_xxxy[i] * tbe_0 + 4.0 * tr_xxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xxxz[i] = -2.0 * tr_xx_xxxz[i] * tbe_0 - 2.0 * tr_xxz_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xxyy[i] = -2.0 * tr_xx_xxyy[i] * tbe_0 + 4.0 * tr_xxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xxyz[i] = -2.0 * tr_xx_xxyz[i] * tbe_0 - 2.0 * tr_xxz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xxzz[i] = -2.0 * tr_xx_xxzz[i] * tbe_0 - 4.0 * tr_xxz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xyyy[i] = -2.0 * tr_xx_xyyy[i] * tbe_0 + 4.0 * tr_xxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xyyz[i] = -2.0 * tr_xx_xyyz[i] * tbe_0 - 2.0 * tr_xxz_xyy[i] * tbe_0 + 4.0 * tr_xxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xyzz[i] = -2.0 * tr_xx_xyzz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tbe_0 + 4.0 * tr_xxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_xzzz[i] = -2.0 * tr_xx_xzzz[i] * tbe_0 - 6.0 * tr_xxz_xzz[i] * tbe_0 + 4.0 * tr_xxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_yyyy[i] = -2.0 * tr_xx_yyyy[i] * tbe_0 + 4.0 * tr_xxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_yyyz[i] = -2.0 * tr_xx_yyyz[i] * tbe_0 - 2.0 * tr_xxz_yyy[i] * tbe_0 + 4.0 * tr_xxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_yyzz[i] = -2.0 * tr_xx_yyzz[i] * tbe_0 - 4.0 * tr_xxz_yyz[i] * tbe_0 + 4.0 * tr_xxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_yzzz[i] = -2.0 * tr_xx_yzzz[i] * tbe_0 - 6.0 * tr_xxz_yzz[i] * tbe_0 + 4.0 * tr_xxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xx_zzzz[i] = -2.0 * tr_xx_zzzz[i] * tbe_0 - 8.0 * tr_xxz_zzz[i] * tbe_0 + 4.0 * tr_xxz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 735-750 components of targeted buffer : DG

    auto tr_z_0_z_xy_xxxx = pbuffer.data(idx_op_geom_110_dg + 735);

    auto tr_z_0_z_xy_xxxy = pbuffer.data(idx_op_geom_110_dg + 736);

    auto tr_z_0_z_xy_xxxz = pbuffer.data(idx_op_geom_110_dg + 737);

    auto tr_z_0_z_xy_xxyy = pbuffer.data(idx_op_geom_110_dg + 738);

    auto tr_z_0_z_xy_xxyz = pbuffer.data(idx_op_geom_110_dg + 739);

    auto tr_z_0_z_xy_xxzz = pbuffer.data(idx_op_geom_110_dg + 740);

    auto tr_z_0_z_xy_xyyy = pbuffer.data(idx_op_geom_110_dg + 741);

    auto tr_z_0_z_xy_xyyz = pbuffer.data(idx_op_geom_110_dg + 742);

    auto tr_z_0_z_xy_xyzz = pbuffer.data(idx_op_geom_110_dg + 743);

    auto tr_z_0_z_xy_xzzz = pbuffer.data(idx_op_geom_110_dg + 744);

    auto tr_z_0_z_xy_yyyy = pbuffer.data(idx_op_geom_110_dg + 745);

    auto tr_z_0_z_xy_yyyz = pbuffer.data(idx_op_geom_110_dg + 746);

    auto tr_z_0_z_xy_yyzz = pbuffer.data(idx_op_geom_110_dg + 747);

    auto tr_z_0_z_xy_yzzz = pbuffer.data(idx_op_geom_110_dg + 748);

    auto tr_z_0_z_xy_zzzz = pbuffer.data(idx_op_geom_110_dg + 749);

    #pragma omp simd aligned(tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_z_0_z_xy_xxxx, tr_z_0_z_xy_xxxy, tr_z_0_z_xy_xxxz, tr_z_0_z_xy_xxyy, tr_z_0_z_xy_xxyz, tr_z_0_z_xy_xxzz, tr_z_0_z_xy_xyyy, tr_z_0_z_xy_xyyz, tr_z_0_z_xy_xyzz, tr_z_0_z_xy_xzzz, tr_z_0_z_xy_yyyy, tr_z_0_z_xy_yyyz, tr_z_0_z_xy_yyzz, tr_z_0_z_xy_yzzz, tr_z_0_z_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xy_xxxx[i] = -2.0 * tr_xy_xxxx[i] * tbe_0 + 4.0 * tr_xyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xxxy[i] = -2.0 * tr_xy_xxxy[i] * tbe_0 + 4.0 * tr_xyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xxxz[i] = -2.0 * tr_xy_xxxz[i] * tbe_0 - 2.0 * tr_xyz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xxyy[i] = -2.0 * tr_xy_xxyy[i] * tbe_0 + 4.0 * tr_xyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xxyz[i] = -2.0 * tr_xy_xxyz[i] * tbe_0 - 2.0 * tr_xyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xxzz[i] = -2.0 * tr_xy_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xyyy[i] = -2.0 * tr_xy_xyyy[i] * tbe_0 + 4.0 * tr_xyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xyyz[i] = -2.0 * tr_xy_xyyz[i] * tbe_0 - 2.0 * tr_xyz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xyzz[i] = -2.0 * tr_xy_xyzz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tbe_0 + 4.0 * tr_xyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_xzzz[i] = -2.0 * tr_xy_xzzz[i] * tbe_0 - 6.0 * tr_xyz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_yyyy[i] = -2.0 * tr_xy_yyyy[i] * tbe_0 + 4.0 * tr_xyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_yyyz[i] = -2.0 * tr_xy_yyyz[i] * tbe_0 - 2.0 * tr_xyz_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_yyzz[i] = -2.0 * tr_xy_yyzz[i] * tbe_0 - 4.0 * tr_xyz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_yzzz[i] = -2.0 * tr_xy_yzzz[i] * tbe_0 - 6.0 * tr_xyz_yzz[i] * tbe_0 + 4.0 * tr_xyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xy_zzzz[i] = -2.0 * tr_xy_zzzz[i] * tbe_0 - 8.0 * tr_xyz_zzz[i] * tbe_0 + 4.0 * tr_xyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 750-765 components of targeted buffer : DG

    auto tr_z_0_z_xz_xxxx = pbuffer.data(idx_op_geom_110_dg + 750);

    auto tr_z_0_z_xz_xxxy = pbuffer.data(idx_op_geom_110_dg + 751);

    auto tr_z_0_z_xz_xxxz = pbuffer.data(idx_op_geom_110_dg + 752);

    auto tr_z_0_z_xz_xxyy = pbuffer.data(idx_op_geom_110_dg + 753);

    auto tr_z_0_z_xz_xxyz = pbuffer.data(idx_op_geom_110_dg + 754);

    auto tr_z_0_z_xz_xxzz = pbuffer.data(idx_op_geom_110_dg + 755);

    auto tr_z_0_z_xz_xyyy = pbuffer.data(idx_op_geom_110_dg + 756);

    auto tr_z_0_z_xz_xyyz = pbuffer.data(idx_op_geom_110_dg + 757);

    auto tr_z_0_z_xz_xyzz = pbuffer.data(idx_op_geom_110_dg + 758);

    auto tr_z_0_z_xz_xzzz = pbuffer.data(idx_op_geom_110_dg + 759);

    auto tr_z_0_z_xz_yyyy = pbuffer.data(idx_op_geom_110_dg + 760);

    auto tr_z_0_z_xz_yyyz = pbuffer.data(idx_op_geom_110_dg + 761);

    auto tr_z_0_z_xz_yyzz = pbuffer.data(idx_op_geom_110_dg + 762);

    auto tr_z_0_z_xz_yzzz = pbuffer.data(idx_op_geom_110_dg + 763);

    auto tr_z_0_z_xz_zzzz = pbuffer.data(idx_op_geom_110_dg + 764);

    #pragma omp simd aligned(tr_x_xxx, tr_x_xxxxz, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_x_zzzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, tr_z_0_z_xz_xxxx, tr_z_0_z_xz_xxxy, tr_z_0_z_xz_xxxz, tr_z_0_z_xz_xxyy, tr_z_0_z_xz_xxyz, tr_z_0_z_xz_xxzz, tr_z_0_z_xz_xyyy, tr_z_0_z_xz_xyyz, tr_z_0_z_xz_xyzz, tr_z_0_z_xz_xzzz, tr_z_0_z_xz_yyyy, tr_z_0_z_xz_yyyz, tr_z_0_z_xz_yyzz, tr_z_0_z_xz_yzzz, tr_z_0_z_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xz_xxxx[i] = -2.0 * tr_x_xxxxz[i] * tke_0 - 6.0 * tr_xz_xxxx[i] * tbe_0 + 4.0 * tr_xzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xxxy[i] = -2.0 * tr_x_xxxyz[i] * tke_0 - 6.0 * tr_xz_xxxy[i] * tbe_0 + 4.0 * tr_xzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xxxz[i] = tr_x_xxx[i] - 2.0 * tr_x_xxxzz[i] * tke_0 - 6.0 * tr_xz_xxxz[i] * tbe_0 - 2.0 * tr_xzz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xxyy[i] = -2.0 * tr_x_xxyyz[i] * tke_0 - 6.0 * tr_xz_xxyy[i] * tbe_0 + 4.0 * tr_xzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xxyz[i] = tr_x_xxy[i] - 2.0 * tr_x_xxyzz[i] * tke_0 - 6.0 * tr_xz_xxyz[i] * tbe_0 - 2.0 * tr_xzz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xxzz[i] = 2.0 * tr_x_xxz[i] - 2.0 * tr_x_xxzzz[i] * tke_0 - 6.0 * tr_xz_xxzz[i] * tbe_0 - 4.0 * tr_xzz_xxz[i] * tbe_0 + 4.0 * tr_xzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xyyy[i] = -2.0 * tr_x_xyyyz[i] * tke_0 - 6.0 * tr_xz_xyyy[i] * tbe_0 + 4.0 * tr_xzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xyyz[i] = tr_x_xyy[i] - 2.0 * tr_x_xyyzz[i] * tke_0 - 6.0 * tr_xz_xyyz[i] * tbe_0 - 2.0 * tr_xzz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xyzz[i] = 2.0 * tr_x_xyz[i] - 2.0 * tr_x_xyzzz[i] * tke_0 - 6.0 * tr_xz_xyzz[i] * tbe_0 - 4.0 * tr_xzz_xyz[i] * tbe_0 + 4.0 * tr_xzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_xzzz[i] = 3.0 * tr_x_xzz[i] - 2.0 * tr_x_xzzzz[i] * tke_0 - 6.0 * tr_xz_xzzz[i] * tbe_0 - 6.0 * tr_xzz_xzz[i] * tbe_0 + 4.0 * tr_xzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_yyyy[i] = -2.0 * tr_x_yyyyz[i] * tke_0 - 6.0 * tr_xz_yyyy[i] * tbe_0 + 4.0 * tr_xzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_yyyz[i] = tr_x_yyy[i] - 2.0 * tr_x_yyyzz[i] * tke_0 - 6.0 * tr_xz_yyyz[i] * tbe_0 - 2.0 * tr_xzz_yyy[i] * tbe_0 + 4.0 * tr_xzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_yyzz[i] = 2.0 * tr_x_yyz[i] - 2.0 * tr_x_yyzzz[i] * tke_0 - 6.0 * tr_xz_yyzz[i] * tbe_0 - 4.0 * tr_xzz_yyz[i] * tbe_0 + 4.0 * tr_xzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_yzzz[i] = 3.0 * tr_x_yzz[i] - 2.0 * tr_x_yzzzz[i] * tke_0 - 6.0 * tr_xz_yzzz[i] * tbe_0 - 6.0 * tr_xzz_yzz[i] * tbe_0 + 4.0 * tr_xzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xz_zzzz[i] = 4.0 * tr_x_zzz[i] - 2.0 * tr_x_zzzzz[i] * tke_0 - 6.0 * tr_xz_zzzz[i] * tbe_0 - 8.0 * tr_xzz_zzz[i] * tbe_0 + 4.0 * tr_xzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 765-780 components of targeted buffer : DG

    auto tr_z_0_z_yy_xxxx = pbuffer.data(idx_op_geom_110_dg + 765);

    auto tr_z_0_z_yy_xxxy = pbuffer.data(idx_op_geom_110_dg + 766);

    auto tr_z_0_z_yy_xxxz = pbuffer.data(idx_op_geom_110_dg + 767);

    auto tr_z_0_z_yy_xxyy = pbuffer.data(idx_op_geom_110_dg + 768);

    auto tr_z_0_z_yy_xxyz = pbuffer.data(idx_op_geom_110_dg + 769);

    auto tr_z_0_z_yy_xxzz = pbuffer.data(idx_op_geom_110_dg + 770);

    auto tr_z_0_z_yy_xyyy = pbuffer.data(idx_op_geom_110_dg + 771);

    auto tr_z_0_z_yy_xyyz = pbuffer.data(idx_op_geom_110_dg + 772);

    auto tr_z_0_z_yy_xyzz = pbuffer.data(idx_op_geom_110_dg + 773);

    auto tr_z_0_z_yy_xzzz = pbuffer.data(idx_op_geom_110_dg + 774);

    auto tr_z_0_z_yy_yyyy = pbuffer.data(idx_op_geom_110_dg + 775);

    auto tr_z_0_z_yy_yyyz = pbuffer.data(idx_op_geom_110_dg + 776);

    auto tr_z_0_z_yy_yyzz = pbuffer.data(idx_op_geom_110_dg + 777);

    auto tr_z_0_z_yy_yzzz = pbuffer.data(idx_op_geom_110_dg + 778);

    auto tr_z_0_z_yy_zzzz = pbuffer.data(idx_op_geom_110_dg + 779);

    #pragma omp simd aligned(tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_z_0_z_yy_xxxx, tr_z_0_z_yy_xxxy, tr_z_0_z_yy_xxxz, tr_z_0_z_yy_xxyy, tr_z_0_z_yy_xxyz, tr_z_0_z_yy_xxzz, tr_z_0_z_yy_xyyy, tr_z_0_z_yy_xyyz, tr_z_0_z_yy_xyzz, tr_z_0_z_yy_xzzz, tr_z_0_z_yy_yyyy, tr_z_0_z_yy_yyyz, tr_z_0_z_yy_yyzz, tr_z_0_z_yy_yzzz, tr_z_0_z_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yy_xxxx[i] = -2.0 * tr_yy_xxxx[i] * tbe_0 + 4.0 * tr_yyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xxxy[i] = -2.0 * tr_yy_xxxy[i] * tbe_0 + 4.0 * tr_yyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xxxz[i] = -2.0 * tr_yy_xxxz[i] * tbe_0 - 2.0 * tr_yyz_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xxyy[i] = -2.0 * tr_yy_xxyy[i] * tbe_0 + 4.0 * tr_yyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xxyz[i] = -2.0 * tr_yy_xxyz[i] * tbe_0 - 2.0 * tr_yyz_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xxzz[i] = -2.0 * tr_yy_xxzz[i] * tbe_0 - 4.0 * tr_yyz_xxz[i] * tbe_0 + 4.0 * tr_yyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xyyy[i] = -2.0 * tr_yy_xyyy[i] * tbe_0 + 4.0 * tr_yyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xyyz[i] = -2.0 * tr_yy_xyyz[i] * tbe_0 - 2.0 * tr_yyz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xyzz[i] = -2.0 * tr_yy_xyzz[i] * tbe_0 - 4.0 * tr_yyz_xyz[i] * tbe_0 + 4.0 * tr_yyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_xzzz[i] = -2.0 * tr_yy_xzzz[i] * tbe_0 - 6.0 * tr_yyz_xzz[i] * tbe_0 + 4.0 * tr_yyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_yyyy[i] = -2.0 * tr_yy_yyyy[i] * tbe_0 + 4.0 * tr_yyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_yyyz[i] = -2.0 * tr_yy_yyyz[i] * tbe_0 - 2.0 * tr_yyz_yyy[i] * tbe_0 + 4.0 * tr_yyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_yyzz[i] = -2.0 * tr_yy_yyzz[i] * tbe_0 - 4.0 * tr_yyz_yyz[i] * tbe_0 + 4.0 * tr_yyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_yzzz[i] = -2.0 * tr_yy_yzzz[i] * tbe_0 - 6.0 * tr_yyz_yzz[i] * tbe_0 + 4.0 * tr_yyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yy_zzzz[i] = -2.0 * tr_yy_zzzz[i] * tbe_0 - 8.0 * tr_yyz_zzz[i] * tbe_0 + 4.0 * tr_yyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 780-795 components of targeted buffer : DG

    auto tr_z_0_z_yz_xxxx = pbuffer.data(idx_op_geom_110_dg + 780);

    auto tr_z_0_z_yz_xxxy = pbuffer.data(idx_op_geom_110_dg + 781);

    auto tr_z_0_z_yz_xxxz = pbuffer.data(idx_op_geom_110_dg + 782);

    auto tr_z_0_z_yz_xxyy = pbuffer.data(idx_op_geom_110_dg + 783);

    auto tr_z_0_z_yz_xxyz = pbuffer.data(idx_op_geom_110_dg + 784);

    auto tr_z_0_z_yz_xxzz = pbuffer.data(idx_op_geom_110_dg + 785);

    auto tr_z_0_z_yz_xyyy = pbuffer.data(idx_op_geom_110_dg + 786);

    auto tr_z_0_z_yz_xyyz = pbuffer.data(idx_op_geom_110_dg + 787);

    auto tr_z_0_z_yz_xyzz = pbuffer.data(idx_op_geom_110_dg + 788);

    auto tr_z_0_z_yz_xzzz = pbuffer.data(idx_op_geom_110_dg + 789);

    auto tr_z_0_z_yz_yyyy = pbuffer.data(idx_op_geom_110_dg + 790);

    auto tr_z_0_z_yz_yyyz = pbuffer.data(idx_op_geom_110_dg + 791);

    auto tr_z_0_z_yz_yyzz = pbuffer.data(idx_op_geom_110_dg + 792);

    auto tr_z_0_z_yz_yzzz = pbuffer.data(idx_op_geom_110_dg + 793);

    auto tr_z_0_z_yz_zzzz = pbuffer.data(idx_op_geom_110_dg + 794);

    #pragma omp simd aligned(tr_y_xxx, tr_y_xxxxz, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_y_zzzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, tr_z_0_z_yz_xxxx, tr_z_0_z_yz_xxxy, tr_z_0_z_yz_xxxz, tr_z_0_z_yz_xxyy, tr_z_0_z_yz_xxyz, tr_z_0_z_yz_xxzz, tr_z_0_z_yz_xyyy, tr_z_0_z_yz_xyyz, tr_z_0_z_yz_xyzz, tr_z_0_z_yz_xzzz, tr_z_0_z_yz_yyyy, tr_z_0_z_yz_yyyz, tr_z_0_z_yz_yyzz, tr_z_0_z_yz_yzzz, tr_z_0_z_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yz_xxxx[i] = -2.0 * tr_y_xxxxz[i] * tke_0 - 6.0 * tr_yz_xxxx[i] * tbe_0 + 4.0 * tr_yzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xxxy[i] = -2.0 * tr_y_xxxyz[i] * tke_0 - 6.0 * tr_yz_xxxy[i] * tbe_0 + 4.0 * tr_yzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xxxz[i] = tr_y_xxx[i] - 2.0 * tr_y_xxxzz[i] * tke_0 - 6.0 * tr_yz_xxxz[i] * tbe_0 - 2.0 * tr_yzz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xxyy[i] = -2.0 * tr_y_xxyyz[i] * tke_0 - 6.0 * tr_yz_xxyy[i] * tbe_0 + 4.0 * tr_yzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xxyz[i] = tr_y_xxy[i] - 2.0 * tr_y_xxyzz[i] * tke_0 - 6.0 * tr_yz_xxyz[i] * tbe_0 - 2.0 * tr_yzz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xxzz[i] = 2.0 * tr_y_xxz[i] - 2.0 * tr_y_xxzzz[i] * tke_0 - 6.0 * tr_yz_xxzz[i] * tbe_0 - 4.0 * tr_yzz_xxz[i] * tbe_0 + 4.0 * tr_yzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xyyy[i] = -2.0 * tr_y_xyyyz[i] * tke_0 - 6.0 * tr_yz_xyyy[i] * tbe_0 + 4.0 * tr_yzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xyyz[i] = tr_y_xyy[i] - 2.0 * tr_y_xyyzz[i] * tke_0 - 6.0 * tr_yz_xyyz[i] * tbe_0 - 2.0 * tr_yzz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xyzz[i] = 2.0 * tr_y_xyz[i] - 2.0 * tr_y_xyzzz[i] * tke_0 - 6.0 * tr_yz_xyzz[i] * tbe_0 - 4.0 * tr_yzz_xyz[i] * tbe_0 + 4.0 * tr_yzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_xzzz[i] = 3.0 * tr_y_xzz[i] - 2.0 * tr_y_xzzzz[i] * tke_0 - 6.0 * tr_yz_xzzz[i] * tbe_0 - 6.0 * tr_yzz_xzz[i] * tbe_0 + 4.0 * tr_yzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_yyyy[i] = -2.0 * tr_y_yyyyz[i] * tke_0 - 6.0 * tr_yz_yyyy[i] * tbe_0 + 4.0 * tr_yzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_yyyz[i] = tr_y_yyy[i] - 2.0 * tr_y_yyyzz[i] * tke_0 - 6.0 * tr_yz_yyyz[i] * tbe_0 - 2.0 * tr_yzz_yyy[i] * tbe_0 + 4.0 * tr_yzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_yyzz[i] = 2.0 * tr_y_yyz[i] - 2.0 * tr_y_yyzzz[i] * tke_0 - 6.0 * tr_yz_yyzz[i] * tbe_0 - 4.0 * tr_yzz_yyz[i] * tbe_0 + 4.0 * tr_yzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_yzzz[i] = 3.0 * tr_y_yzz[i] - 2.0 * tr_y_yzzzz[i] * tke_0 - 6.0 * tr_yz_yzzz[i] * tbe_0 - 6.0 * tr_yzz_yzz[i] * tbe_0 + 4.0 * tr_yzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yz_zzzz[i] = 4.0 * tr_y_zzz[i] - 2.0 * tr_y_zzzzz[i] * tke_0 - 6.0 * tr_yz_zzzz[i] * tbe_0 - 8.0 * tr_yzz_zzz[i] * tbe_0 + 4.0 * tr_yzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 795-810 components of targeted buffer : DG

    auto tr_z_0_z_zz_xxxx = pbuffer.data(idx_op_geom_110_dg + 795);

    auto tr_z_0_z_zz_xxxy = pbuffer.data(idx_op_geom_110_dg + 796);

    auto tr_z_0_z_zz_xxxz = pbuffer.data(idx_op_geom_110_dg + 797);

    auto tr_z_0_z_zz_xxyy = pbuffer.data(idx_op_geom_110_dg + 798);

    auto tr_z_0_z_zz_xxyz = pbuffer.data(idx_op_geom_110_dg + 799);

    auto tr_z_0_z_zz_xxzz = pbuffer.data(idx_op_geom_110_dg + 800);

    auto tr_z_0_z_zz_xyyy = pbuffer.data(idx_op_geom_110_dg + 801);

    auto tr_z_0_z_zz_xyyz = pbuffer.data(idx_op_geom_110_dg + 802);

    auto tr_z_0_z_zz_xyzz = pbuffer.data(idx_op_geom_110_dg + 803);

    auto tr_z_0_z_zz_xzzz = pbuffer.data(idx_op_geom_110_dg + 804);

    auto tr_z_0_z_zz_yyyy = pbuffer.data(idx_op_geom_110_dg + 805);

    auto tr_z_0_z_zz_yyyz = pbuffer.data(idx_op_geom_110_dg + 806);

    auto tr_z_0_z_zz_yyzz = pbuffer.data(idx_op_geom_110_dg + 807);

    auto tr_z_0_z_zz_yzzz = pbuffer.data(idx_op_geom_110_dg + 808);

    auto tr_z_0_z_zz_zzzz = pbuffer.data(idx_op_geom_110_dg + 809);

    #pragma omp simd aligned(tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xzzz, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yzzz, tr_0_zzzz, tr_z_0_z_zz_xxxx, tr_z_0_z_zz_xxxy, tr_z_0_z_zz_xxxz, tr_z_0_z_zz_xxyy, tr_z_0_z_zz_xxyz, tr_z_0_z_zz_xxzz, tr_z_0_z_zz_xyyy, tr_z_0_z_zz_xyyz, tr_z_0_z_zz_xyzz, tr_z_0_z_zz_xzzz, tr_z_0_z_zz_yyyy, tr_z_0_z_zz_yyyz, tr_z_0_z_zz_yyzz, tr_z_0_z_zz_yzzz, tr_z_0_z_zz_zzzz, tr_z_xxx, tr_z_xxxxz, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, tr_z_zzzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxxxz, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, tr_zzz_zzzzz, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xzzz, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yzzz, tr_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zz_xxxx[i] = 2.0 * tr_0_xxxx[i] - 4.0 * tr_z_xxxxz[i] * tke_0 - 10.0 * tr_zz_xxxx[i] * tbe_0 + 4.0 * tr_zzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xxxy[i] = 2.0 * tr_0_xxxy[i] - 4.0 * tr_z_xxxyz[i] * tke_0 - 10.0 * tr_zz_xxxy[i] * tbe_0 + 4.0 * tr_zzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xxxz[i] = 2.0 * tr_0_xxxz[i] + 2.0 * tr_z_xxx[i] - 4.0 * tr_z_xxxzz[i] * tke_0 - 10.0 * tr_zz_xxxz[i] * tbe_0 - 2.0 * tr_zzz_xxx[i] * tbe_0 + 4.0 * tr_zzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xxyy[i] = 2.0 * tr_0_xxyy[i] - 4.0 * tr_z_xxyyz[i] * tke_0 - 10.0 * tr_zz_xxyy[i] * tbe_0 + 4.0 * tr_zzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xxyz[i] = 2.0 * tr_0_xxyz[i] + 2.0 * tr_z_xxy[i] - 4.0 * tr_z_xxyzz[i] * tke_0 - 10.0 * tr_zz_xxyz[i] * tbe_0 - 2.0 * tr_zzz_xxy[i] * tbe_0 + 4.0 * tr_zzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xxzz[i] = 2.0 * tr_0_xxzz[i] + 4.0 * tr_z_xxz[i] - 4.0 * tr_z_xxzzz[i] * tke_0 - 10.0 * tr_zz_xxzz[i] * tbe_0 - 4.0 * tr_zzz_xxz[i] * tbe_0 + 4.0 * tr_zzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xyyy[i] = 2.0 * tr_0_xyyy[i] - 4.0 * tr_z_xyyyz[i] * tke_0 - 10.0 * tr_zz_xyyy[i] * tbe_0 + 4.0 * tr_zzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xyyz[i] = 2.0 * tr_0_xyyz[i] + 2.0 * tr_z_xyy[i] - 4.0 * tr_z_xyyzz[i] * tke_0 - 10.0 * tr_zz_xyyz[i] * tbe_0 - 2.0 * tr_zzz_xyy[i] * tbe_0 + 4.0 * tr_zzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xyzz[i] = 2.0 * tr_0_xyzz[i] + 4.0 * tr_z_xyz[i] - 4.0 * tr_z_xyzzz[i] * tke_0 - 10.0 * tr_zz_xyzz[i] * tbe_0 - 4.0 * tr_zzz_xyz[i] * tbe_0 + 4.0 * tr_zzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_xzzz[i] = 2.0 * tr_0_xzzz[i] + 6.0 * tr_z_xzz[i] - 4.0 * tr_z_xzzzz[i] * tke_0 - 10.0 * tr_zz_xzzz[i] * tbe_0 - 6.0 * tr_zzz_xzz[i] * tbe_0 + 4.0 * tr_zzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_yyyy[i] = 2.0 * tr_0_yyyy[i] - 4.0 * tr_z_yyyyz[i] * tke_0 - 10.0 * tr_zz_yyyy[i] * tbe_0 + 4.0 * tr_zzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_yyyz[i] = 2.0 * tr_0_yyyz[i] + 2.0 * tr_z_yyy[i] - 4.0 * tr_z_yyyzz[i] * tke_0 - 10.0 * tr_zz_yyyz[i] * tbe_0 - 2.0 * tr_zzz_yyy[i] * tbe_0 + 4.0 * tr_zzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_yyzz[i] = 2.0 * tr_0_yyzz[i] + 4.0 * tr_z_yyz[i] - 4.0 * tr_z_yyzzz[i] * tke_0 - 10.0 * tr_zz_yyzz[i] * tbe_0 - 4.0 * tr_zzz_yyz[i] * tbe_0 + 4.0 * tr_zzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_yzzz[i] = 2.0 * tr_0_yzzz[i] + 6.0 * tr_z_yzz[i] - 4.0 * tr_z_yzzzz[i] * tke_0 - 10.0 * tr_zz_yzzz[i] * tbe_0 - 6.0 * tr_zzz_yzz[i] * tbe_0 + 4.0 * tr_zzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zz_zzzz[i] = 2.0 * tr_0_zzzz[i] + 8.0 * tr_z_zzz[i] - 4.0 * tr_z_zzzzz[i] * tke_0 - 10.0 * tr_zz_zzzz[i] * tbe_0 - 8.0 * tr_zzz_zzz[i] * tbe_0 + 4.0 * tr_zzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_zzzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

