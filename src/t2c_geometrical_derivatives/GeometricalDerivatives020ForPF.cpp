#include "GeometricalDerivatives020ForPF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_pf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_pf,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_pp,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dd,
                         const int idx_op_dg,
                         const int idx_op_ff,
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

    // Set up 0-10 components of targeted buffer : PF

    auto tr_0_0_xx_x_xxx = pbuffer.data(idx_op_geom_020_pf);

    auto tr_0_0_xx_x_xxy = pbuffer.data(idx_op_geom_020_pf + 1);

    auto tr_0_0_xx_x_xxz = pbuffer.data(idx_op_geom_020_pf + 2);

    auto tr_0_0_xx_x_xyy = pbuffer.data(idx_op_geom_020_pf + 3);

    auto tr_0_0_xx_x_xyz = pbuffer.data(idx_op_geom_020_pf + 4);

    auto tr_0_0_xx_x_xzz = pbuffer.data(idx_op_geom_020_pf + 5);

    auto tr_0_0_xx_x_yyy = pbuffer.data(idx_op_geom_020_pf + 6);

    auto tr_0_0_xx_x_yyz = pbuffer.data(idx_op_geom_020_pf + 7);

    auto tr_0_0_xx_x_yzz = pbuffer.data(idx_op_geom_020_pf + 8);

    auto tr_0_0_xx_x_zzz = pbuffer.data(idx_op_geom_020_pf + 9);

    #pragma omp simd aligned(tr_0_0_xx_x_xxx, tr_0_0_xx_x_xxy, tr_0_0_xx_x_xxz, tr_0_0_xx_x_xyy, tr_0_0_xx_x_xyz, tr_0_0_xx_x_xzz, tr_0_0_xx_x_yyy, tr_0_0_xx_x_yyz, tr_0_0_xx_x_yzz, tr_0_0_xx_x_zzz, tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_x, tr_x_xxx, tr_x_xxxxx, tr_x_xxxxy, tr_x_xxxxz, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xx_xx, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_x_xxx[i] = 6.0 * tr_0_xx[i] - 4.0 * tr_0_xxxx[i] * tke_0 + 6.0 * tr_x_x[i] - 6.0 * tr_x_xxx[i] * tbe_0 - 14.0 * tr_x_xxx[i] * tke_0 + 4.0 * tr_x_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xx_xx[i] * tbe_0 + 8.0 * tr_xx_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxy[i] = 4.0 * tr_0_xy[i] - 4.0 * tr_0_xxxy[i] * tke_0 + 2.0 * tr_x_y[i] - 6.0 * tr_x_xxy[i] * tbe_0 - 10.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xx_xy[i] * tbe_0 + 8.0 * tr_xx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxz[i] = 4.0 * tr_0_xz[i] - 4.0 * tr_0_xxxz[i] * tke_0 + 2.0 * tr_x_z[i] - 6.0 * tr_x_xxz[i] * tbe_0 - 10.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xx_xz[i] * tbe_0 + 8.0 * tr_xx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xyy[i] = 2.0 * tr_0_yy[i] - 4.0 * tr_0_xxyy[i] * tke_0 - 6.0 * tr_x_xyy[i] * tbe_0 - 6.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_x_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xx_yy[i] * tbe_0 + 8.0 * tr_xx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xyz[i] = 2.0 * tr_0_yz[i] - 4.0 * tr_0_xxyz[i] * tke_0 - 6.0 * tr_x_xyz[i] * tbe_0 - 6.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xx_yz[i] * tbe_0 + 8.0 * tr_xx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xzz[i] = 2.0 * tr_0_zz[i] - 4.0 * tr_0_xxzz[i] * tke_0 - 6.0 * tr_x_xzz[i] * tbe_0 - 6.0 * tr_x_xzz[i] * tke_0 + 4.0 * tr_x_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xx_zz[i] * tbe_0 + 8.0 * tr_xx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yyy[i] = -4.0 * tr_0_xyyy[i] * tke_0 - 6.0 * tr_x_yyy[i] * tbe_0 - 2.0 * tr_x_yyy[i] * tke_0 + 4.0 * tr_x_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yyz[i] = -4.0 * tr_0_xyyz[i] * tke_0 - 6.0 * tr_x_yyz[i] * tbe_0 - 2.0 * tr_x_yyz[i] * tke_0 + 4.0 * tr_x_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yzz[i] = -4.0 * tr_0_xyzz[i] * tke_0 - 6.0 * tr_x_yzz[i] * tbe_0 - 2.0 * tr_x_yzz[i] * tke_0 + 4.0 * tr_x_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_zzz[i] = -4.0 * tr_0_xzzz[i] * tke_0 - 6.0 * tr_x_zzz[i] * tbe_0 - 2.0 * tr_x_zzz[i] * tke_0 + 4.0 * tr_x_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 10-20 components of targeted buffer : PF

    auto tr_0_0_xx_y_xxx = pbuffer.data(idx_op_geom_020_pf + 10);

    auto tr_0_0_xx_y_xxy = pbuffer.data(idx_op_geom_020_pf + 11);

    auto tr_0_0_xx_y_xxz = pbuffer.data(idx_op_geom_020_pf + 12);

    auto tr_0_0_xx_y_xyy = pbuffer.data(idx_op_geom_020_pf + 13);

    auto tr_0_0_xx_y_xyz = pbuffer.data(idx_op_geom_020_pf + 14);

    auto tr_0_0_xx_y_xzz = pbuffer.data(idx_op_geom_020_pf + 15);

    auto tr_0_0_xx_y_yyy = pbuffer.data(idx_op_geom_020_pf + 16);

    auto tr_0_0_xx_y_yyz = pbuffer.data(idx_op_geom_020_pf + 17);

    auto tr_0_0_xx_y_yzz = pbuffer.data(idx_op_geom_020_pf + 18);

    auto tr_0_0_xx_y_zzz = pbuffer.data(idx_op_geom_020_pf + 19);

    #pragma omp simd aligned(tr_0_0_xx_y_xxx, tr_0_0_xx_y_xxy, tr_0_0_xx_y_xxz, tr_0_0_xx_y_xyy, tr_0_0_xx_y_xyz, tr_0_0_xx_y_xzz, tr_0_0_xx_y_yyy, tr_0_0_xx_y_yyz, tr_0_0_xx_y_yzz, tr_0_0_xx_y_zzz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_x, tr_y_xxx, tr_y_xxxxx, tr_y_xxxxy, tr_y_xxxxz, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_y_xxx[i] = 6.0 * tr_y_x[i] - 2.0 * tr_y_xxx[i] * tbe_0 - 14.0 * tr_y_xxx[i] * tke_0 + 4.0 * tr_y_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xy_xx[i] * tbe_0 + 8.0 * tr_xy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxy[i] = 2.0 * tr_y_y[i] - 2.0 * tr_y_xxy[i] * tbe_0 - 10.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xy_xy[i] * tbe_0 + 8.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxz[i] = 2.0 * tr_y_z[i] - 2.0 * tr_y_xxz[i] * tbe_0 - 10.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xy_xz[i] * tbe_0 + 8.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xyy[i] = -2.0 * tr_y_xyy[i] * tbe_0 - 6.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_y_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xy_yy[i] * tbe_0 + 8.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xyz[i] = -2.0 * tr_y_xyz[i] * tbe_0 - 6.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xy_yz[i] * tbe_0 + 8.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xzz[i] = -2.0 * tr_y_xzz[i] * tbe_0 - 6.0 * tr_y_xzz[i] * tke_0 + 4.0 * tr_y_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_zz[i] * tbe_0 + 8.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yyy[i] = -2.0 * tr_y_yyy[i] * tbe_0 - 2.0 * tr_y_yyy[i] * tke_0 + 4.0 * tr_y_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yyz[i] = -2.0 * tr_y_yyz[i] * tbe_0 - 2.0 * tr_y_yyz[i] * tke_0 + 4.0 * tr_y_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yzz[i] = -2.0 * tr_y_yzz[i] * tbe_0 - 2.0 * tr_y_yzz[i] * tke_0 + 4.0 * tr_y_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_zzz[i] = -2.0 * tr_y_zzz[i] * tbe_0 - 2.0 * tr_y_zzz[i] * tke_0 + 4.0 * tr_y_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 20-30 components of targeted buffer : PF

    auto tr_0_0_xx_z_xxx = pbuffer.data(idx_op_geom_020_pf + 20);

    auto tr_0_0_xx_z_xxy = pbuffer.data(idx_op_geom_020_pf + 21);

    auto tr_0_0_xx_z_xxz = pbuffer.data(idx_op_geom_020_pf + 22);

    auto tr_0_0_xx_z_xyy = pbuffer.data(idx_op_geom_020_pf + 23);

    auto tr_0_0_xx_z_xyz = pbuffer.data(idx_op_geom_020_pf + 24);

    auto tr_0_0_xx_z_xzz = pbuffer.data(idx_op_geom_020_pf + 25);

    auto tr_0_0_xx_z_yyy = pbuffer.data(idx_op_geom_020_pf + 26);

    auto tr_0_0_xx_z_yyz = pbuffer.data(idx_op_geom_020_pf + 27);

    auto tr_0_0_xx_z_yzz = pbuffer.data(idx_op_geom_020_pf + 28);

    auto tr_0_0_xx_z_zzz = pbuffer.data(idx_op_geom_020_pf + 29);

    #pragma omp simd aligned(tr_0_0_xx_z_xxx, tr_0_0_xx_z_xxy, tr_0_0_xx_z_xxz, tr_0_0_xx_z_xyy, tr_0_0_xx_z_xyz, tr_0_0_xx_z_xzz, tr_0_0_xx_z_yyy, tr_0_0_xx_z_yyz, tr_0_0_xx_z_yzz, tr_0_0_xx_z_zzz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_x, tr_z_xxx, tr_z_xxxxx, tr_z_xxxxy, tr_z_xxxxz, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_z_xxx[i] = 6.0 * tr_z_x[i] - 2.0 * tr_z_xxx[i] * tbe_0 - 14.0 * tr_z_xxx[i] * tke_0 + 4.0 * tr_z_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xz_xx[i] * tbe_0 + 8.0 * tr_xz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxy[i] = 2.0 * tr_z_y[i] - 2.0 * tr_z_xxy[i] * tbe_0 - 10.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xz_xy[i] * tbe_0 + 8.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxz[i] = 2.0 * tr_z_z[i] - 2.0 * tr_z_xxz[i] * tbe_0 - 10.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xz[i] * tbe_0 + 8.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xyy[i] = -2.0 * tr_z_xyy[i] * tbe_0 - 6.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_z_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xz_yy[i] * tbe_0 + 8.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xyz[i] = -2.0 * tr_z_xyz[i] * tbe_0 - 6.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yz[i] * tbe_0 + 8.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xzz[i] = -2.0 * tr_z_xzz[i] * tbe_0 - 6.0 * tr_z_xzz[i] * tke_0 + 4.0 * tr_z_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_zz[i] * tbe_0 + 8.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yyy[i] = -2.0 * tr_z_yyy[i] * tbe_0 - 2.0 * tr_z_yyy[i] * tke_0 + 4.0 * tr_z_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yyz[i] = -2.0 * tr_z_yyz[i] * tbe_0 - 2.0 * tr_z_yyz[i] * tke_0 + 4.0 * tr_z_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yzz[i] = -2.0 * tr_z_yzz[i] * tbe_0 - 2.0 * tr_z_yzz[i] * tke_0 + 4.0 * tr_z_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_zzz[i] = -2.0 * tr_z_zzz[i] * tbe_0 - 2.0 * tr_z_zzz[i] * tke_0 + 4.0 * tr_z_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-40 components of targeted buffer : PF

    auto tr_0_0_xy_x_xxx = pbuffer.data(idx_op_geom_020_pf + 30);

    auto tr_0_0_xy_x_xxy = pbuffer.data(idx_op_geom_020_pf + 31);

    auto tr_0_0_xy_x_xxz = pbuffer.data(idx_op_geom_020_pf + 32);

    auto tr_0_0_xy_x_xyy = pbuffer.data(idx_op_geom_020_pf + 33);

    auto tr_0_0_xy_x_xyz = pbuffer.data(idx_op_geom_020_pf + 34);

    auto tr_0_0_xy_x_xzz = pbuffer.data(idx_op_geom_020_pf + 35);

    auto tr_0_0_xy_x_yyy = pbuffer.data(idx_op_geom_020_pf + 36);

    auto tr_0_0_xy_x_yyz = pbuffer.data(idx_op_geom_020_pf + 37);

    auto tr_0_0_xy_x_yzz = pbuffer.data(idx_op_geom_020_pf + 38);

    auto tr_0_0_xy_x_zzz = pbuffer.data(idx_op_geom_020_pf + 39);

    #pragma omp simd aligned(tr_0_0_xy_x_xxx, tr_0_0_xy_x_xxy, tr_0_0_xy_x_xxz, tr_0_0_xy_x_xyy, tr_0_0_xy_x_xyz, tr_0_0_xy_x_xzz, tr_0_0_xy_x_yyy, tr_0_0_xy_x_yyz, tr_0_0_xy_x_yzz, tr_0_0_xy_x_zzz, tr_0_xx, tr_0_xxxy, tr_0_xxyy, tr_0_xxyz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_x_x, tr_x_xxx, tr_x_xxxxy, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_y, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_z, tr_xx_xx, tr_xx_xxxy, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xy, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_yy, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xx, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_x_xxx[i] = -2.0 * tr_0_xxxy[i] * tke_0 - 2.0 * tr_y_xxx[i] * tbe_0 - 6.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxy[i] = tr_0_xx[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_y_xxy[i] * tbe_0 + 2.0 * tr_x_x[i] - 4.0 * tr_x_xyy[i] * tke_0 - 2.0 * tr_x_xxx[i] * tke_0 + 4.0 * tr_x_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxz[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_y_xxz[i] * tbe_0 - 4.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xyy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xyyy[i] * tke_0 - 2.0 * tr_y_xyy[i] * tbe_0 + 2.0 * tr_x_y[i] - 2.0 * tr_x_yyy[i] * tke_0 - 4.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xx_xy[i] * tbe_0 + 4.0 * tr_xx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xyz[i] = tr_0_xz[i] - 2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_y_xyz[i] * tbe_0 + tr_x_z[i] - 2.0 * tr_x_yyz[i] * tke_0 - 2.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xz[i] * tbe_0 + 4.0 * tr_xx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_y_xzz[i] * tbe_0 - 2.0 * tr_x_yzz[i] * tke_0 + 4.0 * tr_x_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yyy[i] = 3.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tke_0 - 2.0 * tr_y_yyy[i] * tbe_0 - 6.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_x_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xx_yy[i] * tbe_0 + 4.0 * tr_xx_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_y_yyz[i] * tbe_0 - 4.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_yz[i] * tbe_0 + 4.0 * tr_xx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yzz[i] = tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_y_yzz[i] * tbe_0 - 2.0 * tr_x_xzz[i] * tke_0 + 4.0 * tr_x_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zz[i] * tbe_0 + 4.0 * tr_xx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_zzz[i] = -2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_x_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 40-50 components of targeted buffer : PF

    auto tr_0_0_xy_y_xxx = pbuffer.data(idx_op_geom_020_pf + 40);

    auto tr_0_0_xy_y_xxy = pbuffer.data(idx_op_geom_020_pf + 41);

    auto tr_0_0_xy_y_xxz = pbuffer.data(idx_op_geom_020_pf + 42);

    auto tr_0_0_xy_y_xyy = pbuffer.data(idx_op_geom_020_pf + 43);

    auto tr_0_0_xy_y_xyz = pbuffer.data(idx_op_geom_020_pf + 44);

    auto tr_0_0_xy_y_xzz = pbuffer.data(idx_op_geom_020_pf + 45);

    auto tr_0_0_xy_y_yyy = pbuffer.data(idx_op_geom_020_pf + 46);

    auto tr_0_0_xy_y_yyz = pbuffer.data(idx_op_geom_020_pf + 47);

    auto tr_0_0_xy_y_yzz = pbuffer.data(idx_op_geom_020_pf + 48);

    auto tr_0_0_xy_y_zzz = pbuffer.data(idx_op_geom_020_pf + 49);

    #pragma omp simd aligned(tr_0_0_xy_y_xxx, tr_0_0_xy_y_xxy, tr_0_0_xy_y_xxz, tr_0_0_xy_y_xyy, tr_0_0_xy_y_xyz, tr_0_0_xy_y_xzz, tr_0_0_xy_y_yyy, tr_0_0_xy_y_yyz, tr_0_0_xy_y_yzz, tr_0_0_xy_y_zzz, tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_y_x, tr_y_xxx, tr_y_xxxxy, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_y, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_z, tr_yy_xx, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_y_xxx[i] = 3.0 * tr_0_xx[i] - 2.0 * tr_0_xxxx[i] * tke_0 - 6.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yy_xx[i] * tbe_0 + 4.0 * tr_yy_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xxxy[i] * tke_0 + 2.0 * tr_y_x[i] - 4.0 * tr_y_xyy[i] * tke_0 - 2.0 * tr_y_xxx[i] * tke_0 + 4.0 * tr_y_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yy_xy[i] * tbe_0 + 4.0 * tr_yy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxy[i] * tbe_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xxxz[i] * tke_0 - 4.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yy_xz[i] * tbe_0 + 4.0 * tr_yy_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xyy[i] = tr_0_yy[i] - 2.0 * tr_0_xxyy[i] * tke_0 + 2.0 * tr_y_y[i] - 2.0 * tr_y_yyy[i] * tke_0 - 4.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yy_yy[i] * tbe_0 + 4.0 * tr_yy_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyy[i] * tbe_0 - 4.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xyz[i] = tr_0_yz[i] - 2.0 * tr_0_xxyz[i] * tke_0 + tr_y_z[i] - 2.0 * tr_y_yyz[i] * tke_0 - 2.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yy_yz[i] * tbe_0 + 4.0 * tr_yy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xzz[i] = tr_0_zz[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_y_yzz[i] * tke_0 + 4.0 * tr_y_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yy_zz[i] * tbe_0 + 4.0 * tr_yy_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yyy[i] = -2.0 * tr_0_xyyy[i] * tke_0 - 6.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_y_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyy[i] * tbe_0 - 6.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yyz[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 4.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyz[i] * tbe_0 - 4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_y_xzz[i] * tke_0 + 4.0 * tr_y_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzz[i] * tbe_0 - 2.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_zzz[i] = -2.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_y_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 50-60 components of targeted buffer : PF

    auto tr_0_0_xy_z_xxx = pbuffer.data(idx_op_geom_020_pf + 50);

    auto tr_0_0_xy_z_xxy = pbuffer.data(idx_op_geom_020_pf + 51);

    auto tr_0_0_xy_z_xxz = pbuffer.data(idx_op_geom_020_pf + 52);

    auto tr_0_0_xy_z_xyy = pbuffer.data(idx_op_geom_020_pf + 53);

    auto tr_0_0_xy_z_xyz = pbuffer.data(idx_op_geom_020_pf + 54);

    auto tr_0_0_xy_z_xzz = pbuffer.data(idx_op_geom_020_pf + 55);

    auto tr_0_0_xy_z_yyy = pbuffer.data(idx_op_geom_020_pf + 56);

    auto tr_0_0_xy_z_yyz = pbuffer.data(idx_op_geom_020_pf + 57);

    auto tr_0_0_xy_z_yzz = pbuffer.data(idx_op_geom_020_pf + 58);

    auto tr_0_0_xy_z_zzz = pbuffer.data(idx_op_geom_020_pf + 59);

    #pragma omp simd aligned(tr_0_0_xy_z_xxx, tr_0_0_xy_z_xxy, tr_0_0_xy_z_xxz, tr_0_0_xy_z_xyy, tr_0_0_xy_z_xyz, tr_0_0_xy_z_xzz, tr_0_0_xy_z_yyy, tr_0_0_xy_z_yyz, tr_0_0_xy_z_yzz, tr_0_0_xy_z_zzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_z_x, tr_z_xxx, tr_z_xxxxy, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_y, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_z_xxx[i] = -6.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxy[i] = 2.0 * tr_z_x[i] - 4.0 * tr_z_xyy[i] * tke_0 - 2.0 * tr_z_xxx[i] * tke_0 + 4.0 * tr_z_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxz[i] = -4.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xyy[i] = 2.0 * tr_z_y[i] - 2.0 * tr_z_yyy[i] * tke_0 - 4.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xyz[i] = tr_z_z[i] - 2.0 * tr_z_yyz[i] * tke_0 - 2.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xzz[i] = -2.0 * tr_z_yzz[i] * tke_0 + 4.0 * tr_z_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yyy[i] = -6.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_z_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yyz[i] = -4.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yzz[i] = -2.0 * tr_z_xzz[i] * tke_0 + 4.0 * tr_z_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_zzz[i] = 4.0 * tr_z_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-70 components of targeted buffer : PF

    auto tr_0_0_xz_x_xxx = pbuffer.data(idx_op_geom_020_pf + 60);

    auto tr_0_0_xz_x_xxy = pbuffer.data(idx_op_geom_020_pf + 61);

    auto tr_0_0_xz_x_xxz = pbuffer.data(idx_op_geom_020_pf + 62);

    auto tr_0_0_xz_x_xyy = pbuffer.data(idx_op_geom_020_pf + 63);

    auto tr_0_0_xz_x_xyz = pbuffer.data(idx_op_geom_020_pf + 64);

    auto tr_0_0_xz_x_xzz = pbuffer.data(idx_op_geom_020_pf + 65);

    auto tr_0_0_xz_x_yyy = pbuffer.data(idx_op_geom_020_pf + 66);

    auto tr_0_0_xz_x_yyz = pbuffer.data(idx_op_geom_020_pf + 67);

    auto tr_0_0_xz_x_yzz = pbuffer.data(idx_op_geom_020_pf + 68);

    auto tr_0_0_xz_x_zzz = pbuffer.data(idx_op_geom_020_pf + 69);

    #pragma omp simd aligned(tr_0_0_xz_x_xxx, tr_0_0_xz_x_xxy, tr_0_0_xz_x_xxz, tr_0_0_xz_x_xyy, tr_0_0_xz_x_xyz, tr_0_0_xz_x_xzz, tr_0_0_xz_x_yyy, tr_0_0_xz_x_yyz, tr_0_0_xz_x_yzz, tr_0_0_xz_x_zzz, tr_0_xx, tr_0_xxxz, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_x_x, tr_x_xxx, tr_x_xxxxz, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_y, tr_x_yyz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xx_xx, tr_xx_xxxz, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xz, tr_xx_xzzz, tr_xx_yy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yz, tr_xx_yzzz, tr_xx_zz, tr_xx_zzzz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xx, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_x_xxx[i] = -2.0 * tr_0_xxxz[i] * tke_0 - 2.0 * tr_z_xxx[i] * tbe_0 - 6.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxy[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_z_xxy[i] * tbe_0 - 4.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxz[i] = tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_z_xxz[i] * tbe_0 + 2.0 * tr_x_x[i] - 4.0 * tr_x_xzz[i] * tke_0 - 2.0 * tr_x_xxx[i] * tke_0 + 4.0 * tr_x_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xx[i] * tbe_0 + 4.0 * tr_xx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xyy[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_z_xyy[i] * tbe_0 - 2.0 * tr_x_yyz[i] * tke_0 + 4.0 * tr_x_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xyz[i] = tr_0_xy[i] - 2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_z_xyz[i] * tbe_0 + tr_x_y[i] - 2.0 * tr_x_yzz[i] * tke_0 - 2.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xy[i] * tbe_0 + 4.0 * tr_xx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xzz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tke_0 - 2.0 * tr_z_xzz[i] * tbe_0 + 2.0 * tr_x_z[i] - 2.0 * tr_x_zzz[i] * tke_0 - 4.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_xz[i] * tbe_0 + 4.0 * tr_xx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yyy[i] = -2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_x_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yyz[i] = tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_z_yyz[i] * tbe_0 - 2.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_x_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yy[i] * tbe_0 + 4.0 * tr_xx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yzz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_z_yzz[i] * tbe_0 - 4.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_yz[i] * tbe_0 + 4.0 * tr_xx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_zzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 2.0 * tr_z_zzz[i] * tbe_0 - 6.0 * tr_x_xzz[i] * tke_0 + 4.0 * tr_x_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xx_zz[i] * tbe_0 + 4.0 * tr_xx_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 70-80 components of targeted buffer : PF

    auto tr_0_0_xz_y_xxx = pbuffer.data(idx_op_geom_020_pf + 70);

    auto tr_0_0_xz_y_xxy = pbuffer.data(idx_op_geom_020_pf + 71);

    auto tr_0_0_xz_y_xxz = pbuffer.data(idx_op_geom_020_pf + 72);

    auto tr_0_0_xz_y_xyy = pbuffer.data(idx_op_geom_020_pf + 73);

    auto tr_0_0_xz_y_xyz = pbuffer.data(idx_op_geom_020_pf + 74);

    auto tr_0_0_xz_y_xzz = pbuffer.data(idx_op_geom_020_pf + 75);

    auto tr_0_0_xz_y_yyy = pbuffer.data(idx_op_geom_020_pf + 76);

    auto tr_0_0_xz_y_yyz = pbuffer.data(idx_op_geom_020_pf + 77);

    auto tr_0_0_xz_y_yzz = pbuffer.data(idx_op_geom_020_pf + 78);

    auto tr_0_0_xz_y_zzz = pbuffer.data(idx_op_geom_020_pf + 79);

    #pragma omp simd aligned(tr_0_0_xz_y_xxx, tr_0_0_xz_y_xxy, tr_0_0_xz_y_xxz, tr_0_0_xz_y_xyy, tr_0_0_xz_y_xyz, tr_0_0_xz_y_xzz, tr_0_0_xz_y_yyy, tr_0_0_xz_y_yyz, tr_0_0_xz_y_yzz, tr_0_0_xz_y_zzz, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_y_x, tr_y_xxx, tr_y_xxxxz, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_y, tr_y_yyz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yz_xx, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_y_xxx[i] = -6.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxy[i] = -4.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxz[i] = 2.0 * tr_y_x[i] - 4.0 * tr_y_xzz[i] * tke_0 - 2.0 * tr_y_xxx[i] * tke_0 + 4.0 * tr_y_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xyy[i] = -2.0 * tr_y_yyz[i] * tke_0 + 4.0 * tr_y_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xyz[i] = tr_y_y[i] - 2.0 * tr_y_yzz[i] * tke_0 - 2.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xzz[i] = 2.0 * tr_y_z[i] - 2.0 * tr_y_zzz[i] * tke_0 - 4.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yyy[i] = 4.0 * tr_y_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yyz[i] = -2.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_y_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yzz[i] = -4.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_zzz[i] = -6.0 * tr_y_xzz[i] * tke_0 + 4.0 * tr_y_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 80-90 components of targeted buffer : PF

    auto tr_0_0_xz_z_xxx = pbuffer.data(idx_op_geom_020_pf + 80);

    auto tr_0_0_xz_z_xxy = pbuffer.data(idx_op_geom_020_pf + 81);

    auto tr_0_0_xz_z_xxz = pbuffer.data(idx_op_geom_020_pf + 82);

    auto tr_0_0_xz_z_xyy = pbuffer.data(idx_op_geom_020_pf + 83);

    auto tr_0_0_xz_z_xyz = pbuffer.data(idx_op_geom_020_pf + 84);

    auto tr_0_0_xz_z_xzz = pbuffer.data(idx_op_geom_020_pf + 85);

    auto tr_0_0_xz_z_yyy = pbuffer.data(idx_op_geom_020_pf + 86);

    auto tr_0_0_xz_z_yyz = pbuffer.data(idx_op_geom_020_pf + 87);

    auto tr_0_0_xz_z_yzz = pbuffer.data(idx_op_geom_020_pf + 88);

    auto tr_0_0_xz_z_zzz = pbuffer.data(idx_op_geom_020_pf + 89);

    #pragma omp simd aligned(tr_0_0_xz_z_xxx, tr_0_0_xz_z_xxy, tr_0_0_xz_z_xxz, tr_0_0_xz_z_xyy, tr_0_0_xz_z_xyz, tr_0_0_xz_z_xzz, tr_0_0_xz_z_yyy, tr_0_0_xz_z_yyz, tr_0_0_xz_z_yzz, tr_0_0_xz_z_zzz, tr_0_xx, tr_0_xxxx, tr_0_xxxy, tr_0_xxxz, tr_0_xxyy, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yz, tr_0_zz, tr_x_xxx, tr_x_xxy, tr_x_xxz, tr_x_xyy, tr_x_xyz, tr_x_xzz, tr_x_yyy, tr_x_yyz, tr_x_yzz, tr_x_zzz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_z_x, tr_z_xxx, tr_z_xxxxz, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_y, tr_z_yyz, tr_z_yzz, tr_z_z, tr_z_zzz, tr_zz_xx, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_z_xxx[i] = 3.0 * tr_0_xx[i] - 2.0 * tr_0_xxxx[i] * tke_0 - 6.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_zz_xx[i] * tbe_0 + 4.0 * tr_zz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xxxy[i] * tke_0 - 4.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xy[i] * tbe_0 + 4.0 * tr_zz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xxxz[i] * tke_0 + 2.0 * tr_z_x[i] - 4.0 * tr_z_xzz[i] * tke_0 - 2.0 * tr_z_xxx[i] * tke_0 + 4.0 * tr_z_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xz[i] * tbe_0 + 4.0 * tr_zz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxz[i] * tbe_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xyy[i] = tr_0_yy[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_z_yyz[i] * tke_0 + 4.0 * tr_z_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zz_yy[i] * tbe_0 + 4.0 * tr_zz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xyz[i] = tr_0_yz[i] - 2.0 * tr_0_xxyz[i] * tke_0 + tr_z_y[i] - 2.0 * tr_z_yzz[i] * tke_0 - 2.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_yz[i] * tbe_0 + 4.0 * tr_zz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyz[i] * tbe_0 - 2.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xzz[i] = tr_0_zz[i] - 2.0 * tr_0_xxzz[i] * tke_0 + 2.0 * tr_z_z[i] - 2.0 * tr_z_zzz[i] * tke_0 - 4.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_zz[i] * tbe_0 + 4.0 * tr_zz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzz[i] * tbe_0 - 4.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yyy[i] = -2.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_z_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yyz[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_z_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyz[i] * tbe_0 - 2.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 4.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzz[i] * tbe_0 - 4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_zzz[i] = -2.0 * tr_0_xzzz[i] * tke_0 - 6.0 * tr_z_xzz[i] * tke_0 + 4.0 * tr_z_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzz[i] * tbe_0 - 6.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-100 components of targeted buffer : PF

    auto tr_0_0_yy_x_xxx = pbuffer.data(idx_op_geom_020_pf + 90);

    auto tr_0_0_yy_x_xxy = pbuffer.data(idx_op_geom_020_pf + 91);

    auto tr_0_0_yy_x_xxz = pbuffer.data(idx_op_geom_020_pf + 92);

    auto tr_0_0_yy_x_xyy = pbuffer.data(idx_op_geom_020_pf + 93);

    auto tr_0_0_yy_x_xyz = pbuffer.data(idx_op_geom_020_pf + 94);

    auto tr_0_0_yy_x_xzz = pbuffer.data(idx_op_geom_020_pf + 95);

    auto tr_0_0_yy_x_yyy = pbuffer.data(idx_op_geom_020_pf + 96);

    auto tr_0_0_yy_x_yyz = pbuffer.data(idx_op_geom_020_pf + 97);

    auto tr_0_0_yy_x_yzz = pbuffer.data(idx_op_geom_020_pf + 98);

    auto tr_0_0_yy_x_zzz = pbuffer.data(idx_op_geom_020_pf + 99);

    #pragma omp simd aligned(tr_0_0_yy_x_xxx, tr_0_0_yy_x_xxy, tr_0_0_yy_x_xxz, tr_0_0_yy_x_xyy, tr_0_0_yy_x_xyz, tr_0_0_yy_x_xzz, tr_0_0_yy_x_yyy, tr_0_0_yy_x_yyz, tr_0_0_yy_x_yzz, tr_0_0_yy_x_zzz, tr_x_x, tr_x_xxx, tr_x_xxxyy, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xzz, tr_x_y, tr_x_yyy, tr_x_yyyyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_z, tr_x_zzz, tr_xy_xx, tr_xy_xxxy, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xy, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_yy, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_x_xxx[i] = -2.0 * tr_x_xxx[i] * tbe_0 - 2.0 * tr_x_xxx[i] * tke_0 + 4.0 * tr_x_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxy[i] = -2.0 * tr_x_xxy[i] * tbe_0 - 6.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xy_xx[i] * tbe_0 + 8.0 * tr_xy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxz[i] = -2.0 * tr_x_xxz[i] * tbe_0 - 2.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xyy[i] = 2.0 * tr_x_x[i] - 2.0 * tr_x_xyy[i] * tbe_0 - 10.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_x_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xy_xy[i] * tbe_0 + 8.0 * tr_xy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xyz[i] = -2.0 * tr_x_xyz[i] * tbe_0 - 6.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 8.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xzz[i] = -2.0 * tr_x_xzz[i] * tbe_0 - 2.0 * tr_x_xzz[i] * tke_0 + 4.0 * tr_x_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yyy[i] = 6.0 * tr_x_y[i] - 2.0 * tr_x_yyy[i] * tbe_0 - 14.0 * tr_x_yyy[i] * tke_0 + 4.0 * tr_x_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xy_yy[i] * tbe_0 + 8.0 * tr_xy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yyz[i] = 2.0 * tr_x_z[i] - 2.0 * tr_x_yyz[i] * tbe_0 - 10.0 * tr_x_yyz[i] * tke_0 + 4.0 * tr_x_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xy_yz[i] * tbe_0 + 8.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yzz[i] = -2.0 * tr_x_yzz[i] * tbe_0 - 6.0 * tr_x_yzz[i] * tke_0 + 4.0 * tr_x_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_zz[i] * tbe_0 + 8.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_zzz[i] = -2.0 * tr_x_zzz[i] * tbe_0 - 2.0 * tr_x_zzz[i] * tke_0 + 4.0 * tr_x_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 100-110 components of targeted buffer : PF

    auto tr_0_0_yy_y_xxx = pbuffer.data(idx_op_geom_020_pf + 100);

    auto tr_0_0_yy_y_xxy = pbuffer.data(idx_op_geom_020_pf + 101);

    auto tr_0_0_yy_y_xxz = pbuffer.data(idx_op_geom_020_pf + 102);

    auto tr_0_0_yy_y_xyy = pbuffer.data(idx_op_geom_020_pf + 103);

    auto tr_0_0_yy_y_xyz = pbuffer.data(idx_op_geom_020_pf + 104);

    auto tr_0_0_yy_y_xzz = pbuffer.data(idx_op_geom_020_pf + 105);

    auto tr_0_0_yy_y_yyy = pbuffer.data(idx_op_geom_020_pf + 106);

    auto tr_0_0_yy_y_yyz = pbuffer.data(idx_op_geom_020_pf + 107);

    auto tr_0_0_yy_y_yzz = pbuffer.data(idx_op_geom_020_pf + 108);

    auto tr_0_0_yy_y_zzz = pbuffer.data(idx_op_geom_020_pf + 109);

    #pragma omp simd aligned(tr_0_0_yy_y_xxx, tr_0_0_yy_y_xxy, tr_0_0_yy_y_xxz, tr_0_0_yy_y_xyy, tr_0_0_yy_y_xyz, tr_0_0_yy_y_xzz, tr_0_0_yy_y_yyy, tr_0_0_yy_y_yyz, tr_0_0_yy_y_yzz, tr_0_0_yy_y_zzz, tr_0_xx, tr_0_xxxy, tr_0_xxyy, tr_0_xxyz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_y_x, tr_y_xxx, tr_y_xxxyy, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xzz, tr_y_y, tr_y_yyy, tr_y_yyyyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_z, tr_y_zzz, tr_yy_xx, tr_yy_xxxy, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xy, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_yy, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_y_xxx[i] = -4.0 * tr_0_xxxy[i] * tke_0 - 6.0 * tr_y_xxx[i] * tbe_0 - 2.0 * tr_y_xxx[i] * tke_0 + 4.0 * tr_y_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxy[i] = 2.0 * tr_0_xx[i] - 4.0 * tr_0_xxyy[i] * tke_0 - 6.0 * tr_y_xxy[i] * tbe_0 - 6.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yy_xx[i] * tbe_0 + 8.0 * tr_yy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxz[i] = -4.0 * tr_0_xxyz[i] * tke_0 - 6.0 * tr_y_xxz[i] * tbe_0 - 2.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xyy[i] = 4.0 * tr_0_xy[i] - 4.0 * tr_0_xyyy[i] * tke_0 + 2.0 * tr_y_x[i] - 6.0 * tr_y_xyy[i] * tbe_0 - 10.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_y_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yy_xy[i] * tbe_0 + 8.0 * tr_yy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xyz[i] = 2.0 * tr_0_xz[i] - 4.0 * tr_0_xyyz[i] * tke_0 - 6.0 * tr_y_xyz[i] * tbe_0 - 6.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yy_xz[i] * tbe_0 + 8.0 * tr_yy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xzz[i] = -4.0 * tr_0_xyzz[i] * tke_0 - 6.0 * tr_y_xzz[i] * tbe_0 - 2.0 * tr_y_xzz[i] * tke_0 + 4.0 * tr_y_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yyy[i] = 6.0 * tr_0_yy[i] - 4.0 * tr_0_yyyy[i] * tke_0 + 6.0 * tr_y_y[i] - 6.0 * tr_y_yyy[i] * tbe_0 - 14.0 * tr_y_yyy[i] * tke_0 + 4.0 * tr_y_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yy_yy[i] * tbe_0 + 8.0 * tr_yy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yyz[i] = 4.0 * tr_0_yz[i] - 4.0 * tr_0_yyyz[i] * tke_0 + 2.0 * tr_y_z[i] - 6.0 * tr_y_yyz[i] * tbe_0 - 10.0 * tr_y_yyz[i] * tke_0 + 4.0 * tr_y_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yy_yz[i] * tbe_0 + 8.0 * tr_yy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yzz[i] = 2.0 * tr_0_zz[i] - 4.0 * tr_0_yyzz[i] * tke_0 - 6.0 * tr_y_yzz[i] * tbe_0 - 6.0 * tr_y_yzz[i] * tke_0 + 4.0 * tr_y_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yy_zz[i] * tbe_0 + 8.0 * tr_yy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_zzz[i] = -4.0 * tr_0_yzzz[i] * tke_0 - 6.0 * tr_y_zzz[i] * tbe_0 - 2.0 * tr_y_zzz[i] * tke_0 + 4.0 * tr_y_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 110-120 components of targeted buffer : PF

    auto tr_0_0_yy_z_xxx = pbuffer.data(idx_op_geom_020_pf + 110);

    auto tr_0_0_yy_z_xxy = pbuffer.data(idx_op_geom_020_pf + 111);

    auto tr_0_0_yy_z_xxz = pbuffer.data(idx_op_geom_020_pf + 112);

    auto tr_0_0_yy_z_xyy = pbuffer.data(idx_op_geom_020_pf + 113);

    auto tr_0_0_yy_z_xyz = pbuffer.data(idx_op_geom_020_pf + 114);

    auto tr_0_0_yy_z_xzz = pbuffer.data(idx_op_geom_020_pf + 115);

    auto tr_0_0_yy_z_yyy = pbuffer.data(idx_op_geom_020_pf + 116);

    auto tr_0_0_yy_z_yyz = pbuffer.data(idx_op_geom_020_pf + 117);

    auto tr_0_0_yy_z_yzz = pbuffer.data(idx_op_geom_020_pf + 118);

    auto tr_0_0_yy_z_zzz = pbuffer.data(idx_op_geom_020_pf + 119);

    #pragma omp simd aligned(tr_0_0_yy_z_xxx, tr_0_0_yy_z_xxy, tr_0_0_yy_z_xxz, tr_0_0_yy_z_xyy, tr_0_0_yy_z_xyz, tr_0_0_yy_z_xzz, tr_0_0_yy_z_yyy, tr_0_0_yy_z_yyz, tr_0_0_yy_z_yzz, tr_0_0_yy_z_zzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_z_x, tr_z_xxx, tr_z_xxxyy, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xzz, tr_z_y, tr_z_yyy, tr_z_yyyyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_z, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_z_xxx[i] = -2.0 * tr_z_xxx[i] * tbe_0 - 2.0 * tr_z_xxx[i] * tke_0 + 4.0 * tr_z_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxy[i] = -2.0 * tr_z_xxy[i] * tbe_0 - 6.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yz_xx[i] * tbe_0 + 8.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxz[i] = -2.0 * tr_z_xxz[i] * tbe_0 - 2.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xyy[i] = 2.0 * tr_z_x[i] - 2.0 * tr_z_xyy[i] * tbe_0 - 10.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_z_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yz_xy[i] * tbe_0 + 8.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xyz[i] = -2.0 * tr_z_xyz[i] * tbe_0 - 6.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xz[i] * tbe_0 + 8.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xzz[i] = -2.0 * tr_z_xzz[i] * tbe_0 - 2.0 * tr_z_xzz[i] * tke_0 + 4.0 * tr_z_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yyy[i] = 6.0 * tr_z_y[i] - 2.0 * tr_z_yyy[i] * tbe_0 - 14.0 * tr_z_yyy[i] * tke_0 + 4.0 * tr_z_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yz_yy[i] * tbe_0 + 8.0 * tr_yz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yyz[i] = 2.0 * tr_z_z[i] - 2.0 * tr_z_yyz[i] * tbe_0 - 10.0 * tr_z_yyz[i] * tke_0 + 4.0 * tr_z_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yz_yz[i] * tbe_0 + 8.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yzz[i] = -2.0 * tr_z_yzz[i] * tbe_0 - 6.0 * tr_z_yzz[i] * tke_0 + 4.0 * tr_z_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_zz[i] * tbe_0 + 8.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_zzz[i] = -2.0 * tr_z_zzz[i] * tbe_0 - 2.0 * tr_z_zzz[i] * tke_0 + 4.0 * tr_z_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-130 components of targeted buffer : PF

    auto tr_0_0_yz_x_xxx = pbuffer.data(idx_op_geom_020_pf + 120);

    auto tr_0_0_yz_x_xxy = pbuffer.data(idx_op_geom_020_pf + 121);

    auto tr_0_0_yz_x_xxz = pbuffer.data(idx_op_geom_020_pf + 122);

    auto tr_0_0_yz_x_xyy = pbuffer.data(idx_op_geom_020_pf + 123);

    auto tr_0_0_yz_x_xyz = pbuffer.data(idx_op_geom_020_pf + 124);

    auto tr_0_0_yz_x_xzz = pbuffer.data(idx_op_geom_020_pf + 125);

    auto tr_0_0_yz_x_yyy = pbuffer.data(idx_op_geom_020_pf + 126);

    auto tr_0_0_yz_x_yyz = pbuffer.data(idx_op_geom_020_pf + 127);

    auto tr_0_0_yz_x_yzz = pbuffer.data(idx_op_geom_020_pf + 128);

    auto tr_0_0_yz_x_zzz = pbuffer.data(idx_op_geom_020_pf + 129);

    #pragma omp simd aligned(tr_0_0_yz_x_xxx, tr_0_0_yz_x_xxy, tr_0_0_yz_x_xxz, tr_0_0_yz_x_xyy, tr_0_0_yz_x_xyz, tr_0_0_yz_x_xzz, tr_0_0_yz_x_yyy, tr_0_0_yz_x_yyz, tr_0_0_yz_x_yzz, tr_0_0_yz_x_zzz, tr_x_x, tr_x_xxxyz, tr_x_xxy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_y, tr_x_yyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_z, tr_x_zzz, tr_xy_xx, tr_xy_xxxz, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xz, tr_xy_xzzz, tr_xy_yy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yz, tr_xy_yzzz, tr_xy_zz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xz_xx, tr_xz_xxxy, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xy, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_yy, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_x_xxx[i] = 4.0 * tr_x_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxy[i] = -2.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xz_xx[i] * tbe_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxz[i] = -2.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xx[i] * tbe_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xyy[i] = -4.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xyz[i] = tr_x_x[i] - 2.0 * tr_x_xzz[i] * tke_0 - 2.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_x_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_xz[i] * tbe_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xy[i] * tbe_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xzz[i] = -4.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xz[i] * tbe_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yyy[i] = -6.0 * tr_x_yyz[i] * tke_0 + 4.0 * tr_x_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xz_yy[i] * tbe_0 + 4.0 * tr_xz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yyz[i] = 2.0 * tr_x_y[i] - 4.0 * tr_x_yzz[i] * tke_0 - 2.0 * tr_x_yyy[i] * tke_0 + 4.0 * tr_x_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yz[i] * tbe_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yy[i] * tbe_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yzz[i] = 2.0 * tr_x_z[i] - 2.0 * tr_x_zzz[i] * tke_0 - 4.0 * tr_x_yyz[i] * tke_0 + 4.0 * tr_x_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_zz[i] * tbe_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yz[i] * tbe_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_zzz[i] = -6.0 * tr_x_yzz[i] * tke_0 + 4.0 * tr_x_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xy_zz[i] * tbe_0 + 4.0 * tr_xy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 130-140 components of targeted buffer : PF

    auto tr_0_0_yz_y_xxx = pbuffer.data(idx_op_geom_020_pf + 130);

    auto tr_0_0_yz_y_xxy = pbuffer.data(idx_op_geom_020_pf + 131);

    auto tr_0_0_yz_y_xxz = pbuffer.data(idx_op_geom_020_pf + 132);

    auto tr_0_0_yz_y_xyy = pbuffer.data(idx_op_geom_020_pf + 133);

    auto tr_0_0_yz_y_xyz = pbuffer.data(idx_op_geom_020_pf + 134);

    auto tr_0_0_yz_y_xzz = pbuffer.data(idx_op_geom_020_pf + 135);

    auto tr_0_0_yz_y_yyy = pbuffer.data(idx_op_geom_020_pf + 136);

    auto tr_0_0_yz_y_yyz = pbuffer.data(idx_op_geom_020_pf + 137);

    auto tr_0_0_yz_y_yzz = pbuffer.data(idx_op_geom_020_pf + 138);

    auto tr_0_0_yz_y_zzz = pbuffer.data(idx_op_geom_020_pf + 139);

    #pragma omp simd aligned(tr_0_0_yz_y_xxx, tr_0_0_yz_y_xxy, tr_0_0_yz_y_xxz, tr_0_0_yz_y_xyy, tr_0_0_yz_y_xyz, tr_0_0_yz_y_xzz, tr_0_0_yz_y_yyy, tr_0_0_yz_y_yyz, tr_0_0_yz_y_yzz, tr_0_0_yz_y_zzz, tr_0_xx, tr_0_xxxz, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_y_x, tr_y_xxxyz, tr_y_xxy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_y, tr_y_yyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_z, tr_y_zzz, tr_yy_xx, tr_yy_xxxz, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xz, tr_yy_xzzz, tr_yy_yy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yz, tr_yy_yzzz, tr_yy_zz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yz_xx, tr_yz_xxxy, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xy, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_yy, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_z_xxx, tr_z_xxy, tr_z_xxz, tr_z_xyy, tr_z_xyz, tr_z_xzz, tr_z_yyy, tr_z_yyz, tr_z_yzz, tr_z_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_y_xxx[i] = -2.0 * tr_0_xxxz[i] * tke_0 - 2.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_y_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxy[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_z_xxy[i] * tbe_0 - 2.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxz[i] = tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_z_xxz[i] * tbe_0 - 2.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xx[i] * tbe_0 + 4.0 * tr_yy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xyy[i] = -2.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_z_xyy[i] * tbe_0 - 4.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xyz[i] = tr_0_xy[i] - 2.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_z_xyz[i] * tbe_0 + tr_y_x[i] - 2.0 * tr_y_xzz[i] * tke_0 - 2.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_y_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xy[i] * tbe_0 + 4.0 * tr_yy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xzz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tke_0 - 2.0 * tr_z_xzz[i] * tbe_0 - 4.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yy_xz[i] * tbe_0 + 4.0 * tr_yy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yyy[i] = -2.0 * tr_0_yyyz[i] * tke_0 - 2.0 * tr_z_yyy[i] * tbe_0 - 6.0 * tr_y_yyz[i] * tke_0 + 4.0 * tr_y_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yyz[i] = tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_z_yyz[i] * tbe_0 + 2.0 * tr_y_y[i] - 4.0 * tr_y_yzz[i] * tke_0 - 2.0 * tr_y_yyy[i] * tke_0 + 4.0 * tr_y_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yy[i] * tbe_0 + 4.0 * tr_yy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yzz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tke_0 - 2.0 * tr_z_yzz[i] * tbe_0 + 2.0 * tr_y_z[i] - 2.0 * tr_y_zzz[i] * tke_0 - 4.0 * tr_y_yyz[i] * tke_0 + 4.0 * tr_y_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yy_yz[i] * tbe_0 + 4.0 * tr_yy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_zzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 2.0 * tr_z_zzz[i] * tbe_0 - 6.0 * tr_y_yzz[i] * tke_0 + 4.0 * tr_y_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yy_zz[i] * tbe_0 + 4.0 * tr_yy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 140-150 components of targeted buffer : PF

    auto tr_0_0_yz_z_xxx = pbuffer.data(idx_op_geom_020_pf + 140);

    auto tr_0_0_yz_z_xxy = pbuffer.data(idx_op_geom_020_pf + 141);

    auto tr_0_0_yz_z_xxz = pbuffer.data(idx_op_geom_020_pf + 142);

    auto tr_0_0_yz_z_xyy = pbuffer.data(idx_op_geom_020_pf + 143);

    auto tr_0_0_yz_z_xyz = pbuffer.data(idx_op_geom_020_pf + 144);

    auto tr_0_0_yz_z_xzz = pbuffer.data(idx_op_geom_020_pf + 145);

    auto tr_0_0_yz_z_yyy = pbuffer.data(idx_op_geom_020_pf + 146);

    auto tr_0_0_yz_z_yyz = pbuffer.data(idx_op_geom_020_pf + 147);

    auto tr_0_0_yz_z_yzz = pbuffer.data(idx_op_geom_020_pf + 148);

    auto tr_0_0_yz_z_zzz = pbuffer.data(idx_op_geom_020_pf + 149);

    #pragma omp simd aligned(tr_0_0_yz_z_xxx, tr_0_0_yz_z_xxy, tr_0_0_yz_z_xxz, tr_0_0_yz_z_xyy, tr_0_0_yz_z_xyz, tr_0_0_yz_z_xzz, tr_0_0_yz_z_yyy, tr_0_0_yz_z_yyz, tr_0_0_yz_z_yzz, tr_0_0_yz_z_zzz, tr_0_xx, tr_0_xxxy, tr_0_xxyy, tr_0_xxyz, tr_0_xy, tr_0_xyyy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_yy, tr_0_yyyy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_y_xxx, tr_y_xxy, tr_y_xxz, tr_y_xyy, tr_y_xyz, tr_y_xzz, tr_y_yyy, tr_y_yyz, tr_y_yzz, tr_y_zzz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_z_x, tr_z_xxxyz, tr_z_xxy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_y, tr_z_yyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_z, tr_z_zzz, tr_zz_xx, tr_zz_xxxy, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xy, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_yy, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_z_xxx[i] = -2.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_z_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxy[i] = tr_0_xx[i] - 2.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zz_xx[i] * tbe_0 + 4.0 * tr_zz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxz[i] = -2.0 * tr_0_xxyz[i] * tke_0 - 2.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxz[i] * tbe_0 - 2.0 * tr_yz_xx[i] * tbe_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xyy[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xyyy[i] * tke_0 - 4.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xy[i] * tbe_0 + 4.0 * tr_zz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xyz[i] = tr_0_xz[i] - 2.0 * tr_0_xyyz[i] * tke_0 + tr_z_x[i] - 2.0 * tr_z_xzz[i] * tke_0 - 2.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_z_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_xz[i] * tbe_0 + 4.0 * tr_zz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyz[i] * tbe_0 - 2.0 * tr_yz_xy[i] * tbe_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xzz[i] = -2.0 * tr_0_xyzz[i] * tke_0 - 4.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xzz[i] * tbe_0 - 4.0 * tr_yz_xz[i] * tbe_0 + 4.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yyy[i] = 3.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tke_0 - 6.0 * tr_z_yyz[i] * tke_0 + 4.0 * tr_z_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_zz_yy[i] * tbe_0 + 4.0 * tr_zz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tke_0 + 2.0 * tr_z_y[i] - 4.0 * tr_z_yzz[i] * tke_0 - 2.0 * tr_z_yyy[i] * tke_0 + 4.0 * tr_z_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_yz[i] * tbe_0 + 4.0 * tr_zz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyz[i] * tbe_0 - 2.0 * tr_yz_yy[i] * tbe_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yzz[i] = tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tke_0 + 2.0 * tr_z_z[i] - 2.0 * tr_z_zzz[i] * tke_0 - 4.0 * tr_z_yyz[i] * tke_0 + 4.0 * tr_z_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_zz[i] * tbe_0 + 4.0 * tr_zz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yzz[i] * tbe_0 - 4.0 * tr_yz_yz[i] * tbe_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_zzz[i] = -2.0 * tr_0_yzzz[i] * tke_0 - 6.0 * tr_z_yzz[i] * tke_0 + 4.0 * tr_z_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_zzz[i] * tbe_0 - 6.0 * tr_yz_zz[i] * tbe_0 + 4.0 * tr_yz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-160 components of targeted buffer : PF

    auto tr_0_0_zz_x_xxx = pbuffer.data(idx_op_geom_020_pf + 150);

    auto tr_0_0_zz_x_xxy = pbuffer.data(idx_op_geom_020_pf + 151);

    auto tr_0_0_zz_x_xxz = pbuffer.data(idx_op_geom_020_pf + 152);

    auto tr_0_0_zz_x_xyy = pbuffer.data(idx_op_geom_020_pf + 153);

    auto tr_0_0_zz_x_xyz = pbuffer.data(idx_op_geom_020_pf + 154);

    auto tr_0_0_zz_x_xzz = pbuffer.data(idx_op_geom_020_pf + 155);

    auto tr_0_0_zz_x_yyy = pbuffer.data(idx_op_geom_020_pf + 156);

    auto tr_0_0_zz_x_yyz = pbuffer.data(idx_op_geom_020_pf + 157);

    auto tr_0_0_zz_x_yzz = pbuffer.data(idx_op_geom_020_pf + 158);

    auto tr_0_0_zz_x_zzz = pbuffer.data(idx_op_geom_020_pf + 159);

    #pragma omp simd aligned(tr_0_0_zz_x_xxx, tr_0_0_zz_x_xxy, tr_0_0_zz_x_xxz, tr_0_0_zz_x_xyy, tr_0_0_zz_x_xyz, tr_0_0_zz_x_xzz, tr_0_0_zz_x_yyy, tr_0_0_zz_x_yyz, tr_0_0_zz_x_yzz, tr_0_0_zz_x_zzz, tr_x_x, tr_x_xxx, tr_x_xxxzz, tr_x_xxy, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_y, tr_x_yyy, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_z, tr_x_zzz, tr_x_zzzzz, tr_xz_xx, tr_xz_xxxz, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xz, tr_xz_xzzz, tr_xz_yy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yz, tr_xz_yzzz, tr_xz_zz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_x_xxx[i] = -2.0 * tr_x_xxx[i] * tbe_0 - 2.0 * tr_x_xxx[i] * tke_0 + 4.0 * tr_x_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxy[i] = -2.0 * tr_x_xxy[i] * tbe_0 - 2.0 * tr_x_xxy[i] * tke_0 + 4.0 * tr_x_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxz[i] = -2.0 * tr_x_xxz[i] * tbe_0 - 6.0 * tr_x_xxz[i] * tke_0 + 4.0 * tr_x_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xx[i] * tbe_0 + 8.0 * tr_xz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xyy[i] = -2.0 * tr_x_xyy[i] * tbe_0 - 2.0 * tr_x_xyy[i] * tke_0 + 4.0 * tr_x_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xyz[i] = -2.0 * tr_x_xyz[i] * tbe_0 - 6.0 * tr_x_xyz[i] * tke_0 + 4.0 * tr_x_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xy[i] * tbe_0 + 8.0 * tr_xz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xzz[i] = 2.0 * tr_x_x[i] - 2.0 * tr_x_xzz[i] * tbe_0 - 10.0 * tr_x_xzz[i] * tke_0 + 4.0 * tr_x_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xz[i] * tbe_0 + 8.0 * tr_xz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yyy[i] = -2.0 * tr_x_yyy[i] * tbe_0 - 2.0 * tr_x_yyy[i] * tke_0 + 4.0 * tr_x_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yyz[i] = -2.0 * tr_x_yyz[i] * tbe_0 - 6.0 * tr_x_yyz[i] * tke_0 + 4.0 * tr_x_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yy[i] * tbe_0 + 8.0 * tr_xz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yzz[i] = 2.0 * tr_x_y[i] - 2.0 * tr_x_yzz[i] * tbe_0 - 10.0 * tr_x_yzz[i] * tke_0 + 4.0 * tr_x_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xz_yz[i] * tbe_0 + 8.0 * tr_xz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_zzz[i] = 6.0 * tr_x_z[i] - 2.0 * tr_x_zzz[i] * tbe_0 - 14.0 * tr_x_zzz[i] * tke_0 + 4.0 * tr_x_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xz_zz[i] * tbe_0 + 8.0 * tr_xz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 160-170 components of targeted buffer : PF

    auto tr_0_0_zz_y_xxx = pbuffer.data(idx_op_geom_020_pf + 160);

    auto tr_0_0_zz_y_xxy = pbuffer.data(idx_op_geom_020_pf + 161);

    auto tr_0_0_zz_y_xxz = pbuffer.data(idx_op_geom_020_pf + 162);

    auto tr_0_0_zz_y_xyy = pbuffer.data(idx_op_geom_020_pf + 163);

    auto tr_0_0_zz_y_xyz = pbuffer.data(idx_op_geom_020_pf + 164);

    auto tr_0_0_zz_y_xzz = pbuffer.data(idx_op_geom_020_pf + 165);

    auto tr_0_0_zz_y_yyy = pbuffer.data(idx_op_geom_020_pf + 166);

    auto tr_0_0_zz_y_yyz = pbuffer.data(idx_op_geom_020_pf + 167);

    auto tr_0_0_zz_y_yzz = pbuffer.data(idx_op_geom_020_pf + 168);

    auto tr_0_0_zz_y_zzz = pbuffer.data(idx_op_geom_020_pf + 169);

    #pragma omp simd aligned(tr_0_0_zz_y_xxx, tr_0_0_zz_y_xxy, tr_0_0_zz_y_xxz, tr_0_0_zz_y_xyy, tr_0_0_zz_y_xyz, tr_0_0_zz_y_xzz, tr_0_0_zz_y_yyy, tr_0_0_zz_y_yyz, tr_0_0_zz_y_yzz, tr_0_0_zz_y_zzz, tr_y_x, tr_y_xxx, tr_y_xxxzz, tr_y_xxy, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_y, tr_y_yyy, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_z, tr_y_zzz, tr_y_zzzzz, tr_yz_xx, tr_yz_xxxz, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xz, tr_yz_xzzz, tr_yz_yy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yz, tr_yz_yzzz, tr_yz_zz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_y_xxx[i] = -2.0 * tr_y_xxx[i] * tbe_0 - 2.0 * tr_y_xxx[i] * tke_0 + 4.0 * tr_y_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxy[i] = -2.0 * tr_y_xxy[i] * tbe_0 - 2.0 * tr_y_xxy[i] * tke_0 + 4.0 * tr_y_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxz[i] = -2.0 * tr_y_xxz[i] * tbe_0 - 6.0 * tr_y_xxz[i] * tke_0 + 4.0 * tr_y_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xx[i] * tbe_0 + 8.0 * tr_yz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xyy[i] = -2.0 * tr_y_xyy[i] * tbe_0 - 2.0 * tr_y_xyy[i] * tke_0 + 4.0 * tr_y_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xyz[i] = -2.0 * tr_y_xyz[i] * tbe_0 - 6.0 * tr_y_xyz[i] * tke_0 + 4.0 * tr_y_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xy[i] * tbe_0 + 8.0 * tr_yz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xzz[i] = 2.0 * tr_y_x[i] - 2.0 * tr_y_xzz[i] * tbe_0 - 10.0 * tr_y_xzz[i] * tke_0 + 4.0 * tr_y_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yz_xz[i] * tbe_0 + 8.0 * tr_yz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yyy[i] = -2.0 * tr_y_yyy[i] * tbe_0 - 2.0 * tr_y_yyy[i] * tke_0 + 4.0 * tr_y_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yyz[i] = -2.0 * tr_y_yyz[i] * tbe_0 - 6.0 * tr_y_yyz[i] * tke_0 + 4.0 * tr_y_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_yy[i] * tbe_0 + 8.0 * tr_yz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yzz[i] = 2.0 * tr_y_y[i] - 2.0 * tr_y_yzz[i] * tbe_0 - 10.0 * tr_y_yzz[i] * tke_0 + 4.0 * tr_y_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yz_yz[i] * tbe_0 + 8.0 * tr_yz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_zzz[i] = 6.0 * tr_y_z[i] - 2.0 * tr_y_zzz[i] * tbe_0 - 14.0 * tr_y_zzz[i] * tke_0 + 4.0 * tr_y_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yz_zz[i] * tbe_0 + 8.0 * tr_yz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 170-180 components of targeted buffer : PF

    auto tr_0_0_zz_z_xxx = pbuffer.data(idx_op_geom_020_pf + 170);

    auto tr_0_0_zz_z_xxy = pbuffer.data(idx_op_geom_020_pf + 171);

    auto tr_0_0_zz_z_xxz = pbuffer.data(idx_op_geom_020_pf + 172);

    auto tr_0_0_zz_z_xyy = pbuffer.data(idx_op_geom_020_pf + 173);

    auto tr_0_0_zz_z_xyz = pbuffer.data(idx_op_geom_020_pf + 174);

    auto tr_0_0_zz_z_xzz = pbuffer.data(idx_op_geom_020_pf + 175);

    auto tr_0_0_zz_z_yyy = pbuffer.data(idx_op_geom_020_pf + 176);

    auto tr_0_0_zz_z_yyz = pbuffer.data(idx_op_geom_020_pf + 177);

    auto tr_0_0_zz_z_yzz = pbuffer.data(idx_op_geom_020_pf + 178);

    auto tr_0_0_zz_z_zzz = pbuffer.data(idx_op_geom_020_pf + 179);

    #pragma omp simd aligned(tr_0_0_zz_z_xxx, tr_0_0_zz_z_xxy, tr_0_0_zz_z_xxz, tr_0_0_zz_z_xyy, tr_0_0_zz_z_xyz, tr_0_0_zz_z_xzz, tr_0_0_zz_z_yyy, tr_0_0_zz_z_yyz, tr_0_0_zz_z_yzz, tr_0_0_zz_z_zzz, tr_0_xx, tr_0_xxxz, tr_0_xxyz, tr_0_xxzz, tr_0_xy, tr_0_xyyz, tr_0_xyzz, tr_0_xz, tr_0_xzzz, tr_0_yy, tr_0_yyyz, tr_0_yyzz, tr_0_yz, tr_0_yzzz, tr_0_zz, tr_0_zzzz, tr_z_x, tr_z_xxx, tr_z_xxxzz, tr_z_xxy, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_y, tr_z_yyy, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_z, tr_z_zzz, tr_z_zzzzz, tr_zz_xx, tr_zz_xxxz, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xz, tr_zz_xzzz, tr_zz_yy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yz, tr_zz_yzzz, tr_zz_zz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_z_xxx[i] = -4.0 * tr_0_xxxz[i] * tke_0 - 6.0 * tr_z_xxx[i] * tbe_0 - 2.0 * tr_z_xxx[i] * tke_0 + 4.0 * tr_z_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxy[i] = -4.0 * tr_0_xxyz[i] * tke_0 - 6.0 * tr_z_xxy[i] * tbe_0 - 2.0 * tr_z_xxy[i] * tke_0 + 4.0 * tr_z_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxz[i] = 2.0 * tr_0_xx[i] - 4.0 * tr_0_xxzz[i] * tke_0 - 6.0 * tr_z_xxz[i] * tbe_0 - 6.0 * tr_z_xxz[i] * tke_0 + 4.0 * tr_z_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xx[i] * tbe_0 + 8.0 * tr_zz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xyy[i] = -4.0 * tr_0_xyyz[i] * tke_0 - 6.0 * tr_z_xyy[i] * tbe_0 - 2.0 * tr_z_xyy[i] * tke_0 + 4.0 * tr_z_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xyz[i] = 2.0 * tr_0_xy[i] - 4.0 * tr_0_xyzz[i] * tke_0 - 6.0 * tr_z_xyz[i] * tbe_0 - 6.0 * tr_z_xyz[i] * tke_0 + 4.0 * tr_z_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xy[i] * tbe_0 + 8.0 * tr_zz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xzz[i] = 4.0 * tr_0_xz[i] - 4.0 * tr_0_xzzz[i] * tke_0 + 2.0 * tr_z_x[i] - 6.0 * tr_z_xzz[i] * tbe_0 - 10.0 * tr_z_xzz[i] * tke_0 + 4.0 * tr_z_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zz_xz[i] * tbe_0 + 8.0 * tr_zz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yyy[i] = -4.0 * tr_0_yyyz[i] * tke_0 - 6.0 * tr_z_yyy[i] * tbe_0 - 2.0 * tr_z_yyy[i] * tke_0 + 4.0 * tr_z_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yyz[i] = 2.0 * tr_0_yy[i] - 4.0 * tr_0_yyzz[i] * tke_0 - 6.0 * tr_z_yyz[i] * tbe_0 - 6.0 * tr_z_yyz[i] * tke_0 + 4.0 * tr_z_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_yy[i] * tbe_0 + 8.0 * tr_zz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yzz[i] = 4.0 * tr_0_yz[i] - 4.0 * tr_0_yzzz[i] * tke_0 + 2.0 * tr_z_y[i] - 6.0 * tr_z_yzz[i] * tbe_0 - 10.0 * tr_z_yzz[i] * tke_0 + 4.0 * tr_z_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zz_yz[i] * tbe_0 + 8.0 * tr_zz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_zzz[i] = 6.0 * tr_0_zz[i] - 4.0 * tr_0_zzzz[i] * tke_0 + 6.0 * tr_z_z[i] - 6.0 * tr_z_zzz[i] * tbe_0 - 14.0 * tr_z_zzz[i] * tke_0 + 4.0 * tr_z_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zz_zz[i] * tbe_0 + 8.0 * tr_zz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_zzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

