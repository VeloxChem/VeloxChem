#include "GeometricalDerivatives010ForDG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_dg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_dg,
                         const int idx_op_pg,
                         const int idx_op_df,
                         const int idx_op_dh,
                         const int idx_op_fg,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : PG

    auto tr_x_xxxx = pbuffer.data(idx_op_pg);

    auto tr_x_xxxy = pbuffer.data(idx_op_pg + 1);

    auto tr_x_xxxz = pbuffer.data(idx_op_pg + 2);

    auto tr_x_xxyy = pbuffer.data(idx_op_pg + 3);

    auto tr_x_xxyz = pbuffer.data(idx_op_pg + 4);

    auto tr_x_xxzz = pbuffer.data(idx_op_pg + 5);

    auto tr_x_xyyy = pbuffer.data(idx_op_pg + 6);

    auto tr_x_xyyz = pbuffer.data(idx_op_pg + 7);

    auto tr_x_xyzz = pbuffer.data(idx_op_pg + 8);

    auto tr_x_xzzz = pbuffer.data(idx_op_pg + 9);

    auto tr_x_yyyy = pbuffer.data(idx_op_pg + 10);

    auto tr_x_yyyz = pbuffer.data(idx_op_pg + 11);

    auto tr_x_yyzz = pbuffer.data(idx_op_pg + 12);

    auto tr_x_yzzz = pbuffer.data(idx_op_pg + 13);

    auto tr_x_zzzz = pbuffer.data(idx_op_pg + 14);

    auto tr_y_xxxx = pbuffer.data(idx_op_pg + 15);

    auto tr_y_xxxy = pbuffer.data(idx_op_pg + 16);

    auto tr_y_xxxz = pbuffer.data(idx_op_pg + 17);

    auto tr_y_xxyy = pbuffer.data(idx_op_pg + 18);

    auto tr_y_xxyz = pbuffer.data(idx_op_pg + 19);

    auto tr_y_xxzz = pbuffer.data(idx_op_pg + 20);

    auto tr_y_xyyy = pbuffer.data(idx_op_pg + 21);

    auto tr_y_xyyz = pbuffer.data(idx_op_pg + 22);

    auto tr_y_xyzz = pbuffer.data(idx_op_pg + 23);

    auto tr_y_xzzz = pbuffer.data(idx_op_pg + 24);

    auto tr_y_yyyy = pbuffer.data(idx_op_pg + 25);

    auto tr_y_yyyz = pbuffer.data(idx_op_pg + 26);

    auto tr_y_yyzz = pbuffer.data(idx_op_pg + 27);

    auto tr_y_yzzz = pbuffer.data(idx_op_pg + 28);

    auto tr_y_zzzz = pbuffer.data(idx_op_pg + 29);

    auto tr_z_xxxx = pbuffer.data(idx_op_pg + 30);

    auto tr_z_xxxy = pbuffer.data(idx_op_pg + 31);

    auto tr_z_xxxz = pbuffer.data(idx_op_pg + 32);

    auto tr_z_xxyy = pbuffer.data(idx_op_pg + 33);

    auto tr_z_xxyz = pbuffer.data(idx_op_pg + 34);

    auto tr_z_xxzz = pbuffer.data(idx_op_pg + 35);

    auto tr_z_xyyy = pbuffer.data(idx_op_pg + 36);

    auto tr_z_xyyz = pbuffer.data(idx_op_pg + 37);

    auto tr_z_xyzz = pbuffer.data(idx_op_pg + 38);

    auto tr_z_xzzz = pbuffer.data(idx_op_pg + 39);

    auto tr_z_yyyy = pbuffer.data(idx_op_pg + 40);

    auto tr_z_yyyz = pbuffer.data(idx_op_pg + 41);

    auto tr_z_yyzz = pbuffer.data(idx_op_pg + 42);

    auto tr_z_yzzz = pbuffer.data(idx_op_pg + 43);

    auto tr_z_zzzz = pbuffer.data(idx_op_pg + 44);

    // Set up components of auxiliary buffer : DF

    auto tr_xx_xxx = pbuffer.data(idx_op_df);

    auto tr_xx_xxy = pbuffer.data(idx_op_df + 1);

    auto tr_xx_xxz = pbuffer.data(idx_op_df + 2);

    auto tr_xx_xyy = pbuffer.data(idx_op_df + 3);

    auto tr_xx_xyz = pbuffer.data(idx_op_df + 4);

    auto tr_xx_xzz = pbuffer.data(idx_op_df + 5);

    auto tr_xx_yyy = pbuffer.data(idx_op_df + 6);

    auto tr_xx_yyz = pbuffer.data(idx_op_df + 7);

    auto tr_xx_yzz = pbuffer.data(idx_op_df + 8);

    auto tr_xx_zzz = pbuffer.data(idx_op_df + 9);

    auto tr_xy_xxx = pbuffer.data(idx_op_df + 10);

    auto tr_xy_xxy = pbuffer.data(idx_op_df + 11);

    auto tr_xy_xxz = pbuffer.data(idx_op_df + 12);

    auto tr_xy_xyy = pbuffer.data(idx_op_df + 13);

    auto tr_xy_xyz = pbuffer.data(idx_op_df + 14);

    auto tr_xy_xzz = pbuffer.data(idx_op_df + 15);

    auto tr_xy_yyy = pbuffer.data(idx_op_df + 16);

    auto tr_xy_yyz = pbuffer.data(idx_op_df + 17);

    auto tr_xy_yzz = pbuffer.data(idx_op_df + 18);

    auto tr_xy_zzz = pbuffer.data(idx_op_df + 19);

    auto tr_xz_xxx = pbuffer.data(idx_op_df + 20);

    auto tr_xz_xxy = pbuffer.data(idx_op_df + 21);

    auto tr_xz_xxz = pbuffer.data(idx_op_df + 22);

    auto tr_xz_xyy = pbuffer.data(idx_op_df + 23);

    auto tr_xz_xyz = pbuffer.data(idx_op_df + 24);

    auto tr_xz_xzz = pbuffer.data(idx_op_df + 25);

    auto tr_xz_yyy = pbuffer.data(idx_op_df + 26);

    auto tr_xz_yyz = pbuffer.data(idx_op_df + 27);

    auto tr_xz_yzz = pbuffer.data(idx_op_df + 28);

    auto tr_xz_zzz = pbuffer.data(idx_op_df + 29);

    auto tr_yy_xxx = pbuffer.data(idx_op_df + 30);

    auto tr_yy_xxy = pbuffer.data(idx_op_df + 31);

    auto tr_yy_xxz = pbuffer.data(idx_op_df + 32);

    auto tr_yy_xyy = pbuffer.data(idx_op_df + 33);

    auto tr_yy_xyz = pbuffer.data(idx_op_df + 34);

    auto tr_yy_xzz = pbuffer.data(idx_op_df + 35);

    auto tr_yy_yyy = pbuffer.data(idx_op_df + 36);

    auto tr_yy_yyz = pbuffer.data(idx_op_df + 37);

    auto tr_yy_yzz = pbuffer.data(idx_op_df + 38);

    auto tr_yy_zzz = pbuffer.data(idx_op_df + 39);

    auto tr_yz_xxx = pbuffer.data(idx_op_df + 40);

    auto tr_yz_xxy = pbuffer.data(idx_op_df + 41);

    auto tr_yz_xxz = pbuffer.data(idx_op_df + 42);

    auto tr_yz_xyy = pbuffer.data(idx_op_df + 43);

    auto tr_yz_xyz = pbuffer.data(idx_op_df + 44);

    auto tr_yz_xzz = pbuffer.data(idx_op_df + 45);

    auto tr_yz_yyy = pbuffer.data(idx_op_df + 46);

    auto tr_yz_yyz = pbuffer.data(idx_op_df + 47);

    auto tr_yz_yzz = pbuffer.data(idx_op_df + 48);

    auto tr_yz_zzz = pbuffer.data(idx_op_df + 49);

    auto tr_zz_xxx = pbuffer.data(idx_op_df + 50);

    auto tr_zz_xxy = pbuffer.data(idx_op_df + 51);

    auto tr_zz_xxz = pbuffer.data(idx_op_df + 52);

    auto tr_zz_xyy = pbuffer.data(idx_op_df + 53);

    auto tr_zz_xyz = pbuffer.data(idx_op_df + 54);

    auto tr_zz_xzz = pbuffer.data(idx_op_df + 55);

    auto tr_zz_yyy = pbuffer.data(idx_op_df + 56);

    auto tr_zz_yyz = pbuffer.data(idx_op_df + 57);

    auto tr_zz_yzz = pbuffer.data(idx_op_df + 58);

    auto tr_zz_zzz = pbuffer.data(idx_op_df + 59);

    // Set up components of auxiliary buffer : DH

    auto tr_xx_xxxxx = pbuffer.data(idx_op_dh);

    auto tr_xx_xxxxy = pbuffer.data(idx_op_dh + 1);

    auto tr_xx_xxxxz = pbuffer.data(idx_op_dh + 2);

    auto tr_xx_xxxyy = pbuffer.data(idx_op_dh + 3);

    auto tr_xx_xxxyz = pbuffer.data(idx_op_dh + 4);

    auto tr_xx_xxxzz = pbuffer.data(idx_op_dh + 5);

    auto tr_xx_xxyyy = pbuffer.data(idx_op_dh + 6);

    auto tr_xx_xxyyz = pbuffer.data(idx_op_dh + 7);

    auto tr_xx_xxyzz = pbuffer.data(idx_op_dh + 8);

    auto tr_xx_xxzzz = pbuffer.data(idx_op_dh + 9);

    auto tr_xx_xyyyy = pbuffer.data(idx_op_dh + 10);

    auto tr_xx_xyyyz = pbuffer.data(idx_op_dh + 11);

    auto tr_xx_xyyzz = pbuffer.data(idx_op_dh + 12);

    auto tr_xx_xyzzz = pbuffer.data(idx_op_dh + 13);

    auto tr_xx_xzzzz = pbuffer.data(idx_op_dh + 14);

    auto tr_xx_yyyyy = pbuffer.data(idx_op_dh + 15);

    auto tr_xx_yyyyz = pbuffer.data(idx_op_dh + 16);

    auto tr_xx_yyyzz = pbuffer.data(idx_op_dh + 17);

    auto tr_xx_yyzzz = pbuffer.data(idx_op_dh + 18);

    auto tr_xx_yzzzz = pbuffer.data(idx_op_dh + 19);

    auto tr_xx_zzzzz = pbuffer.data(idx_op_dh + 20);

    auto tr_xy_xxxxx = pbuffer.data(idx_op_dh + 21);

    auto tr_xy_xxxxy = pbuffer.data(idx_op_dh + 22);

    auto tr_xy_xxxxz = pbuffer.data(idx_op_dh + 23);

    auto tr_xy_xxxyy = pbuffer.data(idx_op_dh + 24);

    auto tr_xy_xxxyz = pbuffer.data(idx_op_dh + 25);

    auto tr_xy_xxxzz = pbuffer.data(idx_op_dh + 26);

    auto tr_xy_xxyyy = pbuffer.data(idx_op_dh + 27);

    auto tr_xy_xxyyz = pbuffer.data(idx_op_dh + 28);

    auto tr_xy_xxyzz = pbuffer.data(idx_op_dh + 29);

    auto tr_xy_xxzzz = pbuffer.data(idx_op_dh + 30);

    auto tr_xy_xyyyy = pbuffer.data(idx_op_dh + 31);

    auto tr_xy_xyyyz = pbuffer.data(idx_op_dh + 32);

    auto tr_xy_xyyzz = pbuffer.data(idx_op_dh + 33);

    auto tr_xy_xyzzz = pbuffer.data(idx_op_dh + 34);

    auto tr_xy_xzzzz = pbuffer.data(idx_op_dh + 35);

    auto tr_xy_yyyyy = pbuffer.data(idx_op_dh + 36);

    auto tr_xy_yyyyz = pbuffer.data(idx_op_dh + 37);

    auto tr_xy_yyyzz = pbuffer.data(idx_op_dh + 38);

    auto tr_xy_yyzzz = pbuffer.data(idx_op_dh + 39);

    auto tr_xy_yzzzz = pbuffer.data(idx_op_dh + 40);

    auto tr_xy_zzzzz = pbuffer.data(idx_op_dh + 41);

    auto tr_xz_xxxxx = pbuffer.data(idx_op_dh + 42);

    auto tr_xz_xxxxy = pbuffer.data(idx_op_dh + 43);

    auto tr_xz_xxxxz = pbuffer.data(idx_op_dh + 44);

    auto tr_xz_xxxyy = pbuffer.data(idx_op_dh + 45);

    auto tr_xz_xxxyz = pbuffer.data(idx_op_dh + 46);

    auto tr_xz_xxxzz = pbuffer.data(idx_op_dh + 47);

    auto tr_xz_xxyyy = pbuffer.data(idx_op_dh + 48);

    auto tr_xz_xxyyz = pbuffer.data(idx_op_dh + 49);

    auto tr_xz_xxyzz = pbuffer.data(idx_op_dh + 50);

    auto tr_xz_xxzzz = pbuffer.data(idx_op_dh + 51);

    auto tr_xz_xyyyy = pbuffer.data(idx_op_dh + 52);

    auto tr_xz_xyyyz = pbuffer.data(idx_op_dh + 53);

    auto tr_xz_xyyzz = pbuffer.data(idx_op_dh + 54);

    auto tr_xz_xyzzz = pbuffer.data(idx_op_dh + 55);

    auto tr_xz_xzzzz = pbuffer.data(idx_op_dh + 56);

    auto tr_xz_yyyyy = pbuffer.data(idx_op_dh + 57);

    auto tr_xz_yyyyz = pbuffer.data(idx_op_dh + 58);

    auto tr_xz_yyyzz = pbuffer.data(idx_op_dh + 59);

    auto tr_xz_yyzzz = pbuffer.data(idx_op_dh + 60);

    auto tr_xz_yzzzz = pbuffer.data(idx_op_dh + 61);

    auto tr_xz_zzzzz = pbuffer.data(idx_op_dh + 62);

    auto tr_yy_xxxxx = pbuffer.data(idx_op_dh + 63);

    auto tr_yy_xxxxy = pbuffer.data(idx_op_dh + 64);

    auto tr_yy_xxxxz = pbuffer.data(idx_op_dh + 65);

    auto tr_yy_xxxyy = pbuffer.data(idx_op_dh + 66);

    auto tr_yy_xxxyz = pbuffer.data(idx_op_dh + 67);

    auto tr_yy_xxxzz = pbuffer.data(idx_op_dh + 68);

    auto tr_yy_xxyyy = pbuffer.data(idx_op_dh + 69);

    auto tr_yy_xxyyz = pbuffer.data(idx_op_dh + 70);

    auto tr_yy_xxyzz = pbuffer.data(idx_op_dh + 71);

    auto tr_yy_xxzzz = pbuffer.data(idx_op_dh + 72);

    auto tr_yy_xyyyy = pbuffer.data(idx_op_dh + 73);

    auto tr_yy_xyyyz = pbuffer.data(idx_op_dh + 74);

    auto tr_yy_xyyzz = pbuffer.data(idx_op_dh + 75);

    auto tr_yy_xyzzz = pbuffer.data(idx_op_dh + 76);

    auto tr_yy_xzzzz = pbuffer.data(idx_op_dh + 77);

    auto tr_yy_yyyyy = pbuffer.data(idx_op_dh + 78);

    auto tr_yy_yyyyz = pbuffer.data(idx_op_dh + 79);

    auto tr_yy_yyyzz = pbuffer.data(idx_op_dh + 80);

    auto tr_yy_yyzzz = pbuffer.data(idx_op_dh + 81);

    auto tr_yy_yzzzz = pbuffer.data(idx_op_dh + 82);

    auto tr_yy_zzzzz = pbuffer.data(idx_op_dh + 83);

    auto tr_yz_xxxxx = pbuffer.data(idx_op_dh + 84);

    auto tr_yz_xxxxy = pbuffer.data(idx_op_dh + 85);

    auto tr_yz_xxxxz = pbuffer.data(idx_op_dh + 86);

    auto tr_yz_xxxyy = pbuffer.data(idx_op_dh + 87);

    auto tr_yz_xxxyz = pbuffer.data(idx_op_dh + 88);

    auto tr_yz_xxxzz = pbuffer.data(idx_op_dh + 89);

    auto tr_yz_xxyyy = pbuffer.data(idx_op_dh + 90);

    auto tr_yz_xxyyz = pbuffer.data(idx_op_dh + 91);

    auto tr_yz_xxyzz = pbuffer.data(idx_op_dh + 92);

    auto tr_yz_xxzzz = pbuffer.data(idx_op_dh + 93);

    auto tr_yz_xyyyy = pbuffer.data(idx_op_dh + 94);

    auto tr_yz_xyyyz = pbuffer.data(idx_op_dh + 95);

    auto tr_yz_xyyzz = pbuffer.data(idx_op_dh + 96);

    auto tr_yz_xyzzz = pbuffer.data(idx_op_dh + 97);

    auto tr_yz_xzzzz = pbuffer.data(idx_op_dh + 98);

    auto tr_yz_yyyyy = pbuffer.data(idx_op_dh + 99);

    auto tr_yz_yyyyz = pbuffer.data(idx_op_dh + 100);

    auto tr_yz_yyyzz = pbuffer.data(idx_op_dh + 101);

    auto tr_yz_yyzzz = pbuffer.data(idx_op_dh + 102);

    auto tr_yz_yzzzz = pbuffer.data(idx_op_dh + 103);

    auto tr_yz_zzzzz = pbuffer.data(idx_op_dh + 104);

    auto tr_zz_xxxxx = pbuffer.data(idx_op_dh + 105);

    auto tr_zz_xxxxy = pbuffer.data(idx_op_dh + 106);

    auto tr_zz_xxxxz = pbuffer.data(idx_op_dh + 107);

    auto tr_zz_xxxyy = pbuffer.data(idx_op_dh + 108);

    auto tr_zz_xxxyz = pbuffer.data(idx_op_dh + 109);

    auto tr_zz_xxxzz = pbuffer.data(idx_op_dh + 110);

    auto tr_zz_xxyyy = pbuffer.data(idx_op_dh + 111);

    auto tr_zz_xxyyz = pbuffer.data(idx_op_dh + 112);

    auto tr_zz_xxyzz = pbuffer.data(idx_op_dh + 113);

    auto tr_zz_xxzzz = pbuffer.data(idx_op_dh + 114);

    auto tr_zz_xyyyy = pbuffer.data(idx_op_dh + 115);

    auto tr_zz_xyyyz = pbuffer.data(idx_op_dh + 116);

    auto tr_zz_xyyzz = pbuffer.data(idx_op_dh + 117);

    auto tr_zz_xyzzz = pbuffer.data(idx_op_dh + 118);

    auto tr_zz_xzzzz = pbuffer.data(idx_op_dh + 119);

    auto tr_zz_yyyyy = pbuffer.data(idx_op_dh + 120);

    auto tr_zz_yyyyz = pbuffer.data(idx_op_dh + 121);

    auto tr_zz_yyyzz = pbuffer.data(idx_op_dh + 122);

    auto tr_zz_yyzzz = pbuffer.data(idx_op_dh + 123);

    auto tr_zz_yzzzz = pbuffer.data(idx_op_dh + 124);

    auto tr_zz_zzzzz = pbuffer.data(idx_op_dh + 125);

    // Set up components of auxiliary buffer : FG

    auto tr_xxx_xxxx = pbuffer.data(idx_op_fg);

    auto tr_xxx_xxxy = pbuffer.data(idx_op_fg + 1);

    auto tr_xxx_xxxz = pbuffer.data(idx_op_fg + 2);

    auto tr_xxx_xxyy = pbuffer.data(idx_op_fg + 3);

    auto tr_xxx_xxyz = pbuffer.data(idx_op_fg + 4);

    auto tr_xxx_xxzz = pbuffer.data(idx_op_fg + 5);

    auto tr_xxx_xyyy = pbuffer.data(idx_op_fg + 6);

    auto tr_xxx_xyyz = pbuffer.data(idx_op_fg + 7);

    auto tr_xxx_xyzz = pbuffer.data(idx_op_fg + 8);

    auto tr_xxx_xzzz = pbuffer.data(idx_op_fg + 9);

    auto tr_xxx_yyyy = pbuffer.data(idx_op_fg + 10);

    auto tr_xxx_yyyz = pbuffer.data(idx_op_fg + 11);

    auto tr_xxx_yyzz = pbuffer.data(idx_op_fg + 12);

    auto tr_xxx_yzzz = pbuffer.data(idx_op_fg + 13);

    auto tr_xxx_zzzz = pbuffer.data(idx_op_fg + 14);

    auto tr_xxy_xxxx = pbuffer.data(idx_op_fg + 15);

    auto tr_xxy_xxxy = pbuffer.data(idx_op_fg + 16);

    auto tr_xxy_xxxz = pbuffer.data(idx_op_fg + 17);

    auto tr_xxy_xxyy = pbuffer.data(idx_op_fg + 18);

    auto tr_xxy_xxyz = pbuffer.data(idx_op_fg + 19);

    auto tr_xxy_xxzz = pbuffer.data(idx_op_fg + 20);

    auto tr_xxy_xyyy = pbuffer.data(idx_op_fg + 21);

    auto tr_xxy_xyyz = pbuffer.data(idx_op_fg + 22);

    auto tr_xxy_xyzz = pbuffer.data(idx_op_fg + 23);

    auto tr_xxy_xzzz = pbuffer.data(idx_op_fg + 24);

    auto tr_xxy_yyyy = pbuffer.data(idx_op_fg + 25);

    auto tr_xxy_yyyz = pbuffer.data(idx_op_fg + 26);

    auto tr_xxy_yyzz = pbuffer.data(idx_op_fg + 27);

    auto tr_xxy_yzzz = pbuffer.data(idx_op_fg + 28);

    auto tr_xxy_zzzz = pbuffer.data(idx_op_fg + 29);

    auto tr_xxz_xxxx = pbuffer.data(idx_op_fg + 30);

    auto tr_xxz_xxxy = pbuffer.data(idx_op_fg + 31);

    auto tr_xxz_xxxz = pbuffer.data(idx_op_fg + 32);

    auto tr_xxz_xxyy = pbuffer.data(idx_op_fg + 33);

    auto tr_xxz_xxyz = pbuffer.data(idx_op_fg + 34);

    auto tr_xxz_xxzz = pbuffer.data(idx_op_fg + 35);

    auto tr_xxz_xyyy = pbuffer.data(idx_op_fg + 36);

    auto tr_xxz_xyyz = pbuffer.data(idx_op_fg + 37);

    auto tr_xxz_xyzz = pbuffer.data(idx_op_fg + 38);

    auto tr_xxz_xzzz = pbuffer.data(idx_op_fg + 39);

    auto tr_xxz_yyyy = pbuffer.data(idx_op_fg + 40);

    auto tr_xxz_yyyz = pbuffer.data(idx_op_fg + 41);

    auto tr_xxz_yyzz = pbuffer.data(idx_op_fg + 42);

    auto tr_xxz_yzzz = pbuffer.data(idx_op_fg + 43);

    auto tr_xxz_zzzz = pbuffer.data(idx_op_fg + 44);

    auto tr_xyy_xxxx = pbuffer.data(idx_op_fg + 45);

    auto tr_xyy_xxxy = pbuffer.data(idx_op_fg + 46);

    auto tr_xyy_xxxz = pbuffer.data(idx_op_fg + 47);

    auto tr_xyy_xxyy = pbuffer.data(idx_op_fg + 48);

    auto tr_xyy_xxyz = pbuffer.data(idx_op_fg + 49);

    auto tr_xyy_xxzz = pbuffer.data(idx_op_fg + 50);

    auto tr_xyy_xyyy = pbuffer.data(idx_op_fg + 51);

    auto tr_xyy_xyyz = pbuffer.data(idx_op_fg + 52);

    auto tr_xyy_xyzz = pbuffer.data(idx_op_fg + 53);

    auto tr_xyy_xzzz = pbuffer.data(idx_op_fg + 54);

    auto tr_xyy_yyyy = pbuffer.data(idx_op_fg + 55);

    auto tr_xyy_yyyz = pbuffer.data(idx_op_fg + 56);

    auto tr_xyy_yyzz = pbuffer.data(idx_op_fg + 57);

    auto tr_xyy_yzzz = pbuffer.data(idx_op_fg + 58);

    auto tr_xyy_zzzz = pbuffer.data(idx_op_fg + 59);

    auto tr_xyz_xxxx = pbuffer.data(idx_op_fg + 60);

    auto tr_xyz_xxxy = pbuffer.data(idx_op_fg + 61);

    auto tr_xyz_xxxz = pbuffer.data(idx_op_fg + 62);

    auto tr_xyz_xxyy = pbuffer.data(idx_op_fg + 63);

    auto tr_xyz_xxyz = pbuffer.data(idx_op_fg + 64);

    auto tr_xyz_xxzz = pbuffer.data(idx_op_fg + 65);

    auto tr_xyz_xyyy = pbuffer.data(idx_op_fg + 66);

    auto tr_xyz_xyyz = pbuffer.data(idx_op_fg + 67);

    auto tr_xyz_xyzz = pbuffer.data(idx_op_fg + 68);

    auto tr_xyz_xzzz = pbuffer.data(idx_op_fg + 69);

    auto tr_xyz_yyyy = pbuffer.data(idx_op_fg + 70);

    auto tr_xyz_yyyz = pbuffer.data(idx_op_fg + 71);

    auto tr_xyz_yyzz = pbuffer.data(idx_op_fg + 72);

    auto tr_xyz_yzzz = pbuffer.data(idx_op_fg + 73);

    auto tr_xyz_zzzz = pbuffer.data(idx_op_fg + 74);

    auto tr_xzz_xxxx = pbuffer.data(idx_op_fg + 75);

    auto tr_xzz_xxxy = pbuffer.data(idx_op_fg + 76);

    auto tr_xzz_xxxz = pbuffer.data(idx_op_fg + 77);

    auto tr_xzz_xxyy = pbuffer.data(idx_op_fg + 78);

    auto tr_xzz_xxyz = pbuffer.data(idx_op_fg + 79);

    auto tr_xzz_xxzz = pbuffer.data(idx_op_fg + 80);

    auto tr_xzz_xyyy = pbuffer.data(idx_op_fg + 81);

    auto tr_xzz_xyyz = pbuffer.data(idx_op_fg + 82);

    auto tr_xzz_xyzz = pbuffer.data(idx_op_fg + 83);

    auto tr_xzz_xzzz = pbuffer.data(idx_op_fg + 84);

    auto tr_xzz_yyyy = pbuffer.data(idx_op_fg + 85);

    auto tr_xzz_yyyz = pbuffer.data(idx_op_fg + 86);

    auto tr_xzz_yyzz = pbuffer.data(idx_op_fg + 87);

    auto tr_xzz_yzzz = pbuffer.data(idx_op_fg + 88);

    auto tr_xzz_zzzz = pbuffer.data(idx_op_fg + 89);

    auto tr_yyy_xxxx = pbuffer.data(idx_op_fg + 90);

    auto tr_yyy_xxxy = pbuffer.data(idx_op_fg + 91);

    auto tr_yyy_xxxz = pbuffer.data(idx_op_fg + 92);

    auto tr_yyy_xxyy = pbuffer.data(idx_op_fg + 93);

    auto tr_yyy_xxyz = pbuffer.data(idx_op_fg + 94);

    auto tr_yyy_xxzz = pbuffer.data(idx_op_fg + 95);

    auto tr_yyy_xyyy = pbuffer.data(idx_op_fg + 96);

    auto tr_yyy_xyyz = pbuffer.data(idx_op_fg + 97);

    auto tr_yyy_xyzz = pbuffer.data(idx_op_fg + 98);

    auto tr_yyy_xzzz = pbuffer.data(idx_op_fg + 99);

    auto tr_yyy_yyyy = pbuffer.data(idx_op_fg + 100);

    auto tr_yyy_yyyz = pbuffer.data(idx_op_fg + 101);

    auto tr_yyy_yyzz = pbuffer.data(idx_op_fg + 102);

    auto tr_yyy_yzzz = pbuffer.data(idx_op_fg + 103);

    auto tr_yyy_zzzz = pbuffer.data(idx_op_fg + 104);

    auto tr_yyz_xxxx = pbuffer.data(idx_op_fg + 105);

    auto tr_yyz_xxxy = pbuffer.data(idx_op_fg + 106);

    auto tr_yyz_xxxz = pbuffer.data(idx_op_fg + 107);

    auto tr_yyz_xxyy = pbuffer.data(idx_op_fg + 108);

    auto tr_yyz_xxyz = pbuffer.data(idx_op_fg + 109);

    auto tr_yyz_xxzz = pbuffer.data(idx_op_fg + 110);

    auto tr_yyz_xyyy = pbuffer.data(idx_op_fg + 111);

    auto tr_yyz_xyyz = pbuffer.data(idx_op_fg + 112);

    auto tr_yyz_xyzz = pbuffer.data(idx_op_fg + 113);

    auto tr_yyz_xzzz = pbuffer.data(idx_op_fg + 114);

    auto tr_yyz_yyyy = pbuffer.data(idx_op_fg + 115);

    auto tr_yyz_yyyz = pbuffer.data(idx_op_fg + 116);

    auto tr_yyz_yyzz = pbuffer.data(idx_op_fg + 117);

    auto tr_yyz_yzzz = pbuffer.data(idx_op_fg + 118);

    auto tr_yyz_zzzz = pbuffer.data(idx_op_fg + 119);

    auto tr_yzz_xxxx = pbuffer.data(idx_op_fg + 120);

    auto tr_yzz_xxxy = pbuffer.data(idx_op_fg + 121);

    auto tr_yzz_xxxz = pbuffer.data(idx_op_fg + 122);

    auto tr_yzz_xxyy = pbuffer.data(idx_op_fg + 123);

    auto tr_yzz_xxyz = pbuffer.data(idx_op_fg + 124);

    auto tr_yzz_xxzz = pbuffer.data(idx_op_fg + 125);

    auto tr_yzz_xyyy = pbuffer.data(idx_op_fg + 126);

    auto tr_yzz_xyyz = pbuffer.data(idx_op_fg + 127);

    auto tr_yzz_xyzz = pbuffer.data(idx_op_fg + 128);

    auto tr_yzz_xzzz = pbuffer.data(idx_op_fg + 129);

    auto tr_yzz_yyyy = pbuffer.data(idx_op_fg + 130);

    auto tr_yzz_yyyz = pbuffer.data(idx_op_fg + 131);

    auto tr_yzz_yyzz = pbuffer.data(idx_op_fg + 132);

    auto tr_yzz_yzzz = pbuffer.data(idx_op_fg + 133);

    auto tr_yzz_zzzz = pbuffer.data(idx_op_fg + 134);

    auto tr_zzz_xxxx = pbuffer.data(idx_op_fg + 135);

    auto tr_zzz_xxxy = pbuffer.data(idx_op_fg + 136);

    auto tr_zzz_xxxz = pbuffer.data(idx_op_fg + 137);

    auto tr_zzz_xxyy = pbuffer.data(idx_op_fg + 138);

    auto tr_zzz_xxyz = pbuffer.data(idx_op_fg + 139);

    auto tr_zzz_xxzz = pbuffer.data(idx_op_fg + 140);

    auto tr_zzz_xyyy = pbuffer.data(idx_op_fg + 141);

    auto tr_zzz_xyyz = pbuffer.data(idx_op_fg + 142);

    auto tr_zzz_xyzz = pbuffer.data(idx_op_fg + 143);

    auto tr_zzz_xzzz = pbuffer.data(idx_op_fg + 144);

    auto tr_zzz_yyyy = pbuffer.data(idx_op_fg + 145);

    auto tr_zzz_yyyz = pbuffer.data(idx_op_fg + 146);

    auto tr_zzz_yyzz = pbuffer.data(idx_op_fg + 147);

    auto tr_zzz_yzzz = pbuffer.data(idx_op_fg + 148);

    auto tr_zzz_zzzz = pbuffer.data(idx_op_fg + 149);

    // Set up 0-15 components of targeted buffer : DG

    auto tr_0_0_x_xx_xxxx = pbuffer.data(idx_op_geom_010_dg);

    auto tr_0_0_x_xx_xxxy = pbuffer.data(idx_op_geom_010_dg + 1);

    auto tr_0_0_x_xx_xxxz = pbuffer.data(idx_op_geom_010_dg + 2);

    auto tr_0_0_x_xx_xxyy = pbuffer.data(idx_op_geom_010_dg + 3);

    auto tr_0_0_x_xx_xxyz = pbuffer.data(idx_op_geom_010_dg + 4);

    auto tr_0_0_x_xx_xxzz = pbuffer.data(idx_op_geom_010_dg + 5);

    auto tr_0_0_x_xx_xyyy = pbuffer.data(idx_op_geom_010_dg + 6);

    auto tr_0_0_x_xx_xyyz = pbuffer.data(idx_op_geom_010_dg + 7);

    auto tr_0_0_x_xx_xyzz = pbuffer.data(idx_op_geom_010_dg + 8);

    auto tr_0_0_x_xx_xzzz = pbuffer.data(idx_op_geom_010_dg + 9);

    auto tr_0_0_x_xx_yyyy = pbuffer.data(idx_op_geom_010_dg + 10);

    auto tr_0_0_x_xx_yyyz = pbuffer.data(idx_op_geom_010_dg + 11);

    auto tr_0_0_x_xx_yyzz = pbuffer.data(idx_op_geom_010_dg + 12);

    auto tr_0_0_x_xx_yzzz = pbuffer.data(idx_op_geom_010_dg + 13);

    auto tr_0_0_x_xx_zzzz = pbuffer.data(idx_op_geom_010_dg + 14);

    #pragma omp simd aligned(tr_0_0_x_xx_xxxx, tr_0_0_x_xx_xxxy, tr_0_0_x_xx_xxxz, tr_0_0_x_xx_xxyy, tr_0_0_x_xx_xxyz, tr_0_0_x_xx_xxzz, tr_0_0_x_xx_xyyy, tr_0_0_x_xx_xyyz, tr_0_0_x_xx_xyzz, tr_0_0_x_xx_xzzz, tr_0_0_x_xx_yyyy, tr_0_0_x_xx_yyyz, tr_0_0_x_xx_yyzz, tr_0_0_x_xx_yzzz, tr_0_0_x_xx_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xx_xxxx[i] = 2.0 * tr_xxx_xxxx[i] * tbe_0 + 2.0 * tr_xx_xxxxx[i] * tke_0 - 2.0 * tr_x_xxxx[i] - 4.0 * tr_xx_xxx[i];

        tr_0_0_x_xx_xxxy[i] = 2.0 * tr_xxx_xxxy[i] * tbe_0 + 2.0 * tr_xx_xxxxy[i] * tke_0 - 2.0 * tr_x_xxxy[i] - 3.0 * tr_xx_xxy[i];

        tr_0_0_x_xx_xxxz[i] = 2.0 * tr_xxx_xxxz[i] * tbe_0 + 2.0 * tr_xx_xxxxz[i] * tke_0 - 2.0 * tr_x_xxxz[i] - 3.0 * tr_xx_xxz[i];

        tr_0_0_x_xx_xxyy[i] = 2.0 * tr_xxx_xxyy[i] * tbe_0 + 2.0 * tr_xx_xxxyy[i] * tke_0 - 2.0 * tr_x_xxyy[i] - 2.0 * tr_xx_xyy[i];

        tr_0_0_x_xx_xxyz[i] = 2.0 * tr_xxx_xxyz[i] * tbe_0 + 2.0 * tr_xx_xxxyz[i] * tke_0 - 2.0 * tr_x_xxyz[i] - 2.0 * tr_xx_xyz[i];

        tr_0_0_x_xx_xxzz[i] = 2.0 * tr_xxx_xxzz[i] * tbe_0 + 2.0 * tr_xx_xxxzz[i] * tke_0 - 2.0 * tr_x_xxzz[i] - 2.0 * tr_xx_xzz[i];

        tr_0_0_x_xx_xyyy[i] = 2.0 * tr_xxx_xyyy[i] * tbe_0 + 2.0 * tr_xx_xxyyy[i] * tke_0 - 2.0 * tr_x_xyyy[i] - tr_xx_yyy[i];

        tr_0_0_x_xx_xyyz[i] = 2.0 * tr_xxx_xyyz[i] * tbe_0 + 2.0 * tr_xx_xxyyz[i] * tke_0 - 2.0 * tr_x_xyyz[i] - tr_xx_yyz[i];

        tr_0_0_x_xx_xyzz[i] = 2.0 * tr_xxx_xyzz[i] * tbe_0 + 2.0 * tr_xx_xxyzz[i] * tke_0 - 2.0 * tr_x_xyzz[i] - tr_xx_yzz[i];

        tr_0_0_x_xx_xzzz[i] = 2.0 * tr_xxx_xzzz[i] * tbe_0 + 2.0 * tr_xx_xxzzz[i] * tke_0 - 2.0 * tr_x_xzzz[i] - tr_xx_zzz[i];

        tr_0_0_x_xx_yyyy[i] = 2.0 * tr_xxx_yyyy[i] * tbe_0 + 2.0 * tr_xx_xyyyy[i] * tke_0 - 2.0 * tr_x_yyyy[i];

        tr_0_0_x_xx_yyyz[i] = 2.0 * tr_xxx_yyyz[i] * tbe_0 + 2.0 * tr_xx_xyyyz[i] * tke_0 - 2.0 * tr_x_yyyz[i];

        tr_0_0_x_xx_yyzz[i] = 2.0 * tr_xxx_yyzz[i] * tbe_0 + 2.0 * tr_xx_xyyzz[i] * tke_0 - 2.0 * tr_x_yyzz[i];

        tr_0_0_x_xx_yzzz[i] = 2.0 * tr_xxx_yzzz[i] * tbe_0 + 2.0 * tr_xx_xyzzz[i] * tke_0 - 2.0 * tr_x_yzzz[i];

        tr_0_0_x_xx_zzzz[i] = 2.0 * tr_xxx_zzzz[i] * tbe_0 + 2.0 * tr_xx_xzzzz[i] * tke_0 - 2.0 * tr_x_zzzz[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto tr_0_0_x_xy_xxxx = pbuffer.data(idx_op_geom_010_dg + 15);

    auto tr_0_0_x_xy_xxxy = pbuffer.data(idx_op_geom_010_dg + 16);

    auto tr_0_0_x_xy_xxxz = pbuffer.data(idx_op_geom_010_dg + 17);

    auto tr_0_0_x_xy_xxyy = pbuffer.data(idx_op_geom_010_dg + 18);

    auto tr_0_0_x_xy_xxyz = pbuffer.data(idx_op_geom_010_dg + 19);

    auto tr_0_0_x_xy_xxzz = pbuffer.data(idx_op_geom_010_dg + 20);

    auto tr_0_0_x_xy_xyyy = pbuffer.data(idx_op_geom_010_dg + 21);

    auto tr_0_0_x_xy_xyyz = pbuffer.data(idx_op_geom_010_dg + 22);

    auto tr_0_0_x_xy_xyzz = pbuffer.data(idx_op_geom_010_dg + 23);

    auto tr_0_0_x_xy_xzzz = pbuffer.data(idx_op_geom_010_dg + 24);

    auto tr_0_0_x_xy_yyyy = pbuffer.data(idx_op_geom_010_dg + 25);

    auto tr_0_0_x_xy_yyyz = pbuffer.data(idx_op_geom_010_dg + 26);

    auto tr_0_0_x_xy_yyzz = pbuffer.data(idx_op_geom_010_dg + 27);

    auto tr_0_0_x_xy_yzzz = pbuffer.data(idx_op_geom_010_dg + 28);

    auto tr_0_0_x_xy_zzzz = pbuffer.data(idx_op_geom_010_dg + 29);

    #pragma omp simd aligned(tr_0_0_x_xy_xxxx, tr_0_0_x_xy_xxxy, tr_0_0_x_xy_xxxz, tr_0_0_x_xy_xxyy, tr_0_0_x_xy_xxyz, tr_0_0_x_xy_xxzz, tr_0_0_x_xy_xyyy, tr_0_0_x_xy_xyyz, tr_0_0_x_xy_xyzz, tr_0_0_x_xy_xzzz, tr_0_0_x_xy_yyyy, tr_0_0_x_xy_yyyz, tr_0_0_x_xy_yyzz, tr_0_0_x_xy_yzzz, tr_0_0_x_xy_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xy_xxxx[i] = 2.0 * tr_xxy_xxxx[i] * tbe_0 + 2.0 * tr_xy_xxxxx[i] * tke_0 - tr_y_xxxx[i] - 4.0 * tr_xy_xxx[i];

        tr_0_0_x_xy_xxxy[i] = 2.0 * tr_xxy_xxxy[i] * tbe_0 + 2.0 * tr_xy_xxxxy[i] * tke_0 - tr_y_xxxy[i] - 3.0 * tr_xy_xxy[i];

        tr_0_0_x_xy_xxxz[i] = 2.0 * tr_xxy_xxxz[i] * tbe_0 + 2.0 * tr_xy_xxxxz[i] * tke_0 - tr_y_xxxz[i] - 3.0 * tr_xy_xxz[i];

        tr_0_0_x_xy_xxyy[i] = 2.0 * tr_xxy_xxyy[i] * tbe_0 + 2.0 * tr_xy_xxxyy[i] * tke_0 - tr_y_xxyy[i] - 2.0 * tr_xy_xyy[i];

        tr_0_0_x_xy_xxyz[i] = 2.0 * tr_xxy_xxyz[i] * tbe_0 + 2.0 * tr_xy_xxxyz[i] * tke_0 - tr_y_xxyz[i] - 2.0 * tr_xy_xyz[i];

        tr_0_0_x_xy_xxzz[i] = 2.0 * tr_xxy_xxzz[i] * tbe_0 + 2.0 * tr_xy_xxxzz[i] * tke_0 - tr_y_xxzz[i] - 2.0 * tr_xy_xzz[i];

        tr_0_0_x_xy_xyyy[i] = 2.0 * tr_xxy_xyyy[i] * tbe_0 + 2.0 * tr_xy_xxyyy[i] * tke_0 - tr_y_xyyy[i] - tr_xy_yyy[i];

        tr_0_0_x_xy_xyyz[i] = 2.0 * tr_xxy_xyyz[i] * tbe_0 + 2.0 * tr_xy_xxyyz[i] * tke_0 - tr_y_xyyz[i] - tr_xy_yyz[i];

        tr_0_0_x_xy_xyzz[i] = 2.0 * tr_xxy_xyzz[i] * tbe_0 + 2.0 * tr_xy_xxyzz[i] * tke_0 - tr_y_xyzz[i] - tr_xy_yzz[i];

        tr_0_0_x_xy_xzzz[i] = 2.0 * tr_xxy_xzzz[i] * tbe_0 + 2.0 * tr_xy_xxzzz[i] * tke_0 - tr_y_xzzz[i] - tr_xy_zzz[i];

        tr_0_0_x_xy_yyyy[i] = 2.0 * tr_xxy_yyyy[i] * tbe_0 + 2.0 * tr_xy_xyyyy[i] * tke_0 - tr_y_yyyy[i];

        tr_0_0_x_xy_yyyz[i] = 2.0 * tr_xxy_yyyz[i] * tbe_0 + 2.0 * tr_xy_xyyyz[i] * tke_0 - tr_y_yyyz[i];

        tr_0_0_x_xy_yyzz[i] = 2.0 * tr_xxy_yyzz[i] * tbe_0 + 2.0 * tr_xy_xyyzz[i] * tke_0 - tr_y_yyzz[i];

        tr_0_0_x_xy_yzzz[i] = 2.0 * tr_xxy_yzzz[i] * tbe_0 + 2.0 * tr_xy_xyzzz[i] * tke_0 - tr_y_yzzz[i];

        tr_0_0_x_xy_zzzz[i] = 2.0 * tr_xxy_zzzz[i] * tbe_0 + 2.0 * tr_xy_xzzzz[i] * tke_0 - tr_y_zzzz[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto tr_0_0_x_xz_xxxx = pbuffer.data(idx_op_geom_010_dg + 30);

    auto tr_0_0_x_xz_xxxy = pbuffer.data(idx_op_geom_010_dg + 31);

    auto tr_0_0_x_xz_xxxz = pbuffer.data(idx_op_geom_010_dg + 32);

    auto tr_0_0_x_xz_xxyy = pbuffer.data(idx_op_geom_010_dg + 33);

    auto tr_0_0_x_xz_xxyz = pbuffer.data(idx_op_geom_010_dg + 34);

    auto tr_0_0_x_xz_xxzz = pbuffer.data(idx_op_geom_010_dg + 35);

    auto tr_0_0_x_xz_xyyy = pbuffer.data(idx_op_geom_010_dg + 36);

    auto tr_0_0_x_xz_xyyz = pbuffer.data(idx_op_geom_010_dg + 37);

    auto tr_0_0_x_xz_xyzz = pbuffer.data(idx_op_geom_010_dg + 38);

    auto tr_0_0_x_xz_xzzz = pbuffer.data(idx_op_geom_010_dg + 39);

    auto tr_0_0_x_xz_yyyy = pbuffer.data(idx_op_geom_010_dg + 40);

    auto tr_0_0_x_xz_yyyz = pbuffer.data(idx_op_geom_010_dg + 41);

    auto tr_0_0_x_xz_yyzz = pbuffer.data(idx_op_geom_010_dg + 42);

    auto tr_0_0_x_xz_yzzz = pbuffer.data(idx_op_geom_010_dg + 43);

    auto tr_0_0_x_xz_zzzz = pbuffer.data(idx_op_geom_010_dg + 44);

    #pragma omp simd aligned(tr_0_0_x_xz_xxxx, tr_0_0_x_xz_xxxy, tr_0_0_x_xz_xxxz, tr_0_0_x_xz_xxyy, tr_0_0_x_xz_xxyz, tr_0_0_x_xz_xxzz, tr_0_0_x_xz_xyyy, tr_0_0_x_xz_xyyz, tr_0_0_x_xz_xyzz, tr_0_0_x_xz_xzzz, tr_0_0_x_xz_yyyy, tr_0_0_x_xz_yyyz, tr_0_0_x_xz_yyzz, tr_0_0_x_xz_yzzz, tr_0_0_x_xz_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xz_xxxx[i] = 2.0 * tr_xxz_xxxx[i] * tbe_0 + 2.0 * tr_xz_xxxxx[i] * tke_0 - tr_z_xxxx[i] - 4.0 * tr_xz_xxx[i];

        tr_0_0_x_xz_xxxy[i] = 2.0 * tr_xxz_xxxy[i] * tbe_0 + 2.0 * tr_xz_xxxxy[i] * tke_0 - tr_z_xxxy[i] - 3.0 * tr_xz_xxy[i];

        tr_0_0_x_xz_xxxz[i] = 2.0 * tr_xxz_xxxz[i] * tbe_0 + 2.0 * tr_xz_xxxxz[i] * tke_0 - tr_z_xxxz[i] - 3.0 * tr_xz_xxz[i];

        tr_0_0_x_xz_xxyy[i] = 2.0 * tr_xxz_xxyy[i] * tbe_0 + 2.0 * tr_xz_xxxyy[i] * tke_0 - tr_z_xxyy[i] - 2.0 * tr_xz_xyy[i];

        tr_0_0_x_xz_xxyz[i] = 2.0 * tr_xxz_xxyz[i] * tbe_0 + 2.0 * tr_xz_xxxyz[i] * tke_0 - tr_z_xxyz[i] - 2.0 * tr_xz_xyz[i];

        tr_0_0_x_xz_xxzz[i] = 2.0 * tr_xxz_xxzz[i] * tbe_0 + 2.0 * tr_xz_xxxzz[i] * tke_0 - tr_z_xxzz[i] - 2.0 * tr_xz_xzz[i];

        tr_0_0_x_xz_xyyy[i] = 2.0 * tr_xxz_xyyy[i] * tbe_0 + 2.0 * tr_xz_xxyyy[i] * tke_0 - tr_z_xyyy[i] - tr_xz_yyy[i];

        tr_0_0_x_xz_xyyz[i] = 2.0 * tr_xxz_xyyz[i] * tbe_0 + 2.0 * tr_xz_xxyyz[i] * tke_0 - tr_z_xyyz[i] - tr_xz_yyz[i];

        tr_0_0_x_xz_xyzz[i] = 2.0 * tr_xxz_xyzz[i] * tbe_0 + 2.0 * tr_xz_xxyzz[i] * tke_0 - tr_z_xyzz[i] - tr_xz_yzz[i];

        tr_0_0_x_xz_xzzz[i] = 2.0 * tr_xxz_xzzz[i] * tbe_0 + 2.0 * tr_xz_xxzzz[i] * tke_0 - tr_z_xzzz[i] - tr_xz_zzz[i];

        tr_0_0_x_xz_yyyy[i] = 2.0 * tr_xxz_yyyy[i] * tbe_0 + 2.0 * tr_xz_xyyyy[i] * tke_0 - tr_z_yyyy[i];

        tr_0_0_x_xz_yyyz[i] = 2.0 * tr_xxz_yyyz[i] * tbe_0 + 2.0 * tr_xz_xyyyz[i] * tke_0 - tr_z_yyyz[i];

        tr_0_0_x_xz_yyzz[i] = 2.0 * tr_xxz_yyzz[i] * tbe_0 + 2.0 * tr_xz_xyyzz[i] * tke_0 - tr_z_yyzz[i];

        tr_0_0_x_xz_yzzz[i] = 2.0 * tr_xxz_yzzz[i] * tbe_0 + 2.0 * tr_xz_xyzzz[i] * tke_0 - tr_z_yzzz[i];

        tr_0_0_x_xz_zzzz[i] = 2.0 * tr_xxz_zzzz[i] * tbe_0 + 2.0 * tr_xz_xzzzz[i] * tke_0 - tr_z_zzzz[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto tr_0_0_x_yy_xxxx = pbuffer.data(idx_op_geom_010_dg + 45);

    auto tr_0_0_x_yy_xxxy = pbuffer.data(idx_op_geom_010_dg + 46);

    auto tr_0_0_x_yy_xxxz = pbuffer.data(idx_op_geom_010_dg + 47);

    auto tr_0_0_x_yy_xxyy = pbuffer.data(idx_op_geom_010_dg + 48);

    auto tr_0_0_x_yy_xxyz = pbuffer.data(idx_op_geom_010_dg + 49);

    auto tr_0_0_x_yy_xxzz = pbuffer.data(idx_op_geom_010_dg + 50);

    auto tr_0_0_x_yy_xyyy = pbuffer.data(idx_op_geom_010_dg + 51);

    auto tr_0_0_x_yy_xyyz = pbuffer.data(idx_op_geom_010_dg + 52);

    auto tr_0_0_x_yy_xyzz = pbuffer.data(idx_op_geom_010_dg + 53);

    auto tr_0_0_x_yy_xzzz = pbuffer.data(idx_op_geom_010_dg + 54);

    auto tr_0_0_x_yy_yyyy = pbuffer.data(idx_op_geom_010_dg + 55);

    auto tr_0_0_x_yy_yyyz = pbuffer.data(idx_op_geom_010_dg + 56);

    auto tr_0_0_x_yy_yyzz = pbuffer.data(idx_op_geom_010_dg + 57);

    auto tr_0_0_x_yy_yzzz = pbuffer.data(idx_op_geom_010_dg + 58);

    auto tr_0_0_x_yy_zzzz = pbuffer.data(idx_op_geom_010_dg + 59);

    #pragma omp simd aligned(tr_0_0_x_yy_xxxx, tr_0_0_x_yy_xxxy, tr_0_0_x_yy_xxxz, tr_0_0_x_yy_xxyy, tr_0_0_x_yy_xxyz, tr_0_0_x_yy_xxzz, tr_0_0_x_yy_xyyy, tr_0_0_x_yy_xyyz, tr_0_0_x_yy_xyzz, tr_0_0_x_yy_xzzz, tr_0_0_x_yy_yyyy, tr_0_0_x_yy_yyyz, tr_0_0_x_yy_yyzz, tr_0_0_x_yy_yzzz, tr_0_0_x_yy_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yy_xxxx[i] = 2.0 * tr_xyy_xxxx[i] * tbe_0 + 2.0 * tr_yy_xxxxx[i] * tke_0 - 4.0 * tr_yy_xxx[i];

        tr_0_0_x_yy_xxxy[i] = 2.0 * tr_xyy_xxxy[i] * tbe_0 + 2.0 * tr_yy_xxxxy[i] * tke_0 - 3.0 * tr_yy_xxy[i];

        tr_0_0_x_yy_xxxz[i] = 2.0 * tr_xyy_xxxz[i] * tbe_0 + 2.0 * tr_yy_xxxxz[i] * tke_0 - 3.0 * tr_yy_xxz[i];

        tr_0_0_x_yy_xxyy[i] = 2.0 * tr_xyy_xxyy[i] * tbe_0 + 2.0 * tr_yy_xxxyy[i] * tke_0 - 2.0 * tr_yy_xyy[i];

        tr_0_0_x_yy_xxyz[i] = 2.0 * tr_xyy_xxyz[i] * tbe_0 + 2.0 * tr_yy_xxxyz[i] * tke_0 - 2.0 * tr_yy_xyz[i];

        tr_0_0_x_yy_xxzz[i] = 2.0 * tr_xyy_xxzz[i] * tbe_0 + 2.0 * tr_yy_xxxzz[i] * tke_0 - 2.0 * tr_yy_xzz[i];

        tr_0_0_x_yy_xyyy[i] = 2.0 * tr_xyy_xyyy[i] * tbe_0 + 2.0 * tr_yy_xxyyy[i] * tke_0 - tr_yy_yyy[i];

        tr_0_0_x_yy_xyyz[i] = 2.0 * tr_xyy_xyyz[i] * tbe_0 + 2.0 * tr_yy_xxyyz[i] * tke_0 - tr_yy_yyz[i];

        tr_0_0_x_yy_xyzz[i] = 2.0 * tr_xyy_xyzz[i] * tbe_0 + 2.0 * tr_yy_xxyzz[i] * tke_0 - tr_yy_yzz[i];

        tr_0_0_x_yy_xzzz[i] = 2.0 * tr_xyy_xzzz[i] * tbe_0 + 2.0 * tr_yy_xxzzz[i] * tke_0 - tr_yy_zzz[i];

        tr_0_0_x_yy_yyyy[i] = 2.0 * tr_xyy_yyyy[i] * tbe_0 + 2.0 * tr_yy_xyyyy[i] * tke_0;

        tr_0_0_x_yy_yyyz[i] = 2.0 * tr_xyy_yyyz[i] * tbe_0 + 2.0 * tr_yy_xyyyz[i] * tke_0;

        tr_0_0_x_yy_yyzz[i] = 2.0 * tr_xyy_yyzz[i] * tbe_0 + 2.0 * tr_yy_xyyzz[i] * tke_0;

        tr_0_0_x_yy_yzzz[i] = 2.0 * tr_xyy_yzzz[i] * tbe_0 + 2.0 * tr_yy_xyzzz[i] * tke_0;

        tr_0_0_x_yy_zzzz[i] = 2.0 * tr_xyy_zzzz[i] * tbe_0 + 2.0 * tr_yy_xzzzz[i] * tke_0;
    }

    // Set up 60-75 components of targeted buffer : DG

    auto tr_0_0_x_yz_xxxx = pbuffer.data(idx_op_geom_010_dg + 60);

    auto tr_0_0_x_yz_xxxy = pbuffer.data(idx_op_geom_010_dg + 61);

    auto tr_0_0_x_yz_xxxz = pbuffer.data(idx_op_geom_010_dg + 62);

    auto tr_0_0_x_yz_xxyy = pbuffer.data(idx_op_geom_010_dg + 63);

    auto tr_0_0_x_yz_xxyz = pbuffer.data(idx_op_geom_010_dg + 64);

    auto tr_0_0_x_yz_xxzz = pbuffer.data(idx_op_geom_010_dg + 65);

    auto tr_0_0_x_yz_xyyy = pbuffer.data(idx_op_geom_010_dg + 66);

    auto tr_0_0_x_yz_xyyz = pbuffer.data(idx_op_geom_010_dg + 67);

    auto tr_0_0_x_yz_xyzz = pbuffer.data(idx_op_geom_010_dg + 68);

    auto tr_0_0_x_yz_xzzz = pbuffer.data(idx_op_geom_010_dg + 69);

    auto tr_0_0_x_yz_yyyy = pbuffer.data(idx_op_geom_010_dg + 70);

    auto tr_0_0_x_yz_yyyz = pbuffer.data(idx_op_geom_010_dg + 71);

    auto tr_0_0_x_yz_yyzz = pbuffer.data(idx_op_geom_010_dg + 72);

    auto tr_0_0_x_yz_yzzz = pbuffer.data(idx_op_geom_010_dg + 73);

    auto tr_0_0_x_yz_zzzz = pbuffer.data(idx_op_geom_010_dg + 74);

    #pragma omp simd aligned(tr_0_0_x_yz_xxxx, tr_0_0_x_yz_xxxy, tr_0_0_x_yz_xxxz, tr_0_0_x_yz_xxyy, tr_0_0_x_yz_xxyz, tr_0_0_x_yz_xxzz, tr_0_0_x_yz_xyyy, tr_0_0_x_yz_xyyz, tr_0_0_x_yz_xyzz, tr_0_0_x_yz_xzzz, tr_0_0_x_yz_yyyy, tr_0_0_x_yz_yyyz, tr_0_0_x_yz_yyzz, tr_0_0_x_yz_yzzz, tr_0_0_x_yz_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yz_xxxx[i] = 2.0 * tr_xyz_xxxx[i] * tbe_0 + 2.0 * tr_yz_xxxxx[i] * tke_0 - 4.0 * tr_yz_xxx[i];

        tr_0_0_x_yz_xxxy[i] = 2.0 * tr_xyz_xxxy[i] * tbe_0 + 2.0 * tr_yz_xxxxy[i] * tke_0 - 3.0 * tr_yz_xxy[i];

        tr_0_0_x_yz_xxxz[i] = 2.0 * tr_xyz_xxxz[i] * tbe_0 + 2.0 * tr_yz_xxxxz[i] * tke_0 - 3.0 * tr_yz_xxz[i];

        tr_0_0_x_yz_xxyy[i] = 2.0 * tr_xyz_xxyy[i] * tbe_0 + 2.0 * tr_yz_xxxyy[i] * tke_0 - 2.0 * tr_yz_xyy[i];

        tr_0_0_x_yz_xxyz[i] = 2.0 * tr_xyz_xxyz[i] * tbe_0 + 2.0 * tr_yz_xxxyz[i] * tke_0 - 2.0 * tr_yz_xyz[i];

        tr_0_0_x_yz_xxzz[i] = 2.0 * tr_xyz_xxzz[i] * tbe_0 + 2.0 * tr_yz_xxxzz[i] * tke_0 - 2.0 * tr_yz_xzz[i];

        tr_0_0_x_yz_xyyy[i] = 2.0 * tr_xyz_xyyy[i] * tbe_0 + 2.0 * tr_yz_xxyyy[i] * tke_0 - tr_yz_yyy[i];

        tr_0_0_x_yz_xyyz[i] = 2.0 * tr_xyz_xyyz[i] * tbe_0 + 2.0 * tr_yz_xxyyz[i] * tke_0 - tr_yz_yyz[i];

        tr_0_0_x_yz_xyzz[i] = 2.0 * tr_xyz_xyzz[i] * tbe_0 + 2.0 * tr_yz_xxyzz[i] * tke_0 - tr_yz_yzz[i];

        tr_0_0_x_yz_xzzz[i] = 2.0 * tr_xyz_xzzz[i] * tbe_0 + 2.0 * tr_yz_xxzzz[i] * tke_0 - tr_yz_zzz[i];

        tr_0_0_x_yz_yyyy[i] = 2.0 * tr_xyz_yyyy[i] * tbe_0 + 2.0 * tr_yz_xyyyy[i] * tke_0;

        tr_0_0_x_yz_yyyz[i] = 2.0 * tr_xyz_yyyz[i] * tbe_0 + 2.0 * tr_yz_xyyyz[i] * tke_0;

        tr_0_0_x_yz_yyzz[i] = 2.0 * tr_xyz_yyzz[i] * tbe_0 + 2.0 * tr_yz_xyyzz[i] * tke_0;

        tr_0_0_x_yz_yzzz[i] = 2.0 * tr_xyz_yzzz[i] * tbe_0 + 2.0 * tr_yz_xyzzz[i] * tke_0;

        tr_0_0_x_yz_zzzz[i] = 2.0 * tr_xyz_zzzz[i] * tbe_0 + 2.0 * tr_yz_xzzzz[i] * tke_0;
    }

    // Set up 75-90 components of targeted buffer : DG

    auto tr_0_0_x_zz_xxxx = pbuffer.data(idx_op_geom_010_dg + 75);

    auto tr_0_0_x_zz_xxxy = pbuffer.data(idx_op_geom_010_dg + 76);

    auto tr_0_0_x_zz_xxxz = pbuffer.data(idx_op_geom_010_dg + 77);

    auto tr_0_0_x_zz_xxyy = pbuffer.data(idx_op_geom_010_dg + 78);

    auto tr_0_0_x_zz_xxyz = pbuffer.data(idx_op_geom_010_dg + 79);

    auto tr_0_0_x_zz_xxzz = pbuffer.data(idx_op_geom_010_dg + 80);

    auto tr_0_0_x_zz_xyyy = pbuffer.data(idx_op_geom_010_dg + 81);

    auto tr_0_0_x_zz_xyyz = pbuffer.data(idx_op_geom_010_dg + 82);

    auto tr_0_0_x_zz_xyzz = pbuffer.data(idx_op_geom_010_dg + 83);

    auto tr_0_0_x_zz_xzzz = pbuffer.data(idx_op_geom_010_dg + 84);

    auto tr_0_0_x_zz_yyyy = pbuffer.data(idx_op_geom_010_dg + 85);

    auto tr_0_0_x_zz_yyyz = pbuffer.data(idx_op_geom_010_dg + 86);

    auto tr_0_0_x_zz_yyzz = pbuffer.data(idx_op_geom_010_dg + 87);

    auto tr_0_0_x_zz_yzzz = pbuffer.data(idx_op_geom_010_dg + 88);

    auto tr_0_0_x_zz_zzzz = pbuffer.data(idx_op_geom_010_dg + 89);

    #pragma omp simd aligned(tr_0_0_x_zz_xxxx, tr_0_0_x_zz_xxxy, tr_0_0_x_zz_xxxz, tr_0_0_x_zz_xxyy, tr_0_0_x_zz_xxyz, tr_0_0_x_zz_xxzz, tr_0_0_x_zz_xyyy, tr_0_0_x_zz_xyyz, tr_0_0_x_zz_xyzz, tr_0_0_x_zz_xzzz, tr_0_0_x_zz_yyyy, tr_0_0_x_zz_yyyz, tr_0_0_x_zz_yyzz, tr_0_0_x_zz_yzzz, tr_0_0_x_zz_zzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zz_xxxx[i] = 2.0 * tr_xzz_xxxx[i] * tbe_0 + 2.0 * tr_zz_xxxxx[i] * tke_0 - 4.0 * tr_zz_xxx[i];

        tr_0_0_x_zz_xxxy[i] = 2.0 * tr_xzz_xxxy[i] * tbe_0 + 2.0 * tr_zz_xxxxy[i] * tke_0 - 3.0 * tr_zz_xxy[i];

        tr_0_0_x_zz_xxxz[i] = 2.0 * tr_xzz_xxxz[i] * tbe_0 + 2.0 * tr_zz_xxxxz[i] * tke_0 - 3.0 * tr_zz_xxz[i];

        tr_0_0_x_zz_xxyy[i] = 2.0 * tr_xzz_xxyy[i] * tbe_0 + 2.0 * tr_zz_xxxyy[i] * tke_0 - 2.0 * tr_zz_xyy[i];

        tr_0_0_x_zz_xxyz[i] = 2.0 * tr_xzz_xxyz[i] * tbe_0 + 2.0 * tr_zz_xxxyz[i] * tke_0 - 2.0 * tr_zz_xyz[i];

        tr_0_0_x_zz_xxzz[i] = 2.0 * tr_xzz_xxzz[i] * tbe_0 + 2.0 * tr_zz_xxxzz[i] * tke_0 - 2.0 * tr_zz_xzz[i];

        tr_0_0_x_zz_xyyy[i] = 2.0 * tr_xzz_xyyy[i] * tbe_0 + 2.0 * tr_zz_xxyyy[i] * tke_0 - tr_zz_yyy[i];

        tr_0_0_x_zz_xyyz[i] = 2.0 * tr_xzz_xyyz[i] * tbe_0 + 2.0 * tr_zz_xxyyz[i] * tke_0 - tr_zz_yyz[i];

        tr_0_0_x_zz_xyzz[i] = 2.0 * tr_xzz_xyzz[i] * tbe_0 + 2.0 * tr_zz_xxyzz[i] * tke_0 - tr_zz_yzz[i];

        tr_0_0_x_zz_xzzz[i] = 2.0 * tr_xzz_xzzz[i] * tbe_0 + 2.0 * tr_zz_xxzzz[i] * tke_0 - tr_zz_zzz[i];

        tr_0_0_x_zz_yyyy[i] = 2.0 * tr_xzz_yyyy[i] * tbe_0 + 2.0 * tr_zz_xyyyy[i] * tke_0;

        tr_0_0_x_zz_yyyz[i] = 2.0 * tr_xzz_yyyz[i] * tbe_0 + 2.0 * tr_zz_xyyyz[i] * tke_0;

        tr_0_0_x_zz_yyzz[i] = 2.0 * tr_xzz_yyzz[i] * tbe_0 + 2.0 * tr_zz_xyyzz[i] * tke_0;

        tr_0_0_x_zz_yzzz[i] = 2.0 * tr_xzz_yzzz[i] * tbe_0 + 2.0 * tr_zz_xyzzz[i] * tke_0;

        tr_0_0_x_zz_zzzz[i] = 2.0 * tr_xzz_zzzz[i] * tbe_0 + 2.0 * tr_zz_xzzzz[i] * tke_0;
    }

    // Set up 90-105 components of targeted buffer : DG

    auto tr_0_0_y_xx_xxxx = pbuffer.data(idx_op_geom_010_dg + 90);

    auto tr_0_0_y_xx_xxxy = pbuffer.data(idx_op_geom_010_dg + 91);

    auto tr_0_0_y_xx_xxxz = pbuffer.data(idx_op_geom_010_dg + 92);

    auto tr_0_0_y_xx_xxyy = pbuffer.data(idx_op_geom_010_dg + 93);

    auto tr_0_0_y_xx_xxyz = pbuffer.data(idx_op_geom_010_dg + 94);

    auto tr_0_0_y_xx_xxzz = pbuffer.data(idx_op_geom_010_dg + 95);

    auto tr_0_0_y_xx_xyyy = pbuffer.data(idx_op_geom_010_dg + 96);

    auto tr_0_0_y_xx_xyyz = pbuffer.data(idx_op_geom_010_dg + 97);

    auto tr_0_0_y_xx_xyzz = pbuffer.data(idx_op_geom_010_dg + 98);

    auto tr_0_0_y_xx_xzzz = pbuffer.data(idx_op_geom_010_dg + 99);

    auto tr_0_0_y_xx_yyyy = pbuffer.data(idx_op_geom_010_dg + 100);

    auto tr_0_0_y_xx_yyyz = pbuffer.data(idx_op_geom_010_dg + 101);

    auto tr_0_0_y_xx_yyzz = pbuffer.data(idx_op_geom_010_dg + 102);

    auto tr_0_0_y_xx_yzzz = pbuffer.data(idx_op_geom_010_dg + 103);

    auto tr_0_0_y_xx_zzzz = pbuffer.data(idx_op_geom_010_dg + 104);

    #pragma omp simd aligned(tr_0_0_y_xx_xxxx, tr_0_0_y_xx_xxxy, tr_0_0_y_xx_xxxz, tr_0_0_y_xx_xxyy, tr_0_0_y_xx_xxyz, tr_0_0_y_xx_xxzz, tr_0_0_y_xx_xyyy, tr_0_0_y_xx_xyyz, tr_0_0_y_xx_xyzz, tr_0_0_y_xx_xzzz, tr_0_0_y_xx_yyyy, tr_0_0_y_xx_yyyz, tr_0_0_y_xx_yyzz, tr_0_0_y_xx_yzzz, tr_0_0_y_xx_zzzz, tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xx_xxxx[i] = 2.0 * tr_xxy_xxxx[i] * tbe_0 + 2.0 * tr_xx_xxxxy[i] * tke_0;

        tr_0_0_y_xx_xxxy[i] = 2.0 * tr_xxy_xxxy[i] * tbe_0 + 2.0 * tr_xx_xxxyy[i] * tke_0 - tr_xx_xxx[i];

        tr_0_0_y_xx_xxxz[i] = 2.0 * tr_xxy_xxxz[i] * tbe_0 + 2.0 * tr_xx_xxxyz[i] * tke_0;

        tr_0_0_y_xx_xxyy[i] = 2.0 * tr_xxy_xxyy[i] * tbe_0 + 2.0 * tr_xx_xxyyy[i] * tke_0 - 2.0 * tr_xx_xxy[i];

        tr_0_0_y_xx_xxyz[i] = 2.0 * tr_xxy_xxyz[i] * tbe_0 + 2.0 * tr_xx_xxyyz[i] * tke_0 - tr_xx_xxz[i];

        tr_0_0_y_xx_xxzz[i] = 2.0 * tr_xxy_xxzz[i] * tbe_0 + 2.0 * tr_xx_xxyzz[i] * tke_0;

        tr_0_0_y_xx_xyyy[i] = 2.0 * tr_xxy_xyyy[i] * tbe_0 + 2.0 * tr_xx_xyyyy[i] * tke_0 - 3.0 * tr_xx_xyy[i];

        tr_0_0_y_xx_xyyz[i] = 2.0 * tr_xxy_xyyz[i] * tbe_0 + 2.0 * tr_xx_xyyyz[i] * tke_0 - 2.0 * tr_xx_xyz[i];

        tr_0_0_y_xx_xyzz[i] = 2.0 * tr_xxy_xyzz[i] * tbe_0 + 2.0 * tr_xx_xyyzz[i] * tke_0 - tr_xx_xzz[i];

        tr_0_0_y_xx_xzzz[i] = 2.0 * tr_xxy_xzzz[i] * tbe_0 + 2.0 * tr_xx_xyzzz[i] * tke_0;

        tr_0_0_y_xx_yyyy[i] = 2.0 * tr_xxy_yyyy[i] * tbe_0 + 2.0 * tr_xx_yyyyy[i] * tke_0 - 4.0 * tr_xx_yyy[i];

        tr_0_0_y_xx_yyyz[i] = 2.0 * tr_xxy_yyyz[i] * tbe_0 + 2.0 * tr_xx_yyyyz[i] * tke_0 - 3.0 * tr_xx_yyz[i];

        tr_0_0_y_xx_yyzz[i] = 2.0 * tr_xxy_yyzz[i] * tbe_0 + 2.0 * tr_xx_yyyzz[i] * tke_0 - 2.0 * tr_xx_yzz[i];

        tr_0_0_y_xx_yzzz[i] = 2.0 * tr_xxy_yzzz[i] * tbe_0 + 2.0 * tr_xx_yyzzz[i] * tke_0 - tr_xx_zzz[i];

        tr_0_0_y_xx_zzzz[i] = 2.0 * tr_xxy_zzzz[i] * tbe_0 + 2.0 * tr_xx_yzzzz[i] * tke_0;
    }

    // Set up 105-120 components of targeted buffer : DG

    auto tr_0_0_y_xy_xxxx = pbuffer.data(idx_op_geom_010_dg + 105);

    auto tr_0_0_y_xy_xxxy = pbuffer.data(idx_op_geom_010_dg + 106);

    auto tr_0_0_y_xy_xxxz = pbuffer.data(idx_op_geom_010_dg + 107);

    auto tr_0_0_y_xy_xxyy = pbuffer.data(idx_op_geom_010_dg + 108);

    auto tr_0_0_y_xy_xxyz = pbuffer.data(idx_op_geom_010_dg + 109);

    auto tr_0_0_y_xy_xxzz = pbuffer.data(idx_op_geom_010_dg + 110);

    auto tr_0_0_y_xy_xyyy = pbuffer.data(idx_op_geom_010_dg + 111);

    auto tr_0_0_y_xy_xyyz = pbuffer.data(idx_op_geom_010_dg + 112);

    auto tr_0_0_y_xy_xyzz = pbuffer.data(idx_op_geom_010_dg + 113);

    auto tr_0_0_y_xy_xzzz = pbuffer.data(idx_op_geom_010_dg + 114);

    auto tr_0_0_y_xy_yyyy = pbuffer.data(idx_op_geom_010_dg + 115);

    auto tr_0_0_y_xy_yyyz = pbuffer.data(idx_op_geom_010_dg + 116);

    auto tr_0_0_y_xy_yyzz = pbuffer.data(idx_op_geom_010_dg + 117);

    auto tr_0_0_y_xy_yzzz = pbuffer.data(idx_op_geom_010_dg + 118);

    auto tr_0_0_y_xy_zzzz = pbuffer.data(idx_op_geom_010_dg + 119);

    #pragma omp simd aligned(tr_0_0_y_xy_xxxx, tr_0_0_y_xy_xxxy, tr_0_0_y_xy_xxxz, tr_0_0_y_xy_xxyy, tr_0_0_y_xy_xxyz, tr_0_0_y_xy_xxzz, tr_0_0_y_xy_xyyy, tr_0_0_y_xy_xyyz, tr_0_0_y_xy_xyzz, tr_0_0_y_xy_xzzz, tr_0_0_y_xy_yyyy, tr_0_0_y_xy_yyyz, tr_0_0_y_xy_yyzz, tr_0_0_y_xy_yzzz, tr_0_0_y_xy_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xy_xxxx[i] = 2.0 * tr_xyy_xxxx[i] * tbe_0 + 2.0 * tr_xy_xxxxy[i] * tke_0 - tr_x_xxxx[i];

        tr_0_0_y_xy_xxxy[i] = 2.0 * tr_xyy_xxxy[i] * tbe_0 + 2.0 * tr_xy_xxxyy[i] * tke_0 - tr_x_xxxy[i] - tr_xy_xxx[i];

        tr_0_0_y_xy_xxxz[i] = 2.0 * tr_xyy_xxxz[i] * tbe_0 + 2.0 * tr_xy_xxxyz[i] * tke_0 - tr_x_xxxz[i];

        tr_0_0_y_xy_xxyy[i] = 2.0 * tr_xyy_xxyy[i] * tbe_0 + 2.0 * tr_xy_xxyyy[i] * tke_0 - tr_x_xxyy[i] - 2.0 * tr_xy_xxy[i];

        tr_0_0_y_xy_xxyz[i] = 2.0 * tr_xyy_xxyz[i] * tbe_0 + 2.0 * tr_xy_xxyyz[i] * tke_0 - tr_x_xxyz[i] - tr_xy_xxz[i];

        tr_0_0_y_xy_xxzz[i] = 2.0 * tr_xyy_xxzz[i] * tbe_0 + 2.0 * tr_xy_xxyzz[i] * tke_0 - tr_x_xxzz[i];

        tr_0_0_y_xy_xyyy[i] = 2.0 * tr_xyy_xyyy[i] * tbe_0 + 2.0 * tr_xy_xyyyy[i] * tke_0 - tr_x_xyyy[i] - 3.0 * tr_xy_xyy[i];

        tr_0_0_y_xy_xyyz[i] = 2.0 * tr_xyy_xyyz[i] * tbe_0 + 2.0 * tr_xy_xyyyz[i] * tke_0 - tr_x_xyyz[i] - 2.0 * tr_xy_xyz[i];

        tr_0_0_y_xy_xyzz[i] = 2.0 * tr_xyy_xyzz[i] * tbe_0 + 2.0 * tr_xy_xyyzz[i] * tke_0 - tr_x_xyzz[i] - tr_xy_xzz[i];

        tr_0_0_y_xy_xzzz[i] = 2.0 * tr_xyy_xzzz[i] * tbe_0 + 2.0 * tr_xy_xyzzz[i] * tke_0 - tr_x_xzzz[i];

        tr_0_0_y_xy_yyyy[i] = 2.0 * tr_xyy_yyyy[i] * tbe_0 + 2.0 * tr_xy_yyyyy[i] * tke_0 - tr_x_yyyy[i] - 4.0 * tr_xy_yyy[i];

        tr_0_0_y_xy_yyyz[i] = 2.0 * tr_xyy_yyyz[i] * tbe_0 + 2.0 * tr_xy_yyyyz[i] * tke_0 - tr_x_yyyz[i] - 3.0 * tr_xy_yyz[i];

        tr_0_0_y_xy_yyzz[i] = 2.0 * tr_xyy_yyzz[i] * tbe_0 + 2.0 * tr_xy_yyyzz[i] * tke_0 - tr_x_yyzz[i] - 2.0 * tr_xy_yzz[i];

        tr_0_0_y_xy_yzzz[i] = 2.0 * tr_xyy_yzzz[i] * tbe_0 + 2.0 * tr_xy_yyzzz[i] * tke_0 - tr_x_yzzz[i] - tr_xy_zzz[i];

        tr_0_0_y_xy_zzzz[i] = 2.0 * tr_xyy_zzzz[i] * tbe_0 + 2.0 * tr_xy_yzzzz[i] * tke_0 - tr_x_zzzz[i];
    }

    // Set up 120-135 components of targeted buffer : DG

    auto tr_0_0_y_xz_xxxx = pbuffer.data(idx_op_geom_010_dg + 120);

    auto tr_0_0_y_xz_xxxy = pbuffer.data(idx_op_geom_010_dg + 121);

    auto tr_0_0_y_xz_xxxz = pbuffer.data(idx_op_geom_010_dg + 122);

    auto tr_0_0_y_xz_xxyy = pbuffer.data(idx_op_geom_010_dg + 123);

    auto tr_0_0_y_xz_xxyz = pbuffer.data(idx_op_geom_010_dg + 124);

    auto tr_0_0_y_xz_xxzz = pbuffer.data(idx_op_geom_010_dg + 125);

    auto tr_0_0_y_xz_xyyy = pbuffer.data(idx_op_geom_010_dg + 126);

    auto tr_0_0_y_xz_xyyz = pbuffer.data(idx_op_geom_010_dg + 127);

    auto tr_0_0_y_xz_xyzz = pbuffer.data(idx_op_geom_010_dg + 128);

    auto tr_0_0_y_xz_xzzz = pbuffer.data(idx_op_geom_010_dg + 129);

    auto tr_0_0_y_xz_yyyy = pbuffer.data(idx_op_geom_010_dg + 130);

    auto tr_0_0_y_xz_yyyz = pbuffer.data(idx_op_geom_010_dg + 131);

    auto tr_0_0_y_xz_yyzz = pbuffer.data(idx_op_geom_010_dg + 132);

    auto tr_0_0_y_xz_yzzz = pbuffer.data(idx_op_geom_010_dg + 133);

    auto tr_0_0_y_xz_zzzz = pbuffer.data(idx_op_geom_010_dg + 134);

    #pragma omp simd aligned(tr_0_0_y_xz_xxxx, tr_0_0_y_xz_xxxy, tr_0_0_y_xz_xxxz, tr_0_0_y_xz_xxyy, tr_0_0_y_xz_xxyz, tr_0_0_y_xz_xxzz, tr_0_0_y_xz_xyyy, tr_0_0_y_xz_xyyz, tr_0_0_y_xz_xyzz, tr_0_0_y_xz_xzzz, tr_0_0_y_xz_yyyy, tr_0_0_y_xz_yyyz, tr_0_0_y_xz_yyzz, tr_0_0_y_xz_yzzz, tr_0_0_y_xz_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xz_xxxx[i] = 2.0 * tr_xyz_xxxx[i] * tbe_0 + 2.0 * tr_xz_xxxxy[i] * tke_0;

        tr_0_0_y_xz_xxxy[i] = 2.0 * tr_xyz_xxxy[i] * tbe_0 + 2.0 * tr_xz_xxxyy[i] * tke_0 - tr_xz_xxx[i];

        tr_0_0_y_xz_xxxz[i] = 2.0 * tr_xyz_xxxz[i] * tbe_0 + 2.0 * tr_xz_xxxyz[i] * tke_0;

        tr_0_0_y_xz_xxyy[i] = 2.0 * tr_xyz_xxyy[i] * tbe_0 + 2.0 * tr_xz_xxyyy[i] * tke_0 - 2.0 * tr_xz_xxy[i];

        tr_0_0_y_xz_xxyz[i] = 2.0 * tr_xyz_xxyz[i] * tbe_0 + 2.0 * tr_xz_xxyyz[i] * tke_0 - tr_xz_xxz[i];

        tr_0_0_y_xz_xxzz[i] = 2.0 * tr_xyz_xxzz[i] * tbe_0 + 2.0 * tr_xz_xxyzz[i] * tke_0;

        tr_0_0_y_xz_xyyy[i] = 2.0 * tr_xyz_xyyy[i] * tbe_0 + 2.0 * tr_xz_xyyyy[i] * tke_0 - 3.0 * tr_xz_xyy[i];

        tr_0_0_y_xz_xyyz[i] = 2.0 * tr_xyz_xyyz[i] * tbe_0 + 2.0 * tr_xz_xyyyz[i] * tke_0 - 2.0 * tr_xz_xyz[i];

        tr_0_0_y_xz_xyzz[i] = 2.0 * tr_xyz_xyzz[i] * tbe_0 + 2.0 * tr_xz_xyyzz[i] * tke_0 - tr_xz_xzz[i];

        tr_0_0_y_xz_xzzz[i] = 2.0 * tr_xyz_xzzz[i] * tbe_0 + 2.0 * tr_xz_xyzzz[i] * tke_0;

        tr_0_0_y_xz_yyyy[i] = 2.0 * tr_xyz_yyyy[i] * tbe_0 + 2.0 * tr_xz_yyyyy[i] * tke_0 - 4.0 * tr_xz_yyy[i];

        tr_0_0_y_xz_yyyz[i] = 2.0 * tr_xyz_yyyz[i] * tbe_0 + 2.0 * tr_xz_yyyyz[i] * tke_0 - 3.0 * tr_xz_yyz[i];

        tr_0_0_y_xz_yyzz[i] = 2.0 * tr_xyz_yyzz[i] * tbe_0 + 2.0 * tr_xz_yyyzz[i] * tke_0 - 2.0 * tr_xz_yzz[i];

        tr_0_0_y_xz_yzzz[i] = 2.0 * tr_xyz_yzzz[i] * tbe_0 + 2.0 * tr_xz_yyzzz[i] * tke_0 - tr_xz_zzz[i];

        tr_0_0_y_xz_zzzz[i] = 2.0 * tr_xyz_zzzz[i] * tbe_0 + 2.0 * tr_xz_yzzzz[i] * tke_0;
    }

    // Set up 135-150 components of targeted buffer : DG

    auto tr_0_0_y_yy_xxxx = pbuffer.data(idx_op_geom_010_dg + 135);

    auto tr_0_0_y_yy_xxxy = pbuffer.data(idx_op_geom_010_dg + 136);

    auto tr_0_0_y_yy_xxxz = pbuffer.data(idx_op_geom_010_dg + 137);

    auto tr_0_0_y_yy_xxyy = pbuffer.data(idx_op_geom_010_dg + 138);

    auto tr_0_0_y_yy_xxyz = pbuffer.data(idx_op_geom_010_dg + 139);

    auto tr_0_0_y_yy_xxzz = pbuffer.data(idx_op_geom_010_dg + 140);

    auto tr_0_0_y_yy_xyyy = pbuffer.data(idx_op_geom_010_dg + 141);

    auto tr_0_0_y_yy_xyyz = pbuffer.data(idx_op_geom_010_dg + 142);

    auto tr_0_0_y_yy_xyzz = pbuffer.data(idx_op_geom_010_dg + 143);

    auto tr_0_0_y_yy_xzzz = pbuffer.data(idx_op_geom_010_dg + 144);

    auto tr_0_0_y_yy_yyyy = pbuffer.data(idx_op_geom_010_dg + 145);

    auto tr_0_0_y_yy_yyyz = pbuffer.data(idx_op_geom_010_dg + 146);

    auto tr_0_0_y_yy_yyzz = pbuffer.data(idx_op_geom_010_dg + 147);

    auto tr_0_0_y_yy_yzzz = pbuffer.data(idx_op_geom_010_dg + 148);

    auto tr_0_0_y_yy_zzzz = pbuffer.data(idx_op_geom_010_dg + 149);

    #pragma omp simd aligned(tr_0_0_y_yy_xxxx, tr_0_0_y_yy_xxxy, tr_0_0_y_yy_xxxz, tr_0_0_y_yy_xxyy, tr_0_0_y_yy_xxyz, tr_0_0_y_yy_xxzz, tr_0_0_y_yy_xyyy, tr_0_0_y_yy_xyyz, tr_0_0_y_yy_xyzz, tr_0_0_y_yy_xzzz, tr_0_0_y_yy_yyyy, tr_0_0_y_yy_yyyz, tr_0_0_y_yy_yyzz, tr_0_0_y_yy_yzzz, tr_0_0_y_yy_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yy_xxxx[i] = 2.0 * tr_yyy_xxxx[i] * tbe_0 + 2.0 * tr_yy_xxxxy[i] * tke_0 - 2.0 * tr_y_xxxx[i];

        tr_0_0_y_yy_xxxy[i] = 2.0 * tr_yyy_xxxy[i] * tbe_0 + 2.0 * tr_yy_xxxyy[i] * tke_0 - 2.0 * tr_y_xxxy[i] - tr_yy_xxx[i];

        tr_0_0_y_yy_xxxz[i] = 2.0 * tr_yyy_xxxz[i] * tbe_0 + 2.0 * tr_yy_xxxyz[i] * tke_0 - 2.0 * tr_y_xxxz[i];

        tr_0_0_y_yy_xxyy[i] = 2.0 * tr_yyy_xxyy[i] * tbe_0 + 2.0 * tr_yy_xxyyy[i] * tke_0 - 2.0 * tr_y_xxyy[i] - 2.0 * tr_yy_xxy[i];

        tr_0_0_y_yy_xxyz[i] = 2.0 * tr_yyy_xxyz[i] * tbe_0 + 2.0 * tr_yy_xxyyz[i] * tke_0 - 2.0 * tr_y_xxyz[i] - tr_yy_xxz[i];

        tr_0_0_y_yy_xxzz[i] = 2.0 * tr_yyy_xxzz[i] * tbe_0 + 2.0 * tr_yy_xxyzz[i] * tke_0 - 2.0 * tr_y_xxzz[i];

        tr_0_0_y_yy_xyyy[i] = 2.0 * tr_yyy_xyyy[i] * tbe_0 + 2.0 * tr_yy_xyyyy[i] * tke_0 - 2.0 * tr_y_xyyy[i] - 3.0 * tr_yy_xyy[i];

        tr_0_0_y_yy_xyyz[i] = 2.0 * tr_yyy_xyyz[i] * tbe_0 + 2.0 * tr_yy_xyyyz[i] * tke_0 - 2.0 * tr_y_xyyz[i] - 2.0 * tr_yy_xyz[i];

        tr_0_0_y_yy_xyzz[i] = 2.0 * tr_yyy_xyzz[i] * tbe_0 + 2.0 * tr_yy_xyyzz[i] * tke_0 - 2.0 * tr_y_xyzz[i] - tr_yy_xzz[i];

        tr_0_0_y_yy_xzzz[i] = 2.0 * tr_yyy_xzzz[i] * tbe_0 + 2.0 * tr_yy_xyzzz[i] * tke_0 - 2.0 * tr_y_xzzz[i];

        tr_0_0_y_yy_yyyy[i] = 2.0 * tr_yyy_yyyy[i] * tbe_0 + 2.0 * tr_yy_yyyyy[i] * tke_0 - 2.0 * tr_y_yyyy[i] - 4.0 * tr_yy_yyy[i];

        tr_0_0_y_yy_yyyz[i] = 2.0 * tr_yyy_yyyz[i] * tbe_0 + 2.0 * tr_yy_yyyyz[i] * tke_0 - 2.0 * tr_y_yyyz[i] - 3.0 * tr_yy_yyz[i];

        tr_0_0_y_yy_yyzz[i] = 2.0 * tr_yyy_yyzz[i] * tbe_0 + 2.0 * tr_yy_yyyzz[i] * tke_0 - 2.0 * tr_y_yyzz[i] - 2.0 * tr_yy_yzz[i];

        tr_0_0_y_yy_yzzz[i] = 2.0 * tr_yyy_yzzz[i] * tbe_0 + 2.0 * tr_yy_yyzzz[i] * tke_0 - 2.0 * tr_y_yzzz[i] - tr_yy_zzz[i];

        tr_0_0_y_yy_zzzz[i] = 2.0 * tr_yyy_zzzz[i] * tbe_0 + 2.0 * tr_yy_yzzzz[i] * tke_0 - 2.0 * tr_y_zzzz[i];
    }

    // Set up 150-165 components of targeted buffer : DG

    auto tr_0_0_y_yz_xxxx = pbuffer.data(idx_op_geom_010_dg + 150);

    auto tr_0_0_y_yz_xxxy = pbuffer.data(idx_op_geom_010_dg + 151);

    auto tr_0_0_y_yz_xxxz = pbuffer.data(idx_op_geom_010_dg + 152);

    auto tr_0_0_y_yz_xxyy = pbuffer.data(idx_op_geom_010_dg + 153);

    auto tr_0_0_y_yz_xxyz = pbuffer.data(idx_op_geom_010_dg + 154);

    auto tr_0_0_y_yz_xxzz = pbuffer.data(idx_op_geom_010_dg + 155);

    auto tr_0_0_y_yz_xyyy = pbuffer.data(idx_op_geom_010_dg + 156);

    auto tr_0_0_y_yz_xyyz = pbuffer.data(idx_op_geom_010_dg + 157);

    auto tr_0_0_y_yz_xyzz = pbuffer.data(idx_op_geom_010_dg + 158);

    auto tr_0_0_y_yz_xzzz = pbuffer.data(idx_op_geom_010_dg + 159);

    auto tr_0_0_y_yz_yyyy = pbuffer.data(idx_op_geom_010_dg + 160);

    auto tr_0_0_y_yz_yyyz = pbuffer.data(idx_op_geom_010_dg + 161);

    auto tr_0_0_y_yz_yyzz = pbuffer.data(idx_op_geom_010_dg + 162);

    auto tr_0_0_y_yz_yzzz = pbuffer.data(idx_op_geom_010_dg + 163);

    auto tr_0_0_y_yz_zzzz = pbuffer.data(idx_op_geom_010_dg + 164);

    #pragma omp simd aligned(tr_0_0_y_yz_xxxx, tr_0_0_y_yz_xxxy, tr_0_0_y_yz_xxxz, tr_0_0_y_yz_xxyy, tr_0_0_y_yz_xxyz, tr_0_0_y_yz_xxzz, tr_0_0_y_yz_xyyy, tr_0_0_y_yz_xyyz, tr_0_0_y_yz_xyzz, tr_0_0_y_yz_xzzz, tr_0_0_y_yz_yyyy, tr_0_0_y_yz_yyyz, tr_0_0_y_yz_yyzz, tr_0_0_y_yz_yzzz, tr_0_0_y_yz_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yz_xxxx[i] = 2.0 * tr_yyz_xxxx[i] * tbe_0 + 2.0 * tr_yz_xxxxy[i] * tke_0 - tr_z_xxxx[i];

        tr_0_0_y_yz_xxxy[i] = 2.0 * tr_yyz_xxxy[i] * tbe_0 + 2.0 * tr_yz_xxxyy[i] * tke_0 - tr_z_xxxy[i] - tr_yz_xxx[i];

        tr_0_0_y_yz_xxxz[i] = 2.0 * tr_yyz_xxxz[i] * tbe_0 + 2.0 * tr_yz_xxxyz[i] * tke_0 - tr_z_xxxz[i];

        tr_0_0_y_yz_xxyy[i] = 2.0 * tr_yyz_xxyy[i] * tbe_0 + 2.0 * tr_yz_xxyyy[i] * tke_0 - tr_z_xxyy[i] - 2.0 * tr_yz_xxy[i];

        tr_0_0_y_yz_xxyz[i] = 2.0 * tr_yyz_xxyz[i] * tbe_0 + 2.0 * tr_yz_xxyyz[i] * tke_0 - tr_z_xxyz[i] - tr_yz_xxz[i];

        tr_0_0_y_yz_xxzz[i] = 2.0 * tr_yyz_xxzz[i] * tbe_0 + 2.0 * tr_yz_xxyzz[i] * tke_0 - tr_z_xxzz[i];

        tr_0_0_y_yz_xyyy[i] = 2.0 * tr_yyz_xyyy[i] * tbe_0 + 2.0 * tr_yz_xyyyy[i] * tke_0 - tr_z_xyyy[i] - 3.0 * tr_yz_xyy[i];

        tr_0_0_y_yz_xyyz[i] = 2.0 * tr_yyz_xyyz[i] * tbe_0 + 2.0 * tr_yz_xyyyz[i] * tke_0 - tr_z_xyyz[i] - 2.0 * tr_yz_xyz[i];

        tr_0_0_y_yz_xyzz[i] = 2.0 * tr_yyz_xyzz[i] * tbe_0 + 2.0 * tr_yz_xyyzz[i] * tke_0 - tr_z_xyzz[i] - tr_yz_xzz[i];

        tr_0_0_y_yz_xzzz[i] = 2.0 * tr_yyz_xzzz[i] * tbe_0 + 2.0 * tr_yz_xyzzz[i] * tke_0 - tr_z_xzzz[i];

        tr_0_0_y_yz_yyyy[i] = 2.0 * tr_yyz_yyyy[i] * tbe_0 + 2.0 * tr_yz_yyyyy[i] * tke_0 - tr_z_yyyy[i] - 4.0 * tr_yz_yyy[i];

        tr_0_0_y_yz_yyyz[i] = 2.0 * tr_yyz_yyyz[i] * tbe_0 + 2.0 * tr_yz_yyyyz[i] * tke_0 - tr_z_yyyz[i] - 3.0 * tr_yz_yyz[i];

        tr_0_0_y_yz_yyzz[i] = 2.0 * tr_yyz_yyzz[i] * tbe_0 + 2.0 * tr_yz_yyyzz[i] * tke_0 - tr_z_yyzz[i] - 2.0 * tr_yz_yzz[i];

        tr_0_0_y_yz_yzzz[i] = 2.0 * tr_yyz_yzzz[i] * tbe_0 + 2.0 * tr_yz_yyzzz[i] * tke_0 - tr_z_yzzz[i] - tr_yz_zzz[i];

        tr_0_0_y_yz_zzzz[i] = 2.0 * tr_yyz_zzzz[i] * tbe_0 + 2.0 * tr_yz_yzzzz[i] * tke_0 - tr_z_zzzz[i];
    }

    // Set up 165-180 components of targeted buffer : DG

    auto tr_0_0_y_zz_xxxx = pbuffer.data(idx_op_geom_010_dg + 165);

    auto tr_0_0_y_zz_xxxy = pbuffer.data(idx_op_geom_010_dg + 166);

    auto tr_0_0_y_zz_xxxz = pbuffer.data(idx_op_geom_010_dg + 167);

    auto tr_0_0_y_zz_xxyy = pbuffer.data(idx_op_geom_010_dg + 168);

    auto tr_0_0_y_zz_xxyz = pbuffer.data(idx_op_geom_010_dg + 169);

    auto tr_0_0_y_zz_xxzz = pbuffer.data(idx_op_geom_010_dg + 170);

    auto tr_0_0_y_zz_xyyy = pbuffer.data(idx_op_geom_010_dg + 171);

    auto tr_0_0_y_zz_xyyz = pbuffer.data(idx_op_geom_010_dg + 172);

    auto tr_0_0_y_zz_xyzz = pbuffer.data(idx_op_geom_010_dg + 173);

    auto tr_0_0_y_zz_xzzz = pbuffer.data(idx_op_geom_010_dg + 174);

    auto tr_0_0_y_zz_yyyy = pbuffer.data(idx_op_geom_010_dg + 175);

    auto tr_0_0_y_zz_yyyz = pbuffer.data(idx_op_geom_010_dg + 176);

    auto tr_0_0_y_zz_yyzz = pbuffer.data(idx_op_geom_010_dg + 177);

    auto tr_0_0_y_zz_yzzz = pbuffer.data(idx_op_geom_010_dg + 178);

    auto tr_0_0_y_zz_zzzz = pbuffer.data(idx_op_geom_010_dg + 179);

    #pragma omp simd aligned(tr_0_0_y_zz_xxxx, tr_0_0_y_zz_xxxy, tr_0_0_y_zz_xxxz, tr_0_0_y_zz_xxyy, tr_0_0_y_zz_xxyz, tr_0_0_y_zz_xxzz, tr_0_0_y_zz_xyyy, tr_0_0_y_zz_xyyz, tr_0_0_y_zz_xyzz, tr_0_0_y_zz_xzzz, tr_0_0_y_zz_yyyy, tr_0_0_y_zz_yyyz, tr_0_0_y_zz_yyzz, tr_0_0_y_zz_yzzz, tr_0_0_y_zz_zzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zz_xxxx[i] = 2.0 * tr_yzz_xxxx[i] * tbe_0 + 2.0 * tr_zz_xxxxy[i] * tke_0;

        tr_0_0_y_zz_xxxy[i] = 2.0 * tr_yzz_xxxy[i] * tbe_0 + 2.0 * tr_zz_xxxyy[i] * tke_0 - tr_zz_xxx[i];

        tr_0_0_y_zz_xxxz[i] = 2.0 * tr_yzz_xxxz[i] * tbe_0 + 2.0 * tr_zz_xxxyz[i] * tke_0;

        tr_0_0_y_zz_xxyy[i] = 2.0 * tr_yzz_xxyy[i] * tbe_0 + 2.0 * tr_zz_xxyyy[i] * tke_0 - 2.0 * tr_zz_xxy[i];

        tr_0_0_y_zz_xxyz[i] = 2.0 * tr_yzz_xxyz[i] * tbe_0 + 2.0 * tr_zz_xxyyz[i] * tke_0 - tr_zz_xxz[i];

        tr_0_0_y_zz_xxzz[i] = 2.0 * tr_yzz_xxzz[i] * tbe_0 + 2.0 * tr_zz_xxyzz[i] * tke_0;

        tr_0_0_y_zz_xyyy[i] = 2.0 * tr_yzz_xyyy[i] * tbe_0 + 2.0 * tr_zz_xyyyy[i] * tke_0 - 3.0 * tr_zz_xyy[i];

        tr_0_0_y_zz_xyyz[i] = 2.0 * tr_yzz_xyyz[i] * tbe_0 + 2.0 * tr_zz_xyyyz[i] * tke_0 - 2.0 * tr_zz_xyz[i];

        tr_0_0_y_zz_xyzz[i] = 2.0 * tr_yzz_xyzz[i] * tbe_0 + 2.0 * tr_zz_xyyzz[i] * tke_0 - tr_zz_xzz[i];

        tr_0_0_y_zz_xzzz[i] = 2.0 * tr_yzz_xzzz[i] * tbe_0 + 2.0 * tr_zz_xyzzz[i] * tke_0;

        tr_0_0_y_zz_yyyy[i] = 2.0 * tr_yzz_yyyy[i] * tbe_0 + 2.0 * tr_zz_yyyyy[i] * tke_0 - 4.0 * tr_zz_yyy[i];

        tr_0_0_y_zz_yyyz[i] = 2.0 * tr_yzz_yyyz[i] * tbe_0 + 2.0 * tr_zz_yyyyz[i] * tke_0 - 3.0 * tr_zz_yyz[i];

        tr_0_0_y_zz_yyzz[i] = 2.0 * tr_yzz_yyzz[i] * tbe_0 + 2.0 * tr_zz_yyyzz[i] * tke_0 - 2.0 * tr_zz_yzz[i];

        tr_0_0_y_zz_yzzz[i] = 2.0 * tr_yzz_yzzz[i] * tbe_0 + 2.0 * tr_zz_yyzzz[i] * tke_0 - tr_zz_zzz[i];

        tr_0_0_y_zz_zzzz[i] = 2.0 * tr_yzz_zzzz[i] * tbe_0 + 2.0 * tr_zz_yzzzz[i] * tke_0;
    }

    // Set up 180-195 components of targeted buffer : DG

    auto tr_0_0_z_xx_xxxx = pbuffer.data(idx_op_geom_010_dg + 180);

    auto tr_0_0_z_xx_xxxy = pbuffer.data(idx_op_geom_010_dg + 181);

    auto tr_0_0_z_xx_xxxz = pbuffer.data(idx_op_geom_010_dg + 182);

    auto tr_0_0_z_xx_xxyy = pbuffer.data(idx_op_geom_010_dg + 183);

    auto tr_0_0_z_xx_xxyz = pbuffer.data(idx_op_geom_010_dg + 184);

    auto tr_0_0_z_xx_xxzz = pbuffer.data(idx_op_geom_010_dg + 185);

    auto tr_0_0_z_xx_xyyy = pbuffer.data(idx_op_geom_010_dg + 186);

    auto tr_0_0_z_xx_xyyz = pbuffer.data(idx_op_geom_010_dg + 187);

    auto tr_0_0_z_xx_xyzz = pbuffer.data(idx_op_geom_010_dg + 188);

    auto tr_0_0_z_xx_xzzz = pbuffer.data(idx_op_geom_010_dg + 189);

    auto tr_0_0_z_xx_yyyy = pbuffer.data(idx_op_geom_010_dg + 190);

    auto tr_0_0_z_xx_yyyz = pbuffer.data(idx_op_geom_010_dg + 191);

    auto tr_0_0_z_xx_yyzz = pbuffer.data(idx_op_geom_010_dg + 192);

    auto tr_0_0_z_xx_yzzz = pbuffer.data(idx_op_geom_010_dg + 193);

    auto tr_0_0_z_xx_zzzz = pbuffer.data(idx_op_geom_010_dg + 194);

    #pragma omp simd aligned(tr_0_0_z_xx_xxxx, tr_0_0_z_xx_xxxy, tr_0_0_z_xx_xxxz, tr_0_0_z_xx_xxyy, tr_0_0_z_xx_xxyz, tr_0_0_z_xx_xxzz, tr_0_0_z_xx_xyyy, tr_0_0_z_xx_xyyz, tr_0_0_z_xx_xyzz, tr_0_0_z_xx_xzzz, tr_0_0_z_xx_yyyy, tr_0_0_z_xx_yyyz, tr_0_0_z_xx_yyzz, tr_0_0_z_xx_yzzz, tr_0_0_z_xx_zzzz, tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xx_zzzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xx_xxxx[i] = 2.0 * tr_xxz_xxxx[i] * tbe_0 + 2.0 * tr_xx_xxxxz[i] * tke_0;

        tr_0_0_z_xx_xxxy[i] = 2.0 * tr_xxz_xxxy[i] * tbe_0 + 2.0 * tr_xx_xxxyz[i] * tke_0;

        tr_0_0_z_xx_xxxz[i] = 2.0 * tr_xxz_xxxz[i] * tbe_0 + 2.0 * tr_xx_xxxzz[i] * tke_0 - tr_xx_xxx[i];

        tr_0_0_z_xx_xxyy[i] = 2.0 * tr_xxz_xxyy[i] * tbe_0 + 2.0 * tr_xx_xxyyz[i] * tke_0;

        tr_0_0_z_xx_xxyz[i] = 2.0 * tr_xxz_xxyz[i] * tbe_0 + 2.0 * tr_xx_xxyzz[i] * tke_0 - tr_xx_xxy[i];

        tr_0_0_z_xx_xxzz[i] = 2.0 * tr_xxz_xxzz[i] * tbe_0 + 2.0 * tr_xx_xxzzz[i] * tke_0 - 2.0 * tr_xx_xxz[i];

        tr_0_0_z_xx_xyyy[i] = 2.0 * tr_xxz_xyyy[i] * tbe_0 + 2.0 * tr_xx_xyyyz[i] * tke_0;

        tr_0_0_z_xx_xyyz[i] = 2.0 * tr_xxz_xyyz[i] * tbe_0 + 2.0 * tr_xx_xyyzz[i] * tke_0 - tr_xx_xyy[i];

        tr_0_0_z_xx_xyzz[i] = 2.0 * tr_xxz_xyzz[i] * tbe_0 + 2.0 * tr_xx_xyzzz[i] * tke_0 - 2.0 * tr_xx_xyz[i];

        tr_0_0_z_xx_xzzz[i] = 2.0 * tr_xxz_xzzz[i] * tbe_0 + 2.0 * tr_xx_xzzzz[i] * tke_0 - 3.0 * tr_xx_xzz[i];

        tr_0_0_z_xx_yyyy[i] = 2.0 * tr_xxz_yyyy[i] * tbe_0 + 2.0 * tr_xx_yyyyz[i] * tke_0;

        tr_0_0_z_xx_yyyz[i] = 2.0 * tr_xxz_yyyz[i] * tbe_0 + 2.0 * tr_xx_yyyzz[i] * tke_0 - tr_xx_yyy[i];

        tr_0_0_z_xx_yyzz[i] = 2.0 * tr_xxz_yyzz[i] * tbe_0 + 2.0 * tr_xx_yyzzz[i] * tke_0 - 2.0 * tr_xx_yyz[i];

        tr_0_0_z_xx_yzzz[i] = 2.0 * tr_xxz_yzzz[i] * tbe_0 + 2.0 * tr_xx_yzzzz[i] * tke_0 - 3.0 * tr_xx_yzz[i];

        tr_0_0_z_xx_zzzz[i] = 2.0 * tr_xxz_zzzz[i] * tbe_0 + 2.0 * tr_xx_zzzzz[i] * tke_0 - 4.0 * tr_xx_zzz[i];
    }

    // Set up 195-210 components of targeted buffer : DG

    auto tr_0_0_z_xy_xxxx = pbuffer.data(idx_op_geom_010_dg + 195);

    auto tr_0_0_z_xy_xxxy = pbuffer.data(idx_op_geom_010_dg + 196);

    auto tr_0_0_z_xy_xxxz = pbuffer.data(idx_op_geom_010_dg + 197);

    auto tr_0_0_z_xy_xxyy = pbuffer.data(idx_op_geom_010_dg + 198);

    auto tr_0_0_z_xy_xxyz = pbuffer.data(idx_op_geom_010_dg + 199);

    auto tr_0_0_z_xy_xxzz = pbuffer.data(idx_op_geom_010_dg + 200);

    auto tr_0_0_z_xy_xyyy = pbuffer.data(idx_op_geom_010_dg + 201);

    auto tr_0_0_z_xy_xyyz = pbuffer.data(idx_op_geom_010_dg + 202);

    auto tr_0_0_z_xy_xyzz = pbuffer.data(idx_op_geom_010_dg + 203);

    auto tr_0_0_z_xy_xzzz = pbuffer.data(idx_op_geom_010_dg + 204);

    auto tr_0_0_z_xy_yyyy = pbuffer.data(idx_op_geom_010_dg + 205);

    auto tr_0_0_z_xy_yyyz = pbuffer.data(idx_op_geom_010_dg + 206);

    auto tr_0_0_z_xy_yyzz = pbuffer.data(idx_op_geom_010_dg + 207);

    auto tr_0_0_z_xy_yzzz = pbuffer.data(idx_op_geom_010_dg + 208);

    auto tr_0_0_z_xy_zzzz = pbuffer.data(idx_op_geom_010_dg + 209);

    #pragma omp simd aligned(tr_0_0_z_xy_xxxx, tr_0_0_z_xy_xxxy, tr_0_0_z_xy_xxxz, tr_0_0_z_xy_xxyy, tr_0_0_z_xy_xxyz, tr_0_0_z_xy_xxzz, tr_0_0_z_xy_xyyy, tr_0_0_z_xy_xyyz, tr_0_0_z_xy_xyzz, tr_0_0_z_xy_xzzz, tr_0_0_z_xy_yyyy, tr_0_0_z_xy_yyyz, tr_0_0_z_xy_yyzz, tr_0_0_z_xy_yzzz, tr_0_0_z_xy_zzzz, tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xy_xxxx[i] = 2.0 * tr_xyz_xxxx[i] * tbe_0 + 2.0 * tr_xy_xxxxz[i] * tke_0;

        tr_0_0_z_xy_xxxy[i] = 2.0 * tr_xyz_xxxy[i] * tbe_0 + 2.0 * tr_xy_xxxyz[i] * tke_0;

        tr_0_0_z_xy_xxxz[i] = 2.0 * tr_xyz_xxxz[i] * tbe_0 + 2.0 * tr_xy_xxxzz[i] * tke_0 - tr_xy_xxx[i];

        tr_0_0_z_xy_xxyy[i] = 2.0 * tr_xyz_xxyy[i] * tbe_0 + 2.0 * tr_xy_xxyyz[i] * tke_0;

        tr_0_0_z_xy_xxyz[i] = 2.0 * tr_xyz_xxyz[i] * tbe_0 + 2.0 * tr_xy_xxyzz[i] * tke_0 - tr_xy_xxy[i];

        tr_0_0_z_xy_xxzz[i] = 2.0 * tr_xyz_xxzz[i] * tbe_0 + 2.0 * tr_xy_xxzzz[i] * tke_0 - 2.0 * tr_xy_xxz[i];

        tr_0_0_z_xy_xyyy[i] = 2.0 * tr_xyz_xyyy[i] * tbe_0 + 2.0 * tr_xy_xyyyz[i] * tke_0;

        tr_0_0_z_xy_xyyz[i] = 2.0 * tr_xyz_xyyz[i] * tbe_0 + 2.0 * tr_xy_xyyzz[i] * tke_0 - tr_xy_xyy[i];

        tr_0_0_z_xy_xyzz[i] = 2.0 * tr_xyz_xyzz[i] * tbe_0 + 2.0 * tr_xy_xyzzz[i] * tke_0 - 2.0 * tr_xy_xyz[i];

        tr_0_0_z_xy_xzzz[i] = 2.0 * tr_xyz_xzzz[i] * tbe_0 + 2.0 * tr_xy_xzzzz[i] * tke_0 - 3.0 * tr_xy_xzz[i];

        tr_0_0_z_xy_yyyy[i] = 2.0 * tr_xyz_yyyy[i] * tbe_0 + 2.0 * tr_xy_yyyyz[i] * tke_0;

        tr_0_0_z_xy_yyyz[i] = 2.0 * tr_xyz_yyyz[i] * tbe_0 + 2.0 * tr_xy_yyyzz[i] * tke_0 - tr_xy_yyy[i];

        tr_0_0_z_xy_yyzz[i] = 2.0 * tr_xyz_yyzz[i] * tbe_0 + 2.0 * tr_xy_yyzzz[i] * tke_0 - 2.0 * tr_xy_yyz[i];

        tr_0_0_z_xy_yzzz[i] = 2.0 * tr_xyz_yzzz[i] * tbe_0 + 2.0 * tr_xy_yzzzz[i] * tke_0 - 3.0 * tr_xy_yzz[i];

        tr_0_0_z_xy_zzzz[i] = 2.0 * tr_xyz_zzzz[i] * tbe_0 + 2.0 * tr_xy_zzzzz[i] * tke_0 - 4.0 * tr_xy_zzz[i];
    }

    // Set up 210-225 components of targeted buffer : DG

    auto tr_0_0_z_xz_xxxx = pbuffer.data(idx_op_geom_010_dg + 210);

    auto tr_0_0_z_xz_xxxy = pbuffer.data(idx_op_geom_010_dg + 211);

    auto tr_0_0_z_xz_xxxz = pbuffer.data(idx_op_geom_010_dg + 212);

    auto tr_0_0_z_xz_xxyy = pbuffer.data(idx_op_geom_010_dg + 213);

    auto tr_0_0_z_xz_xxyz = pbuffer.data(idx_op_geom_010_dg + 214);

    auto tr_0_0_z_xz_xxzz = pbuffer.data(idx_op_geom_010_dg + 215);

    auto tr_0_0_z_xz_xyyy = pbuffer.data(idx_op_geom_010_dg + 216);

    auto tr_0_0_z_xz_xyyz = pbuffer.data(idx_op_geom_010_dg + 217);

    auto tr_0_0_z_xz_xyzz = pbuffer.data(idx_op_geom_010_dg + 218);

    auto tr_0_0_z_xz_xzzz = pbuffer.data(idx_op_geom_010_dg + 219);

    auto tr_0_0_z_xz_yyyy = pbuffer.data(idx_op_geom_010_dg + 220);

    auto tr_0_0_z_xz_yyyz = pbuffer.data(idx_op_geom_010_dg + 221);

    auto tr_0_0_z_xz_yyzz = pbuffer.data(idx_op_geom_010_dg + 222);

    auto tr_0_0_z_xz_yzzz = pbuffer.data(idx_op_geom_010_dg + 223);

    auto tr_0_0_z_xz_zzzz = pbuffer.data(idx_op_geom_010_dg + 224);

    #pragma omp simd aligned(tr_0_0_z_xz_xxxx, tr_0_0_z_xz_xxxy, tr_0_0_z_xz_xxxz, tr_0_0_z_xz_xxyy, tr_0_0_z_xz_xxyz, tr_0_0_z_xz_xxzz, tr_0_0_z_xz_xyyy, tr_0_0_z_xz_xyyz, tr_0_0_z_xz_xyzz, tr_0_0_z_xz_xzzz, tr_0_0_z_xz_yyyy, tr_0_0_z_xz_yyyz, tr_0_0_z_xz_yyzz, tr_0_0_z_xz_yzzz, tr_0_0_z_xz_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xz_xxxx[i] = 2.0 * tr_xzz_xxxx[i] * tbe_0 + 2.0 * tr_xz_xxxxz[i] * tke_0 - tr_x_xxxx[i];

        tr_0_0_z_xz_xxxy[i] = 2.0 * tr_xzz_xxxy[i] * tbe_0 + 2.0 * tr_xz_xxxyz[i] * tke_0 - tr_x_xxxy[i];

        tr_0_0_z_xz_xxxz[i] = 2.0 * tr_xzz_xxxz[i] * tbe_0 + 2.0 * tr_xz_xxxzz[i] * tke_0 - tr_x_xxxz[i] - tr_xz_xxx[i];

        tr_0_0_z_xz_xxyy[i] = 2.0 * tr_xzz_xxyy[i] * tbe_0 + 2.0 * tr_xz_xxyyz[i] * tke_0 - tr_x_xxyy[i];

        tr_0_0_z_xz_xxyz[i] = 2.0 * tr_xzz_xxyz[i] * tbe_0 + 2.0 * tr_xz_xxyzz[i] * tke_0 - tr_x_xxyz[i] - tr_xz_xxy[i];

        tr_0_0_z_xz_xxzz[i] = 2.0 * tr_xzz_xxzz[i] * tbe_0 + 2.0 * tr_xz_xxzzz[i] * tke_0 - tr_x_xxzz[i] - 2.0 * tr_xz_xxz[i];

        tr_0_0_z_xz_xyyy[i] = 2.0 * tr_xzz_xyyy[i] * tbe_0 + 2.0 * tr_xz_xyyyz[i] * tke_0 - tr_x_xyyy[i];

        tr_0_0_z_xz_xyyz[i] = 2.0 * tr_xzz_xyyz[i] * tbe_0 + 2.0 * tr_xz_xyyzz[i] * tke_0 - tr_x_xyyz[i] - tr_xz_xyy[i];

        tr_0_0_z_xz_xyzz[i] = 2.0 * tr_xzz_xyzz[i] * tbe_0 + 2.0 * tr_xz_xyzzz[i] * tke_0 - tr_x_xyzz[i] - 2.0 * tr_xz_xyz[i];

        tr_0_0_z_xz_xzzz[i] = 2.0 * tr_xzz_xzzz[i] * tbe_0 + 2.0 * tr_xz_xzzzz[i] * tke_0 - tr_x_xzzz[i] - 3.0 * tr_xz_xzz[i];

        tr_0_0_z_xz_yyyy[i] = 2.0 * tr_xzz_yyyy[i] * tbe_0 + 2.0 * tr_xz_yyyyz[i] * tke_0 - tr_x_yyyy[i];

        tr_0_0_z_xz_yyyz[i] = 2.0 * tr_xzz_yyyz[i] * tbe_0 + 2.0 * tr_xz_yyyzz[i] * tke_0 - tr_x_yyyz[i] - tr_xz_yyy[i];

        tr_0_0_z_xz_yyzz[i] = 2.0 * tr_xzz_yyzz[i] * tbe_0 + 2.0 * tr_xz_yyzzz[i] * tke_0 - tr_x_yyzz[i] - 2.0 * tr_xz_yyz[i];

        tr_0_0_z_xz_yzzz[i] = 2.0 * tr_xzz_yzzz[i] * tbe_0 + 2.0 * tr_xz_yzzzz[i] * tke_0 - tr_x_yzzz[i] - 3.0 * tr_xz_yzz[i];

        tr_0_0_z_xz_zzzz[i] = 2.0 * tr_xzz_zzzz[i] * tbe_0 + 2.0 * tr_xz_zzzzz[i] * tke_0 - tr_x_zzzz[i] - 4.0 * tr_xz_zzz[i];
    }

    // Set up 225-240 components of targeted buffer : DG

    auto tr_0_0_z_yy_xxxx = pbuffer.data(idx_op_geom_010_dg + 225);

    auto tr_0_0_z_yy_xxxy = pbuffer.data(idx_op_geom_010_dg + 226);

    auto tr_0_0_z_yy_xxxz = pbuffer.data(idx_op_geom_010_dg + 227);

    auto tr_0_0_z_yy_xxyy = pbuffer.data(idx_op_geom_010_dg + 228);

    auto tr_0_0_z_yy_xxyz = pbuffer.data(idx_op_geom_010_dg + 229);

    auto tr_0_0_z_yy_xxzz = pbuffer.data(idx_op_geom_010_dg + 230);

    auto tr_0_0_z_yy_xyyy = pbuffer.data(idx_op_geom_010_dg + 231);

    auto tr_0_0_z_yy_xyyz = pbuffer.data(idx_op_geom_010_dg + 232);

    auto tr_0_0_z_yy_xyzz = pbuffer.data(idx_op_geom_010_dg + 233);

    auto tr_0_0_z_yy_xzzz = pbuffer.data(idx_op_geom_010_dg + 234);

    auto tr_0_0_z_yy_yyyy = pbuffer.data(idx_op_geom_010_dg + 235);

    auto tr_0_0_z_yy_yyyz = pbuffer.data(idx_op_geom_010_dg + 236);

    auto tr_0_0_z_yy_yyzz = pbuffer.data(idx_op_geom_010_dg + 237);

    auto tr_0_0_z_yy_yzzz = pbuffer.data(idx_op_geom_010_dg + 238);

    auto tr_0_0_z_yy_zzzz = pbuffer.data(idx_op_geom_010_dg + 239);

    #pragma omp simd aligned(tr_0_0_z_yy_xxxx, tr_0_0_z_yy_xxxy, tr_0_0_z_yy_xxxz, tr_0_0_z_yy_xxyy, tr_0_0_z_yy_xxyz, tr_0_0_z_yy_xxzz, tr_0_0_z_yy_xyyy, tr_0_0_z_yy_xyyz, tr_0_0_z_yy_xyzz, tr_0_0_z_yy_xzzz, tr_0_0_z_yy_yyyy, tr_0_0_z_yy_yyyz, tr_0_0_z_yy_yyzz, tr_0_0_z_yy_yzzz, tr_0_0_z_yy_zzzz, tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yy_zzzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yy_xxxx[i] = 2.0 * tr_yyz_xxxx[i] * tbe_0 + 2.0 * tr_yy_xxxxz[i] * tke_0;

        tr_0_0_z_yy_xxxy[i] = 2.0 * tr_yyz_xxxy[i] * tbe_0 + 2.0 * tr_yy_xxxyz[i] * tke_0;

        tr_0_0_z_yy_xxxz[i] = 2.0 * tr_yyz_xxxz[i] * tbe_0 + 2.0 * tr_yy_xxxzz[i] * tke_0 - tr_yy_xxx[i];

        tr_0_0_z_yy_xxyy[i] = 2.0 * tr_yyz_xxyy[i] * tbe_0 + 2.0 * tr_yy_xxyyz[i] * tke_0;

        tr_0_0_z_yy_xxyz[i] = 2.0 * tr_yyz_xxyz[i] * tbe_0 + 2.0 * tr_yy_xxyzz[i] * tke_0 - tr_yy_xxy[i];

        tr_0_0_z_yy_xxzz[i] = 2.0 * tr_yyz_xxzz[i] * tbe_0 + 2.0 * tr_yy_xxzzz[i] * tke_0 - 2.0 * tr_yy_xxz[i];

        tr_0_0_z_yy_xyyy[i] = 2.0 * tr_yyz_xyyy[i] * tbe_0 + 2.0 * tr_yy_xyyyz[i] * tke_0;

        tr_0_0_z_yy_xyyz[i] = 2.0 * tr_yyz_xyyz[i] * tbe_0 + 2.0 * tr_yy_xyyzz[i] * tke_0 - tr_yy_xyy[i];

        tr_0_0_z_yy_xyzz[i] = 2.0 * tr_yyz_xyzz[i] * tbe_0 + 2.0 * tr_yy_xyzzz[i] * tke_0 - 2.0 * tr_yy_xyz[i];

        tr_0_0_z_yy_xzzz[i] = 2.0 * tr_yyz_xzzz[i] * tbe_0 + 2.0 * tr_yy_xzzzz[i] * tke_0 - 3.0 * tr_yy_xzz[i];

        tr_0_0_z_yy_yyyy[i] = 2.0 * tr_yyz_yyyy[i] * tbe_0 + 2.0 * tr_yy_yyyyz[i] * tke_0;

        tr_0_0_z_yy_yyyz[i] = 2.0 * tr_yyz_yyyz[i] * tbe_0 + 2.0 * tr_yy_yyyzz[i] * tke_0 - tr_yy_yyy[i];

        tr_0_0_z_yy_yyzz[i] = 2.0 * tr_yyz_yyzz[i] * tbe_0 + 2.0 * tr_yy_yyzzz[i] * tke_0 - 2.0 * tr_yy_yyz[i];

        tr_0_0_z_yy_yzzz[i] = 2.0 * tr_yyz_yzzz[i] * tbe_0 + 2.0 * tr_yy_yzzzz[i] * tke_0 - 3.0 * tr_yy_yzz[i];

        tr_0_0_z_yy_zzzz[i] = 2.0 * tr_yyz_zzzz[i] * tbe_0 + 2.0 * tr_yy_zzzzz[i] * tke_0 - 4.0 * tr_yy_zzz[i];
    }

    // Set up 240-255 components of targeted buffer : DG

    auto tr_0_0_z_yz_xxxx = pbuffer.data(idx_op_geom_010_dg + 240);

    auto tr_0_0_z_yz_xxxy = pbuffer.data(idx_op_geom_010_dg + 241);

    auto tr_0_0_z_yz_xxxz = pbuffer.data(idx_op_geom_010_dg + 242);

    auto tr_0_0_z_yz_xxyy = pbuffer.data(idx_op_geom_010_dg + 243);

    auto tr_0_0_z_yz_xxyz = pbuffer.data(idx_op_geom_010_dg + 244);

    auto tr_0_0_z_yz_xxzz = pbuffer.data(idx_op_geom_010_dg + 245);

    auto tr_0_0_z_yz_xyyy = pbuffer.data(idx_op_geom_010_dg + 246);

    auto tr_0_0_z_yz_xyyz = pbuffer.data(idx_op_geom_010_dg + 247);

    auto tr_0_0_z_yz_xyzz = pbuffer.data(idx_op_geom_010_dg + 248);

    auto tr_0_0_z_yz_xzzz = pbuffer.data(idx_op_geom_010_dg + 249);

    auto tr_0_0_z_yz_yyyy = pbuffer.data(idx_op_geom_010_dg + 250);

    auto tr_0_0_z_yz_yyyz = pbuffer.data(idx_op_geom_010_dg + 251);

    auto tr_0_0_z_yz_yyzz = pbuffer.data(idx_op_geom_010_dg + 252);

    auto tr_0_0_z_yz_yzzz = pbuffer.data(idx_op_geom_010_dg + 253);

    auto tr_0_0_z_yz_zzzz = pbuffer.data(idx_op_geom_010_dg + 254);

    #pragma omp simd aligned(tr_0_0_z_yz_xxxx, tr_0_0_z_yz_xxxy, tr_0_0_z_yz_xxxz, tr_0_0_z_yz_xxyy, tr_0_0_z_yz_xxyz, tr_0_0_z_yz_xxzz, tr_0_0_z_yz_xyyy, tr_0_0_z_yz_xyyz, tr_0_0_z_yz_xyzz, tr_0_0_z_yz_xzzz, tr_0_0_z_yz_yyyy, tr_0_0_z_yz_yyyz, tr_0_0_z_yz_yyzz, tr_0_0_z_yz_yzzz, tr_0_0_z_yz_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yz_xxxx[i] = 2.0 * tr_yzz_xxxx[i] * tbe_0 + 2.0 * tr_yz_xxxxz[i] * tke_0 - tr_y_xxxx[i];

        tr_0_0_z_yz_xxxy[i] = 2.0 * tr_yzz_xxxy[i] * tbe_0 + 2.0 * tr_yz_xxxyz[i] * tke_0 - tr_y_xxxy[i];

        tr_0_0_z_yz_xxxz[i] = 2.0 * tr_yzz_xxxz[i] * tbe_0 + 2.0 * tr_yz_xxxzz[i] * tke_0 - tr_y_xxxz[i] - tr_yz_xxx[i];

        tr_0_0_z_yz_xxyy[i] = 2.0 * tr_yzz_xxyy[i] * tbe_0 + 2.0 * tr_yz_xxyyz[i] * tke_0 - tr_y_xxyy[i];

        tr_0_0_z_yz_xxyz[i] = 2.0 * tr_yzz_xxyz[i] * tbe_0 + 2.0 * tr_yz_xxyzz[i] * tke_0 - tr_y_xxyz[i] - tr_yz_xxy[i];

        tr_0_0_z_yz_xxzz[i] = 2.0 * tr_yzz_xxzz[i] * tbe_0 + 2.0 * tr_yz_xxzzz[i] * tke_0 - tr_y_xxzz[i] - 2.0 * tr_yz_xxz[i];

        tr_0_0_z_yz_xyyy[i] = 2.0 * tr_yzz_xyyy[i] * tbe_0 + 2.0 * tr_yz_xyyyz[i] * tke_0 - tr_y_xyyy[i];

        tr_0_0_z_yz_xyyz[i] = 2.0 * tr_yzz_xyyz[i] * tbe_0 + 2.0 * tr_yz_xyyzz[i] * tke_0 - tr_y_xyyz[i] - tr_yz_xyy[i];

        tr_0_0_z_yz_xyzz[i] = 2.0 * tr_yzz_xyzz[i] * tbe_0 + 2.0 * tr_yz_xyzzz[i] * tke_0 - tr_y_xyzz[i] - 2.0 * tr_yz_xyz[i];

        tr_0_0_z_yz_xzzz[i] = 2.0 * tr_yzz_xzzz[i] * tbe_0 + 2.0 * tr_yz_xzzzz[i] * tke_0 - tr_y_xzzz[i] - 3.0 * tr_yz_xzz[i];

        tr_0_0_z_yz_yyyy[i] = 2.0 * tr_yzz_yyyy[i] * tbe_0 + 2.0 * tr_yz_yyyyz[i] * tke_0 - tr_y_yyyy[i];

        tr_0_0_z_yz_yyyz[i] = 2.0 * tr_yzz_yyyz[i] * tbe_0 + 2.0 * tr_yz_yyyzz[i] * tke_0 - tr_y_yyyz[i] - tr_yz_yyy[i];

        tr_0_0_z_yz_yyzz[i] = 2.0 * tr_yzz_yyzz[i] * tbe_0 + 2.0 * tr_yz_yyzzz[i] * tke_0 - tr_y_yyzz[i] - 2.0 * tr_yz_yyz[i];

        tr_0_0_z_yz_yzzz[i] = 2.0 * tr_yzz_yzzz[i] * tbe_0 + 2.0 * tr_yz_yzzzz[i] * tke_0 - tr_y_yzzz[i] - 3.0 * tr_yz_yzz[i];

        tr_0_0_z_yz_zzzz[i] = 2.0 * tr_yzz_zzzz[i] * tbe_0 + 2.0 * tr_yz_zzzzz[i] * tke_0 - tr_y_zzzz[i] - 4.0 * tr_yz_zzz[i];
    }

    // Set up 255-270 components of targeted buffer : DG

    auto tr_0_0_z_zz_xxxx = pbuffer.data(idx_op_geom_010_dg + 255);

    auto tr_0_0_z_zz_xxxy = pbuffer.data(idx_op_geom_010_dg + 256);

    auto tr_0_0_z_zz_xxxz = pbuffer.data(idx_op_geom_010_dg + 257);

    auto tr_0_0_z_zz_xxyy = pbuffer.data(idx_op_geom_010_dg + 258);

    auto tr_0_0_z_zz_xxyz = pbuffer.data(idx_op_geom_010_dg + 259);

    auto tr_0_0_z_zz_xxzz = pbuffer.data(idx_op_geom_010_dg + 260);

    auto tr_0_0_z_zz_xyyy = pbuffer.data(idx_op_geom_010_dg + 261);

    auto tr_0_0_z_zz_xyyz = pbuffer.data(idx_op_geom_010_dg + 262);

    auto tr_0_0_z_zz_xyzz = pbuffer.data(idx_op_geom_010_dg + 263);

    auto tr_0_0_z_zz_xzzz = pbuffer.data(idx_op_geom_010_dg + 264);

    auto tr_0_0_z_zz_yyyy = pbuffer.data(idx_op_geom_010_dg + 265);

    auto tr_0_0_z_zz_yyyz = pbuffer.data(idx_op_geom_010_dg + 266);

    auto tr_0_0_z_zz_yyzz = pbuffer.data(idx_op_geom_010_dg + 267);

    auto tr_0_0_z_zz_yzzz = pbuffer.data(idx_op_geom_010_dg + 268);

    auto tr_0_0_z_zz_zzzz = pbuffer.data(idx_op_geom_010_dg + 269);

    #pragma omp simd aligned(tr_0_0_z_zz_xxxx, tr_0_0_z_zz_xxxy, tr_0_0_z_zz_xxxz, tr_0_0_z_zz_xxyy, tr_0_0_z_zz_xxyz, tr_0_0_z_zz_xxzz, tr_0_0_z_zz_xyyy, tr_0_0_z_zz_xyyz, tr_0_0_z_zz_xyzz, tr_0_0_z_zz_xzzz, tr_0_0_z_zz_yyyy, tr_0_0_z_zz_yyyz, tr_0_0_z_zz_yyzz, tr_0_0_z_zz_yzzz, tr_0_0_z_zz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zz_xxxx[i] = 2.0 * tr_zzz_xxxx[i] * tbe_0 + 2.0 * tr_zz_xxxxz[i] * tke_0 - 2.0 * tr_z_xxxx[i];

        tr_0_0_z_zz_xxxy[i] = 2.0 * tr_zzz_xxxy[i] * tbe_0 + 2.0 * tr_zz_xxxyz[i] * tke_0 - 2.0 * tr_z_xxxy[i];

        tr_0_0_z_zz_xxxz[i] = 2.0 * tr_zzz_xxxz[i] * tbe_0 + 2.0 * tr_zz_xxxzz[i] * tke_0 - 2.0 * tr_z_xxxz[i] - tr_zz_xxx[i];

        tr_0_0_z_zz_xxyy[i] = 2.0 * tr_zzz_xxyy[i] * tbe_0 + 2.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_z_xxyy[i];

        tr_0_0_z_zz_xxyz[i] = 2.0 * tr_zzz_xxyz[i] * tbe_0 + 2.0 * tr_zz_xxyzz[i] * tke_0 - 2.0 * tr_z_xxyz[i] - tr_zz_xxy[i];

        tr_0_0_z_zz_xxzz[i] = 2.0 * tr_zzz_xxzz[i] * tbe_0 + 2.0 * tr_zz_xxzzz[i] * tke_0 - 2.0 * tr_z_xxzz[i] - 2.0 * tr_zz_xxz[i];

        tr_0_0_z_zz_xyyy[i] = 2.0 * tr_zzz_xyyy[i] * tbe_0 + 2.0 * tr_zz_xyyyz[i] * tke_0 - 2.0 * tr_z_xyyy[i];

        tr_0_0_z_zz_xyyz[i] = 2.0 * tr_zzz_xyyz[i] * tbe_0 + 2.0 * tr_zz_xyyzz[i] * tke_0 - 2.0 * tr_z_xyyz[i] - tr_zz_xyy[i];

        tr_0_0_z_zz_xyzz[i] = 2.0 * tr_zzz_xyzz[i] * tbe_0 + 2.0 * tr_zz_xyzzz[i] * tke_0 - 2.0 * tr_z_xyzz[i] - 2.0 * tr_zz_xyz[i];

        tr_0_0_z_zz_xzzz[i] = 2.0 * tr_zzz_xzzz[i] * tbe_0 + 2.0 * tr_zz_xzzzz[i] * tke_0 - 2.0 * tr_z_xzzz[i] - 3.0 * tr_zz_xzz[i];

        tr_0_0_z_zz_yyyy[i] = 2.0 * tr_zzz_yyyy[i] * tbe_0 + 2.0 * tr_zz_yyyyz[i] * tke_0 - 2.0 * tr_z_yyyy[i];

        tr_0_0_z_zz_yyyz[i] = 2.0 * tr_zzz_yyyz[i] * tbe_0 + 2.0 * tr_zz_yyyzz[i] * tke_0 - 2.0 * tr_z_yyyz[i] - tr_zz_yyy[i];

        tr_0_0_z_zz_yyzz[i] = 2.0 * tr_zzz_yyzz[i] * tbe_0 + 2.0 * tr_zz_yyzzz[i] * tke_0 - 2.0 * tr_z_yyzz[i] - 2.0 * tr_zz_yyz[i];

        tr_0_0_z_zz_yzzz[i] = 2.0 * tr_zzz_yzzz[i] * tbe_0 + 2.0 * tr_zz_yzzzz[i] * tke_0 - 2.0 * tr_z_yzzz[i] - 3.0 * tr_zz_yzz[i];

        tr_0_0_z_zz_zzzz[i] = 2.0 * tr_zzz_zzzz[i] * tbe_0 + 2.0 * tr_zz_zzzzz[i] * tke_0 - 2.0 * tr_z_zzzz[i] - 4.0 * tr_zz_zzz[i];
    }

}

} // t2cgeom namespace

