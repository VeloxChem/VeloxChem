#include "GeometricalDerivatives110ForFG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_fg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_fg,
                         const int idx_op_pg,
                         const int idx_op_df,
                         const int idx_op_dh,
                         const int idx_op_fg,
                         const int idx_op_gf,
                         const int idx_op_gh,
                         const int idx_op_hg,
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

    // Set up components of auxiliary buffer : GF

    auto tr_xxxx_xxx = pbuffer.data(idx_op_gf);

    auto tr_xxxx_xxy = pbuffer.data(idx_op_gf + 1);

    auto tr_xxxx_xxz = pbuffer.data(idx_op_gf + 2);

    auto tr_xxxx_xyy = pbuffer.data(idx_op_gf + 3);

    auto tr_xxxx_xyz = pbuffer.data(idx_op_gf + 4);

    auto tr_xxxx_xzz = pbuffer.data(idx_op_gf + 5);

    auto tr_xxxx_yyy = pbuffer.data(idx_op_gf + 6);

    auto tr_xxxx_yyz = pbuffer.data(idx_op_gf + 7);

    auto tr_xxxx_yzz = pbuffer.data(idx_op_gf + 8);

    auto tr_xxxx_zzz = pbuffer.data(idx_op_gf + 9);

    auto tr_xxxy_xxx = pbuffer.data(idx_op_gf + 10);

    auto tr_xxxy_xxy = pbuffer.data(idx_op_gf + 11);

    auto tr_xxxy_xxz = pbuffer.data(idx_op_gf + 12);

    auto tr_xxxy_xyy = pbuffer.data(idx_op_gf + 13);

    auto tr_xxxy_xyz = pbuffer.data(idx_op_gf + 14);

    auto tr_xxxy_xzz = pbuffer.data(idx_op_gf + 15);

    auto tr_xxxy_yyy = pbuffer.data(idx_op_gf + 16);

    auto tr_xxxy_yyz = pbuffer.data(idx_op_gf + 17);

    auto tr_xxxy_yzz = pbuffer.data(idx_op_gf + 18);

    auto tr_xxxy_zzz = pbuffer.data(idx_op_gf + 19);

    auto tr_xxxz_xxx = pbuffer.data(idx_op_gf + 20);

    auto tr_xxxz_xxy = pbuffer.data(idx_op_gf + 21);

    auto tr_xxxz_xxz = pbuffer.data(idx_op_gf + 22);

    auto tr_xxxz_xyy = pbuffer.data(idx_op_gf + 23);

    auto tr_xxxz_xyz = pbuffer.data(idx_op_gf + 24);

    auto tr_xxxz_xzz = pbuffer.data(idx_op_gf + 25);

    auto tr_xxxz_yyy = pbuffer.data(idx_op_gf + 26);

    auto tr_xxxz_yyz = pbuffer.data(idx_op_gf + 27);

    auto tr_xxxz_yzz = pbuffer.data(idx_op_gf + 28);

    auto tr_xxxz_zzz = pbuffer.data(idx_op_gf + 29);

    auto tr_xxyy_xxx = pbuffer.data(idx_op_gf + 30);

    auto tr_xxyy_xxy = pbuffer.data(idx_op_gf + 31);

    auto tr_xxyy_xxz = pbuffer.data(idx_op_gf + 32);

    auto tr_xxyy_xyy = pbuffer.data(idx_op_gf + 33);

    auto tr_xxyy_xyz = pbuffer.data(idx_op_gf + 34);

    auto tr_xxyy_xzz = pbuffer.data(idx_op_gf + 35);

    auto tr_xxyy_yyy = pbuffer.data(idx_op_gf + 36);

    auto tr_xxyy_yyz = pbuffer.data(idx_op_gf + 37);

    auto tr_xxyy_yzz = pbuffer.data(idx_op_gf + 38);

    auto tr_xxyy_zzz = pbuffer.data(idx_op_gf + 39);

    auto tr_xxyz_xxx = pbuffer.data(idx_op_gf + 40);

    auto tr_xxyz_xxy = pbuffer.data(idx_op_gf + 41);

    auto tr_xxyz_xxz = pbuffer.data(idx_op_gf + 42);

    auto tr_xxyz_xyy = pbuffer.data(idx_op_gf + 43);

    auto tr_xxyz_xyz = pbuffer.data(idx_op_gf + 44);

    auto tr_xxyz_xzz = pbuffer.data(idx_op_gf + 45);

    auto tr_xxyz_yyy = pbuffer.data(idx_op_gf + 46);

    auto tr_xxyz_yyz = pbuffer.data(idx_op_gf + 47);

    auto tr_xxyz_yzz = pbuffer.data(idx_op_gf + 48);

    auto tr_xxyz_zzz = pbuffer.data(idx_op_gf + 49);

    auto tr_xxzz_xxx = pbuffer.data(idx_op_gf + 50);

    auto tr_xxzz_xxy = pbuffer.data(idx_op_gf + 51);

    auto tr_xxzz_xxz = pbuffer.data(idx_op_gf + 52);

    auto tr_xxzz_xyy = pbuffer.data(idx_op_gf + 53);

    auto tr_xxzz_xyz = pbuffer.data(idx_op_gf + 54);

    auto tr_xxzz_xzz = pbuffer.data(idx_op_gf + 55);

    auto tr_xxzz_yyy = pbuffer.data(idx_op_gf + 56);

    auto tr_xxzz_yyz = pbuffer.data(idx_op_gf + 57);

    auto tr_xxzz_yzz = pbuffer.data(idx_op_gf + 58);

    auto tr_xxzz_zzz = pbuffer.data(idx_op_gf + 59);

    auto tr_xyyy_xxx = pbuffer.data(idx_op_gf + 60);

    auto tr_xyyy_xxy = pbuffer.data(idx_op_gf + 61);

    auto tr_xyyy_xxz = pbuffer.data(idx_op_gf + 62);

    auto tr_xyyy_xyy = pbuffer.data(idx_op_gf + 63);

    auto tr_xyyy_xyz = pbuffer.data(idx_op_gf + 64);

    auto tr_xyyy_xzz = pbuffer.data(idx_op_gf + 65);

    auto tr_xyyy_yyy = pbuffer.data(idx_op_gf + 66);

    auto tr_xyyy_yyz = pbuffer.data(idx_op_gf + 67);

    auto tr_xyyy_yzz = pbuffer.data(idx_op_gf + 68);

    auto tr_xyyy_zzz = pbuffer.data(idx_op_gf + 69);

    auto tr_xyyz_xxx = pbuffer.data(idx_op_gf + 70);

    auto tr_xyyz_xxy = pbuffer.data(idx_op_gf + 71);

    auto tr_xyyz_xxz = pbuffer.data(idx_op_gf + 72);

    auto tr_xyyz_xyy = pbuffer.data(idx_op_gf + 73);

    auto tr_xyyz_xyz = pbuffer.data(idx_op_gf + 74);

    auto tr_xyyz_xzz = pbuffer.data(idx_op_gf + 75);

    auto tr_xyyz_yyy = pbuffer.data(idx_op_gf + 76);

    auto tr_xyyz_yyz = pbuffer.data(idx_op_gf + 77);

    auto tr_xyyz_yzz = pbuffer.data(idx_op_gf + 78);

    auto tr_xyyz_zzz = pbuffer.data(idx_op_gf + 79);

    auto tr_xyzz_xxx = pbuffer.data(idx_op_gf + 80);

    auto tr_xyzz_xxy = pbuffer.data(idx_op_gf + 81);

    auto tr_xyzz_xxz = pbuffer.data(idx_op_gf + 82);

    auto tr_xyzz_xyy = pbuffer.data(idx_op_gf + 83);

    auto tr_xyzz_xyz = pbuffer.data(idx_op_gf + 84);

    auto tr_xyzz_xzz = pbuffer.data(idx_op_gf + 85);

    auto tr_xyzz_yyy = pbuffer.data(idx_op_gf + 86);

    auto tr_xyzz_yyz = pbuffer.data(idx_op_gf + 87);

    auto tr_xyzz_yzz = pbuffer.data(idx_op_gf + 88);

    auto tr_xyzz_zzz = pbuffer.data(idx_op_gf + 89);

    auto tr_xzzz_xxx = pbuffer.data(idx_op_gf + 90);

    auto tr_xzzz_xxy = pbuffer.data(idx_op_gf + 91);

    auto tr_xzzz_xxz = pbuffer.data(idx_op_gf + 92);

    auto tr_xzzz_xyy = pbuffer.data(idx_op_gf + 93);

    auto tr_xzzz_xyz = pbuffer.data(idx_op_gf + 94);

    auto tr_xzzz_xzz = pbuffer.data(idx_op_gf + 95);

    auto tr_xzzz_yyy = pbuffer.data(idx_op_gf + 96);

    auto tr_xzzz_yyz = pbuffer.data(idx_op_gf + 97);

    auto tr_xzzz_yzz = pbuffer.data(idx_op_gf + 98);

    auto tr_xzzz_zzz = pbuffer.data(idx_op_gf + 99);

    auto tr_yyyy_xxx = pbuffer.data(idx_op_gf + 100);

    auto tr_yyyy_xxy = pbuffer.data(idx_op_gf + 101);

    auto tr_yyyy_xxz = pbuffer.data(idx_op_gf + 102);

    auto tr_yyyy_xyy = pbuffer.data(idx_op_gf + 103);

    auto tr_yyyy_xyz = pbuffer.data(idx_op_gf + 104);

    auto tr_yyyy_xzz = pbuffer.data(idx_op_gf + 105);

    auto tr_yyyy_yyy = pbuffer.data(idx_op_gf + 106);

    auto tr_yyyy_yyz = pbuffer.data(idx_op_gf + 107);

    auto tr_yyyy_yzz = pbuffer.data(idx_op_gf + 108);

    auto tr_yyyy_zzz = pbuffer.data(idx_op_gf + 109);

    auto tr_yyyz_xxx = pbuffer.data(idx_op_gf + 110);

    auto tr_yyyz_xxy = pbuffer.data(idx_op_gf + 111);

    auto tr_yyyz_xxz = pbuffer.data(idx_op_gf + 112);

    auto tr_yyyz_xyy = pbuffer.data(idx_op_gf + 113);

    auto tr_yyyz_xyz = pbuffer.data(idx_op_gf + 114);

    auto tr_yyyz_xzz = pbuffer.data(idx_op_gf + 115);

    auto tr_yyyz_yyy = pbuffer.data(idx_op_gf + 116);

    auto tr_yyyz_yyz = pbuffer.data(idx_op_gf + 117);

    auto tr_yyyz_yzz = pbuffer.data(idx_op_gf + 118);

    auto tr_yyyz_zzz = pbuffer.data(idx_op_gf + 119);

    auto tr_yyzz_xxx = pbuffer.data(idx_op_gf + 120);

    auto tr_yyzz_xxy = pbuffer.data(idx_op_gf + 121);

    auto tr_yyzz_xxz = pbuffer.data(idx_op_gf + 122);

    auto tr_yyzz_xyy = pbuffer.data(idx_op_gf + 123);

    auto tr_yyzz_xyz = pbuffer.data(idx_op_gf + 124);

    auto tr_yyzz_xzz = pbuffer.data(idx_op_gf + 125);

    auto tr_yyzz_yyy = pbuffer.data(idx_op_gf + 126);

    auto tr_yyzz_yyz = pbuffer.data(idx_op_gf + 127);

    auto tr_yyzz_yzz = pbuffer.data(idx_op_gf + 128);

    auto tr_yyzz_zzz = pbuffer.data(idx_op_gf + 129);

    auto tr_yzzz_xxx = pbuffer.data(idx_op_gf + 130);

    auto tr_yzzz_xxy = pbuffer.data(idx_op_gf + 131);

    auto tr_yzzz_xxz = pbuffer.data(idx_op_gf + 132);

    auto tr_yzzz_xyy = pbuffer.data(idx_op_gf + 133);

    auto tr_yzzz_xyz = pbuffer.data(idx_op_gf + 134);

    auto tr_yzzz_xzz = pbuffer.data(idx_op_gf + 135);

    auto tr_yzzz_yyy = pbuffer.data(idx_op_gf + 136);

    auto tr_yzzz_yyz = pbuffer.data(idx_op_gf + 137);

    auto tr_yzzz_yzz = pbuffer.data(idx_op_gf + 138);

    auto tr_yzzz_zzz = pbuffer.data(idx_op_gf + 139);

    auto tr_zzzz_xxx = pbuffer.data(idx_op_gf + 140);

    auto tr_zzzz_xxy = pbuffer.data(idx_op_gf + 141);

    auto tr_zzzz_xxz = pbuffer.data(idx_op_gf + 142);

    auto tr_zzzz_xyy = pbuffer.data(idx_op_gf + 143);

    auto tr_zzzz_xyz = pbuffer.data(idx_op_gf + 144);

    auto tr_zzzz_xzz = pbuffer.data(idx_op_gf + 145);

    auto tr_zzzz_yyy = pbuffer.data(idx_op_gf + 146);

    auto tr_zzzz_yyz = pbuffer.data(idx_op_gf + 147);

    auto tr_zzzz_yzz = pbuffer.data(idx_op_gf + 148);

    auto tr_zzzz_zzz = pbuffer.data(idx_op_gf + 149);

    // Set up components of auxiliary buffer : GH

    auto tr_xxxx_xxxxx = pbuffer.data(idx_op_gh);

    auto tr_xxxx_xxxxy = pbuffer.data(idx_op_gh + 1);

    auto tr_xxxx_xxxxz = pbuffer.data(idx_op_gh + 2);

    auto tr_xxxx_xxxyy = pbuffer.data(idx_op_gh + 3);

    auto tr_xxxx_xxxyz = pbuffer.data(idx_op_gh + 4);

    auto tr_xxxx_xxxzz = pbuffer.data(idx_op_gh + 5);

    auto tr_xxxx_xxyyy = pbuffer.data(idx_op_gh + 6);

    auto tr_xxxx_xxyyz = pbuffer.data(idx_op_gh + 7);

    auto tr_xxxx_xxyzz = pbuffer.data(idx_op_gh + 8);

    auto tr_xxxx_xxzzz = pbuffer.data(idx_op_gh + 9);

    auto tr_xxxx_xyyyy = pbuffer.data(idx_op_gh + 10);

    auto tr_xxxx_xyyyz = pbuffer.data(idx_op_gh + 11);

    auto tr_xxxx_xyyzz = pbuffer.data(idx_op_gh + 12);

    auto tr_xxxx_xyzzz = pbuffer.data(idx_op_gh + 13);

    auto tr_xxxx_xzzzz = pbuffer.data(idx_op_gh + 14);

    auto tr_xxxx_yyyyy = pbuffer.data(idx_op_gh + 15);

    auto tr_xxxx_yyyyz = pbuffer.data(idx_op_gh + 16);

    auto tr_xxxx_yyyzz = pbuffer.data(idx_op_gh + 17);

    auto tr_xxxx_yyzzz = pbuffer.data(idx_op_gh + 18);

    auto tr_xxxx_yzzzz = pbuffer.data(idx_op_gh + 19);

    auto tr_xxxx_zzzzz = pbuffer.data(idx_op_gh + 20);

    auto tr_xxxy_xxxxx = pbuffer.data(idx_op_gh + 21);

    auto tr_xxxy_xxxxy = pbuffer.data(idx_op_gh + 22);

    auto tr_xxxy_xxxxz = pbuffer.data(idx_op_gh + 23);

    auto tr_xxxy_xxxyy = pbuffer.data(idx_op_gh + 24);

    auto tr_xxxy_xxxyz = pbuffer.data(idx_op_gh + 25);

    auto tr_xxxy_xxxzz = pbuffer.data(idx_op_gh + 26);

    auto tr_xxxy_xxyyy = pbuffer.data(idx_op_gh + 27);

    auto tr_xxxy_xxyyz = pbuffer.data(idx_op_gh + 28);

    auto tr_xxxy_xxyzz = pbuffer.data(idx_op_gh + 29);

    auto tr_xxxy_xxzzz = pbuffer.data(idx_op_gh + 30);

    auto tr_xxxy_xyyyy = pbuffer.data(idx_op_gh + 31);

    auto tr_xxxy_xyyyz = pbuffer.data(idx_op_gh + 32);

    auto tr_xxxy_xyyzz = pbuffer.data(idx_op_gh + 33);

    auto tr_xxxy_xyzzz = pbuffer.data(idx_op_gh + 34);

    auto tr_xxxy_xzzzz = pbuffer.data(idx_op_gh + 35);

    auto tr_xxxy_yyyyy = pbuffer.data(idx_op_gh + 36);

    auto tr_xxxy_yyyyz = pbuffer.data(idx_op_gh + 37);

    auto tr_xxxy_yyyzz = pbuffer.data(idx_op_gh + 38);

    auto tr_xxxy_yyzzz = pbuffer.data(idx_op_gh + 39);

    auto tr_xxxy_yzzzz = pbuffer.data(idx_op_gh + 40);

    auto tr_xxxy_zzzzz = pbuffer.data(idx_op_gh + 41);

    auto tr_xxxz_xxxxx = pbuffer.data(idx_op_gh + 42);

    auto tr_xxxz_xxxxy = pbuffer.data(idx_op_gh + 43);

    auto tr_xxxz_xxxxz = pbuffer.data(idx_op_gh + 44);

    auto tr_xxxz_xxxyy = pbuffer.data(idx_op_gh + 45);

    auto tr_xxxz_xxxyz = pbuffer.data(idx_op_gh + 46);

    auto tr_xxxz_xxxzz = pbuffer.data(idx_op_gh + 47);

    auto tr_xxxz_xxyyy = pbuffer.data(idx_op_gh + 48);

    auto tr_xxxz_xxyyz = pbuffer.data(idx_op_gh + 49);

    auto tr_xxxz_xxyzz = pbuffer.data(idx_op_gh + 50);

    auto tr_xxxz_xxzzz = pbuffer.data(idx_op_gh + 51);

    auto tr_xxxz_xyyyy = pbuffer.data(idx_op_gh + 52);

    auto tr_xxxz_xyyyz = pbuffer.data(idx_op_gh + 53);

    auto tr_xxxz_xyyzz = pbuffer.data(idx_op_gh + 54);

    auto tr_xxxz_xyzzz = pbuffer.data(idx_op_gh + 55);

    auto tr_xxxz_xzzzz = pbuffer.data(idx_op_gh + 56);

    auto tr_xxxz_yyyyy = pbuffer.data(idx_op_gh + 57);

    auto tr_xxxz_yyyyz = pbuffer.data(idx_op_gh + 58);

    auto tr_xxxz_yyyzz = pbuffer.data(idx_op_gh + 59);

    auto tr_xxxz_yyzzz = pbuffer.data(idx_op_gh + 60);

    auto tr_xxxz_yzzzz = pbuffer.data(idx_op_gh + 61);

    auto tr_xxxz_zzzzz = pbuffer.data(idx_op_gh + 62);

    auto tr_xxyy_xxxxx = pbuffer.data(idx_op_gh + 63);

    auto tr_xxyy_xxxxy = pbuffer.data(idx_op_gh + 64);

    auto tr_xxyy_xxxxz = pbuffer.data(idx_op_gh + 65);

    auto tr_xxyy_xxxyy = pbuffer.data(idx_op_gh + 66);

    auto tr_xxyy_xxxyz = pbuffer.data(idx_op_gh + 67);

    auto tr_xxyy_xxxzz = pbuffer.data(idx_op_gh + 68);

    auto tr_xxyy_xxyyy = pbuffer.data(idx_op_gh + 69);

    auto tr_xxyy_xxyyz = pbuffer.data(idx_op_gh + 70);

    auto tr_xxyy_xxyzz = pbuffer.data(idx_op_gh + 71);

    auto tr_xxyy_xxzzz = pbuffer.data(idx_op_gh + 72);

    auto tr_xxyy_xyyyy = pbuffer.data(idx_op_gh + 73);

    auto tr_xxyy_xyyyz = pbuffer.data(idx_op_gh + 74);

    auto tr_xxyy_xyyzz = pbuffer.data(idx_op_gh + 75);

    auto tr_xxyy_xyzzz = pbuffer.data(idx_op_gh + 76);

    auto tr_xxyy_xzzzz = pbuffer.data(idx_op_gh + 77);

    auto tr_xxyy_yyyyy = pbuffer.data(idx_op_gh + 78);

    auto tr_xxyy_yyyyz = pbuffer.data(idx_op_gh + 79);

    auto tr_xxyy_yyyzz = pbuffer.data(idx_op_gh + 80);

    auto tr_xxyy_yyzzz = pbuffer.data(idx_op_gh + 81);

    auto tr_xxyy_yzzzz = pbuffer.data(idx_op_gh + 82);

    auto tr_xxyy_zzzzz = pbuffer.data(idx_op_gh + 83);

    auto tr_xxyz_xxxxx = pbuffer.data(idx_op_gh + 84);

    auto tr_xxyz_xxxxy = pbuffer.data(idx_op_gh + 85);

    auto tr_xxyz_xxxxz = pbuffer.data(idx_op_gh + 86);

    auto tr_xxyz_xxxyy = pbuffer.data(idx_op_gh + 87);

    auto tr_xxyz_xxxyz = pbuffer.data(idx_op_gh + 88);

    auto tr_xxyz_xxxzz = pbuffer.data(idx_op_gh + 89);

    auto tr_xxyz_xxyyy = pbuffer.data(idx_op_gh + 90);

    auto tr_xxyz_xxyyz = pbuffer.data(idx_op_gh + 91);

    auto tr_xxyz_xxyzz = pbuffer.data(idx_op_gh + 92);

    auto tr_xxyz_xxzzz = pbuffer.data(idx_op_gh + 93);

    auto tr_xxyz_xyyyy = pbuffer.data(idx_op_gh + 94);

    auto tr_xxyz_xyyyz = pbuffer.data(idx_op_gh + 95);

    auto tr_xxyz_xyyzz = pbuffer.data(idx_op_gh + 96);

    auto tr_xxyz_xyzzz = pbuffer.data(idx_op_gh + 97);

    auto tr_xxyz_xzzzz = pbuffer.data(idx_op_gh + 98);

    auto tr_xxyz_yyyyy = pbuffer.data(idx_op_gh + 99);

    auto tr_xxyz_yyyyz = pbuffer.data(idx_op_gh + 100);

    auto tr_xxyz_yyyzz = pbuffer.data(idx_op_gh + 101);

    auto tr_xxyz_yyzzz = pbuffer.data(idx_op_gh + 102);

    auto tr_xxyz_yzzzz = pbuffer.data(idx_op_gh + 103);

    auto tr_xxyz_zzzzz = pbuffer.data(idx_op_gh + 104);

    auto tr_xxzz_xxxxx = pbuffer.data(idx_op_gh + 105);

    auto tr_xxzz_xxxxy = pbuffer.data(idx_op_gh + 106);

    auto tr_xxzz_xxxxz = pbuffer.data(idx_op_gh + 107);

    auto tr_xxzz_xxxyy = pbuffer.data(idx_op_gh + 108);

    auto tr_xxzz_xxxyz = pbuffer.data(idx_op_gh + 109);

    auto tr_xxzz_xxxzz = pbuffer.data(idx_op_gh + 110);

    auto tr_xxzz_xxyyy = pbuffer.data(idx_op_gh + 111);

    auto tr_xxzz_xxyyz = pbuffer.data(idx_op_gh + 112);

    auto tr_xxzz_xxyzz = pbuffer.data(idx_op_gh + 113);

    auto tr_xxzz_xxzzz = pbuffer.data(idx_op_gh + 114);

    auto tr_xxzz_xyyyy = pbuffer.data(idx_op_gh + 115);

    auto tr_xxzz_xyyyz = pbuffer.data(idx_op_gh + 116);

    auto tr_xxzz_xyyzz = pbuffer.data(idx_op_gh + 117);

    auto tr_xxzz_xyzzz = pbuffer.data(idx_op_gh + 118);

    auto tr_xxzz_xzzzz = pbuffer.data(idx_op_gh + 119);

    auto tr_xxzz_yyyyy = pbuffer.data(idx_op_gh + 120);

    auto tr_xxzz_yyyyz = pbuffer.data(idx_op_gh + 121);

    auto tr_xxzz_yyyzz = pbuffer.data(idx_op_gh + 122);

    auto tr_xxzz_yyzzz = pbuffer.data(idx_op_gh + 123);

    auto tr_xxzz_yzzzz = pbuffer.data(idx_op_gh + 124);

    auto tr_xxzz_zzzzz = pbuffer.data(idx_op_gh + 125);

    auto tr_xyyy_xxxxx = pbuffer.data(idx_op_gh + 126);

    auto tr_xyyy_xxxxy = pbuffer.data(idx_op_gh + 127);

    auto tr_xyyy_xxxxz = pbuffer.data(idx_op_gh + 128);

    auto tr_xyyy_xxxyy = pbuffer.data(idx_op_gh + 129);

    auto tr_xyyy_xxxyz = pbuffer.data(idx_op_gh + 130);

    auto tr_xyyy_xxxzz = pbuffer.data(idx_op_gh + 131);

    auto tr_xyyy_xxyyy = pbuffer.data(idx_op_gh + 132);

    auto tr_xyyy_xxyyz = pbuffer.data(idx_op_gh + 133);

    auto tr_xyyy_xxyzz = pbuffer.data(idx_op_gh + 134);

    auto tr_xyyy_xxzzz = pbuffer.data(idx_op_gh + 135);

    auto tr_xyyy_xyyyy = pbuffer.data(idx_op_gh + 136);

    auto tr_xyyy_xyyyz = pbuffer.data(idx_op_gh + 137);

    auto tr_xyyy_xyyzz = pbuffer.data(idx_op_gh + 138);

    auto tr_xyyy_xyzzz = pbuffer.data(idx_op_gh + 139);

    auto tr_xyyy_xzzzz = pbuffer.data(idx_op_gh + 140);

    auto tr_xyyy_yyyyy = pbuffer.data(idx_op_gh + 141);

    auto tr_xyyy_yyyyz = pbuffer.data(idx_op_gh + 142);

    auto tr_xyyy_yyyzz = pbuffer.data(idx_op_gh + 143);

    auto tr_xyyy_yyzzz = pbuffer.data(idx_op_gh + 144);

    auto tr_xyyy_yzzzz = pbuffer.data(idx_op_gh + 145);

    auto tr_xyyy_zzzzz = pbuffer.data(idx_op_gh + 146);

    auto tr_xyyz_xxxxx = pbuffer.data(idx_op_gh + 147);

    auto tr_xyyz_xxxxy = pbuffer.data(idx_op_gh + 148);

    auto tr_xyyz_xxxxz = pbuffer.data(idx_op_gh + 149);

    auto tr_xyyz_xxxyy = pbuffer.data(idx_op_gh + 150);

    auto tr_xyyz_xxxyz = pbuffer.data(idx_op_gh + 151);

    auto tr_xyyz_xxxzz = pbuffer.data(idx_op_gh + 152);

    auto tr_xyyz_xxyyy = pbuffer.data(idx_op_gh + 153);

    auto tr_xyyz_xxyyz = pbuffer.data(idx_op_gh + 154);

    auto tr_xyyz_xxyzz = pbuffer.data(idx_op_gh + 155);

    auto tr_xyyz_xxzzz = pbuffer.data(idx_op_gh + 156);

    auto tr_xyyz_xyyyy = pbuffer.data(idx_op_gh + 157);

    auto tr_xyyz_xyyyz = pbuffer.data(idx_op_gh + 158);

    auto tr_xyyz_xyyzz = pbuffer.data(idx_op_gh + 159);

    auto tr_xyyz_xyzzz = pbuffer.data(idx_op_gh + 160);

    auto tr_xyyz_xzzzz = pbuffer.data(idx_op_gh + 161);

    auto tr_xyyz_yyyyy = pbuffer.data(idx_op_gh + 162);

    auto tr_xyyz_yyyyz = pbuffer.data(idx_op_gh + 163);

    auto tr_xyyz_yyyzz = pbuffer.data(idx_op_gh + 164);

    auto tr_xyyz_yyzzz = pbuffer.data(idx_op_gh + 165);

    auto tr_xyyz_yzzzz = pbuffer.data(idx_op_gh + 166);

    auto tr_xyyz_zzzzz = pbuffer.data(idx_op_gh + 167);

    auto tr_xyzz_xxxxx = pbuffer.data(idx_op_gh + 168);

    auto tr_xyzz_xxxxy = pbuffer.data(idx_op_gh + 169);

    auto tr_xyzz_xxxxz = pbuffer.data(idx_op_gh + 170);

    auto tr_xyzz_xxxyy = pbuffer.data(idx_op_gh + 171);

    auto tr_xyzz_xxxyz = pbuffer.data(idx_op_gh + 172);

    auto tr_xyzz_xxxzz = pbuffer.data(idx_op_gh + 173);

    auto tr_xyzz_xxyyy = pbuffer.data(idx_op_gh + 174);

    auto tr_xyzz_xxyyz = pbuffer.data(idx_op_gh + 175);

    auto tr_xyzz_xxyzz = pbuffer.data(idx_op_gh + 176);

    auto tr_xyzz_xxzzz = pbuffer.data(idx_op_gh + 177);

    auto tr_xyzz_xyyyy = pbuffer.data(idx_op_gh + 178);

    auto tr_xyzz_xyyyz = pbuffer.data(idx_op_gh + 179);

    auto tr_xyzz_xyyzz = pbuffer.data(idx_op_gh + 180);

    auto tr_xyzz_xyzzz = pbuffer.data(idx_op_gh + 181);

    auto tr_xyzz_xzzzz = pbuffer.data(idx_op_gh + 182);

    auto tr_xyzz_yyyyy = pbuffer.data(idx_op_gh + 183);

    auto tr_xyzz_yyyyz = pbuffer.data(idx_op_gh + 184);

    auto tr_xyzz_yyyzz = pbuffer.data(idx_op_gh + 185);

    auto tr_xyzz_yyzzz = pbuffer.data(idx_op_gh + 186);

    auto tr_xyzz_yzzzz = pbuffer.data(idx_op_gh + 187);

    auto tr_xyzz_zzzzz = pbuffer.data(idx_op_gh + 188);

    auto tr_xzzz_xxxxx = pbuffer.data(idx_op_gh + 189);

    auto tr_xzzz_xxxxy = pbuffer.data(idx_op_gh + 190);

    auto tr_xzzz_xxxxz = pbuffer.data(idx_op_gh + 191);

    auto tr_xzzz_xxxyy = pbuffer.data(idx_op_gh + 192);

    auto tr_xzzz_xxxyz = pbuffer.data(idx_op_gh + 193);

    auto tr_xzzz_xxxzz = pbuffer.data(idx_op_gh + 194);

    auto tr_xzzz_xxyyy = pbuffer.data(idx_op_gh + 195);

    auto tr_xzzz_xxyyz = pbuffer.data(idx_op_gh + 196);

    auto tr_xzzz_xxyzz = pbuffer.data(idx_op_gh + 197);

    auto tr_xzzz_xxzzz = pbuffer.data(idx_op_gh + 198);

    auto tr_xzzz_xyyyy = pbuffer.data(idx_op_gh + 199);

    auto tr_xzzz_xyyyz = pbuffer.data(idx_op_gh + 200);

    auto tr_xzzz_xyyzz = pbuffer.data(idx_op_gh + 201);

    auto tr_xzzz_xyzzz = pbuffer.data(idx_op_gh + 202);

    auto tr_xzzz_xzzzz = pbuffer.data(idx_op_gh + 203);

    auto tr_xzzz_yyyyy = pbuffer.data(idx_op_gh + 204);

    auto tr_xzzz_yyyyz = pbuffer.data(idx_op_gh + 205);

    auto tr_xzzz_yyyzz = pbuffer.data(idx_op_gh + 206);

    auto tr_xzzz_yyzzz = pbuffer.data(idx_op_gh + 207);

    auto tr_xzzz_yzzzz = pbuffer.data(idx_op_gh + 208);

    auto tr_xzzz_zzzzz = pbuffer.data(idx_op_gh + 209);

    auto tr_yyyy_xxxxx = pbuffer.data(idx_op_gh + 210);

    auto tr_yyyy_xxxxy = pbuffer.data(idx_op_gh + 211);

    auto tr_yyyy_xxxxz = pbuffer.data(idx_op_gh + 212);

    auto tr_yyyy_xxxyy = pbuffer.data(idx_op_gh + 213);

    auto tr_yyyy_xxxyz = pbuffer.data(idx_op_gh + 214);

    auto tr_yyyy_xxxzz = pbuffer.data(idx_op_gh + 215);

    auto tr_yyyy_xxyyy = pbuffer.data(idx_op_gh + 216);

    auto tr_yyyy_xxyyz = pbuffer.data(idx_op_gh + 217);

    auto tr_yyyy_xxyzz = pbuffer.data(idx_op_gh + 218);

    auto tr_yyyy_xxzzz = pbuffer.data(idx_op_gh + 219);

    auto tr_yyyy_xyyyy = pbuffer.data(idx_op_gh + 220);

    auto tr_yyyy_xyyyz = pbuffer.data(idx_op_gh + 221);

    auto tr_yyyy_xyyzz = pbuffer.data(idx_op_gh + 222);

    auto tr_yyyy_xyzzz = pbuffer.data(idx_op_gh + 223);

    auto tr_yyyy_xzzzz = pbuffer.data(idx_op_gh + 224);

    auto tr_yyyy_yyyyy = pbuffer.data(idx_op_gh + 225);

    auto tr_yyyy_yyyyz = pbuffer.data(idx_op_gh + 226);

    auto tr_yyyy_yyyzz = pbuffer.data(idx_op_gh + 227);

    auto tr_yyyy_yyzzz = pbuffer.data(idx_op_gh + 228);

    auto tr_yyyy_yzzzz = pbuffer.data(idx_op_gh + 229);

    auto tr_yyyy_zzzzz = pbuffer.data(idx_op_gh + 230);

    auto tr_yyyz_xxxxx = pbuffer.data(idx_op_gh + 231);

    auto tr_yyyz_xxxxy = pbuffer.data(idx_op_gh + 232);

    auto tr_yyyz_xxxxz = pbuffer.data(idx_op_gh + 233);

    auto tr_yyyz_xxxyy = pbuffer.data(idx_op_gh + 234);

    auto tr_yyyz_xxxyz = pbuffer.data(idx_op_gh + 235);

    auto tr_yyyz_xxxzz = pbuffer.data(idx_op_gh + 236);

    auto tr_yyyz_xxyyy = pbuffer.data(idx_op_gh + 237);

    auto tr_yyyz_xxyyz = pbuffer.data(idx_op_gh + 238);

    auto tr_yyyz_xxyzz = pbuffer.data(idx_op_gh + 239);

    auto tr_yyyz_xxzzz = pbuffer.data(idx_op_gh + 240);

    auto tr_yyyz_xyyyy = pbuffer.data(idx_op_gh + 241);

    auto tr_yyyz_xyyyz = pbuffer.data(idx_op_gh + 242);

    auto tr_yyyz_xyyzz = pbuffer.data(idx_op_gh + 243);

    auto tr_yyyz_xyzzz = pbuffer.data(idx_op_gh + 244);

    auto tr_yyyz_xzzzz = pbuffer.data(idx_op_gh + 245);

    auto tr_yyyz_yyyyy = pbuffer.data(idx_op_gh + 246);

    auto tr_yyyz_yyyyz = pbuffer.data(idx_op_gh + 247);

    auto tr_yyyz_yyyzz = pbuffer.data(idx_op_gh + 248);

    auto tr_yyyz_yyzzz = pbuffer.data(idx_op_gh + 249);

    auto tr_yyyz_yzzzz = pbuffer.data(idx_op_gh + 250);

    auto tr_yyyz_zzzzz = pbuffer.data(idx_op_gh + 251);

    auto tr_yyzz_xxxxx = pbuffer.data(idx_op_gh + 252);

    auto tr_yyzz_xxxxy = pbuffer.data(idx_op_gh + 253);

    auto tr_yyzz_xxxxz = pbuffer.data(idx_op_gh + 254);

    auto tr_yyzz_xxxyy = pbuffer.data(idx_op_gh + 255);

    auto tr_yyzz_xxxyz = pbuffer.data(idx_op_gh + 256);

    auto tr_yyzz_xxxzz = pbuffer.data(idx_op_gh + 257);

    auto tr_yyzz_xxyyy = pbuffer.data(idx_op_gh + 258);

    auto tr_yyzz_xxyyz = pbuffer.data(idx_op_gh + 259);

    auto tr_yyzz_xxyzz = pbuffer.data(idx_op_gh + 260);

    auto tr_yyzz_xxzzz = pbuffer.data(idx_op_gh + 261);

    auto tr_yyzz_xyyyy = pbuffer.data(idx_op_gh + 262);

    auto tr_yyzz_xyyyz = pbuffer.data(idx_op_gh + 263);

    auto tr_yyzz_xyyzz = pbuffer.data(idx_op_gh + 264);

    auto tr_yyzz_xyzzz = pbuffer.data(idx_op_gh + 265);

    auto tr_yyzz_xzzzz = pbuffer.data(idx_op_gh + 266);

    auto tr_yyzz_yyyyy = pbuffer.data(idx_op_gh + 267);

    auto tr_yyzz_yyyyz = pbuffer.data(idx_op_gh + 268);

    auto tr_yyzz_yyyzz = pbuffer.data(idx_op_gh + 269);

    auto tr_yyzz_yyzzz = pbuffer.data(idx_op_gh + 270);

    auto tr_yyzz_yzzzz = pbuffer.data(idx_op_gh + 271);

    auto tr_yyzz_zzzzz = pbuffer.data(idx_op_gh + 272);

    auto tr_yzzz_xxxxx = pbuffer.data(idx_op_gh + 273);

    auto tr_yzzz_xxxxy = pbuffer.data(idx_op_gh + 274);

    auto tr_yzzz_xxxxz = pbuffer.data(idx_op_gh + 275);

    auto tr_yzzz_xxxyy = pbuffer.data(idx_op_gh + 276);

    auto tr_yzzz_xxxyz = pbuffer.data(idx_op_gh + 277);

    auto tr_yzzz_xxxzz = pbuffer.data(idx_op_gh + 278);

    auto tr_yzzz_xxyyy = pbuffer.data(idx_op_gh + 279);

    auto tr_yzzz_xxyyz = pbuffer.data(idx_op_gh + 280);

    auto tr_yzzz_xxyzz = pbuffer.data(idx_op_gh + 281);

    auto tr_yzzz_xxzzz = pbuffer.data(idx_op_gh + 282);

    auto tr_yzzz_xyyyy = pbuffer.data(idx_op_gh + 283);

    auto tr_yzzz_xyyyz = pbuffer.data(idx_op_gh + 284);

    auto tr_yzzz_xyyzz = pbuffer.data(idx_op_gh + 285);

    auto tr_yzzz_xyzzz = pbuffer.data(idx_op_gh + 286);

    auto tr_yzzz_xzzzz = pbuffer.data(idx_op_gh + 287);

    auto tr_yzzz_yyyyy = pbuffer.data(idx_op_gh + 288);

    auto tr_yzzz_yyyyz = pbuffer.data(idx_op_gh + 289);

    auto tr_yzzz_yyyzz = pbuffer.data(idx_op_gh + 290);

    auto tr_yzzz_yyzzz = pbuffer.data(idx_op_gh + 291);

    auto tr_yzzz_yzzzz = pbuffer.data(idx_op_gh + 292);

    auto tr_yzzz_zzzzz = pbuffer.data(idx_op_gh + 293);

    auto tr_zzzz_xxxxx = pbuffer.data(idx_op_gh + 294);

    auto tr_zzzz_xxxxy = pbuffer.data(idx_op_gh + 295);

    auto tr_zzzz_xxxxz = pbuffer.data(idx_op_gh + 296);

    auto tr_zzzz_xxxyy = pbuffer.data(idx_op_gh + 297);

    auto tr_zzzz_xxxyz = pbuffer.data(idx_op_gh + 298);

    auto tr_zzzz_xxxzz = pbuffer.data(idx_op_gh + 299);

    auto tr_zzzz_xxyyy = pbuffer.data(idx_op_gh + 300);

    auto tr_zzzz_xxyyz = pbuffer.data(idx_op_gh + 301);

    auto tr_zzzz_xxyzz = pbuffer.data(idx_op_gh + 302);

    auto tr_zzzz_xxzzz = pbuffer.data(idx_op_gh + 303);

    auto tr_zzzz_xyyyy = pbuffer.data(idx_op_gh + 304);

    auto tr_zzzz_xyyyz = pbuffer.data(idx_op_gh + 305);

    auto tr_zzzz_xyyzz = pbuffer.data(idx_op_gh + 306);

    auto tr_zzzz_xyzzz = pbuffer.data(idx_op_gh + 307);

    auto tr_zzzz_xzzzz = pbuffer.data(idx_op_gh + 308);

    auto tr_zzzz_yyyyy = pbuffer.data(idx_op_gh + 309);

    auto tr_zzzz_yyyyz = pbuffer.data(idx_op_gh + 310);

    auto tr_zzzz_yyyzz = pbuffer.data(idx_op_gh + 311);

    auto tr_zzzz_yyzzz = pbuffer.data(idx_op_gh + 312);

    auto tr_zzzz_yzzzz = pbuffer.data(idx_op_gh + 313);

    auto tr_zzzz_zzzzz = pbuffer.data(idx_op_gh + 314);

    // Set up components of auxiliary buffer : HG

    auto tr_xxxxx_xxxx = pbuffer.data(idx_op_hg);

    auto tr_xxxxx_xxxy = pbuffer.data(idx_op_hg + 1);

    auto tr_xxxxx_xxxz = pbuffer.data(idx_op_hg + 2);

    auto tr_xxxxx_xxyy = pbuffer.data(idx_op_hg + 3);

    auto tr_xxxxx_xxyz = pbuffer.data(idx_op_hg + 4);

    auto tr_xxxxx_xxzz = pbuffer.data(idx_op_hg + 5);

    auto tr_xxxxx_xyyy = pbuffer.data(idx_op_hg + 6);

    auto tr_xxxxx_xyyz = pbuffer.data(idx_op_hg + 7);

    auto tr_xxxxx_xyzz = pbuffer.data(idx_op_hg + 8);

    auto tr_xxxxx_xzzz = pbuffer.data(idx_op_hg + 9);

    auto tr_xxxxx_yyyy = pbuffer.data(idx_op_hg + 10);

    auto tr_xxxxx_yyyz = pbuffer.data(idx_op_hg + 11);

    auto tr_xxxxx_yyzz = pbuffer.data(idx_op_hg + 12);

    auto tr_xxxxx_yzzz = pbuffer.data(idx_op_hg + 13);

    auto tr_xxxxx_zzzz = pbuffer.data(idx_op_hg + 14);

    auto tr_xxxxy_xxxx = pbuffer.data(idx_op_hg + 15);

    auto tr_xxxxy_xxxy = pbuffer.data(idx_op_hg + 16);

    auto tr_xxxxy_xxxz = pbuffer.data(idx_op_hg + 17);

    auto tr_xxxxy_xxyy = pbuffer.data(idx_op_hg + 18);

    auto tr_xxxxy_xxyz = pbuffer.data(idx_op_hg + 19);

    auto tr_xxxxy_xxzz = pbuffer.data(idx_op_hg + 20);

    auto tr_xxxxy_xyyy = pbuffer.data(idx_op_hg + 21);

    auto tr_xxxxy_xyyz = pbuffer.data(idx_op_hg + 22);

    auto tr_xxxxy_xyzz = pbuffer.data(idx_op_hg + 23);

    auto tr_xxxxy_xzzz = pbuffer.data(idx_op_hg + 24);

    auto tr_xxxxy_yyyy = pbuffer.data(idx_op_hg + 25);

    auto tr_xxxxy_yyyz = pbuffer.data(idx_op_hg + 26);

    auto tr_xxxxy_yyzz = pbuffer.data(idx_op_hg + 27);

    auto tr_xxxxy_yzzz = pbuffer.data(idx_op_hg + 28);

    auto tr_xxxxy_zzzz = pbuffer.data(idx_op_hg + 29);

    auto tr_xxxxz_xxxx = pbuffer.data(idx_op_hg + 30);

    auto tr_xxxxz_xxxy = pbuffer.data(idx_op_hg + 31);

    auto tr_xxxxz_xxxz = pbuffer.data(idx_op_hg + 32);

    auto tr_xxxxz_xxyy = pbuffer.data(idx_op_hg + 33);

    auto tr_xxxxz_xxyz = pbuffer.data(idx_op_hg + 34);

    auto tr_xxxxz_xxzz = pbuffer.data(idx_op_hg + 35);

    auto tr_xxxxz_xyyy = pbuffer.data(idx_op_hg + 36);

    auto tr_xxxxz_xyyz = pbuffer.data(idx_op_hg + 37);

    auto tr_xxxxz_xyzz = pbuffer.data(idx_op_hg + 38);

    auto tr_xxxxz_xzzz = pbuffer.data(idx_op_hg + 39);

    auto tr_xxxxz_yyyy = pbuffer.data(idx_op_hg + 40);

    auto tr_xxxxz_yyyz = pbuffer.data(idx_op_hg + 41);

    auto tr_xxxxz_yyzz = pbuffer.data(idx_op_hg + 42);

    auto tr_xxxxz_yzzz = pbuffer.data(idx_op_hg + 43);

    auto tr_xxxxz_zzzz = pbuffer.data(idx_op_hg + 44);

    auto tr_xxxyy_xxxx = pbuffer.data(idx_op_hg + 45);

    auto tr_xxxyy_xxxy = pbuffer.data(idx_op_hg + 46);

    auto tr_xxxyy_xxxz = pbuffer.data(idx_op_hg + 47);

    auto tr_xxxyy_xxyy = pbuffer.data(idx_op_hg + 48);

    auto tr_xxxyy_xxyz = pbuffer.data(idx_op_hg + 49);

    auto tr_xxxyy_xxzz = pbuffer.data(idx_op_hg + 50);

    auto tr_xxxyy_xyyy = pbuffer.data(idx_op_hg + 51);

    auto tr_xxxyy_xyyz = pbuffer.data(idx_op_hg + 52);

    auto tr_xxxyy_xyzz = pbuffer.data(idx_op_hg + 53);

    auto tr_xxxyy_xzzz = pbuffer.data(idx_op_hg + 54);

    auto tr_xxxyy_yyyy = pbuffer.data(idx_op_hg + 55);

    auto tr_xxxyy_yyyz = pbuffer.data(idx_op_hg + 56);

    auto tr_xxxyy_yyzz = pbuffer.data(idx_op_hg + 57);

    auto tr_xxxyy_yzzz = pbuffer.data(idx_op_hg + 58);

    auto tr_xxxyy_zzzz = pbuffer.data(idx_op_hg + 59);

    auto tr_xxxyz_xxxx = pbuffer.data(idx_op_hg + 60);

    auto tr_xxxyz_xxxy = pbuffer.data(idx_op_hg + 61);

    auto tr_xxxyz_xxxz = pbuffer.data(idx_op_hg + 62);

    auto tr_xxxyz_xxyy = pbuffer.data(idx_op_hg + 63);

    auto tr_xxxyz_xxyz = pbuffer.data(idx_op_hg + 64);

    auto tr_xxxyz_xxzz = pbuffer.data(idx_op_hg + 65);

    auto tr_xxxyz_xyyy = pbuffer.data(idx_op_hg + 66);

    auto tr_xxxyz_xyyz = pbuffer.data(idx_op_hg + 67);

    auto tr_xxxyz_xyzz = pbuffer.data(idx_op_hg + 68);

    auto tr_xxxyz_xzzz = pbuffer.data(idx_op_hg + 69);

    auto tr_xxxyz_yyyy = pbuffer.data(idx_op_hg + 70);

    auto tr_xxxyz_yyyz = pbuffer.data(idx_op_hg + 71);

    auto tr_xxxyz_yyzz = pbuffer.data(idx_op_hg + 72);

    auto tr_xxxyz_yzzz = pbuffer.data(idx_op_hg + 73);

    auto tr_xxxyz_zzzz = pbuffer.data(idx_op_hg + 74);

    auto tr_xxxzz_xxxx = pbuffer.data(idx_op_hg + 75);

    auto tr_xxxzz_xxxy = pbuffer.data(idx_op_hg + 76);

    auto tr_xxxzz_xxxz = pbuffer.data(idx_op_hg + 77);

    auto tr_xxxzz_xxyy = pbuffer.data(idx_op_hg + 78);

    auto tr_xxxzz_xxyz = pbuffer.data(idx_op_hg + 79);

    auto tr_xxxzz_xxzz = pbuffer.data(idx_op_hg + 80);

    auto tr_xxxzz_xyyy = pbuffer.data(idx_op_hg + 81);

    auto tr_xxxzz_xyyz = pbuffer.data(idx_op_hg + 82);

    auto tr_xxxzz_xyzz = pbuffer.data(idx_op_hg + 83);

    auto tr_xxxzz_xzzz = pbuffer.data(idx_op_hg + 84);

    auto tr_xxxzz_yyyy = pbuffer.data(idx_op_hg + 85);

    auto tr_xxxzz_yyyz = pbuffer.data(idx_op_hg + 86);

    auto tr_xxxzz_yyzz = pbuffer.data(idx_op_hg + 87);

    auto tr_xxxzz_yzzz = pbuffer.data(idx_op_hg + 88);

    auto tr_xxxzz_zzzz = pbuffer.data(idx_op_hg + 89);

    auto tr_xxyyy_xxxx = pbuffer.data(idx_op_hg + 90);

    auto tr_xxyyy_xxxy = pbuffer.data(idx_op_hg + 91);

    auto tr_xxyyy_xxxz = pbuffer.data(idx_op_hg + 92);

    auto tr_xxyyy_xxyy = pbuffer.data(idx_op_hg + 93);

    auto tr_xxyyy_xxyz = pbuffer.data(idx_op_hg + 94);

    auto tr_xxyyy_xxzz = pbuffer.data(idx_op_hg + 95);

    auto tr_xxyyy_xyyy = pbuffer.data(idx_op_hg + 96);

    auto tr_xxyyy_xyyz = pbuffer.data(idx_op_hg + 97);

    auto tr_xxyyy_xyzz = pbuffer.data(idx_op_hg + 98);

    auto tr_xxyyy_xzzz = pbuffer.data(idx_op_hg + 99);

    auto tr_xxyyy_yyyy = pbuffer.data(idx_op_hg + 100);

    auto tr_xxyyy_yyyz = pbuffer.data(idx_op_hg + 101);

    auto tr_xxyyy_yyzz = pbuffer.data(idx_op_hg + 102);

    auto tr_xxyyy_yzzz = pbuffer.data(idx_op_hg + 103);

    auto tr_xxyyy_zzzz = pbuffer.data(idx_op_hg + 104);

    auto tr_xxyyz_xxxx = pbuffer.data(idx_op_hg + 105);

    auto tr_xxyyz_xxxy = pbuffer.data(idx_op_hg + 106);

    auto tr_xxyyz_xxxz = pbuffer.data(idx_op_hg + 107);

    auto tr_xxyyz_xxyy = pbuffer.data(idx_op_hg + 108);

    auto tr_xxyyz_xxyz = pbuffer.data(idx_op_hg + 109);

    auto tr_xxyyz_xxzz = pbuffer.data(idx_op_hg + 110);

    auto tr_xxyyz_xyyy = pbuffer.data(idx_op_hg + 111);

    auto tr_xxyyz_xyyz = pbuffer.data(idx_op_hg + 112);

    auto tr_xxyyz_xyzz = pbuffer.data(idx_op_hg + 113);

    auto tr_xxyyz_xzzz = pbuffer.data(idx_op_hg + 114);

    auto tr_xxyyz_yyyy = pbuffer.data(idx_op_hg + 115);

    auto tr_xxyyz_yyyz = pbuffer.data(idx_op_hg + 116);

    auto tr_xxyyz_yyzz = pbuffer.data(idx_op_hg + 117);

    auto tr_xxyyz_yzzz = pbuffer.data(idx_op_hg + 118);

    auto tr_xxyyz_zzzz = pbuffer.data(idx_op_hg + 119);

    auto tr_xxyzz_xxxx = pbuffer.data(idx_op_hg + 120);

    auto tr_xxyzz_xxxy = pbuffer.data(idx_op_hg + 121);

    auto tr_xxyzz_xxxz = pbuffer.data(idx_op_hg + 122);

    auto tr_xxyzz_xxyy = pbuffer.data(idx_op_hg + 123);

    auto tr_xxyzz_xxyz = pbuffer.data(idx_op_hg + 124);

    auto tr_xxyzz_xxzz = pbuffer.data(idx_op_hg + 125);

    auto tr_xxyzz_xyyy = pbuffer.data(idx_op_hg + 126);

    auto tr_xxyzz_xyyz = pbuffer.data(idx_op_hg + 127);

    auto tr_xxyzz_xyzz = pbuffer.data(idx_op_hg + 128);

    auto tr_xxyzz_xzzz = pbuffer.data(idx_op_hg + 129);

    auto tr_xxyzz_yyyy = pbuffer.data(idx_op_hg + 130);

    auto tr_xxyzz_yyyz = pbuffer.data(idx_op_hg + 131);

    auto tr_xxyzz_yyzz = pbuffer.data(idx_op_hg + 132);

    auto tr_xxyzz_yzzz = pbuffer.data(idx_op_hg + 133);

    auto tr_xxyzz_zzzz = pbuffer.data(idx_op_hg + 134);

    auto tr_xxzzz_xxxx = pbuffer.data(idx_op_hg + 135);

    auto tr_xxzzz_xxxy = pbuffer.data(idx_op_hg + 136);

    auto tr_xxzzz_xxxz = pbuffer.data(idx_op_hg + 137);

    auto tr_xxzzz_xxyy = pbuffer.data(idx_op_hg + 138);

    auto tr_xxzzz_xxyz = pbuffer.data(idx_op_hg + 139);

    auto tr_xxzzz_xxzz = pbuffer.data(idx_op_hg + 140);

    auto tr_xxzzz_xyyy = pbuffer.data(idx_op_hg + 141);

    auto tr_xxzzz_xyyz = pbuffer.data(idx_op_hg + 142);

    auto tr_xxzzz_xyzz = pbuffer.data(idx_op_hg + 143);

    auto tr_xxzzz_xzzz = pbuffer.data(idx_op_hg + 144);

    auto tr_xxzzz_yyyy = pbuffer.data(idx_op_hg + 145);

    auto tr_xxzzz_yyyz = pbuffer.data(idx_op_hg + 146);

    auto tr_xxzzz_yyzz = pbuffer.data(idx_op_hg + 147);

    auto tr_xxzzz_yzzz = pbuffer.data(idx_op_hg + 148);

    auto tr_xxzzz_zzzz = pbuffer.data(idx_op_hg + 149);

    auto tr_xyyyy_xxxx = pbuffer.data(idx_op_hg + 150);

    auto tr_xyyyy_xxxy = pbuffer.data(idx_op_hg + 151);

    auto tr_xyyyy_xxxz = pbuffer.data(idx_op_hg + 152);

    auto tr_xyyyy_xxyy = pbuffer.data(idx_op_hg + 153);

    auto tr_xyyyy_xxyz = pbuffer.data(idx_op_hg + 154);

    auto tr_xyyyy_xxzz = pbuffer.data(idx_op_hg + 155);

    auto tr_xyyyy_xyyy = pbuffer.data(idx_op_hg + 156);

    auto tr_xyyyy_xyyz = pbuffer.data(idx_op_hg + 157);

    auto tr_xyyyy_xyzz = pbuffer.data(idx_op_hg + 158);

    auto tr_xyyyy_xzzz = pbuffer.data(idx_op_hg + 159);

    auto tr_xyyyy_yyyy = pbuffer.data(idx_op_hg + 160);

    auto tr_xyyyy_yyyz = pbuffer.data(idx_op_hg + 161);

    auto tr_xyyyy_yyzz = pbuffer.data(idx_op_hg + 162);

    auto tr_xyyyy_yzzz = pbuffer.data(idx_op_hg + 163);

    auto tr_xyyyy_zzzz = pbuffer.data(idx_op_hg + 164);

    auto tr_xyyyz_xxxx = pbuffer.data(idx_op_hg + 165);

    auto tr_xyyyz_xxxy = pbuffer.data(idx_op_hg + 166);

    auto tr_xyyyz_xxxz = pbuffer.data(idx_op_hg + 167);

    auto tr_xyyyz_xxyy = pbuffer.data(idx_op_hg + 168);

    auto tr_xyyyz_xxyz = pbuffer.data(idx_op_hg + 169);

    auto tr_xyyyz_xxzz = pbuffer.data(idx_op_hg + 170);

    auto tr_xyyyz_xyyy = pbuffer.data(idx_op_hg + 171);

    auto tr_xyyyz_xyyz = pbuffer.data(idx_op_hg + 172);

    auto tr_xyyyz_xyzz = pbuffer.data(idx_op_hg + 173);

    auto tr_xyyyz_xzzz = pbuffer.data(idx_op_hg + 174);

    auto tr_xyyyz_yyyy = pbuffer.data(idx_op_hg + 175);

    auto tr_xyyyz_yyyz = pbuffer.data(idx_op_hg + 176);

    auto tr_xyyyz_yyzz = pbuffer.data(idx_op_hg + 177);

    auto tr_xyyyz_yzzz = pbuffer.data(idx_op_hg + 178);

    auto tr_xyyyz_zzzz = pbuffer.data(idx_op_hg + 179);

    auto tr_xyyzz_xxxx = pbuffer.data(idx_op_hg + 180);

    auto tr_xyyzz_xxxy = pbuffer.data(idx_op_hg + 181);

    auto tr_xyyzz_xxxz = pbuffer.data(idx_op_hg + 182);

    auto tr_xyyzz_xxyy = pbuffer.data(idx_op_hg + 183);

    auto tr_xyyzz_xxyz = pbuffer.data(idx_op_hg + 184);

    auto tr_xyyzz_xxzz = pbuffer.data(idx_op_hg + 185);

    auto tr_xyyzz_xyyy = pbuffer.data(idx_op_hg + 186);

    auto tr_xyyzz_xyyz = pbuffer.data(idx_op_hg + 187);

    auto tr_xyyzz_xyzz = pbuffer.data(idx_op_hg + 188);

    auto tr_xyyzz_xzzz = pbuffer.data(idx_op_hg + 189);

    auto tr_xyyzz_yyyy = pbuffer.data(idx_op_hg + 190);

    auto tr_xyyzz_yyyz = pbuffer.data(idx_op_hg + 191);

    auto tr_xyyzz_yyzz = pbuffer.data(idx_op_hg + 192);

    auto tr_xyyzz_yzzz = pbuffer.data(idx_op_hg + 193);

    auto tr_xyyzz_zzzz = pbuffer.data(idx_op_hg + 194);

    auto tr_xyzzz_xxxx = pbuffer.data(idx_op_hg + 195);

    auto tr_xyzzz_xxxy = pbuffer.data(idx_op_hg + 196);

    auto tr_xyzzz_xxxz = pbuffer.data(idx_op_hg + 197);

    auto tr_xyzzz_xxyy = pbuffer.data(idx_op_hg + 198);

    auto tr_xyzzz_xxyz = pbuffer.data(idx_op_hg + 199);

    auto tr_xyzzz_xxzz = pbuffer.data(idx_op_hg + 200);

    auto tr_xyzzz_xyyy = pbuffer.data(idx_op_hg + 201);

    auto tr_xyzzz_xyyz = pbuffer.data(idx_op_hg + 202);

    auto tr_xyzzz_xyzz = pbuffer.data(idx_op_hg + 203);

    auto tr_xyzzz_xzzz = pbuffer.data(idx_op_hg + 204);

    auto tr_xyzzz_yyyy = pbuffer.data(idx_op_hg + 205);

    auto tr_xyzzz_yyyz = pbuffer.data(idx_op_hg + 206);

    auto tr_xyzzz_yyzz = pbuffer.data(idx_op_hg + 207);

    auto tr_xyzzz_yzzz = pbuffer.data(idx_op_hg + 208);

    auto tr_xyzzz_zzzz = pbuffer.data(idx_op_hg + 209);

    auto tr_xzzzz_xxxx = pbuffer.data(idx_op_hg + 210);

    auto tr_xzzzz_xxxy = pbuffer.data(idx_op_hg + 211);

    auto tr_xzzzz_xxxz = pbuffer.data(idx_op_hg + 212);

    auto tr_xzzzz_xxyy = pbuffer.data(idx_op_hg + 213);

    auto tr_xzzzz_xxyz = pbuffer.data(idx_op_hg + 214);

    auto tr_xzzzz_xxzz = pbuffer.data(idx_op_hg + 215);

    auto tr_xzzzz_xyyy = pbuffer.data(idx_op_hg + 216);

    auto tr_xzzzz_xyyz = pbuffer.data(idx_op_hg + 217);

    auto tr_xzzzz_xyzz = pbuffer.data(idx_op_hg + 218);

    auto tr_xzzzz_xzzz = pbuffer.data(idx_op_hg + 219);

    auto tr_xzzzz_yyyy = pbuffer.data(idx_op_hg + 220);

    auto tr_xzzzz_yyyz = pbuffer.data(idx_op_hg + 221);

    auto tr_xzzzz_yyzz = pbuffer.data(idx_op_hg + 222);

    auto tr_xzzzz_yzzz = pbuffer.data(idx_op_hg + 223);

    auto tr_xzzzz_zzzz = pbuffer.data(idx_op_hg + 224);

    auto tr_yyyyy_xxxx = pbuffer.data(idx_op_hg + 225);

    auto tr_yyyyy_xxxy = pbuffer.data(idx_op_hg + 226);

    auto tr_yyyyy_xxxz = pbuffer.data(idx_op_hg + 227);

    auto tr_yyyyy_xxyy = pbuffer.data(idx_op_hg + 228);

    auto tr_yyyyy_xxyz = pbuffer.data(idx_op_hg + 229);

    auto tr_yyyyy_xxzz = pbuffer.data(idx_op_hg + 230);

    auto tr_yyyyy_xyyy = pbuffer.data(idx_op_hg + 231);

    auto tr_yyyyy_xyyz = pbuffer.data(idx_op_hg + 232);

    auto tr_yyyyy_xyzz = pbuffer.data(idx_op_hg + 233);

    auto tr_yyyyy_xzzz = pbuffer.data(idx_op_hg + 234);

    auto tr_yyyyy_yyyy = pbuffer.data(idx_op_hg + 235);

    auto tr_yyyyy_yyyz = pbuffer.data(idx_op_hg + 236);

    auto tr_yyyyy_yyzz = pbuffer.data(idx_op_hg + 237);

    auto tr_yyyyy_yzzz = pbuffer.data(idx_op_hg + 238);

    auto tr_yyyyy_zzzz = pbuffer.data(idx_op_hg + 239);

    auto tr_yyyyz_xxxx = pbuffer.data(idx_op_hg + 240);

    auto tr_yyyyz_xxxy = pbuffer.data(idx_op_hg + 241);

    auto tr_yyyyz_xxxz = pbuffer.data(idx_op_hg + 242);

    auto tr_yyyyz_xxyy = pbuffer.data(idx_op_hg + 243);

    auto tr_yyyyz_xxyz = pbuffer.data(idx_op_hg + 244);

    auto tr_yyyyz_xxzz = pbuffer.data(idx_op_hg + 245);

    auto tr_yyyyz_xyyy = pbuffer.data(idx_op_hg + 246);

    auto tr_yyyyz_xyyz = pbuffer.data(idx_op_hg + 247);

    auto tr_yyyyz_xyzz = pbuffer.data(idx_op_hg + 248);

    auto tr_yyyyz_xzzz = pbuffer.data(idx_op_hg + 249);

    auto tr_yyyyz_yyyy = pbuffer.data(idx_op_hg + 250);

    auto tr_yyyyz_yyyz = pbuffer.data(idx_op_hg + 251);

    auto tr_yyyyz_yyzz = pbuffer.data(idx_op_hg + 252);

    auto tr_yyyyz_yzzz = pbuffer.data(idx_op_hg + 253);

    auto tr_yyyyz_zzzz = pbuffer.data(idx_op_hg + 254);

    auto tr_yyyzz_xxxx = pbuffer.data(idx_op_hg + 255);

    auto tr_yyyzz_xxxy = pbuffer.data(idx_op_hg + 256);

    auto tr_yyyzz_xxxz = pbuffer.data(idx_op_hg + 257);

    auto tr_yyyzz_xxyy = pbuffer.data(idx_op_hg + 258);

    auto tr_yyyzz_xxyz = pbuffer.data(idx_op_hg + 259);

    auto tr_yyyzz_xxzz = pbuffer.data(idx_op_hg + 260);

    auto tr_yyyzz_xyyy = pbuffer.data(idx_op_hg + 261);

    auto tr_yyyzz_xyyz = pbuffer.data(idx_op_hg + 262);

    auto tr_yyyzz_xyzz = pbuffer.data(idx_op_hg + 263);

    auto tr_yyyzz_xzzz = pbuffer.data(idx_op_hg + 264);

    auto tr_yyyzz_yyyy = pbuffer.data(idx_op_hg + 265);

    auto tr_yyyzz_yyyz = pbuffer.data(idx_op_hg + 266);

    auto tr_yyyzz_yyzz = pbuffer.data(idx_op_hg + 267);

    auto tr_yyyzz_yzzz = pbuffer.data(idx_op_hg + 268);

    auto tr_yyyzz_zzzz = pbuffer.data(idx_op_hg + 269);

    auto tr_yyzzz_xxxx = pbuffer.data(idx_op_hg + 270);

    auto tr_yyzzz_xxxy = pbuffer.data(idx_op_hg + 271);

    auto tr_yyzzz_xxxz = pbuffer.data(idx_op_hg + 272);

    auto tr_yyzzz_xxyy = pbuffer.data(idx_op_hg + 273);

    auto tr_yyzzz_xxyz = pbuffer.data(idx_op_hg + 274);

    auto tr_yyzzz_xxzz = pbuffer.data(idx_op_hg + 275);

    auto tr_yyzzz_xyyy = pbuffer.data(idx_op_hg + 276);

    auto tr_yyzzz_xyyz = pbuffer.data(idx_op_hg + 277);

    auto tr_yyzzz_xyzz = pbuffer.data(idx_op_hg + 278);

    auto tr_yyzzz_xzzz = pbuffer.data(idx_op_hg + 279);

    auto tr_yyzzz_yyyy = pbuffer.data(idx_op_hg + 280);

    auto tr_yyzzz_yyyz = pbuffer.data(idx_op_hg + 281);

    auto tr_yyzzz_yyzz = pbuffer.data(idx_op_hg + 282);

    auto tr_yyzzz_yzzz = pbuffer.data(idx_op_hg + 283);

    auto tr_yyzzz_zzzz = pbuffer.data(idx_op_hg + 284);

    auto tr_yzzzz_xxxx = pbuffer.data(idx_op_hg + 285);

    auto tr_yzzzz_xxxy = pbuffer.data(idx_op_hg + 286);

    auto tr_yzzzz_xxxz = pbuffer.data(idx_op_hg + 287);

    auto tr_yzzzz_xxyy = pbuffer.data(idx_op_hg + 288);

    auto tr_yzzzz_xxyz = pbuffer.data(idx_op_hg + 289);

    auto tr_yzzzz_xxzz = pbuffer.data(idx_op_hg + 290);

    auto tr_yzzzz_xyyy = pbuffer.data(idx_op_hg + 291);

    auto tr_yzzzz_xyyz = pbuffer.data(idx_op_hg + 292);

    auto tr_yzzzz_xyzz = pbuffer.data(idx_op_hg + 293);

    auto tr_yzzzz_xzzz = pbuffer.data(idx_op_hg + 294);

    auto tr_yzzzz_yyyy = pbuffer.data(idx_op_hg + 295);

    auto tr_yzzzz_yyyz = pbuffer.data(idx_op_hg + 296);

    auto tr_yzzzz_yyzz = pbuffer.data(idx_op_hg + 297);

    auto tr_yzzzz_yzzz = pbuffer.data(idx_op_hg + 298);

    auto tr_yzzzz_zzzz = pbuffer.data(idx_op_hg + 299);

    auto tr_zzzzz_xxxx = pbuffer.data(idx_op_hg + 300);

    auto tr_zzzzz_xxxy = pbuffer.data(idx_op_hg + 301);

    auto tr_zzzzz_xxxz = pbuffer.data(idx_op_hg + 302);

    auto tr_zzzzz_xxyy = pbuffer.data(idx_op_hg + 303);

    auto tr_zzzzz_xxyz = pbuffer.data(idx_op_hg + 304);

    auto tr_zzzzz_xxzz = pbuffer.data(idx_op_hg + 305);

    auto tr_zzzzz_xyyy = pbuffer.data(idx_op_hg + 306);

    auto tr_zzzzz_xyyz = pbuffer.data(idx_op_hg + 307);

    auto tr_zzzzz_xyzz = pbuffer.data(idx_op_hg + 308);

    auto tr_zzzzz_xzzz = pbuffer.data(idx_op_hg + 309);

    auto tr_zzzzz_yyyy = pbuffer.data(idx_op_hg + 310);

    auto tr_zzzzz_yyyz = pbuffer.data(idx_op_hg + 311);

    auto tr_zzzzz_yyzz = pbuffer.data(idx_op_hg + 312);

    auto tr_zzzzz_yzzz = pbuffer.data(idx_op_hg + 313);

    auto tr_zzzzz_zzzz = pbuffer.data(idx_op_hg + 314);

    // Set up 0-15 components of targeted buffer : FG

    auto tr_x_0_x_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg);

    auto tr_x_0_x_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 1);

    auto tr_x_0_x_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 2);

    auto tr_x_0_x_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 3);

    auto tr_x_0_x_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 4);

    auto tr_x_0_x_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 5);

    auto tr_x_0_x_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 6);

    auto tr_x_0_x_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 7);

    auto tr_x_0_x_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 8);

    auto tr_x_0_x_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 9);

    auto tr_x_0_x_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 10);

    auto tr_x_0_x_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 11);

    auto tr_x_0_x_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 12);

    auto tr_x_0_x_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 13);

    auto tr_x_0_x_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 14);

    #pragma omp simd aligned(tr_x_0_x_xxx_xxxx, tr_x_0_x_xxx_xxxy, tr_x_0_x_xxx_xxxz, tr_x_0_x_xxx_xxyy, tr_x_0_x_xxx_xxyz, tr_x_0_x_xxx_xxzz, tr_x_0_x_xxx_xyyy, tr_x_0_x_xxx_xyyz, tr_x_0_x_xxx_xyzz, tr_x_0_x_xxx_xzzz, tr_x_0_x_xxx_yyyy, tr_x_0_x_xxx_yyyz, tr_x_0_x_xxx_yyzz, tr_x_0_x_xxx_yzzz, tr_x_0_x_xxx_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxx_xxx, tr_xxxx_xxxxx, tr_xxxx_xxxxy, tr_xxxx_xxxxz, tr_xxxx_xxxyy, tr_xxxx_xxxyz, tr_xxxx_xxxzz, tr_xxxx_xxy, tr_xxxx_xxyyy, tr_xxxx_xxyyz, tr_xxxx_xxyzz, tr_xxxx_xxz, tr_xxxx_xxzzz, tr_xxxx_xyy, tr_xxxx_xyyyy, tr_xxxx_xyyyz, tr_xxxx_xyyzz, tr_xxxx_xyz, tr_xxxx_xyzzz, tr_xxxx_xzz, tr_xxxx_xzzzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxx_xxxx, tr_xxxxx_xxxy, tr_xxxxx_xxxz, tr_xxxxx_xxyy, tr_xxxxx_xxyz, tr_xxxxx_xxzz, tr_xxxxx_xyyy, tr_xxxxx_xyyz, tr_xxxxx_xyzz, tr_xxxxx_xzzz, tr_xxxxx_yyyy, tr_xxxxx_yyyz, tr_xxxxx_yyzz, tr_xxxxx_yzzz, tr_xxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxx_xxxx[i] = 6.0 * tr_x_xxxx[i] + 12.0 * tr_xx_xxx[i] - 6.0 * tr_xx_xxxxx[i] * tke_0 - 14.0 * tr_xxx_xxxx[i] * tbe_0 - 8.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxx_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xxxy[i] = 6.0 * tr_x_xxxy[i] + 9.0 * tr_xx_xxy[i] - 6.0 * tr_xx_xxxxy[i] * tke_0 - 14.0 * tr_xxx_xxxy[i] * tbe_0 - 6.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xxxz[i] = 6.0 * tr_x_xxxz[i] + 9.0 * tr_xx_xxz[i] - 6.0 * tr_xx_xxxxz[i] * tke_0 - 14.0 * tr_xxx_xxxz[i] * tbe_0 - 6.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xxyy[i] = 6.0 * tr_x_xxyy[i] + 6.0 * tr_xx_xyy[i] - 6.0 * tr_xx_xxxyy[i] * tke_0 - 14.0 * tr_xxx_xxyy[i] * tbe_0 - 4.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xxyz[i] = 6.0 * tr_x_xxyz[i] + 6.0 * tr_xx_xyz[i] - 6.0 * tr_xx_xxxyz[i] * tke_0 - 14.0 * tr_xxx_xxyz[i] * tbe_0 - 4.0 * tr_xxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xxzz[i] = 6.0 * tr_x_xxzz[i] + 6.0 * tr_xx_xzz[i] - 6.0 * tr_xx_xxxzz[i] * tke_0 - 14.0 * tr_xxx_xxzz[i] * tbe_0 - 4.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xyyy[i] = 6.0 * tr_x_xyyy[i] + 3.0 * tr_xx_yyy[i] - 6.0 * tr_xx_xxyyy[i] * tke_0 - 14.0 * tr_xxx_xyyy[i] * tbe_0 - 2.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xyyz[i] = 6.0 * tr_x_xyyz[i] + 3.0 * tr_xx_yyz[i] - 6.0 * tr_xx_xxyyz[i] * tke_0 - 14.0 * tr_xxx_xyyz[i] * tbe_0 - 2.0 * tr_xxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xyzz[i] = 6.0 * tr_x_xyzz[i] + 3.0 * tr_xx_yzz[i] - 6.0 * tr_xx_xxyzz[i] * tke_0 - 14.0 * tr_xxx_xyzz[i] * tbe_0 - 2.0 * tr_xxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xzzz[i] = 6.0 * tr_x_xzzz[i] + 3.0 * tr_xx_zzz[i] - 6.0 * tr_xx_xxzzz[i] * tke_0 - 14.0 * tr_xxx_xzzz[i] * tbe_0 - 2.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_yyyy[i] = 6.0 * tr_x_yyyy[i] - 6.0 * tr_xx_xyyyy[i] * tke_0 - 14.0 * tr_xxx_yyyy[i] * tbe_0 + 4.0 * tr_xxxx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_yyyz[i] = 6.0 * tr_x_yyyz[i] - 6.0 * tr_xx_xyyyz[i] * tke_0 - 14.0 * tr_xxx_yyyz[i] * tbe_0 + 4.0 * tr_xxxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_yyzz[i] = 6.0 * tr_x_yyzz[i] - 6.0 * tr_xx_xyyzz[i] * tke_0 - 14.0 * tr_xxx_yyzz[i] * tbe_0 + 4.0 * tr_xxxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_yzzz[i] = 6.0 * tr_x_yzzz[i] - 6.0 * tr_xx_xyzzz[i] * tke_0 - 14.0 * tr_xxx_yzzz[i] * tbe_0 + 4.0 * tr_xxxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_zzzz[i] = 6.0 * tr_x_zzzz[i] - 6.0 * tr_xx_xzzzz[i] * tke_0 - 14.0 * tr_xxx_zzzz[i] * tbe_0 + 4.0 * tr_xxxx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 15-30 components of targeted buffer : FG

    auto tr_x_0_x_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 15);

    auto tr_x_0_x_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 16);

    auto tr_x_0_x_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 17);

    auto tr_x_0_x_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 18);

    auto tr_x_0_x_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 19);

    auto tr_x_0_x_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 20);

    auto tr_x_0_x_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 21);

    auto tr_x_0_x_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 22);

    auto tr_x_0_x_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 23);

    auto tr_x_0_x_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 24);

    auto tr_x_0_x_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 25);

    auto tr_x_0_x_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 26);

    auto tr_x_0_x_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 27);

    auto tr_x_0_x_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 28);

    auto tr_x_0_x_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 29);

    #pragma omp simd aligned(tr_x_0_x_xxy_xxxx, tr_x_0_x_xxy_xxxy, tr_x_0_x_xxy_xxxz, tr_x_0_x_xxy_xxyy, tr_x_0_x_xxy_xxyz, tr_x_0_x_xxy_xxzz, tr_x_0_x_xxy_xyyy, tr_x_0_x_xxy_xyyz, tr_x_0_x_xxy_xyzz, tr_x_0_x_xxy_xzzz, tr_x_0_x_xxy_yyyy, tr_x_0_x_xxy_yyyz, tr_x_0_x_xxy_yyzz, tr_x_0_x_xxy_yzzz, tr_x_0_x_xxy_zzzz, tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xzzz, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yzzz, tr_xxxxy_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxx, tr_xxxy_xxxxy, tr_xxxy_xxxxz, tr_xxxy_xxxyy, tr_xxxy_xxxyz, tr_xxxy_xxxzz, tr_xxxy_xxy, tr_xxxy_xxyyy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xxzzz, tr_xxxy_xyy, tr_xxxy_xyyyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_xzzzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxy_xxxx[i] = 2.0 * tr_y_xxxx[i] + 8.0 * tr_xy_xxx[i] - 4.0 * tr_xy_xxxxx[i] * tke_0 - 10.0 * tr_xxy_xxxx[i] * tbe_0 - 8.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xxxy[i] = 2.0 * tr_y_xxxy[i] + 6.0 * tr_xy_xxy[i] - 4.0 * tr_xy_xxxxy[i] * tke_0 - 10.0 * tr_xxy_xxxy[i] * tbe_0 - 6.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xxxz[i] = 2.0 * tr_y_xxxz[i] + 6.0 * tr_xy_xxz[i] - 4.0 * tr_xy_xxxxz[i] * tke_0 - 10.0 * tr_xxy_xxxz[i] * tbe_0 - 6.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xxyy[i] = 2.0 * tr_y_xxyy[i] + 4.0 * tr_xy_xyy[i] - 4.0 * tr_xy_xxxyy[i] * tke_0 - 10.0 * tr_xxy_xxyy[i] * tbe_0 - 4.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xxyz[i] = 2.0 * tr_y_xxyz[i] + 4.0 * tr_xy_xyz[i] - 4.0 * tr_xy_xxxyz[i] * tke_0 - 10.0 * tr_xxy_xxyz[i] * tbe_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xxzz[i] = 2.0 * tr_y_xxzz[i] + 4.0 * tr_xy_xzz[i] - 4.0 * tr_xy_xxxzz[i] * tke_0 - 10.0 * tr_xxy_xxzz[i] * tbe_0 - 4.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xyyy[i] = 2.0 * tr_y_xyyy[i] + 2.0 * tr_xy_yyy[i] - 4.0 * tr_xy_xxyyy[i] * tke_0 - 10.0 * tr_xxy_xyyy[i] * tbe_0 - 2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xyyz[i] = 2.0 * tr_y_xyyz[i] + 2.0 * tr_xy_yyz[i] - 4.0 * tr_xy_xxyyz[i] * tke_0 - 10.0 * tr_xxy_xyyz[i] * tbe_0 - 2.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xyzz[i] = 2.0 * tr_y_xyzz[i] + 2.0 * tr_xy_yzz[i] - 4.0 * tr_xy_xxyzz[i] * tke_0 - 10.0 * tr_xxy_xyzz[i] * tbe_0 - 2.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xzzz[i] = 2.0 * tr_y_xzzz[i] + 2.0 * tr_xy_zzz[i] - 4.0 * tr_xy_xxzzz[i] * tke_0 - 10.0 * tr_xxy_xzzz[i] * tbe_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_yyyy[i] = 2.0 * tr_y_yyyy[i] - 4.0 * tr_xy_xyyyy[i] * tke_0 - 10.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_yyyz[i] = 2.0 * tr_y_yyyz[i] - 4.0 * tr_xy_xyyyz[i] * tke_0 - 10.0 * tr_xxy_yyyz[i] * tbe_0 + 4.0 * tr_xxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_yyzz[i] = 2.0 * tr_y_yyzz[i] - 4.0 * tr_xy_xyyzz[i] * tke_0 - 10.0 * tr_xxy_yyzz[i] * tbe_0 + 4.0 * tr_xxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_yzzz[i] = 2.0 * tr_y_yzzz[i] - 4.0 * tr_xy_xyzzz[i] * tke_0 - 10.0 * tr_xxy_yzzz[i] * tbe_0 + 4.0 * tr_xxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_zzzz[i] = 2.0 * tr_y_zzzz[i] - 4.0 * tr_xy_xzzzz[i] * tke_0 - 10.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-45 components of targeted buffer : FG

    auto tr_x_0_x_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 30);

    auto tr_x_0_x_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 31);

    auto tr_x_0_x_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 32);

    auto tr_x_0_x_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 33);

    auto tr_x_0_x_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 34);

    auto tr_x_0_x_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 35);

    auto tr_x_0_x_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 36);

    auto tr_x_0_x_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 37);

    auto tr_x_0_x_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 38);

    auto tr_x_0_x_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 39);

    auto tr_x_0_x_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 40);

    auto tr_x_0_x_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 41);

    auto tr_x_0_x_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 42);

    auto tr_x_0_x_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 43);

    auto tr_x_0_x_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 44);

    #pragma omp simd aligned(tr_x_0_x_xxz_xxxx, tr_x_0_x_xxz_xxxy, tr_x_0_x_xxz_xxxz, tr_x_0_x_xxz_xxyy, tr_x_0_x_xxz_xxyz, tr_x_0_x_xxz_xxzz, tr_x_0_x_xxz_xyyy, tr_x_0_x_xxz_xyyz, tr_x_0_x_xxz_xyzz, tr_x_0_x_xxz_xzzz, tr_x_0_x_xxz_yyyy, tr_x_0_x_xxz_yyyz, tr_x_0_x_xxz_yyzz, tr_x_0_x_xxz_yzzz, tr_x_0_x_xxz_zzzz, tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xzzz, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yzzz, tr_xxxxz_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxx, tr_xxxz_xxxxy, tr_xxxz_xxxxz, tr_xxxz_xxxyy, tr_xxxz_xxxyz, tr_xxxz_xxxzz, tr_xxxz_xxy, tr_xxxz_xxyyy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xxzzz, tr_xxxz_xyy, tr_xxxz_xyyyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_xzzzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxz_xxxx[i] = 2.0 * tr_z_xxxx[i] + 8.0 * tr_xz_xxx[i] - 4.0 * tr_xz_xxxxx[i] * tke_0 - 10.0 * tr_xxz_xxxx[i] * tbe_0 - 8.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xxxy[i] = 2.0 * tr_z_xxxy[i] + 6.0 * tr_xz_xxy[i] - 4.0 * tr_xz_xxxxy[i] * tke_0 - 10.0 * tr_xxz_xxxy[i] * tbe_0 - 6.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xxxz[i] = 2.0 * tr_z_xxxz[i] + 6.0 * tr_xz_xxz[i] - 4.0 * tr_xz_xxxxz[i] * tke_0 - 10.0 * tr_xxz_xxxz[i] * tbe_0 - 6.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xxyy[i] = 2.0 * tr_z_xxyy[i] + 4.0 * tr_xz_xyy[i] - 4.0 * tr_xz_xxxyy[i] * tke_0 - 10.0 * tr_xxz_xxyy[i] * tbe_0 - 4.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xxyz[i] = 2.0 * tr_z_xxyz[i] + 4.0 * tr_xz_xyz[i] - 4.0 * tr_xz_xxxyz[i] * tke_0 - 10.0 * tr_xxz_xxyz[i] * tbe_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xxzz[i] = 2.0 * tr_z_xxzz[i] + 4.0 * tr_xz_xzz[i] - 4.0 * tr_xz_xxxzz[i] * tke_0 - 10.0 * tr_xxz_xxzz[i] * tbe_0 - 4.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xyyy[i] = 2.0 * tr_z_xyyy[i] + 2.0 * tr_xz_yyy[i] - 4.0 * tr_xz_xxyyy[i] * tke_0 - 10.0 * tr_xxz_xyyy[i] * tbe_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xyyz[i] = 2.0 * tr_z_xyyz[i] + 2.0 * tr_xz_yyz[i] - 4.0 * tr_xz_xxyyz[i] * tke_0 - 10.0 * tr_xxz_xyyz[i] * tbe_0 - 2.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xyzz[i] = 2.0 * tr_z_xyzz[i] + 2.0 * tr_xz_yzz[i] - 4.0 * tr_xz_xxyzz[i] * tke_0 - 10.0 * tr_xxz_xyzz[i] * tbe_0 - 2.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xzzz[i] = 2.0 * tr_z_xzzz[i] + 2.0 * tr_xz_zzz[i] - 4.0 * tr_xz_xxzzz[i] * tke_0 - 10.0 * tr_xxz_xzzz[i] * tbe_0 - 2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_yyyy[i] = 2.0 * tr_z_yyyy[i] - 4.0 * tr_xz_xyyyy[i] * tke_0 - 10.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_yyyz[i] = 2.0 * tr_z_yyyz[i] - 4.0 * tr_xz_xyyyz[i] * tke_0 - 10.0 * tr_xxz_yyyz[i] * tbe_0 + 4.0 * tr_xxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_yyzz[i] = 2.0 * tr_z_yyzz[i] - 4.0 * tr_xz_xyyzz[i] * tke_0 - 10.0 * tr_xxz_yyzz[i] * tbe_0 + 4.0 * tr_xxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_yzzz[i] = 2.0 * tr_z_yzzz[i] - 4.0 * tr_xz_xyzzz[i] * tke_0 - 10.0 * tr_xxz_yzzz[i] * tbe_0 + 4.0 * tr_xxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_zzzz[i] = 2.0 * tr_z_zzzz[i] - 4.0 * tr_xz_xzzzz[i] * tke_0 - 10.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 45-60 components of targeted buffer : FG

    auto tr_x_0_x_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 45);

    auto tr_x_0_x_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 46);

    auto tr_x_0_x_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 47);

    auto tr_x_0_x_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 48);

    auto tr_x_0_x_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 49);

    auto tr_x_0_x_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 50);

    auto tr_x_0_x_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 51);

    auto tr_x_0_x_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 52);

    auto tr_x_0_x_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 53);

    auto tr_x_0_x_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 54);

    auto tr_x_0_x_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 55);

    auto tr_x_0_x_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 56);

    auto tr_x_0_x_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 57);

    auto tr_x_0_x_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 58);

    auto tr_x_0_x_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 59);

    #pragma omp simd aligned(tr_x_0_x_xyy_xxxx, tr_x_0_x_xyy_xxxy, tr_x_0_x_xyy_xxxz, tr_x_0_x_xyy_xxyy, tr_x_0_x_xyy_xxyz, tr_x_0_x_xyy_xxzz, tr_x_0_x_xyy_xyyy, tr_x_0_x_xyy_xyyz, tr_x_0_x_xyy_xyzz, tr_x_0_x_xyy_xzzz, tr_x_0_x_xyy_yyyy, tr_x_0_x_xyy_yyyz, tr_x_0_x_xyy_yyzz, tr_x_0_x_xyy_yzzz, tr_x_0_x_xyy_zzzz, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xzzz, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yzzz, tr_xxxyy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxx, tr_xxyy_xxxxy, tr_xxyy_xxxxz, tr_xxyy_xxxyy, tr_xxyy_xxxyz, tr_xxyy_xxxzz, tr_xxyy_xxy, tr_xxyy_xxyyy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xxzzz, tr_xxyy_xyy, tr_xxyy_xyyyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_xzzzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyy_xxxx[i] = 4.0 * tr_yy_xxx[i] - 2.0 * tr_yy_xxxxx[i] * tke_0 - 6.0 * tr_xyy_xxxx[i] * tbe_0 - 8.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xxxy[i] = 3.0 * tr_yy_xxy[i] - 2.0 * tr_yy_xxxxy[i] * tke_0 - 6.0 * tr_xyy_xxxy[i] * tbe_0 - 6.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xxxz[i] = 3.0 * tr_yy_xxz[i] - 2.0 * tr_yy_xxxxz[i] * tke_0 - 6.0 * tr_xyy_xxxz[i] * tbe_0 - 6.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xxyy[i] = 2.0 * tr_yy_xyy[i] - 2.0 * tr_yy_xxxyy[i] * tke_0 - 6.0 * tr_xyy_xxyy[i] * tbe_0 - 4.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xxyz[i] = 2.0 * tr_yy_xyz[i] - 2.0 * tr_yy_xxxyz[i] * tke_0 - 6.0 * tr_xyy_xxyz[i] * tbe_0 - 4.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xxzz[i] = 2.0 * tr_yy_xzz[i] - 2.0 * tr_yy_xxxzz[i] * tke_0 - 6.0 * tr_xyy_xxzz[i] * tbe_0 - 4.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xyyy[i] = tr_yy_yyy[i] - 2.0 * tr_yy_xxyyy[i] * tke_0 - 6.0 * tr_xyy_xyyy[i] * tbe_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xyyz[i] = tr_yy_yyz[i] - 2.0 * tr_yy_xxyyz[i] * tke_0 - 6.0 * tr_xyy_xyyz[i] * tbe_0 - 2.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xyzz[i] = tr_yy_yzz[i] - 2.0 * tr_yy_xxyzz[i] * tke_0 - 6.0 * tr_xyy_xyzz[i] * tbe_0 - 2.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xzzz[i] = tr_yy_zzz[i] - 2.0 * tr_yy_xxzzz[i] * tke_0 - 6.0 * tr_xyy_xzzz[i] * tbe_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_yyyy[i] = -2.0 * tr_yy_xyyyy[i] * tke_0 - 6.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_yyyz[i] = -2.0 * tr_yy_xyyyz[i] * tke_0 - 6.0 * tr_xyy_yyyz[i] * tbe_0 + 4.0 * tr_xxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_yyzz[i] = -2.0 * tr_yy_xyyzz[i] * tke_0 - 6.0 * tr_xyy_yyzz[i] * tbe_0 + 4.0 * tr_xxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_yzzz[i] = -2.0 * tr_yy_xyzzz[i] * tke_0 - 6.0 * tr_xyy_yzzz[i] * tbe_0 + 4.0 * tr_xxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_zzzz[i] = -2.0 * tr_yy_xzzzz[i] * tke_0 - 6.0 * tr_xyy_zzzz[i] * tbe_0 + 4.0 * tr_xxyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-75 components of targeted buffer : FG

    auto tr_x_0_x_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 60);

    auto tr_x_0_x_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 61);

    auto tr_x_0_x_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 62);

    auto tr_x_0_x_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 63);

    auto tr_x_0_x_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 64);

    auto tr_x_0_x_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 65);

    auto tr_x_0_x_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 66);

    auto tr_x_0_x_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 67);

    auto tr_x_0_x_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 68);

    auto tr_x_0_x_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 69);

    auto tr_x_0_x_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 70);

    auto tr_x_0_x_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 71);

    auto tr_x_0_x_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 72);

    auto tr_x_0_x_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 73);

    auto tr_x_0_x_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 74);

    #pragma omp simd aligned(tr_x_0_x_xyz_xxxx, tr_x_0_x_xyz_xxxy, tr_x_0_x_xyz_xxxz, tr_x_0_x_xyz_xxyy, tr_x_0_x_xyz_xxyz, tr_x_0_x_xyz_xxzz, tr_x_0_x_xyz_xyyy, tr_x_0_x_xyz_xyyz, tr_x_0_x_xyz_xyzz, tr_x_0_x_xyz_xzzz, tr_x_0_x_xyz_yyyy, tr_x_0_x_xyz_yyyz, tr_x_0_x_xyz_yyzz, tr_x_0_x_xyz_yzzz, tr_x_0_x_xyz_zzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxx, tr_xxyz_xxxxy, tr_xxyz_xxxxz, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyz_xxxx[i] = 4.0 * tr_yz_xxx[i] - 2.0 * tr_yz_xxxxx[i] * tke_0 - 6.0 * tr_xyz_xxxx[i] * tbe_0 - 8.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xxxy[i] = 3.0 * tr_yz_xxy[i] - 2.0 * tr_yz_xxxxy[i] * tke_0 - 6.0 * tr_xyz_xxxy[i] * tbe_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xxxz[i] = 3.0 * tr_yz_xxz[i] - 2.0 * tr_yz_xxxxz[i] * tke_0 - 6.0 * tr_xyz_xxxz[i] * tbe_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xxyy[i] = 2.0 * tr_yz_xyy[i] - 2.0 * tr_yz_xxxyy[i] * tke_0 - 6.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xxyz[i] = 2.0 * tr_yz_xyz[i] - 2.0 * tr_yz_xxxyz[i] * tke_0 - 6.0 * tr_xyz_xxyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xxzz[i] = 2.0 * tr_yz_xzz[i] - 2.0 * tr_yz_xxxzz[i] * tke_0 - 6.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xyyy[i] = tr_yz_yyy[i] - 2.0 * tr_yz_xxyyy[i] * tke_0 - 6.0 * tr_xyz_xyyy[i] * tbe_0 - 2.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xyyz[i] = tr_yz_yyz[i] - 2.0 * tr_yz_xxyyz[i] * tke_0 - 6.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xyzz[i] = tr_yz_yzz[i] - 2.0 * tr_yz_xxyzz[i] * tke_0 - 6.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xzzz[i] = tr_yz_zzz[i] - 2.0 * tr_yz_xxzzz[i] * tke_0 - 6.0 * tr_xyz_xzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_yyyy[i] = -2.0 * tr_yz_xyyyy[i] * tke_0 - 6.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_yyyz[i] = -2.0 * tr_yz_xyyyz[i] * tke_0 - 6.0 * tr_xyz_yyyz[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_yyzz[i] = -2.0 * tr_yz_xyyzz[i] * tke_0 - 6.0 * tr_xyz_yyzz[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_yzzz[i] = -2.0 * tr_yz_xyzzz[i] * tke_0 - 6.0 * tr_xyz_yzzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_zzzz[i] = -2.0 * tr_yz_xzzzz[i] * tke_0 - 6.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 75-90 components of targeted buffer : FG

    auto tr_x_0_x_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 75);

    auto tr_x_0_x_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 76);

    auto tr_x_0_x_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 77);

    auto tr_x_0_x_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 78);

    auto tr_x_0_x_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 79);

    auto tr_x_0_x_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 80);

    auto tr_x_0_x_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 81);

    auto tr_x_0_x_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 82);

    auto tr_x_0_x_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 83);

    auto tr_x_0_x_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 84);

    auto tr_x_0_x_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 85);

    auto tr_x_0_x_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 86);

    auto tr_x_0_x_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 87);

    auto tr_x_0_x_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 88);

    auto tr_x_0_x_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 89);

    #pragma omp simd aligned(tr_x_0_x_xzz_xxxx, tr_x_0_x_xzz_xxxy, tr_x_0_x_xzz_xxxz, tr_x_0_x_xzz_xxyy, tr_x_0_x_xzz_xxyz, tr_x_0_x_xzz_xxzz, tr_x_0_x_xzz_xyyy, tr_x_0_x_xzz_xyyz, tr_x_0_x_xzz_xyzz, tr_x_0_x_xzz_xzzz, tr_x_0_x_xzz_yyyy, tr_x_0_x_xzz_yyyz, tr_x_0_x_xzz_yyzz, tr_x_0_x_xzz_yzzz, tr_x_0_x_xzz_zzzz, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xzzz, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yzzz, tr_xxxzz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxx, tr_xxzz_xxxxy, tr_xxzz_xxxxz, tr_xxzz_xxxyy, tr_xxzz_xxxyz, tr_xxzz_xxxzz, tr_xxzz_xxy, tr_xxzz_xxyyy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xxzzz, tr_xxzz_xyy, tr_xxzz_xyyyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_xzzzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xzz_xxxx[i] = 4.0 * tr_zz_xxx[i] - 2.0 * tr_zz_xxxxx[i] * tke_0 - 6.0 * tr_xzz_xxxx[i] * tbe_0 - 8.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xxxy[i] = 3.0 * tr_zz_xxy[i] - 2.0 * tr_zz_xxxxy[i] * tke_0 - 6.0 * tr_xzz_xxxy[i] * tbe_0 - 6.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xxxz[i] = 3.0 * tr_zz_xxz[i] - 2.0 * tr_zz_xxxxz[i] * tke_0 - 6.0 * tr_xzz_xxxz[i] * tbe_0 - 6.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xxyy[i] = 2.0 * tr_zz_xyy[i] - 2.0 * tr_zz_xxxyy[i] * tke_0 - 6.0 * tr_xzz_xxyy[i] * tbe_0 - 4.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xxyz[i] = 2.0 * tr_zz_xyz[i] - 2.0 * tr_zz_xxxyz[i] * tke_0 - 6.0 * tr_xzz_xxyz[i] * tbe_0 - 4.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xxzz[i] = 2.0 * tr_zz_xzz[i] - 2.0 * tr_zz_xxxzz[i] * tke_0 - 6.0 * tr_xzz_xxzz[i] * tbe_0 - 4.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xyyy[i] = tr_zz_yyy[i] - 2.0 * tr_zz_xxyyy[i] * tke_0 - 6.0 * tr_xzz_xyyy[i] * tbe_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xyyz[i] = tr_zz_yyz[i] - 2.0 * tr_zz_xxyyz[i] * tke_0 - 6.0 * tr_xzz_xyyz[i] * tbe_0 - 2.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xyzz[i] = tr_zz_yzz[i] - 2.0 * tr_zz_xxyzz[i] * tke_0 - 6.0 * tr_xzz_xyzz[i] * tbe_0 - 2.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xzzz[i] = tr_zz_zzz[i] - 2.0 * tr_zz_xxzzz[i] * tke_0 - 6.0 * tr_xzz_xzzz[i] * tbe_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_yyyy[i] = -2.0 * tr_zz_xyyyy[i] * tke_0 - 6.0 * tr_xzz_yyyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_yyyz[i] = -2.0 * tr_zz_xyyyz[i] * tke_0 - 6.0 * tr_xzz_yyyz[i] * tbe_0 + 4.0 * tr_xxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_yyzz[i] = -2.0 * tr_zz_xyyzz[i] * tke_0 - 6.0 * tr_xzz_yyzz[i] * tbe_0 + 4.0 * tr_xxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_yzzz[i] = -2.0 * tr_zz_xyzzz[i] * tke_0 - 6.0 * tr_xzz_yzzz[i] * tbe_0 + 4.0 * tr_xxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_zzzz[i] = -2.0 * tr_zz_xzzzz[i] * tke_0 - 6.0 * tr_xzz_zzzz[i] * tbe_0 + 4.0 * tr_xxzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-105 components of targeted buffer : FG

    auto tr_x_0_x_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 90);

    auto tr_x_0_x_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 91);

    auto tr_x_0_x_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 92);

    auto tr_x_0_x_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 93);

    auto tr_x_0_x_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 94);

    auto tr_x_0_x_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 95);

    auto tr_x_0_x_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 96);

    auto tr_x_0_x_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 97);

    auto tr_x_0_x_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 98);

    auto tr_x_0_x_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 99);

    auto tr_x_0_x_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 100);

    auto tr_x_0_x_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 101);

    auto tr_x_0_x_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 102);

    auto tr_x_0_x_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 103);

    auto tr_x_0_x_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 104);

    #pragma omp simd aligned(tr_x_0_x_yyy_xxxx, tr_x_0_x_yyy_xxxy, tr_x_0_x_yyy_xxxz, tr_x_0_x_yyy_xxyy, tr_x_0_x_yyy_xxyz, tr_x_0_x_yyy_xxzz, tr_x_0_x_yyy_xyyy, tr_x_0_x_yyy_xyyz, tr_x_0_x_yyy_xyzz, tr_x_0_x_yyy_xzzz, tr_x_0_x_yyy_yyyy, tr_x_0_x_yyy_yyyz, tr_x_0_x_yyy_yyzz, tr_x_0_x_yyy_yzzz, tr_x_0_x_yyy_zzzz, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xzzz, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yzzz, tr_xxyyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxx, tr_xyyy_xxxxy, tr_xyyy_xxxxz, tr_xyyy_xxxyy, tr_xyyy_xxxyz, tr_xyyy_xxxzz, tr_xyyy_xxy, tr_xyyy_xxyyy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xxzzz, tr_xyyy_xyy, tr_xyyy_xyyyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_xzzzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyy_xxxx[i] = -2.0 * tr_yyy_xxxx[i] * tbe_0 - 8.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xxxy[i] = -2.0 * tr_yyy_xxxy[i] * tbe_0 - 6.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xxxz[i] = -2.0 * tr_yyy_xxxz[i] * tbe_0 - 6.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xxyy[i] = -2.0 * tr_yyy_xxyy[i] * tbe_0 - 4.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xxyz[i] = -2.0 * tr_yyy_xxyz[i] * tbe_0 - 4.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xxzz[i] = -2.0 * tr_yyy_xxzz[i] * tbe_0 - 4.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xyyy[i] = -2.0 * tr_yyy_xyyy[i] * tbe_0 - 2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xyyz[i] = -2.0 * tr_yyy_xyyz[i] * tbe_0 - 2.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xyzz[i] = -2.0 * tr_yyy_xyzz[i] * tbe_0 - 2.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xzzz[i] = -2.0 * tr_yyy_xzzz[i] * tbe_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_yyyy[i] = -2.0 * tr_yyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_yyyz[i] = -2.0 * tr_yyy_yyyz[i] * tbe_0 + 4.0 * tr_xyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_yyzz[i] = -2.0 * tr_yyy_yyzz[i] * tbe_0 + 4.0 * tr_xyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_yzzz[i] = -2.0 * tr_yyy_yzzz[i] * tbe_0 + 4.0 * tr_xyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_zzzz[i] = -2.0 * tr_yyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 105-120 components of targeted buffer : FG

    auto tr_x_0_x_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 105);

    auto tr_x_0_x_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 106);

    auto tr_x_0_x_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 107);

    auto tr_x_0_x_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 108);

    auto tr_x_0_x_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 109);

    auto tr_x_0_x_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 110);

    auto tr_x_0_x_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 111);

    auto tr_x_0_x_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 112);

    auto tr_x_0_x_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 113);

    auto tr_x_0_x_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 114);

    auto tr_x_0_x_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 115);

    auto tr_x_0_x_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 116);

    auto tr_x_0_x_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 117);

    auto tr_x_0_x_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 118);

    auto tr_x_0_x_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 119);

    #pragma omp simd aligned(tr_x_0_x_yyz_xxxx, tr_x_0_x_yyz_xxxy, tr_x_0_x_yyz_xxxz, tr_x_0_x_yyz_xxyy, tr_x_0_x_yyz_xxyz, tr_x_0_x_yyz_xxzz, tr_x_0_x_yyz_xyyy, tr_x_0_x_yyz_xyyz, tr_x_0_x_yyz_xyzz, tr_x_0_x_yyz_xzzz, tr_x_0_x_yyz_yyyy, tr_x_0_x_yyz_yyyz, tr_x_0_x_yyz_yyzz, tr_x_0_x_yyz_yzzz, tr_x_0_x_yyz_zzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxx, tr_xyyz_xxxxy, tr_xyyz_xxxxz, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyz_xxxx[i] = -2.0 * tr_yyz_xxxx[i] * tbe_0 - 8.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xxxy[i] = -2.0 * tr_yyz_xxxy[i] * tbe_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xxxz[i] = -2.0 * tr_yyz_xxxz[i] * tbe_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xxyy[i] = -2.0 * tr_yyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xxyz[i] = -2.0 * tr_yyz_xxyz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xxzz[i] = -2.0 * tr_yyz_xxzz[i] * tbe_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xyyy[i] = -2.0 * tr_yyz_xyyy[i] * tbe_0 - 2.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xyyz[i] = -2.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xyzz[i] = -2.0 * tr_yyz_xyzz[i] * tbe_0 - 2.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xzzz[i] = -2.0 * tr_yyz_xzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_yyyy[i] = -2.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_yyyz[i] = -2.0 * tr_yyz_yyyz[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_yyzz[i] = -2.0 * tr_yyz_yyzz[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_yzzz[i] = -2.0 * tr_yyz_yzzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_zzzz[i] = -2.0 * tr_yyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-135 components of targeted buffer : FG

    auto tr_x_0_x_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 120);

    auto tr_x_0_x_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 121);

    auto tr_x_0_x_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 122);

    auto tr_x_0_x_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 123);

    auto tr_x_0_x_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 124);

    auto tr_x_0_x_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 125);

    auto tr_x_0_x_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 126);

    auto tr_x_0_x_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 127);

    auto tr_x_0_x_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 128);

    auto tr_x_0_x_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 129);

    auto tr_x_0_x_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 130);

    auto tr_x_0_x_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 131);

    auto tr_x_0_x_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 132);

    auto tr_x_0_x_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 133);

    auto tr_x_0_x_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 134);

    #pragma omp simd aligned(tr_x_0_x_yzz_xxxx, tr_x_0_x_yzz_xxxy, tr_x_0_x_yzz_xxxz, tr_x_0_x_yzz_xxyy, tr_x_0_x_yzz_xxyz, tr_x_0_x_yzz_xxzz, tr_x_0_x_yzz_xyyy, tr_x_0_x_yzz_xyyz, tr_x_0_x_yzz_xyzz, tr_x_0_x_yzz_xzzz, tr_x_0_x_yzz_yyyy, tr_x_0_x_yzz_yyyz, tr_x_0_x_yzz_yyzz, tr_x_0_x_yzz_yzzz, tr_x_0_x_yzz_zzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxx, tr_xyzz_xxxxy, tr_xyzz_xxxxz, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yzz_xxxx[i] = -2.0 * tr_yzz_xxxx[i] * tbe_0 - 8.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xxxy[i] = -2.0 * tr_yzz_xxxy[i] * tbe_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xxxz[i] = -2.0 * tr_yzz_xxxz[i] * tbe_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xxyy[i] = -2.0 * tr_yzz_xxyy[i] * tbe_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xxyz[i] = -2.0 * tr_yzz_xxyz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xxzz[i] = -2.0 * tr_yzz_xxzz[i] * tbe_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xyyy[i] = -2.0 * tr_yzz_xyyy[i] * tbe_0 - 2.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xyyz[i] = -2.0 * tr_yzz_xyyz[i] * tbe_0 - 2.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xyzz[i] = -2.0 * tr_yzz_xyzz[i] * tbe_0 - 2.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xzzz[i] = -2.0 * tr_yzz_xzzz[i] * tbe_0 - 2.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_yyyy[i] = -2.0 * tr_yzz_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_yyyz[i] = -2.0 * tr_yzz_yyyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_yyzz[i] = -2.0 * tr_yzz_yyzz[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_yzzz[i] = -2.0 * tr_yzz_yzzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_zzzz[i] = -2.0 * tr_yzz_zzzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 135-150 components of targeted buffer : FG

    auto tr_x_0_x_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 135);

    auto tr_x_0_x_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 136);

    auto tr_x_0_x_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 137);

    auto tr_x_0_x_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 138);

    auto tr_x_0_x_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 139);

    auto tr_x_0_x_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 140);

    auto tr_x_0_x_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 141);

    auto tr_x_0_x_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 142);

    auto tr_x_0_x_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 143);

    auto tr_x_0_x_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 144);

    auto tr_x_0_x_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 145);

    auto tr_x_0_x_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 146);

    auto tr_x_0_x_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 147);

    auto tr_x_0_x_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 148);

    auto tr_x_0_x_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 149);

    #pragma omp simd aligned(tr_x_0_x_zzz_xxxx, tr_x_0_x_zzz_xxxy, tr_x_0_x_zzz_xxxz, tr_x_0_x_zzz_xxyy, tr_x_0_x_zzz_xxyz, tr_x_0_x_zzz_xxzz, tr_x_0_x_zzz_xyyy, tr_x_0_x_zzz_xyyz, tr_x_0_x_zzz_xyzz, tr_x_0_x_zzz_xzzz, tr_x_0_x_zzz_yyyy, tr_x_0_x_zzz_yyyz, tr_x_0_x_zzz_yyzz, tr_x_0_x_zzz_yzzz, tr_x_0_x_zzz_zzzz, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xzzz, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yzzz, tr_xxzzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxx, tr_xzzz_xxxxy, tr_xzzz_xxxxz, tr_xzzz_xxxyy, tr_xzzz_xxxyz, tr_xzzz_xxxzz, tr_xzzz_xxy, tr_xzzz_xxyyy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xxzzz, tr_xzzz_xyy, tr_xzzz_xyyyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_xzzzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zzz_xxxx[i] = -2.0 * tr_zzz_xxxx[i] * tbe_0 - 8.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xxxy[i] = -2.0 * tr_zzz_xxxy[i] * tbe_0 - 6.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xxxz[i] = -2.0 * tr_zzz_xxxz[i] * tbe_0 - 6.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xxyy[i] = -2.0 * tr_zzz_xxyy[i] * tbe_0 - 4.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xxyz[i] = -2.0 * tr_zzz_xxyz[i] * tbe_0 - 4.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xxzz[i] = -2.0 * tr_zzz_xxzz[i] * tbe_0 - 4.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xyyy[i] = -2.0 * tr_zzz_xyyy[i] * tbe_0 - 2.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xyyz[i] = -2.0 * tr_zzz_xyyz[i] * tbe_0 - 2.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xyzz[i] = -2.0 * tr_zzz_xyzz[i] * tbe_0 - 2.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xzzz[i] = -2.0 * tr_zzz_xzzz[i] * tbe_0 - 2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_yyyy[i] = -2.0 * tr_zzz_yyyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_yyyz[i] = -2.0 * tr_zzz_yyyz[i] * tbe_0 + 4.0 * tr_xzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_yyzz[i] = -2.0 * tr_zzz_yyzz[i] * tbe_0 + 4.0 * tr_xzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_yzzz[i] = -2.0 * tr_zzz_yzzz[i] * tbe_0 + 4.0 * tr_xzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_zzzz[i] = -2.0 * tr_zzz_zzzz[i] * tbe_0 + 4.0 * tr_xzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-165 components of targeted buffer : FG

    auto tr_x_0_y_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 150);

    auto tr_x_0_y_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 151);

    auto tr_x_0_y_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 152);

    auto tr_x_0_y_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 153);

    auto tr_x_0_y_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 154);

    auto tr_x_0_y_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 155);

    auto tr_x_0_y_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 156);

    auto tr_x_0_y_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 157);

    auto tr_x_0_y_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 158);

    auto tr_x_0_y_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 159);

    auto tr_x_0_y_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 160);

    auto tr_x_0_y_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 161);

    auto tr_x_0_y_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 162);

    auto tr_x_0_y_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 163);

    auto tr_x_0_y_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 164);

    #pragma omp simd aligned(tr_x_0_y_xxx_xxxx, tr_x_0_y_xxx_xxxy, tr_x_0_y_xxx_xxxz, tr_x_0_y_xxx_xxyy, tr_x_0_y_xxx_xxyz, tr_x_0_y_xxx_xxzz, tr_x_0_y_xxx_xyyy, tr_x_0_y_xxx_xyyz, tr_x_0_y_xxx_xyzz, tr_x_0_y_xxx_xzzz, tr_x_0_y_xxx_yyyy, tr_x_0_y_xxx_yyyz, tr_x_0_y_xxx_yyzz, tr_x_0_y_xxx_yzzz, tr_x_0_y_xxx_zzzz, tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xxxx_xxx, tr_xxxx_xxxxy, tr_xxxx_xxxyy, tr_xxxx_xxxyz, tr_xxxx_xxy, tr_xxxx_xxyyy, tr_xxxx_xxyyz, tr_xxxx_xxyzz, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyyyy, tr_xxxx_xyyyz, tr_xxxx_xyyzz, tr_xxxx_xyz, tr_xxxx_xyzzz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyyyy, tr_xxxx_yyyyz, tr_xxxx_yyyzz, tr_xxxx_yyz, tr_xxxx_yyzzz, tr_xxxx_yzz, tr_xxxx_yzzzz, tr_xxxx_zzz, tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xzzz, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yzzz, tr_xxxxy_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxx_xxxx[i] = -6.0 * tr_xx_xxxxy[i] * tke_0 - 6.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxxx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xxxy[i] = 3.0 * tr_xx_xxx[i] - 6.0 * tr_xx_xxxyy[i] * tke_0 - 6.0 * tr_xxy_xxxy[i] * tbe_0 - 2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xxxz[i] = -6.0 * tr_xx_xxxyz[i] * tke_0 - 6.0 * tr_xxy_xxxz[i] * tbe_0 + 4.0 * tr_xxxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xxyy[i] = 6.0 * tr_xx_xxy[i] - 6.0 * tr_xx_xxyyy[i] * tke_0 - 6.0 * tr_xxy_xxyy[i] * tbe_0 - 4.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xxyz[i] = 3.0 * tr_xx_xxz[i] - 6.0 * tr_xx_xxyyz[i] * tke_0 - 6.0 * tr_xxy_xxyz[i] * tbe_0 - 2.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xxzz[i] = -6.0 * tr_xx_xxyzz[i] * tke_0 - 6.0 * tr_xxy_xxzz[i] * tbe_0 + 4.0 * tr_xxxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xyyy[i] = 9.0 * tr_xx_xyy[i] - 6.0 * tr_xx_xyyyy[i] * tke_0 - 6.0 * tr_xxy_xyyy[i] * tbe_0 - 6.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xyyz[i] = 6.0 * tr_xx_xyz[i] - 6.0 * tr_xx_xyyyz[i] * tke_0 - 6.0 * tr_xxy_xyyz[i] * tbe_0 - 4.0 * tr_xxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xyzz[i] = 3.0 * tr_xx_xzz[i] - 6.0 * tr_xx_xyyzz[i] * tke_0 - 6.0 * tr_xxy_xyzz[i] * tbe_0 - 2.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xzzz[i] = -6.0 * tr_xx_xyzzz[i] * tke_0 - 6.0 * tr_xxy_xzzz[i] * tbe_0 + 4.0 * tr_xxxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_yyyy[i] = 12.0 * tr_xx_yyy[i] - 6.0 * tr_xx_yyyyy[i] * tke_0 - 6.0 * tr_xxy_yyyy[i] * tbe_0 - 8.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxx_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_yyyz[i] = 9.0 * tr_xx_yyz[i] - 6.0 * tr_xx_yyyyz[i] * tke_0 - 6.0 * tr_xxy_yyyz[i] * tbe_0 - 6.0 * tr_xxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_yyzz[i] = 6.0 * tr_xx_yzz[i] - 6.0 * tr_xx_yyyzz[i] * tke_0 - 6.0 * tr_xxy_yyzz[i] * tbe_0 - 4.0 * tr_xxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_yzzz[i] = 3.0 * tr_xx_zzz[i] - 6.0 * tr_xx_yyzzz[i] * tke_0 - 6.0 * tr_xxy_yzzz[i] * tbe_0 - 2.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_zzzz[i] = -6.0 * tr_xx_yzzzz[i] * tke_0 - 6.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxxx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 165-180 components of targeted buffer : FG

    auto tr_x_0_y_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 165);

    auto tr_x_0_y_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 166);

    auto tr_x_0_y_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 167);

    auto tr_x_0_y_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 168);

    auto tr_x_0_y_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 169);

    auto tr_x_0_y_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 170);

    auto tr_x_0_y_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 171);

    auto tr_x_0_y_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 172);

    auto tr_x_0_y_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 173);

    auto tr_x_0_y_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 174);

    auto tr_x_0_y_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 175);

    auto tr_x_0_y_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 176);

    auto tr_x_0_y_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 177);

    auto tr_x_0_y_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 178);

    auto tr_x_0_y_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 179);

    #pragma omp simd aligned(tr_x_0_y_xxy_xxxx, tr_x_0_y_xxy_xxxy, tr_x_0_y_xxy_xxxz, tr_x_0_y_xxy_xxyy, tr_x_0_y_xxy_xxyz, tr_x_0_y_xxy_xxzz, tr_x_0_y_xxy_xyyy, tr_x_0_y_xxy_xyyz, tr_x_0_y_xxy_xyzz, tr_x_0_y_xxy_xzzz, tr_x_0_y_xxy_yyyy, tr_x_0_y_xxy_yyyz, tr_x_0_y_xxy_yyzz, tr_x_0_y_xxy_yzzz, tr_x_0_y_xxy_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxy, tr_xxxy_xxxyy, tr_xxxy_xxxyz, tr_xxxy_xxy, tr_xxxy_xxyyy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyyyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyyyy, tr_xxxy_yyyyz, tr_xxxy_yyyzz, tr_xxxy_yyz, tr_xxxy_yyzzz, tr_xxxy_yzz, tr_xxxy_yzzzz, tr_xxxy_zzz, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xzzz, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yzzz, tr_xxxyy_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxy_xxxx[i] = 2.0 * tr_x_xxxx[i] - 4.0 * tr_xy_xxxxy[i] * tke_0 - 4.0 * tr_xyy_xxxx[i] * tbe_0 - 2.0 * tr_xxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xxxy[i] = 2.0 * tr_x_xxxy[i] + 2.0 * tr_xy_xxx[i] - 4.0 * tr_xy_xxxyy[i] * tke_0 - 4.0 * tr_xyy_xxxy[i] * tbe_0 - 2.0 * tr_xxx_xxxy[i] * tbe_0 - 2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xxxz[i] = 2.0 * tr_x_xxxz[i] - 4.0 * tr_xy_xxxyz[i] * tke_0 - 4.0 * tr_xyy_xxxz[i] * tbe_0 - 2.0 * tr_xxx_xxxz[i] * tbe_0 + 4.0 * tr_xxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xxyy[i] = 2.0 * tr_x_xxyy[i] + 4.0 * tr_xy_xxy[i] - 4.0 * tr_xy_xxyyy[i] * tke_0 - 4.0 * tr_xyy_xxyy[i] * tbe_0 - 2.0 * tr_xxx_xxyy[i] * tbe_0 - 4.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xxyz[i] = 2.0 * tr_x_xxyz[i] + 2.0 * tr_xy_xxz[i] - 4.0 * tr_xy_xxyyz[i] * tke_0 - 4.0 * tr_xyy_xxyz[i] * tbe_0 - 2.0 * tr_xxx_xxyz[i] * tbe_0 - 2.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xxzz[i] = 2.0 * tr_x_xxzz[i] - 4.0 * tr_xy_xxyzz[i] * tke_0 - 4.0 * tr_xyy_xxzz[i] * tbe_0 - 2.0 * tr_xxx_xxzz[i] * tbe_0 + 4.0 * tr_xxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xyyy[i] = 2.0 * tr_x_xyyy[i] + 6.0 * tr_xy_xyy[i] - 4.0 * tr_xy_xyyyy[i] * tke_0 - 4.0 * tr_xyy_xyyy[i] * tbe_0 - 2.0 * tr_xxx_xyyy[i] * tbe_0 - 6.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xyyz[i] = 2.0 * tr_x_xyyz[i] + 4.0 * tr_xy_xyz[i] - 4.0 * tr_xy_xyyyz[i] * tke_0 - 4.0 * tr_xyy_xyyz[i] * tbe_0 - 2.0 * tr_xxx_xyyz[i] * tbe_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xyzz[i] = 2.0 * tr_x_xyzz[i] + 2.0 * tr_xy_xzz[i] - 4.0 * tr_xy_xyyzz[i] * tke_0 - 4.0 * tr_xyy_xyzz[i] * tbe_0 - 2.0 * tr_xxx_xyzz[i] * tbe_0 - 2.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xzzz[i] = 2.0 * tr_x_xzzz[i] - 4.0 * tr_xy_xyzzz[i] * tke_0 - 4.0 * tr_xyy_xzzz[i] * tbe_0 - 2.0 * tr_xxx_xzzz[i] * tbe_0 + 4.0 * tr_xxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_yyyy[i] = 2.0 * tr_x_yyyy[i] + 8.0 * tr_xy_yyy[i] - 4.0 * tr_xy_yyyyy[i] * tke_0 - 4.0 * tr_xyy_yyyy[i] * tbe_0 - 2.0 * tr_xxx_yyyy[i] * tbe_0 - 8.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_yyyz[i] = 2.0 * tr_x_yyyz[i] + 6.0 * tr_xy_yyz[i] - 4.0 * tr_xy_yyyyz[i] * tke_0 - 4.0 * tr_xyy_yyyz[i] * tbe_0 - 2.0 * tr_xxx_yyyz[i] * tbe_0 - 6.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_yyzz[i] = 2.0 * tr_x_yyzz[i] + 4.0 * tr_xy_yzz[i] - 4.0 * tr_xy_yyyzz[i] * tke_0 - 4.0 * tr_xyy_yyzz[i] * tbe_0 - 2.0 * tr_xxx_yyzz[i] * tbe_0 - 4.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_yzzz[i] = 2.0 * tr_x_yzzz[i] + 2.0 * tr_xy_zzz[i] - 4.0 * tr_xy_yyzzz[i] * tke_0 - 4.0 * tr_xyy_yzzz[i] * tbe_0 - 2.0 * tr_xxx_yzzz[i] * tbe_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_zzzz[i] = 2.0 * tr_x_zzzz[i] - 4.0 * tr_xy_yzzzz[i] * tke_0 - 4.0 * tr_xyy_zzzz[i] * tbe_0 - 2.0 * tr_xxx_zzzz[i] * tbe_0 + 4.0 * tr_xxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-195 components of targeted buffer : FG

    auto tr_x_0_y_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 180);

    auto tr_x_0_y_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 181);

    auto tr_x_0_y_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 182);

    auto tr_x_0_y_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 183);

    auto tr_x_0_y_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 184);

    auto tr_x_0_y_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 185);

    auto tr_x_0_y_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 186);

    auto tr_x_0_y_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 187);

    auto tr_x_0_y_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 188);

    auto tr_x_0_y_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 189);

    auto tr_x_0_y_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 190);

    auto tr_x_0_y_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 191);

    auto tr_x_0_y_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 192);

    auto tr_x_0_y_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 193);

    auto tr_x_0_y_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 194);

    #pragma omp simd aligned(tr_x_0_y_xxz_xxxx, tr_x_0_y_xxz_xxxy, tr_x_0_y_xxz_xxxz, tr_x_0_y_xxz_xxyy, tr_x_0_y_xxz_xxyz, tr_x_0_y_xxz_xxzz, tr_x_0_y_xxz_xyyy, tr_x_0_y_xxz_xyyz, tr_x_0_y_xxz_xyzz, tr_x_0_y_xxz_xzzz, tr_x_0_y_xxz_yyyy, tr_x_0_y_xxz_yyyz, tr_x_0_y_xxz_yyzz, tr_x_0_y_xxz_yzzz, tr_x_0_y_xxz_zzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxy, tr_xxxz_xxxyy, tr_xxxz_xxxyz, tr_xxxz_xxy, tr_xxxz_xxyyy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyyyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyyyy, tr_xxxz_yyyyz, tr_xxxz_yyyzz, tr_xxxz_yyz, tr_xxxz_yyzzz, tr_xxxz_yzz, tr_xxxz_yzzzz, tr_xxxz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxz_xxxx[i] = -4.0 * tr_xz_xxxxy[i] * tke_0 - 4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xxxy[i] = 2.0 * tr_xz_xxx[i] - 4.0 * tr_xz_xxxyy[i] * tke_0 - 4.0 * tr_xyz_xxxy[i] * tbe_0 - 2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xxxz[i] = -4.0 * tr_xz_xxxyz[i] * tke_0 - 4.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xxyy[i] = 4.0 * tr_xz_xxy[i] - 4.0 * tr_xz_xxyyy[i] * tke_0 - 4.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xxyz[i] = 2.0 * tr_xz_xxz[i] - 4.0 * tr_xz_xxyyz[i] * tke_0 - 4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xxzz[i] = -4.0 * tr_xz_xxyzz[i] * tke_0 - 4.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xyyy[i] = 6.0 * tr_xz_xyy[i] - 4.0 * tr_xz_xyyyy[i] * tke_0 - 4.0 * tr_xyz_xyyy[i] * tbe_0 - 6.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xyyz[i] = 4.0 * tr_xz_xyz[i] - 4.0 * tr_xz_xyyyz[i] * tke_0 - 4.0 * tr_xyz_xyyz[i] * tbe_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xyzz[i] = 2.0 * tr_xz_xzz[i] - 4.0 * tr_xz_xyyzz[i] * tke_0 - 4.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xzzz[i] = -4.0 * tr_xz_xyzzz[i] * tke_0 - 4.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_yyyy[i] = 8.0 * tr_xz_yyy[i] - 4.0 * tr_xz_yyyyy[i] * tke_0 - 4.0 * tr_xyz_yyyy[i] * tbe_0 - 8.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_yyyz[i] = 6.0 * tr_xz_yyz[i] - 4.0 * tr_xz_yyyyz[i] * tke_0 - 4.0 * tr_xyz_yyyz[i] * tbe_0 - 6.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_yyzz[i] = 4.0 * tr_xz_yzz[i] - 4.0 * tr_xz_yyyzz[i] * tke_0 - 4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_yzzz[i] = 2.0 * tr_xz_zzz[i] - 4.0 * tr_xz_yyzzz[i] * tke_0 - 4.0 * tr_xyz_yzzz[i] * tbe_0 - 2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_zzzz[i] = -4.0 * tr_xz_yzzzz[i] * tke_0 - 4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 195-210 components of targeted buffer : FG

    auto tr_x_0_y_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 195);

    auto tr_x_0_y_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 196);

    auto tr_x_0_y_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 197);

    auto tr_x_0_y_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 198);

    auto tr_x_0_y_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 199);

    auto tr_x_0_y_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 200);

    auto tr_x_0_y_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 201);

    auto tr_x_0_y_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 202);

    auto tr_x_0_y_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 203);

    auto tr_x_0_y_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 204);

    auto tr_x_0_y_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 205);

    auto tr_x_0_y_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 206);

    auto tr_x_0_y_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 207);

    auto tr_x_0_y_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 208);

    auto tr_x_0_y_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 209);

    #pragma omp simd aligned(tr_x_0_y_xyy_xxxx, tr_x_0_y_xyy_xxxy, tr_x_0_y_xyy_xxxz, tr_x_0_y_xyy_xxyy, tr_x_0_y_xyy_xxyz, tr_x_0_y_xyy_xxzz, tr_x_0_y_xyy_xyyy, tr_x_0_y_xyy_xyyz, tr_x_0_y_xyy_xyzz, tr_x_0_y_xyy_xzzz, tr_x_0_y_xyy_yyyy, tr_x_0_y_xyy_yyyz, tr_x_0_y_xyy_yyzz, tr_x_0_y_xyy_yzzz, tr_x_0_y_xyy_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxy, tr_xxyy_xxxyy, tr_xxyy_xxxyz, tr_xxyy_xxy, tr_xxyy_xxyyy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyyyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyyyy, tr_xxyy_yyyyz, tr_xxyy_yyyzz, tr_xxyy_yyz, tr_xxyy_yyzzz, tr_xxyy_yzz, tr_xxyy_yzzzz, tr_xxyy_zzz, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xzzz, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yzzz, tr_xxyyy_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyy_xxxx[i] = 2.0 * tr_y_xxxx[i] - 2.0 * tr_yy_xxxxy[i] * tke_0 - 2.0 * tr_yyy_xxxx[i] * tbe_0 - 4.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xxxy[i] = 2.0 * tr_y_xxxy[i] + tr_yy_xxx[i] - 2.0 * tr_yy_xxxyy[i] * tke_0 - 2.0 * tr_yyy_xxxy[i] * tbe_0 - 4.0 * tr_xxy_xxxy[i] * tbe_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xxxz[i] = 2.0 * tr_y_xxxz[i] - 2.0 * tr_yy_xxxyz[i] * tke_0 - 2.0 * tr_yyy_xxxz[i] * tbe_0 - 4.0 * tr_xxy_xxxz[i] * tbe_0 + 4.0 * tr_xxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xxyy[i] = 2.0 * tr_y_xxyy[i] + 2.0 * tr_yy_xxy[i] - 2.0 * tr_yy_xxyyy[i] * tke_0 - 2.0 * tr_yyy_xxyy[i] * tbe_0 - 4.0 * tr_xxy_xxyy[i] * tbe_0 - 4.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xxyz[i] = 2.0 * tr_y_xxyz[i] + tr_yy_xxz[i] - 2.0 * tr_yy_xxyyz[i] * tke_0 - 2.0 * tr_yyy_xxyz[i] * tbe_0 - 4.0 * tr_xxy_xxyz[i] * tbe_0 - 2.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xxzz[i] = 2.0 * tr_y_xxzz[i] - 2.0 * tr_yy_xxyzz[i] * tke_0 - 2.0 * tr_yyy_xxzz[i] * tbe_0 - 4.0 * tr_xxy_xxzz[i] * tbe_0 + 4.0 * tr_xxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xyyy[i] = 2.0 * tr_y_xyyy[i] + 3.0 * tr_yy_xyy[i] - 2.0 * tr_yy_xyyyy[i] * tke_0 - 2.0 * tr_yyy_xyyy[i] * tbe_0 - 4.0 * tr_xxy_xyyy[i] * tbe_0 - 6.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xyyz[i] = 2.0 * tr_y_xyyz[i] + 2.0 * tr_yy_xyz[i] - 2.0 * tr_yy_xyyyz[i] * tke_0 - 2.0 * tr_yyy_xyyz[i] * tbe_0 - 4.0 * tr_xxy_xyyz[i] * tbe_0 - 4.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xyzz[i] = 2.0 * tr_y_xyzz[i] + tr_yy_xzz[i] - 2.0 * tr_yy_xyyzz[i] * tke_0 - 2.0 * tr_yyy_xyzz[i] * tbe_0 - 4.0 * tr_xxy_xyzz[i] * tbe_0 - 2.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xzzz[i] = 2.0 * tr_y_xzzz[i] - 2.0 * tr_yy_xyzzz[i] * tke_0 - 2.0 * tr_yyy_xzzz[i] * tbe_0 - 4.0 * tr_xxy_xzzz[i] * tbe_0 + 4.0 * tr_xxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_yyyy[i] = 2.0 * tr_y_yyyy[i] + 4.0 * tr_yy_yyy[i] - 2.0 * tr_yy_yyyyy[i] * tke_0 - 2.0 * tr_yyy_yyyy[i] * tbe_0 - 4.0 * tr_xxy_yyyy[i] * tbe_0 - 8.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_yyyz[i] = 2.0 * tr_y_yyyz[i] + 3.0 * tr_yy_yyz[i] - 2.0 * tr_yy_yyyyz[i] * tke_0 - 2.0 * tr_yyy_yyyz[i] * tbe_0 - 4.0 * tr_xxy_yyyz[i] * tbe_0 - 6.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_yyzz[i] = 2.0 * tr_y_yyzz[i] + 2.0 * tr_yy_yzz[i] - 2.0 * tr_yy_yyyzz[i] * tke_0 - 2.0 * tr_yyy_yyzz[i] * tbe_0 - 4.0 * tr_xxy_yyzz[i] * tbe_0 - 4.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_yzzz[i] = 2.0 * tr_y_yzzz[i] + tr_yy_zzz[i] - 2.0 * tr_yy_yyzzz[i] * tke_0 - 2.0 * tr_yyy_yzzz[i] * tbe_0 - 4.0 * tr_xxy_yzzz[i] * tbe_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_zzzz[i] = 2.0 * tr_y_zzzz[i] - 2.0 * tr_yy_yzzzz[i] * tke_0 - 2.0 * tr_yyy_zzzz[i] * tbe_0 - 4.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-225 components of targeted buffer : FG

    auto tr_x_0_y_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 210);

    auto tr_x_0_y_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 211);

    auto tr_x_0_y_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 212);

    auto tr_x_0_y_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 213);

    auto tr_x_0_y_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 214);

    auto tr_x_0_y_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 215);

    auto tr_x_0_y_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 216);

    auto tr_x_0_y_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 217);

    auto tr_x_0_y_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 218);

    auto tr_x_0_y_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 219);

    auto tr_x_0_y_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 220);

    auto tr_x_0_y_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 221);

    auto tr_x_0_y_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 222);

    auto tr_x_0_y_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 223);

    auto tr_x_0_y_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 224);

    #pragma omp simd aligned(tr_x_0_y_xyz_xxxx, tr_x_0_y_xyz_xxxy, tr_x_0_y_xyz_xxxz, tr_x_0_y_xyz_xxyy, tr_x_0_y_xyz_xxyz, tr_x_0_y_xyz_xxzz, tr_x_0_y_xyz_xyyy, tr_x_0_y_xyz_xyyz, tr_x_0_y_xyz_xyzz, tr_x_0_y_xyz_xzzz, tr_x_0_y_xyz_yyyy, tr_x_0_y_xyz_yyyz, tr_x_0_y_xyz_yyzz, tr_x_0_y_xyz_yzzz, tr_x_0_y_xyz_zzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxy, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyyyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyz_xxxx[i] = tr_z_xxxx[i] - 2.0 * tr_yz_xxxxy[i] * tke_0 - 2.0 * tr_yyz_xxxx[i] * tbe_0 - 2.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xxxy[i] = tr_z_xxxy[i] + tr_yz_xxx[i] - 2.0 * tr_yz_xxxyy[i] * tke_0 - 2.0 * tr_yyz_xxxy[i] * tbe_0 - 2.0 * tr_xxz_xxxy[i] * tbe_0 - 2.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xxxz[i] = tr_z_xxxz[i] - 2.0 * tr_yz_xxxyz[i] * tke_0 - 2.0 * tr_yyz_xxxz[i] * tbe_0 - 2.0 * tr_xxz_xxxz[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xxyy[i] = tr_z_xxyy[i] + 2.0 * tr_yz_xxy[i] - 2.0 * tr_yz_xxyyy[i] * tke_0 - 2.0 * tr_yyz_xxyy[i] * tbe_0 - 2.0 * tr_xxz_xxyy[i] * tbe_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xxyz[i] = tr_z_xxyz[i] + tr_yz_xxz[i] - 2.0 * tr_yz_xxyyz[i] * tke_0 - 2.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xxzz[i] = tr_z_xxzz[i] - 2.0 * tr_yz_xxyzz[i] * tke_0 - 2.0 * tr_yyz_xxzz[i] * tbe_0 - 2.0 * tr_xxz_xxzz[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xyyy[i] = tr_z_xyyy[i] + 3.0 * tr_yz_xyy[i] - 2.0 * tr_yz_xyyyy[i] * tke_0 - 2.0 * tr_yyz_xyyy[i] * tbe_0 - 2.0 * tr_xxz_xyyy[i] * tbe_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xyyz[i] = tr_z_xyyz[i] + 2.0 * tr_yz_xyz[i] - 2.0 * tr_yz_xyyyz[i] * tke_0 - 2.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_xxz_xyyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xyzz[i] = tr_z_xyzz[i] + tr_yz_xzz[i] - 2.0 * tr_yz_xyyzz[i] * tke_0 - 2.0 * tr_yyz_xyzz[i] * tbe_0 - 2.0 * tr_xxz_xyzz[i] * tbe_0 - 2.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xzzz[i] = tr_z_xzzz[i] - 2.0 * tr_yz_xyzzz[i] * tke_0 - 2.0 * tr_yyz_xzzz[i] * tbe_0 - 2.0 * tr_xxz_xzzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_yyyy[i] = tr_z_yyyy[i] + 4.0 * tr_yz_yyy[i] - 2.0 * tr_yz_yyyyy[i] * tke_0 - 2.0 * tr_yyz_yyyy[i] * tbe_0 - 2.0 * tr_xxz_yyyy[i] * tbe_0 - 8.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_yyyz[i] = tr_z_yyyz[i] + 3.0 * tr_yz_yyz[i] - 2.0 * tr_yz_yyyyz[i] * tke_0 - 2.0 * tr_yyz_yyyz[i] * tbe_0 - 2.0 * tr_xxz_yyyz[i] * tbe_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_yyzz[i] = tr_z_yyzz[i] + 2.0 * tr_yz_yzz[i] - 2.0 * tr_yz_yyyzz[i] * tke_0 - 2.0 * tr_yyz_yyzz[i] * tbe_0 - 2.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_yzzz[i] = tr_z_yzzz[i] + tr_yz_zzz[i] - 2.0 * tr_yz_yyzzz[i] * tke_0 - 2.0 * tr_yyz_yzzz[i] * tbe_0 - 2.0 * tr_xxz_yzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_zzzz[i] = tr_z_zzzz[i] - 2.0 * tr_yz_yzzzz[i] * tke_0 - 2.0 * tr_yyz_zzzz[i] * tbe_0 - 2.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 225-240 components of targeted buffer : FG

    auto tr_x_0_y_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 225);

    auto tr_x_0_y_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 226);

    auto tr_x_0_y_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 227);

    auto tr_x_0_y_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 228);

    auto tr_x_0_y_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 229);

    auto tr_x_0_y_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 230);

    auto tr_x_0_y_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 231);

    auto tr_x_0_y_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 232);

    auto tr_x_0_y_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 233);

    auto tr_x_0_y_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 234);

    auto tr_x_0_y_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 235);

    auto tr_x_0_y_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 236);

    auto tr_x_0_y_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 237);

    auto tr_x_0_y_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 238);

    auto tr_x_0_y_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 239);

    #pragma omp simd aligned(tr_x_0_y_xzz_xxxx, tr_x_0_y_xzz_xxxy, tr_x_0_y_xzz_xxxz, tr_x_0_y_xzz_xxyy, tr_x_0_y_xzz_xxyz, tr_x_0_y_xzz_xxzz, tr_x_0_y_xzz_xyyy, tr_x_0_y_xzz_xyyz, tr_x_0_y_xzz_xyzz, tr_x_0_y_xzz_xzzz, tr_x_0_y_xzz_yyyy, tr_x_0_y_xzz_yyyz, tr_x_0_y_xzz_yyzz, tr_x_0_y_xzz_yzzz, tr_x_0_y_xzz_zzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxy, tr_xxzz_xxxyy, tr_xxzz_xxxyz, tr_xxzz_xxy, tr_xxzz_xxyyy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyyyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyyyy, tr_xxzz_yyyyz, tr_xxzz_yyyzz, tr_xxzz_yyz, tr_xxzz_yyzzz, tr_xxzz_yzz, tr_xxzz_yzzzz, tr_xxzz_zzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xzz_xxxx[i] = -2.0 * tr_zz_xxxxy[i] * tke_0 - 2.0 * tr_yzz_xxxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xxxy[i] = tr_zz_xxx[i] - 2.0 * tr_zz_xxxyy[i] * tke_0 - 2.0 * tr_yzz_xxxy[i] * tbe_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xxxz[i] = -2.0 * tr_zz_xxxyz[i] * tke_0 - 2.0 * tr_yzz_xxxz[i] * tbe_0 + 4.0 * tr_xxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xxyy[i] = 2.0 * tr_zz_xxy[i] - 2.0 * tr_zz_xxyyy[i] * tke_0 - 2.0 * tr_yzz_xxyy[i] * tbe_0 - 4.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xxyz[i] = tr_zz_xxz[i] - 2.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_yzz_xxyz[i] * tbe_0 - 2.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xxzz[i] = -2.0 * tr_zz_xxyzz[i] * tke_0 - 2.0 * tr_yzz_xxzz[i] * tbe_0 + 4.0 * tr_xxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xyyy[i] = 3.0 * tr_zz_xyy[i] - 2.0 * tr_zz_xyyyy[i] * tke_0 - 2.0 * tr_yzz_xyyy[i] * tbe_0 - 6.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xyyz[i] = 2.0 * tr_zz_xyz[i] - 2.0 * tr_zz_xyyyz[i] * tke_0 - 2.0 * tr_yzz_xyyz[i] * tbe_0 - 4.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xyzz[i] = tr_zz_xzz[i] - 2.0 * tr_zz_xyyzz[i] * tke_0 - 2.0 * tr_yzz_xyzz[i] * tbe_0 - 2.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xzzz[i] = -2.0 * tr_zz_xyzzz[i] * tke_0 - 2.0 * tr_yzz_xzzz[i] * tbe_0 + 4.0 * tr_xxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_yyyy[i] = 4.0 * tr_zz_yyy[i] - 2.0 * tr_zz_yyyyy[i] * tke_0 - 2.0 * tr_yzz_yyyy[i] * tbe_0 - 8.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_yyyz[i] = 3.0 * tr_zz_yyz[i] - 2.0 * tr_zz_yyyyz[i] * tke_0 - 2.0 * tr_yzz_yyyz[i] * tbe_0 - 6.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_yyzz[i] = 2.0 * tr_zz_yzz[i] - 2.0 * tr_zz_yyyzz[i] * tke_0 - 2.0 * tr_yzz_yyzz[i] * tbe_0 - 4.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_yzzz[i] = tr_zz_zzz[i] - 2.0 * tr_zz_yyzzz[i] * tke_0 - 2.0 * tr_yzz_yzzz[i] * tbe_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_zzzz[i] = -2.0 * tr_zz_yzzzz[i] * tke_0 - 2.0 * tr_yzz_zzzz[i] * tbe_0 + 4.0 * tr_xxzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-255 components of targeted buffer : FG

    auto tr_x_0_y_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 240);

    auto tr_x_0_y_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 241);

    auto tr_x_0_y_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 242);

    auto tr_x_0_y_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 243);

    auto tr_x_0_y_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 244);

    auto tr_x_0_y_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 245);

    auto tr_x_0_y_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 246);

    auto tr_x_0_y_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 247);

    auto tr_x_0_y_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 248);

    auto tr_x_0_y_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 249);

    auto tr_x_0_y_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 250);

    auto tr_x_0_y_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 251);

    auto tr_x_0_y_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 252);

    auto tr_x_0_y_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 253);

    auto tr_x_0_y_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 254);

    #pragma omp simd aligned(tr_x_0_y_yyy_xxxx, tr_x_0_y_yyy_xxxy, tr_x_0_y_yyy_xxxz, tr_x_0_y_yyy_xxyy, tr_x_0_y_yyy_xxyz, tr_x_0_y_yyy_xxzz, tr_x_0_y_yyy_xyyy, tr_x_0_y_yyy_xyyz, tr_x_0_y_yyy_xyzz, tr_x_0_y_yyy_xzzz, tr_x_0_y_yyy_yyyy, tr_x_0_y_yyy_yyyz, tr_x_0_y_yyy_yyzz, tr_x_0_y_yyy_yzzz, tr_x_0_y_yyy_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxy, tr_xyyy_xxxyy, tr_xyyy_xxxyz, tr_xyyy_xxy, tr_xyyy_xxyyy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyyyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyyyy, tr_xyyy_yyyyz, tr_xyyy_yyyzz, tr_xyyy_yyz, tr_xyyy_yyzzz, tr_xyyy_yzz, tr_xyyy_yzzzz, tr_xyyy_zzz, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xzzz, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yzzz, tr_xyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyy_xxxx[i] = -6.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xxxy[i] = -6.0 * tr_xyy_xxxy[i] * tbe_0 - 2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xxxz[i] = -6.0 * tr_xyy_xxxz[i] * tbe_0 + 4.0 * tr_xyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xxyy[i] = -6.0 * tr_xyy_xxyy[i] * tbe_0 - 4.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xxyz[i] = -6.0 * tr_xyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xxzz[i] = -6.0 * tr_xyy_xxzz[i] * tbe_0 + 4.0 * tr_xyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xyyy[i] = -6.0 * tr_xyy_xyyy[i] * tbe_0 - 6.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xyyz[i] = -6.0 * tr_xyy_xyyz[i] * tbe_0 - 4.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xyzz[i] = -6.0 * tr_xyy_xyzz[i] * tbe_0 - 2.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xzzz[i] = -6.0 * tr_xyy_xzzz[i] * tbe_0 + 4.0 * tr_xyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_yyyy[i] = -6.0 * tr_xyy_yyyy[i] * tbe_0 - 8.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_yyyz[i] = -6.0 * tr_xyy_yyyz[i] * tbe_0 - 6.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_yyzz[i] = -6.0 * tr_xyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_yzzz[i] = -6.0 * tr_xyy_yzzz[i] * tbe_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_zzzz[i] = -6.0 * tr_xyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 255-270 components of targeted buffer : FG

    auto tr_x_0_y_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 255);

    auto tr_x_0_y_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 256);

    auto tr_x_0_y_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 257);

    auto tr_x_0_y_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 258);

    auto tr_x_0_y_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 259);

    auto tr_x_0_y_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 260);

    auto tr_x_0_y_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 261);

    auto tr_x_0_y_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 262);

    auto tr_x_0_y_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 263);

    auto tr_x_0_y_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 264);

    auto tr_x_0_y_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 265);

    auto tr_x_0_y_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 266);

    auto tr_x_0_y_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 267);

    auto tr_x_0_y_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 268);

    auto tr_x_0_y_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 269);

    #pragma omp simd aligned(tr_x_0_y_yyz_xxxx, tr_x_0_y_yyz_xxxy, tr_x_0_y_yyz_xxxz, tr_x_0_y_yyz_xxyy, tr_x_0_y_yyz_xxyz, tr_x_0_y_yyz_xxzz, tr_x_0_y_yyz_xyyy, tr_x_0_y_yyz_xyyz, tr_x_0_y_yyz_xyzz, tr_x_0_y_yyz_xzzz, tr_x_0_y_yyz_yyyy, tr_x_0_y_yyz_yyyz, tr_x_0_y_yyz_yyzz, tr_x_0_y_yyz_yzzz, tr_x_0_y_yyz_zzzz, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxy, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyyyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyz_xxxx[i] = -4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xxxy[i] = -4.0 * tr_xyz_xxxy[i] * tbe_0 - 2.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xxxz[i] = -4.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xxyy[i] = -4.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xxyz[i] = -4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xxzz[i] = -4.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xyyy[i] = -4.0 * tr_xyz_xyyy[i] * tbe_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xyyz[i] = -4.0 * tr_xyz_xyyz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xyzz[i] = -4.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xzzz[i] = -4.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_yyyy[i] = -4.0 * tr_xyz_yyyy[i] * tbe_0 - 8.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_yyyz[i] = -4.0 * tr_xyz_yyyz[i] * tbe_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_yyzz[i] = -4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_yzzz[i] = -4.0 * tr_xyz_yzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_zzzz[i] = -4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-285 components of targeted buffer : FG

    auto tr_x_0_y_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 270);

    auto tr_x_0_y_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 271);

    auto tr_x_0_y_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 272);

    auto tr_x_0_y_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 273);

    auto tr_x_0_y_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 274);

    auto tr_x_0_y_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 275);

    auto tr_x_0_y_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 276);

    auto tr_x_0_y_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 277);

    auto tr_x_0_y_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 278);

    auto tr_x_0_y_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 279);

    auto tr_x_0_y_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 280);

    auto tr_x_0_y_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 281);

    auto tr_x_0_y_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 282);

    auto tr_x_0_y_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 283);

    auto tr_x_0_y_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 284);

    #pragma omp simd aligned(tr_x_0_y_yzz_xxxx, tr_x_0_y_yzz_xxxy, tr_x_0_y_yzz_xxxz, tr_x_0_y_yzz_xxyy, tr_x_0_y_yzz_xxyz, tr_x_0_y_yzz_xxzz, tr_x_0_y_yzz_xyyy, tr_x_0_y_yzz_xyyz, tr_x_0_y_yzz_xyzz, tr_x_0_y_yzz_xzzz, tr_x_0_y_yzz_yyyy, tr_x_0_y_yzz_yyyz, tr_x_0_y_yzz_yyzz, tr_x_0_y_yzz_yzzz, tr_x_0_y_yzz_zzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxy, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyyyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yzz_xxxx[i] = -2.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xxxy[i] = -2.0 * tr_xzz_xxxy[i] * tbe_0 - 2.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xxxz[i] = -2.0 * tr_xzz_xxxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xxyy[i] = -2.0 * tr_xzz_xxyy[i] * tbe_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xxyz[i] = -2.0 * tr_xzz_xxyz[i] * tbe_0 - 2.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xxzz[i] = -2.0 * tr_xzz_xxzz[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xyyy[i] = -2.0 * tr_xzz_xyyy[i] * tbe_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xyyz[i] = -2.0 * tr_xzz_xyyz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xyzz[i] = -2.0 * tr_xzz_xyzz[i] * tbe_0 - 2.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xzzz[i] = -2.0 * tr_xzz_xzzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_yyyy[i] = -2.0 * tr_xzz_yyyy[i] * tbe_0 - 8.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_yyyz[i] = -2.0 * tr_xzz_yyyz[i] * tbe_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_yyzz[i] = -2.0 * tr_xzz_yyzz[i] * tbe_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_yzzz[i] = -2.0 * tr_xzz_yzzz[i] * tbe_0 - 2.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_zzzz[i] = -2.0 * tr_xzz_zzzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 285-300 components of targeted buffer : FG

    auto tr_x_0_y_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 285);

    auto tr_x_0_y_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 286);

    auto tr_x_0_y_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 287);

    auto tr_x_0_y_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 288);

    auto tr_x_0_y_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 289);

    auto tr_x_0_y_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 290);

    auto tr_x_0_y_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 291);

    auto tr_x_0_y_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 292);

    auto tr_x_0_y_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 293);

    auto tr_x_0_y_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 294);

    auto tr_x_0_y_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 295);

    auto tr_x_0_y_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 296);

    auto tr_x_0_y_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 297);

    auto tr_x_0_y_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 298);

    auto tr_x_0_y_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 299);

    #pragma omp simd aligned(tr_x_0_y_zzz_xxxx, tr_x_0_y_zzz_xxxy, tr_x_0_y_zzz_xxxz, tr_x_0_y_zzz_xxyy, tr_x_0_y_zzz_xxyz, tr_x_0_y_zzz_xxzz, tr_x_0_y_zzz_xyyy, tr_x_0_y_zzz_xyyz, tr_x_0_y_zzz_xyzz, tr_x_0_y_zzz_xzzz, tr_x_0_y_zzz_yyyy, tr_x_0_y_zzz_yyyz, tr_x_0_y_zzz_yyzz, tr_x_0_y_zzz_yzzz, tr_x_0_y_zzz_zzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxy, tr_xzzz_xxxyy, tr_xzzz_xxxyz, tr_xzzz_xxy, tr_xzzz_xxyyy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyyyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyyyy, tr_xzzz_yyyyz, tr_xzzz_yyyzz, tr_xzzz_yyz, tr_xzzz_yyzzz, tr_xzzz_yzz, tr_xzzz_yzzzz, tr_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zzz_xxxx[i] = 4.0 * tr_xzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xxxy[i] = -2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xxxz[i] = 4.0 * tr_xzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xxyy[i] = -4.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xxyz[i] = -2.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xxzz[i] = 4.0 * tr_xzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xyyy[i] = -6.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xyyz[i] = -4.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xyzz[i] = -2.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xzzz[i] = 4.0 * tr_xzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_yyyy[i] = -8.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_yyyz[i] = -6.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_yyzz[i] = -4.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_yzzz[i] = -2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_zzzz[i] = 4.0 * tr_xzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-315 components of targeted buffer : FG

    auto tr_x_0_z_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 300);

    auto tr_x_0_z_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 301);

    auto tr_x_0_z_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 302);

    auto tr_x_0_z_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 303);

    auto tr_x_0_z_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 304);

    auto tr_x_0_z_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 305);

    auto tr_x_0_z_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 306);

    auto tr_x_0_z_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 307);

    auto tr_x_0_z_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 308);

    auto tr_x_0_z_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 309);

    auto tr_x_0_z_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 310);

    auto tr_x_0_z_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 311);

    auto tr_x_0_z_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 312);

    auto tr_x_0_z_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 313);

    auto tr_x_0_z_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 314);

    #pragma omp simd aligned(tr_x_0_z_xxx_xxxx, tr_x_0_z_xxx_xxxy, tr_x_0_z_xxx_xxxz, tr_x_0_z_xxx_xxyy, tr_x_0_z_xxx_xxyz, tr_x_0_z_xxx_xxzz, tr_x_0_z_xxx_xyyy, tr_x_0_z_xxx_xyyz, tr_x_0_z_xxx_xyzz, tr_x_0_z_xxx_xzzz, tr_x_0_z_xxx_yyyy, tr_x_0_z_xxx_yyyz, tr_x_0_z_xxx_yyzz, tr_x_0_z_xxx_yzzz, tr_x_0_z_xxx_zzzz, tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xx_zzzzz, tr_xxxx_xxx, tr_xxxx_xxxxz, tr_xxxx_xxxyz, tr_xxxx_xxxzz, tr_xxxx_xxy, tr_xxxx_xxyyz, tr_xxxx_xxyzz, tr_xxxx_xxz, tr_xxxx_xxzzz, tr_xxxx_xyy, tr_xxxx_xyyyz, tr_xxxx_xyyzz, tr_xxxx_xyz, tr_xxxx_xyzzz, tr_xxxx_xzz, tr_xxxx_xzzzz, tr_xxxx_yyy, tr_xxxx_yyyyz, tr_xxxx_yyyzz, tr_xxxx_yyz, tr_xxxx_yyzzz, tr_xxxx_yzz, tr_xxxx_yzzzz, tr_xxxx_zzz, tr_xxxx_zzzzz, tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xzzz, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yzzz, tr_xxxxz_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxx_xxxx[i] = -6.0 * tr_xx_xxxxz[i] * tke_0 - 6.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxxx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xxxy[i] = -6.0 * tr_xx_xxxyz[i] * tke_0 - 6.0 * tr_xxz_xxxy[i] * tbe_0 + 4.0 * tr_xxxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xxxz[i] = 3.0 * tr_xx_xxx[i] - 6.0 * tr_xx_xxxzz[i] * tke_0 - 6.0 * tr_xxz_xxxz[i] * tbe_0 - 2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xxyy[i] = -6.0 * tr_xx_xxyyz[i] * tke_0 - 6.0 * tr_xxz_xxyy[i] * tbe_0 + 4.0 * tr_xxxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xxyz[i] = 3.0 * tr_xx_xxy[i] - 6.0 * tr_xx_xxyzz[i] * tke_0 - 6.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xxzz[i] = 6.0 * tr_xx_xxz[i] - 6.0 * tr_xx_xxzzz[i] * tke_0 - 6.0 * tr_xxz_xxzz[i] * tbe_0 - 4.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xyyy[i] = -6.0 * tr_xx_xyyyz[i] * tke_0 - 6.0 * tr_xxz_xyyy[i] * tbe_0 + 4.0 * tr_xxxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xyyz[i] = 3.0 * tr_xx_xyy[i] - 6.0 * tr_xx_xyyzz[i] * tke_0 - 6.0 * tr_xxz_xyyz[i] * tbe_0 - 2.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xyzz[i] = 6.0 * tr_xx_xyz[i] - 6.0 * tr_xx_xyzzz[i] * tke_0 - 6.0 * tr_xxz_xyzz[i] * tbe_0 - 4.0 * tr_xxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xzzz[i] = 9.0 * tr_xx_xzz[i] - 6.0 * tr_xx_xzzzz[i] * tke_0 - 6.0 * tr_xxz_xzzz[i] * tbe_0 - 6.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_yyyy[i] = -6.0 * tr_xx_yyyyz[i] * tke_0 - 6.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxxx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_yyyz[i] = 3.0 * tr_xx_yyy[i] - 6.0 * tr_xx_yyyzz[i] * tke_0 - 6.0 * tr_xxz_yyyz[i] * tbe_0 - 2.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_yyzz[i] = 6.0 * tr_xx_yyz[i] - 6.0 * tr_xx_yyzzz[i] * tke_0 - 6.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_yzzz[i] = 9.0 * tr_xx_yzz[i] - 6.0 * tr_xx_yzzzz[i] * tke_0 - 6.0 * tr_xxz_yzzz[i] * tbe_0 - 6.0 * tr_xxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_zzzz[i] = 12.0 * tr_xx_zzz[i] - 6.0 * tr_xx_zzzzz[i] * tke_0 - 6.0 * tr_xxz_zzzz[i] * tbe_0 - 8.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxx_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 315-330 components of targeted buffer : FG

    auto tr_x_0_z_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 315);

    auto tr_x_0_z_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 316);

    auto tr_x_0_z_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 317);

    auto tr_x_0_z_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 318);

    auto tr_x_0_z_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 319);

    auto tr_x_0_z_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 320);

    auto tr_x_0_z_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 321);

    auto tr_x_0_z_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 322);

    auto tr_x_0_z_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 323);

    auto tr_x_0_z_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 324);

    auto tr_x_0_z_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 325);

    auto tr_x_0_z_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 326);

    auto tr_x_0_z_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 327);

    auto tr_x_0_z_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 328);

    auto tr_x_0_z_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 329);

    #pragma omp simd aligned(tr_x_0_z_xxy_xxxx, tr_x_0_z_xxy_xxxy, tr_x_0_z_xxy_xxxz, tr_x_0_z_xxy_xxyy, tr_x_0_z_xxy_xxyz, tr_x_0_z_xxy_xxzz, tr_x_0_z_xxy_xyyy, tr_x_0_z_xxy_xyyz, tr_x_0_z_xxy_xyzz, tr_x_0_z_xxy_xzzz, tr_x_0_z_xxy_yyyy, tr_x_0_z_xxy_yyyz, tr_x_0_z_xxy_yyzz, tr_x_0_z_xxy_yzzz, tr_x_0_z_xxy_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxz, tr_xxxy_xxxyz, tr_xxxy_xxxzz, tr_xxxy_xxy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xxzzz, tr_xxxy_xyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_xzzzz, tr_xxxy_yyy, tr_xxxy_yyyyz, tr_xxxy_yyyzz, tr_xxxy_yyz, tr_xxxy_yyzzz, tr_xxxy_yzz, tr_xxxy_yzzzz, tr_xxxy_zzz, tr_xxxy_zzzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxy_xxxx[i] = -4.0 * tr_xy_xxxxz[i] * tke_0 - 4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xxxy[i] = -4.0 * tr_xy_xxxyz[i] * tke_0 - 4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xxxz[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xy_xxxzz[i] * tke_0 - 4.0 * tr_xyz_xxxz[i] * tbe_0 - 2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xxyy[i] = -4.0 * tr_xy_xxyyz[i] * tke_0 - 4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xxyz[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xy_xxyzz[i] * tke_0 - 4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xxzz[i] = 4.0 * tr_xy_xxz[i] - 4.0 * tr_xy_xxzzz[i] * tke_0 - 4.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xyyy[i] = -4.0 * tr_xy_xyyyz[i] * tke_0 - 4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xyyz[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xy_xyyzz[i] * tke_0 - 4.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xyzz[i] = 4.0 * tr_xy_xyz[i] - 4.0 * tr_xy_xyzzz[i] * tke_0 - 4.0 * tr_xyz_xyzz[i] * tbe_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xzzz[i] = 6.0 * tr_xy_xzz[i] - 4.0 * tr_xy_xzzzz[i] * tke_0 - 4.0 * tr_xyz_xzzz[i] * tbe_0 - 6.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_yyyy[i] = -4.0 * tr_xy_yyyyz[i] * tke_0 - 4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_yyyz[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xy_yyyzz[i] * tke_0 - 4.0 * tr_xyz_yyyz[i] * tbe_0 - 2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_yyzz[i] = 4.0 * tr_xy_yyz[i] - 4.0 * tr_xy_yyzzz[i] * tke_0 - 4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_yzzz[i] = 6.0 * tr_xy_yzz[i] - 4.0 * tr_xy_yzzzz[i] * tke_0 - 4.0 * tr_xyz_yzzz[i] * tbe_0 - 6.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_zzzz[i] = 8.0 * tr_xy_zzz[i] - 4.0 * tr_xy_zzzzz[i] * tke_0 - 4.0 * tr_xyz_zzzz[i] * tbe_0 - 8.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-345 components of targeted buffer : FG

    auto tr_x_0_z_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 330);

    auto tr_x_0_z_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 331);

    auto tr_x_0_z_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 332);

    auto tr_x_0_z_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 333);

    auto tr_x_0_z_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 334);

    auto tr_x_0_z_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 335);

    auto tr_x_0_z_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 336);

    auto tr_x_0_z_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 337);

    auto tr_x_0_z_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 338);

    auto tr_x_0_z_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 339);

    auto tr_x_0_z_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 340);

    auto tr_x_0_z_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 341);

    auto tr_x_0_z_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 342);

    auto tr_x_0_z_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 343);

    auto tr_x_0_z_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 344);

    #pragma omp simd aligned(tr_x_0_z_xxz_xxxx, tr_x_0_z_xxz_xxxy, tr_x_0_z_xxz_xxxz, tr_x_0_z_xxz_xxyy, tr_x_0_z_xxz_xxyz, tr_x_0_z_xxz_xxzz, tr_x_0_z_xxz_xyyy, tr_x_0_z_xxz_xyyz, tr_x_0_z_xxz_xyzz, tr_x_0_z_xxz_xzzz, tr_x_0_z_xxz_yyyy, tr_x_0_z_xxz_yyyz, tr_x_0_z_xxz_yyzz, tr_x_0_z_xxz_yzzz, tr_x_0_z_xxz_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxz, tr_xxxz_xxxyz, tr_xxxz_xxxzz, tr_xxxz_xxy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xxzzz, tr_xxxz_xyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_xzzzz, tr_xxxz_yyy, tr_xxxz_yyyyz, tr_xxxz_yyyzz, tr_xxxz_yyz, tr_xxxz_yyzzz, tr_xxxz_yzz, tr_xxxz_yzzzz, tr_xxxz_zzz, tr_xxxz_zzzzz, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xzzz, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yzzz, tr_xxxzz_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxz_xxxx[i] = 2.0 * tr_x_xxxx[i] - 4.0 * tr_xz_xxxxz[i] * tke_0 - 4.0 * tr_xzz_xxxx[i] * tbe_0 - 2.0 * tr_xxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xxxy[i] = 2.0 * tr_x_xxxy[i] - 4.0 * tr_xz_xxxyz[i] * tke_0 - 4.0 * tr_xzz_xxxy[i] * tbe_0 - 2.0 * tr_xxx_xxxy[i] * tbe_0 + 4.0 * tr_xxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xxxz[i] = 2.0 * tr_x_xxxz[i] + 2.0 * tr_xz_xxx[i] - 4.0 * tr_xz_xxxzz[i] * tke_0 - 4.0 * tr_xzz_xxxz[i] * tbe_0 - 2.0 * tr_xxx_xxxz[i] * tbe_0 - 2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xxyy[i] = 2.0 * tr_x_xxyy[i] - 4.0 * tr_xz_xxyyz[i] * tke_0 - 4.0 * tr_xzz_xxyy[i] * tbe_0 - 2.0 * tr_xxx_xxyy[i] * tbe_0 + 4.0 * tr_xxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xxyz[i] = 2.0 * tr_x_xxyz[i] + 2.0 * tr_xz_xxy[i] - 4.0 * tr_xz_xxyzz[i] * tke_0 - 4.0 * tr_xzz_xxyz[i] * tbe_0 - 2.0 * tr_xxx_xxyz[i] * tbe_0 - 2.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xxzz[i] = 2.0 * tr_x_xxzz[i] + 4.0 * tr_xz_xxz[i] - 4.0 * tr_xz_xxzzz[i] * tke_0 - 4.0 * tr_xzz_xxzz[i] * tbe_0 - 2.0 * tr_xxx_xxzz[i] * tbe_0 - 4.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xyyy[i] = 2.0 * tr_x_xyyy[i] - 4.0 * tr_xz_xyyyz[i] * tke_0 - 4.0 * tr_xzz_xyyy[i] * tbe_0 - 2.0 * tr_xxx_xyyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xyyz[i] = 2.0 * tr_x_xyyz[i] + 2.0 * tr_xz_xyy[i] - 4.0 * tr_xz_xyyzz[i] * tke_0 - 4.0 * tr_xzz_xyyz[i] * tbe_0 - 2.0 * tr_xxx_xyyz[i] * tbe_0 - 2.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xyzz[i] = 2.0 * tr_x_xyzz[i] + 4.0 * tr_xz_xyz[i] - 4.0 * tr_xz_xyzzz[i] * tke_0 - 4.0 * tr_xzz_xyzz[i] * tbe_0 - 2.0 * tr_xxx_xyzz[i] * tbe_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xzzz[i] = 2.0 * tr_x_xzzz[i] + 6.0 * tr_xz_xzz[i] - 4.0 * tr_xz_xzzzz[i] * tke_0 - 4.0 * tr_xzz_xzzz[i] * tbe_0 - 2.0 * tr_xxx_xzzz[i] * tbe_0 - 6.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_yyyy[i] = 2.0 * tr_x_yyyy[i] - 4.0 * tr_xz_yyyyz[i] * tke_0 - 4.0 * tr_xzz_yyyy[i] * tbe_0 - 2.0 * tr_xxx_yyyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_yyyz[i] = 2.0 * tr_x_yyyz[i] + 2.0 * tr_xz_yyy[i] - 4.0 * tr_xz_yyyzz[i] * tke_0 - 4.0 * tr_xzz_yyyz[i] * tbe_0 - 2.0 * tr_xxx_yyyz[i] * tbe_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_yyzz[i] = 2.0 * tr_x_yyzz[i] + 4.0 * tr_xz_yyz[i] - 4.0 * tr_xz_yyzzz[i] * tke_0 - 4.0 * tr_xzz_yyzz[i] * tbe_0 - 2.0 * tr_xxx_yyzz[i] * tbe_0 - 4.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_yzzz[i] = 2.0 * tr_x_yzzz[i] + 6.0 * tr_xz_yzz[i] - 4.0 * tr_xz_yzzzz[i] * tke_0 - 4.0 * tr_xzz_yzzz[i] * tbe_0 - 2.0 * tr_xxx_yzzz[i] * tbe_0 - 6.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_zzzz[i] = 2.0 * tr_x_zzzz[i] + 8.0 * tr_xz_zzz[i] - 4.0 * tr_xz_zzzzz[i] * tke_0 - 4.0 * tr_xzz_zzzz[i] * tbe_0 - 2.0 * tr_xxx_zzzz[i] * tbe_0 - 8.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 345-360 components of targeted buffer : FG

    auto tr_x_0_z_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 345);

    auto tr_x_0_z_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 346);

    auto tr_x_0_z_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 347);

    auto tr_x_0_z_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 348);

    auto tr_x_0_z_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 349);

    auto tr_x_0_z_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 350);

    auto tr_x_0_z_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 351);

    auto tr_x_0_z_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 352);

    auto tr_x_0_z_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 353);

    auto tr_x_0_z_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 354);

    auto tr_x_0_z_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 355);

    auto tr_x_0_z_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 356);

    auto tr_x_0_z_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 357);

    auto tr_x_0_z_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 358);

    auto tr_x_0_z_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 359);

    #pragma omp simd aligned(tr_x_0_z_xyy_xxxx, tr_x_0_z_xyy_xxxy, tr_x_0_z_xyy_xxxz, tr_x_0_z_xyy_xxyy, tr_x_0_z_xyy_xxyz, tr_x_0_z_xyy_xxzz, tr_x_0_z_xyy_xyyy, tr_x_0_z_xyy_xyyz, tr_x_0_z_xyy_xyzz, tr_x_0_z_xyy_xzzz, tr_x_0_z_xyy_yyyy, tr_x_0_z_xyy_yyyz, tr_x_0_z_xyy_yyzz, tr_x_0_z_xyy_yzzz, tr_x_0_z_xyy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxz, tr_xxyy_xxxyz, tr_xxyy_xxxzz, tr_xxyy_xxy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xxzzz, tr_xxyy_xyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_xzzzz, tr_xxyy_yyy, tr_xxyy_yyyyz, tr_xxyy_yyyzz, tr_xxyy_yyz, tr_xxyy_yyzzz, tr_xxyy_yzz, tr_xxyy_yzzzz, tr_xxyy_zzz, tr_xxyy_zzzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yy_zzzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyy_xxxx[i] = -2.0 * tr_yy_xxxxz[i] * tke_0 - 2.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xxxy[i] = -2.0 * tr_yy_xxxyz[i] * tke_0 - 2.0 * tr_yyz_xxxy[i] * tbe_0 + 4.0 * tr_xxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xxxz[i] = tr_yy_xxx[i] - 2.0 * tr_yy_xxxzz[i] * tke_0 - 2.0 * tr_yyz_xxxz[i] * tbe_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xxyy[i] = -2.0 * tr_yy_xxyyz[i] * tke_0 - 2.0 * tr_yyz_xxyy[i] * tbe_0 + 4.0 * tr_xxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xxyz[i] = tr_yy_xxy[i] - 2.0 * tr_yy_xxyzz[i] * tke_0 - 2.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xxzz[i] = 2.0 * tr_yy_xxz[i] - 2.0 * tr_yy_xxzzz[i] * tke_0 - 2.0 * tr_yyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xyyy[i] = -2.0 * tr_yy_xyyyz[i] * tke_0 - 2.0 * tr_yyz_xyyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xyyz[i] = tr_yy_xyy[i] - 2.0 * tr_yy_xyyzz[i] * tke_0 - 2.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xyzz[i] = 2.0 * tr_yy_xyz[i] - 2.0 * tr_yy_xyzzz[i] * tke_0 - 2.0 * tr_yyz_xyzz[i] * tbe_0 - 4.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xzzz[i] = 3.0 * tr_yy_xzz[i] - 2.0 * tr_yy_xzzzz[i] * tke_0 - 2.0 * tr_yyz_xzzz[i] * tbe_0 - 6.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_yyyy[i] = -2.0 * tr_yy_yyyyz[i] * tke_0 - 2.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_yyyz[i] = tr_yy_yyy[i] - 2.0 * tr_yy_yyyzz[i] * tke_0 - 2.0 * tr_yyz_yyyz[i] * tbe_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_yyzz[i] = 2.0 * tr_yy_yyz[i] - 2.0 * tr_yy_yyzzz[i] * tke_0 - 2.0 * tr_yyz_yyzz[i] * tbe_0 - 4.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_yzzz[i] = 3.0 * tr_yy_yzz[i] - 2.0 * tr_yy_yzzzz[i] * tke_0 - 2.0 * tr_yyz_yzzz[i] * tbe_0 - 6.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_zzzz[i] = 4.0 * tr_yy_zzz[i] - 2.0 * tr_yy_zzzzz[i] * tke_0 - 2.0 * tr_yyz_zzzz[i] * tbe_0 - 8.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-375 components of targeted buffer : FG

    auto tr_x_0_z_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 360);

    auto tr_x_0_z_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 361);

    auto tr_x_0_z_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 362);

    auto tr_x_0_z_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 363);

    auto tr_x_0_z_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 364);

    auto tr_x_0_z_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 365);

    auto tr_x_0_z_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 366);

    auto tr_x_0_z_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 367);

    auto tr_x_0_z_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 368);

    auto tr_x_0_z_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 369);

    auto tr_x_0_z_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 370);

    auto tr_x_0_z_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 371);

    auto tr_x_0_z_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 372);

    auto tr_x_0_z_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 373);

    auto tr_x_0_z_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 374);

    #pragma omp simd aligned(tr_x_0_z_xyz_xxxx, tr_x_0_z_xyz_xxxy, tr_x_0_z_xyz_xxxz, tr_x_0_z_xyz_xxyy, tr_x_0_z_xyz_xxyz, tr_x_0_z_xyz_xxzz, tr_x_0_z_xyz_xyyy, tr_x_0_z_xyz_xyyz, tr_x_0_z_xyz_xyzz, tr_x_0_z_xyz_xzzz, tr_x_0_z_xyz_yyyy, tr_x_0_z_xyz_yyyz, tr_x_0_z_xyz_yyzz, tr_x_0_z_xyz_yzzz, tr_x_0_z_xyz_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxz, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxyz_zzzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyz_xxxx[i] = tr_y_xxxx[i] - 2.0 * tr_yz_xxxxz[i] * tke_0 - 2.0 * tr_yzz_xxxx[i] * tbe_0 - 2.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xxxy[i] = tr_y_xxxy[i] - 2.0 * tr_yz_xxxyz[i] * tke_0 - 2.0 * tr_yzz_xxxy[i] * tbe_0 - 2.0 * tr_xxy_xxxy[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xxxz[i] = tr_y_xxxz[i] + tr_yz_xxx[i] - 2.0 * tr_yz_xxxzz[i] * tke_0 - 2.0 * tr_yzz_xxxz[i] * tbe_0 - 2.0 * tr_xxy_xxxz[i] * tbe_0 - 2.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xxyy[i] = tr_y_xxyy[i] - 2.0 * tr_yz_xxyyz[i] * tke_0 - 2.0 * tr_yzz_xxyy[i] * tbe_0 - 2.0 * tr_xxy_xxyy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xxyz[i] = tr_y_xxyz[i] + tr_yz_xxy[i] - 2.0 * tr_yz_xxyzz[i] * tke_0 - 2.0 * tr_yzz_xxyz[i] * tbe_0 - 2.0 * tr_xxy_xxyz[i] * tbe_0 - 2.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xxzz[i] = tr_y_xxzz[i] + 2.0 * tr_yz_xxz[i] - 2.0 * tr_yz_xxzzz[i] * tke_0 - 2.0 * tr_yzz_xxzz[i] * tbe_0 - 2.0 * tr_xxy_xxzz[i] * tbe_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xyyy[i] = tr_y_xyyy[i] - 2.0 * tr_yz_xyyyz[i] * tke_0 - 2.0 * tr_yzz_xyyy[i] * tbe_0 - 2.0 * tr_xxy_xyyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xyyz[i] = tr_y_xyyz[i] + tr_yz_xyy[i] - 2.0 * tr_yz_xyyzz[i] * tke_0 - 2.0 * tr_yzz_xyyz[i] * tbe_0 - 2.0 * tr_xxy_xyyz[i] * tbe_0 - 2.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xyzz[i] = tr_y_xyzz[i] + 2.0 * tr_yz_xyz[i] - 2.0 * tr_yz_xyzzz[i] * tke_0 - 2.0 * tr_yzz_xyzz[i] * tbe_0 - 2.0 * tr_xxy_xyzz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xzzz[i] = tr_y_xzzz[i] + 3.0 * tr_yz_xzz[i] - 2.0 * tr_yz_xzzzz[i] * tke_0 - 2.0 * tr_yzz_xzzz[i] * tbe_0 - 2.0 * tr_xxy_xzzz[i] * tbe_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_yyyy[i] = tr_y_yyyy[i] - 2.0 * tr_yz_yyyyz[i] * tke_0 - 2.0 * tr_yzz_yyyy[i] * tbe_0 - 2.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_yyyz[i] = tr_y_yyyz[i] + tr_yz_yyy[i] - 2.0 * tr_yz_yyyzz[i] * tke_0 - 2.0 * tr_yzz_yyyz[i] * tbe_0 - 2.0 * tr_xxy_yyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_yyzz[i] = tr_y_yyzz[i] + 2.0 * tr_yz_yyz[i] - 2.0 * tr_yz_yyzzz[i] * tke_0 - 2.0 * tr_yzz_yyzz[i] * tbe_0 - 2.0 * tr_xxy_yyzz[i] * tbe_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_yzzz[i] = tr_y_yzzz[i] + 3.0 * tr_yz_yzz[i] - 2.0 * tr_yz_yzzzz[i] * tke_0 - 2.0 * tr_yzz_yzzz[i] * tbe_0 - 2.0 * tr_xxy_yzzz[i] * tbe_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_zzzz[i] = tr_y_zzzz[i] + 4.0 * tr_yz_zzz[i] - 2.0 * tr_yz_zzzzz[i] * tke_0 - 2.0 * tr_yzz_zzzz[i] * tbe_0 - 2.0 * tr_xxy_zzzz[i] * tbe_0 - 8.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 375-390 components of targeted buffer : FG

    auto tr_x_0_z_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 375);

    auto tr_x_0_z_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 376);

    auto tr_x_0_z_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 377);

    auto tr_x_0_z_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 378);

    auto tr_x_0_z_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 379);

    auto tr_x_0_z_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 380);

    auto tr_x_0_z_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 381);

    auto tr_x_0_z_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 382);

    auto tr_x_0_z_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 383);

    auto tr_x_0_z_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 384);

    auto tr_x_0_z_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 385);

    auto tr_x_0_z_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 386);

    auto tr_x_0_z_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 387);

    auto tr_x_0_z_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 388);

    auto tr_x_0_z_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 389);

    #pragma omp simd aligned(tr_x_0_z_xzz_xxxx, tr_x_0_z_xzz_xxxy, tr_x_0_z_xzz_xxxz, tr_x_0_z_xzz_xxyy, tr_x_0_z_xzz_xxyz, tr_x_0_z_xzz_xxzz, tr_x_0_z_xzz_xyyy, tr_x_0_z_xzz_xyyz, tr_x_0_z_xzz_xyzz, tr_x_0_z_xzz_xzzz, tr_x_0_z_xzz_yyyy, tr_x_0_z_xzz_yyyz, tr_x_0_z_xzz_yyzz, tr_x_0_z_xzz_yzzz, tr_x_0_z_xzz_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxz, tr_xxzz_xxxyz, tr_xxzz_xxxzz, tr_xxzz_xxy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xxzzz, tr_xxzz_xyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_xzzzz, tr_xxzz_yyy, tr_xxzz_yyyyz, tr_xxzz_yyyzz, tr_xxzz_yyz, tr_xxzz_yyzzz, tr_xxzz_yzz, tr_xxzz_yzzzz, tr_xxzz_zzz, tr_xxzz_zzzzz, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xzzz, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yzzz, tr_xxzzz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xzz_xxxx[i] = 2.0 * tr_z_xxxx[i] - 2.0 * tr_zz_xxxxz[i] * tke_0 - 2.0 * tr_zzz_xxxx[i] * tbe_0 - 4.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xxxy[i] = 2.0 * tr_z_xxxy[i] - 2.0 * tr_zz_xxxyz[i] * tke_0 - 2.0 * tr_zzz_xxxy[i] * tbe_0 - 4.0 * tr_xxz_xxxy[i] * tbe_0 + 4.0 * tr_xxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xxxz[i] = 2.0 * tr_z_xxxz[i] + tr_zz_xxx[i] - 2.0 * tr_zz_xxxzz[i] * tke_0 - 2.0 * tr_zzz_xxxz[i] * tbe_0 - 4.0 * tr_xxz_xxxz[i] * tbe_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xxyy[i] = 2.0 * tr_z_xxyy[i] - 2.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_zzz_xxyy[i] * tbe_0 - 4.0 * tr_xxz_xxyy[i] * tbe_0 + 4.0 * tr_xxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xxyz[i] = 2.0 * tr_z_xxyz[i] + tr_zz_xxy[i] - 2.0 * tr_zz_xxyzz[i] * tke_0 - 2.0 * tr_zzz_xxyz[i] * tbe_0 - 4.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xxzz[i] = 2.0 * tr_z_xxzz[i] + 2.0 * tr_zz_xxz[i] - 2.0 * tr_zz_xxzzz[i] * tke_0 - 2.0 * tr_zzz_xxzz[i] * tbe_0 - 4.0 * tr_xxz_xxzz[i] * tbe_0 - 4.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xyyy[i] = 2.0 * tr_z_xyyy[i] - 2.0 * tr_zz_xyyyz[i] * tke_0 - 2.0 * tr_zzz_xyyy[i] * tbe_0 - 4.0 * tr_xxz_xyyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xyyz[i] = 2.0 * tr_z_xyyz[i] + tr_zz_xyy[i] - 2.0 * tr_zz_xyyzz[i] * tke_0 - 2.0 * tr_zzz_xyyz[i] * tbe_0 - 4.0 * tr_xxz_xyyz[i] * tbe_0 - 2.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xyzz[i] = 2.0 * tr_z_xyzz[i] + 2.0 * tr_zz_xyz[i] - 2.0 * tr_zz_xyzzz[i] * tke_0 - 2.0 * tr_zzz_xyzz[i] * tbe_0 - 4.0 * tr_xxz_xyzz[i] * tbe_0 - 4.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xzzz[i] = 2.0 * tr_z_xzzz[i] + 3.0 * tr_zz_xzz[i] - 2.0 * tr_zz_xzzzz[i] * tke_0 - 2.0 * tr_zzz_xzzz[i] * tbe_0 - 4.0 * tr_xxz_xzzz[i] * tbe_0 - 6.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_yyyy[i] = 2.0 * tr_z_yyyy[i] - 2.0 * tr_zz_yyyyz[i] * tke_0 - 2.0 * tr_zzz_yyyy[i] * tbe_0 - 4.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_yyyz[i] = 2.0 * tr_z_yyyz[i] + tr_zz_yyy[i] - 2.0 * tr_zz_yyyzz[i] * tke_0 - 2.0 * tr_zzz_yyyz[i] * tbe_0 - 4.0 * tr_xxz_yyyz[i] * tbe_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_yyzz[i] = 2.0 * tr_z_yyzz[i] + 2.0 * tr_zz_yyz[i] - 2.0 * tr_zz_yyzzz[i] * tke_0 - 2.0 * tr_zzz_yyzz[i] * tbe_0 - 4.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_yzzz[i] = 2.0 * tr_z_yzzz[i] + 3.0 * tr_zz_yzz[i] - 2.0 * tr_zz_yzzzz[i] * tke_0 - 2.0 * tr_zzz_yzzz[i] * tbe_0 - 4.0 * tr_xxz_yzzz[i] * tbe_0 - 6.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_zzzz[i] = 2.0 * tr_z_zzzz[i] + 4.0 * tr_zz_zzz[i] - 2.0 * tr_zz_zzzzz[i] * tke_0 - 2.0 * tr_zzz_zzzz[i] * tbe_0 - 4.0 * tr_xxz_zzzz[i] * tbe_0 - 8.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-405 components of targeted buffer : FG

    auto tr_x_0_z_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 390);

    auto tr_x_0_z_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 391);

    auto tr_x_0_z_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 392);

    auto tr_x_0_z_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 393);

    auto tr_x_0_z_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 394);

    auto tr_x_0_z_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 395);

    auto tr_x_0_z_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 396);

    auto tr_x_0_z_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 397);

    auto tr_x_0_z_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 398);

    auto tr_x_0_z_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 399);

    auto tr_x_0_z_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 400);

    auto tr_x_0_z_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 401);

    auto tr_x_0_z_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 402);

    auto tr_x_0_z_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 403);

    auto tr_x_0_z_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 404);

    #pragma omp simd aligned(tr_x_0_z_yyy_xxxx, tr_x_0_z_yyy_xxxy, tr_x_0_z_yyy_xxxz, tr_x_0_z_yyy_xxyy, tr_x_0_z_yyy_xxyz, tr_x_0_z_yyy_xxzz, tr_x_0_z_yyy_xyyy, tr_x_0_z_yyy_xyyz, tr_x_0_z_yyy_xyzz, tr_x_0_z_yyy_xzzz, tr_x_0_z_yyy_yyyy, tr_x_0_z_yyy_yyyz, tr_x_0_z_yyy_yyzz, tr_x_0_z_yyy_yzzz, tr_x_0_z_yyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxz, tr_xyyy_xxxyz, tr_xyyy_xxxzz, tr_xyyy_xxy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xxzzz, tr_xyyy_xyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_xzzzz, tr_xyyy_yyy, tr_xyyy_yyyyz, tr_xyyy_yyyzz, tr_xyyy_yyz, tr_xyyy_yyzzz, tr_xyyy_yzz, tr_xyyy_yzzzz, tr_xyyy_zzz, tr_xyyy_zzzzz, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyy_xxxx[i] = 4.0 * tr_xyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xxxy[i] = 4.0 * tr_xyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xxxz[i] = -2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xxyy[i] = 4.0 * tr_xyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xxyz[i] = -2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xxzz[i] = -4.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xyyy[i] = 4.0 * tr_xyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xyyz[i] = -2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xyzz[i] = -4.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xzzz[i] = -6.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_yyyy[i] = 4.0 * tr_xyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_yyyz[i] = -2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_yyzz[i] = -4.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_yzzz[i] = -6.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_zzzz[i] = -8.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 405-420 components of targeted buffer : FG

    auto tr_x_0_z_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 405);

    auto tr_x_0_z_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 406);

    auto tr_x_0_z_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 407);

    auto tr_x_0_z_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 408);

    auto tr_x_0_z_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 409);

    auto tr_x_0_z_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 410);

    auto tr_x_0_z_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 411);

    auto tr_x_0_z_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 412);

    auto tr_x_0_z_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 413);

    auto tr_x_0_z_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 414);

    auto tr_x_0_z_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 415);

    auto tr_x_0_z_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 416);

    auto tr_x_0_z_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 417);

    auto tr_x_0_z_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 418);

    auto tr_x_0_z_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 419);

    #pragma omp simd aligned(tr_x_0_z_yyz_xxxx, tr_x_0_z_yyz_xxxy, tr_x_0_z_yyz_xxxz, tr_x_0_z_yyz_xxyy, tr_x_0_z_yyz_xxyz, tr_x_0_z_yyz_xxzz, tr_x_0_z_yyz_xyyy, tr_x_0_z_yyz_xyyz, tr_x_0_z_yyz_xyzz, tr_x_0_z_yyz_xzzz, tr_x_0_z_yyz_yyyy, tr_x_0_z_yyz_yyyz, tr_x_0_z_yyz_yyzz, tr_x_0_z_yyz_yzzz, tr_x_0_z_yyz_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxz, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyyz_zzzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyz_xxxx[i] = -2.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xxxy[i] = -2.0 * tr_xyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xxxz[i] = -2.0 * tr_xyy_xxxz[i] * tbe_0 - 2.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xxyy[i] = -2.0 * tr_xyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xxyz[i] = -2.0 * tr_xyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xxzz[i] = -2.0 * tr_xyy_xxzz[i] * tbe_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xyyy[i] = -2.0 * tr_xyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xyyz[i] = -2.0 * tr_xyy_xyyz[i] * tbe_0 - 2.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xyzz[i] = -2.0 * tr_xyy_xyzz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xzzz[i] = -2.0 * tr_xyy_xzzz[i] * tbe_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_yyyy[i] = -2.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_yyyz[i] = -2.0 * tr_xyy_yyyz[i] * tbe_0 - 2.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_yyzz[i] = -2.0 * tr_xyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_yzzz[i] = -2.0 * tr_xyy_yzzz[i] * tbe_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_zzzz[i] = -2.0 * tr_xyy_zzzz[i] * tbe_0 - 8.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-435 components of targeted buffer : FG

    auto tr_x_0_z_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 420);

    auto tr_x_0_z_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 421);

    auto tr_x_0_z_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 422);

    auto tr_x_0_z_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 423);

    auto tr_x_0_z_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 424);

    auto tr_x_0_z_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 425);

    auto tr_x_0_z_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 426);

    auto tr_x_0_z_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 427);

    auto tr_x_0_z_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 428);

    auto tr_x_0_z_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 429);

    auto tr_x_0_z_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 430);

    auto tr_x_0_z_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 431);

    auto tr_x_0_z_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 432);

    auto tr_x_0_z_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 433);

    auto tr_x_0_z_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 434);

    #pragma omp simd aligned(tr_x_0_z_yzz_xxxx, tr_x_0_z_yzz_xxxy, tr_x_0_z_yzz_xxxz, tr_x_0_z_yzz_xxyy, tr_x_0_z_yzz_xxyz, tr_x_0_z_yzz_xxzz, tr_x_0_z_yzz_xyyy, tr_x_0_z_yzz_xyyz, tr_x_0_z_yzz_xyzz, tr_x_0_z_yzz_xzzz, tr_x_0_z_yzz_yyyy, tr_x_0_z_yzz_yyyz, tr_x_0_z_yzz_yyzz, tr_x_0_z_yzz_yzzz, tr_x_0_z_yzz_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxz, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xyzz_zzzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yzz_xxxx[i] = -4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xxxy[i] = -4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xxxz[i] = -4.0 * tr_xyz_xxxz[i] * tbe_0 - 2.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xxyy[i] = -4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xxyz[i] = -4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xxzz[i] = -4.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xyyy[i] = -4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xyyz[i] = -4.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xyzz[i] = -4.0 * tr_xyz_xyzz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xzzz[i] = -4.0 * tr_xyz_xzzz[i] * tbe_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_yyyy[i] = -4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_yyyz[i] = -4.0 * tr_xyz_yyyz[i] * tbe_0 - 2.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_yyzz[i] = -4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_yzzz[i] = -4.0 * tr_xyz_yzzz[i] * tbe_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_zzzz[i] = -4.0 * tr_xyz_zzzz[i] * tbe_0 - 8.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 435-450 components of targeted buffer : FG

    auto tr_x_0_z_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 435);

    auto tr_x_0_z_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 436);

    auto tr_x_0_z_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 437);

    auto tr_x_0_z_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 438);

    auto tr_x_0_z_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 439);

    auto tr_x_0_z_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 440);

    auto tr_x_0_z_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 441);

    auto tr_x_0_z_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 442);

    auto tr_x_0_z_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 443);

    auto tr_x_0_z_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 444);

    auto tr_x_0_z_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 445);

    auto tr_x_0_z_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 446);

    auto tr_x_0_z_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 447);

    auto tr_x_0_z_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 448);

    auto tr_x_0_z_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 449);

    #pragma omp simd aligned(tr_x_0_z_zzz_xxxx, tr_x_0_z_zzz_xxxy, tr_x_0_z_zzz_xxxz, tr_x_0_z_zzz_xxyy, tr_x_0_z_zzz_xxyz, tr_x_0_z_zzz_xxzz, tr_x_0_z_zzz_xyyy, tr_x_0_z_zzz_xyyz, tr_x_0_z_zzz_xyzz, tr_x_0_z_zzz_xzzz, tr_x_0_z_zzz_yyyy, tr_x_0_z_zzz_yyyz, tr_x_0_z_zzz_yyzz, tr_x_0_z_zzz_yzzz, tr_x_0_z_zzz_zzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxz, tr_xzzz_xxxyz, tr_xzzz_xxxzz, tr_xzzz_xxy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xxzzz, tr_xzzz_xyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_xzzzz, tr_xzzz_yyy, tr_xzzz_yyyyz, tr_xzzz_yyyzz, tr_xzzz_yyz, tr_xzzz_yyzzz, tr_xzzz_yzz, tr_xzzz_yzzzz, tr_xzzz_zzz, tr_xzzz_zzzzz, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xzzz, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yzzz, tr_xzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zzz_xxxx[i] = -6.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xxxy[i] = -6.0 * tr_xzz_xxxy[i] * tbe_0 + 4.0 * tr_xzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xxxz[i] = -6.0 * tr_xzz_xxxz[i] * tbe_0 - 2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xxyy[i] = -6.0 * tr_xzz_xxyy[i] * tbe_0 + 4.0 * tr_xzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xxyz[i] = -6.0 * tr_xzz_xxyz[i] * tbe_0 - 2.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xxzz[i] = -6.0 * tr_xzz_xxzz[i] * tbe_0 - 4.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xyyy[i] = -6.0 * tr_xzz_xyyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xyyz[i] = -6.0 * tr_xzz_xyyz[i] * tbe_0 - 2.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xyzz[i] = -6.0 * tr_xzz_xyzz[i] * tbe_0 - 4.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xzzz[i] = -6.0 * tr_xzz_xzzz[i] * tbe_0 - 6.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_yyyy[i] = -6.0 * tr_xzz_yyyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_yyyz[i] = -6.0 * tr_xzz_yyyz[i] * tbe_0 - 2.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_yyzz[i] = -6.0 * tr_xzz_yyzz[i] * tbe_0 - 4.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_yzzz[i] = -6.0 * tr_xzz_yzzz[i] * tbe_0 - 6.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_zzzz[i] = -6.0 * tr_xzz_zzzz[i] * tbe_0 - 8.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-465 components of targeted buffer : FG

    auto tr_y_0_x_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 450);

    auto tr_y_0_x_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 451);

    auto tr_y_0_x_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 452);

    auto tr_y_0_x_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 453);

    auto tr_y_0_x_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 454);

    auto tr_y_0_x_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 455);

    auto tr_y_0_x_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 456);

    auto tr_y_0_x_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 457);

    auto tr_y_0_x_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 458);

    auto tr_y_0_x_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 459);

    auto tr_y_0_x_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 460);

    auto tr_y_0_x_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 461);

    auto tr_y_0_x_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 462);

    auto tr_y_0_x_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 463);

    auto tr_y_0_x_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 464);

    #pragma omp simd aligned(tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xzzz, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yzzz, tr_xxxxy_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxx, tr_xxxy_xxxxy, tr_xxxy_xxxxz, tr_xxxy_xxxyy, tr_xxxy_xxxyz, tr_xxxy_xxxzz, tr_xxxy_xxy, tr_xxxy_xxyyy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xxzzz, tr_xxxy_xyy, tr_xxxy_xyyyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_xzzzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_y_0_x_xxx_xxxx, tr_y_0_x_xxx_xxxy, tr_y_0_x_xxx_xxxz, tr_y_0_x_xxx_xxyy, tr_y_0_x_xxx_xxyz, tr_y_0_x_xxx_xxzz, tr_y_0_x_xxx_xyyy, tr_y_0_x_xxx_xyyz, tr_y_0_x_xxx_xyzz, tr_y_0_x_xxx_xzzz, tr_y_0_x_xxx_yyyy, tr_y_0_x_xxx_yyyz, tr_y_0_x_xxx_yyzz, tr_y_0_x_xxx_yzzz, tr_y_0_x_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxx_xxxx[i] = -6.0 * tr_xxy_xxxx[i] * tbe_0 - 8.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xxxy[i] = -6.0 * tr_xxy_xxxy[i] * tbe_0 - 6.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xxxz[i] = -6.0 * tr_xxy_xxxz[i] * tbe_0 - 6.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xxyy[i] = -6.0 * tr_xxy_xxyy[i] * tbe_0 - 4.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xxyz[i] = -6.0 * tr_xxy_xxyz[i] * tbe_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xxzz[i] = -6.0 * tr_xxy_xxzz[i] * tbe_0 - 4.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xyyy[i] = -6.0 * tr_xxy_xyyy[i] * tbe_0 - 2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xyyz[i] = -6.0 * tr_xxy_xyyz[i] * tbe_0 - 2.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xyzz[i] = -6.0 * tr_xxy_xyzz[i] * tbe_0 - 2.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xzzz[i] = -6.0 * tr_xxy_xzzz[i] * tbe_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_yyyy[i] = -6.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_yyyz[i] = -6.0 * tr_xxy_yyyz[i] * tbe_0 + 4.0 * tr_xxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_yyzz[i] = -6.0 * tr_xxy_yyzz[i] * tbe_0 + 4.0 * tr_xxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_yzzz[i] = -6.0 * tr_xxy_yzzz[i] * tbe_0 + 4.0 * tr_xxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_zzzz[i] = -6.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 465-480 components of targeted buffer : FG

    auto tr_y_0_x_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 465);

    auto tr_y_0_x_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 466);

    auto tr_y_0_x_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 467);

    auto tr_y_0_x_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 468);

    auto tr_y_0_x_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 469);

    auto tr_y_0_x_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 470);

    auto tr_y_0_x_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 471);

    auto tr_y_0_x_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 472);

    auto tr_y_0_x_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 473);

    auto tr_y_0_x_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 474);

    auto tr_y_0_x_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 475);

    auto tr_y_0_x_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 476);

    auto tr_y_0_x_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 477);

    auto tr_y_0_x_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 478);

    auto tr_y_0_x_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 479);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xzzz, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yzzz, tr_xxxyy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxx, tr_xxyy_xxxxy, tr_xxyy_xxxxz, tr_xxyy_xxxyy, tr_xxyy_xxxyz, tr_xxyy_xxxzz, tr_xxyy_xxy, tr_xxyy_xxyyy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xxzzz, tr_xxyy_xyy, tr_xxyy_xyyyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_xzzzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_y_0_x_xxy_xxxx, tr_y_0_x_xxy_xxxy, tr_y_0_x_xxy_xxxz, tr_y_0_x_xxy_xxyy, tr_y_0_x_xxy_xxyz, tr_y_0_x_xxy_xxzz, tr_y_0_x_xxy_xyyy, tr_y_0_x_xxy_xyyz, tr_y_0_x_xxy_xyzz, tr_y_0_x_xxy_xzzz, tr_y_0_x_xxy_yyyy, tr_y_0_x_xxy_yyyz, tr_y_0_x_xxy_yyzz, tr_y_0_x_xxy_yzzz, tr_y_0_x_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxy_xxxx[i] = 2.0 * tr_x_xxxx[i] - 4.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xx_xxx[i] - 2.0 * tr_xx_xxxxx[i] * tke_0 - 8.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xxxy[i] = 2.0 * tr_x_xxxy[i] - 4.0 * tr_xyy_xxxy[i] * tbe_0 + 3.0 * tr_xx_xxy[i] - 2.0 * tr_xx_xxxxy[i] * tke_0 - 6.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xxxz[i] = 2.0 * tr_x_xxxz[i] - 4.0 * tr_xyy_xxxz[i] * tbe_0 + 3.0 * tr_xx_xxz[i] - 2.0 * tr_xx_xxxxz[i] * tke_0 - 6.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xxyy[i] = 2.0 * tr_x_xxyy[i] - 4.0 * tr_xyy_xxyy[i] * tbe_0 + 2.0 * tr_xx_xyy[i] - 2.0 * tr_xx_xxxyy[i] * tke_0 - 4.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyy_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxyy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xxyz[i] = 2.0 * tr_x_xxyz[i] - 4.0 * tr_xyy_xxyz[i] * tbe_0 + 2.0 * tr_xx_xyz[i] - 2.0 * tr_xx_xxxyz[i] * tke_0 - 4.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxyz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xxzz[i] = 2.0 * tr_x_xxzz[i] - 4.0 * tr_xyy_xxzz[i] * tbe_0 + 2.0 * tr_xx_xzz[i] - 2.0 * tr_xx_xxxzz[i] * tke_0 - 4.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxzz[i] * tbe_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xyyy[i] = 2.0 * tr_x_xyyy[i] - 4.0 * tr_xyy_xyyy[i] * tbe_0 + tr_xx_yyy[i] - 2.0 * tr_xx_xxyyy[i] * tke_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyy_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xyyz[i] = 2.0 * tr_x_xyyz[i] - 4.0 * tr_xyy_xyyz[i] * tbe_0 + tr_xx_yyz[i] - 2.0 * tr_xx_xxyyz[i] * tke_0 - 2.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyy_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyyz[i] * tbe_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xyzz[i] = 2.0 * tr_x_xyzz[i] - 4.0 * tr_xyy_xyzz[i] * tbe_0 + tr_xx_yzz[i] - 2.0 * tr_xx_xxyzz[i] * tke_0 - 2.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyzz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xzzz[i] = 2.0 * tr_x_xzzz[i] - 4.0 * tr_xyy_xzzz[i] * tbe_0 + tr_xx_zzz[i] - 2.0 * tr_xx_xxzzz[i] * tke_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xzzz[i] * tbe_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_yyyy[i] = 2.0 * tr_x_yyyy[i] - 4.0 * tr_xyy_yyyy[i] * tbe_0 - 2.0 * tr_xx_xyyyy[i] * tke_0 + 4.0 * tr_xxyy_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyyy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_yyyz[i] = 2.0 * tr_x_yyyz[i] - 4.0 * tr_xyy_yyyz[i] * tbe_0 - 2.0 * tr_xx_xyyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyyz[i] * tbe_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_yyzz[i] = 2.0 * tr_x_yyzz[i] - 4.0 * tr_xyy_yyzz[i] * tbe_0 - 2.0 * tr_xx_xyyzz[i] * tke_0 + 4.0 * tr_xxyy_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyzz[i] * tbe_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_yzzz[i] = 2.0 * tr_x_yzzz[i] - 4.0 * tr_xyy_yzzz[i] * tbe_0 - 2.0 * tr_xx_xyzzz[i] * tke_0 + 4.0 * tr_xxyy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yzzz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_zzzz[i] = 2.0 * tr_x_zzzz[i] - 4.0 * tr_xyy_zzzz[i] * tbe_0 - 2.0 * tr_xx_xzzzz[i] * tke_0 + 4.0 * tr_xxyy_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zzzz[i] * tbe_0 + 4.0 * tr_xxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-495 components of targeted buffer : FG

    auto tr_y_0_x_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 480);

    auto tr_y_0_x_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 481);

    auto tr_y_0_x_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 482);

    auto tr_y_0_x_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 483);

    auto tr_y_0_x_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 484);

    auto tr_y_0_x_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 485);

    auto tr_y_0_x_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 486);

    auto tr_y_0_x_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 487);

    auto tr_y_0_x_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 488);

    auto tr_y_0_x_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 489);

    auto tr_y_0_x_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 490);

    auto tr_y_0_x_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 491);

    auto tr_y_0_x_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 492);

    auto tr_y_0_x_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 493);

    auto tr_y_0_x_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 494);

    #pragma omp simd aligned(tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxx, tr_xxyz_xxxxy, tr_xxyz_xxxxz, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_y_0_x_xxz_xxxx, tr_y_0_x_xxz_xxxy, tr_y_0_x_xxz_xxxz, tr_y_0_x_xxz_xxyy, tr_y_0_x_xxz_xxyz, tr_y_0_x_xxz_xxzz, tr_y_0_x_xxz_xyyy, tr_y_0_x_xxz_xyyz, tr_y_0_x_xxz_xyzz, tr_y_0_x_xxz_xzzz, tr_y_0_x_xxz_yyyy, tr_y_0_x_xxz_yyyz, tr_y_0_x_xxz_yyzz, tr_y_0_x_xxz_yzzz, tr_y_0_x_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxz_xxxx[i] = -4.0 * tr_xyz_xxxx[i] * tbe_0 - 8.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xxxy[i] = -4.0 * tr_xyz_xxxy[i] * tbe_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xxxz[i] = -4.0 * tr_xyz_xxxz[i] * tbe_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xxyy[i] = -4.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xxyz[i] = -4.0 * tr_xyz_xxyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xxzz[i] = -4.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xyyy[i] = -4.0 * tr_xyz_xyyy[i] * tbe_0 - 2.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xyyz[i] = -4.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xyzz[i] = -4.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xzzz[i] = -4.0 * tr_xyz_xzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_yyyy[i] = -4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_yyyz[i] = -4.0 * tr_xyz_yyyz[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_yyzz[i] = -4.0 * tr_xyz_yyzz[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_yzzz[i] = -4.0 * tr_xyz_yzzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_zzzz[i] = -4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 495-510 components of targeted buffer : FG

    auto tr_y_0_x_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 495);

    auto tr_y_0_x_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 496);

    auto tr_y_0_x_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 497);

    auto tr_y_0_x_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 498);

    auto tr_y_0_x_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 499);

    auto tr_y_0_x_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 500);

    auto tr_y_0_x_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 501);

    auto tr_y_0_x_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 502);

    auto tr_y_0_x_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 503);

    auto tr_y_0_x_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 504);

    auto tr_y_0_x_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 505);

    auto tr_y_0_x_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 506);

    auto tr_y_0_x_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 507);

    auto tr_y_0_x_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 508);

    auto tr_y_0_x_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 509);

    #pragma omp simd aligned(tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xzzz, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yzzz, tr_xxyyy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyyy_xxx, tr_xyyy_xxxxx, tr_xyyy_xxxxy, tr_xyyy_xxxxz, tr_xyyy_xxxyy, tr_xyyy_xxxyz, tr_xyyy_xxxzz, tr_xyyy_xxy, tr_xyyy_xxyyy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xxzzz, tr_xyyy_xyy, tr_xyyy_xyyyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_xzzzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_y_0_x_xyy_xxxx, tr_y_0_x_xyy_xxxy, tr_y_0_x_xyy_xxxz, tr_y_0_x_xyy_xxyy, tr_y_0_x_xyy_xxyz, tr_y_0_x_xyy_xxzz, tr_y_0_x_xyy_xyyy, tr_y_0_x_xyy_xyyz, tr_y_0_x_xyy_xyzz, tr_y_0_x_xyy_xzzz, tr_y_0_x_xyy_yyyy, tr_y_0_x_xyy_yyyz, tr_y_0_x_xyy_yyzz, tr_y_0_x_xyy_yzzz, tr_y_0_x_xyy_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyy_xxxx[i] = 2.0 * tr_y_xxxx[i] - 2.0 * tr_yyy_xxxx[i] * tbe_0 + 8.0 * tr_xy_xxx[i] - 4.0 * tr_xy_xxxxx[i] * tke_0 - 8.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xxxy[i] = 2.0 * tr_y_xxxy[i] - 2.0 * tr_yyy_xxxy[i] * tbe_0 + 6.0 * tr_xy_xxy[i] - 4.0 * tr_xy_xxxxy[i] * tke_0 - 6.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyy_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xxxz[i] = 2.0 * tr_y_xxxz[i] - 2.0 * tr_yyy_xxxz[i] * tbe_0 + 6.0 * tr_xy_xxz[i] - 4.0 * tr_xy_xxxxz[i] * tke_0 - 6.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xxyy[i] = 2.0 * tr_y_xxyy[i] - 2.0 * tr_yyy_xxyy[i] * tbe_0 + 4.0 * tr_xy_xyy[i] - 4.0 * tr_xy_xxxyy[i] * tke_0 - 4.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyy_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxyy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xxyz[i] = 2.0 * tr_y_xxyz[i] - 2.0 * tr_yyy_xxyz[i] * tbe_0 + 4.0 * tr_xy_xyz[i] - 4.0 * tr_xy_xxxyz[i] * tke_0 - 4.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyy_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxyz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xxzz[i] = 2.0 * tr_y_xxzz[i] - 2.0 * tr_yyy_xxzz[i] * tbe_0 + 4.0 * tr_xy_xzz[i] - 4.0 * tr_xy_xxxzz[i] * tke_0 - 4.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xxzz[i] * tbe_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xyyy[i] = 2.0 * tr_y_xyyy[i] - 2.0 * tr_yyy_xyyy[i] * tbe_0 + 2.0 * tr_xy_yyy[i] - 4.0 * tr_xy_xxyyy[i] * tke_0 - 2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyy_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xyyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xyyz[i] = 2.0 * tr_y_xyyz[i] - 2.0 * tr_yyy_xyyz[i] * tbe_0 + 2.0 * tr_xy_yyz[i] - 4.0 * tr_xy_xxyyz[i] * tke_0 - 2.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyy_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xyyz[i] * tbe_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xyzz[i] = 2.0 * tr_y_xyzz[i] - 2.0 * tr_yyy_xyzz[i] * tbe_0 + 2.0 * tr_xy_yzz[i] - 4.0 * tr_xy_xxyzz[i] * tke_0 - 2.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyy_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xyzz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xzzz[i] = 2.0 * tr_y_xzzz[i] - 2.0 * tr_yyy_xzzz[i] * tbe_0 + 2.0 * tr_xy_zzz[i] - 4.0 * tr_xy_xxzzz[i] * tke_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xzzz[i] * tbe_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_yyyy[i] = 2.0 * tr_y_yyyy[i] - 2.0 * tr_yyy_yyyy[i] * tbe_0 - 4.0 * tr_xy_xyyyy[i] * tke_0 + 4.0 * tr_xyyy_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_yyyz[i] = 2.0 * tr_y_yyyz[i] - 2.0 * tr_yyy_yyyz[i] * tbe_0 - 4.0 * tr_xy_xyyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yyyz[i] * tbe_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_yyzz[i] = 2.0 * tr_y_yyzz[i] - 2.0 * tr_yyy_yyzz[i] * tbe_0 - 4.0 * tr_xy_xyyzz[i] * tke_0 + 4.0 * tr_xyyy_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yyzz[i] * tbe_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_yzzz[i] = 2.0 * tr_y_yzzz[i] - 2.0 * tr_yyy_yzzz[i] * tbe_0 - 4.0 * tr_xy_xyzzz[i] * tke_0 + 4.0 * tr_xyyy_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yzzz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_zzzz[i] = 2.0 * tr_y_zzzz[i] - 2.0 * tr_yyy_zzzz[i] * tbe_0 - 4.0 * tr_xy_xzzzz[i] * tke_0 + 4.0 * tr_xyyy_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-525 components of targeted buffer : FG

    auto tr_y_0_x_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 510);

    auto tr_y_0_x_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 511);

    auto tr_y_0_x_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 512);

    auto tr_y_0_x_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 513);

    auto tr_y_0_x_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 514);

    auto tr_y_0_x_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 515);

    auto tr_y_0_x_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 516);

    auto tr_y_0_x_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 517);

    auto tr_y_0_x_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 518);

    auto tr_y_0_x_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 519);

    auto tr_y_0_x_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 520);

    auto tr_y_0_x_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 521);

    auto tr_y_0_x_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 522);

    auto tr_y_0_x_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 523);

    auto tr_y_0_x_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 524);

    #pragma omp simd aligned(tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxx, tr_xyyz_xxxxy, tr_xyyz_xxxxz, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_y_0_x_xyz_xxxx, tr_y_0_x_xyz_xxxy, tr_y_0_x_xyz_xxxz, tr_y_0_x_xyz_xxyy, tr_y_0_x_xyz_xxyz, tr_y_0_x_xyz_xxzz, tr_y_0_x_xyz_xyyy, tr_y_0_x_xyz_xyyz, tr_y_0_x_xyz_xyzz, tr_y_0_x_xyz_xzzz, tr_y_0_x_xyz_yyyy, tr_y_0_x_xyz_yyyz, tr_y_0_x_xyz_yyzz, tr_y_0_x_xyz_yzzz, tr_y_0_x_xyz_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyz_xxxx[i] = tr_z_xxxx[i] - 2.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_xz_xxx[i] - 2.0 * tr_xz_xxxxx[i] * tke_0 - 8.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xxxy[i] = tr_z_xxxy[i] - 2.0 * tr_yyz_xxxy[i] * tbe_0 + 3.0 * tr_xz_xxy[i] - 2.0 * tr_xz_xxxxy[i] * tke_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xxxz[i] = tr_z_xxxz[i] - 2.0 * tr_yyz_xxxz[i] * tbe_0 + 3.0 * tr_xz_xxz[i] - 2.0 * tr_xz_xxxxz[i] * tke_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xxyy[i] = tr_z_xxyy[i] - 2.0 * tr_yyz_xxyy[i] * tbe_0 + 2.0 * tr_xz_xyy[i] - 2.0 * tr_xz_xxxyy[i] * tke_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xxyz[i] = tr_z_xxyz[i] - 2.0 * tr_yyz_xxyz[i] * tbe_0 + 2.0 * tr_xz_xyz[i] - 2.0 * tr_xz_xxxyz[i] * tke_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxyz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xxzz[i] = tr_z_xxzz[i] - 2.0 * tr_yyz_xxzz[i] * tbe_0 + 2.0 * tr_xz_xzz[i] - 2.0 * tr_xz_xxxzz[i] * tke_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xxzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xyyy[i] = tr_z_xyyy[i] - 2.0 * tr_yyz_xyyy[i] * tbe_0 + tr_xz_yyy[i] - 2.0 * tr_xz_xxyyy[i] * tke_0 - 2.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xyyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xyyz[i] = tr_z_xyyz[i] - 2.0 * tr_yyz_xyyz[i] * tbe_0 + tr_xz_yyz[i] - 2.0 * tr_xz_xxyyz[i] * tke_0 - 2.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xyyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xyzz[i] = tr_z_xyzz[i] - 2.0 * tr_yyz_xyzz[i] * tbe_0 + tr_xz_yzz[i] - 2.0 * tr_xz_xxyzz[i] * tke_0 - 2.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xyzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xzzz[i] = tr_z_xzzz[i] - 2.0 * tr_yyz_xzzz[i] * tbe_0 + tr_xz_zzz[i] - 2.0 * tr_xz_xxzzz[i] * tke_0 - 2.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xzzz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_yyyy[i] = tr_z_yyyy[i] - 2.0 * tr_yyz_yyyy[i] * tbe_0 - 2.0 * tr_xz_xyyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_yyyz[i] = tr_z_yyyz[i] - 2.0 * tr_yyz_yyyz[i] * tbe_0 - 2.0 * tr_xz_xyyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yyyz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_yyzz[i] = tr_z_yyzz[i] - 2.0 * tr_yyz_yyzz[i] * tbe_0 - 2.0 * tr_xz_xyyzz[i] * tke_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yyzz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_yzzz[i] = tr_z_yzzz[i] - 2.0 * tr_yyz_yzzz[i] * tbe_0 - 2.0 * tr_xz_xyzzz[i] * tke_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yzzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_zzzz[i] = tr_z_zzzz[i] - 2.0 * tr_yyz_zzzz[i] * tbe_0 - 2.0 * tr_xz_xzzzz[i] * tke_0 + 4.0 * tr_xyyz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 525-540 components of targeted buffer : FG

    auto tr_y_0_x_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 525);

    auto tr_y_0_x_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 526);

    auto tr_y_0_x_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 527);

    auto tr_y_0_x_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 528);

    auto tr_y_0_x_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 529);

    auto tr_y_0_x_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 530);

    auto tr_y_0_x_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 531);

    auto tr_y_0_x_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 532);

    auto tr_y_0_x_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 533);

    auto tr_y_0_x_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 534);

    auto tr_y_0_x_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 535);

    auto tr_y_0_x_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 536);

    auto tr_y_0_x_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 537);

    auto tr_y_0_x_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 538);

    auto tr_y_0_x_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 539);

    #pragma omp simd aligned(tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxx, tr_xyzz_xxxxy, tr_xyzz_xxxxz, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_y_0_x_xzz_xxxx, tr_y_0_x_xzz_xxxy, tr_y_0_x_xzz_xxxz, tr_y_0_x_xzz_xxyy, tr_y_0_x_xzz_xxyz, tr_y_0_x_xzz_xxzz, tr_y_0_x_xzz_xyyy, tr_y_0_x_xzz_xyyz, tr_y_0_x_xzz_xyzz, tr_y_0_x_xzz_xzzz, tr_y_0_x_xzz_yyyy, tr_y_0_x_xzz_yyyz, tr_y_0_x_xzz_yyzz, tr_y_0_x_xzz_yzzz, tr_y_0_x_xzz_zzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xzz_xxxx[i] = -2.0 * tr_yzz_xxxx[i] * tbe_0 - 8.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xxxy[i] = -2.0 * tr_yzz_xxxy[i] * tbe_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xxxz[i] = -2.0 * tr_yzz_xxxz[i] * tbe_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xxyy[i] = -2.0 * tr_yzz_xxyy[i] * tbe_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xxyz[i] = -2.0 * tr_yzz_xxyz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xxzz[i] = -2.0 * tr_yzz_xxzz[i] * tbe_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xyyy[i] = -2.0 * tr_yzz_xyyy[i] * tbe_0 - 2.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xyyz[i] = -2.0 * tr_yzz_xyyz[i] * tbe_0 - 2.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xyzz[i] = -2.0 * tr_yzz_xyzz[i] * tbe_0 - 2.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xzzz[i] = -2.0 * tr_yzz_xzzz[i] * tbe_0 - 2.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_yyyy[i] = -2.0 * tr_yzz_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_yyyz[i] = -2.0 * tr_yzz_yyyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_yyzz[i] = -2.0 * tr_yzz_yyzz[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_yzzz[i] = -2.0 * tr_yzz_yzzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_zzzz[i] = -2.0 * tr_yzz_zzzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 540-555 components of targeted buffer : FG

    auto tr_y_0_x_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 540);

    auto tr_y_0_x_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 541);

    auto tr_y_0_x_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 542);

    auto tr_y_0_x_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 543);

    auto tr_y_0_x_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 544);

    auto tr_y_0_x_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 545);

    auto tr_y_0_x_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 546);

    auto tr_y_0_x_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 547);

    auto tr_y_0_x_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 548);

    auto tr_y_0_x_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 549);

    auto tr_y_0_x_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 550);

    auto tr_y_0_x_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 551);

    auto tr_y_0_x_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 552);

    auto tr_y_0_x_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 553);

    auto tr_y_0_x_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 554);

    #pragma omp simd aligned(tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xzzz, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yzzz, tr_xyyyy_zzzz, tr_y_0_x_yyy_xxxx, tr_y_0_x_yyy_xxxy, tr_y_0_x_yyy_xxxz, tr_y_0_x_yyy_xxyy, tr_y_0_x_yyy_xxyz, tr_y_0_x_yyy_xxzz, tr_y_0_x_yyy_xyyy, tr_y_0_x_yyy_xyyz, tr_y_0_x_yyy_xyzz, tr_y_0_x_yyy_xzzz, tr_y_0_x_yyy_yyyy, tr_y_0_x_yyy_yyyz, tr_y_0_x_yyy_yyzz, tr_y_0_x_yyy_yzzz, tr_y_0_x_yyy_zzzz, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyyy_xxx, tr_yyyy_xxxxx, tr_yyyy_xxxxy, tr_yyyy_xxxxz, tr_yyyy_xxxyy, tr_yyyy_xxxyz, tr_yyyy_xxxzz, tr_yyyy_xxy, tr_yyyy_xxyyy, tr_yyyy_xxyyz, tr_yyyy_xxyzz, tr_yyyy_xxz, tr_yyyy_xxzzz, tr_yyyy_xyy, tr_yyyy_xyyyy, tr_yyyy_xyyyz, tr_yyyy_xyyzz, tr_yyyy_xyz, tr_yyyy_xyzzz, tr_yyyy_xzz, tr_yyyy_xzzzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyy_xxxx[i] = 12.0 * tr_yy_xxx[i] - 6.0 * tr_yy_xxxxx[i] * tke_0 - 8.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyy_xxxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xxxy[i] = 9.0 * tr_yy_xxy[i] - 6.0 * tr_yy_xxxxy[i] * tke_0 - 6.0 * tr_yyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyy_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xxxz[i] = 9.0 * tr_yy_xxz[i] - 6.0 * tr_yy_xxxxz[i] * tke_0 - 6.0 * tr_yyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyy_xxxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xxyy[i] = 6.0 * tr_yy_xyy[i] - 6.0 * tr_yy_xxxyy[i] * tke_0 - 4.0 * tr_yyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyy_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xxyz[i] = 6.0 * tr_yy_xyz[i] - 6.0 * tr_yy_xxxyz[i] * tke_0 - 4.0 * tr_yyyy_xyz[i] * tbe_0 + 4.0 * tr_yyyy_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxyz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xxzz[i] = 6.0 * tr_yy_xzz[i] - 6.0 * tr_yy_xxxzz[i] * tke_0 - 4.0 * tr_yyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyy_xxxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xxzz[i] * tbe_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xyyy[i] = 3.0 * tr_yy_yyy[i] - 6.0 * tr_yy_xxyyy[i] * tke_0 - 2.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyy_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xyyz[i] = 3.0 * tr_yy_yyz[i] - 6.0 * tr_yy_xxyyz[i] * tke_0 - 2.0 * tr_yyyy_yyz[i] * tbe_0 + 4.0 * tr_yyyy_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xyyz[i] * tbe_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xyzz[i] = 3.0 * tr_yy_yzz[i] - 6.0 * tr_yy_xxyzz[i] * tke_0 - 2.0 * tr_yyyy_yzz[i] * tbe_0 + 4.0 * tr_yyyy_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xyzz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xzzz[i] = 3.0 * tr_yy_zzz[i] - 6.0 * tr_yy_xxzzz[i] * tke_0 - 2.0 * tr_yyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyy_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xzzz[i] * tbe_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_yyyy[i] = -6.0 * tr_yy_xyyyy[i] * tke_0 + 4.0 * tr_yyyy_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_yyyz[i] = -6.0 * tr_yy_xyyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yyyz[i] * tbe_0 + 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_yyzz[i] = -6.0 * tr_yy_xyyzz[i] * tke_0 + 4.0 * tr_yyyy_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yyzz[i] * tbe_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_yzzz[i] = -6.0 * tr_yy_xyzzz[i] * tke_0 + 4.0 * tr_yyyy_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yzzz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_zzzz[i] = -6.0 * tr_yy_xzzzz[i] * tke_0 + 4.0 * tr_yyyy_xzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 555-570 components of targeted buffer : FG

    auto tr_y_0_x_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 555);

    auto tr_y_0_x_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 556);

    auto tr_y_0_x_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 557);

    auto tr_y_0_x_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 558);

    auto tr_y_0_x_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 559);

    auto tr_y_0_x_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 560);

    auto tr_y_0_x_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 561);

    auto tr_y_0_x_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 562);

    auto tr_y_0_x_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 563);

    auto tr_y_0_x_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 564);

    auto tr_y_0_x_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 565);

    auto tr_y_0_x_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 566);

    auto tr_y_0_x_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 567);

    auto tr_y_0_x_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 568);

    auto tr_y_0_x_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 569);

    #pragma omp simd aligned(tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_y_0_x_yyz_xxxx, tr_y_0_x_yyz_xxxy, tr_y_0_x_yyz_xxxz, tr_y_0_x_yyz_xxyy, tr_y_0_x_yyz_xxyz, tr_y_0_x_yyz_xxzz, tr_y_0_x_yyz_xyyy, tr_y_0_x_yyz_xyyz, tr_y_0_x_yyz_xyzz, tr_y_0_x_yyz_xzzz, tr_y_0_x_yyz_yyyy, tr_y_0_x_yyz_yyyz, tr_y_0_x_yyz_yyzz, tr_y_0_x_yyz_yzzz, tr_y_0_x_yyz_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxx, tr_yyyz_xxxxy, tr_yyyz_xxxxz, tr_yyyz_xxxyy, tr_yyyz_xxxyz, tr_yyyz_xxxzz, tr_yyyz_xxy, tr_yyyz_xxyyy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xxzzz, tr_yyyz_xyy, tr_yyyz_xyyyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_xzzzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyz_xxxx[i] = 8.0 * tr_yz_xxx[i] - 4.0 * tr_yz_xxxxx[i] * tke_0 - 8.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xxxy[i] = 6.0 * tr_yz_xxy[i] - 4.0 * tr_yz_xxxxy[i] * tke_0 - 6.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xxxz[i] = 6.0 * tr_yz_xxz[i] - 4.0 * tr_yz_xxxxz[i] * tke_0 - 6.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xxyy[i] = 4.0 * tr_yz_xyy[i] - 4.0 * tr_yz_xxxyy[i] * tke_0 - 4.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xxyz[i] = 4.0 * tr_yz_xyz[i] - 4.0 * tr_yz_xxxyz[i] * tke_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxyz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xxzz[i] = 4.0 * tr_yz_xzz[i] - 4.0 * tr_yz_xxxzz[i] * tke_0 - 4.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xyyy[i] = 2.0 * tr_yz_yyy[i] - 4.0 * tr_yz_xxyyy[i] * tke_0 - 2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xyyz[i] = 2.0 * tr_yz_yyz[i] - 4.0 * tr_yz_xxyyz[i] * tke_0 - 2.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xyzz[i] = 2.0 * tr_yz_yzz[i] - 4.0 * tr_yz_xxyzz[i] * tke_0 - 2.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xzzz[i] = 2.0 * tr_yz_zzz[i] - 4.0 * tr_yz_xxzzz[i] * tke_0 - 2.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_yyyy[i] = -4.0 * tr_yz_xyyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_yyyz[i] = -4.0 * tr_yz_xyyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyyz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_yyzz[i] = -4.0 * tr_yz_xyyzz[i] * tke_0 + 4.0 * tr_yyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyzz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_yzzz[i] = -4.0 * tr_yz_xyzzz[i] * tke_0 + 4.0 * tr_yyyz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yzzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_zzzz[i] = -4.0 * tr_yz_xzzzz[i] * tke_0 + 4.0 * tr_yyyz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 570-585 components of targeted buffer : FG

    auto tr_y_0_x_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 570);

    auto tr_y_0_x_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 571);

    auto tr_y_0_x_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 572);

    auto tr_y_0_x_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 573);

    auto tr_y_0_x_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 574);

    auto tr_y_0_x_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 575);

    auto tr_y_0_x_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 576);

    auto tr_y_0_x_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 577);

    auto tr_y_0_x_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 578);

    auto tr_y_0_x_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 579);

    auto tr_y_0_x_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 580);

    auto tr_y_0_x_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 581);

    auto tr_y_0_x_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 582);

    auto tr_y_0_x_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 583);

    auto tr_y_0_x_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 584);

    #pragma omp simd aligned(tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_y_0_x_yzz_xxxx, tr_y_0_x_yzz_xxxy, tr_y_0_x_yzz_xxxz, tr_y_0_x_yzz_xxyy, tr_y_0_x_yzz_xxyz, tr_y_0_x_yzz_xxzz, tr_y_0_x_yzz_xyyy, tr_y_0_x_yzz_xyyz, tr_y_0_x_yzz_xyzz, tr_y_0_x_yzz_xzzz, tr_y_0_x_yzz_yyyy, tr_y_0_x_yzz_yyyz, tr_y_0_x_yzz_yyzz, tr_y_0_x_yzz_yzzz, tr_y_0_x_yzz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxx, tr_yyzz_xxxxy, tr_yyzz_xxxxz, tr_yyzz_xxxyy, tr_yyzz_xxxyz, tr_yyzz_xxxzz, tr_yyzz_xxy, tr_yyzz_xxyyy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xxzzz, tr_yyzz_xyy, tr_yyzz_xyyyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_xzzzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yzz_xxxx[i] = 4.0 * tr_zz_xxx[i] - 2.0 * tr_zz_xxxxx[i] * tke_0 - 8.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xxxy[i] = 3.0 * tr_zz_xxy[i] - 2.0 * tr_zz_xxxxy[i] * tke_0 - 6.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xxxz[i] = 3.0 * tr_zz_xxz[i] - 2.0 * tr_zz_xxxxz[i] * tke_0 - 6.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xxyy[i] = 2.0 * tr_zz_xyy[i] - 2.0 * tr_zz_xxxyy[i] * tke_0 - 4.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xxyz[i] = 2.0 * tr_zz_xyz[i] - 2.0 * tr_zz_xxxyz[i] * tke_0 - 4.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xxzz[i] = 2.0 * tr_zz_xzz[i] - 2.0 * tr_zz_xxxzz[i] * tke_0 - 4.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xxzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xyyy[i] = tr_zz_yyy[i] - 2.0 * tr_zz_xxyyy[i] * tke_0 - 2.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xyyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xyyz[i] = tr_zz_yyz[i] - 2.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xyyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xyzz[i] = tr_zz_yzz[i] - 2.0 * tr_zz_xxyzz[i] * tke_0 - 2.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xyzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xzzz[i] = tr_zz_zzz[i] - 2.0 * tr_zz_xxzzz[i] * tke_0 - 2.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xzzz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_yyyy[i] = -2.0 * tr_zz_xyyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yyyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_yyyz[i] = -2.0 * tr_zz_xyyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yyyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_yyzz[i] = -2.0 * tr_zz_xyyzz[i] * tke_0 + 4.0 * tr_yyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yyzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_yzzz[i] = -2.0 * tr_zz_xyzzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yzzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_zzzz[i] = -2.0 * tr_zz_xzzzz[i] * tke_0 + 4.0 * tr_yyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_zzzz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 585-600 components of targeted buffer : FG

    auto tr_y_0_x_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 585);

    auto tr_y_0_x_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 586);

    auto tr_y_0_x_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 587);

    auto tr_y_0_x_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 588);

    auto tr_y_0_x_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 589);

    auto tr_y_0_x_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 590);

    auto tr_y_0_x_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 591);

    auto tr_y_0_x_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 592);

    auto tr_y_0_x_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 593);

    auto tr_y_0_x_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 594);

    auto tr_y_0_x_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 595);

    auto tr_y_0_x_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 596);

    auto tr_y_0_x_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 597);

    auto tr_y_0_x_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 598);

    auto tr_y_0_x_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 599);

    #pragma omp simd aligned(tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_y_0_x_zzz_xxxx, tr_y_0_x_zzz_xxxy, tr_y_0_x_zzz_xxxz, tr_y_0_x_zzz_xxyy, tr_y_0_x_zzz_xxyz, tr_y_0_x_zzz_xxzz, tr_y_0_x_zzz_xyyy, tr_y_0_x_zzz_xyyz, tr_y_0_x_zzz_xyzz, tr_y_0_x_zzz_xzzz, tr_y_0_x_zzz_yyyy, tr_y_0_x_zzz_yyyz, tr_y_0_x_zzz_yyzz, tr_y_0_x_zzz_yzzz, tr_y_0_x_zzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxx, tr_yzzz_xxxxy, tr_yzzz_xxxxz, tr_yzzz_xxxyy, tr_yzzz_xxxyz, tr_yzzz_xxxzz, tr_yzzz_xxy, tr_yzzz_xxyyy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xxzzz, tr_yzzz_xyy, tr_yzzz_xyyyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_xzzzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zzz_xxxx[i] = -8.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xxxy[i] = -6.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xxxz[i] = -6.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xxyy[i] = -4.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xxyz[i] = -4.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xxzz[i] = -4.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xyyy[i] = -2.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xyyz[i] = -2.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xyzz[i] = -2.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xzzz[i] = -2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_yyyy[i] = 4.0 * tr_yzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_yyyz[i] = 4.0 * tr_yzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_yyzz[i] = 4.0 * tr_yzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_yzzz[i] = 4.0 * tr_yzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_zzzz[i] = 4.0 * tr_yzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 600-615 components of targeted buffer : FG

    auto tr_y_0_y_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 600);

    auto tr_y_0_y_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 601);

    auto tr_y_0_y_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 602);

    auto tr_y_0_y_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 603);

    auto tr_y_0_y_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 604);

    auto tr_y_0_y_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 605);

    auto tr_y_0_y_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 606);

    auto tr_y_0_y_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 607);

    auto tr_y_0_y_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 608);

    auto tr_y_0_y_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 609);

    auto tr_y_0_y_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 610);

    auto tr_y_0_y_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 611);

    auto tr_y_0_y_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 612);

    auto tr_y_0_y_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 613);

    auto tr_y_0_y_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 614);

    #pragma omp simd aligned(tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxy_xxx, tr_xxxy_xxxxy, tr_xxxy_xxxyy, tr_xxxy_xxxyz, tr_xxxy_xxy, tr_xxxy_xxyyy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyyyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyyyy, tr_xxxy_yyyyz, tr_xxxy_yyyzz, tr_xxxy_yyz, tr_xxxy_yyzzz, tr_xxxy_yzz, tr_xxxy_yzzzz, tr_xxxy_zzz, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xzzz, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yzzz, tr_xxxyy_zzzz, tr_y_0_y_xxx_xxxx, tr_y_0_y_xxx_xxxy, tr_y_0_y_xxx_xxxz, tr_y_0_y_xxx_xxyy, tr_y_0_y_xxx_xxyz, tr_y_0_y_xxx_xxzz, tr_y_0_y_xxx_xyyy, tr_y_0_y_xxx_xyyz, tr_y_0_y_xxx_xyzz, tr_y_0_y_xxx_xzzz, tr_y_0_y_xxx_yyyy, tr_y_0_y_xxx_yyyz, tr_y_0_y_xxx_yyzz, tr_y_0_y_xxx_yzzz, tr_y_0_y_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxx_xxxx[i] = -2.0 * tr_xxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xxxy[i] = -2.0 * tr_xxx_xxxy[i] * tbe_0 - 2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xxxz[i] = -2.0 * tr_xxx_xxxz[i] * tbe_0 + 4.0 * tr_xxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xxyy[i] = -2.0 * tr_xxx_xxyy[i] * tbe_0 - 4.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xxyz[i] = -2.0 * tr_xxx_xxyz[i] * tbe_0 - 2.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xxzz[i] = -2.0 * tr_xxx_xxzz[i] * tbe_0 + 4.0 * tr_xxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xyyy[i] = -2.0 * tr_xxx_xyyy[i] * tbe_0 - 6.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xyyz[i] = -2.0 * tr_xxx_xyyz[i] * tbe_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xyzz[i] = -2.0 * tr_xxx_xyzz[i] * tbe_0 - 2.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xzzz[i] = -2.0 * tr_xxx_xzzz[i] * tbe_0 + 4.0 * tr_xxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_yyyy[i] = -2.0 * tr_xxx_yyyy[i] * tbe_0 - 8.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_yyyz[i] = -2.0 * tr_xxx_yyyz[i] * tbe_0 - 6.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_yyzz[i] = -2.0 * tr_xxx_yyzz[i] * tbe_0 - 4.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_yzzz[i] = -2.0 * tr_xxx_yzzz[i] * tbe_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_zzzz[i] = -2.0 * tr_xxx_zzzz[i] * tbe_0 + 4.0 * tr_xxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 615-630 components of targeted buffer : FG

    auto tr_y_0_y_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 615);

    auto tr_y_0_y_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 616);

    auto tr_y_0_y_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 617);

    auto tr_y_0_y_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 618);

    auto tr_y_0_y_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 619);

    auto tr_y_0_y_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 620);

    auto tr_y_0_y_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 621);

    auto tr_y_0_y_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 622);

    auto tr_y_0_y_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 623);

    auto tr_y_0_y_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 624);

    auto tr_y_0_y_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 625);

    auto tr_y_0_y_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 626);

    auto tr_y_0_y_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 627);

    auto tr_y_0_y_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 628);

    auto tr_y_0_y_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 629);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyy_xxx, tr_xxyy_xxxxy, tr_xxyy_xxxyy, tr_xxyy_xxxyz, tr_xxyy_xxy, tr_xxyy_xxyyy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyyyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyyyy, tr_xxyy_yyyyz, tr_xxyy_yyyzz, tr_xxyy_yyz, tr_xxyy_yyzzz, tr_xxyy_yzz, tr_xxyy_yzzzz, tr_xxyy_zzz, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xzzz, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yzzz, tr_xxyyy_zzzz, tr_y_0_y_xxy_xxxx, tr_y_0_y_xxy_xxxy, tr_y_0_y_xxy_xxxz, tr_y_0_y_xxy_xxyy, tr_y_0_y_xxy_xxyz, tr_y_0_y_xxy_xxzz, tr_y_0_y_xxy_xyyy, tr_y_0_y_xxy_xyyz, tr_y_0_y_xxy_xyzz, tr_y_0_y_xxy_xzzz, tr_y_0_y_xxy_yyyy, tr_y_0_y_xxy_yyyz, tr_y_0_y_xxy_yyzz, tr_y_0_y_xxy_yzzz, tr_y_0_y_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxy_xxxx[i] = -2.0 * tr_xx_xxxxy[i] * tke_0 - 6.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xxxy[i] = tr_xx_xxx[i] - 2.0 * tr_xx_xxxyy[i] * tke_0 - 6.0 * tr_xxy_xxxy[i] * tbe_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xxxz[i] = -2.0 * tr_xx_xxxyz[i] * tke_0 - 6.0 * tr_xxy_xxxz[i] * tbe_0 + 4.0 * tr_xxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xxyy[i] = 2.0 * tr_xx_xxy[i] - 2.0 * tr_xx_xxyyy[i] * tke_0 - 6.0 * tr_xxy_xxyy[i] * tbe_0 - 4.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xxyz[i] = tr_xx_xxz[i] - 2.0 * tr_xx_xxyyz[i] * tke_0 - 6.0 * tr_xxy_xxyz[i] * tbe_0 - 2.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xxzz[i] = -2.0 * tr_xx_xxyzz[i] * tke_0 - 6.0 * tr_xxy_xxzz[i] * tbe_0 + 4.0 * tr_xxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xyyy[i] = 3.0 * tr_xx_xyy[i] - 2.0 * tr_xx_xyyyy[i] * tke_0 - 6.0 * tr_xxy_xyyy[i] * tbe_0 - 6.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xyyz[i] = 2.0 * tr_xx_xyz[i] - 2.0 * tr_xx_xyyyz[i] * tke_0 - 6.0 * tr_xxy_xyyz[i] * tbe_0 - 4.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xyzz[i] = tr_xx_xzz[i] - 2.0 * tr_xx_xyyzz[i] * tke_0 - 6.0 * tr_xxy_xyzz[i] * tbe_0 - 2.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xzzz[i] = -2.0 * tr_xx_xyzzz[i] * tke_0 - 6.0 * tr_xxy_xzzz[i] * tbe_0 + 4.0 * tr_xxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_yyyy[i] = 4.0 * tr_xx_yyy[i] - 2.0 * tr_xx_yyyyy[i] * tke_0 - 6.0 * tr_xxy_yyyy[i] * tbe_0 - 8.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_yyyz[i] = 3.0 * tr_xx_yyz[i] - 2.0 * tr_xx_yyyyz[i] * tke_0 - 6.0 * tr_xxy_yyyz[i] * tbe_0 - 6.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_yyzz[i] = 2.0 * tr_xx_yzz[i] - 2.0 * tr_xx_yyyzz[i] * tke_0 - 6.0 * tr_xxy_yyzz[i] * tbe_0 - 4.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_yzzz[i] = tr_xx_zzz[i] - 2.0 * tr_xx_yyzzz[i] * tke_0 - 6.0 * tr_xxy_yzzz[i] * tbe_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_zzzz[i] = -2.0 * tr_xx_yzzzz[i] * tke_0 - 6.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 630-645 components of targeted buffer : FG

    auto tr_y_0_y_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 630);

    auto tr_y_0_y_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 631);

    auto tr_y_0_y_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 632);

    auto tr_y_0_y_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 633);

    auto tr_y_0_y_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 634);

    auto tr_y_0_y_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 635);

    auto tr_y_0_y_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 636);

    auto tr_y_0_y_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 637);

    auto tr_y_0_y_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 638);

    auto tr_y_0_y_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 639);

    auto tr_y_0_y_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 640);

    auto tr_y_0_y_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 641);

    auto tr_y_0_y_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 642);

    auto tr_y_0_y_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 643);

    auto tr_y_0_y_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 644);

    #pragma omp simd aligned(tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxy, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyyyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_y_0_y_xxz_xxxx, tr_y_0_y_xxz_xxxy, tr_y_0_y_xxz_xxxz, tr_y_0_y_xxz_xxyy, tr_y_0_y_xxz_xxyz, tr_y_0_y_xxz_xxzz, tr_y_0_y_xxz_xyyy, tr_y_0_y_xxz_xyyz, tr_y_0_y_xxz_xyzz, tr_y_0_y_xxz_xzzz, tr_y_0_y_xxz_yyyy, tr_y_0_y_xxz_yyyz, tr_y_0_y_xxz_yyzz, tr_y_0_y_xxz_yzzz, tr_y_0_y_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxz_xxxx[i] = -2.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xxxy[i] = -2.0 * tr_xxz_xxxy[i] * tbe_0 - 2.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xxxz[i] = -2.0 * tr_xxz_xxxz[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xxyy[i] = -2.0 * tr_xxz_xxyy[i] * tbe_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xxyz[i] = -2.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xxzz[i] = -2.0 * tr_xxz_xxzz[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xyyy[i] = -2.0 * tr_xxz_xyyy[i] * tbe_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xyyz[i] = -2.0 * tr_xxz_xyyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xyzz[i] = -2.0 * tr_xxz_xyzz[i] * tbe_0 - 2.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xzzz[i] = -2.0 * tr_xxz_xzzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_yyyy[i] = -2.0 * tr_xxz_yyyy[i] * tbe_0 - 8.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_yyyz[i] = -2.0 * tr_xxz_yyyz[i] * tbe_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_yyzz[i] = -2.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_yzzz[i] = -2.0 * tr_xxz_yzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_zzzz[i] = -2.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 645-660 components of targeted buffer : FG

    auto tr_y_0_y_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 645);

    auto tr_y_0_y_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 646);

    auto tr_y_0_y_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 647);

    auto tr_y_0_y_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 648);

    auto tr_y_0_y_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 649);

    auto tr_y_0_y_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 650);

    auto tr_y_0_y_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 651);

    auto tr_y_0_y_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 652);

    auto tr_y_0_y_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 653);

    auto tr_y_0_y_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 654);

    auto tr_y_0_y_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 655);

    auto tr_y_0_y_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 656);

    auto tr_y_0_y_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 657);

    auto tr_y_0_y_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 658);

    auto tr_y_0_y_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 659);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyy_xxx, tr_xyyy_xxxxy, tr_xyyy_xxxyy, tr_xyyy_xxxyz, tr_xyyy_xxy, tr_xyyy_xxyyy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyyyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyyyy, tr_xyyy_yyyyz, tr_xyyy_yyyzz, tr_xyyy_yyz, tr_xyyy_yyzzz, tr_xyyy_yzz, tr_xyyy_yzzzz, tr_xyyy_zzz, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xzzz, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yzzz, tr_xyyyy_zzzz, tr_y_0_y_xyy_xxxx, tr_y_0_y_xyy_xxxy, tr_y_0_y_xyy_xxxz, tr_y_0_y_xyy_xxyy, tr_y_0_y_xyy_xxyz, tr_y_0_y_xyy_xxzz, tr_y_0_y_xyy_xyyy, tr_y_0_y_xyy_xyyz, tr_y_0_y_xyy_xyzz, tr_y_0_y_xyy_xzzz, tr_y_0_y_xyy_yyyy, tr_y_0_y_xyy_yyyz, tr_y_0_y_xyy_yyzz, tr_y_0_y_xyy_yzzz, tr_y_0_y_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyy_xxxx[i] = 2.0 * tr_x_xxxx[i] - 4.0 * tr_xy_xxxxy[i] * tke_0 - 10.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xxxy[i] = 2.0 * tr_x_xxxy[i] + 2.0 * tr_xy_xxx[i] - 4.0 * tr_xy_xxxyy[i] * tke_0 - 10.0 * tr_xyy_xxxy[i] * tbe_0 - 2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xxxz[i] = 2.0 * tr_x_xxxz[i] - 4.0 * tr_xy_xxxyz[i] * tke_0 - 10.0 * tr_xyy_xxxz[i] * tbe_0 + 4.0 * tr_xyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xxyy[i] = 2.0 * tr_x_xxyy[i] + 4.0 * tr_xy_xxy[i] - 4.0 * tr_xy_xxyyy[i] * tke_0 - 10.0 * tr_xyy_xxyy[i] * tbe_0 - 4.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xxyz[i] = 2.0 * tr_x_xxyz[i] + 2.0 * tr_xy_xxz[i] - 4.0 * tr_xy_xxyyz[i] * tke_0 - 10.0 * tr_xyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xxzz[i] = 2.0 * tr_x_xxzz[i] - 4.0 * tr_xy_xxyzz[i] * tke_0 - 10.0 * tr_xyy_xxzz[i] * tbe_0 + 4.0 * tr_xyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xyyy[i] = 2.0 * tr_x_xyyy[i] + 6.0 * tr_xy_xyy[i] - 4.0 * tr_xy_xyyyy[i] * tke_0 - 10.0 * tr_xyy_xyyy[i] * tbe_0 - 6.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xyyz[i] = 2.0 * tr_x_xyyz[i] + 4.0 * tr_xy_xyz[i] - 4.0 * tr_xy_xyyyz[i] * tke_0 - 10.0 * tr_xyy_xyyz[i] * tbe_0 - 4.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xyzz[i] = 2.0 * tr_x_xyzz[i] + 2.0 * tr_xy_xzz[i] - 4.0 * tr_xy_xyyzz[i] * tke_0 - 10.0 * tr_xyy_xyzz[i] * tbe_0 - 2.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xzzz[i] = 2.0 * tr_x_xzzz[i] - 4.0 * tr_xy_xyzzz[i] * tke_0 - 10.0 * tr_xyy_xzzz[i] * tbe_0 + 4.0 * tr_xyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_yyyy[i] = 2.0 * tr_x_yyyy[i] + 8.0 * tr_xy_yyy[i] - 4.0 * tr_xy_yyyyy[i] * tke_0 - 10.0 * tr_xyy_yyyy[i] * tbe_0 - 8.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_yyyz[i] = 2.0 * tr_x_yyyz[i] + 6.0 * tr_xy_yyz[i] - 4.0 * tr_xy_yyyyz[i] * tke_0 - 10.0 * tr_xyy_yyyz[i] * tbe_0 - 6.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_yyzz[i] = 2.0 * tr_x_yyzz[i] + 4.0 * tr_xy_yzz[i] - 4.0 * tr_xy_yyyzz[i] * tke_0 - 10.0 * tr_xyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_yzzz[i] = 2.0 * tr_x_yzzz[i] + 2.0 * tr_xy_zzz[i] - 4.0 * tr_xy_yyzzz[i] * tke_0 - 10.0 * tr_xyy_yzzz[i] * tbe_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_zzzz[i] = 2.0 * tr_x_zzzz[i] - 4.0 * tr_xy_yzzzz[i] * tke_0 - 10.0 * tr_xyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 660-675 components of targeted buffer : FG

    auto tr_y_0_y_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 660);

    auto tr_y_0_y_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 661);

    auto tr_y_0_y_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 662);

    auto tr_y_0_y_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 663);

    auto tr_y_0_y_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 664);

    auto tr_y_0_y_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 665);

    auto tr_y_0_y_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 666);

    auto tr_y_0_y_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 667);

    auto tr_y_0_y_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 668);

    auto tr_y_0_y_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 669);

    auto tr_y_0_y_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 670);

    auto tr_y_0_y_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 671);

    auto tr_y_0_y_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 672);

    auto tr_y_0_y_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 673);

    auto tr_y_0_y_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 674);

    #pragma omp simd aligned(tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxy, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyyyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_y_0_y_xyz_xxxx, tr_y_0_y_xyz_xxxy, tr_y_0_y_xyz_xxxz, tr_y_0_y_xyz_xxyy, tr_y_0_y_xyz_xxyz, tr_y_0_y_xyz_xxzz, tr_y_0_y_xyz_xyyy, tr_y_0_y_xyz_xyyz, tr_y_0_y_xyz_xyzz, tr_y_0_y_xyz_xzzz, tr_y_0_y_xyz_yyyy, tr_y_0_y_xyz_yyyz, tr_y_0_y_xyz_yyzz, tr_y_0_y_xyz_yzzz, tr_y_0_y_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyz_xxxx[i] = -2.0 * tr_xz_xxxxy[i] * tke_0 - 6.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xxxy[i] = tr_xz_xxx[i] - 2.0 * tr_xz_xxxyy[i] * tke_0 - 6.0 * tr_xyz_xxxy[i] * tbe_0 - 2.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xxxz[i] = -2.0 * tr_xz_xxxyz[i] * tke_0 - 6.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xxyy[i] = 2.0 * tr_xz_xxy[i] - 2.0 * tr_xz_xxyyy[i] * tke_0 - 6.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xxyz[i] = tr_xz_xxz[i] - 2.0 * tr_xz_xxyyz[i] * tke_0 - 6.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xxzz[i] = -2.0 * tr_xz_xxyzz[i] * tke_0 - 6.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xyyy[i] = 3.0 * tr_xz_xyy[i] - 2.0 * tr_xz_xyyyy[i] * tke_0 - 6.0 * tr_xyz_xyyy[i] * tbe_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xyyz[i] = 2.0 * tr_xz_xyz[i] - 2.0 * tr_xz_xyyyz[i] * tke_0 - 6.0 * tr_xyz_xyyz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xyzz[i] = tr_xz_xzz[i] - 2.0 * tr_xz_xyyzz[i] * tke_0 - 6.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xzzz[i] = -2.0 * tr_xz_xyzzz[i] * tke_0 - 6.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_yyyy[i] = 4.0 * tr_xz_yyy[i] - 2.0 * tr_xz_yyyyy[i] * tke_0 - 6.0 * tr_xyz_yyyy[i] * tbe_0 - 8.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_yyyz[i] = 3.0 * tr_xz_yyz[i] - 2.0 * tr_xz_yyyyz[i] * tke_0 - 6.0 * tr_xyz_yyyz[i] * tbe_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_yyzz[i] = 2.0 * tr_xz_yzz[i] - 2.0 * tr_xz_yyyzz[i] * tke_0 - 6.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_yzzz[i] = tr_xz_zzz[i] - 2.0 * tr_xz_yyzzz[i] * tke_0 - 6.0 * tr_xyz_yzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_zzzz[i] = -2.0 * tr_xz_yzzzz[i] * tke_0 - 6.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 675-690 components of targeted buffer : FG

    auto tr_y_0_y_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 675);

    auto tr_y_0_y_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 676);

    auto tr_y_0_y_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 677);

    auto tr_y_0_y_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 678);

    auto tr_y_0_y_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 679);

    auto tr_y_0_y_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 680);

    auto tr_y_0_y_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 681);

    auto tr_y_0_y_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 682);

    auto tr_y_0_y_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 683);

    auto tr_y_0_y_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 684);

    auto tr_y_0_y_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 685);

    auto tr_y_0_y_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 686);

    auto tr_y_0_y_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 687);

    auto tr_y_0_y_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 688);

    auto tr_y_0_y_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 689);

    #pragma omp simd aligned(tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxy, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyyyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_y_0_y_xzz_xxxx, tr_y_0_y_xzz_xxxy, tr_y_0_y_xzz_xxxz, tr_y_0_y_xzz_xxyy, tr_y_0_y_xzz_xxyz, tr_y_0_y_xzz_xxzz, tr_y_0_y_xzz_xyyy, tr_y_0_y_xzz_xyyz, tr_y_0_y_xzz_xyzz, tr_y_0_y_xzz_xzzz, tr_y_0_y_xzz_yyyy, tr_y_0_y_xzz_yyyz, tr_y_0_y_xzz_yyzz, tr_y_0_y_xzz_yzzz, tr_y_0_y_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xzz_xxxx[i] = -2.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xxxy[i] = -2.0 * tr_xzz_xxxy[i] * tbe_0 - 2.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xxxz[i] = -2.0 * tr_xzz_xxxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xxyy[i] = -2.0 * tr_xzz_xxyy[i] * tbe_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xxyz[i] = -2.0 * tr_xzz_xxyz[i] * tbe_0 - 2.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xxzz[i] = -2.0 * tr_xzz_xxzz[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xyyy[i] = -2.0 * tr_xzz_xyyy[i] * tbe_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xyyz[i] = -2.0 * tr_xzz_xyyz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xyzz[i] = -2.0 * tr_xzz_xyzz[i] * tbe_0 - 2.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xzzz[i] = -2.0 * tr_xzz_xzzz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_yyyy[i] = -2.0 * tr_xzz_yyyy[i] * tbe_0 - 8.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_yyyz[i] = -2.0 * tr_xzz_yyyz[i] * tbe_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_yyzz[i] = -2.0 * tr_xzz_yyzz[i] * tbe_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_yzzz[i] = -2.0 * tr_xzz_yzzz[i] * tbe_0 - 2.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_zzzz[i] = -2.0 * tr_xzz_zzzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 690-705 components of targeted buffer : FG

    auto tr_y_0_y_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 690);

    auto tr_y_0_y_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 691);

    auto tr_y_0_y_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 692);

    auto tr_y_0_y_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 693);

    auto tr_y_0_y_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 694);

    auto tr_y_0_y_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 695);

    auto tr_y_0_y_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 696);

    auto tr_y_0_y_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 697);

    auto tr_y_0_y_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 698);

    auto tr_y_0_y_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 699);

    auto tr_y_0_y_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 700);

    auto tr_y_0_y_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 701);

    auto tr_y_0_y_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 702);

    auto tr_y_0_y_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 703);

    auto tr_y_0_y_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 704);

    #pragma omp simd aligned(tr_y_0_y_yyy_xxxx, tr_y_0_y_yyy_xxxy, tr_y_0_y_yyy_xxxz, tr_y_0_y_yyy_xxyy, tr_y_0_y_yyy_xxyz, tr_y_0_y_yyy_xxzz, tr_y_0_y_yyy_xyyy, tr_y_0_y_yyy_xyyz, tr_y_0_y_yyy_xyzz, tr_y_0_y_yyy_xzzz, tr_y_0_y_yyy_yyyy, tr_y_0_y_yyy_yyyz, tr_y_0_y_yyy_yyzz, tr_y_0_y_yyy_yzzz, tr_y_0_y_yyy_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, tr_yyyy_xxx, tr_yyyy_xxxxy, tr_yyyy_xxxyy, tr_yyyy_xxxyz, tr_yyyy_xxy, tr_yyyy_xxyyy, tr_yyyy_xxyyz, tr_yyyy_xxyzz, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyyyy, tr_yyyy_xyyyz, tr_yyyy_xyyzz, tr_yyyy_xyz, tr_yyyy_xyzzz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyyyy, tr_yyyy_yyyyz, tr_yyyy_yyyzz, tr_yyyy_yyz, tr_yyyy_yyzzz, tr_yyyy_yzz, tr_yyyy_yzzzz, tr_yyyy_zzz, tr_yyyyy_xxxx, tr_yyyyy_xxxy, tr_yyyyy_xxxz, tr_yyyyy_xxyy, tr_yyyyy_xxyz, tr_yyyyy_xxzz, tr_yyyyy_xyyy, tr_yyyyy_xyyz, tr_yyyyy_xyzz, tr_yyyyy_xzzz, tr_yyyyy_yyyy, tr_yyyyy_yyyz, tr_yyyyy_yyzz, tr_yyyyy_yzzz, tr_yyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyy_xxxx[i] = 6.0 * tr_y_xxxx[i] - 6.0 * tr_yy_xxxxy[i] * tke_0 - 14.0 * tr_yyy_xxxx[i] * tbe_0 + 4.0 * tr_yyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xxxy[i] = 6.0 * tr_y_xxxy[i] + 3.0 * tr_yy_xxx[i] - 6.0 * tr_yy_xxxyy[i] * tke_0 - 14.0 * tr_yyy_xxxy[i] * tbe_0 - 2.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xxxz[i] = 6.0 * tr_y_xxxz[i] - 6.0 * tr_yy_xxxyz[i] * tke_0 - 14.0 * tr_yyy_xxxz[i] * tbe_0 + 4.0 * tr_yyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xxyy[i] = 6.0 * tr_y_xxyy[i] + 6.0 * tr_yy_xxy[i] - 6.0 * tr_yy_xxyyy[i] * tke_0 - 14.0 * tr_yyy_xxyy[i] * tbe_0 - 4.0 * tr_yyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xxyz[i] = 6.0 * tr_y_xxyz[i] + 3.0 * tr_yy_xxz[i] - 6.0 * tr_yy_xxyyz[i] * tke_0 - 14.0 * tr_yyy_xxyz[i] * tbe_0 - 2.0 * tr_yyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xxzz[i] = 6.0 * tr_y_xxzz[i] - 6.0 * tr_yy_xxyzz[i] * tke_0 - 14.0 * tr_yyy_xxzz[i] * tbe_0 + 4.0 * tr_yyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xyyy[i] = 6.0 * tr_y_xyyy[i] + 9.0 * tr_yy_xyy[i] - 6.0 * tr_yy_xyyyy[i] * tke_0 - 14.0 * tr_yyy_xyyy[i] * tbe_0 - 6.0 * tr_yyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xyyz[i] = 6.0 * tr_y_xyyz[i] + 6.0 * tr_yy_xyz[i] - 6.0 * tr_yy_xyyyz[i] * tke_0 - 14.0 * tr_yyy_xyyz[i] * tbe_0 - 4.0 * tr_yyyy_xyz[i] * tbe_0 + 4.0 * tr_yyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xyzz[i] = 6.0 * tr_y_xyzz[i] + 3.0 * tr_yy_xzz[i] - 6.0 * tr_yy_xyyzz[i] * tke_0 - 14.0 * tr_yyy_xyzz[i] * tbe_0 - 2.0 * tr_yyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xzzz[i] = 6.0 * tr_y_xzzz[i] - 6.0 * tr_yy_xyzzz[i] * tke_0 - 14.0 * tr_yyy_xzzz[i] * tbe_0 + 4.0 * tr_yyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_yyyy[i] = 6.0 * tr_y_yyyy[i] + 12.0 * tr_yy_yyy[i] - 6.0 * tr_yy_yyyyy[i] * tke_0 - 14.0 * tr_yyy_yyyy[i] * tbe_0 - 8.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_yyyz[i] = 6.0 * tr_y_yyyz[i] + 9.0 * tr_yy_yyz[i] - 6.0 * tr_yy_yyyyz[i] * tke_0 - 14.0 * tr_yyy_yyyz[i] * tbe_0 - 6.0 * tr_yyyy_yyz[i] * tbe_0 + 4.0 * tr_yyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_yyzz[i] = 6.0 * tr_y_yyzz[i] + 6.0 * tr_yy_yzz[i] - 6.0 * tr_yy_yyyzz[i] * tke_0 - 14.0 * tr_yyy_yyzz[i] * tbe_0 - 4.0 * tr_yyyy_yzz[i] * tbe_0 + 4.0 * tr_yyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_yzzz[i] = 6.0 * tr_y_yzzz[i] + 3.0 * tr_yy_zzz[i] - 6.0 * tr_yy_yyzzz[i] * tke_0 - 14.0 * tr_yyy_yzzz[i] * tbe_0 - 2.0 * tr_yyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_zzzz[i] = 6.0 * tr_y_zzzz[i] - 6.0 * tr_yy_yzzzz[i] * tke_0 - 14.0 * tr_yyy_zzzz[i] * tbe_0 + 4.0 * tr_yyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 705-720 components of targeted buffer : FG

    auto tr_y_0_y_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 705);

    auto tr_y_0_y_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 706);

    auto tr_y_0_y_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 707);

    auto tr_y_0_y_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 708);

    auto tr_y_0_y_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 709);

    auto tr_y_0_y_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 710);

    auto tr_y_0_y_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 711);

    auto tr_y_0_y_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 712);

    auto tr_y_0_y_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 713);

    auto tr_y_0_y_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 714);

    auto tr_y_0_y_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 715);

    auto tr_y_0_y_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 716);

    auto tr_y_0_y_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 717);

    auto tr_y_0_y_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 718);

    auto tr_y_0_y_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 719);

    #pragma omp simd aligned(tr_y_0_y_yyz_xxxx, tr_y_0_y_yyz_xxxy, tr_y_0_y_yyz_xxxz, tr_y_0_y_yyz_xxyy, tr_y_0_y_yyz_xxyz, tr_y_0_y_yyz_xxzz, tr_y_0_y_yyz_xyyy, tr_y_0_y_yyz_xyyz, tr_y_0_y_yyz_xyzz, tr_y_0_y_yyz_xzzz, tr_y_0_y_yyz_yyyy, tr_y_0_y_yyz_yyyz, tr_y_0_y_yyz_yyzz, tr_y_0_y_yyz_yzzz, tr_y_0_y_yyz_zzzz, tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xzzz, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yzzz, tr_yyyyz_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxy, tr_yyyz_xxxyy, tr_yyyz_xxxyz, tr_yyyz_xxy, tr_yyyz_xxyyy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyyyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyyyy, tr_yyyz_yyyyz, tr_yyyz_yyyzz, tr_yyyz_yyz, tr_yyyz_yyzzz, tr_yyyz_yzz, tr_yyyz_yzzzz, tr_yyyz_zzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyz_xxxx[i] = 2.0 * tr_z_xxxx[i] - 4.0 * tr_yz_xxxxy[i] * tke_0 - 10.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xxxy[i] = 2.0 * tr_z_xxxy[i] + 2.0 * tr_yz_xxx[i] - 4.0 * tr_yz_xxxyy[i] * tke_0 - 10.0 * tr_yyz_xxxy[i] * tbe_0 - 2.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xxxz[i] = 2.0 * tr_z_xxxz[i] - 4.0 * tr_yz_xxxyz[i] * tke_0 - 10.0 * tr_yyz_xxxz[i] * tbe_0 + 4.0 * tr_yyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xxyy[i] = 2.0 * tr_z_xxyy[i] + 4.0 * tr_yz_xxy[i] - 4.0 * tr_yz_xxyyy[i] * tke_0 - 10.0 * tr_yyz_xxyy[i] * tbe_0 - 4.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xxyz[i] = 2.0 * tr_z_xxyz[i] + 2.0 * tr_yz_xxz[i] - 4.0 * tr_yz_xxyyz[i] * tke_0 - 10.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xxzz[i] = 2.0 * tr_z_xxzz[i] - 4.0 * tr_yz_xxyzz[i] * tke_0 - 10.0 * tr_yyz_xxzz[i] * tbe_0 + 4.0 * tr_yyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xyyy[i] = 2.0 * tr_z_xyyy[i] + 6.0 * tr_yz_xyy[i] - 4.0 * tr_yz_xyyyy[i] * tke_0 - 10.0 * tr_yyz_xyyy[i] * tbe_0 - 6.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xyyz[i] = 2.0 * tr_z_xyyz[i] + 4.0 * tr_yz_xyz[i] - 4.0 * tr_yz_xyyyz[i] * tke_0 - 10.0 * tr_yyz_xyyz[i] * tbe_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xyzz[i] = 2.0 * tr_z_xyzz[i] + 2.0 * tr_yz_xzz[i] - 4.0 * tr_yz_xyyzz[i] * tke_0 - 10.0 * tr_yyz_xyzz[i] * tbe_0 - 2.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xzzz[i] = 2.0 * tr_z_xzzz[i] - 4.0 * tr_yz_xyzzz[i] * tke_0 - 10.0 * tr_yyz_xzzz[i] * tbe_0 + 4.0 * tr_yyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_yyyy[i] = 2.0 * tr_z_yyyy[i] + 8.0 * tr_yz_yyy[i] - 4.0 * tr_yz_yyyyy[i] * tke_0 - 10.0 * tr_yyz_yyyy[i] * tbe_0 - 8.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_yyyz[i] = 2.0 * tr_z_yyyz[i] + 6.0 * tr_yz_yyz[i] - 4.0 * tr_yz_yyyyz[i] * tke_0 - 10.0 * tr_yyz_yyyz[i] * tbe_0 - 6.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_yyzz[i] = 2.0 * tr_z_yyzz[i] + 4.0 * tr_yz_yzz[i] - 4.0 * tr_yz_yyyzz[i] * tke_0 - 10.0 * tr_yyz_yyzz[i] * tbe_0 - 4.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_yzzz[i] = 2.0 * tr_z_yzzz[i] + 2.0 * tr_yz_zzz[i] - 4.0 * tr_yz_yyzzz[i] * tke_0 - 10.0 * tr_yyz_yzzz[i] * tbe_0 - 2.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_zzzz[i] = 2.0 * tr_z_zzzz[i] - 4.0 * tr_yz_yzzzz[i] * tke_0 - 10.0 * tr_yyz_zzzz[i] * tbe_0 + 4.0 * tr_yyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 720-735 components of targeted buffer : FG

    auto tr_y_0_y_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 720);

    auto tr_y_0_y_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 721);

    auto tr_y_0_y_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 722);

    auto tr_y_0_y_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 723);

    auto tr_y_0_y_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 724);

    auto tr_y_0_y_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 725);

    auto tr_y_0_y_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 726);

    auto tr_y_0_y_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 727);

    auto tr_y_0_y_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 728);

    auto tr_y_0_y_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 729);

    auto tr_y_0_y_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 730);

    auto tr_y_0_y_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 731);

    auto tr_y_0_y_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 732);

    auto tr_y_0_y_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 733);

    auto tr_y_0_y_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 734);

    #pragma omp simd aligned(tr_y_0_y_yzz_xxxx, tr_y_0_y_yzz_xxxy, tr_y_0_y_yzz_xxxz, tr_y_0_y_yzz_xxyy, tr_y_0_y_yzz_xxyz, tr_y_0_y_yzz_xxzz, tr_y_0_y_yzz_xyyy, tr_y_0_y_yzz_xyyz, tr_y_0_y_yzz_xyzz, tr_y_0_y_yzz_xzzz, tr_y_0_y_yzz_yyyy, tr_y_0_y_yzz_yyyz, tr_y_0_y_yzz_yyzz, tr_y_0_y_yzz_yzzz, tr_y_0_y_yzz_zzzz, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xzzz, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yzzz, tr_yyyzz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxy, tr_yyzz_xxxyy, tr_yyzz_xxxyz, tr_yyzz_xxy, tr_yyzz_xxyyy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyyyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyyyy, tr_yyzz_yyyyz, tr_yyzz_yyyzz, tr_yyzz_yyz, tr_yyzz_yyzzz, tr_yyzz_yzz, tr_yyzz_yzzzz, tr_yyzz_zzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yzz_xxxx[i] = -2.0 * tr_zz_xxxxy[i] * tke_0 - 6.0 * tr_yzz_xxxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xxxy[i] = tr_zz_xxx[i] - 2.0 * tr_zz_xxxyy[i] * tke_0 - 6.0 * tr_yzz_xxxy[i] * tbe_0 - 2.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xxxz[i] = -2.0 * tr_zz_xxxyz[i] * tke_0 - 6.0 * tr_yzz_xxxz[i] * tbe_0 + 4.0 * tr_yyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xxyy[i] = 2.0 * tr_zz_xxy[i] - 2.0 * tr_zz_xxyyy[i] * tke_0 - 6.0 * tr_yzz_xxyy[i] * tbe_0 - 4.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xxyz[i] = tr_zz_xxz[i] - 2.0 * tr_zz_xxyyz[i] * tke_0 - 6.0 * tr_yzz_xxyz[i] * tbe_0 - 2.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xxzz[i] = -2.0 * tr_zz_xxyzz[i] * tke_0 - 6.0 * tr_yzz_xxzz[i] * tbe_0 + 4.0 * tr_yyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xyyy[i] = 3.0 * tr_zz_xyy[i] - 2.0 * tr_zz_xyyyy[i] * tke_0 - 6.0 * tr_yzz_xyyy[i] * tbe_0 - 6.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xyyz[i] = 2.0 * tr_zz_xyz[i] - 2.0 * tr_zz_xyyyz[i] * tke_0 - 6.0 * tr_yzz_xyyz[i] * tbe_0 - 4.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xyzz[i] = tr_zz_xzz[i] - 2.0 * tr_zz_xyyzz[i] * tke_0 - 6.0 * tr_yzz_xyzz[i] * tbe_0 - 2.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xzzz[i] = -2.0 * tr_zz_xyzzz[i] * tke_0 - 6.0 * tr_yzz_xzzz[i] * tbe_0 + 4.0 * tr_yyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_yyyy[i] = 4.0 * tr_zz_yyy[i] - 2.0 * tr_zz_yyyyy[i] * tke_0 - 6.0 * tr_yzz_yyyy[i] * tbe_0 - 8.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_yyyz[i] = 3.0 * tr_zz_yyz[i] - 2.0 * tr_zz_yyyyz[i] * tke_0 - 6.0 * tr_yzz_yyyz[i] * tbe_0 - 6.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_yyzz[i] = 2.0 * tr_zz_yzz[i] - 2.0 * tr_zz_yyyzz[i] * tke_0 - 6.0 * tr_yzz_yyzz[i] * tbe_0 - 4.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_yzzz[i] = tr_zz_zzz[i] - 2.0 * tr_zz_yyzzz[i] * tke_0 - 6.0 * tr_yzz_yzzz[i] * tbe_0 - 2.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_zzzz[i] = -2.0 * tr_zz_yzzzz[i] * tke_0 - 6.0 * tr_yzz_zzzz[i] * tbe_0 + 4.0 * tr_yyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 735-750 components of targeted buffer : FG

    auto tr_y_0_y_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 735);

    auto tr_y_0_y_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 736);

    auto tr_y_0_y_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 737);

    auto tr_y_0_y_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 738);

    auto tr_y_0_y_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 739);

    auto tr_y_0_y_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 740);

    auto tr_y_0_y_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 741);

    auto tr_y_0_y_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 742);

    auto tr_y_0_y_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 743);

    auto tr_y_0_y_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 744);

    auto tr_y_0_y_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 745);

    auto tr_y_0_y_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 746);

    auto tr_y_0_y_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 747);

    auto tr_y_0_y_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 748);

    auto tr_y_0_y_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 749);

    #pragma omp simd aligned(tr_y_0_y_zzz_xxxx, tr_y_0_y_zzz_xxxy, tr_y_0_y_zzz_xxxz, tr_y_0_y_zzz_xxyy, tr_y_0_y_zzz_xxyz, tr_y_0_y_zzz_xxzz, tr_y_0_y_zzz_xyyy, tr_y_0_y_zzz_xyyz, tr_y_0_y_zzz_xyzz, tr_y_0_y_zzz_xzzz, tr_y_0_y_zzz_yyyy, tr_y_0_y_zzz_yyyz, tr_y_0_y_zzz_yyzz, tr_y_0_y_zzz_yzzz, tr_y_0_y_zzz_zzzz, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xzzz, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yzzz, tr_yyzzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxy, tr_yzzz_xxxyy, tr_yzzz_xxxyz, tr_yzzz_xxy, tr_yzzz_xxyyy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyyyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyyyy, tr_yzzz_yyyyz, tr_yzzz_yyyzz, tr_yzzz_yyz, tr_yzzz_yyzzz, tr_yzzz_yzz, tr_yzzz_yzzzz, tr_yzzz_zzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zzz_xxxx[i] = -2.0 * tr_zzz_xxxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xxxy[i] = -2.0 * tr_zzz_xxxy[i] * tbe_0 - 2.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xxxz[i] = -2.0 * tr_zzz_xxxz[i] * tbe_0 + 4.0 * tr_yzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xxyy[i] = -2.0 * tr_zzz_xxyy[i] * tbe_0 - 4.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xxyz[i] = -2.0 * tr_zzz_xxyz[i] * tbe_0 - 2.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xxzz[i] = -2.0 * tr_zzz_xxzz[i] * tbe_0 + 4.0 * tr_yzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xyyy[i] = -2.0 * tr_zzz_xyyy[i] * tbe_0 - 6.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xyyz[i] = -2.0 * tr_zzz_xyyz[i] * tbe_0 - 4.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xyzz[i] = -2.0 * tr_zzz_xyzz[i] * tbe_0 - 2.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xzzz[i] = -2.0 * tr_zzz_xzzz[i] * tbe_0 + 4.0 * tr_yzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_yyyy[i] = -2.0 * tr_zzz_yyyy[i] * tbe_0 - 8.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_yyyz[i] = -2.0 * tr_zzz_yyyz[i] * tbe_0 - 6.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_yyzz[i] = -2.0 * tr_zzz_yyzz[i] * tbe_0 - 4.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_yzzz[i] = -2.0 * tr_zzz_yzzz[i] * tbe_0 - 2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_zzzz[i] = -2.0 * tr_zzz_zzzz[i] * tbe_0 + 4.0 * tr_yzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 750-765 components of targeted buffer : FG

    auto tr_y_0_z_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 750);

    auto tr_y_0_z_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 751);

    auto tr_y_0_z_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 752);

    auto tr_y_0_z_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 753);

    auto tr_y_0_z_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 754);

    auto tr_y_0_z_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 755);

    auto tr_y_0_z_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 756);

    auto tr_y_0_z_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 757);

    auto tr_y_0_z_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 758);

    auto tr_y_0_z_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 759);

    auto tr_y_0_z_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 760);

    auto tr_y_0_z_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 761);

    auto tr_y_0_z_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 762);

    auto tr_y_0_z_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 763);

    auto tr_y_0_z_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 764);

    #pragma omp simd aligned(tr_xxxy_xxx, tr_xxxy_xxxxz, tr_xxxy_xxxyz, tr_xxxy_xxxzz, tr_xxxy_xxy, tr_xxxy_xxyyz, tr_xxxy_xxyzz, tr_xxxy_xxz, tr_xxxy_xxzzz, tr_xxxy_xyy, tr_xxxy_xyyyz, tr_xxxy_xyyzz, tr_xxxy_xyz, tr_xxxy_xyzzz, tr_xxxy_xzz, tr_xxxy_xzzzz, tr_xxxy_yyy, tr_xxxy_yyyyz, tr_xxxy_yyyzz, tr_xxxy_yyz, tr_xxxy_yyzzz, tr_xxxy_yzz, tr_xxxy_yzzzz, tr_xxxy_zzz, tr_xxxy_zzzzz, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_y_0_z_xxx_xxxx, tr_y_0_z_xxx_xxxy, tr_y_0_z_xxx_xxxz, tr_y_0_z_xxx_xxyy, tr_y_0_z_xxx_xxyz, tr_y_0_z_xxx_xxzz, tr_y_0_z_xxx_xyyy, tr_y_0_z_xxx_xyyz, tr_y_0_z_xxx_xyzz, tr_y_0_z_xxx_xzzz, tr_y_0_z_xxx_yyyy, tr_y_0_z_xxx_yyyz, tr_y_0_z_xxx_yyzz, tr_y_0_z_xxx_yzzz, tr_y_0_z_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxx_xxxx[i] = 4.0 * tr_xxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xxxy[i] = 4.0 * tr_xxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xxxz[i] = -2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xxyy[i] = 4.0 * tr_xxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xxyz[i] = -2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xxzz[i] = -4.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xyyy[i] = 4.0 * tr_xxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xyyz[i] = -2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xyzz[i] = -4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xzzz[i] = -6.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_yyyy[i] = 4.0 * tr_xxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_yyyz[i] = -2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_yyzz[i] = -4.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_yzzz[i] = -6.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_zzzz[i] = -8.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 765-780 components of targeted buffer : FG

    auto tr_y_0_z_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 765);

    auto tr_y_0_z_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 766);

    auto tr_y_0_z_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 767);

    auto tr_y_0_z_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 768);

    auto tr_y_0_z_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 769);

    auto tr_y_0_z_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 770);

    auto tr_y_0_z_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 771);

    auto tr_y_0_z_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 772);

    auto tr_y_0_z_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 773);

    auto tr_y_0_z_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 774);

    auto tr_y_0_z_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 775);

    auto tr_y_0_z_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 776);

    auto tr_y_0_z_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 777);

    auto tr_y_0_z_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 778);

    auto tr_y_0_z_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 779);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xx_zzzzz, tr_xxyy_xxx, tr_xxyy_xxxxz, tr_xxyy_xxxyz, tr_xxyy_xxxzz, tr_xxyy_xxy, tr_xxyy_xxyyz, tr_xxyy_xxyzz, tr_xxyy_xxz, tr_xxyy_xxzzz, tr_xxyy_xyy, tr_xxyy_xyyyz, tr_xxyy_xyyzz, tr_xxyy_xyz, tr_xxyy_xyzzz, tr_xxyy_xzz, tr_xxyy_xzzzz, tr_xxyy_yyy, tr_xxyy_yyyyz, tr_xxyy_yyyzz, tr_xxyy_yyz, tr_xxyy_yyzzz, tr_xxyy_yzz, tr_xxyy_yzzzz, tr_xxyy_zzz, tr_xxyy_zzzzz, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_y_0_z_xxy_xxxx, tr_y_0_z_xxy_xxxy, tr_y_0_z_xxy_xxxz, tr_y_0_z_xxy_xxyy, tr_y_0_z_xxy_xxyz, tr_y_0_z_xxy_xxzz, tr_y_0_z_xxy_xyyy, tr_y_0_z_xxy_xyyz, tr_y_0_z_xxy_xyzz, tr_y_0_z_xxy_xzzz, tr_y_0_z_xxy_yyyy, tr_y_0_z_xxy_yyyz, tr_y_0_z_xxy_yyzz, tr_y_0_z_xxy_yzzz, tr_y_0_z_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxy_xxxx[i] = -2.0 * tr_xx_xxxxz[i] * tke_0 - 2.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xxxy[i] = -2.0 * tr_xx_xxxyz[i] * tke_0 - 2.0 * tr_xxz_xxxy[i] * tbe_0 + 4.0 * tr_xxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xxxz[i] = tr_xx_xxx[i] - 2.0 * tr_xx_xxxzz[i] * tke_0 - 2.0 * tr_xxz_xxxz[i] * tbe_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xxyy[i] = -2.0 * tr_xx_xxyyz[i] * tke_0 - 2.0 * tr_xxz_xxyy[i] * tbe_0 + 4.0 * tr_xxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xxyz[i] = tr_xx_xxy[i] - 2.0 * tr_xx_xxyzz[i] * tke_0 - 2.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xxzz[i] = 2.0 * tr_xx_xxz[i] - 2.0 * tr_xx_xxzzz[i] * tke_0 - 2.0 * tr_xxz_xxzz[i] * tbe_0 - 4.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xyyy[i] = -2.0 * tr_xx_xyyyz[i] * tke_0 - 2.0 * tr_xxz_xyyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xyyz[i] = tr_xx_xyy[i] - 2.0 * tr_xx_xyyzz[i] * tke_0 - 2.0 * tr_xxz_xyyz[i] * tbe_0 - 2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xyzz[i] = 2.0 * tr_xx_xyz[i] - 2.0 * tr_xx_xyzzz[i] * tke_0 - 2.0 * tr_xxz_xyzz[i] * tbe_0 - 4.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xzzz[i] = 3.0 * tr_xx_xzz[i] - 2.0 * tr_xx_xzzzz[i] * tke_0 - 2.0 * tr_xxz_xzzz[i] * tbe_0 - 6.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_yyyy[i] = -2.0 * tr_xx_yyyyz[i] * tke_0 - 2.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_yyyz[i] = tr_xx_yyy[i] - 2.0 * tr_xx_yyyzz[i] * tke_0 - 2.0 * tr_xxz_yyyz[i] * tbe_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_yyzz[i] = 2.0 * tr_xx_yyz[i] - 2.0 * tr_xx_yyzzz[i] * tke_0 - 2.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_yzzz[i] = 3.0 * tr_xx_yzz[i] - 2.0 * tr_xx_yzzzz[i] * tke_0 - 2.0 * tr_xxz_yzzz[i] * tbe_0 - 6.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_zzzz[i] = 4.0 * tr_xx_zzz[i] - 2.0 * tr_xx_zzzzz[i] * tke_0 - 2.0 * tr_xxz_zzzz[i] * tbe_0 - 8.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 780-795 components of targeted buffer : FG

    auto tr_y_0_z_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 780);

    auto tr_y_0_z_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 781);

    auto tr_y_0_z_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 782);

    auto tr_y_0_z_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 783);

    auto tr_y_0_z_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 784);

    auto tr_y_0_z_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 785);

    auto tr_y_0_z_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 786);

    auto tr_y_0_z_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 787);

    auto tr_y_0_z_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 788);

    auto tr_y_0_z_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 789);

    auto tr_y_0_z_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 790);

    auto tr_y_0_z_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 791);

    auto tr_y_0_z_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 792);

    auto tr_y_0_z_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 793);

    auto tr_y_0_z_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 794);

    #pragma omp simd aligned(tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxz, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxyz_zzzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_y_0_z_xxz_xxxx, tr_y_0_z_xxz_xxxy, tr_y_0_z_xxz_xxxz, tr_y_0_z_xxz_xxyy, tr_y_0_z_xxz_xxyz, tr_y_0_z_xxz_xxzz, tr_y_0_z_xxz_xyyy, tr_y_0_z_xxz_xyyz, tr_y_0_z_xxz_xyzz, tr_y_0_z_xxz_xzzz, tr_y_0_z_xxz_yyyy, tr_y_0_z_xxz_yyyz, tr_y_0_z_xxz_yyzz, tr_y_0_z_xxz_yzzz, tr_y_0_z_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxz_xxxx[i] = -2.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xxxy[i] = -2.0 * tr_xxy_xxxy[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xxxz[i] = -2.0 * tr_xxy_xxxz[i] * tbe_0 - 2.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xxyy[i] = -2.0 * tr_xxy_xxyy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xxyz[i] = -2.0 * tr_xxy_xxyz[i] * tbe_0 - 2.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xxzz[i] = -2.0 * tr_xxy_xxzz[i] * tbe_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xyyy[i] = -2.0 * tr_xxy_xyyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xyyz[i] = -2.0 * tr_xxy_xyyz[i] * tbe_0 - 2.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xyzz[i] = -2.0 * tr_xxy_xyzz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xzzz[i] = -2.0 * tr_xxy_xzzz[i] * tbe_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_yyyy[i] = -2.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_yyyz[i] = -2.0 * tr_xxy_yyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_yyzz[i] = -2.0 * tr_xxy_yyzz[i] * tbe_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_yzzz[i] = -2.0 * tr_xxy_yzzz[i] * tbe_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_zzzz[i] = -2.0 * tr_xxy_zzzz[i] * tbe_0 - 8.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 795-810 components of targeted buffer : FG

    auto tr_y_0_z_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 795);

    auto tr_y_0_z_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 796);

    auto tr_y_0_z_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 797);

    auto tr_y_0_z_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 798);

    auto tr_y_0_z_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 799);

    auto tr_y_0_z_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 800);

    auto tr_y_0_z_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 801);

    auto tr_y_0_z_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 802);

    auto tr_y_0_z_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 803);

    auto tr_y_0_z_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 804);

    auto tr_y_0_z_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 805);

    auto tr_y_0_z_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 806);

    auto tr_y_0_z_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 807);

    auto tr_y_0_z_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 808);

    auto tr_y_0_z_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 809);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyyy_xxx, tr_xyyy_xxxxz, tr_xyyy_xxxyz, tr_xyyy_xxxzz, tr_xyyy_xxy, tr_xyyy_xxyyz, tr_xyyy_xxyzz, tr_xyyy_xxz, tr_xyyy_xxzzz, tr_xyyy_xyy, tr_xyyy_xyyyz, tr_xyyy_xyyzz, tr_xyyy_xyz, tr_xyyy_xyzzz, tr_xyyy_xzz, tr_xyyy_xzzzz, tr_xyyy_yyy, tr_xyyy_yyyyz, tr_xyyy_yyyzz, tr_xyyy_yyz, tr_xyyy_yyzzz, tr_xyyy_yzz, tr_xyyy_yzzzz, tr_xyyy_zzz, tr_xyyy_zzzzz, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_y_0_z_xyy_xxxx, tr_y_0_z_xyy_xxxy, tr_y_0_z_xyy_xxxz, tr_y_0_z_xyy_xxyy, tr_y_0_z_xyy_xxyz, tr_y_0_z_xyy_xxzz, tr_y_0_z_xyy_xyyy, tr_y_0_z_xyy_xyyz, tr_y_0_z_xyy_xyzz, tr_y_0_z_xyy_xzzz, tr_y_0_z_xyy_yyyy, tr_y_0_z_xyy_yyyz, tr_y_0_z_xyy_yyzz, tr_y_0_z_xyy_yzzz, tr_y_0_z_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyy_xxxx[i] = -4.0 * tr_xy_xxxxz[i] * tke_0 - 4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xxxy[i] = -4.0 * tr_xy_xxxyz[i] * tke_0 - 4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xxxz[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xy_xxxzz[i] * tke_0 - 4.0 * tr_xyz_xxxz[i] * tbe_0 - 2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xxyy[i] = -4.0 * tr_xy_xxyyz[i] * tke_0 - 4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xxyz[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xy_xxyzz[i] * tke_0 - 4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xxzz[i] = 4.0 * tr_xy_xxz[i] - 4.0 * tr_xy_xxzzz[i] * tke_0 - 4.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xyyy[i] = -4.0 * tr_xy_xyyyz[i] * tke_0 - 4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xyyz[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xy_xyyzz[i] * tke_0 - 4.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xyzz[i] = 4.0 * tr_xy_xyz[i] - 4.0 * tr_xy_xyzzz[i] * tke_0 - 4.0 * tr_xyz_xyzz[i] * tbe_0 - 4.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xzzz[i] = 6.0 * tr_xy_xzz[i] - 4.0 * tr_xy_xzzzz[i] * tke_0 - 4.0 * tr_xyz_xzzz[i] * tbe_0 - 6.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_yyyy[i] = -4.0 * tr_xy_yyyyz[i] * tke_0 - 4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_yyyz[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xy_yyyzz[i] * tke_0 - 4.0 * tr_xyz_yyyz[i] * tbe_0 - 2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_yyzz[i] = 4.0 * tr_xy_yyz[i] - 4.0 * tr_xy_yyzzz[i] * tke_0 - 4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_yzzz[i] = 6.0 * tr_xy_yzz[i] - 4.0 * tr_xy_yzzzz[i] * tke_0 - 4.0 * tr_xyz_yzzz[i] * tbe_0 - 6.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_zzzz[i] = 8.0 * tr_xy_zzz[i] - 4.0 * tr_xy_zzzzz[i] * tke_0 - 4.0 * tr_xyz_zzzz[i] * tbe_0 - 8.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 810-825 components of targeted buffer : FG

    auto tr_y_0_z_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 810);

    auto tr_y_0_z_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 811);

    auto tr_y_0_z_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 812);

    auto tr_y_0_z_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 813);

    auto tr_y_0_z_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 814);

    auto tr_y_0_z_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 815);

    auto tr_y_0_z_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 816);

    auto tr_y_0_z_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 817);

    auto tr_y_0_z_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 818);

    auto tr_y_0_z_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 819);

    auto tr_y_0_z_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 820);

    auto tr_y_0_z_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 821);

    auto tr_y_0_z_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 822);

    auto tr_y_0_z_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 823);

    auto tr_y_0_z_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 824);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxz, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyyz_zzzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_y_0_z_xyz_xxxx, tr_y_0_z_xyz_xxxy, tr_y_0_z_xyz_xxxz, tr_y_0_z_xyz_xxyy, tr_y_0_z_xyz_xxyz, tr_y_0_z_xyz_xxzz, tr_y_0_z_xyz_xyyy, tr_y_0_z_xyz_xyyz, tr_y_0_z_xyz_xyzz, tr_y_0_z_xyz_xzzz, tr_y_0_z_xyz_yyyy, tr_y_0_z_xyz_yyyz, tr_y_0_z_xyz_yyzz, tr_y_0_z_xyz_yzzz, tr_y_0_z_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyz_xxxx[i] = tr_x_xxxx[i] - 2.0 * tr_xz_xxxxz[i] * tke_0 - 2.0 * tr_xzz_xxxx[i] * tbe_0 - 2.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xxxy[i] = tr_x_xxxy[i] - 2.0 * tr_xz_xxxyz[i] * tke_0 - 2.0 * tr_xzz_xxxy[i] * tbe_0 - 2.0 * tr_xyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xxxz[i] = tr_x_xxxz[i] + tr_xz_xxx[i] - 2.0 * tr_xz_xxxzz[i] * tke_0 - 2.0 * tr_xzz_xxxz[i] * tbe_0 - 2.0 * tr_xyy_xxxz[i] * tbe_0 - 2.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xxyy[i] = tr_x_xxyy[i] - 2.0 * tr_xz_xxyyz[i] * tke_0 - 2.0 * tr_xzz_xxyy[i] * tbe_0 - 2.0 * tr_xyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xxyz[i] = tr_x_xxyz[i] + tr_xz_xxy[i] - 2.0 * tr_xz_xxyzz[i] * tke_0 - 2.0 * tr_xzz_xxyz[i] * tbe_0 - 2.0 * tr_xyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xxzz[i] = tr_x_xxzz[i] + 2.0 * tr_xz_xxz[i] - 2.0 * tr_xz_xxzzz[i] * tke_0 - 2.0 * tr_xzz_xxzz[i] * tbe_0 - 2.0 * tr_xyy_xxzz[i] * tbe_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xyyy[i] = tr_x_xyyy[i] - 2.0 * tr_xz_xyyyz[i] * tke_0 - 2.0 * tr_xzz_xyyy[i] * tbe_0 - 2.0 * tr_xyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xyyz[i] = tr_x_xyyz[i] + tr_xz_xyy[i] - 2.0 * tr_xz_xyyzz[i] * tke_0 - 2.0 * tr_xzz_xyyz[i] * tbe_0 - 2.0 * tr_xyy_xyyz[i] * tbe_0 - 2.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xyzz[i] = tr_x_xyzz[i] + 2.0 * tr_xz_xyz[i] - 2.0 * tr_xz_xyzzz[i] * tke_0 - 2.0 * tr_xzz_xyzz[i] * tbe_0 - 2.0 * tr_xyy_xyzz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xzzz[i] = tr_x_xzzz[i] + 3.0 * tr_xz_xzz[i] - 2.0 * tr_xz_xzzzz[i] * tke_0 - 2.0 * tr_xzz_xzzz[i] * tbe_0 - 2.0 * tr_xyy_xzzz[i] * tbe_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_yyyy[i] = tr_x_yyyy[i] - 2.0 * tr_xz_yyyyz[i] * tke_0 - 2.0 * tr_xzz_yyyy[i] * tbe_0 - 2.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_yyyz[i] = tr_x_yyyz[i] + tr_xz_yyy[i] - 2.0 * tr_xz_yyyzz[i] * tke_0 - 2.0 * tr_xzz_yyyz[i] * tbe_0 - 2.0 * tr_xyy_yyyz[i] * tbe_0 - 2.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_yyzz[i] = tr_x_yyzz[i] + 2.0 * tr_xz_yyz[i] - 2.0 * tr_xz_yyzzz[i] * tke_0 - 2.0 * tr_xzz_yyzz[i] * tbe_0 - 2.0 * tr_xyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_yzzz[i] = tr_x_yzzz[i] + 3.0 * tr_xz_yzz[i] - 2.0 * tr_xz_yzzzz[i] * tke_0 - 2.0 * tr_xzz_yzzz[i] * tbe_0 - 2.0 * tr_xyy_yzzz[i] * tbe_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_zzzz[i] = tr_x_zzzz[i] + 4.0 * tr_xz_zzz[i] - 2.0 * tr_xz_zzzzz[i] * tke_0 - 2.0 * tr_xzz_zzzz[i] * tbe_0 - 2.0 * tr_xyy_zzzz[i] * tbe_0 - 8.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 825-840 components of targeted buffer : FG

    auto tr_y_0_z_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 825);

    auto tr_y_0_z_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 826);

    auto tr_y_0_z_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 827);

    auto tr_y_0_z_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 828);

    auto tr_y_0_z_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 829);

    auto tr_y_0_z_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 830);

    auto tr_y_0_z_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 831);

    auto tr_y_0_z_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 832);

    auto tr_y_0_z_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 833);

    auto tr_y_0_z_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 834);

    auto tr_y_0_z_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 835);

    auto tr_y_0_z_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 836);

    auto tr_y_0_z_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 837);

    auto tr_y_0_z_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 838);

    auto tr_y_0_z_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 839);

    #pragma omp simd aligned(tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxz, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xyzz_zzzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_y_0_z_xzz_xxxx, tr_y_0_z_xzz_xxxy, tr_y_0_z_xzz_xxxz, tr_y_0_z_xzz_xxyy, tr_y_0_z_xzz_xxyz, tr_y_0_z_xzz_xxzz, tr_y_0_z_xzz_xyyy, tr_y_0_z_xzz_xyyz, tr_y_0_z_xzz_xyzz, tr_y_0_z_xzz_xzzz, tr_y_0_z_xzz_yyyy, tr_y_0_z_xzz_yyyz, tr_y_0_z_xzz_yyzz, tr_y_0_z_xzz_yzzz, tr_y_0_z_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xzz_xxxx[i] = -4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xxxy[i] = -4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xxxz[i] = -4.0 * tr_xyz_xxxz[i] * tbe_0 - 2.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xxyy[i] = -4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xxyz[i] = -4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xxzz[i] = -4.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xyyy[i] = -4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xyyz[i] = -4.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xyzz[i] = -4.0 * tr_xyz_xyzz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xzzz[i] = -4.0 * tr_xyz_xzzz[i] * tbe_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_yyyy[i] = -4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_yyyz[i] = -4.0 * tr_xyz_yyyz[i] * tbe_0 - 2.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_yyzz[i] = -4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_yzzz[i] = -4.0 * tr_xyz_yzzz[i] * tbe_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_zzzz[i] = -4.0 * tr_xyz_zzzz[i] * tbe_0 - 8.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 840-855 components of targeted buffer : FG

    auto tr_y_0_z_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 840);

    auto tr_y_0_z_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 841);

    auto tr_y_0_z_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 842);

    auto tr_y_0_z_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 843);

    auto tr_y_0_z_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 844);

    auto tr_y_0_z_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 845);

    auto tr_y_0_z_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 846);

    auto tr_y_0_z_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 847);

    auto tr_y_0_z_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 848);

    auto tr_y_0_z_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 849);

    auto tr_y_0_z_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 850);

    auto tr_y_0_z_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 851);

    auto tr_y_0_z_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 852);

    auto tr_y_0_z_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 853);

    auto tr_y_0_z_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 854);

    #pragma omp simd aligned(tr_y_0_z_yyy_xxxx, tr_y_0_z_yyy_xxxy, tr_y_0_z_yyy_xxxz, tr_y_0_z_yyy_xxyy, tr_y_0_z_yyy_xxyz, tr_y_0_z_yyy_xxzz, tr_y_0_z_yyy_xyyy, tr_y_0_z_yyy_xyyz, tr_y_0_z_yyy_xyzz, tr_y_0_z_yyy_xzzz, tr_y_0_z_yyy_yyyy, tr_y_0_z_yyy_yyyz, tr_y_0_z_yyy_yyzz, tr_y_0_z_yyy_yzzz, tr_y_0_z_yyy_zzzz, tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yy_zzzzz, tr_yyyy_xxx, tr_yyyy_xxxxz, tr_yyyy_xxxyz, tr_yyyy_xxxzz, tr_yyyy_xxy, tr_yyyy_xxyyz, tr_yyyy_xxyzz, tr_yyyy_xxz, tr_yyyy_xxzzz, tr_yyyy_xyy, tr_yyyy_xyyyz, tr_yyyy_xyyzz, tr_yyyy_xyz, tr_yyyy_xyzzz, tr_yyyy_xzz, tr_yyyy_xzzzz, tr_yyyy_yyy, tr_yyyy_yyyyz, tr_yyyy_yyyzz, tr_yyyy_yyz, tr_yyyy_yyzzz, tr_yyyy_yzz, tr_yyyy_yzzzz, tr_yyyy_zzz, tr_yyyy_zzzzz, tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xzzz, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yzzz, tr_yyyyz_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyy_xxxx[i] = -6.0 * tr_yy_xxxxz[i] * tke_0 - 6.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_yyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xxxy[i] = -6.0 * tr_yy_xxxyz[i] * tke_0 - 6.0 * tr_yyz_xxxy[i] * tbe_0 + 4.0 * tr_yyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xxxz[i] = 3.0 * tr_yy_xxx[i] - 6.0 * tr_yy_xxxzz[i] * tke_0 - 6.0 * tr_yyz_xxxz[i] * tbe_0 - 2.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xxyy[i] = -6.0 * tr_yy_xxyyz[i] * tke_0 - 6.0 * tr_yyz_xxyy[i] * tbe_0 + 4.0 * tr_yyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xxyz[i] = 3.0 * tr_yy_xxy[i] - 6.0 * tr_yy_xxyzz[i] * tke_0 - 6.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_yyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xxzz[i] = 6.0 * tr_yy_xxz[i] - 6.0 * tr_yy_xxzzz[i] * tke_0 - 6.0 * tr_yyz_xxzz[i] * tbe_0 - 4.0 * tr_yyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xyyy[i] = -6.0 * tr_yy_xyyyz[i] * tke_0 - 6.0 * tr_yyz_xyyy[i] * tbe_0 + 4.0 * tr_yyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xyyz[i] = 3.0 * tr_yy_xyy[i] - 6.0 * tr_yy_xyyzz[i] * tke_0 - 6.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_yyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xyzz[i] = 6.0 * tr_yy_xyz[i] - 6.0 * tr_yy_xyzzz[i] * tke_0 - 6.0 * tr_yyz_xyzz[i] * tbe_0 - 4.0 * tr_yyyy_xyz[i] * tbe_0 + 4.0 * tr_yyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xzzz[i] = 9.0 * tr_yy_xzz[i] - 6.0 * tr_yy_xzzzz[i] * tke_0 - 6.0 * tr_yyz_xzzz[i] * tbe_0 - 6.0 * tr_yyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_yyyy[i] = -6.0 * tr_yy_yyyyz[i] * tke_0 - 6.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_yyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_yyyz[i] = 3.0 * tr_yy_yyy[i] - 6.0 * tr_yy_yyyzz[i] * tke_0 - 6.0 * tr_yyz_yyyz[i] * tbe_0 - 2.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_yyzz[i] = 6.0 * tr_yy_yyz[i] - 6.0 * tr_yy_yyzzz[i] * tke_0 - 6.0 * tr_yyz_yyzz[i] * tbe_0 - 4.0 * tr_yyyy_yyz[i] * tbe_0 + 4.0 * tr_yyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_yzzz[i] = 9.0 * tr_yy_yzz[i] - 6.0 * tr_yy_yzzzz[i] * tke_0 - 6.0 * tr_yyz_yzzz[i] * tbe_0 - 6.0 * tr_yyyy_yzz[i] * tbe_0 + 4.0 * tr_yyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_zzzz[i] = 12.0 * tr_yy_zzz[i] - 6.0 * tr_yy_zzzzz[i] * tke_0 - 6.0 * tr_yyz_zzzz[i] * tbe_0 - 8.0 * tr_yyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 855-870 components of targeted buffer : FG

    auto tr_y_0_z_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 855);

    auto tr_y_0_z_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 856);

    auto tr_y_0_z_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 857);

    auto tr_y_0_z_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 858);

    auto tr_y_0_z_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 859);

    auto tr_y_0_z_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 860);

    auto tr_y_0_z_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 861);

    auto tr_y_0_z_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 862);

    auto tr_y_0_z_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 863);

    auto tr_y_0_z_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 864);

    auto tr_y_0_z_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 865);

    auto tr_y_0_z_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 866);

    auto tr_y_0_z_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 867);

    auto tr_y_0_z_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 868);

    auto tr_y_0_z_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 869);

    #pragma omp simd aligned(tr_y_0_z_yyz_xxxx, tr_y_0_z_yyz_xxxy, tr_y_0_z_yyz_xxxz, tr_y_0_z_yyz_xxyy, tr_y_0_z_yyz_xxyz, tr_y_0_z_yyz_xxzz, tr_y_0_z_yyz_xyyy, tr_y_0_z_yyz_xyyz, tr_y_0_z_yyz_xyzz, tr_y_0_z_yyz_xzzz, tr_y_0_z_yyz_yyyy, tr_y_0_z_yyz_yyyz, tr_y_0_z_yyz_yyzz, tr_y_0_z_yyz_yzzz, tr_y_0_z_yyz_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxz, tr_yyyz_xxxyz, tr_yyyz_xxxzz, tr_yyyz_xxy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xxzzz, tr_yyyz_xyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_xzzzz, tr_yyyz_yyy, tr_yyyz_yyyyz, tr_yyyz_yyyzz, tr_yyyz_yyz, tr_yyyz_yyzzz, tr_yyyz_yzz, tr_yyyz_yzzzz, tr_yyyz_zzz, tr_yyyz_zzzzz, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xzzz, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yzzz, tr_yyyzz_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyz_xxxx[i] = 2.0 * tr_y_xxxx[i] - 4.0 * tr_yz_xxxxz[i] * tke_0 - 4.0 * tr_yzz_xxxx[i] * tbe_0 - 2.0 * tr_yyy_xxxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xxxy[i] = 2.0 * tr_y_xxxy[i] - 4.0 * tr_yz_xxxyz[i] * tke_0 - 4.0 * tr_yzz_xxxy[i] * tbe_0 - 2.0 * tr_yyy_xxxy[i] * tbe_0 + 4.0 * tr_yyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xxxz[i] = 2.0 * tr_y_xxxz[i] + 2.0 * tr_yz_xxx[i] - 4.0 * tr_yz_xxxzz[i] * tke_0 - 4.0 * tr_yzz_xxxz[i] * tbe_0 - 2.0 * tr_yyy_xxxz[i] * tbe_0 - 2.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xxyy[i] = 2.0 * tr_y_xxyy[i] - 4.0 * tr_yz_xxyyz[i] * tke_0 - 4.0 * tr_yzz_xxyy[i] * tbe_0 - 2.0 * tr_yyy_xxyy[i] * tbe_0 + 4.0 * tr_yyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xxyz[i] = 2.0 * tr_y_xxyz[i] + 2.0 * tr_yz_xxy[i] - 4.0 * tr_yz_xxyzz[i] * tke_0 - 4.0 * tr_yzz_xxyz[i] * tbe_0 - 2.0 * tr_yyy_xxyz[i] * tbe_0 - 2.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xxzz[i] = 2.0 * tr_y_xxzz[i] + 4.0 * tr_yz_xxz[i] - 4.0 * tr_yz_xxzzz[i] * tke_0 - 4.0 * tr_yzz_xxzz[i] * tbe_0 - 2.0 * tr_yyy_xxzz[i] * tbe_0 - 4.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xyyy[i] = 2.0 * tr_y_xyyy[i] - 4.0 * tr_yz_xyyyz[i] * tke_0 - 4.0 * tr_yzz_xyyy[i] * tbe_0 - 2.0 * tr_yyy_xyyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xyyz[i] = 2.0 * tr_y_xyyz[i] + 2.0 * tr_yz_xyy[i] - 4.0 * tr_yz_xyyzz[i] * tke_0 - 4.0 * tr_yzz_xyyz[i] * tbe_0 - 2.0 * tr_yyy_xyyz[i] * tbe_0 - 2.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xyzz[i] = 2.0 * tr_y_xyzz[i] + 4.0 * tr_yz_xyz[i] - 4.0 * tr_yz_xyzzz[i] * tke_0 - 4.0 * tr_yzz_xyzz[i] * tbe_0 - 2.0 * tr_yyy_xyzz[i] * tbe_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xzzz[i] = 2.0 * tr_y_xzzz[i] + 6.0 * tr_yz_xzz[i] - 4.0 * tr_yz_xzzzz[i] * tke_0 - 4.0 * tr_yzz_xzzz[i] * tbe_0 - 2.0 * tr_yyy_xzzz[i] * tbe_0 - 6.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_yyyy[i] = 2.0 * tr_y_yyyy[i] - 4.0 * tr_yz_yyyyz[i] * tke_0 - 4.0 * tr_yzz_yyyy[i] * tbe_0 - 2.0 * tr_yyy_yyyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_yyyz[i] = 2.0 * tr_y_yyyz[i] + 2.0 * tr_yz_yyy[i] - 4.0 * tr_yz_yyyzz[i] * tke_0 - 4.0 * tr_yzz_yyyz[i] * tbe_0 - 2.0 * tr_yyy_yyyz[i] * tbe_0 - 2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_yyzz[i] = 2.0 * tr_y_yyzz[i] + 4.0 * tr_yz_yyz[i] - 4.0 * tr_yz_yyzzz[i] * tke_0 - 4.0 * tr_yzz_yyzz[i] * tbe_0 - 2.0 * tr_yyy_yyzz[i] * tbe_0 - 4.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_yzzz[i] = 2.0 * tr_y_yzzz[i] + 6.0 * tr_yz_yzz[i] - 4.0 * tr_yz_yzzzz[i] * tke_0 - 4.0 * tr_yzz_yzzz[i] * tbe_0 - 2.0 * tr_yyy_yzzz[i] * tbe_0 - 6.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_zzzz[i] = 2.0 * tr_y_zzzz[i] + 8.0 * tr_yz_zzz[i] - 4.0 * tr_yz_zzzzz[i] * tke_0 - 4.0 * tr_yzz_zzzz[i] * tbe_0 - 2.0 * tr_yyy_zzzz[i] * tbe_0 - 8.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 870-885 components of targeted buffer : FG

    auto tr_y_0_z_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 870);

    auto tr_y_0_z_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 871);

    auto tr_y_0_z_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 872);

    auto tr_y_0_z_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 873);

    auto tr_y_0_z_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 874);

    auto tr_y_0_z_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 875);

    auto tr_y_0_z_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 876);

    auto tr_y_0_z_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 877);

    auto tr_y_0_z_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 878);

    auto tr_y_0_z_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 879);

    auto tr_y_0_z_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 880);

    auto tr_y_0_z_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 881);

    auto tr_y_0_z_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 882);

    auto tr_y_0_z_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 883);

    auto tr_y_0_z_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 884);

    #pragma omp simd aligned(tr_y_0_z_yzz_xxxx, tr_y_0_z_yzz_xxxy, tr_y_0_z_yzz_xxxz, tr_y_0_z_yzz_xxyy, tr_y_0_z_yzz_xxyz, tr_y_0_z_yzz_xxzz, tr_y_0_z_yzz_xyyy, tr_y_0_z_yzz_xyyz, tr_y_0_z_yzz_xyzz, tr_y_0_z_yzz_xzzz, tr_y_0_z_yzz_yyyy, tr_y_0_z_yzz_yyyz, tr_y_0_z_yzz_yyzz, tr_y_0_z_yzz_yzzz, tr_y_0_z_yzz_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxz, tr_yyzz_xxxyz, tr_yyzz_xxxzz, tr_yyzz_xxy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xxzzz, tr_yyzz_xyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_xzzzz, tr_yyzz_yyy, tr_yyzz_yyyyz, tr_yyzz_yyyzz, tr_yyzz_yyz, tr_yyzz_yyzzz, tr_yyzz_yzz, tr_yyzz_yzzzz, tr_yyzz_zzz, tr_yyzz_zzzzz, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xzzz, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yzzz, tr_yyzzz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yzz_xxxx[i] = 2.0 * tr_z_xxxx[i] - 2.0 * tr_zz_xxxxz[i] * tke_0 - 2.0 * tr_zzz_xxxx[i] * tbe_0 - 4.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xxxy[i] = 2.0 * tr_z_xxxy[i] - 2.0 * tr_zz_xxxyz[i] * tke_0 - 2.0 * tr_zzz_xxxy[i] * tbe_0 - 4.0 * tr_yyz_xxxy[i] * tbe_0 + 4.0 * tr_yyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xxxz[i] = 2.0 * tr_z_xxxz[i] + tr_zz_xxx[i] - 2.0 * tr_zz_xxxzz[i] * tke_0 - 2.0 * tr_zzz_xxxz[i] * tbe_0 - 4.0 * tr_yyz_xxxz[i] * tbe_0 - 2.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xxyy[i] = 2.0 * tr_z_xxyy[i] - 2.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_zzz_xxyy[i] * tbe_0 - 4.0 * tr_yyz_xxyy[i] * tbe_0 + 4.0 * tr_yyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xxyz[i] = 2.0 * tr_z_xxyz[i] + tr_zz_xxy[i] - 2.0 * tr_zz_xxyzz[i] * tke_0 - 2.0 * tr_zzz_xxyz[i] * tbe_0 - 4.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xxzz[i] = 2.0 * tr_z_xxzz[i] + 2.0 * tr_zz_xxz[i] - 2.0 * tr_zz_xxzzz[i] * tke_0 - 2.0 * tr_zzz_xxzz[i] * tbe_0 - 4.0 * tr_yyz_xxzz[i] * tbe_0 - 4.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xyyy[i] = 2.0 * tr_z_xyyy[i] - 2.0 * tr_zz_xyyyz[i] * tke_0 - 2.0 * tr_zzz_xyyy[i] * tbe_0 - 4.0 * tr_yyz_xyyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xyyz[i] = 2.0 * tr_z_xyyz[i] + tr_zz_xyy[i] - 2.0 * tr_zz_xyyzz[i] * tke_0 - 2.0 * tr_zzz_xyyz[i] * tbe_0 - 4.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xyzz[i] = 2.0 * tr_z_xyzz[i] + 2.0 * tr_zz_xyz[i] - 2.0 * tr_zz_xyzzz[i] * tke_0 - 2.0 * tr_zzz_xyzz[i] * tbe_0 - 4.0 * tr_yyz_xyzz[i] * tbe_0 - 4.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xzzz[i] = 2.0 * tr_z_xzzz[i] + 3.0 * tr_zz_xzz[i] - 2.0 * tr_zz_xzzzz[i] * tke_0 - 2.0 * tr_zzz_xzzz[i] * tbe_0 - 4.0 * tr_yyz_xzzz[i] * tbe_0 - 6.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_yyyy[i] = 2.0 * tr_z_yyyy[i] - 2.0 * tr_zz_yyyyz[i] * tke_0 - 2.0 * tr_zzz_yyyy[i] * tbe_0 - 4.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_yyyz[i] = 2.0 * tr_z_yyyz[i] + tr_zz_yyy[i] - 2.0 * tr_zz_yyyzz[i] * tke_0 - 2.0 * tr_zzz_yyyz[i] * tbe_0 - 4.0 * tr_yyz_yyyz[i] * tbe_0 - 2.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_yyzz[i] = 2.0 * tr_z_yyzz[i] + 2.0 * tr_zz_yyz[i] - 2.0 * tr_zz_yyzzz[i] * tke_0 - 2.0 * tr_zzz_yyzz[i] * tbe_0 - 4.0 * tr_yyz_yyzz[i] * tbe_0 - 4.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_yzzz[i] = 2.0 * tr_z_yzzz[i] + 3.0 * tr_zz_yzz[i] - 2.0 * tr_zz_yzzzz[i] * tke_0 - 2.0 * tr_zzz_yzzz[i] * tbe_0 - 4.0 * tr_yyz_yzzz[i] * tbe_0 - 6.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_zzzz[i] = 2.0 * tr_z_zzzz[i] + 4.0 * tr_zz_zzz[i] - 2.0 * tr_zz_zzzzz[i] * tke_0 - 2.0 * tr_zzz_zzzz[i] * tbe_0 - 4.0 * tr_yyz_zzzz[i] * tbe_0 - 8.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 885-900 components of targeted buffer : FG

    auto tr_y_0_z_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 885);

    auto tr_y_0_z_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 886);

    auto tr_y_0_z_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 887);

    auto tr_y_0_z_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 888);

    auto tr_y_0_z_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 889);

    auto tr_y_0_z_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 890);

    auto tr_y_0_z_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 891);

    auto tr_y_0_z_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 892);

    auto tr_y_0_z_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 893);

    auto tr_y_0_z_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 894);

    auto tr_y_0_z_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 895);

    auto tr_y_0_z_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 896);

    auto tr_y_0_z_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 897);

    auto tr_y_0_z_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 898);

    auto tr_y_0_z_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 899);

    #pragma omp simd aligned(tr_y_0_z_zzz_xxxx, tr_y_0_z_zzz_xxxy, tr_y_0_z_zzz_xxxz, tr_y_0_z_zzz_xxyy, tr_y_0_z_zzz_xxyz, tr_y_0_z_zzz_xxzz, tr_y_0_z_zzz_xyyy, tr_y_0_z_zzz_xyyz, tr_y_0_z_zzz_xyzz, tr_y_0_z_zzz_xzzz, tr_y_0_z_zzz_yyyy, tr_y_0_z_zzz_yyyz, tr_y_0_z_zzz_yyzz, tr_y_0_z_zzz_yzzz, tr_y_0_z_zzz_zzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxz, tr_yzzz_xxxyz, tr_yzzz_xxxzz, tr_yzzz_xxy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xxzzz, tr_yzzz_xyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_xzzzz, tr_yzzz_yyy, tr_yzzz_yyyyz, tr_yzzz_yyyzz, tr_yzzz_yyz, tr_yzzz_yyzzz, tr_yzzz_yzz, tr_yzzz_yzzzz, tr_yzzz_zzz, tr_yzzz_zzzzz, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xzzz, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yzzz, tr_yzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zzz_xxxx[i] = -6.0 * tr_yzz_xxxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xxxy[i] = -6.0 * tr_yzz_xxxy[i] * tbe_0 + 4.0 * tr_yzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xxxz[i] = -6.0 * tr_yzz_xxxz[i] * tbe_0 - 2.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xxyy[i] = -6.0 * tr_yzz_xxyy[i] * tbe_0 + 4.0 * tr_yzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xxyz[i] = -6.0 * tr_yzz_xxyz[i] * tbe_0 - 2.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xxzz[i] = -6.0 * tr_yzz_xxzz[i] * tbe_0 - 4.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xyyy[i] = -6.0 * tr_yzz_xyyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xyyz[i] = -6.0 * tr_yzz_xyyz[i] * tbe_0 - 2.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xyzz[i] = -6.0 * tr_yzz_xyzz[i] * tbe_0 - 4.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xzzz[i] = -6.0 * tr_yzz_xzzz[i] * tbe_0 - 6.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_yyyy[i] = -6.0 * tr_yzz_yyyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_yyyz[i] = -6.0 * tr_yzz_yyyz[i] * tbe_0 - 2.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_yyzz[i] = -6.0 * tr_yzz_yyzz[i] * tbe_0 - 4.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_yzzz[i] = -6.0 * tr_yzz_yzzz[i] * tbe_0 - 6.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_zzzz[i] = -6.0 * tr_yzz_zzzz[i] * tbe_0 - 8.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 900-915 components of targeted buffer : FG

    auto tr_z_0_x_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 900);

    auto tr_z_0_x_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 901);

    auto tr_z_0_x_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 902);

    auto tr_z_0_x_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 903);

    auto tr_z_0_x_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 904);

    auto tr_z_0_x_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 905);

    auto tr_z_0_x_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 906);

    auto tr_z_0_x_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 907);

    auto tr_z_0_x_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 908);

    auto tr_z_0_x_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 909);

    auto tr_z_0_x_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 910);

    auto tr_z_0_x_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 911);

    auto tr_z_0_x_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 912);

    auto tr_z_0_x_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 913);

    auto tr_z_0_x_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 914);

    #pragma omp simd aligned(tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xzzz, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yzzz, tr_xxxxz_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxx, tr_xxxz_xxxxy, tr_xxxz_xxxxz, tr_xxxz_xxxyy, tr_xxxz_xxxyz, tr_xxxz_xxxzz, tr_xxxz_xxy, tr_xxxz_xxyyy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xxzzz, tr_xxxz_xyy, tr_xxxz_xyyyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_xzzzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_z_0_x_xxx_xxxx, tr_z_0_x_xxx_xxxy, tr_z_0_x_xxx_xxxz, tr_z_0_x_xxx_xxyy, tr_z_0_x_xxx_xxyz, tr_z_0_x_xxx_xxzz, tr_z_0_x_xxx_xyyy, tr_z_0_x_xxx_xyyz, tr_z_0_x_xxx_xyzz, tr_z_0_x_xxx_xzzz, tr_z_0_x_xxx_yyyy, tr_z_0_x_xxx_yyyz, tr_z_0_x_xxx_yyzz, tr_z_0_x_xxx_yzzz, tr_z_0_x_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxx_xxxx[i] = -6.0 * tr_xxz_xxxx[i] * tbe_0 - 8.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xxxy[i] = -6.0 * tr_xxz_xxxy[i] * tbe_0 - 6.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xxxz[i] = -6.0 * tr_xxz_xxxz[i] * tbe_0 - 6.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xxyy[i] = -6.0 * tr_xxz_xxyy[i] * tbe_0 - 4.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xxyz[i] = -6.0 * tr_xxz_xxyz[i] * tbe_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xxzz[i] = -6.0 * tr_xxz_xxzz[i] * tbe_0 - 4.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xyyy[i] = -6.0 * tr_xxz_xyyy[i] * tbe_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xyyz[i] = -6.0 * tr_xxz_xyyz[i] * tbe_0 - 2.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xyzz[i] = -6.0 * tr_xxz_xyzz[i] * tbe_0 - 2.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xzzz[i] = -6.0 * tr_xxz_xzzz[i] * tbe_0 - 2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_yyyy[i] = -6.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_yyyz[i] = -6.0 * tr_xxz_yyyz[i] * tbe_0 + 4.0 * tr_xxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_yyzz[i] = -6.0 * tr_xxz_yyzz[i] * tbe_0 + 4.0 * tr_xxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_yzzz[i] = -6.0 * tr_xxz_yzzz[i] * tbe_0 + 4.0 * tr_xxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_zzzz[i] = -6.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 915-930 components of targeted buffer : FG

    auto tr_z_0_x_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 915);

    auto tr_z_0_x_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 916);

    auto tr_z_0_x_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 917);

    auto tr_z_0_x_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 918);

    auto tr_z_0_x_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 919);

    auto tr_z_0_x_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 920);

    auto tr_z_0_x_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 921);

    auto tr_z_0_x_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 922);

    auto tr_z_0_x_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 923);

    auto tr_z_0_x_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 924);

    auto tr_z_0_x_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 925);

    auto tr_z_0_x_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 926);

    auto tr_z_0_x_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 927);

    auto tr_z_0_x_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 928);

    auto tr_z_0_x_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 929);

    #pragma omp simd aligned(tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxx, tr_xxyz_xxxxy, tr_xxyz_xxxxz, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_z_0_x_xxy_xxxx, tr_z_0_x_xxy_xxxy, tr_z_0_x_xxy_xxxz, tr_z_0_x_xxy_xxyy, tr_z_0_x_xxy_xxyz, tr_z_0_x_xxy_xxzz, tr_z_0_x_xxy_xyyy, tr_z_0_x_xxy_xyyz, tr_z_0_x_xxy_xyzz, tr_z_0_x_xxy_xzzz, tr_z_0_x_xxy_yyyy, tr_z_0_x_xxy_yyyz, tr_z_0_x_xxy_yyzz, tr_z_0_x_xxy_yzzz, tr_z_0_x_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxy_xxxx[i] = -4.0 * tr_xyz_xxxx[i] * tbe_0 - 8.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xxxy[i] = -4.0 * tr_xyz_xxxy[i] * tbe_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xxxz[i] = -4.0 * tr_xyz_xxxz[i] * tbe_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xxyy[i] = -4.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xxyz[i] = -4.0 * tr_xyz_xxyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xxzz[i] = -4.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xyyy[i] = -4.0 * tr_xyz_xyyy[i] * tbe_0 - 2.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xyyz[i] = -4.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xyzz[i] = -4.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xzzz[i] = -4.0 * tr_xyz_xzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_yyyy[i] = -4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_yyyz[i] = -4.0 * tr_xyz_yyyz[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_yyzz[i] = -4.0 * tr_xyz_yyzz[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_yzzz[i] = -4.0 * tr_xyz_yzzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_zzzz[i] = -4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 930-945 components of targeted buffer : FG

    auto tr_z_0_x_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 930);

    auto tr_z_0_x_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 931);

    auto tr_z_0_x_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 932);

    auto tr_z_0_x_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 933);

    auto tr_z_0_x_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 934);

    auto tr_z_0_x_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 935);

    auto tr_z_0_x_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 936);

    auto tr_z_0_x_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 937);

    auto tr_z_0_x_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 938);

    auto tr_z_0_x_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 939);

    auto tr_z_0_x_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 940);

    auto tr_z_0_x_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 941);

    auto tr_z_0_x_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 942);

    auto tr_z_0_x_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 943);

    auto tr_z_0_x_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 944);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xzzz, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yzzz, tr_xxxzz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxx, tr_xxzz_xxxxy, tr_xxzz_xxxxz, tr_xxzz_xxxyy, tr_xxzz_xxxyz, tr_xxzz_xxxzz, tr_xxzz_xxy, tr_xxzz_xxyyy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xxzzz, tr_xxzz_xyy, tr_xxzz_xyyyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_xzzzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_z_0_x_xxz_xxxx, tr_z_0_x_xxz_xxxy, tr_z_0_x_xxz_xxxz, tr_z_0_x_xxz_xxyy, tr_z_0_x_xxz_xxyz, tr_z_0_x_xxz_xxzz, tr_z_0_x_xxz_xyyy, tr_z_0_x_xxz_xyyz, tr_z_0_x_xxz_xyzz, tr_z_0_x_xxz_xzzz, tr_z_0_x_xxz_yyyy, tr_z_0_x_xxz_yyyz, tr_z_0_x_xxz_yyzz, tr_z_0_x_xxz_yzzz, tr_z_0_x_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxz_xxxx[i] = 2.0 * tr_x_xxxx[i] - 4.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xx_xxx[i] - 2.0 * tr_xx_xxxxx[i] * tke_0 - 8.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xxxy[i] = 2.0 * tr_x_xxxy[i] - 4.0 * tr_xzz_xxxy[i] * tbe_0 + 3.0 * tr_xx_xxy[i] - 2.0 * tr_xx_xxxxy[i] * tke_0 - 6.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xxxz[i] = 2.0 * tr_x_xxxz[i] - 4.0 * tr_xzz_xxxz[i] * tbe_0 + 3.0 * tr_xx_xxz[i] - 2.0 * tr_xx_xxxxz[i] * tke_0 - 6.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxxz[i] * tbe_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xxyy[i] = 2.0 * tr_x_xxyy[i] - 4.0 * tr_xzz_xxyy[i] * tbe_0 + 2.0 * tr_xx_xyy[i] - 2.0 * tr_xx_xxxyy[i] * tke_0 - 4.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxyy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xxyz[i] = 2.0 * tr_x_xxyz[i] - 4.0 * tr_xzz_xxyz[i] * tbe_0 + 2.0 * tr_xx_xyz[i] - 2.0 * tr_xx_xxxyz[i] * tke_0 - 4.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxyz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xxzz[i] = 2.0 * tr_x_xxzz[i] - 4.0 * tr_xzz_xxzz[i] * tbe_0 + 2.0 * tr_xx_xzz[i] - 2.0 * tr_xx_xxxzz[i] * tke_0 - 4.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xxzz[i] * tbe_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xyyy[i] = 2.0 * tr_x_xyyy[i] - 4.0 * tr_xzz_xyyy[i] * tbe_0 + tr_xx_yyy[i] - 2.0 * tr_xx_xxyyy[i] * tke_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xyyz[i] = 2.0 * tr_x_xyyz[i] - 4.0 * tr_xzz_xyyz[i] * tbe_0 + tr_xx_yyz[i] - 2.0 * tr_xx_xxyyz[i] * tke_0 - 2.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyyz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xyzz[i] = 2.0 * tr_x_xyzz[i] - 4.0 * tr_xzz_xyzz[i] * tbe_0 + tr_xx_yzz[i] - 2.0 * tr_xx_xxyzz[i] * tke_0 - 2.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xyzz[i] * tbe_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xzzz[i] = 2.0 * tr_x_xzzz[i] - 4.0 * tr_xzz_xzzz[i] * tbe_0 + tr_xx_zzz[i] - 2.0 * tr_xx_xxzzz[i] * tke_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xzzz[i] * tbe_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_yyyy[i] = 2.0 * tr_x_yyyy[i] - 4.0 * tr_xzz_yyyy[i] * tbe_0 - 2.0 * tr_xx_xyyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_yyyz[i] = 2.0 * tr_x_yyyz[i] - 4.0 * tr_xzz_yyyz[i] * tbe_0 - 2.0 * tr_xx_xyyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyyz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_yyzz[i] = 2.0 * tr_x_yyzz[i] - 4.0 * tr_xzz_yyzz[i] * tbe_0 - 2.0 * tr_xx_xyyzz[i] * tke_0 + 4.0 * tr_xxzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yyzz[i] * tbe_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_yzzz[i] = 2.0 * tr_x_yzzz[i] - 4.0 * tr_xzz_yzzz[i] * tbe_0 - 2.0 * tr_xx_xyzzz[i] * tke_0 + 4.0 * tr_xxzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yzzz[i] * tbe_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_zzzz[i] = 2.0 * tr_x_zzzz[i] - 4.0 * tr_xzz_zzzz[i] * tbe_0 - 2.0 * tr_xx_xzzzz[i] * tke_0 + 4.0 * tr_xxzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zzzz[i] * tbe_0 + 4.0 * tr_xxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 945-960 components of targeted buffer : FG

    auto tr_z_0_x_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 945);

    auto tr_z_0_x_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 946);

    auto tr_z_0_x_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 947);

    auto tr_z_0_x_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 948);

    auto tr_z_0_x_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 949);

    auto tr_z_0_x_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 950);

    auto tr_z_0_x_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 951);

    auto tr_z_0_x_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 952);

    auto tr_z_0_x_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 953);

    auto tr_z_0_x_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 954);

    auto tr_z_0_x_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 955);

    auto tr_z_0_x_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 956);

    auto tr_z_0_x_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 957);

    auto tr_z_0_x_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 958);

    auto tr_z_0_x_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 959);

    #pragma omp simd aligned(tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxx, tr_xyyz_xxxxy, tr_xyyz_xxxxz, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_z_0_x_xyy_xxxx, tr_z_0_x_xyy_xxxy, tr_z_0_x_xyy_xxxz, tr_z_0_x_xyy_xxyy, tr_z_0_x_xyy_xxyz, tr_z_0_x_xyy_xxzz, tr_z_0_x_xyy_xyyy, tr_z_0_x_xyy_xyyz, tr_z_0_x_xyy_xyzz, tr_z_0_x_xyy_xzzz, tr_z_0_x_xyy_yyyy, tr_z_0_x_xyy_yyyz, tr_z_0_x_xyy_yyzz, tr_z_0_x_xyy_yzzz, tr_z_0_x_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyy_xxxx[i] = -2.0 * tr_yyz_xxxx[i] * tbe_0 - 8.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xxxy[i] = -2.0 * tr_yyz_xxxy[i] * tbe_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xxxz[i] = -2.0 * tr_yyz_xxxz[i] * tbe_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xxyy[i] = -2.0 * tr_yyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xxyz[i] = -2.0 * tr_yyz_xxyz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xxzz[i] = -2.0 * tr_yyz_xxzz[i] * tbe_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xyyy[i] = -2.0 * tr_yyz_xyyy[i] * tbe_0 - 2.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xyyz[i] = -2.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xyzz[i] = -2.0 * tr_yyz_xyzz[i] * tbe_0 - 2.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xzzz[i] = -2.0 * tr_yyz_xzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_yyyy[i] = -2.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_yyyz[i] = -2.0 * tr_yyz_yyyz[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_yyzz[i] = -2.0 * tr_yyz_yyzz[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_yzzz[i] = -2.0 * tr_yyz_yzzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_zzzz[i] = -2.0 * tr_yyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 960-975 components of targeted buffer : FG

    auto tr_z_0_x_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 960);

    auto tr_z_0_x_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 961);

    auto tr_z_0_x_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 962);

    auto tr_z_0_x_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 963);

    auto tr_z_0_x_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 964);

    auto tr_z_0_x_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 965);

    auto tr_z_0_x_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 966);

    auto tr_z_0_x_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 967);

    auto tr_z_0_x_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 968);

    auto tr_z_0_x_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 969);

    auto tr_z_0_x_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 970);

    auto tr_z_0_x_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 971);

    auto tr_z_0_x_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 972);

    auto tr_z_0_x_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 973);

    auto tr_z_0_x_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 974);

    #pragma omp simd aligned(tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyzz_xxx, tr_xyzz_xxxxx, tr_xyzz_xxxxy, tr_xyzz_xxxxz, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_z_0_x_xyz_xxxx, tr_z_0_x_xyz_xxxy, tr_z_0_x_xyz_xxxz, tr_z_0_x_xyz_xxyy, tr_z_0_x_xyz_xxyz, tr_z_0_x_xyz_xxzz, tr_z_0_x_xyz_xyyy, tr_z_0_x_xyz_xyyz, tr_z_0_x_xyz_xyzz, tr_z_0_x_xyz_xzzz, tr_z_0_x_xyz_yyyy, tr_z_0_x_xyz_yyyz, tr_z_0_x_xyz_yyzz, tr_z_0_x_xyz_yzzz, tr_z_0_x_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyz_xxxx[i] = tr_y_xxxx[i] - 2.0 * tr_yzz_xxxx[i] * tbe_0 + 4.0 * tr_xy_xxx[i] - 2.0 * tr_xy_xxxxx[i] * tke_0 - 8.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xxxy[i] = tr_y_xxxy[i] - 2.0 * tr_yzz_xxxy[i] * tbe_0 + 3.0 * tr_xy_xxy[i] - 2.0 * tr_xy_xxxxy[i] * tke_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xxxz[i] = tr_y_xxxz[i] - 2.0 * tr_yzz_xxxz[i] * tbe_0 + 3.0 * tr_xy_xxz[i] - 2.0 * tr_xy_xxxxz[i] * tke_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xxyy[i] = tr_y_xxyy[i] - 2.0 * tr_yzz_xxyy[i] * tbe_0 + 2.0 * tr_xy_xyy[i] - 2.0 * tr_xy_xxxyy[i] * tke_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xxyz[i] = tr_y_xxyz[i] - 2.0 * tr_yzz_xxyz[i] * tbe_0 + 2.0 * tr_xy_xyz[i] - 2.0 * tr_xy_xxxyz[i] * tke_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxyz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xxzz[i] = tr_y_xxzz[i] - 2.0 * tr_yzz_xxzz[i] * tbe_0 + 2.0 * tr_xy_xzz[i] - 2.0 * tr_xy_xxxzz[i] * tke_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xyyy[i] = tr_y_xyyy[i] - 2.0 * tr_yzz_xyyy[i] * tbe_0 + tr_xy_yyy[i] - 2.0 * tr_xy_xxyyy[i] * tke_0 - 2.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xyyz[i] = tr_y_xyyz[i] - 2.0 * tr_yzz_xyyz[i] * tbe_0 + tr_xy_yyz[i] - 2.0 * tr_xy_xxyyz[i] * tke_0 - 2.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xyzz[i] = tr_y_xyzz[i] - 2.0 * tr_yzz_xyzz[i] * tbe_0 + tr_xy_yzz[i] - 2.0 * tr_xy_xxyzz[i] * tke_0 - 2.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xzzz[i] = tr_y_xzzz[i] - 2.0 * tr_yzz_xzzz[i] * tbe_0 + tr_xy_zzz[i] - 2.0 * tr_xy_xxzzz[i] * tke_0 - 2.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xzzz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_yyyy[i] = tr_y_yyyy[i] - 2.0 * tr_yzz_yyyy[i] * tbe_0 - 2.0 * tr_xy_xyyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_yyyz[i] = tr_y_yyyz[i] - 2.0 * tr_yzz_yyyz[i] * tbe_0 - 2.0 * tr_xy_xyyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyyz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_yyzz[i] = tr_y_yyzz[i] - 2.0 * tr_yzz_yyzz[i] * tbe_0 - 2.0 * tr_xy_xyyzz[i] * tke_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyzz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_yzzz[i] = tr_y_yzzz[i] - 2.0 * tr_yzz_yzzz[i] * tbe_0 - 2.0 * tr_xy_xyzzz[i] * tke_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yzzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_zzzz[i] = tr_y_zzzz[i] - 2.0 * tr_yzz_zzzz[i] * tbe_0 - 2.0 * tr_xy_xzzzz[i] * tke_0 + 4.0 * tr_xyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 975-990 components of targeted buffer : FG

    auto tr_z_0_x_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 975);

    auto tr_z_0_x_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 976);

    auto tr_z_0_x_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 977);

    auto tr_z_0_x_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 978);

    auto tr_z_0_x_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 979);

    auto tr_z_0_x_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 980);

    auto tr_z_0_x_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 981);

    auto tr_z_0_x_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 982);

    auto tr_z_0_x_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 983);

    auto tr_z_0_x_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 984);

    auto tr_z_0_x_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 985);

    auto tr_z_0_x_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 986);

    auto tr_z_0_x_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 987);

    auto tr_z_0_x_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 988);

    auto tr_z_0_x_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 989);

    #pragma omp simd aligned(tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xzzz, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yzzz, tr_xxzzz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzzz_xxx, tr_xzzz_xxxxx, tr_xzzz_xxxxy, tr_xzzz_xxxxz, tr_xzzz_xxxyy, tr_xzzz_xxxyz, tr_xzzz_xxxzz, tr_xzzz_xxy, tr_xzzz_xxyyy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xxzzz, tr_xzzz_xyy, tr_xzzz_xyyyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_xzzzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_z_0_x_xzz_xxxx, tr_z_0_x_xzz_xxxy, tr_z_0_x_xzz_xxxz, tr_z_0_x_xzz_xxyy, tr_z_0_x_xzz_xxyz, tr_z_0_x_xzz_xxzz, tr_z_0_x_xzz_xyyy, tr_z_0_x_xzz_xyyz, tr_z_0_x_xzz_xyzz, tr_z_0_x_xzz_xzzz, tr_z_0_x_xzz_yyyy, tr_z_0_x_xzz_yyyz, tr_z_0_x_xzz_yyzz, tr_z_0_x_xzz_yzzz, tr_z_0_x_xzz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xzz_xxxx[i] = 2.0 * tr_z_xxxx[i] - 2.0 * tr_zzz_xxxx[i] * tbe_0 + 8.0 * tr_xz_xxx[i] - 4.0 * tr_xz_xxxxx[i] * tke_0 - 8.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xxxy[i] = 2.0 * tr_z_xxxy[i] - 2.0 * tr_zzz_xxxy[i] * tbe_0 + 6.0 * tr_xz_xxy[i] - 4.0 * tr_xz_xxxxy[i] * tke_0 - 6.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xxxz[i] = 2.0 * tr_z_xxxz[i] - 2.0 * tr_zzz_xxxz[i] * tbe_0 + 6.0 * tr_xz_xxz[i] - 4.0 * tr_xz_xxxxz[i] * tke_0 - 6.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxxz[i] * tbe_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xxyy[i] = 2.0 * tr_z_xxyy[i] - 2.0 * tr_zzz_xxyy[i] * tbe_0 + 4.0 * tr_xz_xyy[i] - 4.0 * tr_xz_xxxyy[i] * tke_0 - 4.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxyy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xxyz[i] = 2.0 * tr_z_xxyz[i] - 2.0 * tr_zzz_xxyz[i] * tbe_0 + 4.0 * tr_xz_xyz[i] - 4.0 * tr_xz_xxxyz[i] * tke_0 - 4.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxyz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xxzz[i] = 2.0 * tr_z_xxzz[i] - 2.0 * tr_zzz_xxzz[i] * tbe_0 + 4.0 * tr_xz_xzz[i] - 4.0 * tr_xz_xxxzz[i] * tke_0 - 4.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xxzz[i] * tbe_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xyyy[i] = 2.0 * tr_z_xyyy[i] - 2.0 * tr_zzz_xyyy[i] * tbe_0 + 2.0 * tr_xz_yyy[i] - 4.0 * tr_xz_xxyyy[i] * tke_0 - 2.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xyyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xyyz[i] = 2.0 * tr_z_xyyz[i] - 2.0 * tr_zzz_xyyz[i] * tbe_0 + 2.0 * tr_xz_yyz[i] - 4.0 * tr_xz_xxyyz[i] * tke_0 - 2.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xyyz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xyzz[i] = 2.0 * tr_z_xyzz[i] - 2.0 * tr_zzz_xyzz[i] * tbe_0 + 2.0 * tr_xz_yzz[i] - 4.0 * tr_xz_xxyzz[i] * tke_0 - 2.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xyzz[i] * tbe_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xzzz[i] = 2.0 * tr_z_xzzz[i] - 2.0 * tr_zzz_xzzz[i] * tbe_0 + 2.0 * tr_xz_zzz[i] - 4.0 * tr_xz_xxzzz[i] * tke_0 - 2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xzzz[i] * tbe_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_yyyy[i] = 2.0 * tr_z_yyyy[i] - 2.0 * tr_zzz_yyyy[i] * tbe_0 - 4.0 * tr_xz_xyyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_yyyz[i] = 2.0 * tr_z_yyyz[i] - 2.0 * tr_zzz_yyyz[i] * tbe_0 - 4.0 * tr_xz_xyyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yyyz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_yyzz[i] = 2.0 * tr_z_yyzz[i] - 2.0 * tr_zzz_yyzz[i] * tbe_0 - 4.0 * tr_xz_xyyzz[i] * tke_0 + 4.0 * tr_xzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yyzz[i] * tbe_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_yzzz[i] = 2.0 * tr_z_yzzz[i] - 2.0 * tr_zzz_yzzz[i] * tbe_0 - 4.0 * tr_xz_xyzzz[i] * tke_0 + 4.0 * tr_xzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yzzz[i] * tbe_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_zzzz[i] = 2.0 * tr_z_zzzz[i] - 2.0 * tr_zzz_zzzz[i] * tbe_0 - 4.0 * tr_xz_xzzzz[i] * tke_0 + 4.0 * tr_xzzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 990-1005 components of targeted buffer : FG

    auto tr_z_0_x_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 990);

    auto tr_z_0_x_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 991);

    auto tr_z_0_x_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 992);

    auto tr_z_0_x_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 993);

    auto tr_z_0_x_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 994);

    auto tr_z_0_x_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 995);

    auto tr_z_0_x_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 996);

    auto tr_z_0_x_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 997);

    auto tr_z_0_x_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 998);

    auto tr_z_0_x_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 999);

    auto tr_z_0_x_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1000);

    auto tr_z_0_x_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1001);

    auto tr_z_0_x_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1002);

    auto tr_z_0_x_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1003);

    auto tr_z_0_x_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1004);

    #pragma omp simd aligned(tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxx, tr_yyyz_xxxxy, tr_yyyz_xxxxz, tr_yyyz_xxxyy, tr_yyyz_xxxyz, tr_yyyz_xxxzz, tr_yyyz_xxy, tr_yyyz_xxyyy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xxzzz, tr_yyyz_xyy, tr_yyyz_xyyyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_xzzzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_z_0_x_yyy_xxxx, tr_z_0_x_yyy_xxxy, tr_z_0_x_yyy_xxxz, tr_z_0_x_yyy_xxyy, tr_z_0_x_yyy_xxyz, tr_z_0_x_yyy_xxzz, tr_z_0_x_yyy_xyyy, tr_z_0_x_yyy_xyyz, tr_z_0_x_yyy_xyzz, tr_z_0_x_yyy_xzzz, tr_z_0_x_yyy_yyyy, tr_z_0_x_yyy_yyyz, tr_z_0_x_yyy_yyzz, tr_z_0_x_yyy_yzzz, tr_z_0_x_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyy_xxxx[i] = -8.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xxxy[i] = -6.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xxxz[i] = -6.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xxyy[i] = -4.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xxyz[i] = -4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xxzz[i] = -4.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xyyy[i] = -2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xyyz[i] = -2.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xyzz[i] = -2.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xzzz[i] = -2.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_yyyy[i] = 4.0 * tr_yyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_yyyz[i] = 4.0 * tr_yyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_yyzz[i] = 4.0 * tr_yyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_yzzz[i] = 4.0 * tr_yyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_zzzz[i] = 4.0 * tr_yyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1005-1020 components of targeted buffer : FG

    auto tr_z_0_x_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1005);

    auto tr_z_0_x_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1006);

    auto tr_z_0_x_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1007);

    auto tr_z_0_x_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1008);

    auto tr_z_0_x_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1009);

    auto tr_z_0_x_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1010);

    auto tr_z_0_x_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1011);

    auto tr_z_0_x_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1012);

    auto tr_z_0_x_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1013);

    auto tr_z_0_x_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1014);

    auto tr_z_0_x_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1015);

    auto tr_z_0_x_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1016);

    auto tr_z_0_x_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1017);

    auto tr_z_0_x_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1018);

    auto tr_z_0_x_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1019);

    #pragma omp simd aligned(tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyzz_xxx, tr_yyzz_xxxxx, tr_yyzz_xxxxy, tr_yyzz_xxxxz, tr_yyzz_xxxyy, tr_yyzz_xxxyz, tr_yyzz_xxxzz, tr_yyzz_xxy, tr_yyzz_xxyyy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xxzzz, tr_yyzz_xyy, tr_yyzz_xyyyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_xzzzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_z_0_x_yyz_xxxx, tr_z_0_x_yyz_xxxy, tr_z_0_x_yyz_xxxz, tr_z_0_x_yyz_xxyy, tr_z_0_x_yyz_xxyz, tr_z_0_x_yyz_xxzz, tr_z_0_x_yyz_xyyy, tr_z_0_x_yyz_xyyz, tr_z_0_x_yyz_xyzz, tr_z_0_x_yyz_xzzz, tr_z_0_x_yyz_yyyy, tr_z_0_x_yyz_yyyz, tr_z_0_x_yyz_yyzz, tr_z_0_x_yyz_yzzz, tr_z_0_x_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyz_xxxx[i] = 4.0 * tr_yy_xxx[i] - 2.0 * tr_yy_xxxxx[i] * tke_0 - 8.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xxxy[i] = 3.0 * tr_yy_xxy[i] - 2.0 * tr_yy_xxxxy[i] * tke_0 - 6.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xxxz[i] = 3.0 * tr_yy_xxz[i] - 2.0 * tr_yy_xxxxz[i] * tke_0 - 6.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xxyy[i] = 2.0 * tr_yy_xyy[i] - 2.0 * tr_yy_xxxyy[i] * tke_0 - 4.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xxyz[i] = 2.0 * tr_yy_xyz[i] - 2.0 * tr_yy_xxxyz[i] * tke_0 - 4.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xxzz[i] = 2.0 * tr_yy_xzz[i] - 2.0 * tr_yy_xxxzz[i] * tke_0 - 4.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xyyy[i] = tr_yy_yyy[i] - 2.0 * tr_yy_xxyyy[i] * tke_0 - 2.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xyyz[i] = tr_yy_yyz[i] - 2.0 * tr_yy_xxyyz[i] * tke_0 - 2.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xyzz[i] = tr_yy_yzz[i] - 2.0 * tr_yy_xxyzz[i] * tke_0 - 2.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xzzz[i] = tr_yy_zzz[i] - 2.0 * tr_yy_xxzzz[i] * tke_0 - 2.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xzzz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_yyyy[i] = -2.0 * tr_yy_xyyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_yyyz[i] = -2.0 * tr_yy_xyyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_yyzz[i] = -2.0 * tr_yy_xyyzz[i] * tke_0 + 4.0 * tr_yyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_yzzz[i] = -2.0 * tr_yy_xyzzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yzzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_zzzz[i] = -2.0 * tr_yy_xzzzz[i] * tke_0 + 4.0 * tr_yyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1020-1035 components of targeted buffer : FG

    auto tr_z_0_x_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1020);

    auto tr_z_0_x_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1021);

    auto tr_z_0_x_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1022);

    auto tr_z_0_x_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1023);

    auto tr_z_0_x_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1024);

    auto tr_z_0_x_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1025);

    auto tr_z_0_x_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1026);

    auto tr_z_0_x_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1027);

    auto tr_z_0_x_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1028);

    auto tr_z_0_x_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1029);

    auto tr_z_0_x_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1030);

    auto tr_z_0_x_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1031);

    auto tr_z_0_x_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1032);

    auto tr_z_0_x_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1033);

    auto tr_z_0_x_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1034);

    #pragma omp simd aligned(tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzzz_xxx, tr_yzzz_xxxxx, tr_yzzz_xxxxy, tr_yzzz_xxxxz, tr_yzzz_xxxyy, tr_yzzz_xxxyz, tr_yzzz_xxxzz, tr_yzzz_xxy, tr_yzzz_xxyyy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xxzzz, tr_yzzz_xyy, tr_yzzz_xyyyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_xzzzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_z_0_x_yzz_xxxx, tr_z_0_x_yzz_xxxy, tr_z_0_x_yzz_xxxz, tr_z_0_x_yzz_xxyy, tr_z_0_x_yzz_xxyz, tr_z_0_x_yzz_xxzz, tr_z_0_x_yzz_xyyy, tr_z_0_x_yzz_xyyz, tr_z_0_x_yzz_xyzz, tr_z_0_x_yzz_xzzz, tr_z_0_x_yzz_yyyy, tr_z_0_x_yzz_yyyz, tr_z_0_x_yzz_yyzz, tr_z_0_x_yzz_yzzz, tr_z_0_x_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yzz_xxxx[i] = 8.0 * tr_yz_xxx[i] - 4.0 * tr_yz_xxxxx[i] * tke_0 - 8.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xxxy[i] = 6.0 * tr_yz_xxy[i] - 4.0 * tr_yz_xxxxy[i] * tke_0 - 6.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xxxz[i] = 6.0 * tr_yz_xxz[i] - 4.0 * tr_yz_xxxxz[i] * tke_0 - 6.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xxyy[i] = 4.0 * tr_yz_xyy[i] - 4.0 * tr_yz_xxxyy[i] * tke_0 - 4.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xxyz[i] = 4.0 * tr_yz_xyz[i] - 4.0 * tr_yz_xxxyz[i] * tke_0 - 4.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxyz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xxzz[i] = 4.0 * tr_yz_xzz[i] - 4.0 * tr_yz_xxxzz[i] * tke_0 - 4.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xyyy[i] = 2.0 * tr_yz_yyy[i] - 4.0 * tr_yz_xxyyy[i] * tke_0 - 2.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xyyz[i] = 2.0 * tr_yz_yyz[i] - 4.0 * tr_yz_xxyyz[i] * tke_0 - 2.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xyzz[i] = 2.0 * tr_yz_yzz[i] - 4.0 * tr_yz_xxyzz[i] * tke_0 - 2.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xzzz[i] = 2.0 * tr_yz_zzz[i] - 4.0 * tr_yz_xxzzz[i] * tke_0 - 2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_yyyy[i] = -4.0 * tr_yz_xyyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_yyyz[i] = -4.0 * tr_yz_xyyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyyz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_yyzz[i] = -4.0 * tr_yz_xyyzz[i] * tke_0 + 4.0 * tr_yzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyzz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_yzzz[i] = -4.0 * tr_yz_xyzzz[i] * tke_0 + 4.0 * tr_yzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yzzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_zzzz[i] = -4.0 * tr_yz_xzzzz[i] * tke_0 + 4.0 * tr_yzzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1035-1050 components of targeted buffer : FG

    auto tr_z_0_x_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1035);

    auto tr_z_0_x_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1036);

    auto tr_z_0_x_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1037);

    auto tr_z_0_x_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1038);

    auto tr_z_0_x_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1039);

    auto tr_z_0_x_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1040);

    auto tr_z_0_x_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1041);

    auto tr_z_0_x_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1042);

    auto tr_z_0_x_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1043);

    auto tr_z_0_x_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1044);

    auto tr_z_0_x_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1045);

    auto tr_z_0_x_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1046);

    auto tr_z_0_x_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1047);

    auto tr_z_0_x_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1048);

    auto tr_z_0_x_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1049);

    #pragma omp simd aligned(tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xzzz, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yzzz, tr_xzzzz_zzzz, tr_z_0_x_zzz_xxxx, tr_z_0_x_zzz_xxxy, tr_z_0_x_zzz_xxxz, tr_z_0_x_zzz_xxyy, tr_z_0_x_zzz_xxyz, tr_z_0_x_zzz_xxzz, tr_z_0_x_zzz_xyyy, tr_z_0_x_zzz_xyyz, tr_z_0_x_zzz_xyzz, tr_z_0_x_zzz_xzzz, tr_z_0_x_zzz_yyyy, tr_z_0_x_zzz_yyyz, tr_z_0_x_zzz_yyzz, tr_z_0_x_zzz_yzzz, tr_z_0_x_zzz_zzzz, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzzz_xxx, tr_zzzz_xxxxx, tr_zzzz_xxxxy, tr_zzzz_xxxxz, tr_zzzz_xxxyy, tr_zzzz_xxxyz, tr_zzzz_xxxzz, tr_zzzz_xxy, tr_zzzz_xxyyy, tr_zzzz_xxyyz, tr_zzzz_xxyzz, tr_zzzz_xxz, tr_zzzz_xxzzz, tr_zzzz_xyy, tr_zzzz_xyyyy, tr_zzzz_xyyyz, tr_zzzz_xyyzz, tr_zzzz_xyz, tr_zzzz_xyzzz, tr_zzzz_xzz, tr_zzzz_xzzzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zzz_xxxx[i] = 12.0 * tr_zz_xxx[i] - 6.0 * tr_zz_xxxxx[i] * tke_0 - 8.0 * tr_zzzz_xxx[i] * tbe_0 + 4.0 * tr_zzzz_xxxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xxxy[i] = 9.0 * tr_zz_xxy[i] - 6.0 * tr_zz_xxxxy[i] * tke_0 - 6.0 * tr_zzzz_xxy[i] * tbe_0 + 4.0 * tr_zzzz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xxxz[i] = 9.0 * tr_zz_xxz[i] - 6.0 * tr_zz_xxxxz[i] * tke_0 - 6.0 * tr_zzzz_xxz[i] * tbe_0 + 4.0 * tr_zzzz_xxxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxxz[i] * tbe_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xxyy[i] = 6.0 * tr_zz_xyy[i] - 6.0 * tr_zz_xxxyy[i] * tke_0 - 4.0 * tr_zzzz_xyy[i] * tbe_0 + 4.0 * tr_zzzz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxyy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xxyz[i] = 6.0 * tr_zz_xyz[i] - 6.0 * tr_zz_xxxyz[i] * tke_0 - 4.0 * tr_zzzz_xyz[i] * tbe_0 + 4.0 * tr_zzzz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxyz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xxzz[i] = 6.0 * tr_zz_xzz[i] - 6.0 * tr_zz_xxxzz[i] * tke_0 - 4.0 * tr_zzzz_xzz[i] * tbe_0 + 4.0 * tr_zzzz_xxxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xxzz[i] * tbe_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xyyy[i] = 3.0 * tr_zz_yyy[i] - 6.0 * tr_zz_xxyyy[i] * tke_0 - 2.0 * tr_zzzz_yyy[i] * tbe_0 + 4.0 * tr_zzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xyyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xyyz[i] = 3.0 * tr_zz_yyz[i] - 6.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_zzzz_yyz[i] * tbe_0 + 4.0 * tr_zzzz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xyyz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xyzz[i] = 3.0 * tr_zz_yzz[i] - 6.0 * tr_zz_xxyzz[i] * tke_0 - 2.0 * tr_zzzz_yzz[i] * tbe_0 + 4.0 * tr_zzzz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xyzz[i] * tbe_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xzzz[i] = 3.0 * tr_zz_zzz[i] - 6.0 * tr_zz_xxzzz[i] * tke_0 - 2.0 * tr_zzzz_zzz[i] * tbe_0 + 4.0 * tr_zzzz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xzzz[i] * tbe_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_yyyy[i] = -6.0 * tr_zz_xyyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yyyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_yyyz[i] = -6.0 * tr_zz_xyyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yyyz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_yyzz[i] = -6.0 * tr_zz_xyyzz[i] * tke_0 + 4.0 * tr_zzzz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yyzz[i] * tbe_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_yzzz[i] = -6.0 * tr_zz_xyzzz[i] * tke_0 + 4.0 * tr_zzzz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yzzz[i] * tbe_0 + 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_zzzz[i] = -6.0 * tr_zz_xzzzz[i] * tke_0 + 4.0 * tr_zzzz_xzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_zzzz[i] * tbe_0 + 4.0 * tr_xzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1050-1065 components of targeted buffer : FG

    auto tr_z_0_y_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 1050);

    auto tr_z_0_y_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 1051);

    auto tr_z_0_y_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 1052);

    auto tr_z_0_y_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 1053);

    auto tr_z_0_y_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 1054);

    auto tr_z_0_y_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 1055);

    auto tr_z_0_y_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 1056);

    auto tr_z_0_y_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 1057);

    auto tr_z_0_y_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 1058);

    auto tr_z_0_y_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 1059);

    auto tr_z_0_y_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 1060);

    auto tr_z_0_y_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 1061);

    auto tr_z_0_y_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 1062);

    auto tr_z_0_y_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 1063);

    auto tr_z_0_y_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 1064);

    #pragma omp simd aligned(tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xzzz, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yzzz, tr_xxxyz_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxy, tr_xxxz_xxxyy, tr_xxxz_xxxyz, tr_xxxz_xxy, tr_xxxz_xxyyy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyyyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyyyy, tr_xxxz_yyyyz, tr_xxxz_yyyzz, tr_xxxz_yyz, tr_xxxz_yyzzz, tr_xxxz_yzz, tr_xxxz_yzzzz, tr_xxxz_zzz, tr_z_0_y_xxx_xxxx, tr_z_0_y_xxx_xxxy, tr_z_0_y_xxx_xxxz, tr_z_0_y_xxx_xxyy, tr_z_0_y_xxx_xxyz, tr_z_0_y_xxx_xxzz, tr_z_0_y_xxx_xyyy, tr_z_0_y_xxx_xyyz, tr_z_0_y_xxx_xyzz, tr_z_0_y_xxx_xzzz, tr_z_0_y_xxx_yyyy, tr_z_0_y_xxx_yyyz, tr_z_0_y_xxx_yyzz, tr_z_0_y_xxx_yzzz, tr_z_0_y_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxx_xxxx[i] = 4.0 * tr_xxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xxxy[i] = -2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xxxz[i] = 4.0 * tr_xxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xxyy[i] = -4.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xxyz[i] = -2.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xxzz[i] = 4.0 * tr_xxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xyyy[i] = -6.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xyyz[i] = -4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xyzz[i] = -2.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xzzz[i] = 4.0 * tr_xxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_yyyy[i] = -8.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_yyyz[i] = -6.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_yyzz[i] = -4.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_yzzz[i] = -2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_zzzz[i] = 4.0 * tr_xxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1065-1080 components of targeted buffer : FG

    auto tr_z_0_y_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 1065);

    auto tr_z_0_y_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 1066);

    auto tr_z_0_y_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 1067);

    auto tr_z_0_y_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 1068);

    auto tr_z_0_y_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 1069);

    auto tr_z_0_y_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 1070);

    auto tr_z_0_y_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 1071);

    auto tr_z_0_y_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 1072);

    auto tr_z_0_y_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 1073);

    auto tr_z_0_y_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 1074);

    auto tr_z_0_y_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1075);

    auto tr_z_0_y_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1076);

    auto tr_z_0_y_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1077);

    auto tr_z_0_y_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1078);

    auto tr_z_0_y_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1079);

    #pragma omp simd aligned(tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xzzz, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yzzz, tr_xxyyz_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxy, tr_xxyz_xxxyy, tr_xxyz_xxxyz, tr_xxyz_xxy, tr_xxyz_xxyyy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyyyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyyyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_z_0_y_xxy_xxxx, tr_z_0_y_xxy_xxxy, tr_z_0_y_xxy_xxxz, tr_z_0_y_xxy_xxyy, tr_z_0_y_xxy_xxyz, tr_z_0_y_xxy_xxzz, tr_z_0_y_xxy_xyyy, tr_z_0_y_xxy_xyyz, tr_z_0_y_xxy_xyzz, tr_z_0_y_xxy_xzzz, tr_z_0_y_xxy_yyyy, tr_z_0_y_xxy_yyyz, tr_z_0_y_xxy_yyzz, tr_z_0_y_xxy_yzzz, tr_z_0_y_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxy_xxxx[i] = -2.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xxxy[i] = -2.0 * tr_xxz_xxxy[i] * tbe_0 - 2.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xxxz[i] = -2.0 * tr_xxz_xxxz[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xxyy[i] = -2.0 * tr_xxz_xxyy[i] * tbe_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xxyz[i] = -2.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xxzz[i] = -2.0 * tr_xxz_xxzz[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xyyy[i] = -2.0 * tr_xxz_xyyy[i] * tbe_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xyyz[i] = -2.0 * tr_xxz_xyyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xyzz[i] = -2.0 * tr_xxz_xyzz[i] * tbe_0 - 2.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xzzz[i] = -2.0 * tr_xxz_xzzz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_yyyy[i] = -2.0 * tr_xxz_yyyy[i] * tbe_0 - 8.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_yyyz[i] = -2.0 * tr_xxz_yyyz[i] * tbe_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_yyzz[i] = -2.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_yzzz[i] = -2.0 * tr_xxz_yzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_zzzz[i] = -2.0 * tr_xxz_zzzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1080-1095 components of targeted buffer : FG

    auto tr_z_0_y_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1080);

    auto tr_z_0_y_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1081);

    auto tr_z_0_y_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1082);

    auto tr_z_0_y_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1083);

    auto tr_z_0_y_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1084);

    auto tr_z_0_y_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1085);

    auto tr_z_0_y_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1086);

    auto tr_z_0_y_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1087);

    auto tr_z_0_y_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1088);

    auto tr_z_0_y_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1089);

    auto tr_z_0_y_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1090);

    auto tr_z_0_y_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1091);

    auto tr_z_0_y_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1092);

    auto tr_z_0_y_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1093);

    auto tr_z_0_y_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1094);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxy, tr_xxzz_xxxyy, tr_xxzz_xxxyz, tr_xxzz_xxy, tr_xxzz_xxyyy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyyyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyyyy, tr_xxzz_yyyyz, tr_xxzz_yyyzz, tr_xxzz_yyz, tr_xxzz_yyzzz, tr_xxzz_yzz, tr_xxzz_yzzzz, tr_xxzz_zzz, tr_z_0_y_xxz_xxxx, tr_z_0_y_xxz_xxxy, tr_z_0_y_xxz_xxxz, tr_z_0_y_xxz_xxyy, tr_z_0_y_xxz_xxyz, tr_z_0_y_xxz_xxzz, tr_z_0_y_xxz_xyyy, tr_z_0_y_xxz_xyyz, tr_z_0_y_xxz_xyzz, tr_z_0_y_xxz_xzzz, tr_z_0_y_xxz_yyyy, tr_z_0_y_xxz_yyyz, tr_z_0_y_xxz_yyzz, tr_z_0_y_xxz_yzzz, tr_z_0_y_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxz_xxxx[i] = -2.0 * tr_xx_xxxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xxxy[i] = tr_xx_xxx[i] - 2.0 * tr_xx_xxxyy[i] * tke_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xxxz[i] = -2.0 * tr_xx_xxxyz[i] * tke_0 + 4.0 * tr_xxzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xxyy[i] = 2.0 * tr_xx_xxy[i] - 2.0 * tr_xx_xxyyy[i] * tke_0 - 4.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xxyz[i] = tr_xx_xxz[i] - 2.0 * tr_xx_xxyyz[i] * tke_0 - 2.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxyz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xxzz[i] = -2.0 * tr_xx_xxyzz[i] * tke_0 + 4.0 * tr_xxzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xxzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xyyy[i] = 3.0 * tr_xx_xyy[i] - 2.0 * tr_xx_xyyyy[i] * tke_0 - 6.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xyyz[i] = 2.0 * tr_xx_xyz[i] - 2.0 * tr_xx_xyyyz[i] * tke_0 - 4.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xyzz[i] = tr_xx_xzz[i] - 2.0 * tr_xx_xyyzz[i] * tke_0 - 2.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xyzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xzzz[i] = -2.0 * tr_xx_xyzzz[i] * tke_0 + 4.0 * tr_xxzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xzzz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_yyyy[i] = 4.0 * tr_xx_yyy[i] - 2.0 * tr_xx_yyyyy[i] * tke_0 - 8.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_yyyz[i] = 3.0 * tr_xx_yyz[i] - 2.0 * tr_xx_yyyyz[i] * tke_0 - 6.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyyz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_yyzz[i] = 2.0 * tr_xx_yzz[i] - 2.0 * tr_xx_yyyzz[i] * tke_0 - 4.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yyzz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_yzzz[i] = tr_xx_zzz[i] - 2.0 * tr_xx_yyzzz[i] * tke_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yzzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_zzzz[i] = -2.0 * tr_xx_yzzzz[i] * tke_0 + 4.0 * tr_xxzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_zzzz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1095-1110 components of targeted buffer : FG

    auto tr_z_0_y_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 1095);

    auto tr_z_0_y_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 1096);

    auto tr_z_0_y_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 1097);

    auto tr_z_0_y_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 1098);

    auto tr_z_0_y_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 1099);

    auto tr_z_0_y_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 1100);

    auto tr_z_0_y_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 1101);

    auto tr_z_0_y_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 1102);

    auto tr_z_0_y_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 1103);

    auto tr_z_0_y_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 1104);

    auto tr_z_0_y_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1105);

    auto tr_z_0_y_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1106);

    auto tr_z_0_y_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1107);

    auto tr_z_0_y_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1108);

    auto tr_z_0_y_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1109);

    #pragma omp simd aligned(tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xzzz, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yzzz, tr_xyyyz_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxy, tr_xyyz_xxxyy, tr_xyyz_xxxyz, tr_xyyz_xxy, tr_xyyz_xxyyy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyyyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyyyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_z_0_y_xyy_xxxx, tr_z_0_y_xyy_xxxy, tr_z_0_y_xyy_xxxz, tr_z_0_y_xyy_xxyy, tr_z_0_y_xyy_xxyz, tr_z_0_y_xyy_xxzz, tr_z_0_y_xyy_xyyy, tr_z_0_y_xyy_xyyz, tr_z_0_y_xyy_xyzz, tr_z_0_y_xyy_xzzz, tr_z_0_y_xyy_yyyy, tr_z_0_y_xyy_yyyz, tr_z_0_y_xyy_yyzz, tr_z_0_y_xyy_yzzz, tr_z_0_y_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyy_xxxx[i] = -4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xxxy[i] = -4.0 * tr_xyz_xxxy[i] * tbe_0 - 2.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xxxz[i] = -4.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xxyy[i] = -4.0 * tr_xyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xxyz[i] = -4.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xxzz[i] = -4.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xyyy[i] = -4.0 * tr_xyz_xyyy[i] * tbe_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xyyz[i] = -4.0 * tr_xyz_xyyz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xyzz[i] = -4.0 * tr_xyz_xyzz[i] * tbe_0 - 2.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xzzz[i] = -4.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_yyyy[i] = -4.0 * tr_xyz_yyyy[i] * tbe_0 - 8.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_yyyz[i] = -4.0 * tr_xyz_yyyz[i] * tbe_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_yyzz[i] = -4.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_yzzz[i] = -4.0 * tr_xyz_yzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_zzzz[i] = -4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1110-1125 components of targeted buffer : FG

    auto tr_z_0_y_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1110);

    auto tr_z_0_y_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1111);

    auto tr_z_0_y_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1112);

    auto tr_z_0_y_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1113);

    auto tr_z_0_y_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1114);

    auto tr_z_0_y_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1115);

    auto tr_z_0_y_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1116);

    auto tr_z_0_y_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1117);

    auto tr_z_0_y_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1118);

    auto tr_z_0_y_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1119);

    auto tr_z_0_y_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1120);

    auto tr_z_0_y_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1121);

    auto tr_z_0_y_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1122);

    auto tr_z_0_y_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1123);

    auto tr_z_0_y_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1124);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxy, tr_xyzz_xxxyy, tr_xyzz_xxxyz, tr_xyzz_xxy, tr_xyzz_xxyyy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyyyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyyyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_z_0_y_xyz_xxxx, tr_z_0_y_xyz_xxxy, tr_z_0_y_xyz_xxxz, tr_z_0_y_xyz_xxyy, tr_z_0_y_xyz_xxyz, tr_z_0_y_xyz_xxzz, tr_z_0_y_xyz_xyyy, tr_z_0_y_xyz_xyyz, tr_z_0_y_xyz_xyzz, tr_z_0_y_xyz_xzzz, tr_z_0_y_xyz_yyyy, tr_z_0_y_xyz_yyyz, tr_z_0_y_xyz_yyzz, tr_z_0_y_xyz_yzzz, tr_z_0_y_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyz_xxxx[i] = tr_x_xxxx[i] - 2.0 * tr_xzz_xxxx[i] * tbe_0 - 2.0 * tr_xy_xxxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xxxy[i] = tr_x_xxxy[i] - 2.0 * tr_xzz_xxxy[i] * tbe_0 + tr_xy_xxx[i] - 2.0 * tr_xy_xxxyy[i] * tke_0 - 2.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xxxz[i] = tr_x_xxxz[i] - 2.0 * tr_xzz_xxxz[i] * tbe_0 - 2.0 * tr_xy_xxxyz[i] * tke_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xxyy[i] = tr_x_xxyy[i] - 2.0 * tr_xzz_xxyy[i] * tbe_0 + 2.0 * tr_xy_xxy[i] - 2.0 * tr_xy_xxyyy[i] * tke_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xxyz[i] = tr_x_xxyz[i] - 2.0 * tr_xzz_xxyz[i] * tbe_0 + tr_xy_xxz[i] - 2.0 * tr_xy_xxyyz[i] * tke_0 - 2.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xxzz[i] = tr_x_xxzz[i] - 2.0 * tr_xzz_xxzz[i] * tbe_0 - 2.0 * tr_xy_xxyzz[i] * tke_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xxzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xyyy[i] = tr_x_xyyy[i] - 2.0 * tr_xzz_xyyy[i] * tbe_0 + 3.0 * tr_xy_xyy[i] - 2.0 * tr_xy_xyyyy[i] * tke_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xyyz[i] = tr_x_xyyz[i] - 2.0 * tr_xzz_xyyz[i] * tbe_0 + 2.0 * tr_xy_xyz[i] - 2.0 * tr_xy_xyyyz[i] * tke_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xyzz[i] = tr_x_xyzz[i] - 2.0 * tr_xzz_xyzz[i] * tbe_0 + tr_xy_xzz[i] - 2.0 * tr_xy_xyyzz[i] * tke_0 - 2.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xyzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xzzz[i] = tr_x_xzzz[i] - 2.0 * tr_xzz_xzzz[i] * tbe_0 - 2.0 * tr_xy_xyzzz[i] * tke_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xzzz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_yyyy[i] = tr_x_yyyy[i] - 2.0 * tr_xzz_yyyy[i] * tbe_0 + 4.0 * tr_xy_yyy[i] - 2.0 * tr_xy_yyyyy[i] * tke_0 - 8.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_yyyz[i] = tr_x_yyyz[i] - 2.0 * tr_xzz_yyyz[i] * tbe_0 + 3.0 * tr_xy_yyz[i] - 2.0 * tr_xy_yyyyz[i] * tke_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_yyzz[i] = tr_x_yyzz[i] - 2.0 * tr_xzz_yyzz[i] * tbe_0 + 2.0 * tr_xy_yzz[i] - 2.0 * tr_xy_yyyzz[i] * tke_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yyzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_yzzz[i] = tr_x_yzzz[i] - 2.0 * tr_xzz_yzzz[i] * tbe_0 + tr_xy_zzz[i] - 2.0 * tr_xy_yyzzz[i] * tke_0 - 2.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yzzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_zzzz[i] = tr_x_zzzz[i] - 2.0 * tr_xzz_zzzz[i] * tbe_0 - 2.0 * tr_xy_yzzzz[i] * tke_0 + 4.0 * tr_xyzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1125-1140 components of targeted buffer : FG

    auto tr_z_0_y_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1125);

    auto tr_z_0_y_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1126);

    auto tr_z_0_y_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1127);

    auto tr_z_0_y_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1128);

    auto tr_z_0_y_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1129);

    auto tr_z_0_y_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1130);

    auto tr_z_0_y_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1131);

    auto tr_z_0_y_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1132);

    auto tr_z_0_y_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1133);

    auto tr_z_0_y_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1134);

    auto tr_z_0_y_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1135);

    auto tr_z_0_y_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1136);

    auto tr_z_0_y_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1137);

    auto tr_z_0_y_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1138);

    auto tr_z_0_y_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1139);

    #pragma omp simd aligned(tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xzzz_xxx, tr_xzzz_xxxxy, tr_xzzz_xxxyy, tr_xzzz_xxxyz, tr_xzzz_xxy, tr_xzzz_xxyyy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyyyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyyyy, tr_xzzz_yyyyz, tr_xzzz_yyyzz, tr_xzzz_yyz, tr_xzzz_yyzzz, tr_xzzz_yzz, tr_xzzz_yzzzz, tr_xzzz_zzz, tr_z_0_y_xzz_xxxx, tr_z_0_y_xzz_xxxy, tr_z_0_y_xzz_xxxz, tr_z_0_y_xzz_xxyy, tr_z_0_y_xzz_xxyz, tr_z_0_y_xzz_xxzz, tr_z_0_y_xzz_xyyy, tr_z_0_y_xzz_xyyz, tr_z_0_y_xzz_xyzz, tr_z_0_y_xzz_xzzz, tr_z_0_y_xzz_yyyy, tr_z_0_y_xzz_yyyz, tr_z_0_y_xzz_yyzz, tr_z_0_y_xzz_yzzz, tr_z_0_y_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xzz_xxxx[i] = -4.0 * tr_xz_xxxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xxxy[i] = 2.0 * tr_xz_xxx[i] - 4.0 * tr_xz_xxxyy[i] * tke_0 - 2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xxxz[i] = -4.0 * tr_xz_xxxyz[i] * tke_0 + 4.0 * tr_xzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xxyy[i] = 4.0 * tr_xz_xxy[i] - 4.0 * tr_xz_xxyyy[i] * tke_0 - 4.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xxyz[i] = 2.0 * tr_xz_xxz[i] - 4.0 * tr_xz_xxyyz[i] * tke_0 - 2.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxyz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xxzz[i] = -4.0 * tr_xz_xxyzz[i] * tke_0 + 4.0 * tr_xzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xxzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xyyy[i] = 6.0 * tr_xz_xyy[i] - 4.0 * tr_xz_xyyyy[i] * tke_0 - 6.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xyyz[i] = 4.0 * tr_xz_xyz[i] - 4.0 * tr_xz_xyyyz[i] * tke_0 - 4.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xyzz[i] = 2.0 * tr_xz_xzz[i] - 4.0 * tr_xz_xyyzz[i] * tke_0 - 2.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xyzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xzzz[i] = -4.0 * tr_xz_xyzzz[i] * tke_0 + 4.0 * tr_xzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xzzz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_yyyy[i] = 8.0 * tr_xz_yyy[i] - 4.0 * tr_xz_yyyyy[i] * tke_0 - 8.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_yyyz[i] = 6.0 * tr_xz_yyz[i] - 4.0 * tr_xz_yyyyz[i] * tke_0 - 6.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzz_yyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyyz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_yyzz[i] = 4.0 * tr_xz_yzz[i] - 4.0 * tr_xz_yyyzz[i] * tke_0 - 4.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yyzz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_yzzz[i] = 2.0 * tr_xz_zzz[i] - 4.0 * tr_xz_yyzzz[i] * tke_0 - 2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzz_yyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yzzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_zzzz[i] = -4.0 * tr_xz_yzzzz[i] * tke_0 + 4.0 * tr_xzzz_yzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zzzz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1140-1155 components of targeted buffer : FG

    auto tr_z_0_y_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 1140);

    auto tr_z_0_y_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 1141);

    auto tr_z_0_y_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 1142);

    auto tr_z_0_y_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 1143);

    auto tr_z_0_y_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 1144);

    auto tr_z_0_y_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 1145);

    auto tr_z_0_y_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 1146);

    auto tr_z_0_y_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 1147);

    auto tr_z_0_y_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 1148);

    auto tr_z_0_y_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 1149);

    auto tr_z_0_y_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1150);

    auto tr_z_0_y_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1151);

    auto tr_z_0_y_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1152);

    auto tr_z_0_y_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1153);

    auto tr_z_0_y_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1154);

    #pragma omp simd aligned(tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xzzz, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yzzz, tr_yyyyz_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxy, tr_yyyz_xxxyy, tr_yyyz_xxxyz, tr_yyyz_xxy, tr_yyyz_xxyyy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyyyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyyyy, tr_yyyz_yyyyz, tr_yyyz_yyyzz, tr_yyyz_yyz, tr_yyyz_yyzzz, tr_yyyz_yzz, tr_yyyz_yzzzz, tr_yyyz_zzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_z_0_y_yyy_xxxx, tr_z_0_y_yyy_xxxy, tr_z_0_y_yyy_xxxz, tr_z_0_y_yyy_xxyy, tr_z_0_y_yyy_xxyz, tr_z_0_y_yyy_xxzz, tr_z_0_y_yyy_xyyy, tr_z_0_y_yyy_xyyz, tr_z_0_y_yyy_xyzz, tr_z_0_y_yyy_xzzz, tr_z_0_y_yyy_yyyy, tr_z_0_y_yyy_yyyz, tr_z_0_y_yyy_yyzz, tr_z_0_y_yyy_yzzz, tr_z_0_y_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyy_xxxx[i] = -6.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xxxy[i] = -6.0 * tr_yyz_xxxy[i] * tbe_0 - 2.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xxxz[i] = -6.0 * tr_yyz_xxxz[i] * tbe_0 + 4.0 * tr_yyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xxyy[i] = -6.0 * tr_yyz_xxyy[i] * tbe_0 - 4.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xxyz[i] = -6.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xxzz[i] = -6.0 * tr_yyz_xxzz[i] * tbe_0 + 4.0 * tr_yyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xyyy[i] = -6.0 * tr_yyz_xyyy[i] * tbe_0 - 6.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xyyz[i] = -6.0 * tr_yyz_xyyz[i] * tbe_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xyzz[i] = -6.0 * tr_yyz_xyzz[i] * tbe_0 - 2.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xzzz[i] = -6.0 * tr_yyz_xzzz[i] * tbe_0 + 4.0 * tr_yyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_yyyy[i] = -6.0 * tr_yyz_yyyy[i] * tbe_0 - 8.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_yyyz[i] = -6.0 * tr_yyz_yyyz[i] * tbe_0 - 6.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_yyzz[i] = -6.0 * tr_yyz_yyzz[i] * tbe_0 - 4.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_yzzz[i] = -6.0 * tr_yyz_yzzz[i] * tbe_0 - 2.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_zzzz[i] = -6.0 * tr_yyz_zzzz[i] * tbe_0 + 4.0 * tr_yyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1155-1170 components of targeted buffer : FG

    auto tr_z_0_y_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1155);

    auto tr_z_0_y_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1156);

    auto tr_z_0_y_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1157);

    auto tr_z_0_y_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1158);

    auto tr_z_0_y_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1159);

    auto tr_z_0_y_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1160);

    auto tr_z_0_y_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1161);

    auto tr_z_0_y_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1162);

    auto tr_z_0_y_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1163);

    auto tr_z_0_y_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1164);

    auto tr_z_0_y_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1165);

    auto tr_z_0_y_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1166);

    auto tr_z_0_y_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1167);

    auto tr_z_0_y_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1168);

    auto tr_z_0_y_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1169);

    #pragma omp simd aligned(tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xzzz, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yzzz, tr_yyyzz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxy, tr_yyzz_xxxyy, tr_yyzz_xxxyz, tr_yyzz_xxy, tr_yyzz_xxyyy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyyyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyyyy, tr_yyzz_yyyyz, tr_yyzz_yyyzz, tr_yyzz_yyz, tr_yyzz_yyzzz, tr_yyzz_yzz, tr_yyzz_yzzzz, tr_yyzz_zzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_z_0_y_yyz_xxxx, tr_z_0_y_yyz_xxxy, tr_z_0_y_yyz_xxxz, tr_z_0_y_yyz_xxyy, tr_z_0_y_yyz_xxyz, tr_z_0_y_yyz_xxzz, tr_z_0_y_yyz_xyyy, tr_z_0_y_yyz_xyyz, tr_z_0_y_yyz_xyzz, tr_z_0_y_yyz_xzzz, tr_z_0_y_yyz_yyyy, tr_z_0_y_yyz_yyyz, tr_z_0_y_yyz_yyzz, tr_z_0_y_yyz_yzzz, tr_z_0_y_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyz_xxxx[i] = 2.0 * tr_y_xxxx[i] - 4.0 * tr_yzz_xxxx[i] * tbe_0 - 2.0 * tr_yy_xxxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xxxy[i] = 2.0 * tr_y_xxxy[i] - 4.0 * tr_yzz_xxxy[i] * tbe_0 + tr_yy_xxx[i] - 2.0 * tr_yy_xxxyy[i] * tke_0 - 2.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xxxz[i] = 2.0 * tr_y_xxxz[i] - 4.0 * tr_yzz_xxxz[i] * tbe_0 - 2.0 * tr_yy_xxxyz[i] * tke_0 + 4.0 * tr_yyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxxz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xxyy[i] = 2.0 * tr_y_xxyy[i] - 4.0 * tr_yzz_xxyy[i] * tbe_0 + 2.0 * tr_yy_xxy[i] - 2.0 * tr_yy_xxyyy[i] * tke_0 - 4.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxyy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xxyz[i] = 2.0 * tr_y_xxyz[i] - 4.0 * tr_yzz_xxyz[i] * tbe_0 + tr_yy_xxz[i] - 2.0 * tr_yy_xxyyz[i] * tke_0 - 2.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxyz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xxzz[i] = 2.0 * tr_y_xxzz[i] - 4.0 * tr_yzz_xxzz[i] * tbe_0 - 2.0 * tr_yy_xxyzz[i] * tke_0 + 4.0 * tr_yyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xxzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xyyy[i] = 2.0 * tr_y_xyyy[i] - 4.0 * tr_yzz_xyyy[i] * tbe_0 + 3.0 * tr_yy_xyy[i] - 2.0 * tr_yy_xyyyy[i] * tke_0 - 6.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xyyy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xyyz[i] = 2.0 * tr_y_xyyz[i] - 4.0 * tr_yzz_xyyz[i] * tbe_0 + 2.0 * tr_yy_xyz[i] - 2.0 * tr_yy_xyyyz[i] * tke_0 - 4.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xyyz[i] * tbe_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xyzz[i] = 2.0 * tr_y_xyzz[i] - 4.0 * tr_yzz_xyzz[i] * tbe_0 + tr_yy_xzz[i] - 2.0 * tr_yy_xyyzz[i] * tke_0 - 2.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xyzz[i] * tbe_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xzzz[i] = 2.0 * tr_y_xzzz[i] - 4.0 * tr_yzz_xzzz[i] * tbe_0 - 2.0 * tr_yy_xyzzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xzzz[i] * tbe_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_yyyy[i] = 2.0 * tr_y_yyyy[i] - 4.0 * tr_yzz_yyyy[i] * tbe_0 + 4.0 * tr_yy_yyy[i] - 2.0 * tr_yy_yyyyy[i] * tke_0 - 8.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yyyy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_yyyz[i] = 2.0 * tr_y_yyyz[i] - 4.0 * tr_yzz_yyyz[i] * tbe_0 + 3.0 * tr_yy_yyz[i] - 2.0 * tr_yy_yyyyz[i] * tke_0 - 6.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yyyz[i] * tbe_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_yyzz[i] = 2.0 * tr_y_yyzz[i] - 4.0 * tr_yzz_yyzz[i] * tbe_0 + 2.0 * tr_yy_yzz[i] - 2.0 * tr_yy_yyyzz[i] * tke_0 - 4.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yyzz[i] * tbe_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_yzzz[i] = 2.0 * tr_y_yzzz[i] - 4.0 * tr_yzz_yzzz[i] * tbe_0 + tr_yy_zzz[i] - 2.0 * tr_yy_yyzzz[i] * tke_0 - 2.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yzzz[i] * tbe_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_zzzz[i] = 2.0 * tr_y_zzzz[i] - 4.0 * tr_yzz_zzzz[i] * tbe_0 - 2.0 * tr_yy_yzzzz[i] * tke_0 + 4.0 * tr_yyzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_zzzz[i] * tbe_0 + 4.0 * tr_yyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1170-1185 components of targeted buffer : FG

    auto tr_z_0_y_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1170);

    auto tr_z_0_y_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1171);

    auto tr_z_0_y_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1172);

    auto tr_z_0_y_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1173);

    auto tr_z_0_y_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1174);

    auto tr_z_0_y_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1175);

    auto tr_z_0_y_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1176);

    auto tr_z_0_y_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1177);

    auto tr_z_0_y_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1178);

    auto tr_z_0_y_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1179);

    auto tr_z_0_y_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1180);

    auto tr_z_0_y_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1181);

    auto tr_z_0_y_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1182);

    auto tr_z_0_y_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1183);

    auto tr_z_0_y_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1184);

    #pragma omp simd aligned(tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xzzz, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yzzz, tr_yyzzz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yzzz_xxx, tr_yzzz_xxxxy, tr_yzzz_xxxyy, tr_yzzz_xxxyz, tr_yzzz_xxy, tr_yzzz_xxyyy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyyyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyyyy, tr_yzzz_yyyyz, tr_yzzz_yyyzz, tr_yzzz_yyz, tr_yzzz_yyzzz, tr_yzzz_yzz, tr_yzzz_yzzzz, tr_yzzz_zzz, tr_z_0_y_yzz_xxxx, tr_z_0_y_yzz_xxxy, tr_z_0_y_yzz_xxxz, tr_z_0_y_yzz_xxyy, tr_z_0_y_yzz_xxyz, tr_z_0_y_yzz_xxzz, tr_z_0_y_yzz_xyyy, tr_z_0_y_yzz_xyyz, tr_z_0_y_yzz_xyzz, tr_z_0_y_yzz_xzzz, tr_z_0_y_yzz_yyyy, tr_z_0_y_yzz_yyyz, tr_z_0_y_yzz_yyzz, tr_z_0_y_yzz_yzzz, tr_z_0_y_yzz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yzz_xxxx[i] = 2.0 * tr_z_xxxx[i] - 2.0 * tr_zzz_xxxx[i] * tbe_0 - 4.0 * tr_yz_xxxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xxxy[i] = 2.0 * tr_z_xxxy[i] - 2.0 * tr_zzz_xxxy[i] * tbe_0 + 2.0 * tr_yz_xxx[i] - 4.0 * tr_yz_xxxyy[i] * tke_0 - 2.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xxxz[i] = 2.0 * tr_z_xxxz[i] - 2.0 * tr_zzz_xxxz[i] * tbe_0 - 4.0 * tr_yz_xxxyz[i] * tke_0 + 4.0 * tr_yzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxxz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xxyy[i] = 2.0 * tr_z_xxyy[i] - 2.0 * tr_zzz_xxyy[i] * tbe_0 + 4.0 * tr_yz_xxy[i] - 4.0 * tr_yz_xxyyy[i] * tke_0 - 4.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxyy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xxyz[i] = 2.0 * tr_z_xxyz[i] - 2.0 * tr_zzz_xxyz[i] * tbe_0 + 2.0 * tr_yz_xxz[i] - 4.0 * tr_yz_xxyyz[i] * tke_0 - 2.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxyz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xxzz[i] = 2.0 * tr_z_xxzz[i] - 2.0 * tr_zzz_xxzz[i] * tbe_0 - 4.0 * tr_yz_xxyzz[i] * tke_0 + 4.0 * tr_yzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xxzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xyyy[i] = 2.0 * tr_z_xyyy[i] - 2.0 * tr_zzz_xyyy[i] * tbe_0 + 6.0 * tr_yz_xyy[i] - 4.0 * tr_yz_xyyyy[i] * tke_0 - 6.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xyyy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xyyz[i] = 2.0 * tr_z_xyyz[i] - 2.0 * tr_zzz_xyyz[i] * tbe_0 + 4.0 * tr_yz_xyz[i] - 4.0 * tr_yz_xyyyz[i] * tke_0 - 4.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xyyz[i] * tbe_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xyzz[i] = 2.0 * tr_z_xyzz[i] - 2.0 * tr_zzz_xyzz[i] * tbe_0 + 2.0 * tr_yz_xzz[i] - 4.0 * tr_yz_xyyzz[i] * tke_0 - 2.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xyzz[i] * tbe_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xzzz[i] = 2.0 * tr_z_xzzz[i] - 2.0 * tr_zzz_xzzz[i] * tbe_0 - 4.0 * tr_yz_xyzzz[i] * tke_0 + 4.0 * tr_yzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xzzz[i] * tbe_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_yyyy[i] = 2.0 * tr_z_yyyy[i] - 2.0 * tr_zzz_yyyy[i] * tbe_0 + 8.0 * tr_yz_yyy[i] - 4.0 * tr_yz_yyyyy[i] * tke_0 - 8.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_yyyz[i] = 2.0 * tr_z_yyyz[i] - 2.0 * tr_zzz_yyyz[i] * tbe_0 + 6.0 * tr_yz_yyz[i] - 4.0 * tr_yz_yyyyz[i] * tke_0 - 6.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzz_yyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yyyz[i] * tbe_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_yyzz[i] = 2.0 * tr_z_yyzz[i] - 2.0 * tr_zzz_yyzz[i] * tbe_0 + 4.0 * tr_yz_yzz[i] - 4.0 * tr_yz_yyyzz[i] * tke_0 - 4.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yyzz[i] * tbe_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_yzzz[i] = 2.0 * tr_z_yzzz[i] - 2.0 * tr_zzz_yzzz[i] * tbe_0 + 2.0 * tr_yz_zzz[i] - 4.0 * tr_yz_yyzzz[i] * tke_0 - 2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzz_yyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yzzz[i] * tbe_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_zzzz[i] = 2.0 * tr_z_zzzz[i] - 2.0 * tr_zzz_zzzz[i] * tbe_0 - 4.0 * tr_yz_yzzzz[i] * tke_0 + 4.0 * tr_yzzz_yzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_zzzz[i] * tbe_0 + 4.0 * tr_yyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1185-1200 components of targeted buffer : FG

    auto tr_z_0_y_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1185);

    auto tr_z_0_y_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1186);

    auto tr_z_0_y_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1187);

    auto tr_z_0_y_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1188);

    auto tr_z_0_y_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1189);

    auto tr_z_0_y_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1190);

    auto tr_z_0_y_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1191);

    auto tr_z_0_y_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1192);

    auto tr_z_0_y_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1193);

    auto tr_z_0_y_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1194);

    auto tr_z_0_y_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1195);

    auto tr_z_0_y_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1196);

    auto tr_z_0_y_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1197);

    auto tr_z_0_y_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1198);

    auto tr_z_0_y_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1199);

    #pragma omp simd aligned(tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xzzz, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yzzz, tr_yzzzz_zzzz, tr_z_0_y_zzz_xxxx, tr_z_0_y_zzz_xxxy, tr_z_0_y_zzz_xxxz, tr_z_0_y_zzz_xxyy, tr_z_0_y_zzz_xxyz, tr_z_0_y_zzz_xxzz, tr_z_0_y_zzz_xyyy, tr_z_0_y_zzz_xyyz, tr_z_0_y_zzz_xyzz, tr_z_0_y_zzz_xzzz, tr_z_0_y_zzz_yyyy, tr_z_0_y_zzz_yyyz, tr_z_0_y_zzz_yyzz, tr_z_0_y_zzz_yzzz, tr_z_0_y_zzz_zzzz, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zzzz_xxx, tr_zzzz_xxxxy, tr_zzzz_xxxyy, tr_zzzz_xxxyz, tr_zzzz_xxy, tr_zzzz_xxyyy, tr_zzzz_xxyyz, tr_zzzz_xxyzz, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyyyy, tr_zzzz_xyyyz, tr_zzzz_xyyzz, tr_zzzz_xyz, tr_zzzz_xyzzz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyyyy, tr_zzzz_yyyyz, tr_zzzz_yyyzz, tr_zzzz_yyz, tr_zzzz_yyzzz, tr_zzzz_yzz, tr_zzzz_yzzzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zzz_xxxx[i] = -6.0 * tr_zz_xxxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xxxy[i] = 3.0 * tr_zz_xxx[i] - 6.0 * tr_zz_xxxyy[i] * tke_0 - 2.0 * tr_zzzz_xxx[i] * tbe_0 + 4.0 * tr_zzzz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xxxz[i] = -6.0 * tr_zz_xxxyz[i] * tke_0 + 4.0 * tr_zzzz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxxz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xxyy[i] = 6.0 * tr_zz_xxy[i] - 6.0 * tr_zz_xxyyy[i] * tke_0 - 4.0 * tr_zzzz_xxy[i] * tbe_0 + 4.0 * tr_zzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxyy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xxyz[i] = 3.0 * tr_zz_xxz[i] - 6.0 * tr_zz_xxyyz[i] * tke_0 - 2.0 * tr_zzzz_xxz[i] * tbe_0 + 4.0 * tr_zzzz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxyz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xxzz[i] = -6.0 * tr_zz_xxyzz[i] * tke_0 + 4.0 * tr_zzzz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xxzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xyyy[i] = 9.0 * tr_zz_xyy[i] - 6.0 * tr_zz_xyyyy[i] * tke_0 - 6.0 * tr_zzzz_xyy[i] * tbe_0 + 4.0 * tr_zzzz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xyyy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xyyz[i] = 6.0 * tr_zz_xyz[i] - 6.0 * tr_zz_xyyyz[i] * tke_0 - 4.0 * tr_zzzz_xyz[i] * tbe_0 + 4.0 * tr_zzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xyyz[i] * tbe_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xyzz[i] = 3.0 * tr_zz_xzz[i] - 6.0 * tr_zz_xyyzz[i] * tke_0 - 2.0 * tr_zzzz_xzz[i] * tbe_0 + 4.0 * tr_zzzz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xyzz[i] * tbe_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xzzz[i] = -6.0 * tr_zz_xyzzz[i] * tke_0 + 4.0 * tr_zzzz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xzzz[i] * tbe_0 + 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_yyyy[i] = 12.0 * tr_zz_yyy[i] - 6.0 * tr_zz_yyyyy[i] * tke_0 - 8.0 * tr_zzzz_yyy[i] * tbe_0 + 4.0 * tr_zzzz_yyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yyyy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_yyyz[i] = 9.0 * tr_zz_yyz[i] - 6.0 * tr_zz_yyyyz[i] * tke_0 - 6.0 * tr_zzzz_yyz[i] * tbe_0 + 4.0 * tr_zzzz_yyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yyyz[i] * tbe_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_yyzz[i] = 6.0 * tr_zz_yzz[i] - 6.0 * tr_zz_yyyzz[i] * tke_0 - 4.0 * tr_zzzz_yzz[i] * tbe_0 + 4.0 * tr_zzzz_yyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yyzz[i] * tbe_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_yzzz[i] = 3.0 * tr_zz_zzz[i] - 6.0 * tr_zz_yyzzz[i] * tke_0 - 2.0 * tr_zzzz_zzz[i] * tbe_0 + 4.0 * tr_zzzz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yzzz[i] * tbe_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_zzzz[i] = -6.0 * tr_zz_yzzzz[i] * tke_0 + 4.0 * tr_zzzz_yzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_zzzz[i] * tbe_0 + 4.0 * tr_yzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1200-1215 components of targeted buffer : FG

    auto tr_z_0_z_xxx_xxxx = pbuffer.data(idx_op_geom_110_fg + 1200);

    auto tr_z_0_z_xxx_xxxy = pbuffer.data(idx_op_geom_110_fg + 1201);

    auto tr_z_0_z_xxx_xxxz = pbuffer.data(idx_op_geom_110_fg + 1202);

    auto tr_z_0_z_xxx_xxyy = pbuffer.data(idx_op_geom_110_fg + 1203);

    auto tr_z_0_z_xxx_xxyz = pbuffer.data(idx_op_geom_110_fg + 1204);

    auto tr_z_0_z_xxx_xxzz = pbuffer.data(idx_op_geom_110_fg + 1205);

    auto tr_z_0_z_xxx_xyyy = pbuffer.data(idx_op_geom_110_fg + 1206);

    auto tr_z_0_z_xxx_xyyz = pbuffer.data(idx_op_geom_110_fg + 1207);

    auto tr_z_0_z_xxx_xyzz = pbuffer.data(idx_op_geom_110_fg + 1208);

    auto tr_z_0_z_xxx_xzzz = pbuffer.data(idx_op_geom_110_fg + 1209);

    auto tr_z_0_z_xxx_yyyy = pbuffer.data(idx_op_geom_110_fg + 1210);

    auto tr_z_0_z_xxx_yyyz = pbuffer.data(idx_op_geom_110_fg + 1211);

    auto tr_z_0_z_xxx_yyzz = pbuffer.data(idx_op_geom_110_fg + 1212);

    auto tr_z_0_z_xxx_yzzz = pbuffer.data(idx_op_geom_110_fg + 1213);

    auto tr_z_0_z_xxx_zzzz = pbuffer.data(idx_op_geom_110_fg + 1214);

    #pragma omp simd aligned(tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, tr_xxxz_xxx, tr_xxxz_xxxxz, tr_xxxz_xxxyz, tr_xxxz_xxxzz, tr_xxxz_xxy, tr_xxxz_xxyyz, tr_xxxz_xxyzz, tr_xxxz_xxz, tr_xxxz_xxzzz, tr_xxxz_xyy, tr_xxxz_xyyyz, tr_xxxz_xyyzz, tr_xxxz_xyz, tr_xxxz_xyzzz, tr_xxxz_xzz, tr_xxxz_xzzzz, tr_xxxz_yyy, tr_xxxz_yyyyz, tr_xxxz_yyyzz, tr_xxxz_yyz, tr_xxxz_yyzzz, tr_xxxz_yzz, tr_xxxz_yzzzz, tr_xxxz_zzz, tr_xxxz_zzzzz, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xzzz, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yzzz, tr_xxxzz_zzzz, tr_z_0_z_xxx_xxxx, tr_z_0_z_xxx_xxxy, tr_z_0_z_xxx_xxxz, tr_z_0_z_xxx_xxyy, tr_z_0_z_xxx_xxyz, tr_z_0_z_xxx_xxzz, tr_z_0_z_xxx_xyyy, tr_z_0_z_xxx_xyyz, tr_z_0_z_xxx_xyzz, tr_z_0_z_xxx_xzzz, tr_z_0_z_xxx_yyyy, tr_z_0_z_xxx_yyyz, tr_z_0_z_xxx_yyzz, tr_z_0_z_xxx_yzzz, tr_z_0_z_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxx_xxxx[i] = -2.0 * tr_xxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xxxy[i] = -2.0 * tr_xxx_xxxy[i] * tbe_0 + 4.0 * tr_xxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xxxz[i] = -2.0 * tr_xxx_xxxz[i] * tbe_0 - 2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xxyy[i] = -2.0 * tr_xxx_xxyy[i] * tbe_0 + 4.0 * tr_xxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xxyz[i] = -2.0 * tr_xxx_xxyz[i] * tbe_0 - 2.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xxzz[i] = -2.0 * tr_xxx_xxzz[i] * tbe_0 - 4.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xyyy[i] = -2.0 * tr_xxx_xyyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xyyz[i] = -2.0 * tr_xxx_xyyz[i] * tbe_0 - 2.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xyzz[i] = -2.0 * tr_xxx_xyzz[i] * tbe_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xzzz[i] = -2.0 * tr_xxx_xzzz[i] * tbe_0 - 6.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_yyyy[i] = -2.0 * tr_xxx_yyyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_yyyz[i] = -2.0 * tr_xxx_yyyz[i] * tbe_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_yyzz[i] = -2.0 * tr_xxx_yyzz[i] * tbe_0 - 4.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_yzzz[i] = -2.0 * tr_xxx_yzzz[i] * tbe_0 - 6.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_zzzz[i] = -2.0 * tr_xxx_zzzz[i] * tbe_0 - 8.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1215-1230 components of targeted buffer : FG

    auto tr_z_0_z_xxy_xxxx = pbuffer.data(idx_op_geom_110_fg + 1215);

    auto tr_z_0_z_xxy_xxxy = pbuffer.data(idx_op_geom_110_fg + 1216);

    auto tr_z_0_z_xxy_xxxz = pbuffer.data(idx_op_geom_110_fg + 1217);

    auto tr_z_0_z_xxy_xxyy = pbuffer.data(idx_op_geom_110_fg + 1218);

    auto tr_z_0_z_xxy_xxyz = pbuffer.data(idx_op_geom_110_fg + 1219);

    auto tr_z_0_z_xxy_xxzz = pbuffer.data(idx_op_geom_110_fg + 1220);

    auto tr_z_0_z_xxy_xyyy = pbuffer.data(idx_op_geom_110_fg + 1221);

    auto tr_z_0_z_xxy_xyyz = pbuffer.data(idx_op_geom_110_fg + 1222);

    auto tr_z_0_z_xxy_xyzz = pbuffer.data(idx_op_geom_110_fg + 1223);

    auto tr_z_0_z_xxy_xzzz = pbuffer.data(idx_op_geom_110_fg + 1224);

    auto tr_z_0_z_xxy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1225);

    auto tr_z_0_z_xxy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1226);

    auto tr_z_0_z_xxy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1227);

    auto tr_z_0_z_xxy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1228);

    auto tr_z_0_z_xxy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1229);

    #pragma omp simd aligned(tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxxxz, tr_xxyz_xxxyz, tr_xxyz_xxxzz, tr_xxyz_xxy, tr_xxyz_xxyyz, tr_xxyz_xxyzz, tr_xxyz_xxz, tr_xxyz_xxzzz, tr_xxyz_xyy, tr_xxyz_xyyyz, tr_xxyz_xyyzz, tr_xxyz_xyz, tr_xxyz_xyzzz, tr_xxyz_xzz, tr_xxyz_xzzzz, tr_xxyz_yyy, tr_xxyz_yyyyz, tr_xxyz_yyyzz, tr_xxyz_yyz, tr_xxyz_yyzzz, tr_xxyz_yzz, tr_xxyz_yzzzz, tr_xxyz_zzz, tr_xxyz_zzzzz, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xzzz, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yzzz, tr_xxyzz_zzzz, tr_z_0_z_xxy_xxxx, tr_z_0_z_xxy_xxxy, tr_z_0_z_xxy_xxxz, tr_z_0_z_xxy_xxyy, tr_z_0_z_xxy_xxyz, tr_z_0_z_xxy_xxzz, tr_z_0_z_xxy_xyyy, tr_z_0_z_xxy_xyyz, tr_z_0_z_xxy_xyzz, tr_z_0_z_xxy_xzzz, tr_z_0_z_xxy_yyyy, tr_z_0_z_xxy_yyyz, tr_z_0_z_xxy_yyzz, tr_z_0_z_xxy_yzzz, tr_z_0_z_xxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxy_xxxx[i] = -2.0 * tr_xxy_xxxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xxxy[i] = -2.0 * tr_xxy_xxxy[i] * tbe_0 + 4.0 * tr_xxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xxxz[i] = -2.0 * tr_xxy_xxxz[i] * tbe_0 - 2.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xxyy[i] = -2.0 * tr_xxy_xxyy[i] * tbe_0 + 4.0 * tr_xxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xxyz[i] = -2.0 * tr_xxy_xxyz[i] * tbe_0 - 2.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xxzz[i] = -2.0 * tr_xxy_xxzz[i] * tbe_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xyyy[i] = -2.0 * tr_xxy_xyyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xyyz[i] = -2.0 * tr_xxy_xyyz[i] * tbe_0 - 2.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xyzz[i] = -2.0 * tr_xxy_xyzz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xzzz[i] = -2.0 * tr_xxy_xzzz[i] * tbe_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_yyyy[i] = -2.0 * tr_xxy_yyyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_yyyz[i] = -2.0 * tr_xxy_yyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_yyzz[i] = -2.0 * tr_xxy_yyzz[i] * tbe_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_yzzz[i] = -2.0 * tr_xxy_yzzz[i] * tbe_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_zzzz[i] = -2.0 * tr_xxy_zzzz[i] * tbe_0 - 8.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1230-1245 components of targeted buffer : FG

    auto tr_z_0_z_xxz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1230);

    auto tr_z_0_z_xxz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1231);

    auto tr_z_0_z_xxz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1232);

    auto tr_z_0_z_xxz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1233);

    auto tr_z_0_z_xxz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1234);

    auto tr_z_0_z_xxz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1235);

    auto tr_z_0_z_xxz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1236);

    auto tr_z_0_z_xxz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1237);

    auto tr_z_0_z_xxz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1238);

    auto tr_z_0_z_xxz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1239);

    auto tr_z_0_z_xxz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1240);

    auto tr_z_0_z_xxz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1241);

    auto tr_z_0_z_xxz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1242);

    auto tr_z_0_z_xxz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1243);

    auto tr_z_0_z_xxz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1244);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xx_zzzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxxxz, tr_xxzz_xxxyz, tr_xxzz_xxxzz, tr_xxzz_xxy, tr_xxzz_xxyyz, tr_xxzz_xxyzz, tr_xxzz_xxz, tr_xxzz_xxzzz, tr_xxzz_xyy, tr_xxzz_xyyyz, tr_xxzz_xyyzz, tr_xxzz_xyz, tr_xxzz_xyzzz, tr_xxzz_xzz, tr_xxzz_xzzzz, tr_xxzz_yyy, tr_xxzz_yyyyz, tr_xxzz_yyyzz, tr_xxzz_yyz, tr_xxzz_yyzzz, tr_xxzz_yzz, tr_xxzz_yzzzz, tr_xxzz_zzz, tr_xxzz_zzzzz, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xzzz, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yzzz, tr_xxzzz_zzzz, tr_z_0_z_xxz_xxxx, tr_z_0_z_xxz_xxxy, tr_z_0_z_xxz_xxxz, tr_z_0_z_xxz_xxyy, tr_z_0_z_xxz_xxyz, tr_z_0_z_xxz_xxzz, tr_z_0_z_xxz_xyyy, tr_z_0_z_xxz_xyyz, tr_z_0_z_xxz_xyzz, tr_z_0_z_xxz_xzzz, tr_z_0_z_xxz_yyyy, tr_z_0_z_xxz_yyyz, tr_z_0_z_xxz_yyzz, tr_z_0_z_xxz_yzzz, tr_z_0_z_xxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxz_xxxx[i] = -2.0 * tr_xx_xxxxz[i] * tke_0 - 6.0 * tr_xxz_xxxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xxxy[i] = -2.0 * tr_xx_xxxyz[i] * tke_0 - 6.0 * tr_xxz_xxxy[i] * tbe_0 + 4.0 * tr_xxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xxxz[i] = tr_xx_xxx[i] - 2.0 * tr_xx_xxxzz[i] * tke_0 - 6.0 * tr_xxz_xxxz[i] * tbe_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xxyy[i] = -2.0 * tr_xx_xxyyz[i] * tke_0 - 6.0 * tr_xxz_xxyy[i] * tbe_0 + 4.0 * tr_xxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xxyz[i] = tr_xx_xxy[i] - 2.0 * tr_xx_xxyzz[i] * tke_0 - 6.0 * tr_xxz_xxyz[i] * tbe_0 - 2.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xxzz[i] = 2.0 * tr_xx_xxz[i] - 2.0 * tr_xx_xxzzz[i] * tke_0 - 6.0 * tr_xxz_xxzz[i] * tbe_0 - 4.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xyyy[i] = -2.0 * tr_xx_xyyyz[i] * tke_0 - 6.0 * tr_xxz_xyyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xyyz[i] = tr_xx_xyy[i] - 2.0 * tr_xx_xyyzz[i] * tke_0 - 6.0 * tr_xxz_xyyz[i] * tbe_0 - 2.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xyzz[i] = 2.0 * tr_xx_xyz[i] - 2.0 * tr_xx_xyzzz[i] * tke_0 - 6.0 * tr_xxz_xyzz[i] * tbe_0 - 4.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xzzz[i] = 3.0 * tr_xx_xzz[i] - 2.0 * tr_xx_xzzzz[i] * tke_0 - 6.0 * tr_xxz_xzzz[i] * tbe_0 - 6.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_yyyy[i] = -2.0 * tr_xx_yyyyz[i] * tke_0 - 6.0 * tr_xxz_yyyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_yyyz[i] = tr_xx_yyy[i] - 2.0 * tr_xx_yyyzz[i] * tke_0 - 6.0 * tr_xxz_yyyz[i] * tbe_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_yyzz[i] = 2.0 * tr_xx_yyz[i] - 2.0 * tr_xx_yyzzz[i] * tke_0 - 6.0 * tr_xxz_yyzz[i] * tbe_0 - 4.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_yzzz[i] = 3.0 * tr_xx_yzz[i] - 2.0 * tr_xx_yzzzz[i] * tke_0 - 6.0 * tr_xxz_yzzz[i] * tbe_0 - 6.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_zzzz[i] = 4.0 * tr_xx_zzz[i] - 2.0 * tr_xx_zzzzz[i] * tke_0 - 6.0 * tr_xxz_zzzz[i] * tbe_0 - 8.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1245-1260 components of targeted buffer : FG

    auto tr_z_0_z_xyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 1245);

    auto tr_z_0_z_xyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 1246);

    auto tr_z_0_z_xyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 1247);

    auto tr_z_0_z_xyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 1248);

    auto tr_z_0_z_xyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 1249);

    auto tr_z_0_z_xyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 1250);

    auto tr_z_0_z_xyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 1251);

    auto tr_z_0_z_xyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 1252);

    auto tr_z_0_z_xyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 1253);

    auto tr_z_0_z_xyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 1254);

    auto tr_z_0_z_xyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1255);

    auto tr_z_0_z_xyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1256);

    auto tr_z_0_z_xyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1257);

    auto tr_z_0_z_xyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1258);

    auto tr_z_0_z_xyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1259);

    #pragma omp simd aligned(tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxxxz, tr_xyyz_xxxyz, tr_xyyz_xxxzz, tr_xyyz_xxy, tr_xyyz_xxyyz, tr_xyyz_xxyzz, tr_xyyz_xxz, tr_xyyz_xxzzz, tr_xyyz_xyy, tr_xyyz_xyyyz, tr_xyyz_xyyzz, tr_xyyz_xyz, tr_xyyz_xyzzz, tr_xyyz_xzz, tr_xyyz_xzzzz, tr_xyyz_yyy, tr_xyyz_yyyyz, tr_xyyz_yyyzz, tr_xyyz_yyz, tr_xyyz_yyzzz, tr_xyyz_yzz, tr_xyyz_yzzzz, tr_xyyz_zzz, tr_xyyz_zzzzz, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xzzz, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yzzz, tr_xyyzz_zzzz, tr_z_0_z_xyy_xxxx, tr_z_0_z_xyy_xxxy, tr_z_0_z_xyy_xxxz, tr_z_0_z_xyy_xxyy, tr_z_0_z_xyy_xxyz, tr_z_0_z_xyy_xxzz, tr_z_0_z_xyy_xyyy, tr_z_0_z_xyy_xyyz, tr_z_0_z_xyy_xyzz, tr_z_0_z_xyy_xzzz, tr_z_0_z_xyy_yyyy, tr_z_0_z_xyy_yyyz, tr_z_0_z_xyy_yyzz, tr_z_0_z_xyy_yzzz, tr_z_0_z_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyy_xxxx[i] = -2.0 * tr_xyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xxxy[i] = -2.0 * tr_xyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xxxz[i] = -2.0 * tr_xyy_xxxz[i] * tbe_0 - 2.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xxyy[i] = -2.0 * tr_xyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xxyz[i] = -2.0 * tr_xyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xxzz[i] = -2.0 * tr_xyy_xxzz[i] * tbe_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xyyy[i] = -2.0 * tr_xyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xyyz[i] = -2.0 * tr_xyy_xyyz[i] * tbe_0 - 2.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xyzz[i] = -2.0 * tr_xyy_xyzz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xzzz[i] = -2.0 * tr_xyy_xzzz[i] * tbe_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_yyyy[i] = -2.0 * tr_xyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_yyyz[i] = -2.0 * tr_xyy_yyyz[i] * tbe_0 - 2.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_yyzz[i] = -2.0 * tr_xyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_yzzz[i] = -2.0 * tr_xyy_yzzz[i] * tbe_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_zzzz[i] = -2.0 * tr_xyy_zzzz[i] * tbe_0 - 8.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1260-1275 components of targeted buffer : FG

    auto tr_z_0_z_xyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1260);

    auto tr_z_0_z_xyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1261);

    auto tr_z_0_z_xyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1262);

    auto tr_z_0_z_xyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1263);

    auto tr_z_0_z_xyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1264);

    auto tr_z_0_z_xyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1265);

    auto tr_z_0_z_xyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1266);

    auto tr_z_0_z_xyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1267);

    auto tr_z_0_z_xyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1268);

    auto tr_z_0_z_xyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1269);

    auto tr_z_0_z_xyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1270);

    auto tr_z_0_z_xyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1271);

    auto tr_z_0_z_xyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1272);

    auto tr_z_0_z_xyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1273);

    auto tr_z_0_z_xyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1274);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxxxz, tr_xyzz_xxxyz, tr_xyzz_xxxzz, tr_xyzz_xxy, tr_xyzz_xxyyz, tr_xyzz_xxyzz, tr_xyzz_xxz, tr_xyzz_xxzzz, tr_xyzz_xyy, tr_xyzz_xyyyz, tr_xyzz_xyyzz, tr_xyzz_xyz, tr_xyzz_xyzzz, tr_xyzz_xzz, tr_xyzz_xzzzz, tr_xyzz_yyy, tr_xyzz_yyyyz, tr_xyzz_yyyzz, tr_xyzz_yyz, tr_xyzz_yyzzz, tr_xyzz_yzz, tr_xyzz_yzzzz, tr_xyzz_zzz, tr_xyzz_zzzzz, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xzzz, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yzzz, tr_xyzzz_zzzz, tr_z_0_z_xyz_xxxx, tr_z_0_z_xyz_xxxy, tr_z_0_z_xyz_xxxz, tr_z_0_z_xyz_xxyy, tr_z_0_z_xyz_xxyz, tr_z_0_z_xyz_xxzz, tr_z_0_z_xyz_xyyy, tr_z_0_z_xyz_xyyz, tr_z_0_z_xyz_xyzz, tr_z_0_z_xyz_xzzz, tr_z_0_z_xyz_yyyy, tr_z_0_z_xyz_yyyz, tr_z_0_z_xyz_yyzz, tr_z_0_z_xyz_yzzz, tr_z_0_z_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyz_xxxx[i] = -2.0 * tr_xy_xxxxz[i] * tke_0 - 6.0 * tr_xyz_xxxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xxxy[i] = -2.0 * tr_xy_xxxyz[i] * tke_0 - 6.0 * tr_xyz_xxxy[i] * tbe_0 + 4.0 * tr_xyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xxxz[i] = tr_xy_xxx[i] - 2.0 * tr_xy_xxxzz[i] * tke_0 - 6.0 * tr_xyz_xxxz[i] * tbe_0 - 2.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xxyy[i] = -2.0 * tr_xy_xxyyz[i] * tke_0 - 6.0 * tr_xyz_xxyy[i] * tbe_0 + 4.0 * tr_xyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xxyz[i] = tr_xy_xxy[i] - 2.0 * tr_xy_xxyzz[i] * tke_0 - 6.0 * tr_xyz_xxyz[i] * tbe_0 - 2.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xxzz[i] = 2.0 * tr_xy_xxz[i] - 2.0 * tr_xy_xxzzz[i] * tke_0 - 6.0 * tr_xyz_xxzz[i] * tbe_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xyyy[i] = -2.0 * tr_xy_xyyyz[i] * tke_0 - 6.0 * tr_xyz_xyyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xyyz[i] = tr_xy_xyy[i] - 2.0 * tr_xy_xyyzz[i] * tke_0 - 6.0 * tr_xyz_xyyz[i] * tbe_0 - 2.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xyzz[i] = 2.0 * tr_xy_xyz[i] - 2.0 * tr_xy_xyzzz[i] * tke_0 - 6.0 * tr_xyz_xyzz[i] * tbe_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xzzz[i] = 3.0 * tr_xy_xzz[i] - 2.0 * tr_xy_xzzzz[i] * tke_0 - 6.0 * tr_xyz_xzzz[i] * tbe_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_yyyy[i] = -2.0 * tr_xy_yyyyz[i] * tke_0 - 6.0 * tr_xyz_yyyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_yyyz[i] = tr_xy_yyy[i] - 2.0 * tr_xy_yyyzz[i] * tke_0 - 6.0 * tr_xyz_yyyz[i] * tbe_0 - 2.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_yyzz[i] = 2.0 * tr_xy_yyz[i] - 2.0 * tr_xy_yyzzz[i] * tke_0 - 6.0 * tr_xyz_yyzz[i] * tbe_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_yzzz[i] = 3.0 * tr_xy_yzz[i] - 2.0 * tr_xy_yzzzz[i] * tke_0 - 6.0 * tr_xyz_yzzz[i] * tbe_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_zzzz[i] = 4.0 * tr_xy_zzz[i] - 2.0 * tr_xy_zzzzz[i] * tke_0 - 6.0 * tr_xyz_zzzz[i] * tbe_0 - 8.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1275-1290 components of targeted buffer : FG

    auto tr_z_0_z_xzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1275);

    auto tr_z_0_z_xzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1276);

    auto tr_z_0_z_xzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1277);

    auto tr_z_0_z_xzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1278);

    auto tr_z_0_z_xzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1279);

    auto tr_z_0_z_xzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1280);

    auto tr_z_0_z_xzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1281);

    auto tr_z_0_z_xzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1282);

    auto tr_z_0_z_xzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1283);

    auto tr_z_0_z_xzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1284);

    auto tr_z_0_z_xzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1285);

    auto tr_z_0_z_xzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1286);

    auto tr_z_0_z_xzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1287);

    auto tr_z_0_z_xzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1288);

    auto tr_z_0_z_xzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1289);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxxxz, tr_xzzz_xxxyz, tr_xzzz_xxxzz, tr_xzzz_xxy, tr_xzzz_xxyyz, tr_xzzz_xxyzz, tr_xzzz_xxz, tr_xzzz_xxzzz, tr_xzzz_xyy, tr_xzzz_xyyyz, tr_xzzz_xyyzz, tr_xzzz_xyz, tr_xzzz_xyzzz, tr_xzzz_xzz, tr_xzzz_xzzzz, tr_xzzz_yyy, tr_xzzz_yyyyz, tr_xzzz_yyyzz, tr_xzzz_yyz, tr_xzzz_yyzzz, tr_xzzz_yzz, tr_xzzz_yzzzz, tr_xzzz_zzz, tr_xzzz_zzzzz, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xzzz, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yzzz, tr_xzzzz_zzzz, tr_z_0_z_xzz_xxxx, tr_z_0_z_xzz_xxxy, tr_z_0_z_xzz_xxxz, tr_z_0_z_xzz_xxyy, tr_z_0_z_xzz_xxyz, tr_z_0_z_xzz_xxzz, tr_z_0_z_xzz_xyyy, tr_z_0_z_xzz_xyyz, tr_z_0_z_xzz_xyzz, tr_z_0_z_xzz_xzzz, tr_z_0_z_xzz_yyyy, tr_z_0_z_xzz_yyyz, tr_z_0_z_xzz_yyzz, tr_z_0_z_xzz_yzzz, tr_z_0_z_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xzz_xxxx[i] = 2.0 * tr_x_xxxx[i] - 4.0 * tr_xz_xxxxz[i] * tke_0 - 10.0 * tr_xzz_xxxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xxxy[i] = 2.0 * tr_x_xxxy[i] - 4.0 * tr_xz_xxxyz[i] * tke_0 - 10.0 * tr_xzz_xxxy[i] * tbe_0 + 4.0 * tr_xzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xxxz[i] = 2.0 * tr_x_xxxz[i] + 2.0 * tr_xz_xxx[i] - 4.0 * tr_xz_xxxzz[i] * tke_0 - 10.0 * tr_xzz_xxxz[i] * tbe_0 - 2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xxyy[i] = 2.0 * tr_x_xxyy[i] - 4.0 * tr_xz_xxyyz[i] * tke_0 - 10.0 * tr_xzz_xxyy[i] * tbe_0 + 4.0 * tr_xzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xxyz[i] = 2.0 * tr_x_xxyz[i] + 2.0 * tr_xz_xxy[i] - 4.0 * tr_xz_xxyzz[i] * tke_0 - 10.0 * tr_xzz_xxyz[i] * tbe_0 - 2.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xxzz[i] = 2.0 * tr_x_xxzz[i] + 4.0 * tr_xz_xxz[i] - 4.0 * tr_xz_xxzzz[i] * tke_0 - 10.0 * tr_xzz_xxzz[i] * tbe_0 - 4.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xyyy[i] = 2.0 * tr_x_xyyy[i] - 4.0 * tr_xz_xyyyz[i] * tke_0 - 10.0 * tr_xzz_xyyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xyyz[i] = 2.0 * tr_x_xyyz[i] + 2.0 * tr_xz_xyy[i] - 4.0 * tr_xz_xyyzz[i] * tke_0 - 10.0 * tr_xzz_xyyz[i] * tbe_0 - 2.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xyzz[i] = 2.0 * tr_x_xyzz[i] + 4.0 * tr_xz_xyz[i] - 4.0 * tr_xz_xyzzz[i] * tke_0 - 10.0 * tr_xzz_xyzz[i] * tbe_0 - 4.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xzzz[i] = 2.0 * tr_x_xzzz[i] + 6.0 * tr_xz_xzz[i] - 4.0 * tr_xz_xzzzz[i] * tke_0 - 10.0 * tr_xzz_xzzz[i] * tbe_0 - 6.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_yyyy[i] = 2.0 * tr_x_yyyy[i] - 4.0 * tr_xz_yyyyz[i] * tke_0 - 10.0 * tr_xzz_yyyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_yyyz[i] = 2.0 * tr_x_yyyz[i] + 2.0 * tr_xz_yyy[i] - 4.0 * tr_xz_yyyzz[i] * tke_0 - 10.0 * tr_xzz_yyyz[i] * tbe_0 - 2.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_yyzz[i] = 2.0 * tr_x_yyzz[i] + 4.0 * tr_xz_yyz[i] - 4.0 * tr_xz_yyzzz[i] * tke_0 - 10.0 * tr_xzz_yyzz[i] * tbe_0 - 4.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_yzzz[i] = 2.0 * tr_x_yzzz[i] + 6.0 * tr_xz_yzz[i] - 4.0 * tr_xz_yzzzz[i] * tke_0 - 10.0 * tr_xzz_yzzz[i] * tbe_0 - 6.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_zzzz[i] = 2.0 * tr_x_zzzz[i] + 8.0 * tr_xz_zzz[i] - 4.0 * tr_xz_zzzzz[i] * tke_0 - 10.0 * tr_xzz_zzzz[i] * tbe_0 - 8.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1290-1305 components of targeted buffer : FG

    auto tr_z_0_z_yyy_xxxx = pbuffer.data(idx_op_geom_110_fg + 1290);

    auto tr_z_0_z_yyy_xxxy = pbuffer.data(idx_op_geom_110_fg + 1291);

    auto tr_z_0_z_yyy_xxxz = pbuffer.data(idx_op_geom_110_fg + 1292);

    auto tr_z_0_z_yyy_xxyy = pbuffer.data(idx_op_geom_110_fg + 1293);

    auto tr_z_0_z_yyy_xxyz = pbuffer.data(idx_op_geom_110_fg + 1294);

    auto tr_z_0_z_yyy_xxzz = pbuffer.data(idx_op_geom_110_fg + 1295);

    auto tr_z_0_z_yyy_xyyy = pbuffer.data(idx_op_geom_110_fg + 1296);

    auto tr_z_0_z_yyy_xyyz = pbuffer.data(idx_op_geom_110_fg + 1297);

    auto tr_z_0_z_yyy_xyzz = pbuffer.data(idx_op_geom_110_fg + 1298);

    auto tr_z_0_z_yyy_xzzz = pbuffer.data(idx_op_geom_110_fg + 1299);

    auto tr_z_0_z_yyy_yyyy = pbuffer.data(idx_op_geom_110_fg + 1300);

    auto tr_z_0_z_yyy_yyyz = pbuffer.data(idx_op_geom_110_fg + 1301);

    auto tr_z_0_z_yyy_yyzz = pbuffer.data(idx_op_geom_110_fg + 1302);

    auto tr_z_0_z_yyy_yzzz = pbuffer.data(idx_op_geom_110_fg + 1303);

    auto tr_z_0_z_yyy_zzzz = pbuffer.data(idx_op_geom_110_fg + 1304);

    #pragma omp simd aligned(tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxxxz, tr_yyyz_xxxyz, tr_yyyz_xxxzz, tr_yyyz_xxy, tr_yyyz_xxyyz, tr_yyyz_xxyzz, tr_yyyz_xxz, tr_yyyz_xxzzz, tr_yyyz_xyy, tr_yyyz_xyyyz, tr_yyyz_xyyzz, tr_yyyz_xyz, tr_yyyz_xyzzz, tr_yyyz_xzz, tr_yyyz_xzzzz, tr_yyyz_yyy, tr_yyyz_yyyyz, tr_yyyz_yyyzz, tr_yyyz_yyz, tr_yyyz_yyzzz, tr_yyyz_yzz, tr_yyyz_yzzzz, tr_yyyz_zzz, tr_yyyz_zzzzz, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xzzz, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yzzz, tr_yyyzz_zzzz, tr_z_0_z_yyy_xxxx, tr_z_0_z_yyy_xxxy, tr_z_0_z_yyy_xxxz, tr_z_0_z_yyy_xxyy, tr_z_0_z_yyy_xxyz, tr_z_0_z_yyy_xxzz, tr_z_0_z_yyy_xyyy, tr_z_0_z_yyy_xyyz, tr_z_0_z_yyy_xyzz, tr_z_0_z_yyy_xzzz, tr_z_0_z_yyy_yyyy, tr_z_0_z_yyy_yyyz, tr_z_0_z_yyy_yyzz, tr_z_0_z_yyy_yzzz, tr_z_0_z_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyy_xxxx[i] = -2.0 * tr_yyy_xxxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xxxy[i] = -2.0 * tr_yyy_xxxy[i] * tbe_0 + 4.0 * tr_yyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xxxz[i] = -2.0 * tr_yyy_xxxz[i] * tbe_0 - 2.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xxyy[i] = -2.0 * tr_yyy_xxyy[i] * tbe_0 + 4.0 * tr_yyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xxyz[i] = -2.0 * tr_yyy_xxyz[i] * tbe_0 - 2.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xxzz[i] = -2.0 * tr_yyy_xxzz[i] * tbe_0 - 4.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xyyy[i] = -2.0 * tr_yyy_xyyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xyyz[i] = -2.0 * tr_yyy_xyyz[i] * tbe_0 - 2.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xyzz[i] = -2.0 * tr_yyy_xyzz[i] * tbe_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xzzz[i] = -2.0 * tr_yyy_xzzz[i] * tbe_0 - 6.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_yyyy[i] = -2.0 * tr_yyy_yyyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_yyyz[i] = -2.0 * tr_yyy_yyyz[i] * tbe_0 - 2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_yyzz[i] = -2.0 * tr_yyy_yyzz[i] * tbe_0 - 4.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_yzzz[i] = -2.0 * tr_yyy_yzzz[i] * tbe_0 - 6.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_zzzz[i] = -2.0 * tr_yyy_zzzz[i] * tbe_0 - 8.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1305-1320 components of targeted buffer : FG

    auto tr_z_0_z_yyz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1305);

    auto tr_z_0_z_yyz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1306);

    auto tr_z_0_z_yyz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1307);

    auto tr_z_0_z_yyz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1308);

    auto tr_z_0_z_yyz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1309);

    auto tr_z_0_z_yyz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1310);

    auto tr_z_0_z_yyz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1311);

    auto tr_z_0_z_yyz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1312);

    auto tr_z_0_z_yyz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1313);

    auto tr_z_0_z_yyz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1314);

    auto tr_z_0_z_yyz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1315);

    auto tr_z_0_z_yyz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1316);

    auto tr_z_0_z_yyz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1317);

    auto tr_z_0_z_yyz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1318);

    auto tr_z_0_z_yyz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1319);

    #pragma omp simd aligned(tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yy_zzzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxxxz, tr_yyzz_xxxyz, tr_yyzz_xxxzz, tr_yyzz_xxy, tr_yyzz_xxyyz, tr_yyzz_xxyzz, tr_yyzz_xxz, tr_yyzz_xxzzz, tr_yyzz_xyy, tr_yyzz_xyyyz, tr_yyzz_xyyzz, tr_yyzz_xyz, tr_yyzz_xyzzz, tr_yyzz_xzz, tr_yyzz_xzzzz, tr_yyzz_yyy, tr_yyzz_yyyyz, tr_yyzz_yyyzz, tr_yyzz_yyz, tr_yyzz_yyzzz, tr_yyzz_yzz, tr_yyzz_yzzzz, tr_yyzz_zzz, tr_yyzz_zzzzz, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xzzz, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yzzz, tr_yyzzz_zzzz, tr_z_0_z_yyz_xxxx, tr_z_0_z_yyz_xxxy, tr_z_0_z_yyz_xxxz, tr_z_0_z_yyz_xxyy, tr_z_0_z_yyz_xxyz, tr_z_0_z_yyz_xxzz, tr_z_0_z_yyz_xyyy, tr_z_0_z_yyz_xyyz, tr_z_0_z_yyz_xyzz, tr_z_0_z_yyz_xzzz, tr_z_0_z_yyz_yyyy, tr_z_0_z_yyz_yyyz, tr_z_0_z_yyz_yyzz, tr_z_0_z_yyz_yzzz, tr_z_0_z_yyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyz_xxxx[i] = -2.0 * tr_yy_xxxxz[i] * tke_0 - 6.0 * tr_yyz_xxxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xxxy[i] = -2.0 * tr_yy_xxxyz[i] * tke_0 - 6.0 * tr_yyz_xxxy[i] * tbe_0 + 4.0 * tr_yyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xxxz[i] = tr_yy_xxx[i] - 2.0 * tr_yy_xxxzz[i] * tke_0 - 6.0 * tr_yyz_xxxz[i] * tbe_0 - 2.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xxyy[i] = -2.0 * tr_yy_xxyyz[i] * tke_0 - 6.0 * tr_yyz_xxyy[i] * tbe_0 + 4.0 * tr_yyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xxyz[i] = tr_yy_xxy[i] - 2.0 * tr_yy_xxyzz[i] * tke_0 - 6.0 * tr_yyz_xxyz[i] * tbe_0 - 2.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xxzz[i] = 2.0 * tr_yy_xxz[i] - 2.0 * tr_yy_xxzzz[i] * tke_0 - 6.0 * tr_yyz_xxzz[i] * tbe_0 - 4.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xyyy[i] = -2.0 * tr_yy_xyyyz[i] * tke_0 - 6.0 * tr_yyz_xyyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xyyz[i] = tr_yy_xyy[i] - 2.0 * tr_yy_xyyzz[i] * tke_0 - 6.0 * tr_yyz_xyyz[i] * tbe_0 - 2.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xyzz[i] = 2.0 * tr_yy_xyz[i] - 2.0 * tr_yy_xyzzz[i] * tke_0 - 6.0 * tr_yyz_xyzz[i] * tbe_0 - 4.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xzzz[i] = 3.0 * tr_yy_xzz[i] - 2.0 * tr_yy_xzzzz[i] * tke_0 - 6.0 * tr_yyz_xzzz[i] * tbe_0 - 6.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_yyyy[i] = -2.0 * tr_yy_yyyyz[i] * tke_0 - 6.0 * tr_yyz_yyyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_yyyz[i] = tr_yy_yyy[i] - 2.0 * tr_yy_yyyzz[i] * tke_0 - 6.0 * tr_yyz_yyyz[i] * tbe_0 - 2.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_yyzz[i] = 2.0 * tr_yy_yyz[i] - 2.0 * tr_yy_yyzzz[i] * tke_0 - 6.0 * tr_yyz_yyzz[i] * tbe_0 - 4.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_yzzz[i] = 3.0 * tr_yy_yzz[i] - 2.0 * tr_yy_yzzzz[i] * tke_0 - 6.0 * tr_yyz_yzzz[i] * tbe_0 - 6.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_zzzz[i] = 4.0 * tr_yy_zzz[i] - 2.0 * tr_yy_zzzzz[i] * tke_0 - 6.0 * tr_yyz_zzzz[i] * tbe_0 - 8.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1320-1335 components of targeted buffer : FG

    auto tr_z_0_z_yzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1320);

    auto tr_z_0_z_yzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1321);

    auto tr_z_0_z_yzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1322);

    auto tr_z_0_z_yzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1323);

    auto tr_z_0_z_yzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1324);

    auto tr_z_0_z_yzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1325);

    auto tr_z_0_z_yzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1326);

    auto tr_z_0_z_yzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1327);

    auto tr_z_0_z_yzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1328);

    auto tr_z_0_z_yzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1329);

    auto tr_z_0_z_yzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1330);

    auto tr_z_0_z_yzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1331);

    auto tr_z_0_z_yzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1332);

    auto tr_z_0_z_yzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1333);

    auto tr_z_0_z_yzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1334);

    #pragma omp simd aligned(tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxxxz, tr_yzzz_xxxyz, tr_yzzz_xxxzz, tr_yzzz_xxy, tr_yzzz_xxyyz, tr_yzzz_xxyzz, tr_yzzz_xxz, tr_yzzz_xxzzz, tr_yzzz_xyy, tr_yzzz_xyyyz, tr_yzzz_xyyzz, tr_yzzz_xyz, tr_yzzz_xyzzz, tr_yzzz_xzz, tr_yzzz_xzzzz, tr_yzzz_yyy, tr_yzzz_yyyyz, tr_yzzz_yyyzz, tr_yzzz_yyz, tr_yzzz_yyzzz, tr_yzzz_yzz, tr_yzzz_yzzzz, tr_yzzz_zzz, tr_yzzz_zzzzz, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xzzz, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yzzz, tr_yzzzz_zzzz, tr_z_0_z_yzz_xxxx, tr_z_0_z_yzz_xxxy, tr_z_0_z_yzz_xxxz, tr_z_0_z_yzz_xxyy, tr_z_0_z_yzz_xxyz, tr_z_0_z_yzz_xxzz, tr_z_0_z_yzz_xyyy, tr_z_0_z_yzz_xyyz, tr_z_0_z_yzz_xyzz, tr_z_0_z_yzz_xzzz, tr_z_0_z_yzz_yyyy, tr_z_0_z_yzz_yyyz, tr_z_0_z_yzz_yyzz, tr_z_0_z_yzz_yzzz, tr_z_0_z_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yzz_xxxx[i] = 2.0 * tr_y_xxxx[i] - 4.0 * tr_yz_xxxxz[i] * tke_0 - 10.0 * tr_yzz_xxxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xxxy[i] = 2.0 * tr_y_xxxy[i] - 4.0 * tr_yz_xxxyz[i] * tke_0 - 10.0 * tr_yzz_xxxy[i] * tbe_0 + 4.0 * tr_yzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xxxz[i] = 2.0 * tr_y_xxxz[i] + 2.0 * tr_yz_xxx[i] - 4.0 * tr_yz_xxxzz[i] * tke_0 - 10.0 * tr_yzz_xxxz[i] * tbe_0 - 2.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xxyy[i] = 2.0 * tr_y_xxyy[i] - 4.0 * tr_yz_xxyyz[i] * tke_0 - 10.0 * tr_yzz_xxyy[i] * tbe_0 + 4.0 * tr_yzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xxyz[i] = 2.0 * tr_y_xxyz[i] + 2.0 * tr_yz_xxy[i] - 4.0 * tr_yz_xxyzz[i] * tke_0 - 10.0 * tr_yzz_xxyz[i] * tbe_0 - 2.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xxzz[i] = 2.0 * tr_y_xxzz[i] + 4.0 * tr_yz_xxz[i] - 4.0 * tr_yz_xxzzz[i] * tke_0 - 10.0 * tr_yzz_xxzz[i] * tbe_0 - 4.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xyyy[i] = 2.0 * tr_y_xyyy[i] - 4.0 * tr_yz_xyyyz[i] * tke_0 - 10.0 * tr_yzz_xyyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xyyz[i] = 2.0 * tr_y_xyyz[i] + 2.0 * tr_yz_xyy[i] - 4.0 * tr_yz_xyyzz[i] * tke_0 - 10.0 * tr_yzz_xyyz[i] * tbe_0 - 2.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xyzz[i] = 2.0 * tr_y_xyzz[i] + 4.0 * tr_yz_xyz[i] - 4.0 * tr_yz_xyzzz[i] * tke_0 - 10.0 * tr_yzz_xyzz[i] * tbe_0 - 4.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xzzz[i] = 2.0 * tr_y_xzzz[i] + 6.0 * tr_yz_xzz[i] - 4.0 * tr_yz_xzzzz[i] * tke_0 - 10.0 * tr_yzz_xzzz[i] * tbe_0 - 6.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_yyyy[i] = 2.0 * tr_y_yyyy[i] - 4.0 * tr_yz_yyyyz[i] * tke_0 - 10.0 * tr_yzz_yyyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_yyyz[i] = 2.0 * tr_y_yyyz[i] + 2.0 * tr_yz_yyy[i] - 4.0 * tr_yz_yyyzz[i] * tke_0 - 10.0 * tr_yzz_yyyz[i] * tbe_0 - 2.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_yyzz[i] = 2.0 * tr_y_yyzz[i] + 4.0 * tr_yz_yyz[i] - 4.0 * tr_yz_yyzzz[i] * tke_0 - 10.0 * tr_yzz_yyzz[i] * tbe_0 - 4.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_yzzz[i] = 2.0 * tr_y_yzzz[i] + 6.0 * tr_yz_yzz[i] - 4.0 * tr_yz_yzzzz[i] * tke_0 - 10.0 * tr_yzz_yzzz[i] * tbe_0 - 6.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_zzzz[i] = 2.0 * tr_y_zzzz[i] + 8.0 * tr_yz_zzz[i] - 4.0 * tr_yz_zzzzz[i] * tke_0 - 10.0 * tr_yzz_zzzz[i] * tbe_0 - 8.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1335-1350 components of targeted buffer : FG

    auto tr_z_0_z_zzz_xxxx = pbuffer.data(idx_op_geom_110_fg + 1335);

    auto tr_z_0_z_zzz_xxxy = pbuffer.data(idx_op_geom_110_fg + 1336);

    auto tr_z_0_z_zzz_xxxz = pbuffer.data(idx_op_geom_110_fg + 1337);

    auto tr_z_0_z_zzz_xxyy = pbuffer.data(idx_op_geom_110_fg + 1338);

    auto tr_z_0_z_zzz_xxyz = pbuffer.data(idx_op_geom_110_fg + 1339);

    auto tr_z_0_z_zzz_xxzz = pbuffer.data(idx_op_geom_110_fg + 1340);

    auto tr_z_0_z_zzz_xyyy = pbuffer.data(idx_op_geom_110_fg + 1341);

    auto tr_z_0_z_zzz_xyyz = pbuffer.data(idx_op_geom_110_fg + 1342);

    auto tr_z_0_z_zzz_xyzz = pbuffer.data(idx_op_geom_110_fg + 1343);

    auto tr_z_0_z_zzz_xzzz = pbuffer.data(idx_op_geom_110_fg + 1344);

    auto tr_z_0_z_zzz_yyyy = pbuffer.data(idx_op_geom_110_fg + 1345);

    auto tr_z_0_z_zzz_yyyz = pbuffer.data(idx_op_geom_110_fg + 1346);

    auto tr_z_0_z_zzz_yyzz = pbuffer.data(idx_op_geom_110_fg + 1347);

    auto tr_z_0_z_zzz_yzzz = pbuffer.data(idx_op_geom_110_fg + 1348);

    auto tr_z_0_z_zzz_zzzz = pbuffer.data(idx_op_geom_110_fg + 1349);

    #pragma omp simd aligned(tr_z_0_z_zzz_xxxx, tr_z_0_z_zzz_xxxy, tr_z_0_z_zzz_xxxz, tr_z_0_z_zzz_xxyy, tr_z_0_z_zzz_xxyz, tr_z_0_z_zzz_xxzz, tr_z_0_z_zzz_xyyy, tr_z_0_z_zzz_xyyz, tr_z_0_z_zzz_xyzz, tr_z_0_z_zzz_xzzz, tr_z_0_z_zzz_yyyy, tr_z_0_z_zzz_yyyz, tr_z_0_z_zzz_yyzz, tr_z_0_z_zzz_yzzz, tr_z_0_z_zzz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxxxz, tr_zzzz_xxxyz, tr_zzzz_xxxzz, tr_zzzz_xxy, tr_zzzz_xxyyz, tr_zzzz_xxyzz, tr_zzzz_xxz, tr_zzzz_xxzzz, tr_zzzz_xyy, tr_zzzz_xyyyz, tr_zzzz_xyyzz, tr_zzzz_xyz, tr_zzzz_xyzzz, tr_zzzz_xzz, tr_zzzz_xzzzz, tr_zzzz_yyy, tr_zzzz_yyyyz, tr_zzzz_yyyzz, tr_zzzz_yyz, tr_zzzz_yyzzz, tr_zzzz_yzz, tr_zzzz_yzzzz, tr_zzzz_zzz, tr_zzzz_zzzzz, tr_zzzzz_xxxx, tr_zzzzz_xxxy, tr_zzzzz_xxxz, tr_zzzzz_xxyy, tr_zzzzz_xxyz, tr_zzzzz_xxzz, tr_zzzzz_xyyy, tr_zzzzz_xyyz, tr_zzzzz_xyzz, tr_zzzzz_xzzz, tr_zzzzz_yyyy, tr_zzzzz_yyyz, tr_zzzzz_yyzz, tr_zzzzz_yzzz, tr_zzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zzz_xxxx[i] = 6.0 * tr_z_xxxx[i] - 6.0 * tr_zz_xxxxz[i] * tke_0 - 14.0 * tr_zzz_xxxx[i] * tbe_0 + 4.0 * tr_zzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xxxy[i] = 6.0 * tr_z_xxxy[i] - 6.0 * tr_zz_xxxyz[i] * tke_0 - 14.0 * tr_zzz_xxxy[i] * tbe_0 + 4.0 * tr_zzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xxxz[i] = 6.0 * tr_z_xxxz[i] + 3.0 * tr_zz_xxx[i] - 6.0 * tr_zz_xxxzz[i] * tke_0 - 14.0 * tr_zzz_xxxz[i] * tbe_0 - 2.0 * tr_zzzz_xxx[i] * tbe_0 + 4.0 * tr_zzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xxyy[i] = 6.0 * tr_z_xxyy[i] - 6.0 * tr_zz_xxyyz[i] * tke_0 - 14.0 * tr_zzz_xxyy[i] * tbe_0 + 4.0 * tr_zzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xxyz[i] = 6.0 * tr_z_xxyz[i] + 3.0 * tr_zz_xxy[i] - 6.0 * tr_zz_xxyzz[i] * tke_0 - 14.0 * tr_zzz_xxyz[i] * tbe_0 - 2.0 * tr_zzzz_xxy[i] * tbe_0 + 4.0 * tr_zzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xxzz[i] = 6.0 * tr_z_xxzz[i] + 6.0 * tr_zz_xxz[i] - 6.0 * tr_zz_xxzzz[i] * tke_0 - 14.0 * tr_zzz_xxzz[i] * tbe_0 - 4.0 * tr_zzzz_xxz[i] * tbe_0 + 4.0 * tr_zzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xyyy[i] = 6.0 * tr_z_xyyy[i] - 6.0 * tr_zz_xyyyz[i] * tke_0 - 14.0 * tr_zzz_xyyy[i] * tbe_0 + 4.0 * tr_zzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xyyz[i] = 6.0 * tr_z_xyyz[i] + 3.0 * tr_zz_xyy[i] - 6.0 * tr_zz_xyyzz[i] * tke_0 - 14.0 * tr_zzz_xyyz[i] * tbe_0 - 2.0 * tr_zzzz_xyy[i] * tbe_0 + 4.0 * tr_zzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xyzz[i] = 6.0 * tr_z_xyzz[i] + 6.0 * tr_zz_xyz[i] - 6.0 * tr_zz_xyzzz[i] * tke_0 - 14.0 * tr_zzz_xyzz[i] * tbe_0 - 4.0 * tr_zzzz_xyz[i] * tbe_0 + 4.0 * tr_zzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xzzz[i] = 6.0 * tr_z_xzzz[i] + 9.0 * tr_zz_xzz[i] - 6.0 * tr_zz_xzzzz[i] * tke_0 - 14.0 * tr_zzz_xzzz[i] * tbe_0 - 6.0 * tr_zzzz_xzz[i] * tbe_0 + 4.0 * tr_zzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_yyyy[i] = 6.0 * tr_z_yyyy[i] - 6.0 * tr_zz_yyyyz[i] * tke_0 - 14.0 * tr_zzz_yyyy[i] * tbe_0 + 4.0 * tr_zzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_yyyz[i] = 6.0 * tr_z_yyyz[i] + 3.0 * tr_zz_yyy[i] - 6.0 * tr_zz_yyyzz[i] * tke_0 - 14.0 * tr_zzz_yyyz[i] * tbe_0 - 2.0 * tr_zzzz_yyy[i] * tbe_0 + 4.0 * tr_zzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_yyzz[i] = 6.0 * tr_z_yyzz[i] + 6.0 * tr_zz_yyz[i] - 6.0 * tr_zz_yyzzz[i] * tke_0 - 14.0 * tr_zzz_yyzz[i] * tbe_0 - 4.0 * tr_zzzz_yyz[i] * tbe_0 + 4.0 * tr_zzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_yzzz[i] = 6.0 * tr_z_yzzz[i] + 9.0 * tr_zz_yzz[i] - 6.0 * tr_zz_yzzzz[i] * tke_0 - 14.0 * tr_zzz_yzzz[i] * tbe_0 - 6.0 * tr_zzzz_yzz[i] * tbe_0 + 4.0 * tr_zzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_zzzz[i] = 6.0 * tr_z_zzzz[i] + 12.0 * tr_zz_zzz[i] - 6.0 * tr_zz_zzzzz[i] * tke_0 - 14.0 * tr_zzz_zzzz[i] * tbe_0 - 8.0 * tr_zzzz_zzz[i] * tbe_0 + 4.0 * tr_zzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

