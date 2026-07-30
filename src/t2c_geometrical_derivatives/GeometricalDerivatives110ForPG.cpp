#include "GeometricalDerivatives110ForPG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_pg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_pg,
                         const int idx_op_sf,
                         const int idx_op_sh,
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

    // Set up components of auxiliary buffer : SF

    auto tr_0_xxx = pbuffer.data(idx_op_sf);

    auto tr_0_xxy = pbuffer.data(idx_op_sf + 1);

    auto tr_0_xxz = pbuffer.data(idx_op_sf + 2);

    auto tr_0_xyy = pbuffer.data(idx_op_sf + 3);

    auto tr_0_xyz = pbuffer.data(idx_op_sf + 4);

    auto tr_0_xzz = pbuffer.data(idx_op_sf + 5);

    auto tr_0_yyy = pbuffer.data(idx_op_sf + 6);

    auto tr_0_yyz = pbuffer.data(idx_op_sf + 7);

    auto tr_0_yzz = pbuffer.data(idx_op_sf + 8);

    auto tr_0_zzz = pbuffer.data(idx_op_sf + 9);

    // Set up components of auxiliary buffer : SH

    auto tr_0_xxxxx = pbuffer.data(idx_op_sh);

    auto tr_0_xxxxy = pbuffer.data(idx_op_sh + 1);

    auto tr_0_xxxxz = pbuffer.data(idx_op_sh + 2);

    auto tr_0_xxxyy = pbuffer.data(idx_op_sh + 3);

    auto tr_0_xxxyz = pbuffer.data(idx_op_sh + 4);

    auto tr_0_xxxzz = pbuffer.data(idx_op_sh + 5);

    auto tr_0_xxyyy = pbuffer.data(idx_op_sh + 6);

    auto tr_0_xxyyz = pbuffer.data(idx_op_sh + 7);

    auto tr_0_xxyzz = pbuffer.data(idx_op_sh + 8);

    auto tr_0_xxzzz = pbuffer.data(idx_op_sh + 9);

    auto tr_0_xyyyy = pbuffer.data(idx_op_sh + 10);

    auto tr_0_xyyyz = pbuffer.data(idx_op_sh + 11);

    auto tr_0_xyyzz = pbuffer.data(idx_op_sh + 12);

    auto tr_0_xyzzz = pbuffer.data(idx_op_sh + 13);

    auto tr_0_xzzzz = pbuffer.data(idx_op_sh + 14);

    auto tr_0_yyyyy = pbuffer.data(idx_op_sh + 15);

    auto tr_0_yyyyz = pbuffer.data(idx_op_sh + 16);

    auto tr_0_yyyzz = pbuffer.data(idx_op_sh + 17);

    auto tr_0_yyzzz = pbuffer.data(idx_op_sh + 18);

    auto tr_0_yzzzz = pbuffer.data(idx_op_sh + 19);

    auto tr_0_zzzzz = pbuffer.data(idx_op_sh + 20);

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

    // Set up 0-15 components of targeted buffer : PG

    auto tr_x_0_x_x_xxxx = pbuffer.data(idx_op_geom_110_pg);

    auto tr_x_0_x_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 1);

    auto tr_x_0_x_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 2);

    auto tr_x_0_x_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 3);

    auto tr_x_0_x_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 4);

    auto tr_x_0_x_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 5);

    auto tr_x_0_x_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 6);

    auto tr_x_0_x_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 7);

    auto tr_x_0_x_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 8);

    auto tr_x_0_x_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 9);

    auto tr_x_0_x_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 10);

    auto tr_x_0_x_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 11);

    auto tr_x_0_x_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 12);

    auto tr_x_0_x_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 13);

    auto tr_x_0_x_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 14);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_0_x_x_xxxx, tr_x_0_x_x_xxxy, tr_x_0_x_x_xxxz, tr_x_0_x_x_xxyy, tr_x_0_x_x_xxyz, tr_x_0_x_x_xxzz, tr_x_0_x_x_xyyy, tr_x_0_x_x_xyyz, tr_x_0_x_x_xyzz, tr_x_0_x_x_xzzz, tr_x_0_x_x_yyyy, tr_x_0_x_x_yyyz, tr_x_0_x_x_yyzz, tr_x_0_x_x_yzzz, tr_x_0_x_x_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_x_xxxx[i] = 4.0 * tr_0_xxx[i] - 2.0 * tr_0_xxxxx[i] * tke_0 - 6.0 * tr_x_xxxx[i] * tbe_0 - 8.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxxy[i] = 3.0 * tr_0_xxy[i] - 2.0 * tr_0_xxxxy[i] * tke_0 - 6.0 * tr_x_xxxy[i] * tbe_0 - 6.0 * tr_xx_xxy[i] * tbe_0 + 4.0 * tr_xx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxxz[i] = 3.0 * tr_0_xxz[i] - 2.0 * tr_0_xxxxz[i] * tke_0 - 6.0 * tr_x_xxxz[i] * tbe_0 - 6.0 * tr_xx_xxz[i] * tbe_0 + 4.0 * tr_xx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxyy[i] = 2.0 * tr_0_xyy[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 6.0 * tr_x_xxyy[i] * tbe_0 - 4.0 * tr_xx_xyy[i] * tbe_0 + 4.0 * tr_xx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xxxyz[i] * tke_0 - 6.0 * tr_x_xxyz[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tbe_0 + 4.0 * tr_xx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xxzz[i] = 2.0 * tr_0_xzz[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 6.0 * tr_x_xxzz[i] * tbe_0 - 4.0 * tr_xx_xzz[i] * tbe_0 + 4.0 * tr_xx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xyyy[i] = tr_0_yyy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 6.0 * tr_x_xyyy[i] * tbe_0 - 2.0 * tr_xx_yyy[i] * tbe_0 + 4.0 * tr_xx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xyyz[i] = tr_0_yyz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 6.0 * tr_x_xyyz[i] * tbe_0 - 2.0 * tr_xx_yyz[i] * tbe_0 + 4.0 * tr_xx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xyzz[i] = tr_0_yzz[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 6.0 * tr_x_xyzz[i] * tbe_0 - 2.0 * tr_xx_yzz[i] * tbe_0 + 4.0 * tr_xx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_xzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 6.0 * tr_x_xzzz[i] * tbe_0 - 2.0 * tr_xx_zzz[i] * tbe_0 + 4.0 * tr_xx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yyyy[i] = -2.0 * tr_0_xyyyy[i] * tke_0 - 6.0 * tr_x_yyyy[i] * tbe_0 + 4.0 * tr_xx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yyyz[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 6.0 * tr_x_yyyz[i] * tbe_0 + 4.0 * tr_xx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yyzz[i] = -2.0 * tr_0_xyyzz[i] * tke_0 - 6.0 * tr_x_yyzz[i] * tbe_0 + 4.0 * tr_xx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_yzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_x_yzzz[i] * tbe_0 + 4.0 * tr_xx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_x_zzzz[i] = -2.0 * tr_0_xzzzz[i] * tke_0 - 6.0 * tr_x_zzzz[i] * tbe_0 + 4.0 * tr_xx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 15-30 components of targeted buffer : PG

    auto tr_x_0_x_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 15);

    auto tr_x_0_x_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 16);

    auto tr_x_0_x_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 17);

    auto tr_x_0_x_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 18);

    auto tr_x_0_x_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 19);

    auto tr_x_0_x_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 20);

    auto tr_x_0_x_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 21);

    auto tr_x_0_x_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 22);

    auto tr_x_0_x_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 23);

    auto tr_x_0_x_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 24);

    auto tr_x_0_x_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 25);

    auto tr_x_0_x_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 26);

    auto tr_x_0_x_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 27);

    auto tr_x_0_x_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 28);

    auto tr_x_0_x_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 29);

    #pragma omp simd aligned(tr_x_0_x_y_xxxx, tr_x_0_x_y_xxxy, tr_x_0_x_y_xxxz, tr_x_0_x_y_xxyy, tr_x_0_x_y_xxyz, tr_x_0_x_y_xxzz, tr_x_0_x_y_xyyy, tr_x_0_x_y_xyyz, tr_x_0_x_y_xyzz, tr_x_0_x_y_xzzz, tr_x_0_x_y_yyyy, tr_x_0_x_y_yyyz, tr_x_0_x_y_yyzz, tr_x_0_x_y_yzzz, tr_x_0_x_y_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_y_xxxx[i] = -2.0 * tr_y_xxxx[i] * tbe_0 - 8.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxxy[i] = -2.0 * tr_y_xxxy[i] * tbe_0 - 6.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxxz[i] = -2.0 * tr_y_xxxz[i] * tbe_0 - 6.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxyy[i] = -2.0 * tr_y_xxyy[i] * tbe_0 - 4.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxyz[i] = -2.0 * tr_y_xxyz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xxzz[i] = -2.0 * tr_y_xxzz[i] * tbe_0 - 4.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xyyy[i] = -2.0 * tr_y_xyyy[i] * tbe_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xyyz[i] = -2.0 * tr_y_xyyz[i] * tbe_0 - 2.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xyzz[i] = -2.0 * tr_y_xyzz[i] * tbe_0 - 2.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_xzzz[i] = -2.0 * tr_y_xzzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yyyy[i] = -2.0 * tr_y_yyyy[i] * tbe_0 + 4.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yyyz[i] = -2.0 * tr_y_yyyz[i] * tbe_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yyzz[i] = -2.0 * tr_y_yyzz[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_yzzz[i] = -2.0 * tr_y_yzzz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_y_zzzz[i] = -2.0 * tr_y_zzzz[i] * tbe_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-45 components of targeted buffer : PG

    auto tr_x_0_x_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 30);

    auto tr_x_0_x_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 31);

    auto tr_x_0_x_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 32);

    auto tr_x_0_x_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 33);

    auto tr_x_0_x_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 34);

    auto tr_x_0_x_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 35);

    auto tr_x_0_x_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 36);

    auto tr_x_0_x_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 37);

    auto tr_x_0_x_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 38);

    auto tr_x_0_x_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 39);

    auto tr_x_0_x_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 40);

    auto tr_x_0_x_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 41);

    auto tr_x_0_x_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 42);

    auto tr_x_0_x_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 43);

    auto tr_x_0_x_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 44);

    #pragma omp simd aligned(tr_x_0_x_z_xxxx, tr_x_0_x_z_xxxy, tr_x_0_x_z_xxxz, tr_x_0_x_z_xxyy, tr_x_0_x_z_xxyz, tr_x_0_x_z_xxzz, tr_x_0_x_z_xyyy, tr_x_0_x_z_xyyz, tr_x_0_x_z_xyzz, tr_x_0_x_z_xzzz, tr_x_0_x_z_yyyy, tr_x_0_x_z_yyyz, tr_x_0_x_z_yyzz, tr_x_0_x_z_yzzz, tr_x_0_x_z_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_z_xxxx[i] = -2.0 * tr_z_xxxx[i] * tbe_0 - 8.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxxy[i] = -2.0 * tr_z_xxxy[i] * tbe_0 - 6.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxxz[i] = -2.0 * tr_z_xxxz[i] * tbe_0 - 6.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxyy[i] = -2.0 * tr_z_xxyy[i] * tbe_0 - 4.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxyz[i] = -2.0 * tr_z_xxyz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xxzz[i] = -2.0 * tr_z_xxzz[i] * tbe_0 - 4.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xyyy[i] = -2.0 * tr_z_xyyy[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xyyz[i] = -2.0 * tr_z_xyyz[i] * tbe_0 - 2.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xyzz[i] = -2.0 * tr_z_xyzz[i] * tbe_0 - 2.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_xzzz[i] = -2.0 * tr_z_xzzz[i] * tbe_0 - 2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yyyy[i] = -2.0 * tr_z_yyyy[i] * tbe_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yyyz[i] = -2.0 * tr_z_yyyz[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yyzz[i] = -2.0 * tr_z_yyzz[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_yzzz[i] = -2.0 * tr_z_yzzz[i] * tbe_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_z_zzzz[i] = -2.0 * tr_z_zzzz[i] * tbe_0 + 4.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 45-60 components of targeted buffer : PG

    auto tr_x_0_y_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 45);

    auto tr_x_0_y_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 46);

    auto tr_x_0_y_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 47);

    auto tr_x_0_y_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 48);

    auto tr_x_0_y_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 49);

    auto tr_x_0_y_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 50);

    auto tr_x_0_y_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 51);

    auto tr_x_0_y_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 52);

    auto tr_x_0_y_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 53);

    auto tr_x_0_y_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 54);

    auto tr_x_0_y_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 55);

    auto tr_x_0_y_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 56);

    auto tr_x_0_y_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 57);

    auto tr_x_0_y_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 58);

    auto tr_x_0_y_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 59);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxy, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_x_0_y_x_xxxx, tr_x_0_y_x_xxxy, tr_x_0_y_x_xxxz, tr_x_0_y_x_xxyy, tr_x_0_y_x_xxyz, tr_x_0_y_x_xxzz, tr_x_0_y_x_xyyy, tr_x_0_y_x_xyyz, tr_x_0_y_x_xyzz, tr_x_0_y_x_xzzz, tr_x_0_y_x_yyyy, tr_x_0_y_x_yyyz, tr_x_0_y_x_yyzz, tr_x_0_y_x_yzzz, tr_x_0_y_x_zzzz, tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_x_xxxx[i] = -2.0 * tr_0_xxxxy[i] * tke_0 - 2.0 * tr_y_xxxx[i] * tbe_0 + 4.0 * tr_xx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxxy[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 2.0 * tr_y_xxxy[i] * tbe_0 - 2.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxxz[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_y_xxxz[i] * tbe_0 + 4.0 * tr_xx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxyy[i] = 2.0 * tr_0_xxy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 2.0 * tr_y_xxyy[i] * tbe_0 - 4.0 * tr_xx_xxy[i] * tbe_0 + 4.0 * tr_xx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxyz[i] = tr_0_xxz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_y_xxyz[i] * tbe_0 - 2.0 * tr_xx_xxz[i] * tbe_0 + 4.0 * tr_xx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xxzz[i] = -2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_y_xxzz[i] * tbe_0 + 4.0 * tr_xx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xyyy[i] = 3.0 * tr_0_xyy[i] - 2.0 * tr_0_xyyyy[i] * tke_0 - 2.0 * tr_y_xyyy[i] * tbe_0 - 6.0 * tr_xx_xyy[i] * tbe_0 + 4.0 * tr_xx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xyyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_y_xyyz[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tbe_0 + 4.0 * tr_xx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xyzz[i] = tr_0_xzz[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_y_xyzz[i] * tbe_0 - 2.0 * tr_xx_xzz[i] * tbe_0 + 4.0 * tr_xx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_xzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_y_xzzz[i] * tbe_0 + 4.0 * tr_xx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yyyy[i] = 4.0 * tr_0_yyy[i] - 2.0 * tr_0_yyyyy[i] * tke_0 - 2.0 * tr_y_yyyy[i] * tbe_0 - 8.0 * tr_xx_yyy[i] * tbe_0 + 4.0 * tr_xx_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yyyz[i] = 3.0 * tr_0_yyz[i] - 2.0 * tr_0_yyyyz[i] * tke_0 - 2.0 * tr_y_yyyz[i] * tbe_0 - 6.0 * tr_xx_yyz[i] * tbe_0 + 4.0 * tr_xx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yyzz[i] = 2.0 * tr_0_yzz[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_y_yyzz[i] * tbe_0 - 4.0 * tr_xx_yzz[i] * tbe_0 + 4.0 * tr_xx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_yzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_y_yzzz[i] * tbe_0 - 2.0 * tr_xx_zzz[i] * tbe_0 + 4.0 * tr_xx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_x_zzzz[i] = -2.0 * tr_0_yzzzz[i] * tke_0 - 2.0 * tr_y_zzzz[i] * tbe_0 + 4.0 * tr_xx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-75 components of targeted buffer : PG

    auto tr_x_0_y_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 60);

    auto tr_x_0_y_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 61);

    auto tr_x_0_y_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 62);

    auto tr_x_0_y_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 63);

    auto tr_x_0_y_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 64);

    auto tr_x_0_y_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 65);

    auto tr_x_0_y_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 66);

    auto tr_x_0_y_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 67);

    auto tr_x_0_y_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 68);

    auto tr_x_0_y_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 69);

    auto tr_x_0_y_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 70);

    auto tr_x_0_y_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 71);

    auto tr_x_0_y_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 72);

    auto tr_x_0_y_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 73);

    auto tr_x_0_y_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 74);

    #pragma omp simd aligned(tr_x_0_y_y_xxxx, tr_x_0_y_y_xxxy, tr_x_0_y_y_xxxz, tr_x_0_y_y_xxyy, tr_x_0_y_y_xxyz, tr_x_0_y_y_xxzz, tr_x_0_y_y_xyyy, tr_x_0_y_y_xyyz, tr_x_0_y_y_xyzz, tr_x_0_y_y_xzzz, tr_x_0_y_y_yyyy, tr_x_0_y_y_yyyz, tr_x_0_y_y_yyzz, tr_x_0_y_y_yzzz, tr_x_0_y_y_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_y_xxxx[i] = -2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxxy[i] = -2.0 * tr_x_xxxy[i] * tbe_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxxz[i] = -2.0 * tr_x_xxxz[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxyy[i] = -2.0 * tr_x_xxyy[i] * tbe_0 - 4.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxyz[i] = -2.0 * tr_x_xxyz[i] * tbe_0 - 2.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xxzz[i] = -2.0 * tr_x_xxzz[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xyyy[i] = -2.0 * tr_x_xyyy[i] * tbe_0 - 6.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xyyz[i] = -2.0 * tr_x_xyyz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xyzz[i] = -2.0 * tr_x_xyzz[i] * tbe_0 - 2.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_xzzz[i] = -2.0 * tr_x_xzzz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yyyy[i] = -2.0 * tr_x_yyyy[i] * tbe_0 - 8.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yyyz[i] = -2.0 * tr_x_yyyz[i] * tbe_0 - 6.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yyzz[i] = -2.0 * tr_x_yyzz[i] * tbe_0 - 4.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_yzzz[i] = -2.0 * tr_x_yzzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_y_zzzz[i] = -2.0 * tr_x_zzzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 75-90 components of targeted buffer : PG

    auto tr_x_0_y_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 75);

    auto tr_x_0_y_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 76);

    auto tr_x_0_y_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 77);

    auto tr_x_0_y_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 78);

    auto tr_x_0_y_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 79);

    auto tr_x_0_y_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 80);

    auto tr_x_0_y_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 81);

    auto tr_x_0_y_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 82);

    auto tr_x_0_y_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 83);

    auto tr_x_0_y_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 84);

    auto tr_x_0_y_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 85);

    auto tr_x_0_y_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 86);

    auto tr_x_0_y_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 87);

    auto tr_x_0_y_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 88);

    auto tr_x_0_y_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 89);

    #pragma omp simd aligned(tr_x_0_y_z_xxxx, tr_x_0_y_z_xxxy, tr_x_0_y_z_xxxz, tr_x_0_y_z_xxyy, tr_x_0_y_z_xxyz, tr_x_0_y_z_xxzz, tr_x_0_y_z_xyyy, tr_x_0_y_z_xyyz, tr_x_0_y_z_xyzz, tr_x_0_y_z_xzzz, tr_x_0_y_z_yyyy, tr_x_0_y_z_yyyz, tr_x_0_y_z_yyzz, tr_x_0_y_z_yzzz, tr_x_0_y_z_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_z_xxxx[i] = 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxxy[i] = -2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxxz[i] = 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxyy[i] = -4.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxyz[i] = -2.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xxzz[i] = 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xyyy[i] = -6.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xyyz[i] = -4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xyzz[i] = -2.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_xzzz[i] = 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yyyy[i] = -8.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yyyz[i] = -6.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yyzz[i] = -4.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_yzzz[i] = -2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_z_zzzz[i] = 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-105 components of targeted buffer : PG

    auto tr_x_0_z_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 90);

    auto tr_x_0_z_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 91);

    auto tr_x_0_z_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 92);

    auto tr_x_0_z_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 93);

    auto tr_x_0_z_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 94);

    auto tr_x_0_z_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 95);

    auto tr_x_0_z_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 96);

    auto tr_x_0_z_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 97);

    auto tr_x_0_z_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 98);

    auto tr_x_0_z_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 99);

    auto tr_x_0_z_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 100);

    auto tr_x_0_z_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 101);

    auto tr_x_0_z_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 102);

    auto tr_x_0_z_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 103);

    auto tr_x_0_z_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 104);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxz, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_x_0_z_x_xxxx, tr_x_0_z_x_xxxy, tr_x_0_z_x_xxxz, tr_x_0_z_x_xxyy, tr_x_0_z_x_xxyz, tr_x_0_z_x_xxzz, tr_x_0_z_x_xyyy, tr_x_0_z_x_xyyz, tr_x_0_z_x_xyzz, tr_x_0_z_x_xzzz, tr_x_0_z_x_yyyy, tr_x_0_z_x_yyyz, tr_x_0_z_x_yyzz, tr_x_0_z_x_yzzz, tr_x_0_z_x_zzzz, tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xx_zzzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_x_xxxx[i] = -2.0 * tr_0_xxxxz[i] * tke_0 - 2.0 * tr_z_xxxx[i] * tbe_0 + 4.0 * tr_xx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxxy[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_z_xxxy[i] * tbe_0 + 4.0 * tr_xx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxxz[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 2.0 * tr_z_xxxz[i] * tbe_0 - 2.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxyy[i] = -2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_z_xxyy[i] * tbe_0 + 4.0 * tr_xx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxyz[i] = tr_0_xxy[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_z_xxyz[i] * tbe_0 - 2.0 * tr_xx_xxy[i] * tbe_0 + 4.0 * tr_xx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xxzz[i] = 2.0 * tr_0_xxz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_z_xxzz[i] * tbe_0 - 4.0 * tr_xx_xxz[i] * tbe_0 + 4.0 * tr_xx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xyyy[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_z_xyyy[i] * tbe_0 + 4.0 * tr_xx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xyyz[i] = tr_0_xyy[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_z_xyyz[i] * tbe_0 - 2.0 * tr_xx_xyy[i] * tbe_0 + 4.0 * tr_xx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xyzz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_z_xyzz[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tbe_0 + 4.0 * tr_xx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_xzzz[i] = 3.0 * tr_0_xzz[i] - 2.0 * tr_0_xzzzz[i] * tke_0 - 2.0 * tr_z_xzzz[i] * tbe_0 - 6.0 * tr_xx_xzz[i] * tbe_0 + 4.0 * tr_xx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yyyy[i] = -2.0 * tr_0_yyyyz[i] * tke_0 - 2.0 * tr_z_yyyy[i] * tbe_0 + 4.0 * tr_xx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yyyz[i] = tr_0_yyy[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_z_yyyz[i] * tbe_0 - 2.0 * tr_xx_yyy[i] * tbe_0 + 4.0 * tr_xx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yyzz[i] = 2.0 * tr_0_yyz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_z_yyzz[i] * tbe_0 - 4.0 * tr_xx_yyz[i] * tbe_0 + 4.0 * tr_xx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_yzzz[i] = 3.0 * tr_0_yzz[i] - 2.0 * tr_0_yzzzz[i] * tke_0 - 2.0 * tr_z_yzzz[i] * tbe_0 - 6.0 * tr_xx_yzz[i] * tbe_0 + 4.0 * tr_xx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_x_zzzz[i] = 4.0 * tr_0_zzz[i] - 2.0 * tr_0_zzzzz[i] * tke_0 - 2.0 * tr_z_zzzz[i] * tbe_0 - 8.0 * tr_xx_zzz[i] * tbe_0 + 4.0 * tr_xx_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 105-120 components of targeted buffer : PG

    auto tr_x_0_z_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 105);

    auto tr_x_0_z_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 106);

    auto tr_x_0_z_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 107);

    auto tr_x_0_z_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 108);

    auto tr_x_0_z_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 109);

    auto tr_x_0_z_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 110);

    auto tr_x_0_z_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 111);

    auto tr_x_0_z_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 112);

    auto tr_x_0_z_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 113);

    auto tr_x_0_z_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 114);

    auto tr_x_0_z_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 115);

    auto tr_x_0_z_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 116);

    auto tr_x_0_z_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 117);

    auto tr_x_0_z_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 118);

    auto tr_x_0_z_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 119);

    #pragma omp simd aligned(tr_x_0_z_y_xxxx, tr_x_0_z_y_xxxy, tr_x_0_z_y_xxxz, tr_x_0_z_y_xxyy, tr_x_0_z_y_xxyz, tr_x_0_z_y_xxzz, tr_x_0_z_y_xyyy, tr_x_0_z_y_xyyz, tr_x_0_z_y_xyzz, tr_x_0_z_y_xzzz, tr_x_0_z_y_yyyy, tr_x_0_z_y_yyyz, tr_x_0_z_y_yyzz, tr_x_0_z_y_yzzz, tr_x_0_z_y_zzzz, tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_y_xxxx[i] = 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxxy[i] = 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxxz[i] = -2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxyy[i] = 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxyz[i] = -2.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xxzz[i] = -4.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xyyy[i] = 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xyyz[i] = -2.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xyzz[i] = -4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_xzzz[i] = -6.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yyyy[i] = 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yyyz[i] = -2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yyzz[i] = -4.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_yzzz[i] = -6.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_y_zzzz[i] = -8.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-135 components of targeted buffer : PG

    auto tr_x_0_z_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 120);

    auto tr_x_0_z_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 121);

    auto tr_x_0_z_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 122);

    auto tr_x_0_z_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 123);

    auto tr_x_0_z_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 124);

    auto tr_x_0_z_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 125);

    auto tr_x_0_z_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 126);

    auto tr_x_0_z_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 127);

    auto tr_x_0_z_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 128);

    auto tr_x_0_z_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 129);

    auto tr_x_0_z_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 130);

    auto tr_x_0_z_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 131);

    auto tr_x_0_z_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 132);

    auto tr_x_0_z_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 133);

    auto tr_x_0_z_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 134);

    #pragma omp simd aligned(tr_x_0_z_z_xxxx, tr_x_0_z_z_xxxy, tr_x_0_z_z_xxxz, tr_x_0_z_z_xxyy, tr_x_0_z_z_xxyz, tr_x_0_z_z_xxzz, tr_x_0_z_z_xyyy, tr_x_0_z_z_xyyz, tr_x_0_z_z_xyzz, tr_x_0_z_z_xzzz, tr_x_0_z_z_yyyy, tr_x_0_z_z_yyyz, tr_x_0_z_z_yyzz, tr_x_0_z_z_yzzz, tr_x_0_z_z_zzzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_z_xxxx[i] = -2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxxy[i] = -2.0 * tr_x_xxxy[i] * tbe_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxxz[i] = -2.0 * tr_x_xxxz[i] * tbe_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxyy[i] = -2.0 * tr_x_xxyy[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxyz[i] = -2.0 * tr_x_xxyz[i] * tbe_0 - 2.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xxzz[i] = -2.0 * tr_x_xxzz[i] * tbe_0 - 4.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xyyy[i] = -2.0 * tr_x_xyyy[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xyyz[i] = -2.0 * tr_x_xyyz[i] * tbe_0 - 2.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xyzz[i] = -2.0 * tr_x_xyzz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_xzzz[i] = -2.0 * tr_x_xzzz[i] * tbe_0 - 6.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yyyy[i] = -2.0 * tr_x_yyyy[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yyyz[i] = -2.0 * tr_x_yyyz[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yyzz[i] = -2.0 * tr_x_yyzz[i] * tbe_0 - 4.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_yzzz[i] = -2.0 * tr_x_yzzz[i] * tbe_0 - 6.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_z_zzzz[i] = -2.0 * tr_x_zzzz[i] * tbe_0 - 8.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 135-150 components of targeted buffer : PG

    auto tr_y_0_x_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 135);

    auto tr_y_0_x_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 136);

    auto tr_y_0_x_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 137);

    auto tr_y_0_x_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 138);

    auto tr_y_0_x_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 139);

    auto tr_y_0_x_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 140);

    auto tr_y_0_x_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 141);

    auto tr_y_0_x_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 142);

    auto tr_y_0_x_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 143);

    auto tr_y_0_x_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 144);

    auto tr_y_0_x_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 145);

    auto tr_y_0_x_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 146);

    auto tr_y_0_x_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 147);

    auto tr_y_0_x_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 148);

    auto tr_y_0_x_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 149);

    #pragma omp simd aligned(tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_0_x_x_xxxx, tr_y_0_x_x_xxxy, tr_y_0_x_x_xxxz, tr_y_0_x_x_xxyy, tr_y_0_x_x_xxyz, tr_y_0_x_x_xxzz, tr_y_0_x_x_xyyy, tr_y_0_x_x_xyyz, tr_y_0_x_x_xyzz, tr_y_0_x_x_xzzz, tr_y_0_x_x_yyyy, tr_y_0_x_x_yyyz, tr_y_0_x_x_yyzz, tr_y_0_x_x_yzzz, tr_y_0_x_x_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_x_xxxx[i] = -2.0 * tr_y_xxxx[i] * tbe_0 - 8.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxxy[i] = -2.0 * tr_y_xxxy[i] * tbe_0 - 6.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxxz[i] = -2.0 * tr_y_xxxz[i] * tbe_0 - 6.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxyy[i] = -2.0 * tr_y_xxyy[i] * tbe_0 - 4.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxyz[i] = -2.0 * tr_y_xxyz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xxzz[i] = -2.0 * tr_y_xxzz[i] * tbe_0 - 4.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xyyy[i] = -2.0 * tr_y_xyyy[i] * tbe_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xyyz[i] = -2.0 * tr_y_xyyz[i] * tbe_0 - 2.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xyzz[i] = -2.0 * tr_y_xyzz[i] * tbe_0 - 2.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_xzzz[i] = -2.0 * tr_y_xzzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yyyy[i] = -2.0 * tr_y_yyyy[i] * tbe_0 + 4.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yyyz[i] = -2.0 * tr_y_yyyz[i] * tbe_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yyzz[i] = -2.0 * tr_y_yyzz[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_yzzz[i] = -2.0 * tr_y_yzzz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_x_zzzz[i] = -2.0 * tr_y_zzzz[i] * tbe_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-165 components of targeted buffer : PG

    auto tr_y_0_x_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 150);

    auto tr_y_0_x_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 151);

    auto tr_y_0_x_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 152);

    auto tr_y_0_x_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 153);

    auto tr_y_0_x_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 154);

    auto tr_y_0_x_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 155);

    auto tr_y_0_x_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 156);

    auto tr_y_0_x_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 157);

    auto tr_y_0_x_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 158);

    auto tr_y_0_x_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 159);

    auto tr_y_0_x_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 160);

    auto tr_y_0_x_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 161);

    auto tr_y_0_x_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 162);

    auto tr_y_0_x_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 163);

    auto tr_y_0_x_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 164);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_y_0_x_y_xxxx, tr_y_0_x_y_xxxy, tr_y_0_x_y_xxxz, tr_y_0_x_y_xxyy, tr_y_0_x_y_xxyz, tr_y_0_x_y_xxzz, tr_y_0_x_y_xyyy, tr_y_0_x_y_xyyz, tr_y_0_x_y_xyzz, tr_y_0_x_y_xzzz, tr_y_0_x_y_yyyy, tr_y_0_x_y_yyyz, tr_y_0_x_y_yyzz, tr_y_0_x_y_yzzz, tr_y_0_x_y_zzzz, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_y_xxxx[i] = 4.0 * tr_0_xxx[i] - 2.0 * tr_0_xxxxx[i] * tke_0 - 8.0 * tr_yy_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxxy[i] = 3.0 * tr_0_xxy[i] - 2.0 * tr_0_xxxxy[i] * tke_0 - 6.0 * tr_yy_xxy[i] * tbe_0 + 4.0 * tr_yy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxy[i] * tbe_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxxz[i] = 3.0 * tr_0_xxz[i] - 2.0 * tr_0_xxxxz[i] * tke_0 - 6.0 * tr_yy_xxz[i] * tbe_0 + 4.0 * tr_yy_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxz[i] * tbe_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxyy[i] = 2.0 * tr_0_xyy[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 4.0 * tr_yy_xyy[i] * tbe_0 + 4.0 * tr_yy_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyy[i] * tbe_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xxxyz[i] * tke_0 - 4.0 * tr_yy_xyz[i] * tbe_0 + 4.0 * tr_yy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyz[i] * tbe_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xxzz[i] = 2.0 * tr_0_xzz[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 4.0 * tr_yy_xzz[i] * tbe_0 + 4.0 * tr_yy_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxzz[i] * tbe_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xyyy[i] = tr_0_yyy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 2.0 * tr_yy_yyy[i] * tbe_0 + 4.0 * tr_yy_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyy[i] * tbe_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xyyz[i] = tr_0_yyz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_yy_yyz[i] * tbe_0 + 4.0 * tr_yy_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyz[i] * tbe_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xyzz[i] = tr_0_yzz[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_yy_yzz[i] * tbe_0 + 4.0 * tr_yy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyzz[i] * tbe_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_xzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_yy_zzz[i] * tbe_0 + 4.0 * tr_yy_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzzz[i] * tbe_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yyyy[i] = -2.0 * tr_0_xyyyy[i] * tke_0 + 4.0 * tr_yy_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyy[i] * tbe_0 + 4.0 * tr_xyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yyyz[i] = -2.0 * tr_0_xyyyz[i] * tke_0 + 4.0 * tr_yy_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyz[i] * tbe_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yyzz[i] = -2.0 * tr_0_xyyzz[i] * tke_0 + 4.0 * tr_yy_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyzz[i] * tbe_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_yzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 + 4.0 * tr_yy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzzz[i] * tbe_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_y_zzzz[i] = -2.0 * tr_0_xzzzz[i] * tke_0 + 4.0 * tr_yy_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzzz[i] * tbe_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 165-180 components of targeted buffer : PG

    auto tr_y_0_x_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 165);

    auto tr_y_0_x_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 166);

    auto tr_y_0_x_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 167);

    auto tr_y_0_x_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 168);

    auto tr_y_0_x_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 169);

    auto tr_y_0_x_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 170);

    auto tr_y_0_x_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 171);

    auto tr_y_0_x_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 172);

    auto tr_y_0_x_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 173);

    auto tr_y_0_x_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 174);

    auto tr_y_0_x_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 175);

    auto tr_y_0_x_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 176);

    auto tr_y_0_x_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 177);

    auto tr_y_0_x_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 178);

    auto tr_y_0_x_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 179);

    #pragma omp simd aligned(tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_y_0_x_z_xxxx, tr_y_0_x_z_xxxy, tr_y_0_x_z_xxxz, tr_y_0_x_z_xxyy, tr_y_0_x_z_xxyz, tr_y_0_x_z_xxzz, tr_y_0_x_z_xyyy, tr_y_0_x_z_xyyz, tr_y_0_x_z_xyzz, tr_y_0_x_z_xzzz, tr_y_0_x_z_yyyy, tr_y_0_x_z_yyyz, tr_y_0_x_z_yyzz, tr_y_0_x_z_yzzz, tr_y_0_x_z_zzzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_z_xxxx[i] = -8.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxxy[i] = -6.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxxz[i] = -6.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxyy[i] = -4.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxyz[i] = -4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xxzz[i] = -4.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xyyy[i] = -2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xyyz[i] = -2.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xyzz[i] = -2.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_xzzz[i] = -2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yyyy[i] = 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yyyz[i] = 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yyzz[i] = 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_yzzz[i] = 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_z_zzzz[i] = 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-195 components of targeted buffer : PG

    auto tr_y_0_y_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 180);

    auto tr_y_0_y_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 181);

    auto tr_y_0_y_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 182);

    auto tr_y_0_y_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 183);

    auto tr_y_0_y_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 184);

    auto tr_y_0_y_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 185);

    auto tr_y_0_y_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 186);

    auto tr_y_0_y_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 187);

    auto tr_y_0_y_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 188);

    auto tr_y_0_y_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 189);

    auto tr_y_0_y_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 190);

    auto tr_y_0_y_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 191);

    auto tr_y_0_y_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 192);

    auto tr_y_0_y_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 193);

    auto tr_y_0_y_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 194);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_y_0_y_x_xxxx, tr_y_0_y_x_xxxy, tr_y_0_y_x_xxxz, tr_y_0_y_x_xxyy, tr_y_0_y_x_xxyz, tr_y_0_y_x_xxzz, tr_y_0_y_x_xyyy, tr_y_0_y_x_xyyz, tr_y_0_y_x_xyzz, tr_y_0_y_x_xzzz, tr_y_0_y_x_yyyy, tr_y_0_y_x_yyyz, tr_y_0_y_x_yyzz, tr_y_0_y_x_yzzz, tr_y_0_y_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_x_xxxx[i] = -2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxxy[i] = -2.0 * tr_x_xxxy[i] * tbe_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxxz[i] = -2.0 * tr_x_xxxz[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxyy[i] = -2.0 * tr_x_xxyy[i] * tbe_0 - 4.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxyz[i] = -2.0 * tr_x_xxyz[i] * tbe_0 - 2.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xxzz[i] = -2.0 * tr_x_xxzz[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xyyy[i] = -2.0 * tr_x_xyyy[i] * tbe_0 - 6.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xyyz[i] = -2.0 * tr_x_xyyz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xyzz[i] = -2.0 * tr_x_xyzz[i] * tbe_0 - 2.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_xzzz[i] = -2.0 * tr_x_xzzz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yyyy[i] = -2.0 * tr_x_yyyy[i] * tbe_0 - 8.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yyyz[i] = -2.0 * tr_x_yyyz[i] * tbe_0 - 6.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yyzz[i] = -2.0 * tr_x_yyzz[i] * tbe_0 - 4.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_yzzz[i] = -2.0 * tr_x_yzzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_x_zzzz[i] = -2.0 * tr_x_zzzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 195-210 components of targeted buffer : PG

    auto tr_y_0_y_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 195);

    auto tr_y_0_y_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 196);

    auto tr_y_0_y_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 197);

    auto tr_y_0_y_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 198);

    auto tr_y_0_y_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 199);

    auto tr_y_0_y_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 200);

    auto tr_y_0_y_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 201);

    auto tr_y_0_y_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 202);

    auto tr_y_0_y_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 203);

    auto tr_y_0_y_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 204);

    auto tr_y_0_y_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 205);

    auto tr_y_0_y_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 206);

    auto tr_y_0_y_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 207);

    auto tr_y_0_y_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 208);

    auto tr_y_0_y_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 209);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxy, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_y_0_y_y_xxxx, tr_y_0_y_y_xxxy, tr_y_0_y_y_xxxz, tr_y_0_y_y_xxyy, tr_y_0_y_y_xxyz, tr_y_0_y_y_xxzz, tr_y_0_y_y_xyyy, tr_y_0_y_y_xyyz, tr_y_0_y_y_xyzz, tr_y_0_y_y_xzzz, tr_y_0_y_y_yyyy, tr_y_0_y_y_yyyz, tr_y_0_y_y_yyzz, tr_y_0_y_y_yzzz, tr_y_0_y_y_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_y_xxxx[i] = -2.0 * tr_0_xxxxy[i] * tke_0 - 6.0 * tr_y_xxxx[i] * tbe_0 + 4.0 * tr_yy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxxy[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 6.0 * tr_y_xxxy[i] * tbe_0 - 2.0 * tr_yy_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxxz[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 6.0 * tr_y_xxxz[i] * tbe_0 + 4.0 * tr_yy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxyy[i] = 2.0 * tr_0_xxy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 6.0 * tr_y_xxyy[i] * tbe_0 - 4.0 * tr_yy_xxy[i] * tbe_0 + 4.0 * tr_yy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxyz[i] = tr_0_xxz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 6.0 * tr_y_xxyz[i] * tbe_0 - 2.0 * tr_yy_xxz[i] * tbe_0 + 4.0 * tr_yy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xxzz[i] = -2.0 * tr_0_xxyzz[i] * tke_0 - 6.0 * tr_y_xxzz[i] * tbe_0 + 4.0 * tr_yy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xyyy[i] = 3.0 * tr_0_xyy[i] - 2.0 * tr_0_xyyyy[i] * tke_0 - 6.0 * tr_y_xyyy[i] * tbe_0 - 6.0 * tr_yy_xyy[i] * tbe_0 + 4.0 * tr_yy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xyyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyyyz[i] * tke_0 - 6.0 * tr_y_xyyz[i] * tbe_0 - 4.0 * tr_yy_xyz[i] * tbe_0 + 4.0 * tr_yy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xyzz[i] = tr_0_xzz[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 6.0 * tr_y_xyzz[i] * tbe_0 - 2.0 * tr_yy_xzz[i] * tbe_0 + 4.0 * tr_yy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_xzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_y_xzzz[i] * tbe_0 + 4.0 * tr_yy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yyyy[i] = 4.0 * tr_0_yyy[i] - 2.0 * tr_0_yyyyy[i] * tke_0 - 6.0 * tr_y_yyyy[i] * tbe_0 - 8.0 * tr_yy_yyy[i] * tbe_0 + 4.0 * tr_yy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yyyz[i] = 3.0 * tr_0_yyz[i] - 2.0 * tr_0_yyyyz[i] * tke_0 - 6.0 * tr_y_yyyz[i] * tbe_0 - 6.0 * tr_yy_yyz[i] * tbe_0 + 4.0 * tr_yy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yyzz[i] = 2.0 * tr_0_yzz[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 6.0 * tr_y_yyzz[i] * tbe_0 - 4.0 * tr_yy_yzz[i] * tbe_0 + 4.0 * tr_yy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_yzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 6.0 * tr_y_yzzz[i] * tbe_0 - 2.0 * tr_yy_zzz[i] * tbe_0 + 4.0 * tr_yy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_y_zzzz[i] = -2.0 * tr_0_yzzzz[i] * tke_0 - 6.0 * tr_y_zzzz[i] * tbe_0 + 4.0 * tr_yy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-225 components of targeted buffer : PG

    auto tr_y_0_y_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 210);

    auto tr_y_0_y_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 211);

    auto tr_y_0_y_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 212);

    auto tr_y_0_y_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 213);

    auto tr_y_0_y_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 214);

    auto tr_y_0_y_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 215);

    auto tr_y_0_y_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 216);

    auto tr_y_0_y_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 217);

    auto tr_y_0_y_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 218);

    auto tr_y_0_y_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 219);

    auto tr_y_0_y_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 220);

    auto tr_y_0_y_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 221);

    auto tr_y_0_y_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 222);

    auto tr_y_0_y_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 223);

    auto tr_y_0_y_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 224);

    #pragma omp simd aligned(tr_y_0_y_z_xxxx, tr_y_0_y_z_xxxy, tr_y_0_y_z_xxxz, tr_y_0_y_z_xxyy, tr_y_0_y_z_xxyz, tr_y_0_y_z_xxzz, tr_y_0_y_z_xyyy, tr_y_0_y_z_xyyz, tr_y_0_y_z_xyzz, tr_y_0_y_z_xzzz, tr_y_0_y_z_yyyy, tr_y_0_y_z_yyyz, tr_y_0_y_z_yyzz, tr_y_0_y_z_yzzz, tr_y_0_y_z_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_z_xxxx[i] = -2.0 * tr_z_xxxx[i] * tbe_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxxy[i] = -2.0 * tr_z_xxxy[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxxz[i] = -2.0 * tr_z_xxxz[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxyy[i] = -2.0 * tr_z_xxyy[i] * tbe_0 - 4.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxyz[i] = -2.0 * tr_z_xxyz[i] * tbe_0 - 2.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xxzz[i] = -2.0 * tr_z_xxzz[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xyyy[i] = -2.0 * tr_z_xyyy[i] * tbe_0 - 6.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xyyz[i] = -2.0 * tr_z_xyyz[i] * tbe_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xyzz[i] = -2.0 * tr_z_xyzz[i] * tbe_0 - 2.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_xzzz[i] = -2.0 * tr_z_xzzz[i] * tbe_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yyyy[i] = -2.0 * tr_z_yyyy[i] * tbe_0 - 8.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yyyz[i] = -2.0 * tr_z_yyyz[i] * tbe_0 - 6.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yyzz[i] = -2.0 * tr_z_yyzz[i] * tbe_0 - 4.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_yzzz[i] = -2.0 * tr_z_yzzz[i] * tbe_0 - 2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_z_zzzz[i] = -2.0 * tr_z_zzzz[i] * tbe_0 + 4.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 225-240 components of targeted buffer : PG

    auto tr_y_0_z_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 225);

    auto tr_y_0_z_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 226);

    auto tr_y_0_z_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 227);

    auto tr_y_0_z_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 228);

    auto tr_y_0_z_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 229);

    auto tr_y_0_z_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 230);

    auto tr_y_0_z_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 231);

    auto tr_y_0_z_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 232);

    auto tr_y_0_z_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 233);

    auto tr_y_0_z_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 234);

    auto tr_y_0_z_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 235);

    auto tr_y_0_z_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 236);

    auto tr_y_0_z_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 237);

    auto tr_y_0_z_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 238);

    auto tr_y_0_z_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 239);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_y_0_z_x_xxxx, tr_y_0_z_x_xxxy, tr_y_0_z_x_xxxz, tr_y_0_z_x_xxyy, tr_y_0_z_x_xxyz, tr_y_0_z_x_xxzz, tr_y_0_z_x_xyyy, tr_y_0_z_x_xyyz, tr_y_0_z_x_xyzz, tr_y_0_z_x_xzzz, tr_y_0_z_x_yyyy, tr_y_0_z_x_yyyz, tr_y_0_z_x_yyzz, tr_y_0_z_x_yzzz, tr_y_0_z_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_x_xxxx[i] = 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxxy[i] = 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxxz[i] = -2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxyy[i] = 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxyz[i] = -2.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xxzz[i] = -4.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xyyy[i] = 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xyyz[i] = -2.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xyzz[i] = -4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_xzzz[i] = -6.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yyyy[i] = 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yyyz[i] = -2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yyzz[i] = -4.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_yzzz[i] = -6.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_x_zzzz[i] = -8.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-255 components of targeted buffer : PG

    auto tr_y_0_z_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 240);

    auto tr_y_0_z_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 241);

    auto tr_y_0_z_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 242);

    auto tr_y_0_z_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 243);

    auto tr_y_0_z_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 244);

    auto tr_y_0_z_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 245);

    auto tr_y_0_z_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 246);

    auto tr_y_0_z_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 247);

    auto tr_y_0_z_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 248);

    auto tr_y_0_z_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 249);

    auto tr_y_0_z_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 250);

    auto tr_y_0_z_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 251);

    auto tr_y_0_z_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 252);

    auto tr_y_0_z_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 253);

    auto tr_y_0_z_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 254);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxz, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_y_0_z_y_xxxx, tr_y_0_z_y_xxxy, tr_y_0_z_y_xxxz, tr_y_0_z_y_xxyy, tr_y_0_z_y_xxyz, tr_y_0_z_y_xxzz, tr_y_0_z_y_xyyy, tr_y_0_z_y_xyyz, tr_y_0_z_y_xyzz, tr_y_0_z_y_xzzz, tr_y_0_z_y_yyyy, tr_y_0_z_y_yyyz, tr_y_0_z_y_yyzz, tr_y_0_z_y_yzzz, tr_y_0_z_y_zzzz, tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yy_zzzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_y_xxxx[i] = -2.0 * tr_0_xxxxz[i] * tke_0 - 2.0 * tr_z_xxxx[i] * tbe_0 + 4.0 * tr_yy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxxy[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_z_xxxy[i] * tbe_0 + 4.0 * tr_yy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxxz[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 2.0 * tr_z_xxxz[i] * tbe_0 - 2.0 * tr_yy_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxyy[i] = -2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_z_xxyy[i] * tbe_0 + 4.0 * tr_yy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxyz[i] = tr_0_xxy[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_z_xxyz[i] * tbe_0 - 2.0 * tr_yy_xxy[i] * tbe_0 + 4.0 * tr_yy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xxzz[i] = 2.0 * tr_0_xxz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_z_xxzz[i] * tbe_0 - 4.0 * tr_yy_xxz[i] * tbe_0 + 4.0 * tr_yy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xyyy[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_z_xyyy[i] * tbe_0 + 4.0 * tr_yy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xyyz[i] = tr_0_xyy[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_z_xyyz[i] * tbe_0 - 2.0 * tr_yy_xyy[i] * tbe_0 + 4.0 * tr_yy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xyzz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_z_xyzz[i] * tbe_0 - 4.0 * tr_yy_xyz[i] * tbe_0 + 4.0 * tr_yy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_xzzz[i] = 3.0 * tr_0_xzz[i] - 2.0 * tr_0_xzzzz[i] * tke_0 - 2.0 * tr_z_xzzz[i] * tbe_0 - 6.0 * tr_yy_xzz[i] * tbe_0 + 4.0 * tr_yy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yyyy[i] = -2.0 * tr_0_yyyyz[i] * tke_0 - 2.0 * tr_z_yyyy[i] * tbe_0 + 4.0 * tr_yy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yyyz[i] = tr_0_yyy[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_z_yyyz[i] * tbe_0 - 2.0 * tr_yy_yyy[i] * tbe_0 + 4.0 * tr_yy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yyzz[i] = 2.0 * tr_0_yyz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_z_yyzz[i] * tbe_0 - 4.0 * tr_yy_yyz[i] * tbe_0 + 4.0 * tr_yy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_yzzz[i] = 3.0 * tr_0_yzz[i] - 2.0 * tr_0_yzzzz[i] * tke_0 - 2.0 * tr_z_yzzz[i] * tbe_0 - 6.0 * tr_yy_yzz[i] * tbe_0 + 4.0 * tr_yy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_y_zzzz[i] = 4.0 * tr_0_zzz[i] - 2.0 * tr_0_zzzzz[i] * tke_0 - 2.0 * tr_z_zzzz[i] * tbe_0 - 8.0 * tr_yy_zzz[i] * tbe_0 + 4.0 * tr_yy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 255-270 components of targeted buffer : PG

    auto tr_y_0_z_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 255);

    auto tr_y_0_z_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 256);

    auto tr_y_0_z_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 257);

    auto tr_y_0_z_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 258);

    auto tr_y_0_z_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 259);

    auto tr_y_0_z_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 260);

    auto tr_y_0_z_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 261);

    auto tr_y_0_z_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 262);

    auto tr_y_0_z_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 263);

    auto tr_y_0_z_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 264);

    auto tr_y_0_z_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 265);

    auto tr_y_0_z_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 266);

    auto tr_y_0_z_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 267);

    auto tr_y_0_z_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 268);

    auto tr_y_0_z_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 269);

    #pragma omp simd aligned(tr_y_0_z_z_xxxx, tr_y_0_z_z_xxxy, tr_y_0_z_z_xxxz, tr_y_0_z_z_xxyy, tr_y_0_z_z_xxyz, tr_y_0_z_z_xxzz, tr_y_0_z_z_xyyy, tr_y_0_z_z_xyyz, tr_y_0_z_z_xyzz, tr_y_0_z_z_xzzz, tr_y_0_z_z_yyyy, tr_y_0_z_z_yyyz, tr_y_0_z_z_yyzz, tr_y_0_z_z_yzzz, tr_y_0_z_z_zzzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_z_xxxx[i] = -2.0 * tr_y_xxxx[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxxy[i] = -2.0 * tr_y_xxxy[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxxz[i] = -2.0 * tr_y_xxxz[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxyy[i] = -2.0 * tr_y_xxyy[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxyz[i] = -2.0 * tr_y_xxyz[i] * tbe_0 - 2.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xxzz[i] = -2.0 * tr_y_xxzz[i] * tbe_0 - 4.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xyyy[i] = -2.0 * tr_y_xyyy[i] * tbe_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xyyz[i] = -2.0 * tr_y_xyyz[i] * tbe_0 - 2.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xyzz[i] = -2.0 * tr_y_xyzz[i] * tbe_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_xzzz[i] = -2.0 * tr_y_xzzz[i] * tbe_0 - 6.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yyyy[i] = -2.0 * tr_y_yyyy[i] * tbe_0 + 4.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yyyz[i] = -2.0 * tr_y_yyyz[i] * tbe_0 - 2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yyzz[i] = -2.0 * tr_y_yyzz[i] * tbe_0 - 4.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_yzzz[i] = -2.0 * tr_y_yzzz[i] * tbe_0 - 6.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_z_zzzz[i] = -2.0 * tr_y_zzzz[i] * tbe_0 - 8.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-285 components of targeted buffer : PG

    auto tr_z_0_x_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 270);

    auto tr_z_0_x_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 271);

    auto tr_z_0_x_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 272);

    auto tr_z_0_x_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 273);

    auto tr_z_0_x_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 274);

    auto tr_z_0_x_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 275);

    auto tr_z_0_x_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 276);

    auto tr_z_0_x_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 277);

    auto tr_z_0_x_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 278);

    auto tr_z_0_x_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 279);

    auto tr_z_0_x_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 280);

    auto tr_z_0_x_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 281);

    auto tr_z_0_x_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 282);

    auto tr_z_0_x_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 283);

    auto tr_z_0_x_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 284);

    #pragma omp simd aligned(tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_0_x_x_xxxx, tr_z_0_x_x_xxxy, tr_z_0_x_x_xxxz, tr_z_0_x_x_xxyy, tr_z_0_x_x_xxyz, tr_z_0_x_x_xxzz, tr_z_0_x_x_xyyy, tr_z_0_x_x_xyyz, tr_z_0_x_x_xyzz, tr_z_0_x_x_xzzz, tr_z_0_x_x_yyyy, tr_z_0_x_x_yyyz, tr_z_0_x_x_yyzz, tr_z_0_x_x_yzzz, tr_z_0_x_x_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_x_xxxx[i] = -2.0 * tr_z_xxxx[i] * tbe_0 - 8.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxxy[i] = -2.0 * tr_z_xxxy[i] * tbe_0 - 6.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxxz[i] = -2.0 * tr_z_xxxz[i] * tbe_0 - 6.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxyy[i] = -2.0 * tr_z_xxyy[i] * tbe_0 - 4.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxyz[i] = -2.0 * tr_z_xxyz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xxzz[i] = -2.0 * tr_z_xxzz[i] * tbe_0 - 4.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xyyy[i] = -2.0 * tr_z_xyyy[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xyyz[i] = -2.0 * tr_z_xyyz[i] * tbe_0 - 2.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xyzz[i] = -2.0 * tr_z_xyzz[i] * tbe_0 - 2.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_xzzz[i] = -2.0 * tr_z_xzzz[i] * tbe_0 - 2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yyyy[i] = -2.0 * tr_z_yyyy[i] * tbe_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yyyz[i] = -2.0 * tr_z_yyyz[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yyzz[i] = -2.0 * tr_z_yyzz[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_yzzz[i] = -2.0 * tr_z_yzzz[i] * tbe_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_x_zzzz[i] = -2.0 * tr_z_zzzz[i] * tbe_0 + 4.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 285-300 components of targeted buffer : PG

    auto tr_z_0_x_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 285);

    auto tr_z_0_x_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 286);

    auto tr_z_0_x_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 287);

    auto tr_z_0_x_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 288);

    auto tr_z_0_x_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 289);

    auto tr_z_0_x_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 290);

    auto tr_z_0_x_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 291);

    auto tr_z_0_x_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 292);

    auto tr_z_0_x_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 293);

    auto tr_z_0_x_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 294);

    auto tr_z_0_x_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 295);

    auto tr_z_0_x_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 296);

    auto tr_z_0_x_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 297);

    auto tr_z_0_x_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 298);

    auto tr_z_0_x_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 299);

    #pragma omp simd aligned(tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_z_0_x_y_xxxx, tr_z_0_x_y_xxxy, tr_z_0_x_y_xxxz, tr_z_0_x_y_xxyy, tr_z_0_x_y_xxyz, tr_z_0_x_y_xxzz, tr_z_0_x_y_xyyy, tr_z_0_x_y_xyyz, tr_z_0_x_y_xyzz, tr_z_0_x_y_xzzz, tr_z_0_x_y_yyyy, tr_z_0_x_y_yyyz, tr_z_0_x_y_yyzz, tr_z_0_x_y_yzzz, tr_z_0_x_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_y_xxxx[i] = -8.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxxy[i] = -6.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxxz[i] = -6.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxyy[i] = -4.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxyz[i] = -4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xxzz[i] = -4.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xyyy[i] = -2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xyyz[i] = -2.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xyzz[i] = -2.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_xzzz[i] = -2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yyyy[i] = 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yyyz[i] = 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yyzz[i] = 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_yzzz[i] = 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_y_zzzz[i] = 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-315 components of targeted buffer : PG

    auto tr_z_0_x_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 300);

    auto tr_z_0_x_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 301);

    auto tr_z_0_x_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 302);

    auto tr_z_0_x_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 303);

    auto tr_z_0_x_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 304);

    auto tr_z_0_x_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 305);

    auto tr_z_0_x_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 306);

    auto tr_z_0_x_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 307);

    auto tr_z_0_x_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 308);

    auto tr_z_0_x_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 309);

    auto tr_z_0_x_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 310);

    auto tr_z_0_x_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 311);

    auto tr_z_0_x_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 312);

    auto tr_z_0_x_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 313);

    auto tr_z_0_x_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 314);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_z_0_x_z_xxxx, tr_z_0_x_z_xxxy, tr_z_0_x_z_xxxz, tr_z_0_x_z_xxyy, tr_z_0_x_z_xxyz, tr_z_0_x_z_xxzz, tr_z_0_x_z_xyyy, tr_z_0_x_z_xyyz, tr_z_0_x_z_xyzz, tr_z_0_x_z_xzzz, tr_z_0_x_z_yyyy, tr_z_0_x_z_yyyz, tr_z_0_x_z_yyzz, tr_z_0_x_z_yzzz, tr_z_0_x_z_zzzz, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_z_xxxx[i] = 4.0 * tr_0_xxx[i] - 2.0 * tr_0_xxxxx[i] * tke_0 - 8.0 * tr_zz_xxx[i] * tbe_0 + 4.0 * tr_zz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxxy[i] = 3.0 * tr_0_xxy[i] - 2.0 * tr_0_xxxxy[i] * tke_0 - 6.0 * tr_zz_xxy[i] * tbe_0 + 4.0 * tr_zz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxy[i] * tbe_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxxz[i] = 3.0 * tr_0_xxz[i] - 2.0 * tr_0_xxxxz[i] * tke_0 - 6.0 * tr_zz_xxz[i] * tbe_0 + 4.0 * tr_zz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxz[i] * tbe_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxyy[i] = 2.0 * tr_0_xyy[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 4.0 * tr_zz_xyy[i] * tbe_0 + 4.0 * tr_zz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyy[i] * tbe_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xxxyz[i] * tke_0 - 4.0 * tr_zz_xyz[i] * tbe_0 + 4.0 * tr_zz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyz[i] * tbe_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xxzz[i] = 2.0 * tr_0_xzz[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 4.0 * tr_zz_xzz[i] * tbe_0 + 4.0 * tr_zz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxzz[i] * tbe_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xyyy[i] = tr_0_yyy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 2.0 * tr_zz_yyy[i] * tbe_0 + 4.0 * tr_zz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyy[i] * tbe_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xyyz[i] = tr_0_yyz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_zz_yyz[i] * tbe_0 + 4.0 * tr_zz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyz[i] * tbe_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xyzz[i] = tr_0_yzz[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_zz_yzz[i] * tbe_0 + 4.0 * tr_zz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyzz[i] * tbe_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_xzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_zz_zzz[i] * tbe_0 + 4.0 * tr_zz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzzz[i] * tbe_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yyyy[i] = -2.0 * tr_0_xyyyy[i] * tke_0 + 4.0 * tr_zz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyy[i] * tbe_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yyyz[i] = -2.0 * tr_0_xyyyz[i] * tke_0 + 4.0 * tr_zz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyz[i] * tbe_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yyzz[i] = -2.0 * tr_0_xyyzz[i] * tke_0 + 4.0 * tr_zz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyzz[i] * tbe_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_yzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 + 4.0 * tr_zz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzzz[i] * tbe_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_z_zzzz[i] = -2.0 * tr_0_xzzzz[i] * tke_0 + 4.0 * tr_zz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzzz[i] * tbe_0 + 4.0 * tr_xzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 315-330 components of targeted buffer : PG

    auto tr_z_0_y_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 315);

    auto tr_z_0_y_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 316);

    auto tr_z_0_y_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 317);

    auto tr_z_0_y_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 318);

    auto tr_z_0_y_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 319);

    auto tr_z_0_y_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 320);

    auto tr_z_0_y_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 321);

    auto tr_z_0_y_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 322);

    auto tr_z_0_y_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 323);

    auto tr_z_0_y_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 324);

    auto tr_z_0_y_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 325);

    auto tr_z_0_y_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 326);

    auto tr_z_0_y_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 327);

    auto tr_z_0_y_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 328);

    auto tr_z_0_y_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 329);

    #pragma omp simd aligned(tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_z_0_y_x_xxxx, tr_z_0_y_x_xxxy, tr_z_0_y_x_xxxz, tr_z_0_y_x_xxyy, tr_z_0_y_x_xxyz, tr_z_0_y_x_xxzz, tr_z_0_y_x_xyyy, tr_z_0_y_x_xyyz, tr_z_0_y_x_xyzz, tr_z_0_y_x_xzzz, tr_z_0_y_x_yyyy, tr_z_0_y_x_yyyz, tr_z_0_y_x_yyzz, tr_z_0_y_x_yzzz, tr_z_0_y_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_x_xxxx[i] = 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxxy[i] = -2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxxz[i] = 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxyy[i] = -4.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxyz[i] = -2.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xxzz[i] = 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xyyy[i] = -6.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xyyz[i] = -4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xyzz[i] = -2.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_xzzz[i] = 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yyyy[i] = -8.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yyyz[i] = -6.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yyzz[i] = -4.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_yzzz[i] = -2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_x_zzzz[i] = 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-345 components of targeted buffer : PG

    auto tr_z_0_y_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 330);

    auto tr_z_0_y_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 331);

    auto tr_z_0_y_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 332);

    auto tr_z_0_y_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 333);

    auto tr_z_0_y_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 334);

    auto tr_z_0_y_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 335);

    auto tr_z_0_y_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 336);

    auto tr_z_0_y_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 337);

    auto tr_z_0_y_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 338);

    auto tr_z_0_y_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 339);

    auto tr_z_0_y_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 340);

    auto tr_z_0_y_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 341);

    auto tr_z_0_y_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 342);

    auto tr_z_0_y_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 343);

    auto tr_z_0_y_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 344);

    #pragma omp simd aligned(tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_0_y_y_xxxx, tr_z_0_y_y_xxxy, tr_z_0_y_y_xxxz, tr_z_0_y_y_xxyy, tr_z_0_y_y_xxyz, tr_z_0_y_y_xxzz, tr_z_0_y_y_xyyy, tr_z_0_y_y_xyyz, tr_z_0_y_y_xyzz, tr_z_0_y_y_xzzz, tr_z_0_y_y_yyyy, tr_z_0_y_y_yyyz, tr_z_0_y_y_yyzz, tr_z_0_y_y_yzzz, tr_z_0_y_y_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_y_xxxx[i] = -2.0 * tr_z_xxxx[i] * tbe_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxxy[i] = -2.0 * tr_z_xxxy[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxxz[i] = -2.0 * tr_z_xxxz[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxyy[i] = -2.0 * tr_z_xxyy[i] * tbe_0 - 4.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxyz[i] = -2.0 * tr_z_xxyz[i] * tbe_0 - 2.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xxzz[i] = -2.0 * tr_z_xxzz[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xyyy[i] = -2.0 * tr_z_xyyy[i] * tbe_0 - 6.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xyyz[i] = -2.0 * tr_z_xyyz[i] * tbe_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xyzz[i] = -2.0 * tr_z_xyzz[i] * tbe_0 - 2.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_xzzz[i] = -2.0 * tr_z_xzzz[i] * tbe_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yyyy[i] = -2.0 * tr_z_yyyy[i] * tbe_0 - 8.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yyyz[i] = -2.0 * tr_z_yyyz[i] * tbe_0 - 6.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yyzz[i] = -2.0 * tr_z_yyzz[i] * tbe_0 - 4.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_yzzz[i] = -2.0 * tr_z_yzzz[i] * tbe_0 - 2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_y_zzzz[i] = -2.0 * tr_z_zzzz[i] * tbe_0 + 4.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 345-360 components of targeted buffer : PG

    auto tr_z_0_y_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 345);

    auto tr_z_0_y_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 346);

    auto tr_z_0_y_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 347);

    auto tr_z_0_y_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 348);

    auto tr_z_0_y_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 349);

    auto tr_z_0_y_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 350);

    auto tr_z_0_y_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 351);

    auto tr_z_0_y_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 352);

    auto tr_z_0_y_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 353);

    auto tr_z_0_y_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 354);

    auto tr_z_0_y_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 355);

    auto tr_z_0_y_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 356);

    auto tr_z_0_y_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 357);

    auto tr_z_0_y_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 358);

    auto tr_z_0_y_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 359);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxy, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_z_0_y_z_xxxx, tr_z_0_y_z_xxxy, tr_z_0_y_z_xxxz, tr_z_0_y_z_xxyy, tr_z_0_y_z_xxyz, tr_z_0_y_z_xxzz, tr_z_0_y_z_xyyy, tr_z_0_y_z_xyyz, tr_z_0_y_z_xyzz, tr_z_0_y_z_xzzz, tr_z_0_y_z_yyyy, tr_z_0_y_z_yyyz, tr_z_0_y_z_yyzz, tr_z_0_y_z_yzzz, tr_z_0_y_z_zzzz, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_z_xxxx[i] = -2.0 * tr_0_xxxxy[i] * tke_0 + 4.0 * tr_zz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxxx[i] * tbe_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxxy[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 2.0 * tr_zz_xxx[i] * tbe_0 + 4.0 * tr_zz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxxy[i] * tbe_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxxz[i] = -2.0 * tr_0_xxxyz[i] * tke_0 + 4.0 * tr_zz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxxz[i] * tbe_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxyy[i] = 2.0 * tr_0_xxy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 4.0 * tr_zz_xxy[i] * tbe_0 + 4.0 * tr_zz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxyy[i] * tbe_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxyz[i] = tr_0_xxz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_zz_xxz[i] * tbe_0 + 4.0 * tr_zz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxyz[i] * tbe_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xxzz[i] = -2.0 * tr_0_xxyzz[i] * tke_0 + 4.0 * tr_zz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxzz[i] * tbe_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xyyy[i] = 3.0 * tr_0_xyy[i] - 2.0 * tr_0_xyyyy[i] * tke_0 - 6.0 * tr_zz_xyy[i] * tbe_0 + 4.0 * tr_zz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyyy[i] * tbe_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xyyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyyyz[i] * tke_0 - 4.0 * tr_zz_xyz[i] * tbe_0 + 4.0 * tr_zz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyyz[i] * tbe_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xyzz[i] = tr_0_xzz[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_zz_xzz[i] * tbe_0 + 4.0 * tr_zz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyzz[i] * tbe_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_xzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 + 4.0 * tr_zz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xzzz[i] * tbe_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yyyy[i] = 4.0 * tr_0_yyy[i] - 2.0 * tr_0_yyyyy[i] * tke_0 - 8.0 * tr_zz_yyy[i] * tbe_0 + 4.0 * tr_zz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyyy[i] * tbe_0 + 4.0 * tr_yzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yyyz[i] = 3.0 * tr_0_yyz[i] - 2.0 * tr_0_yyyyz[i] * tke_0 - 6.0 * tr_zz_yyz[i] * tbe_0 + 4.0 * tr_zz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyyz[i] * tbe_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yyzz[i] = 2.0 * tr_0_yzz[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 4.0 * tr_zz_yzz[i] * tbe_0 + 4.0 * tr_zz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyzz[i] * tbe_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_yzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_zz_zzz[i] * tbe_0 + 4.0 * tr_zz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yzzz[i] * tbe_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_z_zzzz[i] = -2.0 * tr_0_yzzzz[i] * tke_0 + 4.0 * tr_zz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_zzzz[i] * tbe_0 + 4.0 * tr_yzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-375 components of targeted buffer : PG

    auto tr_z_0_z_x_xxxx = pbuffer.data(idx_op_geom_110_pg + 360);

    auto tr_z_0_z_x_xxxy = pbuffer.data(idx_op_geom_110_pg + 361);

    auto tr_z_0_z_x_xxxz = pbuffer.data(idx_op_geom_110_pg + 362);

    auto tr_z_0_z_x_xxyy = pbuffer.data(idx_op_geom_110_pg + 363);

    auto tr_z_0_z_x_xxyz = pbuffer.data(idx_op_geom_110_pg + 364);

    auto tr_z_0_z_x_xxzz = pbuffer.data(idx_op_geom_110_pg + 365);

    auto tr_z_0_z_x_xyyy = pbuffer.data(idx_op_geom_110_pg + 366);

    auto tr_z_0_z_x_xyyz = pbuffer.data(idx_op_geom_110_pg + 367);

    auto tr_z_0_z_x_xyzz = pbuffer.data(idx_op_geom_110_pg + 368);

    auto tr_z_0_z_x_xzzz = pbuffer.data(idx_op_geom_110_pg + 369);

    auto tr_z_0_z_x_yyyy = pbuffer.data(idx_op_geom_110_pg + 370);

    auto tr_z_0_z_x_yyyz = pbuffer.data(idx_op_geom_110_pg + 371);

    auto tr_z_0_z_x_yyzz = pbuffer.data(idx_op_geom_110_pg + 372);

    auto tr_z_0_z_x_yzzz = pbuffer.data(idx_op_geom_110_pg + 373);

    auto tr_z_0_z_x_zzzz = pbuffer.data(idx_op_geom_110_pg + 374);

    #pragma omp simd aligned(tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_z_0_z_x_xxxx, tr_z_0_z_x_xxxy, tr_z_0_z_x_xxxz, tr_z_0_z_x_xxyy, tr_z_0_z_x_xxyz, tr_z_0_z_x_xxzz, tr_z_0_z_x_xyyy, tr_z_0_z_x_xyyz, tr_z_0_z_x_xyzz, tr_z_0_z_x_xzzz, tr_z_0_z_x_yyyy, tr_z_0_z_x_yyyz, tr_z_0_z_x_yyzz, tr_z_0_z_x_yzzz, tr_z_0_z_x_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_x_xxxx[i] = -2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxxy[i] = -2.0 * tr_x_xxxy[i] * tbe_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxxz[i] = -2.0 * tr_x_xxxz[i] * tbe_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxyy[i] = -2.0 * tr_x_xxyy[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxyz[i] = -2.0 * tr_x_xxyz[i] * tbe_0 - 2.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xxzz[i] = -2.0 * tr_x_xxzz[i] * tbe_0 - 4.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xyyy[i] = -2.0 * tr_x_xyyy[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xyyz[i] = -2.0 * tr_x_xyyz[i] * tbe_0 - 2.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xyzz[i] = -2.0 * tr_x_xyzz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_xzzz[i] = -2.0 * tr_x_xzzz[i] * tbe_0 - 6.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yyyy[i] = -2.0 * tr_x_yyyy[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yyyz[i] = -2.0 * tr_x_yyyz[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yyzz[i] = -2.0 * tr_x_yyzz[i] * tbe_0 - 4.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_yzzz[i] = -2.0 * tr_x_yzzz[i] * tbe_0 - 6.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_x_zzzz[i] = -2.0 * tr_x_zzzz[i] * tbe_0 - 8.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 375-390 components of targeted buffer : PG

    auto tr_z_0_z_y_xxxx = pbuffer.data(idx_op_geom_110_pg + 375);

    auto tr_z_0_z_y_xxxy = pbuffer.data(idx_op_geom_110_pg + 376);

    auto tr_z_0_z_y_xxxz = pbuffer.data(idx_op_geom_110_pg + 377);

    auto tr_z_0_z_y_xxyy = pbuffer.data(idx_op_geom_110_pg + 378);

    auto tr_z_0_z_y_xxyz = pbuffer.data(idx_op_geom_110_pg + 379);

    auto tr_z_0_z_y_xxzz = pbuffer.data(idx_op_geom_110_pg + 380);

    auto tr_z_0_z_y_xyyy = pbuffer.data(idx_op_geom_110_pg + 381);

    auto tr_z_0_z_y_xyyz = pbuffer.data(idx_op_geom_110_pg + 382);

    auto tr_z_0_z_y_xyzz = pbuffer.data(idx_op_geom_110_pg + 383);

    auto tr_z_0_z_y_xzzz = pbuffer.data(idx_op_geom_110_pg + 384);

    auto tr_z_0_z_y_yyyy = pbuffer.data(idx_op_geom_110_pg + 385);

    auto tr_z_0_z_y_yyyz = pbuffer.data(idx_op_geom_110_pg + 386);

    auto tr_z_0_z_y_yyzz = pbuffer.data(idx_op_geom_110_pg + 387);

    auto tr_z_0_z_y_yzzz = pbuffer.data(idx_op_geom_110_pg + 388);

    auto tr_z_0_z_y_zzzz = pbuffer.data(idx_op_geom_110_pg + 389);

    #pragma omp simd aligned(tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_z_0_z_y_xxxx, tr_z_0_z_y_xxxy, tr_z_0_z_y_xxxz, tr_z_0_z_y_xxyy, tr_z_0_z_y_xxyz, tr_z_0_z_y_xxzz, tr_z_0_z_y_xyyy, tr_z_0_z_y_xyyz, tr_z_0_z_y_xyzz, tr_z_0_z_y_xzzz, tr_z_0_z_y_yyyy, tr_z_0_z_y_yyyz, tr_z_0_z_y_yyzz, tr_z_0_z_y_yzzz, tr_z_0_z_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_y_xxxx[i] = -2.0 * tr_y_xxxx[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxxy[i] = -2.0 * tr_y_xxxy[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxxz[i] = -2.0 * tr_y_xxxz[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxyy[i] = -2.0 * tr_y_xxyy[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxyz[i] = -2.0 * tr_y_xxyz[i] * tbe_0 - 2.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xxzz[i] = -2.0 * tr_y_xxzz[i] * tbe_0 - 4.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xyyy[i] = -2.0 * tr_y_xyyy[i] * tbe_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xyyz[i] = -2.0 * tr_y_xyyz[i] * tbe_0 - 2.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xyzz[i] = -2.0 * tr_y_xyzz[i] * tbe_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_xzzz[i] = -2.0 * tr_y_xzzz[i] * tbe_0 - 6.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yyyy[i] = -2.0 * tr_y_yyyy[i] * tbe_0 + 4.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yyyz[i] = -2.0 * tr_y_yyyz[i] * tbe_0 - 2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yyzz[i] = -2.0 * tr_y_yyzz[i] * tbe_0 - 4.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_yzzz[i] = -2.0 * tr_y_yzzz[i] * tbe_0 - 6.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_y_zzzz[i] = -2.0 * tr_y_zzzz[i] * tbe_0 - 8.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-405 components of targeted buffer : PG

    auto tr_z_0_z_z_xxxx = pbuffer.data(idx_op_geom_110_pg + 390);

    auto tr_z_0_z_z_xxxy = pbuffer.data(idx_op_geom_110_pg + 391);

    auto tr_z_0_z_z_xxxz = pbuffer.data(idx_op_geom_110_pg + 392);

    auto tr_z_0_z_z_xxyy = pbuffer.data(idx_op_geom_110_pg + 393);

    auto tr_z_0_z_z_xxyz = pbuffer.data(idx_op_geom_110_pg + 394);

    auto tr_z_0_z_z_xxzz = pbuffer.data(idx_op_geom_110_pg + 395);

    auto tr_z_0_z_z_xyyy = pbuffer.data(idx_op_geom_110_pg + 396);

    auto tr_z_0_z_z_xyyz = pbuffer.data(idx_op_geom_110_pg + 397);

    auto tr_z_0_z_z_xyzz = pbuffer.data(idx_op_geom_110_pg + 398);

    auto tr_z_0_z_z_xzzz = pbuffer.data(idx_op_geom_110_pg + 399);

    auto tr_z_0_z_z_yyyy = pbuffer.data(idx_op_geom_110_pg + 400);

    auto tr_z_0_z_z_yyyz = pbuffer.data(idx_op_geom_110_pg + 401);

    auto tr_z_0_z_z_yyzz = pbuffer.data(idx_op_geom_110_pg + 402);

    auto tr_z_0_z_z_yzzz = pbuffer.data(idx_op_geom_110_pg + 403);

    auto tr_z_0_z_z_zzzz = pbuffer.data(idx_op_geom_110_pg + 404);

    #pragma omp simd aligned(tr_0_xxx, tr_0_xxxxz, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_z_0_z_z_xxxx, tr_z_0_z_z_xxxy, tr_z_0_z_z_xxxz, tr_z_0_z_z_xxyy, tr_z_0_z_z_xxyz, tr_z_0_z_z_xxzz, tr_z_0_z_z_xyyy, tr_z_0_z_z_xyyz, tr_z_0_z_z_xyzz, tr_z_0_z_z_xzzz, tr_z_0_z_z_yyyy, tr_z_0_z_z_yyyz, tr_z_0_z_z_yyzz, tr_z_0_z_z_yzzz, tr_z_0_z_z_zzzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_z_xxxx[i] = -2.0 * tr_0_xxxxz[i] * tke_0 - 6.0 * tr_z_xxxx[i] * tbe_0 + 4.0 * tr_zz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxxy[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 6.0 * tr_z_xxxy[i] * tbe_0 + 4.0 * tr_zz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxxz[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 6.0 * tr_z_xxxz[i] * tbe_0 - 2.0 * tr_zz_xxx[i] * tbe_0 + 4.0 * tr_zz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxyy[i] = -2.0 * tr_0_xxyyz[i] * tke_0 - 6.0 * tr_z_xxyy[i] * tbe_0 + 4.0 * tr_zz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxyz[i] = tr_0_xxy[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 6.0 * tr_z_xxyz[i] * tbe_0 - 2.0 * tr_zz_xxy[i] * tbe_0 + 4.0 * tr_zz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xxzz[i] = 2.0 * tr_0_xxz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 6.0 * tr_z_xxzz[i] * tbe_0 - 4.0 * tr_zz_xxz[i] * tbe_0 + 4.0 * tr_zz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xyyy[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 6.0 * tr_z_xyyy[i] * tbe_0 + 4.0 * tr_zz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xyyz[i] = tr_0_xyy[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 6.0 * tr_z_xyyz[i] * tbe_0 - 2.0 * tr_zz_xyy[i] * tbe_0 + 4.0 * tr_zz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xyzz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_z_xyzz[i] * tbe_0 - 4.0 * tr_zz_xyz[i] * tbe_0 + 4.0 * tr_zz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_xzzz[i] = 3.0 * tr_0_xzz[i] - 2.0 * tr_0_xzzzz[i] * tke_0 - 6.0 * tr_z_xzzz[i] * tbe_0 - 6.0 * tr_zz_xzz[i] * tbe_0 + 4.0 * tr_zz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yyyy[i] = -2.0 * tr_0_yyyyz[i] * tke_0 - 6.0 * tr_z_yyyy[i] * tbe_0 + 4.0 * tr_zz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yyyz[i] = tr_0_yyy[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 6.0 * tr_z_yyyz[i] * tbe_0 - 2.0 * tr_zz_yyy[i] * tbe_0 + 4.0 * tr_zz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yyzz[i] = 2.0 * tr_0_yyz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 6.0 * tr_z_yyzz[i] * tbe_0 - 4.0 * tr_zz_yyz[i] * tbe_0 + 4.0 * tr_zz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_yzzz[i] = 3.0 * tr_0_yzz[i] - 2.0 * tr_0_yzzzz[i] * tke_0 - 6.0 * tr_z_yzzz[i] * tbe_0 - 6.0 * tr_zz_yzz[i] * tbe_0 + 4.0 * tr_zz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yzzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_z_zzzz[i] = 4.0 * tr_0_zzz[i] - 2.0 * tr_0_zzzzz[i] * tke_0 - 6.0 * tr_z_zzzz[i] * tbe_0 - 8.0 * tr_zz_zzz[i] * tbe_0 + 4.0 * tr_zz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_zzzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

