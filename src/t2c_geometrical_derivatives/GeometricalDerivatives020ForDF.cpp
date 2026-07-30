#include "GeometricalDerivatives020ForDF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_df(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_df,
                         const int idx_op_sf,
                         const int idx_op_pd,
                         const int idx_op_pg,
                         const int idx_op_dp,
                         const int idx_op_df,
                         const int idx_op_dh,
                         const int idx_op_fd,
                         const int idx_op_fg,
                         const int idx_op_gf,
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

    // Set up components of auxiliary buffer : PD

    auto tr_x_xx = pbuffer.data(idx_op_pd);

    auto tr_x_xy = pbuffer.data(idx_op_pd + 1);

    auto tr_x_xz = pbuffer.data(idx_op_pd + 2);

    auto tr_x_yy = pbuffer.data(idx_op_pd + 3);

    auto tr_x_yz = pbuffer.data(idx_op_pd + 4);

    auto tr_x_zz = pbuffer.data(idx_op_pd + 5);

    auto tr_y_xx = pbuffer.data(idx_op_pd + 6);

    auto tr_y_xy = pbuffer.data(idx_op_pd + 7);

    auto tr_y_xz = pbuffer.data(idx_op_pd + 8);

    auto tr_y_yy = pbuffer.data(idx_op_pd + 9);

    auto tr_y_yz = pbuffer.data(idx_op_pd + 10);

    auto tr_y_zz = pbuffer.data(idx_op_pd + 11);

    auto tr_z_xx = pbuffer.data(idx_op_pd + 12);

    auto tr_z_xy = pbuffer.data(idx_op_pd + 13);

    auto tr_z_xz = pbuffer.data(idx_op_pd + 14);

    auto tr_z_yy = pbuffer.data(idx_op_pd + 15);

    auto tr_z_yz = pbuffer.data(idx_op_pd + 16);

    auto tr_z_zz = pbuffer.data(idx_op_pd + 17);

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

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

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

    // Set up components of auxiliary buffer : FD

    auto tr_xxx_xx = pbuffer.data(idx_op_fd);

    auto tr_xxx_xy = pbuffer.data(idx_op_fd + 1);

    auto tr_xxx_xz = pbuffer.data(idx_op_fd + 2);

    auto tr_xxx_yy = pbuffer.data(idx_op_fd + 3);

    auto tr_xxx_yz = pbuffer.data(idx_op_fd + 4);

    auto tr_xxx_zz = pbuffer.data(idx_op_fd + 5);

    auto tr_xxy_xx = pbuffer.data(idx_op_fd + 6);

    auto tr_xxy_xy = pbuffer.data(idx_op_fd + 7);

    auto tr_xxy_xz = pbuffer.data(idx_op_fd + 8);

    auto tr_xxy_yy = pbuffer.data(idx_op_fd + 9);

    auto tr_xxy_yz = pbuffer.data(idx_op_fd + 10);

    auto tr_xxy_zz = pbuffer.data(idx_op_fd + 11);

    auto tr_xxz_xx = pbuffer.data(idx_op_fd + 12);

    auto tr_xxz_xy = pbuffer.data(idx_op_fd + 13);

    auto tr_xxz_xz = pbuffer.data(idx_op_fd + 14);

    auto tr_xxz_yy = pbuffer.data(idx_op_fd + 15);

    auto tr_xxz_yz = pbuffer.data(idx_op_fd + 16);

    auto tr_xxz_zz = pbuffer.data(idx_op_fd + 17);

    auto tr_xyy_xx = pbuffer.data(idx_op_fd + 18);

    auto tr_xyy_xy = pbuffer.data(idx_op_fd + 19);

    auto tr_xyy_xz = pbuffer.data(idx_op_fd + 20);

    auto tr_xyy_yy = pbuffer.data(idx_op_fd + 21);

    auto tr_xyy_yz = pbuffer.data(idx_op_fd + 22);

    auto tr_xyy_zz = pbuffer.data(idx_op_fd + 23);

    auto tr_xyz_xx = pbuffer.data(idx_op_fd + 24);

    auto tr_xyz_xy = pbuffer.data(idx_op_fd + 25);

    auto tr_xyz_xz = pbuffer.data(idx_op_fd + 26);

    auto tr_xyz_yy = pbuffer.data(idx_op_fd + 27);

    auto tr_xyz_yz = pbuffer.data(idx_op_fd + 28);

    auto tr_xyz_zz = pbuffer.data(idx_op_fd + 29);

    auto tr_xzz_xx = pbuffer.data(idx_op_fd + 30);

    auto tr_xzz_xy = pbuffer.data(idx_op_fd + 31);

    auto tr_xzz_xz = pbuffer.data(idx_op_fd + 32);

    auto tr_xzz_yy = pbuffer.data(idx_op_fd + 33);

    auto tr_xzz_yz = pbuffer.data(idx_op_fd + 34);

    auto tr_xzz_zz = pbuffer.data(idx_op_fd + 35);

    auto tr_yyy_xx = pbuffer.data(idx_op_fd + 36);

    auto tr_yyy_xy = pbuffer.data(idx_op_fd + 37);

    auto tr_yyy_xz = pbuffer.data(idx_op_fd + 38);

    auto tr_yyy_yy = pbuffer.data(idx_op_fd + 39);

    auto tr_yyy_yz = pbuffer.data(idx_op_fd + 40);

    auto tr_yyy_zz = pbuffer.data(idx_op_fd + 41);

    auto tr_yyz_xx = pbuffer.data(idx_op_fd + 42);

    auto tr_yyz_xy = pbuffer.data(idx_op_fd + 43);

    auto tr_yyz_xz = pbuffer.data(idx_op_fd + 44);

    auto tr_yyz_yy = pbuffer.data(idx_op_fd + 45);

    auto tr_yyz_yz = pbuffer.data(idx_op_fd + 46);

    auto tr_yyz_zz = pbuffer.data(idx_op_fd + 47);

    auto tr_yzz_xx = pbuffer.data(idx_op_fd + 48);

    auto tr_yzz_xy = pbuffer.data(idx_op_fd + 49);

    auto tr_yzz_xz = pbuffer.data(idx_op_fd + 50);

    auto tr_yzz_yy = pbuffer.data(idx_op_fd + 51);

    auto tr_yzz_yz = pbuffer.data(idx_op_fd + 52);

    auto tr_yzz_zz = pbuffer.data(idx_op_fd + 53);

    auto tr_zzz_xx = pbuffer.data(idx_op_fd + 54);

    auto tr_zzz_xy = pbuffer.data(idx_op_fd + 55);

    auto tr_zzz_xz = pbuffer.data(idx_op_fd + 56);

    auto tr_zzz_yy = pbuffer.data(idx_op_fd + 57);

    auto tr_zzz_yz = pbuffer.data(idx_op_fd + 58);

    auto tr_zzz_zz = pbuffer.data(idx_op_fd + 59);

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

    // Set up 0-10 components of targeted buffer : DF

    auto tr_0_0_xx_xx_xxx = pbuffer.data(idx_op_geom_020_df);

    auto tr_0_0_xx_xx_xxy = pbuffer.data(idx_op_geom_020_df + 1);

    auto tr_0_0_xx_xx_xxz = pbuffer.data(idx_op_geom_020_df + 2);

    auto tr_0_0_xx_xx_xyy = pbuffer.data(idx_op_geom_020_df + 3);

    auto tr_0_0_xx_xx_xyz = pbuffer.data(idx_op_geom_020_df + 4);

    auto tr_0_0_xx_xx_xzz = pbuffer.data(idx_op_geom_020_df + 5);

    auto tr_0_0_xx_xx_yyy = pbuffer.data(idx_op_geom_020_df + 6);

    auto tr_0_0_xx_xx_yyz = pbuffer.data(idx_op_geom_020_df + 7);

    auto tr_0_0_xx_xx_yzz = pbuffer.data(idx_op_geom_020_df + 8);

    auto tr_0_0_xx_xx_zzz = pbuffer.data(idx_op_geom_020_df + 9);

    #pragma omp simd aligned(tr_0_0_xx_xx_xxx, tr_0_0_xx_xx_xxy, tr_0_0_xx_xx_xxz, tr_0_0_xx_xx_xyy, tr_0_0_xx_xx_xyz, tr_0_0_xx_xx_xzz, tr_0_0_xx_xx_yyy, tr_0_0_xx_xx_yyz, tr_0_0_xx_xx_yzz, tr_0_0_xx_xx_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xx, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxx_xx, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xx_xxx[i] = 2.0 * tr_0_xxx[i] + 12.0 * tr_x_xx[i] - 8.0 * tr_x_xxxx[i] * tke_0 + 6.0 * tr_xx_x[i] - 10.0 * tr_xx_xxx[i] * tbe_0 - 14.0 * tr_xx_xxx[i] * tke_0 + 4.0 * tr_xx_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxx_xx[i] * tbe_0 + 8.0 * tr_xxx_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xxy[i] = 2.0 * tr_0_xxy[i] + 8.0 * tr_x_xy[i] - 8.0 * tr_x_xxxy[i] * tke_0 + 2.0 * tr_xx_y[i] - 10.0 * tr_xx_xxy[i] * tbe_0 - 10.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxx_xy[i] * tbe_0 + 8.0 * tr_xxx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xxz[i] = 2.0 * tr_0_xxz[i] + 8.0 * tr_x_xz[i] - 8.0 * tr_x_xxxz[i] * tke_0 + 2.0 * tr_xx_z[i] - 10.0 * tr_xx_xxz[i] * tbe_0 - 10.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxx_xz[i] * tbe_0 + 8.0 * tr_xxx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xyy[i] = 2.0 * tr_0_xyy[i] + 4.0 * tr_x_yy[i] - 8.0 * tr_x_xxyy[i] * tke_0 - 10.0 * tr_xx_xyy[i] * tbe_0 - 6.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xx_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxx_yy[i] * tbe_0 + 8.0 * tr_xxx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xyz[i] = 2.0 * tr_0_xyz[i] + 4.0 * tr_x_yz[i] - 8.0 * tr_x_xxyz[i] * tke_0 - 10.0 * tr_xx_xyz[i] * tbe_0 - 6.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxx_yz[i] * tbe_0 + 8.0 * tr_xxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_xzz[i] = 2.0 * tr_0_xzz[i] + 4.0 * tr_x_zz[i] - 8.0 * tr_x_xxzz[i] * tke_0 - 10.0 * tr_xx_xzz[i] * tbe_0 - 6.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xx_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxx_zz[i] * tbe_0 + 8.0 * tr_xxx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_yyy[i] = 2.0 * tr_0_yyy[i] - 8.0 * tr_x_xyyy[i] * tke_0 - 10.0 * tr_xx_yyy[i] * tbe_0 - 2.0 * tr_xx_yyy[i] * tke_0 + 4.0 * tr_xx_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_yyz[i] = 2.0 * tr_0_yyz[i] - 8.0 * tr_x_xyyz[i] * tke_0 - 10.0 * tr_xx_yyz[i] * tbe_0 - 2.0 * tr_xx_yyz[i] * tke_0 + 4.0 * tr_xx_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_yzz[i] = 2.0 * tr_0_yzz[i] - 8.0 * tr_x_xyzz[i] * tke_0 - 10.0 * tr_xx_yzz[i] * tbe_0 - 2.0 * tr_xx_yzz[i] * tke_0 + 4.0 * tr_xx_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xx_zzz[i] = 2.0 * tr_0_zzz[i] - 8.0 * tr_x_xzzz[i] * tke_0 - 10.0 * tr_xx_zzz[i] * tbe_0 - 2.0 * tr_xx_zzz[i] * tke_0 + 4.0 * tr_xx_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxx_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 10-20 components of targeted buffer : DF

    auto tr_0_0_xx_xy_xxx = pbuffer.data(idx_op_geom_020_df + 10);

    auto tr_0_0_xx_xy_xxy = pbuffer.data(idx_op_geom_020_df + 11);

    auto tr_0_0_xx_xy_xxz = pbuffer.data(idx_op_geom_020_df + 12);

    auto tr_0_0_xx_xy_xyy = pbuffer.data(idx_op_geom_020_df + 13);

    auto tr_0_0_xx_xy_xyz = pbuffer.data(idx_op_geom_020_df + 14);

    auto tr_0_0_xx_xy_xzz = pbuffer.data(idx_op_geom_020_df + 15);

    auto tr_0_0_xx_xy_yyy = pbuffer.data(idx_op_geom_020_df + 16);

    auto tr_0_0_xx_xy_yyz = pbuffer.data(idx_op_geom_020_df + 17);

    auto tr_0_0_xx_xy_yzz = pbuffer.data(idx_op_geom_020_df + 18);

    auto tr_0_0_xx_xy_zzz = pbuffer.data(idx_op_geom_020_df + 19);

    #pragma omp simd aligned(tr_0_0_xx_xy_xxx, tr_0_0_xx_xy_xxy, tr_0_0_xx_xy_xxz, tr_0_0_xx_xy_xyy, tr_0_0_xx_xy_xyz, tr_0_0_xx_xy_xzz, tr_0_0_xx_xy_yyy, tr_0_0_xx_xy_yyz, tr_0_0_xx_xy_yzz, tr_0_0_xx_xy_zzz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xx, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xy_x, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_y_xx, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yz, tr_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xy_xxx[i] = 6.0 * tr_y_xx[i] - 4.0 * tr_y_xxxx[i] * tke_0 + 6.0 * tr_xy_x[i] - 6.0 * tr_xy_xxx[i] * tbe_0 - 14.0 * tr_xy_xxx[i] * tke_0 + 4.0 * tr_xy_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxy_xx[i] * tbe_0 + 8.0 * tr_xxy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xxy[i] = 4.0 * tr_y_xy[i] - 4.0 * tr_y_xxxy[i] * tke_0 + 2.0 * tr_xy_y[i] - 6.0 * tr_xy_xxy[i] * tbe_0 - 10.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxy_xy[i] * tbe_0 + 8.0 * tr_xxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xxz[i] = 4.0 * tr_y_xz[i] - 4.0 * tr_y_xxxz[i] * tke_0 + 2.0 * tr_xy_z[i] - 6.0 * tr_xy_xxz[i] * tbe_0 - 10.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxy_xz[i] * tbe_0 + 8.0 * tr_xxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xyy[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_y_xxyy[i] * tke_0 - 6.0 * tr_xy_xyy[i] * tbe_0 - 6.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxy_yy[i] * tbe_0 + 8.0 * tr_xxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xyz[i] = 2.0 * tr_y_yz[i] - 4.0 * tr_y_xxyz[i] * tke_0 - 6.0 * tr_xy_xyz[i] * tbe_0 - 6.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_yz[i] * tbe_0 + 8.0 * tr_xxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_xzz[i] = 2.0 * tr_y_zz[i] - 4.0 * tr_y_xxzz[i] * tke_0 - 6.0 * tr_xy_xzz[i] * tbe_0 - 6.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_zz[i] * tbe_0 + 8.0 * tr_xxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_yyy[i] = -4.0 * tr_y_xyyy[i] * tke_0 - 6.0 * tr_xy_yyy[i] * tbe_0 - 2.0 * tr_xy_yyy[i] * tke_0 + 4.0 * tr_xy_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_yyz[i] = -4.0 * tr_y_xyyz[i] * tke_0 - 6.0 * tr_xy_yyz[i] * tbe_0 - 2.0 * tr_xy_yyz[i] * tke_0 + 4.0 * tr_xy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_yzz[i] = -4.0 * tr_y_xyzz[i] * tke_0 - 6.0 * tr_xy_yzz[i] * tbe_0 - 2.0 * tr_xy_yzz[i] * tke_0 + 4.0 * tr_xy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xy_zzz[i] = -4.0 * tr_y_xzzz[i] * tke_0 - 6.0 * tr_xy_zzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tke_0 + 4.0 * tr_xy_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 20-30 components of targeted buffer : DF

    auto tr_0_0_xx_xz_xxx = pbuffer.data(idx_op_geom_020_df + 20);

    auto tr_0_0_xx_xz_xxy = pbuffer.data(idx_op_geom_020_df + 21);

    auto tr_0_0_xx_xz_xxz = pbuffer.data(idx_op_geom_020_df + 22);

    auto tr_0_0_xx_xz_xyy = pbuffer.data(idx_op_geom_020_df + 23);

    auto tr_0_0_xx_xz_xyz = pbuffer.data(idx_op_geom_020_df + 24);

    auto tr_0_0_xx_xz_xzz = pbuffer.data(idx_op_geom_020_df + 25);

    auto tr_0_0_xx_xz_yyy = pbuffer.data(idx_op_geom_020_df + 26);

    auto tr_0_0_xx_xz_yyz = pbuffer.data(idx_op_geom_020_df + 27);

    auto tr_0_0_xx_xz_yzz = pbuffer.data(idx_op_geom_020_df + 28);

    auto tr_0_0_xx_xz_zzz = pbuffer.data(idx_op_geom_020_df + 29);

    #pragma omp simd aligned(tr_0_0_xx_xz_xxx, tr_0_0_xx_xz_xxy, tr_0_0_xx_xz_xxz, tr_0_0_xx_xz_xyy, tr_0_0_xx_xz_xyz, tr_0_0_xx_xz_xzz, tr_0_0_xx_xz_yyy, tr_0_0_xx_xz_yyz, tr_0_0_xx_xz_yzz, tr_0_0_xx_xz_zzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xx, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xz_x, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, tr_z_xx, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xz_xxx[i] = 6.0 * tr_z_xx[i] - 4.0 * tr_z_xxxx[i] * tke_0 + 6.0 * tr_xz_x[i] - 6.0 * tr_xz_xxx[i] * tbe_0 - 14.0 * tr_xz_xxx[i] * tke_0 + 4.0 * tr_xz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xxz_xx[i] * tbe_0 + 8.0 * tr_xxz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xxy[i] = 4.0 * tr_z_xy[i] - 4.0 * tr_z_xxxy[i] * tke_0 + 2.0 * tr_xz_y[i] - 6.0 * tr_xz_xxy[i] * tbe_0 - 10.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxz_xy[i] * tbe_0 + 8.0 * tr_xxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xxz[i] = 4.0 * tr_z_xz[i] - 4.0 * tr_z_xxxz[i] * tke_0 + 2.0 * tr_xz_z[i] - 6.0 * tr_xz_xxz[i] * tbe_0 - 10.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxz_xz[i] * tbe_0 + 8.0 * tr_xxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xyy[i] = 2.0 * tr_z_yy[i] - 4.0 * tr_z_xxyy[i] * tke_0 - 6.0 * tr_xz_xyy[i] * tbe_0 - 6.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxz_yy[i] * tbe_0 + 8.0 * tr_xxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xyz[i] = 2.0 * tr_z_yz[i] - 4.0 * tr_z_xxyz[i] * tke_0 - 6.0 * tr_xz_xyz[i] * tbe_0 - 6.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_yz[i] * tbe_0 + 8.0 * tr_xxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_xzz[i] = 2.0 * tr_z_zz[i] - 4.0 * tr_z_xxzz[i] * tke_0 - 6.0 * tr_xz_xzz[i] * tbe_0 - 6.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_zz[i] * tbe_0 + 8.0 * tr_xxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_yyy[i] = -4.0 * tr_z_xyyy[i] * tke_0 - 6.0 * tr_xz_yyy[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tke_0 + 4.0 * tr_xz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_yyz[i] = -4.0 * tr_z_xyyz[i] * tke_0 - 6.0 * tr_xz_yyz[i] * tbe_0 - 2.0 * tr_xz_yyz[i] * tke_0 + 4.0 * tr_xz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_yzz[i] = -4.0 * tr_z_xyzz[i] * tke_0 - 6.0 * tr_xz_yzz[i] * tbe_0 - 2.0 * tr_xz_yzz[i] * tke_0 + 4.0 * tr_xz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xz_zzz[i] = -4.0 * tr_z_xzzz[i] * tke_0 - 6.0 * tr_xz_zzz[i] * tbe_0 - 2.0 * tr_xz_zzz[i] * tke_0 + 4.0 * tr_xz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-40 components of targeted buffer : DF

    auto tr_0_0_xx_yy_xxx = pbuffer.data(idx_op_geom_020_df + 30);

    auto tr_0_0_xx_yy_xxy = pbuffer.data(idx_op_geom_020_df + 31);

    auto tr_0_0_xx_yy_xxz = pbuffer.data(idx_op_geom_020_df + 32);

    auto tr_0_0_xx_yy_xyy = pbuffer.data(idx_op_geom_020_df + 33);

    auto tr_0_0_xx_yy_xyz = pbuffer.data(idx_op_geom_020_df + 34);

    auto tr_0_0_xx_yy_xzz = pbuffer.data(idx_op_geom_020_df + 35);

    auto tr_0_0_xx_yy_yyy = pbuffer.data(idx_op_geom_020_df + 36);

    auto tr_0_0_xx_yy_yyz = pbuffer.data(idx_op_geom_020_df + 37);

    auto tr_0_0_xx_yy_yzz = pbuffer.data(idx_op_geom_020_df + 38);

    auto tr_0_0_xx_yy_zzz = pbuffer.data(idx_op_geom_020_df + 39);

    #pragma omp simd aligned(tr_0_0_xx_yy_xxx, tr_0_0_xx_yy_xxy, tr_0_0_xx_yy_xxz, tr_0_0_xx_yy_xyy, tr_0_0_xx_yy_xyz, tr_0_0_xx_yy_xzz, tr_0_0_xx_yy_yyy, tr_0_0_xx_yy_yyz, tr_0_0_xx_yy_yzz, tr_0_0_xx_yy_zzz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xyy_xx, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yy_xxx[i] = 6.0 * tr_yy_x[i] - 2.0 * tr_yy_xxx[i] * tbe_0 - 14.0 * tr_yy_xxx[i] * tke_0 + 4.0 * tr_yy_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xyy_xx[i] * tbe_0 + 8.0 * tr_xyy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xxy[i] = 2.0 * tr_yy_y[i] - 2.0 * tr_yy_xxy[i] * tbe_0 - 10.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyy_xy[i] * tbe_0 + 8.0 * tr_xyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xxz[i] = 2.0 * tr_yy_z[i] - 2.0 * tr_yy_xxz[i] * tbe_0 - 10.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyy_xz[i] * tbe_0 + 8.0 * tr_xyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xyy[i] = -2.0 * tr_yy_xyy[i] * tbe_0 - 6.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyy_yy[i] * tbe_0 + 8.0 * tr_xyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xyz[i] = -2.0 * tr_yy_xyz[i] * tbe_0 - 6.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_yz[i] * tbe_0 + 8.0 * tr_xyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_xzz[i] = -2.0 * tr_yy_xzz[i] * tbe_0 - 6.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_zz[i] * tbe_0 + 8.0 * tr_xyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_yyy[i] = -2.0 * tr_yy_yyy[i] * tbe_0 - 2.0 * tr_yy_yyy[i] * tke_0 + 4.0 * tr_yy_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_yyz[i] = -2.0 * tr_yy_yyz[i] * tbe_0 - 2.0 * tr_yy_yyz[i] * tke_0 + 4.0 * tr_yy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_yzz[i] = -2.0 * tr_yy_yzz[i] * tbe_0 - 2.0 * tr_yy_yzz[i] * tke_0 + 4.0 * tr_yy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yy_zzz[i] = -2.0 * tr_yy_zzz[i] * tbe_0 - 2.0 * tr_yy_zzz[i] * tke_0 + 4.0 * tr_yy_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 40-50 components of targeted buffer : DF

    auto tr_0_0_xx_yz_xxx = pbuffer.data(idx_op_geom_020_df + 40);

    auto tr_0_0_xx_yz_xxy = pbuffer.data(idx_op_geom_020_df + 41);

    auto tr_0_0_xx_yz_xxz = pbuffer.data(idx_op_geom_020_df + 42);

    auto tr_0_0_xx_yz_xyy = pbuffer.data(idx_op_geom_020_df + 43);

    auto tr_0_0_xx_yz_xyz = pbuffer.data(idx_op_geom_020_df + 44);

    auto tr_0_0_xx_yz_xzz = pbuffer.data(idx_op_geom_020_df + 45);

    auto tr_0_0_xx_yz_yyy = pbuffer.data(idx_op_geom_020_df + 46);

    auto tr_0_0_xx_yz_yyz = pbuffer.data(idx_op_geom_020_df + 47);

    auto tr_0_0_xx_yz_yzz = pbuffer.data(idx_op_geom_020_df + 48);

    auto tr_0_0_xx_yz_zzz = pbuffer.data(idx_op_geom_020_df + 49);

    #pragma omp simd aligned(tr_0_0_xx_yz_xxx, tr_0_0_xx_yz_xxy, tr_0_0_xx_yz_xxz, tr_0_0_xx_yz_xyy, tr_0_0_xx_yz_xyz, tr_0_0_xx_yz_xzz, tr_0_0_xx_yz_yyy, tr_0_0_xx_yz_yyz, tr_0_0_xx_yz_yzz, tr_0_0_xx_yz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yz_xxx[i] = 6.0 * tr_yz_x[i] - 2.0 * tr_yz_xxx[i] * tbe_0 - 14.0 * tr_yz_xxx[i] * tke_0 + 4.0 * tr_yz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xyz_xx[i] * tbe_0 + 8.0 * tr_xyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xxy[i] = 2.0 * tr_yz_y[i] - 2.0 * tr_yz_xxy[i] * tbe_0 - 10.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyz_xy[i] * tbe_0 + 8.0 * tr_xyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xxz[i] = 2.0 * tr_yz_z[i] - 2.0 * tr_yz_xxz[i] * tbe_0 - 10.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyz_xz[i] * tbe_0 + 8.0 * tr_xyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xyy[i] = -2.0 * tr_yz_xyy[i] * tbe_0 - 6.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 8.0 * tr_xyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xyz[i] = -2.0 * tr_yz_xyz[i] * tbe_0 - 6.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 8.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_xzz[i] = -2.0 * tr_yz_xzz[i] * tbe_0 - 6.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 + 8.0 * tr_xyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_yyy[i] = -2.0 * tr_yz_yyy[i] * tbe_0 - 2.0 * tr_yz_yyy[i] * tke_0 + 4.0 * tr_yz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_yyz[i] = -2.0 * tr_yz_yyz[i] * tbe_0 - 2.0 * tr_yz_yyz[i] * tke_0 + 4.0 * tr_yz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_yzz[i] = -2.0 * tr_yz_yzz[i] * tbe_0 - 2.0 * tr_yz_yzz[i] * tke_0 + 4.0 * tr_yz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yz_zzz[i] = -2.0 * tr_yz_zzz[i] * tbe_0 - 2.0 * tr_yz_zzz[i] * tke_0 + 4.0 * tr_yz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 50-60 components of targeted buffer : DF

    auto tr_0_0_xx_zz_xxx = pbuffer.data(idx_op_geom_020_df + 50);

    auto tr_0_0_xx_zz_xxy = pbuffer.data(idx_op_geom_020_df + 51);

    auto tr_0_0_xx_zz_xxz = pbuffer.data(idx_op_geom_020_df + 52);

    auto tr_0_0_xx_zz_xyy = pbuffer.data(idx_op_geom_020_df + 53);

    auto tr_0_0_xx_zz_xyz = pbuffer.data(idx_op_geom_020_df + 54);

    auto tr_0_0_xx_zz_xzz = pbuffer.data(idx_op_geom_020_df + 55);

    auto tr_0_0_xx_zz_yyy = pbuffer.data(idx_op_geom_020_df + 56);

    auto tr_0_0_xx_zz_yyz = pbuffer.data(idx_op_geom_020_df + 57);

    auto tr_0_0_xx_zz_yzz = pbuffer.data(idx_op_geom_020_df + 58);

    auto tr_0_0_xx_zz_zzz = pbuffer.data(idx_op_geom_020_df + 59);

    #pragma omp simd aligned(tr_0_0_xx_zz_xxx, tr_0_0_xx_zz_xxy, tr_0_0_xx_zz_xxz, tr_0_0_xx_zz_xyy, tr_0_0_xx_zz_xyz, tr_0_0_xx_zz_xzz, tr_0_0_xx_zz_yyy, tr_0_0_xx_zz_yyz, tr_0_0_xx_zz_yzz, tr_0_0_xx_zz_zzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xzz_xx, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_zz_xxx[i] = 6.0 * tr_zz_x[i] - 2.0 * tr_zz_xxx[i] * tbe_0 - 14.0 * tr_zz_xxx[i] * tke_0 + 4.0 * tr_zz_xxxxx[i] * tke_0 * tke_0 - 12.0 * tr_xzz_xx[i] * tbe_0 + 8.0 * tr_xzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xxy[i] = 2.0 * tr_zz_y[i] - 2.0 * tr_zz_xxy[i] * tbe_0 - 10.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xzz_xy[i] * tbe_0 + 8.0 * tr_xzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xxz[i] = 2.0 * tr_zz_z[i] - 2.0 * tr_zz_xxz[i] * tbe_0 - 10.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xzz_xz[i] * tbe_0 + 8.0 * tr_xzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xyy[i] = -2.0 * tr_zz_xyy[i] * tbe_0 - 6.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xzz_yy[i] * tbe_0 + 8.0 * tr_xzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xyz[i] = -2.0 * tr_zz_xyz[i] * tbe_0 - 6.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_yz[i] * tbe_0 + 8.0 * tr_xzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_xzz[i] = -2.0 * tr_zz_xzz[i] * tbe_0 - 6.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_zz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_zz[i] * tbe_0 + 8.0 * tr_xzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_yyy[i] = -2.0 * tr_zz_yyy[i] * tbe_0 - 2.0 * tr_zz_yyy[i] * tke_0 + 4.0 * tr_zz_xxyyy[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_yyz[i] = -2.0 * tr_zz_yyz[i] * tbe_0 - 2.0 * tr_zz_yyz[i] * tke_0 + 4.0 * tr_zz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_yzz[i] = -2.0 * tr_zz_yzz[i] * tbe_0 - 2.0 * tr_zz_yzz[i] * tke_0 + 4.0 * tr_zz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zz_zzz[i] = -2.0 * tr_zz_zzz[i] * tbe_0 - 2.0 * tr_zz_zzz[i] * tke_0 + 4.0 * tr_zz_xxzzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-70 components of targeted buffer : DF

    auto tr_0_0_xy_xx_xxx = pbuffer.data(idx_op_geom_020_df + 60);

    auto tr_0_0_xy_xx_xxy = pbuffer.data(idx_op_geom_020_df + 61);

    auto tr_0_0_xy_xx_xxz = pbuffer.data(idx_op_geom_020_df + 62);

    auto tr_0_0_xy_xx_xyy = pbuffer.data(idx_op_geom_020_df + 63);

    auto tr_0_0_xy_xx_xyz = pbuffer.data(idx_op_geom_020_df + 64);

    auto tr_0_0_xy_xx_xzz = pbuffer.data(idx_op_geom_020_df + 65);

    auto tr_0_0_xy_xx_yyy = pbuffer.data(idx_op_geom_020_df + 66);

    auto tr_0_0_xy_xx_yyz = pbuffer.data(idx_op_geom_020_df + 67);

    auto tr_0_0_xy_xx_yzz = pbuffer.data(idx_op_geom_020_df + 68);

    auto tr_0_0_xy_xx_zzz = pbuffer.data(idx_op_geom_020_df + 69);

    #pragma omp simd aligned(tr_0_0_xy_xx_xxx, tr_0_0_xy_xx_xxy, tr_0_0_xy_xx_xxz, tr_0_0_xy_xx_xyy, tr_0_0_xy_xx_xyz, tr_0_0_xy_xx_xzz, tr_0_0_xy_xx_yyy, tr_0_0_xy_xx_yyz, tr_0_0_xy_xx_yzz, tr_0_0_xy_xx_zzz, tr_x_xx, tr_x_xxxy, tr_x_xxyy, tr_x_xxyz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_yy, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_xx_x, tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_y, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xxx_xx, tr_xxx_xxxy, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_yy, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xx, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xx_xxx[i] = -4.0 * tr_x_xxxy[i] * tke_0 - 4.0 * tr_xy_xxx[i] * tbe_0 - 6.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xxy[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_x_xxyy[i] * tke_0 - 4.0 * tr_xy_xxy[i] * tbe_0 + 2.0 * tr_xx_x[i] - 4.0 * tr_xx_xyy[i] * tke_0 - 2.0 * tr_xx_xxx[i] * tke_0 + 4.0 * tr_xx_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xxz[i] = -4.0 * tr_x_xxyz[i] * tke_0 - 4.0 * tr_xy_xxz[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xyy[i] = 4.0 * tr_x_xy[i] - 4.0 * tr_x_xyyy[i] * tke_0 - 4.0 * tr_xy_xyy[i] * tbe_0 + 2.0 * tr_xx_y[i] - 2.0 * tr_xx_yyy[i] * tke_0 - 4.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxx_xy[i] * tbe_0 + 4.0 * tr_xxx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xyz[i] = 2.0 * tr_x_xz[i] - 4.0 * tr_x_xyyz[i] * tke_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + tr_xx_z[i] - 2.0 * tr_xx_yyz[i] * tke_0 - 2.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xz[i] * tbe_0 + 4.0 * tr_xxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_xzz[i] = -4.0 * tr_x_xyzz[i] * tke_0 - 4.0 * tr_xy_xzz[i] * tbe_0 - 2.0 * tr_xx_yzz[i] * tke_0 + 4.0 * tr_xx_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_yyy[i] = 6.0 * tr_x_yy[i] - 4.0 * tr_x_yyyy[i] * tke_0 - 4.0 * tr_xy_yyy[i] * tbe_0 - 6.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xx_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxx_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_yyz[i] = 4.0 * tr_x_yz[i] - 4.0 * tr_x_yyyz[i] * tke_0 - 4.0 * tr_xy_yyz[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxx_yz[i] * tbe_0 + 4.0 * tr_xxx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_yzz[i] = 2.0 * tr_x_zz[i] - 4.0 * tr_x_yyzz[i] * tke_0 - 4.0 * tr_xy_yzz[i] * tbe_0 - 2.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xx_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xx_zzz[i] = -4.0 * tr_x_yzzz[i] * tke_0 - 4.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xx_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 70-80 components of targeted buffer : DF

    auto tr_0_0_xy_xy_xxx = pbuffer.data(idx_op_geom_020_df + 70);

    auto tr_0_0_xy_xy_xxy = pbuffer.data(idx_op_geom_020_df + 71);

    auto tr_0_0_xy_xy_xxz = pbuffer.data(idx_op_geom_020_df + 72);

    auto tr_0_0_xy_xy_xyy = pbuffer.data(idx_op_geom_020_df + 73);

    auto tr_0_0_xy_xy_xyz = pbuffer.data(idx_op_geom_020_df + 74);

    auto tr_0_0_xy_xy_xzz = pbuffer.data(idx_op_geom_020_df + 75);

    auto tr_0_0_xy_xy_yyy = pbuffer.data(idx_op_geom_020_df + 76);

    auto tr_0_0_xy_xy_yyz = pbuffer.data(idx_op_geom_020_df + 77);

    auto tr_0_0_xy_xy_yzz = pbuffer.data(idx_op_geom_020_df + 78);

    auto tr_0_0_xy_xy_zzz = pbuffer.data(idx_op_geom_020_df + 79);

    #pragma omp simd aligned(tr_0_0_xy_xy_xxx, tr_0_0_xy_xy_xxy, tr_0_0_xy_xy_xxz, tr_0_0_xy_xy_xyy, tr_0_0_xy_xy_xyz, tr_0_0_xy_xy_xzz, tr_0_0_xy_xy_yyy, tr_0_0_xy_xy_yyz, tr_0_0_xy_xy_yzz, tr_0_0_xy_xy_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xx, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxy_xx, tr_xxy_xxxy, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_yy, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xy_x, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_y, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xyy_xx, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_y_xx, tr_y_xxxy, tr_y_xxyy, tr_y_xxyz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_yy, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xy_xxx[i] = tr_0_xxx[i] - 2.0 * tr_y_xxxy[i] * tke_0 - 2.0 * tr_yy_xxx[i] * tbe_0 + 3.0 * tr_x_xx[i] - 2.0 * tr_x_xxxx[i] * tke_0 - 6.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyy_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xxy[i] = tr_0_xxy[i] + tr_y_xx[i] - 2.0 * tr_y_xxyy[i] * tke_0 - 2.0 * tr_yy_xxy[i] * tbe_0 + 2.0 * tr_x_xy[i] - 2.0 * tr_x_xxxy[i] * tke_0 + 2.0 * tr_xy_x[i] - 4.0 * tr_xy_xyy[i] * tke_0 - 2.0 * tr_xy_xxx[i] * tke_0 + 4.0 * tr_xy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxy[i] * tbe_0 - 2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xxz[i] = tr_0_xxz[i] - 2.0 * tr_y_xxyz[i] * tke_0 - 2.0 * tr_yy_xxz[i] * tbe_0 + 2.0 * tr_x_xz[i] - 2.0 * tr_x_xxxz[i] * tke_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxz[i] * tbe_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xyy[i] = tr_0_xyy[i] + 2.0 * tr_y_xy[i] - 2.0 * tr_y_xyyy[i] * tke_0 - 2.0 * tr_yy_xyy[i] * tbe_0 + tr_x_yy[i] - 2.0 * tr_x_xxyy[i] * tke_0 + 2.0 * tr_xy_y[i] - 2.0 * tr_xy_yyy[i] * tke_0 - 4.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyy[i] * tbe_0 - 4.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xyz[i] = tr_0_xyz[i] + tr_y_xz[i] - 2.0 * tr_y_xyyz[i] * tke_0 - 2.0 * tr_yy_xyz[i] * tbe_0 + tr_x_yz[i] - 2.0 * tr_x_xxyz[i] * tke_0 + tr_xy_z[i] - 2.0 * tr_xy_yyz[i] * tke_0 - 2.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyz[i] * tbe_0 - 2.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_xzz[i] = tr_0_xzz[i] - 2.0 * tr_y_xyzz[i] * tke_0 - 2.0 * tr_yy_xzz[i] * tbe_0 + tr_x_zz[i] - 2.0 * tr_x_xxzz[i] * tke_0 - 2.0 * tr_xy_yzz[i] * tke_0 + 4.0 * tr_xy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xzz[i] * tbe_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_yyy[i] = tr_0_yyy[i] + 3.0 * tr_y_yy[i] - 2.0 * tr_y_yyyy[i] * tke_0 - 2.0 * tr_yy_yyy[i] * tbe_0 - 2.0 * tr_x_xyyy[i] * tke_0 - 6.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xy_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyy[i] * tbe_0 - 6.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_yyz[i] = tr_0_yyz[i] + 2.0 * tr_y_yz[i] - 2.0 * tr_y_yyyz[i] * tke_0 - 2.0 * tr_yy_yyz[i] * tbe_0 - 2.0 * tr_x_xyyz[i] * tke_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyz[i] * tbe_0 - 4.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_yzz[i] = tr_0_yzz[i] + tr_y_zz[i] - 2.0 * tr_y_yyzz[i] * tke_0 - 2.0 * tr_yy_yzz[i] * tbe_0 - 2.0 * tr_x_xyzz[i] * tke_0 - 2.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yzz[i] * tbe_0 - 2.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xy_zzz[i] = tr_0_zzz[i] - 2.0 * tr_y_yzzz[i] * tke_0 - 2.0 * tr_yy_zzz[i] * tbe_0 - 2.0 * tr_x_xzzz[i] * tke_0 + 4.0 * tr_xy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zzz[i] * tbe_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 80-90 components of targeted buffer : DF

    auto tr_0_0_xy_xz_xxx = pbuffer.data(idx_op_geom_020_df + 80);

    auto tr_0_0_xy_xz_xxy = pbuffer.data(idx_op_geom_020_df + 81);

    auto tr_0_0_xy_xz_xxz = pbuffer.data(idx_op_geom_020_df + 82);

    auto tr_0_0_xy_xz_xyy = pbuffer.data(idx_op_geom_020_df + 83);

    auto tr_0_0_xy_xz_xyz = pbuffer.data(idx_op_geom_020_df + 84);

    auto tr_0_0_xy_xz_xzz = pbuffer.data(idx_op_geom_020_df + 85);

    auto tr_0_0_xy_xz_yyy = pbuffer.data(idx_op_geom_020_df + 86);

    auto tr_0_0_xy_xz_yyz = pbuffer.data(idx_op_geom_020_df + 87);

    auto tr_0_0_xy_xz_yzz = pbuffer.data(idx_op_geom_020_df + 88);

    auto tr_0_0_xy_xz_zzz = pbuffer.data(idx_op_geom_020_df + 89);

    #pragma omp simd aligned(tr_0_0_xy_xz_xxx, tr_0_0_xy_xz_xxy, tr_0_0_xy_xz_xxz, tr_0_0_xy_xz_xyy, tr_0_0_xy_xz_xyz, tr_0_0_xy_xz_xzz, tr_0_0_xy_xz_yyy, tr_0_0_xy_xz_yyz, tr_0_0_xy_xz_yzz, tr_0_0_xy_xz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxz_xx, tr_xxz_xxxy, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_yy, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xz_x, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_y, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_z_xx, tr_z_xxxy, tr_z_xxyy, tr_z_xxyz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_yy, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xz_xxx[i] = -2.0 * tr_z_xxxy[i] * tke_0 - 2.0 * tr_yz_xxx[i] * tbe_0 - 6.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xxy[i] = tr_z_xx[i] - 2.0 * tr_z_xxyy[i] * tke_0 - 2.0 * tr_yz_xxy[i] * tbe_0 + 2.0 * tr_xz_x[i] - 4.0 * tr_xz_xyy[i] * tke_0 - 2.0 * tr_xz_xxx[i] * tke_0 + 4.0 * tr_xz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xxz[i] = -2.0 * tr_z_xxyz[i] * tke_0 - 2.0 * tr_yz_xxz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xyy[i] = 2.0 * tr_z_xy[i] - 2.0 * tr_z_xyyy[i] * tke_0 - 2.0 * tr_yz_xyy[i] * tbe_0 + 2.0 * tr_xz_y[i] - 2.0 * tr_xz_yyy[i] * tke_0 - 4.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xyz[i] = tr_z_xz[i] - 2.0 * tr_z_xyyz[i] * tke_0 - 2.0 * tr_yz_xyz[i] * tbe_0 + tr_xz_z[i] - 2.0 * tr_xz_yyz[i] * tke_0 - 2.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_xzz[i] = -2.0 * tr_z_xyzz[i] * tke_0 - 2.0 * tr_yz_xzz[i] * tbe_0 - 2.0 * tr_xz_yzz[i] * tke_0 + 4.0 * tr_xz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_yyy[i] = 3.0 * tr_z_yy[i] - 2.0 * tr_z_yyyy[i] * tke_0 - 2.0 * tr_yz_yyy[i] * tbe_0 - 6.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_yyz[i] = 2.0 * tr_z_yz[i] - 2.0 * tr_z_yyyz[i] * tke_0 - 2.0 * tr_yz_yyz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_yzz[i] = tr_z_zz[i] - 2.0 * tr_z_yyzz[i] * tke_0 - 2.0 * tr_yz_yzz[i] * tbe_0 - 2.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xz_zzz[i] = -2.0 * tr_z_yzzz[i] * tke_0 - 2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_xz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-100 components of targeted buffer : DF

    auto tr_0_0_xy_yy_xxx = pbuffer.data(idx_op_geom_020_df + 90);

    auto tr_0_0_xy_yy_xxy = pbuffer.data(idx_op_geom_020_df + 91);

    auto tr_0_0_xy_yy_xxz = pbuffer.data(idx_op_geom_020_df + 92);

    auto tr_0_0_xy_yy_xyy = pbuffer.data(idx_op_geom_020_df + 93);

    auto tr_0_0_xy_yy_xyz = pbuffer.data(idx_op_geom_020_df + 94);

    auto tr_0_0_xy_yy_xzz = pbuffer.data(idx_op_geom_020_df + 95);

    auto tr_0_0_xy_yy_yyy = pbuffer.data(idx_op_geom_020_df + 96);

    auto tr_0_0_xy_yy_yyz = pbuffer.data(idx_op_geom_020_df + 97);

    auto tr_0_0_xy_yy_yzz = pbuffer.data(idx_op_geom_020_df + 98);

    auto tr_0_0_xy_yy_zzz = pbuffer.data(idx_op_geom_020_df + 99);

    #pragma omp simd aligned(tr_0_0_xy_yy_xxx, tr_0_0_xy_yy_xxy, tr_0_0_xy_yy_xxz, tr_0_0_xy_yy_xyy, tr_0_0_xy_yy_xyz, tr_0_0_xy_yy_xzz, tr_0_0_xy_yy_yyy, tr_0_0_xy_yy_yyz, tr_0_0_xy_yy_yzz, tr_0_0_xy_yy_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxy, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_yy, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_y_xx, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_y, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yyy_xx, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yy_xxx[i] = 6.0 * tr_y_xx[i] - 4.0 * tr_y_xxxx[i] * tke_0 - 6.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_yyy_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xxy[i] = 4.0 * tr_y_xy[i] - 4.0 * tr_y_xxxy[i] * tke_0 + 2.0 * tr_yy_x[i] - 4.0 * tr_yy_xyy[i] * tke_0 - 2.0 * tr_yy_xxx[i] * tke_0 + 4.0 * tr_yy_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yyy_xy[i] * tbe_0 + 4.0 * tr_yyy_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxy[i] * tbe_0 - 2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xxz[i] = 4.0 * tr_y_xz[i] - 4.0 * tr_y_xxxz[i] * tke_0 - 4.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyy_xz[i] * tbe_0 + 4.0 * tr_yyy_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xyy[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_y_xxyy[i] * tke_0 + 2.0 * tr_yy_y[i] - 2.0 * tr_yy_yyy[i] * tke_0 - 4.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyy_yy[i] * tbe_0 + 4.0 * tr_yyy_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyy[i] * tbe_0 - 4.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xyz[i] = 2.0 * tr_y_yz[i] - 4.0 * tr_y_xxyz[i] * tke_0 + tr_yy_z[i] - 2.0 * tr_yy_yyz[i] * tke_0 - 2.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyy_yz[i] * tbe_0 + 4.0 * tr_yyy_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyz[i] * tbe_0 - 2.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_xzz[i] = 2.0 * tr_y_zz[i] - 4.0 * tr_y_xxzz[i] * tke_0 - 2.0 * tr_yy_yzz[i] * tke_0 + 4.0 * tr_yy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyy_zz[i] * tbe_0 + 4.0 * tr_yyy_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_yyy[i] = -4.0 * tr_y_xyyy[i] * tke_0 - 6.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yy_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyy[i] * tbe_0 - 6.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_yyz[i] = -4.0 * tr_y_xyyz[i] * tke_0 - 4.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyz[i] * tbe_0 - 4.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_yzz[i] = -4.0 * tr_y_xyzz[i] * tke_0 - 2.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yzz[i] * tbe_0 - 2.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yy_zzz[i] = -4.0 * tr_y_xzzz[i] * tke_0 + 4.0 * tr_yy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyy_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 100-110 components of targeted buffer : DF

    auto tr_0_0_xy_yz_xxx = pbuffer.data(idx_op_geom_020_df + 100);

    auto tr_0_0_xy_yz_xxy = pbuffer.data(idx_op_geom_020_df + 101);

    auto tr_0_0_xy_yz_xxz = pbuffer.data(idx_op_geom_020_df + 102);

    auto tr_0_0_xy_yz_xyy = pbuffer.data(idx_op_geom_020_df + 103);

    auto tr_0_0_xy_yz_xyz = pbuffer.data(idx_op_geom_020_df + 104);

    auto tr_0_0_xy_yz_xzz = pbuffer.data(idx_op_geom_020_df + 105);

    auto tr_0_0_xy_yz_yyy = pbuffer.data(idx_op_geom_020_df + 106);

    auto tr_0_0_xy_yz_yyz = pbuffer.data(idx_op_geom_020_df + 107);

    auto tr_0_0_xy_yz_yzz = pbuffer.data(idx_op_geom_020_df + 108);

    auto tr_0_0_xy_yz_zzz = pbuffer.data(idx_op_geom_020_df + 109);

    #pragma omp simd aligned(tr_0_0_xy_yz_xxx, tr_0_0_xy_yz_xxy, tr_0_0_xy_yz_xxz, tr_0_0_xy_yz_xyy, tr_0_0_xy_yz_xyz, tr_0_0_xy_yz_xzz, tr_0_0_xy_yz_yyy, tr_0_0_xy_yz_yyz, tr_0_0_xy_yz_yzz, tr_0_0_xy_yz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_yyz_xx, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_y, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_z_xx, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yz_xxx[i] = 3.0 * tr_z_xx[i] - 2.0 * tr_z_xxxx[i] * tke_0 - 6.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xxy[i] = 2.0 * tr_z_xy[i] - 2.0 * tr_z_xxxy[i] * tke_0 + 2.0 * tr_yz_x[i] - 4.0 * tr_yz_xyy[i] * tke_0 - 2.0 * tr_yz_xxx[i] * tke_0 + 4.0 * tr_yz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxy[i] * tbe_0 - 2.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xxz[i] = 2.0 * tr_z_xz[i] - 2.0 * tr_z_xxxz[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xyy[i] = tr_z_yy[i] - 2.0 * tr_z_xxyy[i] * tke_0 + 2.0 * tr_yz_y[i] - 2.0 * tr_yz_yyy[i] * tke_0 - 4.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xyy[i] * tbe_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xyz[i] = tr_z_yz[i] - 2.0 * tr_z_xxyz[i] * tke_0 + tr_yz_z[i] - 2.0 * tr_yz_yyz[i] * tke_0 - 2.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xyz[i] * tbe_0 - 2.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_xzz[i] = tr_z_zz[i] - 2.0 * tr_z_xxzz[i] * tke_0 - 2.0 * tr_yz_yzz[i] * tke_0 + 4.0 * tr_yz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_yyy[i] = -2.0 * tr_z_xyyy[i] * tke_0 - 6.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yyy[i] * tbe_0 - 6.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_yyz[i] = -2.0 * tr_z_xyyz[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yyz[i] * tbe_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_yzz[i] = -2.0 * tr_z_xyzz[i] * tke_0 - 2.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_yzz[i] * tbe_0 - 2.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yz_zzz[i] = -2.0 * tr_z_xzzz[i] * tke_0 + 4.0 * tr_yz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 110-120 components of targeted buffer : DF

    auto tr_0_0_xy_zz_xxx = pbuffer.data(idx_op_geom_020_df + 110);

    auto tr_0_0_xy_zz_xxy = pbuffer.data(idx_op_geom_020_df + 111);

    auto tr_0_0_xy_zz_xxz = pbuffer.data(idx_op_geom_020_df + 112);

    auto tr_0_0_xy_zz_xyy = pbuffer.data(idx_op_geom_020_df + 113);

    auto tr_0_0_xy_zz_xyz = pbuffer.data(idx_op_geom_020_df + 114);

    auto tr_0_0_xy_zz_xzz = pbuffer.data(idx_op_geom_020_df + 115);

    auto tr_0_0_xy_zz_yyy = pbuffer.data(idx_op_geom_020_df + 116);

    auto tr_0_0_xy_zz_yyz = pbuffer.data(idx_op_geom_020_df + 117);

    auto tr_0_0_xy_zz_yzz = pbuffer.data(idx_op_geom_020_df + 118);

    auto tr_0_0_xy_zz_zzz = pbuffer.data(idx_op_geom_020_df + 119);

    #pragma omp simd aligned(tr_0_0_xy_zz_xxx, tr_0_0_xy_zz_xxy, tr_0_0_xy_zz_xxz, tr_0_0_xy_zz_xyy, tr_0_0_xy_zz_xyz, tr_0_0_xy_zz_xzz, tr_0_0_xy_zz_yyy, tr_0_0_xy_zz_yyz, tr_0_0_xy_zz_yzz, tr_0_0_xy_zz_zzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xzz_xx, tr_xzz_xxxy, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_yy, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_yzz_xx, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_y, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_zz_xxx[i] = -6.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxxxy[i] * tke_0 * tke_0 - 6.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xxy[i] = 2.0 * tr_zz_x[i] - 4.0 * tr_zz_xyy[i] * tke_0 - 2.0 * tr_zz_xxx[i] * tke_0 + 4.0 * tr_zz_xxxyy[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xxz[i] = -4.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xyy[i] = 2.0 * tr_zz_y[i] - 2.0 * tr_zz_yyy[i] * tke_0 - 4.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxyyy[i] * tke_0 * tke_0 - 2.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xyz[i] = tr_zz_z[i] - 2.0 * tr_zz_yyz[i] * tke_0 - 2.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_xzz[i] = -2.0 * tr_zz_yzz[i] * tke_0 + 4.0 * tr_zz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_yyy[i] = -6.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zz_xyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_yyz[i] = -4.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_yzz[i] = -2.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_zz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zz_zzz[i] = 4.0 * tr_zz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-130 components of targeted buffer : DF

    auto tr_0_0_xz_xx_xxx = pbuffer.data(idx_op_geom_020_df + 120);

    auto tr_0_0_xz_xx_xxy = pbuffer.data(idx_op_geom_020_df + 121);

    auto tr_0_0_xz_xx_xxz = pbuffer.data(idx_op_geom_020_df + 122);

    auto tr_0_0_xz_xx_xyy = pbuffer.data(idx_op_geom_020_df + 123);

    auto tr_0_0_xz_xx_xyz = pbuffer.data(idx_op_geom_020_df + 124);

    auto tr_0_0_xz_xx_xzz = pbuffer.data(idx_op_geom_020_df + 125);

    auto tr_0_0_xz_xx_yyy = pbuffer.data(idx_op_geom_020_df + 126);

    auto tr_0_0_xz_xx_yyz = pbuffer.data(idx_op_geom_020_df + 127);

    auto tr_0_0_xz_xx_yzz = pbuffer.data(idx_op_geom_020_df + 128);

    auto tr_0_0_xz_xx_zzz = pbuffer.data(idx_op_geom_020_df + 129);

    #pragma omp simd aligned(tr_0_0_xz_xx_xxx, tr_0_0_xz_xx_xxy, tr_0_0_xz_xx_xxz, tr_0_0_xz_xx_xyy, tr_0_0_xz_xx_xyz, tr_0_0_xz_xx_xzz, tr_0_0_xz_xx_yyy, tr_0_0_xz_xx_yyz, tr_0_0_xz_xx_yzz, tr_0_0_xz_xx_zzz, tr_x_xx, tr_x_xxxz, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xx_x, tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_y, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxx_xx, tr_xxx_xxxz, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxx_zzzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xx, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xx_xxx[i] = -4.0 * tr_x_xxxz[i] * tke_0 - 4.0 * tr_xz_xxx[i] * tbe_0 - 6.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xxy[i] = -4.0 * tr_x_xxyz[i] * tke_0 - 4.0 * tr_xz_xxy[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xxz[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_x_xxzz[i] * tke_0 - 4.0 * tr_xz_xxz[i] * tbe_0 + 2.0 * tr_xx_x[i] - 4.0 * tr_xx_xzz[i] * tke_0 - 2.0 * tr_xx_xxx[i] * tke_0 + 4.0 * tr_xx_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xyy[i] = -4.0 * tr_x_xyyz[i] * tke_0 - 4.0 * tr_xz_xyy[i] * tbe_0 - 2.0 * tr_xx_yyz[i] * tke_0 + 4.0 * tr_xx_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xyz[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_x_xyzz[i] * tke_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + tr_xx_y[i] - 2.0 * tr_xx_yzz[i] * tke_0 - 2.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xy[i] * tbe_0 + 4.0 * tr_xxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_xzz[i] = 4.0 * tr_x_xz[i] - 4.0 * tr_x_xzzz[i] * tke_0 - 4.0 * tr_xz_xzz[i] * tbe_0 + 2.0 * tr_xx_z[i] - 2.0 * tr_xx_zzz[i] * tke_0 - 4.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxx_xz[i] * tbe_0 + 4.0 * tr_xxx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_yyy[i] = -4.0 * tr_x_yyyz[i] * tke_0 - 4.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xx_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_yyz[i] = 2.0 * tr_x_yy[i] - 4.0 * tr_x_yyzz[i] * tke_0 - 4.0 * tr_xz_yyz[i] * tbe_0 - 2.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xx_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_yzz[i] = 4.0 * tr_x_yz[i] - 4.0 * tr_x_yzzz[i] * tke_0 - 4.0 * tr_xz_yzz[i] * tbe_0 - 4.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxx_yz[i] * tbe_0 + 4.0 * tr_xxx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xx_zzz[i] = 6.0 * tr_x_zz[i] - 4.0 * tr_x_zzzz[i] * tke_0 - 4.0 * tr_xz_zzz[i] * tbe_0 - 6.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xx_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxx_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 130-140 components of targeted buffer : DF

    auto tr_0_0_xz_xy_xxx = pbuffer.data(idx_op_geom_020_df + 130);

    auto tr_0_0_xz_xy_xxy = pbuffer.data(idx_op_geom_020_df + 131);

    auto tr_0_0_xz_xy_xxz = pbuffer.data(idx_op_geom_020_df + 132);

    auto tr_0_0_xz_xy_xyy = pbuffer.data(idx_op_geom_020_df + 133);

    auto tr_0_0_xz_xy_xyz = pbuffer.data(idx_op_geom_020_df + 134);

    auto tr_0_0_xz_xy_xzz = pbuffer.data(idx_op_geom_020_df + 135);

    auto tr_0_0_xz_xy_yyy = pbuffer.data(idx_op_geom_020_df + 136);

    auto tr_0_0_xz_xy_yyz = pbuffer.data(idx_op_geom_020_df + 137);

    auto tr_0_0_xz_xy_yzz = pbuffer.data(idx_op_geom_020_df + 138);

    auto tr_0_0_xz_xy_zzz = pbuffer.data(idx_op_geom_020_df + 139);

    #pragma omp simd aligned(tr_0_0_xz_xy_xxx, tr_0_0_xz_xy_xxy, tr_0_0_xz_xy_xxz, tr_0_0_xz_xy_xyy, tr_0_0_xz_xy_xyz, tr_0_0_xz_xy_xzz, tr_0_0_xz_xy_yyy, tr_0_0_xz_xy_yyz, tr_0_0_xz_xy_yzz, tr_0_0_xz_xy_zzz, tr_xxy_xx, tr_xxy_xxxz, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xy_x, tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_y, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_y_xx, tr_y_xxxz, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xy_xxx[i] = -2.0 * tr_y_xxxz[i] * tke_0 - 2.0 * tr_yz_xxx[i] * tbe_0 - 6.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xxy[i] = -2.0 * tr_y_xxyz[i] * tke_0 - 2.0 * tr_yz_xxy[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xxz[i] = tr_y_xx[i] - 2.0 * tr_y_xxzz[i] * tke_0 - 2.0 * tr_yz_xxz[i] * tbe_0 + 2.0 * tr_xy_x[i] - 4.0 * tr_xy_xzz[i] * tke_0 - 2.0 * tr_xy_xxx[i] * tke_0 + 4.0 * tr_xy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xyy[i] = -2.0 * tr_y_xyyz[i] * tke_0 - 2.0 * tr_yz_xyy[i] * tbe_0 - 2.0 * tr_xy_yyz[i] * tke_0 + 4.0 * tr_xy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xyz[i] = tr_y_xy[i] - 2.0 * tr_y_xyzz[i] * tke_0 - 2.0 * tr_yz_xyz[i] * tbe_0 + tr_xy_y[i] - 2.0 * tr_xy_yzz[i] * tke_0 - 2.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_xzz[i] = 2.0 * tr_y_xz[i] - 2.0 * tr_y_xzzz[i] * tke_0 - 2.0 * tr_yz_xzz[i] * tbe_0 + 2.0 * tr_xy_z[i] - 2.0 * tr_xy_zzz[i] * tke_0 - 4.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_yyy[i] = -2.0 * tr_y_yyyz[i] * tke_0 - 2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_xy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_yyz[i] = tr_y_yy[i] - 2.0 * tr_y_yyzz[i] * tke_0 - 2.0 * tr_yz_yyz[i] * tbe_0 - 2.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_yzz[i] = 2.0 * tr_y_yz[i] - 2.0 * tr_y_yzzz[i] * tke_0 - 2.0 * tr_yz_yzz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xy_zzz[i] = 3.0 * tr_y_zz[i] - 2.0 * tr_y_zzzz[i] * tke_0 - 2.0 * tr_yz_zzz[i] * tbe_0 - 6.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xy_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 140-150 components of targeted buffer : DF

    auto tr_0_0_xz_xz_xxx = pbuffer.data(idx_op_geom_020_df + 140);

    auto tr_0_0_xz_xz_xxy = pbuffer.data(idx_op_geom_020_df + 141);

    auto tr_0_0_xz_xz_xxz = pbuffer.data(idx_op_geom_020_df + 142);

    auto tr_0_0_xz_xz_xyy = pbuffer.data(idx_op_geom_020_df + 143);

    auto tr_0_0_xz_xz_xyz = pbuffer.data(idx_op_geom_020_df + 144);

    auto tr_0_0_xz_xz_xzz = pbuffer.data(idx_op_geom_020_df + 145);

    auto tr_0_0_xz_xz_yyy = pbuffer.data(idx_op_geom_020_df + 146);

    auto tr_0_0_xz_xz_yyz = pbuffer.data(idx_op_geom_020_df + 147);

    auto tr_0_0_xz_xz_yzz = pbuffer.data(idx_op_geom_020_df + 148);

    auto tr_0_0_xz_xz_zzz = pbuffer.data(idx_op_geom_020_df + 149);

    #pragma omp simd aligned(tr_0_0_xz_xz_xxx, tr_0_0_xz_xz_xxy, tr_0_0_xz_xz_xxz, tr_0_0_xz_xz_xyy, tr_0_0_xz_xz_xyz, tr_0_0_xz_xz_xzz, tr_0_0_xz_xz_yyy, tr_0_0_xz_xz_yyz, tr_0_0_xz_xz_yzz, tr_0_0_xz_xz_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xx, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxz_xx, tr_xxz_xxxz, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xz_x, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_y, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_z_xx, tr_z_xxxz, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xz_xxx[i] = tr_0_xxx[i] - 2.0 * tr_z_xxxz[i] * tke_0 - 2.0 * tr_zz_xxx[i] * tbe_0 + 3.0 * tr_x_xx[i] - 2.0 * tr_x_xxxx[i] * tke_0 - 6.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xxy[i] = tr_0_xxy[i] - 2.0 * tr_z_xxyz[i] * tke_0 - 2.0 * tr_zz_xxy[i] * tbe_0 + 2.0 * tr_x_xy[i] - 2.0 * tr_x_xxxy[i] * tke_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxy[i] * tbe_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xxz[i] = tr_0_xxz[i] + tr_z_xx[i] - 2.0 * tr_z_xxzz[i] * tke_0 - 2.0 * tr_zz_xxz[i] * tbe_0 + 2.0 * tr_x_xz[i] - 2.0 * tr_x_xxxz[i] * tke_0 + 2.0 * tr_xz_x[i] - 4.0 * tr_xz_xzz[i] * tke_0 - 2.0 * tr_xz_xxx[i] * tke_0 + 4.0 * tr_xz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxz[i] * tbe_0 - 2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xyy[i] = tr_0_xyy[i] - 2.0 * tr_z_xyyz[i] * tke_0 - 2.0 * tr_zz_xyy[i] * tbe_0 + tr_x_yy[i] - 2.0 * tr_x_xxyy[i] * tke_0 - 2.0 * tr_xz_yyz[i] * tke_0 + 4.0 * tr_xz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyy[i] * tbe_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xyz[i] = tr_0_xyz[i] + tr_z_xy[i] - 2.0 * tr_z_xyzz[i] * tke_0 - 2.0 * tr_zz_xyz[i] * tbe_0 + tr_x_yz[i] - 2.0 * tr_x_xxyz[i] * tke_0 + tr_xz_y[i] - 2.0 * tr_xz_yzz[i] * tke_0 - 2.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyz[i] * tbe_0 - 2.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_xzz[i] = tr_0_xzz[i] + 2.0 * tr_z_xz[i] - 2.0 * tr_z_xzzz[i] * tke_0 - 2.0 * tr_zz_xzz[i] * tbe_0 + tr_x_zz[i] - 2.0 * tr_x_xxzz[i] * tke_0 + 2.0 * tr_xz_z[i] - 2.0 * tr_xz_zzz[i] * tke_0 - 4.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xzz[i] * tbe_0 - 4.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_yyy[i] = tr_0_yyy[i] - 2.0 * tr_z_yyyz[i] * tke_0 - 2.0 * tr_zz_yyy[i] * tbe_0 - 2.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_xz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyy[i] * tbe_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_yyz[i] = tr_0_yyz[i] + tr_z_yy[i] - 2.0 * tr_z_yyzz[i] * tke_0 - 2.0 * tr_zz_yyz[i] * tbe_0 - 2.0 * tr_x_xyyz[i] * tke_0 - 2.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyz[i] * tbe_0 - 2.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_yzz[i] = tr_0_yzz[i] + 2.0 * tr_z_yz[i] - 2.0 * tr_z_yzzz[i] * tke_0 - 2.0 * tr_zz_yzz[i] * tbe_0 - 2.0 * tr_x_xyzz[i] * tke_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yzz[i] * tbe_0 - 4.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xz_zzz[i] = tr_0_zzz[i] + 3.0 * tr_z_zz[i] - 2.0 * tr_z_zzzz[i] * tke_0 - 2.0 * tr_zz_zzz[i] * tbe_0 - 2.0 * tr_x_xzzz[i] * tke_0 - 6.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zzz[i] * tbe_0 - 6.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-160 components of targeted buffer : DF

    auto tr_0_0_xz_yy_xxx = pbuffer.data(idx_op_geom_020_df + 150);

    auto tr_0_0_xz_yy_xxy = pbuffer.data(idx_op_geom_020_df + 151);

    auto tr_0_0_xz_yy_xxz = pbuffer.data(idx_op_geom_020_df + 152);

    auto tr_0_0_xz_yy_xyy = pbuffer.data(idx_op_geom_020_df + 153);

    auto tr_0_0_xz_yy_xyz = pbuffer.data(idx_op_geom_020_df + 154);

    auto tr_0_0_xz_yy_xzz = pbuffer.data(idx_op_geom_020_df + 155);

    auto tr_0_0_xz_yy_yyy = pbuffer.data(idx_op_geom_020_df + 156);

    auto tr_0_0_xz_yy_yyz = pbuffer.data(idx_op_geom_020_df + 157);

    auto tr_0_0_xz_yy_yzz = pbuffer.data(idx_op_geom_020_df + 158);

    auto tr_0_0_xz_yy_zzz = pbuffer.data(idx_op_geom_020_df + 159);

    #pragma omp simd aligned(tr_0_0_xz_yy_xxx, tr_0_0_xz_yy_xxy, tr_0_0_xz_yy_xxz, tr_0_0_xz_yy_xyy, tr_0_0_xz_yy_xyz, tr_0_0_xz_yy_xzz, tr_0_0_xz_yy_yyy, tr_0_0_xz_yy_yyz, tr_0_0_xz_yy_yzz, tr_0_0_xz_yy_zzz, tr_xyy_xx, tr_xyy_xxxz, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_yy_x, tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_y, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, tr_yyz_xx, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yy_xxx[i] = -6.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xxy[i] = -4.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xxz[i] = 2.0 * tr_yy_x[i] - 4.0 * tr_yy_xzz[i] * tke_0 - 2.0 * tr_yy_xxx[i] * tke_0 + 4.0 * tr_yy_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xyy[i] = -2.0 * tr_yy_yyz[i] * tke_0 + 4.0 * tr_yy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xyz[i] = tr_yy_y[i] - 2.0 * tr_yy_yzz[i] * tke_0 - 2.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_xzz[i] = 2.0 * tr_yy_z[i] - 2.0 * tr_yy_zzz[i] * tke_0 - 4.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_yyy[i] = 4.0 * tr_yy_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_yyz[i] = -2.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yy_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_yzz[i] = -4.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yy_zzz[i] = -6.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yy_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 160-170 components of targeted buffer : DF

    auto tr_0_0_xz_yz_xxx = pbuffer.data(idx_op_geom_020_df + 160);

    auto tr_0_0_xz_yz_xxy = pbuffer.data(idx_op_geom_020_df + 161);

    auto tr_0_0_xz_yz_xxz = pbuffer.data(idx_op_geom_020_df + 162);

    auto tr_0_0_xz_yz_xyy = pbuffer.data(idx_op_geom_020_df + 163);

    auto tr_0_0_xz_yz_xyz = pbuffer.data(idx_op_geom_020_df + 164);

    auto tr_0_0_xz_yz_xzz = pbuffer.data(idx_op_geom_020_df + 165);

    auto tr_0_0_xz_yz_yyy = pbuffer.data(idx_op_geom_020_df + 166);

    auto tr_0_0_xz_yz_yyz = pbuffer.data(idx_op_geom_020_df + 167);

    auto tr_0_0_xz_yz_yzz = pbuffer.data(idx_op_geom_020_df + 168);

    auto tr_0_0_xz_yz_zzz = pbuffer.data(idx_op_geom_020_df + 169);

    #pragma omp simd aligned(tr_0_0_xz_yz_xxx, tr_0_0_xz_yz_xxy, tr_0_0_xz_yz_xxz, tr_0_0_xz_yz_xyy, tr_0_0_xz_yz_xyz, tr_0_0_xz_yz_xzz, tr_0_0_xz_yz_yyy, tr_0_0_xz_yz_yyz, tr_0_0_xz_yz_yzz, tr_0_0_xz_yz_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_y_xx, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_y, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yz_xxx[i] = 3.0 * tr_y_xx[i] - 2.0 * tr_y_xxxx[i] * tke_0 - 6.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xxy[i] = 2.0 * tr_y_xy[i] - 2.0 * tr_y_xxxy[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xxz[i] = 2.0 * tr_y_xz[i] - 2.0 * tr_y_xxxz[i] * tke_0 + 2.0 * tr_yz_x[i] - 4.0 * tr_yz_xzz[i] * tke_0 - 2.0 * tr_yz_xxx[i] * tke_0 + 4.0 * tr_yz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxz[i] * tbe_0 - 2.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xyy[i] = tr_y_yy[i] - 2.0 * tr_y_xxyy[i] * tke_0 - 2.0 * tr_yz_yyz[i] * tke_0 + 4.0 * tr_yz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xyz[i] = tr_y_yz[i] - 2.0 * tr_y_xxyz[i] * tke_0 + tr_yz_y[i] - 2.0 * tr_yz_yzz[i] * tke_0 - 2.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyz[i] * tbe_0 - 2.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_xzz[i] = tr_y_zz[i] - 2.0 * tr_y_xxzz[i] * tke_0 + 2.0 * tr_yz_z[i] - 2.0 * tr_yz_zzz[i] * tke_0 - 4.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xzz[i] * tbe_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_yyy[i] = -2.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_yz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_yyz[i] = -2.0 * tr_y_xyyz[i] * tke_0 - 2.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyz[i] * tbe_0 - 2.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_yzz[i] = -2.0 * tr_y_xyzz[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yzz[i] * tbe_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yz_zzz[i] = -2.0 * tr_y_xzzz[i] * tke_0 - 6.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_zzz[i] * tbe_0 - 6.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 170-180 components of targeted buffer : DF

    auto tr_0_0_xz_zz_xxx = pbuffer.data(idx_op_geom_020_df + 170);

    auto tr_0_0_xz_zz_xxy = pbuffer.data(idx_op_geom_020_df + 171);

    auto tr_0_0_xz_zz_xxz = pbuffer.data(idx_op_geom_020_df + 172);

    auto tr_0_0_xz_zz_xyy = pbuffer.data(idx_op_geom_020_df + 173);

    auto tr_0_0_xz_zz_xyz = pbuffer.data(idx_op_geom_020_df + 174);

    auto tr_0_0_xz_zz_xzz = pbuffer.data(idx_op_geom_020_df + 175);

    auto tr_0_0_xz_zz_yyy = pbuffer.data(idx_op_geom_020_df + 176);

    auto tr_0_0_xz_zz_yyz = pbuffer.data(idx_op_geom_020_df + 177);

    auto tr_0_0_xz_zz_yzz = pbuffer.data(idx_op_geom_020_df + 178);

    auto tr_0_0_xz_zz_zzz = pbuffer.data(idx_op_geom_020_df + 179);

    #pragma omp simd aligned(tr_0_0_xz_zz_xxx, tr_0_0_xz_zz_xxy, tr_0_0_xz_zz_xxz, tr_0_0_xz_zz_xyy, tr_0_0_xz_zz_xyz, tr_0_0_xz_zz_xzz, tr_0_0_xz_zz_yyy, tr_0_0_xz_zz_yyz, tr_0_0_xz_zz_yzz, tr_0_0_xz_zz_zzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxz, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_z_xx, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_y, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, tr_zzz_xx, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_zz_xxx[i] = 6.0 * tr_z_xx[i] - 4.0 * tr_z_xxxx[i] * tke_0 - 6.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxxxz[i] * tke_0 * tke_0 - 6.0 * tr_zzz_xx[i] * tbe_0 + 4.0 * tr_zzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xxy[i] = 4.0 * tr_z_xy[i] - 4.0 * tr_z_xxxy[i] * tke_0 - 4.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xxxyz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_xy[i] * tbe_0 + 4.0 * tr_zzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xxz[i] = 4.0 * tr_z_xz[i] - 4.0 * tr_z_xxxz[i] * tke_0 + 2.0 * tr_zz_x[i] - 4.0 * tr_zz_xzz[i] * tke_0 - 2.0 * tr_zz_xxx[i] * tke_0 + 4.0 * tr_zz_xxxzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_xz[i] * tbe_0 + 4.0 * tr_zzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxz[i] * tbe_0 - 2.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xyy[i] = 2.0 * tr_z_yy[i] - 4.0 * tr_z_xxyy[i] * tke_0 - 2.0 * tr_zz_yyz[i] * tke_0 + 4.0 * tr_zz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_yy[i] * tbe_0 + 4.0 * tr_zzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xyz[i] = 2.0 * tr_z_yz[i] - 4.0 * tr_z_xxyz[i] * tke_0 + tr_zz_y[i] - 2.0 * tr_zz_yzz[i] * tke_0 - 2.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_yz[i] * tbe_0 + 4.0 * tr_zzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xyz[i] * tbe_0 - 2.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_xzz[i] = 2.0 * tr_z_zz[i] - 4.0 * tr_z_xxzz[i] * tke_0 + 2.0 * tr_zz_z[i] - 2.0 * tr_zz_zzz[i] * tke_0 - 4.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_zz[i] * tbe_0 + 4.0 * tr_zzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xzz[i] * tbe_0 - 4.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_yyy[i] = -4.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_zz_xyyyz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_yyz[i] = -4.0 * tr_z_xyyz[i] * tke_0 - 2.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zz_xyyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yyz[i] * tbe_0 - 2.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_yzz[i] = -4.0 * tr_z_xyzz[i] * tke_0 - 4.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yzz[i] * tbe_0 - 4.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zz_zzz[i] = -4.0 * tr_z_xzzz[i] * tke_0 - 6.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_zz_xzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_zzz[i] * tbe_0 - 6.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-190 components of targeted buffer : DF

    auto tr_0_0_yy_xx_xxx = pbuffer.data(idx_op_geom_020_df + 180);

    auto tr_0_0_yy_xx_xxy = pbuffer.data(idx_op_geom_020_df + 181);

    auto tr_0_0_yy_xx_xxz = pbuffer.data(idx_op_geom_020_df + 182);

    auto tr_0_0_yy_xx_xyy = pbuffer.data(idx_op_geom_020_df + 183);

    auto tr_0_0_yy_xx_xyz = pbuffer.data(idx_op_geom_020_df + 184);

    auto tr_0_0_yy_xx_xzz = pbuffer.data(idx_op_geom_020_df + 185);

    auto tr_0_0_yy_xx_yyy = pbuffer.data(idx_op_geom_020_df + 186);

    auto tr_0_0_yy_xx_yyz = pbuffer.data(idx_op_geom_020_df + 187);

    auto tr_0_0_yy_xx_yzz = pbuffer.data(idx_op_geom_020_df + 188);

    auto tr_0_0_yy_xx_zzz = pbuffer.data(idx_op_geom_020_df + 189);

    #pragma omp simd aligned(tr_0_0_yy_xx_xxx, tr_0_0_yy_xx_xxy, tr_0_0_yy_xx_xxz, tr_0_0_yy_xx_xyy, tr_0_0_yy_xx_xyz, tr_0_0_yy_xx_xzz, tr_0_0_yy_xx_yyy, tr_0_0_yy_xx_yyz, tr_0_0_yy_xx_yzz, tr_0_0_yy_xx_zzz, tr_xx_x, tr_xx_xxx, tr_xx_xxxyy, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxy_xx, tr_xxy_xxxy, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_yy, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xx_xxx[i] = -2.0 * tr_xx_xxx[i] * tbe_0 - 2.0 * tr_xx_xxx[i] * tke_0 + 4.0 * tr_xx_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xxy[i] = -2.0 * tr_xx_xxy[i] * tbe_0 - 6.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxy_xx[i] * tbe_0 + 8.0 * tr_xxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xxz[i] = -2.0 * tr_xx_xxz[i] * tbe_0 - 2.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xyy[i] = 2.0 * tr_xx_x[i] - 2.0 * tr_xx_xyy[i] * tbe_0 - 10.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xx_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxy_xy[i] * tbe_0 + 8.0 * tr_xxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xyz[i] = -2.0 * tr_xx_xyz[i] * tbe_0 - 6.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_xz[i] * tbe_0 + 8.0 * tr_xxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_xzz[i] = -2.0 * tr_xx_xzz[i] * tbe_0 - 2.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xx_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_yyy[i] = 6.0 * tr_xx_y[i] - 2.0 * tr_xx_yyy[i] * tbe_0 - 14.0 * tr_xx_yyy[i] * tke_0 + 4.0 * tr_xx_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxy_yy[i] * tbe_0 + 8.0 * tr_xxy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_yyz[i] = 2.0 * tr_xx_z[i] - 2.0 * tr_xx_yyz[i] * tbe_0 - 10.0 * tr_xx_yyz[i] * tke_0 + 4.0 * tr_xx_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxy_yz[i] * tbe_0 + 8.0 * tr_xxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_yzz[i] = -2.0 * tr_xx_yzz[i] * tbe_0 - 6.0 * tr_xx_yzz[i] * tke_0 + 4.0 * tr_xx_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxy_zz[i] * tbe_0 + 8.0 * tr_xxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xx_zzz[i] = -2.0 * tr_xx_zzz[i] * tbe_0 - 2.0 * tr_xx_zzz[i] * tke_0 + 4.0 * tr_xx_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 190-200 components of targeted buffer : DF

    auto tr_0_0_yy_xy_xxx = pbuffer.data(idx_op_geom_020_df + 190);

    auto tr_0_0_yy_xy_xxy = pbuffer.data(idx_op_geom_020_df + 191);

    auto tr_0_0_yy_xy_xxz = pbuffer.data(idx_op_geom_020_df + 192);

    auto tr_0_0_yy_xy_xyy = pbuffer.data(idx_op_geom_020_df + 193);

    auto tr_0_0_yy_xy_xyz = pbuffer.data(idx_op_geom_020_df + 194);

    auto tr_0_0_yy_xy_xzz = pbuffer.data(idx_op_geom_020_df + 195);

    auto tr_0_0_yy_xy_yyy = pbuffer.data(idx_op_geom_020_df + 196);

    auto tr_0_0_yy_xy_yyz = pbuffer.data(idx_op_geom_020_df + 197);

    auto tr_0_0_yy_xy_yzz = pbuffer.data(idx_op_geom_020_df + 198);

    auto tr_0_0_yy_xy_zzz = pbuffer.data(idx_op_geom_020_df + 199);

    #pragma omp simd aligned(tr_0_0_yy_xy_xxx, tr_0_0_yy_xy_xxy, tr_0_0_yy_xy_xxz, tr_0_0_yy_xy_xyy, tr_0_0_yy_xy_xyz, tr_0_0_yy_xy_xzz, tr_0_0_yy_xy_yyy, tr_0_0_yy_xy_yyz, tr_0_0_yy_xy_yzz, tr_0_0_yy_xy_zzz, tr_x_xx, tr_x_xxxy, tr_x_xxyy, tr_x_xxyz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_yy, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_xy_x, tr_xy_xxx, tr_xy_xxxyy, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxy, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_yy, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xy_xxx[i] = -4.0 * tr_x_xxxy[i] * tke_0 - 6.0 * tr_xy_xxx[i] * tbe_0 - 2.0 * tr_xy_xxx[i] * tke_0 + 4.0 * tr_xy_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xxy[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_x_xxyy[i] * tke_0 - 6.0 * tr_xy_xxy[i] * tbe_0 - 6.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyy_xx[i] * tbe_0 + 8.0 * tr_xyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xxz[i] = -4.0 * tr_x_xxyz[i] * tke_0 - 6.0 * tr_xy_xxz[i] * tbe_0 - 2.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xyy[i] = 4.0 * tr_x_xy[i] - 4.0 * tr_x_xyyy[i] * tke_0 + 2.0 * tr_xy_x[i] - 6.0 * tr_xy_xyy[i] * tbe_0 - 10.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xy_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyy_xy[i] * tbe_0 + 8.0 * tr_xyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xyz[i] = 2.0 * tr_x_xz[i] - 4.0 * tr_x_xyyz[i] * tke_0 - 6.0 * tr_xy_xyz[i] * tbe_0 - 6.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_xz[i] * tbe_0 + 8.0 * tr_xyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_xzz[i] = -4.0 * tr_x_xyzz[i] * tke_0 - 6.0 * tr_xy_xzz[i] * tbe_0 - 2.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_yyy[i] = 6.0 * tr_x_yy[i] - 4.0 * tr_x_yyyy[i] * tke_0 + 6.0 * tr_xy_y[i] - 6.0 * tr_xy_yyy[i] * tbe_0 - 14.0 * tr_xy_yyy[i] * tke_0 + 4.0 * tr_xy_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyy_yy[i] * tbe_0 + 8.0 * tr_xyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_yyz[i] = 4.0 * tr_x_yz[i] - 4.0 * tr_x_yyyz[i] * tke_0 + 2.0 * tr_xy_z[i] - 6.0 * tr_xy_yyz[i] * tbe_0 - 10.0 * tr_xy_yyz[i] * tke_0 + 4.0 * tr_xy_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyy_yz[i] * tbe_0 + 8.0 * tr_xyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_yzz[i] = 2.0 * tr_x_zz[i] - 4.0 * tr_x_yyzz[i] * tke_0 - 6.0 * tr_xy_yzz[i] * tbe_0 - 6.0 * tr_xy_yzz[i] * tke_0 + 4.0 * tr_xy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyy_zz[i] * tbe_0 + 8.0 * tr_xyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xy_zzz[i] = -4.0 * tr_x_yzzz[i] * tke_0 - 6.0 * tr_xy_zzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tke_0 + 4.0 * tr_xy_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 200-210 components of targeted buffer : DF

    auto tr_0_0_yy_xz_xxx = pbuffer.data(idx_op_geom_020_df + 200);

    auto tr_0_0_yy_xz_xxy = pbuffer.data(idx_op_geom_020_df + 201);

    auto tr_0_0_yy_xz_xxz = pbuffer.data(idx_op_geom_020_df + 202);

    auto tr_0_0_yy_xz_xyy = pbuffer.data(idx_op_geom_020_df + 203);

    auto tr_0_0_yy_xz_xyz = pbuffer.data(idx_op_geom_020_df + 204);

    auto tr_0_0_yy_xz_xzz = pbuffer.data(idx_op_geom_020_df + 205);

    auto tr_0_0_yy_xz_yyy = pbuffer.data(idx_op_geom_020_df + 206);

    auto tr_0_0_yy_xz_yyz = pbuffer.data(idx_op_geom_020_df + 207);

    auto tr_0_0_yy_xz_yzz = pbuffer.data(idx_op_geom_020_df + 208);

    auto tr_0_0_yy_xz_zzz = pbuffer.data(idx_op_geom_020_df + 209);

    #pragma omp simd aligned(tr_0_0_yy_xz_xxx, tr_0_0_yy_xz_xxy, tr_0_0_yy_xz_xxz, tr_0_0_yy_xz_xyy, tr_0_0_yy_xz_xyz, tr_0_0_yy_xz_xzz, tr_0_0_yy_xz_yyy, tr_0_0_yy_xz_yyz, tr_0_0_yy_xz_yzz, tr_0_0_yy_xz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xz_x, tr_xz_xxx, tr_xz_xxxyy, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xz_xxx[i] = -2.0 * tr_xz_xxx[i] * tbe_0 - 2.0 * tr_xz_xxx[i] * tke_0 + 4.0 * tr_xz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xxy[i] = -2.0 * tr_xz_xxy[i] * tbe_0 - 6.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 8.0 * tr_xyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xxz[i] = -2.0 * tr_xz_xxz[i] * tbe_0 - 2.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xyy[i] = 2.0 * tr_xz_x[i] - 2.0 * tr_xz_xyy[i] * tbe_0 - 10.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyz_xy[i] * tbe_0 + 8.0 * tr_xyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xyz[i] = -2.0 * tr_xz_xyz[i] * tbe_0 - 6.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 8.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_xzz[i] = -2.0 * tr_xz_xzz[i] * tbe_0 - 2.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_yyy[i] = 6.0 * tr_xz_y[i] - 2.0 * tr_xz_yyy[i] * tbe_0 - 14.0 * tr_xz_yyy[i] * tke_0 + 4.0 * tr_xz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyz_yy[i] * tbe_0 + 8.0 * tr_xyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_yyz[i] = 2.0 * tr_xz_z[i] - 2.0 * tr_xz_yyz[i] * tbe_0 - 10.0 * tr_xz_yyz[i] * tke_0 + 4.0 * tr_xz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyz_yz[i] * tbe_0 + 8.0 * tr_xyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_yzz[i] = -2.0 * tr_xz_yzz[i] * tbe_0 - 6.0 * tr_xz_yzz[i] * tke_0 + 4.0 * tr_xz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 + 8.0 * tr_xyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xz_zzz[i] = -2.0 * tr_xz_zzz[i] * tbe_0 - 2.0 * tr_xz_zzz[i] * tke_0 + 4.0 * tr_xz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-220 components of targeted buffer : DF

    auto tr_0_0_yy_yy_xxx = pbuffer.data(idx_op_geom_020_df + 210);

    auto tr_0_0_yy_yy_xxy = pbuffer.data(idx_op_geom_020_df + 211);

    auto tr_0_0_yy_yy_xxz = pbuffer.data(idx_op_geom_020_df + 212);

    auto tr_0_0_yy_yy_xyy = pbuffer.data(idx_op_geom_020_df + 213);

    auto tr_0_0_yy_yy_xyz = pbuffer.data(idx_op_geom_020_df + 214);

    auto tr_0_0_yy_yy_xzz = pbuffer.data(idx_op_geom_020_df + 215);

    auto tr_0_0_yy_yy_yyy = pbuffer.data(idx_op_geom_020_df + 216);

    auto tr_0_0_yy_yy_yyz = pbuffer.data(idx_op_geom_020_df + 217);

    auto tr_0_0_yy_yy_yzz = pbuffer.data(idx_op_geom_020_df + 218);

    auto tr_0_0_yy_yy_zzz = pbuffer.data(idx_op_geom_020_df + 219);

    #pragma omp simd aligned(tr_0_0_yy_yy_xxx, tr_0_0_yy_yy_xxy, tr_0_0_yy_yy_xxz, tr_0_0_yy_yy_xyy, tr_0_0_yy_yy_xyz, tr_0_0_yy_yy_xzz, tr_0_0_yy_yy_yyy, tr_0_0_yy_yy_yyz, tr_0_0_yy_yy_yzz, tr_0_0_yy_yy_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_y_xx, tr_y_xxxy, tr_y_xxyy, tr_y_xxyz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_yy, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxxyy, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, tr_yyy_xx, tr_yyy_xxxy, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_yy, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yy_xxx[i] = 2.0 * tr_0_xxx[i] - 8.0 * tr_y_xxxy[i] * tke_0 - 10.0 * tr_yy_xxx[i] * tbe_0 - 2.0 * tr_yy_xxx[i] * tke_0 + 4.0 * tr_yy_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xxy[i] = 2.0 * tr_0_xxy[i] + 4.0 * tr_y_xx[i] - 8.0 * tr_y_xxyy[i] * tke_0 - 10.0 * tr_yy_xxy[i] * tbe_0 - 6.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyy_xx[i] * tbe_0 + 8.0 * tr_yyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xxz[i] = 2.0 * tr_0_xxz[i] - 8.0 * tr_y_xxyz[i] * tke_0 - 10.0 * tr_yy_xxz[i] * tbe_0 - 2.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xyy[i] = 2.0 * tr_0_xyy[i] + 8.0 * tr_y_xy[i] - 8.0 * tr_y_xyyy[i] * tke_0 + 2.0 * tr_yy_x[i] - 10.0 * tr_yy_xyy[i] * tbe_0 - 10.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yy_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyy_xy[i] * tbe_0 + 8.0 * tr_yyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xyz[i] = 2.0 * tr_0_xyz[i] + 4.0 * tr_y_xz[i] - 8.0 * tr_y_xyyz[i] * tke_0 - 10.0 * tr_yy_xyz[i] * tbe_0 - 6.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyy_xz[i] * tbe_0 + 8.0 * tr_yyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_xzz[i] = 2.0 * tr_0_xzz[i] - 8.0 * tr_y_xyzz[i] * tke_0 - 10.0 * tr_yy_xzz[i] * tbe_0 - 2.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_yyy[i] = 2.0 * tr_0_yyy[i] + 12.0 * tr_y_yy[i] - 8.0 * tr_y_yyyy[i] * tke_0 + 6.0 * tr_yy_y[i] - 10.0 * tr_yy_yyy[i] * tbe_0 - 14.0 * tr_yy_yyy[i] * tke_0 + 4.0 * tr_yy_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyy_yy[i] * tbe_0 + 8.0 * tr_yyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_yyz[i] = 2.0 * tr_0_yyz[i] + 8.0 * tr_y_yz[i] - 8.0 * tr_y_yyyz[i] * tke_0 + 2.0 * tr_yy_z[i] - 10.0 * tr_yy_yyz[i] * tbe_0 - 10.0 * tr_yy_yyz[i] * tke_0 + 4.0 * tr_yy_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyy_yz[i] * tbe_0 + 8.0 * tr_yyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_yzz[i] = 2.0 * tr_0_yzz[i] + 4.0 * tr_y_zz[i] - 8.0 * tr_y_yyzz[i] * tke_0 - 10.0 * tr_yy_yzz[i] * tbe_0 - 6.0 * tr_yy_yzz[i] * tke_0 + 4.0 * tr_yy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyy_zz[i] * tbe_0 + 8.0 * tr_yyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yy_zzz[i] = 2.0 * tr_0_zzz[i] - 8.0 * tr_y_yzzz[i] * tke_0 - 10.0 * tr_yy_zzz[i] * tbe_0 - 2.0 * tr_yy_zzz[i] * tke_0 + 4.0 * tr_yy_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 220-230 components of targeted buffer : DF

    auto tr_0_0_yy_yz_xxx = pbuffer.data(idx_op_geom_020_df + 220);

    auto tr_0_0_yy_yz_xxy = pbuffer.data(idx_op_geom_020_df + 221);

    auto tr_0_0_yy_yz_xxz = pbuffer.data(idx_op_geom_020_df + 222);

    auto tr_0_0_yy_yz_xyy = pbuffer.data(idx_op_geom_020_df + 223);

    auto tr_0_0_yy_yz_xyz = pbuffer.data(idx_op_geom_020_df + 224);

    auto tr_0_0_yy_yz_xzz = pbuffer.data(idx_op_geom_020_df + 225);

    auto tr_0_0_yy_yz_yyy = pbuffer.data(idx_op_geom_020_df + 226);

    auto tr_0_0_yy_yz_yyz = pbuffer.data(idx_op_geom_020_df + 227);

    auto tr_0_0_yy_yz_yzz = pbuffer.data(idx_op_geom_020_df + 228);

    auto tr_0_0_yy_yz_zzz = pbuffer.data(idx_op_geom_020_df + 229);

    #pragma omp simd aligned(tr_0_0_yy_yz_xxx, tr_0_0_yy_yz_xxy, tr_0_0_yy_yz_xxz, tr_0_0_yy_yz_xyy, tr_0_0_yy_yz_xyz, tr_0_0_yy_yz_xzz, tr_0_0_yy_yz_yyy, tr_0_0_yy_yz_yyz, tr_0_0_yy_yz_yzz, tr_0_0_yy_yz_zzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyz_xx, tr_yyz_xxxy, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_yy, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxxyy, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, tr_z_xx, tr_z_xxxy, tr_z_xxyy, tr_z_xxyz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_yy, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yz_xxx[i] = -4.0 * tr_z_xxxy[i] * tke_0 - 6.0 * tr_yz_xxx[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tke_0 + 4.0 * tr_yz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xxy[i] = 2.0 * tr_z_xx[i] - 4.0 * tr_z_xxyy[i] * tke_0 - 6.0 * tr_yz_xxy[i] * tbe_0 - 6.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xx[i] * tbe_0 + 8.0 * tr_yyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xxz[i] = -4.0 * tr_z_xxyz[i] * tke_0 - 6.0 * tr_yz_xxz[i] * tbe_0 - 2.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xyy[i] = 4.0 * tr_z_xy[i] - 4.0 * tr_z_xyyy[i] * tke_0 + 2.0 * tr_yz_x[i] - 6.0 * tr_yz_xyy[i] * tbe_0 - 10.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyz_xy[i] * tbe_0 + 8.0 * tr_yyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xyz[i] = 2.0 * tr_z_xz[i] - 4.0 * tr_z_xyyz[i] * tke_0 - 6.0 * tr_yz_xyz[i] * tbe_0 - 6.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xz[i] * tbe_0 + 8.0 * tr_yyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_xzz[i] = -4.0 * tr_z_xyzz[i] * tke_0 - 6.0 * tr_yz_xzz[i] * tbe_0 - 2.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_yyy[i] = 6.0 * tr_z_yy[i] - 4.0 * tr_z_yyyy[i] * tke_0 + 6.0 * tr_yz_y[i] - 6.0 * tr_yz_yyy[i] * tbe_0 - 14.0 * tr_yz_yyy[i] * tke_0 + 4.0 * tr_yz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyz_yy[i] * tbe_0 + 8.0 * tr_yyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_yyz[i] = 4.0 * tr_z_yz[i] - 4.0 * tr_z_yyyz[i] * tke_0 + 2.0 * tr_yz_z[i] - 6.0 * tr_yz_yyz[i] * tbe_0 - 10.0 * tr_yz_yyz[i] * tke_0 + 4.0 * tr_yz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyz_yz[i] * tbe_0 + 8.0 * tr_yyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_yzz[i] = 2.0 * tr_z_zz[i] - 4.0 * tr_z_yyzz[i] * tke_0 - 6.0 * tr_yz_yzz[i] * tbe_0 - 6.0 * tr_yz_yzz[i] * tke_0 + 4.0 * tr_yz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_zz[i] * tbe_0 + 8.0 * tr_yyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yz_zzz[i] = -4.0 * tr_z_yzzz[i] * tke_0 - 6.0 * tr_yz_zzz[i] * tbe_0 - 2.0 * tr_yz_zzz[i] * tke_0 + 4.0 * tr_yz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 230-240 components of targeted buffer : DF

    auto tr_0_0_yy_zz_xxx = pbuffer.data(idx_op_geom_020_df + 230);

    auto tr_0_0_yy_zz_xxy = pbuffer.data(idx_op_geom_020_df + 231);

    auto tr_0_0_yy_zz_xxz = pbuffer.data(idx_op_geom_020_df + 232);

    auto tr_0_0_yy_zz_xyy = pbuffer.data(idx_op_geom_020_df + 233);

    auto tr_0_0_yy_zz_xyz = pbuffer.data(idx_op_geom_020_df + 234);

    auto tr_0_0_yy_zz_xzz = pbuffer.data(idx_op_geom_020_df + 235);

    auto tr_0_0_yy_zz_yyy = pbuffer.data(idx_op_geom_020_df + 236);

    auto tr_0_0_yy_zz_yyz = pbuffer.data(idx_op_geom_020_df + 237);

    auto tr_0_0_yy_zz_yzz = pbuffer.data(idx_op_geom_020_df + 238);

    auto tr_0_0_yy_zz_zzz = pbuffer.data(idx_op_geom_020_df + 239);

    #pragma omp simd aligned(tr_0_0_yy_zz_xxx, tr_0_0_yy_zz_xxy, tr_0_0_yy_zz_xxz, tr_0_0_yy_zz_xyy, tr_0_0_yy_zz_xyz, tr_0_0_yy_zz_xzz, tr_0_0_yy_zz_yyy, tr_0_0_yy_zz_yyz, tr_0_0_yy_zz_yzz, tr_0_0_yy_zz_zzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yzz_xx, tr_yzz_xxxy, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_yy, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxxyy, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_zz_xxx[i] = -2.0 * tr_zz_xxx[i] * tbe_0 - 2.0 * tr_zz_xxx[i] * tke_0 + 4.0 * tr_zz_xxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xxy[i] = -2.0 * tr_zz_xxy[i] * tbe_0 - 6.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xx[i] * tbe_0 + 8.0 * tr_yzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xxz[i] = -2.0 * tr_zz_xxz[i] * tbe_0 - 2.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xyy[i] = 2.0 * tr_zz_x[i] - 2.0 * tr_zz_xyy[i] * tbe_0 - 10.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zz_xyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yzz_xy[i] * tbe_0 + 8.0 * tr_yzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xyz[i] = -2.0 * tr_zz_xyz[i] * tbe_0 - 6.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xz[i] * tbe_0 + 8.0 * tr_yzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_xzz[i] = -2.0 * tr_zz_xzz[i] * tbe_0 - 2.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_zz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_yyy[i] = 6.0 * tr_zz_y[i] - 2.0 * tr_zz_yyy[i] * tbe_0 - 14.0 * tr_zz_yyy[i] * tke_0 + 4.0 * tr_zz_yyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yzz_yy[i] * tbe_0 + 8.0 * tr_yzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_yyz[i] = 2.0 * tr_zz_z[i] - 2.0 * tr_zz_yyz[i] * tbe_0 - 10.0 * tr_zz_yyz[i] * tke_0 + 4.0 * tr_zz_yyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yzz_yz[i] * tbe_0 + 8.0 * tr_yzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_yzz[i] = -2.0 * tr_zz_yzz[i] * tbe_0 - 6.0 * tr_zz_yzz[i] * tke_0 + 4.0 * tr_zz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_zz[i] * tbe_0 + 8.0 * tr_yzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zz_zzz[i] = -2.0 * tr_zz_zzz[i] * tbe_0 - 2.0 * tr_zz_zzz[i] * tke_0 + 4.0 * tr_zz_yyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-250 components of targeted buffer : DF

    auto tr_0_0_yz_xx_xxx = pbuffer.data(idx_op_geom_020_df + 240);

    auto tr_0_0_yz_xx_xxy = pbuffer.data(idx_op_geom_020_df + 241);

    auto tr_0_0_yz_xx_xxz = pbuffer.data(idx_op_geom_020_df + 242);

    auto tr_0_0_yz_xx_xyy = pbuffer.data(idx_op_geom_020_df + 243);

    auto tr_0_0_yz_xx_xyz = pbuffer.data(idx_op_geom_020_df + 244);

    auto tr_0_0_yz_xx_xzz = pbuffer.data(idx_op_geom_020_df + 245);

    auto tr_0_0_yz_xx_yyy = pbuffer.data(idx_op_geom_020_df + 246);

    auto tr_0_0_yz_xx_yyz = pbuffer.data(idx_op_geom_020_df + 247);

    auto tr_0_0_yz_xx_yzz = pbuffer.data(idx_op_geom_020_df + 248);

    auto tr_0_0_yz_xx_zzz = pbuffer.data(idx_op_geom_020_df + 249);

    #pragma omp simd aligned(tr_0_0_yz_xx_xxx, tr_0_0_yz_xx_xxy, tr_0_0_yz_xx_xxz, tr_0_0_yz_xx_xyy, tr_0_0_yz_xx_xyz, tr_0_0_yz_xx_xzz, tr_0_0_yz_xx_yyy, tr_0_0_yz_xx_yyz, tr_0_0_yz_xx_yzz, tr_0_0_yz_xx_zzz, tr_xx_x, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_y, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_z, tr_xx_zzz, tr_xxy_xx, tr_xxy_xxxz, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxz_xx, tr_xxz_xxxy, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_yy, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xx_xxx[i] = 4.0 * tr_xx_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xxy[i] = -2.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xxz[i] = -2.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xyy[i] = -4.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xyz[i] = tr_xx_x[i] - 2.0 * tr_xx_xzz[i] * tke_0 - 2.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xx_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_xzz[i] = -4.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_yyy[i] = -6.0 * tr_xx_yyz[i] * tke_0 + 4.0 * tr_xx_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_yyz[i] = 2.0 * tr_xx_y[i] - 4.0 * tr_xx_yzz[i] * tke_0 - 2.0 * tr_xx_yyy[i] * tke_0 + 4.0 * tr_xx_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_yzz[i] = 2.0 * tr_xx_z[i] - 2.0 * tr_xx_zzz[i] * tke_0 - 4.0 * tr_xx_yyz[i] * tke_0 + 4.0 * tr_xx_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xx_zzz[i] = -6.0 * tr_xx_yzz[i] * tke_0 + 4.0 * tr_xx_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 250-260 components of targeted buffer : DF

    auto tr_0_0_yz_xy_xxx = pbuffer.data(idx_op_geom_020_df + 250);

    auto tr_0_0_yz_xy_xxy = pbuffer.data(idx_op_geom_020_df + 251);

    auto tr_0_0_yz_xy_xxz = pbuffer.data(idx_op_geom_020_df + 252);

    auto tr_0_0_yz_xy_xyy = pbuffer.data(idx_op_geom_020_df + 253);

    auto tr_0_0_yz_xy_xyz = pbuffer.data(idx_op_geom_020_df + 254);

    auto tr_0_0_yz_xy_xzz = pbuffer.data(idx_op_geom_020_df + 255);

    auto tr_0_0_yz_xy_yyy = pbuffer.data(idx_op_geom_020_df + 256);

    auto tr_0_0_yz_xy_yyz = pbuffer.data(idx_op_geom_020_df + 257);

    auto tr_0_0_yz_xy_yzz = pbuffer.data(idx_op_geom_020_df + 258);

    auto tr_0_0_yz_xy_zzz = pbuffer.data(idx_op_geom_020_df + 259);

    #pragma omp simd aligned(tr_0_0_yz_xy_xxx, tr_0_0_yz_xy_xxy, tr_0_0_yz_xy_xxz, tr_0_0_yz_xy_xyy, tr_0_0_yz_xy_xyz, tr_0_0_yz_xy_xzz, tr_0_0_yz_xy_yyy, tr_0_0_yz_xy_yyz, tr_0_0_yz_xy_yzz, tr_0_0_yz_xy_zzz, tr_x_xx, tr_x_xxxz, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xy_x, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_y, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_z, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxz, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xy_xxx[i] = -2.0 * tr_x_xxxz[i] * tke_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xxy[i] = -2.0 * tr_x_xxyz[i] * tke_0 - 2.0 * tr_xz_xxy[i] * tbe_0 - 2.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xxz[i] = tr_x_xx[i] - 2.0 * tr_x_xxzz[i] * tke_0 - 2.0 * tr_xz_xxz[i] * tbe_0 - 2.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xyy[i] = -2.0 * tr_x_xyyz[i] * tke_0 - 2.0 * tr_xz_xyy[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xyz[i] = tr_x_xy[i] - 2.0 * tr_x_xyzz[i] * tke_0 - 2.0 * tr_xz_xyz[i] * tbe_0 + tr_xy_x[i] - 2.0 * tr_xy_xzz[i] * tke_0 - 2.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xy_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_xzz[i] = 2.0 * tr_x_xz[i] - 2.0 * tr_x_xzzz[i] * tke_0 - 2.0 * tr_xz_xzz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_yyy[i] = -2.0 * tr_x_yyyz[i] * tke_0 - 2.0 * tr_xz_yyy[i] * tbe_0 - 6.0 * tr_xy_yyz[i] * tke_0 + 4.0 * tr_xy_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_yyz[i] = tr_x_yy[i] - 2.0 * tr_x_yyzz[i] * tke_0 - 2.0 * tr_xz_yyz[i] * tbe_0 + 2.0 * tr_xy_y[i] - 4.0 * tr_xy_yzz[i] * tke_0 - 2.0 * tr_xy_yyy[i] * tke_0 + 4.0 * tr_xy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_yzz[i] = 2.0 * tr_x_yz[i] - 2.0 * tr_x_yzzz[i] * tke_0 - 2.0 * tr_xz_yzz[i] * tbe_0 + 2.0 * tr_xy_z[i] - 2.0 * tr_xy_zzz[i] * tke_0 - 4.0 * tr_xy_yyz[i] * tke_0 + 4.0 * tr_xy_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xy_zzz[i] = 3.0 * tr_x_zz[i] - 2.0 * tr_x_zzzz[i] * tke_0 - 2.0 * tr_xz_zzz[i] * tbe_0 - 6.0 * tr_xy_yzz[i] * tke_0 + 4.0 * tr_xy_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 260-270 components of targeted buffer : DF

    auto tr_0_0_yz_xz_xxx = pbuffer.data(idx_op_geom_020_df + 260);

    auto tr_0_0_yz_xz_xxy = pbuffer.data(idx_op_geom_020_df + 261);

    auto tr_0_0_yz_xz_xxz = pbuffer.data(idx_op_geom_020_df + 262);

    auto tr_0_0_yz_xz_xyy = pbuffer.data(idx_op_geom_020_df + 263);

    auto tr_0_0_yz_xz_xyz = pbuffer.data(idx_op_geom_020_df + 264);

    auto tr_0_0_yz_xz_xzz = pbuffer.data(idx_op_geom_020_df + 265);

    auto tr_0_0_yz_xz_yyy = pbuffer.data(idx_op_geom_020_df + 266);

    auto tr_0_0_yz_xz_yyz = pbuffer.data(idx_op_geom_020_df + 267);

    auto tr_0_0_yz_xz_yzz = pbuffer.data(idx_op_geom_020_df + 268);

    auto tr_0_0_yz_xz_zzz = pbuffer.data(idx_op_geom_020_df + 269);

    #pragma omp simd aligned(tr_0_0_yz_xz_xxx, tr_0_0_yz_xz_xxy, tr_0_0_yz_xz_xxz, tr_0_0_yz_xz_xyy, tr_0_0_yz_xz_xyz, tr_0_0_yz_xz_xzz, tr_0_0_yz_xz_yyy, tr_0_0_yz_xz_yyz, tr_0_0_yz_xz_yzz, tr_0_0_yz_xz_zzz, tr_x_xx, tr_x_xxxy, tr_x_xxyy, tr_x_xxyz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_yy, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xz_x, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_y, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_z, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxy, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_yy, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xz_xxx[i] = -2.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_xz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xxy[i] = tr_x_xx[i] - 2.0 * tr_x_xxyy[i] * tke_0 - 2.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xxz[i] = -2.0 * tr_x_xxyz[i] * tke_0 - 2.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxz[i] * tbe_0 - 2.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xyy[i] = 2.0 * tr_x_xy[i] - 2.0 * tr_x_xyyy[i] * tke_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xyz[i] = tr_x_xz[i] - 2.0 * tr_x_xyyz[i] * tke_0 + tr_xz_x[i] - 2.0 * tr_xz_xzz[i] * tke_0 - 2.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyz[i] * tbe_0 - 2.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_xzz[i] = -2.0 * tr_x_xyzz[i] * tke_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xzz[i] * tbe_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_yyy[i] = 3.0 * tr_x_yy[i] - 2.0 * tr_x_yyyy[i] * tke_0 - 6.0 * tr_xz_yyz[i] * tke_0 + 4.0 * tr_xz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_yyz[i] = 2.0 * tr_x_yz[i] - 2.0 * tr_x_yyyz[i] * tke_0 + 2.0 * tr_xz_y[i] - 4.0 * tr_xz_yzz[i] * tke_0 - 2.0 * tr_xz_yyy[i] * tke_0 + 4.0 * tr_xz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyz[i] * tbe_0 - 2.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_yzz[i] = tr_x_zz[i] - 2.0 * tr_x_yyzz[i] * tke_0 + 2.0 * tr_xz_z[i] - 2.0 * tr_xz_zzz[i] * tke_0 - 4.0 * tr_xz_yyz[i] * tke_0 + 4.0 * tr_xz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yzz[i] * tbe_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xz_zzz[i] = -2.0 * tr_x_yzzz[i] * tke_0 - 6.0 * tr_xz_yzz[i] * tke_0 + 4.0 * tr_xz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_zzz[i] * tbe_0 - 6.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-280 components of targeted buffer : DF

    auto tr_0_0_yz_yy_xxx = pbuffer.data(idx_op_geom_020_df + 270);

    auto tr_0_0_yz_yy_xxy = pbuffer.data(idx_op_geom_020_df + 271);

    auto tr_0_0_yz_yy_xxz = pbuffer.data(idx_op_geom_020_df + 272);

    auto tr_0_0_yz_yy_xyy = pbuffer.data(idx_op_geom_020_df + 273);

    auto tr_0_0_yz_yy_xyz = pbuffer.data(idx_op_geom_020_df + 274);

    auto tr_0_0_yz_yy_xzz = pbuffer.data(idx_op_geom_020_df + 275);

    auto tr_0_0_yz_yy_yyy = pbuffer.data(idx_op_geom_020_df + 276);

    auto tr_0_0_yz_yy_yyz = pbuffer.data(idx_op_geom_020_df + 277);

    auto tr_0_0_yz_yy_yzz = pbuffer.data(idx_op_geom_020_df + 278);

    auto tr_0_0_yz_yy_zzz = pbuffer.data(idx_op_geom_020_df + 279);

    #pragma omp simd aligned(tr_0_0_yz_yy_xxx, tr_0_0_yz_yy_xxy, tr_0_0_yz_yy_xxz, tr_0_0_yz_yy_xyy, tr_0_0_yz_yy_xyz, tr_0_0_yz_yy_xzz, tr_0_0_yz_yy_yyy, tr_0_0_yz_yy_yyz, tr_0_0_yz_yy_yzz, tr_0_0_yz_yy_zzz, tr_y_xx, tr_y_xxxz, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, tr_yy_x, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_y, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_z, tr_yy_zzz, tr_yyy_xx, tr_yyy_xxxz, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyz_xx, tr_yyz_xxxy, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_yy, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yy_xxx[i] = -4.0 * tr_y_xxxz[i] * tke_0 - 4.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xxy[i] = -4.0 * tr_y_xxyz[i] * tke_0 - 4.0 * tr_yz_xxy[i] * tbe_0 - 2.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xxz[i] = 2.0 * tr_y_xx[i] - 4.0 * tr_y_xxzz[i] * tke_0 - 4.0 * tr_yz_xxz[i] * tbe_0 - 2.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_yyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xyy[i] = -4.0 * tr_y_xyyz[i] * tke_0 - 4.0 * tr_yz_xyy[i] * tbe_0 - 4.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xyz[i] = 2.0 * tr_y_xy[i] - 4.0 * tr_y_xyzz[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + tr_yy_x[i] - 2.0 * tr_yy_xzz[i] * tke_0 - 2.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yy_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xy[i] * tbe_0 + 4.0 * tr_yyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_xzz[i] = 4.0 * tr_y_xz[i] - 4.0 * tr_y_xzzz[i] * tke_0 - 4.0 * tr_yz_xzz[i] * tbe_0 - 4.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyy_xz[i] * tbe_0 + 4.0 * tr_yyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_yyy[i] = -4.0 * tr_y_yyyz[i] * tke_0 - 4.0 * tr_yz_yyy[i] * tbe_0 - 6.0 * tr_yy_yyz[i] * tke_0 + 4.0 * tr_yy_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_yyz[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_y_yyzz[i] * tke_0 - 4.0 * tr_yz_yyz[i] * tbe_0 + 2.0 * tr_yy_y[i] - 4.0 * tr_yy_yzz[i] * tke_0 - 2.0 * tr_yy_yyy[i] * tke_0 + 4.0 * tr_yy_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yy[i] * tbe_0 + 4.0 * tr_yyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_yzz[i] = 4.0 * tr_y_yz[i] - 4.0 * tr_y_yzzz[i] * tke_0 - 4.0 * tr_yz_yzz[i] * tbe_0 + 2.0 * tr_yy_z[i] - 2.0 * tr_yy_zzz[i] * tke_0 - 4.0 * tr_yy_yyz[i] * tke_0 + 4.0 * tr_yy_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyy_yz[i] * tbe_0 + 4.0 * tr_yyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yy_zzz[i] = 6.0 * tr_y_zz[i] - 4.0 * tr_y_zzzz[i] * tke_0 - 4.0 * tr_yz_zzz[i] * tbe_0 - 6.0 * tr_yy_yzz[i] * tke_0 + 4.0 * tr_yy_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyy_zz[i] * tbe_0 + 4.0 * tr_yyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 280-290 components of targeted buffer : DF

    auto tr_0_0_yz_yz_xxx = pbuffer.data(idx_op_geom_020_df + 280);

    auto tr_0_0_yz_yz_xxy = pbuffer.data(idx_op_geom_020_df + 281);

    auto tr_0_0_yz_yz_xxz = pbuffer.data(idx_op_geom_020_df + 282);

    auto tr_0_0_yz_yz_xyy = pbuffer.data(idx_op_geom_020_df + 283);

    auto tr_0_0_yz_yz_xyz = pbuffer.data(idx_op_geom_020_df + 284);

    auto tr_0_0_yz_yz_xzz = pbuffer.data(idx_op_geom_020_df + 285);

    auto tr_0_0_yz_yz_yyy = pbuffer.data(idx_op_geom_020_df + 286);

    auto tr_0_0_yz_yz_yyz = pbuffer.data(idx_op_geom_020_df + 287);

    auto tr_0_0_yz_yz_yzz = pbuffer.data(idx_op_geom_020_df + 288);

    auto tr_0_0_yz_yz_zzz = pbuffer.data(idx_op_geom_020_df + 289);

    #pragma omp simd aligned(tr_0_0_yz_yz_xxx, tr_0_0_yz_yz_xxy, tr_0_0_yz_yz_xxz, tr_0_0_yz_yz_xyy, tr_0_0_yz_yz_xyz, tr_0_0_yz_yz_xzz, tr_0_0_yz_yz_yyy, tr_0_0_yz_yz_yyz, tr_0_0_yz_yz_yzz, tr_0_0_yz_yz_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_y_xx, tr_y_xxxy, tr_y_xxyy, tr_y_xxyz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_yy, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyz_xx, tr_yyz_xxxz, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yz_x, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_y, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_z, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxy, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_yy, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_z_xx, tr_z_xxxz, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yz_xxx[i] = tr_0_xxx[i] - 2.0 * tr_z_xxxz[i] * tke_0 - 2.0 * tr_zz_xxx[i] * tbe_0 - 2.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_yz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxx[i] * tbe_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xxy[i] = tr_0_xxy[i] - 2.0 * tr_z_xxyz[i] * tke_0 - 2.0 * tr_zz_xxy[i] * tbe_0 + tr_y_xx[i] - 2.0 * tr_y_xxyy[i] * tke_0 - 2.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxy[i] * tbe_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xxz[i] = tr_0_xxz[i] + tr_z_xx[i] - 2.0 * tr_z_xxzz[i] * tke_0 - 2.0 * tr_zz_xxz[i] * tbe_0 - 2.0 * tr_y_xxyz[i] * tke_0 - 2.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxz[i] * tbe_0 - 2.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xyy[i] = tr_0_xyy[i] - 2.0 * tr_z_xyyz[i] * tke_0 - 2.0 * tr_zz_xyy[i] * tbe_0 + 2.0 * tr_y_xy[i] - 2.0 * tr_y_xyyy[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xyy[i] * tbe_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xyz[i] = tr_0_xyz[i] + tr_z_xy[i] - 2.0 * tr_z_xyzz[i] * tke_0 - 2.0 * tr_zz_xyz[i] * tbe_0 + tr_y_xz[i] - 2.0 * tr_y_xyyz[i] * tke_0 + tr_yz_x[i] - 2.0 * tr_yz_xzz[i] * tke_0 - 2.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xyz[i] * tbe_0 - 2.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_xzz[i] = tr_0_xzz[i] + 2.0 * tr_z_xz[i] - 2.0 * tr_z_xzzz[i] * tke_0 - 2.0 * tr_zz_xzz[i] * tbe_0 - 2.0 * tr_y_xyzz[i] * tke_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xzz[i] * tbe_0 - 4.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_yyy[i] = tr_0_yyy[i] - 2.0 * tr_z_yyyz[i] * tke_0 - 2.0 * tr_zz_yyy[i] * tbe_0 + 3.0 * tr_y_yy[i] - 2.0 * tr_y_yyyy[i] * tke_0 - 6.0 * tr_yz_yyz[i] * tke_0 + 4.0 * tr_yz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yyy[i] * tbe_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_yyz[i] = tr_0_yyz[i] + tr_z_yy[i] - 2.0 * tr_z_yyzz[i] * tke_0 - 2.0 * tr_zz_yyz[i] * tbe_0 + 2.0 * tr_y_yz[i] - 2.0 * tr_y_yyyz[i] * tke_0 + 2.0 * tr_yz_y[i] - 4.0 * tr_yz_yzz[i] * tke_0 - 2.0 * tr_yz_yyy[i] * tke_0 + 4.0 * tr_yz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yyz[i] * tbe_0 - 2.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_yzz[i] = tr_0_yzz[i] + 2.0 * tr_z_yz[i] - 2.0 * tr_z_yzzz[i] * tke_0 - 2.0 * tr_zz_yzz[i] * tbe_0 + tr_y_zz[i] - 2.0 * tr_y_yyzz[i] * tke_0 + 2.0 * tr_yz_z[i] - 2.0 * tr_yz_zzz[i] * tke_0 - 4.0 * tr_yz_yyz[i] * tke_0 + 4.0 * tr_yz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yzz[i] * tbe_0 - 4.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yz_zzz[i] = tr_0_zzz[i] + 3.0 * tr_z_zz[i] - 2.0 * tr_z_zzzz[i] * tke_0 - 2.0 * tr_zz_zzz[i] * tbe_0 - 2.0 * tr_y_yzzz[i] * tke_0 - 6.0 * tr_yz_yzz[i] * tke_0 + 4.0 * tr_yz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_zzz[i] * tbe_0 - 6.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 290-300 components of targeted buffer : DF

    auto tr_0_0_yz_zz_xxx = pbuffer.data(idx_op_geom_020_df + 290);

    auto tr_0_0_yz_zz_xxy = pbuffer.data(idx_op_geom_020_df + 291);

    auto tr_0_0_yz_zz_xxz = pbuffer.data(idx_op_geom_020_df + 292);

    auto tr_0_0_yz_zz_xyy = pbuffer.data(idx_op_geom_020_df + 293);

    auto tr_0_0_yz_zz_xyz = pbuffer.data(idx_op_geom_020_df + 294);

    auto tr_0_0_yz_zz_xzz = pbuffer.data(idx_op_geom_020_df + 295);

    auto tr_0_0_yz_zz_yyy = pbuffer.data(idx_op_geom_020_df + 296);

    auto tr_0_0_yz_zz_yyz = pbuffer.data(idx_op_geom_020_df + 297);

    auto tr_0_0_yz_zz_yzz = pbuffer.data(idx_op_geom_020_df + 298);

    auto tr_0_0_yz_zz_zzz = pbuffer.data(idx_op_geom_020_df + 299);

    #pragma omp simd aligned(tr_0_0_yz_zz_xxx, tr_0_0_yz_zz_xxy, tr_0_0_yz_zz_xxz, tr_0_0_yz_zz_xyy, tr_0_0_yz_zz_xyz, tr_0_0_yz_zz_xzz, tr_0_0_yz_zz_yyy, tr_0_0_yz_zz_yyz, tr_0_0_yz_zz_yzz, tr_0_0_yz_zz_zzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxz, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_z_xx, tr_z_xxxy, tr_z_xxyy, tr_z_xxyz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_yy, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_zz_x, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_y, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_z, tr_zz_zzz, tr_zzz_xx, tr_zzz_xxxy, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_yy, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_zz_xxx[i] = -4.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_zz_xxxyz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xxy[i] = 2.0 * tr_z_xx[i] - 4.0 * tr_z_xxyy[i] * tke_0 - 2.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_xx[i] * tbe_0 + 4.0 * tr_zzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xxz[i] = -4.0 * tr_z_xxyz[i] * tke_0 - 2.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xxz[i] * tbe_0 - 2.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xyy[i] = 4.0 * tr_z_xy[i] - 4.0 * tr_z_xyyy[i] * tke_0 - 4.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xyyyz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_xy[i] * tbe_0 + 4.0 * tr_zzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xyz[i] = 2.0 * tr_z_xz[i] - 4.0 * tr_z_xyyz[i] * tke_0 + tr_zz_x[i] - 2.0 * tr_zz_xzz[i] * tke_0 - 2.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zz_xyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_xz[i] * tbe_0 + 4.0 * tr_zzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xyz[i] * tbe_0 - 2.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_xzz[i] = -4.0 * tr_z_xyzz[i] * tke_0 - 4.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_xzz[i] * tbe_0 - 4.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_yyy[i] = 6.0 * tr_z_yy[i] - 4.0 * tr_z_yyyy[i] * tke_0 - 6.0 * tr_zz_yyz[i] * tke_0 + 4.0 * tr_zz_yyyyz[i] * tke_0 * tke_0 - 6.0 * tr_zzz_yy[i] * tbe_0 + 4.0 * tr_zzz_yyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_yyz[i] = 4.0 * tr_z_yz[i] - 4.0 * tr_z_yyyz[i] * tke_0 + 2.0 * tr_zz_y[i] - 4.0 * tr_zz_yzz[i] * tke_0 - 2.0 * tr_zz_yyy[i] * tke_0 + 4.0 * tr_zz_yyyzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_yz[i] * tbe_0 + 4.0 * tr_zzz_yyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yyz[i] * tbe_0 - 2.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_yzz[i] = 2.0 * tr_z_zz[i] - 4.0 * tr_z_yyzz[i] * tke_0 + 2.0 * tr_zz_z[i] - 2.0 * tr_zz_zzz[i] * tke_0 - 4.0 * tr_zz_yyz[i] * tke_0 + 4.0 * tr_zz_yyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzz_zz[i] * tbe_0 + 4.0 * tr_zzz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_yzz[i] * tbe_0 - 4.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zz_zzz[i] = -4.0 * tr_z_yzzz[i] * tke_0 - 6.0 * tr_zz_yzz[i] * tke_0 + 4.0 * tr_zz_yzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzz_yzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yz_zzz[i] * tbe_0 - 6.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_yzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-310 components of targeted buffer : DF

    auto tr_0_0_zz_xx_xxx = pbuffer.data(idx_op_geom_020_df + 300);

    auto tr_0_0_zz_xx_xxy = pbuffer.data(idx_op_geom_020_df + 301);

    auto tr_0_0_zz_xx_xxz = pbuffer.data(idx_op_geom_020_df + 302);

    auto tr_0_0_zz_xx_xyy = pbuffer.data(idx_op_geom_020_df + 303);

    auto tr_0_0_zz_xx_xyz = pbuffer.data(idx_op_geom_020_df + 304);

    auto tr_0_0_zz_xx_xzz = pbuffer.data(idx_op_geom_020_df + 305);

    auto tr_0_0_zz_xx_yyy = pbuffer.data(idx_op_geom_020_df + 306);

    auto tr_0_0_zz_xx_yyz = pbuffer.data(idx_op_geom_020_df + 307);

    auto tr_0_0_zz_xx_yzz = pbuffer.data(idx_op_geom_020_df + 308);

    auto tr_0_0_zz_xx_zzz = pbuffer.data(idx_op_geom_020_df + 309);

    #pragma omp simd aligned(tr_0_0_zz_xx_xxx, tr_0_0_zz_xx_xxy, tr_0_0_zz_xx_xxz, tr_0_0_zz_xx_xyy, tr_0_0_zz_xx_xyz, tr_0_0_zz_xx_xzz, tr_0_0_zz_xx_yyy, tr_0_0_zz_xx_yyz, tr_0_0_zz_xx_yzz, tr_0_0_zz_xx_zzz, tr_xx_x, tr_xx_xxx, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_y, tr_xx_yyy, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_z, tr_xx_zzz, tr_xx_zzzzz, tr_xxz_xx, tr_xxz_xxxz, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xx_xxx[i] = -2.0 * tr_xx_xxx[i] * tbe_0 - 2.0 * tr_xx_xxx[i] * tke_0 + 4.0 * tr_xx_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xxy[i] = -2.0 * tr_xx_xxy[i] * tbe_0 - 2.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xx_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xxz[i] = -2.0 * tr_xx_xxz[i] * tbe_0 - 6.0 * tr_xx_xxz[i] * tke_0 + 4.0 * tr_xx_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_xx[i] * tbe_0 + 8.0 * tr_xxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xyy[i] = -2.0 * tr_xx_xyy[i] * tbe_0 - 2.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xx_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xyz[i] = -2.0 * tr_xx_xyz[i] * tbe_0 - 6.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xx_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_xy[i] * tbe_0 + 8.0 * tr_xxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_xzz[i] = 2.0 * tr_xx_x[i] - 2.0 * tr_xx_xzz[i] * tbe_0 - 10.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xx_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxz_xz[i] * tbe_0 + 8.0 * tr_xxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_yyy[i] = -2.0 * tr_xx_yyy[i] * tbe_0 - 2.0 * tr_xx_yyy[i] * tke_0 + 4.0 * tr_xx_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_yyz[i] = -2.0 * tr_xx_yyz[i] * tbe_0 - 6.0 * tr_xx_yyz[i] * tke_0 + 4.0 * tr_xx_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxz_yy[i] * tbe_0 + 8.0 * tr_xxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_yzz[i] = 2.0 * tr_xx_y[i] - 2.0 * tr_xx_yzz[i] * tbe_0 - 10.0 * tr_xx_yzz[i] * tke_0 + 4.0 * tr_xx_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxz_yz[i] * tbe_0 + 8.0 * tr_xxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xx_zzz[i] = 6.0 * tr_xx_z[i] - 2.0 * tr_xx_zzz[i] * tbe_0 - 14.0 * tr_xx_zzz[i] * tke_0 + 4.0 * tr_xx_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxz_zz[i] * tbe_0 + 8.0 * tr_xxz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 310-320 components of targeted buffer : DF

    auto tr_0_0_zz_xy_xxx = pbuffer.data(idx_op_geom_020_df + 310);

    auto tr_0_0_zz_xy_xxy = pbuffer.data(idx_op_geom_020_df + 311);

    auto tr_0_0_zz_xy_xxz = pbuffer.data(idx_op_geom_020_df + 312);

    auto tr_0_0_zz_xy_xyy = pbuffer.data(idx_op_geom_020_df + 313);

    auto tr_0_0_zz_xy_xyz = pbuffer.data(idx_op_geom_020_df + 314);

    auto tr_0_0_zz_xy_xzz = pbuffer.data(idx_op_geom_020_df + 315);

    auto tr_0_0_zz_xy_yyy = pbuffer.data(idx_op_geom_020_df + 316);

    auto tr_0_0_zz_xy_yyz = pbuffer.data(idx_op_geom_020_df + 317);

    auto tr_0_0_zz_xy_yzz = pbuffer.data(idx_op_geom_020_df + 318);

    auto tr_0_0_zz_xy_zzz = pbuffer.data(idx_op_geom_020_df + 319);

    #pragma omp simd aligned(tr_0_0_zz_xy_xxx, tr_0_0_zz_xy_xxy, tr_0_0_zz_xy_xxz, tr_0_0_zz_xy_xyy, tr_0_0_zz_xy_xyz, tr_0_0_zz_xy_xzz, tr_0_0_zz_xy_yyy, tr_0_0_zz_xy_yyz, tr_0_0_zz_xy_yzz, tr_0_0_zz_xy_zzz, tr_xy_x, tr_xy_xxx, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_y, tr_xy_yyy, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_z, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xy_xxx[i] = -2.0 * tr_xy_xxx[i] * tbe_0 - 2.0 * tr_xy_xxx[i] * tke_0 + 4.0 * tr_xy_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xxy[i] = -2.0 * tr_xy_xxy[i] * tbe_0 - 2.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xxz[i] = -2.0 * tr_xy_xxz[i] * tbe_0 - 6.0 * tr_xy_xxz[i] * tke_0 + 4.0 * tr_xy_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 8.0 * tr_xyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xyy[i] = -2.0 * tr_xy_xyy[i] * tbe_0 - 2.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xyz[i] = -2.0 * tr_xy_xyz[i] * tbe_0 - 6.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xy_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 8.0 * tr_xyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_xzz[i] = 2.0 * tr_xy_x[i] - 2.0 * tr_xy_xzz[i] * tbe_0 - 10.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xy_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyz_xz[i] * tbe_0 + 8.0 * tr_xyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_yyy[i] = -2.0 * tr_xy_yyy[i] * tbe_0 - 2.0 * tr_xy_yyy[i] * tke_0 + 4.0 * tr_xy_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_yyz[i] = -2.0 * tr_xy_yyz[i] * tbe_0 - 6.0 * tr_xy_yyz[i] * tke_0 + 4.0 * tr_xy_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 8.0 * tr_xyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_yzz[i] = 2.0 * tr_xy_y[i] - 2.0 * tr_xy_yzz[i] * tbe_0 - 10.0 * tr_xy_yzz[i] * tke_0 + 4.0 * tr_xy_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyz_yz[i] * tbe_0 + 8.0 * tr_xyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xy_zzz[i] = 6.0 * tr_xy_z[i] - 2.0 * tr_xy_zzz[i] * tbe_0 - 14.0 * tr_xy_zzz[i] * tke_0 + 4.0 * tr_xy_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyz_zz[i] * tbe_0 + 8.0 * tr_xyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 320-330 components of targeted buffer : DF

    auto tr_0_0_zz_xz_xxx = pbuffer.data(idx_op_geom_020_df + 320);

    auto tr_0_0_zz_xz_xxy = pbuffer.data(idx_op_geom_020_df + 321);

    auto tr_0_0_zz_xz_xxz = pbuffer.data(idx_op_geom_020_df + 322);

    auto tr_0_0_zz_xz_xyy = pbuffer.data(idx_op_geom_020_df + 323);

    auto tr_0_0_zz_xz_xyz = pbuffer.data(idx_op_geom_020_df + 324);

    auto tr_0_0_zz_xz_xzz = pbuffer.data(idx_op_geom_020_df + 325);

    auto tr_0_0_zz_xz_yyy = pbuffer.data(idx_op_geom_020_df + 326);

    auto tr_0_0_zz_xz_yyz = pbuffer.data(idx_op_geom_020_df + 327);

    auto tr_0_0_zz_xz_yzz = pbuffer.data(idx_op_geom_020_df + 328);

    auto tr_0_0_zz_xz_zzz = pbuffer.data(idx_op_geom_020_df + 329);

    #pragma omp simd aligned(tr_0_0_zz_xz_xxx, tr_0_0_zz_xz_xxy, tr_0_0_zz_xz_xxz, tr_0_0_zz_xz_xyy, tr_0_0_zz_xz_xyz, tr_0_0_zz_xz_xzz, tr_0_0_zz_xz_yyy, tr_0_0_zz_xz_yyz, tr_0_0_zz_xz_yzz, tr_0_0_zz_xz_zzz, tr_x_xx, tr_x_xxxz, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xz_x, tr_xz_xxx, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_y, tr_xz_yyy, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_z, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xx, tr_xzz_xxxz, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xz_xxx[i] = -4.0 * tr_x_xxxz[i] * tke_0 - 6.0 * tr_xz_xxx[i] * tbe_0 - 2.0 * tr_xz_xxx[i] * tke_0 + 4.0 * tr_xz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xxy[i] = -4.0 * tr_x_xxyz[i] * tke_0 - 6.0 * tr_xz_xxy[i] * tbe_0 - 2.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xxz[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_x_xxzz[i] * tke_0 - 6.0 * tr_xz_xxz[i] * tbe_0 - 6.0 * tr_xz_xxz[i] * tke_0 + 4.0 * tr_xz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_xx[i] * tbe_0 + 8.0 * tr_xzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xyy[i] = -4.0 * tr_x_xyyz[i] * tke_0 - 6.0 * tr_xz_xyy[i] * tbe_0 - 2.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xyz[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_x_xyzz[i] * tke_0 - 6.0 * tr_xz_xyz[i] * tbe_0 - 6.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_xy[i] * tbe_0 + 8.0 * tr_xzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_xzz[i] = 4.0 * tr_x_xz[i] - 4.0 * tr_x_xzzz[i] * tke_0 + 2.0 * tr_xz_x[i] - 6.0 * tr_xz_xzz[i] * tbe_0 - 10.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzz_xz[i] * tbe_0 + 8.0 * tr_xzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_yyy[i] = -4.0 * tr_x_yyyz[i] * tke_0 - 6.0 * tr_xz_yyy[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tke_0 + 4.0 * tr_xz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_yyz[i] = 2.0 * tr_x_yy[i] - 4.0 * tr_x_yyzz[i] * tke_0 - 6.0 * tr_xz_yyz[i] * tbe_0 - 6.0 * tr_xz_yyz[i] * tke_0 + 4.0 * tr_xz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzz_yy[i] * tbe_0 + 8.0 * tr_xzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_yzz[i] = 4.0 * tr_x_yz[i] - 4.0 * tr_x_yzzz[i] * tke_0 + 2.0 * tr_xz_y[i] - 6.0 * tr_xz_yzz[i] * tbe_0 - 10.0 * tr_xz_yzz[i] * tke_0 + 4.0 * tr_xz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzz_yz[i] * tbe_0 + 8.0 * tr_xzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xz_zzz[i] = 6.0 * tr_x_zz[i] - 4.0 * tr_x_zzzz[i] * tke_0 + 6.0 * tr_xz_z[i] - 6.0 * tr_xz_zzz[i] * tbe_0 - 14.0 * tr_xz_zzz[i] * tke_0 + 4.0 * tr_xz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xzz_zz[i] * tbe_0 + 8.0 * tr_xzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-340 components of targeted buffer : DF

    auto tr_0_0_zz_yy_xxx = pbuffer.data(idx_op_geom_020_df + 330);

    auto tr_0_0_zz_yy_xxy = pbuffer.data(idx_op_geom_020_df + 331);

    auto tr_0_0_zz_yy_xxz = pbuffer.data(idx_op_geom_020_df + 332);

    auto tr_0_0_zz_yy_xyy = pbuffer.data(idx_op_geom_020_df + 333);

    auto tr_0_0_zz_yy_xyz = pbuffer.data(idx_op_geom_020_df + 334);

    auto tr_0_0_zz_yy_xzz = pbuffer.data(idx_op_geom_020_df + 335);

    auto tr_0_0_zz_yy_yyy = pbuffer.data(idx_op_geom_020_df + 336);

    auto tr_0_0_zz_yy_yyz = pbuffer.data(idx_op_geom_020_df + 337);

    auto tr_0_0_zz_yy_yzz = pbuffer.data(idx_op_geom_020_df + 338);

    auto tr_0_0_zz_yy_zzz = pbuffer.data(idx_op_geom_020_df + 339);

    #pragma omp simd aligned(tr_0_0_zz_yy_xxx, tr_0_0_zz_yy_xxy, tr_0_0_zz_yy_xxz, tr_0_0_zz_yy_xyy, tr_0_0_zz_yy_xyz, tr_0_0_zz_yy_xzz, tr_0_0_zz_yy_yyy, tr_0_0_zz_yy_yyz, tr_0_0_zz_yy_yzz, tr_0_0_zz_yy_zzz, tr_yy_x, tr_yy_xxx, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_y, tr_yy_yyy, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_z, tr_yy_zzz, tr_yy_zzzzz, tr_yyz_xx, tr_yyz_xxxz, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yy_xxx[i] = -2.0 * tr_yy_xxx[i] * tbe_0 - 2.0 * tr_yy_xxx[i] * tke_0 + 4.0 * tr_yy_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xxy[i] = -2.0 * tr_yy_xxy[i] * tbe_0 - 2.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yy_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xxz[i] = -2.0 * tr_yy_xxz[i] * tbe_0 - 6.0 * tr_yy_xxz[i] * tke_0 + 4.0 * tr_yy_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xx[i] * tbe_0 + 8.0 * tr_yyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xyy[i] = -2.0 * tr_yy_xyy[i] * tbe_0 - 2.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yy_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xyz[i] = -2.0 * tr_yy_xyz[i] * tbe_0 - 6.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yy_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_xy[i] * tbe_0 + 8.0 * tr_yyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_xzz[i] = 2.0 * tr_yy_x[i] - 2.0 * tr_yy_xzz[i] * tbe_0 - 10.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yy_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyz_xz[i] * tbe_0 + 8.0 * tr_yyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_yyy[i] = -2.0 * tr_yy_yyy[i] * tbe_0 - 2.0 * tr_yy_yyy[i] * tke_0 + 4.0 * tr_yy_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_yyz[i] = -2.0 * tr_yy_yyz[i] * tbe_0 - 6.0 * tr_yy_yyz[i] * tke_0 + 4.0 * tr_yy_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyz_yy[i] * tbe_0 + 8.0 * tr_yyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_yzz[i] = 2.0 * tr_yy_y[i] - 2.0 * tr_yy_yzz[i] * tbe_0 - 10.0 * tr_yy_yzz[i] * tke_0 + 4.0 * tr_yy_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyz_yz[i] * tbe_0 + 8.0 * tr_yyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yy_zzz[i] = 6.0 * tr_yy_z[i] - 2.0 * tr_yy_zzz[i] * tbe_0 - 14.0 * tr_yy_zzz[i] * tke_0 + 4.0 * tr_yy_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyz_zz[i] * tbe_0 + 8.0 * tr_yyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 340-350 components of targeted buffer : DF

    auto tr_0_0_zz_yz_xxx = pbuffer.data(idx_op_geom_020_df + 340);

    auto tr_0_0_zz_yz_xxy = pbuffer.data(idx_op_geom_020_df + 341);

    auto tr_0_0_zz_yz_xxz = pbuffer.data(idx_op_geom_020_df + 342);

    auto tr_0_0_zz_yz_xyy = pbuffer.data(idx_op_geom_020_df + 343);

    auto tr_0_0_zz_yz_xyz = pbuffer.data(idx_op_geom_020_df + 344);

    auto tr_0_0_zz_yz_xzz = pbuffer.data(idx_op_geom_020_df + 345);

    auto tr_0_0_zz_yz_yyy = pbuffer.data(idx_op_geom_020_df + 346);

    auto tr_0_0_zz_yz_yyz = pbuffer.data(idx_op_geom_020_df + 347);

    auto tr_0_0_zz_yz_yzz = pbuffer.data(idx_op_geom_020_df + 348);

    auto tr_0_0_zz_yz_zzz = pbuffer.data(idx_op_geom_020_df + 349);

    #pragma omp simd aligned(tr_0_0_zz_yz_xxx, tr_0_0_zz_yz_xxy, tr_0_0_zz_yz_xxz, tr_0_0_zz_yz_xyy, tr_0_0_zz_yz_xyz, tr_0_0_zz_yz_xzz, tr_0_0_zz_yz_yyy, tr_0_0_zz_yz_yyz, tr_0_0_zz_yz_yzz, tr_0_0_zz_yz_zzz, tr_y_xx, tr_y_xxxz, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, tr_yz_x, tr_yz_xxx, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_y, tr_yz_yyy, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_z, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xx, tr_yzz_xxxz, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yz_xxx[i] = -4.0 * tr_y_xxxz[i] * tke_0 - 6.0 * tr_yz_xxx[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tke_0 + 4.0 * tr_yz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xxy[i] = -4.0 * tr_y_xxyz[i] * tke_0 - 6.0 * tr_yz_xxy[i] * tbe_0 - 2.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xxz[i] = 2.0 * tr_y_xx[i] - 4.0 * tr_y_xxzz[i] * tke_0 - 6.0 * tr_yz_xxz[i] * tbe_0 - 6.0 * tr_yz_xxz[i] * tke_0 + 4.0 * tr_yz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xx[i] * tbe_0 + 8.0 * tr_yzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xyy[i] = -4.0 * tr_y_xyyz[i] * tke_0 - 6.0 * tr_yz_xyy[i] * tbe_0 - 2.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xyz[i] = 2.0 * tr_y_xy[i] - 4.0 * tr_y_xyzz[i] * tke_0 - 6.0 * tr_yz_xyz[i] * tbe_0 - 6.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_xy[i] * tbe_0 + 8.0 * tr_yzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_xzz[i] = 4.0 * tr_y_xz[i] - 4.0 * tr_y_xzzz[i] * tke_0 + 2.0 * tr_yz_x[i] - 6.0 * tr_yz_xzz[i] * tbe_0 - 10.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzz_xz[i] * tbe_0 + 8.0 * tr_yzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_yyy[i] = -4.0 * tr_y_yyyz[i] * tke_0 - 6.0 * tr_yz_yyy[i] * tbe_0 - 2.0 * tr_yz_yyy[i] * tke_0 + 4.0 * tr_yz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_yyz[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_y_yyzz[i] * tke_0 - 6.0 * tr_yz_yyz[i] * tbe_0 - 6.0 * tr_yz_yyz[i] * tke_0 + 4.0 * tr_yz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzz_yy[i] * tbe_0 + 8.0 * tr_yzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_yzz[i] = 4.0 * tr_y_yz[i] - 4.0 * tr_y_yzzz[i] * tke_0 + 2.0 * tr_yz_y[i] - 6.0 * tr_yz_yzz[i] * tbe_0 - 10.0 * tr_yz_yzz[i] * tke_0 + 4.0 * tr_yz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzz_yz[i] * tbe_0 + 8.0 * tr_yzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yz_zzz[i] = 6.0 * tr_y_zz[i] - 4.0 * tr_y_zzzz[i] * tke_0 + 6.0 * tr_yz_z[i] - 6.0 * tr_yz_zzz[i] * tbe_0 - 14.0 * tr_yz_zzz[i] * tke_0 + 4.0 * tr_yz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yzz_zz[i] * tbe_0 + 8.0 * tr_yzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 350-360 components of targeted buffer : DF

    auto tr_0_0_zz_zz_xxx = pbuffer.data(idx_op_geom_020_df + 350);

    auto tr_0_0_zz_zz_xxy = pbuffer.data(idx_op_geom_020_df + 351);

    auto tr_0_0_zz_zz_xxz = pbuffer.data(idx_op_geom_020_df + 352);

    auto tr_0_0_zz_zz_xyy = pbuffer.data(idx_op_geom_020_df + 353);

    auto tr_0_0_zz_zz_xyz = pbuffer.data(idx_op_geom_020_df + 354);

    auto tr_0_0_zz_zz_xzz = pbuffer.data(idx_op_geom_020_df + 355);

    auto tr_0_0_zz_zz_yyy = pbuffer.data(idx_op_geom_020_df + 356);

    auto tr_0_0_zz_zz_yyz = pbuffer.data(idx_op_geom_020_df + 357);

    auto tr_0_0_zz_zz_yzz = pbuffer.data(idx_op_geom_020_df + 358);

    auto tr_0_0_zz_zz_zzz = pbuffer.data(idx_op_geom_020_df + 359);

    #pragma omp simd aligned(tr_0_0_zz_zz_xxx, tr_0_0_zz_zz_xxy, tr_0_0_zz_zz_xxz, tr_0_0_zz_zz_xyy, tr_0_0_zz_zz_xyz, tr_0_0_zz_zz_xzz, tr_0_0_zz_zz_yyy, tr_0_0_zz_zz_yyz, tr_0_0_zz_zz_yzz, tr_0_0_zz_zz_zzz, tr_0_xxx, tr_0_xxy, tr_0_xxz, tr_0_xyy, tr_0_xyz, tr_0_xzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_z_xx, tr_z_xxxz, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, tr_zz_x, tr_zz_xxx, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_y, tr_zz_yyy, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_z, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xx, tr_zzz_xxxz, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_zz_xxx[i] = 2.0 * tr_0_xxx[i] - 8.0 * tr_z_xxxz[i] * tke_0 - 10.0 * tr_zz_xxx[i] * tbe_0 - 2.0 * tr_zz_xxx[i] * tke_0 + 4.0 * tr_zz_xxxzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xxy[i] = 2.0 * tr_0_xxy[i] - 8.0 * tr_z_xxyz[i] * tke_0 - 10.0 * tr_zz_xxy[i] * tbe_0 - 2.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zz_xxyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xxz[i] = 2.0 * tr_0_xxz[i] + 4.0 * tr_z_xx[i] - 8.0 * tr_z_xxzz[i] * tke_0 - 10.0 * tr_zz_xxz[i] * tbe_0 - 6.0 * tr_zz_xxz[i] * tke_0 + 4.0 * tr_zz_xxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_xx[i] * tbe_0 + 8.0 * tr_zzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xyy[i] = 2.0 * tr_0_xyy[i] - 8.0 * tr_z_xyyz[i] * tke_0 - 10.0 * tr_zz_xyy[i] * tbe_0 - 2.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zz_xyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xyz[i] = 2.0 * tr_0_xyz[i] + 4.0 * tr_z_xy[i] - 8.0 * tr_z_xyzz[i] * tke_0 - 10.0 * tr_zz_xyz[i] * tbe_0 - 6.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zz_xyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_xy[i] * tbe_0 + 8.0 * tr_zzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_xzz[i] = 2.0 * tr_0_xzz[i] + 8.0 * tr_z_xz[i] - 8.0 * tr_z_xzzz[i] * tke_0 + 2.0 * tr_zz_x[i] - 10.0 * tr_zz_xzz[i] * tbe_0 - 10.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_zz_xzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzz_xz[i] * tbe_0 + 8.0 * tr_zzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_xzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_yyy[i] = 2.0 * tr_0_yyy[i] - 8.0 * tr_z_yyyz[i] * tke_0 - 10.0 * tr_zz_yyy[i] * tbe_0 - 2.0 * tr_zz_yyy[i] * tke_0 + 4.0 * tr_zz_yyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_yyz[i] = 2.0 * tr_0_yyz[i] + 4.0 * tr_z_yy[i] - 8.0 * tr_z_yyzz[i] * tke_0 - 10.0 * tr_zz_yyz[i] * tbe_0 - 6.0 * tr_zz_yyz[i] * tke_0 + 4.0 * tr_zz_yyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzz_yy[i] * tbe_0 + 8.0 * tr_zzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_yzz[i] = 2.0 * tr_0_yzz[i] + 8.0 * tr_z_yz[i] - 8.0 * tr_z_yzzz[i] * tke_0 + 2.0 * tr_zz_y[i] - 10.0 * tr_zz_yzz[i] * tbe_0 - 10.0 * tr_zz_yzz[i] * tke_0 + 4.0 * tr_zz_yzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzz_yz[i] * tbe_0 + 8.0 * tr_zzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_yzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zz_zzz[i] = 2.0 * tr_0_zzz[i] + 12.0 * tr_z_zz[i] - 8.0 * tr_z_zzzz[i] * tke_0 + 6.0 * tr_zz_z[i] - 10.0 * tr_zz_zzz[i] * tbe_0 - 14.0 * tr_zz_zzz[i] * tke_0 + 4.0 * tr_zz_zzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zzz_zz[i] * tbe_0 + 8.0 * tr_zzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzz_zzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

