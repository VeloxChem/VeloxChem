#include "GeometricalDerivatives020ForPG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_pg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_pg,
                         const int idx_op_sf,
                         const int idx_op_sh,
                         const int idx_op_pd,
                         const int idx_op_pg,
                         const int idx_op_pi,
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

    // Set up components of auxiliary buffer : PI

    auto tr_x_xxxxxx = pbuffer.data(idx_op_pi);

    auto tr_x_xxxxxy = pbuffer.data(idx_op_pi + 1);

    auto tr_x_xxxxxz = pbuffer.data(idx_op_pi + 2);

    auto tr_x_xxxxyy = pbuffer.data(idx_op_pi + 3);

    auto tr_x_xxxxyz = pbuffer.data(idx_op_pi + 4);

    auto tr_x_xxxxzz = pbuffer.data(idx_op_pi + 5);

    auto tr_x_xxxyyy = pbuffer.data(idx_op_pi + 6);

    auto tr_x_xxxyyz = pbuffer.data(idx_op_pi + 7);

    auto tr_x_xxxyzz = pbuffer.data(idx_op_pi + 8);

    auto tr_x_xxxzzz = pbuffer.data(idx_op_pi + 9);

    auto tr_x_xxyyyy = pbuffer.data(idx_op_pi + 10);

    auto tr_x_xxyyyz = pbuffer.data(idx_op_pi + 11);

    auto tr_x_xxyyzz = pbuffer.data(idx_op_pi + 12);

    auto tr_x_xxyzzz = pbuffer.data(idx_op_pi + 13);

    auto tr_x_xxzzzz = pbuffer.data(idx_op_pi + 14);

    auto tr_x_xyyyyy = pbuffer.data(idx_op_pi + 15);

    auto tr_x_xyyyyz = pbuffer.data(idx_op_pi + 16);

    auto tr_x_xyyyzz = pbuffer.data(idx_op_pi + 17);

    auto tr_x_xyyzzz = pbuffer.data(idx_op_pi + 18);

    auto tr_x_xyzzzz = pbuffer.data(idx_op_pi + 19);

    auto tr_x_xzzzzz = pbuffer.data(idx_op_pi + 20);

    auto tr_x_yyyyyy = pbuffer.data(idx_op_pi + 21);

    auto tr_x_yyyyyz = pbuffer.data(idx_op_pi + 22);

    auto tr_x_yyyyzz = pbuffer.data(idx_op_pi + 23);

    auto tr_x_yyyzzz = pbuffer.data(idx_op_pi + 24);

    auto tr_x_yyzzzz = pbuffer.data(idx_op_pi + 25);

    auto tr_x_yzzzzz = pbuffer.data(idx_op_pi + 26);

    auto tr_x_zzzzzz = pbuffer.data(idx_op_pi + 27);

    auto tr_y_xxxxxx = pbuffer.data(idx_op_pi + 28);

    auto tr_y_xxxxxy = pbuffer.data(idx_op_pi + 29);

    auto tr_y_xxxxxz = pbuffer.data(idx_op_pi + 30);

    auto tr_y_xxxxyy = pbuffer.data(idx_op_pi + 31);

    auto tr_y_xxxxyz = pbuffer.data(idx_op_pi + 32);

    auto tr_y_xxxxzz = pbuffer.data(idx_op_pi + 33);

    auto tr_y_xxxyyy = pbuffer.data(idx_op_pi + 34);

    auto tr_y_xxxyyz = pbuffer.data(idx_op_pi + 35);

    auto tr_y_xxxyzz = pbuffer.data(idx_op_pi + 36);

    auto tr_y_xxxzzz = pbuffer.data(idx_op_pi + 37);

    auto tr_y_xxyyyy = pbuffer.data(idx_op_pi + 38);

    auto tr_y_xxyyyz = pbuffer.data(idx_op_pi + 39);

    auto tr_y_xxyyzz = pbuffer.data(idx_op_pi + 40);

    auto tr_y_xxyzzz = pbuffer.data(idx_op_pi + 41);

    auto tr_y_xxzzzz = pbuffer.data(idx_op_pi + 42);

    auto tr_y_xyyyyy = pbuffer.data(idx_op_pi + 43);

    auto tr_y_xyyyyz = pbuffer.data(idx_op_pi + 44);

    auto tr_y_xyyyzz = pbuffer.data(idx_op_pi + 45);

    auto tr_y_xyyzzz = pbuffer.data(idx_op_pi + 46);

    auto tr_y_xyzzzz = pbuffer.data(idx_op_pi + 47);

    auto tr_y_xzzzzz = pbuffer.data(idx_op_pi + 48);

    auto tr_y_yyyyyy = pbuffer.data(idx_op_pi + 49);

    auto tr_y_yyyyyz = pbuffer.data(idx_op_pi + 50);

    auto tr_y_yyyyzz = pbuffer.data(idx_op_pi + 51);

    auto tr_y_yyyzzz = pbuffer.data(idx_op_pi + 52);

    auto tr_y_yyzzzz = pbuffer.data(idx_op_pi + 53);

    auto tr_y_yzzzzz = pbuffer.data(idx_op_pi + 54);

    auto tr_y_zzzzzz = pbuffer.data(idx_op_pi + 55);

    auto tr_z_xxxxxx = pbuffer.data(idx_op_pi + 56);

    auto tr_z_xxxxxy = pbuffer.data(idx_op_pi + 57);

    auto tr_z_xxxxxz = pbuffer.data(idx_op_pi + 58);

    auto tr_z_xxxxyy = pbuffer.data(idx_op_pi + 59);

    auto tr_z_xxxxyz = pbuffer.data(idx_op_pi + 60);

    auto tr_z_xxxxzz = pbuffer.data(idx_op_pi + 61);

    auto tr_z_xxxyyy = pbuffer.data(idx_op_pi + 62);

    auto tr_z_xxxyyz = pbuffer.data(idx_op_pi + 63);

    auto tr_z_xxxyzz = pbuffer.data(idx_op_pi + 64);

    auto tr_z_xxxzzz = pbuffer.data(idx_op_pi + 65);

    auto tr_z_xxyyyy = pbuffer.data(idx_op_pi + 66);

    auto tr_z_xxyyyz = pbuffer.data(idx_op_pi + 67);

    auto tr_z_xxyyzz = pbuffer.data(idx_op_pi + 68);

    auto tr_z_xxyzzz = pbuffer.data(idx_op_pi + 69);

    auto tr_z_xxzzzz = pbuffer.data(idx_op_pi + 70);

    auto tr_z_xyyyyy = pbuffer.data(idx_op_pi + 71);

    auto tr_z_xyyyyz = pbuffer.data(idx_op_pi + 72);

    auto tr_z_xyyyzz = pbuffer.data(idx_op_pi + 73);

    auto tr_z_xyyzzz = pbuffer.data(idx_op_pi + 74);

    auto tr_z_xyzzzz = pbuffer.data(idx_op_pi + 75);

    auto tr_z_xzzzzz = pbuffer.data(idx_op_pi + 76);

    auto tr_z_yyyyyy = pbuffer.data(idx_op_pi + 77);

    auto tr_z_yyyyyz = pbuffer.data(idx_op_pi + 78);

    auto tr_z_yyyyzz = pbuffer.data(idx_op_pi + 79);

    auto tr_z_yyyzzz = pbuffer.data(idx_op_pi + 80);

    auto tr_z_yyzzzz = pbuffer.data(idx_op_pi + 81);

    auto tr_z_yzzzzz = pbuffer.data(idx_op_pi + 82);

    auto tr_z_zzzzzz = pbuffer.data(idx_op_pi + 83);

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

    auto tr_0_0_xx_x_xxxx = pbuffer.data(idx_op_geom_020_pg);

    auto tr_0_0_xx_x_xxxy = pbuffer.data(idx_op_geom_020_pg + 1);

    auto tr_0_0_xx_x_xxxz = pbuffer.data(idx_op_geom_020_pg + 2);

    auto tr_0_0_xx_x_xxyy = pbuffer.data(idx_op_geom_020_pg + 3);

    auto tr_0_0_xx_x_xxyz = pbuffer.data(idx_op_geom_020_pg + 4);

    auto tr_0_0_xx_x_xxzz = pbuffer.data(idx_op_geom_020_pg + 5);

    auto tr_0_0_xx_x_xyyy = pbuffer.data(idx_op_geom_020_pg + 6);

    auto tr_0_0_xx_x_xyyz = pbuffer.data(idx_op_geom_020_pg + 7);

    auto tr_0_0_xx_x_xyzz = pbuffer.data(idx_op_geom_020_pg + 8);

    auto tr_0_0_xx_x_xzzz = pbuffer.data(idx_op_geom_020_pg + 9);

    auto tr_0_0_xx_x_yyyy = pbuffer.data(idx_op_geom_020_pg + 10);

    auto tr_0_0_xx_x_yyyz = pbuffer.data(idx_op_geom_020_pg + 11);

    auto tr_0_0_xx_x_yyzz = pbuffer.data(idx_op_geom_020_pg + 12);

    auto tr_0_0_xx_x_yzzz = pbuffer.data(idx_op_geom_020_pg + 13);

    auto tr_0_0_xx_x_zzzz = pbuffer.data(idx_op_geom_020_pg + 14);

    #pragma omp simd aligned(tr_0_0_xx_x_xxxx, tr_0_0_xx_x_xxxy, tr_0_0_xx_x_xxxz, tr_0_0_xx_x_xxyy, tr_0_0_xx_x_xxyz, tr_0_0_xx_x_xxzz, tr_0_0_xx_x_xyyy, tr_0_0_xx_x_xyyz, tr_0_0_xx_x_xyzz, tr_0_0_xx_x_xzzz, tr_0_0_xx_x_yyyy, tr_0_0_xx_x_yyyz, tr_0_0_xx_x_yyzz, tr_0_0_xx_x_yzzz, tr_0_0_xx_x_zzzz, tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xx, tr_x_xxxx, tr_x_xxxxxx, tr_x_xxxxxy, tr_x_xxxxxz, tr_x_xxxxyy, tr_x_xxxxyz, tr_x_xxxxzz, tr_x_xxxy, tr_x_xxxyyy, tr_x_xxxyyz, tr_x_xxxyzz, tr_x_xxxz, tr_x_xxxzzz, tr_x_xxyy, tr_x_xxyyyy, tr_x_xxyyyz, tr_x_xxyyzz, tr_x_xxyz, tr_x_xxyzzz, tr_x_xxzz, tr_x_xxzzzz, tr_x_xy, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxx, tr_xx_xxxxy, tr_xx_xxxxz, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xzzz, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yzzz, tr_xxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_x_xxxx[i] = 8.0 * tr_0_xxx[i] - 4.0 * tr_0_xxxxx[i] * tke_0 + 12.0 * tr_x_xx[i] - 6.0 * tr_x_xxxx[i] * tbe_0 - 18.0 * tr_x_xxxx[i] * tke_0 + 4.0 * tr_x_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xx_xxx[i] * tbe_0 + 8.0 * tr_xx_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxxy[i] = 6.0 * tr_0_xxy[i] - 4.0 * tr_0_xxxxy[i] * tke_0 + 6.0 * tr_x_xy[i] - 6.0 * tr_x_xxxy[i] * tbe_0 - 14.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xx_xxy[i] * tbe_0 + 8.0 * tr_xx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxxz[i] = 6.0 * tr_0_xxz[i] - 4.0 * tr_0_xxxxz[i] * tke_0 + 6.0 * tr_x_xz[i] - 6.0 * tr_x_xxxz[i] * tbe_0 - 14.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xx_xxz[i] * tbe_0 + 8.0 * tr_xx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxyy[i] = 4.0 * tr_0_xyy[i] - 4.0 * tr_0_xxxyy[i] * tke_0 + 2.0 * tr_x_yy[i] - 6.0 * tr_x_xxyy[i] * tbe_0 - 10.0 * tr_x_xxyy[i] * tke_0 + 4.0 * tr_x_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xx_xyy[i] * tbe_0 + 8.0 * tr_xx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxyz[i] = 4.0 * tr_0_xyz[i] - 4.0 * tr_0_xxxyz[i] * tke_0 + 2.0 * tr_x_yz[i] - 6.0 * tr_x_xxyz[i] * tbe_0 - 10.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xx_xyz[i] * tbe_0 + 8.0 * tr_xx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xxzz[i] = 4.0 * tr_0_xzz[i] - 4.0 * tr_0_xxxzz[i] * tke_0 + 2.0 * tr_x_zz[i] - 6.0 * tr_x_xxzz[i] * tbe_0 - 10.0 * tr_x_xxzz[i] * tke_0 + 4.0 * tr_x_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xx_xzz[i] * tbe_0 + 8.0 * tr_xx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xyyy[i] = 2.0 * tr_0_yyy[i] - 4.0 * tr_0_xxyyy[i] * tke_0 - 6.0 * tr_x_xyyy[i] * tbe_0 - 6.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_x_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xx_yyy[i] * tbe_0 + 8.0 * tr_xx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xyyz[i] = 2.0 * tr_0_yyz[i] - 4.0 * tr_0_xxyyz[i] * tke_0 - 6.0 * tr_x_xyyz[i] * tbe_0 - 6.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xx_yyz[i] * tbe_0 + 8.0 * tr_xx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xyzz[i] = 2.0 * tr_0_yzz[i] - 4.0 * tr_0_xxyzz[i] * tke_0 - 6.0 * tr_x_xyzz[i] * tbe_0 - 6.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xx_yzz[i] * tbe_0 + 8.0 * tr_xx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_xzzz[i] = 2.0 * tr_0_zzz[i] - 4.0 * tr_0_xxzzz[i] * tke_0 - 6.0 * tr_x_xzzz[i] * tbe_0 - 6.0 * tr_x_xzzz[i] * tke_0 + 4.0 * tr_x_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xx_zzz[i] * tbe_0 + 8.0 * tr_xx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yyyy[i] = -4.0 * tr_0_xyyyy[i] * tke_0 - 6.0 * tr_x_yyyy[i] * tbe_0 - 2.0 * tr_x_yyyy[i] * tke_0 + 4.0 * tr_x_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yyyz[i] = -4.0 * tr_0_xyyyz[i] * tke_0 - 6.0 * tr_x_yyyz[i] * tbe_0 - 2.0 * tr_x_yyyz[i] * tke_0 + 4.0 * tr_x_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yyzz[i] = -4.0 * tr_0_xyyzz[i] * tke_0 - 6.0 * tr_x_yyzz[i] * tbe_0 - 2.0 * tr_x_yyzz[i] * tke_0 + 4.0 * tr_x_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_yzzz[i] = -4.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_x_yzzz[i] * tbe_0 - 2.0 * tr_x_yzzz[i] * tke_0 + 4.0 * tr_x_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_x_zzzz[i] = -4.0 * tr_0_xzzzz[i] * tke_0 - 6.0 * tr_x_zzzz[i] * tbe_0 - 2.0 * tr_x_zzzz[i] * tke_0 + 4.0 * tr_x_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxx_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 15-30 components of targeted buffer : PG

    auto tr_0_0_xx_y_xxxx = pbuffer.data(idx_op_geom_020_pg + 15);

    auto tr_0_0_xx_y_xxxy = pbuffer.data(idx_op_geom_020_pg + 16);

    auto tr_0_0_xx_y_xxxz = pbuffer.data(idx_op_geom_020_pg + 17);

    auto tr_0_0_xx_y_xxyy = pbuffer.data(idx_op_geom_020_pg + 18);

    auto tr_0_0_xx_y_xxyz = pbuffer.data(idx_op_geom_020_pg + 19);

    auto tr_0_0_xx_y_xxzz = pbuffer.data(idx_op_geom_020_pg + 20);

    auto tr_0_0_xx_y_xyyy = pbuffer.data(idx_op_geom_020_pg + 21);

    auto tr_0_0_xx_y_xyyz = pbuffer.data(idx_op_geom_020_pg + 22);

    auto tr_0_0_xx_y_xyzz = pbuffer.data(idx_op_geom_020_pg + 23);

    auto tr_0_0_xx_y_xzzz = pbuffer.data(idx_op_geom_020_pg + 24);

    auto tr_0_0_xx_y_yyyy = pbuffer.data(idx_op_geom_020_pg + 25);

    auto tr_0_0_xx_y_yyyz = pbuffer.data(idx_op_geom_020_pg + 26);

    auto tr_0_0_xx_y_yyzz = pbuffer.data(idx_op_geom_020_pg + 27);

    auto tr_0_0_xx_y_yzzz = pbuffer.data(idx_op_geom_020_pg + 28);

    auto tr_0_0_xx_y_zzzz = pbuffer.data(idx_op_geom_020_pg + 29);

    #pragma omp simd aligned(tr_0_0_xx_y_xxxx, tr_0_0_xx_y_xxxy, tr_0_0_xx_y_xxxz, tr_0_0_xx_y_xxyy, tr_0_0_xx_y_xxyz, tr_0_0_xx_y_xxzz, tr_0_0_xx_y_xyyy, tr_0_0_xx_y_xyyz, tr_0_0_xx_y_xyzz, tr_0_0_xx_y_xzzz, tr_0_0_xx_y_yyyy, tr_0_0_xx_y_yyyz, tr_0_0_xx_y_yyzz, tr_0_0_xx_y_yzzz, tr_0_0_xx_y_zzzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_xx, tr_y_xxxx, tr_y_xxxxxx, tr_y_xxxxxy, tr_y_xxxxxz, tr_y_xxxxyy, tr_y_xxxxyz, tr_y_xxxxzz, tr_y_xxxy, tr_y_xxxyyy, tr_y_xxxyyz, tr_y_xxxyzz, tr_y_xxxz, tr_y_xxxzzz, tr_y_xxyy, tr_y_xxyyyy, tr_y_xxyyyz, tr_y_xxyyzz, tr_y_xxyz, tr_y_xxyzzz, tr_y_xxzz, tr_y_xxzzzz, tr_y_xy, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_y_xxxx[i] = 12.0 * tr_y_xx[i] - 2.0 * tr_y_xxxx[i] * tbe_0 - 18.0 * tr_y_xxxx[i] * tke_0 + 4.0 * tr_y_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xy_xxx[i] * tbe_0 + 8.0 * tr_xy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxxy[i] = 6.0 * tr_y_xy[i] - 2.0 * tr_y_xxxy[i] * tbe_0 - 14.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xy_xxy[i] * tbe_0 + 8.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxxz[i] = 6.0 * tr_y_xz[i] - 2.0 * tr_y_xxxz[i] * tbe_0 - 14.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xy_xxz[i] * tbe_0 + 8.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxyy[i] = 2.0 * tr_y_yy[i] - 2.0 * tr_y_xxyy[i] * tbe_0 - 10.0 * tr_y_xxyy[i] * tke_0 + 4.0 * tr_y_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xy_xyy[i] * tbe_0 + 8.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxyz[i] = 2.0 * tr_y_yz[i] - 2.0 * tr_y_xxyz[i] * tbe_0 - 10.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xy_xyz[i] * tbe_0 + 8.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xxzz[i] = 2.0 * tr_y_zz[i] - 2.0 * tr_y_xxzz[i] * tbe_0 - 10.0 * tr_y_xxzz[i] * tke_0 + 4.0 * tr_y_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xy_xzz[i] * tbe_0 + 8.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xyyy[i] = -2.0 * tr_y_xyyy[i] * tbe_0 - 6.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_y_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xy_yyy[i] * tbe_0 + 8.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xyyz[i] = -2.0 * tr_y_xyyz[i] * tbe_0 - 6.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xy_yyz[i] * tbe_0 + 8.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xyzz[i] = -2.0 * tr_y_xyzz[i] * tbe_0 - 6.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_yzz[i] * tbe_0 + 8.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_xzzz[i] = -2.0 * tr_y_xzzz[i] * tbe_0 - 6.0 * tr_y_xzzz[i] * tke_0 + 4.0 * tr_y_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_zzz[i] * tbe_0 + 8.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yyyy[i] = -2.0 * tr_y_yyyy[i] * tbe_0 - 2.0 * tr_y_yyyy[i] * tke_0 + 4.0 * tr_y_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yyyz[i] = -2.0 * tr_y_yyyz[i] * tbe_0 - 2.0 * tr_y_yyyz[i] * tke_0 + 4.0 * tr_y_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yyzz[i] = -2.0 * tr_y_yyzz[i] * tbe_0 - 2.0 * tr_y_yyzz[i] * tke_0 + 4.0 * tr_y_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_yzzz[i] = -2.0 * tr_y_yzzz[i] * tbe_0 - 2.0 * tr_y_yzzz[i] * tke_0 + 4.0 * tr_y_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_y_zzzz[i] = -2.0 * tr_y_zzzz[i] * tbe_0 - 2.0 * tr_y_zzzz[i] * tke_0 + 4.0 * tr_y_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-45 components of targeted buffer : PG

    auto tr_0_0_xx_z_xxxx = pbuffer.data(idx_op_geom_020_pg + 30);

    auto tr_0_0_xx_z_xxxy = pbuffer.data(idx_op_geom_020_pg + 31);

    auto tr_0_0_xx_z_xxxz = pbuffer.data(idx_op_geom_020_pg + 32);

    auto tr_0_0_xx_z_xxyy = pbuffer.data(idx_op_geom_020_pg + 33);

    auto tr_0_0_xx_z_xxyz = pbuffer.data(idx_op_geom_020_pg + 34);

    auto tr_0_0_xx_z_xxzz = pbuffer.data(idx_op_geom_020_pg + 35);

    auto tr_0_0_xx_z_xyyy = pbuffer.data(idx_op_geom_020_pg + 36);

    auto tr_0_0_xx_z_xyyz = pbuffer.data(idx_op_geom_020_pg + 37);

    auto tr_0_0_xx_z_xyzz = pbuffer.data(idx_op_geom_020_pg + 38);

    auto tr_0_0_xx_z_xzzz = pbuffer.data(idx_op_geom_020_pg + 39);

    auto tr_0_0_xx_z_yyyy = pbuffer.data(idx_op_geom_020_pg + 40);

    auto tr_0_0_xx_z_yyyz = pbuffer.data(idx_op_geom_020_pg + 41);

    auto tr_0_0_xx_z_yyzz = pbuffer.data(idx_op_geom_020_pg + 42);

    auto tr_0_0_xx_z_yzzz = pbuffer.data(idx_op_geom_020_pg + 43);

    auto tr_0_0_xx_z_zzzz = pbuffer.data(idx_op_geom_020_pg + 44);

    #pragma omp simd aligned(tr_0_0_xx_z_xxxx, tr_0_0_xx_z_xxxy, tr_0_0_xx_z_xxxz, tr_0_0_xx_z_xxyy, tr_0_0_xx_z_xxyz, tr_0_0_xx_z_xxzz, tr_0_0_xx_z_xyyy, tr_0_0_xx_z_xyyz, tr_0_0_xx_z_xyzz, tr_0_0_xx_z_xzzz, tr_0_0_xx_z_yyyy, tr_0_0_xx_z_yyyz, tr_0_0_xx_z_yyzz, tr_0_0_xx_z_yzzz, tr_0_0_xx_z_zzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_xx, tr_z_xxxx, tr_z_xxxxxx, tr_z_xxxxxy, tr_z_xxxxxz, tr_z_xxxxyy, tr_z_xxxxyz, tr_z_xxxxzz, tr_z_xxxy, tr_z_xxxyyy, tr_z_xxxyyz, tr_z_xxxyzz, tr_z_xxxz, tr_z_xxxzzz, tr_z_xxyy, tr_z_xxyyyy, tr_z_xxyyyz, tr_z_xxyyzz, tr_z_xxyz, tr_z_xxyzzz, tr_z_xxzz, tr_z_xxzzzz, tr_z_xy, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_z_xxxx[i] = 12.0 * tr_z_xx[i] - 2.0 * tr_z_xxxx[i] * tbe_0 - 18.0 * tr_z_xxxx[i] * tke_0 + 4.0 * tr_z_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xz_xxx[i] * tbe_0 + 8.0 * tr_xz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxxy[i] = 6.0 * tr_z_xy[i] - 2.0 * tr_z_xxxy[i] * tbe_0 - 14.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xz_xxy[i] * tbe_0 + 8.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxxz[i] = 6.0 * tr_z_xz[i] - 2.0 * tr_z_xxxz[i] * tbe_0 - 14.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xz_xxz[i] * tbe_0 + 8.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxyy[i] = 2.0 * tr_z_yy[i] - 2.0 * tr_z_xxyy[i] * tbe_0 - 10.0 * tr_z_xxyy[i] * tke_0 + 4.0 * tr_z_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xz_xyy[i] * tbe_0 + 8.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxyz[i] = 2.0 * tr_z_yz[i] - 2.0 * tr_z_xxyz[i] * tbe_0 - 10.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xyz[i] * tbe_0 + 8.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xxzz[i] = 2.0 * tr_z_zz[i] - 2.0 * tr_z_xxzz[i] * tbe_0 - 10.0 * tr_z_xxzz[i] * tke_0 + 4.0 * tr_z_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xzz[i] * tbe_0 + 8.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xyyy[i] = -2.0 * tr_z_xyyy[i] * tbe_0 - 6.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_z_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xz_yyy[i] * tbe_0 + 8.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xyyz[i] = -2.0 * tr_z_xyyz[i] * tbe_0 - 6.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yyz[i] * tbe_0 + 8.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xyzz[i] = -2.0 * tr_z_xyzz[i] * tbe_0 - 6.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yzz[i] * tbe_0 + 8.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_xzzz[i] = -2.0 * tr_z_xzzz[i] * tbe_0 - 6.0 * tr_z_xzzz[i] * tke_0 + 4.0 * tr_z_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_zzz[i] * tbe_0 + 8.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yyyy[i] = -2.0 * tr_z_yyyy[i] * tbe_0 - 2.0 * tr_z_yyyy[i] * tke_0 + 4.0 * tr_z_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yyyz[i] = -2.0 * tr_z_yyyz[i] * tbe_0 - 2.0 * tr_z_yyyz[i] * tke_0 + 4.0 * tr_z_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yyzz[i] = -2.0 * tr_z_yyzz[i] * tbe_0 - 2.0 * tr_z_yyzz[i] * tke_0 + 4.0 * tr_z_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_yzzz[i] = -2.0 * tr_z_yzzz[i] * tbe_0 - 2.0 * tr_z_yzzz[i] * tke_0 + 4.0 * tr_z_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_z_zzzz[i] = -2.0 * tr_z_zzzz[i] * tbe_0 - 2.0 * tr_z_zzzz[i] * tke_0 + 4.0 * tr_z_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 45-60 components of targeted buffer : PG

    auto tr_0_0_xy_x_xxxx = pbuffer.data(idx_op_geom_020_pg + 45);

    auto tr_0_0_xy_x_xxxy = pbuffer.data(idx_op_geom_020_pg + 46);

    auto tr_0_0_xy_x_xxxz = pbuffer.data(idx_op_geom_020_pg + 47);

    auto tr_0_0_xy_x_xxyy = pbuffer.data(idx_op_geom_020_pg + 48);

    auto tr_0_0_xy_x_xxyz = pbuffer.data(idx_op_geom_020_pg + 49);

    auto tr_0_0_xy_x_xxzz = pbuffer.data(idx_op_geom_020_pg + 50);

    auto tr_0_0_xy_x_xyyy = pbuffer.data(idx_op_geom_020_pg + 51);

    auto tr_0_0_xy_x_xyyz = pbuffer.data(idx_op_geom_020_pg + 52);

    auto tr_0_0_xy_x_xyzz = pbuffer.data(idx_op_geom_020_pg + 53);

    auto tr_0_0_xy_x_xzzz = pbuffer.data(idx_op_geom_020_pg + 54);

    auto tr_0_0_xy_x_yyyy = pbuffer.data(idx_op_geom_020_pg + 55);

    auto tr_0_0_xy_x_yyyz = pbuffer.data(idx_op_geom_020_pg + 56);

    auto tr_0_0_xy_x_yyzz = pbuffer.data(idx_op_geom_020_pg + 57);

    auto tr_0_0_xy_x_yzzz = pbuffer.data(idx_op_geom_020_pg + 58);

    auto tr_0_0_xy_x_zzzz = pbuffer.data(idx_op_geom_020_pg + 59);

    #pragma omp simd aligned(tr_0_0_xy_x_xxxx, tr_0_0_xy_x_xxxy, tr_0_0_xy_x_xxxz, tr_0_0_xy_x_xxyy, tr_0_0_xy_x_xxyz, tr_0_0_xy_x_xxzz, tr_0_0_xy_x_xyyy, tr_0_0_xy_x_xyyz, tr_0_0_xy_x_xyzz, tr_0_0_xy_x_xzzz, tr_0_0_xy_x_yyyy, tr_0_0_xy_x_yyyz, tr_0_0_xy_x_yyzz, tr_0_0_xy_x_yzzz, tr_0_0_xy_x_zzzz, tr_0_xxx, tr_0_xxxxy, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_x_xx, tr_x_xxxx, tr_x_xxxxxy, tr_x_xxxxyy, tr_x_xxxxyz, tr_x_xxxy, tr_x_xxxyyy, tr_x_xxxyyz, tr_x_xxxyzz, tr_x_xxxz, tr_x_xxyy, tr_x_xxyyyy, tr_x_xxyyyz, tr_x_xxyyzz, tr_x_xxyz, tr_x_xxyzzz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyyyy, tr_x_xyyyyz, tr_x_xyyyzz, tr_x_xyyz, tr_x_xyyzzz, tr_x_xyzz, tr_x_xyzzzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_xx_xxx, tr_xx_xxxxy, tr_xx_xxxyy, tr_xx_xxxyz, tr_xx_xxy, tr_xx_xxyyy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xyy, tr_xx_xyyyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyyyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xzzz, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yzzz, tr_xxy_zzzz, tr_xy_xxx, tr_xy_xxxxx, tr_xy_xxxxy, tr_xy_xxxxz, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_x_xxxx[i] = -2.0 * tr_0_xxxxy[i] * tke_0 - 2.0 * tr_y_xxxx[i] * tbe_0 - 8.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxxy[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 2.0 * tr_y_xxxy[i] * tbe_0 + 3.0 * tr_x_xx[i] - 6.0 * tr_x_xxyy[i] * tke_0 - 2.0 * tr_x_xxxx[i] * tke_0 + 4.0 * tr_x_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxxz[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_y_xxxz[i] * tbe_0 - 6.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxyy[i] = 2.0 * tr_0_xxy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 2.0 * tr_y_xxyy[i] * tbe_0 + 4.0 * tr_x_xy[i] - 4.0 * tr_x_xyyy[i] * tke_0 - 4.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xx_xxy[i] * tbe_0 + 4.0 * tr_xx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxyz[i] = tr_0_xxz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_y_xxyz[i] * tbe_0 + 2.0 * tr_x_xz[i] - 4.0 * tr_x_xyyz[i] * tke_0 - 2.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxz[i] * tbe_0 + 4.0 * tr_xx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xxzz[i] = -2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_y_xxzz[i] * tbe_0 - 4.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xyyy[i] = 3.0 * tr_0_xyy[i] - 2.0 * tr_0_xyyyy[i] * tke_0 - 2.0 * tr_y_xyyy[i] * tbe_0 + 3.0 * tr_x_yy[i] - 2.0 * tr_x_yyyy[i] * tke_0 - 6.0 * tr_x_xxyy[i] * tke_0 + 4.0 * tr_x_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xx_xyy[i] * tbe_0 + 4.0 * tr_xx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xyyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_y_xyyz[i] * tbe_0 + 2.0 * tr_x_yz[i] - 2.0 * tr_x_yyyz[i] * tke_0 - 4.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_xyz[i] * tbe_0 + 4.0 * tr_xx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xyzz[i] = tr_0_xzz[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_y_xyzz[i] * tbe_0 + tr_x_zz[i] - 2.0 * tr_x_yyzz[i] * tke_0 - 2.0 * tr_x_xxzz[i] * tke_0 + 4.0 * tr_x_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xzz[i] * tbe_0 + 4.0 * tr_xx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_xzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_y_xzzz[i] * tbe_0 - 2.0 * tr_x_yzzz[i] * tke_0 + 4.0 * tr_x_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yyyy[i] = 4.0 * tr_0_yyy[i] - 2.0 * tr_0_yyyyy[i] * tke_0 - 2.0 * tr_y_yyyy[i] * tbe_0 - 8.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_x_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xx_yyy[i] * tbe_0 + 4.0 * tr_xx_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yyyz[i] = 3.0 * tr_0_yyz[i] - 2.0 * tr_0_yyyyz[i] * tke_0 - 2.0 * tr_y_yyyz[i] * tbe_0 - 6.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xx_yyz[i] * tbe_0 + 4.0 * tr_xx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yyzz[i] = 2.0 * tr_0_yzz[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_y_yyzz[i] * tbe_0 - 4.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_yzz[i] * tbe_0 + 4.0 * tr_xx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_yzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_y_yzzz[i] * tbe_0 - 2.0 * tr_x_xzzz[i] * tke_0 + 4.0 * tr_x_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_zzz[i] * tbe_0 + 4.0 * tr_xx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_x_zzzz[i] = -2.0 * tr_0_yzzzz[i] * tke_0 - 2.0 * tr_y_zzzz[i] * tbe_0 + 4.0 * tr_x_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-75 components of targeted buffer : PG

    auto tr_0_0_xy_y_xxxx = pbuffer.data(idx_op_geom_020_pg + 60);

    auto tr_0_0_xy_y_xxxy = pbuffer.data(idx_op_geom_020_pg + 61);

    auto tr_0_0_xy_y_xxxz = pbuffer.data(idx_op_geom_020_pg + 62);

    auto tr_0_0_xy_y_xxyy = pbuffer.data(idx_op_geom_020_pg + 63);

    auto tr_0_0_xy_y_xxyz = pbuffer.data(idx_op_geom_020_pg + 64);

    auto tr_0_0_xy_y_xxzz = pbuffer.data(idx_op_geom_020_pg + 65);

    auto tr_0_0_xy_y_xyyy = pbuffer.data(idx_op_geom_020_pg + 66);

    auto tr_0_0_xy_y_xyyz = pbuffer.data(idx_op_geom_020_pg + 67);

    auto tr_0_0_xy_y_xyzz = pbuffer.data(idx_op_geom_020_pg + 68);

    auto tr_0_0_xy_y_xzzz = pbuffer.data(idx_op_geom_020_pg + 69);

    auto tr_0_0_xy_y_yyyy = pbuffer.data(idx_op_geom_020_pg + 70);

    auto tr_0_0_xy_y_yyyz = pbuffer.data(idx_op_geom_020_pg + 71);

    auto tr_0_0_xy_y_yyzz = pbuffer.data(idx_op_geom_020_pg + 72);

    auto tr_0_0_xy_y_yzzz = pbuffer.data(idx_op_geom_020_pg + 73);

    auto tr_0_0_xy_y_zzzz = pbuffer.data(idx_op_geom_020_pg + 74);

    #pragma omp simd aligned(tr_0_0_xy_y_xxxx, tr_0_0_xy_y_xxxy, tr_0_0_xy_y_xxxz, tr_0_0_xy_y_xxyy, tr_0_0_xy_y_xxyz, tr_0_0_xy_y_xxzz, tr_0_0_xy_y_xyyy, tr_0_0_xy_y_xyyz, tr_0_0_xy_y_xyzz, tr_0_0_xy_y_xzzz, tr_0_0_xy_y_yyyy, tr_0_0_xy_y_yyyz, tr_0_0_xy_y_yyzz, tr_0_0_xy_y_yzzz, tr_0_0_xy_y_zzzz, tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, tr_y_xx, tr_y_xxxx, tr_y_xxxxxy, tr_y_xxxxyy, tr_y_xxxxyz, tr_y_xxxy, tr_y_xxxyyy, tr_y_xxxyyz, tr_y_xxxyzz, tr_y_xxxz, tr_y_xxyy, tr_y_xxyyyy, tr_y_xxyyyz, tr_y_xxyyzz, tr_y_xxyz, tr_y_xxyzzz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyyyy, tr_y_xyyyyz, tr_y_xyyyzz, tr_y_xyyz, tr_y_xyyzzz, tr_y_xyzz, tr_y_xyzzzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_yy_xxx, tr_yy_xxxxx, tr_yy_xxxxy, tr_yy_xxxxz, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_y_xxxx[i] = 4.0 * tr_0_xxx[i] - 2.0 * tr_0_xxxxx[i] * tke_0 - 8.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yy_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxxy[i] = 3.0 * tr_0_xxy[i] - 2.0 * tr_0_xxxxy[i] * tke_0 + 3.0 * tr_y_xx[i] - 6.0 * tr_y_xxyy[i] * tke_0 - 2.0 * tr_y_xxxx[i] * tke_0 + 4.0 * tr_y_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yy_xxy[i] * tbe_0 + 4.0 * tr_yy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxy[i] * tbe_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxxz[i] = 3.0 * tr_0_xxz[i] - 2.0 * tr_0_xxxxz[i] * tke_0 - 6.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yy_xxz[i] * tbe_0 + 4.0 * tr_yy_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxz[i] * tbe_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxyy[i] = 2.0 * tr_0_xyy[i] - 2.0 * tr_0_xxxyy[i] * tke_0 + 4.0 * tr_y_xy[i] - 4.0 * tr_y_xyyy[i] * tke_0 - 4.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yy_xyy[i] * tbe_0 + 4.0 * tr_yy_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyy[i] * tbe_0 - 4.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xxxyz[i] * tke_0 + 2.0 * tr_y_xz[i] - 4.0 * tr_y_xyyz[i] * tke_0 - 2.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yy_xyz[i] * tbe_0 + 4.0 * tr_yy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyz[i] * tbe_0 - 2.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xxzz[i] = 2.0 * tr_0_xzz[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 4.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yy_xzz[i] * tbe_0 + 4.0 * tr_yy_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxzz[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xyyy[i] = tr_0_yyy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 + 3.0 * tr_y_yy[i] - 2.0 * tr_y_yyyy[i] * tke_0 - 6.0 * tr_y_xxyy[i] * tke_0 + 4.0 * tr_y_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yy_yyy[i] * tbe_0 + 4.0 * tr_yy_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyy[i] * tbe_0 - 6.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xyyz[i] = tr_0_yyz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 + 2.0 * tr_y_yz[i] - 2.0 * tr_y_yyyz[i] * tke_0 - 4.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yy_yyz[i] * tbe_0 + 4.0 * tr_yy_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xyzz[i] = tr_0_yzz[i] - 2.0 * tr_0_xxyzz[i] * tke_0 + tr_y_zz[i] - 2.0 * tr_y_yyzz[i] * tke_0 - 2.0 * tr_y_xxzz[i] * tke_0 + 4.0 * tr_y_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yy_yzz[i] * tbe_0 + 4.0 * tr_yy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyzz[i] * tbe_0 - 2.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_xzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_y_yzzz[i] * tke_0 + 4.0 * tr_y_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yy_zzz[i] * tbe_0 + 4.0 * tr_yy_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzzz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yyyy[i] = -2.0 * tr_0_xyyyy[i] * tke_0 - 8.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_y_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyy[i] * tbe_0 - 8.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yyyz[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 6.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyz[i] * tbe_0 - 6.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yyzz[i] = -2.0 * tr_0_xyyzz[i] * tke_0 - 4.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyzz[i] * tbe_0 - 4.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_yzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_y_xzzz[i] * tke_0 + 4.0 * tr_y_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzzz[i] * tbe_0 - 2.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_y_zzzz[i] = -2.0 * tr_0_xzzzz[i] * tke_0 + 4.0 * tr_y_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yy_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 75-90 components of targeted buffer : PG

    auto tr_0_0_xy_z_xxxx = pbuffer.data(idx_op_geom_020_pg + 75);

    auto tr_0_0_xy_z_xxxy = pbuffer.data(idx_op_geom_020_pg + 76);

    auto tr_0_0_xy_z_xxxz = pbuffer.data(idx_op_geom_020_pg + 77);

    auto tr_0_0_xy_z_xxyy = pbuffer.data(idx_op_geom_020_pg + 78);

    auto tr_0_0_xy_z_xxyz = pbuffer.data(idx_op_geom_020_pg + 79);

    auto tr_0_0_xy_z_xxzz = pbuffer.data(idx_op_geom_020_pg + 80);

    auto tr_0_0_xy_z_xyyy = pbuffer.data(idx_op_geom_020_pg + 81);

    auto tr_0_0_xy_z_xyyz = pbuffer.data(idx_op_geom_020_pg + 82);

    auto tr_0_0_xy_z_xyzz = pbuffer.data(idx_op_geom_020_pg + 83);

    auto tr_0_0_xy_z_xzzz = pbuffer.data(idx_op_geom_020_pg + 84);

    auto tr_0_0_xy_z_yyyy = pbuffer.data(idx_op_geom_020_pg + 85);

    auto tr_0_0_xy_z_yyyz = pbuffer.data(idx_op_geom_020_pg + 86);

    auto tr_0_0_xy_z_yyzz = pbuffer.data(idx_op_geom_020_pg + 87);

    auto tr_0_0_xy_z_yzzz = pbuffer.data(idx_op_geom_020_pg + 88);

    auto tr_0_0_xy_z_zzzz = pbuffer.data(idx_op_geom_020_pg + 89);

    #pragma omp simd aligned(tr_0_0_xy_z_xxxx, tr_0_0_xy_z_xxxy, tr_0_0_xy_z_xxxz, tr_0_0_xy_z_xxyy, tr_0_0_xy_z_xxyz, tr_0_0_xy_z_xxzz, tr_0_0_xy_z_xyyy, tr_0_0_xy_z_xyyz, tr_0_0_xy_z_xyzz, tr_0_0_xy_z_xzzz, tr_0_0_xy_z_yyyy, tr_0_0_xy_z_yyyz, tr_0_0_xy_z_yyzz, tr_0_0_xy_z_yzzz, tr_0_0_xy_z_zzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_z_xx, tr_z_xxxx, tr_z_xxxxxy, tr_z_xxxxyy, tr_z_xxxxyz, tr_z_xxxy, tr_z_xxxyyy, tr_z_xxxyyz, tr_z_xxxyzz, tr_z_xxxz, tr_z_xxyy, tr_z_xxyyyy, tr_z_xxyyyz, tr_z_xxyyzz, tr_z_xxyz, tr_z_xxyzzz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyyyy, tr_z_xyyyyz, tr_z_xyyyzz, tr_z_xyyz, tr_z_xyyzzz, tr_z_xyzz, tr_z_xyzzzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_z_xxxx[i] = -8.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxxy[i] = 3.0 * tr_z_xx[i] - 6.0 * tr_z_xxyy[i] * tke_0 - 2.0 * tr_z_xxxx[i] * tke_0 + 4.0 * tr_z_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxxz[i] = -6.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxyy[i] = 4.0 * tr_z_xy[i] - 4.0 * tr_z_xyyy[i] * tke_0 - 4.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxyz[i] = 2.0 * tr_z_xz[i] - 4.0 * tr_z_xyyz[i] * tke_0 - 2.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xxzz[i] = -4.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xyyy[i] = 3.0 * tr_z_yy[i] - 2.0 * tr_z_yyyy[i] * tke_0 - 6.0 * tr_z_xxyy[i] * tke_0 + 4.0 * tr_z_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xyyz[i] = 2.0 * tr_z_yz[i] - 2.0 * tr_z_yyyz[i] * tke_0 - 4.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xyzz[i] = tr_z_zz[i] - 2.0 * tr_z_yyzz[i] * tke_0 - 2.0 * tr_z_xxzz[i] * tke_0 + 4.0 * tr_z_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_xzzz[i] = -2.0 * tr_z_yzzz[i] * tke_0 + 4.0 * tr_z_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yyyy[i] = -8.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_z_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yyyz[i] = -6.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yyzz[i] = -4.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_yzzz[i] = -2.0 * tr_z_xzzz[i] * tke_0 + 4.0 * tr_z_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_z_zzzz[i] = 4.0 * tr_z_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-105 components of targeted buffer : PG

    auto tr_0_0_xz_x_xxxx = pbuffer.data(idx_op_geom_020_pg + 90);

    auto tr_0_0_xz_x_xxxy = pbuffer.data(idx_op_geom_020_pg + 91);

    auto tr_0_0_xz_x_xxxz = pbuffer.data(idx_op_geom_020_pg + 92);

    auto tr_0_0_xz_x_xxyy = pbuffer.data(idx_op_geom_020_pg + 93);

    auto tr_0_0_xz_x_xxyz = pbuffer.data(idx_op_geom_020_pg + 94);

    auto tr_0_0_xz_x_xxzz = pbuffer.data(idx_op_geom_020_pg + 95);

    auto tr_0_0_xz_x_xyyy = pbuffer.data(idx_op_geom_020_pg + 96);

    auto tr_0_0_xz_x_xyyz = pbuffer.data(idx_op_geom_020_pg + 97);

    auto tr_0_0_xz_x_xyzz = pbuffer.data(idx_op_geom_020_pg + 98);

    auto tr_0_0_xz_x_xzzz = pbuffer.data(idx_op_geom_020_pg + 99);

    auto tr_0_0_xz_x_yyyy = pbuffer.data(idx_op_geom_020_pg + 100);

    auto tr_0_0_xz_x_yyyz = pbuffer.data(idx_op_geom_020_pg + 101);

    auto tr_0_0_xz_x_yyzz = pbuffer.data(idx_op_geom_020_pg + 102);

    auto tr_0_0_xz_x_yzzz = pbuffer.data(idx_op_geom_020_pg + 103);

    auto tr_0_0_xz_x_zzzz = pbuffer.data(idx_op_geom_020_pg + 104);

    #pragma omp simd aligned(tr_0_0_xz_x_xxxx, tr_0_0_xz_x_xxxy, tr_0_0_xz_x_xxxz, tr_0_0_xz_x_xxyy, tr_0_0_xz_x_xxyz, tr_0_0_xz_x_xxzz, tr_0_0_xz_x_xyyy, tr_0_0_xz_x_xyyz, tr_0_0_xz_x_xyzz, tr_0_0_xz_x_xzzz, tr_0_0_xz_x_yyyy, tr_0_0_xz_x_yyyz, tr_0_0_xz_x_yyzz, tr_0_0_xz_x_yzzz, tr_0_0_xz_x_zzzz, tr_0_xxx, tr_0_xxxxz, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_x_xx, tr_x_xxxx, tr_x_xxxxxz, tr_x_xxxxyz, tr_x_xxxxzz, tr_x_xxxy, tr_x_xxxyyz, tr_x_xxxyzz, tr_x_xxxz, tr_x_xxxzzz, tr_x_xxyy, tr_x_xxyyyz, tr_x_xxyyzz, tr_x_xxyz, tr_x_xxyzzz, tr_x_xxzz, tr_x_xxzzzz, tr_x_xy, tr_x_xyyy, tr_x_xyyyyz, tr_x_xyyyzz, tr_x_xyyz, tr_x_xyyzzz, tr_x_xyzz, tr_x_xyzzzz, tr_x_xz, tr_x_xzzz, tr_x_xzzzzz, tr_x_yy, tr_x_yyyz, tr_x_yyzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xx_xxx, tr_xx_xxxxz, tr_xx_xxxyz, tr_xx_xxxzz, tr_xx_xxy, tr_xx_xxyyz, tr_xx_xxyzz, tr_xx_xxz, tr_xx_xxzzz, tr_xx_xyy, tr_xx_xyyyz, tr_xx_xyyzz, tr_xx_xyz, tr_xx_xyzzz, tr_xx_xzz, tr_xx_xzzzz, tr_xx_yyy, tr_xx_yyyyz, tr_xx_yyyzz, tr_xx_yyz, tr_xx_yyzzz, tr_xx_yzz, tr_xx_yzzzz, tr_xx_zzz, tr_xx_zzzzz, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xzzz, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yzzz, tr_xxz_zzzz, tr_xz_xxx, tr_xz_xxxxx, tr_xz_xxxxy, tr_xz_xxxxz, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_x_xxxx[i] = -2.0 * tr_0_xxxxz[i] * tke_0 - 2.0 * tr_z_xxxx[i] * tbe_0 - 8.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxxy[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_z_xxxy[i] * tbe_0 - 6.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxxz[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 2.0 * tr_z_xxxz[i] * tbe_0 + 3.0 * tr_x_xx[i] - 6.0 * tr_x_xxzz[i] * tke_0 - 2.0 * tr_x_xxxx[i] * tke_0 + 4.0 * tr_x_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxx[i] * tbe_0 + 4.0 * tr_xx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxyy[i] = -2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_z_xxyy[i] * tbe_0 - 4.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxyz[i] = tr_0_xxy[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_z_xxyz[i] * tbe_0 + 2.0 * tr_x_xy[i] - 4.0 * tr_x_xyzz[i] * tke_0 - 2.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xxy[i] * tbe_0 + 4.0 * tr_xx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xxzz[i] = 2.0 * tr_0_xxz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_z_xxzz[i] * tbe_0 + 4.0 * tr_x_xz[i] - 4.0 * tr_x_xzzz[i] * tke_0 - 4.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_xxz[i] * tbe_0 + 4.0 * tr_xx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xyyy[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_z_xyyy[i] * tbe_0 - 2.0 * tr_x_yyyz[i] * tke_0 + 4.0 * tr_x_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xyyz[i] = tr_0_xyy[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_z_xyyz[i] * tbe_0 + tr_x_yy[i] - 2.0 * tr_x_yyzz[i] * tke_0 - 2.0 * tr_x_xxyy[i] * tke_0 + 4.0 * tr_x_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_xyy[i] * tbe_0 + 4.0 * tr_xx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xyzz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_z_xyzz[i] * tbe_0 + 2.0 * tr_x_yz[i] - 2.0 * tr_x_yzzz[i] * tke_0 - 4.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_xyz[i] * tbe_0 + 4.0 * tr_xx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_xzzz[i] = 3.0 * tr_0_xzz[i] - 2.0 * tr_0_xzzzz[i] * tke_0 - 2.0 * tr_z_xzzz[i] * tbe_0 + 3.0 * tr_x_zz[i] - 2.0 * tr_x_zzzz[i] * tke_0 - 6.0 * tr_x_xxzz[i] * tke_0 + 4.0 * tr_x_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xx_xzz[i] * tbe_0 + 4.0 * tr_xx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yyyy[i] = -2.0 * tr_0_yyyyz[i] * tke_0 - 2.0 * tr_z_yyyy[i] * tbe_0 + 4.0 * tr_x_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yyyz[i] = tr_0_yyy[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_z_yyyz[i] * tbe_0 - 2.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_x_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xx_yyy[i] * tbe_0 + 4.0 * tr_xx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yyzz[i] = 2.0 * tr_0_yyz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_z_yyzz[i] * tbe_0 - 4.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xx_yyz[i] * tbe_0 + 4.0 * tr_xx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_yzzz[i] = 3.0 * tr_0_yzz[i] - 2.0 * tr_0_yzzzz[i] * tke_0 - 2.0 * tr_z_yzzz[i] * tbe_0 - 6.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xx_yzz[i] * tbe_0 + 4.0 * tr_xx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_x_zzzz[i] = 4.0 * tr_0_zzz[i] - 2.0 * tr_0_zzzzz[i] * tke_0 - 2.0 * tr_z_zzzz[i] * tbe_0 - 8.0 * tr_x_xzzz[i] * tke_0 + 4.0 * tr_x_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xx_zzz[i] * tbe_0 + 4.0 * tr_xx_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 105-120 components of targeted buffer : PG

    auto tr_0_0_xz_y_xxxx = pbuffer.data(idx_op_geom_020_pg + 105);

    auto tr_0_0_xz_y_xxxy = pbuffer.data(idx_op_geom_020_pg + 106);

    auto tr_0_0_xz_y_xxxz = pbuffer.data(idx_op_geom_020_pg + 107);

    auto tr_0_0_xz_y_xxyy = pbuffer.data(idx_op_geom_020_pg + 108);

    auto tr_0_0_xz_y_xxyz = pbuffer.data(idx_op_geom_020_pg + 109);

    auto tr_0_0_xz_y_xxzz = pbuffer.data(idx_op_geom_020_pg + 110);

    auto tr_0_0_xz_y_xyyy = pbuffer.data(idx_op_geom_020_pg + 111);

    auto tr_0_0_xz_y_xyyz = pbuffer.data(idx_op_geom_020_pg + 112);

    auto tr_0_0_xz_y_xyzz = pbuffer.data(idx_op_geom_020_pg + 113);

    auto tr_0_0_xz_y_xzzz = pbuffer.data(idx_op_geom_020_pg + 114);

    auto tr_0_0_xz_y_yyyy = pbuffer.data(idx_op_geom_020_pg + 115);

    auto tr_0_0_xz_y_yyyz = pbuffer.data(idx_op_geom_020_pg + 116);

    auto tr_0_0_xz_y_yyzz = pbuffer.data(idx_op_geom_020_pg + 117);

    auto tr_0_0_xz_y_yzzz = pbuffer.data(idx_op_geom_020_pg + 118);

    auto tr_0_0_xz_y_zzzz = pbuffer.data(idx_op_geom_020_pg + 119);

    #pragma omp simd aligned(tr_0_0_xz_y_xxxx, tr_0_0_xz_y_xxxy, tr_0_0_xz_y_xxxz, tr_0_0_xz_y_xxyy, tr_0_0_xz_y_xxyz, tr_0_0_xz_y_xxzz, tr_0_0_xz_y_xyyy, tr_0_0_xz_y_xyyz, tr_0_0_xz_y_xyzz, tr_0_0_xz_y_xzzz, tr_0_0_xz_y_yyyy, tr_0_0_xz_y_yyyz, tr_0_0_xz_y_yyzz, tr_0_0_xz_y_yzzz, tr_0_0_xz_y_zzzz, tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_y_xx, tr_y_xxxx, tr_y_xxxxxz, tr_y_xxxxyz, tr_y_xxxxzz, tr_y_xxxy, tr_y_xxxyyz, tr_y_xxxyzz, tr_y_xxxz, tr_y_xxxzzz, tr_y_xxyy, tr_y_xxyyyz, tr_y_xxyyzz, tr_y_xxyz, tr_y_xxyzzz, tr_y_xxzz, tr_y_xxzzzz, tr_y_xy, tr_y_xyyy, tr_y_xyyyyz, tr_y_xyyyzz, tr_y_xyyz, tr_y_xyyzzz, tr_y_xyzz, tr_y_xyzzzz, tr_y_xz, tr_y_xzzz, tr_y_xzzzzz, tr_y_yy, tr_y_yyyz, tr_y_yyzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxx, tr_yz_xxxxy, tr_yz_xxxxz, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_y_xxxx[i] = -8.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxxy[i] = -6.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxxz[i] = 3.0 * tr_y_xx[i] - 6.0 * tr_y_xxzz[i] * tke_0 - 2.0 * tr_y_xxxx[i] * tke_0 + 4.0 * tr_y_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxyy[i] = -4.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxyz[i] = 2.0 * tr_y_xy[i] - 4.0 * tr_y_xyzz[i] * tke_0 - 2.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xxzz[i] = 4.0 * tr_y_xz[i] - 4.0 * tr_y_xzzz[i] * tke_0 - 4.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xyyy[i] = -2.0 * tr_y_yyyz[i] * tke_0 + 4.0 * tr_y_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xyyz[i] = tr_y_yy[i] - 2.0 * tr_y_yyzz[i] * tke_0 - 2.0 * tr_y_xxyy[i] * tke_0 + 4.0 * tr_y_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xyzz[i] = 2.0 * tr_y_yz[i] - 2.0 * tr_y_yzzz[i] * tke_0 - 4.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_xzzz[i] = 3.0 * tr_y_zz[i] - 2.0 * tr_y_zzzz[i] * tke_0 - 6.0 * tr_y_xxzz[i] * tke_0 + 4.0 * tr_y_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yyyy[i] = 4.0 * tr_y_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yyyz[i] = -2.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_y_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yyzz[i] = -4.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_yzzz[i] = -6.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_y_zzzz[i] = -8.0 * tr_y_xzzz[i] * tke_0 + 4.0 * tr_y_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-135 components of targeted buffer : PG

    auto tr_0_0_xz_z_xxxx = pbuffer.data(idx_op_geom_020_pg + 120);

    auto tr_0_0_xz_z_xxxy = pbuffer.data(idx_op_geom_020_pg + 121);

    auto tr_0_0_xz_z_xxxz = pbuffer.data(idx_op_geom_020_pg + 122);

    auto tr_0_0_xz_z_xxyy = pbuffer.data(idx_op_geom_020_pg + 123);

    auto tr_0_0_xz_z_xxyz = pbuffer.data(idx_op_geom_020_pg + 124);

    auto tr_0_0_xz_z_xxzz = pbuffer.data(idx_op_geom_020_pg + 125);

    auto tr_0_0_xz_z_xyyy = pbuffer.data(idx_op_geom_020_pg + 126);

    auto tr_0_0_xz_z_xyyz = pbuffer.data(idx_op_geom_020_pg + 127);

    auto tr_0_0_xz_z_xyzz = pbuffer.data(idx_op_geom_020_pg + 128);

    auto tr_0_0_xz_z_xzzz = pbuffer.data(idx_op_geom_020_pg + 129);

    auto tr_0_0_xz_z_yyyy = pbuffer.data(idx_op_geom_020_pg + 130);

    auto tr_0_0_xz_z_yyyz = pbuffer.data(idx_op_geom_020_pg + 131);

    auto tr_0_0_xz_z_yyzz = pbuffer.data(idx_op_geom_020_pg + 132);

    auto tr_0_0_xz_z_yzzz = pbuffer.data(idx_op_geom_020_pg + 133);

    auto tr_0_0_xz_z_zzzz = pbuffer.data(idx_op_geom_020_pg + 134);

    #pragma omp simd aligned(tr_0_0_xz_z_xxxx, tr_0_0_xz_z_xxxy, tr_0_0_xz_z_xxxz, tr_0_0_xz_z_xxyy, tr_0_0_xz_z_xxyz, tr_0_0_xz_z_xxzz, tr_0_0_xz_z_xyyy, tr_0_0_xz_z_xyyz, tr_0_0_xz_z_xyzz, tr_0_0_xz_z_xzzz, tr_0_0_xz_z_yyyy, tr_0_0_xz_z_yyyz, tr_0_0_xz_z_yyzz, tr_0_0_xz_z_yzzz, tr_0_0_xz_z_zzzz, tr_0_xxx, tr_0_xxxxx, tr_0_xxxxy, tr_0_xxxxz, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyz, tr_0_yzz, tr_0_zzz, tr_x_xxxx, tr_x_xxxy, tr_x_xxxz, tr_x_xxyy, tr_x_xxyz, tr_x_xxzz, tr_x_xyyy, tr_x_xyyz, tr_x_xyzz, tr_x_xzzz, tr_x_yyyy, tr_x_yyyz, tr_x_yyzz, tr_x_yzzz, tr_x_zzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, tr_z_xx, tr_z_xxxx, tr_z_xxxxxz, tr_z_xxxxyz, tr_z_xxxxzz, tr_z_xxxy, tr_z_xxxyyz, tr_z_xxxyzz, tr_z_xxxz, tr_z_xxxzzz, tr_z_xxyy, tr_z_xxyyyz, tr_z_xxyyzz, tr_z_xxyz, tr_z_xxyzzz, tr_z_xxzz, tr_z_xxzzzz, tr_z_xy, tr_z_xyyy, tr_z_xyyyyz, tr_z_xyyyzz, tr_z_xyyz, tr_z_xyyzzz, tr_z_xyzz, tr_z_xyzzzz, tr_z_xz, tr_z_xzzz, tr_z_xzzzzz, tr_z_yy, tr_z_yyyz, tr_z_yyzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxx, tr_zz_xxxxy, tr_zz_xxxxz, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_z_xxxx[i] = 4.0 * tr_0_xxx[i] - 2.0 * tr_0_xxxxx[i] * tke_0 - 8.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_zz_xxx[i] * tbe_0 + 4.0 * tr_zz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxx[i] * tbe_0 + 4.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxxy[i] = 3.0 * tr_0_xxy[i] - 2.0 * tr_0_xxxxy[i] * tke_0 - 6.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_zz_xxy[i] * tbe_0 + 4.0 * tr_zz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxy[i] * tbe_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxxz[i] = 3.0 * tr_0_xxz[i] - 2.0 * tr_0_xxxxz[i] * tke_0 + 3.0 * tr_z_xx[i] - 6.0 * tr_z_xxzz[i] * tke_0 - 2.0 * tr_z_xxxx[i] * tke_0 + 4.0 * tr_z_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_zz_xxz[i] * tbe_0 + 4.0 * tr_zz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxxz[i] * tbe_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxyy[i] = 2.0 * tr_0_xyy[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 4.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xyy[i] * tbe_0 + 4.0 * tr_zz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyy[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xxxyz[i] * tke_0 + 2.0 * tr_z_xy[i] - 4.0 * tr_z_xyzz[i] * tke_0 - 2.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xyz[i] * tbe_0 + 4.0 * tr_zz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxyz[i] * tbe_0 - 2.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xxzz[i] = 2.0 * tr_0_xzz[i] - 2.0 * tr_0_xxxzz[i] * tke_0 + 4.0 * tr_z_xz[i] - 4.0 * tr_z_xzzz[i] * tke_0 - 4.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xzz[i] * tbe_0 + 4.0 * tr_zz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxzz[i] * tbe_0 - 4.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xyyy[i] = tr_0_yyy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 2.0 * tr_z_yyyz[i] * tke_0 + 4.0 * tr_z_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_zz_yyy[i] * tbe_0 + 4.0 * tr_zz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyy[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xyyz[i] = tr_0_yyz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 + tr_z_yy[i] - 2.0 * tr_z_yyzz[i] * tke_0 - 2.0 * tr_z_xxyy[i] * tke_0 + 4.0 * tr_z_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_yyz[i] * tbe_0 + 4.0 * tr_zz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyyz[i] * tbe_0 - 2.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xyzz[i] = tr_0_yzz[i] - 2.0 * tr_0_xxyzz[i] * tke_0 + 2.0 * tr_z_yz[i] - 2.0 * tr_z_yzzz[i] * tke_0 - 4.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_yzz[i] * tbe_0 + 4.0 * tr_zz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyzz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_xzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 + 3.0 * tr_z_zz[i] - 2.0 * tr_z_zzzz[i] * tke_0 - 6.0 * tr_z_xxzz[i] * tke_0 + 4.0 * tr_z_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_zzz[i] * tbe_0 + 4.0 * tr_zz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzzz[i] * tbe_0 - 6.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yyyy[i] = -2.0 * tr_0_xyyyy[i] * tke_0 + 4.0 * tr_z_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyy[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yyyz[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_z_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyyz[i] * tbe_0 - 2.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yyzz[i] = -2.0 * tr_0_xyyzz[i] * tke_0 - 4.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyzz[i] * tbe_0 - 4.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_yzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yzzz[i] * tbe_0 - 6.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_z_zzzz[i] = -2.0 * tr_0_xzzzz[i] * tke_0 - 8.0 * tr_z_xzzz[i] * tke_0 + 4.0 * tr_z_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzzz[i] * tbe_0 - 8.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 135-150 components of targeted buffer : PG

    auto tr_0_0_yy_x_xxxx = pbuffer.data(idx_op_geom_020_pg + 135);

    auto tr_0_0_yy_x_xxxy = pbuffer.data(idx_op_geom_020_pg + 136);

    auto tr_0_0_yy_x_xxxz = pbuffer.data(idx_op_geom_020_pg + 137);

    auto tr_0_0_yy_x_xxyy = pbuffer.data(idx_op_geom_020_pg + 138);

    auto tr_0_0_yy_x_xxyz = pbuffer.data(idx_op_geom_020_pg + 139);

    auto tr_0_0_yy_x_xxzz = pbuffer.data(idx_op_geom_020_pg + 140);

    auto tr_0_0_yy_x_xyyy = pbuffer.data(idx_op_geom_020_pg + 141);

    auto tr_0_0_yy_x_xyyz = pbuffer.data(idx_op_geom_020_pg + 142);

    auto tr_0_0_yy_x_xyzz = pbuffer.data(idx_op_geom_020_pg + 143);

    auto tr_0_0_yy_x_xzzz = pbuffer.data(idx_op_geom_020_pg + 144);

    auto tr_0_0_yy_x_yyyy = pbuffer.data(idx_op_geom_020_pg + 145);

    auto tr_0_0_yy_x_yyyz = pbuffer.data(idx_op_geom_020_pg + 146);

    auto tr_0_0_yy_x_yyzz = pbuffer.data(idx_op_geom_020_pg + 147);

    auto tr_0_0_yy_x_yzzz = pbuffer.data(idx_op_geom_020_pg + 148);

    auto tr_0_0_yy_x_zzzz = pbuffer.data(idx_op_geom_020_pg + 149);

    #pragma omp simd aligned(tr_0_0_yy_x_xxxx, tr_0_0_yy_x_xxxy, tr_0_0_yy_x_xxxz, tr_0_0_yy_x_xxyy, tr_0_0_yy_x_xxyz, tr_0_0_yy_x_xxzz, tr_0_0_yy_x_xyyy, tr_0_0_yy_x_xyyz, tr_0_0_yy_x_xyzz, tr_0_0_yy_x_xzzz, tr_0_0_yy_x_yyyy, tr_0_0_yy_x_yyyz, tr_0_0_yy_x_yyzz, tr_0_0_yy_x_yzzz, tr_0_0_yy_x_zzzz, tr_x_xx, tr_x_xxxx, tr_x_xxxxyy, tr_x_xxxy, tr_x_xxxyyy, tr_x_xxxyyz, tr_x_xxxz, tr_x_xxyy, tr_x_xxyyyy, tr_x_xxyyyz, tr_x_xxyyzz, tr_x_xxyz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyyyy, tr_x_xyyyyz, tr_x_xyyyzz, tr_x_xyyz, tr_x_xyyzzz, tr_x_xyzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyy, tr_x_yyyyyy, tr_x_yyyyyz, tr_x_yyyyzz, tr_x_yyyz, tr_x_yyyzzz, tr_x_yyzz, tr_x_yyzzzz, tr_x_yz, tr_x_yzzz, tr_x_zz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxy, tr_xy_xxxyy, tr_xy_xxxyz, tr_xy_xxy, tr_xy_xxyyy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xyy, tr_xy_xyyyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyyyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xzzz, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yzzz, tr_xyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_x_xxxx[i] = -2.0 * tr_x_xxxx[i] * tbe_0 - 2.0 * tr_x_xxxx[i] * tke_0 + 4.0 * tr_x_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxxy[i] = -2.0 * tr_x_xxxy[i] * tbe_0 - 6.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xy_xxx[i] * tbe_0 + 8.0 * tr_xy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxxz[i] = -2.0 * tr_x_xxxz[i] * tbe_0 - 2.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxyy[i] = 2.0 * tr_x_xx[i] - 2.0 * tr_x_xxyy[i] * tbe_0 - 10.0 * tr_x_xxyy[i] * tke_0 + 4.0 * tr_x_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xy_xxy[i] * tbe_0 + 8.0 * tr_xy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxyz[i] = -2.0 * tr_x_xxyz[i] * tbe_0 - 6.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xy_xxz[i] * tbe_0 + 8.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xxzz[i] = -2.0 * tr_x_xxzz[i] * tbe_0 - 2.0 * tr_x_xxzz[i] * tke_0 + 4.0 * tr_x_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xyyy[i] = 6.0 * tr_x_xy[i] - 2.0 * tr_x_xyyy[i] * tbe_0 - 14.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_x_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xy_xyy[i] * tbe_0 + 8.0 * tr_xy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xyyz[i] = 2.0 * tr_x_xz[i] - 2.0 * tr_x_xyyz[i] * tbe_0 - 10.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xy_xyz[i] * tbe_0 + 8.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xyzz[i] = -2.0 * tr_x_xyzz[i] * tbe_0 - 6.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_xzz[i] * tbe_0 + 8.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_xzzz[i] = -2.0 * tr_x_xzzz[i] * tbe_0 - 2.0 * tr_x_xzzz[i] * tke_0 + 4.0 * tr_x_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yyyy[i] = 12.0 * tr_x_yy[i] - 2.0 * tr_x_yyyy[i] * tbe_0 - 18.0 * tr_x_yyyy[i] * tke_0 + 4.0 * tr_x_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xy_yyy[i] * tbe_0 + 8.0 * tr_xy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yyyz[i] = 6.0 * tr_x_yz[i] - 2.0 * tr_x_yyyz[i] * tbe_0 - 14.0 * tr_x_yyyz[i] * tke_0 + 4.0 * tr_x_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xy_yyz[i] * tbe_0 + 8.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yyzz[i] = 2.0 * tr_x_zz[i] - 2.0 * tr_x_yyzz[i] * tbe_0 - 10.0 * tr_x_yyzz[i] * tke_0 + 4.0 * tr_x_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xy_yzz[i] * tbe_0 + 8.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_yzzz[i] = -2.0 * tr_x_yzzz[i] * tbe_0 - 6.0 * tr_x_yzzz[i] * tke_0 + 4.0 * tr_x_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xy_zzz[i] * tbe_0 + 8.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_x_zzzz[i] = -2.0 * tr_x_zzzz[i] * tbe_0 - 2.0 * tr_x_zzzz[i] * tke_0 + 4.0 * tr_x_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-165 components of targeted buffer : PG

    auto tr_0_0_yy_y_xxxx = pbuffer.data(idx_op_geom_020_pg + 150);

    auto tr_0_0_yy_y_xxxy = pbuffer.data(idx_op_geom_020_pg + 151);

    auto tr_0_0_yy_y_xxxz = pbuffer.data(idx_op_geom_020_pg + 152);

    auto tr_0_0_yy_y_xxyy = pbuffer.data(idx_op_geom_020_pg + 153);

    auto tr_0_0_yy_y_xxyz = pbuffer.data(idx_op_geom_020_pg + 154);

    auto tr_0_0_yy_y_xxzz = pbuffer.data(idx_op_geom_020_pg + 155);

    auto tr_0_0_yy_y_xyyy = pbuffer.data(idx_op_geom_020_pg + 156);

    auto tr_0_0_yy_y_xyyz = pbuffer.data(idx_op_geom_020_pg + 157);

    auto tr_0_0_yy_y_xyzz = pbuffer.data(idx_op_geom_020_pg + 158);

    auto tr_0_0_yy_y_xzzz = pbuffer.data(idx_op_geom_020_pg + 159);

    auto tr_0_0_yy_y_yyyy = pbuffer.data(idx_op_geom_020_pg + 160);

    auto tr_0_0_yy_y_yyyz = pbuffer.data(idx_op_geom_020_pg + 161);

    auto tr_0_0_yy_y_yyzz = pbuffer.data(idx_op_geom_020_pg + 162);

    auto tr_0_0_yy_y_yzzz = pbuffer.data(idx_op_geom_020_pg + 163);

    auto tr_0_0_yy_y_zzzz = pbuffer.data(idx_op_geom_020_pg + 164);

    #pragma omp simd aligned(tr_0_0_yy_y_xxxx, tr_0_0_yy_y_xxxy, tr_0_0_yy_y_xxxz, tr_0_0_yy_y_xxyy, tr_0_0_yy_y_xxyz, tr_0_0_yy_y_xxzz, tr_0_0_yy_y_xyyy, tr_0_0_yy_y_xyyz, tr_0_0_yy_y_xyzz, tr_0_0_yy_y_xzzz, tr_0_0_yy_y_yyyy, tr_0_0_yy_y_yyyz, tr_0_0_yy_y_yyzz, tr_0_0_yy_y_yzzz, tr_0_0_yy_y_zzzz, tr_0_xxx, tr_0_xxxxy, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_y_xx, tr_y_xxxx, tr_y_xxxxyy, tr_y_xxxy, tr_y_xxxyyy, tr_y_xxxyyz, tr_y_xxxz, tr_y_xxyy, tr_y_xxyyyy, tr_y_xxyyyz, tr_y_xxyyzz, tr_y_xxyz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyyyy, tr_y_xyyyyz, tr_y_xyyyzz, tr_y_xyyz, tr_y_xyyzzz, tr_y_xyzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyy, tr_y_yyyyyy, tr_y_yyyyyz, tr_y_yyyyzz, tr_y_yyyz, tr_y_yyyzzz, tr_y_yyzz, tr_y_yyzzzz, tr_y_yz, tr_y_yzzz, tr_y_zz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxy, tr_yy_xxxyy, tr_yy_xxxyz, tr_yy_xxy, tr_yy_xxyyy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xyy, tr_yy_xyyyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyyyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xzzz, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yzzz, tr_yyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_y_xxxx[i] = -4.0 * tr_0_xxxxy[i] * tke_0 - 6.0 * tr_y_xxxx[i] * tbe_0 - 2.0 * tr_y_xxxx[i] * tke_0 + 4.0 * tr_y_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxxy[i] = 2.0 * tr_0_xxx[i] - 4.0 * tr_0_xxxyy[i] * tke_0 - 6.0 * tr_y_xxxy[i] * tbe_0 - 6.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yy_xxx[i] * tbe_0 + 8.0 * tr_yy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxxz[i] = -4.0 * tr_0_xxxyz[i] * tke_0 - 6.0 * tr_y_xxxz[i] * tbe_0 - 2.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxyy[i] = 4.0 * tr_0_xxy[i] - 4.0 * tr_0_xxyyy[i] * tke_0 + 2.0 * tr_y_xx[i] - 6.0 * tr_y_xxyy[i] * tbe_0 - 10.0 * tr_y_xxyy[i] * tke_0 + 4.0 * tr_y_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yy_xxy[i] * tbe_0 + 8.0 * tr_yy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxyz[i] = 2.0 * tr_0_xxz[i] - 4.0 * tr_0_xxyyz[i] * tke_0 - 6.0 * tr_y_xxyz[i] * tbe_0 - 6.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yy_xxz[i] * tbe_0 + 8.0 * tr_yy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xxzz[i] = -4.0 * tr_0_xxyzz[i] * tke_0 - 6.0 * tr_y_xxzz[i] * tbe_0 - 2.0 * tr_y_xxzz[i] * tke_0 + 4.0 * tr_y_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xyyy[i] = 6.0 * tr_0_xyy[i] - 4.0 * tr_0_xyyyy[i] * tke_0 + 6.0 * tr_y_xy[i] - 6.0 * tr_y_xyyy[i] * tbe_0 - 14.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_y_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yy_xyy[i] * tbe_0 + 8.0 * tr_yy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xyyz[i] = 4.0 * tr_0_xyz[i] - 4.0 * tr_0_xyyyz[i] * tke_0 + 2.0 * tr_y_xz[i] - 6.0 * tr_y_xyyz[i] * tbe_0 - 10.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yy_xyz[i] * tbe_0 + 8.0 * tr_yy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xyzz[i] = 2.0 * tr_0_xzz[i] - 4.0 * tr_0_xyyzz[i] * tke_0 - 6.0 * tr_y_xyzz[i] * tbe_0 - 6.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yy_xzz[i] * tbe_0 + 8.0 * tr_yy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_xzzz[i] = -4.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_y_xzzz[i] * tbe_0 - 2.0 * tr_y_xzzz[i] * tke_0 + 4.0 * tr_y_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yyyy[i] = 8.0 * tr_0_yyy[i] - 4.0 * tr_0_yyyyy[i] * tke_0 + 12.0 * tr_y_yy[i] - 6.0 * tr_y_yyyy[i] * tbe_0 - 18.0 * tr_y_yyyy[i] * tke_0 + 4.0 * tr_y_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yy_yyy[i] * tbe_0 + 8.0 * tr_yy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yyyz[i] = 6.0 * tr_0_yyz[i] - 4.0 * tr_0_yyyyz[i] * tke_0 + 6.0 * tr_y_yz[i] - 6.0 * tr_y_yyyz[i] * tbe_0 - 14.0 * tr_y_yyyz[i] * tke_0 + 4.0 * tr_y_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yy_yyz[i] * tbe_0 + 8.0 * tr_yy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yyzz[i] = 4.0 * tr_0_yzz[i] - 4.0 * tr_0_yyyzz[i] * tke_0 + 2.0 * tr_y_zz[i] - 6.0 * tr_y_yyzz[i] * tbe_0 - 10.0 * tr_y_yyzz[i] * tke_0 + 4.0 * tr_y_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yy_yzz[i] * tbe_0 + 8.0 * tr_yy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_yzzz[i] = 2.0 * tr_0_zzz[i] - 4.0 * tr_0_yyzzz[i] * tke_0 - 6.0 * tr_y_yzzz[i] * tbe_0 - 6.0 * tr_y_yzzz[i] * tke_0 + 4.0 * tr_y_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yy_zzz[i] * tbe_0 + 8.0 * tr_yy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_y_zzzz[i] = -4.0 * tr_0_yzzzz[i] * tke_0 - 6.0 * tr_y_zzzz[i] * tbe_0 - 2.0 * tr_y_zzzz[i] * tke_0 + 4.0 * tr_y_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 165-180 components of targeted buffer : PG

    auto tr_0_0_yy_z_xxxx = pbuffer.data(idx_op_geom_020_pg + 165);

    auto tr_0_0_yy_z_xxxy = pbuffer.data(idx_op_geom_020_pg + 166);

    auto tr_0_0_yy_z_xxxz = pbuffer.data(idx_op_geom_020_pg + 167);

    auto tr_0_0_yy_z_xxyy = pbuffer.data(idx_op_geom_020_pg + 168);

    auto tr_0_0_yy_z_xxyz = pbuffer.data(idx_op_geom_020_pg + 169);

    auto tr_0_0_yy_z_xxzz = pbuffer.data(idx_op_geom_020_pg + 170);

    auto tr_0_0_yy_z_xyyy = pbuffer.data(idx_op_geom_020_pg + 171);

    auto tr_0_0_yy_z_xyyz = pbuffer.data(idx_op_geom_020_pg + 172);

    auto tr_0_0_yy_z_xyzz = pbuffer.data(idx_op_geom_020_pg + 173);

    auto tr_0_0_yy_z_xzzz = pbuffer.data(idx_op_geom_020_pg + 174);

    auto tr_0_0_yy_z_yyyy = pbuffer.data(idx_op_geom_020_pg + 175);

    auto tr_0_0_yy_z_yyyz = pbuffer.data(idx_op_geom_020_pg + 176);

    auto tr_0_0_yy_z_yyzz = pbuffer.data(idx_op_geom_020_pg + 177);

    auto tr_0_0_yy_z_yzzz = pbuffer.data(idx_op_geom_020_pg + 178);

    auto tr_0_0_yy_z_zzzz = pbuffer.data(idx_op_geom_020_pg + 179);

    #pragma omp simd aligned(tr_0_0_yy_z_xxxx, tr_0_0_yy_z_xxxy, tr_0_0_yy_z_xxxz, tr_0_0_yy_z_xxyy, tr_0_0_yy_z_xxyz, tr_0_0_yy_z_xxzz, tr_0_0_yy_z_xyyy, tr_0_0_yy_z_xyyz, tr_0_0_yy_z_xyzz, tr_0_0_yy_z_xzzz, tr_0_0_yy_z_yyyy, tr_0_0_yy_z_yyyz, tr_0_0_yy_z_yyzz, tr_0_0_yy_z_yzzz, tr_0_0_yy_z_zzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_xx, tr_z_xxxx, tr_z_xxxxyy, tr_z_xxxy, tr_z_xxxyyy, tr_z_xxxyyz, tr_z_xxxz, tr_z_xxyy, tr_z_xxyyyy, tr_z_xxyyyz, tr_z_xxyyzz, tr_z_xxyz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyyyy, tr_z_xyyyyz, tr_z_xyyyzz, tr_z_xyyz, tr_z_xyyzzz, tr_z_xyzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyy, tr_z_yyyyyy, tr_z_yyyyyz, tr_z_yyyyzz, tr_z_yyyz, tr_z_yyyzzz, tr_z_yyzz, tr_z_yyzzzz, tr_z_yz, tr_z_yzzz, tr_z_zz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_z_xxxx[i] = -2.0 * tr_z_xxxx[i] * tbe_0 - 2.0 * tr_z_xxxx[i] * tke_0 + 4.0 * tr_z_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxxy[i] = -2.0 * tr_z_xxxy[i] * tbe_0 - 6.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yz_xxx[i] * tbe_0 + 8.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxxz[i] = -2.0 * tr_z_xxxz[i] * tbe_0 - 2.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxyy[i] = 2.0 * tr_z_xx[i] - 2.0 * tr_z_xxyy[i] * tbe_0 - 10.0 * tr_z_xxyy[i] * tke_0 + 4.0 * tr_z_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yz_xxy[i] * tbe_0 + 8.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxyz[i] = -2.0 * tr_z_xxyz[i] * tbe_0 - 6.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xxz[i] * tbe_0 + 8.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xxzz[i] = -2.0 * tr_z_xxzz[i] * tbe_0 - 2.0 * tr_z_xxzz[i] * tke_0 + 4.0 * tr_z_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xyyy[i] = 6.0 * tr_z_xy[i] - 2.0 * tr_z_xyyy[i] * tbe_0 - 14.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_z_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yz_xyy[i] * tbe_0 + 8.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xyyz[i] = 2.0 * tr_z_xz[i] - 2.0 * tr_z_xyyz[i] * tbe_0 - 10.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yz_xyz[i] * tbe_0 + 8.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xyzz[i] = -2.0 * tr_z_xyzz[i] * tbe_0 - 6.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xzz[i] * tbe_0 + 8.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_xzzz[i] = -2.0 * tr_z_xzzz[i] * tbe_0 - 2.0 * tr_z_xzzz[i] * tke_0 + 4.0 * tr_z_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yyyy[i] = 12.0 * tr_z_yy[i] - 2.0 * tr_z_yyyy[i] * tbe_0 - 18.0 * tr_z_yyyy[i] * tke_0 + 4.0 * tr_z_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yz_yyy[i] * tbe_0 + 8.0 * tr_yz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yyyz[i] = 6.0 * tr_z_yz[i] - 2.0 * tr_z_yyyz[i] * tbe_0 - 14.0 * tr_z_yyyz[i] * tke_0 + 4.0 * tr_z_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yz_yyz[i] * tbe_0 + 8.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yyzz[i] = 2.0 * tr_z_zz[i] - 2.0 * tr_z_yyzz[i] * tbe_0 - 10.0 * tr_z_yyzz[i] * tke_0 + 4.0 * tr_z_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yz_yzz[i] * tbe_0 + 8.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_yzzz[i] = -2.0 * tr_z_yzzz[i] * tbe_0 - 6.0 * tr_z_yzzz[i] * tke_0 + 4.0 * tr_z_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_zzz[i] * tbe_0 + 8.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_z_zzzz[i] = -2.0 * tr_z_zzzz[i] * tbe_0 - 2.0 * tr_z_zzzz[i] * tke_0 + 4.0 * tr_z_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-195 components of targeted buffer : PG

    auto tr_0_0_yz_x_xxxx = pbuffer.data(idx_op_geom_020_pg + 180);

    auto tr_0_0_yz_x_xxxy = pbuffer.data(idx_op_geom_020_pg + 181);

    auto tr_0_0_yz_x_xxxz = pbuffer.data(idx_op_geom_020_pg + 182);

    auto tr_0_0_yz_x_xxyy = pbuffer.data(idx_op_geom_020_pg + 183);

    auto tr_0_0_yz_x_xxyz = pbuffer.data(idx_op_geom_020_pg + 184);

    auto tr_0_0_yz_x_xxzz = pbuffer.data(idx_op_geom_020_pg + 185);

    auto tr_0_0_yz_x_xyyy = pbuffer.data(idx_op_geom_020_pg + 186);

    auto tr_0_0_yz_x_xyyz = pbuffer.data(idx_op_geom_020_pg + 187);

    auto tr_0_0_yz_x_xyzz = pbuffer.data(idx_op_geom_020_pg + 188);

    auto tr_0_0_yz_x_xzzz = pbuffer.data(idx_op_geom_020_pg + 189);

    auto tr_0_0_yz_x_yyyy = pbuffer.data(idx_op_geom_020_pg + 190);

    auto tr_0_0_yz_x_yyyz = pbuffer.data(idx_op_geom_020_pg + 191);

    auto tr_0_0_yz_x_yyzz = pbuffer.data(idx_op_geom_020_pg + 192);

    auto tr_0_0_yz_x_yzzz = pbuffer.data(idx_op_geom_020_pg + 193);

    auto tr_0_0_yz_x_zzzz = pbuffer.data(idx_op_geom_020_pg + 194);

    #pragma omp simd aligned(tr_0_0_yz_x_xxxx, tr_0_0_yz_x_xxxy, tr_0_0_yz_x_xxxz, tr_0_0_yz_x_xxyy, tr_0_0_yz_x_xxyz, tr_0_0_yz_x_xxzz, tr_0_0_yz_x_xyyy, tr_0_0_yz_x_xyyz, tr_0_0_yz_x_xyzz, tr_0_0_yz_x_xzzz, tr_0_0_yz_x_yyyy, tr_0_0_yz_x_yyyz, tr_0_0_yz_x_yyzz, tr_0_0_yz_x_yzzz, tr_0_0_yz_x_zzzz, tr_x_xx, tr_x_xxxxyz, tr_x_xxxy, tr_x_xxxyyz, tr_x_xxxyzz, tr_x_xxxz, tr_x_xxyy, tr_x_xxyyyz, tr_x_xxyyzz, tr_x_xxyz, tr_x_xxyzzz, tr_x_xxzz, tr_x_xy, tr_x_xyyy, tr_x_xyyyyz, tr_x_xyyyzz, tr_x_xyyz, tr_x_xyyzzz, tr_x_xyzz, tr_x_xyzzzz, tr_x_xz, tr_x_xzzz, tr_x_yy, tr_x_yyyy, tr_x_yyyyyz, tr_x_yyyyzz, tr_x_yyyz, tr_x_yyyzzz, tr_x_yyzz, tr_x_yyzzzz, tr_x_yz, tr_x_yzzz, tr_x_yzzzzz, tr_x_zz, tr_x_zzzz, tr_xy_xxx, tr_xy_xxxxz, tr_xy_xxxyz, tr_xy_xxxzz, tr_xy_xxy, tr_xy_xxyyz, tr_xy_xxyzz, tr_xy_xxz, tr_xy_xxzzz, tr_xy_xyy, tr_xy_xyyyz, tr_xy_xyyzz, tr_xy_xyz, tr_xy_xyzzz, tr_xy_xzz, tr_xy_xzzzz, tr_xy_yyy, tr_xy_yyyyz, tr_xy_yyyzz, tr_xy_yyz, tr_xy_yyzzz, tr_xy_yzz, tr_xy_yzzzz, tr_xy_zzz, tr_xy_zzzzz, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xzzz, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yzzz, tr_xyz_zzzz, tr_xz_xxx, tr_xz_xxxxy, tr_xz_xxxyy, tr_xz_xxxyz, tr_xz_xxy, tr_xz_xxyyy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xyy, tr_xz_xyyyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyyyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_x_xxxx[i] = 4.0 * tr_x_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxxy[i] = -2.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xz_xxx[i] * tbe_0 + 4.0 * tr_xz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxxz[i] = -2.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxx[i] * tbe_0 + 4.0 * tr_xy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxyy[i] = -4.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xxy[i] * tbe_0 + 4.0 * tr_xz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxyz[i] = tr_x_xx[i] - 2.0 * tr_x_xxzz[i] * tke_0 - 2.0 * tr_x_xxyy[i] * tke_0 + 4.0 * tr_x_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_xxz[i] * tbe_0 + 4.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xxy[i] * tbe_0 + 4.0 * tr_xy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xxzz[i] = -4.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xxz[i] * tbe_0 + 4.0 * tr_xy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xyyy[i] = -6.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xz_xyy[i] * tbe_0 + 4.0 * tr_xz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xyyz[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_x_xyzz[i] * tke_0 - 2.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_x_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xyz[i] * tbe_0 + 4.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_xyy[i] * tbe_0 + 4.0 * tr_xy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xyzz[i] = 2.0 * tr_x_xz[i] - 2.0 * tr_x_xzzz[i] * tke_0 - 4.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_xzz[i] * tbe_0 + 4.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_xyz[i] * tbe_0 + 4.0 * tr_xy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_xzzz[i] = -6.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xy_xzz[i] * tbe_0 + 4.0 * tr_xy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yyyy[i] = -8.0 * tr_x_yyyz[i] * tke_0 + 4.0 * tr_x_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xz_yyy[i] * tbe_0 + 4.0 * tr_xz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yyyz[i] = 3.0 * tr_x_yy[i] - 6.0 * tr_x_yyzz[i] * tke_0 - 2.0 * tr_x_yyyy[i] * tke_0 + 4.0 * tr_x_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xz_yyz[i] * tbe_0 + 4.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xy_yyy[i] * tbe_0 + 4.0 * tr_xy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yyzz[i] = 4.0 * tr_x_yz[i] - 4.0 * tr_x_yzzz[i] * tke_0 - 4.0 * tr_x_yyyz[i] * tke_0 + 4.0 * tr_x_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yzz[i] * tbe_0 + 4.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xy_yyz[i] * tbe_0 + 4.0 * tr_xy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_yzzz[i] = 3.0 * tr_x_zz[i] - 2.0 * tr_x_zzzz[i] * tke_0 - 6.0 * tr_x_yyzz[i] * tke_0 + 4.0 * tr_x_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xz_zzz[i] * tbe_0 + 4.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xy_yzz[i] * tbe_0 + 4.0 * tr_xy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_x_zzzz[i] = -8.0 * tr_x_yzzz[i] * tke_0 + 4.0 * tr_x_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xy_zzz[i] * tbe_0 + 4.0 * tr_xy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 195-210 components of targeted buffer : PG

    auto tr_0_0_yz_y_xxxx = pbuffer.data(idx_op_geom_020_pg + 195);

    auto tr_0_0_yz_y_xxxy = pbuffer.data(idx_op_geom_020_pg + 196);

    auto tr_0_0_yz_y_xxxz = pbuffer.data(idx_op_geom_020_pg + 197);

    auto tr_0_0_yz_y_xxyy = pbuffer.data(idx_op_geom_020_pg + 198);

    auto tr_0_0_yz_y_xxyz = pbuffer.data(idx_op_geom_020_pg + 199);

    auto tr_0_0_yz_y_xxzz = pbuffer.data(idx_op_geom_020_pg + 200);

    auto tr_0_0_yz_y_xyyy = pbuffer.data(idx_op_geom_020_pg + 201);

    auto tr_0_0_yz_y_xyyz = pbuffer.data(idx_op_geom_020_pg + 202);

    auto tr_0_0_yz_y_xyzz = pbuffer.data(idx_op_geom_020_pg + 203);

    auto tr_0_0_yz_y_xzzz = pbuffer.data(idx_op_geom_020_pg + 204);

    auto tr_0_0_yz_y_yyyy = pbuffer.data(idx_op_geom_020_pg + 205);

    auto tr_0_0_yz_y_yyyz = pbuffer.data(idx_op_geom_020_pg + 206);

    auto tr_0_0_yz_y_yyzz = pbuffer.data(idx_op_geom_020_pg + 207);

    auto tr_0_0_yz_y_yzzz = pbuffer.data(idx_op_geom_020_pg + 208);

    auto tr_0_0_yz_y_zzzz = pbuffer.data(idx_op_geom_020_pg + 209);

    #pragma omp simd aligned(tr_0_0_yz_y_xxxx, tr_0_0_yz_y_xxxy, tr_0_0_yz_y_xxxz, tr_0_0_yz_y_xxyy, tr_0_0_yz_y_xxyz, tr_0_0_yz_y_xxzz, tr_0_0_yz_y_xyyy, tr_0_0_yz_y_xyyz, tr_0_0_yz_y_xyzz, tr_0_0_yz_y_xzzz, tr_0_0_yz_y_yyyy, tr_0_0_yz_y_yyyz, tr_0_0_yz_y_yyzz, tr_0_0_yz_y_yzzz, tr_0_0_yz_y_zzzz, tr_0_xxx, tr_0_xxxxz, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_y_xx, tr_y_xxxxyz, tr_y_xxxy, tr_y_xxxyyz, tr_y_xxxyzz, tr_y_xxxz, tr_y_xxyy, tr_y_xxyyyz, tr_y_xxyyzz, tr_y_xxyz, tr_y_xxyzzz, tr_y_xxzz, tr_y_xy, tr_y_xyyy, tr_y_xyyyyz, tr_y_xyyyzz, tr_y_xyyz, tr_y_xyyzzz, tr_y_xyzz, tr_y_xyzzzz, tr_y_xz, tr_y_xzzz, tr_y_yy, tr_y_yyyy, tr_y_yyyyyz, tr_y_yyyyzz, tr_y_yyyz, tr_y_yyyzzz, tr_y_yyzz, tr_y_yyzzzz, tr_y_yz, tr_y_yzzz, tr_y_yzzzzz, tr_y_zz, tr_y_zzzz, tr_yy_xxx, tr_yy_xxxxz, tr_yy_xxxyz, tr_yy_xxxzz, tr_yy_xxy, tr_yy_xxyyz, tr_yy_xxyzz, tr_yy_xxz, tr_yy_xxzzz, tr_yy_xyy, tr_yy_xyyyz, tr_yy_xyyzz, tr_yy_xyz, tr_yy_xyzzz, tr_yy_xzz, tr_yy_xzzzz, tr_yy_yyy, tr_yy_yyyyz, tr_yy_yyyzz, tr_yy_yyz, tr_yy_yyzzz, tr_yy_yzz, tr_yy_yzzzz, tr_yy_zzz, tr_yy_zzzzz, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xzzz, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yzzz, tr_yyz_zzzz, tr_yz_xxx, tr_yz_xxxxy, tr_yz_xxxyy, tr_yz_xxxyz, tr_yz_xxy, tr_yz_xxyyy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xyy, tr_yz_xyyyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyyyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_z_xxxx, tr_z_xxxy, tr_z_xxxz, tr_z_xxyy, tr_z_xxyz, tr_z_xxzz, tr_z_xyyy, tr_z_xyyz, tr_z_xyzz, tr_z_xzzz, tr_z_yyyy, tr_z_yyyz, tr_z_yyzz, tr_z_yzzz, tr_z_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_y_xxxx[i] = -2.0 * tr_0_xxxxz[i] * tke_0 - 2.0 * tr_z_xxxx[i] * tbe_0 + 4.0 * tr_y_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxxy[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_z_xxxy[i] * tbe_0 - 2.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxxz[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxzz[i] * tke_0 - 2.0 * tr_z_xxxz[i] * tbe_0 - 2.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxx[i] * tbe_0 + 4.0 * tr_yy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxyy[i] = -2.0 * tr_0_xxyyz[i] * tke_0 - 2.0 * tr_z_xxyy[i] * tbe_0 - 4.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxyz[i] = tr_0_xxy[i] - 2.0 * tr_0_xxyzz[i] * tke_0 - 2.0 * tr_z_xxyz[i] * tbe_0 + tr_y_xx[i] - 2.0 * tr_y_xxzz[i] * tke_0 - 2.0 * tr_y_xxyy[i] * tke_0 + 4.0 * tr_y_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xxy[i] * tbe_0 + 4.0 * tr_yy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xxzz[i] = 2.0 * tr_0_xxz[i] - 2.0 * tr_0_xxzzz[i] * tke_0 - 2.0 * tr_z_xxzz[i] * tbe_0 - 4.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yy_xxz[i] * tbe_0 + 4.0 * tr_yy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xyyy[i] = -2.0 * tr_0_xyyyz[i] * tke_0 - 2.0 * tr_z_xyyy[i] * tbe_0 - 6.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xyyz[i] = tr_0_xyy[i] - 2.0 * tr_0_xyyzz[i] * tke_0 - 2.0 * tr_z_xyyz[i] * tbe_0 + 2.0 * tr_y_xy[i] - 4.0 * tr_y_xyzz[i] * tke_0 - 2.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_y_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_xyy[i] * tbe_0 + 4.0 * tr_yy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xyzz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyzzz[i] * tke_0 - 2.0 * tr_z_xyzz[i] * tbe_0 + 2.0 * tr_y_xz[i] - 2.0 * tr_y_xzzz[i] * tke_0 - 4.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yy_xyz[i] * tbe_0 + 4.0 * tr_yy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_xzzz[i] = 3.0 * tr_0_xzz[i] - 2.0 * tr_0_xzzzz[i] * tke_0 - 2.0 * tr_z_xzzz[i] * tbe_0 - 6.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yy_xzz[i] * tbe_0 + 4.0 * tr_yy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yyyy[i] = -2.0 * tr_0_yyyyz[i] * tke_0 - 2.0 * tr_z_yyyy[i] * tbe_0 - 8.0 * tr_y_yyyz[i] * tke_0 + 4.0 * tr_y_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yyyz[i] = tr_0_yyy[i] - 2.0 * tr_0_yyyzz[i] * tke_0 - 2.0 * tr_z_yyyz[i] * tbe_0 + 3.0 * tr_y_yy[i] - 6.0 * tr_y_yyzz[i] * tke_0 - 2.0 * tr_y_yyyy[i] * tke_0 + 4.0 * tr_y_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yy_yyy[i] * tbe_0 + 4.0 * tr_yy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yyzz[i] = 2.0 * tr_0_yyz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 - 2.0 * tr_z_yyzz[i] * tbe_0 + 4.0 * tr_y_yz[i] - 4.0 * tr_y_yzzz[i] * tke_0 - 4.0 * tr_y_yyyz[i] * tke_0 + 4.0 * tr_y_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yy_yyz[i] * tbe_0 + 4.0 * tr_yy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_yzzz[i] = 3.0 * tr_0_yzz[i] - 2.0 * tr_0_yzzzz[i] * tke_0 - 2.0 * tr_z_yzzz[i] * tbe_0 + 3.0 * tr_y_zz[i] - 2.0 * tr_y_zzzz[i] * tke_0 - 6.0 * tr_y_yyzz[i] * tke_0 + 4.0 * tr_y_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yy_yzz[i] * tbe_0 + 4.0 * tr_yy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_y_zzzz[i] = 4.0 * tr_0_zzz[i] - 2.0 * tr_0_zzzzz[i] * tke_0 - 2.0 * tr_z_zzzz[i] * tbe_0 - 8.0 * tr_y_yzzz[i] * tke_0 + 4.0 * tr_y_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_yy_zzz[i] * tbe_0 + 4.0 * tr_yy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-225 components of targeted buffer : PG

    auto tr_0_0_yz_z_xxxx = pbuffer.data(idx_op_geom_020_pg + 210);

    auto tr_0_0_yz_z_xxxy = pbuffer.data(idx_op_geom_020_pg + 211);

    auto tr_0_0_yz_z_xxxz = pbuffer.data(idx_op_geom_020_pg + 212);

    auto tr_0_0_yz_z_xxyy = pbuffer.data(idx_op_geom_020_pg + 213);

    auto tr_0_0_yz_z_xxyz = pbuffer.data(idx_op_geom_020_pg + 214);

    auto tr_0_0_yz_z_xxzz = pbuffer.data(idx_op_geom_020_pg + 215);

    auto tr_0_0_yz_z_xyyy = pbuffer.data(idx_op_geom_020_pg + 216);

    auto tr_0_0_yz_z_xyyz = pbuffer.data(idx_op_geom_020_pg + 217);

    auto tr_0_0_yz_z_xyzz = pbuffer.data(idx_op_geom_020_pg + 218);

    auto tr_0_0_yz_z_xzzz = pbuffer.data(idx_op_geom_020_pg + 219);

    auto tr_0_0_yz_z_yyyy = pbuffer.data(idx_op_geom_020_pg + 220);

    auto tr_0_0_yz_z_yyyz = pbuffer.data(idx_op_geom_020_pg + 221);

    auto tr_0_0_yz_z_yyzz = pbuffer.data(idx_op_geom_020_pg + 222);

    auto tr_0_0_yz_z_yzzz = pbuffer.data(idx_op_geom_020_pg + 223);

    auto tr_0_0_yz_z_zzzz = pbuffer.data(idx_op_geom_020_pg + 224);

    #pragma omp simd aligned(tr_0_0_yz_z_xxxx, tr_0_0_yz_z_xxxy, tr_0_0_yz_z_xxxz, tr_0_0_yz_z_xxyy, tr_0_0_yz_z_xxyz, tr_0_0_yz_z_xxzz, tr_0_0_yz_z_xyyy, tr_0_0_yz_z_xyyz, tr_0_0_yz_z_xyzz, tr_0_0_yz_z_xzzz, tr_0_0_yz_z_yyyy, tr_0_0_yz_z_yyyz, tr_0_0_yz_z_yyzz, tr_0_0_yz_z_yzzz, tr_0_0_yz_z_zzzz, tr_0_xxx, tr_0_xxxxy, tr_0_xxxyy, tr_0_xxxyz, tr_0_xxy, tr_0_xxyyy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xyy, tr_0_xyyyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_yyy, tr_0_yyyyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_y_xxxx, tr_y_xxxy, tr_y_xxxz, tr_y_xxyy, tr_y_xxyz, tr_y_xxzz, tr_y_xyyy, tr_y_xyyz, tr_y_xyzz, tr_y_xzzz, tr_y_yyyy, tr_y_yyyz, tr_y_yyzz, tr_y_yzzz, tr_y_zzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, tr_z_xx, tr_z_xxxxyz, tr_z_xxxy, tr_z_xxxyyz, tr_z_xxxyzz, tr_z_xxxz, tr_z_xxyy, tr_z_xxyyyz, tr_z_xxyyzz, tr_z_xxyz, tr_z_xxyzzz, tr_z_xxzz, tr_z_xy, tr_z_xyyy, tr_z_xyyyyz, tr_z_xyyyzz, tr_z_xyyz, tr_z_xyyzzz, tr_z_xyzz, tr_z_xyzzzz, tr_z_xz, tr_z_xzzz, tr_z_yy, tr_z_yyyy, tr_z_yyyyyz, tr_z_yyyyzz, tr_z_yyyz, tr_z_yyyzzz, tr_z_yyzz, tr_z_yyzzzz, tr_z_yz, tr_z_yzzz, tr_z_yzzzzz, tr_z_zz, tr_z_zzzz, tr_zz_xxx, tr_zz_xxxxy, tr_zz_xxxyy, tr_zz_xxxyz, tr_zz_xxy, tr_zz_xxyyy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xyy, tr_zz_xyyyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyyyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_z_xxxx[i] = -2.0 * tr_0_xxxxy[i] * tke_0 + 4.0 * tr_z_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxxx[i] * tbe_0 + 4.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxxy[i] = tr_0_xxx[i] - 2.0 * tr_0_xxxyy[i] * tke_0 - 2.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zz_xxx[i] * tbe_0 + 4.0 * tr_zz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxxy[i] * tbe_0 + 4.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxxz[i] = -2.0 * tr_0_xxxyz[i] * tke_0 - 2.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxxz[i] * tbe_0 - 2.0 * tr_yz_xxx[i] * tbe_0 + 4.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxyy[i] = 2.0 * tr_0_xxy[i] - 2.0 * tr_0_xxyyy[i] * tke_0 - 4.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xxy[i] * tbe_0 + 4.0 * tr_zz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxyy[i] * tbe_0 + 4.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxyz[i] = tr_0_xxz[i] - 2.0 * tr_0_xxyyz[i] * tke_0 + tr_z_xx[i] - 2.0 * tr_z_xxzz[i] * tke_0 - 2.0 * tr_z_xxyy[i] * tke_0 + 4.0 * tr_z_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_xxz[i] * tbe_0 + 4.0 * tr_zz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxyz[i] * tbe_0 - 2.0 * tr_yz_xxy[i] * tbe_0 + 4.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xxzz[i] = -2.0 * tr_0_xxyzz[i] * tke_0 - 4.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxzz[i] * tbe_0 - 4.0 * tr_yz_xxz[i] * tbe_0 + 4.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xyyy[i] = 3.0 * tr_0_xyy[i] - 2.0 * tr_0_xyyyy[i] * tke_0 - 6.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_zz_xyy[i] * tbe_0 + 4.0 * tr_zz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyyy[i] * tbe_0 + 4.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xyyz[i] = 2.0 * tr_0_xyz[i] - 2.0 * tr_0_xyyyz[i] * tke_0 + 2.0 * tr_z_xy[i] - 4.0 * tr_z_xyzz[i] * tke_0 - 2.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_z_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xyz[i] * tbe_0 + 4.0 * tr_zz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyyz[i] * tbe_0 - 2.0 * tr_yz_xyy[i] * tbe_0 + 4.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xyzz[i] = tr_0_xzz[i] - 2.0 * tr_0_xyyzz[i] * tke_0 + 2.0 * tr_z_xz[i] - 2.0 * tr_z_xzzz[i] * tke_0 - 4.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_xzz[i] * tbe_0 + 4.0 * tr_zz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyzz[i] * tbe_0 - 4.0 * tr_yz_xyz[i] * tbe_0 + 4.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_xzzz[i] = -2.0 * tr_0_xyzzz[i] * tke_0 - 6.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xzzz[i] * tbe_0 - 6.0 * tr_yz_xzz[i] * tbe_0 + 4.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yyyy[i] = 4.0 * tr_0_yyy[i] - 2.0 * tr_0_yyyyy[i] * tke_0 - 8.0 * tr_z_yyyz[i] * tke_0 + 4.0 * tr_z_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_zz_yyy[i] * tbe_0 + 4.0 * tr_zz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyyy[i] * tbe_0 + 4.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yyyz[i] = 3.0 * tr_0_yyz[i] - 2.0 * tr_0_yyyyz[i] * tke_0 + 3.0 * tr_z_yy[i] - 6.0 * tr_z_yyzz[i] * tke_0 - 2.0 * tr_z_yyyy[i] * tke_0 + 4.0 * tr_z_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_zz_yyz[i] * tbe_0 + 4.0 * tr_zz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyyz[i] * tbe_0 - 2.0 * tr_yz_yyy[i] * tbe_0 + 4.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yyzz[i] = 2.0 * tr_0_yzz[i] - 2.0 * tr_0_yyyzz[i] * tke_0 + 4.0 * tr_z_yz[i] - 4.0 * tr_z_yzzz[i] * tke_0 - 4.0 * tr_z_yyyz[i] * tke_0 + 4.0 * tr_z_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_yzz[i] * tbe_0 + 4.0 * tr_zz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyzz[i] * tbe_0 - 4.0 * tr_yz_yyz[i] * tbe_0 + 4.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_yzzz[i] = tr_0_zzz[i] - 2.0 * tr_0_yyzzz[i] * tke_0 + 3.0 * tr_z_zz[i] - 2.0 * tr_z_zzzz[i] * tke_0 - 6.0 * tr_z_yyzz[i] * tke_0 + 4.0 * tr_z_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_zz_zzz[i] * tbe_0 + 4.0 * tr_zz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yzzz[i] * tbe_0 - 6.0 * tr_yz_yzz[i] * tbe_0 + 4.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_z_zzzz[i] = -2.0 * tr_0_yzzzz[i] * tke_0 - 8.0 * tr_z_yzzz[i] * tke_0 + 4.0 * tr_z_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_y_zzzz[i] * tbe_0 - 8.0 * tr_yz_zzz[i] * tbe_0 + 4.0 * tr_yz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 225-240 components of targeted buffer : PG

    auto tr_0_0_zz_x_xxxx = pbuffer.data(idx_op_geom_020_pg + 225);

    auto tr_0_0_zz_x_xxxy = pbuffer.data(idx_op_geom_020_pg + 226);

    auto tr_0_0_zz_x_xxxz = pbuffer.data(idx_op_geom_020_pg + 227);

    auto tr_0_0_zz_x_xxyy = pbuffer.data(idx_op_geom_020_pg + 228);

    auto tr_0_0_zz_x_xxyz = pbuffer.data(idx_op_geom_020_pg + 229);

    auto tr_0_0_zz_x_xxzz = pbuffer.data(idx_op_geom_020_pg + 230);

    auto tr_0_0_zz_x_xyyy = pbuffer.data(idx_op_geom_020_pg + 231);

    auto tr_0_0_zz_x_xyyz = pbuffer.data(idx_op_geom_020_pg + 232);

    auto tr_0_0_zz_x_xyzz = pbuffer.data(idx_op_geom_020_pg + 233);

    auto tr_0_0_zz_x_xzzz = pbuffer.data(idx_op_geom_020_pg + 234);

    auto tr_0_0_zz_x_yyyy = pbuffer.data(idx_op_geom_020_pg + 235);

    auto tr_0_0_zz_x_yyyz = pbuffer.data(idx_op_geom_020_pg + 236);

    auto tr_0_0_zz_x_yyzz = pbuffer.data(idx_op_geom_020_pg + 237);

    auto tr_0_0_zz_x_yzzz = pbuffer.data(idx_op_geom_020_pg + 238);

    auto tr_0_0_zz_x_zzzz = pbuffer.data(idx_op_geom_020_pg + 239);

    #pragma omp simd aligned(tr_0_0_zz_x_xxxx, tr_0_0_zz_x_xxxy, tr_0_0_zz_x_xxxz, tr_0_0_zz_x_xxyy, tr_0_0_zz_x_xxyz, tr_0_0_zz_x_xxzz, tr_0_0_zz_x_xyyy, tr_0_0_zz_x_xyyz, tr_0_0_zz_x_xyzz, tr_0_0_zz_x_xzzz, tr_0_0_zz_x_yyyy, tr_0_0_zz_x_yyyz, tr_0_0_zz_x_yyzz, tr_0_0_zz_x_yzzz, tr_0_0_zz_x_zzzz, tr_x_xx, tr_x_xxxx, tr_x_xxxxzz, tr_x_xxxy, tr_x_xxxyzz, tr_x_xxxz, tr_x_xxxzzz, tr_x_xxyy, tr_x_xxyyzz, tr_x_xxyz, tr_x_xxyzzz, tr_x_xxzz, tr_x_xxzzzz, tr_x_xy, tr_x_xyyy, tr_x_xyyyzz, tr_x_xyyz, tr_x_xyyzzz, tr_x_xyzz, tr_x_xyzzzz, tr_x_xz, tr_x_xzzz, tr_x_xzzzzz, tr_x_yy, tr_x_yyyy, tr_x_yyyyzz, tr_x_yyyz, tr_x_yyyzzz, tr_x_yyzz, tr_x_yyzzzz, tr_x_yz, tr_x_yzzz, tr_x_yzzzzz, tr_x_zz, tr_x_zzzz, tr_x_zzzzzz, tr_xz_xxx, tr_xz_xxxxz, tr_xz_xxxyz, tr_xz_xxxzz, tr_xz_xxy, tr_xz_xxyyz, tr_xz_xxyzz, tr_xz_xxz, tr_xz_xxzzz, tr_xz_xyy, tr_xz_xyyyz, tr_xz_xyyzz, tr_xz_xyz, tr_xz_xyzzz, tr_xz_xzz, tr_xz_xzzzz, tr_xz_yyy, tr_xz_yyyyz, tr_xz_yyyzz, tr_xz_yyz, tr_xz_yyzzz, tr_xz_yzz, tr_xz_yzzzz, tr_xz_zzz, tr_xz_zzzzz, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xzzz, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yzzz, tr_xzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_x_xxxx[i] = -2.0 * tr_x_xxxx[i] * tbe_0 - 2.0 * tr_x_xxxx[i] * tke_0 + 4.0 * tr_x_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxxy[i] = -2.0 * tr_x_xxxy[i] * tbe_0 - 2.0 * tr_x_xxxy[i] * tke_0 + 4.0 * tr_x_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxxz[i] = -2.0 * tr_x_xxxz[i] * tbe_0 - 6.0 * tr_x_xxxz[i] * tke_0 + 4.0 * tr_x_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xxx[i] * tbe_0 + 8.0 * tr_xz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxyy[i] = -2.0 * tr_x_xxyy[i] * tbe_0 - 2.0 * tr_x_xxyy[i] * tke_0 + 4.0 * tr_x_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxyz[i] = -2.0 * tr_x_xxyz[i] * tbe_0 - 6.0 * tr_x_xxyz[i] * tke_0 + 4.0 * tr_x_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xxy[i] * tbe_0 + 8.0 * tr_xz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xxzz[i] = 2.0 * tr_x_xx[i] - 2.0 * tr_x_xxzz[i] * tbe_0 - 10.0 * tr_x_xxzz[i] * tke_0 + 4.0 * tr_x_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xxz[i] * tbe_0 + 8.0 * tr_xz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xyyy[i] = -2.0 * tr_x_xyyy[i] * tbe_0 - 2.0 * tr_x_xyyy[i] * tke_0 + 4.0 * tr_x_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xyyz[i] = -2.0 * tr_x_xyyz[i] * tbe_0 - 6.0 * tr_x_xyyz[i] * tke_0 + 4.0 * tr_x_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_xyy[i] * tbe_0 + 8.0 * tr_xz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xyzz[i] = 2.0 * tr_x_xy[i] - 2.0 * tr_x_xyzz[i] * tbe_0 - 10.0 * tr_x_xyzz[i] * tke_0 + 4.0 * tr_x_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xz_xyz[i] * tbe_0 + 8.0 * tr_xz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_xzzz[i] = 6.0 * tr_x_xz[i] - 2.0 * tr_x_xzzz[i] * tbe_0 - 14.0 * tr_x_xzzz[i] * tke_0 + 4.0 * tr_x_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xz_xzz[i] * tbe_0 + 8.0 * tr_xz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yyyy[i] = -2.0 * tr_x_yyyy[i] * tbe_0 - 2.0 * tr_x_yyyy[i] * tke_0 + 4.0 * tr_x_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yyyz[i] = -2.0 * tr_x_yyyz[i] * tbe_0 - 6.0 * tr_x_yyyz[i] * tke_0 + 4.0 * tr_x_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xz_yyy[i] * tbe_0 + 8.0 * tr_xz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yyzz[i] = 2.0 * tr_x_yy[i] - 2.0 * tr_x_yyzz[i] * tbe_0 - 10.0 * tr_x_yyzz[i] * tke_0 + 4.0 * tr_x_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xz_yyz[i] * tbe_0 + 8.0 * tr_xz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_yzzz[i] = 6.0 * tr_x_yz[i] - 2.0 * tr_x_yzzz[i] * tbe_0 - 14.0 * tr_x_yzzz[i] * tke_0 + 4.0 * tr_x_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xz_yzz[i] * tbe_0 + 8.0 * tr_xz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_x_zzzz[i] = 12.0 * tr_x_zz[i] - 2.0 * tr_x_zzzz[i] * tbe_0 - 18.0 * tr_x_zzzz[i] * tke_0 + 4.0 * tr_x_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xz_zzz[i] * tbe_0 + 8.0 * tr_xz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-255 components of targeted buffer : PG

    auto tr_0_0_zz_y_xxxx = pbuffer.data(idx_op_geom_020_pg + 240);

    auto tr_0_0_zz_y_xxxy = pbuffer.data(idx_op_geom_020_pg + 241);

    auto tr_0_0_zz_y_xxxz = pbuffer.data(idx_op_geom_020_pg + 242);

    auto tr_0_0_zz_y_xxyy = pbuffer.data(idx_op_geom_020_pg + 243);

    auto tr_0_0_zz_y_xxyz = pbuffer.data(idx_op_geom_020_pg + 244);

    auto tr_0_0_zz_y_xxzz = pbuffer.data(idx_op_geom_020_pg + 245);

    auto tr_0_0_zz_y_xyyy = pbuffer.data(idx_op_geom_020_pg + 246);

    auto tr_0_0_zz_y_xyyz = pbuffer.data(idx_op_geom_020_pg + 247);

    auto tr_0_0_zz_y_xyzz = pbuffer.data(idx_op_geom_020_pg + 248);

    auto tr_0_0_zz_y_xzzz = pbuffer.data(idx_op_geom_020_pg + 249);

    auto tr_0_0_zz_y_yyyy = pbuffer.data(idx_op_geom_020_pg + 250);

    auto tr_0_0_zz_y_yyyz = pbuffer.data(idx_op_geom_020_pg + 251);

    auto tr_0_0_zz_y_yyzz = pbuffer.data(idx_op_geom_020_pg + 252);

    auto tr_0_0_zz_y_yzzz = pbuffer.data(idx_op_geom_020_pg + 253);

    auto tr_0_0_zz_y_zzzz = pbuffer.data(idx_op_geom_020_pg + 254);

    #pragma omp simd aligned(tr_0_0_zz_y_xxxx, tr_0_0_zz_y_xxxy, tr_0_0_zz_y_xxxz, tr_0_0_zz_y_xxyy, tr_0_0_zz_y_xxyz, tr_0_0_zz_y_xxzz, tr_0_0_zz_y_xyyy, tr_0_0_zz_y_xyyz, tr_0_0_zz_y_xyzz, tr_0_0_zz_y_xzzz, tr_0_0_zz_y_yyyy, tr_0_0_zz_y_yyyz, tr_0_0_zz_y_yyzz, tr_0_0_zz_y_yzzz, tr_0_0_zz_y_zzzz, tr_y_xx, tr_y_xxxx, tr_y_xxxxzz, tr_y_xxxy, tr_y_xxxyzz, tr_y_xxxz, tr_y_xxxzzz, tr_y_xxyy, tr_y_xxyyzz, tr_y_xxyz, tr_y_xxyzzz, tr_y_xxzz, tr_y_xxzzzz, tr_y_xy, tr_y_xyyy, tr_y_xyyyzz, tr_y_xyyz, tr_y_xyyzzz, tr_y_xyzz, tr_y_xyzzzz, tr_y_xz, tr_y_xzzz, tr_y_xzzzzz, tr_y_yy, tr_y_yyyy, tr_y_yyyyzz, tr_y_yyyz, tr_y_yyyzzz, tr_y_yyzz, tr_y_yyzzzz, tr_y_yz, tr_y_yzzz, tr_y_yzzzzz, tr_y_zz, tr_y_zzzz, tr_y_zzzzzz, tr_yz_xxx, tr_yz_xxxxz, tr_yz_xxxyz, tr_yz_xxxzz, tr_yz_xxy, tr_yz_xxyyz, tr_yz_xxyzz, tr_yz_xxz, tr_yz_xxzzz, tr_yz_xyy, tr_yz_xyyyz, tr_yz_xyyzz, tr_yz_xyz, tr_yz_xyzzz, tr_yz_xzz, tr_yz_xzzzz, tr_yz_yyy, tr_yz_yyyyz, tr_yz_yyyzz, tr_yz_yyz, tr_yz_yyzzz, tr_yz_yzz, tr_yz_yzzzz, tr_yz_zzz, tr_yz_zzzzz, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xzzz, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yzzz, tr_yzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_y_xxxx[i] = -2.0 * tr_y_xxxx[i] * tbe_0 - 2.0 * tr_y_xxxx[i] * tke_0 + 4.0 * tr_y_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxxy[i] = -2.0 * tr_y_xxxy[i] * tbe_0 - 2.0 * tr_y_xxxy[i] * tke_0 + 4.0 * tr_y_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxxz[i] = -2.0 * tr_y_xxxz[i] * tbe_0 - 6.0 * tr_y_xxxz[i] * tke_0 + 4.0 * tr_y_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xxx[i] * tbe_0 + 8.0 * tr_yz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxyy[i] = -2.0 * tr_y_xxyy[i] * tbe_0 - 2.0 * tr_y_xxyy[i] * tke_0 + 4.0 * tr_y_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxyz[i] = -2.0 * tr_y_xxyz[i] * tbe_0 - 6.0 * tr_y_xxyz[i] * tke_0 + 4.0 * tr_y_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xxy[i] * tbe_0 + 8.0 * tr_yz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xxzz[i] = 2.0 * tr_y_xx[i] - 2.0 * tr_y_xxzz[i] * tbe_0 - 10.0 * tr_y_xxzz[i] * tke_0 + 4.0 * tr_y_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yz_xxz[i] * tbe_0 + 8.0 * tr_yz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xyyy[i] = -2.0 * tr_y_xyyy[i] * tbe_0 - 2.0 * tr_y_xyyy[i] * tke_0 + 4.0 * tr_y_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xyyz[i] = -2.0 * tr_y_xyyz[i] * tbe_0 - 6.0 * tr_y_xyyz[i] * tke_0 + 4.0 * tr_y_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_xyy[i] * tbe_0 + 8.0 * tr_yz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xyzz[i] = 2.0 * tr_y_xy[i] - 2.0 * tr_y_xyzz[i] * tbe_0 - 10.0 * tr_y_xyzz[i] * tke_0 + 4.0 * tr_y_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yz_xyz[i] * tbe_0 + 8.0 * tr_yz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_xzzz[i] = 6.0 * tr_y_xz[i] - 2.0 * tr_y_xzzz[i] * tbe_0 - 14.0 * tr_y_xzzz[i] * tke_0 + 4.0 * tr_y_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yz_xzz[i] * tbe_0 + 8.0 * tr_yz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yyyy[i] = -2.0 * tr_y_yyyy[i] * tbe_0 - 2.0 * tr_y_yyyy[i] * tke_0 + 4.0 * tr_y_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yyyz[i] = -2.0 * tr_y_yyyz[i] * tbe_0 - 6.0 * tr_y_yyyz[i] * tke_0 + 4.0 * tr_y_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yz_yyy[i] * tbe_0 + 8.0 * tr_yz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yyzz[i] = 2.0 * tr_y_yy[i] - 2.0 * tr_y_yyzz[i] * tbe_0 - 10.0 * tr_y_yyzz[i] * tke_0 + 4.0 * tr_y_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yz_yyz[i] * tbe_0 + 8.0 * tr_yz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_yzzz[i] = 6.0 * tr_y_yz[i] - 2.0 * tr_y_yzzz[i] * tbe_0 - 14.0 * tr_y_yzzz[i] * tke_0 + 4.0 * tr_y_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yz_yzz[i] * tbe_0 + 8.0 * tr_yz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_y_zzzz[i] = 12.0 * tr_y_zz[i] - 2.0 * tr_y_zzzz[i] * tbe_0 - 18.0 * tr_y_zzzz[i] * tke_0 + 4.0 * tr_y_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_yz_zzz[i] * tbe_0 + 8.0 * tr_yz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 255-270 components of targeted buffer : PG

    auto tr_0_0_zz_z_xxxx = pbuffer.data(idx_op_geom_020_pg + 255);

    auto tr_0_0_zz_z_xxxy = pbuffer.data(idx_op_geom_020_pg + 256);

    auto tr_0_0_zz_z_xxxz = pbuffer.data(idx_op_geom_020_pg + 257);

    auto tr_0_0_zz_z_xxyy = pbuffer.data(idx_op_geom_020_pg + 258);

    auto tr_0_0_zz_z_xxyz = pbuffer.data(idx_op_geom_020_pg + 259);

    auto tr_0_0_zz_z_xxzz = pbuffer.data(idx_op_geom_020_pg + 260);

    auto tr_0_0_zz_z_xyyy = pbuffer.data(idx_op_geom_020_pg + 261);

    auto tr_0_0_zz_z_xyyz = pbuffer.data(idx_op_geom_020_pg + 262);

    auto tr_0_0_zz_z_xyzz = pbuffer.data(idx_op_geom_020_pg + 263);

    auto tr_0_0_zz_z_xzzz = pbuffer.data(idx_op_geom_020_pg + 264);

    auto tr_0_0_zz_z_yyyy = pbuffer.data(idx_op_geom_020_pg + 265);

    auto tr_0_0_zz_z_yyyz = pbuffer.data(idx_op_geom_020_pg + 266);

    auto tr_0_0_zz_z_yyzz = pbuffer.data(idx_op_geom_020_pg + 267);

    auto tr_0_0_zz_z_yzzz = pbuffer.data(idx_op_geom_020_pg + 268);

    auto tr_0_0_zz_z_zzzz = pbuffer.data(idx_op_geom_020_pg + 269);

    #pragma omp simd aligned(tr_0_0_zz_z_xxxx, tr_0_0_zz_z_xxxy, tr_0_0_zz_z_xxxz, tr_0_0_zz_z_xxyy, tr_0_0_zz_z_xxyz, tr_0_0_zz_z_xxzz, tr_0_0_zz_z_xyyy, tr_0_0_zz_z_xyyz, tr_0_0_zz_z_xyzz, tr_0_0_zz_z_xzzz, tr_0_0_zz_z_yyyy, tr_0_0_zz_z_yyyz, tr_0_0_zz_z_yyzz, tr_0_0_zz_z_yzzz, tr_0_0_zz_z_zzzz, tr_0_xxx, tr_0_xxxxz, tr_0_xxxyz, tr_0_xxxzz, tr_0_xxy, tr_0_xxyyz, tr_0_xxyzz, tr_0_xxz, tr_0_xxzzz, tr_0_xyy, tr_0_xyyyz, tr_0_xyyzz, tr_0_xyz, tr_0_xyzzz, tr_0_xzz, tr_0_xzzzz, tr_0_yyy, tr_0_yyyyz, tr_0_yyyzz, tr_0_yyz, tr_0_yyzzz, tr_0_yzz, tr_0_yzzzz, tr_0_zzz, tr_0_zzzzz, tr_z_xx, tr_z_xxxx, tr_z_xxxxzz, tr_z_xxxy, tr_z_xxxyzz, tr_z_xxxz, tr_z_xxxzzz, tr_z_xxyy, tr_z_xxyyzz, tr_z_xxyz, tr_z_xxyzzz, tr_z_xxzz, tr_z_xxzzzz, tr_z_xy, tr_z_xyyy, tr_z_xyyyzz, tr_z_xyyz, tr_z_xyyzzz, tr_z_xyzz, tr_z_xyzzzz, tr_z_xz, tr_z_xzzz, tr_z_xzzzzz, tr_z_yy, tr_z_yyyy, tr_z_yyyyzz, tr_z_yyyz, tr_z_yyyzzz, tr_z_yyzz, tr_z_yyzzzz, tr_z_yz, tr_z_yzzz, tr_z_yzzzzz, tr_z_zz, tr_z_zzzz, tr_z_zzzzzz, tr_zz_xxx, tr_zz_xxxxz, tr_zz_xxxyz, tr_zz_xxxzz, tr_zz_xxy, tr_zz_xxyyz, tr_zz_xxyzz, tr_zz_xxz, tr_zz_xxzzz, tr_zz_xyy, tr_zz_xyyyz, tr_zz_xyyzz, tr_zz_xyz, tr_zz_xyzzz, tr_zz_xzz, tr_zz_xzzzz, tr_zz_yyy, tr_zz_yyyyz, tr_zz_yyyzz, tr_zz_yyz, tr_zz_yyzzz, tr_zz_yzz, tr_zz_yzzzz, tr_zz_zzz, tr_zz_zzzzz, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xzzz, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yzzz, tr_zzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_z_xxxx[i] = -4.0 * tr_0_xxxxz[i] * tke_0 - 6.0 * tr_z_xxxx[i] * tbe_0 - 2.0 * tr_z_xxxx[i] * tke_0 + 4.0 * tr_z_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxxy[i] = -4.0 * tr_0_xxxyz[i] * tke_0 - 6.0 * tr_z_xxxy[i] * tbe_0 - 2.0 * tr_z_xxxy[i] * tke_0 + 4.0 * tr_z_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxxz[i] = 2.0 * tr_0_xxx[i] - 4.0 * tr_0_xxxzz[i] * tke_0 - 6.0 * tr_z_xxxz[i] * tbe_0 - 6.0 * tr_z_xxxz[i] * tke_0 + 4.0 * tr_z_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xxx[i] * tbe_0 + 8.0 * tr_zz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxyy[i] = -4.0 * tr_0_xxyyz[i] * tke_0 - 6.0 * tr_z_xxyy[i] * tbe_0 - 2.0 * tr_z_xxyy[i] * tke_0 + 4.0 * tr_z_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxyz[i] = 2.0 * tr_0_xxy[i] - 4.0 * tr_0_xxyzz[i] * tke_0 - 6.0 * tr_z_xxyz[i] * tbe_0 - 6.0 * tr_z_xxyz[i] * tke_0 + 4.0 * tr_z_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xxy[i] * tbe_0 + 8.0 * tr_zz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xxzz[i] = 4.0 * tr_0_xxz[i] - 4.0 * tr_0_xxzzz[i] * tke_0 + 2.0 * tr_z_xx[i] - 6.0 * tr_z_xxzz[i] * tbe_0 - 10.0 * tr_z_xxzz[i] * tke_0 + 4.0 * tr_z_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zz_xxz[i] * tbe_0 + 8.0 * tr_zz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xyyy[i] = -4.0 * tr_0_xyyyz[i] * tke_0 - 6.0 * tr_z_xyyy[i] * tbe_0 - 2.0 * tr_z_xyyy[i] * tke_0 + 4.0 * tr_z_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xyyz[i] = 2.0 * tr_0_xyy[i] - 4.0 * tr_0_xyyzz[i] * tke_0 - 6.0 * tr_z_xyyz[i] * tbe_0 - 6.0 * tr_z_xyyz[i] * tke_0 + 4.0 * tr_z_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_xyy[i] * tbe_0 + 8.0 * tr_zz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xyzz[i] = 4.0 * tr_0_xyz[i] - 4.0 * tr_0_xyzzz[i] * tke_0 + 2.0 * tr_z_xy[i] - 6.0 * tr_z_xyzz[i] * tbe_0 - 10.0 * tr_z_xyzz[i] * tke_0 + 4.0 * tr_z_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zz_xyz[i] * tbe_0 + 8.0 * tr_zz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_xzzz[i] = 6.0 * tr_0_xzz[i] - 4.0 * tr_0_xzzzz[i] * tke_0 + 6.0 * tr_z_xz[i] - 6.0 * tr_z_xzzz[i] * tbe_0 - 14.0 * tr_z_xzzz[i] * tke_0 + 4.0 * tr_z_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zz_xzz[i] * tbe_0 + 8.0 * tr_zz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yyyy[i] = -4.0 * tr_0_yyyyz[i] * tke_0 - 6.0 * tr_z_yyyy[i] * tbe_0 - 2.0 * tr_z_yyyy[i] * tke_0 + 4.0 * tr_z_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yyyz[i] = 2.0 * tr_0_yyy[i] - 4.0 * tr_0_yyyzz[i] * tke_0 - 6.0 * tr_z_yyyz[i] * tbe_0 - 6.0 * tr_z_yyyz[i] * tke_0 + 4.0 * tr_z_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zz_yyy[i] * tbe_0 + 8.0 * tr_zz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yyzz[i] = 4.0 * tr_0_yyz[i] - 4.0 * tr_0_yyzzz[i] * tke_0 + 2.0 * tr_z_yy[i] - 6.0 * tr_z_yyzz[i] * tbe_0 - 10.0 * tr_z_yyzz[i] * tke_0 + 4.0 * tr_z_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zz_yyz[i] * tbe_0 + 8.0 * tr_zz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_yzzz[i] = 6.0 * tr_0_yzz[i] - 4.0 * tr_0_yzzzz[i] * tke_0 + 6.0 * tr_z_yz[i] - 6.0 * tr_z_yzzz[i] * tbe_0 - 14.0 * tr_z_yzzz[i] * tke_0 + 4.0 * tr_z_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zz_yzz[i] * tbe_0 + 8.0 * tr_zz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_z_zzzz[i] = 8.0 * tr_0_zzz[i] - 4.0 * tr_0_zzzzz[i] * tke_0 + 12.0 * tr_z_zz[i] - 6.0 * tr_z_zzzz[i] * tbe_0 - 18.0 * tr_z_zzzz[i] * tke_0 + 4.0 * tr_z_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_zz_zzz[i] * tbe_0 + 8.0 * tr_zz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzz_zzzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

