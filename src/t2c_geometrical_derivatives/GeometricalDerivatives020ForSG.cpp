#include "GeometricalDerivatives020ForSG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_sg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_sg,
                         const int idx_op_sd,
                         const int idx_op_sg,
                         const int idx_op_si,
                         const int idx_op_pf,
                         const int idx_op_ph,
                         const int idx_op_dg,
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

    // Set up components of auxiliary buffer : SI

    auto tr_0_xxxxxx = pbuffer.data(idx_op_si);

    auto tr_0_xxxxxy = pbuffer.data(idx_op_si + 1);

    auto tr_0_xxxxxz = pbuffer.data(idx_op_si + 2);

    auto tr_0_xxxxyy = pbuffer.data(idx_op_si + 3);

    auto tr_0_xxxxyz = pbuffer.data(idx_op_si + 4);

    auto tr_0_xxxxzz = pbuffer.data(idx_op_si + 5);

    auto tr_0_xxxyyy = pbuffer.data(idx_op_si + 6);

    auto tr_0_xxxyyz = pbuffer.data(idx_op_si + 7);

    auto tr_0_xxxyzz = pbuffer.data(idx_op_si + 8);

    auto tr_0_xxxzzz = pbuffer.data(idx_op_si + 9);

    auto tr_0_xxyyyy = pbuffer.data(idx_op_si + 10);

    auto tr_0_xxyyyz = pbuffer.data(idx_op_si + 11);

    auto tr_0_xxyyzz = pbuffer.data(idx_op_si + 12);

    auto tr_0_xxyzzz = pbuffer.data(idx_op_si + 13);

    auto tr_0_xxzzzz = pbuffer.data(idx_op_si + 14);

    auto tr_0_xyyyyy = pbuffer.data(idx_op_si + 15);

    auto tr_0_xyyyyz = pbuffer.data(idx_op_si + 16);

    auto tr_0_xyyyzz = pbuffer.data(idx_op_si + 17);

    auto tr_0_xyyzzz = pbuffer.data(idx_op_si + 18);

    auto tr_0_xyzzzz = pbuffer.data(idx_op_si + 19);

    auto tr_0_xzzzzz = pbuffer.data(idx_op_si + 20);

    auto tr_0_yyyyyy = pbuffer.data(idx_op_si + 21);

    auto tr_0_yyyyyz = pbuffer.data(idx_op_si + 22);

    auto tr_0_yyyyzz = pbuffer.data(idx_op_si + 23);

    auto tr_0_yyyzzz = pbuffer.data(idx_op_si + 24);

    auto tr_0_yyzzzz = pbuffer.data(idx_op_si + 25);

    auto tr_0_yzzzzz = pbuffer.data(idx_op_si + 26);

    auto tr_0_zzzzzz = pbuffer.data(idx_op_si + 27);

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

    // Set up components of targeted buffer : SG

    auto tr_0_0_xx_0_xxxx = pbuffer.data(idx_op_geom_020_sg);

    auto tr_0_0_xx_0_xxxy = pbuffer.data(idx_op_geom_020_sg + 1);

    auto tr_0_0_xx_0_xxxz = pbuffer.data(idx_op_geom_020_sg + 2);

    auto tr_0_0_xx_0_xxyy = pbuffer.data(idx_op_geom_020_sg + 3);

    auto tr_0_0_xx_0_xxyz = pbuffer.data(idx_op_geom_020_sg + 4);

    auto tr_0_0_xx_0_xxzz = pbuffer.data(idx_op_geom_020_sg + 5);

    auto tr_0_0_xx_0_xyyy = pbuffer.data(idx_op_geom_020_sg + 6);

    auto tr_0_0_xx_0_xyyz = pbuffer.data(idx_op_geom_020_sg + 7);

    auto tr_0_0_xx_0_xyzz = pbuffer.data(idx_op_geom_020_sg + 8);

    auto tr_0_0_xx_0_xzzz = pbuffer.data(idx_op_geom_020_sg + 9);

    auto tr_0_0_xx_0_yyyy = pbuffer.data(idx_op_geom_020_sg + 10);

    auto tr_0_0_xx_0_yyyz = pbuffer.data(idx_op_geom_020_sg + 11);

    auto tr_0_0_xx_0_yyzz = pbuffer.data(idx_op_geom_020_sg + 12);

    auto tr_0_0_xx_0_yzzz = pbuffer.data(idx_op_geom_020_sg + 13);

    auto tr_0_0_xx_0_zzzz = pbuffer.data(idx_op_geom_020_sg + 14);

    auto tr_0_0_xy_0_xxxx = pbuffer.data(idx_op_geom_020_sg + 15);

    auto tr_0_0_xy_0_xxxy = pbuffer.data(idx_op_geom_020_sg + 16);

    auto tr_0_0_xy_0_xxxz = pbuffer.data(idx_op_geom_020_sg + 17);

    auto tr_0_0_xy_0_xxyy = pbuffer.data(idx_op_geom_020_sg + 18);

    auto tr_0_0_xy_0_xxyz = pbuffer.data(idx_op_geom_020_sg + 19);

    auto tr_0_0_xy_0_xxzz = pbuffer.data(idx_op_geom_020_sg + 20);

    auto tr_0_0_xy_0_xyyy = pbuffer.data(idx_op_geom_020_sg + 21);

    auto tr_0_0_xy_0_xyyz = pbuffer.data(idx_op_geom_020_sg + 22);

    auto tr_0_0_xy_0_xyzz = pbuffer.data(idx_op_geom_020_sg + 23);

    auto tr_0_0_xy_0_xzzz = pbuffer.data(idx_op_geom_020_sg + 24);

    auto tr_0_0_xy_0_yyyy = pbuffer.data(idx_op_geom_020_sg + 25);

    auto tr_0_0_xy_0_yyyz = pbuffer.data(idx_op_geom_020_sg + 26);

    auto tr_0_0_xy_0_yyzz = pbuffer.data(idx_op_geom_020_sg + 27);

    auto tr_0_0_xy_0_yzzz = pbuffer.data(idx_op_geom_020_sg + 28);

    auto tr_0_0_xy_0_zzzz = pbuffer.data(idx_op_geom_020_sg + 29);

    auto tr_0_0_xz_0_xxxx = pbuffer.data(idx_op_geom_020_sg + 30);

    auto tr_0_0_xz_0_xxxy = pbuffer.data(idx_op_geom_020_sg + 31);

    auto tr_0_0_xz_0_xxxz = pbuffer.data(idx_op_geom_020_sg + 32);

    auto tr_0_0_xz_0_xxyy = pbuffer.data(idx_op_geom_020_sg + 33);

    auto tr_0_0_xz_0_xxyz = pbuffer.data(idx_op_geom_020_sg + 34);

    auto tr_0_0_xz_0_xxzz = pbuffer.data(idx_op_geom_020_sg + 35);

    auto tr_0_0_xz_0_xyyy = pbuffer.data(idx_op_geom_020_sg + 36);

    auto tr_0_0_xz_0_xyyz = pbuffer.data(idx_op_geom_020_sg + 37);

    auto tr_0_0_xz_0_xyzz = pbuffer.data(idx_op_geom_020_sg + 38);

    auto tr_0_0_xz_0_xzzz = pbuffer.data(idx_op_geom_020_sg + 39);

    auto tr_0_0_xz_0_yyyy = pbuffer.data(idx_op_geom_020_sg + 40);

    auto tr_0_0_xz_0_yyyz = pbuffer.data(idx_op_geom_020_sg + 41);

    auto tr_0_0_xz_0_yyzz = pbuffer.data(idx_op_geom_020_sg + 42);

    auto tr_0_0_xz_0_yzzz = pbuffer.data(idx_op_geom_020_sg + 43);

    auto tr_0_0_xz_0_zzzz = pbuffer.data(idx_op_geom_020_sg + 44);

    auto tr_0_0_yy_0_xxxx = pbuffer.data(idx_op_geom_020_sg + 45);

    auto tr_0_0_yy_0_xxxy = pbuffer.data(idx_op_geom_020_sg + 46);

    auto tr_0_0_yy_0_xxxz = pbuffer.data(idx_op_geom_020_sg + 47);

    auto tr_0_0_yy_0_xxyy = pbuffer.data(idx_op_geom_020_sg + 48);

    auto tr_0_0_yy_0_xxyz = pbuffer.data(idx_op_geom_020_sg + 49);

    auto tr_0_0_yy_0_xxzz = pbuffer.data(idx_op_geom_020_sg + 50);

    auto tr_0_0_yy_0_xyyy = pbuffer.data(idx_op_geom_020_sg + 51);

    auto tr_0_0_yy_0_xyyz = pbuffer.data(idx_op_geom_020_sg + 52);

    auto tr_0_0_yy_0_xyzz = pbuffer.data(idx_op_geom_020_sg + 53);

    auto tr_0_0_yy_0_xzzz = pbuffer.data(idx_op_geom_020_sg + 54);

    auto tr_0_0_yy_0_yyyy = pbuffer.data(idx_op_geom_020_sg + 55);

    auto tr_0_0_yy_0_yyyz = pbuffer.data(idx_op_geom_020_sg + 56);

    auto tr_0_0_yy_0_yyzz = pbuffer.data(idx_op_geom_020_sg + 57);

    auto tr_0_0_yy_0_yzzz = pbuffer.data(idx_op_geom_020_sg + 58);

    auto tr_0_0_yy_0_zzzz = pbuffer.data(idx_op_geom_020_sg + 59);

    auto tr_0_0_yz_0_xxxx = pbuffer.data(idx_op_geom_020_sg + 60);

    auto tr_0_0_yz_0_xxxy = pbuffer.data(idx_op_geom_020_sg + 61);

    auto tr_0_0_yz_0_xxxz = pbuffer.data(idx_op_geom_020_sg + 62);

    auto tr_0_0_yz_0_xxyy = pbuffer.data(idx_op_geom_020_sg + 63);

    auto tr_0_0_yz_0_xxyz = pbuffer.data(idx_op_geom_020_sg + 64);

    auto tr_0_0_yz_0_xxzz = pbuffer.data(idx_op_geom_020_sg + 65);

    auto tr_0_0_yz_0_xyyy = pbuffer.data(idx_op_geom_020_sg + 66);

    auto tr_0_0_yz_0_xyyz = pbuffer.data(idx_op_geom_020_sg + 67);

    auto tr_0_0_yz_0_xyzz = pbuffer.data(idx_op_geom_020_sg + 68);

    auto tr_0_0_yz_0_xzzz = pbuffer.data(idx_op_geom_020_sg + 69);

    auto tr_0_0_yz_0_yyyy = pbuffer.data(idx_op_geom_020_sg + 70);

    auto tr_0_0_yz_0_yyyz = pbuffer.data(idx_op_geom_020_sg + 71);

    auto tr_0_0_yz_0_yyzz = pbuffer.data(idx_op_geom_020_sg + 72);

    auto tr_0_0_yz_0_yzzz = pbuffer.data(idx_op_geom_020_sg + 73);

    auto tr_0_0_yz_0_zzzz = pbuffer.data(idx_op_geom_020_sg + 74);

    auto tr_0_0_zz_0_xxxx = pbuffer.data(idx_op_geom_020_sg + 75);

    auto tr_0_0_zz_0_xxxy = pbuffer.data(idx_op_geom_020_sg + 76);

    auto tr_0_0_zz_0_xxxz = pbuffer.data(idx_op_geom_020_sg + 77);

    auto tr_0_0_zz_0_xxyy = pbuffer.data(idx_op_geom_020_sg + 78);

    auto tr_0_0_zz_0_xxyz = pbuffer.data(idx_op_geom_020_sg + 79);

    auto tr_0_0_zz_0_xxzz = pbuffer.data(idx_op_geom_020_sg + 80);

    auto tr_0_0_zz_0_xyyy = pbuffer.data(idx_op_geom_020_sg + 81);

    auto tr_0_0_zz_0_xyyz = pbuffer.data(idx_op_geom_020_sg + 82);

    auto tr_0_0_zz_0_xyzz = pbuffer.data(idx_op_geom_020_sg + 83);

    auto tr_0_0_zz_0_xzzz = pbuffer.data(idx_op_geom_020_sg + 84);

    auto tr_0_0_zz_0_yyyy = pbuffer.data(idx_op_geom_020_sg + 85);

    auto tr_0_0_zz_0_yyyz = pbuffer.data(idx_op_geom_020_sg + 86);

    auto tr_0_0_zz_0_yyzz = pbuffer.data(idx_op_geom_020_sg + 87);

    auto tr_0_0_zz_0_yzzz = pbuffer.data(idx_op_geom_020_sg + 88);

    auto tr_0_0_zz_0_zzzz = pbuffer.data(idx_op_geom_020_sg + 89);

    #pragma omp simd aligned(tr_0_0_xx_0_xxxx, tr_0_0_xx_0_xxxy, tr_0_0_xx_0_xxxz, tr_0_0_xx_0_xxyy, tr_0_0_xx_0_xxyz, tr_0_0_xx_0_xxzz, tr_0_0_xx_0_xyyy, tr_0_0_xx_0_xyyz, tr_0_0_xx_0_xyzz, tr_0_0_xx_0_xzzz, tr_0_0_xx_0_yyyy, tr_0_0_xx_0_yyyz, tr_0_0_xx_0_yyzz, tr_0_0_xx_0_yzzz, tr_0_0_xx_0_zzzz, tr_0_0_xy_0_xxxx, tr_0_0_xy_0_xxxy, tr_0_0_xy_0_xxxz, tr_0_0_xy_0_xxyy, tr_0_0_xy_0_xxyz, tr_0_0_xy_0_xxzz, tr_0_0_xy_0_xyyy, tr_0_0_xy_0_xyyz, tr_0_0_xy_0_xyzz, tr_0_0_xy_0_xzzz, tr_0_0_xy_0_yyyy, tr_0_0_xy_0_yyyz, tr_0_0_xy_0_yyzz, tr_0_0_xy_0_yzzz, tr_0_0_xy_0_zzzz, tr_0_0_xz_0_xxxx, tr_0_0_xz_0_xxxy, tr_0_0_xz_0_xxxz, tr_0_0_xz_0_xxyy, tr_0_0_xz_0_xxyz, tr_0_0_xz_0_xxzz, tr_0_0_xz_0_xyyy, tr_0_0_xz_0_xyyz, tr_0_0_xz_0_xyzz, tr_0_0_xz_0_xzzz, tr_0_0_xz_0_yyyy, tr_0_0_xz_0_yyyz, tr_0_0_xz_0_yyzz, tr_0_0_xz_0_yzzz, tr_0_0_xz_0_zzzz, tr_0_0_yy_0_xxxx, tr_0_0_yy_0_xxxy, tr_0_0_yy_0_xxxz, tr_0_0_yy_0_xxyy, tr_0_0_yy_0_xxyz, tr_0_0_yy_0_xxzz, tr_0_0_yy_0_xyyy, tr_0_0_yy_0_xyyz, tr_0_0_yy_0_xyzz, tr_0_0_yy_0_xzzz, tr_0_0_yy_0_yyyy, tr_0_0_yy_0_yyyz, tr_0_0_yy_0_yyzz, tr_0_0_yy_0_yzzz, tr_0_0_yy_0_zzzz, tr_0_0_yz_0_xxxx, tr_0_0_yz_0_xxxy, tr_0_0_yz_0_xxxz, tr_0_0_yz_0_xxyy, tr_0_0_yz_0_xxyz, tr_0_0_yz_0_xxzz, tr_0_0_yz_0_xyyy, tr_0_0_yz_0_xyyz, tr_0_0_yz_0_xyzz, tr_0_0_yz_0_xzzz, tr_0_0_yz_0_yyyy, tr_0_0_yz_0_yyyz, tr_0_0_yz_0_yyzz, tr_0_0_yz_0_yzzz, tr_0_0_yz_0_zzzz, tr_0_0_zz_0_xxxx, tr_0_0_zz_0_xxxy, tr_0_0_zz_0_xxxz, tr_0_0_zz_0_xxyy, tr_0_0_zz_0_xxyz, tr_0_0_zz_0_xxzz, tr_0_0_zz_0_xyyy, tr_0_0_zz_0_xyyz, tr_0_0_zz_0_xyzz, tr_0_0_zz_0_xzzz, tr_0_0_zz_0_yyyy, tr_0_0_zz_0_yyyz, tr_0_0_zz_0_yyzz, tr_0_0_zz_0_yzzz, tr_0_0_zz_0_zzzz, tr_0_xx, tr_0_xxxx, tr_0_xxxxxx, tr_0_xxxxxy, tr_0_xxxxxz, tr_0_xxxxyy, tr_0_xxxxyz, tr_0_xxxxzz, tr_0_xxxy, tr_0_xxxyyy, tr_0_xxxyyz, tr_0_xxxyzz, tr_0_xxxz, tr_0_xxxzzz, tr_0_xxyy, tr_0_xxyyyy, tr_0_xxyyyz, tr_0_xxyyzz, tr_0_xxyz, tr_0_xxyzzz, tr_0_xxzz, tr_0_xxzzzz, tr_0_xy, tr_0_xyyy, tr_0_xyyyyy, tr_0_xyyyyz, tr_0_xyyyzz, tr_0_xyyz, tr_0_xyyzzz, tr_0_xyzz, tr_0_xyzzzz, tr_0_xz, tr_0_xzzz, tr_0_xzzzzz, tr_0_yy, tr_0_yyyy, tr_0_yyyyyy, tr_0_yyyyyz, tr_0_yyyyzz, tr_0_yyyz, tr_0_yyyzzz, tr_0_yyzz, tr_0_yyzzzz, tr_0_yz, tr_0_yzzz, tr_0_yzzzzz, tr_0_zz, tr_0_zzzz, tr_0_zzzzzz, tr_x_xxx, tr_x_xxxxx, tr_x_xxxxy, tr_x_xxxxz, tr_x_xxxyy, tr_x_xxxyz, tr_x_xxxzz, tr_x_xxy, tr_x_xxyyy, tr_x_xxyyz, tr_x_xxyzz, tr_x_xxz, tr_x_xxzzz, tr_x_xyy, tr_x_xyyyy, tr_x_xyyyz, tr_x_xyyzz, tr_x_xyz, tr_x_xyzzz, tr_x_xzz, tr_x_xzzzz, tr_x_yyy, tr_x_yyyyy, tr_x_yyyyz, tr_x_yyyzz, tr_x_yyz, tr_x_yyzzz, tr_x_yzz, tr_x_yzzzz, tr_x_zzz, tr_x_zzzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_y_xxx, tr_y_xxxxx, tr_y_xxxxy, tr_y_xxxxz, tr_y_xxxyy, tr_y_xxxyz, tr_y_xxxzz, tr_y_xxy, tr_y_xxyyy, tr_y_xxyyz, tr_y_xxyzz, tr_y_xxz, tr_y_xxzzz, tr_y_xyy, tr_y_xyyyy, tr_y_xyyyz, tr_y_xyyzz, tr_y_xyz, tr_y_xyzzz, tr_y_xzz, tr_y_xzzzz, tr_y_yyy, tr_y_yyyyy, tr_y_yyyyz, tr_y_yyyzz, tr_y_yyz, tr_y_yyzzz, tr_y_yzz, tr_y_yzzzz, tr_y_zzz, tr_y_zzzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_z_xxx, tr_z_xxxxx, tr_z_xxxxy, tr_z_xxxxz, tr_z_xxxyy, tr_z_xxxyz, tr_z_xxxzz, tr_z_xxy, tr_z_xxyyy, tr_z_xxyyz, tr_z_xxyzz, tr_z_xxz, tr_z_xxzzz, tr_z_xyy, tr_z_xyyyy, tr_z_xyyyz, tr_z_xyyzz, tr_z_xyz, tr_z_xyzzz, tr_z_xzz, tr_z_xzzzz, tr_z_yyy, tr_z_yyyyy, tr_z_yyyyz, tr_z_yyyzz, tr_z_yyz, tr_z_yyzzz, tr_z_yzz, tr_z_yzzzz, tr_z_zzz, tr_z_zzzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_0_xxxx[i] = 12.0 * tr_0_xx[i] - 2.0 * tr_0_xxxx[i] * tbe_0 - 18.0 * tr_0_xxxx[i] * tke_0 + 4.0 * tr_0_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_x_xxx[i] * tbe_0 + 8.0 * tr_x_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xxxy[i] = 6.0 * tr_0_xy[i] - 2.0 * tr_0_xxxy[i] * tbe_0 - 14.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_x_xxy[i] * tbe_0 + 8.0 * tr_x_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xxxz[i] = 6.0 * tr_0_xz[i] - 2.0 * tr_0_xxxz[i] * tbe_0 - 14.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_x_xxz[i] * tbe_0 + 8.0 * tr_x_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xxyy[i] = 2.0 * tr_0_yy[i] - 2.0 * tr_0_xxyy[i] * tbe_0 - 10.0 * tr_0_xxyy[i] * tke_0 + 4.0 * tr_0_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_x_xyy[i] * tbe_0 + 8.0 * tr_x_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xxyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_xxyz[i] * tbe_0 - 10.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_x_xyz[i] * tbe_0 + 8.0 * tr_x_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xxzz[i] = 2.0 * tr_0_zz[i] - 2.0 * tr_0_xxzz[i] * tbe_0 - 10.0 * tr_0_xxzz[i] * tke_0 + 4.0 * tr_0_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_x_xzz[i] * tbe_0 + 8.0 * tr_x_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xyyy[i] = -2.0 * tr_0_xyyy[i] * tbe_0 - 6.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_0_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_x_yyy[i] * tbe_0 + 8.0 * tr_x_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xyyz[i] = -2.0 * tr_0_xyyz[i] * tbe_0 - 6.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_x_yyz[i] * tbe_0 + 8.0 * tr_x_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xyzz[i] = -2.0 * tr_0_xyzz[i] * tbe_0 - 6.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_x_yzz[i] * tbe_0 + 8.0 * tr_x_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_xzzz[i] = -2.0 * tr_0_xzzz[i] * tbe_0 - 6.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_0_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_x_zzz[i] * tbe_0 + 8.0 * tr_x_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_yyyy[i] = -2.0 * tr_0_yyyy[i] * tbe_0 - 2.0 * tr_0_yyyy[i] * tke_0 + 4.0 * tr_0_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_x_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_yyyz[i] = -2.0 * tr_0_yyyz[i] * tbe_0 - 2.0 * tr_0_yyyz[i] * tke_0 + 4.0 * tr_0_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_x_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_yyzz[i] = -2.0 * tr_0_yyzz[i] * tbe_0 - 2.0 * tr_0_yyzz[i] * tke_0 + 4.0 * tr_0_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_x_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_yzzz[i] = -2.0 * tr_0_yzzz[i] * tbe_0 - 2.0 * tr_0_yzzz[i] * tke_0 + 4.0 * tr_0_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_x_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_0_zzzz[i] = -2.0 * tr_0_zzzz[i] * tbe_0 - 2.0 * tr_0_zzzz[i] * tke_0 + 4.0 * tr_0_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_x_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xx_zzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xxxx[i] = -8.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_y_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xxxy[i] = 3.0 * tr_0_xx[i] - 6.0 * tr_0_xxyy[i] * tke_0 - 2.0 * tr_0_xxxx[i] * tke_0 + 4.0 * tr_0_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_y_xxy[i] * tbe_0 + 4.0 * tr_y_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_x_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xxxz[i] = -6.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_y_xxz[i] * tbe_0 + 4.0 * tr_y_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xxyy[i] = 4.0 * tr_0_xy[i] - 4.0 * tr_0_xyyy[i] * tke_0 - 4.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_y_xyy[i] * tbe_0 + 4.0 * tr_y_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_x_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xxyz[i] = 2.0 * tr_0_xz[i] - 4.0 * tr_0_xyyz[i] * tke_0 - 2.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_y_xyz[i] * tbe_0 + 4.0 * tr_y_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_x_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xxzz[i] = -4.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_y_xzz[i] * tbe_0 + 4.0 * tr_y_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xyyy[i] = 3.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tke_0 - 6.0 * tr_0_xxyy[i] * tke_0 + 4.0 * tr_0_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_y_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_x_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xyyz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tke_0 - 4.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_y_yyz[i] * tbe_0 + 4.0 * tr_y_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_x_xyz[i] * tbe_0 + 4.0 * tr_x_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xyzz[i] = tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_0_xxzz[i] * tke_0 + 4.0 * tr_0_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_y_yzz[i] * tbe_0 + 4.0 * tr_y_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_x_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_xzzz[i] = -2.0 * tr_0_yzzz[i] * tke_0 + 4.0 * tr_0_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_y_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_x_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_yyyy[i] = -8.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_0_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_y_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_x_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_yyyz[i] = -6.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_y_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_x_yyz[i] * tbe_0 + 4.0 * tr_x_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_yyzz[i] = -4.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_y_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_x_yzz[i] * tbe_0 + 4.0 * tr_x_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_yzzz[i] = -2.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_0_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_y_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_x_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_0_zzzz[i] = 4.0 * tr_0_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_y_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_x_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xy_zzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xxxx[i] = -8.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_z_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xxxy[i] = -6.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_z_xxy[i] * tbe_0 + 4.0 * tr_z_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xxxz[i] = 3.0 * tr_0_xx[i] - 6.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_0_xxxx[i] * tke_0 + 4.0 * tr_0_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_z_xxz[i] * tbe_0 + 4.0 * tr_z_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxx[i] * tbe_0 + 4.0 * tr_x_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xxyy[i] = -4.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_z_xyy[i] * tbe_0 + 4.0 * tr_z_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_x_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xxyz[i] = 2.0 * tr_0_xy[i] - 4.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_z_xyz[i] * tbe_0 + 4.0 * tr_z_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xxy[i] * tbe_0 + 4.0 * tr_x_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xxzz[i] = 4.0 * tr_0_xz[i] - 4.0 * tr_0_xzzz[i] * tke_0 - 4.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_xzz[i] * tbe_0 + 4.0 * tr_z_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_x_xxz[i] * tbe_0 + 4.0 * tr_x_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xyyy[i] = -2.0 * tr_0_yyyz[i] * tke_0 + 4.0 * tr_0_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_z_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_x_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xyyz[i] = tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_0_xxyy[i] * tke_0 + 4.0 * tr_0_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_z_yyz[i] * tbe_0 + 4.0 * tr_z_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_xyy[i] * tbe_0 + 4.0 * tr_x_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xyzz[i] = 2.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tke_0 - 4.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_z_yzz[i] * tbe_0 + 4.0 * tr_z_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_x_xyz[i] * tbe_0 + 4.0 * tr_x_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_xzzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 6.0 * tr_0_xxzz[i] * tke_0 + 4.0 * tr_0_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_z_zzz[i] * tbe_0 + 4.0 * tr_z_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_x_xzz[i] * tbe_0 + 4.0 * tr_x_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_yyyy[i] = 4.0 * tr_0_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_x_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_yyyz[i] = -2.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_0_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_x_yyy[i] * tbe_0 + 4.0 * tr_x_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_yyzz[i] = -4.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_x_yyz[i] * tbe_0 + 4.0 * tr_x_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_yzzz[i] = -6.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_x_yzz[i] * tbe_0 + 4.0 * tr_x_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_0_zzzz[i] = -8.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_0_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_x_zzz[i] * tbe_0 + 4.0 * tr_x_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xz_zzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xxxx[i] = -2.0 * tr_0_xxxx[i] * tbe_0 - 2.0 * tr_0_xxxx[i] * tke_0 + 4.0 * tr_0_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_y_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xxxy[i] = -2.0 * tr_0_xxxy[i] * tbe_0 - 6.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_y_xxx[i] * tbe_0 + 8.0 * tr_y_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xxxz[i] = -2.0 * tr_0_xxxz[i] * tbe_0 - 2.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_y_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xxyy[i] = 2.0 * tr_0_xx[i] - 2.0 * tr_0_xxyy[i] * tbe_0 - 10.0 * tr_0_xxyy[i] * tke_0 + 4.0 * tr_0_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_y_xxy[i] * tbe_0 + 8.0 * tr_y_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xxyz[i] = -2.0 * tr_0_xxyz[i] * tbe_0 - 6.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_y_xxz[i] * tbe_0 + 8.0 * tr_y_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xxzz[i] = -2.0 * tr_0_xxzz[i] * tbe_0 - 2.0 * tr_0_xxzz[i] * tke_0 + 4.0 * tr_0_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_y_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xyyy[i] = 6.0 * tr_0_xy[i] - 2.0 * tr_0_xyyy[i] * tbe_0 - 14.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_0_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_y_xyy[i] * tbe_0 + 8.0 * tr_y_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xyyz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xyyz[i] * tbe_0 - 10.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_y_xyz[i] * tbe_0 + 8.0 * tr_y_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xyzz[i] = -2.0 * tr_0_xyzz[i] * tbe_0 - 6.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_y_xzz[i] * tbe_0 + 8.0 * tr_y_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_xzzz[i] = -2.0 * tr_0_xzzz[i] * tbe_0 - 2.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_0_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_y_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_yyyy[i] = 12.0 * tr_0_yy[i] - 2.0 * tr_0_yyyy[i] * tbe_0 - 18.0 * tr_0_yyyy[i] * tke_0 + 4.0 * tr_0_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_y_yyy[i] * tbe_0 + 8.0 * tr_y_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_yyyz[i] = 6.0 * tr_0_yz[i] - 2.0 * tr_0_yyyz[i] * tbe_0 - 14.0 * tr_0_yyyz[i] * tke_0 + 4.0 * tr_0_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_y_yyz[i] * tbe_0 + 8.0 * tr_y_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_yyzz[i] = 2.0 * tr_0_zz[i] - 2.0 * tr_0_yyzz[i] * tbe_0 - 10.0 * tr_0_yyzz[i] * tke_0 + 4.0 * tr_0_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_y_yzz[i] * tbe_0 + 8.0 * tr_y_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_yzzz[i] = -2.0 * tr_0_yzzz[i] * tbe_0 - 6.0 * tr_0_yzzz[i] * tke_0 + 4.0 * tr_0_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_y_zzz[i] * tbe_0 + 8.0 * tr_y_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_0_zzzz[i] = -2.0 * tr_0_zzzz[i] * tbe_0 - 2.0 * tr_0_zzzz[i] * tke_0 + 4.0 * tr_0_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_y_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yy_zzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xxxx[i] = 4.0 * tr_0_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_z_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_y_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xxxy[i] = -2.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_z_xxx[i] * tbe_0 + 4.0 * tr_z_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_y_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xxxz[i] = -2.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxx[i] * tbe_0 + 4.0 * tr_y_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xxyy[i] = -4.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_z_xxy[i] * tbe_0 + 4.0 * tr_z_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_y_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xxyz[i] = tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tke_0 - 2.0 * tr_0_xxyy[i] * tke_0 + 4.0 * tr_0_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_z_xxz[i] * tbe_0 + 4.0 * tr_z_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xxy[i] * tbe_0 + 4.0 * tr_y_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xxzz[i] = -4.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_y_xxz[i] * tbe_0 + 4.0 * tr_y_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xyyy[i] = -6.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_z_xyy[i] * tbe_0 + 4.0 * tr_z_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_y_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xyyz[i] = 2.0 * tr_0_xy[i] - 4.0 * tr_0_xyzz[i] * tke_0 - 2.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_0_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_z_xyz[i] * tbe_0 + 4.0 * tr_z_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_xyy[i] * tbe_0 + 4.0 * tr_y_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xyzz[i] = 2.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tke_0 - 4.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_z_xzz[i] * tbe_0 + 4.0 * tr_z_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_y_xyz[i] * tbe_0 + 4.0 * tr_y_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_xzzz[i] = -6.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_y_xzz[i] * tbe_0 + 4.0 * tr_y_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_yyyy[i] = -8.0 * tr_0_yyyz[i] * tke_0 + 4.0 * tr_0_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_z_yyy[i] * tbe_0 + 4.0 * tr_z_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_y_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_yyyz[i] = 3.0 * tr_0_yy[i] - 6.0 * tr_0_yyzz[i] * tke_0 - 2.0 * tr_0_yyyy[i] * tke_0 + 4.0 * tr_0_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_z_yyz[i] * tbe_0 + 4.0 * tr_z_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_y_yyy[i] * tbe_0 + 4.0 * tr_y_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_yyzz[i] = 4.0 * tr_0_yz[i] - 4.0 * tr_0_yzzz[i] * tke_0 - 4.0 * tr_0_yyyz[i] * tke_0 + 4.0 * tr_0_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_yzz[i] * tbe_0 + 4.0 * tr_z_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_y_yyz[i] * tbe_0 + 4.0 * tr_y_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_yzzz[i] = 3.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tke_0 - 6.0 * tr_0_yyzz[i] * tke_0 + 4.0 * tr_0_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_z_zzz[i] * tbe_0 + 4.0 * tr_z_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_y_yzz[i] * tbe_0 + 4.0 * tr_y_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_0_zzzz[i] = -8.0 * tr_0_yzzz[i] * tke_0 + 4.0 * tr_0_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_z_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_y_zzz[i] * tbe_0 + 4.0 * tr_y_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yz_zzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xxxx[i] = -2.0 * tr_0_xxxx[i] * tbe_0 - 2.0 * tr_0_xxxx[i] * tke_0 + 4.0 * tr_0_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_z_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xxxy[i] = -2.0 * tr_0_xxxy[i] * tbe_0 - 2.0 * tr_0_xxxy[i] * tke_0 + 4.0 * tr_0_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_z_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xxxz[i] = -2.0 * tr_0_xxxz[i] * tbe_0 - 6.0 * tr_0_xxxz[i] * tke_0 + 4.0 * tr_0_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_xxx[i] * tbe_0 + 8.0 * tr_z_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xxyy[i] = -2.0 * tr_0_xxyy[i] * tbe_0 - 2.0 * tr_0_xxyy[i] * tke_0 + 4.0 * tr_0_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_z_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xxyz[i] = -2.0 * tr_0_xxyz[i] * tbe_0 - 6.0 * tr_0_xxyz[i] * tke_0 + 4.0 * tr_0_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_xxy[i] * tbe_0 + 8.0 * tr_z_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xxzz[i] = 2.0 * tr_0_xx[i] - 2.0 * tr_0_xxzz[i] * tbe_0 - 10.0 * tr_0_xxzz[i] * tke_0 + 4.0 * tr_0_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_z_xxz[i] * tbe_0 + 8.0 * tr_z_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xyyy[i] = -2.0 * tr_0_xyyy[i] * tbe_0 - 2.0 * tr_0_xyyy[i] * tke_0 + 4.0 * tr_0_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_z_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xyyz[i] = -2.0 * tr_0_xyyz[i] * tbe_0 - 6.0 * tr_0_xyyz[i] * tke_0 + 4.0 * tr_0_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_xyy[i] * tbe_0 + 8.0 * tr_z_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xyzz[i] = 2.0 * tr_0_xy[i] - 2.0 * tr_0_xyzz[i] * tbe_0 - 10.0 * tr_0_xyzz[i] * tke_0 + 4.0 * tr_0_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_z_xyz[i] * tbe_0 + 8.0 * tr_z_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_xzzz[i] = 6.0 * tr_0_xz[i] - 2.0 * tr_0_xzzz[i] * tbe_0 - 14.0 * tr_0_xzzz[i] * tke_0 + 4.0 * tr_0_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_z_xzz[i] * tbe_0 + 8.0 * tr_z_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_yyyy[i] = -2.0 * tr_0_yyyy[i] * tbe_0 - 2.0 * tr_0_yyyy[i] * tke_0 + 4.0 * tr_0_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_z_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_yyyz[i] = -2.0 * tr_0_yyyz[i] * tbe_0 - 6.0 * tr_0_yyyz[i] * tke_0 + 4.0 * tr_0_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_z_yyy[i] * tbe_0 + 8.0 * tr_z_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_yyzz[i] = 2.0 * tr_0_yy[i] - 2.0 * tr_0_yyzz[i] * tbe_0 - 10.0 * tr_0_yyzz[i] * tke_0 + 4.0 * tr_0_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_z_yyz[i] * tbe_0 + 8.0 * tr_z_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_yzzz[i] = 6.0 * tr_0_yz[i] - 2.0 * tr_0_yzzz[i] * tbe_0 - 14.0 * tr_0_yzzz[i] * tke_0 + 4.0 * tr_0_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_z_yzz[i] * tbe_0 + 8.0 * tr_z_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_0_zzzz[i] = 12.0 * tr_0_zz[i] - 2.0 * tr_0_zzzz[i] * tbe_0 - 18.0 * tr_0_zzzz[i] * tke_0 + 4.0 * tr_0_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_z_zzz[i] * tbe_0 + 8.0 * tr_z_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zz_zzzz[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

