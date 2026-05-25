#include "GeometricalDerivatives020ForGS.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_gs(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_gs,
                         const int idx_op_ds,
                         const int idx_op_fp,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_hp,
                         const int idx_op_is,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up components of auxiliary buffer : DS

    auto tr_xx_0 = pbuffer.data(idx_op_ds);

    auto tr_xy_0 = pbuffer.data(idx_op_ds + 1);

    auto tr_xz_0 = pbuffer.data(idx_op_ds + 2);

    auto tr_yy_0 = pbuffer.data(idx_op_ds + 3);

    auto tr_yz_0 = pbuffer.data(idx_op_ds + 4);

    auto tr_zz_0 = pbuffer.data(idx_op_ds + 5);

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

    // Set up components of auxiliary buffer : IS

    auto tr_xxxxxx_0 = pbuffer.data(idx_op_is);

    auto tr_xxxxxy_0 = pbuffer.data(idx_op_is + 1);

    auto tr_xxxxxz_0 = pbuffer.data(idx_op_is + 2);

    auto tr_xxxxyy_0 = pbuffer.data(idx_op_is + 3);

    auto tr_xxxxyz_0 = pbuffer.data(idx_op_is + 4);

    auto tr_xxxxzz_0 = pbuffer.data(idx_op_is + 5);

    auto tr_xxxyyy_0 = pbuffer.data(idx_op_is + 6);

    auto tr_xxxyyz_0 = pbuffer.data(idx_op_is + 7);

    auto tr_xxxyzz_0 = pbuffer.data(idx_op_is + 8);

    auto tr_xxxzzz_0 = pbuffer.data(idx_op_is + 9);

    auto tr_xxyyyy_0 = pbuffer.data(idx_op_is + 10);

    auto tr_xxyyyz_0 = pbuffer.data(idx_op_is + 11);

    auto tr_xxyyzz_0 = pbuffer.data(idx_op_is + 12);

    auto tr_xxyzzz_0 = pbuffer.data(idx_op_is + 13);

    auto tr_xxzzzz_0 = pbuffer.data(idx_op_is + 14);

    auto tr_xyyyyy_0 = pbuffer.data(idx_op_is + 15);

    auto tr_xyyyyz_0 = pbuffer.data(idx_op_is + 16);

    auto tr_xyyyzz_0 = pbuffer.data(idx_op_is + 17);

    auto tr_xyyzzz_0 = pbuffer.data(idx_op_is + 18);

    auto tr_xyzzzz_0 = pbuffer.data(idx_op_is + 19);

    auto tr_xzzzzz_0 = pbuffer.data(idx_op_is + 20);

    auto tr_yyyyyy_0 = pbuffer.data(idx_op_is + 21);

    auto tr_yyyyyz_0 = pbuffer.data(idx_op_is + 22);

    auto tr_yyyyzz_0 = pbuffer.data(idx_op_is + 23);

    auto tr_yyyzzz_0 = pbuffer.data(idx_op_is + 24);

    auto tr_yyzzzz_0 = pbuffer.data(idx_op_is + 25);

    auto tr_yzzzzz_0 = pbuffer.data(idx_op_is + 26);

    auto tr_zzzzzz_0 = pbuffer.data(idx_op_is + 27);

    // Set up components of targeted buffer : GS

    auto tr_0_0_xx_xxxx_0 = pbuffer.data(idx_op_geom_020_gs);

    auto tr_0_0_xx_xxxy_0 = pbuffer.data(idx_op_geom_020_gs + 1);

    auto tr_0_0_xx_xxxz_0 = pbuffer.data(idx_op_geom_020_gs + 2);

    auto tr_0_0_xx_xxyy_0 = pbuffer.data(idx_op_geom_020_gs + 3);

    auto tr_0_0_xx_xxyz_0 = pbuffer.data(idx_op_geom_020_gs + 4);

    auto tr_0_0_xx_xxzz_0 = pbuffer.data(idx_op_geom_020_gs + 5);

    auto tr_0_0_xx_xyyy_0 = pbuffer.data(idx_op_geom_020_gs + 6);

    auto tr_0_0_xx_xyyz_0 = pbuffer.data(idx_op_geom_020_gs + 7);

    auto tr_0_0_xx_xyzz_0 = pbuffer.data(idx_op_geom_020_gs + 8);

    auto tr_0_0_xx_xzzz_0 = pbuffer.data(idx_op_geom_020_gs + 9);

    auto tr_0_0_xx_yyyy_0 = pbuffer.data(idx_op_geom_020_gs + 10);

    auto tr_0_0_xx_yyyz_0 = pbuffer.data(idx_op_geom_020_gs + 11);

    auto tr_0_0_xx_yyzz_0 = pbuffer.data(idx_op_geom_020_gs + 12);

    auto tr_0_0_xx_yzzz_0 = pbuffer.data(idx_op_geom_020_gs + 13);

    auto tr_0_0_xx_zzzz_0 = pbuffer.data(idx_op_geom_020_gs + 14);

    auto tr_0_0_xy_xxxx_0 = pbuffer.data(idx_op_geom_020_gs + 15);

    auto tr_0_0_xy_xxxy_0 = pbuffer.data(idx_op_geom_020_gs + 16);

    auto tr_0_0_xy_xxxz_0 = pbuffer.data(idx_op_geom_020_gs + 17);

    auto tr_0_0_xy_xxyy_0 = pbuffer.data(idx_op_geom_020_gs + 18);

    auto tr_0_0_xy_xxyz_0 = pbuffer.data(idx_op_geom_020_gs + 19);

    auto tr_0_0_xy_xxzz_0 = pbuffer.data(idx_op_geom_020_gs + 20);

    auto tr_0_0_xy_xyyy_0 = pbuffer.data(idx_op_geom_020_gs + 21);

    auto tr_0_0_xy_xyyz_0 = pbuffer.data(idx_op_geom_020_gs + 22);

    auto tr_0_0_xy_xyzz_0 = pbuffer.data(idx_op_geom_020_gs + 23);

    auto tr_0_0_xy_xzzz_0 = pbuffer.data(idx_op_geom_020_gs + 24);

    auto tr_0_0_xy_yyyy_0 = pbuffer.data(idx_op_geom_020_gs + 25);

    auto tr_0_0_xy_yyyz_0 = pbuffer.data(idx_op_geom_020_gs + 26);

    auto tr_0_0_xy_yyzz_0 = pbuffer.data(idx_op_geom_020_gs + 27);

    auto tr_0_0_xy_yzzz_0 = pbuffer.data(idx_op_geom_020_gs + 28);

    auto tr_0_0_xy_zzzz_0 = pbuffer.data(idx_op_geom_020_gs + 29);

    auto tr_0_0_xz_xxxx_0 = pbuffer.data(idx_op_geom_020_gs + 30);

    auto tr_0_0_xz_xxxy_0 = pbuffer.data(idx_op_geom_020_gs + 31);

    auto tr_0_0_xz_xxxz_0 = pbuffer.data(idx_op_geom_020_gs + 32);

    auto tr_0_0_xz_xxyy_0 = pbuffer.data(idx_op_geom_020_gs + 33);

    auto tr_0_0_xz_xxyz_0 = pbuffer.data(idx_op_geom_020_gs + 34);

    auto tr_0_0_xz_xxzz_0 = pbuffer.data(idx_op_geom_020_gs + 35);

    auto tr_0_0_xz_xyyy_0 = pbuffer.data(idx_op_geom_020_gs + 36);

    auto tr_0_0_xz_xyyz_0 = pbuffer.data(idx_op_geom_020_gs + 37);

    auto tr_0_0_xz_xyzz_0 = pbuffer.data(idx_op_geom_020_gs + 38);

    auto tr_0_0_xz_xzzz_0 = pbuffer.data(idx_op_geom_020_gs + 39);

    auto tr_0_0_xz_yyyy_0 = pbuffer.data(idx_op_geom_020_gs + 40);

    auto tr_0_0_xz_yyyz_0 = pbuffer.data(idx_op_geom_020_gs + 41);

    auto tr_0_0_xz_yyzz_0 = pbuffer.data(idx_op_geom_020_gs + 42);

    auto tr_0_0_xz_yzzz_0 = pbuffer.data(idx_op_geom_020_gs + 43);

    auto tr_0_0_xz_zzzz_0 = pbuffer.data(idx_op_geom_020_gs + 44);

    auto tr_0_0_yy_xxxx_0 = pbuffer.data(idx_op_geom_020_gs + 45);

    auto tr_0_0_yy_xxxy_0 = pbuffer.data(idx_op_geom_020_gs + 46);

    auto tr_0_0_yy_xxxz_0 = pbuffer.data(idx_op_geom_020_gs + 47);

    auto tr_0_0_yy_xxyy_0 = pbuffer.data(idx_op_geom_020_gs + 48);

    auto tr_0_0_yy_xxyz_0 = pbuffer.data(idx_op_geom_020_gs + 49);

    auto tr_0_0_yy_xxzz_0 = pbuffer.data(idx_op_geom_020_gs + 50);

    auto tr_0_0_yy_xyyy_0 = pbuffer.data(idx_op_geom_020_gs + 51);

    auto tr_0_0_yy_xyyz_0 = pbuffer.data(idx_op_geom_020_gs + 52);

    auto tr_0_0_yy_xyzz_0 = pbuffer.data(idx_op_geom_020_gs + 53);

    auto tr_0_0_yy_xzzz_0 = pbuffer.data(idx_op_geom_020_gs + 54);

    auto tr_0_0_yy_yyyy_0 = pbuffer.data(idx_op_geom_020_gs + 55);

    auto tr_0_0_yy_yyyz_0 = pbuffer.data(idx_op_geom_020_gs + 56);

    auto tr_0_0_yy_yyzz_0 = pbuffer.data(idx_op_geom_020_gs + 57);

    auto tr_0_0_yy_yzzz_0 = pbuffer.data(idx_op_geom_020_gs + 58);

    auto tr_0_0_yy_zzzz_0 = pbuffer.data(idx_op_geom_020_gs + 59);

    auto tr_0_0_yz_xxxx_0 = pbuffer.data(idx_op_geom_020_gs + 60);

    auto tr_0_0_yz_xxxy_0 = pbuffer.data(idx_op_geom_020_gs + 61);

    auto tr_0_0_yz_xxxz_0 = pbuffer.data(idx_op_geom_020_gs + 62);

    auto tr_0_0_yz_xxyy_0 = pbuffer.data(idx_op_geom_020_gs + 63);

    auto tr_0_0_yz_xxyz_0 = pbuffer.data(idx_op_geom_020_gs + 64);

    auto tr_0_0_yz_xxzz_0 = pbuffer.data(idx_op_geom_020_gs + 65);

    auto tr_0_0_yz_xyyy_0 = pbuffer.data(idx_op_geom_020_gs + 66);

    auto tr_0_0_yz_xyyz_0 = pbuffer.data(idx_op_geom_020_gs + 67);

    auto tr_0_0_yz_xyzz_0 = pbuffer.data(idx_op_geom_020_gs + 68);

    auto tr_0_0_yz_xzzz_0 = pbuffer.data(idx_op_geom_020_gs + 69);

    auto tr_0_0_yz_yyyy_0 = pbuffer.data(idx_op_geom_020_gs + 70);

    auto tr_0_0_yz_yyyz_0 = pbuffer.data(idx_op_geom_020_gs + 71);

    auto tr_0_0_yz_yyzz_0 = pbuffer.data(idx_op_geom_020_gs + 72);

    auto tr_0_0_yz_yzzz_0 = pbuffer.data(idx_op_geom_020_gs + 73);

    auto tr_0_0_yz_zzzz_0 = pbuffer.data(idx_op_geom_020_gs + 74);

    auto tr_0_0_zz_xxxx_0 = pbuffer.data(idx_op_geom_020_gs + 75);

    auto tr_0_0_zz_xxxy_0 = pbuffer.data(idx_op_geom_020_gs + 76);

    auto tr_0_0_zz_xxxz_0 = pbuffer.data(idx_op_geom_020_gs + 77);

    auto tr_0_0_zz_xxyy_0 = pbuffer.data(idx_op_geom_020_gs + 78);

    auto tr_0_0_zz_xxyz_0 = pbuffer.data(idx_op_geom_020_gs + 79);

    auto tr_0_0_zz_xxzz_0 = pbuffer.data(idx_op_geom_020_gs + 80);

    auto tr_0_0_zz_xyyy_0 = pbuffer.data(idx_op_geom_020_gs + 81);

    auto tr_0_0_zz_xyyz_0 = pbuffer.data(idx_op_geom_020_gs + 82);

    auto tr_0_0_zz_xyzz_0 = pbuffer.data(idx_op_geom_020_gs + 83);

    auto tr_0_0_zz_xzzz_0 = pbuffer.data(idx_op_geom_020_gs + 84);

    auto tr_0_0_zz_yyyy_0 = pbuffer.data(idx_op_geom_020_gs + 85);

    auto tr_0_0_zz_yyyz_0 = pbuffer.data(idx_op_geom_020_gs + 86);

    auto tr_0_0_zz_yyzz_0 = pbuffer.data(idx_op_geom_020_gs + 87);

    auto tr_0_0_zz_yzzz_0 = pbuffer.data(idx_op_geom_020_gs + 88);

    auto tr_0_0_zz_zzzz_0 = pbuffer.data(idx_op_geom_020_gs + 89);

    #pragma omp simd aligned(tr_0_0_xx_xxxx_0, tr_0_0_xx_xxxy_0, tr_0_0_xx_xxxz_0, tr_0_0_xx_xxyy_0, tr_0_0_xx_xxyz_0, tr_0_0_xx_xxzz_0, tr_0_0_xx_xyyy_0, tr_0_0_xx_xyyz_0, tr_0_0_xx_xyzz_0, tr_0_0_xx_xzzz_0, tr_0_0_xx_yyyy_0, tr_0_0_xx_yyyz_0, tr_0_0_xx_yyzz_0, tr_0_0_xx_yzzz_0, tr_0_0_xx_zzzz_0, tr_0_0_xy_xxxx_0, tr_0_0_xy_xxxy_0, tr_0_0_xy_xxxz_0, tr_0_0_xy_xxyy_0, tr_0_0_xy_xxyz_0, tr_0_0_xy_xxzz_0, tr_0_0_xy_xyyy_0, tr_0_0_xy_xyyz_0, tr_0_0_xy_xyzz_0, tr_0_0_xy_xzzz_0, tr_0_0_xy_yyyy_0, tr_0_0_xy_yyyz_0, tr_0_0_xy_yyzz_0, tr_0_0_xy_yzzz_0, tr_0_0_xy_zzzz_0, tr_0_0_xz_xxxx_0, tr_0_0_xz_xxxy_0, tr_0_0_xz_xxxz_0, tr_0_0_xz_xxyy_0, tr_0_0_xz_xxyz_0, tr_0_0_xz_xxzz_0, tr_0_0_xz_xyyy_0, tr_0_0_xz_xyyz_0, tr_0_0_xz_xyzz_0, tr_0_0_xz_xzzz_0, tr_0_0_xz_yyyy_0, tr_0_0_xz_yyyz_0, tr_0_0_xz_yyzz_0, tr_0_0_xz_yzzz_0, tr_0_0_xz_zzzz_0, tr_0_0_yy_xxxx_0, tr_0_0_yy_xxxy_0, tr_0_0_yy_xxxz_0, tr_0_0_yy_xxyy_0, tr_0_0_yy_xxyz_0, tr_0_0_yy_xxzz_0, tr_0_0_yy_xyyy_0, tr_0_0_yy_xyyz_0, tr_0_0_yy_xyzz_0, tr_0_0_yy_xzzz_0, tr_0_0_yy_yyyy_0, tr_0_0_yy_yyyz_0, tr_0_0_yy_yyzz_0, tr_0_0_yy_yzzz_0, tr_0_0_yy_zzzz_0, tr_0_0_yz_xxxx_0, tr_0_0_yz_xxxy_0, tr_0_0_yz_xxxz_0, tr_0_0_yz_xxyy_0, tr_0_0_yz_xxyz_0, tr_0_0_yz_xxzz_0, tr_0_0_yz_xyyy_0, tr_0_0_yz_xyyz_0, tr_0_0_yz_xyzz_0, tr_0_0_yz_xzzz_0, tr_0_0_yz_yyyy_0, tr_0_0_yz_yyyz_0, tr_0_0_yz_yyzz_0, tr_0_0_yz_yzzz_0, tr_0_0_yz_zzzz_0, tr_0_0_zz_xxxx_0, tr_0_0_zz_xxxy_0, tr_0_0_zz_xxxz_0, tr_0_0_zz_xxyy_0, tr_0_0_zz_xxyz_0, tr_0_0_zz_xxzz_0, tr_0_0_zz_xyyy_0, tr_0_0_zz_xyyz_0, tr_0_0_zz_xyzz_0, tr_0_0_zz_xzzz_0, tr_0_0_zz_yyyy_0, tr_0_0_zz_yyyz_0, tr_0_0_zz_yyzz_0, tr_0_0_zz_yzzz_0, tr_0_0_zz_zzzz_0, tr_xx_0, tr_xxx_x, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxx_x, tr_xxxxx_y, tr_xxxxx_z, tr_xxxxxx_0, tr_xxxxxy_0, tr_xxxxxz_0, tr_xxxxy_x, tr_xxxxy_y, tr_xxxxy_z, tr_xxxxyy_0, tr_xxxxyz_0, tr_xxxxz_x, tr_xxxxz_y, tr_xxxxz_z, tr_xxxxzz_0, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyy_x, tr_xxxyy_y, tr_xxxyy_z, tr_xxxyyy_0, tr_xxxyyz_0, tr_xxxyz_x, tr_xxxyz_y, tr_xxxyz_z, tr_xxxyzz_0, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_y, tr_xxxzz_z, tr_xxxzzz_0, tr_xxy_x, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyy_x, tr_xxyyy_y, tr_xxyyy_z, tr_xxyyyy_0, tr_xxyyyz_0, tr_xxyyz_x, tr_xxyyz_y, tr_xxyyz_z, tr_xxyyzz_0, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_y, tr_xxyzz_z, tr_xxyzzz_0, tr_xxz_x, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_y, tr_xxzzz_z, tr_xxzzzz_0, tr_xy_0, tr_xyy_x, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyy_x, tr_xyyyy_y, tr_xyyyy_z, tr_xyyyyy_0, tr_xyyyyz_0, tr_xyyyz_x, tr_xyyyz_y, tr_xyyyz_z, tr_xyyyzz_0, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_y, tr_xyyzz_z, tr_xyyzzz_0, tr_xyz_x, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_y, tr_xyzzz_z, tr_xyzzzz_0, tr_xz_0, tr_xzz_x, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_y, tr_xzzzz_z, tr_xzzzzz_0, tr_yy_0, tr_yyy_x, tr_yyy_y, tr_yyy_z, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyy_x, tr_yyyyy_y, tr_yyyyy_z, tr_yyyyyy_0, tr_yyyyyz_0, tr_yyyyz_x, tr_yyyyz_y, tr_yyyyz_z, tr_yyyyzz_0, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_y, tr_yyyzz_z, tr_yyyzzz_0, tr_yyz_x, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_y, tr_yyzzz_z, tr_yyzzzz_0, tr_yz_0, tr_yzz_x, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_y, tr_yzzzz_z, tr_yzzzzz_0, tr_zz_0, tr_zzz_x, tr_zzz_y, tr_zzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, tr_zzzzz_x, tr_zzzzz_y, tr_zzzzz_z, tr_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxx_0[i] = 12.0 * tr_xx_0[i] - 16.0 * tr_xxx_x[i] * tke_0 - 18.0 * tr_xxxx_0[i] * tbe_0 - 2.0 * tr_xxxx_0[i] * tke_0 + 4.0 * tr_xxxx_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_0[i] = 6.0 * tr_xy_0[i] - 12.0 * tr_xxy_x[i] * tke_0 - 14.0 * tr_xxxy_0[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tke_0 + 4.0 * tr_xxxy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_0[i] = 6.0 * tr_xz_0[i] - 12.0 * tr_xxz_x[i] * tke_0 - 14.0 * tr_xxxz_0[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tke_0 + 4.0 * tr_xxxz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_0[i] = 2.0 * tr_yy_0[i] - 8.0 * tr_xyy_x[i] * tke_0 - 10.0 * tr_xxyy_0[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tke_0 + 4.0 * tr_xxyy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_0[i] = 2.0 * tr_yz_0[i] - 8.0 * tr_xyz_x[i] * tke_0 - 10.0 * tr_xxyz_0[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tke_0 + 4.0 * tr_xxyz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_0[i] = 2.0 * tr_zz_0[i] - 8.0 * tr_xzz_x[i] * tke_0 - 10.0 * tr_xxzz_0[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tke_0 + 4.0 * tr_xxzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_0[i] = -4.0 * tr_yyy_x[i] * tke_0 - 6.0 * tr_xyyy_0[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tke_0 + 4.0 * tr_xyyy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_0[i] = -4.0 * tr_yyz_x[i] * tke_0 - 6.0 * tr_xyyz_0[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tke_0 + 4.0 * tr_xyyz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_0[i] = -4.0 * tr_yzz_x[i] * tke_0 - 6.0 * tr_xyzz_0[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tke_0 + 4.0 * tr_xyzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_0[i] = -4.0 * tr_zzz_x[i] * tke_0 - 6.0 * tr_xzzz_0[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tke_0 + 4.0 * tr_xzzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_0[i] = -2.0 * tr_yyyy_0[i] * tbe_0 - 2.0 * tr_yyyy_0[i] * tke_0 + 4.0 * tr_yyyy_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_0[i] = -2.0 * tr_yyyz_0[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tke_0 + 4.0 * tr_yyyz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_0[i] = -2.0 * tr_yyzz_0[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tke_0 + 4.0 * tr_yyzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_0[i] = -2.0 * tr_yzzz_0[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tke_0 + 4.0 * tr_yzzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_0[i] = -2.0 * tr_zzzz_0[i] * tbe_0 - 2.0 * tr_zzzz_0[i] * tke_0 + 4.0 * tr_zzzz_xx[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_0[i] = -8.0 * tr_xxx_y[i] * tke_0 - 8.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxx_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_0[i] = 3.0 * tr_xx_0[i] - 6.0 * tr_xxy_y[i] * tke_0 - 6.0 * tr_xxyy_0[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxxy_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxxy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_0[i] = -6.0 * tr_xxz_y[i] * tke_0 - 6.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxxz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_0[i] = 4.0 * tr_xy_0[i] - 4.0 * tr_xyy_y[i] * tke_0 - 4.0 * tr_xyyy_0[i] * tbe_0 - 4.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxyy_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_x[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_0[i] = 2.0 * tr_xz_0[i] - 4.0 * tr_xyz_y[i] * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 - 2.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxyz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_0[i] = -4.0 * tr_xzz_y[i] * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xxzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_0[i] = 3.0 * tr_yy_0[i] - 2.0 * tr_yyy_y[i] * tke_0 - 2.0 * tr_yyyy_0[i] * tbe_0 - 6.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyyy_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_x[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_0[i] = 2.0 * tr_yz_0[i] - 2.0 * tr_yyz_y[i] * tke_0 - 2.0 * tr_yyyz_0[i] * tbe_0 - 4.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyyz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_0[i] = tr_zz_0[i] - 2.0 * tr_yzz_y[i] * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 - 2.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xyzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_0[i] = -2.0 * tr_zzz_y[i] * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_xzzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_0[i] = -8.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyyy_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_x[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_0[i] = -6.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyyz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_0[i] = -4.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yyzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_0[i] = -2.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_yzzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_0[i] = 4.0 * tr_zzzz_xy[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_0[i] = -8.0 * tr_xxx_z[i] * tke_0 - 8.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxx_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_0[i] = -6.0 * tr_xxy_z[i] * tke_0 - 6.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxxy_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_0[i] = 3.0 * tr_xx_0[i] - 6.0 * tr_xxz_z[i] * tke_0 - 6.0 * tr_xxzz_0[i] * tbe_0 - 2.0 * tr_xxx_x[i] * tke_0 + 4.0 * tr_xxxz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_0[i] * tbe_0 + 4.0 * tr_xxxxz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_0[i] = -4.0 * tr_xyy_z[i] * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xxyy_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_0[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xyz_z[i] * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 - 2.0 * tr_xxy_x[i] * tke_0 + 4.0 * tr_xxyz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_0[i] = 4.0 * tr_xz_0[i] - 4.0 * tr_xzz_z[i] * tke_0 - 4.0 * tr_xzzz_0[i] * tbe_0 - 4.0 * tr_xxz_x[i] * tke_0 + 4.0 * tr_xxzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_0[i] = -2.0 * tr_yyy_z[i] * tke_0 - 2.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_xyyy_xz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_0[i] = tr_yy_0[i] - 2.0 * tr_yyz_z[i] * tke_0 - 2.0 * tr_yyzz_0[i] * tbe_0 - 2.0 * tr_xyy_x[i] * tke_0 + 4.0 * tr_xyyz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_0[i] = 2.0 * tr_yz_0[i] - 2.0 * tr_yzz_z[i] * tke_0 - 2.0 * tr_yzzz_0[i] * tbe_0 - 4.0 * tr_xyz_x[i] * tke_0 + 4.0 * tr_xyzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_0[i] = 3.0 * tr_zz_0[i] - 2.0 * tr_zzz_z[i] * tke_0 - 2.0 * tr_zzzz_0[i] * tbe_0 - 6.0 * tr_xzz_x[i] * tke_0 + 4.0 * tr_xzzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_x[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_0[i] * tbe_0 + 4.0 * tr_xxzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_0[i] = 4.0 * tr_yyyy_xz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_x[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_0[i] = -2.0 * tr_yyy_x[i] * tke_0 + 4.0 * tr_yyyz_xz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_x[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_0[i] = -4.0 * tr_yyz_x[i] * tke_0 + 4.0 * tr_yyzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_x[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_0[i] = -6.0 * tr_yzz_x[i] * tke_0 + 4.0 * tr_yzzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_x[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_0[i] = -8.0 * tr_zzz_x[i] * tke_0 + 4.0 * tr_zzzz_xz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_x[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_0[i] * tbe_0 + 4.0 * tr_xzzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_0[i] = -2.0 * tr_xxxx_0[i] * tbe_0 - 2.0 * tr_xxxx_0[i] * tke_0 + 4.0 * tr_xxxx_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_0[i] = -4.0 * tr_xxx_y[i] * tke_0 - 6.0 * tr_xxxy_0[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tke_0 + 4.0 * tr_xxxy_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_0[i] = -2.0 * tr_xxxz_0[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tke_0 + 4.0 * tr_xxxz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_0[i] = 2.0 * tr_xx_0[i] - 8.0 * tr_xxy_y[i] * tke_0 - 10.0 * tr_xxyy_0[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tke_0 + 4.0 * tr_xxyy_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_0[i] = -4.0 * tr_xxz_y[i] * tke_0 - 6.0 * tr_xxyz_0[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tke_0 + 4.0 * tr_xxyz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_0[i] = -2.0 * tr_xxzz_0[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tke_0 + 4.0 * tr_xxzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_0[i] = 6.0 * tr_xy_0[i] - 12.0 * tr_xyy_y[i] * tke_0 - 14.0 * tr_xyyy_0[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tke_0 + 4.0 * tr_xyyy_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_0[i] = 2.0 * tr_xz_0[i] - 8.0 * tr_xyz_y[i] * tke_0 - 10.0 * tr_xyyz_0[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tke_0 + 4.0 * tr_xyyz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_0[i] = -4.0 * tr_xzz_y[i] * tke_0 - 6.0 * tr_xyzz_0[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tke_0 + 4.0 * tr_xyzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_0[i] = -2.0 * tr_xzzz_0[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tke_0 + 4.0 * tr_xzzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_0[i] = 12.0 * tr_yy_0[i] - 16.0 * tr_yyy_y[i] * tke_0 - 18.0 * tr_yyyy_0[i] * tbe_0 - 2.0 * tr_yyyy_0[i] * tke_0 + 4.0 * tr_yyyy_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_0[i] = 6.0 * tr_yz_0[i] - 12.0 * tr_yyz_y[i] * tke_0 - 14.0 * tr_yyyz_0[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tke_0 + 4.0 * tr_yyyz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_0[i] = 2.0 * tr_zz_0[i] - 8.0 * tr_yzz_y[i] * tke_0 - 10.0 * tr_yyzz_0[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tke_0 + 4.0 * tr_yyzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_0[i] = -4.0 * tr_zzz_y[i] * tke_0 - 6.0 * tr_yzzz_0[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tke_0 + 4.0 * tr_yzzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_0[i] = -2.0 * tr_zzzz_0[i] * tbe_0 - 2.0 * tr_zzzz_0[i] * tke_0 + 4.0 * tr_zzzz_yy[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_0[i] = 4.0 * tr_xxxx_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_0[i] = -2.0 * tr_xxx_z[i] * tke_0 - 2.0 * tr_xxxz_0[i] * tbe_0 + 4.0 * tr_xxxy_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_0[i] = -2.0 * tr_xxx_y[i] * tke_0 + 4.0 * tr_xxxz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_0[i] * tbe_0 + 4.0 * tr_xxxyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_0[i] = -4.0 * tr_xxy_z[i] * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_0[i] = tr_xx_0[i] - 2.0 * tr_xxz_z[i] * tke_0 - 2.0 * tr_xxzz_0[i] * tbe_0 - 2.0 * tr_xxy_y[i] * tke_0 + 4.0 * tr_xxyz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_0[i] * tbe_0 + 4.0 * tr_xxyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_0[i] = -4.0 * tr_xxz_y[i] * tke_0 + 4.0 * tr_xxzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_y[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_0[i] * tbe_0 + 4.0 * tr_xxyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_0[i] = -6.0 * tr_xyy_z[i] * tke_0 - 6.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyy_yz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_0[i] = 2.0 * tr_xy_0[i] - 4.0 * tr_xyz_z[i] * tke_0 - 4.0 * tr_xyzz_0[i] * tbe_0 - 2.0 * tr_xyy_y[i] * tke_0 + 4.0 * tr_xyyz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_0[i] * tbe_0 + 4.0 * tr_xyyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_0[i] = 2.0 * tr_xz_0[i] - 2.0 * tr_xzz_z[i] * tke_0 - 2.0 * tr_xzzz_0[i] * tbe_0 - 4.0 * tr_xyz_y[i] * tke_0 + 4.0 * tr_xyzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_y[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_0[i] * tbe_0 + 4.0 * tr_xyyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_0[i] = -6.0 * tr_xzz_y[i] * tke_0 + 4.0 * tr_xzzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_y[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_0[i] * tbe_0 + 4.0 * tr_xyzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_0[i] = -8.0 * tr_yyy_z[i] * tke_0 - 8.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyy_yz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_y[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_0[i] = 3.0 * tr_yy_0[i] - 6.0 * tr_yyz_z[i] * tke_0 - 6.0 * tr_yyzz_0[i] * tbe_0 - 2.0 * tr_yyy_y[i] * tke_0 + 4.0 * tr_yyyz_yz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_y[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_0[i] * tbe_0 + 4.0 * tr_yyyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_0[i] = 4.0 * tr_yz_0[i] - 4.0 * tr_yzz_z[i] * tke_0 - 4.0 * tr_yzzz_0[i] * tbe_0 - 4.0 * tr_yyz_y[i] * tke_0 + 4.0 * tr_yyzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_y[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_0[i] * tbe_0 + 4.0 * tr_yyyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_0[i] = 3.0 * tr_zz_0[i] - 2.0 * tr_zzz_z[i] * tke_0 - 2.0 * tr_zzzz_0[i] * tbe_0 - 6.0 * tr_yzz_y[i] * tke_0 + 4.0 * tr_yzzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_y[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_0[i] * tbe_0 + 4.0 * tr_yyzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_0[i] = -8.0 * tr_zzz_y[i] * tke_0 + 4.0 * tr_zzzz_yz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_y[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_0[i] * tbe_0 + 4.0 * tr_yzzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_0[i] = -2.0 * tr_xxxx_0[i] * tbe_0 - 2.0 * tr_xxxx_0[i] * tke_0 + 4.0 * tr_xxxx_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_0[i] = -2.0 * tr_xxxy_0[i] * tbe_0 - 2.0 * tr_xxxy_0[i] * tke_0 + 4.0 * tr_xxxy_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_0[i] = -4.0 * tr_xxx_z[i] * tke_0 - 6.0 * tr_xxxz_0[i] * tbe_0 - 2.0 * tr_xxxz_0[i] * tke_0 + 4.0 * tr_xxxz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_0[i] = -2.0 * tr_xxyy_0[i] * tbe_0 - 2.0 * tr_xxyy_0[i] * tke_0 + 4.0 * tr_xxyy_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_0[i] = -4.0 * tr_xxy_z[i] * tke_0 - 6.0 * tr_xxyz_0[i] * tbe_0 - 2.0 * tr_xxyz_0[i] * tke_0 + 4.0 * tr_xxyz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_0[i] = 2.0 * tr_xx_0[i] - 8.0 * tr_xxz_z[i] * tke_0 - 10.0 * tr_xxzz_0[i] * tbe_0 - 2.0 * tr_xxzz_0[i] * tke_0 + 4.0 * tr_xxzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_0[i] = -2.0 * tr_xyyy_0[i] * tbe_0 - 2.0 * tr_xyyy_0[i] * tke_0 + 4.0 * tr_xyyy_zz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_0[i] = -4.0 * tr_xyy_z[i] * tke_0 - 6.0 * tr_xyyz_0[i] * tbe_0 - 2.0 * tr_xyyz_0[i] * tke_0 + 4.0 * tr_xyyz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_0[i] = 2.0 * tr_xy_0[i] - 8.0 * tr_xyz_z[i] * tke_0 - 10.0 * tr_xyzz_0[i] * tbe_0 - 2.0 * tr_xyzz_0[i] * tke_0 + 4.0 * tr_xyzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_0[i] = 6.0 * tr_xz_0[i] - 12.0 * tr_xzz_z[i] * tke_0 - 14.0 * tr_xzzz_0[i] * tbe_0 - 2.0 * tr_xzzz_0[i] * tke_0 + 4.0 * tr_xzzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_0[i] = -2.0 * tr_yyyy_0[i] * tbe_0 - 2.0 * tr_yyyy_0[i] * tke_0 + 4.0 * tr_yyyy_zz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_0[i] = -4.0 * tr_yyy_z[i] * tke_0 - 6.0 * tr_yyyz_0[i] * tbe_0 - 2.0 * tr_yyyz_0[i] * tke_0 + 4.0 * tr_yyyz_zz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_0[i] = 2.0 * tr_yy_0[i] - 8.0 * tr_yyz_z[i] * tke_0 - 10.0 * tr_yyzz_0[i] * tbe_0 - 2.0 * tr_yyzz_0[i] * tke_0 + 4.0 * tr_yyzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_0[i] = 6.0 * tr_yz_0[i] - 12.0 * tr_yzz_z[i] * tke_0 - 14.0 * tr_yzzz_0[i] * tbe_0 - 2.0 * tr_yzzz_0[i] * tke_0 + 4.0 * tr_yzzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_0[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_0[i] = 12.0 * tr_zz_0[i] - 16.0 * tr_zzz_z[i] * tke_0 - 18.0 * tr_zzzz_0[i] * tbe_0 - 2.0 * tr_zzzz_0[i] * tke_0 + 4.0 * tr_zzzz_zz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_z[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_0[i] * tbe_0 * tbe_0;
    }
}

} // t2cgeom namespace

