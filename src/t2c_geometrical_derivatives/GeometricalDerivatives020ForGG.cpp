#include "GeometricalDerivatives020ForGG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_gg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_gg,
                         const int idx_op_dg,
                         const int idx_op_ff,
                         const int idx_op_fh,
                         const int idx_op_gd,
                         const int idx_op_gg,
                         const int idx_op_gi,
                         const int idx_op_hf,
                         const int idx_op_hh,
                         const int idx_op_ig,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up components of auxiliary buffer : GI

    auto tr_xxxx_xxxxxx = pbuffer.data(idx_op_gi);

    auto tr_xxxx_xxxxxy = pbuffer.data(idx_op_gi + 1);

    auto tr_xxxx_xxxxxz = pbuffer.data(idx_op_gi + 2);

    auto tr_xxxx_xxxxyy = pbuffer.data(idx_op_gi + 3);

    auto tr_xxxx_xxxxyz = pbuffer.data(idx_op_gi + 4);

    auto tr_xxxx_xxxxzz = pbuffer.data(idx_op_gi + 5);

    auto tr_xxxx_xxxyyy = pbuffer.data(idx_op_gi + 6);

    auto tr_xxxx_xxxyyz = pbuffer.data(idx_op_gi + 7);

    auto tr_xxxx_xxxyzz = pbuffer.data(idx_op_gi + 8);

    auto tr_xxxx_xxxzzz = pbuffer.data(idx_op_gi + 9);

    auto tr_xxxx_xxyyyy = pbuffer.data(idx_op_gi + 10);

    auto tr_xxxx_xxyyyz = pbuffer.data(idx_op_gi + 11);

    auto tr_xxxx_xxyyzz = pbuffer.data(idx_op_gi + 12);

    auto tr_xxxx_xxyzzz = pbuffer.data(idx_op_gi + 13);

    auto tr_xxxx_xxzzzz = pbuffer.data(idx_op_gi + 14);

    auto tr_xxxx_xyyyyy = pbuffer.data(idx_op_gi + 15);

    auto tr_xxxx_xyyyyz = pbuffer.data(idx_op_gi + 16);

    auto tr_xxxx_xyyyzz = pbuffer.data(idx_op_gi + 17);

    auto tr_xxxx_xyyzzz = pbuffer.data(idx_op_gi + 18);

    auto tr_xxxx_xyzzzz = pbuffer.data(idx_op_gi + 19);

    auto tr_xxxx_xzzzzz = pbuffer.data(idx_op_gi + 20);

    auto tr_xxxx_yyyyyy = pbuffer.data(idx_op_gi + 21);

    auto tr_xxxx_yyyyyz = pbuffer.data(idx_op_gi + 22);

    auto tr_xxxx_yyyyzz = pbuffer.data(idx_op_gi + 23);

    auto tr_xxxx_yyyzzz = pbuffer.data(idx_op_gi + 24);

    auto tr_xxxx_yyzzzz = pbuffer.data(idx_op_gi + 25);

    auto tr_xxxx_yzzzzz = pbuffer.data(idx_op_gi + 26);

    auto tr_xxxx_zzzzzz = pbuffer.data(idx_op_gi + 27);

    auto tr_xxxy_xxxxxx = pbuffer.data(idx_op_gi + 28);

    auto tr_xxxy_xxxxxy = pbuffer.data(idx_op_gi + 29);

    auto tr_xxxy_xxxxxz = pbuffer.data(idx_op_gi + 30);

    auto tr_xxxy_xxxxyy = pbuffer.data(idx_op_gi + 31);

    auto tr_xxxy_xxxxyz = pbuffer.data(idx_op_gi + 32);

    auto tr_xxxy_xxxxzz = pbuffer.data(idx_op_gi + 33);

    auto tr_xxxy_xxxyyy = pbuffer.data(idx_op_gi + 34);

    auto tr_xxxy_xxxyyz = pbuffer.data(idx_op_gi + 35);

    auto tr_xxxy_xxxyzz = pbuffer.data(idx_op_gi + 36);

    auto tr_xxxy_xxxzzz = pbuffer.data(idx_op_gi + 37);

    auto tr_xxxy_xxyyyy = pbuffer.data(idx_op_gi + 38);

    auto tr_xxxy_xxyyyz = pbuffer.data(idx_op_gi + 39);

    auto tr_xxxy_xxyyzz = pbuffer.data(idx_op_gi + 40);

    auto tr_xxxy_xxyzzz = pbuffer.data(idx_op_gi + 41);

    auto tr_xxxy_xxzzzz = pbuffer.data(idx_op_gi + 42);

    auto tr_xxxy_xyyyyy = pbuffer.data(idx_op_gi + 43);

    auto tr_xxxy_xyyyyz = pbuffer.data(idx_op_gi + 44);

    auto tr_xxxy_xyyyzz = pbuffer.data(idx_op_gi + 45);

    auto tr_xxxy_xyyzzz = pbuffer.data(idx_op_gi + 46);

    auto tr_xxxy_xyzzzz = pbuffer.data(idx_op_gi + 47);

    auto tr_xxxy_xzzzzz = pbuffer.data(idx_op_gi + 48);

    auto tr_xxxy_yyyyyy = pbuffer.data(idx_op_gi + 49);

    auto tr_xxxy_yyyyyz = pbuffer.data(idx_op_gi + 50);

    auto tr_xxxy_yyyyzz = pbuffer.data(idx_op_gi + 51);

    auto tr_xxxy_yyyzzz = pbuffer.data(idx_op_gi + 52);

    auto tr_xxxy_yyzzzz = pbuffer.data(idx_op_gi + 53);

    auto tr_xxxy_yzzzzz = pbuffer.data(idx_op_gi + 54);

    auto tr_xxxy_zzzzzz = pbuffer.data(idx_op_gi + 55);

    auto tr_xxxz_xxxxxx = pbuffer.data(idx_op_gi + 56);

    auto tr_xxxz_xxxxxy = pbuffer.data(idx_op_gi + 57);

    auto tr_xxxz_xxxxxz = pbuffer.data(idx_op_gi + 58);

    auto tr_xxxz_xxxxyy = pbuffer.data(idx_op_gi + 59);

    auto tr_xxxz_xxxxyz = pbuffer.data(idx_op_gi + 60);

    auto tr_xxxz_xxxxzz = pbuffer.data(idx_op_gi + 61);

    auto tr_xxxz_xxxyyy = pbuffer.data(idx_op_gi + 62);

    auto tr_xxxz_xxxyyz = pbuffer.data(idx_op_gi + 63);

    auto tr_xxxz_xxxyzz = pbuffer.data(idx_op_gi + 64);

    auto tr_xxxz_xxxzzz = pbuffer.data(idx_op_gi + 65);

    auto tr_xxxz_xxyyyy = pbuffer.data(idx_op_gi + 66);

    auto tr_xxxz_xxyyyz = pbuffer.data(idx_op_gi + 67);

    auto tr_xxxz_xxyyzz = pbuffer.data(idx_op_gi + 68);

    auto tr_xxxz_xxyzzz = pbuffer.data(idx_op_gi + 69);

    auto tr_xxxz_xxzzzz = pbuffer.data(idx_op_gi + 70);

    auto tr_xxxz_xyyyyy = pbuffer.data(idx_op_gi + 71);

    auto tr_xxxz_xyyyyz = pbuffer.data(idx_op_gi + 72);

    auto tr_xxxz_xyyyzz = pbuffer.data(idx_op_gi + 73);

    auto tr_xxxz_xyyzzz = pbuffer.data(idx_op_gi + 74);

    auto tr_xxxz_xyzzzz = pbuffer.data(idx_op_gi + 75);

    auto tr_xxxz_xzzzzz = pbuffer.data(idx_op_gi + 76);

    auto tr_xxxz_yyyyyy = pbuffer.data(idx_op_gi + 77);

    auto tr_xxxz_yyyyyz = pbuffer.data(idx_op_gi + 78);

    auto tr_xxxz_yyyyzz = pbuffer.data(idx_op_gi + 79);

    auto tr_xxxz_yyyzzz = pbuffer.data(idx_op_gi + 80);

    auto tr_xxxz_yyzzzz = pbuffer.data(idx_op_gi + 81);

    auto tr_xxxz_yzzzzz = pbuffer.data(idx_op_gi + 82);

    auto tr_xxxz_zzzzzz = pbuffer.data(idx_op_gi + 83);

    auto tr_xxyy_xxxxxx = pbuffer.data(idx_op_gi + 84);

    auto tr_xxyy_xxxxxy = pbuffer.data(idx_op_gi + 85);

    auto tr_xxyy_xxxxxz = pbuffer.data(idx_op_gi + 86);

    auto tr_xxyy_xxxxyy = pbuffer.data(idx_op_gi + 87);

    auto tr_xxyy_xxxxyz = pbuffer.data(idx_op_gi + 88);

    auto tr_xxyy_xxxxzz = pbuffer.data(idx_op_gi + 89);

    auto tr_xxyy_xxxyyy = pbuffer.data(idx_op_gi + 90);

    auto tr_xxyy_xxxyyz = pbuffer.data(idx_op_gi + 91);

    auto tr_xxyy_xxxyzz = pbuffer.data(idx_op_gi + 92);

    auto tr_xxyy_xxxzzz = pbuffer.data(idx_op_gi + 93);

    auto tr_xxyy_xxyyyy = pbuffer.data(idx_op_gi + 94);

    auto tr_xxyy_xxyyyz = pbuffer.data(idx_op_gi + 95);

    auto tr_xxyy_xxyyzz = pbuffer.data(idx_op_gi + 96);

    auto tr_xxyy_xxyzzz = pbuffer.data(idx_op_gi + 97);

    auto tr_xxyy_xxzzzz = pbuffer.data(idx_op_gi + 98);

    auto tr_xxyy_xyyyyy = pbuffer.data(idx_op_gi + 99);

    auto tr_xxyy_xyyyyz = pbuffer.data(idx_op_gi + 100);

    auto tr_xxyy_xyyyzz = pbuffer.data(idx_op_gi + 101);

    auto tr_xxyy_xyyzzz = pbuffer.data(idx_op_gi + 102);

    auto tr_xxyy_xyzzzz = pbuffer.data(idx_op_gi + 103);

    auto tr_xxyy_xzzzzz = pbuffer.data(idx_op_gi + 104);

    auto tr_xxyy_yyyyyy = pbuffer.data(idx_op_gi + 105);

    auto tr_xxyy_yyyyyz = pbuffer.data(idx_op_gi + 106);

    auto tr_xxyy_yyyyzz = pbuffer.data(idx_op_gi + 107);

    auto tr_xxyy_yyyzzz = pbuffer.data(idx_op_gi + 108);

    auto tr_xxyy_yyzzzz = pbuffer.data(idx_op_gi + 109);

    auto tr_xxyy_yzzzzz = pbuffer.data(idx_op_gi + 110);

    auto tr_xxyy_zzzzzz = pbuffer.data(idx_op_gi + 111);

    auto tr_xxyz_xxxxxx = pbuffer.data(idx_op_gi + 112);

    auto tr_xxyz_xxxxxy = pbuffer.data(idx_op_gi + 113);

    auto tr_xxyz_xxxxxz = pbuffer.data(idx_op_gi + 114);

    auto tr_xxyz_xxxxyy = pbuffer.data(idx_op_gi + 115);

    auto tr_xxyz_xxxxyz = pbuffer.data(idx_op_gi + 116);

    auto tr_xxyz_xxxxzz = pbuffer.data(idx_op_gi + 117);

    auto tr_xxyz_xxxyyy = pbuffer.data(idx_op_gi + 118);

    auto tr_xxyz_xxxyyz = pbuffer.data(idx_op_gi + 119);

    auto tr_xxyz_xxxyzz = pbuffer.data(idx_op_gi + 120);

    auto tr_xxyz_xxxzzz = pbuffer.data(idx_op_gi + 121);

    auto tr_xxyz_xxyyyy = pbuffer.data(idx_op_gi + 122);

    auto tr_xxyz_xxyyyz = pbuffer.data(idx_op_gi + 123);

    auto tr_xxyz_xxyyzz = pbuffer.data(idx_op_gi + 124);

    auto tr_xxyz_xxyzzz = pbuffer.data(idx_op_gi + 125);

    auto tr_xxyz_xxzzzz = pbuffer.data(idx_op_gi + 126);

    auto tr_xxyz_xyyyyy = pbuffer.data(idx_op_gi + 127);

    auto tr_xxyz_xyyyyz = pbuffer.data(idx_op_gi + 128);

    auto tr_xxyz_xyyyzz = pbuffer.data(idx_op_gi + 129);

    auto tr_xxyz_xyyzzz = pbuffer.data(idx_op_gi + 130);

    auto tr_xxyz_xyzzzz = pbuffer.data(idx_op_gi + 131);

    auto tr_xxyz_xzzzzz = pbuffer.data(idx_op_gi + 132);

    auto tr_xxyz_yyyyyy = pbuffer.data(idx_op_gi + 133);

    auto tr_xxyz_yyyyyz = pbuffer.data(idx_op_gi + 134);

    auto tr_xxyz_yyyyzz = pbuffer.data(idx_op_gi + 135);

    auto tr_xxyz_yyyzzz = pbuffer.data(idx_op_gi + 136);

    auto tr_xxyz_yyzzzz = pbuffer.data(idx_op_gi + 137);

    auto tr_xxyz_yzzzzz = pbuffer.data(idx_op_gi + 138);

    auto tr_xxyz_zzzzzz = pbuffer.data(idx_op_gi + 139);

    auto tr_xxzz_xxxxxx = pbuffer.data(idx_op_gi + 140);

    auto tr_xxzz_xxxxxy = pbuffer.data(idx_op_gi + 141);

    auto tr_xxzz_xxxxxz = pbuffer.data(idx_op_gi + 142);

    auto tr_xxzz_xxxxyy = pbuffer.data(idx_op_gi + 143);

    auto tr_xxzz_xxxxyz = pbuffer.data(idx_op_gi + 144);

    auto tr_xxzz_xxxxzz = pbuffer.data(idx_op_gi + 145);

    auto tr_xxzz_xxxyyy = pbuffer.data(idx_op_gi + 146);

    auto tr_xxzz_xxxyyz = pbuffer.data(idx_op_gi + 147);

    auto tr_xxzz_xxxyzz = pbuffer.data(idx_op_gi + 148);

    auto tr_xxzz_xxxzzz = pbuffer.data(idx_op_gi + 149);

    auto tr_xxzz_xxyyyy = pbuffer.data(idx_op_gi + 150);

    auto tr_xxzz_xxyyyz = pbuffer.data(idx_op_gi + 151);

    auto tr_xxzz_xxyyzz = pbuffer.data(idx_op_gi + 152);

    auto tr_xxzz_xxyzzz = pbuffer.data(idx_op_gi + 153);

    auto tr_xxzz_xxzzzz = pbuffer.data(idx_op_gi + 154);

    auto tr_xxzz_xyyyyy = pbuffer.data(idx_op_gi + 155);

    auto tr_xxzz_xyyyyz = pbuffer.data(idx_op_gi + 156);

    auto tr_xxzz_xyyyzz = pbuffer.data(idx_op_gi + 157);

    auto tr_xxzz_xyyzzz = pbuffer.data(idx_op_gi + 158);

    auto tr_xxzz_xyzzzz = pbuffer.data(idx_op_gi + 159);

    auto tr_xxzz_xzzzzz = pbuffer.data(idx_op_gi + 160);

    auto tr_xxzz_yyyyyy = pbuffer.data(idx_op_gi + 161);

    auto tr_xxzz_yyyyyz = pbuffer.data(idx_op_gi + 162);

    auto tr_xxzz_yyyyzz = pbuffer.data(idx_op_gi + 163);

    auto tr_xxzz_yyyzzz = pbuffer.data(idx_op_gi + 164);

    auto tr_xxzz_yyzzzz = pbuffer.data(idx_op_gi + 165);

    auto tr_xxzz_yzzzzz = pbuffer.data(idx_op_gi + 166);

    auto tr_xxzz_zzzzzz = pbuffer.data(idx_op_gi + 167);

    auto tr_xyyy_xxxxxx = pbuffer.data(idx_op_gi + 168);

    auto tr_xyyy_xxxxxy = pbuffer.data(idx_op_gi + 169);

    auto tr_xyyy_xxxxxz = pbuffer.data(idx_op_gi + 170);

    auto tr_xyyy_xxxxyy = pbuffer.data(idx_op_gi + 171);

    auto tr_xyyy_xxxxyz = pbuffer.data(idx_op_gi + 172);

    auto tr_xyyy_xxxxzz = pbuffer.data(idx_op_gi + 173);

    auto tr_xyyy_xxxyyy = pbuffer.data(idx_op_gi + 174);

    auto tr_xyyy_xxxyyz = pbuffer.data(idx_op_gi + 175);

    auto tr_xyyy_xxxyzz = pbuffer.data(idx_op_gi + 176);

    auto tr_xyyy_xxxzzz = pbuffer.data(idx_op_gi + 177);

    auto tr_xyyy_xxyyyy = pbuffer.data(idx_op_gi + 178);

    auto tr_xyyy_xxyyyz = pbuffer.data(idx_op_gi + 179);

    auto tr_xyyy_xxyyzz = pbuffer.data(idx_op_gi + 180);

    auto tr_xyyy_xxyzzz = pbuffer.data(idx_op_gi + 181);

    auto tr_xyyy_xxzzzz = pbuffer.data(idx_op_gi + 182);

    auto tr_xyyy_xyyyyy = pbuffer.data(idx_op_gi + 183);

    auto tr_xyyy_xyyyyz = pbuffer.data(idx_op_gi + 184);

    auto tr_xyyy_xyyyzz = pbuffer.data(idx_op_gi + 185);

    auto tr_xyyy_xyyzzz = pbuffer.data(idx_op_gi + 186);

    auto tr_xyyy_xyzzzz = pbuffer.data(idx_op_gi + 187);

    auto tr_xyyy_xzzzzz = pbuffer.data(idx_op_gi + 188);

    auto tr_xyyy_yyyyyy = pbuffer.data(idx_op_gi + 189);

    auto tr_xyyy_yyyyyz = pbuffer.data(idx_op_gi + 190);

    auto tr_xyyy_yyyyzz = pbuffer.data(idx_op_gi + 191);

    auto tr_xyyy_yyyzzz = pbuffer.data(idx_op_gi + 192);

    auto tr_xyyy_yyzzzz = pbuffer.data(idx_op_gi + 193);

    auto tr_xyyy_yzzzzz = pbuffer.data(idx_op_gi + 194);

    auto tr_xyyy_zzzzzz = pbuffer.data(idx_op_gi + 195);

    auto tr_xyyz_xxxxxx = pbuffer.data(idx_op_gi + 196);

    auto tr_xyyz_xxxxxy = pbuffer.data(idx_op_gi + 197);

    auto tr_xyyz_xxxxxz = pbuffer.data(idx_op_gi + 198);

    auto tr_xyyz_xxxxyy = pbuffer.data(idx_op_gi + 199);

    auto tr_xyyz_xxxxyz = pbuffer.data(idx_op_gi + 200);

    auto tr_xyyz_xxxxzz = pbuffer.data(idx_op_gi + 201);

    auto tr_xyyz_xxxyyy = pbuffer.data(idx_op_gi + 202);

    auto tr_xyyz_xxxyyz = pbuffer.data(idx_op_gi + 203);

    auto tr_xyyz_xxxyzz = pbuffer.data(idx_op_gi + 204);

    auto tr_xyyz_xxxzzz = pbuffer.data(idx_op_gi + 205);

    auto tr_xyyz_xxyyyy = pbuffer.data(idx_op_gi + 206);

    auto tr_xyyz_xxyyyz = pbuffer.data(idx_op_gi + 207);

    auto tr_xyyz_xxyyzz = pbuffer.data(idx_op_gi + 208);

    auto tr_xyyz_xxyzzz = pbuffer.data(idx_op_gi + 209);

    auto tr_xyyz_xxzzzz = pbuffer.data(idx_op_gi + 210);

    auto tr_xyyz_xyyyyy = pbuffer.data(idx_op_gi + 211);

    auto tr_xyyz_xyyyyz = pbuffer.data(idx_op_gi + 212);

    auto tr_xyyz_xyyyzz = pbuffer.data(idx_op_gi + 213);

    auto tr_xyyz_xyyzzz = pbuffer.data(idx_op_gi + 214);

    auto tr_xyyz_xyzzzz = pbuffer.data(idx_op_gi + 215);

    auto tr_xyyz_xzzzzz = pbuffer.data(idx_op_gi + 216);

    auto tr_xyyz_yyyyyy = pbuffer.data(idx_op_gi + 217);

    auto tr_xyyz_yyyyyz = pbuffer.data(idx_op_gi + 218);

    auto tr_xyyz_yyyyzz = pbuffer.data(idx_op_gi + 219);

    auto tr_xyyz_yyyzzz = pbuffer.data(idx_op_gi + 220);

    auto tr_xyyz_yyzzzz = pbuffer.data(idx_op_gi + 221);

    auto tr_xyyz_yzzzzz = pbuffer.data(idx_op_gi + 222);

    auto tr_xyyz_zzzzzz = pbuffer.data(idx_op_gi + 223);

    auto tr_xyzz_xxxxxx = pbuffer.data(idx_op_gi + 224);

    auto tr_xyzz_xxxxxy = pbuffer.data(idx_op_gi + 225);

    auto tr_xyzz_xxxxxz = pbuffer.data(idx_op_gi + 226);

    auto tr_xyzz_xxxxyy = pbuffer.data(idx_op_gi + 227);

    auto tr_xyzz_xxxxyz = pbuffer.data(idx_op_gi + 228);

    auto tr_xyzz_xxxxzz = pbuffer.data(idx_op_gi + 229);

    auto tr_xyzz_xxxyyy = pbuffer.data(idx_op_gi + 230);

    auto tr_xyzz_xxxyyz = pbuffer.data(idx_op_gi + 231);

    auto tr_xyzz_xxxyzz = pbuffer.data(idx_op_gi + 232);

    auto tr_xyzz_xxxzzz = pbuffer.data(idx_op_gi + 233);

    auto tr_xyzz_xxyyyy = pbuffer.data(idx_op_gi + 234);

    auto tr_xyzz_xxyyyz = pbuffer.data(idx_op_gi + 235);

    auto tr_xyzz_xxyyzz = pbuffer.data(idx_op_gi + 236);

    auto tr_xyzz_xxyzzz = pbuffer.data(idx_op_gi + 237);

    auto tr_xyzz_xxzzzz = pbuffer.data(idx_op_gi + 238);

    auto tr_xyzz_xyyyyy = pbuffer.data(idx_op_gi + 239);

    auto tr_xyzz_xyyyyz = pbuffer.data(idx_op_gi + 240);

    auto tr_xyzz_xyyyzz = pbuffer.data(idx_op_gi + 241);

    auto tr_xyzz_xyyzzz = pbuffer.data(idx_op_gi + 242);

    auto tr_xyzz_xyzzzz = pbuffer.data(idx_op_gi + 243);

    auto tr_xyzz_xzzzzz = pbuffer.data(idx_op_gi + 244);

    auto tr_xyzz_yyyyyy = pbuffer.data(idx_op_gi + 245);

    auto tr_xyzz_yyyyyz = pbuffer.data(idx_op_gi + 246);

    auto tr_xyzz_yyyyzz = pbuffer.data(idx_op_gi + 247);

    auto tr_xyzz_yyyzzz = pbuffer.data(idx_op_gi + 248);

    auto tr_xyzz_yyzzzz = pbuffer.data(idx_op_gi + 249);

    auto tr_xyzz_yzzzzz = pbuffer.data(idx_op_gi + 250);

    auto tr_xyzz_zzzzzz = pbuffer.data(idx_op_gi + 251);

    auto tr_xzzz_xxxxxx = pbuffer.data(idx_op_gi + 252);

    auto tr_xzzz_xxxxxy = pbuffer.data(idx_op_gi + 253);

    auto tr_xzzz_xxxxxz = pbuffer.data(idx_op_gi + 254);

    auto tr_xzzz_xxxxyy = pbuffer.data(idx_op_gi + 255);

    auto tr_xzzz_xxxxyz = pbuffer.data(idx_op_gi + 256);

    auto tr_xzzz_xxxxzz = pbuffer.data(idx_op_gi + 257);

    auto tr_xzzz_xxxyyy = pbuffer.data(idx_op_gi + 258);

    auto tr_xzzz_xxxyyz = pbuffer.data(idx_op_gi + 259);

    auto tr_xzzz_xxxyzz = pbuffer.data(idx_op_gi + 260);

    auto tr_xzzz_xxxzzz = pbuffer.data(idx_op_gi + 261);

    auto tr_xzzz_xxyyyy = pbuffer.data(idx_op_gi + 262);

    auto tr_xzzz_xxyyyz = pbuffer.data(idx_op_gi + 263);

    auto tr_xzzz_xxyyzz = pbuffer.data(idx_op_gi + 264);

    auto tr_xzzz_xxyzzz = pbuffer.data(idx_op_gi + 265);

    auto tr_xzzz_xxzzzz = pbuffer.data(idx_op_gi + 266);

    auto tr_xzzz_xyyyyy = pbuffer.data(idx_op_gi + 267);

    auto tr_xzzz_xyyyyz = pbuffer.data(idx_op_gi + 268);

    auto tr_xzzz_xyyyzz = pbuffer.data(idx_op_gi + 269);

    auto tr_xzzz_xyyzzz = pbuffer.data(idx_op_gi + 270);

    auto tr_xzzz_xyzzzz = pbuffer.data(idx_op_gi + 271);

    auto tr_xzzz_xzzzzz = pbuffer.data(idx_op_gi + 272);

    auto tr_xzzz_yyyyyy = pbuffer.data(idx_op_gi + 273);

    auto tr_xzzz_yyyyyz = pbuffer.data(idx_op_gi + 274);

    auto tr_xzzz_yyyyzz = pbuffer.data(idx_op_gi + 275);

    auto tr_xzzz_yyyzzz = pbuffer.data(idx_op_gi + 276);

    auto tr_xzzz_yyzzzz = pbuffer.data(idx_op_gi + 277);

    auto tr_xzzz_yzzzzz = pbuffer.data(idx_op_gi + 278);

    auto tr_xzzz_zzzzzz = pbuffer.data(idx_op_gi + 279);

    auto tr_yyyy_xxxxxx = pbuffer.data(idx_op_gi + 280);

    auto tr_yyyy_xxxxxy = pbuffer.data(idx_op_gi + 281);

    auto tr_yyyy_xxxxxz = pbuffer.data(idx_op_gi + 282);

    auto tr_yyyy_xxxxyy = pbuffer.data(idx_op_gi + 283);

    auto tr_yyyy_xxxxyz = pbuffer.data(idx_op_gi + 284);

    auto tr_yyyy_xxxxzz = pbuffer.data(idx_op_gi + 285);

    auto tr_yyyy_xxxyyy = pbuffer.data(idx_op_gi + 286);

    auto tr_yyyy_xxxyyz = pbuffer.data(idx_op_gi + 287);

    auto tr_yyyy_xxxyzz = pbuffer.data(idx_op_gi + 288);

    auto tr_yyyy_xxxzzz = pbuffer.data(idx_op_gi + 289);

    auto tr_yyyy_xxyyyy = pbuffer.data(idx_op_gi + 290);

    auto tr_yyyy_xxyyyz = pbuffer.data(idx_op_gi + 291);

    auto tr_yyyy_xxyyzz = pbuffer.data(idx_op_gi + 292);

    auto tr_yyyy_xxyzzz = pbuffer.data(idx_op_gi + 293);

    auto tr_yyyy_xxzzzz = pbuffer.data(idx_op_gi + 294);

    auto tr_yyyy_xyyyyy = pbuffer.data(idx_op_gi + 295);

    auto tr_yyyy_xyyyyz = pbuffer.data(idx_op_gi + 296);

    auto tr_yyyy_xyyyzz = pbuffer.data(idx_op_gi + 297);

    auto tr_yyyy_xyyzzz = pbuffer.data(idx_op_gi + 298);

    auto tr_yyyy_xyzzzz = pbuffer.data(idx_op_gi + 299);

    auto tr_yyyy_xzzzzz = pbuffer.data(idx_op_gi + 300);

    auto tr_yyyy_yyyyyy = pbuffer.data(idx_op_gi + 301);

    auto tr_yyyy_yyyyyz = pbuffer.data(idx_op_gi + 302);

    auto tr_yyyy_yyyyzz = pbuffer.data(idx_op_gi + 303);

    auto tr_yyyy_yyyzzz = pbuffer.data(idx_op_gi + 304);

    auto tr_yyyy_yyzzzz = pbuffer.data(idx_op_gi + 305);

    auto tr_yyyy_yzzzzz = pbuffer.data(idx_op_gi + 306);

    auto tr_yyyy_zzzzzz = pbuffer.data(idx_op_gi + 307);

    auto tr_yyyz_xxxxxx = pbuffer.data(idx_op_gi + 308);

    auto tr_yyyz_xxxxxy = pbuffer.data(idx_op_gi + 309);

    auto tr_yyyz_xxxxxz = pbuffer.data(idx_op_gi + 310);

    auto tr_yyyz_xxxxyy = pbuffer.data(idx_op_gi + 311);

    auto tr_yyyz_xxxxyz = pbuffer.data(idx_op_gi + 312);

    auto tr_yyyz_xxxxzz = pbuffer.data(idx_op_gi + 313);

    auto tr_yyyz_xxxyyy = pbuffer.data(idx_op_gi + 314);

    auto tr_yyyz_xxxyyz = pbuffer.data(idx_op_gi + 315);

    auto tr_yyyz_xxxyzz = pbuffer.data(idx_op_gi + 316);

    auto tr_yyyz_xxxzzz = pbuffer.data(idx_op_gi + 317);

    auto tr_yyyz_xxyyyy = pbuffer.data(idx_op_gi + 318);

    auto tr_yyyz_xxyyyz = pbuffer.data(idx_op_gi + 319);

    auto tr_yyyz_xxyyzz = pbuffer.data(idx_op_gi + 320);

    auto tr_yyyz_xxyzzz = pbuffer.data(idx_op_gi + 321);

    auto tr_yyyz_xxzzzz = pbuffer.data(idx_op_gi + 322);

    auto tr_yyyz_xyyyyy = pbuffer.data(idx_op_gi + 323);

    auto tr_yyyz_xyyyyz = pbuffer.data(idx_op_gi + 324);

    auto tr_yyyz_xyyyzz = pbuffer.data(idx_op_gi + 325);

    auto tr_yyyz_xyyzzz = pbuffer.data(idx_op_gi + 326);

    auto tr_yyyz_xyzzzz = pbuffer.data(idx_op_gi + 327);

    auto tr_yyyz_xzzzzz = pbuffer.data(idx_op_gi + 328);

    auto tr_yyyz_yyyyyy = pbuffer.data(idx_op_gi + 329);

    auto tr_yyyz_yyyyyz = pbuffer.data(idx_op_gi + 330);

    auto tr_yyyz_yyyyzz = pbuffer.data(idx_op_gi + 331);

    auto tr_yyyz_yyyzzz = pbuffer.data(idx_op_gi + 332);

    auto tr_yyyz_yyzzzz = pbuffer.data(idx_op_gi + 333);

    auto tr_yyyz_yzzzzz = pbuffer.data(idx_op_gi + 334);

    auto tr_yyyz_zzzzzz = pbuffer.data(idx_op_gi + 335);

    auto tr_yyzz_xxxxxx = pbuffer.data(idx_op_gi + 336);

    auto tr_yyzz_xxxxxy = pbuffer.data(idx_op_gi + 337);

    auto tr_yyzz_xxxxxz = pbuffer.data(idx_op_gi + 338);

    auto tr_yyzz_xxxxyy = pbuffer.data(idx_op_gi + 339);

    auto tr_yyzz_xxxxyz = pbuffer.data(idx_op_gi + 340);

    auto tr_yyzz_xxxxzz = pbuffer.data(idx_op_gi + 341);

    auto tr_yyzz_xxxyyy = pbuffer.data(idx_op_gi + 342);

    auto tr_yyzz_xxxyyz = pbuffer.data(idx_op_gi + 343);

    auto tr_yyzz_xxxyzz = pbuffer.data(idx_op_gi + 344);

    auto tr_yyzz_xxxzzz = pbuffer.data(idx_op_gi + 345);

    auto tr_yyzz_xxyyyy = pbuffer.data(idx_op_gi + 346);

    auto tr_yyzz_xxyyyz = pbuffer.data(idx_op_gi + 347);

    auto tr_yyzz_xxyyzz = pbuffer.data(idx_op_gi + 348);

    auto tr_yyzz_xxyzzz = pbuffer.data(idx_op_gi + 349);

    auto tr_yyzz_xxzzzz = pbuffer.data(idx_op_gi + 350);

    auto tr_yyzz_xyyyyy = pbuffer.data(idx_op_gi + 351);

    auto tr_yyzz_xyyyyz = pbuffer.data(idx_op_gi + 352);

    auto tr_yyzz_xyyyzz = pbuffer.data(idx_op_gi + 353);

    auto tr_yyzz_xyyzzz = pbuffer.data(idx_op_gi + 354);

    auto tr_yyzz_xyzzzz = pbuffer.data(idx_op_gi + 355);

    auto tr_yyzz_xzzzzz = pbuffer.data(idx_op_gi + 356);

    auto tr_yyzz_yyyyyy = pbuffer.data(idx_op_gi + 357);

    auto tr_yyzz_yyyyyz = pbuffer.data(idx_op_gi + 358);

    auto tr_yyzz_yyyyzz = pbuffer.data(idx_op_gi + 359);

    auto tr_yyzz_yyyzzz = pbuffer.data(idx_op_gi + 360);

    auto tr_yyzz_yyzzzz = pbuffer.data(idx_op_gi + 361);

    auto tr_yyzz_yzzzzz = pbuffer.data(idx_op_gi + 362);

    auto tr_yyzz_zzzzzz = pbuffer.data(idx_op_gi + 363);

    auto tr_yzzz_xxxxxx = pbuffer.data(idx_op_gi + 364);

    auto tr_yzzz_xxxxxy = pbuffer.data(idx_op_gi + 365);

    auto tr_yzzz_xxxxxz = pbuffer.data(idx_op_gi + 366);

    auto tr_yzzz_xxxxyy = pbuffer.data(idx_op_gi + 367);

    auto tr_yzzz_xxxxyz = pbuffer.data(idx_op_gi + 368);

    auto tr_yzzz_xxxxzz = pbuffer.data(idx_op_gi + 369);

    auto tr_yzzz_xxxyyy = pbuffer.data(idx_op_gi + 370);

    auto tr_yzzz_xxxyyz = pbuffer.data(idx_op_gi + 371);

    auto tr_yzzz_xxxyzz = pbuffer.data(idx_op_gi + 372);

    auto tr_yzzz_xxxzzz = pbuffer.data(idx_op_gi + 373);

    auto tr_yzzz_xxyyyy = pbuffer.data(idx_op_gi + 374);

    auto tr_yzzz_xxyyyz = pbuffer.data(idx_op_gi + 375);

    auto tr_yzzz_xxyyzz = pbuffer.data(idx_op_gi + 376);

    auto tr_yzzz_xxyzzz = pbuffer.data(idx_op_gi + 377);

    auto tr_yzzz_xxzzzz = pbuffer.data(idx_op_gi + 378);

    auto tr_yzzz_xyyyyy = pbuffer.data(idx_op_gi + 379);

    auto tr_yzzz_xyyyyz = pbuffer.data(idx_op_gi + 380);

    auto tr_yzzz_xyyyzz = pbuffer.data(idx_op_gi + 381);

    auto tr_yzzz_xyyzzz = pbuffer.data(idx_op_gi + 382);

    auto tr_yzzz_xyzzzz = pbuffer.data(idx_op_gi + 383);

    auto tr_yzzz_xzzzzz = pbuffer.data(idx_op_gi + 384);

    auto tr_yzzz_yyyyyy = pbuffer.data(idx_op_gi + 385);

    auto tr_yzzz_yyyyyz = pbuffer.data(idx_op_gi + 386);

    auto tr_yzzz_yyyyzz = pbuffer.data(idx_op_gi + 387);

    auto tr_yzzz_yyyzzz = pbuffer.data(idx_op_gi + 388);

    auto tr_yzzz_yyzzzz = pbuffer.data(idx_op_gi + 389);

    auto tr_yzzz_yzzzzz = pbuffer.data(idx_op_gi + 390);

    auto tr_yzzz_zzzzzz = pbuffer.data(idx_op_gi + 391);

    auto tr_zzzz_xxxxxx = pbuffer.data(idx_op_gi + 392);

    auto tr_zzzz_xxxxxy = pbuffer.data(idx_op_gi + 393);

    auto tr_zzzz_xxxxxz = pbuffer.data(idx_op_gi + 394);

    auto tr_zzzz_xxxxyy = pbuffer.data(idx_op_gi + 395);

    auto tr_zzzz_xxxxyz = pbuffer.data(idx_op_gi + 396);

    auto tr_zzzz_xxxxzz = pbuffer.data(idx_op_gi + 397);

    auto tr_zzzz_xxxyyy = pbuffer.data(idx_op_gi + 398);

    auto tr_zzzz_xxxyyz = pbuffer.data(idx_op_gi + 399);

    auto tr_zzzz_xxxyzz = pbuffer.data(idx_op_gi + 400);

    auto tr_zzzz_xxxzzz = pbuffer.data(idx_op_gi + 401);

    auto tr_zzzz_xxyyyy = pbuffer.data(idx_op_gi + 402);

    auto tr_zzzz_xxyyyz = pbuffer.data(idx_op_gi + 403);

    auto tr_zzzz_xxyyzz = pbuffer.data(idx_op_gi + 404);

    auto tr_zzzz_xxyzzz = pbuffer.data(idx_op_gi + 405);

    auto tr_zzzz_xxzzzz = pbuffer.data(idx_op_gi + 406);

    auto tr_zzzz_xyyyyy = pbuffer.data(idx_op_gi + 407);

    auto tr_zzzz_xyyyyz = pbuffer.data(idx_op_gi + 408);

    auto tr_zzzz_xyyyzz = pbuffer.data(idx_op_gi + 409);

    auto tr_zzzz_xyyzzz = pbuffer.data(idx_op_gi + 410);

    auto tr_zzzz_xyzzzz = pbuffer.data(idx_op_gi + 411);

    auto tr_zzzz_xzzzzz = pbuffer.data(idx_op_gi + 412);

    auto tr_zzzz_yyyyyy = pbuffer.data(idx_op_gi + 413);

    auto tr_zzzz_yyyyyz = pbuffer.data(idx_op_gi + 414);

    auto tr_zzzz_yyyyzz = pbuffer.data(idx_op_gi + 415);

    auto tr_zzzz_yyyzzz = pbuffer.data(idx_op_gi + 416);

    auto tr_zzzz_yyzzzz = pbuffer.data(idx_op_gi + 417);

    auto tr_zzzz_yzzzzz = pbuffer.data(idx_op_gi + 418);

    auto tr_zzzz_zzzzzz = pbuffer.data(idx_op_gi + 419);

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

    // Set up components of auxiliary buffer : HH

    auto tr_xxxxx_xxxxx = pbuffer.data(idx_op_hh);

    auto tr_xxxxx_xxxxy = pbuffer.data(idx_op_hh + 1);

    auto tr_xxxxx_xxxxz = pbuffer.data(idx_op_hh + 2);

    auto tr_xxxxx_xxxyy = pbuffer.data(idx_op_hh + 3);

    auto tr_xxxxx_xxxyz = pbuffer.data(idx_op_hh + 4);

    auto tr_xxxxx_xxxzz = pbuffer.data(idx_op_hh + 5);

    auto tr_xxxxx_xxyyy = pbuffer.data(idx_op_hh + 6);

    auto tr_xxxxx_xxyyz = pbuffer.data(idx_op_hh + 7);

    auto tr_xxxxx_xxyzz = pbuffer.data(idx_op_hh + 8);

    auto tr_xxxxx_xxzzz = pbuffer.data(idx_op_hh + 9);

    auto tr_xxxxx_xyyyy = pbuffer.data(idx_op_hh + 10);

    auto tr_xxxxx_xyyyz = pbuffer.data(idx_op_hh + 11);

    auto tr_xxxxx_xyyzz = pbuffer.data(idx_op_hh + 12);

    auto tr_xxxxx_xyzzz = pbuffer.data(idx_op_hh + 13);

    auto tr_xxxxx_xzzzz = pbuffer.data(idx_op_hh + 14);

    auto tr_xxxxx_yyyyy = pbuffer.data(idx_op_hh + 15);

    auto tr_xxxxx_yyyyz = pbuffer.data(idx_op_hh + 16);

    auto tr_xxxxx_yyyzz = pbuffer.data(idx_op_hh + 17);

    auto tr_xxxxx_yyzzz = pbuffer.data(idx_op_hh + 18);

    auto tr_xxxxx_yzzzz = pbuffer.data(idx_op_hh + 19);

    auto tr_xxxxx_zzzzz = pbuffer.data(idx_op_hh + 20);

    auto tr_xxxxy_xxxxx = pbuffer.data(idx_op_hh + 21);

    auto tr_xxxxy_xxxxy = pbuffer.data(idx_op_hh + 22);

    auto tr_xxxxy_xxxxz = pbuffer.data(idx_op_hh + 23);

    auto tr_xxxxy_xxxyy = pbuffer.data(idx_op_hh + 24);

    auto tr_xxxxy_xxxyz = pbuffer.data(idx_op_hh + 25);

    auto tr_xxxxy_xxxzz = pbuffer.data(idx_op_hh + 26);

    auto tr_xxxxy_xxyyy = pbuffer.data(idx_op_hh + 27);

    auto tr_xxxxy_xxyyz = pbuffer.data(idx_op_hh + 28);

    auto tr_xxxxy_xxyzz = pbuffer.data(idx_op_hh + 29);

    auto tr_xxxxy_xxzzz = pbuffer.data(idx_op_hh + 30);

    auto tr_xxxxy_xyyyy = pbuffer.data(idx_op_hh + 31);

    auto tr_xxxxy_xyyyz = pbuffer.data(idx_op_hh + 32);

    auto tr_xxxxy_xyyzz = pbuffer.data(idx_op_hh + 33);

    auto tr_xxxxy_xyzzz = pbuffer.data(idx_op_hh + 34);

    auto tr_xxxxy_xzzzz = pbuffer.data(idx_op_hh + 35);

    auto tr_xxxxy_yyyyy = pbuffer.data(idx_op_hh + 36);

    auto tr_xxxxy_yyyyz = pbuffer.data(idx_op_hh + 37);

    auto tr_xxxxy_yyyzz = pbuffer.data(idx_op_hh + 38);

    auto tr_xxxxy_yyzzz = pbuffer.data(idx_op_hh + 39);

    auto tr_xxxxy_yzzzz = pbuffer.data(idx_op_hh + 40);

    auto tr_xxxxy_zzzzz = pbuffer.data(idx_op_hh + 41);

    auto tr_xxxxz_xxxxx = pbuffer.data(idx_op_hh + 42);

    auto tr_xxxxz_xxxxy = pbuffer.data(idx_op_hh + 43);

    auto tr_xxxxz_xxxxz = pbuffer.data(idx_op_hh + 44);

    auto tr_xxxxz_xxxyy = pbuffer.data(idx_op_hh + 45);

    auto tr_xxxxz_xxxyz = pbuffer.data(idx_op_hh + 46);

    auto tr_xxxxz_xxxzz = pbuffer.data(idx_op_hh + 47);

    auto tr_xxxxz_xxyyy = pbuffer.data(idx_op_hh + 48);

    auto tr_xxxxz_xxyyz = pbuffer.data(idx_op_hh + 49);

    auto tr_xxxxz_xxyzz = pbuffer.data(idx_op_hh + 50);

    auto tr_xxxxz_xxzzz = pbuffer.data(idx_op_hh + 51);

    auto tr_xxxxz_xyyyy = pbuffer.data(idx_op_hh + 52);

    auto tr_xxxxz_xyyyz = pbuffer.data(idx_op_hh + 53);

    auto tr_xxxxz_xyyzz = pbuffer.data(idx_op_hh + 54);

    auto tr_xxxxz_xyzzz = pbuffer.data(idx_op_hh + 55);

    auto tr_xxxxz_xzzzz = pbuffer.data(idx_op_hh + 56);

    auto tr_xxxxz_yyyyy = pbuffer.data(idx_op_hh + 57);

    auto tr_xxxxz_yyyyz = pbuffer.data(idx_op_hh + 58);

    auto tr_xxxxz_yyyzz = pbuffer.data(idx_op_hh + 59);

    auto tr_xxxxz_yyzzz = pbuffer.data(idx_op_hh + 60);

    auto tr_xxxxz_yzzzz = pbuffer.data(idx_op_hh + 61);

    auto tr_xxxxz_zzzzz = pbuffer.data(idx_op_hh + 62);

    auto tr_xxxyy_xxxxx = pbuffer.data(idx_op_hh + 63);

    auto tr_xxxyy_xxxxy = pbuffer.data(idx_op_hh + 64);

    auto tr_xxxyy_xxxxz = pbuffer.data(idx_op_hh + 65);

    auto tr_xxxyy_xxxyy = pbuffer.data(idx_op_hh + 66);

    auto tr_xxxyy_xxxyz = pbuffer.data(idx_op_hh + 67);

    auto tr_xxxyy_xxxzz = pbuffer.data(idx_op_hh + 68);

    auto tr_xxxyy_xxyyy = pbuffer.data(idx_op_hh + 69);

    auto tr_xxxyy_xxyyz = pbuffer.data(idx_op_hh + 70);

    auto tr_xxxyy_xxyzz = pbuffer.data(idx_op_hh + 71);

    auto tr_xxxyy_xxzzz = pbuffer.data(idx_op_hh + 72);

    auto tr_xxxyy_xyyyy = pbuffer.data(idx_op_hh + 73);

    auto tr_xxxyy_xyyyz = pbuffer.data(idx_op_hh + 74);

    auto tr_xxxyy_xyyzz = pbuffer.data(idx_op_hh + 75);

    auto tr_xxxyy_xyzzz = pbuffer.data(idx_op_hh + 76);

    auto tr_xxxyy_xzzzz = pbuffer.data(idx_op_hh + 77);

    auto tr_xxxyy_yyyyy = pbuffer.data(idx_op_hh + 78);

    auto tr_xxxyy_yyyyz = pbuffer.data(idx_op_hh + 79);

    auto tr_xxxyy_yyyzz = pbuffer.data(idx_op_hh + 80);

    auto tr_xxxyy_yyzzz = pbuffer.data(idx_op_hh + 81);

    auto tr_xxxyy_yzzzz = pbuffer.data(idx_op_hh + 82);

    auto tr_xxxyy_zzzzz = pbuffer.data(idx_op_hh + 83);

    auto tr_xxxyz_xxxxx = pbuffer.data(idx_op_hh + 84);

    auto tr_xxxyz_xxxxy = pbuffer.data(idx_op_hh + 85);

    auto tr_xxxyz_xxxxz = pbuffer.data(idx_op_hh + 86);

    auto tr_xxxyz_xxxyy = pbuffer.data(idx_op_hh + 87);

    auto tr_xxxyz_xxxyz = pbuffer.data(idx_op_hh + 88);

    auto tr_xxxyz_xxxzz = pbuffer.data(idx_op_hh + 89);

    auto tr_xxxyz_xxyyy = pbuffer.data(idx_op_hh + 90);

    auto tr_xxxyz_xxyyz = pbuffer.data(idx_op_hh + 91);

    auto tr_xxxyz_xxyzz = pbuffer.data(idx_op_hh + 92);

    auto tr_xxxyz_xxzzz = pbuffer.data(idx_op_hh + 93);

    auto tr_xxxyz_xyyyy = pbuffer.data(idx_op_hh + 94);

    auto tr_xxxyz_xyyyz = pbuffer.data(idx_op_hh + 95);

    auto tr_xxxyz_xyyzz = pbuffer.data(idx_op_hh + 96);

    auto tr_xxxyz_xyzzz = pbuffer.data(idx_op_hh + 97);

    auto tr_xxxyz_xzzzz = pbuffer.data(idx_op_hh + 98);

    auto tr_xxxyz_yyyyy = pbuffer.data(idx_op_hh + 99);

    auto tr_xxxyz_yyyyz = pbuffer.data(idx_op_hh + 100);

    auto tr_xxxyz_yyyzz = pbuffer.data(idx_op_hh + 101);

    auto tr_xxxyz_yyzzz = pbuffer.data(idx_op_hh + 102);

    auto tr_xxxyz_yzzzz = pbuffer.data(idx_op_hh + 103);

    auto tr_xxxyz_zzzzz = pbuffer.data(idx_op_hh + 104);

    auto tr_xxxzz_xxxxx = pbuffer.data(idx_op_hh + 105);

    auto tr_xxxzz_xxxxy = pbuffer.data(idx_op_hh + 106);

    auto tr_xxxzz_xxxxz = pbuffer.data(idx_op_hh + 107);

    auto tr_xxxzz_xxxyy = pbuffer.data(idx_op_hh + 108);

    auto tr_xxxzz_xxxyz = pbuffer.data(idx_op_hh + 109);

    auto tr_xxxzz_xxxzz = pbuffer.data(idx_op_hh + 110);

    auto tr_xxxzz_xxyyy = pbuffer.data(idx_op_hh + 111);

    auto tr_xxxzz_xxyyz = pbuffer.data(idx_op_hh + 112);

    auto tr_xxxzz_xxyzz = pbuffer.data(idx_op_hh + 113);

    auto tr_xxxzz_xxzzz = pbuffer.data(idx_op_hh + 114);

    auto tr_xxxzz_xyyyy = pbuffer.data(idx_op_hh + 115);

    auto tr_xxxzz_xyyyz = pbuffer.data(idx_op_hh + 116);

    auto tr_xxxzz_xyyzz = pbuffer.data(idx_op_hh + 117);

    auto tr_xxxzz_xyzzz = pbuffer.data(idx_op_hh + 118);

    auto tr_xxxzz_xzzzz = pbuffer.data(idx_op_hh + 119);

    auto tr_xxxzz_yyyyy = pbuffer.data(idx_op_hh + 120);

    auto tr_xxxzz_yyyyz = pbuffer.data(idx_op_hh + 121);

    auto tr_xxxzz_yyyzz = pbuffer.data(idx_op_hh + 122);

    auto tr_xxxzz_yyzzz = pbuffer.data(idx_op_hh + 123);

    auto tr_xxxzz_yzzzz = pbuffer.data(idx_op_hh + 124);

    auto tr_xxxzz_zzzzz = pbuffer.data(idx_op_hh + 125);

    auto tr_xxyyy_xxxxx = pbuffer.data(idx_op_hh + 126);

    auto tr_xxyyy_xxxxy = pbuffer.data(idx_op_hh + 127);

    auto tr_xxyyy_xxxxz = pbuffer.data(idx_op_hh + 128);

    auto tr_xxyyy_xxxyy = pbuffer.data(idx_op_hh + 129);

    auto tr_xxyyy_xxxyz = pbuffer.data(idx_op_hh + 130);

    auto tr_xxyyy_xxxzz = pbuffer.data(idx_op_hh + 131);

    auto tr_xxyyy_xxyyy = pbuffer.data(idx_op_hh + 132);

    auto tr_xxyyy_xxyyz = pbuffer.data(idx_op_hh + 133);

    auto tr_xxyyy_xxyzz = pbuffer.data(idx_op_hh + 134);

    auto tr_xxyyy_xxzzz = pbuffer.data(idx_op_hh + 135);

    auto tr_xxyyy_xyyyy = pbuffer.data(idx_op_hh + 136);

    auto tr_xxyyy_xyyyz = pbuffer.data(idx_op_hh + 137);

    auto tr_xxyyy_xyyzz = pbuffer.data(idx_op_hh + 138);

    auto tr_xxyyy_xyzzz = pbuffer.data(idx_op_hh + 139);

    auto tr_xxyyy_xzzzz = pbuffer.data(idx_op_hh + 140);

    auto tr_xxyyy_yyyyy = pbuffer.data(idx_op_hh + 141);

    auto tr_xxyyy_yyyyz = pbuffer.data(idx_op_hh + 142);

    auto tr_xxyyy_yyyzz = pbuffer.data(idx_op_hh + 143);

    auto tr_xxyyy_yyzzz = pbuffer.data(idx_op_hh + 144);

    auto tr_xxyyy_yzzzz = pbuffer.data(idx_op_hh + 145);

    auto tr_xxyyy_zzzzz = pbuffer.data(idx_op_hh + 146);

    auto tr_xxyyz_xxxxx = pbuffer.data(idx_op_hh + 147);

    auto tr_xxyyz_xxxxy = pbuffer.data(idx_op_hh + 148);

    auto tr_xxyyz_xxxxz = pbuffer.data(idx_op_hh + 149);

    auto tr_xxyyz_xxxyy = pbuffer.data(idx_op_hh + 150);

    auto tr_xxyyz_xxxyz = pbuffer.data(idx_op_hh + 151);

    auto tr_xxyyz_xxxzz = pbuffer.data(idx_op_hh + 152);

    auto tr_xxyyz_xxyyy = pbuffer.data(idx_op_hh + 153);

    auto tr_xxyyz_xxyyz = pbuffer.data(idx_op_hh + 154);

    auto tr_xxyyz_xxyzz = pbuffer.data(idx_op_hh + 155);

    auto tr_xxyyz_xxzzz = pbuffer.data(idx_op_hh + 156);

    auto tr_xxyyz_xyyyy = pbuffer.data(idx_op_hh + 157);

    auto tr_xxyyz_xyyyz = pbuffer.data(idx_op_hh + 158);

    auto tr_xxyyz_xyyzz = pbuffer.data(idx_op_hh + 159);

    auto tr_xxyyz_xyzzz = pbuffer.data(idx_op_hh + 160);

    auto tr_xxyyz_xzzzz = pbuffer.data(idx_op_hh + 161);

    auto tr_xxyyz_yyyyy = pbuffer.data(idx_op_hh + 162);

    auto tr_xxyyz_yyyyz = pbuffer.data(idx_op_hh + 163);

    auto tr_xxyyz_yyyzz = pbuffer.data(idx_op_hh + 164);

    auto tr_xxyyz_yyzzz = pbuffer.data(idx_op_hh + 165);

    auto tr_xxyyz_yzzzz = pbuffer.data(idx_op_hh + 166);

    auto tr_xxyyz_zzzzz = pbuffer.data(idx_op_hh + 167);

    auto tr_xxyzz_xxxxx = pbuffer.data(idx_op_hh + 168);

    auto tr_xxyzz_xxxxy = pbuffer.data(idx_op_hh + 169);

    auto tr_xxyzz_xxxxz = pbuffer.data(idx_op_hh + 170);

    auto tr_xxyzz_xxxyy = pbuffer.data(idx_op_hh + 171);

    auto tr_xxyzz_xxxyz = pbuffer.data(idx_op_hh + 172);

    auto tr_xxyzz_xxxzz = pbuffer.data(idx_op_hh + 173);

    auto tr_xxyzz_xxyyy = pbuffer.data(idx_op_hh + 174);

    auto tr_xxyzz_xxyyz = pbuffer.data(idx_op_hh + 175);

    auto tr_xxyzz_xxyzz = pbuffer.data(idx_op_hh + 176);

    auto tr_xxyzz_xxzzz = pbuffer.data(idx_op_hh + 177);

    auto tr_xxyzz_xyyyy = pbuffer.data(idx_op_hh + 178);

    auto tr_xxyzz_xyyyz = pbuffer.data(idx_op_hh + 179);

    auto tr_xxyzz_xyyzz = pbuffer.data(idx_op_hh + 180);

    auto tr_xxyzz_xyzzz = pbuffer.data(idx_op_hh + 181);

    auto tr_xxyzz_xzzzz = pbuffer.data(idx_op_hh + 182);

    auto tr_xxyzz_yyyyy = pbuffer.data(idx_op_hh + 183);

    auto tr_xxyzz_yyyyz = pbuffer.data(idx_op_hh + 184);

    auto tr_xxyzz_yyyzz = pbuffer.data(idx_op_hh + 185);

    auto tr_xxyzz_yyzzz = pbuffer.data(idx_op_hh + 186);

    auto tr_xxyzz_yzzzz = pbuffer.data(idx_op_hh + 187);

    auto tr_xxyzz_zzzzz = pbuffer.data(idx_op_hh + 188);

    auto tr_xxzzz_xxxxx = pbuffer.data(idx_op_hh + 189);

    auto tr_xxzzz_xxxxy = pbuffer.data(idx_op_hh + 190);

    auto tr_xxzzz_xxxxz = pbuffer.data(idx_op_hh + 191);

    auto tr_xxzzz_xxxyy = pbuffer.data(idx_op_hh + 192);

    auto tr_xxzzz_xxxyz = pbuffer.data(idx_op_hh + 193);

    auto tr_xxzzz_xxxzz = pbuffer.data(idx_op_hh + 194);

    auto tr_xxzzz_xxyyy = pbuffer.data(idx_op_hh + 195);

    auto tr_xxzzz_xxyyz = pbuffer.data(idx_op_hh + 196);

    auto tr_xxzzz_xxyzz = pbuffer.data(idx_op_hh + 197);

    auto tr_xxzzz_xxzzz = pbuffer.data(idx_op_hh + 198);

    auto tr_xxzzz_xyyyy = pbuffer.data(idx_op_hh + 199);

    auto tr_xxzzz_xyyyz = pbuffer.data(idx_op_hh + 200);

    auto tr_xxzzz_xyyzz = pbuffer.data(idx_op_hh + 201);

    auto tr_xxzzz_xyzzz = pbuffer.data(idx_op_hh + 202);

    auto tr_xxzzz_xzzzz = pbuffer.data(idx_op_hh + 203);

    auto tr_xxzzz_yyyyy = pbuffer.data(idx_op_hh + 204);

    auto tr_xxzzz_yyyyz = pbuffer.data(idx_op_hh + 205);

    auto tr_xxzzz_yyyzz = pbuffer.data(idx_op_hh + 206);

    auto tr_xxzzz_yyzzz = pbuffer.data(idx_op_hh + 207);

    auto tr_xxzzz_yzzzz = pbuffer.data(idx_op_hh + 208);

    auto tr_xxzzz_zzzzz = pbuffer.data(idx_op_hh + 209);

    auto tr_xyyyy_xxxxx = pbuffer.data(idx_op_hh + 210);

    auto tr_xyyyy_xxxxy = pbuffer.data(idx_op_hh + 211);

    auto tr_xyyyy_xxxxz = pbuffer.data(idx_op_hh + 212);

    auto tr_xyyyy_xxxyy = pbuffer.data(idx_op_hh + 213);

    auto tr_xyyyy_xxxyz = pbuffer.data(idx_op_hh + 214);

    auto tr_xyyyy_xxxzz = pbuffer.data(idx_op_hh + 215);

    auto tr_xyyyy_xxyyy = pbuffer.data(idx_op_hh + 216);

    auto tr_xyyyy_xxyyz = pbuffer.data(idx_op_hh + 217);

    auto tr_xyyyy_xxyzz = pbuffer.data(idx_op_hh + 218);

    auto tr_xyyyy_xxzzz = pbuffer.data(idx_op_hh + 219);

    auto tr_xyyyy_xyyyy = pbuffer.data(idx_op_hh + 220);

    auto tr_xyyyy_xyyyz = pbuffer.data(idx_op_hh + 221);

    auto tr_xyyyy_xyyzz = pbuffer.data(idx_op_hh + 222);

    auto tr_xyyyy_xyzzz = pbuffer.data(idx_op_hh + 223);

    auto tr_xyyyy_xzzzz = pbuffer.data(idx_op_hh + 224);

    auto tr_xyyyy_yyyyy = pbuffer.data(idx_op_hh + 225);

    auto tr_xyyyy_yyyyz = pbuffer.data(idx_op_hh + 226);

    auto tr_xyyyy_yyyzz = pbuffer.data(idx_op_hh + 227);

    auto tr_xyyyy_yyzzz = pbuffer.data(idx_op_hh + 228);

    auto tr_xyyyy_yzzzz = pbuffer.data(idx_op_hh + 229);

    auto tr_xyyyy_zzzzz = pbuffer.data(idx_op_hh + 230);

    auto tr_xyyyz_xxxxx = pbuffer.data(idx_op_hh + 231);

    auto tr_xyyyz_xxxxy = pbuffer.data(idx_op_hh + 232);

    auto tr_xyyyz_xxxxz = pbuffer.data(idx_op_hh + 233);

    auto tr_xyyyz_xxxyy = pbuffer.data(idx_op_hh + 234);

    auto tr_xyyyz_xxxyz = pbuffer.data(idx_op_hh + 235);

    auto tr_xyyyz_xxxzz = pbuffer.data(idx_op_hh + 236);

    auto tr_xyyyz_xxyyy = pbuffer.data(idx_op_hh + 237);

    auto tr_xyyyz_xxyyz = pbuffer.data(idx_op_hh + 238);

    auto tr_xyyyz_xxyzz = pbuffer.data(idx_op_hh + 239);

    auto tr_xyyyz_xxzzz = pbuffer.data(idx_op_hh + 240);

    auto tr_xyyyz_xyyyy = pbuffer.data(idx_op_hh + 241);

    auto tr_xyyyz_xyyyz = pbuffer.data(idx_op_hh + 242);

    auto tr_xyyyz_xyyzz = pbuffer.data(idx_op_hh + 243);

    auto tr_xyyyz_xyzzz = pbuffer.data(idx_op_hh + 244);

    auto tr_xyyyz_xzzzz = pbuffer.data(idx_op_hh + 245);

    auto tr_xyyyz_yyyyy = pbuffer.data(idx_op_hh + 246);

    auto tr_xyyyz_yyyyz = pbuffer.data(idx_op_hh + 247);

    auto tr_xyyyz_yyyzz = pbuffer.data(idx_op_hh + 248);

    auto tr_xyyyz_yyzzz = pbuffer.data(idx_op_hh + 249);

    auto tr_xyyyz_yzzzz = pbuffer.data(idx_op_hh + 250);

    auto tr_xyyyz_zzzzz = pbuffer.data(idx_op_hh + 251);

    auto tr_xyyzz_xxxxx = pbuffer.data(idx_op_hh + 252);

    auto tr_xyyzz_xxxxy = pbuffer.data(idx_op_hh + 253);

    auto tr_xyyzz_xxxxz = pbuffer.data(idx_op_hh + 254);

    auto tr_xyyzz_xxxyy = pbuffer.data(idx_op_hh + 255);

    auto tr_xyyzz_xxxyz = pbuffer.data(idx_op_hh + 256);

    auto tr_xyyzz_xxxzz = pbuffer.data(idx_op_hh + 257);

    auto tr_xyyzz_xxyyy = pbuffer.data(idx_op_hh + 258);

    auto tr_xyyzz_xxyyz = pbuffer.data(idx_op_hh + 259);

    auto tr_xyyzz_xxyzz = pbuffer.data(idx_op_hh + 260);

    auto tr_xyyzz_xxzzz = pbuffer.data(idx_op_hh + 261);

    auto tr_xyyzz_xyyyy = pbuffer.data(idx_op_hh + 262);

    auto tr_xyyzz_xyyyz = pbuffer.data(idx_op_hh + 263);

    auto tr_xyyzz_xyyzz = pbuffer.data(idx_op_hh + 264);

    auto tr_xyyzz_xyzzz = pbuffer.data(idx_op_hh + 265);

    auto tr_xyyzz_xzzzz = pbuffer.data(idx_op_hh + 266);

    auto tr_xyyzz_yyyyy = pbuffer.data(idx_op_hh + 267);

    auto tr_xyyzz_yyyyz = pbuffer.data(idx_op_hh + 268);

    auto tr_xyyzz_yyyzz = pbuffer.data(idx_op_hh + 269);

    auto tr_xyyzz_yyzzz = pbuffer.data(idx_op_hh + 270);

    auto tr_xyyzz_yzzzz = pbuffer.data(idx_op_hh + 271);

    auto tr_xyyzz_zzzzz = pbuffer.data(idx_op_hh + 272);

    auto tr_xyzzz_xxxxx = pbuffer.data(idx_op_hh + 273);

    auto tr_xyzzz_xxxxy = pbuffer.data(idx_op_hh + 274);

    auto tr_xyzzz_xxxxz = pbuffer.data(idx_op_hh + 275);

    auto tr_xyzzz_xxxyy = pbuffer.data(idx_op_hh + 276);

    auto tr_xyzzz_xxxyz = pbuffer.data(idx_op_hh + 277);

    auto tr_xyzzz_xxxzz = pbuffer.data(idx_op_hh + 278);

    auto tr_xyzzz_xxyyy = pbuffer.data(idx_op_hh + 279);

    auto tr_xyzzz_xxyyz = pbuffer.data(idx_op_hh + 280);

    auto tr_xyzzz_xxyzz = pbuffer.data(idx_op_hh + 281);

    auto tr_xyzzz_xxzzz = pbuffer.data(idx_op_hh + 282);

    auto tr_xyzzz_xyyyy = pbuffer.data(idx_op_hh + 283);

    auto tr_xyzzz_xyyyz = pbuffer.data(idx_op_hh + 284);

    auto tr_xyzzz_xyyzz = pbuffer.data(idx_op_hh + 285);

    auto tr_xyzzz_xyzzz = pbuffer.data(idx_op_hh + 286);

    auto tr_xyzzz_xzzzz = pbuffer.data(idx_op_hh + 287);

    auto tr_xyzzz_yyyyy = pbuffer.data(idx_op_hh + 288);

    auto tr_xyzzz_yyyyz = pbuffer.data(idx_op_hh + 289);

    auto tr_xyzzz_yyyzz = pbuffer.data(idx_op_hh + 290);

    auto tr_xyzzz_yyzzz = pbuffer.data(idx_op_hh + 291);

    auto tr_xyzzz_yzzzz = pbuffer.data(idx_op_hh + 292);

    auto tr_xyzzz_zzzzz = pbuffer.data(idx_op_hh + 293);

    auto tr_xzzzz_xxxxx = pbuffer.data(idx_op_hh + 294);

    auto tr_xzzzz_xxxxy = pbuffer.data(idx_op_hh + 295);

    auto tr_xzzzz_xxxxz = pbuffer.data(idx_op_hh + 296);

    auto tr_xzzzz_xxxyy = pbuffer.data(idx_op_hh + 297);

    auto tr_xzzzz_xxxyz = pbuffer.data(idx_op_hh + 298);

    auto tr_xzzzz_xxxzz = pbuffer.data(idx_op_hh + 299);

    auto tr_xzzzz_xxyyy = pbuffer.data(idx_op_hh + 300);

    auto tr_xzzzz_xxyyz = pbuffer.data(idx_op_hh + 301);

    auto tr_xzzzz_xxyzz = pbuffer.data(idx_op_hh + 302);

    auto tr_xzzzz_xxzzz = pbuffer.data(idx_op_hh + 303);

    auto tr_xzzzz_xyyyy = pbuffer.data(idx_op_hh + 304);

    auto tr_xzzzz_xyyyz = pbuffer.data(idx_op_hh + 305);

    auto tr_xzzzz_xyyzz = pbuffer.data(idx_op_hh + 306);

    auto tr_xzzzz_xyzzz = pbuffer.data(idx_op_hh + 307);

    auto tr_xzzzz_xzzzz = pbuffer.data(idx_op_hh + 308);

    auto tr_xzzzz_yyyyy = pbuffer.data(idx_op_hh + 309);

    auto tr_xzzzz_yyyyz = pbuffer.data(idx_op_hh + 310);

    auto tr_xzzzz_yyyzz = pbuffer.data(idx_op_hh + 311);

    auto tr_xzzzz_yyzzz = pbuffer.data(idx_op_hh + 312);

    auto tr_xzzzz_yzzzz = pbuffer.data(idx_op_hh + 313);

    auto tr_xzzzz_zzzzz = pbuffer.data(idx_op_hh + 314);

    auto tr_yyyyy_xxxxx = pbuffer.data(idx_op_hh + 315);

    auto tr_yyyyy_xxxxy = pbuffer.data(idx_op_hh + 316);

    auto tr_yyyyy_xxxxz = pbuffer.data(idx_op_hh + 317);

    auto tr_yyyyy_xxxyy = pbuffer.data(idx_op_hh + 318);

    auto tr_yyyyy_xxxyz = pbuffer.data(idx_op_hh + 319);

    auto tr_yyyyy_xxxzz = pbuffer.data(idx_op_hh + 320);

    auto tr_yyyyy_xxyyy = pbuffer.data(idx_op_hh + 321);

    auto tr_yyyyy_xxyyz = pbuffer.data(idx_op_hh + 322);

    auto tr_yyyyy_xxyzz = pbuffer.data(idx_op_hh + 323);

    auto tr_yyyyy_xxzzz = pbuffer.data(idx_op_hh + 324);

    auto tr_yyyyy_xyyyy = pbuffer.data(idx_op_hh + 325);

    auto tr_yyyyy_xyyyz = pbuffer.data(idx_op_hh + 326);

    auto tr_yyyyy_xyyzz = pbuffer.data(idx_op_hh + 327);

    auto tr_yyyyy_xyzzz = pbuffer.data(idx_op_hh + 328);

    auto tr_yyyyy_xzzzz = pbuffer.data(idx_op_hh + 329);

    auto tr_yyyyy_yyyyy = pbuffer.data(idx_op_hh + 330);

    auto tr_yyyyy_yyyyz = pbuffer.data(idx_op_hh + 331);

    auto tr_yyyyy_yyyzz = pbuffer.data(idx_op_hh + 332);

    auto tr_yyyyy_yyzzz = pbuffer.data(idx_op_hh + 333);

    auto tr_yyyyy_yzzzz = pbuffer.data(idx_op_hh + 334);

    auto tr_yyyyy_zzzzz = pbuffer.data(idx_op_hh + 335);

    auto tr_yyyyz_xxxxx = pbuffer.data(idx_op_hh + 336);

    auto tr_yyyyz_xxxxy = pbuffer.data(idx_op_hh + 337);

    auto tr_yyyyz_xxxxz = pbuffer.data(idx_op_hh + 338);

    auto tr_yyyyz_xxxyy = pbuffer.data(idx_op_hh + 339);

    auto tr_yyyyz_xxxyz = pbuffer.data(idx_op_hh + 340);

    auto tr_yyyyz_xxxzz = pbuffer.data(idx_op_hh + 341);

    auto tr_yyyyz_xxyyy = pbuffer.data(idx_op_hh + 342);

    auto tr_yyyyz_xxyyz = pbuffer.data(idx_op_hh + 343);

    auto tr_yyyyz_xxyzz = pbuffer.data(idx_op_hh + 344);

    auto tr_yyyyz_xxzzz = pbuffer.data(idx_op_hh + 345);

    auto tr_yyyyz_xyyyy = pbuffer.data(idx_op_hh + 346);

    auto tr_yyyyz_xyyyz = pbuffer.data(idx_op_hh + 347);

    auto tr_yyyyz_xyyzz = pbuffer.data(idx_op_hh + 348);

    auto tr_yyyyz_xyzzz = pbuffer.data(idx_op_hh + 349);

    auto tr_yyyyz_xzzzz = pbuffer.data(idx_op_hh + 350);

    auto tr_yyyyz_yyyyy = pbuffer.data(idx_op_hh + 351);

    auto tr_yyyyz_yyyyz = pbuffer.data(idx_op_hh + 352);

    auto tr_yyyyz_yyyzz = pbuffer.data(idx_op_hh + 353);

    auto tr_yyyyz_yyzzz = pbuffer.data(idx_op_hh + 354);

    auto tr_yyyyz_yzzzz = pbuffer.data(idx_op_hh + 355);

    auto tr_yyyyz_zzzzz = pbuffer.data(idx_op_hh + 356);

    auto tr_yyyzz_xxxxx = pbuffer.data(idx_op_hh + 357);

    auto tr_yyyzz_xxxxy = pbuffer.data(idx_op_hh + 358);

    auto tr_yyyzz_xxxxz = pbuffer.data(idx_op_hh + 359);

    auto tr_yyyzz_xxxyy = pbuffer.data(idx_op_hh + 360);

    auto tr_yyyzz_xxxyz = pbuffer.data(idx_op_hh + 361);

    auto tr_yyyzz_xxxzz = pbuffer.data(idx_op_hh + 362);

    auto tr_yyyzz_xxyyy = pbuffer.data(idx_op_hh + 363);

    auto tr_yyyzz_xxyyz = pbuffer.data(idx_op_hh + 364);

    auto tr_yyyzz_xxyzz = pbuffer.data(idx_op_hh + 365);

    auto tr_yyyzz_xxzzz = pbuffer.data(idx_op_hh + 366);

    auto tr_yyyzz_xyyyy = pbuffer.data(idx_op_hh + 367);

    auto tr_yyyzz_xyyyz = pbuffer.data(idx_op_hh + 368);

    auto tr_yyyzz_xyyzz = pbuffer.data(idx_op_hh + 369);

    auto tr_yyyzz_xyzzz = pbuffer.data(idx_op_hh + 370);

    auto tr_yyyzz_xzzzz = pbuffer.data(idx_op_hh + 371);

    auto tr_yyyzz_yyyyy = pbuffer.data(idx_op_hh + 372);

    auto tr_yyyzz_yyyyz = pbuffer.data(idx_op_hh + 373);

    auto tr_yyyzz_yyyzz = pbuffer.data(idx_op_hh + 374);

    auto tr_yyyzz_yyzzz = pbuffer.data(idx_op_hh + 375);

    auto tr_yyyzz_yzzzz = pbuffer.data(idx_op_hh + 376);

    auto tr_yyyzz_zzzzz = pbuffer.data(idx_op_hh + 377);

    auto tr_yyzzz_xxxxx = pbuffer.data(idx_op_hh + 378);

    auto tr_yyzzz_xxxxy = pbuffer.data(idx_op_hh + 379);

    auto tr_yyzzz_xxxxz = pbuffer.data(idx_op_hh + 380);

    auto tr_yyzzz_xxxyy = pbuffer.data(idx_op_hh + 381);

    auto tr_yyzzz_xxxyz = pbuffer.data(idx_op_hh + 382);

    auto tr_yyzzz_xxxzz = pbuffer.data(idx_op_hh + 383);

    auto tr_yyzzz_xxyyy = pbuffer.data(idx_op_hh + 384);

    auto tr_yyzzz_xxyyz = pbuffer.data(idx_op_hh + 385);

    auto tr_yyzzz_xxyzz = pbuffer.data(idx_op_hh + 386);

    auto tr_yyzzz_xxzzz = pbuffer.data(idx_op_hh + 387);

    auto tr_yyzzz_xyyyy = pbuffer.data(idx_op_hh + 388);

    auto tr_yyzzz_xyyyz = pbuffer.data(idx_op_hh + 389);

    auto tr_yyzzz_xyyzz = pbuffer.data(idx_op_hh + 390);

    auto tr_yyzzz_xyzzz = pbuffer.data(idx_op_hh + 391);

    auto tr_yyzzz_xzzzz = pbuffer.data(idx_op_hh + 392);

    auto tr_yyzzz_yyyyy = pbuffer.data(idx_op_hh + 393);

    auto tr_yyzzz_yyyyz = pbuffer.data(idx_op_hh + 394);

    auto tr_yyzzz_yyyzz = pbuffer.data(idx_op_hh + 395);

    auto tr_yyzzz_yyzzz = pbuffer.data(idx_op_hh + 396);

    auto tr_yyzzz_yzzzz = pbuffer.data(idx_op_hh + 397);

    auto tr_yyzzz_zzzzz = pbuffer.data(idx_op_hh + 398);

    auto tr_yzzzz_xxxxx = pbuffer.data(idx_op_hh + 399);

    auto tr_yzzzz_xxxxy = pbuffer.data(idx_op_hh + 400);

    auto tr_yzzzz_xxxxz = pbuffer.data(idx_op_hh + 401);

    auto tr_yzzzz_xxxyy = pbuffer.data(idx_op_hh + 402);

    auto tr_yzzzz_xxxyz = pbuffer.data(idx_op_hh + 403);

    auto tr_yzzzz_xxxzz = pbuffer.data(idx_op_hh + 404);

    auto tr_yzzzz_xxyyy = pbuffer.data(idx_op_hh + 405);

    auto tr_yzzzz_xxyyz = pbuffer.data(idx_op_hh + 406);

    auto tr_yzzzz_xxyzz = pbuffer.data(idx_op_hh + 407);

    auto tr_yzzzz_xxzzz = pbuffer.data(idx_op_hh + 408);

    auto tr_yzzzz_xyyyy = pbuffer.data(idx_op_hh + 409);

    auto tr_yzzzz_xyyyz = pbuffer.data(idx_op_hh + 410);

    auto tr_yzzzz_xyyzz = pbuffer.data(idx_op_hh + 411);

    auto tr_yzzzz_xyzzz = pbuffer.data(idx_op_hh + 412);

    auto tr_yzzzz_xzzzz = pbuffer.data(idx_op_hh + 413);

    auto tr_yzzzz_yyyyy = pbuffer.data(idx_op_hh + 414);

    auto tr_yzzzz_yyyyz = pbuffer.data(idx_op_hh + 415);

    auto tr_yzzzz_yyyzz = pbuffer.data(idx_op_hh + 416);

    auto tr_yzzzz_yyzzz = pbuffer.data(idx_op_hh + 417);

    auto tr_yzzzz_yzzzz = pbuffer.data(idx_op_hh + 418);

    auto tr_yzzzz_zzzzz = pbuffer.data(idx_op_hh + 419);

    auto tr_zzzzz_xxxxx = pbuffer.data(idx_op_hh + 420);

    auto tr_zzzzz_xxxxy = pbuffer.data(idx_op_hh + 421);

    auto tr_zzzzz_xxxxz = pbuffer.data(idx_op_hh + 422);

    auto tr_zzzzz_xxxyy = pbuffer.data(idx_op_hh + 423);

    auto tr_zzzzz_xxxyz = pbuffer.data(idx_op_hh + 424);

    auto tr_zzzzz_xxxzz = pbuffer.data(idx_op_hh + 425);

    auto tr_zzzzz_xxyyy = pbuffer.data(idx_op_hh + 426);

    auto tr_zzzzz_xxyyz = pbuffer.data(idx_op_hh + 427);

    auto tr_zzzzz_xxyzz = pbuffer.data(idx_op_hh + 428);

    auto tr_zzzzz_xxzzz = pbuffer.data(idx_op_hh + 429);

    auto tr_zzzzz_xyyyy = pbuffer.data(idx_op_hh + 430);

    auto tr_zzzzz_xyyyz = pbuffer.data(idx_op_hh + 431);

    auto tr_zzzzz_xyyzz = pbuffer.data(idx_op_hh + 432);

    auto tr_zzzzz_xyzzz = pbuffer.data(idx_op_hh + 433);

    auto tr_zzzzz_xzzzz = pbuffer.data(idx_op_hh + 434);

    auto tr_zzzzz_yyyyy = pbuffer.data(idx_op_hh + 435);

    auto tr_zzzzz_yyyyz = pbuffer.data(idx_op_hh + 436);

    auto tr_zzzzz_yyyzz = pbuffer.data(idx_op_hh + 437);

    auto tr_zzzzz_yyzzz = pbuffer.data(idx_op_hh + 438);

    auto tr_zzzzz_yzzzz = pbuffer.data(idx_op_hh + 439);

    auto tr_zzzzz_zzzzz = pbuffer.data(idx_op_hh + 440);

    // Set up components of auxiliary buffer : IG

    auto tr_xxxxxx_xxxx = pbuffer.data(idx_op_ig);

    auto tr_xxxxxx_xxxy = pbuffer.data(idx_op_ig + 1);

    auto tr_xxxxxx_xxxz = pbuffer.data(idx_op_ig + 2);

    auto tr_xxxxxx_xxyy = pbuffer.data(idx_op_ig + 3);

    auto tr_xxxxxx_xxyz = pbuffer.data(idx_op_ig + 4);

    auto tr_xxxxxx_xxzz = pbuffer.data(idx_op_ig + 5);

    auto tr_xxxxxx_xyyy = pbuffer.data(idx_op_ig + 6);

    auto tr_xxxxxx_xyyz = pbuffer.data(idx_op_ig + 7);

    auto tr_xxxxxx_xyzz = pbuffer.data(idx_op_ig + 8);

    auto tr_xxxxxx_xzzz = pbuffer.data(idx_op_ig + 9);

    auto tr_xxxxxx_yyyy = pbuffer.data(idx_op_ig + 10);

    auto tr_xxxxxx_yyyz = pbuffer.data(idx_op_ig + 11);

    auto tr_xxxxxx_yyzz = pbuffer.data(idx_op_ig + 12);

    auto tr_xxxxxx_yzzz = pbuffer.data(idx_op_ig + 13);

    auto tr_xxxxxx_zzzz = pbuffer.data(idx_op_ig + 14);

    auto tr_xxxxxy_xxxx = pbuffer.data(idx_op_ig + 15);

    auto tr_xxxxxy_xxxy = pbuffer.data(idx_op_ig + 16);

    auto tr_xxxxxy_xxxz = pbuffer.data(idx_op_ig + 17);

    auto tr_xxxxxy_xxyy = pbuffer.data(idx_op_ig + 18);

    auto tr_xxxxxy_xxyz = pbuffer.data(idx_op_ig + 19);

    auto tr_xxxxxy_xxzz = pbuffer.data(idx_op_ig + 20);

    auto tr_xxxxxy_xyyy = pbuffer.data(idx_op_ig + 21);

    auto tr_xxxxxy_xyyz = pbuffer.data(idx_op_ig + 22);

    auto tr_xxxxxy_xyzz = pbuffer.data(idx_op_ig + 23);

    auto tr_xxxxxy_xzzz = pbuffer.data(idx_op_ig + 24);

    auto tr_xxxxxy_yyyy = pbuffer.data(idx_op_ig + 25);

    auto tr_xxxxxy_yyyz = pbuffer.data(idx_op_ig + 26);

    auto tr_xxxxxy_yyzz = pbuffer.data(idx_op_ig + 27);

    auto tr_xxxxxy_yzzz = pbuffer.data(idx_op_ig + 28);

    auto tr_xxxxxy_zzzz = pbuffer.data(idx_op_ig + 29);

    auto tr_xxxxxz_xxxx = pbuffer.data(idx_op_ig + 30);

    auto tr_xxxxxz_xxxy = pbuffer.data(idx_op_ig + 31);

    auto tr_xxxxxz_xxxz = pbuffer.data(idx_op_ig + 32);

    auto tr_xxxxxz_xxyy = pbuffer.data(idx_op_ig + 33);

    auto tr_xxxxxz_xxyz = pbuffer.data(idx_op_ig + 34);

    auto tr_xxxxxz_xxzz = pbuffer.data(idx_op_ig + 35);

    auto tr_xxxxxz_xyyy = pbuffer.data(idx_op_ig + 36);

    auto tr_xxxxxz_xyyz = pbuffer.data(idx_op_ig + 37);

    auto tr_xxxxxz_xyzz = pbuffer.data(idx_op_ig + 38);

    auto tr_xxxxxz_xzzz = pbuffer.data(idx_op_ig + 39);

    auto tr_xxxxxz_yyyy = pbuffer.data(idx_op_ig + 40);

    auto tr_xxxxxz_yyyz = pbuffer.data(idx_op_ig + 41);

    auto tr_xxxxxz_yyzz = pbuffer.data(idx_op_ig + 42);

    auto tr_xxxxxz_yzzz = pbuffer.data(idx_op_ig + 43);

    auto tr_xxxxxz_zzzz = pbuffer.data(idx_op_ig + 44);

    auto tr_xxxxyy_xxxx = pbuffer.data(idx_op_ig + 45);

    auto tr_xxxxyy_xxxy = pbuffer.data(idx_op_ig + 46);

    auto tr_xxxxyy_xxxz = pbuffer.data(idx_op_ig + 47);

    auto tr_xxxxyy_xxyy = pbuffer.data(idx_op_ig + 48);

    auto tr_xxxxyy_xxyz = pbuffer.data(idx_op_ig + 49);

    auto tr_xxxxyy_xxzz = pbuffer.data(idx_op_ig + 50);

    auto tr_xxxxyy_xyyy = pbuffer.data(idx_op_ig + 51);

    auto tr_xxxxyy_xyyz = pbuffer.data(idx_op_ig + 52);

    auto tr_xxxxyy_xyzz = pbuffer.data(idx_op_ig + 53);

    auto tr_xxxxyy_xzzz = pbuffer.data(idx_op_ig + 54);

    auto tr_xxxxyy_yyyy = pbuffer.data(idx_op_ig + 55);

    auto tr_xxxxyy_yyyz = pbuffer.data(idx_op_ig + 56);

    auto tr_xxxxyy_yyzz = pbuffer.data(idx_op_ig + 57);

    auto tr_xxxxyy_yzzz = pbuffer.data(idx_op_ig + 58);

    auto tr_xxxxyy_zzzz = pbuffer.data(idx_op_ig + 59);

    auto tr_xxxxyz_xxxx = pbuffer.data(idx_op_ig + 60);

    auto tr_xxxxyz_xxxy = pbuffer.data(idx_op_ig + 61);

    auto tr_xxxxyz_xxxz = pbuffer.data(idx_op_ig + 62);

    auto tr_xxxxyz_xxyy = pbuffer.data(idx_op_ig + 63);

    auto tr_xxxxyz_xxyz = pbuffer.data(idx_op_ig + 64);

    auto tr_xxxxyz_xxzz = pbuffer.data(idx_op_ig + 65);

    auto tr_xxxxyz_xyyy = pbuffer.data(idx_op_ig + 66);

    auto tr_xxxxyz_xyyz = pbuffer.data(idx_op_ig + 67);

    auto tr_xxxxyz_xyzz = pbuffer.data(idx_op_ig + 68);

    auto tr_xxxxyz_xzzz = pbuffer.data(idx_op_ig + 69);

    auto tr_xxxxyz_yyyy = pbuffer.data(idx_op_ig + 70);

    auto tr_xxxxyz_yyyz = pbuffer.data(idx_op_ig + 71);

    auto tr_xxxxyz_yyzz = pbuffer.data(idx_op_ig + 72);

    auto tr_xxxxyz_yzzz = pbuffer.data(idx_op_ig + 73);

    auto tr_xxxxyz_zzzz = pbuffer.data(idx_op_ig + 74);

    auto tr_xxxxzz_xxxx = pbuffer.data(idx_op_ig + 75);

    auto tr_xxxxzz_xxxy = pbuffer.data(idx_op_ig + 76);

    auto tr_xxxxzz_xxxz = pbuffer.data(idx_op_ig + 77);

    auto tr_xxxxzz_xxyy = pbuffer.data(idx_op_ig + 78);

    auto tr_xxxxzz_xxyz = pbuffer.data(idx_op_ig + 79);

    auto tr_xxxxzz_xxzz = pbuffer.data(idx_op_ig + 80);

    auto tr_xxxxzz_xyyy = pbuffer.data(idx_op_ig + 81);

    auto tr_xxxxzz_xyyz = pbuffer.data(idx_op_ig + 82);

    auto tr_xxxxzz_xyzz = pbuffer.data(idx_op_ig + 83);

    auto tr_xxxxzz_xzzz = pbuffer.data(idx_op_ig + 84);

    auto tr_xxxxzz_yyyy = pbuffer.data(idx_op_ig + 85);

    auto tr_xxxxzz_yyyz = pbuffer.data(idx_op_ig + 86);

    auto tr_xxxxzz_yyzz = pbuffer.data(idx_op_ig + 87);

    auto tr_xxxxzz_yzzz = pbuffer.data(idx_op_ig + 88);

    auto tr_xxxxzz_zzzz = pbuffer.data(idx_op_ig + 89);

    auto tr_xxxyyy_xxxx = pbuffer.data(idx_op_ig + 90);

    auto tr_xxxyyy_xxxy = pbuffer.data(idx_op_ig + 91);

    auto tr_xxxyyy_xxxz = pbuffer.data(idx_op_ig + 92);

    auto tr_xxxyyy_xxyy = pbuffer.data(idx_op_ig + 93);

    auto tr_xxxyyy_xxyz = pbuffer.data(idx_op_ig + 94);

    auto tr_xxxyyy_xxzz = pbuffer.data(idx_op_ig + 95);

    auto tr_xxxyyy_xyyy = pbuffer.data(idx_op_ig + 96);

    auto tr_xxxyyy_xyyz = pbuffer.data(idx_op_ig + 97);

    auto tr_xxxyyy_xyzz = pbuffer.data(idx_op_ig + 98);

    auto tr_xxxyyy_xzzz = pbuffer.data(idx_op_ig + 99);

    auto tr_xxxyyy_yyyy = pbuffer.data(idx_op_ig + 100);

    auto tr_xxxyyy_yyyz = pbuffer.data(idx_op_ig + 101);

    auto tr_xxxyyy_yyzz = pbuffer.data(idx_op_ig + 102);

    auto tr_xxxyyy_yzzz = pbuffer.data(idx_op_ig + 103);

    auto tr_xxxyyy_zzzz = pbuffer.data(idx_op_ig + 104);

    auto tr_xxxyyz_xxxx = pbuffer.data(idx_op_ig + 105);

    auto tr_xxxyyz_xxxy = pbuffer.data(idx_op_ig + 106);

    auto tr_xxxyyz_xxxz = pbuffer.data(idx_op_ig + 107);

    auto tr_xxxyyz_xxyy = pbuffer.data(idx_op_ig + 108);

    auto tr_xxxyyz_xxyz = pbuffer.data(idx_op_ig + 109);

    auto tr_xxxyyz_xxzz = pbuffer.data(idx_op_ig + 110);

    auto tr_xxxyyz_xyyy = pbuffer.data(idx_op_ig + 111);

    auto tr_xxxyyz_xyyz = pbuffer.data(idx_op_ig + 112);

    auto tr_xxxyyz_xyzz = pbuffer.data(idx_op_ig + 113);

    auto tr_xxxyyz_xzzz = pbuffer.data(idx_op_ig + 114);

    auto tr_xxxyyz_yyyy = pbuffer.data(idx_op_ig + 115);

    auto tr_xxxyyz_yyyz = pbuffer.data(idx_op_ig + 116);

    auto tr_xxxyyz_yyzz = pbuffer.data(idx_op_ig + 117);

    auto tr_xxxyyz_yzzz = pbuffer.data(idx_op_ig + 118);

    auto tr_xxxyyz_zzzz = pbuffer.data(idx_op_ig + 119);

    auto tr_xxxyzz_xxxx = pbuffer.data(idx_op_ig + 120);

    auto tr_xxxyzz_xxxy = pbuffer.data(idx_op_ig + 121);

    auto tr_xxxyzz_xxxz = pbuffer.data(idx_op_ig + 122);

    auto tr_xxxyzz_xxyy = pbuffer.data(idx_op_ig + 123);

    auto tr_xxxyzz_xxyz = pbuffer.data(idx_op_ig + 124);

    auto tr_xxxyzz_xxzz = pbuffer.data(idx_op_ig + 125);

    auto tr_xxxyzz_xyyy = pbuffer.data(idx_op_ig + 126);

    auto tr_xxxyzz_xyyz = pbuffer.data(idx_op_ig + 127);

    auto tr_xxxyzz_xyzz = pbuffer.data(idx_op_ig + 128);

    auto tr_xxxyzz_xzzz = pbuffer.data(idx_op_ig + 129);

    auto tr_xxxyzz_yyyy = pbuffer.data(idx_op_ig + 130);

    auto tr_xxxyzz_yyyz = pbuffer.data(idx_op_ig + 131);

    auto tr_xxxyzz_yyzz = pbuffer.data(idx_op_ig + 132);

    auto tr_xxxyzz_yzzz = pbuffer.data(idx_op_ig + 133);

    auto tr_xxxyzz_zzzz = pbuffer.data(idx_op_ig + 134);

    auto tr_xxxzzz_xxxx = pbuffer.data(idx_op_ig + 135);

    auto tr_xxxzzz_xxxy = pbuffer.data(idx_op_ig + 136);

    auto tr_xxxzzz_xxxz = pbuffer.data(idx_op_ig + 137);

    auto tr_xxxzzz_xxyy = pbuffer.data(idx_op_ig + 138);

    auto tr_xxxzzz_xxyz = pbuffer.data(idx_op_ig + 139);

    auto tr_xxxzzz_xxzz = pbuffer.data(idx_op_ig + 140);

    auto tr_xxxzzz_xyyy = pbuffer.data(idx_op_ig + 141);

    auto tr_xxxzzz_xyyz = pbuffer.data(idx_op_ig + 142);

    auto tr_xxxzzz_xyzz = pbuffer.data(idx_op_ig + 143);

    auto tr_xxxzzz_xzzz = pbuffer.data(idx_op_ig + 144);

    auto tr_xxxzzz_yyyy = pbuffer.data(idx_op_ig + 145);

    auto tr_xxxzzz_yyyz = pbuffer.data(idx_op_ig + 146);

    auto tr_xxxzzz_yyzz = pbuffer.data(idx_op_ig + 147);

    auto tr_xxxzzz_yzzz = pbuffer.data(idx_op_ig + 148);

    auto tr_xxxzzz_zzzz = pbuffer.data(idx_op_ig + 149);

    auto tr_xxyyyy_xxxx = pbuffer.data(idx_op_ig + 150);

    auto tr_xxyyyy_xxxy = pbuffer.data(idx_op_ig + 151);

    auto tr_xxyyyy_xxxz = pbuffer.data(idx_op_ig + 152);

    auto tr_xxyyyy_xxyy = pbuffer.data(idx_op_ig + 153);

    auto tr_xxyyyy_xxyz = pbuffer.data(idx_op_ig + 154);

    auto tr_xxyyyy_xxzz = pbuffer.data(idx_op_ig + 155);

    auto tr_xxyyyy_xyyy = pbuffer.data(idx_op_ig + 156);

    auto tr_xxyyyy_xyyz = pbuffer.data(idx_op_ig + 157);

    auto tr_xxyyyy_xyzz = pbuffer.data(idx_op_ig + 158);

    auto tr_xxyyyy_xzzz = pbuffer.data(idx_op_ig + 159);

    auto tr_xxyyyy_yyyy = pbuffer.data(idx_op_ig + 160);

    auto tr_xxyyyy_yyyz = pbuffer.data(idx_op_ig + 161);

    auto tr_xxyyyy_yyzz = pbuffer.data(idx_op_ig + 162);

    auto tr_xxyyyy_yzzz = pbuffer.data(idx_op_ig + 163);

    auto tr_xxyyyy_zzzz = pbuffer.data(idx_op_ig + 164);

    auto tr_xxyyyz_xxxx = pbuffer.data(idx_op_ig + 165);

    auto tr_xxyyyz_xxxy = pbuffer.data(idx_op_ig + 166);

    auto tr_xxyyyz_xxxz = pbuffer.data(idx_op_ig + 167);

    auto tr_xxyyyz_xxyy = pbuffer.data(idx_op_ig + 168);

    auto tr_xxyyyz_xxyz = pbuffer.data(idx_op_ig + 169);

    auto tr_xxyyyz_xxzz = pbuffer.data(idx_op_ig + 170);

    auto tr_xxyyyz_xyyy = pbuffer.data(idx_op_ig + 171);

    auto tr_xxyyyz_xyyz = pbuffer.data(idx_op_ig + 172);

    auto tr_xxyyyz_xyzz = pbuffer.data(idx_op_ig + 173);

    auto tr_xxyyyz_xzzz = pbuffer.data(idx_op_ig + 174);

    auto tr_xxyyyz_yyyy = pbuffer.data(idx_op_ig + 175);

    auto tr_xxyyyz_yyyz = pbuffer.data(idx_op_ig + 176);

    auto tr_xxyyyz_yyzz = pbuffer.data(idx_op_ig + 177);

    auto tr_xxyyyz_yzzz = pbuffer.data(idx_op_ig + 178);

    auto tr_xxyyyz_zzzz = pbuffer.data(idx_op_ig + 179);

    auto tr_xxyyzz_xxxx = pbuffer.data(idx_op_ig + 180);

    auto tr_xxyyzz_xxxy = pbuffer.data(idx_op_ig + 181);

    auto tr_xxyyzz_xxxz = pbuffer.data(idx_op_ig + 182);

    auto tr_xxyyzz_xxyy = pbuffer.data(idx_op_ig + 183);

    auto tr_xxyyzz_xxyz = pbuffer.data(idx_op_ig + 184);

    auto tr_xxyyzz_xxzz = pbuffer.data(idx_op_ig + 185);

    auto tr_xxyyzz_xyyy = pbuffer.data(idx_op_ig + 186);

    auto tr_xxyyzz_xyyz = pbuffer.data(idx_op_ig + 187);

    auto tr_xxyyzz_xyzz = pbuffer.data(idx_op_ig + 188);

    auto tr_xxyyzz_xzzz = pbuffer.data(idx_op_ig + 189);

    auto tr_xxyyzz_yyyy = pbuffer.data(idx_op_ig + 190);

    auto tr_xxyyzz_yyyz = pbuffer.data(idx_op_ig + 191);

    auto tr_xxyyzz_yyzz = pbuffer.data(idx_op_ig + 192);

    auto tr_xxyyzz_yzzz = pbuffer.data(idx_op_ig + 193);

    auto tr_xxyyzz_zzzz = pbuffer.data(idx_op_ig + 194);

    auto tr_xxyzzz_xxxx = pbuffer.data(idx_op_ig + 195);

    auto tr_xxyzzz_xxxy = pbuffer.data(idx_op_ig + 196);

    auto tr_xxyzzz_xxxz = pbuffer.data(idx_op_ig + 197);

    auto tr_xxyzzz_xxyy = pbuffer.data(idx_op_ig + 198);

    auto tr_xxyzzz_xxyz = pbuffer.data(idx_op_ig + 199);

    auto tr_xxyzzz_xxzz = pbuffer.data(idx_op_ig + 200);

    auto tr_xxyzzz_xyyy = pbuffer.data(idx_op_ig + 201);

    auto tr_xxyzzz_xyyz = pbuffer.data(idx_op_ig + 202);

    auto tr_xxyzzz_xyzz = pbuffer.data(idx_op_ig + 203);

    auto tr_xxyzzz_xzzz = pbuffer.data(idx_op_ig + 204);

    auto tr_xxyzzz_yyyy = pbuffer.data(idx_op_ig + 205);

    auto tr_xxyzzz_yyyz = pbuffer.data(idx_op_ig + 206);

    auto tr_xxyzzz_yyzz = pbuffer.data(idx_op_ig + 207);

    auto tr_xxyzzz_yzzz = pbuffer.data(idx_op_ig + 208);

    auto tr_xxyzzz_zzzz = pbuffer.data(idx_op_ig + 209);

    auto tr_xxzzzz_xxxx = pbuffer.data(idx_op_ig + 210);

    auto tr_xxzzzz_xxxy = pbuffer.data(idx_op_ig + 211);

    auto tr_xxzzzz_xxxz = pbuffer.data(idx_op_ig + 212);

    auto tr_xxzzzz_xxyy = pbuffer.data(idx_op_ig + 213);

    auto tr_xxzzzz_xxyz = pbuffer.data(idx_op_ig + 214);

    auto tr_xxzzzz_xxzz = pbuffer.data(idx_op_ig + 215);

    auto tr_xxzzzz_xyyy = pbuffer.data(idx_op_ig + 216);

    auto tr_xxzzzz_xyyz = pbuffer.data(idx_op_ig + 217);

    auto tr_xxzzzz_xyzz = pbuffer.data(idx_op_ig + 218);

    auto tr_xxzzzz_xzzz = pbuffer.data(idx_op_ig + 219);

    auto tr_xxzzzz_yyyy = pbuffer.data(idx_op_ig + 220);

    auto tr_xxzzzz_yyyz = pbuffer.data(idx_op_ig + 221);

    auto tr_xxzzzz_yyzz = pbuffer.data(idx_op_ig + 222);

    auto tr_xxzzzz_yzzz = pbuffer.data(idx_op_ig + 223);

    auto tr_xxzzzz_zzzz = pbuffer.data(idx_op_ig + 224);

    auto tr_xyyyyy_xxxx = pbuffer.data(idx_op_ig + 225);

    auto tr_xyyyyy_xxxy = pbuffer.data(idx_op_ig + 226);

    auto tr_xyyyyy_xxxz = pbuffer.data(idx_op_ig + 227);

    auto tr_xyyyyy_xxyy = pbuffer.data(idx_op_ig + 228);

    auto tr_xyyyyy_xxyz = pbuffer.data(idx_op_ig + 229);

    auto tr_xyyyyy_xxzz = pbuffer.data(idx_op_ig + 230);

    auto tr_xyyyyy_xyyy = pbuffer.data(idx_op_ig + 231);

    auto tr_xyyyyy_xyyz = pbuffer.data(idx_op_ig + 232);

    auto tr_xyyyyy_xyzz = pbuffer.data(idx_op_ig + 233);

    auto tr_xyyyyy_xzzz = pbuffer.data(idx_op_ig + 234);

    auto tr_xyyyyy_yyyy = pbuffer.data(idx_op_ig + 235);

    auto tr_xyyyyy_yyyz = pbuffer.data(idx_op_ig + 236);

    auto tr_xyyyyy_yyzz = pbuffer.data(idx_op_ig + 237);

    auto tr_xyyyyy_yzzz = pbuffer.data(idx_op_ig + 238);

    auto tr_xyyyyy_zzzz = pbuffer.data(idx_op_ig + 239);

    auto tr_xyyyyz_xxxx = pbuffer.data(idx_op_ig + 240);

    auto tr_xyyyyz_xxxy = pbuffer.data(idx_op_ig + 241);

    auto tr_xyyyyz_xxxz = pbuffer.data(idx_op_ig + 242);

    auto tr_xyyyyz_xxyy = pbuffer.data(idx_op_ig + 243);

    auto tr_xyyyyz_xxyz = pbuffer.data(idx_op_ig + 244);

    auto tr_xyyyyz_xxzz = pbuffer.data(idx_op_ig + 245);

    auto tr_xyyyyz_xyyy = pbuffer.data(idx_op_ig + 246);

    auto tr_xyyyyz_xyyz = pbuffer.data(idx_op_ig + 247);

    auto tr_xyyyyz_xyzz = pbuffer.data(idx_op_ig + 248);

    auto tr_xyyyyz_xzzz = pbuffer.data(idx_op_ig + 249);

    auto tr_xyyyyz_yyyy = pbuffer.data(idx_op_ig + 250);

    auto tr_xyyyyz_yyyz = pbuffer.data(idx_op_ig + 251);

    auto tr_xyyyyz_yyzz = pbuffer.data(idx_op_ig + 252);

    auto tr_xyyyyz_yzzz = pbuffer.data(idx_op_ig + 253);

    auto tr_xyyyyz_zzzz = pbuffer.data(idx_op_ig + 254);

    auto tr_xyyyzz_xxxx = pbuffer.data(idx_op_ig + 255);

    auto tr_xyyyzz_xxxy = pbuffer.data(idx_op_ig + 256);

    auto tr_xyyyzz_xxxz = pbuffer.data(idx_op_ig + 257);

    auto tr_xyyyzz_xxyy = pbuffer.data(idx_op_ig + 258);

    auto tr_xyyyzz_xxyz = pbuffer.data(idx_op_ig + 259);

    auto tr_xyyyzz_xxzz = pbuffer.data(idx_op_ig + 260);

    auto tr_xyyyzz_xyyy = pbuffer.data(idx_op_ig + 261);

    auto tr_xyyyzz_xyyz = pbuffer.data(idx_op_ig + 262);

    auto tr_xyyyzz_xyzz = pbuffer.data(idx_op_ig + 263);

    auto tr_xyyyzz_xzzz = pbuffer.data(idx_op_ig + 264);

    auto tr_xyyyzz_yyyy = pbuffer.data(idx_op_ig + 265);

    auto tr_xyyyzz_yyyz = pbuffer.data(idx_op_ig + 266);

    auto tr_xyyyzz_yyzz = pbuffer.data(idx_op_ig + 267);

    auto tr_xyyyzz_yzzz = pbuffer.data(idx_op_ig + 268);

    auto tr_xyyyzz_zzzz = pbuffer.data(idx_op_ig + 269);

    auto tr_xyyzzz_xxxx = pbuffer.data(idx_op_ig + 270);

    auto tr_xyyzzz_xxxy = pbuffer.data(idx_op_ig + 271);

    auto tr_xyyzzz_xxxz = pbuffer.data(idx_op_ig + 272);

    auto tr_xyyzzz_xxyy = pbuffer.data(idx_op_ig + 273);

    auto tr_xyyzzz_xxyz = pbuffer.data(idx_op_ig + 274);

    auto tr_xyyzzz_xxzz = pbuffer.data(idx_op_ig + 275);

    auto tr_xyyzzz_xyyy = pbuffer.data(idx_op_ig + 276);

    auto tr_xyyzzz_xyyz = pbuffer.data(idx_op_ig + 277);

    auto tr_xyyzzz_xyzz = pbuffer.data(idx_op_ig + 278);

    auto tr_xyyzzz_xzzz = pbuffer.data(idx_op_ig + 279);

    auto tr_xyyzzz_yyyy = pbuffer.data(idx_op_ig + 280);

    auto tr_xyyzzz_yyyz = pbuffer.data(idx_op_ig + 281);

    auto tr_xyyzzz_yyzz = pbuffer.data(idx_op_ig + 282);

    auto tr_xyyzzz_yzzz = pbuffer.data(idx_op_ig + 283);

    auto tr_xyyzzz_zzzz = pbuffer.data(idx_op_ig + 284);

    auto tr_xyzzzz_xxxx = pbuffer.data(idx_op_ig + 285);

    auto tr_xyzzzz_xxxy = pbuffer.data(idx_op_ig + 286);

    auto tr_xyzzzz_xxxz = pbuffer.data(idx_op_ig + 287);

    auto tr_xyzzzz_xxyy = pbuffer.data(idx_op_ig + 288);

    auto tr_xyzzzz_xxyz = pbuffer.data(idx_op_ig + 289);

    auto tr_xyzzzz_xxzz = pbuffer.data(idx_op_ig + 290);

    auto tr_xyzzzz_xyyy = pbuffer.data(idx_op_ig + 291);

    auto tr_xyzzzz_xyyz = pbuffer.data(idx_op_ig + 292);

    auto tr_xyzzzz_xyzz = pbuffer.data(idx_op_ig + 293);

    auto tr_xyzzzz_xzzz = pbuffer.data(idx_op_ig + 294);

    auto tr_xyzzzz_yyyy = pbuffer.data(idx_op_ig + 295);

    auto tr_xyzzzz_yyyz = pbuffer.data(idx_op_ig + 296);

    auto tr_xyzzzz_yyzz = pbuffer.data(idx_op_ig + 297);

    auto tr_xyzzzz_yzzz = pbuffer.data(idx_op_ig + 298);

    auto tr_xyzzzz_zzzz = pbuffer.data(idx_op_ig + 299);

    auto tr_xzzzzz_xxxx = pbuffer.data(idx_op_ig + 300);

    auto tr_xzzzzz_xxxy = pbuffer.data(idx_op_ig + 301);

    auto tr_xzzzzz_xxxz = pbuffer.data(idx_op_ig + 302);

    auto tr_xzzzzz_xxyy = pbuffer.data(idx_op_ig + 303);

    auto tr_xzzzzz_xxyz = pbuffer.data(idx_op_ig + 304);

    auto tr_xzzzzz_xxzz = pbuffer.data(idx_op_ig + 305);

    auto tr_xzzzzz_xyyy = pbuffer.data(idx_op_ig + 306);

    auto tr_xzzzzz_xyyz = pbuffer.data(idx_op_ig + 307);

    auto tr_xzzzzz_xyzz = pbuffer.data(idx_op_ig + 308);

    auto tr_xzzzzz_xzzz = pbuffer.data(idx_op_ig + 309);

    auto tr_xzzzzz_yyyy = pbuffer.data(idx_op_ig + 310);

    auto tr_xzzzzz_yyyz = pbuffer.data(idx_op_ig + 311);

    auto tr_xzzzzz_yyzz = pbuffer.data(idx_op_ig + 312);

    auto tr_xzzzzz_yzzz = pbuffer.data(idx_op_ig + 313);

    auto tr_xzzzzz_zzzz = pbuffer.data(idx_op_ig + 314);

    auto tr_yyyyyy_xxxx = pbuffer.data(idx_op_ig + 315);

    auto tr_yyyyyy_xxxy = pbuffer.data(idx_op_ig + 316);

    auto tr_yyyyyy_xxxz = pbuffer.data(idx_op_ig + 317);

    auto tr_yyyyyy_xxyy = pbuffer.data(idx_op_ig + 318);

    auto tr_yyyyyy_xxyz = pbuffer.data(idx_op_ig + 319);

    auto tr_yyyyyy_xxzz = pbuffer.data(idx_op_ig + 320);

    auto tr_yyyyyy_xyyy = pbuffer.data(idx_op_ig + 321);

    auto tr_yyyyyy_xyyz = pbuffer.data(idx_op_ig + 322);

    auto tr_yyyyyy_xyzz = pbuffer.data(idx_op_ig + 323);

    auto tr_yyyyyy_xzzz = pbuffer.data(idx_op_ig + 324);

    auto tr_yyyyyy_yyyy = pbuffer.data(idx_op_ig + 325);

    auto tr_yyyyyy_yyyz = pbuffer.data(idx_op_ig + 326);

    auto tr_yyyyyy_yyzz = pbuffer.data(idx_op_ig + 327);

    auto tr_yyyyyy_yzzz = pbuffer.data(idx_op_ig + 328);

    auto tr_yyyyyy_zzzz = pbuffer.data(idx_op_ig + 329);

    auto tr_yyyyyz_xxxx = pbuffer.data(idx_op_ig + 330);

    auto tr_yyyyyz_xxxy = pbuffer.data(idx_op_ig + 331);

    auto tr_yyyyyz_xxxz = pbuffer.data(idx_op_ig + 332);

    auto tr_yyyyyz_xxyy = pbuffer.data(idx_op_ig + 333);

    auto tr_yyyyyz_xxyz = pbuffer.data(idx_op_ig + 334);

    auto tr_yyyyyz_xxzz = pbuffer.data(idx_op_ig + 335);

    auto tr_yyyyyz_xyyy = pbuffer.data(idx_op_ig + 336);

    auto tr_yyyyyz_xyyz = pbuffer.data(idx_op_ig + 337);

    auto tr_yyyyyz_xyzz = pbuffer.data(idx_op_ig + 338);

    auto tr_yyyyyz_xzzz = pbuffer.data(idx_op_ig + 339);

    auto tr_yyyyyz_yyyy = pbuffer.data(idx_op_ig + 340);

    auto tr_yyyyyz_yyyz = pbuffer.data(idx_op_ig + 341);

    auto tr_yyyyyz_yyzz = pbuffer.data(idx_op_ig + 342);

    auto tr_yyyyyz_yzzz = pbuffer.data(idx_op_ig + 343);

    auto tr_yyyyyz_zzzz = pbuffer.data(idx_op_ig + 344);

    auto tr_yyyyzz_xxxx = pbuffer.data(idx_op_ig + 345);

    auto tr_yyyyzz_xxxy = pbuffer.data(idx_op_ig + 346);

    auto tr_yyyyzz_xxxz = pbuffer.data(idx_op_ig + 347);

    auto tr_yyyyzz_xxyy = pbuffer.data(idx_op_ig + 348);

    auto tr_yyyyzz_xxyz = pbuffer.data(idx_op_ig + 349);

    auto tr_yyyyzz_xxzz = pbuffer.data(idx_op_ig + 350);

    auto tr_yyyyzz_xyyy = pbuffer.data(idx_op_ig + 351);

    auto tr_yyyyzz_xyyz = pbuffer.data(idx_op_ig + 352);

    auto tr_yyyyzz_xyzz = pbuffer.data(idx_op_ig + 353);

    auto tr_yyyyzz_xzzz = pbuffer.data(idx_op_ig + 354);

    auto tr_yyyyzz_yyyy = pbuffer.data(idx_op_ig + 355);

    auto tr_yyyyzz_yyyz = pbuffer.data(idx_op_ig + 356);

    auto tr_yyyyzz_yyzz = pbuffer.data(idx_op_ig + 357);

    auto tr_yyyyzz_yzzz = pbuffer.data(idx_op_ig + 358);

    auto tr_yyyyzz_zzzz = pbuffer.data(idx_op_ig + 359);

    auto tr_yyyzzz_xxxx = pbuffer.data(idx_op_ig + 360);

    auto tr_yyyzzz_xxxy = pbuffer.data(idx_op_ig + 361);

    auto tr_yyyzzz_xxxz = pbuffer.data(idx_op_ig + 362);

    auto tr_yyyzzz_xxyy = pbuffer.data(idx_op_ig + 363);

    auto tr_yyyzzz_xxyz = pbuffer.data(idx_op_ig + 364);

    auto tr_yyyzzz_xxzz = pbuffer.data(idx_op_ig + 365);

    auto tr_yyyzzz_xyyy = pbuffer.data(idx_op_ig + 366);

    auto tr_yyyzzz_xyyz = pbuffer.data(idx_op_ig + 367);

    auto tr_yyyzzz_xyzz = pbuffer.data(idx_op_ig + 368);

    auto tr_yyyzzz_xzzz = pbuffer.data(idx_op_ig + 369);

    auto tr_yyyzzz_yyyy = pbuffer.data(idx_op_ig + 370);

    auto tr_yyyzzz_yyyz = pbuffer.data(idx_op_ig + 371);

    auto tr_yyyzzz_yyzz = pbuffer.data(idx_op_ig + 372);

    auto tr_yyyzzz_yzzz = pbuffer.data(idx_op_ig + 373);

    auto tr_yyyzzz_zzzz = pbuffer.data(idx_op_ig + 374);

    auto tr_yyzzzz_xxxx = pbuffer.data(idx_op_ig + 375);

    auto tr_yyzzzz_xxxy = pbuffer.data(idx_op_ig + 376);

    auto tr_yyzzzz_xxxz = pbuffer.data(idx_op_ig + 377);

    auto tr_yyzzzz_xxyy = pbuffer.data(idx_op_ig + 378);

    auto tr_yyzzzz_xxyz = pbuffer.data(idx_op_ig + 379);

    auto tr_yyzzzz_xxzz = pbuffer.data(idx_op_ig + 380);

    auto tr_yyzzzz_xyyy = pbuffer.data(idx_op_ig + 381);

    auto tr_yyzzzz_xyyz = pbuffer.data(idx_op_ig + 382);

    auto tr_yyzzzz_xyzz = pbuffer.data(idx_op_ig + 383);

    auto tr_yyzzzz_xzzz = pbuffer.data(idx_op_ig + 384);

    auto tr_yyzzzz_yyyy = pbuffer.data(idx_op_ig + 385);

    auto tr_yyzzzz_yyyz = pbuffer.data(idx_op_ig + 386);

    auto tr_yyzzzz_yyzz = pbuffer.data(idx_op_ig + 387);

    auto tr_yyzzzz_yzzz = pbuffer.data(idx_op_ig + 388);

    auto tr_yyzzzz_zzzz = pbuffer.data(idx_op_ig + 389);

    auto tr_yzzzzz_xxxx = pbuffer.data(idx_op_ig + 390);

    auto tr_yzzzzz_xxxy = pbuffer.data(idx_op_ig + 391);

    auto tr_yzzzzz_xxxz = pbuffer.data(idx_op_ig + 392);

    auto tr_yzzzzz_xxyy = pbuffer.data(idx_op_ig + 393);

    auto tr_yzzzzz_xxyz = pbuffer.data(idx_op_ig + 394);

    auto tr_yzzzzz_xxzz = pbuffer.data(idx_op_ig + 395);

    auto tr_yzzzzz_xyyy = pbuffer.data(idx_op_ig + 396);

    auto tr_yzzzzz_xyyz = pbuffer.data(idx_op_ig + 397);

    auto tr_yzzzzz_xyzz = pbuffer.data(idx_op_ig + 398);

    auto tr_yzzzzz_xzzz = pbuffer.data(idx_op_ig + 399);

    auto tr_yzzzzz_yyyy = pbuffer.data(idx_op_ig + 400);

    auto tr_yzzzzz_yyyz = pbuffer.data(idx_op_ig + 401);

    auto tr_yzzzzz_yyzz = pbuffer.data(idx_op_ig + 402);

    auto tr_yzzzzz_yzzz = pbuffer.data(idx_op_ig + 403);

    auto tr_yzzzzz_zzzz = pbuffer.data(idx_op_ig + 404);

    auto tr_zzzzzz_xxxx = pbuffer.data(idx_op_ig + 405);

    auto tr_zzzzzz_xxxy = pbuffer.data(idx_op_ig + 406);

    auto tr_zzzzzz_xxxz = pbuffer.data(idx_op_ig + 407);

    auto tr_zzzzzz_xxyy = pbuffer.data(idx_op_ig + 408);

    auto tr_zzzzzz_xxyz = pbuffer.data(idx_op_ig + 409);

    auto tr_zzzzzz_xxzz = pbuffer.data(idx_op_ig + 410);

    auto tr_zzzzzz_xyyy = pbuffer.data(idx_op_ig + 411);

    auto tr_zzzzzz_xyyz = pbuffer.data(idx_op_ig + 412);

    auto tr_zzzzzz_xyzz = pbuffer.data(idx_op_ig + 413);

    auto tr_zzzzzz_xzzz = pbuffer.data(idx_op_ig + 414);

    auto tr_zzzzzz_yyyy = pbuffer.data(idx_op_ig + 415);

    auto tr_zzzzzz_yyyz = pbuffer.data(idx_op_ig + 416);

    auto tr_zzzzzz_yyzz = pbuffer.data(idx_op_ig + 417);

    auto tr_zzzzzz_yzzz = pbuffer.data(idx_op_ig + 418);

    auto tr_zzzzzz_zzzz = pbuffer.data(idx_op_ig + 419);

    // Set up 0-15 components of targeted buffer : GG

    auto tr_0_0_xx_xxxx_xxxx = pbuffer.data(idx_op_geom_020_gg);

    auto tr_0_0_xx_xxxx_xxxy = pbuffer.data(idx_op_geom_020_gg + 1);

    auto tr_0_0_xx_xxxx_xxxz = pbuffer.data(idx_op_geom_020_gg + 2);

    auto tr_0_0_xx_xxxx_xxyy = pbuffer.data(idx_op_geom_020_gg + 3);

    auto tr_0_0_xx_xxxx_xxyz = pbuffer.data(idx_op_geom_020_gg + 4);

    auto tr_0_0_xx_xxxx_xxzz = pbuffer.data(idx_op_geom_020_gg + 5);

    auto tr_0_0_xx_xxxx_xyyy = pbuffer.data(idx_op_geom_020_gg + 6);

    auto tr_0_0_xx_xxxx_xyyz = pbuffer.data(idx_op_geom_020_gg + 7);

    auto tr_0_0_xx_xxxx_xyzz = pbuffer.data(idx_op_geom_020_gg + 8);

    auto tr_0_0_xx_xxxx_xzzz = pbuffer.data(idx_op_geom_020_gg + 9);

    auto tr_0_0_xx_xxxx_yyyy = pbuffer.data(idx_op_geom_020_gg + 10);

    auto tr_0_0_xx_xxxx_yyyz = pbuffer.data(idx_op_geom_020_gg + 11);

    auto tr_0_0_xx_xxxx_yyzz = pbuffer.data(idx_op_geom_020_gg + 12);

    auto tr_0_0_xx_xxxx_yzzz = pbuffer.data(idx_op_geom_020_gg + 13);

    auto tr_0_0_xx_xxxx_zzzz = pbuffer.data(idx_op_geom_020_gg + 14);

    #pragma omp simd aligned(tr_0_0_xx_xxxx_xxxx, tr_0_0_xx_xxxx_xxxy, tr_0_0_xx_xxxx_xxxz, tr_0_0_xx_xxxx_xxyy, tr_0_0_xx_xxxx_xxyz, tr_0_0_xx_xxxx_xxzz, tr_0_0_xx_xxxx_xyyy, tr_0_0_xx_xxxx_xyyz, tr_0_0_xx_xxxx_xyzz, tr_0_0_xx_xxxx_xzzz, tr_0_0_xx_xxxx_yyyy, tr_0_0_xx_xxxx_yyyz, tr_0_0_xx_xxxx_yyzz, tr_0_0_xx_xxxx_yzzz, tr_0_0_xx_xxxx_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxx_xxx, tr_xxx_xxxxx, tr_xxx_xxxxy, tr_xxx_xxxxz, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxxxx, tr_xxxx_xxxxxy, tr_xxxx_xxxxxz, tr_xxxx_xxxxyy, tr_xxxx_xxxxyz, tr_xxxx_xxxxzz, tr_xxxx_xxxy, tr_xxxx_xxxyyy, tr_xxxx_xxxyyz, tr_xxxx_xxxyzz, tr_xxxx_xxxz, tr_xxxx_xxxzzz, tr_xxxx_xxyy, tr_xxxx_xxyyyy, tr_xxxx_xxyyyz, tr_xxxx_xxyyzz, tr_xxxx_xxyz, tr_xxxx_xxyzzz, tr_xxxx_xxzz, tr_xxxx_xxzzzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxx_xxx, tr_xxxxx_xxxxx, tr_xxxxx_xxxxy, tr_xxxxx_xxxxz, tr_xxxxx_xxxyy, tr_xxxxx_xxxyz, tr_xxxxx_xxxzz, tr_xxxxx_xxy, tr_xxxxx_xxyyy, tr_xxxxx_xxyyz, tr_xxxxx_xxyzz, tr_xxxxx_xxz, tr_xxxxx_xxzzz, tr_xxxxx_xyy, tr_xxxxx_xyyyy, tr_xxxxx_xyyyz, tr_xxxxx_xyyzz, tr_xxxxx_xyz, tr_xxxxx_xyzzz, tr_xxxxx_xzz, tr_xxxxx_xzzzz, tr_xxxxx_yyy, tr_xxxxx_yyz, tr_xxxxx_yzz, tr_xxxxx_zzz, tr_xxxxxx_xxxx, tr_xxxxxx_xxxy, tr_xxxxxx_xxxz, tr_xxxxxx_xxyy, tr_xxxxxx_xxyz, tr_xxxxxx_xxzz, tr_xxxxxx_xyyy, tr_xxxxxx_xyyz, tr_xxxxxx_xyzz, tr_xxxxxx_xzzz, tr_xxxxxx_yyyy, tr_xxxxxx_yyyz, tr_xxxxxx_yyzz, tr_xxxxxx_yzzz, tr_xxxxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxx_xxxx[i] = 12.0 * tr_xx_xxxx[i] + 32.0 * tr_xxx_xxx[i] - 16.0 * tr_xxx_xxxxx[i] * tke_0 + 12.0 * tr_xxxx_xx[i] - 18.0 * tr_xxxx_xxxx[i] * tbe_0 - 18.0 * tr_xxxx_xxxx[i] * tke_0 + 4.0 * tr_xxxx_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxxxx_xxx[i] * tbe_0 + 8.0 * tr_xxxxx_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xxxy[i] = 12.0 * tr_xx_xxxy[i] + 24.0 * tr_xxx_xxy[i] - 16.0 * tr_xxx_xxxxy[i] * tke_0 + 6.0 * tr_xxxx_xy[i] - 18.0 * tr_xxxx_xxxy[i] * tbe_0 - 14.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxxxx_xxy[i] * tbe_0 + 8.0 * tr_xxxxx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xxxz[i] = 12.0 * tr_xx_xxxz[i] + 24.0 * tr_xxx_xxz[i] - 16.0 * tr_xxx_xxxxz[i] * tke_0 + 6.0 * tr_xxxx_xz[i] - 18.0 * tr_xxxx_xxxz[i] * tbe_0 - 14.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxxxx_xxz[i] * tbe_0 + 8.0 * tr_xxxxx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xxyy[i] = 12.0 * tr_xx_xxyy[i] + 16.0 * tr_xxx_xyy[i] - 16.0 * tr_xxx_xxxyy[i] * tke_0 + 2.0 * tr_xxxx_yy[i] - 18.0 * tr_xxxx_xxyy[i] * tbe_0 - 10.0 * tr_xxxx_xxyy[i] * tke_0 + 4.0 * tr_xxxx_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxxx_xyy[i] * tbe_0 + 8.0 * tr_xxxxx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xxyz[i] = 12.0 * tr_xx_xxyz[i] + 16.0 * tr_xxx_xyz[i] - 16.0 * tr_xxx_xxxyz[i] * tke_0 + 2.0 * tr_xxxx_yz[i] - 18.0 * tr_xxxx_xxyz[i] * tbe_0 - 10.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxx_xyz[i] * tbe_0 + 8.0 * tr_xxxxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xxzz[i] = 12.0 * tr_xx_xxzz[i] + 16.0 * tr_xxx_xzz[i] - 16.0 * tr_xxx_xxxzz[i] * tke_0 + 2.0 * tr_xxxx_zz[i] - 18.0 * tr_xxxx_xxzz[i] * tbe_0 - 10.0 * tr_xxxx_xxzz[i] * tke_0 + 4.0 * tr_xxxx_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxx_xzz[i] * tbe_0 + 8.0 * tr_xxxxx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xyyy[i] = 12.0 * tr_xx_xyyy[i] + 8.0 * tr_xxx_yyy[i] - 16.0 * tr_xxx_xxyyy[i] * tke_0 - 18.0 * tr_xxxx_xyyy[i] * tbe_0 - 6.0 * tr_xxxx_xyyy[i] * tke_0 + 4.0 * tr_xxxx_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxx_yyy[i] * tbe_0 + 8.0 * tr_xxxxx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xyyz[i] = 12.0 * tr_xx_xyyz[i] + 8.0 * tr_xxx_yyz[i] - 16.0 * tr_xxx_xxyyz[i] * tke_0 - 18.0 * tr_xxxx_xyyz[i] * tbe_0 - 6.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxx_yyz[i] * tbe_0 + 8.0 * tr_xxxxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xyzz[i] = 12.0 * tr_xx_xyzz[i] + 8.0 * tr_xxx_yzz[i] - 16.0 * tr_xxx_xxyzz[i] * tke_0 - 18.0 * tr_xxxx_xyzz[i] * tbe_0 - 6.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxx_yzz[i] * tbe_0 + 8.0 * tr_xxxxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xzzz[i] = 12.0 * tr_xx_xzzz[i] + 8.0 * tr_xxx_zzz[i] - 16.0 * tr_xxx_xxzzz[i] * tke_0 - 18.0 * tr_xxxx_xzzz[i] * tbe_0 - 6.0 * tr_xxxx_xzzz[i] * tke_0 + 4.0 * tr_xxxx_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxx_zzz[i] * tbe_0 + 8.0 * tr_xxxxx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_yyyy[i] = 12.0 * tr_xx_yyyy[i] - 16.0 * tr_xxx_xyyyy[i] * tke_0 - 18.0 * tr_xxxx_yyyy[i] * tbe_0 - 2.0 * tr_xxxx_yyyy[i] * tke_0 + 4.0 * tr_xxxx_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_yyyz[i] = 12.0 * tr_xx_yyyz[i] - 16.0 * tr_xxx_xyyyz[i] * tke_0 - 18.0 * tr_xxxx_yyyz[i] * tbe_0 - 2.0 * tr_xxxx_yyyz[i] * tke_0 + 4.0 * tr_xxxx_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_yyzz[i] = 12.0 * tr_xx_yyzz[i] - 16.0 * tr_xxx_xyyzz[i] * tke_0 - 18.0 * tr_xxxx_yyzz[i] * tbe_0 - 2.0 * tr_xxxx_yyzz[i] * tke_0 + 4.0 * tr_xxxx_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_yzzz[i] = 12.0 * tr_xx_yzzz[i] - 16.0 * tr_xxx_xyzzz[i] * tke_0 - 18.0 * tr_xxxx_yzzz[i] * tbe_0 - 2.0 * tr_xxxx_yzzz[i] * tke_0 + 4.0 * tr_xxxx_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_zzzz[i] = 12.0 * tr_xx_zzzz[i] - 16.0 * tr_xxx_xzzzz[i] * tke_0 - 18.0 * tr_xxxx_zzzz[i] * tbe_0 - 2.0 * tr_xxxx_zzzz[i] * tke_0 + 4.0 * tr_xxxx_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 15-30 components of targeted buffer : GG

    auto tr_0_0_xx_xxxy_xxxx = pbuffer.data(idx_op_geom_020_gg + 15);

    auto tr_0_0_xx_xxxy_xxxy = pbuffer.data(idx_op_geom_020_gg + 16);

    auto tr_0_0_xx_xxxy_xxxz = pbuffer.data(idx_op_geom_020_gg + 17);

    auto tr_0_0_xx_xxxy_xxyy = pbuffer.data(idx_op_geom_020_gg + 18);

    auto tr_0_0_xx_xxxy_xxyz = pbuffer.data(idx_op_geom_020_gg + 19);

    auto tr_0_0_xx_xxxy_xxzz = pbuffer.data(idx_op_geom_020_gg + 20);

    auto tr_0_0_xx_xxxy_xyyy = pbuffer.data(idx_op_geom_020_gg + 21);

    auto tr_0_0_xx_xxxy_xyyz = pbuffer.data(idx_op_geom_020_gg + 22);

    auto tr_0_0_xx_xxxy_xyzz = pbuffer.data(idx_op_geom_020_gg + 23);

    auto tr_0_0_xx_xxxy_xzzz = pbuffer.data(idx_op_geom_020_gg + 24);

    auto tr_0_0_xx_xxxy_yyyy = pbuffer.data(idx_op_geom_020_gg + 25);

    auto tr_0_0_xx_xxxy_yyyz = pbuffer.data(idx_op_geom_020_gg + 26);

    auto tr_0_0_xx_xxxy_yyzz = pbuffer.data(idx_op_geom_020_gg + 27);

    auto tr_0_0_xx_xxxy_yzzz = pbuffer.data(idx_op_geom_020_gg + 28);

    auto tr_0_0_xx_xxxy_zzzz = pbuffer.data(idx_op_geom_020_gg + 29);

    #pragma omp simd aligned(tr_0_0_xx_xxxy_xxxx, tr_0_0_xx_xxxy_xxxy, tr_0_0_xx_xxxy_xxxz, tr_0_0_xx_xxxy_xxyy, tr_0_0_xx_xxxy_xxyz, tr_0_0_xx_xxxy_xxzz, tr_0_0_xx_xxxy_xyyy, tr_0_0_xx_xxxy_xyyz, tr_0_0_xx_xxxy_xyzz, tr_0_0_xx_xxxy_xzzz, tr_0_0_xx_xxxy_yyyy, tr_0_0_xx_xxxy_yyyz, tr_0_0_xx_xxxy_yyzz, tr_0_0_xx_xxxy_yzzz, tr_0_0_xx_xxxy_zzzz, tr_xxxxxy_xxxx, tr_xxxxxy_xxxy, tr_xxxxxy_xxxz, tr_xxxxxy_xxyy, tr_xxxxxy_xxyz, tr_xxxxxy_xxzz, tr_xxxxxy_xyyy, tr_xxxxxy_xyyz, tr_xxxxxy_xyzz, tr_xxxxxy_xzzz, tr_xxxxxy_yyyy, tr_xxxxxy_yyyz, tr_xxxxxy_yyzz, tr_xxxxxy_yzzz, tr_xxxxxy_zzzz, tr_xxxxy_xxx, tr_xxxxy_xxxxx, tr_xxxxy_xxxxy, tr_xxxxy_xxxxz, tr_xxxxy_xxxyy, tr_xxxxy_xxxyz, tr_xxxxy_xxxzz, tr_xxxxy_xxy, tr_xxxxy_xxyyy, tr_xxxxy_xxyyz, tr_xxxxy_xxyzz, tr_xxxxy_xxz, tr_xxxxy_xxzzz, tr_xxxxy_xyy, tr_xxxxy_xyyyy, tr_xxxxy_xyyyz, tr_xxxxy_xyyzz, tr_xxxxy_xyz, tr_xxxxy_xyzzz, tr_xxxxy_xzz, tr_xxxxy_xzzzz, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_zzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxxxx, tr_xxxy_xxxxxy, tr_xxxy_xxxxxz, tr_xxxy_xxxxyy, tr_xxxy_xxxxyz, tr_xxxy_xxxxzz, tr_xxxy_xxxy, tr_xxxy_xxxyyy, tr_xxxy_xxxyyz, tr_xxxy_xxxyzz, tr_xxxy_xxxz, tr_xxxy_xxxzzz, tr_xxxy_xxyy, tr_xxxy_xxyyyy, tr_xxxy_xxyyyz, tr_xxxy_xxyyzz, tr_xxxy_xxyz, tr_xxxy_xxyzzz, tr_xxxy_xxzz, tr_xxxy_xxzzzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxy_xxxx[i] = 6.0 * tr_xy_xxxx[i] + 24.0 * tr_xxy_xxx[i] - 12.0 * tr_xxy_xxxxx[i] * tke_0 + 12.0 * tr_xxxy_xx[i] - 14.0 * tr_xxxy_xxxx[i] * tbe_0 - 18.0 * tr_xxxy_xxxx[i] * tke_0 + 4.0 * tr_xxxy_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxxxy_xxx[i] * tbe_0 + 8.0 * tr_xxxxy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xxxy[i] = 6.0 * tr_xy_xxxy[i] + 18.0 * tr_xxy_xxy[i] - 12.0 * tr_xxy_xxxxy[i] * tke_0 + 6.0 * tr_xxxy_xy[i] - 14.0 * tr_xxxy_xxxy[i] * tbe_0 - 14.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxxxy_xxy[i] * tbe_0 + 8.0 * tr_xxxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xxxz[i] = 6.0 * tr_xy_xxxz[i] + 18.0 * tr_xxy_xxz[i] - 12.0 * tr_xxy_xxxxz[i] * tke_0 + 6.0 * tr_xxxy_xz[i] - 14.0 * tr_xxxy_xxxz[i] * tbe_0 - 14.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxxxy_xxz[i] * tbe_0 + 8.0 * tr_xxxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xxyy[i] = 6.0 * tr_xy_xxyy[i] + 12.0 * tr_xxy_xyy[i] - 12.0 * tr_xxy_xxxyy[i] * tke_0 + 2.0 * tr_xxxy_yy[i] - 14.0 * tr_xxxy_xxyy[i] * tbe_0 - 10.0 * tr_xxxy_xxyy[i] * tke_0 + 4.0 * tr_xxxy_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_xyy[i] * tbe_0 + 8.0 * tr_xxxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xxyz[i] = 6.0 * tr_xy_xxyz[i] + 12.0 * tr_xxy_xyz[i] - 12.0 * tr_xxy_xxxyz[i] * tke_0 + 2.0 * tr_xxxy_yz[i] - 14.0 * tr_xxxy_xxyz[i] * tbe_0 - 10.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_xyz[i] * tbe_0 + 8.0 * tr_xxxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xxzz[i] = 6.0 * tr_xy_xxzz[i] + 12.0 * tr_xxy_xzz[i] - 12.0 * tr_xxy_xxxzz[i] * tke_0 + 2.0 * tr_xxxy_zz[i] - 14.0 * tr_xxxy_xxzz[i] * tbe_0 - 10.0 * tr_xxxy_xxzz[i] * tke_0 + 4.0 * tr_xxxy_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_xzz[i] * tbe_0 + 8.0 * tr_xxxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xyyy[i] = 6.0 * tr_xy_xyyy[i] + 6.0 * tr_xxy_yyy[i] - 12.0 * tr_xxy_xxyyy[i] * tke_0 - 14.0 * tr_xxxy_xyyy[i] * tbe_0 - 6.0 * tr_xxxy_xyyy[i] * tke_0 + 4.0 * tr_xxxy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_yyy[i] * tbe_0 + 8.0 * tr_xxxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xyyz[i] = 6.0 * tr_xy_xyyz[i] + 6.0 * tr_xxy_yyz[i] - 12.0 * tr_xxy_xxyyz[i] * tke_0 - 14.0 * tr_xxxy_xyyz[i] * tbe_0 - 6.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_yyz[i] * tbe_0 + 8.0 * tr_xxxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xyzz[i] = 6.0 * tr_xy_xyzz[i] + 6.0 * tr_xxy_yzz[i] - 12.0 * tr_xxy_xxyzz[i] * tke_0 - 14.0 * tr_xxxy_xyzz[i] * tbe_0 - 6.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_yzz[i] * tbe_0 + 8.0 * tr_xxxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xzzz[i] = 6.0 * tr_xy_xzzz[i] + 6.0 * tr_xxy_zzz[i] - 12.0 * tr_xxy_xxzzz[i] * tke_0 - 14.0 * tr_xxxy_xzzz[i] * tbe_0 - 6.0 * tr_xxxy_xzzz[i] * tke_0 + 4.0 * tr_xxxy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_zzz[i] * tbe_0 + 8.0 * tr_xxxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_yyyy[i] = 6.0 * tr_xy_yyyy[i] - 12.0 * tr_xxy_xyyyy[i] * tke_0 - 14.0 * tr_xxxy_yyyy[i] * tbe_0 - 2.0 * tr_xxxy_yyyy[i] * tke_0 + 4.0 * tr_xxxy_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_yyyz[i] = 6.0 * tr_xy_yyyz[i] - 12.0 * tr_xxy_xyyyz[i] * tke_0 - 14.0 * tr_xxxy_yyyz[i] * tbe_0 - 2.0 * tr_xxxy_yyyz[i] * tke_0 + 4.0 * tr_xxxy_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_yyzz[i] = 6.0 * tr_xy_yyzz[i] - 12.0 * tr_xxy_xyyzz[i] * tke_0 - 14.0 * tr_xxxy_yyzz[i] * tbe_0 - 2.0 * tr_xxxy_yyzz[i] * tke_0 + 4.0 * tr_xxxy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_yzzz[i] = 6.0 * tr_xy_yzzz[i] - 12.0 * tr_xxy_xyzzz[i] * tke_0 - 14.0 * tr_xxxy_yzzz[i] * tbe_0 - 2.0 * tr_xxxy_yzzz[i] * tke_0 + 4.0 * tr_xxxy_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_zzzz[i] = 6.0 * tr_xy_zzzz[i] - 12.0 * tr_xxy_xzzzz[i] * tke_0 - 14.0 * tr_xxxy_zzzz[i] * tbe_0 - 2.0 * tr_xxxy_zzzz[i] * tke_0 + 4.0 * tr_xxxy_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-45 components of targeted buffer : GG

    auto tr_0_0_xx_xxxz_xxxx = pbuffer.data(idx_op_geom_020_gg + 30);

    auto tr_0_0_xx_xxxz_xxxy = pbuffer.data(idx_op_geom_020_gg + 31);

    auto tr_0_0_xx_xxxz_xxxz = pbuffer.data(idx_op_geom_020_gg + 32);

    auto tr_0_0_xx_xxxz_xxyy = pbuffer.data(idx_op_geom_020_gg + 33);

    auto tr_0_0_xx_xxxz_xxyz = pbuffer.data(idx_op_geom_020_gg + 34);

    auto tr_0_0_xx_xxxz_xxzz = pbuffer.data(idx_op_geom_020_gg + 35);

    auto tr_0_0_xx_xxxz_xyyy = pbuffer.data(idx_op_geom_020_gg + 36);

    auto tr_0_0_xx_xxxz_xyyz = pbuffer.data(idx_op_geom_020_gg + 37);

    auto tr_0_0_xx_xxxz_xyzz = pbuffer.data(idx_op_geom_020_gg + 38);

    auto tr_0_0_xx_xxxz_xzzz = pbuffer.data(idx_op_geom_020_gg + 39);

    auto tr_0_0_xx_xxxz_yyyy = pbuffer.data(idx_op_geom_020_gg + 40);

    auto tr_0_0_xx_xxxz_yyyz = pbuffer.data(idx_op_geom_020_gg + 41);

    auto tr_0_0_xx_xxxz_yyzz = pbuffer.data(idx_op_geom_020_gg + 42);

    auto tr_0_0_xx_xxxz_yzzz = pbuffer.data(idx_op_geom_020_gg + 43);

    auto tr_0_0_xx_xxxz_zzzz = pbuffer.data(idx_op_geom_020_gg + 44);

    #pragma omp simd aligned(tr_0_0_xx_xxxz_xxxx, tr_0_0_xx_xxxz_xxxy, tr_0_0_xx_xxxz_xxxz, tr_0_0_xx_xxxz_xxyy, tr_0_0_xx_xxxz_xxyz, tr_0_0_xx_xxxz_xxzz, tr_0_0_xx_xxxz_xyyy, tr_0_0_xx_xxxz_xyyz, tr_0_0_xx_xxxz_xyzz, tr_0_0_xx_xxxz_xzzz, tr_0_0_xx_xxxz_yyyy, tr_0_0_xx_xxxz_yyyz, tr_0_0_xx_xxxz_yyzz, tr_0_0_xx_xxxz_yzzz, tr_0_0_xx_xxxz_zzzz, tr_xxxxxz_xxxx, tr_xxxxxz_xxxy, tr_xxxxxz_xxxz, tr_xxxxxz_xxyy, tr_xxxxxz_xxyz, tr_xxxxxz_xxzz, tr_xxxxxz_xyyy, tr_xxxxxz_xyyz, tr_xxxxxz_xyzz, tr_xxxxxz_xzzz, tr_xxxxxz_yyyy, tr_xxxxxz_yyyz, tr_xxxxxz_yyzz, tr_xxxxxz_yzzz, tr_xxxxxz_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxxxx, tr_xxxxz_xxxxy, tr_xxxxz_xxxxz, tr_xxxxz_xxxyy, tr_xxxxz_xxxyz, tr_xxxxz_xxxzz, tr_xxxxz_xxy, tr_xxxxz_xxyyy, tr_xxxxz_xxyyz, tr_xxxxz_xxyzz, tr_xxxxz_xxz, tr_xxxxz_xxzzz, tr_xxxxz_xyy, tr_xxxxz_xyyyy, tr_xxxxz_xyyyz, tr_xxxxz_xyyzz, tr_xxxxz_xyz, tr_xxxxz_xyzzz, tr_xxxxz_xzz, tr_xxxxz_xzzzz, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_zzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxxxx, tr_xxxz_xxxxxy, tr_xxxz_xxxxxz, tr_xxxz_xxxxyy, tr_xxxz_xxxxyz, tr_xxxz_xxxxzz, tr_xxxz_xxxy, tr_xxxz_xxxyyy, tr_xxxz_xxxyyz, tr_xxxz_xxxyzz, tr_xxxz_xxxz, tr_xxxz_xxxzzz, tr_xxxz_xxyy, tr_xxxz_xxyyyy, tr_xxxz_xxyyyz, tr_xxxz_xxyyzz, tr_xxxz_xxyz, tr_xxxz_xxyzzz, tr_xxxz_xxzz, tr_xxxz_xxzzzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxz_xxxx[i] = 6.0 * tr_xz_xxxx[i] + 24.0 * tr_xxz_xxx[i] - 12.0 * tr_xxz_xxxxx[i] * tke_0 + 12.0 * tr_xxxz_xx[i] - 14.0 * tr_xxxz_xxxx[i] * tbe_0 - 18.0 * tr_xxxz_xxxx[i] * tke_0 + 4.0 * tr_xxxz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxxxz_xxx[i] * tbe_0 + 8.0 * tr_xxxxz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xxxy[i] = 6.0 * tr_xz_xxxy[i] + 18.0 * tr_xxz_xxy[i] - 12.0 * tr_xxz_xxxxy[i] * tke_0 + 6.0 * tr_xxxz_xy[i] - 14.0 * tr_xxxz_xxxy[i] * tbe_0 - 14.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxxxz_xxy[i] * tbe_0 + 8.0 * tr_xxxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xxxz[i] = 6.0 * tr_xz_xxxz[i] + 18.0 * tr_xxz_xxz[i] - 12.0 * tr_xxz_xxxxz[i] * tke_0 + 6.0 * tr_xxxz_xz[i] - 14.0 * tr_xxxz_xxxz[i] * tbe_0 - 14.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxxxz_xxz[i] * tbe_0 + 8.0 * tr_xxxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xxyy[i] = 6.0 * tr_xz_xxyy[i] + 12.0 * tr_xxz_xyy[i] - 12.0 * tr_xxz_xxxyy[i] * tke_0 + 2.0 * tr_xxxz_yy[i] - 14.0 * tr_xxxz_xxyy[i] * tbe_0 - 10.0 * tr_xxxz_xxyy[i] * tke_0 + 4.0 * tr_xxxz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_xyy[i] * tbe_0 + 8.0 * tr_xxxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xxyz[i] = 6.0 * tr_xz_xxyz[i] + 12.0 * tr_xxz_xyz[i] - 12.0 * tr_xxz_xxxyz[i] * tke_0 + 2.0 * tr_xxxz_yz[i] - 14.0 * tr_xxxz_xxyz[i] * tbe_0 - 10.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_xyz[i] * tbe_0 + 8.0 * tr_xxxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xxzz[i] = 6.0 * tr_xz_xxzz[i] + 12.0 * tr_xxz_xzz[i] - 12.0 * tr_xxz_xxxzz[i] * tke_0 + 2.0 * tr_xxxz_zz[i] - 14.0 * tr_xxxz_xxzz[i] * tbe_0 - 10.0 * tr_xxxz_xxzz[i] * tke_0 + 4.0 * tr_xxxz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_xzz[i] * tbe_0 + 8.0 * tr_xxxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xyyy[i] = 6.0 * tr_xz_xyyy[i] + 6.0 * tr_xxz_yyy[i] - 12.0 * tr_xxz_xxyyy[i] * tke_0 - 14.0 * tr_xxxz_xyyy[i] * tbe_0 - 6.0 * tr_xxxz_xyyy[i] * tke_0 + 4.0 * tr_xxxz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_yyy[i] * tbe_0 + 8.0 * tr_xxxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xyyz[i] = 6.0 * tr_xz_xyyz[i] + 6.0 * tr_xxz_yyz[i] - 12.0 * tr_xxz_xxyyz[i] * tke_0 - 14.0 * tr_xxxz_xyyz[i] * tbe_0 - 6.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_yyz[i] * tbe_0 + 8.0 * tr_xxxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xyzz[i] = 6.0 * tr_xz_xyzz[i] + 6.0 * tr_xxz_yzz[i] - 12.0 * tr_xxz_xxyzz[i] * tke_0 - 14.0 * tr_xxxz_xyzz[i] * tbe_0 - 6.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_yzz[i] * tbe_0 + 8.0 * tr_xxxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xzzz[i] = 6.0 * tr_xz_xzzz[i] + 6.0 * tr_xxz_zzz[i] - 12.0 * tr_xxz_xxzzz[i] * tke_0 - 14.0 * tr_xxxz_xzzz[i] * tbe_0 - 6.0 * tr_xxxz_xzzz[i] * tke_0 + 4.0 * tr_xxxz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_zzz[i] * tbe_0 + 8.0 * tr_xxxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_yyyy[i] = 6.0 * tr_xz_yyyy[i] - 12.0 * tr_xxz_xyyyy[i] * tke_0 - 14.0 * tr_xxxz_yyyy[i] * tbe_0 - 2.0 * tr_xxxz_yyyy[i] * tke_0 + 4.0 * tr_xxxz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_yyyz[i] = 6.0 * tr_xz_yyyz[i] - 12.0 * tr_xxz_xyyyz[i] * tke_0 - 14.0 * tr_xxxz_yyyz[i] * tbe_0 - 2.0 * tr_xxxz_yyyz[i] * tke_0 + 4.0 * tr_xxxz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_yyzz[i] = 6.0 * tr_xz_yyzz[i] - 12.0 * tr_xxz_xyyzz[i] * tke_0 - 14.0 * tr_xxxz_yyzz[i] * tbe_0 - 2.0 * tr_xxxz_yyzz[i] * tke_0 + 4.0 * tr_xxxz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_yzzz[i] = 6.0 * tr_xz_yzzz[i] - 12.0 * tr_xxz_xyzzz[i] * tke_0 - 14.0 * tr_xxxz_yzzz[i] * tbe_0 - 2.0 * tr_xxxz_yzzz[i] * tke_0 + 4.0 * tr_xxxz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_zzzz[i] = 6.0 * tr_xz_zzzz[i] - 12.0 * tr_xxz_xzzzz[i] * tke_0 - 14.0 * tr_xxxz_zzzz[i] * tbe_0 - 2.0 * tr_xxxz_zzzz[i] * tke_0 + 4.0 * tr_xxxz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 45-60 components of targeted buffer : GG

    auto tr_0_0_xx_xxyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 45);

    auto tr_0_0_xx_xxyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 46);

    auto tr_0_0_xx_xxyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 47);

    auto tr_0_0_xx_xxyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 48);

    auto tr_0_0_xx_xxyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 49);

    auto tr_0_0_xx_xxyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 50);

    auto tr_0_0_xx_xxyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 51);

    auto tr_0_0_xx_xxyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 52);

    auto tr_0_0_xx_xxyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 53);

    auto tr_0_0_xx_xxyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 54);

    auto tr_0_0_xx_xxyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 55);

    auto tr_0_0_xx_xxyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 56);

    auto tr_0_0_xx_xxyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 57);

    auto tr_0_0_xx_xxyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 58);

    auto tr_0_0_xx_xxyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 59);

    #pragma omp simd aligned(tr_0_0_xx_xxyy_xxxx, tr_0_0_xx_xxyy_xxxy, tr_0_0_xx_xxyy_xxxz, tr_0_0_xx_xxyy_xxyy, tr_0_0_xx_xxyy_xxyz, tr_0_0_xx_xxyy_xxzz, tr_0_0_xx_xxyy_xyyy, tr_0_0_xx_xxyy_xyyz, tr_0_0_xx_xxyy_xyzz, tr_0_0_xx_xxyy_xzzz, tr_0_0_xx_xxyy_yyyy, tr_0_0_xx_xxyy_yyyz, tr_0_0_xx_xxyy_yyzz, tr_0_0_xx_xxyy_yzzz, tr_0_0_xx_xxyy_zzzz, tr_xxxxyy_xxxx, tr_xxxxyy_xxxy, tr_xxxxyy_xxxz, tr_xxxxyy_xxyy, tr_xxxxyy_xxyz, tr_xxxxyy_xxzz, tr_xxxxyy_xyyy, tr_xxxxyy_xyyz, tr_xxxxyy_xyzz, tr_xxxxyy_xzzz, tr_xxxxyy_yyyy, tr_xxxxyy_yyyz, tr_xxxxyy_yyzz, tr_xxxxyy_yzzz, tr_xxxxyy_zzzz, tr_xxxyy_xxx, tr_xxxyy_xxxxx, tr_xxxyy_xxxxy, tr_xxxyy_xxxxz, tr_xxxyy_xxxyy, tr_xxxyy_xxxyz, tr_xxxyy_xxxzz, tr_xxxyy_xxy, tr_xxxyy_xxyyy, tr_xxxyy_xxyyz, tr_xxxyy_xxyzz, tr_xxxyy_xxz, tr_xxxyy_xxzzz, tr_xxxyy_xyy, tr_xxxyy_xyyyy, tr_xxxyy_xyyyz, tr_xxxyy_xyyzz, tr_xxxyy_xyz, tr_xxxyy_xyzzz, tr_xxxyy_xzz, tr_xxxyy_xzzzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxxxx, tr_xxyy_xxxxxy, tr_xxyy_xxxxxz, tr_xxyy_xxxxyy, tr_xxyy_xxxxyz, tr_xxyy_xxxxzz, tr_xxyy_xxxy, tr_xxyy_xxxyyy, tr_xxyy_xxxyyz, tr_xxyy_xxxyzz, tr_xxyy_xxxz, tr_xxyy_xxxzzz, tr_xxyy_xxyy, tr_xxyy_xxyyyy, tr_xxyy_xxyyyz, tr_xxyy_xxyyzz, tr_xxyy_xxyz, tr_xxyy_xxyzzz, tr_xxyy_xxzz, tr_xxyy_xxzzzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxyy_xxxx[i] = 2.0 * tr_yy_xxxx[i] + 16.0 * tr_xyy_xxx[i] - 8.0 * tr_xyy_xxxxx[i] * tke_0 + 12.0 * tr_xxyy_xx[i] - 10.0 * tr_xxyy_xxxx[i] * tbe_0 - 18.0 * tr_xxyy_xxxx[i] * tke_0 + 4.0 * tr_xxyy_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxxyy_xxx[i] * tbe_0 + 8.0 * tr_xxxyy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xxxy[i] = 2.0 * tr_yy_xxxy[i] + 12.0 * tr_xyy_xxy[i] - 8.0 * tr_xyy_xxxxy[i] * tke_0 + 6.0 * tr_xxyy_xy[i] - 10.0 * tr_xxyy_xxxy[i] * tbe_0 - 14.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxxyy_xxy[i] * tbe_0 + 8.0 * tr_xxxyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xxxz[i] = 2.0 * tr_yy_xxxz[i] + 12.0 * tr_xyy_xxz[i] - 8.0 * tr_xyy_xxxxz[i] * tke_0 + 6.0 * tr_xxyy_xz[i] - 10.0 * tr_xxyy_xxxz[i] * tbe_0 - 14.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxxyy_xxz[i] * tbe_0 + 8.0 * tr_xxxyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xxyy[i] = 2.0 * tr_yy_xxyy[i] + 8.0 * tr_xyy_xyy[i] - 8.0 * tr_xyy_xxxyy[i] * tke_0 + 2.0 * tr_xxyy_yy[i] - 10.0 * tr_xxyy_xxyy[i] * tbe_0 - 10.0 * tr_xxyy_xxyy[i] * tke_0 + 4.0 * tr_xxyy_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_xyy[i] * tbe_0 + 8.0 * tr_xxxyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xxyz[i] = 2.0 * tr_yy_xxyz[i] + 8.0 * tr_xyy_xyz[i] - 8.0 * tr_xyy_xxxyz[i] * tke_0 + 2.0 * tr_xxyy_yz[i] - 10.0 * tr_xxyy_xxyz[i] * tbe_0 - 10.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_xyz[i] * tbe_0 + 8.0 * tr_xxxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xxzz[i] = 2.0 * tr_yy_xxzz[i] + 8.0 * tr_xyy_xzz[i] - 8.0 * tr_xyy_xxxzz[i] * tke_0 + 2.0 * tr_xxyy_zz[i] - 10.0 * tr_xxyy_xxzz[i] * tbe_0 - 10.0 * tr_xxyy_xxzz[i] * tke_0 + 4.0 * tr_xxyy_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_xzz[i] * tbe_0 + 8.0 * tr_xxxyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xyyy[i] = 2.0 * tr_yy_xyyy[i] + 4.0 * tr_xyy_yyy[i] - 8.0 * tr_xyy_xxyyy[i] * tke_0 - 10.0 * tr_xxyy_xyyy[i] * tbe_0 - 6.0 * tr_xxyy_xyyy[i] * tke_0 + 4.0 * tr_xxyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_yyy[i] * tbe_0 + 8.0 * tr_xxxyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xyyz[i] = 2.0 * tr_yy_xyyz[i] + 4.0 * tr_xyy_yyz[i] - 8.0 * tr_xyy_xxyyz[i] * tke_0 - 10.0 * tr_xxyy_xyyz[i] * tbe_0 - 6.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_yyz[i] * tbe_0 + 8.0 * tr_xxxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xyzz[i] = 2.0 * tr_yy_xyzz[i] + 4.0 * tr_xyy_yzz[i] - 8.0 * tr_xyy_xxyzz[i] * tke_0 - 10.0 * tr_xxyy_xyzz[i] * tbe_0 - 6.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_yzz[i] * tbe_0 + 8.0 * tr_xxxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xzzz[i] = 2.0 * tr_yy_xzzz[i] + 4.0 * tr_xyy_zzz[i] - 8.0 * tr_xyy_xxzzz[i] * tke_0 - 10.0 * tr_xxyy_xzzz[i] * tbe_0 - 6.0 * tr_xxyy_xzzz[i] * tke_0 + 4.0 * tr_xxyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_zzz[i] * tbe_0 + 8.0 * tr_xxxyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_yyyy[i] = 2.0 * tr_yy_yyyy[i] - 8.0 * tr_xyy_xyyyy[i] * tke_0 - 10.0 * tr_xxyy_yyyy[i] * tbe_0 - 2.0 * tr_xxyy_yyyy[i] * tke_0 + 4.0 * tr_xxyy_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_yyyz[i] = 2.0 * tr_yy_yyyz[i] - 8.0 * tr_xyy_xyyyz[i] * tke_0 - 10.0 * tr_xxyy_yyyz[i] * tbe_0 - 2.0 * tr_xxyy_yyyz[i] * tke_0 + 4.0 * tr_xxyy_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_yyzz[i] = 2.0 * tr_yy_yyzz[i] - 8.0 * tr_xyy_xyyzz[i] * tke_0 - 10.0 * tr_xxyy_yyzz[i] * tbe_0 - 2.0 * tr_xxyy_yyzz[i] * tke_0 + 4.0 * tr_xxyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_yzzz[i] = 2.0 * tr_yy_yzzz[i] - 8.0 * tr_xyy_xyzzz[i] * tke_0 - 10.0 * tr_xxyy_yzzz[i] * tbe_0 - 2.0 * tr_xxyy_yzzz[i] * tke_0 + 4.0 * tr_xxyy_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_zzzz[i] = 2.0 * tr_yy_zzzz[i] - 8.0 * tr_xyy_xzzzz[i] * tke_0 - 10.0 * tr_xxyy_zzzz[i] * tbe_0 - 2.0 * tr_xxyy_zzzz[i] * tke_0 + 4.0 * tr_xxyy_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-75 components of targeted buffer : GG

    auto tr_0_0_xx_xxyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 60);

    auto tr_0_0_xx_xxyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 61);

    auto tr_0_0_xx_xxyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 62);

    auto tr_0_0_xx_xxyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 63);

    auto tr_0_0_xx_xxyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 64);

    auto tr_0_0_xx_xxyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 65);

    auto tr_0_0_xx_xxyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 66);

    auto tr_0_0_xx_xxyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 67);

    auto tr_0_0_xx_xxyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 68);

    auto tr_0_0_xx_xxyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 69);

    auto tr_0_0_xx_xxyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 70);

    auto tr_0_0_xx_xxyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 71);

    auto tr_0_0_xx_xxyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 72);

    auto tr_0_0_xx_xxyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 73);

    auto tr_0_0_xx_xxyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 74);

    #pragma omp simd aligned(tr_0_0_xx_xxyz_xxxx, tr_0_0_xx_xxyz_xxxy, tr_0_0_xx_xxyz_xxxz, tr_0_0_xx_xxyz_xxyy, tr_0_0_xx_xxyz_xxyz, tr_0_0_xx_xxyz_xxzz, tr_0_0_xx_xxyz_xyyy, tr_0_0_xx_xxyz_xyyz, tr_0_0_xx_xxyz_xyzz, tr_0_0_xx_xxyz_xzzz, tr_0_0_xx_xxyz_yyyy, tr_0_0_xx_xxyz_yyyz, tr_0_0_xx_xxyz_yyzz, tr_0_0_xx_xxyz_yzzz, tr_0_0_xx_xxyz_zzzz, tr_xxxxyz_xxxx, tr_xxxxyz_xxxy, tr_xxxxyz_xxxz, tr_xxxxyz_xxyy, tr_xxxxyz_xxyz, tr_xxxxyz_xxzz, tr_xxxxyz_xyyy, tr_xxxxyz_xyyz, tr_xxxxyz_xyzz, tr_xxxxyz_xzzz, tr_xxxxyz_yyyy, tr_xxxxyz_yyyz, tr_xxxxyz_yyzz, tr_xxxxyz_yzzz, tr_xxxxyz_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxx, tr_xxxyz_xxxxy, tr_xxxyz_xxxxz, tr_xxxyz_xxxyy, tr_xxxyz_xxxyz, tr_xxxyz_xxxzz, tr_xxxyz_xxy, tr_xxxyz_xxyyy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xxzzz, tr_xxxyz_xyy, tr_xxxyz_xyyyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_xzzzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxxxx, tr_xxyz_xxxxxy, tr_xxyz_xxxxxz, tr_xxyz_xxxxyy, tr_xxyz_xxxxyz, tr_xxyz_xxxxzz, tr_xxyz_xxxy, tr_xxyz_xxxyyy, tr_xxyz_xxxyyz, tr_xxyz_xxxyzz, tr_xxyz_xxxz, tr_xxyz_xxxzzz, tr_xxyz_xxyy, tr_xxyz_xxyyyy, tr_xxyz_xxyyyz, tr_xxyz_xxyyzz, tr_xxyz_xxyz, tr_xxyz_xxyzzz, tr_xxyz_xxzz, tr_xxyz_xxzzzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxyz_xxxx[i] = 2.0 * tr_yz_xxxx[i] + 16.0 * tr_xyz_xxx[i] - 8.0 * tr_xyz_xxxxx[i] * tke_0 + 12.0 * tr_xxyz_xx[i] - 10.0 * tr_xxyz_xxxx[i] * tbe_0 - 18.0 * tr_xxyz_xxxx[i] * tke_0 + 4.0 * tr_xxyz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxxyz_xxx[i] * tbe_0 + 8.0 * tr_xxxyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xxxy[i] = 2.0 * tr_yz_xxxy[i] + 12.0 * tr_xyz_xxy[i] - 8.0 * tr_xyz_xxxxy[i] * tke_0 + 6.0 * tr_xxyz_xy[i] - 10.0 * tr_xxyz_xxxy[i] * tbe_0 - 14.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxxyz_xxy[i] * tbe_0 + 8.0 * tr_xxxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xxxz[i] = 2.0 * tr_yz_xxxz[i] + 12.0 * tr_xyz_xxz[i] - 8.0 * tr_xyz_xxxxz[i] * tke_0 + 6.0 * tr_xxyz_xz[i] - 10.0 * tr_xxyz_xxxz[i] * tbe_0 - 14.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxxyz_xxz[i] * tbe_0 + 8.0 * tr_xxxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xxyy[i] = 2.0 * tr_yz_xxyy[i] + 8.0 * tr_xyz_xyy[i] - 8.0 * tr_xyz_xxxyy[i] * tke_0 + 2.0 * tr_xxyz_yy[i] - 10.0 * tr_xxyz_xxyy[i] * tbe_0 - 10.0 * tr_xxyz_xxyy[i] * tke_0 + 4.0 * tr_xxyz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xyy[i] * tbe_0 + 8.0 * tr_xxxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xxyz[i] = 2.0 * tr_yz_xxyz[i] + 8.0 * tr_xyz_xyz[i] - 8.0 * tr_xyz_xxxyz[i] * tke_0 + 2.0 * tr_xxyz_yz[i] - 10.0 * tr_xxyz_xxyz[i] * tbe_0 - 10.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xyz[i] * tbe_0 + 8.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xxzz[i] = 2.0 * tr_yz_xxzz[i] + 8.0 * tr_xyz_xzz[i] - 8.0 * tr_xyz_xxxzz[i] * tke_0 + 2.0 * tr_xxyz_zz[i] - 10.0 * tr_xxyz_xxzz[i] * tbe_0 - 10.0 * tr_xxyz_xxzz[i] * tke_0 + 4.0 * tr_xxyz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xzz[i] * tbe_0 + 8.0 * tr_xxxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xyyy[i] = 2.0 * tr_yz_xyyy[i] + 4.0 * tr_xyz_yyy[i] - 8.0 * tr_xyz_xxyyy[i] * tke_0 - 10.0 * tr_xxyz_xyyy[i] * tbe_0 - 6.0 * tr_xxyz_xyyy[i] * tke_0 + 4.0 * tr_xxyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_yyy[i] * tbe_0 + 8.0 * tr_xxxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xyyz[i] = 2.0 * tr_yz_xyyz[i] + 4.0 * tr_xyz_yyz[i] - 8.0 * tr_xyz_xxyyz[i] * tke_0 - 10.0 * tr_xxyz_xyyz[i] * tbe_0 - 6.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_yyz[i] * tbe_0 + 8.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xyzz[i] = 2.0 * tr_yz_xyzz[i] + 4.0 * tr_xyz_yzz[i] - 8.0 * tr_xyz_xxyzz[i] * tke_0 - 10.0 * tr_xxyz_xyzz[i] * tbe_0 - 6.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_yzz[i] * tbe_0 + 8.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xzzz[i] = 2.0 * tr_yz_xzzz[i] + 4.0 * tr_xyz_zzz[i] - 8.0 * tr_xyz_xxzzz[i] * tke_0 - 10.0 * tr_xxyz_xzzz[i] * tbe_0 - 6.0 * tr_xxyz_xzzz[i] * tke_0 + 4.0 * tr_xxyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_zzz[i] * tbe_0 + 8.0 * tr_xxxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_yyyy[i] = 2.0 * tr_yz_yyyy[i] - 8.0 * tr_xyz_xyyyy[i] * tke_0 - 10.0 * tr_xxyz_yyyy[i] * tbe_0 - 2.0 * tr_xxyz_yyyy[i] * tke_0 + 4.0 * tr_xxyz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_yyyz[i] = 2.0 * tr_yz_yyyz[i] - 8.0 * tr_xyz_xyyyz[i] * tke_0 - 10.0 * tr_xxyz_yyyz[i] * tbe_0 - 2.0 * tr_xxyz_yyyz[i] * tke_0 + 4.0 * tr_xxyz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_yyzz[i] = 2.0 * tr_yz_yyzz[i] - 8.0 * tr_xyz_xyyzz[i] * tke_0 - 10.0 * tr_xxyz_yyzz[i] * tbe_0 - 2.0 * tr_xxyz_yyzz[i] * tke_0 + 4.0 * tr_xxyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_yzzz[i] = 2.0 * tr_yz_yzzz[i] - 8.0 * tr_xyz_xyzzz[i] * tke_0 - 10.0 * tr_xxyz_yzzz[i] * tbe_0 - 2.0 * tr_xxyz_yzzz[i] * tke_0 + 4.0 * tr_xxyz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_zzzz[i] = 2.0 * tr_yz_zzzz[i] - 8.0 * tr_xyz_xzzzz[i] * tke_0 - 10.0 * tr_xxyz_zzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzzz[i] * tke_0 + 4.0 * tr_xxyz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 75-90 components of targeted buffer : GG

    auto tr_0_0_xx_xxzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 75);

    auto tr_0_0_xx_xxzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 76);

    auto tr_0_0_xx_xxzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 77);

    auto tr_0_0_xx_xxzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 78);

    auto tr_0_0_xx_xxzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 79);

    auto tr_0_0_xx_xxzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 80);

    auto tr_0_0_xx_xxzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 81);

    auto tr_0_0_xx_xxzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 82);

    auto tr_0_0_xx_xxzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 83);

    auto tr_0_0_xx_xxzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 84);

    auto tr_0_0_xx_xxzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 85);

    auto tr_0_0_xx_xxzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 86);

    auto tr_0_0_xx_xxzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 87);

    auto tr_0_0_xx_xxzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 88);

    auto tr_0_0_xx_xxzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 89);

    #pragma omp simd aligned(tr_0_0_xx_xxzz_xxxx, tr_0_0_xx_xxzz_xxxy, tr_0_0_xx_xxzz_xxxz, tr_0_0_xx_xxzz_xxyy, tr_0_0_xx_xxzz_xxyz, tr_0_0_xx_xxzz_xxzz, tr_0_0_xx_xxzz_xyyy, tr_0_0_xx_xxzz_xyyz, tr_0_0_xx_xxzz_xyzz, tr_0_0_xx_xxzz_xzzz, tr_0_0_xx_xxzz_yyyy, tr_0_0_xx_xxzz_yyyz, tr_0_0_xx_xxzz_yyzz, tr_0_0_xx_xxzz_yzzz, tr_0_0_xx_xxzz_zzzz, tr_xxxxzz_xxxx, tr_xxxxzz_xxxy, tr_xxxxzz_xxxz, tr_xxxxzz_xxyy, tr_xxxxzz_xxyz, tr_xxxxzz_xxzz, tr_xxxxzz_xyyy, tr_xxxxzz_xyyz, tr_xxxxzz_xyzz, tr_xxxxzz_xzzz, tr_xxxxzz_yyyy, tr_xxxxzz_yyyz, tr_xxxxzz_yyzz, tr_xxxxzz_yzzz, tr_xxxxzz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxxxx, tr_xxxzz_xxxxy, tr_xxxzz_xxxxz, tr_xxxzz_xxxyy, tr_xxxzz_xxxyz, tr_xxxzz_xxxzz, tr_xxxzz_xxy, tr_xxxzz_xxyyy, tr_xxxzz_xxyyz, tr_xxxzz_xxyzz, tr_xxxzz_xxz, tr_xxxzz_xxzzz, tr_xxxzz_xyy, tr_xxxzz_xyyyy, tr_xxxzz_xyyyz, tr_xxxzz_xyyzz, tr_xxxzz_xyz, tr_xxxzz_xyzzz, tr_xxxzz_xzz, tr_xxxzz_xzzzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxxxx, tr_xxzz_xxxxxy, tr_xxzz_xxxxxz, tr_xxzz_xxxxyy, tr_xxzz_xxxxyz, tr_xxzz_xxxxzz, tr_xxzz_xxxy, tr_xxzz_xxxyyy, tr_xxzz_xxxyyz, tr_xxzz_xxxyzz, tr_xxzz_xxxz, tr_xxzz_xxxzzz, tr_xxzz_xxyy, tr_xxzz_xxyyyy, tr_xxzz_xxyyyz, tr_xxzz_xxyyzz, tr_xxzz_xxyz, tr_xxzz_xxyzzz, tr_xxzz_xxzz, tr_xxzz_xxzzzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxzz_xxxx[i] = 2.0 * tr_zz_xxxx[i] + 16.0 * tr_xzz_xxx[i] - 8.0 * tr_xzz_xxxxx[i] * tke_0 + 12.0 * tr_xxzz_xx[i] - 10.0 * tr_xxzz_xxxx[i] * tbe_0 - 18.0 * tr_xxzz_xxxx[i] * tke_0 + 4.0 * tr_xxzz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxxzz_xxx[i] * tbe_0 + 8.0 * tr_xxxzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xxxy[i] = 2.0 * tr_zz_xxxy[i] + 12.0 * tr_xzz_xxy[i] - 8.0 * tr_xzz_xxxxy[i] * tke_0 + 6.0 * tr_xxzz_xy[i] - 10.0 * tr_xxzz_xxxy[i] * tbe_0 - 14.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxxzz_xxy[i] * tbe_0 + 8.0 * tr_xxxzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xxxz[i] = 2.0 * tr_zz_xxxz[i] + 12.0 * tr_xzz_xxz[i] - 8.0 * tr_xzz_xxxxz[i] * tke_0 + 6.0 * tr_xxzz_xz[i] - 10.0 * tr_xxzz_xxxz[i] * tbe_0 - 14.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxxzz_xxz[i] * tbe_0 + 8.0 * tr_xxxzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xxyy[i] = 2.0 * tr_zz_xxyy[i] + 8.0 * tr_xzz_xyy[i] - 8.0 * tr_xzz_xxxyy[i] * tke_0 + 2.0 * tr_xxzz_yy[i] - 10.0 * tr_xxzz_xxyy[i] * tbe_0 - 10.0 * tr_xxzz_xxyy[i] * tke_0 + 4.0 * tr_xxzz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_xyy[i] * tbe_0 + 8.0 * tr_xxxzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xxyz[i] = 2.0 * tr_zz_xxyz[i] + 8.0 * tr_xzz_xyz[i] - 8.0 * tr_xzz_xxxyz[i] * tke_0 + 2.0 * tr_xxzz_yz[i] - 10.0 * tr_xxzz_xxyz[i] * tbe_0 - 10.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_xyz[i] * tbe_0 + 8.0 * tr_xxxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xxzz[i] = 2.0 * tr_zz_xxzz[i] + 8.0 * tr_xzz_xzz[i] - 8.0 * tr_xzz_xxxzz[i] * tke_0 + 2.0 * tr_xxzz_zz[i] - 10.0 * tr_xxzz_xxzz[i] * tbe_0 - 10.0 * tr_xxzz_xxzz[i] * tke_0 + 4.0 * tr_xxzz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_xzz[i] * tbe_0 + 8.0 * tr_xxxzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xyyy[i] = 2.0 * tr_zz_xyyy[i] + 4.0 * tr_xzz_yyy[i] - 8.0 * tr_xzz_xxyyy[i] * tke_0 - 10.0 * tr_xxzz_xyyy[i] * tbe_0 - 6.0 * tr_xxzz_xyyy[i] * tke_0 + 4.0 * tr_xxzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_yyy[i] * tbe_0 + 8.0 * tr_xxxzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xyyz[i] = 2.0 * tr_zz_xyyz[i] + 4.0 * tr_xzz_yyz[i] - 8.0 * tr_xzz_xxyyz[i] * tke_0 - 10.0 * tr_xxzz_xyyz[i] * tbe_0 - 6.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_yyz[i] * tbe_0 + 8.0 * tr_xxxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xyzz[i] = 2.0 * tr_zz_xyzz[i] + 4.0 * tr_xzz_yzz[i] - 8.0 * tr_xzz_xxyzz[i] * tke_0 - 10.0 * tr_xxzz_xyzz[i] * tbe_0 - 6.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_yzz[i] * tbe_0 + 8.0 * tr_xxxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xzzz[i] = 2.0 * tr_zz_xzzz[i] + 4.0 * tr_xzz_zzz[i] - 8.0 * tr_xzz_xxzzz[i] * tke_0 - 10.0 * tr_xxzz_xzzz[i] * tbe_0 - 6.0 * tr_xxzz_xzzz[i] * tke_0 + 4.0 * tr_xxzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_zzz[i] * tbe_0 + 8.0 * tr_xxxzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_yyyy[i] = 2.0 * tr_zz_yyyy[i] - 8.0 * tr_xzz_xyyyy[i] * tke_0 - 10.0 * tr_xxzz_yyyy[i] * tbe_0 - 2.0 * tr_xxzz_yyyy[i] * tke_0 + 4.0 * tr_xxzz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_yyyz[i] = 2.0 * tr_zz_yyyz[i] - 8.0 * tr_xzz_xyyyz[i] * tke_0 - 10.0 * tr_xxzz_yyyz[i] * tbe_0 - 2.0 * tr_xxzz_yyyz[i] * tke_0 + 4.0 * tr_xxzz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_yyzz[i] = 2.0 * tr_zz_yyzz[i] - 8.0 * tr_xzz_xyyzz[i] * tke_0 - 10.0 * tr_xxzz_yyzz[i] * tbe_0 - 2.0 * tr_xxzz_yyzz[i] * tke_0 + 4.0 * tr_xxzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_yzzz[i] = 2.0 * tr_zz_yzzz[i] - 8.0 * tr_xzz_xyzzz[i] * tke_0 - 10.0 * tr_xxzz_yzzz[i] * tbe_0 - 2.0 * tr_xxzz_yzzz[i] * tke_0 + 4.0 * tr_xxzz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_zzzz[i] = 2.0 * tr_zz_zzzz[i] - 8.0 * tr_xzz_xzzzz[i] * tke_0 - 10.0 * tr_xxzz_zzzz[i] * tbe_0 - 2.0 * tr_xxzz_zzzz[i] * tke_0 + 4.0 * tr_xxzz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-105 components of targeted buffer : GG

    auto tr_0_0_xx_xyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 90);

    auto tr_0_0_xx_xyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 91);

    auto tr_0_0_xx_xyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 92);

    auto tr_0_0_xx_xyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 93);

    auto tr_0_0_xx_xyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 94);

    auto tr_0_0_xx_xyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 95);

    auto tr_0_0_xx_xyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 96);

    auto tr_0_0_xx_xyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 97);

    auto tr_0_0_xx_xyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 98);

    auto tr_0_0_xx_xyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 99);

    auto tr_0_0_xx_xyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 100);

    auto tr_0_0_xx_xyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 101);

    auto tr_0_0_xx_xyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 102);

    auto tr_0_0_xx_xyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 103);

    auto tr_0_0_xx_xyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 104);

    #pragma omp simd aligned(tr_0_0_xx_xyyy_xxxx, tr_0_0_xx_xyyy_xxxy, tr_0_0_xx_xyyy_xxxz, tr_0_0_xx_xyyy_xxyy, tr_0_0_xx_xyyy_xxyz, tr_0_0_xx_xyyy_xxzz, tr_0_0_xx_xyyy_xyyy, tr_0_0_xx_xyyy_xyyz, tr_0_0_xx_xyyy_xyzz, tr_0_0_xx_xyyy_xzzz, tr_0_0_xx_xyyy_yyyy, tr_0_0_xx_xyyy_yyyz, tr_0_0_xx_xyyy_yyzz, tr_0_0_xx_xyyy_yzzz, tr_0_0_xx_xyyy_zzzz, tr_xxxyyy_xxxx, tr_xxxyyy_xxxy, tr_xxxyyy_xxxz, tr_xxxyyy_xxyy, tr_xxxyyy_xxyz, tr_xxxyyy_xxzz, tr_xxxyyy_xyyy, tr_xxxyyy_xyyz, tr_xxxyyy_xyzz, tr_xxxyyy_xzzz, tr_xxxyyy_yyyy, tr_xxxyyy_yyyz, tr_xxxyyy_yyzz, tr_xxxyyy_yzzz, tr_xxxyyy_zzzz, tr_xxyyy_xxx, tr_xxyyy_xxxxx, tr_xxyyy_xxxxy, tr_xxyyy_xxxxz, tr_xxyyy_xxxyy, tr_xxyyy_xxxyz, tr_xxyyy_xxxzz, tr_xxyyy_xxy, tr_xxyyy_xxyyy, tr_xxyyy_xxyyz, tr_xxyyy_xxyzz, tr_xxyyy_xxz, tr_xxyyy_xxzzz, tr_xxyyy_xyy, tr_xxyyy_xyyyy, tr_xxyyy_xyyyz, tr_xxyyy_xyyzz, tr_xxyyy_xyz, tr_xxyyy_xyzzz, tr_xxyyy_xzz, tr_xxyyy_xzzzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxxxx, tr_xyyy_xxxxxy, tr_xyyy_xxxxxz, tr_xyyy_xxxxyy, tr_xyyy_xxxxyz, tr_xyyy_xxxxzz, tr_xyyy_xxxy, tr_xyyy_xxxyyy, tr_xyyy_xxxyyz, tr_xyyy_xxxyzz, tr_xyyy_xxxz, tr_xyyy_xxxzzz, tr_xyyy_xxyy, tr_xyyy_xxyyyy, tr_xyyy_xxyyyz, tr_xyyy_xxyyzz, tr_xyyy_xxyz, tr_xyyy_xxyzzz, tr_xyyy_xxzz, tr_xyyy_xxzzzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_yyy_xxx, tr_yyy_xxxxx, tr_yyy_xxxxy, tr_yyy_xxxxz, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyyy_xxxx[i] = 8.0 * tr_yyy_xxx[i] - 4.0 * tr_yyy_xxxxx[i] * tke_0 + 12.0 * tr_xyyy_xx[i] - 6.0 * tr_xyyy_xxxx[i] * tbe_0 - 18.0 * tr_xyyy_xxxx[i] * tke_0 + 4.0 * tr_xyyy_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxyyy_xxx[i] * tbe_0 + 8.0 * tr_xxyyy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xxxy[i] = 6.0 * tr_yyy_xxy[i] - 4.0 * tr_yyy_xxxxy[i] * tke_0 + 6.0 * tr_xyyy_xy[i] - 6.0 * tr_xyyy_xxxy[i] * tbe_0 - 14.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxyyy_xxy[i] * tbe_0 + 8.0 * tr_xxyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xxxz[i] = 6.0 * tr_yyy_xxz[i] - 4.0 * tr_yyy_xxxxz[i] * tke_0 + 6.0 * tr_xyyy_xz[i] - 6.0 * tr_xyyy_xxxz[i] * tbe_0 - 14.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxyyy_xxz[i] * tbe_0 + 8.0 * tr_xxyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xxyy[i] = 4.0 * tr_yyy_xyy[i] - 4.0 * tr_yyy_xxxyy[i] * tke_0 + 2.0 * tr_xyyy_yy[i] - 6.0 * tr_xyyy_xxyy[i] * tbe_0 - 10.0 * tr_xyyy_xxyy[i] * tke_0 + 4.0 * tr_xyyy_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_xyy[i] * tbe_0 + 8.0 * tr_xxyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xxyz[i] = 4.0 * tr_yyy_xyz[i] - 4.0 * tr_yyy_xxxyz[i] * tke_0 + 2.0 * tr_xyyy_yz[i] - 6.0 * tr_xyyy_xxyz[i] * tbe_0 - 10.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_xyz[i] * tbe_0 + 8.0 * tr_xxyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xxzz[i] = 4.0 * tr_yyy_xzz[i] - 4.0 * tr_yyy_xxxzz[i] * tke_0 + 2.0 * tr_xyyy_zz[i] - 6.0 * tr_xyyy_xxzz[i] * tbe_0 - 10.0 * tr_xyyy_xxzz[i] * tke_0 + 4.0 * tr_xyyy_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_xzz[i] * tbe_0 + 8.0 * tr_xxyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xyyy[i] = 2.0 * tr_yyy_yyy[i] - 4.0 * tr_yyy_xxyyy[i] * tke_0 - 6.0 * tr_xyyy_xyyy[i] * tbe_0 - 6.0 * tr_xyyy_xyyy[i] * tke_0 + 4.0 * tr_xyyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_yyy[i] * tbe_0 + 8.0 * tr_xxyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xyyz[i] = 2.0 * tr_yyy_yyz[i] - 4.0 * tr_yyy_xxyyz[i] * tke_0 - 6.0 * tr_xyyy_xyyz[i] * tbe_0 - 6.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_yyz[i] * tbe_0 + 8.0 * tr_xxyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xyzz[i] = 2.0 * tr_yyy_yzz[i] - 4.0 * tr_yyy_xxyzz[i] * tke_0 - 6.0 * tr_xyyy_xyzz[i] * tbe_0 - 6.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_yzz[i] * tbe_0 + 8.0 * tr_xxyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xzzz[i] = 2.0 * tr_yyy_zzz[i] - 4.0 * tr_yyy_xxzzz[i] * tke_0 - 6.0 * tr_xyyy_xzzz[i] * tbe_0 - 6.0 * tr_xyyy_xzzz[i] * tke_0 + 4.0 * tr_xyyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_zzz[i] * tbe_0 + 8.0 * tr_xxyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_yyyy[i] = -4.0 * tr_yyy_xyyyy[i] * tke_0 - 6.0 * tr_xyyy_yyyy[i] * tbe_0 - 2.0 * tr_xyyy_yyyy[i] * tke_0 + 4.0 * tr_xyyy_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_yyyz[i] = -4.0 * tr_yyy_xyyyz[i] * tke_0 - 6.0 * tr_xyyy_yyyz[i] * tbe_0 - 2.0 * tr_xyyy_yyyz[i] * tke_0 + 4.0 * tr_xyyy_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_yyzz[i] = -4.0 * tr_yyy_xyyzz[i] * tke_0 - 6.0 * tr_xyyy_yyzz[i] * tbe_0 - 2.0 * tr_xyyy_yyzz[i] * tke_0 + 4.0 * tr_xyyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_yzzz[i] = -4.0 * tr_yyy_xyzzz[i] * tke_0 - 6.0 * tr_xyyy_yzzz[i] * tbe_0 - 2.0 * tr_xyyy_yzzz[i] * tke_0 + 4.0 * tr_xyyy_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_zzzz[i] = -4.0 * tr_yyy_xzzzz[i] * tke_0 - 6.0 * tr_xyyy_zzzz[i] * tbe_0 - 2.0 * tr_xyyy_zzzz[i] * tke_0 + 4.0 * tr_xyyy_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 105-120 components of targeted buffer : GG

    auto tr_0_0_xx_xyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 105);

    auto tr_0_0_xx_xyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 106);

    auto tr_0_0_xx_xyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 107);

    auto tr_0_0_xx_xyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 108);

    auto tr_0_0_xx_xyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 109);

    auto tr_0_0_xx_xyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 110);

    auto tr_0_0_xx_xyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 111);

    auto tr_0_0_xx_xyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 112);

    auto tr_0_0_xx_xyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 113);

    auto tr_0_0_xx_xyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 114);

    auto tr_0_0_xx_xyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 115);

    auto tr_0_0_xx_xyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 116);

    auto tr_0_0_xx_xyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 117);

    auto tr_0_0_xx_xyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 118);

    auto tr_0_0_xx_xyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 119);

    #pragma omp simd aligned(tr_0_0_xx_xyyz_xxxx, tr_0_0_xx_xyyz_xxxy, tr_0_0_xx_xyyz_xxxz, tr_0_0_xx_xyyz_xxyy, tr_0_0_xx_xyyz_xxyz, tr_0_0_xx_xyyz_xxzz, tr_0_0_xx_xyyz_xyyy, tr_0_0_xx_xyyz_xyyz, tr_0_0_xx_xyyz_xyzz, tr_0_0_xx_xyyz_xzzz, tr_0_0_xx_xyyz_yyyy, tr_0_0_xx_xyyz_yyyz, tr_0_0_xx_xyyz_yyzz, tr_0_0_xx_xyyz_yzzz, tr_0_0_xx_xyyz_zzzz, tr_xxxyyz_xxxx, tr_xxxyyz_xxxy, tr_xxxyyz_xxxz, tr_xxxyyz_xxyy, tr_xxxyyz_xxyz, tr_xxxyyz_xxzz, tr_xxxyyz_xyyy, tr_xxxyyz_xyyz, tr_xxxyyz_xyzz, tr_xxxyyz_xzzz, tr_xxxyyz_yyyy, tr_xxxyyz_yyyz, tr_xxxyyz_yyzz, tr_xxxyyz_yzzz, tr_xxxyyz_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxx, tr_xxyyz_xxxxy, tr_xxyyz_xxxxz, tr_xxyyz_xxxyy, tr_xxyyz_xxxyz, tr_xxyyz_xxxzz, tr_xxyyz_xxy, tr_xxyyz_xxyyy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xxzzz, tr_xxyyz_xyy, tr_xxyyz_xyyyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_xzzzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxxxx, tr_xyyz_xxxxxy, tr_xyyz_xxxxxz, tr_xyyz_xxxxyy, tr_xyyz_xxxxyz, tr_xyyz_xxxxzz, tr_xyyz_xxxy, tr_xyyz_xxxyyy, tr_xyyz_xxxyyz, tr_xyyz_xxxyzz, tr_xyyz_xxxz, tr_xyyz_xxxzzz, tr_xyyz_xxyy, tr_xyyz_xxyyyy, tr_xyyz_xxyyyz, tr_xyyz_xxyyzz, tr_xyyz_xxyz, tr_xyyz_xxyzzz, tr_xyyz_xxzz, tr_xyyz_xxzzzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyyz_xxxx[i] = 8.0 * tr_yyz_xxx[i] - 4.0 * tr_yyz_xxxxx[i] * tke_0 + 12.0 * tr_xyyz_xx[i] - 6.0 * tr_xyyz_xxxx[i] * tbe_0 - 18.0 * tr_xyyz_xxxx[i] * tke_0 + 4.0 * tr_xyyz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxyyz_xxx[i] * tbe_0 + 8.0 * tr_xxyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xxxy[i] = 6.0 * tr_yyz_xxy[i] - 4.0 * tr_yyz_xxxxy[i] * tke_0 + 6.0 * tr_xyyz_xy[i] - 6.0 * tr_xyyz_xxxy[i] * tbe_0 - 14.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxyyz_xxy[i] * tbe_0 + 8.0 * tr_xxyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xxxz[i] = 6.0 * tr_yyz_xxz[i] - 4.0 * tr_yyz_xxxxz[i] * tke_0 + 6.0 * tr_xyyz_xz[i] - 6.0 * tr_xyyz_xxxz[i] * tbe_0 - 14.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxyyz_xxz[i] * tbe_0 + 8.0 * tr_xxyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xxyy[i] = 4.0 * tr_yyz_xyy[i] - 4.0 * tr_yyz_xxxyy[i] * tke_0 + 2.0 * tr_xyyz_yy[i] - 6.0 * tr_xyyz_xxyy[i] * tbe_0 - 10.0 * tr_xyyz_xxyy[i] * tke_0 + 4.0 * tr_xyyz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xyy[i] * tbe_0 + 8.0 * tr_xxyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xxyz[i] = 4.0 * tr_yyz_xyz[i] - 4.0 * tr_yyz_xxxyz[i] * tke_0 + 2.0 * tr_xyyz_yz[i] - 6.0 * tr_xyyz_xxyz[i] * tbe_0 - 10.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xyz[i] * tbe_0 + 8.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xxzz[i] = 4.0 * tr_yyz_xzz[i] - 4.0 * tr_yyz_xxxzz[i] * tke_0 + 2.0 * tr_xyyz_zz[i] - 6.0 * tr_xyyz_xxzz[i] * tbe_0 - 10.0 * tr_xyyz_xxzz[i] * tke_0 + 4.0 * tr_xyyz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xzz[i] * tbe_0 + 8.0 * tr_xxyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xyyy[i] = 2.0 * tr_yyz_yyy[i] - 4.0 * tr_yyz_xxyyy[i] * tke_0 - 6.0 * tr_xyyz_xyyy[i] * tbe_0 - 6.0 * tr_xyyz_xyyy[i] * tke_0 + 4.0 * tr_xyyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_yyy[i] * tbe_0 + 8.0 * tr_xxyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xyyz[i] = 2.0 * tr_yyz_yyz[i] - 4.0 * tr_yyz_xxyyz[i] * tke_0 - 6.0 * tr_xyyz_xyyz[i] * tbe_0 - 6.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_yyz[i] * tbe_0 + 8.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xyzz[i] = 2.0 * tr_yyz_yzz[i] - 4.0 * tr_yyz_xxyzz[i] * tke_0 - 6.0 * tr_xyyz_xyzz[i] * tbe_0 - 6.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_yzz[i] * tbe_0 + 8.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xzzz[i] = 2.0 * tr_yyz_zzz[i] - 4.0 * tr_yyz_xxzzz[i] * tke_0 - 6.0 * tr_xyyz_xzzz[i] * tbe_0 - 6.0 * tr_xyyz_xzzz[i] * tke_0 + 4.0 * tr_xyyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_zzz[i] * tbe_0 + 8.0 * tr_xxyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_yyyy[i] = -4.0 * tr_yyz_xyyyy[i] * tke_0 - 6.0 * tr_xyyz_yyyy[i] * tbe_0 - 2.0 * tr_xyyz_yyyy[i] * tke_0 + 4.0 * tr_xyyz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_yyyz[i] = -4.0 * tr_yyz_xyyyz[i] * tke_0 - 6.0 * tr_xyyz_yyyz[i] * tbe_0 - 2.0 * tr_xyyz_yyyz[i] * tke_0 + 4.0 * tr_xyyz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_yyzz[i] = -4.0 * tr_yyz_xyyzz[i] * tke_0 - 6.0 * tr_xyyz_yyzz[i] * tbe_0 - 2.0 * tr_xyyz_yyzz[i] * tke_0 + 4.0 * tr_xyyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_yzzz[i] = -4.0 * tr_yyz_xyzzz[i] * tke_0 - 6.0 * tr_xyyz_yzzz[i] * tbe_0 - 2.0 * tr_xyyz_yzzz[i] * tke_0 + 4.0 * tr_xyyz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_zzzz[i] = -4.0 * tr_yyz_xzzzz[i] * tke_0 - 6.0 * tr_xyyz_zzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzzz[i] * tke_0 + 4.0 * tr_xyyz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-135 components of targeted buffer : GG

    auto tr_0_0_xx_xyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 120);

    auto tr_0_0_xx_xyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 121);

    auto tr_0_0_xx_xyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 122);

    auto tr_0_0_xx_xyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 123);

    auto tr_0_0_xx_xyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 124);

    auto tr_0_0_xx_xyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 125);

    auto tr_0_0_xx_xyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 126);

    auto tr_0_0_xx_xyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 127);

    auto tr_0_0_xx_xyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 128);

    auto tr_0_0_xx_xyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 129);

    auto tr_0_0_xx_xyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 130);

    auto tr_0_0_xx_xyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 131);

    auto tr_0_0_xx_xyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 132);

    auto tr_0_0_xx_xyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 133);

    auto tr_0_0_xx_xyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 134);

    #pragma omp simd aligned(tr_0_0_xx_xyzz_xxxx, tr_0_0_xx_xyzz_xxxy, tr_0_0_xx_xyzz_xxxz, tr_0_0_xx_xyzz_xxyy, tr_0_0_xx_xyzz_xxyz, tr_0_0_xx_xyzz_xxzz, tr_0_0_xx_xyzz_xyyy, tr_0_0_xx_xyzz_xyyz, tr_0_0_xx_xyzz_xyzz, tr_0_0_xx_xyzz_xzzz, tr_0_0_xx_xyzz_yyyy, tr_0_0_xx_xyzz_yyyz, tr_0_0_xx_xyzz_yyzz, tr_0_0_xx_xyzz_yzzz, tr_0_0_xx_xyzz_zzzz, tr_xxxyzz_xxxx, tr_xxxyzz_xxxy, tr_xxxyzz_xxxz, tr_xxxyzz_xxyy, tr_xxxyzz_xxyz, tr_xxxyzz_xxzz, tr_xxxyzz_xyyy, tr_xxxyzz_xyyz, tr_xxxyzz_xyzz, tr_xxxyzz_xzzz, tr_xxxyzz_yyyy, tr_xxxyzz_yyyz, tr_xxxyzz_yyzz, tr_xxxyzz_yzzz, tr_xxxyzz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxx, tr_xxyzz_xxxxy, tr_xxyzz_xxxxz, tr_xxyzz_xxxyy, tr_xxyzz_xxxyz, tr_xxyzz_xxxzz, tr_xxyzz_xxy, tr_xxyzz_xxyyy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xxzzz, tr_xxyzz_xyy, tr_xxyzz_xyyyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_xzzzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxxxx, tr_xyzz_xxxxxy, tr_xyzz_xxxxxz, tr_xyzz_xxxxyy, tr_xyzz_xxxxyz, tr_xyzz_xxxxzz, tr_xyzz_xxxy, tr_xyzz_xxxyyy, tr_xyzz_xxxyyz, tr_xyzz_xxxyzz, tr_xyzz_xxxz, tr_xyzz_xxxzzz, tr_xyzz_xxyy, tr_xyzz_xxyyyy, tr_xyzz_xxyyyz, tr_xyzz_xxyyzz, tr_xyzz_xxyz, tr_xyzz_xxyzzz, tr_xyzz_xxzz, tr_xyzz_xxzzzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyzz_xxxx[i] = 8.0 * tr_yzz_xxx[i] - 4.0 * tr_yzz_xxxxx[i] * tke_0 + 12.0 * tr_xyzz_xx[i] - 6.0 * tr_xyzz_xxxx[i] * tbe_0 - 18.0 * tr_xyzz_xxxx[i] * tke_0 + 4.0 * tr_xyzz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxyzz_xxx[i] * tbe_0 + 8.0 * tr_xxyzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xxxy[i] = 6.0 * tr_yzz_xxy[i] - 4.0 * tr_yzz_xxxxy[i] * tke_0 + 6.0 * tr_xyzz_xy[i] - 6.0 * tr_xyzz_xxxy[i] * tbe_0 - 14.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxyzz_xxy[i] * tbe_0 + 8.0 * tr_xxyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xxxz[i] = 6.0 * tr_yzz_xxz[i] - 4.0 * tr_yzz_xxxxz[i] * tke_0 + 6.0 * tr_xyzz_xz[i] - 6.0 * tr_xyzz_xxxz[i] * tbe_0 - 14.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxyzz_xxz[i] * tbe_0 + 8.0 * tr_xxyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xxyy[i] = 4.0 * tr_yzz_xyy[i] - 4.0 * tr_yzz_xxxyy[i] * tke_0 + 2.0 * tr_xyzz_yy[i] - 6.0 * tr_xyzz_xxyy[i] * tbe_0 - 10.0 * tr_xyzz_xxyy[i] * tke_0 + 4.0 * tr_xyzz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xyy[i] * tbe_0 + 8.0 * tr_xxyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xxyz[i] = 4.0 * tr_yzz_xyz[i] - 4.0 * tr_yzz_xxxyz[i] * tke_0 + 2.0 * tr_xyzz_yz[i] - 6.0 * tr_xyzz_xxyz[i] * tbe_0 - 10.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xyz[i] * tbe_0 + 8.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xxzz[i] = 4.0 * tr_yzz_xzz[i] - 4.0 * tr_yzz_xxxzz[i] * tke_0 + 2.0 * tr_xyzz_zz[i] - 6.0 * tr_xyzz_xxzz[i] * tbe_0 - 10.0 * tr_xyzz_xxzz[i] * tke_0 + 4.0 * tr_xyzz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xzz[i] * tbe_0 + 8.0 * tr_xxyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xyyy[i] = 2.0 * tr_yzz_yyy[i] - 4.0 * tr_yzz_xxyyy[i] * tke_0 - 6.0 * tr_xyzz_xyyy[i] * tbe_0 - 6.0 * tr_xyzz_xyyy[i] * tke_0 + 4.0 * tr_xyzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_yyy[i] * tbe_0 + 8.0 * tr_xxyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xyyz[i] = 2.0 * tr_yzz_yyz[i] - 4.0 * tr_yzz_xxyyz[i] * tke_0 - 6.0 * tr_xyzz_xyyz[i] * tbe_0 - 6.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_yyz[i] * tbe_0 + 8.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xyzz[i] = 2.0 * tr_yzz_yzz[i] - 4.0 * tr_yzz_xxyzz[i] * tke_0 - 6.0 * tr_xyzz_xyzz[i] * tbe_0 - 6.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_yzz[i] * tbe_0 + 8.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xzzz[i] = 2.0 * tr_yzz_zzz[i] - 4.0 * tr_yzz_xxzzz[i] * tke_0 - 6.0 * tr_xyzz_xzzz[i] * tbe_0 - 6.0 * tr_xyzz_xzzz[i] * tke_0 + 4.0 * tr_xyzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_zzz[i] * tbe_0 + 8.0 * tr_xxyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_yyyy[i] = -4.0 * tr_yzz_xyyyy[i] * tke_0 - 6.0 * tr_xyzz_yyyy[i] * tbe_0 - 2.0 * tr_xyzz_yyyy[i] * tke_0 + 4.0 * tr_xyzz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_yyyz[i] = -4.0 * tr_yzz_xyyyz[i] * tke_0 - 6.0 * tr_xyzz_yyyz[i] * tbe_0 - 2.0 * tr_xyzz_yyyz[i] * tke_0 + 4.0 * tr_xyzz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_yyzz[i] = -4.0 * tr_yzz_xyyzz[i] * tke_0 - 6.0 * tr_xyzz_yyzz[i] * tbe_0 - 2.0 * tr_xyzz_yyzz[i] * tke_0 + 4.0 * tr_xyzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_yzzz[i] = -4.0 * tr_yzz_xyzzz[i] * tke_0 - 6.0 * tr_xyzz_yzzz[i] * tbe_0 - 2.0 * tr_xyzz_yzzz[i] * tke_0 + 4.0 * tr_xyzz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_zzzz[i] = -4.0 * tr_yzz_xzzzz[i] * tke_0 - 6.0 * tr_xyzz_zzzz[i] * tbe_0 - 2.0 * tr_xyzz_zzzz[i] * tke_0 + 4.0 * tr_xyzz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 135-150 components of targeted buffer : GG

    auto tr_0_0_xx_xzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 135);

    auto tr_0_0_xx_xzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 136);

    auto tr_0_0_xx_xzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 137);

    auto tr_0_0_xx_xzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 138);

    auto tr_0_0_xx_xzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 139);

    auto tr_0_0_xx_xzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 140);

    auto tr_0_0_xx_xzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 141);

    auto tr_0_0_xx_xzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 142);

    auto tr_0_0_xx_xzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 143);

    auto tr_0_0_xx_xzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 144);

    auto tr_0_0_xx_xzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 145);

    auto tr_0_0_xx_xzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 146);

    auto tr_0_0_xx_xzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 147);

    auto tr_0_0_xx_xzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 148);

    auto tr_0_0_xx_xzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 149);

    #pragma omp simd aligned(tr_0_0_xx_xzzz_xxxx, tr_0_0_xx_xzzz_xxxy, tr_0_0_xx_xzzz_xxxz, tr_0_0_xx_xzzz_xxyy, tr_0_0_xx_xzzz_xxyz, tr_0_0_xx_xzzz_xxzz, tr_0_0_xx_xzzz_xyyy, tr_0_0_xx_xzzz_xyyz, tr_0_0_xx_xzzz_xyzz, tr_0_0_xx_xzzz_xzzz, tr_0_0_xx_xzzz_yyyy, tr_0_0_xx_xzzz_yyyz, tr_0_0_xx_xzzz_yyzz, tr_0_0_xx_xzzz_yzzz, tr_0_0_xx_xzzz_zzzz, tr_xxxzzz_xxxx, tr_xxxzzz_xxxy, tr_xxxzzz_xxxz, tr_xxxzzz_xxyy, tr_xxxzzz_xxyz, tr_xxxzzz_xxzz, tr_xxxzzz_xyyy, tr_xxxzzz_xyyz, tr_xxxzzz_xyzz, tr_xxxzzz_xzzz, tr_xxxzzz_yyyy, tr_xxxzzz_yyyz, tr_xxxzzz_yyzz, tr_xxxzzz_yzzz, tr_xxxzzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxxxx, tr_xxzzz_xxxxy, tr_xxzzz_xxxxz, tr_xxzzz_xxxyy, tr_xxzzz_xxxyz, tr_xxzzz_xxxzz, tr_xxzzz_xxy, tr_xxzzz_xxyyy, tr_xxzzz_xxyyz, tr_xxzzz_xxyzz, tr_xxzzz_xxz, tr_xxzzz_xxzzz, tr_xxzzz_xyy, tr_xxzzz_xyyyy, tr_xxzzz_xyyyz, tr_xxzzz_xyyzz, tr_xxzzz_xyz, tr_xxzzz_xyzzz, tr_xxzzz_xzz, tr_xxzzz_xzzzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxxxx, tr_xzzz_xxxxxy, tr_xzzz_xxxxxz, tr_xzzz_xxxxyy, tr_xzzz_xxxxyz, tr_xzzz_xxxxzz, tr_xzzz_xxxy, tr_xzzz_xxxyyy, tr_xzzz_xxxyyz, tr_xzzz_xxxyzz, tr_xzzz_xxxz, tr_xzzz_xxxzzz, tr_xzzz_xxyy, tr_xzzz_xxyyyy, tr_xzzz_xxyyyz, tr_xzzz_xxyyzz, tr_xzzz_xxyz, tr_xzzz_xxyzzz, tr_xzzz_xxzz, tr_xzzz_xxzzzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxx, tr_zzz_xxxxy, tr_zzz_xxxxz, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xzzz_xxxx[i] = 8.0 * tr_zzz_xxx[i] - 4.0 * tr_zzz_xxxxx[i] * tke_0 + 12.0 * tr_xzzz_xx[i] - 6.0 * tr_xzzz_xxxx[i] * tbe_0 - 18.0 * tr_xzzz_xxxx[i] * tke_0 + 4.0 * tr_xzzz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xxzzz_xxx[i] * tbe_0 + 8.0 * tr_xxzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xxxy[i] = 6.0 * tr_zzz_xxy[i] - 4.0 * tr_zzz_xxxxy[i] * tke_0 + 6.0 * tr_xzzz_xy[i] - 6.0 * tr_xzzz_xxxy[i] * tbe_0 - 14.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xxzzz_xxy[i] * tbe_0 + 8.0 * tr_xxzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xxxz[i] = 6.0 * tr_zzz_xxz[i] - 4.0 * tr_zzz_xxxxz[i] * tke_0 + 6.0 * tr_xzzz_xz[i] - 6.0 * tr_xzzz_xxxz[i] * tbe_0 - 14.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xxzzz_xxz[i] * tbe_0 + 8.0 * tr_xxzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xxyy[i] = 4.0 * tr_zzz_xyy[i] - 4.0 * tr_zzz_xxxyy[i] * tke_0 + 2.0 * tr_xzzz_yy[i] - 6.0 * tr_xzzz_xxyy[i] * tbe_0 - 10.0 * tr_xzzz_xxyy[i] * tke_0 + 4.0 * tr_xzzz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_xyy[i] * tbe_0 + 8.0 * tr_xxzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xxyz[i] = 4.0 * tr_zzz_xyz[i] - 4.0 * tr_zzz_xxxyz[i] * tke_0 + 2.0 * tr_xzzz_yz[i] - 6.0 * tr_xzzz_xxyz[i] * tbe_0 - 10.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_xyz[i] * tbe_0 + 8.0 * tr_xxzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xxzz[i] = 4.0 * tr_zzz_xzz[i] - 4.0 * tr_zzz_xxxzz[i] * tke_0 + 2.0 * tr_xzzz_zz[i] - 6.0 * tr_xzzz_xxzz[i] * tbe_0 - 10.0 * tr_xzzz_xxzz[i] * tke_0 + 4.0 * tr_xzzz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_xzz[i] * tbe_0 + 8.0 * tr_xxzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xyyy[i] = 2.0 * tr_zzz_yyy[i] - 4.0 * tr_zzz_xxyyy[i] * tke_0 - 6.0 * tr_xzzz_xyyy[i] * tbe_0 - 6.0 * tr_xzzz_xyyy[i] * tke_0 + 4.0 * tr_xzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_yyy[i] * tbe_0 + 8.0 * tr_xxzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xyyz[i] = 2.0 * tr_zzz_yyz[i] - 4.0 * tr_zzz_xxyyz[i] * tke_0 - 6.0 * tr_xzzz_xyyz[i] * tbe_0 - 6.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_yyz[i] * tbe_0 + 8.0 * tr_xxzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xyzz[i] = 2.0 * tr_zzz_yzz[i] - 4.0 * tr_zzz_xxyzz[i] * tke_0 - 6.0 * tr_xzzz_xyzz[i] * tbe_0 - 6.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_yzz[i] * tbe_0 + 8.0 * tr_xxzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xzzz[i] = 2.0 * tr_zzz_zzz[i] - 4.0 * tr_zzz_xxzzz[i] * tke_0 - 6.0 * tr_xzzz_xzzz[i] * tbe_0 - 6.0 * tr_xzzz_xzzz[i] * tke_0 + 4.0 * tr_xzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_zzz[i] * tbe_0 + 8.0 * tr_xxzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_yyyy[i] = -4.0 * tr_zzz_xyyyy[i] * tke_0 - 6.0 * tr_xzzz_yyyy[i] * tbe_0 - 2.0 * tr_xzzz_yyyy[i] * tke_0 + 4.0 * tr_xzzz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_yyyz[i] = -4.0 * tr_zzz_xyyyz[i] * tke_0 - 6.0 * tr_xzzz_yyyz[i] * tbe_0 - 2.0 * tr_xzzz_yyyz[i] * tke_0 + 4.0 * tr_xzzz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_yyzz[i] = -4.0 * tr_zzz_xyyzz[i] * tke_0 - 6.0 * tr_xzzz_yyzz[i] * tbe_0 - 2.0 * tr_xzzz_yyzz[i] * tke_0 + 4.0 * tr_xzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_yzzz[i] = -4.0 * tr_zzz_xyzzz[i] * tke_0 - 6.0 * tr_xzzz_yzzz[i] * tbe_0 - 2.0 * tr_xzzz_yzzz[i] * tke_0 + 4.0 * tr_xzzz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_zzzz[i] = -4.0 * tr_zzz_xzzzz[i] * tke_0 - 6.0 * tr_xzzz_zzzz[i] * tbe_0 - 2.0 * tr_xzzz_zzzz[i] * tke_0 + 4.0 * tr_xzzz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-165 components of targeted buffer : GG

    auto tr_0_0_xx_yyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 150);

    auto tr_0_0_xx_yyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 151);

    auto tr_0_0_xx_yyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 152);

    auto tr_0_0_xx_yyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 153);

    auto tr_0_0_xx_yyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 154);

    auto tr_0_0_xx_yyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 155);

    auto tr_0_0_xx_yyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 156);

    auto tr_0_0_xx_yyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 157);

    auto tr_0_0_xx_yyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 158);

    auto tr_0_0_xx_yyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 159);

    auto tr_0_0_xx_yyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 160);

    auto tr_0_0_xx_yyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 161);

    auto tr_0_0_xx_yyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 162);

    auto tr_0_0_xx_yyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 163);

    auto tr_0_0_xx_yyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 164);

    #pragma omp simd aligned(tr_0_0_xx_yyyy_xxxx, tr_0_0_xx_yyyy_xxxy, tr_0_0_xx_yyyy_xxxz, tr_0_0_xx_yyyy_xxyy, tr_0_0_xx_yyyy_xxyz, tr_0_0_xx_yyyy_xxzz, tr_0_0_xx_yyyy_xyyy, tr_0_0_xx_yyyy_xyyz, tr_0_0_xx_yyyy_xyzz, tr_0_0_xx_yyyy_xzzz, tr_0_0_xx_yyyy_yyyy, tr_0_0_xx_yyyy_yyyz, tr_0_0_xx_yyyy_yyzz, tr_0_0_xx_yyyy_yzzz, tr_0_0_xx_yyyy_zzzz, tr_xxyyyy_xxxx, tr_xxyyyy_xxxy, tr_xxyyyy_xxxz, tr_xxyyyy_xxyy, tr_xxyyyy_xxyz, tr_xxyyyy_xxzz, tr_xxyyyy_xyyy, tr_xxyyyy_xyyz, tr_xxyyyy_xyzz, tr_xxyyyy_xzzz, tr_xxyyyy_yyyy, tr_xxyyyy_yyyz, tr_xxyyyy_yyzz, tr_xxyyyy_yzzz, tr_xxyyyy_zzzz, tr_xyyyy_xxx, tr_xyyyy_xxxxx, tr_xyyyy_xxxxy, tr_xyyyy_xxxxz, tr_xyyyy_xxxyy, tr_xyyyy_xxxyz, tr_xyyyy_xxxzz, tr_xyyyy_xxy, tr_xyyyy_xxyyy, tr_xyyyy_xxyyz, tr_xyyyy_xxyzz, tr_xyyyy_xxz, tr_xyyyy_xxzzz, tr_xyyyy_xyy, tr_xyyyy_xyyyy, tr_xyyyy_xyyyz, tr_xyyyy_xyyzz, tr_xyyyy_xyz, tr_xyyyy_xyzzz, tr_xyyyy_xzz, tr_xyyyy_xzzzz, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_zzz, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxxxx, tr_yyyy_xxxxxy, tr_yyyy_xxxxxz, tr_yyyy_xxxxyy, tr_yyyy_xxxxyz, tr_yyyy_xxxxzz, tr_yyyy_xxxy, tr_yyyy_xxxyyy, tr_yyyy_xxxyyz, tr_yyyy_xxxyzz, tr_yyyy_xxxz, tr_yyyy_xxxzzz, tr_yyyy_xxyy, tr_yyyy_xxyyyy, tr_yyyy_xxyyyz, tr_yyyy_xxyyzz, tr_yyyy_xxyz, tr_yyyy_xxyzzz, tr_yyyy_xxzz, tr_yyyy_xxzzzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyyy_xxxx[i] = 12.0 * tr_yyyy_xx[i] - 2.0 * tr_yyyy_xxxx[i] * tbe_0 - 18.0 * tr_yyyy_xxxx[i] * tke_0 + 4.0 * tr_yyyy_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xyyyy_xxx[i] * tbe_0 + 8.0 * tr_xyyyy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xxxy[i] = 6.0 * tr_yyyy_xy[i] - 2.0 * tr_yyyy_xxxy[i] * tbe_0 - 14.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xyyyy_xxy[i] * tbe_0 + 8.0 * tr_xyyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xxxz[i] = 6.0 * tr_yyyy_xz[i] - 2.0 * tr_yyyy_xxxz[i] * tbe_0 - 14.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xyyyy_xxz[i] * tbe_0 + 8.0 * tr_xyyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xxyy[i] = 2.0 * tr_yyyy_yy[i] - 2.0 * tr_yyyy_xxyy[i] * tbe_0 - 10.0 * tr_yyyy_xxyy[i] * tke_0 + 4.0 * tr_yyyy_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_xyy[i] * tbe_0 + 8.0 * tr_xyyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xxyz[i] = 2.0 * tr_yyyy_yz[i] - 2.0 * tr_yyyy_xxyz[i] * tbe_0 - 10.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_xyz[i] * tbe_0 + 8.0 * tr_xyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xxzz[i] = 2.0 * tr_yyyy_zz[i] - 2.0 * tr_yyyy_xxzz[i] * tbe_0 - 10.0 * tr_yyyy_xxzz[i] * tke_0 + 4.0 * tr_yyyy_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_xzz[i] * tbe_0 + 8.0 * tr_xyyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xyyy[i] = -2.0 * tr_yyyy_xyyy[i] * tbe_0 - 6.0 * tr_yyyy_xyyy[i] * tke_0 + 4.0 * tr_yyyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_yyy[i] * tbe_0 + 8.0 * tr_xyyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xyyz[i] = -2.0 * tr_yyyy_xyyz[i] * tbe_0 - 6.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_yyz[i] * tbe_0 + 8.0 * tr_xyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xyzz[i] = -2.0 * tr_yyyy_xyzz[i] * tbe_0 - 6.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_yzz[i] * tbe_0 + 8.0 * tr_xyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xzzz[i] = -2.0 * tr_yyyy_xzzz[i] * tbe_0 - 6.0 * tr_yyyy_xzzz[i] * tke_0 + 4.0 * tr_yyyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_zzz[i] * tbe_0 + 8.0 * tr_xyyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_yyyy[i] = -2.0 * tr_yyyy_yyyy[i] * tbe_0 - 2.0 * tr_yyyy_yyyy[i] * tke_0 + 4.0 * tr_yyyy_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_yyyz[i] = -2.0 * tr_yyyy_yyyz[i] * tbe_0 - 2.0 * tr_yyyy_yyyz[i] * tke_0 + 4.0 * tr_yyyy_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_yyzz[i] = -2.0 * tr_yyyy_yyzz[i] * tbe_0 - 2.0 * tr_yyyy_yyzz[i] * tke_0 + 4.0 * tr_yyyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_yzzz[i] = -2.0 * tr_yyyy_yzzz[i] * tbe_0 - 2.0 * tr_yyyy_yzzz[i] * tke_0 + 4.0 * tr_yyyy_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_zzzz[i] = -2.0 * tr_yyyy_zzzz[i] * tbe_0 - 2.0 * tr_yyyy_zzzz[i] * tke_0 + 4.0 * tr_yyyy_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 165-180 components of targeted buffer : GG

    auto tr_0_0_xx_yyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 165);

    auto tr_0_0_xx_yyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 166);

    auto tr_0_0_xx_yyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 167);

    auto tr_0_0_xx_yyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 168);

    auto tr_0_0_xx_yyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 169);

    auto tr_0_0_xx_yyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 170);

    auto tr_0_0_xx_yyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 171);

    auto tr_0_0_xx_yyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 172);

    auto tr_0_0_xx_yyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 173);

    auto tr_0_0_xx_yyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 174);

    auto tr_0_0_xx_yyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 175);

    auto tr_0_0_xx_yyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 176);

    auto tr_0_0_xx_yyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 177);

    auto tr_0_0_xx_yyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 178);

    auto tr_0_0_xx_yyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 179);

    #pragma omp simd aligned(tr_0_0_xx_yyyz_xxxx, tr_0_0_xx_yyyz_xxxy, tr_0_0_xx_yyyz_xxxz, tr_0_0_xx_yyyz_xxyy, tr_0_0_xx_yyyz_xxyz, tr_0_0_xx_yyyz_xxzz, tr_0_0_xx_yyyz_xyyy, tr_0_0_xx_yyyz_xyyz, tr_0_0_xx_yyyz_xyzz, tr_0_0_xx_yyyz_xzzz, tr_0_0_xx_yyyz_yyyy, tr_0_0_xx_yyyz_yyyz, tr_0_0_xx_yyyz_yyzz, tr_0_0_xx_yyyz_yzzz, tr_0_0_xx_yyyz_zzzz, tr_xxyyyz_xxxx, tr_xxyyyz_xxxy, tr_xxyyyz_xxxz, tr_xxyyyz_xxyy, tr_xxyyyz_xxyz, tr_xxyyyz_xxzz, tr_xxyyyz_xyyy, tr_xxyyyz_xyyz, tr_xxyyyz_xyzz, tr_xxyyyz_xzzz, tr_xxyyyz_yyyy, tr_xxyyyz_yyyz, tr_xxyyyz_yyzz, tr_xxyyyz_yzzz, tr_xxyyyz_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxx, tr_xyyyz_xxxxy, tr_xyyyz_xxxxz, tr_xyyyz_xxxyy, tr_xyyyz_xxxyz, tr_xyyyz_xxxzz, tr_xyyyz_xxy, tr_xyyyz_xxyyy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xxzzz, tr_xyyyz_xyy, tr_xyyyz_xyyyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_xzzzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxxxx, tr_yyyz_xxxxxy, tr_yyyz_xxxxxz, tr_yyyz_xxxxyy, tr_yyyz_xxxxyz, tr_yyyz_xxxxzz, tr_yyyz_xxxy, tr_yyyz_xxxyyy, tr_yyyz_xxxyyz, tr_yyyz_xxxyzz, tr_yyyz_xxxz, tr_yyyz_xxxzzz, tr_yyyz_xxyy, tr_yyyz_xxyyyy, tr_yyyz_xxyyyz, tr_yyyz_xxyyzz, tr_yyyz_xxyz, tr_yyyz_xxyzzz, tr_yyyz_xxzz, tr_yyyz_xxzzzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyyz_xxxx[i] = 12.0 * tr_yyyz_xx[i] - 2.0 * tr_yyyz_xxxx[i] * tbe_0 - 18.0 * tr_yyyz_xxxx[i] * tke_0 + 4.0 * tr_yyyz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xyyyz_xxx[i] * tbe_0 + 8.0 * tr_xyyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xxxy[i] = 6.0 * tr_yyyz_xy[i] - 2.0 * tr_yyyz_xxxy[i] * tbe_0 - 14.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xyyyz_xxy[i] * tbe_0 + 8.0 * tr_xyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xxxz[i] = 6.0 * tr_yyyz_xz[i] - 2.0 * tr_yyyz_xxxz[i] * tbe_0 - 14.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xyyyz_xxz[i] * tbe_0 + 8.0 * tr_xyyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xxyy[i] = 2.0 * tr_yyyz_yy[i] - 2.0 * tr_yyyz_xxyy[i] * tbe_0 - 10.0 * tr_yyyz_xxyy[i] * tke_0 + 4.0 * tr_yyyz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xyy[i] * tbe_0 + 8.0 * tr_xyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xxyz[i] = 2.0 * tr_yyyz_yz[i] - 2.0 * tr_yyyz_xxyz[i] * tbe_0 - 10.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xyz[i] * tbe_0 + 8.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xxzz[i] = 2.0 * tr_yyyz_zz[i] - 2.0 * tr_yyyz_xxzz[i] * tbe_0 - 10.0 * tr_yyyz_xxzz[i] * tke_0 + 4.0 * tr_yyyz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xzz[i] * tbe_0 + 8.0 * tr_xyyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xyyy[i] = -2.0 * tr_yyyz_xyyy[i] * tbe_0 - 6.0 * tr_yyyz_xyyy[i] * tke_0 + 4.0 * tr_yyyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_yyy[i] * tbe_0 + 8.0 * tr_xyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xyyz[i] = -2.0 * tr_yyyz_xyyz[i] * tbe_0 - 6.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_yyz[i] * tbe_0 + 8.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xyzz[i] = -2.0 * tr_yyyz_xyzz[i] * tbe_0 - 6.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_yzz[i] * tbe_0 + 8.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xzzz[i] = -2.0 * tr_yyyz_xzzz[i] * tbe_0 - 6.0 * tr_yyyz_xzzz[i] * tke_0 + 4.0 * tr_yyyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_zzz[i] * tbe_0 + 8.0 * tr_xyyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_yyyy[i] = -2.0 * tr_yyyz_yyyy[i] * tbe_0 - 2.0 * tr_yyyz_yyyy[i] * tke_0 + 4.0 * tr_yyyz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_yyyz[i] = -2.0 * tr_yyyz_yyyz[i] * tbe_0 - 2.0 * tr_yyyz_yyyz[i] * tke_0 + 4.0 * tr_yyyz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_yyzz[i] = -2.0 * tr_yyyz_yyzz[i] * tbe_0 - 2.0 * tr_yyyz_yyzz[i] * tke_0 + 4.0 * tr_yyyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_yzzz[i] = -2.0 * tr_yyyz_yzzz[i] * tbe_0 - 2.0 * tr_yyyz_yzzz[i] * tke_0 + 4.0 * tr_yyyz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_zzzz[i] = -2.0 * tr_yyyz_zzzz[i] * tbe_0 - 2.0 * tr_yyyz_zzzz[i] * tke_0 + 4.0 * tr_yyyz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-195 components of targeted buffer : GG

    auto tr_0_0_xx_yyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 180);

    auto tr_0_0_xx_yyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 181);

    auto tr_0_0_xx_yyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 182);

    auto tr_0_0_xx_yyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 183);

    auto tr_0_0_xx_yyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 184);

    auto tr_0_0_xx_yyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 185);

    auto tr_0_0_xx_yyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 186);

    auto tr_0_0_xx_yyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 187);

    auto tr_0_0_xx_yyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 188);

    auto tr_0_0_xx_yyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 189);

    auto tr_0_0_xx_yyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 190);

    auto tr_0_0_xx_yyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 191);

    auto tr_0_0_xx_yyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 192);

    auto tr_0_0_xx_yyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 193);

    auto tr_0_0_xx_yyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 194);

    #pragma omp simd aligned(tr_0_0_xx_yyzz_xxxx, tr_0_0_xx_yyzz_xxxy, tr_0_0_xx_yyzz_xxxz, tr_0_0_xx_yyzz_xxyy, tr_0_0_xx_yyzz_xxyz, tr_0_0_xx_yyzz_xxzz, tr_0_0_xx_yyzz_xyyy, tr_0_0_xx_yyzz_xyyz, tr_0_0_xx_yyzz_xyzz, tr_0_0_xx_yyzz_xzzz, tr_0_0_xx_yyzz_yyyy, tr_0_0_xx_yyzz_yyyz, tr_0_0_xx_yyzz_yyzz, tr_0_0_xx_yyzz_yzzz, tr_0_0_xx_yyzz_zzzz, tr_xxyyzz_xxxx, tr_xxyyzz_xxxy, tr_xxyyzz_xxxz, tr_xxyyzz_xxyy, tr_xxyyzz_xxyz, tr_xxyyzz_xxzz, tr_xxyyzz_xyyy, tr_xxyyzz_xyyz, tr_xxyyzz_xyzz, tr_xxyyzz_xzzz, tr_xxyyzz_yyyy, tr_xxyyzz_yyyz, tr_xxyyzz_yyzz, tr_xxyyzz_yzzz, tr_xxyyzz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxx, tr_xyyzz_xxxxy, tr_xyyzz_xxxxz, tr_xyyzz_xxxyy, tr_xyyzz_xxxyz, tr_xyyzz_xxxzz, tr_xyyzz_xxy, tr_xyyzz_xxyyy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xxzzz, tr_xyyzz_xyy, tr_xyyzz_xyyyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_xzzzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxxxx, tr_yyzz_xxxxxy, tr_yyzz_xxxxxz, tr_yyzz_xxxxyy, tr_yyzz_xxxxyz, tr_yyzz_xxxxzz, tr_yyzz_xxxy, tr_yyzz_xxxyyy, tr_yyzz_xxxyyz, tr_yyzz_xxxyzz, tr_yyzz_xxxz, tr_yyzz_xxxzzz, tr_yyzz_xxyy, tr_yyzz_xxyyyy, tr_yyzz_xxyyyz, tr_yyzz_xxyyzz, tr_yyzz_xxyz, tr_yyzz_xxyzzz, tr_yyzz_xxzz, tr_yyzz_xxzzzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyzz_xxxx[i] = 12.0 * tr_yyzz_xx[i] - 2.0 * tr_yyzz_xxxx[i] * tbe_0 - 18.0 * tr_yyzz_xxxx[i] * tke_0 + 4.0 * tr_yyzz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xyyzz_xxx[i] * tbe_0 + 8.0 * tr_xyyzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xxxy[i] = 6.0 * tr_yyzz_xy[i] - 2.0 * tr_yyzz_xxxy[i] * tbe_0 - 14.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xyyzz_xxy[i] * tbe_0 + 8.0 * tr_xyyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xxxz[i] = 6.0 * tr_yyzz_xz[i] - 2.0 * tr_yyzz_xxxz[i] * tbe_0 - 14.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xyyzz_xxz[i] * tbe_0 + 8.0 * tr_xyyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xxyy[i] = 2.0 * tr_yyzz_yy[i] - 2.0 * tr_yyzz_xxyy[i] * tbe_0 - 10.0 * tr_yyzz_xxyy[i] * tke_0 + 4.0 * tr_yyzz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xyy[i] * tbe_0 + 8.0 * tr_xyyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xxyz[i] = 2.0 * tr_yyzz_yz[i] - 2.0 * tr_yyzz_xxyz[i] * tbe_0 - 10.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xyz[i] * tbe_0 + 8.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xxzz[i] = 2.0 * tr_yyzz_zz[i] - 2.0 * tr_yyzz_xxzz[i] * tbe_0 - 10.0 * tr_yyzz_xxzz[i] * tke_0 + 4.0 * tr_yyzz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xzz[i] * tbe_0 + 8.0 * tr_xyyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xyyy[i] = -2.0 * tr_yyzz_xyyy[i] * tbe_0 - 6.0 * tr_yyzz_xyyy[i] * tke_0 + 4.0 * tr_yyzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_yyy[i] * tbe_0 + 8.0 * tr_xyyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xyyz[i] = -2.0 * tr_yyzz_xyyz[i] * tbe_0 - 6.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_yyz[i] * tbe_0 + 8.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xyzz[i] = -2.0 * tr_yyzz_xyzz[i] * tbe_0 - 6.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_yzz[i] * tbe_0 + 8.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xzzz[i] = -2.0 * tr_yyzz_xzzz[i] * tbe_0 - 6.0 * tr_yyzz_xzzz[i] * tke_0 + 4.0 * tr_yyzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_zzz[i] * tbe_0 + 8.0 * tr_xyyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_yyyy[i] = -2.0 * tr_yyzz_yyyy[i] * tbe_0 - 2.0 * tr_yyzz_yyyy[i] * tke_0 + 4.0 * tr_yyzz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_yyyz[i] = -2.0 * tr_yyzz_yyyz[i] * tbe_0 - 2.0 * tr_yyzz_yyyz[i] * tke_0 + 4.0 * tr_yyzz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_yyzz[i] = -2.0 * tr_yyzz_yyzz[i] * tbe_0 - 2.0 * tr_yyzz_yyzz[i] * tke_0 + 4.0 * tr_yyzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_yzzz[i] = -2.0 * tr_yyzz_yzzz[i] * tbe_0 - 2.0 * tr_yyzz_yzzz[i] * tke_0 + 4.0 * tr_yyzz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_zzzz[i] = -2.0 * tr_yyzz_zzzz[i] * tbe_0 - 2.0 * tr_yyzz_zzzz[i] * tke_0 + 4.0 * tr_yyzz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 195-210 components of targeted buffer : GG

    auto tr_0_0_xx_yzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 195);

    auto tr_0_0_xx_yzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 196);

    auto tr_0_0_xx_yzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 197);

    auto tr_0_0_xx_yzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 198);

    auto tr_0_0_xx_yzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 199);

    auto tr_0_0_xx_yzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 200);

    auto tr_0_0_xx_yzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 201);

    auto tr_0_0_xx_yzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 202);

    auto tr_0_0_xx_yzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 203);

    auto tr_0_0_xx_yzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 204);

    auto tr_0_0_xx_yzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 205);

    auto tr_0_0_xx_yzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 206);

    auto tr_0_0_xx_yzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 207);

    auto tr_0_0_xx_yzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 208);

    auto tr_0_0_xx_yzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 209);

    #pragma omp simd aligned(tr_0_0_xx_yzzz_xxxx, tr_0_0_xx_yzzz_xxxy, tr_0_0_xx_yzzz_xxxz, tr_0_0_xx_yzzz_xxyy, tr_0_0_xx_yzzz_xxyz, tr_0_0_xx_yzzz_xxzz, tr_0_0_xx_yzzz_xyyy, tr_0_0_xx_yzzz_xyyz, tr_0_0_xx_yzzz_xyzz, tr_0_0_xx_yzzz_xzzz, tr_0_0_xx_yzzz_yyyy, tr_0_0_xx_yzzz_yyyz, tr_0_0_xx_yzzz_yyzz, tr_0_0_xx_yzzz_yzzz, tr_0_0_xx_yzzz_zzzz, tr_xxyzzz_xxxx, tr_xxyzzz_xxxy, tr_xxyzzz_xxxz, tr_xxyzzz_xxyy, tr_xxyzzz_xxyz, tr_xxyzzz_xxzz, tr_xxyzzz_xyyy, tr_xxyzzz_xyyz, tr_xxyzzz_xyzz, tr_xxyzzz_xzzz, tr_xxyzzz_yyyy, tr_xxyzzz_yyyz, tr_xxyzzz_yyzz, tr_xxyzzz_yzzz, tr_xxyzzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxx, tr_xyzzz_xxxxy, tr_xyzzz_xxxxz, tr_xyzzz_xxxyy, tr_xyzzz_xxxyz, tr_xyzzz_xxxzz, tr_xyzzz_xxy, tr_xyzzz_xxyyy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xxzzz, tr_xyzzz_xyy, tr_xyzzz_xyyyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_xzzzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxxxx, tr_yzzz_xxxxxy, tr_yzzz_xxxxxz, tr_yzzz_xxxxyy, tr_yzzz_xxxxyz, tr_yzzz_xxxxzz, tr_yzzz_xxxy, tr_yzzz_xxxyyy, tr_yzzz_xxxyyz, tr_yzzz_xxxyzz, tr_yzzz_xxxz, tr_yzzz_xxxzzz, tr_yzzz_xxyy, tr_yzzz_xxyyyy, tr_yzzz_xxyyyz, tr_yzzz_xxyyzz, tr_yzzz_xxyz, tr_yzzz_xxyzzz, tr_yzzz_xxzz, tr_yzzz_xxzzzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yzzz_xxxx[i] = 12.0 * tr_yzzz_xx[i] - 2.0 * tr_yzzz_xxxx[i] * tbe_0 - 18.0 * tr_yzzz_xxxx[i] * tke_0 + 4.0 * tr_yzzz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xyzzz_xxx[i] * tbe_0 + 8.0 * tr_xyzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xxxy[i] = 6.0 * tr_yzzz_xy[i] - 2.0 * tr_yzzz_xxxy[i] * tbe_0 - 14.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xyzzz_xxy[i] * tbe_0 + 8.0 * tr_xyzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xxxz[i] = 6.0 * tr_yzzz_xz[i] - 2.0 * tr_yzzz_xxxz[i] * tbe_0 - 14.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xyzzz_xxz[i] * tbe_0 + 8.0 * tr_xyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xxyy[i] = 2.0 * tr_yzzz_yy[i] - 2.0 * tr_yzzz_xxyy[i] * tbe_0 - 10.0 * tr_yzzz_xxyy[i] * tke_0 + 4.0 * tr_yzzz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xyy[i] * tbe_0 + 8.0 * tr_xyzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xxyz[i] = 2.0 * tr_yzzz_yz[i] - 2.0 * tr_yzzz_xxyz[i] * tbe_0 - 10.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xyz[i] * tbe_0 + 8.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xxzz[i] = 2.0 * tr_yzzz_zz[i] - 2.0 * tr_yzzz_xxzz[i] * tbe_0 - 10.0 * tr_yzzz_xxzz[i] * tke_0 + 4.0 * tr_yzzz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xzz[i] * tbe_0 + 8.0 * tr_xyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xyyy[i] = -2.0 * tr_yzzz_xyyy[i] * tbe_0 - 6.0 * tr_yzzz_xyyy[i] * tke_0 + 4.0 * tr_yzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_yyy[i] * tbe_0 + 8.0 * tr_xyzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xyyz[i] = -2.0 * tr_yzzz_xyyz[i] * tbe_0 - 6.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_yyz[i] * tbe_0 + 8.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xyzz[i] = -2.0 * tr_yzzz_xyzz[i] * tbe_0 - 6.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_yzz[i] * tbe_0 + 8.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xzzz[i] = -2.0 * tr_yzzz_xzzz[i] * tbe_0 - 6.0 * tr_yzzz_xzzz[i] * tke_0 + 4.0 * tr_yzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_zzz[i] * tbe_0 + 8.0 * tr_xyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_yyyy[i] = -2.0 * tr_yzzz_yyyy[i] * tbe_0 - 2.0 * tr_yzzz_yyyy[i] * tke_0 + 4.0 * tr_yzzz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_yyyz[i] = -2.0 * tr_yzzz_yyyz[i] * tbe_0 - 2.0 * tr_yzzz_yyyz[i] * tke_0 + 4.0 * tr_yzzz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_yyzz[i] = -2.0 * tr_yzzz_yyzz[i] * tbe_0 - 2.0 * tr_yzzz_yyzz[i] * tke_0 + 4.0 * tr_yzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_yzzz[i] = -2.0 * tr_yzzz_yzzz[i] * tbe_0 - 2.0 * tr_yzzz_yzzz[i] * tke_0 + 4.0 * tr_yzzz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_zzzz[i] = -2.0 * tr_yzzz_zzzz[i] * tbe_0 - 2.0 * tr_yzzz_zzzz[i] * tke_0 + 4.0 * tr_yzzz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-225 components of targeted buffer : GG

    auto tr_0_0_xx_zzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 210);

    auto tr_0_0_xx_zzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 211);

    auto tr_0_0_xx_zzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 212);

    auto tr_0_0_xx_zzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 213);

    auto tr_0_0_xx_zzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 214);

    auto tr_0_0_xx_zzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 215);

    auto tr_0_0_xx_zzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 216);

    auto tr_0_0_xx_zzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 217);

    auto tr_0_0_xx_zzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 218);

    auto tr_0_0_xx_zzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 219);

    auto tr_0_0_xx_zzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 220);

    auto tr_0_0_xx_zzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 221);

    auto tr_0_0_xx_zzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 222);

    auto tr_0_0_xx_zzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 223);

    auto tr_0_0_xx_zzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 224);

    #pragma omp simd aligned(tr_0_0_xx_zzzz_xxxx, tr_0_0_xx_zzzz_xxxy, tr_0_0_xx_zzzz_xxxz, tr_0_0_xx_zzzz_xxyy, tr_0_0_xx_zzzz_xxyz, tr_0_0_xx_zzzz_xxzz, tr_0_0_xx_zzzz_xyyy, tr_0_0_xx_zzzz_xyyz, tr_0_0_xx_zzzz_xyzz, tr_0_0_xx_zzzz_xzzz, tr_0_0_xx_zzzz_yyyy, tr_0_0_xx_zzzz_yyyz, tr_0_0_xx_zzzz_yyzz, tr_0_0_xx_zzzz_yzzz, tr_0_0_xx_zzzz_zzzz, tr_xxzzzz_xxxx, tr_xxzzzz_xxxy, tr_xxzzzz_xxxz, tr_xxzzzz_xxyy, tr_xxzzzz_xxyz, tr_xxzzzz_xxzz, tr_xxzzzz_xyyy, tr_xxzzzz_xyyz, tr_xxzzzz_xyzz, tr_xxzzzz_xzzz, tr_xxzzzz_yyyy, tr_xxzzzz_yyyz, tr_xxzzzz_yyzz, tr_xxzzzz_yzzz, tr_xxzzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxxxx, tr_xzzzz_xxxxy, tr_xzzzz_xxxxz, tr_xzzzz_xxxyy, tr_xzzzz_xxxyz, tr_xzzzz_xxxzz, tr_xzzzz_xxy, tr_xzzzz_xxyyy, tr_xzzzz_xxyyz, tr_xzzzz_xxyzz, tr_xzzzz_xxz, tr_xzzzz_xxzzz, tr_xzzzz_xyy, tr_xzzzz_xyyyy, tr_xzzzz_xyyyz, tr_xzzzz_xyyzz, tr_xzzzz_xyz, tr_xzzzz_xyzzz, tr_xzzzz_xzz, tr_xzzzz_xzzzz, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_zzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxxxx, tr_zzzz_xxxxxy, tr_zzzz_xxxxxz, tr_zzzz_xxxxyy, tr_zzzz_xxxxyz, tr_zzzz_xxxxzz, tr_zzzz_xxxy, tr_zzzz_xxxyyy, tr_zzzz_xxxyyz, tr_zzzz_xxxyzz, tr_zzzz_xxxz, tr_zzzz_xxxzzz, tr_zzzz_xxyy, tr_zzzz_xxyyyy, tr_zzzz_xxyyyz, tr_zzzz_xxyyzz, tr_zzzz_xxyz, tr_zzzz_xxyzzz, tr_zzzz_xxzz, tr_zzzz_xxzzzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_zzzz_xxxx[i] = 12.0 * tr_zzzz_xx[i] - 2.0 * tr_zzzz_xxxx[i] * tbe_0 - 18.0 * tr_zzzz_xxxx[i] * tke_0 + 4.0 * tr_zzzz_xxxxxx[i] * tke_0 * tke_0 - 16.0 * tr_xzzzz_xxx[i] * tbe_0 + 8.0 * tr_xzzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xxxy[i] = 6.0 * tr_zzzz_xy[i] - 2.0 * tr_zzzz_xxxy[i] * tbe_0 - 14.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxxxy[i] * tke_0 * tke_0 - 12.0 * tr_xzzzz_xxy[i] * tbe_0 + 8.0 * tr_xzzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xxxz[i] = 6.0 * tr_zzzz_xz[i] - 2.0 * tr_zzzz_xxxz[i] * tbe_0 - 14.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxxxz[i] * tke_0 * tke_0 - 12.0 * tr_xzzzz_xxz[i] * tbe_0 + 8.0 * tr_xzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xxyy[i] = 2.0 * tr_zzzz_yy[i] - 2.0 * tr_zzzz_xxyy[i] * tbe_0 - 10.0 * tr_zzzz_xxyy[i] * tke_0 + 4.0 * tr_zzzz_xxxxyy[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_xyy[i] * tbe_0 + 8.0 * tr_xzzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xxyz[i] = 2.0 * tr_zzzz_yz[i] - 2.0 * tr_zzzz_xxyz[i] * tbe_0 - 10.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxxxyz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_xyz[i] * tbe_0 + 8.0 * tr_xzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xxzz[i] = 2.0 * tr_zzzz_zz[i] - 2.0 * tr_zzzz_xxzz[i] * tbe_0 - 10.0 * tr_zzzz_xxzz[i] * tke_0 + 4.0 * tr_zzzz_xxxxzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_xzz[i] * tbe_0 + 8.0 * tr_xzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xyyy[i] = -2.0 * tr_zzzz_xyyy[i] * tbe_0 - 6.0 * tr_zzzz_xyyy[i] * tke_0 + 4.0 * tr_zzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_yyy[i] * tbe_0 + 8.0 * tr_xzzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xyyz[i] = -2.0 * tr_zzzz_xyyz[i] * tbe_0 - 6.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_yyz[i] * tbe_0 + 8.0 * tr_xzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xyzz[i] = -2.0 * tr_zzzz_xyzz[i] * tbe_0 - 6.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_yzz[i] * tbe_0 + 8.0 * tr_xzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xzzz[i] = -2.0 * tr_zzzz_xzzz[i] * tbe_0 - 6.0 * tr_zzzz_xzzz[i] * tke_0 + 4.0 * tr_zzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_zzz[i] * tbe_0 + 8.0 * tr_xzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_yyyy[i] = -2.0 * tr_zzzz_yyyy[i] * tbe_0 - 2.0 * tr_zzzz_yyyy[i] * tke_0 + 4.0 * tr_zzzz_xxyyyy[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_yyyz[i] = -2.0 * tr_zzzz_yyyz[i] * tbe_0 - 2.0 * tr_zzzz_yyyz[i] * tke_0 + 4.0 * tr_zzzz_xxyyyz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_yyzz[i] = -2.0 * tr_zzzz_yyzz[i] * tbe_0 - 2.0 * tr_zzzz_yyzz[i] * tke_0 + 4.0 * tr_zzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_yzzz[i] = -2.0 * tr_zzzz_yzzz[i] * tbe_0 - 2.0 * tr_zzzz_yzzz[i] * tke_0 + 4.0 * tr_zzzz_xxyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_zzzz[i] = -2.0 * tr_zzzz_zzzz[i] * tbe_0 - 2.0 * tr_zzzz_zzzz[i] * tke_0 + 4.0 * tr_zzzz_xxzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 225-240 components of targeted buffer : GG

    auto tr_0_0_xy_xxxx_xxxx = pbuffer.data(idx_op_geom_020_gg + 225);

    auto tr_0_0_xy_xxxx_xxxy = pbuffer.data(idx_op_geom_020_gg + 226);

    auto tr_0_0_xy_xxxx_xxxz = pbuffer.data(idx_op_geom_020_gg + 227);

    auto tr_0_0_xy_xxxx_xxyy = pbuffer.data(idx_op_geom_020_gg + 228);

    auto tr_0_0_xy_xxxx_xxyz = pbuffer.data(idx_op_geom_020_gg + 229);

    auto tr_0_0_xy_xxxx_xxzz = pbuffer.data(idx_op_geom_020_gg + 230);

    auto tr_0_0_xy_xxxx_xyyy = pbuffer.data(idx_op_geom_020_gg + 231);

    auto tr_0_0_xy_xxxx_xyyz = pbuffer.data(idx_op_geom_020_gg + 232);

    auto tr_0_0_xy_xxxx_xyzz = pbuffer.data(idx_op_geom_020_gg + 233);

    auto tr_0_0_xy_xxxx_xzzz = pbuffer.data(idx_op_geom_020_gg + 234);

    auto tr_0_0_xy_xxxx_yyyy = pbuffer.data(idx_op_geom_020_gg + 235);

    auto tr_0_0_xy_xxxx_yyyz = pbuffer.data(idx_op_geom_020_gg + 236);

    auto tr_0_0_xy_xxxx_yyzz = pbuffer.data(idx_op_geom_020_gg + 237);

    auto tr_0_0_xy_xxxx_yzzz = pbuffer.data(idx_op_geom_020_gg + 238);

    auto tr_0_0_xy_xxxx_zzzz = pbuffer.data(idx_op_geom_020_gg + 239);

    #pragma omp simd aligned(tr_0_0_xy_xxxx_xxxx, tr_0_0_xy_xxxx_xxxy, tr_0_0_xy_xxxx_xxxz, tr_0_0_xy_xxxx_xxyy, tr_0_0_xy_xxxx_xxyz, tr_0_0_xy_xxxx_xxzz, tr_0_0_xy_xxxx_xyyy, tr_0_0_xy_xxxx_xyyz, tr_0_0_xy_xxxx_xyzz, tr_0_0_xy_xxxx_xzzz, tr_0_0_xy_xxxx_yyyy, tr_0_0_xy_xxxx_yyyz, tr_0_0_xy_xxxx_yyzz, tr_0_0_xy_xxxx_yzzz, tr_0_0_xy_xxxx_zzzz, tr_xxx_xxx, tr_xxx_xxxxy, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyyyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxxxy, tr_xxxx_xxxxyy, tr_xxxx_xxxxyz, tr_xxxx_xxxy, tr_xxxx_xxxyyy, tr_xxxx_xxxyyz, tr_xxxx_xxxyzz, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyyyy, tr_xxxx_xxyyyz, tr_xxxx_xxyyzz, tr_xxxx_xxyz, tr_xxxx_xxyzzz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyyyy, tr_xxxx_xyyyyz, tr_xxxx_xyyyzz, tr_xxxx_xyyz, tr_xxxx_xyyzzz, tr_xxxx_xyzz, tr_xxxx_xyzzzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxxx_xxx, tr_xxxxx_xxxxy, tr_xxxxx_xxxyy, tr_xxxxx_xxxyz, tr_xxxxx_xxy, tr_xxxxx_xxyyy, tr_xxxxx_xxyyz, tr_xxxxx_xxyzz, tr_xxxxx_xxz, tr_xxxxx_xyy, tr_xxxxx_xyyyy, tr_xxxxx_xyyyz, tr_xxxxx_xyyzz, tr_xxxxx_xyz, tr_xxxxx_xyzzz, tr_xxxxx_xzz, tr_xxxxx_yyy, tr_xxxxx_yyyyy, tr_xxxxx_yyyyz, tr_xxxxx_yyyzz, tr_xxxxx_yyz, tr_xxxxx_yyzzz, tr_xxxxx_yzz, tr_xxxxx_yzzzz, tr_xxxxx_zzz, tr_xxxxxy_xxxx, tr_xxxxxy_xxxy, tr_xxxxxy_xxxz, tr_xxxxxy_xxyy, tr_xxxxxy_xxyz, tr_xxxxxy_xxzz, tr_xxxxxy_xyyy, tr_xxxxxy_xyyz, tr_xxxxxy_xyzz, tr_xxxxxy_xzzz, tr_xxxxxy_yyyy, tr_xxxxxy_yyyz, tr_xxxxxy_yyzz, tr_xxxxxy_yzzz, tr_xxxxxy_zzzz, tr_xxxxy_xxx, tr_xxxxy_xxxxx, tr_xxxxy_xxxxy, tr_xxxxy_xxxxz, tr_xxxxy_xxxyy, tr_xxxxy_xxxyz, tr_xxxxy_xxxzz, tr_xxxxy_xxy, tr_xxxxy_xxyyy, tr_xxxxy_xxyyz, tr_xxxxy_xxyzz, tr_xxxxy_xxz, tr_xxxxy_xxzzz, tr_xxxxy_xyy, tr_xxxxy_xyyyy, tr_xxxxy_xyyyz, tr_xxxxy_xyyzz, tr_xxxxy_xyz, tr_xxxxy_xyzzz, tr_xxxxy_xzz, tr_xxxxy_xzzzz, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_zzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxxx_xxxx[i] = -8.0 * tr_xxx_xxxxy[i] * tke_0 - 8.0 * tr_xxxy_xxxx[i] * tbe_0 - 8.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xxxy[i] = 4.0 * tr_xxx_xxx[i] - 8.0 * tr_xxx_xxxyy[i] * tke_0 - 8.0 * tr_xxxy_xxxy[i] * tbe_0 + 3.0 * tr_xxxx_xx[i] - 6.0 * tr_xxxx_xxyy[i] * tke_0 - 2.0 * tr_xxxx_xxxx[i] * tke_0 + 4.0 * tr_xxxx_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xxxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxxy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxx_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xxxz[i] = -8.0 * tr_xxx_xxxyz[i] * tke_0 - 8.0 * tr_xxxy_xxxz[i] * tbe_0 - 6.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xxyy[i] = 8.0 * tr_xxx_xxy[i] - 8.0 * tr_xxx_xxyyy[i] * tke_0 - 8.0 * tr_xxxy_xxyy[i] * tbe_0 + 4.0 * tr_xxxx_xy[i] - 4.0 * tr_xxxx_xyyy[i] * tke_0 - 4.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxxy_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxxx_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xxyz[i] = 4.0 * tr_xxx_xxz[i] - 8.0 * tr_xxx_xxyyz[i] * tke_0 - 8.0 * tr_xxxy_xxyz[i] * tbe_0 + 2.0 * tr_xxxx_xz[i] - 4.0 * tr_xxxx_xyyz[i] * tke_0 - 2.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxxy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xxzz[i] = -8.0 * tr_xxx_xxyzz[i] * tke_0 - 8.0 * tr_xxxy_xxzz[i] * tbe_0 - 4.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xyyy[i] = 12.0 * tr_xxx_xyy[i] - 8.0 * tr_xxx_xyyyy[i] * tke_0 - 8.0 * tr_xxxy_xyyy[i] * tbe_0 + 3.0 * tr_xxxx_yy[i] - 2.0 * tr_xxxx_yyyy[i] * tke_0 - 6.0 * tr_xxxx_xxyy[i] * tke_0 + 4.0 * tr_xxxx_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxxx_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xyyz[i] = 8.0 * tr_xxx_xyz[i] - 8.0 * tr_xxx_xyyyz[i] * tke_0 - 8.0 * tr_xxxy_xyyz[i] * tbe_0 + 2.0 * tr_xxxx_yz[i] - 2.0 * tr_xxxx_yyyz[i] * tke_0 - 4.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xyzz[i] = 4.0 * tr_xxx_xzz[i] - 8.0 * tr_xxx_xyyzz[i] * tke_0 - 8.0 * tr_xxxy_xyzz[i] * tbe_0 + tr_xxxx_zz[i] - 2.0 * tr_xxxx_yyzz[i] * tke_0 - 2.0 * tr_xxxx_xxzz[i] * tke_0 + 4.0 * tr_xxxx_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xzzz[i] = -8.0 * tr_xxx_xyzzz[i] * tke_0 - 8.0 * tr_xxxy_xzzz[i] * tbe_0 - 2.0 * tr_xxxx_yzzz[i] * tke_0 + 4.0 * tr_xxxx_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_yyyy[i] = 16.0 * tr_xxx_yyy[i] - 8.0 * tr_xxx_yyyyy[i] * tke_0 - 8.0 * tr_xxxy_yyyy[i] * tbe_0 - 8.0 * tr_xxxx_xyyy[i] * tke_0 + 4.0 * tr_xxxx_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xxxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxx_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_yyyz[i] = 12.0 * tr_xxx_yyz[i] - 8.0 * tr_xxx_yyyyz[i] * tke_0 - 8.0 * tr_xxxy_yyyz[i] * tbe_0 - 6.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxxx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_yyzz[i] = 8.0 * tr_xxx_yzz[i] - 8.0 * tr_xxx_yyyzz[i] * tke_0 - 8.0 * tr_xxxy_yyzz[i] * tbe_0 - 4.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxxx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_yzzz[i] = 4.0 * tr_xxx_zzz[i] - 8.0 * tr_xxx_yyzzz[i] * tke_0 - 8.0 * tr_xxxy_yzzz[i] * tbe_0 - 2.0 * tr_xxxx_xzzz[i] * tke_0 + 4.0 * tr_xxxx_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_zzzz[i] = -8.0 * tr_xxx_yzzzz[i] * tke_0 - 8.0 * tr_xxxy_zzzz[i] * tbe_0 + 4.0 * tr_xxxx_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-255 components of targeted buffer : GG

    auto tr_0_0_xy_xxxy_xxxx = pbuffer.data(idx_op_geom_020_gg + 240);

    auto tr_0_0_xy_xxxy_xxxy = pbuffer.data(idx_op_geom_020_gg + 241);

    auto tr_0_0_xy_xxxy_xxxz = pbuffer.data(idx_op_geom_020_gg + 242);

    auto tr_0_0_xy_xxxy_xxyy = pbuffer.data(idx_op_geom_020_gg + 243);

    auto tr_0_0_xy_xxxy_xxyz = pbuffer.data(idx_op_geom_020_gg + 244);

    auto tr_0_0_xy_xxxy_xxzz = pbuffer.data(idx_op_geom_020_gg + 245);

    auto tr_0_0_xy_xxxy_xyyy = pbuffer.data(idx_op_geom_020_gg + 246);

    auto tr_0_0_xy_xxxy_xyyz = pbuffer.data(idx_op_geom_020_gg + 247);

    auto tr_0_0_xy_xxxy_xyzz = pbuffer.data(idx_op_geom_020_gg + 248);

    auto tr_0_0_xy_xxxy_xzzz = pbuffer.data(idx_op_geom_020_gg + 249);

    auto tr_0_0_xy_xxxy_yyyy = pbuffer.data(idx_op_geom_020_gg + 250);

    auto tr_0_0_xy_xxxy_yyyz = pbuffer.data(idx_op_geom_020_gg + 251);

    auto tr_0_0_xy_xxxy_yyzz = pbuffer.data(idx_op_geom_020_gg + 252);

    auto tr_0_0_xy_xxxy_yzzz = pbuffer.data(idx_op_geom_020_gg + 253);

    auto tr_0_0_xy_xxxy_zzzz = pbuffer.data(idx_op_geom_020_gg + 254);

    #pragma omp simd aligned(tr_0_0_xy_xxxy_xxxx, tr_0_0_xy_xxxy_xxxy, tr_0_0_xy_xxxy_xxxz, tr_0_0_xy_xxxy_xxyy, tr_0_0_xy_xxxy_xxyz, tr_0_0_xy_xxxy_xxzz, tr_0_0_xy_xxxy_xyyy, tr_0_0_xy_xxxy_xyyz, tr_0_0_xy_xxxy_xyzz, tr_0_0_xy_xxxy_xzzz, tr_0_0_xy_xxxy_yyyy, tr_0_0_xy_xxxy_yyyz, tr_0_0_xy_xxxy_yyzz, tr_0_0_xy_xxxy_yzzz, tr_0_0_xy_xxxy_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxx_xxx, tr_xxx_xxxxx, tr_xxx_xxxxy, tr_xxx_xxxxz, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xzzz, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yzzz, tr_xxxx_zzzz, tr_xxxxy_xxx, tr_xxxxy_xxxxy, tr_xxxxy_xxxyy, tr_xxxxy_xxxyz, tr_xxxxy_xxy, tr_xxxxy_xxyyy, tr_xxxxy_xxyyz, tr_xxxxy_xxyzz, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyyyy, tr_xxxxy_xyyyz, tr_xxxxy_xyyzz, tr_xxxxy_xyz, tr_xxxxy_xyzzz, tr_xxxxy_xzz, tr_xxxxy_yyy, tr_xxxxy_yyyyy, tr_xxxxy_yyyyz, tr_xxxxy_yyyzz, tr_xxxxy_yyz, tr_xxxxy_yyzzz, tr_xxxxy_yzz, tr_xxxxy_yzzzz, tr_xxxxy_zzz, tr_xxxxyy_xxxx, tr_xxxxyy_xxxy, tr_xxxxyy_xxxz, tr_xxxxyy_xxyy, tr_xxxxyy_xxyz, tr_xxxxyy_xxzz, tr_xxxxyy_xyyy, tr_xxxxyy_xyyz, tr_xxxxyy_xyzz, tr_xxxxyy_xzzz, tr_xxxxyy_yyyy, tr_xxxxyy_yyyz, tr_xxxxyy_yyzz, tr_xxxxyy_yzzz, tr_xxxxyy_zzzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxxxy, tr_xxxy_xxxxyy, tr_xxxy_xxxxyz, tr_xxxy_xxxy, tr_xxxy_xxxyyy, tr_xxxy_xxxyyz, tr_xxxy_xxxyzz, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyyyy, tr_xxxy_xxyyyz, tr_xxxy_xxyyzz, tr_xxxy_xxyz, tr_xxxy_xxyzzz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyyyy, tr_xxxy_xyyyyz, tr_xxxy_xyyyzz, tr_xxxy_xyyz, tr_xxxy_xyyzzz, tr_xxxy_xyzz, tr_xxxy_xyzzzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxyy_xxx, tr_xxxyy_xxxxx, tr_xxxyy_xxxxy, tr_xxxyy_xxxxz, tr_xxxyy_xxxyy, tr_xxxyy_xxxyz, tr_xxxyy_xxxzz, tr_xxxyy_xxy, tr_xxxyy_xxyyy, tr_xxxyy_xxyyz, tr_xxxyy_xxyzz, tr_xxxyy_xxz, tr_xxxyy_xxzzz, tr_xxxyy_xyy, tr_xxxyy_xyyyy, tr_xxxyy_xyyyz, tr_xxxyy_xyyzz, tr_xxxyy_xyz, tr_xxxyy_xyzzz, tr_xxxyy_xzz, tr_xxxyy_xzzzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxxy_xxxx[i] = 3.0 * tr_xx_xxxx[i] - 6.0 * tr_xxy_xxxxy[i] * tke_0 - 6.0 * tr_xxyy_xxxx[i] * tbe_0 + 4.0 * tr_xxx_xxx[i] - 2.0 * tr_xxx_xxxxx[i] * tke_0 - 8.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xxxy[i] = 3.0 * tr_xx_xxxy[i] + 3.0 * tr_xxy_xxx[i] - 6.0 * tr_xxy_xxxyy[i] * tke_0 - 6.0 * tr_xxyy_xxxy[i] * tbe_0 + 3.0 * tr_xxx_xxy[i] - 2.0 * tr_xxx_xxxxy[i] * tke_0 + 3.0 * tr_xxxy_xx[i] - 6.0 * tr_xxxy_xxyy[i] * tke_0 - 2.0 * tr_xxxy_xxxx[i] * tke_0 + 4.0 * tr_xxxy_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xxxyy_xxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxxy[i] * tbe_0 - 2.0 * tr_xxxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xxxz[i] = 3.0 * tr_xx_xxxz[i] - 6.0 * tr_xxy_xxxyz[i] * tke_0 - 6.0 * tr_xxyy_xxxz[i] * tbe_0 + 3.0 * tr_xxx_xxz[i] - 2.0 * tr_xxx_xxxxz[i] * tke_0 - 6.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxyy_xxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xxyy[i] = 3.0 * tr_xx_xxyy[i] + 6.0 * tr_xxy_xxy[i] - 6.0 * tr_xxy_xxyyy[i] * tke_0 - 6.0 * tr_xxyy_xxyy[i] * tbe_0 + 2.0 * tr_xxx_xyy[i] - 2.0 * tr_xxx_xxxyy[i] * tke_0 + 4.0 * tr_xxxy_xy[i] - 4.0 * tr_xxxy_xyyy[i] * tke_0 - 4.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_xyy[i] * tbe_0 + 4.0 * tr_xxxyy_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxyy[i] * tbe_0 - 4.0 * tr_xxxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xxyz[i] = 3.0 * tr_xx_xxyz[i] + 3.0 * tr_xxy_xxz[i] - 6.0 * tr_xxy_xxyyz[i] * tke_0 - 6.0 * tr_xxyy_xxyz[i] * tbe_0 + 2.0 * tr_xxx_xyz[i] - 2.0 * tr_xxx_xxxyz[i] * tke_0 + 2.0 * tr_xxxy_xz[i] - 4.0 * tr_xxxy_xyyz[i] * tke_0 - 2.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_xyz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxyz[i] * tbe_0 - 2.0 * tr_xxxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xxzz[i] = 3.0 * tr_xx_xxzz[i] - 6.0 * tr_xxy_xxyzz[i] * tke_0 - 6.0 * tr_xxyy_xxzz[i] * tbe_0 + 2.0 * tr_xxx_xzz[i] - 2.0 * tr_xxx_xxxzz[i] * tke_0 - 4.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_xzz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxzz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xyyy[i] = 3.0 * tr_xx_xyyy[i] + 9.0 * tr_xxy_xyy[i] - 6.0 * tr_xxy_xyyyy[i] * tke_0 - 6.0 * tr_xxyy_xyyy[i] * tbe_0 + tr_xxx_yyy[i] - 2.0 * tr_xxx_xxyyy[i] * tke_0 + 3.0 * tr_xxxy_yy[i] - 2.0 * tr_xxxy_yyyy[i] * tke_0 - 6.0 * tr_xxxy_xxyy[i] * tke_0 + 4.0 * tr_xxxy_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxyy_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyyy[i] * tbe_0 - 6.0 * tr_xxxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xyyz[i] = 3.0 * tr_xx_xyyz[i] + 6.0 * tr_xxy_xyz[i] - 6.0 * tr_xxy_xyyyz[i] * tke_0 - 6.0 * tr_xxyy_xyyz[i] * tbe_0 + tr_xxx_yyz[i] - 2.0 * tr_xxx_xxyyz[i] * tke_0 + 2.0 * tr_xxxy_yz[i] - 2.0 * tr_xxxy_yyyz[i] * tke_0 - 4.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyy_yyz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyyz[i] * tbe_0 - 4.0 * tr_xxxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xyzz[i] = 3.0 * tr_xx_xyzz[i] + 3.0 * tr_xxy_xzz[i] - 6.0 * tr_xxy_xyyzz[i] * tke_0 - 6.0 * tr_xxyy_xyzz[i] * tbe_0 + tr_xxx_yzz[i] - 2.0 * tr_xxx_xxyzz[i] * tke_0 + tr_xxxy_zz[i] - 2.0 * tr_xxxy_yyzz[i] * tke_0 - 2.0 * tr_xxxy_xxzz[i] * tke_0 + 4.0 * tr_xxxy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyy_yzz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyzz[i] * tbe_0 - 2.0 * tr_xxxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xzzz[i] = 3.0 * tr_xx_xzzz[i] - 6.0 * tr_xxy_xyzzz[i] * tke_0 - 6.0 * tr_xxyy_xzzz[i] * tbe_0 + tr_xxx_zzz[i] - 2.0 * tr_xxx_xxzzz[i] * tke_0 - 2.0 * tr_xxxy_yzzz[i] * tke_0 + 4.0 * tr_xxxy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xzzz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_yyyy[i] = 3.0 * tr_xx_yyyy[i] + 12.0 * tr_xxy_yyy[i] - 6.0 * tr_xxy_yyyyy[i] * tke_0 - 6.0 * tr_xxyy_yyyy[i] * tbe_0 - 2.0 * tr_xxx_xyyyy[i] * tke_0 - 8.0 * tr_xxxy_xyyy[i] * tke_0 + 4.0 * tr_xxxy_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyyy[i] * tbe_0 - 8.0 * tr_xxxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_yyyz[i] = 3.0 * tr_xx_yyyz[i] + 9.0 * tr_xxy_yyz[i] - 6.0 * tr_xxy_yyyyz[i] * tke_0 - 6.0 * tr_xxyy_yyyz[i] * tbe_0 - 2.0 * tr_xxx_xyyyz[i] * tke_0 - 6.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyyz[i] * tbe_0 - 6.0 * tr_xxxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_yyzz[i] = 3.0 * tr_xx_yyzz[i] + 6.0 * tr_xxy_yzz[i] - 6.0 * tr_xxy_yyyzz[i] * tke_0 - 6.0 * tr_xxyy_yyzz[i] * tbe_0 - 2.0 * tr_xxx_xyyzz[i] * tke_0 - 4.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyzz[i] * tbe_0 - 4.0 * tr_xxxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_yzzz[i] = 3.0 * tr_xx_yzzz[i] + 3.0 * tr_xxy_zzz[i] - 6.0 * tr_xxy_yyzzz[i] * tke_0 - 6.0 * tr_xxyy_yzzz[i] * tbe_0 - 2.0 * tr_xxx_xyzzz[i] * tke_0 - 2.0 * tr_xxxy_xzzz[i] * tke_0 + 4.0 * tr_xxxy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yzzz[i] * tbe_0 - 2.0 * tr_xxxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_zzzz[i] = 3.0 * tr_xx_zzzz[i] - 6.0 * tr_xxy_yzzzz[i] * tke_0 - 6.0 * tr_xxyy_zzzz[i] * tbe_0 - 2.0 * tr_xxx_xzzzz[i] * tke_0 + 4.0 * tr_xxxy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zzzz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 255-270 components of targeted buffer : GG

    auto tr_0_0_xy_xxxz_xxxx = pbuffer.data(idx_op_geom_020_gg + 255);

    auto tr_0_0_xy_xxxz_xxxy = pbuffer.data(idx_op_geom_020_gg + 256);

    auto tr_0_0_xy_xxxz_xxxz = pbuffer.data(idx_op_geom_020_gg + 257);

    auto tr_0_0_xy_xxxz_xxyy = pbuffer.data(idx_op_geom_020_gg + 258);

    auto tr_0_0_xy_xxxz_xxyz = pbuffer.data(idx_op_geom_020_gg + 259);

    auto tr_0_0_xy_xxxz_xxzz = pbuffer.data(idx_op_geom_020_gg + 260);

    auto tr_0_0_xy_xxxz_xyyy = pbuffer.data(idx_op_geom_020_gg + 261);

    auto tr_0_0_xy_xxxz_xyyz = pbuffer.data(idx_op_geom_020_gg + 262);

    auto tr_0_0_xy_xxxz_xyzz = pbuffer.data(idx_op_geom_020_gg + 263);

    auto tr_0_0_xy_xxxz_xzzz = pbuffer.data(idx_op_geom_020_gg + 264);

    auto tr_0_0_xy_xxxz_yyyy = pbuffer.data(idx_op_geom_020_gg + 265);

    auto tr_0_0_xy_xxxz_yyyz = pbuffer.data(idx_op_geom_020_gg + 266);

    auto tr_0_0_xy_xxxz_yyzz = pbuffer.data(idx_op_geom_020_gg + 267);

    auto tr_0_0_xy_xxxz_yzzz = pbuffer.data(idx_op_geom_020_gg + 268);

    auto tr_0_0_xy_xxxz_zzzz = pbuffer.data(idx_op_geom_020_gg + 269);

    #pragma omp simd aligned(tr_0_0_xy_xxxz_xxxx, tr_0_0_xy_xxxz_xxxy, tr_0_0_xy_xxxz_xxxz, tr_0_0_xy_xxxz_xxyy, tr_0_0_xy_xxxz_xxyz, tr_0_0_xy_xxxz_xxzz, tr_0_0_xy_xxxz_xyyy, tr_0_0_xy_xxxz_xyyz, tr_0_0_xy_xxxz_xyzz, tr_0_0_xy_xxxz_xzzz, tr_0_0_xy_xxxz_yyyy, tr_0_0_xy_xxxz_yyyz, tr_0_0_xy_xxxz_yyzz, tr_0_0_xy_xxxz_yzzz, tr_0_0_xy_xxxz_zzzz, tr_xxxxyz_xxxx, tr_xxxxyz_xxxy, tr_xxxxyz_xxxz, tr_xxxxyz_xxyy, tr_xxxxyz_xxyz, tr_xxxxyz_xxzz, tr_xxxxyz_xyyy, tr_xxxxyz_xyyz, tr_xxxxyz_xyzz, tr_xxxxyz_xzzz, tr_xxxxyz_yyyy, tr_xxxxyz_yyyz, tr_xxxxyz_yyzz, tr_xxxxyz_yzzz, tr_xxxxyz_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxxxy, tr_xxxxz_xxxyy, tr_xxxxz_xxxyz, tr_xxxxz_xxy, tr_xxxxz_xxyyy, tr_xxxxz_xxyyz, tr_xxxxz_xxyzz, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyyyy, tr_xxxxz_xyyyz, tr_xxxxz_xyyzz, tr_xxxxz_xyz, tr_xxxxz_xyzzz, tr_xxxxz_xzz, tr_xxxxz_yyy, tr_xxxxz_yyyyy, tr_xxxxz_yyyyz, tr_xxxxz_yyyzz, tr_xxxxz_yyz, tr_xxxxz_yyzzz, tr_xxxxz_yzz, tr_xxxxz_yzzzz, tr_xxxxz_zzz, tr_xxxyz_xxx, tr_xxxyz_xxxxx, tr_xxxyz_xxxxy, tr_xxxyz_xxxxz, tr_xxxyz_xxxyy, tr_xxxyz_xxxyz, tr_xxxyz_xxxzz, tr_xxxyz_xxy, tr_xxxyz_xxyyy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xxzzz, tr_xxxyz_xyy, tr_xxxyz_xyyyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_xzzzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxxxy, tr_xxxz_xxxxyy, tr_xxxz_xxxxyz, tr_xxxz_xxxy, tr_xxxz_xxxyyy, tr_xxxz_xxxyyz, tr_xxxz_xxxyzz, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyyyy, tr_xxxz_xxyyyz, tr_xxxz_xxyyzz, tr_xxxz_xxyz, tr_xxxz_xxyzzz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyyyy, tr_xxxz_xyyyyz, tr_xxxz_xyyyzz, tr_xxxz_xyyz, tr_xxxz_xyyzzz, tr_xxxz_xyzz, tr_xxxz_xyzzzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxxz_xxxx[i] = -6.0 * tr_xxz_xxxxy[i] * tke_0 - 6.0 * tr_xxyz_xxxx[i] * tbe_0 - 8.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xxxy[i] = 3.0 * tr_xxz_xxx[i] - 6.0 * tr_xxz_xxxyy[i] * tke_0 - 6.0 * tr_xxyz_xxxy[i] * tbe_0 + 3.0 * tr_xxxz_xx[i] - 6.0 * tr_xxxz_xxyy[i] * tke_0 - 2.0 * tr_xxxz_xxxx[i] * tke_0 + 4.0 * tr_xxxz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xxxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xxxz[i] = -6.0 * tr_xxz_xxxyz[i] * tke_0 - 6.0 * tr_xxyz_xxxz[i] * tbe_0 - 6.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xxyy[i] = 6.0 * tr_xxz_xxy[i] - 6.0 * tr_xxz_xxyyy[i] * tke_0 - 6.0 * tr_xxyz_xxyy[i] * tbe_0 + 4.0 * tr_xxxz_xy[i] - 4.0 * tr_xxxz_xyyy[i] * tke_0 - 4.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xxyz[i] = 3.0 * tr_xxz_xxz[i] - 6.0 * tr_xxz_xxyyz[i] * tke_0 - 6.0 * tr_xxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxxz_xz[i] - 4.0 * tr_xxxz_xyyz[i] * tke_0 - 2.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xyz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xxzz[i] = -6.0 * tr_xxz_xxyzz[i] * tke_0 - 6.0 * tr_xxyz_xxzz[i] * tbe_0 - 4.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xyyy[i] = 9.0 * tr_xxz_xyy[i] - 6.0 * tr_xxz_xyyyy[i] * tke_0 - 6.0 * tr_xxyz_xyyy[i] * tbe_0 + 3.0 * tr_xxxz_yy[i] - 2.0 * tr_xxxz_yyyy[i] * tke_0 - 6.0 * tr_xxxz_xxyy[i] * tke_0 + 4.0 * tr_xxxz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xyyz[i] = 6.0 * tr_xxz_xyz[i] - 6.0 * tr_xxz_xyyyz[i] * tke_0 - 6.0 * tr_xxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxxz_yz[i] - 2.0 * tr_xxxz_yyyz[i] * tke_0 - 4.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xyzz[i] = 3.0 * tr_xxz_xzz[i] - 6.0 * tr_xxz_xyyzz[i] * tke_0 - 6.0 * tr_xxyz_xyzz[i] * tbe_0 + tr_xxxz_zz[i] - 2.0 * tr_xxxz_yyzz[i] * tke_0 - 2.0 * tr_xxxz_xxzz[i] * tke_0 + 4.0 * tr_xxxz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xzzz[i] = -6.0 * tr_xxz_xyzzz[i] * tke_0 - 6.0 * tr_xxyz_xzzz[i] * tbe_0 - 2.0 * tr_xxxz_yzzz[i] * tke_0 + 4.0 * tr_xxxz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_yyyy[i] = 12.0 * tr_xxz_yyy[i] - 6.0 * tr_xxz_yyyyy[i] * tke_0 - 6.0 * tr_xxyz_yyyy[i] * tbe_0 - 8.0 * tr_xxxz_xyyy[i] * tke_0 + 4.0 * tr_xxxz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xxxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_yyyz[i] = 9.0 * tr_xxz_yyz[i] - 6.0 * tr_xxz_yyyyz[i] * tke_0 - 6.0 * tr_xxyz_yyyz[i] * tbe_0 - 6.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_yyzz[i] = 6.0 * tr_xxz_yzz[i] - 6.0 * tr_xxz_yyyzz[i] * tke_0 - 6.0 * tr_xxyz_yyzz[i] * tbe_0 - 4.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_yzzz[i] = 3.0 * tr_xxz_zzz[i] - 6.0 * tr_xxz_yyzzz[i] * tke_0 - 6.0 * tr_xxyz_yzzz[i] * tbe_0 - 2.0 * tr_xxxz_xzzz[i] * tke_0 + 4.0 * tr_xxxz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_zzzz[i] = -6.0 * tr_xxz_yzzzz[i] * tke_0 - 6.0 * tr_xxyz_zzzz[i] * tbe_0 + 4.0 * tr_xxxz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-285 components of targeted buffer : GG

    auto tr_0_0_xy_xxyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 270);

    auto tr_0_0_xy_xxyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 271);

    auto tr_0_0_xy_xxyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 272);

    auto tr_0_0_xy_xxyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 273);

    auto tr_0_0_xy_xxyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 274);

    auto tr_0_0_xy_xxyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 275);

    auto tr_0_0_xy_xxyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 276);

    auto tr_0_0_xy_xxyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 277);

    auto tr_0_0_xy_xxyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 278);

    auto tr_0_0_xy_xxyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 279);

    auto tr_0_0_xy_xxyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 280);

    auto tr_0_0_xy_xxyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 281);

    auto tr_0_0_xy_xxyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 282);

    auto tr_0_0_xy_xxyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 283);

    auto tr_0_0_xy_xxyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 284);

    #pragma omp simd aligned(tr_0_0_xy_xxyy_xxxx, tr_0_0_xy_xxyy_xxxy, tr_0_0_xy_xxyy_xxxz, tr_0_0_xy_xxyy_xxyy, tr_0_0_xy_xxyy_xxyz, tr_0_0_xy_xxyy_xxzz, tr_0_0_xy_xxyy_xyyy, tr_0_0_xy_xxyy_xyyz, tr_0_0_xy_xxyy_xyzz, tr_0_0_xy_xxyy_xzzz, tr_0_0_xy_xxyy_yyyy, tr_0_0_xy_xxyy_yyyz, tr_0_0_xy_xxyy_yyzz, tr_0_0_xy_xxyy_yzzz, tr_0_0_xy_xxyy_zzzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xxxyy_xxx, tr_xxxyy_xxxxy, tr_xxxyy_xxxyy, tr_xxxyy_xxxyz, tr_xxxyy_xxy, tr_xxxyy_xxyyy, tr_xxxyy_xxyyz, tr_xxxyy_xxyzz, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyyyy, tr_xxxyy_xyyyz, tr_xxxyy_xyyzz, tr_xxxyy_xyz, tr_xxxyy_xyzzz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyyyy, tr_xxxyy_yyyyz, tr_xxxyy_yyyzz, tr_xxxyy_yyz, tr_xxxyy_yyzzz, tr_xxxyy_yzz, tr_xxxyy_yzzzz, tr_xxxyy_zzz, tr_xxxyyy_xxxx, tr_xxxyyy_xxxy, tr_xxxyyy_xxxz, tr_xxxyyy_xxyy, tr_xxxyyy_xxyz, tr_xxxyyy_xxzz, tr_xxxyyy_xyyy, tr_xxxyyy_xyyz, tr_xxxyyy_xyzz, tr_xxxyyy_xzzz, tr_xxxyyy_yyyy, tr_xxxyyy_yyyz, tr_xxxyyy_yyzz, tr_xxxyyy_yzzz, tr_xxxyyy_zzzz, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxxxy, tr_xxyy_xxxxyy, tr_xxyy_xxxxyz, tr_xxyy_xxxy, tr_xxyy_xxxyyy, tr_xxyy_xxxyyz, tr_xxyy_xxxyzz, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyyyy, tr_xxyy_xxyyyz, tr_xxyy_xxyyzz, tr_xxyy_xxyz, tr_xxyy_xxyzzz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyyyy, tr_xxyy_xyyyyz, tr_xxyy_xyyyzz, tr_xxyy_xyyz, tr_xxyy_xyyzzz, tr_xxyy_xyzz, tr_xxyy_xyzzzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyyy_xxx, tr_xxyyy_xxxxx, tr_xxyyy_xxxxy, tr_xxyyy_xxxxz, tr_xxyyy_xxxyy, tr_xxyyy_xxxyz, tr_xxyyy_xxxzz, tr_xxyyy_xxy, tr_xxyyy_xxyyy, tr_xxyyy_xxyyz, tr_xxyyy_xxyzz, tr_xxyyy_xxz, tr_xxyyy_xxzzz, tr_xxyyy_xyy, tr_xxyyy_xyyyy, tr_xxyyy_xyyyz, tr_xxyyy_xyyzz, tr_xxyyy_xyz, tr_xxyyy_xyzzz, tr_xxyyy_xzz, tr_xxyyy_xzzzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxyy_xxxx[i] = 4.0 * tr_xy_xxxx[i] - 4.0 * tr_xyy_xxxxy[i] * tke_0 - 4.0 * tr_xyyy_xxxx[i] * tbe_0 + 8.0 * tr_xxy_xxx[i] - 4.0 * tr_xxy_xxxxx[i] * tke_0 - 8.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xxxy[i] = 4.0 * tr_xy_xxxy[i] + 2.0 * tr_xyy_xxx[i] - 4.0 * tr_xyy_xxxyy[i] * tke_0 - 4.0 * tr_xyyy_xxxy[i] * tbe_0 + 6.0 * tr_xxy_xxy[i] - 4.0 * tr_xxy_xxxxy[i] * tke_0 + 3.0 * tr_xxyy_xx[i] - 6.0 * tr_xxyy_xxyy[i] * tke_0 - 2.0 * tr_xxyy_xxxx[i] * tke_0 + 4.0 * tr_xxyy_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xxyyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxxy[i] * tbe_0 - 2.0 * tr_xxxyy_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xxxz[i] = 4.0 * tr_xy_xxxz[i] - 4.0 * tr_xyy_xxxyz[i] * tke_0 - 4.0 * tr_xyyy_xxxz[i] * tbe_0 + 6.0 * tr_xxy_xxz[i] - 4.0 * tr_xxy_xxxxz[i] * tke_0 - 6.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xxyy[i] = 4.0 * tr_xy_xxyy[i] + 4.0 * tr_xyy_xxy[i] - 4.0 * tr_xyy_xxyyy[i] * tke_0 - 4.0 * tr_xyyy_xxyy[i] * tbe_0 + 4.0 * tr_xxy_xyy[i] - 4.0 * tr_xxy_xxxyy[i] * tke_0 + 4.0 * tr_xxyy_xy[i] - 4.0 * tr_xxyy_xyyy[i] * tke_0 - 4.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyy_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxyy[i] * tbe_0 - 4.0 * tr_xxxyy_xxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xxyz[i] = 4.0 * tr_xy_xxyz[i] + 2.0 * tr_xyy_xxz[i] - 4.0 * tr_xyy_xxyyz[i] * tke_0 - 4.0 * tr_xyyy_xxyz[i] * tbe_0 + 4.0 * tr_xxy_xyz[i] - 4.0 * tr_xxy_xxxyz[i] * tke_0 + 2.0 * tr_xxyy_xz[i] - 4.0 * tr_xxyy_xyyz[i] * tke_0 - 2.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxyz[i] * tbe_0 - 2.0 * tr_xxxyy_xxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xxzz[i] = 4.0 * tr_xy_xxzz[i] - 4.0 * tr_xyy_xxyzz[i] * tke_0 - 4.0 * tr_xyyy_xxzz[i] * tbe_0 + 4.0 * tr_xxy_xzz[i] - 4.0 * tr_xxy_xxxzz[i] * tke_0 - 4.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxzz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xyyy[i] = 4.0 * tr_xy_xyyy[i] + 6.0 * tr_xyy_xyy[i] - 4.0 * tr_xyy_xyyyy[i] * tke_0 - 4.0 * tr_xyyy_xyyy[i] * tbe_0 + 2.0 * tr_xxy_yyy[i] - 4.0 * tr_xxy_xxyyy[i] * tke_0 + 3.0 * tr_xxyy_yy[i] - 2.0 * tr_xxyy_yyyy[i] * tke_0 - 6.0 * tr_xxyy_xxyy[i] * tke_0 + 4.0 * tr_xxyy_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xyyy[i] * tbe_0 - 6.0 * tr_xxxyy_xyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xyyz[i] = 4.0 * tr_xy_xyyz[i] + 4.0 * tr_xyy_xyz[i] - 4.0 * tr_xyy_xyyyz[i] * tke_0 - 4.0 * tr_xyyy_xyyz[i] * tbe_0 + 2.0 * tr_xxy_yyz[i] - 4.0 * tr_xxy_xxyyz[i] * tke_0 + 2.0 * tr_xxyy_yz[i] - 2.0 * tr_xxyy_yyyz[i] * tke_0 - 4.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xyyz[i] * tbe_0 - 4.0 * tr_xxxyy_xyz[i] * tbe_0 + 4.0 * tr_xxxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xyzz[i] = 4.0 * tr_xy_xyzz[i] + 2.0 * tr_xyy_xzz[i] - 4.0 * tr_xyy_xyyzz[i] * tke_0 - 4.0 * tr_xyyy_xyzz[i] * tbe_0 + 2.0 * tr_xxy_yzz[i] - 4.0 * tr_xxy_xxyzz[i] * tke_0 + tr_xxyy_zz[i] - 2.0 * tr_xxyy_yyzz[i] * tke_0 - 2.0 * tr_xxyy_xxzz[i] * tke_0 + 4.0 * tr_xxyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xyzz[i] * tbe_0 - 2.0 * tr_xxxyy_xzz[i] * tbe_0 + 4.0 * tr_xxxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xzzz[i] = 4.0 * tr_xy_xzzz[i] - 4.0 * tr_xyy_xyzzz[i] * tke_0 - 4.0 * tr_xyyy_xzzz[i] * tbe_0 + 2.0 * tr_xxy_zzz[i] - 4.0 * tr_xxy_xxzzz[i] * tke_0 - 2.0 * tr_xxyy_yzzz[i] * tke_0 + 4.0 * tr_xxyy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xzzz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_yyyy[i] = 4.0 * tr_xy_yyyy[i] + 8.0 * tr_xyy_yyy[i] - 4.0 * tr_xyy_yyyyy[i] * tke_0 - 4.0 * tr_xyyy_yyyy[i] * tbe_0 - 4.0 * tr_xxy_xyyyy[i] * tke_0 - 8.0 * tr_xxyy_xyyy[i] * tke_0 + 4.0 * tr_xxyy_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yyyy[i] * tbe_0 - 8.0 * tr_xxxyy_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_yyyz[i] = 4.0 * tr_xy_yyyz[i] + 6.0 * tr_xyy_yyz[i] - 4.0 * tr_xyy_yyyyz[i] * tke_0 - 4.0 * tr_xyyy_yyyz[i] * tbe_0 - 4.0 * tr_xxy_xyyyz[i] * tke_0 - 6.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yyyz[i] * tbe_0 - 6.0 * tr_xxxyy_yyz[i] * tbe_0 + 4.0 * tr_xxxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_yyzz[i] = 4.0 * tr_xy_yyzz[i] + 4.0 * tr_xyy_yzz[i] - 4.0 * tr_xyy_yyyzz[i] * tke_0 - 4.0 * tr_xyyy_yyzz[i] * tbe_0 - 4.0 * tr_xxy_xyyzz[i] * tke_0 - 4.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yyzz[i] * tbe_0 - 4.0 * tr_xxxyy_yzz[i] * tbe_0 + 4.0 * tr_xxxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_yzzz[i] = 4.0 * tr_xy_yzzz[i] + 2.0 * tr_xyy_zzz[i] - 4.0 * tr_xyy_yyzzz[i] * tke_0 - 4.0 * tr_xyyy_yzzz[i] * tbe_0 - 4.0 * tr_xxy_xyzzz[i] * tke_0 - 2.0 * tr_xxyy_xzzz[i] * tke_0 + 4.0 * tr_xxyy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yzzz[i] * tbe_0 - 2.0 * tr_xxxyy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_zzzz[i] = 4.0 * tr_xy_zzzz[i] - 4.0 * tr_xyy_yzzzz[i] * tke_0 - 4.0 * tr_xyyy_zzzz[i] * tbe_0 - 4.0 * tr_xxy_xzzzz[i] * tke_0 + 4.0 * tr_xxyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_zzzz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 285-300 components of targeted buffer : GG

    auto tr_0_0_xy_xxyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 285);

    auto tr_0_0_xy_xxyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 286);

    auto tr_0_0_xy_xxyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 287);

    auto tr_0_0_xy_xxyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 288);

    auto tr_0_0_xy_xxyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 289);

    auto tr_0_0_xy_xxyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 290);

    auto tr_0_0_xy_xxyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 291);

    auto tr_0_0_xy_xxyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 292);

    auto tr_0_0_xy_xxyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 293);

    auto tr_0_0_xy_xxyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 294);

    auto tr_0_0_xy_xxyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 295);

    auto tr_0_0_xy_xxyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 296);

    auto tr_0_0_xy_xxyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 297);

    auto tr_0_0_xy_xxyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 298);

    auto tr_0_0_xy_xxyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 299);

    #pragma omp simd aligned(tr_0_0_xy_xxyz_xxxx, tr_0_0_xy_xxyz_xxxy, tr_0_0_xy_xxyz_xxxz, tr_0_0_xy_xxyz_xxyy, tr_0_0_xy_xxyz_xxyz, tr_0_0_xy_xxyz_xxzz, tr_0_0_xy_xxyz_xyyy, tr_0_0_xy_xxyz_xyyz, tr_0_0_xy_xxyz_xyzz, tr_0_0_xy_xxyz_xzzz, tr_0_0_xy_xxyz_yyyy, tr_0_0_xy_xxyz_yyyz, tr_0_0_xy_xxyz_yyzz, tr_0_0_xy_xxyz_yzzz, tr_0_0_xy_xxyz_zzzz, tr_xxxyyz_xxxx, tr_xxxyyz_xxxy, tr_xxxyyz_xxxz, tr_xxxyyz_xxyy, tr_xxxyyz_xxyz, tr_xxxyyz_xxzz, tr_xxxyyz_xyyy, tr_xxxyyz_xyyz, tr_xxxyyz_xyzz, tr_xxxyyz_xzzz, tr_xxxyyz_yyyy, tr_xxxyyz_yyyz, tr_xxxyyz_yyzz, tr_xxxyyz_yzzz, tr_xxxyyz_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxy, tr_xxxyz_xxxyy, tr_xxxyz_xxxyz, tr_xxxyz_xxy, tr_xxxyz_xxyyy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyyyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyyyy, tr_xxxyz_yyyyz, tr_xxxyz_yyyzz, tr_xxxyz_yyz, tr_xxxyz_yyzzz, tr_xxxyz_yzz, tr_xxxyz_yzzzz, tr_xxxyz_zzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxx, tr_xxyyz_xxxxy, tr_xxyyz_xxxxz, tr_xxyyz_xxxyy, tr_xxyyz_xxxyz, tr_xxyyz_xxxzz, tr_xxyyz_xxy, tr_xxyyz_xxyyy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xxzzz, tr_xxyyz_xyy, tr_xxyyz_xyyyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_xzzzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxxxy, tr_xxyz_xxxxyy, tr_xxyz_xxxxyz, tr_xxyz_xxxy, tr_xxyz_xxxyyy, tr_xxyz_xxxyyz, tr_xxyz_xxxyzz, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyyyy, tr_xxyz_xxyyyz, tr_xxyz_xxyyzz, tr_xxyz_xxyz, tr_xxyz_xxyzzz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyyyy, tr_xxyz_xyyyyz, tr_xxyz_xyyyzz, tr_xxyz_xyyz, tr_xxyz_xyyzzz, tr_xxyz_xyzz, tr_xxyz_xyzzzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxyz_xxxx[i] = 2.0 * tr_xz_xxxx[i] - 4.0 * tr_xyz_xxxxy[i] * tke_0 - 4.0 * tr_xyyz_xxxx[i] * tbe_0 + 4.0 * tr_xxz_xxx[i] - 2.0 * tr_xxz_xxxxx[i] * tke_0 - 8.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xxxy[i] = 2.0 * tr_xz_xxxy[i] + 2.0 * tr_xyz_xxx[i] - 4.0 * tr_xyz_xxxyy[i] * tke_0 - 4.0 * tr_xyyz_xxxy[i] * tbe_0 + 3.0 * tr_xxz_xxy[i] - 2.0 * tr_xxz_xxxxy[i] * tke_0 + 3.0 * tr_xxyz_xx[i] - 6.0 * tr_xxyz_xxyy[i] * tke_0 - 2.0 * tr_xxyz_xxxx[i] * tke_0 + 4.0 * tr_xxyz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xxyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxxy[i] * tbe_0 - 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xxxz[i] = 2.0 * tr_xz_xxxz[i] - 4.0 * tr_xyz_xxxyz[i] * tke_0 - 4.0 * tr_xyyz_xxxz[i] * tbe_0 + 3.0 * tr_xxz_xxz[i] - 2.0 * tr_xxz_xxxxz[i] * tke_0 - 6.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xxyy[i] = 2.0 * tr_xz_xxyy[i] + 4.0 * tr_xyz_xxy[i] - 4.0 * tr_xyz_xxyyy[i] * tke_0 - 4.0 * tr_xyyz_xxyy[i] * tbe_0 + 2.0 * tr_xxz_xyy[i] - 2.0 * tr_xxz_xxxyy[i] * tke_0 + 4.0 * tr_xxyz_xy[i] - 4.0 * tr_xxyz_xyyy[i] * tke_0 - 4.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxyy[i] * tbe_0 - 4.0 * tr_xxxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xxyz[i] = 2.0 * tr_xz_xxyz[i] + 2.0 * tr_xyz_xxz[i] - 4.0 * tr_xyz_xxyyz[i] * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tbe_0 + 2.0 * tr_xxz_xyz[i] - 2.0 * tr_xxz_xxxyz[i] * tke_0 + 2.0 * tr_xxyz_xz[i] - 4.0 * tr_xxyz_xyyz[i] * tke_0 - 2.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxyz[i] * tbe_0 - 2.0 * tr_xxxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xxzz[i] = 2.0 * tr_xz_xxzz[i] - 4.0 * tr_xyz_xxyzz[i] * tke_0 - 4.0 * tr_xyyz_xxzz[i] * tbe_0 + 2.0 * tr_xxz_xzz[i] - 2.0 * tr_xxz_xxxzz[i] * tke_0 - 4.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xyyy[i] = 2.0 * tr_xz_xyyy[i] + 6.0 * tr_xyz_xyy[i] - 4.0 * tr_xyz_xyyyy[i] * tke_0 - 4.0 * tr_xyyz_xyyy[i] * tbe_0 + tr_xxz_yyy[i] - 2.0 * tr_xxz_xxyyy[i] * tke_0 + 3.0 * tr_xxyz_yy[i] - 2.0 * tr_xxyz_yyyy[i] * tke_0 - 6.0 * tr_xxyz_xxyy[i] * tke_0 + 4.0 * tr_xxyz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xyyy[i] * tbe_0 - 6.0 * tr_xxxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xyyz[i] = 2.0 * tr_xz_xyyz[i] + 4.0 * tr_xyz_xyz[i] - 4.0 * tr_xyz_xyyyz[i] * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tbe_0 + tr_xxz_yyz[i] - 2.0 * tr_xxz_xxyyz[i] * tke_0 + 2.0 * tr_xxyz_yz[i] - 2.0 * tr_xxyz_yyyz[i] * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xyyz[i] * tbe_0 - 4.0 * tr_xxxyz_xyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xyzz[i] = 2.0 * tr_xz_xyzz[i] + 2.0 * tr_xyz_xzz[i] - 4.0 * tr_xyz_xyyzz[i] * tke_0 - 4.0 * tr_xyyz_xyzz[i] * tbe_0 + tr_xxz_yzz[i] - 2.0 * tr_xxz_xxyzz[i] * tke_0 + tr_xxyz_zz[i] - 2.0 * tr_xxyz_yyzz[i] * tke_0 - 2.0 * tr_xxyz_xxzz[i] * tke_0 + 4.0 * tr_xxyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xyzz[i] * tbe_0 - 2.0 * tr_xxxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xzzz[i] = 2.0 * tr_xz_xzzz[i] - 4.0 * tr_xyz_xyzzz[i] * tke_0 - 4.0 * tr_xyyz_xzzz[i] * tbe_0 + tr_xxz_zzz[i] - 2.0 * tr_xxz_xxzzz[i] * tke_0 - 2.0 * tr_xxyz_yzzz[i] * tke_0 + 4.0 * tr_xxyz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xzzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_yyyy[i] = 2.0 * tr_xz_yyyy[i] + 8.0 * tr_xyz_yyy[i] - 4.0 * tr_xyz_yyyyy[i] * tke_0 - 4.0 * tr_xyyz_yyyy[i] * tbe_0 - 2.0 * tr_xxz_xyyyy[i] * tke_0 - 8.0 * tr_xxyz_xyyy[i] * tke_0 + 4.0 * tr_xxyz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yyyy[i] * tbe_0 - 8.0 * tr_xxxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_yyyz[i] = 2.0 * tr_xz_yyyz[i] + 6.0 * tr_xyz_yyz[i] - 4.0 * tr_xyz_yyyyz[i] * tke_0 - 4.0 * tr_xyyz_yyyz[i] * tbe_0 - 2.0 * tr_xxz_xyyyz[i] * tke_0 - 6.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yyyz[i] * tbe_0 - 6.0 * tr_xxxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_yyzz[i] = 2.0 * tr_xz_yyzz[i] + 4.0 * tr_xyz_yzz[i] - 4.0 * tr_xyz_yyyzz[i] * tke_0 - 4.0 * tr_xyyz_yyzz[i] * tbe_0 - 2.0 * tr_xxz_xyyzz[i] * tke_0 - 4.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yyzz[i] * tbe_0 - 4.0 * tr_xxxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_yzzz[i] = 2.0 * tr_xz_yzzz[i] + 2.0 * tr_xyz_zzz[i] - 4.0 * tr_xyz_yyzzz[i] * tke_0 - 4.0 * tr_xyyz_yzzz[i] * tbe_0 - 2.0 * tr_xxz_xyzzz[i] * tke_0 - 2.0 * tr_xxyz_xzzz[i] * tke_0 + 4.0 * tr_xxyz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yzzz[i] * tbe_0 - 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_zzzz[i] = 2.0 * tr_xz_zzzz[i] - 4.0 * tr_xyz_yzzzz[i] * tke_0 - 4.0 * tr_xyyz_zzzz[i] * tbe_0 - 2.0 * tr_xxz_xzzzz[i] * tke_0 + 4.0 * tr_xxyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_zzzz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-315 components of targeted buffer : GG

    auto tr_0_0_xy_xxzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 300);

    auto tr_0_0_xy_xxzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 301);

    auto tr_0_0_xy_xxzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 302);

    auto tr_0_0_xy_xxzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 303);

    auto tr_0_0_xy_xxzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 304);

    auto tr_0_0_xy_xxzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 305);

    auto tr_0_0_xy_xxzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 306);

    auto tr_0_0_xy_xxzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 307);

    auto tr_0_0_xy_xxzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 308);

    auto tr_0_0_xy_xxzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 309);

    auto tr_0_0_xy_xxzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 310);

    auto tr_0_0_xy_xxzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 311);

    auto tr_0_0_xy_xxzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 312);

    auto tr_0_0_xy_xxzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 313);

    auto tr_0_0_xy_xxzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 314);

    #pragma omp simd aligned(tr_0_0_xy_xxzz_xxxx, tr_0_0_xy_xxzz_xxxy, tr_0_0_xy_xxzz_xxxz, tr_0_0_xy_xxzz_xxyy, tr_0_0_xy_xxzz_xxyz, tr_0_0_xy_xxzz_xxzz, tr_0_0_xy_xxzz_xyyy, tr_0_0_xy_xxzz_xyyz, tr_0_0_xy_xxzz_xyzz, tr_0_0_xy_xxzz_xzzz, tr_0_0_xy_xxzz_yyyy, tr_0_0_xy_xxzz_yyyz, tr_0_0_xy_xxzz_yyzz, tr_0_0_xy_xxzz_yzzz, tr_0_0_xy_xxzz_zzzz, tr_xxxyzz_xxxx, tr_xxxyzz_xxxy, tr_xxxyzz_xxxz, tr_xxxyzz_xxyy, tr_xxxyzz_xxyz, tr_xxxyzz_xxzz, tr_xxxyzz_xyyy, tr_xxxyzz_xyyz, tr_xxxyzz_xyzz, tr_xxxyzz_xzzz, tr_xxxyzz_yyyy, tr_xxxyzz_yyyz, tr_xxxyzz_yyzz, tr_xxxyzz_yzzz, tr_xxxyzz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxxxy, tr_xxxzz_xxxyy, tr_xxxzz_xxxyz, tr_xxxzz_xxy, tr_xxxzz_xxyyy, tr_xxxzz_xxyyz, tr_xxxzz_xxyzz, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyyyy, tr_xxxzz_xyyyz, tr_xxxzz_xyyzz, tr_xxxzz_xyz, tr_xxxzz_xyzzz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyyyy, tr_xxxzz_yyyyz, tr_xxxzz_yyyzz, tr_xxxzz_yyz, tr_xxxzz_yyzzz, tr_xxxzz_yzz, tr_xxxzz_yzzzz, tr_xxxzz_zzz, tr_xxyzz_xxx, tr_xxyzz_xxxxx, tr_xxyzz_xxxxy, tr_xxyzz_xxxxz, tr_xxyzz_xxxyy, tr_xxyzz_xxxyz, tr_xxyzz_xxxzz, tr_xxyzz_xxy, tr_xxyzz_xxyyy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xxzzz, tr_xxyzz_xyy, tr_xxyzz_xyyyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_xzzzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxxxy, tr_xxzz_xxxxyy, tr_xxzz_xxxxyz, tr_xxzz_xxxy, tr_xxzz_xxxyyy, tr_xxzz_xxxyyz, tr_xxzz_xxxyzz, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyyyy, tr_xxzz_xxyyyz, tr_xxzz_xxyyzz, tr_xxzz_xxyz, tr_xxzz_xxyzzz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyyyy, tr_xxzz_xyyyyz, tr_xxzz_xyyyzz, tr_xxzz_xyyz, tr_xxzz_xyyzzz, tr_xxzz_xyzz, tr_xxzz_xyzzzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxzz_xxxx[i] = -4.0 * tr_xzz_xxxxy[i] * tke_0 - 4.0 * tr_xyzz_xxxx[i] * tbe_0 - 8.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xxxy[i] = 2.0 * tr_xzz_xxx[i] - 4.0 * tr_xzz_xxxyy[i] * tke_0 - 4.0 * tr_xyzz_xxxy[i] * tbe_0 + 3.0 * tr_xxzz_xx[i] - 6.0 * tr_xxzz_xxyy[i] * tke_0 - 2.0 * tr_xxzz_xxxx[i] * tke_0 + 4.0 * tr_xxzz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xxyzz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxzz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xxxz[i] = -4.0 * tr_xzz_xxxyz[i] * tke_0 - 4.0 * tr_xyzz_xxxz[i] * tbe_0 - 6.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xxyy[i] = 4.0 * tr_xzz_xxy[i] - 4.0 * tr_xzz_xxyyy[i] * tke_0 - 4.0 * tr_xyzz_xxyy[i] * tbe_0 + 4.0 * tr_xxzz_xy[i] - 4.0 * tr_xxzz_xyyy[i] * tke_0 - 4.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxzz_xxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xxyz[i] = 2.0 * tr_xzz_xxz[i] - 4.0 * tr_xzz_xxyyz[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tbe_0 + 2.0 * tr_xxzz_xz[i] - 4.0 * tr_xxzz_xyyz[i] * tke_0 - 2.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xyz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxzz_xxz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xxzz[i] = -4.0 * tr_xzz_xxyzz[i] * tke_0 - 4.0 * tr_xyzz_xxzz[i] * tbe_0 - 4.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xyyy[i] = 6.0 * tr_xzz_xyy[i] - 4.0 * tr_xzz_xyyyy[i] * tke_0 - 4.0 * tr_xyzz_xyyy[i] * tbe_0 + 3.0 * tr_xxzz_yy[i] - 2.0 * tr_xxzz_yyyy[i] * tke_0 - 6.0 * tr_xxzz_xxyy[i] * tke_0 + 4.0 * tr_xxzz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxxzz_xyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xyyz[i] = 4.0 * tr_xzz_xyz[i] - 4.0 * tr_xzz_xyyyz[i] * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tbe_0 + 2.0 * tr_xxzz_yz[i] - 2.0 * tr_xxzz_yyyz[i] * tke_0 - 4.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxzz_xyz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xyzz[i] = 2.0 * tr_xzz_xzz[i] - 4.0 * tr_xzz_xyyzz[i] * tke_0 - 4.0 * tr_xyzz_xyzz[i] * tbe_0 + tr_xxzz_zz[i] - 2.0 * tr_xxzz_yyzz[i] * tke_0 - 2.0 * tr_xxzz_xxzz[i] * tke_0 + 4.0 * tr_xxzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxzz_xzz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xzzz[i] = -4.0 * tr_xzz_xyzzz[i] * tke_0 - 4.0 * tr_xyzz_xzzz[i] * tbe_0 - 2.0 * tr_xxzz_yzzz[i] * tke_0 + 4.0 * tr_xxzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_yyyy[i] = 8.0 * tr_xzz_yyy[i] - 4.0 * tr_xzz_yyyyy[i] * tke_0 - 4.0 * tr_xyzz_yyyy[i] * tbe_0 - 8.0 * tr_xxzz_xyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xxxzz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_yyyz[i] = 6.0 * tr_xzz_yyz[i] - 4.0 * tr_xzz_yyyyz[i] * tke_0 - 4.0 * tr_xyzz_yyyz[i] * tbe_0 - 6.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxzz_yyz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_yyzz[i] = 4.0 * tr_xzz_yzz[i] - 4.0 * tr_xzz_yyyzz[i] * tke_0 - 4.0 * tr_xyzz_yyzz[i] * tbe_0 - 4.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxzz_yzz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_yzzz[i] = 2.0 * tr_xzz_zzz[i] - 4.0 * tr_xzz_yyzzz[i] * tke_0 - 4.0 * tr_xyzz_yzzz[i] * tbe_0 - 2.0 * tr_xxzz_xzzz[i] * tke_0 + 4.0 * tr_xxzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxzz_zzz[i] * tbe_0 + 4.0 * tr_xxxzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_zzzz[i] = -4.0 * tr_xzz_yzzzz[i] * tke_0 - 4.0 * tr_xyzz_zzzz[i] * tbe_0 + 4.0 * tr_xxzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 315-330 components of targeted buffer : GG

    auto tr_0_0_xy_xyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 315);

    auto tr_0_0_xy_xyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 316);

    auto tr_0_0_xy_xyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 317);

    auto tr_0_0_xy_xyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 318);

    auto tr_0_0_xy_xyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 319);

    auto tr_0_0_xy_xyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 320);

    auto tr_0_0_xy_xyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 321);

    auto tr_0_0_xy_xyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 322);

    auto tr_0_0_xy_xyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 323);

    auto tr_0_0_xy_xyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 324);

    auto tr_0_0_xy_xyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 325);

    auto tr_0_0_xy_xyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 326);

    auto tr_0_0_xy_xyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 327);

    auto tr_0_0_xy_xyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 328);

    auto tr_0_0_xy_xyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 329);

    #pragma omp simd aligned(tr_0_0_xy_xyyy_xxxx, tr_0_0_xy_xyyy_xxxy, tr_0_0_xy_xyyy_xxxz, tr_0_0_xy_xyyy_xxyy, tr_0_0_xy_xyyy_xxyz, tr_0_0_xy_xyyy_xxzz, tr_0_0_xy_xyyy_xyyy, tr_0_0_xy_xyyy_xyyz, tr_0_0_xy_xyyy_xyzz, tr_0_0_xy_xyyy_xzzz, tr_0_0_xy_xyyy_yyyy, tr_0_0_xy_xyyy_yyyz, tr_0_0_xy_xyyy_yyzz, tr_0_0_xy_xyyy_yzzz, tr_0_0_xy_xyyy_zzzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_xxyyy_xxx, tr_xxyyy_xxxxy, tr_xxyyy_xxxyy, tr_xxyyy_xxxyz, tr_xxyyy_xxy, tr_xxyyy_xxyyy, tr_xxyyy_xxyyz, tr_xxyyy_xxyzz, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyyyy, tr_xxyyy_xyyyz, tr_xxyyy_xyyzz, tr_xxyyy_xyz, tr_xxyyy_xyzzz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyyyy, tr_xxyyy_yyyyz, tr_xxyyy_yyyzz, tr_xxyyy_yyz, tr_xxyyy_yyzzz, tr_xxyyy_yzz, tr_xxyyy_yzzzz, tr_xxyyy_zzz, tr_xxyyyy_xxxx, tr_xxyyyy_xxxy, tr_xxyyyy_xxxz, tr_xxyyyy_xxyy, tr_xxyyyy_xxyz, tr_xxyyyy_xxzz, tr_xxyyyy_xyyy, tr_xxyyyy_xyyz, tr_xxyyyy_xyzz, tr_xxyyyy_xzzz, tr_xxyyyy_yyyy, tr_xxyyyy_yyyz, tr_xxyyyy_yyzz, tr_xxyyyy_yzzz, tr_xxyyyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxxxy, tr_xyyy_xxxxyy, tr_xyyy_xxxxyz, tr_xyyy_xxxy, tr_xyyy_xxxyyy, tr_xyyy_xxxyyz, tr_xyyy_xxxyzz, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyyyy, tr_xyyy_xxyyyz, tr_xyyy_xxyyzz, tr_xyyy_xxyz, tr_xyyy_xxyzzz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyyyy, tr_xyyy_xyyyyz, tr_xyyy_xyyyzz, tr_xyyy_xyyz, tr_xyyy_xyyzzz, tr_xyyy_xyzz, tr_xyyy_xyzzzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyyy_xxx, tr_xyyyy_xxxxx, tr_xyyyy_xxxxy, tr_xyyyy_xxxxz, tr_xyyyy_xxxyy, tr_xyyyy_xxxyz, tr_xyyyy_xxxzz, tr_xyyyy_xxy, tr_xyyyy_xxyyy, tr_xyyyy_xxyyz, tr_xyyyy_xxyzz, tr_xyyyy_xxz, tr_xyyyy_xxzzz, tr_xyyyy_xyy, tr_xyyyy_xyyyy, tr_xyyyy_xyyyz, tr_xyyyy_xyyzz, tr_xyyyy_xyz, tr_xyyyy_xyzzz, tr_xyyyy_xzz, tr_xyyyy_xzzzz, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyy_xxx, tr_yyy_xxxxy, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyyyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xzzz, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yzzz, tr_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyyy_xxxx[i] = 3.0 * tr_yy_xxxx[i] - 2.0 * tr_yyy_xxxxy[i] * tke_0 - 2.0 * tr_yyyy_xxxx[i] * tbe_0 + 12.0 * tr_xyy_xxx[i] - 6.0 * tr_xyy_xxxxx[i] * tke_0 - 8.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xxxy[i] = 3.0 * tr_yy_xxxy[i] + tr_yyy_xxx[i] - 2.0 * tr_yyy_xxxyy[i] * tke_0 - 2.0 * tr_yyyy_xxxy[i] * tbe_0 + 9.0 * tr_xyy_xxy[i] - 6.0 * tr_xyy_xxxxy[i] * tke_0 + 3.0 * tr_xyyy_xx[i] - 6.0 * tr_xyyy_xxyy[i] * tke_0 - 2.0 * tr_xyyy_xxxx[i] * tke_0 + 4.0 * tr_xyyy_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xyyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyy_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxxy[i] * tbe_0 - 2.0 * tr_xxyyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xxxz[i] = 3.0 * tr_yy_xxxz[i] - 2.0 * tr_yyy_xxxyz[i] * tke_0 - 2.0 * tr_yyyy_xxxz[i] * tbe_0 + 9.0 * tr_xyy_xxz[i] - 6.0 * tr_xyy_xxxxz[i] * tke_0 - 6.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xxyy[i] = 3.0 * tr_yy_xxyy[i] + 2.0 * tr_yyy_xxy[i] - 2.0 * tr_yyy_xxyyy[i] * tke_0 - 2.0 * tr_yyyy_xxyy[i] * tbe_0 + 6.0 * tr_xyy_xyy[i] - 6.0 * tr_xyy_xxxyy[i] * tke_0 + 4.0 * tr_xyyy_xy[i] - 4.0 * tr_xyyy_xyyy[i] * tke_0 - 4.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyy_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxyy[i] * tbe_0 - 4.0 * tr_xxyyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xxyz[i] = 3.0 * tr_yy_xxyz[i] + tr_yyy_xxz[i] - 2.0 * tr_yyy_xxyyz[i] * tke_0 - 2.0 * tr_yyyy_xxyz[i] * tbe_0 + 6.0 * tr_xyy_xyz[i] - 6.0 * tr_xyy_xxxyz[i] * tke_0 + 2.0 * tr_xyyy_xz[i] - 4.0 * tr_xyyy_xyyz[i] * tke_0 - 2.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxyz[i] * tbe_0 - 2.0 * tr_xxyyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xxzz[i] = 3.0 * tr_yy_xxzz[i] - 2.0 * tr_yyy_xxyzz[i] * tke_0 - 2.0 * tr_yyyy_xxzz[i] * tbe_0 + 6.0 * tr_xyy_xzz[i] - 6.0 * tr_xyy_xxxzz[i] * tke_0 - 4.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxzz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xyyy[i] = 3.0 * tr_yy_xyyy[i] + 3.0 * tr_yyy_xyy[i] - 2.0 * tr_yyy_xyyyy[i] * tke_0 - 2.0 * tr_yyyy_xyyy[i] * tbe_0 + 3.0 * tr_xyy_yyy[i] - 6.0 * tr_xyy_xxyyy[i] * tke_0 + 3.0 * tr_xyyy_yy[i] - 2.0 * tr_xyyy_yyyy[i] * tke_0 - 6.0 * tr_xyyy_xxyy[i] * tke_0 + 4.0 * tr_xyyy_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xyyy[i] * tbe_0 - 6.0 * tr_xxyyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xyyz[i] = 3.0 * tr_yy_xyyz[i] + 2.0 * tr_yyy_xyz[i] - 2.0 * tr_yyy_xyyyz[i] * tke_0 - 2.0 * tr_yyyy_xyyz[i] * tbe_0 + 3.0 * tr_xyy_yyz[i] - 6.0 * tr_xyy_xxyyz[i] * tke_0 + 2.0 * tr_xyyy_yz[i] - 2.0 * tr_xyyy_yyyz[i] * tke_0 - 4.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xyyz[i] * tbe_0 - 4.0 * tr_xxyyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xyzz[i] = 3.0 * tr_yy_xyzz[i] + tr_yyy_xzz[i] - 2.0 * tr_yyy_xyyzz[i] * tke_0 - 2.0 * tr_yyyy_xyzz[i] * tbe_0 + 3.0 * tr_xyy_yzz[i] - 6.0 * tr_xyy_xxyzz[i] * tke_0 + tr_xyyy_zz[i] - 2.0 * tr_xyyy_yyzz[i] * tke_0 - 2.0 * tr_xyyy_xxzz[i] * tke_0 + 4.0 * tr_xyyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xyzz[i] * tbe_0 - 2.0 * tr_xxyyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xzzz[i] = 3.0 * tr_yy_xzzz[i] - 2.0 * tr_yyy_xyzzz[i] * tke_0 - 2.0 * tr_yyyy_xzzz[i] * tbe_0 + 3.0 * tr_xyy_zzz[i] - 6.0 * tr_xyy_xxzzz[i] * tke_0 - 2.0 * tr_xyyy_yzzz[i] * tke_0 + 4.0 * tr_xyyy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xzzz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_yyyy[i] = 3.0 * tr_yy_yyyy[i] + 4.0 * tr_yyy_yyy[i] - 2.0 * tr_yyy_yyyyy[i] * tke_0 - 2.0 * tr_yyyy_yyyy[i] * tbe_0 - 6.0 * tr_xyy_xyyyy[i] * tke_0 - 8.0 * tr_xyyy_xyyy[i] * tke_0 + 4.0 * tr_xyyy_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yyyy[i] * tbe_0 - 8.0 * tr_xxyyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_yyyz[i] = 3.0 * tr_yy_yyyz[i] + 3.0 * tr_yyy_yyz[i] - 2.0 * tr_yyy_yyyyz[i] * tke_0 - 2.0 * tr_yyyy_yyyz[i] * tbe_0 - 6.0 * tr_xyy_xyyyz[i] * tke_0 - 6.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yyyz[i] * tbe_0 - 6.0 * tr_xxyyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_yyzz[i] = 3.0 * tr_yy_yyzz[i] + 2.0 * tr_yyy_yzz[i] - 2.0 * tr_yyy_yyyzz[i] * tke_0 - 2.0 * tr_yyyy_yyzz[i] * tbe_0 - 6.0 * tr_xyy_xyyzz[i] * tke_0 - 4.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yyzz[i] * tbe_0 - 4.0 * tr_xxyyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_yzzz[i] = 3.0 * tr_yy_yzzz[i] + tr_yyy_zzz[i] - 2.0 * tr_yyy_yyzzz[i] * tke_0 - 2.0 * tr_yyyy_yzzz[i] * tbe_0 - 6.0 * tr_xyy_xyzzz[i] * tke_0 - 2.0 * tr_xyyy_xzzz[i] * tke_0 + 4.0 * tr_xyyy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yzzz[i] * tbe_0 - 2.0 * tr_xxyyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_zzzz[i] = 3.0 * tr_yy_zzzz[i] - 2.0 * tr_yyy_yzzzz[i] * tke_0 - 2.0 * tr_yyyy_zzzz[i] * tbe_0 - 6.0 * tr_xyy_xzzzz[i] * tke_0 + 4.0 * tr_xyyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_zzzz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-345 components of targeted buffer : GG

    auto tr_0_0_xy_xyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 330);

    auto tr_0_0_xy_xyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 331);

    auto tr_0_0_xy_xyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 332);

    auto tr_0_0_xy_xyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 333);

    auto tr_0_0_xy_xyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 334);

    auto tr_0_0_xy_xyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 335);

    auto tr_0_0_xy_xyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 336);

    auto tr_0_0_xy_xyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 337);

    auto tr_0_0_xy_xyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 338);

    auto tr_0_0_xy_xyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 339);

    auto tr_0_0_xy_xyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 340);

    auto tr_0_0_xy_xyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 341);

    auto tr_0_0_xy_xyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 342);

    auto tr_0_0_xy_xyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 343);

    auto tr_0_0_xy_xyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 344);

    #pragma omp simd aligned(tr_0_0_xy_xyyz_xxxx, tr_0_0_xy_xyyz_xxxy, tr_0_0_xy_xyyz_xxxz, tr_0_0_xy_xyyz_xxyy, tr_0_0_xy_xyyz_xxyz, tr_0_0_xy_xyyz_xxzz, tr_0_0_xy_xyyz_xyyy, tr_0_0_xy_xyyz_xyyz, tr_0_0_xy_xyyz_xyzz, tr_0_0_xy_xyyz_xzzz, tr_0_0_xy_xyyz_yyyy, tr_0_0_xy_xyyz_yyyz, tr_0_0_xy_xyyz_yyzz, tr_0_0_xy_xyyz_yzzz, tr_0_0_xy_xyyz_zzzz, tr_xxyyyz_xxxx, tr_xxyyyz_xxxy, tr_xxyyyz_xxxz, tr_xxyyyz_xxyy, tr_xxyyyz_xxyz, tr_xxyyyz_xxzz, tr_xxyyyz_xyyy, tr_xxyyyz_xyyz, tr_xxyyyz_xyzz, tr_xxyyyz_xzzz, tr_xxyyyz_yyyy, tr_xxyyyz_yyyz, tr_xxyyyz_yyzz, tr_xxyyyz_yzzz, tr_xxyyyz_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxy, tr_xxyyz_xxxyy, tr_xxyyz_xxxyz, tr_xxyyz_xxy, tr_xxyyz_xxyyy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyyyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyyyy, tr_xxyyz_yyyyz, tr_xxyyz_yyyzz, tr_xxyyz_yyz, tr_xxyyz_yyzzz, tr_xxyyz_yzz, tr_xxyyz_yzzzz, tr_xxyyz_zzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxx, tr_xyyyz_xxxxy, tr_xyyyz_xxxxz, tr_xyyyz_xxxyy, tr_xyyyz_xxxyz, tr_xyyyz_xxxzz, tr_xyyyz_xxy, tr_xyyyz_xxyyy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xxzzz, tr_xyyyz_xyy, tr_xyyyz_xyyyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_xzzzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxxxy, tr_xyyz_xxxxyy, tr_xyyz_xxxxyz, tr_xyyz_xxxy, tr_xyyz_xxxyyy, tr_xyyz_xxxyyz, tr_xyyz_xxxyzz, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyyyy, tr_xyyz_xxyyyz, tr_xyyz_xxyyzz, tr_xyyz_xxyz, tr_xyyz_xxyzzz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyyyy, tr_xyyz_xyyyyz, tr_xyyz_xyyyzz, tr_xyyz_xyyz, tr_xyyz_xyyzzz, tr_xyyz_xyzz, tr_xyyz_xyzzzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyyz_xxxx[i] = 2.0 * tr_yz_xxxx[i] - 2.0 * tr_yyz_xxxxy[i] * tke_0 - 2.0 * tr_yyyz_xxxx[i] * tbe_0 + 8.0 * tr_xyz_xxx[i] - 4.0 * tr_xyz_xxxxx[i] * tke_0 - 8.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xxxy[i] = 2.0 * tr_yz_xxxy[i] + tr_yyz_xxx[i] - 2.0 * tr_yyz_xxxyy[i] * tke_0 - 2.0 * tr_yyyz_xxxy[i] * tbe_0 + 6.0 * tr_xyz_xxy[i] - 4.0 * tr_xyz_xxxxy[i] * tke_0 + 3.0 * tr_xyyz_xx[i] - 6.0 * tr_xyyz_xxyy[i] * tke_0 - 2.0 * tr_xyyz_xxxx[i] * tke_0 + 4.0 * tr_xyyz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xyyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxy[i] * tbe_0 - 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xxxz[i] = 2.0 * tr_yz_xxxz[i] - 2.0 * tr_yyz_xxxyz[i] * tke_0 - 2.0 * tr_yyyz_xxxz[i] * tbe_0 + 6.0 * tr_xyz_xxz[i] - 4.0 * tr_xyz_xxxxz[i] * tke_0 - 6.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xxyy[i] = 2.0 * tr_yz_xxyy[i] + 2.0 * tr_yyz_xxy[i] - 2.0 * tr_yyz_xxyyy[i] * tke_0 - 2.0 * tr_yyyz_xxyy[i] * tbe_0 + 4.0 * tr_xyz_xyy[i] - 4.0 * tr_xyz_xxxyy[i] * tke_0 + 4.0 * tr_xyyz_xy[i] - 4.0 * tr_xyyz_xyyy[i] * tke_0 - 4.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxyy[i] * tbe_0 - 4.0 * tr_xxyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xxyz[i] = 2.0 * tr_yz_xxyz[i] + tr_yyz_xxz[i] - 2.0 * tr_yyz_xxyyz[i] * tke_0 - 2.0 * tr_yyyz_xxyz[i] * tbe_0 + 4.0 * tr_xyz_xyz[i] - 4.0 * tr_xyz_xxxyz[i] * tke_0 + 2.0 * tr_xyyz_xz[i] - 4.0 * tr_xyyz_xyyz[i] * tke_0 - 2.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tbe_0 - 2.0 * tr_xxyyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xxzz[i] = 2.0 * tr_yz_xxzz[i] - 2.0 * tr_yyz_xxyzz[i] * tke_0 - 2.0 * tr_yyyz_xxzz[i] * tbe_0 + 4.0 * tr_xyz_xzz[i] - 4.0 * tr_xyz_xxxzz[i] * tke_0 - 4.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xyyy[i] = 2.0 * tr_yz_xyyy[i] + 3.0 * tr_yyz_xyy[i] - 2.0 * tr_yyz_xyyyy[i] * tke_0 - 2.0 * tr_yyyz_xyyy[i] * tbe_0 + 2.0 * tr_xyz_yyy[i] - 4.0 * tr_xyz_xxyyy[i] * tke_0 + 3.0 * tr_xyyz_yy[i] - 2.0 * tr_xyyz_yyyy[i] * tke_0 - 6.0 * tr_xyyz_xxyy[i] * tke_0 + 4.0 * tr_xyyz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyyy[i] * tbe_0 - 6.0 * tr_xxyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xyyz[i] = 2.0 * tr_yz_xyyz[i] + 2.0 * tr_yyz_xyz[i] - 2.0 * tr_yyz_xyyyz[i] * tke_0 - 2.0 * tr_yyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyz_yyz[i] - 4.0 * tr_xyz_xxyyz[i] * tke_0 + 2.0 * tr_xyyz_yz[i] - 2.0 * tr_xyyz_yyyz[i] * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tbe_0 - 4.0 * tr_xxyyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xyzz[i] = 2.0 * tr_yz_xyzz[i] + tr_yyz_xzz[i] - 2.0 * tr_yyz_xyyzz[i] * tke_0 - 2.0 * tr_yyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyz_yzz[i] - 4.0 * tr_xyz_xxyzz[i] * tke_0 + tr_xyyz_zz[i] - 2.0 * tr_xyyz_yyzz[i] * tke_0 - 2.0 * tr_xyyz_xxzz[i] * tke_0 + 4.0 * tr_xyyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyzz[i] * tbe_0 - 2.0 * tr_xxyyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xzzz[i] = 2.0 * tr_yz_xzzz[i] - 2.0 * tr_yyz_xyzzz[i] * tke_0 - 2.0 * tr_yyyz_xzzz[i] * tbe_0 + 2.0 * tr_xyz_zzz[i] - 4.0 * tr_xyz_xxzzz[i] * tke_0 - 2.0 * tr_xyyz_yzzz[i] * tke_0 + 4.0 * tr_xyyz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xzzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_yyyy[i] = 2.0 * tr_yz_yyyy[i] + 4.0 * tr_yyz_yyy[i] - 2.0 * tr_yyz_yyyyy[i] * tke_0 - 2.0 * tr_yyyz_yyyy[i] * tbe_0 - 4.0 * tr_xyz_xyyyy[i] * tke_0 - 8.0 * tr_xyyz_xyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyyy[i] * tbe_0 - 8.0 * tr_xxyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_yyyz[i] = 2.0 * tr_yz_yyyz[i] + 3.0 * tr_yyz_yyz[i] - 2.0 * tr_yyz_yyyyz[i] * tke_0 - 2.0 * tr_yyyz_yyyz[i] * tbe_0 - 4.0 * tr_xyz_xyyyz[i] * tke_0 - 6.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyyz[i] * tbe_0 - 6.0 * tr_xxyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_yyzz[i] = 2.0 * tr_yz_yyzz[i] + 2.0 * tr_yyz_yzz[i] - 2.0 * tr_yyz_yyyzz[i] * tke_0 - 2.0 * tr_yyyz_yyzz[i] * tbe_0 - 4.0 * tr_xyz_xyyzz[i] * tke_0 - 4.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyzz[i] * tbe_0 - 4.0 * tr_xxyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_yzzz[i] = 2.0 * tr_yz_yzzz[i] + tr_yyz_zzz[i] - 2.0 * tr_yyz_yyzzz[i] * tke_0 - 2.0 * tr_yyyz_yzzz[i] * tbe_0 - 4.0 * tr_xyz_xyzzz[i] * tke_0 - 2.0 * tr_xyyz_xzzz[i] * tke_0 + 4.0 * tr_xyyz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yzzz[i] * tbe_0 - 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_zzzz[i] = 2.0 * tr_yz_zzzz[i] - 2.0 * tr_yyz_yzzzz[i] * tke_0 - 2.0 * tr_yyyz_zzzz[i] * tbe_0 - 4.0 * tr_xyz_xzzzz[i] * tke_0 + 4.0 * tr_xyyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zzzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 345-360 components of targeted buffer : GG

    auto tr_0_0_xy_xyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 345);

    auto tr_0_0_xy_xyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 346);

    auto tr_0_0_xy_xyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 347);

    auto tr_0_0_xy_xyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 348);

    auto tr_0_0_xy_xyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 349);

    auto tr_0_0_xy_xyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 350);

    auto tr_0_0_xy_xyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 351);

    auto tr_0_0_xy_xyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 352);

    auto tr_0_0_xy_xyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 353);

    auto tr_0_0_xy_xyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 354);

    auto tr_0_0_xy_xyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 355);

    auto tr_0_0_xy_xyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 356);

    auto tr_0_0_xy_xyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 357);

    auto tr_0_0_xy_xyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 358);

    auto tr_0_0_xy_xyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 359);

    #pragma omp simd aligned(tr_0_0_xy_xyzz_xxxx, tr_0_0_xy_xyzz_xxxy, tr_0_0_xy_xyzz_xxxz, tr_0_0_xy_xyzz_xxyy, tr_0_0_xy_xyzz_xxyz, tr_0_0_xy_xyzz_xxzz, tr_0_0_xy_xyzz_xyyy, tr_0_0_xy_xyzz_xyyz, tr_0_0_xy_xyzz_xyzz, tr_0_0_xy_xyzz_xzzz, tr_0_0_xy_xyzz_yyyy, tr_0_0_xy_xyzz_yyyz, tr_0_0_xy_xyzz_yyzz, tr_0_0_xy_xyzz_yzzz, tr_0_0_xy_xyzz_zzzz, tr_xxyyzz_xxxx, tr_xxyyzz_xxxy, tr_xxyyzz_xxxz, tr_xxyyzz_xxyy, tr_xxyyzz_xxyz, tr_xxyyzz_xxzz, tr_xxyyzz_xyyy, tr_xxyyzz_xyyz, tr_xxyyzz_xyzz, tr_xxyyzz_xzzz, tr_xxyyzz_yyyy, tr_xxyyzz_yyyz, tr_xxyyzz_yyzz, tr_xxyyzz_yzzz, tr_xxyyzz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxy, tr_xxyzz_xxxyy, tr_xxyzz_xxxyz, tr_xxyzz_xxy, tr_xxyzz_xxyyy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyyyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyyyy, tr_xxyzz_yyyyz, tr_xxyzz_yyyzz, tr_xxyzz_yyz, tr_xxyzz_yyzzz, tr_xxyzz_yzz, tr_xxyzz_yzzzz, tr_xxyzz_zzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxx, tr_xyyzz_xxxxy, tr_xyyzz_xxxxz, tr_xyyzz_xxxyy, tr_xyyzz_xxxyz, tr_xyyzz_xxxzz, tr_xyyzz_xxy, tr_xyyzz_xxyyy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xxzzz, tr_xyyzz_xyy, tr_xyyzz_xyyyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_xzzzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxxxy, tr_xyzz_xxxxyy, tr_xyzz_xxxxyz, tr_xyzz_xxxy, tr_xyzz_xxxyyy, tr_xyzz_xxxyyz, tr_xyzz_xxxyzz, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyyyy, tr_xyzz_xxyyyz, tr_xyzz_xxyyzz, tr_xyzz_xxyz, tr_xyzz_xxyzzz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyyyy, tr_xyzz_xyyyyz, tr_xyzz_xyyyzz, tr_xyzz_xyyz, tr_xyzz_xyyzzz, tr_xyzz_xyzz, tr_xyzz_xyzzzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyzz_xxxx[i] = tr_zz_xxxx[i] - 2.0 * tr_yzz_xxxxy[i] * tke_0 - 2.0 * tr_yyzz_xxxx[i] * tbe_0 + 4.0 * tr_xzz_xxx[i] - 2.0 * tr_xzz_xxxxx[i] * tke_0 - 8.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xxxy[i] = tr_zz_xxxy[i] + tr_yzz_xxx[i] - 2.0 * tr_yzz_xxxyy[i] * tke_0 - 2.0 * tr_yyzz_xxxy[i] * tbe_0 + 3.0 * tr_xzz_xxy[i] - 2.0 * tr_xzz_xxxxy[i] * tke_0 + 3.0 * tr_xyzz_xx[i] - 6.0 * tr_xyzz_xxyy[i] * tke_0 - 2.0 * tr_xyzz_xxxx[i] * tke_0 + 4.0 * tr_xyzz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xyyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxxy[i] * tbe_0 - 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xxxz[i] = tr_zz_xxxz[i] - 2.0 * tr_yzz_xxxyz[i] * tke_0 - 2.0 * tr_yyzz_xxxz[i] * tbe_0 + 3.0 * tr_xzz_xxz[i] - 2.0 * tr_xzz_xxxxz[i] * tke_0 - 6.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xxyy[i] = tr_zz_xxyy[i] + 2.0 * tr_yzz_xxy[i] - 2.0 * tr_yzz_xxyyy[i] * tke_0 - 2.0 * tr_yyzz_xxyy[i] * tbe_0 + 2.0 * tr_xzz_xyy[i] - 2.0 * tr_xzz_xxxyy[i] * tke_0 + 4.0 * tr_xyzz_xy[i] - 4.0 * tr_xyzz_xyyy[i] * tke_0 - 4.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxyy[i] * tbe_0 - 4.0 * tr_xxyzz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xxyz[i] = tr_zz_xxyz[i] + tr_yzz_xxz[i] - 2.0 * tr_yzz_xxyyz[i] * tke_0 - 2.0 * tr_yyzz_xxyz[i] * tbe_0 + 2.0 * tr_xzz_xyz[i] - 2.0 * tr_xzz_xxxyz[i] * tke_0 + 2.0 * tr_xyzz_xz[i] - 4.0 * tr_xyzz_xyyz[i] * tke_0 - 2.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxyz[i] * tbe_0 - 2.0 * tr_xxyzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xxzz[i] = tr_zz_xxzz[i] - 2.0 * tr_yzz_xxyzz[i] * tke_0 - 2.0 * tr_yyzz_xxzz[i] * tbe_0 + 2.0 * tr_xzz_xzz[i] - 2.0 * tr_xzz_xxxzz[i] * tke_0 - 4.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xyyy[i] = tr_zz_xyyy[i] + 3.0 * tr_yzz_xyy[i] - 2.0 * tr_yzz_xyyyy[i] * tke_0 - 2.0 * tr_yyzz_xyyy[i] * tbe_0 + tr_xzz_yyy[i] - 2.0 * tr_xzz_xxyyy[i] * tke_0 + 3.0 * tr_xyzz_yy[i] - 2.0 * tr_xyzz_yyyy[i] * tke_0 - 6.0 * tr_xyzz_xxyy[i] * tke_0 + 4.0 * tr_xyzz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xyyy[i] * tbe_0 - 6.0 * tr_xxyzz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xyyz[i] = tr_zz_xyyz[i] + 2.0 * tr_yzz_xyz[i] - 2.0 * tr_yzz_xyyyz[i] * tke_0 - 2.0 * tr_yyzz_xyyz[i] * tbe_0 + tr_xzz_yyz[i] - 2.0 * tr_xzz_xxyyz[i] * tke_0 + 2.0 * tr_xyzz_yz[i] - 2.0 * tr_xyzz_yyyz[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xyyz[i] * tbe_0 - 4.0 * tr_xxyzz_xyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xyzz[i] = tr_zz_xyzz[i] + tr_yzz_xzz[i] - 2.0 * tr_yzz_xyyzz[i] * tke_0 - 2.0 * tr_yyzz_xyzz[i] * tbe_0 + tr_xzz_yzz[i] - 2.0 * tr_xzz_xxyzz[i] * tke_0 + tr_xyzz_zz[i] - 2.0 * tr_xyzz_yyzz[i] * tke_0 - 2.0 * tr_xyzz_xxzz[i] * tke_0 + 4.0 * tr_xyzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xyzz[i] * tbe_0 - 2.0 * tr_xxyzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xzzz[i] = tr_zz_xzzz[i] - 2.0 * tr_yzz_xyzzz[i] * tke_0 - 2.0 * tr_yyzz_xzzz[i] * tbe_0 + tr_xzz_zzz[i] - 2.0 * tr_xzz_xxzzz[i] * tke_0 - 2.0 * tr_xyzz_yzzz[i] * tke_0 + 4.0 * tr_xyzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xzzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_yyyy[i] = tr_zz_yyyy[i] + 4.0 * tr_yzz_yyy[i] - 2.0 * tr_yzz_yyyyy[i] * tke_0 - 2.0 * tr_yyzz_yyyy[i] * tbe_0 - 2.0 * tr_xzz_xyyyy[i] * tke_0 - 8.0 * tr_xyzz_xyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yyyy[i] * tbe_0 - 8.0 * tr_xxyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_yyyz[i] = tr_zz_yyyz[i] + 3.0 * tr_yzz_yyz[i] - 2.0 * tr_yzz_yyyyz[i] * tke_0 - 2.0 * tr_yyzz_yyyz[i] * tbe_0 - 2.0 * tr_xzz_xyyyz[i] * tke_0 - 6.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yyyz[i] * tbe_0 - 6.0 * tr_xxyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_yyzz[i] = tr_zz_yyzz[i] + 2.0 * tr_yzz_yzz[i] - 2.0 * tr_yzz_yyyzz[i] * tke_0 - 2.0 * tr_yyzz_yyzz[i] * tbe_0 - 2.0 * tr_xzz_xyyzz[i] * tke_0 - 4.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yyzz[i] * tbe_0 - 4.0 * tr_xxyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_yzzz[i] = tr_zz_yzzz[i] + tr_yzz_zzz[i] - 2.0 * tr_yzz_yyzzz[i] * tke_0 - 2.0 * tr_yyzz_yzzz[i] * tbe_0 - 2.0 * tr_xzz_xyzzz[i] * tke_0 - 2.0 * tr_xyzz_xzzz[i] * tke_0 + 4.0 * tr_xyzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yzzz[i] * tbe_0 - 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_zzzz[i] = tr_zz_zzzz[i] - 2.0 * tr_yzz_yzzzz[i] * tke_0 - 2.0 * tr_yyzz_zzzz[i] * tbe_0 - 2.0 * tr_xzz_xzzzz[i] * tke_0 + 4.0 * tr_xyzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_zzzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-375 components of targeted buffer : GG

    auto tr_0_0_xy_xzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 360);

    auto tr_0_0_xy_xzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 361);

    auto tr_0_0_xy_xzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 362);

    auto tr_0_0_xy_xzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 363);

    auto tr_0_0_xy_xzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 364);

    auto tr_0_0_xy_xzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 365);

    auto tr_0_0_xy_xzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 366);

    auto tr_0_0_xy_xzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 367);

    auto tr_0_0_xy_xzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 368);

    auto tr_0_0_xy_xzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 369);

    auto tr_0_0_xy_xzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 370);

    auto tr_0_0_xy_xzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 371);

    auto tr_0_0_xy_xzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 372);

    auto tr_0_0_xy_xzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 373);

    auto tr_0_0_xy_xzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 374);

    #pragma omp simd aligned(tr_0_0_xy_xzzz_xxxx, tr_0_0_xy_xzzz_xxxy, tr_0_0_xy_xzzz_xxxz, tr_0_0_xy_xzzz_xxyy, tr_0_0_xy_xzzz_xxyz, tr_0_0_xy_xzzz_xxzz, tr_0_0_xy_xzzz_xyyy, tr_0_0_xy_xzzz_xyyz, tr_0_0_xy_xzzz_xyzz, tr_0_0_xy_xzzz_xzzz, tr_0_0_xy_xzzz_yyyy, tr_0_0_xy_xzzz_yyyz, tr_0_0_xy_xzzz_yyzz, tr_0_0_xy_xzzz_yzzz, tr_0_0_xy_xzzz_zzzz, tr_xxyzzz_xxxx, tr_xxyzzz_xxxy, tr_xxyzzz_xxxz, tr_xxyzzz_xxyy, tr_xxyzzz_xxyz, tr_xxyzzz_xxzz, tr_xxyzzz_xyyy, tr_xxyzzz_xyyz, tr_xxyzzz_xyzz, tr_xxyzzz_xzzz, tr_xxyzzz_yyyy, tr_xxyzzz_yyyz, tr_xxyzzz_yyzz, tr_xxyzzz_yzzz, tr_xxyzzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxxxy, tr_xxzzz_xxxyy, tr_xxzzz_xxxyz, tr_xxzzz_xxy, tr_xxzzz_xxyyy, tr_xxzzz_xxyyz, tr_xxzzz_xxyzz, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyyyy, tr_xxzzz_xyyyz, tr_xxzzz_xyyzz, tr_xxzzz_xyz, tr_xxzzz_xyzzz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyyyy, tr_xxzzz_yyyyz, tr_xxzzz_yyyzz, tr_xxzzz_yyz, tr_xxzzz_yyzzz, tr_xxzzz_yzz, tr_xxzzz_yzzzz, tr_xxzzz_zzz, tr_xyzzz_xxx, tr_xyzzz_xxxxx, tr_xyzzz_xxxxy, tr_xyzzz_xxxxz, tr_xyzzz_xxxyy, tr_xyzzz_xxxyz, tr_xyzzz_xxxzz, tr_xyzzz_xxy, tr_xyzzz_xxyyy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xxzzz, tr_xyzzz_xyy, tr_xyzzz_xyyyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_xzzzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxxxy, tr_xzzz_xxxxyy, tr_xzzz_xxxxyz, tr_xzzz_xxxy, tr_xzzz_xxxyyy, tr_xzzz_xxxyyz, tr_xzzz_xxxyzz, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyyyy, tr_xzzz_xxyyyz, tr_xzzz_xxyyzz, tr_xzzz_xxyz, tr_xzzz_xxyzzz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyyyy, tr_xzzz_xyyyyz, tr_xzzz_xyyyzz, tr_xzzz_xyyz, tr_xzzz_xyyzzz, tr_xzzz_xyzz, tr_xzzz_xyzzzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxy, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyyyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xzzz_xxxx[i] = -2.0 * tr_zzz_xxxxy[i] * tke_0 - 2.0 * tr_yzzz_xxxx[i] * tbe_0 - 8.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xxxy[i] = tr_zzz_xxx[i] - 2.0 * tr_zzz_xxxyy[i] * tke_0 - 2.0 * tr_yzzz_xxxy[i] * tbe_0 + 3.0 * tr_xzzz_xx[i] - 6.0 * tr_xzzz_xxyy[i] * tke_0 - 2.0 * tr_xzzz_xxxx[i] * tke_0 + 4.0 * tr_xzzz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_xyzzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xxxz[i] = -2.0 * tr_zzz_xxxyz[i] * tke_0 - 2.0 * tr_yzzz_xxxz[i] * tbe_0 - 6.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xxyy[i] = 2.0 * tr_zzz_xxy[i] - 2.0 * tr_zzz_xxyyy[i] * tke_0 - 2.0 * tr_yzzz_xxyy[i] * tbe_0 + 4.0 * tr_xzzz_xy[i] - 4.0 * tr_xzzz_xyyy[i] * tke_0 - 4.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxzzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xxyz[i] = tr_zzz_xxz[i] - 2.0 * tr_zzz_xxyyz[i] * tke_0 - 2.0 * tr_yzzz_xxyz[i] * tbe_0 + 2.0 * tr_xzzz_xz[i] - 4.0 * tr_xzzz_xyyz[i] * tke_0 - 2.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzzz_xxz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xxzz[i] = -2.0 * tr_zzz_xxyzz[i] * tke_0 - 2.0 * tr_yzzz_xxzz[i] * tbe_0 - 4.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xyyy[i] = 3.0 * tr_zzz_xyy[i] - 2.0 * tr_zzz_xyyyy[i] * tke_0 - 2.0 * tr_yzzz_xyyy[i] * tbe_0 + 3.0 * tr_xzzz_yy[i] - 2.0 * tr_xzzz_yyyy[i] * tke_0 - 6.0 * tr_xzzz_xxyy[i] * tke_0 + 4.0 * tr_xzzz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xyyz[i] = 2.0 * tr_zzz_xyz[i] - 2.0 * tr_zzz_xyyyz[i] * tke_0 - 2.0 * tr_yzzz_xyyz[i] * tbe_0 + 2.0 * tr_xzzz_yz[i] - 2.0 * tr_xzzz_yyyz[i] * tke_0 - 4.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxzzz_xyz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xyzz[i] = tr_zzz_xzz[i] - 2.0 * tr_zzz_xyyzz[i] * tke_0 - 2.0 * tr_yzzz_xyzz[i] * tbe_0 + tr_xzzz_zz[i] - 2.0 * tr_xzzz_yyzz[i] * tke_0 - 2.0 * tr_xzzz_xxzz[i] * tke_0 + 4.0 * tr_xzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzzz_xzz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xzzz[i] = -2.0 * tr_zzz_xyzzz[i] * tke_0 - 2.0 * tr_yzzz_xzzz[i] * tbe_0 - 2.0 * tr_xzzz_yzzz[i] * tke_0 + 4.0 * tr_xzzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_yyyy[i] = 4.0 * tr_zzz_yyy[i] - 2.0 * tr_zzz_yyyyy[i] * tke_0 - 2.0 * tr_yzzz_yyyy[i] * tbe_0 - 8.0 * tr_xzzz_xyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xxzzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_yyyz[i] = 3.0 * tr_zzz_yyz[i] - 2.0 * tr_zzz_yyyyz[i] * tke_0 - 2.0 * tr_yzzz_yyyz[i] * tbe_0 - 6.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzzz_yyz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_yyzz[i] = 2.0 * tr_zzz_yzz[i] - 2.0 * tr_zzz_yyyzz[i] * tke_0 - 2.0 * tr_yzzz_yyzz[i] * tbe_0 - 4.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxzzz_yzz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_yzzz[i] = tr_zzz_zzz[i] - 2.0 * tr_zzz_yyzzz[i] * tke_0 - 2.0 * tr_yzzz_yzzz[i] * tbe_0 - 2.0 * tr_xzzz_xzzz[i] * tke_0 + 4.0 * tr_xzzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_zzzz[i] = -2.0 * tr_zzz_yzzzz[i] * tke_0 - 2.0 * tr_yzzz_zzzz[i] * tbe_0 + 4.0 * tr_xzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 375-390 components of targeted buffer : GG

    auto tr_0_0_xy_yyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 375);

    auto tr_0_0_xy_yyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 376);

    auto tr_0_0_xy_yyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 377);

    auto tr_0_0_xy_yyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 378);

    auto tr_0_0_xy_yyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 379);

    auto tr_0_0_xy_yyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 380);

    auto tr_0_0_xy_yyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 381);

    auto tr_0_0_xy_yyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 382);

    auto tr_0_0_xy_yyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 383);

    auto tr_0_0_xy_yyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 384);

    auto tr_0_0_xy_yyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 385);

    auto tr_0_0_xy_yyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 386);

    auto tr_0_0_xy_yyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 387);

    auto tr_0_0_xy_yyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 388);

    auto tr_0_0_xy_yyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 389);

    #pragma omp simd aligned(tr_0_0_xy_yyyy_xxxx, tr_0_0_xy_yyyy_xxxy, tr_0_0_xy_yyyy_xxxz, tr_0_0_xy_yyyy_xxyy, tr_0_0_xy_yyyy_xxyz, tr_0_0_xy_yyyy_xxzz, tr_0_0_xy_yyyy_xyyy, tr_0_0_xy_yyyy_xyyz, tr_0_0_xy_yyyy_xyzz, tr_0_0_xy_yyyy_xzzz, tr_0_0_xy_yyyy_yyyy, tr_0_0_xy_yyyy_yyyz, tr_0_0_xy_yyyy_yyzz, tr_0_0_xy_yyyy_yzzz, tr_0_0_xy_yyyy_zzzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, tr_xyyyy_xxx, tr_xyyyy_xxxxy, tr_xyyyy_xxxyy, tr_xyyyy_xxxyz, tr_xyyyy_xxy, tr_xyyyy_xxyyy, tr_xyyyy_xxyyz, tr_xyyyy_xxyzz, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyyyy, tr_xyyyy_xyyyz, tr_xyyyy_xyyzz, tr_xyyyy_xyz, tr_xyyyy_xyzzz, tr_xyyyy_xzz, tr_xyyyy_yyy, tr_xyyyy_yyyyy, tr_xyyyy_yyyyz, tr_xyyyy_yyyzz, tr_xyyyy_yyz, tr_xyyyy_yyzzz, tr_xyyyy_yzz, tr_xyyyy_yzzzz, tr_xyyyy_zzz, tr_xyyyyy_xxxx, tr_xyyyyy_xxxy, tr_xyyyyy_xxxz, tr_xyyyyy_xxyy, tr_xyyyyy_xxyz, tr_xyyyyy_xxzz, tr_xyyyyy_xyyy, tr_xyyyyy_xyyz, tr_xyyyyy_xyzz, tr_xyyyyy_xzzz, tr_xyyyyy_yyyy, tr_xyyyyy_yyyz, tr_xyyyyy_yyzz, tr_xyyyyy_yzzz, tr_xyyyyy_zzzz, tr_yyy_xxx, tr_yyy_xxxxx, tr_yyy_xxxxy, tr_yyy_xxxxz, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxxxy, tr_yyyy_xxxxyy, tr_yyyy_xxxxyz, tr_yyyy_xxxy, tr_yyyy_xxxyyy, tr_yyyy_xxxyyz, tr_yyyy_xxxyzz, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyyyy, tr_yyyy_xxyyyz, tr_yyyy_xxyyzz, tr_yyyy_xxyz, tr_yyyy_xxyzzz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyyyy, tr_yyyy_xyyyyz, tr_yyyy_xyyyzz, tr_yyyy_xyyz, tr_yyyy_xyyzzz, tr_yyyy_xyzz, tr_yyyy_xyzzzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyyy_xxx, tr_yyyyy_xxxxx, tr_yyyyy_xxxxy, tr_yyyyy_xxxxz, tr_yyyyy_xxxyy, tr_yyyyy_xxxyz, tr_yyyyy_xxxzz, tr_yyyyy_xxy, tr_yyyyy_xxyyy, tr_yyyyy_xxyyz, tr_yyyyy_xxyzz, tr_yyyyy_xxz, tr_yyyyy_xxzzz, tr_yyyyy_xyy, tr_yyyyy_xyyyy, tr_yyyyy_xyyyz, tr_yyyyy_xyyzz, tr_yyyyy_xyz, tr_yyyyy_xyzzz, tr_yyyyy_xzz, tr_yyyyy_xzzzz, tr_yyyyy_yyy, tr_yyyyy_yyz, tr_yyyyy_yzz, tr_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyyy_xxxx[i] = 16.0 * tr_yyy_xxx[i] - 8.0 * tr_yyy_xxxxx[i] * tke_0 - 8.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yyyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyyy_xxxxx[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xxxy[i] = 12.0 * tr_yyy_xxy[i] - 8.0 * tr_yyy_xxxxy[i] * tke_0 + 3.0 * tr_yyyy_xx[i] - 6.0 * tr_yyyy_xxyy[i] * tke_0 - 2.0 * tr_yyyy_xxxx[i] * tke_0 + 4.0 * tr_yyyy_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yyyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyyy_xxxxy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxxy[i] * tbe_0 - 2.0 * tr_xyyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xxxz[i] = 12.0 * tr_yyy_xxz[i] - 8.0 * tr_yyy_xxxxz[i] * tke_0 - 6.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyyy_xxxxz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xxyy[i] = 8.0 * tr_yyy_xyy[i] - 8.0 * tr_yyy_xxxyy[i] * tke_0 + 4.0 * tr_yyyy_xy[i] - 4.0 * tr_yyyy_xyyy[i] * tke_0 - 4.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyyy_xxxyy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxyy[i] * tbe_0 - 4.0 * tr_xyyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xxyz[i] = 8.0 * tr_yyy_xyz[i] - 8.0 * tr_yyy_xxxyz[i] * tke_0 + 2.0 * tr_yyyy_xz[i] - 4.0 * tr_yyyy_xyyz[i] * tke_0 - 2.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_xyz[i] * tbe_0 + 4.0 * tr_yyyyy_xxxyz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xxzz[i] = 8.0 * tr_yyy_xzz[i] - 8.0 * tr_yyy_xxxzz[i] * tke_0 - 4.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyyy_xxxzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxzz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xyyy[i] = 4.0 * tr_yyy_yyy[i] - 8.0 * tr_yyy_xxyyy[i] * tke_0 + 3.0 * tr_yyyy_yy[i] - 2.0 * tr_yyyy_yyyy[i] * tke_0 - 6.0 * tr_yyyy_xxyy[i] * tke_0 + 4.0 * tr_yyyy_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyyy_xxyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xyyy[i] * tbe_0 - 6.0 * tr_xyyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xyyz[i] = 4.0 * tr_yyy_yyz[i] - 8.0 * tr_yyy_xxyyz[i] * tke_0 + 2.0 * tr_yyyy_yz[i] - 2.0 * tr_yyyy_yyyz[i] * tke_0 - 4.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyy_yyz[i] * tbe_0 + 4.0 * tr_yyyyy_xxyyz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xyyz[i] * tbe_0 - 4.0 * tr_xyyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xyzz[i] = 4.0 * tr_yyy_yzz[i] - 8.0 * tr_yyy_xxyzz[i] * tke_0 + tr_yyyy_zz[i] - 2.0 * tr_yyyy_yyzz[i] * tke_0 - 2.0 * tr_yyyy_xxzz[i] * tke_0 + 4.0 * tr_yyyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyy_yzz[i] * tbe_0 + 4.0 * tr_yyyyy_xxyzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xyzz[i] * tbe_0 - 2.0 * tr_xyyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xzzz[i] = 4.0 * tr_yyy_zzz[i] - 8.0 * tr_yyy_xxzzz[i] * tke_0 - 2.0 * tr_yyyy_yzzz[i] * tke_0 + 4.0 * tr_yyyy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyyy_xxzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xzzz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_yyyy[i] = -8.0 * tr_yyy_xyyyy[i] * tke_0 - 8.0 * tr_yyyy_xyyy[i] * tke_0 + 4.0 * tr_yyyy_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yyyy[i] * tbe_0 - 8.0 * tr_xyyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_yyyz[i] = -8.0 * tr_yyy_xyyyz[i] * tke_0 - 6.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xyyyz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yyyz[i] * tbe_0 - 6.0 * tr_xyyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_yyzz[i] = -8.0 * tr_yyy_xyyzz[i] * tke_0 - 4.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xyyzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_yzzz[i] = -8.0 * tr_yyy_xyzzz[i] * tke_0 - 2.0 * tr_yyyy_xzzz[i] * tke_0 + 4.0 * tr_yyyy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xyzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yzzz[i] * tbe_0 - 2.0 * tr_xyyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_zzzz[i] = -8.0 * tr_yyy_xzzzz[i] * tke_0 + 4.0 * tr_yyyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_zzzz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-405 components of targeted buffer : GG

    auto tr_0_0_xy_yyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 390);

    auto tr_0_0_xy_yyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 391);

    auto tr_0_0_xy_yyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 392);

    auto tr_0_0_xy_yyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 393);

    auto tr_0_0_xy_yyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 394);

    auto tr_0_0_xy_yyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 395);

    auto tr_0_0_xy_yyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 396);

    auto tr_0_0_xy_yyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 397);

    auto tr_0_0_xy_yyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 398);

    auto tr_0_0_xy_yyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 399);

    auto tr_0_0_xy_yyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 400);

    auto tr_0_0_xy_yyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 401);

    auto tr_0_0_xy_yyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 402);

    auto tr_0_0_xy_yyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 403);

    auto tr_0_0_xy_yyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 404);

    #pragma omp simd aligned(tr_0_0_xy_yyyz_xxxx, tr_0_0_xy_yyyz_xxxy, tr_0_0_xy_yyyz_xxxz, tr_0_0_xy_yyyz_xxyy, tr_0_0_xy_yyyz_xxyz, tr_0_0_xy_yyyz_xxzz, tr_0_0_xy_yyyz_xyyy, tr_0_0_xy_yyyz_xyyz, tr_0_0_xy_yyyz_xyzz, tr_0_0_xy_yyyz_xzzz, tr_0_0_xy_yyyz_yyyy, tr_0_0_xy_yyyz_yyyz, tr_0_0_xy_yyyz_yyzz, tr_0_0_xy_yyyz_yzzz, tr_0_0_xy_yyyz_zzzz, tr_xyyyyz_xxxx, tr_xyyyyz_xxxy, tr_xyyyyz_xxxz, tr_xyyyyz_xxyy, tr_xyyyyz_xxyz, tr_xyyyyz_xxzz, tr_xyyyyz_xyyy, tr_xyyyyz_xyyz, tr_xyyyyz_xyzz, tr_xyyyyz_xzzz, tr_xyyyyz_yyyy, tr_xyyyyz_yyyz, tr_xyyyyz_yyzz, tr_xyyyyz_yzzz, tr_xyyyyz_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxy, tr_xyyyz_xxxyy, tr_xyyyz_xxxyz, tr_xyyyz_xxy, tr_xyyyz_xxyyy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyyyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyyyy, tr_xyyyz_yyyyz, tr_xyyyz_yyyzz, tr_xyyyz_yyz, tr_xyyyz_yyzzz, tr_xyyyz_yzz, tr_xyyyz_yzzzz, tr_xyyyz_zzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxxxx, tr_yyyyz_xxxxy, tr_yyyyz_xxxxz, tr_yyyyz_xxxyy, tr_yyyyz_xxxyz, tr_yyyyz_xxxzz, tr_yyyyz_xxy, tr_yyyyz_xxyyy, tr_yyyyz_xxyyz, tr_yyyyz_xxyzz, tr_yyyyz_xxz, tr_yyyyz_xxzzz, tr_yyyyz_xyy, tr_yyyyz_xyyyy, tr_yyyyz_xyyyz, tr_yyyyz_xyyzz, tr_yyyyz_xyz, tr_yyyyz_xyzzz, tr_yyyyz_xzz, tr_yyyyz_xzzzz, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxxxy, tr_yyyz_xxxxyy, tr_yyyz_xxxxyz, tr_yyyz_xxxy, tr_yyyz_xxxyyy, tr_yyyz_xxxyyz, tr_yyyz_xxxyzz, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyyyy, tr_yyyz_xxyyyz, tr_yyyz_xxyyzz, tr_yyyz_xxyz, tr_yyyz_xxyzzz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyyyy, tr_yyyz_xyyyyz, tr_yyyz_xyyyzz, tr_yyyz_xyyz, tr_yyyz_xyyzzz, tr_yyyz_xyzz, tr_yyyz_xyzzzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyyz_xxxx[i] = 12.0 * tr_yyz_xxx[i] - 6.0 * tr_yyz_xxxxx[i] * tke_0 - 8.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xxxy[i] = 9.0 * tr_yyz_xxy[i] - 6.0 * tr_yyz_xxxxy[i] * tke_0 + 3.0 * tr_yyyz_xx[i] - 6.0 * tr_yyyz_xxyy[i] * tke_0 - 2.0 * tr_yyyz_xxxx[i] * tke_0 + 4.0 * tr_yyyz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yyyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxxy[i] * tbe_0 - 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xxxz[i] = 9.0 * tr_yyz_xxz[i] - 6.0 * tr_yyz_xxxxz[i] * tke_0 - 6.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xxyy[i] = 6.0 * tr_yyz_xyy[i] - 6.0 * tr_yyz_xxxyy[i] * tke_0 + 4.0 * tr_yyyz_xy[i] - 4.0 * tr_yyyz_xyyy[i] * tke_0 - 4.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xxyz[i] = 6.0 * tr_yyz_xyz[i] - 6.0 * tr_yyz_xxxyz[i] * tke_0 + 2.0 * tr_yyyz_xz[i] - 4.0 * tr_yyyz_xyyz[i] * tke_0 - 2.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xxzz[i] = 6.0 * tr_yyz_xzz[i] - 6.0 * tr_yyz_xxxzz[i] * tke_0 - 4.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xyyy[i] = 3.0 * tr_yyz_yyy[i] - 6.0 * tr_yyz_xxyyy[i] * tke_0 + 3.0 * tr_yyyz_yy[i] - 2.0 * tr_yyyz_yyyy[i] * tke_0 - 6.0 * tr_yyyz_xxyy[i] * tke_0 + 4.0 * tr_yyyz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xyyy[i] * tbe_0 - 6.0 * tr_xyyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xyyz[i] = 3.0 * tr_yyz_yyz[i] - 6.0 * tr_yyz_xxyyz[i] * tke_0 + 2.0 * tr_yyyz_yz[i] - 2.0 * tr_yyyz_yyyz[i] * tke_0 - 4.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xyyz[i] * tbe_0 - 4.0 * tr_xyyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xyzz[i] = 3.0 * tr_yyz_yzz[i] - 6.0 * tr_yyz_xxyzz[i] * tke_0 + tr_yyyz_zz[i] - 2.0 * tr_yyyz_yyzz[i] * tke_0 - 2.0 * tr_yyyz_xxzz[i] * tke_0 + 4.0 * tr_yyyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xyzz[i] * tbe_0 - 2.0 * tr_xyyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xzzz[i] = 3.0 * tr_yyz_zzz[i] - 6.0 * tr_yyz_xxzzz[i] * tke_0 - 2.0 * tr_yyyz_yzzz[i] * tke_0 + 4.0 * tr_yyyz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyyz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xzzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_yyyy[i] = -6.0 * tr_yyz_xyyyy[i] * tke_0 - 8.0 * tr_yyyz_xyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yyyy[i] * tbe_0 - 8.0 * tr_xyyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_yyyz[i] = -6.0 * tr_yyz_xyyyz[i] * tke_0 - 6.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yyyz[i] * tbe_0 - 6.0 * tr_xyyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_yyzz[i] = -6.0 * tr_yyz_xyyzz[i] * tke_0 - 4.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_yzzz[i] = -6.0 * tr_yyz_xyzzz[i] * tke_0 - 2.0 * tr_yyyz_xzzz[i] * tke_0 + 4.0 * tr_yyyz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yzzz[i] * tbe_0 - 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_zzzz[i] = -6.0 * tr_yyz_xzzzz[i] * tke_0 + 4.0 * tr_yyyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_zzzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 405-420 components of targeted buffer : GG

    auto tr_0_0_xy_yyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 405);

    auto tr_0_0_xy_yyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 406);

    auto tr_0_0_xy_yyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 407);

    auto tr_0_0_xy_yyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 408);

    auto tr_0_0_xy_yyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 409);

    auto tr_0_0_xy_yyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 410);

    auto tr_0_0_xy_yyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 411);

    auto tr_0_0_xy_yyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 412);

    auto tr_0_0_xy_yyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 413);

    auto tr_0_0_xy_yyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 414);

    auto tr_0_0_xy_yyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 415);

    auto tr_0_0_xy_yyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 416);

    auto tr_0_0_xy_yyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 417);

    auto tr_0_0_xy_yyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 418);

    auto tr_0_0_xy_yyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 419);

    #pragma omp simd aligned(tr_0_0_xy_yyzz_xxxx, tr_0_0_xy_yyzz_xxxy, tr_0_0_xy_yyzz_xxxz, tr_0_0_xy_yyzz_xxyy, tr_0_0_xy_yyzz_xxyz, tr_0_0_xy_yyzz_xxzz, tr_0_0_xy_yyzz_xyyy, tr_0_0_xy_yyzz_xyyz, tr_0_0_xy_yyzz_xyzz, tr_0_0_xy_yyzz_xzzz, tr_0_0_xy_yyzz_yyyy, tr_0_0_xy_yyzz_yyyz, tr_0_0_xy_yyzz_yyzz, tr_0_0_xy_yyzz_yzzz, tr_0_0_xy_yyzz_zzzz, tr_xyyyzz_xxxx, tr_xyyyzz_xxxy, tr_xyyyzz_xxxz, tr_xyyyzz_xxyy, tr_xyyyzz_xxyz, tr_xyyyzz_xxzz, tr_xyyyzz_xyyy, tr_xyyyzz_xyyz, tr_xyyyzz_xyzz, tr_xyyyzz_xzzz, tr_xyyyzz_yyyy, tr_xyyyzz_yyyz, tr_xyyyzz_yyzz, tr_xyyyzz_yzzz, tr_xyyyzz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxy, tr_xyyzz_xxxyy, tr_xyyzz_xxxyz, tr_xyyzz_xxy, tr_xyyzz_xxyyy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyyyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyyyy, tr_xyyzz_yyyyz, tr_xyyzz_yyyzz, tr_xyyzz_yyz, tr_xyyzz_yyzzz, tr_xyyzz_yzz, tr_xyyzz_yzzzz, tr_xyyzz_zzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxxxx, tr_yyyzz_xxxxy, tr_yyyzz_xxxxz, tr_yyyzz_xxxyy, tr_yyyzz_xxxyz, tr_yyyzz_xxxzz, tr_yyyzz_xxy, tr_yyyzz_xxyyy, tr_yyyzz_xxyyz, tr_yyyzz_xxyzz, tr_yyyzz_xxz, tr_yyyzz_xxzzz, tr_yyyzz_xyy, tr_yyyzz_xyyyy, tr_yyyzz_xyyyz, tr_yyyzz_xyyzz, tr_yyyzz_xyz, tr_yyyzz_xyzzz, tr_yyyzz_xzz, tr_yyyzz_xzzzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxxxy, tr_yyzz_xxxxyy, tr_yyzz_xxxxyz, tr_yyzz_xxxy, tr_yyzz_xxxyyy, tr_yyzz_xxxyyz, tr_yyzz_xxxyzz, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyyyy, tr_yyzz_xxyyyz, tr_yyzz_xxyyzz, tr_yyzz_xxyz, tr_yyzz_xxyzzz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyyyy, tr_yyzz_xyyyyz, tr_yyzz_xyyyzz, tr_yyzz_xyyz, tr_yyzz_xyyzzz, tr_yyzz_xyzz, tr_yyzz_xyzzzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyzz_xxxx[i] = 8.0 * tr_yzz_xxx[i] - 4.0 * tr_yzz_xxxxx[i] * tke_0 - 8.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xxxy[i] = 6.0 * tr_yzz_xxy[i] - 4.0 * tr_yzz_xxxxy[i] * tke_0 + 3.0 * tr_yyzz_xx[i] - 6.0 * tr_yyzz_xxyy[i] * tke_0 - 2.0 * tr_yyzz_xxxx[i] * tke_0 + 4.0 * tr_yyzz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yyyzz_xxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxxy[i] * tbe_0 - 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xxxz[i] = 6.0 * tr_yzz_xxz[i] - 4.0 * tr_yzz_xxxxz[i] * tke_0 - 6.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyzz_xxz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xxyy[i] = 4.0 * tr_yzz_xyy[i] - 4.0 * tr_yzz_xxxyy[i] * tke_0 + 4.0 * tr_yyzz_xy[i] - 4.0 * tr_yyzz_xyyy[i] * tke_0 - 4.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xyy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxyy[i] * tbe_0 - 4.0 * tr_xyyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xxyz[i] = 4.0 * tr_yzz_xyz[i] - 4.0 * tr_yzz_xxxyz[i] * tke_0 + 2.0 * tr_yyzz_xz[i] - 4.0 * tr_yyzz_xyyz[i] * tke_0 - 2.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xyz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tbe_0 - 2.0 * tr_xyyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xxzz[i] = 4.0 * tr_yzz_xzz[i] - 4.0 * tr_yzz_xxxzz[i] * tke_0 - 4.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xyyy[i] = 2.0 * tr_yzz_yyy[i] - 4.0 * tr_yzz_xxyyy[i] * tke_0 + 3.0 * tr_yyzz_yy[i] - 2.0 * tr_yyzz_yyyy[i] * tke_0 - 6.0 * tr_yyzz_xxyy[i] * tke_0 + 4.0 * tr_yyzz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_yyy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xyyy[i] * tbe_0 - 6.0 * tr_xyyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xyyz[i] = 2.0 * tr_yzz_yyz[i] - 4.0 * tr_yzz_xxyyz[i] * tke_0 + 2.0 * tr_yyzz_yz[i] - 2.0 * tr_yyzz_yyyz[i] * tke_0 - 4.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_yyz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tbe_0 - 4.0 * tr_xyyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xyzz[i] = 2.0 * tr_yzz_yzz[i] - 4.0 * tr_yzz_xxyzz[i] * tke_0 + tr_yyzz_zz[i] - 2.0 * tr_yyzz_yyzz[i] * tke_0 - 2.0 * tr_yyzz_xxzz[i] * tke_0 + 4.0 * tr_yyzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_yzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xyzz[i] * tbe_0 - 2.0 * tr_xyyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xzzz[i] = 2.0 * tr_yzz_zzz[i] - 4.0 * tr_yzz_xxzzz[i] * tke_0 - 2.0 * tr_yyzz_yzzz[i] * tke_0 + 4.0 * tr_yyzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_zzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xzzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_yyyy[i] = -4.0 * tr_yzz_xyyyy[i] * tke_0 - 8.0 * tr_yyzz_xyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yyyy[i] * tbe_0 - 8.0 * tr_xyyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_yyyz[i] = -4.0 * tr_yzz_xyyyz[i] * tke_0 - 6.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yyyz[i] * tbe_0 - 6.0 * tr_xyyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_yyzz[i] = -4.0 * tr_yzz_xyyzz[i] * tke_0 - 4.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yyzz[i] * tbe_0 - 4.0 * tr_xyyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_yzzz[i] = -4.0 * tr_yzz_xyzzz[i] * tke_0 - 2.0 * tr_yyzz_xzzz[i] * tke_0 + 4.0 * tr_yyzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yzzz[i] * tbe_0 - 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_zzzz[i] = -4.0 * tr_yzz_xzzzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_zzzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-435 components of targeted buffer : GG

    auto tr_0_0_xy_yzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 420);

    auto tr_0_0_xy_yzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 421);

    auto tr_0_0_xy_yzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 422);

    auto tr_0_0_xy_yzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 423);

    auto tr_0_0_xy_yzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 424);

    auto tr_0_0_xy_yzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 425);

    auto tr_0_0_xy_yzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 426);

    auto tr_0_0_xy_yzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 427);

    auto tr_0_0_xy_yzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 428);

    auto tr_0_0_xy_yzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 429);

    auto tr_0_0_xy_yzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 430);

    auto tr_0_0_xy_yzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 431);

    auto tr_0_0_xy_yzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 432);

    auto tr_0_0_xy_yzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 433);

    auto tr_0_0_xy_yzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 434);

    #pragma omp simd aligned(tr_0_0_xy_yzzz_xxxx, tr_0_0_xy_yzzz_xxxy, tr_0_0_xy_yzzz_xxxz, tr_0_0_xy_yzzz_xxyy, tr_0_0_xy_yzzz_xxyz, tr_0_0_xy_yzzz_xxzz, tr_0_0_xy_yzzz_xyyy, tr_0_0_xy_yzzz_xyyz, tr_0_0_xy_yzzz_xyzz, tr_0_0_xy_yzzz_xzzz, tr_0_0_xy_yzzz_yyyy, tr_0_0_xy_yzzz_yyyz, tr_0_0_xy_yzzz_yyzz, tr_0_0_xy_yzzz_yzzz, tr_0_0_xy_yzzz_zzzz, tr_xyyzzz_xxxx, tr_xyyzzz_xxxy, tr_xyyzzz_xxxz, tr_xyyzzz_xxyy, tr_xyyzzz_xxyz, tr_xyyzzz_xxzz, tr_xyyzzz_xyyy, tr_xyyzzz_xyyz, tr_xyyzzz_xyzz, tr_xyyzzz_xzzz, tr_xyyzzz_yyyy, tr_xyyzzz_yyyz, tr_xyyzzz_yyzz, tr_xyyzzz_yzzz, tr_xyyzzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxy, tr_xyzzz_xxxyy, tr_xyzzz_xxxyz, tr_xyzzz_xxy, tr_xyzzz_xxyyy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyyyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyyyy, tr_xyzzz_yyyyz, tr_xyzzz_yyyzz, tr_xyzzz_yyz, tr_xyzzz_yyzzz, tr_xyzzz_yzz, tr_xyzzz_yzzzz, tr_xyzzz_zzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxxxx, tr_yyzzz_xxxxy, tr_yyzzz_xxxxz, tr_yyzzz_xxxyy, tr_yyzzz_xxxyz, tr_yyzzz_xxxzz, tr_yyzzz_xxy, tr_yyzzz_xxyyy, tr_yyzzz_xxyyz, tr_yyzzz_xxyzz, tr_yyzzz_xxz, tr_yyzzz_xxzzz, tr_yyzzz_xyy, tr_yyzzz_xyyyy, tr_yyzzz_xyyyz, tr_yyzzz_xyyzz, tr_yyzzz_xyz, tr_yyzzz_xyzzz, tr_yyzzz_xzz, tr_yyzzz_xzzzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxxxy, tr_yzzz_xxxxyy, tr_yzzz_xxxxyz, tr_yzzz_xxxy, tr_yzzz_xxxyyy, tr_yzzz_xxxyyz, tr_yzzz_xxxyzz, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyyyy, tr_yzzz_xxyyyz, tr_yzzz_xxyyzz, tr_yzzz_xxyz, tr_yzzz_xxyzzz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyyyy, tr_yzzz_xyyyyz, tr_yzzz_xyyyzz, tr_yzzz_xyyz, tr_yzzz_xyyzzz, tr_yzzz_xyzz, tr_yzzz_xyzzzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_zzz_xxx, tr_zzz_xxxxx, tr_zzz_xxxxy, tr_zzz_xxxxz, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yzzz_xxxx[i] = 4.0 * tr_zzz_xxx[i] - 2.0 * tr_zzz_xxxxx[i] * tke_0 - 8.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xxxy[i] = 3.0 * tr_zzz_xxy[i] - 2.0 * tr_zzz_xxxxy[i] * tke_0 + 3.0 * tr_yzzz_xx[i] - 6.0 * tr_yzzz_xxyy[i] * tke_0 - 2.0 * tr_yzzz_xxxx[i] * tke_0 + 4.0 * tr_yzzz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yyzzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxxy[i] * tbe_0 - 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xxxz[i] = 3.0 * tr_zzz_xxz[i] - 2.0 * tr_zzz_xxxxz[i] * tke_0 - 6.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyzzz_xxz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xxyy[i] = 2.0 * tr_zzz_xyy[i] - 2.0 * tr_zzz_xxxyy[i] * tke_0 + 4.0 * tr_yzzz_xy[i] - 4.0 * tr_yzzz_xyyy[i] * tke_0 - 4.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxyy[i] * tbe_0 - 4.0 * tr_xyzzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xxyz[i] = 2.0 * tr_zzz_xyz[i] - 2.0 * tr_zzz_xxxyz[i] * tke_0 + 2.0 * tr_yzzz_xz[i] - 4.0 * tr_yzzz_xyyz[i] * tke_0 - 2.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xyz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxyz[i] * tbe_0 - 2.0 * tr_xyzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xxzz[i] = 2.0 * tr_zzz_xzz[i] - 2.0 * tr_zzz_xxxzz[i] * tke_0 - 4.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xyyy[i] = tr_zzz_yyy[i] - 2.0 * tr_zzz_xxyyy[i] * tke_0 + 3.0 * tr_yzzz_yy[i] - 2.0 * tr_yzzz_yyyy[i] * tke_0 - 6.0 * tr_yzzz_xxyy[i] * tke_0 + 4.0 * tr_yzzz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xyyy[i] * tbe_0 - 6.0 * tr_xyzzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xyyz[i] = tr_zzz_yyz[i] - 2.0 * tr_zzz_xxyyz[i] * tke_0 + 2.0 * tr_yzzz_yz[i] - 2.0 * tr_yzzz_yyyz[i] * tke_0 - 4.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_yyz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xyyz[i] * tbe_0 - 4.0 * tr_xyzzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xyzz[i] = tr_zzz_yzz[i] - 2.0 * tr_zzz_xxyzz[i] * tke_0 + tr_yzzz_zz[i] - 2.0 * tr_yzzz_yyzz[i] * tke_0 - 2.0 * tr_yzzz_xxzz[i] * tke_0 + 4.0 * tr_yzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_yzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xyzz[i] * tbe_0 - 2.0 * tr_xyzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xzzz[i] = tr_zzz_zzz[i] - 2.0 * tr_zzz_xxzzz[i] * tke_0 - 2.0 * tr_yzzz_yzzz[i] * tke_0 + 4.0 * tr_yzzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_zzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xzzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_yyyy[i] = -2.0 * tr_zzz_xyyyy[i] * tke_0 - 8.0 * tr_yzzz_xyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yyyy[i] * tbe_0 - 8.0 * tr_xyzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_yyyz[i] = -2.0 * tr_zzz_xyyyz[i] * tke_0 - 6.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yyyz[i] * tbe_0 - 6.0 * tr_xyzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_yyzz[i] = -2.0 * tr_zzz_xyyzz[i] * tke_0 - 4.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yyzz[i] * tbe_0 - 4.0 * tr_xyzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_yzzz[i] = -2.0 * tr_zzz_xyzzz[i] * tke_0 - 2.0 * tr_yzzz_xzzz[i] * tke_0 + 4.0 * tr_yzzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yzzz[i] * tbe_0 - 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_zzzz[i] = -2.0 * tr_zzz_xzzzz[i] * tke_0 + 4.0 * tr_yzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_zzzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 435-450 components of targeted buffer : GG

    auto tr_0_0_xy_zzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 435);

    auto tr_0_0_xy_zzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 436);

    auto tr_0_0_xy_zzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 437);

    auto tr_0_0_xy_zzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 438);

    auto tr_0_0_xy_zzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 439);

    auto tr_0_0_xy_zzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 440);

    auto tr_0_0_xy_zzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 441);

    auto tr_0_0_xy_zzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 442);

    auto tr_0_0_xy_zzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 443);

    auto tr_0_0_xy_zzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 444);

    auto tr_0_0_xy_zzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 445);

    auto tr_0_0_xy_zzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 446);

    auto tr_0_0_xy_zzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 447);

    auto tr_0_0_xy_zzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 448);

    auto tr_0_0_xy_zzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 449);

    #pragma omp simd aligned(tr_0_0_xy_zzzz_xxxx, tr_0_0_xy_zzzz_xxxy, tr_0_0_xy_zzzz_xxxz, tr_0_0_xy_zzzz_xxyy, tr_0_0_xy_zzzz_xxyz, tr_0_0_xy_zzzz_xxzz, tr_0_0_xy_zzzz_xyyy, tr_0_0_xy_zzzz_xyyz, tr_0_0_xy_zzzz_xyzz, tr_0_0_xy_zzzz_xzzz, tr_0_0_xy_zzzz_yyyy, tr_0_0_xy_zzzz_yyyz, tr_0_0_xy_zzzz_yyzz, tr_0_0_xy_zzzz_yzzz, tr_0_0_xy_zzzz_zzzz, tr_xyzzzz_xxxx, tr_xyzzzz_xxxy, tr_xyzzzz_xxxz, tr_xyzzzz_xxyy, tr_xyzzzz_xxyz, tr_xyzzzz_xxzz, tr_xyzzzz_xyyy, tr_xyzzzz_xyyz, tr_xyzzzz_xyzz, tr_xyzzzz_xzzz, tr_xyzzzz_yyyy, tr_xyzzzz_yyyz, tr_xyzzzz_yyzz, tr_xyzzzz_yzzz, tr_xyzzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxxxy, tr_xzzzz_xxxyy, tr_xzzzz_xxxyz, tr_xzzzz_xxy, tr_xzzzz_xxyyy, tr_xzzzz_xxyyz, tr_xzzzz_xxyzz, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyyyy, tr_xzzzz_xyyyz, tr_xzzzz_xyyzz, tr_xzzzz_xyz, tr_xzzzz_xyzzz, tr_xzzzz_xzz, tr_xzzzz_yyy, tr_xzzzz_yyyyy, tr_xzzzz_yyyyz, tr_xzzzz_yyyzz, tr_xzzzz_yyz, tr_xzzzz_yyzzz, tr_xzzzz_yzz, tr_xzzzz_yzzzz, tr_xzzzz_zzz, tr_yzzzz_xxx, tr_yzzzz_xxxxx, tr_yzzzz_xxxxy, tr_yzzzz_xxxxz, tr_yzzzz_xxxyy, tr_yzzzz_xxxyz, tr_yzzzz_xxxzz, tr_yzzzz_xxy, tr_yzzzz_xxyyy, tr_yzzzz_xxyyz, tr_yzzzz_xxyzz, tr_yzzzz_xxz, tr_yzzzz_xxzzz, tr_yzzzz_xyy, tr_yzzzz_xyyyy, tr_yzzzz_xyyyz, tr_yzzzz_xyyzz, tr_yzzzz_xyz, tr_yzzzz_xyzzz, tr_yzzzz_xzz, tr_yzzzz_xzzzz, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_zzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxxxy, tr_zzzz_xxxxyy, tr_zzzz_xxxxyz, tr_zzzz_xxxy, tr_zzzz_xxxyyy, tr_zzzz_xxxyyz, tr_zzzz_xxxyzz, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyyyy, tr_zzzz_xxyyyz, tr_zzzz_xxyyzz, tr_zzzz_xxyz, tr_zzzz_xxyzzz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyyyy, tr_zzzz_xyyyyz, tr_zzzz_xyyyzz, tr_zzzz_xyyz, tr_zzzz_xyyzzz, tr_zzzz_xyzz, tr_zzzz_xyzzzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_zzzz_xxxx[i] = -8.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxxxy[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xxxy[i] = 3.0 * tr_zzzz_xx[i] - 6.0 * tr_zzzz_xxyy[i] * tke_0 - 2.0 * tr_zzzz_xxxx[i] * tke_0 + 4.0 * tr_zzzz_xxxxyy[i] * tke_0 * tke_0 - 6.0 * tr_yzzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xxxz[i] = -6.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yzzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xxyy[i] = 4.0 * tr_zzzz_xy[i] - 4.0 * tr_zzzz_xyyy[i] * tke_0 - 4.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xzzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xxyz[i] = 2.0 * tr_zzzz_xz[i] - 4.0 * tr_zzzz_xyyz[i] * tke_0 - 2.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xxzz[i] = -4.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xyyy[i] = 3.0 * tr_zzzz_yy[i] - 2.0 * tr_zzzz_yyyy[i] * tke_0 - 6.0 * tr_zzzz_xxyy[i] * tke_0 + 4.0 * tr_zzzz_xxyyyy[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xyyz[i] = 2.0 * tr_zzzz_yz[i] - 2.0 * tr_zzzz_yyyz[i] * tke_0 - 4.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xzzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xyzz[i] = tr_zzzz_zz[i] - 2.0 * tr_zzzz_yyzz[i] * tke_0 - 2.0 * tr_zzzz_xxzz[i] * tke_0 + 4.0 * tr_zzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xzzz[i] = -2.0 * tr_zzzz_yzzz[i] * tke_0 + 4.0 * tr_zzzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_yyyy[i] = -8.0 * tr_zzzz_xyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyyy[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_yyyz[i] = -6.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_yyzz[i] = -4.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xzzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_yzzz[i] = -2.0 * tr_zzzz_xzzz[i] * tke_0 + 4.0 * tr_zzzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_zzzz[i] = 4.0 * tr_zzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-465 components of targeted buffer : GG

    auto tr_0_0_xz_xxxx_xxxx = pbuffer.data(idx_op_geom_020_gg + 450);

    auto tr_0_0_xz_xxxx_xxxy = pbuffer.data(idx_op_geom_020_gg + 451);

    auto tr_0_0_xz_xxxx_xxxz = pbuffer.data(idx_op_geom_020_gg + 452);

    auto tr_0_0_xz_xxxx_xxyy = pbuffer.data(idx_op_geom_020_gg + 453);

    auto tr_0_0_xz_xxxx_xxyz = pbuffer.data(idx_op_geom_020_gg + 454);

    auto tr_0_0_xz_xxxx_xxzz = pbuffer.data(idx_op_geom_020_gg + 455);

    auto tr_0_0_xz_xxxx_xyyy = pbuffer.data(idx_op_geom_020_gg + 456);

    auto tr_0_0_xz_xxxx_xyyz = pbuffer.data(idx_op_geom_020_gg + 457);

    auto tr_0_0_xz_xxxx_xyzz = pbuffer.data(idx_op_geom_020_gg + 458);

    auto tr_0_0_xz_xxxx_xzzz = pbuffer.data(idx_op_geom_020_gg + 459);

    auto tr_0_0_xz_xxxx_yyyy = pbuffer.data(idx_op_geom_020_gg + 460);

    auto tr_0_0_xz_xxxx_yyyz = pbuffer.data(idx_op_geom_020_gg + 461);

    auto tr_0_0_xz_xxxx_yyzz = pbuffer.data(idx_op_geom_020_gg + 462);

    auto tr_0_0_xz_xxxx_yzzz = pbuffer.data(idx_op_geom_020_gg + 463);

    auto tr_0_0_xz_xxxx_zzzz = pbuffer.data(idx_op_geom_020_gg + 464);

    #pragma omp simd aligned(tr_0_0_xz_xxxx_xxxx, tr_0_0_xz_xxxx_xxxy, tr_0_0_xz_xxxx_xxxz, tr_0_0_xz_xxxx_xxyy, tr_0_0_xz_xxxx_xxyz, tr_0_0_xz_xxxx_xxzz, tr_0_0_xz_xxxx_xyyy, tr_0_0_xz_xxxx_xyyz, tr_0_0_xz_xxxx_xyzz, tr_0_0_xz_xxxx_xzzz, tr_0_0_xz_xxxx_yyyy, tr_0_0_xz_xxxx_yyyz, tr_0_0_xz_xxxx_yyzz, tr_0_0_xz_xxxx_yzzz, tr_0_0_xz_xxxx_zzzz, tr_xxx_xxx, tr_xxx_xxxxz, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxx_zzzzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxxxz, tr_xxxx_xxxxyz, tr_xxxx_xxxxzz, tr_xxxx_xxxy, tr_xxxx_xxxyyz, tr_xxxx_xxxyzz, tr_xxxx_xxxz, tr_xxxx_xxxzzz, tr_xxxx_xxyy, tr_xxxx_xxyyyz, tr_xxxx_xxyyzz, tr_xxxx_xxyz, tr_xxxx_xxyzzz, tr_xxxx_xxzz, tr_xxxx_xxzzzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyyyz, tr_xxxx_xyyyzz, tr_xxxx_xyyz, tr_xxxx_xyyzzz, tr_xxxx_xyzz, tr_xxxx_xyzzzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_xzzzzz, tr_xxxx_yy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxx_xxx, tr_xxxxx_xxxxz, tr_xxxxx_xxxyz, tr_xxxxx_xxxzz, tr_xxxxx_xxy, tr_xxxxx_xxyyz, tr_xxxxx_xxyzz, tr_xxxxx_xxz, tr_xxxxx_xxzzz, tr_xxxxx_xyy, tr_xxxxx_xyyyz, tr_xxxxx_xyyzz, tr_xxxxx_xyz, tr_xxxxx_xyzzz, tr_xxxxx_xzz, tr_xxxxx_xzzzz, tr_xxxxx_yyy, tr_xxxxx_yyyyz, tr_xxxxx_yyyzz, tr_xxxxx_yyz, tr_xxxxx_yyzzz, tr_xxxxx_yzz, tr_xxxxx_yzzzz, tr_xxxxx_zzz, tr_xxxxx_zzzzz, tr_xxxxxz_xxxx, tr_xxxxxz_xxxy, tr_xxxxxz_xxxz, tr_xxxxxz_xxyy, tr_xxxxxz_xxyz, tr_xxxxxz_xxzz, tr_xxxxxz_xyyy, tr_xxxxxz_xyyz, tr_xxxxxz_xyzz, tr_xxxxxz_xzzz, tr_xxxxxz_yyyy, tr_xxxxxz_yyyz, tr_xxxxxz_yyzz, tr_xxxxxz_yzzz, tr_xxxxxz_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxxxx, tr_xxxxz_xxxxy, tr_xxxxz_xxxxz, tr_xxxxz_xxxyy, tr_xxxxz_xxxyz, tr_xxxxz_xxxzz, tr_xxxxz_xxy, tr_xxxxz_xxyyy, tr_xxxxz_xxyyz, tr_xxxxz_xxyzz, tr_xxxxz_xxz, tr_xxxxz_xxzzz, tr_xxxxz_xyy, tr_xxxxz_xyyyy, tr_xxxxz_xyyyz, tr_xxxxz_xyyzz, tr_xxxxz_xyz, tr_xxxxz_xyzzz, tr_xxxxz_xzz, tr_xxxxz_xzzzz, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_zzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxxx_xxxx[i] = -8.0 * tr_xxx_xxxxz[i] * tke_0 - 8.0 * tr_xxxz_xxxx[i] * tbe_0 - 8.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xxxy[i] = -8.0 * tr_xxx_xxxyz[i] * tke_0 - 8.0 * tr_xxxz_xxxy[i] * tbe_0 - 6.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xxxz[i] = 4.0 * tr_xxx_xxx[i] - 8.0 * tr_xxx_xxxzz[i] * tke_0 - 8.0 * tr_xxxz_xxxz[i] * tbe_0 + 3.0 * tr_xxxx_xx[i] - 6.0 * tr_xxxx_xxzz[i] * tke_0 - 2.0 * tr_xxxx_xxxx[i] * tke_0 + 4.0 * tr_xxxx_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xxxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxxz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxx_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xxyy[i] = -8.0 * tr_xxx_xxyyz[i] * tke_0 - 8.0 * tr_xxxz_xxyy[i] * tbe_0 - 4.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xxyz[i] = 4.0 * tr_xxx_xxy[i] - 8.0 * tr_xxx_xxyzz[i] * tke_0 - 8.0 * tr_xxxz_xxyz[i] * tbe_0 + 2.0 * tr_xxxx_xy[i] - 4.0 * tr_xxxx_xyzz[i] * tke_0 - 2.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxxz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxxx_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xxzz[i] = 8.0 * tr_xxx_xxz[i] - 8.0 * tr_xxx_xxzzz[i] * tke_0 - 8.0 * tr_xxxz_xxzz[i] * tbe_0 + 4.0 * tr_xxxx_xz[i] - 4.0 * tr_xxxx_xzzz[i] * tke_0 - 4.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxxz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxxx_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xyyy[i] = -8.0 * tr_xxx_xyyyz[i] * tke_0 - 8.0 * tr_xxxz_xyyy[i] * tbe_0 - 2.0 * tr_xxxx_yyyz[i] * tke_0 + 4.0 * tr_xxxx_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xyyz[i] = 4.0 * tr_xxx_xyy[i] - 8.0 * tr_xxx_xyyzz[i] * tke_0 - 8.0 * tr_xxxz_xyyz[i] * tbe_0 + tr_xxxx_yy[i] - 2.0 * tr_xxxx_yyzz[i] * tke_0 - 2.0 * tr_xxxx_xxyy[i] * tke_0 + 4.0 * tr_xxxx_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxxx_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xyzz[i] = 8.0 * tr_xxx_xyz[i] - 8.0 * tr_xxx_xyzzz[i] * tke_0 - 8.0 * tr_xxxz_xyzz[i] * tbe_0 + 2.0 * tr_xxxx_yz[i] - 2.0 * tr_xxxx_yzzz[i] * tke_0 - 4.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxxx_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xzzz[i] = 12.0 * tr_xxx_xzz[i] - 8.0 * tr_xxx_xzzzz[i] * tke_0 - 8.0 * tr_xxxz_xzzz[i] * tbe_0 + 3.0 * tr_xxxx_zz[i] - 2.0 * tr_xxxx_zzzz[i] * tke_0 - 6.0 * tr_xxxx_xxzz[i] * tke_0 + 4.0 * tr_xxxx_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxxx_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_yyyy[i] = -8.0 * tr_xxx_yyyyz[i] * tke_0 - 8.0 * tr_xxxz_yyyy[i] * tbe_0 + 4.0 * tr_xxxx_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_yyyz[i] = 4.0 * tr_xxx_yyy[i] - 8.0 * tr_xxx_yyyzz[i] * tke_0 - 8.0 * tr_xxxz_yyyz[i] * tbe_0 - 2.0 * tr_xxxx_xyyy[i] * tke_0 + 4.0 * tr_xxxx_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxx_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_yyzz[i] = 8.0 * tr_xxx_yyz[i] - 8.0 * tr_xxx_yyzzz[i] * tke_0 - 8.0 * tr_xxxz_yyzz[i] * tbe_0 - 4.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxxx_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_yzzz[i] = 12.0 * tr_xxx_yzz[i] - 8.0 * tr_xxx_yzzzz[i] * tke_0 - 8.0 * tr_xxxz_yzzz[i] * tbe_0 - 6.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxxx_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_zzzz[i] = 16.0 * tr_xxx_zzz[i] - 8.0 * tr_xxx_zzzzz[i] * tke_0 - 8.0 * tr_xxxz_zzzz[i] * tbe_0 - 8.0 * tr_xxxx_xzzz[i] * tke_0 + 4.0 * tr_xxxx_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxx_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 465-480 components of targeted buffer : GG

    auto tr_0_0_xz_xxxy_xxxx = pbuffer.data(idx_op_geom_020_gg + 465);

    auto tr_0_0_xz_xxxy_xxxy = pbuffer.data(idx_op_geom_020_gg + 466);

    auto tr_0_0_xz_xxxy_xxxz = pbuffer.data(idx_op_geom_020_gg + 467);

    auto tr_0_0_xz_xxxy_xxyy = pbuffer.data(idx_op_geom_020_gg + 468);

    auto tr_0_0_xz_xxxy_xxyz = pbuffer.data(idx_op_geom_020_gg + 469);

    auto tr_0_0_xz_xxxy_xxzz = pbuffer.data(idx_op_geom_020_gg + 470);

    auto tr_0_0_xz_xxxy_xyyy = pbuffer.data(idx_op_geom_020_gg + 471);

    auto tr_0_0_xz_xxxy_xyyz = pbuffer.data(idx_op_geom_020_gg + 472);

    auto tr_0_0_xz_xxxy_xyzz = pbuffer.data(idx_op_geom_020_gg + 473);

    auto tr_0_0_xz_xxxy_xzzz = pbuffer.data(idx_op_geom_020_gg + 474);

    auto tr_0_0_xz_xxxy_yyyy = pbuffer.data(idx_op_geom_020_gg + 475);

    auto tr_0_0_xz_xxxy_yyyz = pbuffer.data(idx_op_geom_020_gg + 476);

    auto tr_0_0_xz_xxxy_yyzz = pbuffer.data(idx_op_geom_020_gg + 477);

    auto tr_0_0_xz_xxxy_yzzz = pbuffer.data(idx_op_geom_020_gg + 478);

    auto tr_0_0_xz_xxxy_zzzz = pbuffer.data(idx_op_geom_020_gg + 479);

    #pragma omp simd aligned(tr_0_0_xz_xxxy_xxxx, tr_0_0_xz_xxxy_xxxy, tr_0_0_xz_xxxy_xxxz, tr_0_0_xz_xxxy_xxyy, tr_0_0_xz_xxxy_xxyz, tr_0_0_xz_xxxy_xxzz, tr_0_0_xz_xxxy_xyyy, tr_0_0_xz_xxxy_xyyz, tr_0_0_xz_xxxy_xyzz, tr_0_0_xz_xxxy_xzzz, tr_0_0_xz_xxxy_yyyy, tr_0_0_xz_xxxy_yyyz, tr_0_0_xz_xxxy_yyzz, tr_0_0_xz_xxxy_yzzz, tr_0_0_xz_xxxy_zzzz, tr_xxxxy_xxx, tr_xxxxy_xxxxz, tr_xxxxy_xxxyz, tr_xxxxy_xxxzz, tr_xxxxy_xxy, tr_xxxxy_xxyyz, tr_xxxxy_xxyzz, tr_xxxxy_xxz, tr_xxxxy_xxzzz, tr_xxxxy_xyy, tr_xxxxy_xyyyz, tr_xxxxy_xyyzz, tr_xxxxy_xyz, tr_xxxxy_xyzzz, tr_xxxxy_xzz, tr_xxxxy_xzzzz, tr_xxxxy_yyy, tr_xxxxy_yyyyz, tr_xxxxy_yyyzz, tr_xxxxy_yyz, tr_xxxxy_yyzzz, tr_xxxxy_yzz, tr_xxxxy_yzzzz, tr_xxxxy_zzz, tr_xxxxy_zzzzz, tr_xxxxyz_xxxx, tr_xxxxyz_xxxy, tr_xxxxyz_xxxz, tr_xxxxyz_xxyy, tr_xxxxyz_xxyz, tr_xxxxyz_xxzz, tr_xxxxyz_xyyy, tr_xxxxyz_xyyz, tr_xxxxyz_xyzz, tr_xxxxyz_xzzz, tr_xxxxyz_yyyy, tr_xxxxyz_yyyz, tr_xxxxyz_yyzz, tr_xxxxyz_yzzz, tr_xxxxyz_zzzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxxxz, tr_xxxy_xxxxyz, tr_xxxy_xxxxzz, tr_xxxy_xxxy, tr_xxxy_xxxyyz, tr_xxxy_xxxyzz, tr_xxxy_xxxz, tr_xxxy_xxxzzz, tr_xxxy_xxyy, tr_xxxy_xxyyyz, tr_xxxy_xxyyzz, tr_xxxy_xxyz, tr_xxxy_xxyzzz, tr_xxxy_xxzz, tr_xxxy_xxzzzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyyyz, tr_xxxy_xyyyzz, tr_xxxy_xyyz, tr_xxxy_xyyzzz, tr_xxxy_xyzz, tr_xxxy_xyzzzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_xzzzzz, tr_xxxy_yy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxx, tr_xxxyz_xxxxy, tr_xxxyz_xxxxz, tr_xxxyz_xxxyy, tr_xxxyz_xxxyz, tr_xxxyz_xxxzz, tr_xxxyz_xxy, tr_xxxyz_xxyyy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xxzzz, tr_xxxyz_xyy, tr_xxxyz_xyyyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_xzzzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxxy_xxxx[i] = -6.0 * tr_xxy_xxxxz[i] * tke_0 - 6.0 * tr_xxyz_xxxx[i] * tbe_0 - 8.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xxxy[i] = -6.0 * tr_xxy_xxxyz[i] * tke_0 - 6.0 * tr_xxyz_xxxy[i] * tbe_0 - 6.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xxxz[i] = 3.0 * tr_xxy_xxx[i] - 6.0 * tr_xxy_xxxzz[i] * tke_0 - 6.0 * tr_xxyz_xxxz[i] * tbe_0 + 3.0 * tr_xxxy_xx[i] - 6.0 * tr_xxxy_xxzz[i] * tke_0 - 2.0 * tr_xxxy_xxxx[i] * tke_0 + 4.0 * tr_xxxy_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xxxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xxyy[i] = -6.0 * tr_xxy_xxyyz[i] * tke_0 - 6.0 * tr_xxyz_xxyy[i] * tbe_0 - 4.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xxyz[i] = 3.0 * tr_xxy_xxy[i] - 6.0 * tr_xxy_xxyzz[i] * tke_0 - 6.0 * tr_xxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxxy_xy[i] - 4.0 * tr_xxxy_xyzz[i] * tke_0 - 2.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xyz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xxzz[i] = 6.0 * tr_xxy_xxz[i] - 6.0 * tr_xxy_xxzzz[i] * tke_0 - 6.0 * tr_xxyz_xxzz[i] * tbe_0 + 4.0 * tr_xxxy_xz[i] - 4.0 * tr_xxxy_xzzz[i] * tke_0 - 4.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xyyy[i] = -6.0 * tr_xxy_xyyyz[i] * tke_0 - 6.0 * tr_xxyz_xyyy[i] * tbe_0 - 2.0 * tr_xxxy_yyyz[i] * tke_0 + 4.0 * tr_xxxy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xyyz[i] = 3.0 * tr_xxy_xyy[i] - 6.0 * tr_xxy_xyyzz[i] * tke_0 - 6.0 * tr_xxyz_xyyz[i] * tbe_0 + tr_xxxy_yy[i] - 2.0 * tr_xxxy_yyzz[i] * tke_0 - 2.0 * tr_xxxy_xxyy[i] * tke_0 + 4.0 * tr_xxxy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xyzz[i] = 6.0 * tr_xxy_xyz[i] - 6.0 * tr_xxy_xyzzz[i] * tke_0 - 6.0 * tr_xxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxxy_yz[i] - 2.0 * tr_xxxy_yzzz[i] * tke_0 - 4.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xzzz[i] = 9.0 * tr_xxy_xzz[i] - 6.0 * tr_xxy_xzzzz[i] * tke_0 - 6.0 * tr_xxyz_xzzz[i] * tbe_0 + 3.0 * tr_xxxy_zz[i] - 2.0 * tr_xxxy_zzzz[i] * tke_0 - 6.0 * tr_xxxy_xxzz[i] * tke_0 + 4.0 * tr_xxxy_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_yyyy[i] = -6.0 * tr_xxy_yyyyz[i] * tke_0 - 6.0 * tr_xxyz_yyyy[i] * tbe_0 + 4.0 * tr_xxxy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_yyyz[i] = 3.0 * tr_xxy_yyy[i] - 6.0 * tr_xxy_yyyzz[i] * tke_0 - 6.0 * tr_xxyz_yyyz[i] * tbe_0 - 2.0 * tr_xxxy_xyyy[i] * tke_0 + 4.0 * tr_xxxy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_yyzz[i] = 6.0 * tr_xxy_yyz[i] - 6.0 * tr_xxy_yyzzz[i] * tke_0 - 6.0 * tr_xxyz_yyzz[i] * tbe_0 - 4.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_yzzz[i] = 9.0 * tr_xxy_yzz[i] - 6.0 * tr_xxy_yzzzz[i] * tke_0 - 6.0 * tr_xxyz_yzzz[i] * tbe_0 - 6.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_zzzz[i] = 12.0 * tr_xxy_zzz[i] - 6.0 * tr_xxy_zzzzz[i] * tke_0 - 6.0 * tr_xxyz_zzzz[i] * tbe_0 - 8.0 * tr_xxxy_xzzz[i] * tke_0 + 4.0 * tr_xxxy_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-495 components of targeted buffer : GG

    auto tr_0_0_xz_xxxz_xxxx = pbuffer.data(idx_op_geom_020_gg + 480);

    auto tr_0_0_xz_xxxz_xxxy = pbuffer.data(idx_op_geom_020_gg + 481);

    auto tr_0_0_xz_xxxz_xxxz = pbuffer.data(idx_op_geom_020_gg + 482);

    auto tr_0_0_xz_xxxz_xxyy = pbuffer.data(idx_op_geom_020_gg + 483);

    auto tr_0_0_xz_xxxz_xxyz = pbuffer.data(idx_op_geom_020_gg + 484);

    auto tr_0_0_xz_xxxz_xxzz = pbuffer.data(idx_op_geom_020_gg + 485);

    auto tr_0_0_xz_xxxz_xyyy = pbuffer.data(idx_op_geom_020_gg + 486);

    auto tr_0_0_xz_xxxz_xyyz = pbuffer.data(idx_op_geom_020_gg + 487);

    auto tr_0_0_xz_xxxz_xyzz = pbuffer.data(idx_op_geom_020_gg + 488);

    auto tr_0_0_xz_xxxz_xzzz = pbuffer.data(idx_op_geom_020_gg + 489);

    auto tr_0_0_xz_xxxz_yyyy = pbuffer.data(idx_op_geom_020_gg + 490);

    auto tr_0_0_xz_xxxz_yyyz = pbuffer.data(idx_op_geom_020_gg + 491);

    auto tr_0_0_xz_xxxz_yyzz = pbuffer.data(idx_op_geom_020_gg + 492);

    auto tr_0_0_xz_xxxz_yzzz = pbuffer.data(idx_op_geom_020_gg + 493);

    auto tr_0_0_xz_xxxz_zzzz = pbuffer.data(idx_op_geom_020_gg + 494);

    #pragma omp simd aligned(tr_0_0_xz_xxxz_xxxx, tr_0_0_xz_xxxz_xxxy, tr_0_0_xz_xxxz_xxxz, tr_0_0_xz_xxxz_xxyy, tr_0_0_xz_xxxz_xxyz, tr_0_0_xz_xxxz_xxzz, tr_0_0_xz_xxxz_xyyy, tr_0_0_xz_xxxz_xyyz, tr_0_0_xz_xxxz_xyzz, tr_0_0_xz_xxxz_xzzz, tr_0_0_xz_xxxz_yyyy, tr_0_0_xz_xxxz_yyyz, tr_0_0_xz_xxxz_yyzz, tr_0_0_xz_xxxz_yzzz, tr_0_0_xz_xxxz_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxx_xxx, tr_xxx_xxxxx, tr_xxx_xxxxy, tr_xxx_xxxxz, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xzzz, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yzzz, tr_xxxx_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxxxz, tr_xxxxz_xxxyz, tr_xxxxz_xxxzz, tr_xxxxz_xxy, tr_xxxxz_xxyyz, tr_xxxxz_xxyzz, tr_xxxxz_xxz, tr_xxxxz_xxzzz, tr_xxxxz_xyy, tr_xxxxz_xyyyz, tr_xxxxz_xyyzz, tr_xxxxz_xyz, tr_xxxxz_xyzzz, tr_xxxxz_xzz, tr_xxxxz_xzzzz, tr_xxxxz_yyy, tr_xxxxz_yyyyz, tr_xxxxz_yyyzz, tr_xxxxz_yyz, tr_xxxxz_yyzzz, tr_xxxxz_yzz, tr_xxxxz_yzzzz, tr_xxxxz_zzz, tr_xxxxz_zzzzz, tr_xxxxzz_xxxx, tr_xxxxzz_xxxy, tr_xxxxzz_xxxz, tr_xxxxzz_xxyy, tr_xxxxzz_xxyz, tr_xxxxzz_xxzz, tr_xxxxzz_xyyy, tr_xxxxzz_xyyz, tr_xxxxzz_xyzz, tr_xxxxzz_xzzz, tr_xxxxzz_yyyy, tr_xxxxzz_yyyz, tr_xxxxzz_yyzz, tr_xxxxzz_yzzz, tr_xxxxzz_zzzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxxxz, tr_xxxz_xxxxyz, tr_xxxz_xxxxzz, tr_xxxz_xxxy, tr_xxxz_xxxyyz, tr_xxxz_xxxyzz, tr_xxxz_xxxz, tr_xxxz_xxxzzz, tr_xxxz_xxyy, tr_xxxz_xxyyyz, tr_xxxz_xxyyzz, tr_xxxz_xxyz, tr_xxxz_xxyzzz, tr_xxxz_xxzz, tr_xxxz_xxzzzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyyyz, tr_xxxz_xyyyzz, tr_xxxz_xyyz, tr_xxxz_xyyzzz, tr_xxxz_xyzz, tr_xxxz_xyzzzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_xzzzzz, tr_xxxz_yy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxxxx, tr_xxxzz_xxxxy, tr_xxxzz_xxxxz, tr_xxxzz_xxxyy, tr_xxxzz_xxxyz, tr_xxxzz_xxxzz, tr_xxxzz_xxy, tr_xxxzz_xxyyy, tr_xxxzz_xxyyz, tr_xxxzz_xxyzz, tr_xxxzz_xxz, tr_xxxzz_xxzzz, tr_xxxzz_xyy, tr_xxxzz_xyyyy, tr_xxxzz_xyyyz, tr_xxxzz_xyyzz, tr_xxxzz_xyz, tr_xxxzz_xyzzz, tr_xxxzz_xzz, tr_xxxzz_xzzzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxxz_xxxx[i] = 3.0 * tr_xx_xxxx[i] - 6.0 * tr_xxz_xxxxz[i] * tke_0 - 6.0 * tr_xxzz_xxxx[i] * tbe_0 + 4.0 * tr_xxx_xxx[i] - 2.0 * tr_xxx_xxxxx[i] * tke_0 - 8.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xxxy[i] = 3.0 * tr_xx_xxxy[i] - 6.0 * tr_xxz_xxxyz[i] * tke_0 - 6.0 * tr_xxzz_xxxy[i] * tbe_0 + 3.0 * tr_xxx_xxy[i] - 2.0 * tr_xxx_xxxxy[i] * tke_0 - 6.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxzz_xxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xxxz[i] = 3.0 * tr_xx_xxxz[i] + 3.0 * tr_xxz_xxx[i] - 6.0 * tr_xxz_xxxzz[i] * tke_0 - 6.0 * tr_xxzz_xxxz[i] * tbe_0 + 3.0 * tr_xxx_xxz[i] - 2.0 * tr_xxx_xxxxz[i] * tke_0 + 3.0 * tr_xxxz_xx[i] - 6.0 * tr_xxxz_xxzz[i] * tke_0 - 2.0 * tr_xxxz_xxxx[i] * tke_0 + 4.0 * tr_xxxz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xxxzz_xxz[i] * tbe_0 + 4.0 * tr_xxxzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxxz[i] * tbe_0 - 2.0 * tr_xxxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xxyy[i] = 3.0 * tr_xx_xxyy[i] - 6.0 * tr_xxz_xxyyz[i] * tke_0 - 6.0 * tr_xxzz_xxyy[i] * tbe_0 + 2.0 * tr_xxx_xyy[i] - 2.0 * tr_xxx_xxxyy[i] * tke_0 - 4.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xyy[i] * tbe_0 + 4.0 * tr_xxxzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxyy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xxyz[i] = 3.0 * tr_xx_xxyz[i] + 3.0 * tr_xxz_xxy[i] - 6.0 * tr_xxz_xxyzz[i] * tke_0 - 6.0 * tr_xxzz_xxyz[i] * tbe_0 + 2.0 * tr_xxx_xyz[i] - 2.0 * tr_xxx_xxxyz[i] * tke_0 + 2.0 * tr_xxxz_xy[i] - 4.0 * tr_xxxz_xyzz[i] * tke_0 - 2.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xyz[i] * tbe_0 + 4.0 * tr_xxxzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxyz[i] * tbe_0 - 2.0 * tr_xxxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xxzz[i] = 3.0 * tr_xx_xxzz[i] + 6.0 * tr_xxz_xxz[i] - 6.0 * tr_xxz_xxzzz[i] * tke_0 - 6.0 * tr_xxzz_xxzz[i] * tbe_0 + 2.0 * tr_xxx_xzz[i] - 2.0 * tr_xxx_xxxzz[i] * tke_0 + 4.0 * tr_xxxz_xz[i] - 4.0 * tr_xxxz_xzzz[i] * tke_0 - 4.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xzz[i] * tbe_0 + 4.0 * tr_xxxzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxzz[i] * tbe_0 - 4.0 * tr_xxxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xyyy[i] = 3.0 * tr_xx_xyyy[i] - 6.0 * tr_xxz_xyyyz[i] * tke_0 - 6.0 * tr_xxzz_xyyy[i] * tbe_0 + tr_xxx_yyy[i] - 2.0 * tr_xxx_xxyyy[i] * tke_0 - 2.0 * tr_xxxz_yyyz[i] * tke_0 + 4.0 * tr_xxxz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xyyz[i] = 3.0 * tr_xx_xyyz[i] + 3.0 * tr_xxz_xyy[i] - 6.0 * tr_xxz_xyyzz[i] * tke_0 - 6.0 * tr_xxzz_xyyz[i] * tbe_0 + tr_xxx_yyz[i] - 2.0 * tr_xxx_xxyyz[i] * tke_0 + tr_xxxz_yy[i] - 2.0 * tr_xxxz_yyzz[i] * tke_0 - 2.0 * tr_xxxz_xxyy[i] * tke_0 + 4.0 * tr_xxxz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_yyz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyyz[i] * tbe_0 - 2.0 * tr_xxxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xyzz[i] = 3.0 * tr_xx_xyzz[i] + 6.0 * tr_xxz_xyz[i] - 6.0 * tr_xxz_xyzzz[i] * tke_0 - 6.0 * tr_xxzz_xyzz[i] * tbe_0 + tr_xxx_yzz[i] - 2.0 * tr_xxx_xxyzz[i] * tke_0 + 2.0 * tr_xxxz_yz[i] - 2.0 * tr_xxxz_yzzz[i] * tke_0 - 4.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_yzz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyzz[i] * tbe_0 - 4.0 * tr_xxxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xzzz[i] = 3.0 * tr_xx_xzzz[i] + 9.0 * tr_xxz_xzz[i] - 6.0 * tr_xxz_xzzzz[i] * tke_0 - 6.0 * tr_xxzz_xzzz[i] * tbe_0 + tr_xxx_zzz[i] - 2.0 * tr_xxx_xxzzz[i] * tke_0 + 3.0 * tr_xxxz_zz[i] - 2.0 * tr_xxxz_zzzz[i] * tke_0 - 6.0 * tr_xxxz_xxzz[i] * tke_0 + 4.0 * tr_xxxz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_zzz[i] * tbe_0 + 4.0 * tr_xxxzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xzzz[i] * tbe_0 - 6.0 * tr_xxxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_yyyy[i] = 3.0 * tr_xx_yyyy[i] - 6.0 * tr_xxz_yyyyz[i] * tke_0 - 6.0 * tr_xxzz_yyyy[i] * tbe_0 - 2.0 * tr_xxx_xyyyy[i] * tke_0 + 4.0 * tr_xxxz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyyy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_yyyz[i] = 3.0 * tr_xx_yyyz[i] + 3.0 * tr_xxz_yyy[i] - 6.0 * tr_xxz_yyyzz[i] * tke_0 - 6.0 * tr_xxzz_yyyz[i] * tbe_0 - 2.0 * tr_xxx_xyyyz[i] * tke_0 - 2.0 * tr_xxxz_xyyy[i] * tke_0 + 4.0 * tr_xxxz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyyz[i] * tbe_0 - 2.0 * tr_xxxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_yyzz[i] = 3.0 * tr_xx_yyzz[i] + 6.0 * tr_xxz_yyz[i] - 6.0 * tr_xxz_yyzzz[i] * tke_0 - 6.0 * tr_xxzz_yyzz[i] * tbe_0 - 2.0 * tr_xxx_xyyzz[i] * tke_0 - 4.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyzz[i] * tbe_0 - 4.0 * tr_xxxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_yzzz[i] = 3.0 * tr_xx_yzzz[i] + 9.0 * tr_xxz_yzz[i] - 6.0 * tr_xxz_yzzzz[i] * tke_0 - 6.0 * tr_xxzz_yzzz[i] * tbe_0 - 2.0 * tr_xxx_xyzzz[i] * tke_0 - 6.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yzzz[i] * tbe_0 - 6.0 * tr_xxxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_zzzz[i] = 3.0 * tr_xx_zzzz[i] + 12.0 * tr_xxz_zzz[i] - 6.0 * tr_xxz_zzzzz[i] * tke_0 - 6.0 * tr_xxzz_zzzz[i] * tbe_0 - 2.0 * tr_xxx_xzzzz[i] * tke_0 - 8.0 * tr_xxxz_xzzz[i] * tke_0 + 4.0 * tr_xxxz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zzzz[i] * tbe_0 - 8.0 * tr_xxxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 495-510 components of targeted buffer : GG

    auto tr_0_0_xz_xxyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 495);

    auto tr_0_0_xz_xxyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 496);

    auto tr_0_0_xz_xxyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 497);

    auto tr_0_0_xz_xxyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 498);

    auto tr_0_0_xz_xxyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 499);

    auto tr_0_0_xz_xxyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 500);

    auto tr_0_0_xz_xxyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 501);

    auto tr_0_0_xz_xxyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 502);

    auto tr_0_0_xz_xxyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 503);

    auto tr_0_0_xz_xxyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 504);

    auto tr_0_0_xz_xxyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 505);

    auto tr_0_0_xz_xxyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 506);

    auto tr_0_0_xz_xxyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 507);

    auto tr_0_0_xz_xxyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 508);

    auto tr_0_0_xz_xxyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 509);

    #pragma omp simd aligned(tr_0_0_xz_xxyy_xxxx, tr_0_0_xz_xxyy_xxxy, tr_0_0_xz_xxyy_xxxz, tr_0_0_xz_xxyy_xxyy, tr_0_0_xz_xxyy_xxyz, tr_0_0_xz_xxyy_xxzz, tr_0_0_xz_xxyy_xyyy, tr_0_0_xz_xxyy_xyyz, tr_0_0_xz_xxyy_xyzz, tr_0_0_xz_xxyy_xzzz, tr_0_0_xz_xxyy_yyyy, tr_0_0_xz_xxyy_yyyz, tr_0_0_xz_xxyy_yyzz, tr_0_0_xz_xxyy_yzzz, tr_0_0_xz_xxyy_zzzz, tr_xxxyy_xxx, tr_xxxyy_xxxxz, tr_xxxyy_xxxyz, tr_xxxyy_xxxzz, tr_xxxyy_xxy, tr_xxxyy_xxyyz, tr_xxxyy_xxyzz, tr_xxxyy_xxz, tr_xxxyy_xxzzz, tr_xxxyy_xyy, tr_xxxyy_xyyyz, tr_xxxyy_xyyzz, tr_xxxyy_xyz, tr_xxxyy_xyzzz, tr_xxxyy_xzz, tr_xxxyy_xzzzz, tr_xxxyy_yyy, tr_xxxyy_yyyyz, tr_xxxyy_yyyzz, tr_xxxyy_yyz, tr_xxxyy_yyzzz, tr_xxxyy_yzz, tr_xxxyy_yzzzz, tr_xxxyy_zzz, tr_xxxyy_zzzzz, tr_xxxyyz_xxxx, tr_xxxyyz_xxxy, tr_xxxyyz_xxxz, tr_xxxyyz_xxyy, tr_xxxyyz_xxyz, tr_xxxyyz_xxzz, tr_xxxyyz_xyyy, tr_xxxyyz_xyyz, tr_xxxyyz_xyzz, tr_xxxyyz_xzzz, tr_xxxyyz_yyyy, tr_xxxyyz_yyyz, tr_xxxyyz_yyzz, tr_xxxyyz_yzzz, tr_xxxyyz_zzzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxxxz, tr_xxyy_xxxxyz, tr_xxyy_xxxxzz, tr_xxyy_xxxy, tr_xxyy_xxxyyz, tr_xxyy_xxxyzz, tr_xxyy_xxxz, tr_xxyy_xxxzzz, tr_xxyy_xxyy, tr_xxyy_xxyyyz, tr_xxyy_xxyyzz, tr_xxyy_xxyz, tr_xxyy_xxyzzz, tr_xxyy_xxzz, tr_xxyy_xxzzzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyyyz, tr_xxyy_xyyyzz, tr_xxyy_xyyz, tr_xxyy_xyyzzz, tr_xxyy_xyzz, tr_xxyy_xyzzzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_xzzzzz, tr_xxyy_yy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxx, tr_xxyyz_xxxxy, tr_xxyyz_xxxxz, tr_xxyyz_xxxyy, tr_xxyyz_xxxyz, tr_xxyyz_xxxzz, tr_xxyyz_xxy, tr_xxyyz_xxyyy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xxzzz, tr_xxyyz_xyy, tr_xxyyz_xyyyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_xzzzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxyy_xxxx[i] = -4.0 * tr_xyy_xxxxz[i] * tke_0 - 4.0 * tr_xyyz_xxxx[i] * tbe_0 - 8.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xxxy[i] = -4.0 * tr_xyy_xxxyz[i] * tke_0 - 4.0 * tr_xyyz_xxxy[i] * tbe_0 - 6.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xxxz[i] = 2.0 * tr_xyy_xxx[i] - 4.0 * tr_xyy_xxxzz[i] * tke_0 - 4.0 * tr_xyyz_xxxz[i] * tbe_0 + 3.0 * tr_xxyy_xx[i] - 6.0 * tr_xxyy_xxzz[i] * tke_0 - 2.0 * tr_xxyy_xxxx[i] * tke_0 + 4.0 * tr_xxyy_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xxyyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xxyy[i] = -4.0 * tr_xyy_xxyyz[i] * tke_0 - 4.0 * tr_xyyz_xxyy[i] * tbe_0 - 4.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xxyz[i] = 2.0 * tr_xyy_xxy[i] - 4.0 * tr_xyy_xxyzz[i] * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tbe_0 + 2.0 * tr_xxyy_xy[i] - 4.0 * tr_xxyy_xyzz[i] * tke_0 - 2.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_xxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xxzz[i] = 4.0 * tr_xyy_xxz[i] - 4.0 * tr_xyy_xxzzz[i] * tke_0 - 4.0 * tr_xyyz_xxzz[i] * tbe_0 + 4.0 * tr_xxyy_xz[i] - 4.0 * tr_xxyy_xzzz[i] * tke_0 - 4.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_xxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xyyy[i] = -4.0 * tr_xyy_xyyyz[i] * tke_0 - 4.0 * tr_xyyz_xyyy[i] * tbe_0 - 2.0 * tr_xxyy_yyyz[i] * tke_0 + 4.0 * tr_xxyy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xyyz[i] = 2.0 * tr_xyy_xyy[i] - 4.0 * tr_xyy_xyyzz[i] * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tbe_0 + tr_xxyy_yy[i] - 2.0 * tr_xxyy_yyzz[i] * tke_0 - 2.0 * tr_xxyy_xxyy[i] * tke_0 + 4.0 * tr_xxyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_xyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xyzz[i] = 4.0 * tr_xyy_xyz[i] - 4.0 * tr_xyy_xyzzz[i] * tke_0 - 4.0 * tr_xyyz_xyzz[i] * tbe_0 + 2.0 * tr_xxyy_yz[i] - 2.0 * tr_xxyy_yzzz[i] * tke_0 - 4.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_xyz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xzzz[i] = 6.0 * tr_xyy_xzz[i] - 4.0 * tr_xyy_xzzzz[i] * tke_0 - 4.0 * tr_xyyz_xzzz[i] * tbe_0 + 3.0 * tr_xxyy_zz[i] - 2.0 * tr_xxyy_zzzz[i] * tke_0 - 6.0 * tr_xxyy_xxzz[i] * tke_0 + 4.0 * tr_xxyy_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxyy_xzz[i] * tbe_0 + 4.0 * tr_xxxyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_yyyy[i] = -4.0 * tr_xyy_yyyyz[i] * tke_0 - 4.0 * tr_xyyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_yyyz[i] = 2.0 * tr_xyy_yyy[i] - 4.0 * tr_xyy_yyyzz[i] * tke_0 - 4.0 * tr_xyyz_yyyz[i] * tbe_0 - 2.0 * tr_xxyy_xyyy[i] * tke_0 + 4.0 * tr_xxyy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_yyzz[i] = 4.0 * tr_xyy_yyz[i] - 4.0 * tr_xyy_yyzzz[i] * tke_0 - 4.0 * tr_xyyz_yyzz[i] * tbe_0 - 4.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_yyz[i] * tbe_0 + 4.0 * tr_xxxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_yzzz[i] = 6.0 * tr_xyy_yzz[i] - 4.0 * tr_xyy_yzzzz[i] * tke_0 - 4.0 * tr_xyyz_yzzz[i] * tbe_0 - 6.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxyy_yzz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_zzzz[i] = 8.0 * tr_xyy_zzz[i] - 4.0 * tr_xyy_zzzzz[i] * tke_0 - 4.0 * tr_xyyz_zzzz[i] * tbe_0 - 8.0 * tr_xxyy_xzzz[i] * tke_0 + 4.0 * tr_xxyy_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxxyy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-525 components of targeted buffer : GG

    auto tr_0_0_xz_xxyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 510);

    auto tr_0_0_xz_xxyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 511);

    auto tr_0_0_xz_xxyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 512);

    auto tr_0_0_xz_xxyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 513);

    auto tr_0_0_xz_xxyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 514);

    auto tr_0_0_xz_xxyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 515);

    auto tr_0_0_xz_xxyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 516);

    auto tr_0_0_xz_xxyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 517);

    auto tr_0_0_xz_xxyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 518);

    auto tr_0_0_xz_xxyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 519);

    auto tr_0_0_xz_xxyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 520);

    auto tr_0_0_xz_xxyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 521);

    auto tr_0_0_xz_xxyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 522);

    auto tr_0_0_xz_xxyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 523);

    auto tr_0_0_xz_xxyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 524);

    #pragma omp simd aligned(tr_0_0_xz_xxyz_xxxx, tr_0_0_xz_xxyz_xxxy, tr_0_0_xz_xxyz_xxxz, tr_0_0_xz_xxyz_xxyy, tr_0_0_xz_xxyz_xxyz, tr_0_0_xz_xxyz_xxzz, tr_0_0_xz_xxyz_xyyy, tr_0_0_xz_xxyz_xyyz, tr_0_0_xz_xxyz_xyzz, tr_0_0_xz_xxyz_xzzz, tr_0_0_xz_xxyz_yyyy, tr_0_0_xz_xxyz_yyyz, tr_0_0_xz_xxyz_yyzz, tr_0_0_xz_xxyz_yzzz, tr_0_0_xz_xxyz_zzzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxz, tr_xxxyz_xxxyz, tr_xxxyz_xxxzz, tr_xxxyz_xxy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xxzzz, tr_xxxyz_xyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_xzzzz, tr_xxxyz_yyy, tr_xxxyz_yyyyz, tr_xxxyz_yyyzz, tr_xxxyz_yyz, tr_xxxyz_yyzzz, tr_xxxyz_yzz, tr_xxxyz_yzzzz, tr_xxxyz_zzz, tr_xxxyz_zzzzz, tr_xxxyzz_xxxx, tr_xxxyzz_xxxy, tr_xxxyzz_xxxz, tr_xxxyzz_xxyy, tr_xxxyzz_xxyz, tr_xxxyzz_xxzz, tr_xxxyzz_xyyy, tr_xxxyzz_xyyz, tr_xxxyzz_xyzz, tr_xxxyzz_xzzz, tr_xxxyzz_yyyy, tr_xxxyzz_yyyz, tr_xxxyzz_yyzz, tr_xxxyzz_yzzz, tr_xxxyzz_zzzz, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxxxz, tr_xxyz_xxxxyz, tr_xxyz_xxxxzz, tr_xxyz_xxxy, tr_xxyz_xxxyyz, tr_xxyz_xxxyzz, tr_xxyz_xxxz, tr_xxyz_xxxzzz, tr_xxyz_xxyy, tr_xxyz_xxyyyz, tr_xxyz_xxyyzz, tr_xxyz_xxyz, tr_xxyz_xxyzzz, tr_xxyz_xxzz, tr_xxyz_xxzzzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyyyz, tr_xxyz_xyyyzz, tr_xxyz_xyyz, tr_xxyz_xyyzzz, tr_xxyz_xyzz, tr_xxyz_xyzzzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_xzzzzz, tr_xxyz_yy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxx, tr_xxyzz_xxxxy, tr_xxyzz_xxxxz, tr_xxyzz_xxxyy, tr_xxyzz_xxxyz, tr_xxyzz_xxxzz, tr_xxyzz_xxy, tr_xxyzz_xxyyy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xxzzz, tr_xxyzz_xyy, tr_xxyzz_xyyyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_xzzzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxyz_xxxx[i] = 2.0 * tr_xy_xxxx[i] - 4.0 * tr_xyz_xxxxz[i] * tke_0 - 4.0 * tr_xyzz_xxxx[i] * tbe_0 + 4.0 * tr_xxy_xxx[i] - 2.0 * tr_xxy_xxxxx[i] * tke_0 - 8.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xxxy[i] = 2.0 * tr_xy_xxxy[i] - 4.0 * tr_xyz_xxxyz[i] * tke_0 - 4.0 * tr_xyzz_xxxy[i] * tbe_0 + 3.0 * tr_xxy_xxy[i] - 2.0 * tr_xxy_xxxxy[i] * tke_0 - 6.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyzz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xxxz[i] = 2.0 * tr_xy_xxxz[i] + 2.0 * tr_xyz_xxx[i] - 4.0 * tr_xyz_xxxzz[i] * tke_0 - 4.0 * tr_xyzz_xxxz[i] * tbe_0 + 3.0 * tr_xxy_xxz[i] - 2.0 * tr_xxy_xxxxz[i] * tke_0 + 3.0 * tr_xxyz_xx[i] - 6.0 * tr_xxyz_xxzz[i] * tke_0 - 2.0 * tr_xxyz_xxxx[i] * tke_0 + 4.0 * tr_xxyz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xxyzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxxz[i] * tbe_0 - 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xxyy[i] = 2.0 * tr_xy_xxyy[i] - 4.0 * tr_xyz_xxyyz[i] * tke_0 - 4.0 * tr_xyzz_xxyy[i] * tbe_0 + 2.0 * tr_xxy_xyy[i] - 2.0 * tr_xxy_xxxyy[i] * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxyy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xxyz[i] = 2.0 * tr_xy_xxyz[i] + 2.0 * tr_xyz_xxy[i] - 4.0 * tr_xyz_xxyzz[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tbe_0 + 2.0 * tr_xxy_xyz[i] - 2.0 * tr_xxy_xxxyz[i] * tke_0 + 2.0 * tr_xxyz_xy[i] - 4.0 * tr_xxyz_xyzz[i] * tke_0 - 2.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xyz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxyz[i] * tbe_0 - 2.0 * tr_xxxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xxzz[i] = 2.0 * tr_xy_xxzz[i] + 4.0 * tr_xyz_xxz[i] - 4.0 * tr_xyz_xxzzz[i] * tke_0 - 4.0 * tr_xyzz_xxzz[i] * tbe_0 + 2.0 * tr_xxy_xzz[i] - 2.0 * tr_xxy_xxxzz[i] * tke_0 + 4.0 * tr_xxyz_xz[i] - 4.0 * tr_xxyz_xzzz[i] * tke_0 - 4.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxzz[i] * tbe_0 - 4.0 * tr_xxxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xyyy[i] = 2.0 * tr_xy_xyyy[i] - 4.0 * tr_xyz_xyyyz[i] * tke_0 - 4.0 * tr_xyzz_xyyy[i] * tbe_0 + tr_xxy_yyy[i] - 2.0 * tr_xxy_xxyyy[i] * tke_0 - 2.0 * tr_xxyz_yyyz[i] * tke_0 + 4.0 * tr_xxyz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xyyz[i] = 2.0 * tr_xy_xyyz[i] + 2.0 * tr_xyz_xyy[i] - 4.0 * tr_xyz_xyyzz[i] * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tbe_0 + tr_xxy_yyz[i] - 2.0 * tr_xxy_xxyyz[i] * tke_0 + tr_xxyz_yy[i] - 2.0 * tr_xxyz_yyzz[i] * tke_0 - 2.0 * tr_xxyz_xxyy[i] * tke_0 + 4.0 * tr_xxyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyyz[i] * tbe_0 - 2.0 * tr_xxxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xyzz[i] = 2.0 * tr_xy_xyzz[i] + 4.0 * tr_xyz_xyz[i] - 4.0 * tr_xyz_xyzzz[i] * tke_0 - 4.0 * tr_xyzz_xyzz[i] * tbe_0 + tr_xxy_yzz[i] - 2.0 * tr_xxy_xxyzz[i] * tke_0 + 2.0 * tr_xxyz_yz[i] - 2.0 * tr_xxyz_yzzz[i] * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyzz[i] * tbe_0 - 4.0 * tr_xxxyz_xyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xzzz[i] = 2.0 * tr_xy_xzzz[i] + 6.0 * tr_xyz_xzz[i] - 4.0 * tr_xyz_xzzzz[i] * tke_0 - 4.0 * tr_xyzz_xzzz[i] * tbe_0 + tr_xxy_zzz[i] - 2.0 * tr_xxy_xxzzz[i] * tke_0 + 3.0 * tr_xxyz_zz[i] - 2.0 * tr_xxyz_zzzz[i] * tke_0 - 6.0 * tr_xxyz_xxzz[i] * tke_0 + 4.0 * tr_xxyz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xzzz[i] * tbe_0 - 6.0 * tr_xxxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_yyyy[i] = 2.0 * tr_xy_yyyy[i] - 4.0 * tr_xyz_yyyyz[i] * tke_0 - 4.0 * tr_xyzz_yyyy[i] * tbe_0 - 2.0 * tr_xxy_xyyyy[i] * tke_0 + 4.0 * tr_xxyz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_yyyz[i] = 2.0 * tr_xy_yyyz[i] + 2.0 * tr_xyz_yyy[i] - 4.0 * tr_xyz_yyyzz[i] * tke_0 - 4.0 * tr_xyzz_yyyz[i] * tbe_0 - 2.0 * tr_xxy_xyyyz[i] * tke_0 - 2.0 * tr_xxyz_xyyy[i] * tke_0 + 4.0 * tr_xxyz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyyz[i] * tbe_0 - 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_yyzz[i] = 2.0 * tr_xy_yyzz[i] + 4.0 * tr_xyz_yyz[i] - 4.0 * tr_xyz_yyzzz[i] * tke_0 - 4.0 * tr_xyzz_yyzz[i] * tbe_0 - 2.0 * tr_xxy_xyyzz[i] * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyzz[i] * tbe_0 - 4.0 * tr_xxxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_yzzz[i] = 2.0 * tr_xy_yzzz[i] + 6.0 * tr_xyz_yzz[i] - 4.0 * tr_xyz_yzzzz[i] * tke_0 - 4.0 * tr_xyzz_yzzz[i] * tbe_0 - 2.0 * tr_xxy_xyzzz[i] * tke_0 - 6.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yzzz[i] * tbe_0 - 6.0 * tr_xxxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_zzzz[i] = 2.0 * tr_xy_zzzz[i] + 8.0 * tr_xyz_zzz[i] - 4.0 * tr_xyz_zzzzz[i] * tke_0 - 4.0 * tr_xyzz_zzzz[i] * tbe_0 - 2.0 * tr_xxy_xzzzz[i] * tke_0 - 8.0 * tr_xxyz_xzzz[i] * tke_0 + 4.0 * tr_xxyz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_zzzz[i] * tbe_0 - 8.0 * tr_xxxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 525-540 components of targeted buffer : GG

    auto tr_0_0_xz_xxzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 525);

    auto tr_0_0_xz_xxzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 526);

    auto tr_0_0_xz_xxzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 527);

    auto tr_0_0_xz_xxzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 528);

    auto tr_0_0_xz_xxzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 529);

    auto tr_0_0_xz_xxzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 530);

    auto tr_0_0_xz_xxzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 531);

    auto tr_0_0_xz_xxzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 532);

    auto tr_0_0_xz_xxzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 533);

    auto tr_0_0_xz_xxzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 534);

    auto tr_0_0_xz_xxzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 535);

    auto tr_0_0_xz_xxzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 536);

    auto tr_0_0_xz_xxzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 537);

    auto tr_0_0_xz_xxzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 538);

    auto tr_0_0_xz_xxzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 539);

    #pragma omp simd aligned(tr_0_0_xz_xxzz_xxxx, tr_0_0_xz_xxzz_xxxy, tr_0_0_xz_xxzz_xxxz, tr_0_0_xz_xxzz_xxyy, tr_0_0_xz_xxzz_xxyz, tr_0_0_xz_xxzz_xxzz, tr_0_0_xz_xxzz_xyyy, tr_0_0_xz_xxzz_xyyz, tr_0_0_xz_xxzz_xyzz, tr_0_0_xz_xxzz_xzzz, tr_0_0_xz_xxzz_yyyy, tr_0_0_xz_xxzz_yyyz, tr_0_0_xz_xxzz_yyzz, tr_0_0_xz_xxzz_yzzz, tr_0_0_xz_xxzz_zzzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxxxz, tr_xxxzz_xxxyz, tr_xxxzz_xxxzz, tr_xxxzz_xxy, tr_xxxzz_xxyyz, tr_xxxzz_xxyzz, tr_xxxzz_xxz, tr_xxxzz_xxzzz, tr_xxxzz_xyy, tr_xxxzz_xyyyz, tr_xxxzz_xyyzz, tr_xxxzz_xyz, tr_xxxzz_xyzzz, tr_xxxzz_xzz, tr_xxxzz_xzzzz, tr_xxxzz_yyy, tr_xxxzz_yyyyz, tr_xxxzz_yyyzz, tr_xxxzz_yyz, tr_xxxzz_yyzzz, tr_xxxzz_yzz, tr_xxxzz_yzzzz, tr_xxxzz_zzz, tr_xxxzz_zzzzz, tr_xxxzzz_xxxx, tr_xxxzzz_xxxy, tr_xxxzzz_xxxz, tr_xxxzzz_xxyy, tr_xxxzzz_xxyz, tr_xxxzzz_xxzz, tr_xxxzzz_xyyy, tr_xxxzzz_xyyz, tr_xxxzzz_xyzz, tr_xxxzzz_xzzz, tr_xxxzzz_yyyy, tr_xxxzzz_yyyz, tr_xxxzzz_yyzz, tr_xxxzzz_yzzz, tr_xxxzzz_zzzz, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxxxz, tr_xxzz_xxxxyz, tr_xxzz_xxxxzz, tr_xxzz_xxxy, tr_xxzz_xxxyyz, tr_xxzz_xxxyzz, tr_xxzz_xxxz, tr_xxzz_xxxzzz, tr_xxzz_xxyy, tr_xxzz_xxyyyz, tr_xxzz_xxyyzz, tr_xxzz_xxyz, tr_xxzz_xxyzzz, tr_xxzz_xxzz, tr_xxzz_xxzzzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyyyz, tr_xxzz_xyyyzz, tr_xxzz_xyyz, tr_xxzz_xyyzzz, tr_xxzz_xyzz, tr_xxzz_xyzzzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_xzzzzz, tr_xxzz_yy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxxxx, tr_xxzzz_xxxxy, tr_xxzzz_xxxxz, tr_xxzzz_xxxyy, tr_xxzzz_xxxyz, tr_xxzzz_xxxzz, tr_xxzzz_xxy, tr_xxzzz_xxyyy, tr_xxzzz_xxyyz, tr_xxzzz_xxyzz, tr_xxzzz_xxz, tr_xxzzz_xxzzz, tr_xxzzz_xyy, tr_xxzzz_xyyyy, tr_xxzzz_xyyyz, tr_xxzzz_xyyzz, tr_xxzzz_xyz, tr_xxzzz_xyzzz, tr_xxzzz_xzz, tr_xxzzz_xzzzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxzz_xxxx[i] = 4.0 * tr_xz_xxxx[i] - 4.0 * tr_xzz_xxxxz[i] * tke_0 - 4.0 * tr_xzzz_xxxx[i] * tbe_0 + 8.0 * tr_xxz_xxx[i] - 4.0 * tr_xxz_xxxxx[i] * tke_0 - 8.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xxxy[i] = 4.0 * tr_xz_xxxy[i] - 4.0 * tr_xzz_xxxyz[i] * tke_0 - 4.0 * tr_xzzz_xxxy[i] * tbe_0 + 6.0 * tr_xxz_xxy[i] - 4.0 * tr_xxz_xxxxy[i] * tke_0 - 6.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xxzzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xxxz[i] = 4.0 * tr_xz_xxxz[i] + 2.0 * tr_xzz_xxx[i] - 4.0 * tr_xzz_xxxzz[i] * tke_0 - 4.0 * tr_xzzz_xxxz[i] * tbe_0 + 6.0 * tr_xxz_xxz[i] - 4.0 * tr_xxz_xxxxz[i] * tke_0 + 3.0 * tr_xxzz_xx[i] - 6.0 * tr_xxzz_xxzz[i] * tke_0 - 2.0 * tr_xxzz_xxxx[i] * tke_0 + 4.0 * tr_xxzz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xxzzz_xxz[i] * tbe_0 + 4.0 * tr_xxzzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxxz[i] * tbe_0 - 2.0 * tr_xxxzz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xxyy[i] = 4.0 * tr_xz_xxyy[i] - 4.0 * tr_xzz_xxyyz[i] * tke_0 - 4.0 * tr_xzzz_xxyy[i] * tbe_0 + 4.0 * tr_xxz_xyy[i] - 4.0 * tr_xxz_xxxyy[i] * tke_0 - 4.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxyy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xxyz[i] = 4.0 * tr_xz_xxyz[i] + 2.0 * tr_xzz_xxy[i] - 4.0 * tr_xzz_xxyzz[i] * tke_0 - 4.0 * tr_xzzz_xxyz[i] * tbe_0 + 4.0 * tr_xxz_xyz[i] - 4.0 * tr_xxz_xxxyz[i] * tke_0 + 2.0 * tr_xxzz_xy[i] - 4.0 * tr_xxzz_xyzz[i] * tke_0 - 2.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xyz[i] * tbe_0 + 4.0 * tr_xxzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxyz[i] * tbe_0 - 2.0 * tr_xxxzz_xxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xxzz[i] = 4.0 * tr_xz_xxzz[i] + 4.0 * tr_xzz_xxz[i] - 4.0 * tr_xzz_xxzzz[i] * tke_0 - 4.0 * tr_xzzz_xxzz[i] * tbe_0 + 4.0 * tr_xxz_xzz[i] - 4.0 * tr_xxz_xxxzz[i] * tke_0 + 4.0 * tr_xxzz_xz[i] - 4.0 * tr_xxzz_xzzz[i] * tke_0 - 4.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xzz[i] * tbe_0 + 4.0 * tr_xxzzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxzz[i] * tbe_0 - 4.0 * tr_xxxzz_xxz[i] * tbe_0 + 4.0 * tr_xxxzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xyyy[i] = 4.0 * tr_xz_xyyy[i] - 4.0 * tr_xzz_xyyyz[i] * tke_0 - 4.0 * tr_xzzz_xyyy[i] * tbe_0 + 2.0 * tr_xxz_yyy[i] - 4.0 * tr_xxz_xxyyy[i] * tke_0 - 2.0 * tr_xxzz_yyyz[i] * tke_0 + 4.0 * tr_xxzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xyyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xyyz[i] = 4.0 * tr_xz_xyyz[i] + 2.0 * tr_xzz_xyy[i] - 4.0 * tr_xzz_xyyzz[i] * tke_0 - 4.0 * tr_xzzz_xyyz[i] * tbe_0 + 2.0 * tr_xxz_yyz[i] - 4.0 * tr_xxz_xxyyz[i] * tke_0 + tr_xxzz_yy[i] - 2.0 * tr_xxzz_yyzz[i] * tke_0 - 2.0 * tr_xxzz_xxyy[i] * tke_0 + 4.0 * tr_xxzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_yyz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xyyz[i] * tbe_0 - 2.0 * tr_xxxzz_xyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xyzz[i] = 4.0 * tr_xz_xyzz[i] + 4.0 * tr_xzz_xyz[i] - 4.0 * tr_xzz_xyzzz[i] * tke_0 - 4.0 * tr_xzzz_xyzz[i] * tbe_0 + 2.0 * tr_xxz_yzz[i] - 4.0 * tr_xxz_xxyzz[i] * tke_0 + 2.0 * tr_xxzz_yz[i] - 2.0 * tr_xxzz_yzzz[i] * tke_0 - 4.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_yzz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xyzz[i] * tbe_0 - 4.0 * tr_xxxzz_xyz[i] * tbe_0 + 4.0 * tr_xxxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xzzz[i] = 4.0 * tr_xz_xzzz[i] + 6.0 * tr_xzz_xzz[i] - 4.0 * tr_xzz_xzzzz[i] * tke_0 - 4.0 * tr_xzzz_xzzz[i] * tbe_0 + 2.0 * tr_xxz_zzz[i] - 4.0 * tr_xxz_xxzzz[i] * tke_0 + 3.0 * tr_xxzz_zz[i] - 2.0 * tr_xxzz_zzzz[i] * tke_0 - 6.0 * tr_xxzz_xxzz[i] * tke_0 + 4.0 * tr_xxzz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xzzz[i] * tbe_0 - 6.0 * tr_xxxzz_xzz[i] * tbe_0 + 4.0 * tr_xxxzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_yyyy[i] = 4.0 * tr_xz_yyyy[i] - 4.0 * tr_xzz_yyyyz[i] * tke_0 - 4.0 * tr_xzzz_yyyy[i] * tbe_0 - 4.0 * tr_xxz_xyyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yyyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_yyyz[i] = 4.0 * tr_xz_yyyz[i] + 2.0 * tr_xzz_yyy[i] - 4.0 * tr_xzz_yyyzz[i] * tke_0 - 4.0 * tr_xzzz_yyyz[i] * tbe_0 - 4.0 * tr_xxz_xyyyz[i] * tke_0 - 2.0 * tr_xxzz_xyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yyyz[i] * tbe_0 - 2.0 * tr_xxxzz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_yyzz[i] = 4.0 * tr_xz_yyzz[i] + 4.0 * tr_xzz_yyz[i] - 4.0 * tr_xzz_yyzzz[i] * tke_0 - 4.0 * tr_xzzz_yyzz[i] * tbe_0 - 4.0 * tr_xxz_xyyzz[i] * tke_0 - 4.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yyzz[i] * tbe_0 - 4.0 * tr_xxxzz_yyz[i] * tbe_0 + 4.0 * tr_xxxzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_yzzz[i] = 4.0 * tr_xz_yzzz[i] + 6.0 * tr_xzz_yzz[i] - 4.0 * tr_xzz_yzzzz[i] * tke_0 - 4.0 * tr_xzzz_yzzz[i] * tbe_0 - 4.0 * tr_xxz_xyzzz[i] * tke_0 - 6.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yzzz[i] * tbe_0 - 6.0 * tr_xxxzz_yzz[i] * tbe_0 + 4.0 * tr_xxxzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_zzzz[i] = 4.0 * tr_xz_zzzz[i] + 8.0 * tr_xzz_zzz[i] - 4.0 * tr_xzz_zzzzz[i] * tke_0 - 4.0 * tr_xzzz_zzzz[i] * tbe_0 - 4.0 * tr_xxz_xzzzz[i] * tke_0 - 8.0 * tr_xxzz_xzzz[i] * tke_0 + 4.0 * tr_xxzz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_zzzz[i] * tbe_0 - 8.0 * tr_xxxzz_zzz[i] * tbe_0 + 4.0 * tr_xxxzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 540-555 components of targeted buffer : GG

    auto tr_0_0_xz_xyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 540);

    auto tr_0_0_xz_xyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 541);

    auto tr_0_0_xz_xyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 542);

    auto tr_0_0_xz_xyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 543);

    auto tr_0_0_xz_xyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 544);

    auto tr_0_0_xz_xyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 545);

    auto tr_0_0_xz_xyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 546);

    auto tr_0_0_xz_xyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 547);

    auto tr_0_0_xz_xyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 548);

    auto tr_0_0_xz_xyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 549);

    auto tr_0_0_xz_xyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 550);

    auto tr_0_0_xz_xyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 551);

    auto tr_0_0_xz_xyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 552);

    auto tr_0_0_xz_xyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 553);

    auto tr_0_0_xz_xyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 554);

    #pragma omp simd aligned(tr_0_0_xz_xyyy_xxxx, tr_0_0_xz_xyyy_xxxy, tr_0_0_xz_xyyy_xxxz, tr_0_0_xz_xyyy_xxyy, tr_0_0_xz_xyyy_xxyz, tr_0_0_xz_xyyy_xxzz, tr_0_0_xz_xyyy_xyyy, tr_0_0_xz_xyyy_xyyz, tr_0_0_xz_xyyy_xyzz, tr_0_0_xz_xyyy_xzzz, tr_0_0_xz_xyyy_yyyy, tr_0_0_xz_xyyy_yyyz, tr_0_0_xz_xyyy_yyzz, tr_0_0_xz_xyyy_yzzz, tr_0_0_xz_xyyy_zzzz, tr_xxyyy_xxx, tr_xxyyy_xxxxz, tr_xxyyy_xxxyz, tr_xxyyy_xxxzz, tr_xxyyy_xxy, tr_xxyyy_xxyyz, tr_xxyyy_xxyzz, tr_xxyyy_xxz, tr_xxyyy_xxzzz, tr_xxyyy_xyy, tr_xxyyy_xyyyz, tr_xxyyy_xyyzz, tr_xxyyy_xyz, tr_xxyyy_xyzzz, tr_xxyyy_xzz, tr_xxyyy_xzzzz, tr_xxyyy_yyy, tr_xxyyy_yyyyz, tr_xxyyy_yyyzz, tr_xxyyy_yyz, tr_xxyyy_yyzzz, tr_xxyyy_yzz, tr_xxyyy_yzzzz, tr_xxyyy_zzz, tr_xxyyy_zzzzz, tr_xxyyyz_xxxx, tr_xxyyyz_xxxy, tr_xxyyyz_xxxz, tr_xxyyyz_xxyy, tr_xxyyyz_xxyz, tr_xxyyyz_xxzz, tr_xxyyyz_xyyy, tr_xxyyyz_xyyz, tr_xxyyyz_xyzz, tr_xxyyyz_xzzz, tr_xxyyyz_yyyy, tr_xxyyyz_yyyz, tr_xxyyyz_yyzz, tr_xxyyyz_yzzz, tr_xxyyyz_zzzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxxxz, tr_xyyy_xxxxyz, tr_xyyy_xxxxzz, tr_xyyy_xxxy, tr_xyyy_xxxyyz, tr_xyyy_xxxyzz, tr_xyyy_xxxz, tr_xyyy_xxxzzz, tr_xyyy_xxyy, tr_xyyy_xxyyyz, tr_xyyy_xxyyzz, tr_xyyy_xxyz, tr_xyyy_xxyzzz, tr_xyyy_xxzz, tr_xyyy_xxzzzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyyyz, tr_xyyy_xyyyzz, tr_xyyy_xyyz, tr_xyyy_xyyzzz, tr_xyyy_xyzz, tr_xyyy_xyzzzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_xzzzzz, tr_xyyy_yy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxx, tr_xyyyz_xxxxy, tr_xyyyz_xxxxz, tr_xyyyz_xxxyy, tr_xyyyz_xxxyz, tr_xyyyz_xxxzz, tr_xyyyz_xxy, tr_xyyyz_xxyyy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xxzzz, tr_xyyyz_xyy, tr_xyyyz_xyyyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_xzzzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_yyy_xxx, tr_yyy_xxxxz, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyy_zzzzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyyy_xxxx[i] = -2.0 * tr_yyy_xxxxz[i] * tke_0 - 2.0 * tr_yyyz_xxxx[i] * tbe_0 - 8.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xxxy[i] = -2.0 * tr_yyy_xxxyz[i] * tke_0 - 2.0 * tr_yyyz_xxxy[i] * tbe_0 - 6.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xxxz[i] = tr_yyy_xxx[i] - 2.0 * tr_yyy_xxxzz[i] * tke_0 - 2.0 * tr_yyyz_xxxz[i] * tbe_0 + 3.0 * tr_xyyy_xx[i] - 6.0 * tr_xyyy_xxzz[i] * tke_0 - 2.0 * tr_xyyy_xxxx[i] * tke_0 + 4.0 * tr_xyyy_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xyyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xxyy[i] = -2.0 * tr_yyy_xxyyz[i] * tke_0 - 2.0 * tr_yyyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xxyz[i] = tr_yyy_xxy[i] - 2.0 * tr_yyy_xxyzz[i] * tke_0 - 2.0 * tr_yyyz_xxyz[i] * tbe_0 + 2.0 * tr_xyyy_xy[i] - 4.0 * tr_xyyy_xyzz[i] * tke_0 - 2.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xxzz[i] = 2.0 * tr_yyy_xxz[i] - 2.0 * tr_yyy_xxzzz[i] * tke_0 - 2.0 * tr_yyyz_xxzz[i] * tbe_0 + 4.0 * tr_xyyy_xz[i] - 4.0 * tr_xyyy_xzzz[i] * tke_0 - 4.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xyyy[i] = -2.0 * tr_yyy_xyyyz[i] * tke_0 - 2.0 * tr_yyyz_xyyy[i] * tbe_0 - 2.0 * tr_xyyy_yyyz[i] * tke_0 + 4.0 * tr_xyyy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xyyz[i] = tr_yyy_xyy[i] - 2.0 * tr_yyy_xyyzz[i] * tke_0 - 2.0 * tr_yyyz_xyyz[i] * tbe_0 + tr_xyyy_yy[i] - 2.0 * tr_xyyy_yyzz[i] * tke_0 - 2.0 * tr_xyyy_xxyy[i] * tke_0 + 4.0 * tr_xyyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xyzz[i] = 2.0 * tr_yyy_xyz[i] - 2.0 * tr_yyy_xyzzz[i] * tke_0 - 2.0 * tr_yyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyyy_yz[i] - 2.0 * tr_xyyy_yzzz[i] * tke_0 - 4.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xzzz[i] = 3.0 * tr_yyy_xzz[i] - 2.0 * tr_yyy_xzzzz[i] * tke_0 - 2.0 * tr_yyyz_xzzz[i] * tbe_0 + 3.0 * tr_xyyy_zz[i] - 2.0 * tr_xyyy_zzzz[i] * tke_0 - 6.0 * tr_xyyy_xxzz[i] * tke_0 + 4.0 * tr_xyyy_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_yyyy[i] = -2.0 * tr_yyy_yyyyz[i] * tke_0 - 2.0 * tr_yyyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_yyyz[i] = tr_yyy_yyy[i] - 2.0 * tr_yyy_yyyzz[i] * tke_0 - 2.0 * tr_yyyz_yyyz[i] * tbe_0 - 2.0 * tr_xyyy_xyyy[i] * tke_0 + 4.0 * tr_xyyy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_yyzz[i] = 2.0 * tr_yyy_yyz[i] - 2.0 * tr_yyy_yyzzz[i] * tke_0 - 2.0 * tr_yyyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_yzzz[i] = 3.0 * tr_yyy_yzz[i] - 2.0 * tr_yyy_yzzzz[i] * tke_0 - 2.0 * tr_yyyz_yzzz[i] * tbe_0 - 6.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_zzzz[i] = 4.0 * tr_yyy_zzz[i] - 2.0 * tr_yyy_zzzzz[i] * tke_0 - 2.0 * tr_yyyz_zzzz[i] * tbe_0 - 8.0 * tr_xyyy_xzzz[i] * tke_0 + 4.0 * tr_xyyy_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxyyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 555-570 components of targeted buffer : GG

    auto tr_0_0_xz_xyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 555);

    auto tr_0_0_xz_xyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 556);

    auto tr_0_0_xz_xyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 557);

    auto tr_0_0_xz_xyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 558);

    auto tr_0_0_xz_xyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 559);

    auto tr_0_0_xz_xyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 560);

    auto tr_0_0_xz_xyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 561);

    auto tr_0_0_xz_xyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 562);

    auto tr_0_0_xz_xyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 563);

    auto tr_0_0_xz_xyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 564);

    auto tr_0_0_xz_xyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 565);

    auto tr_0_0_xz_xyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 566);

    auto tr_0_0_xz_xyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 567);

    auto tr_0_0_xz_xyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 568);

    auto tr_0_0_xz_xyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 569);

    #pragma omp simd aligned(tr_0_0_xz_xyyz_xxxx, tr_0_0_xz_xyyz_xxxy, tr_0_0_xz_xyyz_xxxz, tr_0_0_xz_xyyz_xxyy, tr_0_0_xz_xyyz_xxyz, tr_0_0_xz_xyyz_xxzz, tr_0_0_xz_xyyz_xyyy, tr_0_0_xz_xyyz_xyyz, tr_0_0_xz_xyyz_xyzz, tr_0_0_xz_xyyz_xzzz, tr_0_0_xz_xyyz_yyyy, tr_0_0_xz_xyyz_yyyz, tr_0_0_xz_xyyz_yyzz, tr_0_0_xz_xyyz_yzzz, tr_0_0_xz_xyyz_zzzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxz, tr_xxyyz_xxxyz, tr_xxyyz_xxxzz, tr_xxyyz_xxy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xxzzz, tr_xxyyz_xyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_xzzzz, tr_xxyyz_yyy, tr_xxyyz_yyyyz, tr_xxyyz_yyyzz, tr_xxyyz_yyz, tr_xxyyz_yyzzz, tr_xxyyz_yzz, tr_xxyyz_yzzzz, tr_xxyyz_zzz, tr_xxyyz_zzzzz, tr_xxyyzz_xxxx, tr_xxyyzz_xxxy, tr_xxyyzz_xxxz, tr_xxyyzz_xxyy, tr_xxyyzz_xxyz, tr_xxyyzz_xxzz, tr_xxyyzz_xyyy, tr_xxyyzz_xyyz, tr_xxyyzz_xyzz, tr_xxyyzz_xzzz, tr_xxyyzz_yyyy, tr_xxyyzz_yyyz, tr_xxyyzz_yyzz, tr_xxyyzz_yzzz, tr_xxyyzz_zzzz, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxxxz, tr_xyyz_xxxxyz, tr_xyyz_xxxxzz, tr_xyyz_xxxy, tr_xyyz_xxxyyz, tr_xyyz_xxxyzz, tr_xyyz_xxxz, tr_xyyz_xxxzzz, tr_xyyz_xxyy, tr_xyyz_xxyyyz, tr_xyyz_xxyyzz, tr_xyyz_xxyz, tr_xyyz_xxyzzz, tr_xyyz_xxzz, tr_xyyz_xxzzzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyyyz, tr_xyyz_xyyyzz, tr_xyyz_xyyz, tr_xyyz_xyyzzz, tr_xyyz_xyzz, tr_xyyz_xyzzzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_xzzzzz, tr_xyyz_yy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxx, tr_xyyzz_xxxxy, tr_xyyzz_xxxxz, tr_xyyzz_xxxyy, tr_xyyzz_xxxyz, tr_xyyzz_xxxzz, tr_xyyzz_xxy, tr_xyyzz_xxyyy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xxzzz, tr_xyyzz_xyy, tr_xyyzz_xyyyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_xzzzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyyz_xxxx[i] = tr_yy_xxxx[i] - 2.0 * tr_yyz_xxxxz[i] * tke_0 - 2.0 * tr_yyzz_xxxx[i] * tbe_0 + 4.0 * tr_xyy_xxx[i] - 2.0 * tr_xyy_xxxxx[i] * tke_0 - 8.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xxxy[i] = tr_yy_xxxy[i] - 2.0 * tr_yyz_xxxyz[i] * tke_0 - 2.0 * tr_yyzz_xxxy[i] * tbe_0 + 3.0 * tr_xyy_xxy[i] - 2.0 * tr_xyy_xxxxy[i] * tke_0 - 6.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xxxz[i] = tr_yy_xxxz[i] + tr_yyz_xxx[i] - 2.0 * tr_yyz_xxxzz[i] * tke_0 - 2.0 * tr_yyzz_xxxz[i] * tbe_0 + 3.0 * tr_xyy_xxz[i] - 2.0 * tr_xyy_xxxxz[i] * tke_0 + 3.0 * tr_xyyz_xx[i] - 6.0 * tr_xyyz_xxzz[i] * tke_0 - 2.0 * tr_xyyz_xxxx[i] * tke_0 + 4.0 * tr_xyyz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xyyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxxz[i] * tbe_0 - 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xxyy[i] = tr_yy_xxyy[i] - 2.0 * tr_yyz_xxyyz[i] * tke_0 - 2.0 * tr_yyzz_xxyy[i] * tbe_0 + 2.0 * tr_xyy_xyy[i] - 2.0 * tr_xyy_xxxyy[i] * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xxyz[i] = tr_yy_xxyz[i] + tr_yyz_xxy[i] - 2.0 * tr_yyz_xxyzz[i] * tke_0 - 2.0 * tr_yyzz_xxyz[i] * tbe_0 + 2.0 * tr_xyy_xyz[i] - 2.0 * tr_xyy_xxxyz[i] * tke_0 + 2.0 * tr_xyyz_xy[i] - 4.0 * tr_xyyz_xyzz[i] * tke_0 - 2.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxyz[i] * tbe_0 - 2.0 * tr_xxyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xxzz[i] = tr_yy_xxzz[i] + 2.0 * tr_yyz_xxz[i] - 2.0 * tr_yyz_xxzzz[i] * tke_0 - 2.0 * tr_yyzz_xxzz[i] * tbe_0 + 2.0 * tr_xyy_xzz[i] - 2.0 * tr_xyy_xxxzz[i] * tke_0 + 4.0 * tr_xyyz_xz[i] - 4.0 * tr_xyyz_xzzz[i] * tke_0 - 4.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxzz[i] * tbe_0 - 4.0 * tr_xxyyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xyyy[i] = tr_yy_xyyy[i] - 2.0 * tr_yyz_xyyyz[i] * tke_0 - 2.0 * tr_yyzz_xyyy[i] * tbe_0 + tr_xyy_yyy[i] - 2.0 * tr_xyy_xxyyy[i] * tke_0 - 2.0 * tr_xyyz_yyyz[i] * tke_0 + 4.0 * tr_xyyz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xyyz[i] = tr_yy_xyyz[i] + tr_yyz_xyy[i] - 2.0 * tr_yyz_xyyzz[i] * tke_0 - 2.0 * tr_yyzz_xyyz[i] * tbe_0 + tr_xyy_yyz[i] - 2.0 * tr_xyy_xxyyz[i] * tke_0 + tr_xyyz_yy[i] - 2.0 * tr_xyyz_yyzz[i] * tke_0 - 2.0 * tr_xyyz_xxyy[i] * tke_0 + 4.0 * tr_xyyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyyz[i] * tbe_0 - 2.0 * tr_xxyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xyzz[i] = tr_yy_xyzz[i] + 2.0 * tr_yyz_xyz[i] - 2.0 * tr_yyz_xyzzz[i] * tke_0 - 2.0 * tr_yyzz_xyzz[i] * tbe_0 + tr_xyy_yzz[i] - 2.0 * tr_xyy_xxyzz[i] * tke_0 + 2.0 * tr_xyyz_yz[i] - 2.0 * tr_xyyz_yzzz[i] * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyzz[i] * tbe_0 - 4.0 * tr_xxyyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xzzz[i] = tr_yy_xzzz[i] + 3.0 * tr_yyz_xzz[i] - 2.0 * tr_yyz_xzzzz[i] * tke_0 - 2.0 * tr_yyzz_xzzz[i] * tbe_0 + tr_xyy_zzz[i] - 2.0 * tr_xyy_xxzzz[i] * tke_0 + 3.0 * tr_xyyz_zz[i] - 2.0 * tr_xyyz_zzzz[i] * tke_0 - 6.0 * tr_xyyz_xxzz[i] * tke_0 + 4.0 * tr_xyyz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xzzz[i] * tbe_0 - 6.0 * tr_xxyyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_yyyy[i] = tr_yy_yyyy[i] - 2.0 * tr_yyz_yyyyz[i] * tke_0 - 2.0 * tr_yyzz_yyyy[i] * tbe_0 - 2.0 * tr_xyy_xyyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_yyyz[i] = tr_yy_yyyz[i] + tr_yyz_yyy[i] - 2.0 * tr_yyz_yyyzz[i] * tke_0 - 2.0 * tr_yyzz_yyyz[i] * tbe_0 - 2.0 * tr_xyy_xyyyz[i] * tke_0 - 2.0 * tr_xyyz_xyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyyz[i] * tbe_0 - 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_yyzz[i] = tr_yy_yyzz[i] + 2.0 * tr_yyz_yyz[i] - 2.0 * tr_yyz_yyzzz[i] * tke_0 - 2.0 * tr_yyzz_yyzz[i] * tbe_0 - 2.0 * tr_xyy_xyyzz[i] * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyzz[i] * tbe_0 - 4.0 * tr_xxyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_yzzz[i] = tr_yy_yzzz[i] + 3.0 * tr_yyz_yzz[i] - 2.0 * tr_yyz_yzzzz[i] * tke_0 - 2.0 * tr_yyzz_yzzz[i] * tbe_0 - 2.0 * tr_xyy_xyzzz[i] * tke_0 - 6.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yzzz[i] * tbe_0 - 6.0 * tr_xxyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_zzzz[i] = tr_yy_zzzz[i] + 4.0 * tr_yyz_zzz[i] - 2.0 * tr_yyz_zzzzz[i] * tke_0 - 2.0 * tr_yyzz_zzzz[i] * tbe_0 - 2.0 * tr_xyy_xzzzz[i] * tke_0 - 8.0 * tr_xyyz_xzzz[i] * tke_0 + 4.0 * tr_xyyz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_zzzz[i] * tbe_0 - 8.0 * tr_xxyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 570-585 components of targeted buffer : GG

    auto tr_0_0_xz_xyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 570);

    auto tr_0_0_xz_xyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 571);

    auto tr_0_0_xz_xyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 572);

    auto tr_0_0_xz_xyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 573);

    auto tr_0_0_xz_xyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 574);

    auto tr_0_0_xz_xyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 575);

    auto tr_0_0_xz_xyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 576);

    auto tr_0_0_xz_xyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 577);

    auto tr_0_0_xz_xyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 578);

    auto tr_0_0_xz_xyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 579);

    auto tr_0_0_xz_xyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 580);

    auto tr_0_0_xz_xyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 581);

    auto tr_0_0_xz_xyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 582);

    auto tr_0_0_xz_xyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 583);

    auto tr_0_0_xz_xyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 584);

    #pragma omp simd aligned(tr_0_0_xz_xyzz_xxxx, tr_0_0_xz_xyzz_xxxy, tr_0_0_xz_xyzz_xxxz, tr_0_0_xz_xyzz_xxyy, tr_0_0_xz_xyzz_xxyz, tr_0_0_xz_xyzz_xxzz, tr_0_0_xz_xyzz_xyyy, tr_0_0_xz_xyzz_xyyz, tr_0_0_xz_xyzz_xyzz, tr_0_0_xz_xyzz_xzzz, tr_0_0_xz_xyzz_yyyy, tr_0_0_xz_xyzz_yyyz, tr_0_0_xz_xyzz_yyzz, tr_0_0_xz_xyzz_yzzz, tr_0_0_xz_xyzz_zzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxz, tr_xxyzz_xxxyz, tr_xxyzz_xxxzz, tr_xxyzz_xxy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xxzzz, tr_xxyzz_xyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_xzzzz, tr_xxyzz_yyy, tr_xxyzz_yyyyz, tr_xxyzz_yyyzz, tr_xxyzz_yyz, tr_xxyzz_yyzzz, tr_xxyzz_yzz, tr_xxyzz_yzzzz, tr_xxyzz_zzz, tr_xxyzz_zzzzz, tr_xxyzzz_xxxx, tr_xxyzzz_xxxy, tr_xxyzzz_xxxz, tr_xxyzzz_xxyy, tr_xxyzzz_xxyz, tr_xxyzzz_xxzz, tr_xxyzzz_xyyy, tr_xxyzzz_xyyz, tr_xxyzzz_xyzz, tr_xxyzzz_xzzz, tr_xxyzzz_yyyy, tr_xxyzzz_yyyz, tr_xxyzzz_yyzz, tr_xxyzzz_yzzz, tr_xxyzzz_zzzz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxxxz, tr_xyzz_xxxxyz, tr_xyzz_xxxxzz, tr_xyzz_xxxy, tr_xyzz_xxxyyz, tr_xyzz_xxxyzz, tr_xyzz_xxxz, tr_xyzz_xxxzzz, tr_xyzz_xxyy, tr_xyzz_xxyyyz, tr_xyzz_xxyyzz, tr_xyzz_xxyz, tr_xyzz_xxyzzz, tr_xyzz_xxzz, tr_xyzz_xxzzzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyyyz, tr_xyzz_xyyyzz, tr_xyzz_xyyz, tr_xyzz_xyyzzz, tr_xyzz_xyzz, tr_xyzz_xyzzzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_xzzzzz, tr_xyzz_yy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxx, tr_xyzzz_xxxxy, tr_xyzzz_xxxxz, tr_xyzzz_xxxyy, tr_xyzzz_xxxyz, tr_xyzzz_xxxzz, tr_xyzzz_xxy, tr_xyzzz_xxyyy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xxzzz, tr_xyzzz_xyy, tr_xyzzz_xyyyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_xzzzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyzz_xxxx[i] = 2.0 * tr_yz_xxxx[i] - 2.0 * tr_yzz_xxxxz[i] * tke_0 - 2.0 * tr_yzzz_xxxx[i] * tbe_0 + 8.0 * tr_xyz_xxx[i] - 4.0 * tr_xyz_xxxxx[i] * tke_0 - 8.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xxxy[i] = 2.0 * tr_yz_xxxy[i] - 2.0 * tr_yzz_xxxyz[i] * tke_0 - 2.0 * tr_yzzz_xxxy[i] * tbe_0 + 6.0 * tr_xyz_xxy[i] - 4.0 * tr_xyz_xxxxy[i] * tke_0 - 6.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xyzzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xxxz[i] = 2.0 * tr_yz_xxxz[i] + tr_yzz_xxx[i] - 2.0 * tr_yzz_xxxzz[i] * tke_0 - 2.0 * tr_yzzz_xxxz[i] * tbe_0 + 6.0 * tr_xyz_xxz[i] - 4.0 * tr_xyz_xxxxz[i] * tke_0 + 3.0 * tr_xyzz_xx[i] - 6.0 * tr_xyzz_xxzz[i] * tke_0 - 2.0 * tr_xyzz_xxxx[i] * tke_0 + 4.0 * tr_xyzz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xyzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxz[i] * tbe_0 - 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xxyy[i] = 2.0 * tr_yz_xxyy[i] - 2.0 * tr_yzz_xxyyz[i] * tke_0 - 2.0 * tr_yzzz_xxyy[i] * tbe_0 + 4.0 * tr_xyz_xyy[i] - 4.0 * tr_xyz_xxxyy[i] * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xxyz[i] = 2.0 * tr_yz_xxyz[i] + tr_yzz_xxy[i] - 2.0 * tr_yzz_xxyzz[i] * tke_0 - 2.0 * tr_yzzz_xxyz[i] * tbe_0 + 4.0 * tr_xyz_xyz[i] - 4.0 * tr_xyz_xxxyz[i] * tke_0 + 2.0 * tr_xyzz_xy[i] - 4.0 * tr_xyzz_xyzz[i] * tke_0 - 2.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tbe_0 - 2.0 * tr_xxyzz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xxzz[i] = 2.0 * tr_yz_xxzz[i] + 2.0 * tr_yzz_xxz[i] - 2.0 * tr_yzz_xxzzz[i] * tke_0 - 2.0 * tr_yzzz_xxzz[i] * tbe_0 + 4.0 * tr_xyz_xzz[i] - 4.0 * tr_xyz_xxxzz[i] * tke_0 + 4.0 * tr_xyzz_xz[i] - 4.0 * tr_xyzz_xzzz[i] * tke_0 - 4.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xyyy[i] = 2.0 * tr_yz_xyyy[i] - 2.0 * tr_yzz_xyyyz[i] * tke_0 - 2.0 * tr_yzzz_xyyy[i] * tbe_0 + 2.0 * tr_xyz_yyy[i] - 4.0 * tr_xyz_xxyyy[i] * tke_0 - 2.0 * tr_xyzz_yyyz[i] * tke_0 + 4.0 * tr_xyzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xyyz[i] = 2.0 * tr_yz_xyyz[i] + tr_yzz_xyy[i] - 2.0 * tr_yzz_xyyzz[i] * tke_0 - 2.0 * tr_yzzz_xyyz[i] * tbe_0 + 2.0 * tr_xyz_yyz[i] - 4.0 * tr_xyz_xxyyz[i] * tke_0 + tr_xyzz_yy[i] - 2.0 * tr_xyzz_yyzz[i] * tke_0 - 2.0 * tr_xyzz_xxyy[i] * tke_0 + 4.0 * tr_xyzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tbe_0 - 2.0 * tr_xxyzz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xyzz[i] = 2.0 * tr_yz_xyzz[i] + 2.0 * tr_yzz_xyz[i] - 2.0 * tr_yzz_xyzzz[i] * tke_0 - 2.0 * tr_yzzz_xyzz[i] * tbe_0 + 2.0 * tr_xyz_yzz[i] - 4.0 * tr_xyz_xxyzz[i] * tke_0 + 2.0 * tr_xyzz_yz[i] - 2.0 * tr_xyzz_yzzz[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyzz[i] * tbe_0 - 4.0 * tr_xxyzz_xyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xzzz[i] = 2.0 * tr_yz_xzzz[i] + 3.0 * tr_yzz_xzz[i] - 2.0 * tr_yzz_xzzzz[i] * tke_0 - 2.0 * tr_yzzz_xzzz[i] * tbe_0 + 2.0 * tr_xyz_zzz[i] - 4.0 * tr_xyz_xxzzz[i] * tke_0 + 3.0 * tr_xyzz_zz[i] - 2.0 * tr_xyzz_zzzz[i] * tke_0 - 6.0 * tr_xyzz_xxzz[i] * tke_0 + 4.0 * tr_xyzz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xzzz[i] * tbe_0 - 6.0 * tr_xxyzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_yyyy[i] = 2.0 * tr_yz_yyyy[i] - 2.0 * tr_yzz_yyyyz[i] * tke_0 - 2.0 * tr_yzzz_yyyy[i] * tbe_0 - 4.0 * tr_xyz_xyyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_yyyz[i] = 2.0 * tr_yz_yyyz[i] + tr_yzz_yyy[i] - 2.0 * tr_yzz_yyyzz[i] * tke_0 - 2.0 * tr_yzzz_yyyz[i] * tbe_0 - 4.0 * tr_xyz_xyyyz[i] * tke_0 - 2.0 * tr_xyzz_xyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyyz[i] * tbe_0 - 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_yyzz[i] = 2.0 * tr_yz_yyzz[i] + 2.0 * tr_yzz_yyz[i] - 2.0 * tr_yzz_yyzzz[i] * tke_0 - 2.0 * tr_yzzz_yyzz[i] * tbe_0 - 4.0 * tr_xyz_xyyzz[i] * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyzz[i] * tbe_0 - 4.0 * tr_xxyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_yzzz[i] = 2.0 * tr_yz_yzzz[i] + 3.0 * tr_yzz_yzz[i] - 2.0 * tr_yzz_yzzzz[i] * tke_0 - 2.0 * tr_yzzz_yzzz[i] * tbe_0 - 4.0 * tr_xyz_xyzzz[i] * tke_0 - 6.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yzzz[i] * tbe_0 - 6.0 * tr_xxyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_zzzz[i] = 2.0 * tr_yz_zzzz[i] + 4.0 * tr_yzz_zzz[i] - 2.0 * tr_yzz_zzzzz[i] * tke_0 - 2.0 * tr_yzzz_zzzz[i] * tbe_0 - 4.0 * tr_xyz_xzzzz[i] * tke_0 - 8.0 * tr_xyzz_xzzz[i] * tke_0 + 4.0 * tr_xyzz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zzzz[i] * tbe_0 - 8.0 * tr_xxyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 585-600 components of targeted buffer : GG

    auto tr_0_0_xz_xzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 585);

    auto tr_0_0_xz_xzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 586);

    auto tr_0_0_xz_xzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 587);

    auto tr_0_0_xz_xzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 588);

    auto tr_0_0_xz_xzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 589);

    auto tr_0_0_xz_xzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 590);

    auto tr_0_0_xz_xzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 591);

    auto tr_0_0_xz_xzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 592);

    auto tr_0_0_xz_xzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 593);

    auto tr_0_0_xz_xzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 594);

    auto tr_0_0_xz_xzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 595);

    auto tr_0_0_xz_xzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 596);

    auto tr_0_0_xz_xzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 597);

    auto tr_0_0_xz_xzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 598);

    auto tr_0_0_xz_xzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 599);

    #pragma omp simd aligned(tr_0_0_xz_xzzz_xxxx, tr_0_0_xz_xzzz_xxxy, tr_0_0_xz_xzzz_xxxz, tr_0_0_xz_xzzz_xxyy, tr_0_0_xz_xzzz_xxyz, tr_0_0_xz_xzzz_xxzz, tr_0_0_xz_xzzz_xyyy, tr_0_0_xz_xzzz_xyyz, tr_0_0_xz_xzzz_xyzz, tr_0_0_xz_xzzz_xzzz, tr_0_0_xz_xzzz_yyyy, tr_0_0_xz_xzzz_yyyz, tr_0_0_xz_xzzz_yyzz, tr_0_0_xz_xzzz_yzzz, tr_0_0_xz_xzzz_zzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxxxz, tr_xxzzz_xxxyz, tr_xxzzz_xxxzz, tr_xxzzz_xxy, tr_xxzzz_xxyyz, tr_xxzzz_xxyzz, tr_xxzzz_xxz, tr_xxzzz_xxzzz, tr_xxzzz_xyy, tr_xxzzz_xyyyz, tr_xxzzz_xyyzz, tr_xxzzz_xyz, tr_xxzzz_xyzzz, tr_xxzzz_xzz, tr_xxzzz_xzzzz, tr_xxzzz_yyy, tr_xxzzz_yyyyz, tr_xxzzz_yyyzz, tr_xxzzz_yyz, tr_xxzzz_yyzzz, tr_xxzzz_yzz, tr_xxzzz_yzzzz, tr_xxzzz_zzz, tr_xxzzz_zzzzz, tr_xxzzzz_xxxx, tr_xxzzzz_xxxy, tr_xxzzzz_xxxz, tr_xxzzzz_xxyy, tr_xxzzzz_xxyz, tr_xxzzzz_xxzz, tr_xxzzzz_xyyy, tr_xxzzzz_xyyz, tr_xxzzzz_xyzz, tr_xxzzzz_xzzz, tr_xxzzzz_yyyy, tr_xxzzzz_yyyz, tr_xxzzzz_yyzz, tr_xxzzzz_yzzz, tr_xxzzzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxxxz, tr_xzzz_xxxxyz, tr_xzzz_xxxxzz, tr_xzzz_xxxy, tr_xzzz_xxxyyz, tr_xzzz_xxxyzz, tr_xzzz_xxxz, tr_xzzz_xxxzzz, tr_xzzz_xxyy, tr_xzzz_xxyyyz, tr_xzzz_xxyyzz, tr_xzzz_xxyz, tr_xzzz_xxyzzz, tr_xzzz_xxzz, tr_xzzz_xxzzzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyyyz, tr_xzzz_xyyyzz, tr_xzzz_xyyz, tr_xzzz_xyyzzz, tr_xzzz_xyzz, tr_xzzz_xyzzzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_xzzzzz, tr_xzzz_yy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxxxx, tr_xzzzz_xxxxy, tr_xzzzz_xxxxz, tr_xzzzz_xxxyy, tr_xzzzz_xxxyz, tr_xzzzz_xxxzz, tr_xzzzz_xxy, tr_xzzzz_xxyyy, tr_xzzzz_xxyyz, tr_xzzzz_xxyzz, tr_xzzzz_xxz, tr_xzzzz_xxzzz, tr_xzzzz_xyy, tr_xzzzz_xyyyy, tr_xzzzz_xyyyz, tr_xzzzz_xyyzz, tr_xzzzz_xyz, tr_xzzzz_xyzzz, tr_xzzzz_xzz, tr_xzzzz_xzzzz, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxxxz, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, tr_zzz_zzzzz, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xzzz, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yzzz, tr_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xzzz_xxxx[i] = 3.0 * tr_zz_xxxx[i] - 2.0 * tr_zzz_xxxxz[i] * tke_0 - 2.0 * tr_zzzz_xxxx[i] * tbe_0 + 12.0 * tr_xzz_xxx[i] - 6.0 * tr_xzz_xxxxx[i] * tke_0 - 8.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xxxy[i] = 3.0 * tr_zz_xxxy[i] - 2.0 * tr_zzz_xxxyz[i] * tke_0 - 2.0 * tr_zzzz_xxxy[i] * tbe_0 + 9.0 * tr_xzz_xxy[i] - 6.0 * tr_xzz_xxxxy[i] * tke_0 - 6.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_xzzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xxxz[i] = 3.0 * tr_zz_xxxz[i] + tr_zzz_xxx[i] - 2.0 * tr_zzz_xxxzz[i] * tke_0 - 2.0 * tr_zzzz_xxxz[i] * tbe_0 + 9.0 * tr_xzz_xxz[i] - 6.0 * tr_xzz_xxxxz[i] * tke_0 + 3.0 * tr_xzzz_xx[i] - 6.0 * tr_xzzz_xxzz[i] * tke_0 - 2.0 * tr_xzzz_xxxx[i] * tke_0 + 4.0 * tr_xzzz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_xzzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzzz_xxxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxxz[i] * tbe_0 - 2.0 * tr_xxzzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xxyy[i] = 3.0 * tr_zz_xxyy[i] - 2.0 * tr_zzz_xxyyz[i] * tke_0 - 2.0 * tr_zzzz_xxyy[i] * tbe_0 + 6.0 * tr_xzz_xyy[i] - 6.0 * tr_xzz_xxxyy[i] * tke_0 - 4.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxyy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xxyz[i] = 3.0 * tr_zz_xxyz[i] + tr_zzz_xxy[i] - 2.0 * tr_zzz_xxyzz[i] * tke_0 - 2.0 * tr_zzzz_xxyz[i] * tbe_0 + 6.0 * tr_xzz_xyz[i] - 6.0 * tr_xzz_xxxyz[i] * tke_0 + 2.0 * tr_xzzz_xy[i] - 4.0 * tr_xzzz_xyzz[i] * tke_0 - 2.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzzz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxyz[i] * tbe_0 - 2.0 * tr_xxzzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xxzz[i] = 3.0 * tr_zz_xxzz[i] + 2.0 * tr_zzz_xxz[i] - 2.0 * tr_zzz_xxzzz[i] * tke_0 - 2.0 * tr_zzzz_xxzz[i] * tbe_0 + 6.0 * tr_xzz_xzz[i] - 6.0 * tr_xzz_xxxzz[i] * tke_0 + 4.0 * tr_xzzz_xz[i] - 4.0 * tr_xzzz_xzzz[i] * tke_0 - 4.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzzz_xxxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxzz[i] * tbe_0 - 4.0 * tr_xxzzz_xxz[i] * tbe_0 + 4.0 * tr_xxzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xyyy[i] = 3.0 * tr_zz_xyyy[i] - 2.0 * tr_zzz_xyyyz[i] * tke_0 - 2.0 * tr_zzzz_xyyy[i] * tbe_0 + 3.0 * tr_xzz_yyy[i] - 6.0 * tr_xzz_xxyyy[i] * tke_0 - 2.0 * tr_xzzz_yyyz[i] * tke_0 + 4.0 * tr_xzzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xyyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xyyz[i] = 3.0 * tr_zz_xyyz[i] + tr_zzz_xyy[i] - 2.0 * tr_zzz_xyyzz[i] * tke_0 - 2.0 * tr_zzzz_xyyz[i] * tbe_0 + 3.0 * tr_xzz_yyz[i] - 6.0 * tr_xzz_xxyyz[i] * tke_0 + tr_xzzz_yy[i] - 2.0 * tr_xzzz_yyzz[i] * tke_0 - 2.0 * tr_xzzz_xxyy[i] * tke_0 + 4.0 * tr_xzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xyyz[i] * tbe_0 - 2.0 * tr_xxzzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xyzz[i] = 3.0 * tr_zz_xyzz[i] + 2.0 * tr_zzz_xyz[i] - 2.0 * tr_zzz_xyzzz[i] * tke_0 - 2.0 * tr_zzzz_xyzz[i] * tbe_0 + 3.0 * tr_xzz_yzz[i] - 6.0 * tr_xzz_xxyzz[i] * tke_0 + 2.0 * tr_xzzz_yz[i] - 2.0 * tr_xzzz_yzzz[i] * tke_0 - 4.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xyzz[i] * tbe_0 - 4.0 * tr_xxzzz_xyz[i] * tbe_0 + 4.0 * tr_xxzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xzzz[i] = 3.0 * tr_zz_xzzz[i] + 3.0 * tr_zzz_xzz[i] - 2.0 * tr_zzz_xzzzz[i] * tke_0 - 2.0 * tr_zzzz_xzzz[i] * tbe_0 + 3.0 * tr_xzz_zzz[i] - 6.0 * tr_xzz_xxzzz[i] * tke_0 + 3.0 * tr_xzzz_zz[i] - 2.0 * tr_xzzz_zzzz[i] * tke_0 - 6.0 * tr_xzzz_xxzz[i] * tke_0 + 4.0 * tr_xzzz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzzz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xzzz[i] * tbe_0 - 6.0 * tr_xxzzz_xzz[i] * tbe_0 + 4.0 * tr_xxzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_yyyy[i] = 3.0 * tr_zz_yyyy[i] - 2.0 * tr_zzz_yyyyz[i] * tke_0 - 2.0 * tr_zzzz_yyyy[i] * tbe_0 - 6.0 * tr_xzz_xyyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yyyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_yyyz[i] = 3.0 * tr_zz_yyyz[i] + tr_zzz_yyy[i] - 2.0 * tr_zzz_yyyzz[i] * tke_0 - 2.0 * tr_zzzz_yyyz[i] * tbe_0 - 6.0 * tr_xzz_xyyyz[i] * tke_0 - 2.0 * tr_xzzz_xyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yyyz[i] * tbe_0 - 2.0 * tr_xxzzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_yyzz[i] = 3.0 * tr_zz_yyzz[i] + 2.0 * tr_zzz_yyz[i] - 2.0 * tr_zzz_yyzzz[i] * tke_0 - 2.0 * tr_zzzz_yyzz[i] * tbe_0 - 6.0 * tr_xzz_xyyzz[i] * tke_0 - 4.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yyzz[i] * tbe_0 - 4.0 * tr_xxzzz_yyz[i] * tbe_0 + 4.0 * tr_xxzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_yzzz[i] = 3.0 * tr_zz_yzzz[i] + 3.0 * tr_zzz_yzz[i] - 2.0 * tr_zzz_yzzzz[i] * tke_0 - 2.0 * tr_zzzz_yzzz[i] * tbe_0 - 6.0 * tr_xzz_xyzzz[i] * tke_0 - 6.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yzzz[i] * tbe_0 - 6.0 * tr_xxzzz_yzz[i] * tbe_0 + 4.0 * tr_xxzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_zzzz[i] = 3.0 * tr_zz_zzzz[i] + 4.0 * tr_zzz_zzz[i] - 2.0 * tr_zzz_zzzzz[i] * tke_0 - 2.0 * tr_zzzz_zzzz[i] * tbe_0 - 6.0 * tr_xzz_xzzzz[i] * tke_0 - 8.0 * tr_xzzz_xzzz[i] * tke_0 + 4.0 * tr_xzzz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_zzzz[i] * tbe_0 - 8.0 * tr_xxzzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 600-615 components of targeted buffer : GG

    auto tr_0_0_xz_yyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 600);

    auto tr_0_0_xz_yyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 601);

    auto tr_0_0_xz_yyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 602);

    auto tr_0_0_xz_yyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 603);

    auto tr_0_0_xz_yyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 604);

    auto tr_0_0_xz_yyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 605);

    auto tr_0_0_xz_yyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 606);

    auto tr_0_0_xz_yyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 607);

    auto tr_0_0_xz_yyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 608);

    auto tr_0_0_xz_yyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 609);

    auto tr_0_0_xz_yyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 610);

    auto tr_0_0_xz_yyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 611);

    auto tr_0_0_xz_yyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 612);

    auto tr_0_0_xz_yyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 613);

    auto tr_0_0_xz_yyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 614);

    #pragma omp simd aligned(tr_0_0_xz_yyyy_xxxx, tr_0_0_xz_yyyy_xxxy, tr_0_0_xz_yyyy_xxxz, tr_0_0_xz_yyyy_xxyy, tr_0_0_xz_yyyy_xxyz, tr_0_0_xz_yyyy_xxzz, tr_0_0_xz_yyyy_xyyy, tr_0_0_xz_yyyy_xyyz, tr_0_0_xz_yyyy_xyzz, tr_0_0_xz_yyyy_xzzz, tr_0_0_xz_yyyy_yyyy, tr_0_0_xz_yyyy_yyyz, tr_0_0_xz_yyyy_yyzz, tr_0_0_xz_yyyy_yzzz, tr_0_0_xz_yyyy_zzzz, tr_xyyyy_xxx, tr_xyyyy_xxxxz, tr_xyyyy_xxxyz, tr_xyyyy_xxxzz, tr_xyyyy_xxy, tr_xyyyy_xxyyz, tr_xyyyy_xxyzz, tr_xyyyy_xxz, tr_xyyyy_xxzzz, tr_xyyyy_xyy, tr_xyyyy_xyyyz, tr_xyyyy_xyyzz, tr_xyyyy_xyz, tr_xyyyy_xyzzz, tr_xyyyy_xzz, tr_xyyyy_xzzzz, tr_xyyyy_yyy, tr_xyyyy_yyyyz, tr_xyyyy_yyyzz, tr_xyyyy_yyz, tr_xyyyy_yyzzz, tr_xyyyy_yzz, tr_xyyyy_yzzzz, tr_xyyyy_zzz, tr_xyyyy_zzzzz, tr_xyyyyz_xxxx, tr_xyyyyz_xxxy, tr_xyyyyz_xxxz, tr_xyyyyz_xxyy, tr_xyyyyz_xxyz, tr_xyyyyz_xxzz, tr_xyyyyz_xyyy, tr_xyyyyz_xyyz, tr_xyyyyz_xyzz, tr_xyyyyz_xzzz, tr_xyyyyz_yyyy, tr_xyyyyz_yyyz, tr_xyyyyz_yyzz, tr_xyyyyz_yzzz, tr_xyyyyz_zzzz, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxxxz, tr_yyyy_xxxxyz, tr_yyyy_xxxxzz, tr_yyyy_xxxy, tr_yyyy_xxxyyz, tr_yyyy_xxxyzz, tr_yyyy_xxxz, tr_yyyy_xxxzzz, tr_yyyy_xxyy, tr_yyyy_xxyyyz, tr_yyyy_xxyyzz, tr_yyyy_xxyz, tr_yyyy_xxyzzz, tr_yyyy_xxzz, tr_yyyy_xxzzzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyyyz, tr_yyyy_xyyyzz, tr_yyyy_xyyz, tr_yyyy_xyyzzz, tr_yyyy_xyzz, tr_yyyy_xyzzzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_xzzzzz, tr_yyyy_yy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxxxx, tr_yyyyz_xxxxy, tr_yyyyz_xxxxz, tr_yyyyz_xxxyy, tr_yyyyz_xxxyz, tr_yyyyz_xxxzz, tr_yyyyz_xxy, tr_yyyyz_xxyyy, tr_yyyyz_xxyyz, tr_yyyyz_xxyzz, tr_yyyyz_xxz, tr_yyyyz_xxzzz, tr_yyyyz_xyy, tr_yyyyz_xyyyy, tr_yyyyz_xyyyz, tr_yyyyz_xyyzz, tr_yyyyz_xyz, tr_yyyyz_xyzzz, tr_yyyyz_xzz, tr_yyyyz_xzzzz, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyyy_xxxx[i] = -8.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xxxy[i] = -6.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xxxz[i] = 3.0 * tr_yyyy_xx[i] - 6.0 * tr_yyyy_xxzz[i] * tke_0 - 2.0 * tr_yyyy_xxxx[i] * tke_0 + 4.0 * tr_yyyy_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_yyyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xxyy[i] = -4.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xxyz[i] = 2.0 * tr_yyyy_xy[i] - 4.0 * tr_yyyy_xyzz[i] * tke_0 - 2.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xxzz[i] = 4.0 * tr_yyyy_xz[i] - 4.0 * tr_yyyy_xzzz[i] * tke_0 - 4.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xyyy[i] = -2.0 * tr_yyyy_yyyz[i] * tke_0 + 4.0 * tr_yyyy_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xyyz[i] = tr_yyyy_yy[i] - 2.0 * tr_yyyy_yyzz[i] * tke_0 - 2.0 * tr_yyyy_xxyy[i] * tke_0 + 4.0 * tr_yyyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xyzz[i] = 2.0 * tr_yyyy_yz[i] - 2.0 * tr_yyyy_yzzz[i] * tke_0 - 4.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xzzz[i] = 3.0 * tr_yyyy_zz[i] - 2.0 * tr_yyyy_zzzz[i] * tke_0 - 6.0 * tr_yyyy_xxzz[i] * tke_0 + 4.0 * tr_yyyy_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyyz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_yyyy[i] = 4.0 * tr_yyyy_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_yyyz[i] = -2.0 * tr_yyyy_xyyy[i] * tke_0 + 4.0 * tr_yyyy_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_yyzz[i] = -4.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_yzzz[i] = -6.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_zzzz[i] = -8.0 * tr_yyyy_xzzz[i] * tke_0 + 4.0 * tr_yyyy_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 615-630 components of targeted buffer : GG

    auto tr_0_0_xz_yyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 615);

    auto tr_0_0_xz_yyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 616);

    auto tr_0_0_xz_yyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 617);

    auto tr_0_0_xz_yyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 618);

    auto tr_0_0_xz_yyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 619);

    auto tr_0_0_xz_yyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 620);

    auto tr_0_0_xz_yyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 621);

    auto tr_0_0_xz_yyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 622);

    auto tr_0_0_xz_yyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 623);

    auto tr_0_0_xz_yyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 624);

    auto tr_0_0_xz_yyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 625);

    auto tr_0_0_xz_yyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 626);

    auto tr_0_0_xz_yyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 627);

    auto tr_0_0_xz_yyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 628);

    auto tr_0_0_xz_yyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 629);

    #pragma omp simd aligned(tr_0_0_xz_yyyz_xxxx, tr_0_0_xz_yyyz_xxxy, tr_0_0_xz_yyyz_xxxz, tr_0_0_xz_yyyz_xxyy, tr_0_0_xz_yyyz_xxyz, tr_0_0_xz_yyyz_xxzz, tr_0_0_xz_yyyz_xyyy, tr_0_0_xz_yyyz_xyyz, tr_0_0_xz_yyyz_xyzz, tr_0_0_xz_yyyz_xzzz, tr_0_0_xz_yyyz_yyyy, tr_0_0_xz_yyyz_yyyz, tr_0_0_xz_yyyz_yyzz, tr_0_0_xz_yyyz_yzzz, tr_0_0_xz_yyyz_zzzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxz, tr_xyyyz_xxxyz, tr_xyyyz_xxxzz, tr_xyyyz_xxy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xxzzz, tr_xyyyz_xyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_xzzzz, tr_xyyyz_yyy, tr_xyyyz_yyyyz, tr_xyyyz_yyyzz, tr_xyyyz_yyz, tr_xyyyz_yyzzz, tr_xyyyz_yzz, tr_xyyyz_yzzzz, tr_xyyyz_zzz, tr_xyyyz_zzzzz, tr_xyyyzz_xxxx, tr_xyyyzz_xxxy, tr_xyyyzz_xxxz, tr_xyyyzz_xxyy, tr_xyyyzz_xxyz, tr_xyyyzz_xxzz, tr_xyyyzz_xyyy, tr_xyyyzz_xyyz, tr_xyyyzz_xyzz, tr_xyyyzz_xzzz, tr_xyyyzz_yyyy, tr_xyyyzz_yyyz, tr_xyyyzz_yyzz, tr_xyyyzz_yzzz, tr_xyyyzz_zzzz, tr_yyy_xxx, tr_yyy_xxxxx, tr_yyy_xxxxy, tr_yyy_xxxxz, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxxxz, tr_yyyz_xxxxyz, tr_yyyz_xxxxzz, tr_yyyz_xxxy, tr_yyyz_xxxyyz, tr_yyyz_xxxyzz, tr_yyyz_xxxz, tr_yyyz_xxxzzz, tr_yyyz_xxyy, tr_yyyz_xxyyyz, tr_yyyz_xxyyzz, tr_yyyz_xxyz, tr_yyyz_xxyzzz, tr_yyyz_xxzz, tr_yyyz_xxzzzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyyyz, tr_yyyz_xyyyzz, tr_yyyz_xyyz, tr_yyyz_xyyzzz, tr_yyyz_xyzz, tr_yyyz_xyzzzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_xzzzzz, tr_yyyz_yy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxxxx, tr_yyyzz_xxxxy, tr_yyyzz_xxxxz, tr_yyyzz_xxxyy, tr_yyyzz_xxxyz, tr_yyyzz_xxxzz, tr_yyyzz_xxy, tr_yyyzz_xxyyy, tr_yyyzz_xxyyz, tr_yyyzz_xxyzz, tr_yyyzz_xxz, tr_yyyzz_xxzzz, tr_yyyzz_xyy, tr_yyyzz_xyyyy, tr_yyyzz_xyyyz, tr_yyyzz_xyyzz, tr_yyyzz_xyz, tr_yyyzz_xyzzz, tr_yyyzz_xzz, tr_yyyzz_xzzzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyyz_xxxx[i] = 4.0 * tr_yyy_xxx[i] - 2.0 * tr_yyy_xxxxx[i] * tke_0 - 8.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xxxy[i] = 3.0 * tr_yyy_xxy[i] - 2.0 * tr_yyy_xxxxy[i] * tke_0 - 6.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyzz_xxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xxxz[i] = 3.0 * tr_yyy_xxz[i] - 2.0 * tr_yyy_xxxxz[i] * tke_0 + 3.0 * tr_yyyz_xx[i] - 6.0 * tr_yyyz_xxzz[i] * tke_0 - 2.0 * tr_yyyz_xxxx[i] * tke_0 + 4.0 * tr_yyyz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_yyyzz_xxz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxxz[i] * tbe_0 - 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xxyy[i] = 2.0 * tr_yyy_xyy[i] - 2.0 * tr_yyy_xxxyy[i] * tke_0 - 4.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xyy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xxyz[i] = 2.0 * tr_yyy_xyz[i] - 2.0 * tr_yyy_xxxyz[i] * tke_0 + 2.0 * tr_yyyz_xy[i] - 4.0 * tr_yyyz_xyzz[i] * tke_0 - 2.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xyz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xxzz[i] = 2.0 * tr_yyy_xzz[i] - 2.0 * tr_yyy_xxxzz[i] * tke_0 + 4.0 * tr_yyyz_xz[i] - 4.0 * tr_yyyz_xzzz[i] * tke_0 - 4.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxzz[i] * tbe_0 - 4.0 * tr_xyyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xyyy[i] = tr_yyy_yyy[i] - 2.0 * tr_yyy_xxyyy[i] * tke_0 - 2.0 * tr_yyyz_yyyz[i] * tke_0 + 4.0 * tr_yyyz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_yyy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xyyz[i] = tr_yyy_yyz[i] - 2.0 * tr_yyy_xxyyz[i] * tke_0 + tr_yyyz_yy[i] - 2.0 * tr_yyyz_yyzz[i] * tke_0 - 2.0 * tr_yyyz_xxyy[i] * tke_0 + 4.0 * tr_yyyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_yyz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyyz[i] * tbe_0 - 2.0 * tr_xyyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xyzz[i] = tr_yyy_yzz[i] - 2.0 * tr_yyy_xxyzz[i] * tke_0 + 2.0 * tr_yyyz_yz[i] - 2.0 * tr_yyyz_yzzz[i] * tke_0 - 4.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_yzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyzz[i] * tbe_0 - 4.0 * tr_xyyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xzzz[i] = tr_yyy_zzz[i] - 2.0 * tr_yyy_xxzzz[i] * tke_0 + 3.0 * tr_yyyz_zz[i] - 2.0 * tr_yyyz_zzzz[i] * tke_0 - 6.0 * tr_yyyz_xxzz[i] * tke_0 + 4.0 * tr_yyyz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_zzz[i] * tbe_0 + 4.0 * tr_yyyzz_xxzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xzzz[i] * tbe_0 - 6.0 * tr_xyyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_yyyy[i] = -2.0 * tr_yyy_xyyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_yyyz[i] = -2.0 * tr_yyy_xyyyz[i] * tke_0 - 2.0 * tr_yyyz_xyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyyz[i] * tbe_0 - 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_yyzz[i] = -2.0 * tr_yyy_xyyzz[i] * tke_0 - 4.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_yzzz[i] = -2.0 * tr_yyy_xyzzz[i] * tke_0 - 6.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yzzz[i] * tbe_0 - 6.0 * tr_xyyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_zzzz[i] = -2.0 * tr_yyy_xzzzz[i] * tke_0 - 8.0 * tr_yyyz_xzzz[i] * tke_0 + 4.0 * tr_yyyz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_zzzz[i] * tbe_0 - 8.0 * tr_xyyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 630-645 components of targeted buffer : GG

    auto tr_0_0_xz_yyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 630);

    auto tr_0_0_xz_yyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 631);

    auto tr_0_0_xz_yyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 632);

    auto tr_0_0_xz_yyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 633);

    auto tr_0_0_xz_yyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 634);

    auto tr_0_0_xz_yyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 635);

    auto tr_0_0_xz_yyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 636);

    auto tr_0_0_xz_yyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 637);

    auto tr_0_0_xz_yyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 638);

    auto tr_0_0_xz_yyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 639);

    auto tr_0_0_xz_yyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 640);

    auto tr_0_0_xz_yyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 641);

    auto tr_0_0_xz_yyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 642);

    auto tr_0_0_xz_yyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 643);

    auto tr_0_0_xz_yyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 644);

    #pragma omp simd aligned(tr_0_0_xz_yyzz_xxxx, tr_0_0_xz_yyzz_xxxy, tr_0_0_xz_yyzz_xxxz, tr_0_0_xz_yyzz_xxyy, tr_0_0_xz_yyzz_xxyz, tr_0_0_xz_yyzz_xxzz, tr_0_0_xz_yyzz_xyyy, tr_0_0_xz_yyzz_xyyz, tr_0_0_xz_yyzz_xyzz, tr_0_0_xz_yyzz_xzzz, tr_0_0_xz_yyzz_yyyy, tr_0_0_xz_yyzz_yyyz, tr_0_0_xz_yyzz_yyzz, tr_0_0_xz_yyzz_yzzz, tr_0_0_xz_yyzz_zzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxz, tr_xyyzz_xxxyz, tr_xyyzz_xxxzz, tr_xyyzz_xxy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xxzzz, tr_xyyzz_xyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_xzzzz, tr_xyyzz_yyy, tr_xyyzz_yyyyz, tr_xyyzz_yyyzz, tr_xyyzz_yyz, tr_xyyzz_yyzzz, tr_xyyzz_yzz, tr_xyyzz_yzzzz, tr_xyyzz_zzz, tr_xyyzz_zzzzz, tr_xyyzzz_xxxx, tr_xyyzzz_xxxy, tr_xyyzzz_xxxz, tr_xyyzzz_xxyy, tr_xyyzzz_xxyz, tr_xyyzzz_xxzz, tr_xyyzzz_xyyy, tr_xyyzzz_xyyz, tr_xyyzzz_xyzz, tr_xyyzzz_xzzz, tr_xyyzzz_yyyy, tr_xyyzzz_yyyz, tr_xyyzzz_yyzz, tr_xyyzzz_yzzz, tr_xyyzzz_zzzz, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxxxz, tr_yyzz_xxxxyz, tr_yyzz_xxxxzz, tr_yyzz_xxxy, tr_yyzz_xxxyyz, tr_yyzz_xxxyzz, tr_yyzz_xxxz, tr_yyzz_xxxzzz, tr_yyzz_xxyy, tr_yyzz_xxyyyz, tr_yyzz_xxyyzz, tr_yyzz_xxyz, tr_yyzz_xxyzzz, tr_yyzz_xxzz, tr_yyzz_xxzzzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyyyz, tr_yyzz_xyyyzz, tr_yyzz_xyyz, tr_yyzz_xyyzzz, tr_yyzz_xyzz, tr_yyzz_xyzzzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_xzzzzz, tr_yyzz_yy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxxxx, tr_yyzzz_xxxxy, tr_yyzzz_xxxxz, tr_yyzzz_xxxyy, tr_yyzzz_xxxyz, tr_yyzzz_xxxzz, tr_yyzzz_xxy, tr_yyzzz_xxyyy, tr_yyzzz_xxyyz, tr_yyzzz_xxyzz, tr_yyzzz_xxz, tr_yyzzz_xxzzz, tr_yyzzz_xyy, tr_yyzzz_xyyyy, tr_yyzzz_xyyyz, tr_yyzzz_xyyzz, tr_yyzzz_xyz, tr_yyzzz_xyzzz, tr_yyzzz_xzz, tr_yyzzz_xzzzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyzz_xxxx[i] = 8.0 * tr_yyz_xxx[i] - 4.0 * tr_yyz_xxxxx[i] * tke_0 - 8.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xxxy[i] = 6.0 * tr_yyz_xxy[i] - 4.0 * tr_yyz_xxxxy[i] * tke_0 - 6.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yyzzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xxxz[i] = 6.0 * tr_yyz_xxz[i] - 4.0 * tr_yyz_xxxxz[i] * tke_0 + 3.0 * tr_yyzz_xx[i] - 6.0 * tr_yyzz_xxzz[i] * tke_0 - 2.0 * tr_yyzz_xxxx[i] * tke_0 + 4.0 * tr_yyzz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_yyzzz_xxz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxxz[i] * tbe_0 - 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xxyy[i] = 4.0 * tr_yyz_xyy[i] - 4.0 * tr_yyz_xxxyy[i] * tke_0 - 4.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xxyz[i] = 4.0 * tr_yyz_xyz[i] - 4.0 * tr_yyz_xxxyz[i] * tke_0 + 2.0 * tr_yyzz_xy[i] - 4.0 * tr_yyzz_xyzz[i] * tke_0 - 2.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xyz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xxzz[i] = 4.0 * tr_yyz_xzz[i] - 4.0 * tr_yyz_xxxzz[i] * tke_0 + 4.0 * tr_yyzz_xz[i] - 4.0 * tr_yyzz_xzzz[i] * tke_0 - 4.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxzz[i] * tbe_0 - 4.0 * tr_xyyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xyyy[i] = 2.0 * tr_yyz_yyy[i] - 4.0 * tr_yyz_xxyyy[i] * tke_0 - 2.0 * tr_yyzz_yyyz[i] * tke_0 + 4.0 * tr_yyzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xyyz[i] = 2.0 * tr_yyz_yyz[i] - 4.0 * tr_yyz_xxyyz[i] * tke_0 + tr_yyzz_yy[i] - 2.0 * tr_yyzz_yyzz[i] * tke_0 - 2.0 * tr_yyzz_xxyy[i] * tke_0 + 4.0 * tr_yyzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_yyz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tbe_0 - 2.0 * tr_xyyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xyzz[i] = 2.0 * tr_yyz_yzz[i] - 4.0 * tr_yyz_xxyzz[i] * tke_0 + 2.0 * tr_yyzz_yz[i] - 2.0 * tr_yyzz_yzzz[i] * tke_0 - 4.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_yzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyzz[i] * tbe_0 - 4.0 * tr_xyyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xzzz[i] = 2.0 * tr_yyz_zzz[i] - 4.0 * tr_yyz_xxzzz[i] * tke_0 + 3.0 * tr_yyzz_zz[i] - 2.0 * tr_yyzz_zzzz[i] * tke_0 - 6.0 * tr_yyzz_xxzz[i] * tke_0 + 4.0 * tr_yyzz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_zzz[i] * tbe_0 + 4.0 * tr_yyzzz_xxzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xzzz[i] * tbe_0 - 6.0 * tr_xyyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_yyyy[i] = -4.0 * tr_yyz_xyyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_yyyz[i] = -4.0 * tr_yyz_xyyyz[i] * tke_0 - 2.0 * tr_yyzz_xyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyyz[i] * tbe_0 - 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_yyzz[i] = -4.0 * tr_yyz_xyyzz[i] * tke_0 - 4.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_yzzz[i] = -4.0 * tr_yyz_xyzzz[i] * tke_0 - 6.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yzzz[i] * tbe_0 - 6.0 * tr_xyyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_zzzz[i] = -4.0 * tr_yyz_xzzzz[i] * tke_0 - 8.0 * tr_yyzz_xzzz[i] * tke_0 + 4.0 * tr_yyzz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_zzzz[i] * tbe_0 - 8.0 * tr_xyyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 645-660 components of targeted buffer : GG

    auto tr_0_0_xz_yzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 645);

    auto tr_0_0_xz_yzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 646);

    auto tr_0_0_xz_yzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 647);

    auto tr_0_0_xz_yzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 648);

    auto tr_0_0_xz_yzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 649);

    auto tr_0_0_xz_yzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 650);

    auto tr_0_0_xz_yzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 651);

    auto tr_0_0_xz_yzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 652);

    auto tr_0_0_xz_yzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 653);

    auto tr_0_0_xz_yzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 654);

    auto tr_0_0_xz_yzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 655);

    auto tr_0_0_xz_yzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 656);

    auto tr_0_0_xz_yzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 657);

    auto tr_0_0_xz_yzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 658);

    auto tr_0_0_xz_yzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 659);

    #pragma omp simd aligned(tr_0_0_xz_yzzz_xxxx, tr_0_0_xz_yzzz_xxxy, tr_0_0_xz_yzzz_xxxz, tr_0_0_xz_yzzz_xxyy, tr_0_0_xz_yzzz_xxyz, tr_0_0_xz_yzzz_xxzz, tr_0_0_xz_yzzz_xyyy, tr_0_0_xz_yzzz_xyyz, tr_0_0_xz_yzzz_xyzz, tr_0_0_xz_yzzz_xzzz, tr_0_0_xz_yzzz_yyyy, tr_0_0_xz_yzzz_yyyz, tr_0_0_xz_yzzz_yyzz, tr_0_0_xz_yzzz_yzzz, tr_0_0_xz_yzzz_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxz, tr_xyzzz_xxxyz, tr_xyzzz_xxxzz, tr_xyzzz_xxy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xxzzz, tr_xyzzz_xyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_xzzzz, tr_xyzzz_yyy, tr_xyzzz_yyyyz, tr_xyzzz_yyyzz, tr_xyzzz_yyz, tr_xyzzz_yyzzz, tr_xyzzz_yzz, tr_xyzzz_yzzzz, tr_xyzzz_zzz, tr_xyzzz_zzzzz, tr_xyzzzz_xxxx, tr_xyzzzz_xxxy, tr_xyzzzz_xxxz, tr_xyzzzz_xxyy, tr_xyzzzz_xxyz, tr_xyzzzz_xxzz, tr_xyzzzz_xyyy, tr_xyzzzz_xyyz, tr_xyzzzz_xyzz, tr_xyzzzz_xzzz, tr_xyzzzz_yyyy, tr_xyzzzz_yyyz, tr_xyzzzz_yyzz, tr_xyzzzz_yzzz, tr_xyzzzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxxxz, tr_yzzz_xxxxyz, tr_yzzz_xxxxzz, tr_yzzz_xxxy, tr_yzzz_xxxyyz, tr_yzzz_xxxyzz, tr_yzzz_xxxz, tr_yzzz_xxxzzz, tr_yzzz_xxyy, tr_yzzz_xxyyyz, tr_yzzz_xxyyzz, tr_yzzz_xxyz, tr_yzzz_xxyzzz, tr_yzzz_xxzz, tr_yzzz_xxzzzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyyyz, tr_yzzz_xyyyzz, tr_yzzz_xyyz, tr_yzzz_xyyzzz, tr_yzzz_xyzz, tr_yzzz_xyzzzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_xzzzzz, tr_yzzz_yy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxxxx, tr_yzzzz_xxxxy, tr_yzzzz_xxxxz, tr_yzzzz_xxxyy, tr_yzzzz_xxxyz, tr_yzzzz_xxxzz, tr_yzzzz_xxy, tr_yzzzz_xxyyy, tr_yzzzz_xxyyz, tr_yzzzz_xxyzz, tr_yzzzz_xxz, tr_yzzzz_xxzzz, tr_yzzzz_xyy, tr_yzzzz_xyyyy, tr_yzzzz_xyyyz, tr_yzzzz_xyyzz, tr_yzzzz_xyz, tr_yzzzz_xyzzz, tr_yzzzz_xzz, tr_yzzzz_xzzzz, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yzzz_xxxx[i] = 12.0 * tr_yzz_xxx[i] - 6.0 * tr_yzz_xxxxx[i] * tke_0 - 8.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xxxy[i] = 9.0 * tr_yzz_xxy[i] - 6.0 * tr_yzz_xxxxy[i] * tke_0 - 6.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_yzzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xxxz[i] = 9.0 * tr_yzz_xxz[i] - 6.0 * tr_yzz_xxxxz[i] * tke_0 + 3.0 * tr_yzzz_xx[i] - 6.0 * tr_yzzz_xxzz[i] * tke_0 - 2.0 * tr_yzzz_xxxx[i] * tke_0 + 4.0 * tr_yzzz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_yzzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxxz[i] * tbe_0 - 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xxyy[i] = 6.0 * tr_yzz_xyy[i] - 6.0 * tr_yzz_xxxyy[i] * tke_0 - 4.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xxyz[i] = 6.0 * tr_yzz_xyz[i] - 6.0 * tr_yzz_xxxyz[i] * tke_0 + 2.0 * tr_yzzz_xy[i] - 4.0 * tr_yzzz_xyzz[i] * tke_0 - 2.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxyz[i] * tbe_0 - 2.0 * tr_xyzzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xxzz[i] = 6.0 * tr_yzz_xzz[i] - 6.0 * tr_yzz_xxxzz[i] * tke_0 + 4.0 * tr_yzzz_xz[i] - 4.0 * tr_yzzz_xzzz[i] * tke_0 - 4.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxzz[i] * tbe_0 - 4.0 * tr_xyzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xyyy[i] = 3.0 * tr_yzz_yyy[i] - 6.0 * tr_yzz_xxyyy[i] * tke_0 - 2.0 * tr_yzzz_yyyz[i] * tke_0 + 4.0 * tr_yzzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xyyz[i] = 3.0 * tr_yzz_yyz[i] - 6.0 * tr_yzz_xxyyz[i] * tke_0 + tr_yzzz_yy[i] - 2.0 * tr_yzzz_yyzz[i] * tke_0 - 2.0 * tr_yzzz_xxyy[i] * tke_0 + 4.0 * tr_yzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyyz[i] * tbe_0 - 2.0 * tr_xyzzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xyzz[i] = 3.0 * tr_yzz_yzz[i] - 6.0 * tr_yzz_xxyzz[i] * tke_0 + 2.0 * tr_yzzz_yz[i] - 2.0 * tr_yzzz_yzzz[i] * tke_0 - 4.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyzz[i] * tbe_0 - 4.0 * tr_xyzzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xzzz[i] = 3.0 * tr_yzz_zzz[i] - 6.0 * tr_yzz_xxzzz[i] * tke_0 + 3.0 * tr_yzzz_zz[i] - 2.0 * tr_yzzz_zzzz[i] * tke_0 - 6.0 * tr_yzzz_xxzz[i] * tke_0 + 4.0 * tr_yzzz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzzz_xxzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xzzz[i] * tbe_0 - 6.0 * tr_xyzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_yyyy[i] = -6.0 * tr_yzz_xyyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_yyyz[i] = -6.0 * tr_yzz_xyyyz[i] * tke_0 - 2.0 * tr_yzzz_xyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyyz[i] * tbe_0 - 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_yyzz[i] = -6.0 * tr_yzz_xyyzz[i] * tke_0 - 4.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyzz[i] * tbe_0 - 4.0 * tr_xyzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_yzzz[i] = -6.0 * tr_yzz_xyzzz[i] * tke_0 - 6.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yzzz[i] * tbe_0 - 6.0 * tr_xyzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_zzzz[i] = -6.0 * tr_yzz_xzzzz[i] * tke_0 - 8.0 * tr_yzzz_xzzz[i] * tke_0 + 4.0 * tr_yzzz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_zzzz[i] * tbe_0 - 8.0 * tr_xyzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 660-675 components of targeted buffer : GG

    auto tr_0_0_xz_zzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 660);

    auto tr_0_0_xz_zzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 661);

    auto tr_0_0_xz_zzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 662);

    auto tr_0_0_xz_zzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 663);

    auto tr_0_0_xz_zzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 664);

    auto tr_0_0_xz_zzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 665);

    auto tr_0_0_xz_zzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 666);

    auto tr_0_0_xz_zzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 667);

    auto tr_0_0_xz_zzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 668);

    auto tr_0_0_xz_zzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 669);

    auto tr_0_0_xz_zzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 670);

    auto tr_0_0_xz_zzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 671);

    auto tr_0_0_xz_zzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 672);

    auto tr_0_0_xz_zzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 673);

    auto tr_0_0_xz_zzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 674);

    #pragma omp simd aligned(tr_0_0_xz_zzzz_xxxx, tr_0_0_xz_zzzz_xxxy, tr_0_0_xz_zzzz_xxxz, tr_0_0_xz_zzzz_xxyy, tr_0_0_xz_zzzz_xxyz, tr_0_0_xz_zzzz_xxzz, tr_0_0_xz_zzzz_xyyy, tr_0_0_xz_zzzz_xyyz, tr_0_0_xz_zzzz_xyzz, tr_0_0_xz_zzzz_xzzz, tr_0_0_xz_zzzz_yyyy, tr_0_0_xz_zzzz_yyyz, tr_0_0_xz_zzzz_yyzz, tr_0_0_xz_zzzz_yzzz, tr_0_0_xz_zzzz_zzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxxxz, tr_xzzzz_xxxyz, tr_xzzzz_xxxzz, tr_xzzzz_xxy, tr_xzzzz_xxyyz, tr_xzzzz_xxyzz, tr_xzzzz_xxz, tr_xzzzz_xxzzz, tr_xzzzz_xyy, tr_xzzzz_xyyyz, tr_xzzzz_xyyzz, tr_xzzzz_xyz, tr_xzzzz_xyzzz, tr_xzzzz_xzz, tr_xzzzz_xzzzz, tr_xzzzz_yyy, tr_xzzzz_yyyyz, tr_xzzzz_yyyzz, tr_xzzzz_yyz, tr_xzzzz_yyzzz, tr_xzzzz_yzz, tr_xzzzz_yzzzz, tr_xzzzz_zzz, tr_xzzzz_zzzzz, tr_xzzzzz_xxxx, tr_xzzzzz_xxxy, tr_xzzzzz_xxxz, tr_xzzzzz_xxyy, tr_xzzzzz_xxyz, tr_xzzzzz_xxzz, tr_xzzzzz_xyyy, tr_xzzzzz_xyyz, tr_xzzzzz_xyzz, tr_xzzzzz_xzzz, tr_xzzzzz_yyyy, tr_xzzzzz_yyyz, tr_xzzzzz_yyzz, tr_xzzzzz_yzzz, tr_xzzzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxx, tr_zzz_xxxxy, tr_zzz_xxxxz, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxxxz, tr_zzzz_xxxxyz, tr_zzzz_xxxxzz, tr_zzzz_xxxy, tr_zzzz_xxxyyz, tr_zzzz_xxxyzz, tr_zzzz_xxxz, tr_zzzz_xxxzzz, tr_zzzz_xxyy, tr_zzzz_xxyyyz, tr_zzzz_xxyyzz, tr_zzzz_xxyz, tr_zzzz_xxyzzz, tr_zzzz_xxzz, tr_zzzz_xxzzzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyyyz, tr_zzzz_xyyyzz, tr_zzzz_xyyz, tr_zzzz_xyyzzz, tr_zzzz_xyzz, tr_zzzz_xyzzzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_xzzzzz, tr_zzzz_yy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzz_zzzz, tr_zzzzz_xxx, tr_zzzzz_xxxxx, tr_zzzzz_xxxxy, tr_zzzzz_xxxxz, tr_zzzzz_xxxyy, tr_zzzzz_xxxyz, tr_zzzzz_xxxzz, tr_zzzzz_xxy, tr_zzzzz_xxyyy, tr_zzzzz_xxyyz, tr_zzzzz_xxyzz, tr_zzzzz_xxz, tr_zzzzz_xxzzz, tr_zzzzz_xyy, tr_zzzzz_xyyyy, tr_zzzzz_xyyyz, tr_zzzzz_xyyzz, tr_zzzzz_xyz, tr_zzzzz_xyzzz, tr_zzzzz_xzz, tr_zzzzz_xzzzz, tr_zzzzz_yyy, tr_zzzzz_yyz, tr_zzzzz_yzz, tr_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_zzzz_xxxx[i] = 16.0 * tr_zzz_xxx[i] - 8.0 * tr_zzz_xxxxx[i] * tke_0 - 8.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxxxz[i] * tke_0 * tke_0 - 8.0 * tr_zzzzz_xxx[i] * tbe_0 + 4.0 * tr_zzzzz_xxxxx[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xxxy[i] = 12.0 * tr_zzz_xxy[i] - 8.0 * tr_zzz_xxxxy[i] * tke_0 - 6.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxxxyz[i] * tke_0 * tke_0 - 6.0 * tr_zzzzz_xxy[i] * tbe_0 + 4.0 * tr_zzzzz_xxxxy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xxxz[i] = 12.0 * tr_zzz_xxz[i] - 8.0 * tr_zzz_xxxxz[i] * tke_0 + 3.0 * tr_zzzz_xx[i] - 6.0 * tr_zzzz_xxzz[i] * tke_0 - 2.0 * tr_zzzz_xxxx[i] * tke_0 + 4.0 * tr_zzzz_xxxxzz[i] * tke_0 * tke_0 - 6.0 * tr_zzzzz_xxz[i] * tbe_0 + 4.0 * tr_zzzzz_xxxxz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxxz[i] * tbe_0 - 2.0 * tr_xzzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xxyy[i] = 8.0 * tr_zzz_xyy[i] - 8.0 * tr_zzz_xxxyy[i] * tke_0 - 4.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xxxyyz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xyy[i] * tbe_0 + 4.0 * tr_zzzzz_xxxyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxyy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xxyz[i] = 8.0 * tr_zzz_xyz[i] - 8.0 * tr_zzz_xxxyz[i] * tke_0 + 2.0 * tr_zzzz_xy[i] - 4.0 * tr_zzzz_xyzz[i] * tke_0 - 2.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxyzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xyz[i] * tbe_0 + 4.0 * tr_zzzzz_xxxyz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxyz[i] * tbe_0 - 2.0 * tr_xzzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xxzz[i] = 8.0 * tr_zzz_xzz[i] - 8.0 * tr_zzz_xxxzz[i] * tke_0 + 4.0 * tr_zzzz_xz[i] - 4.0 * tr_zzzz_xzzz[i] * tke_0 - 4.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xzz[i] * tbe_0 + 4.0 * tr_zzzzz_xxxzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxzz[i] * tbe_0 - 4.0 * tr_xzzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xyyy[i] = 4.0 * tr_zzz_yyy[i] - 8.0 * tr_zzz_xxyyy[i] * tke_0 - 2.0 * tr_zzzz_yyyz[i] * tke_0 + 4.0 * tr_zzzz_xxyyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_yyy[i] * tbe_0 + 4.0 * tr_zzzzz_xxyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xyyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xyyz[i] = 4.0 * tr_zzz_yyz[i] - 8.0 * tr_zzz_xxyyz[i] * tke_0 + tr_zzzz_yy[i] - 2.0 * tr_zzzz_yyzz[i] * tke_0 - 2.0 * tr_zzzz_xxyy[i] * tke_0 + 4.0 * tr_zzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_yyz[i] * tbe_0 + 4.0 * tr_zzzzz_xxyyz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xyyz[i] * tbe_0 - 2.0 * tr_xzzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xyzz[i] = 4.0 * tr_zzz_yzz[i] - 8.0 * tr_zzz_xxyzz[i] * tke_0 + 2.0 * tr_zzzz_yz[i] - 2.0 * tr_zzzz_yzzz[i] * tke_0 - 4.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_yzz[i] * tbe_0 + 4.0 * tr_zzzzz_xxyzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xyzz[i] * tbe_0 - 4.0 * tr_xzzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xzzz[i] = 4.0 * tr_zzz_zzz[i] - 8.0 * tr_zzz_xxzzz[i] * tke_0 + 3.0 * tr_zzzz_zz[i] - 2.0 * tr_zzzz_zzzz[i] * tke_0 - 6.0 * tr_zzzz_xxzz[i] * tke_0 + 4.0 * tr_zzzz_xxzzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_zzz[i] * tbe_0 + 4.0 * tr_zzzzz_xxzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xzzz[i] * tbe_0 - 6.0 * tr_xzzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_yyyy[i] = -8.0 * tr_zzz_xyyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yyyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_yyyz[i] = -8.0 * tr_zzz_xyyyz[i] * tke_0 - 2.0 * tr_zzzz_xyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyyyz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yyyz[i] * tbe_0 - 2.0 * tr_xzzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_yyzz[i] = -8.0 * tr_zzz_xyyzz[i] * tke_0 - 4.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyyzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yyzz[i] * tbe_0 - 4.0 * tr_xzzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_yzzz[i] = -8.0 * tr_zzz_xyzzz[i] * tke_0 - 6.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yzzz[i] * tbe_0 - 6.0 * tr_xzzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_zzzz[i] = -8.0 * tr_zzz_xzzzz[i] * tke_0 - 8.0 * tr_zzzz_xzzz[i] * tke_0 + 4.0 * tr_zzzz_xzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_zzzz[i] * tbe_0 - 8.0 * tr_xzzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 675-690 components of targeted buffer : GG

    auto tr_0_0_yy_xxxx_xxxx = pbuffer.data(idx_op_geom_020_gg + 675);

    auto tr_0_0_yy_xxxx_xxxy = pbuffer.data(idx_op_geom_020_gg + 676);

    auto tr_0_0_yy_xxxx_xxxz = pbuffer.data(idx_op_geom_020_gg + 677);

    auto tr_0_0_yy_xxxx_xxyy = pbuffer.data(idx_op_geom_020_gg + 678);

    auto tr_0_0_yy_xxxx_xxyz = pbuffer.data(idx_op_geom_020_gg + 679);

    auto tr_0_0_yy_xxxx_xxzz = pbuffer.data(idx_op_geom_020_gg + 680);

    auto tr_0_0_yy_xxxx_xyyy = pbuffer.data(idx_op_geom_020_gg + 681);

    auto tr_0_0_yy_xxxx_xyyz = pbuffer.data(idx_op_geom_020_gg + 682);

    auto tr_0_0_yy_xxxx_xyzz = pbuffer.data(idx_op_geom_020_gg + 683);

    auto tr_0_0_yy_xxxx_xzzz = pbuffer.data(idx_op_geom_020_gg + 684);

    auto tr_0_0_yy_xxxx_yyyy = pbuffer.data(idx_op_geom_020_gg + 685);

    auto tr_0_0_yy_xxxx_yyyz = pbuffer.data(idx_op_geom_020_gg + 686);

    auto tr_0_0_yy_xxxx_yyzz = pbuffer.data(idx_op_geom_020_gg + 687);

    auto tr_0_0_yy_xxxx_yzzz = pbuffer.data(idx_op_geom_020_gg + 688);

    auto tr_0_0_yy_xxxx_zzzz = pbuffer.data(idx_op_geom_020_gg + 689);

    #pragma omp simd aligned(tr_0_0_yy_xxxx_xxxx, tr_0_0_yy_xxxx_xxxy, tr_0_0_yy_xxxx_xxxz, tr_0_0_yy_xxxx_xxyy, tr_0_0_yy_xxxx_xxyz, tr_0_0_yy_xxxx_xxzz, tr_0_0_yy_xxxx_xyyy, tr_0_0_yy_xxxx_xyyz, tr_0_0_yy_xxxx_xyzz, tr_0_0_yy_xxxx_xzzz, tr_0_0_yy_xxxx_yyyy, tr_0_0_yy_xxxx_yyyz, tr_0_0_yy_xxxx_yyzz, tr_0_0_yy_xxxx_yzzz, tr_0_0_yy_xxxx_zzzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxxyy, tr_xxxx_xxxy, tr_xxxx_xxxyyy, tr_xxxx_xxxyyz, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyyyy, tr_xxxx_xxyyyz, tr_xxxx_xxyyzz, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyyyy, tr_xxxx_xyyyyz, tr_xxxx_xyyyzz, tr_xxxx_xyyz, tr_xxxx_xyyzzz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyyyy, tr_xxxx_yyyyyz, tr_xxxx_yyyyzz, tr_xxxx_yyyz, tr_xxxx_yyyzzz, tr_xxxx_yyzz, tr_xxxx_yyzzzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxy_xxx, tr_xxxxy_xxxxy, tr_xxxxy_xxxyy, tr_xxxxy_xxxyz, tr_xxxxy_xxy, tr_xxxxy_xxyyy, tr_xxxxy_xxyyz, tr_xxxxy_xxyzz, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyyyy, tr_xxxxy_xyyyz, tr_xxxxy_xyyzz, tr_xxxxy_xyz, tr_xxxxy_xyzzz, tr_xxxxy_xzz, tr_xxxxy_yyy, tr_xxxxy_yyyyy, tr_xxxxy_yyyyz, tr_xxxxy_yyyzz, tr_xxxxy_yyz, tr_xxxxy_yyzzz, tr_xxxxy_yzz, tr_xxxxy_yzzzz, tr_xxxxy_zzz, tr_xxxxyy_xxxx, tr_xxxxyy_xxxy, tr_xxxxyy_xxxz, tr_xxxxyy_xxyy, tr_xxxxyy_xxyz, tr_xxxxyy_xxzz, tr_xxxxyy_xyyy, tr_xxxxyy_xyyz, tr_xxxxyy_xyzz, tr_xxxxyy_xzzz, tr_xxxxyy_yyyy, tr_xxxxyy_yyyz, tr_xxxxyy_yyzz, tr_xxxxyy_yzzz, tr_xxxxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxxx_xxxx[i] = -2.0 * tr_xxxx_xxxx[i] * tbe_0 - 2.0 * tr_xxxx_xxxx[i] * tke_0 + 4.0 * tr_xxxx_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xxxy[i] = -2.0 * tr_xxxx_xxxy[i] * tbe_0 - 6.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_xxx[i] * tbe_0 + 8.0 * tr_xxxxy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xxxz[i] = -2.0 * tr_xxxx_xxxz[i] * tbe_0 - 2.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xxyy[i] = 2.0 * tr_xxxx_xx[i] - 2.0 * tr_xxxx_xxyy[i] * tbe_0 - 10.0 * tr_xxxx_xxyy[i] * tke_0 + 4.0 * tr_xxxx_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_xxy[i] * tbe_0 + 8.0 * tr_xxxxy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xxyz[i] = -2.0 * tr_xxxx_xxyz[i] * tbe_0 - 6.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_xxz[i] * tbe_0 + 8.0 * tr_xxxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xxzz[i] = -2.0 * tr_xxxx_xxzz[i] * tbe_0 - 2.0 * tr_xxxx_xxzz[i] * tke_0 + 4.0 * tr_xxxx_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xyyy[i] = 6.0 * tr_xxxx_xy[i] - 2.0 * tr_xxxx_xyyy[i] * tbe_0 - 14.0 * tr_xxxx_xyyy[i] * tke_0 + 4.0 * tr_xxxx_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxxxy_xyy[i] * tbe_0 + 8.0 * tr_xxxxy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xyyz[i] = 2.0 * tr_xxxx_xz[i] - 2.0 * tr_xxxx_xyyz[i] * tbe_0 - 10.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_xyz[i] * tbe_0 + 8.0 * tr_xxxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xyzz[i] = -2.0 * tr_xxxx_xyzz[i] * tbe_0 - 6.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_xzz[i] * tbe_0 + 8.0 * tr_xxxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xzzz[i] = -2.0 * tr_xxxx_xzzz[i] * tbe_0 - 2.0 * tr_xxxx_xzzz[i] * tke_0 + 4.0 * tr_xxxx_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_yyyy[i] = 12.0 * tr_xxxx_yy[i] - 2.0 * tr_xxxx_yyyy[i] * tbe_0 - 18.0 * tr_xxxx_yyyy[i] * tke_0 + 4.0 * tr_xxxx_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xxxxy_yyy[i] * tbe_0 + 8.0 * tr_xxxxy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_yyyz[i] = 6.0 * tr_xxxx_yz[i] - 2.0 * tr_xxxx_yyyz[i] * tbe_0 - 14.0 * tr_xxxx_yyyz[i] * tke_0 + 4.0 * tr_xxxx_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xxxxy_yyz[i] * tbe_0 + 8.0 * tr_xxxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_yyzz[i] = 2.0 * tr_xxxx_zz[i] - 2.0 * tr_xxxx_yyzz[i] * tbe_0 - 10.0 * tr_xxxx_yyzz[i] * tke_0 + 4.0 * tr_xxxx_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_yzz[i] * tbe_0 + 8.0 * tr_xxxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_yzzz[i] = -2.0 * tr_xxxx_yzzz[i] * tbe_0 - 6.0 * tr_xxxx_yzzz[i] * tke_0 + 4.0 * tr_xxxx_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_zzz[i] * tbe_0 + 8.0 * tr_xxxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_zzzz[i] = -2.0 * tr_xxxx_zzzz[i] * tbe_0 - 2.0 * tr_xxxx_zzzz[i] * tke_0 + 4.0 * tr_xxxx_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 690-705 components of targeted buffer : GG

    auto tr_0_0_yy_xxxy_xxxx = pbuffer.data(idx_op_geom_020_gg + 690);

    auto tr_0_0_yy_xxxy_xxxy = pbuffer.data(idx_op_geom_020_gg + 691);

    auto tr_0_0_yy_xxxy_xxxz = pbuffer.data(idx_op_geom_020_gg + 692);

    auto tr_0_0_yy_xxxy_xxyy = pbuffer.data(idx_op_geom_020_gg + 693);

    auto tr_0_0_yy_xxxy_xxyz = pbuffer.data(idx_op_geom_020_gg + 694);

    auto tr_0_0_yy_xxxy_xxzz = pbuffer.data(idx_op_geom_020_gg + 695);

    auto tr_0_0_yy_xxxy_xyyy = pbuffer.data(idx_op_geom_020_gg + 696);

    auto tr_0_0_yy_xxxy_xyyz = pbuffer.data(idx_op_geom_020_gg + 697);

    auto tr_0_0_yy_xxxy_xyzz = pbuffer.data(idx_op_geom_020_gg + 698);

    auto tr_0_0_yy_xxxy_xzzz = pbuffer.data(idx_op_geom_020_gg + 699);

    auto tr_0_0_yy_xxxy_yyyy = pbuffer.data(idx_op_geom_020_gg + 700);

    auto tr_0_0_yy_xxxy_yyyz = pbuffer.data(idx_op_geom_020_gg + 701);

    auto tr_0_0_yy_xxxy_yyzz = pbuffer.data(idx_op_geom_020_gg + 702);

    auto tr_0_0_yy_xxxy_yzzz = pbuffer.data(idx_op_geom_020_gg + 703);

    auto tr_0_0_yy_xxxy_zzzz = pbuffer.data(idx_op_geom_020_gg + 704);

    #pragma omp simd aligned(tr_0_0_yy_xxxy_xxxx, tr_0_0_yy_xxxy_xxxy, tr_0_0_yy_xxxy_xxxz, tr_0_0_yy_xxxy_xxyy, tr_0_0_yy_xxxy_xxyz, tr_0_0_yy_xxxy_xxzz, tr_0_0_yy_xxxy_xyyy, tr_0_0_yy_xxxy_xyyz, tr_0_0_yy_xxxy_xyzz, tr_0_0_yy_xxxy_xzzz, tr_0_0_yy_xxxy_yyyy, tr_0_0_yy_xxxy_yyyz, tr_0_0_yy_xxxy_yyzz, tr_0_0_yy_xxxy_yzzz, tr_0_0_yy_xxxy_zzzz, tr_xxx_xxx, tr_xxx_xxxxy, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyyyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxxyy, tr_xxxy_xxxy, tr_xxxy_xxxyyy, tr_xxxy_xxxyyz, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyyyy, tr_xxxy_xxyyyz, tr_xxxy_xxyyzz, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyyyy, tr_xxxy_xyyyyz, tr_xxxy_xyyyzz, tr_xxxy_xyyz, tr_xxxy_xyyzzz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyyyy, tr_xxxy_yyyyyz, tr_xxxy_yyyyzz, tr_xxxy_yyyz, tr_xxxy_yyyzzz, tr_xxxy_yyzz, tr_xxxy_yyzzzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyy_xxx, tr_xxxyy_xxxxy, tr_xxxyy_xxxyy, tr_xxxyy_xxxyz, tr_xxxyy_xxy, tr_xxxyy_xxyyy, tr_xxxyy_xxyyz, tr_xxxyy_xxyzz, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyyyy, tr_xxxyy_xyyyz, tr_xxxyy_xyyzz, tr_xxxyy_xyz, tr_xxxyy_xyzzz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyyyy, tr_xxxyy_yyyyz, tr_xxxyy_yyyzz, tr_xxxyy_yyz, tr_xxxyy_yyzzz, tr_xxxyy_yzz, tr_xxxyy_yzzzz, tr_xxxyy_zzz, tr_xxxyyy_xxxx, tr_xxxyyy_xxxy, tr_xxxyyy_xxxz, tr_xxxyyy_xxyy, tr_xxxyyy_xxyz, tr_xxxyyy_xxzz, tr_xxxyyy_xyyy, tr_xxxyyy_xyyz, tr_xxxyyy_xyzz, tr_xxxyyy_xzzz, tr_xxxyyy_yyyy, tr_xxxyyy_yyyz, tr_xxxyyy_yyzz, tr_xxxyyy_yzzz, tr_xxxyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxxy_xxxx[i] = -4.0 * tr_xxx_xxxxy[i] * tke_0 - 6.0 * tr_xxxy_xxxx[i] * tbe_0 - 2.0 * tr_xxxy_xxxx[i] * tke_0 + 4.0 * tr_xxxy_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xxxy[i] = 2.0 * tr_xxx_xxx[i] - 4.0 * tr_xxx_xxxyy[i] * tke_0 - 6.0 * tr_xxxy_xxxy[i] * tbe_0 - 6.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_xxx[i] * tbe_0 + 8.0 * tr_xxxyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xxxz[i] = -4.0 * tr_xxx_xxxyz[i] * tke_0 - 6.0 * tr_xxxy_xxxz[i] * tbe_0 - 2.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xxyy[i] = 4.0 * tr_xxx_xxy[i] - 4.0 * tr_xxx_xxyyy[i] * tke_0 + 2.0 * tr_xxxy_xx[i] - 6.0 * tr_xxxy_xxyy[i] * tbe_0 - 10.0 * tr_xxxy_xxyy[i] * tke_0 + 4.0 * tr_xxxy_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_xxy[i] * tbe_0 + 8.0 * tr_xxxyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xxyz[i] = 2.0 * tr_xxx_xxz[i] - 4.0 * tr_xxx_xxyyz[i] * tke_0 - 6.0 * tr_xxxy_xxyz[i] * tbe_0 - 6.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_xxz[i] * tbe_0 + 8.0 * tr_xxxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xxzz[i] = -4.0 * tr_xxx_xxyzz[i] * tke_0 - 6.0 * tr_xxxy_xxzz[i] * tbe_0 - 2.0 * tr_xxxy_xxzz[i] * tke_0 + 4.0 * tr_xxxy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xyyy[i] = 6.0 * tr_xxx_xyy[i] - 4.0 * tr_xxx_xyyyy[i] * tke_0 + 6.0 * tr_xxxy_xy[i] - 6.0 * tr_xxxy_xyyy[i] * tbe_0 - 14.0 * tr_xxxy_xyyy[i] * tke_0 + 4.0 * tr_xxxy_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxxyy_xyy[i] * tbe_0 + 8.0 * tr_xxxyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xyyz[i] = 4.0 * tr_xxx_xyz[i] - 4.0 * tr_xxx_xyyyz[i] * tke_0 + 2.0 * tr_xxxy_xz[i] - 6.0 * tr_xxxy_xyyz[i] * tbe_0 - 10.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_xyz[i] * tbe_0 + 8.0 * tr_xxxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xyzz[i] = 2.0 * tr_xxx_xzz[i] - 4.0 * tr_xxx_xyyzz[i] * tke_0 - 6.0 * tr_xxxy_xyzz[i] * tbe_0 - 6.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_xzz[i] * tbe_0 + 8.0 * tr_xxxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xzzz[i] = -4.0 * tr_xxx_xyzzz[i] * tke_0 - 6.0 * tr_xxxy_xzzz[i] * tbe_0 - 2.0 * tr_xxxy_xzzz[i] * tke_0 + 4.0 * tr_xxxy_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_yyyy[i] = 8.0 * tr_xxx_yyy[i] - 4.0 * tr_xxx_yyyyy[i] * tke_0 + 12.0 * tr_xxxy_yy[i] - 6.0 * tr_xxxy_yyyy[i] * tbe_0 - 18.0 * tr_xxxy_yyyy[i] * tke_0 + 4.0 * tr_xxxy_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xxxyy_yyy[i] * tbe_0 + 8.0 * tr_xxxyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_yyyz[i] = 6.0 * tr_xxx_yyz[i] - 4.0 * tr_xxx_yyyyz[i] * tke_0 + 6.0 * tr_xxxy_yz[i] - 6.0 * tr_xxxy_yyyz[i] * tbe_0 - 14.0 * tr_xxxy_yyyz[i] * tke_0 + 4.0 * tr_xxxy_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xxxyy_yyz[i] * tbe_0 + 8.0 * tr_xxxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_yyzz[i] = 4.0 * tr_xxx_yzz[i] - 4.0 * tr_xxx_yyyzz[i] * tke_0 + 2.0 * tr_xxxy_zz[i] - 6.0 * tr_xxxy_yyzz[i] * tbe_0 - 10.0 * tr_xxxy_yyzz[i] * tke_0 + 4.0 * tr_xxxy_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_yzz[i] * tbe_0 + 8.0 * tr_xxxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_yzzz[i] = 2.0 * tr_xxx_zzz[i] - 4.0 * tr_xxx_yyzzz[i] * tke_0 - 6.0 * tr_xxxy_yzzz[i] * tbe_0 - 6.0 * tr_xxxy_yzzz[i] * tke_0 + 4.0 * tr_xxxy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_zzz[i] * tbe_0 + 8.0 * tr_xxxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_zzzz[i] = -4.0 * tr_xxx_yzzzz[i] * tke_0 - 6.0 * tr_xxxy_zzzz[i] * tbe_0 - 2.0 * tr_xxxy_zzzz[i] * tke_0 + 4.0 * tr_xxxy_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 705-720 components of targeted buffer : GG

    auto tr_0_0_yy_xxxz_xxxx = pbuffer.data(idx_op_geom_020_gg + 705);

    auto tr_0_0_yy_xxxz_xxxy = pbuffer.data(idx_op_geom_020_gg + 706);

    auto tr_0_0_yy_xxxz_xxxz = pbuffer.data(idx_op_geom_020_gg + 707);

    auto tr_0_0_yy_xxxz_xxyy = pbuffer.data(idx_op_geom_020_gg + 708);

    auto tr_0_0_yy_xxxz_xxyz = pbuffer.data(idx_op_geom_020_gg + 709);

    auto tr_0_0_yy_xxxz_xxzz = pbuffer.data(idx_op_geom_020_gg + 710);

    auto tr_0_0_yy_xxxz_xyyy = pbuffer.data(idx_op_geom_020_gg + 711);

    auto tr_0_0_yy_xxxz_xyyz = pbuffer.data(idx_op_geom_020_gg + 712);

    auto tr_0_0_yy_xxxz_xyzz = pbuffer.data(idx_op_geom_020_gg + 713);

    auto tr_0_0_yy_xxxz_xzzz = pbuffer.data(idx_op_geom_020_gg + 714);

    auto tr_0_0_yy_xxxz_yyyy = pbuffer.data(idx_op_geom_020_gg + 715);

    auto tr_0_0_yy_xxxz_yyyz = pbuffer.data(idx_op_geom_020_gg + 716);

    auto tr_0_0_yy_xxxz_yyzz = pbuffer.data(idx_op_geom_020_gg + 717);

    auto tr_0_0_yy_xxxz_yzzz = pbuffer.data(idx_op_geom_020_gg + 718);

    auto tr_0_0_yy_xxxz_zzzz = pbuffer.data(idx_op_geom_020_gg + 719);

    #pragma omp simd aligned(tr_0_0_yy_xxxz_xxxx, tr_0_0_yy_xxxz_xxxy, tr_0_0_yy_xxxz_xxxz, tr_0_0_yy_xxxz_xxyy, tr_0_0_yy_xxxz_xxyz, tr_0_0_yy_xxxz_xxzz, tr_0_0_yy_xxxz_xyyy, tr_0_0_yy_xxxz_xyyz, tr_0_0_yy_xxxz_xyzz, tr_0_0_yy_xxxz_xzzz, tr_0_0_yy_xxxz_yyyy, tr_0_0_yy_xxxz_yyyz, tr_0_0_yy_xxxz_yyzz, tr_0_0_yy_xxxz_yzzz, tr_0_0_yy_xxxz_zzzz, tr_xxxyyz_xxxx, tr_xxxyyz_xxxy, tr_xxxyyz_xxxz, tr_xxxyyz_xxyy, tr_xxxyyz_xxyz, tr_xxxyyz_xxzz, tr_xxxyyz_xyyy, tr_xxxyyz_xyyz, tr_xxxyyz_xyzz, tr_xxxyyz_xzzz, tr_xxxyyz_yyyy, tr_xxxyyz_yyyz, tr_xxxyyz_yyzz, tr_xxxyyz_yzzz, tr_xxxyyz_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxy, tr_xxxyz_xxxyy, tr_xxxyz_xxxyz, tr_xxxyz_xxy, tr_xxxyz_xxyyy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyyyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyyyy, tr_xxxyz_yyyyz, tr_xxxyz_yyyzz, tr_xxxyz_yyz, tr_xxxyz_yyzzz, tr_xxxyz_yzz, tr_xxxyz_yzzzz, tr_xxxyz_zzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxxyy, tr_xxxz_xxxy, tr_xxxz_xxxyyy, tr_xxxz_xxxyyz, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyyyy, tr_xxxz_xxyyyz, tr_xxxz_xxyyzz, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyyyy, tr_xxxz_xyyyyz, tr_xxxz_xyyyzz, tr_xxxz_xyyz, tr_xxxz_xyyzzz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyyyy, tr_xxxz_yyyyyz, tr_xxxz_yyyyzz, tr_xxxz_yyyz, tr_xxxz_yyyzzz, tr_xxxz_yyzz, tr_xxxz_yyzzzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxxz_xxxx[i] = -2.0 * tr_xxxz_xxxx[i] * tbe_0 - 2.0 * tr_xxxz_xxxx[i] * tke_0 + 4.0 * tr_xxxz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xxxy[i] = -2.0 * tr_xxxz_xxxy[i] * tbe_0 - 6.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xxx[i] * tbe_0 + 8.0 * tr_xxxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xxxz[i] = -2.0 * tr_xxxz_xxxz[i] * tbe_0 - 2.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xxyy[i] = 2.0 * tr_xxxz_xx[i] - 2.0 * tr_xxxz_xxyy[i] * tbe_0 - 10.0 * tr_xxxz_xxyy[i] * tke_0 + 4.0 * tr_xxxz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xxy[i] * tbe_0 + 8.0 * tr_xxxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xxyz[i] = -2.0 * tr_xxxz_xxyz[i] * tbe_0 - 6.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xxz[i] * tbe_0 + 8.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xxzz[i] = -2.0 * tr_xxxz_xxzz[i] * tbe_0 - 2.0 * tr_xxxz_xxzz[i] * tke_0 + 4.0 * tr_xxxz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xyyy[i] = 6.0 * tr_xxxz_xy[i] - 2.0 * tr_xxxz_xyyy[i] * tbe_0 - 14.0 * tr_xxxz_xyyy[i] * tke_0 + 4.0 * tr_xxxz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxxyz_xyy[i] * tbe_0 + 8.0 * tr_xxxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xyyz[i] = 2.0 * tr_xxxz_xz[i] - 2.0 * tr_xxxz_xyyz[i] * tbe_0 - 10.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xyz[i] * tbe_0 + 8.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xyzz[i] = -2.0 * tr_xxxz_xyzz[i] * tbe_0 - 6.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xzz[i] * tbe_0 + 8.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xzzz[i] = -2.0 * tr_xxxz_xzzz[i] * tbe_0 - 2.0 * tr_xxxz_xzzz[i] * tke_0 + 4.0 * tr_xxxz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_yyyy[i] = 12.0 * tr_xxxz_yy[i] - 2.0 * tr_xxxz_yyyy[i] * tbe_0 - 18.0 * tr_xxxz_yyyy[i] * tke_0 + 4.0 * tr_xxxz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xxxyz_yyy[i] * tbe_0 + 8.0 * tr_xxxyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_yyyz[i] = 6.0 * tr_xxxz_yz[i] - 2.0 * tr_xxxz_yyyz[i] * tbe_0 - 14.0 * tr_xxxz_yyyz[i] * tke_0 + 4.0 * tr_xxxz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xxxyz_yyz[i] * tbe_0 + 8.0 * tr_xxxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_yyzz[i] = 2.0 * tr_xxxz_zz[i] - 2.0 * tr_xxxz_yyzz[i] * tbe_0 - 10.0 * tr_xxxz_yyzz[i] * tke_0 + 4.0 * tr_xxxz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_yzz[i] * tbe_0 + 8.0 * tr_xxxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_yzzz[i] = -2.0 * tr_xxxz_yzzz[i] * tbe_0 - 6.0 * tr_xxxz_yzzz[i] * tke_0 + 4.0 * tr_xxxz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_zzz[i] * tbe_0 + 8.0 * tr_xxxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_zzzz[i] = -2.0 * tr_xxxz_zzzz[i] * tbe_0 - 2.0 * tr_xxxz_zzzz[i] * tke_0 + 4.0 * tr_xxxz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 720-735 components of targeted buffer : GG

    auto tr_0_0_yy_xxyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 720);

    auto tr_0_0_yy_xxyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 721);

    auto tr_0_0_yy_xxyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 722);

    auto tr_0_0_yy_xxyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 723);

    auto tr_0_0_yy_xxyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 724);

    auto tr_0_0_yy_xxyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 725);

    auto tr_0_0_yy_xxyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 726);

    auto tr_0_0_yy_xxyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 727);

    auto tr_0_0_yy_xxyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 728);

    auto tr_0_0_yy_xxyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 729);

    auto tr_0_0_yy_xxyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 730);

    auto tr_0_0_yy_xxyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 731);

    auto tr_0_0_yy_xxyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 732);

    auto tr_0_0_yy_xxyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 733);

    auto tr_0_0_yy_xxyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 734);

    #pragma omp simd aligned(tr_0_0_yy_xxyy_xxxx, tr_0_0_yy_xxyy_xxxy, tr_0_0_yy_xxyy_xxxz, tr_0_0_yy_xxyy_xxyy, tr_0_0_yy_xxyy_xxyz, tr_0_0_yy_xxyy_xxzz, tr_0_0_yy_xxyy_xyyy, tr_0_0_yy_xxyy_xyyz, tr_0_0_yy_xxyy_xyzz, tr_0_0_yy_xxyy_xzzz, tr_0_0_yy_xxyy_yyyy, tr_0_0_yy_xxyy_yyyz, tr_0_0_yy_xxyy_yyzz, tr_0_0_yy_xxyy_yzzz, tr_0_0_yy_xxyy_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxxyy, tr_xxyy_xxxy, tr_xxyy_xxxyyy, tr_xxyy_xxxyyz, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyyyy, tr_xxyy_xxyyyz, tr_xxyy_xxyyzz, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyyyy, tr_xxyy_xyyyyz, tr_xxyy_xyyyzz, tr_xxyy_xyyz, tr_xxyy_xyyzzz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyyyy, tr_xxyy_yyyyyz, tr_xxyy_yyyyzz, tr_xxyy_yyyz, tr_xxyy_yyyzzz, tr_xxyy_yyzz, tr_xxyy_yyzzzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyy_xxx, tr_xxyyy_xxxxy, tr_xxyyy_xxxyy, tr_xxyyy_xxxyz, tr_xxyyy_xxy, tr_xxyyy_xxyyy, tr_xxyyy_xxyyz, tr_xxyyy_xxyzz, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyyyy, tr_xxyyy_xyyyz, tr_xxyyy_xyyzz, tr_xxyyy_xyz, tr_xxyyy_xyzzz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyyyy, tr_xxyyy_yyyyz, tr_xxyyy_yyyzz, tr_xxyyy_yyz, tr_xxyyy_yyzzz, tr_xxyyy_yzz, tr_xxyyy_yzzzz, tr_xxyyy_zzz, tr_xxyyyy_xxxx, tr_xxyyyy_xxxy, tr_xxyyyy_xxxz, tr_xxyyyy_xxyy, tr_xxyyyy_xxyz, tr_xxyyyy_xxzz, tr_xxyyyy_xyyy, tr_xxyyyy_xyyz, tr_xxyyyy_xyzz, tr_xxyyyy_xzzz, tr_xxyyyy_yyyy, tr_xxyyyy_yyyz, tr_xxyyyy_yyzz, tr_xxyyyy_yzzz, tr_xxyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxyy_xxxx[i] = 2.0 * tr_xx_xxxx[i] - 8.0 * tr_xxy_xxxxy[i] * tke_0 - 10.0 * tr_xxyy_xxxx[i] * tbe_0 - 2.0 * tr_xxyy_xxxx[i] * tke_0 + 4.0 * tr_xxyy_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xxxy[i] = 2.0 * tr_xx_xxxy[i] + 4.0 * tr_xxy_xxx[i] - 8.0 * tr_xxy_xxxyy[i] * tke_0 - 10.0 * tr_xxyy_xxxy[i] * tbe_0 - 6.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_xxx[i] * tbe_0 + 8.0 * tr_xxyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xxxz[i] = 2.0 * tr_xx_xxxz[i] - 8.0 * tr_xxy_xxxyz[i] * tke_0 - 10.0 * tr_xxyy_xxxz[i] * tbe_0 - 2.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xxyy[i] = 2.0 * tr_xx_xxyy[i] + 8.0 * tr_xxy_xxy[i] - 8.0 * tr_xxy_xxyyy[i] * tke_0 + 2.0 * tr_xxyy_xx[i] - 10.0 * tr_xxyy_xxyy[i] * tbe_0 - 10.0 * tr_xxyy_xxyy[i] * tke_0 + 4.0 * tr_xxyy_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_xxy[i] * tbe_0 + 8.0 * tr_xxyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xxyz[i] = 2.0 * tr_xx_xxyz[i] + 4.0 * tr_xxy_xxz[i] - 8.0 * tr_xxy_xxyyz[i] * tke_0 - 10.0 * tr_xxyy_xxyz[i] * tbe_0 - 6.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_xxz[i] * tbe_0 + 8.0 * tr_xxyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xxzz[i] = 2.0 * tr_xx_xxzz[i] - 8.0 * tr_xxy_xxyzz[i] * tke_0 - 10.0 * tr_xxyy_xxzz[i] * tbe_0 - 2.0 * tr_xxyy_xxzz[i] * tke_0 + 4.0 * tr_xxyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xyyy[i] = 2.0 * tr_xx_xyyy[i] + 12.0 * tr_xxy_xyy[i] - 8.0 * tr_xxy_xyyyy[i] * tke_0 + 6.0 * tr_xxyy_xy[i] - 10.0 * tr_xxyy_xyyy[i] * tbe_0 - 14.0 * tr_xxyy_xyyy[i] * tke_0 + 4.0 * tr_xxyy_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxyyy_xyy[i] * tbe_0 + 8.0 * tr_xxyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xyyz[i] = 2.0 * tr_xx_xyyz[i] + 8.0 * tr_xxy_xyz[i] - 8.0 * tr_xxy_xyyyz[i] * tke_0 + 2.0 * tr_xxyy_xz[i] - 10.0 * tr_xxyy_xyyz[i] * tbe_0 - 10.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_xyz[i] * tbe_0 + 8.0 * tr_xxyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xyzz[i] = 2.0 * tr_xx_xyzz[i] + 4.0 * tr_xxy_xzz[i] - 8.0 * tr_xxy_xyyzz[i] * tke_0 - 10.0 * tr_xxyy_xyzz[i] * tbe_0 - 6.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_xzz[i] * tbe_0 + 8.0 * tr_xxyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xzzz[i] = 2.0 * tr_xx_xzzz[i] - 8.0 * tr_xxy_xyzzz[i] * tke_0 - 10.0 * tr_xxyy_xzzz[i] * tbe_0 - 2.0 * tr_xxyy_xzzz[i] * tke_0 + 4.0 * tr_xxyy_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_yyyy[i] = 2.0 * tr_xx_yyyy[i] + 16.0 * tr_xxy_yyy[i] - 8.0 * tr_xxy_yyyyy[i] * tke_0 + 12.0 * tr_xxyy_yy[i] - 10.0 * tr_xxyy_yyyy[i] * tbe_0 - 18.0 * tr_xxyy_yyyy[i] * tke_0 + 4.0 * tr_xxyy_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xxyyy_yyy[i] * tbe_0 + 8.0 * tr_xxyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_yyyz[i] = 2.0 * tr_xx_yyyz[i] + 12.0 * tr_xxy_yyz[i] - 8.0 * tr_xxy_yyyyz[i] * tke_0 + 6.0 * tr_xxyy_yz[i] - 10.0 * tr_xxyy_yyyz[i] * tbe_0 - 14.0 * tr_xxyy_yyyz[i] * tke_0 + 4.0 * tr_xxyy_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xxyyy_yyz[i] * tbe_0 + 8.0 * tr_xxyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_yyzz[i] = 2.0 * tr_xx_yyzz[i] + 8.0 * tr_xxy_yzz[i] - 8.0 * tr_xxy_yyyzz[i] * tke_0 + 2.0 * tr_xxyy_zz[i] - 10.0 * tr_xxyy_yyzz[i] * tbe_0 - 10.0 * tr_xxyy_yyzz[i] * tke_0 + 4.0 * tr_xxyy_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_yzz[i] * tbe_0 + 8.0 * tr_xxyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_yzzz[i] = 2.0 * tr_xx_yzzz[i] + 4.0 * tr_xxy_zzz[i] - 8.0 * tr_xxy_yyzzz[i] * tke_0 - 10.0 * tr_xxyy_yzzz[i] * tbe_0 - 6.0 * tr_xxyy_yzzz[i] * tke_0 + 4.0 * tr_xxyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_zzz[i] * tbe_0 + 8.0 * tr_xxyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_zzzz[i] = 2.0 * tr_xx_zzzz[i] - 8.0 * tr_xxy_yzzzz[i] * tke_0 - 10.0 * tr_xxyy_zzzz[i] * tbe_0 - 2.0 * tr_xxyy_zzzz[i] * tke_0 + 4.0 * tr_xxyy_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 735-750 components of targeted buffer : GG

    auto tr_0_0_yy_xxyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 735);

    auto tr_0_0_yy_xxyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 736);

    auto tr_0_0_yy_xxyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 737);

    auto tr_0_0_yy_xxyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 738);

    auto tr_0_0_yy_xxyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 739);

    auto tr_0_0_yy_xxyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 740);

    auto tr_0_0_yy_xxyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 741);

    auto tr_0_0_yy_xxyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 742);

    auto tr_0_0_yy_xxyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 743);

    auto tr_0_0_yy_xxyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 744);

    auto tr_0_0_yy_xxyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 745);

    auto tr_0_0_yy_xxyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 746);

    auto tr_0_0_yy_xxyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 747);

    auto tr_0_0_yy_xxyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 748);

    auto tr_0_0_yy_xxyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 749);

    #pragma omp simd aligned(tr_0_0_yy_xxyz_xxxx, tr_0_0_yy_xxyz_xxxy, tr_0_0_yy_xxyz_xxxz, tr_0_0_yy_xxyz_xxyy, tr_0_0_yy_xxyz_xxyz, tr_0_0_yy_xxyz_xxzz, tr_0_0_yy_xxyz_xyyy, tr_0_0_yy_xxyz_xyyz, tr_0_0_yy_xxyz_xyzz, tr_0_0_yy_xxyz_xzzz, tr_0_0_yy_xxyz_yyyy, tr_0_0_yy_xxyz_yyyz, tr_0_0_yy_xxyz_yyzz, tr_0_0_yy_xxyz_yzzz, tr_0_0_yy_xxyz_zzzz, tr_xxyyyz_xxxx, tr_xxyyyz_xxxy, tr_xxyyyz_xxxz, tr_xxyyyz_xxyy, tr_xxyyyz_xxyz, tr_xxyyyz_xxzz, tr_xxyyyz_xyyy, tr_xxyyyz_xyyz, tr_xxyyyz_xyzz, tr_xxyyyz_xzzz, tr_xxyyyz_yyyy, tr_xxyyyz_yyyz, tr_xxyyyz_yyzz, tr_xxyyyz_yzzz, tr_xxyyyz_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxy, tr_xxyyz_xxxyy, tr_xxyyz_xxxyz, tr_xxyyz_xxy, tr_xxyyz_xxyyy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyyyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyyyy, tr_xxyyz_yyyyz, tr_xxyyz_yyyzz, tr_xxyyz_yyz, tr_xxyyz_yyzzz, tr_xxyyz_yzz, tr_xxyyz_yzzzz, tr_xxyyz_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxxyy, tr_xxyz_xxxy, tr_xxyz_xxxyyy, tr_xxyz_xxxyyz, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyyyy, tr_xxyz_xxyyyz, tr_xxyz_xxyyzz, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyyyy, tr_xxyz_xyyyyz, tr_xxyz_xyyyzz, tr_xxyz_xyyz, tr_xxyz_xyyzzz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyyyy, tr_xxyz_yyyyyz, tr_xxyz_yyyyzz, tr_xxyz_yyyz, tr_xxyz_yyyzzz, tr_xxyz_yyzz, tr_xxyz_yyzzzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxyz_xxxx[i] = -4.0 * tr_xxz_xxxxy[i] * tke_0 - 6.0 * tr_xxyz_xxxx[i] * tbe_0 - 2.0 * tr_xxyz_xxxx[i] * tke_0 + 4.0 * tr_xxyz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xxxy[i] = 2.0 * tr_xxz_xxx[i] - 4.0 * tr_xxz_xxxyy[i] * tke_0 - 6.0 * tr_xxyz_xxxy[i] * tbe_0 - 6.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xxx[i] * tbe_0 + 8.0 * tr_xxyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xxxz[i] = -4.0 * tr_xxz_xxxyz[i] * tke_0 - 6.0 * tr_xxyz_xxxz[i] * tbe_0 - 2.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xxyy[i] = 4.0 * tr_xxz_xxy[i] - 4.0 * tr_xxz_xxyyy[i] * tke_0 + 2.0 * tr_xxyz_xx[i] - 6.0 * tr_xxyz_xxyy[i] * tbe_0 - 10.0 * tr_xxyz_xxyy[i] * tke_0 + 4.0 * tr_xxyz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xxy[i] * tbe_0 + 8.0 * tr_xxyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xxyz[i] = 2.0 * tr_xxz_xxz[i] - 4.0 * tr_xxz_xxyyz[i] * tke_0 - 6.0 * tr_xxyz_xxyz[i] * tbe_0 - 6.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xxz[i] * tbe_0 + 8.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xxzz[i] = -4.0 * tr_xxz_xxyzz[i] * tke_0 - 6.0 * tr_xxyz_xxzz[i] * tbe_0 - 2.0 * tr_xxyz_xxzz[i] * tke_0 + 4.0 * tr_xxyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xyyy[i] = 6.0 * tr_xxz_xyy[i] - 4.0 * tr_xxz_xyyyy[i] * tke_0 + 6.0 * tr_xxyz_xy[i] - 6.0 * tr_xxyz_xyyy[i] * tbe_0 - 14.0 * tr_xxyz_xyyy[i] * tke_0 + 4.0 * tr_xxyz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxyyz_xyy[i] * tbe_0 + 8.0 * tr_xxyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xyyz[i] = 4.0 * tr_xxz_xyz[i] - 4.0 * tr_xxz_xyyyz[i] * tke_0 + 2.0 * tr_xxyz_xz[i] - 6.0 * tr_xxyz_xyyz[i] * tbe_0 - 10.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xyz[i] * tbe_0 + 8.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xyzz[i] = 2.0 * tr_xxz_xzz[i] - 4.0 * tr_xxz_xyyzz[i] * tke_0 - 6.0 * tr_xxyz_xyzz[i] * tbe_0 - 6.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xzz[i] * tbe_0 + 8.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xzzz[i] = -4.0 * tr_xxz_xyzzz[i] * tke_0 - 6.0 * tr_xxyz_xzzz[i] * tbe_0 - 2.0 * tr_xxyz_xzzz[i] * tke_0 + 4.0 * tr_xxyz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_yyyy[i] = 8.0 * tr_xxz_yyy[i] - 4.0 * tr_xxz_yyyyy[i] * tke_0 + 12.0 * tr_xxyz_yy[i] - 6.0 * tr_xxyz_yyyy[i] * tbe_0 - 18.0 * tr_xxyz_yyyy[i] * tke_0 + 4.0 * tr_xxyz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xxyyz_yyy[i] * tbe_0 + 8.0 * tr_xxyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_yyyz[i] = 6.0 * tr_xxz_yyz[i] - 4.0 * tr_xxz_yyyyz[i] * tke_0 + 6.0 * tr_xxyz_yz[i] - 6.0 * tr_xxyz_yyyz[i] * tbe_0 - 14.0 * tr_xxyz_yyyz[i] * tke_0 + 4.0 * tr_xxyz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xxyyz_yyz[i] * tbe_0 + 8.0 * tr_xxyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_yyzz[i] = 4.0 * tr_xxz_yzz[i] - 4.0 * tr_xxz_yyyzz[i] * tke_0 + 2.0 * tr_xxyz_zz[i] - 6.0 * tr_xxyz_yyzz[i] * tbe_0 - 10.0 * tr_xxyz_yyzz[i] * tke_0 + 4.0 * tr_xxyz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_yzz[i] * tbe_0 + 8.0 * tr_xxyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_yzzz[i] = 2.0 * tr_xxz_zzz[i] - 4.0 * tr_xxz_yyzzz[i] * tke_0 - 6.0 * tr_xxyz_yzzz[i] * tbe_0 - 6.0 * tr_xxyz_yzzz[i] * tke_0 + 4.0 * tr_xxyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_zzz[i] * tbe_0 + 8.0 * tr_xxyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_zzzz[i] = -4.0 * tr_xxz_yzzzz[i] * tke_0 - 6.0 * tr_xxyz_zzzz[i] * tbe_0 - 2.0 * tr_xxyz_zzzz[i] * tke_0 + 4.0 * tr_xxyz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 750-765 components of targeted buffer : GG

    auto tr_0_0_yy_xxzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 750);

    auto tr_0_0_yy_xxzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 751);

    auto tr_0_0_yy_xxzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 752);

    auto tr_0_0_yy_xxzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 753);

    auto tr_0_0_yy_xxzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 754);

    auto tr_0_0_yy_xxzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 755);

    auto tr_0_0_yy_xxzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 756);

    auto tr_0_0_yy_xxzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 757);

    auto tr_0_0_yy_xxzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 758);

    auto tr_0_0_yy_xxzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 759);

    auto tr_0_0_yy_xxzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 760);

    auto tr_0_0_yy_xxzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 761);

    auto tr_0_0_yy_xxzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 762);

    auto tr_0_0_yy_xxzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 763);

    auto tr_0_0_yy_xxzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 764);

    #pragma omp simd aligned(tr_0_0_yy_xxzz_xxxx, tr_0_0_yy_xxzz_xxxy, tr_0_0_yy_xxzz_xxxz, tr_0_0_yy_xxzz_xxyy, tr_0_0_yy_xxzz_xxyz, tr_0_0_yy_xxzz_xxzz, tr_0_0_yy_xxzz_xyyy, tr_0_0_yy_xxzz_xyyz, tr_0_0_yy_xxzz_xyzz, tr_0_0_yy_xxzz_xzzz, tr_0_0_yy_xxzz_yyyy, tr_0_0_yy_xxzz_yyyz, tr_0_0_yy_xxzz_yyzz, tr_0_0_yy_xxzz_yzzz, tr_0_0_yy_xxzz_zzzz, tr_xxyyzz_xxxx, tr_xxyyzz_xxxy, tr_xxyyzz_xxxz, tr_xxyyzz_xxyy, tr_xxyyzz_xxyz, tr_xxyyzz_xxzz, tr_xxyyzz_xyyy, tr_xxyyzz_xyyz, tr_xxyyzz_xyzz, tr_xxyyzz_xzzz, tr_xxyyzz_yyyy, tr_xxyyzz_yyyz, tr_xxyyzz_yyzz, tr_xxyyzz_yzzz, tr_xxyyzz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxy, tr_xxyzz_xxxyy, tr_xxyzz_xxxyz, tr_xxyzz_xxy, tr_xxyzz_xxyyy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyyyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyyyy, tr_xxyzz_yyyyz, tr_xxyzz_yyyzz, tr_xxyzz_yyz, tr_xxyzz_yyzzz, tr_xxyzz_yzz, tr_xxyzz_yzzzz, tr_xxyzz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxxyy, tr_xxzz_xxxy, tr_xxzz_xxxyyy, tr_xxzz_xxxyyz, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyyyy, tr_xxzz_xxyyyz, tr_xxzz_xxyyzz, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyyyy, tr_xxzz_xyyyyz, tr_xxzz_xyyyzz, tr_xxzz_xyyz, tr_xxzz_xyyzzz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyyyy, tr_xxzz_yyyyyz, tr_xxzz_yyyyzz, tr_xxzz_yyyz, tr_xxzz_yyyzzz, tr_xxzz_yyzz, tr_xxzz_yyzzzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxzz_xxxx[i] = -2.0 * tr_xxzz_xxxx[i] * tbe_0 - 2.0 * tr_xxzz_xxxx[i] * tke_0 + 4.0 * tr_xxzz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xxxy[i] = -2.0 * tr_xxzz_xxxy[i] * tbe_0 - 6.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xxx[i] * tbe_0 + 8.0 * tr_xxyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xxxz[i] = -2.0 * tr_xxzz_xxxz[i] * tbe_0 - 2.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xxyy[i] = 2.0 * tr_xxzz_xx[i] - 2.0 * tr_xxzz_xxyy[i] * tbe_0 - 10.0 * tr_xxzz_xxyy[i] * tke_0 + 4.0 * tr_xxzz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xxy[i] * tbe_0 + 8.0 * tr_xxyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xxyz[i] = -2.0 * tr_xxzz_xxyz[i] * tbe_0 - 6.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xxz[i] * tbe_0 + 8.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xxzz[i] = -2.0 * tr_xxzz_xxzz[i] * tbe_0 - 2.0 * tr_xxzz_xxzz[i] * tke_0 + 4.0 * tr_xxzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xyyy[i] = 6.0 * tr_xxzz_xy[i] - 2.0 * tr_xxzz_xyyy[i] * tbe_0 - 14.0 * tr_xxzz_xyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xxyzz_xyy[i] * tbe_0 + 8.0 * tr_xxyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xyyz[i] = 2.0 * tr_xxzz_xz[i] - 2.0 * tr_xxzz_xyyz[i] * tbe_0 - 10.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xyz[i] * tbe_0 + 8.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xyzz[i] = -2.0 * tr_xxzz_xyzz[i] * tbe_0 - 6.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xzz[i] * tbe_0 + 8.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xzzz[i] = -2.0 * tr_xxzz_xzzz[i] * tbe_0 - 2.0 * tr_xxzz_xzzz[i] * tke_0 + 4.0 * tr_xxzz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_yyyy[i] = 12.0 * tr_xxzz_yy[i] - 2.0 * tr_xxzz_yyyy[i] * tbe_0 - 18.0 * tr_xxzz_yyyy[i] * tke_0 + 4.0 * tr_xxzz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xxyzz_yyy[i] * tbe_0 + 8.0 * tr_xxyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_yyyz[i] = 6.0 * tr_xxzz_yz[i] - 2.0 * tr_xxzz_yyyz[i] * tbe_0 - 14.0 * tr_xxzz_yyyz[i] * tke_0 + 4.0 * tr_xxzz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xxyzz_yyz[i] * tbe_0 + 8.0 * tr_xxyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_yyzz[i] = 2.0 * tr_xxzz_zz[i] - 2.0 * tr_xxzz_yyzz[i] * tbe_0 - 10.0 * tr_xxzz_yyzz[i] * tke_0 + 4.0 * tr_xxzz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_yzz[i] * tbe_0 + 8.0 * tr_xxyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_yzzz[i] = -2.0 * tr_xxzz_yzzz[i] * tbe_0 - 6.0 * tr_xxzz_yzzz[i] * tke_0 + 4.0 * tr_xxzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_zzz[i] * tbe_0 + 8.0 * tr_xxyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_zzzz[i] = -2.0 * tr_xxzz_zzzz[i] * tbe_0 - 2.0 * tr_xxzz_zzzz[i] * tke_0 + 4.0 * tr_xxzz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 765-780 components of targeted buffer : GG

    auto tr_0_0_yy_xyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 765);

    auto tr_0_0_yy_xyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 766);

    auto tr_0_0_yy_xyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 767);

    auto tr_0_0_yy_xyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 768);

    auto tr_0_0_yy_xyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 769);

    auto tr_0_0_yy_xyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 770);

    auto tr_0_0_yy_xyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 771);

    auto tr_0_0_yy_xyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 772);

    auto tr_0_0_yy_xyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 773);

    auto tr_0_0_yy_xyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 774);

    auto tr_0_0_yy_xyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 775);

    auto tr_0_0_yy_xyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 776);

    auto tr_0_0_yy_xyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 777);

    auto tr_0_0_yy_xyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 778);

    auto tr_0_0_yy_xyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 779);

    #pragma omp simd aligned(tr_0_0_yy_xyyy_xxxx, tr_0_0_yy_xyyy_xxxy, tr_0_0_yy_xyyy_xxxz, tr_0_0_yy_xyyy_xxyy, tr_0_0_yy_xyyy_xxyz, tr_0_0_yy_xyyy_xxzz, tr_0_0_yy_xyyy_xyyy, tr_0_0_yy_xyyy_xyyz, tr_0_0_yy_xyyy_xyzz, tr_0_0_yy_xyyy_xzzz, tr_0_0_yy_xyyy_yyyy, tr_0_0_yy_xyyy_yyyz, tr_0_0_yy_xyyy_yyzz, tr_0_0_yy_xyyy_yzzz, tr_0_0_yy_xyyy_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxxyy, tr_xyyy_xxxy, tr_xyyy_xxxyyy, tr_xyyy_xxxyyz, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyyyy, tr_xyyy_xxyyyz, tr_xyyy_xxyyzz, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyyyy, tr_xyyy_xyyyyz, tr_xyyy_xyyyzz, tr_xyyy_xyyz, tr_xyyy_xyyzzz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyyyy, tr_xyyy_yyyyyz, tr_xyyy_yyyyzz, tr_xyyy_yyyz, tr_xyyy_yyyzzz, tr_xyyy_yyzz, tr_xyyy_yyzzzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyy_xxx, tr_xyyyy_xxxxy, tr_xyyyy_xxxyy, tr_xyyyy_xxxyz, tr_xyyyy_xxy, tr_xyyyy_xxyyy, tr_xyyyy_xxyyz, tr_xyyyy_xxyzz, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyyyy, tr_xyyyy_xyyyz, tr_xyyyy_xyyzz, tr_xyyyy_xyz, tr_xyyyy_xyzzz, tr_xyyyy_xzz, tr_xyyyy_yyy, tr_xyyyy_yyyyy, tr_xyyyy_yyyyz, tr_xyyyy_yyyzz, tr_xyyyy_yyz, tr_xyyyy_yyzzz, tr_xyyyy_yzz, tr_xyyyy_yzzzz, tr_xyyyy_zzz, tr_xyyyyy_xxxx, tr_xyyyyy_xxxy, tr_xyyyyy_xxxz, tr_xyyyyy_xxyy, tr_xyyyyy_xxyz, tr_xyyyyy_xxzz, tr_xyyyyy_xyyy, tr_xyyyyy_xyyz, tr_xyyyyy_xyzz, tr_xyyyyy_xzzz, tr_xyyyyy_yyyy, tr_xyyyyy_yyyz, tr_xyyyyy_yyzz, tr_xyyyyy_yzzz, tr_xyyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyyy_xxxx[i] = 6.0 * tr_xy_xxxx[i] - 12.0 * tr_xyy_xxxxy[i] * tke_0 - 14.0 * tr_xyyy_xxxx[i] * tbe_0 - 2.0 * tr_xyyy_xxxx[i] * tke_0 + 4.0 * tr_xyyy_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xxxy[i] = 6.0 * tr_xy_xxxy[i] + 6.0 * tr_xyy_xxx[i] - 12.0 * tr_xyy_xxxyy[i] * tke_0 - 14.0 * tr_xyyy_xxxy[i] * tbe_0 - 6.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_xxx[i] * tbe_0 + 8.0 * tr_xyyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xxxz[i] = 6.0 * tr_xy_xxxz[i] - 12.0 * tr_xyy_xxxyz[i] * tke_0 - 14.0 * tr_xyyy_xxxz[i] * tbe_0 - 2.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xxyy[i] = 6.0 * tr_xy_xxyy[i] + 12.0 * tr_xyy_xxy[i] - 12.0 * tr_xyy_xxyyy[i] * tke_0 + 2.0 * tr_xyyy_xx[i] - 14.0 * tr_xyyy_xxyy[i] * tbe_0 - 10.0 * tr_xyyy_xxyy[i] * tke_0 + 4.0 * tr_xyyy_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_xxy[i] * tbe_0 + 8.0 * tr_xyyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xxyz[i] = 6.0 * tr_xy_xxyz[i] + 6.0 * tr_xyy_xxz[i] - 12.0 * tr_xyy_xxyyz[i] * tke_0 - 14.0 * tr_xyyy_xxyz[i] * tbe_0 - 6.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_xxz[i] * tbe_0 + 8.0 * tr_xyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xxzz[i] = 6.0 * tr_xy_xxzz[i] - 12.0 * tr_xyy_xxyzz[i] * tke_0 - 14.0 * tr_xyyy_xxzz[i] * tbe_0 - 2.0 * tr_xyyy_xxzz[i] * tke_0 + 4.0 * tr_xyyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xyyy[i] = 6.0 * tr_xy_xyyy[i] + 18.0 * tr_xyy_xyy[i] - 12.0 * tr_xyy_xyyyy[i] * tke_0 + 6.0 * tr_xyyy_xy[i] - 14.0 * tr_xyyy_xyyy[i] * tbe_0 - 14.0 * tr_xyyy_xyyy[i] * tke_0 + 4.0 * tr_xyyy_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyyyy_xyy[i] * tbe_0 + 8.0 * tr_xyyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xyyz[i] = 6.0 * tr_xy_xyyz[i] + 12.0 * tr_xyy_xyz[i] - 12.0 * tr_xyy_xyyyz[i] * tke_0 + 2.0 * tr_xyyy_xz[i] - 14.0 * tr_xyyy_xyyz[i] * tbe_0 - 10.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_xyz[i] * tbe_0 + 8.0 * tr_xyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xyzz[i] = 6.0 * tr_xy_xyzz[i] + 6.0 * tr_xyy_xzz[i] - 12.0 * tr_xyy_xyyzz[i] * tke_0 - 14.0 * tr_xyyy_xyzz[i] * tbe_0 - 6.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_xzz[i] * tbe_0 + 8.0 * tr_xyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xzzz[i] = 6.0 * tr_xy_xzzz[i] - 12.0 * tr_xyy_xyzzz[i] * tke_0 - 14.0 * tr_xyyy_xzzz[i] * tbe_0 - 2.0 * tr_xyyy_xzzz[i] * tke_0 + 4.0 * tr_xyyy_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_yyyy[i] = 6.0 * tr_xy_yyyy[i] + 24.0 * tr_xyy_yyy[i] - 12.0 * tr_xyy_yyyyy[i] * tke_0 + 12.0 * tr_xyyy_yy[i] - 14.0 * tr_xyyy_yyyy[i] * tbe_0 - 18.0 * tr_xyyy_yyyy[i] * tke_0 + 4.0 * tr_xyyy_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xyyyy_yyy[i] * tbe_0 + 8.0 * tr_xyyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_yyyz[i] = 6.0 * tr_xy_yyyz[i] + 18.0 * tr_xyy_yyz[i] - 12.0 * tr_xyy_yyyyz[i] * tke_0 + 6.0 * tr_xyyy_yz[i] - 14.0 * tr_xyyy_yyyz[i] * tbe_0 - 14.0 * tr_xyyy_yyyz[i] * tke_0 + 4.0 * tr_xyyy_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xyyyy_yyz[i] * tbe_0 + 8.0 * tr_xyyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_yyzz[i] = 6.0 * tr_xy_yyzz[i] + 12.0 * tr_xyy_yzz[i] - 12.0 * tr_xyy_yyyzz[i] * tke_0 + 2.0 * tr_xyyy_zz[i] - 14.0 * tr_xyyy_yyzz[i] * tbe_0 - 10.0 * tr_xyyy_yyzz[i] * tke_0 + 4.0 * tr_xyyy_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_yzz[i] * tbe_0 + 8.0 * tr_xyyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_yzzz[i] = 6.0 * tr_xy_yzzz[i] + 6.0 * tr_xyy_zzz[i] - 12.0 * tr_xyy_yyzzz[i] * tke_0 - 14.0 * tr_xyyy_yzzz[i] * tbe_0 - 6.0 * tr_xyyy_yzzz[i] * tke_0 + 4.0 * tr_xyyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_zzz[i] * tbe_0 + 8.0 * tr_xyyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_zzzz[i] = 6.0 * tr_xy_zzzz[i] - 12.0 * tr_xyy_yzzzz[i] * tke_0 - 14.0 * tr_xyyy_zzzz[i] * tbe_0 - 2.0 * tr_xyyy_zzzz[i] * tke_0 + 4.0 * tr_xyyy_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 780-795 components of targeted buffer : GG

    auto tr_0_0_yy_xyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 780);

    auto tr_0_0_yy_xyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 781);

    auto tr_0_0_yy_xyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 782);

    auto tr_0_0_yy_xyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 783);

    auto tr_0_0_yy_xyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 784);

    auto tr_0_0_yy_xyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 785);

    auto tr_0_0_yy_xyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 786);

    auto tr_0_0_yy_xyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 787);

    auto tr_0_0_yy_xyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 788);

    auto tr_0_0_yy_xyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 789);

    auto tr_0_0_yy_xyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 790);

    auto tr_0_0_yy_xyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 791);

    auto tr_0_0_yy_xyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 792);

    auto tr_0_0_yy_xyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 793);

    auto tr_0_0_yy_xyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 794);

    #pragma omp simd aligned(tr_0_0_yy_xyyz_xxxx, tr_0_0_yy_xyyz_xxxy, tr_0_0_yy_xyyz_xxxz, tr_0_0_yy_xyyz_xxyy, tr_0_0_yy_xyyz_xxyz, tr_0_0_yy_xyyz_xxzz, tr_0_0_yy_xyyz_xyyy, tr_0_0_yy_xyyz_xyyz, tr_0_0_yy_xyyz_xyzz, tr_0_0_yy_xyyz_xzzz, tr_0_0_yy_xyyz_yyyy, tr_0_0_yy_xyyz_yyyz, tr_0_0_yy_xyyz_yyzz, tr_0_0_yy_xyyz_yzzz, tr_0_0_yy_xyyz_zzzz, tr_xyyyyz_xxxx, tr_xyyyyz_xxxy, tr_xyyyyz_xxxz, tr_xyyyyz_xxyy, tr_xyyyyz_xxyz, tr_xyyyyz_xxzz, tr_xyyyyz_xyyy, tr_xyyyyz_xyyz, tr_xyyyyz_xyzz, tr_xyyyyz_xzzz, tr_xyyyyz_yyyy, tr_xyyyyz_yyyz, tr_xyyyyz_yyzz, tr_xyyyyz_yzzz, tr_xyyyyz_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxy, tr_xyyyz_xxxyy, tr_xyyyz_xxxyz, tr_xyyyz_xxy, tr_xyyyz_xxyyy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyyyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyyyy, tr_xyyyz_yyyyz, tr_xyyyz_yyyzz, tr_xyyyz_yyz, tr_xyyyz_yyzzz, tr_xyyyz_yzz, tr_xyyyz_yzzzz, tr_xyyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxxyy, tr_xyyz_xxxy, tr_xyyz_xxxyyy, tr_xyyz_xxxyyz, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyyyy, tr_xyyz_xxyyyz, tr_xyyz_xxyyzz, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyyyy, tr_xyyz_xyyyyz, tr_xyyz_xyyyzz, tr_xyyz_xyyz, tr_xyyz_xyyzzz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyyyy, tr_xyyz_yyyyyz, tr_xyyz_yyyyzz, tr_xyyz_yyyz, tr_xyyz_yyyzzz, tr_xyyz_yyzz, tr_xyyz_yyzzzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyyz_xxxx[i] = 2.0 * tr_xz_xxxx[i] - 8.0 * tr_xyz_xxxxy[i] * tke_0 - 10.0 * tr_xyyz_xxxx[i] * tbe_0 - 2.0 * tr_xyyz_xxxx[i] * tke_0 + 4.0 * tr_xyyz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xxxy[i] = 2.0 * tr_xz_xxxy[i] + 4.0 * tr_xyz_xxx[i] - 8.0 * tr_xyz_xxxyy[i] * tke_0 - 10.0 * tr_xyyz_xxxy[i] * tbe_0 - 6.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xxx[i] * tbe_0 + 8.0 * tr_xyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xxxz[i] = 2.0 * tr_xz_xxxz[i] - 8.0 * tr_xyz_xxxyz[i] * tke_0 - 10.0 * tr_xyyz_xxxz[i] * tbe_0 - 2.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xxyy[i] = 2.0 * tr_xz_xxyy[i] + 8.0 * tr_xyz_xxy[i] - 8.0 * tr_xyz_xxyyy[i] * tke_0 + 2.0 * tr_xyyz_xx[i] - 10.0 * tr_xyyz_xxyy[i] * tbe_0 - 10.0 * tr_xyyz_xxyy[i] * tke_0 + 4.0 * tr_xyyz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xxy[i] * tbe_0 + 8.0 * tr_xyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xxyz[i] = 2.0 * tr_xz_xxyz[i] + 4.0 * tr_xyz_xxz[i] - 8.0 * tr_xyz_xxyyz[i] * tke_0 - 10.0 * tr_xyyz_xxyz[i] * tbe_0 - 6.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xxz[i] * tbe_0 + 8.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xxzz[i] = 2.0 * tr_xz_xxzz[i] - 8.0 * tr_xyz_xxyzz[i] * tke_0 - 10.0 * tr_xyyz_xxzz[i] * tbe_0 - 2.0 * tr_xyyz_xxzz[i] * tke_0 + 4.0 * tr_xyyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xyyy[i] = 2.0 * tr_xz_xyyy[i] + 12.0 * tr_xyz_xyy[i] - 8.0 * tr_xyz_xyyyy[i] * tke_0 + 6.0 * tr_xyyz_xy[i] - 10.0 * tr_xyyz_xyyy[i] * tbe_0 - 14.0 * tr_xyyz_xyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyyyz_xyy[i] * tbe_0 + 8.0 * tr_xyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xyyz[i] = 2.0 * tr_xz_xyyz[i] + 8.0 * tr_xyz_xyz[i] - 8.0 * tr_xyz_xyyyz[i] * tke_0 + 2.0 * tr_xyyz_xz[i] - 10.0 * tr_xyyz_xyyz[i] * tbe_0 - 10.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xyz[i] * tbe_0 + 8.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xyzz[i] = 2.0 * tr_xz_xyzz[i] + 4.0 * tr_xyz_xzz[i] - 8.0 * tr_xyz_xyyzz[i] * tke_0 - 10.0 * tr_xyyz_xyzz[i] * tbe_0 - 6.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xzz[i] * tbe_0 + 8.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xzzz[i] = 2.0 * tr_xz_xzzz[i] - 8.0 * tr_xyz_xyzzz[i] * tke_0 - 10.0 * tr_xyyz_xzzz[i] * tbe_0 - 2.0 * tr_xyyz_xzzz[i] * tke_0 + 4.0 * tr_xyyz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_yyyy[i] = 2.0 * tr_xz_yyyy[i] + 16.0 * tr_xyz_yyy[i] - 8.0 * tr_xyz_yyyyy[i] * tke_0 + 12.0 * tr_xyyz_yy[i] - 10.0 * tr_xyyz_yyyy[i] * tbe_0 - 18.0 * tr_xyyz_yyyy[i] * tke_0 + 4.0 * tr_xyyz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xyyyz_yyy[i] * tbe_0 + 8.0 * tr_xyyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_yyyz[i] = 2.0 * tr_xz_yyyz[i] + 12.0 * tr_xyz_yyz[i] - 8.0 * tr_xyz_yyyyz[i] * tke_0 + 6.0 * tr_xyyz_yz[i] - 10.0 * tr_xyyz_yyyz[i] * tbe_0 - 14.0 * tr_xyyz_yyyz[i] * tke_0 + 4.0 * tr_xyyz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xyyyz_yyz[i] * tbe_0 + 8.0 * tr_xyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_yyzz[i] = 2.0 * tr_xz_yyzz[i] + 8.0 * tr_xyz_yzz[i] - 8.0 * tr_xyz_yyyzz[i] * tke_0 + 2.0 * tr_xyyz_zz[i] - 10.0 * tr_xyyz_yyzz[i] * tbe_0 - 10.0 * tr_xyyz_yyzz[i] * tke_0 + 4.0 * tr_xyyz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_yzz[i] * tbe_0 + 8.0 * tr_xyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_yzzz[i] = 2.0 * tr_xz_yzzz[i] + 4.0 * tr_xyz_zzz[i] - 8.0 * tr_xyz_yyzzz[i] * tke_0 - 10.0 * tr_xyyz_yzzz[i] * tbe_0 - 6.0 * tr_xyyz_yzzz[i] * tke_0 + 4.0 * tr_xyyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_zzz[i] * tbe_0 + 8.0 * tr_xyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_zzzz[i] = 2.0 * tr_xz_zzzz[i] - 8.0 * tr_xyz_yzzzz[i] * tke_0 - 10.0 * tr_xyyz_zzzz[i] * tbe_0 - 2.0 * tr_xyyz_zzzz[i] * tke_0 + 4.0 * tr_xyyz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 795-810 components of targeted buffer : GG

    auto tr_0_0_yy_xyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 795);

    auto tr_0_0_yy_xyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 796);

    auto tr_0_0_yy_xyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 797);

    auto tr_0_0_yy_xyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 798);

    auto tr_0_0_yy_xyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 799);

    auto tr_0_0_yy_xyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 800);

    auto tr_0_0_yy_xyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 801);

    auto tr_0_0_yy_xyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 802);

    auto tr_0_0_yy_xyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 803);

    auto tr_0_0_yy_xyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 804);

    auto tr_0_0_yy_xyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 805);

    auto tr_0_0_yy_xyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 806);

    auto tr_0_0_yy_xyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 807);

    auto tr_0_0_yy_xyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 808);

    auto tr_0_0_yy_xyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 809);

    #pragma omp simd aligned(tr_0_0_yy_xyzz_xxxx, tr_0_0_yy_xyzz_xxxy, tr_0_0_yy_xyzz_xxxz, tr_0_0_yy_xyzz_xxyy, tr_0_0_yy_xyzz_xxyz, tr_0_0_yy_xyzz_xxzz, tr_0_0_yy_xyzz_xyyy, tr_0_0_yy_xyzz_xyyz, tr_0_0_yy_xyzz_xyzz, tr_0_0_yy_xyzz_xzzz, tr_0_0_yy_xyzz_yyyy, tr_0_0_yy_xyzz_yyyz, tr_0_0_yy_xyzz_yyzz, tr_0_0_yy_xyzz_yzzz, tr_0_0_yy_xyzz_zzzz, tr_xyyyzz_xxxx, tr_xyyyzz_xxxy, tr_xyyyzz_xxxz, tr_xyyyzz_xxyy, tr_xyyyzz_xxyz, tr_xyyyzz_xxzz, tr_xyyyzz_xyyy, tr_xyyyzz_xyyz, tr_xyyyzz_xyzz, tr_xyyyzz_xzzz, tr_xyyyzz_yyyy, tr_xyyyzz_yyyz, tr_xyyyzz_yyzz, tr_xyyyzz_yzzz, tr_xyyyzz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxy, tr_xyyzz_xxxyy, tr_xyyzz_xxxyz, tr_xyyzz_xxy, tr_xyyzz_xxyyy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyyyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyyyy, tr_xyyzz_yyyyz, tr_xyyzz_yyyzz, tr_xyyzz_yyz, tr_xyyzz_yyzzz, tr_xyyzz_yzz, tr_xyyzz_yzzzz, tr_xyyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxxyy, tr_xyzz_xxxy, tr_xyzz_xxxyyy, tr_xyzz_xxxyyz, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyyyy, tr_xyzz_xxyyyz, tr_xyzz_xxyyzz, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyyyy, tr_xyzz_xyyyyz, tr_xyzz_xyyyzz, tr_xyzz_xyyz, tr_xyzz_xyyzzz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyyyy, tr_xyzz_yyyyyz, tr_xyzz_yyyyzz, tr_xyzz_yyyz, tr_xyzz_yyyzzz, tr_xyzz_yyzz, tr_xyzz_yyzzzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyzz_xxxx[i] = -4.0 * tr_xzz_xxxxy[i] * tke_0 - 6.0 * tr_xyzz_xxxx[i] * tbe_0 - 2.0 * tr_xyzz_xxxx[i] * tke_0 + 4.0 * tr_xyzz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xxxy[i] = 2.0 * tr_xzz_xxx[i] - 4.0 * tr_xzz_xxxyy[i] * tke_0 - 6.0 * tr_xyzz_xxxy[i] * tbe_0 - 6.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xxx[i] * tbe_0 + 8.0 * tr_xyyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xxxz[i] = -4.0 * tr_xzz_xxxyz[i] * tke_0 - 6.0 * tr_xyzz_xxxz[i] * tbe_0 - 2.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xxyy[i] = 4.0 * tr_xzz_xxy[i] - 4.0 * tr_xzz_xxyyy[i] * tke_0 + 2.0 * tr_xyzz_xx[i] - 6.0 * tr_xyzz_xxyy[i] * tbe_0 - 10.0 * tr_xyzz_xxyy[i] * tke_0 + 4.0 * tr_xyzz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xxy[i] * tbe_0 + 8.0 * tr_xyyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xxyz[i] = 2.0 * tr_xzz_xxz[i] - 4.0 * tr_xzz_xxyyz[i] * tke_0 - 6.0 * tr_xyzz_xxyz[i] * tbe_0 - 6.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xxz[i] * tbe_0 + 8.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xxzz[i] = -4.0 * tr_xzz_xxyzz[i] * tke_0 - 6.0 * tr_xyzz_xxzz[i] * tbe_0 - 2.0 * tr_xyzz_xxzz[i] * tke_0 + 4.0 * tr_xyzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xyyy[i] = 6.0 * tr_xzz_xyy[i] - 4.0 * tr_xzz_xyyyy[i] * tke_0 + 6.0 * tr_xyzz_xy[i] - 6.0 * tr_xyzz_xyyy[i] * tbe_0 - 14.0 * tr_xyzz_xyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyyzz_xyy[i] * tbe_0 + 8.0 * tr_xyyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xyyz[i] = 4.0 * tr_xzz_xyz[i] - 4.0 * tr_xzz_xyyyz[i] * tke_0 + 2.0 * tr_xyzz_xz[i] - 6.0 * tr_xyzz_xyyz[i] * tbe_0 - 10.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xyz[i] * tbe_0 + 8.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xyzz[i] = 2.0 * tr_xzz_xzz[i] - 4.0 * tr_xzz_xyyzz[i] * tke_0 - 6.0 * tr_xyzz_xyzz[i] * tbe_0 - 6.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xzz[i] * tbe_0 + 8.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xzzz[i] = -4.0 * tr_xzz_xyzzz[i] * tke_0 - 6.0 * tr_xyzz_xzzz[i] * tbe_0 - 2.0 * tr_xyzz_xzzz[i] * tke_0 + 4.0 * tr_xyzz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_yyyy[i] = 8.0 * tr_xzz_yyy[i] - 4.0 * tr_xzz_yyyyy[i] * tke_0 + 12.0 * tr_xyzz_yy[i] - 6.0 * tr_xyzz_yyyy[i] * tbe_0 - 18.0 * tr_xyzz_yyyy[i] * tke_0 + 4.0 * tr_xyzz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xyyzz_yyy[i] * tbe_0 + 8.0 * tr_xyyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_yyyz[i] = 6.0 * tr_xzz_yyz[i] - 4.0 * tr_xzz_yyyyz[i] * tke_0 + 6.0 * tr_xyzz_yz[i] - 6.0 * tr_xyzz_yyyz[i] * tbe_0 - 14.0 * tr_xyzz_yyyz[i] * tke_0 + 4.0 * tr_xyzz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xyyzz_yyz[i] * tbe_0 + 8.0 * tr_xyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_yyzz[i] = 4.0 * tr_xzz_yzz[i] - 4.0 * tr_xzz_yyyzz[i] * tke_0 + 2.0 * tr_xyzz_zz[i] - 6.0 * tr_xyzz_yyzz[i] * tbe_0 - 10.0 * tr_xyzz_yyzz[i] * tke_0 + 4.0 * tr_xyzz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_yzz[i] * tbe_0 + 8.0 * tr_xyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_yzzz[i] = 2.0 * tr_xzz_zzz[i] - 4.0 * tr_xzz_yyzzz[i] * tke_0 - 6.0 * tr_xyzz_yzzz[i] * tbe_0 - 6.0 * tr_xyzz_yzzz[i] * tke_0 + 4.0 * tr_xyzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_zzz[i] * tbe_0 + 8.0 * tr_xyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_zzzz[i] = -4.0 * tr_xzz_yzzzz[i] * tke_0 - 6.0 * tr_xyzz_zzzz[i] * tbe_0 - 2.0 * tr_xyzz_zzzz[i] * tke_0 + 4.0 * tr_xyzz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 810-825 components of targeted buffer : GG

    auto tr_0_0_yy_xzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 810);

    auto tr_0_0_yy_xzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 811);

    auto tr_0_0_yy_xzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 812);

    auto tr_0_0_yy_xzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 813);

    auto tr_0_0_yy_xzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 814);

    auto tr_0_0_yy_xzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 815);

    auto tr_0_0_yy_xzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 816);

    auto tr_0_0_yy_xzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 817);

    auto tr_0_0_yy_xzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 818);

    auto tr_0_0_yy_xzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 819);

    auto tr_0_0_yy_xzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 820);

    auto tr_0_0_yy_xzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 821);

    auto tr_0_0_yy_xzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 822);

    auto tr_0_0_yy_xzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 823);

    auto tr_0_0_yy_xzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 824);

    #pragma omp simd aligned(tr_0_0_yy_xzzz_xxxx, tr_0_0_yy_xzzz_xxxy, tr_0_0_yy_xzzz_xxxz, tr_0_0_yy_xzzz_xxyy, tr_0_0_yy_xzzz_xxyz, tr_0_0_yy_xzzz_xxzz, tr_0_0_yy_xzzz_xyyy, tr_0_0_yy_xzzz_xyyz, tr_0_0_yy_xzzz_xyzz, tr_0_0_yy_xzzz_xzzz, tr_0_0_yy_xzzz_yyyy, tr_0_0_yy_xzzz_yyyz, tr_0_0_yy_xzzz_yyzz, tr_0_0_yy_xzzz_yzzz, tr_0_0_yy_xzzz_zzzz, tr_xyyzzz_xxxx, tr_xyyzzz_xxxy, tr_xyyzzz_xxxz, tr_xyyzzz_xxyy, tr_xyyzzz_xxyz, tr_xyyzzz_xxzz, tr_xyyzzz_xyyy, tr_xyyzzz_xyyz, tr_xyyzzz_xyzz, tr_xyyzzz_xzzz, tr_xyyzzz_yyyy, tr_xyyzzz_yyyz, tr_xyyzzz_yyzz, tr_xyyzzz_yzzz, tr_xyyzzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxy, tr_xyzzz_xxxyy, tr_xyzzz_xxxyz, tr_xyzzz_xxy, tr_xyzzz_xxyyy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyyyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyyyy, tr_xyzzz_yyyyz, tr_xyzzz_yyyzz, tr_xyzzz_yyz, tr_xyzzz_yyzzz, tr_xyzzz_yzz, tr_xyzzz_yzzzz, tr_xyzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxxyy, tr_xzzz_xxxy, tr_xzzz_xxxyyy, tr_xzzz_xxxyyz, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyyyy, tr_xzzz_xxyyyz, tr_xzzz_xxyyzz, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyyyy, tr_xzzz_xyyyyz, tr_xzzz_xyyyzz, tr_xzzz_xyyz, tr_xzzz_xyyzzz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyyyy, tr_xzzz_yyyyyz, tr_xzzz_yyyyzz, tr_xzzz_yyyz, tr_xzzz_yyyzzz, tr_xzzz_yyzz, tr_xzzz_yyzzzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xzzz_xxxx[i] = -2.0 * tr_xzzz_xxxx[i] * tbe_0 - 2.0 * tr_xzzz_xxxx[i] * tke_0 + 4.0 * tr_xzzz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xxxy[i] = -2.0 * tr_xzzz_xxxy[i] * tbe_0 - 6.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xxx[i] * tbe_0 + 8.0 * tr_xyzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xxxz[i] = -2.0 * tr_xzzz_xxxz[i] * tbe_0 - 2.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xxyy[i] = 2.0 * tr_xzzz_xx[i] - 2.0 * tr_xzzz_xxyy[i] * tbe_0 - 10.0 * tr_xzzz_xxyy[i] * tke_0 + 4.0 * tr_xzzz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xxy[i] * tbe_0 + 8.0 * tr_xyzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xxyz[i] = -2.0 * tr_xzzz_xxyz[i] * tbe_0 - 6.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xxz[i] * tbe_0 + 8.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xxzz[i] = -2.0 * tr_xzzz_xxzz[i] * tbe_0 - 2.0 * tr_xzzz_xxzz[i] * tke_0 + 4.0 * tr_xzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xyyy[i] = 6.0 * tr_xzzz_xy[i] - 2.0 * tr_xzzz_xyyy[i] * tbe_0 - 14.0 * tr_xzzz_xyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_xyzzz_xyy[i] * tbe_0 + 8.0 * tr_xyzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xyyz[i] = 2.0 * tr_xzzz_xz[i] - 2.0 * tr_xzzz_xyyz[i] * tbe_0 - 10.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xyz[i] * tbe_0 + 8.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xyzz[i] = -2.0 * tr_xzzz_xyzz[i] * tbe_0 - 6.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xzz[i] * tbe_0 + 8.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xzzz[i] = -2.0 * tr_xzzz_xzzz[i] * tbe_0 - 2.0 * tr_xzzz_xzzz[i] * tke_0 + 4.0 * tr_xzzz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_yyyy[i] = 12.0 * tr_xzzz_yy[i] - 2.0 * tr_xzzz_yyyy[i] * tbe_0 - 18.0 * tr_xzzz_yyyy[i] * tke_0 + 4.0 * tr_xzzz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_xyzzz_yyy[i] * tbe_0 + 8.0 * tr_xyzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_yyyz[i] = 6.0 * tr_xzzz_yz[i] - 2.0 * tr_xzzz_yyyz[i] * tbe_0 - 14.0 * tr_xzzz_yyyz[i] * tke_0 + 4.0 * tr_xzzz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_xyzzz_yyz[i] * tbe_0 + 8.0 * tr_xyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_yyzz[i] = 2.0 * tr_xzzz_zz[i] - 2.0 * tr_xzzz_yyzz[i] * tbe_0 - 10.0 * tr_xzzz_yyzz[i] * tke_0 + 4.0 * tr_xzzz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_yzz[i] * tbe_0 + 8.0 * tr_xyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_yzzz[i] = -2.0 * tr_xzzz_yzzz[i] * tbe_0 - 6.0 * tr_xzzz_yzzz[i] * tke_0 + 4.0 * tr_xzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_zzz[i] * tbe_0 + 8.0 * tr_xyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_zzzz[i] = -2.0 * tr_xzzz_zzzz[i] * tbe_0 - 2.0 * tr_xzzz_zzzz[i] * tke_0 + 4.0 * tr_xzzz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 825-840 components of targeted buffer : GG

    auto tr_0_0_yy_yyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 825);

    auto tr_0_0_yy_yyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 826);

    auto tr_0_0_yy_yyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 827);

    auto tr_0_0_yy_yyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 828);

    auto tr_0_0_yy_yyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 829);

    auto tr_0_0_yy_yyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 830);

    auto tr_0_0_yy_yyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 831);

    auto tr_0_0_yy_yyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 832);

    auto tr_0_0_yy_yyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 833);

    auto tr_0_0_yy_yyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 834);

    auto tr_0_0_yy_yyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 835);

    auto tr_0_0_yy_yyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 836);

    auto tr_0_0_yy_yyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 837);

    auto tr_0_0_yy_yyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 838);

    auto tr_0_0_yy_yyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 839);

    #pragma omp simd aligned(tr_0_0_yy_yyyy_xxxx, tr_0_0_yy_yyyy_xxxy, tr_0_0_yy_yyyy_xxxz, tr_0_0_yy_yyyy_xxyy, tr_0_0_yy_yyyy_xxyz, tr_0_0_yy_yyyy_xxzz, tr_0_0_yy_yyyy_xyyy, tr_0_0_yy_yyyy_xyyz, tr_0_0_yy_yyyy_xyzz, tr_0_0_yy_yyyy_xzzz, tr_0_0_yy_yyyy_yyyy, tr_0_0_yy_yyyy_yyyz, tr_0_0_yy_yyyy_yyzz, tr_0_0_yy_yyyy_yzzz, tr_0_0_yy_yyyy_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyy_xxx, tr_yyy_xxxxy, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyyyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxxyy, tr_yyyy_xxxy, tr_yyyy_xxxyyy, tr_yyyy_xxxyyz, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyyyy, tr_yyyy_xxyyyz, tr_yyyy_xxyyzz, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyyyy, tr_yyyy_xyyyyz, tr_yyyy_xyyyzz, tr_yyyy_xyyz, tr_yyyy_xyyzzz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyyyy, tr_yyyy_yyyyyz, tr_yyyy_yyyyzz, tr_yyyy_yyyz, tr_yyyy_yyyzzz, tr_yyyy_yyzz, tr_yyyy_yyzzzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyyy_xxx, tr_yyyyy_xxxxy, tr_yyyyy_xxxyy, tr_yyyyy_xxxyz, tr_yyyyy_xxy, tr_yyyyy_xxyyy, tr_yyyyy_xxyyz, tr_yyyyy_xxyzz, tr_yyyyy_xxz, tr_yyyyy_xyy, tr_yyyyy_xyyyy, tr_yyyyy_xyyyz, tr_yyyyy_xyyzz, tr_yyyyy_xyz, tr_yyyyy_xyzzz, tr_yyyyy_xzz, tr_yyyyy_yyy, tr_yyyyy_yyyyy, tr_yyyyy_yyyyz, tr_yyyyy_yyyzz, tr_yyyyy_yyz, tr_yyyyy_yyzzz, tr_yyyyy_yzz, tr_yyyyy_yzzzz, tr_yyyyy_zzz, tr_yyyyyy_xxxx, tr_yyyyyy_xxxy, tr_yyyyyy_xxxz, tr_yyyyyy_xxyy, tr_yyyyyy_xxyz, tr_yyyyyy_xxzz, tr_yyyyyy_xyyy, tr_yyyyyy_xyyz, tr_yyyyyy_xyzz, tr_yyyyyy_xzzz, tr_yyyyyy_yyyy, tr_yyyyyy_yyyz, tr_yyyyyy_yyzz, tr_yyyyyy_yzzz, tr_yyyyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyyy_xxxx[i] = 12.0 * tr_yy_xxxx[i] - 16.0 * tr_yyy_xxxxy[i] * tke_0 - 18.0 * tr_yyyy_xxxx[i] * tbe_0 - 2.0 * tr_yyyy_xxxx[i] * tke_0 + 4.0 * tr_yyyy_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xxxy[i] = 12.0 * tr_yy_xxxy[i] + 8.0 * tr_yyy_xxx[i] - 16.0 * tr_yyy_xxxyy[i] * tke_0 - 18.0 * tr_yyyy_xxxy[i] * tbe_0 - 6.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_xxx[i] * tbe_0 + 8.0 * tr_yyyyy_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xxxz[i] = 12.0 * tr_yy_xxxz[i] - 16.0 * tr_yyy_xxxyz[i] * tke_0 - 18.0 * tr_yyyy_xxxz[i] * tbe_0 - 2.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xxyy[i] = 12.0 * tr_yy_xxyy[i] + 16.0 * tr_yyy_xxy[i] - 16.0 * tr_yyy_xxyyy[i] * tke_0 + 2.0 * tr_yyyy_xx[i] - 18.0 * tr_yyyy_xxyy[i] * tbe_0 - 10.0 * tr_yyyy_xxyy[i] * tke_0 + 4.0 * tr_yyyy_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyyy_xxy[i] * tbe_0 + 8.0 * tr_yyyyy_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xxyz[i] = 12.0 * tr_yy_xxyz[i] + 8.0 * tr_yyy_xxz[i] - 16.0 * tr_yyy_xxyyz[i] * tke_0 - 18.0 * tr_yyyy_xxyz[i] * tbe_0 - 6.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_xxz[i] * tbe_0 + 8.0 * tr_yyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xxzz[i] = 12.0 * tr_yy_xxzz[i] - 16.0 * tr_yyy_xxyzz[i] * tke_0 - 18.0 * tr_yyyy_xxzz[i] * tbe_0 - 2.0 * tr_yyyy_xxzz[i] * tke_0 + 4.0 * tr_yyyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xyyy[i] = 12.0 * tr_yy_xyyy[i] + 24.0 * tr_yyy_xyy[i] - 16.0 * tr_yyy_xyyyy[i] * tke_0 + 6.0 * tr_yyyy_xy[i] - 18.0 * tr_yyyy_xyyy[i] * tbe_0 - 14.0 * tr_yyyy_xyyy[i] * tke_0 + 4.0 * tr_yyyy_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyyyy_xyy[i] * tbe_0 + 8.0 * tr_yyyyy_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xyyz[i] = 12.0 * tr_yy_xyyz[i] + 16.0 * tr_yyy_xyz[i] - 16.0 * tr_yyy_xyyyz[i] * tke_0 + 2.0 * tr_yyyy_xz[i] - 18.0 * tr_yyyy_xyyz[i] * tbe_0 - 10.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyy_xyz[i] * tbe_0 + 8.0 * tr_yyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xyzz[i] = 12.0 * tr_yy_xyzz[i] + 8.0 * tr_yyy_xzz[i] - 16.0 * tr_yyy_xyyzz[i] * tke_0 - 18.0 * tr_yyyy_xyzz[i] * tbe_0 - 6.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_xzz[i] * tbe_0 + 8.0 * tr_yyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xzzz[i] = 12.0 * tr_yy_xzzz[i] - 16.0 * tr_yyy_xyzzz[i] * tke_0 - 18.0 * tr_yyyy_xzzz[i] * tbe_0 - 2.0 * tr_yyyy_xzzz[i] * tke_0 + 4.0 * tr_yyyy_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_yyyy[i] = 12.0 * tr_yy_yyyy[i] + 32.0 * tr_yyy_yyy[i] - 16.0 * tr_yyy_yyyyy[i] * tke_0 + 12.0 * tr_yyyy_yy[i] - 18.0 * tr_yyyy_yyyy[i] * tbe_0 - 18.0 * tr_yyyy_yyyy[i] * tke_0 + 4.0 * tr_yyyy_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yyyyy_yyy[i] * tbe_0 + 8.0 * tr_yyyyy_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_yyyz[i] = 12.0 * tr_yy_yyyz[i] + 24.0 * tr_yyy_yyz[i] - 16.0 * tr_yyy_yyyyz[i] * tke_0 + 6.0 * tr_yyyy_yz[i] - 18.0 * tr_yyyy_yyyz[i] * tbe_0 - 14.0 * tr_yyyy_yyyz[i] * tke_0 + 4.0 * tr_yyyy_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yyyyy_yyz[i] * tbe_0 + 8.0 * tr_yyyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_yyzz[i] = 12.0 * tr_yy_yyzz[i] + 16.0 * tr_yyy_yzz[i] - 16.0 * tr_yyy_yyyzz[i] * tke_0 + 2.0 * tr_yyyy_zz[i] - 18.0 * tr_yyyy_yyzz[i] * tbe_0 - 10.0 * tr_yyyy_yyzz[i] * tke_0 + 4.0 * tr_yyyy_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyy_yzz[i] * tbe_0 + 8.0 * tr_yyyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_yzzz[i] = 12.0 * tr_yy_yzzz[i] + 8.0 * tr_yyy_zzz[i] - 16.0 * tr_yyy_yyzzz[i] * tke_0 - 18.0 * tr_yyyy_yzzz[i] * tbe_0 - 6.0 * tr_yyyy_yzzz[i] * tke_0 + 4.0 * tr_yyyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_zzz[i] * tbe_0 + 8.0 * tr_yyyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_zzzz[i] = 12.0 * tr_yy_zzzz[i] - 16.0 * tr_yyy_yzzzz[i] * tke_0 - 18.0 * tr_yyyy_zzzz[i] * tbe_0 - 2.0 * tr_yyyy_zzzz[i] * tke_0 + 4.0 * tr_yyyy_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 840-855 components of targeted buffer : GG

    auto tr_0_0_yy_yyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 840);

    auto tr_0_0_yy_yyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 841);

    auto tr_0_0_yy_yyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 842);

    auto tr_0_0_yy_yyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 843);

    auto tr_0_0_yy_yyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 844);

    auto tr_0_0_yy_yyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 845);

    auto tr_0_0_yy_yyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 846);

    auto tr_0_0_yy_yyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 847);

    auto tr_0_0_yy_yyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 848);

    auto tr_0_0_yy_yyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 849);

    auto tr_0_0_yy_yyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 850);

    auto tr_0_0_yy_yyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 851);

    auto tr_0_0_yy_yyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 852);

    auto tr_0_0_yy_yyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 853);

    auto tr_0_0_yy_yyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 854);

    #pragma omp simd aligned(tr_0_0_yy_yyyz_xxxx, tr_0_0_yy_yyyz_xxxy, tr_0_0_yy_yyyz_xxxz, tr_0_0_yy_yyyz_xxyy, tr_0_0_yy_yyyz_xxyz, tr_0_0_yy_yyyz_xxzz, tr_0_0_yy_yyyz_xyyy, tr_0_0_yy_yyyz_xyyz, tr_0_0_yy_yyyz_xyzz, tr_0_0_yy_yyyz_xzzz, tr_0_0_yy_yyyz_yyyy, tr_0_0_yy_yyyz_yyyz, tr_0_0_yy_yyyz_yyzz, tr_0_0_yy_yyyz_yzzz, tr_0_0_yy_yyyz_zzzz, tr_yyyyyz_xxxx, tr_yyyyyz_xxxy, tr_yyyyyz_xxxz, tr_yyyyyz_xxyy, tr_yyyyyz_xxyz, tr_yyyyyz_xxzz, tr_yyyyyz_xyyy, tr_yyyyyz_xyyz, tr_yyyyyz_xyzz, tr_yyyyyz_xzzz, tr_yyyyyz_yyyy, tr_yyyyyz_yyyz, tr_yyyyyz_yyzz, tr_yyyyyz_yzzz, tr_yyyyyz_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxxxy, tr_yyyyz_xxxyy, tr_yyyyz_xxxyz, tr_yyyyz_xxy, tr_yyyyz_xxyyy, tr_yyyyz_xxyyz, tr_yyyyz_xxyzz, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyyyy, tr_yyyyz_xyyyz, tr_yyyyz_xyyzz, tr_yyyyz_xyz, tr_yyyyz_xyzzz, tr_yyyyz_xzz, tr_yyyyz_yyy, tr_yyyyz_yyyyy, tr_yyyyz_yyyyz, tr_yyyyz_yyyzz, tr_yyyyz_yyz, tr_yyyyz_yyzzz, tr_yyyyz_yzz, tr_yyyyz_yzzzz, tr_yyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxxyy, tr_yyyz_xxxy, tr_yyyz_xxxyyy, tr_yyyz_xxxyyz, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyyyy, tr_yyyz_xxyyyz, tr_yyyz_xxyyzz, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyyyy, tr_yyyz_xyyyyz, tr_yyyz_xyyyzz, tr_yyyz_xyyz, tr_yyyz_xyyzzz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyyyy, tr_yyyz_yyyyyz, tr_yyyz_yyyyzz, tr_yyyz_yyyz, tr_yyyz_yyyzzz, tr_yyyz_yyzz, tr_yyyz_yyzzzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyyz_xxxx[i] = 6.0 * tr_yz_xxxx[i] - 12.0 * tr_yyz_xxxxy[i] * tke_0 - 14.0 * tr_yyyz_xxxx[i] * tbe_0 - 2.0 * tr_yyyz_xxxx[i] * tke_0 + 4.0 * tr_yyyz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xxxy[i] = 6.0 * tr_yz_xxxy[i] + 6.0 * tr_yyz_xxx[i] - 12.0 * tr_yyz_xxxyy[i] * tke_0 - 14.0 * tr_yyyz_xxxy[i] * tbe_0 - 6.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xxx[i] * tbe_0 + 8.0 * tr_yyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xxxz[i] = 6.0 * tr_yz_xxxz[i] - 12.0 * tr_yyz_xxxyz[i] * tke_0 - 14.0 * tr_yyyz_xxxz[i] * tbe_0 - 2.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xxyy[i] = 6.0 * tr_yz_xxyy[i] + 12.0 * tr_yyz_xxy[i] - 12.0 * tr_yyz_xxyyy[i] * tke_0 + 2.0 * tr_yyyz_xx[i] - 14.0 * tr_yyyz_xxyy[i] * tbe_0 - 10.0 * tr_yyyz_xxyy[i] * tke_0 + 4.0 * tr_yyyz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_xxy[i] * tbe_0 + 8.0 * tr_yyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xxyz[i] = 6.0 * tr_yz_xxyz[i] + 6.0 * tr_yyz_xxz[i] - 12.0 * tr_yyz_xxyyz[i] * tke_0 - 14.0 * tr_yyyz_xxyz[i] * tbe_0 - 6.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xxz[i] * tbe_0 + 8.0 * tr_yyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xxzz[i] = 6.0 * tr_yz_xxzz[i] - 12.0 * tr_yyz_xxyzz[i] * tke_0 - 14.0 * tr_yyyz_xxzz[i] * tbe_0 - 2.0 * tr_yyyz_xxzz[i] * tke_0 + 4.0 * tr_yyyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xyyy[i] = 6.0 * tr_yz_xyyy[i] + 18.0 * tr_yyz_xyy[i] - 12.0 * tr_yyz_xyyyy[i] * tke_0 + 6.0 * tr_yyyz_xy[i] - 14.0 * tr_yyyz_xyyy[i] * tbe_0 - 14.0 * tr_yyyz_xyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyyyz_xyy[i] * tbe_0 + 8.0 * tr_yyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xyyz[i] = 6.0 * tr_yz_xyyz[i] + 12.0 * tr_yyz_xyz[i] - 12.0 * tr_yyz_xyyyz[i] * tke_0 + 2.0 * tr_yyyz_xz[i] - 14.0 * tr_yyyz_xyyz[i] * tbe_0 - 10.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_xyz[i] * tbe_0 + 8.0 * tr_yyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xyzz[i] = 6.0 * tr_yz_xyzz[i] + 6.0 * tr_yyz_xzz[i] - 12.0 * tr_yyz_xyyzz[i] * tke_0 - 14.0 * tr_yyyz_xyzz[i] * tbe_0 - 6.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xzz[i] * tbe_0 + 8.0 * tr_yyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xzzz[i] = 6.0 * tr_yz_xzzz[i] - 12.0 * tr_yyz_xyzzz[i] * tke_0 - 14.0 * tr_yyyz_xzzz[i] * tbe_0 - 2.0 * tr_yyyz_xzzz[i] * tke_0 + 4.0 * tr_yyyz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_yyyy[i] = 6.0 * tr_yz_yyyy[i] + 24.0 * tr_yyz_yyy[i] - 12.0 * tr_yyz_yyyyy[i] * tke_0 + 12.0 * tr_yyyz_yy[i] - 14.0 * tr_yyyz_yyyy[i] * tbe_0 - 18.0 * tr_yyyz_yyyy[i] * tke_0 + 4.0 * tr_yyyz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yyyyz_yyy[i] * tbe_0 + 8.0 * tr_yyyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_yyyz[i] = 6.0 * tr_yz_yyyz[i] + 18.0 * tr_yyz_yyz[i] - 12.0 * tr_yyz_yyyyz[i] * tke_0 + 6.0 * tr_yyyz_yz[i] - 14.0 * tr_yyyz_yyyz[i] * tbe_0 - 14.0 * tr_yyyz_yyyz[i] * tke_0 + 4.0 * tr_yyyz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yyyyz_yyz[i] * tbe_0 + 8.0 * tr_yyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_yyzz[i] = 6.0 * tr_yz_yyzz[i] + 12.0 * tr_yyz_yzz[i] - 12.0 * tr_yyz_yyyzz[i] * tke_0 + 2.0 * tr_yyyz_zz[i] - 14.0 * tr_yyyz_yyzz[i] * tbe_0 - 10.0 * tr_yyyz_yyzz[i] * tke_0 + 4.0 * tr_yyyz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_yzz[i] * tbe_0 + 8.0 * tr_yyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_yzzz[i] = 6.0 * tr_yz_yzzz[i] + 6.0 * tr_yyz_zzz[i] - 12.0 * tr_yyz_yyzzz[i] * tke_0 - 14.0 * tr_yyyz_yzzz[i] * tbe_0 - 6.0 * tr_yyyz_yzzz[i] * tke_0 + 4.0 * tr_yyyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_zzz[i] * tbe_0 + 8.0 * tr_yyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_zzzz[i] = 6.0 * tr_yz_zzzz[i] - 12.0 * tr_yyz_yzzzz[i] * tke_0 - 14.0 * tr_yyyz_zzzz[i] * tbe_0 - 2.0 * tr_yyyz_zzzz[i] * tke_0 + 4.0 * tr_yyyz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 855-870 components of targeted buffer : GG

    auto tr_0_0_yy_yyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 855);

    auto tr_0_0_yy_yyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 856);

    auto tr_0_0_yy_yyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 857);

    auto tr_0_0_yy_yyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 858);

    auto tr_0_0_yy_yyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 859);

    auto tr_0_0_yy_yyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 860);

    auto tr_0_0_yy_yyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 861);

    auto tr_0_0_yy_yyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 862);

    auto tr_0_0_yy_yyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 863);

    auto tr_0_0_yy_yyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 864);

    auto tr_0_0_yy_yyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 865);

    auto tr_0_0_yy_yyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 866);

    auto tr_0_0_yy_yyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 867);

    auto tr_0_0_yy_yyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 868);

    auto tr_0_0_yy_yyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 869);

    #pragma omp simd aligned(tr_0_0_yy_yyzz_xxxx, tr_0_0_yy_yyzz_xxxy, tr_0_0_yy_yyzz_xxxz, tr_0_0_yy_yyzz_xxyy, tr_0_0_yy_yyzz_xxyz, tr_0_0_yy_yyzz_xxzz, tr_0_0_yy_yyzz_xyyy, tr_0_0_yy_yyzz_xyyz, tr_0_0_yy_yyzz_xyzz, tr_0_0_yy_yyzz_xzzz, tr_0_0_yy_yyzz_yyyy, tr_0_0_yy_yyzz_yyyz, tr_0_0_yy_yyzz_yyzz, tr_0_0_yy_yyzz_yzzz, tr_0_0_yy_yyzz_zzzz, tr_yyyyzz_xxxx, tr_yyyyzz_xxxy, tr_yyyyzz_xxxz, tr_yyyyzz_xxyy, tr_yyyyzz_xxyz, tr_yyyyzz_xxzz, tr_yyyyzz_xyyy, tr_yyyyzz_xyyz, tr_yyyyzz_xyzz, tr_yyyyzz_xzzz, tr_yyyyzz_yyyy, tr_yyyyzz_yyyz, tr_yyyyzz_yyzz, tr_yyyyzz_yzzz, tr_yyyyzz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxxxy, tr_yyyzz_xxxyy, tr_yyyzz_xxxyz, tr_yyyzz_xxy, tr_yyyzz_xxyyy, tr_yyyzz_xxyyz, tr_yyyzz_xxyzz, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyyyy, tr_yyyzz_xyyyz, tr_yyyzz_xyyzz, tr_yyyzz_xyz, tr_yyyzz_xyzzz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyyyy, tr_yyyzz_yyyyz, tr_yyyzz_yyyzz, tr_yyyzz_yyz, tr_yyyzz_yyzzz, tr_yyyzz_yzz, tr_yyyzz_yzzzz, tr_yyyzz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxxyy, tr_yyzz_xxxy, tr_yyzz_xxxyyy, tr_yyzz_xxxyyz, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyyyy, tr_yyzz_xxyyyz, tr_yyzz_xxyyzz, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyyyy, tr_yyzz_xyyyyz, tr_yyzz_xyyyzz, tr_yyzz_xyyz, tr_yyzz_xyyzzz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyyyy, tr_yyzz_yyyyyz, tr_yyzz_yyyyzz, tr_yyzz_yyyz, tr_yyzz_yyyzzz, tr_yyzz_yyzz, tr_yyzz_yyzzzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyzz_xxxx[i] = 2.0 * tr_zz_xxxx[i] - 8.0 * tr_yzz_xxxxy[i] * tke_0 - 10.0 * tr_yyzz_xxxx[i] * tbe_0 - 2.0 * tr_yyzz_xxxx[i] * tke_0 + 4.0 * tr_yyzz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xxxy[i] = 2.0 * tr_zz_xxxy[i] + 4.0 * tr_yzz_xxx[i] - 8.0 * tr_yzz_xxxyy[i] * tke_0 - 10.0 * tr_yyzz_xxxy[i] * tbe_0 - 6.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xxx[i] * tbe_0 + 8.0 * tr_yyyzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xxxz[i] = 2.0 * tr_zz_xxxz[i] - 8.0 * tr_yzz_xxxyz[i] * tke_0 - 10.0 * tr_yyzz_xxxz[i] * tbe_0 - 2.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xxyy[i] = 2.0 * tr_zz_xxyy[i] + 8.0 * tr_yzz_xxy[i] - 8.0 * tr_yzz_xxyyy[i] * tke_0 + 2.0 * tr_yyzz_xx[i] - 10.0 * tr_yyzz_xxyy[i] * tbe_0 - 10.0 * tr_yyzz_xxyy[i] * tke_0 + 4.0 * tr_yyzz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_xxy[i] * tbe_0 + 8.0 * tr_yyyzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xxyz[i] = 2.0 * tr_zz_xxyz[i] + 4.0 * tr_yzz_xxz[i] - 8.0 * tr_yzz_xxyyz[i] * tke_0 - 10.0 * tr_yyzz_xxyz[i] * tbe_0 - 6.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xxz[i] * tbe_0 + 8.0 * tr_yyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xxzz[i] = 2.0 * tr_zz_xxzz[i] - 8.0 * tr_yzz_xxyzz[i] * tke_0 - 10.0 * tr_yyzz_xxzz[i] * tbe_0 - 2.0 * tr_yyzz_xxzz[i] * tke_0 + 4.0 * tr_yyzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xyyy[i] = 2.0 * tr_zz_xyyy[i] + 12.0 * tr_yzz_xyy[i] - 8.0 * tr_yzz_xyyyy[i] * tke_0 + 6.0 * tr_yyzz_xy[i] - 10.0 * tr_yyzz_xyyy[i] * tbe_0 - 14.0 * tr_yyzz_xyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyyzz_xyy[i] * tbe_0 + 8.0 * tr_yyyzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xyyz[i] = 2.0 * tr_zz_xyyz[i] + 8.0 * tr_yzz_xyz[i] - 8.0 * tr_yzz_xyyyz[i] * tke_0 + 2.0 * tr_yyzz_xz[i] - 10.0 * tr_yyzz_xyyz[i] * tbe_0 - 10.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_xyz[i] * tbe_0 + 8.0 * tr_yyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xyzz[i] = 2.0 * tr_zz_xyzz[i] + 4.0 * tr_yzz_xzz[i] - 8.0 * tr_yzz_xyyzz[i] * tke_0 - 10.0 * tr_yyzz_xyzz[i] * tbe_0 - 6.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xzz[i] * tbe_0 + 8.0 * tr_yyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xzzz[i] = 2.0 * tr_zz_xzzz[i] - 8.0 * tr_yzz_xyzzz[i] * tke_0 - 10.0 * tr_yyzz_xzzz[i] * tbe_0 - 2.0 * tr_yyzz_xzzz[i] * tke_0 + 4.0 * tr_yyzz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_yyyy[i] = 2.0 * tr_zz_yyyy[i] + 16.0 * tr_yzz_yyy[i] - 8.0 * tr_yzz_yyyyy[i] * tke_0 + 12.0 * tr_yyzz_yy[i] - 10.0 * tr_yyzz_yyyy[i] * tbe_0 - 18.0 * tr_yyzz_yyyy[i] * tke_0 + 4.0 * tr_yyzz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yyyzz_yyy[i] * tbe_0 + 8.0 * tr_yyyzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_yyyz[i] = 2.0 * tr_zz_yyyz[i] + 12.0 * tr_yzz_yyz[i] - 8.0 * tr_yzz_yyyyz[i] * tke_0 + 6.0 * tr_yyzz_yz[i] - 10.0 * tr_yyzz_yyyz[i] * tbe_0 - 14.0 * tr_yyzz_yyyz[i] * tke_0 + 4.0 * tr_yyzz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yyyzz_yyz[i] * tbe_0 + 8.0 * tr_yyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_yyzz[i] = 2.0 * tr_zz_yyzz[i] + 8.0 * tr_yzz_yzz[i] - 8.0 * tr_yzz_yyyzz[i] * tke_0 + 2.0 * tr_yyzz_zz[i] - 10.0 * tr_yyzz_yyzz[i] * tbe_0 - 10.0 * tr_yyzz_yyzz[i] * tke_0 + 4.0 * tr_yyzz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_yzz[i] * tbe_0 + 8.0 * tr_yyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_yzzz[i] = 2.0 * tr_zz_yzzz[i] + 4.0 * tr_yzz_zzz[i] - 8.0 * tr_yzz_yyzzz[i] * tke_0 - 10.0 * tr_yyzz_yzzz[i] * tbe_0 - 6.0 * tr_yyzz_yzzz[i] * tke_0 + 4.0 * tr_yyzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_zzz[i] * tbe_0 + 8.0 * tr_yyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_zzzz[i] = 2.0 * tr_zz_zzzz[i] - 8.0 * tr_yzz_yzzzz[i] * tke_0 - 10.0 * tr_yyzz_zzzz[i] * tbe_0 - 2.0 * tr_yyzz_zzzz[i] * tke_0 + 4.0 * tr_yyzz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 870-885 components of targeted buffer : GG

    auto tr_0_0_yy_yzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 870);

    auto tr_0_0_yy_yzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 871);

    auto tr_0_0_yy_yzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 872);

    auto tr_0_0_yy_yzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 873);

    auto tr_0_0_yy_yzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 874);

    auto tr_0_0_yy_yzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 875);

    auto tr_0_0_yy_yzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 876);

    auto tr_0_0_yy_yzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 877);

    auto tr_0_0_yy_yzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 878);

    auto tr_0_0_yy_yzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 879);

    auto tr_0_0_yy_yzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 880);

    auto tr_0_0_yy_yzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 881);

    auto tr_0_0_yy_yzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 882);

    auto tr_0_0_yy_yzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 883);

    auto tr_0_0_yy_yzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 884);

    #pragma omp simd aligned(tr_0_0_yy_yzzz_xxxx, tr_0_0_yy_yzzz_xxxy, tr_0_0_yy_yzzz_xxxz, tr_0_0_yy_yzzz_xxyy, tr_0_0_yy_yzzz_xxyz, tr_0_0_yy_yzzz_xxzz, tr_0_0_yy_yzzz_xyyy, tr_0_0_yy_yzzz_xyyz, tr_0_0_yy_yzzz_xyzz, tr_0_0_yy_yzzz_xzzz, tr_0_0_yy_yzzz_yyyy, tr_0_0_yy_yzzz_yyyz, tr_0_0_yy_yzzz_yyzz, tr_0_0_yy_yzzz_yzzz, tr_0_0_yy_yzzz_zzzz, tr_yyyzzz_xxxx, tr_yyyzzz_xxxy, tr_yyyzzz_xxxz, tr_yyyzzz_xxyy, tr_yyyzzz_xxyz, tr_yyyzzz_xxzz, tr_yyyzzz_xyyy, tr_yyyzzz_xyyz, tr_yyyzzz_xyzz, tr_yyyzzz_xzzz, tr_yyyzzz_yyyy, tr_yyyzzz_yyyz, tr_yyyzzz_yyzz, tr_yyyzzz_yzzz, tr_yyyzzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxxxy, tr_yyzzz_xxxyy, tr_yyzzz_xxxyz, tr_yyzzz_xxy, tr_yyzzz_xxyyy, tr_yyzzz_xxyyz, tr_yyzzz_xxyzz, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyyyy, tr_yyzzz_xyyyz, tr_yyzzz_xyyzz, tr_yyzzz_xyz, tr_yyzzz_xyzzz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyyyy, tr_yyzzz_yyyyz, tr_yyzzz_yyyzz, tr_yyzzz_yyz, tr_yyzzz_yyzzz, tr_yyzzz_yzz, tr_yyzzz_yzzzz, tr_yyzzz_zzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxxyy, tr_yzzz_xxxy, tr_yzzz_xxxyyy, tr_yzzz_xxxyyz, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyyyy, tr_yzzz_xxyyyz, tr_yzzz_xxyyzz, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyyyy, tr_yzzz_xyyyyz, tr_yzzz_xyyyzz, tr_yzzz_xyyz, tr_yzzz_xyyzzz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyyyy, tr_yzzz_yyyyyz, tr_yzzz_yyyyzz, tr_yzzz_yyyz, tr_yzzz_yyyzzz, tr_yzzz_yyzz, tr_yzzz_yyzzzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxy, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyyyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yzzz_xxxx[i] = -4.0 * tr_zzz_xxxxy[i] * tke_0 - 6.0 * tr_yzzz_xxxx[i] * tbe_0 - 2.0 * tr_yzzz_xxxx[i] * tke_0 + 4.0 * tr_yzzz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xxxy[i] = 2.0 * tr_zzz_xxx[i] - 4.0 * tr_zzz_xxxyy[i] * tke_0 - 6.0 * tr_yzzz_xxxy[i] * tbe_0 - 6.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xxx[i] * tbe_0 + 8.0 * tr_yyzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xxxz[i] = -4.0 * tr_zzz_xxxyz[i] * tke_0 - 6.0 * tr_yzzz_xxxz[i] * tbe_0 - 2.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xxyy[i] = 4.0 * tr_zzz_xxy[i] - 4.0 * tr_zzz_xxyyy[i] * tke_0 + 2.0 * tr_yzzz_xx[i] - 6.0 * tr_yzzz_xxyy[i] * tbe_0 - 10.0 * tr_yzzz_xxyy[i] * tke_0 + 4.0 * tr_yzzz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_xxy[i] * tbe_0 + 8.0 * tr_yyzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xxyz[i] = 2.0 * tr_zzz_xxz[i] - 4.0 * tr_zzz_xxyyz[i] * tke_0 - 6.0 * tr_yzzz_xxyz[i] * tbe_0 - 6.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xxz[i] * tbe_0 + 8.0 * tr_yyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xxzz[i] = -4.0 * tr_zzz_xxyzz[i] * tke_0 - 6.0 * tr_yzzz_xxzz[i] * tbe_0 - 2.0 * tr_yzzz_xxzz[i] * tke_0 + 4.0 * tr_yzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xyyy[i] = 6.0 * tr_zzz_xyy[i] - 4.0 * tr_zzz_xyyyy[i] * tke_0 + 6.0 * tr_yzzz_xy[i] - 6.0 * tr_yzzz_xyyy[i] * tbe_0 - 14.0 * tr_yzzz_xyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yyzzz_xyy[i] * tbe_0 + 8.0 * tr_yyzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xyyz[i] = 4.0 * tr_zzz_xyz[i] - 4.0 * tr_zzz_xyyyz[i] * tke_0 + 2.0 * tr_yzzz_xz[i] - 6.0 * tr_yzzz_xyyz[i] * tbe_0 - 10.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_xyz[i] * tbe_0 + 8.0 * tr_yyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xyzz[i] = 2.0 * tr_zzz_xzz[i] - 4.0 * tr_zzz_xyyzz[i] * tke_0 - 6.0 * tr_yzzz_xyzz[i] * tbe_0 - 6.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xzz[i] * tbe_0 + 8.0 * tr_yyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xzzz[i] = -4.0 * tr_zzz_xyzzz[i] * tke_0 - 6.0 * tr_yzzz_xzzz[i] * tbe_0 - 2.0 * tr_yzzz_xzzz[i] * tke_0 + 4.0 * tr_yzzz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_yyyy[i] = 8.0 * tr_zzz_yyy[i] - 4.0 * tr_zzz_yyyyy[i] * tke_0 + 12.0 * tr_yzzz_yy[i] - 6.0 * tr_yzzz_yyyy[i] * tbe_0 - 18.0 * tr_yzzz_yyyy[i] * tke_0 + 4.0 * tr_yzzz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yyzzz_yyy[i] * tbe_0 + 8.0 * tr_yyzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_yyyz[i] = 6.0 * tr_zzz_yyz[i] - 4.0 * tr_zzz_yyyyz[i] * tke_0 + 6.0 * tr_yzzz_yz[i] - 6.0 * tr_yzzz_yyyz[i] * tbe_0 - 14.0 * tr_yzzz_yyyz[i] * tke_0 + 4.0 * tr_yzzz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yyzzz_yyz[i] * tbe_0 + 8.0 * tr_yyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_yyzz[i] = 4.0 * tr_zzz_yzz[i] - 4.0 * tr_zzz_yyyzz[i] * tke_0 + 2.0 * tr_yzzz_zz[i] - 6.0 * tr_yzzz_yyzz[i] * tbe_0 - 10.0 * tr_yzzz_yyzz[i] * tke_0 + 4.0 * tr_yzzz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_yzz[i] * tbe_0 + 8.0 * tr_yyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_yzzz[i] = 2.0 * tr_zzz_zzz[i] - 4.0 * tr_zzz_yyzzz[i] * tke_0 - 6.0 * tr_yzzz_yzzz[i] * tbe_0 - 6.0 * tr_yzzz_yzzz[i] * tke_0 + 4.0 * tr_yzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_zzz[i] * tbe_0 + 8.0 * tr_yyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_zzzz[i] = -4.0 * tr_zzz_yzzzz[i] * tke_0 - 6.0 * tr_yzzz_zzzz[i] * tbe_0 - 2.0 * tr_yzzz_zzzz[i] * tke_0 + 4.0 * tr_yzzz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 885-900 components of targeted buffer : GG

    auto tr_0_0_yy_zzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 885);

    auto tr_0_0_yy_zzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 886);

    auto tr_0_0_yy_zzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 887);

    auto tr_0_0_yy_zzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 888);

    auto tr_0_0_yy_zzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 889);

    auto tr_0_0_yy_zzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 890);

    auto tr_0_0_yy_zzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 891);

    auto tr_0_0_yy_zzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 892);

    auto tr_0_0_yy_zzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 893);

    auto tr_0_0_yy_zzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 894);

    auto tr_0_0_yy_zzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 895);

    auto tr_0_0_yy_zzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 896);

    auto tr_0_0_yy_zzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 897);

    auto tr_0_0_yy_zzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 898);

    auto tr_0_0_yy_zzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 899);

    #pragma omp simd aligned(tr_0_0_yy_zzzz_xxxx, tr_0_0_yy_zzzz_xxxy, tr_0_0_yy_zzzz_xxxz, tr_0_0_yy_zzzz_xxyy, tr_0_0_yy_zzzz_xxyz, tr_0_0_yy_zzzz_xxzz, tr_0_0_yy_zzzz_xyyy, tr_0_0_yy_zzzz_xyyz, tr_0_0_yy_zzzz_xyzz, tr_0_0_yy_zzzz_xzzz, tr_0_0_yy_zzzz_yyyy, tr_0_0_yy_zzzz_yyyz, tr_0_0_yy_zzzz_yyzz, tr_0_0_yy_zzzz_yzzz, tr_0_0_yy_zzzz_zzzz, tr_yyzzzz_xxxx, tr_yyzzzz_xxxy, tr_yyzzzz_xxxz, tr_yyzzzz_xxyy, tr_yyzzzz_xxyz, tr_yyzzzz_xxzz, tr_yyzzzz_xyyy, tr_yyzzzz_xyyz, tr_yyzzzz_xyzz, tr_yyzzzz_xzzz, tr_yyzzzz_yyyy, tr_yyzzzz_yyyz, tr_yyzzzz_yyzz, tr_yyzzzz_yzzz, tr_yyzzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxxxy, tr_yzzzz_xxxyy, tr_yzzzz_xxxyz, tr_yzzzz_xxy, tr_yzzzz_xxyyy, tr_yzzzz_xxyyz, tr_yzzzz_xxyzz, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyyyy, tr_yzzzz_xyyyz, tr_yzzzz_xyyzz, tr_yzzzz_xyz, tr_yzzzz_xyzzz, tr_yzzzz_xzz, tr_yzzzz_yyy, tr_yzzzz_yyyyy, tr_yzzzz_yyyyz, tr_yzzzz_yyyzz, tr_yzzzz_yyz, tr_yzzzz_yyzzz, tr_yzzzz_yzz, tr_yzzzz_yzzzz, tr_yzzzz_zzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxxyy, tr_zzzz_xxxy, tr_zzzz_xxxyyy, tr_zzzz_xxxyyz, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyyyy, tr_zzzz_xxyyyz, tr_zzzz_xxyyzz, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyyyy, tr_zzzz_xyyyyz, tr_zzzz_xyyyzz, tr_zzzz_xyyz, tr_zzzz_xyyzzz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyyyy, tr_zzzz_yyyyyz, tr_zzzz_yyyyzz, tr_zzzz_yyyz, tr_zzzz_yyyzzz, tr_zzzz_yyzz, tr_zzzz_yyzzzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_zzzz_xxxx[i] = -2.0 * tr_zzzz_xxxx[i] * tbe_0 - 2.0 * tr_zzzz_xxxx[i] * tke_0 + 4.0 * tr_zzzz_xxxxyy[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xxxy[i] = -2.0 * tr_zzzz_xxxy[i] * tbe_0 - 6.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxyyy[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xxx[i] * tbe_0 + 8.0 * tr_yzzzz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xxxz[i] = -2.0 * tr_zzzz_xxxz[i] * tbe_0 - 2.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxyyz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xxyy[i] = 2.0 * tr_zzzz_xx[i] - 2.0 * tr_zzzz_xxyy[i] * tbe_0 - 10.0 * tr_zzzz_xxyy[i] * tke_0 + 4.0 * tr_zzzz_xxyyyy[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_xxy[i] * tbe_0 + 8.0 * tr_yzzzz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xxyz[i] = -2.0 * tr_zzzz_xxyz[i] * tbe_0 - 6.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xxz[i] * tbe_0 + 8.0 * tr_yzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xxzz[i] = -2.0 * tr_zzzz_xxzz[i] * tbe_0 - 2.0 * tr_zzzz_xxzz[i] * tke_0 + 4.0 * tr_zzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xyyy[i] = 6.0 * tr_zzzz_xy[i] - 2.0 * tr_zzzz_xyyy[i] * tbe_0 - 14.0 * tr_zzzz_xyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyyy[i] * tke_0 * tke_0 - 12.0 * tr_yzzzz_xyy[i] * tbe_0 + 8.0 * tr_yzzzz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xyyz[i] = 2.0 * tr_zzzz_xz[i] - 2.0 * tr_zzzz_xyyz[i] * tbe_0 - 10.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_xyz[i] * tbe_0 + 8.0 * tr_yzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xyzz[i] = -2.0 * tr_zzzz_xyzz[i] * tbe_0 - 6.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xzz[i] * tbe_0 + 8.0 * tr_yzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xzzz[i] = -2.0 * tr_zzzz_xzzz[i] * tbe_0 - 2.0 * tr_zzzz_xzzz[i] * tke_0 + 4.0 * tr_zzzz_xyyzzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_yyyy[i] = 12.0 * tr_zzzz_yy[i] - 2.0 * tr_zzzz_yyyy[i] * tbe_0 - 18.0 * tr_zzzz_yyyy[i] * tke_0 + 4.0 * tr_zzzz_yyyyyy[i] * tke_0 * tke_0 - 16.0 * tr_yzzzz_yyy[i] * tbe_0 + 8.0 * tr_yzzzz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_yyyz[i] = 6.0 * tr_zzzz_yz[i] - 2.0 * tr_zzzz_yyyz[i] * tbe_0 - 14.0 * tr_zzzz_yyyz[i] * tke_0 + 4.0 * tr_zzzz_yyyyyz[i] * tke_0 * tke_0 - 12.0 * tr_yzzzz_yyz[i] * tbe_0 + 8.0 * tr_yzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_yyzz[i] = 2.0 * tr_zzzz_zz[i] - 2.0 * tr_zzzz_yyzz[i] * tbe_0 - 10.0 * tr_zzzz_yyzz[i] * tke_0 + 4.0 * tr_zzzz_yyyyzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_yzz[i] * tbe_0 + 8.0 * tr_yzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_yzzz[i] = -2.0 * tr_zzzz_yzzz[i] * tbe_0 - 6.0 * tr_zzzz_yzzz[i] * tke_0 + 4.0 * tr_zzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_zzz[i] * tbe_0 + 8.0 * tr_yzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_zzzz[i] = -2.0 * tr_zzzz_zzzz[i] * tbe_0 - 2.0 * tr_zzzz_zzzz[i] * tke_0 + 4.0 * tr_zzzz_yyzzzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 900-915 components of targeted buffer : GG

    auto tr_0_0_yz_xxxx_xxxx = pbuffer.data(idx_op_geom_020_gg + 900);

    auto tr_0_0_yz_xxxx_xxxy = pbuffer.data(idx_op_geom_020_gg + 901);

    auto tr_0_0_yz_xxxx_xxxz = pbuffer.data(idx_op_geom_020_gg + 902);

    auto tr_0_0_yz_xxxx_xxyy = pbuffer.data(idx_op_geom_020_gg + 903);

    auto tr_0_0_yz_xxxx_xxyz = pbuffer.data(idx_op_geom_020_gg + 904);

    auto tr_0_0_yz_xxxx_xxzz = pbuffer.data(idx_op_geom_020_gg + 905);

    auto tr_0_0_yz_xxxx_xyyy = pbuffer.data(idx_op_geom_020_gg + 906);

    auto tr_0_0_yz_xxxx_xyyz = pbuffer.data(idx_op_geom_020_gg + 907);

    auto tr_0_0_yz_xxxx_xyzz = pbuffer.data(idx_op_geom_020_gg + 908);

    auto tr_0_0_yz_xxxx_xzzz = pbuffer.data(idx_op_geom_020_gg + 909);

    auto tr_0_0_yz_xxxx_yyyy = pbuffer.data(idx_op_geom_020_gg + 910);

    auto tr_0_0_yz_xxxx_yyyz = pbuffer.data(idx_op_geom_020_gg + 911);

    auto tr_0_0_yz_xxxx_yyzz = pbuffer.data(idx_op_geom_020_gg + 912);

    auto tr_0_0_yz_xxxx_yzzz = pbuffer.data(idx_op_geom_020_gg + 913);

    auto tr_0_0_yz_xxxx_zzzz = pbuffer.data(idx_op_geom_020_gg + 914);

    #pragma omp simd aligned(tr_0_0_yz_xxxx_xxxx, tr_0_0_yz_xxxx_xxxy, tr_0_0_yz_xxxx_xxxz, tr_0_0_yz_xxxx_xxyy, tr_0_0_yz_xxxx_xxyz, tr_0_0_yz_xxxx_xxzz, tr_0_0_yz_xxxx_xyyy, tr_0_0_yz_xxxx_xyyz, tr_0_0_yz_xxxx_xyzz, tr_0_0_yz_xxxx_xzzz, tr_0_0_yz_xxxx_yyyy, tr_0_0_yz_xxxx_yyyz, tr_0_0_yz_xxxx_yyzz, tr_0_0_yz_xxxx_yzzz, tr_0_0_yz_xxxx_zzzz, tr_xxxx_xx, tr_xxxx_xxxxyz, tr_xxxx_xxxy, tr_xxxx_xxxyyz, tr_xxxx_xxxyzz, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyyyz, tr_xxxx_xxyyzz, tr_xxxx_xxyz, tr_xxxx_xxyzzz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyyyz, tr_xxxx_xyyyzz, tr_xxxx_xyyz, tr_xxxx_xyyzzz, tr_xxxx_xyzz, tr_xxxx_xyzzzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyyyz, tr_xxxx_yyyyzz, tr_xxxx_yyyz, tr_xxxx_yyyzzz, tr_xxxx_yyzz, tr_xxxx_yyzzzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_yzzzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxy_xxx, tr_xxxxy_xxxxz, tr_xxxxy_xxxyz, tr_xxxxy_xxxzz, tr_xxxxy_xxy, tr_xxxxy_xxyyz, tr_xxxxy_xxyzz, tr_xxxxy_xxz, tr_xxxxy_xxzzz, tr_xxxxy_xyy, tr_xxxxy_xyyyz, tr_xxxxy_xyyzz, tr_xxxxy_xyz, tr_xxxxy_xyzzz, tr_xxxxy_xzz, tr_xxxxy_xzzzz, tr_xxxxy_yyy, tr_xxxxy_yyyyz, tr_xxxxy_yyyzz, tr_xxxxy_yyz, tr_xxxxy_yyzzz, tr_xxxxy_yzz, tr_xxxxy_yzzzz, tr_xxxxy_zzz, tr_xxxxy_zzzzz, tr_xxxxyz_xxxx, tr_xxxxyz_xxxy, tr_xxxxyz_xxxz, tr_xxxxyz_xxyy, tr_xxxxyz_xxyz, tr_xxxxyz_xxzz, tr_xxxxyz_xyyy, tr_xxxxyz_xyyz, tr_xxxxyz_xyzz, tr_xxxxyz_xzzz, tr_xxxxyz_yyyy, tr_xxxxyz_yyyz, tr_xxxxyz_yyzz, tr_xxxxyz_yzzz, tr_xxxxyz_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxxxy, tr_xxxxz_xxxyy, tr_xxxxz_xxxyz, tr_xxxxz_xxy, tr_xxxxz_xxyyy, tr_xxxxz_xxyyz, tr_xxxxz_xxyzz, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyyyy, tr_xxxxz_xyyyz, tr_xxxxz_xyyzz, tr_xxxxz_xyz, tr_xxxxz_xyzzz, tr_xxxxz_xzz, tr_xxxxz_yyy, tr_xxxxz_yyyyy, tr_xxxxz_yyyyz, tr_xxxxz_yyyzz, tr_xxxxz_yyz, tr_xxxxz_yyzzz, tr_xxxxz_yzz, tr_xxxxz_yzzzz, tr_xxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxxx_xxxx[i] = 4.0 * tr_xxxx_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xxxy[i] = -2.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xxxz[i] = -2.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xxyy[i] = -4.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xxyz[i] = tr_xxxx_xx[i] - 2.0 * tr_xxxx_xxzz[i] * tke_0 - 2.0 * tr_xxxx_xxyy[i] * tke_0 + 4.0 * tr_xxxx_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xxzz[i] = -4.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xyyy[i] = -6.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xyyz[i] = 2.0 * tr_xxxx_xy[i] - 4.0 * tr_xxxx_xyzz[i] * tke_0 - 2.0 * tr_xxxx_xyyy[i] * tke_0 + 4.0 * tr_xxxx_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xyzz[i] = 2.0 * tr_xxxx_xz[i] - 2.0 * tr_xxxx_xzzz[i] * tke_0 - 4.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xzzz[i] = -6.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxxy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_yyyy[i] = -8.0 * tr_xxxx_yyyz[i] * tke_0 + 4.0 * tr_xxxx_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_yyyz[i] = 3.0 * tr_xxxx_yy[i] - 6.0 * tr_xxxx_yyzz[i] * tke_0 - 2.0 * tr_xxxx_yyyy[i] * tke_0 + 4.0 * tr_xxxx_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xxxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxxz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_yyzz[i] = 4.0 * tr_xxxx_yz[i] - 4.0 * tr_xxxx_yzzz[i] * tke_0 - 4.0 * tr_xxxx_yyyz[i] * tke_0 + 4.0 * tr_xxxx_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxxz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxxy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_yzzz[i] = 3.0 * tr_xxxx_zz[i] - 2.0 * tr_xxxx_zzzz[i] * tke_0 - 6.0 * tr_xxxx_yyzz[i] * tke_0 + 4.0 * tr_xxxx_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_zzzz[i] = -8.0 * tr_xxxx_yzzz[i] * tke_0 + 4.0 * tr_xxxx_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 915-930 components of targeted buffer : GG

    auto tr_0_0_yz_xxxy_xxxx = pbuffer.data(idx_op_geom_020_gg + 915);

    auto tr_0_0_yz_xxxy_xxxy = pbuffer.data(idx_op_geom_020_gg + 916);

    auto tr_0_0_yz_xxxy_xxxz = pbuffer.data(idx_op_geom_020_gg + 917);

    auto tr_0_0_yz_xxxy_xxyy = pbuffer.data(idx_op_geom_020_gg + 918);

    auto tr_0_0_yz_xxxy_xxyz = pbuffer.data(idx_op_geom_020_gg + 919);

    auto tr_0_0_yz_xxxy_xxzz = pbuffer.data(idx_op_geom_020_gg + 920);

    auto tr_0_0_yz_xxxy_xyyy = pbuffer.data(idx_op_geom_020_gg + 921);

    auto tr_0_0_yz_xxxy_xyyz = pbuffer.data(idx_op_geom_020_gg + 922);

    auto tr_0_0_yz_xxxy_xyzz = pbuffer.data(idx_op_geom_020_gg + 923);

    auto tr_0_0_yz_xxxy_xzzz = pbuffer.data(idx_op_geom_020_gg + 924);

    auto tr_0_0_yz_xxxy_yyyy = pbuffer.data(idx_op_geom_020_gg + 925);

    auto tr_0_0_yz_xxxy_yyyz = pbuffer.data(idx_op_geom_020_gg + 926);

    auto tr_0_0_yz_xxxy_yyzz = pbuffer.data(idx_op_geom_020_gg + 927);

    auto tr_0_0_yz_xxxy_yzzz = pbuffer.data(idx_op_geom_020_gg + 928);

    auto tr_0_0_yz_xxxy_zzzz = pbuffer.data(idx_op_geom_020_gg + 929);

    #pragma omp simd aligned(tr_0_0_yz_xxxy_xxxx, tr_0_0_yz_xxxy_xxxy, tr_0_0_yz_xxxy_xxxz, tr_0_0_yz_xxxy_xxyy, tr_0_0_yz_xxxy_xxyz, tr_0_0_yz_xxxy_xxzz, tr_0_0_yz_xxxy_xyyy, tr_0_0_yz_xxxy_xyyz, tr_0_0_yz_xxxy_xyzz, tr_0_0_yz_xxxy_xzzz, tr_0_0_yz_xxxy_yyyy, tr_0_0_yz_xxxy_yyyz, tr_0_0_yz_xxxy_yyzz, tr_0_0_yz_xxxy_yzzz, tr_0_0_yz_xxxy_zzzz, tr_xxx_xxx, tr_xxx_xxxxz, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxx_zzzzz, tr_xxxy_xx, tr_xxxy_xxxxyz, tr_xxxy_xxxy, tr_xxxy_xxxyyz, tr_xxxy_xxxyzz, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyyyz, tr_xxxy_xxyyzz, tr_xxxy_xxyz, tr_xxxy_xxyzzz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyyyz, tr_xxxy_xyyyzz, tr_xxxy_xyyz, tr_xxxy_xyyzzz, tr_xxxy_xyzz, tr_xxxy_xyzzzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyyyz, tr_xxxy_yyyyzz, tr_xxxy_yyyz, tr_xxxy_yyyzzz, tr_xxxy_yyzz, tr_xxxy_yyzzzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_yzzzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyy_xxx, tr_xxxyy_xxxxz, tr_xxxyy_xxxyz, tr_xxxyy_xxxzz, tr_xxxyy_xxy, tr_xxxyy_xxyyz, tr_xxxyy_xxyzz, tr_xxxyy_xxz, tr_xxxyy_xxzzz, tr_xxxyy_xyy, tr_xxxyy_xyyyz, tr_xxxyy_xyyzz, tr_xxxyy_xyz, tr_xxxyy_xyzzz, tr_xxxyy_xzz, tr_xxxyy_xzzzz, tr_xxxyy_yyy, tr_xxxyy_yyyyz, tr_xxxyy_yyyzz, tr_xxxyy_yyz, tr_xxxyy_yyzzz, tr_xxxyy_yzz, tr_xxxyy_yzzzz, tr_xxxyy_zzz, tr_xxxyy_zzzzz, tr_xxxyyz_xxxx, tr_xxxyyz_xxxy, tr_xxxyyz_xxxz, tr_xxxyyz_xxyy, tr_xxxyyz_xxyz, tr_xxxyyz_xxzz, tr_xxxyyz_xyyy, tr_xxxyyz_xyyz, tr_xxxyyz_xyzz, tr_xxxyyz_xzzz, tr_xxxyyz_yyyy, tr_xxxyyz_yyyz, tr_xxxyyz_yyzz, tr_xxxyyz_yzzz, tr_xxxyyz_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxy, tr_xxxyz_xxxyy, tr_xxxyz_xxxyz, tr_xxxyz_xxy, tr_xxxyz_xxyyy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyyyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyyyy, tr_xxxyz_yyyyz, tr_xxxyz_yyyzz, tr_xxxyz_yyz, tr_xxxyz_yyzzz, tr_xxxyz_yzz, tr_xxxyz_yzzzz, tr_xxxyz_zzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxxy_xxxx[i] = -2.0 * tr_xxx_xxxxz[i] * tke_0 - 2.0 * tr_xxxz_xxxx[i] * tbe_0 + 4.0 * tr_xxxy_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xxxy[i] = -2.0 * tr_xxx_xxxyz[i] * tke_0 - 2.0 * tr_xxxz_xxxy[i] * tbe_0 - 2.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xxxz[i] = tr_xxx_xxx[i] - 2.0 * tr_xxx_xxxzz[i] * tke_0 - 2.0 * tr_xxxz_xxxz[i] * tbe_0 - 2.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xxyy[i] = -2.0 * tr_xxx_xxyyz[i] * tke_0 - 2.0 * tr_xxxz_xxyy[i] * tbe_0 - 4.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xxyz[i] = tr_xxx_xxy[i] - 2.0 * tr_xxx_xxyzz[i] * tke_0 - 2.0 * tr_xxxz_xxyz[i] * tbe_0 + tr_xxxy_xx[i] - 2.0 * tr_xxxy_xxzz[i] * tke_0 - 2.0 * tr_xxxy_xxyy[i] * tke_0 + 4.0 * tr_xxxy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_xxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xxzz[i] = 2.0 * tr_xxx_xxz[i] - 2.0 * tr_xxx_xxzzz[i] * tke_0 - 2.0 * tr_xxxz_xxzz[i] * tbe_0 - 4.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_xxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xyyy[i] = -2.0 * tr_xxx_xyyyz[i] * tke_0 - 2.0 * tr_xxxz_xyyy[i] * tbe_0 - 6.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xyyz[i] = tr_xxx_xyy[i] - 2.0 * tr_xxx_xyyzz[i] * tke_0 - 2.0 * tr_xxxz_xyyz[i] * tbe_0 + 2.0 * tr_xxxy_xy[i] - 4.0 * tr_xxxy_xyzz[i] * tke_0 - 2.0 * tr_xxxy_xyyy[i] * tke_0 + 4.0 * tr_xxxy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_xyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xyzz[i] = 2.0 * tr_xxx_xyz[i] - 2.0 * tr_xxx_xyzzz[i] * tke_0 - 2.0 * tr_xxxz_xyzz[i] * tbe_0 + 2.0 * tr_xxxy_xz[i] - 2.0 * tr_xxxy_xzzz[i] * tke_0 - 4.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_xyz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xzzz[i] = 3.0 * tr_xxx_xzz[i] - 2.0 * tr_xxx_xzzzz[i] * tke_0 - 2.0 * tr_xxxz_xzzz[i] * tbe_0 - 6.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxyy_xzz[i] * tbe_0 + 4.0 * tr_xxxyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_yyyy[i] = -2.0 * tr_xxx_yyyyz[i] * tke_0 - 2.0 * tr_xxxz_yyyy[i] * tbe_0 - 8.0 * tr_xxxy_yyyz[i] * tke_0 + 4.0 * tr_xxxy_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_yyyz[i] = tr_xxx_yyy[i] - 2.0 * tr_xxx_yyyzz[i] * tke_0 - 2.0 * tr_xxxz_yyyz[i] * tbe_0 + 3.0 * tr_xxxy_yy[i] - 6.0 * tr_xxxy_yyzz[i] * tke_0 - 2.0 * tr_xxxy_yyyy[i] * tke_0 + 4.0 * tr_xxxy_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xxxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_yyzz[i] = 2.0 * tr_xxx_yyz[i] - 2.0 * tr_xxx_yyzzz[i] * tke_0 - 2.0 * tr_xxxz_yyzz[i] * tbe_0 + 4.0 * tr_xxxy_yz[i] - 4.0 * tr_xxxy_yzzz[i] * tke_0 - 4.0 * tr_xxxy_yyyz[i] * tke_0 + 4.0 * tr_xxxy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_yyz[i] * tbe_0 + 4.0 * tr_xxxyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_yzzz[i] = 3.0 * tr_xxx_yzz[i] - 2.0 * tr_xxx_yzzzz[i] * tke_0 - 2.0 * tr_xxxz_yzzz[i] * tbe_0 + 3.0 * tr_xxxy_zz[i] - 2.0 * tr_xxxy_zzzz[i] * tke_0 - 6.0 * tr_xxxy_yyzz[i] * tke_0 + 4.0 * tr_xxxy_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxxyy_yzz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_zzzz[i] = 4.0 * tr_xxx_zzz[i] - 2.0 * tr_xxx_zzzzz[i] * tke_0 - 2.0 * tr_xxxz_zzzz[i] * tbe_0 - 8.0 * tr_xxxy_yzzz[i] * tke_0 + 4.0 * tr_xxxy_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxxyy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 930-945 components of targeted buffer : GG

    auto tr_0_0_yz_xxxz_xxxx = pbuffer.data(idx_op_geom_020_gg + 930);

    auto tr_0_0_yz_xxxz_xxxy = pbuffer.data(idx_op_geom_020_gg + 931);

    auto tr_0_0_yz_xxxz_xxxz = pbuffer.data(idx_op_geom_020_gg + 932);

    auto tr_0_0_yz_xxxz_xxyy = pbuffer.data(idx_op_geom_020_gg + 933);

    auto tr_0_0_yz_xxxz_xxyz = pbuffer.data(idx_op_geom_020_gg + 934);

    auto tr_0_0_yz_xxxz_xxzz = pbuffer.data(idx_op_geom_020_gg + 935);

    auto tr_0_0_yz_xxxz_xyyy = pbuffer.data(idx_op_geom_020_gg + 936);

    auto tr_0_0_yz_xxxz_xyyz = pbuffer.data(idx_op_geom_020_gg + 937);

    auto tr_0_0_yz_xxxz_xyzz = pbuffer.data(idx_op_geom_020_gg + 938);

    auto tr_0_0_yz_xxxz_xzzz = pbuffer.data(idx_op_geom_020_gg + 939);

    auto tr_0_0_yz_xxxz_yyyy = pbuffer.data(idx_op_geom_020_gg + 940);

    auto tr_0_0_yz_xxxz_yyyz = pbuffer.data(idx_op_geom_020_gg + 941);

    auto tr_0_0_yz_xxxz_yyzz = pbuffer.data(idx_op_geom_020_gg + 942);

    auto tr_0_0_yz_xxxz_yzzz = pbuffer.data(idx_op_geom_020_gg + 943);

    auto tr_0_0_yz_xxxz_zzzz = pbuffer.data(idx_op_geom_020_gg + 944);

    #pragma omp simd aligned(tr_0_0_yz_xxxz_xxxx, tr_0_0_yz_xxxz_xxxy, tr_0_0_yz_xxxz_xxxz, tr_0_0_yz_xxxz_xxyy, tr_0_0_yz_xxxz_xxyz, tr_0_0_yz_xxxz_xxzz, tr_0_0_yz_xxxz_xyyy, tr_0_0_yz_xxxz_xyyz, tr_0_0_yz_xxxz_xyzz, tr_0_0_yz_xxxz_xzzz, tr_0_0_yz_xxxz_yyyy, tr_0_0_yz_xxxz_yyyz, tr_0_0_yz_xxxz_yyzz, tr_0_0_yz_xxxz_yzzz, tr_0_0_yz_xxxz_zzzz, tr_xxx_xxx, tr_xxx_xxxxy, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyyyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxz, tr_xxxyz_xxxyz, tr_xxxyz_xxxzz, tr_xxxyz_xxy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xxzzz, tr_xxxyz_xyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_xzzzz, tr_xxxyz_yyy, tr_xxxyz_yyyyz, tr_xxxyz_yyyzz, tr_xxxyz_yyz, tr_xxxyz_yyzzz, tr_xxxyz_yzz, tr_xxxyz_yzzzz, tr_xxxyz_zzz, tr_xxxyz_zzzzz, tr_xxxyzz_xxxx, tr_xxxyzz_xxxy, tr_xxxyzz_xxxz, tr_xxxyzz_xxyy, tr_xxxyzz_xxyz, tr_xxxyzz_xxzz, tr_xxxyzz_xyyy, tr_xxxyzz_xyyz, tr_xxxyzz_xyzz, tr_xxxyzz_xzzz, tr_xxxyzz_yyyy, tr_xxxyzz_yyyz, tr_xxxyzz_yyzz, tr_xxxyzz_yzzz, tr_xxxyzz_zzzz, tr_xxxz_xx, tr_xxxz_xxxxyz, tr_xxxz_xxxy, tr_xxxz_xxxyyz, tr_xxxz_xxxyzz, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyyyz, tr_xxxz_xxyyzz, tr_xxxz_xxyz, tr_xxxz_xxyzzz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyyyz, tr_xxxz_xyyyzz, tr_xxxz_xyyz, tr_xxxz_xyyzzz, tr_xxxz_xyzz, tr_xxxz_xyzzzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyyyz, tr_xxxz_yyyyzz, tr_xxxz_yyyz, tr_xxxz_yyyzzz, tr_xxxz_yyzz, tr_xxxz_yyzzzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_yzzzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxxxy, tr_xxxzz_xxxyy, tr_xxxzz_xxxyz, tr_xxxzz_xxy, tr_xxxzz_xxyyy, tr_xxxzz_xxyyz, tr_xxxzz_xxyzz, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyyyy, tr_xxxzz_xyyyz, tr_xxxzz_xyyzz, tr_xxxzz_xyz, tr_xxxzz_xyzzz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyyyy, tr_xxxzz_yyyyz, tr_xxxzz_yyyzz, tr_xxxzz_yyz, tr_xxxzz_yyzzz, tr_xxxzz_yzz, tr_xxxzz_yzzzz, tr_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxxz_xxxx[i] = -2.0 * tr_xxx_xxxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xxxy[i] = tr_xxx_xxx[i] - 2.0 * tr_xxx_xxxyy[i] * tke_0 - 2.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xxxz[i] = -2.0 * tr_xxx_xxxyz[i] * tke_0 - 2.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxxz[i] * tbe_0 - 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xxyy[i] = 2.0 * tr_xxx_xxy[i] - 2.0 * tr_xxx_xxyyy[i] * tke_0 - 4.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxyy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xxyz[i] = tr_xxx_xxz[i] - 2.0 * tr_xxx_xxyyz[i] * tke_0 + tr_xxxz_xx[i] - 2.0 * tr_xxxz_xxzz[i] * tke_0 - 2.0 * tr_xxxz_xxyy[i] * tke_0 + 4.0 * tr_xxxz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_xxz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxyz[i] * tbe_0 - 2.0 * tr_xxxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xxzz[i] = -2.0 * tr_xxx_xxyzz[i] * tke_0 - 4.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxzz[i] * tbe_0 - 4.0 * tr_xxxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xyyy[i] = 3.0 * tr_xxx_xyy[i] - 2.0 * tr_xxx_xyyyy[i] * tke_0 - 6.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxxzz_xyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xyyz[i] = 2.0 * tr_xxx_xyz[i] - 2.0 * tr_xxx_xyyyz[i] * tke_0 + 2.0 * tr_xxxz_xy[i] - 4.0 * tr_xxxz_xyzz[i] * tke_0 - 2.0 * tr_xxxz_xyyy[i] * tke_0 + 4.0 * tr_xxxz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xyz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyyz[i] * tbe_0 - 2.0 * tr_xxxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xyzz[i] = tr_xxx_xzz[i] - 2.0 * tr_xxx_xyyzz[i] * tke_0 + 2.0 * tr_xxxz_xz[i] - 2.0 * tr_xxxz_xzzz[i] * tke_0 - 4.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_xzz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyzz[i] * tbe_0 - 4.0 * tr_xxxyz_xyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xzzz[i] = -2.0 * tr_xxx_xyzzz[i] * tke_0 - 6.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xzzz[i] * tbe_0 - 6.0 * tr_xxxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_yyyy[i] = 4.0 * tr_xxx_yyy[i] - 2.0 * tr_xxx_yyyyy[i] * tke_0 - 8.0 * tr_xxxz_yyyz[i] * tke_0 + 4.0 * tr_xxxz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_yyyz[i] = 3.0 * tr_xxx_yyz[i] - 2.0 * tr_xxx_yyyyz[i] * tke_0 + 3.0 * tr_xxxz_yy[i] - 6.0 * tr_xxxz_yyzz[i] * tke_0 - 2.0 * tr_xxxz_yyyy[i] * tke_0 + 4.0 * tr_xxxz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xxxzz_yyz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyyz[i] * tbe_0 - 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_yyzz[i] = 2.0 * tr_xxx_yzz[i] - 2.0 * tr_xxx_yyyzz[i] * tke_0 + 4.0 * tr_xxxz_yz[i] - 4.0 * tr_xxxz_yzzz[i] * tke_0 - 4.0 * tr_xxxz_yyyz[i] * tke_0 + 4.0 * tr_xxxz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_yzz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyzz[i] * tbe_0 - 4.0 * tr_xxxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_yzzz[i] = tr_xxx_zzz[i] - 2.0 * tr_xxx_yyzzz[i] * tke_0 + 3.0 * tr_xxxz_zz[i] - 2.0 * tr_xxxz_zzzz[i] * tke_0 - 6.0 * tr_xxxz_yyzz[i] * tke_0 + 4.0 * tr_xxxz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_zzz[i] * tbe_0 + 4.0 * tr_xxxzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yzzz[i] * tbe_0 - 6.0 * tr_xxxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_zzzz[i] = -2.0 * tr_xxx_yzzzz[i] * tke_0 - 8.0 * tr_xxxz_yzzz[i] * tke_0 + 4.0 * tr_xxxz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_zzzz[i] * tbe_0 - 8.0 * tr_xxxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 945-960 components of targeted buffer : GG

    auto tr_0_0_yz_xxyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 945);

    auto tr_0_0_yz_xxyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 946);

    auto tr_0_0_yz_xxyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 947);

    auto tr_0_0_yz_xxyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 948);

    auto tr_0_0_yz_xxyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 949);

    auto tr_0_0_yz_xxyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 950);

    auto tr_0_0_yz_xxyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 951);

    auto tr_0_0_yz_xxyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 952);

    auto tr_0_0_yz_xxyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 953);

    auto tr_0_0_yz_xxyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 954);

    auto tr_0_0_yz_xxyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 955);

    auto tr_0_0_yz_xxyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 956);

    auto tr_0_0_yz_xxyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 957);

    auto tr_0_0_yz_xxyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 958);

    auto tr_0_0_yz_xxyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 959);

    #pragma omp simd aligned(tr_0_0_yz_xxyy_xxxx, tr_0_0_yz_xxyy_xxxy, tr_0_0_yz_xxyy_xxxz, tr_0_0_yz_xxyy_xxyy, tr_0_0_yz_xxyy_xxyz, tr_0_0_yz_xxyy_xxzz, tr_0_0_yz_xxyy_xyyy, tr_0_0_yz_xxyy_xyyz, tr_0_0_yz_xxyy_xyzz, tr_0_0_yz_xxyy_xzzz, tr_0_0_yz_xxyy_yyyy, tr_0_0_yz_xxyy_yyyz, tr_0_0_yz_xxyy_yyzz, tr_0_0_yz_xxyy_yzzz, tr_0_0_yz_xxyy_zzzz, tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyy_xx, tr_xxyy_xxxxyz, tr_xxyy_xxxy, tr_xxyy_xxxyyz, tr_xxyy_xxxyzz, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyyyz, tr_xxyy_xxyyzz, tr_xxyy_xxyz, tr_xxyy_xxyzzz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyyyz, tr_xxyy_xyyyzz, tr_xxyy_xyyz, tr_xxyy_xyyzzz, tr_xxyy_xyzz, tr_xxyy_xyzzzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyyyz, tr_xxyy_yyyyzz, tr_xxyy_yyyz, tr_xxyy_yyyzzz, tr_xxyy_yyzz, tr_xxyy_yyzzzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_yzzzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyy_xxx, tr_xxyyy_xxxxz, tr_xxyyy_xxxyz, tr_xxyyy_xxxzz, tr_xxyyy_xxy, tr_xxyyy_xxyyz, tr_xxyyy_xxyzz, tr_xxyyy_xxz, tr_xxyyy_xxzzz, tr_xxyyy_xyy, tr_xxyyy_xyyyz, tr_xxyyy_xyyzz, tr_xxyyy_xyz, tr_xxyyy_xyzzz, tr_xxyyy_xzz, tr_xxyyy_xzzzz, tr_xxyyy_yyy, tr_xxyyy_yyyyz, tr_xxyyy_yyyzz, tr_xxyyy_yyz, tr_xxyyy_yyzzz, tr_xxyyy_yzz, tr_xxyyy_yzzzz, tr_xxyyy_zzz, tr_xxyyy_zzzzz, tr_xxyyyz_xxxx, tr_xxyyyz_xxxy, tr_xxyyyz_xxxz, tr_xxyyyz_xxyy, tr_xxyyyz_xxyz, tr_xxyyyz_xxzz, tr_xxyyyz_xyyy, tr_xxyyyz_xyyz, tr_xxyyyz_xyzz, tr_xxyyyz_xzzz, tr_xxyyyz_yyyy, tr_xxyyyz_yyyz, tr_xxyyyz_yyzz, tr_xxyyyz_yzzz, tr_xxyyyz_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxy, tr_xxyyz_xxxyy, tr_xxyyz_xxxyz, tr_xxyyz_xxy, tr_xxyyz_xxyyy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyyyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyyyy, tr_xxyyz_yyyyz, tr_xxyyz_yyyzz, tr_xxyyz_yyz, tr_xxyyz_yyzzz, tr_xxyyz_yzz, tr_xxyyz_yzzzz, tr_xxyyz_zzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxyy_xxxx[i] = -4.0 * tr_xxy_xxxxz[i] * tke_0 - 4.0 * tr_xxyz_xxxx[i] * tbe_0 + 4.0 * tr_xxyy_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xxxy[i] = -4.0 * tr_xxy_xxxyz[i] * tke_0 - 4.0 * tr_xxyz_xxxy[i] * tbe_0 - 2.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xxxz[i] = 2.0 * tr_xxy_xxx[i] - 4.0 * tr_xxy_xxxzz[i] * tke_0 - 4.0 * tr_xxyz_xxxz[i] * tbe_0 - 2.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xxyy[i] = -4.0 * tr_xxy_xxyyz[i] * tke_0 - 4.0 * tr_xxyz_xxyy[i] * tbe_0 - 4.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xxyz[i] = 2.0 * tr_xxy_xxy[i] - 4.0 * tr_xxy_xxyzz[i] * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tbe_0 + tr_xxyy_xx[i] - 2.0 * tr_xxyy_xxzz[i] * tke_0 - 2.0 * tr_xxyy_xxyy[i] * tke_0 + 4.0 * tr_xxyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xxzz[i] = 4.0 * tr_xxy_xxz[i] - 4.0 * tr_xxy_xxzzz[i] * tke_0 - 4.0 * tr_xxyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xyyy[i] = -4.0 * tr_xxy_xyyyz[i] * tke_0 - 4.0 * tr_xxyz_xyyy[i] * tbe_0 - 6.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xyyz[i] = 2.0 * tr_xxy_xyy[i] - 4.0 * tr_xxy_xyyzz[i] * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxyy_xy[i] - 4.0 * tr_xxyy_xyzz[i] * tke_0 - 2.0 * tr_xxyy_xyyy[i] * tke_0 + 4.0 * tr_xxyy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xyzz[i] = 4.0 * tr_xxy_xyz[i] - 4.0 * tr_xxy_xyzzz[i] * tke_0 - 4.0 * tr_xxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxyy_xz[i] - 2.0 * tr_xxyy_xzzz[i] * tke_0 - 4.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xzzz[i] = 6.0 * tr_xxy_xzz[i] - 4.0 * tr_xxy_xzzzz[i] * tke_0 - 4.0 * tr_xxyz_xzzz[i] * tbe_0 - 6.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_yyyy[i] = -4.0 * tr_xxy_yyyyz[i] * tke_0 - 4.0 * tr_xxyz_yyyy[i] * tbe_0 - 8.0 * tr_xxyy_yyyz[i] * tke_0 + 4.0 * tr_xxyy_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_yyyz[i] = 2.0 * tr_xxy_yyy[i] - 4.0 * tr_xxy_yyyzz[i] * tke_0 - 4.0 * tr_xxyz_yyyz[i] * tbe_0 + 3.0 * tr_xxyy_yy[i] - 6.0 * tr_xxyy_yyzz[i] * tke_0 - 2.0 * tr_xxyy_yyyy[i] * tke_0 + 4.0 * tr_xxyy_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xxyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_yyzz[i] = 4.0 * tr_xxy_yyz[i] - 4.0 * tr_xxy_yyzzz[i] * tke_0 - 4.0 * tr_xxyz_yyzz[i] * tbe_0 + 4.0 * tr_xxyy_yz[i] - 4.0 * tr_xxyy_yzzz[i] * tke_0 - 4.0 * tr_xxyy_yyyz[i] * tke_0 + 4.0 * tr_xxyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_yzzz[i] = 6.0 * tr_xxy_yzz[i] - 4.0 * tr_xxy_yzzzz[i] * tke_0 - 4.0 * tr_xxyz_yzzz[i] * tbe_0 + 3.0 * tr_xxyy_zz[i] - 2.0 * tr_xxyy_zzzz[i] * tke_0 - 6.0 * tr_xxyy_yyzz[i] * tke_0 + 4.0 * tr_xxyy_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_zzzz[i] = 8.0 * tr_xxy_zzz[i] - 4.0 * tr_xxy_zzzzz[i] * tke_0 - 4.0 * tr_xxyz_zzzz[i] * tbe_0 - 8.0 * tr_xxyy_yzzz[i] * tke_0 + 4.0 * tr_xxyy_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xxyyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 960-975 components of targeted buffer : GG

    auto tr_0_0_yz_xxyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 960);

    auto tr_0_0_yz_xxyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 961);

    auto tr_0_0_yz_xxyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 962);

    auto tr_0_0_yz_xxyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 963);

    auto tr_0_0_yz_xxyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 964);

    auto tr_0_0_yz_xxyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 965);

    auto tr_0_0_yz_xxyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 966);

    auto tr_0_0_yz_xxyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 967);

    auto tr_0_0_yz_xxyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 968);

    auto tr_0_0_yz_xxyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 969);

    auto tr_0_0_yz_xxyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 970);

    auto tr_0_0_yz_xxyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 971);

    auto tr_0_0_yz_xxyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 972);

    auto tr_0_0_yz_xxyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 973);

    auto tr_0_0_yz_xxyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 974);

    #pragma omp simd aligned(tr_0_0_yz_xxyz_xxxx, tr_0_0_yz_xxyz_xxxy, tr_0_0_yz_xxyz_xxxz, tr_0_0_yz_xxyz_xxyy, tr_0_0_yz_xxyz_xxyz, tr_0_0_yz_xxyz_xxzz, tr_0_0_yz_xxyz_xyyy, tr_0_0_yz_xxyz_xyyz, tr_0_0_yz_xxyz_xyzz, tr_0_0_yz_xxyz_xzzz, tr_0_0_yz_xxyz_yyyy, tr_0_0_yz_xxyz_yyyz, tr_0_0_yz_xxyz_yyzz, tr_0_0_yz_xxyz_yzzz, tr_0_0_yz_xxyz_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxz, tr_xxyyz_xxxyz, tr_xxyyz_xxxzz, tr_xxyyz_xxy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xxzzz, tr_xxyyz_xyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_xzzzz, tr_xxyyz_yyy, tr_xxyyz_yyyyz, tr_xxyyz_yyyzz, tr_xxyyz_yyz, tr_xxyyz_yyzzz, tr_xxyyz_yzz, tr_xxyyz_yzzzz, tr_xxyyz_zzz, tr_xxyyz_zzzzz, tr_xxyyzz_xxxx, tr_xxyyzz_xxxy, tr_xxyyzz_xxxz, tr_xxyyzz_xxyy, tr_xxyyzz_xxyz, tr_xxyyzz_xxzz, tr_xxyyzz_xyyy, tr_xxyyzz_xyyz, tr_xxyyzz_xyzz, tr_xxyyzz_xzzz, tr_xxyyzz_yyyy, tr_xxyyzz_yyyz, tr_xxyyzz_yyzz, tr_xxyyzz_yzzz, tr_xxyyzz_zzzz, tr_xxyz_xx, tr_xxyz_xxxxyz, tr_xxyz_xxxy, tr_xxyz_xxxyyz, tr_xxyz_xxxyzz, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyyyz, tr_xxyz_xxyyzz, tr_xxyz_xxyz, tr_xxyz_xxyzzz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyyyz, tr_xxyz_xyyyzz, tr_xxyz_xyyz, tr_xxyz_xyyzzz, tr_xxyz_xyzz, tr_xxyz_xyzzzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyyyz, tr_xxyz_yyyyzz, tr_xxyz_yyyz, tr_xxyz_yyyzzz, tr_xxyz_yyzz, tr_xxyz_yyzzzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_yzzzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxy, tr_xxyzz_xxxyy, tr_xxyzz_xxxyz, tr_xxyzz_xxy, tr_xxyzz_xxyyy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyyyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyyyy, tr_xxyzz_yyyyz, tr_xxyzz_yyyzz, tr_xxyzz_yyz, tr_xxyzz_yyzzz, tr_xxyzz_yzz, tr_xxyzz_yzzzz, tr_xxyzz_zzz, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxyz_xxxx[i] = tr_xx_xxxx[i] - 2.0 * tr_xxz_xxxxz[i] * tke_0 - 2.0 * tr_xxzz_xxxx[i] * tbe_0 - 2.0 * tr_xxy_xxxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xxxy[i] = tr_xx_xxxy[i] - 2.0 * tr_xxz_xxxyz[i] * tke_0 - 2.0 * tr_xxzz_xxxy[i] * tbe_0 + tr_xxy_xxx[i] - 2.0 * tr_xxy_xxxyy[i] * tke_0 - 2.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xxxz[i] = tr_xx_xxxz[i] + tr_xxz_xxx[i] - 2.0 * tr_xxz_xxxzz[i] * tke_0 - 2.0 * tr_xxzz_xxxz[i] * tbe_0 - 2.0 * tr_xxy_xxxyz[i] * tke_0 - 2.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxxz[i] * tbe_0 - 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xxyy[i] = tr_xx_xxyy[i] - 2.0 * tr_xxz_xxyyz[i] * tke_0 - 2.0 * tr_xxzz_xxyy[i] * tbe_0 + 2.0 * tr_xxy_xxy[i] - 2.0 * tr_xxy_xxyyy[i] * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxyy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xxyz[i] = tr_xx_xxyz[i] + tr_xxz_xxy[i] - 2.0 * tr_xxz_xxyzz[i] * tke_0 - 2.0 * tr_xxzz_xxyz[i] * tbe_0 + tr_xxy_xxz[i] - 2.0 * tr_xxy_xxyyz[i] * tke_0 + tr_xxyz_xx[i] - 2.0 * tr_xxyz_xxzz[i] * tke_0 - 2.0 * tr_xxyz_xxyy[i] * tke_0 + 4.0 * tr_xxyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxyz[i] * tbe_0 - 2.0 * tr_xxyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xxzz[i] = tr_xx_xxzz[i] + 2.0 * tr_xxz_xxz[i] - 2.0 * tr_xxz_xxzzz[i] * tke_0 - 2.0 * tr_xxzz_xxzz[i] * tbe_0 - 2.0 * tr_xxy_xxyzz[i] * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxzz[i] * tbe_0 - 4.0 * tr_xxyyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xyyy[i] = tr_xx_xyyy[i] - 2.0 * tr_xxz_xyyyz[i] * tke_0 - 2.0 * tr_xxzz_xyyy[i] * tbe_0 + 3.0 * tr_xxy_xyy[i] - 2.0 * tr_xxy_xyyyy[i] * tke_0 - 6.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxyzz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xyyz[i] = tr_xx_xyyz[i] + tr_xxz_xyy[i] - 2.0 * tr_xxz_xyyzz[i] * tke_0 - 2.0 * tr_xxzz_xyyz[i] * tbe_0 + 2.0 * tr_xxy_xyz[i] - 2.0 * tr_xxy_xyyyz[i] * tke_0 + 2.0 * tr_xxyz_xy[i] - 4.0 * tr_xxyz_xyzz[i] * tke_0 - 2.0 * tr_xxyz_xyyy[i] * tke_0 + 4.0 * tr_xxyz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyyz[i] * tbe_0 - 2.0 * tr_xxyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xyzz[i] = tr_xx_xyzz[i] + 2.0 * tr_xxz_xyz[i] - 2.0 * tr_xxz_xyzzz[i] * tke_0 - 2.0 * tr_xxzz_xyzz[i] * tbe_0 + tr_xxy_xzz[i] - 2.0 * tr_xxy_xyyzz[i] * tke_0 + 2.0 * tr_xxyz_xz[i] - 2.0 * tr_xxyz_xzzz[i] * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyzz[i] * tbe_0 - 4.0 * tr_xxyyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xzzz[i] = tr_xx_xzzz[i] + 3.0 * tr_xxz_xzz[i] - 2.0 * tr_xxz_xzzzz[i] * tke_0 - 2.0 * tr_xxzz_xzzz[i] * tbe_0 - 2.0 * tr_xxy_xyzzz[i] * tke_0 - 6.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xzzz[i] * tbe_0 - 6.0 * tr_xxyyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_yyyy[i] = tr_xx_yyyy[i] - 2.0 * tr_xxz_yyyyz[i] * tke_0 - 2.0 * tr_xxzz_yyyy[i] * tbe_0 + 4.0 * tr_xxy_yyy[i] - 2.0 * tr_xxy_yyyyy[i] * tke_0 - 8.0 * tr_xxyz_yyyz[i] * tke_0 + 4.0 * tr_xxyz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_yyyz[i] = tr_xx_yyyz[i] + tr_xxz_yyy[i] - 2.0 * tr_xxz_yyyzz[i] * tke_0 - 2.0 * tr_xxzz_yyyz[i] * tbe_0 + 3.0 * tr_xxy_yyz[i] - 2.0 * tr_xxy_yyyyz[i] * tke_0 + 3.0 * tr_xxyz_yy[i] - 6.0 * tr_xxyz_yyzz[i] * tke_0 - 2.0 * tr_xxyz_yyyy[i] * tke_0 + 4.0 * tr_xxyz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xxyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyyz[i] * tbe_0 - 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_yyzz[i] = tr_xx_yyzz[i] + 2.0 * tr_xxz_yyz[i] - 2.0 * tr_xxz_yyzzz[i] * tke_0 - 2.0 * tr_xxzz_yyzz[i] * tbe_0 + 2.0 * tr_xxy_yzz[i] - 2.0 * tr_xxy_yyyzz[i] * tke_0 + 4.0 * tr_xxyz_yz[i] - 4.0 * tr_xxyz_yzzz[i] * tke_0 - 4.0 * tr_xxyz_yyyz[i] * tke_0 + 4.0 * tr_xxyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyzz[i] * tbe_0 - 4.0 * tr_xxyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_yzzz[i] = tr_xx_yzzz[i] + 3.0 * tr_xxz_yzz[i] - 2.0 * tr_xxz_yzzzz[i] * tke_0 - 2.0 * tr_xxzz_yzzz[i] * tbe_0 + tr_xxy_zzz[i] - 2.0 * tr_xxy_yyzzz[i] * tke_0 + 3.0 * tr_xxyz_zz[i] - 2.0 * tr_xxyz_zzzz[i] * tke_0 - 6.0 * tr_xxyz_yyzz[i] * tke_0 + 4.0 * tr_xxyz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yzzz[i] * tbe_0 - 6.0 * tr_xxyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_zzzz[i] = tr_xx_zzzz[i] + 4.0 * tr_xxz_zzz[i] - 2.0 * tr_xxz_zzzzz[i] * tke_0 - 2.0 * tr_xxzz_zzzz[i] * tbe_0 - 2.0 * tr_xxy_yzzzz[i] * tke_0 - 8.0 * tr_xxyz_yzzz[i] * tke_0 + 4.0 * tr_xxyz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_zzzz[i] * tbe_0 - 8.0 * tr_xxyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 975-990 components of targeted buffer : GG

    auto tr_0_0_yz_xxzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 975);

    auto tr_0_0_yz_xxzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 976);

    auto tr_0_0_yz_xxzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 977);

    auto tr_0_0_yz_xxzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 978);

    auto tr_0_0_yz_xxzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 979);

    auto tr_0_0_yz_xxzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 980);

    auto tr_0_0_yz_xxzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 981);

    auto tr_0_0_yz_xxzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 982);

    auto tr_0_0_yz_xxzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 983);

    auto tr_0_0_yz_xxzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 984);

    auto tr_0_0_yz_xxzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 985);

    auto tr_0_0_yz_xxzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 986);

    auto tr_0_0_yz_xxzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 987);

    auto tr_0_0_yz_xxzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 988);

    auto tr_0_0_yz_xxzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 989);

    #pragma omp simd aligned(tr_0_0_yz_xxzz_xxxx, tr_0_0_yz_xxzz_xxxy, tr_0_0_yz_xxzz_xxxz, tr_0_0_yz_xxzz_xxyy, tr_0_0_yz_xxzz_xxyz, tr_0_0_yz_xxzz_xxzz, tr_0_0_yz_xxzz_xyyy, tr_0_0_yz_xxzz_xyyz, tr_0_0_yz_xxzz_xyzz, tr_0_0_yz_xxzz_xzzz, tr_0_0_yz_xxzz_yyyy, tr_0_0_yz_xxzz_yyyz, tr_0_0_yz_xxzz_yyzz, tr_0_0_yz_xxzz_yzzz, tr_0_0_yz_xxzz_zzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxz, tr_xxyzz_xxxyz, tr_xxyzz_xxxzz, tr_xxyzz_xxy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xxzzz, tr_xxyzz_xyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_xzzzz, tr_xxyzz_yyy, tr_xxyzz_yyyyz, tr_xxyzz_yyyzz, tr_xxyzz_yyz, tr_xxyzz_yyzzz, tr_xxyzz_yzz, tr_xxyzz_yzzzz, tr_xxyzz_zzz, tr_xxyzz_zzzzz, tr_xxyzzz_xxxx, tr_xxyzzz_xxxy, tr_xxyzzz_xxxz, tr_xxyzzz_xxyy, tr_xxyzzz_xxyz, tr_xxyzzz_xxzz, tr_xxyzzz_xyyy, tr_xxyzzz_xyyz, tr_xxyzzz_xyzz, tr_xxyzzz_xzzz, tr_xxyzzz_yyyy, tr_xxyzzz_yyyz, tr_xxyzzz_yyzz, tr_xxyzzz_yzzz, tr_xxyzzz_zzzz, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xxxxyz, tr_xxzz_xxxy, tr_xxzz_xxxyyz, tr_xxzz_xxxyzz, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyyyz, tr_xxzz_xxyyzz, tr_xxzz_xxyz, tr_xxzz_xxyzzz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyyyz, tr_xxzz_xyyyzz, tr_xxzz_xyyz, tr_xxzz_xyyzzz, tr_xxzz_xyzz, tr_xxzz_xyzzzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyyyz, tr_xxzz_yyyyzz, tr_xxzz_yyyz, tr_xxzz_yyyzzz, tr_xxzz_yyzz, tr_xxzz_yyzzzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_yzzzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxxxy, tr_xxzzz_xxxyy, tr_xxzzz_xxxyz, tr_xxzzz_xxy, tr_xxzzz_xxyyy, tr_xxzzz_xxyyz, tr_xxzzz_xxyzz, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyyyy, tr_xxzzz_xyyyz, tr_xxzzz_xyyzz, tr_xxzzz_xyz, tr_xxzzz_xyzzz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyyyy, tr_xxzzz_yyyyz, tr_xxzzz_yyyzz, tr_xxzzz_yyz, tr_xxzzz_yyzzz, tr_xxzzz_yzz, tr_xxzzz_yzzzz, tr_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxzz_xxxx[i] = -4.0 * tr_xxz_xxxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xxxy[i] = 2.0 * tr_xxz_xxx[i] - 4.0 * tr_xxz_xxxyy[i] * tke_0 - 2.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xxxz[i] = -4.0 * tr_xxz_xxxyz[i] * tke_0 - 2.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxxz[i] * tbe_0 - 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xxyy[i] = 4.0 * tr_xxz_xxy[i] - 4.0 * tr_xxz_xxyyy[i] * tke_0 - 4.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxyy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xxyz[i] = 2.0 * tr_xxz_xxz[i] - 4.0 * tr_xxz_xxyyz[i] * tke_0 + tr_xxzz_xx[i] - 2.0 * tr_xxzz_xxzz[i] * tke_0 - 2.0 * tr_xxzz_xxyy[i] * tke_0 + 4.0 * tr_xxzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_xxz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxyz[i] * tbe_0 - 2.0 * tr_xxyzz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xxzz[i] = -4.0 * tr_xxz_xxyzz[i] * tke_0 - 4.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxzz[i] * tbe_0 - 4.0 * tr_xxyzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xyyy[i] = 6.0 * tr_xxz_xyy[i] - 4.0 * tr_xxz_xyyyy[i] * tke_0 - 6.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xxzzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xyyz[i] = 4.0 * tr_xxz_xyz[i] - 4.0 * tr_xxz_xyyyz[i] * tke_0 + 2.0 * tr_xxzz_xy[i] - 4.0 * tr_xxzz_xyzz[i] * tke_0 - 2.0 * tr_xxzz_xyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xyz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyyz[i] * tbe_0 - 2.0 * tr_xxyzz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xyzz[i] = 2.0 * tr_xxz_xzz[i] - 4.0 * tr_xxz_xyyzz[i] * tke_0 + 2.0 * tr_xxzz_xz[i] - 2.0 * tr_xxzz_xzzz[i] * tke_0 - 4.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_xzz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyzz[i] * tbe_0 - 4.0 * tr_xxyzz_xyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xzzz[i] = -4.0 * tr_xxz_xyzzz[i] * tke_0 - 6.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xzzz[i] * tbe_0 - 6.0 * tr_xxyzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_yyyy[i] = 8.0 * tr_xxz_yyy[i] - 4.0 * tr_xxz_yyyyy[i] * tke_0 - 8.0 * tr_xxzz_yyyz[i] * tke_0 + 4.0 * tr_xxzz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_yyyz[i] = 6.0 * tr_xxz_yyz[i] - 4.0 * tr_xxz_yyyyz[i] * tke_0 + 3.0 * tr_xxzz_yy[i] - 6.0 * tr_xxzz_yyzz[i] * tke_0 - 2.0 * tr_xxzz_yyyy[i] * tke_0 + 4.0 * tr_xxzz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xxzzz_yyz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyyz[i] * tbe_0 - 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_yyzz[i] = 4.0 * tr_xxz_yzz[i] - 4.0 * tr_xxz_yyyzz[i] * tke_0 + 4.0 * tr_xxzz_yz[i] - 4.0 * tr_xxzz_yzzz[i] * tke_0 - 4.0 * tr_xxzz_yyyz[i] * tke_0 + 4.0 * tr_xxzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_yzz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyzz[i] * tbe_0 - 4.0 * tr_xxyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_yzzz[i] = 2.0 * tr_xxz_zzz[i] - 4.0 * tr_xxz_yyzzz[i] * tke_0 + 3.0 * tr_xxzz_zz[i] - 2.0 * tr_xxzz_zzzz[i] * tke_0 - 6.0 * tr_xxzz_yyzz[i] * tke_0 + 4.0 * tr_xxzz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzz_yyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yzzz[i] * tbe_0 - 6.0 * tr_xxyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_zzzz[i] = -4.0 * tr_xxz_yzzzz[i] * tke_0 - 8.0 * tr_xxzz_yzzz[i] * tke_0 + 4.0 * tr_xxzz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_yzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zzzz[i] * tbe_0 - 8.0 * tr_xxyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 990-1005 components of targeted buffer : GG

    auto tr_0_0_yz_xyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 990);

    auto tr_0_0_yz_xyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 991);

    auto tr_0_0_yz_xyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 992);

    auto tr_0_0_yz_xyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 993);

    auto tr_0_0_yz_xyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 994);

    auto tr_0_0_yz_xyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 995);

    auto tr_0_0_yz_xyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 996);

    auto tr_0_0_yz_xyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 997);

    auto tr_0_0_yz_xyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 998);

    auto tr_0_0_yz_xyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 999);

    auto tr_0_0_yz_xyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 1000);

    auto tr_0_0_yz_xyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 1001);

    auto tr_0_0_yz_xyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 1002);

    auto tr_0_0_yz_xyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 1003);

    auto tr_0_0_yz_xyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 1004);

    #pragma omp simd aligned(tr_0_0_yz_xyyy_xxxx, tr_0_0_yz_xyyy_xxxy, tr_0_0_yz_xyyy_xxxz, tr_0_0_yz_xyyy_xxyy, tr_0_0_yz_xyyy_xxyz, tr_0_0_yz_xyyy_xxzz, tr_0_0_yz_xyyy_xyyy, tr_0_0_yz_xyyy_xyyz, tr_0_0_yz_xyyy_xyzz, tr_0_0_yz_xyyy_xzzz, tr_0_0_yz_xyyy_yyyy, tr_0_0_yz_xyyy_yyyz, tr_0_0_yz_xyyy_yyzz, tr_0_0_yz_xyyy_yzzz, tr_0_0_yz_xyyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyy_xx, tr_xyyy_xxxxyz, tr_xyyy_xxxy, tr_xyyy_xxxyyz, tr_xyyy_xxxyzz, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyyyz, tr_xyyy_xxyyzz, tr_xyyy_xxyz, tr_xyyy_xxyzzz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyyyz, tr_xyyy_xyyyzz, tr_xyyy_xyyz, tr_xyyy_xyyzzz, tr_xyyy_xyzz, tr_xyyy_xyzzzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyyyz, tr_xyyy_yyyyzz, tr_xyyy_yyyz, tr_xyyy_yyyzzz, tr_xyyy_yyzz, tr_xyyy_yyzzzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_yzzzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyy_xxx, tr_xyyyy_xxxxz, tr_xyyyy_xxxyz, tr_xyyyy_xxxzz, tr_xyyyy_xxy, tr_xyyyy_xxyyz, tr_xyyyy_xxyzz, tr_xyyyy_xxz, tr_xyyyy_xxzzz, tr_xyyyy_xyy, tr_xyyyy_xyyyz, tr_xyyyy_xyyzz, tr_xyyyy_xyz, tr_xyyyy_xyzzz, tr_xyyyy_xzz, tr_xyyyy_xzzzz, tr_xyyyy_yyy, tr_xyyyy_yyyyz, tr_xyyyy_yyyzz, tr_xyyyy_yyz, tr_xyyyy_yyzzz, tr_xyyyy_yzz, tr_xyyyy_yzzzz, tr_xyyyy_zzz, tr_xyyyy_zzzzz, tr_xyyyyz_xxxx, tr_xyyyyz_xxxy, tr_xyyyyz_xxxz, tr_xyyyyz_xxyy, tr_xyyyyz_xxyz, tr_xyyyyz_xxzz, tr_xyyyyz_xyyy, tr_xyyyyz_xyyz, tr_xyyyyz_xyzz, tr_xyyyyz_xzzz, tr_xyyyyz_yyyy, tr_xyyyyz_yyyz, tr_xyyyyz_yyzz, tr_xyyyyz_yzzz, tr_xyyyyz_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxy, tr_xyyyz_xxxyy, tr_xyyyz_xxxyz, tr_xyyyz_xxy, tr_xyyyz_xxyyy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyyyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyyyy, tr_xyyyz_yyyyz, tr_xyyyz_yyyzz, tr_xyyyz_yyz, tr_xyyyz_yyzzz, tr_xyyyz_yzz, tr_xyyyz_yzzzz, tr_xyyyz_zzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyyy_xxxx[i] = -6.0 * tr_xyy_xxxxz[i] * tke_0 - 6.0 * tr_xyyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyy_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xxxy[i] = -6.0 * tr_xyy_xxxyz[i] * tke_0 - 6.0 * tr_xyyz_xxxy[i] * tbe_0 - 2.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xxxz[i] = 3.0 * tr_xyy_xxx[i] - 6.0 * tr_xyy_xxxzz[i] * tke_0 - 6.0 * tr_xyyz_xxxz[i] * tbe_0 - 2.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xxyy[i] = -6.0 * tr_xyy_xxyyz[i] * tke_0 - 6.0 * tr_xyyz_xxyy[i] * tbe_0 - 4.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xxyz[i] = 3.0 * tr_xyy_xxy[i] - 6.0 * tr_xyy_xxyzz[i] * tke_0 - 6.0 * tr_xyyz_xxyz[i] * tbe_0 + tr_xyyy_xx[i] - 2.0 * tr_xyyy_xxzz[i] * tke_0 - 2.0 * tr_xyyy_xxyy[i] * tke_0 + 4.0 * tr_xyyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xxzz[i] = 6.0 * tr_xyy_xxz[i] - 6.0 * tr_xyy_xxzzz[i] * tke_0 - 6.0 * tr_xyyz_xxzz[i] * tbe_0 - 4.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xyyy[i] = -6.0 * tr_xyy_xyyyz[i] * tke_0 - 6.0 * tr_xyyz_xyyy[i] * tbe_0 - 6.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xyyz[i] = 3.0 * tr_xyy_xyy[i] - 6.0 * tr_xyy_xyyzz[i] * tke_0 - 6.0 * tr_xyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyyy_xy[i] - 4.0 * tr_xyyy_xyzz[i] * tke_0 - 2.0 * tr_xyyy_xyyy[i] * tke_0 + 4.0 * tr_xyyy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xyzz[i] = 6.0 * tr_xyy_xyz[i] - 6.0 * tr_xyy_xyzzz[i] * tke_0 - 6.0 * tr_xyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyyy_xz[i] - 2.0 * tr_xyyy_xzzz[i] * tke_0 - 4.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xzzz[i] = 9.0 * tr_xyy_xzz[i] - 6.0 * tr_xyy_xzzzz[i] * tke_0 - 6.0 * tr_xyyz_xzzz[i] * tbe_0 - 6.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_yyyy[i] = -6.0 * tr_xyy_yyyyz[i] * tke_0 - 6.0 * tr_xyyz_yyyy[i] * tbe_0 - 8.0 * tr_xyyy_yyyz[i] * tke_0 + 4.0 * tr_xyyy_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_yyyz[i] = 3.0 * tr_xyy_yyy[i] - 6.0 * tr_xyy_yyyzz[i] * tke_0 - 6.0 * tr_xyyz_yyyz[i] * tbe_0 + 3.0 * tr_xyyy_yy[i] - 6.0 * tr_xyyy_yyzz[i] * tke_0 - 2.0 * tr_xyyy_yyyy[i] * tke_0 + 4.0 * tr_xyyy_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xyyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_yyzz[i] = 6.0 * tr_xyy_yyz[i] - 6.0 * tr_xyy_yyzzz[i] * tke_0 - 6.0 * tr_xyyz_yyzz[i] * tbe_0 + 4.0 * tr_xyyy_yz[i] - 4.0 * tr_xyyy_yzzz[i] * tke_0 - 4.0 * tr_xyyy_yyyz[i] * tke_0 + 4.0 * tr_xyyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_yzzz[i] = 9.0 * tr_xyy_yzz[i] - 6.0 * tr_xyy_yzzzz[i] * tke_0 - 6.0 * tr_xyyz_yzzz[i] * tbe_0 + 3.0 * tr_xyyy_zz[i] - 2.0 * tr_xyyy_zzzz[i] * tke_0 - 6.0 * tr_xyyy_yyzz[i] * tke_0 + 4.0 * tr_xyyy_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_zzzz[i] = 12.0 * tr_xyy_zzz[i] - 6.0 * tr_xyy_zzzzz[i] * tke_0 - 6.0 * tr_xyyz_zzzz[i] * tbe_0 - 8.0 * tr_xyyy_yzzz[i] * tke_0 + 4.0 * tr_xyyy_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1005-1020 components of targeted buffer : GG

    auto tr_0_0_yz_xyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1005);

    auto tr_0_0_yz_xyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1006);

    auto tr_0_0_yz_xyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1007);

    auto tr_0_0_yz_xyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1008);

    auto tr_0_0_yz_xyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1009);

    auto tr_0_0_yz_xyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1010);

    auto tr_0_0_yz_xyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1011);

    auto tr_0_0_yz_xyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1012);

    auto tr_0_0_yz_xyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1013);

    auto tr_0_0_yz_xyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1014);

    auto tr_0_0_yz_xyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1015);

    auto tr_0_0_yz_xyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1016);

    auto tr_0_0_yz_xyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1017);

    auto tr_0_0_yz_xyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1018);

    auto tr_0_0_yz_xyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1019);

    #pragma omp simd aligned(tr_0_0_yz_xyyz_xxxx, tr_0_0_yz_xyyz_xxxy, tr_0_0_yz_xyyz_xxxz, tr_0_0_yz_xyyz_xxyy, tr_0_0_yz_xyyz_xxyz, tr_0_0_yz_xyyz_xxzz, tr_0_0_yz_xyyz_xyyy, tr_0_0_yz_xyyz_xyyz, tr_0_0_yz_xyyz_xyzz, tr_0_0_yz_xyyz_xzzz, tr_0_0_yz_xyyz_yyyy, tr_0_0_yz_xyyz_yyyz, tr_0_0_yz_xyyz_yyzz, tr_0_0_yz_xyyz_yzzz, tr_0_0_yz_xyyz_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxz, tr_xyyyz_xxxyz, tr_xyyyz_xxxzz, tr_xyyyz_xxy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xxzzz, tr_xyyyz_xyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_xzzzz, tr_xyyyz_yyy, tr_xyyyz_yyyyz, tr_xyyyz_yyyzz, tr_xyyyz_yyz, tr_xyyyz_yyzzz, tr_xyyyz_yzz, tr_xyyyz_yzzzz, tr_xyyyz_zzz, tr_xyyyz_zzzzz, tr_xyyyzz_xxxx, tr_xyyyzz_xxxy, tr_xyyyzz_xxxz, tr_xyyyzz_xxyy, tr_xyyyzz_xxyz, tr_xyyyzz_xxzz, tr_xyyyzz_xyyy, tr_xyyyzz_xyyz, tr_xyyyzz_xyzz, tr_xyyyzz_xzzz, tr_xyyyzz_yyyy, tr_xyyyzz_yyyz, tr_xyyyzz_yyzz, tr_xyyyzz_yzzz, tr_xyyyzz_zzzz, tr_xyyz_xx, tr_xyyz_xxxxyz, tr_xyyz_xxxy, tr_xyyz_xxxyyz, tr_xyyz_xxxyzz, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyyyz, tr_xyyz_xxyyzz, tr_xyyz_xxyz, tr_xyyz_xxyzzz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyyyz, tr_xyyz_xyyyzz, tr_xyyz_xyyz, tr_xyyz_xyyzzz, tr_xyyz_xyzz, tr_xyyz_xyzzzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyyyz, tr_xyyz_yyyyzz, tr_xyyz_yyyz, tr_xyyz_yyyzzz, tr_xyyz_yyzz, tr_xyyz_yyzzzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_yzzzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxy, tr_xyyzz_xxxyy, tr_xyyzz_xxxyz, tr_xyyzz_xxy, tr_xyyzz_xxyyy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyyyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyyyy, tr_xyyzz_yyyyz, tr_xyyzz_yyyzz, tr_xyyzz_yyz, tr_xyyzz_yyzzz, tr_xyyzz_yzz, tr_xyyzz_yzzzz, tr_xyyzz_zzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyyz_xxxx[i] = 2.0 * tr_xy_xxxx[i] - 4.0 * tr_xyz_xxxxz[i] * tke_0 - 4.0 * tr_xyzz_xxxx[i] * tbe_0 - 2.0 * tr_xyy_xxxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xxxy[i] = 2.0 * tr_xy_xxxy[i] - 4.0 * tr_xyz_xxxyz[i] * tke_0 - 4.0 * tr_xyzz_xxxy[i] * tbe_0 + tr_xyy_xxx[i] - 2.0 * tr_xyy_xxxyy[i] * tke_0 - 2.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xxxz[i] = 2.0 * tr_xy_xxxz[i] + 2.0 * tr_xyz_xxx[i] - 4.0 * tr_xyz_xxxzz[i] * tke_0 - 4.0 * tr_xyzz_xxxz[i] * tbe_0 - 2.0 * tr_xyy_xxxyz[i] * tke_0 - 2.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxxz[i] * tbe_0 - 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xxyy[i] = 2.0 * tr_xy_xxyy[i] - 4.0 * tr_xyz_xxyyz[i] * tke_0 - 4.0 * tr_xyzz_xxyy[i] * tbe_0 + 2.0 * tr_xyy_xxy[i] - 2.0 * tr_xyy_xxyyy[i] * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxyy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xxyz[i] = 2.0 * tr_xy_xxyz[i] + 2.0 * tr_xyz_xxy[i] - 4.0 * tr_xyz_xxyzz[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tbe_0 + tr_xyy_xxz[i] - 2.0 * tr_xyy_xxyyz[i] * tke_0 + tr_xyyz_xx[i] - 2.0 * tr_xyyz_xxzz[i] * tke_0 - 2.0 * tr_xyyz_xxyy[i] * tke_0 + 4.0 * tr_xyyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxyz[i] * tbe_0 - 2.0 * tr_xyyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xxzz[i] = 2.0 * tr_xy_xxzz[i] + 4.0 * tr_xyz_xxz[i] - 4.0 * tr_xyz_xxzzz[i] * tke_0 - 4.0 * tr_xyzz_xxzz[i] * tbe_0 - 2.0 * tr_xyy_xxyzz[i] * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxzz[i] * tbe_0 - 4.0 * tr_xyyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xyyy[i] = 2.0 * tr_xy_xyyy[i] - 4.0 * tr_xyz_xyyyz[i] * tke_0 - 4.0 * tr_xyzz_xyyy[i] * tbe_0 + 3.0 * tr_xyy_xyy[i] - 2.0 * tr_xyy_xyyyy[i] * tke_0 - 6.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xyyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xyyz[i] = 2.0 * tr_xy_xyyz[i] + 2.0 * tr_xyz_xyy[i] - 4.0 * tr_xyz_xyyzz[i] * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tbe_0 + 2.0 * tr_xyy_xyz[i] - 2.0 * tr_xyy_xyyyz[i] * tke_0 + 2.0 * tr_xyyz_xy[i] - 4.0 * tr_xyyz_xyzz[i] * tke_0 - 2.0 * tr_xyyz_xyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyyz[i] * tbe_0 - 2.0 * tr_xyyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xyzz[i] = 2.0 * tr_xy_xyzz[i] + 4.0 * tr_xyz_xyz[i] - 4.0 * tr_xyz_xyzzz[i] * tke_0 - 4.0 * tr_xyzz_xyzz[i] * tbe_0 + tr_xyy_xzz[i] - 2.0 * tr_xyy_xyyzz[i] * tke_0 + 2.0 * tr_xyyz_xz[i] - 2.0 * tr_xyyz_xzzz[i] * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyzz[i] * tbe_0 - 4.0 * tr_xyyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xzzz[i] = 2.0 * tr_xy_xzzz[i] + 6.0 * tr_xyz_xzz[i] - 4.0 * tr_xyz_xzzzz[i] * tke_0 - 4.0 * tr_xyzz_xzzz[i] * tbe_0 - 2.0 * tr_xyy_xyzzz[i] * tke_0 - 6.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xzzz[i] * tbe_0 - 6.0 * tr_xyyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_yyyy[i] = 2.0 * tr_xy_yyyy[i] - 4.0 * tr_xyz_yyyyz[i] * tke_0 - 4.0 * tr_xyzz_yyyy[i] * tbe_0 + 4.0 * tr_xyy_yyy[i] - 2.0 * tr_xyy_yyyyy[i] * tke_0 - 8.0 * tr_xyyz_yyyz[i] * tke_0 + 4.0 * tr_xyyz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_yyyz[i] = 2.0 * tr_xy_yyyz[i] + 2.0 * tr_xyz_yyy[i] - 4.0 * tr_xyz_yyyzz[i] * tke_0 - 4.0 * tr_xyzz_yyyz[i] * tbe_0 + 3.0 * tr_xyy_yyz[i] - 2.0 * tr_xyy_yyyyz[i] * tke_0 + 3.0 * tr_xyyz_yy[i] - 6.0 * tr_xyyz_yyzz[i] * tke_0 - 2.0 * tr_xyyz_yyyy[i] * tke_0 + 4.0 * tr_xyyz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xyyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyyz[i] * tbe_0 - 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_yyzz[i] = 2.0 * tr_xy_yyzz[i] + 4.0 * tr_xyz_yyz[i] - 4.0 * tr_xyz_yyzzz[i] * tke_0 - 4.0 * tr_xyzz_yyzz[i] * tbe_0 + 2.0 * tr_xyy_yzz[i] - 2.0 * tr_xyy_yyyzz[i] * tke_0 + 4.0 * tr_xyyz_yz[i] - 4.0 * tr_xyyz_yzzz[i] * tke_0 - 4.0 * tr_xyyz_yyyz[i] * tke_0 + 4.0 * tr_xyyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyzz[i] * tbe_0 - 4.0 * tr_xyyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_yzzz[i] = 2.0 * tr_xy_yzzz[i] + 6.0 * tr_xyz_yzz[i] - 4.0 * tr_xyz_yzzzz[i] * tke_0 - 4.0 * tr_xyzz_yzzz[i] * tbe_0 + tr_xyy_zzz[i] - 2.0 * tr_xyy_yyzzz[i] * tke_0 + 3.0 * tr_xyyz_zz[i] - 2.0 * tr_xyyz_zzzz[i] * tke_0 - 6.0 * tr_xyyz_yyzz[i] * tke_0 + 4.0 * tr_xyyz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yzzz[i] * tbe_0 - 6.0 * tr_xyyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_zzzz[i] = 2.0 * tr_xy_zzzz[i] + 8.0 * tr_xyz_zzz[i] - 4.0 * tr_xyz_zzzzz[i] * tke_0 - 4.0 * tr_xyzz_zzzz[i] * tbe_0 - 2.0 * tr_xyy_yzzzz[i] * tke_0 - 8.0 * tr_xyyz_yzzz[i] * tke_0 + 4.0 * tr_xyyz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_zzzz[i] * tbe_0 - 8.0 * tr_xyyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1020-1035 components of targeted buffer : GG

    auto tr_0_0_yz_xyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1020);

    auto tr_0_0_yz_xyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1021);

    auto tr_0_0_yz_xyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1022);

    auto tr_0_0_yz_xyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1023);

    auto tr_0_0_yz_xyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1024);

    auto tr_0_0_yz_xyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1025);

    auto tr_0_0_yz_xyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1026);

    auto tr_0_0_yz_xyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1027);

    auto tr_0_0_yz_xyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1028);

    auto tr_0_0_yz_xyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1029);

    auto tr_0_0_yz_xyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1030);

    auto tr_0_0_yz_xyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1031);

    auto tr_0_0_yz_xyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1032);

    auto tr_0_0_yz_xyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1033);

    auto tr_0_0_yz_xyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1034);

    #pragma omp simd aligned(tr_0_0_yz_xyzz_xxxx, tr_0_0_yz_xyzz_xxxy, tr_0_0_yz_xyzz_xxxz, tr_0_0_yz_xyzz_xxyy, tr_0_0_yz_xyzz_xxyz, tr_0_0_yz_xyzz_xxzz, tr_0_0_yz_xyzz_xyyy, tr_0_0_yz_xyzz_xyyz, tr_0_0_yz_xyzz_xyzz, tr_0_0_yz_xyzz_xzzz, tr_0_0_yz_xyzz_yyyy, tr_0_0_yz_xyzz_yyyz, tr_0_0_yz_xyzz_yyzz, tr_0_0_yz_xyzz_yzzz, tr_0_0_yz_xyzz_zzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxz, tr_xyyzz_xxxyz, tr_xyyzz_xxxzz, tr_xyyzz_xxy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xxzzz, tr_xyyzz_xyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_xzzzz, tr_xyyzz_yyy, tr_xyyzz_yyyyz, tr_xyyzz_yyyzz, tr_xyyzz_yyz, tr_xyyzz_yyzzz, tr_xyyzz_yzz, tr_xyyzz_yzzzz, tr_xyyzz_zzz, tr_xyyzz_zzzzz, tr_xyyzzz_xxxx, tr_xyyzzz_xxxy, tr_xyyzzz_xxxz, tr_xyyzzz_xxyy, tr_xyyzzz_xxyz, tr_xyyzzz_xxzz, tr_xyyzzz_xyyy, tr_xyyzzz_xyyz, tr_xyyzzz_xyzz, tr_xyyzzz_xzzz, tr_xyyzzz_yyyy, tr_xyyzzz_yyyz, tr_xyyzzz_yyzz, tr_xyyzzz_yzzz, tr_xyyzzz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxxyz, tr_xyzz_xxxy, tr_xyzz_xxxyyz, tr_xyzz_xxxyzz, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyyyz, tr_xyzz_xxyyzz, tr_xyzz_xxyz, tr_xyzz_xxyzzz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyyyz, tr_xyzz_xyyyzz, tr_xyzz_xyyz, tr_xyzz_xyyzzz, tr_xyzz_xyzz, tr_xyzz_xyzzzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyyyz, tr_xyzz_yyyyzz, tr_xyzz_yyyz, tr_xyzz_yyyzzz, tr_xyzz_yyzz, tr_xyzz_yyzzzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_yzzzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxy, tr_xyzzz_xxxyy, tr_xyzzz_xxxyz, tr_xyzzz_xxy, tr_xyzzz_xxyyy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyyyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyyyy, tr_xyzzz_yyyyz, tr_xyzzz_yyyzz, tr_xyzzz_yyz, tr_xyzzz_yyzzz, tr_xyzzz_yzz, tr_xyzzz_yzzzz, tr_xyzzz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyzz_xxxx[i] = 2.0 * tr_xz_xxxx[i] - 2.0 * tr_xzz_xxxxz[i] * tke_0 - 2.0 * tr_xzzz_xxxx[i] * tbe_0 - 4.0 * tr_xyz_xxxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xxxy[i] = 2.0 * tr_xz_xxxy[i] - 2.0 * tr_xzz_xxxyz[i] * tke_0 - 2.0 * tr_xzzz_xxxy[i] * tbe_0 + 2.0 * tr_xyz_xxx[i] - 4.0 * tr_xyz_xxxyy[i] * tke_0 - 2.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xxxz[i] = 2.0 * tr_xz_xxxz[i] + tr_xzz_xxx[i] - 2.0 * tr_xzz_xxxzz[i] * tke_0 - 2.0 * tr_xzzz_xxxz[i] * tbe_0 - 4.0 * tr_xyz_xxxyz[i] * tke_0 - 2.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxxz[i] * tbe_0 - 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xxyy[i] = 2.0 * tr_xz_xxyy[i] - 2.0 * tr_xzz_xxyyz[i] * tke_0 - 2.0 * tr_xzzz_xxyy[i] * tbe_0 + 4.0 * tr_xyz_xxy[i] - 4.0 * tr_xyz_xxyyy[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxyy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xxyz[i] = 2.0 * tr_xz_xxyz[i] + tr_xzz_xxy[i] - 2.0 * tr_xzz_xxyzz[i] * tke_0 - 2.0 * tr_xzzz_xxyz[i] * tbe_0 + 2.0 * tr_xyz_xxz[i] - 4.0 * tr_xyz_xxyyz[i] * tke_0 + tr_xyzz_xx[i] - 2.0 * tr_xyzz_xxzz[i] * tke_0 - 2.0 * tr_xyzz_xxyy[i] * tke_0 + 4.0 * tr_xyzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxyz[i] * tbe_0 - 2.0 * tr_xyyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xxzz[i] = 2.0 * tr_xz_xxzz[i] + 2.0 * tr_xzz_xxz[i] - 2.0 * tr_xzz_xxzzz[i] * tke_0 - 2.0 * tr_xzzz_xxzz[i] * tbe_0 - 4.0 * tr_xyz_xxyzz[i] * tke_0 - 4.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxzz[i] * tbe_0 - 4.0 * tr_xyyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xyyy[i] = 2.0 * tr_xz_xyyy[i] - 2.0 * tr_xzz_xyyyz[i] * tke_0 - 2.0 * tr_xzzz_xyyy[i] * tbe_0 + 6.0 * tr_xyz_xyy[i] - 4.0 * tr_xyz_xyyyy[i] * tke_0 - 6.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xyzzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xyyz[i] = 2.0 * tr_xz_xyyz[i] + tr_xzz_xyy[i] - 2.0 * tr_xzz_xyyzz[i] * tke_0 - 2.0 * tr_xzzz_xyyz[i] * tbe_0 + 4.0 * tr_xyz_xyz[i] - 4.0 * tr_xyz_xyyyz[i] * tke_0 + 2.0 * tr_xyzz_xy[i] - 4.0 * tr_xyzz_xyzz[i] * tke_0 - 2.0 * tr_xyzz_xyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyyz[i] * tbe_0 - 2.0 * tr_xyyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xyzz[i] = 2.0 * tr_xz_xyzz[i] + 2.0 * tr_xzz_xyz[i] - 2.0 * tr_xzz_xyzzz[i] * tke_0 - 2.0 * tr_xzzz_xyzz[i] * tbe_0 + 2.0 * tr_xyz_xzz[i] - 4.0 * tr_xyz_xyyzz[i] * tke_0 + 2.0 * tr_xyzz_xz[i] - 2.0 * tr_xyzz_xzzz[i] * tke_0 - 4.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyzz[i] * tbe_0 - 4.0 * tr_xyyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xzzz[i] = 2.0 * tr_xz_xzzz[i] + 3.0 * tr_xzz_xzz[i] - 2.0 * tr_xzz_xzzzz[i] * tke_0 - 2.0 * tr_xzzz_xzzz[i] * tbe_0 - 4.0 * tr_xyz_xyzzz[i] * tke_0 - 6.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xzzz[i] * tbe_0 - 6.0 * tr_xyyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_yyyy[i] = 2.0 * tr_xz_yyyy[i] - 2.0 * tr_xzz_yyyyz[i] * tke_0 - 2.0 * tr_xzzz_yyyy[i] * tbe_0 + 8.0 * tr_xyz_yyy[i] - 4.0 * tr_xyz_yyyyy[i] * tke_0 - 8.0 * tr_xyzz_yyyz[i] * tke_0 + 4.0 * tr_xyzz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_yyyz[i] = 2.0 * tr_xz_yyyz[i] + tr_xzz_yyy[i] - 2.0 * tr_xzz_yyyzz[i] * tke_0 - 2.0 * tr_xzzz_yyyz[i] * tbe_0 + 6.0 * tr_xyz_yyz[i] - 4.0 * tr_xyz_yyyyz[i] * tke_0 + 3.0 * tr_xyzz_yy[i] - 6.0 * tr_xyzz_yyzz[i] * tke_0 - 2.0 * tr_xyzz_yyyy[i] * tke_0 + 4.0 * tr_xyzz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xyzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyyz[i] * tbe_0 - 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_yyzz[i] = 2.0 * tr_xz_yyzz[i] + 2.0 * tr_xzz_yyz[i] - 2.0 * tr_xzz_yyzzz[i] * tke_0 - 2.0 * tr_xzzz_yyzz[i] * tbe_0 + 4.0 * tr_xyz_yzz[i] - 4.0 * tr_xyz_yyyzz[i] * tke_0 + 4.0 * tr_xyzz_yz[i] - 4.0 * tr_xyzz_yzzz[i] * tke_0 - 4.0 * tr_xyzz_yyyz[i] * tke_0 + 4.0 * tr_xyzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyzz[i] * tbe_0 - 4.0 * tr_xyyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_yzzz[i] = 2.0 * tr_xz_yzzz[i] + 3.0 * tr_xzz_yzz[i] - 2.0 * tr_xzz_yzzzz[i] * tke_0 - 2.0 * tr_xzzz_yzzz[i] * tbe_0 + 2.0 * tr_xyz_zzz[i] - 4.0 * tr_xyz_yyzzz[i] * tke_0 + 3.0 * tr_xyzz_zz[i] - 2.0 * tr_xyzz_zzzz[i] * tke_0 - 6.0 * tr_xyzz_yyzz[i] * tke_0 + 4.0 * tr_xyzz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yzzz[i] * tbe_0 - 6.0 * tr_xyyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_zzzz[i] = 2.0 * tr_xz_zzzz[i] + 4.0 * tr_xzz_zzz[i] - 2.0 * tr_xzz_zzzzz[i] * tke_0 - 2.0 * tr_xzzz_zzzz[i] * tbe_0 - 4.0 * tr_xyz_yzzzz[i] * tke_0 - 8.0 * tr_xyzz_yzzz[i] * tke_0 + 4.0 * tr_xyzz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_yzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_zzzz[i] * tbe_0 - 8.0 * tr_xyyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1035-1050 components of targeted buffer : GG

    auto tr_0_0_yz_xzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1035);

    auto tr_0_0_yz_xzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1036);

    auto tr_0_0_yz_xzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1037);

    auto tr_0_0_yz_xzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1038);

    auto tr_0_0_yz_xzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1039);

    auto tr_0_0_yz_xzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1040);

    auto tr_0_0_yz_xzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1041);

    auto tr_0_0_yz_xzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1042);

    auto tr_0_0_yz_xzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1043);

    auto tr_0_0_yz_xzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1044);

    auto tr_0_0_yz_xzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1045);

    auto tr_0_0_yz_xzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1046);

    auto tr_0_0_yz_xzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1047);

    auto tr_0_0_yz_xzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1048);

    auto tr_0_0_yz_xzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1049);

    #pragma omp simd aligned(tr_0_0_yz_xzzz_xxxx, tr_0_0_yz_xzzz_xxxy, tr_0_0_yz_xzzz_xxxz, tr_0_0_yz_xzzz_xxyy, tr_0_0_yz_xzzz_xxyz, tr_0_0_yz_xzzz_xxzz, tr_0_0_yz_xzzz_xyyy, tr_0_0_yz_xzzz_xyyz, tr_0_0_yz_xzzz_xyzz, tr_0_0_yz_xzzz_xzzz, tr_0_0_yz_xzzz_yyyy, tr_0_0_yz_xzzz_yyyz, tr_0_0_yz_xzzz_yyzz, tr_0_0_yz_xzzz_yzzz, tr_0_0_yz_xzzz_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxz, tr_xyzzz_xxxyz, tr_xyzzz_xxxzz, tr_xyzzz_xxy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xxzzz, tr_xyzzz_xyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_xzzzz, tr_xyzzz_yyy, tr_xyzzz_yyyyz, tr_xyzzz_yyyzz, tr_xyzzz_yyz, tr_xyzzz_yyzzz, tr_xyzzz_yzz, tr_xyzzz_yzzzz, tr_xyzzz_zzz, tr_xyzzz_zzzzz, tr_xyzzzz_xxxx, tr_xyzzzz_xxxy, tr_xyzzzz_xxxz, tr_xyzzzz_xxyy, tr_xyzzzz_xxyz, tr_xyzzzz_xxzz, tr_xyzzzz_xyyy, tr_xyzzzz_xyyz, tr_xyzzzz_xyzz, tr_xyzzzz_xzzz, tr_xyzzzz_yyyy, tr_xyzzzz_yyyz, tr_xyzzzz_yyzz, tr_xyzzzz_yzzz, tr_xyzzzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xxxxyz, tr_xzzz_xxxy, tr_xzzz_xxxyyz, tr_xzzz_xxxyzz, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyyyz, tr_xzzz_xxyyzz, tr_xzzz_xxyz, tr_xzzz_xxyzzz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyyyz, tr_xzzz_xyyyzz, tr_xzzz_xyyz, tr_xzzz_xyyzzz, tr_xzzz_xyzz, tr_xzzz_xyzzzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyyyz, tr_xzzz_yyyyzz, tr_xzzz_yyyz, tr_xzzz_yyyzzz, tr_xzzz_yyzz, tr_xzzz_yyzzzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_yzzzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxxxy, tr_xzzzz_xxxyy, tr_xzzzz_xxxyz, tr_xzzzz_xxy, tr_xzzzz_xxyyy, tr_xzzzz_xxyyz, tr_xzzzz_xxyzz, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyyyy, tr_xzzzz_xyyyz, tr_xzzzz_xyyzz, tr_xzzzz_xyz, tr_xzzzz_xyzzz, tr_xzzzz_xzz, tr_xzzzz_yyy, tr_xzzzz_yyyyy, tr_xzzzz_yyyyz, tr_xzzzz_yyyzz, tr_xzzzz_yyz, tr_xzzzz_yyzzz, tr_xzzzz_yzz, tr_xzzzz_yzzzz, tr_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xzzz_xxxx[i] = -6.0 * tr_xzz_xxxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xxxy[i] = 3.0 * tr_xzz_xxx[i] - 6.0 * tr_xzz_xxxyy[i] * tke_0 - 2.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xxxz[i] = -6.0 * tr_xzz_xxxyz[i] * tke_0 - 2.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxxz[i] * tbe_0 - 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xxyy[i] = 6.0 * tr_xzz_xxy[i] - 6.0 * tr_xzz_xxyyy[i] * tke_0 - 4.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxyy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xxyz[i] = 3.0 * tr_xzz_xxz[i] - 6.0 * tr_xzz_xxyyz[i] * tke_0 + tr_xzzz_xx[i] - 2.0 * tr_xzzz_xxzz[i] * tke_0 - 2.0 * tr_xzzz_xxyy[i] * tke_0 + 4.0 * tr_xzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxyz[i] * tbe_0 - 2.0 * tr_xyzzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xxzz[i] = -6.0 * tr_xzz_xxyzz[i] * tke_0 - 4.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxzz[i] * tbe_0 - 4.0 * tr_xyzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xyyy[i] = 9.0 * tr_xzz_xyy[i] - 6.0 * tr_xzz_xyyyy[i] * tke_0 - 6.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_xzzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xyyz[i] = 6.0 * tr_xzz_xyz[i] - 6.0 * tr_xzz_xyyyz[i] * tke_0 + 2.0 * tr_xzzz_xy[i] - 4.0 * tr_xzzz_xyzz[i] * tke_0 - 2.0 * tr_xzzz_xyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyyz[i] * tbe_0 - 2.0 * tr_xyzzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xyzz[i] = 3.0 * tr_xzz_xzz[i] - 6.0 * tr_xzz_xyyzz[i] * tke_0 + 2.0 * tr_xzzz_xz[i] - 2.0 * tr_xzzz_xzzz[i] * tke_0 - 4.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyzz[i] * tbe_0 - 4.0 * tr_xyzzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xzzz[i] = -6.0 * tr_xzz_xyzzz[i] * tke_0 - 6.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xzzz[i] * tbe_0 - 6.0 * tr_xyzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_yyyy[i] = 12.0 * tr_xzz_yyy[i] - 6.0 * tr_xzz_yyyyy[i] * tke_0 - 8.0 * tr_xzzz_yyyz[i] * tke_0 + 4.0 * tr_xzzz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_yyyz[i] = 9.0 * tr_xzz_yyz[i] - 6.0 * tr_xzz_yyyyz[i] * tke_0 + 3.0 * tr_xzzz_yy[i] - 6.0 * tr_xzzz_yyzz[i] * tke_0 - 2.0 * tr_xzzz_yyyy[i] * tke_0 + 4.0 * tr_xzzz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_xzzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyyz[i] * tbe_0 - 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_yyzz[i] = 6.0 * tr_xzz_yzz[i] - 6.0 * tr_xzz_yyyzz[i] * tke_0 + 4.0 * tr_xzzz_yz[i] - 4.0 * tr_xzzz_yzzz[i] * tke_0 - 4.0 * tr_xzzz_yyyz[i] * tke_0 + 4.0 * tr_xzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyzz[i] * tbe_0 - 4.0 * tr_xyzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_yzzz[i] = 3.0 * tr_xzz_zzz[i] - 6.0 * tr_xzz_yyzzz[i] * tke_0 + 3.0 * tr_xzzz_zz[i] - 2.0 * tr_xzzz_zzzz[i] * tke_0 - 6.0 * tr_xzzz_yyzz[i] * tke_0 + 4.0 * tr_xzzz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzzz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yzzz[i] * tbe_0 - 6.0 * tr_xyzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_zzzz[i] = -6.0 * tr_xzz_yzzzz[i] * tke_0 - 8.0 * tr_xzzz_yzzz[i] * tke_0 + 4.0 * tr_xzzz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_yzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_zzzz[i] * tbe_0 - 8.0 * tr_xyzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1050-1065 components of targeted buffer : GG

    auto tr_0_0_yz_yyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 1050);

    auto tr_0_0_yz_yyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 1051);

    auto tr_0_0_yz_yyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 1052);

    auto tr_0_0_yz_yyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 1053);

    auto tr_0_0_yz_yyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 1054);

    auto tr_0_0_yz_yyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 1055);

    auto tr_0_0_yz_yyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 1056);

    auto tr_0_0_yz_yyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 1057);

    auto tr_0_0_yz_yyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 1058);

    auto tr_0_0_yz_yyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 1059);

    auto tr_0_0_yz_yyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 1060);

    auto tr_0_0_yz_yyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 1061);

    auto tr_0_0_yz_yyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 1062);

    auto tr_0_0_yz_yyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 1063);

    auto tr_0_0_yz_yyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 1064);

    #pragma omp simd aligned(tr_0_0_yz_yyyy_xxxx, tr_0_0_yz_yyyy_xxxy, tr_0_0_yz_yyyy_xxxz, tr_0_0_yz_yyyy_xxyy, tr_0_0_yz_yyyy_xxyz, tr_0_0_yz_yyyy_xxzz, tr_0_0_yz_yyyy_xyyy, tr_0_0_yz_yyyy_xyyz, tr_0_0_yz_yyyy_xyzz, tr_0_0_yz_yyyy_xzzz, tr_0_0_yz_yyyy_yyyy, tr_0_0_yz_yyyy_yyyz, tr_0_0_yz_yyyy_yyzz, tr_0_0_yz_yyyy_yzzz, tr_0_0_yz_yyyy_zzzz, tr_yyy_xxx, tr_yyy_xxxxz, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyy_zzzzz, tr_yyyy_xx, tr_yyyy_xxxxyz, tr_yyyy_xxxy, tr_yyyy_xxxyyz, tr_yyyy_xxxyzz, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyyyz, tr_yyyy_xxyyzz, tr_yyyy_xxyz, tr_yyyy_xxyzzz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyyyz, tr_yyyy_xyyyzz, tr_yyyy_xyyz, tr_yyyy_xyyzzz, tr_yyyy_xyzz, tr_yyyy_xyzzzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyyyz, tr_yyyy_yyyyzz, tr_yyyy_yyyz, tr_yyyy_yyyzzz, tr_yyyy_yyzz, tr_yyyy_yyzzzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_yzzzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyyy_xxx, tr_yyyyy_xxxxz, tr_yyyyy_xxxyz, tr_yyyyy_xxxzz, tr_yyyyy_xxy, tr_yyyyy_xxyyz, tr_yyyyy_xxyzz, tr_yyyyy_xxz, tr_yyyyy_xxzzz, tr_yyyyy_xyy, tr_yyyyy_xyyyz, tr_yyyyy_xyyzz, tr_yyyyy_xyz, tr_yyyyy_xyzzz, tr_yyyyy_xzz, tr_yyyyy_xzzzz, tr_yyyyy_yyy, tr_yyyyy_yyyyz, tr_yyyyy_yyyzz, tr_yyyyy_yyz, tr_yyyyy_yyzzz, tr_yyyyy_yzz, tr_yyyyy_yzzzz, tr_yyyyy_zzz, tr_yyyyy_zzzzz, tr_yyyyyz_xxxx, tr_yyyyyz_xxxy, tr_yyyyyz_xxxz, tr_yyyyyz_xxyy, tr_yyyyyz_xxyz, tr_yyyyyz_xxzz, tr_yyyyyz_xyyy, tr_yyyyyz_xyyz, tr_yyyyyz_xyzz, tr_yyyyyz_xzzz, tr_yyyyyz_yyyy, tr_yyyyyz_yyyz, tr_yyyyyz_yyzz, tr_yyyyyz_yzzz, tr_yyyyyz_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxxxy, tr_yyyyz_xxxyy, tr_yyyyz_xxxyz, tr_yyyyz_xxy, tr_yyyyz_xxyyy, tr_yyyyz_xxyyz, tr_yyyyz_xxyzz, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyyyy, tr_yyyyz_xyyyz, tr_yyyyz_xyyzz, tr_yyyyz_xyz, tr_yyyyz_xyzzz, tr_yyyyz_xzz, tr_yyyyz_yyy, tr_yyyyz_yyyyy, tr_yyyyz_yyyyz, tr_yyyyz_yyyzz, tr_yyyyz_yyz, tr_yyyyz_yyzzz, tr_yyyyz_yzz, tr_yyyyz_yzzzz, tr_yyyyz_zzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyyy_xxxx[i] = -8.0 * tr_yyy_xxxxz[i] * tke_0 - 8.0 * tr_yyyz_xxxx[i] * tbe_0 + 4.0 * tr_yyyy_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xxxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xxxy[i] = -8.0 * tr_yyy_xxxyz[i] * tke_0 - 8.0 * tr_yyyz_xxxy[i] * tbe_0 - 2.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xxxz[i] = 4.0 * tr_yyy_xxx[i] - 8.0 * tr_yyy_xxxzz[i] * tke_0 - 8.0 * tr_yyyz_xxxz[i] * tbe_0 - 2.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyyy_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xxyy[i] = -8.0 * tr_yyy_xxyyz[i] * tke_0 - 8.0 * tr_yyyz_xxyy[i] * tbe_0 - 4.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xxyz[i] = 4.0 * tr_yyy_xxy[i] - 8.0 * tr_yyy_xxyzz[i] * tke_0 - 8.0 * tr_yyyz_xxyz[i] * tbe_0 + tr_yyyy_xx[i] - 2.0 * tr_yyyy_xxzz[i] * tke_0 - 2.0 * tr_yyyy_xxyy[i] * tke_0 + 4.0 * tr_yyyy_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyyy_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xxzz[i] = 8.0 * tr_yyy_xxz[i] - 8.0 * tr_yyy_xxzzz[i] * tke_0 - 8.0 * tr_yyyz_xxzz[i] * tbe_0 - 4.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyyy_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xyyy[i] = -8.0 * tr_yyy_xyyyz[i] * tke_0 - 8.0 * tr_yyyz_xyyy[i] * tbe_0 - 6.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xyyz[i] = 4.0 * tr_yyy_xyy[i] - 8.0 * tr_yyy_xyyzz[i] * tke_0 - 8.0 * tr_yyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyyy_xy[i] - 4.0 * tr_yyyy_xyzz[i] * tke_0 - 2.0 * tr_yyyy_xyyy[i] * tke_0 + 4.0 * tr_yyyy_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyyz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyyy_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xyzz[i] = 8.0 * tr_yyy_xyz[i] - 8.0 * tr_yyy_xyzzz[i] * tke_0 - 8.0 * tr_yyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyyy_xz[i] - 2.0 * tr_yyyy_xzzz[i] * tke_0 - 4.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyyz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyyy_xyz[i] * tbe_0 + 4.0 * tr_yyyyy_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xzzz[i] = 12.0 * tr_yyy_xzz[i] - 8.0 * tr_yyy_xzzzz[i] * tke_0 - 8.0 * tr_yyyz_xzzz[i] * tbe_0 - 6.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyyy_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_yyyy[i] = -8.0 * tr_yyy_yyyyz[i] * tke_0 - 8.0 * tr_yyyz_yyyy[i] * tbe_0 - 8.0 * tr_yyyy_yyyz[i] * tke_0 + 4.0 * tr_yyyy_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_yyyz[i] = 4.0 * tr_yyy_yyy[i] - 8.0 * tr_yyy_yyyzz[i] * tke_0 - 8.0 * tr_yyyz_yyyz[i] * tbe_0 + 3.0 * tr_yyyy_yy[i] - 6.0 * tr_yyyy_yyzz[i] * tke_0 - 2.0 * tr_yyyy_yyyy[i] * tke_0 + 4.0 * tr_yyyy_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_yyyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyyz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyyy_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_yyzz[i] = 8.0 * tr_yyy_yyz[i] - 8.0 * tr_yyy_yyzzz[i] * tke_0 - 8.0 * tr_yyyz_yyzz[i] * tbe_0 + 4.0 * tr_yyyy_yz[i] - 4.0 * tr_yyyy_yzzz[i] * tke_0 - 4.0 * tr_yyyy_yyyz[i] * tke_0 + 4.0 * tr_yyyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyyz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyyy_yyz[i] * tbe_0 + 4.0 * tr_yyyyy_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_yzzz[i] = 12.0 * tr_yyy_yzz[i] - 8.0 * tr_yyy_yzzzz[i] * tke_0 - 8.0 * tr_yyyz_yzzz[i] * tbe_0 + 3.0 * tr_yyyy_zz[i] - 2.0 * tr_yyyy_zzzz[i] * tke_0 - 6.0 * tr_yyyy_yyzz[i] * tke_0 + 4.0 * tr_yyyy_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyyz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyyyy_yzz[i] * tbe_0 + 4.0 * tr_yyyyy_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_zzzz[i] = 16.0 * tr_yyy_zzz[i] - 8.0 * tr_yyy_zzzzz[i] * tke_0 - 8.0 * tr_yyyz_zzzz[i] * tbe_0 - 8.0 * tr_yyyy_yzzz[i] * tke_0 + 4.0 * tr_yyyy_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_yyyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyyy_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1065-1080 components of targeted buffer : GG

    auto tr_0_0_yz_yyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1065);

    auto tr_0_0_yz_yyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1066);

    auto tr_0_0_yz_yyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1067);

    auto tr_0_0_yz_yyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1068);

    auto tr_0_0_yz_yyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1069);

    auto tr_0_0_yz_yyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1070);

    auto tr_0_0_yz_yyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1071);

    auto tr_0_0_yz_yyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1072);

    auto tr_0_0_yz_yyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1073);

    auto tr_0_0_yz_yyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1074);

    auto tr_0_0_yz_yyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1075);

    auto tr_0_0_yz_yyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1076);

    auto tr_0_0_yz_yyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1077);

    auto tr_0_0_yz_yyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1078);

    auto tr_0_0_yz_yyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1079);

    #pragma omp simd aligned(tr_0_0_yz_yyyz_xxxx, tr_0_0_yz_yyyz_xxxy, tr_0_0_yz_yyyz_xxxz, tr_0_0_yz_yyyz_xxyy, tr_0_0_yz_yyyz_xxyz, tr_0_0_yz_yyyz_xxzz, tr_0_0_yz_yyyz_xyyy, tr_0_0_yz_yyyz_xyyz, tr_0_0_yz_yyyz_xyzz, tr_0_0_yz_yyyz_xzzz, tr_0_0_yz_yyyz_yyyy, tr_0_0_yz_yyyz_yyyz, tr_0_0_yz_yyyz_yyzz, tr_0_0_yz_yyyz_yzzz, tr_0_0_yz_yyyz_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyy_xxx, tr_yyy_xxxxy, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyyyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xzzz, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yzzz, tr_yyyy_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxxxz, tr_yyyyz_xxxyz, tr_yyyyz_xxxzz, tr_yyyyz_xxy, tr_yyyyz_xxyyz, tr_yyyyz_xxyzz, tr_yyyyz_xxz, tr_yyyyz_xxzzz, tr_yyyyz_xyy, tr_yyyyz_xyyyz, tr_yyyyz_xyyzz, tr_yyyyz_xyz, tr_yyyyz_xyzzz, tr_yyyyz_xzz, tr_yyyyz_xzzzz, tr_yyyyz_yyy, tr_yyyyz_yyyyz, tr_yyyyz_yyyzz, tr_yyyyz_yyz, tr_yyyyz_yyzzz, tr_yyyyz_yzz, tr_yyyyz_yzzzz, tr_yyyyz_zzz, tr_yyyyz_zzzzz, tr_yyyyzz_xxxx, tr_yyyyzz_xxxy, tr_yyyyzz_xxxz, tr_yyyyzz_xxyy, tr_yyyyzz_xxyz, tr_yyyyzz_xxzz, tr_yyyyzz_xyyy, tr_yyyyzz_xyyz, tr_yyyyzz_xyzz, tr_yyyyzz_xzzz, tr_yyyyzz_yyyy, tr_yyyyzz_yyyz, tr_yyyyzz_yyzz, tr_yyyyzz_yzzz, tr_yyyyzz_zzzz, tr_yyyz_xx, tr_yyyz_xxxxyz, tr_yyyz_xxxy, tr_yyyz_xxxyyz, tr_yyyz_xxxyzz, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyyyz, tr_yyyz_xxyyzz, tr_yyyz_xxyz, tr_yyyz_xxyzzz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyyyz, tr_yyyz_xyyyzz, tr_yyyz_xyyz, tr_yyyz_xyyzzz, tr_yyyz_xyzz, tr_yyyz_xyzzzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyyyz, tr_yyyz_yyyyzz, tr_yyyz_yyyz, tr_yyyz_yyyzzz, tr_yyyz_yyzz, tr_yyyz_yyzzzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_yzzzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxxxy, tr_yyyzz_xxxyy, tr_yyyzz_xxxyz, tr_yyyzz_xxy, tr_yyyzz_xxyyy, tr_yyyzz_xxyyz, tr_yyyzz_xxyzz, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyyyy, tr_yyyzz_xyyyz, tr_yyyzz_xyyzz, tr_yyyzz_xyz, tr_yyyzz_xyzzz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyyyy, tr_yyyzz_yyyyz, tr_yyyzz_yyyzz, tr_yyyzz_yyz, tr_yyyzz_yyzzz, tr_yyyzz_yzz, tr_yyyzz_yzzzz, tr_yyyzz_zzz, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyyz_xxxx[i] = 3.0 * tr_yy_xxxx[i] - 6.0 * tr_yyz_xxxxz[i] * tke_0 - 6.0 * tr_yyzz_xxxx[i] * tbe_0 - 2.0 * tr_yyy_xxxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xxxxy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xxxy[i] = 3.0 * tr_yy_xxxy[i] - 6.0 * tr_yyz_xxxyz[i] * tke_0 - 6.0 * tr_yyzz_xxxy[i] * tbe_0 + tr_yyy_xxx[i] - 2.0 * tr_yyy_xxxyy[i] * tke_0 - 2.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xxxz[i] = 3.0 * tr_yy_xxxz[i] + 3.0 * tr_yyz_xxx[i] - 6.0 * tr_yyz_xxxzz[i] * tke_0 - 6.0 * tr_yyzz_xxxz[i] * tbe_0 - 2.0 * tr_yyy_xxxyz[i] * tke_0 - 2.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xxxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxxz[i] * tbe_0 - 2.0 * tr_yyyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xxyy[i] = 3.0 * tr_yy_xxyy[i] - 6.0 * tr_yyz_xxyyz[i] * tke_0 - 6.0 * tr_yyzz_xxyy[i] * tbe_0 + 2.0 * tr_yyy_xxy[i] - 2.0 * tr_yyy_xxyyy[i] * tke_0 - 4.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxyy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xxyz[i] = 3.0 * tr_yy_xxyz[i] + 3.0 * tr_yyz_xxy[i] - 6.0 * tr_yyz_xxyzz[i] * tke_0 - 6.0 * tr_yyzz_xxyz[i] * tbe_0 + tr_yyy_xxz[i] - 2.0 * tr_yyy_xxyyz[i] * tke_0 + tr_yyyz_xx[i] - 2.0 * tr_yyyz_xxzz[i] * tke_0 - 2.0 * tr_yyyz_xxyy[i] * tke_0 + 4.0 * tr_yyyz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_xxz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxyz[i] * tbe_0 - 2.0 * tr_yyyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xxzz[i] = 3.0 * tr_yy_xxzz[i] + 6.0 * tr_yyz_xxz[i] - 6.0 * tr_yyz_xxzzz[i] * tke_0 - 6.0 * tr_yyzz_xxzz[i] * tbe_0 - 2.0 * tr_yyy_xxyzz[i] * tke_0 - 4.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xxyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxzz[i] * tbe_0 - 4.0 * tr_yyyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xyyy[i] = 3.0 * tr_yy_xyyy[i] - 6.0 * tr_yyz_xyyyz[i] * tke_0 - 6.0 * tr_yyzz_xyyy[i] * tbe_0 + 3.0 * tr_yyy_xyy[i] - 2.0 * tr_yyy_xyyyy[i] * tke_0 - 6.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yyyzz_xyy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xyyy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xyyz[i] = 3.0 * tr_yy_xyyz[i] + 3.0 * tr_yyz_xyy[i] - 6.0 * tr_yyz_xyyzz[i] * tke_0 - 6.0 * tr_yyzz_xyyz[i] * tbe_0 + 2.0 * tr_yyy_xyz[i] - 2.0 * tr_yyy_xyyyz[i] * tke_0 + 2.0 * tr_yyyz_xy[i] - 4.0 * tr_yyyz_xyzz[i] * tke_0 - 2.0 * tr_yyyz_xyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xyz[i] * tbe_0 + 4.0 * tr_yyyzz_xyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xyyz[i] * tbe_0 - 2.0 * tr_yyyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xyzz[i] = 3.0 * tr_yy_xyzz[i] + 6.0 * tr_yyz_xyz[i] - 6.0 * tr_yyz_xyzzz[i] * tke_0 - 6.0 * tr_yyzz_xyzz[i] * tbe_0 + tr_yyy_xzz[i] - 2.0 * tr_yyy_xyyzz[i] * tke_0 + 2.0 * tr_yyyz_xz[i] - 2.0 * tr_yyyz_xzzz[i] * tke_0 - 4.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_xzz[i] * tbe_0 + 4.0 * tr_yyyzz_xyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xyzz[i] * tbe_0 - 4.0 * tr_yyyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xzzz[i] = 3.0 * tr_yy_xzzz[i] + 9.0 * tr_yyz_xzz[i] - 6.0 * tr_yyz_xzzzz[i] * tke_0 - 6.0 * tr_yyzz_xzzz[i] * tbe_0 - 2.0 * tr_yyy_xyzzz[i] * tke_0 - 6.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xzzz[i] * tbe_0 - 6.0 * tr_yyyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_yyyy[i] = 3.0 * tr_yy_yyyy[i] - 6.0 * tr_yyz_yyyyz[i] * tke_0 - 6.0 * tr_yyzz_yyyy[i] * tbe_0 + 4.0 * tr_yyy_yyy[i] - 2.0 * tr_yyy_yyyyy[i] * tke_0 - 8.0 * tr_yyyz_yyyz[i] * tke_0 + 4.0 * tr_yyyz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_yyy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yyyy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_yyyz[i] = 3.0 * tr_yy_yyyz[i] + 3.0 * tr_yyz_yyy[i] - 6.0 * tr_yyz_yyyzz[i] * tke_0 - 6.0 * tr_yyzz_yyyz[i] * tbe_0 + 3.0 * tr_yyy_yyz[i] - 2.0 * tr_yyy_yyyyz[i] * tke_0 + 3.0 * tr_yyyz_yy[i] - 6.0 * tr_yyyz_yyzz[i] * tke_0 - 2.0 * tr_yyyz_yyyy[i] * tke_0 + 4.0 * tr_yyyz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_yyyzz_yyz[i] * tbe_0 + 4.0 * tr_yyyzz_yyyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yyyz[i] * tbe_0 - 2.0 * tr_yyyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_yyzz[i] = 3.0 * tr_yy_yyzz[i] + 6.0 * tr_yyz_yyz[i] - 6.0 * tr_yyz_yyzzz[i] * tke_0 - 6.0 * tr_yyzz_yyzz[i] * tbe_0 + 2.0 * tr_yyy_yzz[i] - 2.0 * tr_yyy_yyyzz[i] * tke_0 + 4.0 * tr_yyyz_yz[i] - 4.0 * tr_yyyz_yzzz[i] * tke_0 - 4.0 * tr_yyyz_yyyz[i] * tke_0 + 4.0 * tr_yyyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_yzz[i] * tbe_0 + 4.0 * tr_yyyzz_yyyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yyzz[i] * tbe_0 - 4.0 * tr_yyyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_yzzz[i] = 3.0 * tr_yy_yzzz[i] + 9.0 * tr_yyz_yzz[i] - 6.0 * tr_yyz_yzzzz[i] * tke_0 - 6.0 * tr_yyzz_yzzz[i] * tbe_0 + tr_yyy_zzz[i] - 2.0 * tr_yyy_yyzzz[i] * tke_0 + 3.0 * tr_yyyz_zz[i] - 2.0 * tr_yyyz_zzzz[i] * tke_0 - 6.0 * tr_yyyz_yyzz[i] * tke_0 + 4.0 * tr_yyyz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_zzz[i] * tbe_0 + 4.0 * tr_yyyzz_yyzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yzzz[i] * tbe_0 - 6.0 * tr_yyyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_zzzz[i] = 3.0 * tr_yy_zzzz[i] + 12.0 * tr_yyz_zzz[i] - 6.0 * tr_yyz_zzzzz[i] * tke_0 - 6.0 * tr_yyzz_zzzz[i] * tbe_0 - 2.0 * tr_yyy_yzzzz[i] * tke_0 - 8.0 * tr_yyyz_yzzz[i] * tke_0 + 4.0 * tr_yyyz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_yzzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_zzzz[i] * tbe_0 - 8.0 * tr_yyyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1080-1095 components of targeted buffer : GG

    auto tr_0_0_yz_yyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1080);

    auto tr_0_0_yz_yyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1081);

    auto tr_0_0_yz_yyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1082);

    auto tr_0_0_yz_yyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1083);

    auto tr_0_0_yz_yyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1084);

    auto tr_0_0_yz_yyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1085);

    auto tr_0_0_yz_yyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1086);

    auto tr_0_0_yz_yyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1087);

    auto tr_0_0_yz_yyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1088);

    auto tr_0_0_yz_yyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1089);

    auto tr_0_0_yz_yyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1090);

    auto tr_0_0_yz_yyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1091);

    auto tr_0_0_yz_yyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1092);

    auto tr_0_0_yz_yyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1093);

    auto tr_0_0_yz_yyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1094);

    #pragma omp simd aligned(tr_0_0_yz_yyzz_xxxx, tr_0_0_yz_yyzz_xxxy, tr_0_0_yz_yyzz_xxxz, tr_0_0_yz_yyzz_xxyy, tr_0_0_yz_yyzz_xxyz, tr_0_0_yz_yyzz_xxzz, tr_0_0_yz_yyzz_xyyy, tr_0_0_yz_yyzz_xyyz, tr_0_0_yz_yyzz_xyzz, tr_0_0_yz_yyzz_xzzz, tr_0_0_yz_yyzz_yyyy, tr_0_0_yz_yyzz_yyyz, tr_0_0_yz_yyzz_yyzz, tr_0_0_yz_yyzz_yzzz, tr_0_0_yz_yyzz_zzzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxxxz, tr_yyyzz_xxxyz, tr_yyyzz_xxxzz, tr_yyyzz_xxy, tr_yyyzz_xxyyz, tr_yyyzz_xxyzz, tr_yyyzz_xxz, tr_yyyzz_xxzzz, tr_yyyzz_xyy, tr_yyyzz_xyyyz, tr_yyyzz_xyyzz, tr_yyyzz_xyz, tr_yyyzz_xyzzz, tr_yyyzz_xzz, tr_yyyzz_xzzzz, tr_yyyzz_yyy, tr_yyyzz_yyyyz, tr_yyyzz_yyyzz, tr_yyyzz_yyz, tr_yyyzz_yyzzz, tr_yyyzz_yzz, tr_yyyzz_yzzzz, tr_yyyzz_zzz, tr_yyyzz_zzzzz, tr_yyyzzz_xxxx, tr_yyyzzz_xxxy, tr_yyyzzz_xxxz, tr_yyyzzz_xxyy, tr_yyyzzz_xxyz, tr_yyyzzz_xxzz, tr_yyyzzz_xyyy, tr_yyyzzz_xyyz, tr_yyyzzz_xyzz, tr_yyyzzz_xzzz, tr_yyyzzz_yyyy, tr_yyyzzz_yyyz, tr_yyyzzz_yyzz, tr_yyyzzz_yzzz, tr_yyyzzz_zzzz, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xxxxyz, tr_yyzz_xxxy, tr_yyzz_xxxyyz, tr_yyzz_xxxyzz, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyyyz, tr_yyzz_xxyyzz, tr_yyzz_xxyz, tr_yyzz_xxyzzz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyyyz, tr_yyzz_xyyyzz, tr_yyzz_xyyz, tr_yyzz_xyyzzz, tr_yyzz_xyzz, tr_yyzz_xyzzzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyyyz, tr_yyzz_yyyyzz, tr_yyzz_yyyz, tr_yyzz_yyyzzz, tr_yyzz_yyzz, tr_yyzz_yyzzzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_yzzzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxxxy, tr_yyzzz_xxxyy, tr_yyzzz_xxxyz, tr_yyzzz_xxy, tr_yyzzz_xxyyy, tr_yyzzz_xxyyz, tr_yyzzz_xxyzz, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyyyy, tr_yyzzz_xyyyz, tr_yyzzz_xyyzz, tr_yyzzz_xyz, tr_yyzzz_xyzzz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyyyy, tr_yyzzz_yyyyz, tr_yyzzz_yyyzz, tr_yyzzz_yyz, tr_yyzzz_yyzzz, tr_yyzzz_yzz, tr_yyzzz_yzzzz, tr_yyzzz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyzz_xxxx[i] = 4.0 * tr_yz_xxxx[i] - 4.0 * tr_yzz_xxxxz[i] * tke_0 - 4.0 * tr_yzzz_xxxx[i] * tbe_0 - 4.0 * tr_yyz_xxxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xxxxy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xxxy[i] = 4.0 * tr_yz_xxxy[i] - 4.0 * tr_yzz_xxxyz[i] * tke_0 - 4.0 * tr_yzzz_xxxy[i] * tbe_0 + 2.0 * tr_yyz_xxx[i] - 4.0 * tr_yyz_xxxyy[i] * tke_0 - 2.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xxxz[i] = 4.0 * tr_yz_xxxz[i] + 2.0 * tr_yzz_xxx[i] - 4.0 * tr_yzz_xxxzz[i] * tke_0 - 4.0 * tr_yzzz_xxxz[i] * tbe_0 - 4.0 * tr_yyz_xxxyz[i] * tke_0 - 2.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xxxyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxxz[i] * tbe_0 - 2.0 * tr_yyyzz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xxyy[i] = 4.0 * tr_yz_xxyy[i] - 4.0 * tr_yzz_xxyyz[i] * tke_0 - 4.0 * tr_yzzz_xxyy[i] * tbe_0 + 4.0 * tr_yyz_xxy[i] - 4.0 * tr_yyz_xxyyy[i] * tke_0 - 4.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxyy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xxyz[i] = 4.0 * tr_yz_xxyz[i] + 2.0 * tr_yzz_xxy[i] - 4.0 * tr_yzz_xxyzz[i] * tke_0 - 4.0 * tr_yzzz_xxyz[i] * tbe_0 + 2.0 * tr_yyz_xxz[i] - 4.0 * tr_yyz_xxyyz[i] * tke_0 + tr_yyzz_xx[i] - 2.0 * tr_yyzz_xxzz[i] * tke_0 - 2.0 * tr_yyzz_xxyy[i] * tke_0 + 4.0 * tr_yyzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_xxz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxyz[i] * tbe_0 - 2.0 * tr_yyyzz_xxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xxzz[i] = 4.0 * tr_yz_xxzz[i] + 4.0 * tr_yzz_xxz[i] - 4.0 * tr_yzz_xxzzz[i] * tke_0 - 4.0 * tr_yzzz_xxzz[i] * tbe_0 - 4.0 * tr_yyz_xxyzz[i] * tke_0 - 4.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xxyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxzz[i] * tbe_0 - 4.0 * tr_yyyzz_xxz[i] * tbe_0 + 4.0 * tr_yyyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xyyy[i] = 4.0 * tr_yz_xyyy[i] - 4.0 * tr_yzz_xyyyz[i] * tke_0 - 4.0 * tr_yzzz_xyyy[i] * tbe_0 + 6.0 * tr_yyz_xyy[i] - 4.0 * tr_yyz_xyyyy[i] * tke_0 - 6.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yyzzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xyyy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xyyz[i] = 4.0 * tr_yz_xyyz[i] + 2.0 * tr_yzz_xyy[i] - 4.0 * tr_yzz_xyyzz[i] * tke_0 - 4.0 * tr_yzzz_xyyz[i] * tbe_0 + 4.0 * tr_yyz_xyz[i] - 4.0 * tr_yyz_xyyyz[i] * tke_0 + 2.0 * tr_yyzz_xy[i] - 4.0 * tr_yyzz_xyzz[i] * tke_0 - 2.0 * tr_yyzz_xyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xyz[i] * tbe_0 + 4.0 * tr_yyzzz_xyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xyyz[i] * tbe_0 - 2.0 * tr_yyyzz_xyy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xyzz[i] = 4.0 * tr_yz_xyzz[i] + 4.0 * tr_yzz_xyz[i] - 4.0 * tr_yzz_xyzzz[i] * tke_0 - 4.0 * tr_yzzz_xyzz[i] * tbe_0 + 2.0 * tr_yyz_xzz[i] - 4.0 * tr_yyz_xyyzz[i] * tke_0 + 2.0 * tr_yyzz_xz[i] - 2.0 * tr_yyzz_xzzz[i] * tke_0 - 4.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_xzz[i] * tbe_0 + 4.0 * tr_yyzzz_xyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xyzz[i] * tbe_0 - 4.0 * tr_yyyzz_xyz[i] * tbe_0 + 4.0 * tr_yyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xzzz[i] = 4.0 * tr_yz_xzzz[i] + 6.0 * tr_yzz_xzz[i] - 4.0 * tr_yzz_xzzzz[i] * tke_0 - 4.0 * tr_yzzz_xzzz[i] * tbe_0 - 4.0 * tr_yyz_xyzzz[i] * tke_0 - 6.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xzzz[i] * tbe_0 - 6.0 * tr_yyyzz_xzz[i] * tbe_0 + 4.0 * tr_yyyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_yyyy[i] = 4.0 * tr_yz_yyyy[i] - 4.0 * tr_yzz_yyyyz[i] * tke_0 - 4.0 * tr_yzzz_yyyy[i] * tbe_0 + 8.0 * tr_yyz_yyy[i] - 4.0 * tr_yyz_yyyyy[i] * tke_0 - 8.0 * tr_yyzz_yyyz[i] * tke_0 + 4.0 * tr_yyzz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yyyy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_yyyz[i] = 4.0 * tr_yz_yyyz[i] + 2.0 * tr_yzz_yyy[i] - 4.0 * tr_yzz_yyyzz[i] * tke_0 - 4.0 * tr_yzzz_yyyz[i] * tbe_0 + 6.0 * tr_yyz_yyz[i] - 4.0 * tr_yyz_yyyyz[i] * tke_0 + 3.0 * tr_yyzz_yy[i] - 6.0 * tr_yyzz_yyzz[i] * tke_0 - 2.0 * tr_yyzz_yyyy[i] * tke_0 + 4.0 * tr_yyzz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_yyzzz_yyz[i] * tbe_0 + 4.0 * tr_yyzzz_yyyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yyyz[i] * tbe_0 - 2.0 * tr_yyyzz_yyy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_yyzz[i] = 4.0 * tr_yz_yyzz[i] + 4.0 * tr_yzz_yyz[i] - 4.0 * tr_yzz_yyzzz[i] * tke_0 - 4.0 * tr_yzzz_yyzz[i] * tbe_0 + 4.0 * tr_yyz_yzz[i] - 4.0 * tr_yyz_yyyzz[i] * tke_0 + 4.0 * tr_yyzz_yz[i] - 4.0 * tr_yyzz_yzzz[i] * tke_0 - 4.0 * tr_yyzz_yyyz[i] * tke_0 + 4.0 * tr_yyzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_yzz[i] * tbe_0 + 4.0 * tr_yyzzz_yyyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yyzz[i] * tbe_0 - 4.0 * tr_yyyzz_yyz[i] * tbe_0 + 4.0 * tr_yyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_yzzz[i] = 4.0 * tr_yz_yzzz[i] + 6.0 * tr_yzz_yzz[i] - 4.0 * tr_yzz_yzzzz[i] * tke_0 - 4.0 * tr_yzzz_yzzz[i] * tbe_0 + 2.0 * tr_yyz_zzz[i] - 4.0 * tr_yyz_yyzzz[i] * tke_0 + 3.0 * tr_yyzz_zz[i] - 2.0 * tr_yyzz_zzzz[i] * tke_0 - 6.0 * tr_yyzz_yyzz[i] * tke_0 + 4.0 * tr_yyzz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_zzz[i] * tbe_0 + 4.0 * tr_yyzzz_yyzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yzzz[i] * tbe_0 - 6.0 * tr_yyyzz_yzz[i] * tbe_0 + 4.0 * tr_yyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_zzzz[i] = 4.0 * tr_yz_zzzz[i] + 8.0 * tr_yzz_zzz[i] - 4.0 * tr_yzz_zzzzz[i] * tke_0 - 4.0 * tr_yzzz_zzzz[i] * tbe_0 - 4.0 * tr_yyz_yzzzz[i] * tke_0 - 8.0 * tr_yyzz_yzzz[i] * tke_0 + 4.0 * tr_yyzz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_yzzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_zzzz[i] * tbe_0 - 8.0 * tr_yyyzz_zzz[i] * tbe_0 + 4.0 * tr_yyyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1095-1110 components of targeted buffer : GG

    auto tr_0_0_yz_yzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1095);

    auto tr_0_0_yz_yzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1096);

    auto tr_0_0_yz_yzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1097);

    auto tr_0_0_yz_yzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1098);

    auto tr_0_0_yz_yzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1099);

    auto tr_0_0_yz_yzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1100);

    auto tr_0_0_yz_yzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1101);

    auto tr_0_0_yz_yzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1102);

    auto tr_0_0_yz_yzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1103);

    auto tr_0_0_yz_yzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1104);

    auto tr_0_0_yz_yzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1105);

    auto tr_0_0_yz_yzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1106);

    auto tr_0_0_yz_yzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1107);

    auto tr_0_0_yz_yzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1108);

    auto tr_0_0_yz_yzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1109);

    #pragma omp simd aligned(tr_0_0_yz_yzzz_xxxx, tr_0_0_yz_yzzz_xxxy, tr_0_0_yz_yzzz_xxxz, tr_0_0_yz_yzzz_xxyy, tr_0_0_yz_yzzz_xxyz, tr_0_0_yz_yzzz_xxzz, tr_0_0_yz_yzzz_xyyy, tr_0_0_yz_yzzz_xyyz, tr_0_0_yz_yzzz_xyzz, tr_0_0_yz_yzzz_xzzz, tr_0_0_yz_yzzz_yyyy, tr_0_0_yz_yzzz_yyyz, tr_0_0_yz_yzzz_yyzz, tr_0_0_yz_yzzz_yzzz, tr_0_0_yz_yzzz_zzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxxxz, tr_yyzzz_xxxyz, tr_yyzzz_xxxzz, tr_yyzzz_xxy, tr_yyzzz_xxyyz, tr_yyzzz_xxyzz, tr_yyzzz_xxz, tr_yyzzz_xxzzz, tr_yyzzz_xyy, tr_yyzzz_xyyyz, tr_yyzzz_xyyzz, tr_yyzzz_xyz, tr_yyzzz_xyzzz, tr_yyzzz_xzz, tr_yyzzz_xzzzz, tr_yyzzz_yyy, tr_yyzzz_yyyyz, tr_yyzzz_yyyzz, tr_yyzzz_yyz, tr_yyzzz_yyzzz, tr_yyzzz_yzz, tr_yyzzz_yzzzz, tr_yyzzz_zzz, tr_yyzzz_zzzzz, tr_yyzzzz_xxxx, tr_yyzzzz_xxxy, tr_yyzzzz_xxxz, tr_yyzzzz_xxyy, tr_yyzzzz_xxyz, tr_yyzzzz_xxzz, tr_yyzzzz_xyyy, tr_yyzzzz_xyyz, tr_yyzzzz_xyzz, tr_yyzzzz_xzzz, tr_yyzzzz_yyyy, tr_yyzzzz_yyyz, tr_yyzzzz_yyzz, tr_yyzzzz_yzzz, tr_yyzzzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xxxxyz, tr_yzzz_xxxy, tr_yzzz_xxxyyz, tr_yzzz_xxxyzz, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyyyz, tr_yzzz_xxyyzz, tr_yzzz_xxyz, tr_yzzz_xxyzzz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyyyz, tr_yzzz_xyyyzz, tr_yzzz_xyyz, tr_yzzz_xyyzzz, tr_yzzz_xyzz, tr_yzzz_xyzzzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyyyz, tr_yzzz_yyyyzz, tr_yzzz_yyyz, tr_yzzz_yyyzzz, tr_yzzz_yyzz, tr_yzzz_yyzzzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_yzzzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxxxy, tr_yzzzz_xxxyy, tr_yzzzz_xxxyz, tr_yzzzz_xxy, tr_yzzzz_xxyyy, tr_yzzzz_xxyyz, tr_yzzzz_xxyzz, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyyyy, tr_yzzzz_xyyyz, tr_yzzzz_xyyzz, tr_yzzzz_xyz, tr_yzzzz_xyzzz, tr_yzzzz_xzz, tr_yzzzz_yyy, tr_yzzzz_yyyyy, tr_yzzzz_yyyyz, tr_yzzzz_yyyzz, tr_yzzzz_yyz, tr_yzzzz_yyzzz, tr_yzzzz_yzz, tr_yzzzz_yzzzz, tr_yzzzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxxxz, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, tr_zzz_zzzzz, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xzzz, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yzzz, tr_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yzzz_xxxx[i] = 3.0 * tr_zz_xxxx[i] - 2.0 * tr_zzz_xxxxz[i] * tke_0 - 2.0 * tr_zzzz_xxxx[i] * tbe_0 - 6.0 * tr_yzz_xxxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xxxxy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xxxy[i] = 3.0 * tr_zz_xxxy[i] - 2.0 * tr_zzz_xxxyz[i] * tke_0 - 2.0 * tr_zzzz_xxxy[i] * tbe_0 + 3.0 * tr_yzz_xxx[i] - 6.0 * tr_yzz_xxxyy[i] * tke_0 - 2.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xxxz[i] = 3.0 * tr_zz_xxxz[i] + tr_zzz_xxx[i] - 2.0 * tr_zzz_xxxzz[i] * tke_0 - 2.0 * tr_zzzz_xxxz[i] * tbe_0 - 6.0 * tr_yzz_xxxyz[i] * tke_0 - 2.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xxxyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxxz[i] * tbe_0 - 2.0 * tr_yyzzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xxyy[i] = 3.0 * tr_zz_xxyy[i] - 2.0 * tr_zzz_xxyyz[i] * tke_0 - 2.0 * tr_zzzz_xxyy[i] * tbe_0 + 6.0 * tr_yzz_xxy[i] - 6.0 * tr_yzz_xxyyy[i] * tke_0 - 4.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxyy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xxyz[i] = 3.0 * tr_zz_xxyz[i] + tr_zzz_xxy[i] - 2.0 * tr_zzz_xxyzz[i] * tke_0 - 2.0 * tr_zzzz_xxyz[i] * tbe_0 + 3.0 * tr_yzz_xxz[i] - 6.0 * tr_yzz_xxyyz[i] * tke_0 + tr_yzzz_xx[i] - 2.0 * tr_yzzz_xxzz[i] * tke_0 - 2.0 * tr_yzzz_xxyy[i] * tke_0 + 4.0 * tr_yzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxyz[i] * tbe_0 - 2.0 * tr_yyzzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xxzz[i] = 3.0 * tr_zz_xxzz[i] + 2.0 * tr_zzz_xxz[i] - 2.0 * tr_zzz_xxzzz[i] * tke_0 - 2.0 * tr_zzzz_xxzz[i] * tbe_0 - 6.0 * tr_yzz_xxyzz[i] * tke_0 - 4.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xxyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxzz[i] * tbe_0 - 4.0 * tr_yyzzz_xxz[i] * tbe_0 + 4.0 * tr_yyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xyyy[i] = 3.0 * tr_zz_xyyy[i] - 2.0 * tr_zzz_xyyyz[i] * tke_0 - 2.0 * tr_zzzz_xyyy[i] * tbe_0 + 9.0 * tr_yzz_xyy[i] - 6.0 * tr_yzz_xyyyy[i] * tke_0 - 6.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_yzzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xyyy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xyyz[i] = 3.0 * tr_zz_xyyz[i] + tr_zzz_xyy[i] - 2.0 * tr_zzz_xyyzz[i] * tke_0 - 2.0 * tr_zzzz_xyyz[i] * tbe_0 + 6.0 * tr_yzz_xyz[i] - 6.0 * tr_yzz_xyyyz[i] * tke_0 + 2.0 * tr_yzzz_xy[i] - 4.0 * tr_yzzz_xyzz[i] * tke_0 - 2.0 * tr_yzzz_xyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzzz_xyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xyyz[i] * tbe_0 - 2.0 * tr_yyzzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xyzz[i] = 3.0 * tr_zz_xyzz[i] + 2.0 * tr_zzz_xyz[i] - 2.0 * tr_zzz_xyzzz[i] * tke_0 - 2.0 * tr_zzzz_xyzz[i] * tbe_0 + 3.0 * tr_yzz_xzz[i] - 6.0 * tr_yzz_xyyzz[i] * tke_0 + 2.0 * tr_yzzz_xz[i] - 2.0 * tr_yzzz_xzzz[i] * tke_0 - 4.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzzz_xyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xyzz[i] * tbe_0 - 4.0 * tr_yyzzz_xyz[i] * tbe_0 + 4.0 * tr_yyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xzzz[i] = 3.0 * tr_zz_xzzz[i] + 3.0 * tr_zzz_xzz[i] - 2.0 * tr_zzz_xzzzz[i] * tke_0 - 2.0 * tr_zzzz_xzzz[i] * tbe_0 - 6.0 * tr_yzz_xyzzz[i] * tke_0 - 6.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xzzz[i] * tbe_0 - 6.0 * tr_yyzzz_xzz[i] * tbe_0 + 4.0 * tr_yyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_yyyy[i] = 3.0 * tr_zz_yyyy[i] - 2.0 * tr_zzz_yyyyz[i] * tke_0 - 2.0 * tr_zzzz_yyyy[i] * tbe_0 + 12.0 * tr_yzz_yyy[i] - 6.0 * tr_yzz_yyyyy[i] * tke_0 - 8.0 * tr_yzzz_yyyz[i] * tke_0 + 4.0 * tr_yzzz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yyyy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_yyyz[i] = 3.0 * tr_zz_yyyz[i] + tr_zzz_yyy[i] - 2.0 * tr_zzz_yyyzz[i] * tke_0 - 2.0 * tr_zzzz_yyyz[i] * tbe_0 + 9.0 * tr_yzz_yyz[i] - 6.0 * tr_yzz_yyyyz[i] * tke_0 + 3.0 * tr_yzzz_yy[i] - 6.0 * tr_yzzz_yyzz[i] * tke_0 - 2.0 * tr_yzzz_yyyy[i] * tke_0 + 4.0 * tr_yzzz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_yzzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzzz_yyyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yyyz[i] * tbe_0 - 2.0 * tr_yyzzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_yyzz[i] = 3.0 * tr_zz_yyzz[i] + 2.0 * tr_zzz_yyz[i] - 2.0 * tr_zzz_yyzzz[i] * tke_0 - 2.0 * tr_zzzz_yyzz[i] * tbe_0 + 6.0 * tr_yzz_yzz[i] - 6.0 * tr_yzz_yyyzz[i] * tke_0 + 4.0 * tr_yzzz_yz[i] - 4.0 * tr_yzzz_yzzz[i] * tke_0 - 4.0 * tr_yzzz_yyyz[i] * tke_0 + 4.0 * tr_yzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzzz_yyyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yyzz[i] * tbe_0 - 4.0 * tr_yyzzz_yyz[i] * tbe_0 + 4.0 * tr_yyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_yzzz[i] = 3.0 * tr_zz_yzzz[i] + 3.0 * tr_zzz_yzz[i] - 2.0 * tr_zzz_yzzzz[i] * tke_0 - 2.0 * tr_zzzz_yzzz[i] * tbe_0 + 3.0 * tr_yzz_zzz[i] - 6.0 * tr_yzz_yyzzz[i] * tke_0 + 3.0 * tr_yzzz_zz[i] - 2.0 * tr_yzzz_zzzz[i] * tke_0 - 6.0 * tr_yzzz_yyzz[i] * tke_0 + 4.0 * tr_yzzz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzzz_yyzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yzzz[i] * tbe_0 - 6.0 * tr_yyzzz_yzz[i] * tbe_0 + 4.0 * tr_yyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_zzzz[i] = 3.0 * tr_zz_zzzz[i] + 4.0 * tr_zzz_zzz[i] - 2.0 * tr_zzz_zzzzz[i] * tke_0 - 2.0 * tr_zzzz_zzzz[i] * tbe_0 - 6.0 * tr_yzz_yzzzz[i] * tke_0 - 8.0 * tr_yzzz_yzzz[i] * tke_0 + 4.0 * tr_yzzz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_yzzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_zzzz[i] * tbe_0 - 8.0 * tr_yyzzz_zzz[i] * tbe_0 + 4.0 * tr_yyzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1110-1125 components of targeted buffer : GG

    auto tr_0_0_yz_zzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1110);

    auto tr_0_0_yz_zzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1111);

    auto tr_0_0_yz_zzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1112);

    auto tr_0_0_yz_zzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1113);

    auto tr_0_0_yz_zzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1114);

    auto tr_0_0_yz_zzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1115);

    auto tr_0_0_yz_zzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1116);

    auto tr_0_0_yz_zzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1117);

    auto tr_0_0_yz_zzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1118);

    auto tr_0_0_yz_zzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1119);

    auto tr_0_0_yz_zzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1120);

    auto tr_0_0_yz_zzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1121);

    auto tr_0_0_yz_zzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1122);

    auto tr_0_0_yz_zzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1123);

    auto tr_0_0_yz_zzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1124);

    #pragma omp simd aligned(tr_0_0_yz_zzzz_xxxx, tr_0_0_yz_zzzz_xxxy, tr_0_0_yz_zzzz_xxxz, tr_0_0_yz_zzzz_xxyy, tr_0_0_yz_zzzz_xxyz, tr_0_0_yz_zzzz_xxzz, tr_0_0_yz_zzzz_xyyy, tr_0_0_yz_zzzz_xyyz, tr_0_0_yz_zzzz_xyzz, tr_0_0_yz_zzzz_xzzz, tr_0_0_yz_zzzz_yyyy, tr_0_0_yz_zzzz_yyyz, tr_0_0_yz_zzzz_yyzz, tr_0_0_yz_zzzz_yzzz, tr_0_0_yz_zzzz_zzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxxxz, tr_yzzzz_xxxyz, tr_yzzzz_xxxzz, tr_yzzzz_xxy, tr_yzzzz_xxyyz, tr_yzzzz_xxyzz, tr_yzzzz_xxz, tr_yzzzz_xxzzz, tr_yzzzz_xyy, tr_yzzzz_xyyyz, tr_yzzzz_xyyzz, tr_yzzzz_xyz, tr_yzzzz_xyzzz, tr_yzzzz_xzz, tr_yzzzz_xzzzz, tr_yzzzz_yyy, tr_yzzzz_yyyyz, tr_yzzzz_yyyzz, tr_yzzzz_yyz, tr_yzzzz_yyzzz, tr_yzzzz_yzz, tr_yzzzz_yzzzz, tr_yzzzz_zzz, tr_yzzzz_zzzzz, tr_yzzzzz_xxxx, tr_yzzzzz_xxxy, tr_yzzzzz_xxxz, tr_yzzzzz_xxyy, tr_yzzzzz_xxyz, tr_yzzzzz_xxzz, tr_yzzzzz_xyyy, tr_yzzzzz_xyyz, tr_yzzzzz_xyzz, tr_yzzzzz_xzzz, tr_yzzzzz_yyyy, tr_yzzzzz_yyyz, tr_yzzzzz_yyzz, tr_yzzzzz_yzzz, tr_yzzzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxy, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyyyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xxxxyz, tr_zzzz_xxxy, tr_zzzz_xxxyyz, tr_zzzz_xxxyzz, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyyyz, tr_zzzz_xxyyzz, tr_zzzz_xxyz, tr_zzzz_xxyzzz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyyyz, tr_zzzz_xyyyzz, tr_zzzz_xyyz, tr_zzzz_xyyzzz, tr_zzzz_xyzz, tr_zzzz_xyzzzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyyyz, tr_zzzz_yyyyzz, tr_zzzz_yyyz, tr_zzzz_yyyzzz, tr_zzzz_yyzz, tr_zzzz_yyzzzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_yzzzzz, tr_zzzz_zz, tr_zzzz_zzzz, tr_zzzzz_xxx, tr_zzzzz_xxxxy, tr_zzzzz_xxxyy, tr_zzzzz_xxxyz, tr_zzzzz_xxy, tr_zzzzz_xxyyy, tr_zzzzz_xxyyz, tr_zzzzz_xxyzz, tr_zzzzz_xxz, tr_zzzzz_xyy, tr_zzzzz_xyyyy, tr_zzzzz_xyyyz, tr_zzzzz_xyyzz, tr_zzzzz_xyz, tr_zzzzz_xyzzz, tr_zzzzz_xzz, tr_zzzzz_yyy, tr_zzzzz_yyyyy, tr_zzzzz_yyyyz, tr_zzzzz_yyyzz, tr_zzzzz_yyz, tr_zzzzz_yyzzz, tr_zzzzz_yzz, tr_zzzzz_yzzzz, tr_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_zzzz_xxxx[i] = -8.0 * tr_zzz_xxxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxxyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xxxxy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xxxy[i] = 4.0 * tr_zzz_xxx[i] - 8.0 * tr_zzz_xxxyy[i] * tke_0 - 2.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_xxx[i] * tbe_0 + 4.0 * tr_zzzzz_xxxyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xxxz[i] = -8.0 * tr_zzz_xxxyz[i] * tke_0 - 2.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xxxyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxxz[i] * tbe_0 - 2.0 * tr_yzzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xxyy[i] = 8.0 * tr_zzz_xxy[i] - 8.0 * tr_zzz_xxyyy[i] * tke_0 - 4.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxyyyz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xxy[i] * tbe_0 + 4.0 * tr_zzzzz_xxyyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxyy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xxyz[i] = 4.0 * tr_zzz_xxz[i] - 8.0 * tr_zzz_xxyyz[i] * tke_0 + tr_zzzz_xx[i] - 2.0 * tr_zzzz_xxzz[i] * tke_0 - 2.0 * tr_zzzz_xxyy[i] * tke_0 + 4.0 * tr_zzzz_xxyyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_xxz[i] * tbe_0 + 4.0 * tr_zzzzz_xxyyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxyz[i] * tbe_0 - 2.0 * tr_yzzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xxzz[i] = -8.0 * tr_zzz_xxyzz[i] * tke_0 - 4.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxyzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xxyzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxzz[i] * tbe_0 - 4.0 * tr_yzzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xyyy[i] = 12.0 * tr_zzz_xyy[i] - 8.0 * tr_zzz_xyyyy[i] * tke_0 - 6.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyyyz[i] * tke_0 * tke_0 - 6.0 * tr_zzzzz_xyy[i] * tbe_0 + 4.0 * tr_zzzzz_xyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xyyy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xyyz[i] = 8.0 * tr_zzz_xyz[i] - 8.0 * tr_zzz_xyyyz[i] * tke_0 + 2.0 * tr_zzzz_xy[i] - 4.0 * tr_zzzz_xyzz[i] * tke_0 - 2.0 * tr_zzzz_xyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xyz[i] * tbe_0 + 4.0 * tr_zzzzz_xyyyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xyyz[i] * tbe_0 - 2.0 * tr_yzzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xyzz[i] = 4.0 * tr_zzz_xzz[i] - 8.0 * tr_zzz_xyyzz[i] * tke_0 + 2.0 * tr_zzzz_xz[i] - 2.0 * tr_zzzz_xzzz[i] * tke_0 - 4.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_xzz[i] * tbe_0 + 4.0 * tr_zzzzz_xyyzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xyzz[i] * tbe_0 - 4.0 * tr_yzzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xzzz[i] = -8.0 * tr_zzz_xyzzz[i] * tke_0 - 6.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xyzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyzzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xzzz[i] * tbe_0 - 6.0 * tr_yzzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_yyyy[i] = 16.0 * tr_zzz_yyy[i] - 8.0 * tr_zzz_yyyyy[i] * tke_0 - 8.0 * tr_zzzz_yyyz[i] * tke_0 + 4.0 * tr_zzzz_yyyyyz[i] * tke_0 * tke_0 - 8.0 * tr_zzzzz_yyy[i] * tbe_0 + 4.0 * tr_zzzzz_yyyyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yyyy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_yyyz[i] = 12.0 * tr_zzz_yyz[i] - 8.0 * tr_zzz_yyyyz[i] * tke_0 + 3.0 * tr_zzzz_yy[i] - 6.0 * tr_zzzz_yyzz[i] * tke_0 - 2.0 * tr_zzzz_yyyy[i] * tke_0 + 4.0 * tr_zzzz_yyyyzz[i] * tke_0 * tke_0 - 6.0 * tr_zzzzz_yyz[i] * tbe_0 + 4.0 * tr_zzzzz_yyyyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yyyz[i] * tbe_0 - 2.0 * tr_yzzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_yyzz[i] = 8.0 * tr_zzz_yzz[i] - 8.0 * tr_zzz_yyyzz[i] * tke_0 + 4.0 * tr_zzzz_yz[i] - 4.0 * tr_zzzz_yzzz[i] * tke_0 - 4.0 * tr_zzzz_yyyz[i] * tke_0 + 4.0 * tr_zzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_yzz[i] * tbe_0 + 4.0 * tr_zzzzz_yyyzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yyzz[i] * tbe_0 - 4.0 * tr_yzzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_yzzz[i] = 4.0 * tr_zzz_zzz[i] - 8.0 * tr_zzz_yyzzz[i] * tke_0 + 3.0 * tr_zzzz_zz[i] - 2.0 * tr_zzzz_zzzz[i] * tke_0 - 6.0 * tr_zzzz_yyzz[i] * tke_0 + 4.0 * tr_zzzz_yyzzzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_zzz[i] * tbe_0 + 4.0 * tr_zzzzz_yyzzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yzzz[i] * tbe_0 - 6.0 * tr_yzzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_zzzz[i] = -8.0 * tr_zzz_yzzzz[i] * tke_0 - 8.0 * tr_zzzz_yzzz[i] * tke_0 + 4.0 * tr_zzzz_yzzzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_yzzzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_zzzz[i] * tbe_0 - 8.0 * tr_yzzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1125-1140 components of targeted buffer : GG

    auto tr_0_0_zz_xxxx_xxxx = pbuffer.data(idx_op_geom_020_gg + 1125);

    auto tr_0_0_zz_xxxx_xxxy = pbuffer.data(idx_op_geom_020_gg + 1126);

    auto tr_0_0_zz_xxxx_xxxz = pbuffer.data(idx_op_geom_020_gg + 1127);

    auto tr_0_0_zz_xxxx_xxyy = pbuffer.data(idx_op_geom_020_gg + 1128);

    auto tr_0_0_zz_xxxx_xxyz = pbuffer.data(idx_op_geom_020_gg + 1129);

    auto tr_0_0_zz_xxxx_xxzz = pbuffer.data(idx_op_geom_020_gg + 1130);

    auto tr_0_0_zz_xxxx_xyyy = pbuffer.data(idx_op_geom_020_gg + 1131);

    auto tr_0_0_zz_xxxx_xyyz = pbuffer.data(idx_op_geom_020_gg + 1132);

    auto tr_0_0_zz_xxxx_xyzz = pbuffer.data(idx_op_geom_020_gg + 1133);

    auto tr_0_0_zz_xxxx_xzzz = pbuffer.data(idx_op_geom_020_gg + 1134);

    auto tr_0_0_zz_xxxx_yyyy = pbuffer.data(idx_op_geom_020_gg + 1135);

    auto tr_0_0_zz_xxxx_yyyz = pbuffer.data(idx_op_geom_020_gg + 1136);

    auto tr_0_0_zz_xxxx_yyzz = pbuffer.data(idx_op_geom_020_gg + 1137);

    auto tr_0_0_zz_xxxx_yzzz = pbuffer.data(idx_op_geom_020_gg + 1138);

    auto tr_0_0_zz_xxxx_zzzz = pbuffer.data(idx_op_geom_020_gg + 1139);

    #pragma omp simd aligned(tr_0_0_zz_xxxx_xxxx, tr_0_0_zz_xxxx_xxxy, tr_0_0_zz_xxxx_xxxz, tr_0_0_zz_xxxx_xxyy, tr_0_0_zz_xxxx_xxyz, tr_0_0_zz_xxxx_xxzz, tr_0_0_zz_xxxx_xyyy, tr_0_0_zz_xxxx_xyyz, tr_0_0_zz_xxxx_xyzz, tr_0_0_zz_xxxx_xzzz, tr_0_0_zz_xxxx_yyyy, tr_0_0_zz_xxxx_yyyz, tr_0_0_zz_xxxx_yyzz, tr_0_0_zz_xxxx_yzzz, tr_0_0_zz_xxxx_zzzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxxzz, tr_xxxx_xxxy, tr_xxxx_xxxyzz, tr_xxxx_xxxz, tr_xxxx_xxxzzz, tr_xxxx_xxyy, tr_xxxx_xxyyzz, tr_xxxx_xxyz, tr_xxxx_xxyzzz, tr_xxxx_xxzz, tr_xxxx_xxzzzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyyzz, tr_xxxx_xyyz, tr_xxxx_xyyzzz, tr_xxxx_xyzz, tr_xxxx_xyzzzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_xzzzzz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyyzz, tr_xxxx_yyyz, tr_xxxx_yyyzzz, tr_xxxx_yyzz, tr_xxxx_yyzzzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_yzzzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxx_zzzzzz, tr_xxxxz_xxx, tr_xxxxz_xxxxz, tr_xxxxz_xxxyz, tr_xxxxz_xxxzz, tr_xxxxz_xxy, tr_xxxxz_xxyyz, tr_xxxxz_xxyzz, tr_xxxxz_xxz, tr_xxxxz_xxzzz, tr_xxxxz_xyy, tr_xxxxz_xyyyz, tr_xxxxz_xyyzz, tr_xxxxz_xyz, tr_xxxxz_xyzzz, tr_xxxxz_xzz, tr_xxxxz_xzzzz, tr_xxxxz_yyy, tr_xxxxz_yyyyz, tr_xxxxz_yyyzz, tr_xxxxz_yyz, tr_xxxxz_yyzzz, tr_xxxxz_yzz, tr_xxxxz_yzzzz, tr_xxxxz_zzz, tr_xxxxz_zzzzz, tr_xxxxzz_xxxx, tr_xxxxzz_xxxy, tr_xxxxzz_xxxz, tr_xxxxzz_xxyy, tr_xxxxzz_xxyz, tr_xxxxzz_xxzz, tr_xxxxzz_xyyy, tr_xxxxzz_xyyz, tr_xxxxzz_xyzz, tr_xxxxzz_xzzz, tr_xxxxzz_yyyy, tr_xxxxzz_yyyz, tr_xxxxzz_yyzz, tr_xxxxzz_yzzz, tr_xxxxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxxx_xxxx[i] = -2.0 * tr_xxxx_xxxx[i] * tbe_0 - 2.0 * tr_xxxx_xxxx[i] * tke_0 + 4.0 * tr_xxxx_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xxxy[i] = -2.0 * tr_xxxx_xxxy[i] * tbe_0 - 2.0 * tr_xxxx_xxxy[i] * tke_0 + 4.0 * tr_xxxx_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xxxz[i] = -2.0 * tr_xxxx_xxxz[i] * tbe_0 - 6.0 * tr_xxxx_xxxz[i] * tke_0 + 4.0 * tr_xxxx_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xxx[i] * tbe_0 + 8.0 * tr_xxxxz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xxyy[i] = -2.0 * tr_xxxx_xxyy[i] * tbe_0 - 2.0 * tr_xxxx_xxyy[i] * tke_0 + 4.0 * tr_xxxx_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xxyz[i] = -2.0 * tr_xxxx_xxyz[i] * tbe_0 - 6.0 * tr_xxxx_xxyz[i] * tke_0 + 4.0 * tr_xxxx_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xxy[i] * tbe_0 + 8.0 * tr_xxxxz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xxzz[i] = 2.0 * tr_xxxx_xx[i] - 2.0 * tr_xxxx_xxzz[i] * tbe_0 - 10.0 * tr_xxxx_xxzz[i] * tke_0 + 4.0 * tr_xxxx_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_xxz[i] * tbe_0 + 8.0 * tr_xxxxz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xyyy[i] = -2.0 * tr_xxxx_xyyy[i] * tbe_0 - 2.0 * tr_xxxx_xyyy[i] * tke_0 + 4.0 * tr_xxxx_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xyyz[i] = -2.0 * tr_xxxx_xyyz[i] * tbe_0 - 6.0 * tr_xxxx_xyyz[i] * tke_0 + 4.0 * tr_xxxx_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_xyy[i] * tbe_0 + 8.0 * tr_xxxxz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xyzz[i] = 2.0 * tr_xxxx_xy[i] - 2.0 * tr_xxxx_xyzz[i] * tbe_0 - 10.0 * tr_xxxx_xyzz[i] * tke_0 + 4.0 * tr_xxxx_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_xyz[i] * tbe_0 + 8.0 * tr_xxxxz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xzzz[i] = 6.0 * tr_xxxx_xz[i] - 2.0 * tr_xxxx_xzzz[i] * tbe_0 - 14.0 * tr_xxxx_xzzz[i] * tke_0 + 4.0 * tr_xxxx_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxxz_xzz[i] * tbe_0 + 8.0 * tr_xxxxz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_yyyy[i] = -2.0 * tr_xxxx_yyyy[i] * tbe_0 - 2.0 * tr_xxxx_yyyy[i] * tke_0 + 4.0 * tr_xxxx_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_yyyz[i] = -2.0 * tr_xxxx_yyyz[i] * tbe_0 - 6.0 * tr_xxxx_yyyz[i] * tke_0 + 4.0 * tr_xxxx_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_yyy[i] * tbe_0 + 8.0 * tr_xxxxz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_yyzz[i] = 2.0 * tr_xxxx_yy[i] - 2.0 * tr_xxxx_yyzz[i] * tbe_0 - 10.0 * tr_xxxx_yyzz[i] * tke_0 + 4.0 * tr_xxxx_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_yyz[i] * tbe_0 + 8.0 * tr_xxxxz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_yzzz[i] = 6.0 * tr_xxxx_yz[i] - 2.0 * tr_xxxx_yzzz[i] * tbe_0 - 14.0 * tr_xxxx_yzzz[i] * tke_0 + 4.0 * tr_xxxx_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxxz_yzz[i] * tbe_0 + 8.0 * tr_xxxxz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_zzzz[i] = 12.0 * tr_xxxx_zz[i] - 2.0 * tr_xxxx_zzzz[i] * tbe_0 - 18.0 * tr_xxxx_zzzz[i] * tke_0 + 4.0 * tr_xxxx_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xxxxz_zzz[i] * tbe_0 + 8.0 * tr_xxxxz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1140-1155 components of targeted buffer : GG

    auto tr_0_0_zz_xxxy_xxxx = pbuffer.data(idx_op_geom_020_gg + 1140);

    auto tr_0_0_zz_xxxy_xxxy = pbuffer.data(idx_op_geom_020_gg + 1141);

    auto tr_0_0_zz_xxxy_xxxz = pbuffer.data(idx_op_geom_020_gg + 1142);

    auto tr_0_0_zz_xxxy_xxyy = pbuffer.data(idx_op_geom_020_gg + 1143);

    auto tr_0_0_zz_xxxy_xxyz = pbuffer.data(idx_op_geom_020_gg + 1144);

    auto tr_0_0_zz_xxxy_xxzz = pbuffer.data(idx_op_geom_020_gg + 1145);

    auto tr_0_0_zz_xxxy_xyyy = pbuffer.data(idx_op_geom_020_gg + 1146);

    auto tr_0_0_zz_xxxy_xyyz = pbuffer.data(idx_op_geom_020_gg + 1147);

    auto tr_0_0_zz_xxxy_xyzz = pbuffer.data(idx_op_geom_020_gg + 1148);

    auto tr_0_0_zz_xxxy_xzzz = pbuffer.data(idx_op_geom_020_gg + 1149);

    auto tr_0_0_zz_xxxy_yyyy = pbuffer.data(idx_op_geom_020_gg + 1150);

    auto tr_0_0_zz_xxxy_yyyz = pbuffer.data(idx_op_geom_020_gg + 1151);

    auto tr_0_0_zz_xxxy_yyzz = pbuffer.data(idx_op_geom_020_gg + 1152);

    auto tr_0_0_zz_xxxy_yzzz = pbuffer.data(idx_op_geom_020_gg + 1153);

    auto tr_0_0_zz_xxxy_zzzz = pbuffer.data(idx_op_geom_020_gg + 1154);

    #pragma omp simd aligned(tr_0_0_zz_xxxy_xxxx, tr_0_0_zz_xxxy_xxxy, tr_0_0_zz_xxxy_xxxz, tr_0_0_zz_xxxy_xxyy, tr_0_0_zz_xxxy_xxyz, tr_0_0_zz_xxxy_xxzz, tr_0_0_zz_xxxy_xyyy, tr_0_0_zz_xxxy_xyyz, tr_0_0_zz_xxxy_xyzz, tr_0_0_zz_xxxy_xzzz, tr_0_0_zz_xxxy_yyyy, tr_0_0_zz_xxxy_yyyz, tr_0_0_zz_xxxy_yyzz, tr_0_0_zz_xxxy_yzzz, tr_0_0_zz_xxxy_zzzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxxzz, tr_xxxy_xxxy, tr_xxxy_xxxyzz, tr_xxxy_xxxz, tr_xxxy_xxxzzz, tr_xxxy_xxyy, tr_xxxy_xxyyzz, tr_xxxy_xxyz, tr_xxxy_xxyzzz, tr_xxxy_xxzz, tr_xxxy_xxzzzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyyzz, tr_xxxy_xyyz, tr_xxxy_xyyzzz, tr_xxxy_xyzz, tr_xxxy_xyzzzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_xzzzzz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyyzz, tr_xxxy_yyyz, tr_xxxy_yyyzzz, tr_xxxy_yyzz, tr_xxxy_yyzzzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_yzzzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxy_zzzzzz, tr_xxxyz_xxx, tr_xxxyz_xxxxz, tr_xxxyz_xxxyz, tr_xxxyz_xxxzz, tr_xxxyz_xxy, tr_xxxyz_xxyyz, tr_xxxyz_xxyzz, tr_xxxyz_xxz, tr_xxxyz_xxzzz, tr_xxxyz_xyy, tr_xxxyz_xyyyz, tr_xxxyz_xyyzz, tr_xxxyz_xyz, tr_xxxyz_xyzzz, tr_xxxyz_xzz, tr_xxxyz_xzzzz, tr_xxxyz_yyy, tr_xxxyz_yyyyz, tr_xxxyz_yyyzz, tr_xxxyz_yyz, tr_xxxyz_yyzzz, tr_xxxyz_yzz, tr_xxxyz_yzzzz, tr_xxxyz_zzz, tr_xxxyz_zzzzz, tr_xxxyzz_xxxx, tr_xxxyzz_xxxy, tr_xxxyzz_xxxz, tr_xxxyzz_xxyy, tr_xxxyzz_xxyz, tr_xxxyzz_xxzz, tr_xxxyzz_xyyy, tr_xxxyzz_xyyz, tr_xxxyzz_xyzz, tr_xxxyzz_xzzz, tr_xxxyzz_yyyy, tr_xxxyzz_yyyz, tr_xxxyzz_yyzz, tr_xxxyzz_yzzz, tr_xxxyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxxy_xxxx[i] = -2.0 * tr_xxxy_xxxx[i] * tbe_0 - 2.0 * tr_xxxy_xxxx[i] * tke_0 + 4.0 * tr_xxxy_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xxxy[i] = -2.0 * tr_xxxy_xxxy[i] * tbe_0 - 2.0 * tr_xxxy_xxxy[i] * tke_0 + 4.0 * tr_xxxy_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xxxz[i] = -2.0 * tr_xxxy_xxxz[i] * tbe_0 - 6.0 * tr_xxxy_xxxz[i] * tke_0 + 4.0 * tr_xxxy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xxx[i] * tbe_0 + 8.0 * tr_xxxyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xxyy[i] = -2.0 * tr_xxxy_xxyy[i] * tbe_0 - 2.0 * tr_xxxy_xxyy[i] * tke_0 + 4.0 * tr_xxxy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xxyz[i] = -2.0 * tr_xxxy_xxyz[i] * tbe_0 - 6.0 * tr_xxxy_xxyz[i] * tke_0 + 4.0 * tr_xxxy_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xxy[i] * tbe_0 + 8.0 * tr_xxxyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xxzz[i] = 2.0 * tr_xxxy_xx[i] - 2.0 * tr_xxxy_xxzz[i] * tbe_0 - 10.0 * tr_xxxy_xxzz[i] * tke_0 + 4.0 * tr_xxxy_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xxz[i] * tbe_0 + 8.0 * tr_xxxyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xyyy[i] = -2.0 * tr_xxxy_xyyy[i] * tbe_0 - 2.0 * tr_xxxy_xyyy[i] * tke_0 + 4.0 * tr_xxxy_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xyyz[i] = -2.0 * tr_xxxy_xyyz[i] * tbe_0 - 6.0 * tr_xxxy_xyyz[i] * tke_0 + 4.0 * tr_xxxy_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_xyy[i] * tbe_0 + 8.0 * tr_xxxyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xyzz[i] = 2.0 * tr_xxxy_xy[i] - 2.0 * tr_xxxy_xyzz[i] * tbe_0 - 10.0 * tr_xxxy_xyzz[i] * tke_0 + 4.0 * tr_xxxy_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_xyz[i] * tbe_0 + 8.0 * tr_xxxyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xzzz[i] = 6.0 * tr_xxxy_xz[i] - 2.0 * tr_xxxy_xzzz[i] * tbe_0 - 14.0 * tr_xxxy_xzzz[i] * tke_0 + 4.0 * tr_xxxy_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxyz_xzz[i] * tbe_0 + 8.0 * tr_xxxyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_yyyy[i] = -2.0 * tr_xxxy_yyyy[i] * tbe_0 - 2.0 * tr_xxxy_yyyy[i] * tke_0 + 4.0 * tr_xxxy_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_yyyz[i] = -2.0 * tr_xxxy_yyyz[i] * tbe_0 - 6.0 * tr_xxxy_yyyz[i] * tke_0 + 4.0 * tr_xxxy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_yyy[i] * tbe_0 + 8.0 * tr_xxxyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_yyzz[i] = 2.0 * tr_xxxy_yy[i] - 2.0 * tr_xxxy_yyzz[i] * tbe_0 - 10.0 * tr_xxxy_yyzz[i] * tke_0 + 4.0 * tr_xxxy_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_yyz[i] * tbe_0 + 8.0 * tr_xxxyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_yzzz[i] = 6.0 * tr_xxxy_yz[i] - 2.0 * tr_xxxy_yzzz[i] * tbe_0 - 14.0 * tr_xxxy_yzzz[i] * tke_0 + 4.0 * tr_xxxy_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxyz_yzz[i] * tbe_0 + 8.0 * tr_xxxyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_zzzz[i] = 12.0 * tr_xxxy_zz[i] - 2.0 * tr_xxxy_zzzz[i] * tbe_0 - 18.0 * tr_xxxy_zzzz[i] * tke_0 + 4.0 * tr_xxxy_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xxxyz_zzz[i] * tbe_0 + 8.0 * tr_xxxyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1155-1170 components of targeted buffer : GG

    auto tr_0_0_zz_xxxz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1155);

    auto tr_0_0_zz_xxxz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1156);

    auto tr_0_0_zz_xxxz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1157);

    auto tr_0_0_zz_xxxz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1158);

    auto tr_0_0_zz_xxxz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1159);

    auto tr_0_0_zz_xxxz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1160);

    auto tr_0_0_zz_xxxz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1161);

    auto tr_0_0_zz_xxxz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1162);

    auto tr_0_0_zz_xxxz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1163);

    auto tr_0_0_zz_xxxz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1164);

    auto tr_0_0_zz_xxxz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1165);

    auto tr_0_0_zz_xxxz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1166);

    auto tr_0_0_zz_xxxz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1167);

    auto tr_0_0_zz_xxxz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1168);

    auto tr_0_0_zz_xxxz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1169);

    #pragma omp simd aligned(tr_0_0_zz_xxxz_xxxx, tr_0_0_zz_xxxz_xxxy, tr_0_0_zz_xxxz_xxxz, tr_0_0_zz_xxxz_xxyy, tr_0_0_zz_xxxz_xxyz, tr_0_0_zz_xxxz_xxzz, tr_0_0_zz_xxxz_xyyy, tr_0_0_zz_xxxz_xyyz, tr_0_0_zz_xxxz_xyzz, tr_0_0_zz_xxxz_xzzz, tr_0_0_zz_xxxz_yyyy, tr_0_0_zz_xxxz_yyyz, tr_0_0_zz_xxxz_yyzz, tr_0_0_zz_xxxz_yzzz, tr_0_0_zz_xxxz_zzzz, tr_xxx_xxx, tr_xxx_xxxxz, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxx_zzzzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxxzz, tr_xxxz_xxxy, tr_xxxz_xxxyzz, tr_xxxz_xxxz, tr_xxxz_xxxzzz, tr_xxxz_xxyy, tr_xxxz_xxyyzz, tr_xxxz_xxyz, tr_xxxz_xxyzzz, tr_xxxz_xxzz, tr_xxxz_xxzzzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyyzz, tr_xxxz_xyyz, tr_xxxz_xyyzzz, tr_xxxz_xyzz, tr_xxxz_xyzzzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_xzzzzz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyyzz, tr_xxxz_yyyz, tr_xxxz_yyyzzz, tr_xxxz_yyzz, tr_xxxz_yyzzzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_yzzzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxz_zzzzzz, tr_xxxzz_xxx, tr_xxxzz_xxxxz, tr_xxxzz_xxxyz, tr_xxxzz_xxxzz, tr_xxxzz_xxy, tr_xxxzz_xxyyz, tr_xxxzz_xxyzz, tr_xxxzz_xxz, tr_xxxzz_xxzzz, tr_xxxzz_xyy, tr_xxxzz_xyyyz, tr_xxxzz_xyyzz, tr_xxxzz_xyz, tr_xxxzz_xyzzz, tr_xxxzz_xzz, tr_xxxzz_xzzzz, tr_xxxzz_yyy, tr_xxxzz_yyyyz, tr_xxxzz_yyyzz, tr_xxxzz_yyz, tr_xxxzz_yyzzz, tr_xxxzz_yzz, tr_xxxzz_yzzzz, tr_xxxzz_zzz, tr_xxxzz_zzzzz, tr_xxxzzz_xxxx, tr_xxxzzz_xxxy, tr_xxxzzz_xxxz, tr_xxxzzz_xxyy, tr_xxxzzz_xxyz, tr_xxxzzz_xxzz, tr_xxxzzz_xyyy, tr_xxxzzz_xyyz, tr_xxxzzz_xyzz, tr_xxxzzz_xzzz, tr_xxxzzz_yyyy, tr_xxxzzz_yyyz, tr_xxxzzz_yyzz, tr_xxxzzz_yzzz, tr_xxxzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxxz_xxxx[i] = -4.0 * tr_xxx_xxxxz[i] * tke_0 - 6.0 * tr_xxxz_xxxx[i] * tbe_0 - 2.0 * tr_xxxz_xxxx[i] * tke_0 + 4.0 * tr_xxxz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xxxy[i] = -4.0 * tr_xxx_xxxyz[i] * tke_0 - 6.0 * tr_xxxz_xxxy[i] * tbe_0 - 2.0 * tr_xxxz_xxxy[i] * tke_0 + 4.0 * tr_xxxz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xxxz[i] = 2.0 * tr_xxx_xxx[i] - 4.0 * tr_xxx_xxxzz[i] * tke_0 - 6.0 * tr_xxxz_xxxz[i] * tbe_0 - 6.0 * tr_xxxz_xxxz[i] * tke_0 + 4.0 * tr_xxxz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xxx[i] * tbe_0 + 8.0 * tr_xxxzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xxyy[i] = -4.0 * tr_xxx_xxyyz[i] * tke_0 - 6.0 * tr_xxxz_xxyy[i] * tbe_0 - 2.0 * tr_xxxz_xxyy[i] * tke_0 + 4.0 * tr_xxxz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xxyz[i] = 2.0 * tr_xxx_xxy[i] - 4.0 * tr_xxx_xxyzz[i] * tke_0 - 6.0 * tr_xxxz_xxyz[i] * tbe_0 - 6.0 * tr_xxxz_xxyz[i] * tke_0 + 4.0 * tr_xxxz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xxy[i] * tbe_0 + 8.0 * tr_xxxzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xxzz[i] = 4.0 * tr_xxx_xxz[i] - 4.0 * tr_xxx_xxzzz[i] * tke_0 + 2.0 * tr_xxxz_xx[i] - 6.0 * tr_xxxz_xxzz[i] * tbe_0 - 10.0 * tr_xxxz_xxzz[i] * tke_0 + 4.0 * tr_xxxz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_xxz[i] * tbe_0 + 8.0 * tr_xxxzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xyyy[i] = -4.0 * tr_xxx_xyyyz[i] * tke_0 - 6.0 * tr_xxxz_xyyy[i] * tbe_0 - 2.0 * tr_xxxz_xyyy[i] * tke_0 + 4.0 * tr_xxxz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xyyz[i] = 2.0 * tr_xxx_xyy[i] - 4.0 * tr_xxx_xyyzz[i] * tke_0 - 6.0 * tr_xxxz_xyyz[i] * tbe_0 - 6.0 * tr_xxxz_xyyz[i] * tke_0 + 4.0 * tr_xxxz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_xyy[i] * tbe_0 + 8.0 * tr_xxxzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xyzz[i] = 4.0 * tr_xxx_xyz[i] - 4.0 * tr_xxx_xyzzz[i] * tke_0 + 2.0 * tr_xxxz_xy[i] - 6.0 * tr_xxxz_xyzz[i] * tbe_0 - 10.0 * tr_xxxz_xyzz[i] * tke_0 + 4.0 * tr_xxxz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_xyz[i] * tbe_0 + 8.0 * tr_xxxzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xzzz[i] = 6.0 * tr_xxx_xzz[i] - 4.0 * tr_xxx_xzzzz[i] * tke_0 + 6.0 * tr_xxxz_xz[i] - 6.0 * tr_xxxz_xzzz[i] * tbe_0 - 14.0 * tr_xxxz_xzzz[i] * tke_0 + 4.0 * tr_xxxz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxzz_xzz[i] * tbe_0 + 8.0 * tr_xxxzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_yyyy[i] = -4.0 * tr_xxx_yyyyz[i] * tke_0 - 6.0 * tr_xxxz_yyyy[i] * tbe_0 - 2.0 * tr_xxxz_yyyy[i] * tke_0 + 4.0 * tr_xxxz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_yyyz[i] = 2.0 * tr_xxx_yyy[i] - 4.0 * tr_xxx_yyyzz[i] * tke_0 - 6.0 * tr_xxxz_yyyz[i] * tbe_0 - 6.0 * tr_xxxz_yyyz[i] * tke_0 + 4.0 * tr_xxxz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_yyy[i] * tbe_0 + 8.0 * tr_xxxzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_yyzz[i] = 4.0 * tr_xxx_yyz[i] - 4.0 * tr_xxx_yyzzz[i] * tke_0 + 2.0 * tr_xxxz_yy[i] - 6.0 * tr_xxxz_yyzz[i] * tbe_0 - 10.0 * tr_xxxz_yyzz[i] * tke_0 + 4.0 * tr_xxxz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_yyz[i] * tbe_0 + 8.0 * tr_xxxzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_yzzz[i] = 6.0 * tr_xxx_yzz[i] - 4.0 * tr_xxx_yzzzz[i] * tke_0 + 6.0 * tr_xxxz_yz[i] - 6.0 * tr_xxxz_yzzz[i] * tbe_0 - 14.0 * tr_xxxz_yzzz[i] * tke_0 + 4.0 * tr_xxxz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxxzz_yzz[i] * tbe_0 + 8.0 * tr_xxxzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_zzzz[i] = 8.0 * tr_xxx_zzz[i] - 4.0 * tr_xxx_zzzzz[i] * tke_0 + 12.0 * tr_xxxz_zz[i] - 6.0 * tr_xxxz_zzzz[i] * tbe_0 - 18.0 * tr_xxxz_zzzz[i] * tke_0 + 4.0 * tr_xxxz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xxxzz_zzz[i] * tbe_0 + 8.0 * tr_xxxzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1170-1185 components of targeted buffer : GG

    auto tr_0_0_zz_xxyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 1170);

    auto tr_0_0_zz_xxyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 1171);

    auto tr_0_0_zz_xxyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 1172);

    auto tr_0_0_zz_xxyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 1173);

    auto tr_0_0_zz_xxyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 1174);

    auto tr_0_0_zz_xxyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 1175);

    auto tr_0_0_zz_xxyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 1176);

    auto tr_0_0_zz_xxyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 1177);

    auto tr_0_0_zz_xxyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 1178);

    auto tr_0_0_zz_xxyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 1179);

    auto tr_0_0_zz_xxyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 1180);

    auto tr_0_0_zz_xxyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 1181);

    auto tr_0_0_zz_xxyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 1182);

    auto tr_0_0_zz_xxyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 1183);

    auto tr_0_0_zz_xxyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 1184);

    #pragma omp simd aligned(tr_0_0_zz_xxyy_xxxx, tr_0_0_zz_xxyy_xxxy, tr_0_0_zz_xxyy_xxxz, tr_0_0_zz_xxyy_xxyy, tr_0_0_zz_xxyy_xxyz, tr_0_0_zz_xxyy_xxzz, tr_0_0_zz_xxyy_xyyy, tr_0_0_zz_xxyy_xyyz, tr_0_0_zz_xxyy_xyzz, tr_0_0_zz_xxyy_xzzz, tr_0_0_zz_xxyy_yyyy, tr_0_0_zz_xxyy_yyyz, tr_0_0_zz_xxyy_yyzz, tr_0_0_zz_xxyy_yzzz, tr_0_0_zz_xxyy_zzzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxxzz, tr_xxyy_xxxy, tr_xxyy_xxxyzz, tr_xxyy_xxxz, tr_xxyy_xxxzzz, tr_xxyy_xxyy, tr_xxyy_xxyyzz, tr_xxyy_xxyz, tr_xxyy_xxyzzz, tr_xxyy_xxzz, tr_xxyy_xxzzzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyyzz, tr_xxyy_xyyz, tr_xxyy_xyyzzz, tr_xxyy_xyzz, tr_xxyy_xyzzzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_xzzzzz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyyzz, tr_xxyy_yyyz, tr_xxyy_yyyzzz, tr_xxyy_yyzz, tr_xxyy_yyzzzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_yzzzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyy_zzzzzz, tr_xxyyz_xxx, tr_xxyyz_xxxxz, tr_xxyyz_xxxyz, tr_xxyyz_xxxzz, tr_xxyyz_xxy, tr_xxyyz_xxyyz, tr_xxyyz_xxyzz, tr_xxyyz_xxz, tr_xxyyz_xxzzz, tr_xxyyz_xyy, tr_xxyyz_xyyyz, tr_xxyyz_xyyzz, tr_xxyyz_xyz, tr_xxyyz_xyzzz, tr_xxyyz_xzz, tr_xxyyz_xzzzz, tr_xxyyz_yyy, tr_xxyyz_yyyyz, tr_xxyyz_yyyzz, tr_xxyyz_yyz, tr_xxyyz_yyzzz, tr_xxyyz_yzz, tr_xxyyz_yzzzz, tr_xxyyz_zzz, tr_xxyyz_zzzzz, tr_xxyyzz_xxxx, tr_xxyyzz_xxxy, tr_xxyyzz_xxxz, tr_xxyyzz_xxyy, tr_xxyyzz_xxyz, tr_xxyyzz_xxzz, tr_xxyyzz_xyyy, tr_xxyyzz_xyyz, tr_xxyyzz_xyzz, tr_xxyyzz_xzzz, tr_xxyyzz_yyyy, tr_xxyyzz_yyyz, tr_xxyyzz_yyzz, tr_xxyyzz_yzzz, tr_xxyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxyy_xxxx[i] = -2.0 * tr_xxyy_xxxx[i] * tbe_0 - 2.0 * tr_xxyy_xxxx[i] * tke_0 + 4.0 * tr_xxyy_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xxxy[i] = -2.0 * tr_xxyy_xxxy[i] * tbe_0 - 2.0 * tr_xxyy_xxxy[i] * tke_0 + 4.0 * tr_xxyy_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xxxz[i] = -2.0 * tr_xxyy_xxxz[i] * tbe_0 - 6.0 * tr_xxyy_xxxz[i] * tke_0 + 4.0 * tr_xxyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xxx[i] * tbe_0 + 8.0 * tr_xxyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xxyy[i] = -2.0 * tr_xxyy_xxyy[i] * tbe_0 - 2.0 * tr_xxyy_xxyy[i] * tke_0 + 4.0 * tr_xxyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xxyz[i] = -2.0 * tr_xxyy_xxyz[i] * tbe_0 - 6.0 * tr_xxyy_xxyz[i] * tke_0 + 4.0 * tr_xxyy_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xxy[i] * tbe_0 + 8.0 * tr_xxyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xxzz[i] = 2.0 * tr_xxyy_xx[i] - 2.0 * tr_xxyy_xxzz[i] * tbe_0 - 10.0 * tr_xxyy_xxzz[i] * tke_0 + 4.0 * tr_xxyy_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xxz[i] * tbe_0 + 8.0 * tr_xxyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xyyy[i] = -2.0 * tr_xxyy_xyyy[i] * tbe_0 - 2.0 * tr_xxyy_xyyy[i] * tke_0 + 4.0 * tr_xxyy_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xyyz[i] = -2.0 * tr_xxyy_xyyz[i] * tbe_0 - 6.0 * tr_xxyy_xyyz[i] * tke_0 + 4.0 * tr_xxyy_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_xyy[i] * tbe_0 + 8.0 * tr_xxyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xyzz[i] = 2.0 * tr_xxyy_xy[i] - 2.0 * tr_xxyy_xyzz[i] * tbe_0 - 10.0 * tr_xxyy_xyzz[i] * tke_0 + 4.0 * tr_xxyy_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_xyz[i] * tbe_0 + 8.0 * tr_xxyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xzzz[i] = 6.0 * tr_xxyy_xz[i] - 2.0 * tr_xxyy_xzzz[i] * tbe_0 - 14.0 * tr_xxyy_xzzz[i] * tke_0 + 4.0 * tr_xxyy_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxyyz_xzz[i] * tbe_0 + 8.0 * tr_xxyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_yyyy[i] = -2.0 * tr_xxyy_yyyy[i] * tbe_0 - 2.0 * tr_xxyy_yyyy[i] * tke_0 + 4.0 * tr_xxyy_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_yyyz[i] = -2.0 * tr_xxyy_yyyz[i] * tbe_0 - 6.0 * tr_xxyy_yyyz[i] * tke_0 + 4.0 * tr_xxyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_yyy[i] * tbe_0 + 8.0 * tr_xxyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_yyzz[i] = 2.0 * tr_xxyy_yy[i] - 2.0 * tr_xxyy_yyzz[i] * tbe_0 - 10.0 * tr_xxyy_yyzz[i] * tke_0 + 4.0 * tr_xxyy_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_yyz[i] * tbe_0 + 8.0 * tr_xxyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_yzzz[i] = 6.0 * tr_xxyy_yz[i] - 2.0 * tr_xxyy_yzzz[i] * tbe_0 - 14.0 * tr_xxyy_yzzz[i] * tke_0 + 4.0 * tr_xxyy_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxyyz_yzz[i] * tbe_0 + 8.0 * tr_xxyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_zzzz[i] = 12.0 * tr_xxyy_zz[i] - 2.0 * tr_xxyy_zzzz[i] * tbe_0 - 18.0 * tr_xxyy_zzzz[i] * tke_0 + 4.0 * tr_xxyy_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xxyyz_zzz[i] * tbe_0 + 8.0 * tr_xxyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1185-1200 components of targeted buffer : GG

    auto tr_0_0_zz_xxyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1185);

    auto tr_0_0_zz_xxyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1186);

    auto tr_0_0_zz_xxyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1187);

    auto tr_0_0_zz_xxyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1188);

    auto tr_0_0_zz_xxyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1189);

    auto tr_0_0_zz_xxyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1190);

    auto tr_0_0_zz_xxyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1191);

    auto tr_0_0_zz_xxyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1192);

    auto tr_0_0_zz_xxyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1193);

    auto tr_0_0_zz_xxyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1194);

    auto tr_0_0_zz_xxyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1195);

    auto tr_0_0_zz_xxyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1196);

    auto tr_0_0_zz_xxyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1197);

    auto tr_0_0_zz_xxyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1198);

    auto tr_0_0_zz_xxyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1199);

    #pragma omp simd aligned(tr_0_0_zz_xxyz_xxxx, tr_0_0_zz_xxyz_xxxy, tr_0_0_zz_xxyz_xxxz, tr_0_0_zz_xxyz_xxyy, tr_0_0_zz_xxyz_xxyz, tr_0_0_zz_xxyz_xxzz, tr_0_0_zz_xxyz_xyyy, tr_0_0_zz_xxyz_xyyz, tr_0_0_zz_xxyz_xyzz, tr_0_0_zz_xxyz_xzzz, tr_0_0_zz_xxyz_yyyy, tr_0_0_zz_xxyz_yyyz, tr_0_0_zz_xxyz_yyzz, tr_0_0_zz_xxyz_yzzz, tr_0_0_zz_xxyz_zzzz, tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxxzz, tr_xxyz_xxxy, tr_xxyz_xxxyzz, tr_xxyz_xxxz, tr_xxyz_xxxzzz, tr_xxyz_xxyy, tr_xxyz_xxyyzz, tr_xxyz_xxyz, tr_xxyz_xxyzzz, tr_xxyz_xxzz, tr_xxyz_xxzzzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyyzz, tr_xxyz_xyyz, tr_xxyz_xyyzzz, tr_xxyz_xyzz, tr_xxyz_xyzzzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_xzzzzz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyyzz, tr_xxyz_yyyz, tr_xxyz_yyyzzz, tr_xxyz_yyzz, tr_xxyz_yyzzzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_yzzzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyz_zzzzzz, tr_xxyzz_xxx, tr_xxyzz_xxxxz, tr_xxyzz_xxxyz, tr_xxyzz_xxxzz, tr_xxyzz_xxy, tr_xxyzz_xxyyz, tr_xxyzz_xxyzz, tr_xxyzz_xxz, tr_xxyzz_xxzzz, tr_xxyzz_xyy, tr_xxyzz_xyyyz, tr_xxyzz_xyyzz, tr_xxyzz_xyz, tr_xxyzz_xyzzz, tr_xxyzz_xzz, tr_xxyzz_xzzzz, tr_xxyzz_yyy, tr_xxyzz_yyyyz, tr_xxyzz_yyyzz, tr_xxyzz_yyz, tr_xxyzz_yyzzz, tr_xxyzz_yzz, tr_xxyzz_yzzzz, tr_xxyzz_zzz, tr_xxyzz_zzzzz, tr_xxyzzz_xxxx, tr_xxyzzz_xxxy, tr_xxyzzz_xxxz, tr_xxyzzz_xxyy, tr_xxyzzz_xxyz, tr_xxyzzz_xxzz, tr_xxyzzz_xyyy, tr_xxyzzz_xyyz, tr_xxyzzz_xyzz, tr_xxyzzz_xzzz, tr_xxyzzz_yyyy, tr_xxyzzz_yyyz, tr_xxyzzz_yyzz, tr_xxyzzz_yzzz, tr_xxyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxyz_xxxx[i] = -4.0 * tr_xxy_xxxxz[i] * tke_0 - 6.0 * tr_xxyz_xxxx[i] * tbe_0 - 2.0 * tr_xxyz_xxxx[i] * tke_0 + 4.0 * tr_xxyz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xxxy[i] = -4.0 * tr_xxy_xxxyz[i] * tke_0 - 6.0 * tr_xxyz_xxxy[i] * tbe_0 - 2.0 * tr_xxyz_xxxy[i] * tke_0 + 4.0 * tr_xxyz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xxxz[i] = 2.0 * tr_xxy_xxx[i] - 4.0 * tr_xxy_xxxzz[i] * tke_0 - 6.0 * tr_xxyz_xxxz[i] * tbe_0 - 6.0 * tr_xxyz_xxxz[i] * tke_0 + 4.0 * tr_xxyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xxx[i] * tbe_0 + 8.0 * tr_xxyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xxyy[i] = -4.0 * tr_xxy_xxyyz[i] * tke_0 - 6.0 * tr_xxyz_xxyy[i] * tbe_0 - 2.0 * tr_xxyz_xxyy[i] * tke_0 + 4.0 * tr_xxyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xxyz[i] = 2.0 * tr_xxy_xxy[i] - 4.0 * tr_xxy_xxyzz[i] * tke_0 - 6.0 * tr_xxyz_xxyz[i] * tbe_0 - 6.0 * tr_xxyz_xxyz[i] * tke_0 + 4.0 * tr_xxyz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xxy[i] * tbe_0 + 8.0 * tr_xxyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xxzz[i] = 4.0 * tr_xxy_xxz[i] - 4.0 * tr_xxy_xxzzz[i] * tke_0 + 2.0 * tr_xxyz_xx[i] - 6.0 * tr_xxyz_xxzz[i] * tbe_0 - 10.0 * tr_xxyz_xxzz[i] * tke_0 + 4.0 * tr_xxyz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xxz[i] * tbe_0 + 8.0 * tr_xxyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xyyy[i] = -4.0 * tr_xxy_xyyyz[i] * tke_0 - 6.0 * tr_xxyz_xyyy[i] * tbe_0 - 2.0 * tr_xxyz_xyyy[i] * tke_0 + 4.0 * tr_xxyz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xyyz[i] = 2.0 * tr_xxy_xyy[i] - 4.0 * tr_xxy_xyyzz[i] * tke_0 - 6.0 * tr_xxyz_xyyz[i] * tbe_0 - 6.0 * tr_xxyz_xyyz[i] * tke_0 + 4.0 * tr_xxyz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_xyy[i] * tbe_0 + 8.0 * tr_xxyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xyzz[i] = 4.0 * tr_xxy_xyz[i] - 4.0 * tr_xxy_xyzzz[i] * tke_0 + 2.0 * tr_xxyz_xy[i] - 6.0 * tr_xxyz_xyzz[i] * tbe_0 - 10.0 * tr_xxyz_xyzz[i] * tke_0 + 4.0 * tr_xxyz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_xyz[i] * tbe_0 + 8.0 * tr_xxyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xzzz[i] = 6.0 * tr_xxy_xzz[i] - 4.0 * tr_xxy_xzzzz[i] * tke_0 + 6.0 * tr_xxyz_xz[i] - 6.0 * tr_xxyz_xzzz[i] * tbe_0 - 14.0 * tr_xxyz_xzzz[i] * tke_0 + 4.0 * tr_xxyz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxyzz_xzz[i] * tbe_0 + 8.0 * tr_xxyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_yyyy[i] = -4.0 * tr_xxy_yyyyz[i] * tke_0 - 6.0 * tr_xxyz_yyyy[i] * tbe_0 - 2.0 * tr_xxyz_yyyy[i] * tke_0 + 4.0 * tr_xxyz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_yyyz[i] = 2.0 * tr_xxy_yyy[i] - 4.0 * tr_xxy_yyyzz[i] * tke_0 - 6.0 * tr_xxyz_yyyz[i] * tbe_0 - 6.0 * tr_xxyz_yyyz[i] * tke_0 + 4.0 * tr_xxyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_yyy[i] * tbe_0 + 8.0 * tr_xxyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_yyzz[i] = 4.0 * tr_xxy_yyz[i] - 4.0 * tr_xxy_yyzzz[i] * tke_0 + 2.0 * tr_xxyz_yy[i] - 6.0 * tr_xxyz_yyzz[i] * tbe_0 - 10.0 * tr_xxyz_yyzz[i] * tke_0 + 4.0 * tr_xxyz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_yyz[i] * tbe_0 + 8.0 * tr_xxyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_yzzz[i] = 6.0 * tr_xxy_yzz[i] - 4.0 * tr_xxy_yzzzz[i] * tke_0 + 6.0 * tr_xxyz_yz[i] - 6.0 * tr_xxyz_yzzz[i] * tbe_0 - 14.0 * tr_xxyz_yzzz[i] * tke_0 + 4.0 * tr_xxyz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxyzz_yzz[i] * tbe_0 + 8.0 * tr_xxyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_zzzz[i] = 8.0 * tr_xxy_zzz[i] - 4.0 * tr_xxy_zzzzz[i] * tke_0 + 12.0 * tr_xxyz_zz[i] - 6.0 * tr_xxyz_zzzz[i] * tbe_0 - 18.0 * tr_xxyz_zzzz[i] * tke_0 + 4.0 * tr_xxyz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xxyzz_zzz[i] * tbe_0 + 8.0 * tr_xxyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1200-1215 components of targeted buffer : GG

    auto tr_0_0_zz_xxzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1200);

    auto tr_0_0_zz_xxzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1201);

    auto tr_0_0_zz_xxzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1202);

    auto tr_0_0_zz_xxzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1203);

    auto tr_0_0_zz_xxzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1204);

    auto tr_0_0_zz_xxzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1205);

    auto tr_0_0_zz_xxzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1206);

    auto tr_0_0_zz_xxzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1207);

    auto tr_0_0_zz_xxzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1208);

    auto tr_0_0_zz_xxzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1209);

    auto tr_0_0_zz_xxzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1210);

    auto tr_0_0_zz_xxzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1211);

    auto tr_0_0_zz_xxzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1212);

    auto tr_0_0_zz_xxzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1213);

    auto tr_0_0_zz_xxzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1214);

    #pragma omp simd aligned(tr_0_0_zz_xxzz_xxxx, tr_0_0_zz_xxzz_xxxy, tr_0_0_zz_xxzz_xxxz, tr_0_0_zz_xxzz_xxyy, tr_0_0_zz_xxzz_xxyz, tr_0_0_zz_xxzz_xxzz, tr_0_0_zz_xxzz_xyyy, tr_0_0_zz_xxzz_xyyz, tr_0_0_zz_xxzz_xyzz, tr_0_0_zz_xxzz_xzzz, tr_0_0_zz_xxzz_yyyy, tr_0_0_zz_xxzz_yyyz, tr_0_0_zz_xxzz_yyzz, tr_0_0_zz_xxzz_yzzz, tr_0_0_zz_xxzz_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxxzz, tr_xxzz_xxxy, tr_xxzz_xxxyzz, tr_xxzz_xxxz, tr_xxzz_xxxzzz, tr_xxzz_xxyy, tr_xxzz_xxyyzz, tr_xxzz_xxyz, tr_xxzz_xxyzzz, tr_xxzz_xxzz, tr_xxzz_xxzzzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyyzz, tr_xxzz_xyyz, tr_xxzz_xyyzzz, tr_xxzz_xyzz, tr_xxzz_xyzzzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_xzzzzz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyyzz, tr_xxzz_yyyz, tr_xxzz_yyyzzz, tr_xxzz_yyzz, tr_xxzz_yyzzzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_yzzzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzz_zzzzzz, tr_xxzzz_xxx, tr_xxzzz_xxxxz, tr_xxzzz_xxxyz, tr_xxzzz_xxxzz, tr_xxzzz_xxy, tr_xxzzz_xxyyz, tr_xxzzz_xxyzz, tr_xxzzz_xxz, tr_xxzzz_xxzzz, tr_xxzzz_xyy, tr_xxzzz_xyyyz, tr_xxzzz_xyyzz, tr_xxzzz_xyz, tr_xxzzz_xyzzz, tr_xxzzz_xzz, tr_xxzzz_xzzzz, tr_xxzzz_yyy, tr_xxzzz_yyyyz, tr_xxzzz_yyyzz, tr_xxzzz_yyz, tr_xxzzz_yyzzz, tr_xxzzz_yzz, tr_xxzzz_yzzzz, tr_xxzzz_zzz, tr_xxzzz_zzzzz, tr_xxzzzz_xxxx, tr_xxzzzz_xxxy, tr_xxzzzz_xxxz, tr_xxzzzz_xxyy, tr_xxzzzz_xxyz, tr_xxzzzz_xxzz, tr_xxzzzz_xyyy, tr_xxzzzz_xyyz, tr_xxzzzz_xyzz, tr_xxzzzz_xzzz, tr_xxzzzz_yyyy, tr_xxzzzz_yyyz, tr_xxzzzz_yyzz, tr_xxzzzz_yzzz, tr_xxzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxzz_xxxx[i] = 2.0 * tr_xx_xxxx[i] - 8.0 * tr_xxz_xxxxz[i] * tke_0 - 10.0 * tr_xxzz_xxxx[i] * tbe_0 - 2.0 * tr_xxzz_xxxx[i] * tke_0 + 4.0 * tr_xxzz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xxxy[i] = 2.0 * tr_xx_xxxy[i] - 8.0 * tr_xxz_xxxyz[i] * tke_0 - 10.0 * tr_xxzz_xxxy[i] * tbe_0 - 2.0 * tr_xxzz_xxxy[i] * tke_0 + 4.0 * tr_xxzz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xxxz[i] = 2.0 * tr_xx_xxxz[i] + 4.0 * tr_xxz_xxx[i] - 8.0 * tr_xxz_xxxzz[i] * tke_0 - 10.0 * tr_xxzz_xxxz[i] * tbe_0 - 6.0 * tr_xxzz_xxxz[i] * tke_0 + 4.0 * tr_xxzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xxx[i] * tbe_0 + 8.0 * tr_xxzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xxyy[i] = 2.0 * tr_xx_xxyy[i] - 8.0 * tr_xxz_xxyyz[i] * tke_0 - 10.0 * tr_xxzz_xxyy[i] * tbe_0 - 2.0 * tr_xxzz_xxyy[i] * tke_0 + 4.0 * tr_xxzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xxyz[i] = 2.0 * tr_xx_xxyz[i] + 4.0 * tr_xxz_xxy[i] - 8.0 * tr_xxz_xxyzz[i] * tke_0 - 10.0 * tr_xxzz_xxyz[i] * tbe_0 - 6.0 * tr_xxzz_xxyz[i] * tke_0 + 4.0 * tr_xxzz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xxy[i] * tbe_0 + 8.0 * tr_xxzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xxzz[i] = 2.0 * tr_xx_xxzz[i] + 8.0 * tr_xxz_xxz[i] - 8.0 * tr_xxz_xxzzz[i] * tke_0 + 2.0 * tr_xxzz_xx[i] - 10.0 * tr_xxzz_xxzz[i] * tbe_0 - 10.0 * tr_xxzz_xxzz[i] * tke_0 + 4.0 * tr_xxzz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_xxz[i] * tbe_0 + 8.0 * tr_xxzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xyyy[i] = 2.0 * tr_xx_xyyy[i] - 8.0 * tr_xxz_xyyyz[i] * tke_0 - 10.0 * tr_xxzz_xyyy[i] * tbe_0 - 2.0 * tr_xxzz_xyyy[i] * tke_0 + 4.0 * tr_xxzz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xyyz[i] = 2.0 * tr_xx_xyyz[i] + 4.0 * tr_xxz_xyy[i] - 8.0 * tr_xxz_xyyzz[i] * tke_0 - 10.0 * tr_xxzz_xyyz[i] * tbe_0 - 6.0 * tr_xxzz_xyyz[i] * tke_0 + 4.0 * tr_xxzz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_xyy[i] * tbe_0 + 8.0 * tr_xxzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xyzz[i] = 2.0 * tr_xx_xyzz[i] + 8.0 * tr_xxz_xyz[i] - 8.0 * tr_xxz_xyzzz[i] * tke_0 + 2.0 * tr_xxzz_xy[i] - 10.0 * tr_xxzz_xyzz[i] * tbe_0 - 10.0 * tr_xxzz_xyzz[i] * tke_0 + 4.0 * tr_xxzz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_xyz[i] * tbe_0 + 8.0 * tr_xxzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xzzz[i] = 2.0 * tr_xx_xzzz[i] + 12.0 * tr_xxz_xzz[i] - 8.0 * tr_xxz_xzzzz[i] * tke_0 + 6.0 * tr_xxzz_xz[i] - 10.0 * tr_xxzz_xzzz[i] * tbe_0 - 14.0 * tr_xxzz_xzzz[i] * tke_0 + 4.0 * tr_xxzz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxzzz_xzz[i] * tbe_0 + 8.0 * tr_xxzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_yyyy[i] = 2.0 * tr_xx_yyyy[i] - 8.0 * tr_xxz_yyyyz[i] * tke_0 - 10.0 * tr_xxzz_yyyy[i] * tbe_0 - 2.0 * tr_xxzz_yyyy[i] * tke_0 + 4.0 * tr_xxzz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_yyyz[i] = 2.0 * tr_xx_yyyz[i] + 4.0 * tr_xxz_yyy[i] - 8.0 * tr_xxz_yyyzz[i] * tke_0 - 10.0 * tr_xxzz_yyyz[i] * tbe_0 - 6.0 * tr_xxzz_yyyz[i] * tke_0 + 4.0 * tr_xxzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_yyy[i] * tbe_0 + 8.0 * tr_xxzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_yyzz[i] = 2.0 * tr_xx_yyzz[i] + 8.0 * tr_xxz_yyz[i] - 8.0 * tr_xxz_yyzzz[i] * tke_0 + 2.0 * tr_xxzz_yy[i] - 10.0 * tr_xxzz_yyzz[i] * tbe_0 - 10.0 * tr_xxzz_yyzz[i] * tke_0 + 4.0 * tr_xxzz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_yyz[i] * tbe_0 + 8.0 * tr_xxzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_yzzz[i] = 2.0 * tr_xx_yzzz[i] + 12.0 * tr_xxz_yzz[i] - 8.0 * tr_xxz_yzzzz[i] * tke_0 + 6.0 * tr_xxzz_yz[i] - 10.0 * tr_xxzz_yzzz[i] * tbe_0 - 14.0 * tr_xxzz_yzzz[i] * tke_0 + 4.0 * tr_xxzz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xxzzz_yzz[i] * tbe_0 + 8.0 * tr_xxzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_zzzz[i] = 2.0 * tr_xx_zzzz[i] + 16.0 * tr_xxz_zzz[i] - 8.0 * tr_xxz_zzzzz[i] * tke_0 + 12.0 * tr_xxzz_zz[i] - 10.0 * tr_xxzz_zzzz[i] * tbe_0 - 18.0 * tr_xxzz_zzzz[i] * tke_0 + 4.0 * tr_xxzz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xxzzz_zzz[i] * tbe_0 + 8.0 * tr_xxzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1215-1230 components of targeted buffer : GG

    auto tr_0_0_zz_xyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 1215);

    auto tr_0_0_zz_xyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 1216);

    auto tr_0_0_zz_xyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 1217);

    auto tr_0_0_zz_xyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 1218);

    auto tr_0_0_zz_xyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 1219);

    auto tr_0_0_zz_xyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 1220);

    auto tr_0_0_zz_xyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 1221);

    auto tr_0_0_zz_xyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 1222);

    auto tr_0_0_zz_xyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 1223);

    auto tr_0_0_zz_xyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 1224);

    auto tr_0_0_zz_xyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 1225);

    auto tr_0_0_zz_xyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 1226);

    auto tr_0_0_zz_xyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 1227);

    auto tr_0_0_zz_xyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 1228);

    auto tr_0_0_zz_xyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 1229);

    #pragma omp simd aligned(tr_0_0_zz_xyyy_xxxx, tr_0_0_zz_xyyy_xxxy, tr_0_0_zz_xyyy_xxxz, tr_0_0_zz_xyyy_xxyy, tr_0_0_zz_xyyy_xxyz, tr_0_0_zz_xyyy_xxzz, tr_0_0_zz_xyyy_xyyy, tr_0_0_zz_xyyy_xyyz, tr_0_0_zz_xyyy_xyzz, tr_0_0_zz_xyyy_xzzz, tr_0_0_zz_xyyy_yyyy, tr_0_0_zz_xyyy_yyyz, tr_0_0_zz_xyyy_yyzz, tr_0_0_zz_xyyy_yzzz, tr_0_0_zz_xyyy_zzzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxxzz, tr_xyyy_xxxy, tr_xyyy_xxxyzz, tr_xyyy_xxxz, tr_xyyy_xxxzzz, tr_xyyy_xxyy, tr_xyyy_xxyyzz, tr_xyyy_xxyz, tr_xyyy_xxyzzz, tr_xyyy_xxzz, tr_xyyy_xxzzzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyyzz, tr_xyyy_xyyz, tr_xyyy_xyyzzz, tr_xyyy_xyzz, tr_xyyy_xyzzzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_xzzzzz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyyzz, tr_xyyy_yyyz, tr_xyyy_yyyzzz, tr_xyyy_yyzz, tr_xyyy_yyzzzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_yzzzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyy_zzzzzz, tr_xyyyz_xxx, tr_xyyyz_xxxxz, tr_xyyyz_xxxyz, tr_xyyyz_xxxzz, tr_xyyyz_xxy, tr_xyyyz_xxyyz, tr_xyyyz_xxyzz, tr_xyyyz_xxz, tr_xyyyz_xxzzz, tr_xyyyz_xyy, tr_xyyyz_xyyyz, tr_xyyyz_xyyzz, tr_xyyyz_xyz, tr_xyyyz_xyzzz, tr_xyyyz_xzz, tr_xyyyz_xzzzz, tr_xyyyz_yyy, tr_xyyyz_yyyyz, tr_xyyyz_yyyzz, tr_xyyyz_yyz, tr_xyyyz_yyzzz, tr_xyyyz_yzz, tr_xyyyz_yzzzz, tr_xyyyz_zzz, tr_xyyyz_zzzzz, tr_xyyyzz_xxxx, tr_xyyyzz_xxxy, tr_xyyyzz_xxxz, tr_xyyyzz_xxyy, tr_xyyyzz_xxyz, tr_xyyyzz_xxzz, tr_xyyyzz_xyyy, tr_xyyyzz_xyyz, tr_xyyyzz_xyzz, tr_xyyyzz_xzzz, tr_xyyyzz_yyyy, tr_xyyyzz_yyyz, tr_xyyyzz_yyzz, tr_xyyyzz_yzzz, tr_xyyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyyy_xxxx[i] = -2.0 * tr_xyyy_xxxx[i] * tbe_0 - 2.0 * tr_xyyy_xxxx[i] * tke_0 + 4.0 * tr_xyyy_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xxxy[i] = -2.0 * tr_xyyy_xxxy[i] * tbe_0 - 2.0 * tr_xyyy_xxxy[i] * tke_0 + 4.0 * tr_xyyy_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xxxz[i] = -2.0 * tr_xyyy_xxxz[i] * tbe_0 - 6.0 * tr_xyyy_xxxz[i] * tke_0 + 4.0 * tr_xyyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xxx[i] * tbe_0 + 8.0 * tr_xyyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xxyy[i] = -2.0 * tr_xyyy_xxyy[i] * tbe_0 - 2.0 * tr_xyyy_xxyy[i] * tke_0 + 4.0 * tr_xyyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xxyz[i] = -2.0 * tr_xyyy_xxyz[i] * tbe_0 - 6.0 * tr_xyyy_xxyz[i] * tke_0 + 4.0 * tr_xyyy_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xxy[i] * tbe_0 + 8.0 * tr_xyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xxzz[i] = 2.0 * tr_xyyy_xx[i] - 2.0 * tr_xyyy_xxzz[i] * tbe_0 - 10.0 * tr_xyyy_xxzz[i] * tke_0 + 4.0 * tr_xyyy_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xxz[i] * tbe_0 + 8.0 * tr_xyyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xyyy[i] = -2.0 * tr_xyyy_xyyy[i] * tbe_0 - 2.0 * tr_xyyy_xyyy[i] * tke_0 + 4.0 * tr_xyyy_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xyyz[i] = -2.0 * tr_xyyy_xyyz[i] * tbe_0 - 6.0 * tr_xyyy_xyyz[i] * tke_0 + 4.0 * tr_xyyy_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_xyy[i] * tbe_0 + 8.0 * tr_xyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xyzz[i] = 2.0 * tr_xyyy_xy[i] - 2.0 * tr_xyyy_xyzz[i] * tbe_0 - 10.0 * tr_xyyy_xyzz[i] * tke_0 + 4.0 * tr_xyyy_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_xyz[i] * tbe_0 + 8.0 * tr_xyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xzzz[i] = 6.0 * tr_xyyy_xz[i] - 2.0 * tr_xyyy_xzzz[i] * tbe_0 - 14.0 * tr_xyyy_xzzz[i] * tke_0 + 4.0 * tr_xyyy_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyyyz_xzz[i] * tbe_0 + 8.0 * tr_xyyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_yyyy[i] = -2.0 * tr_xyyy_yyyy[i] * tbe_0 - 2.0 * tr_xyyy_yyyy[i] * tke_0 + 4.0 * tr_xyyy_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_yyyz[i] = -2.0 * tr_xyyy_yyyz[i] * tbe_0 - 6.0 * tr_xyyy_yyyz[i] * tke_0 + 4.0 * tr_xyyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_yyy[i] * tbe_0 + 8.0 * tr_xyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_yyzz[i] = 2.0 * tr_xyyy_yy[i] - 2.0 * tr_xyyy_yyzz[i] * tbe_0 - 10.0 * tr_xyyy_yyzz[i] * tke_0 + 4.0 * tr_xyyy_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_yyz[i] * tbe_0 + 8.0 * tr_xyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_yzzz[i] = 6.0 * tr_xyyy_yz[i] - 2.0 * tr_xyyy_yzzz[i] * tbe_0 - 14.0 * tr_xyyy_yzzz[i] * tke_0 + 4.0 * tr_xyyy_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyyyz_yzz[i] * tbe_0 + 8.0 * tr_xyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_zzzz[i] = 12.0 * tr_xyyy_zz[i] - 2.0 * tr_xyyy_zzzz[i] * tbe_0 - 18.0 * tr_xyyy_zzzz[i] * tke_0 + 4.0 * tr_xyyy_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xyyyz_zzz[i] * tbe_0 + 8.0 * tr_xyyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1230-1245 components of targeted buffer : GG

    auto tr_0_0_zz_xyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1230);

    auto tr_0_0_zz_xyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1231);

    auto tr_0_0_zz_xyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1232);

    auto tr_0_0_zz_xyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1233);

    auto tr_0_0_zz_xyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1234);

    auto tr_0_0_zz_xyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1235);

    auto tr_0_0_zz_xyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1236);

    auto tr_0_0_zz_xyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1237);

    auto tr_0_0_zz_xyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1238);

    auto tr_0_0_zz_xyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1239);

    auto tr_0_0_zz_xyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1240);

    auto tr_0_0_zz_xyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1241);

    auto tr_0_0_zz_xyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1242);

    auto tr_0_0_zz_xyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1243);

    auto tr_0_0_zz_xyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1244);

    #pragma omp simd aligned(tr_0_0_zz_xyyz_xxxx, tr_0_0_zz_xyyz_xxxy, tr_0_0_zz_xyyz_xxxz, tr_0_0_zz_xyyz_xxyy, tr_0_0_zz_xyyz_xxyz, tr_0_0_zz_xyyz_xxzz, tr_0_0_zz_xyyz_xyyy, tr_0_0_zz_xyyz_xyyz, tr_0_0_zz_xyyz_xyzz, tr_0_0_zz_xyyz_xzzz, tr_0_0_zz_xyyz_yyyy, tr_0_0_zz_xyyz_yyyz, tr_0_0_zz_xyyz_yyzz, tr_0_0_zz_xyyz_yzzz, tr_0_0_zz_xyyz_zzzz, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxxzz, tr_xyyz_xxxy, tr_xyyz_xxxyzz, tr_xyyz_xxxz, tr_xyyz_xxxzzz, tr_xyyz_xxyy, tr_xyyz_xxyyzz, tr_xyyz_xxyz, tr_xyyz_xxyzzz, tr_xyyz_xxzz, tr_xyyz_xxzzzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyyzz, tr_xyyz_xyyz, tr_xyyz_xyyzzz, tr_xyyz_xyzz, tr_xyyz_xyzzzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_xzzzzz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyyzz, tr_xyyz_yyyz, tr_xyyz_yyyzzz, tr_xyyz_yyzz, tr_xyyz_yyzzzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_yzzzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyz_zzzzzz, tr_xyyzz_xxx, tr_xyyzz_xxxxz, tr_xyyzz_xxxyz, tr_xyyzz_xxxzz, tr_xyyzz_xxy, tr_xyyzz_xxyyz, tr_xyyzz_xxyzz, tr_xyyzz_xxz, tr_xyyzz_xxzzz, tr_xyyzz_xyy, tr_xyyzz_xyyyz, tr_xyyzz_xyyzz, tr_xyyzz_xyz, tr_xyyzz_xyzzz, tr_xyyzz_xzz, tr_xyyzz_xzzzz, tr_xyyzz_yyy, tr_xyyzz_yyyyz, tr_xyyzz_yyyzz, tr_xyyzz_yyz, tr_xyyzz_yyzzz, tr_xyyzz_yzz, tr_xyyzz_yzzzz, tr_xyyzz_zzz, tr_xyyzz_zzzzz, tr_xyyzzz_xxxx, tr_xyyzzz_xxxy, tr_xyyzzz_xxxz, tr_xyyzzz_xxyy, tr_xyyzzz_xxyz, tr_xyyzzz_xxzz, tr_xyyzzz_xyyy, tr_xyyzzz_xyyz, tr_xyyzzz_xyzz, tr_xyyzzz_xzzz, tr_xyyzzz_yyyy, tr_xyyzzz_yyyz, tr_xyyzzz_yyzz, tr_xyyzzz_yzzz, tr_xyyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyyz_xxxx[i] = -4.0 * tr_xyy_xxxxz[i] * tke_0 - 6.0 * tr_xyyz_xxxx[i] * tbe_0 - 2.0 * tr_xyyz_xxxx[i] * tke_0 + 4.0 * tr_xyyz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xxxy[i] = -4.0 * tr_xyy_xxxyz[i] * tke_0 - 6.0 * tr_xyyz_xxxy[i] * tbe_0 - 2.0 * tr_xyyz_xxxy[i] * tke_0 + 4.0 * tr_xyyz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xxxz[i] = 2.0 * tr_xyy_xxx[i] - 4.0 * tr_xyy_xxxzz[i] * tke_0 - 6.0 * tr_xyyz_xxxz[i] * tbe_0 - 6.0 * tr_xyyz_xxxz[i] * tke_0 + 4.0 * tr_xyyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xxx[i] * tbe_0 + 8.0 * tr_xyyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xxyy[i] = -4.0 * tr_xyy_xxyyz[i] * tke_0 - 6.0 * tr_xyyz_xxyy[i] * tbe_0 - 2.0 * tr_xyyz_xxyy[i] * tke_0 + 4.0 * tr_xyyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xxyz[i] = 2.0 * tr_xyy_xxy[i] - 4.0 * tr_xyy_xxyzz[i] * tke_0 - 6.0 * tr_xyyz_xxyz[i] * tbe_0 - 6.0 * tr_xyyz_xxyz[i] * tke_0 + 4.0 * tr_xyyz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xxy[i] * tbe_0 + 8.0 * tr_xyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xxzz[i] = 4.0 * tr_xyy_xxz[i] - 4.0 * tr_xyy_xxzzz[i] * tke_0 + 2.0 * tr_xyyz_xx[i] - 6.0 * tr_xyyz_xxzz[i] * tbe_0 - 10.0 * tr_xyyz_xxzz[i] * tke_0 + 4.0 * tr_xyyz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xxz[i] * tbe_0 + 8.0 * tr_xyyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xyyy[i] = -4.0 * tr_xyy_xyyyz[i] * tke_0 - 6.0 * tr_xyyz_xyyy[i] * tbe_0 - 2.0 * tr_xyyz_xyyy[i] * tke_0 + 4.0 * tr_xyyz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xyyz[i] = 2.0 * tr_xyy_xyy[i] - 4.0 * tr_xyy_xyyzz[i] * tke_0 - 6.0 * tr_xyyz_xyyz[i] * tbe_0 - 6.0 * tr_xyyz_xyyz[i] * tke_0 + 4.0 * tr_xyyz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_xyy[i] * tbe_0 + 8.0 * tr_xyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xyzz[i] = 4.0 * tr_xyy_xyz[i] - 4.0 * tr_xyy_xyzzz[i] * tke_0 + 2.0 * tr_xyyz_xy[i] - 6.0 * tr_xyyz_xyzz[i] * tbe_0 - 10.0 * tr_xyyz_xyzz[i] * tke_0 + 4.0 * tr_xyyz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_xyz[i] * tbe_0 + 8.0 * tr_xyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xzzz[i] = 6.0 * tr_xyy_xzz[i] - 4.0 * tr_xyy_xzzzz[i] * tke_0 + 6.0 * tr_xyyz_xz[i] - 6.0 * tr_xyyz_xzzz[i] * tbe_0 - 14.0 * tr_xyyz_xzzz[i] * tke_0 + 4.0 * tr_xyyz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyyzz_xzz[i] * tbe_0 + 8.0 * tr_xyyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_yyyy[i] = -4.0 * tr_xyy_yyyyz[i] * tke_0 - 6.0 * tr_xyyz_yyyy[i] * tbe_0 - 2.0 * tr_xyyz_yyyy[i] * tke_0 + 4.0 * tr_xyyz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_yyyz[i] = 2.0 * tr_xyy_yyy[i] - 4.0 * tr_xyy_yyyzz[i] * tke_0 - 6.0 * tr_xyyz_yyyz[i] * tbe_0 - 6.0 * tr_xyyz_yyyz[i] * tke_0 + 4.0 * tr_xyyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_yyy[i] * tbe_0 + 8.0 * tr_xyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_yyzz[i] = 4.0 * tr_xyy_yyz[i] - 4.0 * tr_xyy_yyzzz[i] * tke_0 + 2.0 * tr_xyyz_yy[i] - 6.0 * tr_xyyz_yyzz[i] * tbe_0 - 10.0 * tr_xyyz_yyzz[i] * tke_0 + 4.0 * tr_xyyz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_yyz[i] * tbe_0 + 8.0 * tr_xyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_yzzz[i] = 6.0 * tr_xyy_yzz[i] - 4.0 * tr_xyy_yzzzz[i] * tke_0 + 6.0 * tr_xyyz_yz[i] - 6.0 * tr_xyyz_yzzz[i] * tbe_0 - 14.0 * tr_xyyz_yzzz[i] * tke_0 + 4.0 * tr_xyyz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyyzz_yzz[i] * tbe_0 + 8.0 * tr_xyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_zzzz[i] = 8.0 * tr_xyy_zzz[i] - 4.0 * tr_xyy_zzzzz[i] * tke_0 + 12.0 * tr_xyyz_zz[i] - 6.0 * tr_xyyz_zzzz[i] * tbe_0 - 18.0 * tr_xyyz_zzzz[i] * tke_0 + 4.0 * tr_xyyz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xyyzz_zzz[i] * tbe_0 + 8.0 * tr_xyyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1245-1260 components of targeted buffer : GG

    auto tr_0_0_zz_xyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1245);

    auto tr_0_0_zz_xyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1246);

    auto tr_0_0_zz_xyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1247);

    auto tr_0_0_zz_xyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1248);

    auto tr_0_0_zz_xyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1249);

    auto tr_0_0_zz_xyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1250);

    auto tr_0_0_zz_xyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1251);

    auto tr_0_0_zz_xyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1252);

    auto tr_0_0_zz_xyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1253);

    auto tr_0_0_zz_xyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1254);

    auto tr_0_0_zz_xyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1255);

    auto tr_0_0_zz_xyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1256);

    auto tr_0_0_zz_xyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1257);

    auto tr_0_0_zz_xyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1258);

    auto tr_0_0_zz_xyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1259);

    #pragma omp simd aligned(tr_0_0_zz_xyzz_xxxx, tr_0_0_zz_xyzz_xxxy, tr_0_0_zz_xyzz_xxxz, tr_0_0_zz_xyzz_xxyy, tr_0_0_zz_xyzz_xxyz, tr_0_0_zz_xyzz_xxzz, tr_0_0_zz_xyzz_xyyy, tr_0_0_zz_xyzz_xyyz, tr_0_0_zz_xyzz_xyzz, tr_0_0_zz_xyzz_xzzz, tr_0_0_zz_xyzz_yyyy, tr_0_0_zz_xyzz_yyyz, tr_0_0_zz_xyzz_yyzz, tr_0_0_zz_xyzz_yzzz, tr_0_0_zz_xyzz_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxxzz, tr_xyzz_xxxy, tr_xyzz_xxxyzz, tr_xyzz_xxxz, tr_xyzz_xxxzzz, tr_xyzz_xxyy, tr_xyzz_xxyyzz, tr_xyzz_xxyz, tr_xyzz_xxyzzz, tr_xyzz_xxzz, tr_xyzz_xxzzzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyyzz, tr_xyzz_xyyz, tr_xyzz_xyyzzz, tr_xyzz_xyzz, tr_xyzz_xyzzzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_xzzzzz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyyzz, tr_xyzz_yyyz, tr_xyzz_yyyzzz, tr_xyzz_yyzz, tr_xyzz_yyzzzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_yzzzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzz_zzzzzz, tr_xyzzz_xxx, tr_xyzzz_xxxxz, tr_xyzzz_xxxyz, tr_xyzzz_xxxzz, tr_xyzzz_xxy, tr_xyzzz_xxyyz, tr_xyzzz_xxyzz, tr_xyzzz_xxz, tr_xyzzz_xxzzz, tr_xyzzz_xyy, tr_xyzzz_xyyyz, tr_xyzzz_xyyzz, tr_xyzzz_xyz, tr_xyzzz_xyzzz, tr_xyzzz_xzz, tr_xyzzz_xzzzz, tr_xyzzz_yyy, tr_xyzzz_yyyyz, tr_xyzzz_yyyzz, tr_xyzzz_yyz, tr_xyzzz_yyzzz, tr_xyzzz_yzz, tr_xyzzz_yzzzz, tr_xyzzz_zzz, tr_xyzzz_zzzzz, tr_xyzzzz_xxxx, tr_xyzzzz_xxxy, tr_xyzzzz_xxxz, tr_xyzzzz_xxyy, tr_xyzzzz_xxyz, tr_xyzzzz_xxzz, tr_xyzzzz_xyyy, tr_xyzzzz_xyyz, tr_xyzzzz_xyzz, tr_xyzzzz_xzzz, tr_xyzzzz_yyyy, tr_xyzzzz_yyyz, tr_xyzzzz_yyzz, tr_xyzzzz_yzzz, tr_xyzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyzz_xxxx[i] = 2.0 * tr_xy_xxxx[i] - 8.0 * tr_xyz_xxxxz[i] * tke_0 - 10.0 * tr_xyzz_xxxx[i] * tbe_0 - 2.0 * tr_xyzz_xxxx[i] * tke_0 + 4.0 * tr_xyzz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xxxy[i] = 2.0 * tr_xy_xxxy[i] - 8.0 * tr_xyz_xxxyz[i] * tke_0 - 10.0 * tr_xyzz_xxxy[i] * tbe_0 - 2.0 * tr_xyzz_xxxy[i] * tke_0 + 4.0 * tr_xyzz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xxxz[i] = 2.0 * tr_xy_xxxz[i] + 4.0 * tr_xyz_xxx[i] - 8.0 * tr_xyz_xxxzz[i] * tke_0 - 10.0 * tr_xyzz_xxxz[i] * tbe_0 - 6.0 * tr_xyzz_xxxz[i] * tke_0 + 4.0 * tr_xyzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xxx[i] * tbe_0 + 8.0 * tr_xyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xxyy[i] = 2.0 * tr_xy_xxyy[i] - 8.0 * tr_xyz_xxyyz[i] * tke_0 - 10.0 * tr_xyzz_xxyy[i] * tbe_0 - 2.0 * tr_xyzz_xxyy[i] * tke_0 + 4.0 * tr_xyzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xxyz[i] = 2.0 * tr_xy_xxyz[i] + 4.0 * tr_xyz_xxy[i] - 8.0 * tr_xyz_xxyzz[i] * tke_0 - 10.0 * tr_xyzz_xxyz[i] * tbe_0 - 6.0 * tr_xyzz_xxyz[i] * tke_0 + 4.0 * tr_xyzz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xxy[i] * tbe_0 + 8.0 * tr_xyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xxzz[i] = 2.0 * tr_xy_xxzz[i] + 8.0 * tr_xyz_xxz[i] - 8.0 * tr_xyz_xxzzz[i] * tke_0 + 2.0 * tr_xyzz_xx[i] - 10.0 * tr_xyzz_xxzz[i] * tbe_0 - 10.0 * tr_xyzz_xxzz[i] * tke_0 + 4.0 * tr_xyzz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xxz[i] * tbe_0 + 8.0 * tr_xyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xyyy[i] = 2.0 * tr_xy_xyyy[i] - 8.0 * tr_xyz_xyyyz[i] * tke_0 - 10.0 * tr_xyzz_xyyy[i] * tbe_0 - 2.0 * tr_xyzz_xyyy[i] * tke_0 + 4.0 * tr_xyzz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xyyz[i] = 2.0 * tr_xy_xyyz[i] + 4.0 * tr_xyz_xyy[i] - 8.0 * tr_xyz_xyyzz[i] * tke_0 - 10.0 * tr_xyzz_xyyz[i] * tbe_0 - 6.0 * tr_xyzz_xyyz[i] * tke_0 + 4.0 * tr_xyzz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_xyy[i] * tbe_0 + 8.0 * tr_xyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xyzz[i] = 2.0 * tr_xy_xyzz[i] + 8.0 * tr_xyz_xyz[i] - 8.0 * tr_xyz_xyzzz[i] * tke_0 + 2.0 * tr_xyzz_xy[i] - 10.0 * tr_xyzz_xyzz[i] * tbe_0 - 10.0 * tr_xyzz_xyzz[i] * tke_0 + 4.0 * tr_xyzz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_xyz[i] * tbe_0 + 8.0 * tr_xyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xzzz[i] = 2.0 * tr_xy_xzzz[i] + 12.0 * tr_xyz_xzz[i] - 8.0 * tr_xyz_xzzzz[i] * tke_0 + 6.0 * tr_xyzz_xz[i] - 10.0 * tr_xyzz_xzzz[i] * tbe_0 - 14.0 * tr_xyzz_xzzz[i] * tke_0 + 4.0 * tr_xyzz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyzzz_xzz[i] * tbe_0 + 8.0 * tr_xyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_yyyy[i] = 2.0 * tr_xy_yyyy[i] - 8.0 * tr_xyz_yyyyz[i] * tke_0 - 10.0 * tr_xyzz_yyyy[i] * tbe_0 - 2.0 * tr_xyzz_yyyy[i] * tke_0 + 4.0 * tr_xyzz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_yyyz[i] = 2.0 * tr_xy_yyyz[i] + 4.0 * tr_xyz_yyy[i] - 8.0 * tr_xyz_yyyzz[i] * tke_0 - 10.0 * tr_xyzz_yyyz[i] * tbe_0 - 6.0 * tr_xyzz_yyyz[i] * tke_0 + 4.0 * tr_xyzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_yyy[i] * tbe_0 + 8.0 * tr_xyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_yyzz[i] = 2.0 * tr_xy_yyzz[i] + 8.0 * tr_xyz_yyz[i] - 8.0 * tr_xyz_yyzzz[i] * tke_0 + 2.0 * tr_xyzz_yy[i] - 10.0 * tr_xyzz_yyzz[i] * tbe_0 - 10.0 * tr_xyzz_yyzz[i] * tke_0 + 4.0 * tr_xyzz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_yyz[i] * tbe_0 + 8.0 * tr_xyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_yzzz[i] = 2.0 * tr_xy_yzzz[i] + 12.0 * tr_xyz_yzz[i] - 8.0 * tr_xyz_yzzzz[i] * tke_0 + 6.0 * tr_xyzz_yz[i] - 10.0 * tr_xyzz_yzzz[i] * tbe_0 - 14.0 * tr_xyzz_yzzz[i] * tke_0 + 4.0 * tr_xyzz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xyzzz_yzz[i] * tbe_0 + 8.0 * tr_xyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_zzzz[i] = 2.0 * tr_xy_zzzz[i] + 16.0 * tr_xyz_zzz[i] - 8.0 * tr_xyz_zzzzz[i] * tke_0 + 12.0 * tr_xyzz_zz[i] - 10.0 * tr_xyzz_zzzz[i] * tbe_0 - 18.0 * tr_xyzz_zzzz[i] * tke_0 + 4.0 * tr_xyzz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xyzzz_zzz[i] * tbe_0 + 8.0 * tr_xyzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1260-1275 components of targeted buffer : GG

    auto tr_0_0_zz_xzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1260);

    auto tr_0_0_zz_xzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1261);

    auto tr_0_0_zz_xzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1262);

    auto tr_0_0_zz_xzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1263);

    auto tr_0_0_zz_xzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1264);

    auto tr_0_0_zz_xzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1265);

    auto tr_0_0_zz_xzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1266);

    auto tr_0_0_zz_xzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1267);

    auto tr_0_0_zz_xzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1268);

    auto tr_0_0_zz_xzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1269);

    auto tr_0_0_zz_xzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1270);

    auto tr_0_0_zz_xzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1271);

    auto tr_0_0_zz_xzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1272);

    auto tr_0_0_zz_xzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1273);

    auto tr_0_0_zz_xzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1274);

    #pragma omp simd aligned(tr_0_0_zz_xzzz_xxxx, tr_0_0_zz_xzzz_xxxy, tr_0_0_zz_xzzz_xxxz, tr_0_0_zz_xzzz_xxyy, tr_0_0_zz_xzzz_xxyz, tr_0_0_zz_xzzz_xxzz, tr_0_0_zz_xzzz_xyyy, tr_0_0_zz_xzzz_xyyz, tr_0_0_zz_xzzz_xyzz, tr_0_0_zz_xzzz_xzzz, tr_0_0_zz_xzzz_yyyy, tr_0_0_zz_xzzz_yyyz, tr_0_0_zz_xzzz_yyzz, tr_0_0_zz_xzzz_yzzz, tr_0_0_zz_xzzz_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxxzz, tr_xzzz_xxxy, tr_xzzz_xxxyzz, tr_xzzz_xxxz, tr_xzzz_xxxzzz, tr_xzzz_xxyy, tr_xzzz_xxyyzz, tr_xzzz_xxyz, tr_xzzz_xxyzzz, tr_xzzz_xxzz, tr_xzzz_xxzzzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyyzz, tr_xzzz_xyyz, tr_xzzz_xyyzzz, tr_xzzz_xyzz, tr_xzzz_xyzzzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_xzzzzz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyyzz, tr_xzzz_yyyz, tr_xzzz_yyyzzz, tr_xzzz_yyzz, tr_xzzz_yyzzzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_yzzzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzz_zzzzzz, tr_xzzzz_xxx, tr_xzzzz_xxxxz, tr_xzzzz_xxxyz, tr_xzzzz_xxxzz, tr_xzzzz_xxy, tr_xzzzz_xxyyz, tr_xzzzz_xxyzz, tr_xzzzz_xxz, tr_xzzzz_xxzzz, tr_xzzzz_xyy, tr_xzzzz_xyyyz, tr_xzzzz_xyyzz, tr_xzzzz_xyz, tr_xzzzz_xyzzz, tr_xzzzz_xzz, tr_xzzzz_xzzzz, tr_xzzzz_yyy, tr_xzzzz_yyyyz, tr_xzzzz_yyyzz, tr_xzzzz_yyz, tr_xzzzz_yyzzz, tr_xzzzz_yzz, tr_xzzzz_yzzzz, tr_xzzzz_zzz, tr_xzzzz_zzzzz, tr_xzzzzz_xxxx, tr_xzzzzz_xxxy, tr_xzzzzz_xxxz, tr_xzzzzz_xxyy, tr_xzzzzz_xxyz, tr_xzzzzz_xxzz, tr_xzzzzz_xyyy, tr_xzzzzz_xyyz, tr_xzzzzz_xyzz, tr_xzzzzz_xzzz, tr_xzzzzz_yyyy, tr_xzzzzz_yyyz, tr_xzzzzz_yyzz, tr_xzzzzz_yzzz, tr_xzzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xzzz_xxxx[i] = 6.0 * tr_xz_xxxx[i] - 12.0 * tr_xzz_xxxxz[i] * tke_0 - 14.0 * tr_xzzz_xxxx[i] * tbe_0 - 2.0 * tr_xzzz_xxxx[i] * tke_0 + 4.0 * tr_xzzz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xxxy[i] = 6.0 * tr_xz_xxxy[i] - 12.0 * tr_xzz_xxxyz[i] * tke_0 - 14.0 * tr_xzzz_xxxy[i] * tbe_0 - 2.0 * tr_xzzz_xxxy[i] * tke_0 + 4.0 * tr_xzzz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xxxz[i] = 6.0 * tr_xz_xxxz[i] + 6.0 * tr_xzz_xxx[i] - 12.0 * tr_xzz_xxxzz[i] * tke_0 - 14.0 * tr_xzzz_xxxz[i] * tbe_0 - 6.0 * tr_xzzz_xxxz[i] * tke_0 + 4.0 * tr_xzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xxx[i] * tbe_0 + 8.0 * tr_xzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xxyy[i] = 6.0 * tr_xz_xxyy[i] - 12.0 * tr_xzz_xxyyz[i] * tke_0 - 14.0 * tr_xzzz_xxyy[i] * tbe_0 - 2.0 * tr_xzzz_xxyy[i] * tke_0 + 4.0 * tr_xzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xxyz[i] = 6.0 * tr_xz_xxyz[i] + 6.0 * tr_xzz_xxy[i] - 12.0 * tr_xzz_xxyzz[i] * tke_0 - 14.0 * tr_xzzz_xxyz[i] * tbe_0 - 6.0 * tr_xzzz_xxyz[i] * tke_0 + 4.0 * tr_xzzz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xxy[i] * tbe_0 + 8.0 * tr_xzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xxzz[i] = 6.0 * tr_xz_xxzz[i] + 12.0 * tr_xzz_xxz[i] - 12.0 * tr_xzz_xxzzz[i] * tke_0 + 2.0 * tr_xzzz_xx[i] - 14.0 * tr_xzzz_xxzz[i] * tbe_0 - 10.0 * tr_xzzz_xxzz[i] * tke_0 + 4.0 * tr_xzzz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_xxz[i] * tbe_0 + 8.0 * tr_xzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xyyy[i] = 6.0 * tr_xz_xyyy[i] - 12.0 * tr_xzz_xyyyz[i] * tke_0 - 14.0 * tr_xzzz_xyyy[i] * tbe_0 - 2.0 * tr_xzzz_xyyy[i] * tke_0 + 4.0 * tr_xzzz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xyyz[i] = 6.0 * tr_xz_xyyz[i] + 6.0 * tr_xzz_xyy[i] - 12.0 * tr_xzz_xyyzz[i] * tke_0 - 14.0 * tr_xzzz_xyyz[i] * tbe_0 - 6.0 * tr_xzzz_xyyz[i] * tke_0 + 4.0 * tr_xzzz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_xyy[i] * tbe_0 + 8.0 * tr_xzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xyzz[i] = 6.0 * tr_xz_xyzz[i] + 12.0 * tr_xzz_xyz[i] - 12.0 * tr_xzz_xyzzz[i] * tke_0 + 2.0 * tr_xzzz_xy[i] - 14.0 * tr_xzzz_xyzz[i] * tbe_0 - 10.0 * tr_xzzz_xyzz[i] * tke_0 + 4.0 * tr_xzzz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_xyz[i] * tbe_0 + 8.0 * tr_xzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xzzz[i] = 6.0 * tr_xz_xzzz[i] + 18.0 * tr_xzz_xzz[i] - 12.0 * tr_xzz_xzzzz[i] * tke_0 + 6.0 * tr_xzzz_xz[i] - 14.0 * tr_xzzz_xzzz[i] * tbe_0 - 14.0 * tr_xzzz_xzzz[i] * tke_0 + 4.0 * tr_xzzz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xzzzz_xzz[i] * tbe_0 + 8.0 * tr_xzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_yyyy[i] = 6.0 * tr_xz_yyyy[i] - 12.0 * tr_xzz_yyyyz[i] * tke_0 - 14.0 * tr_xzzz_yyyy[i] * tbe_0 - 2.0 * tr_xzzz_yyyy[i] * tke_0 + 4.0 * tr_xzzz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_yyyz[i] = 6.0 * tr_xz_yyyz[i] + 6.0 * tr_xzz_yyy[i] - 12.0 * tr_xzz_yyyzz[i] * tke_0 - 14.0 * tr_xzzz_yyyz[i] * tbe_0 - 6.0 * tr_xzzz_yyyz[i] * tke_0 + 4.0 * tr_xzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_yyy[i] * tbe_0 + 8.0 * tr_xzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_yyzz[i] = 6.0 * tr_xz_yyzz[i] + 12.0 * tr_xzz_yyz[i] - 12.0 * tr_xzz_yyzzz[i] * tke_0 + 2.0 * tr_xzzz_yy[i] - 14.0 * tr_xzzz_yyzz[i] * tbe_0 - 10.0 * tr_xzzz_yyzz[i] * tke_0 + 4.0 * tr_xzzz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_yyz[i] * tbe_0 + 8.0 * tr_xzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_yzzz[i] = 6.0 * tr_xz_yzzz[i] + 18.0 * tr_xzz_yzz[i] - 12.0 * tr_xzz_yzzzz[i] * tke_0 + 6.0 * tr_xzzz_yz[i] - 14.0 * tr_xzzz_yzzz[i] * tbe_0 - 14.0 * tr_xzzz_yzzz[i] * tke_0 + 4.0 * tr_xzzz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_xzzzz_yzz[i] * tbe_0 + 8.0 * tr_xzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_zzzz[i] = 6.0 * tr_xz_zzzz[i] + 24.0 * tr_xzz_zzz[i] - 12.0 * tr_xzz_zzzzz[i] * tke_0 + 12.0 * tr_xzzz_zz[i] - 14.0 * tr_xzzz_zzzz[i] * tbe_0 - 18.0 * tr_xzzz_zzzz[i] * tke_0 + 4.0 * tr_xzzz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_xzzzz_zzz[i] * tbe_0 + 8.0 * tr_xzzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1275-1290 components of targeted buffer : GG

    auto tr_0_0_zz_yyyy_xxxx = pbuffer.data(idx_op_geom_020_gg + 1275);

    auto tr_0_0_zz_yyyy_xxxy = pbuffer.data(idx_op_geom_020_gg + 1276);

    auto tr_0_0_zz_yyyy_xxxz = pbuffer.data(idx_op_geom_020_gg + 1277);

    auto tr_0_0_zz_yyyy_xxyy = pbuffer.data(idx_op_geom_020_gg + 1278);

    auto tr_0_0_zz_yyyy_xxyz = pbuffer.data(idx_op_geom_020_gg + 1279);

    auto tr_0_0_zz_yyyy_xxzz = pbuffer.data(idx_op_geom_020_gg + 1280);

    auto tr_0_0_zz_yyyy_xyyy = pbuffer.data(idx_op_geom_020_gg + 1281);

    auto tr_0_0_zz_yyyy_xyyz = pbuffer.data(idx_op_geom_020_gg + 1282);

    auto tr_0_0_zz_yyyy_xyzz = pbuffer.data(idx_op_geom_020_gg + 1283);

    auto tr_0_0_zz_yyyy_xzzz = pbuffer.data(idx_op_geom_020_gg + 1284);

    auto tr_0_0_zz_yyyy_yyyy = pbuffer.data(idx_op_geom_020_gg + 1285);

    auto tr_0_0_zz_yyyy_yyyz = pbuffer.data(idx_op_geom_020_gg + 1286);

    auto tr_0_0_zz_yyyy_yyzz = pbuffer.data(idx_op_geom_020_gg + 1287);

    auto tr_0_0_zz_yyyy_yzzz = pbuffer.data(idx_op_geom_020_gg + 1288);

    auto tr_0_0_zz_yyyy_zzzz = pbuffer.data(idx_op_geom_020_gg + 1289);

    #pragma omp simd aligned(tr_0_0_zz_yyyy_xxxx, tr_0_0_zz_yyyy_xxxy, tr_0_0_zz_yyyy_xxxz, tr_0_0_zz_yyyy_xxyy, tr_0_0_zz_yyyy_xxyz, tr_0_0_zz_yyyy_xxzz, tr_0_0_zz_yyyy_xyyy, tr_0_0_zz_yyyy_xyyz, tr_0_0_zz_yyyy_xyzz, tr_0_0_zz_yyyy_xzzz, tr_0_0_zz_yyyy_yyyy, tr_0_0_zz_yyyy_yyyz, tr_0_0_zz_yyyy_yyzz, tr_0_0_zz_yyyy_yzzz, tr_0_0_zz_yyyy_zzzz, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxxzz, tr_yyyy_xxxy, tr_yyyy_xxxyzz, tr_yyyy_xxxz, tr_yyyy_xxxzzz, tr_yyyy_xxyy, tr_yyyy_xxyyzz, tr_yyyy_xxyz, tr_yyyy_xxyzzz, tr_yyyy_xxzz, tr_yyyy_xxzzzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyyzz, tr_yyyy_xyyz, tr_yyyy_xyyzzz, tr_yyyy_xyzz, tr_yyyy_xyzzzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_xzzzzz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyyzz, tr_yyyy_yyyz, tr_yyyy_yyyzzz, tr_yyyy_yyzz, tr_yyyy_yyzzzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_yzzzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyy_zzzzzz, tr_yyyyz_xxx, tr_yyyyz_xxxxz, tr_yyyyz_xxxyz, tr_yyyyz_xxxzz, tr_yyyyz_xxy, tr_yyyyz_xxyyz, tr_yyyyz_xxyzz, tr_yyyyz_xxz, tr_yyyyz_xxzzz, tr_yyyyz_xyy, tr_yyyyz_xyyyz, tr_yyyyz_xyyzz, tr_yyyyz_xyz, tr_yyyyz_xyzzz, tr_yyyyz_xzz, tr_yyyyz_xzzzz, tr_yyyyz_yyy, tr_yyyyz_yyyyz, tr_yyyyz_yyyzz, tr_yyyyz_yyz, tr_yyyyz_yyzzz, tr_yyyyz_yzz, tr_yyyyz_yzzzz, tr_yyyyz_zzz, tr_yyyyz_zzzzz, tr_yyyyzz_xxxx, tr_yyyyzz_xxxy, tr_yyyyzz_xxxz, tr_yyyyzz_xxyy, tr_yyyyzz_xxyz, tr_yyyyzz_xxzz, tr_yyyyzz_xyyy, tr_yyyyzz_xyyz, tr_yyyyzz_xyzz, tr_yyyyzz_xzzz, tr_yyyyzz_yyyy, tr_yyyyzz_yyyz, tr_yyyyzz_yyzz, tr_yyyyzz_yzzz, tr_yyyyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyyy_xxxx[i] = -2.0 * tr_yyyy_xxxx[i] * tbe_0 - 2.0 * tr_yyyy_xxxx[i] * tke_0 + 4.0 * tr_yyyy_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xxxy[i] = -2.0 * tr_yyyy_xxxy[i] * tbe_0 - 2.0 * tr_yyyy_xxxy[i] * tke_0 + 4.0 * tr_yyyy_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xxxz[i] = -2.0 * tr_yyyy_xxxz[i] * tbe_0 - 6.0 * tr_yyyy_xxxz[i] * tke_0 + 4.0 * tr_yyyy_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xxx[i] * tbe_0 + 8.0 * tr_yyyyz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xxyy[i] = -2.0 * tr_yyyy_xxyy[i] * tbe_0 - 2.0 * tr_yyyy_xxyy[i] * tke_0 + 4.0 * tr_yyyy_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xxyz[i] = -2.0 * tr_yyyy_xxyz[i] * tbe_0 - 6.0 * tr_yyyy_xxyz[i] * tke_0 + 4.0 * tr_yyyy_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xxy[i] * tbe_0 + 8.0 * tr_yyyyz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xxzz[i] = 2.0 * tr_yyyy_xx[i] - 2.0 * tr_yyyy_xxzz[i] * tbe_0 - 10.0 * tr_yyyy_xxzz[i] * tke_0 + 4.0 * tr_yyyy_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_xxz[i] * tbe_0 + 8.0 * tr_yyyyz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xyyy[i] = -2.0 * tr_yyyy_xyyy[i] * tbe_0 - 2.0 * tr_yyyy_xyyy[i] * tke_0 + 4.0 * tr_yyyy_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xyyz[i] = -2.0 * tr_yyyy_xyyz[i] * tbe_0 - 6.0 * tr_yyyy_xyyz[i] * tke_0 + 4.0 * tr_yyyy_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_xyy[i] * tbe_0 + 8.0 * tr_yyyyz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xyzz[i] = 2.0 * tr_yyyy_xy[i] - 2.0 * tr_yyyy_xyzz[i] * tbe_0 - 10.0 * tr_yyyy_xyzz[i] * tke_0 + 4.0 * tr_yyyy_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_xyz[i] * tbe_0 + 8.0 * tr_yyyyz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xzzz[i] = 6.0 * tr_yyyy_xz[i] - 2.0 * tr_yyyy_xzzz[i] * tbe_0 - 14.0 * tr_yyyy_xzzz[i] * tke_0 + 4.0 * tr_yyyy_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyyyz_xzz[i] * tbe_0 + 8.0 * tr_yyyyz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_yyyy[i] = -2.0 * tr_yyyy_yyyy[i] * tbe_0 - 2.0 * tr_yyyy_yyyy[i] * tke_0 + 4.0 * tr_yyyy_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_yyyz[i] = -2.0 * tr_yyyy_yyyz[i] * tbe_0 - 6.0 * tr_yyyy_yyyz[i] * tke_0 + 4.0 * tr_yyyy_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_yyy[i] * tbe_0 + 8.0 * tr_yyyyz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_yyzz[i] = 2.0 * tr_yyyy_yy[i] - 2.0 * tr_yyyy_yyzz[i] * tbe_0 - 10.0 * tr_yyyy_yyzz[i] * tke_0 + 4.0 * tr_yyyy_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_yyz[i] * tbe_0 + 8.0 * tr_yyyyz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_yzzz[i] = 6.0 * tr_yyyy_yz[i] - 2.0 * tr_yyyy_yzzz[i] * tbe_0 - 14.0 * tr_yyyy_yzzz[i] * tke_0 + 4.0 * tr_yyyy_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyyyz_yzz[i] * tbe_0 + 8.0 * tr_yyyyz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_zzzz[i] = 12.0 * tr_yyyy_zz[i] - 2.0 * tr_yyyy_zzzz[i] * tbe_0 - 18.0 * tr_yyyy_zzzz[i] * tke_0 + 4.0 * tr_yyyy_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_yyyyz_zzz[i] * tbe_0 + 8.0 * tr_yyyyz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1290-1305 components of targeted buffer : GG

    auto tr_0_0_zz_yyyz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1290);

    auto tr_0_0_zz_yyyz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1291);

    auto tr_0_0_zz_yyyz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1292);

    auto tr_0_0_zz_yyyz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1293);

    auto tr_0_0_zz_yyyz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1294);

    auto tr_0_0_zz_yyyz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1295);

    auto tr_0_0_zz_yyyz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1296);

    auto tr_0_0_zz_yyyz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1297);

    auto tr_0_0_zz_yyyz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1298);

    auto tr_0_0_zz_yyyz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1299);

    auto tr_0_0_zz_yyyz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1300);

    auto tr_0_0_zz_yyyz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1301);

    auto tr_0_0_zz_yyyz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1302);

    auto tr_0_0_zz_yyyz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1303);

    auto tr_0_0_zz_yyyz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1304);

    #pragma omp simd aligned(tr_0_0_zz_yyyz_xxxx, tr_0_0_zz_yyyz_xxxy, tr_0_0_zz_yyyz_xxxz, tr_0_0_zz_yyyz_xxyy, tr_0_0_zz_yyyz_xxyz, tr_0_0_zz_yyyz_xxzz, tr_0_0_zz_yyyz_xyyy, tr_0_0_zz_yyyz_xyyz, tr_0_0_zz_yyyz_xyzz, tr_0_0_zz_yyyz_xzzz, tr_0_0_zz_yyyz_yyyy, tr_0_0_zz_yyyz_yyyz, tr_0_0_zz_yyyz_yyzz, tr_0_0_zz_yyyz_yzzz, tr_0_0_zz_yyyz_zzzz, tr_yyy_xxx, tr_yyy_xxxxz, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyy_zzzzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxxzz, tr_yyyz_xxxy, tr_yyyz_xxxyzz, tr_yyyz_xxxz, tr_yyyz_xxxzzz, tr_yyyz_xxyy, tr_yyyz_xxyyzz, tr_yyyz_xxyz, tr_yyyz_xxyzzz, tr_yyyz_xxzz, tr_yyyz_xxzzzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyyzz, tr_yyyz_xyyz, tr_yyyz_xyyzzz, tr_yyyz_xyzz, tr_yyyz_xyzzzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_xzzzzz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyyzz, tr_yyyz_yyyz, tr_yyyz_yyyzzz, tr_yyyz_yyzz, tr_yyyz_yyzzzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_yzzzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyz_zzzzzz, tr_yyyzz_xxx, tr_yyyzz_xxxxz, tr_yyyzz_xxxyz, tr_yyyzz_xxxzz, tr_yyyzz_xxy, tr_yyyzz_xxyyz, tr_yyyzz_xxyzz, tr_yyyzz_xxz, tr_yyyzz_xxzzz, tr_yyyzz_xyy, tr_yyyzz_xyyyz, tr_yyyzz_xyyzz, tr_yyyzz_xyz, tr_yyyzz_xyzzz, tr_yyyzz_xzz, tr_yyyzz_xzzzz, tr_yyyzz_yyy, tr_yyyzz_yyyyz, tr_yyyzz_yyyzz, tr_yyyzz_yyz, tr_yyyzz_yyzzz, tr_yyyzz_yzz, tr_yyyzz_yzzzz, tr_yyyzz_zzz, tr_yyyzz_zzzzz, tr_yyyzzz_xxxx, tr_yyyzzz_xxxy, tr_yyyzzz_xxxz, tr_yyyzzz_xxyy, tr_yyyzzz_xxyz, tr_yyyzzz_xxzz, tr_yyyzzz_xyyy, tr_yyyzzz_xyyz, tr_yyyzzz_xyzz, tr_yyyzzz_xzzz, tr_yyyzzz_yyyy, tr_yyyzzz_yyyz, tr_yyyzzz_yyzz, tr_yyyzzz_yzzz, tr_yyyzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyyz_xxxx[i] = -4.0 * tr_yyy_xxxxz[i] * tke_0 - 6.0 * tr_yyyz_xxxx[i] * tbe_0 - 2.0 * tr_yyyz_xxxx[i] * tke_0 + 4.0 * tr_yyyz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xxxy[i] = -4.0 * tr_yyy_xxxyz[i] * tke_0 - 6.0 * tr_yyyz_xxxy[i] * tbe_0 - 2.0 * tr_yyyz_xxxy[i] * tke_0 + 4.0 * tr_yyyz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xxxz[i] = 2.0 * tr_yyy_xxx[i] - 4.0 * tr_yyy_xxxzz[i] * tke_0 - 6.0 * tr_yyyz_xxxz[i] * tbe_0 - 6.0 * tr_yyyz_xxxz[i] * tke_0 + 4.0 * tr_yyyz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xxx[i] * tbe_0 + 8.0 * tr_yyyzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xxyy[i] = -4.0 * tr_yyy_xxyyz[i] * tke_0 - 6.0 * tr_yyyz_xxyy[i] * tbe_0 - 2.0 * tr_yyyz_xxyy[i] * tke_0 + 4.0 * tr_yyyz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xxyz[i] = 2.0 * tr_yyy_xxy[i] - 4.0 * tr_yyy_xxyzz[i] * tke_0 - 6.0 * tr_yyyz_xxyz[i] * tbe_0 - 6.0 * tr_yyyz_xxyz[i] * tke_0 + 4.0 * tr_yyyz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xxy[i] * tbe_0 + 8.0 * tr_yyyzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xxzz[i] = 4.0 * tr_yyy_xxz[i] - 4.0 * tr_yyy_xxzzz[i] * tke_0 + 2.0 * tr_yyyz_xx[i] - 6.0 * tr_yyyz_xxzz[i] * tbe_0 - 10.0 * tr_yyyz_xxzz[i] * tke_0 + 4.0 * tr_yyyz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_xxz[i] * tbe_0 + 8.0 * tr_yyyzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xyyy[i] = -4.0 * tr_yyy_xyyyz[i] * tke_0 - 6.0 * tr_yyyz_xyyy[i] * tbe_0 - 2.0 * tr_yyyz_xyyy[i] * tke_0 + 4.0 * tr_yyyz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xyyz[i] = 2.0 * tr_yyy_xyy[i] - 4.0 * tr_yyy_xyyzz[i] * tke_0 - 6.0 * tr_yyyz_xyyz[i] * tbe_0 - 6.0 * tr_yyyz_xyyz[i] * tke_0 + 4.0 * tr_yyyz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_xyy[i] * tbe_0 + 8.0 * tr_yyyzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xyzz[i] = 4.0 * tr_yyy_xyz[i] - 4.0 * tr_yyy_xyzzz[i] * tke_0 + 2.0 * tr_yyyz_xy[i] - 6.0 * tr_yyyz_xyzz[i] * tbe_0 - 10.0 * tr_yyyz_xyzz[i] * tke_0 + 4.0 * tr_yyyz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_xyz[i] * tbe_0 + 8.0 * tr_yyyzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xzzz[i] = 6.0 * tr_yyy_xzz[i] - 4.0 * tr_yyy_xzzzz[i] * tke_0 + 6.0 * tr_yyyz_xz[i] - 6.0 * tr_yyyz_xzzz[i] * tbe_0 - 14.0 * tr_yyyz_xzzz[i] * tke_0 + 4.0 * tr_yyyz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyyzz_xzz[i] * tbe_0 + 8.0 * tr_yyyzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_yyyy[i] = -4.0 * tr_yyy_yyyyz[i] * tke_0 - 6.0 * tr_yyyz_yyyy[i] * tbe_0 - 2.0 * tr_yyyz_yyyy[i] * tke_0 + 4.0 * tr_yyyz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_yyyz[i] = 2.0 * tr_yyy_yyy[i] - 4.0 * tr_yyy_yyyzz[i] * tke_0 - 6.0 * tr_yyyz_yyyz[i] * tbe_0 - 6.0 * tr_yyyz_yyyz[i] * tke_0 + 4.0 * tr_yyyz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_yyy[i] * tbe_0 + 8.0 * tr_yyyzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_yyzz[i] = 4.0 * tr_yyy_yyz[i] - 4.0 * tr_yyy_yyzzz[i] * tke_0 + 2.0 * tr_yyyz_yy[i] - 6.0 * tr_yyyz_yyzz[i] * tbe_0 - 10.0 * tr_yyyz_yyzz[i] * tke_0 + 4.0 * tr_yyyz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_yyz[i] * tbe_0 + 8.0 * tr_yyyzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_yzzz[i] = 6.0 * tr_yyy_yzz[i] - 4.0 * tr_yyy_yzzzz[i] * tke_0 + 6.0 * tr_yyyz_yz[i] - 6.0 * tr_yyyz_yzzz[i] * tbe_0 - 14.0 * tr_yyyz_yzzz[i] * tke_0 + 4.0 * tr_yyyz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyyzz_yzz[i] * tbe_0 + 8.0 * tr_yyyzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_zzzz[i] = 8.0 * tr_yyy_zzz[i] - 4.0 * tr_yyy_zzzzz[i] * tke_0 + 12.0 * tr_yyyz_zz[i] - 6.0 * tr_yyyz_zzzz[i] * tbe_0 - 18.0 * tr_yyyz_zzzz[i] * tke_0 + 4.0 * tr_yyyz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_yyyzz_zzz[i] * tbe_0 + 8.0 * tr_yyyzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1305-1320 components of targeted buffer : GG

    auto tr_0_0_zz_yyzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1305);

    auto tr_0_0_zz_yyzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1306);

    auto tr_0_0_zz_yyzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1307);

    auto tr_0_0_zz_yyzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1308);

    auto tr_0_0_zz_yyzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1309);

    auto tr_0_0_zz_yyzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1310);

    auto tr_0_0_zz_yyzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1311);

    auto tr_0_0_zz_yyzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1312);

    auto tr_0_0_zz_yyzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1313);

    auto tr_0_0_zz_yyzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1314);

    auto tr_0_0_zz_yyzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1315);

    auto tr_0_0_zz_yyzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1316);

    auto tr_0_0_zz_yyzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1317);

    auto tr_0_0_zz_yyzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1318);

    auto tr_0_0_zz_yyzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1319);

    #pragma omp simd aligned(tr_0_0_zz_yyzz_xxxx, tr_0_0_zz_yyzz_xxxy, tr_0_0_zz_yyzz_xxxz, tr_0_0_zz_yyzz_xxyy, tr_0_0_zz_yyzz_xxyz, tr_0_0_zz_yyzz_xxzz, tr_0_0_zz_yyzz_xyyy, tr_0_0_zz_yyzz_xyyz, tr_0_0_zz_yyzz_xyzz, tr_0_0_zz_yyzz_xzzz, tr_0_0_zz_yyzz_yyyy, tr_0_0_zz_yyzz_yyyz, tr_0_0_zz_yyzz_yyzz, tr_0_0_zz_yyzz_yzzz, tr_0_0_zz_yyzz_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxxzz, tr_yyzz_xxxy, tr_yyzz_xxxyzz, tr_yyzz_xxxz, tr_yyzz_xxxzzz, tr_yyzz_xxyy, tr_yyzz_xxyyzz, tr_yyzz_xxyz, tr_yyzz_xxyzzz, tr_yyzz_xxzz, tr_yyzz_xxzzzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyyzz, tr_yyzz_xyyz, tr_yyzz_xyyzzz, tr_yyzz_xyzz, tr_yyzz_xyzzzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_xzzzzz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyyzz, tr_yyzz_yyyz, tr_yyzz_yyyzzz, tr_yyzz_yyzz, tr_yyzz_yyzzzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_yzzzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzz_zzzzzz, tr_yyzzz_xxx, tr_yyzzz_xxxxz, tr_yyzzz_xxxyz, tr_yyzzz_xxxzz, tr_yyzzz_xxy, tr_yyzzz_xxyyz, tr_yyzzz_xxyzz, tr_yyzzz_xxz, tr_yyzzz_xxzzz, tr_yyzzz_xyy, tr_yyzzz_xyyyz, tr_yyzzz_xyyzz, tr_yyzzz_xyz, tr_yyzzz_xyzzz, tr_yyzzz_xzz, tr_yyzzz_xzzzz, tr_yyzzz_yyy, tr_yyzzz_yyyyz, tr_yyzzz_yyyzz, tr_yyzzz_yyz, tr_yyzzz_yyzzz, tr_yyzzz_yzz, tr_yyzzz_yzzzz, tr_yyzzz_zzz, tr_yyzzz_zzzzz, tr_yyzzzz_xxxx, tr_yyzzzz_xxxy, tr_yyzzzz_xxxz, tr_yyzzzz_xxyy, tr_yyzzzz_xxyz, tr_yyzzzz_xxzz, tr_yyzzzz_xyyy, tr_yyzzzz_xyyz, tr_yyzzzz_xyzz, tr_yyzzzz_xzzz, tr_yyzzzz_yyyy, tr_yyzzzz_yyyz, tr_yyzzzz_yyzz, tr_yyzzzz_yzzz, tr_yyzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyzz_xxxx[i] = 2.0 * tr_yy_xxxx[i] - 8.0 * tr_yyz_xxxxz[i] * tke_0 - 10.0 * tr_yyzz_xxxx[i] * tbe_0 - 2.0 * tr_yyzz_xxxx[i] * tke_0 + 4.0 * tr_yyzz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xxxy[i] = 2.0 * tr_yy_xxxy[i] - 8.0 * tr_yyz_xxxyz[i] * tke_0 - 10.0 * tr_yyzz_xxxy[i] * tbe_0 - 2.0 * tr_yyzz_xxxy[i] * tke_0 + 4.0 * tr_yyzz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xxxz[i] = 2.0 * tr_yy_xxxz[i] + 4.0 * tr_yyz_xxx[i] - 8.0 * tr_yyz_xxxzz[i] * tke_0 - 10.0 * tr_yyzz_xxxz[i] * tbe_0 - 6.0 * tr_yyzz_xxxz[i] * tke_0 + 4.0 * tr_yyzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xxx[i] * tbe_0 + 8.0 * tr_yyzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xxyy[i] = 2.0 * tr_yy_xxyy[i] - 8.0 * tr_yyz_xxyyz[i] * tke_0 - 10.0 * tr_yyzz_xxyy[i] * tbe_0 - 2.0 * tr_yyzz_xxyy[i] * tke_0 + 4.0 * tr_yyzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xxyz[i] = 2.0 * tr_yy_xxyz[i] + 4.0 * tr_yyz_xxy[i] - 8.0 * tr_yyz_xxyzz[i] * tke_0 - 10.0 * tr_yyzz_xxyz[i] * tbe_0 - 6.0 * tr_yyzz_xxyz[i] * tke_0 + 4.0 * tr_yyzz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xxy[i] * tbe_0 + 8.0 * tr_yyzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xxzz[i] = 2.0 * tr_yy_xxzz[i] + 8.0 * tr_yyz_xxz[i] - 8.0 * tr_yyz_xxzzz[i] * tke_0 + 2.0 * tr_yyzz_xx[i] - 10.0 * tr_yyzz_xxzz[i] * tbe_0 - 10.0 * tr_yyzz_xxzz[i] * tke_0 + 4.0 * tr_yyzz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_xxz[i] * tbe_0 + 8.0 * tr_yyzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xyyy[i] = 2.0 * tr_yy_xyyy[i] - 8.0 * tr_yyz_xyyyz[i] * tke_0 - 10.0 * tr_yyzz_xyyy[i] * tbe_0 - 2.0 * tr_yyzz_xyyy[i] * tke_0 + 4.0 * tr_yyzz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xyyz[i] = 2.0 * tr_yy_xyyz[i] + 4.0 * tr_yyz_xyy[i] - 8.0 * tr_yyz_xyyzz[i] * tke_0 - 10.0 * tr_yyzz_xyyz[i] * tbe_0 - 6.0 * tr_yyzz_xyyz[i] * tke_0 + 4.0 * tr_yyzz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_xyy[i] * tbe_0 + 8.0 * tr_yyzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xyzz[i] = 2.0 * tr_yy_xyzz[i] + 8.0 * tr_yyz_xyz[i] - 8.0 * tr_yyz_xyzzz[i] * tke_0 + 2.0 * tr_yyzz_xy[i] - 10.0 * tr_yyzz_xyzz[i] * tbe_0 - 10.0 * tr_yyzz_xyzz[i] * tke_0 + 4.0 * tr_yyzz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_xyz[i] * tbe_0 + 8.0 * tr_yyzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xzzz[i] = 2.0 * tr_yy_xzzz[i] + 12.0 * tr_yyz_xzz[i] - 8.0 * tr_yyz_xzzzz[i] * tke_0 + 6.0 * tr_yyzz_xz[i] - 10.0 * tr_yyzz_xzzz[i] * tbe_0 - 14.0 * tr_yyzz_xzzz[i] * tke_0 + 4.0 * tr_yyzz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyzzz_xzz[i] * tbe_0 + 8.0 * tr_yyzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_yyyy[i] = 2.0 * tr_yy_yyyy[i] - 8.0 * tr_yyz_yyyyz[i] * tke_0 - 10.0 * tr_yyzz_yyyy[i] * tbe_0 - 2.0 * tr_yyzz_yyyy[i] * tke_0 + 4.0 * tr_yyzz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_yyyz[i] = 2.0 * tr_yy_yyyz[i] + 4.0 * tr_yyz_yyy[i] - 8.0 * tr_yyz_yyyzz[i] * tke_0 - 10.0 * tr_yyzz_yyyz[i] * tbe_0 - 6.0 * tr_yyzz_yyyz[i] * tke_0 + 4.0 * tr_yyzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_yyy[i] * tbe_0 + 8.0 * tr_yyzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_yyzz[i] = 2.0 * tr_yy_yyzz[i] + 8.0 * tr_yyz_yyz[i] - 8.0 * tr_yyz_yyzzz[i] * tke_0 + 2.0 * tr_yyzz_yy[i] - 10.0 * tr_yyzz_yyzz[i] * tbe_0 - 10.0 * tr_yyzz_yyzz[i] * tke_0 + 4.0 * tr_yyzz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_yyz[i] * tbe_0 + 8.0 * tr_yyzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_yzzz[i] = 2.0 * tr_yy_yzzz[i] + 12.0 * tr_yyz_yzz[i] - 8.0 * tr_yyz_yzzzz[i] * tke_0 + 6.0 * tr_yyzz_yz[i] - 10.0 * tr_yyzz_yzzz[i] * tbe_0 - 14.0 * tr_yyzz_yzzz[i] * tke_0 + 4.0 * tr_yyzz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yyzzz_yzz[i] * tbe_0 + 8.0 * tr_yyzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_zzzz[i] = 2.0 * tr_yy_zzzz[i] + 16.0 * tr_yyz_zzz[i] - 8.0 * tr_yyz_zzzzz[i] * tke_0 + 12.0 * tr_yyzz_zz[i] - 10.0 * tr_yyzz_zzzz[i] * tbe_0 - 18.0 * tr_yyzz_zzzz[i] * tke_0 + 4.0 * tr_yyzz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_yyzzz_zzz[i] * tbe_0 + 8.0 * tr_yyzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1320-1335 components of targeted buffer : GG

    auto tr_0_0_zz_yzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1320);

    auto tr_0_0_zz_yzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1321);

    auto tr_0_0_zz_yzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1322);

    auto tr_0_0_zz_yzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1323);

    auto tr_0_0_zz_yzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1324);

    auto tr_0_0_zz_yzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1325);

    auto tr_0_0_zz_yzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1326);

    auto tr_0_0_zz_yzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1327);

    auto tr_0_0_zz_yzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1328);

    auto tr_0_0_zz_yzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1329);

    auto tr_0_0_zz_yzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1330);

    auto tr_0_0_zz_yzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1331);

    auto tr_0_0_zz_yzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1332);

    auto tr_0_0_zz_yzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1333);

    auto tr_0_0_zz_yzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1334);

    #pragma omp simd aligned(tr_0_0_zz_yzzz_xxxx, tr_0_0_zz_yzzz_xxxy, tr_0_0_zz_yzzz_xxxz, tr_0_0_zz_yzzz_xxyy, tr_0_0_zz_yzzz_xxyz, tr_0_0_zz_yzzz_xxzz, tr_0_0_zz_yzzz_xyyy, tr_0_0_zz_yzzz_xyyz, tr_0_0_zz_yzzz_xyzz, tr_0_0_zz_yzzz_xzzz, tr_0_0_zz_yzzz_yyyy, tr_0_0_zz_yzzz_yyyz, tr_0_0_zz_yzzz_yyzz, tr_0_0_zz_yzzz_yzzz, tr_0_0_zz_yzzz_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxxzz, tr_yzzz_xxxy, tr_yzzz_xxxyzz, tr_yzzz_xxxz, tr_yzzz_xxxzzz, tr_yzzz_xxyy, tr_yzzz_xxyyzz, tr_yzzz_xxyz, tr_yzzz_xxyzzz, tr_yzzz_xxzz, tr_yzzz_xxzzzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyyzz, tr_yzzz_xyyz, tr_yzzz_xyyzzz, tr_yzzz_xyzz, tr_yzzz_xyzzzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_xzzzzz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyyzz, tr_yzzz_yyyz, tr_yzzz_yyyzzz, tr_yzzz_yyzz, tr_yzzz_yyzzzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_yzzzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzz_zzzzzz, tr_yzzzz_xxx, tr_yzzzz_xxxxz, tr_yzzzz_xxxyz, tr_yzzzz_xxxzz, tr_yzzzz_xxy, tr_yzzzz_xxyyz, tr_yzzzz_xxyzz, tr_yzzzz_xxz, tr_yzzzz_xxzzz, tr_yzzzz_xyy, tr_yzzzz_xyyyz, tr_yzzzz_xyyzz, tr_yzzzz_xyz, tr_yzzzz_xyzzz, tr_yzzzz_xzz, tr_yzzzz_xzzzz, tr_yzzzz_yyy, tr_yzzzz_yyyyz, tr_yzzzz_yyyzz, tr_yzzzz_yyz, tr_yzzzz_yyzzz, tr_yzzzz_yzz, tr_yzzzz_yzzzz, tr_yzzzz_zzz, tr_yzzzz_zzzzz, tr_yzzzzz_xxxx, tr_yzzzzz_xxxy, tr_yzzzzz_xxxz, tr_yzzzzz_xxyy, tr_yzzzzz_xxyz, tr_yzzzzz_xxzz, tr_yzzzzz_xyyy, tr_yzzzzz_xyyz, tr_yzzzzz_xyzz, tr_yzzzzz_xzzz, tr_yzzzzz_yyyy, tr_yzzzzz_yyyz, tr_yzzzzz_yyzz, tr_yzzzzz_yzzz, tr_yzzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yzzz_xxxx[i] = 6.0 * tr_yz_xxxx[i] - 12.0 * tr_yzz_xxxxz[i] * tke_0 - 14.0 * tr_yzzz_xxxx[i] * tbe_0 - 2.0 * tr_yzzz_xxxx[i] * tke_0 + 4.0 * tr_yzzz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xxxy[i] = 6.0 * tr_yz_xxxy[i] - 12.0 * tr_yzz_xxxyz[i] * tke_0 - 14.0 * tr_yzzz_xxxy[i] * tbe_0 - 2.0 * tr_yzzz_xxxy[i] * tke_0 + 4.0 * tr_yzzz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xxxz[i] = 6.0 * tr_yz_xxxz[i] + 6.0 * tr_yzz_xxx[i] - 12.0 * tr_yzz_xxxzz[i] * tke_0 - 14.0 * tr_yzzz_xxxz[i] * tbe_0 - 6.0 * tr_yzzz_xxxz[i] * tke_0 + 4.0 * tr_yzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xxx[i] * tbe_0 + 8.0 * tr_yzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xxyy[i] = 6.0 * tr_yz_xxyy[i] - 12.0 * tr_yzz_xxyyz[i] * tke_0 - 14.0 * tr_yzzz_xxyy[i] * tbe_0 - 2.0 * tr_yzzz_xxyy[i] * tke_0 + 4.0 * tr_yzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xxyz[i] = 6.0 * tr_yz_xxyz[i] + 6.0 * tr_yzz_xxy[i] - 12.0 * tr_yzz_xxyzz[i] * tke_0 - 14.0 * tr_yzzz_xxyz[i] * tbe_0 - 6.0 * tr_yzzz_xxyz[i] * tke_0 + 4.0 * tr_yzzz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xxy[i] * tbe_0 + 8.0 * tr_yzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xxzz[i] = 6.0 * tr_yz_xxzz[i] + 12.0 * tr_yzz_xxz[i] - 12.0 * tr_yzz_xxzzz[i] * tke_0 + 2.0 * tr_yzzz_xx[i] - 14.0 * tr_yzzz_xxzz[i] * tbe_0 - 10.0 * tr_yzzz_xxzz[i] * tke_0 + 4.0 * tr_yzzz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_xxz[i] * tbe_0 + 8.0 * tr_yzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xyyy[i] = 6.0 * tr_yz_xyyy[i] - 12.0 * tr_yzz_xyyyz[i] * tke_0 - 14.0 * tr_yzzz_xyyy[i] * tbe_0 - 2.0 * tr_yzzz_xyyy[i] * tke_0 + 4.0 * tr_yzzz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xyyz[i] = 6.0 * tr_yz_xyyz[i] + 6.0 * tr_yzz_xyy[i] - 12.0 * tr_yzz_xyyzz[i] * tke_0 - 14.0 * tr_yzzz_xyyz[i] * tbe_0 - 6.0 * tr_yzzz_xyyz[i] * tke_0 + 4.0 * tr_yzzz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_xyy[i] * tbe_0 + 8.0 * tr_yzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xyzz[i] = 6.0 * tr_yz_xyzz[i] + 12.0 * tr_yzz_xyz[i] - 12.0 * tr_yzz_xyzzz[i] * tke_0 + 2.0 * tr_yzzz_xy[i] - 14.0 * tr_yzzz_xyzz[i] * tbe_0 - 10.0 * tr_yzzz_xyzz[i] * tke_0 + 4.0 * tr_yzzz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_xyz[i] * tbe_0 + 8.0 * tr_yzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xzzz[i] = 6.0 * tr_yz_xzzz[i] + 18.0 * tr_yzz_xzz[i] - 12.0 * tr_yzz_xzzzz[i] * tke_0 + 6.0 * tr_yzzz_xz[i] - 14.0 * tr_yzzz_xzzz[i] * tbe_0 - 14.0 * tr_yzzz_xzzz[i] * tke_0 + 4.0 * tr_yzzz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yzzzz_xzz[i] * tbe_0 + 8.0 * tr_yzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_yyyy[i] = 6.0 * tr_yz_yyyy[i] - 12.0 * tr_yzz_yyyyz[i] * tke_0 - 14.0 * tr_yzzz_yyyy[i] * tbe_0 - 2.0 * tr_yzzz_yyyy[i] * tke_0 + 4.0 * tr_yzzz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_yyyz[i] = 6.0 * tr_yz_yyyz[i] + 6.0 * tr_yzz_yyy[i] - 12.0 * tr_yzz_yyyzz[i] * tke_0 - 14.0 * tr_yzzz_yyyz[i] * tbe_0 - 6.0 * tr_yzzz_yyyz[i] * tke_0 + 4.0 * tr_yzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_yyy[i] * tbe_0 + 8.0 * tr_yzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_yyzz[i] = 6.0 * tr_yz_yyzz[i] + 12.0 * tr_yzz_yyz[i] - 12.0 * tr_yzz_yyzzz[i] * tke_0 + 2.0 * tr_yzzz_yy[i] - 14.0 * tr_yzzz_yyzz[i] * tbe_0 - 10.0 * tr_yzzz_yyzz[i] * tke_0 + 4.0 * tr_yzzz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_yyz[i] * tbe_0 + 8.0 * tr_yzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_yzzz[i] = 6.0 * tr_yz_yzzz[i] + 18.0 * tr_yzz_yzz[i] - 12.0 * tr_yzz_yzzzz[i] * tke_0 + 6.0 * tr_yzzz_yz[i] - 14.0 * tr_yzzz_yzzz[i] * tbe_0 - 14.0 * tr_yzzz_yzzz[i] * tke_0 + 4.0 * tr_yzzz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_yzzzz_yzz[i] * tbe_0 + 8.0 * tr_yzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_zzzz[i] = 6.0 * tr_yz_zzzz[i] + 24.0 * tr_yzz_zzz[i] - 12.0 * tr_yzz_zzzzz[i] * tke_0 + 12.0 * tr_yzzz_zz[i] - 14.0 * tr_yzzz_zzzz[i] * tbe_0 - 18.0 * tr_yzzz_zzzz[i] * tke_0 + 4.0 * tr_yzzz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_yzzzz_zzz[i] * tbe_0 + 8.0 * tr_yzzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1335-1350 components of targeted buffer : GG

    auto tr_0_0_zz_zzzz_xxxx = pbuffer.data(idx_op_geom_020_gg + 1335);

    auto tr_0_0_zz_zzzz_xxxy = pbuffer.data(idx_op_geom_020_gg + 1336);

    auto tr_0_0_zz_zzzz_xxxz = pbuffer.data(idx_op_geom_020_gg + 1337);

    auto tr_0_0_zz_zzzz_xxyy = pbuffer.data(idx_op_geom_020_gg + 1338);

    auto tr_0_0_zz_zzzz_xxyz = pbuffer.data(idx_op_geom_020_gg + 1339);

    auto tr_0_0_zz_zzzz_xxzz = pbuffer.data(idx_op_geom_020_gg + 1340);

    auto tr_0_0_zz_zzzz_xyyy = pbuffer.data(idx_op_geom_020_gg + 1341);

    auto tr_0_0_zz_zzzz_xyyz = pbuffer.data(idx_op_geom_020_gg + 1342);

    auto tr_0_0_zz_zzzz_xyzz = pbuffer.data(idx_op_geom_020_gg + 1343);

    auto tr_0_0_zz_zzzz_xzzz = pbuffer.data(idx_op_geom_020_gg + 1344);

    auto tr_0_0_zz_zzzz_yyyy = pbuffer.data(idx_op_geom_020_gg + 1345);

    auto tr_0_0_zz_zzzz_yyyz = pbuffer.data(idx_op_geom_020_gg + 1346);

    auto tr_0_0_zz_zzzz_yyzz = pbuffer.data(idx_op_geom_020_gg + 1347);

    auto tr_0_0_zz_zzzz_yzzz = pbuffer.data(idx_op_geom_020_gg + 1348);

    auto tr_0_0_zz_zzzz_zzzz = pbuffer.data(idx_op_geom_020_gg + 1349);

    #pragma omp simd aligned(tr_0_0_zz_zzzz_xxxx, tr_0_0_zz_zzzz_xxxy, tr_0_0_zz_zzzz_xxxz, tr_0_0_zz_zzzz_xxyy, tr_0_0_zz_zzzz_xxyz, tr_0_0_zz_zzzz_xxzz, tr_0_0_zz_zzzz_xyyy, tr_0_0_zz_zzzz_xyyz, tr_0_0_zz_zzzz_xyzz, tr_0_0_zz_zzzz_xzzz, tr_0_0_zz_zzzz_yyyy, tr_0_0_zz_zzzz_yyyz, tr_0_0_zz_zzzz_yyzz, tr_0_0_zz_zzzz_yzzz, tr_0_0_zz_zzzz_zzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxxxz, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, tr_zzz_zzzzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxxzz, tr_zzzz_xxxy, tr_zzzz_xxxyzz, tr_zzzz_xxxz, tr_zzzz_xxxzzz, tr_zzzz_xxyy, tr_zzzz_xxyyzz, tr_zzzz_xxyz, tr_zzzz_xxyzzz, tr_zzzz_xxzz, tr_zzzz_xxzzzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyyzz, tr_zzzz_xyyz, tr_zzzz_xyyzzz, tr_zzzz_xyzz, tr_zzzz_xyzzzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_xzzzzz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyyzz, tr_zzzz_yyyz, tr_zzzz_yyyzzz, tr_zzzz_yyzz, tr_zzzz_yyzzzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_yzzzzz, tr_zzzz_zz, tr_zzzz_zzzz, tr_zzzz_zzzzzz, tr_zzzzz_xxx, tr_zzzzz_xxxxz, tr_zzzzz_xxxyz, tr_zzzzz_xxxzz, tr_zzzzz_xxy, tr_zzzzz_xxyyz, tr_zzzzz_xxyzz, tr_zzzzz_xxz, tr_zzzzz_xxzzz, tr_zzzzz_xyy, tr_zzzzz_xyyyz, tr_zzzzz_xyyzz, tr_zzzzz_xyz, tr_zzzzz_xyzzz, tr_zzzzz_xzz, tr_zzzzz_xzzzz, tr_zzzzz_yyy, tr_zzzzz_yyyyz, tr_zzzzz_yyyzz, tr_zzzzz_yyz, tr_zzzzz_yyzzz, tr_zzzzz_yzz, tr_zzzzz_yzzzz, tr_zzzzz_zzz, tr_zzzzz_zzzzz, tr_zzzzzz_xxxx, tr_zzzzzz_xxxy, tr_zzzzzz_xxxz, tr_zzzzzz_xxyy, tr_zzzzzz_xxyz, tr_zzzzzz_xxzz, tr_zzzzzz_xyyy, tr_zzzzzz_xyyz, tr_zzzzzz_xyzz, tr_zzzzzz_xzzz, tr_zzzzzz_yyyy, tr_zzzzzz_yyyz, tr_zzzzzz_yyzz, tr_zzzzzz_yzzz, tr_zzzzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_zzzz_xxxx[i] = 12.0 * tr_zz_xxxx[i] - 16.0 * tr_zzz_xxxxz[i] * tke_0 - 18.0 * tr_zzzz_xxxx[i] * tbe_0 - 2.0 * tr_zzzz_xxxx[i] * tke_0 + 4.0 * tr_zzzz_xxxxzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_xxxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxxx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xxxy[i] = 12.0 * tr_zz_xxxy[i] - 16.0 * tr_zzz_xxxyz[i] * tke_0 - 18.0 * tr_zzzz_xxxy[i] * tbe_0 - 2.0 * tr_zzzz_xxxy[i] * tke_0 + 4.0 * tr_zzzz_xxxyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_xxxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxxy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xxxz[i] = 12.0 * tr_zz_xxxz[i] + 8.0 * tr_zzz_xxx[i] - 16.0 * tr_zzz_xxxzz[i] * tke_0 - 18.0 * tr_zzzz_xxxz[i] * tbe_0 - 6.0 * tr_zzzz_xxxz[i] * tke_0 + 4.0 * tr_zzzz_xxxzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xxx[i] * tbe_0 + 8.0 * tr_zzzzz_xxxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxxz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xxyy[i] = 12.0 * tr_zz_xxyy[i] - 16.0 * tr_zzz_xxyyz[i] * tke_0 - 18.0 * tr_zzzz_xxyy[i] * tbe_0 - 2.0 * tr_zzzz_xxyy[i] * tke_0 + 4.0 * tr_zzzz_xxyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_xxyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xxyz[i] = 12.0 * tr_zz_xxyz[i] + 8.0 * tr_zzz_xxy[i] - 16.0 * tr_zzz_xxyzz[i] * tke_0 - 18.0 * tr_zzzz_xxyz[i] * tbe_0 - 6.0 * tr_zzzz_xxyz[i] * tke_0 + 4.0 * tr_zzzz_xxyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xxy[i] * tbe_0 + 8.0 * tr_zzzzz_xxyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xxzz[i] = 12.0 * tr_zz_xxzz[i] + 16.0 * tr_zzz_xxz[i] - 16.0 * tr_zzz_xxzzz[i] * tke_0 + 2.0 * tr_zzzz_xx[i] - 18.0 * tr_zzzz_xxzz[i] * tbe_0 - 10.0 * tr_zzzz_xxzz[i] * tke_0 + 4.0 * tr_zzzz_xxzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzzzz_xxz[i] * tbe_0 + 8.0 * tr_zzzzz_xxzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xyyy[i] = 12.0 * tr_zz_xyyy[i] - 16.0 * tr_zzz_xyyyz[i] * tke_0 - 18.0 * tr_zzzz_xyyy[i] * tbe_0 - 2.0 * tr_zzzz_xyyy[i] * tke_0 + 4.0 * tr_zzzz_xyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_xyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xyyz[i] = 12.0 * tr_zz_xyyz[i] + 8.0 * tr_zzz_xyy[i] - 16.0 * tr_zzz_xyyzz[i] * tke_0 - 18.0 * tr_zzzz_xyyz[i] * tbe_0 - 6.0 * tr_zzzz_xyyz[i] * tke_0 + 4.0 * tr_zzzz_xyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_xyy[i] * tbe_0 + 8.0 * tr_zzzzz_xyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xyzz[i] = 12.0 * tr_zz_xyzz[i] + 16.0 * tr_zzz_xyz[i] - 16.0 * tr_zzz_xyzzz[i] * tke_0 + 2.0 * tr_zzzz_xy[i] - 18.0 * tr_zzzz_xyzz[i] * tbe_0 - 10.0 * tr_zzzz_xyzz[i] * tke_0 + 4.0 * tr_zzzz_xyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzzzz_xyz[i] * tbe_0 + 8.0 * tr_zzzzz_xyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xzzz[i] = 12.0 * tr_zz_xzzz[i] + 24.0 * tr_zzz_xzz[i] - 16.0 * tr_zzz_xzzzz[i] * tke_0 + 6.0 * tr_zzzz_xz[i] - 18.0 * tr_zzzz_xzzz[i] * tbe_0 - 14.0 * tr_zzzz_xzzz[i] * tke_0 + 4.0 * tr_zzzz_xzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zzzzz_xzz[i] * tbe_0 + 8.0 * tr_zzzzz_xzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_yyyy[i] = 12.0 * tr_zz_yyyy[i] - 16.0 * tr_zzz_yyyyz[i] * tke_0 - 18.0 * tr_zzzz_yyyy[i] * tbe_0 - 2.0 * tr_zzzz_yyyy[i] * tke_0 + 4.0 * tr_zzzz_yyyyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_yyyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yyyy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_yyyz[i] = 12.0 * tr_zz_yyyz[i] + 8.0 * tr_zzz_yyy[i] - 16.0 * tr_zzz_yyyzz[i] * tke_0 - 18.0 * tr_zzzz_yyyz[i] * tbe_0 - 6.0 * tr_zzzz_yyyz[i] * tke_0 + 4.0 * tr_zzzz_yyyzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_yyy[i] * tbe_0 + 8.0 * tr_zzzzz_yyyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yyyz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_yyzz[i] = 12.0 * tr_zz_yyzz[i] + 16.0 * tr_zzz_yyz[i] - 16.0 * tr_zzz_yyzzz[i] * tke_0 + 2.0 * tr_zzzz_yy[i] - 18.0 * tr_zzzz_yyzz[i] * tbe_0 - 10.0 * tr_zzzz_yyzz[i] * tke_0 + 4.0 * tr_zzzz_yyzzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzzzz_yyz[i] * tbe_0 + 8.0 * tr_zzzzz_yyzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yyzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_yzzz[i] = 12.0 * tr_zz_yzzz[i] + 24.0 * tr_zzz_yzz[i] - 16.0 * tr_zzz_yzzzz[i] * tke_0 + 6.0 * tr_zzzz_yz[i] - 18.0 * tr_zzzz_yzzz[i] * tbe_0 - 14.0 * tr_zzzz_yzzz[i] * tke_0 + 4.0 * tr_zzzz_yzzzzz[i] * tke_0 * tke_0 - 12.0 * tr_zzzzz_yzz[i] * tbe_0 + 8.0 * tr_zzzzz_yzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yzzz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_zzzz[i] = 12.0 * tr_zz_zzzz[i] + 32.0 * tr_zzz_zzz[i] - 16.0 * tr_zzz_zzzzz[i] * tke_0 + 12.0 * tr_zzzz_zz[i] - 18.0 * tr_zzzz_zzzz[i] * tbe_0 - 18.0 * tr_zzzz_zzzz[i] * tke_0 + 4.0 * tr_zzzz_zzzzzz[i] * tke_0 * tke_0 - 16.0 * tr_zzzzz_zzz[i] * tbe_0 + 8.0 * tr_zzzzz_zzzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_zzzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

