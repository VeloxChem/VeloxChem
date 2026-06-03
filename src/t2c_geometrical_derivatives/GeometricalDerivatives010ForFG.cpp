#include "GeometricalDerivatives010ForFG.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_fg(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_fg,
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

    // Set up 0-15 components of targeted buffer : FG

    auto tr_0_0_x_xxx_xxxx = pbuffer.data(idx_op_geom_010_fg);

    auto tr_0_0_x_xxx_xxxy = pbuffer.data(idx_op_geom_010_fg + 1);

    auto tr_0_0_x_xxx_xxxz = pbuffer.data(idx_op_geom_010_fg + 2);

    auto tr_0_0_x_xxx_xxyy = pbuffer.data(idx_op_geom_010_fg + 3);

    auto tr_0_0_x_xxx_xxyz = pbuffer.data(idx_op_geom_010_fg + 4);

    auto tr_0_0_x_xxx_xxzz = pbuffer.data(idx_op_geom_010_fg + 5);

    auto tr_0_0_x_xxx_xyyy = pbuffer.data(idx_op_geom_010_fg + 6);

    auto tr_0_0_x_xxx_xyyz = pbuffer.data(idx_op_geom_010_fg + 7);

    auto tr_0_0_x_xxx_xyzz = pbuffer.data(idx_op_geom_010_fg + 8);

    auto tr_0_0_x_xxx_xzzz = pbuffer.data(idx_op_geom_010_fg + 9);

    auto tr_0_0_x_xxx_yyyy = pbuffer.data(idx_op_geom_010_fg + 10);

    auto tr_0_0_x_xxx_yyyz = pbuffer.data(idx_op_geom_010_fg + 11);

    auto tr_0_0_x_xxx_yyzz = pbuffer.data(idx_op_geom_010_fg + 12);

    auto tr_0_0_x_xxx_yzzz = pbuffer.data(idx_op_geom_010_fg + 13);

    auto tr_0_0_x_xxx_zzzz = pbuffer.data(idx_op_geom_010_fg + 14);

    #pragma omp simd aligned(tr_0_0_x_xxx_xxxx, tr_0_0_x_xxx_xxxy, tr_0_0_x_xxx_xxxz, tr_0_0_x_xxx_xxyy, tr_0_0_x_xxx_xxyz, tr_0_0_x_xxx_xxzz, tr_0_0_x_xxx_xyyy, tr_0_0_x_xxx_xyyz, tr_0_0_x_xxx_xyzz, tr_0_0_x_xxx_xzzz, tr_0_0_x_xxx_yyyy, tr_0_0_x_xxx_yyyz, tr_0_0_x_xxx_yyzz, tr_0_0_x_xxx_yzzz, tr_0_0_x_xxx_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxx_xxx, tr_xxx_xxxxx, tr_xxx_xxxxy, tr_xxx_xxxxz, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xzzz, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yzzz, tr_xxxx_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxx_xxxx[i] = 2.0 * tr_xxxx_xxxx[i] * tbe_0 + 2.0 * tr_xxx_xxxxx[i] * tke_0 - 3.0 * tr_xx_xxxx[i] - 4.0 * tr_xxx_xxx[i];

        tr_0_0_x_xxx_xxxy[i] = 2.0 * tr_xxxx_xxxy[i] * tbe_0 + 2.0 * tr_xxx_xxxxy[i] * tke_0 - 3.0 * tr_xx_xxxy[i] - 3.0 * tr_xxx_xxy[i];

        tr_0_0_x_xxx_xxxz[i] = 2.0 * tr_xxxx_xxxz[i] * tbe_0 + 2.0 * tr_xxx_xxxxz[i] * tke_0 - 3.0 * tr_xx_xxxz[i] - 3.0 * tr_xxx_xxz[i];

        tr_0_0_x_xxx_xxyy[i] = 2.0 * tr_xxxx_xxyy[i] * tbe_0 + 2.0 * tr_xxx_xxxyy[i] * tke_0 - 3.0 * tr_xx_xxyy[i] - 2.0 * tr_xxx_xyy[i];

        tr_0_0_x_xxx_xxyz[i] = 2.0 * tr_xxxx_xxyz[i] * tbe_0 + 2.0 * tr_xxx_xxxyz[i] * tke_0 - 3.0 * tr_xx_xxyz[i] - 2.0 * tr_xxx_xyz[i];

        tr_0_0_x_xxx_xxzz[i] = 2.0 * tr_xxxx_xxzz[i] * tbe_0 + 2.0 * tr_xxx_xxxzz[i] * tke_0 - 3.0 * tr_xx_xxzz[i] - 2.0 * tr_xxx_xzz[i];

        tr_0_0_x_xxx_xyyy[i] = 2.0 * tr_xxxx_xyyy[i] * tbe_0 + 2.0 * tr_xxx_xxyyy[i] * tke_0 - 3.0 * tr_xx_xyyy[i] - tr_xxx_yyy[i];

        tr_0_0_x_xxx_xyyz[i] = 2.0 * tr_xxxx_xyyz[i] * tbe_0 + 2.0 * tr_xxx_xxyyz[i] * tke_0 - 3.0 * tr_xx_xyyz[i] - tr_xxx_yyz[i];

        tr_0_0_x_xxx_xyzz[i] = 2.0 * tr_xxxx_xyzz[i] * tbe_0 + 2.0 * tr_xxx_xxyzz[i] * tke_0 - 3.0 * tr_xx_xyzz[i] - tr_xxx_yzz[i];

        tr_0_0_x_xxx_xzzz[i] = 2.0 * tr_xxxx_xzzz[i] * tbe_0 + 2.0 * tr_xxx_xxzzz[i] * tke_0 - 3.0 * tr_xx_xzzz[i] - tr_xxx_zzz[i];

        tr_0_0_x_xxx_yyyy[i] = 2.0 * tr_xxxx_yyyy[i] * tbe_0 + 2.0 * tr_xxx_xyyyy[i] * tke_0 - 3.0 * tr_xx_yyyy[i];

        tr_0_0_x_xxx_yyyz[i] = 2.0 * tr_xxxx_yyyz[i] * tbe_0 + 2.0 * tr_xxx_xyyyz[i] * tke_0 - 3.0 * tr_xx_yyyz[i];

        tr_0_0_x_xxx_yyzz[i] = 2.0 * tr_xxxx_yyzz[i] * tbe_0 + 2.0 * tr_xxx_xyyzz[i] * tke_0 - 3.0 * tr_xx_yyzz[i];

        tr_0_0_x_xxx_yzzz[i] = 2.0 * tr_xxxx_yzzz[i] * tbe_0 + 2.0 * tr_xxx_xyzzz[i] * tke_0 - 3.0 * tr_xx_yzzz[i];

        tr_0_0_x_xxx_zzzz[i] = 2.0 * tr_xxxx_zzzz[i] * tbe_0 + 2.0 * tr_xxx_xzzzz[i] * tke_0 - 3.0 * tr_xx_zzzz[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto tr_0_0_x_xxy_xxxx = pbuffer.data(idx_op_geom_010_fg + 15);

    auto tr_0_0_x_xxy_xxxy = pbuffer.data(idx_op_geom_010_fg + 16);

    auto tr_0_0_x_xxy_xxxz = pbuffer.data(idx_op_geom_010_fg + 17);

    auto tr_0_0_x_xxy_xxyy = pbuffer.data(idx_op_geom_010_fg + 18);

    auto tr_0_0_x_xxy_xxyz = pbuffer.data(idx_op_geom_010_fg + 19);

    auto tr_0_0_x_xxy_xxzz = pbuffer.data(idx_op_geom_010_fg + 20);

    auto tr_0_0_x_xxy_xyyy = pbuffer.data(idx_op_geom_010_fg + 21);

    auto tr_0_0_x_xxy_xyyz = pbuffer.data(idx_op_geom_010_fg + 22);

    auto tr_0_0_x_xxy_xyzz = pbuffer.data(idx_op_geom_010_fg + 23);

    auto tr_0_0_x_xxy_xzzz = pbuffer.data(idx_op_geom_010_fg + 24);

    auto tr_0_0_x_xxy_yyyy = pbuffer.data(idx_op_geom_010_fg + 25);

    auto tr_0_0_x_xxy_yyyz = pbuffer.data(idx_op_geom_010_fg + 26);

    auto tr_0_0_x_xxy_yyzz = pbuffer.data(idx_op_geom_010_fg + 27);

    auto tr_0_0_x_xxy_yzzz = pbuffer.data(idx_op_geom_010_fg + 28);

    auto tr_0_0_x_xxy_zzzz = pbuffer.data(idx_op_geom_010_fg + 29);

    #pragma omp simd aligned(tr_0_0_x_xxy_xxxx, tr_0_0_x_xxy_xxxy, tr_0_0_x_xxy_xxxz, tr_0_0_x_xxy_xxyy, tr_0_0_x_xxy_xxyz, tr_0_0_x_xxy_xxzz, tr_0_0_x_xxy_xyyy, tr_0_0_x_xxy_xyyz, tr_0_0_x_xxy_xyzz, tr_0_0_x_xxy_xzzz, tr_0_0_x_xxy_yyyy, tr_0_0_x_xxy_yyyz, tr_0_0_x_xxy_yyzz, tr_0_0_x_xxy_yzzz, tr_0_0_x_xxy_zzzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, tr_xxy_xxx, tr_xxy_xxxxx, tr_xxy_xxxxy, tr_xxy_xxxxz, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxy_xxxx[i] = 2.0 * tr_xxxy_xxxx[i] * tbe_0 + 2.0 * tr_xxy_xxxxx[i] * tke_0 - 2.0 * tr_xy_xxxx[i] - 4.0 * tr_xxy_xxx[i];

        tr_0_0_x_xxy_xxxy[i] = 2.0 * tr_xxxy_xxxy[i] * tbe_0 + 2.0 * tr_xxy_xxxxy[i] * tke_0 - 2.0 * tr_xy_xxxy[i] - 3.0 * tr_xxy_xxy[i];

        tr_0_0_x_xxy_xxxz[i] = 2.0 * tr_xxxy_xxxz[i] * tbe_0 + 2.0 * tr_xxy_xxxxz[i] * tke_0 - 2.0 * tr_xy_xxxz[i] - 3.0 * tr_xxy_xxz[i];

        tr_0_0_x_xxy_xxyy[i] = 2.0 * tr_xxxy_xxyy[i] * tbe_0 + 2.0 * tr_xxy_xxxyy[i] * tke_0 - 2.0 * tr_xy_xxyy[i] - 2.0 * tr_xxy_xyy[i];

        tr_0_0_x_xxy_xxyz[i] = 2.0 * tr_xxxy_xxyz[i] * tbe_0 + 2.0 * tr_xxy_xxxyz[i] * tke_0 - 2.0 * tr_xy_xxyz[i] - 2.0 * tr_xxy_xyz[i];

        tr_0_0_x_xxy_xxzz[i] = 2.0 * tr_xxxy_xxzz[i] * tbe_0 + 2.0 * tr_xxy_xxxzz[i] * tke_0 - 2.0 * tr_xy_xxzz[i] - 2.0 * tr_xxy_xzz[i];

        tr_0_0_x_xxy_xyyy[i] = 2.0 * tr_xxxy_xyyy[i] * tbe_0 + 2.0 * tr_xxy_xxyyy[i] * tke_0 - 2.0 * tr_xy_xyyy[i] - tr_xxy_yyy[i];

        tr_0_0_x_xxy_xyyz[i] = 2.0 * tr_xxxy_xyyz[i] * tbe_0 + 2.0 * tr_xxy_xxyyz[i] * tke_0 - 2.0 * tr_xy_xyyz[i] - tr_xxy_yyz[i];

        tr_0_0_x_xxy_xyzz[i] = 2.0 * tr_xxxy_xyzz[i] * tbe_0 + 2.0 * tr_xxy_xxyzz[i] * tke_0 - 2.0 * tr_xy_xyzz[i] - tr_xxy_yzz[i];

        tr_0_0_x_xxy_xzzz[i] = 2.0 * tr_xxxy_xzzz[i] * tbe_0 + 2.0 * tr_xxy_xxzzz[i] * tke_0 - 2.0 * tr_xy_xzzz[i] - tr_xxy_zzz[i];

        tr_0_0_x_xxy_yyyy[i] = 2.0 * tr_xxxy_yyyy[i] * tbe_0 + 2.0 * tr_xxy_xyyyy[i] * tke_0 - 2.0 * tr_xy_yyyy[i];

        tr_0_0_x_xxy_yyyz[i] = 2.0 * tr_xxxy_yyyz[i] * tbe_0 + 2.0 * tr_xxy_xyyyz[i] * tke_0 - 2.0 * tr_xy_yyyz[i];

        tr_0_0_x_xxy_yyzz[i] = 2.0 * tr_xxxy_yyzz[i] * tbe_0 + 2.0 * tr_xxy_xyyzz[i] * tke_0 - 2.0 * tr_xy_yyzz[i];

        tr_0_0_x_xxy_yzzz[i] = 2.0 * tr_xxxy_yzzz[i] * tbe_0 + 2.0 * tr_xxy_xyzzz[i] * tke_0 - 2.0 * tr_xy_yzzz[i];

        tr_0_0_x_xxy_zzzz[i] = 2.0 * tr_xxxy_zzzz[i] * tbe_0 + 2.0 * tr_xxy_xzzzz[i] * tke_0 - 2.0 * tr_xy_zzzz[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto tr_0_0_x_xxz_xxxx = pbuffer.data(idx_op_geom_010_fg + 30);

    auto tr_0_0_x_xxz_xxxy = pbuffer.data(idx_op_geom_010_fg + 31);

    auto tr_0_0_x_xxz_xxxz = pbuffer.data(idx_op_geom_010_fg + 32);

    auto tr_0_0_x_xxz_xxyy = pbuffer.data(idx_op_geom_010_fg + 33);

    auto tr_0_0_x_xxz_xxyz = pbuffer.data(idx_op_geom_010_fg + 34);

    auto tr_0_0_x_xxz_xxzz = pbuffer.data(idx_op_geom_010_fg + 35);

    auto tr_0_0_x_xxz_xyyy = pbuffer.data(idx_op_geom_010_fg + 36);

    auto tr_0_0_x_xxz_xyyz = pbuffer.data(idx_op_geom_010_fg + 37);

    auto tr_0_0_x_xxz_xyzz = pbuffer.data(idx_op_geom_010_fg + 38);

    auto tr_0_0_x_xxz_xzzz = pbuffer.data(idx_op_geom_010_fg + 39);

    auto tr_0_0_x_xxz_yyyy = pbuffer.data(idx_op_geom_010_fg + 40);

    auto tr_0_0_x_xxz_yyyz = pbuffer.data(idx_op_geom_010_fg + 41);

    auto tr_0_0_x_xxz_yyzz = pbuffer.data(idx_op_geom_010_fg + 42);

    auto tr_0_0_x_xxz_yzzz = pbuffer.data(idx_op_geom_010_fg + 43);

    auto tr_0_0_x_xxz_zzzz = pbuffer.data(idx_op_geom_010_fg + 44);

    #pragma omp simd aligned(tr_0_0_x_xxz_xxxx, tr_0_0_x_xxz_xxxy, tr_0_0_x_xxz_xxxz, tr_0_0_x_xxz_xxyy, tr_0_0_x_xxz_xxyz, tr_0_0_x_xxz_xxzz, tr_0_0_x_xxz_xyyy, tr_0_0_x_xxz_xyyz, tr_0_0_x_xxz_xyzz, tr_0_0_x_xxz_xzzz, tr_0_0_x_xxz_yyyy, tr_0_0_x_xxz_yyyz, tr_0_0_x_xxz_yyzz, tr_0_0_x_xxz_yzzz, tr_0_0_x_xxz_zzzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, tr_xxz_xxx, tr_xxz_xxxxx, tr_xxz_xxxxy, tr_xxz_xxxxz, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxz_xxxx[i] = 2.0 * tr_xxxz_xxxx[i] * tbe_0 + 2.0 * tr_xxz_xxxxx[i] * tke_0 - 2.0 * tr_xz_xxxx[i] - 4.0 * tr_xxz_xxx[i];

        tr_0_0_x_xxz_xxxy[i] = 2.0 * tr_xxxz_xxxy[i] * tbe_0 + 2.0 * tr_xxz_xxxxy[i] * tke_0 - 2.0 * tr_xz_xxxy[i] - 3.0 * tr_xxz_xxy[i];

        tr_0_0_x_xxz_xxxz[i] = 2.0 * tr_xxxz_xxxz[i] * tbe_0 + 2.0 * tr_xxz_xxxxz[i] * tke_0 - 2.0 * tr_xz_xxxz[i] - 3.0 * tr_xxz_xxz[i];

        tr_0_0_x_xxz_xxyy[i] = 2.0 * tr_xxxz_xxyy[i] * tbe_0 + 2.0 * tr_xxz_xxxyy[i] * tke_0 - 2.0 * tr_xz_xxyy[i] - 2.0 * tr_xxz_xyy[i];

        tr_0_0_x_xxz_xxyz[i] = 2.0 * tr_xxxz_xxyz[i] * tbe_0 + 2.0 * tr_xxz_xxxyz[i] * tke_0 - 2.0 * tr_xz_xxyz[i] - 2.0 * tr_xxz_xyz[i];

        tr_0_0_x_xxz_xxzz[i] = 2.0 * tr_xxxz_xxzz[i] * tbe_0 + 2.0 * tr_xxz_xxxzz[i] * tke_0 - 2.0 * tr_xz_xxzz[i] - 2.0 * tr_xxz_xzz[i];

        tr_0_0_x_xxz_xyyy[i] = 2.0 * tr_xxxz_xyyy[i] * tbe_0 + 2.0 * tr_xxz_xxyyy[i] * tke_0 - 2.0 * tr_xz_xyyy[i] - tr_xxz_yyy[i];

        tr_0_0_x_xxz_xyyz[i] = 2.0 * tr_xxxz_xyyz[i] * tbe_0 + 2.0 * tr_xxz_xxyyz[i] * tke_0 - 2.0 * tr_xz_xyyz[i] - tr_xxz_yyz[i];

        tr_0_0_x_xxz_xyzz[i] = 2.0 * tr_xxxz_xyzz[i] * tbe_0 + 2.0 * tr_xxz_xxyzz[i] * tke_0 - 2.0 * tr_xz_xyzz[i] - tr_xxz_yzz[i];

        tr_0_0_x_xxz_xzzz[i] = 2.0 * tr_xxxz_xzzz[i] * tbe_0 + 2.0 * tr_xxz_xxzzz[i] * tke_0 - 2.0 * tr_xz_xzzz[i] - tr_xxz_zzz[i];

        tr_0_0_x_xxz_yyyy[i] = 2.0 * tr_xxxz_yyyy[i] * tbe_0 + 2.0 * tr_xxz_xyyyy[i] * tke_0 - 2.0 * tr_xz_yyyy[i];

        tr_0_0_x_xxz_yyyz[i] = 2.0 * tr_xxxz_yyyz[i] * tbe_0 + 2.0 * tr_xxz_xyyyz[i] * tke_0 - 2.0 * tr_xz_yyyz[i];

        tr_0_0_x_xxz_yyzz[i] = 2.0 * tr_xxxz_yyzz[i] * tbe_0 + 2.0 * tr_xxz_xyyzz[i] * tke_0 - 2.0 * tr_xz_yyzz[i];

        tr_0_0_x_xxz_yzzz[i] = 2.0 * tr_xxxz_yzzz[i] * tbe_0 + 2.0 * tr_xxz_xyzzz[i] * tke_0 - 2.0 * tr_xz_yzzz[i];

        tr_0_0_x_xxz_zzzz[i] = 2.0 * tr_xxxz_zzzz[i] * tbe_0 + 2.0 * tr_xxz_xzzzz[i] * tke_0 - 2.0 * tr_xz_zzzz[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto tr_0_0_x_xyy_xxxx = pbuffer.data(idx_op_geom_010_fg + 45);

    auto tr_0_0_x_xyy_xxxy = pbuffer.data(idx_op_geom_010_fg + 46);

    auto tr_0_0_x_xyy_xxxz = pbuffer.data(idx_op_geom_010_fg + 47);

    auto tr_0_0_x_xyy_xxyy = pbuffer.data(idx_op_geom_010_fg + 48);

    auto tr_0_0_x_xyy_xxyz = pbuffer.data(idx_op_geom_010_fg + 49);

    auto tr_0_0_x_xyy_xxzz = pbuffer.data(idx_op_geom_010_fg + 50);

    auto tr_0_0_x_xyy_xyyy = pbuffer.data(idx_op_geom_010_fg + 51);

    auto tr_0_0_x_xyy_xyyz = pbuffer.data(idx_op_geom_010_fg + 52);

    auto tr_0_0_x_xyy_xyzz = pbuffer.data(idx_op_geom_010_fg + 53);

    auto tr_0_0_x_xyy_xzzz = pbuffer.data(idx_op_geom_010_fg + 54);

    auto tr_0_0_x_xyy_yyyy = pbuffer.data(idx_op_geom_010_fg + 55);

    auto tr_0_0_x_xyy_yyyz = pbuffer.data(idx_op_geom_010_fg + 56);

    auto tr_0_0_x_xyy_yyzz = pbuffer.data(idx_op_geom_010_fg + 57);

    auto tr_0_0_x_xyy_yzzz = pbuffer.data(idx_op_geom_010_fg + 58);

    auto tr_0_0_x_xyy_zzzz = pbuffer.data(idx_op_geom_010_fg + 59);

    #pragma omp simd aligned(tr_0_0_x_xyy_xxxx, tr_0_0_x_xyy_xxxy, tr_0_0_x_xyy_xxxz, tr_0_0_x_xyy_xxyy, tr_0_0_x_xyy_xxyz, tr_0_0_x_xyy_xxzz, tr_0_0_x_xyy_xyyy, tr_0_0_x_xyy_xyyz, tr_0_0_x_xyy_xyzz, tr_0_0_x_xyy_xzzz, tr_0_0_x_xyy_yyyy, tr_0_0_x_xyy_yyyz, tr_0_0_x_xyy_yyzz, tr_0_0_x_xyy_yzzz, tr_0_0_x_xyy_zzzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxx, tr_xyy_xxxxy, tr_xyy_xxxxz, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyy_xxxx[i] = 2.0 * tr_xxyy_xxxx[i] * tbe_0 + 2.0 * tr_xyy_xxxxx[i] * tke_0 - tr_yy_xxxx[i] - 4.0 * tr_xyy_xxx[i];

        tr_0_0_x_xyy_xxxy[i] = 2.0 * tr_xxyy_xxxy[i] * tbe_0 + 2.0 * tr_xyy_xxxxy[i] * tke_0 - tr_yy_xxxy[i] - 3.0 * tr_xyy_xxy[i];

        tr_0_0_x_xyy_xxxz[i] = 2.0 * tr_xxyy_xxxz[i] * tbe_0 + 2.0 * tr_xyy_xxxxz[i] * tke_0 - tr_yy_xxxz[i] - 3.0 * tr_xyy_xxz[i];

        tr_0_0_x_xyy_xxyy[i] = 2.0 * tr_xxyy_xxyy[i] * tbe_0 + 2.0 * tr_xyy_xxxyy[i] * tke_0 - tr_yy_xxyy[i] - 2.0 * tr_xyy_xyy[i];

        tr_0_0_x_xyy_xxyz[i] = 2.0 * tr_xxyy_xxyz[i] * tbe_0 + 2.0 * tr_xyy_xxxyz[i] * tke_0 - tr_yy_xxyz[i] - 2.0 * tr_xyy_xyz[i];

        tr_0_0_x_xyy_xxzz[i] = 2.0 * tr_xxyy_xxzz[i] * tbe_0 + 2.0 * tr_xyy_xxxzz[i] * tke_0 - tr_yy_xxzz[i] - 2.0 * tr_xyy_xzz[i];

        tr_0_0_x_xyy_xyyy[i] = 2.0 * tr_xxyy_xyyy[i] * tbe_0 + 2.0 * tr_xyy_xxyyy[i] * tke_0 - tr_yy_xyyy[i] - tr_xyy_yyy[i];

        tr_0_0_x_xyy_xyyz[i] = 2.0 * tr_xxyy_xyyz[i] * tbe_0 + 2.0 * tr_xyy_xxyyz[i] * tke_0 - tr_yy_xyyz[i] - tr_xyy_yyz[i];

        tr_0_0_x_xyy_xyzz[i] = 2.0 * tr_xxyy_xyzz[i] * tbe_0 + 2.0 * tr_xyy_xxyzz[i] * tke_0 - tr_yy_xyzz[i] - tr_xyy_yzz[i];

        tr_0_0_x_xyy_xzzz[i] = 2.0 * tr_xxyy_xzzz[i] * tbe_0 + 2.0 * tr_xyy_xxzzz[i] * tke_0 - tr_yy_xzzz[i] - tr_xyy_zzz[i];

        tr_0_0_x_xyy_yyyy[i] = 2.0 * tr_xxyy_yyyy[i] * tbe_0 + 2.0 * tr_xyy_xyyyy[i] * tke_0 - tr_yy_yyyy[i];

        tr_0_0_x_xyy_yyyz[i] = 2.0 * tr_xxyy_yyyz[i] * tbe_0 + 2.0 * tr_xyy_xyyyz[i] * tke_0 - tr_yy_yyyz[i];

        tr_0_0_x_xyy_yyzz[i] = 2.0 * tr_xxyy_yyzz[i] * tbe_0 + 2.0 * tr_xyy_xyyzz[i] * tke_0 - tr_yy_yyzz[i];

        tr_0_0_x_xyy_yzzz[i] = 2.0 * tr_xxyy_yzzz[i] * tbe_0 + 2.0 * tr_xyy_xyzzz[i] * tke_0 - tr_yy_yzzz[i];

        tr_0_0_x_xyy_zzzz[i] = 2.0 * tr_xxyy_zzzz[i] * tbe_0 + 2.0 * tr_xyy_xzzzz[i] * tke_0 - tr_yy_zzzz[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto tr_0_0_x_xyz_xxxx = pbuffer.data(idx_op_geom_010_fg + 60);

    auto tr_0_0_x_xyz_xxxy = pbuffer.data(idx_op_geom_010_fg + 61);

    auto tr_0_0_x_xyz_xxxz = pbuffer.data(idx_op_geom_010_fg + 62);

    auto tr_0_0_x_xyz_xxyy = pbuffer.data(idx_op_geom_010_fg + 63);

    auto tr_0_0_x_xyz_xxyz = pbuffer.data(idx_op_geom_010_fg + 64);

    auto tr_0_0_x_xyz_xxzz = pbuffer.data(idx_op_geom_010_fg + 65);

    auto tr_0_0_x_xyz_xyyy = pbuffer.data(idx_op_geom_010_fg + 66);

    auto tr_0_0_x_xyz_xyyz = pbuffer.data(idx_op_geom_010_fg + 67);

    auto tr_0_0_x_xyz_xyzz = pbuffer.data(idx_op_geom_010_fg + 68);

    auto tr_0_0_x_xyz_xzzz = pbuffer.data(idx_op_geom_010_fg + 69);

    auto tr_0_0_x_xyz_yyyy = pbuffer.data(idx_op_geom_010_fg + 70);

    auto tr_0_0_x_xyz_yyyz = pbuffer.data(idx_op_geom_010_fg + 71);

    auto tr_0_0_x_xyz_yyzz = pbuffer.data(idx_op_geom_010_fg + 72);

    auto tr_0_0_x_xyz_yzzz = pbuffer.data(idx_op_geom_010_fg + 73);

    auto tr_0_0_x_xyz_zzzz = pbuffer.data(idx_op_geom_010_fg + 74);

    #pragma omp simd aligned(tr_0_0_x_xyz_xxxx, tr_0_0_x_xyz_xxxy, tr_0_0_x_xyz_xxxz, tr_0_0_x_xyz_xxyy, tr_0_0_x_xyz_xxyz, tr_0_0_x_xyz_xxzz, tr_0_0_x_xyz_xyyy, tr_0_0_x_xyz_xyyz, tr_0_0_x_xyz_xyzz, tr_0_0_x_xyz_xzzz, tr_0_0_x_xyz_yyyy, tr_0_0_x_xyz_yyyz, tr_0_0_x_xyz_yyzz, tr_0_0_x_xyz_yzzz, tr_0_0_x_xyz_zzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxx, tr_xyz_xxxxy, tr_xyz_xxxxz, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyz_xxxx[i] = 2.0 * tr_xxyz_xxxx[i] * tbe_0 + 2.0 * tr_xyz_xxxxx[i] * tke_0 - tr_yz_xxxx[i] - 4.0 * tr_xyz_xxx[i];

        tr_0_0_x_xyz_xxxy[i] = 2.0 * tr_xxyz_xxxy[i] * tbe_0 + 2.0 * tr_xyz_xxxxy[i] * tke_0 - tr_yz_xxxy[i] - 3.0 * tr_xyz_xxy[i];

        tr_0_0_x_xyz_xxxz[i] = 2.0 * tr_xxyz_xxxz[i] * tbe_0 + 2.0 * tr_xyz_xxxxz[i] * tke_0 - tr_yz_xxxz[i] - 3.0 * tr_xyz_xxz[i];

        tr_0_0_x_xyz_xxyy[i] = 2.0 * tr_xxyz_xxyy[i] * tbe_0 + 2.0 * tr_xyz_xxxyy[i] * tke_0 - tr_yz_xxyy[i] - 2.0 * tr_xyz_xyy[i];

        tr_0_0_x_xyz_xxyz[i] = 2.0 * tr_xxyz_xxyz[i] * tbe_0 + 2.0 * tr_xyz_xxxyz[i] * tke_0 - tr_yz_xxyz[i] - 2.0 * tr_xyz_xyz[i];

        tr_0_0_x_xyz_xxzz[i] = 2.0 * tr_xxyz_xxzz[i] * tbe_0 + 2.0 * tr_xyz_xxxzz[i] * tke_0 - tr_yz_xxzz[i] - 2.0 * tr_xyz_xzz[i];

        tr_0_0_x_xyz_xyyy[i] = 2.0 * tr_xxyz_xyyy[i] * tbe_0 + 2.0 * tr_xyz_xxyyy[i] * tke_0 - tr_yz_xyyy[i] - tr_xyz_yyy[i];

        tr_0_0_x_xyz_xyyz[i] = 2.0 * tr_xxyz_xyyz[i] * tbe_0 + 2.0 * tr_xyz_xxyyz[i] * tke_0 - tr_yz_xyyz[i] - tr_xyz_yyz[i];

        tr_0_0_x_xyz_xyzz[i] = 2.0 * tr_xxyz_xyzz[i] * tbe_0 + 2.0 * tr_xyz_xxyzz[i] * tke_0 - tr_yz_xyzz[i] - tr_xyz_yzz[i];

        tr_0_0_x_xyz_xzzz[i] = 2.0 * tr_xxyz_xzzz[i] * tbe_0 + 2.0 * tr_xyz_xxzzz[i] * tke_0 - tr_yz_xzzz[i] - tr_xyz_zzz[i];

        tr_0_0_x_xyz_yyyy[i] = 2.0 * tr_xxyz_yyyy[i] * tbe_0 + 2.0 * tr_xyz_xyyyy[i] * tke_0 - tr_yz_yyyy[i];

        tr_0_0_x_xyz_yyyz[i] = 2.0 * tr_xxyz_yyyz[i] * tbe_0 + 2.0 * tr_xyz_xyyyz[i] * tke_0 - tr_yz_yyyz[i];

        tr_0_0_x_xyz_yyzz[i] = 2.0 * tr_xxyz_yyzz[i] * tbe_0 + 2.0 * tr_xyz_xyyzz[i] * tke_0 - tr_yz_yyzz[i];

        tr_0_0_x_xyz_yzzz[i] = 2.0 * tr_xxyz_yzzz[i] * tbe_0 + 2.0 * tr_xyz_xyzzz[i] * tke_0 - tr_yz_yzzz[i];

        tr_0_0_x_xyz_zzzz[i] = 2.0 * tr_xxyz_zzzz[i] * tbe_0 + 2.0 * tr_xyz_xzzzz[i] * tke_0 - tr_yz_zzzz[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto tr_0_0_x_xzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 75);

    auto tr_0_0_x_xzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 76);

    auto tr_0_0_x_xzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 77);

    auto tr_0_0_x_xzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 78);

    auto tr_0_0_x_xzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 79);

    auto tr_0_0_x_xzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 80);

    auto tr_0_0_x_xzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 81);

    auto tr_0_0_x_xzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 82);

    auto tr_0_0_x_xzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 83);

    auto tr_0_0_x_xzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 84);

    auto tr_0_0_x_xzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 85);

    auto tr_0_0_x_xzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 86);

    auto tr_0_0_x_xzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 87);

    auto tr_0_0_x_xzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 88);

    auto tr_0_0_x_xzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 89);

    #pragma omp simd aligned(tr_0_0_x_xzz_xxxx, tr_0_0_x_xzz_xxxy, tr_0_0_x_xzz_xxxz, tr_0_0_x_xzz_xxyy, tr_0_0_x_xzz_xxyz, tr_0_0_x_xzz_xxzz, tr_0_0_x_xzz_xyyy, tr_0_0_x_xzz_xyyz, tr_0_0_x_xzz_xyzz, tr_0_0_x_xzz_xzzz, tr_0_0_x_xzz_yyyy, tr_0_0_x_xzz_yyyz, tr_0_0_x_xzz_yyzz, tr_0_0_x_xzz_yzzz, tr_0_0_x_xzz_zzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxx, tr_xzz_xxxxy, tr_xzz_xxxxz, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzz_xxxx[i] = 2.0 * tr_xxzz_xxxx[i] * tbe_0 + 2.0 * tr_xzz_xxxxx[i] * tke_0 - tr_zz_xxxx[i] - 4.0 * tr_xzz_xxx[i];

        tr_0_0_x_xzz_xxxy[i] = 2.0 * tr_xxzz_xxxy[i] * tbe_0 + 2.0 * tr_xzz_xxxxy[i] * tke_0 - tr_zz_xxxy[i] - 3.0 * tr_xzz_xxy[i];

        tr_0_0_x_xzz_xxxz[i] = 2.0 * tr_xxzz_xxxz[i] * tbe_0 + 2.0 * tr_xzz_xxxxz[i] * tke_0 - tr_zz_xxxz[i] - 3.0 * tr_xzz_xxz[i];

        tr_0_0_x_xzz_xxyy[i] = 2.0 * tr_xxzz_xxyy[i] * tbe_0 + 2.0 * tr_xzz_xxxyy[i] * tke_0 - tr_zz_xxyy[i] - 2.0 * tr_xzz_xyy[i];

        tr_0_0_x_xzz_xxyz[i] = 2.0 * tr_xxzz_xxyz[i] * tbe_0 + 2.0 * tr_xzz_xxxyz[i] * tke_0 - tr_zz_xxyz[i] - 2.0 * tr_xzz_xyz[i];

        tr_0_0_x_xzz_xxzz[i] = 2.0 * tr_xxzz_xxzz[i] * tbe_0 + 2.0 * tr_xzz_xxxzz[i] * tke_0 - tr_zz_xxzz[i] - 2.0 * tr_xzz_xzz[i];

        tr_0_0_x_xzz_xyyy[i] = 2.0 * tr_xxzz_xyyy[i] * tbe_0 + 2.0 * tr_xzz_xxyyy[i] * tke_0 - tr_zz_xyyy[i] - tr_xzz_yyy[i];

        tr_0_0_x_xzz_xyyz[i] = 2.0 * tr_xxzz_xyyz[i] * tbe_0 + 2.0 * tr_xzz_xxyyz[i] * tke_0 - tr_zz_xyyz[i] - tr_xzz_yyz[i];

        tr_0_0_x_xzz_xyzz[i] = 2.0 * tr_xxzz_xyzz[i] * tbe_0 + 2.0 * tr_xzz_xxyzz[i] * tke_0 - tr_zz_xyzz[i] - tr_xzz_yzz[i];

        tr_0_0_x_xzz_xzzz[i] = 2.0 * tr_xxzz_xzzz[i] * tbe_0 + 2.0 * tr_xzz_xxzzz[i] * tke_0 - tr_zz_xzzz[i] - tr_xzz_zzz[i];

        tr_0_0_x_xzz_yyyy[i] = 2.0 * tr_xxzz_yyyy[i] * tbe_0 + 2.0 * tr_xzz_xyyyy[i] * tke_0 - tr_zz_yyyy[i];

        tr_0_0_x_xzz_yyyz[i] = 2.0 * tr_xxzz_yyyz[i] * tbe_0 + 2.0 * tr_xzz_xyyyz[i] * tke_0 - tr_zz_yyyz[i];

        tr_0_0_x_xzz_yyzz[i] = 2.0 * tr_xxzz_yyzz[i] * tbe_0 + 2.0 * tr_xzz_xyyzz[i] * tke_0 - tr_zz_yyzz[i];

        tr_0_0_x_xzz_yzzz[i] = 2.0 * tr_xxzz_yzzz[i] * tbe_0 + 2.0 * tr_xzz_xyzzz[i] * tke_0 - tr_zz_yzzz[i];

        tr_0_0_x_xzz_zzzz[i] = 2.0 * tr_xxzz_zzzz[i] * tbe_0 + 2.0 * tr_xzz_xzzzz[i] * tke_0 - tr_zz_zzzz[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto tr_0_0_x_yyy_xxxx = pbuffer.data(idx_op_geom_010_fg + 90);

    auto tr_0_0_x_yyy_xxxy = pbuffer.data(idx_op_geom_010_fg + 91);

    auto tr_0_0_x_yyy_xxxz = pbuffer.data(idx_op_geom_010_fg + 92);

    auto tr_0_0_x_yyy_xxyy = pbuffer.data(idx_op_geom_010_fg + 93);

    auto tr_0_0_x_yyy_xxyz = pbuffer.data(idx_op_geom_010_fg + 94);

    auto tr_0_0_x_yyy_xxzz = pbuffer.data(idx_op_geom_010_fg + 95);

    auto tr_0_0_x_yyy_xyyy = pbuffer.data(idx_op_geom_010_fg + 96);

    auto tr_0_0_x_yyy_xyyz = pbuffer.data(idx_op_geom_010_fg + 97);

    auto tr_0_0_x_yyy_xyzz = pbuffer.data(idx_op_geom_010_fg + 98);

    auto tr_0_0_x_yyy_xzzz = pbuffer.data(idx_op_geom_010_fg + 99);

    auto tr_0_0_x_yyy_yyyy = pbuffer.data(idx_op_geom_010_fg + 100);

    auto tr_0_0_x_yyy_yyyz = pbuffer.data(idx_op_geom_010_fg + 101);

    auto tr_0_0_x_yyy_yyzz = pbuffer.data(idx_op_geom_010_fg + 102);

    auto tr_0_0_x_yyy_yzzz = pbuffer.data(idx_op_geom_010_fg + 103);

    auto tr_0_0_x_yyy_zzzz = pbuffer.data(idx_op_geom_010_fg + 104);

    #pragma omp simd aligned(tr_0_0_x_yyy_xxxx, tr_0_0_x_yyy_xxxy, tr_0_0_x_yyy_xxxz, tr_0_0_x_yyy_xxyy, tr_0_0_x_yyy_xxyz, tr_0_0_x_yyy_xxzz, tr_0_0_x_yyy_xyyy, tr_0_0_x_yyy_xyyz, tr_0_0_x_yyy_xyzz, tr_0_0_x_yyy_xzzz, tr_0_0_x_yyy_yyyy, tr_0_0_x_yyy_yyyz, tr_0_0_x_yyy_yyzz, tr_0_0_x_yyy_yzzz, tr_0_0_x_yyy_zzzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, tr_yyy_xxx, tr_yyy_xxxxx, tr_yyy_xxxxy, tr_yyy_xxxxz, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyy_xxxx[i] = 2.0 * tr_xyyy_xxxx[i] * tbe_0 + 2.0 * tr_yyy_xxxxx[i] * tke_0 - 4.0 * tr_yyy_xxx[i];

        tr_0_0_x_yyy_xxxy[i] = 2.0 * tr_xyyy_xxxy[i] * tbe_0 + 2.0 * tr_yyy_xxxxy[i] * tke_0 - 3.0 * tr_yyy_xxy[i];

        tr_0_0_x_yyy_xxxz[i] = 2.0 * tr_xyyy_xxxz[i] * tbe_0 + 2.0 * tr_yyy_xxxxz[i] * tke_0 - 3.0 * tr_yyy_xxz[i];

        tr_0_0_x_yyy_xxyy[i] = 2.0 * tr_xyyy_xxyy[i] * tbe_0 + 2.0 * tr_yyy_xxxyy[i] * tke_0 - 2.0 * tr_yyy_xyy[i];

        tr_0_0_x_yyy_xxyz[i] = 2.0 * tr_xyyy_xxyz[i] * tbe_0 + 2.0 * tr_yyy_xxxyz[i] * tke_0 - 2.0 * tr_yyy_xyz[i];

        tr_0_0_x_yyy_xxzz[i] = 2.0 * tr_xyyy_xxzz[i] * tbe_0 + 2.0 * tr_yyy_xxxzz[i] * tke_0 - 2.0 * tr_yyy_xzz[i];

        tr_0_0_x_yyy_xyyy[i] = 2.0 * tr_xyyy_xyyy[i] * tbe_0 + 2.0 * tr_yyy_xxyyy[i] * tke_0 - tr_yyy_yyy[i];

        tr_0_0_x_yyy_xyyz[i] = 2.0 * tr_xyyy_xyyz[i] * tbe_0 + 2.0 * tr_yyy_xxyyz[i] * tke_0 - tr_yyy_yyz[i];

        tr_0_0_x_yyy_xyzz[i] = 2.0 * tr_xyyy_xyzz[i] * tbe_0 + 2.0 * tr_yyy_xxyzz[i] * tke_0 - tr_yyy_yzz[i];

        tr_0_0_x_yyy_xzzz[i] = 2.0 * tr_xyyy_xzzz[i] * tbe_0 + 2.0 * tr_yyy_xxzzz[i] * tke_0 - tr_yyy_zzz[i];

        tr_0_0_x_yyy_yyyy[i] = 2.0 * tr_xyyy_yyyy[i] * tbe_0 + 2.0 * tr_yyy_xyyyy[i] * tke_0;

        tr_0_0_x_yyy_yyyz[i] = 2.0 * tr_xyyy_yyyz[i] * tbe_0 + 2.0 * tr_yyy_xyyyz[i] * tke_0;

        tr_0_0_x_yyy_yyzz[i] = 2.0 * tr_xyyy_yyzz[i] * tbe_0 + 2.0 * tr_yyy_xyyzz[i] * tke_0;

        tr_0_0_x_yyy_yzzz[i] = 2.0 * tr_xyyy_yzzz[i] * tbe_0 + 2.0 * tr_yyy_xyzzz[i] * tke_0;

        tr_0_0_x_yyy_zzzz[i] = 2.0 * tr_xyyy_zzzz[i] * tbe_0 + 2.0 * tr_yyy_xzzzz[i] * tke_0;
    }

    // Set up 105-120 components of targeted buffer : FG

    auto tr_0_0_x_yyz_xxxx = pbuffer.data(idx_op_geom_010_fg + 105);

    auto tr_0_0_x_yyz_xxxy = pbuffer.data(idx_op_geom_010_fg + 106);

    auto tr_0_0_x_yyz_xxxz = pbuffer.data(idx_op_geom_010_fg + 107);

    auto tr_0_0_x_yyz_xxyy = pbuffer.data(idx_op_geom_010_fg + 108);

    auto tr_0_0_x_yyz_xxyz = pbuffer.data(idx_op_geom_010_fg + 109);

    auto tr_0_0_x_yyz_xxzz = pbuffer.data(idx_op_geom_010_fg + 110);

    auto tr_0_0_x_yyz_xyyy = pbuffer.data(idx_op_geom_010_fg + 111);

    auto tr_0_0_x_yyz_xyyz = pbuffer.data(idx_op_geom_010_fg + 112);

    auto tr_0_0_x_yyz_xyzz = pbuffer.data(idx_op_geom_010_fg + 113);

    auto tr_0_0_x_yyz_xzzz = pbuffer.data(idx_op_geom_010_fg + 114);

    auto tr_0_0_x_yyz_yyyy = pbuffer.data(idx_op_geom_010_fg + 115);

    auto tr_0_0_x_yyz_yyyz = pbuffer.data(idx_op_geom_010_fg + 116);

    auto tr_0_0_x_yyz_yyzz = pbuffer.data(idx_op_geom_010_fg + 117);

    auto tr_0_0_x_yyz_yzzz = pbuffer.data(idx_op_geom_010_fg + 118);

    auto tr_0_0_x_yyz_zzzz = pbuffer.data(idx_op_geom_010_fg + 119);

    #pragma omp simd aligned(tr_0_0_x_yyz_xxxx, tr_0_0_x_yyz_xxxy, tr_0_0_x_yyz_xxxz, tr_0_0_x_yyz_xxyy, tr_0_0_x_yyz_xxyz, tr_0_0_x_yyz_xxzz, tr_0_0_x_yyz_xyyy, tr_0_0_x_yyz_xyyz, tr_0_0_x_yyz_xyzz, tr_0_0_x_yyz_xzzz, tr_0_0_x_yyz_yyyy, tr_0_0_x_yyz_yyyz, tr_0_0_x_yyz_yyzz, tr_0_0_x_yyz_yzzz, tr_0_0_x_yyz_zzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxx, tr_yyz_xxxxy, tr_yyz_xxxxz, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyz_xxxx[i] = 2.0 * tr_xyyz_xxxx[i] * tbe_0 + 2.0 * tr_yyz_xxxxx[i] * tke_0 - 4.0 * tr_yyz_xxx[i];

        tr_0_0_x_yyz_xxxy[i] = 2.0 * tr_xyyz_xxxy[i] * tbe_0 + 2.0 * tr_yyz_xxxxy[i] * tke_0 - 3.0 * tr_yyz_xxy[i];

        tr_0_0_x_yyz_xxxz[i] = 2.0 * tr_xyyz_xxxz[i] * tbe_0 + 2.0 * tr_yyz_xxxxz[i] * tke_0 - 3.0 * tr_yyz_xxz[i];

        tr_0_0_x_yyz_xxyy[i] = 2.0 * tr_xyyz_xxyy[i] * tbe_0 + 2.0 * tr_yyz_xxxyy[i] * tke_0 - 2.0 * tr_yyz_xyy[i];

        tr_0_0_x_yyz_xxyz[i] = 2.0 * tr_xyyz_xxyz[i] * tbe_0 + 2.0 * tr_yyz_xxxyz[i] * tke_0 - 2.0 * tr_yyz_xyz[i];

        tr_0_0_x_yyz_xxzz[i] = 2.0 * tr_xyyz_xxzz[i] * tbe_0 + 2.0 * tr_yyz_xxxzz[i] * tke_0 - 2.0 * tr_yyz_xzz[i];

        tr_0_0_x_yyz_xyyy[i] = 2.0 * tr_xyyz_xyyy[i] * tbe_0 + 2.0 * tr_yyz_xxyyy[i] * tke_0 - tr_yyz_yyy[i];

        tr_0_0_x_yyz_xyyz[i] = 2.0 * tr_xyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyz_xxyyz[i] * tke_0 - tr_yyz_yyz[i];

        tr_0_0_x_yyz_xyzz[i] = 2.0 * tr_xyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyz_xxyzz[i] * tke_0 - tr_yyz_yzz[i];

        tr_0_0_x_yyz_xzzz[i] = 2.0 * tr_xyyz_xzzz[i] * tbe_0 + 2.0 * tr_yyz_xxzzz[i] * tke_0 - tr_yyz_zzz[i];

        tr_0_0_x_yyz_yyyy[i] = 2.0 * tr_xyyz_yyyy[i] * tbe_0 + 2.0 * tr_yyz_xyyyy[i] * tke_0;

        tr_0_0_x_yyz_yyyz[i] = 2.0 * tr_xyyz_yyyz[i] * tbe_0 + 2.0 * tr_yyz_xyyyz[i] * tke_0;

        tr_0_0_x_yyz_yyzz[i] = 2.0 * tr_xyyz_yyzz[i] * tbe_0 + 2.0 * tr_yyz_xyyzz[i] * tke_0;

        tr_0_0_x_yyz_yzzz[i] = 2.0 * tr_xyyz_yzzz[i] * tbe_0 + 2.0 * tr_yyz_xyzzz[i] * tke_0;

        tr_0_0_x_yyz_zzzz[i] = 2.0 * tr_xyyz_zzzz[i] * tbe_0 + 2.0 * tr_yyz_xzzzz[i] * tke_0;
    }

    // Set up 120-135 components of targeted buffer : FG

    auto tr_0_0_x_yzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 120);

    auto tr_0_0_x_yzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 121);

    auto tr_0_0_x_yzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 122);

    auto tr_0_0_x_yzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 123);

    auto tr_0_0_x_yzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 124);

    auto tr_0_0_x_yzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 125);

    auto tr_0_0_x_yzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 126);

    auto tr_0_0_x_yzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 127);

    auto tr_0_0_x_yzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 128);

    auto tr_0_0_x_yzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 129);

    auto tr_0_0_x_yzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 130);

    auto tr_0_0_x_yzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 131);

    auto tr_0_0_x_yzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 132);

    auto tr_0_0_x_yzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 133);

    auto tr_0_0_x_yzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 134);

    #pragma omp simd aligned(tr_0_0_x_yzz_xxxx, tr_0_0_x_yzz_xxxy, tr_0_0_x_yzz_xxxz, tr_0_0_x_yzz_xxyy, tr_0_0_x_yzz_xxyz, tr_0_0_x_yzz_xxzz, tr_0_0_x_yzz_xyyy, tr_0_0_x_yzz_xyyz, tr_0_0_x_yzz_xyzz, tr_0_0_x_yzz_xzzz, tr_0_0_x_yzz_yyyy, tr_0_0_x_yzz_yyyz, tr_0_0_x_yzz_yyzz, tr_0_0_x_yzz_yzzz, tr_0_0_x_yzz_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxx, tr_yzz_xxxxy, tr_yzz_xxxxz, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzz_xxxx[i] = 2.0 * tr_xyzz_xxxx[i] * tbe_0 + 2.0 * tr_yzz_xxxxx[i] * tke_0 - 4.0 * tr_yzz_xxx[i];

        tr_0_0_x_yzz_xxxy[i] = 2.0 * tr_xyzz_xxxy[i] * tbe_0 + 2.0 * tr_yzz_xxxxy[i] * tke_0 - 3.0 * tr_yzz_xxy[i];

        tr_0_0_x_yzz_xxxz[i] = 2.0 * tr_xyzz_xxxz[i] * tbe_0 + 2.0 * tr_yzz_xxxxz[i] * tke_0 - 3.0 * tr_yzz_xxz[i];

        tr_0_0_x_yzz_xxyy[i] = 2.0 * tr_xyzz_xxyy[i] * tbe_0 + 2.0 * tr_yzz_xxxyy[i] * tke_0 - 2.0 * tr_yzz_xyy[i];

        tr_0_0_x_yzz_xxyz[i] = 2.0 * tr_xyzz_xxyz[i] * tbe_0 + 2.0 * tr_yzz_xxxyz[i] * tke_0 - 2.0 * tr_yzz_xyz[i];

        tr_0_0_x_yzz_xxzz[i] = 2.0 * tr_xyzz_xxzz[i] * tbe_0 + 2.0 * tr_yzz_xxxzz[i] * tke_0 - 2.0 * tr_yzz_xzz[i];

        tr_0_0_x_yzz_xyyy[i] = 2.0 * tr_xyzz_xyyy[i] * tbe_0 + 2.0 * tr_yzz_xxyyy[i] * tke_0 - tr_yzz_yyy[i];

        tr_0_0_x_yzz_xyyz[i] = 2.0 * tr_xyzz_xyyz[i] * tbe_0 + 2.0 * tr_yzz_xxyyz[i] * tke_0 - tr_yzz_yyz[i];

        tr_0_0_x_yzz_xyzz[i] = 2.0 * tr_xyzz_xyzz[i] * tbe_0 + 2.0 * tr_yzz_xxyzz[i] * tke_0 - tr_yzz_yzz[i];

        tr_0_0_x_yzz_xzzz[i] = 2.0 * tr_xyzz_xzzz[i] * tbe_0 + 2.0 * tr_yzz_xxzzz[i] * tke_0 - tr_yzz_zzz[i];

        tr_0_0_x_yzz_yyyy[i] = 2.0 * tr_xyzz_yyyy[i] * tbe_0 + 2.0 * tr_yzz_xyyyy[i] * tke_0;

        tr_0_0_x_yzz_yyyz[i] = 2.0 * tr_xyzz_yyyz[i] * tbe_0 + 2.0 * tr_yzz_xyyyz[i] * tke_0;

        tr_0_0_x_yzz_yyzz[i] = 2.0 * tr_xyzz_yyzz[i] * tbe_0 + 2.0 * tr_yzz_xyyzz[i] * tke_0;

        tr_0_0_x_yzz_yzzz[i] = 2.0 * tr_xyzz_yzzz[i] * tbe_0 + 2.0 * tr_yzz_xyzzz[i] * tke_0;

        tr_0_0_x_yzz_zzzz[i] = 2.0 * tr_xyzz_zzzz[i] * tbe_0 + 2.0 * tr_yzz_xzzzz[i] * tke_0;
    }

    // Set up 135-150 components of targeted buffer : FG

    auto tr_0_0_x_zzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 135);

    auto tr_0_0_x_zzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 136);

    auto tr_0_0_x_zzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 137);

    auto tr_0_0_x_zzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 138);

    auto tr_0_0_x_zzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 139);

    auto tr_0_0_x_zzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 140);

    auto tr_0_0_x_zzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 141);

    auto tr_0_0_x_zzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 142);

    auto tr_0_0_x_zzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 143);

    auto tr_0_0_x_zzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 144);

    auto tr_0_0_x_zzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 145);

    auto tr_0_0_x_zzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 146);

    auto tr_0_0_x_zzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 147);

    auto tr_0_0_x_zzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 148);

    auto tr_0_0_x_zzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 149);

    #pragma omp simd aligned(tr_0_0_x_zzz_xxxx, tr_0_0_x_zzz_xxxy, tr_0_0_x_zzz_xxxz, tr_0_0_x_zzz_xxyy, tr_0_0_x_zzz_xxyz, tr_0_0_x_zzz_xxzz, tr_0_0_x_zzz_xyyy, tr_0_0_x_zzz_xyyz, tr_0_0_x_zzz_xyzz, tr_0_0_x_zzz_xzzz, tr_0_0_x_zzz_yyyy, tr_0_0_x_zzz_yyyz, tr_0_0_x_zzz_yyzz, tr_0_0_x_zzz_yzzz, tr_0_0_x_zzz_zzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxx, tr_zzz_xxxxy, tr_zzz_xxxxz, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzz_xxxx[i] = 2.0 * tr_xzzz_xxxx[i] * tbe_0 + 2.0 * tr_zzz_xxxxx[i] * tke_0 - 4.0 * tr_zzz_xxx[i];

        tr_0_0_x_zzz_xxxy[i] = 2.0 * tr_xzzz_xxxy[i] * tbe_0 + 2.0 * tr_zzz_xxxxy[i] * tke_0 - 3.0 * tr_zzz_xxy[i];

        tr_0_0_x_zzz_xxxz[i] = 2.0 * tr_xzzz_xxxz[i] * tbe_0 + 2.0 * tr_zzz_xxxxz[i] * tke_0 - 3.0 * tr_zzz_xxz[i];

        tr_0_0_x_zzz_xxyy[i] = 2.0 * tr_xzzz_xxyy[i] * tbe_0 + 2.0 * tr_zzz_xxxyy[i] * tke_0 - 2.0 * tr_zzz_xyy[i];

        tr_0_0_x_zzz_xxyz[i] = 2.0 * tr_xzzz_xxyz[i] * tbe_0 + 2.0 * tr_zzz_xxxyz[i] * tke_0 - 2.0 * tr_zzz_xyz[i];

        tr_0_0_x_zzz_xxzz[i] = 2.0 * tr_xzzz_xxzz[i] * tbe_0 + 2.0 * tr_zzz_xxxzz[i] * tke_0 - 2.0 * tr_zzz_xzz[i];

        tr_0_0_x_zzz_xyyy[i] = 2.0 * tr_xzzz_xyyy[i] * tbe_0 + 2.0 * tr_zzz_xxyyy[i] * tke_0 - tr_zzz_yyy[i];

        tr_0_0_x_zzz_xyyz[i] = 2.0 * tr_xzzz_xyyz[i] * tbe_0 + 2.0 * tr_zzz_xxyyz[i] * tke_0 - tr_zzz_yyz[i];

        tr_0_0_x_zzz_xyzz[i] = 2.0 * tr_xzzz_xyzz[i] * tbe_0 + 2.0 * tr_zzz_xxyzz[i] * tke_0 - tr_zzz_yzz[i];

        tr_0_0_x_zzz_xzzz[i] = 2.0 * tr_xzzz_xzzz[i] * tbe_0 + 2.0 * tr_zzz_xxzzz[i] * tke_0 - tr_zzz_zzz[i];

        tr_0_0_x_zzz_yyyy[i] = 2.0 * tr_xzzz_yyyy[i] * tbe_0 + 2.0 * tr_zzz_xyyyy[i] * tke_0;

        tr_0_0_x_zzz_yyyz[i] = 2.0 * tr_xzzz_yyyz[i] * tbe_0 + 2.0 * tr_zzz_xyyyz[i] * tke_0;

        tr_0_0_x_zzz_yyzz[i] = 2.0 * tr_xzzz_yyzz[i] * tbe_0 + 2.0 * tr_zzz_xyyzz[i] * tke_0;

        tr_0_0_x_zzz_yzzz[i] = 2.0 * tr_xzzz_yzzz[i] * tbe_0 + 2.0 * tr_zzz_xyzzz[i] * tke_0;

        tr_0_0_x_zzz_zzzz[i] = 2.0 * tr_xzzz_zzzz[i] * tbe_0 + 2.0 * tr_zzz_xzzzz[i] * tke_0;
    }

    // Set up 150-165 components of targeted buffer : FG

    auto tr_0_0_y_xxx_xxxx = pbuffer.data(idx_op_geom_010_fg + 150);

    auto tr_0_0_y_xxx_xxxy = pbuffer.data(idx_op_geom_010_fg + 151);

    auto tr_0_0_y_xxx_xxxz = pbuffer.data(idx_op_geom_010_fg + 152);

    auto tr_0_0_y_xxx_xxyy = pbuffer.data(idx_op_geom_010_fg + 153);

    auto tr_0_0_y_xxx_xxyz = pbuffer.data(idx_op_geom_010_fg + 154);

    auto tr_0_0_y_xxx_xxzz = pbuffer.data(idx_op_geom_010_fg + 155);

    auto tr_0_0_y_xxx_xyyy = pbuffer.data(idx_op_geom_010_fg + 156);

    auto tr_0_0_y_xxx_xyyz = pbuffer.data(idx_op_geom_010_fg + 157);

    auto tr_0_0_y_xxx_xyzz = pbuffer.data(idx_op_geom_010_fg + 158);

    auto tr_0_0_y_xxx_xzzz = pbuffer.data(idx_op_geom_010_fg + 159);

    auto tr_0_0_y_xxx_yyyy = pbuffer.data(idx_op_geom_010_fg + 160);

    auto tr_0_0_y_xxx_yyyz = pbuffer.data(idx_op_geom_010_fg + 161);

    auto tr_0_0_y_xxx_yyzz = pbuffer.data(idx_op_geom_010_fg + 162);

    auto tr_0_0_y_xxx_yzzz = pbuffer.data(idx_op_geom_010_fg + 163);

    auto tr_0_0_y_xxx_zzzz = pbuffer.data(idx_op_geom_010_fg + 164);

    #pragma omp simd aligned(tr_0_0_y_xxx_xxxx, tr_0_0_y_xxx_xxxy, tr_0_0_y_xxx_xxxz, tr_0_0_y_xxx_xxyy, tr_0_0_y_xxx_xxyz, tr_0_0_y_xxx_xxzz, tr_0_0_y_xxx_xyyy, tr_0_0_y_xxx_xyyz, tr_0_0_y_xxx_xyzz, tr_0_0_y_xxx_xzzz, tr_0_0_y_xxx_yyyy, tr_0_0_y_xxx_yyyz, tr_0_0_y_xxx_yyzz, tr_0_0_y_xxx_yzzz, tr_0_0_y_xxx_zzzz, tr_xxx_xxx, tr_xxx_xxxxy, tr_xxx_xxxyy, tr_xxx_xxxyz, tr_xxx_xxy, tr_xxx_xxyyy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyyyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyyyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xzzz, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yzzz, tr_xxxy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxx_xxxx[i] = 2.0 * tr_xxxy_xxxx[i] * tbe_0 + 2.0 * tr_xxx_xxxxy[i] * tke_0;

        tr_0_0_y_xxx_xxxy[i] = 2.0 * tr_xxxy_xxxy[i] * tbe_0 + 2.0 * tr_xxx_xxxyy[i] * tke_0 - tr_xxx_xxx[i];

        tr_0_0_y_xxx_xxxz[i] = 2.0 * tr_xxxy_xxxz[i] * tbe_0 + 2.0 * tr_xxx_xxxyz[i] * tke_0;

        tr_0_0_y_xxx_xxyy[i] = 2.0 * tr_xxxy_xxyy[i] * tbe_0 + 2.0 * tr_xxx_xxyyy[i] * tke_0 - 2.0 * tr_xxx_xxy[i];

        tr_0_0_y_xxx_xxyz[i] = 2.0 * tr_xxxy_xxyz[i] * tbe_0 + 2.0 * tr_xxx_xxyyz[i] * tke_0 - tr_xxx_xxz[i];

        tr_0_0_y_xxx_xxzz[i] = 2.0 * tr_xxxy_xxzz[i] * tbe_0 + 2.0 * tr_xxx_xxyzz[i] * tke_0;

        tr_0_0_y_xxx_xyyy[i] = 2.0 * tr_xxxy_xyyy[i] * tbe_0 + 2.0 * tr_xxx_xyyyy[i] * tke_0 - 3.0 * tr_xxx_xyy[i];

        tr_0_0_y_xxx_xyyz[i] = 2.0 * tr_xxxy_xyyz[i] * tbe_0 + 2.0 * tr_xxx_xyyyz[i] * tke_0 - 2.0 * tr_xxx_xyz[i];

        tr_0_0_y_xxx_xyzz[i] = 2.0 * tr_xxxy_xyzz[i] * tbe_0 + 2.0 * tr_xxx_xyyzz[i] * tke_0 - tr_xxx_xzz[i];

        tr_0_0_y_xxx_xzzz[i] = 2.0 * tr_xxxy_xzzz[i] * tbe_0 + 2.0 * tr_xxx_xyzzz[i] * tke_0;

        tr_0_0_y_xxx_yyyy[i] = 2.0 * tr_xxxy_yyyy[i] * tbe_0 + 2.0 * tr_xxx_yyyyy[i] * tke_0 - 4.0 * tr_xxx_yyy[i];

        tr_0_0_y_xxx_yyyz[i] = 2.0 * tr_xxxy_yyyz[i] * tbe_0 + 2.0 * tr_xxx_yyyyz[i] * tke_0 - 3.0 * tr_xxx_yyz[i];

        tr_0_0_y_xxx_yyzz[i] = 2.0 * tr_xxxy_yyzz[i] * tbe_0 + 2.0 * tr_xxx_yyyzz[i] * tke_0 - 2.0 * tr_xxx_yzz[i];

        tr_0_0_y_xxx_yzzz[i] = 2.0 * tr_xxxy_yzzz[i] * tbe_0 + 2.0 * tr_xxx_yyzzz[i] * tke_0 - tr_xxx_zzz[i];

        tr_0_0_y_xxx_zzzz[i] = 2.0 * tr_xxxy_zzzz[i] * tbe_0 + 2.0 * tr_xxx_yzzzz[i] * tke_0;
    }

    // Set up 165-180 components of targeted buffer : FG

    auto tr_0_0_y_xxy_xxxx = pbuffer.data(idx_op_geom_010_fg + 165);

    auto tr_0_0_y_xxy_xxxy = pbuffer.data(idx_op_geom_010_fg + 166);

    auto tr_0_0_y_xxy_xxxz = pbuffer.data(idx_op_geom_010_fg + 167);

    auto tr_0_0_y_xxy_xxyy = pbuffer.data(idx_op_geom_010_fg + 168);

    auto tr_0_0_y_xxy_xxyz = pbuffer.data(idx_op_geom_010_fg + 169);

    auto tr_0_0_y_xxy_xxzz = pbuffer.data(idx_op_geom_010_fg + 170);

    auto tr_0_0_y_xxy_xyyy = pbuffer.data(idx_op_geom_010_fg + 171);

    auto tr_0_0_y_xxy_xyyz = pbuffer.data(idx_op_geom_010_fg + 172);

    auto tr_0_0_y_xxy_xyzz = pbuffer.data(idx_op_geom_010_fg + 173);

    auto tr_0_0_y_xxy_xzzz = pbuffer.data(idx_op_geom_010_fg + 174);

    auto tr_0_0_y_xxy_yyyy = pbuffer.data(idx_op_geom_010_fg + 175);

    auto tr_0_0_y_xxy_yyyz = pbuffer.data(idx_op_geom_010_fg + 176);

    auto tr_0_0_y_xxy_yyzz = pbuffer.data(idx_op_geom_010_fg + 177);

    auto tr_0_0_y_xxy_yzzz = pbuffer.data(idx_op_geom_010_fg + 178);

    auto tr_0_0_y_xxy_zzzz = pbuffer.data(idx_op_geom_010_fg + 179);

    #pragma omp simd aligned(tr_0_0_y_xxy_xxxx, tr_0_0_y_xxy_xxxy, tr_0_0_y_xxy_xxxz, tr_0_0_y_xxy_xxyy, tr_0_0_y_xxy_xxyz, tr_0_0_y_xxy_xxzz, tr_0_0_y_xxy_xyyy, tr_0_0_y_xxy_xyyz, tr_0_0_y_xxy_xyzz, tr_0_0_y_xxy_xzzz, tr_0_0_y_xxy_yyyy, tr_0_0_y_xxy_yyyz, tr_0_0_y_xxy_yyzz, tr_0_0_y_xxy_yzzz, tr_0_0_y_xxy_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxy_xxx, tr_xxy_xxxxy, tr_xxy_xxxyy, tr_xxy_xxxyz, tr_xxy_xxy, tr_xxy_xxyyy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyyyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyyyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xzzz, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yzzz, tr_xxyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxy_xxxx[i] = 2.0 * tr_xxyy_xxxx[i] * tbe_0 + 2.0 * tr_xxy_xxxxy[i] * tke_0 - tr_xx_xxxx[i];

        tr_0_0_y_xxy_xxxy[i] = 2.0 * tr_xxyy_xxxy[i] * tbe_0 + 2.0 * tr_xxy_xxxyy[i] * tke_0 - tr_xx_xxxy[i] - tr_xxy_xxx[i];

        tr_0_0_y_xxy_xxxz[i] = 2.0 * tr_xxyy_xxxz[i] * tbe_0 + 2.0 * tr_xxy_xxxyz[i] * tke_0 - tr_xx_xxxz[i];

        tr_0_0_y_xxy_xxyy[i] = 2.0 * tr_xxyy_xxyy[i] * tbe_0 + 2.0 * tr_xxy_xxyyy[i] * tke_0 - tr_xx_xxyy[i] - 2.0 * tr_xxy_xxy[i];

        tr_0_0_y_xxy_xxyz[i] = 2.0 * tr_xxyy_xxyz[i] * tbe_0 + 2.0 * tr_xxy_xxyyz[i] * tke_0 - tr_xx_xxyz[i] - tr_xxy_xxz[i];

        tr_0_0_y_xxy_xxzz[i] = 2.0 * tr_xxyy_xxzz[i] * tbe_0 + 2.0 * tr_xxy_xxyzz[i] * tke_0 - tr_xx_xxzz[i];

        tr_0_0_y_xxy_xyyy[i] = 2.0 * tr_xxyy_xyyy[i] * tbe_0 + 2.0 * tr_xxy_xyyyy[i] * tke_0 - tr_xx_xyyy[i] - 3.0 * tr_xxy_xyy[i];

        tr_0_0_y_xxy_xyyz[i] = 2.0 * tr_xxyy_xyyz[i] * tbe_0 + 2.0 * tr_xxy_xyyyz[i] * tke_0 - tr_xx_xyyz[i] - 2.0 * tr_xxy_xyz[i];

        tr_0_0_y_xxy_xyzz[i] = 2.0 * tr_xxyy_xyzz[i] * tbe_0 + 2.0 * tr_xxy_xyyzz[i] * tke_0 - tr_xx_xyzz[i] - tr_xxy_xzz[i];

        tr_0_0_y_xxy_xzzz[i] = 2.0 * tr_xxyy_xzzz[i] * tbe_0 + 2.0 * tr_xxy_xyzzz[i] * tke_0 - tr_xx_xzzz[i];

        tr_0_0_y_xxy_yyyy[i] = 2.0 * tr_xxyy_yyyy[i] * tbe_0 + 2.0 * tr_xxy_yyyyy[i] * tke_0 - tr_xx_yyyy[i] - 4.0 * tr_xxy_yyy[i];

        tr_0_0_y_xxy_yyyz[i] = 2.0 * tr_xxyy_yyyz[i] * tbe_0 + 2.0 * tr_xxy_yyyyz[i] * tke_0 - tr_xx_yyyz[i] - 3.0 * tr_xxy_yyz[i];

        tr_0_0_y_xxy_yyzz[i] = 2.0 * tr_xxyy_yyzz[i] * tbe_0 + 2.0 * tr_xxy_yyyzz[i] * tke_0 - tr_xx_yyzz[i] - 2.0 * tr_xxy_yzz[i];

        tr_0_0_y_xxy_yzzz[i] = 2.0 * tr_xxyy_yzzz[i] * tbe_0 + 2.0 * tr_xxy_yyzzz[i] * tke_0 - tr_xx_yzzz[i] - tr_xxy_zzz[i];

        tr_0_0_y_xxy_zzzz[i] = 2.0 * tr_xxyy_zzzz[i] * tbe_0 + 2.0 * tr_xxy_yzzzz[i] * tke_0 - tr_xx_zzzz[i];
    }

    // Set up 180-195 components of targeted buffer : FG

    auto tr_0_0_y_xxz_xxxx = pbuffer.data(idx_op_geom_010_fg + 180);

    auto tr_0_0_y_xxz_xxxy = pbuffer.data(idx_op_geom_010_fg + 181);

    auto tr_0_0_y_xxz_xxxz = pbuffer.data(idx_op_geom_010_fg + 182);

    auto tr_0_0_y_xxz_xxyy = pbuffer.data(idx_op_geom_010_fg + 183);

    auto tr_0_0_y_xxz_xxyz = pbuffer.data(idx_op_geom_010_fg + 184);

    auto tr_0_0_y_xxz_xxzz = pbuffer.data(idx_op_geom_010_fg + 185);

    auto tr_0_0_y_xxz_xyyy = pbuffer.data(idx_op_geom_010_fg + 186);

    auto tr_0_0_y_xxz_xyyz = pbuffer.data(idx_op_geom_010_fg + 187);

    auto tr_0_0_y_xxz_xyzz = pbuffer.data(idx_op_geom_010_fg + 188);

    auto tr_0_0_y_xxz_xzzz = pbuffer.data(idx_op_geom_010_fg + 189);

    auto tr_0_0_y_xxz_yyyy = pbuffer.data(idx_op_geom_010_fg + 190);

    auto tr_0_0_y_xxz_yyyz = pbuffer.data(idx_op_geom_010_fg + 191);

    auto tr_0_0_y_xxz_yyzz = pbuffer.data(idx_op_geom_010_fg + 192);

    auto tr_0_0_y_xxz_yzzz = pbuffer.data(idx_op_geom_010_fg + 193);

    auto tr_0_0_y_xxz_zzzz = pbuffer.data(idx_op_geom_010_fg + 194);

    #pragma omp simd aligned(tr_0_0_y_xxz_xxxx, tr_0_0_y_xxz_xxxy, tr_0_0_y_xxz_xxxz, tr_0_0_y_xxz_xxyy, tr_0_0_y_xxz_xxyz, tr_0_0_y_xxz_xxzz, tr_0_0_y_xxz_xyyy, tr_0_0_y_xxz_xyyz, tr_0_0_y_xxz_xyzz, tr_0_0_y_xxz_xzzz, tr_0_0_y_xxz_yyyy, tr_0_0_y_xxz_yyyz, tr_0_0_y_xxz_yyzz, tr_0_0_y_xxz_yzzz, tr_0_0_y_xxz_zzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, tr_xxz_xxx, tr_xxz_xxxxy, tr_xxz_xxxyy, tr_xxz_xxxyz, tr_xxz_xxy, tr_xxz_xxyyy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyyyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyyyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxz_xxxx[i] = 2.0 * tr_xxyz_xxxx[i] * tbe_0 + 2.0 * tr_xxz_xxxxy[i] * tke_0;

        tr_0_0_y_xxz_xxxy[i] = 2.0 * tr_xxyz_xxxy[i] * tbe_0 + 2.0 * tr_xxz_xxxyy[i] * tke_0 - tr_xxz_xxx[i];

        tr_0_0_y_xxz_xxxz[i] = 2.0 * tr_xxyz_xxxz[i] * tbe_0 + 2.0 * tr_xxz_xxxyz[i] * tke_0;

        tr_0_0_y_xxz_xxyy[i] = 2.0 * tr_xxyz_xxyy[i] * tbe_0 + 2.0 * tr_xxz_xxyyy[i] * tke_0 - 2.0 * tr_xxz_xxy[i];

        tr_0_0_y_xxz_xxyz[i] = 2.0 * tr_xxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxz_xxyyz[i] * tke_0 - tr_xxz_xxz[i];

        tr_0_0_y_xxz_xxzz[i] = 2.0 * tr_xxyz_xxzz[i] * tbe_0 + 2.0 * tr_xxz_xxyzz[i] * tke_0;

        tr_0_0_y_xxz_xyyy[i] = 2.0 * tr_xxyz_xyyy[i] * tbe_0 + 2.0 * tr_xxz_xyyyy[i] * tke_0 - 3.0 * tr_xxz_xyy[i];

        tr_0_0_y_xxz_xyyz[i] = 2.0 * tr_xxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxz_xyyyz[i] * tke_0 - 2.0 * tr_xxz_xyz[i];

        tr_0_0_y_xxz_xyzz[i] = 2.0 * tr_xxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxz_xyyzz[i] * tke_0 - tr_xxz_xzz[i];

        tr_0_0_y_xxz_xzzz[i] = 2.0 * tr_xxyz_xzzz[i] * tbe_0 + 2.0 * tr_xxz_xyzzz[i] * tke_0;

        tr_0_0_y_xxz_yyyy[i] = 2.0 * tr_xxyz_yyyy[i] * tbe_0 + 2.0 * tr_xxz_yyyyy[i] * tke_0 - 4.0 * tr_xxz_yyy[i];

        tr_0_0_y_xxz_yyyz[i] = 2.0 * tr_xxyz_yyyz[i] * tbe_0 + 2.0 * tr_xxz_yyyyz[i] * tke_0 - 3.0 * tr_xxz_yyz[i];

        tr_0_0_y_xxz_yyzz[i] = 2.0 * tr_xxyz_yyzz[i] * tbe_0 + 2.0 * tr_xxz_yyyzz[i] * tke_0 - 2.0 * tr_xxz_yzz[i];

        tr_0_0_y_xxz_yzzz[i] = 2.0 * tr_xxyz_yzzz[i] * tbe_0 + 2.0 * tr_xxz_yyzzz[i] * tke_0 - tr_xxz_zzz[i];

        tr_0_0_y_xxz_zzzz[i] = 2.0 * tr_xxyz_zzzz[i] * tbe_0 + 2.0 * tr_xxz_yzzzz[i] * tke_0;
    }

    // Set up 195-210 components of targeted buffer : FG

    auto tr_0_0_y_xyy_xxxx = pbuffer.data(idx_op_geom_010_fg + 195);

    auto tr_0_0_y_xyy_xxxy = pbuffer.data(idx_op_geom_010_fg + 196);

    auto tr_0_0_y_xyy_xxxz = pbuffer.data(idx_op_geom_010_fg + 197);

    auto tr_0_0_y_xyy_xxyy = pbuffer.data(idx_op_geom_010_fg + 198);

    auto tr_0_0_y_xyy_xxyz = pbuffer.data(idx_op_geom_010_fg + 199);

    auto tr_0_0_y_xyy_xxzz = pbuffer.data(idx_op_geom_010_fg + 200);

    auto tr_0_0_y_xyy_xyyy = pbuffer.data(idx_op_geom_010_fg + 201);

    auto tr_0_0_y_xyy_xyyz = pbuffer.data(idx_op_geom_010_fg + 202);

    auto tr_0_0_y_xyy_xyzz = pbuffer.data(idx_op_geom_010_fg + 203);

    auto tr_0_0_y_xyy_xzzz = pbuffer.data(idx_op_geom_010_fg + 204);

    auto tr_0_0_y_xyy_yyyy = pbuffer.data(idx_op_geom_010_fg + 205);

    auto tr_0_0_y_xyy_yyyz = pbuffer.data(idx_op_geom_010_fg + 206);

    auto tr_0_0_y_xyy_yyzz = pbuffer.data(idx_op_geom_010_fg + 207);

    auto tr_0_0_y_xyy_yzzz = pbuffer.data(idx_op_geom_010_fg + 208);

    auto tr_0_0_y_xyy_zzzz = pbuffer.data(idx_op_geom_010_fg + 209);

    #pragma omp simd aligned(tr_0_0_y_xyy_xxxx, tr_0_0_y_xyy_xxxy, tr_0_0_y_xyy_xxxz, tr_0_0_y_xyy_xxyy, tr_0_0_y_xyy_xxyz, tr_0_0_y_xyy_xxzz, tr_0_0_y_xyy_xyyy, tr_0_0_y_xyy_xyyz, tr_0_0_y_xyy_xyzz, tr_0_0_y_xyy_xzzz, tr_0_0_y_xyy_yyyy, tr_0_0_y_xyy_yyyz, tr_0_0_y_xyy_yyzz, tr_0_0_y_xyy_yzzz, tr_0_0_y_xyy_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyy_xxx, tr_xyy_xxxxy, tr_xyy_xxxyy, tr_xyy_xxxyz, tr_xyy_xxy, tr_xyy_xxyyy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyyyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyyyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xzzz, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yzzz, tr_xyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyy_xxxx[i] = 2.0 * tr_xyyy_xxxx[i] * tbe_0 + 2.0 * tr_xyy_xxxxy[i] * tke_0 - 2.0 * tr_xy_xxxx[i];

        tr_0_0_y_xyy_xxxy[i] = 2.0 * tr_xyyy_xxxy[i] * tbe_0 + 2.0 * tr_xyy_xxxyy[i] * tke_0 - 2.0 * tr_xy_xxxy[i] - tr_xyy_xxx[i];

        tr_0_0_y_xyy_xxxz[i] = 2.0 * tr_xyyy_xxxz[i] * tbe_0 + 2.0 * tr_xyy_xxxyz[i] * tke_0 - 2.0 * tr_xy_xxxz[i];

        tr_0_0_y_xyy_xxyy[i] = 2.0 * tr_xyyy_xxyy[i] * tbe_0 + 2.0 * tr_xyy_xxyyy[i] * tke_0 - 2.0 * tr_xy_xxyy[i] - 2.0 * tr_xyy_xxy[i];

        tr_0_0_y_xyy_xxyz[i] = 2.0 * tr_xyyy_xxyz[i] * tbe_0 + 2.0 * tr_xyy_xxyyz[i] * tke_0 - 2.0 * tr_xy_xxyz[i] - tr_xyy_xxz[i];

        tr_0_0_y_xyy_xxzz[i] = 2.0 * tr_xyyy_xxzz[i] * tbe_0 + 2.0 * tr_xyy_xxyzz[i] * tke_0 - 2.0 * tr_xy_xxzz[i];

        tr_0_0_y_xyy_xyyy[i] = 2.0 * tr_xyyy_xyyy[i] * tbe_0 + 2.0 * tr_xyy_xyyyy[i] * tke_0 - 2.0 * tr_xy_xyyy[i] - 3.0 * tr_xyy_xyy[i];

        tr_0_0_y_xyy_xyyz[i] = 2.0 * tr_xyyy_xyyz[i] * tbe_0 + 2.0 * tr_xyy_xyyyz[i] * tke_0 - 2.0 * tr_xy_xyyz[i] - 2.0 * tr_xyy_xyz[i];

        tr_0_0_y_xyy_xyzz[i] = 2.0 * tr_xyyy_xyzz[i] * tbe_0 + 2.0 * tr_xyy_xyyzz[i] * tke_0 - 2.0 * tr_xy_xyzz[i] - tr_xyy_xzz[i];

        tr_0_0_y_xyy_xzzz[i] = 2.0 * tr_xyyy_xzzz[i] * tbe_0 + 2.0 * tr_xyy_xyzzz[i] * tke_0 - 2.0 * tr_xy_xzzz[i];

        tr_0_0_y_xyy_yyyy[i] = 2.0 * tr_xyyy_yyyy[i] * tbe_0 + 2.0 * tr_xyy_yyyyy[i] * tke_0 - 2.0 * tr_xy_yyyy[i] - 4.0 * tr_xyy_yyy[i];

        tr_0_0_y_xyy_yyyz[i] = 2.0 * tr_xyyy_yyyz[i] * tbe_0 + 2.0 * tr_xyy_yyyyz[i] * tke_0 - 2.0 * tr_xy_yyyz[i] - 3.0 * tr_xyy_yyz[i];

        tr_0_0_y_xyy_yyzz[i] = 2.0 * tr_xyyy_yyzz[i] * tbe_0 + 2.0 * tr_xyy_yyyzz[i] * tke_0 - 2.0 * tr_xy_yyzz[i] - 2.0 * tr_xyy_yzz[i];

        tr_0_0_y_xyy_yzzz[i] = 2.0 * tr_xyyy_yzzz[i] * tbe_0 + 2.0 * tr_xyy_yyzzz[i] * tke_0 - 2.0 * tr_xy_yzzz[i] - tr_xyy_zzz[i];

        tr_0_0_y_xyy_zzzz[i] = 2.0 * tr_xyyy_zzzz[i] * tbe_0 + 2.0 * tr_xyy_yzzzz[i] * tke_0 - 2.0 * tr_xy_zzzz[i];
    }

    // Set up 210-225 components of targeted buffer : FG

    auto tr_0_0_y_xyz_xxxx = pbuffer.data(idx_op_geom_010_fg + 210);

    auto tr_0_0_y_xyz_xxxy = pbuffer.data(idx_op_geom_010_fg + 211);

    auto tr_0_0_y_xyz_xxxz = pbuffer.data(idx_op_geom_010_fg + 212);

    auto tr_0_0_y_xyz_xxyy = pbuffer.data(idx_op_geom_010_fg + 213);

    auto tr_0_0_y_xyz_xxyz = pbuffer.data(idx_op_geom_010_fg + 214);

    auto tr_0_0_y_xyz_xxzz = pbuffer.data(idx_op_geom_010_fg + 215);

    auto tr_0_0_y_xyz_xyyy = pbuffer.data(idx_op_geom_010_fg + 216);

    auto tr_0_0_y_xyz_xyyz = pbuffer.data(idx_op_geom_010_fg + 217);

    auto tr_0_0_y_xyz_xyzz = pbuffer.data(idx_op_geom_010_fg + 218);

    auto tr_0_0_y_xyz_xzzz = pbuffer.data(idx_op_geom_010_fg + 219);

    auto tr_0_0_y_xyz_yyyy = pbuffer.data(idx_op_geom_010_fg + 220);

    auto tr_0_0_y_xyz_yyyz = pbuffer.data(idx_op_geom_010_fg + 221);

    auto tr_0_0_y_xyz_yyzz = pbuffer.data(idx_op_geom_010_fg + 222);

    auto tr_0_0_y_xyz_yzzz = pbuffer.data(idx_op_geom_010_fg + 223);

    auto tr_0_0_y_xyz_zzzz = pbuffer.data(idx_op_geom_010_fg + 224);

    #pragma omp simd aligned(tr_0_0_y_xyz_xxxx, tr_0_0_y_xyz_xxxy, tr_0_0_y_xyz_xxxz, tr_0_0_y_xyz_xxyy, tr_0_0_y_xyz_xxyz, tr_0_0_y_xyz_xxzz, tr_0_0_y_xyz_xyyy, tr_0_0_y_xyz_xyyz, tr_0_0_y_xyz_xyzz, tr_0_0_y_xyz_xzzz, tr_0_0_y_xyz_yyyy, tr_0_0_y_xyz_yyyz, tr_0_0_y_xyz_yyzz, tr_0_0_y_xyz_yzzz, tr_0_0_y_xyz_zzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, tr_xyz_xxx, tr_xyz_xxxxy, tr_xyz_xxxyy, tr_xyz_xxxyz, tr_xyz_xxy, tr_xyz_xxyyy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyyyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyyyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyz_xxxx[i] = 2.0 * tr_xyyz_xxxx[i] * tbe_0 + 2.0 * tr_xyz_xxxxy[i] * tke_0 - tr_xz_xxxx[i];

        tr_0_0_y_xyz_xxxy[i] = 2.0 * tr_xyyz_xxxy[i] * tbe_0 + 2.0 * tr_xyz_xxxyy[i] * tke_0 - tr_xz_xxxy[i] - tr_xyz_xxx[i];

        tr_0_0_y_xyz_xxxz[i] = 2.0 * tr_xyyz_xxxz[i] * tbe_0 + 2.0 * tr_xyz_xxxyz[i] * tke_0 - tr_xz_xxxz[i];

        tr_0_0_y_xyz_xxyy[i] = 2.0 * tr_xyyz_xxyy[i] * tbe_0 + 2.0 * tr_xyz_xxyyy[i] * tke_0 - tr_xz_xxyy[i] - 2.0 * tr_xyz_xxy[i];

        tr_0_0_y_xyz_xxyz[i] = 2.0 * tr_xyyz_xxyz[i] * tbe_0 + 2.0 * tr_xyz_xxyyz[i] * tke_0 - tr_xz_xxyz[i] - tr_xyz_xxz[i];

        tr_0_0_y_xyz_xxzz[i] = 2.0 * tr_xyyz_xxzz[i] * tbe_0 + 2.0 * tr_xyz_xxyzz[i] * tke_0 - tr_xz_xxzz[i];

        tr_0_0_y_xyz_xyyy[i] = 2.0 * tr_xyyz_xyyy[i] * tbe_0 + 2.0 * tr_xyz_xyyyy[i] * tke_0 - tr_xz_xyyy[i] - 3.0 * tr_xyz_xyy[i];

        tr_0_0_y_xyz_xyyz[i] = 2.0 * tr_xyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyz_xyyyz[i] * tke_0 - tr_xz_xyyz[i] - 2.0 * tr_xyz_xyz[i];

        tr_0_0_y_xyz_xyzz[i] = 2.0 * tr_xyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyz_xyyzz[i] * tke_0 - tr_xz_xyzz[i] - tr_xyz_xzz[i];

        tr_0_0_y_xyz_xzzz[i] = 2.0 * tr_xyyz_xzzz[i] * tbe_0 + 2.0 * tr_xyz_xyzzz[i] * tke_0 - tr_xz_xzzz[i];

        tr_0_0_y_xyz_yyyy[i] = 2.0 * tr_xyyz_yyyy[i] * tbe_0 + 2.0 * tr_xyz_yyyyy[i] * tke_0 - tr_xz_yyyy[i] - 4.0 * tr_xyz_yyy[i];

        tr_0_0_y_xyz_yyyz[i] = 2.0 * tr_xyyz_yyyz[i] * tbe_0 + 2.0 * tr_xyz_yyyyz[i] * tke_0 - tr_xz_yyyz[i] - 3.0 * tr_xyz_yyz[i];

        tr_0_0_y_xyz_yyzz[i] = 2.0 * tr_xyyz_yyzz[i] * tbe_0 + 2.0 * tr_xyz_yyyzz[i] * tke_0 - tr_xz_yyzz[i] - 2.0 * tr_xyz_yzz[i];

        tr_0_0_y_xyz_yzzz[i] = 2.0 * tr_xyyz_yzzz[i] * tbe_0 + 2.0 * tr_xyz_yyzzz[i] * tke_0 - tr_xz_yzzz[i] - tr_xyz_zzz[i];

        tr_0_0_y_xyz_zzzz[i] = 2.0 * tr_xyyz_zzzz[i] * tbe_0 + 2.0 * tr_xyz_yzzzz[i] * tke_0 - tr_xz_zzzz[i];
    }

    // Set up 225-240 components of targeted buffer : FG

    auto tr_0_0_y_xzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 225);

    auto tr_0_0_y_xzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 226);

    auto tr_0_0_y_xzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 227);

    auto tr_0_0_y_xzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 228);

    auto tr_0_0_y_xzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 229);

    auto tr_0_0_y_xzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 230);

    auto tr_0_0_y_xzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 231);

    auto tr_0_0_y_xzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 232);

    auto tr_0_0_y_xzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 233);

    auto tr_0_0_y_xzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 234);

    auto tr_0_0_y_xzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 235);

    auto tr_0_0_y_xzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 236);

    auto tr_0_0_y_xzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 237);

    auto tr_0_0_y_xzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 238);

    auto tr_0_0_y_xzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 239);

    #pragma omp simd aligned(tr_0_0_y_xzz_xxxx, tr_0_0_y_xzz_xxxy, tr_0_0_y_xzz_xxxz, tr_0_0_y_xzz_xxyy, tr_0_0_y_xzz_xxyz, tr_0_0_y_xzz_xxzz, tr_0_0_y_xzz_xyyy, tr_0_0_y_xzz_xyyz, tr_0_0_y_xzz_xyzz, tr_0_0_y_xzz_xzzz, tr_0_0_y_xzz_yyyy, tr_0_0_y_xzz_yyyz, tr_0_0_y_xzz_yyzz, tr_0_0_y_xzz_yzzz, tr_0_0_y_xzz_zzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, tr_xzz_xxx, tr_xzz_xxxxy, tr_xzz_xxxyy, tr_xzz_xxxyz, tr_xzz_xxy, tr_xzz_xxyyy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyyyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyyyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzz_xxxx[i] = 2.0 * tr_xyzz_xxxx[i] * tbe_0 + 2.0 * tr_xzz_xxxxy[i] * tke_0;

        tr_0_0_y_xzz_xxxy[i] = 2.0 * tr_xyzz_xxxy[i] * tbe_0 + 2.0 * tr_xzz_xxxyy[i] * tke_0 - tr_xzz_xxx[i];

        tr_0_0_y_xzz_xxxz[i] = 2.0 * tr_xyzz_xxxz[i] * tbe_0 + 2.0 * tr_xzz_xxxyz[i] * tke_0;

        tr_0_0_y_xzz_xxyy[i] = 2.0 * tr_xyzz_xxyy[i] * tbe_0 + 2.0 * tr_xzz_xxyyy[i] * tke_0 - 2.0 * tr_xzz_xxy[i];

        tr_0_0_y_xzz_xxyz[i] = 2.0 * tr_xyzz_xxyz[i] * tbe_0 + 2.0 * tr_xzz_xxyyz[i] * tke_0 - tr_xzz_xxz[i];

        tr_0_0_y_xzz_xxzz[i] = 2.0 * tr_xyzz_xxzz[i] * tbe_0 + 2.0 * tr_xzz_xxyzz[i] * tke_0;

        tr_0_0_y_xzz_xyyy[i] = 2.0 * tr_xyzz_xyyy[i] * tbe_0 + 2.0 * tr_xzz_xyyyy[i] * tke_0 - 3.0 * tr_xzz_xyy[i];

        tr_0_0_y_xzz_xyyz[i] = 2.0 * tr_xyzz_xyyz[i] * tbe_0 + 2.0 * tr_xzz_xyyyz[i] * tke_0 - 2.0 * tr_xzz_xyz[i];

        tr_0_0_y_xzz_xyzz[i] = 2.0 * tr_xyzz_xyzz[i] * tbe_0 + 2.0 * tr_xzz_xyyzz[i] * tke_0 - tr_xzz_xzz[i];

        tr_0_0_y_xzz_xzzz[i] = 2.0 * tr_xyzz_xzzz[i] * tbe_0 + 2.0 * tr_xzz_xyzzz[i] * tke_0;

        tr_0_0_y_xzz_yyyy[i] = 2.0 * tr_xyzz_yyyy[i] * tbe_0 + 2.0 * tr_xzz_yyyyy[i] * tke_0 - 4.0 * tr_xzz_yyy[i];

        tr_0_0_y_xzz_yyyz[i] = 2.0 * tr_xyzz_yyyz[i] * tbe_0 + 2.0 * tr_xzz_yyyyz[i] * tke_0 - 3.0 * tr_xzz_yyz[i];

        tr_0_0_y_xzz_yyzz[i] = 2.0 * tr_xyzz_yyzz[i] * tbe_0 + 2.0 * tr_xzz_yyyzz[i] * tke_0 - 2.0 * tr_xzz_yzz[i];

        tr_0_0_y_xzz_yzzz[i] = 2.0 * tr_xyzz_yzzz[i] * tbe_0 + 2.0 * tr_xzz_yyzzz[i] * tke_0 - tr_xzz_zzz[i];

        tr_0_0_y_xzz_zzzz[i] = 2.0 * tr_xyzz_zzzz[i] * tbe_0 + 2.0 * tr_xzz_yzzzz[i] * tke_0;
    }

    // Set up 240-255 components of targeted buffer : FG

    auto tr_0_0_y_yyy_xxxx = pbuffer.data(idx_op_geom_010_fg + 240);

    auto tr_0_0_y_yyy_xxxy = pbuffer.data(idx_op_geom_010_fg + 241);

    auto tr_0_0_y_yyy_xxxz = pbuffer.data(idx_op_geom_010_fg + 242);

    auto tr_0_0_y_yyy_xxyy = pbuffer.data(idx_op_geom_010_fg + 243);

    auto tr_0_0_y_yyy_xxyz = pbuffer.data(idx_op_geom_010_fg + 244);

    auto tr_0_0_y_yyy_xxzz = pbuffer.data(idx_op_geom_010_fg + 245);

    auto tr_0_0_y_yyy_xyyy = pbuffer.data(idx_op_geom_010_fg + 246);

    auto tr_0_0_y_yyy_xyyz = pbuffer.data(idx_op_geom_010_fg + 247);

    auto tr_0_0_y_yyy_xyzz = pbuffer.data(idx_op_geom_010_fg + 248);

    auto tr_0_0_y_yyy_xzzz = pbuffer.data(idx_op_geom_010_fg + 249);

    auto tr_0_0_y_yyy_yyyy = pbuffer.data(idx_op_geom_010_fg + 250);

    auto tr_0_0_y_yyy_yyyz = pbuffer.data(idx_op_geom_010_fg + 251);

    auto tr_0_0_y_yyy_yyzz = pbuffer.data(idx_op_geom_010_fg + 252);

    auto tr_0_0_y_yyy_yzzz = pbuffer.data(idx_op_geom_010_fg + 253);

    auto tr_0_0_y_yyy_zzzz = pbuffer.data(idx_op_geom_010_fg + 254);

    #pragma omp simd aligned(tr_0_0_y_yyy_xxxx, tr_0_0_y_yyy_xxxy, tr_0_0_y_yyy_xxxz, tr_0_0_y_yyy_xxyy, tr_0_0_y_yyy_xxyz, tr_0_0_y_yyy_xxzz, tr_0_0_y_yyy_xyyy, tr_0_0_y_yyy_xyyz, tr_0_0_y_yyy_xyzz, tr_0_0_y_yyy_xzzz, tr_0_0_y_yyy_yyyy, tr_0_0_y_yyy_yyyz, tr_0_0_y_yyy_yyzz, tr_0_0_y_yyy_yzzz, tr_0_0_y_yyy_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyy_xxx, tr_yyy_xxxxy, tr_yyy_xxxyy, tr_yyy_xxxyz, tr_yyy_xxy, tr_yyy_xxyyy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyyyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyyyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xzzz, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yzzz, tr_yyyy_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyy_xxxx[i] = 2.0 * tr_yyyy_xxxx[i] * tbe_0 + 2.0 * tr_yyy_xxxxy[i] * tke_0 - 3.0 * tr_yy_xxxx[i];

        tr_0_0_y_yyy_xxxy[i] = 2.0 * tr_yyyy_xxxy[i] * tbe_0 + 2.0 * tr_yyy_xxxyy[i] * tke_0 - 3.0 * tr_yy_xxxy[i] - tr_yyy_xxx[i];

        tr_0_0_y_yyy_xxxz[i] = 2.0 * tr_yyyy_xxxz[i] * tbe_0 + 2.0 * tr_yyy_xxxyz[i] * tke_0 - 3.0 * tr_yy_xxxz[i];

        tr_0_0_y_yyy_xxyy[i] = 2.0 * tr_yyyy_xxyy[i] * tbe_0 + 2.0 * tr_yyy_xxyyy[i] * tke_0 - 3.0 * tr_yy_xxyy[i] - 2.0 * tr_yyy_xxy[i];

        tr_0_0_y_yyy_xxyz[i] = 2.0 * tr_yyyy_xxyz[i] * tbe_0 + 2.0 * tr_yyy_xxyyz[i] * tke_0 - 3.0 * tr_yy_xxyz[i] - tr_yyy_xxz[i];

        tr_0_0_y_yyy_xxzz[i] = 2.0 * tr_yyyy_xxzz[i] * tbe_0 + 2.0 * tr_yyy_xxyzz[i] * tke_0 - 3.0 * tr_yy_xxzz[i];

        tr_0_0_y_yyy_xyyy[i] = 2.0 * tr_yyyy_xyyy[i] * tbe_0 + 2.0 * tr_yyy_xyyyy[i] * tke_0 - 3.0 * tr_yy_xyyy[i] - 3.0 * tr_yyy_xyy[i];

        tr_0_0_y_yyy_xyyz[i] = 2.0 * tr_yyyy_xyyz[i] * tbe_0 + 2.0 * tr_yyy_xyyyz[i] * tke_0 - 3.0 * tr_yy_xyyz[i] - 2.0 * tr_yyy_xyz[i];

        tr_0_0_y_yyy_xyzz[i] = 2.0 * tr_yyyy_xyzz[i] * tbe_0 + 2.0 * tr_yyy_xyyzz[i] * tke_0 - 3.0 * tr_yy_xyzz[i] - tr_yyy_xzz[i];

        tr_0_0_y_yyy_xzzz[i] = 2.0 * tr_yyyy_xzzz[i] * tbe_0 + 2.0 * tr_yyy_xyzzz[i] * tke_0 - 3.0 * tr_yy_xzzz[i];

        tr_0_0_y_yyy_yyyy[i] = 2.0 * tr_yyyy_yyyy[i] * tbe_0 + 2.0 * tr_yyy_yyyyy[i] * tke_0 - 3.0 * tr_yy_yyyy[i] - 4.0 * tr_yyy_yyy[i];

        tr_0_0_y_yyy_yyyz[i] = 2.0 * tr_yyyy_yyyz[i] * tbe_0 + 2.0 * tr_yyy_yyyyz[i] * tke_0 - 3.0 * tr_yy_yyyz[i] - 3.0 * tr_yyy_yyz[i];

        tr_0_0_y_yyy_yyzz[i] = 2.0 * tr_yyyy_yyzz[i] * tbe_0 + 2.0 * tr_yyy_yyyzz[i] * tke_0 - 3.0 * tr_yy_yyzz[i] - 2.0 * tr_yyy_yzz[i];

        tr_0_0_y_yyy_yzzz[i] = 2.0 * tr_yyyy_yzzz[i] * tbe_0 + 2.0 * tr_yyy_yyzzz[i] * tke_0 - 3.0 * tr_yy_yzzz[i] - tr_yyy_zzz[i];

        tr_0_0_y_yyy_zzzz[i] = 2.0 * tr_yyyy_zzzz[i] * tbe_0 + 2.0 * tr_yyy_yzzzz[i] * tke_0 - 3.0 * tr_yy_zzzz[i];
    }

    // Set up 255-270 components of targeted buffer : FG

    auto tr_0_0_y_yyz_xxxx = pbuffer.data(idx_op_geom_010_fg + 255);

    auto tr_0_0_y_yyz_xxxy = pbuffer.data(idx_op_geom_010_fg + 256);

    auto tr_0_0_y_yyz_xxxz = pbuffer.data(idx_op_geom_010_fg + 257);

    auto tr_0_0_y_yyz_xxyy = pbuffer.data(idx_op_geom_010_fg + 258);

    auto tr_0_0_y_yyz_xxyz = pbuffer.data(idx_op_geom_010_fg + 259);

    auto tr_0_0_y_yyz_xxzz = pbuffer.data(idx_op_geom_010_fg + 260);

    auto tr_0_0_y_yyz_xyyy = pbuffer.data(idx_op_geom_010_fg + 261);

    auto tr_0_0_y_yyz_xyyz = pbuffer.data(idx_op_geom_010_fg + 262);

    auto tr_0_0_y_yyz_xyzz = pbuffer.data(idx_op_geom_010_fg + 263);

    auto tr_0_0_y_yyz_xzzz = pbuffer.data(idx_op_geom_010_fg + 264);

    auto tr_0_0_y_yyz_yyyy = pbuffer.data(idx_op_geom_010_fg + 265);

    auto tr_0_0_y_yyz_yyyz = pbuffer.data(idx_op_geom_010_fg + 266);

    auto tr_0_0_y_yyz_yyzz = pbuffer.data(idx_op_geom_010_fg + 267);

    auto tr_0_0_y_yyz_yzzz = pbuffer.data(idx_op_geom_010_fg + 268);

    auto tr_0_0_y_yyz_zzzz = pbuffer.data(idx_op_geom_010_fg + 269);

    #pragma omp simd aligned(tr_0_0_y_yyz_xxxx, tr_0_0_y_yyz_xxxy, tr_0_0_y_yyz_xxxz, tr_0_0_y_yyz_xxyy, tr_0_0_y_yyz_xxyz, tr_0_0_y_yyz_xxzz, tr_0_0_y_yyz_xyyy, tr_0_0_y_yyz_xyyz, tr_0_0_y_yyz_xyzz, tr_0_0_y_yyz_xzzz, tr_0_0_y_yyz_yyyy, tr_0_0_y_yyz_yyyz, tr_0_0_y_yyz_yyzz, tr_0_0_y_yyz_yzzz, tr_0_0_y_yyz_zzzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, tr_yyz_xxx, tr_yyz_xxxxy, tr_yyz_xxxyy, tr_yyz_xxxyz, tr_yyz_xxy, tr_yyz_xxyyy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyyyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyyyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyz_xxxx[i] = 2.0 * tr_yyyz_xxxx[i] * tbe_0 + 2.0 * tr_yyz_xxxxy[i] * tke_0 - 2.0 * tr_yz_xxxx[i];

        tr_0_0_y_yyz_xxxy[i] = 2.0 * tr_yyyz_xxxy[i] * tbe_0 + 2.0 * tr_yyz_xxxyy[i] * tke_0 - 2.0 * tr_yz_xxxy[i] - tr_yyz_xxx[i];

        tr_0_0_y_yyz_xxxz[i] = 2.0 * tr_yyyz_xxxz[i] * tbe_0 + 2.0 * tr_yyz_xxxyz[i] * tke_0 - 2.0 * tr_yz_xxxz[i];

        tr_0_0_y_yyz_xxyy[i] = 2.0 * tr_yyyz_xxyy[i] * tbe_0 + 2.0 * tr_yyz_xxyyy[i] * tke_0 - 2.0 * tr_yz_xxyy[i] - 2.0 * tr_yyz_xxy[i];

        tr_0_0_y_yyz_xxyz[i] = 2.0 * tr_yyyz_xxyz[i] * tbe_0 + 2.0 * tr_yyz_xxyyz[i] * tke_0 - 2.0 * tr_yz_xxyz[i] - tr_yyz_xxz[i];

        tr_0_0_y_yyz_xxzz[i] = 2.0 * tr_yyyz_xxzz[i] * tbe_0 + 2.0 * tr_yyz_xxyzz[i] * tke_0 - 2.0 * tr_yz_xxzz[i];

        tr_0_0_y_yyz_xyyy[i] = 2.0 * tr_yyyz_xyyy[i] * tbe_0 + 2.0 * tr_yyz_xyyyy[i] * tke_0 - 2.0 * tr_yz_xyyy[i] - 3.0 * tr_yyz_xyy[i];

        tr_0_0_y_yyz_xyyz[i] = 2.0 * tr_yyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyz_xyyyz[i] * tke_0 - 2.0 * tr_yz_xyyz[i] - 2.0 * tr_yyz_xyz[i];

        tr_0_0_y_yyz_xyzz[i] = 2.0 * tr_yyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyz_xyyzz[i] * tke_0 - 2.0 * tr_yz_xyzz[i] - tr_yyz_xzz[i];

        tr_0_0_y_yyz_xzzz[i] = 2.0 * tr_yyyz_xzzz[i] * tbe_0 + 2.0 * tr_yyz_xyzzz[i] * tke_0 - 2.0 * tr_yz_xzzz[i];

        tr_0_0_y_yyz_yyyy[i] = 2.0 * tr_yyyz_yyyy[i] * tbe_0 + 2.0 * tr_yyz_yyyyy[i] * tke_0 - 2.0 * tr_yz_yyyy[i] - 4.0 * tr_yyz_yyy[i];

        tr_0_0_y_yyz_yyyz[i] = 2.0 * tr_yyyz_yyyz[i] * tbe_0 + 2.0 * tr_yyz_yyyyz[i] * tke_0 - 2.0 * tr_yz_yyyz[i] - 3.0 * tr_yyz_yyz[i];

        tr_0_0_y_yyz_yyzz[i] = 2.0 * tr_yyyz_yyzz[i] * tbe_0 + 2.0 * tr_yyz_yyyzz[i] * tke_0 - 2.0 * tr_yz_yyzz[i] - 2.0 * tr_yyz_yzz[i];

        tr_0_0_y_yyz_yzzz[i] = 2.0 * tr_yyyz_yzzz[i] * tbe_0 + 2.0 * tr_yyz_yyzzz[i] * tke_0 - 2.0 * tr_yz_yzzz[i] - tr_yyz_zzz[i];

        tr_0_0_y_yyz_zzzz[i] = 2.0 * tr_yyyz_zzzz[i] * tbe_0 + 2.0 * tr_yyz_yzzzz[i] * tke_0 - 2.0 * tr_yz_zzzz[i];
    }

    // Set up 270-285 components of targeted buffer : FG

    auto tr_0_0_y_yzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 270);

    auto tr_0_0_y_yzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 271);

    auto tr_0_0_y_yzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 272);

    auto tr_0_0_y_yzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 273);

    auto tr_0_0_y_yzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 274);

    auto tr_0_0_y_yzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 275);

    auto tr_0_0_y_yzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 276);

    auto tr_0_0_y_yzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 277);

    auto tr_0_0_y_yzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 278);

    auto tr_0_0_y_yzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 279);

    auto tr_0_0_y_yzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 280);

    auto tr_0_0_y_yzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 281);

    auto tr_0_0_y_yzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 282);

    auto tr_0_0_y_yzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 283);

    auto tr_0_0_y_yzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 284);

    #pragma omp simd aligned(tr_0_0_y_yzz_xxxx, tr_0_0_y_yzz_xxxy, tr_0_0_y_yzz_xxxz, tr_0_0_y_yzz_xxyy, tr_0_0_y_yzz_xxyz, tr_0_0_y_yzz_xxzz, tr_0_0_y_yzz_xyyy, tr_0_0_y_yzz_xyyz, tr_0_0_y_yzz_xyzz, tr_0_0_y_yzz_xzzz, tr_0_0_y_yzz_yyyy, tr_0_0_y_yzz_yyyz, tr_0_0_y_yzz_yyzz, tr_0_0_y_yzz_yzzz, tr_0_0_y_yzz_zzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, tr_yzz_xxx, tr_yzz_xxxxy, tr_yzz_xxxyy, tr_yzz_xxxyz, tr_yzz_xxy, tr_yzz_xxyyy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyyyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyyyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzz_xxxx[i] = 2.0 * tr_yyzz_xxxx[i] * tbe_0 + 2.0 * tr_yzz_xxxxy[i] * tke_0 - tr_zz_xxxx[i];

        tr_0_0_y_yzz_xxxy[i] = 2.0 * tr_yyzz_xxxy[i] * tbe_0 + 2.0 * tr_yzz_xxxyy[i] * tke_0 - tr_zz_xxxy[i] - tr_yzz_xxx[i];

        tr_0_0_y_yzz_xxxz[i] = 2.0 * tr_yyzz_xxxz[i] * tbe_0 + 2.0 * tr_yzz_xxxyz[i] * tke_0 - tr_zz_xxxz[i];

        tr_0_0_y_yzz_xxyy[i] = 2.0 * tr_yyzz_xxyy[i] * tbe_0 + 2.0 * tr_yzz_xxyyy[i] * tke_0 - tr_zz_xxyy[i] - 2.0 * tr_yzz_xxy[i];

        tr_0_0_y_yzz_xxyz[i] = 2.0 * tr_yyzz_xxyz[i] * tbe_0 + 2.0 * tr_yzz_xxyyz[i] * tke_0 - tr_zz_xxyz[i] - tr_yzz_xxz[i];

        tr_0_0_y_yzz_xxzz[i] = 2.0 * tr_yyzz_xxzz[i] * tbe_0 + 2.0 * tr_yzz_xxyzz[i] * tke_0 - tr_zz_xxzz[i];

        tr_0_0_y_yzz_xyyy[i] = 2.0 * tr_yyzz_xyyy[i] * tbe_0 + 2.0 * tr_yzz_xyyyy[i] * tke_0 - tr_zz_xyyy[i] - 3.0 * tr_yzz_xyy[i];

        tr_0_0_y_yzz_xyyz[i] = 2.0 * tr_yyzz_xyyz[i] * tbe_0 + 2.0 * tr_yzz_xyyyz[i] * tke_0 - tr_zz_xyyz[i] - 2.0 * tr_yzz_xyz[i];

        tr_0_0_y_yzz_xyzz[i] = 2.0 * tr_yyzz_xyzz[i] * tbe_0 + 2.0 * tr_yzz_xyyzz[i] * tke_0 - tr_zz_xyzz[i] - tr_yzz_xzz[i];

        tr_0_0_y_yzz_xzzz[i] = 2.0 * tr_yyzz_xzzz[i] * tbe_0 + 2.0 * tr_yzz_xyzzz[i] * tke_0 - tr_zz_xzzz[i];

        tr_0_0_y_yzz_yyyy[i] = 2.0 * tr_yyzz_yyyy[i] * tbe_0 + 2.0 * tr_yzz_yyyyy[i] * tke_0 - tr_zz_yyyy[i] - 4.0 * tr_yzz_yyy[i];

        tr_0_0_y_yzz_yyyz[i] = 2.0 * tr_yyzz_yyyz[i] * tbe_0 + 2.0 * tr_yzz_yyyyz[i] * tke_0 - tr_zz_yyyz[i] - 3.0 * tr_yzz_yyz[i];

        tr_0_0_y_yzz_yyzz[i] = 2.0 * tr_yyzz_yyzz[i] * tbe_0 + 2.0 * tr_yzz_yyyzz[i] * tke_0 - tr_zz_yyzz[i] - 2.0 * tr_yzz_yzz[i];

        tr_0_0_y_yzz_yzzz[i] = 2.0 * tr_yyzz_yzzz[i] * tbe_0 + 2.0 * tr_yzz_yyzzz[i] * tke_0 - tr_zz_yzzz[i] - tr_yzz_zzz[i];

        tr_0_0_y_yzz_zzzz[i] = 2.0 * tr_yyzz_zzzz[i] * tbe_0 + 2.0 * tr_yzz_yzzzz[i] * tke_0 - tr_zz_zzzz[i];
    }

    // Set up 285-300 components of targeted buffer : FG

    auto tr_0_0_y_zzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 285);

    auto tr_0_0_y_zzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 286);

    auto tr_0_0_y_zzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 287);

    auto tr_0_0_y_zzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 288);

    auto tr_0_0_y_zzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 289);

    auto tr_0_0_y_zzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 290);

    auto tr_0_0_y_zzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 291);

    auto tr_0_0_y_zzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 292);

    auto tr_0_0_y_zzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 293);

    auto tr_0_0_y_zzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 294);

    auto tr_0_0_y_zzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 295);

    auto tr_0_0_y_zzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 296);

    auto tr_0_0_y_zzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 297);

    auto tr_0_0_y_zzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 298);

    auto tr_0_0_y_zzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 299);

    #pragma omp simd aligned(tr_0_0_y_zzz_xxxx, tr_0_0_y_zzz_xxxy, tr_0_0_y_zzz_xxxz, tr_0_0_y_zzz_xxyy, tr_0_0_y_zzz_xxyz, tr_0_0_y_zzz_xxzz, tr_0_0_y_zzz_xyyy, tr_0_0_y_zzz_xyyz, tr_0_0_y_zzz_xyzz, tr_0_0_y_zzz_xzzz, tr_0_0_y_zzz_yyyy, tr_0_0_y_zzz_yyyz, tr_0_0_y_zzz_yyzz, tr_0_0_y_zzz_yzzz, tr_0_0_y_zzz_zzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, tr_zzz_xxx, tr_zzz_xxxxy, tr_zzz_xxxyy, tr_zzz_xxxyz, tr_zzz_xxy, tr_zzz_xxyyy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyyyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyyyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzz_xxxx[i] = 2.0 * tr_yzzz_xxxx[i] * tbe_0 + 2.0 * tr_zzz_xxxxy[i] * tke_0;

        tr_0_0_y_zzz_xxxy[i] = 2.0 * tr_yzzz_xxxy[i] * tbe_0 + 2.0 * tr_zzz_xxxyy[i] * tke_0 - tr_zzz_xxx[i];

        tr_0_0_y_zzz_xxxz[i] = 2.0 * tr_yzzz_xxxz[i] * tbe_0 + 2.0 * tr_zzz_xxxyz[i] * tke_0;

        tr_0_0_y_zzz_xxyy[i] = 2.0 * tr_yzzz_xxyy[i] * tbe_0 + 2.0 * tr_zzz_xxyyy[i] * tke_0 - 2.0 * tr_zzz_xxy[i];

        tr_0_0_y_zzz_xxyz[i] = 2.0 * tr_yzzz_xxyz[i] * tbe_0 + 2.0 * tr_zzz_xxyyz[i] * tke_0 - tr_zzz_xxz[i];

        tr_0_0_y_zzz_xxzz[i] = 2.0 * tr_yzzz_xxzz[i] * tbe_0 + 2.0 * tr_zzz_xxyzz[i] * tke_0;

        tr_0_0_y_zzz_xyyy[i] = 2.0 * tr_yzzz_xyyy[i] * tbe_0 + 2.0 * tr_zzz_xyyyy[i] * tke_0 - 3.0 * tr_zzz_xyy[i];

        tr_0_0_y_zzz_xyyz[i] = 2.0 * tr_yzzz_xyyz[i] * tbe_0 + 2.0 * tr_zzz_xyyyz[i] * tke_0 - 2.0 * tr_zzz_xyz[i];

        tr_0_0_y_zzz_xyzz[i] = 2.0 * tr_yzzz_xyzz[i] * tbe_0 + 2.0 * tr_zzz_xyyzz[i] * tke_0 - tr_zzz_xzz[i];

        tr_0_0_y_zzz_xzzz[i] = 2.0 * tr_yzzz_xzzz[i] * tbe_0 + 2.0 * tr_zzz_xyzzz[i] * tke_0;

        tr_0_0_y_zzz_yyyy[i] = 2.0 * tr_yzzz_yyyy[i] * tbe_0 + 2.0 * tr_zzz_yyyyy[i] * tke_0 - 4.0 * tr_zzz_yyy[i];

        tr_0_0_y_zzz_yyyz[i] = 2.0 * tr_yzzz_yyyz[i] * tbe_0 + 2.0 * tr_zzz_yyyyz[i] * tke_0 - 3.0 * tr_zzz_yyz[i];

        tr_0_0_y_zzz_yyzz[i] = 2.0 * tr_yzzz_yyzz[i] * tbe_0 + 2.0 * tr_zzz_yyyzz[i] * tke_0 - 2.0 * tr_zzz_yzz[i];

        tr_0_0_y_zzz_yzzz[i] = 2.0 * tr_yzzz_yzzz[i] * tbe_0 + 2.0 * tr_zzz_yyzzz[i] * tke_0 - tr_zzz_zzz[i];

        tr_0_0_y_zzz_zzzz[i] = 2.0 * tr_yzzz_zzzz[i] * tbe_0 + 2.0 * tr_zzz_yzzzz[i] * tke_0;
    }

    // Set up 300-315 components of targeted buffer : FG

    auto tr_0_0_z_xxx_xxxx = pbuffer.data(idx_op_geom_010_fg + 300);

    auto tr_0_0_z_xxx_xxxy = pbuffer.data(idx_op_geom_010_fg + 301);

    auto tr_0_0_z_xxx_xxxz = pbuffer.data(idx_op_geom_010_fg + 302);

    auto tr_0_0_z_xxx_xxyy = pbuffer.data(idx_op_geom_010_fg + 303);

    auto tr_0_0_z_xxx_xxyz = pbuffer.data(idx_op_geom_010_fg + 304);

    auto tr_0_0_z_xxx_xxzz = pbuffer.data(idx_op_geom_010_fg + 305);

    auto tr_0_0_z_xxx_xyyy = pbuffer.data(idx_op_geom_010_fg + 306);

    auto tr_0_0_z_xxx_xyyz = pbuffer.data(idx_op_geom_010_fg + 307);

    auto tr_0_0_z_xxx_xyzz = pbuffer.data(idx_op_geom_010_fg + 308);

    auto tr_0_0_z_xxx_xzzz = pbuffer.data(idx_op_geom_010_fg + 309);

    auto tr_0_0_z_xxx_yyyy = pbuffer.data(idx_op_geom_010_fg + 310);

    auto tr_0_0_z_xxx_yyyz = pbuffer.data(idx_op_geom_010_fg + 311);

    auto tr_0_0_z_xxx_yyzz = pbuffer.data(idx_op_geom_010_fg + 312);

    auto tr_0_0_z_xxx_yzzz = pbuffer.data(idx_op_geom_010_fg + 313);

    auto tr_0_0_z_xxx_zzzz = pbuffer.data(idx_op_geom_010_fg + 314);

    #pragma omp simd aligned(tr_0_0_z_xxx_xxxx, tr_0_0_z_xxx_xxxy, tr_0_0_z_xxx_xxxz, tr_0_0_z_xxx_xxyy, tr_0_0_z_xxx_xxyz, tr_0_0_z_xxx_xxzz, tr_0_0_z_xxx_xyyy, tr_0_0_z_xxx_xyyz, tr_0_0_z_xxx_xyzz, tr_0_0_z_xxx_xzzz, tr_0_0_z_xxx_yyyy, tr_0_0_z_xxx_yyyz, tr_0_0_z_xxx_yyzz, tr_0_0_z_xxx_yzzz, tr_0_0_z_xxx_zzzz, tr_xxx_xxx, tr_xxx_xxxxz, tr_xxx_xxxyz, tr_xxx_xxxzz, tr_xxx_xxy, tr_xxx_xxyyz, tr_xxx_xxyzz, tr_xxx_xxz, tr_xxx_xxzzz, tr_xxx_xyy, tr_xxx_xyyyz, tr_xxx_xyyzz, tr_xxx_xyz, tr_xxx_xyzzz, tr_xxx_xzz, tr_xxx_xzzzz, tr_xxx_yyy, tr_xxx_yyyyz, tr_xxx_yyyzz, tr_xxx_yyz, tr_xxx_yyzzz, tr_xxx_yzz, tr_xxx_yzzzz, tr_xxx_zzz, tr_xxx_zzzzz, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xzzz, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yzzz, tr_xxxz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxx_xxxx[i] = 2.0 * tr_xxxz_xxxx[i] * tbe_0 + 2.0 * tr_xxx_xxxxz[i] * tke_0;

        tr_0_0_z_xxx_xxxy[i] = 2.0 * tr_xxxz_xxxy[i] * tbe_0 + 2.0 * tr_xxx_xxxyz[i] * tke_0;

        tr_0_0_z_xxx_xxxz[i] = 2.0 * tr_xxxz_xxxz[i] * tbe_0 + 2.0 * tr_xxx_xxxzz[i] * tke_0 - tr_xxx_xxx[i];

        tr_0_0_z_xxx_xxyy[i] = 2.0 * tr_xxxz_xxyy[i] * tbe_0 + 2.0 * tr_xxx_xxyyz[i] * tke_0;

        tr_0_0_z_xxx_xxyz[i] = 2.0 * tr_xxxz_xxyz[i] * tbe_0 + 2.0 * tr_xxx_xxyzz[i] * tke_0 - tr_xxx_xxy[i];

        tr_0_0_z_xxx_xxzz[i] = 2.0 * tr_xxxz_xxzz[i] * tbe_0 + 2.0 * tr_xxx_xxzzz[i] * tke_0 - 2.0 * tr_xxx_xxz[i];

        tr_0_0_z_xxx_xyyy[i] = 2.0 * tr_xxxz_xyyy[i] * tbe_0 + 2.0 * tr_xxx_xyyyz[i] * tke_0;

        tr_0_0_z_xxx_xyyz[i] = 2.0 * tr_xxxz_xyyz[i] * tbe_0 + 2.0 * tr_xxx_xyyzz[i] * tke_0 - tr_xxx_xyy[i];

        tr_0_0_z_xxx_xyzz[i] = 2.0 * tr_xxxz_xyzz[i] * tbe_0 + 2.0 * tr_xxx_xyzzz[i] * tke_0 - 2.0 * tr_xxx_xyz[i];

        tr_0_0_z_xxx_xzzz[i] = 2.0 * tr_xxxz_xzzz[i] * tbe_0 + 2.0 * tr_xxx_xzzzz[i] * tke_0 - 3.0 * tr_xxx_xzz[i];

        tr_0_0_z_xxx_yyyy[i] = 2.0 * tr_xxxz_yyyy[i] * tbe_0 + 2.0 * tr_xxx_yyyyz[i] * tke_0;

        tr_0_0_z_xxx_yyyz[i] = 2.0 * tr_xxxz_yyyz[i] * tbe_0 + 2.0 * tr_xxx_yyyzz[i] * tke_0 - tr_xxx_yyy[i];

        tr_0_0_z_xxx_yyzz[i] = 2.0 * tr_xxxz_yyzz[i] * tbe_0 + 2.0 * tr_xxx_yyzzz[i] * tke_0 - 2.0 * tr_xxx_yyz[i];

        tr_0_0_z_xxx_yzzz[i] = 2.0 * tr_xxxz_yzzz[i] * tbe_0 + 2.0 * tr_xxx_yzzzz[i] * tke_0 - 3.0 * tr_xxx_yzz[i];

        tr_0_0_z_xxx_zzzz[i] = 2.0 * tr_xxxz_zzzz[i] * tbe_0 + 2.0 * tr_xxx_zzzzz[i] * tke_0 - 4.0 * tr_xxx_zzz[i];
    }

    // Set up 315-330 components of targeted buffer : FG

    auto tr_0_0_z_xxy_xxxx = pbuffer.data(idx_op_geom_010_fg + 315);

    auto tr_0_0_z_xxy_xxxy = pbuffer.data(idx_op_geom_010_fg + 316);

    auto tr_0_0_z_xxy_xxxz = pbuffer.data(idx_op_geom_010_fg + 317);

    auto tr_0_0_z_xxy_xxyy = pbuffer.data(idx_op_geom_010_fg + 318);

    auto tr_0_0_z_xxy_xxyz = pbuffer.data(idx_op_geom_010_fg + 319);

    auto tr_0_0_z_xxy_xxzz = pbuffer.data(idx_op_geom_010_fg + 320);

    auto tr_0_0_z_xxy_xyyy = pbuffer.data(idx_op_geom_010_fg + 321);

    auto tr_0_0_z_xxy_xyyz = pbuffer.data(idx_op_geom_010_fg + 322);

    auto tr_0_0_z_xxy_xyzz = pbuffer.data(idx_op_geom_010_fg + 323);

    auto tr_0_0_z_xxy_xzzz = pbuffer.data(idx_op_geom_010_fg + 324);

    auto tr_0_0_z_xxy_yyyy = pbuffer.data(idx_op_geom_010_fg + 325);

    auto tr_0_0_z_xxy_yyyz = pbuffer.data(idx_op_geom_010_fg + 326);

    auto tr_0_0_z_xxy_yyzz = pbuffer.data(idx_op_geom_010_fg + 327);

    auto tr_0_0_z_xxy_yzzz = pbuffer.data(idx_op_geom_010_fg + 328);

    auto tr_0_0_z_xxy_zzzz = pbuffer.data(idx_op_geom_010_fg + 329);

    #pragma omp simd aligned(tr_0_0_z_xxy_xxxx, tr_0_0_z_xxy_xxxy, tr_0_0_z_xxy_xxxz, tr_0_0_z_xxy_xxyy, tr_0_0_z_xxy_xxyz, tr_0_0_z_xxy_xxzz, tr_0_0_z_xxy_xyyy, tr_0_0_z_xxy_xyyz, tr_0_0_z_xxy_xyzz, tr_0_0_z_xxy_xzzz, tr_0_0_z_xxy_yyyy, tr_0_0_z_xxy_yyyz, tr_0_0_z_xxy_yyzz, tr_0_0_z_xxy_yzzz, tr_0_0_z_xxy_zzzz, tr_xxy_xxx, tr_xxy_xxxxz, tr_xxy_xxxyz, tr_xxy_xxxzz, tr_xxy_xxy, tr_xxy_xxyyz, tr_xxy_xxyzz, tr_xxy_xxz, tr_xxy_xxzzz, tr_xxy_xyy, tr_xxy_xyyyz, tr_xxy_xyyzz, tr_xxy_xyz, tr_xxy_xyzzz, tr_xxy_xzz, tr_xxy_xzzzz, tr_xxy_yyy, tr_xxy_yyyyz, tr_xxy_yyyzz, tr_xxy_yyz, tr_xxy_yyzzz, tr_xxy_yzz, tr_xxy_yzzzz, tr_xxy_zzz, tr_xxy_zzzzz, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xzzz, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yzzz, tr_xxyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxy_xxxx[i] = 2.0 * tr_xxyz_xxxx[i] * tbe_0 + 2.0 * tr_xxy_xxxxz[i] * tke_0;

        tr_0_0_z_xxy_xxxy[i] = 2.0 * tr_xxyz_xxxy[i] * tbe_0 + 2.0 * tr_xxy_xxxyz[i] * tke_0;

        tr_0_0_z_xxy_xxxz[i] = 2.0 * tr_xxyz_xxxz[i] * tbe_0 + 2.0 * tr_xxy_xxxzz[i] * tke_0 - tr_xxy_xxx[i];

        tr_0_0_z_xxy_xxyy[i] = 2.0 * tr_xxyz_xxyy[i] * tbe_0 + 2.0 * tr_xxy_xxyyz[i] * tke_0;

        tr_0_0_z_xxy_xxyz[i] = 2.0 * tr_xxyz_xxyz[i] * tbe_0 + 2.0 * tr_xxy_xxyzz[i] * tke_0 - tr_xxy_xxy[i];

        tr_0_0_z_xxy_xxzz[i] = 2.0 * tr_xxyz_xxzz[i] * tbe_0 + 2.0 * tr_xxy_xxzzz[i] * tke_0 - 2.0 * tr_xxy_xxz[i];

        tr_0_0_z_xxy_xyyy[i] = 2.0 * tr_xxyz_xyyy[i] * tbe_0 + 2.0 * tr_xxy_xyyyz[i] * tke_0;

        tr_0_0_z_xxy_xyyz[i] = 2.0 * tr_xxyz_xyyz[i] * tbe_0 + 2.0 * tr_xxy_xyyzz[i] * tke_0 - tr_xxy_xyy[i];

        tr_0_0_z_xxy_xyzz[i] = 2.0 * tr_xxyz_xyzz[i] * tbe_0 + 2.0 * tr_xxy_xyzzz[i] * tke_0 - 2.0 * tr_xxy_xyz[i];

        tr_0_0_z_xxy_xzzz[i] = 2.0 * tr_xxyz_xzzz[i] * tbe_0 + 2.0 * tr_xxy_xzzzz[i] * tke_0 - 3.0 * tr_xxy_xzz[i];

        tr_0_0_z_xxy_yyyy[i] = 2.0 * tr_xxyz_yyyy[i] * tbe_0 + 2.0 * tr_xxy_yyyyz[i] * tke_0;

        tr_0_0_z_xxy_yyyz[i] = 2.0 * tr_xxyz_yyyz[i] * tbe_0 + 2.0 * tr_xxy_yyyzz[i] * tke_0 - tr_xxy_yyy[i];

        tr_0_0_z_xxy_yyzz[i] = 2.0 * tr_xxyz_yyzz[i] * tbe_0 + 2.0 * tr_xxy_yyzzz[i] * tke_0 - 2.0 * tr_xxy_yyz[i];

        tr_0_0_z_xxy_yzzz[i] = 2.0 * tr_xxyz_yzzz[i] * tbe_0 + 2.0 * tr_xxy_yzzzz[i] * tke_0 - 3.0 * tr_xxy_yzz[i];

        tr_0_0_z_xxy_zzzz[i] = 2.0 * tr_xxyz_zzzz[i] * tbe_0 + 2.0 * tr_xxy_zzzzz[i] * tke_0 - 4.0 * tr_xxy_zzz[i];
    }

    // Set up 330-345 components of targeted buffer : FG

    auto tr_0_0_z_xxz_xxxx = pbuffer.data(idx_op_geom_010_fg + 330);

    auto tr_0_0_z_xxz_xxxy = pbuffer.data(idx_op_geom_010_fg + 331);

    auto tr_0_0_z_xxz_xxxz = pbuffer.data(idx_op_geom_010_fg + 332);

    auto tr_0_0_z_xxz_xxyy = pbuffer.data(idx_op_geom_010_fg + 333);

    auto tr_0_0_z_xxz_xxyz = pbuffer.data(idx_op_geom_010_fg + 334);

    auto tr_0_0_z_xxz_xxzz = pbuffer.data(idx_op_geom_010_fg + 335);

    auto tr_0_0_z_xxz_xyyy = pbuffer.data(idx_op_geom_010_fg + 336);

    auto tr_0_0_z_xxz_xyyz = pbuffer.data(idx_op_geom_010_fg + 337);

    auto tr_0_0_z_xxz_xyzz = pbuffer.data(idx_op_geom_010_fg + 338);

    auto tr_0_0_z_xxz_xzzz = pbuffer.data(idx_op_geom_010_fg + 339);

    auto tr_0_0_z_xxz_yyyy = pbuffer.data(idx_op_geom_010_fg + 340);

    auto tr_0_0_z_xxz_yyyz = pbuffer.data(idx_op_geom_010_fg + 341);

    auto tr_0_0_z_xxz_yyzz = pbuffer.data(idx_op_geom_010_fg + 342);

    auto tr_0_0_z_xxz_yzzz = pbuffer.data(idx_op_geom_010_fg + 343);

    auto tr_0_0_z_xxz_zzzz = pbuffer.data(idx_op_geom_010_fg + 344);

    #pragma omp simd aligned(tr_0_0_z_xxz_xxxx, tr_0_0_z_xxz_xxxy, tr_0_0_z_xxz_xxxz, tr_0_0_z_xxz_xxyy, tr_0_0_z_xxz_xxyz, tr_0_0_z_xxz_xxzz, tr_0_0_z_xxz_xyyy, tr_0_0_z_xxz_xyyz, tr_0_0_z_xxz_xyzz, tr_0_0_z_xxz_xzzz, tr_0_0_z_xxz_yyyy, tr_0_0_z_xxz_yyyz, tr_0_0_z_xxz_yyzz, tr_0_0_z_xxz_yzzz, tr_0_0_z_xxz_zzzz, tr_xx_xxxx, tr_xx_xxxy, tr_xx_xxxz, tr_xx_xxyy, tr_xx_xxyz, tr_xx_xxzz, tr_xx_xyyy, tr_xx_xyyz, tr_xx_xyzz, tr_xx_xzzz, tr_xx_yyyy, tr_xx_yyyz, tr_xx_yyzz, tr_xx_yzzz, tr_xx_zzzz, tr_xxz_xxx, tr_xxz_xxxxz, tr_xxz_xxxyz, tr_xxz_xxxzz, tr_xxz_xxy, tr_xxz_xxyyz, tr_xxz_xxyzz, tr_xxz_xxz, tr_xxz_xxzzz, tr_xxz_xyy, tr_xxz_xyyyz, tr_xxz_xyyzz, tr_xxz_xyz, tr_xxz_xyzzz, tr_xxz_xzz, tr_xxz_xzzzz, tr_xxz_yyy, tr_xxz_yyyyz, tr_xxz_yyyzz, tr_xxz_yyz, tr_xxz_yyzzz, tr_xxz_yzz, tr_xxz_yzzzz, tr_xxz_zzz, tr_xxz_zzzzz, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xzzz, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yzzz, tr_xxzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxz_xxxx[i] = 2.0 * tr_xxzz_xxxx[i] * tbe_0 + 2.0 * tr_xxz_xxxxz[i] * tke_0 - tr_xx_xxxx[i];

        tr_0_0_z_xxz_xxxy[i] = 2.0 * tr_xxzz_xxxy[i] * tbe_0 + 2.0 * tr_xxz_xxxyz[i] * tke_0 - tr_xx_xxxy[i];

        tr_0_0_z_xxz_xxxz[i] = 2.0 * tr_xxzz_xxxz[i] * tbe_0 + 2.0 * tr_xxz_xxxzz[i] * tke_0 - tr_xx_xxxz[i] - tr_xxz_xxx[i];

        tr_0_0_z_xxz_xxyy[i] = 2.0 * tr_xxzz_xxyy[i] * tbe_0 + 2.0 * tr_xxz_xxyyz[i] * tke_0 - tr_xx_xxyy[i];

        tr_0_0_z_xxz_xxyz[i] = 2.0 * tr_xxzz_xxyz[i] * tbe_0 + 2.0 * tr_xxz_xxyzz[i] * tke_0 - tr_xx_xxyz[i] - tr_xxz_xxy[i];

        tr_0_0_z_xxz_xxzz[i] = 2.0 * tr_xxzz_xxzz[i] * tbe_0 + 2.0 * tr_xxz_xxzzz[i] * tke_0 - tr_xx_xxzz[i] - 2.0 * tr_xxz_xxz[i];

        tr_0_0_z_xxz_xyyy[i] = 2.0 * tr_xxzz_xyyy[i] * tbe_0 + 2.0 * tr_xxz_xyyyz[i] * tke_0 - tr_xx_xyyy[i];

        tr_0_0_z_xxz_xyyz[i] = 2.0 * tr_xxzz_xyyz[i] * tbe_0 + 2.0 * tr_xxz_xyyzz[i] * tke_0 - tr_xx_xyyz[i] - tr_xxz_xyy[i];

        tr_0_0_z_xxz_xyzz[i] = 2.0 * tr_xxzz_xyzz[i] * tbe_0 + 2.0 * tr_xxz_xyzzz[i] * tke_0 - tr_xx_xyzz[i] - 2.0 * tr_xxz_xyz[i];

        tr_0_0_z_xxz_xzzz[i] = 2.0 * tr_xxzz_xzzz[i] * tbe_0 + 2.0 * tr_xxz_xzzzz[i] * tke_0 - tr_xx_xzzz[i] - 3.0 * tr_xxz_xzz[i];

        tr_0_0_z_xxz_yyyy[i] = 2.0 * tr_xxzz_yyyy[i] * tbe_0 + 2.0 * tr_xxz_yyyyz[i] * tke_0 - tr_xx_yyyy[i];

        tr_0_0_z_xxz_yyyz[i] = 2.0 * tr_xxzz_yyyz[i] * tbe_0 + 2.0 * tr_xxz_yyyzz[i] * tke_0 - tr_xx_yyyz[i] - tr_xxz_yyy[i];

        tr_0_0_z_xxz_yyzz[i] = 2.0 * tr_xxzz_yyzz[i] * tbe_0 + 2.0 * tr_xxz_yyzzz[i] * tke_0 - tr_xx_yyzz[i] - 2.0 * tr_xxz_yyz[i];

        tr_0_0_z_xxz_yzzz[i] = 2.0 * tr_xxzz_yzzz[i] * tbe_0 + 2.0 * tr_xxz_yzzzz[i] * tke_0 - tr_xx_yzzz[i] - 3.0 * tr_xxz_yzz[i];

        tr_0_0_z_xxz_zzzz[i] = 2.0 * tr_xxzz_zzzz[i] * tbe_0 + 2.0 * tr_xxz_zzzzz[i] * tke_0 - tr_xx_zzzz[i] - 4.0 * tr_xxz_zzz[i];
    }

    // Set up 345-360 components of targeted buffer : FG

    auto tr_0_0_z_xyy_xxxx = pbuffer.data(idx_op_geom_010_fg + 345);

    auto tr_0_0_z_xyy_xxxy = pbuffer.data(idx_op_geom_010_fg + 346);

    auto tr_0_0_z_xyy_xxxz = pbuffer.data(idx_op_geom_010_fg + 347);

    auto tr_0_0_z_xyy_xxyy = pbuffer.data(idx_op_geom_010_fg + 348);

    auto tr_0_0_z_xyy_xxyz = pbuffer.data(idx_op_geom_010_fg + 349);

    auto tr_0_0_z_xyy_xxzz = pbuffer.data(idx_op_geom_010_fg + 350);

    auto tr_0_0_z_xyy_xyyy = pbuffer.data(idx_op_geom_010_fg + 351);

    auto tr_0_0_z_xyy_xyyz = pbuffer.data(idx_op_geom_010_fg + 352);

    auto tr_0_0_z_xyy_xyzz = pbuffer.data(idx_op_geom_010_fg + 353);

    auto tr_0_0_z_xyy_xzzz = pbuffer.data(idx_op_geom_010_fg + 354);

    auto tr_0_0_z_xyy_yyyy = pbuffer.data(idx_op_geom_010_fg + 355);

    auto tr_0_0_z_xyy_yyyz = pbuffer.data(idx_op_geom_010_fg + 356);

    auto tr_0_0_z_xyy_yyzz = pbuffer.data(idx_op_geom_010_fg + 357);

    auto tr_0_0_z_xyy_yzzz = pbuffer.data(idx_op_geom_010_fg + 358);

    auto tr_0_0_z_xyy_zzzz = pbuffer.data(idx_op_geom_010_fg + 359);

    #pragma omp simd aligned(tr_0_0_z_xyy_xxxx, tr_0_0_z_xyy_xxxy, tr_0_0_z_xyy_xxxz, tr_0_0_z_xyy_xxyy, tr_0_0_z_xyy_xxyz, tr_0_0_z_xyy_xxzz, tr_0_0_z_xyy_xyyy, tr_0_0_z_xyy_xyyz, tr_0_0_z_xyy_xyzz, tr_0_0_z_xyy_xzzz, tr_0_0_z_xyy_yyyy, tr_0_0_z_xyy_yyyz, tr_0_0_z_xyy_yyzz, tr_0_0_z_xyy_yzzz, tr_0_0_z_xyy_zzzz, tr_xyy_xxx, tr_xyy_xxxxz, tr_xyy_xxxyz, tr_xyy_xxxzz, tr_xyy_xxy, tr_xyy_xxyyz, tr_xyy_xxyzz, tr_xyy_xxz, tr_xyy_xxzzz, tr_xyy_xyy, tr_xyy_xyyyz, tr_xyy_xyyzz, tr_xyy_xyz, tr_xyy_xyzzz, tr_xyy_xzz, tr_xyy_xzzzz, tr_xyy_yyy, tr_xyy_yyyyz, tr_xyy_yyyzz, tr_xyy_yyz, tr_xyy_yyzzz, tr_xyy_yzz, tr_xyy_yzzzz, tr_xyy_zzz, tr_xyy_zzzzz, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xzzz, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yzzz, tr_xyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyy_xxxx[i] = 2.0 * tr_xyyz_xxxx[i] * tbe_0 + 2.0 * tr_xyy_xxxxz[i] * tke_0;

        tr_0_0_z_xyy_xxxy[i] = 2.0 * tr_xyyz_xxxy[i] * tbe_0 + 2.0 * tr_xyy_xxxyz[i] * tke_0;

        tr_0_0_z_xyy_xxxz[i] = 2.0 * tr_xyyz_xxxz[i] * tbe_0 + 2.0 * tr_xyy_xxxzz[i] * tke_0 - tr_xyy_xxx[i];

        tr_0_0_z_xyy_xxyy[i] = 2.0 * tr_xyyz_xxyy[i] * tbe_0 + 2.0 * tr_xyy_xxyyz[i] * tke_0;

        tr_0_0_z_xyy_xxyz[i] = 2.0 * tr_xyyz_xxyz[i] * tbe_0 + 2.0 * tr_xyy_xxyzz[i] * tke_0 - tr_xyy_xxy[i];

        tr_0_0_z_xyy_xxzz[i] = 2.0 * tr_xyyz_xxzz[i] * tbe_0 + 2.0 * tr_xyy_xxzzz[i] * tke_0 - 2.0 * tr_xyy_xxz[i];

        tr_0_0_z_xyy_xyyy[i] = 2.0 * tr_xyyz_xyyy[i] * tbe_0 + 2.0 * tr_xyy_xyyyz[i] * tke_0;

        tr_0_0_z_xyy_xyyz[i] = 2.0 * tr_xyyz_xyyz[i] * tbe_0 + 2.0 * tr_xyy_xyyzz[i] * tke_0 - tr_xyy_xyy[i];

        tr_0_0_z_xyy_xyzz[i] = 2.0 * tr_xyyz_xyzz[i] * tbe_0 + 2.0 * tr_xyy_xyzzz[i] * tke_0 - 2.0 * tr_xyy_xyz[i];

        tr_0_0_z_xyy_xzzz[i] = 2.0 * tr_xyyz_xzzz[i] * tbe_0 + 2.0 * tr_xyy_xzzzz[i] * tke_0 - 3.0 * tr_xyy_xzz[i];

        tr_0_0_z_xyy_yyyy[i] = 2.0 * tr_xyyz_yyyy[i] * tbe_0 + 2.0 * tr_xyy_yyyyz[i] * tke_0;

        tr_0_0_z_xyy_yyyz[i] = 2.0 * tr_xyyz_yyyz[i] * tbe_0 + 2.0 * tr_xyy_yyyzz[i] * tke_0 - tr_xyy_yyy[i];

        tr_0_0_z_xyy_yyzz[i] = 2.0 * tr_xyyz_yyzz[i] * tbe_0 + 2.0 * tr_xyy_yyzzz[i] * tke_0 - 2.0 * tr_xyy_yyz[i];

        tr_0_0_z_xyy_yzzz[i] = 2.0 * tr_xyyz_yzzz[i] * tbe_0 + 2.0 * tr_xyy_yzzzz[i] * tke_0 - 3.0 * tr_xyy_yzz[i];

        tr_0_0_z_xyy_zzzz[i] = 2.0 * tr_xyyz_zzzz[i] * tbe_0 + 2.0 * tr_xyy_zzzzz[i] * tke_0 - 4.0 * tr_xyy_zzz[i];
    }

    // Set up 360-375 components of targeted buffer : FG

    auto tr_0_0_z_xyz_xxxx = pbuffer.data(idx_op_geom_010_fg + 360);

    auto tr_0_0_z_xyz_xxxy = pbuffer.data(idx_op_geom_010_fg + 361);

    auto tr_0_0_z_xyz_xxxz = pbuffer.data(idx_op_geom_010_fg + 362);

    auto tr_0_0_z_xyz_xxyy = pbuffer.data(idx_op_geom_010_fg + 363);

    auto tr_0_0_z_xyz_xxyz = pbuffer.data(idx_op_geom_010_fg + 364);

    auto tr_0_0_z_xyz_xxzz = pbuffer.data(idx_op_geom_010_fg + 365);

    auto tr_0_0_z_xyz_xyyy = pbuffer.data(idx_op_geom_010_fg + 366);

    auto tr_0_0_z_xyz_xyyz = pbuffer.data(idx_op_geom_010_fg + 367);

    auto tr_0_0_z_xyz_xyzz = pbuffer.data(idx_op_geom_010_fg + 368);

    auto tr_0_0_z_xyz_xzzz = pbuffer.data(idx_op_geom_010_fg + 369);

    auto tr_0_0_z_xyz_yyyy = pbuffer.data(idx_op_geom_010_fg + 370);

    auto tr_0_0_z_xyz_yyyz = pbuffer.data(idx_op_geom_010_fg + 371);

    auto tr_0_0_z_xyz_yyzz = pbuffer.data(idx_op_geom_010_fg + 372);

    auto tr_0_0_z_xyz_yzzz = pbuffer.data(idx_op_geom_010_fg + 373);

    auto tr_0_0_z_xyz_zzzz = pbuffer.data(idx_op_geom_010_fg + 374);

    #pragma omp simd aligned(tr_0_0_z_xyz_xxxx, tr_0_0_z_xyz_xxxy, tr_0_0_z_xyz_xxxz, tr_0_0_z_xyz_xxyy, tr_0_0_z_xyz_xxyz, tr_0_0_z_xyz_xxzz, tr_0_0_z_xyz_xyyy, tr_0_0_z_xyz_xyyz, tr_0_0_z_xyz_xyzz, tr_0_0_z_xyz_xzzz, tr_0_0_z_xyz_yyyy, tr_0_0_z_xyz_yyyz, tr_0_0_z_xyz_yyzz, tr_0_0_z_xyz_yzzz, tr_0_0_z_xyz_zzzz, tr_xy_xxxx, tr_xy_xxxy, tr_xy_xxxz, tr_xy_xxyy, tr_xy_xxyz, tr_xy_xxzz, tr_xy_xyyy, tr_xy_xyyz, tr_xy_xyzz, tr_xy_xzzz, tr_xy_yyyy, tr_xy_yyyz, tr_xy_yyzz, tr_xy_yzzz, tr_xy_zzzz, tr_xyz_xxx, tr_xyz_xxxxz, tr_xyz_xxxyz, tr_xyz_xxxzz, tr_xyz_xxy, tr_xyz_xxyyz, tr_xyz_xxyzz, tr_xyz_xxz, tr_xyz_xxzzz, tr_xyz_xyy, tr_xyz_xyyyz, tr_xyz_xyyzz, tr_xyz_xyz, tr_xyz_xyzzz, tr_xyz_xzz, tr_xyz_xzzzz, tr_xyz_yyy, tr_xyz_yyyyz, tr_xyz_yyyzz, tr_xyz_yyz, tr_xyz_yyzzz, tr_xyz_yzz, tr_xyz_yzzzz, tr_xyz_zzz, tr_xyz_zzzzz, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xzzz, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yzzz, tr_xyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyz_xxxx[i] = 2.0 * tr_xyzz_xxxx[i] * tbe_0 + 2.0 * tr_xyz_xxxxz[i] * tke_0 - tr_xy_xxxx[i];

        tr_0_0_z_xyz_xxxy[i] = 2.0 * tr_xyzz_xxxy[i] * tbe_0 + 2.0 * tr_xyz_xxxyz[i] * tke_0 - tr_xy_xxxy[i];

        tr_0_0_z_xyz_xxxz[i] = 2.0 * tr_xyzz_xxxz[i] * tbe_0 + 2.0 * tr_xyz_xxxzz[i] * tke_0 - tr_xy_xxxz[i] - tr_xyz_xxx[i];

        tr_0_0_z_xyz_xxyy[i] = 2.0 * tr_xyzz_xxyy[i] * tbe_0 + 2.0 * tr_xyz_xxyyz[i] * tke_0 - tr_xy_xxyy[i];

        tr_0_0_z_xyz_xxyz[i] = 2.0 * tr_xyzz_xxyz[i] * tbe_0 + 2.0 * tr_xyz_xxyzz[i] * tke_0 - tr_xy_xxyz[i] - tr_xyz_xxy[i];

        tr_0_0_z_xyz_xxzz[i] = 2.0 * tr_xyzz_xxzz[i] * tbe_0 + 2.0 * tr_xyz_xxzzz[i] * tke_0 - tr_xy_xxzz[i] - 2.0 * tr_xyz_xxz[i];

        tr_0_0_z_xyz_xyyy[i] = 2.0 * tr_xyzz_xyyy[i] * tbe_0 + 2.0 * tr_xyz_xyyyz[i] * tke_0 - tr_xy_xyyy[i];

        tr_0_0_z_xyz_xyyz[i] = 2.0 * tr_xyzz_xyyz[i] * tbe_0 + 2.0 * tr_xyz_xyyzz[i] * tke_0 - tr_xy_xyyz[i] - tr_xyz_xyy[i];

        tr_0_0_z_xyz_xyzz[i] = 2.0 * tr_xyzz_xyzz[i] * tbe_0 + 2.0 * tr_xyz_xyzzz[i] * tke_0 - tr_xy_xyzz[i] - 2.0 * tr_xyz_xyz[i];

        tr_0_0_z_xyz_xzzz[i] = 2.0 * tr_xyzz_xzzz[i] * tbe_0 + 2.0 * tr_xyz_xzzzz[i] * tke_0 - tr_xy_xzzz[i] - 3.0 * tr_xyz_xzz[i];

        tr_0_0_z_xyz_yyyy[i] = 2.0 * tr_xyzz_yyyy[i] * tbe_0 + 2.0 * tr_xyz_yyyyz[i] * tke_0 - tr_xy_yyyy[i];

        tr_0_0_z_xyz_yyyz[i] = 2.0 * tr_xyzz_yyyz[i] * tbe_0 + 2.0 * tr_xyz_yyyzz[i] * tke_0 - tr_xy_yyyz[i] - tr_xyz_yyy[i];

        tr_0_0_z_xyz_yyzz[i] = 2.0 * tr_xyzz_yyzz[i] * tbe_0 + 2.0 * tr_xyz_yyzzz[i] * tke_0 - tr_xy_yyzz[i] - 2.0 * tr_xyz_yyz[i];

        tr_0_0_z_xyz_yzzz[i] = 2.0 * tr_xyzz_yzzz[i] * tbe_0 + 2.0 * tr_xyz_yzzzz[i] * tke_0 - tr_xy_yzzz[i] - 3.0 * tr_xyz_yzz[i];

        tr_0_0_z_xyz_zzzz[i] = 2.0 * tr_xyzz_zzzz[i] * tbe_0 + 2.0 * tr_xyz_zzzzz[i] * tke_0 - tr_xy_zzzz[i] - 4.0 * tr_xyz_zzz[i];
    }

    // Set up 375-390 components of targeted buffer : FG

    auto tr_0_0_z_xzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 375);

    auto tr_0_0_z_xzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 376);

    auto tr_0_0_z_xzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 377);

    auto tr_0_0_z_xzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 378);

    auto tr_0_0_z_xzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 379);

    auto tr_0_0_z_xzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 380);

    auto tr_0_0_z_xzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 381);

    auto tr_0_0_z_xzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 382);

    auto tr_0_0_z_xzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 383);

    auto tr_0_0_z_xzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 384);

    auto tr_0_0_z_xzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 385);

    auto tr_0_0_z_xzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 386);

    auto tr_0_0_z_xzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 387);

    auto tr_0_0_z_xzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 388);

    auto tr_0_0_z_xzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 389);

    #pragma omp simd aligned(tr_0_0_z_xzz_xxxx, tr_0_0_z_xzz_xxxy, tr_0_0_z_xzz_xxxz, tr_0_0_z_xzz_xxyy, tr_0_0_z_xzz_xxyz, tr_0_0_z_xzz_xxzz, tr_0_0_z_xzz_xyyy, tr_0_0_z_xzz_xyyz, tr_0_0_z_xzz_xyzz, tr_0_0_z_xzz_xzzz, tr_0_0_z_xzz_yyyy, tr_0_0_z_xzz_yyyz, tr_0_0_z_xzz_yyzz, tr_0_0_z_xzz_yzzz, tr_0_0_z_xzz_zzzz, tr_xz_xxxx, tr_xz_xxxy, tr_xz_xxxz, tr_xz_xxyy, tr_xz_xxyz, tr_xz_xxzz, tr_xz_xyyy, tr_xz_xyyz, tr_xz_xyzz, tr_xz_xzzz, tr_xz_yyyy, tr_xz_yyyz, tr_xz_yyzz, tr_xz_yzzz, tr_xz_zzzz, tr_xzz_xxx, tr_xzz_xxxxz, tr_xzz_xxxyz, tr_xzz_xxxzz, tr_xzz_xxy, tr_xzz_xxyyz, tr_xzz_xxyzz, tr_xzz_xxz, tr_xzz_xxzzz, tr_xzz_xyy, tr_xzz_xyyyz, tr_xzz_xyyzz, tr_xzz_xyz, tr_xzz_xyzzz, tr_xzz_xzz, tr_xzz_xzzzz, tr_xzz_yyy, tr_xzz_yyyyz, tr_xzz_yyyzz, tr_xzz_yyz, tr_xzz_yyzzz, tr_xzz_yzz, tr_xzz_yzzzz, tr_xzz_zzz, tr_xzz_zzzzz, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xzzz, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yzzz, tr_xzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzz_xxxx[i] = 2.0 * tr_xzzz_xxxx[i] * tbe_0 + 2.0 * tr_xzz_xxxxz[i] * tke_0 - 2.0 * tr_xz_xxxx[i];

        tr_0_0_z_xzz_xxxy[i] = 2.0 * tr_xzzz_xxxy[i] * tbe_0 + 2.0 * tr_xzz_xxxyz[i] * tke_0 - 2.0 * tr_xz_xxxy[i];

        tr_0_0_z_xzz_xxxz[i] = 2.0 * tr_xzzz_xxxz[i] * tbe_0 + 2.0 * tr_xzz_xxxzz[i] * tke_0 - 2.0 * tr_xz_xxxz[i] - tr_xzz_xxx[i];

        tr_0_0_z_xzz_xxyy[i] = 2.0 * tr_xzzz_xxyy[i] * tbe_0 + 2.0 * tr_xzz_xxyyz[i] * tke_0 - 2.0 * tr_xz_xxyy[i];

        tr_0_0_z_xzz_xxyz[i] = 2.0 * tr_xzzz_xxyz[i] * tbe_0 + 2.0 * tr_xzz_xxyzz[i] * tke_0 - 2.0 * tr_xz_xxyz[i] - tr_xzz_xxy[i];

        tr_0_0_z_xzz_xxzz[i] = 2.0 * tr_xzzz_xxzz[i] * tbe_0 + 2.0 * tr_xzz_xxzzz[i] * tke_0 - 2.0 * tr_xz_xxzz[i] - 2.0 * tr_xzz_xxz[i];

        tr_0_0_z_xzz_xyyy[i] = 2.0 * tr_xzzz_xyyy[i] * tbe_0 + 2.0 * tr_xzz_xyyyz[i] * tke_0 - 2.0 * tr_xz_xyyy[i];

        tr_0_0_z_xzz_xyyz[i] = 2.0 * tr_xzzz_xyyz[i] * tbe_0 + 2.0 * tr_xzz_xyyzz[i] * tke_0 - 2.0 * tr_xz_xyyz[i] - tr_xzz_xyy[i];

        tr_0_0_z_xzz_xyzz[i] = 2.0 * tr_xzzz_xyzz[i] * tbe_0 + 2.0 * tr_xzz_xyzzz[i] * tke_0 - 2.0 * tr_xz_xyzz[i] - 2.0 * tr_xzz_xyz[i];

        tr_0_0_z_xzz_xzzz[i] = 2.0 * tr_xzzz_xzzz[i] * tbe_0 + 2.0 * tr_xzz_xzzzz[i] * tke_0 - 2.0 * tr_xz_xzzz[i] - 3.0 * tr_xzz_xzz[i];

        tr_0_0_z_xzz_yyyy[i] = 2.0 * tr_xzzz_yyyy[i] * tbe_0 + 2.0 * tr_xzz_yyyyz[i] * tke_0 - 2.0 * tr_xz_yyyy[i];

        tr_0_0_z_xzz_yyyz[i] = 2.0 * tr_xzzz_yyyz[i] * tbe_0 + 2.0 * tr_xzz_yyyzz[i] * tke_0 - 2.0 * tr_xz_yyyz[i] - tr_xzz_yyy[i];

        tr_0_0_z_xzz_yyzz[i] = 2.0 * tr_xzzz_yyzz[i] * tbe_0 + 2.0 * tr_xzz_yyzzz[i] * tke_0 - 2.0 * tr_xz_yyzz[i] - 2.0 * tr_xzz_yyz[i];

        tr_0_0_z_xzz_yzzz[i] = 2.0 * tr_xzzz_yzzz[i] * tbe_0 + 2.0 * tr_xzz_yzzzz[i] * tke_0 - 2.0 * tr_xz_yzzz[i] - 3.0 * tr_xzz_yzz[i];

        tr_0_0_z_xzz_zzzz[i] = 2.0 * tr_xzzz_zzzz[i] * tbe_0 + 2.0 * tr_xzz_zzzzz[i] * tke_0 - 2.0 * tr_xz_zzzz[i] - 4.0 * tr_xzz_zzz[i];
    }

    // Set up 390-405 components of targeted buffer : FG

    auto tr_0_0_z_yyy_xxxx = pbuffer.data(idx_op_geom_010_fg + 390);

    auto tr_0_0_z_yyy_xxxy = pbuffer.data(idx_op_geom_010_fg + 391);

    auto tr_0_0_z_yyy_xxxz = pbuffer.data(idx_op_geom_010_fg + 392);

    auto tr_0_0_z_yyy_xxyy = pbuffer.data(idx_op_geom_010_fg + 393);

    auto tr_0_0_z_yyy_xxyz = pbuffer.data(idx_op_geom_010_fg + 394);

    auto tr_0_0_z_yyy_xxzz = pbuffer.data(idx_op_geom_010_fg + 395);

    auto tr_0_0_z_yyy_xyyy = pbuffer.data(idx_op_geom_010_fg + 396);

    auto tr_0_0_z_yyy_xyyz = pbuffer.data(idx_op_geom_010_fg + 397);

    auto tr_0_0_z_yyy_xyzz = pbuffer.data(idx_op_geom_010_fg + 398);

    auto tr_0_0_z_yyy_xzzz = pbuffer.data(idx_op_geom_010_fg + 399);

    auto tr_0_0_z_yyy_yyyy = pbuffer.data(idx_op_geom_010_fg + 400);

    auto tr_0_0_z_yyy_yyyz = pbuffer.data(idx_op_geom_010_fg + 401);

    auto tr_0_0_z_yyy_yyzz = pbuffer.data(idx_op_geom_010_fg + 402);

    auto tr_0_0_z_yyy_yzzz = pbuffer.data(idx_op_geom_010_fg + 403);

    auto tr_0_0_z_yyy_zzzz = pbuffer.data(idx_op_geom_010_fg + 404);

    #pragma omp simd aligned(tr_0_0_z_yyy_xxxx, tr_0_0_z_yyy_xxxy, tr_0_0_z_yyy_xxxz, tr_0_0_z_yyy_xxyy, tr_0_0_z_yyy_xxyz, tr_0_0_z_yyy_xxzz, tr_0_0_z_yyy_xyyy, tr_0_0_z_yyy_xyyz, tr_0_0_z_yyy_xyzz, tr_0_0_z_yyy_xzzz, tr_0_0_z_yyy_yyyy, tr_0_0_z_yyy_yyyz, tr_0_0_z_yyy_yyzz, tr_0_0_z_yyy_yzzz, tr_0_0_z_yyy_zzzz, tr_yyy_xxx, tr_yyy_xxxxz, tr_yyy_xxxyz, tr_yyy_xxxzz, tr_yyy_xxy, tr_yyy_xxyyz, tr_yyy_xxyzz, tr_yyy_xxz, tr_yyy_xxzzz, tr_yyy_xyy, tr_yyy_xyyyz, tr_yyy_xyyzz, tr_yyy_xyz, tr_yyy_xyzzz, tr_yyy_xzz, tr_yyy_xzzzz, tr_yyy_yyy, tr_yyy_yyyyz, tr_yyy_yyyzz, tr_yyy_yyz, tr_yyy_yyzzz, tr_yyy_yzz, tr_yyy_yzzzz, tr_yyy_zzz, tr_yyy_zzzzz, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xzzz, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yzzz, tr_yyyz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyy_xxxx[i] = 2.0 * tr_yyyz_xxxx[i] * tbe_0 + 2.0 * tr_yyy_xxxxz[i] * tke_0;

        tr_0_0_z_yyy_xxxy[i] = 2.0 * tr_yyyz_xxxy[i] * tbe_0 + 2.0 * tr_yyy_xxxyz[i] * tke_0;

        tr_0_0_z_yyy_xxxz[i] = 2.0 * tr_yyyz_xxxz[i] * tbe_0 + 2.0 * tr_yyy_xxxzz[i] * tke_0 - tr_yyy_xxx[i];

        tr_0_0_z_yyy_xxyy[i] = 2.0 * tr_yyyz_xxyy[i] * tbe_0 + 2.0 * tr_yyy_xxyyz[i] * tke_0;

        tr_0_0_z_yyy_xxyz[i] = 2.0 * tr_yyyz_xxyz[i] * tbe_0 + 2.0 * tr_yyy_xxyzz[i] * tke_0 - tr_yyy_xxy[i];

        tr_0_0_z_yyy_xxzz[i] = 2.0 * tr_yyyz_xxzz[i] * tbe_0 + 2.0 * tr_yyy_xxzzz[i] * tke_0 - 2.0 * tr_yyy_xxz[i];

        tr_0_0_z_yyy_xyyy[i] = 2.0 * tr_yyyz_xyyy[i] * tbe_0 + 2.0 * tr_yyy_xyyyz[i] * tke_0;

        tr_0_0_z_yyy_xyyz[i] = 2.0 * tr_yyyz_xyyz[i] * tbe_0 + 2.0 * tr_yyy_xyyzz[i] * tke_0 - tr_yyy_xyy[i];

        tr_0_0_z_yyy_xyzz[i] = 2.0 * tr_yyyz_xyzz[i] * tbe_0 + 2.0 * tr_yyy_xyzzz[i] * tke_0 - 2.0 * tr_yyy_xyz[i];

        tr_0_0_z_yyy_xzzz[i] = 2.0 * tr_yyyz_xzzz[i] * tbe_0 + 2.0 * tr_yyy_xzzzz[i] * tke_0 - 3.0 * tr_yyy_xzz[i];

        tr_0_0_z_yyy_yyyy[i] = 2.0 * tr_yyyz_yyyy[i] * tbe_0 + 2.0 * tr_yyy_yyyyz[i] * tke_0;

        tr_0_0_z_yyy_yyyz[i] = 2.0 * tr_yyyz_yyyz[i] * tbe_0 + 2.0 * tr_yyy_yyyzz[i] * tke_0 - tr_yyy_yyy[i];

        tr_0_0_z_yyy_yyzz[i] = 2.0 * tr_yyyz_yyzz[i] * tbe_0 + 2.0 * tr_yyy_yyzzz[i] * tke_0 - 2.0 * tr_yyy_yyz[i];

        tr_0_0_z_yyy_yzzz[i] = 2.0 * tr_yyyz_yzzz[i] * tbe_0 + 2.0 * tr_yyy_yzzzz[i] * tke_0 - 3.0 * tr_yyy_yzz[i];

        tr_0_0_z_yyy_zzzz[i] = 2.0 * tr_yyyz_zzzz[i] * tbe_0 + 2.0 * tr_yyy_zzzzz[i] * tke_0 - 4.0 * tr_yyy_zzz[i];
    }

    // Set up 405-420 components of targeted buffer : FG

    auto tr_0_0_z_yyz_xxxx = pbuffer.data(idx_op_geom_010_fg + 405);

    auto tr_0_0_z_yyz_xxxy = pbuffer.data(idx_op_geom_010_fg + 406);

    auto tr_0_0_z_yyz_xxxz = pbuffer.data(idx_op_geom_010_fg + 407);

    auto tr_0_0_z_yyz_xxyy = pbuffer.data(idx_op_geom_010_fg + 408);

    auto tr_0_0_z_yyz_xxyz = pbuffer.data(idx_op_geom_010_fg + 409);

    auto tr_0_0_z_yyz_xxzz = pbuffer.data(idx_op_geom_010_fg + 410);

    auto tr_0_0_z_yyz_xyyy = pbuffer.data(idx_op_geom_010_fg + 411);

    auto tr_0_0_z_yyz_xyyz = pbuffer.data(idx_op_geom_010_fg + 412);

    auto tr_0_0_z_yyz_xyzz = pbuffer.data(idx_op_geom_010_fg + 413);

    auto tr_0_0_z_yyz_xzzz = pbuffer.data(idx_op_geom_010_fg + 414);

    auto tr_0_0_z_yyz_yyyy = pbuffer.data(idx_op_geom_010_fg + 415);

    auto tr_0_0_z_yyz_yyyz = pbuffer.data(idx_op_geom_010_fg + 416);

    auto tr_0_0_z_yyz_yyzz = pbuffer.data(idx_op_geom_010_fg + 417);

    auto tr_0_0_z_yyz_yzzz = pbuffer.data(idx_op_geom_010_fg + 418);

    auto tr_0_0_z_yyz_zzzz = pbuffer.data(idx_op_geom_010_fg + 419);

    #pragma omp simd aligned(tr_0_0_z_yyz_xxxx, tr_0_0_z_yyz_xxxy, tr_0_0_z_yyz_xxxz, tr_0_0_z_yyz_xxyy, tr_0_0_z_yyz_xxyz, tr_0_0_z_yyz_xxzz, tr_0_0_z_yyz_xyyy, tr_0_0_z_yyz_xyyz, tr_0_0_z_yyz_xyzz, tr_0_0_z_yyz_xzzz, tr_0_0_z_yyz_yyyy, tr_0_0_z_yyz_yyyz, tr_0_0_z_yyz_yyzz, tr_0_0_z_yyz_yzzz, tr_0_0_z_yyz_zzzz, tr_yy_xxxx, tr_yy_xxxy, tr_yy_xxxz, tr_yy_xxyy, tr_yy_xxyz, tr_yy_xxzz, tr_yy_xyyy, tr_yy_xyyz, tr_yy_xyzz, tr_yy_xzzz, tr_yy_yyyy, tr_yy_yyyz, tr_yy_yyzz, tr_yy_yzzz, tr_yy_zzzz, tr_yyz_xxx, tr_yyz_xxxxz, tr_yyz_xxxyz, tr_yyz_xxxzz, tr_yyz_xxy, tr_yyz_xxyyz, tr_yyz_xxyzz, tr_yyz_xxz, tr_yyz_xxzzz, tr_yyz_xyy, tr_yyz_xyyyz, tr_yyz_xyyzz, tr_yyz_xyz, tr_yyz_xyzzz, tr_yyz_xzz, tr_yyz_xzzzz, tr_yyz_yyy, tr_yyz_yyyyz, tr_yyz_yyyzz, tr_yyz_yyz, tr_yyz_yyzzz, tr_yyz_yzz, tr_yyz_yzzzz, tr_yyz_zzz, tr_yyz_zzzzz, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xzzz, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yzzz, tr_yyzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyz_xxxx[i] = 2.0 * tr_yyzz_xxxx[i] * tbe_0 + 2.0 * tr_yyz_xxxxz[i] * tke_0 - tr_yy_xxxx[i];

        tr_0_0_z_yyz_xxxy[i] = 2.0 * tr_yyzz_xxxy[i] * tbe_0 + 2.0 * tr_yyz_xxxyz[i] * tke_0 - tr_yy_xxxy[i];

        tr_0_0_z_yyz_xxxz[i] = 2.0 * tr_yyzz_xxxz[i] * tbe_0 + 2.0 * tr_yyz_xxxzz[i] * tke_0 - tr_yy_xxxz[i] - tr_yyz_xxx[i];

        tr_0_0_z_yyz_xxyy[i] = 2.0 * tr_yyzz_xxyy[i] * tbe_0 + 2.0 * tr_yyz_xxyyz[i] * tke_0 - tr_yy_xxyy[i];

        tr_0_0_z_yyz_xxyz[i] = 2.0 * tr_yyzz_xxyz[i] * tbe_0 + 2.0 * tr_yyz_xxyzz[i] * tke_0 - tr_yy_xxyz[i] - tr_yyz_xxy[i];

        tr_0_0_z_yyz_xxzz[i] = 2.0 * tr_yyzz_xxzz[i] * tbe_0 + 2.0 * tr_yyz_xxzzz[i] * tke_0 - tr_yy_xxzz[i] - 2.0 * tr_yyz_xxz[i];

        tr_0_0_z_yyz_xyyy[i] = 2.0 * tr_yyzz_xyyy[i] * tbe_0 + 2.0 * tr_yyz_xyyyz[i] * tke_0 - tr_yy_xyyy[i];

        tr_0_0_z_yyz_xyyz[i] = 2.0 * tr_yyzz_xyyz[i] * tbe_0 + 2.0 * tr_yyz_xyyzz[i] * tke_0 - tr_yy_xyyz[i] - tr_yyz_xyy[i];

        tr_0_0_z_yyz_xyzz[i] = 2.0 * tr_yyzz_xyzz[i] * tbe_0 + 2.0 * tr_yyz_xyzzz[i] * tke_0 - tr_yy_xyzz[i] - 2.0 * tr_yyz_xyz[i];

        tr_0_0_z_yyz_xzzz[i] = 2.0 * tr_yyzz_xzzz[i] * tbe_0 + 2.0 * tr_yyz_xzzzz[i] * tke_0 - tr_yy_xzzz[i] - 3.0 * tr_yyz_xzz[i];

        tr_0_0_z_yyz_yyyy[i] = 2.0 * tr_yyzz_yyyy[i] * tbe_0 + 2.0 * tr_yyz_yyyyz[i] * tke_0 - tr_yy_yyyy[i];

        tr_0_0_z_yyz_yyyz[i] = 2.0 * tr_yyzz_yyyz[i] * tbe_0 + 2.0 * tr_yyz_yyyzz[i] * tke_0 - tr_yy_yyyz[i] - tr_yyz_yyy[i];

        tr_0_0_z_yyz_yyzz[i] = 2.0 * tr_yyzz_yyzz[i] * tbe_0 + 2.0 * tr_yyz_yyzzz[i] * tke_0 - tr_yy_yyzz[i] - 2.0 * tr_yyz_yyz[i];

        tr_0_0_z_yyz_yzzz[i] = 2.0 * tr_yyzz_yzzz[i] * tbe_0 + 2.0 * tr_yyz_yzzzz[i] * tke_0 - tr_yy_yzzz[i] - 3.0 * tr_yyz_yzz[i];

        tr_0_0_z_yyz_zzzz[i] = 2.0 * tr_yyzz_zzzz[i] * tbe_0 + 2.0 * tr_yyz_zzzzz[i] * tke_0 - tr_yy_zzzz[i] - 4.0 * tr_yyz_zzz[i];
    }

    // Set up 420-435 components of targeted buffer : FG

    auto tr_0_0_z_yzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 420);

    auto tr_0_0_z_yzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 421);

    auto tr_0_0_z_yzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 422);

    auto tr_0_0_z_yzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 423);

    auto tr_0_0_z_yzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 424);

    auto tr_0_0_z_yzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 425);

    auto tr_0_0_z_yzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 426);

    auto tr_0_0_z_yzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 427);

    auto tr_0_0_z_yzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 428);

    auto tr_0_0_z_yzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 429);

    auto tr_0_0_z_yzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 430);

    auto tr_0_0_z_yzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 431);

    auto tr_0_0_z_yzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 432);

    auto tr_0_0_z_yzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 433);

    auto tr_0_0_z_yzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 434);

    #pragma omp simd aligned(tr_0_0_z_yzz_xxxx, tr_0_0_z_yzz_xxxy, tr_0_0_z_yzz_xxxz, tr_0_0_z_yzz_xxyy, tr_0_0_z_yzz_xxyz, tr_0_0_z_yzz_xxzz, tr_0_0_z_yzz_xyyy, tr_0_0_z_yzz_xyyz, tr_0_0_z_yzz_xyzz, tr_0_0_z_yzz_xzzz, tr_0_0_z_yzz_yyyy, tr_0_0_z_yzz_yyyz, tr_0_0_z_yzz_yyzz, tr_0_0_z_yzz_yzzz, tr_0_0_z_yzz_zzzz, tr_yz_xxxx, tr_yz_xxxy, tr_yz_xxxz, tr_yz_xxyy, tr_yz_xxyz, tr_yz_xxzz, tr_yz_xyyy, tr_yz_xyyz, tr_yz_xyzz, tr_yz_xzzz, tr_yz_yyyy, tr_yz_yyyz, tr_yz_yyzz, tr_yz_yzzz, tr_yz_zzzz, tr_yzz_xxx, tr_yzz_xxxxz, tr_yzz_xxxyz, tr_yzz_xxxzz, tr_yzz_xxy, tr_yzz_xxyyz, tr_yzz_xxyzz, tr_yzz_xxz, tr_yzz_xxzzz, tr_yzz_xyy, tr_yzz_xyyyz, tr_yzz_xyyzz, tr_yzz_xyz, tr_yzz_xyzzz, tr_yzz_xzz, tr_yzz_xzzzz, tr_yzz_yyy, tr_yzz_yyyyz, tr_yzz_yyyzz, tr_yzz_yyz, tr_yzz_yyzzz, tr_yzz_yzz, tr_yzz_yzzzz, tr_yzz_zzz, tr_yzz_zzzzz, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xzzz, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yzzz, tr_yzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzz_xxxx[i] = 2.0 * tr_yzzz_xxxx[i] * tbe_0 + 2.0 * tr_yzz_xxxxz[i] * tke_0 - 2.0 * tr_yz_xxxx[i];

        tr_0_0_z_yzz_xxxy[i] = 2.0 * tr_yzzz_xxxy[i] * tbe_0 + 2.0 * tr_yzz_xxxyz[i] * tke_0 - 2.0 * tr_yz_xxxy[i];

        tr_0_0_z_yzz_xxxz[i] = 2.0 * tr_yzzz_xxxz[i] * tbe_0 + 2.0 * tr_yzz_xxxzz[i] * tke_0 - 2.0 * tr_yz_xxxz[i] - tr_yzz_xxx[i];

        tr_0_0_z_yzz_xxyy[i] = 2.0 * tr_yzzz_xxyy[i] * tbe_0 + 2.0 * tr_yzz_xxyyz[i] * tke_0 - 2.0 * tr_yz_xxyy[i];

        tr_0_0_z_yzz_xxyz[i] = 2.0 * tr_yzzz_xxyz[i] * tbe_0 + 2.0 * tr_yzz_xxyzz[i] * tke_0 - 2.0 * tr_yz_xxyz[i] - tr_yzz_xxy[i];

        tr_0_0_z_yzz_xxzz[i] = 2.0 * tr_yzzz_xxzz[i] * tbe_0 + 2.0 * tr_yzz_xxzzz[i] * tke_0 - 2.0 * tr_yz_xxzz[i] - 2.0 * tr_yzz_xxz[i];

        tr_0_0_z_yzz_xyyy[i] = 2.0 * tr_yzzz_xyyy[i] * tbe_0 + 2.0 * tr_yzz_xyyyz[i] * tke_0 - 2.0 * tr_yz_xyyy[i];

        tr_0_0_z_yzz_xyyz[i] = 2.0 * tr_yzzz_xyyz[i] * tbe_0 + 2.0 * tr_yzz_xyyzz[i] * tke_0 - 2.0 * tr_yz_xyyz[i] - tr_yzz_xyy[i];

        tr_0_0_z_yzz_xyzz[i] = 2.0 * tr_yzzz_xyzz[i] * tbe_0 + 2.0 * tr_yzz_xyzzz[i] * tke_0 - 2.0 * tr_yz_xyzz[i] - 2.0 * tr_yzz_xyz[i];

        tr_0_0_z_yzz_xzzz[i] = 2.0 * tr_yzzz_xzzz[i] * tbe_0 + 2.0 * tr_yzz_xzzzz[i] * tke_0 - 2.0 * tr_yz_xzzz[i] - 3.0 * tr_yzz_xzz[i];

        tr_0_0_z_yzz_yyyy[i] = 2.0 * tr_yzzz_yyyy[i] * tbe_0 + 2.0 * tr_yzz_yyyyz[i] * tke_0 - 2.0 * tr_yz_yyyy[i];

        tr_0_0_z_yzz_yyyz[i] = 2.0 * tr_yzzz_yyyz[i] * tbe_0 + 2.0 * tr_yzz_yyyzz[i] * tke_0 - 2.0 * tr_yz_yyyz[i] - tr_yzz_yyy[i];

        tr_0_0_z_yzz_yyzz[i] = 2.0 * tr_yzzz_yyzz[i] * tbe_0 + 2.0 * tr_yzz_yyzzz[i] * tke_0 - 2.0 * tr_yz_yyzz[i] - 2.0 * tr_yzz_yyz[i];

        tr_0_0_z_yzz_yzzz[i] = 2.0 * tr_yzzz_yzzz[i] * tbe_0 + 2.0 * tr_yzz_yzzzz[i] * tke_0 - 2.0 * tr_yz_yzzz[i] - 3.0 * tr_yzz_yzz[i];

        tr_0_0_z_yzz_zzzz[i] = 2.0 * tr_yzzz_zzzz[i] * tbe_0 + 2.0 * tr_yzz_zzzzz[i] * tke_0 - 2.0 * tr_yz_zzzz[i] - 4.0 * tr_yzz_zzz[i];
    }

    // Set up 435-450 components of targeted buffer : FG

    auto tr_0_0_z_zzz_xxxx = pbuffer.data(idx_op_geom_010_fg + 435);

    auto tr_0_0_z_zzz_xxxy = pbuffer.data(idx_op_geom_010_fg + 436);

    auto tr_0_0_z_zzz_xxxz = pbuffer.data(idx_op_geom_010_fg + 437);

    auto tr_0_0_z_zzz_xxyy = pbuffer.data(idx_op_geom_010_fg + 438);

    auto tr_0_0_z_zzz_xxyz = pbuffer.data(idx_op_geom_010_fg + 439);

    auto tr_0_0_z_zzz_xxzz = pbuffer.data(idx_op_geom_010_fg + 440);

    auto tr_0_0_z_zzz_xyyy = pbuffer.data(idx_op_geom_010_fg + 441);

    auto tr_0_0_z_zzz_xyyz = pbuffer.data(idx_op_geom_010_fg + 442);

    auto tr_0_0_z_zzz_xyzz = pbuffer.data(idx_op_geom_010_fg + 443);

    auto tr_0_0_z_zzz_xzzz = pbuffer.data(idx_op_geom_010_fg + 444);

    auto tr_0_0_z_zzz_yyyy = pbuffer.data(idx_op_geom_010_fg + 445);

    auto tr_0_0_z_zzz_yyyz = pbuffer.data(idx_op_geom_010_fg + 446);

    auto tr_0_0_z_zzz_yyzz = pbuffer.data(idx_op_geom_010_fg + 447);

    auto tr_0_0_z_zzz_yzzz = pbuffer.data(idx_op_geom_010_fg + 448);

    auto tr_0_0_z_zzz_zzzz = pbuffer.data(idx_op_geom_010_fg + 449);

    #pragma omp simd aligned(tr_0_0_z_zzz_xxxx, tr_0_0_z_zzz_xxxy, tr_0_0_z_zzz_xxxz, tr_0_0_z_zzz_xxyy, tr_0_0_z_zzz_xxyz, tr_0_0_z_zzz_xxzz, tr_0_0_z_zzz_xyyy, tr_0_0_z_zzz_xyyz, tr_0_0_z_zzz_xyzz, tr_0_0_z_zzz_xzzz, tr_0_0_z_zzz_yyyy, tr_0_0_z_zzz_yyyz, tr_0_0_z_zzz_yyzz, tr_0_0_z_zzz_yzzz, tr_0_0_z_zzz_zzzz, tr_zz_xxxx, tr_zz_xxxy, tr_zz_xxxz, tr_zz_xxyy, tr_zz_xxyz, tr_zz_xxzz, tr_zz_xyyy, tr_zz_xyyz, tr_zz_xyzz, tr_zz_xzzz, tr_zz_yyyy, tr_zz_yyyz, tr_zz_yyzz, tr_zz_yzzz, tr_zz_zzzz, tr_zzz_xxx, tr_zzz_xxxxz, tr_zzz_xxxyz, tr_zzz_xxxzz, tr_zzz_xxy, tr_zzz_xxyyz, tr_zzz_xxyzz, tr_zzz_xxz, tr_zzz_xxzzz, tr_zzz_xyy, tr_zzz_xyyyz, tr_zzz_xyyzz, tr_zzz_xyz, tr_zzz_xyzzz, tr_zzz_xzz, tr_zzz_xzzzz, tr_zzz_yyy, tr_zzz_yyyyz, tr_zzz_yyyzz, tr_zzz_yyz, tr_zzz_yyzzz, tr_zzz_yzz, tr_zzz_yzzzz, tr_zzz_zzz, tr_zzz_zzzzz, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xzzz, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yzzz, tr_zzzz_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzz_xxxx[i] = 2.0 * tr_zzzz_xxxx[i] * tbe_0 + 2.0 * tr_zzz_xxxxz[i] * tke_0 - 3.0 * tr_zz_xxxx[i];

        tr_0_0_z_zzz_xxxy[i] = 2.0 * tr_zzzz_xxxy[i] * tbe_0 + 2.0 * tr_zzz_xxxyz[i] * tke_0 - 3.0 * tr_zz_xxxy[i];

        tr_0_0_z_zzz_xxxz[i] = 2.0 * tr_zzzz_xxxz[i] * tbe_0 + 2.0 * tr_zzz_xxxzz[i] * tke_0 - 3.0 * tr_zz_xxxz[i] - tr_zzz_xxx[i];

        tr_0_0_z_zzz_xxyy[i] = 2.0 * tr_zzzz_xxyy[i] * tbe_0 + 2.0 * tr_zzz_xxyyz[i] * tke_0 - 3.0 * tr_zz_xxyy[i];

        tr_0_0_z_zzz_xxyz[i] = 2.0 * tr_zzzz_xxyz[i] * tbe_0 + 2.0 * tr_zzz_xxyzz[i] * tke_0 - 3.0 * tr_zz_xxyz[i] - tr_zzz_xxy[i];

        tr_0_0_z_zzz_xxzz[i] = 2.0 * tr_zzzz_xxzz[i] * tbe_0 + 2.0 * tr_zzz_xxzzz[i] * tke_0 - 3.0 * tr_zz_xxzz[i] - 2.0 * tr_zzz_xxz[i];

        tr_0_0_z_zzz_xyyy[i] = 2.0 * tr_zzzz_xyyy[i] * tbe_0 + 2.0 * tr_zzz_xyyyz[i] * tke_0 - 3.0 * tr_zz_xyyy[i];

        tr_0_0_z_zzz_xyyz[i] = 2.0 * tr_zzzz_xyyz[i] * tbe_0 + 2.0 * tr_zzz_xyyzz[i] * tke_0 - 3.0 * tr_zz_xyyz[i] - tr_zzz_xyy[i];

        tr_0_0_z_zzz_xyzz[i] = 2.0 * tr_zzzz_xyzz[i] * tbe_0 + 2.0 * tr_zzz_xyzzz[i] * tke_0 - 3.0 * tr_zz_xyzz[i] - 2.0 * tr_zzz_xyz[i];

        tr_0_0_z_zzz_xzzz[i] = 2.0 * tr_zzzz_xzzz[i] * tbe_0 + 2.0 * tr_zzz_xzzzz[i] * tke_0 - 3.0 * tr_zz_xzzz[i] - 3.0 * tr_zzz_xzz[i];

        tr_0_0_z_zzz_yyyy[i] = 2.0 * tr_zzzz_yyyy[i] * tbe_0 + 2.0 * tr_zzz_yyyyz[i] * tke_0 - 3.0 * tr_zz_yyyy[i];

        tr_0_0_z_zzz_yyyz[i] = 2.0 * tr_zzzz_yyyz[i] * tbe_0 + 2.0 * tr_zzz_yyyzz[i] * tke_0 - 3.0 * tr_zz_yyyz[i] - tr_zzz_yyy[i];

        tr_0_0_z_zzz_yyzz[i] = 2.0 * tr_zzzz_yyzz[i] * tbe_0 + 2.0 * tr_zzz_yyzzz[i] * tke_0 - 3.0 * tr_zz_yyzz[i] - 2.0 * tr_zzz_yyz[i];

        tr_0_0_z_zzz_yzzz[i] = 2.0 * tr_zzzz_yzzz[i] * tbe_0 + 2.0 * tr_zzz_yzzzz[i] * tke_0 - 3.0 * tr_zz_yzzz[i] - 3.0 * tr_zzz_yzz[i];

        tr_0_0_z_zzz_zzzz[i] = 2.0 * tr_zzzz_zzzz[i] * tbe_0 + 2.0 * tr_zzz_zzzzz[i] * tke_0 - 3.0 * tr_zz_zzzz[i] - 4.0 * tr_zzz_zzz[i];
    }

}

} // t2cgeom namespace

