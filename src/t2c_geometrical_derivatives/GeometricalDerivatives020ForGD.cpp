#include "GeometricalDerivatives020ForGD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_020_gd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_020_gd,
                         const int idx_op_dd,
                         const int idx_op_fp,
                         const int idx_op_ff,
                         const int idx_op_gs,
                         const int idx_op_gd,
                         const int idx_op_gg,
                         const int idx_op_hp,
                         const int idx_op_hf,
                         const int idx_op_id,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up components of auxiliary buffer : ID

    auto tr_xxxxxx_xx = pbuffer.data(idx_op_id);

    auto tr_xxxxxx_xy = pbuffer.data(idx_op_id + 1);

    auto tr_xxxxxx_xz = pbuffer.data(idx_op_id + 2);

    auto tr_xxxxxx_yy = pbuffer.data(idx_op_id + 3);

    auto tr_xxxxxx_yz = pbuffer.data(idx_op_id + 4);

    auto tr_xxxxxx_zz = pbuffer.data(idx_op_id + 5);

    auto tr_xxxxxy_xx = pbuffer.data(idx_op_id + 6);

    auto tr_xxxxxy_xy = pbuffer.data(idx_op_id + 7);

    auto tr_xxxxxy_xz = pbuffer.data(idx_op_id + 8);

    auto tr_xxxxxy_yy = pbuffer.data(idx_op_id + 9);

    auto tr_xxxxxy_yz = pbuffer.data(idx_op_id + 10);

    auto tr_xxxxxy_zz = pbuffer.data(idx_op_id + 11);

    auto tr_xxxxxz_xx = pbuffer.data(idx_op_id + 12);

    auto tr_xxxxxz_xy = pbuffer.data(idx_op_id + 13);

    auto tr_xxxxxz_xz = pbuffer.data(idx_op_id + 14);

    auto tr_xxxxxz_yy = pbuffer.data(idx_op_id + 15);

    auto tr_xxxxxz_yz = pbuffer.data(idx_op_id + 16);

    auto tr_xxxxxz_zz = pbuffer.data(idx_op_id + 17);

    auto tr_xxxxyy_xx = pbuffer.data(idx_op_id + 18);

    auto tr_xxxxyy_xy = pbuffer.data(idx_op_id + 19);

    auto tr_xxxxyy_xz = pbuffer.data(idx_op_id + 20);

    auto tr_xxxxyy_yy = pbuffer.data(idx_op_id + 21);

    auto tr_xxxxyy_yz = pbuffer.data(idx_op_id + 22);

    auto tr_xxxxyy_zz = pbuffer.data(idx_op_id + 23);

    auto tr_xxxxyz_xx = pbuffer.data(idx_op_id + 24);

    auto tr_xxxxyz_xy = pbuffer.data(idx_op_id + 25);

    auto tr_xxxxyz_xz = pbuffer.data(idx_op_id + 26);

    auto tr_xxxxyz_yy = pbuffer.data(idx_op_id + 27);

    auto tr_xxxxyz_yz = pbuffer.data(idx_op_id + 28);

    auto tr_xxxxyz_zz = pbuffer.data(idx_op_id + 29);

    auto tr_xxxxzz_xx = pbuffer.data(idx_op_id + 30);

    auto tr_xxxxzz_xy = pbuffer.data(idx_op_id + 31);

    auto tr_xxxxzz_xz = pbuffer.data(idx_op_id + 32);

    auto tr_xxxxzz_yy = pbuffer.data(idx_op_id + 33);

    auto tr_xxxxzz_yz = pbuffer.data(idx_op_id + 34);

    auto tr_xxxxzz_zz = pbuffer.data(idx_op_id + 35);

    auto tr_xxxyyy_xx = pbuffer.data(idx_op_id + 36);

    auto tr_xxxyyy_xy = pbuffer.data(idx_op_id + 37);

    auto tr_xxxyyy_xz = pbuffer.data(idx_op_id + 38);

    auto tr_xxxyyy_yy = pbuffer.data(idx_op_id + 39);

    auto tr_xxxyyy_yz = pbuffer.data(idx_op_id + 40);

    auto tr_xxxyyy_zz = pbuffer.data(idx_op_id + 41);

    auto tr_xxxyyz_xx = pbuffer.data(idx_op_id + 42);

    auto tr_xxxyyz_xy = pbuffer.data(idx_op_id + 43);

    auto tr_xxxyyz_xz = pbuffer.data(idx_op_id + 44);

    auto tr_xxxyyz_yy = pbuffer.data(idx_op_id + 45);

    auto tr_xxxyyz_yz = pbuffer.data(idx_op_id + 46);

    auto tr_xxxyyz_zz = pbuffer.data(idx_op_id + 47);

    auto tr_xxxyzz_xx = pbuffer.data(idx_op_id + 48);

    auto tr_xxxyzz_xy = pbuffer.data(idx_op_id + 49);

    auto tr_xxxyzz_xz = pbuffer.data(idx_op_id + 50);

    auto tr_xxxyzz_yy = pbuffer.data(idx_op_id + 51);

    auto tr_xxxyzz_yz = pbuffer.data(idx_op_id + 52);

    auto tr_xxxyzz_zz = pbuffer.data(idx_op_id + 53);

    auto tr_xxxzzz_xx = pbuffer.data(idx_op_id + 54);

    auto tr_xxxzzz_xy = pbuffer.data(idx_op_id + 55);

    auto tr_xxxzzz_xz = pbuffer.data(idx_op_id + 56);

    auto tr_xxxzzz_yy = pbuffer.data(idx_op_id + 57);

    auto tr_xxxzzz_yz = pbuffer.data(idx_op_id + 58);

    auto tr_xxxzzz_zz = pbuffer.data(idx_op_id + 59);

    auto tr_xxyyyy_xx = pbuffer.data(idx_op_id + 60);

    auto tr_xxyyyy_xy = pbuffer.data(idx_op_id + 61);

    auto tr_xxyyyy_xz = pbuffer.data(idx_op_id + 62);

    auto tr_xxyyyy_yy = pbuffer.data(idx_op_id + 63);

    auto tr_xxyyyy_yz = pbuffer.data(idx_op_id + 64);

    auto tr_xxyyyy_zz = pbuffer.data(idx_op_id + 65);

    auto tr_xxyyyz_xx = pbuffer.data(idx_op_id + 66);

    auto tr_xxyyyz_xy = pbuffer.data(idx_op_id + 67);

    auto tr_xxyyyz_xz = pbuffer.data(idx_op_id + 68);

    auto tr_xxyyyz_yy = pbuffer.data(idx_op_id + 69);

    auto tr_xxyyyz_yz = pbuffer.data(idx_op_id + 70);

    auto tr_xxyyyz_zz = pbuffer.data(idx_op_id + 71);

    auto tr_xxyyzz_xx = pbuffer.data(idx_op_id + 72);

    auto tr_xxyyzz_xy = pbuffer.data(idx_op_id + 73);

    auto tr_xxyyzz_xz = pbuffer.data(idx_op_id + 74);

    auto tr_xxyyzz_yy = pbuffer.data(idx_op_id + 75);

    auto tr_xxyyzz_yz = pbuffer.data(idx_op_id + 76);

    auto tr_xxyyzz_zz = pbuffer.data(idx_op_id + 77);

    auto tr_xxyzzz_xx = pbuffer.data(idx_op_id + 78);

    auto tr_xxyzzz_xy = pbuffer.data(idx_op_id + 79);

    auto tr_xxyzzz_xz = pbuffer.data(idx_op_id + 80);

    auto tr_xxyzzz_yy = pbuffer.data(idx_op_id + 81);

    auto tr_xxyzzz_yz = pbuffer.data(idx_op_id + 82);

    auto tr_xxyzzz_zz = pbuffer.data(idx_op_id + 83);

    auto tr_xxzzzz_xx = pbuffer.data(idx_op_id + 84);

    auto tr_xxzzzz_xy = pbuffer.data(idx_op_id + 85);

    auto tr_xxzzzz_xz = pbuffer.data(idx_op_id + 86);

    auto tr_xxzzzz_yy = pbuffer.data(idx_op_id + 87);

    auto tr_xxzzzz_yz = pbuffer.data(idx_op_id + 88);

    auto tr_xxzzzz_zz = pbuffer.data(idx_op_id + 89);

    auto tr_xyyyyy_xx = pbuffer.data(idx_op_id + 90);

    auto tr_xyyyyy_xy = pbuffer.data(idx_op_id + 91);

    auto tr_xyyyyy_xz = pbuffer.data(idx_op_id + 92);

    auto tr_xyyyyy_yy = pbuffer.data(idx_op_id + 93);

    auto tr_xyyyyy_yz = pbuffer.data(idx_op_id + 94);

    auto tr_xyyyyy_zz = pbuffer.data(idx_op_id + 95);

    auto tr_xyyyyz_xx = pbuffer.data(idx_op_id + 96);

    auto tr_xyyyyz_xy = pbuffer.data(idx_op_id + 97);

    auto tr_xyyyyz_xz = pbuffer.data(idx_op_id + 98);

    auto tr_xyyyyz_yy = pbuffer.data(idx_op_id + 99);

    auto tr_xyyyyz_yz = pbuffer.data(idx_op_id + 100);

    auto tr_xyyyyz_zz = pbuffer.data(idx_op_id + 101);

    auto tr_xyyyzz_xx = pbuffer.data(idx_op_id + 102);

    auto tr_xyyyzz_xy = pbuffer.data(idx_op_id + 103);

    auto tr_xyyyzz_xz = pbuffer.data(idx_op_id + 104);

    auto tr_xyyyzz_yy = pbuffer.data(idx_op_id + 105);

    auto tr_xyyyzz_yz = pbuffer.data(idx_op_id + 106);

    auto tr_xyyyzz_zz = pbuffer.data(idx_op_id + 107);

    auto tr_xyyzzz_xx = pbuffer.data(idx_op_id + 108);

    auto tr_xyyzzz_xy = pbuffer.data(idx_op_id + 109);

    auto tr_xyyzzz_xz = pbuffer.data(idx_op_id + 110);

    auto tr_xyyzzz_yy = pbuffer.data(idx_op_id + 111);

    auto tr_xyyzzz_yz = pbuffer.data(idx_op_id + 112);

    auto tr_xyyzzz_zz = pbuffer.data(idx_op_id + 113);

    auto tr_xyzzzz_xx = pbuffer.data(idx_op_id + 114);

    auto tr_xyzzzz_xy = pbuffer.data(idx_op_id + 115);

    auto tr_xyzzzz_xz = pbuffer.data(idx_op_id + 116);

    auto tr_xyzzzz_yy = pbuffer.data(idx_op_id + 117);

    auto tr_xyzzzz_yz = pbuffer.data(idx_op_id + 118);

    auto tr_xyzzzz_zz = pbuffer.data(idx_op_id + 119);

    auto tr_xzzzzz_xx = pbuffer.data(idx_op_id + 120);

    auto tr_xzzzzz_xy = pbuffer.data(idx_op_id + 121);

    auto tr_xzzzzz_xz = pbuffer.data(idx_op_id + 122);

    auto tr_xzzzzz_yy = pbuffer.data(idx_op_id + 123);

    auto tr_xzzzzz_yz = pbuffer.data(idx_op_id + 124);

    auto tr_xzzzzz_zz = pbuffer.data(idx_op_id + 125);

    auto tr_yyyyyy_xx = pbuffer.data(idx_op_id + 126);

    auto tr_yyyyyy_xy = pbuffer.data(idx_op_id + 127);

    auto tr_yyyyyy_xz = pbuffer.data(idx_op_id + 128);

    auto tr_yyyyyy_yy = pbuffer.data(idx_op_id + 129);

    auto tr_yyyyyy_yz = pbuffer.data(idx_op_id + 130);

    auto tr_yyyyyy_zz = pbuffer.data(idx_op_id + 131);

    auto tr_yyyyyz_xx = pbuffer.data(idx_op_id + 132);

    auto tr_yyyyyz_xy = pbuffer.data(idx_op_id + 133);

    auto tr_yyyyyz_xz = pbuffer.data(idx_op_id + 134);

    auto tr_yyyyyz_yy = pbuffer.data(idx_op_id + 135);

    auto tr_yyyyyz_yz = pbuffer.data(idx_op_id + 136);

    auto tr_yyyyyz_zz = pbuffer.data(idx_op_id + 137);

    auto tr_yyyyzz_xx = pbuffer.data(idx_op_id + 138);

    auto tr_yyyyzz_xy = pbuffer.data(idx_op_id + 139);

    auto tr_yyyyzz_xz = pbuffer.data(idx_op_id + 140);

    auto tr_yyyyzz_yy = pbuffer.data(idx_op_id + 141);

    auto tr_yyyyzz_yz = pbuffer.data(idx_op_id + 142);

    auto tr_yyyyzz_zz = pbuffer.data(idx_op_id + 143);

    auto tr_yyyzzz_xx = pbuffer.data(idx_op_id + 144);

    auto tr_yyyzzz_xy = pbuffer.data(idx_op_id + 145);

    auto tr_yyyzzz_xz = pbuffer.data(idx_op_id + 146);

    auto tr_yyyzzz_yy = pbuffer.data(idx_op_id + 147);

    auto tr_yyyzzz_yz = pbuffer.data(idx_op_id + 148);

    auto tr_yyyzzz_zz = pbuffer.data(idx_op_id + 149);

    auto tr_yyzzzz_xx = pbuffer.data(idx_op_id + 150);

    auto tr_yyzzzz_xy = pbuffer.data(idx_op_id + 151);

    auto tr_yyzzzz_xz = pbuffer.data(idx_op_id + 152);

    auto tr_yyzzzz_yy = pbuffer.data(idx_op_id + 153);

    auto tr_yyzzzz_yz = pbuffer.data(idx_op_id + 154);

    auto tr_yyzzzz_zz = pbuffer.data(idx_op_id + 155);

    auto tr_yzzzzz_xx = pbuffer.data(idx_op_id + 156);

    auto tr_yzzzzz_xy = pbuffer.data(idx_op_id + 157);

    auto tr_yzzzzz_xz = pbuffer.data(idx_op_id + 158);

    auto tr_yzzzzz_yy = pbuffer.data(idx_op_id + 159);

    auto tr_yzzzzz_yz = pbuffer.data(idx_op_id + 160);

    auto tr_yzzzzz_zz = pbuffer.data(idx_op_id + 161);

    auto tr_zzzzzz_xx = pbuffer.data(idx_op_id + 162);

    auto tr_zzzzzz_xy = pbuffer.data(idx_op_id + 163);

    auto tr_zzzzzz_xz = pbuffer.data(idx_op_id + 164);

    auto tr_zzzzzz_yy = pbuffer.data(idx_op_id + 165);

    auto tr_zzzzzz_yz = pbuffer.data(idx_op_id + 166);

    auto tr_zzzzzz_zz = pbuffer.data(idx_op_id + 167);

    // Set up 0-6 components of targeted buffer : GD

    auto tr_0_0_xx_xxxx_xx = pbuffer.data(idx_op_geom_020_gd);

    auto tr_0_0_xx_xxxx_xy = pbuffer.data(idx_op_geom_020_gd + 1);

    auto tr_0_0_xx_xxxx_xz = pbuffer.data(idx_op_geom_020_gd + 2);

    auto tr_0_0_xx_xxxx_yy = pbuffer.data(idx_op_geom_020_gd + 3);

    auto tr_0_0_xx_xxxx_yz = pbuffer.data(idx_op_geom_020_gd + 4);

    auto tr_0_0_xx_xxxx_zz = pbuffer.data(idx_op_geom_020_gd + 5);

    #pragma omp simd aligned(tr_0_0_xx_xxxx_xx, tr_0_0_xx_xxxx_xy, tr_0_0_xx_xxxx_xz, tr_0_0_xx_xxxx_yy, tr_0_0_xx_xxxx_yz, tr_0_0_xx_xxxx_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_z, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxx_x, tr_xxxxx_xxx, tr_xxxxx_xxy, tr_xxxxx_xxz, tr_xxxxx_xyy, tr_xxxxx_xyz, tr_xxxxx_xzz, tr_xxxxx_y, tr_xxxxx_z, tr_xxxxxx_xx, tr_xxxxxx_xy, tr_xxxxxx_xz, tr_xxxxxx_yy, tr_xxxxxx_yz, tr_xxxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxx_xx[i] = 12.0 * tr_xx_xx[i] + 16.0 * tr_xxx_x[i] - 16.0 * tr_xxx_xxx[i] * tke_0 + 2.0 * tr_xxxx_0[i] - 18.0 * tr_xxxx_xx[i] * tbe_0 - 10.0 * tr_xxxx_xx[i] * tke_0 + 4.0 * tr_xxxx_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxxxx_x[i] * tbe_0 + 8.0 * tr_xxxxx_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xy[i] = 12.0 * tr_xx_xy[i] + 8.0 * tr_xxx_y[i] - 16.0 * tr_xxx_xxy[i] * tke_0 - 18.0 * tr_xxxx_xy[i] * tbe_0 - 6.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxx_y[i] * tbe_0 + 8.0 * tr_xxxxx_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_xz[i] = 12.0 * tr_xx_xz[i] + 8.0 * tr_xxx_z[i] - 16.0 * tr_xxx_xxz[i] * tke_0 - 18.0 * tr_xxxx_xz[i] * tbe_0 - 6.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxx_z[i] * tbe_0 + 8.0 * tr_xxxxx_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_yy[i] = 12.0 * tr_xx_yy[i] - 16.0 * tr_xxx_xyy[i] * tke_0 - 18.0 * tr_xxxx_yy[i] * tbe_0 - 2.0 * tr_xxxx_yy[i] * tke_0 + 4.0 * tr_xxxx_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_yz[i] = 12.0 * tr_xx_yz[i] - 16.0 * tr_xxx_xyz[i] * tke_0 - 18.0 * tr_xxxx_yz[i] * tbe_0 - 2.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxx_zz[i] = 12.0 * tr_xx_zz[i] - 16.0 * tr_xxx_xzz[i] * tke_0 - 18.0 * tr_xxxx_zz[i] * tbe_0 - 2.0 * tr_xxxx_zz[i] * tke_0 + 4.0 * tr_xxxx_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxx_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 6-12 components of targeted buffer : GD

    auto tr_0_0_xx_xxxy_xx = pbuffer.data(idx_op_geom_020_gd + 6);

    auto tr_0_0_xx_xxxy_xy = pbuffer.data(idx_op_geom_020_gd + 7);

    auto tr_0_0_xx_xxxy_xz = pbuffer.data(idx_op_geom_020_gd + 8);

    auto tr_0_0_xx_xxxy_yy = pbuffer.data(idx_op_geom_020_gd + 9);

    auto tr_0_0_xx_xxxy_yz = pbuffer.data(idx_op_geom_020_gd + 10);

    auto tr_0_0_xx_xxxy_zz = pbuffer.data(idx_op_geom_020_gd + 11);

    #pragma omp simd aligned(tr_0_0_xx_xxxy_xx, tr_0_0_xx_xxxy_xy, tr_0_0_xx_xxxy_xz, tr_0_0_xx_xxxy_yy, tr_0_0_xx_xxxy_yz, tr_0_0_xx_xxxy_zz, tr_xxxxxy_xx, tr_xxxxxy_xy, tr_xxxxxy_xz, tr_xxxxxy_yy, tr_xxxxxy_yz, tr_xxxxxy_zz, tr_xxxxy_x, tr_xxxxy_xxx, tr_xxxxy_xxy, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxy_xx[i] = 6.0 * tr_xy_xx[i] + 12.0 * tr_xxy_x[i] - 12.0 * tr_xxy_xxx[i] * tke_0 + 2.0 * tr_xxxy_0[i] - 14.0 * tr_xxxy_xx[i] * tbe_0 - 10.0 * tr_xxxy_xx[i] * tke_0 + 4.0 * tr_xxxy_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_x[i] * tbe_0 + 8.0 * tr_xxxxy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xy[i] = 6.0 * tr_xy_xy[i] + 6.0 * tr_xxy_y[i] - 12.0 * tr_xxy_xxy[i] * tke_0 - 14.0 * tr_xxxy_xy[i] * tbe_0 - 6.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_y[i] * tbe_0 + 8.0 * tr_xxxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_xz[i] = 6.0 * tr_xy_xz[i] + 6.0 * tr_xxy_z[i] - 12.0 * tr_xxy_xxz[i] * tke_0 - 14.0 * tr_xxxy_xz[i] * tbe_0 - 6.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_z[i] * tbe_0 + 8.0 * tr_xxxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_yy[i] = 6.0 * tr_xy_yy[i] - 12.0 * tr_xxy_xyy[i] * tke_0 - 14.0 * tr_xxxy_yy[i] * tbe_0 - 2.0 * tr_xxxy_yy[i] * tke_0 + 4.0 * tr_xxxy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_yz[i] = 6.0 * tr_xy_yz[i] - 12.0 * tr_xxy_xyz[i] * tke_0 - 14.0 * tr_xxxy_yz[i] * tbe_0 - 2.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxy_zz[i] = 6.0 * tr_xy_zz[i] - 12.0 * tr_xxy_xzz[i] * tke_0 - 14.0 * tr_xxxy_zz[i] * tbe_0 - 2.0 * tr_xxxy_zz[i] * tke_0 + 4.0 * tr_xxxy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 12-18 components of targeted buffer : GD

    auto tr_0_0_xx_xxxz_xx = pbuffer.data(idx_op_geom_020_gd + 12);

    auto tr_0_0_xx_xxxz_xy = pbuffer.data(idx_op_geom_020_gd + 13);

    auto tr_0_0_xx_xxxz_xz = pbuffer.data(idx_op_geom_020_gd + 14);

    auto tr_0_0_xx_xxxz_yy = pbuffer.data(idx_op_geom_020_gd + 15);

    auto tr_0_0_xx_xxxz_yz = pbuffer.data(idx_op_geom_020_gd + 16);

    auto tr_0_0_xx_xxxz_zz = pbuffer.data(idx_op_geom_020_gd + 17);

    #pragma omp simd aligned(tr_0_0_xx_xxxz_xx, tr_0_0_xx_xxxz_xy, tr_0_0_xx_xxxz_xz, tr_0_0_xx_xxxz_yy, tr_0_0_xx_xxxz_yz, tr_0_0_xx_xxxz_zz, tr_xxxxxz_xx, tr_xxxxxz_xy, tr_xxxxxz_xz, tr_xxxxxz_yy, tr_xxxxxz_yz, tr_xxxxxz_zz, tr_xxxxz_x, tr_xxxxz_xxx, tr_xxxxz_xxy, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxxz_xx[i] = 6.0 * tr_xz_xx[i] + 12.0 * tr_xxz_x[i] - 12.0 * tr_xxz_xxx[i] * tke_0 + 2.0 * tr_xxxz_0[i] - 14.0 * tr_xxxz_xx[i] * tbe_0 - 10.0 * tr_xxxz_xx[i] * tke_0 + 4.0 * tr_xxxz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_x[i] * tbe_0 + 8.0 * tr_xxxxz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xy[i] = 6.0 * tr_xz_xy[i] + 6.0 * tr_xxz_y[i] - 12.0 * tr_xxz_xxy[i] * tke_0 - 14.0 * tr_xxxz_xy[i] * tbe_0 - 6.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_y[i] * tbe_0 + 8.0 * tr_xxxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_xz[i] = 6.0 * tr_xz_xz[i] + 6.0 * tr_xxz_z[i] - 12.0 * tr_xxz_xxz[i] * tke_0 - 14.0 * tr_xxxz_xz[i] * tbe_0 - 6.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_z[i] * tbe_0 + 8.0 * tr_xxxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_yy[i] = 6.0 * tr_xz_yy[i] - 12.0 * tr_xxz_xyy[i] * tke_0 - 14.0 * tr_xxxz_yy[i] * tbe_0 - 2.0 * tr_xxxz_yy[i] * tke_0 + 4.0 * tr_xxxz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_yz[i] = 6.0 * tr_xz_yz[i] - 12.0 * tr_xxz_xyz[i] * tke_0 - 14.0 * tr_xxxz_yz[i] * tbe_0 - 2.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxxz_zz[i] = 6.0 * tr_xz_zz[i] - 12.0 * tr_xxz_xzz[i] * tke_0 - 14.0 * tr_xxxz_zz[i] * tbe_0 - 2.0 * tr_xxxz_zz[i] * tke_0 + 4.0 * tr_xxxz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 18-24 components of targeted buffer : GD

    auto tr_0_0_xx_xxyy_xx = pbuffer.data(idx_op_geom_020_gd + 18);

    auto tr_0_0_xx_xxyy_xy = pbuffer.data(idx_op_geom_020_gd + 19);

    auto tr_0_0_xx_xxyy_xz = pbuffer.data(idx_op_geom_020_gd + 20);

    auto tr_0_0_xx_xxyy_yy = pbuffer.data(idx_op_geom_020_gd + 21);

    auto tr_0_0_xx_xxyy_yz = pbuffer.data(idx_op_geom_020_gd + 22);

    auto tr_0_0_xx_xxyy_zz = pbuffer.data(idx_op_geom_020_gd + 23);

    #pragma omp simd aligned(tr_0_0_xx_xxyy_xx, tr_0_0_xx_xxyy_xy, tr_0_0_xx_xxyy_xz, tr_0_0_xx_xxyy_yy, tr_0_0_xx_xxyy_yz, tr_0_0_xx_xxyy_zz, tr_xxxxyy_xx, tr_xxxxyy_xy, tr_xxxxyy_xz, tr_xxxxyy_yy, tr_xxxxyy_yz, tr_xxxxyy_zz, tr_xxxyy_x, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_y, tr_xxxyy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxyy_xx[i] = 2.0 * tr_yy_xx[i] + 8.0 * tr_xyy_x[i] - 8.0 * tr_xyy_xxx[i] * tke_0 + 2.0 * tr_xxyy_0[i] - 10.0 * tr_xxyy_xx[i] * tbe_0 - 10.0 * tr_xxyy_xx[i] * tke_0 + 4.0 * tr_xxyy_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_x[i] * tbe_0 + 8.0 * tr_xxxyy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xy[i] = 2.0 * tr_yy_xy[i] + 4.0 * tr_xyy_y[i] - 8.0 * tr_xyy_xxy[i] * tke_0 - 10.0 * tr_xxyy_xy[i] * tbe_0 - 6.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_y[i] * tbe_0 + 8.0 * tr_xxxyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_xz[i] = 2.0 * tr_yy_xz[i] + 4.0 * tr_xyy_z[i] - 8.0 * tr_xyy_xxz[i] * tke_0 - 10.0 * tr_xxyy_xz[i] * tbe_0 - 6.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_z[i] * tbe_0 + 8.0 * tr_xxxyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_yy[i] = 2.0 * tr_yy_yy[i] - 8.0 * tr_xyy_xyy[i] * tke_0 - 10.0 * tr_xxyy_yy[i] * tbe_0 - 2.0 * tr_xxyy_yy[i] * tke_0 + 4.0 * tr_xxyy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_yz[i] = 2.0 * tr_yy_yz[i] - 8.0 * tr_xyy_xyz[i] * tke_0 - 10.0 * tr_xxyy_yz[i] * tbe_0 - 2.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyy_zz[i] = 2.0 * tr_yy_zz[i] - 8.0 * tr_xyy_xzz[i] * tke_0 - 10.0 * tr_xxyy_zz[i] * tbe_0 - 2.0 * tr_xxyy_zz[i] * tke_0 + 4.0 * tr_xxyy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 24-30 components of targeted buffer : GD

    auto tr_0_0_xx_xxyz_xx = pbuffer.data(idx_op_geom_020_gd + 24);

    auto tr_0_0_xx_xxyz_xy = pbuffer.data(idx_op_geom_020_gd + 25);

    auto tr_0_0_xx_xxyz_xz = pbuffer.data(idx_op_geom_020_gd + 26);

    auto tr_0_0_xx_xxyz_yy = pbuffer.data(idx_op_geom_020_gd + 27);

    auto tr_0_0_xx_xxyz_yz = pbuffer.data(idx_op_geom_020_gd + 28);

    auto tr_0_0_xx_xxyz_zz = pbuffer.data(idx_op_geom_020_gd + 29);

    #pragma omp simd aligned(tr_0_0_xx_xxyz_xx, tr_0_0_xx_xxyz_xy, tr_0_0_xx_xxyz_xz, tr_0_0_xx_xxyz_yy, tr_0_0_xx_xxyz_yz, tr_0_0_xx_xxyz_zz, tr_xxxxyz_xx, tr_xxxxyz_xy, tr_xxxxyz_xz, tr_xxxxyz_yy, tr_xxxxyz_yz, tr_xxxxyz_zz, tr_xxxyz_x, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_y, tr_xxxyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxyz_xx[i] = 2.0 * tr_yz_xx[i] + 8.0 * tr_xyz_x[i] - 8.0 * tr_xyz_xxx[i] * tke_0 + 2.0 * tr_xxyz_0[i] - 10.0 * tr_xxyz_xx[i] * tbe_0 - 10.0 * tr_xxyz_xx[i] * tke_0 + 4.0 * tr_xxyz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_x[i] * tbe_0 + 8.0 * tr_xxxyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xy[i] = 2.0 * tr_yz_xy[i] + 4.0 * tr_xyz_y[i] - 8.0 * tr_xyz_xxy[i] * tke_0 - 10.0 * tr_xxyz_xy[i] * tbe_0 - 6.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_y[i] * tbe_0 + 8.0 * tr_xxxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_xz[i] = 2.0 * tr_yz_xz[i] + 4.0 * tr_xyz_z[i] - 8.0 * tr_xyz_xxz[i] * tke_0 - 10.0 * tr_xxyz_xz[i] * tbe_0 - 6.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_z[i] * tbe_0 + 8.0 * tr_xxxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_yy[i] = 2.0 * tr_yz_yy[i] - 8.0 * tr_xyz_xyy[i] * tke_0 - 10.0 * tr_xxyz_yy[i] * tbe_0 - 2.0 * tr_xxyz_yy[i] * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_yz[i] = 2.0 * tr_yz_yz[i] - 8.0 * tr_xyz_xyz[i] * tke_0 - 10.0 * tr_xxyz_yz[i] * tbe_0 - 2.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxyz_zz[i] = 2.0 * tr_yz_zz[i] - 8.0 * tr_xyz_xzz[i] * tke_0 - 10.0 * tr_xxyz_zz[i] * tbe_0 - 2.0 * tr_xxyz_zz[i] * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-36 components of targeted buffer : GD

    auto tr_0_0_xx_xxzz_xx = pbuffer.data(idx_op_geom_020_gd + 30);

    auto tr_0_0_xx_xxzz_xy = pbuffer.data(idx_op_geom_020_gd + 31);

    auto tr_0_0_xx_xxzz_xz = pbuffer.data(idx_op_geom_020_gd + 32);

    auto tr_0_0_xx_xxzz_yy = pbuffer.data(idx_op_geom_020_gd + 33);

    auto tr_0_0_xx_xxzz_yz = pbuffer.data(idx_op_geom_020_gd + 34);

    auto tr_0_0_xx_xxzz_zz = pbuffer.data(idx_op_geom_020_gd + 35);

    #pragma omp simd aligned(tr_0_0_xx_xxzz_xx, tr_0_0_xx_xxzz_xy, tr_0_0_xx_xxzz_xz, tr_0_0_xx_xxzz_yy, tr_0_0_xx_xxzz_yz, tr_0_0_xx_xxzz_zz, tr_xxxxzz_xx, tr_xxxxzz_xy, tr_xxxxzz_xz, tr_xxxxzz_yy, tr_xxxxzz_yz, tr_xxxxzz_zz, tr_xxxzz_x, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_y, tr_xxxzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xxzz_xx[i] = 2.0 * tr_zz_xx[i] + 8.0 * tr_xzz_x[i] - 8.0 * tr_xzz_xxx[i] * tke_0 + 2.0 * tr_xxzz_0[i] - 10.0 * tr_xxzz_xx[i] * tbe_0 - 10.0 * tr_xxzz_xx[i] * tke_0 + 4.0 * tr_xxzz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_x[i] * tbe_0 + 8.0 * tr_xxxzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xy[i] = 2.0 * tr_zz_xy[i] + 4.0 * tr_xzz_y[i] - 8.0 * tr_xzz_xxy[i] * tke_0 - 10.0 * tr_xxzz_xy[i] * tbe_0 - 6.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_y[i] * tbe_0 + 8.0 * tr_xxxzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_xz[i] = 2.0 * tr_zz_xz[i] + 4.0 * tr_xzz_z[i] - 8.0 * tr_xzz_xxz[i] * tke_0 - 10.0 * tr_xxzz_xz[i] * tbe_0 - 6.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_z[i] * tbe_0 + 8.0 * tr_xxxzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_yy[i] = 2.0 * tr_zz_yy[i] - 8.0 * tr_xzz_xyy[i] * tke_0 - 10.0 * tr_xxzz_yy[i] * tbe_0 - 2.0 * tr_xxzz_yy[i] * tke_0 + 4.0 * tr_xxzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_yz[i] = 2.0 * tr_zz_yz[i] - 8.0 * tr_xzz_xyz[i] * tke_0 - 10.0 * tr_xxzz_yz[i] * tbe_0 - 2.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xxzz_zz[i] = 2.0 * tr_zz_zz[i] - 8.0 * tr_xzz_xzz[i] * tke_0 - 10.0 * tr_xxzz_zz[i] * tbe_0 - 2.0 * tr_xxzz_zz[i] * tke_0 + 4.0 * tr_xxzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 36-42 components of targeted buffer : GD

    auto tr_0_0_xx_xyyy_xx = pbuffer.data(idx_op_geom_020_gd + 36);

    auto tr_0_0_xx_xyyy_xy = pbuffer.data(idx_op_geom_020_gd + 37);

    auto tr_0_0_xx_xyyy_xz = pbuffer.data(idx_op_geom_020_gd + 38);

    auto tr_0_0_xx_xyyy_yy = pbuffer.data(idx_op_geom_020_gd + 39);

    auto tr_0_0_xx_xyyy_yz = pbuffer.data(idx_op_geom_020_gd + 40);

    auto tr_0_0_xx_xyyy_zz = pbuffer.data(idx_op_geom_020_gd + 41);

    #pragma omp simd aligned(tr_0_0_xx_xyyy_xx, tr_0_0_xx_xyyy_xy, tr_0_0_xx_xyyy_xz, tr_0_0_xx_xyyy_yy, tr_0_0_xx_xyyy_yz, tr_0_0_xx_xyyy_zz, tr_xxxyyy_xx, tr_xxxyyy_xy, tr_xxxyyy_xz, tr_xxxyyy_yy, tr_xxxyyy_yz, tr_xxxyyy_zz, tr_xxyyy_x, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_y, tr_xxyyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyyy_xx[i] = 4.0 * tr_yyy_x[i] - 4.0 * tr_yyy_xxx[i] * tke_0 + 2.0 * tr_xyyy_0[i] - 6.0 * tr_xyyy_xx[i] * tbe_0 - 10.0 * tr_xyyy_xx[i] * tke_0 + 4.0 * tr_xyyy_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_x[i] * tbe_0 + 8.0 * tr_xxyyy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xy[i] = 2.0 * tr_yyy_y[i] - 4.0 * tr_yyy_xxy[i] * tke_0 - 6.0 * tr_xyyy_xy[i] * tbe_0 - 6.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_y[i] * tbe_0 + 8.0 * tr_xxyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_xz[i] = 2.0 * tr_yyy_z[i] - 4.0 * tr_yyy_xxz[i] * tke_0 - 6.0 * tr_xyyy_xz[i] * tbe_0 - 6.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_z[i] * tbe_0 + 8.0 * tr_xxyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_yy[i] = -4.0 * tr_yyy_xyy[i] * tke_0 - 6.0 * tr_xyyy_yy[i] * tbe_0 - 2.0 * tr_xyyy_yy[i] * tke_0 + 4.0 * tr_xyyy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_yz[i] = -4.0 * tr_yyy_xyz[i] * tke_0 - 6.0 * tr_xyyy_yz[i] * tbe_0 - 2.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyy_zz[i] = -4.0 * tr_yyy_xzz[i] * tke_0 - 6.0 * tr_xyyy_zz[i] * tbe_0 - 2.0 * tr_xyyy_zz[i] * tke_0 + 4.0 * tr_xyyy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 42-48 components of targeted buffer : GD

    auto tr_0_0_xx_xyyz_xx = pbuffer.data(idx_op_geom_020_gd + 42);

    auto tr_0_0_xx_xyyz_xy = pbuffer.data(idx_op_geom_020_gd + 43);

    auto tr_0_0_xx_xyyz_xz = pbuffer.data(idx_op_geom_020_gd + 44);

    auto tr_0_0_xx_xyyz_yy = pbuffer.data(idx_op_geom_020_gd + 45);

    auto tr_0_0_xx_xyyz_yz = pbuffer.data(idx_op_geom_020_gd + 46);

    auto tr_0_0_xx_xyyz_zz = pbuffer.data(idx_op_geom_020_gd + 47);

    #pragma omp simd aligned(tr_0_0_xx_xyyz_xx, tr_0_0_xx_xyyz_xy, tr_0_0_xx_xyyz_xz, tr_0_0_xx_xyyz_yy, tr_0_0_xx_xyyz_yz, tr_0_0_xx_xyyz_zz, tr_xxxyyz_xx, tr_xxxyyz_xy, tr_xxxyyz_xz, tr_xxxyyz_yy, tr_xxxyyz_yz, tr_xxxyyz_zz, tr_xxyyz_x, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_y, tr_xxyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyyz_xx[i] = 4.0 * tr_yyz_x[i] - 4.0 * tr_yyz_xxx[i] * tke_0 + 2.0 * tr_xyyz_0[i] - 6.0 * tr_xyyz_xx[i] * tbe_0 - 10.0 * tr_xyyz_xx[i] * tke_0 + 4.0 * tr_xyyz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_x[i] * tbe_0 + 8.0 * tr_xxyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xy[i] = 2.0 * tr_yyz_y[i] - 4.0 * tr_yyz_xxy[i] * tke_0 - 6.0 * tr_xyyz_xy[i] * tbe_0 - 6.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_y[i] * tbe_0 + 8.0 * tr_xxyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_xz[i] = 2.0 * tr_yyz_z[i] - 4.0 * tr_yyz_xxz[i] * tke_0 - 6.0 * tr_xyyz_xz[i] * tbe_0 - 6.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_z[i] * tbe_0 + 8.0 * tr_xxyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_yy[i] = -4.0 * tr_yyz_xyy[i] * tke_0 - 6.0 * tr_xyyz_yy[i] * tbe_0 - 2.0 * tr_xyyz_yy[i] * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_yz[i] = -4.0 * tr_yyz_xyz[i] * tke_0 - 6.0 * tr_xyyz_yz[i] * tbe_0 - 2.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyyz_zz[i] = -4.0 * tr_yyz_xzz[i] * tke_0 - 6.0 * tr_xyyz_zz[i] * tbe_0 - 2.0 * tr_xyyz_zz[i] * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 48-54 components of targeted buffer : GD

    auto tr_0_0_xx_xyzz_xx = pbuffer.data(idx_op_geom_020_gd + 48);

    auto tr_0_0_xx_xyzz_xy = pbuffer.data(idx_op_geom_020_gd + 49);

    auto tr_0_0_xx_xyzz_xz = pbuffer.data(idx_op_geom_020_gd + 50);

    auto tr_0_0_xx_xyzz_yy = pbuffer.data(idx_op_geom_020_gd + 51);

    auto tr_0_0_xx_xyzz_yz = pbuffer.data(idx_op_geom_020_gd + 52);

    auto tr_0_0_xx_xyzz_zz = pbuffer.data(idx_op_geom_020_gd + 53);

    #pragma omp simd aligned(tr_0_0_xx_xyzz_xx, tr_0_0_xx_xyzz_xy, tr_0_0_xx_xyzz_xz, tr_0_0_xx_xyzz_yy, tr_0_0_xx_xyzz_yz, tr_0_0_xx_xyzz_zz, tr_xxxyzz_xx, tr_xxxyzz_xy, tr_xxxyzz_xz, tr_xxxyzz_yy, tr_xxxyzz_yz, tr_xxxyzz_zz, tr_xxyzz_x, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_y, tr_xxyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xyzz_xx[i] = 4.0 * tr_yzz_x[i] - 4.0 * tr_yzz_xxx[i] * tke_0 + 2.0 * tr_xyzz_0[i] - 6.0 * tr_xyzz_xx[i] * tbe_0 - 10.0 * tr_xyzz_xx[i] * tke_0 + 4.0 * tr_xyzz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_x[i] * tbe_0 + 8.0 * tr_xxyzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xy[i] = 2.0 * tr_yzz_y[i] - 4.0 * tr_yzz_xxy[i] * tke_0 - 6.0 * tr_xyzz_xy[i] * tbe_0 - 6.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_y[i] * tbe_0 + 8.0 * tr_xxyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_xz[i] = 2.0 * tr_yzz_z[i] - 4.0 * tr_yzz_xxz[i] * tke_0 - 6.0 * tr_xyzz_xz[i] * tbe_0 - 6.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_z[i] * tbe_0 + 8.0 * tr_xxyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_yy[i] = -4.0 * tr_yzz_xyy[i] * tke_0 - 6.0 * tr_xyzz_yy[i] * tbe_0 - 2.0 * tr_xyzz_yy[i] * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_yz[i] = -4.0 * tr_yzz_xyz[i] * tke_0 - 6.0 * tr_xyzz_yz[i] * tbe_0 - 2.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xyzz_zz[i] = -4.0 * tr_yzz_xzz[i] * tke_0 - 6.0 * tr_xyzz_zz[i] * tbe_0 - 2.0 * tr_xyzz_zz[i] * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 54-60 components of targeted buffer : GD

    auto tr_0_0_xx_xzzz_xx = pbuffer.data(idx_op_geom_020_gd + 54);

    auto tr_0_0_xx_xzzz_xy = pbuffer.data(idx_op_geom_020_gd + 55);

    auto tr_0_0_xx_xzzz_xz = pbuffer.data(idx_op_geom_020_gd + 56);

    auto tr_0_0_xx_xzzz_yy = pbuffer.data(idx_op_geom_020_gd + 57);

    auto tr_0_0_xx_xzzz_yz = pbuffer.data(idx_op_geom_020_gd + 58);

    auto tr_0_0_xx_xzzz_zz = pbuffer.data(idx_op_geom_020_gd + 59);

    #pragma omp simd aligned(tr_0_0_xx_xzzz_xx, tr_0_0_xx_xzzz_xy, tr_0_0_xx_xzzz_xz, tr_0_0_xx_xzzz_yy, tr_0_0_xx_xzzz_yz, tr_0_0_xx_xzzz_zz, tr_xxxzzz_xx, tr_xxxzzz_xy, tr_xxxzzz_xz, tr_xxxzzz_yy, tr_xxxzzz_yz, tr_xxxzzz_zz, tr_xxzzz_x, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_y, tr_xxzzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_xzzz_xx[i] = 4.0 * tr_zzz_x[i] - 4.0 * tr_zzz_xxx[i] * tke_0 + 2.0 * tr_xzzz_0[i] - 6.0 * tr_xzzz_xx[i] * tbe_0 - 10.0 * tr_xzzz_xx[i] * tke_0 + 4.0 * tr_xzzz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_x[i] * tbe_0 + 8.0 * tr_xxzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xy[i] = 2.0 * tr_zzz_y[i] - 4.0 * tr_zzz_xxy[i] * tke_0 - 6.0 * tr_xzzz_xy[i] * tbe_0 - 6.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_y[i] * tbe_0 + 8.0 * tr_xxzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_xz[i] = 2.0 * tr_zzz_z[i] - 4.0 * tr_zzz_xxz[i] * tke_0 - 6.0 * tr_xzzz_xz[i] * tbe_0 - 6.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_z[i] * tbe_0 + 8.0 * tr_xxzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_yy[i] = -4.0 * tr_zzz_xyy[i] * tke_0 - 6.0 * tr_xzzz_yy[i] * tbe_0 - 2.0 * tr_xzzz_yy[i] * tke_0 + 4.0 * tr_xzzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_yz[i] = -4.0 * tr_zzz_xyz[i] * tke_0 - 6.0 * tr_xzzz_yz[i] * tbe_0 - 2.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_xzzz_zz[i] = -4.0 * tr_zzz_xzz[i] * tke_0 - 6.0 * tr_xzzz_zz[i] * tbe_0 - 2.0 * tr_xzzz_zz[i] * tke_0 + 4.0 * tr_xzzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-66 components of targeted buffer : GD

    auto tr_0_0_xx_yyyy_xx = pbuffer.data(idx_op_geom_020_gd + 60);

    auto tr_0_0_xx_yyyy_xy = pbuffer.data(idx_op_geom_020_gd + 61);

    auto tr_0_0_xx_yyyy_xz = pbuffer.data(idx_op_geom_020_gd + 62);

    auto tr_0_0_xx_yyyy_yy = pbuffer.data(idx_op_geom_020_gd + 63);

    auto tr_0_0_xx_yyyy_yz = pbuffer.data(idx_op_geom_020_gd + 64);

    auto tr_0_0_xx_yyyy_zz = pbuffer.data(idx_op_geom_020_gd + 65);

    #pragma omp simd aligned(tr_0_0_xx_yyyy_xx, tr_0_0_xx_yyyy_xy, tr_0_0_xx_yyyy_xz, tr_0_0_xx_yyyy_yy, tr_0_0_xx_yyyy_yz, tr_0_0_xx_yyyy_zz, tr_xxyyyy_xx, tr_xxyyyy_xy, tr_xxyyyy_xz, tr_xxyyyy_yy, tr_xxyyyy_yz, tr_xxyyyy_zz, tr_xyyyy_x, tr_xyyyy_xxx, tr_xyyyy_xxy, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_y, tr_xyyyy_z, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyyy_xx[i] = 2.0 * tr_yyyy_0[i] - 2.0 * tr_yyyy_xx[i] * tbe_0 - 10.0 * tr_yyyy_xx[i] * tke_0 + 4.0 * tr_yyyy_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_x[i] * tbe_0 + 8.0 * tr_xyyyy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xy[i] = -2.0 * tr_yyyy_xy[i] * tbe_0 - 6.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_y[i] * tbe_0 + 8.0 * tr_xyyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_xz[i] = -2.0 * tr_yyyy_xz[i] * tbe_0 - 6.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_z[i] * tbe_0 + 8.0 * tr_xyyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_yy[i] = -2.0 * tr_yyyy_yy[i] * tbe_0 - 2.0 * tr_yyyy_yy[i] * tke_0 + 4.0 * tr_yyyy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_yz[i] = -2.0 * tr_yyyy_yz[i] * tbe_0 - 2.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyy_zz[i] = -2.0 * tr_yyyy_zz[i] * tbe_0 - 2.0 * tr_yyyy_zz[i] * tke_0 + 4.0 * tr_yyyy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 66-72 components of targeted buffer : GD

    auto tr_0_0_xx_yyyz_xx = pbuffer.data(idx_op_geom_020_gd + 66);

    auto tr_0_0_xx_yyyz_xy = pbuffer.data(idx_op_geom_020_gd + 67);

    auto tr_0_0_xx_yyyz_xz = pbuffer.data(idx_op_geom_020_gd + 68);

    auto tr_0_0_xx_yyyz_yy = pbuffer.data(idx_op_geom_020_gd + 69);

    auto tr_0_0_xx_yyyz_yz = pbuffer.data(idx_op_geom_020_gd + 70);

    auto tr_0_0_xx_yyyz_zz = pbuffer.data(idx_op_geom_020_gd + 71);

    #pragma omp simd aligned(tr_0_0_xx_yyyz_xx, tr_0_0_xx_yyyz_xy, tr_0_0_xx_yyyz_xz, tr_0_0_xx_yyyz_yy, tr_0_0_xx_yyyz_yz, tr_0_0_xx_yyyz_zz, tr_xxyyyz_xx, tr_xxyyyz_xy, tr_xxyyyz_xz, tr_xxyyyz_yy, tr_xxyyyz_yz, tr_xxyyyz_zz, tr_xyyyz_x, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_y, tr_xyyyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyyz_xx[i] = 2.0 * tr_yyyz_0[i] - 2.0 * tr_yyyz_xx[i] * tbe_0 - 10.0 * tr_yyyz_xx[i] * tke_0 + 4.0 * tr_yyyz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_x[i] * tbe_0 + 8.0 * tr_xyyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xy[i] = -2.0 * tr_yyyz_xy[i] * tbe_0 - 6.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_y[i] * tbe_0 + 8.0 * tr_xyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_xz[i] = -2.0 * tr_yyyz_xz[i] * tbe_0 - 6.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_z[i] * tbe_0 + 8.0 * tr_xyyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_yy[i] = -2.0 * tr_yyyz_yy[i] * tbe_0 - 2.0 * tr_yyyz_yy[i] * tke_0 + 4.0 * tr_yyyz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_yz[i] = -2.0 * tr_yyyz_yz[i] * tbe_0 - 2.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyyz_zz[i] = -2.0 * tr_yyyz_zz[i] * tbe_0 - 2.0 * tr_yyyz_zz[i] * tke_0 + 4.0 * tr_yyyz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 72-78 components of targeted buffer : GD

    auto tr_0_0_xx_yyzz_xx = pbuffer.data(idx_op_geom_020_gd + 72);

    auto tr_0_0_xx_yyzz_xy = pbuffer.data(idx_op_geom_020_gd + 73);

    auto tr_0_0_xx_yyzz_xz = pbuffer.data(idx_op_geom_020_gd + 74);

    auto tr_0_0_xx_yyzz_yy = pbuffer.data(idx_op_geom_020_gd + 75);

    auto tr_0_0_xx_yyzz_yz = pbuffer.data(idx_op_geom_020_gd + 76);

    auto tr_0_0_xx_yyzz_zz = pbuffer.data(idx_op_geom_020_gd + 77);

    #pragma omp simd aligned(tr_0_0_xx_yyzz_xx, tr_0_0_xx_yyzz_xy, tr_0_0_xx_yyzz_xz, tr_0_0_xx_yyzz_yy, tr_0_0_xx_yyzz_yz, tr_0_0_xx_yyzz_zz, tr_xxyyzz_xx, tr_xxyyzz_xy, tr_xxyyzz_xz, tr_xxyyzz_yy, tr_xxyyzz_yz, tr_xxyyzz_zz, tr_xyyzz_x, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_y, tr_xyyzz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yyzz_xx[i] = 2.0 * tr_yyzz_0[i] - 2.0 * tr_yyzz_xx[i] * tbe_0 - 10.0 * tr_yyzz_xx[i] * tke_0 + 4.0 * tr_yyzz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_x[i] * tbe_0 + 8.0 * tr_xyyzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xy[i] = -2.0 * tr_yyzz_xy[i] * tbe_0 - 6.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_y[i] * tbe_0 + 8.0 * tr_xyyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_xz[i] = -2.0 * tr_yyzz_xz[i] * tbe_0 - 6.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_z[i] * tbe_0 + 8.0 * tr_xyyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_yy[i] = -2.0 * tr_yyzz_yy[i] * tbe_0 - 2.0 * tr_yyzz_yy[i] * tke_0 + 4.0 * tr_yyzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_yz[i] = -2.0 * tr_yyzz_yz[i] * tbe_0 - 2.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yyzz_zz[i] = -2.0 * tr_yyzz_zz[i] * tbe_0 - 2.0 * tr_yyzz_zz[i] * tke_0 + 4.0 * tr_yyzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 78-84 components of targeted buffer : GD

    auto tr_0_0_xx_yzzz_xx = pbuffer.data(idx_op_geom_020_gd + 78);

    auto tr_0_0_xx_yzzz_xy = pbuffer.data(idx_op_geom_020_gd + 79);

    auto tr_0_0_xx_yzzz_xz = pbuffer.data(idx_op_geom_020_gd + 80);

    auto tr_0_0_xx_yzzz_yy = pbuffer.data(idx_op_geom_020_gd + 81);

    auto tr_0_0_xx_yzzz_yz = pbuffer.data(idx_op_geom_020_gd + 82);

    auto tr_0_0_xx_yzzz_zz = pbuffer.data(idx_op_geom_020_gd + 83);

    #pragma omp simd aligned(tr_0_0_xx_yzzz_xx, tr_0_0_xx_yzzz_xy, tr_0_0_xx_yzzz_xz, tr_0_0_xx_yzzz_yy, tr_0_0_xx_yzzz_yz, tr_0_0_xx_yzzz_zz, tr_xxyzzz_xx, tr_xxyzzz_xy, tr_xxyzzz_xz, tr_xxyzzz_yy, tr_xxyzzz_yz, tr_xxyzzz_zz, tr_xyzzz_x, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_y, tr_xyzzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_yzzz_xx[i] = 2.0 * tr_yzzz_0[i] - 2.0 * tr_yzzz_xx[i] * tbe_0 - 10.0 * tr_yzzz_xx[i] * tke_0 + 4.0 * tr_yzzz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_x[i] * tbe_0 + 8.0 * tr_xyzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xy[i] = -2.0 * tr_yzzz_xy[i] * tbe_0 - 6.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_y[i] * tbe_0 + 8.0 * tr_xyzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_xz[i] = -2.0 * tr_yzzz_xz[i] * tbe_0 - 6.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_z[i] * tbe_0 + 8.0 * tr_xyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_yy[i] = -2.0 * tr_yzzz_yy[i] * tbe_0 - 2.0 * tr_yzzz_yy[i] * tke_0 + 4.0 * tr_yzzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_yz[i] = -2.0 * tr_yzzz_yz[i] * tbe_0 - 2.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_yzzz_zz[i] = -2.0 * tr_yzzz_zz[i] * tbe_0 - 2.0 * tr_yzzz_zz[i] * tke_0 + 4.0 * tr_yzzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 84-90 components of targeted buffer : GD

    auto tr_0_0_xx_zzzz_xx = pbuffer.data(idx_op_geom_020_gd + 84);

    auto tr_0_0_xx_zzzz_xy = pbuffer.data(idx_op_geom_020_gd + 85);

    auto tr_0_0_xx_zzzz_xz = pbuffer.data(idx_op_geom_020_gd + 86);

    auto tr_0_0_xx_zzzz_yy = pbuffer.data(idx_op_geom_020_gd + 87);

    auto tr_0_0_xx_zzzz_yz = pbuffer.data(idx_op_geom_020_gd + 88);

    auto tr_0_0_xx_zzzz_zz = pbuffer.data(idx_op_geom_020_gd + 89);

    #pragma omp simd aligned(tr_0_0_xx_zzzz_xx, tr_0_0_xx_zzzz_xy, tr_0_0_xx_zzzz_xz, tr_0_0_xx_zzzz_yy, tr_0_0_xx_zzzz_yz, tr_0_0_xx_zzzz_zz, tr_xxzzzz_xx, tr_xxzzzz_xy, tr_xxzzzz_xz, tr_xxzzzz_yy, tr_xxzzzz_yz, tr_xxzzzz_zz, tr_xzzzz_x, tr_xzzzz_xxx, tr_xzzzz_xxy, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_y, tr_xzzzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xx_zzzz_xx[i] = 2.0 * tr_zzzz_0[i] - 2.0 * tr_zzzz_xx[i] * tbe_0 - 10.0 * tr_zzzz_xx[i] * tke_0 + 4.0 * tr_zzzz_xxxx[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_x[i] * tbe_0 + 8.0 * tr_xzzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xy[i] = -2.0 * tr_zzzz_xy[i] * tbe_0 - 6.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_y[i] * tbe_0 + 8.0 * tr_xzzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_xz[i] = -2.0 * tr_zzzz_xz[i] * tbe_0 - 6.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_z[i] * tbe_0 + 8.0 * tr_xzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_yy[i] = -2.0 * tr_zzzz_yy[i] * tbe_0 - 2.0 * tr_zzzz_yy[i] * tke_0 + 4.0 * tr_zzzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_yz[i] = -2.0 * tr_zzzz_yz[i] * tbe_0 - 2.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_xxyz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xx_zzzz_zz[i] = -2.0 * tr_zzzz_zz[i] * tbe_0 - 2.0 * tr_zzzz_zz[i] * tke_0 + 4.0 * tr_zzzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-96 components of targeted buffer : GD

    auto tr_0_0_xy_xxxx_xx = pbuffer.data(idx_op_geom_020_gd + 90);

    auto tr_0_0_xy_xxxx_xy = pbuffer.data(idx_op_geom_020_gd + 91);

    auto tr_0_0_xy_xxxx_xz = pbuffer.data(idx_op_geom_020_gd + 92);

    auto tr_0_0_xy_xxxx_yy = pbuffer.data(idx_op_geom_020_gd + 93);

    auto tr_0_0_xy_xxxx_yz = pbuffer.data(idx_op_geom_020_gd + 94);

    auto tr_0_0_xy_xxxx_zz = pbuffer.data(idx_op_geom_020_gd + 95);

    #pragma omp simd aligned(tr_0_0_xy_xxxx_xx, tr_0_0_xy_xxxx_xy, tr_0_0_xy_xxxx_xz, tr_0_0_xy_xxxx_yy, tr_0_0_xy_xxxx_yz, tr_0_0_xy_xxxx_zz, tr_xxx_x, tr_xxx_xxy, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xxxy, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxxx_x, tr_xxxxx_xxy, tr_xxxxx_xyy, tr_xxxxx_xyz, tr_xxxxx_y, tr_xxxxx_yyy, tr_xxxxx_yyz, tr_xxxxx_yzz, tr_xxxxx_z, tr_xxxxxy_xx, tr_xxxxxy_xy, tr_xxxxxy_xz, tr_xxxxxy_yy, tr_xxxxxy_yz, tr_xxxxxy_zz, tr_xxxxy_x, tr_xxxxy_xxx, tr_xxxxy_xxy, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_y, tr_xxxxy_z, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxxx_xx[i] = -8.0 * tr_xxx_xxy[i] * tke_0 - 8.0 * tr_xxxy_xx[i] * tbe_0 - 4.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_x[i] * tbe_0 + 4.0 * tr_xxxxy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xy[i] = 4.0 * tr_xxx_x[i] - 8.0 * tr_xxx_xyy[i] * tke_0 - 8.0 * tr_xxxy_xy[i] * tbe_0 + tr_xxxx_0[i] - 2.0 * tr_xxxx_yy[i] * tke_0 - 2.0 * tr_xxxx_xx[i] * tke_0 + 4.0 * tr_xxxx_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxxy_y[i] * tbe_0 + 4.0 * tr_xxxxy_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_x[i] * tbe_0 + 4.0 * tr_xxxxx_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_xz[i] = -8.0 * tr_xxx_xyz[i] * tke_0 - 8.0 * tr_xxxy_xz[i] * tbe_0 - 2.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxy_z[i] * tbe_0 + 4.0 * tr_xxxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_yy[i] = 8.0 * tr_xxx_y[i] - 8.0 * tr_xxx_yyy[i] * tke_0 - 8.0 * tr_xxxy_yy[i] * tbe_0 - 4.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_y[i] * tbe_0 + 4.0 * tr_xxxxx_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_yz[i] = 4.0 * tr_xxx_z[i] - 8.0 * tr_xxx_yyz[i] * tke_0 - 8.0 * tr_xxxy_yz[i] * tbe_0 - 2.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_z[i] * tbe_0 + 4.0 * tr_xxxxx_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxx_zz[i] = -8.0 * tr_xxx_yzz[i] * tke_0 - 8.0 * tr_xxxy_zz[i] * tbe_0 + 4.0 * tr_xxxx_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 96-102 components of targeted buffer : GD

    auto tr_0_0_xy_xxxy_xx = pbuffer.data(idx_op_geom_020_gd + 96);

    auto tr_0_0_xy_xxxy_xy = pbuffer.data(idx_op_geom_020_gd + 97);

    auto tr_0_0_xy_xxxy_xz = pbuffer.data(idx_op_geom_020_gd + 98);

    auto tr_0_0_xy_xxxy_yy = pbuffer.data(idx_op_geom_020_gd + 99);

    auto tr_0_0_xy_xxxy_yz = pbuffer.data(idx_op_geom_020_gd + 100);

    auto tr_0_0_xy_xxxy_zz = pbuffer.data(idx_op_geom_020_gd + 101);

    #pragma omp simd aligned(tr_0_0_xy_xxxy_xx, tr_0_0_xy_xxxy_xy, tr_0_0_xy_xxxy_xz, tr_0_0_xy_xxxy_yy, tr_0_0_xy_xxxy_yz, tr_0_0_xy_xxxy_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_z, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxy_x, tr_xxxxy_xxy, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_y, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_z, tr_xxxxyy_xx, tr_xxxxyy_xy, tr_xxxxyy_xz, tr_xxxxyy_yy, tr_xxxxyy_yz, tr_xxxxyy_zz, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xxxy, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxyy_x, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_y, tr_xxxyy_z, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxxy_xx[i] = 3.0 * tr_xx_xx[i] - 6.0 * tr_xxy_xxy[i] * tke_0 - 6.0 * tr_xxyy_xx[i] * tbe_0 + 2.0 * tr_xxx_x[i] - 2.0 * tr_xxx_xxx[i] * tke_0 - 4.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_x[i] * tbe_0 + 4.0 * tr_xxxyy_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xy[i] = 3.0 * tr_xx_xy[i] + 3.0 * tr_xxy_x[i] - 6.0 * tr_xxy_xyy[i] * tke_0 - 6.0 * tr_xxyy_xy[i] * tbe_0 + tr_xxx_y[i] - 2.0 * tr_xxx_xxy[i] * tke_0 + tr_xxxy_0[i] - 2.0 * tr_xxxy_yy[i] * tke_0 - 2.0 * tr_xxxy_xx[i] * tke_0 + 4.0 * tr_xxxy_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxyy_y[i] * tbe_0 + 4.0 * tr_xxxyy_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xy[i] * tbe_0 - 2.0 * tr_xxxxy_x[i] * tbe_0 + 4.0 * tr_xxxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_xz[i] = 3.0 * tr_xx_xz[i] - 6.0 * tr_xxy_xyz[i] * tke_0 - 6.0 * tr_xxyy_xz[i] * tbe_0 + tr_xxx_z[i] - 2.0 * tr_xxx_xxz[i] * tke_0 - 2.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyy_z[i] * tbe_0 + 4.0 * tr_xxxyy_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_yy[i] = 3.0 * tr_xx_yy[i] + 6.0 * tr_xxy_y[i] - 6.0 * tr_xxy_yyy[i] * tke_0 - 6.0 * tr_xxyy_yy[i] * tbe_0 - 2.0 * tr_xxx_xyy[i] * tke_0 - 4.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yy[i] * tbe_0 - 4.0 * tr_xxxxy_y[i] * tbe_0 + 4.0 * tr_xxxxy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_yz[i] = 3.0 * tr_xx_yz[i] + 3.0 * tr_xxy_z[i] - 6.0 * tr_xxy_yyz[i] * tke_0 - 6.0 * tr_xxyy_yz[i] * tbe_0 - 2.0 * tr_xxx_xyz[i] * tke_0 - 2.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yz[i] * tbe_0 - 2.0 * tr_xxxxy_z[i] * tbe_0 + 4.0 * tr_xxxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxy_zz[i] = 3.0 * tr_xx_zz[i] - 6.0 * tr_xxy_yzz[i] * tke_0 - 6.0 * tr_xxyy_zz[i] * tbe_0 - 2.0 * tr_xxx_xzz[i] * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyy_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zz[i] * tbe_0 + 4.0 * tr_xxxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 102-108 components of targeted buffer : GD

    auto tr_0_0_xy_xxxz_xx = pbuffer.data(idx_op_geom_020_gd + 102);

    auto tr_0_0_xy_xxxz_xy = pbuffer.data(idx_op_geom_020_gd + 103);

    auto tr_0_0_xy_xxxz_xz = pbuffer.data(idx_op_geom_020_gd + 104);

    auto tr_0_0_xy_xxxz_yy = pbuffer.data(idx_op_geom_020_gd + 105);

    auto tr_0_0_xy_xxxz_yz = pbuffer.data(idx_op_geom_020_gd + 106);

    auto tr_0_0_xy_xxxz_zz = pbuffer.data(idx_op_geom_020_gd + 107);

    #pragma omp simd aligned(tr_0_0_xy_xxxz_xx, tr_0_0_xy_xxxz_xy, tr_0_0_xy_xxxz_xz, tr_0_0_xy_xxxz_yy, tr_0_0_xy_xxxz_yz, tr_0_0_xy_xxxz_zz, tr_xxxxyz_xx, tr_xxxxyz_xy, tr_xxxxyz_xz, tr_xxxxyz_yy, tr_xxxxyz_yz, tr_xxxxyz_zz, tr_xxxxz_x, tr_xxxxz_xxy, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_y, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_z, tr_xxxyz_x, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_y, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xxxy, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxxz_xx[i] = -6.0 * tr_xxz_xxy[i] * tke_0 - 6.0 * tr_xxyz_xx[i] * tbe_0 - 4.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xy[i] = 3.0 * tr_xxz_x[i] - 6.0 * tr_xxz_xyy[i] * tke_0 - 6.0 * tr_xxyz_xy[i] * tbe_0 + tr_xxxz_0[i] - 2.0 * tr_xxxz_yy[i] * tke_0 - 2.0 * tr_xxxz_xx[i] * tke_0 + 4.0 * tr_xxxz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxz_x[i] * tbe_0 + 4.0 * tr_xxxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_xz[i] = -6.0 * tr_xxz_xyz[i] * tke_0 - 6.0 * tr_xxyz_xz[i] * tbe_0 - 2.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_yy[i] = 6.0 * tr_xxz_y[i] - 6.0 * tr_xxz_yyy[i] * tke_0 - 6.0 * tr_xxyz_yy[i] * tbe_0 - 4.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxz_y[i] * tbe_0 + 4.0 * tr_xxxxz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_yz[i] = 3.0 * tr_xxz_z[i] - 6.0 * tr_xxz_yyz[i] * tke_0 - 6.0 * tr_xxyz_yz[i] * tbe_0 - 2.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxz_z[i] * tbe_0 + 4.0 * tr_xxxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxxz_zz[i] = -6.0 * tr_xxz_yzz[i] * tke_0 - 6.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxxz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 108-114 components of targeted buffer : GD

    auto tr_0_0_xy_xxyy_xx = pbuffer.data(idx_op_geom_020_gd + 108);

    auto tr_0_0_xy_xxyy_xy = pbuffer.data(idx_op_geom_020_gd + 109);

    auto tr_0_0_xy_xxyy_xz = pbuffer.data(idx_op_geom_020_gd + 110);

    auto tr_0_0_xy_xxyy_yy = pbuffer.data(idx_op_geom_020_gd + 111);

    auto tr_0_0_xy_xxyy_yz = pbuffer.data(idx_op_geom_020_gd + 112);

    auto tr_0_0_xy_xxyy_zz = pbuffer.data(idx_op_geom_020_gd + 113);

    #pragma omp simd aligned(tr_0_0_xy_xxyy_xx, tr_0_0_xy_xxyy_xy, tr_0_0_xy_xxyy_xz, tr_0_0_xy_xxyy_yy, tr_0_0_xy_xxyy_yz, tr_0_0_xy_xxyy_zz, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyy_x, tr_xxxyy_xxy, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_y, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_z, tr_xxxyyy_xx, tr_xxxyyy_xy, tr_xxxyyy_xz, tr_xxxyyy_yy, tr_xxxyyy_yz, tr_xxxyyy_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xxxy, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyyy_x, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_y, tr_xxyyy_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxyy_xx[i] = 4.0 * tr_xy_xx[i] - 4.0 * tr_xyy_xxy[i] * tke_0 - 4.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xxy_x[i] - 4.0 * tr_xxy_xxx[i] * tke_0 - 4.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_x[i] * tbe_0 + 4.0 * tr_xxyyy_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xy[i] = 4.0 * tr_xy_xy[i] + 2.0 * tr_xyy_x[i] - 4.0 * tr_xyy_xyy[i] * tke_0 - 4.0 * tr_xyyy_xy[i] * tbe_0 + 2.0 * tr_xxy_y[i] - 4.0 * tr_xxy_xxy[i] * tke_0 + tr_xxyy_0[i] - 2.0 * tr_xxyy_yy[i] * tke_0 - 2.0 * tr_xxyy_xx[i] * tke_0 + 4.0 * tr_xxyy_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyyy_y[i] * tbe_0 + 4.0 * tr_xxyyy_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xy[i] * tbe_0 - 2.0 * tr_xxxyy_x[i] * tbe_0 + 4.0 * tr_xxxyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_xz[i] = 4.0 * tr_xy_xz[i] - 4.0 * tr_xyy_xyz[i] * tke_0 - 4.0 * tr_xyyy_xz[i] * tbe_0 + 2.0 * tr_xxy_z[i] - 4.0 * tr_xxy_xxz[i] * tke_0 - 2.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_yy[i] = 4.0 * tr_xy_yy[i] + 4.0 * tr_xyy_y[i] - 4.0 * tr_xyy_yyy[i] * tke_0 - 4.0 * tr_xyyy_yy[i] * tbe_0 - 4.0 * tr_xxy_xyy[i] * tke_0 - 4.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yy[i] * tbe_0 - 4.0 * tr_xxxyy_y[i] * tbe_0 + 4.0 * tr_xxxyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_yz[i] = 4.0 * tr_xy_yz[i] + 2.0 * tr_xyy_z[i] - 4.0 * tr_xyy_yyz[i] * tke_0 - 4.0 * tr_xyyy_yz[i] * tbe_0 - 4.0 * tr_xxy_xyz[i] * tke_0 - 2.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yz[i] * tbe_0 - 2.0 * tr_xxxyy_z[i] * tbe_0 + 4.0 * tr_xxxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyy_zz[i] = 4.0 * tr_xy_zz[i] - 4.0 * tr_xyy_yzz[i] * tke_0 - 4.0 * tr_xyyy_zz[i] * tbe_0 - 4.0 * tr_xxy_xzz[i] * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyy_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 114-120 components of targeted buffer : GD

    auto tr_0_0_xy_xxyz_xx = pbuffer.data(idx_op_geom_020_gd + 114);

    auto tr_0_0_xy_xxyz_xy = pbuffer.data(idx_op_geom_020_gd + 115);

    auto tr_0_0_xy_xxyz_xz = pbuffer.data(idx_op_geom_020_gd + 116);

    auto tr_0_0_xy_xxyz_yy = pbuffer.data(idx_op_geom_020_gd + 117);

    auto tr_0_0_xy_xxyz_yz = pbuffer.data(idx_op_geom_020_gd + 118);

    auto tr_0_0_xy_xxyz_zz = pbuffer.data(idx_op_geom_020_gd + 119);

    #pragma omp simd aligned(tr_0_0_xy_xxyz_xx, tr_0_0_xy_xxyz_xy, tr_0_0_xy_xxyz_xz, tr_0_0_xy_xxyz_yy, tr_0_0_xy_xxyz_yz, tr_0_0_xy_xxyz_zz, tr_xxxyyz_xx, tr_xxxyyz_xy, tr_xxxyyz_xz, tr_xxxyyz_yy, tr_xxxyyz_yz, tr_xxxyyz_zz, tr_xxxyz_x, tr_xxxyz_xxy, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_y, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_z, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxyyz_x, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_y, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xxxy, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_z, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxyz_xx[i] = 2.0 * tr_xz_xx[i] - 4.0 * tr_xyz_xxy[i] * tke_0 - 4.0 * tr_xyyz_xx[i] * tbe_0 + 2.0 * tr_xxz_x[i] - 2.0 * tr_xxz_xxx[i] * tke_0 - 4.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xy[i] = 2.0 * tr_xz_xy[i] + 2.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xyy[i] * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + tr_xxz_y[i] - 2.0 * tr_xxz_xxy[i] * tke_0 + tr_xxyz_0[i] - 2.0 * tr_xxyz_yy[i] * tke_0 - 2.0 * tr_xxyz_xx[i] * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xy[i] * tbe_0 - 2.0 * tr_xxxyz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_xz[i] = 2.0 * tr_xz_xz[i] - 4.0 * tr_xyz_xyz[i] * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + tr_xxz_z[i] - 2.0 * tr_xxz_xxz[i] * tke_0 - 2.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_yy[i] = 2.0 * tr_xz_yy[i] + 4.0 * tr_xyz_y[i] - 4.0 * tr_xyz_yyy[i] * tke_0 - 4.0 * tr_xyyz_yy[i] * tbe_0 - 2.0 * tr_xxz_xyy[i] * tke_0 - 4.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yy[i] * tbe_0 - 4.0 * tr_xxxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_yz[i] = 2.0 * tr_xz_yz[i] + 2.0 * tr_xyz_z[i] - 4.0 * tr_xyz_yyz[i] * tke_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 - 2.0 * tr_xxz_xyz[i] * tke_0 - 2.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yz[i] * tbe_0 - 2.0 * tr_xxxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxyz_zz[i] = 2.0 * tr_xz_zz[i] - 4.0 * tr_xyz_yzz[i] * tke_0 - 4.0 * tr_xyyz_zz[i] * tbe_0 - 2.0 * tr_xxz_xzz[i] * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-126 components of targeted buffer : GD

    auto tr_0_0_xy_xxzz_xx = pbuffer.data(idx_op_geom_020_gd + 120);

    auto tr_0_0_xy_xxzz_xy = pbuffer.data(idx_op_geom_020_gd + 121);

    auto tr_0_0_xy_xxzz_xz = pbuffer.data(idx_op_geom_020_gd + 122);

    auto tr_0_0_xy_xxzz_yy = pbuffer.data(idx_op_geom_020_gd + 123);

    auto tr_0_0_xy_xxzz_yz = pbuffer.data(idx_op_geom_020_gd + 124);

    auto tr_0_0_xy_xxzz_zz = pbuffer.data(idx_op_geom_020_gd + 125);

    #pragma omp simd aligned(tr_0_0_xy_xxzz_xx, tr_0_0_xy_xxzz_xy, tr_0_0_xy_xxzz_xz, tr_0_0_xy_xxzz_yy, tr_0_0_xy_xxzz_yz, tr_0_0_xy_xxzz_zz, tr_xxxyzz_xx, tr_xxxyzz_xy, tr_xxxyzz_xz, tr_xxxyzz_yy, tr_xxxyzz_yz, tr_xxxyzz_zz, tr_xxxzz_x, tr_xxxzz_xxy, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_y, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_z, tr_xxyzz_x, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_y, tr_xxyzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xxxy, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xxzz_xx[i] = -4.0 * tr_xzz_xxy[i] * tke_0 - 4.0 * tr_xyzz_xx[i] * tbe_0 - 4.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xy[i] = 2.0 * tr_xzz_x[i] - 4.0 * tr_xzz_xyy[i] * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + tr_xxzz_0[i] - 2.0 * tr_xxzz_yy[i] * tke_0 - 2.0 * tr_xxzz_xx[i] * tke_0 + 4.0 * tr_xxzz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxzz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_xz[i] = -4.0 * tr_xzz_xyz[i] * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 - 2.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_yy[i] = 4.0 * tr_xzz_y[i] - 4.0 * tr_xzz_yyy[i] * tke_0 - 4.0 * tr_xyzz_yy[i] * tbe_0 - 4.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxzz_y[i] * tbe_0 + 4.0 * tr_xxxzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_yz[i] = 2.0 * tr_xzz_z[i] - 4.0 * tr_xzz_yyz[i] * tke_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 - 2.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxzz_z[i] * tbe_0 + 4.0 * tr_xxxzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xxzz_zz[i] = -4.0 * tr_xzz_yzz[i] * tke_0 - 4.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 126-132 components of targeted buffer : GD

    auto tr_0_0_xy_xyyy_xx = pbuffer.data(idx_op_geom_020_gd + 126);

    auto tr_0_0_xy_xyyy_xy = pbuffer.data(idx_op_geom_020_gd + 127);

    auto tr_0_0_xy_xyyy_xz = pbuffer.data(idx_op_geom_020_gd + 128);

    auto tr_0_0_xy_xyyy_yy = pbuffer.data(idx_op_geom_020_gd + 129);

    auto tr_0_0_xy_xyyy_yz = pbuffer.data(idx_op_geom_020_gd + 130);

    auto tr_0_0_xy_xyyy_zz = pbuffer.data(idx_op_geom_020_gd + 131);

    #pragma omp simd aligned(tr_0_0_xy_xyyy_xx, tr_0_0_xy_xyyy_xy, tr_0_0_xy_xyyy_xz, tr_0_0_xy_xyyy_yy, tr_0_0_xy_xyyy_yz, tr_0_0_xy_xyyy_zz, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyy_x, tr_xxyyy_xxy, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_y, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_z, tr_xxyyyy_xx, tr_xxyyyy_xy, tr_xxyyyy_xz, tr_xxyyyy_yy, tr_xxyyyy_yz, tr_xxyyyy_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xxxy, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyyy_x, tr_xyyyy_xxx, tr_xyyyy_xxy, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_y, tr_xyyyy_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyy_x, tr_yyy_xxy, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyyy_xx[i] = 3.0 * tr_yy_xx[i] - 2.0 * tr_yyy_xxy[i] * tke_0 - 2.0 * tr_yyyy_xx[i] * tbe_0 + 6.0 * tr_xyy_x[i] - 6.0 * tr_xyy_xxx[i] * tke_0 - 4.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_xxx[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xy[i] = 3.0 * tr_yy_xy[i] + tr_yyy_x[i] - 2.0 * tr_yyy_xyy[i] * tke_0 - 2.0 * tr_yyyy_xy[i] * tbe_0 + 3.0 * tr_xyy_y[i] - 6.0 * tr_xyy_xxy[i] * tke_0 + tr_xyyy_0[i] - 2.0 * tr_xyyy_yy[i] * tke_0 - 2.0 * tr_xyyy_xx[i] * tke_0 + 4.0 * tr_xyyy_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyyy_y[i] * tbe_0 + 4.0 * tr_xyyyy_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xy[i] * tbe_0 - 2.0 * tr_xxyyy_x[i] * tbe_0 + 4.0 * tr_xxyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_xz[i] = 3.0 * tr_yy_xz[i] - 2.0 * tr_yyy_xyz[i] * tke_0 - 2.0 * tr_yyyy_xz[i] * tbe_0 + 3.0 * tr_xyy_z[i] - 6.0 * tr_xyy_xxz[i] * tke_0 - 2.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_xxz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_yy[i] = 3.0 * tr_yy_yy[i] + 2.0 * tr_yyy_y[i] - 2.0 * tr_yyy_yyy[i] * tke_0 - 2.0 * tr_yyyy_yy[i] * tbe_0 - 6.0 * tr_xyy_xyy[i] * tke_0 - 4.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yy[i] * tbe_0 - 4.0 * tr_xxyyy_y[i] * tbe_0 + 4.0 * tr_xxyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_yz[i] = 3.0 * tr_yy_yz[i] + tr_yyy_z[i] - 2.0 * tr_yyy_yyz[i] * tke_0 - 2.0 * tr_yyyy_yz[i] * tbe_0 - 6.0 * tr_xyy_xyz[i] * tke_0 - 2.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yz[i] * tbe_0 - 2.0 * tr_xxyyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyy_zz[i] = 3.0 * tr_yy_zz[i] - 2.0 * tr_yyy_yzz[i] * tke_0 - 2.0 * tr_yyyy_zz[i] * tbe_0 - 6.0 * tr_xyy_xzz[i] * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyy_xzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 132-138 components of targeted buffer : GD

    auto tr_0_0_xy_xyyz_xx = pbuffer.data(idx_op_geom_020_gd + 132);

    auto tr_0_0_xy_xyyz_xy = pbuffer.data(idx_op_geom_020_gd + 133);

    auto tr_0_0_xy_xyyz_xz = pbuffer.data(idx_op_geom_020_gd + 134);

    auto tr_0_0_xy_xyyz_yy = pbuffer.data(idx_op_geom_020_gd + 135);

    auto tr_0_0_xy_xyyz_yz = pbuffer.data(idx_op_geom_020_gd + 136);

    auto tr_0_0_xy_xyyz_zz = pbuffer.data(idx_op_geom_020_gd + 137);

    #pragma omp simd aligned(tr_0_0_xy_xyyz_xx, tr_0_0_xy_xyyz_xy, tr_0_0_xy_xyyz_xz, tr_0_0_xy_xyyz_yy, tr_0_0_xy_xyyz_yz, tr_0_0_xy_xyyz_zz, tr_xxyyyz_xx, tr_xxyyyz_xy, tr_xxyyyz_xz, tr_xxyyyz_yy, tr_xxyyyz_yz, tr_xxyyyz_zz, tr_xxyyz_x, tr_xxyyz_xxy, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_y, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_z, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xyyyz_x, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_y, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xxxy, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyyz_xx[i] = 2.0 * tr_yz_xx[i] - 2.0 * tr_yyz_xxy[i] * tke_0 - 2.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xxx[i] * tke_0 - 4.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xy[i] = 2.0 * tr_yz_xy[i] + tr_yyz_x[i] - 2.0 * tr_yyz_xyy[i] * tke_0 - 2.0 * tr_yyyz_xy[i] * tbe_0 + 2.0 * tr_xyz_y[i] - 4.0 * tr_xyz_xxy[i] * tke_0 + tr_xyyz_0[i] - 2.0 * tr_xyyz_yy[i] * tke_0 - 2.0 * tr_xyyz_xx[i] * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 - 2.0 * tr_xxyyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_xz[i] = 2.0 * tr_yz_xz[i] - 2.0 * tr_yyz_xyz[i] * tke_0 - 2.0 * tr_yyyz_xz[i] * tbe_0 + 2.0 * tr_xyz_z[i] - 4.0 * tr_xyz_xxz[i] * tke_0 - 2.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_yy[i] = 2.0 * tr_yz_yy[i] + 2.0 * tr_yyz_y[i] - 2.0 * tr_yyz_yyy[i] * tke_0 - 2.0 * tr_yyyz_yy[i] * tbe_0 - 4.0 * tr_xyz_xyy[i] * tke_0 - 4.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yy[i] * tbe_0 - 4.0 * tr_xxyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_yz[i] = 2.0 * tr_yz_yz[i] + tr_yyz_z[i] - 2.0 * tr_yyz_yyz[i] * tke_0 - 2.0 * tr_yyyz_yz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tke_0 - 2.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 - 2.0 * tr_xxyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyyz_zz[i] = 2.0 * tr_yz_zz[i] - 2.0 * tr_yyz_yzz[i] * tke_0 - 2.0 * tr_yyyz_zz[i] * tbe_0 - 4.0 * tr_xyz_xzz[i] * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 138-144 components of targeted buffer : GD

    auto tr_0_0_xy_xyzz_xx = pbuffer.data(idx_op_geom_020_gd + 138);

    auto tr_0_0_xy_xyzz_xy = pbuffer.data(idx_op_geom_020_gd + 139);

    auto tr_0_0_xy_xyzz_xz = pbuffer.data(idx_op_geom_020_gd + 140);

    auto tr_0_0_xy_xyzz_yy = pbuffer.data(idx_op_geom_020_gd + 141);

    auto tr_0_0_xy_xyzz_yz = pbuffer.data(idx_op_geom_020_gd + 142);

    auto tr_0_0_xy_xyzz_zz = pbuffer.data(idx_op_geom_020_gd + 143);

    #pragma omp simd aligned(tr_0_0_xy_xyzz_xx, tr_0_0_xy_xyzz_xy, tr_0_0_xy_xyzz_xz, tr_0_0_xy_xyzz_yy, tr_0_0_xy_xyzz_yz, tr_0_0_xy_xyzz_zz, tr_xxyyzz_xx, tr_xxyyzz_xy, tr_xxyyzz_xz, tr_xxyyzz_yy, tr_xxyyzz_yz, tr_xxyyzz_zz, tr_xxyzz_x, tr_xxyzz_xxy, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_y, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_z, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xyyzz_x, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_y, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xxxy, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_z, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xyzz_xx[i] = tr_zz_xx[i] - 2.0 * tr_yzz_xxy[i] * tke_0 - 2.0 * tr_yyzz_xx[i] * tbe_0 + 2.0 * tr_xzz_x[i] - 2.0 * tr_xzz_xxx[i] * tke_0 - 4.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xy[i] = tr_zz_xy[i] + tr_yzz_x[i] - 2.0 * tr_yzz_xyy[i] * tke_0 - 2.0 * tr_yyzz_xy[i] * tbe_0 + tr_xzz_y[i] - 2.0 * tr_xzz_xxy[i] * tke_0 + tr_xyzz_0[i] - 2.0 * tr_xyzz_yy[i] * tke_0 - 2.0 * tr_xyzz_xx[i] * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xy[i] * tbe_0 - 2.0 * tr_xxyzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_xz[i] = tr_zz_xz[i] - 2.0 * tr_yzz_xyz[i] * tke_0 - 2.0 * tr_yyzz_xz[i] * tbe_0 + tr_xzz_z[i] - 2.0 * tr_xzz_xxz[i] * tke_0 - 2.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_yy[i] = tr_zz_yy[i] + 2.0 * tr_yzz_y[i] - 2.0 * tr_yzz_yyy[i] * tke_0 - 2.0 * tr_yyzz_yy[i] * tbe_0 - 2.0 * tr_xzz_xyy[i] * tke_0 - 4.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yy[i] * tbe_0 - 4.0 * tr_xxyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_yz[i] = tr_zz_yz[i] + tr_yzz_z[i] - 2.0 * tr_yzz_yyz[i] * tke_0 - 2.0 * tr_yyzz_yz[i] * tbe_0 - 2.0 * tr_xzz_xyz[i] * tke_0 - 2.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yz[i] * tbe_0 - 2.0 * tr_xxyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xyzz_zz[i] = tr_zz_zz[i] - 2.0 * tr_yzz_yzz[i] * tke_0 - 2.0 * tr_yyzz_zz[i] * tbe_0 - 2.0 * tr_xzz_xzz[i] * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 144-150 components of targeted buffer : GD

    auto tr_0_0_xy_xzzz_xx = pbuffer.data(idx_op_geom_020_gd + 144);

    auto tr_0_0_xy_xzzz_xy = pbuffer.data(idx_op_geom_020_gd + 145);

    auto tr_0_0_xy_xzzz_xz = pbuffer.data(idx_op_geom_020_gd + 146);

    auto tr_0_0_xy_xzzz_yy = pbuffer.data(idx_op_geom_020_gd + 147);

    auto tr_0_0_xy_xzzz_yz = pbuffer.data(idx_op_geom_020_gd + 148);

    auto tr_0_0_xy_xzzz_zz = pbuffer.data(idx_op_geom_020_gd + 149);

    #pragma omp simd aligned(tr_0_0_xy_xzzz_xx, tr_0_0_xy_xzzz_xy, tr_0_0_xy_xzzz_xz, tr_0_0_xy_xzzz_yy, tr_0_0_xy_xzzz_yz, tr_0_0_xy_xzzz_zz, tr_xxyzzz_xx, tr_xxyzzz_xy, tr_xxyzzz_xz, tr_xxyzzz_yy, tr_xxyzzz_yz, tr_xxyzzz_zz, tr_xxzzz_x, tr_xxzzz_xxy, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_y, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_z, tr_xyzzz_x, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_y, tr_xyzzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xxxy, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, tr_zzz_x, tr_zzz_xxy, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_xzzz_xx[i] = -2.0 * tr_zzz_xxy[i] * tke_0 - 2.0 * tr_yzzz_xx[i] * tbe_0 - 4.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xy[i] = tr_zzz_x[i] - 2.0 * tr_zzz_xyy[i] * tke_0 - 2.0 * tr_yzzz_xy[i] * tbe_0 + tr_xzzz_0[i] - 2.0 * tr_xzzz_yy[i] * tke_0 - 2.0 * tr_xzzz_xx[i] * tke_0 + 4.0 * tr_xzzz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_xz[i] = -2.0 * tr_zzz_xyz[i] * tke_0 - 2.0 * tr_yzzz_xz[i] * tbe_0 - 2.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_yy[i] = 2.0 * tr_zzz_y[i] - 2.0 * tr_zzz_yyy[i] * tke_0 - 2.0 * tr_yzzz_yy[i] * tbe_0 - 4.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxzzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_yz[i] = tr_zzz_z[i] - 2.0 * tr_zzz_yyz[i] * tke_0 - 2.0 * tr_yzzz_yz[i] * tbe_0 - 2.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzzz_z[i] * tbe_0 + 4.0 * tr_xxzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_xzzz_zz[i] = -2.0 * tr_zzz_yzz[i] * tke_0 - 2.0 * tr_yzzz_zz[i] * tbe_0 + 4.0 * tr_xzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-156 components of targeted buffer : GD

    auto tr_0_0_xy_yyyy_xx = pbuffer.data(idx_op_geom_020_gd + 150);

    auto tr_0_0_xy_yyyy_xy = pbuffer.data(idx_op_geom_020_gd + 151);

    auto tr_0_0_xy_yyyy_xz = pbuffer.data(idx_op_geom_020_gd + 152);

    auto tr_0_0_xy_yyyy_yy = pbuffer.data(idx_op_geom_020_gd + 153);

    auto tr_0_0_xy_yyyy_yz = pbuffer.data(idx_op_geom_020_gd + 154);

    auto tr_0_0_xy_yyyy_zz = pbuffer.data(idx_op_geom_020_gd + 155);

    #pragma omp simd aligned(tr_0_0_xy_yyyy_xx, tr_0_0_xy_yyyy_xy, tr_0_0_xy_yyyy_xz, tr_0_0_xy_yyyy_yy, tr_0_0_xy_yyyy_yz, tr_0_0_xy_yyyy_zz, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyy_x, tr_xyyyy_xxy, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_y, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_z, tr_xyyyyy_xx, tr_xyyyyy_xy, tr_xyyyyy_xz, tr_xyyyyy_yy, tr_xyyyyy_yz, tr_xyyyyy_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_z, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xxxy, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyyy_x, tr_yyyyy_xxx, tr_yyyyy_xxy, tr_yyyyy_xxz, tr_yyyyy_xyy, tr_yyyyy_xyz, tr_yyyyy_xzz, tr_yyyyy_y, tr_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyyy_xx[i] = 8.0 * tr_yyy_x[i] - 8.0 * tr_yyy_xxx[i] * tke_0 - 4.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_x[i] * tbe_0 + 4.0 * tr_yyyyy_xxx[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xy[i] = 4.0 * tr_yyy_y[i] - 8.0 * tr_yyy_xxy[i] * tke_0 + tr_yyyy_0[i] - 2.0 * tr_yyyy_yy[i] * tke_0 - 2.0 * tr_yyyy_xx[i] * tke_0 + 4.0 * tr_yyyy_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyyy_y[i] * tbe_0 + 4.0 * tr_yyyyy_xxy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xy[i] * tbe_0 - 2.0 * tr_xyyyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_xz[i] = 4.0 * tr_yyy_z[i] - 8.0 * tr_yyy_xxz[i] * tke_0 - 2.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyy_z[i] * tbe_0 + 4.0 * tr_yyyyy_xxz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_yy[i] = -8.0 * tr_yyy_xyy[i] * tke_0 - 4.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xyy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yy[i] * tbe_0 - 4.0 * tr_xyyyy_y[i] * tbe_0 + 4.0 * tr_xyyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_yz[i] = -8.0 * tr_yyy_xyz[i] * tke_0 - 2.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xyz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yz[i] * tbe_0 - 2.0 * tr_xyyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyy_zz[i] = -8.0 * tr_yyy_xzz[i] * tke_0 + 4.0 * tr_yyyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyy_xzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 156-162 components of targeted buffer : GD

    auto tr_0_0_xy_yyyz_xx = pbuffer.data(idx_op_geom_020_gd + 156);

    auto tr_0_0_xy_yyyz_xy = pbuffer.data(idx_op_geom_020_gd + 157);

    auto tr_0_0_xy_yyyz_xz = pbuffer.data(idx_op_geom_020_gd + 158);

    auto tr_0_0_xy_yyyz_yy = pbuffer.data(idx_op_geom_020_gd + 159);

    auto tr_0_0_xy_yyyz_yz = pbuffer.data(idx_op_geom_020_gd + 160);

    auto tr_0_0_xy_yyyz_zz = pbuffer.data(idx_op_geom_020_gd + 161);

    #pragma omp simd aligned(tr_0_0_xy_yyyz_xx, tr_0_0_xy_yyyz_xy, tr_0_0_xy_yyyz_xz, tr_0_0_xy_yyyz_yy, tr_0_0_xy_yyyz_yz, tr_0_0_xy_yyyz_zz, tr_xyyyyz_xx, tr_xyyyyz_xy, tr_xyyyyz_xz, tr_xyyyyz_yy, tr_xyyyyz_yz, tr_xyyyyz_zz, tr_xyyyz_x, tr_xyyyz_xxy, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_y, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_z, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yyyyz_x, tr_yyyyz_xxx, tr_yyyyz_xxy, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_y, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xxxy, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyyz_xx[i] = 6.0 * tr_yyz_x[i] - 6.0 * tr_yyz_xxx[i] * tke_0 - 4.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_x[i] * tbe_0 + 4.0 * tr_yyyyz_xxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xy[i] = 3.0 * tr_yyz_y[i] - 6.0 * tr_yyz_xxy[i] * tke_0 + tr_yyyz_0[i] - 2.0 * tr_yyyz_yy[i] * tke_0 - 2.0 * tr_yyyz_xx[i] * tke_0 + 4.0 * tr_yyyz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_y[i] * tbe_0 + 4.0 * tr_yyyyz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xy[i] * tbe_0 - 2.0 * tr_xyyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_xz[i] = 3.0 * tr_yyz_z[i] - 6.0 * tr_yyz_xxz[i] * tke_0 - 2.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_z[i] * tbe_0 + 4.0 * tr_yyyyz_xxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_yy[i] = -6.0 * tr_yyz_xyy[i] * tke_0 - 4.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yy[i] * tbe_0 - 4.0 * tr_xyyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_yz[i] = -6.0 * tr_yyz_xyz[i] * tke_0 - 2.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yz[i] * tbe_0 - 2.0 * tr_xyyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyyz_zz[i] = -6.0 * tr_yyz_xzz[i] * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 162-168 components of targeted buffer : GD

    auto tr_0_0_xy_yyzz_xx = pbuffer.data(idx_op_geom_020_gd + 162);

    auto tr_0_0_xy_yyzz_xy = pbuffer.data(idx_op_geom_020_gd + 163);

    auto tr_0_0_xy_yyzz_xz = pbuffer.data(idx_op_geom_020_gd + 164);

    auto tr_0_0_xy_yyzz_yy = pbuffer.data(idx_op_geom_020_gd + 165);

    auto tr_0_0_xy_yyzz_yz = pbuffer.data(idx_op_geom_020_gd + 166);

    auto tr_0_0_xy_yyzz_zz = pbuffer.data(idx_op_geom_020_gd + 167);

    #pragma omp simd aligned(tr_0_0_xy_yyzz_xx, tr_0_0_xy_yyzz_xy, tr_0_0_xy_yyzz_xz, tr_0_0_xy_yyzz_yy, tr_0_0_xy_yyzz_yz, tr_0_0_xy_yyzz_zz, tr_xyyyzz_xx, tr_xyyyzz_xy, tr_xyyyzz_xz, tr_xyyyzz_yy, tr_xyyyzz_yz, tr_xyyyzz_zz, tr_xyyzz_x, tr_xyyzz_xxy, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_y, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_z, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_yyyzz_x, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_y, tr_yyyzz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xxxy, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yyzz_xx[i] = 4.0 * tr_yzz_x[i] - 4.0 * tr_yzz_xxx[i] * tke_0 - 4.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xy[i] = 2.0 * tr_yzz_y[i] - 4.0 * tr_yzz_xxy[i] * tke_0 + tr_yyzz_0[i] - 2.0 * tr_yyzz_yy[i] * tke_0 - 2.0 * tr_yyzz_xx[i] * tke_0 + 4.0 * tr_yyzz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_y[i] * tbe_0 + 4.0 * tr_yyyzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 - 2.0 * tr_xyyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_xz[i] = 2.0 * tr_yzz_z[i] - 4.0 * tr_yzz_xxz[i] * tke_0 - 2.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_z[i] * tbe_0 + 4.0 * tr_yyyzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_yy[i] = -4.0 * tr_yzz_xyy[i] * tke_0 - 4.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yy[i] * tbe_0 - 4.0 * tr_xyyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_yz[i] = -4.0 * tr_yzz_xyz[i] * tke_0 - 2.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 - 2.0 * tr_xyyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yyzz_zz[i] = -4.0 * tr_yzz_xzz[i] * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 168-174 components of targeted buffer : GD

    auto tr_0_0_xy_yzzz_xx = pbuffer.data(idx_op_geom_020_gd + 168);

    auto tr_0_0_xy_yzzz_xy = pbuffer.data(idx_op_geom_020_gd + 169);

    auto tr_0_0_xy_yzzz_xz = pbuffer.data(idx_op_geom_020_gd + 170);

    auto tr_0_0_xy_yzzz_yy = pbuffer.data(idx_op_geom_020_gd + 171);

    auto tr_0_0_xy_yzzz_yz = pbuffer.data(idx_op_geom_020_gd + 172);

    auto tr_0_0_xy_yzzz_zz = pbuffer.data(idx_op_geom_020_gd + 173);

    #pragma omp simd aligned(tr_0_0_xy_yzzz_xx, tr_0_0_xy_yzzz_xy, tr_0_0_xy_yzzz_xz, tr_0_0_xy_yzzz_yy, tr_0_0_xy_yzzz_yz, tr_0_0_xy_yzzz_zz, tr_xyyzzz_xx, tr_xyyzzz_xy, tr_xyyzzz_xz, tr_xyyzzz_yy, tr_xyyzzz_yz, tr_xyyzzz_zz, tr_xyzzz_x, tr_xyzzz_xxy, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_y, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_z, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_yyzzz_x, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_y, tr_yyzzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xxxy, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_yzzz_xx[i] = 2.0 * tr_zzz_x[i] - 2.0 * tr_zzz_xxx[i] * tke_0 - 4.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xy[i] = tr_zzz_y[i] - 2.0 * tr_zzz_xxy[i] * tke_0 + tr_yzzz_0[i] - 2.0 * tr_yzzz_yy[i] * tke_0 - 2.0 * tr_yzzz_xx[i] * tke_0 + 4.0 * tr_yzzz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_y[i] * tbe_0 + 4.0 * tr_yyzzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xy[i] * tbe_0 - 2.0 * tr_xyzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_xz[i] = tr_zzz_z[i] - 2.0 * tr_zzz_xxz[i] * tke_0 - 2.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_z[i] * tbe_0 + 4.0 * tr_yyzzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_yy[i] = -2.0 * tr_zzz_xyy[i] * tke_0 - 4.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yy[i] * tbe_0 - 4.0 * tr_xyzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_yz[i] = -2.0 * tr_zzz_xyz[i] * tke_0 - 2.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yz[i] * tbe_0 - 2.0 * tr_xyzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_yzzz_zz[i] = -2.0 * tr_zzz_xzz[i] * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 174-180 components of targeted buffer : GD

    auto tr_0_0_xy_zzzz_xx = pbuffer.data(idx_op_geom_020_gd + 174);

    auto tr_0_0_xy_zzzz_xy = pbuffer.data(idx_op_geom_020_gd + 175);

    auto tr_0_0_xy_zzzz_xz = pbuffer.data(idx_op_geom_020_gd + 176);

    auto tr_0_0_xy_zzzz_yy = pbuffer.data(idx_op_geom_020_gd + 177);

    auto tr_0_0_xy_zzzz_yz = pbuffer.data(idx_op_geom_020_gd + 178);

    auto tr_0_0_xy_zzzz_zz = pbuffer.data(idx_op_geom_020_gd + 179);

    #pragma omp simd aligned(tr_0_0_xy_zzzz_xx, tr_0_0_xy_zzzz_xy, tr_0_0_xy_zzzz_xz, tr_0_0_xy_zzzz_yy, tr_0_0_xy_zzzz_yz, tr_0_0_xy_zzzz_zz, tr_xyzzzz_xx, tr_xyzzzz_xy, tr_xyzzzz_xz, tr_xyzzzz_yy, tr_xyzzzz_yz, tr_xyzzzz_zz, tr_xzzzz_x, tr_xzzzz_xxy, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_y, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_z, tr_yzzzz_x, tr_yzzzz_xxx, tr_yzzzz_xxy, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_y, tr_yzzzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xxxy, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xy_zzzz_xx[i] = -4.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xxxy[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xy[i] = tr_zzzz_0[i] - 2.0 * tr_zzzz_yy[i] * tke_0 - 2.0 * tr_zzzz_xx[i] * tke_0 + 4.0 * tr_zzzz_xxyy[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_xz[i] = -2.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_z[i] * tbe_0 + 4.0 * tr_yzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_yy[i] = -4.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xyyy[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xzzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_yz[i] = -2.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzzz_z[i] * tbe_0 + 4.0 * tr_xzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xy_zzzz_zz[i] = 4.0 * tr_zzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-186 components of targeted buffer : GD

    auto tr_0_0_xz_xxxx_xx = pbuffer.data(idx_op_geom_020_gd + 180);

    auto tr_0_0_xz_xxxx_xy = pbuffer.data(idx_op_geom_020_gd + 181);

    auto tr_0_0_xz_xxxx_xz = pbuffer.data(idx_op_geom_020_gd + 182);

    auto tr_0_0_xz_xxxx_yy = pbuffer.data(idx_op_geom_020_gd + 183);

    auto tr_0_0_xz_xxxx_yz = pbuffer.data(idx_op_geom_020_gd + 184);

    auto tr_0_0_xz_xxxx_zz = pbuffer.data(idx_op_geom_020_gd + 185);

    #pragma omp simd aligned(tr_0_0_xz_xxxx_xx, tr_0_0_xz_xxxx_xy, tr_0_0_xz_xxxx_xz, tr_0_0_xz_xxxx_yy, tr_0_0_xz_xxxx_yz, tr_0_0_xz_xxxx_zz, tr_xxx_x, tr_xxx_xxz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xxxz, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxx_x, tr_xxxxx_xxz, tr_xxxxx_xyz, tr_xxxxx_xzz, tr_xxxxx_y, tr_xxxxx_yyz, tr_xxxxx_yzz, tr_xxxxx_z, tr_xxxxx_zzz, tr_xxxxxz_xx, tr_xxxxxz_xy, tr_xxxxxz_xz, tr_xxxxxz_yy, tr_xxxxxz_yz, tr_xxxxxz_zz, tr_xxxxz_x, tr_xxxxz_xxx, tr_xxxxz_xxy, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_y, tr_xxxxz_z, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxxx_xx[i] = -8.0 * tr_xxx_xxz[i] * tke_0 - 8.0 * tr_xxxz_xx[i] * tbe_0 - 4.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_x[i] * tbe_0 + 4.0 * tr_xxxxz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xy[i] = -8.0 * tr_xxx_xyz[i] * tke_0 - 8.0 * tr_xxxz_xy[i] * tbe_0 - 2.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_y[i] * tbe_0 + 4.0 * tr_xxxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_xz[i] = 4.0 * tr_xxx_x[i] - 8.0 * tr_xxx_xzz[i] * tke_0 - 8.0 * tr_xxxz_xz[i] * tbe_0 + tr_xxxx_0[i] - 2.0 * tr_xxxx_zz[i] * tke_0 - 2.0 * tr_xxxx_xx[i] * tke_0 + 4.0 * tr_xxxx_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_z[i] * tbe_0 + 4.0 * tr_xxxxz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_x[i] * tbe_0 + 4.0 * tr_xxxxx_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_yy[i] = -8.0 * tr_xxx_yyz[i] * tke_0 - 8.0 * tr_xxxz_yy[i] * tbe_0 + 4.0 * tr_xxxx_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_yz[i] = 4.0 * tr_xxx_y[i] - 8.0 * tr_xxx_yzz[i] * tke_0 - 8.0 * tr_xxxz_yz[i] * tbe_0 - 2.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxx_y[i] * tbe_0 + 4.0 * tr_xxxxx_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxx_zz[i] = 8.0 * tr_xxx_z[i] - 8.0 * tr_xxx_zzz[i] * tke_0 - 8.0 * tr_xxxz_zz[i] * tbe_0 - 4.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxx_z[i] * tbe_0 + 4.0 * tr_xxxxx_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 186-192 components of targeted buffer : GD

    auto tr_0_0_xz_xxxy_xx = pbuffer.data(idx_op_geom_020_gd + 186);

    auto tr_0_0_xz_xxxy_xy = pbuffer.data(idx_op_geom_020_gd + 187);

    auto tr_0_0_xz_xxxy_xz = pbuffer.data(idx_op_geom_020_gd + 188);

    auto tr_0_0_xz_xxxy_yy = pbuffer.data(idx_op_geom_020_gd + 189);

    auto tr_0_0_xz_xxxy_yz = pbuffer.data(idx_op_geom_020_gd + 190);

    auto tr_0_0_xz_xxxy_zz = pbuffer.data(idx_op_geom_020_gd + 191);

    #pragma omp simd aligned(tr_0_0_xz_xxxy_xx, tr_0_0_xz_xxxy_xy, tr_0_0_xz_xxxy_xz, tr_0_0_xz_xxxy_yy, tr_0_0_xz_xxxy_yz, tr_0_0_xz_xxxy_zz, tr_xxxxy_x, tr_xxxxy_xxz, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_y, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_z, tr_xxxxy_zzz, tr_xxxxyz_xx, tr_xxxxyz_xy, tr_xxxxyz_xz, tr_xxxxyz_yy, tr_xxxxyz_yz, tr_xxxxyz_zz, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xxxz, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_y, tr_xxxyz_z, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxxy_xx[i] = -6.0 * tr_xxy_xxz[i] * tke_0 - 6.0 * tr_xxyz_xx[i] * tbe_0 - 4.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xy[i] = -6.0 * tr_xxy_xyz[i] * tke_0 - 6.0 * tr_xxyz_xy[i] * tbe_0 - 2.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_xz[i] = 3.0 * tr_xxy_x[i] - 6.0 * tr_xxy_xzz[i] * tke_0 - 6.0 * tr_xxyz_xz[i] * tbe_0 + tr_xxxy_0[i] - 2.0 * tr_xxxy_zz[i] * tke_0 - 2.0 * tr_xxxy_xx[i] * tke_0 + 4.0 * tr_xxxy_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_x[i] * tbe_0 + 4.0 * tr_xxxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_yy[i] = -6.0 * tr_xxy_yyz[i] * tke_0 - 6.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxxy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_yz[i] = 3.0 * tr_xxy_y[i] - 6.0 * tr_xxy_yzz[i] * tke_0 - 6.0 * tr_xxyz_yz[i] * tbe_0 - 2.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_y[i] * tbe_0 + 4.0 * tr_xxxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxy_zz[i] = 6.0 * tr_xxy_z[i] - 6.0 * tr_xxy_zzz[i] * tke_0 - 6.0 * tr_xxyz_zz[i] * tbe_0 - 4.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_z[i] * tbe_0 + 4.0 * tr_xxxxy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 192-198 components of targeted buffer : GD

    auto tr_0_0_xz_xxxz_xx = pbuffer.data(idx_op_geom_020_gd + 192);

    auto tr_0_0_xz_xxxz_xy = pbuffer.data(idx_op_geom_020_gd + 193);

    auto tr_0_0_xz_xxxz_xz = pbuffer.data(idx_op_geom_020_gd + 194);

    auto tr_0_0_xz_xxxz_yy = pbuffer.data(idx_op_geom_020_gd + 195);

    auto tr_0_0_xz_xxxz_yz = pbuffer.data(idx_op_geom_020_gd + 196);

    auto tr_0_0_xz_xxxz_zz = pbuffer.data(idx_op_geom_020_gd + 197);

    #pragma omp simd aligned(tr_0_0_xz_xxxz_xx, tr_0_0_xz_xxxz_xy, tr_0_0_xz_xxxz_xz, tr_0_0_xz_xxxz_yy, tr_0_0_xz_xxxz_yz, tr_0_0_xz_xxxz_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxx_x, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_z, tr_xxxx_xx, tr_xxxx_xy, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxz_x, tr_xxxxz_xxz, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_y, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_z, tr_xxxxz_zzz, tr_xxxxzz_xx, tr_xxxxzz_xy, tr_xxxxzz_xz, tr_xxxxzz_yy, tr_xxxxzz_yz, tr_xxxxzz_zz, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xxxz, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_y, tr_xxxzz_z, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxxz_xx[i] = 3.0 * tr_xx_xx[i] - 6.0 * tr_xxz_xxz[i] * tke_0 - 6.0 * tr_xxzz_xx[i] * tbe_0 + 2.0 * tr_xxx_x[i] - 2.0 * tr_xxx_xxx[i] * tke_0 - 4.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xy[i] = 3.0 * tr_xx_xy[i] - 6.0 * tr_xxz_xyz[i] * tke_0 - 6.0 * tr_xxzz_xy[i] * tbe_0 + tr_xxx_y[i] - 2.0 * tr_xxx_xxy[i] * tke_0 - 2.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_y[i] * tbe_0 + 4.0 * tr_xxxzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_xz[i] = 3.0 * tr_xx_xz[i] + 3.0 * tr_xxz_x[i] - 6.0 * tr_xxz_xzz[i] * tke_0 - 6.0 * tr_xxzz_xz[i] * tbe_0 + tr_xxx_z[i] - 2.0 * tr_xxx_xxz[i] * tke_0 + tr_xxxz_0[i] - 2.0 * tr_xxxz_zz[i] * tke_0 - 2.0 * tr_xxxz_xx[i] * tke_0 + 4.0 * tr_xxxz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_z[i] * tbe_0 + 4.0 * tr_xxxzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xz[i] * tbe_0 - 2.0 * tr_xxxxz_x[i] * tbe_0 + 4.0 * tr_xxxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_yy[i] = 3.0 * tr_xx_yy[i] - 6.0 * tr_xxz_yyz[i] * tke_0 - 6.0 * tr_xxzz_yy[i] * tbe_0 - 2.0 * tr_xxx_xyy[i] * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yy[i] * tbe_0 + 4.0 * tr_xxxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_yz[i] = 3.0 * tr_xx_yz[i] + 3.0 * tr_xxz_y[i] - 6.0 * tr_xxz_yzz[i] * tke_0 - 6.0 * tr_xxzz_yz[i] * tbe_0 - 2.0 * tr_xxx_xyz[i] * tke_0 - 2.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yz[i] * tbe_0 - 2.0 * tr_xxxxz_y[i] * tbe_0 + 4.0 * tr_xxxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxxz_zz[i] = 3.0 * tr_xx_zz[i] + 6.0 * tr_xxz_z[i] - 6.0 * tr_xxz_zzz[i] * tke_0 - 6.0 * tr_xxzz_zz[i] * tbe_0 - 2.0 * tr_xxx_xzz[i] * tke_0 - 4.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zz[i] * tbe_0 - 4.0 * tr_xxxxz_z[i] * tbe_0 + 4.0 * tr_xxxxz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 198-204 components of targeted buffer : GD

    auto tr_0_0_xz_xxyy_xx = pbuffer.data(idx_op_geom_020_gd + 198);

    auto tr_0_0_xz_xxyy_xy = pbuffer.data(idx_op_geom_020_gd + 199);

    auto tr_0_0_xz_xxyy_xz = pbuffer.data(idx_op_geom_020_gd + 200);

    auto tr_0_0_xz_xxyy_yy = pbuffer.data(idx_op_geom_020_gd + 201);

    auto tr_0_0_xz_xxyy_yz = pbuffer.data(idx_op_geom_020_gd + 202);

    auto tr_0_0_xz_xxyy_zz = pbuffer.data(idx_op_geom_020_gd + 203);

    #pragma omp simd aligned(tr_0_0_xz_xxyy_xx, tr_0_0_xz_xxyy_xy, tr_0_0_xz_xxyy_xz, tr_0_0_xz_xxyy_yy, tr_0_0_xz_xxyy_yz, tr_0_0_xz_xxyy_zz, tr_xxxyy_x, tr_xxxyy_xxz, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_y, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_z, tr_xxxyy_zzz, tr_xxxyyz_xx, tr_xxxyyz_xy, tr_xxxyyz_xz, tr_xxxyyz_yy, tr_xxxyyz_yz, tr_xxxyyz_zz, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xxxz, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_y, tr_xxyyz_z, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxyy_xx[i] = -4.0 * tr_xyy_xxz[i] * tke_0 - 4.0 * tr_xyyz_xx[i] * tbe_0 - 4.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xy[i] = -4.0 * tr_xyy_xyz[i] * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 - 2.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_xz[i] = 2.0 * tr_xyy_x[i] - 4.0 * tr_xyy_xzz[i] * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 + tr_xxyy_0[i] - 2.0 * tr_xxyy_zz[i] * tke_0 - 2.0 * tr_xxyy_xx[i] * tke_0 + 4.0 * tr_xxyy_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_x[i] * tbe_0 + 4.0 * tr_xxxyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_yy[i] = -4.0 * tr_xyy_yyz[i] * tke_0 - 4.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xxyy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_yz[i] = 2.0 * tr_xyy_y[i] - 4.0 * tr_xyy_yzz[i] * tke_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 - 2.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_y[i] * tbe_0 + 4.0 * tr_xxxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyy_zz[i] = 4.0 * tr_xyy_z[i] - 4.0 * tr_xyy_zzz[i] * tke_0 - 4.0 * tr_xyyz_zz[i] * tbe_0 - 4.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_z[i] * tbe_0 + 4.0 * tr_xxxyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 204-210 components of targeted buffer : GD

    auto tr_0_0_xz_xxyz_xx = pbuffer.data(idx_op_geom_020_gd + 204);

    auto tr_0_0_xz_xxyz_xy = pbuffer.data(idx_op_geom_020_gd + 205);

    auto tr_0_0_xz_xxyz_xz = pbuffer.data(idx_op_geom_020_gd + 206);

    auto tr_0_0_xz_xxyz_yy = pbuffer.data(idx_op_geom_020_gd + 207);

    auto tr_0_0_xz_xxyz_yz = pbuffer.data(idx_op_geom_020_gd + 208);

    auto tr_0_0_xz_xxyz_zz = pbuffer.data(idx_op_geom_020_gd + 209);

    #pragma omp simd aligned(tr_0_0_xz_xxyz_xx, tr_0_0_xz_xxyz_xy, tr_0_0_xz_xxyz_xz, tr_0_0_xz_xxyz_yy, tr_0_0_xz_xxyz_yz, tr_0_0_xz_xxyz_zz, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_xxz, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_y, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_z, tr_xxxyz_zzz, tr_xxxyzz_xx, tr_xxxyzz_xy, tr_xxxyzz_xz, tr_xxxyzz_yy, tr_xxxyzz_yz, tr_xxxyzz_zz, tr_xxy_x, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xxxz, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_y, tr_xxyzz_z, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxyz_xx[i] = 2.0 * tr_xy_xx[i] - 4.0 * tr_xyz_xxz[i] * tke_0 - 4.0 * tr_xyzz_xx[i] * tbe_0 + 2.0 * tr_xxy_x[i] - 2.0 * tr_xxy_xxx[i] * tke_0 - 4.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xy[i] = 2.0 * tr_xy_xy[i] - 4.0 * tr_xyz_xyz[i] * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + tr_xxy_y[i] - 2.0 * tr_xxy_xxy[i] * tke_0 - 2.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_xz[i] = 2.0 * tr_xy_xz[i] + 2.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xzz[i] * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 + tr_xxy_z[i] - 2.0 * tr_xxy_xxz[i] * tke_0 + tr_xxyz_0[i] - 2.0 * tr_xxyz_zz[i] * tke_0 - 2.0 * tr_xxyz_xx[i] * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xz[i] * tbe_0 - 2.0 * tr_xxxyz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_yy[i] = 2.0 * tr_xy_yy[i] - 4.0 * tr_xyz_yyz[i] * tke_0 - 4.0 * tr_xyzz_yy[i] * tbe_0 - 2.0 * tr_xxy_xyy[i] * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_yz[i] = 2.0 * tr_xy_yz[i] + 2.0 * tr_xyz_y[i] - 4.0 * tr_xyz_yzz[i] * tke_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 - 2.0 * tr_xxy_xyz[i] * tke_0 - 2.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yz[i] * tbe_0 - 2.0 * tr_xxxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxyz_zz[i] = 2.0 * tr_xy_zz[i] + 4.0 * tr_xyz_z[i] - 4.0 * tr_xyz_zzz[i] * tke_0 - 4.0 * tr_xyzz_zz[i] * tbe_0 - 2.0 * tr_xxy_xzz[i] * tke_0 - 4.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_zz[i] * tbe_0 - 4.0 * tr_xxxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-216 components of targeted buffer : GD

    auto tr_0_0_xz_xxzz_xx = pbuffer.data(idx_op_geom_020_gd + 210);

    auto tr_0_0_xz_xxzz_xy = pbuffer.data(idx_op_geom_020_gd + 211);

    auto tr_0_0_xz_xxzz_xz = pbuffer.data(idx_op_geom_020_gd + 212);

    auto tr_0_0_xz_xxzz_yy = pbuffer.data(idx_op_geom_020_gd + 213);

    auto tr_0_0_xz_xxzz_yz = pbuffer.data(idx_op_geom_020_gd + 214);

    auto tr_0_0_xz_xxzz_zz = pbuffer.data(idx_op_geom_020_gd + 215);

    #pragma omp simd aligned(tr_0_0_xz_xxzz_xx, tr_0_0_xz_xxzz_xy, tr_0_0_xz_xxzz_xz, tr_0_0_xz_xxzz_yy, tr_0_0_xz_xxzz_yz, tr_0_0_xz_xxzz_zz, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_xxz, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_y, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_z, tr_xxxzz_zzz, tr_xxxzzz_xx, tr_xxxzzz_xy, tr_xxxzzz_xz, tr_xxxzzz_yy, tr_xxxzzz_yz, tr_xxxzzz_zz, tr_xxz_x, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xxxz, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_y, tr_xxzzz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xxzz_xx[i] = 4.0 * tr_xz_xx[i] - 4.0 * tr_xzz_xxz[i] * tke_0 - 4.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xxz_x[i] - 4.0 * tr_xxz_xxx[i] * tke_0 - 4.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xy[i] = 4.0 * tr_xz_xy[i] - 4.0 * tr_xzz_xyz[i] * tke_0 - 4.0 * tr_xzzz_xy[i] * tbe_0 + 2.0 * tr_xxz_y[i] - 4.0 * tr_xxz_xxy[i] * tke_0 - 2.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_xz[i] = 4.0 * tr_xz_xz[i] + 2.0 * tr_xzz_x[i] - 4.0 * tr_xzz_xzz[i] * tke_0 - 4.0 * tr_xzzz_xz[i] * tbe_0 + 2.0 * tr_xxz_z[i] - 4.0 * tr_xxz_xxz[i] * tke_0 + tr_xxzz_0[i] - 2.0 * tr_xxzz_zz[i] * tke_0 - 2.0 * tr_xxzz_xx[i] * tke_0 + 4.0 * tr_xxzz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_z[i] * tbe_0 + 4.0 * tr_xxzzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xz[i] * tbe_0 - 2.0 * tr_xxxzz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_yy[i] = 4.0 * tr_xz_yy[i] - 4.0 * tr_xzz_yyz[i] * tke_0 - 4.0 * tr_xzzz_yy[i] * tbe_0 - 4.0 * tr_xxz_xyy[i] * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_yz[i] = 4.0 * tr_xz_yz[i] + 2.0 * tr_xzz_y[i] - 4.0 * tr_xzz_yzz[i] * tke_0 - 4.0 * tr_xzzz_yz[i] * tbe_0 - 4.0 * tr_xxz_xyz[i] * tke_0 - 2.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yz[i] * tbe_0 - 2.0 * tr_xxxzz_y[i] * tbe_0 + 4.0 * tr_xxxzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xxzz_zz[i] = 4.0 * tr_xz_zz[i] + 4.0 * tr_xzz_z[i] - 4.0 * tr_xzz_zzz[i] * tke_0 - 4.0 * tr_xzzz_zz[i] * tbe_0 - 4.0 * tr_xxz_xzz[i] * tke_0 - 4.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_zz[i] * tbe_0 - 4.0 * tr_xxxzz_z[i] * tbe_0 + 4.0 * tr_xxxzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 216-222 components of targeted buffer : GD

    auto tr_0_0_xz_xyyy_xx = pbuffer.data(idx_op_geom_020_gd + 216);

    auto tr_0_0_xz_xyyy_xy = pbuffer.data(idx_op_geom_020_gd + 217);

    auto tr_0_0_xz_xyyy_xz = pbuffer.data(idx_op_geom_020_gd + 218);

    auto tr_0_0_xz_xyyy_yy = pbuffer.data(idx_op_geom_020_gd + 219);

    auto tr_0_0_xz_xyyy_yz = pbuffer.data(idx_op_geom_020_gd + 220);

    auto tr_0_0_xz_xyyy_zz = pbuffer.data(idx_op_geom_020_gd + 221);

    #pragma omp simd aligned(tr_0_0_xz_xyyy_xx, tr_0_0_xz_xyyy_xy, tr_0_0_xz_xyyy_xz, tr_0_0_xz_xyyy_yy, tr_0_0_xz_xyyy_yz, tr_0_0_xz_xyyy_zz, tr_xxyyy_x, tr_xxyyy_xxz, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_y, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_z, tr_xxyyy_zzz, tr_xxyyyz_xx, tr_xxyyyz_xy, tr_xxyyyz_xz, tr_xxyyyz_yy, tr_xxyyyz_yz, tr_xxyyyz_zz, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xxxz, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_y, tr_xyyyz_z, tr_yyy_x, tr_yyy_xxz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyyy_xx[i] = -2.0 * tr_yyy_xxz[i] * tke_0 - 2.0 * tr_yyyz_xx[i] * tbe_0 - 4.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xy[i] = -2.0 * tr_yyy_xyz[i] * tke_0 - 2.0 * tr_yyyz_xy[i] * tbe_0 - 2.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_xz[i] = tr_yyy_x[i] - 2.0 * tr_yyy_xzz[i] * tke_0 - 2.0 * tr_yyyz_xz[i] * tbe_0 + tr_xyyy_0[i] - 2.0 * tr_xyyy_zz[i] * tke_0 - 2.0 * tr_xyyy_xx[i] * tke_0 + 4.0 * tr_xyyy_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_x[i] * tbe_0 + 4.0 * tr_xxyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_yy[i] = -2.0 * tr_yyy_yyz[i] * tke_0 - 2.0 * tr_yyyz_yy[i] * tbe_0 + 4.0 * tr_xyyy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_yz[i] = tr_yyy_y[i] - 2.0 * tr_yyy_yzz[i] * tke_0 - 2.0 * tr_yyyz_yz[i] * tbe_0 - 2.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_y[i] * tbe_0 + 4.0 * tr_xxyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyy_zz[i] = 2.0 * tr_yyy_z[i] - 2.0 * tr_yyy_zzz[i] * tke_0 - 2.0 * tr_yyyz_zz[i] * tbe_0 - 4.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 222-228 components of targeted buffer : GD

    auto tr_0_0_xz_xyyz_xx = pbuffer.data(idx_op_geom_020_gd + 222);

    auto tr_0_0_xz_xyyz_xy = pbuffer.data(idx_op_geom_020_gd + 223);

    auto tr_0_0_xz_xyyz_xz = pbuffer.data(idx_op_geom_020_gd + 224);

    auto tr_0_0_xz_xyyz_yy = pbuffer.data(idx_op_geom_020_gd + 225);

    auto tr_0_0_xz_xyyz_yz = pbuffer.data(idx_op_geom_020_gd + 226);

    auto tr_0_0_xz_xyyz_zz = pbuffer.data(idx_op_geom_020_gd + 227);

    #pragma omp simd aligned(tr_0_0_xz_xyyz_xx, tr_0_0_xz_xyyz_xy, tr_0_0_xz_xyyz_xz, tr_0_0_xz_xyyz_yy, tr_0_0_xz_xyyz_yz, tr_0_0_xz_xyyz_zz, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_xxz, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_y, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_z, tr_xxyyz_zzz, tr_xxyyzz_xx, tr_xxyyzz_xy, tr_xxyyzz_xz, tr_xxyyzz_yy, tr_xxyyzz_yz, tr_xxyyzz_zz, tr_xyy_x, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xxxz, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_y, tr_xyyzz_z, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyyz_xx[i] = tr_yy_xx[i] - 2.0 * tr_yyz_xxz[i] * tke_0 - 2.0 * tr_yyzz_xx[i] * tbe_0 + 2.0 * tr_xyy_x[i] - 2.0 * tr_xyy_xxx[i] * tke_0 - 4.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xy[i] = tr_yy_xy[i] - 2.0 * tr_yyz_xyz[i] * tke_0 - 2.0 * tr_yyzz_xy[i] * tbe_0 + tr_xyy_y[i] - 2.0 * tr_xyy_xxy[i] * tke_0 - 2.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_xz[i] = tr_yy_xz[i] + tr_yyz_x[i] - 2.0 * tr_yyz_xzz[i] * tke_0 - 2.0 * tr_yyzz_xz[i] * tbe_0 + tr_xyy_z[i] - 2.0 * tr_xyy_xxz[i] * tke_0 + tr_xyyz_0[i] - 2.0 * tr_xyyz_zz[i] * tke_0 - 2.0 * tr_xyyz_xx[i] * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xz[i] * tbe_0 - 2.0 * tr_xxyyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_yy[i] = tr_yy_yy[i] - 2.0 * tr_yyz_yyz[i] * tke_0 - 2.0 * tr_yyzz_yy[i] * tbe_0 - 2.0 * tr_xyy_xyy[i] * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_yz[i] = tr_yy_yz[i] + tr_yyz_y[i] - 2.0 * tr_yyz_yzz[i] * tke_0 - 2.0 * tr_yyzz_yz[i] * tbe_0 - 2.0 * tr_xyy_xyz[i] * tke_0 - 2.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yz[i] * tbe_0 - 2.0 * tr_xxyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyyz_zz[i] = tr_yy_zz[i] + 2.0 * tr_yyz_z[i] - 2.0 * tr_yyz_zzz[i] * tke_0 - 2.0 * tr_yyzz_zz[i] * tbe_0 - 2.0 * tr_xyy_xzz[i] * tke_0 - 4.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_zz[i] * tbe_0 - 4.0 * tr_xxyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 228-234 components of targeted buffer : GD

    auto tr_0_0_xz_xyzz_xx = pbuffer.data(idx_op_geom_020_gd + 228);

    auto tr_0_0_xz_xyzz_xy = pbuffer.data(idx_op_geom_020_gd + 229);

    auto tr_0_0_xz_xyzz_xz = pbuffer.data(idx_op_geom_020_gd + 230);

    auto tr_0_0_xz_xyzz_yy = pbuffer.data(idx_op_geom_020_gd + 231);

    auto tr_0_0_xz_xyzz_yz = pbuffer.data(idx_op_geom_020_gd + 232);

    auto tr_0_0_xz_xyzz_zz = pbuffer.data(idx_op_geom_020_gd + 233);

    #pragma omp simd aligned(tr_0_0_xz_xyzz_xx, tr_0_0_xz_xyzz_xy, tr_0_0_xz_xyzz_xz, tr_0_0_xz_xyzz_yy, tr_0_0_xz_xyzz_yz, tr_0_0_xz_xyzz_zz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_xxz, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_y, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_z, tr_xxyzz_zzz, tr_xxyzzz_xx, tr_xxyzzz_xy, tr_xxyzzz_xz, tr_xxyzzz_yy, tr_xxyzzz_yz, tr_xxyzzz_zz, tr_xyz_x, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xxxz, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_y, tr_xyzzz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xyzz_xx[i] = 2.0 * tr_yz_xx[i] - 2.0 * tr_yzz_xxz[i] * tke_0 - 2.0 * tr_yzzz_xx[i] * tbe_0 + 4.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xxx[i] * tke_0 - 4.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xy[i] = 2.0 * tr_yz_xy[i] - 2.0 * tr_yzz_xyz[i] * tke_0 - 2.0 * tr_yzzz_xy[i] * tbe_0 + 2.0 * tr_xyz_y[i] - 4.0 * tr_xyz_xxy[i] * tke_0 - 2.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_xz[i] = 2.0 * tr_yz_xz[i] + tr_yzz_x[i] - 2.0 * tr_yzz_xzz[i] * tke_0 - 2.0 * tr_yzzz_xz[i] * tbe_0 + 2.0 * tr_xyz_z[i] - 4.0 * tr_xyz_xxz[i] * tke_0 + tr_xyzz_0[i] - 2.0 * tr_xyzz_zz[i] * tke_0 - 2.0 * tr_xyzz_xx[i] * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 - 2.0 * tr_xxyzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_yy[i] = 2.0 * tr_yz_yy[i] - 2.0 * tr_yzz_yyz[i] * tke_0 - 2.0 * tr_yzzz_yy[i] * tbe_0 - 4.0 * tr_xyz_xyy[i] * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_yz[i] = 2.0 * tr_yz_yz[i] + tr_yzz_y[i] - 2.0 * tr_yzz_yzz[i] * tke_0 - 2.0 * tr_yzzz_yz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tke_0 - 2.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 - 2.0 * tr_xxyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xyzz_zz[i] = 2.0 * tr_yz_zz[i] + 2.0 * tr_yzz_z[i] - 2.0 * tr_yzz_zzz[i] * tke_0 - 2.0 * tr_yzzz_zz[i] * tbe_0 - 4.0 * tr_xyz_xzz[i] * tke_0 - 4.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zz[i] * tbe_0 - 4.0 * tr_xxyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 234-240 components of targeted buffer : GD

    auto tr_0_0_xz_xzzz_xx = pbuffer.data(idx_op_geom_020_gd + 234);

    auto tr_0_0_xz_xzzz_xy = pbuffer.data(idx_op_geom_020_gd + 235);

    auto tr_0_0_xz_xzzz_xz = pbuffer.data(idx_op_geom_020_gd + 236);

    auto tr_0_0_xz_xzzz_yy = pbuffer.data(idx_op_geom_020_gd + 237);

    auto tr_0_0_xz_xzzz_yz = pbuffer.data(idx_op_geom_020_gd + 238);

    auto tr_0_0_xz_xzzz_zz = pbuffer.data(idx_op_geom_020_gd + 239);

    #pragma omp simd aligned(tr_0_0_xz_xzzz_xx, tr_0_0_xz_xzzz_xy, tr_0_0_xz_xzzz_xz, tr_0_0_xz_xzzz_yy, tr_0_0_xz_xzzz_yz, tr_0_0_xz_xzzz_zz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_xxz, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_y, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_z, tr_xxzzz_zzz, tr_xxzzzz_xx, tr_xxzzzz_xy, tr_xxzzzz_xz, tr_xxzzzz_yy, tr_xxzzzz_yz, tr_xxzzzz_zz, tr_xzz_x, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xxxz, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_xxx, tr_xzzzz_xxy, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_y, tr_xzzzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xxz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_xzzz_xx[i] = 3.0 * tr_zz_xx[i] - 2.0 * tr_zzz_xxz[i] * tke_0 - 2.0 * tr_zzzz_xx[i] * tbe_0 + 6.0 * tr_xzz_x[i] - 6.0 * tr_xzz_xxx[i] * tke_0 - 4.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_xxx[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xy[i] = 3.0 * tr_zz_xy[i] - 2.0 * tr_zzz_xyz[i] * tke_0 - 2.0 * tr_zzzz_xy[i] * tbe_0 + 3.0 * tr_xzz_y[i] - 6.0 * tr_xzz_xxy[i] * tke_0 - 2.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_xz[i] = 3.0 * tr_zz_xz[i] + tr_zzz_x[i] - 2.0 * tr_zzz_xzz[i] * tke_0 - 2.0 * tr_zzzz_xz[i] * tbe_0 + 3.0 * tr_xzz_z[i] - 6.0 * tr_xzz_xxz[i] * tke_0 + tr_xzzz_0[i] - 2.0 * tr_xzzz_zz[i] * tke_0 - 2.0 * tr_xzzz_xx[i] * tke_0 + 4.0 * tr_xzzz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_z[i] * tbe_0 + 4.0 * tr_xzzzz_xxz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xz[i] * tbe_0 - 2.0 * tr_xxzzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_yy[i] = 3.0 * tr_zz_yy[i] - 2.0 * tr_zzz_yyz[i] * tke_0 - 2.0 * tr_zzzz_yy[i] * tbe_0 - 6.0 * tr_xzz_xyy[i] * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_yz[i] = 3.0 * tr_zz_yz[i] + tr_zzz_y[i] - 2.0 * tr_zzz_yzz[i] * tke_0 - 2.0 * tr_zzzz_yz[i] * tbe_0 - 6.0 * tr_xzz_xyz[i] * tke_0 - 2.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yz[i] * tbe_0 - 2.0 * tr_xxzzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_xzzz_zz[i] = 3.0 * tr_zz_zz[i] + 2.0 * tr_zzz_z[i] - 2.0 * tr_zzz_zzz[i] * tke_0 - 2.0 * tr_zzzz_zz[i] * tbe_0 - 6.0 * tr_xzz_xzz[i] * tke_0 - 4.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_zz[i] * tbe_0 - 4.0 * tr_xxzzz_z[i] * tbe_0 + 4.0 * tr_xxzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-246 components of targeted buffer : GD

    auto tr_0_0_xz_yyyy_xx = pbuffer.data(idx_op_geom_020_gd + 240);

    auto tr_0_0_xz_yyyy_xy = pbuffer.data(idx_op_geom_020_gd + 241);

    auto tr_0_0_xz_yyyy_xz = pbuffer.data(idx_op_geom_020_gd + 242);

    auto tr_0_0_xz_yyyy_yy = pbuffer.data(idx_op_geom_020_gd + 243);

    auto tr_0_0_xz_yyyy_yz = pbuffer.data(idx_op_geom_020_gd + 244);

    auto tr_0_0_xz_yyyy_zz = pbuffer.data(idx_op_geom_020_gd + 245);

    #pragma omp simd aligned(tr_0_0_xz_yyyy_xx, tr_0_0_xz_yyyy_xy, tr_0_0_xz_yyyy_xz, tr_0_0_xz_yyyy_yy, tr_0_0_xz_yyyy_yz, tr_0_0_xz_yyyy_zz, tr_xyyyy_x, tr_xyyyy_xxz, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_y, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_z, tr_xyyyy_zzz, tr_xyyyyz_xx, tr_xyyyyz_xy, tr_xyyyyz_xz, tr_xyyyyz_yy, tr_xyyyyz_yz, tr_xyyyyz_zz, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xxxz, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyz_x, tr_yyyyz_xxx, tr_yyyyz_xxy, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_y, tr_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyyy_xx[i] = -4.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_x[i] * tbe_0 + 4.0 * tr_yyyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xy[i] = -2.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_y[i] * tbe_0 + 4.0 * tr_yyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_xz[i] = tr_yyyy_0[i] - 2.0 * tr_yyyy_zz[i] * tke_0 - 2.0 * tr_yyyy_xx[i] * tke_0 + 4.0 * tr_yyyy_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_z[i] * tbe_0 + 4.0 * tr_yyyyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_yy[i] = 4.0 * tr_yyyy_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_yz[i] = -2.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_y[i] * tbe_0 + 4.0 * tr_xyyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyy_zz[i] = -4.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 246-252 components of targeted buffer : GD

    auto tr_0_0_xz_yyyz_xx = pbuffer.data(idx_op_geom_020_gd + 246);

    auto tr_0_0_xz_yyyz_xy = pbuffer.data(idx_op_geom_020_gd + 247);

    auto tr_0_0_xz_yyyz_xz = pbuffer.data(idx_op_geom_020_gd + 248);

    auto tr_0_0_xz_yyyz_yy = pbuffer.data(idx_op_geom_020_gd + 249);

    auto tr_0_0_xz_yyyz_yz = pbuffer.data(idx_op_geom_020_gd + 250);

    auto tr_0_0_xz_yyyz_zz = pbuffer.data(idx_op_geom_020_gd + 251);

    #pragma omp simd aligned(tr_0_0_xz_yyyz_xx, tr_0_0_xz_yyyz_xy, tr_0_0_xz_yyyz_xz, tr_0_0_xz_yyyz_yy, tr_0_0_xz_yyyz_yz, tr_0_0_xz_yyyz_zz, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_xxz, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_y, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_z, tr_xyyyz_zzz, tr_xyyyzz_xx, tr_xyyyzz_xy, tr_xyyyzz_xz, tr_xyyyzz_yy, tr_xyyyzz_yz, tr_xyyyzz_zz, tr_yyy_x, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xxxz, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_y, tr_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyyz_xx[i] = 2.0 * tr_yyy_x[i] - 2.0 * tr_yyy_xxx[i] * tke_0 - 4.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xy[i] = tr_yyy_y[i] - 2.0 * tr_yyy_xxy[i] * tke_0 - 2.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_y[i] * tbe_0 + 4.0 * tr_yyyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_xz[i] = tr_yyy_z[i] - 2.0 * tr_yyy_xxz[i] * tke_0 + tr_yyyz_0[i] - 2.0 * tr_yyyz_zz[i] * tke_0 - 2.0 * tr_yyyz_xx[i] * tke_0 + 4.0 * tr_yyyz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_z[i] * tbe_0 + 4.0 * tr_yyyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xz[i] * tbe_0 - 2.0 * tr_xyyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_yy[i] = -2.0 * tr_yyy_xyy[i] * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_yz[i] = -2.0 * tr_yyy_xyz[i] * tke_0 - 2.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yz[i] * tbe_0 - 2.0 * tr_xyyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyyz_zz[i] = -2.0 * tr_yyy_xzz[i] * tke_0 - 4.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_zz[i] * tbe_0 - 4.0 * tr_xyyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 252-258 components of targeted buffer : GD

    auto tr_0_0_xz_yyzz_xx = pbuffer.data(idx_op_geom_020_gd + 252);

    auto tr_0_0_xz_yyzz_xy = pbuffer.data(idx_op_geom_020_gd + 253);

    auto tr_0_0_xz_yyzz_xz = pbuffer.data(idx_op_geom_020_gd + 254);

    auto tr_0_0_xz_yyzz_yy = pbuffer.data(idx_op_geom_020_gd + 255);

    auto tr_0_0_xz_yyzz_yz = pbuffer.data(idx_op_geom_020_gd + 256);

    auto tr_0_0_xz_yyzz_zz = pbuffer.data(idx_op_geom_020_gd + 257);

    #pragma omp simd aligned(tr_0_0_xz_yyzz_xx, tr_0_0_xz_yyzz_xy, tr_0_0_xz_yyzz_xz, tr_0_0_xz_yyzz_yy, tr_0_0_xz_yyzz_yz, tr_0_0_xz_yyzz_zz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_xxz, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_y, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_z, tr_xyyzz_zzz, tr_xyyzzz_xx, tr_xyyzzz_xy, tr_xyyzzz_xz, tr_xyyzzz_yy, tr_xyyzzz_yz, tr_xyyzzz_zz, tr_yyz_x, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xxxz, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_y, tr_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yyzz_xx[i] = 4.0 * tr_yyz_x[i] - 4.0 * tr_yyz_xxx[i] * tke_0 - 4.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xy[i] = 2.0 * tr_yyz_y[i] - 4.0 * tr_yyz_xxy[i] * tke_0 - 2.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_y[i] * tbe_0 + 4.0 * tr_yyzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_xz[i] = 2.0 * tr_yyz_z[i] - 4.0 * tr_yyz_xxz[i] * tke_0 + tr_yyzz_0[i] - 2.0 * tr_yyzz_zz[i] * tke_0 - 2.0 * tr_yyzz_xx[i] * tke_0 + 4.0 * tr_yyzz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_z[i] * tbe_0 + 4.0 * tr_yyzzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 - 2.0 * tr_xyyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_yy[i] = -4.0 * tr_yyz_xyy[i] * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_yz[i] = -4.0 * tr_yyz_xyz[i] * tke_0 - 2.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 - 2.0 * tr_xyyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yyzz_zz[i] = -4.0 * tr_yyz_xzz[i] * tke_0 - 4.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_zz[i] * tbe_0 - 4.0 * tr_xyyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 258-264 components of targeted buffer : GD

    auto tr_0_0_xz_yzzz_xx = pbuffer.data(idx_op_geom_020_gd + 258);

    auto tr_0_0_xz_yzzz_xy = pbuffer.data(idx_op_geom_020_gd + 259);

    auto tr_0_0_xz_yzzz_xz = pbuffer.data(idx_op_geom_020_gd + 260);

    auto tr_0_0_xz_yzzz_yy = pbuffer.data(idx_op_geom_020_gd + 261);

    auto tr_0_0_xz_yzzz_yz = pbuffer.data(idx_op_geom_020_gd + 262);

    auto tr_0_0_xz_yzzz_zz = pbuffer.data(idx_op_geom_020_gd + 263);

    #pragma omp simd aligned(tr_0_0_xz_yzzz_xx, tr_0_0_xz_yzzz_xy, tr_0_0_xz_yzzz_xz, tr_0_0_xz_yzzz_yy, tr_0_0_xz_yzzz_yz, tr_0_0_xz_yzzz_zz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_xxz, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_y, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_z, tr_xyzzz_zzz, tr_xyzzzz_xx, tr_xyzzzz_xy, tr_xyzzzz_xz, tr_xyzzzz_yy, tr_xyzzzz_yz, tr_xyzzzz_zz, tr_yzz_x, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xxxz, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_xxx, tr_yzzzz_xxy, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_y, tr_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_yzzz_xx[i] = 6.0 * tr_yzz_x[i] - 6.0 * tr_yzz_xxx[i] * tke_0 - 4.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xy[i] = 3.0 * tr_yzz_y[i] - 6.0 * tr_yzz_xxy[i] * tke_0 - 2.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_xz[i] = 3.0 * tr_yzz_z[i] - 6.0 * tr_yzz_xxz[i] * tke_0 + tr_yzzz_0[i] - 2.0 * tr_yzzz_zz[i] * tke_0 - 2.0 * tr_yzzz_xx[i] * tke_0 + 4.0 * tr_yzzz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_z[i] * tbe_0 + 4.0 * tr_yzzzz_xxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xz[i] * tbe_0 - 2.0 * tr_xyzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_yy[i] = -6.0 * tr_yzz_xyy[i] * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_yz[i] = -6.0 * tr_yzz_xyz[i] * tke_0 - 2.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yz[i] * tbe_0 - 2.0 * tr_xyzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_yzzz_zz[i] = -6.0 * tr_yzz_xzz[i] * tke_0 - 4.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_zz[i] * tbe_0 - 4.0 * tr_xyzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 264-270 components of targeted buffer : GD

    auto tr_0_0_xz_zzzz_xx = pbuffer.data(idx_op_geom_020_gd + 264);

    auto tr_0_0_xz_zzzz_xy = pbuffer.data(idx_op_geom_020_gd + 265);

    auto tr_0_0_xz_zzzz_xz = pbuffer.data(idx_op_geom_020_gd + 266);

    auto tr_0_0_xz_zzzz_yy = pbuffer.data(idx_op_geom_020_gd + 267);

    auto tr_0_0_xz_zzzz_yz = pbuffer.data(idx_op_geom_020_gd + 268);

    auto tr_0_0_xz_zzzz_zz = pbuffer.data(idx_op_geom_020_gd + 269);

    #pragma omp simd aligned(tr_0_0_xz_zzzz_xx, tr_0_0_xz_zzzz_xy, tr_0_0_xz_zzzz_xz, tr_0_0_xz_zzzz_yy, tr_0_0_xz_zzzz_yz, tr_0_0_xz_zzzz_zz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_xxz, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_y, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_z, tr_xzzzz_zzz, tr_xzzzzz_xx, tr_xzzzzz_xy, tr_xzzzzz_xz, tr_xzzzzz_yy, tr_xzzzzz_yz, tr_xzzzzz_zz, tr_zzz_x, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xxxz, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yz, tr_zzzz_zz, tr_zzzzz_x, tr_zzzzz_xxx, tr_zzzzz_xxy, tr_zzzzz_xxz, tr_zzzzz_xyy, tr_zzzzz_xyz, tr_zzzzz_xzz, tr_zzzzz_y, tr_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_xz_zzzz_xx[i] = 8.0 * tr_zzz_x[i] - 8.0 * tr_zzz_xxx[i] * tke_0 - 4.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xxxz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_x[i] * tbe_0 + 4.0 * tr_zzzzz_xxx[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xy[i] = 4.0 * tr_zzz_y[i] - 8.0 * tr_zzz_xxy[i] * tke_0 - 2.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_xxyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_y[i] * tbe_0 + 4.0 * tr_zzzzz_xxy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_xz[i] = 4.0 * tr_zzz_z[i] - 8.0 * tr_zzz_xxz[i] * tke_0 + tr_zzzz_0[i] - 2.0 * tr_zzzz_zz[i] * tke_0 - 2.0 * tr_zzzz_xx[i] * tke_0 + 4.0 * tr_zzzz_xxzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_z[i] * tbe_0 + 4.0 * tr_zzzzz_xxz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xz[i] * tbe_0 - 2.0 * tr_xzzzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_yy[i] = -8.0 * tr_zzz_xyy[i] * tke_0 + 4.0 * tr_zzzz_xyyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_yz[i] = -8.0 * tr_zzz_xyz[i] * tke_0 - 2.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yz[i] * tbe_0 - 2.0 * tr_xzzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_xz_zzzz_zz[i] = -8.0 * tr_zzz_xzz[i] * tke_0 - 4.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_zz[i] * tbe_0 - 4.0 * tr_xzzzz_z[i] * tbe_0 + 4.0 * tr_xzzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-276 components of targeted buffer : GD

    auto tr_0_0_yy_xxxx_xx = pbuffer.data(idx_op_geom_020_gd + 270);

    auto tr_0_0_yy_xxxx_xy = pbuffer.data(idx_op_geom_020_gd + 271);

    auto tr_0_0_yy_xxxx_xz = pbuffer.data(idx_op_geom_020_gd + 272);

    auto tr_0_0_yy_xxxx_yy = pbuffer.data(idx_op_geom_020_gd + 273);

    auto tr_0_0_yy_xxxx_yz = pbuffer.data(idx_op_geom_020_gd + 274);

    auto tr_0_0_yy_xxxx_zz = pbuffer.data(idx_op_geom_020_gd + 275);

    #pragma omp simd aligned(tr_0_0_yy_xxxx_xx, tr_0_0_yy_xxxx_xy, tr_0_0_yy_xxxx_xz, tr_0_0_yy_xxxx_yy, tr_0_0_yy_xxxx_yz, tr_0_0_yy_xxxx_zz, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xxyy, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxy_x, tr_xxxxy_xxy, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_y, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_z, tr_xxxxyy_xx, tr_xxxxyy_xy, tr_xxxxyy_xz, tr_xxxxyy_yy, tr_xxxxyy_yz, tr_xxxxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxxx_xx[i] = -2.0 * tr_xxxx_xx[i] * tbe_0 - 2.0 * tr_xxxx_xx[i] * tke_0 + 4.0 * tr_xxxx_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xy[i] = -2.0 * tr_xxxx_xy[i] * tbe_0 - 6.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_x[i] * tbe_0 + 8.0 * tr_xxxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_xz[i] = -2.0 * tr_xxxx_xz[i] * tbe_0 - 2.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_yy[i] = 2.0 * tr_xxxx_0[i] - 2.0 * tr_xxxx_yy[i] * tbe_0 - 10.0 * tr_xxxx_yy[i] * tke_0 + 4.0 * tr_xxxx_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxxy_y[i] * tbe_0 + 8.0 * tr_xxxxy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_yz[i] = -2.0 * tr_xxxx_yz[i] * tbe_0 - 6.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxy_z[i] * tbe_0 + 8.0 * tr_xxxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxx_zz[i] = -2.0 * tr_xxxx_zz[i] * tbe_0 - 2.0 * tr_xxxx_zz[i] * tke_0 + 4.0 * tr_xxxx_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 276-282 components of targeted buffer : GD

    auto tr_0_0_yy_xxxy_xx = pbuffer.data(idx_op_geom_020_gd + 276);

    auto tr_0_0_yy_xxxy_xy = pbuffer.data(idx_op_geom_020_gd + 277);

    auto tr_0_0_yy_xxxy_xz = pbuffer.data(idx_op_geom_020_gd + 278);

    auto tr_0_0_yy_xxxy_yy = pbuffer.data(idx_op_geom_020_gd + 279);

    auto tr_0_0_yy_xxxy_yz = pbuffer.data(idx_op_geom_020_gd + 280);

    auto tr_0_0_yy_xxxy_zz = pbuffer.data(idx_op_geom_020_gd + 281);

    #pragma omp simd aligned(tr_0_0_yy_xxxy_xx, tr_0_0_yy_xxxy_xy, tr_0_0_yy_xxxy_xz, tr_0_0_yy_xxxy_yy, tr_0_0_yy_xxxy_yz, tr_0_0_yy_xxxy_zz, tr_xxx_x, tr_xxx_xxy, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xxyy, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyy_x, tr_xxxyy_xxy, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_y, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_z, tr_xxxyyy_xx, tr_xxxyyy_xy, tr_xxxyyy_xz, tr_xxxyyy_yy, tr_xxxyyy_yz, tr_xxxyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxxy_xx[i] = -4.0 * tr_xxx_xxy[i] * tke_0 - 6.0 * tr_xxxy_xx[i] * tbe_0 - 2.0 * tr_xxxy_xx[i] * tke_0 + 4.0 * tr_xxxy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xy[i] = 2.0 * tr_xxx_x[i] - 4.0 * tr_xxx_xyy[i] * tke_0 - 6.0 * tr_xxxy_xy[i] * tbe_0 - 6.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_x[i] * tbe_0 + 8.0 * tr_xxxyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_xz[i] = -4.0 * tr_xxx_xyz[i] * tke_0 - 6.0 * tr_xxxy_xz[i] * tbe_0 - 2.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_yy[i] = 4.0 * tr_xxx_y[i] - 4.0 * tr_xxx_yyy[i] * tke_0 + 2.0 * tr_xxxy_0[i] - 6.0 * tr_xxxy_yy[i] * tbe_0 - 10.0 * tr_xxxy_yy[i] * tke_0 + 4.0 * tr_xxxy_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyy_y[i] * tbe_0 + 8.0 * tr_xxxyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_yz[i] = 2.0 * tr_xxx_z[i] - 4.0 * tr_xxx_yyz[i] * tke_0 - 6.0 * tr_xxxy_yz[i] * tbe_0 - 6.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyy_z[i] * tbe_0 + 8.0 * tr_xxxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxy_zz[i] = -4.0 * tr_xxx_yzz[i] * tke_0 - 6.0 * tr_xxxy_zz[i] * tbe_0 - 2.0 * tr_xxxy_zz[i] * tke_0 + 4.0 * tr_xxxy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 282-288 components of targeted buffer : GD

    auto tr_0_0_yy_xxxz_xx = pbuffer.data(idx_op_geom_020_gd + 282);

    auto tr_0_0_yy_xxxz_xy = pbuffer.data(idx_op_geom_020_gd + 283);

    auto tr_0_0_yy_xxxz_xz = pbuffer.data(idx_op_geom_020_gd + 284);

    auto tr_0_0_yy_xxxz_yy = pbuffer.data(idx_op_geom_020_gd + 285);

    auto tr_0_0_yy_xxxz_yz = pbuffer.data(idx_op_geom_020_gd + 286);

    auto tr_0_0_yy_xxxz_zz = pbuffer.data(idx_op_geom_020_gd + 287);

    #pragma omp simd aligned(tr_0_0_yy_xxxz_xx, tr_0_0_yy_xxxz_xy, tr_0_0_yy_xxxz_xz, tr_0_0_yy_xxxz_yy, tr_0_0_yy_xxxz_yz, tr_0_0_yy_xxxz_zz, tr_xxxyyz_xx, tr_xxxyyz_xy, tr_xxxyyz_xz, tr_xxxyyz_yy, tr_xxxyyz_yz, tr_xxxyyz_zz, tr_xxxyz_x, tr_xxxyz_xxy, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_y, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_z, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xxyy, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxxz_xx[i] = -2.0 * tr_xxxz_xx[i] * tbe_0 - 2.0 * tr_xxxz_xx[i] * tke_0 + 4.0 * tr_xxxz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xy[i] = -2.0 * tr_xxxz_xy[i] * tbe_0 - 6.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_x[i] * tbe_0 + 8.0 * tr_xxxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_xz[i] = -2.0 * tr_xxxz_xz[i] * tbe_0 - 2.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_yy[i] = 2.0 * tr_xxxz_0[i] - 2.0 * tr_xxxz_yy[i] * tbe_0 - 10.0 * tr_xxxz_yy[i] * tke_0 + 4.0 * tr_xxxz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_y[i] * tbe_0 + 8.0 * tr_xxxyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_yz[i] = -2.0 * tr_xxxz_yz[i] * tbe_0 - 6.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_z[i] * tbe_0 + 8.0 * tr_xxxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxxz_zz[i] = -2.0 * tr_xxxz_zz[i] * tbe_0 - 2.0 * tr_xxxz_zz[i] * tke_0 + 4.0 * tr_xxxz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 288-294 components of targeted buffer : GD

    auto tr_0_0_yy_xxyy_xx = pbuffer.data(idx_op_geom_020_gd + 288);

    auto tr_0_0_yy_xxyy_xy = pbuffer.data(idx_op_geom_020_gd + 289);

    auto tr_0_0_yy_xxyy_xz = pbuffer.data(idx_op_geom_020_gd + 290);

    auto tr_0_0_yy_xxyy_yy = pbuffer.data(idx_op_geom_020_gd + 291);

    auto tr_0_0_yy_xxyy_yz = pbuffer.data(idx_op_geom_020_gd + 292);

    auto tr_0_0_yy_xxyy_zz = pbuffer.data(idx_op_geom_020_gd + 293);

    #pragma omp simd aligned(tr_0_0_yy_xxyy_xx, tr_0_0_yy_xxyy_xy, tr_0_0_yy_xxyy_xz, tr_0_0_yy_xxyy_yy, tr_0_0_yy_xxyy_yz, tr_0_0_yy_xxyy_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xxyy, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyy_x, tr_xxyyy_xxy, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_y, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_z, tr_xxyyyy_xx, tr_xxyyyy_xy, tr_xxyyyy_xz, tr_xxyyyy_yy, tr_xxyyyy_yz, tr_xxyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxyy_xx[i] = 2.0 * tr_xx_xx[i] - 8.0 * tr_xxy_xxy[i] * tke_0 - 10.0 * tr_xxyy_xx[i] * tbe_0 - 2.0 * tr_xxyy_xx[i] * tke_0 + 4.0 * tr_xxyy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xy[i] = 2.0 * tr_xx_xy[i] + 4.0 * tr_xxy_x[i] - 8.0 * tr_xxy_xyy[i] * tke_0 - 10.0 * tr_xxyy_xy[i] * tbe_0 - 6.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_x[i] * tbe_0 + 8.0 * tr_xxyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_xz[i] = 2.0 * tr_xx_xz[i] - 8.0 * tr_xxy_xyz[i] * tke_0 - 10.0 * tr_xxyy_xz[i] * tbe_0 - 2.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_yy[i] = 2.0 * tr_xx_yy[i] + 8.0 * tr_xxy_y[i] - 8.0 * tr_xxy_yyy[i] * tke_0 + 2.0 * tr_xxyy_0[i] - 10.0 * tr_xxyy_yy[i] * tbe_0 - 10.0 * tr_xxyy_yy[i] * tke_0 + 4.0 * tr_xxyy_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyy_y[i] * tbe_0 + 8.0 * tr_xxyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_yz[i] = 2.0 * tr_xx_yz[i] + 4.0 * tr_xxy_z[i] - 8.0 * tr_xxy_yyz[i] * tke_0 - 10.0 * tr_xxyy_yz[i] * tbe_0 - 6.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyy_z[i] * tbe_0 + 8.0 * tr_xxyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyy_zz[i] = 2.0 * tr_xx_zz[i] - 8.0 * tr_xxy_yzz[i] * tke_0 - 10.0 * tr_xxyy_zz[i] * tbe_0 - 2.0 * tr_xxyy_zz[i] * tke_0 + 4.0 * tr_xxyy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 294-300 components of targeted buffer : GD

    auto tr_0_0_yy_xxyz_xx = pbuffer.data(idx_op_geom_020_gd + 294);

    auto tr_0_0_yy_xxyz_xy = pbuffer.data(idx_op_geom_020_gd + 295);

    auto tr_0_0_yy_xxyz_xz = pbuffer.data(idx_op_geom_020_gd + 296);

    auto tr_0_0_yy_xxyz_yy = pbuffer.data(idx_op_geom_020_gd + 297);

    auto tr_0_0_yy_xxyz_yz = pbuffer.data(idx_op_geom_020_gd + 298);

    auto tr_0_0_yy_xxyz_zz = pbuffer.data(idx_op_geom_020_gd + 299);

    #pragma omp simd aligned(tr_0_0_yy_xxyz_xx, tr_0_0_yy_xxyz_xy, tr_0_0_yy_xxyz_xz, tr_0_0_yy_xxyz_yy, tr_0_0_yy_xxyz_yz, tr_0_0_yy_xxyz_zz, tr_xxyyyz_xx, tr_xxyyyz_xy, tr_xxyyyz_xz, tr_xxyyyz_yy, tr_xxyyyz_yz, tr_xxyyyz_zz, tr_xxyyz_x, tr_xxyyz_xxy, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_y, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_z, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xxyy, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_zz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxyz_xx[i] = -4.0 * tr_xxz_xxy[i] * tke_0 - 6.0 * tr_xxyz_xx[i] * tbe_0 - 2.0 * tr_xxyz_xx[i] * tke_0 + 4.0 * tr_xxyz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xy[i] = 2.0 * tr_xxz_x[i] - 4.0 * tr_xxz_xyy[i] * tke_0 - 6.0 * tr_xxyz_xy[i] * tbe_0 - 6.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_x[i] * tbe_0 + 8.0 * tr_xxyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_xz[i] = -4.0 * tr_xxz_xyz[i] * tke_0 - 6.0 * tr_xxyz_xz[i] * tbe_0 - 2.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_yy[i] = 4.0 * tr_xxz_y[i] - 4.0 * tr_xxz_yyy[i] * tke_0 + 2.0 * tr_xxyz_0[i] - 6.0 * tr_xxyz_yy[i] * tbe_0 - 10.0 * tr_xxyz_yy[i] * tke_0 + 4.0 * tr_xxyz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_y[i] * tbe_0 + 8.0 * tr_xxyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_yz[i] = 2.0 * tr_xxz_z[i] - 4.0 * tr_xxz_yyz[i] * tke_0 - 6.0 * tr_xxyz_yz[i] * tbe_0 - 6.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_z[i] * tbe_0 + 8.0 * tr_xxyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxyz_zz[i] = -4.0 * tr_xxz_yzz[i] * tke_0 - 6.0 * tr_xxyz_zz[i] * tbe_0 - 2.0 * tr_xxyz_zz[i] * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-306 components of targeted buffer : GD

    auto tr_0_0_yy_xxzz_xx = pbuffer.data(idx_op_geom_020_gd + 300);

    auto tr_0_0_yy_xxzz_xy = pbuffer.data(idx_op_geom_020_gd + 301);

    auto tr_0_0_yy_xxzz_xz = pbuffer.data(idx_op_geom_020_gd + 302);

    auto tr_0_0_yy_xxzz_yy = pbuffer.data(idx_op_geom_020_gd + 303);

    auto tr_0_0_yy_xxzz_yz = pbuffer.data(idx_op_geom_020_gd + 304);

    auto tr_0_0_yy_xxzz_zz = pbuffer.data(idx_op_geom_020_gd + 305);

    #pragma omp simd aligned(tr_0_0_yy_xxzz_xx, tr_0_0_yy_xxzz_xy, tr_0_0_yy_xxzz_xz, tr_0_0_yy_xxzz_yy, tr_0_0_yy_xxzz_yz, tr_0_0_yy_xxzz_zz, tr_xxyyzz_xx, tr_xxyyzz_xy, tr_xxyyzz_xz, tr_xxyyzz_yy, tr_xxyyzz_yz, tr_xxyyzz_zz, tr_xxyzz_x, tr_xxyzz_xxy, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_y, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_z, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xxyy, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xxzz_xx[i] = -2.0 * tr_xxzz_xx[i] * tbe_0 - 2.0 * tr_xxzz_xx[i] * tke_0 + 4.0 * tr_xxzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xy[i] = -2.0 * tr_xxzz_xy[i] * tbe_0 - 6.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_x[i] * tbe_0 + 8.0 * tr_xxyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_xz[i] = -2.0 * tr_xxzz_xz[i] * tbe_0 - 2.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_yy[i] = 2.0 * tr_xxzz_0[i] - 2.0 * tr_xxzz_yy[i] * tbe_0 - 10.0 * tr_xxzz_yy[i] * tke_0 + 4.0 * tr_xxzz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_y[i] * tbe_0 + 8.0 * tr_xxyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_yz[i] = -2.0 * tr_xxzz_yz[i] * tbe_0 - 6.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_z[i] * tbe_0 + 8.0 * tr_xxyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xxzz_zz[i] = -2.0 * tr_xxzz_zz[i] * tbe_0 - 2.0 * tr_xxzz_zz[i] * tke_0 + 4.0 * tr_xxzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 306-312 components of targeted buffer : GD

    auto tr_0_0_yy_xyyy_xx = pbuffer.data(idx_op_geom_020_gd + 306);

    auto tr_0_0_yy_xyyy_xy = pbuffer.data(idx_op_geom_020_gd + 307);

    auto tr_0_0_yy_xyyy_xz = pbuffer.data(idx_op_geom_020_gd + 308);

    auto tr_0_0_yy_xyyy_yy = pbuffer.data(idx_op_geom_020_gd + 309);

    auto tr_0_0_yy_xyyy_yz = pbuffer.data(idx_op_geom_020_gd + 310);

    auto tr_0_0_yy_xyyy_zz = pbuffer.data(idx_op_geom_020_gd + 311);

    #pragma omp simd aligned(tr_0_0_yy_xyyy_xx, tr_0_0_yy_xyyy_xy, tr_0_0_yy_xyyy_xz, tr_0_0_yy_xyyy_yy, tr_0_0_yy_xyyy_yz, tr_0_0_yy_xyyy_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xxyy, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyy_x, tr_xyyyy_xxy, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_y, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_z, tr_xyyyyy_xx, tr_xyyyyy_xy, tr_xyyyyy_xz, tr_xyyyyy_yy, tr_xyyyyy_yz, tr_xyyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyyy_xx[i] = 6.0 * tr_xy_xx[i] - 12.0 * tr_xyy_xxy[i] * tke_0 - 14.0 * tr_xyyy_xx[i] * tbe_0 - 2.0 * tr_xyyy_xx[i] * tke_0 + 4.0 * tr_xyyy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xy[i] = 6.0 * tr_xy_xy[i] + 6.0 * tr_xyy_x[i] - 12.0 * tr_xyy_xyy[i] * tke_0 - 14.0 * tr_xyyy_xy[i] * tbe_0 - 6.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_x[i] * tbe_0 + 8.0 * tr_xyyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_xz[i] = 6.0 * tr_xy_xz[i] - 12.0 * tr_xyy_xyz[i] * tke_0 - 14.0 * tr_xyyy_xz[i] * tbe_0 - 2.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_yy[i] = 6.0 * tr_xy_yy[i] + 12.0 * tr_xyy_y[i] - 12.0 * tr_xyy_yyy[i] * tke_0 + 2.0 * tr_xyyy_0[i] - 14.0 * tr_xyyy_yy[i] * tbe_0 - 10.0 * tr_xyyy_yy[i] * tke_0 + 4.0 * tr_xyyy_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyy_y[i] * tbe_0 + 8.0 * tr_xyyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_yz[i] = 6.0 * tr_xy_yz[i] + 6.0 * tr_xyy_z[i] - 12.0 * tr_xyy_yyz[i] * tke_0 - 14.0 * tr_xyyy_yz[i] * tbe_0 - 6.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyy_z[i] * tbe_0 + 8.0 * tr_xyyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyy_zz[i] = 6.0 * tr_xy_zz[i] - 12.0 * tr_xyy_yzz[i] * tke_0 - 14.0 * tr_xyyy_zz[i] * tbe_0 - 2.0 * tr_xyyy_zz[i] * tke_0 + 4.0 * tr_xyyy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 312-318 components of targeted buffer : GD

    auto tr_0_0_yy_xyyz_xx = pbuffer.data(idx_op_geom_020_gd + 312);

    auto tr_0_0_yy_xyyz_xy = pbuffer.data(idx_op_geom_020_gd + 313);

    auto tr_0_0_yy_xyyz_xz = pbuffer.data(idx_op_geom_020_gd + 314);

    auto tr_0_0_yy_xyyz_yy = pbuffer.data(idx_op_geom_020_gd + 315);

    auto tr_0_0_yy_xyyz_yz = pbuffer.data(idx_op_geom_020_gd + 316);

    auto tr_0_0_yy_xyyz_zz = pbuffer.data(idx_op_geom_020_gd + 317);

    #pragma omp simd aligned(tr_0_0_yy_xyyz_xx, tr_0_0_yy_xyyz_xy, tr_0_0_yy_xyyz_xz, tr_0_0_yy_xyyz_yy, tr_0_0_yy_xyyz_yz, tr_0_0_yy_xyyz_zz, tr_xyyyyz_xx, tr_xyyyyz_xy, tr_xyyyyz_xz, tr_xyyyyz_yy, tr_xyyyyz_yz, tr_xyyyyz_zz, tr_xyyyz_x, tr_xyyyz_xxy, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_y, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_z, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xxyy, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyyz_xx[i] = 2.0 * tr_xz_xx[i] - 8.0 * tr_xyz_xxy[i] * tke_0 - 10.0 * tr_xyyz_xx[i] * tbe_0 - 2.0 * tr_xyyz_xx[i] * tke_0 + 4.0 * tr_xyyz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xy[i] = 2.0 * tr_xz_xy[i] + 4.0 * tr_xyz_x[i] - 8.0 * tr_xyz_xyy[i] * tke_0 - 10.0 * tr_xyyz_xy[i] * tbe_0 - 6.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_x[i] * tbe_0 + 8.0 * tr_xyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_xz[i] = 2.0 * tr_xz_xz[i] - 8.0 * tr_xyz_xyz[i] * tke_0 - 10.0 * tr_xyyz_xz[i] * tbe_0 - 2.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_yy[i] = 2.0 * tr_xz_yy[i] + 8.0 * tr_xyz_y[i] - 8.0 * tr_xyz_yyy[i] * tke_0 + 2.0 * tr_xyyz_0[i] - 10.0 * tr_xyyz_yy[i] * tbe_0 - 10.0 * tr_xyyz_yy[i] * tke_0 + 4.0 * tr_xyyz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_y[i] * tbe_0 + 8.0 * tr_xyyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_yz[i] = 2.0 * tr_xz_yz[i] + 4.0 * tr_xyz_z[i] - 8.0 * tr_xyz_yyz[i] * tke_0 - 10.0 * tr_xyyz_yz[i] * tbe_0 - 6.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_z[i] * tbe_0 + 8.0 * tr_xyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyyz_zz[i] = 2.0 * tr_xz_zz[i] - 8.0 * tr_xyz_yzz[i] * tke_0 - 10.0 * tr_xyyz_zz[i] * tbe_0 - 2.0 * tr_xyyz_zz[i] * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 318-324 components of targeted buffer : GD

    auto tr_0_0_yy_xyzz_xx = pbuffer.data(idx_op_geom_020_gd + 318);

    auto tr_0_0_yy_xyzz_xy = pbuffer.data(idx_op_geom_020_gd + 319);

    auto tr_0_0_yy_xyzz_xz = pbuffer.data(idx_op_geom_020_gd + 320);

    auto tr_0_0_yy_xyzz_yy = pbuffer.data(idx_op_geom_020_gd + 321);

    auto tr_0_0_yy_xyzz_yz = pbuffer.data(idx_op_geom_020_gd + 322);

    auto tr_0_0_yy_xyzz_zz = pbuffer.data(idx_op_geom_020_gd + 323);

    #pragma omp simd aligned(tr_0_0_yy_xyzz_xx, tr_0_0_yy_xyzz_xy, tr_0_0_yy_xyzz_xz, tr_0_0_yy_xyzz_yy, tr_0_0_yy_xyzz_yz, tr_0_0_yy_xyzz_zz, tr_xyyyzz_xx, tr_xyyyzz_xy, tr_xyyyzz_xz, tr_xyyyzz_yy, tr_xyyyzz_yz, tr_xyyyzz_zz, tr_xyyzz_x, tr_xyyzz_xxy, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_y, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_z, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xxyy, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_zz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xyzz_xx[i] = -4.0 * tr_xzz_xxy[i] * tke_0 - 6.0 * tr_xyzz_xx[i] * tbe_0 - 2.0 * tr_xyzz_xx[i] * tke_0 + 4.0 * tr_xyzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xy[i] = 2.0 * tr_xzz_x[i] - 4.0 * tr_xzz_xyy[i] * tke_0 - 6.0 * tr_xyzz_xy[i] * tbe_0 - 6.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_x[i] * tbe_0 + 8.0 * tr_xyyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_xz[i] = -4.0 * tr_xzz_xyz[i] * tke_0 - 6.0 * tr_xyzz_xz[i] * tbe_0 - 2.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_yy[i] = 4.0 * tr_xzz_y[i] - 4.0 * tr_xzz_yyy[i] * tke_0 + 2.0 * tr_xyzz_0[i] - 6.0 * tr_xyzz_yy[i] * tbe_0 - 10.0 * tr_xyzz_yy[i] * tke_0 + 4.0 * tr_xyzz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_y[i] * tbe_0 + 8.0 * tr_xyyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_yz[i] = 2.0 * tr_xzz_z[i] - 4.0 * tr_xzz_yyz[i] * tke_0 - 6.0 * tr_xyzz_yz[i] * tbe_0 - 6.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_z[i] * tbe_0 + 8.0 * tr_xyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xyzz_zz[i] = -4.0 * tr_xzz_yzz[i] * tke_0 - 6.0 * tr_xyzz_zz[i] * tbe_0 - 2.0 * tr_xyzz_zz[i] * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 324-330 components of targeted buffer : GD

    auto tr_0_0_yy_xzzz_xx = pbuffer.data(idx_op_geom_020_gd + 324);

    auto tr_0_0_yy_xzzz_xy = pbuffer.data(idx_op_geom_020_gd + 325);

    auto tr_0_0_yy_xzzz_xz = pbuffer.data(idx_op_geom_020_gd + 326);

    auto tr_0_0_yy_xzzz_yy = pbuffer.data(idx_op_geom_020_gd + 327);

    auto tr_0_0_yy_xzzz_yz = pbuffer.data(idx_op_geom_020_gd + 328);

    auto tr_0_0_yy_xzzz_zz = pbuffer.data(idx_op_geom_020_gd + 329);

    #pragma omp simd aligned(tr_0_0_yy_xzzz_xx, tr_0_0_yy_xzzz_xy, tr_0_0_yy_xzzz_xz, tr_0_0_yy_xzzz_yy, tr_0_0_yy_xzzz_yz, tr_0_0_yy_xzzz_zz, tr_xyyzzz_xx, tr_xyyzzz_xy, tr_xyyzzz_xz, tr_xyyzzz_yy, tr_xyyzzz_yz, tr_xyyzzz_zz, tr_xyzzz_x, tr_xyzzz_xxy, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_y, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_z, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xxyy, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_xzzz_xx[i] = -2.0 * tr_xzzz_xx[i] * tbe_0 - 2.0 * tr_xzzz_xx[i] * tke_0 + 4.0 * tr_xzzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xy[i] = -2.0 * tr_xzzz_xy[i] * tbe_0 - 6.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_x[i] * tbe_0 + 8.0 * tr_xyzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_xz[i] = -2.0 * tr_xzzz_xz[i] * tbe_0 - 2.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_yy[i] = 2.0 * tr_xzzz_0[i] - 2.0 * tr_xzzz_yy[i] * tbe_0 - 10.0 * tr_xzzz_yy[i] * tke_0 + 4.0 * tr_xzzz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_y[i] * tbe_0 + 8.0 * tr_xyzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_yz[i] = -2.0 * tr_xzzz_yz[i] * tbe_0 - 6.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_z[i] * tbe_0 + 8.0 * tr_xyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_xzzz_zz[i] = -2.0 * tr_xzzz_zz[i] * tbe_0 - 2.0 * tr_xzzz_zz[i] * tke_0 + 4.0 * tr_xzzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-336 components of targeted buffer : GD

    auto tr_0_0_yy_yyyy_xx = pbuffer.data(idx_op_geom_020_gd + 330);

    auto tr_0_0_yy_yyyy_xy = pbuffer.data(idx_op_geom_020_gd + 331);

    auto tr_0_0_yy_yyyy_xz = pbuffer.data(idx_op_geom_020_gd + 332);

    auto tr_0_0_yy_yyyy_yy = pbuffer.data(idx_op_geom_020_gd + 333);

    auto tr_0_0_yy_yyyy_yz = pbuffer.data(idx_op_geom_020_gd + 334);

    auto tr_0_0_yy_yyyy_zz = pbuffer.data(idx_op_geom_020_gd + 335);

    #pragma omp simd aligned(tr_0_0_yy_yyyy_xx, tr_0_0_yy_yyyy_xy, tr_0_0_yy_yyyy_xz, tr_0_0_yy_yyyy_yy, tr_0_0_yy_yyyy_yz, tr_0_0_yy_yyyy_zz, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyy_x, tr_yyy_xxy, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xxyy, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyy_x, tr_yyyyy_xxy, tr_yyyyy_xyy, tr_yyyyy_xyz, tr_yyyyy_y, tr_yyyyy_yyy, tr_yyyyy_yyz, tr_yyyyy_yzz, tr_yyyyy_z, tr_yyyyyy_xx, tr_yyyyyy_xy, tr_yyyyyy_xz, tr_yyyyyy_yy, tr_yyyyyy_yz, tr_yyyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyyy_xx[i] = 12.0 * tr_yy_xx[i] - 16.0 * tr_yyy_xxy[i] * tke_0 - 18.0 * tr_yyyy_xx[i] * tbe_0 - 2.0 * tr_yyyy_xx[i] * tke_0 + 4.0 * tr_yyyy_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xy[i] = 12.0 * tr_yy_xy[i] + 8.0 * tr_yyy_x[i] - 16.0 * tr_yyy_xyy[i] * tke_0 - 18.0 * tr_yyyy_xy[i] * tbe_0 - 6.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_x[i] * tbe_0 + 8.0 * tr_yyyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_xz[i] = 12.0 * tr_yy_xz[i] - 16.0 * tr_yyy_xyz[i] * tke_0 - 18.0 * tr_yyyy_xz[i] * tbe_0 - 2.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_yy[i] = 12.0 * tr_yy_yy[i] + 16.0 * tr_yyy_y[i] - 16.0 * tr_yyy_yyy[i] * tke_0 + 2.0 * tr_yyyy_0[i] - 18.0 * tr_yyyy_yy[i] * tbe_0 - 10.0 * tr_yyyy_yy[i] * tke_0 + 4.0 * tr_yyyy_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyyy_y[i] * tbe_0 + 8.0 * tr_yyyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_yz[i] = 12.0 * tr_yy_yz[i] + 8.0 * tr_yyy_z[i] - 16.0 * tr_yyy_yyz[i] * tke_0 - 18.0 * tr_yyyy_yz[i] * tbe_0 - 6.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyy_z[i] * tbe_0 + 8.0 * tr_yyyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyy_zz[i] = 12.0 * tr_yy_zz[i] - 16.0 * tr_yyy_yzz[i] * tke_0 - 18.0 * tr_yyyy_zz[i] * tbe_0 - 2.0 * tr_yyyy_zz[i] * tke_0 + 4.0 * tr_yyyy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 336-342 components of targeted buffer : GD

    auto tr_0_0_yy_yyyz_xx = pbuffer.data(idx_op_geom_020_gd + 336);

    auto tr_0_0_yy_yyyz_xy = pbuffer.data(idx_op_geom_020_gd + 337);

    auto tr_0_0_yy_yyyz_xz = pbuffer.data(idx_op_geom_020_gd + 338);

    auto tr_0_0_yy_yyyz_yy = pbuffer.data(idx_op_geom_020_gd + 339);

    auto tr_0_0_yy_yyyz_yz = pbuffer.data(idx_op_geom_020_gd + 340);

    auto tr_0_0_yy_yyyz_zz = pbuffer.data(idx_op_geom_020_gd + 341);

    #pragma omp simd aligned(tr_0_0_yy_yyyz_xx, tr_0_0_yy_yyyz_xy, tr_0_0_yy_yyyz_xz, tr_0_0_yy_yyyz_yy, tr_0_0_yy_yyyz_yz, tr_0_0_yy_yyyz_zz, tr_yyyyyz_xx, tr_yyyyyz_xy, tr_yyyyyz_xz, tr_yyyyyz_yy, tr_yyyyyz_yz, tr_yyyyyz_zz, tr_yyyyz_x, tr_yyyyz_xxy, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_y, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_z, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xxyy, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_zz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyyz_xx[i] = 6.0 * tr_yz_xx[i] - 12.0 * tr_yyz_xxy[i] * tke_0 - 14.0 * tr_yyyz_xx[i] * tbe_0 - 2.0 * tr_yyyz_xx[i] * tke_0 + 4.0 * tr_yyyz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xy[i] = 6.0 * tr_yz_xy[i] + 6.0 * tr_yyz_x[i] - 12.0 * tr_yyz_xyy[i] * tke_0 - 14.0 * tr_yyyz_xy[i] * tbe_0 - 6.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_x[i] * tbe_0 + 8.0 * tr_yyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_xz[i] = 6.0 * tr_yz_xz[i] - 12.0 * tr_yyz_xyz[i] * tke_0 - 14.0 * tr_yyyz_xz[i] * tbe_0 - 2.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_yy[i] = 6.0 * tr_yz_yy[i] + 12.0 * tr_yyz_y[i] - 12.0 * tr_yyz_yyy[i] * tke_0 + 2.0 * tr_yyyz_0[i] - 14.0 * tr_yyyz_yy[i] * tbe_0 - 10.0 * tr_yyyz_yy[i] * tke_0 + 4.0 * tr_yyyz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_y[i] * tbe_0 + 8.0 * tr_yyyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_yz[i] = 6.0 * tr_yz_yz[i] + 6.0 * tr_yyz_z[i] - 12.0 * tr_yyz_yyz[i] * tke_0 - 14.0 * tr_yyyz_yz[i] * tbe_0 - 6.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_z[i] * tbe_0 + 8.0 * tr_yyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyyz_zz[i] = 6.0 * tr_yz_zz[i] - 12.0 * tr_yyz_yzz[i] * tke_0 - 14.0 * tr_yyyz_zz[i] * tbe_0 - 2.0 * tr_yyyz_zz[i] * tke_0 + 4.0 * tr_yyyz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 342-348 components of targeted buffer : GD

    auto tr_0_0_yy_yyzz_xx = pbuffer.data(idx_op_geom_020_gd + 342);

    auto tr_0_0_yy_yyzz_xy = pbuffer.data(idx_op_geom_020_gd + 343);

    auto tr_0_0_yy_yyzz_xz = pbuffer.data(idx_op_geom_020_gd + 344);

    auto tr_0_0_yy_yyzz_yy = pbuffer.data(idx_op_geom_020_gd + 345);

    auto tr_0_0_yy_yyzz_yz = pbuffer.data(idx_op_geom_020_gd + 346);

    auto tr_0_0_yy_yyzz_zz = pbuffer.data(idx_op_geom_020_gd + 347);

    #pragma omp simd aligned(tr_0_0_yy_yyzz_xx, tr_0_0_yy_yyzz_xy, tr_0_0_yy_yyzz_xz, tr_0_0_yy_yyzz_yy, tr_0_0_yy_yyzz_yz, tr_0_0_yy_yyzz_zz, tr_yyyyzz_xx, tr_yyyyzz_xy, tr_yyyyzz_xz, tr_yyyyzz_yy, tr_yyyyzz_yz, tr_yyyyzz_zz, tr_yyyzz_x, tr_yyyzz_xxy, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_y, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_z, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xxyy, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_zz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yyzz_xx[i] = 2.0 * tr_zz_xx[i] - 8.0 * tr_yzz_xxy[i] * tke_0 - 10.0 * tr_yyzz_xx[i] * tbe_0 - 2.0 * tr_yyzz_xx[i] * tke_0 + 4.0 * tr_yyzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xy[i] = 2.0 * tr_zz_xy[i] + 4.0 * tr_yzz_x[i] - 8.0 * tr_yzz_xyy[i] * tke_0 - 10.0 * tr_yyzz_xy[i] * tbe_0 - 6.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_x[i] * tbe_0 + 8.0 * tr_yyyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_xz[i] = 2.0 * tr_zz_xz[i] - 8.0 * tr_yzz_xyz[i] * tke_0 - 10.0 * tr_yyzz_xz[i] * tbe_0 - 2.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_yy[i] = 2.0 * tr_zz_yy[i] + 8.0 * tr_yzz_y[i] - 8.0 * tr_yzz_yyy[i] * tke_0 + 2.0 * tr_yyzz_0[i] - 10.0 * tr_yyzz_yy[i] * tbe_0 - 10.0 * tr_yyzz_yy[i] * tke_0 + 4.0 * tr_yyzz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_y[i] * tbe_0 + 8.0 * tr_yyyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_yz[i] = 2.0 * tr_zz_yz[i] + 4.0 * tr_yzz_z[i] - 8.0 * tr_yzz_yyz[i] * tke_0 - 10.0 * tr_yyzz_yz[i] * tbe_0 - 6.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_z[i] * tbe_0 + 8.0 * tr_yyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yyzz_zz[i] = 2.0 * tr_zz_zz[i] - 8.0 * tr_yzz_yzz[i] * tke_0 - 10.0 * tr_yyzz_zz[i] * tbe_0 - 2.0 * tr_yyzz_zz[i] * tke_0 + 4.0 * tr_yyzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 348-354 components of targeted buffer : GD

    auto tr_0_0_yy_yzzz_xx = pbuffer.data(idx_op_geom_020_gd + 348);

    auto tr_0_0_yy_yzzz_xy = pbuffer.data(idx_op_geom_020_gd + 349);

    auto tr_0_0_yy_yzzz_xz = pbuffer.data(idx_op_geom_020_gd + 350);

    auto tr_0_0_yy_yzzz_yy = pbuffer.data(idx_op_geom_020_gd + 351);

    auto tr_0_0_yy_yzzz_yz = pbuffer.data(idx_op_geom_020_gd + 352);

    auto tr_0_0_yy_yzzz_zz = pbuffer.data(idx_op_geom_020_gd + 353);

    #pragma omp simd aligned(tr_0_0_yy_yzzz_xx, tr_0_0_yy_yzzz_xy, tr_0_0_yy_yzzz_xz, tr_0_0_yy_yzzz_yy, tr_0_0_yy_yzzz_yz, tr_0_0_yy_yzzz_zz, tr_yyyzzz_xx, tr_yyyzzz_xy, tr_yyyzzz_xz, tr_yyyzzz_yy, tr_yyyzzz_yz, tr_yyyzzz_zz, tr_yyzzz_x, tr_yyzzz_xxy, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_y, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_z, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xxyy, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_zz, tr_zzz_x, tr_zzz_xxy, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_yzzz_xx[i] = -4.0 * tr_zzz_xxy[i] * tke_0 - 6.0 * tr_yzzz_xx[i] * tbe_0 - 2.0 * tr_yzzz_xx[i] * tke_0 + 4.0 * tr_yzzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xy[i] = 2.0 * tr_zzz_x[i] - 4.0 * tr_zzz_xyy[i] * tke_0 - 6.0 * tr_yzzz_xy[i] * tbe_0 - 6.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_x[i] * tbe_0 + 8.0 * tr_yyzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_xz[i] = -4.0 * tr_zzz_xyz[i] * tke_0 - 6.0 * tr_yzzz_xz[i] * tbe_0 - 2.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_yy[i] = 4.0 * tr_zzz_y[i] - 4.0 * tr_zzz_yyy[i] * tke_0 + 2.0 * tr_yzzz_0[i] - 6.0 * tr_yzzz_yy[i] * tbe_0 - 10.0 * tr_yzzz_yy[i] * tke_0 + 4.0 * tr_yzzz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_y[i] * tbe_0 + 8.0 * tr_yyzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_yz[i] = 2.0 * tr_zzz_z[i] - 4.0 * tr_zzz_yyz[i] * tke_0 - 6.0 * tr_yzzz_yz[i] * tbe_0 - 6.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_z[i] * tbe_0 + 8.0 * tr_yyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_yzzz_zz[i] = -4.0 * tr_zzz_yzz[i] * tke_0 - 6.0 * tr_yzzz_zz[i] * tbe_0 - 2.0 * tr_yzzz_zz[i] * tke_0 + 4.0 * tr_yzzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 354-360 components of targeted buffer : GD

    auto tr_0_0_yy_zzzz_xx = pbuffer.data(idx_op_geom_020_gd + 354);

    auto tr_0_0_yy_zzzz_xy = pbuffer.data(idx_op_geom_020_gd + 355);

    auto tr_0_0_yy_zzzz_xz = pbuffer.data(idx_op_geom_020_gd + 356);

    auto tr_0_0_yy_zzzz_yy = pbuffer.data(idx_op_geom_020_gd + 357);

    auto tr_0_0_yy_zzzz_yz = pbuffer.data(idx_op_geom_020_gd + 358);

    auto tr_0_0_yy_zzzz_zz = pbuffer.data(idx_op_geom_020_gd + 359);

    #pragma omp simd aligned(tr_0_0_yy_zzzz_xx, tr_0_0_yy_zzzz_xy, tr_0_0_yy_zzzz_xz, tr_0_0_yy_zzzz_yy, tr_0_0_yy_zzzz_yz, tr_0_0_yy_zzzz_zz, tr_yyzzzz_xx, tr_yyzzzz_xy, tr_yyzzzz_xz, tr_yyzzzz_yy, tr_yyzzzz_yz, tr_yyzzzz_zz, tr_yzzzz_x, tr_yzzzz_xxy, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_y, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_z, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xxyy, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yy_zzzz_xx[i] = -2.0 * tr_zzzz_xx[i] * tbe_0 - 2.0 * tr_zzzz_xx[i] * tke_0 + 4.0 * tr_zzzz_xxyy[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xy[i] = -2.0 * tr_zzzz_xy[i] * tbe_0 - 6.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xyyy[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_x[i] * tbe_0 + 8.0 * tr_yzzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_xz[i] = -2.0 * tr_zzzz_xz[i] * tbe_0 - 2.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xyyz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_yy[i] = 2.0 * tr_zzzz_0[i] - 2.0 * tr_zzzz_yy[i] * tbe_0 - 10.0 * tr_zzzz_yy[i] * tke_0 + 4.0 * tr_zzzz_yyyy[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_y[i] * tbe_0 + 8.0 * tr_yzzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_yz[i] = -2.0 * tr_zzzz_yz[i] * tbe_0 - 6.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_z[i] * tbe_0 + 8.0 * tr_yzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yy_zzzz_zz[i] = -2.0 * tr_zzzz_zz[i] * tbe_0 - 2.0 * tr_zzzz_zz[i] * tke_0 + 4.0 * tr_zzzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-366 components of targeted buffer : GD

    auto tr_0_0_yz_xxxx_xx = pbuffer.data(idx_op_geom_020_gd + 360);

    auto tr_0_0_yz_xxxx_xy = pbuffer.data(idx_op_geom_020_gd + 361);

    auto tr_0_0_yz_xxxx_xz = pbuffer.data(idx_op_geom_020_gd + 362);

    auto tr_0_0_yz_xxxx_yy = pbuffer.data(idx_op_geom_020_gd + 363);

    auto tr_0_0_yz_xxxx_yz = pbuffer.data(idx_op_geom_020_gd + 364);

    auto tr_0_0_yz_xxxx_zz = pbuffer.data(idx_op_geom_020_gd + 365);

    #pragma omp simd aligned(tr_0_0_yz_xxxx_xx, tr_0_0_yz_xxxx_xy, tr_0_0_yz_xxxx_xz, tr_0_0_yz_xxxx_yy, tr_0_0_yz_xxxx_yz, tr_0_0_yz_xxxx_zz, tr_xxxx_0, tr_xxxx_xxyz, tr_xxxx_xy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxxy_x, tr_xxxxy_xxz, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_y, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_z, tr_xxxxy_zzz, tr_xxxxyz_xx, tr_xxxxyz_xy, tr_xxxxyz_xz, tr_xxxxyz_yy, tr_xxxxyz_yz, tr_xxxxyz_zz, tr_xxxxz_x, tr_xxxxz_xxy, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_y, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxxx_xx[i] = 4.0 * tr_xxxx_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xy[i] = -2.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_x[i] * tbe_0 + 4.0 * tr_xxxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_xz[i] = -2.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_x[i] * tbe_0 + 4.0 * tr_xxxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_yy[i] = -4.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_y[i] * tbe_0 + 4.0 * tr_xxxxz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_yz[i] = tr_xxxx_0[i] - 2.0 * tr_xxxx_zz[i] * tke_0 - 2.0 * tr_xxxx_yy[i] * tke_0 + 4.0 * tr_xxxx_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxxz_z[i] * tbe_0 + 4.0 * tr_xxxxz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxxy_y[i] * tbe_0 + 4.0 * tr_xxxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxx_zz[i] = -4.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxxz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxxy_z[i] * tbe_0 + 4.0 * tr_xxxxy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 366-372 components of targeted buffer : GD

    auto tr_0_0_yz_xxxy_xx = pbuffer.data(idx_op_geom_020_gd + 366);

    auto tr_0_0_yz_xxxy_xy = pbuffer.data(idx_op_geom_020_gd + 367);

    auto tr_0_0_yz_xxxy_xz = pbuffer.data(idx_op_geom_020_gd + 368);

    auto tr_0_0_yz_xxxy_yy = pbuffer.data(idx_op_geom_020_gd + 369);

    auto tr_0_0_yz_xxxy_yz = pbuffer.data(idx_op_geom_020_gd + 370);

    auto tr_0_0_yz_xxxy_zz = pbuffer.data(idx_op_geom_020_gd + 371);

    #pragma omp simd aligned(tr_0_0_yz_xxxy_xx, tr_0_0_yz_xxxy_xy, tr_0_0_yz_xxxy_xz, tr_0_0_yz_xxxy_yy, tr_0_0_yz_xxxy_yz, tr_0_0_yz_xxxy_zz, tr_xxx_x, tr_xxx_xxz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxy_0, tr_xxxy_xxyz, tr_xxxy_xy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxyy_x, tr_xxxyy_xxz, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_y, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_z, tr_xxxyy_zzz, tr_xxxyyz_xx, tr_xxxyyz_xy, tr_xxxyyz_xz, tr_xxxyyz_yy, tr_xxxyyz_yz, tr_xxxyyz_zz, tr_xxxyz_x, tr_xxxyz_xxy, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_y, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_z, tr_xxxz_xx, tr_xxxz_xy, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxxy_xx[i] = -2.0 * tr_xxx_xxz[i] * tke_0 - 2.0 * tr_xxxz_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xy[i] = -2.0 * tr_xxx_xyz[i] * tke_0 - 2.0 * tr_xxxz_xy[i] * tbe_0 - 2.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_xz[i] = tr_xxx_x[i] - 2.0 * tr_xxx_xzz[i] * tke_0 - 2.0 * tr_xxxz_xz[i] * tbe_0 - 2.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_x[i] * tbe_0 + 4.0 * tr_xxxyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_yy[i] = -2.0 * tr_xxx_yyz[i] * tke_0 - 2.0 * tr_xxxz_yy[i] * tbe_0 - 4.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_yz[i] = tr_xxx_y[i] - 2.0 * tr_xxx_yzz[i] * tke_0 - 2.0 * tr_xxxz_yz[i] * tbe_0 + tr_xxxy_0[i] - 2.0 * tr_xxxy_zz[i] * tke_0 - 2.0 * tr_xxxy_yy[i] * tke_0 + 4.0 * tr_xxxy_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxyy_y[i] * tbe_0 + 4.0 * tr_xxxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxy_zz[i] = 2.0 * tr_xxx_z[i] - 2.0 * tr_xxx_zzz[i] * tke_0 - 2.0 * tr_xxxz_zz[i] * tbe_0 - 4.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxyy_z[i] * tbe_0 + 4.0 * tr_xxxyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 372-378 components of targeted buffer : GD

    auto tr_0_0_yz_xxxz_xx = pbuffer.data(idx_op_geom_020_gd + 372);

    auto tr_0_0_yz_xxxz_xy = pbuffer.data(idx_op_geom_020_gd + 373);

    auto tr_0_0_yz_xxxz_xz = pbuffer.data(idx_op_geom_020_gd + 374);

    auto tr_0_0_yz_xxxz_yy = pbuffer.data(idx_op_geom_020_gd + 375);

    auto tr_0_0_yz_xxxz_yz = pbuffer.data(idx_op_geom_020_gd + 376);

    auto tr_0_0_yz_xxxz_zz = pbuffer.data(idx_op_geom_020_gd + 377);

    #pragma omp simd aligned(tr_0_0_yz_xxxz_xx, tr_0_0_yz_xxxz_xy, tr_0_0_yz_xxxz_xz, tr_0_0_yz_xxxz_yy, tr_0_0_yz_xxxz_yz, tr_0_0_yz_xxxz_zz, tr_xxx_x, tr_xxx_xxy, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_y, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxxy_xx, tr_xxxy_xy, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxxyz_x, tr_xxxyz_xxz, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_y, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_z, tr_xxxyz_zzz, tr_xxxyzz_xx, tr_xxxyzz_xy, tr_xxxyzz_xz, tr_xxxyzz_yy, tr_xxxyzz_yz, tr_xxxyzz_zz, tr_xxxz_0, tr_xxxz_xxyz, tr_xxxz_xy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxzz_x, tr_xxxzz_xxy, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_y, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxxz_xx[i] = -2.0 * tr_xxx_xxy[i] * tke_0 + 4.0 * tr_xxxz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xy[i] = tr_xxx_x[i] - 2.0 * tr_xxx_xyy[i] * tke_0 - 2.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_x[i] * tbe_0 + 4.0 * tr_xxxzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_xz[i] = -2.0 * tr_xxx_xyz[i] * tke_0 - 2.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xz[i] * tbe_0 - 2.0 * tr_xxxyz_x[i] * tbe_0 + 4.0 * tr_xxxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_yy[i] = 2.0 * tr_xxx_y[i] - 2.0 * tr_xxx_yyy[i] * tke_0 - 4.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_y[i] * tbe_0 + 4.0 * tr_xxxzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_yz[i] = tr_xxx_z[i] - 2.0 * tr_xxx_yyz[i] * tke_0 + tr_xxxz_0[i] - 2.0 * tr_xxxz_zz[i] * tke_0 - 2.0 * tr_xxxz_yy[i] * tke_0 + 4.0 * tr_xxxz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxxzz_z[i] * tbe_0 + 4.0 * tr_xxxzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yz[i] * tbe_0 - 2.0 * tr_xxxyz_y[i] * tbe_0 + 4.0 * tr_xxxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxxz_zz[i] = -2.0 * tr_xxx_yzz[i] * tke_0 - 4.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxxzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_zz[i] * tbe_0 - 4.0 * tr_xxxyz_z[i] * tbe_0 + 4.0 * tr_xxxyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 378-384 components of targeted buffer : GD

    auto tr_0_0_yz_xxyy_xx = pbuffer.data(idx_op_geom_020_gd + 378);

    auto tr_0_0_yz_xxyy_xy = pbuffer.data(idx_op_geom_020_gd + 379);

    auto tr_0_0_yz_xxyy_xz = pbuffer.data(idx_op_geom_020_gd + 380);

    auto tr_0_0_yz_xxyy_yy = pbuffer.data(idx_op_geom_020_gd + 381);

    auto tr_0_0_yz_xxyy_yz = pbuffer.data(idx_op_geom_020_gd + 382);

    auto tr_0_0_yz_xxyy_zz = pbuffer.data(idx_op_geom_020_gd + 383);

    #pragma omp simd aligned(tr_0_0_yz_xxyy_xx, tr_0_0_yz_xxyy_xy, tr_0_0_yz_xxyy_xz, tr_0_0_yz_xxyy_yy, tr_0_0_yz_xxyy_yz, tr_0_0_yz_xxyy_zz, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyy_0, tr_xxyy_xxyz, tr_xxyy_xy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyyy_x, tr_xxyyy_xxz, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_y, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_z, tr_xxyyy_zzz, tr_xxyyyz_xx, tr_xxyyyz_xy, tr_xxyyyz_xz, tr_xxyyyz_yy, tr_xxyyyz_yz, tr_xxyyyz_zz, tr_xxyyz_x, tr_xxyyz_xxy, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_y, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_z, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxyy_xx[i] = -4.0 * tr_xxy_xxz[i] * tke_0 - 4.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xy[i] = -4.0 * tr_xxy_xyz[i] * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 - 2.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_xz[i] = 2.0 * tr_xxy_x[i] - 4.0 * tr_xxy_xzz[i] * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 - 2.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_x[i] * tbe_0 + 4.0 * tr_xxyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_yy[i] = -4.0 * tr_xxy_yyz[i] * tke_0 - 4.0 * tr_xxyz_yy[i] * tbe_0 - 4.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_yz[i] = 2.0 * tr_xxy_y[i] - 4.0 * tr_xxy_yzz[i] * tke_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 + tr_xxyy_0[i] - 2.0 * tr_xxyy_zz[i] * tke_0 - 2.0 * tr_xxyy_yy[i] * tke_0 + 4.0 * tr_xxyy_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyyy_y[i] * tbe_0 + 4.0 * tr_xxyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyy_zz[i] = 4.0 * tr_xxy_z[i] - 4.0 * tr_xxy_zzz[i] * tke_0 - 4.0 * tr_xxyz_zz[i] * tbe_0 - 4.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyyy_z[i] * tbe_0 + 4.0 * tr_xxyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 384-390 components of targeted buffer : GD

    auto tr_0_0_yz_xxyz_xx = pbuffer.data(idx_op_geom_020_gd + 384);

    auto tr_0_0_yz_xxyz_xy = pbuffer.data(idx_op_geom_020_gd + 385);

    auto tr_0_0_yz_xxyz_xz = pbuffer.data(idx_op_geom_020_gd + 386);

    auto tr_0_0_yz_xxyz_yy = pbuffer.data(idx_op_geom_020_gd + 387);

    auto tr_0_0_yz_xxyz_yz = pbuffer.data(idx_op_geom_020_gd + 388);

    auto tr_0_0_yz_xxyz_zz = pbuffer.data(idx_op_geom_020_gd + 389);

    #pragma omp simd aligned(tr_0_0_yz_xxyz_xx, tr_0_0_yz_xxyz_xy, tr_0_0_yz_xxyz_xz, tr_0_0_yz_xxyz_yy, tr_0_0_yz_xxyz_yz, tr_0_0_yz_xxyz_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxy_x, tr_xxy_xxy, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_y, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxyy_xx, tr_xxyy_xy, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xxyyz_x, tr_xxyyz_xxz, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_y, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_z, tr_xxyyz_zzz, tr_xxyyzz_xx, tr_xxyyzz_xy, tr_xxyyzz_xz, tr_xxyyzz_yy, tr_xxyyzz_yz, tr_xxyyzz_zz, tr_xxyz_0, tr_xxyz_xxyz, tr_xxyz_xy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_xxy, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_y, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_z, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xy, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxyz_xx[i] = tr_xx_xx[i] - 2.0 * tr_xxz_xxz[i] * tke_0 - 2.0 * tr_xxzz_xx[i] * tbe_0 - 2.0 * tr_xxy_xxy[i] * tke_0 + 4.0 * tr_xxyz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xy[i] = tr_xx_xy[i] - 2.0 * tr_xxz_xyz[i] * tke_0 - 2.0 * tr_xxzz_xy[i] * tbe_0 + tr_xxy_x[i] - 2.0 * tr_xxy_xyy[i] * tke_0 - 2.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_xz[i] = tr_xx_xz[i] + tr_xxz_x[i] - 2.0 * tr_xxz_xzz[i] * tke_0 - 2.0 * tr_xxzz_xz[i] * tbe_0 - 2.0 * tr_xxy_xyz[i] * tke_0 - 2.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xz[i] * tbe_0 - 2.0 * tr_xxyyz_x[i] * tbe_0 + 4.0 * tr_xxyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_yy[i] = tr_xx_yy[i] - 2.0 * tr_xxz_yyz[i] * tke_0 - 2.0 * tr_xxzz_yy[i] * tbe_0 + 2.0 * tr_xxy_y[i] - 2.0 * tr_xxy_yyy[i] * tke_0 - 4.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_yz[i] = tr_xx_yz[i] + tr_xxz_y[i] - 2.0 * tr_xxz_yzz[i] * tke_0 - 2.0 * tr_xxzz_yz[i] * tbe_0 + tr_xxy_z[i] - 2.0 * tr_xxy_yyz[i] * tke_0 + tr_xxyz_0[i] - 2.0 * tr_xxyz_zz[i] * tke_0 - 2.0 * tr_xxyz_yy[i] * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yz[i] * tbe_0 - 2.0 * tr_xxyyz_y[i] * tbe_0 + 4.0 * tr_xxyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxyz_zz[i] = tr_xx_zz[i] + 2.0 * tr_xxz_z[i] - 2.0 * tr_xxz_zzz[i] * tke_0 - 2.0 * tr_xxzz_zz[i] * tbe_0 - 2.0 * tr_xxy_yzz[i] * tke_0 - 4.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_zz[i] * tbe_0 - 4.0 * tr_xxyyz_z[i] * tbe_0 + 4.0 * tr_xxyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-396 components of targeted buffer : GD

    auto tr_0_0_yz_xxzz_xx = pbuffer.data(idx_op_geom_020_gd + 390);

    auto tr_0_0_yz_xxzz_xy = pbuffer.data(idx_op_geom_020_gd + 391);

    auto tr_0_0_yz_xxzz_xz = pbuffer.data(idx_op_geom_020_gd + 392);

    auto tr_0_0_yz_xxzz_yy = pbuffer.data(idx_op_geom_020_gd + 393);

    auto tr_0_0_yz_xxzz_yz = pbuffer.data(idx_op_geom_020_gd + 394);

    auto tr_0_0_yz_xxzz_zz = pbuffer.data(idx_op_geom_020_gd + 395);

    #pragma omp simd aligned(tr_0_0_yz_xxzz_xx, tr_0_0_yz_xxzz_xy, tr_0_0_yz_xxzz_xz, tr_0_0_yz_xxzz_yy, tr_0_0_yz_xxzz_yz, tr_0_0_yz_xxzz_zz, tr_xxyz_xx, tr_xxyz_xy, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xxyzz_x, tr_xxyzz_xxz, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_y, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_z, tr_xxyzz_zzz, tr_xxyzzz_xx, tr_xxyzzz_xy, tr_xxyzzz_xz, tr_xxyzzz_yy, tr_xxyzzz_yz, tr_xxyzzz_zz, tr_xxz_x, tr_xxz_xxy, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_y, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxzz_0, tr_xxzz_xxyz, tr_xxzz_xy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzzz_x, tr_xxzzz_xxy, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_y, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xxzz_xx[i] = -4.0 * tr_xxz_xxy[i] * tke_0 + 4.0 * tr_xxzz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xy[i] = 2.0 * tr_xxz_x[i] - 4.0 * tr_xxz_xyy[i] * tke_0 - 2.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_x[i] * tbe_0 + 4.0 * tr_xxzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_xz[i] = -4.0 * tr_xxz_xyz[i] * tke_0 - 2.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xz[i] * tbe_0 - 2.0 * tr_xxyzz_x[i] * tbe_0 + 4.0 * tr_xxyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_yy[i] = 4.0 * tr_xxz_y[i] - 4.0 * tr_xxz_yyy[i] * tke_0 - 4.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_y[i] * tbe_0 + 4.0 * tr_xxzzz_yyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_yz[i] = 2.0 * tr_xxz_z[i] - 4.0 * tr_xxz_yyz[i] * tke_0 + tr_xxzz_0[i] - 2.0 * tr_xxzz_zz[i] * tke_0 - 2.0 * tr_xxzz_yy[i] * tke_0 + 4.0 * tr_xxzz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xxzzz_z[i] * tbe_0 + 4.0 * tr_xxzzz_yyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yz[i] * tbe_0 - 2.0 * tr_xxyzz_y[i] * tbe_0 + 4.0 * tr_xxyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xxzz_zz[i] = -4.0 * tr_xxz_yzz[i] * tke_0 - 4.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xxzzz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zz[i] * tbe_0 - 4.0 * tr_xxyzz_z[i] * tbe_0 + 4.0 * tr_xxyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 396-402 components of targeted buffer : GD

    auto tr_0_0_yz_xyyy_xx = pbuffer.data(idx_op_geom_020_gd + 396);

    auto tr_0_0_yz_xyyy_xy = pbuffer.data(idx_op_geom_020_gd + 397);

    auto tr_0_0_yz_xyyy_xz = pbuffer.data(idx_op_geom_020_gd + 398);

    auto tr_0_0_yz_xyyy_yy = pbuffer.data(idx_op_geom_020_gd + 399);

    auto tr_0_0_yz_xyyy_yz = pbuffer.data(idx_op_geom_020_gd + 400);

    auto tr_0_0_yz_xyyy_zz = pbuffer.data(idx_op_geom_020_gd + 401);

    #pragma omp simd aligned(tr_0_0_yz_xyyy_xx, tr_0_0_yz_xyyy_xy, tr_0_0_yz_xyyy_xz, tr_0_0_yz_xyyy_yy, tr_0_0_yz_xyyy_yz, tr_0_0_yz_xyyy_zz, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyy_0, tr_xyyy_xxyz, tr_xyyy_xy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyyy_x, tr_xyyyy_xxz, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_y, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_z, tr_xyyyy_zzz, tr_xyyyyz_xx, tr_xyyyyz_xy, tr_xyyyyz_xz, tr_xyyyyz_yy, tr_xyyyyz_yz, tr_xyyyyz_zz, tr_xyyyz_x, tr_xyyyz_xxy, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_y, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_z, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyyy_xx[i] = -6.0 * tr_xyy_xxz[i] * tke_0 - 6.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xy[i] = -6.0 * tr_xyy_xyz[i] * tke_0 - 6.0 * tr_xyyz_xy[i] * tbe_0 - 2.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_xz[i] = 3.0 * tr_xyy_x[i] - 6.0 * tr_xyy_xzz[i] * tke_0 - 6.0 * tr_xyyz_xz[i] * tbe_0 - 2.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_x[i] * tbe_0 + 4.0 * tr_xyyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_yy[i] = -6.0 * tr_xyy_yyz[i] * tke_0 - 6.0 * tr_xyyz_yy[i] * tbe_0 - 4.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_yz[i] = 3.0 * tr_xyy_y[i] - 6.0 * tr_xyy_yzz[i] * tke_0 - 6.0 * tr_xyyz_yz[i] * tbe_0 + tr_xyyy_0[i] - 2.0 * tr_xyyy_zz[i] * tke_0 - 2.0 * tr_xyyy_yy[i] * tke_0 + 4.0 * tr_xyyy_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyyy_y[i] * tbe_0 + 4.0 * tr_xyyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyy_zz[i] = 6.0 * tr_xyy_z[i] - 6.0 * tr_xyy_zzz[i] * tke_0 - 6.0 * tr_xyyz_zz[i] * tbe_0 - 4.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyyy_z[i] * tbe_0 + 4.0 * tr_xyyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 402-408 components of targeted buffer : GD

    auto tr_0_0_yz_xyyz_xx = pbuffer.data(idx_op_geom_020_gd + 402);

    auto tr_0_0_yz_xyyz_xy = pbuffer.data(idx_op_geom_020_gd + 403);

    auto tr_0_0_yz_xyyz_xz = pbuffer.data(idx_op_geom_020_gd + 404);

    auto tr_0_0_yz_xyyz_yy = pbuffer.data(idx_op_geom_020_gd + 405);

    auto tr_0_0_yz_xyyz_yz = pbuffer.data(idx_op_geom_020_gd + 406);

    auto tr_0_0_yz_xyyz_zz = pbuffer.data(idx_op_geom_020_gd + 407);

    #pragma omp simd aligned(tr_0_0_yz_xyyz_xx, tr_0_0_yz_xyyz_xy, tr_0_0_yz_xyyz_xz, tr_0_0_yz_xyyz_yy, tr_0_0_yz_xyyz_yz, tr_0_0_yz_xyyz_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyy_x, tr_xyy_xxy, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_y, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyyy_xx, tr_xyyy_xy, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_xyyyz_x, tr_xyyyz_xxz, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_y, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_z, tr_xyyyz_zzz, tr_xyyyzz_xx, tr_xyyyzz_xy, tr_xyyyzz_xz, tr_xyyyzz_yy, tr_xyyyzz_yz, tr_xyyyzz_zz, tr_xyyz_0, tr_xyyz_xxyz, tr_xyyz_xy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_xxy, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_y, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_z, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyyz_xx[i] = 2.0 * tr_xy_xx[i] - 4.0 * tr_xyz_xxz[i] * tke_0 - 4.0 * tr_xyzz_xx[i] * tbe_0 - 2.0 * tr_xyy_xxy[i] * tke_0 + 4.0 * tr_xyyz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xy[i] = 2.0 * tr_xy_xy[i] - 4.0 * tr_xyz_xyz[i] * tke_0 - 4.0 * tr_xyzz_xy[i] * tbe_0 + tr_xyy_x[i] - 2.0 * tr_xyy_xyy[i] * tke_0 - 2.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_xz[i] = 2.0 * tr_xy_xz[i] + 2.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xzz[i] * tke_0 - 4.0 * tr_xyzz_xz[i] * tbe_0 - 2.0 * tr_xyy_xyz[i] * tke_0 - 2.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xz[i] * tbe_0 - 2.0 * tr_xyyyz_x[i] * tbe_0 + 4.0 * tr_xyyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_yy[i] = 2.0 * tr_xy_yy[i] - 4.0 * tr_xyz_yyz[i] * tke_0 - 4.0 * tr_xyzz_yy[i] * tbe_0 + 2.0 * tr_xyy_y[i] - 2.0 * tr_xyy_yyy[i] * tke_0 - 4.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_yz[i] = 2.0 * tr_xy_yz[i] + 2.0 * tr_xyz_y[i] - 4.0 * tr_xyz_yzz[i] * tke_0 - 4.0 * tr_xyzz_yz[i] * tbe_0 + tr_xyy_z[i] - 2.0 * tr_xyy_yyz[i] * tke_0 + tr_xyyz_0[i] - 2.0 * tr_xyyz_zz[i] * tke_0 - 2.0 * tr_xyyz_yy[i] * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yz[i] * tbe_0 - 2.0 * tr_xyyyz_y[i] * tbe_0 + 4.0 * tr_xyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyyz_zz[i] = 2.0 * tr_xy_zz[i] + 4.0 * tr_xyz_z[i] - 4.0 * tr_xyz_zzz[i] * tke_0 - 4.0 * tr_xyzz_zz[i] * tbe_0 - 2.0 * tr_xyy_yzz[i] * tke_0 - 4.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_zz[i] * tbe_0 - 4.0 * tr_xyyyz_z[i] * tbe_0 + 4.0 * tr_xyyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 408-414 components of targeted buffer : GD

    auto tr_0_0_yz_xyzz_xx = pbuffer.data(idx_op_geom_020_gd + 408);

    auto tr_0_0_yz_xyzz_xy = pbuffer.data(idx_op_geom_020_gd + 409);

    auto tr_0_0_yz_xyzz_xz = pbuffer.data(idx_op_geom_020_gd + 410);

    auto tr_0_0_yz_xyzz_yy = pbuffer.data(idx_op_geom_020_gd + 411);

    auto tr_0_0_yz_xyzz_yz = pbuffer.data(idx_op_geom_020_gd + 412);

    auto tr_0_0_yz_xyzz_zz = pbuffer.data(idx_op_geom_020_gd + 413);

    #pragma omp simd aligned(tr_0_0_yz_xyzz_xx, tr_0_0_yz_xyzz_xy, tr_0_0_yz_xyzz_xz, tr_0_0_yz_xyzz_yy, tr_0_0_yz_xyzz_yz, tr_0_0_yz_xyzz_zz, tr_xyyz_xx, tr_xyyz_xy, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_xyyzz_x, tr_xyyzz_xxz, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_y, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_z, tr_xyyzz_zzz, tr_xyyzzz_xx, tr_xyyzzz_xy, tr_xyyzzz_xz, tr_xyyzzz_yy, tr_xyyzzz_yz, tr_xyyzzz_zz, tr_xyz_x, tr_xyz_xxy, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_y, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyzz_0, tr_xyzz_xxyz, tr_xyzz_xy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_xxy, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_y, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_z, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xy, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xyzz_xx[i] = 2.0 * tr_xz_xx[i] - 2.0 * tr_xzz_xxz[i] * tke_0 - 2.0 * tr_xzzz_xx[i] * tbe_0 - 4.0 * tr_xyz_xxy[i] * tke_0 + 4.0 * tr_xyzz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xy[i] = 2.0 * tr_xz_xy[i] - 2.0 * tr_xzz_xyz[i] * tke_0 - 2.0 * tr_xzzz_xy[i] * tbe_0 + 2.0 * tr_xyz_x[i] - 4.0 * tr_xyz_xyy[i] * tke_0 - 2.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_xz[i] = 2.0 * tr_xz_xz[i] + tr_xzz_x[i] - 2.0 * tr_xzz_xzz[i] * tke_0 - 2.0 * tr_xzzz_xz[i] * tbe_0 - 4.0 * tr_xyz_xyz[i] * tke_0 - 2.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xz[i] * tbe_0 - 2.0 * tr_xyyzz_x[i] * tbe_0 + 4.0 * tr_xyyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_yy[i] = 2.0 * tr_xz_yy[i] - 2.0 * tr_xzz_yyz[i] * tke_0 - 2.0 * tr_xzzz_yy[i] * tbe_0 + 4.0 * tr_xyz_y[i] - 4.0 * tr_xyz_yyy[i] * tke_0 - 4.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_yz[i] = 2.0 * tr_xz_yz[i] + tr_xzz_y[i] - 2.0 * tr_xzz_yzz[i] * tke_0 - 2.0 * tr_xzzz_yz[i] * tbe_0 + 2.0 * tr_xyz_z[i] - 4.0 * tr_xyz_yyz[i] * tke_0 + tr_xyzz_0[i] - 2.0 * tr_xyzz_zz[i] * tke_0 - 2.0 * tr_xyzz_yy[i] * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xyzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yz[i] * tbe_0 - 2.0 * tr_xyyzz_y[i] * tbe_0 + 4.0 * tr_xyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xyzz_zz[i] = 2.0 * tr_xz_zz[i] + 2.0 * tr_xzz_z[i] - 2.0 * tr_xzz_zzz[i] * tke_0 - 2.0 * tr_xzzz_zz[i] * tbe_0 - 4.0 * tr_xyz_yzz[i] * tke_0 - 4.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_zz[i] * tbe_0 - 4.0 * tr_xyyzz_z[i] * tbe_0 + 4.0 * tr_xyyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 414-420 components of targeted buffer : GD

    auto tr_0_0_yz_xzzz_xx = pbuffer.data(idx_op_geom_020_gd + 414);

    auto tr_0_0_yz_xzzz_xy = pbuffer.data(idx_op_geom_020_gd + 415);

    auto tr_0_0_yz_xzzz_xz = pbuffer.data(idx_op_geom_020_gd + 416);

    auto tr_0_0_yz_xzzz_yy = pbuffer.data(idx_op_geom_020_gd + 417);

    auto tr_0_0_yz_xzzz_yz = pbuffer.data(idx_op_geom_020_gd + 418);

    auto tr_0_0_yz_xzzz_zz = pbuffer.data(idx_op_geom_020_gd + 419);

    #pragma omp simd aligned(tr_0_0_yz_xzzz_xx, tr_0_0_yz_xzzz_xy, tr_0_0_yz_xzzz_xz, tr_0_0_yz_xzzz_yy, tr_0_0_yz_xzzz_yz, tr_0_0_yz_xzzz_zz, tr_xyzz_xx, tr_xyzz_xy, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_xyzzz_x, tr_xyzzz_xxz, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_y, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_z, tr_xyzzz_zzz, tr_xyzzzz_xx, tr_xyzzzz_xy, tr_xyzzzz_xz, tr_xyzzzz_yy, tr_xyzzzz_yz, tr_xyzzzz_zz, tr_xzz_x, tr_xzz_xxy, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_y, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzzz_0, tr_xzzz_xxyz, tr_xzzz_xy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzzz_x, tr_xzzzz_xxy, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_y, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_xzzz_xx[i] = -6.0 * tr_xzz_xxy[i] * tke_0 + 4.0 * tr_xzzz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xy[i] = 3.0 * tr_xzz_x[i] - 6.0 * tr_xzz_xyy[i] * tke_0 - 2.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_x[i] * tbe_0 + 4.0 * tr_xzzzz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_xz[i] = -6.0 * tr_xzz_xyz[i] * tke_0 - 2.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xz[i] * tbe_0 - 2.0 * tr_xyzzz_x[i] * tbe_0 + 4.0 * tr_xyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_yy[i] = 6.0 * tr_xzz_y[i] - 6.0 * tr_xzz_yyy[i] * tke_0 - 4.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_y[i] * tbe_0 + 4.0 * tr_xzzzz_yyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_yz[i] = 3.0 * tr_xzz_z[i] - 6.0 * tr_xzz_yyz[i] * tke_0 + tr_xzzz_0[i] - 2.0 * tr_xzzz_zz[i] * tke_0 - 2.0 * tr_xzzz_yy[i] * tke_0 + 4.0 * tr_xzzz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_xzzzz_z[i] * tbe_0 + 4.0 * tr_xzzzz_yyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yz[i] * tbe_0 - 2.0 * tr_xyzzz_y[i] * tbe_0 + 4.0 * tr_xyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_xzzz_zz[i] = -6.0 * tr_xzz_yzz[i] * tke_0 - 4.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_xzzzz_yzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_zz[i] * tbe_0 - 4.0 * tr_xyzzz_z[i] * tbe_0 + 4.0 * tr_xyzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-426 components of targeted buffer : GD

    auto tr_0_0_yz_yyyy_xx = pbuffer.data(idx_op_geom_020_gd + 420);

    auto tr_0_0_yz_yyyy_xy = pbuffer.data(idx_op_geom_020_gd + 421);

    auto tr_0_0_yz_yyyy_xz = pbuffer.data(idx_op_geom_020_gd + 422);

    auto tr_0_0_yz_yyyy_yy = pbuffer.data(idx_op_geom_020_gd + 423);

    auto tr_0_0_yz_yyyy_yz = pbuffer.data(idx_op_geom_020_gd + 424);

    auto tr_0_0_yz_yyyy_zz = pbuffer.data(idx_op_geom_020_gd + 425);

    #pragma omp simd aligned(tr_0_0_yz_yyyy_xx, tr_0_0_yz_yyyy_xy, tr_0_0_yz_yyyy_xz, tr_0_0_yz_yyyy_yy, tr_0_0_yz_yyyy_yz, tr_0_0_yz_yyyy_zz, tr_yyy_x, tr_yyy_xxz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyy_0, tr_yyyy_xxyz, tr_yyyy_xy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyyy_x, tr_yyyyy_xxz, tr_yyyyy_xyz, tr_yyyyy_xzz, tr_yyyyy_y, tr_yyyyy_yyz, tr_yyyyy_yzz, tr_yyyyy_z, tr_yyyyy_zzz, tr_yyyyyz_xx, tr_yyyyyz_xy, tr_yyyyyz_xz, tr_yyyyyz_yy, tr_yyyyyz_yz, tr_yyyyyz_zz, tr_yyyyz_x, tr_yyyyz_xxy, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_y, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_z, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyyy_xx[i] = -8.0 * tr_yyy_xxz[i] * tke_0 - 8.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_yyyy_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xy[i] = -8.0 * tr_yyy_xyz[i] * tke_0 - 8.0 * tr_yyyz_xy[i] * tbe_0 - 2.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_x[i] * tbe_0 + 4.0 * tr_yyyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_xz[i] = 4.0 * tr_yyy_x[i] - 8.0 * tr_yyy_xzz[i] * tke_0 - 8.0 * tr_yyyz_xz[i] * tbe_0 - 2.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyyy_x[i] * tbe_0 + 4.0 * tr_yyyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_yy[i] = -8.0 * tr_yyy_yyz[i] * tke_0 - 8.0 * tr_yyyz_yy[i] * tbe_0 - 4.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_y[i] * tbe_0 + 4.0 * tr_yyyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_yz[i] = 4.0 * tr_yyy_y[i] - 8.0 * tr_yyy_yzz[i] * tke_0 - 8.0 * tr_yyyz_yz[i] * tbe_0 + tr_yyyy_0[i] - 2.0 * tr_yyyy_zz[i] * tke_0 - 2.0 * tr_yyyy_yy[i] * tke_0 + 4.0 * tr_yyyy_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyyz_z[i] * tbe_0 + 4.0 * tr_yyyyz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyyy_y[i] * tbe_0 + 4.0 * tr_yyyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyy_zz[i] = 8.0 * tr_yyy_z[i] - 8.0 * tr_yyy_zzz[i] * tke_0 - 8.0 * tr_yyyz_zz[i] * tbe_0 - 4.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyyz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyyy_z[i] * tbe_0 + 4.0 * tr_yyyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 426-432 components of targeted buffer : GD

    auto tr_0_0_yz_yyyz_xx = pbuffer.data(idx_op_geom_020_gd + 426);

    auto tr_0_0_yz_yyyz_xy = pbuffer.data(idx_op_geom_020_gd + 427);

    auto tr_0_0_yz_yyyz_xz = pbuffer.data(idx_op_geom_020_gd + 428);

    auto tr_0_0_yz_yyyz_yy = pbuffer.data(idx_op_geom_020_gd + 429);

    auto tr_0_0_yz_yyyz_yz = pbuffer.data(idx_op_geom_020_gd + 430);

    auto tr_0_0_yz_yyyz_zz = pbuffer.data(idx_op_geom_020_gd + 431);

    #pragma omp simd aligned(tr_0_0_yz_yyyz_xx, tr_0_0_yz_yyyz_xy, tr_0_0_yz_yyyz_xz, tr_0_0_yz_yyyz_yy, tr_0_0_yz_yyyz_yz, tr_0_0_yz_yyyz_zz, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyy_x, tr_yyy_xxy, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_y, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyyy_xx, tr_yyyy_xy, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, tr_yyyyz_x, tr_yyyyz_xxz, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_y, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_z, tr_yyyyz_zzz, tr_yyyyzz_xx, tr_yyyyzz_xy, tr_yyyyzz_xz, tr_yyyyzz_yy, tr_yyyyzz_yz, tr_yyyyzz_zz, tr_yyyz_0, tr_yyyz_xxyz, tr_yyyz_xy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_xxy, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_y, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_z, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyyz_xx[i] = 3.0 * tr_yy_xx[i] - 6.0 * tr_yyz_xxz[i] * tke_0 - 6.0 * tr_yyzz_xx[i] * tbe_0 - 2.0 * tr_yyy_xxy[i] * tke_0 + 4.0 * tr_yyyz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xy[i] = 3.0 * tr_yy_xy[i] - 6.0 * tr_yyz_xyz[i] * tke_0 - 6.0 * tr_yyzz_xy[i] * tbe_0 + tr_yyy_x[i] - 2.0 * tr_yyy_xyy[i] * tke_0 - 2.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_xz[i] = 3.0 * tr_yy_xz[i] + 3.0 * tr_yyz_x[i] - 6.0 * tr_yyz_xzz[i] * tke_0 - 6.0 * tr_yyzz_xz[i] * tbe_0 - 2.0 * tr_yyy_xyz[i] * tke_0 - 2.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xz[i] * tbe_0 - 2.0 * tr_yyyyz_x[i] * tbe_0 + 4.0 * tr_yyyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_yy[i] = 3.0 * tr_yy_yy[i] - 6.0 * tr_yyz_yyz[i] * tke_0 - 6.0 * tr_yyzz_yy[i] * tbe_0 + 2.0 * tr_yyy_y[i] - 2.0 * tr_yyy_yyy[i] * tke_0 - 4.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_y[i] * tbe_0 + 4.0 * tr_yyyzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yy[i] * tbe_0 + 4.0 * tr_yyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_yz[i] = 3.0 * tr_yy_yz[i] + 3.0 * tr_yyz_y[i] - 6.0 * tr_yyz_yzz[i] * tke_0 - 6.0 * tr_yyzz_yz[i] * tbe_0 + tr_yyy_z[i] - 2.0 * tr_yyy_yyz[i] * tke_0 + tr_yyyz_0[i] - 2.0 * tr_yyyz_zz[i] * tke_0 - 2.0 * tr_yyyz_yy[i] * tke_0 + 4.0 * tr_yyyz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyyzz_z[i] * tbe_0 + 4.0 * tr_yyyzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yz[i] * tbe_0 - 2.0 * tr_yyyyz_y[i] * tbe_0 + 4.0 * tr_yyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyyz_zz[i] = 3.0 * tr_yy_zz[i] + 6.0 * tr_yyz_z[i] - 6.0 * tr_yyz_zzz[i] * tke_0 - 6.0 * tr_yyzz_zz[i] * tbe_0 - 2.0 * tr_yyy_yzz[i] * tke_0 - 4.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyyzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_zz[i] * tbe_0 - 4.0 * tr_yyyyz_z[i] * tbe_0 + 4.0 * tr_yyyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 432-438 components of targeted buffer : GD

    auto tr_0_0_yz_yyzz_xx = pbuffer.data(idx_op_geom_020_gd + 432);

    auto tr_0_0_yz_yyzz_xy = pbuffer.data(idx_op_geom_020_gd + 433);

    auto tr_0_0_yz_yyzz_xz = pbuffer.data(idx_op_geom_020_gd + 434);

    auto tr_0_0_yz_yyzz_yy = pbuffer.data(idx_op_geom_020_gd + 435);

    auto tr_0_0_yz_yyzz_yz = pbuffer.data(idx_op_geom_020_gd + 436);

    auto tr_0_0_yz_yyzz_zz = pbuffer.data(idx_op_geom_020_gd + 437);

    #pragma omp simd aligned(tr_0_0_yz_yyzz_xx, tr_0_0_yz_yyzz_xy, tr_0_0_yz_yyzz_xz, tr_0_0_yz_yyzz_yy, tr_0_0_yz_yyzz_yz, tr_0_0_yz_yyzz_zz, tr_yyyz_xx, tr_yyyz_xy, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, tr_yyyzz_x, tr_yyyzz_xxz, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_y, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_z, tr_yyyzz_zzz, tr_yyyzzz_xx, tr_yyyzzz_xy, tr_yyyzzz_xz, tr_yyyzzz_yy, tr_yyyzzz_yz, tr_yyyzzz_zz, tr_yyz_x, tr_yyz_xxy, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_y, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyzz_0, tr_yyzz_xxyz, tr_yyzz_xy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_xxy, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_y, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_z, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yyzz_xx[i] = 4.0 * tr_yz_xx[i] - 4.0 * tr_yzz_xxz[i] * tke_0 - 4.0 * tr_yzzz_xx[i] * tbe_0 - 4.0 * tr_yyz_xxy[i] * tke_0 + 4.0 * tr_yyzz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xy[i] = 4.0 * tr_yz_xy[i] - 4.0 * tr_yzz_xyz[i] * tke_0 - 4.0 * tr_yzzz_xy[i] * tbe_0 + 2.0 * tr_yyz_x[i] - 4.0 * tr_yyz_xyy[i] * tke_0 - 2.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_xz[i] = 4.0 * tr_yz_xz[i] + 2.0 * tr_yzz_x[i] - 4.0 * tr_yzz_xzz[i] * tke_0 - 4.0 * tr_yzzz_xz[i] * tbe_0 - 4.0 * tr_yyz_xyz[i] * tke_0 - 2.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xz[i] * tbe_0 - 2.0 * tr_yyyzz_x[i] * tbe_0 + 4.0 * tr_yyyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_yy[i] = 4.0 * tr_yz_yy[i] - 4.0 * tr_yzz_yyz[i] * tke_0 - 4.0 * tr_yzzz_yy[i] * tbe_0 + 4.0 * tr_yyz_y[i] - 4.0 * tr_yyz_yyy[i] * tke_0 - 4.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_y[i] * tbe_0 + 4.0 * tr_yyzzz_yyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_yz[i] = 4.0 * tr_yz_yz[i] + 2.0 * tr_yzz_y[i] - 4.0 * tr_yzz_yzz[i] * tke_0 - 4.0 * tr_yzzz_yz[i] * tbe_0 + 2.0 * tr_yyz_z[i] - 4.0 * tr_yyz_yyz[i] * tke_0 + tr_yyzz_0[i] - 2.0 * tr_yyzz_zz[i] * tke_0 - 2.0 * tr_yyzz_yy[i] * tke_0 + 4.0 * tr_yyzz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_yyzzz_z[i] * tbe_0 + 4.0 * tr_yyzzz_yyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yz[i] * tbe_0 - 2.0 * tr_yyyzz_y[i] * tbe_0 + 4.0 * tr_yyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yyzz_zz[i] = 4.0 * tr_yz_zz[i] + 4.0 * tr_yzz_z[i] - 4.0 * tr_yzz_zzz[i] * tke_0 - 4.0 * tr_yzzz_zz[i] * tbe_0 - 4.0 * tr_yyz_yzz[i] * tke_0 - 4.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_yyzzz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_zz[i] * tbe_0 - 4.0 * tr_yyyzz_z[i] * tbe_0 + 4.0 * tr_yyyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 438-444 components of targeted buffer : GD

    auto tr_0_0_yz_yzzz_xx = pbuffer.data(idx_op_geom_020_gd + 438);

    auto tr_0_0_yz_yzzz_xy = pbuffer.data(idx_op_geom_020_gd + 439);

    auto tr_0_0_yz_yzzz_xz = pbuffer.data(idx_op_geom_020_gd + 440);

    auto tr_0_0_yz_yzzz_yy = pbuffer.data(idx_op_geom_020_gd + 441);

    auto tr_0_0_yz_yzzz_yz = pbuffer.data(idx_op_geom_020_gd + 442);

    auto tr_0_0_yz_yzzz_zz = pbuffer.data(idx_op_geom_020_gd + 443);

    #pragma omp simd aligned(tr_0_0_yz_yzzz_xx, tr_0_0_yz_yzzz_xy, tr_0_0_yz_yzzz_xz, tr_0_0_yz_yzzz_yy, tr_0_0_yz_yzzz_yz, tr_0_0_yz_yzzz_zz, tr_yyzz_xx, tr_yyzz_xy, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, tr_yyzzz_x, tr_yyzzz_xxz, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_y, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_z, tr_yyzzz_zzz, tr_yyzzzz_xx, tr_yyzzzz_xy, tr_yyzzzz_xz, tr_yyzzzz_yy, tr_yyzzzz_yz, tr_yyzzzz_zz, tr_yzz_x, tr_yzz_xxy, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_y, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzzz_0, tr_yzzz_xxyz, tr_yzzz_xy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_xxy, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_y, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_z, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xxz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xy, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_yzzz_xx[i] = 3.0 * tr_zz_xx[i] - 2.0 * tr_zzz_xxz[i] * tke_0 - 2.0 * tr_zzzz_xx[i] * tbe_0 - 6.0 * tr_yzz_xxy[i] * tke_0 + 4.0 * tr_yzzz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xy[i] = 3.0 * tr_zz_xy[i] - 2.0 * tr_zzz_xyz[i] * tke_0 - 2.0 * tr_zzzz_xy[i] * tbe_0 + 3.0 * tr_yzz_x[i] - 6.0 * tr_yzz_xyy[i] * tke_0 - 2.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_xz[i] = 3.0 * tr_zz_xz[i] + tr_zzz_x[i] - 2.0 * tr_zzz_xzz[i] * tke_0 - 2.0 * tr_zzzz_xz[i] * tbe_0 - 6.0 * tr_yzz_xyz[i] * tke_0 - 2.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xz[i] * tbe_0 - 2.0 * tr_yyzzz_x[i] * tbe_0 + 4.0 * tr_yyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_yy[i] = 3.0 * tr_zz_yy[i] - 2.0 * tr_zzz_yyz[i] * tke_0 - 2.0 * tr_zzzz_yy[i] * tbe_0 + 6.0 * tr_yzz_y[i] - 6.0 * tr_yzz_yyy[i] * tke_0 - 4.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_yyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_yz[i] = 3.0 * tr_zz_yz[i] + tr_zzz_y[i] - 2.0 * tr_zzz_yzz[i] * tke_0 - 2.0 * tr_zzzz_yz[i] * tbe_0 + 3.0 * tr_yzz_z[i] - 6.0 * tr_yzz_yyz[i] * tke_0 + tr_yzzz_0[i] - 2.0 * tr_yzzz_zz[i] * tke_0 - 2.0 * tr_yzzz_yy[i] * tke_0 + 4.0 * tr_yzzz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_yzzzz_z[i] * tbe_0 + 4.0 * tr_yzzzz_yyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yz[i] * tbe_0 - 2.0 * tr_yyzzz_y[i] * tbe_0 + 4.0 * tr_yyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_yzzz_zz[i] = 3.0 * tr_zz_zz[i] + 2.0 * tr_zzz_z[i] - 2.0 * tr_zzz_zzz[i] * tke_0 - 2.0 * tr_zzzz_zz[i] * tbe_0 - 6.0 * tr_yzz_yzz[i] * tke_0 - 4.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_yzzzz_yzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_zz[i] * tbe_0 - 4.0 * tr_yyzzz_z[i] * tbe_0 + 4.0 * tr_yyzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 444-450 components of targeted buffer : GD

    auto tr_0_0_yz_zzzz_xx = pbuffer.data(idx_op_geom_020_gd + 444);

    auto tr_0_0_yz_zzzz_xy = pbuffer.data(idx_op_geom_020_gd + 445);

    auto tr_0_0_yz_zzzz_xz = pbuffer.data(idx_op_geom_020_gd + 446);

    auto tr_0_0_yz_zzzz_yy = pbuffer.data(idx_op_geom_020_gd + 447);

    auto tr_0_0_yz_zzzz_yz = pbuffer.data(idx_op_geom_020_gd + 448);

    auto tr_0_0_yz_zzzz_zz = pbuffer.data(idx_op_geom_020_gd + 449);

    #pragma omp simd aligned(tr_0_0_yz_zzzz_xx, tr_0_0_yz_zzzz_xy, tr_0_0_yz_zzzz_xz, tr_0_0_yz_zzzz_yy, tr_0_0_yz_zzzz_yz, tr_0_0_yz_zzzz_zz, tr_yzzz_xx, tr_yzzz_xy, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, tr_yzzzz_x, tr_yzzzz_xxz, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_y, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_z, tr_yzzzz_zzz, tr_yzzzzz_xx, tr_yzzzzz_xy, tr_yzzzzz_xz, tr_yzzzzz_yy, tr_yzzzzz_yz, tr_yzzzzz_zz, tr_zzz_x, tr_zzz_xxy, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_y, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzzz_0, tr_zzzz_xxyz, tr_zzzz_xy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzzz_x, tr_zzzzz_xxy, tr_zzzzz_xyy, tr_zzzzz_xyz, tr_zzzzz_y, tr_zzzzz_yyy, tr_zzzzz_yyz, tr_zzzzz_yzz, tr_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_yz_zzzz_xx[i] = -8.0 * tr_zzz_xxy[i] * tke_0 + 4.0 * tr_zzzz_xxyz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xxy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xy[i] = 4.0 * tr_zzz_x[i] - 8.0 * tr_zzz_xyy[i] * tke_0 - 2.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xyyz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_x[i] * tbe_0 + 4.0 * tr_zzzzz_xyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_xz[i] = -8.0 * tr_zzz_xyz[i] * tke_0 - 2.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xyzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_xyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xz[i] * tbe_0 - 2.0 * tr_yzzzz_x[i] * tbe_0 + 4.0 * tr_yzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_yy[i] = 8.0 * tr_zzz_y[i] - 8.0 * tr_zzz_yyy[i] * tke_0 - 4.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_yyyz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_y[i] * tbe_0 + 4.0 * tr_zzzzz_yyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_yz[i] = 4.0 * tr_zzz_z[i] - 8.0 * tr_zzz_yyz[i] * tke_0 + tr_zzzz_0[i] - 2.0 * tr_zzzz_zz[i] * tke_0 - 2.0 * tr_zzzz_yy[i] * tke_0 + 4.0 * tr_zzzz_yyzz[i] * tke_0 * tke_0 - 2.0 * tr_zzzzz_z[i] * tbe_0 + 4.0 * tr_zzzzz_yyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yz[i] * tbe_0 - 2.0 * tr_yzzzz_y[i] * tbe_0 + 4.0 * tr_yzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_yz_zzzz_zz[i] = -8.0 * tr_zzz_yzz[i] * tke_0 - 4.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_yzzz[i] * tke_0 * tke_0 + 4.0 * tr_zzzzz_yzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_zz[i] * tbe_0 - 4.0 * tr_yzzzz_z[i] * tbe_0 + 4.0 * tr_yzzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-456 components of targeted buffer : GD

    auto tr_0_0_zz_xxxx_xx = pbuffer.data(idx_op_geom_020_gd + 450);

    auto tr_0_0_zz_xxxx_xy = pbuffer.data(idx_op_geom_020_gd + 451);

    auto tr_0_0_zz_xxxx_xz = pbuffer.data(idx_op_geom_020_gd + 452);

    auto tr_0_0_zz_xxxx_yy = pbuffer.data(idx_op_geom_020_gd + 453);

    auto tr_0_0_zz_xxxx_yz = pbuffer.data(idx_op_geom_020_gd + 454);

    auto tr_0_0_zz_xxxx_zz = pbuffer.data(idx_op_geom_020_gd + 455);

    #pragma omp simd aligned(tr_0_0_zz_xxxx_xx, tr_0_0_zz_xxxx_xy, tr_0_0_zz_xxxx_xz, tr_0_0_zz_xxxx_yy, tr_0_0_zz_xxxx_yz, tr_0_0_zz_xxxx_zz, tr_xxxx_0, tr_xxxx_xx, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxz_x, tr_xxxxz_xxz, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_y, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_z, tr_xxxxz_zzz, tr_xxxxzz_xx, tr_xxxxzz_xy, tr_xxxxzz_xz, tr_xxxxzz_yy, tr_xxxxzz_yz, tr_xxxxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxxx_xx[i] = -2.0 * tr_xxxx_xx[i] * tbe_0 - 2.0 * tr_xxxx_xx[i] * tke_0 + 4.0 * tr_xxxx_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xy[i] = -2.0 * tr_xxxx_xy[i] * tbe_0 - 2.0 * tr_xxxx_xy[i] * tke_0 + 4.0 * tr_xxxx_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_xz[i] = -2.0 * tr_xxxx_xz[i] * tbe_0 - 6.0 * tr_xxxx_xz[i] * tke_0 + 4.0 * tr_xxxx_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_x[i] * tbe_0 + 8.0 * tr_xxxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_yy[i] = -2.0 * tr_xxxx_yy[i] * tbe_0 - 2.0 * tr_xxxx_yy[i] * tke_0 + 4.0 * tr_xxxx_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_yz[i] = -2.0 * tr_xxxx_yz[i] * tbe_0 - 6.0 * tr_xxxx_yz[i] * tke_0 + 4.0 * tr_xxxx_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxxz_y[i] * tbe_0 + 8.0 * tr_xxxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxx_zz[i] = 2.0 * tr_xxxx_0[i] - 2.0 * tr_xxxx_zz[i] * tbe_0 - 10.0 * tr_xxxx_zz[i] * tke_0 + 4.0 * tr_xxxx_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxxz_z[i] * tbe_0 + 8.0 * tr_xxxxz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 456-462 components of targeted buffer : GD

    auto tr_0_0_zz_xxxy_xx = pbuffer.data(idx_op_geom_020_gd + 456);

    auto tr_0_0_zz_xxxy_xy = pbuffer.data(idx_op_geom_020_gd + 457);

    auto tr_0_0_zz_xxxy_xz = pbuffer.data(idx_op_geom_020_gd + 458);

    auto tr_0_0_zz_xxxy_yy = pbuffer.data(idx_op_geom_020_gd + 459);

    auto tr_0_0_zz_xxxy_yz = pbuffer.data(idx_op_geom_020_gd + 460);

    auto tr_0_0_zz_xxxy_zz = pbuffer.data(idx_op_geom_020_gd + 461);

    #pragma omp simd aligned(tr_0_0_zz_xxxy_xx, tr_0_0_zz_xxxy_xy, tr_0_0_zz_xxxy_xz, tr_0_0_zz_xxxy_yy, tr_0_0_zz_xxxy_yz, tr_0_0_zz_xxxy_zz, tr_xxxy_0, tr_xxxy_xx, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyz_x, tr_xxxyz_xxz, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_y, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_z, tr_xxxyz_zzz, tr_xxxyzz_xx, tr_xxxyzz_xy, tr_xxxyzz_xz, tr_xxxyzz_yy, tr_xxxyzz_yz, tr_xxxyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxxy_xx[i] = -2.0 * tr_xxxy_xx[i] * tbe_0 - 2.0 * tr_xxxy_xx[i] * tke_0 + 4.0 * tr_xxxy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xy[i] = -2.0 * tr_xxxy_xy[i] * tbe_0 - 2.0 * tr_xxxy_xy[i] * tke_0 + 4.0 * tr_xxxy_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_xz[i] = -2.0 * tr_xxxy_xz[i] * tbe_0 - 6.0 * tr_xxxy_xz[i] * tke_0 + 4.0 * tr_xxxy_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_x[i] * tbe_0 + 8.0 * tr_xxxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_yy[i] = -2.0 * tr_xxxy_yy[i] * tbe_0 - 2.0 * tr_xxxy_yy[i] * tke_0 + 4.0 * tr_xxxy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_yz[i] = -2.0 * tr_xxxy_yz[i] * tbe_0 - 6.0 * tr_xxxy_yz[i] * tke_0 + 4.0 * tr_xxxy_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxyz_y[i] * tbe_0 + 8.0 * tr_xxxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxy_zz[i] = 2.0 * tr_xxxy_0[i] - 2.0 * tr_xxxy_zz[i] * tbe_0 - 10.0 * tr_xxxy_zz[i] * tke_0 + 4.0 * tr_xxxy_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxyz_z[i] * tbe_0 + 8.0 * tr_xxxyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 462-468 components of targeted buffer : GD

    auto tr_0_0_zz_xxxz_xx = pbuffer.data(idx_op_geom_020_gd + 462);

    auto tr_0_0_zz_xxxz_xy = pbuffer.data(idx_op_geom_020_gd + 463);

    auto tr_0_0_zz_xxxz_xz = pbuffer.data(idx_op_geom_020_gd + 464);

    auto tr_0_0_zz_xxxz_yy = pbuffer.data(idx_op_geom_020_gd + 465);

    auto tr_0_0_zz_xxxz_yz = pbuffer.data(idx_op_geom_020_gd + 466);

    auto tr_0_0_zz_xxxz_zz = pbuffer.data(idx_op_geom_020_gd + 467);

    #pragma omp simd aligned(tr_0_0_zz_xxxz_xx, tr_0_0_zz_xxxz_xy, tr_0_0_zz_xxxz_xz, tr_0_0_zz_xxxz_yy, tr_0_0_zz_xxxz_yz, tr_0_0_zz_xxxz_zz, tr_xxx_x, tr_xxx_xxz, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_y, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_z, tr_xxx_zzz, tr_xxxz_0, tr_xxxz_xx, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxzz_x, tr_xxxzz_xxz, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_y, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_z, tr_xxxzz_zzz, tr_xxxzzz_xx, tr_xxxzzz_xy, tr_xxxzzz_xz, tr_xxxzzz_yy, tr_xxxzzz_yz, tr_xxxzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxxz_xx[i] = -4.0 * tr_xxx_xxz[i] * tke_0 - 6.0 * tr_xxxz_xx[i] * tbe_0 - 2.0 * tr_xxxz_xx[i] * tke_0 + 4.0 * tr_xxxz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xy[i] = -4.0 * tr_xxx_xyz[i] * tke_0 - 6.0 * tr_xxxz_xy[i] * tbe_0 - 2.0 * tr_xxxz_xy[i] * tke_0 + 4.0 * tr_xxxz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_xz[i] = 2.0 * tr_xxx_x[i] - 4.0 * tr_xxx_xzz[i] * tke_0 - 6.0 * tr_xxxz_xz[i] * tbe_0 - 6.0 * tr_xxxz_xz[i] * tke_0 + 4.0 * tr_xxxz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_x[i] * tbe_0 + 8.0 * tr_xxxzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_yy[i] = -4.0 * tr_xxx_yyz[i] * tke_0 - 6.0 * tr_xxxz_yy[i] * tbe_0 - 2.0 * tr_xxxz_yy[i] * tke_0 + 4.0 * tr_xxxz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxxzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_yz[i] = 2.0 * tr_xxx_y[i] - 4.0 * tr_xxx_yzz[i] * tke_0 - 6.0 * tr_xxxz_yz[i] * tbe_0 - 6.0 * tr_xxxz_yz[i] * tke_0 + 4.0 * tr_xxxz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxxzz_y[i] * tbe_0 + 8.0 * tr_xxxzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxxz_zz[i] = 4.0 * tr_xxx_z[i] - 4.0 * tr_xxx_zzz[i] * tke_0 + 2.0 * tr_xxxz_0[i] - 6.0 * tr_xxxz_zz[i] * tbe_0 - 10.0 * tr_xxxz_zz[i] * tke_0 + 4.0 * tr_xxxz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxxzz_z[i] * tbe_0 + 8.0 * tr_xxxzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 468-474 components of targeted buffer : GD

    auto tr_0_0_zz_xxyy_xx = pbuffer.data(idx_op_geom_020_gd + 468);

    auto tr_0_0_zz_xxyy_xy = pbuffer.data(idx_op_geom_020_gd + 469);

    auto tr_0_0_zz_xxyy_xz = pbuffer.data(idx_op_geom_020_gd + 470);

    auto tr_0_0_zz_xxyy_yy = pbuffer.data(idx_op_geom_020_gd + 471);

    auto tr_0_0_zz_xxyy_yz = pbuffer.data(idx_op_geom_020_gd + 472);

    auto tr_0_0_zz_xxyy_zz = pbuffer.data(idx_op_geom_020_gd + 473);

    #pragma omp simd aligned(tr_0_0_zz_xxyy_xx, tr_0_0_zz_xxyy_xy, tr_0_0_zz_xxyy_xz, tr_0_0_zz_xxyy_yy, tr_0_0_zz_xxyy_yz, tr_0_0_zz_xxyy_zz, tr_xxyy_0, tr_xxyy_xx, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyz_x, tr_xxyyz_xxz, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_y, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_z, tr_xxyyz_zzz, tr_xxyyzz_xx, tr_xxyyzz_xy, tr_xxyyzz_xz, tr_xxyyzz_yy, tr_xxyyzz_yz, tr_xxyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxyy_xx[i] = -2.0 * tr_xxyy_xx[i] * tbe_0 - 2.0 * tr_xxyy_xx[i] * tke_0 + 4.0 * tr_xxyy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xy[i] = -2.0 * tr_xxyy_xy[i] * tbe_0 - 2.0 * tr_xxyy_xy[i] * tke_0 + 4.0 * tr_xxyy_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_xz[i] = -2.0 * tr_xxyy_xz[i] * tbe_0 - 6.0 * tr_xxyy_xz[i] * tke_0 + 4.0 * tr_xxyy_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_x[i] * tbe_0 + 8.0 * tr_xxyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_yy[i] = -2.0 * tr_xxyy_yy[i] * tbe_0 - 2.0 * tr_xxyy_yy[i] * tke_0 + 4.0 * tr_xxyy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_yz[i] = -2.0 * tr_xxyy_yz[i] * tbe_0 - 6.0 * tr_xxyy_yz[i] * tke_0 + 4.0 * tr_xxyy_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyyz_y[i] * tbe_0 + 8.0 * tr_xxyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyy_zz[i] = 2.0 * tr_xxyy_0[i] - 2.0 * tr_xxyy_zz[i] * tbe_0 - 10.0 * tr_xxyy_zz[i] * tke_0 + 4.0 * tr_xxyy_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyyz_z[i] * tbe_0 + 8.0 * tr_xxyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 474-480 components of targeted buffer : GD

    auto tr_0_0_zz_xxyz_xx = pbuffer.data(idx_op_geom_020_gd + 474);

    auto tr_0_0_zz_xxyz_xy = pbuffer.data(idx_op_geom_020_gd + 475);

    auto tr_0_0_zz_xxyz_xz = pbuffer.data(idx_op_geom_020_gd + 476);

    auto tr_0_0_zz_xxyz_yy = pbuffer.data(idx_op_geom_020_gd + 477);

    auto tr_0_0_zz_xxyz_yz = pbuffer.data(idx_op_geom_020_gd + 478);

    auto tr_0_0_zz_xxyz_zz = pbuffer.data(idx_op_geom_020_gd + 479);

    #pragma omp simd aligned(tr_0_0_zz_xxyz_xx, tr_0_0_zz_xxyz_xy, tr_0_0_zz_xxyz_xz, tr_0_0_zz_xxyz_yy, tr_0_0_zz_xxyz_yz, tr_0_0_zz_xxyz_zz, tr_xxy_x, tr_xxy_xxz, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_y, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_z, tr_xxy_zzz, tr_xxyz_0, tr_xxyz_xx, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_x, tr_xxyzz_xxz, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_y, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_z, tr_xxyzz_zzz, tr_xxyzzz_xx, tr_xxyzzz_xy, tr_xxyzzz_xz, tr_xxyzzz_yy, tr_xxyzzz_yz, tr_xxyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxyz_xx[i] = -4.0 * tr_xxy_xxz[i] * tke_0 - 6.0 * tr_xxyz_xx[i] * tbe_0 - 2.0 * tr_xxyz_xx[i] * tke_0 + 4.0 * tr_xxyz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xy[i] = -4.0 * tr_xxy_xyz[i] * tke_0 - 6.0 * tr_xxyz_xy[i] * tbe_0 - 2.0 * tr_xxyz_xy[i] * tke_0 + 4.0 * tr_xxyz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_xz[i] = 2.0 * tr_xxy_x[i] - 4.0 * tr_xxy_xzz[i] * tke_0 - 6.0 * tr_xxyz_xz[i] * tbe_0 - 6.0 * tr_xxyz_xz[i] * tke_0 + 4.0 * tr_xxyz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_x[i] * tbe_0 + 8.0 * tr_xxyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_yy[i] = -4.0 * tr_xxy_yyz[i] * tke_0 - 6.0 * tr_xxyz_yy[i] * tbe_0 - 2.0 * tr_xxyz_yy[i] * tke_0 + 4.0 * tr_xxyz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_yz[i] = 2.0 * tr_xxy_y[i] - 4.0 * tr_xxy_yzz[i] * tke_0 - 6.0 * tr_xxyz_yz[i] * tbe_0 - 6.0 * tr_xxyz_yz[i] * tke_0 + 4.0 * tr_xxyz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxyzz_y[i] * tbe_0 + 8.0 * tr_xxyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxyz_zz[i] = 4.0 * tr_xxy_z[i] - 4.0 * tr_xxy_zzz[i] * tke_0 + 2.0 * tr_xxyz_0[i] - 6.0 * tr_xxyz_zz[i] * tbe_0 - 10.0 * tr_xxyz_zz[i] * tke_0 + 4.0 * tr_xxyz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxyzz_z[i] * tbe_0 + 8.0 * tr_xxyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-486 components of targeted buffer : GD

    auto tr_0_0_zz_xxzz_xx = pbuffer.data(idx_op_geom_020_gd + 480);

    auto tr_0_0_zz_xxzz_xy = pbuffer.data(idx_op_geom_020_gd + 481);

    auto tr_0_0_zz_xxzz_xz = pbuffer.data(idx_op_geom_020_gd + 482);

    auto tr_0_0_zz_xxzz_yy = pbuffer.data(idx_op_geom_020_gd + 483);

    auto tr_0_0_zz_xxzz_yz = pbuffer.data(idx_op_geom_020_gd + 484);

    auto tr_0_0_zz_xxzz_zz = pbuffer.data(idx_op_geom_020_gd + 485);

    #pragma omp simd aligned(tr_0_0_zz_xxzz_xx, tr_0_0_zz_xxzz_xy, tr_0_0_zz_xxzz_xz, tr_0_0_zz_xxzz_yy, tr_0_0_zz_xxzz_yz, tr_0_0_zz_xxzz_zz, tr_xx_xx, tr_xx_xy, tr_xx_xz, tr_xx_yy, tr_xx_yz, tr_xx_zz, tr_xxz_x, tr_xxz_xxz, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_y, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_z, tr_xxz_zzz, tr_xxzz_0, tr_xxzz_xx, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzzz_x, tr_xxzzz_xxz, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_y, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_z, tr_xxzzz_zzz, tr_xxzzzz_xx, tr_xxzzzz_xy, tr_xxzzzz_xz, tr_xxzzzz_yy, tr_xxzzzz_yz, tr_xxzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xxzz_xx[i] = 2.0 * tr_xx_xx[i] - 8.0 * tr_xxz_xxz[i] * tke_0 - 10.0 * tr_xxzz_xx[i] * tbe_0 - 2.0 * tr_xxzz_xx[i] * tke_0 + 4.0 * tr_xxzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xy[i] = 2.0 * tr_xx_xy[i] - 8.0 * tr_xxz_xyz[i] * tke_0 - 10.0 * tr_xxzz_xy[i] * tbe_0 - 2.0 * tr_xxzz_xy[i] * tke_0 + 4.0 * tr_xxzz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_xz[i] = 2.0 * tr_xx_xz[i] + 4.0 * tr_xxz_x[i] - 8.0 * tr_xxz_xzz[i] * tke_0 - 10.0 * tr_xxzz_xz[i] * tbe_0 - 6.0 * tr_xxzz_xz[i] * tke_0 + 4.0 * tr_xxzz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_x[i] * tbe_0 + 8.0 * tr_xxzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_yy[i] = 2.0 * tr_xx_yy[i] - 8.0 * tr_xxz_yyz[i] * tke_0 - 10.0 * tr_xxzz_yy[i] * tbe_0 - 2.0 * tr_xxzz_yy[i] * tke_0 + 4.0 * tr_xxzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xxzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_yz[i] = 2.0 * tr_xx_yz[i] + 4.0 * tr_xxz_y[i] - 8.0 * tr_xxz_yzz[i] * tke_0 - 10.0 * tr_xxzz_yz[i] * tbe_0 - 6.0 * tr_xxzz_yz[i] * tke_0 + 4.0 * tr_xxzz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xxzzz_y[i] * tbe_0 + 8.0 * tr_xxzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xxzz_zz[i] = 2.0 * tr_xx_zz[i] + 8.0 * tr_xxz_z[i] - 8.0 * tr_xxz_zzz[i] * tke_0 + 2.0 * tr_xxzz_0[i] - 10.0 * tr_xxzz_zz[i] * tbe_0 - 10.0 * tr_xxzz_zz[i] * tke_0 + 4.0 * tr_xxzz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xxzzz_z[i] * tbe_0 + 8.0 * tr_xxzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 486-492 components of targeted buffer : GD

    auto tr_0_0_zz_xyyy_xx = pbuffer.data(idx_op_geom_020_gd + 486);

    auto tr_0_0_zz_xyyy_xy = pbuffer.data(idx_op_geom_020_gd + 487);

    auto tr_0_0_zz_xyyy_xz = pbuffer.data(idx_op_geom_020_gd + 488);

    auto tr_0_0_zz_xyyy_yy = pbuffer.data(idx_op_geom_020_gd + 489);

    auto tr_0_0_zz_xyyy_yz = pbuffer.data(idx_op_geom_020_gd + 490);

    auto tr_0_0_zz_xyyy_zz = pbuffer.data(idx_op_geom_020_gd + 491);

    #pragma omp simd aligned(tr_0_0_zz_xyyy_xx, tr_0_0_zz_xyyy_xy, tr_0_0_zz_xyyy_xz, tr_0_0_zz_xyyy_yy, tr_0_0_zz_xyyy_yz, tr_0_0_zz_xyyy_zz, tr_xyyy_0, tr_xyyy_xx, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyz_x, tr_xyyyz_xxz, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_y, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_z, tr_xyyyz_zzz, tr_xyyyzz_xx, tr_xyyyzz_xy, tr_xyyyzz_xz, tr_xyyyzz_yy, tr_xyyyzz_yz, tr_xyyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyyy_xx[i] = -2.0 * tr_xyyy_xx[i] * tbe_0 - 2.0 * tr_xyyy_xx[i] * tke_0 + 4.0 * tr_xyyy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xy[i] = -2.0 * tr_xyyy_xy[i] * tbe_0 - 2.0 * tr_xyyy_xy[i] * tke_0 + 4.0 * tr_xyyy_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_xz[i] = -2.0 * tr_xyyy_xz[i] * tbe_0 - 6.0 * tr_xyyy_xz[i] * tke_0 + 4.0 * tr_xyyy_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_x[i] * tbe_0 + 8.0 * tr_xyyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_yy[i] = -2.0 * tr_xyyy_yy[i] * tbe_0 - 2.0 * tr_xyyy_yy[i] * tke_0 + 4.0 * tr_xyyy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_yz[i] = -2.0 * tr_xyyy_yz[i] * tbe_0 - 6.0 * tr_xyyy_yz[i] * tke_0 + 4.0 * tr_xyyy_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyyz_y[i] * tbe_0 + 8.0 * tr_xyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyy_zz[i] = 2.0 * tr_xyyy_0[i] - 2.0 * tr_xyyy_zz[i] * tbe_0 - 10.0 * tr_xyyy_zz[i] * tke_0 + 4.0 * tr_xyyy_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyyz_z[i] * tbe_0 + 8.0 * tr_xyyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 492-498 components of targeted buffer : GD

    auto tr_0_0_zz_xyyz_xx = pbuffer.data(idx_op_geom_020_gd + 492);

    auto tr_0_0_zz_xyyz_xy = pbuffer.data(idx_op_geom_020_gd + 493);

    auto tr_0_0_zz_xyyz_xz = pbuffer.data(idx_op_geom_020_gd + 494);

    auto tr_0_0_zz_xyyz_yy = pbuffer.data(idx_op_geom_020_gd + 495);

    auto tr_0_0_zz_xyyz_yz = pbuffer.data(idx_op_geom_020_gd + 496);

    auto tr_0_0_zz_xyyz_zz = pbuffer.data(idx_op_geom_020_gd + 497);

    #pragma omp simd aligned(tr_0_0_zz_xyyz_xx, tr_0_0_zz_xyyz_xy, tr_0_0_zz_xyyz_xz, tr_0_0_zz_xyyz_yy, tr_0_0_zz_xyyz_yz, tr_0_0_zz_xyyz_zz, tr_xyy_x, tr_xyy_xxz, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_y, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_z, tr_xyy_zzz, tr_xyyz_0, tr_xyyz_xx, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_x, tr_xyyzz_xxz, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_y, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_z, tr_xyyzz_zzz, tr_xyyzzz_xx, tr_xyyzzz_xy, tr_xyyzzz_xz, tr_xyyzzz_yy, tr_xyyzzz_yz, tr_xyyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyyz_xx[i] = -4.0 * tr_xyy_xxz[i] * tke_0 - 6.0 * tr_xyyz_xx[i] * tbe_0 - 2.0 * tr_xyyz_xx[i] * tke_0 + 4.0 * tr_xyyz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xy[i] = -4.0 * tr_xyy_xyz[i] * tke_0 - 6.0 * tr_xyyz_xy[i] * tbe_0 - 2.0 * tr_xyyz_xy[i] * tke_0 + 4.0 * tr_xyyz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_xz[i] = 2.0 * tr_xyy_x[i] - 4.0 * tr_xyy_xzz[i] * tke_0 - 6.0 * tr_xyyz_xz[i] * tbe_0 - 6.0 * tr_xyyz_xz[i] * tke_0 + 4.0 * tr_xyyz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_x[i] * tbe_0 + 8.0 * tr_xyyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_yy[i] = -4.0 * tr_xyy_yyz[i] * tke_0 - 6.0 * tr_xyyz_yy[i] * tbe_0 - 2.0 * tr_xyyz_yy[i] * tke_0 + 4.0 * tr_xyyz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_yz[i] = 2.0 * tr_xyy_y[i] - 4.0 * tr_xyy_yzz[i] * tke_0 - 6.0 * tr_xyyz_yz[i] * tbe_0 - 6.0 * tr_xyyz_yz[i] * tke_0 + 4.0 * tr_xyyz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyyzz_y[i] * tbe_0 + 8.0 * tr_xyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyyz_zz[i] = 4.0 * tr_xyy_z[i] - 4.0 * tr_xyy_zzz[i] * tke_0 + 2.0 * tr_xyyz_0[i] - 6.0 * tr_xyyz_zz[i] * tbe_0 - 10.0 * tr_xyyz_zz[i] * tke_0 + 4.0 * tr_xyyz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyyzz_z[i] * tbe_0 + 8.0 * tr_xyyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 498-504 components of targeted buffer : GD

    auto tr_0_0_zz_xyzz_xx = pbuffer.data(idx_op_geom_020_gd + 498);

    auto tr_0_0_zz_xyzz_xy = pbuffer.data(idx_op_geom_020_gd + 499);

    auto tr_0_0_zz_xyzz_xz = pbuffer.data(idx_op_geom_020_gd + 500);

    auto tr_0_0_zz_xyzz_yy = pbuffer.data(idx_op_geom_020_gd + 501);

    auto tr_0_0_zz_xyzz_yz = pbuffer.data(idx_op_geom_020_gd + 502);

    auto tr_0_0_zz_xyzz_zz = pbuffer.data(idx_op_geom_020_gd + 503);

    #pragma omp simd aligned(tr_0_0_zz_xyzz_xx, tr_0_0_zz_xyzz_xy, tr_0_0_zz_xyzz_xz, tr_0_0_zz_xyzz_yy, tr_0_0_zz_xyzz_yz, tr_0_0_zz_xyzz_zz, tr_xy_xx, tr_xy_xy, tr_xy_xz, tr_xy_yy, tr_xy_yz, tr_xy_zz, tr_xyz_x, tr_xyz_xxz, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_y, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_z, tr_xyz_zzz, tr_xyzz_0, tr_xyzz_xx, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_x, tr_xyzzz_xxz, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_y, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_z, tr_xyzzz_zzz, tr_xyzzzz_xx, tr_xyzzzz_xy, tr_xyzzzz_xz, tr_xyzzzz_yy, tr_xyzzzz_yz, tr_xyzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xyzz_xx[i] = 2.0 * tr_xy_xx[i] - 8.0 * tr_xyz_xxz[i] * tke_0 - 10.0 * tr_xyzz_xx[i] * tbe_0 - 2.0 * tr_xyzz_xx[i] * tke_0 + 4.0 * tr_xyzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xy[i] = 2.0 * tr_xy_xy[i] - 8.0 * tr_xyz_xyz[i] * tke_0 - 10.0 * tr_xyzz_xy[i] * tbe_0 - 2.0 * tr_xyzz_xy[i] * tke_0 + 4.0 * tr_xyzz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_xz[i] = 2.0 * tr_xy_xz[i] + 4.0 * tr_xyz_x[i] - 8.0 * tr_xyz_xzz[i] * tke_0 - 10.0 * tr_xyzz_xz[i] * tbe_0 - 6.0 * tr_xyzz_xz[i] * tke_0 + 4.0 * tr_xyzz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_x[i] * tbe_0 + 8.0 * tr_xyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_yy[i] = 2.0 * tr_xy_yy[i] - 8.0 * tr_xyz_yyz[i] * tke_0 - 10.0 * tr_xyzz_yy[i] * tbe_0 - 2.0 * tr_xyzz_yy[i] * tke_0 + 4.0 * tr_xyzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_yz[i] = 2.0 * tr_xy_yz[i] + 4.0 * tr_xyz_y[i] - 8.0 * tr_xyz_yzz[i] * tke_0 - 10.0 * tr_xyzz_yz[i] * tbe_0 - 6.0 * tr_xyzz_yz[i] * tke_0 + 4.0 * tr_xyzz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xyzzz_y[i] * tbe_0 + 8.0 * tr_xyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xyzz_zz[i] = 2.0 * tr_xy_zz[i] + 8.0 * tr_xyz_z[i] - 8.0 * tr_xyz_zzz[i] * tke_0 + 2.0 * tr_xyzz_0[i] - 10.0 * tr_xyzz_zz[i] * tbe_0 - 10.0 * tr_xyzz_zz[i] * tke_0 + 4.0 * tr_xyzz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xyzzz_z[i] * tbe_0 + 8.0 * tr_xyzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 504-510 components of targeted buffer : GD

    auto tr_0_0_zz_xzzz_xx = pbuffer.data(idx_op_geom_020_gd + 504);

    auto tr_0_0_zz_xzzz_xy = pbuffer.data(idx_op_geom_020_gd + 505);

    auto tr_0_0_zz_xzzz_xz = pbuffer.data(idx_op_geom_020_gd + 506);

    auto tr_0_0_zz_xzzz_yy = pbuffer.data(idx_op_geom_020_gd + 507);

    auto tr_0_0_zz_xzzz_yz = pbuffer.data(idx_op_geom_020_gd + 508);

    auto tr_0_0_zz_xzzz_zz = pbuffer.data(idx_op_geom_020_gd + 509);

    #pragma omp simd aligned(tr_0_0_zz_xzzz_xx, tr_0_0_zz_xzzz_xy, tr_0_0_zz_xzzz_xz, tr_0_0_zz_xzzz_yy, tr_0_0_zz_xzzz_yz, tr_0_0_zz_xzzz_zz, tr_xz_xx, tr_xz_xy, tr_xz_xz, tr_xz_yy, tr_xz_yz, tr_xz_zz, tr_xzz_x, tr_xzz_xxz, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_y, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_z, tr_xzz_zzz, tr_xzzz_0, tr_xzzz_xx, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzzz_x, tr_xzzzz_xxz, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_y, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_z, tr_xzzzz_zzz, tr_xzzzzz_xx, tr_xzzzzz_xy, tr_xzzzzz_xz, tr_xzzzzz_yy, tr_xzzzzz_yz, tr_xzzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_xzzz_xx[i] = 6.0 * tr_xz_xx[i] - 12.0 * tr_xzz_xxz[i] * tke_0 - 14.0 * tr_xzzz_xx[i] * tbe_0 - 2.0 * tr_xzzz_xx[i] * tke_0 + 4.0 * tr_xzzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xy[i] = 6.0 * tr_xz_xy[i] - 12.0 * tr_xzz_xyz[i] * tke_0 - 14.0 * tr_xzzz_xy[i] * tbe_0 - 2.0 * tr_xzzz_xy[i] * tke_0 + 4.0 * tr_xzzz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_xz[i] = 6.0 * tr_xz_xz[i] + 6.0 * tr_xzz_x[i] - 12.0 * tr_xzz_xzz[i] * tke_0 - 14.0 * tr_xzzz_xz[i] * tbe_0 - 6.0 * tr_xzzz_xz[i] * tke_0 + 4.0 * tr_xzzz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_x[i] * tbe_0 + 8.0 * tr_xzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_yy[i] = 6.0 * tr_xz_yy[i] - 12.0 * tr_xzz_yyz[i] * tke_0 - 14.0 * tr_xzzz_yy[i] * tbe_0 - 2.0 * tr_xzzz_yy[i] * tke_0 + 4.0 * tr_xzzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_xzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_yz[i] = 6.0 * tr_xz_yz[i] + 6.0 * tr_xzz_y[i] - 12.0 * tr_xzz_yzz[i] * tke_0 - 14.0 * tr_xzzz_yz[i] * tbe_0 - 6.0 * tr_xzzz_yz[i] * tke_0 + 4.0 * tr_xzzz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_xzzzz_y[i] * tbe_0 + 8.0 * tr_xzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_xzzz_zz[i] = 6.0 * tr_xz_zz[i] + 12.0 * tr_xzz_z[i] - 12.0 * tr_xzz_zzz[i] * tke_0 + 2.0 * tr_xzzz_0[i] - 14.0 * tr_xzzz_zz[i] * tbe_0 - 10.0 * tr_xzzz_zz[i] * tke_0 + 4.0 * tr_xzzz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_xzzzz_z[i] * tbe_0 + 8.0 * tr_xzzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-516 components of targeted buffer : GD

    auto tr_0_0_zz_yyyy_xx = pbuffer.data(idx_op_geom_020_gd + 510);

    auto tr_0_0_zz_yyyy_xy = pbuffer.data(idx_op_geom_020_gd + 511);

    auto tr_0_0_zz_yyyy_xz = pbuffer.data(idx_op_geom_020_gd + 512);

    auto tr_0_0_zz_yyyy_yy = pbuffer.data(idx_op_geom_020_gd + 513);

    auto tr_0_0_zz_yyyy_yz = pbuffer.data(idx_op_geom_020_gd + 514);

    auto tr_0_0_zz_yyyy_zz = pbuffer.data(idx_op_geom_020_gd + 515);

    #pragma omp simd aligned(tr_0_0_zz_yyyy_xx, tr_0_0_zz_yyyy_xy, tr_0_0_zz_yyyy_xz, tr_0_0_zz_yyyy_yy, tr_0_0_zz_yyyy_yz, tr_0_0_zz_yyyy_zz, tr_yyyy_0, tr_yyyy_xx, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyyz_x, tr_yyyyz_xxz, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_y, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_z, tr_yyyyz_zzz, tr_yyyyzz_xx, tr_yyyyzz_xy, tr_yyyyzz_xz, tr_yyyyzz_yy, tr_yyyyzz_yz, tr_yyyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyyy_xx[i] = -2.0 * tr_yyyy_xx[i] * tbe_0 - 2.0 * tr_yyyy_xx[i] * tke_0 + 4.0 * tr_yyyy_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xy[i] = -2.0 * tr_yyyy_xy[i] * tbe_0 - 2.0 * tr_yyyy_xy[i] * tke_0 + 4.0 * tr_yyyy_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_xz[i] = -2.0 * tr_yyyy_xz[i] * tbe_0 - 6.0 * tr_yyyy_xz[i] * tke_0 + 4.0 * tr_yyyy_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_x[i] * tbe_0 + 8.0 * tr_yyyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_yy[i] = -2.0 * tr_yyyy_yy[i] * tbe_0 - 2.0 * tr_yyyy_yy[i] * tke_0 + 4.0 * tr_yyyy_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_yz[i] = -2.0 * tr_yyyy_yz[i] * tbe_0 - 6.0 * tr_yyyy_yz[i] * tke_0 + 4.0 * tr_yyyy_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyyz_y[i] * tbe_0 + 8.0 * tr_yyyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyy_zz[i] = 2.0 * tr_yyyy_0[i] - 2.0 * tr_yyyy_zz[i] * tbe_0 - 10.0 * tr_yyyy_zz[i] * tke_0 + 4.0 * tr_yyyy_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyyz_z[i] * tbe_0 + 8.0 * tr_yyyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 516-522 components of targeted buffer : GD

    auto tr_0_0_zz_yyyz_xx = pbuffer.data(idx_op_geom_020_gd + 516);

    auto tr_0_0_zz_yyyz_xy = pbuffer.data(idx_op_geom_020_gd + 517);

    auto tr_0_0_zz_yyyz_xz = pbuffer.data(idx_op_geom_020_gd + 518);

    auto tr_0_0_zz_yyyz_yy = pbuffer.data(idx_op_geom_020_gd + 519);

    auto tr_0_0_zz_yyyz_yz = pbuffer.data(idx_op_geom_020_gd + 520);

    auto tr_0_0_zz_yyyz_zz = pbuffer.data(idx_op_geom_020_gd + 521);

    #pragma omp simd aligned(tr_0_0_zz_yyyz_xx, tr_0_0_zz_yyyz_xy, tr_0_0_zz_yyyz_xz, tr_0_0_zz_yyyz_yy, tr_0_0_zz_yyyz_yz, tr_0_0_zz_yyyz_zz, tr_yyy_x, tr_yyy_xxz, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_y, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_z, tr_yyy_zzz, tr_yyyz_0, tr_yyyz_xx, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyzz_x, tr_yyyzz_xxz, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_y, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_z, tr_yyyzz_zzz, tr_yyyzzz_xx, tr_yyyzzz_xy, tr_yyyzzz_xz, tr_yyyzzz_yy, tr_yyyzzz_yz, tr_yyyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyyz_xx[i] = -4.0 * tr_yyy_xxz[i] * tke_0 - 6.0 * tr_yyyz_xx[i] * tbe_0 - 2.0 * tr_yyyz_xx[i] * tke_0 + 4.0 * tr_yyyz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xy[i] = -4.0 * tr_yyy_xyz[i] * tke_0 - 6.0 * tr_yyyz_xy[i] * tbe_0 - 2.0 * tr_yyyz_xy[i] * tke_0 + 4.0 * tr_yyyz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_xz[i] = 2.0 * tr_yyy_x[i] - 4.0 * tr_yyy_xzz[i] * tke_0 - 6.0 * tr_yyyz_xz[i] * tbe_0 - 6.0 * tr_yyyz_xz[i] * tke_0 + 4.0 * tr_yyyz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_x[i] * tbe_0 + 8.0 * tr_yyyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_yy[i] = -4.0 * tr_yyy_yyz[i] * tke_0 - 6.0 * tr_yyyz_yy[i] * tbe_0 - 2.0 * tr_yyyz_yy[i] * tke_0 + 4.0 * tr_yyyz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_yz[i] = 2.0 * tr_yyy_y[i] - 4.0 * tr_yyy_yzz[i] * tke_0 - 6.0 * tr_yyyz_yz[i] * tbe_0 - 6.0 * tr_yyyz_yz[i] * tke_0 + 4.0 * tr_yyyz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyyzz_y[i] * tbe_0 + 8.0 * tr_yyyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyyz_zz[i] = 4.0 * tr_yyy_z[i] - 4.0 * tr_yyy_zzz[i] * tke_0 + 2.0 * tr_yyyz_0[i] - 6.0 * tr_yyyz_zz[i] * tbe_0 - 10.0 * tr_yyyz_zz[i] * tke_0 + 4.0 * tr_yyyz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyyzz_z[i] * tbe_0 + 8.0 * tr_yyyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 522-528 components of targeted buffer : GD

    auto tr_0_0_zz_yyzz_xx = pbuffer.data(idx_op_geom_020_gd + 522);

    auto tr_0_0_zz_yyzz_xy = pbuffer.data(idx_op_geom_020_gd + 523);

    auto tr_0_0_zz_yyzz_xz = pbuffer.data(idx_op_geom_020_gd + 524);

    auto tr_0_0_zz_yyzz_yy = pbuffer.data(idx_op_geom_020_gd + 525);

    auto tr_0_0_zz_yyzz_yz = pbuffer.data(idx_op_geom_020_gd + 526);

    auto tr_0_0_zz_yyzz_zz = pbuffer.data(idx_op_geom_020_gd + 527);

    #pragma omp simd aligned(tr_0_0_zz_yyzz_xx, tr_0_0_zz_yyzz_xy, tr_0_0_zz_yyzz_xz, tr_0_0_zz_yyzz_yy, tr_0_0_zz_yyzz_yz, tr_0_0_zz_yyzz_zz, tr_yy_xx, tr_yy_xy, tr_yy_xz, tr_yy_yy, tr_yy_yz, tr_yy_zz, tr_yyz_x, tr_yyz_xxz, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_y, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_z, tr_yyz_zzz, tr_yyzz_0, tr_yyzz_xx, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzzz_x, tr_yyzzz_xxz, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_y, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_z, tr_yyzzz_zzz, tr_yyzzzz_xx, tr_yyzzzz_xy, tr_yyzzzz_xz, tr_yyzzzz_yy, tr_yyzzzz_yz, tr_yyzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yyzz_xx[i] = 2.0 * tr_yy_xx[i] - 8.0 * tr_yyz_xxz[i] * tke_0 - 10.0 * tr_yyzz_xx[i] * tbe_0 - 2.0 * tr_yyzz_xx[i] * tke_0 + 4.0 * tr_yyzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xy[i] = 2.0 * tr_yy_xy[i] - 8.0 * tr_yyz_xyz[i] * tke_0 - 10.0 * tr_yyzz_xy[i] * tbe_0 - 2.0 * tr_yyzz_xy[i] * tke_0 + 4.0 * tr_yyzz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_xz[i] = 2.0 * tr_yy_xz[i] + 4.0 * tr_yyz_x[i] - 8.0 * tr_yyz_xzz[i] * tke_0 - 10.0 * tr_yyzz_xz[i] * tbe_0 - 6.0 * tr_yyzz_xz[i] * tke_0 + 4.0 * tr_yyzz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_x[i] * tbe_0 + 8.0 * tr_yyzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_yy[i] = 2.0 * tr_yy_yy[i] - 8.0 * tr_yyz_yyz[i] * tke_0 - 10.0 * tr_yyzz_yy[i] * tbe_0 - 2.0 * tr_yyzz_yy[i] * tke_0 + 4.0 * tr_yyzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yyzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_yz[i] = 2.0 * tr_yy_yz[i] + 4.0 * tr_yyz_y[i] - 8.0 * tr_yyz_yzz[i] * tke_0 - 10.0 * tr_yyzz_yz[i] * tbe_0 - 6.0 * tr_yyzz_yz[i] * tke_0 + 4.0 * tr_yyzz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_yyzzz_y[i] * tbe_0 + 8.0 * tr_yyzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yyzz_zz[i] = 2.0 * tr_yy_zz[i] + 8.0 * tr_yyz_z[i] - 8.0 * tr_yyz_zzz[i] * tke_0 + 2.0 * tr_yyzz_0[i] - 10.0 * tr_yyzz_zz[i] * tbe_0 - 10.0 * tr_yyzz_zz[i] * tke_0 + 4.0 * tr_yyzz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_yyzzz_z[i] * tbe_0 + 8.0 * tr_yyzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 528-534 components of targeted buffer : GD

    auto tr_0_0_zz_yzzz_xx = pbuffer.data(idx_op_geom_020_gd + 528);

    auto tr_0_0_zz_yzzz_xy = pbuffer.data(idx_op_geom_020_gd + 529);

    auto tr_0_0_zz_yzzz_xz = pbuffer.data(idx_op_geom_020_gd + 530);

    auto tr_0_0_zz_yzzz_yy = pbuffer.data(idx_op_geom_020_gd + 531);

    auto tr_0_0_zz_yzzz_yz = pbuffer.data(idx_op_geom_020_gd + 532);

    auto tr_0_0_zz_yzzz_zz = pbuffer.data(idx_op_geom_020_gd + 533);

    #pragma omp simd aligned(tr_0_0_zz_yzzz_xx, tr_0_0_zz_yzzz_xy, tr_0_0_zz_yzzz_xz, tr_0_0_zz_yzzz_yy, tr_0_0_zz_yzzz_yz, tr_0_0_zz_yzzz_zz, tr_yz_xx, tr_yz_xy, tr_yz_xz, tr_yz_yy, tr_yz_yz, tr_yz_zz, tr_yzz_x, tr_yzz_xxz, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_y, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_z, tr_yzz_zzz, tr_yzzz_0, tr_yzzz_xx, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzzz_x, tr_yzzzz_xxz, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_y, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_z, tr_yzzzz_zzz, tr_yzzzzz_xx, tr_yzzzzz_xy, tr_yzzzzz_xz, tr_yzzzzz_yy, tr_yzzzzz_yz, tr_yzzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_yzzz_xx[i] = 6.0 * tr_yz_xx[i] - 12.0 * tr_yzz_xxz[i] * tke_0 - 14.0 * tr_yzzz_xx[i] * tbe_0 - 2.0 * tr_yzzz_xx[i] * tke_0 + 4.0 * tr_yzzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xy[i] = 6.0 * tr_yz_xy[i] - 12.0 * tr_yzz_xyz[i] * tke_0 - 14.0 * tr_yzzz_xy[i] * tbe_0 - 2.0 * tr_yzzz_xy[i] * tke_0 + 4.0 * tr_yzzz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_xz[i] = 6.0 * tr_yz_xz[i] + 6.0 * tr_yzz_x[i] - 12.0 * tr_yzz_xzz[i] * tke_0 - 14.0 * tr_yzzz_xz[i] * tbe_0 - 6.0 * tr_yzzz_xz[i] * tke_0 + 4.0 * tr_yzzz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_x[i] * tbe_0 + 8.0 * tr_yzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_yy[i] = 6.0 * tr_yz_yy[i] - 12.0 * tr_yzz_yyz[i] * tke_0 - 14.0 * tr_yzzz_yy[i] * tbe_0 - 2.0 * tr_yzzz_yy[i] * tke_0 + 4.0 * tr_yzzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_yzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_yz[i] = 6.0 * tr_yz_yz[i] + 6.0 * tr_yzz_y[i] - 12.0 * tr_yzz_yzz[i] * tke_0 - 14.0 * tr_yzzz_yz[i] * tbe_0 - 6.0 * tr_yzzz_yz[i] * tke_0 + 4.0 * tr_yzzz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_yzzzz_y[i] * tbe_0 + 8.0 * tr_yzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_yzzz_zz[i] = 6.0 * tr_yz_zz[i] + 12.0 * tr_yzz_z[i] - 12.0 * tr_yzz_zzz[i] * tke_0 + 2.0 * tr_yzzz_0[i] - 14.0 * tr_yzzz_zz[i] * tbe_0 - 10.0 * tr_yzzz_zz[i] * tke_0 + 4.0 * tr_yzzz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_yzzzz_z[i] * tbe_0 + 8.0 * tr_yzzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 534-540 components of targeted buffer : GD

    auto tr_0_0_zz_zzzz_xx = pbuffer.data(idx_op_geom_020_gd + 534);

    auto tr_0_0_zz_zzzz_xy = pbuffer.data(idx_op_geom_020_gd + 535);

    auto tr_0_0_zz_zzzz_xz = pbuffer.data(idx_op_geom_020_gd + 536);

    auto tr_0_0_zz_zzzz_yy = pbuffer.data(idx_op_geom_020_gd + 537);

    auto tr_0_0_zz_zzzz_yz = pbuffer.data(idx_op_geom_020_gd + 538);

    auto tr_0_0_zz_zzzz_zz = pbuffer.data(idx_op_geom_020_gd + 539);

    #pragma omp simd aligned(tr_0_0_zz_zzzz_xx, tr_0_0_zz_zzzz_xy, tr_0_0_zz_zzzz_xz, tr_0_0_zz_zzzz_yy, tr_0_0_zz_zzzz_yz, tr_0_0_zz_zzzz_zz, tr_zz_xx, tr_zz_xy, tr_zz_xz, tr_zz_yy, tr_zz_yz, tr_zz_zz, tr_zzz_x, tr_zzz_xxz, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_y, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_z, tr_zzz_zzz, tr_zzzz_0, tr_zzzz_xx, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzz_zzzz, tr_zzzzz_x, tr_zzzzz_xxz, tr_zzzzz_xyz, tr_zzzzz_xzz, tr_zzzzz_y, tr_zzzzz_yyz, tr_zzzzz_yzz, tr_zzzzz_z, tr_zzzzz_zzz, tr_zzzzzz_xx, tr_zzzzzz_xy, tr_zzzzzz_xz, tr_zzzzzz_yy, tr_zzzzzz_yz, tr_zzzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_zz_zzzz_xx[i] = 12.0 * tr_zz_xx[i] - 16.0 * tr_zzz_xxz[i] * tke_0 - 18.0 * tr_zzzz_xx[i] * tbe_0 - 2.0 * tr_zzzz_xx[i] * tke_0 + 4.0 * tr_zzzz_xxzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xx[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xy[i] = 12.0 * tr_zz_xy[i] - 16.0 * tr_zzz_xyz[i] * tke_0 - 18.0 * tr_zzzz_xy[i] * tbe_0 - 2.0 * tr_zzzz_xy[i] * tke_0 + 4.0 * tr_zzzz_xyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_xz[i] = 12.0 * tr_zz_xz[i] + 8.0 * tr_zzz_x[i] - 16.0 * tr_zzz_xzz[i] * tke_0 - 18.0 * tr_zzzz_xz[i] * tbe_0 - 6.0 * tr_zzzz_xz[i] * tke_0 + 4.0 * tr_zzzz_xzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_x[i] * tbe_0 + 8.0 * tr_zzzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_yy[i] = 12.0 * tr_zz_yy[i] - 16.0 * tr_zzz_yyz[i] * tke_0 - 18.0 * tr_zzzz_yy[i] * tbe_0 - 2.0 * tr_zzzz_yy[i] * tke_0 + 4.0 * tr_zzzz_yyzz[i] * tke_0 * tke_0 + 8.0 * tr_zzzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yy[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_yz[i] = 12.0 * tr_zz_yz[i] + 8.0 * tr_zzz_y[i] - 16.0 * tr_zzz_yzz[i] * tke_0 - 18.0 * tr_zzzz_yz[i] * tbe_0 - 6.0 * tr_zzzz_yz[i] * tke_0 + 4.0 * tr_zzzz_yzzz[i] * tke_0 * tke_0 - 4.0 * tr_zzzzz_y[i] * tbe_0 + 8.0 * tr_zzzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yz[i] * tbe_0 * tbe_0;

        tr_0_0_zz_zzzz_zz[i] = 12.0 * tr_zz_zz[i] + 16.0 * tr_zzz_z[i] - 16.0 * tr_zzz_zzz[i] * tke_0 + 2.0 * tr_zzzz_0[i] - 18.0 * tr_zzzz_zz[i] * tbe_0 - 10.0 * tr_zzzz_zz[i] * tke_0 + 4.0 * tr_zzzz_zzzz[i] * tke_0 * tke_0 - 8.0 * tr_zzzzz_z[i] * tbe_0 + 8.0 * tr_zzzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_zz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

