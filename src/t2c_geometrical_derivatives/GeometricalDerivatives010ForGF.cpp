#include "GeometricalDerivatives010ForGF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_gf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gf,
                         const int idx_op_ff,
                         const int idx_op_gd,
                         const int idx_op_gg,
                         const int idx_op_hf,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up 0-10 components of targeted buffer : GF

    auto tr_0_0_x_xxxx_xxx = pbuffer.data(idx_op_geom_010_gf);

    auto tr_0_0_x_xxxx_xxy = pbuffer.data(idx_op_geom_010_gf + 1);

    auto tr_0_0_x_xxxx_xxz = pbuffer.data(idx_op_geom_010_gf + 2);

    auto tr_0_0_x_xxxx_xyy = pbuffer.data(idx_op_geom_010_gf + 3);

    auto tr_0_0_x_xxxx_xyz = pbuffer.data(idx_op_geom_010_gf + 4);

    auto tr_0_0_x_xxxx_xzz = pbuffer.data(idx_op_geom_010_gf + 5);

    auto tr_0_0_x_xxxx_yyy = pbuffer.data(idx_op_geom_010_gf + 6);

    auto tr_0_0_x_xxxx_yyz = pbuffer.data(idx_op_geom_010_gf + 7);

    auto tr_0_0_x_xxxx_yzz = pbuffer.data(idx_op_geom_010_gf + 8);

    auto tr_0_0_x_xxxx_zzz = pbuffer.data(idx_op_geom_010_gf + 9);

    #pragma omp simd aligned(tr_0_0_x_xxxx_xxx, tr_0_0_x_xxxx_xxy, tr_0_0_x_xxxx_xxz, tr_0_0_x_xxxx_xyy, tr_0_0_x_xxxx_xyz, tr_0_0_x_xxxx_xzz, tr_0_0_x_xxxx_yyy, tr_0_0_x_xxxx_yyz, tr_0_0_x_xxxx_yzz, tr_0_0_x_xxxx_zzz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxx_xx, tr_xxxx_xxxx, tr_xxxx_xxxy, tr_xxxx_xxxz, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yz, tr_xxxx_zz, tr_xxxxx_xxx, tr_xxxxx_xxy, tr_xxxxx_xxz, tr_xxxxx_xyy, tr_xxxxx_xyz, tr_xxxxx_xzz, tr_xxxxx_yyy, tr_xxxxx_yyz, tr_xxxxx_yzz, tr_xxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxx_xxx[i] = 2.0 * tr_xxxxx_xxx[i] * tbe_0 + 2.0 * tr_xxxx_xxxx[i] * tke_0 - 4.0 * tr_xxx_xxx[i] - 3.0 * tr_xxxx_xx[i];

        tr_0_0_x_xxxx_xxy[i] = 2.0 * tr_xxxxx_xxy[i] * tbe_0 + 2.0 * tr_xxxx_xxxy[i] * tke_0 - 4.0 * tr_xxx_xxy[i] - 2.0 * tr_xxxx_xy[i];

        tr_0_0_x_xxxx_xxz[i] = 2.0 * tr_xxxxx_xxz[i] * tbe_0 + 2.0 * tr_xxxx_xxxz[i] * tke_0 - 4.0 * tr_xxx_xxz[i] - 2.0 * tr_xxxx_xz[i];

        tr_0_0_x_xxxx_xyy[i] = 2.0 * tr_xxxxx_xyy[i] * tbe_0 + 2.0 * tr_xxxx_xxyy[i] * tke_0 - 4.0 * tr_xxx_xyy[i] - tr_xxxx_yy[i];

        tr_0_0_x_xxxx_xyz[i] = 2.0 * tr_xxxxx_xyz[i] * tbe_0 + 2.0 * tr_xxxx_xxyz[i] * tke_0 - 4.0 * tr_xxx_xyz[i] - tr_xxxx_yz[i];

        tr_0_0_x_xxxx_xzz[i] = 2.0 * tr_xxxxx_xzz[i] * tbe_0 + 2.0 * tr_xxxx_xxzz[i] * tke_0 - 4.0 * tr_xxx_xzz[i] - tr_xxxx_zz[i];

        tr_0_0_x_xxxx_yyy[i] = 2.0 * tr_xxxxx_yyy[i] * tbe_0 + 2.0 * tr_xxxx_xyyy[i] * tke_0 - 4.0 * tr_xxx_yyy[i];

        tr_0_0_x_xxxx_yyz[i] = 2.0 * tr_xxxxx_yyz[i] * tbe_0 + 2.0 * tr_xxxx_xyyz[i] * tke_0 - 4.0 * tr_xxx_yyz[i];

        tr_0_0_x_xxxx_yzz[i] = 2.0 * tr_xxxxx_yzz[i] * tbe_0 + 2.0 * tr_xxxx_xyzz[i] * tke_0 - 4.0 * tr_xxx_yzz[i];

        tr_0_0_x_xxxx_zzz[i] = 2.0 * tr_xxxxx_zzz[i] * tbe_0 + 2.0 * tr_xxxx_xzzz[i] * tke_0 - 4.0 * tr_xxx_zzz[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto tr_0_0_x_xxxy_xxx = pbuffer.data(idx_op_geom_010_gf + 10);

    auto tr_0_0_x_xxxy_xxy = pbuffer.data(idx_op_geom_010_gf + 11);

    auto tr_0_0_x_xxxy_xxz = pbuffer.data(idx_op_geom_010_gf + 12);

    auto tr_0_0_x_xxxy_xyy = pbuffer.data(idx_op_geom_010_gf + 13);

    auto tr_0_0_x_xxxy_xyz = pbuffer.data(idx_op_geom_010_gf + 14);

    auto tr_0_0_x_xxxy_xzz = pbuffer.data(idx_op_geom_010_gf + 15);

    auto tr_0_0_x_xxxy_yyy = pbuffer.data(idx_op_geom_010_gf + 16);

    auto tr_0_0_x_xxxy_yyz = pbuffer.data(idx_op_geom_010_gf + 17);

    auto tr_0_0_x_xxxy_yzz = pbuffer.data(idx_op_geom_010_gf + 18);

    auto tr_0_0_x_xxxy_zzz = pbuffer.data(idx_op_geom_010_gf + 19);

    #pragma omp simd aligned(tr_0_0_x_xxxy_xxx, tr_0_0_x_xxxy_xxy, tr_0_0_x_xxxy_xxz, tr_0_0_x_xxxy_xyy, tr_0_0_x_xxxy_xyz, tr_0_0_x_xxxy_xzz, tr_0_0_x_xxxy_yyy, tr_0_0_x_xxxy_yyz, tr_0_0_x_xxxy_yzz, tr_0_0_x_xxxy_zzz, tr_xxxxy_xxx, tr_xxxxy_xxy, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_zzz, tr_xxxy_xx, tr_xxxy_xxxx, tr_xxxy_xxxy, tr_xxxy_xxxz, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yz, tr_xxxy_zz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxy_xxx[i] = 2.0 * tr_xxxxy_xxx[i] * tbe_0 + 2.0 * tr_xxxy_xxxx[i] * tke_0 - 3.0 * tr_xxy_xxx[i] - 3.0 * tr_xxxy_xx[i];

        tr_0_0_x_xxxy_xxy[i] = 2.0 * tr_xxxxy_xxy[i] * tbe_0 + 2.0 * tr_xxxy_xxxy[i] * tke_0 - 3.0 * tr_xxy_xxy[i] - 2.0 * tr_xxxy_xy[i];

        tr_0_0_x_xxxy_xxz[i] = 2.0 * tr_xxxxy_xxz[i] * tbe_0 + 2.0 * tr_xxxy_xxxz[i] * tke_0 - 3.0 * tr_xxy_xxz[i] - 2.0 * tr_xxxy_xz[i];

        tr_0_0_x_xxxy_xyy[i] = 2.0 * tr_xxxxy_xyy[i] * tbe_0 + 2.0 * tr_xxxy_xxyy[i] * tke_0 - 3.0 * tr_xxy_xyy[i] - tr_xxxy_yy[i];

        tr_0_0_x_xxxy_xyz[i] = 2.0 * tr_xxxxy_xyz[i] * tbe_0 + 2.0 * tr_xxxy_xxyz[i] * tke_0 - 3.0 * tr_xxy_xyz[i] - tr_xxxy_yz[i];

        tr_0_0_x_xxxy_xzz[i] = 2.0 * tr_xxxxy_xzz[i] * tbe_0 + 2.0 * tr_xxxy_xxzz[i] * tke_0 - 3.0 * tr_xxy_xzz[i] - tr_xxxy_zz[i];

        tr_0_0_x_xxxy_yyy[i] = 2.0 * tr_xxxxy_yyy[i] * tbe_0 + 2.0 * tr_xxxy_xyyy[i] * tke_0 - 3.0 * tr_xxy_yyy[i];

        tr_0_0_x_xxxy_yyz[i] = 2.0 * tr_xxxxy_yyz[i] * tbe_0 + 2.0 * tr_xxxy_xyyz[i] * tke_0 - 3.0 * tr_xxy_yyz[i];

        tr_0_0_x_xxxy_yzz[i] = 2.0 * tr_xxxxy_yzz[i] * tbe_0 + 2.0 * tr_xxxy_xyzz[i] * tke_0 - 3.0 * tr_xxy_yzz[i];

        tr_0_0_x_xxxy_zzz[i] = 2.0 * tr_xxxxy_zzz[i] * tbe_0 + 2.0 * tr_xxxy_xzzz[i] * tke_0 - 3.0 * tr_xxy_zzz[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto tr_0_0_x_xxxz_xxx = pbuffer.data(idx_op_geom_010_gf + 20);

    auto tr_0_0_x_xxxz_xxy = pbuffer.data(idx_op_geom_010_gf + 21);

    auto tr_0_0_x_xxxz_xxz = pbuffer.data(idx_op_geom_010_gf + 22);

    auto tr_0_0_x_xxxz_xyy = pbuffer.data(idx_op_geom_010_gf + 23);

    auto tr_0_0_x_xxxz_xyz = pbuffer.data(idx_op_geom_010_gf + 24);

    auto tr_0_0_x_xxxz_xzz = pbuffer.data(idx_op_geom_010_gf + 25);

    auto tr_0_0_x_xxxz_yyy = pbuffer.data(idx_op_geom_010_gf + 26);

    auto tr_0_0_x_xxxz_yyz = pbuffer.data(idx_op_geom_010_gf + 27);

    auto tr_0_0_x_xxxz_yzz = pbuffer.data(idx_op_geom_010_gf + 28);

    auto tr_0_0_x_xxxz_zzz = pbuffer.data(idx_op_geom_010_gf + 29);

    #pragma omp simd aligned(tr_0_0_x_xxxz_xxx, tr_0_0_x_xxxz_xxy, tr_0_0_x_xxxz_xxz, tr_0_0_x_xxxz_xyy, tr_0_0_x_xxxz_xyz, tr_0_0_x_xxxz_xzz, tr_0_0_x_xxxz_yyy, tr_0_0_x_xxxz_yyz, tr_0_0_x_xxxz_yzz, tr_0_0_x_xxxz_zzz, tr_xxxxz_xxx, tr_xxxxz_xxy, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_zzz, tr_xxxz_xx, tr_xxxz_xxxx, tr_xxxz_xxxy, tr_xxxz_xxxz, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yz, tr_xxxz_zz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxz_xxx[i] = 2.0 * tr_xxxxz_xxx[i] * tbe_0 + 2.0 * tr_xxxz_xxxx[i] * tke_0 - 3.0 * tr_xxz_xxx[i] - 3.0 * tr_xxxz_xx[i];

        tr_0_0_x_xxxz_xxy[i] = 2.0 * tr_xxxxz_xxy[i] * tbe_0 + 2.0 * tr_xxxz_xxxy[i] * tke_0 - 3.0 * tr_xxz_xxy[i] - 2.0 * tr_xxxz_xy[i];

        tr_0_0_x_xxxz_xxz[i] = 2.0 * tr_xxxxz_xxz[i] * tbe_0 + 2.0 * tr_xxxz_xxxz[i] * tke_0 - 3.0 * tr_xxz_xxz[i] - 2.0 * tr_xxxz_xz[i];

        tr_0_0_x_xxxz_xyy[i] = 2.0 * tr_xxxxz_xyy[i] * tbe_0 + 2.0 * tr_xxxz_xxyy[i] * tke_0 - 3.0 * tr_xxz_xyy[i] - tr_xxxz_yy[i];

        tr_0_0_x_xxxz_xyz[i] = 2.0 * tr_xxxxz_xyz[i] * tbe_0 + 2.0 * tr_xxxz_xxyz[i] * tke_0 - 3.0 * tr_xxz_xyz[i] - tr_xxxz_yz[i];

        tr_0_0_x_xxxz_xzz[i] = 2.0 * tr_xxxxz_xzz[i] * tbe_0 + 2.0 * tr_xxxz_xxzz[i] * tke_0 - 3.0 * tr_xxz_xzz[i] - tr_xxxz_zz[i];

        tr_0_0_x_xxxz_yyy[i] = 2.0 * tr_xxxxz_yyy[i] * tbe_0 + 2.0 * tr_xxxz_xyyy[i] * tke_0 - 3.0 * tr_xxz_yyy[i];

        tr_0_0_x_xxxz_yyz[i] = 2.0 * tr_xxxxz_yyz[i] * tbe_0 + 2.0 * tr_xxxz_xyyz[i] * tke_0 - 3.0 * tr_xxz_yyz[i];

        tr_0_0_x_xxxz_yzz[i] = 2.0 * tr_xxxxz_yzz[i] * tbe_0 + 2.0 * tr_xxxz_xyzz[i] * tke_0 - 3.0 * tr_xxz_yzz[i];

        tr_0_0_x_xxxz_zzz[i] = 2.0 * tr_xxxxz_zzz[i] * tbe_0 + 2.0 * tr_xxxz_xzzz[i] * tke_0 - 3.0 * tr_xxz_zzz[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto tr_0_0_x_xxyy_xxx = pbuffer.data(idx_op_geom_010_gf + 30);

    auto tr_0_0_x_xxyy_xxy = pbuffer.data(idx_op_geom_010_gf + 31);

    auto tr_0_0_x_xxyy_xxz = pbuffer.data(idx_op_geom_010_gf + 32);

    auto tr_0_0_x_xxyy_xyy = pbuffer.data(idx_op_geom_010_gf + 33);

    auto tr_0_0_x_xxyy_xyz = pbuffer.data(idx_op_geom_010_gf + 34);

    auto tr_0_0_x_xxyy_xzz = pbuffer.data(idx_op_geom_010_gf + 35);

    auto tr_0_0_x_xxyy_yyy = pbuffer.data(idx_op_geom_010_gf + 36);

    auto tr_0_0_x_xxyy_yyz = pbuffer.data(idx_op_geom_010_gf + 37);

    auto tr_0_0_x_xxyy_yzz = pbuffer.data(idx_op_geom_010_gf + 38);

    auto tr_0_0_x_xxyy_zzz = pbuffer.data(idx_op_geom_010_gf + 39);

    #pragma omp simd aligned(tr_0_0_x_xxyy_xxx, tr_0_0_x_xxyy_xxy, tr_0_0_x_xxyy_xxz, tr_0_0_x_xxyy_xyy, tr_0_0_x_xxyy_xyz, tr_0_0_x_xxyy_xzz, tr_0_0_x_xxyy_yyy, tr_0_0_x_xxyy_yyz, tr_0_0_x_xxyy_yzz, tr_0_0_x_xxyy_zzz, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, tr_xxyy_xx, tr_xxyy_xxxx, tr_xxyy_xxxy, tr_xxyy_xxxz, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yz, tr_xxyy_zz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyy_xxx[i] = 2.0 * tr_xxxyy_xxx[i] * tbe_0 + 2.0 * tr_xxyy_xxxx[i] * tke_0 - 2.0 * tr_xyy_xxx[i] - 3.0 * tr_xxyy_xx[i];

        tr_0_0_x_xxyy_xxy[i] = 2.0 * tr_xxxyy_xxy[i] * tbe_0 + 2.0 * tr_xxyy_xxxy[i] * tke_0 - 2.0 * tr_xyy_xxy[i] - 2.0 * tr_xxyy_xy[i];

        tr_0_0_x_xxyy_xxz[i] = 2.0 * tr_xxxyy_xxz[i] * tbe_0 + 2.0 * tr_xxyy_xxxz[i] * tke_0 - 2.0 * tr_xyy_xxz[i] - 2.0 * tr_xxyy_xz[i];

        tr_0_0_x_xxyy_xyy[i] = 2.0 * tr_xxxyy_xyy[i] * tbe_0 + 2.0 * tr_xxyy_xxyy[i] * tke_0 - 2.0 * tr_xyy_xyy[i] - tr_xxyy_yy[i];

        tr_0_0_x_xxyy_xyz[i] = 2.0 * tr_xxxyy_xyz[i] * tbe_0 + 2.0 * tr_xxyy_xxyz[i] * tke_0 - 2.0 * tr_xyy_xyz[i] - tr_xxyy_yz[i];

        tr_0_0_x_xxyy_xzz[i] = 2.0 * tr_xxxyy_xzz[i] * tbe_0 + 2.0 * tr_xxyy_xxzz[i] * tke_0 - 2.0 * tr_xyy_xzz[i] - tr_xxyy_zz[i];

        tr_0_0_x_xxyy_yyy[i] = 2.0 * tr_xxxyy_yyy[i] * tbe_0 + 2.0 * tr_xxyy_xyyy[i] * tke_0 - 2.0 * tr_xyy_yyy[i];

        tr_0_0_x_xxyy_yyz[i] = 2.0 * tr_xxxyy_yyz[i] * tbe_0 + 2.0 * tr_xxyy_xyyz[i] * tke_0 - 2.0 * tr_xyy_yyz[i];

        tr_0_0_x_xxyy_yzz[i] = 2.0 * tr_xxxyy_yzz[i] * tbe_0 + 2.0 * tr_xxyy_xyzz[i] * tke_0 - 2.0 * tr_xyy_yzz[i];

        tr_0_0_x_xxyy_zzz[i] = 2.0 * tr_xxxyy_zzz[i] * tbe_0 + 2.0 * tr_xxyy_xzzz[i] * tke_0 - 2.0 * tr_xyy_zzz[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto tr_0_0_x_xxyz_xxx = pbuffer.data(idx_op_geom_010_gf + 40);

    auto tr_0_0_x_xxyz_xxy = pbuffer.data(idx_op_geom_010_gf + 41);

    auto tr_0_0_x_xxyz_xxz = pbuffer.data(idx_op_geom_010_gf + 42);

    auto tr_0_0_x_xxyz_xyy = pbuffer.data(idx_op_geom_010_gf + 43);

    auto tr_0_0_x_xxyz_xyz = pbuffer.data(idx_op_geom_010_gf + 44);

    auto tr_0_0_x_xxyz_xzz = pbuffer.data(idx_op_geom_010_gf + 45);

    auto tr_0_0_x_xxyz_yyy = pbuffer.data(idx_op_geom_010_gf + 46);

    auto tr_0_0_x_xxyz_yyz = pbuffer.data(idx_op_geom_010_gf + 47);

    auto tr_0_0_x_xxyz_yzz = pbuffer.data(idx_op_geom_010_gf + 48);

    auto tr_0_0_x_xxyz_zzz = pbuffer.data(idx_op_geom_010_gf + 49);

    #pragma omp simd aligned(tr_0_0_x_xxyz_xxx, tr_0_0_x_xxyz_xxy, tr_0_0_x_xxyz_xxz, tr_0_0_x_xxyz_xyy, tr_0_0_x_xxyz_xyz, tr_0_0_x_xxyz_xzz, tr_0_0_x_xxyz_yyy, tr_0_0_x_xxyz_yyz, tr_0_0_x_xxyz_yzz, tr_0_0_x_xxyz_zzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxyz_xx, tr_xxyz_xxxx, tr_xxyz_xxxy, tr_xxyz_xxxz, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yz, tr_xxyz_zz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyz_xxx[i] = 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 2.0 * tr_xxyz_xxxx[i] * tke_0 - 2.0 * tr_xyz_xxx[i] - 3.0 * tr_xxyz_xx[i];

        tr_0_0_x_xxyz_xxy[i] = 2.0 * tr_xxxyz_xxy[i] * tbe_0 + 2.0 * tr_xxyz_xxxy[i] * tke_0 - 2.0 * tr_xyz_xxy[i] - 2.0 * tr_xxyz_xy[i];

        tr_0_0_x_xxyz_xxz[i] = 2.0 * tr_xxxyz_xxz[i] * tbe_0 + 2.0 * tr_xxyz_xxxz[i] * tke_0 - 2.0 * tr_xyz_xxz[i] - 2.0 * tr_xxyz_xz[i];

        tr_0_0_x_xxyz_xyy[i] = 2.0 * tr_xxxyz_xyy[i] * tbe_0 + 2.0 * tr_xxyz_xxyy[i] * tke_0 - 2.0 * tr_xyz_xyy[i] - tr_xxyz_yy[i];

        tr_0_0_x_xxyz_xyz[i] = 2.0 * tr_xxxyz_xyz[i] * tbe_0 + 2.0 * tr_xxyz_xxyz[i] * tke_0 - 2.0 * tr_xyz_xyz[i] - tr_xxyz_yz[i];

        tr_0_0_x_xxyz_xzz[i] = 2.0 * tr_xxxyz_xzz[i] * tbe_0 + 2.0 * tr_xxyz_xxzz[i] * tke_0 - 2.0 * tr_xyz_xzz[i] - tr_xxyz_zz[i];

        tr_0_0_x_xxyz_yyy[i] = 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 2.0 * tr_xxyz_xyyy[i] * tke_0 - 2.0 * tr_xyz_yyy[i];

        tr_0_0_x_xxyz_yyz[i] = 2.0 * tr_xxxyz_yyz[i] * tbe_0 + 2.0 * tr_xxyz_xyyz[i] * tke_0 - 2.0 * tr_xyz_yyz[i];

        tr_0_0_x_xxyz_yzz[i] = 2.0 * tr_xxxyz_yzz[i] * tbe_0 + 2.0 * tr_xxyz_xyzz[i] * tke_0 - 2.0 * tr_xyz_yzz[i];

        tr_0_0_x_xxyz_zzz[i] = 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 2.0 * tr_xxyz_xzzz[i] * tke_0 - 2.0 * tr_xyz_zzz[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto tr_0_0_x_xxzz_xxx = pbuffer.data(idx_op_geom_010_gf + 50);

    auto tr_0_0_x_xxzz_xxy = pbuffer.data(idx_op_geom_010_gf + 51);

    auto tr_0_0_x_xxzz_xxz = pbuffer.data(idx_op_geom_010_gf + 52);

    auto tr_0_0_x_xxzz_xyy = pbuffer.data(idx_op_geom_010_gf + 53);

    auto tr_0_0_x_xxzz_xyz = pbuffer.data(idx_op_geom_010_gf + 54);

    auto tr_0_0_x_xxzz_xzz = pbuffer.data(idx_op_geom_010_gf + 55);

    auto tr_0_0_x_xxzz_yyy = pbuffer.data(idx_op_geom_010_gf + 56);

    auto tr_0_0_x_xxzz_yyz = pbuffer.data(idx_op_geom_010_gf + 57);

    auto tr_0_0_x_xxzz_yzz = pbuffer.data(idx_op_geom_010_gf + 58);

    auto tr_0_0_x_xxzz_zzz = pbuffer.data(idx_op_geom_010_gf + 59);

    #pragma omp simd aligned(tr_0_0_x_xxzz_xxx, tr_0_0_x_xxzz_xxy, tr_0_0_x_xxzz_xxz, tr_0_0_x_xxzz_xyy, tr_0_0_x_xxzz_xyz, tr_0_0_x_xxzz_xzz, tr_0_0_x_xxzz_yyy, tr_0_0_x_xxzz_yyz, tr_0_0_x_xxzz_yzz, tr_0_0_x_xxzz_zzz, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, tr_xxzz_xx, tr_xxzz_xxxx, tr_xxzz_xxxy, tr_xxzz_xxxz, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yz, tr_xxzz_zz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxzz_xxx[i] = 2.0 * tr_xxxzz_xxx[i] * tbe_0 + 2.0 * tr_xxzz_xxxx[i] * tke_0 - 2.0 * tr_xzz_xxx[i] - 3.0 * tr_xxzz_xx[i];

        tr_0_0_x_xxzz_xxy[i] = 2.0 * tr_xxxzz_xxy[i] * tbe_0 + 2.0 * tr_xxzz_xxxy[i] * tke_0 - 2.0 * tr_xzz_xxy[i] - 2.0 * tr_xxzz_xy[i];

        tr_0_0_x_xxzz_xxz[i] = 2.0 * tr_xxxzz_xxz[i] * tbe_0 + 2.0 * tr_xxzz_xxxz[i] * tke_0 - 2.0 * tr_xzz_xxz[i] - 2.0 * tr_xxzz_xz[i];

        tr_0_0_x_xxzz_xyy[i] = 2.0 * tr_xxxzz_xyy[i] * tbe_0 + 2.0 * tr_xxzz_xxyy[i] * tke_0 - 2.0 * tr_xzz_xyy[i] - tr_xxzz_yy[i];

        tr_0_0_x_xxzz_xyz[i] = 2.0 * tr_xxxzz_xyz[i] * tbe_0 + 2.0 * tr_xxzz_xxyz[i] * tke_0 - 2.0 * tr_xzz_xyz[i] - tr_xxzz_yz[i];

        tr_0_0_x_xxzz_xzz[i] = 2.0 * tr_xxxzz_xzz[i] * tbe_0 + 2.0 * tr_xxzz_xxzz[i] * tke_0 - 2.0 * tr_xzz_xzz[i] - tr_xxzz_zz[i];

        tr_0_0_x_xxzz_yyy[i] = 2.0 * tr_xxxzz_yyy[i] * tbe_0 + 2.0 * tr_xxzz_xyyy[i] * tke_0 - 2.0 * tr_xzz_yyy[i];

        tr_0_0_x_xxzz_yyz[i] = 2.0 * tr_xxxzz_yyz[i] * tbe_0 + 2.0 * tr_xxzz_xyyz[i] * tke_0 - 2.0 * tr_xzz_yyz[i];

        tr_0_0_x_xxzz_yzz[i] = 2.0 * tr_xxxzz_yzz[i] * tbe_0 + 2.0 * tr_xxzz_xyzz[i] * tke_0 - 2.0 * tr_xzz_yzz[i];

        tr_0_0_x_xxzz_zzz[i] = 2.0 * tr_xxxzz_zzz[i] * tbe_0 + 2.0 * tr_xxzz_xzzz[i] * tke_0 - 2.0 * tr_xzz_zzz[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto tr_0_0_x_xyyy_xxx = pbuffer.data(idx_op_geom_010_gf + 60);

    auto tr_0_0_x_xyyy_xxy = pbuffer.data(idx_op_geom_010_gf + 61);

    auto tr_0_0_x_xyyy_xxz = pbuffer.data(idx_op_geom_010_gf + 62);

    auto tr_0_0_x_xyyy_xyy = pbuffer.data(idx_op_geom_010_gf + 63);

    auto tr_0_0_x_xyyy_xyz = pbuffer.data(idx_op_geom_010_gf + 64);

    auto tr_0_0_x_xyyy_xzz = pbuffer.data(idx_op_geom_010_gf + 65);

    auto tr_0_0_x_xyyy_yyy = pbuffer.data(idx_op_geom_010_gf + 66);

    auto tr_0_0_x_xyyy_yyz = pbuffer.data(idx_op_geom_010_gf + 67);

    auto tr_0_0_x_xyyy_yzz = pbuffer.data(idx_op_geom_010_gf + 68);

    auto tr_0_0_x_xyyy_zzz = pbuffer.data(idx_op_geom_010_gf + 69);

    #pragma omp simd aligned(tr_0_0_x_xyyy_xxx, tr_0_0_x_xyyy_xxy, tr_0_0_x_xyyy_xxz, tr_0_0_x_xyyy_xyy, tr_0_0_x_xyyy_xyz, tr_0_0_x_xyyy_xzz, tr_0_0_x_xyyy_yyy, tr_0_0_x_xyyy_yyz, tr_0_0_x_xyyy_yzz, tr_0_0_x_xyyy_zzz, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, tr_xyyy_xx, tr_xyyy_xxxx, tr_xyyy_xxxy, tr_xyyy_xxxz, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yz, tr_xyyy_zz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyy_xxx[i] = 2.0 * tr_xxyyy_xxx[i] * tbe_0 + 2.0 * tr_xyyy_xxxx[i] * tke_0 - tr_yyy_xxx[i] - 3.0 * tr_xyyy_xx[i];

        tr_0_0_x_xyyy_xxy[i] = 2.0 * tr_xxyyy_xxy[i] * tbe_0 + 2.0 * tr_xyyy_xxxy[i] * tke_0 - tr_yyy_xxy[i] - 2.0 * tr_xyyy_xy[i];

        tr_0_0_x_xyyy_xxz[i] = 2.0 * tr_xxyyy_xxz[i] * tbe_0 + 2.0 * tr_xyyy_xxxz[i] * tke_0 - tr_yyy_xxz[i] - 2.0 * tr_xyyy_xz[i];

        tr_0_0_x_xyyy_xyy[i] = 2.0 * tr_xxyyy_xyy[i] * tbe_0 + 2.0 * tr_xyyy_xxyy[i] * tke_0 - tr_yyy_xyy[i] - tr_xyyy_yy[i];

        tr_0_0_x_xyyy_xyz[i] = 2.0 * tr_xxyyy_xyz[i] * tbe_0 + 2.0 * tr_xyyy_xxyz[i] * tke_0 - tr_yyy_xyz[i] - tr_xyyy_yz[i];

        tr_0_0_x_xyyy_xzz[i] = 2.0 * tr_xxyyy_xzz[i] * tbe_0 + 2.0 * tr_xyyy_xxzz[i] * tke_0 - tr_yyy_xzz[i] - tr_xyyy_zz[i];

        tr_0_0_x_xyyy_yyy[i] = 2.0 * tr_xxyyy_yyy[i] * tbe_0 + 2.0 * tr_xyyy_xyyy[i] * tke_0 - tr_yyy_yyy[i];

        tr_0_0_x_xyyy_yyz[i] = 2.0 * tr_xxyyy_yyz[i] * tbe_0 + 2.0 * tr_xyyy_xyyz[i] * tke_0 - tr_yyy_yyz[i];

        tr_0_0_x_xyyy_yzz[i] = 2.0 * tr_xxyyy_yzz[i] * tbe_0 + 2.0 * tr_xyyy_xyzz[i] * tke_0 - tr_yyy_yzz[i];

        tr_0_0_x_xyyy_zzz[i] = 2.0 * tr_xxyyy_zzz[i] * tbe_0 + 2.0 * tr_xyyy_xzzz[i] * tke_0 - tr_yyy_zzz[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto tr_0_0_x_xyyz_xxx = pbuffer.data(idx_op_geom_010_gf + 70);

    auto tr_0_0_x_xyyz_xxy = pbuffer.data(idx_op_geom_010_gf + 71);

    auto tr_0_0_x_xyyz_xxz = pbuffer.data(idx_op_geom_010_gf + 72);

    auto tr_0_0_x_xyyz_xyy = pbuffer.data(idx_op_geom_010_gf + 73);

    auto tr_0_0_x_xyyz_xyz = pbuffer.data(idx_op_geom_010_gf + 74);

    auto tr_0_0_x_xyyz_xzz = pbuffer.data(idx_op_geom_010_gf + 75);

    auto tr_0_0_x_xyyz_yyy = pbuffer.data(idx_op_geom_010_gf + 76);

    auto tr_0_0_x_xyyz_yyz = pbuffer.data(idx_op_geom_010_gf + 77);

    auto tr_0_0_x_xyyz_yzz = pbuffer.data(idx_op_geom_010_gf + 78);

    auto tr_0_0_x_xyyz_zzz = pbuffer.data(idx_op_geom_010_gf + 79);

    #pragma omp simd aligned(tr_0_0_x_xyyz_xxx, tr_0_0_x_xyyz_xxy, tr_0_0_x_xyyz_xxz, tr_0_0_x_xyyz_xyy, tr_0_0_x_xyyz_xyz, tr_0_0_x_xyyz_xzz, tr_0_0_x_xyyz_yyy, tr_0_0_x_xyyz_yyz, tr_0_0_x_xyyz_yzz, tr_0_0_x_xyyz_zzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxx, tr_xyyz_xxxy, tr_xyyz_xxxz, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yz, tr_xyyz_zz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyz_xxx[i] = 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 2.0 * tr_xyyz_xxxx[i] * tke_0 - tr_yyz_xxx[i] - 3.0 * tr_xyyz_xx[i];

        tr_0_0_x_xyyz_xxy[i] = 2.0 * tr_xxyyz_xxy[i] * tbe_0 + 2.0 * tr_xyyz_xxxy[i] * tke_0 - tr_yyz_xxy[i] - 2.0 * tr_xyyz_xy[i];

        tr_0_0_x_xyyz_xxz[i] = 2.0 * tr_xxyyz_xxz[i] * tbe_0 + 2.0 * tr_xyyz_xxxz[i] * tke_0 - tr_yyz_xxz[i] - 2.0 * tr_xyyz_xz[i];

        tr_0_0_x_xyyz_xyy[i] = 2.0 * tr_xxyyz_xyy[i] * tbe_0 + 2.0 * tr_xyyz_xxyy[i] * tke_0 - tr_yyz_xyy[i] - tr_xyyz_yy[i];

        tr_0_0_x_xyyz_xyz[i] = 2.0 * tr_xxyyz_xyz[i] * tbe_0 + 2.0 * tr_xyyz_xxyz[i] * tke_0 - tr_yyz_xyz[i] - tr_xyyz_yz[i];

        tr_0_0_x_xyyz_xzz[i] = 2.0 * tr_xxyyz_xzz[i] * tbe_0 + 2.0 * tr_xyyz_xxzz[i] * tke_0 - tr_yyz_xzz[i] - tr_xyyz_zz[i];

        tr_0_0_x_xyyz_yyy[i] = 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 2.0 * tr_xyyz_xyyy[i] * tke_0 - tr_yyz_yyy[i];

        tr_0_0_x_xyyz_yyz[i] = 2.0 * tr_xxyyz_yyz[i] * tbe_0 + 2.0 * tr_xyyz_xyyz[i] * tke_0 - tr_yyz_yyz[i];

        tr_0_0_x_xyyz_yzz[i] = 2.0 * tr_xxyyz_yzz[i] * tbe_0 + 2.0 * tr_xyyz_xyzz[i] * tke_0 - tr_yyz_yzz[i];

        tr_0_0_x_xyyz_zzz[i] = 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 2.0 * tr_xyyz_xzzz[i] * tke_0 - tr_yyz_zzz[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto tr_0_0_x_xyzz_xxx = pbuffer.data(idx_op_geom_010_gf + 80);

    auto tr_0_0_x_xyzz_xxy = pbuffer.data(idx_op_geom_010_gf + 81);

    auto tr_0_0_x_xyzz_xxz = pbuffer.data(idx_op_geom_010_gf + 82);

    auto tr_0_0_x_xyzz_xyy = pbuffer.data(idx_op_geom_010_gf + 83);

    auto tr_0_0_x_xyzz_xyz = pbuffer.data(idx_op_geom_010_gf + 84);

    auto tr_0_0_x_xyzz_xzz = pbuffer.data(idx_op_geom_010_gf + 85);

    auto tr_0_0_x_xyzz_yyy = pbuffer.data(idx_op_geom_010_gf + 86);

    auto tr_0_0_x_xyzz_yyz = pbuffer.data(idx_op_geom_010_gf + 87);

    auto tr_0_0_x_xyzz_yzz = pbuffer.data(idx_op_geom_010_gf + 88);

    auto tr_0_0_x_xyzz_zzz = pbuffer.data(idx_op_geom_010_gf + 89);

    #pragma omp simd aligned(tr_0_0_x_xyzz_xxx, tr_0_0_x_xyzz_xxy, tr_0_0_x_xyzz_xxz, tr_0_0_x_xyzz_xyy, tr_0_0_x_xyzz_xyz, tr_0_0_x_xyzz_xzz, tr_0_0_x_xyzz_yyy, tr_0_0_x_xyzz_yyz, tr_0_0_x_xyzz_yzz, tr_0_0_x_xyzz_zzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxx, tr_xyzz_xxxy, tr_xyzz_xxxz, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yz, tr_xyzz_zz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyzz_xxx[i] = 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 2.0 * tr_xyzz_xxxx[i] * tke_0 - tr_yzz_xxx[i] - 3.0 * tr_xyzz_xx[i];

        tr_0_0_x_xyzz_xxy[i] = 2.0 * tr_xxyzz_xxy[i] * tbe_0 + 2.0 * tr_xyzz_xxxy[i] * tke_0 - tr_yzz_xxy[i] - 2.0 * tr_xyzz_xy[i];

        tr_0_0_x_xyzz_xxz[i] = 2.0 * tr_xxyzz_xxz[i] * tbe_0 + 2.0 * tr_xyzz_xxxz[i] * tke_0 - tr_yzz_xxz[i] - 2.0 * tr_xyzz_xz[i];

        tr_0_0_x_xyzz_xyy[i] = 2.0 * tr_xxyzz_xyy[i] * tbe_0 + 2.0 * tr_xyzz_xxyy[i] * tke_0 - tr_yzz_xyy[i] - tr_xyzz_yy[i];

        tr_0_0_x_xyzz_xyz[i] = 2.0 * tr_xxyzz_xyz[i] * tbe_0 + 2.0 * tr_xyzz_xxyz[i] * tke_0 - tr_yzz_xyz[i] - tr_xyzz_yz[i];

        tr_0_0_x_xyzz_xzz[i] = 2.0 * tr_xxyzz_xzz[i] * tbe_0 + 2.0 * tr_xyzz_xxzz[i] * tke_0 - tr_yzz_xzz[i] - tr_xyzz_zz[i];

        tr_0_0_x_xyzz_yyy[i] = 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 2.0 * tr_xyzz_xyyy[i] * tke_0 - tr_yzz_yyy[i];

        tr_0_0_x_xyzz_yyz[i] = 2.0 * tr_xxyzz_yyz[i] * tbe_0 + 2.0 * tr_xyzz_xyyz[i] * tke_0 - tr_yzz_yyz[i];

        tr_0_0_x_xyzz_yzz[i] = 2.0 * tr_xxyzz_yzz[i] * tbe_0 + 2.0 * tr_xyzz_xyzz[i] * tke_0 - tr_yzz_yzz[i];

        tr_0_0_x_xyzz_zzz[i] = 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 2.0 * tr_xyzz_xzzz[i] * tke_0 - tr_yzz_zzz[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto tr_0_0_x_xzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 90);

    auto tr_0_0_x_xzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 91);

    auto tr_0_0_x_xzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 92);

    auto tr_0_0_x_xzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 93);

    auto tr_0_0_x_xzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 94);

    auto tr_0_0_x_xzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 95);

    auto tr_0_0_x_xzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 96);

    auto tr_0_0_x_xzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 97);

    auto tr_0_0_x_xzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 98);

    auto tr_0_0_x_xzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 99);

    #pragma omp simd aligned(tr_0_0_x_xzzz_xxx, tr_0_0_x_xzzz_xxy, tr_0_0_x_xzzz_xxz, tr_0_0_x_xzzz_xyy, tr_0_0_x_xzzz_xyz, tr_0_0_x_xzzz_xzz, tr_0_0_x_xzzz_yyy, tr_0_0_x_xzzz_yyz, tr_0_0_x_xzzz_yzz, tr_0_0_x_xzzz_zzz, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxx, tr_xzzz_xxxy, tr_xzzz_xxxz, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yz, tr_xzzz_zz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzzz_xxx[i] = 2.0 * tr_xxzzz_xxx[i] * tbe_0 + 2.0 * tr_xzzz_xxxx[i] * tke_0 - tr_zzz_xxx[i] - 3.0 * tr_xzzz_xx[i];

        tr_0_0_x_xzzz_xxy[i] = 2.0 * tr_xxzzz_xxy[i] * tbe_0 + 2.0 * tr_xzzz_xxxy[i] * tke_0 - tr_zzz_xxy[i] - 2.0 * tr_xzzz_xy[i];

        tr_0_0_x_xzzz_xxz[i] = 2.0 * tr_xxzzz_xxz[i] * tbe_0 + 2.0 * tr_xzzz_xxxz[i] * tke_0 - tr_zzz_xxz[i] - 2.0 * tr_xzzz_xz[i];

        tr_0_0_x_xzzz_xyy[i] = 2.0 * tr_xxzzz_xyy[i] * tbe_0 + 2.0 * tr_xzzz_xxyy[i] * tke_0 - tr_zzz_xyy[i] - tr_xzzz_yy[i];

        tr_0_0_x_xzzz_xyz[i] = 2.0 * tr_xxzzz_xyz[i] * tbe_0 + 2.0 * tr_xzzz_xxyz[i] * tke_0 - tr_zzz_xyz[i] - tr_xzzz_yz[i];

        tr_0_0_x_xzzz_xzz[i] = 2.0 * tr_xxzzz_xzz[i] * tbe_0 + 2.0 * tr_xzzz_xxzz[i] * tke_0 - tr_zzz_xzz[i] - tr_xzzz_zz[i];

        tr_0_0_x_xzzz_yyy[i] = 2.0 * tr_xxzzz_yyy[i] * tbe_0 + 2.0 * tr_xzzz_xyyy[i] * tke_0 - tr_zzz_yyy[i];

        tr_0_0_x_xzzz_yyz[i] = 2.0 * tr_xxzzz_yyz[i] * tbe_0 + 2.0 * tr_xzzz_xyyz[i] * tke_0 - tr_zzz_yyz[i];

        tr_0_0_x_xzzz_yzz[i] = 2.0 * tr_xxzzz_yzz[i] * tbe_0 + 2.0 * tr_xzzz_xyzz[i] * tke_0 - tr_zzz_yzz[i];

        tr_0_0_x_xzzz_zzz[i] = 2.0 * tr_xxzzz_zzz[i] * tbe_0 + 2.0 * tr_xzzz_xzzz[i] * tke_0 - tr_zzz_zzz[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto tr_0_0_x_yyyy_xxx = pbuffer.data(idx_op_geom_010_gf + 100);

    auto tr_0_0_x_yyyy_xxy = pbuffer.data(idx_op_geom_010_gf + 101);

    auto tr_0_0_x_yyyy_xxz = pbuffer.data(idx_op_geom_010_gf + 102);

    auto tr_0_0_x_yyyy_xyy = pbuffer.data(idx_op_geom_010_gf + 103);

    auto tr_0_0_x_yyyy_xyz = pbuffer.data(idx_op_geom_010_gf + 104);

    auto tr_0_0_x_yyyy_xzz = pbuffer.data(idx_op_geom_010_gf + 105);

    auto tr_0_0_x_yyyy_yyy = pbuffer.data(idx_op_geom_010_gf + 106);

    auto tr_0_0_x_yyyy_yyz = pbuffer.data(idx_op_geom_010_gf + 107);

    auto tr_0_0_x_yyyy_yzz = pbuffer.data(idx_op_geom_010_gf + 108);

    auto tr_0_0_x_yyyy_zzz = pbuffer.data(idx_op_geom_010_gf + 109);

    #pragma omp simd aligned(tr_0_0_x_yyyy_xxx, tr_0_0_x_yyyy_xxy, tr_0_0_x_yyyy_xxz, tr_0_0_x_yyyy_xyy, tr_0_0_x_yyyy_xyz, tr_0_0_x_yyyy_xzz, tr_0_0_x_yyyy_yyy, tr_0_0_x_yyyy_yyz, tr_0_0_x_yyyy_yzz, tr_0_0_x_yyyy_zzz, tr_xyyyy_xxx, tr_xyyyy_xxy, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_zzz, tr_yyyy_xx, tr_yyyy_xxxx, tr_yyyy_xxxy, tr_yyyy_xxxz, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yz, tr_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyy_xxx[i] = 2.0 * tr_xyyyy_xxx[i] * tbe_0 + 2.0 * tr_yyyy_xxxx[i] * tke_0 - 3.0 * tr_yyyy_xx[i];

        tr_0_0_x_yyyy_xxy[i] = 2.0 * tr_xyyyy_xxy[i] * tbe_0 + 2.0 * tr_yyyy_xxxy[i] * tke_0 - 2.0 * tr_yyyy_xy[i];

        tr_0_0_x_yyyy_xxz[i] = 2.0 * tr_xyyyy_xxz[i] * tbe_0 + 2.0 * tr_yyyy_xxxz[i] * tke_0 - 2.0 * tr_yyyy_xz[i];

        tr_0_0_x_yyyy_xyy[i] = 2.0 * tr_xyyyy_xyy[i] * tbe_0 + 2.0 * tr_yyyy_xxyy[i] * tke_0 - tr_yyyy_yy[i];

        tr_0_0_x_yyyy_xyz[i] = 2.0 * tr_xyyyy_xyz[i] * tbe_0 + 2.0 * tr_yyyy_xxyz[i] * tke_0 - tr_yyyy_yz[i];

        tr_0_0_x_yyyy_xzz[i] = 2.0 * tr_xyyyy_xzz[i] * tbe_0 + 2.0 * tr_yyyy_xxzz[i] * tke_0 - tr_yyyy_zz[i];

        tr_0_0_x_yyyy_yyy[i] = 2.0 * tr_xyyyy_yyy[i] * tbe_0 + 2.0 * tr_yyyy_xyyy[i] * tke_0;

        tr_0_0_x_yyyy_yyz[i] = 2.0 * tr_xyyyy_yyz[i] * tbe_0 + 2.0 * tr_yyyy_xyyz[i] * tke_0;

        tr_0_0_x_yyyy_yzz[i] = 2.0 * tr_xyyyy_yzz[i] * tbe_0 + 2.0 * tr_yyyy_xyzz[i] * tke_0;

        tr_0_0_x_yyyy_zzz[i] = 2.0 * tr_xyyyy_zzz[i] * tbe_0 + 2.0 * tr_yyyy_xzzz[i] * tke_0;
    }

    // Set up 110-120 components of targeted buffer : GF

    auto tr_0_0_x_yyyz_xxx = pbuffer.data(idx_op_geom_010_gf + 110);

    auto tr_0_0_x_yyyz_xxy = pbuffer.data(idx_op_geom_010_gf + 111);

    auto tr_0_0_x_yyyz_xxz = pbuffer.data(idx_op_geom_010_gf + 112);

    auto tr_0_0_x_yyyz_xyy = pbuffer.data(idx_op_geom_010_gf + 113);

    auto tr_0_0_x_yyyz_xyz = pbuffer.data(idx_op_geom_010_gf + 114);

    auto tr_0_0_x_yyyz_xzz = pbuffer.data(idx_op_geom_010_gf + 115);

    auto tr_0_0_x_yyyz_yyy = pbuffer.data(idx_op_geom_010_gf + 116);

    auto tr_0_0_x_yyyz_yyz = pbuffer.data(idx_op_geom_010_gf + 117);

    auto tr_0_0_x_yyyz_yzz = pbuffer.data(idx_op_geom_010_gf + 118);

    auto tr_0_0_x_yyyz_zzz = pbuffer.data(idx_op_geom_010_gf + 119);

    #pragma omp simd aligned(tr_0_0_x_yyyz_xxx, tr_0_0_x_yyyz_xxy, tr_0_0_x_yyyz_xxz, tr_0_0_x_yyyz_xyy, tr_0_0_x_yyyz_xyz, tr_0_0_x_yyyz_xzz, tr_0_0_x_yyyz_yyy, tr_0_0_x_yyyz_yyz, tr_0_0_x_yyyz_yzz, tr_0_0_x_yyyz_zzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxx, tr_yyyz_xxxy, tr_yyyz_xxxz, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yz, tr_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyz_xxx[i] = 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 2.0 * tr_yyyz_xxxx[i] * tke_0 - 3.0 * tr_yyyz_xx[i];

        tr_0_0_x_yyyz_xxy[i] = 2.0 * tr_xyyyz_xxy[i] * tbe_0 + 2.0 * tr_yyyz_xxxy[i] * tke_0 - 2.0 * tr_yyyz_xy[i];

        tr_0_0_x_yyyz_xxz[i] = 2.0 * tr_xyyyz_xxz[i] * tbe_0 + 2.0 * tr_yyyz_xxxz[i] * tke_0 - 2.0 * tr_yyyz_xz[i];

        tr_0_0_x_yyyz_xyy[i] = 2.0 * tr_xyyyz_xyy[i] * tbe_0 + 2.0 * tr_yyyz_xxyy[i] * tke_0 - tr_yyyz_yy[i];

        tr_0_0_x_yyyz_xyz[i] = 2.0 * tr_xyyyz_xyz[i] * tbe_0 + 2.0 * tr_yyyz_xxyz[i] * tke_0 - tr_yyyz_yz[i];

        tr_0_0_x_yyyz_xzz[i] = 2.0 * tr_xyyyz_xzz[i] * tbe_0 + 2.0 * tr_yyyz_xxzz[i] * tke_0 - tr_yyyz_zz[i];

        tr_0_0_x_yyyz_yyy[i] = 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 2.0 * tr_yyyz_xyyy[i] * tke_0;

        tr_0_0_x_yyyz_yyz[i] = 2.0 * tr_xyyyz_yyz[i] * tbe_0 + 2.0 * tr_yyyz_xyyz[i] * tke_0;

        tr_0_0_x_yyyz_yzz[i] = 2.0 * tr_xyyyz_yzz[i] * tbe_0 + 2.0 * tr_yyyz_xyzz[i] * tke_0;

        tr_0_0_x_yyyz_zzz[i] = 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 2.0 * tr_yyyz_xzzz[i] * tke_0;
    }

    // Set up 120-130 components of targeted buffer : GF

    auto tr_0_0_x_yyzz_xxx = pbuffer.data(idx_op_geom_010_gf + 120);

    auto tr_0_0_x_yyzz_xxy = pbuffer.data(idx_op_geom_010_gf + 121);

    auto tr_0_0_x_yyzz_xxz = pbuffer.data(idx_op_geom_010_gf + 122);

    auto tr_0_0_x_yyzz_xyy = pbuffer.data(idx_op_geom_010_gf + 123);

    auto tr_0_0_x_yyzz_xyz = pbuffer.data(idx_op_geom_010_gf + 124);

    auto tr_0_0_x_yyzz_xzz = pbuffer.data(idx_op_geom_010_gf + 125);

    auto tr_0_0_x_yyzz_yyy = pbuffer.data(idx_op_geom_010_gf + 126);

    auto tr_0_0_x_yyzz_yyz = pbuffer.data(idx_op_geom_010_gf + 127);

    auto tr_0_0_x_yyzz_yzz = pbuffer.data(idx_op_geom_010_gf + 128);

    auto tr_0_0_x_yyzz_zzz = pbuffer.data(idx_op_geom_010_gf + 129);

    #pragma omp simd aligned(tr_0_0_x_yyzz_xxx, tr_0_0_x_yyzz_xxy, tr_0_0_x_yyzz_xxz, tr_0_0_x_yyzz_xyy, tr_0_0_x_yyzz_xyz, tr_0_0_x_yyzz_xzz, tr_0_0_x_yyzz_yyy, tr_0_0_x_yyzz_yyz, tr_0_0_x_yyzz_yzz, tr_0_0_x_yyzz_zzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_yyzz_xx, tr_yyzz_xxxx, tr_yyzz_xxxy, tr_yyzz_xxxz, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yz, tr_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyzz_xxx[i] = 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 2.0 * tr_yyzz_xxxx[i] * tke_0 - 3.0 * tr_yyzz_xx[i];

        tr_0_0_x_yyzz_xxy[i] = 2.0 * tr_xyyzz_xxy[i] * tbe_0 + 2.0 * tr_yyzz_xxxy[i] * tke_0 - 2.0 * tr_yyzz_xy[i];

        tr_0_0_x_yyzz_xxz[i] = 2.0 * tr_xyyzz_xxz[i] * tbe_0 + 2.0 * tr_yyzz_xxxz[i] * tke_0 - 2.0 * tr_yyzz_xz[i];

        tr_0_0_x_yyzz_xyy[i] = 2.0 * tr_xyyzz_xyy[i] * tbe_0 + 2.0 * tr_yyzz_xxyy[i] * tke_0 - tr_yyzz_yy[i];

        tr_0_0_x_yyzz_xyz[i] = 2.0 * tr_xyyzz_xyz[i] * tbe_0 + 2.0 * tr_yyzz_xxyz[i] * tke_0 - tr_yyzz_yz[i];

        tr_0_0_x_yyzz_xzz[i] = 2.0 * tr_xyyzz_xzz[i] * tbe_0 + 2.0 * tr_yyzz_xxzz[i] * tke_0 - tr_yyzz_zz[i];

        tr_0_0_x_yyzz_yyy[i] = 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 2.0 * tr_yyzz_xyyy[i] * tke_0;

        tr_0_0_x_yyzz_yyz[i] = 2.0 * tr_xyyzz_yyz[i] * tbe_0 + 2.0 * tr_yyzz_xyyz[i] * tke_0;

        tr_0_0_x_yyzz_yzz[i] = 2.0 * tr_xyyzz_yzz[i] * tbe_0 + 2.0 * tr_yyzz_xyzz[i] * tke_0;

        tr_0_0_x_yyzz_zzz[i] = 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 2.0 * tr_yyzz_xzzz[i] * tke_0;
    }

    // Set up 130-140 components of targeted buffer : GF

    auto tr_0_0_x_yzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 130);

    auto tr_0_0_x_yzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 131);

    auto tr_0_0_x_yzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 132);

    auto tr_0_0_x_yzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 133);

    auto tr_0_0_x_yzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 134);

    auto tr_0_0_x_yzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 135);

    auto tr_0_0_x_yzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 136);

    auto tr_0_0_x_yzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 137);

    auto tr_0_0_x_yzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 138);

    auto tr_0_0_x_yzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 139);

    #pragma omp simd aligned(tr_0_0_x_yzzz_xxx, tr_0_0_x_yzzz_xxy, tr_0_0_x_yzzz_xxz, tr_0_0_x_yzzz_xyy, tr_0_0_x_yzzz_xyz, tr_0_0_x_yzzz_xzz, tr_0_0_x_yzzz_yyy, tr_0_0_x_yzzz_yyz, tr_0_0_x_yzzz_yzz, tr_0_0_x_yzzz_zzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_yzzz_xx, tr_yzzz_xxxx, tr_yzzz_xxxy, tr_yzzz_xxxz, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yz, tr_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzzz_xxx[i] = 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 2.0 * tr_yzzz_xxxx[i] * tke_0 - 3.0 * tr_yzzz_xx[i];

        tr_0_0_x_yzzz_xxy[i] = 2.0 * tr_xyzzz_xxy[i] * tbe_0 + 2.0 * tr_yzzz_xxxy[i] * tke_0 - 2.0 * tr_yzzz_xy[i];

        tr_0_0_x_yzzz_xxz[i] = 2.0 * tr_xyzzz_xxz[i] * tbe_0 + 2.0 * tr_yzzz_xxxz[i] * tke_0 - 2.0 * tr_yzzz_xz[i];

        tr_0_0_x_yzzz_xyy[i] = 2.0 * tr_xyzzz_xyy[i] * tbe_0 + 2.0 * tr_yzzz_xxyy[i] * tke_0 - tr_yzzz_yy[i];

        tr_0_0_x_yzzz_xyz[i] = 2.0 * tr_xyzzz_xyz[i] * tbe_0 + 2.0 * tr_yzzz_xxyz[i] * tke_0 - tr_yzzz_yz[i];

        tr_0_0_x_yzzz_xzz[i] = 2.0 * tr_xyzzz_xzz[i] * tbe_0 + 2.0 * tr_yzzz_xxzz[i] * tke_0 - tr_yzzz_zz[i];

        tr_0_0_x_yzzz_yyy[i] = 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 2.0 * tr_yzzz_xyyy[i] * tke_0;

        tr_0_0_x_yzzz_yyz[i] = 2.0 * tr_xyzzz_yyz[i] * tbe_0 + 2.0 * tr_yzzz_xyyz[i] * tke_0;

        tr_0_0_x_yzzz_yzz[i] = 2.0 * tr_xyzzz_yzz[i] * tbe_0 + 2.0 * tr_yzzz_xyzz[i] * tke_0;

        tr_0_0_x_yzzz_zzz[i] = 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 2.0 * tr_yzzz_xzzz[i] * tke_0;
    }

    // Set up 140-150 components of targeted buffer : GF

    auto tr_0_0_x_zzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 140);

    auto tr_0_0_x_zzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 141);

    auto tr_0_0_x_zzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 142);

    auto tr_0_0_x_zzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 143);

    auto tr_0_0_x_zzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 144);

    auto tr_0_0_x_zzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 145);

    auto tr_0_0_x_zzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 146);

    auto tr_0_0_x_zzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 147);

    auto tr_0_0_x_zzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 148);

    auto tr_0_0_x_zzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 149);

    #pragma omp simd aligned(tr_0_0_x_zzzz_xxx, tr_0_0_x_zzzz_xxy, tr_0_0_x_zzzz_xxz, tr_0_0_x_zzzz_xyy, tr_0_0_x_zzzz_xyz, tr_0_0_x_zzzz_xzz, tr_0_0_x_zzzz_yyy, tr_0_0_x_zzzz_yyz, tr_0_0_x_zzzz_yzz, tr_0_0_x_zzzz_zzz, tr_xzzzz_xxx, tr_xzzzz_xxy, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_zzz, tr_zzzz_xx, tr_zzzz_xxxx, tr_zzzz_xxxy, tr_zzzz_xxxz, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzzz_xxx[i] = 2.0 * tr_xzzzz_xxx[i] * tbe_0 + 2.0 * tr_zzzz_xxxx[i] * tke_0 - 3.0 * tr_zzzz_xx[i];

        tr_0_0_x_zzzz_xxy[i] = 2.0 * tr_xzzzz_xxy[i] * tbe_0 + 2.0 * tr_zzzz_xxxy[i] * tke_0 - 2.0 * tr_zzzz_xy[i];

        tr_0_0_x_zzzz_xxz[i] = 2.0 * tr_xzzzz_xxz[i] * tbe_0 + 2.0 * tr_zzzz_xxxz[i] * tke_0 - 2.0 * tr_zzzz_xz[i];

        tr_0_0_x_zzzz_xyy[i] = 2.0 * tr_xzzzz_xyy[i] * tbe_0 + 2.0 * tr_zzzz_xxyy[i] * tke_0 - tr_zzzz_yy[i];

        tr_0_0_x_zzzz_xyz[i] = 2.0 * tr_xzzzz_xyz[i] * tbe_0 + 2.0 * tr_zzzz_xxyz[i] * tke_0 - tr_zzzz_yz[i];

        tr_0_0_x_zzzz_xzz[i] = 2.0 * tr_xzzzz_xzz[i] * tbe_0 + 2.0 * tr_zzzz_xxzz[i] * tke_0 - tr_zzzz_zz[i];

        tr_0_0_x_zzzz_yyy[i] = 2.0 * tr_xzzzz_yyy[i] * tbe_0 + 2.0 * tr_zzzz_xyyy[i] * tke_0;

        tr_0_0_x_zzzz_yyz[i] = 2.0 * tr_xzzzz_yyz[i] * tbe_0 + 2.0 * tr_zzzz_xyyz[i] * tke_0;

        tr_0_0_x_zzzz_yzz[i] = 2.0 * tr_xzzzz_yzz[i] * tbe_0 + 2.0 * tr_zzzz_xyzz[i] * tke_0;

        tr_0_0_x_zzzz_zzz[i] = 2.0 * tr_xzzzz_zzz[i] * tbe_0 + 2.0 * tr_zzzz_xzzz[i] * tke_0;
    }

    // Set up 150-160 components of targeted buffer : GF

    auto tr_0_0_y_xxxx_xxx = pbuffer.data(idx_op_geom_010_gf + 150);

    auto tr_0_0_y_xxxx_xxy = pbuffer.data(idx_op_geom_010_gf + 151);

    auto tr_0_0_y_xxxx_xxz = pbuffer.data(idx_op_geom_010_gf + 152);

    auto tr_0_0_y_xxxx_xyy = pbuffer.data(idx_op_geom_010_gf + 153);

    auto tr_0_0_y_xxxx_xyz = pbuffer.data(idx_op_geom_010_gf + 154);

    auto tr_0_0_y_xxxx_xzz = pbuffer.data(idx_op_geom_010_gf + 155);

    auto tr_0_0_y_xxxx_yyy = pbuffer.data(idx_op_geom_010_gf + 156);

    auto tr_0_0_y_xxxx_yyz = pbuffer.data(idx_op_geom_010_gf + 157);

    auto tr_0_0_y_xxxx_yzz = pbuffer.data(idx_op_geom_010_gf + 158);

    auto tr_0_0_y_xxxx_zzz = pbuffer.data(idx_op_geom_010_gf + 159);

    #pragma omp simd aligned(tr_0_0_y_xxxx_xxx, tr_0_0_y_xxxx_xxy, tr_0_0_y_xxxx_xxz, tr_0_0_y_xxxx_xyy, tr_0_0_y_xxxx_xyz, tr_0_0_y_xxxx_xzz, tr_0_0_y_xxxx_yyy, tr_0_0_y_xxxx_yyz, tr_0_0_y_xxxx_yzz, tr_0_0_y_xxxx_zzz, tr_xxxx_xx, tr_xxxx_xxxy, tr_xxxx_xxyy, tr_xxxx_xxyz, tr_xxxx_xy, tr_xxxx_xyyy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_yy, tr_xxxx_yyyy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxxy_xxx, tr_xxxxy_xxy, tr_xxxxy_xxz, tr_xxxxy_xyy, tr_xxxxy_xyz, tr_xxxxy_xzz, tr_xxxxy_yyy, tr_xxxxy_yyz, tr_xxxxy_yzz, tr_xxxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxx_xxx[i] = 2.0 * tr_xxxxy_xxx[i] * tbe_0 + 2.0 * tr_xxxx_xxxy[i] * tke_0;

        tr_0_0_y_xxxx_xxy[i] = 2.0 * tr_xxxxy_xxy[i] * tbe_0 + 2.0 * tr_xxxx_xxyy[i] * tke_0 - tr_xxxx_xx[i];

        tr_0_0_y_xxxx_xxz[i] = 2.0 * tr_xxxxy_xxz[i] * tbe_0 + 2.0 * tr_xxxx_xxyz[i] * tke_0;

        tr_0_0_y_xxxx_xyy[i] = 2.0 * tr_xxxxy_xyy[i] * tbe_0 + 2.0 * tr_xxxx_xyyy[i] * tke_0 - 2.0 * tr_xxxx_xy[i];

        tr_0_0_y_xxxx_xyz[i] = 2.0 * tr_xxxxy_xyz[i] * tbe_0 + 2.0 * tr_xxxx_xyyz[i] * tke_0 - tr_xxxx_xz[i];

        tr_0_0_y_xxxx_xzz[i] = 2.0 * tr_xxxxy_xzz[i] * tbe_0 + 2.0 * tr_xxxx_xyzz[i] * tke_0;

        tr_0_0_y_xxxx_yyy[i] = 2.0 * tr_xxxxy_yyy[i] * tbe_0 + 2.0 * tr_xxxx_yyyy[i] * tke_0 - 3.0 * tr_xxxx_yy[i];

        tr_0_0_y_xxxx_yyz[i] = 2.0 * tr_xxxxy_yyz[i] * tbe_0 + 2.0 * tr_xxxx_yyyz[i] * tke_0 - 2.0 * tr_xxxx_yz[i];

        tr_0_0_y_xxxx_yzz[i] = 2.0 * tr_xxxxy_yzz[i] * tbe_0 + 2.0 * tr_xxxx_yyzz[i] * tke_0 - tr_xxxx_zz[i];

        tr_0_0_y_xxxx_zzz[i] = 2.0 * tr_xxxxy_zzz[i] * tbe_0 + 2.0 * tr_xxxx_yzzz[i] * tke_0;
    }

    // Set up 160-170 components of targeted buffer : GF

    auto tr_0_0_y_xxxy_xxx = pbuffer.data(idx_op_geom_010_gf + 160);

    auto tr_0_0_y_xxxy_xxy = pbuffer.data(idx_op_geom_010_gf + 161);

    auto tr_0_0_y_xxxy_xxz = pbuffer.data(idx_op_geom_010_gf + 162);

    auto tr_0_0_y_xxxy_xyy = pbuffer.data(idx_op_geom_010_gf + 163);

    auto tr_0_0_y_xxxy_xyz = pbuffer.data(idx_op_geom_010_gf + 164);

    auto tr_0_0_y_xxxy_xzz = pbuffer.data(idx_op_geom_010_gf + 165);

    auto tr_0_0_y_xxxy_yyy = pbuffer.data(idx_op_geom_010_gf + 166);

    auto tr_0_0_y_xxxy_yyz = pbuffer.data(idx_op_geom_010_gf + 167);

    auto tr_0_0_y_xxxy_yzz = pbuffer.data(idx_op_geom_010_gf + 168);

    auto tr_0_0_y_xxxy_zzz = pbuffer.data(idx_op_geom_010_gf + 169);

    #pragma omp simd aligned(tr_0_0_y_xxxy_xxx, tr_0_0_y_xxxy_xxy, tr_0_0_y_xxxy_xxz, tr_0_0_y_xxxy_xyy, tr_0_0_y_xxxy_xyz, tr_0_0_y_xxxy_xzz, tr_0_0_y_xxxy_yyy, tr_0_0_y_xxxy_yyz, tr_0_0_y_xxxy_yzz, tr_0_0_y_xxxy_zzz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxy_xx, tr_xxxy_xxxy, tr_xxxy_xxyy, tr_xxxy_xxyz, tr_xxxy_xy, tr_xxxy_xyyy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_yy, tr_xxxy_yyyy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxyy_xxx, tr_xxxyy_xxy, tr_xxxyy_xxz, tr_xxxyy_xyy, tr_xxxyy_xyz, tr_xxxyy_xzz, tr_xxxyy_yyy, tr_xxxyy_yyz, tr_xxxyy_yzz, tr_xxxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxy_xxx[i] = 2.0 * tr_xxxyy_xxx[i] * tbe_0 + 2.0 * tr_xxxy_xxxy[i] * tke_0 - tr_xxx_xxx[i];

        tr_0_0_y_xxxy_xxy[i] = 2.0 * tr_xxxyy_xxy[i] * tbe_0 + 2.0 * tr_xxxy_xxyy[i] * tke_0 - tr_xxx_xxy[i] - tr_xxxy_xx[i];

        tr_0_0_y_xxxy_xxz[i] = 2.0 * tr_xxxyy_xxz[i] * tbe_0 + 2.0 * tr_xxxy_xxyz[i] * tke_0 - tr_xxx_xxz[i];

        tr_0_0_y_xxxy_xyy[i] = 2.0 * tr_xxxyy_xyy[i] * tbe_0 + 2.0 * tr_xxxy_xyyy[i] * tke_0 - tr_xxx_xyy[i] - 2.0 * tr_xxxy_xy[i];

        tr_0_0_y_xxxy_xyz[i] = 2.0 * tr_xxxyy_xyz[i] * tbe_0 + 2.0 * tr_xxxy_xyyz[i] * tke_0 - tr_xxx_xyz[i] - tr_xxxy_xz[i];

        tr_0_0_y_xxxy_xzz[i] = 2.0 * tr_xxxyy_xzz[i] * tbe_0 + 2.0 * tr_xxxy_xyzz[i] * tke_0 - tr_xxx_xzz[i];

        tr_0_0_y_xxxy_yyy[i] = 2.0 * tr_xxxyy_yyy[i] * tbe_0 + 2.0 * tr_xxxy_yyyy[i] * tke_0 - tr_xxx_yyy[i] - 3.0 * tr_xxxy_yy[i];

        tr_0_0_y_xxxy_yyz[i] = 2.0 * tr_xxxyy_yyz[i] * tbe_0 + 2.0 * tr_xxxy_yyyz[i] * tke_0 - tr_xxx_yyz[i] - 2.0 * tr_xxxy_yz[i];

        tr_0_0_y_xxxy_yzz[i] = 2.0 * tr_xxxyy_yzz[i] * tbe_0 + 2.0 * tr_xxxy_yyzz[i] * tke_0 - tr_xxx_yzz[i] - tr_xxxy_zz[i];

        tr_0_0_y_xxxy_zzz[i] = 2.0 * tr_xxxyy_zzz[i] * tbe_0 + 2.0 * tr_xxxy_yzzz[i] * tke_0 - tr_xxx_zzz[i];
    }

    // Set up 170-180 components of targeted buffer : GF

    auto tr_0_0_y_xxxz_xxx = pbuffer.data(idx_op_geom_010_gf + 170);

    auto tr_0_0_y_xxxz_xxy = pbuffer.data(idx_op_geom_010_gf + 171);

    auto tr_0_0_y_xxxz_xxz = pbuffer.data(idx_op_geom_010_gf + 172);

    auto tr_0_0_y_xxxz_xyy = pbuffer.data(idx_op_geom_010_gf + 173);

    auto tr_0_0_y_xxxz_xyz = pbuffer.data(idx_op_geom_010_gf + 174);

    auto tr_0_0_y_xxxz_xzz = pbuffer.data(idx_op_geom_010_gf + 175);

    auto tr_0_0_y_xxxz_yyy = pbuffer.data(idx_op_geom_010_gf + 176);

    auto tr_0_0_y_xxxz_yyz = pbuffer.data(idx_op_geom_010_gf + 177);

    auto tr_0_0_y_xxxz_yzz = pbuffer.data(idx_op_geom_010_gf + 178);

    auto tr_0_0_y_xxxz_zzz = pbuffer.data(idx_op_geom_010_gf + 179);

    #pragma omp simd aligned(tr_0_0_y_xxxz_xxx, tr_0_0_y_xxxz_xxy, tr_0_0_y_xxxz_xxz, tr_0_0_y_xxxz_xyy, tr_0_0_y_xxxz_xyz, tr_0_0_y_xxxz_xzz, tr_0_0_y_xxxz_yyy, tr_0_0_y_xxxz_yyz, tr_0_0_y_xxxz_yzz, tr_0_0_y_xxxz_zzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, tr_xxxz_xx, tr_xxxz_xxxy, tr_xxxz_xxyy, tr_xxxz_xxyz, tr_xxxz_xy, tr_xxxz_xyyy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_yy, tr_xxxz_yyyy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxz_xxx[i] = 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 2.0 * tr_xxxz_xxxy[i] * tke_0;

        tr_0_0_y_xxxz_xxy[i] = 2.0 * tr_xxxyz_xxy[i] * tbe_0 + 2.0 * tr_xxxz_xxyy[i] * tke_0 - tr_xxxz_xx[i];

        tr_0_0_y_xxxz_xxz[i] = 2.0 * tr_xxxyz_xxz[i] * tbe_0 + 2.0 * tr_xxxz_xxyz[i] * tke_0;

        tr_0_0_y_xxxz_xyy[i] = 2.0 * tr_xxxyz_xyy[i] * tbe_0 + 2.0 * tr_xxxz_xyyy[i] * tke_0 - 2.0 * tr_xxxz_xy[i];

        tr_0_0_y_xxxz_xyz[i] = 2.0 * tr_xxxyz_xyz[i] * tbe_0 + 2.0 * tr_xxxz_xyyz[i] * tke_0 - tr_xxxz_xz[i];

        tr_0_0_y_xxxz_xzz[i] = 2.0 * tr_xxxyz_xzz[i] * tbe_0 + 2.0 * tr_xxxz_xyzz[i] * tke_0;

        tr_0_0_y_xxxz_yyy[i] = 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 2.0 * tr_xxxz_yyyy[i] * tke_0 - 3.0 * tr_xxxz_yy[i];

        tr_0_0_y_xxxz_yyz[i] = 2.0 * tr_xxxyz_yyz[i] * tbe_0 + 2.0 * tr_xxxz_yyyz[i] * tke_0 - 2.0 * tr_xxxz_yz[i];

        tr_0_0_y_xxxz_yzz[i] = 2.0 * tr_xxxyz_yzz[i] * tbe_0 + 2.0 * tr_xxxz_yyzz[i] * tke_0 - tr_xxxz_zz[i];

        tr_0_0_y_xxxz_zzz[i] = 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 2.0 * tr_xxxz_yzzz[i] * tke_0;
    }

    // Set up 180-190 components of targeted buffer : GF

    auto tr_0_0_y_xxyy_xxx = pbuffer.data(idx_op_geom_010_gf + 180);

    auto tr_0_0_y_xxyy_xxy = pbuffer.data(idx_op_geom_010_gf + 181);

    auto tr_0_0_y_xxyy_xxz = pbuffer.data(idx_op_geom_010_gf + 182);

    auto tr_0_0_y_xxyy_xyy = pbuffer.data(idx_op_geom_010_gf + 183);

    auto tr_0_0_y_xxyy_xyz = pbuffer.data(idx_op_geom_010_gf + 184);

    auto tr_0_0_y_xxyy_xzz = pbuffer.data(idx_op_geom_010_gf + 185);

    auto tr_0_0_y_xxyy_yyy = pbuffer.data(idx_op_geom_010_gf + 186);

    auto tr_0_0_y_xxyy_yyz = pbuffer.data(idx_op_geom_010_gf + 187);

    auto tr_0_0_y_xxyy_yzz = pbuffer.data(idx_op_geom_010_gf + 188);

    auto tr_0_0_y_xxyy_zzz = pbuffer.data(idx_op_geom_010_gf + 189);

    #pragma omp simd aligned(tr_0_0_y_xxyy_xxx, tr_0_0_y_xxyy_xxy, tr_0_0_y_xxyy_xxz, tr_0_0_y_xxyy_xyy, tr_0_0_y_xxyy_xyz, tr_0_0_y_xxyy_xzz, tr_0_0_y_xxyy_yyy, tr_0_0_y_xxyy_yyz, tr_0_0_y_xxyy_yzz, tr_0_0_y_xxyy_zzz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyy_xx, tr_xxyy_xxxy, tr_xxyy_xxyy, tr_xxyy_xxyz, tr_xxyy_xy, tr_xxyy_xyyy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_yy, tr_xxyy_yyyy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyyy_xxx, tr_xxyyy_xxy, tr_xxyyy_xxz, tr_xxyyy_xyy, tr_xxyyy_xyz, tr_xxyyy_xzz, tr_xxyyy_yyy, tr_xxyyy_yyz, tr_xxyyy_yzz, tr_xxyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyy_xxx[i] = 2.0 * tr_xxyyy_xxx[i] * tbe_0 + 2.0 * tr_xxyy_xxxy[i] * tke_0 - 2.0 * tr_xxy_xxx[i];

        tr_0_0_y_xxyy_xxy[i] = 2.0 * tr_xxyyy_xxy[i] * tbe_0 + 2.0 * tr_xxyy_xxyy[i] * tke_0 - 2.0 * tr_xxy_xxy[i] - tr_xxyy_xx[i];

        tr_0_0_y_xxyy_xxz[i] = 2.0 * tr_xxyyy_xxz[i] * tbe_0 + 2.0 * tr_xxyy_xxyz[i] * tke_0 - 2.0 * tr_xxy_xxz[i];

        tr_0_0_y_xxyy_xyy[i] = 2.0 * tr_xxyyy_xyy[i] * tbe_0 + 2.0 * tr_xxyy_xyyy[i] * tke_0 - 2.0 * tr_xxy_xyy[i] - 2.0 * tr_xxyy_xy[i];

        tr_0_0_y_xxyy_xyz[i] = 2.0 * tr_xxyyy_xyz[i] * tbe_0 + 2.0 * tr_xxyy_xyyz[i] * tke_0 - 2.0 * tr_xxy_xyz[i] - tr_xxyy_xz[i];

        tr_0_0_y_xxyy_xzz[i] = 2.0 * tr_xxyyy_xzz[i] * tbe_0 + 2.0 * tr_xxyy_xyzz[i] * tke_0 - 2.0 * tr_xxy_xzz[i];

        tr_0_0_y_xxyy_yyy[i] = 2.0 * tr_xxyyy_yyy[i] * tbe_0 + 2.0 * tr_xxyy_yyyy[i] * tke_0 - 2.0 * tr_xxy_yyy[i] - 3.0 * tr_xxyy_yy[i];

        tr_0_0_y_xxyy_yyz[i] = 2.0 * tr_xxyyy_yyz[i] * tbe_0 + 2.0 * tr_xxyy_yyyz[i] * tke_0 - 2.0 * tr_xxy_yyz[i] - 2.0 * tr_xxyy_yz[i];

        tr_0_0_y_xxyy_yzz[i] = 2.0 * tr_xxyyy_yzz[i] * tbe_0 + 2.0 * tr_xxyy_yyzz[i] * tke_0 - 2.0 * tr_xxy_yzz[i] - tr_xxyy_zz[i];

        tr_0_0_y_xxyy_zzz[i] = 2.0 * tr_xxyyy_zzz[i] * tbe_0 + 2.0 * tr_xxyy_yzzz[i] * tke_0 - 2.0 * tr_xxy_zzz[i];
    }

    // Set up 190-200 components of targeted buffer : GF

    auto tr_0_0_y_xxyz_xxx = pbuffer.data(idx_op_geom_010_gf + 190);

    auto tr_0_0_y_xxyz_xxy = pbuffer.data(idx_op_geom_010_gf + 191);

    auto tr_0_0_y_xxyz_xxz = pbuffer.data(idx_op_geom_010_gf + 192);

    auto tr_0_0_y_xxyz_xyy = pbuffer.data(idx_op_geom_010_gf + 193);

    auto tr_0_0_y_xxyz_xyz = pbuffer.data(idx_op_geom_010_gf + 194);

    auto tr_0_0_y_xxyz_xzz = pbuffer.data(idx_op_geom_010_gf + 195);

    auto tr_0_0_y_xxyz_yyy = pbuffer.data(idx_op_geom_010_gf + 196);

    auto tr_0_0_y_xxyz_yyz = pbuffer.data(idx_op_geom_010_gf + 197);

    auto tr_0_0_y_xxyz_yzz = pbuffer.data(idx_op_geom_010_gf + 198);

    auto tr_0_0_y_xxyz_zzz = pbuffer.data(idx_op_geom_010_gf + 199);

    #pragma omp simd aligned(tr_0_0_y_xxyz_xxx, tr_0_0_y_xxyz_xxy, tr_0_0_y_xxyz_xxz, tr_0_0_y_xxyz_xyy, tr_0_0_y_xxyz_xyz, tr_0_0_y_xxyz_xzz, tr_0_0_y_xxyz_yyy, tr_0_0_y_xxyz_yyz, tr_0_0_y_xxyz_yzz, tr_0_0_y_xxyz_zzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, tr_xxyz_xx, tr_xxyz_xxxy, tr_xxyz_xxyy, tr_xxyz_xxyz, tr_xxyz_xy, tr_xxyz_xyyy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_yy, tr_xxyz_yyyy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyz_xxx[i] = 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 2.0 * tr_xxyz_xxxy[i] * tke_0 - tr_xxz_xxx[i];

        tr_0_0_y_xxyz_xxy[i] = 2.0 * tr_xxyyz_xxy[i] * tbe_0 + 2.0 * tr_xxyz_xxyy[i] * tke_0 - tr_xxz_xxy[i] - tr_xxyz_xx[i];

        tr_0_0_y_xxyz_xxz[i] = 2.0 * tr_xxyyz_xxz[i] * tbe_0 + 2.0 * tr_xxyz_xxyz[i] * tke_0 - tr_xxz_xxz[i];

        tr_0_0_y_xxyz_xyy[i] = 2.0 * tr_xxyyz_xyy[i] * tbe_0 + 2.0 * tr_xxyz_xyyy[i] * tke_0 - tr_xxz_xyy[i] - 2.0 * tr_xxyz_xy[i];

        tr_0_0_y_xxyz_xyz[i] = 2.0 * tr_xxyyz_xyz[i] * tbe_0 + 2.0 * tr_xxyz_xyyz[i] * tke_0 - tr_xxz_xyz[i] - tr_xxyz_xz[i];

        tr_0_0_y_xxyz_xzz[i] = 2.0 * tr_xxyyz_xzz[i] * tbe_0 + 2.0 * tr_xxyz_xyzz[i] * tke_0 - tr_xxz_xzz[i];

        tr_0_0_y_xxyz_yyy[i] = 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 2.0 * tr_xxyz_yyyy[i] * tke_0 - tr_xxz_yyy[i] - 3.0 * tr_xxyz_yy[i];

        tr_0_0_y_xxyz_yyz[i] = 2.0 * tr_xxyyz_yyz[i] * tbe_0 + 2.0 * tr_xxyz_yyyz[i] * tke_0 - tr_xxz_yyz[i] - 2.0 * tr_xxyz_yz[i];

        tr_0_0_y_xxyz_yzz[i] = 2.0 * tr_xxyyz_yzz[i] * tbe_0 + 2.0 * tr_xxyz_yyzz[i] * tke_0 - tr_xxz_yzz[i] - tr_xxyz_zz[i];

        tr_0_0_y_xxyz_zzz[i] = 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 2.0 * tr_xxyz_yzzz[i] * tke_0 - tr_xxz_zzz[i];
    }

    // Set up 200-210 components of targeted buffer : GF

    auto tr_0_0_y_xxzz_xxx = pbuffer.data(idx_op_geom_010_gf + 200);

    auto tr_0_0_y_xxzz_xxy = pbuffer.data(idx_op_geom_010_gf + 201);

    auto tr_0_0_y_xxzz_xxz = pbuffer.data(idx_op_geom_010_gf + 202);

    auto tr_0_0_y_xxzz_xyy = pbuffer.data(idx_op_geom_010_gf + 203);

    auto tr_0_0_y_xxzz_xyz = pbuffer.data(idx_op_geom_010_gf + 204);

    auto tr_0_0_y_xxzz_xzz = pbuffer.data(idx_op_geom_010_gf + 205);

    auto tr_0_0_y_xxzz_yyy = pbuffer.data(idx_op_geom_010_gf + 206);

    auto tr_0_0_y_xxzz_yyz = pbuffer.data(idx_op_geom_010_gf + 207);

    auto tr_0_0_y_xxzz_yzz = pbuffer.data(idx_op_geom_010_gf + 208);

    auto tr_0_0_y_xxzz_zzz = pbuffer.data(idx_op_geom_010_gf + 209);

    #pragma omp simd aligned(tr_0_0_y_xxzz_xxx, tr_0_0_y_xxzz_xxy, tr_0_0_y_xxzz_xxz, tr_0_0_y_xxzz_xyy, tr_0_0_y_xxzz_xyz, tr_0_0_y_xxzz_xzz, tr_0_0_y_xxzz_yyy, tr_0_0_y_xxzz_yyz, tr_0_0_y_xxzz_yzz, tr_0_0_y_xxzz_zzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, tr_xxzz_xx, tr_xxzz_xxxy, tr_xxzz_xxyy, tr_xxzz_xxyz, tr_xxzz_xy, tr_xxzz_xyyy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_yy, tr_xxzz_yyyy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxzz_xxx[i] = 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 2.0 * tr_xxzz_xxxy[i] * tke_0;

        tr_0_0_y_xxzz_xxy[i] = 2.0 * tr_xxyzz_xxy[i] * tbe_0 + 2.0 * tr_xxzz_xxyy[i] * tke_0 - tr_xxzz_xx[i];

        tr_0_0_y_xxzz_xxz[i] = 2.0 * tr_xxyzz_xxz[i] * tbe_0 + 2.0 * tr_xxzz_xxyz[i] * tke_0;

        tr_0_0_y_xxzz_xyy[i] = 2.0 * tr_xxyzz_xyy[i] * tbe_0 + 2.0 * tr_xxzz_xyyy[i] * tke_0 - 2.0 * tr_xxzz_xy[i];

        tr_0_0_y_xxzz_xyz[i] = 2.0 * tr_xxyzz_xyz[i] * tbe_0 + 2.0 * tr_xxzz_xyyz[i] * tke_0 - tr_xxzz_xz[i];

        tr_0_0_y_xxzz_xzz[i] = 2.0 * tr_xxyzz_xzz[i] * tbe_0 + 2.0 * tr_xxzz_xyzz[i] * tke_0;

        tr_0_0_y_xxzz_yyy[i] = 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 2.0 * tr_xxzz_yyyy[i] * tke_0 - 3.0 * tr_xxzz_yy[i];

        tr_0_0_y_xxzz_yyz[i] = 2.0 * tr_xxyzz_yyz[i] * tbe_0 + 2.0 * tr_xxzz_yyyz[i] * tke_0 - 2.0 * tr_xxzz_yz[i];

        tr_0_0_y_xxzz_yzz[i] = 2.0 * tr_xxyzz_yzz[i] * tbe_0 + 2.0 * tr_xxzz_yyzz[i] * tke_0 - tr_xxzz_zz[i];

        tr_0_0_y_xxzz_zzz[i] = 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 2.0 * tr_xxzz_yzzz[i] * tke_0;
    }

    // Set up 210-220 components of targeted buffer : GF

    auto tr_0_0_y_xyyy_xxx = pbuffer.data(idx_op_geom_010_gf + 210);

    auto tr_0_0_y_xyyy_xxy = pbuffer.data(idx_op_geom_010_gf + 211);

    auto tr_0_0_y_xyyy_xxz = pbuffer.data(idx_op_geom_010_gf + 212);

    auto tr_0_0_y_xyyy_xyy = pbuffer.data(idx_op_geom_010_gf + 213);

    auto tr_0_0_y_xyyy_xyz = pbuffer.data(idx_op_geom_010_gf + 214);

    auto tr_0_0_y_xyyy_xzz = pbuffer.data(idx_op_geom_010_gf + 215);

    auto tr_0_0_y_xyyy_yyy = pbuffer.data(idx_op_geom_010_gf + 216);

    auto tr_0_0_y_xyyy_yyz = pbuffer.data(idx_op_geom_010_gf + 217);

    auto tr_0_0_y_xyyy_yzz = pbuffer.data(idx_op_geom_010_gf + 218);

    auto tr_0_0_y_xyyy_zzz = pbuffer.data(idx_op_geom_010_gf + 219);

    #pragma omp simd aligned(tr_0_0_y_xyyy_xxx, tr_0_0_y_xyyy_xxy, tr_0_0_y_xyyy_xxz, tr_0_0_y_xyyy_xyy, tr_0_0_y_xyyy_xyz, tr_0_0_y_xyyy_xzz, tr_0_0_y_xyyy_yyy, tr_0_0_y_xyyy_yyz, tr_0_0_y_xyyy_yzz, tr_0_0_y_xyyy_zzz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyy_xx, tr_xyyy_xxxy, tr_xyyy_xxyy, tr_xyyy_xxyz, tr_xyyy_xy, tr_xyyy_xyyy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_yy, tr_xyyy_yyyy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyyy_xxx, tr_xyyyy_xxy, tr_xyyyy_xxz, tr_xyyyy_xyy, tr_xyyyy_xyz, tr_xyyyy_xzz, tr_xyyyy_yyy, tr_xyyyy_yyz, tr_xyyyy_yzz, tr_xyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyy_xxx[i] = 2.0 * tr_xyyyy_xxx[i] * tbe_0 + 2.0 * tr_xyyy_xxxy[i] * tke_0 - 3.0 * tr_xyy_xxx[i];

        tr_0_0_y_xyyy_xxy[i] = 2.0 * tr_xyyyy_xxy[i] * tbe_0 + 2.0 * tr_xyyy_xxyy[i] * tke_0 - 3.0 * tr_xyy_xxy[i] - tr_xyyy_xx[i];

        tr_0_0_y_xyyy_xxz[i] = 2.0 * tr_xyyyy_xxz[i] * tbe_0 + 2.0 * tr_xyyy_xxyz[i] * tke_0 - 3.0 * tr_xyy_xxz[i];

        tr_0_0_y_xyyy_xyy[i] = 2.0 * tr_xyyyy_xyy[i] * tbe_0 + 2.0 * tr_xyyy_xyyy[i] * tke_0 - 3.0 * tr_xyy_xyy[i] - 2.0 * tr_xyyy_xy[i];

        tr_0_0_y_xyyy_xyz[i] = 2.0 * tr_xyyyy_xyz[i] * tbe_0 + 2.0 * tr_xyyy_xyyz[i] * tke_0 - 3.0 * tr_xyy_xyz[i] - tr_xyyy_xz[i];

        tr_0_0_y_xyyy_xzz[i] = 2.0 * tr_xyyyy_xzz[i] * tbe_0 + 2.0 * tr_xyyy_xyzz[i] * tke_0 - 3.0 * tr_xyy_xzz[i];

        tr_0_0_y_xyyy_yyy[i] = 2.0 * tr_xyyyy_yyy[i] * tbe_0 + 2.0 * tr_xyyy_yyyy[i] * tke_0 - 3.0 * tr_xyy_yyy[i] - 3.0 * tr_xyyy_yy[i];

        tr_0_0_y_xyyy_yyz[i] = 2.0 * tr_xyyyy_yyz[i] * tbe_0 + 2.0 * tr_xyyy_yyyz[i] * tke_0 - 3.0 * tr_xyy_yyz[i] - 2.0 * tr_xyyy_yz[i];

        tr_0_0_y_xyyy_yzz[i] = 2.0 * tr_xyyyy_yzz[i] * tbe_0 + 2.0 * tr_xyyy_yyzz[i] * tke_0 - 3.0 * tr_xyy_yzz[i] - tr_xyyy_zz[i];

        tr_0_0_y_xyyy_zzz[i] = 2.0 * tr_xyyyy_zzz[i] * tbe_0 + 2.0 * tr_xyyy_yzzz[i] * tke_0 - 3.0 * tr_xyy_zzz[i];
    }

    // Set up 220-230 components of targeted buffer : GF

    auto tr_0_0_y_xyyz_xxx = pbuffer.data(idx_op_geom_010_gf + 220);

    auto tr_0_0_y_xyyz_xxy = pbuffer.data(idx_op_geom_010_gf + 221);

    auto tr_0_0_y_xyyz_xxz = pbuffer.data(idx_op_geom_010_gf + 222);

    auto tr_0_0_y_xyyz_xyy = pbuffer.data(idx_op_geom_010_gf + 223);

    auto tr_0_0_y_xyyz_xyz = pbuffer.data(idx_op_geom_010_gf + 224);

    auto tr_0_0_y_xyyz_xzz = pbuffer.data(idx_op_geom_010_gf + 225);

    auto tr_0_0_y_xyyz_yyy = pbuffer.data(idx_op_geom_010_gf + 226);

    auto tr_0_0_y_xyyz_yyz = pbuffer.data(idx_op_geom_010_gf + 227);

    auto tr_0_0_y_xyyz_yzz = pbuffer.data(idx_op_geom_010_gf + 228);

    auto tr_0_0_y_xyyz_zzz = pbuffer.data(idx_op_geom_010_gf + 229);

    #pragma omp simd aligned(tr_0_0_y_xyyz_xxx, tr_0_0_y_xyyz_xxy, tr_0_0_y_xyyz_xxz, tr_0_0_y_xyyz_xyy, tr_0_0_y_xyyz_xyz, tr_0_0_y_xyyz_xzz, tr_0_0_y_xyyz_yyy, tr_0_0_y_xyyz_yyz, tr_0_0_y_xyyz_yzz, tr_0_0_y_xyyz_zzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, tr_xyyz_xx, tr_xyyz_xxxy, tr_xyyz_xxyy, tr_xyyz_xxyz, tr_xyyz_xy, tr_xyyz_xyyy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_yy, tr_xyyz_yyyy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyz_xxx[i] = 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 2.0 * tr_xyyz_xxxy[i] * tke_0 - 2.0 * tr_xyz_xxx[i];

        tr_0_0_y_xyyz_xxy[i] = 2.0 * tr_xyyyz_xxy[i] * tbe_0 + 2.0 * tr_xyyz_xxyy[i] * tke_0 - 2.0 * tr_xyz_xxy[i] - tr_xyyz_xx[i];

        tr_0_0_y_xyyz_xxz[i] = 2.0 * tr_xyyyz_xxz[i] * tbe_0 + 2.0 * tr_xyyz_xxyz[i] * tke_0 - 2.0 * tr_xyz_xxz[i];

        tr_0_0_y_xyyz_xyy[i] = 2.0 * tr_xyyyz_xyy[i] * tbe_0 + 2.0 * tr_xyyz_xyyy[i] * tke_0 - 2.0 * tr_xyz_xyy[i] - 2.0 * tr_xyyz_xy[i];

        tr_0_0_y_xyyz_xyz[i] = 2.0 * tr_xyyyz_xyz[i] * tbe_0 + 2.0 * tr_xyyz_xyyz[i] * tke_0 - 2.0 * tr_xyz_xyz[i] - tr_xyyz_xz[i];

        tr_0_0_y_xyyz_xzz[i] = 2.0 * tr_xyyyz_xzz[i] * tbe_0 + 2.0 * tr_xyyz_xyzz[i] * tke_0 - 2.0 * tr_xyz_xzz[i];

        tr_0_0_y_xyyz_yyy[i] = 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 2.0 * tr_xyyz_yyyy[i] * tke_0 - 2.0 * tr_xyz_yyy[i] - 3.0 * tr_xyyz_yy[i];

        tr_0_0_y_xyyz_yyz[i] = 2.0 * tr_xyyyz_yyz[i] * tbe_0 + 2.0 * tr_xyyz_yyyz[i] * tke_0 - 2.0 * tr_xyz_yyz[i] - 2.0 * tr_xyyz_yz[i];

        tr_0_0_y_xyyz_yzz[i] = 2.0 * tr_xyyyz_yzz[i] * tbe_0 + 2.0 * tr_xyyz_yyzz[i] * tke_0 - 2.0 * tr_xyz_yzz[i] - tr_xyyz_zz[i];

        tr_0_0_y_xyyz_zzz[i] = 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 2.0 * tr_xyyz_yzzz[i] * tke_0 - 2.0 * tr_xyz_zzz[i];
    }

    // Set up 230-240 components of targeted buffer : GF

    auto tr_0_0_y_xyzz_xxx = pbuffer.data(idx_op_geom_010_gf + 230);

    auto tr_0_0_y_xyzz_xxy = pbuffer.data(idx_op_geom_010_gf + 231);

    auto tr_0_0_y_xyzz_xxz = pbuffer.data(idx_op_geom_010_gf + 232);

    auto tr_0_0_y_xyzz_xyy = pbuffer.data(idx_op_geom_010_gf + 233);

    auto tr_0_0_y_xyzz_xyz = pbuffer.data(idx_op_geom_010_gf + 234);

    auto tr_0_0_y_xyzz_xzz = pbuffer.data(idx_op_geom_010_gf + 235);

    auto tr_0_0_y_xyzz_yyy = pbuffer.data(idx_op_geom_010_gf + 236);

    auto tr_0_0_y_xyzz_yyz = pbuffer.data(idx_op_geom_010_gf + 237);

    auto tr_0_0_y_xyzz_yzz = pbuffer.data(idx_op_geom_010_gf + 238);

    auto tr_0_0_y_xyzz_zzz = pbuffer.data(idx_op_geom_010_gf + 239);

    #pragma omp simd aligned(tr_0_0_y_xyzz_xxx, tr_0_0_y_xyzz_xxy, tr_0_0_y_xyzz_xxz, tr_0_0_y_xyzz_xyy, tr_0_0_y_xyzz_xyz, tr_0_0_y_xyzz_xzz, tr_0_0_y_xyzz_yyy, tr_0_0_y_xyzz_yyz, tr_0_0_y_xyzz_yzz, tr_0_0_y_xyzz_zzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, tr_xyzz_xx, tr_xyzz_xxxy, tr_xyzz_xxyy, tr_xyzz_xxyz, tr_xyzz_xy, tr_xyzz_xyyy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_yy, tr_xyzz_yyyy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyzz_xxx[i] = 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 2.0 * tr_xyzz_xxxy[i] * tke_0 - tr_xzz_xxx[i];

        tr_0_0_y_xyzz_xxy[i] = 2.0 * tr_xyyzz_xxy[i] * tbe_0 + 2.0 * tr_xyzz_xxyy[i] * tke_0 - tr_xzz_xxy[i] - tr_xyzz_xx[i];

        tr_0_0_y_xyzz_xxz[i] = 2.0 * tr_xyyzz_xxz[i] * tbe_0 + 2.0 * tr_xyzz_xxyz[i] * tke_0 - tr_xzz_xxz[i];

        tr_0_0_y_xyzz_xyy[i] = 2.0 * tr_xyyzz_xyy[i] * tbe_0 + 2.0 * tr_xyzz_xyyy[i] * tke_0 - tr_xzz_xyy[i] - 2.0 * tr_xyzz_xy[i];

        tr_0_0_y_xyzz_xyz[i] = 2.0 * tr_xyyzz_xyz[i] * tbe_0 + 2.0 * tr_xyzz_xyyz[i] * tke_0 - tr_xzz_xyz[i] - tr_xyzz_xz[i];

        tr_0_0_y_xyzz_xzz[i] = 2.0 * tr_xyyzz_xzz[i] * tbe_0 + 2.0 * tr_xyzz_xyzz[i] * tke_0 - tr_xzz_xzz[i];

        tr_0_0_y_xyzz_yyy[i] = 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 2.0 * tr_xyzz_yyyy[i] * tke_0 - tr_xzz_yyy[i] - 3.0 * tr_xyzz_yy[i];

        tr_0_0_y_xyzz_yyz[i] = 2.0 * tr_xyyzz_yyz[i] * tbe_0 + 2.0 * tr_xyzz_yyyz[i] * tke_0 - tr_xzz_yyz[i] - 2.0 * tr_xyzz_yz[i];

        tr_0_0_y_xyzz_yzz[i] = 2.0 * tr_xyyzz_yzz[i] * tbe_0 + 2.0 * tr_xyzz_yyzz[i] * tke_0 - tr_xzz_yzz[i] - tr_xyzz_zz[i];

        tr_0_0_y_xyzz_zzz[i] = 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 2.0 * tr_xyzz_yzzz[i] * tke_0 - tr_xzz_zzz[i];
    }

    // Set up 240-250 components of targeted buffer : GF

    auto tr_0_0_y_xzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 240);

    auto tr_0_0_y_xzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 241);

    auto tr_0_0_y_xzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 242);

    auto tr_0_0_y_xzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 243);

    auto tr_0_0_y_xzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 244);

    auto tr_0_0_y_xzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 245);

    auto tr_0_0_y_xzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 246);

    auto tr_0_0_y_xzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 247);

    auto tr_0_0_y_xzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 248);

    auto tr_0_0_y_xzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 249);

    #pragma omp simd aligned(tr_0_0_y_xzzz_xxx, tr_0_0_y_xzzz_xxy, tr_0_0_y_xzzz_xxz, tr_0_0_y_xzzz_xyy, tr_0_0_y_xzzz_xyz, tr_0_0_y_xzzz_xzz, tr_0_0_y_xzzz_yyy, tr_0_0_y_xzzz_yyz, tr_0_0_y_xzzz_yzz, tr_0_0_y_xzzz_zzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, tr_xzzz_xx, tr_xzzz_xxxy, tr_xzzz_xxyy, tr_xzzz_xxyz, tr_xzzz_xy, tr_xzzz_xyyy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_yy, tr_xzzz_yyyy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzzz_xxx[i] = 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 2.0 * tr_xzzz_xxxy[i] * tke_0;

        tr_0_0_y_xzzz_xxy[i] = 2.0 * tr_xyzzz_xxy[i] * tbe_0 + 2.0 * tr_xzzz_xxyy[i] * tke_0 - tr_xzzz_xx[i];

        tr_0_0_y_xzzz_xxz[i] = 2.0 * tr_xyzzz_xxz[i] * tbe_0 + 2.0 * tr_xzzz_xxyz[i] * tke_0;

        tr_0_0_y_xzzz_xyy[i] = 2.0 * tr_xyzzz_xyy[i] * tbe_0 + 2.0 * tr_xzzz_xyyy[i] * tke_0 - 2.0 * tr_xzzz_xy[i];

        tr_0_0_y_xzzz_xyz[i] = 2.0 * tr_xyzzz_xyz[i] * tbe_0 + 2.0 * tr_xzzz_xyyz[i] * tke_0 - tr_xzzz_xz[i];

        tr_0_0_y_xzzz_xzz[i] = 2.0 * tr_xyzzz_xzz[i] * tbe_0 + 2.0 * tr_xzzz_xyzz[i] * tke_0;

        tr_0_0_y_xzzz_yyy[i] = 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 2.0 * tr_xzzz_yyyy[i] * tke_0 - 3.0 * tr_xzzz_yy[i];

        tr_0_0_y_xzzz_yyz[i] = 2.0 * tr_xyzzz_yyz[i] * tbe_0 + 2.0 * tr_xzzz_yyyz[i] * tke_0 - 2.0 * tr_xzzz_yz[i];

        tr_0_0_y_xzzz_yzz[i] = 2.0 * tr_xyzzz_yzz[i] * tbe_0 + 2.0 * tr_xzzz_yyzz[i] * tke_0 - tr_xzzz_zz[i];

        tr_0_0_y_xzzz_zzz[i] = 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 2.0 * tr_xzzz_yzzz[i] * tke_0;
    }

    // Set up 250-260 components of targeted buffer : GF

    auto tr_0_0_y_yyyy_xxx = pbuffer.data(idx_op_geom_010_gf + 250);

    auto tr_0_0_y_yyyy_xxy = pbuffer.data(idx_op_geom_010_gf + 251);

    auto tr_0_0_y_yyyy_xxz = pbuffer.data(idx_op_geom_010_gf + 252);

    auto tr_0_0_y_yyyy_xyy = pbuffer.data(idx_op_geom_010_gf + 253);

    auto tr_0_0_y_yyyy_xyz = pbuffer.data(idx_op_geom_010_gf + 254);

    auto tr_0_0_y_yyyy_xzz = pbuffer.data(idx_op_geom_010_gf + 255);

    auto tr_0_0_y_yyyy_yyy = pbuffer.data(idx_op_geom_010_gf + 256);

    auto tr_0_0_y_yyyy_yyz = pbuffer.data(idx_op_geom_010_gf + 257);

    auto tr_0_0_y_yyyy_yzz = pbuffer.data(idx_op_geom_010_gf + 258);

    auto tr_0_0_y_yyyy_zzz = pbuffer.data(idx_op_geom_010_gf + 259);

    #pragma omp simd aligned(tr_0_0_y_yyyy_xxx, tr_0_0_y_yyyy_xxy, tr_0_0_y_yyyy_xxz, tr_0_0_y_yyyy_xyy, tr_0_0_y_yyyy_xyz, tr_0_0_y_yyyy_xzz, tr_0_0_y_yyyy_yyy, tr_0_0_y_yyyy_yyz, tr_0_0_y_yyyy_yzz, tr_0_0_y_yyyy_zzz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, tr_yyyy_xx, tr_yyyy_xxxy, tr_yyyy_xxyy, tr_yyyy_xxyz, tr_yyyy_xy, tr_yyyy_xyyy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_yy, tr_yyyy_yyyy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyyy_xxx, tr_yyyyy_xxy, tr_yyyyy_xxz, tr_yyyyy_xyy, tr_yyyyy_xyz, tr_yyyyy_xzz, tr_yyyyy_yyy, tr_yyyyy_yyz, tr_yyyyy_yzz, tr_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyy_xxx[i] = 2.0 * tr_yyyyy_xxx[i] * tbe_0 + 2.0 * tr_yyyy_xxxy[i] * tke_0 - 4.0 * tr_yyy_xxx[i];

        tr_0_0_y_yyyy_xxy[i] = 2.0 * tr_yyyyy_xxy[i] * tbe_0 + 2.0 * tr_yyyy_xxyy[i] * tke_0 - 4.0 * tr_yyy_xxy[i] - tr_yyyy_xx[i];

        tr_0_0_y_yyyy_xxz[i] = 2.0 * tr_yyyyy_xxz[i] * tbe_0 + 2.0 * tr_yyyy_xxyz[i] * tke_0 - 4.0 * tr_yyy_xxz[i];

        tr_0_0_y_yyyy_xyy[i] = 2.0 * tr_yyyyy_xyy[i] * tbe_0 + 2.0 * tr_yyyy_xyyy[i] * tke_0 - 4.0 * tr_yyy_xyy[i] - 2.0 * tr_yyyy_xy[i];

        tr_0_0_y_yyyy_xyz[i] = 2.0 * tr_yyyyy_xyz[i] * tbe_0 + 2.0 * tr_yyyy_xyyz[i] * tke_0 - 4.0 * tr_yyy_xyz[i] - tr_yyyy_xz[i];

        tr_0_0_y_yyyy_xzz[i] = 2.0 * tr_yyyyy_xzz[i] * tbe_0 + 2.0 * tr_yyyy_xyzz[i] * tke_0 - 4.0 * tr_yyy_xzz[i];

        tr_0_0_y_yyyy_yyy[i] = 2.0 * tr_yyyyy_yyy[i] * tbe_0 + 2.0 * tr_yyyy_yyyy[i] * tke_0 - 4.0 * tr_yyy_yyy[i] - 3.0 * tr_yyyy_yy[i];

        tr_0_0_y_yyyy_yyz[i] = 2.0 * tr_yyyyy_yyz[i] * tbe_0 + 2.0 * tr_yyyy_yyyz[i] * tke_0 - 4.0 * tr_yyy_yyz[i] - 2.0 * tr_yyyy_yz[i];

        tr_0_0_y_yyyy_yzz[i] = 2.0 * tr_yyyyy_yzz[i] * tbe_0 + 2.0 * tr_yyyy_yyzz[i] * tke_0 - 4.0 * tr_yyy_yzz[i] - tr_yyyy_zz[i];

        tr_0_0_y_yyyy_zzz[i] = 2.0 * tr_yyyyy_zzz[i] * tbe_0 + 2.0 * tr_yyyy_yzzz[i] * tke_0 - 4.0 * tr_yyy_zzz[i];
    }

    // Set up 260-270 components of targeted buffer : GF

    auto tr_0_0_y_yyyz_xxx = pbuffer.data(idx_op_geom_010_gf + 260);

    auto tr_0_0_y_yyyz_xxy = pbuffer.data(idx_op_geom_010_gf + 261);

    auto tr_0_0_y_yyyz_xxz = pbuffer.data(idx_op_geom_010_gf + 262);

    auto tr_0_0_y_yyyz_xyy = pbuffer.data(idx_op_geom_010_gf + 263);

    auto tr_0_0_y_yyyz_xyz = pbuffer.data(idx_op_geom_010_gf + 264);

    auto tr_0_0_y_yyyz_xzz = pbuffer.data(idx_op_geom_010_gf + 265);

    auto tr_0_0_y_yyyz_yyy = pbuffer.data(idx_op_geom_010_gf + 266);

    auto tr_0_0_y_yyyz_yyz = pbuffer.data(idx_op_geom_010_gf + 267);

    auto tr_0_0_y_yyyz_yzz = pbuffer.data(idx_op_geom_010_gf + 268);

    auto tr_0_0_y_yyyz_zzz = pbuffer.data(idx_op_geom_010_gf + 269);

    #pragma omp simd aligned(tr_0_0_y_yyyz_xxx, tr_0_0_y_yyyz_xxy, tr_0_0_y_yyyz_xxz, tr_0_0_y_yyyz_xyy, tr_0_0_y_yyyz_xyz, tr_0_0_y_yyyz_xzz, tr_0_0_y_yyyz_yyy, tr_0_0_y_yyyz_yyz, tr_0_0_y_yyyz_yzz, tr_0_0_y_yyyz_zzz, tr_yyyyz_xxx, tr_yyyyz_xxy, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_zzz, tr_yyyz_xx, tr_yyyz_xxxy, tr_yyyz_xxyy, tr_yyyz_xxyz, tr_yyyz_xy, tr_yyyz_xyyy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_yy, tr_yyyz_yyyy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyz_xxx[i] = 2.0 * tr_yyyyz_xxx[i] * tbe_0 + 2.0 * tr_yyyz_xxxy[i] * tke_0 - 3.0 * tr_yyz_xxx[i];

        tr_0_0_y_yyyz_xxy[i] = 2.0 * tr_yyyyz_xxy[i] * tbe_0 + 2.0 * tr_yyyz_xxyy[i] * tke_0 - 3.0 * tr_yyz_xxy[i] - tr_yyyz_xx[i];

        tr_0_0_y_yyyz_xxz[i] = 2.0 * tr_yyyyz_xxz[i] * tbe_0 + 2.0 * tr_yyyz_xxyz[i] * tke_0 - 3.0 * tr_yyz_xxz[i];

        tr_0_0_y_yyyz_xyy[i] = 2.0 * tr_yyyyz_xyy[i] * tbe_0 + 2.0 * tr_yyyz_xyyy[i] * tke_0 - 3.0 * tr_yyz_xyy[i] - 2.0 * tr_yyyz_xy[i];

        tr_0_0_y_yyyz_xyz[i] = 2.0 * tr_yyyyz_xyz[i] * tbe_0 + 2.0 * tr_yyyz_xyyz[i] * tke_0 - 3.0 * tr_yyz_xyz[i] - tr_yyyz_xz[i];

        tr_0_0_y_yyyz_xzz[i] = 2.0 * tr_yyyyz_xzz[i] * tbe_0 + 2.0 * tr_yyyz_xyzz[i] * tke_0 - 3.0 * tr_yyz_xzz[i];

        tr_0_0_y_yyyz_yyy[i] = 2.0 * tr_yyyyz_yyy[i] * tbe_0 + 2.0 * tr_yyyz_yyyy[i] * tke_0 - 3.0 * tr_yyz_yyy[i] - 3.0 * tr_yyyz_yy[i];

        tr_0_0_y_yyyz_yyz[i] = 2.0 * tr_yyyyz_yyz[i] * tbe_0 + 2.0 * tr_yyyz_yyyz[i] * tke_0 - 3.0 * tr_yyz_yyz[i] - 2.0 * tr_yyyz_yz[i];

        tr_0_0_y_yyyz_yzz[i] = 2.0 * tr_yyyyz_yzz[i] * tbe_0 + 2.0 * tr_yyyz_yyzz[i] * tke_0 - 3.0 * tr_yyz_yzz[i] - tr_yyyz_zz[i];

        tr_0_0_y_yyyz_zzz[i] = 2.0 * tr_yyyyz_zzz[i] * tbe_0 + 2.0 * tr_yyyz_yzzz[i] * tke_0 - 3.0 * tr_yyz_zzz[i];
    }

    // Set up 270-280 components of targeted buffer : GF

    auto tr_0_0_y_yyzz_xxx = pbuffer.data(idx_op_geom_010_gf + 270);

    auto tr_0_0_y_yyzz_xxy = pbuffer.data(idx_op_geom_010_gf + 271);

    auto tr_0_0_y_yyzz_xxz = pbuffer.data(idx_op_geom_010_gf + 272);

    auto tr_0_0_y_yyzz_xyy = pbuffer.data(idx_op_geom_010_gf + 273);

    auto tr_0_0_y_yyzz_xyz = pbuffer.data(idx_op_geom_010_gf + 274);

    auto tr_0_0_y_yyzz_xzz = pbuffer.data(idx_op_geom_010_gf + 275);

    auto tr_0_0_y_yyzz_yyy = pbuffer.data(idx_op_geom_010_gf + 276);

    auto tr_0_0_y_yyzz_yyz = pbuffer.data(idx_op_geom_010_gf + 277);

    auto tr_0_0_y_yyzz_yzz = pbuffer.data(idx_op_geom_010_gf + 278);

    auto tr_0_0_y_yyzz_zzz = pbuffer.data(idx_op_geom_010_gf + 279);

    #pragma omp simd aligned(tr_0_0_y_yyzz_xxx, tr_0_0_y_yyzz_xxy, tr_0_0_y_yyzz_xxz, tr_0_0_y_yyzz_xyy, tr_0_0_y_yyzz_xyz, tr_0_0_y_yyzz_xzz, tr_0_0_y_yyzz_yyy, tr_0_0_y_yyzz_yyz, tr_0_0_y_yyzz_yzz, tr_0_0_y_yyzz_zzz, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, tr_yyzz_xx, tr_yyzz_xxxy, tr_yyzz_xxyy, tr_yyzz_xxyz, tr_yyzz_xy, tr_yyzz_xyyy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_yy, tr_yyzz_yyyy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyzz_xxx[i] = 2.0 * tr_yyyzz_xxx[i] * tbe_0 + 2.0 * tr_yyzz_xxxy[i] * tke_0 - 2.0 * tr_yzz_xxx[i];

        tr_0_0_y_yyzz_xxy[i] = 2.0 * tr_yyyzz_xxy[i] * tbe_0 + 2.0 * tr_yyzz_xxyy[i] * tke_0 - 2.0 * tr_yzz_xxy[i] - tr_yyzz_xx[i];

        tr_0_0_y_yyzz_xxz[i] = 2.0 * tr_yyyzz_xxz[i] * tbe_0 + 2.0 * tr_yyzz_xxyz[i] * tke_0 - 2.0 * tr_yzz_xxz[i];

        tr_0_0_y_yyzz_xyy[i] = 2.0 * tr_yyyzz_xyy[i] * tbe_0 + 2.0 * tr_yyzz_xyyy[i] * tke_0 - 2.0 * tr_yzz_xyy[i] - 2.0 * tr_yyzz_xy[i];

        tr_0_0_y_yyzz_xyz[i] = 2.0 * tr_yyyzz_xyz[i] * tbe_0 + 2.0 * tr_yyzz_xyyz[i] * tke_0 - 2.0 * tr_yzz_xyz[i] - tr_yyzz_xz[i];

        tr_0_0_y_yyzz_xzz[i] = 2.0 * tr_yyyzz_xzz[i] * tbe_0 + 2.0 * tr_yyzz_xyzz[i] * tke_0 - 2.0 * tr_yzz_xzz[i];

        tr_0_0_y_yyzz_yyy[i] = 2.0 * tr_yyyzz_yyy[i] * tbe_0 + 2.0 * tr_yyzz_yyyy[i] * tke_0 - 2.0 * tr_yzz_yyy[i] - 3.0 * tr_yyzz_yy[i];

        tr_0_0_y_yyzz_yyz[i] = 2.0 * tr_yyyzz_yyz[i] * tbe_0 + 2.0 * tr_yyzz_yyyz[i] * tke_0 - 2.0 * tr_yzz_yyz[i] - 2.0 * tr_yyzz_yz[i];

        tr_0_0_y_yyzz_yzz[i] = 2.0 * tr_yyyzz_yzz[i] * tbe_0 + 2.0 * tr_yyzz_yyzz[i] * tke_0 - 2.0 * tr_yzz_yzz[i] - tr_yyzz_zz[i];

        tr_0_0_y_yyzz_zzz[i] = 2.0 * tr_yyyzz_zzz[i] * tbe_0 + 2.0 * tr_yyzz_yzzz[i] * tke_0 - 2.0 * tr_yzz_zzz[i];
    }

    // Set up 280-290 components of targeted buffer : GF

    auto tr_0_0_y_yzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 280);

    auto tr_0_0_y_yzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 281);

    auto tr_0_0_y_yzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 282);

    auto tr_0_0_y_yzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 283);

    auto tr_0_0_y_yzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 284);

    auto tr_0_0_y_yzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 285);

    auto tr_0_0_y_yzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 286);

    auto tr_0_0_y_yzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 287);

    auto tr_0_0_y_yzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 288);

    auto tr_0_0_y_yzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 289);

    #pragma omp simd aligned(tr_0_0_y_yzzz_xxx, tr_0_0_y_yzzz_xxy, tr_0_0_y_yzzz_xxz, tr_0_0_y_yzzz_xyy, tr_0_0_y_yzzz_xyz, tr_0_0_y_yzzz_xzz, tr_0_0_y_yzzz_yyy, tr_0_0_y_yzzz_yyz, tr_0_0_y_yzzz_yzz, tr_0_0_y_yzzz_zzz, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, tr_yzzz_xx, tr_yzzz_xxxy, tr_yzzz_xxyy, tr_yzzz_xxyz, tr_yzzz_xy, tr_yzzz_xyyy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_yy, tr_yzzz_yyyy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzzz_xxx[i] = 2.0 * tr_yyzzz_xxx[i] * tbe_0 + 2.0 * tr_yzzz_xxxy[i] * tke_0 - tr_zzz_xxx[i];

        tr_0_0_y_yzzz_xxy[i] = 2.0 * tr_yyzzz_xxy[i] * tbe_0 + 2.0 * tr_yzzz_xxyy[i] * tke_0 - tr_zzz_xxy[i] - tr_yzzz_xx[i];

        tr_0_0_y_yzzz_xxz[i] = 2.0 * tr_yyzzz_xxz[i] * tbe_0 + 2.0 * tr_yzzz_xxyz[i] * tke_0 - tr_zzz_xxz[i];

        tr_0_0_y_yzzz_xyy[i] = 2.0 * tr_yyzzz_xyy[i] * tbe_0 + 2.0 * tr_yzzz_xyyy[i] * tke_0 - tr_zzz_xyy[i] - 2.0 * tr_yzzz_xy[i];

        tr_0_0_y_yzzz_xyz[i] = 2.0 * tr_yyzzz_xyz[i] * tbe_0 + 2.0 * tr_yzzz_xyyz[i] * tke_0 - tr_zzz_xyz[i] - tr_yzzz_xz[i];

        tr_0_0_y_yzzz_xzz[i] = 2.0 * tr_yyzzz_xzz[i] * tbe_0 + 2.0 * tr_yzzz_xyzz[i] * tke_0 - tr_zzz_xzz[i];

        tr_0_0_y_yzzz_yyy[i] = 2.0 * tr_yyzzz_yyy[i] * tbe_0 + 2.0 * tr_yzzz_yyyy[i] * tke_0 - tr_zzz_yyy[i] - 3.0 * tr_yzzz_yy[i];

        tr_0_0_y_yzzz_yyz[i] = 2.0 * tr_yyzzz_yyz[i] * tbe_0 + 2.0 * tr_yzzz_yyyz[i] * tke_0 - tr_zzz_yyz[i] - 2.0 * tr_yzzz_yz[i];

        tr_0_0_y_yzzz_yzz[i] = 2.0 * tr_yyzzz_yzz[i] * tbe_0 + 2.0 * tr_yzzz_yyzz[i] * tke_0 - tr_zzz_yzz[i] - tr_yzzz_zz[i];

        tr_0_0_y_yzzz_zzz[i] = 2.0 * tr_yyzzz_zzz[i] * tbe_0 + 2.0 * tr_yzzz_yzzz[i] * tke_0 - tr_zzz_zzz[i];
    }

    // Set up 290-300 components of targeted buffer : GF

    auto tr_0_0_y_zzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 290);

    auto tr_0_0_y_zzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 291);

    auto tr_0_0_y_zzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 292);

    auto tr_0_0_y_zzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 293);

    auto tr_0_0_y_zzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 294);

    auto tr_0_0_y_zzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 295);

    auto tr_0_0_y_zzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 296);

    auto tr_0_0_y_zzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 297);

    auto tr_0_0_y_zzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 298);

    auto tr_0_0_y_zzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 299);

    #pragma omp simd aligned(tr_0_0_y_zzzz_xxx, tr_0_0_y_zzzz_xxy, tr_0_0_y_zzzz_xxz, tr_0_0_y_zzzz_xyy, tr_0_0_y_zzzz_xyz, tr_0_0_y_zzzz_xzz, tr_0_0_y_zzzz_yyy, tr_0_0_y_zzzz_yyz, tr_0_0_y_zzzz_yzz, tr_0_0_y_zzzz_zzz, tr_yzzzz_xxx, tr_yzzzz_xxy, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_zzz, tr_zzzz_xx, tr_zzzz_xxxy, tr_zzzz_xxyy, tr_zzzz_xxyz, tr_zzzz_xy, tr_zzzz_xyyy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_yy, tr_zzzz_yyyy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzzz_xxx[i] = 2.0 * tr_yzzzz_xxx[i] * tbe_0 + 2.0 * tr_zzzz_xxxy[i] * tke_0;

        tr_0_0_y_zzzz_xxy[i] = 2.0 * tr_yzzzz_xxy[i] * tbe_0 + 2.0 * tr_zzzz_xxyy[i] * tke_0 - tr_zzzz_xx[i];

        tr_0_0_y_zzzz_xxz[i] = 2.0 * tr_yzzzz_xxz[i] * tbe_0 + 2.0 * tr_zzzz_xxyz[i] * tke_0;

        tr_0_0_y_zzzz_xyy[i] = 2.0 * tr_yzzzz_xyy[i] * tbe_0 + 2.0 * tr_zzzz_xyyy[i] * tke_0 - 2.0 * tr_zzzz_xy[i];

        tr_0_0_y_zzzz_xyz[i] = 2.0 * tr_yzzzz_xyz[i] * tbe_0 + 2.0 * tr_zzzz_xyyz[i] * tke_0 - tr_zzzz_xz[i];

        tr_0_0_y_zzzz_xzz[i] = 2.0 * tr_yzzzz_xzz[i] * tbe_0 + 2.0 * tr_zzzz_xyzz[i] * tke_0;

        tr_0_0_y_zzzz_yyy[i] = 2.0 * tr_yzzzz_yyy[i] * tbe_0 + 2.0 * tr_zzzz_yyyy[i] * tke_0 - 3.0 * tr_zzzz_yy[i];

        tr_0_0_y_zzzz_yyz[i] = 2.0 * tr_yzzzz_yyz[i] * tbe_0 + 2.0 * tr_zzzz_yyyz[i] * tke_0 - 2.0 * tr_zzzz_yz[i];

        tr_0_0_y_zzzz_yzz[i] = 2.0 * tr_yzzzz_yzz[i] * tbe_0 + 2.0 * tr_zzzz_yyzz[i] * tke_0 - tr_zzzz_zz[i];

        tr_0_0_y_zzzz_zzz[i] = 2.0 * tr_yzzzz_zzz[i] * tbe_0 + 2.0 * tr_zzzz_yzzz[i] * tke_0;
    }

    // Set up 300-310 components of targeted buffer : GF

    auto tr_0_0_z_xxxx_xxx = pbuffer.data(idx_op_geom_010_gf + 300);

    auto tr_0_0_z_xxxx_xxy = pbuffer.data(idx_op_geom_010_gf + 301);

    auto tr_0_0_z_xxxx_xxz = pbuffer.data(idx_op_geom_010_gf + 302);

    auto tr_0_0_z_xxxx_xyy = pbuffer.data(idx_op_geom_010_gf + 303);

    auto tr_0_0_z_xxxx_xyz = pbuffer.data(idx_op_geom_010_gf + 304);

    auto tr_0_0_z_xxxx_xzz = pbuffer.data(idx_op_geom_010_gf + 305);

    auto tr_0_0_z_xxxx_yyy = pbuffer.data(idx_op_geom_010_gf + 306);

    auto tr_0_0_z_xxxx_yyz = pbuffer.data(idx_op_geom_010_gf + 307);

    auto tr_0_0_z_xxxx_yzz = pbuffer.data(idx_op_geom_010_gf + 308);

    auto tr_0_0_z_xxxx_zzz = pbuffer.data(idx_op_geom_010_gf + 309);

    #pragma omp simd aligned(tr_0_0_z_xxxx_xxx, tr_0_0_z_xxxx_xxy, tr_0_0_z_xxxx_xxz, tr_0_0_z_xxxx_xyy, tr_0_0_z_xxxx_xyz, tr_0_0_z_xxxx_xzz, tr_0_0_z_xxxx_yyy, tr_0_0_z_xxxx_yyz, tr_0_0_z_xxxx_yzz, tr_0_0_z_xxxx_zzz, tr_xxxx_xx, tr_xxxx_xxxz, tr_xxxx_xxyz, tr_xxxx_xxzz, tr_xxxx_xy, tr_xxxx_xyyz, tr_xxxx_xyzz, tr_xxxx_xz, tr_xxxx_xzzz, tr_xxxx_yy, tr_xxxx_yyyz, tr_xxxx_yyzz, tr_xxxx_yz, tr_xxxx_yzzz, tr_xxxx_zz, tr_xxxx_zzzz, tr_xxxxz_xxx, tr_xxxxz_xxy, tr_xxxxz_xxz, tr_xxxxz_xyy, tr_xxxxz_xyz, tr_xxxxz_xzz, tr_xxxxz_yyy, tr_xxxxz_yyz, tr_xxxxz_yzz, tr_xxxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxx_xxx[i] = 2.0 * tr_xxxxz_xxx[i] * tbe_0 + 2.0 * tr_xxxx_xxxz[i] * tke_0;

        tr_0_0_z_xxxx_xxy[i] = 2.0 * tr_xxxxz_xxy[i] * tbe_0 + 2.0 * tr_xxxx_xxyz[i] * tke_0;

        tr_0_0_z_xxxx_xxz[i] = 2.0 * tr_xxxxz_xxz[i] * tbe_0 + 2.0 * tr_xxxx_xxzz[i] * tke_0 - tr_xxxx_xx[i];

        tr_0_0_z_xxxx_xyy[i] = 2.0 * tr_xxxxz_xyy[i] * tbe_0 + 2.0 * tr_xxxx_xyyz[i] * tke_0;

        tr_0_0_z_xxxx_xyz[i] = 2.0 * tr_xxxxz_xyz[i] * tbe_0 + 2.0 * tr_xxxx_xyzz[i] * tke_0 - tr_xxxx_xy[i];

        tr_0_0_z_xxxx_xzz[i] = 2.0 * tr_xxxxz_xzz[i] * tbe_0 + 2.0 * tr_xxxx_xzzz[i] * tke_0 - 2.0 * tr_xxxx_xz[i];

        tr_0_0_z_xxxx_yyy[i] = 2.0 * tr_xxxxz_yyy[i] * tbe_0 + 2.0 * tr_xxxx_yyyz[i] * tke_0;

        tr_0_0_z_xxxx_yyz[i] = 2.0 * tr_xxxxz_yyz[i] * tbe_0 + 2.0 * tr_xxxx_yyzz[i] * tke_0 - tr_xxxx_yy[i];

        tr_0_0_z_xxxx_yzz[i] = 2.0 * tr_xxxxz_yzz[i] * tbe_0 + 2.0 * tr_xxxx_yzzz[i] * tke_0 - 2.0 * tr_xxxx_yz[i];

        tr_0_0_z_xxxx_zzz[i] = 2.0 * tr_xxxxz_zzz[i] * tbe_0 + 2.0 * tr_xxxx_zzzz[i] * tke_0 - 3.0 * tr_xxxx_zz[i];
    }

    // Set up 310-320 components of targeted buffer : GF

    auto tr_0_0_z_xxxy_xxx = pbuffer.data(idx_op_geom_010_gf + 310);

    auto tr_0_0_z_xxxy_xxy = pbuffer.data(idx_op_geom_010_gf + 311);

    auto tr_0_0_z_xxxy_xxz = pbuffer.data(idx_op_geom_010_gf + 312);

    auto tr_0_0_z_xxxy_xyy = pbuffer.data(idx_op_geom_010_gf + 313);

    auto tr_0_0_z_xxxy_xyz = pbuffer.data(idx_op_geom_010_gf + 314);

    auto tr_0_0_z_xxxy_xzz = pbuffer.data(idx_op_geom_010_gf + 315);

    auto tr_0_0_z_xxxy_yyy = pbuffer.data(idx_op_geom_010_gf + 316);

    auto tr_0_0_z_xxxy_yyz = pbuffer.data(idx_op_geom_010_gf + 317);

    auto tr_0_0_z_xxxy_yzz = pbuffer.data(idx_op_geom_010_gf + 318);

    auto tr_0_0_z_xxxy_zzz = pbuffer.data(idx_op_geom_010_gf + 319);

    #pragma omp simd aligned(tr_0_0_z_xxxy_xxx, tr_0_0_z_xxxy_xxy, tr_0_0_z_xxxy_xxz, tr_0_0_z_xxxy_xyy, tr_0_0_z_xxxy_xyz, tr_0_0_z_xxxy_xzz, tr_0_0_z_xxxy_yyy, tr_0_0_z_xxxy_yyz, tr_0_0_z_xxxy_yzz, tr_0_0_z_xxxy_zzz, tr_xxxy_xx, tr_xxxy_xxxz, tr_xxxy_xxyz, tr_xxxy_xxzz, tr_xxxy_xy, tr_xxxy_xyyz, tr_xxxy_xyzz, tr_xxxy_xz, tr_xxxy_xzzz, tr_xxxy_yy, tr_xxxy_yyyz, tr_xxxy_yyzz, tr_xxxy_yz, tr_xxxy_yzzz, tr_xxxy_zz, tr_xxxy_zzzz, tr_xxxyz_xxx, tr_xxxyz_xxy, tr_xxxyz_xxz, tr_xxxyz_xyy, tr_xxxyz_xyz, tr_xxxyz_xzz, tr_xxxyz_yyy, tr_xxxyz_yyz, tr_xxxyz_yzz, tr_xxxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxy_xxx[i] = 2.0 * tr_xxxyz_xxx[i] * tbe_0 + 2.0 * tr_xxxy_xxxz[i] * tke_0;

        tr_0_0_z_xxxy_xxy[i] = 2.0 * tr_xxxyz_xxy[i] * tbe_0 + 2.0 * tr_xxxy_xxyz[i] * tke_0;

        tr_0_0_z_xxxy_xxz[i] = 2.0 * tr_xxxyz_xxz[i] * tbe_0 + 2.0 * tr_xxxy_xxzz[i] * tke_0 - tr_xxxy_xx[i];

        tr_0_0_z_xxxy_xyy[i] = 2.0 * tr_xxxyz_xyy[i] * tbe_0 + 2.0 * tr_xxxy_xyyz[i] * tke_0;

        tr_0_0_z_xxxy_xyz[i] = 2.0 * tr_xxxyz_xyz[i] * tbe_0 + 2.0 * tr_xxxy_xyzz[i] * tke_0 - tr_xxxy_xy[i];

        tr_0_0_z_xxxy_xzz[i] = 2.0 * tr_xxxyz_xzz[i] * tbe_0 + 2.0 * tr_xxxy_xzzz[i] * tke_0 - 2.0 * tr_xxxy_xz[i];

        tr_0_0_z_xxxy_yyy[i] = 2.0 * tr_xxxyz_yyy[i] * tbe_0 + 2.0 * tr_xxxy_yyyz[i] * tke_0;

        tr_0_0_z_xxxy_yyz[i] = 2.0 * tr_xxxyz_yyz[i] * tbe_0 + 2.0 * tr_xxxy_yyzz[i] * tke_0 - tr_xxxy_yy[i];

        tr_0_0_z_xxxy_yzz[i] = 2.0 * tr_xxxyz_yzz[i] * tbe_0 + 2.0 * tr_xxxy_yzzz[i] * tke_0 - 2.0 * tr_xxxy_yz[i];

        tr_0_0_z_xxxy_zzz[i] = 2.0 * tr_xxxyz_zzz[i] * tbe_0 + 2.0 * tr_xxxy_zzzz[i] * tke_0 - 3.0 * tr_xxxy_zz[i];
    }

    // Set up 320-330 components of targeted buffer : GF

    auto tr_0_0_z_xxxz_xxx = pbuffer.data(idx_op_geom_010_gf + 320);

    auto tr_0_0_z_xxxz_xxy = pbuffer.data(idx_op_geom_010_gf + 321);

    auto tr_0_0_z_xxxz_xxz = pbuffer.data(idx_op_geom_010_gf + 322);

    auto tr_0_0_z_xxxz_xyy = pbuffer.data(idx_op_geom_010_gf + 323);

    auto tr_0_0_z_xxxz_xyz = pbuffer.data(idx_op_geom_010_gf + 324);

    auto tr_0_0_z_xxxz_xzz = pbuffer.data(idx_op_geom_010_gf + 325);

    auto tr_0_0_z_xxxz_yyy = pbuffer.data(idx_op_geom_010_gf + 326);

    auto tr_0_0_z_xxxz_yyz = pbuffer.data(idx_op_geom_010_gf + 327);

    auto tr_0_0_z_xxxz_yzz = pbuffer.data(idx_op_geom_010_gf + 328);

    auto tr_0_0_z_xxxz_zzz = pbuffer.data(idx_op_geom_010_gf + 329);

    #pragma omp simd aligned(tr_0_0_z_xxxz_xxx, tr_0_0_z_xxxz_xxy, tr_0_0_z_xxxz_xxz, tr_0_0_z_xxxz_xyy, tr_0_0_z_xxxz_xyz, tr_0_0_z_xxxz_xzz, tr_0_0_z_xxxz_yyy, tr_0_0_z_xxxz_yyz, tr_0_0_z_xxxz_yzz, tr_0_0_z_xxxz_zzz, tr_xxx_xxx, tr_xxx_xxy, tr_xxx_xxz, tr_xxx_xyy, tr_xxx_xyz, tr_xxx_xzz, tr_xxx_yyy, tr_xxx_yyz, tr_xxx_yzz, tr_xxx_zzz, tr_xxxz_xx, tr_xxxz_xxxz, tr_xxxz_xxyz, tr_xxxz_xxzz, tr_xxxz_xy, tr_xxxz_xyyz, tr_xxxz_xyzz, tr_xxxz_xz, tr_xxxz_xzzz, tr_xxxz_yy, tr_xxxz_yyyz, tr_xxxz_yyzz, tr_xxxz_yz, tr_xxxz_yzzz, tr_xxxz_zz, tr_xxxz_zzzz, tr_xxxzz_xxx, tr_xxxzz_xxy, tr_xxxzz_xxz, tr_xxxzz_xyy, tr_xxxzz_xyz, tr_xxxzz_xzz, tr_xxxzz_yyy, tr_xxxzz_yyz, tr_xxxzz_yzz, tr_xxxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxz_xxx[i] = 2.0 * tr_xxxzz_xxx[i] * tbe_0 + 2.0 * tr_xxxz_xxxz[i] * tke_0 - tr_xxx_xxx[i];

        tr_0_0_z_xxxz_xxy[i] = 2.0 * tr_xxxzz_xxy[i] * tbe_0 + 2.0 * tr_xxxz_xxyz[i] * tke_0 - tr_xxx_xxy[i];

        tr_0_0_z_xxxz_xxz[i] = 2.0 * tr_xxxzz_xxz[i] * tbe_0 + 2.0 * tr_xxxz_xxzz[i] * tke_0 - tr_xxx_xxz[i] - tr_xxxz_xx[i];

        tr_0_0_z_xxxz_xyy[i] = 2.0 * tr_xxxzz_xyy[i] * tbe_0 + 2.0 * tr_xxxz_xyyz[i] * tke_0 - tr_xxx_xyy[i];

        tr_0_0_z_xxxz_xyz[i] = 2.0 * tr_xxxzz_xyz[i] * tbe_0 + 2.0 * tr_xxxz_xyzz[i] * tke_0 - tr_xxx_xyz[i] - tr_xxxz_xy[i];

        tr_0_0_z_xxxz_xzz[i] = 2.0 * tr_xxxzz_xzz[i] * tbe_0 + 2.0 * tr_xxxz_xzzz[i] * tke_0 - tr_xxx_xzz[i] - 2.0 * tr_xxxz_xz[i];

        tr_0_0_z_xxxz_yyy[i] = 2.0 * tr_xxxzz_yyy[i] * tbe_0 + 2.0 * tr_xxxz_yyyz[i] * tke_0 - tr_xxx_yyy[i];

        tr_0_0_z_xxxz_yyz[i] = 2.0 * tr_xxxzz_yyz[i] * tbe_0 + 2.0 * tr_xxxz_yyzz[i] * tke_0 - tr_xxx_yyz[i] - tr_xxxz_yy[i];

        tr_0_0_z_xxxz_yzz[i] = 2.0 * tr_xxxzz_yzz[i] * tbe_0 + 2.0 * tr_xxxz_yzzz[i] * tke_0 - tr_xxx_yzz[i] - 2.0 * tr_xxxz_yz[i];

        tr_0_0_z_xxxz_zzz[i] = 2.0 * tr_xxxzz_zzz[i] * tbe_0 + 2.0 * tr_xxxz_zzzz[i] * tke_0 - tr_xxx_zzz[i] - 3.0 * tr_xxxz_zz[i];
    }

    // Set up 330-340 components of targeted buffer : GF

    auto tr_0_0_z_xxyy_xxx = pbuffer.data(idx_op_geom_010_gf + 330);

    auto tr_0_0_z_xxyy_xxy = pbuffer.data(idx_op_geom_010_gf + 331);

    auto tr_0_0_z_xxyy_xxz = pbuffer.data(idx_op_geom_010_gf + 332);

    auto tr_0_0_z_xxyy_xyy = pbuffer.data(idx_op_geom_010_gf + 333);

    auto tr_0_0_z_xxyy_xyz = pbuffer.data(idx_op_geom_010_gf + 334);

    auto tr_0_0_z_xxyy_xzz = pbuffer.data(idx_op_geom_010_gf + 335);

    auto tr_0_0_z_xxyy_yyy = pbuffer.data(idx_op_geom_010_gf + 336);

    auto tr_0_0_z_xxyy_yyz = pbuffer.data(idx_op_geom_010_gf + 337);

    auto tr_0_0_z_xxyy_yzz = pbuffer.data(idx_op_geom_010_gf + 338);

    auto tr_0_0_z_xxyy_zzz = pbuffer.data(idx_op_geom_010_gf + 339);

    #pragma omp simd aligned(tr_0_0_z_xxyy_xxx, tr_0_0_z_xxyy_xxy, tr_0_0_z_xxyy_xxz, tr_0_0_z_xxyy_xyy, tr_0_0_z_xxyy_xyz, tr_0_0_z_xxyy_xzz, tr_0_0_z_xxyy_yyy, tr_0_0_z_xxyy_yyz, tr_0_0_z_xxyy_yzz, tr_0_0_z_xxyy_zzz, tr_xxyy_xx, tr_xxyy_xxxz, tr_xxyy_xxyz, tr_xxyy_xxzz, tr_xxyy_xy, tr_xxyy_xyyz, tr_xxyy_xyzz, tr_xxyy_xz, tr_xxyy_xzzz, tr_xxyy_yy, tr_xxyy_yyyz, tr_xxyy_yyzz, tr_xxyy_yz, tr_xxyy_yzzz, tr_xxyy_zz, tr_xxyy_zzzz, tr_xxyyz_xxx, tr_xxyyz_xxy, tr_xxyyz_xxz, tr_xxyyz_xyy, tr_xxyyz_xyz, tr_xxyyz_xzz, tr_xxyyz_yyy, tr_xxyyz_yyz, tr_xxyyz_yzz, tr_xxyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyy_xxx[i] = 2.0 * tr_xxyyz_xxx[i] * tbe_0 + 2.0 * tr_xxyy_xxxz[i] * tke_0;

        tr_0_0_z_xxyy_xxy[i] = 2.0 * tr_xxyyz_xxy[i] * tbe_0 + 2.0 * tr_xxyy_xxyz[i] * tke_0;

        tr_0_0_z_xxyy_xxz[i] = 2.0 * tr_xxyyz_xxz[i] * tbe_0 + 2.0 * tr_xxyy_xxzz[i] * tke_0 - tr_xxyy_xx[i];

        tr_0_0_z_xxyy_xyy[i] = 2.0 * tr_xxyyz_xyy[i] * tbe_0 + 2.0 * tr_xxyy_xyyz[i] * tke_0;

        tr_0_0_z_xxyy_xyz[i] = 2.0 * tr_xxyyz_xyz[i] * tbe_0 + 2.0 * tr_xxyy_xyzz[i] * tke_0 - tr_xxyy_xy[i];

        tr_0_0_z_xxyy_xzz[i] = 2.0 * tr_xxyyz_xzz[i] * tbe_0 + 2.0 * tr_xxyy_xzzz[i] * tke_0 - 2.0 * tr_xxyy_xz[i];

        tr_0_0_z_xxyy_yyy[i] = 2.0 * tr_xxyyz_yyy[i] * tbe_0 + 2.0 * tr_xxyy_yyyz[i] * tke_0;

        tr_0_0_z_xxyy_yyz[i] = 2.0 * tr_xxyyz_yyz[i] * tbe_0 + 2.0 * tr_xxyy_yyzz[i] * tke_0 - tr_xxyy_yy[i];

        tr_0_0_z_xxyy_yzz[i] = 2.0 * tr_xxyyz_yzz[i] * tbe_0 + 2.0 * tr_xxyy_yzzz[i] * tke_0 - 2.0 * tr_xxyy_yz[i];

        tr_0_0_z_xxyy_zzz[i] = 2.0 * tr_xxyyz_zzz[i] * tbe_0 + 2.0 * tr_xxyy_zzzz[i] * tke_0 - 3.0 * tr_xxyy_zz[i];
    }

    // Set up 340-350 components of targeted buffer : GF

    auto tr_0_0_z_xxyz_xxx = pbuffer.data(idx_op_geom_010_gf + 340);

    auto tr_0_0_z_xxyz_xxy = pbuffer.data(idx_op_geom_010_gf + 341);

    auto tr_0_0_z_xxyz_xxz = pbuffer.data(idx_op_geom_010_gf + 342);

    auto tr_0_0_z_xxyz_xyy = pbuffer.data(idx_op_geom_010_gf + 343);

    auto tr_0_0_z_xxyz_xyz = pbuffer.data(idx_op_geom_010_gf + 344);

    auto tr_0_0_z_xxyz_xzz = pbuffer.data(idx_op_geom_010_gf + 345);

    auto tr_0_0_z_xxyz_yyy = pbuffer.data(idx_op_geom_010_gf + 346);

    auto tr_0_0_z_xxyz_yyz = pbuffer.data(idx_op_geom_010_gf + 347);

    auto tr_0_0_z_xxyz_yzz = pbuffer.data(idx_op_geom_010_gf + 348);

    auto tr_0_0_z_xxyz_zzz = pbuffer.data(idx_op_geom_010_gf + 349);

    #pragma omp simd aligned(tr_0_0_z_xxyz_xxx, tr_0_0_z_xxyz_xxy, tr_0_0_z_xxyz_xxz, tr_0_0_z_xxyz_xyy, tr_0_0_z_xxyz_xyz, tr_0_0_z_xxyz_xzz, tr_0_0_z_xxyz_yyy, tr_0_0_z_xxyz_yyz, tr_0_0_z_xxyz_yzz, tr_0_0_z_xxyz_zzz, tr_xxy_xxx, tr_xxy_xxy, tr_xxy_xxz, tr_xxy_xyy, tr_xxy_xyz, tr_xxy_xzz, tr_xxy_yyy, tr_xxy_yyz, tr_xxy_yzz, tr_xxy_zzz, tr_xxyz_xx, tr_xxyz_xxxz, tr_xxyz_xxyz, tr_xxyz_xxzz, tr_xxyz_xy, tr_xxyz_xyyz, tr_xxyz_xyzz, tr_xxyz_xz, tr_xxyz_xzzz, tr_xxyz_yy, tr_xxyz_yyyz, tr_xxyz_yyzz, tr_xxyz_yz, tr_xxyz_yzzz, tr_xxyz_zz, tr_xxyz_zzzz, tr_xxyzz_xxx, tr_xxyzz_xxy, tr_xxyzz_xxz, tr_xxyzz_xyy, tr_xxyzz_xyz, tr_xxyzz_xzz, tr_xxyzz_yyy, tr_xxyzz_yyz, tr_xxyzz_yzz, tr_xxyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyz_xxx[i] = 2.0 * tr_xxyzz_xxx[i] * tbe_0 + 2.0 * tr_xxyz_xxxz[i] * tke_0 - tr_xxy_xxx[i];

        tr_0_0_z_xxyz_xxy[i] = 2.0 * tr_xxyzz_xxy[i] * tbe_0 + 2.0 * tr_xxyz_xxyz[i] * tke_0 - tr_xxy_xxy[i];

        tr_0_0_z_xxyz_xxz[i] = 2.0 * tr_xxyzz_xxz[i] * tbe_0 + 2.0 * tr_xxyz_xxzz[i] * tke_0 - tr_xxy_xxz[i] - tr_xxyz_xx[i];

        tr_0_0_z_xxyz_xyy[i] = 2.0 * tr_xxyzz_xyy[i] * tbe_0 + 2.0 * tr_xxyz_xyyz[i] * tke_0 - tr_xxy_xyy[i];

        tr_0_0_z_xxyz_xyz[i] = 2.0 * tr_xxyzz_xyz[i] * tbe_0 + 2.0 * tr_xxyz_xyzz[i] * tke_0 - tr_xxy_xyz[i] - tr_xxyz_xy[i];

        tr_0_0_z_xxyz_xzz[i] = 2.0 * tr_xxyzz_xzz[i] * tbe_0 + 2.0 * tr_xxyz_xzzz[i] * tke_0 - tr_xxy_xzz[i] - 2.0 * tr_xxyz_xz[i];

        tr_0_0_z_xxyz_yyy[i] = 2.0 * tr_xxyzz_yyy[i] * tbe_0 + 2.0 * tr_xxyz_yyyz[i] * tke_0 - tr_xxy_yyy[i];

        tr_0_0_z_xxyz_yyz[i] = 2.0 * tr_xxyzz_yyz[i] * tbe_0 + 2.0 * tr_xxyz_yyzz[i] * tke_0 - tr_xxy_yyz[i] - tr_xxyz_yy[i];

        tr_0_0_z_xxyz_yzz[i] = 2.0 * tr_xxyzz_yzz[i] * tbe_0 + 2.0 * tr_xxyz_yzzz[i] * tke_0 - tr_xxy_yzz[i] - 2.0 * tr_xxyz_yz[i];

        tr_0_0_z_xxyz_zzz[i] = 2.0 * tr_xxyzz_zzz[i] * tbe_0 + 2.0 * tr_xxyz_zzzz[i] * tke_0 - tr_xxy_zzz[i] - 3.0 * tr_xxyz_zz[i];
    }

    // Set up 350-360 components of targeted buffer : GF

    auto tr_0_0_z_xxzz_xxx = pbuffer.data(idx_op_geom_010_gf + 350);

    auto tr_0_0_z_xxzz_xxy = pbuffer.data(idx_op_geom_010_gf + 351);

    auto tr_0_0_z_xxzz_xxz = pbuffer.data(idx_op_geom_010_gf + 352);

    auto tr_0_0_z_xxzz_xyy = pbuffer.data(idx_op_geom_010_gf + 353);

    auto tr_0_0_z_xxzz_xyz = pbuffer.data(idx_op_geom_010_gf + 354);

    auto tr_0_0_z_xxzz_xzz = pbuffer.data(idx_op_geom_010_gf + 355);

    auto tr_0_0_z_xxzz_yyy = pbuffer.data(idx_op_geom_010_gf + 356);

    auto tr_0_0_z_xxzz_yyz = pbuffer.data(idx_op_geom_010_gf + 357);

    auto tr_0_0_z_xxzz_yzz = pbuffer.data(idx_op_geom_010_gf + 358);

    auto tr_0_0_z_xxzz_zzz = pbuffer.data(idx_op_geom_010_gf + 359);

    #pragma omp simd aligned(tr_0_0_z_xxzz_xxx, tr_0_0_z_xxzz_xxy, tr_0_0_z_xxzz_xxz, tr_0_0_z_xxzz_xyy, tr_0_0_z_xxzz_xyz, tr_0_0_z_xxzz_xzz, tr_0_0_z_xxzz_yyy, tr_0_0_z_xxzz_yyz, tr_0_0_z_xxzz_yzz, tr_0_0_z_xxzz_zzz, tr_xxz_xxx, tr_xxz_xxy, tr_xxz_xxz, tr_xxz_xyy, tr_xxz_xyz, tr_xxz_xzz, tr_xxz_yyy, tr_xxz_yyz, tr_xxz_yzz, tr_xxz_zzz, tr_xxzz_xx, tr_xxzz_xxxz, tr_xxzz_xxyz, tr_xxzz_xxzz, tr_xxzz_xy, tr_xxzz_xyyz, tr_xxzz_xyzz, tr_xxzz_xz, tr_xxzz_xzzz, tr_xxzz_yy, tr_xxzz_yyyz, tr_xxzz_yyzz, tr_xxzz_yz, tr_xxzz_yzzz, tr_xxzz_zz, tr_xxzz_zzzz, tr_xxzzz_xxx, tr_xxzzz_xxy, tr_xxzzz_xxz, tr_xxzzz_xyy, tr_xxzzz_xyz, tr_xxzzz_xzz, tr_xxzzz_yyy, tr_xxzzz_yyz, tr_xxzzz_yzz, tr_xxzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxzz_xxx[i] = 2.0 * tr_xxzzz_xxx[i] * tbe_0 + 2.0 * tr_xxzz_xxxz[i] * tke_0 - 2.0 * tr_xxz_xxx[i];

        tr_0_0_z_xxzz_xxy[i] = 2.0 * tr_xxzzz_xxy[i] * tbe_0 + 2.0 * tr_xxzz_xxyz[i] * tke_0 - 2.0 * tr_xxz_xxy[i];

        tr_0_0_z_xxzz_xxz[i] = 2.0 * tr_xxzzz_xxz[i] * tbe_0 + 2.0 * tr_xxzz_xxzz[i] * tke_0 - 2.0 * tr_xxz_xxz[i] - tr_xxzz_xx[i];

        tr_0_0_z_xxzz_xyy[i] = 2.0 * tr_xxzzz_xyy[i] * tbe_0 + 2.0 * tr_xxzz_xyyz[i] * tke_0 - 2.0 * tr_xxz_xyy[i];

        tr_0_0_z_xxzz_xyz[i] = 2.0 * tr_xxzzz_xyz[i] * tbe_0 + 2.0 * tr_xxzz_xyzz[i] * tke_0 - 2.0 * tr_xxz_xyz[i] - tr_xxzz_xy[i];

        tr_0_0_z_xxzz_xzz[i] = 2.0 * tr_xxzzz_xzz[i] * tbe_0 + 2.0 * tr_xxzz_xzzz[i] * tke_0 - 2.0 * tr_xxz_xzz[i] - 2.0 * tr_xxzz_xz[i];

        tr_0_0_z_xxzz_yyy[i] = 2.0 * tr_xxzzz_yyy[i] * tbe_0 + 2.0 * tr_xxzz_yyyz[i] * tke_0 - 2.0 * tr_xxz_yyy[i];

        tr_0_0_z_xxzz_yyz[i] = 2.0 * tr_xxzzz_yyz[i] * tbe_0 + 2.0 * tr_xxzz_yyzz[i] * tke_0 - 2.0 * tr_xxz_yyz[i] - tr_xxzz_yy[i];

        tr_0_0_z_xxzz_yzz[i] = 2.0 * tr_xxzzz_yzz[i] * tbe_0 + 2.0 * tr_xxzz_yzzz[i] * tke_0 - 2.0 * tr_xxz_yzz[i] - 2.0 * tr_xxzz_yz[i];

        tr_0_0_z_xxzz_zzz[i] = 2.0 * tr_xxzzz_zzz[i] * tbe_0 + 2.0 * tr_xxzz_zzzz[i] * tke_0 - 2.0 * tr_xxz_zzz[i] - 3.0 * tr_xxzz_zz[i];
    }

    // Set up 360-370 components of targeted buffer : GF

    auto tr_0_0_z_xyyy_xxx = pbuffer.data(idx_op_geom_010_gf + 360);

    auto tr_0_0_z_xyyy_xxy = pbuffer.data(idx_op_geom_010_gf + 361);

    auto tr_0_0_z_xyyy_xxz = pbuffer.data(idx_op_geom_010_gf + 362);

    auto tr_0_0_z_xyyy_xyy = pbuffer.data(idx_op_geom_010_gf + 363);

    auto tr_0_0_z_xyyy_xyz = pbuffer.data(idx_op_geom_010_gf + 364);

    auto tr_0_0_z_xyyy_xzz = pbuffer.data(idx_op_geom_010_gf + 365);

    auto tr_0_0_z_xyyy_yyy = pbuffer.data(idx_op_geom_010_gf + 366);

    auto tr_0_0_z_xyyy_yyz = pbuffer.data(idx_op_geom_010_gf + 367);

    auto tr_0_0_z_xyyy_yzz = pbuffer.data(idx_op_geom_010_gf + 368);

    auto tr_0_0_z_xyyy_zzz = pbuffer.data(idx_op_geom_010_gf + 369);

    #pragma omp simd aligned(tr_0_0_z_xyyy_xxx, tr_0_0_z_xyyy_xxy, tr_0_0_z_xyyy_xxz, tr_0_0_z_xyyy_xyy, tr_0_0_z_xyyy_xyz, tr_0_0_z_xyyy_xzz, tr_0_0_z_xyyy_yyy, tr_0_0_z_xyyy_yyz, tr_0_0_z_xyyy_yzz, tr_0_0_z_xyyy_zzz, tr_xyyy_xx, tr_xyyy_xxxz, tr_xyyy_xxyz, tr_xyyy_xxzz, tr_xyyy_xy, tr_xyyy_xyyz, tr_xyyy_xyzz, tr_xyyy_xz, tr_xyyy_xzzz, tr_xyyy_yy, tr_xyyy_yyyz, tr_xyyy_yyzz, tr_xyyy_yz, tr_xyyy_yzzz, tr_xyyy_zz, tr_xyyy_zzzz, tr_xyyyz_xxx, tr_xyyyz_xxy, tr_xyyyz_xxz, tr_xyyyz_xyy, tr_xyyyz_xyz, tr_xyyyz_xzz, tr_xyyyz_yyy, tr_xyyyz_yyz, tr_xyyyz_yzz, tr_xyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyy_xxx[i] = 2.0 * tr_xyyyz_xxx[i] * tbe_0 + 2.0 * tr_xyyy_xxxz[i] * tke_0;

        tr_0_0_z_xyyy_xxy[i] = 2.0 * tr_xyyyz_xxy[i] * tbe_0 + 2.0 * tr_xyyy_xxyz[i] * tke_0;

        tr_0_0_z_xyyy_xxz[i] = 2.0 * tr_xyyyz_xxz[i] * tbe_0 + 2.0 * tr_xyyy_xxzz[i] * tke_0 - tr_xyyy_xx[i];

        tr_0_0_z_xyyy_xyy[i] = 2.0 * tr_xyyyz_xyy[i] * tbe_0 + 2.0 * tr_xyyy_xyyz[i] * tke_0;

        tr_0_0_z_xyyy_xyz[i] = 2.0 * tr_xyyyz_xyz[i] * tbe_0 + 2.0 * tr_xyyy_xyzz[i] * tke_0 - tr_xyyy_xy[i];

        tr_0_0_z_xyyy_xzz[i] = 2.0 * tr_xyyyz_xzz[i] * tbe_0 + 2.0 * tr_xyyy_xzzz[i] * tke_0 - 2.0 * tr_xyyy_xz[i];

        tr_0_0_z_xyyy_yyy[i] = 2.0 * tr_xyyyz_yyy[i] * tbe_0 + 2.0 * tr_xyyy_yyyz[i] * tke_0;

        tr_0_0_z_xyyy_yyz[i] = 2.0 * tr_xyyyz_yyz[i] * tbe_0 + 2.0 * tr_xyyy_yyzz[i] * tke_0 - tr_xyyy_yy[i];

        tr_0_0_z_xyyy_yzz[i] = 2.0 * tr_xyyyz_yzz[i] * tbe_0 + 2.0 * tr_xyyy_yzzz[i] * tke_0 - 2.0 * tr_xyyy_yz[i];

        tr_0_0_z_xyyy_zzz[i] = 2.0 * tr_xyyyz_zzz[i] * tbe_0 + 2.0 * tr_xyyy_zzzz[i] * tke_0 - 3.0 * tr_xyyy_zz[i];
    }

    // Set up 370-380 components of targeted buffer : GF

    auto tr_0_0_z_xyyz_xxx = pbuffer.data(idx_op_geom_010_gf + 370);

    auto tr_0_0_z_xyyz_xxy = pbuffer.data(idx_op_geom_010_gf + 371);

    auto tr_0_0_z_xyyz_xxz = pbuffer.data(idx_op_geom_010_gf + 372);

    auto tr_0_0_z_xyyz_xyy = pbuffer.data(idx_op_geom_010_gf + 373);

    auto tr_0_0_z_xyyz_xyz = pbuffer.data(idx_op_geom_010_gf + 374);

    auto tr_0_0_z_xyyz_xzz = pbuffer.data(idx_op_geom_010_gf + 375);

    auto tr_0_0_z_xyyz_yyy = pbuffer.data(idx_op_geom_010_gf + 376);

    auto tr_0_0_z_xyyz_yyz = pbuffer.data(idx_op_geom_010_gf + 377);

    auto tr_0_0_z_xyyz_yzz = pbuffer.data(idx_op_geom_010_gf + 378);

    auto tr_0_0_z_xyyz_zzz = pbuffer.data(idx_op_geom_010_gf + 379);

    #pragma omp simd aligned(tr_0_0_z_xyyz_xxx, tr_0_0_z_xyyz_xxy, tr_0_0_z_xyyz_xxz, tr_0_0_z_xyyz_xyy, tr_0_0_z_xyyz_xyz, tr_0_0_z_xyyz_xzz, tr_0_0_z_xyyz_yyy, tr_0_0_z_xyyz_yyz, tr_0_0_z_xyyz_yzz, tr_0_0_z_xyyz_zzz, tr_xyy_xxx, tr_xyy_xxy, tr_xyy_xxz, tr_xyy_xyy, tr_xyy_xyz, tr_xyy_xzz, tr_xyy_yyy, tr_xyy_yyz, tr_xyy_yzz, tr_xyy_zzz, tr_xyyz_xx, tr_xyyz_xxxz, tr_xyyz_xxyz, tr_xyyz_xxzz, tr_xyyz_xy, tr_xyyz_xyyz, tr_xyyz_xyzz, tr_xyyz_xz, tr_xyyz_xzzz, tr_xyyz_yy, tr_xyyz_yyyz, tr_xyyz_yyzz, tr_xyyz_yz, tr_xyyz_yzzz, tr_xyyz_zz, tr_xyyz_zzzz, tr_xyyzz_xxx, tr_xyyzz_xxy, tr_xyyzz_xxz, tr_xyyzz_xyy, tr_xyyzz_xyz, tr_xyyzz_xzz, tr_xyyzz_yyy, tr_xyyzz_yyz, tr_xyyzz_yzz, tr_xyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyz_xxx[i] = 2.0 * tr_xyyzz_xxx[i] * tbe_0 + 2.0 * tr_xyyz_xxxz[i] * tke_0 - tr_xyy_xxx[i];

        tr_0_0_z_xyyz_xxy[i] = 2.0 * tr_xyyzz_xxy[i] * tbe_0 + 2.0 * tr_xyyz_xxyz[i] * tke_0 - tr_xyy_xxy[i];

        tr_0_0_z_xyyz_xxz[i] = 2.0 * tr_xyyzz_xxz[i] * tbe_0 + 2.0 * tr_xyyz_xxzz[i] * tke_0 - tr_xyy_xxz[i] - tr_xyyz_xx[i];

        tr_0_0_z_xyyz_xyy[i] = 2.0 * tr_xyyzz_xyy[i] * tbe_0 + 2.0 * tr_xyyz_xyyz[i] * tke_0 - tr_xyy_xyy[i];

        tr_0_0_z_xyyz_xyz[i] = 2.0 * tr_xyyzz_xyz[i] * tbe_0 + 2.0 * tr_xyyz_xyzz[i] * tke_0 - tr_xyy_xyz[i] - tr_xyyz_xy[i];

        tr_0_0_z_xyyz_xzz[i] = 2.0 * tr_xyyzz_xzz[i] * tbe_0 + 2.0 * tr_xyyz_xzzz[i] * tke_0 - tr_xyy_xzz[i] - 2.0 * tr_xyyz_xz[i];

        tr_0_0_z_xyyz_yyy[i] = 2.0 * tr_xyyzz_yyy[i] * tbe_0 + 2.0 * tr_xyyz_yyyz[i] * tke_0 - tr_xyy_yyy[i];

        tr_0_0_z_xyyz_yyz[i] = 2.0 * tr_xyyzz_yyz[i] * tbe_0 + 2.0 * tr_xyyz_yyzz[i] * tke_0 - tr_xyy_yyz[i] - tr_xyyz_yy[i];

        tr_0_0_z_xyyz_yzz[i] = 2.0 * tr_xyyzz_yzz[i] * tbe_0 + 2.0 * tr_xyyz_yzzz[i] * tke_0 - tr_xyy_yzz[i] - 2.0 * tr_xyyz_yz[i];

        tr_0_0_z_xyyz_zzz[i] = 2.0 * tr_xyyzz_zzz[i] * tbe_0 + 2.0 * tr_xyyz_zzzz[i] * tke_0 - tr_xyy_zzz[i] - 3.0 * tr_xyyz_zz[i];
    }

    // Set up 380-390 components of targeted buffer : GF

    auto tr_0_0_z_xyzz_xxx = pbuffer.data(idx_op_geom_010_gf + 380);

    auto tr_0_0_z_xyzz_xxy = pbuffer.data(idx_op_geom_010_gf + 381);

    auto tr_0_0_z_xyzz_xxz = pbuffer.data(idx_op_geom_010_gf + 382);

    auto tr_0_0_z_xyzz_xyy = pbuffer.data(idx_op_geom_010_gf + 383);

    auto tr_0_0_z_xyzz_xyz = pbuffer.data(idx_op_geom_010_gf + 384);

    auto tr_0_0_z_xyzz_xzz = pbuffer.data(idx_op_geom_010_gf + 385);

    auto tr_0_0_z_xyzz_yyy = pbuffer.data(idx_op_geom_010_gf + 386);

    auto tr_0_0_z_xyzz_yyz = pbuffer.data(idx_op_geom_010_gf + 387);

    auto tr_0_0_z_xyzz_yzz = pbuffer.data(idx_op_geom_010_gf + 388);

    auto tr_0_0_z_xyzz_zzz = pbuffer.data(idx_op_geom_010_gf + 389);

    #pragma omp simd aligned(tr_0_0_z_xyzz_xxx, tr_0_0_z_xyzz_xxy, tr_0_0_z_xyzz_xxz, tr_0_0_z_xyzz_xyy, tr_0_0_z_xyzz_xyz, tr_0_0_z_xyzz_xzz, tr_0_0_z_xyzz_yyy, tr_0_0_z_xyzz_yyz, tr_0_0_z_xyzz_yzz, tr_0_0_z_xyzz_zzz, tr_xyz_xxx, tr_xyz_xxy, tr_xyz_xxz, tr_xyz_xyy, tr_xyz_xyz, tr_xyz_xzz, tr_xyz_yyy, tr_xyz_yyz, tr_xyz_yzz, tr_xyz_zzz, tr_xyzz_xx, tr_xyzz_xxxz, tr_xyzz_xxyz, tr_xyzz_xxzz, tr_xyzz_xy, tr_xyzz_xyyz, tr_xyzz_xyzz, tr_xyzz_xz, tr_xyzz_xzzz, tr_xyzz_yy, tr_xyzz_yyyz, tr_xyzz_yyzz, tr_xyzz_yz, tr_xyzz_yzzz, tr_xyzz_zz, tr_xyzz_zzzz, tr_xyzzz_xxx, tr_xyzzz_xxy, tr_xyzzz_xxz, tr_xyzzz_xyy, tr_xyzzz_xyz, tr_xyzzz_xzz, tr_xyzzz_yyy, tr_xyzzz_yyz, tr_xyzzz_yzz, tr_xyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyzz_xxx[i] = 2.0 * tr_xyzzz_xxx[i] * tbe_0 + 2.0 * tr_xyzz_xxxz[i] * tke_0 - 2.0 * tr_xyz_xxx[i];

        tr_0_0_z_xyzz_xxy[i] = 2.0 * tr_xyzzz_xxy[i] * tbe_0 + 2.0 * tr_xyzz_xxyz[i] * tke_0 - 2.0 * tr_xyz_xxy[i];

        tr_0_0_z_xyzz_xxz[i] = 2.0 * tr_xyzzz_xxz[i] * tbe_0 + 2.0 * tr_xyzz_xxzz[i] * tke_0 - 2.0 * tr_xyz_xxz[i] - tr_xyzz_xx[i];

        tr_0_0_z_xyzz_xyy[i] = 2.0 * tr_xyzzz_xyy[i] * tbe_0 + 2.0 * tr_xyzz_xyyz[i] * tke_0 - 2.0 * tr_xyz_xyy[i];

        tr_0_0_z_xyzz_xyz[i] = 2.0 * tr_xyzzz_xyz[i] * tbe_0 + 2.0 * tr_xyzz_xyzz[i] * tke_0 - 2.0 * tr_xyz_xyz[i] - tr_xyzz_xy[i];

        tr_0_0_z_xyzz_xzz[i] = 2.0 * tr_xyzzz_xzz[i] * tbe_0 + 2.0 * tr_xyzz_xzzz[i] * tke_0 - 2.0 * tr_xyz_xzz[i] - 2.0 * tr_xyzz_xz[i];

        tr_0_0_z_xyzz_yyy[i] = 2.0 * tr_xyzzz_yyy[i] * tbe_0 + 2.0 * tr_xyzz_yyyz[i] * tke_0 - 2.0 * tr_xyz_yyy[i];

        tr_0_0_z_xyzz_yyz[i] = 2.0 * tr_xyzzz_yyz[i] * tbe_0 + 2.0 * tr_xyzz_yyzz[i] * tke_0 - 2.0 * tr_xyz_yyz[i] - tr_xyzz_yy[i];

        tr_0_0_z_xyzz_yzz[i] = 2.0 * tr_xyzzz_yzz[i] * tbe_0 + 2.0 * tr_xyzz_yzzz[i] * tke_0 - 2.0 * tr_xyz_yzz[i] - 2.0 * tr_xyzz_yz[i];

        tr_0_0_z_xyzz_zzz[i] = 2.0 * tr_xyzzz_zzz[i] * tbe_0 + 2.0 * tr_xyzz_zzzz[i] * tke_0 - 2.0 * tr_xyz_zzz[i] - 3.0 * tr_xyzz_zz[i];
    }

    // Set up 390-400 components of targeted buffer : GF

    auto tr_0_0_z_xzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 390);

    auto tr_0_0_z_xzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 391);

    auto tr_0_0_z_xzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 392);

    auto tr_0_0_z_xzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 393);

    auto tr_0_0_z_xzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 394);

    auto tr_0_0_z_xzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 395);

    auto tr_0_0_z_xzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 396);

    auto tr_0_0_z_xzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 397);

    auto tr_0_0_z_xzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 398);

    auto tr_0_0_z_xzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 399);

    #pragma omp simd aligned(tr_0_0_z_xzzz_xxx, tr_0_0_z_xzzz_xxy, tr_0_0_z_xzzz_xxz, tr_0_0_z_xzzz_xyy, tr_0_0_z_xzzz_xyz, tr_0_0_z_xzzz_xzz, tr_0_0_z_xzzz_yyy, tr_0_0_z_xzzz_yyz, tr_0_0_z_xzzz_yzz, tr_0_0_z_xzzz_zzz, tr_xzz_xxx, tr_xzz_xxy, tr_xzz_xxz, tr_xzz_xyy, tr_xzz_xyz, tr_xzz_xzz, tr_xzz_yyy, tr_xzz_yyz, tr_xzz_yzz, tr_xzz_zzz, tr_xzzz_xx, tr_xzzz_xxxz, tr_xzzz_xxyz, tr_xzzz_xxzz, tr_xzzz_xy, tr_xzzz_xyyz, tr_xzzz_xyzz, tr_xzzz_xz, tr_xzzz_xzzz, tr_xzzz_yy, tr_xzzz_yyyz, tr_xzzz_yyzz, tr_xzzz_yz, tr_xzzz_yzzz, tr_xzzz_zz, tr_xzzz_zzzz, tr_xzzzz_xxx, tr_xzzzz_xxy, tr_xzzzz_xxz, tr_xzzzz_xyy, tr_xzzzz_xyz, tr_xzzzz_xzz, tr_xzzzz_yyy, tr_xzzzz_yyz, tr_xzzzz_yzz, tr_xzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzzz_xxx[i] = 2.0 * tr_xzzzz_xxx[i] * tbe_0 + 2.0 * tr_xzzz_xxxz[i] * tke_0 - 3.0 * tr_xzz_xxx[i];

        tr_0_0_z_xzzz_xxy[i] = 2.0 * tr_xzzzz_xxy[i] * tbe_0 + 2.0 * tr_xzzz_xxyz[i] * tke_0 - 3.0 * tr_xzz_xxy[i];

        tr_0_0_z_xzzz_xxz[i] = 2.0 * tr_xzzzz_xxz[i] * tbe_0 + 2.0 * tr_xzzz_xxzz[i] * tke_0 - 3.0 * tr_xzz_xxz[i] - tr_xzzz_xx[i];

        tr_0_0_z_xzzz_xyy[i] = 2.0 * tr_xzzzz_xyy[i] * tbe_0 + 2.0 * tr_xzzz_xyyz[i] * tke_0 - 3.0 * tr_xzz_xyy[i];

        tr_0_0_z_xzzz_xyz[i] = 2.0 * tr_xzzzz_xyz[i] * tbe_0 + 2.0 * tr_xzzz_xyzz[i] * tke_0 - 3.0 * tr_xzz_xyz[i] - tr_xzzz_xy[i];

        tr_0_0_z_xzzz_xzz[i] = 2.0 * tr_xzzzz_xzz[i] * tbe_0 + 2.0 * tr_xzzz_xzzz[i] * tke_0 - 3.0 * tr_xzz_xzz[i] - 2.0 * tr_xzzz_xz[i];

        tr_0_0_z_xzzz_yyy[i] = 2.0 * tr_xzzzz_yyy[i] * tbe_0 + 2.0 * tr_xzzz_yyyz[i] * tke_0 - 3.0 * tr_xzz_yyy[i];

        tr_0_0_z_xzzz_yyz[i] = 2.0 * tr_xzzzz_yyz[i] * tbe_0 + 2.0 * tr_xzzz_yyzz[i] * tke_0 - 3.0 * tr_xzz_yyz[i] - tr_xzzz_yy[i];

        tr_0_0_z_xzzz_yzz[i] = 2.0 * tr_xzzzz_yzz[i] * tbe_0 + 2.0 * tr_xzzz_yzzz[i] * tke_0 - 3.0 * tr_xzz_yzz[i] - 2.0 * tr_xzzz_yz[i];

        tr_0_0_z_xzzz_zzz[i] = 2.0 * tr_xzzzz_zzz[i] * tbe_0 + 2.0 * tr_xzzz_zzzz[i] * tke_0 - 3.0 * tr_xzz_zzz[i] - 3.0 * tr_xzzz_zz[i];
    }

    // Set up 400-410 components of targeted buffer : GF

    auto tr_0_0_z_yyyy_xxx = pbuffer.data(idx_op_geom_010_gf + 400);

    auto tr_0_0_z_yyyy_xxy = pbuffer.data(idx_op_geom_010_gf + 401);

    auto tr_0_0_z_yyyy_xxz = pbuffer.data(idx_op_geom_010_gf + 402);

    auto tr_0_0_z_yyyy_xyy = pbuffer.data(idx_op_geom_010_gf + 403);

    auto tr_0_0_z_yyyy_xyz = pbuffer.data(idx_op_geom_010_gf + 404);

    auto tr_0_0_z_yyyy_xzz = pbuffer.data(idx_op_geom_010_gf + 405);

    auto tr_0_0_z_yyyy_yyy = pbuffer.data(idx_op_geom_010_gf + 406);

    auto tr_0_0_z_yyyy_yyz = pbuffer.data(idx_op_geom_010_gf + 407);

    auto tr_0_0_z_yyyy_yzz = pbuffer.data(idx_op_geom_010_gf + 408);

    auto tr_0_0_z_yyyy_zzz = pbuffer.data(idx_op_geom_010_gf + 409);

    #pragma omp simd aligned(tr_0_0_z_yyyy_xxx, tr_0_0_z_yyyy_xxy, tr_0_0_z_yyyy_xxz, tr_0_0_z_yyyy_xyy, tr_0_0_z_yyyy_xyz, tr_0_0_z_yyyy_xzz, tr_0_0_z_yyyy_yyy, tr_0_0_z_yyyy_yyz, tr_0_0_z_yyyy_yzz, tr_0_0_z_yyyy_zzz, tr_yyyy_xx, tr_yyyy_xxxz, tr_yyyy_xxyz, tr_yyyy_xxzz, tr_yyyy_xy, tr_yyyy_xyyz, tr_yyyy_xyzz, tr_yyyy_xz, tr_yyyy_xzzz, tr_yyyy_yy, tr_yyyy_yyyz, tr_yyyy_yyzz, tr_yyyy_yz, tr_yyyy_yzzz, tr_yyyy_zz, tr_yyyy_zzzz, tr_yyyyz_xxx, tr_yyyyz_xxy, tr_yyyyz_xxz, tr_yyyyz_xyy, tr_yyyyz_xyz, tr_yyyyz_xzz, tr_yyyyz_yyy, tr_yyyyz_yyz, tr_yyyyz_yzz, tr_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyy_xxx[i] = 2.0 * tr_yyyyz_xxx[i] * tbe_0 + 2.0 * tr_yyyy_xxxz[i] * tke_0;

        tr_0_0_z_yyyy_xxy[i] = 2.0 * tr_yyyyz_xxy[i] * tbe_0 + 2.0 * tr_yyyy_xxyz[i] * tke_0;

        tr_0_0_z_yyyy_xxz[i] = 2.0 * tr_yyyyz_xxz[i] * tbe_0 + 2.0 * tr_yyyy_xxzz[i] * tke_0 - tr_yyyy_xx[i];

        tr_0_0_z_yyyy_xyy[i] = 2.0 * tr_yyyyz_xyy[i] * tbe_0 + 2.0 * tr_yyyy_xyyz[i] * tke_0;

        tr_0_0_z_yyyy_xyz[i] = 2.0 * tr_yyyyz_xyz[i] * tbe_0 + 2.0 * tr_yyyy_xyzz[i] * tke_0 - tr_yyyy_xy[i];

        tr_0_0_z_yyyy_xzz[i] = 2.0 * tr_yyyyz_xzz[i] * tbe_0 + 2.0 * tr_yyyy_xzzz[i] * tke_0 - 2.0 * tr_yyyy_xz[i];

        tr_0_0_z_yyyy_yyy[i] = 2.0 * tr_yyyyz_yyy[i] * tbe_0 + 2.0 * tr_yyyy_yyyz[i] * tke_0;

        tr_0_0_z_yyyy_yyz[i] = 2.0 * tr_yyyyz_yyz[i] * tbe_0 + 2.0 * tr_yyyy_yyzz[i] * tke_0 - tr_yyyy_yy[i];

        tr_0_0_z_yyyy_yzz[i] = 2.0 * tr_yyyyz_yzz[i] * tbe_0 + 2.0 * tr_yyyy_yzzz[i] * tke_0 - 2.0 * tr_yyyy_yz[i];

        tr_0_0_z_yyyy_zzz[i] = 2.0 * tr_yyyyz_zzz[i] * tbe_0 + 2.0 * tr_yyyy_zzzz[i] * tke_0 - 3.0 * tr_yyyy_zz[i];
    }

    // Set up 410-420 components of targeted buffer : GF

    auto tr_0_0_z_yyyz_xxx = pbuffer.data(idx_op_geom_010_gf + 410);

    auto tr_0_0_z_yyyz_xxy = pbuffer.data(idx_op_geom_010_gf + 411);

    auto tr_0_0_z_yyyz_xxz = pbuffer.data(idx_op_geom_010_gf + 412);

    auto tr_0_0_z_yyyz_xyy = pbuffer.data(idx_op_geom_010_gf + 413);

    auto tr_0_0_z_yyyz_xyz = pbuffer.data(idx_op_geom_010_gf + 414);

    auto tr_0_0_z_yyyz_xzz = pbuffer.data(idx_op_geom_010_gf + 415);

    auto tr_0_0_z_yyyz_yyy = pbuffer.data(idx_op_geom_010_gf + 416);

    auto tr_0_0_z_yyyz_yyz = pbuffer.data(idx_op_geom_010_gf + 417);

    auto tr_0_0_z_yyyz_yzz = pbuffer.data(idx_op_geom_010_gf + 418);

    auto tr_0_0_z_yyyz_zzz = pbuffer.data(idx_op_geom_010_gf + 419);

    #pragma omp simd aligned(tr_0_0_z_yyyz_xxx, tr_0_0_z_yyyz_xxy, tr_0_0_z_yyyz_xxz, tr_0_0_z_yyyz_xyy, tr_0_0_z_yyyz_xyz, tr_0_0_z_yyyz_xzz, tr_0_0_z_yyyz_yyy, tr_0_0_z_yyyz_yyz, tr_0_0_z_yyyz_yzz, tr_0_0_z_yyyz_zzz, tr_yyy_xxx, tr_yyy_xxy, tr_yyy_xxz, tr_yyy_xyy, tr_yyy_xyz, tr_yyy_xzz, tr_yyy_yyy, tr_yyy_yyz, tr_yyy_yzz, tr_yyy_zzz, tr_yyyz_xx, tr_yyyz_xxxz, tr_yyyz_xxyz, tr_yyyz_xxzz, tr_yyyz_xy, tr_yyyz_xyyz, tr_yyyz_xyzz, tr_yyyz_xz, tr_yyyz_xzzz, tr_yyyz_yy, tr_yyyz_yyyz, tr_yyyz_yyzz, tr_yyyz_yz, tr_yyyz_yzzz, tr_yyyz_zz, tr_yyyz_zzzz, tr_yyyzz_xxx, tr_yyyzz_xxy, tr_yyyzz_xxz, tr_yyyzz_xyy, tr_yyyzz_xyz, tr_yyyzz_xzz, tr_yyyzz_yyy, tr_yyyzz_yyz, tr_yyyzz_yzz, tr_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyz_xxx[i] = 2.0 * tr_yyyzz_xxx[i] * tbe_0 + 2.0 * tr_yyyz_xxxz[i] * tke_0 - tr_yyy_xxx[i];

        tr_0_0_z_yyyz_xxy[i] = 2.0 * tr_yyyzz_xxy[i] * tbe_0 + 2.0 * tr_yyyz_xxyz[i] * tke_0 - tr_yyy_xxy[i];

        tr_0_0_z_yyyz_xxz[i] = 2.0 * tr_yyyzz_xxz[i] * tbe_0 + 2.0 * tr_yyyz_xxzz[i] * tke_0 - tr_yyy_xxz[i] - tr_yyyz_xx[i];

        tr_0_0_z_yyyz_xyy[i] = 2.0 * tr_yyyzz_xyy[i] * tbe_0 + 2.0 * tr_yyyz_xyyz[i] * tke_0 - tr_yyy_xyy[i];

        tr_0_0_z_yyyz_xyz[i] = 2.0 * tr_yyyzz_xyz[i] * tbe_0 + 2.0 * tr_yyyz_xyzz[i] * tke_0 - tr_yyy_xyz[i] - tr_yyyz_xy[i];

        tr_0_0_z_yyyz_xzz[i] = 2.0 * tr_yyyzz_xzz[i] * tbe_0 + 2.0 * tr_yyyz_xzzz[i] * tke_0 - tr_yyy_xzz[i] - 2.0 * tr_yyyz_xz[i];

        tr_0_0_z_yyyz_yyy[i] = 2.0 * tr_yyyzz_yyy[i] * tbe_0 + 2.0 * tr_yyyz_yyyz[i] * tke_0 - tr_yyy_yyy[i];

        tr_0_0_z_yyyz_yyz[i] = 2.0 * tr_yyyzz_yyz[i] * tbe_0 + 2.0 * tr_yyyz_yyzz[i] * tke_0 - tr_yyy_yyz[i] - tr_yyyz_yy[i];

        tr_0_0_z_yyyz_yzz[i] = 2.0 * tr_yyyzz_yzz[i] * tbe_0 + 2.0 * tr_yyyz_yzzz[i] * tke_0 - tr_yyy_yzz[i] - 2.0 * tr_yyyz_yz[i];

        tr_0_0_z_yyyz_zzz[i] = 2.0 * tr_yyyzz_zzz[i] * tbe_0 + 2.0 * tr_yyyz_zzzz[i] * tke_0 - tr_yyy_zzz[i] - 3.0 * tr_yyyz_zz[i];
    }

    // Set up 420-430 components of targeted buffer : GF

    auto tr_0_0_z_yyzz_xxx = pbuffer.data(idx_op_geom_010_gf + 420);

    auto tr_0_0_z_yyzz_xxy = pbuffer.data(idx_op_geom_010_gf + 421);

    auto tr_0_0_z_yyzz_xxz = pbuffer.data(idx_op_geom_010_gf + 422);

    auto tr_0_0_z_yyzz_xyy = pbuffer.data(idx_op_geom_010_gf + 423);

    auto tr_0_0_z_yyzz_xyz = pbuffer.data(idx_op_geom_010_gf + 424);

    auto tr_0_0_z_yyzz_xzz = pbuffer.data(idx_op_geom_010_gf + 425);

    auto tr_0_0_z_yyzz_yyy = pbuffer.data(idx_op_geom_010_gf + 426);

    auto tr_0_0_z_yyzz_yyz = pbuffer.data(idx_op_geom_010_gf + 427);

    auto tr_0_0_z_yyzz_yzz = pbuffer.data(idx_op_geom_010_gf + 428);

    auto tr_0_0_z_yyzz_zzz = pbuffer.data(idx_op_geom_010_gf + 429);

    #pragma omp simd aligned(tr_0_0_z_yyzz_xxx, tr_0_0_z_yyzz_xxy, tr_0_0_z_yyzz_xxz, tr_0_0_z_yyzz_xyy, tr_0_0_z_yyzz_xyz, tr_0_0_z_yyzz_xzz, tr_0_0_z_yyzz_yyy, tr_0_0_z_yyzz_yyz, tr_0_0_z_yyzz_yzz, tr_0_0_z_yyzz_zzz, tr_yyz_xxx, tr_yyz_xxy, tr_yyz_xxz, tr_yyz_xyy, tr_yyz_xyz, tr_yyz_xzz, tr_yyz_yyy, tr_yyz_yyz, tr_yyz_yzz, tr_yyz_zzz, tr_yyzz_xx, tr_yyzz_xxxz, tr_yyzz_xxyz, tr_yyzz_xxzz, tr_yyzz_xy, tr_yyzz_xyyz, tr_yyzz_xyzz, tr_yyzz_xz, tr_yyzz_xzzz, tr_yyzz_yy, tr_yyzz_yyyz, tr_yyzz_yyzz, tr_yyzz_yz, tr_yyzz_yzzz, tr_yyzz_zz, tr_yyzz_zzzz, tr_yyzzz_xxx, tr_yyzzz_xxy, tr_yyzzz_xxz, tr_yyzzz_xyy, tr_yyzzz_xyz, tr_yyzzz_xzz, tr_yyzzz_yyy, tr_yyzzz_yyz, tr_yyzzz_yzz, tr_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyzz_xxx[i] = 2.0 * tr_yyzzz_xxx[i] * tbe_0 + 2.0 * tr_yyzz_xxxz[i] * tke_0 - 2.0 * tr_yyz_xxx[i];

        tr_0_0_z_yyzz_xxy[i] = 2.0 * tr_yyzzz_xxy[i] * tbe_0 + 2.0 * tr_yyzz_xxyz[i] * tke_0 - 2.0 * tr_yyz_xxy[i];

        tr_0_0_z_yyzz_xxz[i] = 2.0 * tr_yyzzz_xxz[i] * tbe_0 + 2.0 * tr_yyzz_xxzz[i] * tke_0 - 2.0 * tr_yyz_xxz[i] - tr_yyzz_xx[i];

        tr_0_0_z_yyzz_xyy[i] = 2.0 * tr_yyzzz_xyy[i] * tbe_0 + 2.0 * tr_yyzz_xyyz[i] * tke_0 - 2.0 * tr_yyz_xyy[i];

        tr_0_0_z_yyzz_xyz[i] = 2.0 * tr_yyzzz_xyz[i] * tbe_0 + 2.0 * tr_yyzz_xyzz[i] * tke_0 - 2.0 * tr_yyz_xyz[i] - tr_yyzz_xy[i];

        tr_0_0_z_yyzz_xzz[i] = 2.0 * tr_yyzzz_xzz[i] * tbe_0 + 2.0 * tr_yyzz_xzzz[i] * tke_0 - 2.0 * tr_yyz_xzz[i] - 2.0 * tr_yyzz_xz[i];

        tr_0_0_z_yyzz_yyy[i] = 2.0 * tr_yyzzz_yyy[i] * tbe_0 + 2.0 * tr_yyzz_yyyz[i] * tke_0 - 2.0 * tr_yyz_yyy[i];

        tr_0_0_z_yyzz_yyz[i] = 2.0 * tr_yyzzz_yyz[i] * tbe_0 + 2.0 * tr_yyzz_yyzz[i] * tke_0 - 2.0 * tr_yyz_yyz[i] - tr_yyzz_yy[i];

        tr_0_0_z_yyzz_yzz[i] = 2.0 * tr_yyzzz_yzz[i] * tbe_0 + 2.0 * tr_yyzz_yzzz[i] * tke_0 - 2.0 * tr_yyz_yzz[i] - 2.0 * tr_yyzz_yz[i];

        tr_0_0_z_yyzz_zzz[i] = 2.0 * tr_yyzzz_zzz[i] * tbe_0 + 2.0 * tr_yyzz_zzzz[i] * tke_0 - 2.0 * tr_yyz_zzz[i] - 3.0 * tr_yyzz_zz[i];
    }

    // Set up 430-440 components of targeted buffer : GF

    auto tr_0_0_z_yzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 430);

    auto tr_0_0_z_yzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 431);

    auto tr_0_0_z_yzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 432);

    auto tr_0_0_z_yzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 433);

    auto tr_0_0_z_yzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 434);

    auto tr_0_0_z_yzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 435);

    auto tr_0_0_z_yzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 436);

    auto tr_0_0_z_yzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 437);

    auto tr_0_0_z_yzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 438);

    auto tr_0_0_z_yzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 439);

    #pragma omp simd aligned(tr_0_0_z_yzzz_xxx, tr_0_0_z_yzzz_xxy, tr_0_0_z_yzzz_xxz, tr_0_0_z_yzzz_xyy, tr_0_0_z_yzzz_xyz, tr_0_0_z_yzzz_xzz, tr_0_0_z_yzzz_yyy, tr_0_0_z_yzzz_yyz, tr_0_0_z_yzzz_yzz, tr_0_0_z_yzzz_zzz, tr_yzz_xxx, tr_yzz_xxy, tr_yzz_xxz, tr_yzz_xyy, tr_yzz_xyz, tr_yzz_xzz, tr_yzz_yyy, tr_yzz_yyz, tr_yzz_yzz, tr_yzz_zzz, tr_yzzz_xx, tr_yzzz_xxxz, tr_yzzz_xxyz, tr_yzzz_xxzz, tr_yzzz_xy, tr_yzzz_xyyz, tr_yzzz_xyzz, tr_yzzz_xz, tr_yzzz_xzzz, tr_yzzz_yy, tr_yzzz_yyyz, tr_yzzz_yyzz, tr_yzzz_yz, tr_yzzz_yzzz, tr_yzzz_zz, tr_yzzz_zzzz, tr_yzzzz_xxx, tr_yzzzz_xxy, tr_yzzzz_xxz, tr_yzzzz_xyy, tr_yzzzz_xyz, tr_yzzzz_xzz, tr_yzzzz_yyy, tr_yzzzz_yyz, tr_yzzzz_yzz, tr_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzzz_xxx[i] = 2.0 * tr_yzzzz_xxx[i] * tbe_0 + 2.0 * tr_yzzz_xxxz[i] * tke_0 - 3.0 * tr_yzz_xxx[i];

        tr_0_0_z_yzzz_xxy[i] = 2.0 * tr_yzzzz_xxy[i] * tbe_0 + 2.0 * tr_yzzz_xxyz[i] * tke_0 - 3.0 * tr_yzz_xxy[i];

        tr_0_0_z_yzzz_xxz[i] = 2.0 * tr_yzzzz_xxz[i] * tbe_0 + 2.0 * tr_yzzz_xxzz[i] * tke_0 - 3.0 * tr_yzz_xxz[i] - tr_yzzz_xx[i];

        tr_0_0_z_yzzz_xyy[i] = 2.0 * tr_yzzzz_xyy[i] * tbe_0 + 2.0 * tr_yzzz_xyyz[i] * tke_0 - 3.0 * tr_yzz_xyy[i];

        tr_0_0_z_yzzz_xyz[i] = 2.0 * tr_yzzzz_xyz[i] * tbe_0 + 2.0 * tr_yzzz_xyzz[i] * tke_0 - 3.0 * tr_yzz_xyz[i] - tr_yzzz_xy[i];

        tr_0_0_z_yzzz_xzz[i] = 2.0 * tr_yzzzz_xzz[i] * tbe_0 + 2.0 * tr_yzzz_xzzz[i] * tke_0 - 3.0 * tr_yzz_xzz[i] - 2.0 * tr_yzzz_xz[i];

        tr_0_0_z_yzzz_yyy[i] = 2.0 * tr_yzzzz_yyy[i] * tbe_0 + 2.0 * tr_yzzz_yyyz[i] * tke_0 - 3.0 * tr_yzz_yyy[i];

        tr_0_0_z_yzzz_yyz[i] = 2.0 * tr_yzzzz_yyz[i] * tbe_0 + 2.0 * tr_yzzz_yyzz[i] * tke_0 - 3.0 * tr_yzz_yyz[i] - tr_yzzz_yy[i];

        tr_0_0_z_yzzz_yzz[i] = 2.0 * tr_yzzzz_yzz[i] * tbe_0 + 2.0 * tr_yzzz_yzzz[i] * tke_0 - 3.0 * tr_yzz_yzz[i] - 2.0 * tr_yzzz_yz[i];

        tr_0_0_z_yzzz_zzz[i] = 2.0 * tr_yzzzz_zzz[i] * tbe_0 + 2.0 * tr_yzzz_zzzz[i] * tke_0 - 3.0 * tr_yzz_zzz[i] - 3.0 * tr_yzzz_zz[i];
    }

    // Set up 440-450 components of targeted buffer : GF

    auto tr_0_0_z_zzzz_xxx = pbuffer.data(idx_op_geom_010_gf + 440);

    auto tr_0_0_z_zzzz_xxy = pbuffer.data(idx_op_geom_010_gf + 441);

    auto tr_0_0_z_zzzz_xxz = pbuffer.data(idx_op_geom_010_gf + 442);

    auto tr_0_0_z_zzzz_xyy = pbuffer.data(idx_op_geom_010_gf + 443);

    auto tr_0_0_z_zzzz_xyz = pbuffer.data(idx_op_geom_010_gf + 444);

    auto tr_0_0_z_zzzz_xzz = pbuffer.data(idx_op_geom_010_gf + 445);

    auto tr_0_0_z_zzzz_yyy = pbuffer.data(idx_op_geom_010_gf + 446);

    auto tr_0_0_z_zzzz_yyz = pbuffer.data(idx_op_geom_010_gf + 447);

    auto tr_0_0_z_zzzz_yzz = pbuffer.data(idx_op_geom_010_gf + 448);

    auto tr_0_0_z_zzzz_zzz = pbuffer.data(idx_op_geom_010_gf + 449);

    #pragma omp simd aligned(tr_0_0_z_zzzz_xxx, tr_0_0_z_zzzz_xxy, tr_0_0_z_zzzz_xxz, tr_0_0_z_zzzz_xyy, tr_0_0_z_zzzz_xyz, tr_0_0_z_zzzz_xzz, tr_0_0_z_zzzz_yyy, tr_0_0_z_zzzz_yyz, tr_0_0_z_zzzz_yzz, tr_0_0_z_zzzz_zzz, tr_zzz_xxx, tr_zzz_xxy, tr_zzz_xxz, tr_zzz_xyy, tr_zzz_xyz, tr_zzz_xzz, tr_zzz_yyy, tr_zzz_yyz, tr_zzz_yzz, tr_zzz_zzz, tr_zzzz_xx, tr_zzzz_xxxz, tr_zzzz_xxyz, tr_zzzz_xxzz, tr_zzzz_xy, tr_zzzz_xyyz, tr_zzzz_xyzz, tr_zzzz_xz, tr_zzzz_xzzz, tr_zzzz_yy, tr_zzzz_yyyz, tr_zzzz_yyzz, tr_zzzz_yz, tr_zzzz_yzzz, tr_zzzz_zz, tr_zzzz_zzzz, tr_zzzzz_xxx, tr_zzzzz_xxy, tr_zzzzz_xxz, tr_zzzzz_xyy, tr_zzzzz_xyz, tr_zzzzz_xzz, tr_zzzzz_yyy, tr_zzzzz_yyz, tr_zzzzz_yzz, tr_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzzz_xxx[i] = 2.0 * tr_zzzzz_xxx[i] * tbe_0 + 2.0 * tr_zzzz_xxxz[i] * tke_0 - 4.0 * tr_zzz_xxx[i];

        tr_0_0_z_zzzz_xxy[i] = 2.0 * tr_zzzzz_xxy[i] * tbe_0 + 2.0 * tr_zzzz_xxyz[i] * tke_0 - 4.0 * tr_zzz_xxy[i];

        tr_0_0_z_zzzz_xxz[i] = 2.0 * tr_zzzzz_xxz[i] * tbe_0 + 2.0 * tr_zzzz_xxzz[i] * tke_0 - 4.0 * tr_zzz_xxz[i] - tr_zzzz_xx[i];

        tr_0_0_z_zzzz_xyy[i] = 2.0 * tr_zzzzz_xyy[i] * tbe_0 + 2.0 * tr_zzzz_xyyz[i] * tke_0 - 4.0 * tr_zzz_xyy[i];

        tr_0_0_z_zzzz_xyz[i] = 2.0 * tr_zzzzz_xyz[i] * tbe_0 + 2.0 * tr_zzzz_xyzz[i] * tke_0 - 4.0 * tr_zzz_xyz[i] - tr_zzzz_xy[i];

        tr_0_0_z_zzzz_xzz[i] = 2.0 * tr_zzzzz_xzz[i] * tbe_0 + 2.0 * tr_zzzz_xzzz[i] * tke_0 - 4.0 * tr_zzz_xzz[i] - 2.0 * tr_zzzz_xz[i];

        tr_0_0_z_zzzz_yyy[i] = 2.0 * tr_zzzzz_yyy[i] * tbe_0 + 2.0 * tr_zzzz_yyyz[i] * tke_0 - 4.0 * tr_zzz_yyy[i];

        tr_0_0_z_zzzz_yyz[i] = 2.0 * tr_zzzzz_yyz[i] * tbe_0 + 2.0 * tr_zzzz_yyzz[i] * tke_0 - 4.0 * tr_zzz_yyz[i] - tr_zzzz_yy[i];

        tr_0_0_z_zzzz_yzz[i] = 2.0 * tr_zzzzz_yzz[i] * tbe_0 + 2.0 * tr_zzzz_yzzz[i] * tke_0 - 4.0 * tr_zzz_yzz[i] - 2.0 * tr_zzzz_yz[i];

        tr_0_0_z_zzzz_zzz[i] = 2.0 * tr_zzzzz_zzz[i] * tbe_0 + 2.0 * tr_zzzz_zzzz[i] * tke_0 - 4.0 * tr_zzz_zzz[i] - 3.0 * tr_zzzz_zz[i];
    }

}

} // t2cgeom namespace

