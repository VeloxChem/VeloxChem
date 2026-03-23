#include "GeometricalDerivatives110ForGF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_gf(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_gf,
                         const int idx_op_df,
                         const int idx_op_fd,
                         const int idx_op_fg,
                         const int idx_op_gf,
                         const int idx_op_hd,
                         const int idx_op_hg,
                         const int idx_op_if,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up components of auxiliary buffer : HD

    auto tr_xxxxx_xx = pbuffer.data(idx_op_hd);

    auto tr_xxxxx_xy = pbuffer.data(idx_op_hd + 1);

    auto tr_xxxxx_xz = pbuffer.data(idx_op_hd + 2);

    auto tr_xxxxx_yy = pbuffer.data(idx_op_hd + 3);

    auto tr_xxxxx_yz = pbuffer.data(idx_op_hd + 4);

    auto tr_xxxxx_zz = pbuffer.data(idx_op_hd + 5);

    auto tr_xxxxy_xx = pbuffer.data(idx_op_hd + 6);

    auto tr_xxxxy_xy = pbuffer.data(idx_op_hd + 7);

    auto tr_xxxxy_xz = pbuffer.data(idx_op_hd + 8);

    auto tr_xxxxy_yy = pbuffer.data(idx_op_hd + 9);

    auto tr_xxxxy_yz = pbuffer.data(idx_op_hd + 10);

    auto tr_xxxxy_zz = pbuffer.data(idx_op_hd + 11);

    auto tr_xxxxz_xx = pbuffer.data(idx_op_hd + 12);

    auto tr_xxxxz_xy = pbuffer.data(idx_op_hd + 13);

    auto tr_xxxxz_xz = pbuffer.data(idx_op_hd + 14);

    auto tr_xxxxz_yy = pbuffer.data(idx_op_hd + 15);

    auto tr_xxxxz_yz = pbuffer.data(idx_op_hd + 16);

    auto tr_xxxxz_zz = pbuffer.data(idx_op_hd + 17);

    auto tr_xxxyy_xx = pbuffer.data(idx_op_hd + 18);

    auto tr_xxxyy_xy = pbuffer.data(idx_op_hd + 19);

    auto tr_xxxyy_xz = pbuffer.data(idx_op_hd + 20);

    auto tr_xxxyy_yy = pbuffer.data(idx_op_hd + 21);

    auto tr_xxxyy_yz = pbuffer.data(idx_op_hd + 22);

    auto tr_xxxyy_zz = pbuffer.data(idx_op_hd + 23);

    auto tr_xxxyz_xx = pbuffer.data(idx_op_hd + 24);

    auto tr_xxxyz_xy = pbuffer.data(idx_op_hd + 25);

    auto tr_xxxyz_xz = pbuffer.data(idx_op_hd + 26);

    auto tr_xxxyz_yy = pbuffer.data(idx_op_hd + 27);

    auto tr_xxxyz_yz = pbuffer.data(idx_op_hd + 28);

    auto tr_xxxyz_zz = pbuffer.data(idx_op_hd + 29);

    auto tr_xxxzz_xx = pbuffer.data(idx_op_hd + 30);

    auto tr_xxxzz_xy = pbuffer.data(idx_op_hd + 31);

    auto tr_xxxzz_xz = pbuffer.data(idx_op_hd + 32);

    auto tr_xxxzz_yy = pbuffer.data(idx_op_hd + 33);

    auto tr_xxxzz_yz = pbuffer.data(idx_op_hd + 34);

    auto tr_xxxzz_zz = pbuffer.data(idx_op_hd + 35);

    auto tr_xxyyy_xx = pbuffer.data(idx_op_hd + 36);

    auto tr_xxyyy_xy = pbuffer.data(idx_op_hd + 37);

    auto tr_xxyyy_xz = pbuffer.data(idx_op_hd + 38);

    auto tr_xxyyy_yy = pbuffer.data(idx_op_hd + 39);

    auto tr_xxyyy_yz = pbuffer.data(idx_op_hd + 40);

    auto tr_xxyyy_zz = pbuffer.data(idx_op_hd + 41);

    auto tr_xxyyz_xx = pbuffer.data(idx_op_hd + 42);

    auto tr_xxyyz_xy = pbuffer.data(idx_op_hd + 43);

    auto tr_xxyyz_xz = pbuffer.data(idx_op_hd + 44);

    auto tr_xxyyz_yy = pbuffer.data(idx_op_hd + 45);

    auto tr_xxyyz_yz = pbuffer.data(idx_op_hd + 46);

    auto tr_xxyyz_zz = pbuffer.data(idx_op_hd + 47);

    auto tr_xxyzz_xx = pbuffer.data(idx_op_hd + 48);

    auto tr_xxyzz_xy = pbuffer.data(idx_op_hd + 49);

    auto tr_xxyzz_xz = pbuffer.data(idx_op_hd + 50);

    auto tr_xxyzz_yy = pbuffer.data(idx_op_hd + 51);

    auto tr_xxyzz_yz = pbuffer.data(idx_op_hd + 52);

    auto tr_xxyzz_zz = pbuffer.data(idx_op_hd + 53);

    auto tr_xxzzz_xx = pbuffer.data(idx_op_hd + 54);

    auto tr_xxzzz_xy = pbuffer.data(idx_op_hd + 55);

    auto tr_xxzzz_xz = pbuffer.data(idx_op_hd + 56);

    auto tr_xxzzz_yy = pbuffer.data(idx_op_hd + 57);

    auto tr_xxzzz_yz = pbuffer.data(idx_op_hd + 58);

    auto tr_xxzzz_zz = pbuffer.data(idx_op_hd + 59);

    auto tr_xyyyy_xx = pbuffer.data(idx_op_hd + 60);

    auto tr_xyyyy_xy = pbuffer.data(idx_op_hd + 61);

    auto tr_xyyyy_xz = pbuffer.data(idx_op_hd + 62);

    auto tr_xyyyy_yy = pbuffer.data(idx_op_hd + 63);

    auto tr_xyyyy_yz = pbuffer.data(idx_op_hd + 64);

    auto tr_xyyyy_zz = pbuffer.data(idx_op_hd + 65);

    auto tr_xyyyz_xx = pbuffer.data(idx_op_hd + 66);

    auto tr_xyyyz_xy = pbuffer.data(idx_op_hd + 67);

    auto tr_xyyyz_xz = pbuffer.data(idx_op_hd + 68);

    auto tr_xyyyz_yy = pbuffer.data(idx_op_hd + 69);

    auto tr_xyyyz_yz = pbuffer.data(idx_op_hd + 70);

    auto tr_xyyyz_zz = pbuffer.data(idx_op_hd + 71);

    auto tr_xyyzz_xx = pbuffer.data(idx_op_hd + 72);

    auto tr_xyyzz_xy = pbuffer.data(idx_op_hd + 73);

    auto tr_xyyzz_xz = pbuffer.data(idx_op_hd + 74);

    auto tr_xyyzz_yy = pbuffer.data(idx_op_hd + 75);

    auto tr_xyyzz_yz = pbuffer.data(idx_op_hd + 76);

    auto tr_xyyzz_zz = pbuffer.data(idx_op_hd + 77);

    auto tr_xyzzz_xx = pbuffer.data(idx_op_hd + 78);

    auto tr_xyzzz_xy = pbuffer.data(idx_op_hd + 79);

    auto tr_xyzzz_xz = pbuffer.data(idx_op_hd + 80);

    auto tr_xyzzz_yy = pbuffer.data(idx_op_hd + 81);

    auto tr_xyzzz_yz = pbuffer.data(idx_op_hd + 82);

    auto tr_xyzzz_zz = pbuffer.data(idx_op_hd + 83);

    auto tr_xzzzz_xx = pbuffer.data(idx_op_hd + 84);

    auto tr_xzzzz_xy = pbuffer.data(idx_op_hd + 85);

    auto tr_xzzzz_xz = pbuffer.data(idx_op_hd + 86);

    auto tr_xzzzz_yy = pbuffer.data(idx_op_hd + 87);

    auto tr_xzzzz_yz = pbuffer.data(idx_op_hd + 88);

    auto tr_xzzzz_zz = pbuffer.data(idx_op_hd + 89);

    auto tr_yyyyy_xx = pbuffer.data(idx_op_hd + 90);

    auto tr_yyyyy_xy = pbuffer.data(idx_op_hd + 91);

    auto tr_yyyyy_xz = pbuffer.data(idx_op_hd + 92);

    auto tr_yyyyy_yy = pbuffer.data(idx_op_hd + 93);

    auto tr_yyyyy_yz = pbuffer.data(idx_op_hd + 94);

    auto tr_yyyyy_zz = pbuffer.data(idx_op_hd + 95);

    auto tr_yyyyz_xx = pbuffer.data(idx_op_hd + 96);

    auto tr_yyyyz_xy = pbuffer.data(idx_op_hd + 97);

    auto tr_yyyyz_xz = pbuffer.data(idx_op_hd + 98);

    auto tr_yyyyz_yy = pbuffer.data(idx_op_hd + 99);

    auto tr_yyyyz_yz = pbuffer.data(idx_op_hd + 100);

    auto tr_yyyyz_zz = pbuffer.data(idx_op_hd + 101);

    auto tr_yyyzz_xx = pbuffer.data(idx_op_hd + 102);

    auto tr_yyyzz_xy = pbuffer.data(idx_op_hd + 103);

    auto tr_yyyzz_xz = pbuffer.data(idx_op_hd + 104);

    auto tr_yyyzz_yy = pbuffer.data(idx_op_hd + 105);

    auto tr_yyyzz_yz = pbuffer.data(idx_op_hd + 106);

    auto tr_yyyzz_zz = pbuffer.data(idx_op_hd + 107);

    auto tr_yyzzz_xx = pbuffer.data(idx_op_hd + 108);

    auto tr_yyzzz_xy = pbuffer.data(idx_op_hd + 109);

    auto tr_yyzzz_xz = pbuffer.data(idx_op_hd + 110);

    auto tr_yyzzz_yy = pbuffer.data(idx_op_hd + 111);

    auto tr_yyzzz_yz = pbuffer.data(idx_op_hd + 112);

    auto tr_yyzzz_zz = pbuffer.data(idx_op_hd + 113);

    auto tr_yzzzz_xx = pbuffer.data(idx_op_hd + 114);

    auto tr_yzzzz_xy = pbuffer.data(idx_op_hd + 115);

    auto tr_yzzzz_xz = pbuffer.data(idx_op_hd + 116);

    auto tr_yzzzz_yy = pbuffer.data(idx_op_hd + 117);

    auto tr_yzzzz_yz = pbuffer.data(idx_op_hd + 118);

    auto tr_yzzzz_zz = pbuffer.data(idx_op_hd + 119);

    auto tr_zzzzz_xx = pbuffer.data(idx_op_hd + 120);

    auto tr_zzzzz_xy = pbuffer.data(idx_op_hd + 121);

    auto tr_zzzzz_xz = pbuffer.data(idx_op_hd + 122);

    auto tr_zzzzz_yy = pbuffer.data(idx_op_hd + 123);

    auto tr_zzzzz_yz = pbuffer.data(idx_op_hd + 124);

    auto tr_zzzzz_zz = pbuffer.data(idx_op_hd + 125);

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

    // Set up components of auxiliary buffer : IF

    auto tr_xxxxxx_xxx = pbuffer.data(idx_op_if);

    auto tr_xxxxxx_xxy = pbuffer.data(idx_op_if + 1);

    auto tr_xxxxxx_xxz = pbuffer.data(idx_op_if + 2);

    auto tr_xxxxxx_xyy = pbuffer.data(idx_op_if + 3);

    auto tr_xxxxxx_xyz = pbuffer.data(idx_op_if + 4);

    auto tr_xxxxxx_xzz = pbuffer.data(idx_op_if + 5);

    auto tr_xxxxxx_yyy = pbuffer.data(idx_op_if + 6);

    auto tr_xxxxxx_yyz = pbuffer.data(idx_op_if + 7);

    auto tr_xxxxxx_yzz = pbuffer.data(idx_op_if + 8);

    auto tr_xxxxxx_zzz = pbuffer.data(idx_op_if + 9);

    auto tr_xxxxxy_xxx = pbuffer.data(idx_op_if + 10);

    auto tr_xxxxxy_xxy = pbuffer.data(idx_op_if + 11);

    auto tr_xxxxxy_xxz = pbuffer.data(idx_op_if + 12);

    auto tr_xxxxxy_xyy = pbuffer.data(idx_op_if + 13);

    auto tr_xxxxxy_xyz = pbuffer.data(idx_op_if + 14);

    auto tr_xxxxxy_xzz = pbuffer.data(idx_op_if + 15);

    auto tr_xxxxxy_yyy = pbuffer.data(idx_op_if + 16);

    auto tr_xxxxxy_yyz = pbuffer.data(idx_op_if + 17);

    auto tr_xxxxxy_yzz = pbuffer.data(idx_op_if + 18);

    auto tr_xxxxxy_zzz = pbuffer.data(idx_op_if + 19);

    auto tr_xxxxxz_xxx = pbuffer.data(idx_op_if + 20);

    auto tr_xxxxxz_xxy = pbuffer.data(idx_op_if + 21);

    auto tr_xxxxxz_xxz = pbuffer.data(idx_op_if + 22);

    auto tr_xxxxxz_xyy = pbuffer.data(idx_op_if + 23);

    auto tr_xxxxxz_xyz = pbuffer.data(idx_op_if + 24);

    auto tr_xxxxxz_xzz = pbuffer.data(idx_op_if + 25);

    auto tr_xxxxxz_yyy = pbuffer.data(idx_op_if + 26);

    auto tr_xxxxxz_yyz = pbuffer.data(idx_op_if + 27);

    auto tr_xxxxxz_yzz = pbuffer.data(idx_op_if + 28);

    auto tr_xxxxxz_zzz = pbuffer.data(idx_op_if + 29);

    auto tr_xxxxyy_xxx = pbuffer.data(idx_op_if + 30);

    auto tr_xxxxyy_xxy = pbuffer.data(idx_op_if + 31);

    auto tr_xxxxyy_xxz = pbuffer.data(idx_op_if + 32);

    auto tr_xxxxyy_xyy = pbuffer.data(idx_op_if + 33);

    auto tr_xxxxyy_xyz = pbuffer.data(idx_op_if + 34);

    auto tr_xxxxyy_xzz = pbuffer.data(idx_op_if + 35);

    auto tr_xxxxyy_yyy = pbuffer.data(idx_op_if + 36);

    auto tr_xxxxyy_yyz = pbuffer.data(idx_op_if + 37);

    auto tr_xxxxyy_yzz = pbuffer.data(idx_op_if + 38);

    auto tr_xxxxyy_zzz = pbuffer.data(idx_op_if + 39);

    auto tr_xxxxyz_xxx = pbuffer.data(idx_op_if + 40);

    auto tr_xxxxyz_xxy = pbuffer.data(idx_op_if + 41);

    auto tr_xxxxyz_xxz = pbuffer.data(idx_op_if + 42);

    auto tr_xxxxyz_xyy = pbuffer.data(idx_op_if + 43);

    auto tr_xxxxyz_xyz = pbuffer.data(idx_op_if + 44);

    auto tr_xxxxyz_xzz = pbuffer.data(idx_op_if + 45);

    auto tr_xxxxyz_yyy = pbuffer.data(idx_op_if + 46);

    auto tr_xxxxyz_yyz = pbuffer.data(idx_op_if + 47);

    auto tr_xxxxyz_yzz = pbuffer.data(idx_op_if + 48);

    auto tr_xxxxyz_zzz = pbuffer.data(idx_op_if + 49);

    auto tr_xxxxzz_xxx = pbuffer.data(idx_op_if + 50);

    auto tr_xxxxzz_xxy = pbuffer.data(idx_op_if + 51);

    auto tr_xxxxzz_xxz = pbuffer.data(idx_op_if + 52);

    auto tr_xxxxzz_xyy = pbuffer.data(idx_op_if + 53);

    auto tr_xxxxzz_xyz = pbuffer.data(idx_op_if + 54);

    auto tr_xxxxzz_xzz = pbuffer.data(idx_op_if + 55);

    auto tr_xxxxzz_yyy = pbuffer.data(idx_op_if + 56);

    auto tr_xxxxzz_yyz = pbuffer.data(idx_op_if + 57);

    auto tr_xxxxzz_yzz = pbuffer.data(idx_op_if + 58);

    auto tr_xxxxzz_zzz = pbuffer.data(idx_op_if + 59);

    auto tr_xxxyyy_xxx = pbuffer.data(idx_op_if + 60);

    auto tr_xxxyyy_xxy = pbuffer.data(idx_op_if + 61);

    auto tr_xxxyyy_xxz = pbuffer.data(idx_op_if + 62);

    auto tr_xxxyyy_xyy = pbuffer.data(idx_op_if + 63);

    auto tr_xxxyyy_xyz = pbuffer.data(idx_op_if + 64);

    auto tr_xxxyyy_xzz = pbuffer.data(idx_op_if + 65);

    auto tr_xxxyyy_yyy = pbuffer.data(idx_op_if + 66);

    auto tr_xxxyyy_yyz = pbuffer.data(idx_op_if + 67);

    auto tr_xxxyyy_yzz = pbuffer.data(idx_op_if + 68);

    auto tr_xxxyyy_zzz = pbuffer.data(idx_op_if + 69);

    auto tr_xxxyyz_xxx = pbuffer.data(idx_op_if + 70);

    auto tr_xxxyyz_xxy = pbuffer.data(idx_op_if + 71);

    auto tr_xxxyyz_xxz = pbuffer.data(idx_op_if + 72);

    auto tr_xxxyyz_xyy = pbuffer.data(idx_op_if + 73);

    auto tr_xxxyyz_xyz = pbuffer.data(idx_op_if + 74);

    auto tr_xxxyyz_xzz = pbuffer.data(idx_op_if + 75);

    auto tr_xxxyyz_yyy = pbuffer.data(idx_op_if + 76);

    auto tr_xxxyyz_yyz = pbuffer.data(idx_op_if + 77);

    auto tr_xxxyyz_yzz = pbuffer.data(idx_op_if + 78);

    auto tr_xxxyyz_zzz = pbuffer.data(idx_op_if + 79);

    auto tr_xxxyzz_xxx = pbuffer.data(idx_op_if + 80);

    auto tr_xxxyzz_xxy = pbuffer.data(idx_op_if + 81);

    auto tr_xxxyzz_xxz = pbuffer.data(idx_op_if + 82);

    auto tr_xxxyzz_xyy = pbuffer.data(idx_op_if + 83);

    auto tr_xxxyzz_xyz = pbuffer.data(idx_op_if + 84);

    auto tr_xxxyzz_xzz = pbuffer.data(idx_op_if + 85);

    auto tr_xxxyzz_yyy = pbuffer.data(idx_op_if + 86);

    auto tr_xxxyzz_yyz = pbuffer.data(idx_op_if + 87);

    auto tr_xxxyzz_yzz = pbuffer.data(idx_op_if + 88);

    auto tr_xxxyzz_zzz = pbuffer.data(idx_op_if + 89);

    auto tr_xxxzzz_xxx = pbuffer.data(idx_op_if + 90);

    auto tr_xxxzzz_xxy = pbuffer.data(idx_op_if + 91);

    auto tr_xxxzzz_xxz = pbuffer.data(idx_op_if + 92);

    auto tr_xxxzzz_xyy = pbuffer.data(idx_op_if + 93);

    auto tr_xxxzzz_xyz = pbuffer.data(idx_op_if + 94);

    auto tr_xxxzzz_xzz = pbuffer.data(idx_op_if + 95);

    auto tr_xxxzzz_yyy = pbuffer.data(idx_op_if + 96);

    auto tr_xxxzzz_yyz = pbuffer.data(idx_op_if + 97);

    auto tr_xxxzzz_yzz = pbuffer.data(idx_op_if + 98);

    auto tr_xxxzzz_zzz = pbuffer.data(idx_op_if + 99);

    auto tr_xxyyyy_xxx = pbuffer.data(idx_op_if + 100);

    auto tr_xxyyyy_xxy = pbuffer.data(idx_op_if + 101);

    auto tr_xxyyyy_xxz = pbuffer.data(idx_op_if + 102);

    auto tr_xxyyyy_xyy = pbuffer.data(idx_op_if + 103);

    auto tr_xxyyyy_xyz = pbuffer.data(idx_op_if + 104);

    auto tr_xxyyyy_xzz = pbuffer.data(idx_op_if + 105);

    auto tr_xxyyyy_yyy = pbuffer.data(idx_op_if + 106);

    auto tr_xxyyyy_yyz = pbuffer.data(idx_op_if + 107);

    auto tr_xxyyyy_yzz = pbuffer.data(idx_op_if + 108);

    auto tr_xxyyyy_zzz = pbuffer.data(idx_op_if + 109);

    auto tr_xxyyyz_xxx = pbuffer.data(idx_op_if + 110);

    auto tr_xxyyyz_xxy = pbuffer.data(idx_op_if + 111);

    auto tr_xxyyyz_xxz = pbuffer.data(idx_op_if + 112);

    auto tr_xxyyyz_xyy = pbuffer.data(idx_op_if + 113);

    auto tr_xxyyyz_xyz = pbuffer.data(idx_op_if + 114);

    auto tr_xxyyyz_xzz = pbuffer.data(idx_op_if + 115);

    auto tr_xxyyyz_yyy = pbuffer.data(idx_op_if + 116);

    auto tr_xxyyyz_yyz = pbuffer.data(idx_op_if + 117);

    auto tr_xxyyyz_yzz = pbuffer.data(idx_op_if + 118);

    auto tr_xxyyyz_zzz = pbuffer.data(idx_op_if + 119);

    auto tr_xxyyzz_xxx = pbuffer.data(idx_op_if + 120);

    auto tr_xxyyzz_xxy = pbuffer.data(idx_op_if + 121);

    auto tr_xxyyzz_xxz = pbuffer.data(idx_op_if + 122);

    auto tr_xxyyzz_xyy = pbuffer.data(idx_op_if + 123);

    auto tr_xxyyzz_xyz = pbuffer.data(idx_op_if + 124);

    auto tr_xxyyzz_xzz = pbuffer.data(idx_op_if + 125);

    auto tr_xxyyzz_yyy = pbuffer.data(idx_op_if + 126);

    auto tr_xxyyzz_yyz = pbuffer.data(idx_op_if + 127);

    auto tr_xxyyzz_yzz = pbuffer.data(idx_op_if + 128);

    auto tr_xxyyzz_zzz = pbuffer.data(idx_op_if + 129);

    auto tr_xxyzzz_xxx = pbuffer.data(idx_op_if + 130);

    auto tr_xxyzzz_xxy = pbuffer.data(idx_op_if + 131);

    auto tr_xxyzzz_xxz = pbuffer.data(idx_op_if + 132);

    auto tr_xxyzzz_xyy = pbuffer.data(idx_op_if + 133);

    auto tr_xxyzzz_xyz = pbuffer.data(idx_op_if + 134);

    auto tr_xxyzzz_xzz = pbuffer.data(idx_op_if + 135);

    auto tr_xxyzzz_yyy = pbuffer.data(idx_op_if + 136);

    auto tr_xxyzzz_yyz = pbuffer.data(idx_op_if + 137);

    auto tr_xxyzzz_yzz = pbuffer.data(idx_op_if + 138);

    auto tr_xxyzzz_zzz = pbuffer.data(idx_op_if + 139);

    auto tr_xxzzzz_xxx = pbuffer.data(idx_op_if + 140);

    auto tr_xxzzzz_xxy = pbuffer.data(idx_op_if + 141);

    auto tr_xxzzzz_xxz = pbuffer.data(idx_op_if + 142);

    auto tr_xxzzzz_xyy = pbuffer.data(idx_op_if + 143);

    auto tr_xxzzzz_xyz = pbuffer.data(idx_op_if + 144);

    auto tr_xxzzzz_xzz = pbuffer.data(idx_op_if + 145);

    auto tr_xxzzzz_yyy = pbuffer.data(idx_op_if + 146);

    auto tr_xxzzzz_yyz = pbuffer.data(idx_op_if + 147);

    auto tr_xxzzzz_yzz = pbuffer.data(idx_op_if + 148);

    auto tr_xxzzzz_zzz = pbuffer.data(idx_op_if + 149);

    auto tr_xyyyyy_xxx = pbuffer.data(idx_op_if + 150);

    auto tr_xyyyyy_xxy = pbuffer.data(idx_op_if + 151);

    auto tr_xyyyyy_xxz = pbuffer.data(idx_op_if + 152);

    auto tr_xyyyyy_xyy = pbuffer.data(idx_op_if + 153);

    auto tr_xyyyyy_xyz = pbuffer.data(idx_op_if + 154);

    auto tr_xyyyyy_xzz = pbuffer.data(idx_op_if + 155);

    auto tr_xyyyyy_yyy = pbuffer.data(idx_op_if + 156);

    auto tr_xyyyyy_yyz = pbuffer.data(idx_op_if + 157);

    auto tr_xyyyyy_yzz = pbuffer.data(idx_op_if + 158);

    auto tr_xyyyyy_zzz = pbuffer.data(idx_op_if + 159);

    auto tr_xyyyyz_xxx = pbuffer.data(idx_op_if + 160);

    auto tr_xyyyyz_xxy = pbuffer.data(idx_op_if + 161);

    auto tr_xyyyyz_xxz = pbuffer.data(idx_op_if + 162);

    auto tr_xyyyyz_xyy = pbuffer.data(idx_op_if + 163);

    auto tr_xyyyyz_xyz = pbuffer.data(idx_op_if + 164);

    auto tr_xyyyyz_xzz = pbuffer.data(idx_op_if + 165);

    auto tr_xyyyyz_yyy = pbuffer.data(idx_op_if + 166);

    auto tr_xyyyyz_yyz = pbuffer.data(idx_op_if + 167);

    auto tr_xyyyyz_yzz = pbuffer.data(idx_op_if + 168);

    auto tr_xyyyyz_zzz = pbuffer.data(idx_op_if + 169);

    auto tr_xyyyzz_xxx = pbuffer.data(idx_op_if + 170);

    auto tr_xyyyzz_xxy = pbuffer.data(idx_op_if + 171);

    auto tr_xyyyzz_xxz = pbuffer.data(idx_op_if + 172);

    auto tr_xyyyzz_xyy = pbuffer.data(idx_op_if + 173);

    auto tr_xyyyzz_xyz = pbuffer.data(idx_op_if + 174);

    auto tr_xyyyzz_xzz = pbuffer.data(idx_op_if + 175);

    auto tr_xyyyzz_yyy = pbuffer.data(idx_op_if + 176);

    auto tr_xyyyzz_yyz = pbuffer.data(idx_op_if + 177);

    auto tr_xyyyzz_yzz = pbuffer.data(idx_op_if + 178);

    auto tr_xyyyzz_zzz = pbuffer.data(idx_op_if + 179);

    auto tr_xyyzzz_xxx = pbuffer.data(idx_op_if + 180);

    auto tr_xyyzzz_xxy = pbuffer.data(idx_op_if + 181);

    auto tr_xyyzzz_xxz = pbuffer.data(idx_op_if + 182);

    auto tr_xyyzzz_xyy = pbuffer.data(idx_op_if + 183);

    auto tr_xyyzzz_xyz = pbuffer.data(idx_op_if + 184);

    auto tr_xyyzzz_xzz = pbuffer.data(idx_op_if + 185);

    auto tr_xyyzzz_yyy = pbuffer.data(idx_op_if + 186);

    auto tr_xyyzzz_yyz = pbuffer.data(idx_op_if + 187);

    auto tr_xyyzzz_yzz = pbuffer.data(idx_op_if + 188);

    auto tr_xyyzzz_zzz = pbuffer.data(idx_op_if + 189);

    auto tr_xyzzzz_xxx = pbuffer.data(idx_op_if + 190);

    auto tr_xyzzzz_xxy = pbuffer.data(idx_op_if + 191);

    auto tr_xyzzzz_xxz = pbuffer.data(idx_op_if + 192);

    auto tr_xyzzzz_xyy = pbuffer.data(idx_op_if + 193);

    auto tr_xyzzzz_xyz = pbuffer.data(idx_op_if + 194);

    auto tr_xyzzzz_xzz = pbuffer.data(idx_op_if + 195);

    auto tr_xyzzzz_yyy = pbuffer.data(idx_op_if + 196);

    auto tr_xyzzzz_yyz = pbuffer.data(idx_op_if + 197);

    auto tr_xyzzzz_yzz = pbuffer.data(idx_op_if + 198);

    auto tr_xyzzzz_zzz = pbuffer.data(idx_op_if + 199);

    auto tr_xzzzzz_xxx = pbuffer.data(idx_op_if + 200);

    auto tr_xzzzzz_xxy = pbuffer.data(idx_op_if + 201);

    auto tr_xzzzzz_xxz = pbuffer.data(idx_op_if + 202);

    auto tr_xzzzzz_xyy = pbuffer.data(idx_op_if + 203);

    auto tr_xzzzzz_xyz = pbuffer.data(idx_op_if + 204);

    auto tr_xzzzzz_xzz = pbuffer.data(idx_op_if + 205);

    auto tr_xzzzzz_yyy = pbuffer.data(idx_op_if + 206);

    auto tr_xzzzzz_yyz = pbuffer.data(idx_op_if + 207);

    auto tr_xzzzzz_yzz = pbuffer.data(idx_op_if + 208);

    auto tr_xzzzzz_zzz = pbuffer.data(idx_op_if + 209);

    auto tr_yyyyyy_xxx = pbuffer.data(idx_op_if + 210);

    auto tr_yyyyyy_xxy = pbuffer.data(idx_op_if + 211);

    auto tr_yyyyyy_xxz = pbuffer.data(idx_op_if + 212);

    auto tr_yyyyyy_xyy = pbuffer.data(idx_op_if + 213);

    auto tr_yyyyyy_xyz = pbuffer.data(idx_op_if + 214);

    auto tr_yyyyyy_xzz = pbuffer.data(idx_op_if + 215);

    auto tr_yyyyyy_yyy = pbuffer.data(idx_op_if + 216);

    auto tr_yyyyyy_yyz = pbuffer.data(idx_op_if + 217);

    auto tr_yyyyyy_yzz = pbuffer.data(idx_op_if + 218);

    auto tr_yyyyyy_zzz = pbuffer.data(idx_op_if + 219);

    auto tr_yyyyyz_xxx = pbuffer.data(idx_op_if + 220);

    auto tr_yyyyyz_xxy = pbuffer.data(idx_op_if + 221);

    auto tr_yyyyyz_xxz = pbuffer.data(idx_op_if + 222);

    auto tr_yyyyyz_xyy = pbuffer.data(idx_op_if + 223);

    auto tr_yyyyyz_xyz = pbuffer.data(idx_op_if + 224);

    auto tr_yyyyyz_xzz = pbuffer.data(idx_op_if + 225);

    auto tr_yyyyyz_yyy = pbuffer.data(idx_op_if + 226);

    auto tr_yyyyyz_yyz = pbuffer.data(idx_op_if + 227);

    auto tr_yyyyyz_yzz = pbuffer.data(idx_op_if + 228);

    auto tr_yyyyyz_zzz = pbuffer.data(idx_op_if + 229);

    auto tr_yyyyzz_xxx = pbuffer.data(idx_op_if + 230);

    auto tr_yyyyzz_xxy = pbuffer.data(idx_op_if + 231);

    auto tr_yyyyzz_xxz = pbuffer.data(idx_op_if + 232);

    auto tr_yyyyzz_xyy = pbuffer.data(idx_op_if + 233);

    auto tr_yyyyzz_xyz = pbuffer.data(idx_op_if + 234);

    auto tr_yyyyzz_xzz = pbuffer.data(idx_op_if + 235);

    auto tr_yyyyzz_yyy = pbuffer.data(idx_op_if + 236);

    auto tr_yyyyzz_yyz = pbuffer.data(idx_op_if + 237);

    auto tr_yyyyzz_yzz = pbuffer.data(idx_op_if + 238);

    auto tr_yyyyzz_zzz = pbuffer.data(idx_op_if + 239);

    auto tr_yyyzzz_xxx = pbuffer.data(idx_op_if + 240);

    auto tr_yyyzzz_xxy = pbuffer.data(idx_op_if + 241);

    auto tr_yyyzzz_xxz = pbuffer.data(idx_op_if + 242);

    auto tr_yyyzzz_xyy = pbuffer.data(idx_op_if + 243);

    auto tr_yyyzzz_xyz = pbuffer.data(idx_op_if + 244);

    auto tr_yyyzzz_xzz = pbuffer.data(idx_op_if + 245);

    auto tr_yyyzzz_yyy = pbuffer.data(idx_op_if + 246);

    auto tr_yyyzzz_yyz = pbuffer.data(idx_op_if + 247);

    auto tr_yyyzzz_yzz = pbuffer.data(idx_op_if + 248);

    auto tr_yyyzzz_zzz = pbuffer.data(idx_op_if + 249);

    auto tr_yyzzzz_xxx = pbuffer.data(idx_op_if + 250);

    auto tr_yyzzzz_xxy = pbuffer.data(idx_op_if + 251);

    auto tr_yyzzzz_xxz = pbuffer.data(idx_op_if + 252);

    auto tr_yyzzzz_xyy = pbuffer.data(idx_op_if + 253);

    auto tr_yyzzzz_xyz = pbuffer.data(idx_op_if + 254);

    auto tr_yyzzzz_xzz = pbuffer.data(idx_op_if + 255);

    auto tr_yyzzzz_yyy = pbuffer.data(idx_op_if + 256);

    auto tr_yyzzzz_yyz = pbuffer.data(idx_op_if + 257);

    auto tr_yyzzzz_yzz = pbuffer.data(idx_op_if + 258);

    auto tr_yyzzzz_zzz = pbuffer.data(idx_op_if + 259);

    auto tr_yzzzzz_xxx = pbuffer.data(idx_op_if + 260);

    auto tr_yzzzzz_xxy = pbuffer.data(idx_op_if + 261);

    auto tr_yzzzzz_xxz = pbuffer.data(idx_op_if + 262);

    auto tr_yzzzzz_xyy = pbuffer.data(idx_op_if + 263);

    auto tr_yzzzzz_xyz = pbuffer.data(idx_op_if + 264);

    auto tr_yzzzzz_xzz = pbuffer.data(idx_op_if + 265);

    auto tr_yzzzzz_yyy = pbuffer.data(idx_op_if + 266);

    auto tr_yzzzzz_yyz = pbuffer.data(idx_op_if + 267);

    auto tr_yzzzzz_yzz = pbuffer.data(idx_op_if + 268);

    auto tr_yzzzzz_zzz = pbuffer.data(idx_op_if + 269);

    auto tr_zzzzzz_xxx = pbuffer.data(idx_op_if + 270);

    auto tr_zzzzzz_xxy = pbuffer.data(idx_op_if + 271);

    auto tr_zzzzzz_xxz = pbuffer.data(idx_op_if + 272);

    auto tr_zzzzzz_xyy = pbuffer.data(idx_op_if + 273);

    auto tr_zzzzzz_xyz = pbuffer.data(idx_op_if + 274);

    auto tr_zzzzzz_xzz = pbuffer.data(idx_op_if + 275);

    auto tr_zzzzzz_yyy = pbuffer.data(idx_op_if + 276);

    auto tr_zzzzzz_yyz = pbuffer.data(idx_op_if + 277);

    auto tr_zzzzzz_yzz = pbuffer.data(idx_op_if + 278);

    auto tr_zzzzzz_zzz = pbuffer.data(idx_op_if + 279);

    // Set up 0-10 components of targeted buffer : GF

    auto tr_x_0_x_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf);

    auto tr_x_0_x_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 1);

    auto tr_x_0_x_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 2);

    auto tr_x_0_x_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 3);

    auto tr_x_0_x_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 4);

    auto tr_x_0_x_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 5);

    auto tr_x_0_x_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 6);

    auto tr_x_0_x_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 7);

    auto tr_x_0_x_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 8);

    auto tr_x_0_x_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 9);

    #pragma omp simd aligned(tr_x_0_x_xxxx_xxx, tr_x_0_x_xxxx_xxy, tr_x_0_x_xxxx_xxz, tr_x_0_x_xxxx_xyy, tr_x_0_x_xxxx_xyz, tr_x_0_x_xxxx_xzz, tr_x_0_x_xxxx_yyy, tr_x_0_x_xxxx_yyz, tr_x_0_x_xxxx_yzz, tr_x_0_x_xxxx_zzz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xx, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxx_xx, tr_xxxxx_xxxx, tr_xxxxx_xxxy, tr_xxxxx_xxxz, tr_xxxxx_xxyy, tr_xxxxx_xxyz, tr_xxxxx_xxzz, tr_xxxxx_xy, tr_xxxxx_xyyy, tr_xxxxx_xyyz, tr_xxxxx_xyzz, tr_xxxxx_xz, tr_xxxxx_xzzz, tr_xxxxx_yy, tr_xxxxx_yz, tr_xxxxx_zz, tr_xxxxxx_xxx, tr_xxxxxx_xxy, tr_xxxxxx_xxz, tr_xxxxxx_xyy, tr_xxxxxx_xyz, tr_xxxxxx_xzz, tr_xxxxxx_yyy, tr_xxxxxx_yyz, tr_xxxxxx_yzz, tr_xxxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxxx_xxx[i] = 12.0 * tr_xx_xxx[i] + 12.0 * tr_xxx_xx[i] - 8.0 * tr_xxx_xxxx[i] * tke_0 - 18.0 * tr_xxxx_xxx[i] * tbe_0 - 6.0 * tr_xxxxx_xx[i] * tbe_0 + 4.0 * tr_xxxxx_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_xxy[i] = 12.0 * tr_xx_xxy[i] + 8.0 * tr_xxx_xy[i] - 8.0 * tr_xxx_xxxy[i] * tke_0 - 18.0 * tr_xxxx_xxy[i] * tbe_0 - 4.0 * tr_xxxxx_xy[i] * tbe_0 + 4.0 * tr_xxxxx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_xxz[i] = 12.0 * tr_xx_xxz[i] + 8.0 * tr_xxx_xz[i] - 8.0 * tr_xxx_xxxz[i] * tke_0 - 18.0 * tr_xxxx_xxz[i] * tbe_0 - 4.0 * tr_xxxxx_xz[i] * tbe_0 + 4.0 * tr_xxxxx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_xyy[i] = 12.0 * tr_xx_xyy[i] + 4.0 * tr_xxx_yy[i] - 8.0 * tr_xxx_xxyy[i] * tke_0 - 18.0 * tr_xxxx_xyy[i] * tbe_0 - 2.0 * tr_xxxxx_yy[i] * tbe_0 + 4.0 * tr_xxxxx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_xyz[i] = 12.0 * tr_xx_xyz[i] + 4.0 * tr_xxx_yz[i] - 8.0 * tr_xxx_xxyz[i] * tke_0 - 18.0 * tr_xxxx_xyz[i] * tbe_0 - 2.0 * tr_xxxxx_yz[i] * tbe_0 + 4.0 * tr_xxxxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_xzz[i] = 12.0 * tr_xx_xzz[i] + 4.0 * tr_xxx_zz[i] - 8.0 * tr_xxx_xxzz[i] * tke_0 - 18.0 * tr_xxxx_xzz[i] * tbe_0 - 2.0 * tr_xxxxx_zz[i] * tbe_0 + 4.0 * tr_xxxxx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_yyy[i] = 12.0 * tr_xx_yyy[i] - 8.0 * tr_xxx_xyyy[i] * tke_0 - 18.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_yyz[i] = 12.0 * tr_xx_yyz[i] - 8.0 * tr_xxx_xyyz[i] * tke_0 - 18.0 * tr_xxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_yzz[i] = 12.0 * tr_xx_yzz[i] - 8.0 * tr_xxx_xyzz[i] * tke_0 - 18.0 * tr_xxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxx_zzz[i] = 12.0 * tr_xx_zzz[i] - 8.0 * tr_xxx_xzzz[i] * tke_0 - 18.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxx_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 10-20 components of targeted buffer : GF

    auto tr_x_0_x_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 10);

    auto tr_x_0_x_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 11);

    auto tr_x_0_x_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 12);

    auto tr_x_0_x_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 13);

    auto tr_x_0_x_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 14);

    auto tr_x_0_x_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 15);

    auto tr_x_0_x_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 16);

    auto tr_x_0_x_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 17);

    auto tr_x_0_x_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 18);

    auto tr_x_0_x_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 19);

    #pragma omp simd aligned(tr_x_0_x_xxxy_xxx, tr_x_0_x_xxxy_xxy, tr_x_0_x_xxxy_xxz, tr_x_0_x_xxxy_xyy, tr_x_0_x_xxxy_xyz, tr_x_0_x_xxxy_xzz, tr_x_0_x_xxxy_yyy, tr_x_0_x_xxxy_yyz, tr_x_0_x_xxxy_yzz, tr_x_0_x_xxxy_zzz, tr_xxxxxy_xxx, tr_xxxxxy_xxy, tr_xxxxxy_xxz, tr_xxxxxy_xyy, tr_xxxxxy_xyz, tr_xxxxxy_xzz, tr_xxxxxy_yyy, tr_xxxxxy_yyz, tr_xxxxxy_yzz, tr_xxxxxy_zzz, tr_xxxxy_xx, tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xy, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xz, tr_xxxxy_xzzz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xx, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxxy_xxx[i] = 6.0 * tr_xy_xxx[i] + 9.0 * tr_xxy_xx[i] - 6.0 * tr_xxy_xxxx[i] * tke_0 - 14.0 * tr_xxxy_xxx[i] * tbe_0 - 6.0 * tr_xxxxy_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_xxy[i] = 6.0 * tr_xy_xxy[i] + 6.0 * tr_xxy_xy[i] - 6.0 * tr_xxy_xxxy[i] * tke_0 - 14.0 * tr_xxxy_xxy[i] * tbe_0 - 4.0 * tr_xxxxy_xy[i] * tbe_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_xxz[i] = 6.0 * tr_xy_xxz[i] + 6.0 * tr_xxy_xz[i] - 6.0 * tr_xxy_xxxz[i] * tke_0 - 14.0 * tr_xxxy_xxz[i] * tbe_0 - 4.0 * tr_xxxxy_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_xyy[i] = 6.0 * tr_xy_xyy[i] + 3.0 * tr_xxy_yy[i] - 6.0 * tr_xxy_xxyy[i] * tke_0 - 14.0 * tr_xxxy_xyy[i] * tbe_0 - 2.0 * tr_xxxxy_yy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_xyz[i] = 6.0 * tr_xy_xyz[i] + 3.0 * tr_xxy_yz[i] - 6.0 * tr_xxy_xxyz[i] * tke_0 - 14.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxxy_yz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_xzz[i] = 6.0 * tr_xy_xzz[i] + 3.0 * tr_xxy_zz[i] - 6.0 * tr_xxy_xxzz[i] * tke_0 - 14.0 * tr_xxxy_xzz[i] * tbe_0 - 2.0 * tr_xxxxy_zz[i] * tbe_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_yyy[i] = 6.0 * tr_xy_yyy[i] - 6.0 * tr_xxy_xyyy[i] * tke_0 - 14.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_yyz[i] = 6.0 * tr_xy_yyz[i] - 6.0 * tr_xxy_xyyz[i] * tke_0 - 14.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_yzz[i] = 6.0 * tr_xy_yzz[i] - 6.0 * tr_xxy_xyzz[i] * tke_0 - 14.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxy_zzz[i] = 6.0 * tr_xy_zzz[i] - 6.0 * tr_xxy_xzzz[i] * tke_0 - 14.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 20-30 components of targeted buffer : GF

    auto tr_x_0_x_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 20);

    auto tr_x_0_x_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 21);

    auto tr_x_0_x_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 22);

    auto tr_x_0_x_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 23);

    auto tr_x_0_x_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 24);

    auto tr_x_0_x_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 25);

    auto tr_x_0_x_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 26);

    auto tr_x_0_x_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 27);

    auto tr_x_0_x_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 28);

    auto tr_x_0_x_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 29);

    #pragma omp simd aligned(tr_x_0_x_xxxz_xxx, tr_x_0_x_xxxz_xxy, tr_x_0_x_xxxz_xxz, tr_x_0_x_xxxz_xyy, tr_x_0_x_xxxz_xyz, tr_x_0_x_xxxz_xzz, tr_x_0_x_xxxz_yyy, tr_x_0_x_xxxz_yyz, tr_x_0_x_xxxz_yzz, tr_x_0_x_xxxz_zzz, tr_xxxxxz_xxx, tr_xxxxxz_xxy, tr_xxxxxz_xxz, tr_xxxxxz_xyy, tr_xxxxxz_xyz, tr_xxxxxz_xzz, tr_xxxxxz_yyy, tr_xxxxxz_yyz, tr_xxxxxz_yzz, tr_xxxxxz_zzz, tr_xxxxz_xx, tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xy, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xz, tr_xxxxz_xzzz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xx, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxxz_xxx[i] = 6.0 * tr_xz_xxx[i] + 9.0 * tr_xxz_xx[i] - 6.0 * tr_xxz_xxxx[i] * tke_0 - 14.0 * tr_xxxz_xxx[i] * tbe_0 - 6.0 * tr_xxxxz_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_xxy[i] = 6.0 * tr_xz_xxy[i] + 6.0 * tr_xxz_xy[i] - 6.0 * tr_xxz_xxxy[i] * tke_0 - 14.0 * tr_xxxz_xxy[i] * tbe_0 - 4.0 * tr_xxxxz_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_xxz[i] = 6.0 * tr_xz_xxz[i] + 6.0 * tr_xxz_xz[i] - 6.0 * tr_xxz_xxxz[i] * tke_0 - 14.0 * tr_xxxz_xxz[i] * tbe_0 - 4.0 * tr_xxxxz_xz[i] * tbe_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_xyy[i] = 6.0 * tr_xz_xyy[i] + 3.0 * tr_xxz_yy[i] - 6.0 * tr_xxz_xxyy[i] * tke_0 - 14.0 * tr_xxxz_xyy[i] * tbe_0 - 2.0 * tr_xxxxz_yy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_xyz[i] = 6.0 * tr_xz_xyz[i] + 3.0 * tr_xxz_yz[i] - 6.0 * tr_xxz_xxyz[i] * tke_0 - 14.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxxz_yz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_xzz[i] = 6.0 * tr_xz_xzz[i] + 3.0 * tr_xxz_zz[i] - 6.0 * tr_xxz_xxzz[i] * tke_0 - 14.0 * tr_xxxz_xzz[i] * tbe_0 - 2.0 * tr_xxxxz_zz[i] * tbe_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_yyy[i] = 6.0 * tr_xz_yyy[i] - 6.0 * tr_xxz_xyyy[i] * tke_0 - 14.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_yyz[i] = 6.0 * tr_xz_yyz[i] - 6.0 * tr_xxz_xyyz[i] * tke_0 - 14.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_yzz[i] = 6.0 * tr_xz_yzz[i] - 6.0 * tr_xxz_xyzz[i] * tke_0 - 14.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxxz_zzz[i] = 6.0 * tr_xz_zzz[i] - 6.0 * tr_xxz_xzzz[i] * tke_0 - 14.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-40 components of targeted buffer : GF

    auto tr_x_0_x_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 30);

    auto tr_x_0_x_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 31);

    auto tr_x_0_x_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 32);

    auto tr_x_0_x_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 33);

    auto tr_x_0_x_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 34);

    auto tr_x_0_x_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 35);

    auto tr_x_0_x_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 36);

    auto tr_x_0_x_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 37);

    auto tr_x_0_x_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 38);

    auto tr_x_0_x_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 39);

    #pragma omp simd aligned(tr_x_0_x_xxyy_xxx, tr_x_0_x_xxyy_xxy, tr_x_0_x_xxyy_xxz, tr_x_0_x_xxyy_xyy, tr_x_0_x_xxyy_xyz, tr_x_0_x_xxyy_xzz, tr_x_0_x_xxyy_yyy, tr_x_0_x_xxyy_yyz, tr_x_0_x_xxyy_yzz, tr_x_0_x_xxyy_zzz, tr_xxxxyy_xxx, tr_xxxxyy_xxy, tr_xxxxyy_xxz, tr_xxxxyy_xyy, tr_xxxxyy_xyz, tr_xxxxyy_xzz, tr_xxxxyy_yyy, tr_xxxxyy_yyz, tr_xxxxyy_yzz, tr_xxxxyy_zzz, tr_xxxyy_xx, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xy, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xz, tr_xxxyy_xzzz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xyy_xx, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxyy_xxx[i] = 2.0 * tr_yy_xxx[i] + 6.0 * tr_xyy_xx[i] - 4.0 * tr_xyy_xxxx[i] * tke_0 - 10.0 * tr_xxyy_xxx[i] * tbe_0 - 6.0 * tr_xxxyy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_xxy[i] = 2.0 * tr_yy_xxy[i] + 4.0 * tr_xyy_xy[i] - 4.0 * tr_xyy_xxxy[i] * tke_0 - 10.0 * tr_xxyy_xxy[i] * tbe_0 - 4.0 * tr_xxxyy_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_xxz[i] = 2.0 * tr_yy_xxz[i] + 4.0 * tr_xyy_xz[i] - 4.0 * tr_xyy_xxxz[i] * tke_0 - 10.0 * tr_xxyy_xxz[i] * tbe_0 - 4.0 * tr_xxxyy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_xyy[i] = 2.0 * tr_yy_xyy[i] + 2.0 * tr_xyy_yy[i] - 4.0 * tr_xyy_xxyy[i] * tke_0 - 10.0 * tr_xxyy_xyy[i] * tbe_0 - 2.0 * tr_xxxyy_yy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_xyz[i] = 2.0 * tr_yy_xyz[i] + 2.0 * tr_xyy_yz[i] - 4.0 * tr_xyy_xxyz[i] * tke_0 - 10.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxxyy_yz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_xzz[i] = 2.0 * tr_yy_xzz[i] + 2.0 * tr_xyy_zz[i] - 4.0 * tr_xyy_xxzz[i] * tke_0 - 10.0 * tr_xxyy_xzz[i] * tbe_0 - 2.0 * tr_xxxyy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_yyy[i] = 2.0 * tr_yy_yyy[i] - 4.0 * tr_xyy_xyyy[i] * tke_0 - 10.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_yyz[i] = 2.0 * tr_yy_yyz[i] - 4.0 * tr_xyy_xyyz[i] * tke_0 - 10.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_yzz[i] = 2.0 * tr_yy_yzz[i] - 4.0 * tr_xyy_xyzz[i] * tke_0 - 10.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyy_zzz[i] = 2.0 * tr_yy_zzz[i] - 4.0 * tr_xyy_xzzz[i] * tke_0 - 10.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 40-50 components of targeted buffer : GF

    auto tr_x_0_x_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 40);

    auto tr_x_0_x_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 41);

    auto tr_x_0_x_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 42);

    auto tr_x_0_x_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 43);

    auto tr_x_0_x_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 44);

    auto tr_x_0_x_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 45);

    auto tr_x_0_x_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 46);

    auto tr_x_0_x_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 47);

    auto tr_x_0_x_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 48);

    auto tr_x_0_x_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 49);

    #pragma omp simd aligned(tr_x_0_x_xxyz_xxx, tr_x_0_x_xxyz_xxy, tr_x_0_x_xxyz_xxz, tr_x_0_x_xxyz_xyy, tr_x_0_x_xxyz_xyz, tr_x_0_x_xxyz_xzz, tr_x_0_x_xxyz_yyy, tr_x_0_x_xxyz_yyz, tr_x_0_x_xxyz_yzz, tr_x_0_x_xxyz_zzz, tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_xxxyz_xx, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xy, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_xzzz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxyz_xxx[i] = 2.0 * tr_yz_xxx[i] + 6.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxxx[i] * tke_0 - 10.0 * tr_xxyz_xxx[i] * tbe_0 - 6.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_xxy[i] = 2.0 * tr_yz_xxy[i] + 4.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xxxy[i] * tke_0 - 10.0 * tr_xxyz_xxy[i] * tbe_0 - 4.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_xxz[i] = 2.0 * tr_yz_xxz[i] + 4.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xxxz[i] * tke_0 - 10.0 * tr_xxyz_xxz[i] * tbe_0 - 4.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_xyy[i] = 2.0 * tr_yz_xyy[i] + 2.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_xxyy[i] * tke_0 - 10.0 * tr_xxyz_xyy[i] * tbe_0 - 2.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_xyz[i] = 2.0 * tr_yz_xyz[i] + 2.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 10.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_xzz[i] = 2.0 * tr_yz_xzz[i] + 2.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_xxzz[i] * tke_0 - 10.0 * tr_xxyz_xzz[i] * tbe_0 - 2.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_yyy[i] = 2.0 * tr_yz_yyy[i] - 4.0 * tr_xyz_xyyy[i] * tke_0 - 10.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_yyz[i] = 2.0 * tr_yz_yyz[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 10.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_yzz[i] = 2.0 * tr_yz_yzz[i] - 4.0 * tr_xyz_xyzz[i] * tke_0 - 10.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxyz_zzz[i] = 2.0 * tr_yz_zzz[i] - 4.0 * tr_xyz_xzzz[i] * tke_0 - 10.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 50-60 components of targeted buffer : GF

    auto tr_x_0_x_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 50);

    auto tr_x_0_x_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 51);

    auto tr_x_0_x_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 52);

    auto tr_x_0_x_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 53);

    auto tr_x_0_x_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 54);

    auto tr_x_0_x_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 55);

    auto tr_x_0_x_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 56);

    auto tr_x_0_x_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 57);

    auto tr_x_0_x_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 58);

    auto tr_x_0_x_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 59);

    #pragma omp simd aligned(tr_x_0_x_xxzz_xxx, tr_x_0_x_xxzz_xxy, tr_x_0_x_xxzz_xxz, tr_x_0_x_xxzz_xyy, tr_x_0_x_xxzz_xyz, tr_x_0_x_xxzz_xzz, tr_x_0_x_xxzz_yyy, tr_x_0_x_xxzz_yyz, tr_x_0_x_xxzz_yzz, tr_x_0_x_xxzz_zzz, tr_xxxxzz_xxx, tr_xxxxzz_xxy, tr_xxxxzz_xxz, tr_xxxxzz_xyy, tr_xxxxzz_xyz, tr_xxxxzz_xzz, tr_xxxxzz_yyy, tr_xxxxzz_yyz, tr_xxxxzz_yzz, tr_xxxxzz_zzz, tr_xxxzz_xx, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xy, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xz, tr_xxxzz_xzzz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xzz_xx, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxzz_xxx[i] = 2.0 * tr_zz_xxx[i] + 6.0 * tr_xzz_xx[i] - 4.0 * tr_xzz_xxxx[i] * tke_0 - 10.0 * tr_xxzz_xxx[i] * tbe_0 - 6.0 * tr_xxxzz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_xxy[i] = 2.0 * tr_zz_xxy[i] + 4.0 * tr_xzz_xy[i] - 4.0 * tr_xzz_xxxy[i] * tke_0 - 10.0 * tr_xxzz_xxy[i] * tbe_0 - 4.0 * tr_xxxzz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_xxz[i] = 2.0 * tr_zz_xxz[i] + 4.0 * tr_xzz_xz[i] - 4.0 * tr_xzz_xxxz[i] * tke_0 - 10.0 * tr_xxzz_xxz[i] * tbe_0 - 4.0 * tr_xxxzz_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_xyy[i] = 2.0 * tr_zz_xyy[i] + 2.0 * tr_xzz_yy[i] - 4.0 * tr_xzz_xxyy[i] * tke_0 - 10.0 * tr_xxzz_xyy[i] * tbe_0 - 2.0 * tr_xxxzz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_xyz[i] = 2.0 * tr_zz_xyz[i] + 2.0 * tr_xzz_yz[i] - 4.0 * tr_xzz_xxyz[i] * tke_0 - 10.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxxzz_yz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_xzz[i] = 2.0 * tr_zz_xzz[i] + 2.0 * tr_xzz_zz[i] - 4.0 * tr_xzz_xxzz[i] * tke_0 - 10.0 * tr_xxzz_xzz[i] * tbe_0 - 2.0 * tr_xxxzz_zz[i] * tbe_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_yyy[i] = 2.0 * tr_zz_yyy[i] - 4.0 * tr_xzz_xyyy[i] * tke_0 - 10.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_yyz[i] = 2.0 * tr_zz_yyz[i] - 4.0 * tr_xzz_xyyz[i] * tke_0 - 10.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_yzz[i] = 2.0 * tr_zz_yzz[i] - 4.0 * tr_xzz_xyzz[i] * tke_0 - 10.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxzz_zzz[i] = 2.0 * tr_zz_zzz[i] - 4.0 * tr_xzz_xzzz[i] * tke_0 - 10.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-70 components of targeted buffer : GF

    auto tr_x_0_x_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 60);

    auto tr_x_0_x_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 61);

    auto tr_x_0_x_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 62);

    auto tr_x_0_x_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 63);

    auto tr_x_0_x_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 64);

    auto tr_x_0_x_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 65);

    auto tr_x_0_x_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 66);

    auto tr_x_0_x_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 67);

    auto tr_x_0_x_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 68);

    auto tr_x_0_x_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 69);

    #pragma omp simd aligned(tr_x_0_x_xyyy_xxx, tr_x_0_x_xyyy_xxy, tr_x_0_x_xyyy_xxz, tr_x_0_x_xyyy_xyy, tr_x_0_x_xyyy_xyz, tr_x_0_x_xyyy_xzz, tr_x_0_x_xyyy_yyy, tr_x_0_x_xyyy_yyz, tr_x_0_x_xyyy_yzz, tr_x_0_x_xyyy_zzz, tr_xxxyyy_xxx, tr_xxxyyy_xxy, tr_xxxyyy_xxz, tr_xxxyyy_xyy, tr_xxxyyy_xyz, tr_xxxyyy_xzz, tr_xxxyyy_yyy, tr_xxxyyy_yyz, tr_xxxyyy_yzz, tr_xxxyyy_zzz, tr_xxyyy_xx, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xy, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xz, tr_xxyyy_xzzz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_yyy_xx, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyyy_xxx[i] = 3.0 * tr_yyy_xx[i] - 2.0 * tr_yyy_xxxx[i] * tke_0 - 6.0 * tr_xyyy_xxx[i] * tbe_0 - 6.0 * tr_xxyyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_xxy[i] = 2.0 * tr_yyy_xy[i] - 2.0 * tr_yyy_xxxy[i] * tke_0 - 6.0 * tr_xyyy_xxy[i] * tbe_0 - 4.0 * tr_xxyyy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_xxz[i] = 2.0 * tr_yyy_xz[i] - 2.0 * tr_yyy_xxxz[i] * tke_0 - 6.0 * tr_xyyy_xxz[i] * tbe_0 - 4.0 * tr_xxyyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_xyy[i] = tr_yyy_yy[i] - 2.0 * tr_yyy_xxyy[i] * tke_0 - 6.0 * tr_xyyy_xyy[i] * tbe_0 - 2.0 * tr_xxyyy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_xyz[i] = tr_yyy_yz[i] - 2.0 * tr_yyy_xxyz[i] * tke_0 - 6.0 * tr_xyyy_xyz[i] * tbe_0 - 2.0 * tr_xxyyy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_xzz[i] = tr_yyy_zz[i] - 2.0 * tr_yyy_xxzz[i] * tke_0 - 6.0 * tr_xyyy_xzz[i] * tbe_0 - 2.0 * tr_xxyyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_yyy[i] = -2.0 * tr_yyy_xyyy[i] * tke_0 - 6.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_yyz[i] = -2.0 * tr_yyy_xyyz[i] * tke_0 - 6.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_yzz[i] = -2.0 * tr_yyy_xyzz[i] * tke_0 - 6.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyy_zzz[i] = -2.0 * tr_yyy_xzzz[i] * tke_0 - 6.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 70-80 components of targeted buffer : GF

    auto tr_x_0_x_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 70);

    auto tr_x_0_x_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 71);

    auto tr_x_0_x_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 72);

    auto tr_x_0_x_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 73);

    auto tr_x_0_x_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 74);

    auto tr_x_0_x_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 75);

    auto tr_x_0_x_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 76);

    auto tr_x_0_x_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 77);

    auto tr_x_0_x_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 78);

    auto tr_x_0_x_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 79);

    #pragma omp simd aligned(tr_x_0_x_xyyz_xxx, tr_x_0_x_xyyz_xxy, tr_x_0_x_xyyz_xxz, tr_x_0_x_xyyz_xyy, tr_x_0_x_xyyz_xyz, tr_x_0_x_xyyz_xzz, tr_x_0_x_xyyz_yyy, tr_x_0_x_xyyz_yyz, tr_x_0_x_xyyz_yzz, tr_x_0_x_xyyz_zzz, tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxyyz_xx, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xy, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_xzzz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_yyz_xx, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyyz_xxx[i] = 3.0 * tr_yyz_xx[i] - 2.0 * tr_yyz_xxxx[i] * tke_0 - 6.0 * tr_xyyz_xxx[i] * tbe_0 - 6.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_xxy[i] = 2.0 * tr_yyz_xy[i] - 2.0 * tr_yyz_xxxy[i] * tke_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 - 4.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_xxz[i] = 2.0 * tr_yyz_xz[i] - 2.0 * tr_yyz_xxxz[i] * tke_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 - 4.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_xyy[i] = tr_yyz_yy[i] - 2.0 * tr_yyz_xxyy[i] * tke_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 - 2.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_xyz[i] = tr_yyz_yz[i] - 2.0 * tr_yyz_xxyz[i] * tke_0 - 6.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_xzz[i] = tr_yyz_zz[i] - 2.0 * tr_yyz_xxzz[i] * tke_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 - 2.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_yyy[i] = -2.0 * tr_yyz_xyyy[i] * tke_0 - 6.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_yyz[i] = -2.0 * tr_yyz_xyyz[i] * tke_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_yzz[i] = -2.0 * tr_yyz_xyzz[i] * tke_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyyz_zzz[i] = -2.0 * tr_yyz_xzzz[i] * tke_0 - 6.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 80-90 components of targeted buffer : GF

    auto tr_x_0_x_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 80);

    auto tr_x_0_x_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 81);

    auto tr_x_0_x_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 82);

    auto tr_x_0_x_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 83);

    auto tr_x_0_x_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 84);

    auto tr_x_0_x_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 85);

    auto tr_x_0_x_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 86);

    auto tr_x_0_x_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 87);

    auto tr_x_0_x_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 88);

    auto tr_x_0_x_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 89);

    #pragma omp simd aligned(tr_x_0_x_xyzz_xxx, tr_x_0_x_xyzz_xxy, tr_x_0_x_xyzz_xxz, tr_x_0_x_xyzz_xyy, tr_x_0_x_xyzz_xyz, tr_x_0_x_xyzz_xzz, tr_x_0_x_xyzz_yyy, tr_x_0_x_xyzz_yyz, tr_x_0_x_xyzz_yzz, tr_x_0_x_xyzz_zzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xy, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_xzzz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_yzz_xx, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyzz_xxx[i] = 3.0 * tr_yzz_xx[i] - 2.0 * tr_yzz_xxxx[i] * tke_0 - 6.0 * tr_xyzz_xxx[i] * tbe_0 - 6.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_xxy[i] = 2.0 * tr_yzz_xy[i] - 2.0 * tr_yzz_xxxy[i] * tke_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 - 4.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_xxz[i] = 2.0 * tr_yzz_xz[i] - 2.0 * tr_yzz_xxxz[i] * tke_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 - 4.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_xyy[i] = tr_yzz_yy[i] - 2.0 * tr_yzz_xxyy[i] * tke_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 - 2.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_xyz[i] = tr_yzz_yz[i] - 2.0 * tr_yzz_xxyz[i] * tke_0 - 6.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_xzz[i] = tr_yzz_zz[i] - 2.0 * tr_yzz_xxzz[i] * tke_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 - 2.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_yyy[i] = -2.0 * tr_yzz_xyyy[i] * tke_0 - 6.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_yyz[i] = -2.0 * tr_yzz_xyyz[i] * tke_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_yzz[i] = -2.0 * tr_yzz_xyzz[i] * tke_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyzz_zzz[i] = -2.0 * tr_yzz_xzzz[i] * tke_0 - 6.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-100 components of targeted buffer : GF

    auto tr_x_0_x_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 90);

    auto tr_x_0_x_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 91);

    auto tr_x_0_x_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 92);

    auto tr_x_0_x_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 93);

    auto tr_x_0_x_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 94);

    auto tr_x_0_x_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 95);

    auto tr_x_0_x_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 96);

    auto tr_x_0_x_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 97);

    auto tr_x_0_x_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 98);

    auto tr_x_0_x_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 99);

    #pragma omp simd aligned(tr_x_0_x_xzzz_xxx, tr_x_0_x_xzzz_xxy, tr_x_0_x_xzzz_xxz, tr_x_0_x_xzzz_xyy, tr_x_0_x_xzzz_xyz, tr_x_0_x_xzzz_xzz, tr_x_0_x_xzzz_yyy, tr_x_0_x_xzzz_yyz, tr_x_0_x_xzzz_yzz, tr_x_0_x_xzzz_zzz, tr_xxxzzz_xxx, tr_xxxzzz_xxy, tr_xxxzzz_xxz, tr_xxxzzz_xyy, tr_xxxzzz_xyz, tr_xxxzzz_xzz, tr_xxxzzz_yyy, tr_xxxzzz_yyz, tr_xxxzzz_yzz, tr_xxxzzz_zzz, tr_xxzzz_xx, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xy, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xz, tr_xxzzz_xzzz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_zzz_xx, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xzzz_xxx[i] = 3.0 * tr_zzz_xx[i] - 2.0 * tr_zzz_xxxx[i] * tke_0 - 6.0 * tr_xzzz_xxx[i] * tbe_0 - 6.0 * tr_xxzzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_xxy[i] = 2.0 * tr_zzz_xy[i] - 2.0 * tr_zzz_xxxy[i] * tke_0 - 6.0 * tr_xzzz_xxy[i] * tbe_0 - 4.0 * tr_xxzzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_xxz[i] = 2.0 * tr_zzz_xz[i] - 2.0 * tr_zzz_xxxz[i] * tke_0 - 6.0 * tr_xzzz_xxz[i] * tbe_0 - 4.0 * tr_xxzzz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_xyy[i] = tr_zzz_yy[i] - 2.0 * tr_zzz_xxyy[i] * tke_0 - 6.0 * tr_xzzz_xyy[i] * tbe_0 - 2.0 * tr_xxzzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_xyz[i] = tr_zzz_yz[i] - 2.0 * tr_zzz_xxyz[i] * tke_0 - 6.0 * tr_xzzz_xyz[i] * tbe_0 - 2.0 * tr_xxzzz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_xzz[i] = tr_zzz_zz[i] - 2.0 * tr_zzz_xxzz[i] * tke_0 - 6.0 * tr_xzzz_xzz[i] * tbe_0 - 2.0 * tr_xxzzz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_yyy[i] = -2.0 * tr_zzz_xyyy[i] * tke_0 - 6.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_yyz[i] = -2.0 * tr_zzz_xyyz[i] * tke_0 - 6.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_yzz[i] = -2.0 * tr_zzz_xyzz[i] * tke_0 - 6.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzzz_zzz[i] = -2.0 * tr_zzz_xzzz[i] * tke_0 - 6.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 100-110 components of targeted buffer : GF

    auto tr_x_0_x_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 100);

    auto tr_x_0_x_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 101);

    auto tr_x_0_x_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 102);

    auto tr_x_0_x_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 103);

    auto tr_x_0_x_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 104);

    auto tr_x_0_x_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 105);

    auto tr_x_0_x_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 106);

    auto tr_x_0_x_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 107);

    auto tr_x_0_x_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 108);

    auto tr_x_0_x_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 109);

    #pragma omp simd aligned(tr_x_0_x_yyyy_xxx, tr_x_0_x_yyyy_xxy, tr_x_0_x_yyyy_xxz, tr_x_0_x_yyyy_xyy, tr_x_0_x_yyyy_xyz, tr_x_0_x_yyyy_xzz, tr_x_0_x_yyyy_yyy, tr_x_0_x_yyyy_yyz, tr_x_0_x_yyyy_yzz, tr_x_0_x_yyyy_zzz, tr_xxyyyy_xxx, tr_xxyyyy_xxy, tr_xxyyyy_xxz, tr_xxyyyy_xyy, tr_xxyyyy_xyz, tr_xxyyyy_xzz, tr_xxyyyy_yyy, tr_xxyyyy_yyz, tr_xxyyyy_yzz, tr_xxyyyy_zzz, tr_xyyyy_xx, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xy, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xz, tr_xyyyy_xzzz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyyy_xxx[i] = -2.0 * tr_yyyy_xxx[i] * tbe_0 - 6.0 * tr_xyyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_xxy[i] = -2.0 * tr_yyyy_xxy[i] * tbe_0 - 4.0 * tr_xyyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_xxz[i] = -2.0 * tr_yyyy_xxz[i] * tbe_0 - 4.0 * tr_xyyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_xyy[i] = -2.0 * tr_yyyy_xyy[i] * tbe_0 - 2.0 * tr_xyyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_xyz[i] = -2.0 * tr_yyyy_xyz[i] * tbe_0 - 2.0 * tr_xyyyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_xzz[i] = -2.0 * tr_yyyy_xzz[i] * tbe_0 - 2.0 * tr_xyyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_yyy[i] = -2.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_yyz[i] = -2.0 * tr_yyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_yzz[i] = -2.0 * tr_yyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyy_zzz[i] = -2.0 * tr_yyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 110-120 components of targeted buffer : GF

    auto tr_x_0_x_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 110);

    auto tr_x_0_x_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 111);

    auto tr_x_0_x_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 112);

    auto tr_x_0_x_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 113);

    auto tr_x_0_x_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 114);

    auto tr_x_0_x_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 115);

    auto tr_x_0_x_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 116);

    auto tr_x_0_x_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 117);

    auto tr_x_0_x_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 118);

    auto tr_x_0_x_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 119);

    #pragma omp simd aligned(tr_x_0_x_yyyz_xxx, tr_x_0_x_yyyz_xxy, tr_x_0_x_yyyz_xxz, tr_x_0_x_yyyz_xyy, tr_x_0_x_yyyz_xyz, tr_x_0_x_yyyz_xzz, tr_x_0_x_yyyz_yyy, tr_x_0_x_yyyz_yyz, tr_x_0_x_yyyz_yzz, tr_x_0_x_yyyz_zzz, tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xyyyz_xx, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xy, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_xzzz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyyz_xxx[i] = -2.0 * tr_yyyz_xxx[i] * tbe_0 - 6.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_xxy[i] = -2.0 * tr_yyyz_xxy[i] * tbe_0 - 4.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_xxz[i] = -2.0 * tr_yyyz_xxz[i] * tbe_0 - 4.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_xyy[i] = -2.0 * tr_yyyz_xyy[i] * tbe_0 - 2.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_xyz[i] = -2.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_xzz[i] = -2.0 * tr_yyyz_xzz[i] * tbe_0 - 2.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_yyy[i] = -2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_yyz[i] = -2.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_yzz[i] = -2.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyyz_zzz[i] = -2.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-130 components of targeted buffer : GF

    auto tr_x_0_x_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 120);

    auto tr_x_0_x_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 121);

    auto tr_x_0_x_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 122);

    auto tr_x_0_x_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 123);

    auto tr_x_0_x_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 124);

    auto tr_x_0_x_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 125);

    auto tr_x_0_x_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 126);

    auto tr_x_0_x_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 127);

    auto tr_x_0_x_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 128);

    auto tr_x_0_x_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 129);

    #pragma omp simd aligned(tr_x_0_x_yyzz_xxx, tr_x_0_x_yyzz_xxy, tr_x_0_x_yyzz_xxz, tr_x_0_x_yyzz_xyy, tr_x_0_x_yyzz_xyz, tr_x_0_x_yyzz_xzz, tr_x_0_x_yyzz_yyy, tr_x_0_x_yyzz_yyz, tr_x_0_x_yyzz_yzz, tr_x_0_x_yyzz_zzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xy, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_xzzz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyzz_xxx[i] = -2.0 * tr_yyzz_xxx[i] * tbe_0 - 6.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_xxy[i] = -2.0 * tr_yyzz_xxy[i] * tbe_0 - 4.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_xxz[i] = -2.0 * tr_yyzz_xxz[i] * tbe_0 - 4.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_xyy[i] = -2.0 * tr_yyzz_xyy[i] * tbe_0 - 2.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_xyz[i] = -2.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_xzz[i] = -2.0 * tr_yyzz_xzz[i] * tbe_0 - 2.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_yyy[i] = -2.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_yyz[i] = -2.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_yzz[i] = -2.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyzz_zzz[i] = -2.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 130-140 components of targeted buffer : GF

    auto tr_x_0_x_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 130);

    auto tr_x_0_x_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 131);

    auto tr_x_0_x_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 132);

    auto tr_x_0_x_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 133);

    auto tr_x_0_x_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 134);

    auto tr_x_0_x_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 135);

    auto tr_x_0_x_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 136);

    auto tr_x_0_x_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 137);

    auto tr_x_0_x_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 138);

    auto tr_x_0_x_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 139);

    #pragma omp simd aligned(tr_x_0_x_yzzz_xxx, tr_x_0_x_yzzz_xxy, tr_x_0_x_yzzz_xxz, tr_x_0_x_yzzz_xyy, tr_x_0_x_yzzz_xyz, tr_x_0_x_yzzz_xzz, tr_x_0_x_yzzz_yyy, tr_x_0_x_yzzz_yyz, tr_x_0_x_yzzz_yzz, tr_x_0_x_yzzz_zzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xy, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_xzzz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yzzz_xxx[i] = -2.0 * tr_yzzz_xxx[i] * tbe_0 - 6.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_xxy[i] = -2.0 * tr_yzzz_xxy[i] * tbe_0 - 4.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_xxz[i] = -2.0 * tr_yzzz_xxz[i] * tbe_0 - 4.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_xyy[i] = -2.0 * tr_yzzz_xyy[i] * tbe_0 - 2.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_xyz[i] = -2.0 * tr_yzzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_xzz[i] = -2.0 * tr_yzzz_xzz[i] * tbe_0 - 2.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_yyy[i] = -2.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_yyz[i] = -2.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_yzz[i] = -2.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzzz_zzz[i] = -2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 140-150 components of targeted buffer : GF

    auto tr_x_0_x_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 140);

    auto tr_x_0_x_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 141);

    auto tr_x_0_x_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 142);

    auto tr_x_0_x_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 143);

    auto tr_x_0_x_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 144);

    auto tr_x_0_x_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 145);

    auto tr_x_0_x_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 146);

    auto tr_x_0_x_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 147);

    auto tr_x_0_x_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 148);

    auto tr_x_0_x_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 149);

    #pragma omp simd aligned(tr_x_0_x_zzzz_xxx, tr_x_0_x_zzzz_xxy, tr_x_0_x_zzzz_xxz, tr_x_0_x_zzzz_xyy, tr_x_0_x_zzzz_xyz, tr_x_0_x_zzzz_xzz, tr_x_0_x_zzzz_yyy, tr_x_0_x_zzzz_yyz, tr_x_0_x_zzzz_yzz, tr_x_0_x_zzzz_zzz, tr_xxzzzz_xxx, tr_xxzzzz_xxy, tr_xxzzzz_xxz, tr_xxzzzz_xyy, tr_xxzzzz_xyz, tr_xxzzzz_xzz, tr_xxzzzz_yyy, tr_xxzzzz_yyz, tr_xxzzzz_yzz, tr_xxzzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xy, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xz, tr_xzzzz_xzzz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zzzz_xxx[i] = -2.0 * tr_zzzz_xxx[i] * tbe_0 - 6.0 * tr_xzzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_xxy[i] = -2.0 * tr_zzzz_xxy[i] * tbe_0 - 4.0 * tr_xzzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_xxz[i] = -2.0 * tr_zzzz_xxz[i] * tbe_0 - 4.0 * tr_xzzzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_xyy[i] = -2.0 * tr_zzzz_xyy[i] * tbe_0 - 2.0 * tr_xzzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_xyz[i] = -2.0 * tr_zzzz_xyz[i] * tbe_0 - 2.0 * tr_xzzzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_xzz[i] = -2.0 * tr_zzzz_xzz[i] * tbe_0 - 2.0 * tr_xzzzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_yyy[i] = -2.0 * tr_zzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_yyz[i] = -2.0 * tr_zzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_yzz[i] = -2.0 * tr_zzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzzz_zzz[i] = -2.0 * tr_zzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-160 components of targeted buffer : GF

    auto tr_x_0_y_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 150);

    auto tr_x_0_y_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 151);

    auto tr_x_0_y_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 152);

    auto tr_x_0_y_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 153);

    auto tr_x_0_y_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 154);

    auto tr_x_0_y_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 155);

    auto tr_x_0_y_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 156);

    auto tr_x_0_y_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 157);

    auto tr_x_0_y_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 158);

    auto tr_x_0_y_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 159);

    #pragma omp simd aligned(tr_x_0_y_xxxx_xxx, tr_x_0_y_xxxx_xxy, tr_x_0_y_xxxx_xxz, tr_x_0_y_xxxx_xyy, tr_x_0_y_xxxx_xyz, tr_x_0_y_xxxx_xzz, tr_x_0_y_xxxx_yyy, tr_x_0_y_xxxx_yyz, tr_x_0_y_xxxx_yzz, tr_x_0_y_xxxx_zzz, tr_xxx_xx, tr_xxx_xxxy, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_yy, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxxxx_xx, tr_xxxxx_xxxy, tr_xxxxx_xxyy, tr_xxxxx_xxyz, tr_xxxxx_xy, tr_xxxxx_xyyy, tr_xxxxx_xyyz, tr_xxxxx_xyzz, tr_xxxxx_xz, tr_xxxxx_yy, tr_xxxxx_yyyy, tr_xxxxx_yyyz, tr_xxxxx_yyzz, tr_xxxxx_yz, tr_xxxxx_yzzz, tr_xxxxx_zz, tr_xxxxxy_xxx, tr_xxxxxy_xxy, tr_xxxxxy_xxz, tr_xxxxxy_xyy, tr_xxxxxy_xyz, tr_xxxxxy_xzz, tr_xxxxxy_yyy, tr_xxxxxy_yyz, tr_xxxxxy_yzz, tr_xxxxxy_zzz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxxx_xxx[i] = -8.0 * tr_xxx_xxxy[i] * tke_0 - 8.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxxx_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_xxy[i] = 4.0 * tr_xxx_xx[i] - 8.0 * tr_xxx_xxyy[i] * tke_0 - 8.0 * tr_xxxy_xxy[i] * tbe_0 - 2.0 * tr_xxxxx_xx[i] * tbe_0 + 4.0 * tr_xxxxx_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_xxz[i] = -8.0 * tr_xxx_xxyz[i] * tke_0 - 8.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_xyy[i] = 8.0 * tr_xxx_xy[i] - 8.0 * tr_xxx_xyyy[i] * tke_0 - 8.0 * tr_xxxy_xyy[i] * tbe_0 - 4.0 * tr_xxxxx_xy[i] * tbe_0 + 4.0 * tr_xxxxx_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_xyz[i] = 4.0 * tr_xxx_xz[i] - 8.0 * tr_xxx_xyyz[i] * tke_0 - 8.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxxx_xz[i] * tbe_0 + 4.0 * tr_xxxxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_xzz[i] = -8.0 * tr_xxx_xyzz[i] * tke_0 - 8.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_yyy[i] = 12.0 * tr_xxx_yy[i] - 8.0 * tr_xxx_yyyy[i] * tke_0 - 8.0 * tr_xxxy_yyy[i] * tbe_0 - 6.0 * tr_xxxxx_yy[i] * tbe_0 + 4.0 * tr_xxxxx_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_yyz[i] = 8.0 * tr_xxx_yz[i] - 8.0 * tr_xxx_yyyz[i] * tke_0 - 8.0 * tr_xxxy_yyz[i] * tbe_0 - 4.0 * tr_xxxxx_yz[i] * tbe_0 + 4.0 * tr_xxxxx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_yzz[i] = 4.0 * tr_xxx_zz[i] - 8.0 * tr_xxx_yyzz[i] * tke_0 - 8.0 * tr_xxxy_yzz[i] * tbe_0 - 2.0 * tr_xxxxx_zz[i] * tbe_0 + 4.0 * tr_xxxxx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxx_zzz[i] = -8.0 * tr_xxx_yzzz[i] * tke_0 - 8.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 160-170 components of targeted buffer : GF

    auto tr_x_0_y_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 160);

    auto tr_x_0_y_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 161);

    auto tr_x_0_y_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 162);

    auto tr_x_0_y_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 163);

    auto tr_x_0_y_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 164);

    auto tr_x_0_y_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 165);

    auto tr_x_0_y_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 166);

    auto tr_x_0_y_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 167);

    auto tr_x_0_y_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 168);

    auto tr_x_0_y_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 169);

    #pragma omp simd aligned(tr_x_0_y_xxxy_xxx, tr_x_0_y_xxxy_xxy, tr_x_0_y_xxxy_xxz, tr_x_0_y_xxxy_xyy, tr_x_0_y_xxxy_xyz, tr_x_0_y_xxxy_xzz, tr_x_0_y_xxxy_yyy, tr_x_0_y_xxxy_yyz, tr_x_0_y_xxxy_yzz, tr_x_0_y_xxxy_zzz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxy_xx, tr_xxxxy_xxxy, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xy, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yz, tr_xxxxy_yzzz, tr_xxxxy_zz, tr_xxxxyy_xxx, tr_xxxxyy_xxy, tr_xxxxyy_xxz, tr_xxxxyy_xyy, tr_xxxxyy_xyz, tr_xxxxyy_xzz, tr_xxxxyy_yyy, tr_xxxxyy_yyz, tr_xxxxyy_yzz, tr_xxxxyy_zzz, tr_xxy_xx, tr_xxy_xxxy, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_yy, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxxy_xxx[i] = 3.0 * tr_xx_xxx[i] - 6.0 * tr_xxy_xxxy[i] * tke_0 - 6.0 * tr_xxyy_xxx[i] * tbe_0 - 2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_xxy[i] = 3.0 * tr_xx_xxy[i] + 3.0 * tr_xxy_xx[i] - 6.0 * tr_xxy_xxyy[i] * tke_0 - 6.0 * tr_xxyy_xxy[i] * tbe_0 - 2.0 * tr_xxxx_xxy[i] * tbe_0 - 2.0 * tr_xxxxy_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_xxz[i] = 3.0 * tr_xx_xxz[i] - 6.0 * tr_xxy_xxyz[i] * tke_0 - 6.0 * tr_xxyy_xxz[i] * tbe_0 - 2.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_xyy[i] = 3.0 * tr_xx_xyy[i] + 6.0 * tr_xxy_xy[i] - 6.0 * tr_xxy_xyyy[i] * tke_0 - 6.0 * tr_xxyy_xyy[i] * tbe_0 - 2.0 * tr_xxxx_xyy[i] * tbe_0 - 4.0 * tr_xxxxy_xy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_xyz[i] = 3.0 * tr_xx_xyz[i] + 3.0 * tr_xxy_xz[i] - 6.0 * tr_xxy_xyyz[i] * tke_0 - 6.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxxx_xyz[i] * tbe_0 - 2.0 * tr_xxxxy_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_xzz[i] = 3.0 * tr_xx_xzz[i] - 6.0 * tr_xxy_xyzz[i] * tke_0 - 6.0 * tr_xxyy_xzz[i] * tbe_0 - 2.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_yyy[i] = 3.0 * tr_xx_yyy[i] + 9.0 * tr_xxy_yy[i] - 6.0 * tr_xxy_yyyy[i] * tke_0 - 6.0 * tr_xxyy_yyy[i] * tbe_0 - 2.0 * tr_xxxx_yyy[i] * tbe_0 - 6.0 * tr_xxxxy_yy[i] * tbe_0 + 4.0 * tr_xxxxy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_yyz[i] = 3.0 * tr_xx_yyz[i] + 6.0 * tr_xxy_yz[i] - 6.0 * tr_xxy_yyyz[i] * tke_0 - 6.0 * tr_xxyy_yyz[i] * tbe_0 - 2.0 * tr_xxxx_yyz[i] * tbe_0 - 4.0 * tr_xxxxy_yz[i] * tbe_0 + 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_yzz[i] = 3.0 * tr_xx_yzz[i] + 3.0 * tr_xxy_zz[i] - 6.0 * tr_xxy_yyzz[i] * tke_0 - 6.0 * tr_xxyy_yzz[i] * tbe_0 - 2.0 * tr_xxxx_yzz[i] * tbe_0 - 2.0 * tr_xxxxy_zz[i] * tbe_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxy_zzz[i] = 3.0 * tr_xx_zzz[i] - 6.0 * tr_xxy_yzzz[i] * tke_0 - 6.0 * tr_xxyy_zzz[i] * tbe_0 - 2.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 170-180 components of targeted buffer : GF

    auto tr_x_0_y_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 170);

    auto tr_x_0_y_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 171);

    auto tr_x_0_y_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 172);

    auto tr_x_0_y_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 173);

    auto tr_x_0_y_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 174);

    auto tr_x_0_y_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 175);

    auto tr_x_0_y_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 176);

    auto tr_x_0_y_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 177);

    auto tr_x_0_y_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 178);

    auto tr_x_0_y_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 179);

    #pragma omp simd aligned(tr_x_0_y_xxxz_xxx, tr_x_0_y_xxxz_xxy, tr_x_0_y_xxxz_xxz, tr_x_0_y_xxxz_xyy, tr_x_0_y_xxxz_xyz, tr_x_0_y_xxxz_xzz, tr_x_0_y_xxxz_yyy, tr_x_0_y_xxxz_yyz, tr_x_0_y_xxxz_yzz, tr_x_0_y_xxxz_zzz, tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_xxxxz_xx, tr_xxxxz_xxxy, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xy, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yz, tr_xxxxz_yzzz, tr_xxxxz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxz_xx, tr_xxz_xxxy, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_yy, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxxz_xxx[i] = -6.0 * tr_xxz_xxxy[i] * tke_0 - 6.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_xxy[i] = 3.0 * tr_xxz_xx[i] - 6.0 * tr_xxz_xxyy[i] * tke_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 - 2.0 * tr_xxxxz_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_xxz[i] = -6.0 * tr_xxz_xxyz[i] * tke_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_xyy[i] = 6.0 * tr_xxz_xy[i] - 6.0 * tr_xxz_xyyy[i] * tke_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 - 4.0 * tr_xxxxz_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_xyz[i] = 3.0 * tr_xxz_xz[i] - 6.0 * tr_xxz_xyyz[i] * tke_0 - 6.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxxxz_xz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_xzz[i] = -6.0 * tr_xxz_xyzz[i] * tke_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_yyy[i] = 9.0 * tr_xxz_yy[i] - 6.0 * tr_xxz_yyyy[i] * tke_0 - 6.0 * tr_xxyz_yyy[i] * tbe_0 - 6.0 * tr_xxxxz_yy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_yyz[i] = 6.0 * tr_xxz_yz[i] - 6.0 * tr_xxz_yyyz[i] * tke_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 - 4.0 * tr_xxxxz_yz[i] * tbe_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_yzz[i] = 3.0 * tr_xxz_zz[i] - 6.0 * tr_xxz_yyzz[i] * tke_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 - 2.0 * tr_xxxxz_zz[i] * tbe_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxxz_zzz[i] = -6.0 * tr_xxz_yzzz[i] * tke_0 - 6.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-190 components of targeted buffer : GF

    auto tr_x_0_y_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 180);

    auto tr_x_0_y_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 181);

    auto tr_x_0_y_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 182);

    auto tr_x_0_y_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 183);

    auto tr_x_0_y_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 184);

    auto tr_x_0_y_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 185);

    auto tr_x_0_y_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 186);

    auto tr_x_0_y_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 187);

    auto tr_x_0_y_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 188);

    auto tr_x_0_y_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 189);

    #pragma omp simd aligned(tr_x_0_y_xxyy_xxx, tr_x_0_y_xxyy_xxy, tr_x_0_y_xxyy_xxz, tr_x_0_y_xxyy_xyy, tr_x_0_y_xxyy_xyz, tr_x_0_y_xxyy_xzz, tr_x_0_y_xxyy_yyy, tr_x_0_y_xxyy_yyz, tr_x_0_y_xxyy_yzz, tr_x_0_y_xxyy_zzz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyy_xx, tr_xxxyy_xxxy, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xy, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yz, tr_xxxyy_yzzz, tr_xxxyy_zz, tr_xxxyyy_xxx, tr_xxxyyy_xxy, tr_xxxyyy_xxz, tr_xxxyyy_xyy, tr_xxxyyy_xyz, tr_xxxyyy_xzz, tr_xxxyyy_yyy, tr_xxxyyy_yyz, tr_xxxyyy_yzz, tr_xxxyyy_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxy, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_yy, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxyy_xxx[i] = 4.0 * tr_xy_xxx[i] - 4.0 * tr_xyy_xxxy[i] * tke_0 - 4.0 * tr_xyyy_xxx[i] * tbe_0 - 4.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_xxy[i] = 4.0 * tr_xy_xxy[i] + 2.0 * tr_xyy_xx[i] - 4.0 * tr_xyy_xxyy[i] * tke_0 - 4.0 * tr_xyyy_xxy[i] * tbe_0 - 4.0 * tr_xxxy_xxy[i] * tbe_0 - 2.0 * tr_xxxyy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_xxz[i] = 4.0 * tr_xy_xxz[i] - 4.0 * tr_xyy_xxyz[i] * tke_0 - 4.0 * tr_xyyy_xxz[i] * tbe_0 - 4.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_xyy[i] = 4.0 * tr_xy_xyy[i] + 4.0 * tr_xyy_xy[i] - 4.0 * tr_xyy_xyyy[i] * tke_0 - 4.0 * tr_xyyy_xyy[i] * tbe_0 - 4.0 * tr_xxxy_xyy[i] * tbe_0 - 4.0 * tr_xxxyy_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_xyz[i] = 4.0 * tr_xy_xyz[i] + 2.0 * tr_xyy_xz[i] - 4.0 * tr_xyy_xyyz[i] * tke_0 - 4.0 * tr_xyyy_xyz[i] * tbe_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxyy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_xzz[i] = 4.0 * tr_xy_xzz[i] - 4.0 * tr_xyy_xyzz[i] * tke_0 - 4.0 * tr_xyyy_xzz[i] * tbe_0 - 4.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_yyy[i] = 4.0 * tr_xy_yyy[i] + 6.0 * tr_xyy_yy[i] - 4.0 * tr_xyy_yyyy[i] * tke_0 - 4.0 * tr_xyyy_yyy[i] * tbe_0 - 4.0 * tr_xxxy_yyy[i] * tbe_0 - 6.0 * tr_xxxyy_yy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_yyz[i] = 4.0 * tr_xy_yyz[i] + 4.0 * tr_xyy_yz[i] - 4.0 * tr_xyy_yyyz[i] * tke_0 - 4.0 * tr_xyyy_yyz[i] * tbe_0 - 4.0 * tr_xxxy_yyz[i] * tbe_0 - 4.0 * tr_xxxyy_yz[i] * tbe_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_yzz[i] = 4.0 * tr_xy_yzz[i] + 2.0 * tr_xyy_zz[i] - 4.0 * tr_xyy_yyzz[i] * tke_0 - 4.0 * tr_xyyy_yzz[i] * tbe_0 - 4.0 * tr_xxxy_yzz[i] * tbe_0 - 2.0 * tr_xxxyy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyy_zzz[i] = 4.0 * tr_xy_zzz[i] - 4.0 * tr_xyy_yzzz[i] * tke_0 - 4.0 * tr_xyyy_zzz[i] * tbe_0 - 4.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 190-200 components of targeted buffer : GF

    auto tr_x_0_y_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 190);

    auto tr_x_0_y_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 191);

    auto tr_x_0_y_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 192);

    auto tr_x_0_y_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 193);

    auto tr_x_0_y_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 194);

    auto tr_x_0_y_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 195);

    auto tr_x_0_y_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 196);

    auto tr_x_0_y_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 197);

    auto tr_x_0_y_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 198);

    auto tr_x_0_y_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 199);

    #pragma omp simd aligned(tr_x_0_y_xxyz_xxx, tr_x_0_y_xxyz_xxy, tr_x_0_y_xxyz_xxz, tr_x_0_y_xxyz_xyy, tr_x_0_y_xxyz_xyz, tr_x_0_y_xxyz_xzz, tr_x_0_y_xxyz_yyy, tr_x_0_y_xxyz_yyz, tr_x_0_y_xxyz_yzz, tr_x_0_y_xxyz_zzz, tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxxyz_xx, tr_xxxyz_xxxy, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xy, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yz, tr_xxxyz_yzzz, tr_xxxyz_zz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxyz_xxx[i] = 2.0 * tr_xz_xxx[i] - 4.0 * tr_xyz_xxxy[i] * tke_0 - 4.0 * tr_xyyz_xxx[i] * tbe_0 - 2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_xxy[i] = 2.0 * tr_xz_xxy[i] + 2.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxyy[i] * tke_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 - 2.0 * tr_xxxz_xxy[i] * tbe_0 - 2.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_xxz[i] = 2.0 * tr_xz_xxz[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 - 2.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_xyy[i] = 2.0 * tr_xz_xyy[i] + 4.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xyyy[i] * tke_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 - 2.0 * tr_xxxz_xyy[i] * tbe_0 - 4.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_xyz[i] = 2.0 * tr_xz_xyz[i] + 2.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_xzz[i] = 2.0 * tr_xz_xzz[i] - 4.0 * tr_xyz_xyzz[i] * tke_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 - 2.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_yyy[i] = 2.0 * tr_xz_yyy[i] + 6.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_yyyy[i] * tke_0 - 4.0 * tr_xyyz_yyy[i] * tbe_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 - 6.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_yyz[i] = 2.0 * tr_xz_yyz[i] + 4.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_yyyz[i] * tke_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xxxz_yyz[i] * tbe_0 - 4.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_yzz[i] = 2.0 * tr_xz_yzz[i] + 2.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_yyzz[i] * tke_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 - 2.0 * tr_xxxz_yzz[i] * tbe_0 - 2.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxyz_zzz[i] = 2.0 * tr_xz_zzz[i] - 4.0 * tr_xyz_yzzz[i] * tke_0 - 4.0 * tr_xyyz_zzz[i] * tbe_0 - 2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 200-210 components of targeted buffer : GF

    auto tr_x_0_y_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 200);

    auto tr_x_0_y_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 201);

    auto tr_x_0_y_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 202);

    auto tr_x_0_y_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 203);

    auto tr_x_0_y_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 204);

    auto tr_x_0_y_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 205);

    auto tr_x_0_y_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 206);

    auto tr_x_0_y_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 207);

    auto tr_x_0_y_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 208);

    auto tr_x_0_y_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 209);

    #pragma omp simd aligned(tr_x_0_y_xxzz_xxx, tr_x_0_y_xxzz_xxy, tr_x_0_y_xxzz_xxz, tr_x_0_y_xxzz_xyy, tr_x_0_y_xxzz_xyz, tr_x_0_y_xxzz_xzz, tr_x_0_y_xxzz_yyy, tr_x_0_y_xxzz_yyz, tr_x_0_y_xxzz_yzz, tr_x_0_y_xxzz_zzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_xxxzz_xx, tr_xxxzz_xxxy, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xy, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yz, tr_xxxzz_yzzz, tr_xxxzz_zz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xzz_xx, tr_xzz_xxxy, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_yy, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxzz_xxx[i] = -4.0 * tr_xzz_xxxy[i] * tke_0 - 4.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_xxy[i] = 2.0 * tr_xzz_xx[i] - 4.0 * tr_xzz_xxyy[i] * tke_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 - 2.0 * tr_xxxzz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_xxz[i] = -4.0 * tr_xzz_xxyz[i] * tke_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_xyy[i] = 4.0 * tr_xzz_xy[i] - 4.0 * tr_xzz_xyyy[i] * tke_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 - 4.0 * tr_xxxzz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_xyz[i] = 2.0 * tr_xzz_xz[i] - 4.0 * tr_xzz_xyyz[i] * tke_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xxxzz_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_xzz[i] = -4.0 * tr_xzz_xyzz[i] * tke_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_yyy[i] = 6.0 * tr_xzz_yy[i] - 4.0 * tr_xzz_yyyy[i] * tke_0 - 4.0 * tr_xyzz_yyy[i] * tbe_0 - 6.0 * tr_xxxzz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_yyz[i] = 4.0 * tr_xzz_yz[i] - 4.0 * tr_xzz_yyyz[i] * tke_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 - 4.0 * tr_xxxzz_yz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_yzz[i] = 2.0 * tr_xzz_zz[i] - 4.0 * tr_xzz_yyzz[i] * tke_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 - 2.0 * tr_xxxzz_zz[i] * tbe_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxzz_zzz[i] = -4.0 * tr_xzz_yzzz[i] * tke_0 - 4.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-220 components of targeted buffer : GF

    auto tr_x_0_y_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 210);

    auto tr_x_0_y_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 211);

    auto tr_x_0_y_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 212);

    auto tr_x_0_y_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 213);

    auto tr_x_0_y_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 214);

    auto tr_x_0_y_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 215);

    auto tr_x_0_y_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 216);

    auto tr_x_0_y_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 217);

    auto tr_x_0_y_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 218);

    auto tr_x_0_y_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 219);

    #pragma omp simd aligned(tr_x_0_y_xyyy_xxx, tr_x_0_y_xyyy_xxy, tr_x_0_y_xyyy_xxz, tr_x_0_y_xyyy_xyy, tr_x_0_y_xyyy_xyz, tr_x_0_y_xyyy_xzz, tr_x_0_y_xyyy_yyy, tr_x_0_y_xyyy_yyz, tr_x_0_y_xyyy_yzz, tr_x_0_y_xyyy_zzz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyy_xx, tr_xxyyy_xxxy, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xy, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yz, tr_xxyyy_yzzz, tr_xxyyy_zz, tr_xxyyyy_xxx, tr_xxyyyy_xxy, tr_xxyyyy_xxz, tr_xxyyyy_xyy, tr_xxyyyy_xyz, tr_xxyyyy_xzz, tr_xxyyyy_yyy, tr_xxyyyy_yyz, tr_xxyyyy_yzz, tr_xxyyyy_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyy_xx, tr_yyy_xxxy, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_yy, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyyy_xxx[i] = 3.0 * tr_yy_xxx[i] - 2.0 * tr_yyy_xxxy[i] * tke_0 - 2.0 * tr_yyyy_xxx[i] * tbe_0 - 6.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_xxy[i] = 3.0 * tr_yy_xxy[i] + tr_yyy_xx[i] - 2.0 * tr_yyy_xxyy[i] * tke_0 - 2.0 * tr_yyyy_xxy[i] * tbe_0 - 6.0 * tr_xxyy_xxy[i] * tbe_0 - 2.0 * tr_xxyyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_xxz[i] = 3.0 * tr_yy_xxz[i] - 2.0 * tr_yyy_xxyz[i] * tke_0 - 2.0 * tr_yyyy_xxz[i] * tbe_0 - 6.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_xyy[i] = 3.0 * tr_yy_xyy[i] + 2.0 * tr_yyy_xy[i] - 2.0 * tr_yyy_xyyy[i] * tke_0 - 2.0 * tr_yyyy_xyy[i] * tbe_0 - 6.0 * tr_xxyy_xyy[i] * tbe_0 - 4.0 * tr_xxyyy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_xyz[i] = 3.0 * tr_yy_xyz[i] + tr_yyy_xz[i] - 2.0 * tr_yyy_xyyz[i] * tke_0 - 2.0 * tr_yyyy_xyz[i] * tbe_0 - 6.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxyyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_xzz[i] = 3.0 * tr_yy_xzz[i] - 2.0 * tr_yyy_xyzz[i] * tke_0 - 2.0 * tr_yyyy_xzz[i] * tbe_0 - 6.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_yyy[i] = 3.0 * tr_yy_yyy[i] + 3.0 * tr_yyy_yy[i] - 2.0 * tr_yyy_yyyy[i] * tke_0 - 2.0 * tr_yyyy_yyy[i] * tbe_0 - 6.0 * tr_xxyy_yyy[i] * tbe_0 - 6.0 * tr_xxyyy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_yyz[i] = 3.0 * tr_yy_yyz[i] + 2.0 * tr_yyy_yz[i] - 2.0 * tr_yyy_yyyz[i] * tke_0 - 2.0 * tr_yyyy_yyz[i] * tbe_0 - 6.0 * tr_xxyy_yyz[i] * tbe_0 - 4.0 * tr_xxyyy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_yzz[i] = 3.0 * tr_yy_yzz[i] + tr_yyy_zz[i] - 2.0 * tr_yyy_yyzz[i] * tke_0 - 2.0 * tr_yyyy_yzz[i] * tbe_0 - 6.0 * tr_xxyy_yzz[i] * tbe_0 - 2.0 * tr_xxyyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyy_zzz[i] = 3.0 * tr_yy_zzz[i] - 2.0 * tr_yyy_yzzz[i] * tke_0 - 2.0 * tr_yyyy_zzz[i] * tbe_0 - 6.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 220-230 components of targeted buffer : GF

    auto tr_x_0_y_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 220);

    auto tr_x_0_y_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 221);

    auto tr_x_0_y_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 222);

    auto tr_x_0_y_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 223);

    auto tr_x_0_y_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 224);

    auto tr_x_0_y_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 225);

    auto tr_x_0_y_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 226);

    auto tr_x_0_y_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 227);

    auto tr_x_0_y_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 228);

    auto tr_x_0_y_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 229);

    #pragma omp simd aligned(tr_x_0_y_xyyz_xxx, tr_x_0_y_xyyz_xxy, tr_x_0_y_xyyz_xxz, tr_x_0_y_xyyz_xyy, tr_x_0_y_xyyz_xyz, tr_x_0_y_xyyz_xzz, tr_x_0_y_xyyz_yyy, tr_x_0_y_xyyz_yyz, tr_x_0_y_xyyz_yzz, tr_x_0_y_xyyz_zzz, tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xxyyz_xx, tr_xxyyz_xxxy, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xy, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yz, tr_xxyyz_yzzz, tr_xxyyz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyz_xx, tr_yyz_xxxy, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_yy, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyyz_xxx[i] = 2.0 * tr_yz_xxx[i] - 2.0 * tr_yyz_xxxy[i] * tke_0 - 2.0 * tr_yyyz_xxx[i] * tbe_0 - 4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_xxy[i] = 2.0 * tr_yz_xxy[i] + tr_yyz_xx[i] - 2.0 * tr_yyz_xxyy[i] * tke_0 - 2.0 * tr_yyyz_xxy[i] * tbe_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 - 2.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_xxz[i] = 2.0 * tr_yz_xxz[i] - 2.0 * tr_yyz_xxyz[i] * tke_0 - 2.0 * tr_yyyz_xxz[i] * tbe_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_xyy[i] = 2.0 * tr_yz_xyy[i] + 2.0 * tr_yyz_xy[i] - 2.0 * tr_yyz_xyyy[i] * tke_0 - 2.0 * tr_yyyz_xyy[i] * tbe_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 - 4.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_xyz[i] = 2.0 * tr_yz_xyz[i] + tr_yyz_xz[i] - 2.0 * tr_yyz_xyyz[i] * tke_0 - 2.0 * tr_yyyz_xyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_xzz[i] = 2.0 * tr_yz_xzz[i] - 2.0 * tr_yyz_xyzz[i] * tke_0 - 2.0 * tr_yyyz_xzz[i] * tbe_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_yyy[i] = 2.0 * tr_yz_yyy[i] + 3.0 * tr_yyz_yy[i] - 2.0 * tr_yyz_yyyy[i] * tke_0 - 2.0 * tr_yyyz_yyy[i] * tbe_0 - 4.0 * tr_xxyz_yyy[i] * tbe_0 - 6.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_yyz[i] = 2.0 * tr_yz_yyz[i] + 2.0 * tr_yyz_yz[i] - 2.0 * tr_yyz_yyyz[i] * tke_0 - 2.0 * tr_yyyz_yyz[i] * tbe_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 - 4.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_yzz[i] = 2.0 * tr_yz_yzz[i] + tr_yyz_zz[i] - 2.0 * tr_yyz_yyzz[i] * tke_0 - 2.0 * tr_yyyz_yzz[i] * tbe_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 - 2.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyyz_zzz[i] = 2.0 * tr_yz_zzz[i] - 2.0 * tr_yyz_yzzz[i] * tke_0 - 2.0 * tr_yyyz_zzz[i] * tbe_0 - 4.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 230-240 components of targeted buffer : GF

    auto tr_x_0_y_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 230);

    auto tr_x_0_y_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 231);

    auto tr_x_0_y_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 232);

    auto tr_x_0_y_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 233);

    auto tr_x_0_y_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 234);

    auto tr_x_0_y_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 235);

    auto tr_x_0_y_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 236);

    auto tr_x_0_y_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 237);

    auto tr_x_0_y_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 238);

    auto tr_x_0_y_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 239);

    #pragma omp simd aligned(tr_x_0_y_xyzz_xxx, tr_x_0_y_xyzz_xxy, tr_x_0_y_xyzz_xxz, tr_x_0_y_xyzz_xyy, tr_x_0_y_xyzz_xyz, tr_x_0_y_xyzz_xzz, tr_x_0_y_xyzz_yyy, tr_x_0_y_xyzz_yyz, tr_x_0_y_xyzz_yzz, tr_x_0_y_xyzz_zzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxy, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xy, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yz, tr_xxyzz_yzzz, tr_xxyzz_zz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yzz_xx, tr_yzz_xxxy, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_yy, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyzz_xxx[i] = tr_zz_xxx[i] - 2.0 * tr_yzz_xxxy[i] * tke_0 - 2.0 * tr_yyzz_xxx[i] * tbe_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_xxy[i] = tr_zz_xxy[i] + tr_yzz_xx[i] - 2.0 * tr_yzz_xxyy[i] * tke_0 - 2.0 * tr_yyzz_xxy[i] * tbe_0 - 2.0 * tr_xxzz_xxy[i] * tbe_0 - 2.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_xxz[i] = tr_zz_xxz[i] - 2.0 * tr_yzz_xxyz[i] * tke_0 - 2.0 * tr_yyzz_xxz[i] * tbe_0 - 2.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_xyy[i] = tr_zz_xyy[i] + 2.0 * tr_yzz_xy[i] - 2.0 * tr_yzz_xyyy[i] * tke_0 - 2.0 * tr_yyzz_xyy[i] * tbe_0 - 2.0 * tr_xxzz_xyy[i] * tbe_0 - 4.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_xyz[i] = tr_zz_xyz[i] + tr_yzz_xz[i] - 2.0 * tr_yzz_xyyz[i] * tke_0 - 2.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_xzz[i] = tr_zz_xzz[i] - 2.0 * tr_yzz_xyzz[i] * tke_0 - 2.0 * tr_yyzz_xzz[i] * tbe_0 - 2.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_yyy[i] = tr_zz_yyy[i] + 3.0 * tr_yzz_yy[i] - 2.0 * tr_yzz_yyyy[i] * tke_0 - 2.0 * tr_yyzz_yyy[i] * tbe_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 - 6.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_yyz[i] = tr_zz_yyz[i] + 2.0 * tr_yzz_yz[i] - 2.0 * tr_yzz_yyyz[i] * tke_0 - 2.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_xxzz_yyz[i] * tbe_0 - 4.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_yzz[i] = tr_zz_yzz[i] + tr_yzz_zz[i] - 2.0 * tr_yzz_yyzz[i] * tke_0 - 2.0 * tr_yyzz_yzz[i] * tbe_0 - 2.0 * tr_xxzz_yzz[i] * tbe_0 - 2.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyzz_zzz[i] = tr_zz_zzz[i] - 2.0 * tr_yzz_yzzz[i] * tke_0 - 2.0 * tr_yyzz_zzz[i] * tbe_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-250 components of targeted buffer : GF

    auto tr_x_0_y_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 240);

    auto tr_x_0_y_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 241);

    auto tr_x_0_y_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 242);

    auto tr_x_0_y_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 243);

    auto tr_x_0_y_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 244);

    auto tr_x_0_y_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 245);

    auto tr_x_0_y_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 246);

    auto tr_x_0_y_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 247);

    auto tr_x_0_y_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 248);

    auto tr_x_0_y_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 249);

    #pragma omp simd aligned(tr_x_0_y_xzzz_xxx, tr_x_0_y_xzzz_xxy, tr_x_0_y_xzzz_xxz, tr_x_0_y_xzzz_xyy, tr_x_0_y_xzzz_xyz, tr_x_0_y_xzzz_xzz, tr_x_0_y_xzzz_yyy, tr_x_0_y_xzzz_yyz, tr_x_0_y_xzzz_yzz, tr_x_0_y_xzzz_zzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_xxzzz_xx, tr_xxzzz_xxxy, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xy, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yz, tr_xxzzz_yzzz, tr_xxzzz_zz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_zzz_xx, tr_zzz_xxxy, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_yy, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xzzz_xxx[i] = -2.0 * tr_zzz_xxxy[i] * tke_0 - 2.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_xxy[i] = tr_zzz_xx[i] - 2.0 * tr_zzz_xxyy[i] * tke_0 - 2.0 * tr_yzzz_xxy[i] * tbe_0 - 2.0 * tr_xxzzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_xxz[i] = -2.0 * tr_zzz_xxyz[i] * tke_0 - 2.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_xyy[i] = 2.0 * tr_zzz_xy[i] - 2.0 * tr_zzz_xyyy[i] * tke_0 - 2.0 * tr_yzzz_xyy[i] * tbe_0 - 4.0 * tr_xxzzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_xyz[i] = tr_zzz_xz[i] - 2.0 * tr_zzz_xyyz[i] * tke_0 - 2.0 * tr_yzzz_xyz[i] * tbe_0 - 2.0 * tr_xxzzz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_xzz[i] = -2.0 * tr_zzz_xyzz[i] * tke_0 - 2.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_yyy[i] = 3.0 * tr_zzz_yy[i] - 2.0 * tr_zzz_yyyy[i] * tke_0 - 2.0 * tr_yzzz_yyy[i] * tbe_0 - 6.0 * tr_xxzzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_yyz[i] = 2.0 * tr_zzz_yz[i] - 2.0 * tr_zzz_yyyz[i] * tke_0 - 2.0 * tr_yzzz_yyz[i] * tbe_0 - 4.0 * tr_xxzzz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_yzz[i] = tr_zzz_zz[i] - 2.0 * tr_zzz_yyzz[i] * tke_0 - 2.0 * tr_yzzz_yzz[i] * tbe_0 - 2.0 * tr_xxzzz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzzz_zzz[i] = -2.0 * tr_zzz_yzzz[i] * tke_0 - 2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 250-260 components of targeted buffer : GF

    auto tr_x_0_y_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 250);

    auto tr_x_0_y_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 251);

    auto tr_x_0_y_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 252);

    auto tr_x_0_y_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 253);

    auto tr_x_0_y_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 254);

    auto tr_x_0_y_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 255);

    auto tr_x_0_y_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 256);

    auto tr_x_0_y_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 257);

    auto tr_x_0_y_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 258);

    auto tr_x_0_y_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 259);

    #pragma omp simd aligned(tr_x_0_y_yyyy_xxx, tr_x_0_y_yyyy_xxy, tr_x_0_y_yyyy_xxz, tr_x_0_y_yyyy_xyy, tr_x_0_y_yyyy_xyz, tr_x_0_y_yyyy_xzz, tr_x_0_y_yyyy_yyy, tr_x_0_y_yyyy_yyz, tr_x_0_y_yyyy_yzz, tr_x_0_y_yyyy_zzz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyy_xx, tr_xyyyy_xxxy, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xy, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yz, tr_xyyyy_yzzz, tr_xyyyy_zz, tr_xyyyyy_xxx, tr_xyyyyy_xxy, tr_xyyyyy_xxz, tr_xyyyyy_xyy, tr_xyyyyy_xyz, tr_xyyyyy_xzz, tr_xyyyyy_yyy, tr_xyyyyy_yyz, tr_xyyyyy_yzz, tr_xyyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyyy_xxx[i] = -8.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_xxy[i] = -8.0 * tr_xyyy_xxy[i] * tbe_0 - 2.0 * tr_xyyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_xxz[i] = -8.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_xyy[i] = -8.0 * tr_xyyy_xyy[i] * tbe_0 - 4.0 * tr_xyyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_xyz[i] = -8.0 * tr_xyyy_xyz[i] * tbe_0 - 2.0 * tr_xyyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_xzz[i] = -8.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_yyy[i] = -8.0 * tr_xyyy_yyy[i] * tbe_0 - 6.0 * tr_xyyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_yyz[i] = -8.0 * tr_xyyy_yyz[i] * tbe_0 - 4.0 * tr_xyyyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_yzz[i] = -8.0 * tr_xyyy_yzz[i] * tbe_0 - 2.0 * tr_xyyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyy_zzz[i] = -8.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 260-270 components of targeted buffer : GF

    auto tr_x_0_y_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 260);

    auto tr_x_0_y_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 261);

    auto tr_x_0_y_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 262);

    auto tr_x_0_y_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 263);

    auto tr_x_0_y_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 264);

    auto tr_x_0_y_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 265);

    auto tr_x_0_y_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 266);

    auto tr_x_0_y_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 267);

    auto tr_x_0_y_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 268);

    auto tr_x_0_y_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 269);

    #pragma omp simd aligned(tr_x_0_y_yyyz_xxx, tr_x_0_y_yyyz_xxy, tr_x_0_y_yyyz_xxz, tr_x_0_y_yyyz_xyy, tr_x_0_y_yyyz_xyz, tr_x_0_y_yyyz_xzz, tr_x_0_y_yyyz_yyy, tr_x_0_y_yyyz_yyz, tr_x_0_y_yyyz_yzz, tr_x_0_y_yyyz_zzz, tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, tr_xyyyz_xx, tr_xyyyz_xxxy, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xy, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yz, tr_xyyyz_yzzz, tr_xyyyz_zz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyyz_xxx[i] = -6.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_xxy[i] = -6.0 * tr_xyyz_xxy[i] * tbe_0 - 2.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_xxz[i] = -6.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_xyy[i] = -6.0 * tr_xyyz_xyy[i] * tbe_0 - 4.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_xyz[i] = -6.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_xzz[i] = -6.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_yyy[i] = -6.0 * tr_xyyz_yyy[i] * tbe_0 - 6.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_yyz[i] = -6.0 * tr_xyyz_yyz[i] * tbe_0 - 4.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_yzz[i] = -6.0 * tr_xyyz_yzz[i] * tbe_0 - 2.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyyz_zzz[i] = -6.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-280 components of targeted buffer : GF

    auto tr_x_0_y_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 270);

    auto tr_x_0_y_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 271);

    auto tr_x_0_y_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 272);

    auto tr_x_0_y_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 273);

    auto tr_x_0_y_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 274);

    auto tr_x_0_y_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 275);

    auto tr_x_0_y_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 276);

    auto tr_x_0_y_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 277);

    auto tr_x_0_y_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 278);

    auto tr_x_0_y_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 279);

    #pragma omp simd aligned(tr_x_0_y_yyzz_xxx, tr_x_0_y_yyzz_xxy, tr_x_0_y_yyzz_xxz, tr_x_0_y_yyzz_xyy, tr_x_0_y_yyzz_xyz, tr_x_0_y_yyzz_xzz, tr_x_0_y_yyzz_yyy, tr_x_0_y_yyzz_yyz, tr_x_0_y_yyzz_yzz, tr_x_0_y_yyzz_zzz, tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxy, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xy, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yz, tr_xyyzz_yzzz, tr_xyyzz_zz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyzz_xxx[i] = -4.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_xxy[i] = -4.0 * tr_xyzz_xxy[i] * tbe_0 - 2.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_xxz[i] = -4.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_xyy[i] = -4.0 * tr_xyzz_xyy[i] * tbe_0 - 4.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_xyz[i] = -4.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_xzz[i] = -4.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_yyy[i] = -4.0 * tr_xyzz_yyy[i] * tbe_0 - 6.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_yyz[i] = -4.0 * tr_xyzz_yyz[i] * tbe_0 - 4.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_yzz[i] = -4.0 * tr_xyzz_yzz[i] * tbe_0 - 2.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyzz_zzz[i] = -4.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 280-290 components of targeted buffer : GF

    auto tr_x_0_y_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 280);

    auto tr_x_0_y_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 281);

    auto tr_x_0_y_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 282);

    auto tr_x_0_y_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 283);

    auto tr_x_0_y_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 284);

    auto tr_x_0_y_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 285);

    auto tr_x_0_y_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 286);

    auto tr_x_0_y_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 287);

    auto tr_x_0_y_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 288);

    auto tr_x_0_y_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 289);

    #pragma omp simd aligned(tr_x_0_y_yzzz_xxx, tr_x_0_y_yzzz_xxy, tr_x_0_y_yzzz_xxz, tr_x_0_y_yzzz_xyy, tr_x_0_y_yzzz_xyz, tr_x_0_y_yzzz_xzz, tr_x_0_y_yzzz_yyy, tr_x_0_y_yzzz_yyz, tr_x_0_y_yzzz_yzz, tr_x_0_y_yzzz_zzz, tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxy, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xy, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yz, tr_xyzzz_yzzz, tr_xyzzz_zz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yzzz_xxx[i] = -2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_xxy[i] = -2.0 * tr_xzzz_xxy[i] * tbe_0 - 2.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_xxz[i] = -2.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_xyy[i] = -2.0 * tr_xzzz_xyy[i] * tbe_0 - 4.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_xyz[i] = -2.0 * tr_xzzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_xzz[i] = -2.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_yyy[i] = -2.0 * tr_xzzz_yyy[i] * tbe_0 - 6.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_yyz[i] = -2.0 * tr_xzzz_yyz[i] * tbe_0 - 4.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_yzz[i] = -2.0 * tr_xzzz_yzz[i] * tbe_0 - 2.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzzz_zzz[i] = -2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 290-300 components of targeted buffer : GF

    auto tr_x_0_y_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 290);

    auto tr_x_0_y_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 291);

    auto tr_x_0_y_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 292);

    auto tr_x_0_y_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 293);

    auto tr_x_0_y_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 294);

    auto tr_x_0_y_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 295);

    auto tr_x_0_y_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 296);

    auto tr_x_0_y_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 297);

    auto tr_x_0_y_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 298);

    auto tr_x_0_y_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 299);

    #pragma omp simd aligned(tr_x_0_y_zzzz_xxx, tr_x_0_y_zzzz_xxy, tr_x_0_y_zzzz_xxz, tr_x_0_y_zzzz_xyy, tr_x_0_y_zzzz_xyz, tr_x_0_y_zzzz_xzz, tr_x_0_y_zzzz_yyy, tr_x_0_y_zzzz_yyz, tr_x_0_y_zzzz_yzz, tr_x_0_y_zzzz_zzz, tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xxxy, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xy, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yz, tr_xzzzz_yzzz, tr_xzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zzzz_xxx[i] = 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_xxy[i] = -2.0 * tr_xzzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_xxz[i] = 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_xyy[i] = -4.0 * tr_xzzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_xyz[i] = -2.0 * tr_xzzzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_xzz[i] = 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_yyy[i] = -6.0 * tr_xzzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_yyz[i] = -4.0 * tr_xzzzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_yzz[i] = -2.0 * tr_xzzzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzzz_zzz[i] = 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-310 components of targeted buffer : GF

    auto tr_x_0_z_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 300);

    auto tr_x_0_z_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 301);

    auto tr_x_0_z_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 302);

    auto tr_x_0_z_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 303);

    auto tr_x_0_z_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 304);

    auto tr_x_0_z_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 305);

    auto tr_x_0_z_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 306);

    auto tr_x_0_z_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 307);

    auto tr_x_0_z_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 308);

    auto tr_x_0_z_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 309);

    #pragma omp simd aligned(tr_x_0_z_xxxx_xxx, tr_x_0_z_xxxx_xxy, tr_x_0_z_xxxx_xxz, tr_x_0_z_xxxx_xyy, tr_x_0_z_xxxx_xyz, tr_x_0_z_xxxx_xzz, tr_x_0_z_xxxx_yyy, tr_x_0_z_xxxx_yyz, tr_x_0_z_xxxx_yzz, tr_x_0_z_xxxx_zzz, tr_xxx_xx, tr_xxx_xxxz, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxx_zzzz, tr_xxxxx_xx, tr_xxxxx_xxxz, tr_xxxxx_xxyz, tr_xxxxx_xxzz, tr_xxxxx_xy, tr_xxxxx_xyyz, tr_xxxxx_xyzz, tr_xxxxx_xz, tr_xxxxx_xzzz, tr_xxxxx_yy, tr_xxxxx_yyyz, tr_xxxxx_yyzz, tr_xxxxx_yz, tr_xxxxx_yzzz, tr_xxxxx_zz, tr_xxxxx_zzzz, tr_xxxxxz_xxx, tr_xxxxxz_xxy, tr_xxxxxz_xxz, tr_xxxxxz_xyy, tr_xxxxxz_xyz, tr_xxxxxz_xzz, tr_xxxxxz_yyy, tr_xxxxxz_yyz, tr_xxxxxz_yzz, tr_xxxxxz_zzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxxx_xxx[i] = -8.0 * tr_xxx_xxxz[i] * tke_0 - 8.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxxx_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_xxy[i] = -8.0 * tr_xxx_xxyz[i] * tke_0 - 8.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxxx_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_xxz[i] = 4.0 * tr_xxx_xx[i] - 8.0 * tr_xxx_xxzz[i] * tke_0 - 8.0 * tr_xxxz_xxz[i] * tbe_0 - 2.0 * tr_xxxxx_xx[i] * tbe_0 + 4.0 * tr_xxxxx_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_xyy[i] = -8.0 * tr_xxx_xyyz[i] * tke_0 - 8.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxxx_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_xyz[i] = 4.0 * tr_xxx_xy[i] - 8.0 * tr_xxx_xyzz[i] * tke_0 - 8.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxxx_xy[i] * tbe_0 + 4.0 * tr_xxxxx_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_xzz[i] = 8.0 * tr_xxx_xz[i] - 8.0 * tr_xxx_xzzz[i] * tke_0 - 8.0 * tr_xxxz_xzz[i] * tbe_0 - 4.0 * tr_xxxxx_xz[i] * tbe_0 + 4.0 * tr_xxxxx_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_yyy[i] = -8.0 * tr_xxx_yyyz[i] * tke_0 - 8.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxx_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_yyz[i] = 4.0 * tr_xxx_yy[i] - 8.0 * tr_xxx_yyzz[i] * tke_0 - 8.0 * tr_xxxz_yyz[i] * tbe_0 - 2.0 * tr_xxxxx_yy[i] * tbe_0 + 4.0 * tr_xxxxx_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_yzz[i] = 8.0 * tr_xxx_yz[i] - 8.0 * tr_xxx_yzzz[i] * tke_0 - 8.0 * tr_xxxz_yzz[i] * tbe_0 - 4.0 * tr_xxxxx_yz[i] * tbe_0 + 4.0 * tr_xxxxx_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxx_zzz[i] = 12.0 * tr_xxx_zz[i] - 8.0 * tr_xxx_zzzz[i] * tke_0 - 8.0 * tr_xxxz_zzz[i] * tbe_0 - 6.0 * tr_xxxxx_zz[i] * tbe_0 + 4.0 * tr_xxxxx_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 310-320 components of targeted buffer : GF

    auto tr_x_0_z_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 310);

    auto tr_x_0_z_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 311);

    auto tr_x_0_z_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 312);

    auto tr_x_0_z_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 313);

    auto tr_x_0_z_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 314);

    auto tr_x_0_z_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 315);

    auto tr_x_0_z_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 316);

    auto tr_x_0_z_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 317);

    auto tr_x_0_z_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 318);

    auto tr_x_0_z_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 319);

    #pragma omp simd aligned(tr_x_0_z_xxxy_xxx, tr_x_0_z_xxxy_xxy, tr_x_0_z_xxxy_xxz, tr_x_0_z_xxxy_xyy, tr_x_0_z_xxxy_xyz, tr_x_0_z_xxxy_xzz, tr_x_0_z_xxxy_yyy, tr_x_0_z_xxxy_yyz, tr_x_0_z_xxxy_yzz, tr_x_0_z_xxxy_zzz, tr_xxxxy_xx, tr_xxxxy_xxxz, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xz, tr_xxxxy_xzzz, tr_xxxxy_yy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yz, tr_xxxxy_yzzz, tr_xxxxy_zz, tr_xxxxy_zzzz, tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_xxy_xx, tr_xxy_xxxz, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxxy_xxx[i] = -6.0 * tr_xxy_xxxz[i] * tke_0 - 6.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_xxy[i] = -6.0 * tr_xxy_xxyz[i] * tke_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_xxz[i] = 3.0 * tr_xxy_xx[i] - 6.0 * tr_xxy_xxzz[i] * tke_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 - 2.0 * tr_xxxxy_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_xyy[i] = -6.0 * tr_xxy_xyyz[i] * tke_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_xyz[i] = 3.0 * tr_xxy_xy[i] - 6.0 * tr_xxy_xyzz[i] * tke_0 - 6.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxxxy_xy[i] * tbe_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_xzz[i] = 6.0 * tr_xxy_xz[i] - 6.0 * tr_xxy_xzzz[i] * tke_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 - 4.0 * tr_xxxxy_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_yyy[i] = -6.0 * tr_xxy_yyyz[i] * tke_0 - 6.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_yyz[i] = 3.0 * tr_xxy_yy[i] - 6.0 * tr_xxy_yyzz[i] * tke_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 - 2.0 * tr_xxxxy_yy[i] * tbe_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_yzz[i] = 6.0 * tr_xxy_yz[i] - 6.0 * tr_xxy_yzzz[i] * tke_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 - 4.0 * tr_xxxxy_yz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxy_zzz[i] = 9.0 * tr_xxy_zz[i] - 6.0 * tr_xxy_zzzz[i] * tke_0 - 6.0 * tr_xxyz_zzz[i] * tbe_0 - 6.0 * tr_xxxxy_zz[i] * tbe_0 + 4.0 * tr_xxxxy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 320-330 components of targeted buffer : GF

    auto tr_x_0_z_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 320);

    auto tr_x_0_z_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 321);

    auto tr_x_0_z_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 322);

    auto tr_x_0_z_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 323);

    auto tr_x_0_z_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 324);

    auto tr_x_0_z_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 325);

    auto tr_x_0_z_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 326);

    auto tr_x_0_z_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 327);

    auto tr_x_0_z_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 328);

    auto tr_x_0_z_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 329);

    #pragma omp simd aligned(tr_x_0_z_xxxz_xxx, tr_x_0_z_xxxz_xxy, tr_x_0_z_xxxz_xxz, tr_x_0_z_xxxz_xyy, tr_x_0_z_xxxz_xyz, tr_x_0_z_xxxz_xzz, tr_x_0_z_xxxz_yyy, tr_x_0_z_xxxz_yyz, tr_x_0_z_xxxz_yzz, tr_x_0_z_xxxz_zzz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxz_xx, tr_xxxxz_xxxz, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xz, tr_xxxxz_xzzz, tr_xxxxz_yy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yz, tr_xxxxz_yzzz, tr_xxxxz_zz, tr_xxxxz_zzzz, tr_xxxxzz_xxx, tr_xxxxzz_xxy, tr_xxxxzz_xxz, tr_xxxxzz_xyy, tr_xxxxzz_xyz, tr_xxxxzz_xzz, tr_xxxxzz_yyy, tr_xxxxzz_yyz, tr_xxxxzz_yzz, tr_xxxxzz_zzz, tr_xxz_xx, tr_xxz_xxxz, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxxz_xxx[i] = 3.0 * tr_xx_xxx[i] - 6.0 * tr_xxz_xxxz[i] * tke_0 - 6.0 * tr_xxzz_xxx[i] * tbe_0 - 2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_xxy[i] = 3.0 * tr_xx_xxy[i] - 6.0 * tr_xxz_xxyz[i] * tke_0 - 6.0 * tr_xxzz_xxy[i] * tbe_0 - 2.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_xxz[i] = 3.0 * tr_xx_xxz[i] + 3.0 * tr_xxz_xx[i] - 6.0 * tr_xxz_xxzz[i] * tke_0 - 6.0 * tr_xxzz_xxz[i] * tbe_0 - 2.0 * tr_xxxx_xxz[i] * tbe_0 - 2.0 * tr_xxxxz_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_xyy[i] = 3.0 * tr_xx_xyy[i] - 6.0 * tr_xxz_xyyz[i] * tke_0 - 6.0 * tr_xxzz_xyy[i] * tbe_0 - 2.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_xyz[i] = 3.0 * tr_xx_xyz[i] + 3.0 * tr_xxz_xy[i] - 6.0 * tr_xxz_xyzz[i] * tke_0 - 6.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxxx_xyz[i] * tbe_0 - 2.0 * tr_xxxxz_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_xzz[i] = 3.0 * tr_xx_xzz[i] + 6.0 * tr_xxz_xz[i] - 6.0 * tr_xxz_xzzz[i] * tke_0 - 6.0 * tr_xxzz_xzz[i] * tbe_0 - 2.0 * tr_xxxx_xzz[i] * tbe_0 - 4.0 * tr_xxxxz_xz[i] * tbe_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_yyy[i] = 3.0 * tr_xx_yyy[i] - 6.0 * tr_xxz_yyyz[i] * tke_0 - 6.0 * tr_xxzz_yyy[i] * tbe_0 - 2.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_yyz[i] = 3.0 * tr_xx_yyz[i] + 3.0 * tr_xxz_yy[i] - 6.0 * tr_xxz_yyzz[i] * tke_0 - 6.0 * tr_xxzz_yyz[i] * tbe_0 - 2.0 * tr_xxxx_yyz[i] * tbe_0 - 2.0 * tr_xxxxz_yy[i] * tbe_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_yzz[i] = 3.0 * tr_xx_yzz[i] + 6.0 * tr_xxz_yz[i] - 6.0 * tr_xxz_yzzz[i] * tke_0 - 6.0 * tr_xxzz_yzz[i] * tbe_0 - 2.0 * tr_xxxx_yzz[i] * tbe_0 - 4.0 * tr_xxxxz_yz[i] * tbe_0 + 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxxz_zzz[i] = 3.0 * tr_xx_zzz[i] + 9.0 * tr_xxz_zz[i] - 6.0 * tr_xxz_zzzz[i] * tke_0 - 6.0 * tr_xxzz_zzz[i] * tbe_0 - 2.0 * tr_xxxx_zzz[i] * tbe_0 - 6.0 * tr_xxxxz_zz[i] * tbe_0 + 4.0 * tr_xxxxz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-340 components of targeted buffer : GF

    auto tr_x_0_z_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 330);

    auto tr_x_0_z_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 331);

    auto tr_x_0_z_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 332);

    auto tr_x_0_z_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 333);

    auto tr_x_0_z_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 334);

    auto tr_x_0_z_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 335);

    auto tr_x_0_z_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 336);

    auto tr_x_0_z_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 337);

    auto tr_x_0_z_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 338);

    auto tr_x_0_z_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 339);

    #pragma omp simd aligned(tr_x_0_z_xxyy_xxx, tr_x_0_z_xxyy_xxy, tr_x_0_z_xxyy_xxz, tr_x_0_z_xxyy_xyy, tr_x_0_z_xxyy_xyz, tr_x_0_z_xxyy_xzz, tr_x_0_z_xxyy_yyy, tr_x_0_z_xxyy_yyz, tr_x_0_z_xxyy_yzz, tr_x_0_z_xxyy_zzz, tr_xxxyy_xx, tr_xxxyy_xxxz, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xz, tr_xxxyy_xzzz, tr_xxxyy_yy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yz, tr_xxxyy_yzzz, tr_xxxyy_zz, tr_xxxyy_zzzz, tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xyy_xx, tr_xyy_xxxz, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxyy_xxx[i] = -4.0 * tr_xyy_xxxz[i] * tke_0 - 4.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_xxy[i] = -4.0 * tr_xyy_xxyz[i] * tke_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_xxz[i] = 2.0 * tr_xyy_xx[i] - 4.0 * tr_xyy_xxzz[i] * tke_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 - 2.0 * tr_xxxyy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_xyy[i] = -4.0 * tr_xyy_xyyz[i] * tke_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_xyz[i] = 2.0 * tr_xyy_xy[i] - 4.0 * tr_xyy_xyzz[i] * tke_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xxxyy_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_xzz[i] = 4.0 * tr_xyy_xz[i] - 4.0 * tr_xyy_xzzz[i] * tke_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 - 4.0 * tr_xxxyy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_yyy[i] = -4.0 * tr_xyy_yyyz[i] * tke_0 - 4.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_yyz[i] = 2.0 * tr_xyy_yy[i] - 4.0 * tr_xyy_yyzz[i] * tke_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xxxyy_yy[i] * tbe_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_yzz[i] = 4.0 * tr_xyy_yz[i] - 4.0 * tr_xyy_yzzz[i] * tke_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 - 4.0 * tr_xxxyy_yz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyy_zzz[i] = 6.0 * tr_xyy_zz[i] - 4.0 * tr_xyy_zzzz[i] * tke_0 - 4.0 * tr_xyyz_zzz[i] * tbe_0 - 6.0 * tr_xxxyy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 340-350 components of targeted buffer : GF

    auto tr_x_0_z_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 340);

    auto tr_x_0_z_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 341);

    auto tr_x_0_z_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 342);

    auto tr_x_0_z_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 343);

    auto tr_x_0_z_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 344);

    auto tr_x_0_z_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 345);

    auto tr_x_0_z_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 346);

    auto tr_x_0_z_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 347);

    auto tr_x_0_z_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 348);

    auto tr_x_0_z_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 349);

    #pragma omp simd aligned(tr_x_0_z_xxyz_xxx, tr_x_0_z_xxyz_xxy, tr_x_0_z_xxyz_xxz, tr_x_0_z_xxyz_xyy, tr_x_0_z_xxyz_xyz, tr_x_0_z_xxyz_xzz, tr_x_0_z_xxyz_yyy, tr_x_0_z_xxyz_yyz, tr_x_0_z_xxyz_yzz, tr_x_0_z_xxyz_zzz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyz_xx, tr_xxxyz_xxxz, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_xzzz, tr_xxxyz_yy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yz, tr_xxxyz_yzzz, tr_xxxyz_zz, tr_xxxyz_zzzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxyz_xxx[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xyz_xxxz[i] * tke_0 - 4.0 * tr_xyzz_xxx[i] * tbe_0 - 2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_xxy[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 - 2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_xxz[i] = 2.0 * tr_xy_xxz[i] + 2.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxzz[i] * tke_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 - 2.0 * tr_xxxy_xxz[i] * tbe_0 - 2.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_xyy[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 - 2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_xyz[i] = 2.0 * tr_xy_xyz[i] + 2.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xyzz[i] * tke_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_xzz[i] = 2.0 * tr_xy_xzz[i] + 4.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xzzz[i] * tke_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 - 2.0 * tr_xxxy_xzz[i] * tbe_0 - 4.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_yyy[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xyz_yyyz[i] * tke_0 - 4.0 * tr_xyzz_yyy[i] * tbe_0 - 2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_yyz[i] = 2.0 * tr_xy_yyz[i] + 2.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_yyzz[i] * tke_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 - 2.0 * tr_xxxy_yyz[i] * tbe_0 - 2.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_yzz[i] = 2.0 * tr_xy_yzz[i] + 4.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_yzzz[i] * tke_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 - 2.0 * tr_xxxy_yzz[i] * tbe_0 - 4.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxyz_zzz[i] = 2.0 * tr_xy_zzz[i] + 6.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_zzzz[i] * tke_0 - 4.0 * tr_xyzz_zzz[i] * tbe_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 - 6.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 350-360 components of targeted buffer : GF

    auto tr_x_0_z_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 350);

    auto tr_x_0_z_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 351);

    auto tr_x_0_z_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 352);

    auto tr_x_0_z_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 353);

    auto tr_x_0_z_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 354);

    auto tr_x_0_z_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 355);

    auto tr_x_0_z_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 356);

    auto tr_x_0_z_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 357);

    auto tr_x_0_z_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 358);

    auto tr_x_0_z_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 359);

    #pragma omp simd aligned(tr_x_0_z_xxzz_xxx, tr_x_0_z_xxzz_xxy, tr_x_0_z_xxzz_xxz, tr_x_0_z_xxzz_xyy, tr_x_0_z_xxzz_xyz, tr_x_0_z_xxzz_xzz, tr_x_0_z_xxzz_yyy, tr_x_0_z_xxzz_yyz, tr_x_0_z_xxzz_yzz, tr_x_0_z_xxzz_zzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxxzz_xx, tr_xxxzz_xxxz, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xz, tr_xxxzz_xzzz, tr_xxxzz_yy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yz, tr_xxxzz_yzzz, tr_xxxzz_zz, tr_xxxzz_zzzz, tr_xxxzzz_xxx, tr_xxxzzz_xxy, tr_xxxzzz_xxz, tr_xxxzzz_xyy, tr_xxxzzz_xyz, tr_xxxzzz_xzz, tr_xxxzzz_yyy, tr_xxxzzz_yyz, tr_xxxzzz_yzz, tr_xxxzzz_zzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxz, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxzz_xxx[i] = 4.0 * tr_xz_xxx[i] - 4.0 * tr_xzz_xxxz[i] * tke_0 - 4.0 * tr_xzzz_xxx[i] * tbe_0 - 4.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_xxy[i] = 4.0 * tr_xz_xxy[i] - 4.0 * tr_xzz_xxyz[i] * tke_0 - 4.0 * tr_xzzz_xxy[i] * tbe_0 - 4.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_xxz[i] = 4.0 * tr_xz_xxz[i] + 2.0 * tr_xzz_xx[i] - 4.0 * tr_xzz_xxzz[i] * tke_0 - 4.0 * tr_xzzz_xxz[i] * tbe_0 - 4.0 * tr_xxxz_xxz[i] * tbe_0 - 2.0 * tr_xxxzz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_xyy[i] = 4.0 * tr_xz_xyy[i] - 4.0 * tr_xzz_xyyz[i] * tke_0 - 4.0 * tr_xzzz_xyy[i] * tbe_0 - 4.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_xyz[i] = 4.0 * tr_xz_xyz[i] + 2.0 * tr_xzz_xy[i] - 4.0 * tr_xzz_xyzz[i] * tke_0 - 4.0 * tr_xzzz_xyz[i] * tbe_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxzz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_xzz[i] = 4.0 * tr_xz_xzz[i] + 4.0 * tr_xzz_xz[i] - 4.0 * tr_xzz_xzzz[i] * tke_0 - 4.0 * tr_xzzz_xzz[i] * tbe_0 - 4.0 * tr_xxxz_xzz[i] * tbe_0 - 4.0 * tr_xxxzz_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_yyy[i] = 4.0 * tr_xz_yyy[i] - 4.0 * tr_xzz_yyyz[i] * tke_0 - 4.0 * tr_xzzz_yyy[i] * tbe_0 - 4.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_yyz[i] = 4.0 * tr_xz_yyz[i] + 2.0 * tr_xzz_yy[i] - 4.0 * tr_xzz_yyzz[i] * tke_0 - 4.0 * tr_xzzz_yyz[i] * tbe_0 - 4.0 * tr_xxxz_yyz[i] * tbe_0 - 2.0 * tr_xxxzz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_yzz[i] = 4.0 * tr_xz_yzz[i] + 4.0 * tr_xzz_yz[i] - 4.0 * tr_xzz_yzzz[i] * tke_0 - 4.0 * tr_xzzz_yzz[i] * tbe_0 - 4.0 * tr_xxxz_yzz[i] * tbe_0 - 4.0 * tr_xxxzz_yz[i] * tbe_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxzz_zzz[i] = 4.0 * tr_xz_zzz[i] + 6.0 * tr_xzz_zz[i] - 4.0 * tr_xzz_zzzz[i] * tke_0 - 4.0 * tr_xzzz_zzz[i] * tbe_0 - 4.0 * tr_xxxz_zzz[i] * tbe_0 - 6.0 * tr_xxxzz_zz[i] * tbe_0 + 4.0 * tr_xxxzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-370 components of targeted buffer : GF

    auto tr_x_0_z_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 360);

    auto tr_x_0_z_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 361);

    auto tr_x_0_z_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 362);

    auto tr_x_0_z_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 363);

    auto tr_x_0_z_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 364);

    auto tr_x_0_z_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 365);

    auto tr_x_0_z_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 366);

    auto tr_x_0_z_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 367);

    auto tr_x_0_z_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 368);

    auto tr_x_0_z_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 369);

    #pragma omp simd aligned(tr_x_0_z_xyyy_xxx, tr_x_0_z_xyyy_xxy, tr_x_0_z_xyyy_xxz, tr_x_0_z_xyyy_xyy, tr_x_0_z_xyyy_xyz, tr_x_0_z_xyyy_xzz, tr_x_0_z_xyyy_yyy, tr_x_0_z_xyyy_yyz, tr_x_0_z_xyyy_yzz, tr_x_0_z_xyyy_zzz, tr_xxyyy_xx, tr_xxyyy_xxxz, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xz, tr_xxyyy_xzzz, tr_xxyyy_yy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yz, tr_xxyyy_yzzz, tr_xxyyy_zz, tr_xxyyy_zzzz, tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_yyy_xx, tr_yyy_xxxz, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyyy_xxx[i] = -2.0 * tr_yyy_xxxz[i] * tke_0 - 2.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_xxy[i] = -2.0 * tr_yyy_xxyz[i] * tke_0 - 2.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_xxz[i] = tr_yyy_xx[i] - 2.0 * tr_yyy_xxzz[i] * tke_0 - 2.0 * tr_yyyz_xxz[i] * tbe_0 - 2.0 * tr_xxyyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_xyy[i] = -2.0 * tr_yyy_xyyz[i] * tke_0 - 2.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_xyz[i] = tr_yyy_xy[i] - 2.0 * tr_yyy_xyzz[i] * tke_0 - 2.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_xzz[i] = 2.0 * tr_yyy_xz[i] - 2.0 * tr_yyy_xzzz[i] * tke_0 - 2.0 * tr_yyyz_xzz[i] * tbe_0 - 4.0 * tr_xxyyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_yyy[i] = -2.0 * tr_yyy_yyyz[i] * tke_0 - 2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_yyz[i] = tr_yyy_yy[i] - 2.0 * tr_yyy_yyzz[i] * tke_0 - 2.0 * tr_yyyz_yyz[i] * tbe_0 - 2.0 * tr_xxyyy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_yzz[i] = 2.0 * tr_yyy_yz[i] - 2.0 * tr_yyy_yzzz[i] * tke_0 - 2.0 * tr_yyyz_yzz[i] * tbe_0 - 4.0 * tr_xxyyy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyy_zzz[i] = 3.0 * tr_yyy_zz[i] - 2.0 * tr_yyy_zzzz[i] * tke_0 - 2.0 * tr_yyyz_zzz[i] * tbe_0 - 6.0 * tr_xxyyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 370-380 components of targeted buffer : GF

    auto tr_x_0_z_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 370);

    auto tr_x_0_z_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 371);

    auto tr_x_0_z_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 372);

    auto tr_x_0_z_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 373);

    auto tr_x_0_z_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 374);

    auto tr_x_0_z_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 375);

    auto tr_x_0_z_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 376);

    auto tr_x_0_z_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 377);

    auto tr_x_0_z_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 378);

    auto tr_x_0_z_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 379);

    #pragma omp simd aligned(tr_x_0_z_xyyz_xxx, tr_x_0_z_xyyz_xxy, tr_x_0_z_xyyz_xxz, tr_x_0_z_xyyz_xyy, tr_x_0_z_xyyz_xyz, tr_x_0_z_xyyz_xzz, tr_x_0_z_xyyz_yyy, tr_x_0_z_xyyz_yyz, tr_x_0_z_xyyz_yzz, tr_x_0_z_xyyz_zzz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyz_xx, tr_xxyyz_xxxz, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_xzzz, tr_xxyyz_yy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yz, tr_xxyyz_yzzz, tr_xxyyz_zz, tr_xxyyz_zzzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyz_xx, tr_yyz_xxxz, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyyz_xxx[i] = tr_yy_xxx[i] - 2.0 * tr_yyz_xxxz[i] * tke_0 - 2.0 * tr_yyzz_xxx[i] * tbe_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_xxy[i] = tr_yy_xxy[i] - 2.0 * tr_yyz_xxyz[i] * tke_0 - 2.0 * tr_yyzz_xxy[i] * tbe_0 - 2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_xxz[i] = tr_yy_xxz[i] + tr_yyz_xx[i] - 2.0 * tr_yyz_xxzz[i] * tke_0 - 2.0 * tr_yyzz_xxz[i] * tbe_0 - 2.0 * tr_xxyy_xxz[i] * tbe_0 - 2.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_xyy[i] = tr_yy_xyy[i] - 2.0 * tr_yyz_xyyz[i] * tke_0 - 2.0 * tr_yyzz_xyy[i] * tbe_0 - 2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_xyz[i] = tr_yy_xyz[i] + tr_yyz_xy[i] - 2.0 * tr_yyz_xyzz[i] * tke_0 - 2.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_xzz[i] = tr_yy_xzz[i] + 2.0 * tr_yyz_xz[i] - 2.0 * tr_yyz_xzzz[i] * tke_0 - 2.0 * tr_yyzz_xzz[i] * tbe_0 - 2.0 * tr_xxyy_xzz[i] * tbe_0 - 4.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_yyy[i] = tr_yy_yyy[i] - 2.0 * tr_yyz_yyyz[i] * tke_0 - 2.0 * tr_yyzz_yyy[i] * tbe_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_yyz[i] = tr_yy_yyz[i] + tr_yyz_yy[i] - 2.0 * tr_yyz_yyzz[i] * tke_0 - 2.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_xxyy_yyz[i] * tbe_0 - 2.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_yzz[i] = tr_yy_yzz[i] + 2.0 * tr_yyz_yz[i] - 2.0 * tr_yyz_yzzz[i] * tke_0 - 2.0 * tr_yyzz_yzz[i] * tbe_0 - 2.0 * tr_xxyy_yzz[i] * tbe_0 - 4.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyyz_zzz[i] = tr_yy_zzz[i] + 3.0 * tr_yyz_zz[i] - 2.0 * tr_yyz_zzzz[i] * tke_0 - 2.0 * tr_yyzz_zzz[i] * tbe_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 - 6.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 380-390 components of targeted buffer : GF

    auto tr_x_0_z_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 380);

    auto tr_x_0_z_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 381);

    auto tr_x_0_z_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 382);

    auto tr_x_0_z_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 383);

    auto tr_x_0_z_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 384);

    auto tr_x_0_z_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 385);

    auto tr_x_0_z_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 386);

    auto tr_x_0_z_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 387);

    auto tr_x_0_z_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 388);

    auto tr_x_0_z_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 389);

    #pragma omp simd aligned(tr_x_0_z_xyzz_xxx, tr_x_0_z_xyzz_xxy, tr_x_0_z_xyzz_xxz, tr_x_0_z_xyzz_xyy, tr_x_0_z_xyzz_xyz, tr_x_0_z_xyzz_xzz, tr_x_0_z_xyzz_yyy, tr_x_0_z_xyzz_yyz, tr_x_0_z_xyzz_yzz, tr_x_0_z_xyzz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxz, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_xzzz, tr_xxyzz_yy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yz, tr_xxyzz_yzzz, tr_xxyzz_zz, tr_xxyzz_zzzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxz, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyzz_xxx[i] = 2.0 * tr_yz_xxx[i] - 2.0 * tr_yzz_xxxz[i] * tke_0 - 2.0 * tr_yzzz_xxx[i] * tbe_0 - 4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_xxy[i] = 2.0 * tr_yz_xxy[i] - 2.0 * tr_yzz_xxyz[i] * tke_0 - 2.0 * tr_yzzz_xxy[i] * tbe_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_xxz[i] = 2.0 * tr_yz_xxz[i] + tr_yzz_xx[i] - 2.0 * tr_yzz_xxzz[i] * tke_0 - 2.0 * tr_yzzz_xxz[i] * tbe_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 - 2.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_xyy[i] = 2.0 * tr_yz_xyy[i] - 2.0 * tr_yzz_xyyz[i] * tke_0 - 2.0 * tr_yzzz_xyy[i] * tbe_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_xyz[i] = 2.0 * tr_yz_xyz[i] + tr_yzz_xy[i] - 2.0 * tr_yzz_xyzz[i] * tke_0 - 2.0 * tr_yzzz_xyz[i] * tbe_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_xzz[i] = 2.0 * tr_yz_xzz[i] + 2.0 * tr_yzz_xz[i] - 2.0 * tr_yzz_xzzz[i] * tke_0 - 2.0 * tr_yzzz_xzz[i] * tbe_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 - 4.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_yyy[i] = 2.0 * tr_yz_yyy[i] - 2.0 * tr_yzz_yyyz[i] * tke_0 - 2.0 * tr_yzzz_yyy[i] * tbe_0 - 4.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_yyz[i] = 2.0 * tr_yz_yyz[i] + tr_yzz_yy[i] - 2.0 * tr_yzz_yyzz[i] * tke_0 - 2.0 * tr_yzzz_yyz[i] * tbe_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 - 2.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_yzz[i] = 2.0 * tr_yz_yzz[i] + 2.0 * tr_yzz_yz[i] - 2.0 * tr_yzz_yzzz[i] * tke_0 - 2.0 * tr_yzzz_yzz[i] * tbe_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 - 4.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyzz_zzz[i] = 2.0 * tr_yz_zzz[i] + 3.0 * tr_yzz_zz[i] - 2.0 * tr_yzz_zzzz[i] * tke_0 - 2.0 * tr_yzzz_zzz[i] * tbe_0 - 4.0 * tr_xxyz_zzz[i] * tbe_0 - 6.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-400 components of targeted buffer : GF

    auto tr_x_0_z_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 390);

    auto tr_x_0_z_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 391);

    auto tr_x_0_z_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 392);

    auto tr_x_0_z_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 393);

    auto tr_x_0_z_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 394);

    auto tr_x_0_z_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 395);

    auto tr_x_0_z_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 396);

    auto tr_x_0_z_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 397);

    auto tr_x_0_z_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 398);

    auto tr_x_0_z_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 399);

    #pragma omp simd aligned(tr_x_0_z_xzzz_xxx, tr_x_0_z_xzzz_xxy, tr_x_0_z_xzzz_xxz, tr_x_0_z_xzzz_xyy, tr_x_0_z_xzzz_xyz, tr_x_0_z_xzzz_xzz, tr_x_0_z_xzzz_yyy, tr_x_0_z_xzzz_yyz, tr_x_0_z_xzzz_yzz, tr_x_0_z_xzzz_zzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xxzzz_xx, tr_xxzzz_xxxz, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xz, tr_xxzzz_xzzz, tr_xxzzz_yy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yz, tr_xxzzz_yzzz, tr_xxzzz_zz, tr_xxzzz_zzzz, tr_xxzzzz_xxx, tr_xxzzzz_xxy, tr_xxzzzz_xxz, tr_xxzzzz_xyy, tr_xxzzzz_xyz, tr_xxzzzz_xzz, tr_xxzzzz_yyy, tr_xxzzzz_yyz, tr_xxzzzz_yzz, tr_xxzzzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzz_xx, tr_zzz_xxxz, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xzzz_xxx[i] = 3.0 * tr_zz_xxx[i] - 2.0 * tr_zzz_xxxz[i] * tke_0 - 2.0 * tr_zzzz_xxx[i] * tbe_0 - 6.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_xxy[i] = 3.0 * tr_zz_xxy[i] - 2.0 * tr_zzz_xxyz[i] * tke_0 - 2.0 * tr_zzzz_xxy[i] * tbe_0 - 6.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_xxz[i] = 3.0 * tr_zz_xxz[i] + tr_zzz_xx[i] - 2.0 * tr_zzz_xxzz[i] * tke_0 - 2.0 * tr_zzzz_xxz[i] * tbe_0 - 6.0 * tr_xxzz_xxz[i] * tbe_0 - 2.0 * tr_xxzzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_xyy[i] = 3.0 * tr_zz_xyy[i] - 2.0 * tr_zzz_xyyz[i] * tke_0 - 2.0 * tr_zzzz_xyy[i] * tbe_0 - 6.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_xyz[i] = 3.0 * tr_zz_xyz[i] + tr_zzz_xy[i] - 2.0 * tr_zzz_xyzz[i] * tke_0 - 2.0 * tr_zzzz_xyz[i] * tbe_0 - 6.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxzzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_xzz[i] = 3.0 * tr_zz_xzz[i] + 2.0 * tr_zzz_xz[i] - 2.0 * tr_zzz_xzzz[i] * tke_0 - 2.0 * tr_zzzz_xzz[i] * tbe_0 - 6.0 * tr_xxzz_xzz[i] * tbe_0 - 4.0 * tr_xxzzz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_yyy[i] = 3.0 * tr_zz_yyy[i] - 2.0 * tr_zzz_yyyz[i] * tke_0 - 2.0 * tr_zzzz_yyy[i] * tbe_0 - 6.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_yyz[i] = 3.0 * tr_zz_yyz[i] + tr_zzz_yy[i] - 2.0 * tr_zzz_yyzz[i] * tke_0 - 2.0 * tr_zzzz_yyz[i] * tbe_0 - 6.0 * tr_xxzz_yyz[i] * tbe_0 - 2.0 * tr_xxzzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_yzz[i] = 3.0 * tr_zz_yzz[i] + 2.0 * tr_zzz_yz[i] - 2.0 * tr_zzz_yzzz[i] * tke_0 - 2.0 * tr_zzzz_yzz[i] * tbe_0 - 6.0 * tr_xxzz_yzz[i] * tbe_0 - 4.0 * tr_xxzzz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzzz_zzz[i] = 3.0 * tr_zz_zzz[i] + 3.0 * tr_zzz_zz[i] - 2.0 * tr_zzz_zzzz[i] * tke_0 - 2.0 * tr_zzzz_zzz[i] * tbe_0 - 6.0 * tr_xxzz_zzz[i] * tbe_0 - 6.0 * tr_xxzzz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 400-410 components of targeted buffer : GF

    auto tr_x_0_z_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 400);

    auto tr_x_0_z_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 401);

    auto tr_x_0_z_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 402);

    auto tr_x_0_z_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 403);

    auto tr_x_0_z_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 404);

    auto tr_x_0_z_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 405);

    auto tr_x_0_z_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 406);

    auto tr_x_0_z_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 407);

    auto tr_x_0_z_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 408);

    auto tr_x_0_z_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 409);

    #pragma omp simd aligned(tr_x_0_z_yyyy_xxx, tr_x_0_z_yyyy_xxy, tr_x_0_z_yyyy_xxz, tr_x_0_z_yyyy_xyy, tr_x_0_z_yyyy_xyz, tr_x_0_z_yyyy_xzz, tr_x_0_z_yyyy_yyy, tr_x_0_z_yyyy_yyz, tr_x_0_z_yyyy_yzz, tr_x_0_z_yyyy_zzz, tr_xyyyy_xx, tr_xyyyy_xxxz, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xz, tr_xyyyy_xzzz, tr_xyyyy_yy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yz, tr_xyyyy_yzzz, tr_xyyyy_zz, tr_xyyyy_zzzz, tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyyy_xxx[i] = 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_xxy[i] = 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_xxz[i] = -2.0 * tr_xyyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_xyy[i] = 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_xyz[i] = -2.0 * tr_xyyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_xzz[i] = -4.0 * tr_xyyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_yyy[i] = 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_yyz[i] = -2.0 * tr_xyyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_yzz[i] = -4.0 * tr_xyyyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyy_zzz[i] = -6.0 * tr_xyyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 410-420 components of targeted buffer : GF

    auto tr_x_0_z_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 410);

    auto tr_x_0_z_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 411);

    auto tr_x_0_z_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 412);

    auto tr_x_0_z_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 413);

    auto tr_x_0_z_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 414);

    auto tr_x_0_z_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 415);

    auto tr_x_0_z_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 416);

    auto tr_x_0_z_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 417);

    auto tr_x_0_z_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 418);

    auto tr_x_0_z_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 419);

    #pragma omp simd aligned(tr_x_0_z_yyyz_xxx, tr_x_0_z_yyyz_xxy, tr_x_0_z_yyyz_xxz, tr_x_0_z_yyyz_xyy, tr_x_0_z_yyyz_xyz, tr_x_0_z_yyyz_xzz, tr_x_0_z_yyyz_yyy, tr_x_0_z_yyyz_yyz, tr_x_0_z_yyyz_yzz, tr_x_0_z_yyyz_zzz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyz_xx, tr_xyyyz_xxxz, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_xzzz, tr_xyyyz_yy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yz, tr_xyyyz_yzzz, tr_xyyyz_zz, tr_xyyyz_zzzz, tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyyz_xxx[i] = -2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_xxy[i] = -2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_xxz[i] = -2.0 * tr_xyyy_xxz[i] * tbe_0 - 2.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_xyy[i] = -2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_xyz[i] = -2.0 * tr_xyyy_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_xzz[i] = -2.0 * tr_xyyy_xzz[i] * tbe_0 - 4.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_yyy[i] = -2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_yyz[i] = -2.0 * tr_xyyy_yyz[i] * tbe_0 - 2.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_yzz[i] = -2.0 * tr_xyyy_yzz[i] * tbe_0 - 4.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyyz_zzz[i] = -2.0 * tr_xyyy_zzz[i] * tbe_0 - 6.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-430 components of targeted buffer : GF

    auto tr_x_0_z_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 420);

    auto tr_x_0_z_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 421);

    auto tr_x_0_z_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 422);

    auto tr_x_0_z_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 423);

    auto tr_x_0_z_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 424);

    auto tr_x_0_z_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 425);

    auto tr_x_0_z_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 426);

    auto tr_x_0_z_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 427);

    auto tr_x_0_z_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 428);

    auto tr_x_0_z_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 429);

    #pragma omp simd aligned(tr_x_0_z_yyzz_xxx, tr_x_0_z_yyzz_xxy, tr_x_0_z_yyzz_xxz, tr_x_0_z_yyzz_xyy, tr_x_0_z_yyzz_xyz, tr_x_0_z_yyzz_xzz, tr_x_0_z_yyzz_yyy, tr_x_0_z_yyzz_yyz, tr_x_0_z_yyzz_yzz, tr_x_0_z_yyzz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxz, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_xzzz, tr_xyyzz_yy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yz, tr_xyyzz_yzzz, tr_xyyzz_zz, tr_xyyzz_zzzz, tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyzz_xxx[i] = -4.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_xxy[i] = -4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_xxz[i] = -4.0 * tr_xyyz_xxz[i] * tbe_0 - 2.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_xyy[i] = -4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_xyz[i] = -4.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_xzz[i] = -4.0 * tr_xyyz_xzz[i] * tbe_0 - 4.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_yyy[i] = -4.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_yyz[i] = -4.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_yzz[i] = -4.0 * tr_xyyz_yzz[i] * tbe_0 - 4.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyzz_zzz[i] = -4.0 * tr_xyyz_zzz[i] * tbe_0 - 6.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 430-440 components of targeted buffer : GF

    auto tr_x_0_z_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 430);

    auto tr_x_0_z_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 431);

    auto tr_x_0_z_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 432);

    auto tr_x_0_z_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 433);

    auto tr_x_0_z_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 434);

    auto tr_x_0_z_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 435);

    auto tr_x_0_z_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 436);

    auto tr_x_0_z_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 437);

    auto tr_x_0_z_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 438);

    auto tr_x_0_z_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 439);

    #pragma omp simd aligned(tr_x_0_z_yzzz_xxx, tr_x_0_z_yzzz_xxy, tr_x_0_z_yzzz_xxz, tr_x_0_z_yzzz_xyy, tr_x_0_z_yzzz_xyz, tr_x_0_z_yzzz_xzz, tr_x_0_z_yzzz_yyy, tr_x_0_z_yzzz_yyz, tr_x_0_z_yzzz_yzz, tr_x_0_z_yzzz_zzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxz, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_xzzz, tr_xyzzz_yy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yz, tr_xyzzz_yzzz, tr_xyzzz_zz, tr_xyzzz_zzzz, tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yzzz_xxx[i] = -6.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_xxy[i] = -6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_xxz[i] = -6.0 * tr_xyzz_xxz[i] * tbe_0 - 2.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_xyy[i] = -6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_xyz[i] = -6.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_xzz[i] = -6.0 * tr_xyzz_xzz[i] * tbe_0 - 4.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_yyy[i] = -6.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_yyz[i] = -6.0 * tr_xyzz_yyz[i] * tbe_0 - 2.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_yzz[i] = -6.0 * tr_xyzz_yzz[i] * tbe_0 - 4.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzzz_zzz[i] = -6.0 * tr_xyzz_zzz[i] * tbe_0 - 6.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 440-450 components of targeted buffer : GF

    auto tr_x_0_z_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 440);

    auto tr_x_0_z_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 441);

    auto tr_x_0_z_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 442);

    auto tr_x_0_z_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 443);

    auto tr_x_0_z_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 444);

    auto tr_x_0_z_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 445);

    auto tr_x_0_z_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 446);

    auto tr_x_0_z_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 447);

    auto tr_x_0_z_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 448);

    auto tr_x_0_z_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 449);

    #pragma omp simd aligned(tr_x_0_z_zzzz_xxx, tr_x_0_z_zzzz_xxy, tr_x_0_z_zzzz_xxz, tr_x_0_z_zzzz_xyy, tr_x_0_z_zzzz_xyz, tr_x_0_z_zzzz_xzz, tr_x_0_z_zzzz_yyy, tr_x_0_z_zzzz_yyz, tr_x_0_z_zzzz_yzz, tr_x_0_z_zzzz_zzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xxxz, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xz, tr_xzzzz_xzzz, tr_xzzzz_yy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yz, tr_xzzzz_yzzz, tr_xzzzz_zz, tr_xzzzz_zzzz, tr_xzzzzz_xxx, tr_xzzzzz_xxy, tr_xzzzzz_xxz, tr_xzzzzz_xyy, tr_xzzzzz_xyz, tr_xzzzzz_xzz, tr_xzzzzz_yyy, tr_xzzzzz_yyz, tr_xzzzzz_yzz, tr_xzzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zzzz_xxx[i] = -8.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_xxy[i] = -8.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_xxz[i] = -8.0 * tr_xzzz_xxz[i] * tbe_0 - 2.0 * tr_xzzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_xyy[i] = -8.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_xyz[i] = -8.0 * tr_xzzz_xyz[i] * tbe_0 - 2.0 * tr_xzzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_xzz[i] = -8.0 * tr_xzzz_xzz[i] * tbe_0 - 4.0 * tr_xzzzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_yyy[i] = -8.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_yyz[i] = -8.0 * tr_xzzz_yyz[i] * tbe_0 - 2.0 * tr_xzzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_yzz[i] = -8.0 * tr_xzzz_yzz[i] * tbe_0 - 4.0 * tr_xzzzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzzz_zzz[i] = -8.0 * tr_xzzz_zzz[i] * tbe_0 - 6.0 * tr_xzzzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-460 components of targeted buffer : GF

    auto tr_y_0_x_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 450);

    auto tr_y_0_x_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 451);

    auto tr_y_0_x_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 452);

    auto tr_y_0_x_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 453);

    auto tr_y_0_x_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 454);

    auto tr_y_0_x_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 455);

    auto tr_y_0_x_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 456);

    auto tr_y_0_x_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 457);

    auto tr_y_0_x_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 458);

    auto tr_y_0_x_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 459);

    #pragma omp simd aligned(tr_xxxxxy_xxx, tr_xxxxxy_xxy, tr_xxxxxy_xxz, tr_xxxxxy_xyy, tr_xxxxxy_xyz, tr_xxxxxy_xzz, tr_xxxxxy_yyy, tr_xxxxxy_yyz, tr_xxxxxy_yzz, tr_xxxxxy_zzz, tr_xxxxy_xx, tr_xxxxy_xxxx, tr_xxxxy_xxxy, tr_xxxxy_xxxz, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xy, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xz, tr_xxxxy_xzzz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_y_0_x_xxxx_xxx, tr_y_0_x_xxxx_xxy, tr_y_0_x_xxxx_xxz, tr_y_0_x_xxxx_xyy, tr_y_0_x_xxxx_xyz, tr_y_0_x_xxxx_xzz, tr_y_0_x_xxxx_yyy, tr_y_0_x_xxxx_yyz, tr_y_0_x_xxxx_yzz, tr_y_0_x_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxxx_xxx[i] = -8.0 * tr_xxxy_xxx[i] * tbe_0 - 6.0 * tr_xxxxy_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_xxy[i] = -8.0 * tr_xxxy_xxy[i] * tbe_0 - 4.0 * tr_xxxxy_xy[i] * tbe_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_xxz[i] = -8.0 * tr_xxxy_xxz[i] * tbe_0 - 4.0 * tr_xxxxy_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_xyy[i] = -8.0 * tr_xxxy_xyy[i] * tbe_0 - 2.0 * tr_xxxxy_yy[i] * tbe_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_xyz[i] = -8.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxxy_yz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_xzz[i] = -8.0 * tr_xxxy_xzz[i] * tbe_0 - 2.0 * tr_xxxxy_zz[i] * tbe_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_yyy[i] = -8.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_yyz[i] = -8.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_yzz[i] = -8.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxx_zzz[i] = -8.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 460-470 components of targeted buffer : GF

    auto tr_y_0_x_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 460);

    auto tr_y_0_x_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 461);

    auto tr_y_0_x_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 462);

    auto tr_y_0_x_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 463);

    auto tr_y_0_x_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 464);

    auto tr_y_0_x_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 465);

    auto tr_y_0_x_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 466);

    auto tr_y_0_x_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 467);

    auto tr_y_0_x_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 468);

    auto tr_y_0_x_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 469);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xx, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxyy_xxx, tr_xxxxyy_xxy, tr_xxxxyy_xxz, tr_xxxxyy_xyy, tr_xxxxyy_xyz, tr_xxxxyy_xzz, tr_xxxxyy_yyy, tr_xxxxyy_yyz, tr_xxxxyy_yzz, tr_xxxxyy_zzz, tr_xxxyy_xx, tr_xxxyy_xxxx, tr_xxxyy_xxxy, tr_xxxyy_xxxz, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xy, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xz, tr_xxxyy_xzzz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_y_0_x_xxxy_xxx, tr_y_0_x_xxxy_xxy, tr_y_0_x_xxxy_xxz, tr_y_0_x_xxxy_xyy, tr_y_0_x_xxxy_xyz, tr_y_0_x_xxxy_xzz, tr_y_0_x_xxxy_yyy, tr_y_0_x_xxxy_yyz, tr_y_0_x_xxxy_yzz, tr_y_0_x_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxxy_xxx[i] = 3.0 * tr_xx_xxx[i] - 6.0 * tr_xxyy_xxx[i] * tbe_0 + 3.0 * tr_xxx_xx[i] - 2.0 * tr_xxx_xxxx[i] * tke_0 - 6.0 * tr_xxxyy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_xxy[i] = 3.0 * tr_xx_xxy[i] - 6.0 * tr_xxyy_xxy[i] * tbe_0 + 2.0 * tr_xxx_xy[i] - 2.0 * tr_xxx_xxxy[i] * tke_0 - 4.0 * tr_xxxyy_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_xxz[i] = 3.0 * tr_xx_xxz[i] - 6.0 * tr_xxyy_xxz[i] * tbe_0 + 2.0 * tr_xxx_xz[i] - 2.0 * tr_xxx_xxxz[i] * tke_0 - 4.0 * tr_xxxyy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_xyy[i] = 3.0 * tr_xx_xyy[i] - 6.0 * tr_xxyy_xyy[i] * tbe_0 + tr_xxx_yy[i] - 2.0 * tr_xxx_xxyy[i] * tke_0 - 2.0 * tr_xxxyy_yy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_xyz[i] = 3.0 * tr_xx_xyz[i] - 6.0 * tr_xxyy_xyz[i] * tbe_0 + tr_xxx_yz[i] - 2.0 * tr_xxx_xxyz[i] * tke_0 - 2.0 * tr_xxxyy_yz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_xzz[i] = 3.0 * tr_xx_xzz[i] - 6.0 * tr_xxyy_xzz[i] * tbe_0 + tr_xxx_zz[i] - 2.0 * tr_xxx_xxzz[i] * tke_0 - 2.0 * tr_xxxyy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_yyy[i] = 3.0 * tr_xx_yyy[i] - 6.0 * tr_xxyy_yyy[i] * tbe_0 - 2.0 * tr_xxx_xyyy[i] * tke_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_yyz[i] = 3.0 * tr_xx_yyz[i] - 6.0 * tr_xxyy_yyz[i] * tbe_0 - 2.0 * tr_xxx_xyyz[i] * tke_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_yzz[i] = 3.0 * tr_xx_yzz[i] - 6.0 * tr_xxyy_yzz[i] * tbe_0 - 2.0 * tr_xxx_xyzz[i] * tke_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxy_zzz[i] = 3.0 * tr_xx_zzz[i] - 6.0 * tr_xxyy_zzz[i] * tbe_0 - 2.0 * tr_xxx_xzzz[i] * tke_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 470-480 components of targeted buffer : GF

    auto tr_y_0_x_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 470);

    auto tr_y_0_x_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 471);

    auto tr_y_0_x_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 472);

    auto tr_y_0_x_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 473);

    auto tr_y_0_x_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 474);

    auto tr_y_0_x_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 475);

    auto tr_y_0_x_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 476);

    auto tr_y_0_x_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 477);

    auto tr_y_0_x_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 478);

    auto tr_y_0_x_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 479);

    #pragma omp simd aligned(tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_xxxyz_xx, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xy, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_xzzz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_y_0_x_xxxz_xxx, tr_y_0_x_xxxz_xxy, tr_y_0_x_xxxz_xxz, tr_y_0_x_xxxz_xyy, tr_y_0_x_xxxz_xyz, tr_y_0_x_xxxz_xzz, tr_y_0_x_xxxz_yyy, tr_y_0_x_xxxz_yyz, tr_y_0_x_xxxz_yzz, tr_y_0_x_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxxz_xxx[i] = -6.0 * tr_xxyz_xxx[i] * tbe_0 - 6.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_xxy[i] = -6.0 * tr_xxyz_xxy[i] * tbe_0 - 4.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_xxz[i] = -6.0 * tr_xxyz_xxz[i] * tbe_0 - 4.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_xyy[i] = -6.0 * tr_xxyz_xyy[i] * tbe_0 - 2.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_xyz[i] = -6.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_xzz[i] = -6.0 * tr_xxyz_xzz[i] * tbe_0 - 2.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_yyy[i] = -6.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_yyz[i] = -6.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_yzz[i] = -6.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxxz_zzz[i] = -6.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-490 components of targeted buffer : GF

    auto tr_y_0_x_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 480);

    auto tr_y_0_x_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 481);

    auto tr_y_0_x_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 482);

    auto tr_y_0_x_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 483);

    auto tr_y_0_x_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 484);

    auto tr_y_0_x_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 485);

    auto tr_y_0_x_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 486);

    auto tr_y_0_x_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 487);

    auto tr_y_0_x_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 488);

    auto tr_y_0_x_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 489);

    #pragma omp simd aligned(tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyyy_xxx, tr_xxxyyy_xxy, tr_xxxyyy_xxz, tr_xxxyyy_xyy, tr_xxxyyy_xyz, tr_xxxyyy_xzz, tr_xxxyyy_yyy, tr_xxxyyy_yyz, tr_xxxyyy_yzz, tr_xxxyyy_zzz, tr_xxy_xx, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyyy_xx, tr_xxyyy_xxxx, tr_xxyyy_xxxy, tr_xxyyy_xxxz, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xy, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xz, tr_xxyyy_xzzz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_y_0_x_xxyy_xxx, tr_y_0_x_xxyy_xxy, tr_y_0_x_xxyy_xxz, tr_y_0_x_xxyy_xyy, tr_y_0_x_xxyy_xyz, tr_y_0_x_xxyy_xzz, tr_y_0_x_xxyy_yyy, tr_y_0_x_xxyy_yyz, tr_y_0_x_xxyy_yzz, tr_y_0_x_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxyy_xxx[i] = 4.0 * tr_xy_xxx[i] - 4.0 * tr_xyyy_xxx[i] * tbe_0 + 6.0 * tr_xxy_xx[i] - 4.0 * tr_xxy_xxxx[i] * tke_0 - 6.0 * tr_xxyyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_xxy[i] = 4.0 * tr_xy_xxy[i] - 4.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xxy_xy[i] - 4.0 * tr_xxy_xxxy[i] * tke_0 - 4.0 * tr_xxyyy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_xxz[i] = 4.0 * tr_xy_xxz[i] - 4.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xxy_xz[i] - 4.0 * tr_xxy_xxxz[i] * tke_0 - 4.0 * tr_xxyyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_xyy[i] = 4.0 * tr_xy_xyy[i] - 4.0 * tr_xyyy_xyy[i] * tbe_0 + 2.0 * tr_xxy_yy[i] - 4.0 * tr_xxy_xxyy[i] * tke_0 - 2.0 * tr_xxyyy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_xyz[i] = 4.0 * tr_xy_xyz[i] - 4.0 * tr_xyyy_xyz[i] * tbe_0 + 2.0 * tr_xxy_yz[i] - 4.0 * tr_xxy_xxyz[i] * tke_0 - 2.0 * tr_xxyyy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_xzz[i] = 4.0 * tr_xy_xzz[i] - 4.0 * tr_xyyy_xzz[i] * tbe_0 + 2.0 * tr_xxy_zz[i] - 4.0 * tr_xxy_xxzz[i] * tke_0 - 2.0 * tr_xxyyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_yyy[i] = 4.0 * tr_xy_yyy[i] - 4.0 * tr_xyyy_yyy[i] * tbe_0 - 4.0 * tr_xxy_xyyy[i] * tke_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_yyz[i] = 4.0 * tr_xy_yyz[i] - 4.0 * tr_xyyy_yyz[i] * tbe_0 - 4.0 * tr_xxy_xyyz[i] * tke_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_yzz[i] = 4.0 * tr_xy_yzz[i] - 4.0 * tr_xyyy_yzz[i] * tbe_0 - 4.0 * tr_xxy_xyzz[i] * tke_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyy_zzz[i] = 4.0 * tr_xy_zzz[i] - 4.0 * tr_xyyy_zzz[i] * tbe_0 - 4.0 * tr_xxy_xzzz[i] * tke_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 490-500 components of targeted buffer : GF

    auto tr_y_0_x_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 490);

    auto tr_y_0_x_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 491);

    auto tr_y_0_x_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 492);

    auto tr_y_0_x_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 493);

    auto tr_y_0_x_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 494);

    auto tr_y_0_x_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 495);

    auto tr_y_0_x_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 496);

    auto tr_y_0_x_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 497);

    auto tr_y_0_x_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 498);

    auto tr_y_0_x_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 499);

    #pragma omp simd aligned(tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxyyz_xx, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xy, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_xzzz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxz_xx, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_y_0_x_xxyz_xxx, tr_y_0_x_xxyz_xxy, tr_y_0_x_xxyz_xxz, tr_y_0_x_xxyz_xyy, tr_y_0_x_xxyz_xyz, tr_y_0_x_xxyz_xzz, tr_y_0_x_xxyz_yyy, tr_y_0_x_xxyz_yyz, tr_y_0_x_xxyz_yzz, tr_y_0_x_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxyz_xxx[i] = 2.0 * tr_xz_xxx[i] - 4.0 * tr_xyyz_xxx[i] * tbe_0 + 3.0 * tr_xxz_xx[i] - 2.0 * tr_xxz_xxxx[i] * tke_0 - 6.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_xxy[i] = 2.0 * tr_xz_xxy[i] - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 2.0 * tr_xxz_xy[i] - 2.0 * tr_xxz_xxxy[i] * tke_0 - 4.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_xxz[i] = 2.0 * tr_xz_xxz[i] - 4.0 * tr_xyyz_xxz[i] * tbe_0 + 2.0 * tr_xxz_xz[i] - 2.0 * tr_xxz_xxxz[i] * tke_0 - 4.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_xyy[i] = 2.0 * tr_xz_xyy[i] - 4.0 * tr_xyyz_xyy[i] * tbe_0 + tr_xxz_yy[i] - 2.0 * tr_xxz_xxyy[i] * tke_0 - 2.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_xyz[i] = 2.0 * tr_xz_xyz[i] - 4.0 * tr_xyyz_xyz[i] * tbe_0 + tr_xxz_yz[i] - 2.0 * tr_xxz_xxyz[i] * tke_0 - 2.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_xzz[i] = 2.0 * tr_xz_xzz[i] - 4.0 * tr_xyyz_xzz[i] * tbe_0 + tr_xxz_zz[i] - 2.0 * tr_xxz_xxzz[i] * tke_0 - 2.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_yyy[i] = 2.0 * tr_xz_yyy[i] - 4.0 * tr_xyyz_yyy[i] * tbe_0 - 2.0 * tr_xxz_xyyy[i] * tke_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_yyz[i] = 2.0 * tr_xz_yyz[i] - 4.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xxz_xyyz[i] * tke_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_yzz[i] = 2.0 * tr_xz_yzz[i] - 4.0 * tr_xyyz_yzz[i] * tbe_0 - 2.0 * tr_xxz_xyzz[i] * tke_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxyz_zzz[i] = 2.0 * tr_xz_zzz[i] - 4.0 * tr_xyyz_zzz[i] * tbe_0 - 2.0 * tr_xxz_xzzz[i] * tke_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 500-510 components of targeted buffer : GF

    auto tr_y_0_x_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 500);

    auto tr_y_0_x_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 501);

    auto tr_y_0_x_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 502);

    auto tr_y_0_x_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 503);

    auto tr_y_0_x_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 504);

    auto tr_y_0_x_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 505);

    auto tr_y_0_x_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 506);

    auto tr_y_0_x_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 507);

    auto tr_y_0_x_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 508);

    auto tr_y_0_x_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 509);

    #pragma omp simd aligned(tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xy, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_xzzz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_y_0_x_xxzz_xxx, tr_y_0_x_xxzz_xxy, tr_y_0_x_xxzz_xxz, tr_y_0_x_xxzz_xyy, tr_y_0_x_xxzz_xyz, tr_y_0_x_xxzz_xzz, tr_y_0_x_xxzz_yyy, tr_y_0_x_xxzz_yyz, tr_y_0_x_xxzz_yzz, tr_y_0_x_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxzz_xxx[i] = -4.0 * tr_xyzz_xxx[i] * tbe_0 - 6.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_xxy[i] = -4.0 * tr_xyzz_xxy[i] * tbe_0 - 4.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_xxz[i] = -4.0 * tr_xyzz_xxz[i] * tbe_0 - 4.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_xyy[i] = -4.0 * tr_xyzz_xyy[i] * tbe_0 - 2.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_xyz[i] = -4.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_xzz[i] = -4.0 * tr_xyzz_xzz[i] * tbe_0 - 2.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_yyy[i] = -4.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_yyz[i] = -4.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_yzz[i] = -4.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxzz_zzz[i] = -4.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-520 components of targeted buffer : GF

    auto tr_y_0_x_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 510);

    auto tr_y_0_x_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 511);

    auto tr_y_0_x_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 512);

    auto tr_y_0_x_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 513);

    auto tr_y_0_x_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 514);

    auto tr_y_0_x_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 515);

    auto tr_y_0_x_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 516);

    auto tr_y_0_x_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 517);

    auto tr_y_0_x_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 518);

    auto tr_y_0_x_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 519);

    #pragma omp simd aligned(tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyyy_xxx, tr_xxyyyy_xxy, tr_xxyyyy_xxz, tr_xxyyyy_xyy, tr_xxyyyy_xyz, tr_xxyyyy_xzz, tr_xxyyyy_yyy, tr_xxyyyy_yyz, tr_xxyyyy_yzz, tr_xxyyyy_zzz, tr_xyy_xx, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyyy_xx, tr_xyyyy_xxxx, tr_xyyyy_xxxy, tr_xyyyy_xxxz, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xy, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xz, tr_xyyyy_xzzz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, tr_y_0_x_xyyy_xxx, tr_y_0_x_xyyy_xxy, tr_y_0_x_xyyy_xxz, tr_y_0_x_xyyy_xyy, tr_y_0_x_xyyy_xyz, tr_y_0_x_xyyy_xzz, tr_y_0_x_xyyy_yyy, tr_y_0_x_xyyy_yyz, tr_y_0_x_xyyy_yzz, tr_y_0_x_xyyy_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyyy_xxx[i] = 3.0 * tr_yy_xxx[i] - 2.0 * tr_yyyy_xxx[i] * tbe_0 + 9.0 * tr_xyy_xx[i] - 6.0 * tr_xyy_xxxx[i] * tke_0 - 6.0 * tr_xyyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_xxy[i] = 3.0 * tr_yy_xxy[i] - 2.0 * tr_yyyy_xxy[i] * tbe_0 + 6.0 * tr_xyy_xy[i] - 6.0 * tr_xyy_xxxy[i] * tke_0 - 4.0 * tr_xyyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_xxz[i] = 3.0 * tr_yy_xxz[i] - 2.0 * tr_yyyy_xxz[i] * tbe_0 + 6.0 * tr_xyy_xz[i] - 6.0 * tr_xyy_xxxz[i] * tke_0 - 4.0 * tr_xyyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_xyy[i] = 3.0 * tr_yy_xyy[i] - 2.0 * tr_yyyy_xyy[i] * tbe_0 + 3.0 * tr_xyy_yy[i] - 6.0 * tr_xyy_xxyy[i] * tke_0 - 2.0 * tr_xyyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_xyz[i] = 3.0 * tr_yy_xyz[i] - 2.0 * tr_yyyy_xyz[i] * tbe_0 + 3.0 * tr_xyy_yz[i] - 6.0 * tr_xyy_xxyz[i] * tke_0 - 2.0 * tr_xyyyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_xzz[i] = 3.0 * tr_yy_xzz[i] - 2.0 * tr_yyyy_xzz[i] * tbe_0 + 3.0 * tr_xyy_zz[i] - 6.0 * tr_xyy_xxzz[i] * tke_0 - 2.0 * tr_xyyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_yyy[i] = 3.0 * tr_yy_yyy[i] - 2.0 * tr_yyyy_yyy[i] * tbe_0 - 6.0 * tr_xyy_xyyy[i] * tke_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_yyz[i] = 3.0 * tr_yy_yyz[i] - 2.0 * tr_yyyy_yyz[i] * tbe_0 - 6.0 * tr_xyy_xyyz[i] * tke_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_yzz[i] = 3.0 * tr_yy_yzz[i] - 2.0 * tr_yyyy_yzz[i] * tbe_0 - 6.0 * tr_xyy_xyzz[i] * tke_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyy_zzz[i] = 3.0 * tr_yy_zzz[i] - 2.0 * tr_yyyy_zzz[i] * tbe_0 - 6.0 * tr_xyy_xzzz[i] * tke_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 520-530 components of targeted buffer : GF

    auto tr_y_0_x_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 520);

    auto tr_y_0_x_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 521);

    auto tr_y_0_x_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 522);

    auto tr_y_0_x_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 523);

    auto tr_y_0_x_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 524);

    auto tr_y_0_x_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 525);

    auto tr_y_0_x_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 526);

    auto tr_y_0_x_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 527);

    auto tr_y_0_x_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 528);

    auto tr_y_0_x_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 529);

    #pragma omp simd aligned(tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyyyz_xx, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xy, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_xzzz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_y_0_x_xyyz_xxx, tr_y_0_x_xyyz_xxy, tr_y_0_x_xyyz_xxz, tr_y_0_x_xyyz_xyy, tr_y_0_x_xyyz_xyz, tr_y_0_x_xyyz_xzz, tr_y_0_x_xyyz_yyy, tr_y_0_x_xyyz_yyz, tr_y_0_x_xyyz_yzz, tr_y_0_x_xyyz_zzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyyz_xxx[i] = 2.0 * tr_yz_xxx[i] - 2.0 * tr_yyyz_xxx[i] * tbe_0 + 6.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxxx[i] * tke_0 - 6.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_xxy[i] = 2.0 * tr_yz_xxy[i] - 2.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xxxy[i] * tke_0 - 4.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_xxz[i] = 2.0 * tr_yz_xxz[i] - 2.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xxxz[i] * tke_0 - 4.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_xyy[i] = 2.0 * tr_yz_xyy[i] - 2.0 * tr_yyyz_xyy[i] * tbe_0 + 2.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_xxyy[i] * tke_0 - 2.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_xyz[i] = 2.0 * tr_yz_xyz[i] - 2.0 * tr_yyyz_xyz[i] * tbe_0 + 2.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 2.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_xzz[i] = 2.0 * tr_yz_xzz[i] - 2.0 * tr_yyyz_xzz[i] * tbe_0 + 2.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_xxzz[i] * tke_0 - 2.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_yyy[i] = 2.0 * tr_yz_yyy[i] - 2.0 * tr_yyyz_yyy[i] * tbe_0 - 4.0 * tr_xyz_xyyy[i] * tke_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_yyz[i] = 2.0 * tr_yz_yyz[i] - 2.0 * tr_yyyz_yyz[i] * tbe_0 - 4.0 * tr_xyz_xyyz[i] * tke_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_yzz[i] = 2.0 * tr_yz_yzz[i] - 2.0 * tr_yyyz_yzz[i] * tbe_0 - 4.0 * tr_xyz_xyzz[i] * tke_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyyz_zzz[i] = 2.0 * tr_yz_zzz[i] - 2.0 * tr_yyyz_zzz[i] * tbe_0 - 4.0 * tr_xyz_xzzz[i] * tke_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 530-540 components of targeted buffer : GF

    auto tr_y_0_x_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 530);

    auto tr_y_0_x_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 531);

    auto tr_y_0_x_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 532);

    auto tr_y_0_x_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 533);

    auto tr_y_0_x_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 534);

    auto tr_y_0_x_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 535);

    auto tr_y_0_x_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 536);

    auto tr_y_0_x_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 537);

    auto tr_y_0_x_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 538);

    auto tr_y_0_x_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 539);

    #pragma omp simd aligned(tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xy, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_xzzz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xzz_xx, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_y_0_x_xyzz_xxx, tr_y_0_x_xyzz_xxy, tr_y_0_x_xyzz_xxz, tr_y_0_x_xyzz_xyy, tr_y_0_x_xyzz_xyz, tr_y_0_x_xyzz_xzz, tr_y_0_x_xyzz_yyy, tr_y_0_x_xyzz_yyz, tr_y_0_x_xyzz_yzz, tr_y_0_x_xyzz_zzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyzz_xxx[i] = tr_zz_xxx[i] - 2.0 * tr_yyzz_xxx[i] * tbe_0 + 3.0 * tr_xzz_xx[i] - 2.0 * tr_xzz_xxxx[i] * tke_0 - 6.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_xxy[i] = tr_zz_xxy[i] - 2.0 * tr_yyzz_xxy[i] * tbe_0 + 2.0 * tr_xzz_xy[i] - 2.0 * tr_xzz_xxxy[i] * tke_0 - 4.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_xxz[i] = tr_zz_xxz[i] - 2.0 * tr_yyzz_xxz[i] * tbe_0 + 2.0 * tr_xzz_xz[i] - 2.0 * tr_xzz_xxxz[i] * tke_0 - 4.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_xyy[i] = tr_zz_xyy[i] - 2.0 * tr_yyzz_xyy[i] * tbe_0 + tr_xzz_yy[i] - 2.0 * tr_xzz_xxyy[i] * tke_0 - 2.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_xyz[i] = tr_zz_xyz[i] - 2.0 * tr_yyzz_xyz[i] * tbe_0 + tr_xzz_yz[i] - 2.0 * tr_xzz_xxyz[i] * tke_0 - 2.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_xzz[i] = tr_zz_xzz[i] - 2.0 * tr_yyzz_xzz[i] * tbe_0 + tr_xzz_zz[i] - 2.0 * tr_xzz_xxzz[i] * tke_0 - 2.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_yyy[i] = tr_zz_yyy[i] - 2.0 * tr_yyzz_yyy[i] * tbe_0 - 2.0 * tr_xzz_xyyy[i] * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_yyz[i] = tr_zz_yyz[i] - 2.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_xzz_xyyz[i] * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_yzz[i] = tr_zz_yzz[i] - 2.0 * tr_yyzz_yzz[i] * tbe_0 - 2.0 * tr_xzz_xyzz[i] * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyzz_zzz[i] = tr_zz_zzz[i] - 2.0 * tr_yyzz_zzz[i] * tbe_0 - 2.0 * tr_xzz_xzzz[i] * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 540-550 components of targeted buffer : GF

    auto tr_y_0_x_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 540);

    auto tr_y_0_x_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 541);

    auto tr_y_0_x_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 542);

    auto tr_y_0_x_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 543);

    auto tr_y_0_x_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 544);

    auto tr_y_0_x_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 545);

    auto tr_y_0_x_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 546);

    auto tr_y_0_x_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 547);

    auto tr_y_0_x_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 548);

    auto tr_y_0_x_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 549);

    #pragma omp simd aligned(tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xy, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_xzzz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_y_0_x_xzzz_xxx, tr_y_0_x_xzzz_xxy, tr_y_0_x_xzzz_xxz, tr_y_0_x_xzzz_xyy, tr_y_0_x_xzzz_xyz, tr_y_0_x_xzzz_xzz, tr_y_0_x_xzzz_yyy, tr_y_0_x_xzzz_yyz, tr_y_0_x_xzzz_yzz, tr_y_0_x_xzzz_zzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xzzz_xxx[i] = -2.0 * tr_yzzz_xxx[i] * tbe_0 - 6.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_xxy[i] = -2.0 * tr_yzzz_xxy[i] * tbe_0 - 4.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_xxz[i] = -2.0 * tr_yzzz_xxz[i] * tbe_0 - 4.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_xyy[i] = -2.0 * tr_yzzz_xyy[i] * tbe_0 - 2.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_xyz[i] = -2.0 * tr_yzzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_xzz[i] = -2.0 * tr_yzzz_xzz[i] * tbe_0 - 2.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_yyy[i] = -2.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_yyz[i] = -2.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_yzz[i] = -2.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzzz_zzz[i] = -2.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 550-560 components of targeted buffer : GF

    auto tr_y_0_x_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 550);

    auto tr_y_0_x_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 551);

    auto tr_y_0_x_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 552);

    auto tr_y_0_x_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 553);

    auto tr_y_0_x_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 554);

    auto tr_y_0_x_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 555);

    auto tr_y_0_x_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 556);

    auto tr_y_0_x_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 557);

    auto tr_y_0_x_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 558);

    auto tr_y_0_x_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 559);

    #pragma omp simd aligned(tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyyy_xxx, tr_xyyyyy_xxy, tr_xyyyyy_xxz, tr_xyyyyy_xyy, tr_xyyyyy_xyz, tr_xyyyyy_xzz, tr_xyyyyy_yyy, tr_xyyyyy_yyz, tr_xyyyyy_yzz, tr_xyyyyy_zzz, tr_y_0_x_yyyy_xxx, tr_y_0_x_yyyy_xxy, tr_y_0_x_yyyy_xxz, tr_y_0_x_yyyy_xyy, tr_y_0_x_yyyy_xyz, tr_y_0_x_yyyy_xzz, tr_y_0_x_yyyy_yyy, tr_y_0_x_yyyy_yyz, tr_y_0_x_yyyy_yzz, tr_y_0_x_yyyy_zzz, tr_yyy_xx, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyyy_xx, tr_yyyyy_xxxx, tr_yyyyy_xxxy, tr_yyyyy_xxxz, tr_yyyyy_xxyy, tr_yyyyy_xxyz, tr_yyyyy_xxzz, tr_yyyyy_xy, tr_yyyyy_xyyy, tr_yyyyy_xyyz, tr_yyyyy_xyzz, tr_yyyyy_xz, tr_yyyyy_xzzz, tr_yyyyy_yy, tr_yyyyy_yz, tr_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyyy_xxx[i] = 12.0 * tr_yyy_xx[i] - 8.0 * tr_yyy_xxxx[i] * tke_0 - 6.0 * tr_yyyyy_xx[i] * tbe_0 + 4.0 * tr_yyyyy_xxxx[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_xxy[i] = 8.0 * tr_yyy_xy[i] - 8.0 * tr_yyy_xxxy[i] * tke_0 - 4.0 * tr_yyyyy_xy[i] * tbe_0 + 4.0 * tr_yyyyy_xxxy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_xxz[i] = 8.0 * tr_yyy_xz[i] - 8.0 * tr_yyy_xxxz[i] * tke_0 - 4.0 * tr_yyyyy_xz[i] * tbe_0 + 4.0 * tr_yyyyy_xxxz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_xyy[i] = 4.0 * tr_yyy_yy[i] - 8.0 * tr_yyy_xxyy[i] * tke_0 - 2.0 * tr_yyyyy_yy[i] * tbe_0 + 4.0 * tr_yyyyy_xxyy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_xyz[i] = 4.0 * tr_yyy_yz[i] - 8.0 * tr_yyy_xxyz[i] * tke_0 - 2.0 * tr_yyyyy_yz[i] * tbe_0 + 4.0 * tr_yyyyy_xxyz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_xzz[i] = 4.0 * tr_yyy_zz[i] - 8.0 * tr_yyy_xxzz[i] * tke_0 - 2.0 * tr_yyyyy_zz[i] * tbe_0 + 4.0 * tr_yyyyy_xxzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_yyy[i] = -8.0 * tr_yyy_xyyy[i] * tke_0 + 4.0 * tr_yyyyy_xyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_yyz[i] = -8.0 * tr_yyy_xyyz[i] * tke_0 + 4.0 * tr_yyyyy_xyyz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_yzz[i] = -8.0 * tr_yyy_xyzz[i] * tke_0 + 4.0 * tr_yyyyy_xyzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyy_zzz[i] = -8.0 * tr_yyy_xzzz[i] * tke_0 + 4.0 * tr_yyyyy_xzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 560-570 components of targeted buffer : GF

    auto tr_y_0_x_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 560);

    auto tr_y_0_x_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 561);

    auto tr_y_0_x_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 562);

    auto tr_y_0_x_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 563);

    auto tr_y_0_x_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 564);

    auto tr_y_0_x_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 565);

    auto tr_y_0_x_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 566);

    auto tr_y_0_x_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 567);

    auto tr_y_0_x_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 568);

    auto tr_y_0_x_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 569);

    #pragma omp simd aligned(tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_y_0_x_yyyz_xxx, tr_y_0_x_yyyz_xxy, tr_y_0_x_yyyz_xxz, tr_y_0_x_yyyz_xyy, tr_y_0_x_yyyz_xyz, tr_y_0_x_yyyz_xzz, tr_y_0_x_yyyz_yyy, tr_y_0_x_yyyz_yyz, tr_y_0_x_yyyz_yzz, tr_y_0_x_yyyz_zzz, tr_yyyyz_xx, tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xy, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xz, tr_yyyyz_xzzz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyz_xx, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyyz_xxx[i] = 9.0 * tr_yyz_xx[i] - 6.0 * tr_yyz_xxxx[i] * tke_0 - 6.0 * tr_yyyyz_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_xxy[i] = 6.0 * tr_yyz_xy[i] - 6.0 * tr_yyz_xxxy[i] * tke_0 - 4.0 * tr_yyyyz_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_xxz[i] = 6.0 * tr_yyz_xz[i] - 6.0 * tr_yyz_xxxz[i] * tke_0 - 4.0 * tr_yyyyz_xz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_xyy[i] = 3.0 * tr_yyz_yy[i] - 6.0 * tr_yyz_xxyy[i] * tke_0 - 2.0 * tr_yyyyz_yy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_xyz[i] = 3.0 * tr_yyz_yz[i] - 6.0 * tr_yyz_xxyz[i] * tke_0 - 2.0 * tr_yyyyz_yz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_xzz[i] = 3.0 * tr_yyz_zz[i] - 6.0 * tr_yyz_xxzz[i] * tke_0 - 2.0 * tr_yyyyz_zz[i] * tbe_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_yyy[i] = -6.0 * tr_yyz_xyyy[i] * tke_0 + 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_yyz[i] = -6.0 * tr_yyz_xyyz[i] * tke_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_yzz[i] = -6.0 * tr_yyz_xyzz[i] * tke_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyyz_zzz[i] = -6.0 * tr_yyz_xzzz[i] * tke_0 + 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 570-580 components of targeted buffer : GF

    auto tr_y_0_x_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 570);

    auto tr_y_0_x_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 571);

    auto tr_y_0_x_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 572);

    auto tr_y_0_x_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 573);

    auto tr_y_0_x_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 574);

    auto tr_y_0_x_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 575);

    auto tr_y_0_x_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 576);

    auto tr_y_0_x_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 577);

    auto tr_y_0_x_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 578);

    auto tr_y_0_x_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 579);

    #pragma omp simd aligned(tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_y_0_x_yyzz_xxx, tr_y_0_x_yyzz_xxy, tr_y_0_x_yyzz_xxz, tr_y_0_x_yyzz_xyy, tr_y_0_x_yyzz_xyz, tr_y_0_x_yyzz_xzz, tr_y_0_x_yyzz_yyy, tr_y_0_x_yyzz_yyz, tr_y_0_x_yyzz_yzz, tr_y_0_x_yyzz_zzz, tr_yyyzz_xx, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xy, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xz, tr_yyyzz_xzzz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_yzz_xx, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyzz_xxx[i] = 6.0 * tr_yzz_xx[i] - 4.0 * tr_yzz_xxxx[i] * tke_0 - 6.0 * tr_yyyzz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_xxy[i] = 4.0 * tr_yzz_xy[i] - 4.0 * tr_yzz_xxxy[i] * tke_0 - 4.0 * tr_yyyzz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_xxz[i] = 4.0 * tr_yzz_xz[i] - 4.0 * tr_yzz_xxxz[i] * tke_0 - 4.0 * tr_yyyzz_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_xyy[i] = 2.0 * tr_yzz_yy[i] - 4.0 * tr_yzz_xxyy[i] * tke_0 - 2.0 * tr_yyyzz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_xyz[i] = 2.0 * tr_yzz_yz[i] - 4.0 * tr_yzz_xxyz[i] * tke_0 - 2.0 * tr_yyyzz_yz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_xzz[i] = 2.0 * tr_yzz_zz[i] - 4.0 * tr_yzz_xxzz[i] * tke_0 - 2.0 * tr_yyyzz_zz[i] * tbe_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_yyy[i] = -4.0 * tr_yzz_xyyy[i] * tke_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_yyz[i] = -4.0 * tr_yzz_xyyz[i] * tke_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_yzz[i] = -4.0 * tr_yzz_xyzz[i] * tke_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyzz_zzz[i] = -4.0 * tr_yzz_xzzz[i] * tke_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 580-590 components of targeted buffer : GF

    auto tr_y_0_x_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 580);

    auto tr_y_0_x_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 581);

    auto tr_y_0_x_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 582);

    auto tr_y_0_x_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 583);

    auto tr_y_0_x_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 584);

    auto tr_y_0_x_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 585);

    auto tr_y_0_x_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 586);

    auto tr_y_0_x_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 587);

    auto tr_y_0_x_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 588);

    auto tr_y_0_x_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 589);

    #pragma omp simd aligned(tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_y_0_x_yzzz_xxx, tr_y_0_x_yzzz_xxy, tr_y_0_x_yzzz_xxz, tr_y_0_x_yzzz_xyy, tr_y_0_x_yzzz_xyz, tr_y_0_x_yzzz_xzz, tr_y_0_x_yzzz_yyy, tr_y_0_x_yzzz_yyz, tr_y_0_x_yzzz_yzz, tr_y_0_x_yzzz_zzz, tr_yyzzz_xx, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xy, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xz, tr_yyzzz_xzzz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_zzz_xx, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yzzz_xxx[i] = 3.0 * tr_zzz_xx[i] - 2.0 * tr_zzz_xxxx[i] * tke_0 - 6.0 * tr_yyzzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_xxy[i] = 2.0 * tr_zzz_xy[i] - 2.0 * tr_zzz_xxxy[i] * tke_0 - 4.0 * tr_yyzzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_xxz[i] = 2.0 * tr_zzz_xz[i] - 2.0 * tr_zzz_xxxz[i] * tke_0 - 4.0 * tr_yyzzz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_xyy[i] = tr_zzz_yy[i] - 2.0 * tr_zzz_xxyy[i] * tke_0 - 2.0 * tr_yyzzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_xyz[i] = tr_zzz_yz[i] - 2.0 * tr_zzz_xxyz[i] * tke_0 - 2.0 * tr_yyzzz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_xzz[i] = tr_zzz_zz[i] - 2.0 * tr_zzz_xxzz[i] * tke_0 - 2.0 * tr_yyzzz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_yyy[i] = -2.0 * tr_zzz_xyyy[i] * tke_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_yyz[i] = -2.0 * tr_zzz_xyyz[i] * tke_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_yzz[i] = -2.0 * tr_zzz_xyzz[i] * tke_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzzz_zzz[i] = -2.0 * tr_zzz_xzzz[i] * tke_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 590-600 components of targeted buffer : GF

    auto tr_y_0_x_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 590);

    auto tr_y_0_x_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 591);

    auto tr_y_0_x_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 592);

    auto tr_y_0_x_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 593);

    auto tr_y_0_x_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 594);

    auto tr_y_0_x_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 595);

    auto tr_y_0_x_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 596);

    auto tr_y_0_x_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 597);

    auto tr_y_0_x_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 598);

    auto tr_y_0_x_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 599);

    #pragma omp simd aligned(tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, tr_y_0_x_zzzz_xxx, tr_y_0_x_zzzz_xxy, tr_y_0_x_zzzz_xxz, tr_y_0_x_zzzz_xyy, tr_y_0_x_zzzz_xyz, tr_y_0_x_zzzz_xzz, tr_y_0_x_zzzz_yyy, tr_y_0_x_zzzz_yyz, tr_y_0_x_zzzz_yzz, tr_y_0_x_zzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xy, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xz, tr_yzzzz_xzzz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zzzz_xxx[i] = -6.0 * tr_yzzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_xxy[i] = -4.0 * tr_yzzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_xxz[i] = -4.0 * tr_yzzzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_xyy[i] = -2.0 * tr_yzzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_xyz[i] = -2.0 * tr_yzzzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_xzz[i] = -2.0 * tr_yzzzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_yyy[i] = 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_yyz[i] = 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_yzz[i] = 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzzz_zzz[i] = 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 600-610 components of targeted buffer : GF

    auto tr_y_0_y_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 600);

    auto tr_y_0_y_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 601);

    auto tr_y_0_y_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 602);

    auto tr_y_0_y_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 603);

    auto tr_y_0_y_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 604);

    auto tr_y_0_y_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 605);

    auto tr_y_0_y_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 606);

    auto tr_y_0_y_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 607);

    auto tr_y_0_y_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 608);

    auto tr_y_0_y_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 609);

    #pragma omp simd aligned(tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxy_xx, tr_xxxxy_xxxy, tr_xxxxy_xxyy, tr_xxxxy_xxyz, tr_xxxxy_xy, tr_xxxxy_xyyy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yyyy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yz, tr_xxxxy_yzzz, tr_xxxxy_zz, tr_xxxxyy_xxx, tr_xxxxyy_xxy, tr_xxxxyy_xxz, tr_xxxxyy_xyy, tr_xxxxyy_xyz, tr_xxxxyy_xzz, tr_xxxxyy_yyy, tr_xxxxyy_yyz, tr_xxxxyy_yzz, tr_xxxxyy_zzz, tr_y_0_y_xxxx_xxx, tr_y_0_y_xxxx_xxy, tr_y_0_y_xxxx_xxz, tr_y_0_y_xxxx_xyy, tr_y_0_y_xxxx_xyz, tr_y_0_y_xxxx_xzz, tr_y_0_y_xxxx_yyy, tr_y_0_y_xxxx_yyz, tr_y_0_y_xxxx_yzz, tr_y_0_y_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxxx_xxx[i] = -2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_xxy[i] = -2.0 * tr_xxxx_xxy[i] * tbe_0 - 2.0 * tr_xxxxy_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_xxz[i] = -2.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_xyy[i] = -2.0 * tr_xxxx_xyy[i] * tbe_0 - 4.0 * tr_xxxxy_xy[i] * tbe_0 + 4.0 * tr_xxxxy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_xyz[i] = -2.0 * tr_xxxx_xyz[i] * tbe_0 - 2.0 * tr_xxxxy_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_xzz[i] = -2.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_yyy[i] = -2.0 * tr_xxxx_yyy[i] * tbe_0 - 6.0 * tr_xxxxy_yy[i] * tbe_0 + 4.0 * tr_xxxxy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_yyz[i] = -2.0 * tr_xxxx_yyz[i] * tbe_0 - 4.0 * tr_xxxxy_yz[i] * tbe_0 + 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_yzz[i] = -2.0 * tr_xxxx_yzz[i] * tbe_0 - 2.0 * tr_xxxxy_zz[i] * tbe_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxx_zzz[i] = -2.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 610-620 components of targeted buffer : GF

    auto tr_y_0_y_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 610);

    auto tr_y_0_y_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 611);

    auto tr_y_0_y_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 612);

    auto tr_y_0_y_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 613);

    auto tr_y_0_y_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 614);

    auto tr_y_0_y_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 615);

    auto tr_y_0_y_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 616);

    auto tr_y_0_y_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 617);

    auto tr_y_0_y_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 618);

    auto tr_y_0_y_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 619);

    #pragma omp simd aligned(tr_xxx_xx, tr_xxx_xxxy, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_yy, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyy_xx, tr_xxxyy_xxxy, tr_xxxyy_xxyy, tr_xxxyy_xxyz, tr_xxxyy_xy, tr_xxxyy_xyyy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yyyy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yz, tr_xxxyy_yzzz, tr_xxxyy_zz, tr_xxxyyy_xxx, tr_xxxyyy_xxy, tr_xxxyyy_xxz, tr_xxxyyy_xyy, tr_xxxyyy_xyz, tr_xxxyyy_xzz, tr_xxxyyy_yyy, tr_xxxyyy_yyz, tr_xxxyyy_yzz, tr_xxxyyy_zzz, tr_y_0_y_xxxy_xxx, tr_y_0_y_xxxy_xxy, tr_y_0_y_xxxy_xxz, tr_y_0_y_xxxy_xyy, tr_y_0_y_xxxy_xyz, tr_y_0_y_xxxy_xzz, tr_y_0_y_xxxy_yyy, tr_y_0_y_xxxy_yyz, tr_y_0_y_xxxy_yzz, tr_y_0_y_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxxy_xxx[i] = -2.0 * tr_xxx_xxxy[i] * tke_0 - 6.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_xxy[i] = tr_xxx_xx[i] - 2.0 * tr_xxx_xxyy[i] * tke_0 - 6.0 * tr_xxxy_xxy[i] * tbe_0 - 2.0 * tr_xxxyy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_xxz[i] = -2.0 * tr_xxx_xxyz[i] * tke_0 - 6.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_xyy[i] = 2.0 * tr_xxx_xy[i] - 2.0 * tr_xxx_xyyy[i] * tke_0 - 6.0 * tr_xxxy_xyy[i] * tbe_0 - 4.0 * tr_xxxyy_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_xyz[i] = tr_xxx_xz[i] - 2.0 * tr_xxx_xyyz[i] * tke_0 - 6.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxyy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_xzz[i] = -2.0 * tr_xxx_xyzz[i] * tke_0 - 6.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_yyy[i] = 3.0 * tr_xxx_yy[i] - 2.0 * tr_xxx_yyyy[i] * tke_0 - 6.0 * tr_xxxy_yyy[i] * tbe_0 - 6.0 * tr_xxxyy_yy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_yyz[i] = 2.0 * tr_xxx_yz[i] - 2.0 * tr_xxx_yyyz[i] * tke_0 - 6.0 * tr_xxxy_yyz[i] * tbe_0 - 4.0 * tr_xxxyy_yz[i] * tbe_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_yzz[i] = tr_xxx_zz[i] - 2.0 * tr_xxx_yyzz[i] * tke_0 - 6.0 * tr_xxxy_yzz[i] * tbe_0 - 2.0 * tr_xxxyy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxy_zzz[i] = -2.0 * tr_xxx_yzzz[i] * tke_0 - 6.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 620-630 components of targeted buffer : GF

    auto tr_y_0_y_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 620);

    auto tr_y_0_y_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 621);

    auto tr_y_0_y_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 622);

    auto tr_y_0_y_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 623);

    auto tr_y_0_y_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 624);

    auto tr_y_0_y_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 625);

    auto tr_y_0_y_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 626);

    auto tr_y_0_y_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 627);

    auto tr_y_0_y_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 628);

    auto tr_y_0_y_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 629);

    #pragma omp simd aligned(tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxxyz_xx, tr_xxxyz_xxxy, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xy, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yz, tr_xxxyz_yzzz, tr_xxxyz_zz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_y_0_y_xxxz_xxx, tr_y_0_y_xxxz_xxy, tr_y_0_y_xxxz_xxz, tr_y_0_y_xxxz_xyy, tr_y_0_y_xxxz_xyz, tr_y_0_y_xxxz_xzz, tr_y_0_y_xxxz_yyy, tr_y_0_y_xxxz_yyz, tr_y_0_y_xxxz_yzz, tr_y_0_y_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxxz_xxx[i] = -2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_xxy[i] = -2.0 * tr_xxxz_xxy[i] * tbe_0 - 2.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_xxz[i] = -2.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_xyy[i] = -2.0 * tr_xxxz_xyy[i] * tbe_0 - 4.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_xyz[i] = -2.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_xzz[i] = -2.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_yyy[i] = -2.0 * tr_xxxz_yyy[i] * tbe_0 - 6.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_yyz[i] = -2.0 * tr_xxxz_yyz[i] * tbe_0 - 4.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_yzz[i] = -2.0 * tr_xxxz_yzz[i] * tbe_0 - 2.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxxz_zzz[i] = -2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 630-640 components of targeted buffer : GF

    auto tr_y_0_y_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 630);

    auto tr_y_0_y_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 631);

    auto tr_y_0_y_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 632);

    auto tr_y_0_y_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 633);

    auto tr_y_0_y_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 634);

    auto tr_y_0_y_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 635);

    auto tr_y_0_y_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 636);

    auto tr_y_0_y_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 637);

    auto tr_y_0_y_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 638);

    auto tr_y_0_y_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 639);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxy_xx, tr_xxy_xxxy, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_yy, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyy_xx, tr_xxyyy_xxxy, tr_xxyyy_xxyy, tr_xxyyy_xxyz, tr_xxyyy_xy, tr_xxyyy_xyyy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yyyy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yz, tr_xxyyy_yzzz, tr_xxyyy_zz, tr_xxyyyy_xxx, tr_xxyyyy_xxy, tr_xxyyyy_xxz, tr_xxyyyy_xyy, tr_xxyyyy_xyz, tr_xxyyyy_xzz, tr_xxyyyy_yyy, tr_xxyyyy_yyz, tr_xxyyyy_yzz, tr_xxyyyy_zzz, tr_y_0_y_xxyy_xxx, tr_y_0_y_xxyy_xxy, tr_y_0_y_xxyy_xxz, tr_y_0_y_xxyy_xyy, tr_y_0_y_xxyy_xyz, tr_y_0_y_xxyy_xzz, tr_y_0_y_xxyy_yyy, tr_y_0_y_xxyy_yyz, tr_y_0_y_xxyy_yzz, tr_y_0_y_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxyy_xxx[i] = 2.0 * tr_xx_xxx[i] - 4.0 * tr_xxy_xxxy[i] * tke_0 - 10.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_xxy[i] = 2.0 * tr_xx_xxy[i] + 2.0 * tr_xxy_xx[i] - 4.0 * tr_xxy_xxyy[i] * tke_0 - 10.0 * tr_xxyy_xxy[i] * tbe_0 - 2.0 * tr_xxyyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_xxz[i] = 2.0 * tr_xx_xxz[i] - 4.0 * tr_xxy_xxyz[i] * tke_0 - 10.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_xyy[i] = 2.0 * tr_xx_xyy[i] + 4.0 * tr_xxy_xy[i] - 4.0 * tr_xxy_xyyy[i] * tke_0 - 10.0 * tr_xxyy_xyy[i] * tbe_0 - 4.0 * tr_xxyyy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_xyz[i] = 2.0 * tr_xx_xyz[i] + 2.0 * tr_xxy_xz[i] - 4.0 * tr_xxy_xyyz[i] * tke_0 - 10.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxyyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_xzz[i] = 2.0 * tr_xx_xzz[i] - 4.0 * tr_xxy_xyzz[i] * tke_0 - 10.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_yyy[i] = 2.0 * tr_xx_yyy[i] + 6.0 * tr_xxy_yy[i] - 4.0 * tr_xxy_yyyy[i] * tke_0 - 10.0 * tr_xxyy_yyy[i] * tbe_0 - 6.0 * tr_xxyyy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_yyz[i] = 2.0 * tr_xx_yyz[i] + 4.0 * tr_xxy_yz[i] - 4.0 * tr_xxy_yyyz[i] * tke_0 - 10.0 * tr_xxyy_yyz[i] * tbe_0 - 4.0 * tr_xxyyy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_yzz[i] = 2.0 * tr_xx_yzz[i] + 2.0 * tr_xxy_zz[i] - 4.0 * tr_xxy_yyzz[i] * tke_0 - 10.0 * tr_xxyy_yzz[i] * tbe_0 - 2.0 * tr_xxyyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyy_zzz[i] = 2.0 * tr_xx_zzz[i] - 4.0 * tr_xxy_yzzz[i] * tke_0 - 10.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 640-650 components of targeted buffer : GF

    auto tr_y_0_y_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 640);

    auto tr_y_0_y_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 641);

    auto tr_y_0_y_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 642);

    auto tr_y_0_y_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 643);

    auto tr_y_0_y_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 644);

    auto tr_y_0_y_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 645);

    auto tr_y_0_y_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 646);

    auto tr_y_0_y_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 647);

    auto tr_y_0_y_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 648);

    auto tr_y_0_y_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 649);

    #pragma omp simd aligned(tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xxyyz_xx, tr_xxyyz_xxxy, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xy, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yz, tr_xxyyz_yzzz, tr_xxyyz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxz_xx, tr_xxz_xxxy, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_yy, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_y_0_y_xxyz_xxx, tr_y_0_y_xxyz_xxy, tr_y_0_y_xxyz_xxz, tr_y_0_y_xxyz_xyy, tr_y_0_y_xxyz_xyz, tr_y_0_y_xxyz_xzz, tr_y_0_y_xxyz_yyy, tr_y_0_y_xxyz_yyz, tr_y_0_y_xxyz_yzz, tr_y_0_y_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxyz_xxx[i] = -2.0 * tr_xxz_xxxy[i] * tke_0 - 6.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_xxy[i] = tr_xxz_xx[i] - 2.0 * tr_xxz_xxyy[i] * tke_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 - 2.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_xxz[i] = -2.0 * tr_xxz_xxyz[i] * tke_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_xyy[i] = 2.0 * tr_xxz_xy[i] - 2.0 * tr_xxz_xyyy[i] * tke_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 - 4.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_xyz[i] = tr_xxz_xz[i] - 2.0 * tr_xxz_xyyz[i] * tke_0 - 6.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_xzz[i] = -2.0 * tr_xxz_xyzz[i] * tke_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_yyy[i] = 3.0 * tr_xxz_yy[i] - 2.0 * tr_xxz_yyyy[i] * tke_0 - 6.0 * tr_xxyz_yyy[i] * tbe_0 - 6.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_yyz[i] = 2.0 * tr_xxz_yz[i] - 2.0 * tr_xxz_yyyz[i] * tke_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 - 4.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_yzz[i] = tr_xxz_zz[i] - 2.0 * tr_xxz_yyzz[i] * tke_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 - 2.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxyz_zzz[i] = -2.0 * tr_xxz_yzzz[i] * tke_0 - 6.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 650-660 components of targeted buffer : GF

    auto tr_y_0_y_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 650);

    auto tr_y_0_y_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 651);

    auto tr_y_0_y_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 652);

    auto tr_y_0_y_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 653);

    auto tr_y_0_y_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 654);

    auto tr_y_0_y_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 655);

    auto tr_y_0_y_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 656);

    auto tr_y_0_y_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 657);

    auto tr_y_0_y_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 658);

    auto tr_y_0_y_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 659);

    #pragma omp simd aligned(tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxy, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xy, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yz, tr_xxyzz_yzzz, tr_xxyzz_zz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_y_0_y_xxzz_xxx, tr_y_0_y_xxzz_xxy, tr_y_0_y_xxzz_xxz, tr_y_0_y_xxzz_xyy, tr_y_0_y_xxzz_xyz, tr_y_0_y_xxzz_xzz, tr_y_0_y_xxzz_yyy, tr_y_0_y_xxzz_yyz, tr_y_0_y_xxzz_yzz, tr_y_0_y_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxzz_xxx[i] = -2.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_xxy[i] = -2.0 * tr_xxzz_xxy[i] * tbe_0 - 2.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_xxz[i] = -2.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_xyy[i] = -2.0 * tr_xxzz_xyy[i] * tbe_0 - 4.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_xyz[i] = -2.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_xzz[i] = -2.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_yyy[i] = -2.0 * tr_xxzz_yyy[i] * tbe_0 - 6.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_yyz[i] = -2.0 * tr_xxzz_yyz[i] * tbe_0 - 4.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_yzz[i] = -2.0 * tr_xxzz_yzz[i] * tbe_0 - 2.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxzz_zzz[i] = -2.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 660-670 components of targeted buffer : GF

    auto tr_y_0_y_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 660);

    auto tr_y_0_y_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 661);

    auto tr_y_0_y_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 662);

    auto tr_y_0_y_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 663);

    auto tr_y_0_y_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 664);

    auto tr_y_0_y_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 665);

    auto tr_y_0_y_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 666);

    auto tr_y_0_y_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 667);

    auto tr_y_0_y_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 668);

    auto tr_y_0_y_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 669);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxy, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_yy, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyy_xx, tr_xyyyy_xxxy, tr_xyyyy_xxyy, tr_xyyyy_xxyz, tr_xyyyy_xy, tr_xyyyy_xyyy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yyyy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yz, tr_xyyyy_yzzz, tr_xyyyy_zz, tr_xyyyyy_xxx, tr_xyyyyy_xxy, tr_xyyyyy_xxz, tr_xyyyyy_xyy, tr_xyyyyy_xyz, tr_xyyyyy_xzz, tr_xyyyyy_yyy, tr_xyyyyy_yyz, tr_xyyyyy_yzz, tr_xyyyyy_zzz, tr_y_0_y_xyyy_xxx, tr_y_0_y_xyyy_xxy, tr_y_0_y_xyyy_xxz, tr_y_0_y_xyyy_xyy, tr_y_0_y_xyyy_xyz, tr_y_0_y_xyyy_xzz, tr_y_0_y_xyyy_yyy, tr_y_0_y_xyyy_yyz, tr_y_0_y_xyyy_yzz, tr_y_0_y_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyyy_xxx[i] = 6.0 * tr_xy_xxx[i] - 6.0 * tr_xyy_xxxy[i] * tke_0 - 14.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_xxy[i] = 6.0 * tr_xy_xxy[i] + 3.0 * tr_xyy_xx[i] - 6.0 * tr_xyy_xxyy[i] * tke_0 - 14.0 * tr_xyyy_xxy[i] * tbe_0 - 2.0 * tr_xyyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_xxz[i] = 6.0 * tr_xy_xxz[i] - 6.0 * tr_xyy_xxyz[i] * tke_0 - 14.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_xyy[i] = 6.0 * tr_xy_xyy[i] + 6.0 * tr_xyy_xy[i] - 6.0 * tr_xyy_xyyy[i] * tke_0 - 14.0 * tr_xyyy_xyy[i] * tbe_0 - 4.0 * tr_xyyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_xyz[i] = 6.0 * tr_xy_xyz[i] + 3.0 * tr_xyy_xz[i] - 6.0 * tr_xyy_xyyz[i] * tke_0 - 14.0 * tr_xyyy_xyz[i] * tbe_0 - 2.0 * tr_xyyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_xzz[i] = 6.0 * tr_xy_xzz[i] - 6.0 * tr_xyy_xyzz[i] * tke_0 - 14.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_yyy[i] = 6.0 * tr_xy_yyy[i] + 9.0 * tr_xyy_yy[i] - 6.0 * tr_xyy_yyyy[i] * tke_0 - 14.0 * tr_xyyy_yyy[i] * tbe_0 - 6.0 * tr_xyyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_yyz[i] = 6.0 * tr_xy_yyz[i] + 6.0 * tr_xyy_yz[i] - 6.0 * tr_xyy_yyyz[i] * tke_0 - 14.0 * tr_xyyy_yyz[i] * tbe_0 - 4.0 * tr_xyyyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_yzz[i] = 6.0 * tr_xy_yzz[i] + 3.0 * tr_xyy_zz[i] - 6.0 * tr_xyy_yyzz[i] * tke_0 - 14.0 * tr_xyyy_yzz[i] * tbe_0 - 2.0 * tr_xyyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyy_zzz[i] = 6.0 * tr_xy_zzz[i] - 6.0 * tr_xyy_yzzz[i] * tke_0 - 14.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 670-680 components of targeted buffer : GF

    auto tr_y_0_y_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 670);

    auto tr_y_0_y_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 671);

    auto tr_y_0_y_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 672);

    auto tr_y_0_y_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 673);

    auto tr_y_0_y_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 674);

    auto tr_y_0_y_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 675);

    auto tr_y_0_y_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 676);

    auto tr_y_0_y_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 677);

    auto tr_y_0_y_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 678);

    auto tr_y_0_y_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 679);

    #pragma omp simd aligned(tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, tr_xyyyz_xx, tr_xyyyz_xxxy, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xy, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yz, tr_xyyyz_yzzz, tr_xyyyz_zz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_y_0_y_xyyz_xxx, tr_y_0_y_xyyz_xxy, tr_y_0_y_xyyz_xxz, tr_y_0_y_xyyz_xyy, tr_y_0_y_xyyz_xyz, tr_y_0_y_xyyz_xzz, tr_y_0_y_xyyz_yyy, tr_y_0_y_xyyz_yyz, tr_y_0_y_xyyz_yzz, tr_y_0_y_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyyz_xxx[i] = 2.0 * tr_xz_xxx[i] - 4.0 * tr_xyz_xxxy[i] * tke_0 - 10.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_xxy[i] = 2.0 * tr_xz_xxy[i] + 2.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxyy[i] * tke_0 - 10.0 * tr_xyyz_xxy[i] * tbe_0 - 2.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_xxz[i] = 2.0 * tr_xz_xxz[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 10.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_xyy[i] = 2.0 * tr_xz_xyy[i] + 4.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xyyy[i] * tke_0 - 10.0 * tr_xyyz_xyy[i] * tbe_0 - 4.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_xyz[i] = 2.0 * tr_xz_xyz[i] + 2.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 10.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_xzz[i] = 2.0 * tr_xz_xzz[i] - 4.0 * tr_xyz_xyzz[i] * tke_0 - 10.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_yyy[i] = 2.0 * tr_xz_yyy[i] + 6.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_yyyy[i] * tke_0 - 10.0 * tr_xyyz_yyy[i] * tbe_0 - 6.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_yyz[i] = 2.0 * tr_xz_yyz[i] + 4.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_yyyz[i] * tke_0 - 10.0 * tr_xyyz_yyz[i] * tbe_0 - 4.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_yzz[i] = 2.0 * tr_xz_yzz[i] + 2.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_yyzz[i] * tke_0 - 10.0 * tr_xyyz_yzz[i] * tbe_0 - 2.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyyz_zzz[i] = 2.0 * tr_xz_zzz[i] - 4.0 * tr_xyz_yzzz[i] * tke_0 - 10.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 680-690 components of targeted buffer : GF

    auto tr_y_0_y_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 680);

    auto tr_y_0_y_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 681);

    auto tr_y_0_y_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 682);

    auto tr_y_0_y_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 683);

    auto tr_y_0_y_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 684);

    auto tr_y_0_y_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 685);

    auto tr_y_0_y_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 686);

    auto tr_y_0_y_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 687);

    auto tr_y_0_y_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 688);

    auto tr_y_0_y_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 689);

    #pragma omp simd aligned(tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxy, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xy, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yz, tr_xyyzz_yzzz, tr_xyyzz_zz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xzz_xx, tr_xzz_xxxy, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_yy, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_y_0_y_xyzz_xxx, tr_y_0_y_xyzz_xxy, tr_y_0_y_xyzz_xxz, tr_y_0_y_xyzz_xyy, tr_y_0_y_xyzz_xyz, tr_y_0_y_xyzz_xzz, tr_y_0_y_xyzz_yyy, tr_y_0_y_xyzz_yyz, tr_y_0_y_xyzz_yzz, tr_y_0_y_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyzz_xxx[i] = -2.0 * tr_xzz_xxxy[i] * tke_0 - 6.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_xxy[i] = tr_xzz_xx[i] - 2.0 * tr_xzz_xxyy[i] * tke_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 - 2.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_xxz[i] = -2.0 * tr_xzz_xxyz[i] * tke_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_xyy[i] = 2.0 * tr_xzz_xy[i] - 2.0 * tr_xzz_xyyy[i] * tke_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 - 4.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_xyz[i] = tr_xzz_xz[i] - 2.0 * tr_xzz_xyyz[i] * tke_0 - 6.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_xzz[i] = -2.0 * tr_xzz_xyzz[i] * tke_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_yyy[i] = 3.0 * tr_xzz_yy[i] - 2.0 * tr_xzz_yyyy[i] * tke_0 - 6.0 * tr_xyzz_yyy[i] * tbe_0 - 6.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_yyz[i] = 2.0 * tr_xzz_yz[i] - 2.0 * tr_xzz_yyyz[i] * tke_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 - 4.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_yzz[i] = tr_xzz_zz[i] - 2.0 * tr_xzz_yyzz[i] * tke_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 - 2.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyzz_zzz[i] = -2.0 * tr_xzz_yzzz[i] * tke_0 - 6.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 690-700 components of targeted buffer : GF

    auto tr_y_0_y_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 690);

    auto tr_y_0_y_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 691);

    auto tr_y_0_y_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 692);

    auto tr_y_0_y_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 693);

    auto tr_y_0_y_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 694);

    auto tr_y_0_y_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 695);

    auto tr_y_0_y_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 696);

    auto tr_y_0_y_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 697);

    auto tr_y_0_y_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 698);

    auto tr_y_0_y_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 699);

    #pragma omp simd aligned(tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxy, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xy, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yz, tr_xyzzz_yzzz, tr_xyzzz_zz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_y_0_y_xzzz_xxx, tr_y_0_y_xzzz_xxy, tr_y_0_y_xzzz_xxz, tr_y_0_y_xzzz_xyy, tr_y_0_y_xzzz_xyz, tr_y_0_y_xzzz_xzz, tr_y_0_y_xzzz_yyy, tr_y_0_y_xzzz_yyz, tr_y_0_y_xzzz_yzz, tr_y_0_y_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xzzz_xxx[i] = -2.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_xxy[i] = -2.0 * tr_xzzz_xxy[i] * tbe_0 - 2.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_xxz[i] = -2.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_xyy[i] = -2.0 * tr_xzzz_xyy[i] * tbe_0 - 4.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_xyz[i] = -2.0 * tr_xzzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_xzz[i] = -2.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_yyy[i] = -2.0 * tr_xzzz_yyy[i] * tbe_0 - 6.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_yyz[i] = -2.0 * tr_xzzz_yyz[i] * tbe_0 - 4.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_yzz[i] = -2.0 * tr_xzzz_yzz[i] * tbe_0 - 2.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzzz_zzz[i] = -2.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 700-710 components of targeted buffer : GF

    auto tr_y_0_y_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 700);

    auto tr_y_0_y_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 701);

    auto tr_y_0_y_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 702);

    auto tr_y_0_y_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 703);

    auto tr_y_0_y_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 704);

    auto tr_y_0_y_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 705);

    auto tr_y_0_y_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 706);

    auto tr_y_0_y_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 707);

    auto tr_y_0_y_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 708);

    auto tr_y_0_y_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 709);

    #pragma omp simd aligned(tr_y_0_y_yyyy_xxx, tr_y_0_y_yyyy_xxy, tr_y_0_y_yyyy_xxz, tr_y_0_y_yyyy_xyy, tr_y_0_y_yyyy_xyz, tr_y_0_y_yyyy_xzz, tr_y_0_y_yyyy_yyy, tr_y_0_y_yyyy_yyz, tr_y_0_y_yyyy_yzz, tr_y_0_y_yyyy_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyy_xx, tr_yyy_xxxy, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_yy, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, tr_yyyyy_xx, tr_yyyyy_xxxy, tr_yyyyy_xxyy, tr_yyyyy_xxyz, tr_yyyyy_xy, tr_yyyyy_xyyy, tr_yyyyy_xyyz, tr_yyyyy_xyzz, tr_yyyyy_xz, tr_yyyyy_yy, tr_yyyyy_yyyy, tr_yyyyy_yyyz, tr_yyyyy_yyzz, tr_yyyyy_yz, tr_yyyyy_yzzz, tr_yyyyy_zz, tr_yyyyyy_xxx, tr_yyyyyy_xxy, tr_yyyyyy_xxz, tr_yyyyyy_xyy, tr_yyyyyy_xyz, tr_yyyyyy_xzz, tr_yyyyyy_yyy, tr_yyyyyy_yyz, tr_yyyyyy_yzz, tr_yyyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyyy_xxx[i] = 12.0 * tr_yy_xxx[i] - 8.0 * tr_yyy_xxxy[i] * tke_0 - 18.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyyy_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_xxy[i] = 12.0 * tr_yy_xxy[i] + 4.0 * tr_yyy_xx[i] - 8.0 * tr_yyy_xxyy[i] * tke_0 - 18.0 * tr_yyyy_xxy[i] * tbe_0 - 2.0 * tr_yyyyy_xx[i] * tbe_0 + 4.0 * tr_yyyyy_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_xxz[i] = 12.0 * tr_yy_xxz[i] - 8.0 * tr_yyy_xxyz[i] * tke_0 - 18.0 * tr_yyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_xyy[i] = 12.0 * tr_yy_xyy[i] + 8.0 * tr_yyy_xy[i] - 8.0 * tr_yyy_xyyy[i] * tke_0 - 18.0 * tr_yyyy_xyy[i] * tbe_0 - 4.0 * tr_yyyyy_xy[i] * tbe_0 + 4.0 * tr_yyyyy_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_xyz[i] = 12.0 * tr_yy_xyz[i] + 4.0 * tr_yyy_xz[i] - 8.0 * tr_yyy_xyyz[i] * tke_0 - 18.0 * tr_yyyy_xyz[i] * tbe_0 - 2.0 * tr_yyyyy_xz[i] * tbe_0 + 4.0 * tr_yyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_xzz[i] = 12.0 * tr_yy_xzz[i] - 8.0 * tr_yyy_xyzz[i] * tke_0 - 18.0 * tr_yyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_yyy[i] = 12.0 * tr_yy_yyy[i] + 12.0 * tr_yyy_yy[i] - 8.0 * tr_yyy_yyyy[i] * tke_0 - 18.0 * tr_yyyy_yyy[i] * tbe_0 - 6.0 * tr_yyyyy_yy[i] * tbe_0 + 4.0 * tr_yyyyy_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_yyz[i] = 12.0 * tr_yy_yyz[i] + 8.0 * tr_yyy_yz[i] - 8.0 * tr_yyy_yyyz[i] * tke_0 - 18.0 * tr_yyyy_yyz[i] * tbe_0 - 4.0 * tr_yyyyy_yz[i] * tbe_0 + 4.0 * tr_yyyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_yzz[i] = 12.0 * tr_yy_yzz[i] + 4.0 * tr_yyy_zz[i] - 8.0 * tr_yyy_yyzz[i] * tke_0 - 18.0 * tr_yyyy_yzz[i] * tbe_0 - 2.0 * tr_yyyyy_zz[i] * tbe_0 + 4.0 * tr_yyyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyy_zzz[i] = 12.0 * tr_yy_zzz[i] - 8.0 * tr_yyy_yzzz[i] * tke_0 - 18.0 * tr_yyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyy_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 710-720 components of targeted buffer : GF

    auto tr_y_0_y_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 710);

    auto tr_y_0_y_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 711);

    auto tr_y_0_y_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 712);

    auto tr_y_0_y_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 713);

    auto tr_y_0_y_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 714);

    auto tr_y_0_y_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 715);

    auto tr_y_0_y_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 716);

    auto tr_y_0_y_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 717);

    auto tr_y_0_y_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 718);

    auto tr_y_0_y_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 719);

    #pragma omp simd aligned(tr_y_0_y_yyyz_xxx, tr_y_0_y_yyyz_xxy, tr_y_0_y_yyyz_xxz, tr_y_0_y_yyyz_xyy, tr_y_0_y_yyyz_xyz, tr_y_0_y_yyyz_xzz, tr_y_0_y_yyyz_yyy, tr_y_0_y_yyyz_yyz, tr_y_0_y_yyyz_yzz, tr_y_0_y_yyyz_zzz, tr_yyyyyz_xxx, tr_yyyyyz_xxy, tr_yyyyyz_xxz, tr_yyyyyz_xyy, tr_yyyyyz_xyz, tr_yyyyyz_xzz, tr_yyyyyz_yyy, tr_yyyyyz_yyz, tr_yyyyyz_yzz, tr_yyyyyz_zzz, tr_yyyyz_xx, tr_yyyyz_xxxy, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xy, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yz, tr_yyyyz_yzzz, tr_yyyyz_zz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyz_xx, tr_yyz_xxxy, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_yy, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyyz_xxx[i] = 6.0 * tr_yz_xxx[i] - 6.0 * tr_yyz_xxxy[i] * tke_0 - 14.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_xxy[i] = 6.0 * tr_yz_xxy[i] + 3.0 * tr_yyz_xx[i] - 6.0 * tr_yyz_xxyy[i] * tke_0 - 14.0 * tr_yyyz_xxy[i] * tbe_0 - 2.0 * tr_yyyyz_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_xxz[i] = 6.0 * tr_yz_xxz[i] - 6.0 * tr_yyz_xxyz[i] * tke_0 - 14.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_xyy[i] = 6.0 * tr_yz_xyy[i] + 6.0 * tr_yyz_xy[i] - 6.0 * tr_yyz_xyyy[i] * tke_0 - 14.0 * tr_yyyz_xyy[i] * tbe_0 - 4.0 * tr_yyyyz_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_xyz[i] = 6.0 * tr_yz_xyz[i] + 3.0 * tr_yyz_xz[i] - 6.0 * tr_yyz_xyyz[i] * tke_0 - 14.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_yyyyz_xz[i] * tbe_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_xzz[i] = 6.0 * tr_yz_xzz[i] - 6.0 * tr_yyz_xyzz[i] * tke_0 - 14.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_yyy[i] = 6.0 * tr_yz_yyy[i] + 9.0 * tr_yyz_yy[i] - 6.0 * tr_yyz_yyyy[i] * tke_0 - 14.0 * tr_yyyz_yyy[i] * tbe_0 - 6.0 * tr_yyyyz_yy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_yyz[i] = 6.0 * tr_yz_yyz[i] + 6.0 * tr_yyz_yz[i] - 6.0 * tr_yyz_yyyz[i] * tke_0 - 14.0 * tr_yyyz_yyz[i] * tbe_0 - 4.0 * tr_yyyyz_yz[i] * tbe_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_yzz[i] = 6.0 * tr_yz_yzz[i] + 3.0 * tr_yyz_zz[i] - 6.0 * tr_yyz_yyzz[i] * tke_0 - 14.0 * tr_yyyz_yzz[i] * tbe_0 - 2.0 * tr_yyyyz_zz[i] * tbe_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyyz_zzz[i] = 6.0 * tr_yz_zzz[i] - 6.0 * tr_yyz_yzzz[i] * tke_0 - 14.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 720-730 components of targeted buffer : GF

    auto tr_y_0_y_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 720);

    auto tr_y_0_y_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 721);

    auto tr_y_0_y_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 722);

    auto tr_y_0_y_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 723);

    auto tr_y_0_y_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 724);

    auto tr_y_0_y_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 725);

    auto tr_y_0_y_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 726);

    auto tr_y_0_y_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 727);

    auto tr_y_0_y_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 728);

    auto tr_y_0_y_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 729);

    #pragma omp simd aligned(tr_y_0_y_yyzz_xxx, tr_y_0_y_yyzz_xxy, tr_y_0_y_yyzz_xxz, tr_y_0_y_yyzz_xyy, tr_y_0_y_yyzz_xyz, tr_y_0_y_yyzz_xzz, tr_y_0_y_yyzz_yyy, tr_y_0_y_yyzz_yyz, tr_y_0_y_yyzz_yzz, tr_y_0_y_yyzz_zzz, tr_yyyyzz_xxx, tr_yyyyzz_xxy, tr_yyyyzz_xxz, tr_yyyyzz_xyy, tr_yyyyzz_xyz, tr_yyyyzz_xzz, tr_yyyyzz_yyy, tr_yyyyzz_yyz, tr_yyyyzz_yzz, tr_yyyyzz_zzz, tr_yyyzz_xx, tr_yyyzz_xxxy, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xy, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yz, tr_yyyzz_yzzz, tr_yyyzz_zz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yzz_xx, tr_yzz_xxxy, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_yy, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyzz_xxx[i] = 2.0 * tr_zz_xxx[i] - 4.0 * tr_yzz_xxxy[i] * tke_0 - 10.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_xxy[i] = 2.0 * tr_zz_xxy[i] + 2.0 * tr_yzz_xx[i] - 4.0 * tr_yzz_xxyy[i] * tke_0 - 10.0 * tr_yyzz_xxy[i] * tbe_0 - 2.0 * tr_yyyzz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_xxz[i] = 2.0 * tr_zz_xxz[i] - 4.0 * tr_yzz_xxyz[i] * tke_0 - 10.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_xyy[i] = 2.0 * tr_zz_xyy[i] + 4.0 * tr_yzz_xy[i] - 4.0 * tr_yzz_xyyy[i] * tke_0 - 10.0 * tr_yyzz_xyy[i] * tbe_0 - 4.0 * tr_yyyzz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_xyz[i] = 2.0 * tr_zz_xyz[i] + 2.0 * tr_yzz_xz[i] - 4.0 * tr_yzz_xyyz[i] * tke_0 - 10.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_yyyzz_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_xzz[i] = 2.0 * tr_zz_xzz[i] - 4.0 * tr_yzz_xyzz[i] * tke_0 - 10.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_yyy[i] = 2.0 * tr_zz_yyy[i] + 6.0 * tr_yzz_yy[i] - 4.0 * tr_yzz_yyyy[i] * tke_0 - 10.0 * tr_yyzz_yyy[i] * tbe_0 - 6.0 * tr_yyyzz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_yyz[i] = 2.0 * tr_zz_yyz[i] + 4.0 * tr_yzz_yz[i] - 4.0 * tr_yzz_yyyz[i] * tke_0 - 10.0 * tr_yyzz_yyz[i] * tbe_0 - 4.0 * tr_yyyzz_yz[i] * tbe_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_yzz[i] = 2.0 * tr_zz_yzz[i] + 2.0 * tr_yzz_zz[i] - 4.0 * tr_yzz_yyzz[i] * tke_0 - 10.0 * tr_yyzz_yzz[i] * tbe_0 - 2.0 * tr_yyyzz_zz[i] * tbe_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyzz_zzz[i] = 2.0 * tr_zz_zzz[i] - 4.0 * tr_yzz_yzzz[i] * tke_0 - 10.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 730-740 components of targeted buffer : GF

    auto tr_y_0_y_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 730);

    auto tr_y_0_y_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 731);

    auto tr_y_0_y_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 732);

    auto tr_y_0_y_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 733);

    auto tr_y_0_y_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 734);

    auto tr_y_0_y_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 735);

    auto tr_y_0_y_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 736);

    auto tr_y_0_y_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 737);

    auto tr_y_0_y_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 738);

    auto tr_y_0_y_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 739);

    #pragma omp simd aligned(tr_y_0_y_yzzz_xxx, tr_y_0_y_yzzz_xxy, tr_y_0_y_yzzz_xxz, tr_y_0_y_yzzz_xyy, tr_y_0_y_yzzz_xyz, tr_y_0_y_yzzz_xzz, tr_y_0_y_yzzz_yyy, tr_y_0_y_yzzz_yyz, tr_y_0_y_yzzz_yzz, tr_y_0_y_yzzz_zzz, tr_yyyzzz_xxx, tr_yyyzzz_xxy, tr_yyyzzz_xxz, tr_yyyzzz_xyy, tr_yyyzzz_xyz, tr_yyyzzz_xzz, tr_yyyzzz_yyy, tr_yyyzzz_yyz, tr_yyyzzz_yzz, tr_yyyzzz_zzz, tr_yyzzz_xx, tr_yyzzz_xxxy, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xy, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yz, tr_yyzzz_yzzz, tr_yyzzz_zz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_zzz_xx, tr_zzz_xxxy, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_yy, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yzzz_xxx[i] = -2.0 * tr_zzz_xxxy[i] * tke_0 - 6.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_xxy[i] = tr_zzz_xx[i] - 2.0 * tr_zzz_xxyy[i] * tke_0 - 6.0 * tr_yzzz_xxy[i] * tbe_0 - 2.0 * tr_yyzzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_xxz[i] = -2.0 * tr_zzz_xxyz[i] * tke_0 - 6.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_xyy[i] = 2.0 * tr_zzz_xy[i] - 2.0 * tr_zzz_xyyy[i] * tke_0 - 6.0 * tr_yzzz_xyy[i] * tbe_0 - 4.0 * tr_yyzzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_xyz[i] = tr_zzz_xz[i] - 2.0 * tr_zzz_xyyz[i] * tke_0 - 6.0 * tr_yzzz_xyz[i] * tbe_0 - 2.0 * tr_yyzzz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_xzz[i] = -2.0 * tr_zzz_xyzz[i] * tke_0 - 6.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_yyy[i] = 3.0 * tr_zzz_yy[i] - 2.0 * tr_zzz_yyyy[i] * tke_0 - 6.0 * tr_yzzz_yyy[i] * tbe_0 - 6.0 * tr_yyzzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_yyz[i] = 2.0 * tr_zzz_yz[i] - 2.0 * tr_zzz_yyyz[i] * tke_0 - 6.0 * tr_yzzz_yyz[i] * tbe_0 - 4.0 * tr_yyzzz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_yzz[i] = tr_zzz_zz[i] - 2.0 * tr_zzz_yyzz[i] * tke_0 - 6.0 * tr_yzzz_yzz[i] * tbe_0 - 2.0 * tr_yyzzz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzzz_zzz[i] = -2.0 * tr_zzz_yzzz[i] * tke_0 - 6.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 740-750 components of targeted buffer : GF

    auto tr_y_0_y_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 740);

    auto tr_y_0_y_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 741);

    auto tr_y_0_y_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 742);

    auto tr_y_0_y_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 743);

    auto tr_y_0_y_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 744);

    auto tr_y_0_y_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 745);

    auto tr_y_0_y_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 746);

    auto tr_y_0_y_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 747);

    auto tr_y_0_y_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 748);

    auto tr_y_0_y_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 749);

    #pragma omp simd aligned(tr_y_0_y_zzzz_xxx, tr_y_0_y_zzzz_xxy, tr_y_0_y_zzzz_xxz, tr_y_0_y_zzzz_xyy, tr_y_0_y_zzzz_xyz, tr_y_0_y_zzzz_xzz, tr_y_0_y_zzzz_yyy, tr_y_0_y_zzzz_yyz, tr_y_0_y_zzzz_yzz, tr_y_0_y_zzzz_zzz, tr_yyzzzz_xxx, tr_yyzzzz_xxy, tr_yyzzzz_xxz, tr_yyzzzz_xyy, tr_yyzzzz_xyz, tr_yyzzzz_xzz, tr_yyzzzz_yyy, tr_yyzzzz_yyz, tr_yyzzzz_yzz, tr_yyzzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xxxy, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xy, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yz, tr_yzzzz_yzzz, tr_yzzzz_zz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zzzz_xxx[i] = -2.0 * tr_zzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_xxy[i] = -2.0 * tr_zzzz_xxy[i] * tbe_0 - 2.0 * tr_yzzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_xxz[i] = -2.0 * tr_zzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_xyy[i] = -2.0 * tr_zzzz_xyy[i] * tbe_0 - 4.0 * tr_yzzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_xyz[i] = -2.0 * tr_zzzz_xyz[i] * tbe_0 - 2.0 * tr_yzzzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_xzz[i] = -2.0 * tr_zzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_yyy[i] = -2.0 * tr_zzzz_yyy[i] * tbe_0 - 6.0 * tr_yzzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_yyz[i] = -2.0 * tr_zzzz_yyz[i] * tbe_0 - 4.0 * tr_yzzzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_yzz[i] = -2.0 * tr_zzzz_yzz[i] * tbe_0 - 2.0 * tr_yzzzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzzz_zzz[i] = -2.0 * tr_zzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 750-760 components of targeted buffer : GF

    auto tr_y_0_z_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 750);

    auto tr_y_0_z_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 751);

    auto tr_y_0_z_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 752);

    auto tr_y_0_z_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 753);

    auto tr_y_0_z_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 754);

    auto tr_y_0_z_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 755);

    auto tr_y_0_z_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 756);

    auto tr_y_0_z_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 757);

    auto tr_y_0_z_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 758);

    auto tr_y_0_z_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 759);

    #pragma omp simd aligned(tr_xxxxy_xx, tr_xxxxy_xxxz, tr_xxxxy_xxyz, tr_xxxxy_xxzz, tr_xxxxy_xy, tr_xxxxy_xyyz, tr_xxxxy_xyzz, tr_xxxxy_xz, tr_xxxxy_xzzz, tr_xxxxy_yy, tr_xxxxy_yyyz, tr_xxxxy_yyzz, tr_xxxxy_yz, tr_xxxxy_yzzz, tr_xxxxy_zz, tr_xxxxy_zzzz, tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_y_0_z_xxxx_xxx, tr_y_0_z_xxxx_xxy, tr_y_0_z_xxxx_xxz, tr_y_0_z_xxxx_xyy, tr_y_0_z_xxxx_xyz, tr_y_0_z_xxxx_xzz, tr_y_0_z_xxxx_yyy, tr_y_0_z_xxxx_yyz, tr_y_0_z_xxxx_yzz, tr_y_0_z_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxxx_xxx[i] = 4.0 * tr_xxxxy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_xxy[i] = 4.0 * tr_xxxxy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_xxz[i] = -2.0 * tr_xxxxy_xx[i] * tbe_0 + 4.0 * tr_xxxxy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_xyy[i] = 4.0 * tr_xxxxy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_xyz[i] = -2.0 * tr_xxxxy_xy[i] * tbe_0 + 4.0 * tr_xxxxy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_xzz[i] = -4.0 * tr_xxxxy_xz[i] * tbe_0 + 4.0 * tr_xxxxy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_yyy[i] = 4.0 * tr_xxxxy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_yyz[i] = -2.0 * tr_xxxxy_yy[i] * tbe_0 + 4.0 * tr_xxxxy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_yzz[i] = -4.0 * tr_xxxxy_yz[i] * tbe_0 + 4.0 * tr_xxxxy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxx_zzz[i] = -6.0 * tr_xxxxy_zz[i] * tbe_0 + 4.0 * tr_xxxxy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 760-770 components of targeted buffer : GF

    auto tr_y_0_z_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 760);

    auto tr_y_0_z_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 761);

    auto tr_y_0_z_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 762);

    auto tr_y_0_z_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 763);

    auto tr_y_0_z_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 764);

    auto tr_y_0_z_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 765);

    auto tr_y_0_z_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 766);

    auto tr_y_0_z_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 767);

    auto tr_y_0_z_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 768);

    auto tr_y_0_z_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 769);

    #pragma omp simd aligned(tr_xxx_xx, tr_xxx_xxxz, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxx_zzzz, tr_xxxyy_xx, tr_xxxyy_xxxz, tr_xxxyy_xxyz, tr_xxxyy_xxzz, tr_xxxyy_xy, tr_xxxyy_xyyz, tr_xxxyy_xyzz, tr_xxxyy_xz, tr_xxxyy_xzzz, tr_xxxyy_yy, tr_xxxyy_yyyz, tr_xxxyy_yyzz, tr_xxxyy_yz, tr_xxxyy_yzzz, tr_xxxyy_zz, tr_xxxyy_zzzz, tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_y_0_z_xxxy_xxx, tr_y_0_z_xxxy_xxy, tr_y_0_z_xxxy_xxz, tr_y_0_z_xxxy_xyy, tr_y_0_z_xxxy_xyz, tr_y_0_z_xxxy_xzz, tr_y_0_z_xxxy_yyy, tr_y_0_z_xxxy_yyz, tr_y_0_z_xxxy_yzz, tr_y_0_z_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxxy_xxx[i] = -2.0 * tr_xxx_xxxz[i] * tke_0 - 2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_xxy[i] = -2.0 * tr_xxx_xxyz[i] * tke_0 - 2.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_xxz[i] = tr_xxx_xx[i] - 2.0 * tr_xxx_xxzz[i] * tke_0 - 2.0 * tr_xxxz_xxz[i] * tbe_0 - 2.0 * tr_xxxyy_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_xyy[i] = -2.0 * tr_xxx_xyyz[i] * tke_0 - 2.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_xyz[i] = tr_xxx_xy[i] - 2.0 * tr_xxx_xyzz[i] * tke_0 - 2.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxyy_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_xzz[i] = 2.0 * tr_xxx_xz[i] - 2.0 * tr_xxx_xzzz[i] * tke_0 - 2.0 * tr_xxxz_xzz[i] * tbe_0 - 4.0 * tr_xxxyy_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_yyy[i] = -2.0 * tr_xxx_yyyz[i] * tke_0 - 2.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_yyz[i] = tr_xxx_yy[i] - 2.0 * tr_xxx_yyzz[i] * tke_0 - 2.0 * tr_xxxz_yyz[i] * tbe_0 - 2.0 * tr_xxxyy_yy[i] * tbe_0 + 4.0 * tr_xxxyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_yzz[i] = 2.0 * tr_xxx_yz[i] - 2.0 * tr_xxx_yzzz[i] * tke_0 - 2.0 * tr_xxxz_yzz[i] * tbe_0 - 4.0 * tr_xxxyy_yz[i] * tbe_0 + 4.0 * tr_xxxyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxy_zzz[i] = 3.0 * tr_xxx_zz[i] - 2.0 * tr_xxx_zzzz[i] * tke_0 - 2.0 * tr_xxxz_zzz[i] * tbe_0 - 6.0 * tr_xxxyy_zz[i] * tbe_0 + 4.0 * tr_xxxyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 770-780 components of targeted buffer : GF

    auto tr_y_0_z_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 770);

    auto tr_y_0_z_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 771);

    auto tr_y_0_z_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 772);

    auto tr_y_0_z_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 773);

    auto tr_y_0_z_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 774);

    auto tr_y_0_z_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 775);

    auto tr_y_0_z_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 776);

    auto tr_y_0_z_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 777);

    auto tr_y_0_z_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 778);

    auto tr_y_0_z_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 779);

    #pragma omp simd aligned(tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyz_xx, tr_xxxyz_xxxz, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_xzzz, tr_xxxyz_yy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yz, tr_xxxyz_yzzz, tr_xxxyz_zz, tr_xxxyz_zzzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_y_0_z_xxxz_xxx, tr_y_0_z_xxxz_xxy, tr_y_0_z_xxxz_xxz, tr_y_0_z_xxxz_xyy, tr_y_0_z_xxxz_xyz, tr_y_0_z_xxxz_xzz, tr_y_0_z_xxxz_yyy, tr_y_0_z_xxxz_yyz, tr_y_0_z_xxxz_yzz, tr_y_0_z_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxxz_xxx[i] = -2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_xxy[i] = -2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_xxz[i] = -2.0 * tr_xxxy_xxz[i] * tbe_0 - 2.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_xyy[i] = -2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_xyz[i] = -2.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_xzz[i] = -2.0 * tr_xxxy_xzz[i] * tbe_0 - 4.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_yyy[i] = -2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_yyz[i] = -2.0 * tr_xxxy_yyz[i] * tbe_0 - 2.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_yzz[i] = -2.0 * tr_xxxy_yzz[i] * tbe_0 - 4.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxxz_zzz[i] = -2.0 * tr_xxxy_zzz[i] * tbe_0 - 6.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 780-790 components of targeted buffer : GF

    auto tr_y_0_z_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 780);

    auto tr_y_0_z_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 781);

    auto tr_y_0_z_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 782);

    auto tr_y_0_z_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 783);

    auto tr_y_0_z_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 784);

    auto tr_y_0_z_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 785);

    auto tr_y_0_z_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 786);

    auto tr_y_0_z_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 787);

    auto tr_y_0_z_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 788);

    auto tr_y_0_z_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 789);

    #pragma omp simd aligned(tr_xxy_xx, tr_xxy_xxxz, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxy_zzzz, tr_xxyyy_xx, tr_xxyyy_xxxz, tr_xxyyy_xxyz, tr_xxyyy_xxzz, tr_xxyyy_xy, tr_xxyyy_xyyz, tr_xxyyy_xyzz, tr_xxyyy_xz, tr_xxyyy_xzzz, tr_xxyyy_yy, tr_xxyyy_yyyz, tr_xxyyy_yyzz, tr_xxyyy_yz, tr_xxyyy_yzzz, tr_xxyyy_zz, tr_xxyyy_zzzz, tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_y_0_z_xxyy_xxx, tr_y_0_z_xxyy_xxy, tr_y_0_z_xxyy_xxz, tr_y_0_z_xxyy_xyy, tr_y_0_z_xxyy_xyz, tr_y_0_z_xxyy_xzz, tr_y_0_z_xxyy_yyy, tr_y_0_z_xxyy_yyz, tr_y_0_z_xxyy_yzz, tr_y_0_z_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxyy_xxx[i] = -4.0 * tr_xxy_xxxz[i] * tke_0 - 4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_xxy[i] = -4.0 * tr_xxy_xxyz[i] * tke_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_xxz[i] = 2.0 * tr_xxy_xx[i] - 4.0 * tr_xxy_xxzz[i] * tke_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 - 2.0 * tr_xxyyy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_xyy[i] = -4.0 * tr_xxy_xyyz[i] * tke_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_xyz[i] = 2.0 * tr_xxy_xy[i] - 4.0 * tr_xxy_xyzz[i] * tke_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_xzz[i] = 4.0 * tr_xxy_xz[i] - 4.0 * tr_xxy_xzzz[i] * tke_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 - 4.0 * tr_xxyyy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_yyy[i] = -4.0 * tr_xxy_yyyz[i] * tke_0 - 4.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_yyz[i] = 2.0 * tr_xxy_yy[i] - 4.0 * tr_xxy_yyzz[i] * tke_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 - 2.0 * tr_xxyyy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_yzz[i] = 4.0 * tr_xxy_yz[i] - 4.0 * tr_xxy_yzzz[i] * tke_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 - 4.0 * tr_xxyyy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyy_zzz[i] = 6.0 * tr_xxy_zz[i] - 4.0 * tr_xxy_zzzz[i] * tke_0 - 4.0 * tr_xxyz_zzz[i] * tbe_0 - 6.0 * tr_xxyyy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 790-800 components of targeted buffer : GF

    auto tr_y_0_z_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 790);

    auto tr_y_0_z_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 791);

    auto tr_y_0_z_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 792);

    auto tr_y_0_z_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 793);

    auto tr_y_0_z_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 794);

    auto tr_y_0_z_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 795);

    auto tr_y_0_z_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 796);

    auto tr_y_0_z_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 797);

    auto tr_y_0_z_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 798);

    auto tr_y_0_z_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 799);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyz_xx, tr_xxyyz_xxxz, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_xzzz, tr_xxyyz_yy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yz, tr_xxyyz_yzzz, tr_xxyyz_zz, tr_xxyyz_zzzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xxz_xx, tr_xxz_xxxz, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_y_0_z_xxyz_xxx, tr_y_0_z_xxyz_xxy, tr_y_0_z_xxyz_xxz, tr_y_0_z_xxyz_xyy, tr_y_0_z_xxyz_xyz, tr_y_0_z_xxyz_xzz, tr_y_0_z_xxyz_yyy, tr_y_0_z_xxyz_yyz, tr_y_0_z_xxyz_yzz, tr_y_0_z_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxyz_xxx[i] = tr_xx_xxx[i] - 2.0 * tr_xxz_xxxz[i] * tke_0 - 2.0 * tr_xxzz_xxx[i] * tbe_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_xxy[i] = tr_xx_xxy[i] - 2.0 * tr_xxz_xxyz[i] * tke_0 - 2.0 * tr_xxzz_xxy[i] * tbe_0 - 2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_xxz[i] = tr_xx_xxz[i] + tr_xxz_xx[i] - 2.0 * tr_xxz_xxzz[i] * tke_0 - 2.0 * tr_xxzz_xxz[i] * tbe_0 - 2.0 * tr_xxyy_xxz[i] * tbe_0 - 2.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_xyy[i] = tr_xx_xyy[i] - 2.0 * tr_xxz_xyyz[i] * tke_0 - 2.0 * tr_xxzz_xyy[i] * tbe_0 - 2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_xyz[i] = tr_xx_xyz[i] + tr_xxz_xy[i] - 2.0 * tr_xxz_xyzz[i] * tke_0 - 2.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_xzz[i] = tr_xx_xzz[i] + 2.0 * tr_xxz_xz[i] - 2.0 * tr_xxz_xzzz[i] * tke_0 - 2.0 * tr_xxzz_xzz[i] * tbe_0 - 2.0 * tr_xxyy_xzz[i] * tbe_0 - 4.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_yyy[i] = tr_xx_yyy[i] - 2.0 * tr_xxz_yyyz[i] * tke_0 - 2.0 * tr_xxzz_yyy[i] * tbe_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_yyz[i] = tr_xx_yyz[i] + tr_xxz_yy[i] - 2.0 * tr_xxz_yyzz[i] * tke_0 - 2.0 * tr_xxzz_yyz[i] * tbe_0 - 2.0 * tr_xxyy_yyz[i] * tbe_0 - 2.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_yzz[i] = tr_xx_yzz[i] + 2.0 * tr_xxz_yz[i] - 2.0 * tr_xxz_yzzz[i] * tke_0 - 2.0 * tr_xxzz_yzz[i] * tbe_0 - 2.0 * tr_xxyy_yzz[i] * tbe_0 - 4.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxyz_zzz[i] = tr_xx_zzz[i] + 3.0 * tr_xxz_zz[i] - 2.0 * tr_xxz_zzzz[i] * tke_0 - 2.0 * tr_xxzz_zzz[i] * tbe_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 - 6.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 800-810 components of targeted buffer : GF

    auto tr_y_0_z_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 800);

    auto tr_y_0_z_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 801);

    auto tr_y_0_z_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 802);

    auto tr_y_0_z_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 803);

    auto tr_y_0_z_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 804);

    auto tr_y_0_z_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 805);

    auto tr_y_0_z_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 806);

    auto tr_y_0_z_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 807);

    auto tr_y_0_z_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 808);

    auto tr_y_0_z_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 809);

    #pragma omp simd aligned(tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxz, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_xzzz, tr_xxyzz_yy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yz, tr_xxyzz_yzzz, tr_xxyzz_zz, tr_xxyzz_zzzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_y_0_z_xxzz_xxx, tr_y_0_z_xxzz_xxy, tr_y_0_z_xxzz_xxz, tr_y_0_z_xxzz_xyy, tr_y_0_z_xxzz_xyz, tr_y_0_z_xxzz_xzz, tr_y_0_z_xxzz_yyy, tr_y_0_z_xxzz_yyz, tr_y_0_z_xxzz_yzz, tr_y_0_z_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxzz_xxx[i] = -4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_xxy[i] = -4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_xxz[i] = -4.0 * tr_xxyz_xxz[i] * tbe_0 - 2.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_xyy[i] = -4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_xyz[i] = -4.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_xzz[i] = -4.0 * tr_xxyz_xzz[i] * tbe_0 - 4.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_yyy[i] = -4.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_yyz[i] = -4.0 * tr_xxyz_yyz[i] * tbe_0 - 2.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_yzz[i] = -4.0 * tr_xxyz_yzz[i] * tbe_0 - 4.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxzz_zzz[i] = -4.0 * tr_xxyz_zzz[i] * tbe_0 - 6.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 810-820 components of targeted buffer : GF

    auto tr_y_0_z_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 810);

    auto tr_y_0_z_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 811);

    auto tr_y_0_z_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 812);

    auto tr_y_0_z_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 813);

    auto tr_y_0_z_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 814);

    auto tr_y_0_z_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 815);

    auto tr_y_0_z_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 816);

    auto tr_y_0_z_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 817);

    auto tr_y_0_z_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 818);

    auto tr_y_0_z_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 819);

    #pragma omp simd aligned(tr_xyy_xx, tr_xyy_xxxz, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyy_zzzz, tr_xyyyy_xx, tr_xyyyy_xxxz, tr_xyyyy_xxyz, tr_xyyyy_xxzz, tr_xyyyy_xy, tr_xyyyy_xyyz, tr_xyyyy_xyzz, tr_xyyyy_xz, tr_xyyyy_xzzz, tr_xyyyy_yy, tr_xyyyy_yyyz, tr_xyyyy_yyzz, tr_xyyyy_yz, tr_xyyyy_yzzz, tr_xyyyy_zz, tr_xyyyy_zzzz, tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_y_0_z_xyyy_xxx, tr_y_0_z_xyyy_xxy, tr_y_0_z_xyyy_xxz, tr_y_0_z_xyyy_xyy, tr_y_0_z_xyyy_xyz, tr_y_0_z_xyyy_xzz, tr_y_0_z_xyyy_yyy, tr_y_0_z_xyyy_yyz, tr_y_0_z_xyyy_yzz, tr_y_0_z_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyyy_xxx[i] = -6.0 * tr_xyy_xxxz[i] * tke_0 - 6.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_xxy[i] = -6.0 * tr_xyy_xxyz[i] * tke_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_xxz[i] = 3.0 * tr_xyy_xx[i] - 6.0 * tr_xyy_xxzz[i] * tke_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 - 2.0 * tr_xyyyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_xyy[i] = -6.0 * tr_xyy_xyyz[i] * tke_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_xyz[i] = 3.0 * tr_xyy_xy[i] - 6.0 * tr_xyy_xyzz[i] * tke_0 - 6.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_xzz[i] = 6.0 * tr_xyy_xz[i] - 6.0 * tr_xyy_xzzz[i] * tke_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 - 4.0 * tr_xyyyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_yyy[i] = -6.0 * tr_xyy_yyyz[i] * tke_0 - 6.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_yyz[i] = 3.0 * tr_xyy_yy[i] - 6.0 * tr_xyy_yyzz[i] * tke_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xyyyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_yzz[i] = 6.0 * tr_xyy_yz[i] - 6.0 * tr_xyy_yzzz[i] * tke_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 - 4.0 * tr_xyyyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyy_zzz[i] = 9.0 * tr_xyy_zz[i] - 6.0 * tr_xyy_zzzz[i] * tke_0 - 6.0 * tr_xyyz_zzz[i] * tbe_0 - 6.0 * tr_xyyyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 820-830 components of targeted buffer : GF

    auto tr_y_0_z_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 820);

    auto tr_y_0_z_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 821);

    auto tr_y_0_z_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 822);

    auto tr_y_0_z_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 823);

    auto tr_y_0_z_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 824);

    auto tr_y_0_z_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 825);

    auto tr_y_0_z_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 826);

    auto tr_y_0_z_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 827);

    auto tr_y_0_z_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 828);

    auto tr_y_0_z_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 829);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyz_xx, tr_xyyyz_xxxz, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_xzzz, tr_xyyyz_yy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yz, tr_xyyyz_yzzz, tr_xyyyz_zz, tr_xyyyz_zzzz, tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_y_0_z_xyyz_xxx, tr_y_0_z_xyyz_xxy, tr_y_0_z_xyyz_xxz, tr_y_0_z_xyyz_xyy, tr_y_0_z_xyyz_xyz, tr_y_0_z_xyyz_xzz, tr_y_0_z_xyyz_yyy, tr_y_0_z_xyyz_yyz, tr_y_0_z_xyyz_yzz, tr_y_0_z_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyyz_xxx[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xyz_xxxz[i] * tke_0 - 4.0 * tr_xyzz_xxx[i] * tbe_0 - 2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_xxy[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 4.0 * tr_xyzz_xxy[i] * tbe_0 - 2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_xxz[i] = 2.0 * tr_xy_xxz[i] + 2.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxzz[i] * tke_0 - 4.0 * tr_xyzz_xxz[i] * tbe_0 - 2.0 * tr_xyyy_xxz[i] * tbe_0 - 2.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_xyy[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 4.0 * tr_xyzz_xyy[i] * tbe_0 - 2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_xyz[i] = 2.0 * tr_xy_xyz[i] + 2.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xyzz[i] * tke_0 - 4.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xyyy_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_xzz[i] = 2.0 * tr_xy_xzz[i] + 4.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xzzz[i] * tke_0 - 4.0 * tr_xyzz_xzz[i] * tbe_0 - 2.0 * tr_xyyy_xzz[i] * tbe_0 - 4.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_yyy[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xyz_yyyz[i] * tke_0 - 4.0 * tr_xyzz_yyy[i] * tbe_0 - 2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_yyz[i] = 2.0 * tr_xy_yyz[i] + 2.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_yyzz[i] * tke_0 - 4.0 * tr_xyzz_yyz[i] * tbe_0 - 2.0 * tr_xyyy_yyz[i] * tbe_0 - 2.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_yzz[i] = 2.0 * tr_xy_yzz[i] + 4.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_yzzz[i] * tke_0 - 4.0 * tr_xyzz_yzz[i] * tbe_0 - 2.0 * tr_xyyy_yzz[i] * tbe_0 - 4.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyyz_zzz[i] = 2.0 * tr_xy_zzz[i] + 6.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_zzzz[i] * tke_0 - 4.0 * tr_xyzz_zzz[i] * tbe_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 - 6.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 830-840 components of targeted buffer : GF

    auto tr_y_0_z_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 830);

    auto tr_y_0_z_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 831);

    auto tr_y_0_z_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 832);

    auto tr_y_0_z_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 833);

    auto tr_y_0_z_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 834);

    auto tr_y_0_z_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 835);

    auto tr_y_0_z_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 836);

    auto tr_y_0_z_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 837);

    auto tr_y_0_z_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 838);

    auto tr_y_0_z_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 839);

    #pragma omp simd aligned(tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxz, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_xzzz, tr_xyyzz_yy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yz, tr_xyyzz_yzzz, tr_xyyzz_zz, tr_xyyzz_zzzz, tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxz, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_y_0_z_xyzz_xxx, tr_y_0_z_xyzz_xxy, tr_y_0_z_xyzz_xxz, tr_y_0_z_xyzz_xyy, tr_y_0_z_xyzz_xyz, tr_y_0_z_xyzz_xzz, tr_y_0_z_xyzz_yyy, tr_y_0_z_xyzz_yyz, tr_y_0_z_xyzz_yzz, tr_y_0_z_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyzz_xxx[i] = 2.0 * tr_xz_xxx[i] - 2.0 * tr_xzz_xxxz[i] * tke_0 - 2.0 * tr_xzzz_xxx[i] * tbe_0 - 4.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_xxy[i] = 2.0 * tr_xz_xxy[i] - 2.0 * tr_xzz_xxyz[i] * tke_0 - 2.0 * tr_xzzz_xxy[i] * tbe_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_xxz[i] = 2.0 * tr_xz_xxz[i] + tr_xzz_xx[i] - 2.0 * tr_xzz_xxzz[i] * tke_0 - 2.0 * tr_xzzz_xxz[i] * tbe_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 - 2.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_xyy[i] = 2.0 * tr_xz_xyy[i] - 2.0 * tr_xzz_xyyz[i] * tke_0 - 2.0 * tr_xzzz_xyy[i] * tbe_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_xyz[i] = 2.0 * tr_xz_xyz[i] + tr_xzz_xy[i] - 2.0 * tr_xzz_xyzz[i] * tke_0 - 2.0 * tr_xzzz_xyz[i] * tbe_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_xzz[i] = 2.0 * tr_xz_xzz[i] + 2.0 * tr_xzz_xz[i] - 2.0 * tr_xzz_xzzz[i] * tke_0 - 2.0 * tr_xzzz_xzz[i] * tbe_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 - 4.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_yyy[i] = 2.0 * tr_xz_yyy[i] - 2.0 * tr_xzz_yyyz[i] * tke_0 - 2.0 * tr_xzzz_yyy[i] * tbe_0 - 4.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_yyz[i] = 2.0 * tr_xz_yyz[i] + tr_xzz_yy[i] - 2.0 * tr_xzz_yyzz[i] * tke_0 - 2.0 * tr_xzzz_yyz[i] * tbe_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_yzz[i] = 2.0 * tr_xz_yzz[i] + 2.0 * tr_xzz_yz[i] - 2.0 * tr_xzz_yzzz[i] * tke_0 - 2.0 * tr_xzzz_yzz[i] * tbe_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 - 4.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyzz_zzz[i] = 2.0 * tr_xz_zzz[i] + 3.0 * tr_xzz_zz[i] - 2.0 * tr_xzz_zzzz[i] * tke_0 - 2.0 * tr_xzzz_zzz[i] * tbe_0 - 4.0 * tr_xyyz_zzz[i] * tbe_0 - 6.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 840-850 components of targeted buffer : GF

    auto tr_y_0_z_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 840);

    auto tr_y_0_z_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 841);

    auto tr_y_0_z_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 842);

    auto tr_y_0_z_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 843);

    auto tr_y_0_z_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 844);

    auto tr_y_0_z_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 845);

    auto tr_y_0_z_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 846);

    auto tr_y_0_z_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 847);

    auto tr_y_0_z_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 848);

    auto tr_y_0_z_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 849);

    #pragma omp simd aligned(tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxz, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_xzzz, tr_xyzzz_yy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yz, tr_xyzzz_yzzz, tr_xyzzz_zz, tr_xyzzz_zzzz, tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, tr_y_0_z_xzzz_xxx, tr_y_0_z_xzzz_xxy, tr_y_0_z_xzzz_xxz, tr_y_0_z_xzzz_xyy, tr_y_0_z_xzzz_xyz, tr_y_0_z_xzzz_xzz, tr_y_0_z_xzzz_yyy, tr_y_0_z_xzzz_yyz, tr_y_0_z_xzzz_yzz, tr_y_0_z_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xzzz_xxx[i] = -6.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_xxy[i] = -6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_xxz[i] = -6.0 * tr_xyzz_xxz[i] * tbe_0 - 2.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_xyy[i] = -6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_xyz[i] = -6.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_xzz[i] = -6.0 * tr_xyzz_xzz[i] * tbe_0 - 4.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_yyy[i] = -6.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_yyz[i] = -6.0 * tr_xyzz_yyz[i] * tbe_0 - 2.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_yzz[i] = -6.0 * tr_xyzz_yzz[i] * tbe_0 - 4.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzzz_zzz[i] = -6.0 * tr_xyzz_zzz[i] * tbe_0 - 6.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 850-860 components of targeted buffer : GF

    auto tr_y_0_z_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 850);

    auto tr_y_0_z_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 851);

    auto tr_y_0_z_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 852);

    auto tr_y_0_z_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 853);

    auto tr_y_0_z_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 854);

    auto tr_y_0_z_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 855);

    auto tr_y_0_z_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 856);

    auto tr_y_0_z_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 857);

    auto tr_y_0_z_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 858);

    auto tr_y_0_z_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 859);

    #pragma omp simd aligned(tr_y_0_z_yyyy_xxx, tr_y_0_z_yyyy_xxy, tr_y_0_z_yyyy_xxz, tr_y_0_z_yyyy_xyy, tr_y_0_z_yyyy_xyz, tr_y_0_z_yyyy_xzz, tr_y_0_z_yyyy_yyy, tr_y_0_z_yyyy_yyz, tr_y_0_z_yyyy_yzz, tr_y_0_z_yyyy_zzz, tr_yyy_xx, tr_yyy_xxxz, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyy_zzzz, tr_yyyyy_xx, tr_yyyyy_xxxz, tr_yyyyy_xxyz, tr_yyyyy_xxzz, tr_yyyyy_xy, tr_yyyyy_xyyz, tr_yyyyy_xyzz, tr_yyyyy_xz, tr_yyyyy_xzzz, tr_yyyyy_yy, tr_yyyyy_yyyz, tr_yyyyy_yyzz, tr_yyyyy_yz, tr_yyyyy_yzzz, tr_yyyyy_zz, tr_yyyyy_zzzz, tr_yyyyyz_xxx, tr_yyyyyz_xxy, tr_yyyyyz_xxz, tr_yyyyyz_xyy, tr_yyyyyz_xyz, tr_yyyyyz_xzz, tr_yyyyyz_yyy, tr_yyyyyz_yyz, tr_yyyyyz_yzz, tr_yyyyyz_zzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyyy_xxx[i] = -8.0 * tr_yyy_xxxz[i] * tke_0 - 8.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyy_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_xxy[i] = -8.0 * tr_yyy_xxyz[i] * tke_0 - 8.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyyy_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_xxz[i] = 4.0 * tr_yyy_xx[i] - 8.0 * tr_yyy_xxzz[i] * tke_0 - 8.0 * tr_yyyz_xxz[i] * tbe_0 - 2.0 * tr_yyyyy_xx[i] * tbe_0 + 4.0 * tr_yyyyy_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_xyy[i] = -8.0 * tr_yyy_xyyz[i] * tke_0 - 8.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyyy_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_xyz[i] = 4.0 * tr_yyy_xy[i] - 8.0 * tr_yyy_xyzz[i] * tke_0 - 8.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_yyyyy_xy[i] * tbe_0 + 4.0 * tr_yyyyy_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_xzz[i] = 8.0 * tr_yyy_xz[i] - 8.0 * tr_yyy_xzzz[i] * tke_0 - 8.0 * tr_yyyz_xzz[i] * tbe_0 - 4.0 * tr_yyyyy_xz[i] * tbe_0 + 4.0 * tr_yyyyy_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_yyy[i] = -8.0 * tr_yyy_yyyz[i] * tke_0 - 8.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyyy_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_yyz[i] = 4.0 * tr_yyy_yy[i] - 8.0 * tr_yyy_yyzz[i] * tke_0 - 8.0 * tr_yyyz_yyz[i] * tbe_0 - 2.0 * tr_yyyyy_yy[i] * tbe_0 + 4.0 * tr_yyyyy_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_yzz[i] = 8.0 * tr_yyy_yz[i] - 8.0 * tr_yyy_yzzz[i] * tke_0 - 8.0 * tr_yyyz_yzz[i] * tbe_0 - 4.0 * tr_yyyyy_yz[i] * tbe_0 + 4.0 * tr_yyyyy_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyy_zzz[i] = 12.0 * tr_yyy_zz[i] - 8.0 * tr_yyy_zzzz[i] * tke_0 - 8.0 * tr_yyyz_zzz[i] * tbe_0 - 6.0 * tr_yyyyy_zz[i] * tbe_0 + 4.0 * tr_yyyyy_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 860-870 components of targeted buffer : GF

    auto tr_y_0_z_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 860);

    auto tr_y_0_z_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 861);

    auto tr_y_0_z_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 862);

    auto tr_y_0_z_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 863);

    auto tr_y_0_z_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 864);

    auto tr_y_0_z_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 865);

    auto tr_y_0_z_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 866);

    auto tr_y_0_z_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 867);

    auto tr_y_0_z_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 868);

    auto tr_y_0_z_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 869);

    #pragma omp simd aligned(tr_y_0_z_yyyz_xxx, tr_y_0_z_yyyz_xxy, tr_y_0_z_yyyz_xxz, tr_y_0_z_yyyz_xyy, tr_y_0_z_yyyz_xyz, tr_y_0_z_yyyz_xzz, tr_y_0_z_yyyz_yyy, tr_y_0_z_yyyz_yyz, tr_y_0_z_yyyz_yzz, tr_y_0_z_yyyz_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, tr_yyyyz_xx, tr_yyyyz_xxxz, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xz, tr_yyyyz_xzzz, tr_yyyyz_yy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yz, tr_yyyyz_yzzz, tr_yyyyz_zz, tr_yyyyz_zzzz, tr_yyyyzz_xxx, tr_yyyyzz_xxy, tr_yyyyzz_xxz, tr_yyyyzz_xyy, tr_yyyyzz_xyz, tr_yyyyzz_xzz, tr_yyyyzz_yyy, tr_yyyyzz_yyz, tr_yyyyzz_yzz, tr_yyyyzz_zzz, tr_yyz_xx, tr_yyz_xxxz, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyyz_xxx[i] = 3.0 * tr_yy_xxx[i] - 6.0 * tr_yyz_xxxz[i] * tke_0 - 6.0 * tr_yyzz_xxx[i] * tbe_0 - 2.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_xxy[i] = 3.0 * tr_yy_xxy[i] - 6.0 * tr_yyz_xxyz[i] * tke_0 - 6.0 * tr_yyzz_xxy[i] * tbe_0 - 2.0 * tr_yyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_xxz[i] = 3.0 * tr_yy_xxz[i] + 3.0 * tr_yyz_xx[i] - 6.0 * tr_yyz_xxzz[i] * tke_0 - 6.0 * tr_yyzz_xxz[i] * tbe_0 - 2.0 * tr_yyyy_xxz[i] * tbe_0 - 2.0 * tr_yyyyz_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_xyy[i] = 3.0 * tr_yy_xyy[i] - 6.0 * tr_yyz_xyyz[i] * tke_0 - 6.0 * tr_yyzz_xyy[i] * tbe_0 - 2.0 * tr_yyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_xyz[i] = 3.0 * tr_yy_xyz[i] + 3.0 * tr_yyz_xy[i] - 6.0 * tr_yyz_xyzz[i] * tke_0 - 6.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_yyyy_xyz[i] * tbe_0 - 2.0 * tr_yyyyz_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_xzz[i] = 3.0 * tr_yy_xzz[i] + 6.0 * tr_yyz_xz[i] - 6.0 * tr_yyz_xzzz[i] * tke_0 - 6.0 * tr_yyzz_xzz[i] * tbe_0 - 2.0 * tr_yyyy_xzz[i] * tbe_0 - 4.0 * tr_yyyyz_xz[i] * tbe_0 + 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_yyy[i] = 3.0 * tr_yy_yyy[i] - 6.0 * tr_yyz_yyyz[i] * tke_0 - 6.0 * tr_yyzz_yyy[i] * tbe_0 - 2.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_yyz[i] = 3.0 * tr_yy_yyz[i] + 3.0 * tr_yyz_yy[i] - 6.0 * tr_yyz_yyzz[i] * tke_0 - 6.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_yyyy_yyz[i] * tbe_0 - 2.0 * tr_yyyyz_yy[i] * tbe_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_yzz[i] = 3.0 * tr_yy_yzz[i] + 6.0 * tr_yyz_yz[i] - 6.0 * tr_yyz_yzzz[i] * tke_0 - 6.0 * tr_yyzz_yzz[i] * tbe_0 - 2.0 * tr_yyyy_yzz[i] * tbe_0 - 4.0 * tr_yyyyz_yz[i] * tbe_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyyz_zzz[i] = 3.0 * tr_yy_zzz[i] + 9.0 * tr_yyz_zz[i] - 6.0 * tr_yyz_zzzz[i] * tke_0 - 6.0 * tr_yyzz_zzz[i] * tbe_0 - 2.0 * tr_yyyy_zzz[i] * tbe_0 - 6.0 * tr_yyyyz_zz[i] * tbe_0 + 4.0 * tr_yyyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 870-880 components of targeted buffer : GF

    auto tr_y_0_z_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 870);

    auto tr_y_0_z_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 871);

    auto tr_y_0_z_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 872);

    auto tr_y_0_z_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 873);

    auto tr_y_0_z_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 874);

    auto tr_y_0_z_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 875);

    auto tr_y_0_z_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 876);

    auto tr_y_0_z_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 877);

    auto tr_y_0_z_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 878);

    auto tr_y_0_z_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 879);

    #pragma omp simd aligned(tr_y_0_z_yyzz_xxx, tr_y_0_z_yyzz_xxy, tr_y_0_z_yyzz_xxz, tr_y_0_z_yyzz_xyy, tr_y_0_z_yyzz_xyz, tr_y_0_z_yyzz_xzz, tr_y_0_z_yyzz_yyy, tr_y_0_z_yyzz_yyz, tr_y_0_z_yyzz_yzz, tr_y_0_z_yyzz_zzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyyzz_xx, tr_yyyzz_xxxz, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xz, tr_yyyzz_xzzz, tr_yyyzz_yy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yz, tr_yyyzz_yzzz, tr_yyyzz_zz, tr_yyyzz_zzzz, tr_yyyzzz_xxx, tr_yyyzzz_xxy, tr_yyyzzz_xxz, tr_yyyzzz_xyy, tr_yyyzzz_xyz, tr_yyyzzz_xzz, tr_yyyzzz_yyy, tr_yyyzzz_yyz, tr_yyyzzz_yzz, tr_yyyzzz_zzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxz, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyzz_xxx[i] = 4.0 * tr_yz_xxx[i] - 4.0 * tr_yzz_xxxz[i] * tke_0 - 4.0 * tr_yzzz_xxx[i] * tbe_0 - 4.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_xxy[i] = 4.0 * tr_yz_xxy[i] - 4.0 * tr_yzz_xxyz[i] * tke_0 - 4.0 * tr_yzzz_xxy[i] * tbe_0 - 4.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_xxz[i] = 4.0 * tr_yz_xxz[i] + 2.0 * tr_yzz_xx[i] - 4.0 * tr_yzz_xxzz[i] * tke_0 - 4.0 * tr_yzzz_xxz[i] * tbe_0 - 4.0 * tr_yyyz_xxz[i] * tbe_0 - 2.0 * tr_yyyzz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_xyy[i] = 4.0 * tr_yz_xyy[i] - 4.0 * tr_yzz_xyyz[i] * tke_0 - 4.0 * tr_yzzz_xyy[i] * tbe_0 - 4.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_xyz[i] = 4.0 * tr_yz_xyz[i] + 2.0 * tr_yzz_xy[i] - 4.0 * tr_yzz_xyzz[i] * tke_0 - 4.0 * tr_yzzz_xyz[i] * tbe_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_yyyzz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_xzz[i] = 4.0 * tr_yz_xzz[i] + 4.0 * tr_yzz_xz[i] - 4.0 * tr_yzz_xzzz[i] * tke_0 - 4.0 * tr_yzzz_xzz[i] * tbe_0 - 4.0 * tr_yyyz_xzz[i] * tbe_0 - 4.0 * tr_yyyzz_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_yyy[i] = 4.0 * tr_yz_yyy[i] - 4.0 * tr_yzz_yyyz[i] * tke_0 - 4.0 * tr_yzzz_yyy[i] * tbe_0 - 4.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_yyz[i] = 4.0 * tr_yz_yyz[i] + 2.0 * tr_yzz_yy[i] - 4.0 * tr_yzz_yyzz[i] * tke_0 - 4.0 * tr_yzzz_yyz[i] * tbe_0 - 4.0 * tr_yyyz_yyz[i] * tbe_0 - 2.0 * tr_yyyzz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_yzz[i] = 4.0 * tr_yz_yzz[i] + 4.0 * tr_yzz_yz[i] - 4.0 * tr_yzz_yzzz[i] * tke_0 - 4.0 * tr_yzzz_yzz[i] * tbe_0 - 4.0 * tr_yyyz_yzz[i] * tbe_0 - 4.0 * tr_yyyzz_yz[i] * tbe_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyzz_zzz[i] = 4.0 * tr_yz_zzz[i] + 6.0 * tr_yzz_zz[i] - 4.0 * tr_yzz_zzzz[i] * tke_0 - 4.0 * tr_yzzz_zzz[i] * tbe_0 - 4.0 * tr_yyyz_zzz[i] * tbe_0 - 6.0 * tr_yyyzz_zz[i] * tbe_0 + 4.0 * tr_yyyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 880-890 components of targeted buffer : GF

    auto tr_y_0_z_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 880);

    auto tr_y_0_z_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 881);

    auto tr_y_0_z_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 882);

    auto tr_y_0_z_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 883);

    auto tr_y_0_z_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 884);

    auto tr_y_0_z_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 885);

    auto tr_y_0_z_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 886);

    auto tr_y_0_z_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 887);

    auto tr_y_0_z_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 888);

    auto tr_y_0_z_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 889);

    #pragma omp simd aligned(tr_y_0_z_yzzz_xxx, tr_y_0_z_yzzz_xxy, tr_y_0_z_yzzz_xxz, tr_y_0_z_yzzz_xyy, tr_y_0_z_yzzz_xyz, tr_y_0_z_yzzz_xzz, tr_y_0_z_yzzz_yyy, tr_y_0_z_yzzz_yyz, tr_y_0_z_yzzz_yzz, tr_y_0_z_yzzz_zzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yyzzz_xx, tr_yyzzz_xxxz, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xz, tr_yyzzz_xzzz, tr_yyzzz_yy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yz, tr_yyzzz_yzzz, tr_yyzzz_zz, tr_yyzzz_zzzz, tr_yyzzzz_xxx, tr_yyzzzz_xxy, tr_yyzzzz_xxz, tr_yyzzzz_xyy, tr_yyzzzz_xyz, tr_yyzzzz_xzz, tr_yyzzzz_yyy, tr_yyzzzz_yyz, tr_yyzzzz_yzz, tr_yyzzzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzz_xx, tr_zzz_xxxz, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yzzz_xxx[i] = 3.0 * tr_zz_xxx[i] - 2.0 * tr_zzz_xxxz[i] * tke_0 - 2.0 * tr_zzzz_xxx[i] * tbe_0 - 6.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_xxy[i] = 3.0 * tr_zz_xxy[i] - 2.0 * tr_zzz_xxyz[i] * tke_0 - 2.0 * tr_zzzz_xxy[i] * tbe_0 - 6.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_xxz[i] = 3.0 * tr_zz_xxz[i] + tr_zzz_xx[i] - 2.0 * tr_zzz_xxzz[i] * tke_0 - 2.0 * tr_zzzz_xxz[i] * tbe_0 - 6.0 * tr_yyzz_xxz[i] * tbe_0 - 2.0 * tr_yyzzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_xyy[i] = 3.0 * tr_zz_xyy[i] - 2.0 * tr_zzz_xyyz[i] * tke_0 - 2.0 * tr_zzzz_xyy[i] * tbe_0 - 6.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_xyz[i] = 3.0 * tr_zz_xyz[i] + tr_zzz_xy[i] - 2.0 * tr_zzz_xyzz[i] * tke_0 - 2.0 * tr_zzzz_xyz[i] * tbe_0 - 6.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_yyzzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_xzz[i] = 3.0 * tr_zz_xzz[i] + 2.0 * tr_zzz_xz[i] - 2.0 * tr_zzz_xzzz[i] * tke_0 - 2.0 * tr_zzzz_xzz[i] * tbe_0 - 6.0 * tr_yyzz_xzz[i] * tbe_0 - 4.0 * tr_yyzzz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_yyy[i] = 3.0 * tr_zz_yyy[i] - 2.0 * tr_zzz_yyyz[i] * tke_0 - 2.0 * tr_zzzz_yyy[i] * tbe_0 - 6.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_yyz[i] = 3.0 * tr_zz_yyz[i] + tr_zzz_yy[i] - 2.0 * tr_zzz_yyzz[i] * tke_0 - 2.0 * tr_zzzz_yyz[i] * tbe_0 - 6.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_yyzzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_yzz[i] = 3.0 * tr_zz_yzz[i] + 2.0 * tr_zzz_yz[i] - 2.0 * tr_zzz_yzzz[i] * tke_0 - 2.0 * tr_zzzz_yzz[i] * tbe_0 - 6.0 * tr_yyzz_yzz[i] * tbe_0 - 4.0 * tr_yyzzz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzzz_zzz[i] = 3.0 * tr_zz_zzz[i] + 3.0 * tr_zzz_zz[i] - 2.0 * tr_zzz_zzzz[i] * tke_0 - 2.0 * tr_zzzz_zzz[i] * tbe_0 - 6.0 * tr_yyzz_zzz[i] * tbe_0 - 6.0 * tr_yyzzz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 890-900 components of targeted buffer : GF

    auto tr_y_0_z_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 890);

    auto tr_y_0_z_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 891);

    auto tr_y_0_z_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 892);

    auto tr_y_0_z_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 893);

    auto tr_y_0_z_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 894);

    auto tr_y_0_z_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 895);

    auto tr_y_0_z_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 896);

    auto tr_y_0_z_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 897);

    auto tr_y_0_z_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 898);

    auto tr_y_0_z_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 899);

    #pragma omp simd aligned(tr_y_0_z_zzzz_xxx, tr_y_0_z_zzzz_xxy, tr_y_0_z_zzzz_xxz, tr_y_0_z_zzzz_xyy, tr_y_0_z_zzzz_xyz, tr_y_0_z_zzzz_xzz, tr_y_0_z_zzzz_yyy, tr_y_0_z_zzzz_yyz, tr_y_0_z_zzzz_yzz, tr_y_0_z_zzzz_zzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xxxz, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xz, tr_yzzzz_xzzz, tr_yzzzz_yy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yz, tr_yzzzz_yzzz, tr_yzzzz_zz, tr_yzzzz_zzzz, tr_yzzzzz_xxx, tr_yzzzzz_xxy, tr_yzzzzz_xxz, tr_yzzzzz_xyy, tr_yzzzzz_xyz, tr_yzzzzz_xzz, tr_yzzzzz_yyy, tr_yzzzzz_yyz, tr_yzzzzz_yzz, tr_yzzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zzzz_xxx[i] = -8.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_xxy[i] = -8.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_xxz[i] = -8.0 * tr_yzzz_xxz[i] * tbe_0 - 2.0 * tr_yzzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_xyy[i] = -8.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_xyz[i] = -8.0 * tr_yzzz_xyz[i] * tbe_0 - 2.0 * tr_yzzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_xzz[i] = -8.0 * tr_yzzz_xzz[i] * tbe_0 - 4.0 * tr_yzzzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_yyy[i] = -8.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_yyz[i] = -8.0 * tr_yzzz_yyz[i] * tbe_0 - 2.0 * tr_yzzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_yzz[i] = -8.0 * tr_yzzz_yzz[i] * tbe_0 - 4.0 * tr_yzzzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzzz_zzz[i] = -8.0 * tr_yzzz_zzz[i] * tbe_0 - 6.0 * tr_yzzzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 900-910 components of targeted buffer : GF

    auto tr_z_0_x_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 900);

    auto tr_z_0_x_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 901);

    auto tr_z_0_x_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 902);

    auto tr_z_0_x_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 903);

    auto tr_z_0_x_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 904);

    auto tr_z_0_x_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 905);

    auto tr_z_0_x_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 906);

    auto tr_z_0_x_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 907);

    auto tr_z_0_x_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 908);

    auto tr_z_0_x_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 909);

    #pragma omp simd aligned(tr_xxxxxz_xxx, tr_xxxxxz_xxy, tr_xxxxxz_xxz, tr_xxxxxz_xyy, tr_xxxxxz_xyz, tr_xxxxxz_xzz, tr_xxxxxz_yyy, tr_xxxxxz_yyz, tr_xxxxxz_yzz, tr_xxxxxz_zzz, tr_xxxxz_xx, tr_xxxxz_xxxx, tr_xxxxz_xxxy, tr_xxxxz_xxxz, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xy, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xz, tr_xxxxz_xzzz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_z_0_x_xxxx_xxx, tr_z_0_x_xxxx_xxy, tr_z_0_x_xxxx_xxz, tr_z_0_x_xxxx_xyy, tr_z_0_x_xxxx_xyz, tr_z_0_x_xxxx_xzz, tr_z_0_x_xxxx_yyy, tr_z_0_x_xxxx_yyz, tr_z_0_x_xxxx_yzz, tr_z_0_x_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxxx_xxx[i] = -8.0 * tr_xxxz_xxx[i] * tbe_0 - 6.0 * tr_xxxxz_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_xxy[i] = -8.0 * tr_xxxz_xxy[i] * tbe_0 - 4.0 * tr_xxxxz_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_xxz[i] = -8.0 * tr_xxxz_xxz[i] * tbe_0 - 4.0 * tr_xxxxz_xz[i] * tbe_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_xyy[i] = -8.0 * tr_xxxz_xyy[i] * tbe_0 - 2.0 * tr_xxxxz_yy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_xyz[i] = -8.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxxz_yz[i] * tbe_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_xzz[i] = -8.0 * tr_xxxz_xzz[i] * tbe_0 - 2.0 * tr_xxxxz_zz[i] * tbe_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_yyy[i] = -8.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_yyz[i] = -8.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_yzz[i] = -8.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxx_zzz[i] = -8.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxxz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 910-920 components of targeted buffer : GF

    auto tr_z_0_x_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 910);

    auto tr_z_0_x_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 911);

    auto tr_z_0_x_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 912);

    auto tr_z_0_x_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 913);

    auto tr_z_0_x_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 914);

    auto tr_z_0_x_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 915);

    auto tr_z_0_x_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 916);

    auto tr_z_0_x_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 917);

    auto tr_z_0_x_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 918);

    auto tr_z_0_x_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 919);

    #pragma omp simd aligned(tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_xxxyz_xx, tr_xxxyz_xxxx, tr_xxxyz_xxxy, tr_xxxyz_xxxz, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xy, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_xzzz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_z_0_x_xxxy_xxx, tr_z_0_x_xxxy_xxy, tr_z_0_x_xxxy_xxz, tr_z_0_x_xxxy_xyy, tr_z_0_x_xxxy_xyz, tr_z_0_x_xxxy_xzz, tr_z_0_x_xxxy_yyy, tr_z_0_x_xxxy_yyz, tr_z_0_x_xxxy_yzz, tr_z_0_x_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxxy_xxx[i] = -6.0 * tr_xxyz_xxx[i] * tbe_0 - 6.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_xxy[i] = -6.0 * tr_xxyz_xxy[i] * tbe_0 - 4.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_xxz[i] = -6.0 * tr_xxyz_xxz[i] * tbe_0 - 4.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_xyy[i] = -6.0 * tr_xxyz_xyy[i] * tbe_0 - 2.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_xyz[i] = -6.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_xzz[i] = -6.0 * tr_xxyz_xzz[i] * tbe_0 - 2.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_yyy[i] = -6.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_yyz[i] = -6.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_yzz[i] = -6.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxy_zzz[i] = -6.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 920-930 components of targeted buffer : GF

    auto tr_z_0_x_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 920);

    auto tr_z_0_x_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 921);

    auto tr_z_0_x_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 922);

    auto tr_z_0_x_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 923);

    auto tr_z_0_x_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 924);

    auto tr_z_0_x_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 925);

    auto tr_z_0_x_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 926);

    auto tr_z_0_x_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 927);

    auto tr_z_0_x_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 928);

    auto tr_z_0_x_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 929);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xx, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxzz_xxx, tr_xxxxzz_xxy, tr_xxxxzz_xxz, tr_xxxxzz_xyy, tr_xxxxzz_xyz, tr_xxxxzz_xzz, tr_xxxxzz_yyy, tr_xxxxzz_yyz, tr_xxxxzz_yzz, tr_xxxxzz_zzz, tr_xxxzz_xx, tr_xxxzz_xxxx, tr_xxxzz_xxxy, tr_xxxzz_xxxz, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xy, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xz, tr_xxxzz_xzzz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_z_0_x_xxxz_xxx, tr_z_0_x_xxxz_xxy, tr_z_0_x_xxxz_xxz, tr_z_0_x_xxxz_xyy, tr_z_0_x_xxxz_xyz, tr_z_0_x_xxxz_xzz, tr_z_0_x_xxxz_yyy, tr_z_0_x_xxxz_yyz, tr_z_0_x_xxxz_yzz, tr_z_0_x_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxxz_xxx[i] = 3.0 * tr_xx_xxx[i] - 6.0 * tr_xxzz_xxx[i] * tbe_0 + 3.0 * tr_xxx_xx[i] - 2.0 * tr_xxx_xxxx[i] * tke_0 - 6.0 * tr_xxxzz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_xxy[i] = 3.0 * tr_xx_xxy[i] - 6.0 * tr_xxzz_xxy[i] * tbe_0 + 2.0 * tr_xxx_xy[i] - 2.0 * tr_xxx_xxxy[i] * tke_0 - 4.0 * tr_xxxzz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_xxz[i] = 3.0 * tr_xx_xxz[i] - 6.0 * tr_xxzz_xxz[i] * tbe_0 + 2.0 * tr_xxx_xz[i] - 2.0 * tr_xxx_xxxz[i] * tke_0 - 4.0 * tr_xxxzz_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xxz[i] * tbe_0 + 4.0 * tr_xxxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_xyy[i] = 3.0 * tr_xx_xyy[i] - 6.0 * tr_xxzz_xyy[i] * tbe_0 + tr_xxx_yy[i] - 2.0 * tr_xxx_xxyy[i] * tke_0 - 2.0 * tr_xxxzz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_xyz[i] = 3.0 * tr_xx_xyz[i] - 6.0 * tr_xxzz_xyz[i] * tbe_0 + tr_xxx_yz[i] - 2.0 * tr_xxx_xxyz[i] * tke_0 - 2.0 * tr_xxxzz_yz[i] * tbe_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xyz[i] * tbe_0 + 4.0 * tr_xxxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_xzz[i] = 3.0 * tr_xx_xzz[i] - 6.0 * tr_xxzz_xzz[i] * tbe_0 + tr_xxx_zz[i] - 2.0 * tr_xxx_xxzz[i] * tke_0 - 2.0 * tr_xxxzz_zz[i] * tbe_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_xzz[i] * tbe_0 + 4.0 * tr_xxxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_yyy[i] = 3.0 * tr_xx_yyy[i] - 6.0 * tr_xxzz_yyy[i] * tbe_0 - 2.0 * tr_xxx_xyyy[i] * tke_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_yyz[i] = 3.0 * tr_xx_yyz[i] - 6.0 * tr_xxzz_yyz[i] * tbe_0 - 2.0 * tr_xxx_xyyz[i] * tke_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yyz[i] * tbe_0 + 4.0 * tr_xxxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_yzz[i] = 3.0 * tr_xx_yzz[i] - 6.0 * tr_xxzz_yzz[i] * tbe_0 - 2.0 * tr_xxx_xyzz[i] * tke_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_yzz[i] * tbe_0 + 4.0 * tr_xxxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxxz_zzz[i] = 3.0 * tr_xx_zzz[i] - 6.0 * tr_xxzz_zzz[i] * tbe_0 - 2.0 * tr_xxx_xzzz[i] * tke_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxx_zzz[i] * tbe_0 + 4.0 * tr_xxxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 930-940 components of targeted buffer : GF

    auto tr_z_0_x_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 930);

    auto tr_z_0_x_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 931);

    auto tr_z_0_x_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 932);

    auto tr_z_0_x_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 933);

    auto tr_z_0_x_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 934);

    auto tr_z_0_x_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 935);

    auto tr_z_0_x_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 936);

    auto tr_z_0_x_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 937);

    auto tr_z_0_x_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 938);

    auto tr_z_0_x_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 939);

    #pragma omp simd aligned(tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxyyz_xx, tr_xxyyz_xxxx, tr_xxyyz_xxxy, tr_xxyyz_xxxz, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xy, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_xzzz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_z_0_x_xxyy_xxx, tr_z_0_x_xxyy_xxy, tr_z_0_x_xxyy_xxz, tr_z_0_x_xxyy_xyy, tr_z_0_x_xxyy_xyz, tr_z_0_x_xxyy_xzz, tr_z_0_x_xxyy_yyy, tr_z_0_x_xxyy_yyz, tr_z_0_x_xxyy_yzz, tr_z_0_x_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxyy_xxx[i] = -4.0 * tr_xyyz_xxx[i] * tbe_0 - 6.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_xxy[i] = -4.0 * tr_xyyz_xxy[i] * tbe_0 - 4.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_xxz[i] = -4.0 * tr_xyyz_xxz[i] * tbe_0 - 4.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_xyy[i] = -4.0 * tr_xyyz_xyy[i] * tbe_0 - 2.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_xyz[i] = -4.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_xzz[i] = -4.0 * tr_xyyz_xzz[i] * tbe_0 - 2.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_yyy[i] = -4.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_yyz[i] = -4.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_yzz[i] = -4.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyy_zzz[i] = -4.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 940-950 components of targeted buffer : GF

    auto tr_z_0_x_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 940);

    auto tr_z_0_x_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 941);

    auto tr_z_0_x_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 942);

    auto tr_z_0_x_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 943);

    auto tr_z_0_x_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 944);

    auto tr_z_0_x_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 945);

    auto tr_z_0_x_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 946);

    auto tr_z_0_x_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 947);

    auto tr_z_0_x_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 948);

    auto tr_z_0_x_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 949);

    #pragma omp simd aligned(tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_xxy_xx, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyzz_xx, tr_xxyzz_xxxx, tr_xxyzz_xxxy, tr_xxyzz_xxxz, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xy, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_xzzz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_z_0_x_xxyz_xxx, tr_z_0_x_xxyz_xxy, tr_z_0_x_xxyz_xxz, tr_z_0_x_xxyz_xyy, tr_z_0_x_xxyz_xyz, tr_z_0_x_xxyz_xzz, tr_z_0_x_xxyz_yyy, tr_z_0_x_xxyz_yyz, tr_z_0_x_xxyz_yzz, tr_z_0_x_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxyz_xxx[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xyzz_xxx[i] * tbe_0 + 3.0 * tr_xxy_xx[i] - 2.0 * tr_xxy_xxxx[i] * tke_0 - 6.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_xxy[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xyzz_xxy[i] * tbe_0 + 2.0 * tr_xxy_xy[i] - 2.0 * tr_xxy_xxxy[i] * tke_0 - 4.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_xxz[i] = 2.0 * tr_xy_xxz[i] - 4.0 * tr_xyzz_xxz[i] * tbe_0 + 2.0 * tr_xxy_xz[i] - 2.0 * tr_xxy_xxxz[i] * tke_0 - 4.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_xyy[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xyzz_xyy[i] * tbe_0 + tr_xxy_yy[i] - 2.0 * tr_xxy_xxyy[i] * tke_0 - 2.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_xyz[i] = 2.0 * tr_xy_xyz[i] - 4.0 * tr_xyzz_xyz[i] * tbe_0 + tr_xxy_yz[i] - 2.0 * tr_xxy_xxyz[i] * tke_0 - 2.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_xzz[i] = 2.0 * tr_xy_xzz[i] - 4.0 * tr_xyzz_xzz[i] * tbe_0 + tr_xxy_zz[i] - 2.0 * tr_xxy_xxzz[i] * tke_0 - 2.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_yyy[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xyzz_yyy[i] * tbe_0 - 2.0 * tr_xxy_xyyy[i] * tke_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_yyz[i] = 2.0 * tr_xy_yyz[i] - 4.0 * tr_xyzz_yyz[i] * tbe_0 - 2.0 * tr_xxy_xyyz[i] * tke_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_yzz[i] = 2.0 * tr_xy_yzz[i] - 4.0 * tr_xyzz_yzz[i] * tbe_0 - 2.0 * tr_xxy_xyzz[i] * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxyz_zzz[i] = 2.0 * tr_xy_zzz[i] - 4.0 * tr_xyzz_zzz[i] * tbe_0 - 2.0 * tr_xxy_xzzz[i] * tke_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 950-960 components of targeted buffer : GF

    auto tr_z_0_x_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 950);

    auto tr_z_0_x_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 951);

    auto tr_z_0_x_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 952);

    auto tr_z_0_x_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 953);

    auto tr_z_0_x_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 954);

    auto tr_z_0_x_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 955);

    auto tr_z_0_x_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 956);

    auto tr_z_0_x_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 957);

    auto tr_z_0_x_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 958);

    auto tr_z_0_x_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 959);

    #pragma omp simd aligned(tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxxzzz_xxx, tr_xxxzzz_xxy, tr_xxxzzz_xxz, tr_xxxzzz_xyy, tr_xxxzzz_xyz, tr_xxxzzz_xzz, tr_xxxzzz_yyy, tr_xxxzzz_yyz, tr_xxxzzz_yzz, tr_xxxzzz_zzz, tr_xxz_xx, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xxzzz_xx, tr_xxzzz_xxxx, tr_xxzzz_xxxy, tr_xxzzz_xxxz, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xy, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xz, tr_xxzzz_xzzz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_z_0_x_xxzz_xxx, tr_z_0_x_xxzz_xxy, tr_z_0_x_xxzz_xxz, tr_z_0_x_xxzz_xyy, tr_z_0_x_xxzz_xyz, tr_z_0_x_xxzz_xzz, tr_z_0_x_xxzz_yyy, tr_z_0_x_xxzz_yyz, tr_z_0_x_xxzz_yzz, tr_z_0_x_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxzz_xxx[i] = 4.0 * tr_xz_xxx[i] - 4.0 * tr_xzzz_xxx[i] * tbe_0 + 6.0 * tr_xxz_xx[i] - 4.0 * tr_xxz_xxxx[i] * tke_0 - 6.0 * tr_xxzzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_xxy[i] = 4.0 * tr_xz_xxy[i] - 4.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xxz_xy[i] - 4.0 * tr_xxz_xxxy[i] * tke_0 - 4.0 * tr_xxzzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_xxz[i] = 4.0 * tr_xz_xxz[i] - 4.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xxz_xz[i] - 4.0 * tr_xxz_xxxz[i] * tke_0 - 4.0 * tr_xxzzz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_xyy[i] = 4.0 * tr_xz_xyy[i] - 4.0 * tr_xzzz_xyy[i] * tbe_0 + 2.0 * tr_xxz_yy[i] - 4.0 * tr_xxz_xxyy[i] * tke_0 - 2.0 * tr_xxzzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_xyz[i] = 4.0 * tr_xz_xyz[i] - 4.0 * tr_xzzz_xyz[i] * tbe_0 + 2.0 * tr_xxz_yz[i] - 4.0 * tr_xxz_xxyz[i] * tke_0 - 2.0 * tr_xxzzz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xyz[i] * tbe_0 + 4.0 * tr_xxxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_xzz[i] = 4.0 * tr_xz_xzz[i] - 4.0 * tr_xzzz_xzz[i] * tbe_0 + 2.0 * tr_xxz_zz[i] - 4.0 * tr_xxz_xxzz[i] * tke_0 - 2.0 * tr_xxzzz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_yyy[i] = 4.0 * tr_xz_yyy[i] - 4.0 * tr_xzzz_yyy[i] * tbe_0 - 4.0 * tr_xxz_xyyy[i] * tke_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_yyz[i] = 4.0 * tr_xz_yyz[i] - 4.0 * tr_xzzz_yyz[i] * tbe_0 - 4.0 * tr_xxz_xyyz[i] * tke_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yyz[i] * tbe_0 + 4.0 * tr_xxxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_yzz[i] = 4.0 * tr_xz_yzz[i] - 4.0 * tr_xzzz_yzz[i] * tbe_0 - 4.0 * tr_xxz_xyzz[i] * tke_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_yzz[i] * tbe_0 + 4.0 * tr_xxxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxzz_zzz[i] = 4.0 * tr_xz_zzz[i] - 4.0 * tr_xzzz_zzz[i] * tbe_0 - 4.0 * tr_xxz_xzzz[i] * tke_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 960-970 components of targeted buffer : GF

    auto tr_z_0_x_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 960);

    auto tr_z_0_x_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 961);

    auto tr_z_0_x_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 962);

    auto tr_z_0_x_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 963);

    auto tr_z_0_x_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 964);

    auto tr_z_0_x_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 965);

    auto tr_z_0_x_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 966);

    auto tr_z_0_x_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 967);

    auto tr_z_0_x_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 968);

    auto tr_z_0_x_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 969);

    #pragma omp simd aligned(tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xyyyz_xx, tr_xyyyz_xxxx, tr_xyyyz_xxxy, tr_xyyyz_xxxz, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xy, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_xzzz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_z_0_x_xyyy_xxx, tr_z_0_x_xyyy_xxy, tr_z_0_x_xyyy_xxz, tr_z_0_x_xyyy_xyy, tr_z_0_x_xyyy_xyz, tr_z_0_x_xyyy_xzz, tr_z_0_x_xyyy_yyy, tr_z_0_x_xyyy_yyz, tr_z_0_x_xyyy_yzz, tr_z_0_x_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyyy_xxx[i] = -2.0 * tr_yyyz_xxx[i] * tbe_0 - 6.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_xxy[i] = -2.0 * tr_yyyz_xxy[i] * tbe_0 - 4.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_xxz[i] = -2.0 * tr_yyyz_xxz[i] * tbe_0 - 4.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_xyy[i] = -2.0 * tr_yyyz_xyy[i] * tbe_0 - 2.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_xyz[i] = -2.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_xzz[i] = -2.0 * tr_yyyz_xzz[i] * tbe_0 - 2.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_yyy[i] = -2.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_yyz[i] = -2.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_yzz[i] = -2.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyy_zzz[i] = -2.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 970-980 components of targeted buffer : GF

    auto tr_z_0_x_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 970);

    auto tr_z_0_x_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 971);

    auto tr_z_0_x_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 972);

    auto tr_z_0_x_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 973);

    auto tr_z_0_x_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 974);

    auto tr_z_0_x_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 975);

    auto tr_z_0_x_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 976);

    auto tr_z_0_x_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 977);

    auto tr_z_0_x_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 978);

    auto tr_z_0_x_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 979);

    #pragma omp simd aligned(tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xyy_xx, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyzz_xx, tr_xyyzz_xxxx, tr_xyyzz_xxxy, tr_xyyzz_xxxz, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xy, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_xzzz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_z_0_x_xyyz_xxx, tr_z_0_x_xyyz_xxy, tr_z_0_x_xyyz_xxz, tr_z_0_x_xyyz_xyy, tr_z_0_x_xyyz_xyz, tr_z_0_x_xyyz_xzz, tr_z_0_x_xyyz_yyy, tr_z_0_x_xyyz_yyz, tr_z_0_x_xyyz_yzz, tr_z_0_x_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyyz_xxx[i] = tr_yy_xxx[i] - 2.0 * tr_yyzz_xxx[i] * tbe_0 + 3.0 * tr_xyy_xx[i] - 2.0 * tr_xyy_xxxx[i] * tke_0 - 6.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_xxy[i] = tr_yy_xxy[i] - 2.0 * tr_yyzz_xxy[i] * tbe_0 + 2.0 * tr_xyy_xy[i] - 2.0 * tr_xyy_xxxy[i] * tke_0 - 4.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_xxz[i] = tr_yy_xxz[i] - 2.0 * tr_yyzz_xxz[i] * tbe_0 + 2.0 * tr_xyy_xz[i] - 2.0 * tr_xyy_xxxz[i] * tke_0 - 4.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_xyy[i] = tr_yy_xyy[i] - 2.0 * tr_yyzz_xyy[i] * tbe_0 + tr_xyy_yy[i] - 2.0 * tr_xyy_xxyy[i] * tke_0 - 2.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_xyz[i] = tr_yy_xyz[i] - 2.0 * tr_yyzz_xyz[i] * tbe_0 + tr_xyy_yz[i] - 2.0 * tr_xyy_xxyz[i] * tke_0 - 2.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_xzz[i] = tr_yy_xzz[i] - 2.0 * tr_yyzz_xzz[i] * tbe_0 + tr_xyy_zz[i] - 2.0 * tr_xyy_xxzz[i] * tke_0 - 2.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_yyy[i] = tr_yy_yyy[i] - 2.0 * tr_yyzz_yyy[i] * tbe_0 - 2.0 * tr_xyy_xyyy[i] * tke_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_yyz[i] = tr_yy_yyz[i] - 2.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_xyy_xyyz[i] * tke_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_yzz[i] = tr_yy_yzz[i] - 2.0 * tr_yyzz_yzz[i] * tbe_0 - 2.0 * tr_xyy_xyzz[i] * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyyz_zzz[i] = tr_yy_zzz[i] - 2.0 * tr_yyzz_zzz[i] * tbe_0 - 2.0 * tr_xyy_xzzz[i] * tke_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 980-990 components of targeted buffer : GF

    auto tr_z_0_x_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 980);

    auto tr_z_0_x_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 981);

    auto tr_z_0_x_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 982);

    auto tr_z_0_x_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 983);

    auto tr_z_0_x_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 984);

    auto tr_z_0_x_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 985);

    auto tr_z_0_x_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 986);

    auto tr_z_0_x_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 987);

    auto tr_z_0_x_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 988);

    auto tr_z_0_x_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 989);

    #pragma omp simd aligned(tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzzz_xx, tr_xyzzz_xxxx, tr_xyzzz_xxxy, tr_xyzzz_xxxz, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xy, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_xzzz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_z_0_x_xyzz_xxx, tr_z_0_x_xyzz_xxy, tr_z_0_x_xyzz_xxz, tr_z_0_x_xyzz_xyy, tr_z_0_x_xyzz_xyz, tr_z_0_x_xyzz_xzz, tr_z_0_x_xyzz_yyy, tr_z_0_x_xyzz_yyz, tr_z_0_x_xyzz_yzz, tr_z_0_x_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyzz_xxx[i] = 2.0 * tr_yz_xxx[i] - 2.0 * tr_yzzz_xxx[i] * tbe_0 + 6.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxxx[i] * tke_0 - 6.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_xxy[i] = 2.0 * tr_yz_xxy[i] - 2.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xxxy[i] * tke_0 - 4.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_xxz[i] = 2.0 * tr_yz_xxz[i] - 2.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xxxz[i] * tke_0 - 4.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_xyy[i] = 2.0 * tr_yz_xyy[i] - 2.0 * tr_yzzz_xyy[i] * tbe_0 + 2.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_xxyy[i] * tke_0 - 2.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_xyz[i] = 2.0 * tr_yz_xyz[i] - 2.0 * tr_yzzz_xyz[i] * tbe_0 + 2.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 2.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_xzz[i] = 2.0 * tr_yz_xzz[i] - 2.0 * tr_yzzz_xzz[i] * tbe_0 + 2.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_xxzz[i] * tke_0 - 2.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_yyy[i] = 2.0 * tr_yz_yyy[i] - 2.0 * tr_yzzz_yyy[i] * tbe_0 - 4.0 * tr_xyz_xyyy[i] * tke_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_yyz[i] = 2.0 * tr_yz_yyz[i] - 2.0 * tr_yzzz_yyz[i] * tbe_0 - 4.0 * tr_xyz_xyyz[i] * tke_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_yzz[i] = 2.0 * tr_yz_yzz[i] - 2.0 * tr_yzzz_yzz[i] * tbe_0 - 4.0 * tr_xyz_xyzz[i] * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyzz_zzz[i] = 2.0 * tr_yz_zzz[i] - 2.0 * tr_yzzz_zzz[i] * tbe_0 - 4.0 * tr_xyz_xzzz[i] * tke_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 990-1000 components of targeted buffer : GF

    auto tr_z_0_x_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 990);

    auto tr_z_0_x_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 991);

    auto tr_z_0_x_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 992);

    auto tr_z_0_x_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 993);

    auto tr_z_0_x_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 994);

    auto tr_z_0_x_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 995);

    auto tr_z_0_x_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 996);

    auto tr_z_0_x_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 997);

    auto tr_z_0_x_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 998);

    auto tr_z_0_x_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 999);

    #pragma omp simd aligned(tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xxzzzz_xxx, tr_xxzzzz_xxy, tr_xxzzzz_xxz, tr_xxzzzz_xyy, tr_xxzzzz_xyz, tr_xxzzzz_xzz, tr_xxzzzz_yyy, tr_xxzzzz_yyz, tr_xxzzzz_yzz, tr_xxzzzz_zzz, tr_xzz_xx, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_xzzzz_xx, tr_xzzzz_xxxx, tr_xzzzz_xxxy, tr_xzzzz_xxxz, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xy, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xz, tr_xzzzz_xzzz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, tr_z_0_x_xzzz_xxx, tr_z_0_x_xzzz_xxy, tr_z_0_x_xzzz_xxz, tr_z_0_x_xzzz_xyy, tr_z_0_x_xzzz_xyz, tr_z_0_x_xzzz_xzz, tr_z_0_x_xzzz_yyy, tr_z_0_x_xzzz_yyz, tr_z_0_x_xzzz_yzz, tr_z_0_x_xzzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xzzz_xxx[i] = 3.0 * tr_zz_xxx[i] - 2.0 * tr_zzzz_xxx[i] * tbe_0 + 9.0 * tr_xzz_xx[i] - 6.0 * tr_xzz_xxxx[i] * tke_0 - 6.0 * tr_xzzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_xxy[i] = 3.0 * tr_zz_xxy[i] - 2.0 * tr_zzzz_xxy[i] * tbe_0 + 6.0 * tr_xzz_xy[i] - 6.0 * tr_xzz_xxxy[i] * tke_0 - 4.0 * tr_xzzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_xxz[i] = 3.0 * tr_zz_xxz[i] - 2.0 * tr_zzzz_xxz[i] * tbe_0 + 6.0 * tr_xzz_xz[i] - 6.0 * tr_xzz_xxxz[i] * tke_0 - 4.0 * tr_xzzzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xxz[i] * tbe_0 + 4.0 * tr_xxzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_xyy[i] = 3.0 * tr_zz_xyy[i] - 2.0 * tr_zzzz_xyy[i] * tbe_0 + 3.0 * tr_xzz_yy[i] - 6.0 * tr_xzz_xxyy[i] * tke_0 - 2.0 * tr_xzzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_xyz[i] = 3.0 * tr_zz_xyz[i] - 2.0 * tr_zzzz_xyz[i] * tbe_0 + 3.0 * tr_xzz_yz[i] - 6.0 * tr_xzz_xxyz[i] * tke_0 - 2.0 * tr_xzzzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xyz[i] * tbe_0 + 4.0 * tr_xxzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_xzz[i] = 3.0 * tr_zz_xzz[i] - 2.0 * tr_zzzz_xzz[i] * tbe_0 + 3.0 * tr_xzz_zz[i] - 6.0 * tr_xzz_xxzz[i] * tke_0 - 2.0 * tr_xzzzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_xzz[i] * tbe_0 + 4.0 * tr_xxzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_yyy[i] = 3.0 * tr_zz_yyy[i] - 2.0 * tr_zzzz_yyy[i] * tbe_0 - 6.0 * tr_xzz_xyyy[i] * tke_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_yyz[i] = 3.0 * tr_zz_yyz[i] - 2.0 * tr_zzzz_yyz[i] * tbe_0 - 6.0 * tr_xzz_xyyz[i] * tke_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yyz[i] * tbe_0 + 4.0 * tr_xxzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_yzz[i] = 3.0 * tr_zz_yzz[i] - 2.0 * tr_zzzz_yzz[i] * tbe_0 - 6.0 * tr_xzz_xyzz[i] * tke_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_yzz[i] * tbe_0 + 4.0 * tr_xxzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzzz_zzz[i] = 3.0 * tr_zz_zzz[i] - 2.0 * tr_zzzz_zzz[i] * tbe_0 - 6.0 * tr_xzz_xzzz[i] * tke_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xxzz_zzz[i] * tbe_0 + 4.0 * tr_xxzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1000-1010 components of targeted buffer : GF

    auto tr_z_0_x_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1000);

    auto tr_z_0_x_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1001);

    auto tr_z_0_x_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1002);

    auto tr_z_0_x_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1003);

    auto tr_z_0_x_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1004);

    auto tr_z_0_x_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1005);

    auto tr_z_0_x_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1006);

    auto tr_z_0_x_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1007);

    auto tr_z_0_x_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1008);

    auto tr_z_0_x_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1009);

    #pragma omp simd aligned(tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, tr_yyyyz_xx, tr_yyyyz_xxxx, tr_yyyyz_xxxy, tr_yyyyz_xxxz, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xy, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xz, tr_yyyyz_xzzz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, tr_z_0_x_yyyy_xxx, tr_z_0_x_yyyy_xxy, tr_z_0_x_yyyy_xxz, tr_z_0_x_yyyy_xyy, tr_z_0_x_yyyy_xyz, tr_z_0_x_yyyy_xzz, tr_z_0_x_yyyy_yyy, tr_z_0_x_yyyy_yyz, tr_z_0_x_yyyy_yzz, tr_z_0_x_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyyy_xxx[i] = -6.0 * tr_yyyyz_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_xxy[i] = -4.0 * tr_yyyyz_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_xxz[i] = -4.0 * tr_yyyyz_xz[i] * tbe_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_xyy[i] = -2.0 * tr_yyyyz_yy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_xyz[i] = -2.0 * tr_yyyyz_yz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_xzz[i] = -2.0 * tr_yyyyz_zz[i] * tbe_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_yyy[i] = 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_yyz[i] = 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_yzz[i] = 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyy_zzz[i] = 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1010-1020 components of targeted buffer : GF

    auto tr_z_0_x_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1010);

    auto tr_z_0_x_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1011);

    auto tr_z_0_x_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1012);

    auto tr_z_0_x_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1013);

    auto tr_z_0_x_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1014);

    auto tr_z_0_x_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1015);

    auto tr_z_0_x_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1016);

    auto tr_z_0_x_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1017);

    auto tr_z_0_x_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1018);

    auto tr_z_0_x_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1019);

    #pragma omp simd aligned(tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_yyy_xx, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyzz_xx, tr_yyyzz_xxxx, tr_yyyzz_xxxy, tr_yyyzz_xxxz, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xy, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xz, tr_yyyzz_xzzz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_z_0_x_yyyz_xxx, tr_z_0_x_yyyz_xxy, tr_z_0_x_yyyz_xxz, tr_z_0_x_yyyz_xyy, tr_z_0_x_yyyz_xyz, tr_z_0_x_yyyz_xzz, tr_z_0_x_yyyz_yyy, tr_z_0_x_yyyz_yyz, tr_z_0_x_yyyz_yzz, tr_z_0_x_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyyz_xxx[i] = 3.0 * tr_yyy_xx[i] - 2.0 * tr_yyy_xxxx[i] * tke_0 - 6.0 * tr_yyyzz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxx[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_xxy[i] = 2.0 * tr_yyy_xy[i] - 2.0 * tr_yyy_xxxy[i] * tke_0 - 4.0 * tr_yyyzz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_xxz[i] = 2.0 * tr_yyy_xz[i] - 2.0 * tr_yyy_xxxz[i] * tke_0 - 4.0 * tr_yyyzz_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_xyy[i] = tr_yyy_yy[i] - 2.0 * tr_yyy_xxyy[i] * tke_0 - 2.0 * tr_yyyzz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_xyz[i] = tr_yyy_yz[i] - 2.0 * tr_yyy_xxyz[i] * tke_0 - 2.0 * tr_yyyzz_yz[i] * tbe_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_xzz[i] = tr_yyy_zz[i] - 2.0 * tr_yyy_xxzz[i] * tke_0 - 2.0 * tr_yyyzz_zz[i] * tbe_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_yyy[i] = -2.0 * tr_yyy_xyyy[i] * tke_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_yyz[i] = -2.0 * tr_yyy_xyyz[i] * tke_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_yzz[i] = -2.0 * tr_yyy_xyzz[i] * tke_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyyz_zzz[i] = -2.0 * tr_yyy_xzzz[i] * tke_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1020-1030 components of targeted buffer : GF

    auto tr_z_0_x_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1020);

    auto tr_z_0_x_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1021);

    auto tr_z_0_x_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1022);

    auto tr_z_0_x_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1023);

    auto tr_z_0_x_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1024);

    auto tr_z_0_x_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1025);

    auto tr_z_0_x_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1026);

    auto tr_z_0_x_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1027);

    auto tr_z_0_x_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1028);

    auto tr_z_0_x_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1029);

    #pragma omp simd aligned(tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_yyz_xx, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yyzzz_xx, tr_yyzzz_xxxx, tr_yyzzz_xxxy, tr_yyzzz_xxxz, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xy, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xz, tr_yyzzz_xzzz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_z_0_x_yyzz_xxx, tr_z_0_x_yyzz_xxy, tr_z_0_x_yyzz_xxz, tr_z_0_x_yyzz_xyy, tr_z_0_x_yyzz_xyz, tr_z_0_x_yyzz_xzz, tr_z_0_x_yyzz_yyy, tr_z_0_x_yyzz_yyz, tr_z_0_x_yyzz_yzz, tr_z_0_x_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyzz_xxx[i] = 6.0 * tr_yyz_xx[i] - 4.0 * tr_yyz_xxxx[i] * tke_0 - 6.0 * tr_yyzzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_xxy[i] = 4.0 * tr_yyz_xy[i] - 4.0 * tr_yyz_xxxy[i] * tke_0 - 4.0 * tr_yyzzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_xxz[i] = 4.0 * tr_yyz_xz[i] - 4.0 * tr_yyz_xxxz[i] * tke_0 - 4.0 * tr_yyzzz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_xyy[i] = 2.0 * tr_yyz_yy[i] - 4.0 * tr_yyz_xxyy[i] * tke_0 - 2.0 * tr_yyzzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_xyz[i] = 2.0 * tr_yyz_yz[i] - 4.0 * tr_yyz_xxyz[i] * tke_0 - 2.0 * tr_yyzzz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_xzz[i] = 2.0 * tr_yyz_zz[i] - 4.0 * tr_yyz_xxzz[i] * tke_0 - 2.0 * tr_yyzzz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_yyy[i] = -4.0 * tr_yyz_xyyy[i] * tke_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_yyz[i] = -4.0 * tr_yyz_xyyz[i] * tke_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_yzz[i] = -4.0 * tr_yyz_xyzz[i] * tke_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyzz_zzz[i] = -4.0 * tr_yyz_xzzz[i] * tke_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1030-1040 components of targeted buffer : GF

    auto tr_z_0_x_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1030);

    auto tr_z_0_x_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1031);

    auto tr_z_0_x_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1032);

    auto tr_z_0_x_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1033);

    auto tr_z_0_x_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1034);

    auto tr_z_0_x_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1035);

    auto tr_z_0_x_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1036);

    auto tr_z_0_x_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1037);

    auto tr_z_0_x_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1038);

    auto tr_z_0_x_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1039);

    #pragma omp simd aligned(tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, tr_yzz_xx, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_yzzzz_xx, tr_yzzzz_xxxx, tr_yzzzz_xxxy, tr_yzzzz_xxxz, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xy, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xz, tr_yzzzz_xzzz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, tr_z_0_x_yzzz_xxx, tr_z_0_x_yzzz_xxy, tr_z_0_x_yzzz_xxz, tr_z_0_x_yzzz_xyy, tr_z_0_x_yzzz_xyz, tr_z_0_x_yzzz_xzz, tr_z_0_x_yzzz_yyy, tr_z_0_x_yzzz_yyz, tr_z_0_x_yzzz_yzz, tr_z_0_x_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yzzz_xxx[i] = 9.0 * tr_yzz_xx[i] - 6.0 * tr_yzz_xxxx[i] * tke_0 - 6.0 * tr_yzzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_xxy[i] = 6.0 * tr_yzz_xy[i] - 6.0 * tr_yzz_xxxy[i] * tke_0 - 4.0 * tr_yzzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_xxz[i] = 6.0 * tr_yzz_xz[i] - 6.0 * tr_yzz_xxxz[i] * tke_0 - 4.0 * tr_yzzzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_xyy[i] = 3.0 * tr_yzz_yy[i] - 6.0 * tr_yzz_xxyy[i] * tke_0 - 2.0 * tr_yzzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_xyz[i] = 3.0 * tr_yzz_yz[i] - 6.0 * tr_yzz_xxyz[i] * tke_0 - 2.0 * tr_yzzzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_xzz[i] = 3.0 * tr_yzz_zz[i] - 6.0 * tr_yzz_xxzz[i] * tke_0 - 2.0 * tr_yzzzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_yyy[i] = -6.0 * tr_yzz_xyyy[i] * tke_0 + 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_yyz[i] = -6.0 * tr_yzz_xyyz[i] * tke_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_yzz[i] = -6.0 * tr_yzz_xyzz[i] * tke_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzzz_zzz[i] = -6.0 * tr_yzz_xzzz[i] * tke_0 + 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1040-1050 components of targeted buffer : GF

    auto tr_z_0_x_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1040);

    auto tr_z_0_x_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1041);

    auto tr_z_0_x_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1042);

    auto tr_z_0_x_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1043);

    auto tr_z_0_x_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1044);

    auto tr_z_0_x_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1045);

    auto tr_z_0_x_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1046);

    auto tr_z_0_x_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1047);

    auto tr_z_0_x_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1048);

    auto tr_z_0_x_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1049);

    #pragma omp simd aligned(tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_xzzzzz_xxx, tr_xzzzzz_xxy, tr_xzzzzz_xxz, tr_xzzzzz_xyy, tr_xzzzzz_xyz, tr_xzzzzz_xzz, tr_xzzzzz_yyy, tr_xzzzzz_yyz, tr_xzzzzz_yzz, tr_xzzzzz_zzz, tr_z_0_x_zzzz_xxx, tr_z_0_x_zzzz_xxy, tr_z_0_x_zzzz_xxz, tr_z_0_x_zzzz_xyy, tr_z_0_x_zzzz_xyz, tr_z_0_x_zzzz_xzz, tr_z_0_x_zzzz_yyy, tr_z_0_x_zzzz_yyz, tr_z_0_x_zzzz_yzz, tr_z_0_x_zzzz_zzz, tr_zzz_xx, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, tr_zzzzz_xx, tr_zzzzz_xxxx, tr_zzzzz_xxxy, tr_zzzzz_xxxz, tr_zzzzz_xxyy, tr_zzzzz_xxyz, tr_zzzzz_xxzz, tr_zzzzz_xy, tr_zzzzz_xyyy, tr_zzzzz_xyyz, tr_zzzzz_xyzz, tr_zzzzz_xz, tr_zzzzz_xzzz, tr_zzzzz_yy, tr_zzzzz_yz, tr_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zzzz_xxx[i] = 12.0 * tr_zzz_xx[i] - 8.0 * tr_zzz_xxxx[i] * tke_0 - 6.0 * tr_zzzzz_xx[i] * tbe_0 + 4.0 * tr_zzzzz_xxxx[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_xxy[i] = 8.0 * tr_zzz_xy[i] - 8.0 * tr_zzz_xxxy[i] * tke_0 - 4.0 * tr_zzzzz_xy[i] * tbe_0 + 4.0 * tr_zzzzz_xxxy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_xxz[i] = 8.0 * tr_zzz_xz[i] - 8.0 * tr_zzz_xxxz[i] * tke_0 - 4.0 * tr_zzzzz_xz[i] * tbe_0 + 4.0 * tr_zzzzz_xxxz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xxz[i] * tbe_0 + 4.0 * tr_xzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_xyy[i] = 4.0 * tr_zzz_yy[i] - 8.0 * tr_zzz_xxyy[i] * tke_0 - 2.0 * tr_zzzzz_yy[i] * tbe_0 + 4.0 * tr_zzzzz_xxyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_xyz[i] = 4.0 * tr_zzz_yz[i] - 8.0 * tr_zzz_xxyz[i] * tke_0 - 2.0 * tr_zzzzz_yz[i] * tbe_0 + 4.0 * tr_zzzzz_xxyz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xyz[i] * tbe_0 + 4.0 * tr_xzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_xzz[i] = 4.0 * tr_zzz_zz[i] - 8.0 * tr_zzz_xxzz[i] * tke_0 - 2.0 * tr_zzzzz_zz[i] * tbe_0 + 4.0 * tr_zzzzz_xxzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_xzz[i] * tbe_0 + 4.0 * tr_xzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_yyy[i] = -8.0 * tr_zzz_xyyy[i] * tke_0 + 4.0 * tr_zzzzz_xyyy[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_yyz[i] = -8.0 * tr_zzz_xyyz[i] * tke_0 + 4.0 * tr_zzzzz_xyyz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_yzz[i] = -8.0 * tr_zzz_xyzz[i] * tke_0 + 4.0 * tr_zzzzz_xyzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_yzz[i] * tbe_0 + 4.0 * tr_xzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzzz_zzz[i] = -8.0 * tr_zzz_xzzz[i] * tke_0 + 4.0 * tr_zzzzz_xzzz[i] * tbe_0 * tke_0 - 8.0 * tr_xzzz_zzz[i] * tbe_0 + 4.0 * tr_xzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1050-1060 components of targeted buffer : GF

    auto tr_z_0_y_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 1050);

    auto tr_z_0_y_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 1051);

    auto tr_z_0_y_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 1052);

    auto tr_z_0_y_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 1053);

    auto tr_z_0_y_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 1054);

    auto tr_z_0_y_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 1055);

    auto tr_z_0_y_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 1056);

    auto tr_z_0_y_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 1057);

    auto tr_z_0_y_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 1058);

    auto tr_z_0_y_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 1059);

    #pragma omp simd aligned(tr_xxxxyz_xxx, tr_xxxxyz_xxy, tr_xxxxyz_xxz, tr_xxxxyz_xyy, tr_xxxxyz_xyz, tr_xxxxyz_xzz, tr_xxxxyz_yyy, tr_xxxxyz_yyz, tr_xxxxyz_yzz, tr_xxxxyz_zzz, tr_xxxxz_xx, tr_xxxxz_xxxy, tr_xxxxz_xxyy, tr_xxxxz_xxyz, tr_xxxxz_xy, tr_xxxxz_xyyy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yyyy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yz, tr_xxxxz_yzzz, tr_xxxxz_zz, tr_z_0_y_xxxx_xxx, tr_z_0_y_xxxx_xxy, tr_z_0_y_xxxx_xxz, tr_z_0_y_xxxx_xyy, tr_z_0_y_xxxx_xyz, tr_z_0_y_xxxx_xzz, tr_z_0_y_xxxx_yyy, tr_z_0_y_xxxx_yyz, tr_z_0_y_xxxx_yzz, tr_z_0_y_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxxx_xxx[i] = 4.0 * tr_xxxxz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_xxy[i] = -2.0 * tr_xxxxz_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_xxz[i] = 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_xyy[i] = -4.0 * tr_xxxxz_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_xyz[i] = -2.0 * tr_xxxxz_xz[i] * tbe_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_xzz[i] = 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_yyy[i] = -6.0 * tr_xxxxz_yy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_yyz[i] = -4.0 * tr_xxxxz_yz[i] * tbe_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_yzz[i] = -2.0 * tr_xxxxz_zz[i] * tbe_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxx_zzz[i] = 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1060-1070 components of targeted buffer : GF

    auto tr_z_0_y_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 1060);

    auto tr_z_0_y_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 1061);

    auto tr_z_0_y_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 1062);

    auto tr_z_0_y_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 1063);

    auto tr_z_0_y_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 1064);

    auto tr_z_0_y_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 1065);

    auto tr_z_0_y_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 1066);

    auto tr_z_0_y_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 1067);

    auto tr_z_0_y_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 1068);

    auto tr_z_0_y_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 1069);

    #pragma omp simd aligned(tr_xxxyyz_xxx, tr_xxxyyz_xxy, tr_xxxyyz_xxz, tr_xxxyyz_xyy, tr_xxxyyz_xyz, tr_xxxyyz_xzz, tr_xxxyyz_yyy, tr_xxxyyz_yyz, tr_xxxyyz_yzz, tr_xxxyyz_zzz, tr_xxxyz_xx, tr_xxxyz_xxxy, tr_xxxyz_xxyy, tr_xxxyz_xxyz, tr_xxxyz_xy, tr_xxxyz_xyyy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yyyy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yz, tr_xxxyz_yzzz, tr_xxxyz_zz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_z_0_y_xxxy_xxx, tr_z_0_y_xxxy_xxy, tr_z_0_y_xxxy_xxz, tr_z_0_y_xxxy_xyy, tr_z_0_y_xxxy_xyz, tr_z_0_y_xxxy_xzz, tr_z_0_y_xxxy_yyy, tr_z_0_y_xxxy_yyz, tr_z_0_y_xxxy_yzz, tr_z_0_y_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxxy_xxx[i] = -2.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_xxy[i] = -2.0 * tr_xxxz_xxy[i] * tbe_0 - 2.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_xxz[i] = -2.0 * tr_xxxz_xxz[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_xyy[i] = -2.0 * tr_xxxz_xyy[i] * tbe_0 - 4.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_xyz[i] = -2.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_xzz[i] = -2.0 * tr_xxxz_xzz[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_yyy[i] = -2.0 * tr_xxxz_yyy[i] * tbe_0 - 6.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_yyz[i] = -2.0 * tr_xxxz_yyz[i] * tbe_0 - 4.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_yzz[i] = -2.0 * tr_xxxz_yzz[i] * tbe_0 - 2.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxy_zzz[i] = -2.0 * tr_xxxz_zzz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1070-1080 components of targeted buffer : GF

    auto tr_z_0_y_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 1070);

    auto tr_z_0_y_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 1071);

    auto tr_z_0_y_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 1072);

    auto tr_z_0_y_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 1073);

    auto tr_z_0_y_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 1074);

    auto tr_z_0_y_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 1075);

    auto tr_z_0_y_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 1076);

    auto tr_z_0_y_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 1077);

    auto tr_z_0_y_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 1078);

    auto tr_z_0_y_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 1079);

    #pragma omp simd aligned(tr_xxx_xx, tr_xxx_xxxy, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_yy, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_xxxzz_xx, tr_xxxzz_xxxy, tr_xxxzz_xxyy, tr_xxxzz_xxyz, tr_xxxzz_xy, tr_xxxzz_xyyy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yyyy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yz, tr_xxxzz_yzzz, tr_xxxzz_zz, tr_z_0_y_xxxz_xxx, tr_z_0_y_xxxz_xxy, tr_z_0_y_xxxz_xxz, tr_z_0_y_xxxz_xyy, tr_z_0_y_xxxz_xyz, tr_z_0_y_xxxz_xzz, tr_z_0_y_xxxz_yyy, tr_z_0_y_xxxz_yyz, tr_z_0_y_xxxz_yzz, tr_z_0_y_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxxz_xxx[i] = -2.0 * tr_xxx_xxxy[i] * tke_0 + 4.0 * tr_xxxzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_xxy[i] = tr_xxx_xx[i] - 2.0 * tr_xxx_xxyy[i] * tke_0 - 2.0 * tr_xxxzz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_xxz[i] = -2.0 * tr_xxx_xxyz[i] * tke_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xxz[i] * tbe_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_xyy[i] = 2.0 * tr_xxx_xy[i] - 2.0 * tr_xxx_xyyy[i] * tke_0 - 4.0 * tr_xxxzz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_xyz[i] = tr_xxx_xz[i] - 2.0 * tr_xxx_xyyz[i] * tke_0 - 2.0 * tr_xxxzz_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xyz[i] * tbe_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_xzz[i] = -2.0 * tr_xxx_xyzz[i] * tke_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_xzz[i] * tbe_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_yyy[i] = 3.0 * tr_xxx_yy[i] - 2.0 * tr_xxx_yyyy[i] * tke_0 - 6.0 * tr_xxxzz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_yyz[i] = 2.0 * tr_xxx_yz[i] - 2.0 * tr_xxx_yyyz[i] * tke_0 - 4.0 * tr_xxxzz_yz[i] * tbe_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yyz[i] * tbe_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_yzz[i] = tr_xxx_zz[i] - 2.0 * tr_xxx_yyzz[i] * tke_0 - 2.0 * tr_xxxzz_zz[i] * tbe_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_yzz[i] * tbe_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxxz_zzz[i] = -2.0 * tr_xxx_yzzz[i] * tke_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxxy_zzz[i] * tbe_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1080-1090 components of targeted buffer : GF

    auto tr_z_0_y_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1080);

    auto tr_z_0_y_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1081);

    auto tr_z_0_y_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1082);

    auto tr_z_0_y_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1083);

    auto tr_z_0_y_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1084);

    auto tr_z_0_y_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1085);

    auto tr_z_0_y_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1086);

    auto tr_z_0_y_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1087);

    auto tr_z_0_y_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1088);

    auto tr_z_0_y_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1089);

    #pragma omp simd aligned(tr_xxyyyz_xxx, tr_xxyyyz_xxy, tr_xxyyyz_xxz, tr_xxyyyz_xyy, tr_xxyyyz_xyz, tr_xxyyyz_xzz, tr_xxyyyz_yyy, tr_xxyyyz_yyz, tr_xxyyyz_yzz, tr_xxyyyz_zzz, tr_xxyyz_xx, tr_xxyyz_xxxy, tr_xxyyz_xxyy, tr_xxyyz_xxyz, tr_xxyyz_xy, tr_xxyyz_xyyy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yyyy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yz, tr_xxyyz_yzzz, tr_xxyyz_zz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_z_0_y_xxyy_xxx, tr_z_0_y_xxyy_xxy, tr_z_0_y_xxyy_xxz, tr_z_0_y_xxyy_xyy, tr_z_0_y_xxyy_xyz, tr_z_0_y_xxyy_xzz, tr_z_0_y_xxyy_yyy, tr_z_0_y_xxyy_yyz, tr_z_0_y_xxyy_yzz, tr_z_0_y_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxyy_xxx[i] = -4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_xxy[i] = -4.0 * tr_xxyz_xxy[i] * tbe_0 - 2.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_xxz[i] = -4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_xyy[i] = -4.0 * tr_xxyz_xyy[i] * tbe_0 - 4.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_xyz[i] = -4.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_xzz[i] = -4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_yyy[i] = -4.0 * tr_xxyz_yyy[i] * tbe_0 - 6.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_yyz[i] = -4.0 * tr_xxyz_yyz[i] * tbe_0 - 4.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_yzz[i] = -4.0 * tr_xxyz_yzz[i] * tbe_0 - 2.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyy_zzz[i] = -4.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1090-1100 components of targeted buffer : GF

    auto tr_z_0_y_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1090);

    auto tr_z_0_y_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1091);

    auto tr_z_0_y_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1092);

    auto tr_z_0_y_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1093);

    auto tr_z_0_y_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1094);

    auto tr_z_0_y_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1095);

    auto tr_z_0_y_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1096);

    auto tr_z_0_y_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1097);

    auto tr_z_0_y_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1098);

    auto tr_z_0_y_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1099);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxy_xx, tr_xxy_xxxy, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_yy, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxy, tr_xxyzz_xxyy, tr_xxyzz_xxyz, tr_xxyzz_xy, tr_xxyzz_xyyy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yyyy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yz, tr_xxyzz_yzzz, tr_xxyzz_zz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_z_0_y_xxyz_xxx, tr_z_0_y_xxyz_xxy, tr_z_0_y_xxyz_xxz, tr_z_0_y_xxyz_xyy, tr_z_0_y_xxyz_xyz, tr_z_0_y_xxyz_xzz, tr_z_0_y_xxyz_yyy, tr_z_0_y_xxyz_yyz, tr_z_0_y_xxyz_yzz, tr_z_0_y_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxyz_xxx[i] = tr_xx_xxx[i] - 2.0 * tr_xxzz_xxx[i] * tbe_0 - 2.0 * tr_xxy_xxxy[i] * tke_0 + 4.0 * tr_xxyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_xxy[i] = tr_xx_xxy[i] - 2.0 * tr_xxzz_xxy[i] * tbe_0 + tr_xxy_xx[i] - 2.0 * tr_xxy_xxyy[i] * tke_0 - 2.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_xxz[i] = tr_xx_xxz[i] - 2.0 * tr_xxzz_xxz[i] * tbe_0 - 2.0 * tr_xxy_xxyz[i] * tke_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xxz[i] * tbe_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_xyy[i] = tr_xx_xyy[i] - 2.0 * tr_xxzz_xyy[i] * tbe_0 + 2.0 * tr_xxy_xy[i] - 2.0 * tr_xxy_xyyy[i] * tke_0 - 4.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_xyz[i] = tr_xx_xyz[i] - 2.0 * tr_xxzz_xyz[i] * tbe_0 + tr_xxy_xz[i] - 2.0 * tr_xxy_xyyz[i] * tke_0 - 2.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xyz[i] * tbe_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_xzz[i] = tr_xx_xzz[i] - 2.0 * tr_xxzz_xzz[i] * tbe_0 - 2.0 * tr_xxy_xyzz[i] * tke_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_xzz[i] * tbe_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_yyy[i] = tr_xx_yyy[i] - 2.0 * tr_xxzz_yyy[i] * tbe_0 + 3.0 * tr_xxy_yy[i] - 2.0 * tr_xxy_yyyy[i] * tke_0 - 6.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_yyz[i] = tr_xx_yyz[i] - 2.0 * tr_xxzz_yyz[i] * tbe_0 + 2.0 * tr_xxy_yz[i] - 2.0 * tr_xxy_yyyz[i] * tke_0 - 4.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yyz[i] * tbe_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_yzz[i] = tr_xx_yzz[i] - 2.0 * tr_xxzz_yzz[i] * tbe_0 + tr_xxy_zz[i] - 2.0 * tr_xxy_yyzz[i] * tke_0 - 2.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_yzz[i] * tbe_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxyz_zzz[i] = tr_xx_zzz[i] - 2.0 * tr_xxzz_zzz[i] * tbe_0 - 2.0 * tr_xxy_yzzz[i] * tke_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxyy_zzz[i] * tbe_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1100-1110 components of targeted buffer : GF

    auto tr_z_0_y_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1100);

    auto tr_z_0_y_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1101);

    auto tr_z_0_y_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1102);

    auto tr_z_0_y_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1103);

    auto tr_z_0_y_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1104);

    auto tr_z_0_y_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1105);

    auto tr_z_0_y_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1106);

    auto tr_z_0_y_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1107);

    auto tr_z_0_y_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1108);

    auto tr_z_0_y_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1109);

    #pragma omp simd aligned(tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_xxz_xx, tr_xxz_xxxy, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_yy, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxzzz_xx, tr_xxzzz_xxxy, tr_xxzzz_xxyy, tr_xxzzz_xxyz, tr_xxzzz_xy, tr_xxzzz_xyyy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yyyy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yz, tr_xxzzz_yzzz, tr_xxzzz_zz, tr_z_0_y_xxzz_xxx, tr_z_0_y_xxzz_xxy, tr_z_0_y_xxzz_xxz, tr_z_0_y_xxzz_xyy, tr_z_0_y_xxzz_xyz, tr_z_0_y_xxzz_xzz, tr_z_0_y_xxzz_yyy, tr_z_0_y_xxzz_yyz, tr_z_0_y_xxzz_yzz, tr_z_0_y_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxzz_xxx[i] = -4.0 * tr_xxz_xxxy[i] * tke_0 + 4.0 * tr_xxzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_xxy[i] = 2.0 * tr_xxz_xx[i] - 4.0 * tr_xxz_xxyy[i] * tke_0 - 2.0 * tr_xxzzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_xxz[i] = -4.0 * tr_xxz_xxyz[i] * tke_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xxz[i] * tbe_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_xyy[i] = 4.0 * tr_xxz_xy[i] - 4.0 * tr_xxz_xyyy[i] * tke_0 - 4.0 * tr_xxzzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_xyz[i] = 2.0 * tr_xxz_xz[i] - 4.0 * tr_xxz_xyyz[i] * tke_0 - 2.0 * tr_xxzzz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xyz[i] * tbe_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_xzz[i] = -4.0 * tr_xxz_xyzz[i] * tke_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_xzz[i] * tbe_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_yyy[i] = 6.0 * tr_xxz_yy[i] - 4.0 * tr_xxz_yyyy[i] * tke_0 - 6.0 * tr_xxzzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_yyz[i] = 4.0 * tr_xxz_yz[i] - 4.0 * tr_xxz_yyyz[i] * tke_0 - 4.0 * tr_xxzzz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yyz[i] * tbe_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_yzz[i] = 2.0 * tr_xxz_zz[i] - 4.0 * tr_xxz_yyzz[i] * tke_0 - 2.0 * tr_xxzzz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_yzz[i] * tbe_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxzz_zzz[i] = -4.0 * tr_xxz_yzzz[i] * tke_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxyz_zzz[i] * tbe_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1110-1120 components of targeted buffer : GF

    auto tr_z_0_y_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1110);

    auto tr_z_0_y_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1111);

    auto tr_z_0_y_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1112);

    auto tr_z_0_y_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1113);

    auto tr_z_0_y_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1114);

    auto tr_z_0_y_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1115);

    auto tr_z_0_y_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1116);

    auto tr_z_0_y_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1117);

    auto tr_z_0_y_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1118);

    auto tr_z_0_y_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1119);

    #pragma omp simd aligned(tr_xyyyyz_xxx, tr_xyyyyz_xxy, tr_xyyyyz_xxz, tr_xyyyyz_xyy, tr_xyyyyz_xyz, tr_xyyyyz_xzz, tr_xyyyyz_yyy, tr_xyyyyz_yyz, tr_xyyyyz_yzz, tr_xyyyyz_zzz, tr_xyyyz_xx, tr_xyyyz_xxxy, tr_xyyyz_xxyy, tr_xyyyz_xxyz, tr_xyyyz_xy, tr_xyyyz_xyyy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yyyy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yz, tr_xyyyz_yzzz, tr_xyyyz_zz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_z_0_y_xyyy_xxx, tr_z_0_y_xyyy_xxy, tr_z_0_y_xyyy_xxz, tr_z_0_y_xyyy_xyy, tr_z_0_y_xyyy_xyz, tr_z_0_y_xyyy_xzz, tr_z_0_y_xyyy_yyy, tr_z_0_y_xyyy_yyz, tr_z_0_y_xyyy_yzz, tr_z_0_y_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyyy_xxx[i] = -6.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_xxy[i] = -6.0 * tr_xyyz_xxy[i] * tbe_0 - 2.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_xxz[i] = -6.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_xyy[i] = -6.0 * tr_xyyz_xyy[i] * tbe_0 - 4.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_xyz[i] = -6.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_xzz[i] = -6.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_yyy[i] = -6.0 * tr_xyyz_yyy[i] * tbe_0 - 6.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_yyz[i] = -6.0 * tr_xyyz_yyz[i] * tbe_0 - 4.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_yzz[i] = -6.0 * tr_xyyz_yzz[i] * tbe_0 - 2.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyy_zzz[i] = -6.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1120-1130 components of targeted buffer : GF

    auto tr_z_0_y_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1120);

    auto tr_z_0_y_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1121);

    auto tr_z_0_y_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1122);

    auto tr_z_0_y_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1123);

    auto tr_z_0_y_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1124);

    auto tr_z_0_y_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1125);

    auto tr_z_0_y_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1126);

    auto tr_z_0_y_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1127);

    auto tr_z_0_y_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1128);

    auto tr_z_0_y_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1129);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxy, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_yy, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxy, tr_xyyzz_xxyy, tr_xyyzz_xxyz, tr_xyyzz_xy, tr_xyyzz_xyyy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yyyy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yz, tr_xyyzz_yzzz, tr_xyyzz_zz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_z_0_y_xyyz_xxx, tr_z_0_y_xyyz_xxy, tr_z_0_y_xyyz_xxz, tr_z_0_y_xyyz_xyy, tr_z_0_y_xyyz_xyz, tr_z_0_y_xyyz_xzz, tr_z_0_y_xyyz_yyy, tr_z_0_y_xyyz_yyz, tr_z_0_y_xyyz_yzz, tr_z_0_y_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyyz_xxx[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xyzz_xxx[i] * tbe_0 - 2.0 * tr_xyy_xxxy[i] * tke_0 + 4.0 * tr_xyyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_xxy[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xyzz_xxy[i] * tbe_0 + tr_xyy_xx[i] - 2.0 * tr_xyy_xxyy[i] * tke_0 - 2.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_xxz[i] = 2.0 * tr_xy_xxz[i] - 4.0 * tr_xyzz_xxz[i] * tbe_0 - 2.0 * tr_xyy_xxyz[i] * tke_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xxz[i] * tbe_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_xyy[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xyzz_xyy[i] * tbe_0 + 2.0 * tr_xyy_xy[i] - 2.0 * tr_xyy_xyyy[i] * tke_0 - 4.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_xyz[i] = 2.0 * tr_xy_xyz[i] - 4.0 * tr_xyzz_xyz[i] * tbe_0 + tr_xyy_xz[i] - 2.0 * tr_xyy_xyyz[i] * tke_0 - 2.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xyz[i] * tbe_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_xzz[i] = 2.0 * tr_xy_xzz[i] - 4.0 * tr_xyzz_xzz[i] * tbe_0 - 2.0 * tr_xyy_xyzz[i] * tke_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_xzz[i] * tbe_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_yyy[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xyzz_yyy[i] * tbe_0 + 3.0 * tr_xyy_yy[i] - 2.0 * tr_xyy_yyyy[i] * tke_0 - 6.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_yyz[i] = 2.0 * tr_xy_yyz[i] - 4.0 * tr_xyzz_yyz[i] * tbe_0 + 2.0 * tr_xyy_yz[i] - 2.0 * tr_xyy_yyyz[i] * tke_0 - 4.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yyz[i] * tbe_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_yzz[i] = 2.0 * tr_xy_yzz[i] - 4.0 * tr_xyzz_yzz[i] * tbe_0 + tr_xyy_zz[i] - 2.0 * tr_xyy_yyzz[i] * tke_0 - 2.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_yzz[i] * tbe_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyyz_zzz[i] = 2.0 * tr_xy_zzz[i] - 4.0 * tr_xyzz_zzz[i] * tbe_0 - 2.0 * tr_xyy_yzzz[i] * tke_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyyy_zzz[i] * tbe_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1130-1140 components of targeted buffer : GF

    auto tr_z_0_y_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1130);

    auto tr_z_0_y_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1131);

    auto tr_z_0_y_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1132);

    auto tr_z_0_y_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1133);

    auto tr_z_0_y_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1134);

    auto tr_z_0_y_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1135);

    auto tr_z_0_y_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1136);

    auto tr_z_0_y_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1137);

    auto tr_z_0_y_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1138);

    auto tr_z_0_y_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1139);

    #pragma omp simd aligned(tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyzzz_xx, tr_xyzzz_xxxy, tr_xyzzz_xxyy, tr_xyzzz_xxyz, tr_xyzzz_xy, tr_xyzzz_xyyy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yyyy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yz, tr_xyzzz_yzzz, tr_xyzzz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_z_0_y_xyzz_xxx, tr_z_0_y_xyzz_xxy, tr_z_0_y_xyzz_xxz, tr_z_0_y_xyzz_xyy, tr_z_0_y_xyzz_xyz, tr_z_0_y_xyzz_xzz, tr_z_0_y_xyzz_yyy, tr_z_0_y_xyzz_yyz, tr_z_0_y_xyzz_yzz, tr_z_0_y_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyzz_xxx[i] = 2.0 * tr_xz_xxx[i] - 2.0 * tr_xzzz_xxx[i] * tbe_0 - 4.0 * tr_xyz_xxxy[i] * tke_0 + 4.0 * tr_xyzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_xxy[i] = 2.0 * tr_xz_xxy[i] - 2.0 * tr_xzzz_xxy[i] * tbe_0 + 2.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxyy[i] * tke_0 - 2.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_xxz[i] = 2.0 * tr_xz_xxz[i] - 2.0 * tr_xzzz_xxz[i] * tbe_0 - 4.0 * tr_xyz_xxyz[i] * tke_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xxz[i] * tbe_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_xyy[i] = 2.0 * tr_xz_xyy[i] - 2.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xyyy[i] * tke_0 - 4.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_xyz[i] = 2.0 * tr_xz_xyz[i] - 2.0 * tr_xzzz_xyz[i] * tbe_0 + 2.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 2.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xyz[i] * tbe_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_xzz[i] = 2.0 * tr_xz_xzz[i] - 2.0 * tr_xzzz_xzz[i] * tbe_0 - 4.0 * tr_xyz_xyzz[i] * tke_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_xzz[i] * tbe_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_yyy[i] = 2.0 * tr_xz_yyy[i] - 2.0 * tr_xzzz_yyy[i] * tbe_0 + 6.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_yyyy[i] * tke_0 - 6.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_yyz[i] = 2.0 * tr_xz_yyz[i] - 2.0 * tr_xzzz_yyz[i] * tbe_0 + 4.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_yyyz[i] * tke_0 - 4.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yyz[i] * tbe_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_yzz[i] = 2.0 * tr_xz_yzz[i] - 2.0 * tr_xzzz_yzz[i] * tbe_0 + 2.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_yyzz[i] * tke_0 - 2.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_yzz[i] * tbe_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyzz_zzz[i] = 2.0 * tr_xz_zzz[i] - 2.0 * tr_xzzz_zzz[i] * tbe_0 - 4.0 * tr_xyz_yzzz[i] * tke_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyyz_zzz[i] * tbe_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1140-1150 components of targeted buffer : GF

    auto tr_z_0_y_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1140);

    auto tr_z_0_y_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1141);

    auto tr_z_0_y_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1142);

    auto tr_z_0_y_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1143);

    auto tr_z_0_y_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1144);

    auto tr_z_0_y_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1145);

    auto tr_z_0_y_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1146);

    auto tr_z_0_y_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1147);

    auto tr_z_0_y_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1148);

    auto tr_z_0_y_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1149);

    #pragma omp simd aligned(tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, tr_xzz_xx, tr_xzz_xxxy, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_yy, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzzzz_xx, tr_xzzzz_xxxy, tr_xzzzz_xxyy, tr_xzzzz_xxyz, tr_xzzzz_xy, tr_xzzzz_xyyy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yyyy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yz, tr_xzzzz_yzzz, tr_xzzzz_zz, tr_z_0_y_xzzz_xxx, tr_z_0_y_xzzz_xxy, tr_z_0_y_xzzz_xxz, tr_z_0_y_xzzz_xyy, tr_z_0_y_xzzz_xyz, tr_z_0_y_xzzz_xzz, tr_z_0_y_xzzz_yyy, tr_z_0_y_xzzz_yyz, tr_z_0_y_xzzz_yzz, tr_z_0_y_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xzzz_xxx[i] = -6.0 * tr_xzz_xxxy[i] * tke_0 + 4.0 * tr_xzzzz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_xxy[i] = 3.0 * tr_xzz_xx[i] - 6.0 * tr_xzz_xxyy[i] * tke_0 - 2.0 * tr_xzzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_xxz[i] = -6.0 * tr_xzz_xxyz[i] * tke_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xxz[i] * tbe_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_xyy[i] = 6.0 * tr_xzz_xy[i] - 6.0 * tr_xzz_xyyy[i] * tke_0 - 4.0 * tr_xzzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_xyz[i] = 3.0 * tr_xzz_xz[i] - 6.0 * tr_xzz_xyyz[i] * tke_0 - 2.0 * tr_xzzzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xyz[i] * tbe_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_xzz[i] = -6.0 * tr_xzz_xyzz[i] * tke_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_xzz[i] * tbe_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_yyy[i] = 9.0 * tr_xzz_yy[i] - 6.0 * tr_xzz_yyyy[i] * tke_0 - 6.0 * tr_xzzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_yyz[i] = 6.0 * tr_xzz_yz[i] - 6.0 * tr_xzz_yyyz[i] * tke_0 - 4.0 * tr_xzzzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yyz[i] * tbe_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_yzz[i] = 3.0 * tr_xzz_zz[i] - 6.0 * tr_xzz_yyzz[i] * tke_0 - 2.0 * tr_xzzzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_yzz[i] * tbe_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzzz_zzz[i] = -6.0 * tr_xzz_yzzz[i] * tke_0 + 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyzz_zzz[i] * tbe_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1150-1160 components of targeted buffer : GF

    auto tr_z_0_y_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1150);

    auto tr_z_0_y_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1151);

    auto tr_z_0_y_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1152);

    auto tr_z_0_y_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1153);

    auto tr_z_0_y_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1154);

    auto tr_z_0_y_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1155);

    auto tr_z_0_y_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1156);

    auto tr_z_0_y_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1157);

    auto tr_z_0_y_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1158);

    auto tr_z_0_y_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1159);

    #pragma omp simd aligned(tr_yyyyyz_xxx, tr_yyyyyz_xxy, tr_yyyyyz_xxz, tr_yyyyyz_xyy, tr_yyyyyz_xyz, tr_yyyyyz_xzz, tr_yyyyyz_yyy, tr_yyyyyz_yyz, tr_yyyyyz_yzz, tr_yyyyyz_zzz, tr_yyyyz_xx, tr_yyyyz_xxxy, tr_yyyyz_xxyy, tr_yyyyz_xxyz, tr_yyyyz_xy, tr_yyyyz_xyyy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yyyy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yz, tr_yyyyz_yzzz, tr_yyyyz_zz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_z_0_y_yyyy_xxx, tr_z_0_y_yyyy_xxy, tr_z_0_y_yyyy_xxz, tr_z_0_y_yyyy_xyy, tr_z_0_y_yyyy_xyz, tr_z_0_y_yyyy_xzz, tr_z_0_y_yyyy_yyy, tr_z_0_y_yyyy_yyz, tr_z_0_y_yyyy_yzz, tr_z_0_y_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyyy_xxx[i] = -8.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_xxy[i] = -8.0 * tr_yyyz_xxy[i] * tbe_0 - 2.0 * tr_yyyyz_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_xxz[i] = -8.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_xyy[i] = -8.0 * tr_yyyz_xyy[i] * tbe_0 - 4.0 * tr_yyyyz_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_xyz[i] = -8.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_yyyyz_xz[i] * tbe_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_xzz[i] = -8.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_yyy[i] = -8.0 * tr_yyyz_yyy[i] * tbe_0 - 6.0 * tr_yyyyz_yy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_yyz[i] = -8.0 * tr_yyyz_yyz[i] * tbe_0 - 4.0 * tr_yyyyz_yz[i] * tbe_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_yzz[i] = -8.0 * tr_yyyz_yzz[i] * tbe_0 - 2.0 * tr_yyyyz_zz[i] * tbe_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyy_zzz[i] = -8.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyyz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1160-1170 components of targeted buffer : GF

    auto tr_z_0_y_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1160);

    auto tr_z_0_y_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1161);

    auto tr_z_0_y_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1162);

    auto tr_z_0_y_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1163);

    auto tr_z_0_y_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1164);

    auto tr_z_0_y_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1165);

    auto tr_z_0_y_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1166);

    auto tr_z_0_y_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1167);

    auto tr_z_0_y_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1168);

    auto tr_z_0_y_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1169);

    #pragma omp simd aligned(tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyy_xx, tr_yyy_xxxy, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_yy, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, tr_yyyyzz_xxx, tr_yyyyzz_xxy, tr_yyyyzz_xxz, tr_yyyyzz_xyy, tr_yyyyzz_xyz, tr_yyyyzz_xzz, tr_yyyyzz_yyy, tr_yyyyzz_yyz, tr_yyyyzz_yzz, tr_yyyyzz_zzz, tr_yyyzz_xx, tr_yyyzz_xxxy, tr_yyyzz_xxyy, tr_yyyzz_xxyz, tr_yyyzz_xy, tr_yyyzz_xyyy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yyyy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yz, tr_yyyzz_yzzz, tr_yyyzz_zz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_z_0_y_yyyz_xxx, tr_z_0_y_yyyz_xxy, tr_z_0_y_yyyz_xxz, tr_z_0_y_yyyz_xyy, tr_z_0_y_yyyz_xyz, tr_z_0_y_yyyz_xzz, tr_z_0_y_yyyz_yyy, tr_z_0_y_yyyz_yyz, tr_z_0_y_yyyz_yzz, tr_z_0_y_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyyz_xxx[i] = 3.0 * tr_yy_xxx[i] - 6.0 * tr_yyzz_xxx[i] * tbe_0 - 2.0 * tr_yyy_xxxy[i] * tke_0 + 4.0 * tr_yyyzz_xxxy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_xxy[i] = 3.0 * tr_yy_xxy[i] - 6.0 * tr_yyzz_xxy[i] * tbe_0 + tr_yyy_xx[i] - 2.0 * tr_yyy_xxyy[i] * tke_0 - 2.0 * tr_yyyzz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_xxz[i] = 3.0 * tr_yy_xxz[i] - 6.0 * tr_yyzz_xxz[i] * tbe_0 - 2.0 * tr_yyy_xxyz[i] * tke_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xxz[i] * tbe_0 + 4.0 * tr_yyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_xyy[i] = 3.0 * tr_yy_xyy[i] - 6.0 * tr_yyzz_xyy[i] * tbe_0 + 2.0 * tr_yyy_xy[i] - 2.0 * tr_yyy_xyyy[i] * tke_0 - 4.0 * tr_yyyzz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_xyz[i] = 3.0 * tr_yy_xyz[i] - 6.0 * tr_yyzz_xyz[i] * tbe_0 + tr_yyy_xz[i] - 2.0 * tr_yyy_xyyz[i] * tke_0 - 2.0 * tr_yyyzz_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xyz[i] * tbe_0 + 4.0 * tr_yyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_xzz[i] = 3.0 * tr_yy_xzz[i] - 6.0 * tr_yyzz_xzz[i] * tbe_0 - 2.0 * tr_yyy_xyzz[i] * tke_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_xzz[i] * tbe_0 + 4.0 * tr_yyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_yyy[i] = 3.0 * tr_yy_yyy[i] - 6.0 * tr_yyzz_yyy[i] * tbe_0 + 3.0 * tr_yyy_yy[i] - 2.0 * tr_yyy_yyyy[i] * tke_0 - 6.0 * tr_yyyzz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_yyz[i] = 3.0 * tr_yy_yyz[i] - 6.0 * tr_yyzz_yyz[i] * tbe_0 + 2.0 * tr_yyy_yz[i] - 2.0 * tr_yyy_yyyz[i] * tke_0 - 4.0 * tr_yyyzz_yz[i] * tbe_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yyz[i] * tbe_0 + 4.0 * tr_yyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_yzz[i] = 3.0 * tr_yy_yzz[i] - 6.0 * tr_yyzz_yzz[i] * tbe_0 + tr_yyy_zz[i] - 2.0 * tr_yyy_yyzz[i] * tke_0 - 2.0 * tr_yyyzz_zz[i] * tbe_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_yzz[i] * tbe_0 + 4.0 * tr_yyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyyz_zzz[i] = 3.0 * tr_yy_zzz[i] - 6.0 * tr_yyzz_zzz[i] * tbe_0 - 2.0 * tr_yyy_yzzz[i] * tke_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyyy_zzz[i] * tbe_0 + 4.0 * tr_yyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1170-1180 components of targeted buffer : GF

    auto tr_z_0_y_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1170);

    auto tr_z_0_y_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1171);

    auto tr_z_0_y_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1172);

    auto tr_z_0_y_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1173);

    auto tr_z_0_y_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1174);

    auto tr_z_0_y_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1175);

    auto tr_z_0_y_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1176);

    auto tr_z_0_y_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1177);

    auto tr_z_0_y_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1178);

    auto tr_z_0_y_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1179);

    #pragma omp simd aligned(tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyyzzz_xxx, tr_yyyzzz_xxy, tr_yyyzzz_xxz, tr_yyyzzz_xyy, tr_yyyzzz_xyz, tr_yyyzzz_xzz, tr_yyyzzz_yyy, tr_yyyzzz_yyz, tr_yyyzzz_yzz, tr_yyyzzz_zzz, tr_yyz_xx, tr_yyz_xxxy, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_yy, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyzzz_xx, tr_yyzzz_xxxy, tr_yyzzz_xxyy, tr_yyzzz_xxyz, tr_yyzzz_xy, tr_yyzzz_xyyy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yyyy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yz, tr_yyzzz_yzzz, tr_yyzzz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_z_0_y_yyzz_xxx, tr_z_0_y_yyzz_xxy, tr_z_0_y_yyzz_xxz, tr_z_0_y_yyzz_xyy, tr_z_0_y_yyzz_xyz, tr_z_0_y_yyzz_xzz, tr_z_0_y_yyzz_yyy, tr_z_0_y_yyzz_yyz, tr_z_0_y_yyzz_yzz, tr_z_0_y_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyzz_xxx[i] = 4.0 * tr_yz_xxx[i] - 4.0 * tr_yzzz_xxx[i] * tbe_0 - 4.0 * tr_yyz_xxxy[i] * tke_0 + 4.0 * tr_yyzzz_xxxy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_xxy[i] = 4.0 * tr_yz_xxy[i] - 4.0 * tr_yzzz_xxy[i] * tbe_0 + 2.0 * tr_yyz_xx[i] - 4.0 * tr_yyz_xxyy[i] * tke_0 - 2.0 * tr_yyzzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_xxz[i] = 4.0 * tr_yz_xxz[i] - 4.0 * tr_yzzz_xxz[i] * tbe_0 - 4.0 * tr_yyz_xxyz[i] * tke_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xxz[i] * tbe_0 + 4.0 * tr_yyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_xyy[i] = 4.0 * tr_yz_xyy[i] - 4.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yyz_xy[i] - 4.0 * tr_yyz_xyyy[i] * tke_0 - 4.0 * tr_yyzzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_xyz[i] = 4.0 * tr_yz_xyz[i] - 4.0 * tr_yzzz_xyz[i] * tbe_0 + 2.0 * tr_yyz_xz[i] - 4.0 * tr_yyz_xyyz[i] * tke_0 - 2.0 * tr_yyzzz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xyz[i] * tbe_0 + 4.0 * tr_yyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_xzz[i] = 4.0 * tr_yz_xzz[i] - 4.0 * tr_yzzz_xzz[i] * tbe_0 - 4.0 * tr_yyz_xyzz[i] * tke_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_xzz[i] * tbe_0 + 4.0 * tr_yyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_yyy[i] = 4.0 * tr_yz_yyy[i] - 4.0 * tr_yzzz_yyy[i] * tbe_0 + 6.0 * tr_yyz_yy[i] - 4.0 * tr_yyz_yyyy[i] * tke_0 - 6.0 * tr_yyzzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_yyz[i] = 4.0 * tr_yz_yyz[i] - 4.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yyz_yz[i] - 4.0 * tr_yyz_yyyz[i] * tke_0 - 4.0 * tr_yyzzz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yyz[i] * tbe_0 + 4.0 * tr_yyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_yzz[i] = 4.0 * tr_yz_yzz[i] - 4.0 * tr_yzzz_yzz[i] * tbe_0 + 2.0 * tr_yyz_zz[i] - 4.0 * tr_yyz_yyzz[i] * tke_0 - 2.0 * tr_yyzzz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_yzz[i] * tbe_0 + 4.0 * tr_yyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyzz_zzz[i] = 4.0 * tr_yz_zzz[i] - 4.0 * tr_yzzz_zzz[i] * tbe_0 - 4.0 * tr_yyz_yzzz[i] * tke_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyyz_zzz[i] * tbe_0 + 4.0 * tr_yyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1180-1190 components of targeted buffer : GF

    auto tr_z_0_y_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1180);

    auto tr_z_0_y_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1181);

    auto tr_z_0_y_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1182);

    auto tr_z_0_y_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1183);

    auto tr_z_0_y_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1184);

    auto tr_z_0_y_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1185);

    auto tr_z_0_y_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1186);

    auto tr_z_0_y_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1187);

    auto tr_z_0_y_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1188);

    auto tr_z_0_y_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1189);

    #pragma omp simd aligned(tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yyzzzz_xxx, tr_yyzzzz_xxy, tr_yyzzzz_xxz, tr_yyzzzz_xyy, tr_yyzzzz_xyz, tr_yyzzzz_xzz, tr_yyzzzz_yyy, tr_yyzzzz_yyz, tr_yyzzzz_yzz, tr_yyzzzz_zzz, tr_yzz_xx, tr_yzz_xxxy, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_yy, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzzzz_xx, tr_yzzzz_xxxy, tr_yzzzz_xxyy, tr_yzzzz_xxyz, tr_yzzzz_xy, tr_yzzzz_xyyy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yyyy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yz, tr_yzzzz_yzzz, tr_yzzzz_zz, tr_z_0_y_yzzz_xxx, tr_z_0_y_yzzz_xxy, tr_z_0_y_yzzz_xxz, tr_z_0_y_yzzz_xyy, tr_z_0_y_yzzz_xyz, tr_z_0_y_yzzz_xzz, tr_z_0_y_yzzz_yyy, tr_z_0_y_yzzz_yyz, tr_z_0_y_yzzz_yzz, tr_z_0_y_yzzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yzzz_xxx[i] = 3.0 * tr_zz_xxx[i] - 2.0 * tr_zzzz_xxx[i] * tbe_0 - 6.0 * tr_yzz_xxxy[i] * tke_0 + 4.0 * tr_yzzzz_xxxy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_xxy[i] = 3.0 * tr_zz_xxy[i] - 2.0 * tr_zzzz_xxy[i] * tbe_0 + 3.0 * tr_yzz_xx[i] - 6.0 * tr_yzz_xxyy[i] * tke_0 - 2.0 * tr_yzzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_xxz[i] = 3.0 * tr_zz_xxz[i] - 2.0 * tr_zzzz_xxz[i] * tbe_0 - 6.0 * tr_yzz_xxyz[i] * tke_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xxz[i] * tbe_0 + 4.0 * tr_yyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_xyy[i] = 3.0 * tr_zz_xyy[i] - 2.0 * tr_zzzz_xyy[i] * tbe_0 + 6.0 * tr_yzz_xy[i] - 6.0 * tr_yzz_xyyy[i] * tke_0 - 4.0 * tr_yzzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_xyz[i] = 3.0 * tr_zz_xyz[i] - 2.0 * tr_zzzz_xyz[i] * tbe_0 + 3.0 * tr_yzz_xz[i] - 6.0 * tr_yzz_xyyz[i] * tke_0 - 2.0 * tr_yzzzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xyz[i] * tbe_0 + 4.0 * tr_yyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_xzz[i] = 3.0 * tr_zz_xzz[i] - 2.0 * tr_zzzz_xzz[i] * tbe_0 - 6.0 * tr_yzz_xyzz[i] * tke_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_xzz[i] * tbe_0 + 4.0 * tr_yyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_yyy[i] = 3.0 * tr_zz_yyy[i] - 2.0 * tr_zzzz_yyy[i] * tbe_0 + 9.0 * tr_yzz_yy[i] - 6.0 * tr_yzz_yyyy[i] * tke_0 - 6.0 * tr_yzzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyy[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_yyz[i] = 3.0 * tr_zz_yyz[i] - 2.0 * tr_zzzz_yyz[i] * tbe_0 + 6.0 * tr_yzz_yz[i] - 6.0 * tr_yzz_yyyz[i] * tke_0 - 4.0 * tr_yzzzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yyz[i] * tbe_0 + 4.0 * tr_yyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_yzz[i] = 3.0 * tr_zz_yzz[i] - 2.0 * tr_zzzz_yzz[i] * tbe_0 + 3.0 * tr_yzz_zz[i] - 6.0 * tr_yzz_yyzz[i] * tke_0 - 2.0 * tr_yzzzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_yzz[i] * tbe_0 + 4.0 * tr_yyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzzz_zzz[i] = 3.0 * tr_zz_zzz[i] - 2.0 * tr_zzzz_zzz[i] * tbe_0 - 6.0 * tr_yzz_yzzz[i] * tke_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tke_0 - 6.0 * tr_yyzz_zzz[i] * tbe_0 + 4.0 * tr_yyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1190-1200 components of targeted buffer : GF

    auto tr_z_0_y_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1190);

    auto tr_z_0_y_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1191);

    auto tr_z_0_y_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1192);

    auto tr_z_0_y_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1193);

    auto tr_z_0_y_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1194);

    auto tr_z_0_y_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1195);

    auto tr_z_0_y_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1196);

    auto tr_z_0_y_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1197);

    auto tr_z_0_y_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1198);

    auto tr_z_0_y_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1199);

    #pragma omp simd aligned(tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_yzzzzz_xxx, tr_yzzzzz_xxy, tr_yzzzzz_xxz, tr_yzzzzz_xyy, tr_yzzzzz_xyz, tr_yzzzzz_xzz, tr_yzzzzz_yyy, tr_yzzzzz_yyz, tr_yzzzzz_yzz, tr_yzzzzz_zzz, tr_z_0_y_zzzz_xxx, tr_z_0_y_zzzz_xxy, tr_z_0_y_zzzz_xxz, tr_z_0_y_zzzz_xyy, tr_z_0_y_zzzz_xyz, tr_z_0_y_zzzz_xzz, tr_z_0_y_zzzz_yyy, tr_z_0_y_zzzz_yyz, tr_z_0_y_zzzz_yzz, tr_z_0_y_zzzz_zzz, tr_zzz_xx, tr_zzz_xxxy, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_yy, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, tr_zzzzz_xx, tr_zzzzz_xxxy, tr_zzzzz_xxyy, tr_zzzzz_xxyz, tr_zzzzz_xy, tr_zzzzz_xyyy, tr_zzzzz_xyyz, tr_zzzzz_xyzz, tr_zzzzz_xz, tr_zzzzz_yy, tr_zzzzz_yyyy, tr_zzzzz_yyyz, tr_zzzzz_yyzz, tr_zzzzz_yz, tr_zzzzz_yzzz, tr_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zzzz_xxx[i] = -8.0 * tr_zzz_xxxy[i] * tke_0 + 4.0 * tr_zzzzz_xxxy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_xxy[i] = 4.0 * tr_zzz_xx[i] - 8.0 * tr_zzz_xxyy[i] * tke_0 - 2.0 * tr_zzzzz_xx[i] * tbe_0 + 4.0 * tr_zzzzz_xxyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_xxz[i] = -8.0 * tr_zzz_xxyz[i] * tke_0 + 4.0 * tr_zzzzz_xxyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xxz[i] * tbe_0 + 4.0 * tr_yzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_xyy[i] = 8.0 * tr_zzz_xy[i] - 8.0 * tr_zzz_xyyy[i] * tke_0 - 4.0 * tr_zzzzz_xy[i] * tbe_0 + 4.0 * tr_zzzzz_xyyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_xyz[i] = 4.0 * tr_zzz_xz[i] - 8.0 * tr_zzz_xyyz[i] * tke_0 - 2.0 * tr_zzzzz_xz[i] * tbe_0 + 4.0 * tr_zzzzz_xyyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xyz[i] * tbe_0 + 4.0 * tr_yzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_xzz[i] = -8.0 * tr_zzz_xyzz[i] * tke_0 + 4.0 * tr_zzzzz_xyzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_xzz[i] * tbe_0 + 4.0 * tr_yzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_yyy[i] = 12.0 * tr_zzz_yy[i] - 8.0 * tr_zzz_yyyy[i] * tke_0 - 6.0 * tr_zzzzz_yy[i] * tbe_0 + 4.0 * tr_zzzzz_yyyy[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_yyz[i] = 8.0 * tr_zzz_yz[i] - 8.0 * tr_zzz_yyyz[i] * tke_0 - 4.0 * tr_zzzzz_yz[i] * tbe_0 + 4.0 * tr_zzzzz_yyyz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yyz[i] * tbe_0 + 4.0 * tr_yzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_yzz[i] = 4.0 * tr_zzz_zz[i] - 8.0 * tr_zzz_yyzz[i] * tke_0 - 2.0 * tr_zzzzz_zz[i] * tbe_0 + 4.0 * tr_zzzzz_yyzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_yzz[i] * tbe_0 + 4.0 * tr_yzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzzz_zzz[i] = -8.0 * tr_zzz_yzzz[i] * tke_0 + 4.0 * tr_zzzzz_yzzz[i] * tbe_0 * tke_0 - 8.0 * tr_yzzz_zzz[i] * tbe_0 + 4.0 * tr_yzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1200-1210 components of targeted buffer : GF

    auto tr_z_0_z_xxxx_xxx = pbuffer.data(idx_op_geom_110_gf + 1200);

    auto tr_z_0_z_xxxx_xxy = pbuffer.data(idx_op_geom_110_gf + 1201);

    auto tr_z_0_z_xxxx_xxz = pbuffer.data(idx_op_geom_110_gf + 1202);

    auto tr_z_0_z_xxxx_xyy = pbuffer.data(idx_op_geom_110_gf + 1203);

    auto tr_z_0_z_xxxx_xyz = pbuffer.data(idx_op_geom_110_gf + 1204);

    auto tr_z_0_z_xxxx_xzz = pbuffer.data(idx_op_geom_110_gf + 1205);

    auto tr_z_0_z_xxxx_yyy = pbuffer.data(idx_op_geom_110_gf + 1206);

    auto tr_z_0_z_xxxx_yyz = pbuffer.data(idx_op_geom_110_gf + 1207);

    auto tr_z_0_z_xxxx_yzz = pbuffer.data(idx_op_geom_110_gf + 1208);

    auto tr_z_0_z_xxxx_zzz = pbuffer.data(idx_op_geom_110_gf + 1209);

    #pragma omp simd aligned(tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, tr_xxxxz_xx, tr_xxxxz_xxxz, tr_xxxxz_xxyz, tr_xxxxz_xxzz, tr_xxxxz_xy, tr_xxxxz_xyyz, tr_xxxxz_xyzz, tr_xxxxz_xz, tr_xxxxz_xzzz, tr_xxxxz_yy, tr_xxxxz_yyyz, tr_xxxxz_yyzz, tr_xxxxz_yz, tr_xxxxz_yzzz, tr_xxxxz_zz, tr_xxxxz_zzzz, tr_xxxxzz_xxx, tr_xxxxzz_xxy, tr_xxxxzz_xxz, tr_xxxxzz_xyy, tr_xxxxzz_xyz, tr_xxxxzz_xzz, tr_xxxxzz_yyy, tr_xxxxzz_yyz, tr_xxxxzz_yzz, tr_xxxxzz_zzz, tr_z_0_z_xxxx_xxx, tr_z_0_z_xxxx_xxy, tr_z_0_z_xxxx_xxz, tr_z_0_z_xxxx_xyy, tr_z_0_z_xxxx_xyz, tr_z_0_z_xxxx_xzz, tr_z_0_z_xxxx_yyy, tr_z_0_z_xxxx_yyz, tr_z_0_z_xxxx_yzz, tr_z_0_z_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxxx_xxx[i] = -2.0 * tr_xxxx_xxx[i] * tbe_0 + 4.0 * tr_xxxxz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_xxy[i] = -2.0 * tr_xxxx_xxy[i] * tbe_0 + 4.0 * tr_xxxxz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_xxz[i] = -2.0 * tr_xxxx_xxz[i] * tbe_0 - 2.0 * tr_xxxxz_xx[i] * tbe_0 + 4.0 * tr_xxxxz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_xyy[i] = -2.0 * tr_xxxx_xyy[i] * tbe_0 + 4.0 * tr_xxxxz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_xyz[i] = -2.0 * tr_xxxx_xyz[i] * tbe_0 - 2.0 * tr_xxxxz_xy[i] * tbe_0 + 4.0 * tr_xxxxz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_xzz[i] = -2.0 * tr_xxxx_xzz[i] * tbe_0 - 4.0 * tr_xxxxz_xz[i] * tbe_0 + 4.0 * tr_xxxxz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_yyy[i] = -2.0 * tr_xxxx_yyy[i] * tbe_0 + 4.0 * tr_xxxxz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_yyz[i] = -2.0 * tr_xxxx_yyz[i] * tbe_0 - 2.0 * tr_xxxxz_yy[i] * tbe_0 + 4.0 * tr_xxxxz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_yzz[i] = -2.0 * tr_xxxx_yzz[i] * tbe_0 - 4.0 * tr_xxxxz_yz[i] * tbe_0 + 4.0 * tr_xxxxz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxx_zzz[i] = -2.0 * tr_xxxx_zzz[i] * tbe_0 - 6.0 * tr_xxxxz_zz[i] * tbe_0 + 4.0 * tr_xxxxz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1210-1220 components of targeted buffer : GF

    auto tr_z_0_z_xxxy_xxx = pbuffer.data(idx_op_geom_110_gf + 1210);

    auto tr_z_0_z_xxxy_xxy = pbuffer.data(idx_op_geom_110_gf + 1211);

    auto tr_z_0_z_xxxy_xxz = pbuffer.data(idx_op_geom_110_gf + 1212);

    auto tr_z_0_z_xxxy_xyy = pbuffer.data(idx_op_geom_110_gf + 1213);

    auto tr_z_0_z_xxxy_xyz = pbuffer.data(idx_op_geom_110_gf + 1214);

    auto tr_z_0_z_xxxy_xzz = pbuffer.data(idx_op_geom_110_gf + 1215);

    auto tr_z_0_z_xxxy_yyy = pbuffer.data(idx_op_geom_110_gf + 1216);

    auto tr_z_0_z_xxxy_yyz = pbuffer.data(idx_op_geom_110_gf + 1217);

    auto tr_z_0_z_xxxy_yzz = pbuffer.data(idx_op_geom_110_gf + 1218);

    auto tr_z_0_z_xxxy_zzz = pbuffer.data(idx_op_geom_110_gf + 1219);

    #pragma omp simd aligned(tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxxyz_xx, tr_xxxyz_xxxz, tr_xxxyz_xxyz, tr_xxxyz_xxzz, tr_xxxyz_xy, tr_xxxyz_xyyz, tr_xxxyz_xyzz, tr_xxxyz_xz, tr_xxxyz_xzzz, tr_xxxyz_yy, tr_xxxyz_yyyz, tr_xxxyz_yyzz, tr_xxxyz_yz, tr_xxxyz_yzzz, tr_xxxyz_zz, tr_xxxyz_zzzz, tr_xxxyzz_xxx, tr_xxxyzz_xxy, tr_xxxyzz_xxz, tr_xxxyzz_xyy, tr_xxxyzz_xyz, tr_xxxyzz_xzz, tr_xxxyzz_yyy, tr_xxxyzz_yyz, tr_xxxyzz_yzz, tr_xxxyzz_zzz, tr_z_0_z_xxxy_xxx, tr_z_0_z_xxxy_xxy, tr_z_0_z_xxxy_xxz, tr_z_0_z_xxxy_xyy, tr_z_0_z_xxxy_xyz, tr_z_0_z_xxxy_xzz, tr_z_0_z_xxxy_yyy, tr_z_0_z_xxxy_yyz, tr_z_0_z_xxxy_yzz, tr_z_0_z_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxxy_xxx[i] = -2.0 * tr_xxxy_xxx[i] * tbe_0 + 4.0 * tr_xxxyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_xxy[i] = -2.0 * tr_xxxy_xxy[i] * tbe_0 + 4.0 * tr_xxxyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_xxz[i] = -2.0 * tr_xxxy_xxz[i] * tbe_0 - 2.0 * tr_xxxyz_xx[i] * tbe_0 + 4.0 * tr_xxxyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_xyy[i] = -2.0 * tr_xxxy_xyy[i] * tbe_0 + 4.0 * tr_xxxyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_xyz[i] = -2.0 * tr_xxxy_xyz[i] * tbe_0 - 2.0 * tr_xxxyz_xy[i] * tbe_0 + 4.0 * tr_xxxyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_xzz[i] = -2.0 * tr_xxxy_xzz[i] * tbe_0 - 4.0 * tr_xxxyz_xz[i] * tbe_0 + 4.0 * tr_xxxyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_yyy[i] = -2.0 * tr_xxxy_yyy[i] * tbe_0 + 4.0 * tr_xxxyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_yyz[i] = -2.0 * tr_xxxy_yyz[i] * tbe_0 - 2.0 * tr_xxxyz_yy[i] * tbe_0 + 4.0 * tr_xxxyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_yzz[i] = -2.0 * tr_xxxy_yzz[i] * tbe_0 - 4.0 * tr_xxxyz_yz[i] * tbe_0 + 4.0 * tr_xxxyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxy_zzz[i] = -2.0 * tr_xxxy_zzz[i] * tbe_0 - 6.0 * tr_xxxyz_zz[i] * tbe_0 + 4.0 * tr_xxxyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1220-1230 components of targeted buffer : GF

    auto tr_z_0_z_xxxz_xxx = pbuffer.data(idx_op_geom_110_gf + 1220);

    auto tr_z_0_z_xxxz_xxy = pbuffer.data(idx_op_geom_110_gf + 1221);

    auto tr_z_0_z_xxxz_xxz = pbuffer.data(idx_op_geom_110_gf + 1222);

    auto tr_z_0_z_xxxz_xyy = pbuffer.data(idx_op_geom_110_gf + 1223);

    auto tr_z_0_z_xxxz_xyz = pbuffer.data(idx_op_geom_110_gf + 1224);

    auto tr_z_0_z_xxxz_xzz = pbuffer.data(idx_op_geom_110_gf + 1225);

    auto tr_z_0_z_xxxz_yyy = pbuffer.data(idx_op_geom_110_gf + 1226);

    auto tr_z_0_z_xxxz_yyz = pbuffer.data(idx_op_geom_110_gf + 1227);

    auto tr_z_0_z_xxxz_yzz = pbuffer.data(idx_op_geom_110_gf + 1228);

    auto tr_z_0_z_xxxz_zzz = pbuffer.data(idx_op_geom_110_gf + 1229);

    #pragma omp simd aligned(tr_xxx_xx, tr_xxx_xxxz, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxx_zzzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxxzz_xx, tr_xxxzz_xxxz, tr_xxxzz_xxyz, tr_xxxzz_xxzz, tr_xxxzz_xy, tr_xxxzz_xyyz, tr_xxxzz_xyzz, tr_xxxzz_xz, tr_xxxzz_xzzz, tr_xxxzz_yy, tr_xxxzz_yyyz, tr_xxxzz_yyzz, tr_xxxzz_yz, tr_xxxzz_yzzz, tr_xxxzz_zz, tr_xxxzz_zzzz, tr_xxxzzz_xxx, tr_xxxzzz_xxy, tr_xxxzzz_xxz, tr_xxxzzz_xyy, tr_xxxzzz_xyz, tr_xxxzzz_xzz, tr_xxxzzz_yyy, tr_xxxzzz_yyz, tr_xxxzzz_yzz, tr_xxxzzz_zzz, tr_z_0_z_xxxz_xxx, tr_z_0_z_xxxz_xxy, tr_z_0_z_xxxz_xxz, tr_z_0_z_xxxz_xyy, tr_z_0_z_xxxz_xyz, tr_z_0_z_xxxz_xzz, tr_z_0_z_xxxz_yyy, tr_z_0_z_xxxz_yyz, tr_z_0_z_xxxz_yzz, tr_z_0_z_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxxz_xxx[i] = -2.0 * tr_xxx_xxxz[i] * tke_0 - 6.0 * tr_xxxz_xxx[i] * tbe_0 + 4.0 * tr_xxxzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_xxy[i] = -2.0 * tr_xxx_xxyz[i] * tke_0 - 6.0 * tr_xxxz_xxy[i] * tbe_0 + 4.0 * tr_xxxzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_xxz[i] = tr_xxx_xx[i] - 2.0 * tr_xxx_xxzz[i] * tke_0 - 6.0 * tr_xxxz_xxz[i] * tbe_0 - 2.0 * tr_xxxzz_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_xyy[i] = -2.0 * tr_xxx_xyyz[i] * tke_0 - 6.0 * tr_xxxz_xyy[i] * tbe_0 + 4.0 * tr_xxxzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_xyz[i] = tr_xxx_xy[i] - 2.0 * tr_xxx_xyzz[i] * tke_0 - 6.0 * tr_xxxz_xyz[i] * tbe_0 - 2.0 * tr_xxxzz_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_xzz[i] = 2.0 * tr_xxx_xz[i] - 2.0 * tr_xxx_xzzz[i] * tke_0 - 6.0 * tr_xxxz_xzz[i] * tbe_0 - 4.0 * tr_xxxzz_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_yyy[i] = -2.0 * tr_xxx_yyyz[i] * tke_0 - 6.0 * tr_xxxz_yyy[i] * tbe_0 + 4.0 * tr_xxxzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_yyz[i] = tr_xxx_yy[i] - 2.0 * tr_xxx_yyzz[i] * tke_0 - 6.0 * tr_xxxz_yyz[i] * tbe_0 - 2.0 * tr_xxxzz_yy[i] * tbe_0 + 4.0 * tr_xxxzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_yzz[i] = 2.0 * tr_xxx_yz[i] - 2.0 * tr_xxx_yzzz[i] * tke_0 - 6.0 * tr_xxxz_yzz[i] * tbe_0 - 4.0 * tr_xxxzz_yz[i] * tbe_0 + 4.0 * tr_xxxzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxxz_zzz[i] = 3.0 * tr_xxx_zz[i] - 2.0 * tr_xxx_zzzz[i] * tke_0 - 6.0 * tr_xxxz_zzz[i] * tbe_0 - 6.0 * tr_xxxzz_zz[i] * tbe_0 + 4.0 * tr_xxxzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1230-1240 components of targeted buffer : GF

    auto tr_z_0_z_xxyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1230);

    auto tr_z_0_z_xxyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1231);

    auto tr_z_0_z_xxyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1232);

    auto tr_z_0_z_xxyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1233);

    auto tr_z_0_z_xxyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1234);

    auto tr_z_0_z_xxyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1235);

    auto tr_z_0_z_xxyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1236);

    auto tr_z_0_z_xxyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1237);

    auto tr_z_0_z_xxyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1238);

    auto tr_z_0_z_xxyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1239);

    #pragma omp simd aligned(tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xxyyz_xx, tr_xxyyz_xxxz, tr_xxyyz_xxyz, tr_xxyyz_xxzz, tr_xxyyz_xy, tr_xxyyz_xyyz, tr_xxyyz_xyzz, tr_xxyyz_xz, tr_xxyyz_xzzz, tr_xxyyz_yy, tr_xxyyz_yyyz, tr_xxyyz_yyzz, tr_xxyyz_yz, tr_xxyyz_yzzz, tr_xxyyz_zz, tr_xxyyz_zzzz, tr_xxyyzz_xxx, tr_xxyyzz_xxy, tr_xxyyzz_xxz, tr_xxyyzz_xyy, tr_xxyyzz_xyz, tr_xxyyzz_xzz, tr_xxyyzz_yyy, tr_xxyyzz_yyz, tr_xxyyzz_yzz, tr_xxyyzz_zzz, tr_z_0_z_xxyy_xxx, tr_z_0_z_xxyy_xxy, tr_z_0_z_xxyy_xxz, tr_z_0_z_xxyy_xyy, tr_z_0_z_xxyy_xyz, tr_z_0_z_xxyy_xzz, tr_z_0_z_xxyy_yyy, tr_z_0_z_xxyy_yyz, tr_z_0_z_xxyy_yzz, tr_z_0_z_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxyy_xxx[i] = -2.0 * tr_xxyy_xxx[i] * tbe_0 + 4.0 * tr_xxyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_xxy[i] = -2.0 * tr_xxyy_xxy[i] * tbe_0 + 4.0 * tr_xxyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_xxz[i] = -2.0 * tr_xxyy_xxz[i] * tbe_0 - 2.0 * tr_xxyyz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_xyy[i] = -2.0 * tr_xxyy_xyy[i] * tbe_0 + 4.0 * tr_xxyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_xyz[i] = -2.0 * tr_xxyy_xyz[i] * tbe_0 - 2.0 * tr_xxyyz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_xzz[i] = -2.0 * tr_xxyy_xzz[i] * tbe_0 - 4.0 * tr_xxyyz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_yyy[i] = -2.0 * tr_xxyy_yyy[i] * tbe_0 + 4.0 * tr_xxyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_yyz[i] = -2.0 * tr_xxyy_yyz[i] * tbe_0 - 2.0 * tr_xxyyz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_yzz[i] = -2.0 * tr_xxyy_yzz[i] * tbe_0 - 4.0 * tr_xxyyz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyy_zzz[i] = -2.0 * tr_xxyy_zzz[i] * tbe_0 - 6.0 * tr_xxyyz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1240-1250 components of targeted buffer : GF

    auto tr_z_0_z_xxyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1240);

    auto tr_z_0_z_xxyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1241);

    auto tr_z_0_z_xxyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1242);

    auto tr_z_0_z_xxyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1243);

    auto tr_z_0_z_xxyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1244);

    auto tr_z_0_z_xxyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1245);

    auto tr_z_0_z_xxyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1246);

    auto tr_z_0_z_xxyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1247);

    auto tr_z_0_z_xxyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1248);

    auto tr_z_0_z_xxyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1249);

    #pragma omp simd aligned(tr_xxy_xx, tr_xxy_xxxz, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xxxz, tr_xxyzz_xxyz, tr_xxyzz_xxzz, tr_xxyzz_xy, tr_xxyzz_xyyz, tr_xxyzz_xyzz, tr_xxyzz_xz, tr_xxyzz_xzzz, tr_xxyzz_yy, tr_xxyzz_yyyz, tr_xxyzz_yyzz, tr_xxyzz_yz, tr_xxyzz_yzzz, tr_xxyzz_zz, tr_xxyzz_zzzz, tr_xxyzzz_xxx, tr_xxyzzz_xxy, tr_xxyzzz_xxz, tr_xxyzzz_xyy, tr_xxyzzz_xyz, tr_xxyzzz_xzz, tr_xxyzzz_yyy, tr_xxyzzz_yyz, tr_xxyzzz_yzz, tr_xxyzzz_zzz, tr_z_0_z_xxyz_xxx, tr_z_0_z_xxyz_xxy, tr_z_0_z_xxyz_xxz, tr_z_0_z_xxyz_xyy, tr_z_0_z_xxyz_xyz, tr_z_0_z_xxyz_xzz, tr_z_0_z_xxyz_yyy, tr_z_0_z_xxyz_yyz, tr_z_0_z_xxyz_yzz, tr_z_0_z_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxyz_xxx[i] = -2.0 * tr_xxy_xxxz[i] * tke_0 - 6.0 * tr_xxyz_xxx[i] * tbe_0 + 4.0 * tr_xxyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_xxy[i] = -2.0 * tr_xxy_xxyz[i] * tke_0 - 6.0 * tr_xxyz_xxy[i] * tbe_0 + 4.0 * tr_xxyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_xxz[i] = tr_xxy_xx[i] - 2.0 * tr_xxy_xxzz[i] * tke_0 - 6.0 * tr_xxyz_xxz[i] * tbe_0 - 2.0 * tr_xxyzz_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_xyy[i] = -2.0 * tr_xxy_xyyz[i] * tke_0 - 6.0 * tr_xxyz_xyy[i] * tbe_0 + 4.0 * tr_xxyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_xyz[i] = tr_xxy_xy[i] - 2.0 * tr_xxy_xyzz[i] * tke_0 - 6.0 * tr_xxyz_xyz[i] * tbe_0 - 2.0 * tr_xxyzz_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_xzz[i] = 2.0 * tr_xxy_xz[i] - 2.0 * tr_xxy_xzzz[i] * tke_0 - 6.0 * tr_xxyz_xzz[i] * tbe_0 - 4.0 * tr_xxyzz_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_yyy[i] = -2.0 * tr_xxy_yyyz[i] * tke_0 - 6.0 * tr_xxyz_yyy[i] * tbe_0 + 4.0 * tr_xxyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_yyz[i] = tr_xxy_yy[i] - 2.0 * tr_xxy_yyzz[i] * tke_0 - 6.0 * tr_xxyz_yyz[i] * tbe_0 - 2.0 * tr_xxyzz_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_yzz[i] = 2.0 * tr_xxy_yz[i] - 2.0 * tr_xxy_yzzz[i] * tke_0 - 6.0 * tr_xxyz_yzz[i] * tbe_0 - 4.0 * tr_xxyzz_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxyz_zzz[i] = 3.0 * tr_xxy_zz[i] - 2.0 * tr_xxy_zzzz[i] * tke_0 - 6.0 * tr_xxyz_zzz[i] * tbe_0 - 6.0 * tr_xxyzz_zz[i] * tbe_0 + 4.0 * tr_xxyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1250-1260 components of targeted buffer : GF

    auto tr_z_0_z_xxzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1250);

    auto tr_z_0_z_xxzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1251);

    auto tr_z_0_z_xxzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1252);

    auto tr_z_0_z_xxzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1253);

    auto tr_z_0_z_xxzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1254);

    auto tr_z_0_z_xxzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1255);

    auto tr_z_0_z_xxzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1256);

    auto tr_z_0_z_xxzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1257);

    auto tr_z_0_z_xxzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1258);

    auto tr_z_0_z_xxzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1259);

    #pragma omp simd aligned(tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxz_xx, tr_xxz_xxxz, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xxzzz_xx, tr_xxzzz_xxxz, tr_xxzzz_xxyz, tr_xxzzz_xxzz, tr_xxzzz_xy, tr_xxzzz_xyyz, tr_xxzzz_xyzz, tr_xxzzz_xz, tr_xxzzz_xzzz, tr_xxzzz_yy, tr_xxzzz_yyyz, tr_xxzzz_yyzz, tr_xxzzz_yz, tr_xxzzz_yzzz, tr_xxzzz_zz, tr_xxzzz_zzzz, tr_xxzzzz_xxx, tr_xxzzzz_xxy, tr_xxzzzz_xxz, tr_xxzzzz_xyy, tr_xxzzzz_xyz, tr_xxzzzz_xzz, tr_xxzzzz_yyy, tr_xxzzzz_yyz, tr_xxzzzz_yzz, tr_xxzzzz_zzz, tr_z_0_z_xxzz_xxx, tr_z_0_z_xxzz_xxy, tr_z_0_z_xxzz_xxz, tr_z_0_z_xxzz_xyy, tr_z_0_z_xxzz_xyz, tr_z_0_z_xxzz_xzz, tr_z_0_z_xxzz_yyy, tr_z_0_z_xxzz_yyz, tr_z_0_z_xxzz_yzz, tr_z_0_z_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxzz_xxx[i] = 2.0 * tr_xx_xxx[i] - 4.0 * tr_xxz_xxxz[i] * tke_0 - 10.0 * tr_xxzz_xxx[i] * tbe_0 + 4.0 * tr_xxzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_xxy[i] = 2.0 * tr_xx_xxy[i] - 4.0 * tr_xxz_xxyz[i] * tke_0 - 10.0 * tr_xxzz_xxy[i] * tbe_0 + 4.0 * tr_xxzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_xxz[i] = 2.0 * tr_xx_xxz[i] + 2.0 * tr_xxz_xx[i] - 4.0 * tr_xxz_xxzz[i] * tke_0 - 10.0 * tr_xxzz_xxz[i] * tbe_0 - 2.0 * tr_xxzzz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_xyy[i] = 2.0 * tr_xx_xyy[i] - 4.0 * tr_xxz_xyyz[i] * tke_0 - 10.0 * tr_xxzz_xyy[i] * tbe_0 + 4.0 * tr_xxzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_xyz[i] = 2.0 * tr_xx_xyz[i] + 2.0 * tr_xxz_xy[i] - 4.0 * tr_xxz_xyzz[i] * tke_0 - 10.0 * tr_xxzz_xyz[i] * tbe_0 - 2.0 * tr_xxzzz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_xzz[i] = 2.0 * tr_xx_xzz[i] + 4.0 * tr_xxz_xz[i] - 4.0 * tr_xxz_xzzz[i] * tke_0 - 10.0 * tr_xxzz_xzz[i] * tbe_0 - 4.0 * tr_xxzzz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_yyy[i] = 2.0 * tr_xx_yyy[i] - 4.0 * tr_xxz_yyyz[i] * tke_0 - 10.0 * tr_xxzz_yyy[i] * tbe_0 + 4.0 * tr_xxzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_yyz[i] = 2.0 * tr_xx_yyz[i] + 2.0 * tr_xxz_yy[i] - 4.0 * tr_xxz_yyzz[i] * tke_0 - 10.0 * tr_xxzz_yyz[i] * tbe_0 - 2.0 * tr_xxzzz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_yzz[i] = 2.0 * tr_xx_yzz[i] + 4.0 * tr_xxz_yz[i] - 4.0 * tr_xxz_yzzz[i] * tke_0 - 10.0 * tr_xxzz_yzz[i] * tbe_0 - 4.0 * tr_xxzzz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxzz_zzz[i] = 2.0 * tr_xx_zzz[i] + 6.0 * tr_xxz_zz[i] - 4.0 * tr_xxz_zzzz[i] * tke_0 - 10.0 * tr_xxzz_zzz[i] * tbe_0 - 6.0 * tr_xxzzz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1260-1270 components of targeted buffer : GF

    auto tr_z_0_z_xyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1260);

    auto tr_z_0_z_xyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1261);

    auto tr_z_0_z_xyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1262);

    auto tr_z_0_z_xyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1263);

    auto tr_z_0_z_xyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1264);

    auto tr_z_0_z_xyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1265);

    auto tr_z_0_z_xyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1266);

    auto tr_z_0_z_xyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1267);

    auto tr_z_0_z_xyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1268);

    auto tr_z_0_z_xyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1269);

    #pragma omp simd aligned(tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_xyyyz_xx, tr_xyyyz_xxxz, tr_xyyyz_xxyz, tr_xyyyz_xxzz, tr_xyyyz_xy, tr_xyyyz_xyyz, tr_xyyyz_xyzz, tr_xyyyz_xz, tr_xyyyz_xzzz, tr_xyyyz_yy, tr_xyyyz_yyyz, tr_xyyyz_yyzz, tr_xyyyz_yz, tr_xyyyz_yzzz, tr_xyyyz_zz, tr_xyyyz_zzzz, tr_xyyyzz_xxx, tr_xyyyzz_xxy, tr_xyyyzz_xxz, tr_xyyyzz_xyy, tr_xyyyzz_xyz, tr_xyyyzz_xzz, tr_xyyyzz_yyy, tr_xyyyzz_yyz, tr_xyyyzz_yzz, tr_xyyyzz_zzz, tr_z_0_z_xyyy_xxx, tr_z_0_z_xyyy_xxy, tr_z_0_z_xyyy_xxz, tr_z_0_z_xyyy_xyy, tr_z_0_z_xyyy_xyz, tr_z_0_z_xyyy_xzz, tr_z_0_z_xyyy_yyy, tr_z_0_z_xyyy_yyz, tr_z_0_z_xyyy_yzz, tr_z_0_z_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyyy_xxx[i] = -2.0 * tr_xyyy_xxx[i] * tbe_0 + 4.0 * tr_xyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_xxy[i] = -2.0 * tr_xyyy_xxy[i] * tbe_0 + 4.0 * tr_xyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_xxz[i] = -2.0 * tr_xyyy_xxz[i] * tbe_0 - 2.0 * tr_xyyyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_xyy[i] = -2.0 * tr_xyyy_xyy[i] * tbe_0 + 4.0 * tr_xyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_xyz[i] = -2.0 * tr_xyyy_xyz[i] * tbe_0 - 2.0 * tr_xyyyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_xzz[i] = -2.0 * tr_xyyy_xzz[i] * tbe_0 - 4.0 * tr_xyyyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_yyy[i] = -2.0 * tr_xyyy_yyy[i] * tbe_0 + 4.0 * tr_xyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_yyz[i] = -2.0 * tr_xyyy_yyz[i] * tbe_0 - 2.0 * tr_xyyyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_yzz[i] = -2.0 * tr_xyyy_yzz[i] * tbe_0 - 4.0 * tr_xyyyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyy_zzz[i] = -2.0 * tr_xyyy_zzz[i] * tbe_0 - 6.0 * tr_xyyyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1270-1280 components of targeted buffer : GF

    auto tr_z_0_z_xyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1270);

    auto tr_z_0_z_xyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1271);

    auto tr_z_0_z_xyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1272);

    auto tr_z_0_z_xyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1273);

    auto tr_z_0_z_xyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1274);

    auto tr_z_0_z_xyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1275);

    auto tr_z_0_z_xyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1276);

    auto tr_z_0_z_xyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1277);

    auto tr_z_0_z_xyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1278);

    auto tr_z_0_z_xyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1279);

    #pragma omp simd aligned(tr_xyy_xx, tr_xyy_xxxz, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xxxz, tr_xyyzz_xxyz, tr_xyyzz_xxzz, tr_xyyzz_xy, tr_xyyzz_xyyz, tr_xyyzz_xyzz, tr_xyyzz_xz, tr_xyyzz_xzzz, tr_xyyzz_yy, tr_xyyzz_yyyz, tr_xyyzz_yyzz, tr_xyyzz_yz, tr_xyyzz_yzzz, tr_xyyzz_zz, tr_xyyzz_zzzz, tr_xyyzzz_xxx, tr_xyyzzz_xxy, tr_xyyzzz_xxz, tr_xyyzzz_xyy, tr_xyyzzz_xyz, tr_xyyzzz_xzz, tr_xyyzzz_yyy, tr_xyyzzz_yyz, tr_xyyzzz_yzz, tr_xyyzzz_zzz, tr_z_0_z_xyyz_xxx, tr_z_0_z_xyyz_xxy, tr_z_0_z_xyyz_xxz, tr_z_0_z_xyyz_xyy, tr_z_0_z_xyyz_xyz, tr_z_0_z_xyyz_xzz, tr_z_0_z_xyyz_yyy, tr_z_0_z_xyyz_yyz, tr_z_0_z_xyyz_yzz, tr_z_0_z_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyyz_xxx[i] = -2.0 * tr_xyy_xxxz[i] * tke_0 - 6.0 * tr_xyyz_xxx[i] * tbe_0 + 4.0 * tr_xyyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_xxy[i] = -2.0 * tr_xyy_xxyz[i] * tke_0 - 6.0 * tr_xyyz_xxy[i] * tbe_0 + 4.0 * tr_xyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_xxz[i] = tr_xyy_xx[i] - 2.0 * tr_xyy_xxzz[i] * tke_0 - 6.0 * tr_xyyz_xxz[i] * tbe_0 - 2.0 * tr_xyyzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_xyy[i] = -2.0 * tr_xyy_xyyz[i] * tke_0 - 6.0 * tr_xyyz_xyy[i] * tbe_0 + 4.0 * tr_xyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_xyz[i] = tr_xyy_xy[i] - 2.0 * tr_xyy_xyzz[i] * tke_0 - 6.0 * tr_xyyz_xyz[i] * tbe_0 - 2.0 * tr_xyyzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_xzz[i] = 2.0 * tr_xyy_xz[i] - 2.0 * tr_xyy_xzzz[i] * tke_0 - 6.0 * tr_xyyz_xzz[i] * tbe_0 - 4.0 * tr_xyyzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_yyy[i] = -2.0 * tr_xyy_yyyz[i] * tke_0 - 6.0 * tr_xyyz_yyy[i] * tbe_0 + 4.0 * tr_xyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_yyz[i] = tr_xyy_yy[i] - 2.0 * tr_xyy_yyzz[i] * tke_0 - 6.0 * tr_xyyz_yyz[i] * tbe_0 - 2.0 * tr_xyyzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_yzz[i] = 2.0 * tr_xyy_yz[i] - 2.0 * tr_xyy_yzzz[i] * tke_0 - 6.0 * tr_xyyz_yzz[i] * tbe_0 - 4.0 * tr_xyyzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyyz_zzz[i] = 3.0 * tr_xyy_zz[i] - 2.0 * tr_xyy_zzzz[i] * tke_0 - 6.0 * tr_xyyz_zzz[i] * tbe_0 - 6.0 * tr_xyyzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1280-1290 components of targeted buffer : GF

    auto tr_z_0_z_xyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1280);

    auto tr_z_0_z_xyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1281);

    auto tr_z_0_z_xyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1282);

    auto tr_z_0_z_xyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1283);

    auto tr_z_0_z_xyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1284);

    auto tr_z_0_z_xyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1285);

    auto tr_z_0_z_xyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1286);

    auto tr_z_0_z_xyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1287);

    auto tr_z_0_z_xyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1288);

    auto tr_z_0_z_xyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1289);

    #pragma omp simd aligned(tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xxxz, tr_xyzzz_xxyz, tr_xyzzz_xxzz, tr_xyzzz_xy, tr_xyzzz_xyyz, tr_xyzzz_xyzz, tr_xyzzz_xz, tr_xyzzz_xzzz, tr_xyzzz_yy, tr_xyzzz_yyyz, tr_xyzzz_yyzz, tr_xyzzz_yz, tr_xyzzz_yzzz, tr_xyzzz_zz, tr_xyzzz_zzzz, tr_xyzzzz_xxx, tr_xyzzzz_xxy, tr_xyzzzz_xxz, tr_xyzzzz_xyy, tr_xyzzzz_xyz, tr_xyzzzz_xzz, tr_xyzzzz_yyy, tr_xyzzzz_yyz, tr_xyzzzz_yzz, tr_xyzzzz_zzz, tr_z_0_z_xyzz_xxx, tr_z_0_z_xyzz_xxy, tr_z_0_z_xyzz_xxz, tr_z_0_z_xyzz_xyy, tr_z_0_z_xyzz_xyz, tr_z_0_z_xyzz_xzz, tr_z_0_z_xyzz_yyy, tr_z_0_z_xyzz_yyz, tr_z_0_z_xyzz_yzz, tr_z_0_z_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyzz_xxx[i] = 2.0 * tr_xy_xxx[i] - 4.0 * tr_xyz_xxxz[i] * tke_0 - 10.0 * tr_xyzz_xxx[i] * tbe_0 + 4.0 * tr_xyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_xxy[i] = 2.0 * tr_xy_xxy[i] - 4.0 * tr_xyz_xxyz[i] * tke_0 - 10.0 * tr_xyzz_xxy[i] * tbe_0 + 4.0 * tr_xyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_xxz[i] = 2.0 * tr_xy_xxz[i] + 2.0 * tr_xyz_xx[i] - 4.0 * tr_xyz_xxzz[i] * tke_0 - 10.0 * tr_xyzz_xxz[i] * tbe_0 - 2.0 * tr_xyzzz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_xyy[i] = 2.0 * tr_xy_xyy[i] - 4.0 * tr_xyz_xyyz[i] * tke_0 - 10.0 * tr_xyzz_xyy[i] * tbe_0 + 4.0 * tr_xyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_xyz[i] = 2.0 * tr_xy_xyz[i] + 2.0 * tr_xyz_xy[i] - 4.0 * tr_xyz_xyzz[i] * tke_0 - 10.0 * tr_xyzz_xyz[i] * tbe_0 - 2.0 * tr_xyzzz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_xzz[i] = 2.0 * tr_xy_xzz[i] + 4.0 * tr_xyz_xz[i] - 4.0 * tr_xyz_xzzz[i] * tke_0 - 10.0 * tr_xyzz_xzz[i] * tbe_0 - 4.0 * tr_xyzzz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_yyy[i] = 2.0 * tr_xy_yyy[i] - 4.0 * tr_xyz_yyyz[i] * tke_0 - 10.0 * tr_xyzz_yyy[i] * tbe_0 + 4.0 * tr_xyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_yyz[i] = 2.0 * tr_xy_yyz[i] + 2.0 * tr_xyz_yy[i] - 4.0 * tr_xyz_yyzz[i] * tke_0 - 10.0 * tr_xyzz_yyz[i] * tbe_0 - 2.0 * tr_xyzzz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_yzz[i] = 2.0 * tr_xy_yzz[i] + 4.0 * tr_xyz_yz[i] - 4.0 * tr_xyz_yzzz[i] * tke_0 - 10.0 * tr_xyzz_yzz[i] * tbe_0 - 4.0 * tr_xyzzz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyzz_zzz[i] = 2.0 * tr_xy_zzz[i] + 6.0 * tr_xyz_zz[i] - 4.0 * tr_xyz_zzzz[i] * tke_0 - 10.0 * tr_xyzz_zzz[i] * tbe_0 - 6.0 * tr_xyzzz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1290-1300 components of targeted buffer : GF

    auto tr_z_0_z_xzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1290);

    auto tr_z_0_z_xzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1291);

    auto tr_z_0_z_xzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1292);

    auto tr_z_0_z_xzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1293);

    auto tr_z_0_z_xzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1294);

    auto tr_z_0_z_xzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1295);

    auto tr_z_0_z_xzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1296);

    auto tr_z_0_z_xzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1297);

    auto tr_z_0_z_xzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1298);

    auto tr_z_0_z_xzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1299);

    #pragma omp simd aligned(tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxz, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xxxz, tr_xzzzz_xxyz, tr_xzzzz_xxzz, tr_xzzzz_xy, tr_xzzzz_xyyz, tr_xzzzz_xyzz, tr_xzzzz_xz, tr_xzzzz_xzzz, tr_xzzzz_yy, tr_xzzzz_yyyz, tr_xzzzz_yyzz, tr_xzzzz_yz, tr_xzzzz_yzzz, tr_xzzzz_zz, tr_xzzzz_zzzz, tr_xzzzzz_xxx, tr_xzzzzz_xxy, tr_xzzzzz_xxz, tr_xzzzzz_xyy, tr_xzzzzz_xyz, tr_xzzzzz_xzz, tr_xzzzzz_yyy, tr_xzzzzz_yyz, tr_xzzzzz_yzz, tr_xzzzzz_zzz, tr_z_0_z_xzzz_xxx, tr_z_0_z_xzzz_xxy, tr_z_0_z_xzzz_xxz, tr_z_0_z_xzzz_xyy, tr_z_0_z_xzzz_xyz, tr_z_0_z_xzzz_xzz, tr_z_0_z_xzzz_yyy, tr_z_0_z_xzzz_yyz, tr_z_0_z_xzzz_yzz, tr_z_0_z_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xzzz_xxx[i] = 6.0 * tr_xz_xxx[i] - 6.0 * tr_xzz_xxxz[i] * tke_0 - 14.0 * tr_xzzz_xxx[i] * tbe_0 + 4.0 * tr_xzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_xxy[i] = 6.0 * tr_xz_xxy[i] - 6.0 * tr_xzz_xxyz[i] * tke_0 - 14.0 * tr_xzzz_xxy[i] * tbe_0 + 4.0 * tr_xzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_xxz[i] = 6.0 * tr_xz_xxz[i] + 3.0 * tr_xzz_xx[i] - 6.0 * tr_xzz_xxzz[i] * tke_0 - 14.0 * tr_xzzz_xxz[i] * tbe_0 - 2.0 * tr_xzzzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_xyy[i] = 6.0 * tr_xz_xyy[i] - 6.0 * tr_xzz_xyyz[i] * tke_0 - 14.0 * tr_xzzz_xyy[i] * tbe_0 + 4.0 * tr_xzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_xyz[i] = 6.0 * tr_xz_xyz[i] + 3.0 * tr_xzz_xy[i] - 6.0 * tr_xzz_xyzz[i] * tke_0 - 14.0 * tr_xzzz_xyz[i] * tbe_0 - 2.0 * tr_xzzzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_xzz[i] = 6.0 * tr_xz_xzz[i] + 6.0 * tr_xzz_xz[i] - 6.0 * tr_xzz_xzzz[i] * tke_0 - 14.0 * tr_xzzz_xzz[i] * tbe_0 - 4.0 * tr_xzzzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_yyy[i] = 6.0 * tr_xz_yyy[i] - 6.0 * tr_xzz_yyyz[i] * tke_0 - 14.0 * tr_xzzz_yyy[i] * tbe_0 + 4.0 * tr_xzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_yyz[i] = 6.0 * tr_xz_yyz[i] + 3.0 * tr_xzz_yy[i] - 6.0 * tr_xzz_yyzz[i] * tke_0 - 14.0 * tr_xzzz_yyz[i] * tbe_0 - 2.0 * tr_xzzzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_yzz[i] = 6.0 * tr_xz_yzz[i] + 6.0 * tr_xzz_yz[i] - 6.0 * tr_xzz_yzzz[i] * tke_0 - 14.0 * tr_xzzz_yzz[i] * tbe_0 - 4.0 * tr_xzzzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzzz_zzz[i] = 6.0 * tr_xz_zzz[i] + 9.0 * tr_xzz_zz[i] - 6.0 * tr_xzz_zzzz[i] * tke_0 - 14.0 * tr_xzzz_zzz[i] * tbe_0 - 6.0 * tr_xzzzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1300-1310 components of targeted buffer : GF

    auto tr_z_0_z_yyyy_xxx = pbuffer.data(idx_op_geom_110_gf + 1300);

    auto tr_z_0_z_yyyy_xxy = pbuffer.data(idx_op_geom_110_gf + 1301);

    auto tr_z_0_z_yyyy_xxz = pbuffer.data(idx_op_geom_110_gf + 1302);

    auto tr_z_0_z_yyyy_xyy = pbuffer.data(idx_op_geom_110_gf + 1303);

    auto tr_z_0_z_yyyy_xyz = pbuffer.data(idx_op_geom_110_gf + 1304);

    auto tr_z_0_z_yyyy_xzz = pbuffer.data(idx_op_geom_110_gf + 1305);

    auto tr_z_0_z_yyyy_yyy = pbuffer.data(idx_op_geom_110_gf + 1306);

    auto tr_z_0_z_yyyy_yyz = pbuffer.data(idx_op_geom_110_gf + 1307);

    auto tr_z_0_z_yyyy_yzz = pbuffer.data(idx_op_geom_110_gf + 1308);

    auto tr_z_0_z_yyyy_zzz = pbuffer.data(idx_op_geom_110_gf + 1309);

    #pragma omp simd aligned(tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, tr_yyyyz_xx, tr_yyyyz_xxxz, tr_yyyyz_xxyz, tr_yyyyz_xxzz, tr_yyyyz_xy, tr_yyyyz_xyyz, tr_yyyyz_xyzz, tr_yyyyz_xz, tr_yyyyz_xzzz, tr_yyyyz_yy, tr_yyyyz_yyyz, tr_yyyyz_yyzz, tr_yyyyz_yz, tr_yyyyz_yzzz, tr_yyyyz_zz, tr_yyyyz_zzzz, tr_yyyyzz_xxx, tr_yyyyzz_xxy, tr_yyyyzz_xxz, tr_yyyyzz_xyy, tr_yyyyzz_xyz, tr_yyyyzz_xzz, tr_yyyyzz_yyy, tr_yyyyzz_yyz, tr_yyyyzz_yzz, tr_yyyyzz_zzz, tr_z_0_z_yyyy_xxx, tr_z_0_z_yyyy_xxy, tr_z_0_z_yyyy_xxz, tr_z_0_z_yyyy_xyy, tr_z_0_z_yyyy_xyz, tr_z_0_z_yyyy_xzz, tr_z_0_z_yyyy_yyy, tr_z_0_z_yyyy_yyz, tr_z_0_z_yyyy_yzz, tr_z_0_z_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyyy_xxx[i] = -2.0 * tr_yyyy_xxx[i] * tbe_0 + 4.0 * tr_yyyyz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_xxy[i] = -2.0 * tr_yyyy_xxy[i] * tbe_0 + 4.0 * tr_yyyyz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_xxz[i] = -2.0 * tr_yyyy_xxz[i] * tbe_0 - 2.0 * tr_yyyyz_xx[i] * tbe_0 + 4.0 * tr_yyyyz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_xyy[i] = -2.0 * tr_yyyy_xyy[i] * tbe_0 + 4.0 * tr_yyyyz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_xyz[i] = -2.0 * tr_yyyy_xyz[i] * tbe_0 - 2.0 * tr_yyyyz_xy[i] * tbe_0 + 4.0 * tr_yyyyz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_xzz[i] = -2.0 * tr_yyyy_xzz[i] * tbe_0 - 4.0 * tr_yyyyz_xz[i] * tbe_0 + 4.0 * tr_yyyyz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_yyy[i] = -2.0 * tr_yyyy_yyy[i] * tbe_0 + 4.0 * tr_yyyyz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_yyz[i] = -2.0 * tr_yyyy_yyz[i] * tbe_0 - 2.0 * tr_yyyyz_yy[i] * tbe_0 + 4.0 * tr_yyyyz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_yzz[i] = -2.0 * tr_yyyy_yzz[i] * tbe_0 - 4.0 * tr_yyyyz_yz[i] * tbe_0 + 4.0 * tr_yyyyz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyy_zzz[i] = -2.0 * tr_yyyy_zzz[i] * tbe_0 - 6.0 * tr_yyyyz_zz[i] * tbe_0 + 4.0 * tr_yyyyz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1310-1320 components of targeted buffer : GF

    auto tr_z_0_z_yyyz_xxx = pbuffer.data(idx_op_geom_110_gf + 1310);

    auto tr_z_0_z_yyyz_xxy = pbuffer.data(idx_op_geom_110_gf + 1311);

    auto tr_z_0_z_yyyz_xxz = pbuffer.data(idx_op_geom_110_gf + 1312);

    auto tr_z_0_z_yyyz_xyy = pbuffer.data(idx_op_geom_110_gf + 1313);

    auto tr_z_0_z_yyyz_xyz = pbuffer.data(idx_op_geom_110_gf + 1314);

    auto tr_z_0_z_yyyz_xzz = pbuffer.data(idx_op_geom_110_gf + 1315);

    auto tr_z_0_z_yyyz_yyy = pbuffer.data(idx_op_geom_110_gf + 1316);

    auto tr_z_0_z_yyyz_yyz = pbuffer.data(idx_op_geom_110_gf + 1317);

    auto tr_z_0_z_yyyz_yzz = pbuffer.data(idx_op_geom_110_gf + 1318);

    auto tr_z_0_z_yyyz_zzz = pbuffer.data(idx_op_geom_110_gf + 1319);

    #pragma omp simd aligned(tr_yyy_xx, tr_yyy_xxxz, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyyzz_xx, tr_yyyzz_xxxz, tr_yyyzz_xxyz, tr_yyyzz_xxzz, tr_yyyzz_xy, tr_yyyzz_xyyz, tr_yyyzz_xyzz, tr_yyyzz_xz, tr_yyyzz_xzzz, tr_yyyzz_yy, tr_yyyzz_yyyz, tr_yyyzz_yyzz, tr_yyyzz_yz, tr_yyyzz_yzzz, tr_yyyzz_zz, tr_yyyzz_zzzz, tr_yyyzzz_xxx, tr_yyyzzz_xxy, tr_yyyzzz_xxz, tr_yyyzzz_xyy, tr_yyyzzz_xyz, tr_yyyzzz_xzz, tr_yyyzzz_yyy, tr_yyyzzz_yyz, tr_yyyzzz_yzz, tr_yyyzzz_zzz, tr_z_0_z_yyyz_xxx, tr_z_0_z_yyyz_xxy, tr_z_0_z_yyyz_xxz, tr_z_0_z_yyyz_xyy, tr_z_0_z_yyyz_xyz, tr_z_0_z_yyyz_xzz, tr_z_0_z_yyyz_yyy, tr_z_0_z_yyyz_yyz, tr_z_0_z_yyyz_yzz, tr_z_0_z_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyyz_xxx[i] = -2.0 * tr_yyy_xxxz[i] * tke_0 - 6.0 * tr_yyyz_xxx[i] * tbe_0 + 4.0 * tr_yyyzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_xxy[i] = -2.0 * tr_yyy_xxyz[i] * tke_0 - 6.0 * tr_yyyz_xxy[i] * tbe_0 + 4.0 * tr_yyyzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_xxz[i] = tr_yyy_xx[i] - 2.0 * tr_yyy_xxzz[i] * tke_0 - 6.0 * tr_yyyz_xxz[i] * tbe_0 - 2.0 * tr_yyyzz_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_xyy[i] = -2.0 * tr_yyy_xyyz[i] * tke_0 - 6.0 * tr_yyyz_xyy[i] * tbe_0 + 4.0 * tr_yyyzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_xyz[i] = tr_yyy_xy[i] - 2.0 * tr_yyy_xyzz[i] * tke_0 - 6.0 * tr_yyyz_xyz[i] * tbe_0 - 2.0 * tr_yyyzz_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_xzz[i] = 2.0 * tr_yyy_xz[i] - 2.0 * tr_yyy_xzzz[i] * tke_0 - 6.0 * tr_yyyz_xzz[i] * tbe_0 - 4.0 * tr_yyyzz_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_yyy[i] = -2.0 * tr_yyy_yyyz[i] * tke_0 - 6.0 * tr_yyyz_yyy[i] * tbe_0 + 4.0 * tr_yyyzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_yyz[i] = tr_yyy_yy[i] - 2.0 * tr_yyy_yyzz[i] * tke_0 - 6.0 * tr_yyyz_yyz[i] * tbe_0 - 2.0 * tr_yyyzz_yy[i] * tbe_0 + 4.0 * tr_yyyzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_yzz[i] = 2.0 * tr_yyy_yz[i] - 2.0 * tr_yyy_yzzz[i] * tke_0 - 6.0 * tr_yyyz_yzz[i] * tbe_0 - 4.0 * tr_yyyzz_yz[i] * tbe_0 + 4.0 * tr_yyyzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyyz_zzz[i] = 3.0 * tr_yyy_zz[i] - 2.0 * tr_yyy_zzzz[i] * tke_0 - 6.0 * tr_yyyz_zzz[i] * tbe_0 - 6.0 * tr_yyyzz_zz[i] * tbe_0 + 4.0 * tr_yyyzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1320-1330 components of targeted buffer : GF

    auto tr_z_0_z_yyzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1320);

    auto tr_z_0_z_yyzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1321);

    auto tr_z_0_z_yyzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1322);

    auto tr_z_0_z_yyzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1323);

    auto tr_z_0_z_yyzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1324);

    auto tr_z_0_z_yyzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1325);

    auto tr_z_0_z_yyzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1326);

    auto tr_z_0_z_yyzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1327);

    auto tr_z_0_z_yyzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1328);

    auto tr_z_0_z_yyzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1329);

    #pragma omp simd aligned(tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyz_xx, tr_yyz_xxxz, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yyzzz_xx, tr_yyzzz_xxxz, tr_yyzzz_xxyz, tr_yyzzz_xxzz, tr_yyzzz_xy, tr_yyzzz_xyyz, tr_yyzzz_xyzz, tr_yyzzz_xz, tr_yyzzz_xzzz, tr_yyzzz_yy, tr_yyzzz_yyyz, tr_yyzzz_yyzz, tr_yyzzz_yz, tr_yyzzz_yzzz, tr_yyzzz_zz, tr_yyzzz_zzzz, tr_yyzzzz_xxx, tr_yyzzzz_xxy, tr_yyzzzz_xxz, tr_yyzzzz_xyy, tr_yyzzzz_xyz, tr_yyzzzz_xzz, tr_yyzzzz_yyy, tr_yyzzzz_yyz, tr_yyzzzz_yzz, tr_yyzzzz_zzz, tr_z_0_z_yyzz_xxx, tr_z_0_z_yyzz_xxy, tr_z_0_z_yyzz_xxz, tr_z_0_z_yyzz_xyy, tr_z_0_z_yyzz_xyz, tr_z_0_z_yyzz_xzz, tr_z_0_z_yyzz_yyy, tr_z_0_z_yyzz_yyz, tr_z_0_z_yyzz_yzz, tr_z_0_z_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyzz_xxx[i] = 2.0 * tr_yy_xxx[i] - 4.0 * tr_yyz_xxxz[i] * tke_0 - 10.0 * tr_yyzz_xxx[i] * tbe_0 + 4.0 * tr_yyzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_xxy[i] = 2.0 * tr_yy_xxy[i] - 4.0 * tr_yyz_xxyz[i] * tke_0 - 10.0 * tr_yyzz_xxy[i] * tbe_0 + 4.0 * tr_yyzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_xxz[i] = 2.0 * tr_yy_xxz[i] + 2.0 * tr_yyz_xx[i] - 4.0 * tr_yyz_xxzz[i] * tke_0 - 10.0 * tr_yyzz_xxz[i] * tbe_0 - 2.0 * tr_yyzzz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_xyy[i] = 2.0 * tr_yy_xyy[i] - 4.0 * tr_yyz_xyyz[i] * tke_0 - 10.0 * tr_yyzz_xyy[i] * tbe_0 + 4.0 * tr_yyzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_xyz[i] = 2.0 * tr_yy_xyz[i] + 2.0 * tr_yyz_xy[i] - 4.0 * tr_yyz_xyzz[i] * tke_0 - 10.0 * tr_yyzz_xyz[i] * tbe_0 - 2.0 * tr_yyzzz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_xzz[i] = 2.0 * tr_yy_xzz[i] + 4.0 * tr_yyz_xz[i] - 4.0 * tr_yyz_xzzz[i] * tke_0 - 10.0 * tr_yyzz_xzz[i] * tbe_0 - 4.0 * tr_yyzzz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_yyy[i] = 2.0 * tr_yy_yyy[i] - 4.0 * tr_yyz_yyyz[i] * tke_0 - 10.0 * tr_yyzz_yyy[i] * tbe_0 + 4.0 * tr_yyzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_yyz[i] = 2.0 * tr_yy_yyz[i] + 2.0 * tr_yyz_yy[i] - 4.0 * tr_yyz_yyzz[i] * tke_0 - 10.0 * tr_yyzz_yyz[i] * tbe_0 - 2.0 * tr_yyzzz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_yzz[i] = 2.0 * tr_yy_yzz[i] + 4.0 * tr_yyz_yz[i] - 4.0 * tr_yyz_yzzz[i] * tke_0 - 10.0 * tr_yyzz_yzz[i] * tbe_0 - 4.0 * tr_yyzzz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyzz_zzz[i] = 2.0 * tr_yy_zzz[i] + 6.0 * tr_yyz_zz[i] - 4.0 * tr_yyz_zzzz[i] * tke_0 - 10.0 * tr_yyzz_zzz[i] * tbe_0 - 6.0 * tr_yyzzz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1330-1340 components of targeted buffer : GF

    auto tr_z_0_z_yzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1330);

    auto tr_z_0_z_yzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1331);

    auto tr_z_0_z_yzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1332);

    auto tr_z_0_z_yzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1333);

    auto tr_z_0_z_yzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1334);

    auto tr_z_0_z_yzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1335);

    auto tr_z_0_z_yzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1336);

    auto tr_z_0_z_yzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1337);

    auto tr_z_0_z_yzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1338);

    auto tr_z_0_z_yzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1339);

    #pragma omp simd aligned(tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxz, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xxxz, tr_yzzzz_xxyz, tr_yzzzz_xxzz, tr_yzzzz_xy, tr_yzzzz_xyyz, tr_yzzzz_xyzz, tr_yzzzz_xz, tr_yzzzz_xzzz, tr_yzzzz_yy, tr_yzzzz_yyyz, tr_yzzzz_yyzz, tr_yzzzz_yz, tr_yzzzz_yzzz, tr_yzzzz_zz, tr_yzzzz_zzzz, tr_yzzzzz_xxx, tr_yzzzzz_xxy, tr_yzzzzz_xxz, tr_yzzzzz_xyy, tr_yzzzzz_xyz, tr_yzzzzz_xzz, tr_yzzzzz_yyy, tr_yzzzzz_yyz, tr_yzzzzz_yzz, tr_yzzzzz_zzz, tr_z_0_z_yzzz_xxx, tr_z_0_z_yzzz_xxy, tr_z_0_z_yzzz_xxz, tr_z_0_z_yzzz_xyy, tr_z_0_z_yzzz_xyz, tr_z_0_z_yzzz_xzz, tr_z_0_z_yzzz_yyy, tr_z_0_z_yzzz_yyz, tr_z_0_z_yzzz_yzz, tr_z_0_z_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yzzz_xxx[i] = 6.0 * tr_yz_xxx[i] - 6.0 * tr_yzz_xxxz[i] * tke_0 - 14.0 * tr_yzzz_xxx[i] * tbe_0 + 4.0 * tr_yzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_xxy[i] = 6.0 * tr_yz_xxy[i] - 6.0 * tr_yzz_xxyz[i] * tke_0 - 14.0 * tr_yzzz_xxy[i] * tbe_0 + 4.0 * tr_yzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_xxz[i] = 6.0 * tr_yz_xxz[i] + 3.0 * tr_yzz_xx[i] - 6.0 * tr_yzz_xxzz[i] * tke_0 - 14.0 * tr_yzzz_xxz[i] * tbe_0 - 2.0 * tr_yzzzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_xyy[i] = 6.0 * tr_yz_xyy[i] - 6.0 * tr_yzz_xyyz[i] * tke_0 - 14.0 * tr_yzzz_xyy[i] * tbe_0 + 4.0 * tr_yzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_xyz[i] = 6.0 * tr_yz_xyz[i] + 3.0 * tr_yzz_xy[i] - 6.0 * tr_yzz_xyzz[i] * tke_0 - 14.0 * tr_yzzz_xyz[i] * tbe_0 - 2.0 * tr_yzzzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_xzz[i] = 6.0 * tr_yz_xzz[i] + 6.0 * tr_yzz_xz[i] - 6.0 * tr_yzz_xzzz[i] * tke_0 - 14.0 * tr_yzzz_xzz[i] * tbe_0 - 4.0 * tr_yzzzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_yyy[i] = 6.0 * tr_yz_yyy[i] - 6.0 * tr_yzz_yyyz[i] * tke_0 - 14.0 * tr_yzzz_yyy[i] * tbe_0 + 4.0 * tr_yzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_yyz[i] = 6.0 * tr_yz_yyz[i] + 3.0 * tr_yzz_yy[i] - 6.0 * tr_yzz_yyzz[i] * tke_0 - 14.0 * tr_yzzz_yyz[i] * tbe_0 - 2.0 * tr_yzzzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_yzz[i] = 6.0 * tr_yz_yzz[i] + 6.0 * tr_yzz_yz[i] - 6.0 * tr_yzz_yzzz[i] * tke_0 - 14.0 * tr_yzzz_yzz[i] * tbe_0 - 4.0 * tr_yzzzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzzz_zzz[i] = 6.0 * tr_yz_zzz[i] + 9.0 * tr_yzz_zz[i] - 6.0 * tr_yzz_zzzz[i] * tke_0 - 14.0 * tr_yzzz_zzz[i] * tbe_0 - 6.0 * tr_yzzzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

    // Set up 1340-1350 components of targeted buffer : GF

    auto tr_z_0_z_zzzz_xxx = pbuffer.data(idx_op_geom_110_gf + 1340);

    auto tr_z_0_z_zzzz_xxy = pbuffer.data(idx_op_geom_110_gf + 1341);

    auto tr_z_0_z_zzzz_xxz = pbuffer.data(idx_op_geom_110_gf + 1342);

    auto tr_z_0_z_zzzz_xyy = pbuffer.data(idx_op_geom_110_gf + 1343);

    auto tr_z_0_z_zzzz_xyz = pbuffer.data(idx_op_geom_110_gf + 1344);

    auto tr_z_0_z_zzzz_xzz = pbuffer.data(idx_op_geom_110_gf + 1345);

    auto tr_z_0_z_zzzz_yyy = pbuffer.data(idx_op_geom_110_gf + 1346);

    auto tr_z_0_z_zzzz_yyz = pbuffer.data(idx_op_geom_110_gf + 1347);

    auto tr_z_0_z_zzzz_yzz = pbuffer.data(idx_op_geom_110_gf + 1348);

    auto tr_z_0_z_zzzz_zzz = pbuffer.data(idx_op_geom_110_gf + 1349);

    #pragma omp simd aligned(tr_z_0_z_zzzz_xxx, tr_z_0_z_zzzz_xxy, tr_z_0_z_zzzz_xxz, tr_z_0_z_zzzz_xyy, tr_z_0_z_zzzz_xyz, tr_z_0_z_zzzz_xzz, tr_z_0_z_zzzz_yyy, tr_z_0_z_zzzz_yyz, tr_z_0_z_zzzz_yzz, tr_z_0_z_zzzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzz_xx, tr_zzz_xxxz, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, tr_zzzzz_xx, tr_zzzzz_xxxz, tr_zzzzz_xxyz, tr_zzzzz_xxzz, tr_zzzzz_xy, tr_zzzzz_xyyz, tr_zzzzz_xyzz, tr_zzzzz_xz, tr_zzzzz_xzzz, tr_zzzzz_yy, tr_zzzzz_yyyz, tr_zzzzz_yyzz, tr_zzzzz_yz, tr_zzzzz_yzzz, tr_zzzzz_zz, tr_zzzzz_zzzz, tr_zzzzzz_xxx, tr_zzzzzz_xxy, tr_zzzzzz_xxz, tr_zzzzzz_xyy, tr_zzzzzz_xyz, tr_zzzzzz_xzz, tr_zzzzzz_yyy, tr_zzzzzz_yyz, tr_zzzzzz_yzz, tr_zzzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zzzz_xxx[i] = 12.0 * tr_zz_xxx[i] - 8.0 * tr_zzz_xxxz[i] * tke_0 - 18.0 * tr_zzzz_xxx[i] * tbe_0 + 4.0 * tr_zzzzz_xxxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxx[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_xxy[i] = 12.0 * tr_zz_xxy[i] - 8.0 * tr_zzz_xxyz[i] * tke_0 - 18.0 * tr_zzzz_xxy[i] * tbe_0 + 4.0 * tr_zzzzz_xxyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_xxz[i] = 12.0 * tr_zz_xxz[i] + 4.0 * tr_zzz_xx[i] - 8.0 * tr_zzz_xxzz[i] * tke_0 - 18.0 * tr_zzzz_xxz[i] * tbe_0 - 2.0 * tr_zzzzz_xx[i] * tbe_0 + 4.0 * tr_zzzzz_xxzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xxz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_xyy[i] = 12.0 * tr_zz_xyy[i] - 8.0 * tr_zzz_xyyz[i] * tke_0 - 18.0 * tr_zzzz_xyy[i] * tbe_0 + 4.0 * tr_zzzzz_xyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_xyz[i] = 12.0 * tr_zz_xyz[i] + 4.0 * tr_zzz_xy[i] - 8.0 * tr_zzz_xyzz[i] * tke_0 - 18.0 * tr_zzzz_xyz[i] * tbe_0 - 2.0 * tr_zzzzz_xy[i] * tbe_0 + 4.0 * tr_zzzzz_xyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_xzz[i] = 12.0 * tr_zz_xzz[i] + 8.0 * tr_zzz_xz[i] - 8.0 * tr_zzz_xzzz[i] * tke_0 - 18.0 * tr_zzzz_xzz[i] * tbe_0 - 4.0 * tr_zzzzz_xz[i] * tbe_0 + 4.0 * tr_zzzzz_xzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_xzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_yyy[i] = 12.0 * tr_zz_yyy[i] - 8.0 * tr_zzz_yyyz[i] * tke_0 - 18.0 * tr_zzzz_yyy[i] * tbe_0 + 4.0 * tr_zzzzz_yyyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yyy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_yyz[i] = 12.0 * tr_zz_yyz[i] + 4.0 * tr_zzz_yy[i] - 8.0 * tr_zzz_yyzz[i] * tke_0 - 18.0 * tr_zzzz_yyz[i] * tbe_0 - 2.0 * tr_zzzzz_yy[i] * tbe_0 + 4.0 * tr_zzzzz_yyzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yyz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_yzz[i] = 12.0 * tr_zz_yzz[i] + 8.0 * tr_zzz_yz[i] - 8.0 * tr_zzz_yzzz[i] * tke_0 - 18.0 * tr_zzzz_yzz[i] * tbe_0 - 4.0 * tr_zzzzz_yz[i] * tbe_0 + 4.0 * tr_zzzzz_yzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_yzz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzzz_zzz[i] = 12.0 * tr_zz_zzz[i] + 12.0 * tr_zzz_zz[i] - 8.0 * tr_zzz_zzzz[i] * tke_0 - 18.0 * tr_zzzz_zzz[i] * tbe_0 - 6.0 * tr_zzzzz_zz[i] * tbe_0 + 4.0 * tr_zzzzz_zzzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzzz_zzz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

