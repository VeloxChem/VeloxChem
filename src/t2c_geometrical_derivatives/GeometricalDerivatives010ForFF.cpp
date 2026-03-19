#include "GeometricalDerivatives010ForFF.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_ff(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_ff,
                         const int idx_op_df,
                         const int idx_op_fd,
                         const int idx_op_fg,
                         const int idx_op_gf,
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

    // Set up 0-10 components of targeted buffer : FF

    auto tr_0_0_x_xxx_xxx = pbuffer.data(idx_op_geom_010_ff);

    auto tr_0_0_x_xxx_xxy = pbuffer.data(idx_op_geom_010_ff + 1);

    auto tr_0_0_x_xxx_xxz = pbuffer.data(idx_op_geom_010_ff + 2);

    auto tr_0_0_x_xxx_xyy = pbuffer.data(idx_op_geom_010_ff + 3);

    auto tr_0_0_x_xxx_xyz = pbuffer.data(idx_op_geom_010_ff + 4);

    auto tr_0_0_x_xxx_xzz = pbuffer.data(idx_op_geom_010_ff + 5);

    auto tr_0_0_x_xxx_yyy = pbuffer.data(idx_op_geom_010_ff + 6);

    auto tr_0_0_x_xxx_yyz = pbuffer.data(idx_op_geom_010_ff + 7);

    auto tr_0_0_x_xxx_yzz = pbuffer.data(idx_op_geom_010_ff + 8);

    auto tr_0_0_x_xxx_zzz = pbuffer.data(idx_op_geom_010_ff + 9);

    #pragma omp simd aligned(tr_0_0_x_xxx_xxx, tr_0_0_x_xxx_xxy, tr_0_0_x_xxx_xxz, tr_0_0_x_xxx_xyy, tr_0_0_x_xxx_xyz, tr_0_0_x_xxx_xzz, tr_0_0_x_xxx_yyy, tr_0_0_x_xxx_yyz, tr_0_0_x_xxx_yzz, tr_0_0_x_xxx_zzz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxx_xx, tr_xxx_xxxx, tr_xxx_xxxy, tr_xxx_xxxz, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxx_xxx[i] = 2.0 * tr_xxxx_xxx[i] * tbe_0 + 2.0 * tr_xxx_xxxx[i] * tke_0 - 3.0 * tr_xx_xxx[i] - 3.0 * tr_xxx_xx[i];

        tr_0_0_x_xxx_xxy[i] = 2.0 * tr_xxxx_xxy[i] * tbe_0 + 2.0 * tr_xxx_xxxy[i] * tke_0 - 3.0 * tr_xx_xxy[i] - 2.0 * tr_xxx_xy[i];

        tr_0_0_x_xxx_xxz[i] = 2.0 * tr_xxxx_xxz[i] * tbe_0 + 2.0 * tr_xxx_xxxz[i] * tke_0 - 3.0 * tr_xx_xxz[i] - 2.0 * tr_xxx_xz[i];

        tr_0_0_x_xxx_xyy[i] = 2.0 * tr_xxxx_xyy[i] * tbe_0 + 2.0 * tr_xxx_xxyy[i] * tke_0 - 3.0 * tr_xx_xyy[i] - tr_xxx_yy[i];

        tr_0_0_x_xxx_xyz[i] = 2.0 * tr_xxxx_xyz[i] * tbe_0 + 2.0 * tr_xxx_xxyz[i] * tke_0 - 3.0 * tr_xx_xyz[i] - tr_xxx_yz[i];

        tr_0_0_x_xxx_xzz[i] = 2.0 * tr_xxxx_xzz[i] * tbe_0 + 2.0 * tr_xxx_xxzz[i] * tke_0 - 3.0 * tr_xx_xzz[i] - tr_xxx_zz[i];

        tr_0_0_x_xxx_yyy[i] = 2.0 * tr_xxxx_yyy[i] * tbe_0 + 2.0 * tr_xxx_xyyy[i] * tke_0 - 3.0 * tr_xx_yyy[i];

        tr_0_0_x_xxx_yyz[i] = 2.0 * tr_xxxx_yyz[i] * tbe_0 + 2.0 * tr_xxx_xyyz[i] * tke_0 - 3.0 * tr_xx_yyz[i];

        tr_0_0_x_xxx_yzz[i] = 2.0 * tr_xxxx_yzz[i] * tbe_0 + 2.0 * tr_xxx_xyzz[i] * tke_0 - 3.0 * tr_xx_yzz[i];

        tr_0_0_x_xxx_zzz[i] = 2.0 * tr_xxxx_zzz[i] * tbe_0 + 2.0 * tr_xxx_xzzz[i] * tke_0 - 3.0 * tr_xx_zzz[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto tr_0_0_x_xxy_xxx = pbuffer.data(idx_op_geom_010_ff + 10);

    auto tr_0_0_x_xxy_xxy = pbuffer.data(idx_op_geom_010_ff + 11);

    auto tr_0_0_x_xxy_xxz = pbuffer.data(idx_op_geom_010_ff + 12);

    auto tr_0_0_x_xxy_xyy = pbuffer.data(idx_op_geom_010_ff + 13);

    auto tr_0_0_x_xxy_xyz = pbuffer.data(idx_op_geom_010_ff + 14);

    auto tr_0_0_x_xxy_xzz = pbuffer.data(idx_op_geom_010_ff + 15);

    auto tr_0_0_x_xxy_yyy = pbuffer.data(idx_op_geom_010_ff + 16);

    auto tr_0_0_x_xxy_yyz = pbuffer.data(idx_op_geom_010_ff + 17);

    auto tr_0_0_x_xxy_yzz = pbuffer.data(idx_op_geom_010_ff + 18);

    auto tr_0_0_x_xxy_zzz = pbuffer.data(idx_op_geom_010_ff + 19);

    #pragma omp simd aligned(tr_0_0_x_xxy_xxx, tr_0_0_x_xxy_xxy, tr_0_0_x_xxy_xxz, tr_0_0_x_xxy_xyy, tr_0_0_x_xxy_xyz, tr_0_0_x_xxy_xzz, tr_0_0_x_xxy_yyy, tr_0_0_x_xxy_yyz, tr_0_0_x_xxy_yzz, tr_0_0_x_xxy_zzz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, tr_xxy_xx, tr_xxy_xxxx, tr_xxy_xxxy, tr_xxy_xxxz, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxy_xxx[i] = 2.0 * tr_xxxy_xxx[i] * tbe_0 + 2.0 * tr_xxy_xxxx[i] * tke_0 - 2.0 * tr_xy_xxx[i] - 3.0 * tr_xxy_xx[i];

        tr_0_0_x_xxy_xxy[i] = 2.0 * tr_xxxy_xxy[i] * tbe_0 + 2.0 * tr_xxy_xxxy[i] * tke_0 - 2.0 * tr_xy_xxy[i] - 2.0 * tr_xxy_xy[i];

        tr_0_0_x_xxy_xxz[i] = 2.0 * tr_xxxy_xxz[i] * tbe_0 + 2.0 * tr_xxy_xxxz[i] * tke_0 - 2.0 * tr_xy_xxz[i] - 2.0 * tr_xxy_xz[i];

        tr_0_0_x_xxy_xyy[i] = 2.0 * tr_xxxy_xyy[i] * tbe_0 + 2.0 * tr_xxy_xxyy[i] * tke_0 - 2.0 * tr_xy_xyy[i] - tr_xxy_yy[i];

        tr_0_0_x_xxy_xyz[i] = 2.0 * tr_xxxy_xyz[i] * tbe_0 + 2.0 * tr_xxy_xxyz[i] * tke_0 - 2.0 * tr_xy_xyz[i] - tr_xxy_yz[i];

        tr_0_0_x_xxy_xzz[i] = 2.0 * tr_xxxy_xzz[i] * tbe_0 + 2.0 * tr_xxy_xxzz[i] * tke_0 - 2.0 * tr_xy_xzz[i] - tr_xxy_zz[i];

        tr_0_0_x_xxy_yyy[i] = 2.0 * tr_xxxy_yyy[i] * tbe_0 + 2.0 * tr_xxy_xyyy[i] * tke_0 - 2.0 * tr_xy_yyy[i];

        tr_0_0_x_xxy_yyz[i] = 2.0 * tr_xxxy_yyz[i] * tbe_0 + 2.0 * tr_xxy_xyyz[i] * tke_0 - 2.0 * tr_xy_yyz[i];

        tr_0_0_x_xxy_yzz[i] = 2.0 * tr_xxxy_yzz[i] * tbe_0 + 2.0 * tr_xxy_xyzz[i] * tke_0 - 2.0 * tr_xy_yzz[i];

        tr_0_0_x_xxy_zzz[i] = 2.0 * tr_xxxy_zzz[i] * tbe_0 + 2.0 * tr_xxy_xzzz[i] * tke_0 - 2.0 * tr_xy_zzz[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto tr_0_0_x_xxz_xxx = pbuffer.data(idx_op_geom_010_ff + 20);

    auto tr_0_0_x_xxz_xxy = pbuffer.data(idx_op_geom_010_ff + 21);

    auto tr_0_0_x_xxz_xxz = pbuffer.data(idx_op_geom_010_ff + 22);

    auto tr_0_0_x_xxz_xyy = pbuffer.data(idx_op_geom_010_ff + 23);

    auto tr_0_0_x_xxz_xyz = pbuffer.data(idx_op_geom_010_ff + 24);

    auto tr_0_0_x_xxz_xzz = pbuffer.data(idx_op_geom_010_ff + 25);

    auto tr_0_0_x_xxz_yyy = pbuffer.data(idx_op_geom_010_ff + 26);

    auto tr_0_0_x_xxz_yyz = pbuffer.data(idx_op_geom_010_ff + 27);

    auto tr_0_0_x_xxz_yzz = pbuffer.data(idx_op_geom_010_ff + 28);

    auto tr_0_0_x_xxz_zzz = pbuffer.data(idx_op_geom_010_ff + 29);

    #pragma omp simd aligned(tr_0_0_x_xxz_xxx, tr_0_0_x_xxz_xxy, tr_0_0_x_xxz_xxz, tr_0_0_x_xxz_xyy, tr_0_0_x_xxz_xyz, tr_0_0_x_xxz_xzz, tr_0_0_x_xxz_yyy, tr_0_0_x_xxz_yyz, tr_0_0_x_xxz_yzz, tr_0_0_x_xxz_zzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, tr_xxz_xx, tr_xxz_xxxx, tr_xxz_xxxy, tr_xxz_xxxz, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxz_xxx[i] = 2.0 * tr_xxxz_xxx[i] * tbe_0 + 2.0 * tr_xxz_xxxx[i] * tke_0 - 2.0 * tr_xz_xxx[i] - 3.0 * tr_xxz_xx[i];

        tr_0_0_x_xxz_xxy[i] = 2.0 * tr_xxxz_xxy[i] * tbe_0 + 2.0 * tr_xxz_xxxy[i] * tke_0 - 2.0 * tr_xz_xxy[i] - 2.0 * tr_xxz_xy[i];

        tr_0_0_x_xxz_xxz[i] = 2.0 * tr_xxxz_xxz[i] * tbe_0 + 2.0 * tr_xxz_xxxz[i] * tke_0 - 2.0 * tr_xz_xxz[i] - 2.0 * tr_xxz_xz[i];

        tr_0_0_x_xxz_xyy[i] = 2.0 * tr_xxxz_xyy[i] * tbe_0 + 2.0 * tr_xxz_xxyy[i] * tke_0 - 2.0 * tr_xz_xyy[i] - tr_xxz_yy[i];

        tr_0_0_x_xxz_xyz[i] = 2.0 * tr_xxxz_xyz[i] * tbe_0 + 2.0 * tr_xxz_xxyz[i] * tke_0 - 2.0 * tr_xz_xyz[i] - tr_xxz_yz[i];

        tr_0_0_x_xxz_xzz[i] = 2.0 * tr_xxxz_xzz[i] * tbe_0 + 2.0 * tr_xxz_xxzz[i] * tke_0 - 2.0 * tr_xz_xzz[i] - tr_xxz_zz[i];

        tr_0_0_x_xxz_yyy[i] = 2.0 * tr_xxxz_yyy[i] * tbe_0 + 2.0 * tr_xxz_xyyy[i] * tke_0 - 2.0 * tr_xz_yyy[i];

        tr_0_0_x_xxz_yyz[i] = 2.0 * tr_xxxz_yyz[i] * tbe_0 + 2.0 * tr_xxz_xyyz[i] * tke_0 - 2.0 * tr_xz_yyz[i];

        tr_0_0_x_xxz_yzz[i] = 2.0 * tr_xxxz_yzz[i] * tbe_0 + 2.0 * tr_xxz_xyzz[i] * tke_0 - 2.0 * tr_xz_yzz[i];

        tr_0_0_x_xxz_zzz[i] = 2.0 * tr_xxxz_zzz[i] * tbe_0 + 2.0 * tr_xxz_xzzz[i] * tke_0 - 2.0 * tr_xz_zzz[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto tr_0_0_x_xyy_xxx = pbuffer.data(idx_op_geom_010_ff + 30);

    auto tr_0_0_x_xyy_xxy = pbuffer.data(idx_op_geom_010_ff + 31);

    auto tr_0_0_x_xyy_xxz = pbuffer.data(idx_op_geom_010_ff + 32);

    auto tr_0_0_x_xyy_xyy = pbuffer.data(idx_op_geom_010_ff + 33);

    auto tr_0_0_x_xyy_xyz = pbuffer.data(idx_op_geom_010_ff + 34);

    auto tr_0_0_x_xyy_xzz = pbuffer.data(idx_op_geom_010_ff + 35);

    auto tr_0_0_x_xyy_yyy = pbuffer.data(idx_op_geom_010_ff + 36);

    auto tr_0_0_x_xyy_yyz = pbuffer.data(idx_op_geom_010_ff + 37);

    auto tr_0_0_x_xyy_yzz = pbuffer.data(idx_op_geom_010_ff + 38);

    auto tr_0_0_x_xyy_zzz = pbuffer.data(idx_op_geom_010_ff + 39);

    #pragma omp simd aligned(tr_0_0_x_xyy_xxx, tr_0_0_x_xyy_xxy, tr_0_0_x_xyy_xxz, tr_0_0_x_xyy_xyy, tr_0_0_x_xyy_xyz, tr_0_0_x_xyy_xzz, tr_0_0_x_xyy_yyy, tr_0_0_x_xyy_yyz, tr_0_0_x_xyy_yzz, tr_0_0_x_xyy_zzz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, tr_xyy_xx, tr_xyy_xxxx, tr_xyy_xxxy, tr_xyy_xxxz, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyy_xxx[i] = 2.0 * tr_xxyy_xxx[i] * tbe_0 + 2.0 * tr_xyy_xxxx[i] * tke_0 - tr_yy_xxx[i] - 3.0 * tr_xyy_xx[i];

        tr_0_0_x_xyy_xxy[i] = 2.0 * tr_xxyy_xxy[i] * tbe_0 + 2.0 * tr_xyy_xxxy[i] * tke_0 - tr_yy_xxy[i] - 2.0 * tr_xyy_xy[i];

        tr_0_0_x_xyy_xxz[i] = 2.0 * tr_xxyy_xxz[i] * tbe_0 + 2.0 * tr_xyy_xxxz[i] * tke_0 - tr_yy_xxz[i] - 2.0 * tr_xyy_xz[i];

        tr_0_0_x_xyy_xyy[i] = 2.0 * tr_xxyy_xyy[i] * tbe_0 + 2.0 * tr_xyy_xxyy[i] * tke_0 - tr_yy_xyy[i] - tr_xyy_yy[i];

        tr_0_0_x_xyy_xyz[i] = 2.0 * tr_xxyy_xyz[i] * tbe_0 + 2.0 * tr_xyy_xxyz[i] * tke_0 - tr_yy_xyz[i] - tr_xyy_yz[i];

        tr_0_0_x_xyy_xzz[i] = 2.0 * tr_xxyy_xzz[i] * tbe_0 + 2.0 * tr_xyy_xxzz[i] * tke_0 - tr_yy_xzz[i] - tr_xyy_zz[i];

        tr_0_0_x_xyy_yyy[i] = 2.0 * tr_xxyy_yyy[i] * tbe_0 + 2.0 * tr_xyy_xyyy[i] * tke_0 - tr_yy_yyy[i];

        tr_0_0_x_xyy_yyz[i] = 2.0 * tr_xxyy_yyz[i] * tbe_0 + 2.0 * tr_xyy_xyyz[i] * tke_0 - tr_yy_yyz[i];

        tr_0_0_x_xyy_yzz[i] = 2.0 * tr_xxyy_yzz[i] * tbe_0 + 2.0 * tr_xyy_xyzz[i] * tke_0 - tr_yy_yzz[i];

        tr_0_0_x_xyy_zzz[i] = 2.0 * tr_xxyy_zzz[i] * tbe_0 + 2.0 * tr_xyy_xzzz[i] * tke_0 - tr_yy_zzz[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto tr_0_0_x_xyz_xxx = pbuffer.data(idx_op_geom_010_ff + 40);

    auto tr_0_0_x_xyz_xxy = pbuffer.data(idx_op_geom_010_ff + 41);

    auto tr_0_0_x_xyz_xxz = pbuffer.data(idx_op_geom_010_ff + 42);

    auto tr_0_0_x_xyz_xyy = pbuffer.data(idx_op_geom_010_ff + 43);

    auto tr_0_0_x_xyz_xyz = pbuffer.data(idx_op_geom_010_ff + 44);

    auto tr_0_0_x_xyz_xzz = pbuffer.data(idx_op_geom_010_ff + 45);

    auto tr_0_0_x_xyz_yyy = pbuffer.data(idx_op_geom_010_ff + 46);

    auto tr_0_0_x_xyz_yyz = pbuffer.data(idx_op_geom_010_ff + 47);

    auto tr_0_0_x_xyz_yzz = pbuffer.data(idx_op_geom_010_ff + 48);

    auto tr_0_0_x_xyz_zzz = pbuffer.data(idx_op_geom_010_ff + 49);

    #pragma omp simd aligned(tr_0_0_x_xyz_xxx, tr_0_0_x_xyz_xxy, tr_0_0_x_xyz_xxz, tr_0_0_x_xyz_xyy, tr_0_0_x_xyz_xyz, tr_0_0_x_xyz_xzz, tr_0_0_x_xyz_yyy, tr_0_0_x_xyz_yyz, tr_0_0_x_xyz_yzz, tr_0_0_x_xyz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xyz_xx, tr_xyz_xxxx, tr_xyz_xxxy, tr_xyz_xxxz, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyz_xxx[i] = 2.0 * tr_xxyz_xxx[i] * tbe_0 + 2.0 * tr_xyz_xxxx[i] * tke_0 - tr_yz_xxx[i] - 3.0 * tr_xyz_xx[i];

        tr_0_0_x_xyz_xxy[i] = 2.0 * tr_xxyz_xxy[i] * tbe_0 + 2.0 * tr_xyz_xxxy[i] * tke_0 - tr_yz_xxy[i] - 2.0 * tr_xyz_xy[i];

        tr_0_0_x_xyz_xxz[i] = 2.0 * tr_xxyz_xxz[i] * tbe_0 + 2.0 * tr_xyz_xxxz[i] * tke_0 - tr_yz_xxz[i] - 2.0 * tr_xyz_xz[i];

        tr_0_0_x_xyz_xyy[i] = 2.0 * tr_xxyz_xyy[i] * tbe_0 + 2.0 * tr_xyz_xxyy[i] * tke_0 - tr_yz_xyy[i] - tr_xyz_yy[i];

        tr_0_0_x_xyz_xyz[i] = 2.0 * tr_xxyz_xyz[i] * tbe_0 + 2.0 * tr_xyz_xxyz[i] * tke_0 - tr_yz_xyz[i] - tr_xyz_yz[i];

        tr_0_0_x_xyz_xzz[i] = 2.0 * tr_xxyz_xzz[i] * tbe_0 + 2.0 * tr_xyz_xxzz[i] * tke_0 - tr_yz_xzz[i] - tr_xyz_zz[i];

        tr_0_0_x_xyz_yyy[i] = 2.0 * tr_xxyz_yyy[i] * tbe_0 + 2.0 * tr_xyz_xyyy[i] * tke_0 - tr_yz_yyy[i];

        tr_0_0_x_xyz_yyz[i] = 2.0 * tr_xxyz_yyz[i] * tbe_0 + 2.0 * tr_xyz_xyyz[i] * tke_0 - tr_yz_yyz[i];

        tr_0_0_x_xyz_yzz[i] = 2.0 * tr_xxyz_yzz[i] * tbe_0 + 2.0 * tr_xyz_xyzz[i] * tke_0 - tr_yz_yzz[i];

        tr_0_0_x_xyz_zzz[i] = 2.0 * tr_xxyz_zzz[i] * tbe_0 + 2.0 * tr_xyz_xzzz[i] * tke_0 - tr_yz_zzz[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto tr_0_0_x_xzz_xxx = pbuffer.data(idx_op_geom_010_ff + 50);

    auto tr_0_0_x_xzz_xxy = pbuffer.data(idx_op_geom_010_ff + 51);

    auto tr_0_0_x_xzz_xxz = pbuffer.data(idx_op_geom_010_ff + 52);

    auto tr_0_0_x_xzz_xyy = pbuffer.data(idx_op_geom_010_ff + 53);

    auto tr_0_0_x_xzz_xyz = pbuffer.data(idx_op_geom_010_ff + 54);

    auto tr_0_0_x_xzz_xzz = pbuffer.data(idx_op_geom_010_ff + 55);

    auto tr_0_0_x_xzz_yyy = pbuffer.data(idx_op_geom_010_ff + 56);

    auto tr_0_0_x_xzz_yyz = pbuffer.data(idx_op_geom_010_ff + 57);

    auto tr_0_0_x_xzz_yzz = pbuffer.data(idx_op_geom_010_ff + 58);

    auto tr_0_0_x_xzz_zzz = pbuffer.data(idx_op_geom_010_ff + 59);

    #pragma omp simd aligned(tr_0_0_x_xzz_xxx, tr_0_0_x_xzz_xxy, tr_0_0_x_xzz_xxz, tr_0_0_x_xzz_xyy, tr_0_0_x_xzz_xyz, tr_0_0_x_xzz_xzz, tr_0_0_x_xzz_yyy, tr_0_0_x_xzz_yyz, tr_0_0_x_xzz_yzz, tr_0_0_x_xzz_zzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, tr_xzz_xx, tr_xzz_xxxx, tr_xzz_xxxy, tr_xzz_xxxz, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzz_xxx[i] = 2.0 * tr_xxzz_xxx[i] * tbe_0 + 2.0 * tr_xzz_xxxx[i] * tke_0 - tr_zz_xxx[i] - 3.0 * tr_xzz_xx[i];

        tr_0_0_x_xzz_xxy[i] = 2.0 * tr_xxzz_xxy[i] * tbe_0 + 2.0 * tr_xzz_xxxy[i] * tke_0 - tr_zz_xxy[i] - 2.0 * tr_xzz_xy[i];

        tr_0_0_x_xzz_xxz[i] = 2.0 * tr_xxzz_xxz[i] * tbe_0 + 2.0 * tr_xzz_xxxz[i] * tke_0 - tr_zz_xxz[i] - 2.0 * tr_xzz_xz[i];

        tr_0_0_x_xzz_xyy[i] = 2.0 * tr_xxzz_xyy[i] * tbe_0 + 2.0 * tr_xzz_xxyy[i] * tke_0 - tr_zz_xyy[i] - tr_xzz_yy[i];

        tr_0_0_x_xzz_xyz[i] = 2.0 * tr_xxzz_xyz[i] * tbe_0 + 2.0 * tr_xzz_xxyz[i] * tke_0 - tr_zz_xyz[i] - tr_xzz_yz[i];

        tr_0_0_x_xzz_xzz[i] = 2.0 * tr_xxzz_xzz[i] * tbe_0 + 2.0 * tr_xzz_xxzz[i] * tke_0 - tr_zz_xzz[i] - tr_xzz_zz[i];

        tr_0_0_x_xzz_yyy[i] = 2.0 * tr_xxzz_yyy[i] * tbe_0 + 2.0 * tr_xzz_xyyy[i] * tke_0 - tr_zz_yyy[i];

        tr_0_0_x_xzz_yyz[i] = 2.0 * tr_xxzz_yyz[i] * tbe_0 + 2.0 * tr_xzz_xyyz[i] * tke_0 - tr_zz_yyz[i];

        tr_0_0_x_xzz_yzz[i] = 2.0 * tr_xxzz_yzz[i] * tbe_0 + 2.0 * tr_xzz_xyzz[i] * tke_0 - tr_zz_yzz[i];

        tr_0_0_x_xzz_zzz[i] = 2.0 * tr_xxzz_zzz[i] * tbe_0 + 2.0 * tr_xzz_xzzz[i] * tke_0 - tr_zz_zzz[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto tr_0_0_x_yyy_xxx = pbuffer.data(idx_op_geom_010_ff + 60);

    auto tr_0_0_x_yyy_xxy = pbuffer.data(idx_op_geom_010_ff + 61);

    auto tr_0_0_x_yyy_xxz = pbuffer.data(idx_op_geom_010_ff + 62);

    auto tr_0_0_x_yyy_xyy = pbuffer.data(idx_op_geom_010_ff + 63);

    auto tr_0_0_x_yyy_xyz = pbuffer.data(idx_op_geom_010_ff + 64);

    auto tr_0_0_x_yyy_xzz = pbuffer.data(idx_op_geom_010_ff + 65);

    auto tr_0_0_x_yyy_yyy = pbuffer.data(idx_op_geom_010_ff + 66);

    auto tr_0_0_x_yyy_yyz = pbuffer.data(idx_op_geom_010_ff + 67);

    auto tr_0_0_x_yyy_yzz = pbuffer.data(idx_op_geom_010_ff + 68);

    auto tr_0_0_x_yyy_zzz = pbuffer.data(idx_op_geom_010_ff + 69);

    #pragma omp simd aligned(tr_0_0_x_yyy_xxx, tr_0_0_x_yyy_xxy, tr_0_0_x_yyy_xxz, tr_0_0_x_yyy_xyy, tr_0_0_x_yyy_xyz, tr_0_0_x_yyy_xzz, tr_0_0_x_yyy_yyy, tr_0_0_x_yyy_yyz, tr_0_0_x_yyy_yzz, tr_0_0_x_yyy_zzz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, tr_yyy_xx, tr_yyy_xxxx, tr_yyy_xxxy, tr_yyy_xxxz, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyy_xxx[i] = 2.0 * tr_xyyy_xxx[i] * tbe_0 + 2.0 * tr_yyy_xxxx[i] * tke_0 - 3.0 * tr_yyy_xx[i];

        tr_0_0_x_yyy_xxy[i] = 2.0 * tr_xyyy_xxy[i] * tbe_0 + 2.0 * tr_yyy_xxxy[i] * tke_0 - 2.0 * tr_yyy_xy[i];

        tr_0_0_x_yyy_xxz[i] = 2.0 * tr_xyyy_xxz[i] * tbe_0 + 2.0 * tr_yyy_xxxz[i] * tke_0 - 2.0 * tr_yyy_xz[i];

        tr_0_0_x_yyy_xyy[i] = 2.0 * tr_xyyy_xyy[i] * tbe_0 + 2.0 * tr_yyy_xxyy[i] * tke_0 - tr_yyy_yy[i];

        tr_0_0_x_yyy_xyz[i] = 2.0 * tr_xyyy_xyz[i] * tbe_0 + 2.0 * tr_yyy_xxyz[i] * tke_0 - tr_yyy_yz[i];

        tr_0_0_x_yyy_xzz[i] = 2.0 * tr_xyyy_xzz[i] * tbe_0 + 2.0 * tr_yyy_xxzz[i] * tke_0 - tr_yyy_zz[i];

        tr_0_0_x_yyy_yyy[i] = 2.0 * tr_xyyy_yyy[i] * tbe_0 + 2.0 * tr_yyy_xyyy[i] * tke_0;

        tr_0_0_x_yyy_yyz[i] = 2.0 * tr_xyyy_yyz[i] * tbe_0 + 2.0 * tr_yyy_xyyz[i] * tke_0;

        tr_0_0_x_yyy_yzz[i] = 2.0 * tr_xyyy_yzz[i] * tbe_0 + 2.0 * tr_yyy_xyzz[i] * tke_0;

        tr_0_0_x_yyy_zzz[i] = 2.0 * tr_xyyy_zzz[i] * tbe_0 + 2.0 * tr_yyy_xzzz[i] * tke_0;
    }

    // Set up 70-80 components of targeted buffer : FF

    auto tr_0_0_x_yyz_xxx = pbuffer.data(idx_op_geom_010_ff + 70);

    auto tr_0_0_x_yyz_xxy = pbuffer.data(idx_op_geom_010_ff + 71);

    auto tr_0_0_x_yyz_xxz = pbuffer.data(idx_op_geom_010_ff + 72);

    auto tr_0_0_x_yyz_xyy = pbuffer.data(idx_op_geom_010_ff + 73);

    auto tr_0_0_x_yyz_xyz = pbuffer.data(idx_op_geom_010_ff + 74);

    auto tr_0_0_x_yyz_xzz = pbuffer.data(idx_op_geom_010_ff + 75);

    auto tr_0_0_x_yyz_yyy = pbuffer.data(idx_op_geom_010_ff + 76);

    auto tr_0_0_x_yyz_yyz = pbuffer.data(idx_op_geom_010_ff + 77);

    auto tr_0_0_x_yyz_yzz = pbuffer.data(idx_op_geom_010_ff + 78);

    auto tr_0_0_x_yyz_zzz = pbuffer.data(idx_op_geom_010_ff + 79);

    #pragma omp simd aligned(tr_0_0_x_yyz_xxx, tr_0_0_x_yyz_xxy, tr_0_0_x_yyz_xxz, tr_0_0_x_yyz_xyy, tr_0_0_x_yyz_xyz, tr_0_0_x_yyz_xzz, tr_0_0_x_yyz_yyy, tr_0_0_x_yyz_yyz, tr_0_0_x_yyz_yzz, tr_0_0_x_yyz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_yyz_xx, tr_yyz_xxxx, tr_yyz_xxxy, tr_yyz_xxxz, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyz_xxx[i] = 2.0 * tr_xyyz_xxx[i] * tbe_0 + 2.0 * tr_yyz_xxxx[i] * tke_0 - 3.0 * tr_yyz_xx[i];

        tr_0_0_x_yyz_xxy[i] = 2.0 * tr_xyyz_xxy[i] * tbe_0 + 2.0 * tr_yyz_xxxy[i] * tke_0 - 2.0 * tr_yyz_xy[i];

        tr_0_0_x_yyz_xxz[i] = 2.0 * tr_xyyz_xxz[i] * tbe_0 + 2.0 * tr_yyz_xxxz[i] * tke_0 - 2.0 * tr_yyz_xz[i];

        tr_0_0_x_yyz_xyy[i] = 2.0 * tr_xyyz_xyy[i] * tbe_0 + 2.0 * tr_yyz_xxyy[i] * tke_0 - tr_yyz_yy[i];

        tr_0_0_x_yyz_xyz[i] = 2.0 * tr_xyyz_xyz[i] * tbe_0 + 2.0 * tr_yyz_xxyz[i] * tke_0 - tr_yyz_yz[i];

        tr_0_0_x_yyz_xzz[i] = 2.0 * tr_xyyz_xzz[i] * tbe_0 + 2.0 * tr_yyz_xxzz[i] * tke_0 - tr_yyz_zz[i];

        tr_0_0_x_yyz_yyy[i] = 2.0 * tr_xyyz_yyy[i] * tbe_0 + 2.0 * tr_yyz_xyyy[i] * tke_0;

        tr_0_0_x_yyz_yyz[i] = 2.0 * tr_xyyz_yyz[i] * tbe_0 + 2.0 * tr_yyz_xyyz[i] * tke_0;

        tr_0_0_x_yyz_yzz[i] = 2.0 * tr_xyyz_yzz[i] * tbe_0 + 2.0 * tr_yyz_xyzz[i] * tke_0;

        tr_0_0_x_yyz_zzz[i] = 2.0 * tr_xyyz_zzz[i] * tbe_0 + 2.0 * tr_yyz_xzzz[i] * tke_0;
    }

    // Set up 80-90 components of targeted buffer : FF

    auto tr_0_0_x_yzz_xxx = pbuffer.data(idx_op_geom_010_ff + 80);

    auto tr_0_0_x_yzz_xxy = pbuffer.data(idx_op_geom_010_ff + 81);

    auto tr_0_0_x_yzz_xxz = pbuffer.data(idx_op_geom_010_ff + 82);

    auto tr_0_0_x_yzz_xyy = pbuffer.data(idx_op_geom_010_ff + 83);

    auto tr_0_0_x_yzz_xyz = pbuffer.data(idx_op_geom_010_ff + 84);

    auto tr_0_0_x_yzz_xzz = pbuffer.data(idx_op_geom_010_ff + 85);

    auto tr_0_0_x_yzz_yyy = pbuffer.data(idx_op_geom_010_ff + 86);

    auto tr_0_0_x_yzz_yyz = pbuffer.data(idx_op_geom_010_ff + 87);

    auto tr_0_0_x_yzz_yzz = pbuffer.data(idx_op_geom_010_ff + 88);

    auto tr_0_0_x_yzz_zzz = pbuffer.data(idx_op_geom_010_ff + 89);

    #pragma omp simd aligned(tr_0_0_x_yzz_xxx, tr_0_0_x_yzz_xxy, tr_0_0_x_yzz_xxz, tr_0_0_x_yzz_xyy, tr_0_0_x_yzz_xyz, tr_0_0_x_yzz_xzz, tr_0_0_x_yzz_yyy, tr_0_0_x_yzz_yyz, tr_0_0_x_yzz_yzz, tr_0_0_x_yzz_zzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_yzz_xx, tr_yzz_xxxx, tr_yzz_xxxy, tr_yzz_xxxz, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzz_xxx[i] = 2.0 * tr_xyzz_xxx[i] * tbe_0 + 2.0 * tr_yzz_xxxx[i] * tke_0 - 3.0 * tr_yzz_xx[i];

        tr_0_0_x_yzz_xxy[i] = 2.0 * tr_xyzz_xxy[i] * tbe_0 + 2.0 * tr_yzz_xxxy[i] * tke_0 - 2.0 * tr_yzz_xy[i];

        tr_0_0_x_yzz_xxz[i] = 2.0 * tr_xyzz_xxz[i] * tbe_0 + 2.0 * tr_yzz_xxxz[i] * tke_0 - 2.0 * tr_yzz_xz[i];

        tr_0_0_x_yzz_xyy[i] = 2.0 * tr_xyzz_xyy[i] * tbe_0 + 2.0 * tr_yzz_xxyy[i] * tke_0 - tr_yzz_yy[i];

        tr_0_0_x_yzz_xyz[i] = 2.0 * tr_xyzz_xyz[i] * tbe_0 + 2.0 * tr_yzz_xxyz[i] * tke_0 - tr_yzz_yz[i];

        tr_0_0_x_yzz_xzz[i] = 2.0 * tr_xyzz_xzz[i] * tbe_0 + 2.0 * tr_yzz_xxzz[i] * tke_0 - tr_yzz_zz[i];

        tr_0_0_x_yzz_yyy[i] = 2.0 * tr_xyzz_yyy[i] * tbe_0 + 2.0 * tr_yzz_xyyy[i] * tke_0;

        tr_0_0_x_yzz_yyz[i] = 2.0 * tr_xyzz_yyz[i] * tbe_0 + 2.0 * tr_yzz_xyyz[i] * tke_0;

        tr_0_0_x_yzz_yzz[i] = 2.0 * tr_xyzz_yzz[i] * tbe_0 + 2.0 * tr_yzz_xyzz[i] * tke_0;

        tr_0_0_x_yzz_zzz[i] = 2.0 * tr_xyzz_zzz[i] * tbe_0 + 2.0 * tr_yzz_xzzz[i] * tke_0;
    }

    // Set up 90-100 components of targeted buffer : FF

    auto tr_0_0_x_zzz_xxx = pbuffer.data(idx_op_geom_010_ff + 90);

    auto tr_0_0_x_zzz_xxy = pbuffer.data(idx_op_geom_010_ff + 91);

    auto tr_0_0_x_zzz_xxz = pbuffer.data(idx_op_geom_010_ff + 92);

    auto tr_0_0_x_zzz_xyy = pbuffer.data(idx_op_geom_010_ff + 93);

    auto tr_0_0_x_zzz_xyz = pbuffer.data(idx_op_geom_010_ff + 94);

    auto tr_0_0_x_zzz_xzz = pbuffer.data(idx_op_geom_010_ff + 95);

    auto tr_0_0_x_zzz_yyy = pbuffer.data(idx_op_geom_010_ff + 96);

    auto tr_0_0_x_zzz_yyz = pbuffer.data(idx_op_geom_010_ff + 97);

    auto tr_0_0_x_zzz_yzz = pbuffer.data(idx_op_geom_010_ff + 98);

    auto tr_0_0_x_zzz_zzz = pbuffer.data(idx_op_geom_010_ff + 99);

    #pragma omp simd aligned(tr_0_0_x_zzz_xxx, tr_0_0_x_zzz_xxy, tr_0_0_x_zzz_xxz, tr_0_0_x_zzz_xyy, tr_0_0_x_zzz_xyz, tr_0_0_x_zzz_xzz, tr_0_0_x_zzz_yyy, tr_0_0_x_zzz_yyz, tr_0_0_x_zzz_yzz, tr_0_0_x_zzz_zzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, tr_zzz_xx, tr_zzz_xxxx, tr_zzz_xxxy, tr_zzz_xxxz, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzz_xxx[i] = 2.0 * tr_xzzz_xxx[i] * tbe_0 + 2.0 * tr_zzz_xxxx[i] * tke_0 - 3.0 * tr_zzz_xx[i];

        tr_0_0_x_zzz_xxy[i] = 2.0 * tr_xzzz_xxy[i] * tbe_0 + 2.0 * tr_zzz_xxxy[i] * tke_0 - 2.0 * tr_zzz_xy[i];

        tr_0_0_x_zzz_xxz[i] = 2.0 * tr_xzzz_xxz[i] * tbe_0 + 2.0 * tr_zzz_xxxz[i] * tke_0 - 2.0 * tr_zzz_xz[i];

        tr_0_0_x_zzz_xyy[i] = 2.0 * tr_xzzz_xyy[i] * tbe_0 + 2.0 * tr_zzz_xxyy[i] * tke_0 - tr_zzz_yy[i];

        tr_0_0_x_zzz_xyz[i] = 2.0 * tr_xzzz_xyz[i] * tbe_0 + 2.0 * tr_zzz_xxyz[i] * tke_0 - tr_zzz_yz[i];

        tr_0_0_x_zzz_xzz[i] = 2.0 * tr_xzzz_xzz[i] * tbe_0 + 2.0 * tr_zzz_xxzz[i] * tke_0 - tr_zzz_zz[i];

        tr_0_0_x_zzz_yyy[i] = 2.0 * tr_xzzz_yyy[i] * tbe_0 + 2.0 * tr_zzz_xyyy[i] * tke_0;

        tr_0_0_x_zzz_yyz[i] = 2.0 * tr_xzzz_yyz[i] * tbe_0 + 2.0 * tr_zzz_xyyz[i] * tke_0;

        tr_0_0_x_zzz_yzz[i] = 2.0 * tr_xzzz_yzz[i] * tbe_0 + 2.0 * tr_zzz_xyzz[i] * tke_0;

        tr_0_0_x_zzz_zzz[i] = 2.0 * tr_xzzz_zzz[i] * tbe_0 + 2.0 * tr_zzz_xzzz[i] * tke_0;
    }

    // Set up 100-110 components of targeted buffer : FF

    auto tr_0_0_y_xxx_xxx = pbuffer.data(idx_op_geom_010_ff + 100);

    auto tr_0_0_y_xxx_xxy = pbuffer.data(idx_op_geom_010_ff + 101);

    auto tr_0_0_y_xxx_xxz = pbuffer.data(idx_op_geom_010_ff + 102);

    auto tr_0_0_y_xxx_xyy = pbuffer.data(idx_op_geom_010_ff + 103);

    auto tr_0_0_y_xxx_xyz = pbuffer.data(idx_op_geom_010_ff + 104);

    auto tr_0_0_y_xxx_xzz = pbuffer.data(idx_op_geom_010_ff + 105);

    auto tr_0_0_y_xxx_yyy = pbuffer.data(idx_op_geom_010_ff + 106);

    auto tr_0_0_y_xxx_yyz = pbuffer.data(idx_op_geom_010_ff + 107);

    auto tr_0_0_y_xxx_yzz = pbuffer.data(idx_op_geom_010_ff + 108);

    auto tr_0_0_y_xxx_zzz = pbuffer.data(idx_op_geom_010_ff + 109);

    #pragma omp simd aligned(tr_0_0_y_xxx_xxx, tr_0_0_y_xxx_xxy, tr_0_0_y_xxx_xxz, tr_0_0_y_xxx_xyy, tr_0_0_y_xxx_xyz, tr_0_0_y_xxx_xzz, tr_0_0_y_xxx_yyy, tr_0_0_y_xxx_yyz, tr_0_0_y_xxx_yzz, tr_0_0_y_xxx_zzz, tr_xxx_xx, tr_xxx_xxxy, tr_xxx_xxyy, tr_xxx_xxyz, tr_xxx_xy, tr_xxx_xyyy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_yy, tr_xxx_yyyy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxx_xxx[i] = 2.0 * tr_xxxy_xxx[i] * tbe_0 + 2.0 * tr_xxx_xxxy[i] * tke_0;

        tr_0_0_y_xxx_xxy[i] = 2.0 * tr_xxxy_xxy[i] * tbe_0 + 2.0 * tr_xxx_xxyy[i] * tke_0 - tr_xxx_xx[i];

        tr_0_0_y_xxx_xxz[i] = 2.0 * tr_xxxy_xxz[i] * tbe_0 + 2.0 * tr_xxx_xxyz[i] * tke_0;

        tr_0_0_y_xxx_xyy[i] = 2.0 * tr_xxxy_xyy[i] * tbe_0 + 2.0 * tr_xxx_xyyy[i] * tke_0 - 2.0 * tr_xxx_xy[i];

        tr_0_0_y_xxx_xyz[i] = 2.0 * tr_xxxy_xyz[i] * tbe_0 + 2.0 * tr_xxx_xyyz[i] * tke_0 - tr_xxx_xz[i];

        tr_0_0_y_xxx_xzz[i] = 2.0 * tr_xxxy_xzz[i] * tbe_0 + 2.0 * tr_xxx_xyzz[i] * tke_0;

        tr_0_0_y_xxx_yyy[i] = 2.0 * tr_xxxy_yyy[i] * tbe_0 + 2.0 * tr_xxx_yyyy[i] * tke_0 - 3.0 * tr_xxx_yy[i];

        tr_0_0_y_xxx_yyz[i] = 2.0 * tr_xxxy_yyz[i] * tbe_0 + 2.0 * tr_xxx_yyyz[i] * tke_0 - 2.0 * tr_xxx_yz[i];

        tr_0_0_y_xxx_yzz[i] = 2.0 * tr_xxxy_yzz[i] * tbe_0 + 2.0 * tr_xxx_yyzz[i] * tke_0 - tr_xxx_zz[i];

        tr_0_0_y_xxx_zzz[i] = 2.0 * tr_xxxy_zzz[i] * tbe_0 + 2.0 * tr_xxx_yzzz[i] * tke_0;
    }

    // Set up 110-120 components of targeted buffer : FF

    auto tr_0_0_y_xxy_xxx = pbuffer.data(idx_op_geom_010_ff + 110);

    auto tr_0_0_y_xxy_xxy = pbuffer.data(idx_op_geom_010_ff + 111);

    auto tr_0_0_y_xxy_xxz = pbuffer.data(idx_op_geom_010_ff + 112);

    auto tr_0_0_y_xxy_xyy = pbuffer.data(idx_op_geom_010_ff + 113);

    auto tr_0_0_y_xxy_xyz = pbuffer.data(idx_op_geom_010_ff + 114);

    auto tr_0_0_y_xxy_xzz = pbuffer.data(idx_op_geom_010_ff + 115);

    auto tr_0_0_y_xxy_yyy = pbuffer.data(idx_op_geom_010_ff + 116);

    auto tr_0_0_y_xxy_yyz = pbuffer.data(idx_op_geom_010_ff + 117);

    auto tr_0_0_y_xxy_yzz = pbuffer.data(idx_op_geom_010_ff + 118);

    auto tr_0_0_y_xxy_zzz = pbuffer.data(idx_op_geom_010_ff + 119);

    #pragma omp simd aligned(tr_0_0_y_xxy_xxx, tr_0_0_y_xxy_xxy, tr_0_0_y_xxy_xxz, tr_0_0_y_xxy_xyy, tr_0_0_y_xxy_xyz, tr_0_0_y_xxy_xzz, tr_0_0_y_xxy_yyy, tr_0_0_y_xxy_yyz, tr_0_0_y_xxy_yzz, tr_0_0_y_xxy_zzz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxy_xx, tr_xxy_xxxy, tr_xxy_xxyy, tr_xxy_xxyz, tr_xxy_xy, tr_xxy_xyyy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_yy, tr_xxy_yyyy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxy_xxx[i] = 2.0 * tr_xxyy_xxx[i] * tbe_0 + 2.0 * tr_xxy_xxxy[i] * tke_0 - tr_xx_xxx[i];

        tr_0_0_y_xxy_xxy[i] = 2.0 * tr_xxyy_xxy[i] * tbe_0 + 2.0 * tr_xxy_xxyy[i] * tke_0 - tr_xx_xxy[i] - tr_xxy_xx[i];

        tr_0_0_y_xxy_xxz[i] = 2.0 * tr_xxyy_xxz[i] * tbe_0 + 2.0 * tr_xxy_xxyz[i] * tke_0 - tr_xx_xxz[i];

        tr_0_0_y_xxy_xyy[i] = 2.0 * tr_xxyy_xyy[i] * tbe_0 + 2.0 * tr_xxy_xyyy[i] * tke_0 - tr_xx_xyy[i] - 2.0 * tr_xxy_xy[i];

        tr_0_0_y_xxy_xyz[i] = 2.0 * tr_xxyy_xyz[i] * tbe_0 + 2.0 * tr_xxy_xyyz[i] * tke_0 - tr_xx_xyz[i] - tr_xxy_xz[i];

        tr_0_0_y_xxy_xzz[i] = 2.0 * tr_xxyy_xzz[i] * tbe_0 + 2.0 * tr_xxy_xyzz[i] * tke_0 - tr_xx_xzz[i];

        tr_0_0_y_xxy_yyy[i] = 2.0 * tr_xxyy_yyy[i] * tbe_0 + 2.0 * tr_xxy_yyyy[i] * tke_0 - tr_xx_yyy[i] - 3.0 * tr_xxy_yy[i];

        tr_0_0_y_xxy_yyz[i] = 2.0 * tr_xxyy_yyz[i] * tbe_0 + 2.0 * tr_xxy_yyyz[i] * tke_0 - tr_xx_yyz[i] - 2.0 * tr_xxy_yz[i];

        tr_0_0_y_xxy_yzz[i] = 2.0 * tr_xxyy_yzz[i] * tbe_0 + 2.0 * tr_xxy_yyzz[i] * tke_0 - tr_xx_yzz[i] - tr_xxy_zz[i];

        tr_0_0_y_xxy_zzz[i] = 2.0 * tr_xxyy_zzz[i] * tbe_0 + 2.0 * tr_xxy_yzzz[i] * tke_0 - tr_xx_zzz[i];
    }

    // Set up 120-130 components of targeted buffer : FF

    auto tr_0_0_y_xxz_xxx = pbuffer.data(idx_op_geom_010_ff + 120);

    auto tr_0_0_y_xxz_xxy = pbuffer.data(idx_op_geom_010_ff + 121);

    auto tr_0_0_y_xxz_xxz = pbuffer.data(idx_op_geom_010_ff + 122);

    auto tr_0_0_y_xxz_xyy = pbuffer.data(idx_op_geom_010_ff + 123);

    auto tr_0_0_y_xxz_xyz = pbuffer.data(idx_op_geom_010_ff + 124);

    auto tr_0_0_y_xxz_xzz = pbuffer.data(idx_op_geom_010_ff + 125);

    auto tr_0_0_y_xxz_yyy = pbuffer.data(idx_op_geom_010_ff + 126);

    auto tr_0_0_y_xxz_yyz = pbuffer.data(idx_op_geom_010_ff + 127);

    auto tr_0_0_y_xxz_yzz = pbuffer.data(idx_op_geom_010_ff + 128);

    auto tr_0_0_y_xxz_zzz = pbuffer.data(idx_op_geom_010_ff + 129);

    #pragma omp simd aligned(tr_0_0_y_xxz_xxx, tr_0_0_y_xxz_xxy, tr_0_0_y_xxz_xxz, tr_0_0_y_xxz_xyy, tr_0_0_y_xxz_xyz, tr_0_0_y_xxz_xzz, tr_0_0_y_xxz_yyy, tr_0_0_y_xxz_yyz, tr_0_0_y_xxz_yzz, tr_0_0_y_xxz_zzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, tr_xxz_xx, tr_xxz_xxxy, tr_xxz_xxyy, tr_xxz_xxyz, tr_xxz_xy, tr_xxz_xyyy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_yy, tr_xxz_yyyy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxz_xxx[i] = 2.0 * tr_xxyz_xxx[i] * tbe_0 + 2.0 * tr_xxz_xxxy[i] * tke_0;

        tr_0_0_y_xxz_xxy[i] = 2.0 * tr_xxyz_xxy[i] * tbe_0 + 2.0 * tr_xxz_xxyy[i] * tke_0 - tr_xxz_xx[i];

        tr_0_0_y_xxz_xxz[i] = 2.0 * tr_xxyz_xxz[i] * tbe_0 + 2.0 * tr_xxz_xxyz[i] * tke_0;

        tr_0_0_y_xxz_xyy[i] = 2.0 * tr_xxyz_xyy[i] * tbe_0 + 2.0 * tr_xxz_xyyy[i] * tke_0 - 2.0 * tr_xxz_xy[i];

        tr_0_0_y_xxz_xyz[i] = 2.0 * tr_xxyz_xyz[i] * tbe_0 + 2.0 * tr_xxz_xyyz[i] * tke_0 - tr_xxz_xz[i];

        tr_0_0_y_xxz_xzz[i] = 2.0 * tr_xxyz_xzz[i] * tbe_0 + 2.0 * tr_xxz_xyzz[i] * tke_0;

        tr_0_0_y_xxz_yyy[i] = 2.0 * tr_xxyz_yyy[i] * tbe_0 + 2.0 * tr_xxz_yyyy[i] * tke_0 - 3.0 * tr_xxz_yy[i];

        tr_0_0_y_xxz_yyz[i] = 2.0 * tr_xxyz_yyz[i] * tbe_0 + 2.0 * tr_xxz_yyyz[i] * tke_0 - 2.0 * tr_xxz_yz[i];

        tr_0_0_y_xxz_yzz[i] = 2.0 * tr_xxyz_yzz[i] * tbe_0 + 2.0 * tr_xxz_yyzz[i] * tke_0 - tr_xxz_zz[i];

        tr_0_0_y_xxz_zzz[i] = 2.0 * tr_xxyz_zzz[i] * tbe_0 + 2.0 * tr_xxz_yzzz[i] * tke_0;
    }

    // Set up 130-140 components of targeted buffer : FF

    auto tr_0_0_y_xyy_xxx = pbuffer.data(idx_op_geom_010_ff + 130);

    auto tr_0_0_y_xyy_xxy = pbuffer.data(idx_op_geom_010_ff + 131);

    auto tr_0_0_y_xyy_xxz = pbuffer.data(idx_op_geom_010_ff + 132);

    auto tr_0_0_y_xyy_xyy = pbuffer.data(idx_op_geom_010_ff + 133);

    auto tr_0_0_y_xyy_xyz = pbuffer.data(idx_op_geom_010_ff + 134);

    auto tr_0_0_y_xyy_xzz = pbuffer.data(idx_op_geom_010_ff + 135);

    auto tr_0_0_y_xyy_yyy = pbuffer.data(idx_op_geom_010_ff + 136);

    auto tr_0_0_y_xyy_yyz = pbuffer.data(idx_op_geom_010_ff + 137);

    auto tr_0_0_y_xyy_yzz = pbuffer.data(idx_op_geom_010_ff + 138);

    auto tr_0_0_y_xyy_zzz = pbuffer.data(idx_op_geom_010_ff + 139);

    #pragma omp simd aligned(tr_0_0_y_xyy_xxx, tr_0_0_y_xyy_xxy, tr_0_0_y_xyy_xxz, tr_0_0_y_xyy_xyy, tr_0_0_y_xyy_xyz, tr_0_0_y_xyy_xzz, tr_0_0_y_xyy_yyy, tr_0_0_y_xyy_yyz, tr_0_0_y_xyy_yzz, tr_0_0_y_xyy_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyy_xx, tr_xyy_xxxy, tr_xyy_xxyy, tr_xyy_xxyz, tr_xyy_xy, tr_xyy_xyyy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_yy, tr_xyy_yyyy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyy_xxx[i] = 2.0 * tr_xyyy_xxx[i] * tbe_0 + 2.0 * tr_xyy_xxxy[i] * tke_0 - 2.0 * tr_xy_xxx[i];

        tr_0_0_y_xyy_xxy[i] = 2.0 * tr_xyyy_xxy[i] * tbe_0 + 2.0 * tr_xyy_xxyy[i] * tke_0 - 2.0 * tr_xy_xxy[i] - tr_xyy_xx[i];

        tr_0_0_y_xyy_xxz[i] = 2.0 * tr_xyyy_xxz[i] * tbe_0 + 2.0 * tr_xyy_xxyz[i] * tke_0 - 2.0 * tr_xy_xxz[i];

        tr_0_0_y_xyy_xyy[i] = 2.0 * tr_xyyy_xyy[i] * tbe_0 + 2.0 * tr_xyy_xyyy[i] * tke_0 - 2.0 * tr_xy_xyy[i] - 2.0 * tr_xyy_xy[i];

        tr_0_0_y_xyy_xyz[i] = 2.0 * tr_xyyy_xyz[i] * tbe_0 + 2.0 * tr_xyy_xyyz[i] * tke_0 - 2.0 * tr_xy_xyz[i] - tr_xyy_xz[i];

        tr_0_0_y_xyy_xzz[i] = 2.0 * tr_xyyy_xzz[i] * tbe_0 + 2.0 * tr_xyy_xyzz[i] * tke_0 - 2.0 * tr_xy_xzz[i];

        tr_0_0_y_xyy_yyy[i] = 2.0 * tr_xyyy_yyy[i] * tbe_0 + 2.0 * tr_xyy_yyyy[i] * tke_0 - 2.0 * tr_xy_yyy[i] - 3.0 * tr_xyy_yy[i];

        tr_0_0_y_xyy_yyz[i] = 2.0 * tr_xyyy_yyz[i] * tbe_0 + 2.0 * tr_xyy_yyyz[i] * tke_0 - 2.0 * tr_xy_yyz[i] - 2.0 * tr_xyy_yz[i];

        tr_0_0_y_xyy_yzz[i] = 2.0 * tr_xyyy_yzz[i] * tbe_0 + 2.0 * tr_xyy_yyzz[i] * tke_0 - 2.0 * tr_xy_yzz[i] - tr_xyy_zz[i];

        tr_0_0_y_xyy_zzz[i] = 2.0 * tr_xyyy_zzz[i] * tbe_0 + 2.0 * tr_xyy_yzzz[i] * tke_0 - 2.0 * tr_xy_zzz[i];
    }

    // Set up 140-150 components of targeted buffer : FF

    auto tr_0_0_y_xyz_xxx = pbuffer.data(idx_op_geom_010_ff + 140);

    auto tr_0_0_y_xyz_xxy = pbuffer.data(idx_op_geom_010_ff + 141);

    auto tr_0_0_y_xyz_xxz = pbuffer.data(idx_op_geom_010_ff + 142);

    auto tr_0_0_y_xyz_xyy = pbuffer.data(idx_op_geom_010_ff + 143);

    auto tr_0_0_y_xyz_xyz = pbuffer.data(idx_op_geom_010_ff + 144);

    auto tr_0_0_y_xyz_xzz = pbuffer.data(idx_op_geom_010_ff + 145);

    auto tr_0_0_y_xyz_yyy = pbuffer.data(idx_op_geom_010_ff + 146);

    auto tr_0_0_y_xyz_yyz = pbuffer.data(idx_op_geom_010_ff + 147);

    auto tr_0_0_y_xyz_yzz = pbuffer.data(idx_op_geom_010_ff + 148);

    auto tr_0_0_y_xyz_zzz = pbuffer.data(idx_op_geom_010_ff + 149);

    #pragma omp simd aligned(tr_0_0_y_xyz_xxx, tr_0_0_y_xyz_xxy, tr_0_0_y_xyz_xxz, tr_0_0_y_xyz_xyy, tr_0_0_y_xyz_xyz, tr_0_0_y_xyz_xzz, tr_0_0_y_xyz_yyy, tr_0_0_y_xyz_yyz, tr_0_0_y_xyz_yzz, tr_0_0_y_xyz_zzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, tr_xyz_xx, tr_xyz_xxxy, tr_xyz_xxyy, tr_xyz_xxyz, tr_xyz_xy, tr_xyz_xyyy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_yy, tr_xyz_yyyy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyz_xxx[i] = 2.0 * tr_xyyz_xxx[i] * tbe_0 + 2.0 * tr_xyz_xxxy[i] * tke_0 - tr_xz_xxx[i];

        tr_0_0_y_xyz_xxy[i] = 2.0 * tr_xyyz_xxy[i] * tbe_0 + 2.0 * tr_xyz_xxyy[i] * tke_0 - tr_xz_xxy[i] - tr_xyz_xx[i];

        tr_0_0_y_xyz_xxz[i] = 2.0 * tr_xyyz_xxz[i] * tbe_0 + 2.0 * tr_xyz_xxyz[i] * tke_0 - tr_xz_xxz[i];

        tr_0_0_y_xyz_xyy[i] = 2.0 * tr_xyyz_xyy[i] * tbe_0 + 2.0 * tr_xyz_xyyy[i] * tke_0 - tr_xz_xyy[i] - 2.0 * tr_xyz_xy[i];

        tr_0_0_y_xyz_xyz[i] = 2.0 * tr_xyyz_xyz[i] * tbe_0 + 2.0 * tr_xyz_xyyz[i] * tke_0 - tr_xz_xyz[i] - tr_xyz_xz[i];

        tr_0_0_y_xyz_xzz[i] = 2.0 * tr_xyyz_xzz[i] * tbe_0 + 2.0 * tr_xyz_xyzz[i] * tke_0 - tr_xz_xzz[i];

        tr_0_0_y_xyz_yyy[i] = 2.0 * tr_xyyz_yyy[i] * tbe_0 + 2.0 * tr_xyz_yyyy[i] * tke_0 - tr_xz_yyy[i] - 3.0 * tr_xyz_yy[i];

        tr_0_0_y_xyz_yyz[i] = 2.0 * tr_xyyz_yyz[i] * tbe_0 + 2.0 * tr_xyz_yyyz[i] * tke_0 - tr_xz_yyz[i] - 2.0 * tr_xyz_yz[i];

        tr_0_0_y_xyz_yzz[i] = 2.0 * tr_xyyz_yzz[i] * tbe_0 + 2.0 * tr_xyz_yyzz[i] * tke_0 - tr_xz_yzz[i] - tr_xyz_zz[i];

        tr_0_0_y_xyz_zzz[i] = 2.0 * tr_xyyz_zzz[i] * tbe_0 + 2.0 * tr_xyz_yzzz[i] * tke_0 - tr_xz_zzz[i];
    }

    // Set up 150-160 components of targeted buffer : FF

    auto tr_0_0_y_xzz_xxx = pbuffer.data(idx_op_geom_010_ff + 150);

    auto tr_0_0_y_xzz_xxy = pbuffer.data(idx_op_geom_010_ff + 151);

    auto tr_0_0_y_xzz_xxz = pbuffer.data(idx_op_geom_010_ff + 152);

    auto tr_0_0_y_xzz_xyy = pbuffer.data(idx_op_geom_010_ff + 153);

    auto tr_0_0_y_xzz_xyz = pbuffer.data(idx_op_geom_010_ff + 154);

    auto tr_0_0_y_xzz_xzz = pbuffer.data(idx_op_geom_010_ff + 155);

    auto tr_0_0_y_xzz_yyy = pbuffer.data(idx_op_geom_010_ff + 156);

    auto tr_0_0_y_xzz_yyz = pbuffer.data(idx_op_geom_010_ff + 157);

    auto tr_0_0_y_xzz_yzz = pbuffer.data(idx_op_geom_010_ff + 158);

    auto tr_0_0_y_xzz_zzz = pbuffer.data(idx_op_geom_010_ff + 159);

    #pragma omp simd aligned(tr_0_0_y_xzz_xxx, tr_0_0_y_xzz_xxy, tr_0_0_y_xzz_xxz, tr_0_0_y_xzz_xyy, tr_0_0_y_xzz_xyz, tr_0_0_y_xzz_xzz, tr_0_0_y_xzz_yyy, tr_0_0_y_xzz_yyz, tr_0_0_y_xzz_yzz, tr_0_0_y_xzz_zzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, tr_xzz_xx, tr_xzz_xxxy, tr_xzz_xxyy, tr_xzz_xxyz, tr_xzz_xy, tr_xzz_xyyy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_yy, tr_xzz_yyyy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzz_xxx[i] = 2.0 * tr_xyzz_xxx[i] * tbe_0 + 2.0 * tr_xzz_xxxy[i] * tke_0;

        tr_0_0_y_xzz_xxy[i] = 2.0 * tr_xyzz_xxy[i] * tbe_0 + 2.0 * tr_xzz_xxyy[i] * tke_0 - tr_xzz_xx[i];

        tr_0_0_y_xzz_xxz[i] = 2.0 * tr_xyzz_xxz[i] * tbe_0 + 2.0 * tr_xzz_xxyz[i] * tke_0;

        tr_0_0_y_xzz_xyy[i] = 2.0 * tr_xyzz_xyy[i] * tbe_0 + 2.0 * tr_xzz_xyyy[i] * tke_0 - 2.0 * tr_xzz_xy[i];

        tr_0_0_y_xzz_xyz[i] = 2.0 * tr_xyzz_xyz[i] * tbe_0 + 2.0 * tr_xzz_xyyz[i] * tke_0 - tr_xzz_xz[i];

        tr_0_0_y_xzz_xzz[i] = 2.0 * tr_xyzz_xzz[i] * tbe_0 + 2.0 * tr_xzz_xyzz[i] * tke_0;

        tr_0_0_y_xzz_yyy[i] = 2.0 * tr_xyzz_yyy[i] * tbe_0 + 2.0 * tr_xzz_yyyy[i] * tke_0 - 3.0 * tr_xzz_yy[i];

        tr_0_0_y_xzz_yyz[i] = 2.0 * tr_xyzz_yyz[i] * tbe_0 + 2.0 * tr_xzz_yyyz[i] * tke_0 - 2.0 * tr_xzz_yz[i];

        tr_0_0_y_xzz_yzz[i] = 2.0 * tr_xyzz_yzz[i] * tbe_0 + 2.0 * tr_xzz_yyzz[i] * tke_0 - tr_xzz_zz[i];

        tr_0_0_y_xzz_zzz[i] = 2.0 * tr_xyzz_zzz[i] * tbe_0 + 2.0 * tr_xzz_yzzz[i] * tke_0;
    }

    // Set up 160-170 components of targeted buffer : FF

    auto tr_0_0_y_yyy_xxx = pbuffer.data(idx_op_geom_010_ff + 160);

    auto tr_0_0_y_yyy_xxy = pbuffer.data(idx_op_geom_010_ff + 161);

    auto tr_0_0_y_yyy_xxz = pbuffer.data(idx_op_geom_010_ff + 162);

    auto tr_0_0_y_yyy_xyy = pbuffer.data(idx_op_geom_010_ff + 163);

    auto tr_0_0_y_yyy_xyz = pbuffer.data(idx_op_geom_010_ff + 164);

    auto tr_0_0_y_yyy_xzz = pbuffer.data(idx_op_geom_010_ff + 165);

    auto tr_0_0_y_yyy_yyy = pbuffer.data(idx_op_geom_010_ff + 166);

    auto tr_0_0_y_yyy_yyz = pbuffer.data(idx_op_geom_010_ff + 167);

    auto tr_0_0_y_yyy_yzz = pbuffer.data(idx_op_geom_010_ff + 168);

    auto tr_0_0_y_yyy_zzz = pbuffer.data(idx_op_geom_010_ff + 169);

    #pragma omp simd aligned(tr_0_0_y_yyy_xxx, tr_0_0_y_yyy_xxy, tr_0_0_y_yyy_xxz, tr_0_0_y_yyy_xyy, tr_0_0_y_yyy_xyz, tr_0_0_y_yyy_xzz, tr_0_0_y_yyy_yyy, tr_0_0_y_yyy_yyz, tr_0_0_y_yyy_yzz, tr_0_0_y_yyy_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyy_xx, tr_yyy_xxxy, tr_yyy_xxyy, tr_yyy_xxyz, tr_yyy_xy, tr_yyy_xyyy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_yy, tr_yyy_yyyy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyy_xxx[i] = 2.0 * tr_yyyy_xxx[i] * tbe_0 + 2.0 * tr_yyy_xxxy[i] * tke_0 - 3.0 * tr_yy_xxx[i];

        tr_0_0_y_yyy_xxy[i] = 2.0 * tr_yyyy_xxy[i] * tbe_0 + 2.0 * tr_yyy_xxyy[i] * tke_0 - 3.0 * tr_yy_xxy[i] - tr_yyy_xx[i];

        tr_0_0_y_yyy_xxz[i] = 2.0 * tr_yyyy_xxz[i] * tbe_0 + 2.0 * tr_yyy_xxyz[i] * tke_0 - 3.0 * tr_yy_xxz[i];

        tr_0_0_y_yyy_xyy[i] = 2.0 * tr_yyyy_xyy[i] * tbe_0 + 2.0 * tr_yyy_xyyy[i] * tke_0 - 3.0 * tr_yy_xyy[i] - 2.0 * tr_yyy_xy[i];

        tr_0_0_y_yyy_xyz[i] = 2.0 * tr_yyyy_xyz[i] * tbe_0 + 2.0 * tr_yyy_xyyz[i] * tke_0 - 3.0 * tr_yy_xyz[i] - tr_yyy_xz[i];

        tr_0_0_y_yyy_xzz[i] = 2.0 * tr_yyyy_xzz[i] * tbe_0 + 2.0 * tr_yyy_xyzz[i] * tke_0 - 3.0 * tr_yy_xzz[i];

        tr_0_0_y_yyy_yyy[i] = 2.0 * tr_yyyy_yyy[i] * tbe_0 + 2.0 * tr_yyy_yyyy[i] * tke_0 - 3.0 * tr_yy_yyy[i] - 3.0 * tr_yyy_yy[i];

        tr_0_0_y_yyy_yyz[i] = 2.0 * tr_yyyy_yyz[i] * tbe_0 + 2.0 * tr_yyy_yyyz[i] * tke_0 - 3.0 * tr_yy_yyz[i] - 2.0 * tr_yyy_yz[i];

        tr_0_0_y_yyy_yzz[i] = 2.0 * tr_yyyy_yzz[i] * tbe_0 + 2.0 * tr_yyy_yyzz[i] * tke_0 - 3.0 * tr_yy_yzz[i] - tr_yyy_zz[i];

        tr_0_0_y_yyy_zzz[i] = 2.0 * tr_yyyy_zzz[i] * tbe_0 + 2.0 * tr_yyy_yzzz[i] * tke_0 - 3.0 * tr_yy_zzz[i];
    }

    // Set up 170-180 components of targeted buffer : FF

    auto tr_0_0_y_yyz_xxx = pbuffer.data(idx_op_geom_010_ff + 170);

    auto tr_0_0_y_yyz_xxy = pbuffer.data(idx_op_geom_010_ff + 171);

    auto tr_0_0_y_yyz_xxz = pbuffer.data(idx_op_geom_010_ff + 172);

    auto tr_0_0_y_yyz_xyy = pbuffer.data(idx_op_geom_010_ff + 173);

    auto tr_0_0_y_yyz_xyz = pbuffer.data(idx_op_geom_010_ff + 174);

    auto tr_0_0_y_yyz_xzz = pbuffer.data(idx_op_geom_010_ff + 175);

    auto tr_0_0_y_yyz_yyy = pbuffer.data(idx_op_geom_010_ff + 176);

    auto tr_0_0_y_yyz_yyz = pbuffer.data(idx_op_geom_010_ff + 177);

    auto tr_0_0_y_yyz_yzz = pbuffer.data(idx_op_geom_010_ff + 178);

    auto tr_0_0_y_yyz_zzz = pbuffer.data(idx_op_geom_010_ff + 179);

    #pragma omp simd aligned(tr_0_0_y_yyz_xxx, tr_0_0_y_yyz_xxy, tr_0_0_y_yyz_xxz, tr_0_0_y_yyz_xyy, tr_0_0_y_yyz_xyz, tr_0_0_y_yyz_xzz, tr_0_0_y_yyz_yyy, tr_0_0_y_yyz_yyz, tr_0_0_y_yyz_yzz, tr_0_0_y_yyz_zzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, tr_yyz_xx, tr_yyz_xxxy, tr_yyz_xxyy, tr_yyz_xxyz, tr_yyz_xy, tr_yyz_xyyy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_yy, tr_yyz_yyyy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyz_xxx[i] = 2.0 * tr_yyyz_xxx[i] * tbe_0 + 2.0 * tr_yyz_xxxy[i] * tke_0 - 2.0 * tr_yz_xxx[i];

        tr_0_0_y_yyz_xxy[i] = 2.0 * tr_yyyz_xxy[i] * tbe_0 + 2.0 * tr_yyz_xxyy[i] * tke_0 - 2.0 * tr_yz_xxy[i] - tr_yyz_xx[i];

        tr_0_0_y_yyz_xxz[i] = 2.0 * tr_yyyz_xxz[i] * tbe_0 + 2.0 * tr_yyz_xxyz[i] * tke_0 - 2.0 * tr_yz_xxz[i];

        tr_0_0_y_yyz_xyy[i] = 2.0 * tr_yyyz_xyy[i] * tbe_0 + 2.0 * tr_yyz_xyyy[i] * tke_0 - 2.0 * tr_yz_xyy[i] - 2.0 * tr_yyz_xy[i];

        tr_0_0_y_yyz_xyz[i] = 2.0 * tr_yyyz_xyz[i] * tbe_0 + 2.0 * tr_yyz_xyyz[i] * tke_0 - 2.0 * tr_yz_xyz[i] - tr_yyz_xz[i];

        tr_0_0_y_yyz_xzz[i] = 2.0 * tr_yyyz_xzz[i] * tbe_0 + 2.0 * tr_yyz_xyzz[i] * tke_0 - 2.0 * tr_yz_xzz[i];

        tr_0_0_y_yyz_yyy[i] = 2.0 * tr_yyyz_yyy[i] * tbe_0 + 2.0 * tr_yyz_yyyy[i] * tke_0 - 2.0 * tr_yz_yyy[i] - 3.0 * tr_yyz_yy[i];

        tr_0_0_y_yyz_yyz[i] = 2.0 * tr_yyyz_yyz[i] * tbe_0 + 2.0 * tr_yyz_yyyz[i] * tke_0 - 2.0 * tr_yz_yyz[i] - 2.0 * tr_yyz_yz[i];

        tr_0_0_y_yyz_yzz[i] = 2.0 * tr_yyyz_yzz[i] * tbe_0 + 2.0 * tr_yyz_yyzz[i] * tke_0 - 2.0 * tr_yz_yzz[i] - tr_yyz_zz[i];

        tr_0_0_y_yyz_zzz[i] = 2.0 * tr_yyyz_zzz[i] * tbe_0 + 2.0 * tr_yyz_yzzz[i] * tke_0 - 2.0 * tr_yz_zzz[i];
    }

    // Set up 180-190 components of targeted buffer : FF

    auto tr_0_0_y_yzz_xxx = pbuffer.data(idx_op_geom_010_ff + 180);

    auto tr_0_0_y_yzz_xxy = pbuffer.data(idx_op_geom_010_ff + 181);

    auto tr_0_0_y_yzz_xxz = pbuffer.data(idx_op_geom_010_ff + 182);

    auto tr_0_0_y_yzz_xyy = pbuffer.data(idx_op_geom_010_ff + 183);

    auto tr_0_0_y_yzz_xyz = pbuffer.data(idx_op_geom_010_ff + 184);

    auto tr_0_0_y_yzz_xzz = pbuffer.data(idx_op_geom_010_ff + 185);

    auto tr_0_0_y_yzz_yyy = pbuffer.data(idx_op_geom_010_ff + 186);

    auto tr_0_0_y_yzz_yyz = pbuffer.data(idx_op_geom_010_ff + 187);

    auto tr_0_0_y_yzz_yzz = pbuffer.data(idx_op_geom_010_ff + 188);

    auto tr_0_0_y_yzz_zzz = pbuffer.data(idx_op_geom_010_ff + 189);

    #pragma omp simd aligned(tr_0_0_y_yzz_xxx, tr_0_0_y_yzz_xxy, tr_0_0_y_yzz_xxz, tr_0_0_y_yzz_xyy, tr_0_0_y_yzz_xyz, tr_0_0_y_yzz_xzz, tr_0_0_y_yzz_yyy, tr_0_0_y_yzz_yyz, tr_0_0_y_yzz_yzz, tr_0_0_y_yzz_zzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, tr_yzz_xx, tr_yzz_xxxy, tr_yzz_xxyy, tr_yzz_xxyz, tr_yzz_xy, tr_yzz_xyyy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_yy, tr_yzz_yyyy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzz_xxx[i] = 2.0 * tr_yyzz_xxx[i] * tbe_0 + 2.0 * tr_yzz_xxxy[i] * tke_0 - tr_zz_xxx[i];

        tr_0_0_y_yzz_xxy[i] = 2.0 * tr_yyzz_xxy[i] * tbe_0 + 2.0 * tr_yzz_xxyy[i] * tke_0 - tr_zz_xxy[i] - tr_yzz_xx[i];

        tr_0_0_y_yzz_xxz[i] = 2.0 * tr_yyzz_xxz[i] * tbe_0 + 2.0 * tr_yzz_xxyz[i] * tke_0 - tr_zz_xxz[i];

        tr_0_0_y_yzz_xyy[i] = 2.0 * tr_yyzz_xyy[i] * tbe_0 + 2.0 * tr_yzz_xyyy[i] * tke_0 - tr_zz_xyy[i] - 2.0 * tr_yzz_xy[i];

        tr_0_0_y_yzz_xyz[i] = 2.0 * tr_yyzz_xyz[i] * tbe_0 + 2.0 * tr_yzz_xyyz[i] * tke_0 - tr_zz_xyz[i] - tr_yzz_xz[i];

        tr_0_0_y_yzz_xzz[i] = 2.0 * tr_yyzz_xzz[i] * tbe_0 + 2.0 * tr_yzz_xyzz[i] * tke_0 - tr_zz_xzz[i];

        tr_0_0_y_yzz_yyy[i] = 2.0 * tr_yyzz_yyy[i] * tbe_0 + 2.0 * tr_yzz_yyyy[i] * tke_0 - tr_zz_yyy[i] - 3.0 * tr_yzz_yy[i];

        tr_0_0_y_yzz_yyz[i] = 2.0 * tr_yyzz_yyz[i] * tbe_0 + 2.0 * tr_yzz_yyyz[i] * tke_0 - tr_zz_yyz[i] - 2.0 * tr_yzz_yz[i];

        tr_0_0_y_yzz_yzz[i] = 2.0 * tr_yyzz_yzz[i] * tbe_0 + 2.0 * tr_yzz_yyzz[i] * tke_0 - tr_zz_yzz[i] - tr_yzz_zz[i];

        tr_0_0_y_yzz_zzz[i] = 2.0 * tr_yyzz_zzz[i] * tbe_0 + 2.0 * tr_yzz_yzzz[i] * tke_0 - tr_zz_zzz[i];
    }

    // Set up 190-200 components of targeted buffer : FF

    auto tr_0_0_y_zzz_xxx = pbuffer.data(idx_op_geom_010_ff + 190);

    auto tr_0_0_y_zzz_xxy = pbuffer.data(idx_op_geom_010_ff + 191);

    auto tr_0_0_y_zzz_xxz = pbuffer.data(idx_op_geom_010_ff + 192);

    auto tr_0_0_y_zzz_xyy = pbuffer.data(idx_op_geom_010_ff + 193);

    auto tr_0_0_y_zzz_xyz = pbuffer.data(idx_op_geom_010_ff + 194);

    auto tr_0_0_y_zzz_xzz = pbuffer.data(idx_op_geom_010_ff + 195);

    auto tr_0_0_y_zzz_yyy = pbuffer.data(idx_op_geom_010_ff + 196);

    auto tr_0_0_y_zzz_yyz = pbuffer.data(idx_op_geom_010_ff + 197);

    auto tr_0_0_y_zzz_yzz = pbuffer.data(idx_op_geom_010_ff + 198);

    auto tr_0_0_y_zzz_zzz = pbuffer.data(idx_op_geom_010_ff + 199);

    #pragma omp simd aligned(tr_0_0_y_zzz_xxx, tr_0_0_y_zzz_xxy, tr_0_0_y_zzz_xxz, tr_0_0_y_zzz_xyy, tr_0_0_y_zzz_xyz, tr_0_0_y_zzz_xzz, tr_0_0_y_zzz_yyy, tr_0_0_y_zzz_yyz, tr_0_0_y_zzz_yzz, tr_0_0_y_zzz_zzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, tr_zzz_xx, tr_zzz_xxxy, tr_zzz_xxyy, tr_zzz_xxyz, tr_zzz_xy, tr_zzz_xyyy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_yy, tr_zzz_yyyy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzz_xxx[i] = 2.0 * tr_yzzz_xxx[i] * tbe_0 + 2.0 * tr_zzz_xxxy[i] * tke_0;

        tr_0_0_y_zzz_xxy[i] = 2.0 * tr_yzzz_xxy[i] * tbe_0 + 2.0 * tr_zzz_xxyy[i] * tke_0 - tr_zzz_xx[i];

        tr_0_0_y_zzz_xxz[i] = 2.0 * tr_yzzz_xxz[i] * tbe_0 + 2.0 * tr_zzz_xxyz[i] * tke_0;

        tr_0_0_y_zzz_xyy[i] = 2.0 * tr_yzzz_xyy[i] * tbe_0 + 2.0 * tr_zzz_xyyy[i] * tke_0 - 2.0 * tr_zzz_xy[i];

        tr_0_0_y_zzz_xyz[i] = 2.0 * tr_yzzz_xyz[i] * tbe_0 + 2.0 * tr_zzz_xyyz[i] * tke_0 - tr_zzz_xz[i];

        tr_0_0_y_zzz_xzz[i] = 2.0 * tr_yzzz_xzz[i] * tbe_0 + 2.0 * tr_zzz_xyzz[i] * tke_0;

        tr_0_0_y_zzz_yyy[i] = 2.0 * tr_yzzz_yyy[i] * tbe_0 + 2.0 * tr_zzz_yyyy[i] * tke_0 - 3.0 * tr_zzz_yy[i];

        tr_0_0_y_zzz_yyz[i] = 2.0 * tr_yzzz_yyz[i] * tbe_0 + 2.0 * tr_zzz_yyyz[i] * tke_0 - 2.0 * tr_zzz_yz[i];

        tr_0_0_y_zzz_yzz[i] = 2.0 * tr_yzzz_yzz[i] * tbe_0 + 2.0 * tr_zzz_yyzz[i] * tke_0 - tr_zzz_zz[i];

        tr_0_0_y_zzz_zzz[i] = 2.0 * tr_yzzz_zzz[i] * tbe_0 + 2.0 * tr_zzz_yzzz[i] * tke_0;
    }

    // Set up 200-210 components of targeted buffer : FF

    auto tr_0_0_z_xxx_xxx = pbuffer.data(idx_op_geom_010_ff + 200);

    auto tr_0_0_z_xxx_xxy = pbuffer.data(idx_op_geom_010_ff + 201);

    auto tr_0_0_z_xxx_xxz = pbuffer.data(idx_op_geom_010_ff + 202);

    auto tr_0_0_z_xxx_xyy = pbuffer.data(idx_op_geom_010_ff + 203);

    auto tr_0_0_z_xxx_xyz = pbuffer.data(idx_op_geom_010_ff + 204);

    auto tr_0_0_z_xxx_xzz = pbuffer.data(idx_op_geom_010_ff + 205);

    auto tr_0_0_z_xxx_yyy = pbuffer.data(idx_op_geom_010_ff + 206);

    auto tr_0_0_z_xxx_yyz = pbuffer.data(idx_op_geom_010_ff + 207);

    auto tr_0_0_z_xxx_yzz = pbuffer.data(idx_op_geom_010_ff + 208);

    auto tr_0_0_z_xxx_zzz = pbuffer.data(idx_op_geom_010_ff + 209);

    #pragma omp simd aligned(tr_0_0_z_xxx_xxx, tr_0_0_z_xxx_xxy, tr_0_0_z_xxx_xxz, tr_0_0_z_xxx_xyy, tr_0_0_z_xxx_xyz, tr_0_0_z_xxx_xzz, tr_0_0_z_xxx_yyy, tr_0_0_z_xxx_yyz, tr_0_0_z_xxx_yzz, tr_0_0_z_xxx_zzz, tr_xxx_xx, tr_xxx_xxxz, tr_xxx_xxyz, tr_xxx_xxzz, tr_xxx_xy, tr_xxx_xyyz, tr_xxx_xyzz, tr_xxx_xz, tr_xxx_xzzz, tr_xxx_yy, tr_xxx_yyyz, tr_xxx_yyzz, tr_xxx_yz, tr_xxx_yzzz, tr_xxx_zz, tr_xxx_zzzz, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxx_xxx[i] = 2.0 * tr_xxxz_xxx[i] * tbe_0 + 2.0 * tr_xxx_xxxz[i] * tke_0;

        tr_0_0_z_xxx_xxy[i] = 2.0 * tr_xxxz_xxy[i] * tbe_0 + 2.0 * tr_xxx_xxyz[i] * tke_0;

        tr_0_0_z_xxx_xxz[i] = 2.0 * tr_xxxz_xxz[i] * tbe_0 + 2.0 * tr_xxx_xxzz[i] * tke_0 - tr_xxx_xx[i];

        tr_0_0_z_xxx_xyy[i] = 2.0 * tr_xxxz_xyy[i] * tbe_0 + 2.0 * tr_xxx_xyyz[i] * tke_0;

        tr_0_0_z_xxx_xyz[i] = 2.0 * tr_xxxz_xyz[i] * tbe_0 + 2.0 * tr_xxx_xyzz[i] * tke_0 - tr_xxx_xy[i];

        tr_0_0_z_xxx_xzz[i] = 2.0 * tr_xxxz_xzz[i] * tbe_0 + 2.0 * tr_xxx_xzzz[i] * tke_0 - 2.0 * tr_xxx_xz[i];

        tr_0_0_z_xxx_yyy[i] = 2.0 * tr_xxxz_yyy[i] * tbe_0 + 2.0 * tr_xxx_yyyz[i] * tke_0;

        tr_0_0_z_xxx_yyz[i] = 2.0 * tr_xxxz_yyz[i] * tbe_0 + 2.0 * tr_xxx_yyzz[i] * tke_0 - tr_xxx_yy[i];

        tr_0_0_z_xxx_yzz[i] = 2.0 * tr_xxxz_yzz[i] * tbe_0 + 2.0 * tr_xxx_yzzz[i] * tke_0 - 2.0 * tr_xxx_yz[i];

        tr_0_0_z_xxx_zzz[i] = 2.0 * tr_xxxz_zzz[i] * tbe_0 + 2.0 * tr_xxx_zzzz[i] * tke_0 - 3.0 * tr_xxx_zz[i];
    }

    // Set up 210-220 components of targeted buffer : FF

    auto tr_0_0_z_xxy_xxx = pbuffer.data(idx_op_geom_010_ff + 210);

    auto tr_0_0_z_xxy_xxy = pbuffer.data(idx_op_geom_010_ff + 211);

    auto tr_0_0_z_xxy_xxz = pbuffer.data(idx_op_geom_010_ff + 212);

    auto tr_0_0_z_xxy_xyy = pbuffer.data(idx_op_geom_010_ff + 213);

    auto tr_0_0_z_xxy_xyz = pbuffer.data(idx_op_geom_010_ff + 214);

    auto tr_0_0_z_xxy_xzz = pbuffer.data(idx_op_geom_010_ff + 215);

    auto tr_0_0_z_xxy_yyy = pbuffer.data(idx_op_geom_010_ff + 216);

    auto tr_0_0_z_xxy_yyz = pbuffer.data(idx_op_geom_010_ff + 217);

    auto tr_0_0_z_xxy_yzz = pbuffer.data(idx_op_geom_010_ff + 218);

    auto tr_0_0_z_xxy_zzz = pbuffer.data(idx_op_geom_010_ff + 219);

    #pragma omp simd aligned(tr_0_0_z_xxy_xxx, tr_0_0_z_xxy_xxy, tr_0_0_z_xxy_xxz, tr_0_0_z_xxy_xyy, tr_0_0_z_xxy_xyz, tr_0_0_z_xxy_xzz, tr_0_0_z_xxy_yyy, tr_0_0_z_xxy_yyz, tr_0_0_z_xxy_yzz, tr_0_0_z_xxy_zzz, tr_xxy_xx, tr_xxy_xxxz, tr_xxy_xxyz, tr_xxy_xxzz, tr_xxy_xy, tr_xxy_xyyz, tr_xxy_xyzz, tr_xxy_xz, tr_xxy_xzzz, tr_xxy_yy, tr_xxy_yyyz, tr_xxy_yyzz, tr_xxy_yz, tr_xxy_yzzz, tr_xxy_zz, tr_xxy_zzzz, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxy_xxx[i] = 2.0 * tr_xxyz_xxx[i] * tbe_0 + 2.0 * tr_xxy_xxxz[i] * tke_0;

        tr_0_0_z_xxy_xxy[i] = 2.0 * tr_xxyz_xxy[i] * tbe_0 + 2.0 * tr_xxy_xxyz[i] * tke_0;

        tr_0_0_z_xxy_xxz[i] = 2.0 * tr_xxyz_xxz[i] * tbe_0 + 2.0 * tr_xxy_xxzz[i] * tke_0 - tr_xxy_xx[i];

        tr_0_0_z_xxy_xyy[i] = 2.0 * tr_xxyz_xyy[i] * tbe_0 + 2.0 * tr_xxy_xyyz[i] * tke_0;

        tr_0_0_z_xxy_xyz[i] = 2.0 * tr_xxyz_xyz[i] * tbe_0 + 2.0 * tr_xxy_xyzz[i] * tke_0 - tr_xxy_xy[i];

        tr_0_0_z_xxy_xzz[i] = 2.0 * tr_xxyz_xzz[i] * tbe_0 + 2.0 * tr_xxy_xzzz[i] * tke_0 - 2.0 * tr_xxy_xz[i];

        tr_0_0_z_xxy_yyy[i] = 2.0 * tr_xxyz_yyy[i] * tbe_0 + 2.0 * tr_xxy_yyyz[i] * tke_0;

        tr_0_0_z_xxy_yyz[i] = 2.0 * tr_xxyz_yyz[i] * tbe_0 + 2.0 * tr_xxy_yyzz[i] * tke_0 - tr_xxy_yy[i];

        tr_0_0_z_xxy_yzz[i] = 2.0 * tr_xxyz_yzz[i] * tbe_0 + 2.0 * tr_xxy_yzzz[i] * tke_0 - 2.0 * tr_xxy_yz[i];

        tr_0_0_z_xxy_zzz[i] = 2.0 * tr_xxyz_zzz[i] * tbe_0 + 2.0 * tr_xxy_zzzz[i] * tke_0 - 3.0 * tr_xxy_zz[i];
    }

    // Set up 220-230 components of targeted buffer : FF

    auto tr_0_0_z_xxz_xxx = pbuffer.data(idx_op_geom_010_ff + 220);

    auto tr_0_0_z_xxz_xxy = pbuffer.data(idx_op_geom_010_ff + 221);

    auto tr_0_0_z_xxz_xxz = pbuffer.data(idx_op_geom_010_ff + 222);

    auto tr_0_0_z_xxz_xyy = pbuffer.data(idx_op_geom_010_ff + 223);

    auto tr_0_0_z_xxz_xyz = pbuffer.data(idx_op_geom_010_ff + 224);

    auto tr_0_0_z_xxz_xzz = pbuffer.data(idx_op_geom_010_ff + 225);

    auto tr_0_0_z_xxz_yyy = pbuffer.data(idx_op_geom_010_ff + 226);

    auto tr_0_0_z_xxz_yyz = pbuffer.data(idx_op_geom_010_ff + 227);

    auto tr_0_0_z_xxz_yzz = pbuffer.data(idx_op_geom_010_ff + 228);

    auto tr_0_0_z_xxz_zzz = pbuffer.data(idx_op_geom_010_ff + 229);

    #pragma omp simd aligned(tr_0_0_z_xxz_xxx, tr_0_0_z_xxz_xxy, tr_0_0_z_xxz_xxz, tr_0_0_z_xxz_xyy, tr_0_0_z_xxz_xyz, tr_0_0_z_xxz_xzz, tr_0_0_z_xxz_yyy, tr_0_0_z_xxz_yyz, tr_0_0_z_xxz_yzz, tr_0_0_z_xxz_zzz, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_zzz, tr_xxz_xx, tr_xxz_xxxz, tr_xxz_xxyz, tr_xxz_xxzz, tr_xxz_xy, tr_xxz_xyyz, tr_xxz_xyzz, tr_xxz_xz, tr_xxz_xzzz, tr_xxz_yy, tr_xxz_yyyz, tr_xxz_yyzz, tr_xxz_yz, tr_xxz_yzzz, tr_xxz_zz, tr_xxz_zzzz, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxz_xxx[i] = 2.0 * tr_xxzz_xxx[i] * tbe_0 + 2.0 * tr_xxz_xxxz[i] * tke_0 - tr_xx_xxx[i];

        tr_0_0_z_xxz_xxy[i] = 2.0 * tr_xxzz_xxy[i] * tbe_0 + 2.0 * tr_xxz_xxyz[i] * tke_0 - tr_xx_xxy[i];

        tr_0_0_z_xxz_xxz[i] = 2.0 * tr_xxzz_xxz[i] * tbe_0 + 2.0 * tr_xxz_xxzz[i] * tke_0 - tr_xx_xxz[i] - tr_xxz_xx[i];

        tr_0_0_z_xxz_xyy[i] = 2.0 * tr_xxzz_xyy[i] * tbe_0 + 2.0 * tr_xxz_xyyz[i] * tke_0 - tr_xx_xyy[i];

        tr_0_0_z_xxz_xyz[i] = 2.0 * tr_xxzz_xyz[i] * tbe_0 + 2.0 * tr_xxz_xyzz[i] * tke_0 - tr_xx_xyz[i] - tr_xxz_xy[i];

        tr_0_0_z_xxz_xzz[i] = 2.0 * tr_xxzz_xzz[i] * tbe_0 + 2.0 * tr_xxz_xzzz[i] * tke_0 - tr_xx_xzz[i] - 2.0 * tr_xxz_xz[i];

        tr_0_0_z_xxz_yyy[i] = 2.0 * tr_xxzz_yyy[i] * tbe_0 + 2.0 * tr_xxz_yyyz[i] * tke_0 - tr_xx_yyy[i];

        tr_0_0_z_xxz_yyz[i] = 2.0 * tr_xxzz_yyz[i] * tbe_0 + 2.0 * tr_xxz_yyzz[i] * tke_0 - tr_xx_yyz[i] - tr_xxz_yy[i];

        tr_0_0_z_xxz_yzz[i] = 2.0 * tr_xxzz_yzz[i] * tbe_0 + 2.0 * tr_xxz_yzzz[i] * tke_0 - tr_xx_yzz[i] - 2.0 * tr_xxz_yz[i];

        tr_0_0_z_xxz_zzz[i] = 2.0 * tr_xxzz_zzz[i] * tbe_0 + 2.0 * tr_xxz_zzzz[i] * tke_0 - tr_xx_zzz[i] - 3.0 * tr_xxz_zz[i];
    }

    // Set up 230-240 components of targeted buffer : FF

    auto tr_0_0_z_xyy_xxx = pbuffer.data(idx_op_geom_010_ff + 230);

    auto tr_0_0_z_xyy_xxy = pbuffer.data(idx_op_geom_010_ff + 231);

    auto tr_0_0_z_xyy_xxz = pbuffer.data(idx_op_geom_010_ff + 232);

    auto tr_0_0_z_xyy_xyy = pbuffer.data(idx_op_geom_010_ff + 233);

    auto tr_0_0_z_xyy_xyz = pbuffer.data(idx_op_geom_010_ff + 234);

    auto tr_0_0_z_xyy_xzz = pbuffer.data(idx_op_geom_010_ff + 235);

    auto tr_0_0_z_xyy_yyy = pbuffer.data(idx_op_geom_010_ff + 236);

    auto tr_0_0_z_xyy_yyz = pbuffer.data(idx_op_geom_010_ff + 237);

    auto tr_0_0_z_xyy_yzz = pbuffer.data(idx_op_geom_010_ff + 238);

    auto tr_0_0_z_xyy_zzz = pbuffer.data(idx_op_geom_010_ff + 239);

    #pragma omp simd aligned(tr_0_0_z_xyy_xxx, tr_0_0_z_xyy_xxy, tr_0_0_z_xyy_xxz, tr_0_0_z_xyy_xyy, tr_0_0_z_xyy_xyz, tr_0_0_z_xyy_xzz, tr_0_0_z_xyy_yyy, tr_0_0_z_xyy_yyz, tr_0_0_z_xyy_yzz, tr_0_0_z_xyy_zzz, tr_xyy_xx, tr_xyy_xxxz, tr_xyy_xxyz, tr_xyy_xxzz, tr_xyy_xy, tr_xyy_xyyz, tr_xyy_xyzz, tr_xyy_xz, tr_xyy_xzzz, tr_xyy_yy, tr_xyy_yyyz, tr_xyy_yyzz, tr_xyy_yz, tr_xyy_yzzz, tr_xyy_zz, tr_xyy_zzzz, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyy_xxx[i] = 2.0 * tr_xyyz_xxx[i] * tbe_0 + 2.0 * tr_xyy_xxxz[i] * tke_0;

        tr_0_0_z_xyy_xxy[i] = 2.0 * tr_xyyz_xxy[i] * tbe_0 + 2.0 * tr_xyy_xxyz[i] * tke_0;

        tr_0_0_z_xyy_xxz[i] = 2.0 * tr_xyyz_xxz[i] * tbe_0 + 2.0 * tr_xyy_xxzz[i] * tke_0 - tr_xyy_xx[i];

        tr_0_0_z_xyy_xyy[i] = 2.0 * tr_xyyz_xyy[i] * tbe_0 + 2.0 * tr_xyy_xyyz[i] * tke_0;

        tr_0_0_z_xyy_xyz[i] = 2.0 * tr_xyyz_xyz[i] * tbe_0 + 2.0 * tr_xyy_xyzz[i] * tke_0 - tr_xyy_xy[i];

        tr_0_0_z_xyy_xzz[i] = 2.0 * tr_xyyz_xzz[i] * tbe_0 + 2.0 * tr_xyy_xzzz[i] * tke_0 - 2.0 * tr_xyy_xz[i];

        tr_0_0_z_xyy_yyy[i] = 2.0 * tr_xyyz_yyy[i] * tbe_0 + 2.0 * tr_xyy_yyyz[i] * tke_0;

        tr_0_0_z_xyy_yyz[i] = 2.0 * tr_xyyz_yyz[i] * tbe_0 + 2.0 * tr_xyy_yyzz[i] * tke_0 - tr_xyy_yy[i];

        tr_0_0_z_xyy_yzz[i] = 2.0 * tr_xyyz_yzz[i] * tbe_0 + 2.0 * tr_xyy_yzzz[i] * tke_0 - 2.0 * tr_xyy_yz[i];

        tr_0_0_z_xyy_zzz[i] = 2.0 * tr_xyyz_zzz[i] * tbe_0 + 2.0 * tr_xyy_zzzz[i] * tke_0 - 3.0 * tr_xyy_zz[i];
    }

    // Set up 240-250 components of targeted buffer : FF

    auto tr_0_0_z_xyz_xxx = pbuffer.data(idx_op_geom_010_ff + 240);

    auto tr_0_0_z_xyz_xxy = pbuffer.data(idx_op_geom_010_ff + 241);

    auto tr_0_0_z_xyz_xxz = pbuffer.data(idx_op_geom_010_ff + 242);

    auto tr_0_0_z_xyz_xyy = pbuffer.data(idx_op_geom_010_ff + 243);

    auto tr_0_0_z_xyz_xyz = pbuffer.data(idx_op_geom_010_ff + 244);

    auto tr_0_0_z_xyz_xzz = pbuffer.data(idx_op_geom_010_ff + 245);

    auto tr_0_0_z_xyz_yyy = pbuffer.data(idx_op_geom_010_ff + 246);

    auto tr_0_0_z_xyz_yyz = pbuffer.data(idx_op_geom_010_ff + 247);

    auto tr_0_0_z_xyz_yzz = pbuffer.data(idx_op_geom_010_ff + 248);

    auto tr_0_0_z_xyz_zzz = pbuffer.data(idx_op_geom_010_ff + 249);

    #pragma omp simd aligned(tr_0_0_z_xyz_xxx, tr_0_0_z_xyz_xxy, tr_0_0_z_xyz_xxz, tr_0_0_z_xyz_xyy, tr_0_0_z_xyz_xyz, tr_0_0_z_xyz_xzz, tr_0_0_z_xyz_yyy, tr_0_0_z_xyz_yyz, tr_0_0_z_xyz_yzz, tr_0_0_z_xyz_zzz, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_zzz, tr_xyz_xx, tr_xyz_xxxz, tr_xyz_xxyz, tr_xyz_xxzz, tr_xyz_xy, tr_xyz_xyyz, tr_xyz_xyzz, tr_xyz_xz, tr_xyz_xzzz, tr_xyz_yy, tr_xyz_yyyz, tr_xyz_yyzz, tr_xyz_yz, tr_xyz_yzzz, tr_xyz_zz, tr_xyz_zzzz, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyz_xxx[i] = 2.0 * tr_xyzz_xxx[i] * tbe_0 + 2.0 * tr_xyz_xxxz[i] * tke_0 - tr_xy_xxx[i];

        tr_0_0_z_xyz_xxy[i] = 2.0 * tr_xyzz_xxy[i] * tbe_0 + 2.0 * tr_xyz_xxyz[i] * tke_0 - tr_xy_xxy[i];

        tr_0_0_z_xyz_xxz[i] = 2.0 * tr_xyzz_xxz[i] * tbe_0 + 2.0 * tr_xyz_xxzz[i] * tke_0 - tr_xy_xxz[i] - tr_xyz_xx[i];

        tr_0_0_z_xyz_xyy[i] = 2.0 * tr_xyzz_xyy[i] * tbe_0 + 2.0 * tr_xyz_xyyz[i] * tke_0 - tr_xy_xyy[i];

        tr_0_0_z_xyz_xyz[i] = 2.0 * tr_xyzz_xyz[i] * tbe_0 + 2.0 * tr_xyz_xyzz[i] * tke_0 - tr_xy_xyz[i] - tr_xyz_xy[i];

        tr_0_0_z_xyz_xzz[i] = 2.0 * tr_xyzz_xzz[i] * tbe_0 + 2.0 * tr_xyz_xzzz[i] * tke_0 - tr_xy_xzz[i] - 2.0 * tr_xyz_xz[i];

        tr_0_0_z_xyz_yyy[i] = 2.0 * tr_xyzz_yyy[i] * tbe_0 + 2.0 * tr_xyz_yyyz[i] * tke_0 - tr_xy_yyy[i];

        tr_0_0_z_xyz_yyz[i] = 2.0 * tr_xyzz_yyz[i] * tbe_0 + 2.0 * tr_xyz_yyzz[i] * tke_0 - tr_xy_yyz[i] - tr_xyz_yy[i];

        tr_0_0_z_xyz_yzz[i] = 2.0 * tr_xyzz_yzz[i] * tbe_0 + 2.0 * tr_xyz_yzzz[i] * tke_0 - tr_xy_yzz[i] - 2.0 * tr_xyz_yz[i];

        tr_0_0_z_xyz_zzz[i] = 2.0 * tr_xyzz_zzz[i] * tbe_0 + 2.0 * tr_xyz_zzzz[i] * tke_0 - tr_xy_zzz[i] - 3.0 * tr_xyz_zz[i];
    }

    // Set up 250-260 components of targeted buffer : FF

    auto tr_0_0_z_xzz_xxx = pbuffer.data(idx_op_geom_010_ff + 250);

    auto tr_0_0_z_xzz_xxy = pbuffer.data(idx_op_geom_010_ff + 251);

    auto tr_0_0_z_xzz_xxz = pbuffer.data(idx_op_geom_010_ff + 252);

    auto tr_0_0_z_xzz_xyy = pbuffer.data(idx_op_geom_010_ff + 253);

    auto tr_0_0_z_xzz_xyz = pbuffer.data(idx_op_geom_010_ff + 254);

    auto tr_0_0_z_xzz_xzz = pbuffer.data(idx_op_geom_010_ff + 255);

    auto tr_0_0_z_xzz_yyy = pbuffer.data(idx_op_geom_010_ff + 256);

    auto tr_0_0_z_xzz_yyz = pbuffer.data(idx_op_geom_010_ff + 257);

    auto tr_0_0_z_xzz_yzz = pbuffer.data(idx_op_geom_010_ff + 258);

    auto tr_0_0_z_xzz_zzz = pbuffer.data(idx_op_geom_010_ff + 259);

    #pragma omp simd aligned(tr_0_0_z_xzz_xxx, tr_0_0_z_xzz_xxy, tr_0_0_z_xzz_xxz, tr_0_0_z_xzz_xyy, tr_0_0_z_xzz_xyz, tr_0_0_z_xzz_xzz, tr_0_0_z_xzz_yyy, tr_0_0_z_xzz_yyz, tr_0_0_z_xzz_yzz, tr_0_0_z_xzz_zzz, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_zzz, tr_xzz_xx, tr_xzz_xxxz, tr_xzz_xxyz, tr_xzz_xxzz, tr_xzz_xy, tr_xzz_xyyz, tr_xzz_xyzz, tr_xzz_xz, tr_xzz_xzzz, tr_xzz_yy, tr_xzz_yyyz, tr_xzz_yyzz, tr_xzz_yz, tr_xzz_yzzz, tr_xzz_zz, tr_xzz_zzzz, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzz_xxx[i] = 2.0 * tr_xzzz_xxx[i] * tbe_0 + 2.0 * tr_xzz_xxxz[i] * tke_0 - 2.0 * tr_xz_xxx[i];

        tr_0_0_z_xzz_xxy[i] = 2.0 * tr_xzzz_xxy[i] * tbe_0 + 2.0 * tr_xzz_xxyz[i] * tke_0 - 2.0 * tr_xz_xxy[i];

        tr_0_0_z_xzz_xxz[i] = 2.0 * tr_xzzz_xxz[i] * tbe_0 + 2.0 * tr_xzz_xxzz[i] * tke_0 - 2.0 * tr_xz_xxz[i] - tr_xzz_xx[i];

        tr_0_0_z_xzz_xyy[i] = 2.0 * tr_xzzz_xyy[i] * tbe_0 + 2.0 * tr_xzz_xyyz[i] * tke_0 - 2.0 * tr_xz_xyy[i];

        tr_0_0_z_xzz_xyz[i] = 2.0 * tr_xzzz_xyz[i] * tbe_0 + 2.0 * tr_xzz_xyzz[i] * tke_0 - 2.0 * tr_xz_xyz[i] - tr_xzz_xy[i];

        tr_0_0_z_xzz_xzz[i] = 2.0 * tr_xzzz_xzz[i] * tbe_0 + 2.0 * tr_xzz_xzzz[i] * tke_0 - 2.0 * tr_xz_xzz[i] - 2.0 * tr_xzz_xz[i];

        tr_0_0_z_xzz_yyy[i] = 2.0 * tr_xzzz_yyy[i] * tbe_0 + 2.0 * tr_xzz_yyyz[i] * tke_0 - 2.0 * tr_xz_yyy[i];

        tr_0_0_z_xzz_yyz[i] = 2.0 * tr_xzzz_yyz[i] * tbe_0 + 2.0 * tr_xzz_yyzz[i] * tke_0 - 2.0 * tr_xz_yyz[i] - tr_xzz_yy[i];

        tr_0_0_z_xzz_yzz[i] = 2.0 * tr_xzzz_yzz[i] * tbe_0 + 2.0 * tr_xzz_yzzz[i] * tke_0 - 2.0 * tr_xz_yzz[i] - 2.0 * tr_xzz_yz[i];

        tr_0_0_z_xzz_zzz[i] = 2.0 * tr_xzzz_zzz[i] * tbe_0 + 2.0 * tr_xzz_zzzz[i] * tke_0 - 2.0 * tr_xz_zzz[i] - 3.0 * tr_xzz_zz[i];
    }

    // Set up 260-270 components of targeted buffer : FF

    auto tr_0_0_z_yyy_xxx = pbuffer.data(idx_op_geom_010_ff + 260);

    auto tr_0_0_z_yyy_xxy = pbuffer.data(idx_op_geom_010_ff + 261);

    auto tr_0_0_z_yyy_xxz = pbuffer.data(idx_op_geom_010_ff + 262);

    auto tr_0_0_z_yyy_xyy = pbuffer.data(idx_op_geom_010_ff + 263);

    auto tr_0_0_z_yyy_xyz = pbuffer.data(idx_op_geom_010_ff + 264);

    auto tr_0_0_z_yyy_xzz = pbuffer.data(idx_op_geom_010_ff + 265);

    auto tr_0_0_z_yyy_yyy = pbuffer.data(idx_op_geom_010_ff + 266);

    auto tr_0_0_z_yyy_yyz = pbuffer.data(idx_op_geom_010_ff + 267);

    auto tr_0_0_z_yyy_yzz = pbuffer.data(idx_op_geom_010_ff + 268);

    auto tr_0_0_z_yyy_zzz = pbuffer.data(idx_op_geom_010_ff + 269);

    #pragma omp simd aligned(tr_0_0_z_yyy_xxx, tr_0_0_z_yyy_xxy, tr_0_0_z_yyy_xxz, tr_0_0_z_yyy_xyy, tr_0_0_z_yyy_xyz, tr_0_0_z_yyy_xzz, tr_0_0_z_yyy_yyy, tr_0_0_z_yyy_yyz, tr_0_0_z_yyy_yzz, tr_0_0_z_yyy_zzz, tr_yyy_xx, tr_yyy_xxxz, tr_yyy_xxyz, tr_yyy_xxzz, tr_yyy_xy, tr_yyy_xyyz, tr_yyy_xyzz, tr_yyy_xz, tr_yyy_xzzz, tr_yyy_yy, tr_yyy_yyyz, tr_yyy_yyzz, tr_yyy_yz, tr_yyy_yzzz, tr_yyy_zz, tr_yyy_zzzz, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyy_xxx[i] = 2.0 * tr_yyyz_xxx[i] * tbe_0 + 2.0 * tr_yyy_xxxz[i] * tke_0;

        tr_0_0_z_yyy_xxy[i] = 2.0 * tr_yyyz_xxy[i] * tbe_0 + 2.0 * tr_yyy_xxyz[i] * tke_0;

        tr_0_0_z_yyy_xxz[i] = 2.0 * tr_yyyz_xxz[i] * tbe_0 + 2.0 * tr_yyy_xxzz[i] * tke_0 - tr_yyy_xx[i];

        tr_0_0_z_yyy_xyy[i] = 2.0 * tr_yyyz_xyy[i] * tbe_0 + 2.0 * tr_yyy_xyyz[i] * tke_0;

        tr_0_0_z_yyy_xyz[i] = 2.0 * tr_yyyz_xyz[i] * tbe_0 + 2.0 * tr_yyy_xyzz[i] * tke_0 - tr_yyy_xy[i];

        tr_0_0_z_yyy_xzz[i] = 2.0 * tr_yyyz_xzz[i] * tbe_0 + 2.0 * tr_yyy_xzzz[i] * tke_0 - 2.0 * tr_yyy_xz[i];

        tr_0_0_z_yyy_yyy[i] = 2.0 * tr_yyyz_yyy[i] * tbe_0 + 2.0 * tr_yyy_yyyz[i] * tke_0;

        tr_0_0_z_yyy_yyz[i] = 2.0 * tr_yyyz_yyz[i] * tbe_0 + 2.0 * tr_yyy_yyzz[i] * tke_0 - tr_yyy_yy[i];

        tr_0_0_z_yyy_yzz[i] = 2.0 * tr_yyyz_yzz[i] * tbe_0 + 2.0 * tr_yyy_yzzz[i] * tke_0 - 2.0 * tr_yyy_yz[i];

        tr_0_0_z_yyy_zzz[i] = 2.0 * tr_yyyz_zzz[i] * tbe_0 + 2.0 * tr_yyy_zzzz[i] * tke_0 - 3.0 * tr_yyy_zz[i];
    }

    // Set up 270-280 components of targeted buffer : FF

    auto tr_0_0_z_yyz_xxx = pbuffer.data(idx_op_geom_010_ff + 270);

    auto tr_0_0_z_yyz_xxy = pbuffer.data(idx_op_geom_010_ff + 271);

    auto tr_0_0_z_yyz_xxz = pbuffer.data(idx_op_geom_010_ff + 272);

    auto tr_0_0_z_yyz_xyy = pbuffer.data(idx_op_geom_010_ff + 273);

    auto tr_0_0_z_yyz_xyz = pbuffer.data(idx_op_geom_010_ff + 274);

    auto tr_0_0_z_yyz_xzz = pbuffer.data(idx_op_geom_010_ff + 275);

    auto tr_0_0_z_yyz_yyy = pbuffer.data(idx_op_geom_010_ff + 276);

    auto tr_0_0_z_yyz_yyz = pbuffer.data(idx_op_geom_010_ff + 277);

    auto tr_0_0_z_yyz_yzz = pbuffer.data(idx_op_geom_010_ff + 278);

    auto tr_0_0_z_yyz_zzz = pbuffer.data(idx_op_geom_010_ff + 279);

    #pragma omp simd aligned(tr_0_0_z_yyz_xxx, tr_0_0_z_yyz_xxy, tr_0_0_z_yyz_xxz, tr_0_0_z_yyz_xyy, tr_0_0_z_yyz_xyz, tr_0_0_z_yyz_xzz, tr_0_0_z_yyz_yyy, tr_0_0_z_yyz_yyz, tr_0_0_z_yyz_yzz, tr_0_0_z_yyz_zzz, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_zzz, tr_yyz_xx, tr_yyz_xxxz, tr_yyz_xxyz, tr_yyz_xxzz, tr_yyz_xy, tr_yyz_xyyz, tr_yyz_xyzz, tr_yyz_xz, tr_yyz_xzzz, tr_yyz_yy, tr_yyz_yyyz, tr_yyz_yyzz, tr_yyz_yz, tr_yyz_yzzz, tr_yyz_zz, tr_yyz_zzzz, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyz_xxx[i] = 2.0 * tr_yyzz_xxx[i] * tbe_0 + 2.0 * tr_yyz_xxxz[i] * tke_0 - tr_yy_xxx[i];

        tr_0_0_z_yyz_xxy[i] = 2.0 * tr_yyzz_xxy[i] * tbe_0 + 2.0 * tr_yyz_xxyz[i] * tke_0 - tr_yy_xxy[i];

        tr_0_0_z_yyz_xxz[i] = 2.0 * tr_yyzz_xxz[i] * tbe_0 + 2.0 * tr_yyz_xxzz[i] * tke_0 - tr_yy_xxz[i] - tr_yyz_xx[i];

        tr_0_0_z_yyz_xyy[i] = 2.0 * tr_yyzz_xyy[i] * tbe_0 + 2.0 * tr_yyz_xyyz[i] * tke_0 - tr_yy_xyy[i];

        tr_0_0_z_yyz_xyz[i] = 2.0 * tr_yyzz_xyz[i] * tbe_0 + 2.0 * tr_yyz_xyzz[i] * tke_0 - tr_yy_xyz[i] - tr_yyz_xy[i];

        tr_0_0_z_yyz_xzz[i] = 2.0 * tr_yyzz_xzz[i] * tbe_0 + 2.0 * tr_yyz_xzzz[i] * tke_0 - tr_yy_xzz[i] - 2.0 * tr_yyz_xz[i];

        tr_0_0_z_yyz_yyy[i] = 2.0 * tr_yyzz_yyy[i] * tbe_0 + 2.0 * tr_yyz_yyyz[i] * tke_0 - tr_yy_yyy[i];

        tr_0_0_z_yyz_yyz[i] = 2.0 * tr_yyzz_yyz[i] * tbe_0 + 2.0 * tr_yyz_yyzz[i] * tke_0 - tr_yy_yyz[i] - tr_yyz_yy[i];

        tr_0_0_z_yyz_yzz[i] = 2.0 * tr_yyzz_yzz[i] * tbe_0 + 2.0 * tr_yyz_yzzz[i] * tke_0 - tr_yy_yzz[i] - 2.0 * tr_yyz_yz[i];

        tr_0_0_z_yyz_zzz[i] = 2.0 * tr_yyzz_zzz[i] * tbe_0 + 2.0 * tr_yyz_zzzz[i] * tke_0 - tr_yy_zzz[i] - 3.0 * tr_yyz_zz[i];
    }

    // Set up 280-290 components of targeted buffer : FF

    auto tr_0_0_z_yzz_xxx = pbuffer.data(idx_op_geom_010_ff + 280);

    auto tr_0_0_z_yzz_xxy = pbuffer.data(idx_op_geom_010_ff + 281);

    auto tr_0_0_z_yzz_xxz = pbuffer.data(idx_op_geom_010_ff + 282);

    auto tr_0_0_z_yzz_xyy = pbuffer.data(idx_op_geom_010_ff + 283);

    auto tr_0_0_z_yzz_xyz = pbuffer.data(idx_op_geom_010_ff + 284);

    auto tr_0_0_z_yzz_xzz = pbuffer.data(idx_op_geom_010_ff + 285);

    auto tr_0_0_z_yzz_yyy = pbuffer.data(idx_op_geom_010_ff + 286);

    auto tr_0_0_z_yzz_yyz = pbuffer.data(idx_op_geom_010_ff + 287);

    auto tr_0_0_z_yzz_yzz = pbuffer.data(idx_op_geom_010_ff + 288);

    auto tr_0_0_z_yzz_zzz = pbuffer.data(idx_op_geom_010_ff + 289);

    #pragma omp simd aligned(tr_0_0_z_yzz_xxx, tr_0_0_z_yzz_xxy, tr_0_0_z_yzz_xxz, tr_0_0_z_yzz_xyy, tr_0_0_z_yzz_xyz, tr_0_0_z_yzz_xzz, tr_0_0_z_yzz_yyy, tr_0_0_z_yzz_yyz, tr_0_0_z_yzz_yzz, tr_0_0_z_yzz_zzz, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_zzz, tr_yzz_xx, tr_yzz_xxxz, tr_yzz_xxyz, tr_yzz_xxzz, tr_yzz_xy, tr_yzz_xyyz, tr_yzz_xyzz, tr_yzz_xz, tr_yzz_xzzz, tr_yzz_yy, tr_yzz_yyyz, tr_yzz_yyzz, tr_yzz_yz, tr_yzz_yzzz, tr_yzz_zz, tr_yzz_zzzz, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzz_xxx[i] = 2.0 * tr_yzzz_xxx[i] * tbe_0 + 2.0 * tr_yzz_xxxz[i] * tke_0 - 2.0 * tr_yz_xxx[i];

        tr_0_0_z_yzz_xxy[i] = 2.0 * tr_yzzz_xxy[i] * tbe_0 + 2.0 * tr_yzz_xxyz[i] * tke_0 - 2.0 * tr_yz_xxy[i];

        tr_0_0_z_yzz_xxz[i] = 2.0 * tr_yzzz_xxz[i] * tbe_0 + 2.0 * tr_yzz_xxzz[i] * tke_0 - 2.0 * tr_yz_xxz[i] - tr_yzz_xx[i];

        tr_0_0_z_yzz_xyy[i] = 2.0 * tr_yzzz_xyy[i] * tbe_0 + 2.0 * tr_yzz_xyyz[i] * tke_0 - 2.0 * tr_yz_xyy[i];

        tr_0_0_z_yzz_xyz[i] = 2.0 * tr_yzzz_xyz[i] * tbe_0 + 2.0 * tr_yzz_xyzz[i] * tke_0 - 2.0 * tr_yz_xyz[i] - tr_yzz_xy[i];

        tr_0_0_z_yzz_xzz[i] = 2.0 * tr_yzzz_xzz[i] * tbe_0 + 2.0 * tr_yzz_xzzz[i] * tke_0 - 2.0 * tr_yz_xzz[i] - 2.0 * tr_yzz_xz[i];

        tr_0_0_z_yzz_yyy[i] = 2.0 * tr_yzzz_yyy[i] * tbe_0 + 2.0 * tr_yzz_yyyz[i] * tke_0 - 2.0 * tr_yz_yyy[i];

        tr_0_0_z_yzz_yyz[i] = 2.0 * tr_yzzz_yyz[i] * tbe_0 + 2.0 * tr_yzz_yyzz[i] * tke_0 - 2.0 * tr_yz_yyz[i] - tr_yzz_yy[i];

        tr_0_0_z_yzz_yzz[i] = 2.0 * tr_yzzz_yzz[i] * tbe_0 + 2.0 * tr_yzz_yzzz[i] * tke_0 - 2.0 * tr_yz_yzz[i] - 2.0 * tr_yzz_yz[i];

        tr_0_0_z_yzz_zzz[i] = 2.0 * tr_yzzz_zzz[i] * tbe_0 + 2.0 * tr_yzz_zzzz[i] * tke_0 - 2.0 * tr_yz_zzz[i] - 3.0 * tr_yzz_zz[i];
    }

    // Set up 290-300 components of targeted buffer : FF

    auto tr_0_0_z_zzz_xxx = pbuffer.data(idx_op_geom_010_ff + 290);

    auto tr_0_0_z_zzz_xxy = pbuffer.data(idx_op_geom_010_ff + 291);

    auto tr_0_0_z_zzz_xxz = pbuffer.data(idx_op_geom_010_ff + 292);

    auto tr_0_0_z_zzz_xyy = pbuffer.data(idx_op_geom_010_ff + 293);

    auto tr_0_0_z_zzz_xyz = pbuffer.data(idx_op_geom_010_ff + 294);

    auto tr_0_0_z_zzz_xzz = pbuffer.data(idx_op_geom_010_ff + 295);

    auto tr_0_0_z_zzz_yyy = pbuffer.data(idx_op_geom_010_ff + 296);

    auto tr_0_0_z_zzz_yyz = pbuffer.data(idx_op_geom_010_ff + 297);

    auto tr_0_0_z_zzz_yzz = pbuffer.data(idx_op_geom_010_ff + 298);

    auto tr_0_0_z_zzz_zzz = pbuffer.data(idx_op_geom_010_ff + 299);

    #pragma omp simd aligned(tr_0_0_z_zzz_xxx, tr_0_0_z_zzz_xxy, tr_0_0_z_zzz_xxz, tr_0_0_z_zzz_xyy, tr_0_0_z_zzz_xyz, tr_0_0_z_zzz_xzz, tr_0_0_z_zzz_yyy, tr_0_0_z_zzz_yyz, tr_0_0_z_zzz_yzz, tr_0_0_z_zzz_zzz, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_zzz, tr_zzz_xx, tr_zzz_xxxz, tr_zzz_xxyz, tr_zzz_xxzz, tr_zzz_xy, tr_zzz_xyyz, tr_zzz_xyzz, tr_zzz_xz, tr_zzz_xzzz, tr_zzz_yy, tr_zzz_yyyz, tr_zzz_yyzz, tr_zzz_yz, tr_zzz_yzzz, tr_zzz_zz, tr_zzz_zzzz, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzz_xxx[i] = 2.0 * tr_zzzz_xxx[i] * tbe_0 + 2.0 * tr_zzz_xxxz[i] * tke_0 - 3.0 * tr_zz_xxx[i];

        tr_0_0_z_zzz_xxy[i] = 2.0 * tr_zzzz_xxy[i] * tbe_0 + 2.0 * tr_zzz_xxyz[i] * tke_0 - 3.0 * tr_zz_xxy[i];

        tr_0_0_z_zzz_xxz[i] = 2.0 * tr_zzzz_xxz[i] * tbe_0 + 2.0 * tr_zzz_xxzz[i] * tke_0 - 3.0 * tr_zz_xxz[i] - tr_zzz_xx[i];

        tr_0_0_z_zzz_xyy[i] = 2.0 * tr_zzzz_xyy[i] * tbe_0 + 2.0 * tr_zzz_xyyz[i] * tke_0 - 3.0 * tr_zz_xyy[i];

        tr_0_0_z_zzz_xyz[i] = 2.0 * tr_zzzz_xyz[i] * tbe_0 + 2.0 * tr_zzz_xyzz[i] * tke_0 - 3.0 * tr_zz_xyz[i] - tr_zzz_xy[i];

        tr_0_0_z_zzz_xzz[i] = 2.0 * tr_zzzz_xzz[i] * tbe_0 + 2.0 * tr_zzz_xzzz[i] * tke_0 - 3.0 * tr_zz_xzz[i] - 2.0 * tr_zzz_xz[i];

        tr_0_0_z_zzz_yyy[i] = 2.0 * tr_zzzz_yyy[i] * tbe_0 + 2.0 * tr_zzz_yyyz[i] * tke_0 - 3.0 * tr_zz_yyy[i];

        tr_0_0_z_zzz_yyz[i] = 2.0 * tr_zzzz_yyz[i] * tbe_0 + 2.0 * tr_zzz_yyzz[i] * tke_0 - 3.0 * tr_zz_yyz[i] - tr_zzz_yy[i];

        tr_0_0_z_zzz_yzz[i] = 2.0 * tr_zzzz_yzz[i] * tbe_0 + 2.0 * tr_zzz_yzzz[i] * tke_0 - 3.0 * tr_zz_yzz[i] - 2.0 * tr_zzz_yz[i];

        tr_0_0_z_zzz_zzz[i] = 2.0 * tr_zzzz_zzz[i] * tbe_0 + 2.0 * tr_zzz_zzzz[i] * tke_0 - 3.0 * tr_zz_zzz[i] - 3.0 * tr_zzz_zz[i];
    }

}

} // t2cgeom namespace

