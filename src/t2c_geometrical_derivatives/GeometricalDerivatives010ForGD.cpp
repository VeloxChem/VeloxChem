#include "GeometricalDerivatives010ForGD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_010_gd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_010_gd,
                         const int idx_op_fd,
                         const int idx_op_gp,
                         const int idx_op_gf,
                         const int idx_op_hd,
                         const CSimdArray<double>& factors,
                         const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

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

    // Set up components of auxiliary buffer : GP

    auto tr_xxxx_x = pbuffer.data(idx_op_gp);

    auto tr_xxxx_y = pbuffer.data(idx_op_gp + 1);

    auto tr_xxxx_z = pbuffer.data(idx_op_gp + 2);

    auto tr_xxxy_x = pbuffer.data(idx_op_gp + 3);

    auto tr_xxxy_y = pbuffer.data(idx_op_gp + 4);

    auto tr_xxxy_z = pbuffer.data(idx_op_gp + 5);

    auto tr_xxxz_x = pbuffer.data(idx_op_gp + 6);

    auto tr_xxxz_y = pbuffer.data(idx_op_gp + 7);

    auto tr_xxxz_z = pbuffer.data(idx_op_gp + 8);

    auto tr_xxyy_x = pbuffer.data(idx_op_gp + 9);

    auto tr_xxyy_y = pbuffer.data(idx_op_gp + 10);

    auto tr_xxyy_z = pbuffer.data(idx_op_gp + 11);

    auto tr_xxyz_x = pbuffer.data(idx_op_gp + 12);

    auto tr_xxyz_y = pbuffer.data(idx_op_gp + 13);

    auto tr_xxyz_z = pbuffer.data(idx_op_gp + 14);

    auto tr_xxzz_x = pbuffer.data(idx_op_gp + 15);

    auto tr_xxzz_y = pbuffer.data(idx_op_gp + 16);

    auto tr_xxzz_z = pbuffer.data(idx_op_gp + 17);

    auto tr_xyyy_x = pbuffer.data(idx_op_gp + 18);

    auto tr_xyyy_y = pbuffer.data(idx_op_gp + 19);

    auto tr_xyyy_z = pbuffer.data(idx_op_gp + 20);

    auto tr_xyyz_x = pbuffer.data(idx_op_gp + 21);

    auto tr_xyyz_y = pbuffer.data(idx_op_gp + 22);

    auto tr_xyyz_z = pbuffer.data(idx_op_gp + 23);

    auto tr_xyzz_x = pbuffer.data(idx_op_gp + 24);

    auto tr_xyzz_y = pbuffer.data(idx_op_gp + 25);

    auto tr_xyzz_z = pbuffer.data(idx_op_gp + 26);

    auto tr_xzzz_x = pbuffer.data(idx_op_gp + 27);

    auto tr_xzzz_y = pbuffer.data(idx_op_gp + 28);

    auto tr_xzzz_z = pbuffer.data(idx_op_gp + 29);

    auto tr_yyyy_x = pbuffer.data(idx_op_gp + 30);

    auto tr_yyyy_y = pbuffer.data(idx_op_gp + 31);

    auto tr_yyyy_z = pbuffer.data(idx_op_gp + 32);

    auto tr_yyyz_x = pbuffer.data(idx_op_gp + 33);

    auto tr_yyyz_y = pbuffer.data(idx_op_gp + 34);

    auto tr_yyyz_z = pbuffer.data(idx_op_gp + 35);

    auto tr_yyzz_x = pbuffer.data(idx_op_gp + 36);

    auto tr_yyzz_y = pbuffer.data(idx_op_gp + 37);

    auto tr_yyzz_z = pbuffer.data(idx_op_gp + 38);

    auto tr_yzzz_x = pbuffer.data(idx_op_gp + 39);

    auto tr_yzzz_y = pbuffer.data(idx_op_gp + 40);

    auto tr_yzzz_z = pbuffer.data(idx_op_gp + 41);

    auto tr_zzzz_x = pbuffer.data(idx_op_gp + 42);

    auto tr_zzzz_y = pbuffer.data(idx_op_gp + 43);

    auto tr_zzzz_z = pbuffer.data(idx_op_gp + 44);

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

    // Set up 0-6 components of targeted buffer : GD

    auto tr_0_0_x_xxxx_xx = pbuffer.data(idx_op_geom_010_gd);

    auto tr_0_0_x_xxxx_xy = pbuffer.data(idx_op_geom_010_gd + 1);

    auto tr_0_0_x_xxxx_xz = pbuffer.data(idx_op_geom_010_gd + 2);

    auto tr_0_0_x_xxxx_yy = pbuffer.data(idx_op_geom_010_gd + 3);

    auto tr_0_0_x_xxxx_yz = pbuffer.data(idx_op_geom_010_gd + 4);

    auto tr_0_0_x_xxxx_zz = pbuffer.data(idx_op_geom_010_gd + 5);

    #pragma omp simd aligned(tr_0_0_x_xxxx_xx, tr_0_0_x_xxxx_xy, tr_0_0_x_xxxx_xz, tr_0_0_x_xxxx_yy, tr_0_0_x_xxxx_yz, tr_0_0_x_xxxx_zz, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_x, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_y, tr_xxxx_z, tr_xxxxx_xx, tr_xxxxx_xy, tr_xxxxx_xz, tr_xxxxx_yy, tr_xxxxx_yz, tr_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxx_xx[i] = 2.0 * tr_xxxxx_xx[i] * tbe_0 + 2.0 * tr_xxxx_xxx[i] * tke_0 - 4.0 * tr_xxx_xx[i] - 2.0 * tr_xxxx_x[i];

        tr_0_0_x_xxxx_xy[i] = 2.0 * tr_xxxxx_xy[i] * tbe_0 + 2.0 * tr_xxxx_xxy[i] * tke_0 - 4.0 * tr_xxx_xy[i] - tr_xxxx_y[i];

        tr_0_0_x_xxxx_xz[i] = 2.0 * tr_xxxxx_xz[i] * tbe_0 + 2.0 * tr_xxxx_xxz[i] * tke_0 - 4.0 * tr_xxx_xz[i] - tr_xxxx_z[i];

        tr_0_0_x_xxxx_yy[i] = 2.0 * tr_xxxxx_yy[i] * tbe_0 + 2.0 * tr_xxxx_xyy[i] * tke_0 - 4.0 * tr_xxx_yy[i];

        tr_0_0_x_xxxx_yz[i] = 2.0 * tr_xxxxx_yz[i] * tbe_0 + 2.0 * tr_xxxx_xyz[i] * tke_0 - 4.0 * tr_xxx_yz[i];

        tr_0_0_x_xxxx_zz[i] = 2.0 * tr_xxxxx_zz[i] * tbe_0 + 2.0 * tr_xxxx_xzz[i] * tke_0 - 4.0 * tr_xxx_zz[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto tr_0_0_x_xxxy_xx = pbuffer.data(idx_op_geom_010_gd + 6);

    auto tr_0_0_x_xxxy_xy = pbuffer.data(idx_op_geom_010_gd + 7);

    auto tr_0_0_x_xxxy_xz = pbuffer.data(idx_op_geom_010_gd + 8);

    auto tr_0_0_x_xxxy_yy = pbuffer.data(idx_op_geom_010_gd + 9);

    auto tr_0_0_x_xxxy_yz = pbuffer.data(idx_op_geom_010_gd + 10);

    auto tr_0_0_x_xxxy_zz = pbuffer.data(idx_op_geom_010_gd + 11);

    #pragma omp simd aligned(tr_0_0_x_xxxy_xx, tr_0_0_x_xxxy_xy, tr_0_0_x_xxxy_xz, tr_0_0_x_xxxy_yy, tr_0_0_x_xxxy_yz, tr_0_0_x_xxxy_zz, tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxy_x, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_y, tr_xxxy_z, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxy_xx[i] = 2.0 * tr_xxxxy_xx[i] * tbe_0 + 2.0 * tr_xxxy_xxx[i] * tke_0 - 3.0 * tr_xxy_xx[i] - 2.0 * tr_xxxy_x[i];

        tr_0_0_x_xxxy_xy[i] = 2.0 * tr_xxxxy_xy[i] * tbe_0 + 2.0 * tr_xxxy_xxy[i] * tke_0 - 3.0 * tr_xxy_xy[i] - tr_xxxy_y[i];

        tr_0_0_x_xxxy_xz[i] = 2.0 * tr_xxxxy_xz[i] * tbe_0 + 2.0 * tr_xxxy_xxz[i] * tke_0 - 3.0 * tr_xxy_xz[i] - tr_xxxy_z[i];

        tr_0_0_x_xxxy_yy[i] = 2.0 * tr_xxxxy_yy[i] * tbe_0 + 2.0 * tr_xxxy_xyy[i] * tke_0 - 3.0 * tr_xxy_yy[i];

        tr_0_0_x_xxxy_yz[i] = 2.0 * tr_xxxxy_yz[i] * tbe_0 + 2.0 * tr_xxxy_xyz[i] * tke_0 - 3.0 * tr_xxy_yz[i];

        tr_0_0_x_xxxy_zz[i] = 2.0 * tr_xxxxy_zz[i] * tbe_0 + 2.0 * tr_xxxy_xzz[i] * tke_0 - 3.0 * tr_xxy_zz[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto tr_0_0_x_xxxz_xx = pbuffer.data(idx_op_geom_010_gd + 12);

    auto tr_0_0_x_xxxz_xy = pbuffer.data(idx_op_geom_010_gd + 13);

    auto tr_0_0_x_xxxz_xz = pbuffer.data(idx_op_geom_010_gd + 14);

    auto tr_0_0_x_xxxz_yy = pbuffer.data(idx_op_geom_010_gd + 15);

    auto tr_0_0_x_xxxz_yz = pbuffer.data(idx_op_geom_010_gd + 16);

    auto tr_0_0_x_xxxz_zz = pbuffer.data(idx_op_geom_010_gd + 17);

    #pragma omp simd aligned(tr_0_0_x_xxxz_xx, tr_0_0_x_xxxz_xy, tr_0_0_x_xxxz_xz, tr_0_0_x_xxxz_yy, tr_0_0_x_xxxz_yz, tr_0_0_x_xxxz_zz, tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxz_x, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_y, tr_xxxz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxxz_xx[i] = 2.0 * tr_xxxxz_xx[i] * tbe_0 + 2.0 * tr_xxxz_xxx[i] * tke_0 - 3.0 * tr_xxz_xx[i] - 2.0 * tr_xxxz_x[i];

        tr_0_0_x_xxxz_xy[i] = 2.0 * tr_xxxxz_xy[i] * tbe_0 + 2.0 * tr_xxxz_xxy[i] * tke_0 - 3.0 * tr_xxz_xy[i] - tr_xxxz_y[i];

        tr_0_0_x_xxxz_xz[i] = 2.0 * tr_xxxxz_xz[i] * tbe_0 + 2.0 * tr_xxxz_xxz[i] * tke_0 - 3.0 * tr_xxz_xz[i] - tr_xxxz_z[i];

        tr_0_0_x_xxxz_yy[i] = 2.0 * tr_xxxxz_yy[i] * tbe_0 + 2.0 * tr_xxxz_xyy[i] * tke_0 - 3.0 * tr_xxz_yy[i];

        tr_0_0_x_xxxz_yz[i] = 2.0 * tr_xxxxz_yz[i] * tbe_0 + 2.0 * tr_xxxz_xyz[i] * tke_0 - 3.0 * tr_xxz_yz[i];

        tr_0_0_x_xxxz_zz[i] = 2.0 * tr_xxxxz_zz[i] * tbe_0 + 2.0 * tr_xxxz_xzz[i] * tke_0 - 3.0 * tr_xxz_zz[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto tr_0_0_x_xxyy_xx = pbuffer.data(idx_op_geom_010_gd + 18);

    auto tr_0_0_x_xxyy_xy = pbuffer.data(idx_op_geom_010_gd + 19);

    auto tr_0_0_x_xxyy_xz = pbuffer.data(idx_op_geom_010_gd + 20);

    auto tr_0_0_x_xxyy_yy = pbuffer.data(idx_op_geom_010_gd + 21);

    auto tr_0_0_x_xxyy_yz = pbuffer.data(idx_op_geom_010_gd + 22);

    auto tr_0_0_x_xxyy_zz = pbuffer.data(idx_op_geom_010_gd + 23);

    #pragma omp simd aligned(tr_0_0_x_xxyy_xx, tr_0_0_x_xxyy_xy, tr_0_0_x_xxyy_xz, tr_0_0_x_xxyy_yy, tr_0_0_x_xxyy_yz, tr_0_0_x_xxyy_zz, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxyy_x, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_y, tr_xxyy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyy_xx[i] = 2.0 * tr_xxxyy_xx[i] * tbe_0 + 2.0 * tr_xxyy_xxx[i] * tke_0 - 2.0 * tr_xyy_xx[i] - 2.0 * tr_xxyy_x[i];

        tr_0_0_x_xxyy_xy[i] = 2.0 * tr_xxxyy_xy[i] * tbe_0 + 2.0 * tr_xxyy_xxy[i] * tke_0 - 2.0 * tr_xyy_xy[i] - tr_xxyy_y[i];

        tr_0_0_x_xxyy_xz[i] = 2.0 * tr_xxxyy_xz[i] * tbe_0 + 2.0 * tr_xxyy_xxz[i] * tke_0 - 2.0 * tr_xyy_xz[i] - tr_xxyy_z[i];

        tr_0_0_x_xxyy_yy[i] = 2.0 * tr_xxxyy_yy[i] * tbe_0 + 2.0 * tr_xxyy_xyy[i] * tke_0 - 2.0 * tr_xyy_yy[i];

        tr_0_0_x_xxyy_yz[i] = 2.0 * tr_xxxyy_yz[i] * tbe_0 + 2.0 * tr_xxyy_xyz[i] * tke_0 - 2.0 * tr_xyy_yz[i];

        tr_0_0_x_xxyy_zz[i] = 2.0 * tr_xxxyy_zz[i] * tbe_0 + 2.0 * tr_xxyy_xzz[i] * tke_0 - 2.0 * tr_xyy_zz[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto tr_0_0_x_xxyz_xx = pbuffer.data(idx_op_geom_010_gd + 24);

    auto tr_0_0_x_xxyz_xy = pbuffer.data(idx_op_geom_010_gd + 25);

    auto tr_0_0_x_xxyz_xz = pbuffer.data(idx_op_geom_010_gd + 26);

    auto tr_0_0_x_xxyz_yy = pbuffer.data(idx_op_geom_010_gd + 27);

    auto tr_0_0_x_xxyz_yz = pbuffer.data(idx_op_geom_010_gd + 28);

    auto tr_0_0_x_xxyz_zz = pbuffer.data(idx_op_geom_010_gd + 29);

    #pragma omp simd aligned(tr_0_0_x_xxyz_xx, tr_0_0_x_xxyz_xy, tr_0_0_x_xxyz_xz, tr_0_0_x_xxyz_yy, tr_0_0_x_xxyz_yz, tr_0_0_x_xxyz_zz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_x, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxyz_xx[i] = 2.0 * tr_xxxyz_xx[i] * tbe_0 + 2.0 * tr_xxyz_xxx[i] * tke_0 - 2.0 * tr_xyz_xx[i] - 2.0 * tr_xxyz_x[i];

        tr_0_0_x_xxyz_xy[i] = 2.0 * tr_xxxyz_xy[i] * tbe_0 + 2.0 * tr_xxyz_xxy[i] * tke_0 - 2.0 * tr_xyz_xy[i] - tr_xxyz_y[i];

        tr_0_0_x_xxyz_xz[i] = 2.0 * tr_xxxyz_xz[i] * tbe_0 + 2.0 * tr_xxyz_xxz[i] * tke_0 - 2.0 * tr_xyz_xz[i] - tr_xxyz_z[i];

        tr_0_0_x_xxyz_yy[i] = 2.0 * tr_xxxyz_yy[i] * tbe_0 + 2.0 * tr_xxyz_xyy[i] * tke_0 - 2.0 * tr_xyz_yy[i];

        tr_0_0_x_xxyz_yz[i] = 2.0 * tr_xxxyz_yz[i] * tbe_0 + 2.0 * tr_xxyz_xyz[i] * tke_0 - 2.0 * tr_xyz_yz[i];

        tr_0_0_x_xxyz_zz[i] = 2.0 * tr_xxxyz_zz[i] * tbe_0 + 2.0 * tr_xxyz_xzz[i] * tke_0 - 2.0 * tr_xyz_zz[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto tr_0_0_x_xxzz_xx = pbuffer.data(idx_op_geom_010_gd + 30);

    auto tr_0_0_x_xxzz_xy = pbuffer.data(idx_op_geom_010_gd + 31);

    auto tr_0_0_x_xxzz_xz = pbuffer.data(idx_op_geom_010_gd + 32);

    auto tr_0_0_x_xxzz_yy = pbuffer.data(idx_op_geom_010_gd + 33);

    auto tr_0_0_x_xxzz_yz = pbuffer.data(idx_op_geom_010_gd + 34);

    auto tr_0_0_x_xxzz_zz = pbuffer.data(idx_op_geom_010_gd + 35);

    #pragma omp simd aligned(tr_0_0_x_xxzz_xx, tr_0_0_x_xxzz_xy, tr_0_0_x_xxzz_xz, tr_0_0_x_xxzz_yy, tr_0_0_x_xxzz_yz, tr_0_0_x_xxzz_zz, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxzz_x, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_y, tr_xxzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xxzz_xx[i] = 2.0 * tr_xxxzz_xx[i] * tbe_0 + 2.0 * tr_xxzz_xxx[i] * tke_0 - 2.0 * tr_xzz_xx[i] - 2.0 * tr_xxzz_x[i];

        tr_0_0_x_xxzz_xy[i] = 2.0 * tr_xxxzz_xy[i] * tbe_0 + 2.0 * tr_xxzz_xxy[i] * tke_0 - 2.0 * tr_xzz_xy[i] - tr_xxzz_y[i];

        tr_0_0_x_xxzz_xz[i] = 2.0 * tr_xxxzz_xz[i] * tbe_0 + 2.0 * tr_xxzz_xxz[i] * tke_0 - 2.0 * tr_xzz_xz[i] - tr_xxzz_z[i];

        tr_0_0_x_xxzz_yy[i] = 2.0 * tr_xxxzz_yy[i] * tbe_0 + 2.0 * tr_xxzz_xyy[i] * tke_0 - 2.0 * tr_xzz_yy[i];

        tr_0_0_x_xxzz_yz[i] = 2.0 * tr_xxxzz_yz[i] * tbe_0 + 2.0 * tr_xxzz_xyz[i] * tke_0 - 2.0 * tr_xzz_yz[i];

        tr_0_0_x_xxzz_zz[i] = 2.0 * tr_xxxzz_zz[i] * tbe_0 + 2.0 * tr_xxzz_xzz[i] * tke_0 - 2.0 * tr_xzz_zz[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto tr_0_0_x_xyyy_xx = pbuffer.data(idx_op_geom_010_gd + 36);

    auto tr_0_0_x_xyyy_xy = pbuffer.data(idx_op_geom_010_gd + 37);

    auto tr_0_0_x_xyyy_xz = pbuffer.data(idx_op_geom_010_gd + 38);

    auto tr_0_0_x_xyyy_yy = pbuffer.data(idx_op_geom_010_gd + 39);

    auto tr_0_0_x_xyyy_yz = pbuffer.data(idx_op_geom_010_gd + 40);

    auto tr_0_0_x_xyyy_zz = pbuffer.data(idx_op_geom_010_gd + 41);

    #pragma omp simd aligned(tr_0_0_x_xyyy_xx, tr_0_0_x_xyyy_xy, tr_0_0_x_xyyy_xz, tr_0_0_x_xyyy_yy, tr_0_0_x_xyyy_yz, tr_0_0_x_xyyy_zz, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_xyyy_x, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_y, tr_xyyy_z, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyy_xx[i] = 2.0 * tr_xxyyy_xx[i] * tbe_0 + 2.0 * tr_xyyy_xxx[i] * tke_0 - tr_yyy_xx[i] - 2.0 * tr_xyyy_x[i];

        tr_0_0_x_xyyy_xy[i] = 2.0 * tr_xxyyy_xy[i] * tbe_0 + 2.0 * tr_xyyy_xxy[i] * tke_0 - tr_yyy_xy[i] - tr_xyyy_y[i];

        tr_0_0_x_xyyy_xz[i] = 2.0 * tr_xxyyy_xz[i] * tbe_0 + 2.0 * tr_xyyy_xxz[i] * tke_0 - tr_yyy_xz[i] - tr_xyyy_z[i];

        tr_0_0_x_xyyy_yy[i] = 2.0 * tr_xxyyy_yy[i] * tbe_0 + 2.0 * tr_xyyy_xyy[i] * tke_0 - tr_yyy_yy[i];

        tr_0_0_x_xyyy_yz[i] = 2.0 * tr_xxyyy_yz[i] * tbe_0 + 2.0 * tr_xyyy_xyz[i] * tke_0 - tr_yyy_yz[i];

        tr_0_0_x_xyyy_zz[i] = 2.0 * tr_xxyyy_zz[i] * tbe_0 + 2.0 * tr_xyyy_xzz[i] * tke_0 - tr_yyy_zz[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto tr_0_0_x_xyyz_xx = pbuffer.data(idx_op_geom_010_gd + 42);

    auto tr_0_0_x_xyyz_xy = pbuffer.data(idx_op_geom_010_gd + 43);

    auto tr_0_0_x_xyyz_xz = pbuffer.data(idx_op_geom_010_gd + 44);

    auto tr_0_0_x_xyyz_yy = pbuffer.data(idx_op_geom_010_gd + 45);

    auto tr_0_0_x_xyyz_yz = pbuffer.data(idx_op_geom_010_gd + 46);

    auto tr_0_0_x_xyyz_zz = pbuffer.data(idx_op_geom_010_gd + 47);

    #pragma omp simd aligned(tr_0_0_x_xyyz_xx, tr_0_0_x_xyyz_xy, tr_0_0_x_xyyz_xz, tr_0_0_x_xyyz_yy, tr_0_0_x_xyyz_yz, tr_0_0_x_xyyz_zz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xyyz_x, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_z, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyyz_xx[i] = 2.0 * tr_xxyyz_xx[i] * tbe_0 + 2.0 * tr_xyyz_xxx[i] * tke_0 - tr_yyz_xx[i] - 2.0 * tr_xyyz_x[i];

        tr_0_0_x_xyyz_xy[i] = 2.0 * tr_xxyyz_xy[i] * tbe_0 + 2.0 * tr_xyyz_xxy[i] * tke_0 - tr_yyz_xy[i] - tr_xyyz_y[i];

        tr_0_0_x_xyyz_xz[i] = 2.0 * tr_xxyyz_xz[i] * tbe_0 + 2.0 * tr_xyyz_xxz[i] * tke_0 - tr_yyz_xz[i] - tr_xyyz_z[i];

        tr_0_0_x_xyyz_yy[i] = 2.0 * tr_xxyyz_yy[i] * tbe_0 + 2.0 * tr_xyyz_xyy[i] * tke_0 - tr_yyz_yy[i];

        tr_0_0_x_xyyz_yz[i] = 2.0 * tr_xxyyz_yz[i] * tbe_0 + 2.0 * tr_xyyz_xyz[i] * tke_0 - tr_yyz_yz[i];

        tr_0_0_x_xyyz_zz[i] = 2.0 * tr_xxyyz_zz[i] * tbe_0 + 2.0 * tr_xyyz_xzz[i] * tke_0 - tr_yyz_zz[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto tr_0_0_x_xyzz_xx = pbuffer.data(idx_op_geom_010_gd + 48);

    auto tr_0_0_x_xyzz_xy = pbuffer.data(idx_op_geom_010_gd + 49);

    auto tr_0_0_x_xyzz_xz = pbuffer.data(idx_op_geom_010_gd + 50);

    auto tr_0_0_x_xyzz_yy = pbuffer.data(idx_op_geom_010_gd + 51);

    auto tr_0_0_x_xyzz_yz = pbuffer.data(idx_op_geom_010_gd + 52);

    auto tr_0_0_x_xyzz_zz = pbuffer.data(idx_op_geom_010_gd + 53);

    #pragma omp simd aligned(tr_0_0_x_xyzz_xx, tr_0_0_x_xyzz_xy, tr_0_0_x_xyzz_xz, tr_0_0_x_xyzz_yy, tr_0_0_x_xyzz_yz, tr_0_0_x_xyzz_zz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xyzz_x, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_z, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xyzz_xx[i] = 2.0 * tr_xxyzz_xx[i] * tbe_0 + 2.0 * tr_xyzz_xxx[i] * tke_0 - tr_yzz_xx[i] - 2.0 * tr_xyzz_x[i];

        tr_0_0_x_xyzz_xy[i] = 2.0 * tr_xxyzz_xy[i] * tbe_0 + 2.0 * tr_xyzz_xxy[i] * tke_0 - tr_yzz_xy[i] - tr_xyzz_y[i];

        tr_0_0_x_xyzz_xz[i] = 2.0 * tr_xxyzz_xz[i] * tbe_0 + 2.0 * tr_xyzz_xxz[i] * tke_0 - tr_yzz_xz[i] - tr_xyzz_z[i];

        tr_0_0_x_xyzz_yy[i] = 2.0 * tr_xxyzz_yy[i] * tbe_0 + 2.0 * tr_xyzz_xyy[i] * tke_0 - tr_yzz_yy[i];

        tr_0_0_x_xyzz_yz[i] = 2.0 * tr_xxyzz_yz[i] * tbe_0 + 2.0 * tr_xyzz_xyz[i] * tke_0 - tr_yzz_yz[i];

        tr_0_0_x_xyzz_zz[i] = 2.0 * tr_xxyzz_zz[i] * tbe_0 + 2.0 * tr_xyzz_xzz[i] * tke_0 - tr_yzz_zz[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto tr_0_0_x_xzzz_xx = pbuffer.data(idx_op_geom_010_gd + 54);

    auto tr_0_0_x_xzzz_xy = pbuffer.data(idx_op_geom_010_gd + 55);

    auto tr_0_0_x_xzzz_xz = pbuffer.data(idx_op_geom_010_gd + 56);

    auto tr_0_0_x_xzzz_yy = pbuffer.data(idx_op_geom_010_gd + 57);

    auto tr_0_0_x_xzzz_yz = pbuffer.data(idx_op_geom_010_gd + 58);

    auto tr_0_0_x_xzzz_zz = pbuffer.data(idx_op_geom_010_gd + 59);

    #pragma omp simd aligned(tr_0_0_x_xzzz_xx, tr_0_0_x_xzzz_xy, tr_0_0_x_xzzz_xz, tr_0_0_x_xzzz_yy, tr_0_0_x_xzzz_yz, tr_0_0_x_xzzz_zz, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_xzzz_x, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_y, tr_xzzz_z, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_xzzz_xx[i] = 2.0 * tr_xxzzz_xx[i] * tbe_0 + 2.0 * tr_xzzz_xxx[i] * tke_0 - tr_zzz_xx[i] - 2.0 * tr_xzzz_x[i];

        tr_0_0_x_xzzz_xy[i] = 2.0 * tr_xxzzz_xy[i] * tbe_0 + 2.0 * tr_xzzz_xxy[i] * tke_0 - tr_zzz_xy[i] - tr_xzzz_y[i];

        tr_0_0_x_xzzz_xz[i] = 2.0 * tr_xxzzz_xz[i] * tbe_0 + 2.0 * tr_xzzz_xxz[i] * tke_0 - tr_zzz_xz[i] - tr_xzzz_z[i];

        tr_0_0_x_xzzz_yy[i] = 2.0 * tr_xxzzz_yy[i] * tbe_0 + 2.0 * tr_xzzz_xyy[i] * tke_0 - tr_zzz_yy[i];

        tr_0_0_x_xzzz_yz[i] = 2.0 * tr_xxzzz_yz[i] * tbe_0 + 2.0 * tr_xzzz_xyz[i] * tke_0 - tr_zzz_yz[i];

        tr_0_0_x_xzzz_zz[i] = 2.0 * tr_xxzzz_zz[i] * tbe_0 + 2.0 * tr_xzzz_xzz[i] * tke_0 - tr_zzz_zz[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto tr_0_0_x_yyyy_xx = pbuffer.data(idx_op_geom_010_gd + 60);

    auto tr_0_0_x_yyyy_xy = pbuffer.data(idx_op_geom_010_gd + 61);

    auto tr_0_0_x_yyyy_xz = pbuffer.data(idx_op_geom_010_gd + 62);

    auto tr_0_0_x_yyyy_yy = pbuffer.data(idx_op_geom_010_gd + 63);

    auto tr_0_0_x_yyyy_yz = pbuffer.data(idx_op_geom_010_gd + 64);

    auto tr_0_0_x_yyyy_zz = pbuffer.data(idx_op_geom_010_gd + 65);

    #pragma omp simd aligned(tr_0_0_x_yyyy_xx, tr_0_0_x_yyyy_xy, tr_0_0_x_yyyy_xz, tr_0_0_x_yyyy_yy, tr_0_0_x_yyyy_yz, tr_0_0_x_yyyy_zz, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, tr_yyyy_x, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyy_xx[i] = 2.0 * tr_xyyyy_xx[i] * tbe_0 + 2.0 * tr_yyyy_xxx[i] * tke_0 - 2.0 * tr_yyyy_x[i];

        tr_0_0_x_yyyy_xy[i] = 2.0 * tr_xyyyy_xy[i] * tbe_0 + 2.0 * tr_yyyy_xxy[i] * tke_0 - tr_yyyy_y[i];

        tr_0_0_x_yyyy_xz[i] = 2.0 * tr_xyyyy_xz[i] * tbe_0 + 2.0 * tr_yyyy_xxz[i] * tke_0 - tr_yyyy_z[i];

        tr_0_0_x_yyyy_yy[i] = 2.0 * tr_xyyyy_yy[i] * tbe_0 + 2.0 * tr_yyyy_xyy[i] * tke_0;

        tr_0_0_x_yyyy_yz[i] = 2.0 * tr_xyyyy_yz[i] * tbe_0 + 2.0 * tr_yyyy_xyz[i] * tke_0;

        tr_0_0_x_yyyy_zz[i] = 2.0 * tr_xyyyy_zz[i] * tbe_0 + 2.0 * tr_yyyy_xzz[i] * tke_0;
    }

    // Set up 66-72 components of targeted buffer : GD

    auto tr_0_0_x_yyyz_xx = pbuffer.data(idx_op_geom_010_gd + 66);

    auto tr_0_0_x_yyyz_xy = pbuffer.data(idx_op_geom_010_gd + 67);

    auto tr_0_0_x_yyyz_xz = pbuffer.data(idx_op_geom_010_gd + 68);

    auto tr_0_0_x_yyyz_yy = pbuffer.data(idx_op_geom_010_gd + 69);

    auto tr_0_0_x_yyyz_yz = pbuffer.data(idx_op_geom_010_gd + 70);

    auto tr_0_0_x_yyyz_zz = pbuffer.data(idx_op_geom_010_gd + 71);

    #pragma omp simd aligned(tr_0_0_x_yyyz_xx, tr_0_0_x_yyyz_xy, tr_0_0_x_yyyz_xz, tr_0_0_x_yyyz_yy, tr_0_0_x_yyyz_yz, tr_0_0_x_yyyz_zz, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_yyyz_x, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_y, tr_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyyz_xx[i] = 2.0 * tr_xyyyz_xx[i] * tbe_0 + 2.0 * tr_yyyz_xxx[i] * tke_0 - 2.0 * tr_yyyz_x[i];

        tr_0_0_x_yyyz_xy[i] = 2.0 * tr_xyyyz_xy[i] * tbe_0 + 2.0 * tr_yyyz_xxy[i] * tke_0 - tr_yyyz_y[i];

        tr_0_0_x_yyyz_xz[i] = 2.0 * tr_xyyyz_xz[i] * tbe_0 + 2.0 * tr_yyyz_xxz[i] * tke_0 - tr_yyyz_z[i];

        tr_0_0_x_yyyz_yy[i] = 2.0 * tr_xyyyz_yy[i] * tbe_0 + 2.0 * tr_yyyz_xyy[i] * tke_0;

        tr_0_0_x_yyyz_yz[i] = 2.0 * tr_xyyyz_yz[i] * tbe_0 + 2.0 * tr_yyyz_xyz[i] * tke_0;

        tr_0_0_x_yyyz_zz[i] = 2.0 * tr_xyyyz_zz[i] * tbe_0 + 2.0 * tr_yyyz_xzz[i] * tke_0;
    }

    // Set up 72-78 components of targeted buffer : GD

    auto tr_0_0_x_yyzz_xx = pbuffer.data(idx_op_geom_010_gd + 72);

    auto tr_0_0_x_yyzz_xy = pbuffer.data(idx_op_geom_010_gd + 73);

    auto tr_0_0_x_yyzz_xz = pbuffer.data(idx_op_geom_010_gd + 74);

    auto tr_0_0_x_yyzz_yy = pbuffer.data(idx_op_geom_010_gd + 75);

    auto tr_0_0_x_yyzz_yz = pbuffer.data(idx_op_geom_010_gd + 76);

    auto tr_0_0_x_yyzz_zz = pbuffer.data(idx_op_geom_010_gd + 77);

    #pragma omp simd aligned(tr_0_0_x_yyzz_xx, tr_0_0_x_yyzz_xy, tr_0_0_x_yyzz_xz, tr_0_0_x_yyzz_yy, tr_0_0_x_yyzz_yz, tr_0_0_x_yyzz_zz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_yyzz_x, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_y, tr_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yyzz_xx[i] = 2.0 * tr_xyyzz_xx[i] * tbe_0 + 2.0 * tr_yyzz_xxx[i] * tke_0 - 2.0 * tr_yyzz_x[i];

        tr_0_0_x_yyzz_xy[i] = 2.0 * tr_xyyzz_xy[i] * tbe_0 + 2.0 * tr_yyzz_xxy[i] * tke_0 - tr_yyzz_y[i];

        tr_0_0_x_yyzz_xz[i] = 2.0 * tr_xyyzz_xz[i] * tbe_0 + 2.0 * tr_yyzz_xxz[i] * tke_0 - tr_yyzz_z[i];

        tr_0_0_x_yyzz_yy[i] = 2.0 * tr_xyyzz_yy[i] * tbe_0 + 2.0 * tr_yyzz_xyy[i] * tke_0;

        tr_0_0_x_yyzz_yz[i] = 2.0 * tr_xyyzz_yz[i] * tbe_0 + 2.0 * tr_yyzz_xyz[i] * tke_0;

        tr_0_0_x_yyzz_zz[i] = 2.0 * tr_xyyzz_zz[i] * tbe_0 + 2.0 * tr_yyzz_xzz[i] * tke_0;
    }

    // Set up 78-84 components of targeted buffer : GD

    auto tr_0_0_x_yzzz_xx = pbuffer.data(idx_op_geom_010_gd + 78);

    auto tr_0_0_x_yzzz_xy = pbuffer.data(idx_op_geom_010_gd + 79);

    auto tr_0_0_x_yzzz_xz = pbuffer.data(idx_op_geom_010_gd + 80);

    auto tr_0_0_x_yzzz_yy = pbuffer.data(idx_op_geom_010_gd + 81);

    auto tr_0_0_x_yzzz_yz = pbuffer.data(idx_op_geom_010_gd + 82);

    auto tr_0_0_x_yzzz_zz = pbuffer.data(idx_op_geom_010_gd + 83);

    #pragma omp simd aligned(tr_0_0_x_yzzz_xx, tr_0_0_x_yzzz_xy, tr_0_0_x_yzzz_xz, tr_0_0_x_yzzz_yy, tr_0_0_x_yzzz_yz, tr_0_0_x_yzzz_zz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_yzzz_x, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_yzzz_xx[i] = 2.0 * tr_xyzzz_xx[i] * tbe_0 + 2.0 * tr_yzzz_xxx[i] * tke_0 - 2.0 * tr_yzzz_x[i];

        tr_0_0_x_yzzz_xy[i] = 2.0 * tr_xyzzz_xy[i] * tbe_0 + 2.0 * tr_yzzz_xxy[i] * tke_0 - tr_yzzz_y[i];

        tr_0_0_x_yzzz_xz[i] = 2.0 * tr_xyzzz_xz[i] * tbe_0 + 2.0 * tr_yzzz_xxz[i] * tke_0 - tr_yzzz_z[i];

        tr_0_0_x_yzzz_yy[i] = 2.0 * tr_xyzzz_yy[i] * tbe_0 + 2.0 * tr_yzzz_xyy[i] * tke_0;

        tr_0_0_x_yzzz_yz[i] = 2.0 * tr_xyzzz_yz[i] * tbe_0 + 2.0 * tr_yzzz_xyz[i] * tke_0;

        tr_0_0_x_yzzz_zz[i] = 2.0 * tr_xyzzz_zz[i] * tbe_0 + 2.0 * tr_yzzz_xzz[i] * tke_0;
    }

    // Set up 84-90 components of targeted buffer : GD

    auto tr_0_0_x_zzzz_xx = pbuffer.data(idx_op_geom_010_gd + 84);

    auto tr_0_0_x_zzzz_xy = pbuffer.data(idx_op_geom_010_gd + 85);

    auto tr_0_0_x_zzzz_xz = pbuffer.data(idx_op_geom_010_gd + 86);

    auto tr_0_0_x_zzzz_yy = pbuffer.data(idx_op_geom_010_gd + 87);

    auto tr_0_0_x_zzzz_yz = pbuffer.data(idx_op_geom_010_gd + 88);

    auto tr_0_0_x_zzzz_zz = pbuffer.data(idx_op_geom_010_gd + 89);

    #pragma omp simd aligned(tr_0_0_x_zzzz_xx, tr_0_0_x_zzzz_xy, tr_0_0_x_zzzz_xz, tr_0_0_x_zzzz_yy, tr_0_0_x_zzzz_yz, tr_0_0_x_zzzz_zz, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, tr_zzzz_x, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_x_zzzz_xx[i] = 2.0 * tr_xzzzz_xx[i] * tbe_0 + 2.0 * tr_zzzz_xxx[i] * tke_0 - 2.0 * tr_zzzz_x[i];

        tr_0_0_x_zzzz_xy[i] = 2.0 * tr_xzzzz_xy[i] * tbe_0 + 2.0 * tr_zzzz_xxy[i] * tke_0 - tr_zzzz_y[i];

        tr_0_0_x_zzzz_xz[i] = 2.0 * tr_xzzzz_xz[i] * tbe_0 + 2.0 * tr_zzzz_xxz[i] * tke_0 - tr_zzzz_z[i];

        tr_0_0_x_zzzz_yy[i] = 2.0 * tr_xzzzz_yy[i] * tbe_0 + 2.0 * tr_zzzz_xyy[i] * tke_0;

        tr_0_0_x_zzzz_yz[i] = 2.0 * tr_xzzzz_yz[i] * tbe_0 + 2.0 * tr_zzzz_xyz[i] * tke_0;

        tr_0_0_x_zzzz_zz[i] = 2.0 * tr_xzzzz_zz[i] * tbe_0 + 2.0 * tr_zzzz_xzz[i] * tke_0;
    }

    // Set up 90-96 components of targeted buffer : GD

    auto tr_0_0_y_xxxx_xx = pbuffer.data(idx_op_geom_010_gd + 90);

    auto tr_0_0_y_xxxx_xy = pbuffer.data(idx_op_geom_010_gd + 91);

    auto tr_0_0_y_xxxx_xz = pbuffer.data(idx_op_geom_010_gd + 92);

    auto tr_0_0_y_xxxx_yy = pbuffer.data(idx_op_geom_010_gd + 93);

    auto tr_0_0_y_xxxx_yz = pbuffer.data(idx_op_geom_010_gd + 94);

    auto tr_0_0_y_xxxx_zz = pbuffer.data(idx_op_geom_010_gd + 95);

    #pragma omp simd aligned(tr_0_0_y_xxxx_xx, tr_0_0_y_xxxx_xy, tr_0_0_y_xxxx_xz, tr_0_0_y_xxxx_yy, tr_0_0_y_xxxx_yz, tr_0_0_y_xxxx_zz, tr_xxxx_x, tr_xxxx_xxy, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_y, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_z, tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxx_xx[i] = 2.0 * tr_xxxxy_xx[i] * tbe_0 + 2.0 * tr_xxxx_xxy[i] * tke_0;

        tr_0_0_y_xxxx_xy[i] = 2.0 * tr_xxxxy_xy[i] * tbe_0 + 2.0 * tr_xxxx_xyy[i] * tke_0 - tr_xxxx_x[i];

        tr_0_0_y_xxxx_xz[i] = 2.0 * tr_xxxxy_xz[i] * tbe_0 + 2.0 * tr_xxxx_xyz[i] * tke_0;

        tr_0_0_y_xxxx_yy[i] = 2.0 * tr_xxxxy_yy[i] * tbe_0 + 2.0 * tr_xxxx_yyy[i] * tke_0 - 2.0 * tr_xxxx_y[i];

        tr_0_0_y_xxxx_yz[i] = 2.0 * tr_xxxxy_yz[i] * tbe_0 + 2.0 * tr_xxxx_yyz[i] * tke_0 - tr_xxxx_z[i];

        tr_0_0_y_xxxx_zz[i] = 2.0 * tr_xxxxy_zz[i] * tbe_0 + 2.0 * tr_xxxx_yzz[i] * tke_0;
    }

    // Set up 96-102 components of targeted buffer : GD

    auto tr_0_0_y_xxxy_xx = pbuffer.data(idx_op_geom_010_gd + 96);

    auto tr_0_0_y_xxxy_xy = pbuffer.data(idx_op_geom_010_gd + 97);

    auto tr_0_0_y_xxxy_xz = pbuffer.data(idx_op_geom_010_gd + 98);

    auto tr_0_0_y_xxxy_yy = pbuffer.data(idx_op_geom_010_gd + 99);

    auto tr_0_0_y_xxxy_yz = pbuffer.data(idx_op_geom_010_gd + 100);

    auto tr_0_0_y_xxxy_zz = pbuffer.data(idx_op_geom_010_gd + 101);

    #pragma omp simd aligned(tr_0_0_y_xxxy_xx, tr_0_0_y_xxxy_xy, tr_0_0_y_xxxy_xz, tr_0_0_y_xxxy_yy, tr_0_0_y_xxxy_yz, tr_0_0_y_xxxy_zz, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxy_x, tr_xxxy_xxy, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_y, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_z, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxy_xx[i] = 2.0 * tr_xxxyy_xx[i] * tbe_0 + 2.0 * tr_xxxy_xxy[i] * tke_0 - tr_xxx_xx[i];

        tr_0_0_y_xxxy_xy[i] = 2.0 * tr_xxxyy_xy[i] * tbe_0 + 2.0 * tr_xxxy_xyy[i] * tke_0 - tr_xxx_xy[i] - tr_xxxy_x[i];

        tr_0_0_y_xxxy_xz[i] = 2.0 * tr_xxxyy_xz[i] * tbe_0 + 2.0 * tr_xxxy_xyz[i] * tke_0 - tr_xxx_xz[i];

        tr_0_0_y_xxxy_yy[i] = 2.0 * tr_xxxyy_yy[i] * tbe_0 + 2.0 * tr_xxxy_yyy[i] * tke_0 - tr_xxx_yy[i] - 2.0 * tr_xxxy_y[i];

        tr_0_0_y_xxxy_yz[i] = 2.0 * tr_xxxyy_yz[i] * tbe_0 + 2.0 * tr_xxxy_yyz[i] * tke_0 - tr_xxx_yz[i] - tr_xxxy_z[i];

        tr_0_0_y_xxxy_zz[i] = 2.0 * tr_xxxyy_zz[i] * tbe_0 + 2.0 * tr_xxxy_yzz[i] * tke_0 - tr_xxx_zz[i];
    }

    // Set up 102-108 components of targeted buffer : GD

    auto tr_0_0_y_xxxz_xx = pbuffer.data(idx_op_geom_010_gd + 102);

    auto tr_0_0_y_xxxz_xy = pbuffer.data(idx_op_geom_010_gd + 103);

    auto tr_0_0_y_xxxz_xz = pbuffer.data(idx_op_geom_010_gd + 104);

    auto tr_0_0_y_xxxz_yy = pbuffer.data(idx_op_geom_010_gd + 105);

    auto tr_0_0_y_xxxz_yz = pbuffer.data(idx_op_geom_010_gd + 106);

    auto tr_0_0_y_xxxz_zz = pbuffer.data(idx_op_geom_010_gd + 107);

    #pragma omp simd aligned(tr_0_0_y_xxxz_xx, tr_0_0_y_xxxz_xy, tr_0_0_y_xxxz_xz, tr_0_0_y_xxxz_yy, tr_0_0_y_xxxz_yz, tr_0_0_y_xxxz_zz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxxz_x, tr_xxxz_xxy, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_y, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxxz_xx[i] = 2.0 * tr_xxxyz_xx[i] * tbe_0 + 2.0 * tr_xxxz_xxy[i] * tke_0;

        tr_0_0_y_xxxz_xy[i] = 2.0 * tr_xxxyz_xy[i] * tbe_0 + 2.0 * tr_xxxz_xyy[i] * tke_0 - tr_xxxz_x[i];

        tr_0_0_y_xxxz_xz[i] = 2.0 * tr_xxxyz_xz[i] * tbe_0 + 2.0 * tr_xxxz_xyz[i] * tke_0;

        tr_0_0_y_xxxz_yy[i] = 2.0 * tr_xxxyz_yy[i] * tbe_0 + 2.0 * tr_xxxz_yyy[i] * tke_0 - 2.0 * tr_xxxz_y[i];

        tr_0_0_y_xxxz_yz[i] = 2.0 * tr_xxxyz_yz[i] * tbe_0 + 2.0 * tr_xxxz_yyz[i] * tke_0 - tr_xxxz_z[i];

        tr_0_0_y_xxxz_zz[i] = 2.0 * tr_xxxyz_zz[i] * tbe_0 + 2.0 * tr_xxxz_yzz[i] * tke_0;
    }

    // Set up 108-114 components of targeted buffer : GD

    auto tr_0_0_y_xxyy_xx = pbuffer.data(idx_op_geom_010_gd + 108);

    auto tr_0_0_y_xxyy_xy = pbuffer.data(idx_op_geom_010_gd + 109);

    auto tr_0_0_y_xxyy_xz = pbuffer.data(idx_op_geom_010_gd + 110);

    auto tr_0_0_y_xxyy_yy = pbuffer.data(idx_op_geom_010_gd + 111);

    auto tr_0_0_y_xxyy_yz = pbuffer.data(idx_op_geom_010_gd + 112);

    auto tr_0_0_y_xxyy_zz = pbuffer.data(idx_op_geom_010_gd + 113);

    #pragma omp simd aligned(tr_0_0_y_xxyy_xx, tr_0_0_y_xxyy_xy, tr_0_0_y_xxyy_xz, tr_0_0_y_xxyy_yy, tr_0_0_y_xxyy_yz, tr_0_0_y_xxyy_zz, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyy_x, tr_xxyy_xxy, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_y, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_z, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyy_xx[i] = 2.0 * tr_xxyyy_xx[i] * tbe_0 + 2.0 * tr_xxyy_xxy[i] * tke_0 - 2.0 * tr_xxy_xx[i];

        tr_0_0_y_xxyy_xy[i] = 2.0 * tr_xxyyy_xy[i] * tbe_0 + 2.0 * tr_xxyy_xyy[i] * tke_0 - 2.0 * tr_xxy_xy[i] - tr_xxyy_x[i];

        tr_0_0_y_xxyy_xz[i] = 2.0 * tr_xxyyy_xz[i] * tbe_0 + 2.0 * tr_xxyy_xyz[i] * tke_0 - 2.0 * tr_xxy_xz[i];

        tr_0_0_y_xxyy_yy[i] = 2.0 * tr_xxyyy_yy[i] * tbe_0 + 2.0 * tr_xxyy_yyy[i] * tke_0 - 2.0 * tr_xxy_yy[i] - 2.0 * tr_xxyy_y[i];

        tr_0_0_y_xxyy_yz[i] = 2.0 * tr_xxyyy_yz[i] * tbe_0 + 2.0 * tr_xxyy_yyz[i] * tke_0 - 2.0 * tr_xxy_yz[i] - tr_xxyy_z[i];

        tr_0_0_y_xxyy_zz[i] = 2.0 * tr_xxyyy_zz[i] * tbe_0 + 2.0 * tr_xxyy_yzz[i] * tke_0 - 2.0 * tr_xxy_zz[i];
    }

    // Set up 114-120 components of targeted buffer : GD

    auto tr_0_0_y_xxyz_xx = pbuffer.data(idx_op_geom_010_gd + 114);

    auto tr_0_0_y_xxyz_xy = pbuffer.data(idx_op_geom_010_gd + 115);

    auto tr_0_0_y_xxyz_xz = pbuffer.data(idx_op_geom_010_gd + 116);

    auto tr_0_0_y_xxyz_yy = pbuffer.data(idx_op_geom_010_gd + 117);

    auto tr_0_0_y_xxyz_yz = pbuffer.data(idx_op_geom_010_gd + 118);

    auto tr_0_0_y_xxyz_zz = pbuffer.data(idx_op_geom_010_gd + 119);

    #pragma omp simd aligned(tr_0_0_y_xxyz_xx, tr_0_0_y_xxyz_xy, tr_0_0_y_xxyz_xz, tr_0_0_y_xxyz_yy, tr_0_0_y_xxyz_yz, tr_0_0_y_xxyz_zz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyz_x, tr_xxyz_xxy, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_y, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxyz_xx[i] = 2.0 * tr_xxyyz_xx[i] * tbe_0 + 2.0 * tr_xxyz_xxy[i] * tke_0 - tr_xxz_xx[i];

        tr_0_0_y_xxyz_xy[i] = 2.0 * tr_xxyyz_xy[i] * tbe_0 + 2.0 * tr_xxyz_xyy[i] * tke_0 - tr_xxz_xy[i] - tr_xxyz_x[i];

        tr_0_0_y_xxyz_xz[i] = 2.0 * tr_xxyyz_xz[i] * tbe_0 + 2.0 * tr_xxyz_xyz[i] * tke_0 - tr_xxz_xz[i];

        tr_0_0_y_xxyz_yy[i] = 2.0 * tr_xxyyz_yy[i] * tbe_0 + 2.0 * tr_xxyz_yyy[i] * tke_0 - tr_xxz_yy[i] - 2.0 * tr_xxyz_y[i];

        tr_0_0_y_xxyz_yz[i] = 2.0 * tr_xxyyz_yz[i] * tbe_0 + 2.0 * tr_xxyz_yyz[i] * tke_0 - tr_xxz_yz[i] - tr_xxyz_z[i];

        tr_0_0_y_xxyz_zz[i] = 2.0 * tr_xxyyz_zz[i] * tbe_0 + 2.0 * tr_xxyz_yzz[i] * tke_0 - tr_xxz_zz[i];
    }

    // Set up 120-126 components of targeted buffer : GD

    auto tr_0_0_y_xxzz_xx = pbuffer.data(idx_op_geom_010_gd + 120);

    auto tr_0_0_y_xxzz_xy = pbuffer.data(idx_op_geom_010_gd + 121);

    auto tr_0_0_y_xxzz_xz = pbuffer.data(idx_op_geom_010_gd + 122);

    auto tr_0_0_y_xxzz_yy = pbuffer.data(idx_op_geom_010_gd + 123);

    auto tr_0_0_y_xxzz_yz = pbuffer.data(idx_op_geom_010_gd + 124);

    auto tr_0_0_y_xxzz_zz = pbuffer.data(idx_op_geom_010_gd + 125);

    #pragma omp simd aligned(tr_0_0_y_xxzz_xx, tr_0_0_y_xxzz_xy, tr_0_0_y_xxzz_xz, tr_0_0_y_xxzz_yy, tr_0_0_y_xxzz_yz, tr_0_0_y_xxzz_zz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xxzz_x, tr_xxzz_xxy, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_y, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xxzz_xx[i] = 2.0 * tr_xxyzz_xx[i] * tbe_0 + 2.0 * tr_xxzz_xxy[i] * tke_0;

        tr_0_0_y_xxzz_xy[i] = 2.0 * tr_xxyzz_xy[i] * tbe_0 + 2.0 * tr_xxzz_xyy[i] * tke_0 - tr_xxzz_x[i];

        tr_0_0_y_xxzz_xz[i] = 2.0 * tr_xxyzz_xz[i] * tbe_0 + 2.0 * tr_xxzz_xyz[i] * tke_0;

        tr_0_0_y_xxzz_yy[i] = 2.0 * tr_xxyzz_yy[i] * tbe_0 + 2.0 * tr_xxzz_yyy[i] * tke_0 - 2.0 * tr_xxzz_y[i];

        tr_0_0_y_xxzz_yz[i] = 2.0 * tr_xxyzz_yz[i] * tbe_0 + 2.0 * tr_xxzz_yyz[i] * tke_0 - tr_xxzz_z[i];

        tr_0_0_y_xxzz_zz[i] = 2.0 * tr_xxyzz_zz[i] * tbe_0 + 2.0 * tr_xxzz_yzz[i] * tke_0;
    }

    // Set up 126-132 components of targeted buffer : GD

    auto tr_0_0_y_xyyy_xx = pbuffer.data(idx_op_geom_010_gd + 126);

    auto tr_0_0_y_xyyy_xy = pbuffer.data(idx_op_geom_010_gd + 127);

    auto tr_0_0_y_xyyy_xz = pbuffer.data(idx_op_geom_010_gd + 128);

    auto tr_0_0_y_xyyy_yy = pbuffer.data(idx_op_geom_010_gd + 129);

    auto tr_0_0_y_xyyy_yz = pbuffer.data(idx_op_geom_010_gd + 130);

    auto tr_0_0_y_xyyy_zz = pbuffer.data(idx_op_geom_010_gd + 131);

    #pragma omp simd aligned(tr_0_0_y_xyyy_xx, tr_0_0_y_xyyy_xy, tr_0_0_y_xyyy_xz, tr_0_0_y_xyyy_yy, tr_0_0_y_xyyy_yz, tr_0_0_y_xyyy_zz, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyy_x, tr_xyyy_xxy, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_y, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_z, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyy_xx[i] = 2.0 * tr_xyyyy_xx[i] * tbe_0 + 2.0 * tr_xyyy_xxy[i] * tke_0 - 3.0 * tr_xyy_xx[i];

        tr_0_0_y_xyyy_xy[i] = 2.0 * tr_xyyyy_xy[i] * tbe_0 + 2.0 * tr_xyyy_xyy[i] * tke_0 - 3.0 * tr_xyy_xy[i] - tr_xyyy_x[i];

        tr_0_0_y_xyyy_xz[i] = 2.0 * tr_xyyyy_xz[i] * tbe_0 + 2.0 * tr_xyyy_xyz[i] * tke_0 - 3.0 * tr_xyy_xz[i];

        tr_0_0_y_xyyy_yy[i] = 2.0 * tr_xyyyy_yy[i] * tbe_0 + 2.0 * tr_xyyy_yyy[i] * tke_0 - 3.0 * tr_xyy_yy[i] - 2.0 * tr_xyyy_y[i];

        tr_0_0_y_xyyy_yz[i] = 2.0 * tr_xyyyy_yz[i] * tbe_0 + 2.0 * tr_xyyy_yyz[i] * tke_0 - 3.0 * tr_xyy_yz[i] - tr_xyyy_z[i];

        tr_0_0_y_xyyy_zz[i] = 2.0 * tr_xyyyy_zz[i] * tbe_0 + 2.0 * tr_xyyy_yzz[i] * tke_0 - 3.0 * tr_xyy_zz[i];
    }

    // Set up 132-138 components of targeted buffer : GD

    auto tr_0_0_y_xyyz_xx = pbuffer.data(idx_op_geom_010_gd + 132);

    auto tr_0_0_y_xyyz_xy = pbuffer.data(idx_op_geom_010_gd + 133);

    auto tr_0_0_y_xyyz_xz = pbuffer.data(idx_op_geom_010_gd + 134);

    auto tr_0_0_y_xyyz_yy = pbuffer.data(idx_op_geom_010_gd + 135);

    auto tr_0_0_y_xyyz_yz = pbuffer.data(idx_op_geom_010_gd + 136);

    auto tr_0_0_y_xyyz_zz = pbuffer.data(idx_op_geom_010_gd + 137);

    #pragma omp simd aligned(tr_0_0_y_xyyz_xx, tr_0_0_y_xyyz_xy, tr_0_0_y_xyyz_xz, tr_0_0_y_xyyz_yy, tr_0_0_y_xyyz_yz, tr_0_0_y_xyyz_zz, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyz_x, tr_xyyz_xxy, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_y, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyyz_xx[i] = 2.0 * tr_xyyyz_xx[i] * tbe_0 + 2.0 * tr_xyyz_xxy[i] * tke_0 - 2.0 * tr_xyz_xx[i];

        tr_0_0_y_xyyz_xy[i] = 2.0 * tr_xyyyz_xy[i] * tbe_0 + 2.0 * tr_xyyz_xyy[i] * tke_0 - 2.0 * tr_xyz_xy[i] - tr_xyyz_x[i];

        tr_0_0_y_xyyz_xz[i] = 2.0 * tr_xyyyz_xz[i] * tbe_0 + 2.0 * tr_xyyz_xyz[i] * tke_0 - 2.0 * tr_xyz_xz[i];

        tr_0_0_y_xyyz_yy[i] = 2.0 * tr_xyyyz_yy[i] * tbe_0 + 2.0 * tr_xyyz_yyy[i] * tke_0 - 2.0 * tr_xyz_yy[i] - 2.0 * tr_xyyz_y[i];

        tr_0_0_y_xyyz_yz[i] = 2.0 * tr_xyyyz_yz[i] * tbe_0 + 2.0 * tr_xyyz_yyz[i] * tke_0 - 2.0 * tr_xyz_yz[i] - tr_xyyz_z[i];

        tr_0_0_y_xyyz_zz[i] = 2.0 * tr_xyyyz_zz[i] * tbe_0 + 2.0 * tr_xyyz_yzz[i] * tke_0 - 2.0 * tr_xyz_zz[i];
    }

    // Set up 138-144 components of targeted buffer : GD

    auto tr_0_0_y_xyzz_xx = pbuffer.data(idx_op_geom_010_gd + 138);

    auto tr_0_0_y_xyzz_xy = pbuffer.data(idx_op_geom_010_gd + 139);

    auto tr_0_0_y_xyzz_xz = pbuffer.data(idx_op_geom_010_gd + 140);

    auto tr_0_0_y_xyzz_yy = pbuffer.data(idx_op_geom_010_gd + 141);

    auto tr_0_0_y_xyzz_yz = pbuffer.data(idx_op_geom_010_gd + 142);

    auto tr_0_0_y_xyzz_zz = pbuffer.data(idx_op_geom_010_gd + 143);

    #pragma omp simd aligned(tr_0_0_y_xyzz_xx, tr_0_0_y_xyzz_xy, tr_0_0_y_xyzz_xz, tr_0_0_y_xyzz_yy, tr_0_0_y_xyzz_yz, tr_0_0_y_xyzz_zz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyzz_x, tr_xyzz_xxy, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_y, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xyzz_xx[i] = 2.0 * tr_xyyzz_xx[i] * tbe_0 + 2.0 * tr_xyzz_xxy[i] * tke_0 - tr_xzz_xx[i];

        tr_0_0_y_xyzz_xy[i] = 2.0 * tr_xyyzz_xy[i] * tbe_0 + 2.0 * tr_xyzz_xyy[i] * tke_0 - tr_xzz_xy[i] - tr_xyzz_x[i];

        tr_0_0_y_xyzz_xz[i] = 2.0 * tr_xyyzz_xz[i] * tbe_0 + 2.0 * tr_xyzz_xyz[i] * tke_0 - tr_xzz_xz[i];

        tr_0_0_y_xyzz_yy[i] = 2.0 * tr_xyyzz_yy[i] * tbe_0 + 2.0 * tr_xyzz_yyy[i] * tke_0 - tr_xzz_yy[i] - 2.0 * tr_xyzz_y[i];

        tr_0_0_y_xyzz_yz[i] = 2.0 * tr_xyyzz_yz[i] * tbe_0 + 2.0 * tr_xyzz_yyz[i] * tke_0 - tr_xzz_yz[i] - tr_xyzz_z[i];

        tr_0_0_y_xyzz_zz[i] = 2.0 * tr_xyyzz_zz[i] * tbe_0 + 2.0 * tr_xyzz_yzz[i] * tke_0 - tr_xzz_zz[i];
    }

    // Set up 144-150 components of targeted buffer : GD

    auto tr_0_0_y_xzzz_xx = pbuffer.data(idx_op_geom_010_gd + 144);

    auto tr_0_0_y_xzzz_xy = pbuffer.data(idx_op_geom_010_gd + 145);

    auto tr_0_0_y_xzzz_xz = pbuffer.data(idx_op_geom_010_gd + 146);

    auto tr_0_0_y_xzzz_yy = pbuffer.data(idx_op_geom_010_gd + 147);

    auto tr_0_0_y_xzzz_yz = pbuffer.data(idx_op_geom_010_gd + 148);

    auto tr_0_0_y_xzzz_zz = pbuffer.data(idx_op_geom_010_gd + 149);

    #pragma omp simd aligned(tr_0_0_y_xzzz_xx, tr_0_0_y_xzzz_xy, tr_0_0_y_xzzz_xz, tr_0_0_y_xzzz_yy, tr_0_0_y_xzzz_yz, tr_0_0_y_xzzz_zz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_xzzz_x, tr_xzzz_xxy, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_y, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_xzzz_xx[i] = 2.0 * tr_xyzzz_xx[i] * tbe_0 + 2.0 * tr_xzzz_xxy[i] * tke_0;

        tr_0_0_y_xzzz_xy[i] = 2.0 * tr_xyzzz_xy[i] * tbe_0 + 2.0 * tr_xzzz_xyy[i] * tke_0 - tr_xzzz_x[i];

        tr_0_0_y_xzzz_xz[i] = 2.0 * tr_xyzzz_xz[i] * tbe_0 + 2.0 * tr_xzzz_xyz[i] * tke_0;

        tr_0_0_y_xzzz_yy[i] = 2.0 * tr_xyzzz_yy[i] * tbe_0 + 2.0 * tr_xzzz_yyy[i] * tke_0 - 2.0 * tr_xzzz_y[i];

        tr_0_0_y_xzzz_yz[i] = 2.0 * tr_xyzzz_yz[i] * tbe_0 + 2.0 * tr_xzzz_yyz[i] * tke_0 - tr_xzzz_z[i];

        tr_0_0_y_xzzz_zz[i] = 2.0 * tr_xyzzz_zz[i] * tbe_0 + 2.0 * tr_xzzz_yzz[i] * tke_0;
    }

    // Set up 150-156 components of targeted buffer : GD

    auto tr_0_0_y_yyyy_xx = pbuffer.data(idx_op_geom_010_gd + 150);

    auto tr_0_0_y_yyyy_xy = pbuffer.data(idx_op_geom_010_gd + 151);

    auto tr_0_0_y_yyyy_xz = pbuffer.data(idx_op_geom_010_gd + 152);

    auto tr_0_0_y_yyyy_yy = pbuffer.data(idx_op_geom_010_gd + 153);

    auto tr_0_0_y_yyyy_yz = pbuffer.data(idx_op_geom_010_gd + 154);

    auto tr_0_0_y_yyyy_zz = pbuffer.data(idx_op_geom_010_gd + 155);

    #pragma omp simd aligned(tr_0_0_y_yyyy_xx, tr_0_0_y_yyyy_xy, tr_0_0_y_yyyy_xz, tr_0_0_y_yyyy_yy, tr_0_0_y_yyyy_yz, tr_0_0_y_yyyy_zz, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyy_x, tr_yyyy_xxy, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_y, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_z, tr_yyyyy_xx, tr_yyyyy_xy, tr_yyyyy_xz, tr_yyyyy_yy, tr_yyyyy_yz, tr_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyy_xx[i] = 2.0 * tr_yyyyy_xx[i] * tbe_0 + 2.0 * tr_yyyy_xxy[i] * tke_0 - 4.0 * tr_yyy_xx[i];

        tr_0_0_y_yyyy_xy[i] = 2.0 * tr_yyyyy_xy[i] * tbe_0 + 2.0 * tr_yyyy_xyy[i] * tke_0 - 4.0 * tr_yyy_xy[i] - tr_yyyy_x[i];

        tr_0_0_y_yyyy_xz[i] = 2.0 * tr_yyyyy_xz[i] * tbe_0 + 2.0 * tr_yyyy_xyz[i] * tke_0 - 4.0 * tr_yyy_xz[i];

        tr_0_0_y_yyyy_yy[i] = 2.0 * tr_yyyyy_yy[i] * tbe_0 + 2.0 * tr_yyyy_yyy[i] * tke_0 - 4.0 * tr_yyy_yy[i] - 2.0 * tr_yyyy_y[i];

        tr_0_0_y_yyyy_yz[i] = 2.0 * tr_yyyyy_yz[i] * tbe_0 + 2.0 * tr_yyyy_yyz[i] * tke_0 - 4.0 * tr_yyy_yz[i] - tr_yyyy_z[i];

        tr_0_0_y_yyyy_zz[i] = 2.0 * tr_yyyyy_zz[i] * tbe_0 + 2.0 * tr_yyyy_yzz[i] * tke_0 - 4.0 * tr_yyy_zz[i];
    }

    // Set up 156-162 components of targeted buffer : GD

    auto tr_0_0_y_yyyz_xx = pbuffer.data(idx_op_geom_010_gd + 156);

    auto tr_0_0_y_yyyz_xy = pbuffer.data(idx_op_geom_010_gd + 157);

    auto tr_0_0_y_yyyz_xz = pbuffer.data(idx_op_geom_010_gd + 158);

    auto tr_0_0_y_yyyz_yy = pbuffer.data(idx_op_geom_010_gd + 159);

    auto tr_0_0_y_yyyz_yz = pbuffer.data(idx_op_geom_010_gd + 160);

    auto tr_0_0_y_yyyz_zz = pbuffer.data(idx_op_geom_010_gd + 161);

    #pragma omp simd aligned(tr_0_0_y_yyyz_xx, tr_0_0_y_yyyz_xy, tr_0_0_y_yyyz_xz, tr_0_0_y_yyyz_yy, tr_0_0_y_yyyz_yz, tr_0_0_y_yyyz_zz, tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyyz_x, tr_yyyz_xxy, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_y, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_z, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyyz_xx[i] = 2.0 * tr_yyyyz_xx[i] * tbe_0 + 2.0 * tr_yyyz_xxy[i] * tke_0 - 3.0 * tr_yyz_xx[i];

        tr_0_0_y_yyyz_xy[i] = 2.0 * tr_yyyyz_xy[i] * tbe_0 + 2.0 * tr_yyyz_xyy[i] * tke_0 - 3.0 * tr_yyz_xy[i] - tr_yyyz_x[i];

        tr_0_0_y_yyyz_xz[i] = 2.0 * tr_yyyyz_xz[i] * tbe_0 + 2.0 * tr_yyyz_xyz[i] * tke_0 - 3.0 * tr_yyz_xz[i];

        tr_0_0_y_yyyz_yy[i] = 2.0 * tr_yyyyz_yy[i] * tbe_0 + 2.0 * tr_yyyz_yyy[i] * tke_0 - 3.0 * tr_yyz_yy[i] - 2.0 * tr_yyyz_y[i];

        tr_0_0_y_yyyz_yz[i] = 2.0 * tr_yyyyz_yz[i] * tbe_0 + 2.0 * tr_yyyz_yyz[i] * tke_0 - 3.0 * tr_yyz_yz[i] - tr_yyyz_z[i];

        tr_0_0_y_yyyz_zz[i] = 2.0 * tr_yyyyz_zz[i] * tbe_0 + 2.0 * tr_yyyz_yzz[i] * tke_0 - 3.0 * tr_yyz_zz[i];
    }

    // Set up 162-168 components of targeted buffer : GD

    auto tr_0_0_y_yyzz_xx = pbuffer.data(idx_op_geom_010_gd + 162);

    auto tr_0_0_y_yyzz_xy = pbuffer.data(idx_op_geom_010_gd + 163);

    auto tr_0_0_y_yyzz_xz = pbuffer.data(idx_op_geom_010_gd + 164);

    auto tr_0_0_y_yyzz_yy = pbuffer.data(idx_op_geom_010_gd + 165);

    auto tr_0_0_y_yyzz_yz = pbuffer.data(idx_op_geom_010_gd + 166);

    auto tr_0_0_y_yyzz_zz = pbuffer.data(idx_op_geom_010_gd + 167);

    #pragma omp simd aligned(tr_0_0_y_yyzz_xx, tr_0_0_y_yyzz_xy, tr_0_0_y_yyzz_xz, tr_0_0_y_yyzz_yy, tr_0_0_y_yyzz_yz, tr_0_0_y_yyzz_zz, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_yyzz_x, tr_yyzz_xxy, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_y, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_z, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yyzz_xx[i] = 2.0 * tr_yyyzz_xx[i] * tbe_0 + 2.0 * tr_yyzz_xxy[i] * tke_0 - 2.0 * tr_yzz_xx[i];

        tr_0_0_y_yyzz_xy[i] = 2.0 * tr_yyyzz_xy[i] * tbe_0 + 2.0 * tr_yyzz_xyy[i] * tke_0 - 2.0 * tr_yzz_xy[i] - tr_yyzz_x[i];

        tr_0_0_y_yyzz_xz[i] = 2.0 * tr_yyyzz_xz[i] * tbe_0 + 2.0 * tr_yyzz_xyz[i] * tke_0 - 2.0 * tr_yzz_xz[i];

        tr_0_0_y_yyzz_yy[i] = 2.0 * tr_yyyzz_yy[i] * tbe_0 + 2.0 * tr_yyzz_yyy[i] * tke_0 - 2.0 * tr_yzz_yy[i] - 2.0 * tr_yyzz_y[i];

        tr_0_0_y_yyzz_yz[i] = 2.0 * tr_yyyzz_yz[i] * tbe_0 + 2.0 * tr_yyzz_yyz[i] * tke_0 - 2.0 * tr_yzz_yz[i] - tr_yyzz_z[i];

        tr_0_0_y_yyzz_zz[i] = 2.0 * tr_yyyzz_zz[i] * tbe_0 + 2.0 * tr_yyzz_yzz[i] * tke_0 - 2.0 * tr_yzz_zz[i];
    }

    // Set up 168-174 components of targeted buffer : GD

    auto tr_0_0_y_yzzz_xx = pbuffer.data(idx_op_geom_010_gd + 168);

    auto tr_0_0_y_yzzz_xy = pbuffer.data(idx_op_geom_010_gd + 169);

    auto tr_0_0_y_yzzz_xz = pbuffer.data(idx_op_geom_010_gd + 170);

    auto tr_0_0_y_yzzz_yy = pbuffer.data(idx_op_geom_010_gd + 171);

    auto tr_0_0_y_yzzz_yz = pbuffer.data(idx_op_geom_010_gd + 172);

    auto tr_0_0_y_yzzz_zz = pbuffer.data(idx_op_geom_010_gd + 173);

    #pragma omp simd aligned(tr_0_0_y_yzzz_xx, tr_0_0_y_yzzz_xy, tr_0_0_y_yzzz_xz, tr_0_0_y_yzzz_yy, tr_0_0_y_yzzz_yz, tr_0_0_y_yzzz_zz, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_yzzz_x, tr_yzzz_xxy, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_y, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_z, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_yzzz_xx[i] = 2.0 * tr_yyzzz_xx[i] * tbe_0 + 2.0 * tr_yzzz_xxy[i] * tke_0 - tr_zzz_xx[i];

        tr_0_0_y_yzzz_xy[i] = 2.0 * tr_yyzzz_xy[i] * tbe_0 + 2.0 * tr_yzzz_xyy[i] * tke_0 - tr_zzz_xy[i] - tr_yzzz_x[i];

        tr_0_0_y_yzzz_xz[i] = 2.0 * tr_yyzzz_xz[i] * tbe_0 + 2.0 * tr_yzzz_xyz[i] * tke_0 - tr_zzz_xz[i];

        tr_0_0_y_yzzz_yy[i] = 2.0 * tr_yyzzz_yy[i] * tbe_0 + 2.0 * tr_yzzz_yyy[i] * tke_0 - tr_zzz_yy[i] - 2.0 * tr_yzzz_y[i];

        tr_0_0_y_yzzz_yz[i] = 2.0 * tr_yyzzz_yz[i] * tbe_0 + 2.0 * tr_yzzz_yyz[i] * tke_0 - tr_zzz_yz[i] - tr_yzzz_z[i];

        tr_0_0_y_yzzz_zz[i] = 2.0 * tr_yyzzz_zz[i] * tbe_0 + 2.0 * tr_yzzz_yzz[i] * tke_0 - tr_zzz_zz[i];
    }

    // Set up 174-180 components of targeted buffer : GD

    auto tr_0_0_y_zzzz_xx = pbuffer.data(idx_op_geom_010_gd + 174);

    auto tr_0_0_y_zzzz_xy = pbuffer.data(idx_op_geom_010_gd + 175);

    auto tr_0_0_y_zzzz_xz = pbuffer.data(idx_op_geom_010_gd + 176);

    auto tr_0_0_y_zzzz_yy = pbuffer.data(idx_op_geom_010_gd + 177);

    auto tr_0_0_y_zzzz_yz = pbuffer.data(idx_op_geom_010_gd + 178);

    auto tr_0_0_y_zzzz_zz = pbuffer.data(idx_op_geom_010_gd + 179);

    #pragma omp simd aligned(tr_0_0_y_zzzz_xx, tr_0_0_y_zzzz_xy, tr_0_0_y_zzzz_xz, tr_0_0_y_zzzz_yy, tr_0_0_y_zzzz_yz, tr_0_0_y_zzzz_zz, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, tr_zzzz_x, tr_zzzz_xxy, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_y, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_y_zzzz_xx[i] = 2.0 * tr_yzzzz_xx[i] * tbe_0 + 2.0 * tr_zzzz_xxy[i] * tke_0;

        tr_0_0_y_zzzz_xy[i] = 2.0 * tr_yzzzz_xy[i] * tbe_0 + 2.0 * tr_zzzz_xyy[i] * tke_0 - tr_zzzz_x[i];

        tr_0_0_y_zzzz_xz[i] = 2.0 * tr_yzzzz_xz[i] * tbe_0 + 2.0 * tr_zzzz_xyz[i] * tke_0;

        tr_0_0_y_zzzz_yy[i] = 2.0 * tr_yzzzz_yy[i] * tbe_0 + 2.0 * tr_zzzz_yyy[i] * tke_0 - 2.0 * tr_zzzz_y[i];

        tr_0_0_y_zzzz_yz[i] = 2.0 * tr_yzzzz_yz[i] * tbe_0 + 2.0 * tr_zzzz_yyz[i] * tke_0 - tr_zzzz_z[i];

        tr_0_0_y_zzzz_zz[i] = 2.0 * tr_yzzzz_zz[i] * tbe_0 + 2.0 * tr_zzzz_yzz[i] * tke_0;
    }

    // Set up 180-186 components of targeted buffer : GD

    auto tr_0_0_z_xxxx_xx = pbuffer.data(idx_op_geom_010_gd + 180);

    auto tr_0_0_z_xxxx_xy = pbuffer.data(idx_op_geom_010_gd + 181);

    auto tr_0_0_z_xxxx_xz = pbuffer.data(idx_op_geom_010_gd + 182);

    auto tr_0_0_z_xxxx_yy = pbuffer.data(idx_op_geom_010_gd + 183);

    auto tr_0_0_z_xxxx_yz = pbuffer.data(idx_op_geom_010_gd + 184);

    auto tr_0_0_z_xxxx_zz = pbuffer.data(idx_op_geom_010_gd + 185);

    #pragma omp simd aligned(tr_0_0_z_xxxx_xx, tr_0_0_z_xxxx_xy, tr_0_0_z_xxxx_xz, tr_0_0_z_xxxx_yy, tr_0_0_z_xxxx_yz, tr_0_0_z_xxxx_zz, tr_xxxx_x, tr_xxxx_xxz, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_y, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_z, tr_xxxx_zzz, tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxx_xx[i] = 2.0 * tr_xxxxz_xx[i] * tbe_0 + 2.0 * tr_xxxx_xxz[i] * tke_0;

        tr_0_0_z_xxxx_xy[i] = 2.0 * tr_xxxxz_xy[i] * tbe_0 + 2.0 * tr_xxxx_xyz[i] * tke_0;

        tr_0_0_z_xxxx_xz[i] = 2.0 * tr_xxxxz_xz[i] * tbe_0 + 2.0 * tr_xxxx_xzz[i] * tke_0 - tr_xxxx_x[i];

        tr_0_0_z_xxxx_yy[i] = 2.0 * tr_xxxxz_yy[i] * tbe_0 + 2.0 * tr_xxxx_yyz[i] * tke_0;

        tr_0_0_z_xxxx_yz[i] = 2.0 * tr_xxxxz_yz[i] * tbe_0 + 2.0 * tr_xxxx_yzz[i] * tke_0 - tr_xxxx_y[i];

        tr_0_0_z_xxxx_zz[i] = 2.0 * tr_xxxxz_zz[i] * tbe_0 + 2.0 * tr_xxxx_zzz[i] * tke_0 - 2.0 * tr_xxxx_z[i];
    }

    // Set up 186-192 components of targeted buffer : GD

    auto tr_0_0_z_xxxy_xx = pbuffer.data(idx_op_geom_010_gd + 186);

    auto tr_0_0_z_xxxy_xy = pbuffer.data(idx_op_geom_010_gd + 187);

    auto tr_0_0_z_xxxy_xz = pbuffer.data(idx_op_geom_010_gd + 188);

    auto tr_0_0_z_xxxy_yy = pbuffer.data(idx_op_geom_010_gd + 189);

    auto tr_0_0_z_xxxy_yz = pbuffer.data(idx_op_geom_010_gd + 190);

    auto tr_0_0_z_xxxy_zz = pbuffer.data(idx_op_geom_010_gd + 191);

    #pragma omp simd aligned(tr_0_0_z_xxxy_xx, tr_0_0_z_xxxy_xy, tr_0_0_z_xxxy_xz, tr_0_0_z_xxxy_yy, tr_0_0_z_xxxy_yz, tr_0_0_z_xxxy_zz, tr_xxxy_x, tr_xxxy_xxz, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_y, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_z, tr_xxxy_zzz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxy_xx[i] = 2.0 * tr_xxxyz_xx[i] * tbe_0 + 2.0 * tr_xxxy_xxz[i] * tke_0;

        tr_0_0_z_xxxy_xy[i] = 2.0 * tr_xxxyz_xy[i] * tbe_0 + 2.0 * tr_xxxy_xyz[i] * tke_0;

        tr_0_0_z_xxxy_xz[i] = 2.0 * tr_xxxyz_xz[i] * tbe_0 + 2.0 * tr_xxxy_xzz[i] * tke_0 - tr_xxxy_x[i];

        tr_0_0_z_xxxy_yy[i] = 2.0 * tr_xxxyz_yy[i] * tbe_0 + 2.0 * tr_xxxy_yyz[i] * tke_0;

        tr_0_0_z_xxxy_yz[i] = 2.0 * tr_xxxyz_yz[i] * tbe_0 + 2.0 * tr_xxxy_yzz[i] * tke_0 - tr_xxxy_y[i];

        tr_0_0_z_xxxy_zz[i] = 2.0 * tr_xxxyz_zz[i] * tbe_0 + 2.0 * tr_xxxy_zzz[i] * tke_0 - 2.0 * tr_xxxy_z[i];
    }

    // Set up 192-198 components of targeted buffer : GD

    auto tr_0_0_z_xxxz_xx = pbuffer.data(idx_op_geom_010_gd + 192);

    auto tr_0_0_z_xxxz_xy = pbuffer.data(idx_op_geom_010_gd + 193);

    auto tr_0_0_z_xxxz_xz = pbuffer.data(idx_op_geom_010_gd + 194);

    auto tr_0_0_z_xxxz_yy = pbuffer.data(idx_op_geom_010_gd + 195);

    auto tr_0_0_z_xxxz_yz = pbuffer.data(idx_op_geom_010_gd + 196);

    auto tr_0_0_z_xxxz_zz = pbuffer.data(idx_op_geom_010_gd + 197);

    #pragma omp simd aligned(tr_0_0_z_xxxz_xx, tr_0_0_z_xxxz_xy, tr_0_0_z_xxxz_xz, tr_0_0_z_xxxz_yy, tr_0_0_z_xxxz_yz, tr_0_0_z_xxxz_zz, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxz_x, tr_xxxz_xxz, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_y, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_z, tr_xxxz_zzz, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxxz_xx[i] = 2.0 * tr_xxxzz_xx[i] * tbe_0 + 2.0 * tr_xxxz_xxz[i] * tke_0 - tr_xxx_xx[i];

        tr_0_0_z_xxxz_xy[i] = 2.0 * tr_xxxzz_xy[i] * tbe_0 + 2.0 * tr_xxxz_xyz[i] * tke_0 - tr_xxx_xy[i];

        tr_0_0_z_xxxz_xz[i] = 2.0 * tr_xxxzz_xz[i] * tbe_0 + 2.0 * tr_xxxz_xzz[i] * tke_0 - tr_xxx_xz[i] - tr_xxxz_x[i];

        tr_0_0_z_xxxz_yy[i] = 2.0 * tr_xxxzz_yy[i] * tbe_0 + 2.0 * tr_xxxz_yyz[i] * tke_0 - tr_xxx_yy[i];

        tr_0_0_z_xxxz_yz[i] = 2.0 * tr_xxxzz_yz[i] * tbe_0 + 2.0 * tr_xxxz_yzz[i] * tke_0 - tr_xxx_yz[i] - tr_xxxz_y[i];

        tr_0_0_z_xxxz_zz[i] = 2.0 * tr_xxxzz_zz[i] * tbe_0 + 2.0 * tr_xxxz_zzz[i] * tke_0 - tr_xxx_zz[i] - 2.0 * tr_xxxz_z[i];
    }

    // Set up 198-204 components of targeted buffer : GD

    auto tr_0_0_z_xxyy_xx = pbuffer.data(idx_op_geom_010_gd + 198);

    auto tr_0_0_z_xxyy_xy = pbuffer.data(idx_op_geom_010_gd + 199);

    auto tr_0_0_z_xxyy_xz = pbuffer.data(idx_op_geom_010_gd + 200);

    auto tr_0_0_z_xxyy_yy = pbuffer.data(idx_op_geom_010_gd + 201);

    auto tr_0_0_z_xxyy_yz = pbuffer.data(idx_op_geom_010_gd + 202);

    auto tr_0_0_z_xxyy_zz = pbuffer.data(idx_op_geom_010_gd + 203);

    #pragma omp simd aligned(tr_0_0_z_xxyy_xx, tr_0_0_z_xxyy_xy, tr_0_0_z_xxyy_xz, tr_0_0_z_xxyy_yy, tr_0_0_z_xxyy_yz, tr_0_0_z_xxyy_zz, tr_xxyy_x, tr_xxyy_xxz, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_y, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_z, tr_xxyy_zzz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyy_xx[i] = 2.0 * tr_xxyyz_xx[i] * tbe_0 + 2.0 * tr_xxyy_xxz[i] * tke_0;

        tr_0_0_z_xxyy_xy[i] = 2.0 * tr_xxyyz_xy[i] * tbe_0 + 2.0 * tr_xxyy_xyz[i] * tke_0;

        tr_0_0_z_xxyy_xz[i] = 2.0 * tr_xxyyz_xz[i] * tbe_0 + 2.0 * tr_xxyy_xzz[i] * tke_0 - tr_xxyy_x[i];

        tr_0_0_z_xxyy_yy[i] = 2.0 * tr_xxyyz_yy[i] * tbe_0 + 2.0 * tr_xxyy_yyz[i] * tke_0;

        tr_0_0_z_xxyy_yz[i] = 2.0 * tr_xxyyz_yz[i] * tbe_0 + 2.0 * tr_xxyy_yzz[i] * tke_0 - tr_xxyy_y[i];

        tr_0_0_z_xxyy_zz[i] = 2.0 * tr_xxyyz_zz[i] * tbe_0 + 2.0 * tr_xxyy_zzz[i] * tke_0 - 2.0 * tr_xxyy_z[i];
    }

    // Set up 204-210 components of targeted buffer : GD

    auto tr_0_0_z_xxyz_xx = pbuffer.data(idx_op_geom_010_gd + 204);

    auto tr_0_0_z_xxyz_xy = pbuffer.data(idx_op_geom_010_gd + 205);

    auto tr_0_0_z_xxyz_xz = pbuffer.data(idx_op_geom_010_gd + 206);

    auto tr_0_0_z_xxyz_yy = pbuffer.data(idx_op_geom_010_gd + 207);

    auto tr_0_0_z_xxyz_yz = pbuffer.data(idx_op_geom_010_gd + 208);

    auto tr_0_0_z_xxyz_zz = pbuffer.data(idx_op_geom_010_gd + 209);

    #pragma omp simd aligned(tr_0_0_z_xxyz_xx, tr_0_0_z_xxyz_xy, tr_0_0_z_xxyz_xz, tr_0_0_z_xxyz_yy, tr_0_0_z_xxyz_yz, tr_0_0_z_xxyz_zz, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_xxz, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxyz_xx[i] = 2.0 * tr_xxyzz_xx[i] * tbe_0 + 2.0 * tr_xxyz_xxz[i] * tke_0 - tr_xxy_xx[i];

        tr_0_0_z_xxyz_xy[i] = 2.0 * tr_xxyzz_xy[i] * tbe_0 + 2.0 * tr_xxyz_xyz[i] * tke_0 - tr_xxy_xy[i];

        tr_0_0_z_xxyz_xz[i] = 2.0 * tr_xxyzz_xz[i] * tbe_0 + 2.0 * tr_xxyz_xzz[i] * tke_0 - tr_xxy_xz[i] - tr_xxyz_x[i];

        tr_0_0_z_xxyz_yy[i] = 2.0 * tr_xxyzz_yy[i] * tbe_0 + 2.0 * tr_xxyz_yyz[i] * tke_0 - tr_xxy_yy[i];

        tr_0_0_z_xxyz_yz[i] = 2.0 * tr_xxyzz_yz[i] * tbe_0 + 2.0 * tr_xxyz_yzz[i] * tke_0 - tr_xxy_yz[i] - tr_xxyz_y[i];

        tr_0_0_z_xxyz_zz[i] = 2.0 * tr_xxyzz_zz[i] * tbe_0 + 2.0 * tr_xxyz_zzz[i] * tke_0 - tr_xxy_zz[i] - 2.0 * tr_xxyz_z[i];
    }

    // Set up 210-216 components of targeted buffer : GD

    auto tr_0_0_z_xxzz_xx = pbuffer.data(idx_op_geom_010_gd + 210);

    auto tr_0_0_z_xxzz_xy = pbuffer.data(idx_op_geom_010_gd + 211);

    auto tr_0_0_z_xxzz_xz = pbuffer.data(idx_op_geom_010_gd + 212);

    auto tr_0_0_z_xxzz_yy = pbuffer.data(idx_op_geom_010_gd + 213);

    auto tr_0_0_z_xxzz_yz = pbuffer.data(idx_op_geom_010_gd + 214);

    auto tr_0_0_z_xxzz_zz = pbuffer.data(idx_op_geom_010_gd + 215);

    #pragma omp simd aligned(tr_0_0_z_xxzz_xx, tr_0_0_z_xxzz_xy, tr_0_0_z_xxzz_xz, tr_0_0_z_xxzz_yy, tr_0_0_z_xxzz_yz, tr_0_0_z_xxzz_zz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_xxz, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_y, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_z, tr_xxzz_zzz, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xxzz_xx[i] = 2.0 * tr_xxzzz_xx[i] * tbe_0 + 2.0 * tr_xxzz_xxz[i] * tke_0 - 2.0 * tr_xxz_xx[i];

        tr_0_0_z_xxzz_xy[i] = 2.0 * tr_xxzzz_xy[i] * tbe_0 + 2.0 * tr_xxzz_xyz[i] * tke_0 - 2.0 * tr_xxz_xy[i];

        tr_0_0_z_xxzz_xz[i] = 2.0 * tr_xxzzz_xz[i] * tbe_0 + 2.0 * tr_xxzz_xzz[i] * tke_0 - 2.0 * tr_xxz_xz[i] - tr_xxzz_x[i];

        tr_0_0_z_xxzz_yy[i] = 2.0 * tr_xxzzz_yy[i] * tbe_0 + 2.0 * tr_xxzz_yyz[i] * tke_0 - 2.0 * tr_xxz_yy[i];

        tr_0_0_z_xxzz_yz[i] = 2.0 * tr_xxzzz_yz[i] * tbe_0 + 2.0 * tr_xxzz_yzz[i] * tke_0 - 2.0 * tr_xxz_yz[i] - tr_xxzz_y[i];

        tr_0_0_z_xxzz_zz[i] = 2.0 * tr_xxzzz_zz[i] * tbe_0 + 2.0 * tr_xxzz_zzz[i] * tke_0 - 2.0 * tr_xxz_zz[i] - 2.0 * tr_xxzz_z[i];
    }

    // Set up 216-222 components of targeted buffer : GD

    auto tr_0_0_z_xyyy_xx = pbuffer.data(idx_op_geom_010_gd + 216);

    auto tr_0_0_z_xyyy_xy = pbuffer.data(idx_op_geom_010_gd + 217);

    auto tr_0_0_z_xyyy_xz = pbuffer.data(idx_op_geom_010_gd + 218);

    auto tr_0_0_z_xyyy_yy = pbuffer.data(idx_op_geom_010_gd + 219);

    auto tr_0_0_z_xyyy_yz = pbuffer.data(idx_op_geom_010_gd + 220);

    auto tr_0_0_z_xyyy_zz = pbuffer.data(idx_op_geom_010_gd + 221);

    #pragma omp simd aligned(tr_0_0_z_xyyy_xx, tr_0_0_z_xyyy_xy, tr_0_0_z_xyyy_xz, tr_0_0_z_xyyy_yy, tr_0_0_z_xyyy_yz, tr_0_0_z_xyyy_zz, tr_xyyy_x, tr_xyyy_xxz, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_y, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_z, tr_xyyy_zzz, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyy_xx[i] = 2.0 * tr_xyyyz_xx[i] * tbe_0 + 2.0 * tr_xyyy_xxz[i] * tke_0;

        tr_0_0_z_xyyy_xy[i] = 2.0 * tr_xyyyz_xy[i] * tbe_0 + 2.0 * tr_xyyy_xyz[i] * tke_0;

        tr_0_0_z_xyyy_xz[i] = 2.0 * tr_xyyyz_xz[i] * tbe_0 + 2.0 * tr_xyyy_xzz[i] * tke_0 - tr_xyyy_x[i];

        tr_0_0_z_xyyy_yy[i] = 2.0 * tr_xyyyz_yy[i] * tbe_0 + 2.0 * tr_xyyy_yyz[i] * tke_0;

        tr_0_0_z_xyyy_yz[i] = 2.0 * tr_xyyyz_yz[i] * tbe_0 + 2.0 * tr_xyyy_yzz[i] * tke_0 - tr_xyyy_y[i];

        tr_0_0_z_xyyy_zz[i] = 2.0 * tr_xyyyz_zz[i] * tbe_0 + 2.0 * tr_xyyy_zzz[i] * tke_0 - 2.0 * tr_xyyy_z[i];
    }

    // Set up 222-228 components of targeted buffer : GD

    auto tr_0_0_z_xyyz_xx = pbuffer.data(idx_op_geom_010_gd + 222);

    auto tr_0_0_z_xyyz_xy = pbuffer.data(idx_op_geom_010_gd + 223);

    auto tr_0_0_z_xyyz_xz = pbuffer.data(idx_op_geom_010_gd + 224);

    auto tr_0_0_z_xyyz_yy = pbuffer.data(idx_op_geom_010_gd + 225);

    auto tr_0_0_z_xyyz_yz = pbuffer.data(idx_op_geom_010_gd + 226);

    auto tr_0_0_z_xyyz_zz = pbuffer.data(idx_op_geom_010_gd + 227);

    #pragma omp simd aligned(tr_0_0_z_xyyz_xx, tr_0_0_z_xyyz_xy, tr_0_0_z_xyyz_xz, tr_0_0_z_xyyz_yy, tr_0_0_z_xyyz_yz, tr_0_0_z_xyyz_zz, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_xxz, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyyz_xx[i] = 2.0 * tr_xyyzz_xx[i] * tbe_0 + 2.0 * tr_xyyz_xxz[i] * tke_0 - tr_xyy_xx[i];

        tr_0_0_z_xyyz_xy[i] = 2.0 * tr_xyyzz_xy[i] * tbe_0 + 2.0 * tr_xyyz_xyz[i] * tke_0 - tr_xyy_xy[i];

        tr_0_0_z_xyyz_xz[i] = 2.0 * tr_xyyzz_xz[i] * tbe_0 + 2.0 * tr_xyyz_xzz[i] * tke_0 - tr_xyy_xz[i] - tr_xyyz_x[i];

        tr_0_0_z_xyyz_yy[i] = 2.0 * tr_xyyzz_yy[i] * tbe_0 + 2.0 * tr_xyyz_yyz[i] * tke_0 - tr_xyy_yy[i];

        tr_0_0_z_xyyz_yz[i] = 2.0 * tr_xyyzz_yz[i] * tbe_0 + 2.0 * tr_xyyz_yzz[i] * tke_0 - tr_xyy_yz[i] - tr_xyyz_y[i];

        tr_0_0_z_xyyz_zz[i] = 2.0 * tr_xyyzz_zz[i] * tbe_0 + 2.0 * tr_xyyz_zzz[i] * tke_0 - tr_xyy_zz[i] - 2.0 * tr_xyyz_z[i];
    }

    // Set up 228-234 components of targeted buffer : GD

    auto tr_0_0_z_xyzz_xx = pbuffer.data(idx_op_geom_010_gd + 228);

    auto tr_0_0_z_xyzz_xy = pbuffer.data(idx_op_geom_010_gd + 229);

    auto tr_0_0_z_xyzz_xz = pbuffer.data(idx_op_geom_010_gd + 230);

    auto tr_0_0_z_xyzz_yy = pbuffer.data(idx_op_geom_010_gd + 231);

    auto tr_0_0_z_xyzz_yz = pbuffer.data(idx_op_geom_010_gd + 232);

    auto tr_0_0_z_xyzz_zz = pbuffer.data(idx_op_geom_010_gd + 233);

    #pragma omp simd aligned(tr_0_0_z_xyzz_xx, tr_0_0_z_xyzz_xy, tr_0_0_z_xyzz_xz, tr_0_0_z_xyzz_yy, tr_0_0_z_xyzz_yz, tr_0_0_z_xyzz_zz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_xxz, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xyzz_xx[i] = 2.0 * tr_xyzzz_xx[i] * tbe_0 + 2.0 * tr_xyzz_xxz[i] * tke_0 - 2.0 * tr_xyz_xx[i];

        tr_0_0_z_xyzz_xy[i] = 2.0 * tr_xyzzz_xy[i] * tbe_0 + 2.0 * tr_xyzz_xyz[i] * tke_0 - 2.0 * tr_xyz_xy[i];

        tr_0_0_z_xyzz_xz[i] = 2.0 * tr_xyzzz_xz[i] * tbe_0 + 2.0 * tr_xyzz_xzz[i] * tke_0 - 2.0 * tr_xyz_xz[i] - tr_xyzz_x[i];

        tr_0_0_z_xyzz_yy[i] = 2.0 * tr_xyzzz_yy[i] * tbe_0 + 2.0 * tr_xyzz_yyz[i] * tke_0 - 2.0 * tr_xyz_yy[i];

        tr_0_0_z_xyzz_yz[i] = 2.0 * tr_xyzzz_yz[i] * tbe_0 + 2.0 * tr_xyzz_yzz[i] * tke_0 - 2.0 * tr_xyz_yz[i] - tr_xyzz_y[i];

        tr_0_0_z_xyzz_zz[i] = 2.0 * tr_xyzzz_zz[i] * tbe_0 + 2.0 * tr_xyzz_zzz[i] * tke_0 - 2.0 * tr_xyz_zz[i] - 2.0 * tr_xyzz_z[i];
    }

    // Set up 234-240 components of targeted buffer : GD

    auto tr_0_0_z_xzzz_xx = pbuffer.data(idx_op_geom_010_gd + 234);

    auto tr_0_0_z_xzzz_xy = pbuffer.data(idx_op_geom_010_gd + 235);

    auto tr_0_0_z_xzzz_xz = pbuffer.data(idx_op_geom_010_gd + 236);

    auto tr_0_0_z_xzzz_yy = pbuffer.data(idx_op_geom_010_gd + 237);

    auto tr_0_0_z_xzzz_yz = pbuffer.data(idx_op_geom_010_gd + 238);

    auto tr_0_0_z_xzzz_zz = pbuffer.data(idx_op_geom_010_gd + 239);

    #pragma omp simd aligned(tr_0_0_z_xzzz_xx, tr_0_0_z_xzzz_xy, tr_0_0_z_xzzz_xz, tr_0_0_z_xzzz_yy, tr_0_0_z_xzzz_yz, tr_0_0_z_xzzz_zz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_xxz, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_y, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_z, tr_xzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_xzzz_xx[i] = 2.0 * tr_xzzzz_xx[i] * tbe_0 + 2.0 * tr_xzzz_xxz[i] * tke_0 - 3.0 * tr_xzz_xx[i];

        tr_0_0_z_xzzz_xy[i] = 2.0 * tr_xzzzz_xy[i] * tbe_0 + 2.0 * tr_xzzz_xyz[i] * tke_0 - 3.0 * tr_xzz_xy[i];

        tr_0_0_z_xzzz_xz[i] = 2.0 * tr_xzzzz_xz[i] * tbe_0 + 2.0 * tr_xzzz_xzz[i] * tke_0 - 3.0 * tr_xzz_xz[i] - tr_xzzz_x[i];

        tr_0_0_z_xzzz_yy[i] = 2.0 * tr_xzzzz_yy[i] * tbe_0 + 2.0 * tr_xzzz_yyz[i] * tke_0 - 3.0 * tr_xzz_yy[i];

        tr_0_0_z_xzzz_yz[i] = 2.0 * tr_xzzzz_yz[i] * tbe_0 + 2.0 * tr_xzzz_yzz[i] * tke_0 - 3.0 * tr_xzz_yz[i] - tr_xzzz_y[i];

        tr_0_0_z_xzzz_zz[i] = 2.0 * tr_xzzzz_zz[i] * tbe_0 + 2.0 * tr_xzzz_zzz[i] * tke_0 - 3.0 * tr_xzz_zz[i] - 2.0 * tr_xzzz_z[i];
    }

    // Set up 240-246 components of targeted buffer : GD

    auto tr_0_0_z_yyyy_xx = pbuffer.data(idx_op_geom_010_gd + 240);

    auto tr_0_0_z_yyyy_xy = pbuffer.data(idx_op_geom_010_gd + 241);

    auto tr_0_0_z_yyyy_xz = pbuffer.data(idx_op_geom_010_gd + 242);

    auto tr_0_0_z_yyyy_yy = pbuffer.data(idx_op_geom_010_gd + 243);

    auto tr_0_0_z_yyyy_yz = pbuffer.data(idx_op_geom_010_gd + 244);

    auto tr_0_0_z_yyyy_zz = pbuffer.data(idx_op_geom_010_gd + 245);

    #pragma omp simd aligned(tr_0_0_z_yyyy_xx, tr_0_0_z_yyyy_xy, tr_0_0_z_yyyy_xz, tr_0_0_z_yyyy_yy, tr_0_0_z_yyyy_yz, tr_0_0_z_yyyy_zz, tr_yyyy_x, tr_yyyy_xxz, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_y, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_z, tr_yyyy_zzz, tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyy_xx[i] = 2.0 * tr_yyyyz_xx[i] * tbe_0 + 2.0 * tr_yyyy_xxz[i] * tke_0;

        tr_0_0_z_yyyy_xy[i] = 2.0 * tr_yyyyz_xy[i] * tbe_0 + 2.0 * tr_yyyy_xyz[i] * tke_0;

        tr_0_0_z_yyyy_xz[i] = 2.0 * tr_yyyyz_xz[i] * tbe_0 + 2.0 * tr_yyyy_xzz[i] * tke_0 - tr_yyyy_x[i];

        tr_0_0_z_yyyy_yy[i] = 2.0 * tr_yyyyz_yy[i] * tbe_0 + 2.0 * tr_yyyy_yyz[i] * tke_0;

        tr_0_0_z_yyyy_yz[i] = 2.0 * tr_yyyyz_yz[i] * tbe_0 + 2.0 * tr_yyyy_yzz[i] * tke_0 - tr_yyyy_y[i];

        tr_0_0_z_yyyy_zz[i] = 2.0 * tr_yyyyz_zz[i] * tbe_0 + 2.0 * tr_yyyy_zzz[i] * tke_0 - 2.0 * tr_yyyy_z[i];
    }

    // Set up 246-252 components of targeted buffer : GD

    auto tr_0_0_z_yyyz_xx = pbuffer.data(idx_op_geom_010_gd + 246);

    auto tr_0_0_z_yyyz_xy = pbuffer.data(idx_op_geom_010_gd + 247);

    auto tr_0_0_z_yyyz_xz = pbuffer.data(idx_op_geom_010_gd + 248);

    auto tr_0_0_z_yyyz_yy = pbuffer.data(idx_op_geom_010_gd + 249);

    auto tr_0_0_z_yyyz_yz = pbuffer.data(idx_op_geom_010_gd + 250);

    auto tr_0_0_z_yyyz_zz = pbuffer.data(idx_op_geom_010_gd + 251);

    #pragma omp simd aligned(tr_0_0_z_yyyz_xx, tr_0_0_z_yyyz_xy, tr_0_0_z_yyyz_xz, tr_0_0_z_yyyz_yy, tr_0_0_z_yyyz_yz, tr_0_0_z_yyyz_zz, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_xxz, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_y, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_z, tr_yyyz_zzz, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyyz_xx[i] = 2.0 * tr_yyyzz_xx[i] * tbe_0 + 2.0 * tr_yyyz_xxz[i] * tke_0 - tr_yyy_xx[i];

        tr_0_0_z_yyyz_xy[i] = 2.0 * tr_yyyzz_xy[i] * tbe_0 + 2.0 * tr_yyyz_xyz[i] * tke_0 - tr_yyy_xy[i];

        tr_0_0_z_yyyz_xz[i] = 2.0 * tr_yyyzz_xz[i] * tbe_0 + 2.0 * tr_yyyz_xzz[i] * tke_0 - tr_yyy_xz[i] - tr_yyyz_x[i];

        tr_0_0_z_yyyz_yy[i] = 2.0 * tr_yyyzz_yy[i] * tbe_0 + 2.0 * tr_yyyz_yyz[i] * tke_0 - tr_yyy_yy[i];

        tr_0_0_z_yyyz_yz[i] = 2.0 * tr_yyyzz_yz[i] * tbe_0 + 2.0 * tr_yyyz_yzz[i] * tke_0 - tr_yyy_yz[i] - tr_yyyz_y[i];

        tr_0_0_z_yyyz_zz[i] = 2.0 * tr_yyyzz_zz[i] * tbe_0 + 2.0 * tr_yyyz_zzz[i] * tke_0 - tr_yyy_zz[i] - 2.0 * tr_yyyz_z[i];
    }

    // Set up 252-258 components of targeted buffer : GD

    auto tr_0_0_z_yyzz_xx = pbuffer.data(idx_op_geom_010_gd + 252);

    auto tr_0_0_z_yyzz_xy = pbuffer.data(idx_op_geom_010_gd + 253);

    auto tr_0_0_z_yyzz_xz = pbuffer.data(idx_op_geom_010_gd + 254);

    auto tr_0_0_z_yyzz_yy = pbuffer.data(idx_op_geom_010_gd + 255);

    auto tr_0_0_z_yyzz_yz = pbuffer.data(idx_op_geom_010_gd + 256);

    auto tr_0_0_z_yyzz_zz = pbuffer.data(idx_op_geom_010_gd + 257);

    #pragma omp simd aligned(tr_0_0_z_yyzz_xx, tr_0_0_z_yyzz_xy, tr_0_0_z_yyzz_xz, tr_0_0_z_yyzz_yy, tr_0_0_z_yyzz_yz, tr_0_0_z_yyzz_zz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_xxz, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_y, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_z, tr_yyzz_zzz, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yyzz_xx[i] = 2.0 * tr_yyzzz_xx[i] * tbe_0 + 2.0 * tr_yyzz_xxz[i] * tke_0 - 2.0 * tr_yyz_xx[i];

        tr_0_0_z_yyzz_xy[i] = 2.0 * tr_yyzzz_xy[i] * tbe_0 + 2.0 * tr_yyzz_xyz[i] * tke_0 - 2.0 * tr_yyz_xy[i];

        tr_0_0_z_yyzz_xz[i] = 2.0 * tr_yyzzz_xz[i] * tbe_0 + 2.0 * tr_yyzz_xzz[i] * tke_0 - 2.0 * tr_yyz_xz[i] - tr_yyzz_x[i];

        tr_0_0_z_yyzz_yy[i] = 2.0 * tr_yyzzz_yy[i] * tbe_0 + 2.0 * tr_yyzz_yyz[i] * tke_0 - 2.0 * tr_yyz_yy[i];

        tr_0_0_z_yyzz_yz[i] = 2.0 * tr_yyzzz_yz[i] * tbe_0 + 2.0 * tr_yyzz_yzz[i] * tke_0 - 2.0 * tr_yyz_yz[i] - tr_yyzz_y[i];

        tr_0_0_z_yyzz_zz[i] = 2.0 * tr_yyzzz_zz[i] * tbe_0 + 2.0 * tr_yyzz_zzz[i] * tke_0 - 2.0 * tr_yyz_zz[i] - 2.0 * tr_yyzz_z[i];
    }

    // Set up 258-264 components of targeted buffer : GD

    auto tr_0_0_z_yzzz_xx = pbuffer.data(idx_op_geom_010_gd + 258);

    auto tr_0_0_z_yzzz_xy = pbuffer.data(idx_op_geom_010_gd + 259);

    auto tr_0_0_z_yzzz_xz = pbuffer.data(idx_op_geom_010_gd + 260);

    auto tr_0_0_z_yzzz_yy = pbuffer.data(idx_op_geom_010_gd + 261);

    auto tr_0_0_z_yzzz_yz = pbuffer.data(idx_op_geom_010_gd + 262);

    auto tr_0_0_z_yzzz_zz = pbuffer.data(idx_op_geom_010_gd + 263);

    #pragma omp simd aligned(tr_0_0_z_yzzz_xx, tr_0_0_z_yzzz_xy, tr_0_0_z_yzzz_xz, tr_0_0_z_yzzz_yy, tr_0_0_z_yzzz_yz, tr_0_0_z_yzzz_zz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_xxz, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_y, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_z, tr_yzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_yzzz_xx[i] = 2.0 * tr_yzzzz_xx[i] * tbe_0 + 2.0 * tr_yzzz_xxz[i] * tke_0 - 3.0 * tr_yzz_xx[i];

        tr_0_0_z_yzzz_xy[i] = 2.0 * tr_yzzzz_xy[i] * tbe_0 + 2.0 * tr_yzzz_xyz[i] * tke_0 - 3.0 * tr_yzz_xy[i];

        tr_0_0_z_yzzz_xz[i] = 2.0 * tr_yzzzz_xz[i] * tbe_0 + 2.0 * tr_yzzz_xzz[i] * tke_0 - 3.0 * tr_yzz_xz[i] - tr_yzzz_x[i];

        tr_0_0_z_yzzz_yy[i] = 2.0 * tr_yzzzz_yy[i] * tbe_0 + 2.0 * tr_yzzz_yyz[i] * tke_0 - 3.0 * tr_yzz_yy[i];

        tr_0_0_z_yzzz_yz[i] = 2.0 * tr_yzzzz_yz[i] * tbe_0 + 2.0 * tr_yzzz_yzz[i] * tke_0 - 3.0 * tr_yzz_yz[i] - tr_yzzz_y[i];

        tr_0_0_z_yzzz_zz[i] = 2.0 * tr_yzzzz_zz[i] * tbe_0 + 2.0 * tr_yzzz_zzz[i] * tke_0 - 3.0 * tr_yzz_zz[i] - 2.0 * tr_yzzz_z[i];
    }

    // Set up 264-270 components of targeted buffer : GD

    auto tr_0_0_z_zzzz_xx = pbuffer.data(idx_op_geom_010_gd + 264);

    auto tr_0_0_z_zzzz_xy = pbuffer.data(idx_op_geom_010_gd + 265);

    auto tr_0_0_z_zzzz_xz = pbuffer.data(idx_op_geom_010_gd + 266);

    auto tr_0_0_z_zzzz_yy = pbuffer.data(idx_op_geom_010_gd + 267);

    auto tr_0_0_z_zzzz_yz = pbuffer.data(idx_op_geom_010_gd + 268);

    auto tr_0_0_z_zzzz_zz = pbuffer.data(idx_op_geom_010_gd + 269);

    #pragma omp simd aligned(tr_0_0_z_zzzz_xx, tr_0_0_z_zzzz_xy, tr_0_0_z_zzzz_xz, tr_0_0_z_zzzz_yy, tr_0_0_z_zzzz_yz, tr_0_0_z_zzzz_zz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_xxz, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_y, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_z, tr_zzzz_zzz, tr_zzzzz_xx, tr_zzzzz_xy, tr_zzzzz_xz, tr_zzzzz_yy, tr_zzzzz_yz, tr_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_0_0_z_zzzz_xx[i] = 2.0 * tr_zzzzz_xx[i] * tbe_0 + 2.0 * tr_zzzz_xxz[i] * tke_0 - 4.0 * tr_zzz_xx[i];

        tr_0_0_z_zzzz_xy[i] = 2.0 * tr_zzzzz_xy[i] * tbe_0 + 2.0 * tr_zzzz_xyz[i] * tke_0 - 4.0 * tr_zzz_xy[i];

        tr_0_0_z_zzzz_xz[i] = 2.0 * tr_zzzzz_xz[i] * tbe_0 + 2.0 * tr_zzzz_xzz[i] * tke_0 - 4.0 * tr_zzz_xz[i] - tr_zzzz_x[i];

        tr_0_0_z_zzzz_yy[i] = 2.0 * tr_zzzzz_yy[i] * tbe_0 + 2.0 * tr_zzzz_yyz[i] * tke_0 - 4.0 * tr_zzz_yy[i];

        tr_0_0_z_zzzz_yz[i] = 2.0 * tr_zzzzz_yz[i] * tbe_0 + 2.0 * tr_zzzz_yzz[i] * tke_0 - 4.0 * tr_zzz_yz[i] - tr_zzzz_y[i];

        tr_0_0_z_zzzz_zz[i] = 2.0 * tr_zzzzz_zz[i] * tbe_0 + 2.0 * tr_zzzz_zzz[i] * tke_0 - 4.0 * tr_zzz_zz[i] - 2.0 * tr_zzzz_z[i];
    }

}

} // t2cgeom namespace

