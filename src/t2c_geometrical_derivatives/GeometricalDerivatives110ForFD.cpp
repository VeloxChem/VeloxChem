#include "GeometricalDerivatives110ForFD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_110_fd(CSimdArray<double>& pbuffer,
                         const int idx_op_geom_110_fd,
                         const int idx_op_pd,
                         const int idx_op_dp,
                         const int idx_op_df,
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

    // Set up components of auxiliary buffer : DP

    auto tr_xx_x = pbuffer.data(idx_op_dp);

    auto tr_xx_y = pbuffer.data(idx_op_dp + 1);

    auto tr_xx_z = pbuffer.data(idx_op_dp + 2);

    auto tr_xy_x = pbuffer.data(idx_op_dp + 3);

    auto tr_xy_y = pbuffer.data(idx_op_dp + 4);

    auto tr_xy_z = pbuffer.data(idx_op_dp + 5);

    auto tr_xz_x = pbuffer.data(idx_op_dp + 6);

    auto tr_xz_y = pbuffer.data(idx_op_dp + 7);

    auto tr_xz_z = pbuffer.data(idx_op_dp + 8);

    auto tr_yy_x = pbuffer.data(idx_op_dp + 9);

    auto tr_yy_y = pbuffer.data(idx_op_dp + 10);

    auto tr_yy_z = pbuffer.data(idx_op_dp + 11);

    auto tr_yz_x = pbuffer.data(idx_op_dp + 12);

    auto tr_yz_y = pbuffer.data(idx_op_dp + 13);

    auto tr_yz_z = pbuffer.data(idx_op_dp + 14);

    auto tr_zz_x = pbuffer.data(idx_op_dp + 15);

    auto tr_zz_y = pbuffer.data(idx_op_dp + 16);

    auto tr_zz_z = pbuffer.data(idx_op_dp + 17);

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

    // Set up 0-6 components of targeted buffer : FD

    auto tr_x_0_x_xxx_xx = pbuffer.data(idx_op_geom_110_fd);

    auto tr_x_0_x_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 1);

    auto tr_x_0_x_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 2);

    auto tr_x_0_x_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 3);

    auto tr_x_0_x_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 4);

    auto tr_x_0_x_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 5);

    #pragma omp simd aligned(tr_x_0_x_xxx_xx, tr_x_0_x_xxx_xy, tr_x_0_x_xxx_xz, tr_x_0_x_xxx_yy, tr_x_0_x_xxx_yz, tr_x_0_x_xxx_zz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_z, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxx_x, tr_xxxx_xxx, tr_xxxx_xxy, tr_xxxx_xxz, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_y, tr_xxxx_z, tr_xxxxx_xx, tr_xxxxx_xy, tr_xxxxx_xz, tr_xxxxx_yy, tr_xxxxx_yz, tr_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxx_xx[i] = 6.0 * tr_x_xx[i] + 6.0 * tr_xx_x[i] - 6.0 * tr_xx_xxx[i] * tke_0 - 14.0 * tr_xxx_xx[i] * tbe_0 - 4.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxx_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xy[i] = 6.0 * tr_x_xy[i] + 3.0 * tr_xx_y[i] - 6.0 * tr_xx_xxy[i] * tke_0 - 14.0 * tr_xxx_xy[i] * tbe_0 - 2.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxx_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_xz[i] = 6.0 * tr_x_xz[i] + 3.0 * tr_xx_z[i] - 6.0 * tr_xx_xxz[i] * tke_0 - 14.0 * tr_xxx_xz[i] * tbe_0 - 2.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxx_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_yy[i] = 6.0 * tr_x_yy[i] - 6.0 * tr_xx_xyy[i] * tke_0 - 14.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxxx_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_yz[i] = 6.0 * tr_x_yz[i] - 6.0 * tr_xx_xyz[i] * tke_0 - 14.0 * tr_xxx_yz[i] * tbe_0 + 4.0 * tr_xxxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxx_zz[i] = 6.0 * tr_x_zz[i] - 6.0 * tr_xx_xzz[i] * tke_0 - 14.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxxx_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxx_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 6-12 components of targeted buffer : FD

    auto tr_x_0_x_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 6);

    auto tr_x_0_x_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 7);

    auto tr_x_0_x_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 8);

    auto tr_x_0_x_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 9);

    auto tr_x_0_x_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 10);

    auto tr_x_0_x_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 11);

    #pragma omp simd aligned(tr_x_0_x_xxy_xx, tr_x_0_x_xxy_xy, tr_x_0_x_xxy_xz, tr_x_0_x_xxy_yy, tr_x_0_x_xxy_yz, tr_x_0_x_xxy_zz, tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxy_x, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_y, tr_xxxy_z, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xy_x, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_z, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxy_xx[i] = 2.0 * tr_y_xx[i] + 4.0 * tr_xy_x[i] - 4.0 * tr_xy_xxx[i] * tke_0 - 10.0 * tr_xxy_xx[i] * tbe_0 - 4.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xy[i] = 2.0 * tr_y_xy[i] + 2.0 * tr_xy_y[i] - 4.0 * tr_xy_xxy[i] * tke_0 - 10.0 * tr_xxy_xy[i] * tbe_0 - 2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_xz[i] = 2.0 * tr_y_xz[i] + 2.0 * tr_xy_z[i] - 4.0 * tr_xy_xxz[i] * tke_0 - 10.0 * tr_xxy_xz[i] * tbe_0 - 2.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_yy[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_xy_xyy[i] * tke_0 - 10.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_yz[i] = 2.0 * tr_y_yz[i] - 4.0 * tr_xy_xyz[i] * tke_0 - 10.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxy_zz[i] = 2.0 * tr_y_zz[i] - 4.0 * tr_xy_xzz[i] * tke_0 - 10.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 12-18 components of targeted buffer : FD

    auto tr_x_0_x_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 12);

    auto tr_x_0_x_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 13);

    auto tr_x_0_x_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 14);

    auto tr_x_0_x_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 15);

    auto tr_x_0_x_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 16);

    auto tr_x_0_x_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 17);

    #pragma omp simd aligned(tr_x_0_x_xxz_xx, tr_x_0_x_xxz_xy, tr_x_0_x_xxz_xz, tr_x_0_x_xxz_yy, tr_x_0_x_xxz_yz, tr_x_0_x_xxz_zz, tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxz_x, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_y, tr_xxxz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xz_x, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xxz_xx[i] = 2.0 * tr_z_xx[i] + 4.0 * tr_xz_x[i] - 4.0 * tr_xz_xxx[i] * tke_0 - 10.0 * tr_xxz_xx[i] * tbe_0 - 4.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xy[i] = 2.0 * tr_z_xy[i] + 2.0 * tr_xz_y[i] - 4.0 * tr_xz_xxy[i] * tke_0 - 10.0 * tr_xxz_xy[i] * tbe_0 - 2.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_xz[i] = 2.0 * tr_z_xz[i] + 2.0 * tr_xz_z[i] - 4.0 * tr_xz_xxz[i] * tke_0 - 10.0 * tr_xxz_xz[i] * tbe_0 - 2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_yy[i] = 2.0 * tr_z_yy[i] - 4.0 * tr_xz_xyy[i] * tke_0 - 10.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_yz[i] = 2.0 * tr_z_yz[i] - 4.0 * tr_xz_xyz[i] * tke_0 - 10.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xxz_zz[i] = 2.0 * tr_z_zz[i] - 4.0 * tr_xz_xzz[i] * tke_0 - 10.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 18-24 components of targeted buffer : FD

    auto tr_x_0_x_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 18);

    auto tr_x_0_x_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 19);

    auto tr_x_0_x_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 20);

    auto tr_x_0_x_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 21);

    auto tr_x_0_x_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 22);

    auto tr_x_0_x_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 23);

    #pragma omp simd aligned(tr_x_0_x_xyy_xx, tr_x_0_x_xyy_xy, tr_x_0_x_xyy_xz, tr_x_0_x_xyy_yy, tr_x_0_x_xyy_yz, tr_x_0_x_xyy_zz, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxyy_x, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_y, tr_xxyy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyy_xx[i] = 2.0 * tr_yy_x[i] - 2.0 * tr_yy_xxx[i] * tke_0 - 6.0 * tr_xyy_xx[i] * tbe_0 - 4.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xy[i] = tr_yy_y[i] - 2.0 * tr_yy_xxy[i] * tke_0 - 6.0 * tr_xyy_xy[i] * tbe_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_xz[i] = tr_yy_z[i] - 2.0 * tr_yy_xxz[i] * tke_0 - 6.0 * tr_xyy_xz[i] * tbe_0 - 2.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_yy[i] = -2.0 * tr_yy_xyy[i] * tke_0 - 6.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_yz[i] = -2.0 * tr_yy_xyz[i] * tke_0 - 6.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyy_zz[i] = -2.0 * tr_yy_xzz[i] * tke_0 - 6.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 24-30 components of targeted buffer : FD

    auto tr_x_0_x_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 24);

    auto tr_x_0_x_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 25);

    auto tr_x_0_x_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 26);

    auto tr_x_0_x_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 27);

    auto tr_x_0_x_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 28);

    auto tr_x_0_x_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 29);

    #pragma omp simd aligned(tr_x_0_x_xyz_xx, tr_x_0_x_xyz_xy, tr_x_0_x_xyz_xz, tr_x_0_x_xyz_yy, tr_x_0_x_xyz_yz, tr_x_0_x_xyz_zz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_x, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xyz_xx[i] = 2.0 * tr_yz_x[i] - 2.0 * tr_yz_xxx[i] * tke_0 - 6.0 * tr_xyz_xx[i] * tbe_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xy[i] = tr_yz_y[i] - 2.0 * tr_yz_xxy[i] * tke_0 - 6.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_xz[i] = tr_yz_z[i] - 2.0 * tr_yz_xxz[i] * tke_0 - 6.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_yy[i] = -2.0 * tr_yz_xyy[i] * tke_0 - 6.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_yz[i] = -2.0 * tr_yz_xyz[i] * tke_0 - 6.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xyz_zz[i] = -2.0 * tr_yz_xzz[i] * tke_0 - 6.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 30-36 components of targeted buffer : FD

    auto tr_x_0_x_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 30);

    auto tr_x_0_x_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 31);

    auto tr_x_0_x_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 32);

    auto tr_x_0_x_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 33);

    auto tr_x_0_x_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 34);

    auto tr_x_0_x_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 35);

    #pragma omp simd aligned(tr_x_0_x_xzz_xx, tr_x_0_x_xzz_xy, tr_x_0_x_xzz_xz, tr_x_0_x_xzz_yy, tr_x_0_x_xzz_yz, tr_x_0_x_xzz_zz, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxzz_x, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_y, tr_xxzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_xzz_xx[i] = 2.0 * tr_zz_x[i] - 2.0 * tr_zz_xxx[i] * tke_0 - 6.0 * tr_xzz_xx[i] * tbe_0 - 4.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xy[i] = tr_zz_y[i] - 2.0 * tr_zz_xxy[i] * tke_0 - 6.0 * tr_xzz_xy[i] * tbe_0 - 2.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_xz[i] = tr_zz_z[i] - 2.0 * tr_zz_xxz[i] * tke_0 - 6.0 * tr_xzz_xz[i] * tbe_0 - 2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_yy[i] = -2.0 * tr_zz_xyy[i] * tke_0 - 6.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_yz[i] = -2.0 * tr_zz_xyz[i] * tke_0 - 6.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_xzz_zz[i] = -2.0 * tr_zz_xzz[i] * tke_0 - 6.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 36-42 components of targeted buffer : FD

    auto tr_x_0_x_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 36);

    auto tr_x_0_x_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 37);

    auto tr_x_0_x_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 38);

    auto tr_x_0_x_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 39);

    auto tr_x_0_x_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 40);

    auto tr_x_0_x_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 41);

    #pragma omp simd aligned(tr_x_0_x_yyy_xx, tr_x_0_x_yyy_xy, tr_x_0_x_yyy_xz, tr_x_0_x_yyy_yy, tr_x_0_x_yyy_yz, tr_x_0_x_yyy_zz, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_xyyy_x, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_y, tr_xyyy_z, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyy_xx[i] = -2.0 * tr_yyy_xx[i] * tbe_0 - 4.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xy[i] = -2.0 * tr_yyy_xy[i] * tbe_0 - 2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_xz[i] = -2.0 * tr_yyy_xz[i] * tbe_0 - 2.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_yy[i] = -2.0 * tr_yyy_yy[i] * tbe_0 + 4.0 * tr_xyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_yz[i] = -2.0 * tr_yyy_yz[i] * tbe_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyy_zz[i] = -2.0 * tr_yyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 42-48 components of targeted buffer : FD

    auto tr_x_0_x_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 42);

    auto tr_x_0_x_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 43);

    auto tr_x_0_x_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 44);

    auto tr_x_0_x_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 45);

    auto tr_x_0_x_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 46);

    auto tr_x_0_x_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 47);

    #pragma omp simd aligned(tr_x_0_x_yyz_xx, tr_x_0_x_yyz_xy, tr_x_0_x_yyz_xz, tr_x_0_x_yyz_yy, tr_x_0_x_yyz_yz, tr_x_0_x_yyz_zz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xyyz_x, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_z, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yyz_xx[i] = -2.0 * tr_yyz_xx[i] * tbe_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xy[i] = -2.0 * tr_yyz_xy[i] * tbe_0 - 2.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_xz[i] = -2.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_yy[i] = -2.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_yz[i] = -2.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yyz_zz[i] = -2.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 48-54 components of targeted buffer : FD

    auto tr_x_0_x_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 48);

    auto tr_x_0_x_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 49);

    auto tr_x_0_x_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 50);

    auto tr_x_0_x_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 51);

    auto tr_x_0_x_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 52);

    auto tr_x_0_x_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 53);

    #pragma omp simd aligned(tr_x_0_x_yzz_xx, tr_x_0_x_yzz_xy, tr_x_0_x_yzz_xz, tr_x_0_x_yzz_yy, tr_x_0_x_yzz_yz, tr_x_0_x_yzz_zz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xyzz_x, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_z, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_yzz_xx[i] = -2.0 * tr_yzz_xx[i] * tbe_0 - 4.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xy[i] = -2.0 * tr_yzz_xy[i] * tbe_0 - 2.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_xz[i] = -2.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_yy[i] = -2.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_yz[i] = -2.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_yzz_zz[i] = -2.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 54-60 components of targeted buffer : FD

    auto tr_x_0_x_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 54);

    auto tr_x_0_x_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 55);

    auto tr_x_0_x_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 56);

    auto tr_x_0_x_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 57);

    auto tr_x_0_x_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 58);

    auto tr_x_0_x_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 59);

    #pragma omp simd aligned(tr_x_0_x_zzz_xx, tr_x_0_x_zzz_xy, tr_x_0_x_zzz_xz, tr_x_0_x_zzz_yy, tr_x_0_x_zzz_yz, tr_x_0_x_zzz_zz, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_xzzz_x, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_y, tr_xzzz_z, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_x_zzz_xx[i] = -2.0 * tr_zzz_xx[i] * tbe_0 - 4.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xy[i] = -2.0 * tr_zzz_xy[i] * tbe_0 - 2.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_xz[i] = -2.0 * tr_zzz_xz[i] * tbe_0 - 2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_yy[i] = -2.0 * tr_zzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_yz[i] = -2.0 * tr_zzz_yz[i] * tbe_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_x_zzz_zz[i] = -2.0 * tr_zzz_zz[i] * tbe_0 + 4.0 * tr_xzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 60-66 components of targeted buffer : FD

    auto tr_x_0_y_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 60);

    auto tr_x_0_y_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 61);

    auto tr_x_0_y_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 62);

    auto tr_x_0_y_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 63);

    auto tr_x_0_y_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 64);

    auto tr_x_0_y_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 65);

    #pragma omp simd aligned(tr_x_0_y_xxx_xx, tr_x_0_y_xxx_xy, tr_x_0_y_xxx_xz, tr_x_0_y_xxx_yy, tr_x_0_y_xxx_yz, tr_x_0_y_xxx_zz, tr_xx_x, tr_xx_xxy, tr_xx_xyy, tr_xx_xyz, tr_xx_y, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xxxx_x, tr_xxxx_xxy, tr_xxxx_xyy, tr_xxxx_xyz, tr_xxxx_y, tr_xxxx_yyy, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_z, tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxx_xx[i] = -6.0 * tr_xx_xxy[i] * tke_0 - 6.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxxx_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xy[i] = 3.0 * tr_xx_x[i] - 6.0 * tr_xx_xyy[i] * tke_0 - 6.0 * tr_xxy_xy[i] * tbe_0 - 2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxx_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_xz[i] = -6.0 * tr_xx_xyz[i] * tke_0 - 6.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_yy[i] = 6.0 * tr_xx_y[i] - 6.0 * tr_xx_yyy[i] * tke_0 - 6.0 * tr_xxy_yy[i] * tbe_0 - 4.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxx_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_yz[i] = 3.0 * tr_xx_z[i] - 6.0 * tr_xx_yyz[i] * tke_0 - 6.0 * tr_xxy_yz[i] * tbe_0 - 2.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxx_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxx_zz[i] = -6.0 * tr_xx_yzz[i] * tke_0 - 6.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxxx_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 66-72 components of targeted buffer : FD

    auto tr_x_0_y_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 66);

    auto tr_x_0_y_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 67);

    auto tr_x_0_y_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 68);

    auto tr_x_0_y_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 69);

    auto tr_x_0_y_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 70);

    auto tr_x_0_y_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 71);

    #pragma omp simd aligned(tr_x_0_y_xxy_xx, tr_x_0_y_xxy_xy, tr_x_0_y_xxy_xz, tr_x_0_y_xxy_yy, tr_x_0_y_xxy_yz, tr_x_0_y_xxy_zz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxy_x, tr_xxxy_xxy, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_y, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_z, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_xy_x, tr_xy_xxy, tr_xy_xyy, tr_xy_xyz, tr_xy_y, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxy_xx[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_xy_xxy[i] * tke_0 - 4.0 * tr_xyy_xx[i] * tbe_0 - 2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xy[i] = 2.0 * tr_x_xy[i] + 2.0 * tr_xy_x[i] - 4.0 * tr_xy_xyy[i] * tke_0 - 4.0 * tr_xyy_xy[i] * tbe_0 - 2.0 * tr_xxx_xy[i] * tbe_0 - 2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_xz[i] = 2.0 * tr_x_xz[i] - 4.0 * tr_xy_xyz[i] * tke_0 - 4.0 * tr_xyy_xz[i] * tbe_0 - 2.0 * tr_xxx_xz[i] * tbe_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_yy[i] = 2.0 * tr_x_yy[i] + 4.0 * tr_xy_y[i] - 4.0 * tr_xy_yyy[i] * tke_0 - 4.0 * tr_xyy_yy[i] * tbe_0 - 2.0 * tr_xxx_yy[i] * tbe_0 - 4.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_yz[i] = 2.0 * tr_x_yz[i] + 2.0 * tr_xy_z[i] - 4.0 * tr_xy_yyz[i] * tke_0 - 4.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xxx_yz[i] * tbe_0 - 2.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxy_zz[i] = 2.0 * tr_x_zz[i] - 4.0 * tr_xy_yzz[i] * tke_0 - 4.0 * tr_xyy_zz[i] * tbe_0 - 2.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 72-78 components of targeted buffer : FD

    auto tr_x_0_y_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 72);

    auto tr_x_0_y_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 73);

    auto tr_x_0_y_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 74);

    auto tr_x_0_y_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 75);

    auto tr_x_0_y_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 76);

    auto tr_x_0_y_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 77);

    #pragma omp simd aligned(tr_x_0_y_xxz_xx, tr_x_0_y_xxz_xy, tr_x_0_y_xxz_xz, tr_x_0_y_xxz_yy, tr_x_0_y_xxz_yz, tr_x_0_y_xxz_zz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxxz_x, tr_xxxz_xxy, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_y, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xz_x, tr_xz_xxy, tr_xz_xyy, tr_xz_xyz, tr_xz_y, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xxz_xx[i] = -4.0 * tr_xz_xxy[i] * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xy[i] = 2.0 * tr_xz_x[i] - 4.0 * tr_xz_xyy[i] * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_xz[i] = -4.0 * tr_xz_xyz[i] * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_yy[i] = 4.0 * tr_xz_y[i] - 4.0 * tr_xz_yyy[i] * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 - 4.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_yz[i] = 2.0 * tr_xz_z[i] - 4.0 * tr_xz_yyz[i] * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xxz_zz[i] = -4.0 * tr_xz_yzz[i] * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 78-84 components of targeted buffer : FD

    auto tr_x_0_y_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 78);

    auto tr_x_0_y_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 79);

    auto tr_x_0_y_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 80);

    auto tr_x_0_y_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 81);

    auto tr_x_0_y_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 82);

    auto tr_x_0_y_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 83);

    #pragma omp simd aligned(tr_x_0_y_xyy_xx, tr_x_0_y_xyy_xy, tr_x_0_y_xyy_xz, tr_x_0_y_xyy_yy, tr_x_0_y_xyy_yz, tr_x_0_y_xyy_zz, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyy_x, tr_xxyy_xxy, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_y, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_z, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_xxy, tr_yy_xyy, tr_yy_xyz, tr_yy_y, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyy_xx[i] = 2.0 * tr_y_xx[i] - 2.0 * tr_yy_xxy[i] * tke_0 - 2.0 * tr_yyy_xx[i] * tbe_0 - 4.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xy[i] = 2.0 * tr_y_xy[i] + tr_yy_x[i] - 2.0 * tr_yy_xyy[i] * tke_0 - 2.0 * tr_yyy_xy[i] * tbe_0 - 4.0 * tr_xxy_xy[i] * tbe_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_xz[i] = 2.0 * tr_y_xz[i] - 2.0 * tr_yy_xyz[i] * tke_0 - 2.0 * tr_yyy_xz[i] * tbe_0 - 4.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_yy[i] = 2.0 * tr_y_yy[i] + 2.0 * tr_yy_y[i] - 2.0 * tr_yy_yyy[i] * tke_0 - 2.0 * tr_yyy_yy[i] * tbe_0 - 4.0 * tr_xxy_yy[i] * tbe_0 - 4.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_yz[i] = 2.0 * tr_y_yz[i] + tr_yy_z[i] - 2.0 * tr_yy_yyz[i] * tke_0 - 2.0 * tr_yyy_yz[i] * tbe_0 - 4.0 * tr_xxy_yz[i] * tbe_0 - 2.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyy_zz[i] = 2.0 * tr_y_zz[i] - 2.0 * tr_yy_yzz[i] * tke_0 - 2.0 * tr_yyy_zz[i] * tbe_0 - 4.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 84-90 components of targeted buffer : FD

    auto tr_x_0_y_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 84);

    auto tr_x_0_y_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 85);

    auto tr_x_0_y_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 86);

    auto tr_x_0_y_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 87);

    auto tr_x_0_y_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 88);

    auto tr_x_0_y_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 89);

    #pragma omp simd aligned(tr_x_0_y_xyz_xx, tr_x_0_y_xyz_xy, tr_x_0_y_xyz_xz, tr_x_0_y_xyz_yy, tr_x_0_y_xyz_yz, tr_x_0_y_xyz_zz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyz_x, tr_xxyz_xxy, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_y, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yz_x, tr_yz_xxy, tr_yz_xyy, tr_yz_xyz, tr_yz_y, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xyz_xx[i] = tr_z_xx[i] - 2.0 * tr_yz_xxy[i] * tke_0 - 2.0 * tr_yyz_xx[i] * tbe_0 - 2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xy[i] = tr_z_xy[i] + tr_yz_x[i] - 2.0 * tr_yz_xyy[i] * tke_0 - 2.0 * tr_yyz_xy[i] * tbe_0 - 2.0 * tr_xxz_xy[i] * tbe_0 - 2.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_xz[i] = tr_z_xz[i] - 2.0 * tr_yz_xyz[i] * tke_0 - 2.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_yy[i] = tr_z_yy[i] + 2.0 * tr_yz_y[i] - 2.0 * tr_yz_yyy[i] * tke_0 - 2.0 * tr_yyz_yy[i] * tbe_0 - 2.0 * tr_xxz_yy[i] * tbe_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_yz[i] = tr_z_yz[i] + tr_yz_z[i] - 2.0 * tr_yz_yyz[i] * tke_0 - 2.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xyz_zz[i] = tr_z_zz[i] - 2.0 * tr_yz_yzz[i] * tke_0 - 2.0 * tr_yyz_zz[i] * tbe_0 - 2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 90-96 components of targeted buffer : FD

    auto tr_x_0_y_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 90);

    auto tr_x_0_y_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 91);

    auto tr_x_0_y_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 92);

    auto tr_x_0_y_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 93);

    auto tr_x_0_y_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 94);

    auto tr_x_0_y_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 95);

    #pragma omp simd aligned(tr_x_0_y_xzz_xx, tr_x_0_y_xzz_xy, tr_x_0_y_xzz_xz, tr_x_0_y_xzz_yy, tr_x_0_y_xzz_yz, tr_x_0_y_xzz_zz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xxzz_x, tr_xxzz_xxy, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_y, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_z, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_zz_x, tr_zz_xxy, tr_zz_xyy, tr_zz_xyz, tr_zz_y, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_xzz_xx[i] = -2.0 * tr_zz_xxy[i] * tke_0 - 2.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xy[i] = tr_zz_x[i] - 2.0 * tr_zz_xyy[i] * tke_0 - 2.0 * tr_yzz_xy[i] * tbe_0 - 2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_xz[i] = -2.0 * tr_zz_xyz[i] * tke_0 - 2.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_yy[i] = 2.0 * tr_zz_y[i] - 2.0 * tr_zz_yyy[i] * tke_0 - 2.0 * tr_yzz_yy[i] * tbe_0 - 4.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_yz[i] = tr_zz_z[i] - 2.0 * tr_zz_yyz[i] * tke_0 - 2.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_xzz_zz[i] = -2.0 * tr_zz_yzz[i] * tke_0 - 2.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 96-102 components of targeted buffer : FD

    auto tr_x_0_y_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 96);

    auto tr_x_0_y_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 97);

    auto tr_x_0_y_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 98);

    auto tr_x_0_y_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 99);

    auto tr_x_0_y_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 100);

    auto tr_x_0_y_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 101);

    #pragma omp simd aligned(tr_x_0_y_yyy_xx, tr_x_0_y_yyy_xy, tr_x_0_y_yyy_xz, tr_x_0_y_yyy_yy, tr_x_0_y_yyy_yz, tr_x_0_y_yyy_zz, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyy_x, tr_xyyy_xxy, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_y, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_z, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyy_xx[i] = -6.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xy[i] = -6.0 * tr_xyy_xy[i] * tbe_0 - 2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_xz[i] = -6.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_yy[i] = -6.0 * tr_xyy_yy[i] * tbe_0 - 4.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_yz[i] = -6.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyy_zz[i] = -6.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 102-108 components of targeted buffer : FD

    auto tr_x_0_y_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 102);

    auto tr_x_0_y_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 103);

    auto tr_x_0_y_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 104);

    auto tr_x_0_y_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 105);

    auto tr_x_0_y_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 106);

    auto tr_x_0_y_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 107);

    #pragma omp simd aligned(tr_x_0_y_yyz_xx, tr_x_0_y_yyz_xy, tr_x_0_y_yyz_xz, tr_x_0_y_yyz_yy, tr_x_0_y_yyz_yz, tr_x_0_y_yyz_zz, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyz_x, tr_xyyz_xxy, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_y, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yyz_xx[i] = -4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xy[i] = -4.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_xz[i] = -4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_yy[i] = -4.0 * tr_xyz_yy[i] * tbe_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_yz[i] = -4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yyz_zz[i] = -4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 108-114 components of targeted buffer : FD

    auto tr_x_0_y_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 108);

    auto tr_x_0_y_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 109);

    auto tr_x_0_y_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 110);

    auto tr_x_0_y_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 111);

    auto tr_x_0_y_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 112);

    auto tr_x_0_y_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 113);

    #pragma omp simd aligned(tr_x_0_y_yzz_xx, tr_x_0_y_yzz_xy, tr_x_0_y_yzz_xz, tr_x_0_y_yzz_yy, tr_x_0_y_yzz_yz, tr_x_0_y_yzz_zz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyzz_x, tr_xyzz_xxy, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_y, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_yzz_xx[i] = -2.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xy[i] = -2.0 * tr_xzz_xy[i] * tbe_0 - 2.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_xz[i] = -2.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_yy[i] = -2.0 * tr_xzz_yy[i] * tbe_0 - 4.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_yz[i] = -2.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_yzz_zz[i] = -2.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 114-120 components of targeted buffer : FD

    auto tr_x_0_y_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 114);

    auto tr_x_0_y_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 115);

    auto tr_x_0_y_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 116);

    auto tr_x_0_y_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 117);

    auto tr_x_0_y_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 118);

    auto tr_x_0_y_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 119);

    #pragma omp simd aligned(tr_x_0_y_zzz_xx, tr_x_0_y_zzz_xy, tr_x_0_y_zzz_xz, tr_x_0_y_zzz_yy, tr_x_0_y_zzz_yz, tr_x_0_y_zzz_zz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_xzzz_x, tr_xzzz_xxy, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_y, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_y_zzz_xx[i] = 4.0 * tr_xzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xy[i] = -2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_xz[i] = 4.0 * tr_xzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_yy[i] = -4.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_yz[i] = -2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_y_zzz_zz[i] = 4.0 * tr_xzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 120-126 components of targeted buffer : FD

    auto tr_x_0_z_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 120);

    auto tr_x_0_z_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 121);

    auto tr_x_0_z_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 122);

    auto tr_x_0_z_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 123);

    auto tr_x_0_z_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 124);

    auto tr_x_0_z_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 125);

    #pragma omp simd aligned(tr_x_0_z_xxx_xx, tr_x_0_z_xxx_xy, tr_x_0_z_xxx_xz, tr_x_0_z_xxx_yy, tr_x_0_z_xxx_yz, tr_x_0_z_xxx_zz, tr_xx_x, tr_xx_xxz, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxxx_x, tr_xxxx_xxz, tr_xxxx_xyz, tr_xxxx_xzz, tr_xxxx_y, tr_xxxx_yyz, tr_xxxx_yzz, tr_xxxx_z, tr_xxxx_zzz, tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxx_xx[i] = -6.0 * tr_xx_xxz[i] * tke_0 - 6.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxxx_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xy[i] = -6.0 * tr_xx_xyz[i] * tke_0 - 6.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxxx_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_xz[i] = 3.0 * tr_xx_x[i] - 6.0 * tr_xx_xzz[i] * tke_0 - 6.0 * tr_xxz_xz[i] * tbe_0 - 2.0 * tr_xxxx_x[i] * tbe_0 + 4.0 * tr_xxxx_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_yy[i] = -6.0 * tr_xx_yyz[i] * tke_0 - 6.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxxx_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_yz[i] = 3.0 * tr_xx_y[i] - 6.0 * tr_xx_yzz[i] * tke_0 - 6.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxxx_y[i] * tbe_0 + 4.0 * tr_xxxx_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxx_zz[i] = 6.0 * tr_xx_z[i] - 6.0 * tr_xx_zzz[i] * tke_0 - 6.0 * tr_xxz_zz[i] * tbe_0 - 4.0 * tr_xxxx_z[i] * tbe_0 + 4.0 * tr_xxxx_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 126-132 components of targeted buffer : FD

    auto tr_x_0_z_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 126);

    auto tr_x_0_z_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 127);

    auto tr_x_0_z_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 128);

    auto tr_x_0_z_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 129);

    auto tr_x_0_z_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 130);

    auto tr_x_0_z_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 131);

    #pragma omp simd aligned(tr_x_0_z_xxy_xx, tr_x_0_z_xxy_xy, tr_x_0_z_xxy_xz, tr_x_0_z_xxy_yy, tr_x_0_z_xxy_yz, tr_x_0_z_xxy_zz, tr_xxxy_x, tr_xxxy_xxz, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_y, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_z, tr_xxxy_zzz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xy_x, tr_xy_xxz, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxy_xx[i] = -4.0 * tr_xy_xxz[i] * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xy[i] = -4.0 * tr_xy_xyz[i] * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_xz[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xy_xzz[i] * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_yy[i] = -4.0 * tr_xy_yyz[i] * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_yz[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xy_yzz[i] * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxy_zz[i] = 4.0 * tr_xy_z[i] - 4.0 * tr_xy_zzz[i] * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 - 4.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 132-138 components of targeted buffer : FD

    auto tr_x_0_z_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 132);

    auto tr_x_0_z_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 133);

    auto tr_x_0_z_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 134);

    auto tr_x_0_z_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 135);

    auto tr_x_0_z_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 136);

    auto tr_x_0_z_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 137);

    #pragma omp simd aligned(tr_x_0_z_xxz_xx, tr_x_0_z_xxz_xy, tr_x_0_z_xxz_xz, tr_x_0_z_xxz_yy, tr_x_0_z_xxz_yz, tr_x_0_z_xxz_zz, tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxz_x, tr_xxxz_xxz, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_y, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_z, tr_xxxz_zzz, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_xz_x, tr_xz_xxz, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xxz_xx[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_xz_xxz[i] * tke_0 - 4.0 * tr_xzz_xx[i] * tbe_0 - 2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xy[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_xz_xyz[i] * tke_0 - 4.0 * tr_xzz_xy[i] * tbe_0 - 2.0 * tr_xxx_xy[i] * tbe_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_xz[i] = 2.0 * tr_x_xz[i] + 2.0 * tr_xz_x[i] - 4.0 * tr_xz_xzz[i] * tke_0 - 4.0 * tr_xzz_xz[i] * tbe_0 - 2.0 * tr_xxx_xz[i] * tbe_0 - 2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_yy[i] = 2.0 * tr_x_yy[i] - 4.0 * tr_xz_yyz[i] * tke_0 - 4.0 * tr_xzz_yy[i] * tbe_0 - 2.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_yz[i] = 2.0 * tr_x_yz[i] + 2.0 * tr_xz_y[i] - 4.0 * tr_xz_yzz[i] * tke_0 - 4.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xxx_yz[i] * tbe_0 - 2.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xxz_zz[i] = 2.0 * tr_x_zz[i] + 4.0 * tr_xz_z[i] - 4.0 * tr_xz_zzz[i] * tke_0 - 4.0 * tr_xzz_zz[i] * tbe_0 - 2.0 * tr_xxx_zz[i] * tbe_0 - 4.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 138-144 components of targeted buffer : FD

    auto tr_x_0_z_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 138);

    auto tr_x_0_z_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 139);

    auto tr_x_0_z_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 140);

    auto tr_x_0_z_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 141);

    auto tr_x_0_z_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 142);

    auto tr_x_0_z_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 143);

    #pragma omp simd aligned(tr_x_0_z_xyy_xx, tr_x_0_z_xyy_xy, tr_x_0_z_xyy_xz, tr_x_0_z_xyy_yy, tr_x_0_z_xyy_yz, tr_x_0_z_xyy_zz, tr_xxyy_x, tr_xxyy_xxz, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_y, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_z, tr_xxyy_zzz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_yy_x, tr_yy_xxz, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyy_xx[i] = -2.0 * tr_yy_xxz[i] * tke_0 - 2.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xy[i] = -2.0 * tr_yy_xyz[i] * tke_0 - 2.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_xz[i] = tr_yy_x[i] - 2.0 * tr_yy_xzz[i] * tke_0 - 2.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_yy[i] = -2.0 * tr_yy_yyz[i] * tke_0 - 2.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_yz[i] = tr_yy_y[i] - 2.0 * tr_yy_yzz[i] * tke_0 - 2.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyy_zz[i] = 2.0 * tr_yy_z[i] - 2.0 * tr_yy_zzz[i] * tke_0 - 2.0 * tr_yyz_zz[i] * tbe_0 - 4.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 144-150 components of targeted buffer : FD

    auto tr_x_0_z_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 144);

    auto tr_x_0_z_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 145);

    auto tr_x_0_z_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 146);

    auto tr_x_0_z_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 147);

    auto tr_x_0_z_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 148);

    auto tr_x_0_z_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 149);

    #pragma omp simd aligned(tr_x_0_z_xyz_xx, tr_x_0_z_xyz_xy, tr_x_0_z_xyz_xz, tr_x_0_z_xyz_yy, tr_x_0_z_xyz_yz, tr_x_0_z_xyz_zz, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_xxz, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_xxz, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xyz_xx[i] = tr_y_xx[i] - 2.0 * tr_yz_xxz[i] * tke_0 - 2.0 * tr_yzz_xx[i] * tbe_0 - 2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xy[i] = tr_y_xy[i] - 2.0 * tr_yz_xyz[i] * tke_0 - 2.0 * tr_yzz_xy[i] * tbe_0 - 2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_xz[i] = tr_y_xz[i] + tr_yz_x[i] - 2.0 * tr_yz_xzz[i] * tke_0 - 2.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_xxy_xz[i] * tbe_0 - 2.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_yy[i] = tr_y_yy[i] - 2.0 * tr_yz_yyz[i] * tke_0 - 2.0 * tr_yzz_yy[i] * tbe_0 - 2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_yz[i] = tr_y_yz[i] + tr_yz_y[i] - 2.0 * tr_yz_yzz[i] * tke_0 - 2.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_xxy_yz[i] * tbe_0 - 2.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xyz_zz[i] = tr_y_zz[i] + 2.0 * tr_yz_z[i] - 2.0 * tr_yz_zzz[i] * tke_0 - 2.0 * tr_yzz_zz[i] * tbe_0 - 2.0 * tr_xxy_zz[i] * tbe_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 150-156 components of targeted buffer : FD

    auto tr_x_0_z_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 150);

    auto tr_x_0_z_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 151);

    auto tr_x_0_z_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 152);

    auto tr_x_0_z_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 153);

    auto tr_x_0_z_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 154);

    auto tr_x_0_z_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 155);

    #pragma omp simd aligned(tr_x_0_z_xzz_xx, tr_x_0_z_xzz_xy, tr_x_0_z_xzz_xz, tr_x_0_z_xzz_yy, tr_x_0_z_xzz_yz, tr_x_0_z_xzz_zz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_xxz, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_y, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_z, tr_xxzz_zzz, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_xxz, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_xzz_xx[i] = 2.0 * tr_z_xx[i] - 2.0 * tr_zz_xxz[i] * tke_0 - 2.0 * tr_zzz_xx[i] * tbe_0 - 4.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xy[i] = 2.0 * tr_z_xy[i] - 2.0 * tr_zz_xyz[i] * tke_0 - 2.0 * tr_zzz_xy[i] * tbe_0 - 4.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_xz[i] = 2.0 * tr_z_xz[i] + tr_zz_x[i] - 2.0 * tr_zz_xzz[i] * tke_0 - 2.0 * tr_zzz_xz[i] * tbe_0 - 4.0 * tr_xxz_xz[i] * tbe_0 - 2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_yy[i] = 2.0 * tr_z_yy[i] - 2.0 * tr_zz_yyz[i] * tke_0 - 2.0 * tr_zzz_yy[i] * tbe_0 - 4.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_yz[i] = 2.0 * tr_z_yz[i] + tr_zz_y[i] - 2.0 * tr_zz_yzz[i] * tke_0 - 2.0 * tr_zzz_yz[i] * tbe_0 - 4.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_xzz_zz[i] = 2.0 * tr_z_zz[i] + 2.0 * tr_zz_z[i] - 2.0 * tr_zz_zzz[i] * tke_0 - 2.0 * tr_zzz_zz[i] * tbe_0 - 4.0 * tr_xxz_zz[i] * tbe_0 - 4.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 156-162 components of targeted buffer : FD

    auto tr_x_0_z_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 156);

    auto tr_x_0_z_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 157);

    auto tr_x_0_z_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 158);

    auto tr_x_0_z_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 159);

    auto tr_x_0_z_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 160);

    auto tr_x_0_z_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 161);

    #pragma omp simd aligned(tr_x_0_z_yyy_xx, tr_x_0_z_yyy_xy, tr_x_0_z_yyy_xz, tr_x_0_z_yyy_yy, tr_x_0_z_yyy_yz, tr_x_0_z_yyy_zz, tr_xyyy_x, tr_xyyy_xxz, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_y, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_z, tr_xyyy_zzz, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyy_xx[i] = 4.0 * tr_xyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xy[i] = 4.0 * tr_xyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_xz[i] = -2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_yy[i] = 4.0 * tr_xyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_yz[i] = -2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyy_zz[i] = -4.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 162-168 components of targeted buffer : FD

    auto tr_x_0_z_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 162);

    auto tr_x_0_z_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 163);

    auto tr_x_0_z_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 164);

    auto tr_x_0_z_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 165);

    auto tr_x_0_z_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 166);

    auto tr_x_0_z_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 167);

    #pragma omp simd aligned(tr_x_0_z_yyz_xx, tr_x_0_z_yyz_xy, tr_x_0_z_yyz_xz, tr_x_0_z_yyz_yy, tr_x_0_z_yyz_yz, tr_x_0_z_yyz_zz, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_xxz, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yyz_xx[i] = -2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xy[i] = -2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_xz[i] = -2.0 * tr_xyy_xz[i] * tbe_0 - 2.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_yy[i] = -2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_yz[i] = -2.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yyz_zz[i] = -2.0 * tr_xyy_zz[i] * tbe_0 - 4.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 168-174 components of targeted buffer : FD

    auto tr_x_0_z_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 168);

    auto tr_x_0_z_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 169);

    auto tr_x_0_z_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 170);

    auto tr_x_0_z_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 171);

    auto tr_x_0_z_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 172);

    auto tr_x_0_z_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 173);

    #pragma omp simd aligned(tr_x_0_z_yzz_xx, tr_x_0_z_yzz_xy, tr_x_0_z_yzz_xz, tr_x_0_z_yzz_yy, tr_x_0_z_yzz_yz, tr_x_0_z_yzz_zz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_xxz, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_yzz_xx[i] = -4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xy[i] = -4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_xz[i] = -4.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_yy[i] = -4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_yz[i] = -4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_yzz_zz[i] = -4.0 * tr_xyz_zz[i] * tbe_0 - 4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 174-180 components of targeted buffer : FD

    auto tr_x_0_z_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 174);

    auto tr_x_0_z_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 175);

    auto tr_x_0_z_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 176);

    auto tr_x_0_z_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 177);

    auto tr_x_0_z_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 178);

    auto tr_x_0_z_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 179);

    #pragma omp simd aligned(tr_x_0_z_zzz_xx, tr_x_0_z_zzz_xy, tr_x_0_z_zzz_xz, tr_x_0_z_zzz_yy, tr_x_0_z_zzz_yz, tr_x_0_z_zzz_zz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_xxz, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_y, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_z, tr_xzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_x_0_z_zzz_xx[i] = -6.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xx[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xy[i] = -6.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_xz[i] = -6.0 * tr_xzz_xz[i] * tbe_0 - 2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_yy[i] = -6.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yy[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_yz[i] = -6.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yz[i] * tbe_0 * tbe_0;

        tr_x_0_z_zzz_zz[i] = -6.0 * tr_xzz_zz[i] * tbe_0 - 4.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 180-186 components of targeted buffer : FD

    auto tr_y_0_x_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 180);

    auto tr_y_0_x_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 181);

    auto tr_y_0_x_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 182);

    auto tr_y_0_x_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 183);

    auto tr_y_0_x_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 184);

    auto tr_y_0_x_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 185);

    #pragma omp simd aligned(tr_xxxxy_xx, tr_xxxxy_xy, tr_xxxxy_xz, tr_xxxxy_yy, tr_xxxxy_yz, tr_xxxxy_zz, tr_xxxy_x, tr_xxxy_xxx, tr_xxxy_xxy, tr_xxxy_xxz, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_y, tr_xxxy_z, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_y_0_x_xxx_xx, tr_y_0_x_xxx_xy, tr_y_0_x_xxx_xz, tr_y_0_x_xxx_yy, tr_y_0_x_xxx_yz, tr_y_0_x_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxx_xx[i] = -6.0 * tr_xxy_xx[i] * tbe_0 - 4.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxy_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xy[i] = -6.0 * tr_xxy_xy[i] * tbe_0 - 2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_xz[i] = -6.0 * tr_xxy_xz[i] * tbe_0 - 2.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_yy[i] = -6.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_yz[i] = -6.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxx_zz[i] = -6.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 186-192 components of targeted buffer : FD

    auto tr_y_0_x_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 186);

    auto tr_y_0_x_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 187);

    auto tr_y_0_x_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 188);

    auto tr_y_0_x_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 189);

    auto tr_y_0_x_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 190);

    auto tr_y_0_x_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 191);

    #pragma omp simd aligned(tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_z, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_xxyy_x, tr_xxyy_xxx, tr_xxyy_xxy, tr_xxyy_xxz, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_y, tr_xxyy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_y_0_x_xxy_xx, tr_y_0_x_xxy_xy, tr_y_0_x_xxy_xz, tr_y_0_x_xxy_yy, tr_y_0_x_xxy_yz, tr_y_0_x_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxy_xx[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_xyy_xx[i] * tbe_0 + 2.0 * tr_xx_x[i] - 2.0 * tr_xx_xxx[i] * tke_0 - 4.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyy_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxxyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xy[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_xyy_xy[i] * tbe_0 + tr_xx_y[i] - 2.0 * tr_xx_xxy[i] * tke_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xy[i] * tbe_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_xz[i] = 2.0 * tr_x_xz[i] - 4.0 * tr_xyy_xz[i] * tbe_0 + tr_xx_z[i] - 2.0 * tr_xx_xxz[i] * tke_0 - 2.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xz[i] * tbe_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_yy[i] = 2.0 * tr_x_yy[i] - 4.0 * tr_xyy_yy[i] * tbe_0 - 2.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxxyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_yz[i] = 2.0 * tr_x_yz[i] - 4.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yz[i] * tbe_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxy_zz[i] = 2.0 * tr_x_zz[i] - 4.0 * tr_xyy_zz[i] * tbe_0 - 2.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 192-198 components of targeted buffer : FD

    auto tr_y_0_x_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 192);

    auto tr_y_0_x_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 193);

    auto tr_y_0_x_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 194);

    auto tr_y_0_x_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 195);

    auto tr_y_0_x_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 196);

    auto tr_y_0_x_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 197);

    #pragma omp simd aligned(tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_x, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_y_0_x_xxz_xx, tr_y_0_x_xxz_xy, tr_y_0_x_xxz_xz, tr_y_0_x_xxz_yy, tr_y_0_x_xxz_yz, tr_y_0_x_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xxz_xx[i] = -4.0 * tr_xyz_xx[i] * tbe_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xy[i] = -4.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_xz[i] = -4.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_yy[i] = -4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_yz[i] = -4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xxz_zz[i] = -4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 198-204 components of targeted buffer : FD

    auto tr_y_0_x_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 198);

    auto tr_y_0_x_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 199);

    auto tr_y_0_x_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 200);

    auto tr_y_0_x_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 201);

    auto tr_y_0_x_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 202);

    auto tr_y_0_x_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 203);

    #pragma omp simd aligned(tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_xy_x, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_z, tr_xyyy_x, tr_xyyy_xxx, tr_xyyy_xxy, tr_xyyy_xxz, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_y, tr_xyyy_z, tr_y_0_x_xyy_xx, tr_y_0_x_xyy_xy, tr_y_0_x_xyy_xz, tr_y_0_x_xyy_yy, tr_y_0_x_xyy_yz, tr_y_0_x_xyy_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyy_xx[i] = 2.0 * tr_y_xx[i] - 2.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_xy_x[i] - 4.0 * tr_xy_xxx[i] * tke_0 - 4.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xy[i] = 2.0 * tr_y_xy[i] - 2.0 * tr_yyy_xy[i] * tbe_0 + 2.0 * tr_xy_y[i] - 4.0 * tr_xy_xxy[i] * tke_0 - 2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyy_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_xz[i] = 2.0 * tr_y_xz[i] - 2.0 * tr_yyy_xz[i] * tbe_0 + 2.0 * tr_xy_z[i] - 4.0 * tr_xy_xxz[i] * tke_0 - 2.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyy_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_yy[i] = 2.0 * tr_y_yy[i] - 2.0 * tr_yyy_yy[i] * tbe_0 - 4.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xyyy_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxyyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_yz[i] = 2.0 * tr_y_yz[i] - 2.0 * tr_yyy_yz[i] * tbe_0 - 4.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyy_zz[i] = 2.0 * tr_y_zz[i] - 2.0 * tr_yyy_zz[i] * tbe_0 - 4.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xyyy_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 204-210 components of targeted buffer : FD

    auto tr_y_0_x_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 204);

    auto tr_y_0_x_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 205);

    auto tr_y_0_x_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 206);

    auto tr_y_0_x_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 207);

    auto tr_y_0_x_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 208);

    auto tr_y_0_x_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 209);

    #pragma omp simd aligned(tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xyyz_x, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_z, tr_xz_x, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_z, tr_y_0_x_xyz_xx, tr_y_0_x_xyz_xy, tr_y_0_x_xyz_xz, tr_y_0_x_xyz_yy, tr_y_0_x_xyz_yz, tr_y_0_x_xyz_zz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xyz_xx[i] = tr_z_xx[i] - 2.0 * tr_yyz_xx[i] * tbe_0 + 2.0 * tr_xz_x[i] - 2.0 * tr_xz_xxx[i] * tke_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xy[i] = tr_z_xy[i] - 2.0 * tr_yyz_xy[i] * tbe_0 + tr_xz_y[i] - 2.0 * tr_xz_xxy[i] * tke_0 - 2.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_xz[i] = tr_z_xz[i] - 2.0 * tr_yyz_xz[i] * tbe_0 + tr_xz_z[i] - 2.0 * tr_xz_xxz[i] * tke_0 - 2.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_yy[i] = tr_z_yy[i] - 2.0 * tr_yyz_yy[i] * tbe_0 - 2.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_yz[i] = tr_z_yz[i] - 2.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xyz_zz[i] = tr_z_zz[i] - 2.0 * tr_yyz_zz[i] * tbe_0 - 2.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 210-216 components of targeted buffer : FD

    auto tr_y_0_x_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 210);

    auto tr_y_0_x_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 211);

    auto tr_y_0_x_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 212);

    auto tr_y_0_x_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 213);

    auto tr_y_0_x_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 214);

    auto tr_y_0_x_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 215);

    #pragma omp simd aligned(tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xyzz_x, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_z, tr_y_0_x_xzz_xx, tr_y_0_x_xzz_xy, tr_y_0_x_xzz_xz, tr_y_0_x_xzz_yy, tr_y_0_x_xzz_yz, tr_y_0_x_xzz_zz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_xzz_xx[i] = -2.0 * tr_yzz_xx[i] * tbe_0 - 4.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xy[i] = -2.0 * tr_yzz_xy[i] * tbe_0 - 2.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_xz[i] = -2.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_yy[i] = -2.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_yz[i] = -2.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_xzz_zz[i] = -2.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 216-222 components of targeted buffer : FD

    auto tr_y_0_x_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 216);

    auto tr_y_0_x_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 217);

    auto tr_y_0_x_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 218);

    auto tr_y_0_x_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 219);

    auto tr_y_0_x_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 220);

    auto tr_y_0_x_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 221);

    #pragma omp simd aligned(tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, tr_y_0_x_yyy_xx, tr_y_0_x_yyy_xy, tr_y_0_x_yyy_xz, tr_y_0_x_yyy_yy, tr_y_0_x_yyy_yz, tr_y_0_x_yyy_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_z, tr_yyyy_x, tr_yyyy_xxx, tr_yyyy_xxy, tr_yyyy_xxz, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_y, tr_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyy_xx[i] = 6.0 * tr_yy_x[i] - 6.0 * tr_yy_xxx[i] * tke_0 - 4.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyy_xxx[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xy[i] = 3.0 * tr_yy_y[i] - 6.0 * tr_yy_xxy[i] * tke_0 - 2.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_yyyy_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_xz[i] = 3.0 * tr_yy_z[i] - 6.0 * tr_yy_xxz[i] * tke_0 - 2.0 * tr_yyyy_z[i] * tbe_0 + 4.0 * tr_yyyy_xxz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyyyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_yy[i] = -6.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yyyy_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyyyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_yz[i] = -6.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yyyy_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyyyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyy_zz[i] = -6.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yyyy_xzz[i] * tbe_0 * tke_0 - 6.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 222-228 components of targeted buffer : FD

    auto tr_y_0_x_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 222);

    auto tr_y_0_x_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 223);

    auto tr_y_0_x_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 224);

    auto tr_y_0_x_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 225);

    auto tr_y_0_x_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 226);

    auto tr_y_0_x_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 227);

    #pragma omp simd aligned(tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_y_0_x_yyz_xx, tr_y_0_x_yyz_xy, tr_y_0_x_yyz_xz, tr_y_0_x_yyz_yy, tr_y_0_x_yyz_yz, tr_y_0_x_yyz_zz, tr_yyyz_x, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_y, tr_yyyz_z, tr_yz_x, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yyz_xx[i] = 4.0 * tr_yz_x[i] - 4.0 * tr_yz_xxx[i] * tke_0 - 4.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xy[i] = 2.0 * tr_yz_y[i] - 4.0 * tr_yz_xxy[i] * tke_0 - 2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_xz[i] = 2.0 * tr_yz_z[i] - 4.0 * tr_yz_xxz[i] * tke_0 - 2.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_yy[i] = -4.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yyyz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_yz[i] = -4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yyz_zz[i] = -4.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yyyz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 228-234 components of targeted buffer : FD

    auto tr_y_0_x_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 228);

    auto tr_y_0_x_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 229);

    auto tr_y_0_x_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 230);

    auto tr_y_0_x_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 231);

    auto tr_y_0_x_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 232);

    auto tr_y_0_x_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 233);

    #pragma omp simd aligned(tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_y_0_x_yzz_xx, tr_y_0_x_yzz_xy, tr_y_0_x_yzz_xz, tr_y_0_x_yzz_yy, tr_y_0_x_yzz_yz, tr_y_0_x_yzz_zz, tr_yyzz_x, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_y, tr_yyzz_z, tr_zz_x, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_yzz_xx[i] = 2.0 * tr_zz_x[i] - 2.0 * tr_zz_xxx[i] * tke_0 - 4.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xy[i] = tr_zz_y[i] - 2.0 * tr_zz_xxy[i] * tke_0 - 2.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_xz[i] = tr_zz_z[i] - 2.0 * tr_zz_xxz[i] * tke_0 - 2.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_yy[i] = -2.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_yz[i] = -2.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_yzz_zz[i] = -2.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 234-240 components of targeted buffer : FD

    auto tr_y_0_x_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 234);

    auto tr_y_0_x_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 235);

    auto tr_y_0_x_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 236);

    auto tr_y_0_x_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 237);

    auto tr_y_0_x_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 238);

    auto tr_y_0_x_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 239);

    #pragma omp simd aligned(tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_y_0_x_zzz_xx, tr_y_0_x_zzz_xy, tr_y_0_x_zzz_xz, tr_y_0_x_zzz_yy, tr_y_0_x_zzz_yz, tr_y_0_x_zzz_zz, tr_yzzz_x, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_y, tr_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_x_zzz_xx[i] = -4.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xy[i] = -2.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_xz[i] = -2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_yy[i] = 4.0 * tr_yzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_yz[i] = 4.0 * tr_yzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_x_zzz_zz[i] = 4.0 * tr_yzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 240-246 components of targeted buffer : FD

    auto tr_y_0_y_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 240);

    auto tr_y_0_y_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 241);

    auto tr_y_0_y_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 242);

    auto tr_y_0_y_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 243);

    auto tr_y_0_y_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 244);

    auto tr_y_0_y_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 245);

    #pragma omp simd aligned(tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxy_x, tr_xxxy_xxy, tr_xxxy_xyy, tr_xxxy_xyz, tr_xxxy_y, tr_xxxy_yyy, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_z, tr_xxxyy_xx, tr_xxxyy_xy, tr_xxxyy_xz, tr_xxxyy_yy, tr_xxxyy_yz, tr_xxxyy_zz, tr_y_0_y_xxx_xx, tr_y_0_y_xxx_xy, tr_y_0_y_xxx_xz, tr_y_0_y_xxx_yy, tr_y_0_y_xxx_yz, tr_y_0_y_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxx_xx[i] = -2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxxy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xy[i] = -2.0 * tr_xxx_xy[i] * tbe_0 - 2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_xz[i] = -2.0 * tr_xxx_xz[i] * tbe_0 + 4.0 * tr_xxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_yy[i] = -2.0 * tr_xxx_yy[i] * tbe_0 - 4.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_yz[i] = -2.0 * tr_xxx_yz[i] * tbe_0 - 2.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxx_zz[i] = -2.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 246-252 components of targeted buffer : FD

    auto tr_y_0_y_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 246);

    auto tr_y_0_y_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 247);

    auto tr_y_0_y_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 248);

    auto tr_y_0_y_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 249);

    auto tr_y_0_y_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 250);

    auto tr_y_0_y_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 251);

    #pragma omp simd aligned(tr_xx_x, tr_xx_xxy, tr_xx_xyy, tr_xx_xyz, tr_xx_y, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyy_x, tr_xxyy_xxy, tr_xxyy_xyy, tr_xxyy_xyz, tr_xxyy_y, tr_xxyy_yyy, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_z, tr_xxyyy_xx, tr_xxyyy_xy, tr_xxyyy_xz, tr_xxyyy_yy, tr_xxyyy_yz, tr_xxyyy_zz, tr_y_0_y_xxy_xx, tr_y_0_y_xxy_xy, tr_y_0_y_xxy_xz, tr_y_0_y_xxy_yy, tr_y_0_y_xxy_yz, tr_y_0_y_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxy_xx[i] = -2.0 * tr_xx_xxy[i] * tke_0 - 6.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xy[i] = tr_xx_x[i] - 2.0 * tr_xx_xyy[i] * tke_0 - 6.0 * tr_xxy_xy[i] * tbe_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_xz[i] = -2.0 * tr_xx_xyz[i] * tke_0 - 6.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_yy[i] = 2.0 * tr_xx_y[i] - 2.0 * tr_xx_yyy[i] * tke_0 - 6.0 * tr_xxy_yy[i] * tbe_0 - 4.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_yz[i] = tr_xx_z[i] - 2.0 * tr_xx_yyz[i] * tke_0 - 6.0 * tr_xxy_yz[i] * tbe_0 - 2.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxy_zz[i] = -2.0 * tr_xx_yzz[i] * tke_0 - 6.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 252-258 components of targeted buffer : FD

    auto tr_y_0_y_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 252);

    auto tr_y_0_y_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 253);

    auto tr_y_0_y_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 254);

    auto tr_y_0_y_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 255);

    auto tr_y_0_y_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 256);

    auto tr_y_0_y_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 257);

    #pragma omp simd aligned(tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyz_x, tr_xxyz_xxy, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_y, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_y_0_y_xxz_xx, tr_y_0_y_xxz_xy, tr_y_0_y_xxz_xz, tr_y_0_y_xxz_yy, tr_y_0_y_xxz_yz, tr_y_0_y_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xxz_xx[i] = -2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xy[i] = -2.0 * tr_xxz_xy[i] * tbe_0 - 2.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_xz[i] = -2.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_yy[i] = -2.0 * tr_xxz_yy[i] * tbe_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_yz[i] = -2.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xxz_zz[i] = -2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 258-264 components of targeted buffer : FD

    auto tr_y_0_y_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 258);

    auto tr_y_0_y_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 259);

    auto tr_y_0_y_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 260);

    auto tr_y_0_y_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 261);

    auto tr_y_0_y_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 262);

    auto tr_y_0_y_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 263);

    #pragma omp simd aligned(tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xy_x, tr_xy_xxy, tr_xy_xyy, tr_xy_xyz, tr_xy_y, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyy_x, tr_xyyy_xxy, tr_xyyy_xyy, tr_xyyy_xyz, tr_xyyy_y, tr_xyyy_yyy, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_z, tr_xyyyy_xx, tr_xyyyy_xy, tr_xyyyy_xz, tr_xyyyy_yy, tr_xyyyy_yz, tr_xyyyy_zz, tr_y_0_y_xyy_xx, tr_y_0_y_xyy_xy, tr_y_0_y_xyy_xz, tr_y_0_y_xyy_yy, tr_y_0_y_xyy_yz, tr_y_0_y_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyy_xx[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_xy_xxy[i] * tke_0 - 10.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xy[i] = 2.0 * tr_x_xy[i] + 2.0 * tr_xy_x[i] - 4.0 * tr_xy_xyy[i] * tke_0 - 10.0 * tr_xyy_xy[i] * tbe_0 - 2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_xz[i] = 2.0 * tr_x_xz[i] - 4.0 * tr_xy_xyz[i] * tke_0 - 10.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_yy[i] = 2.0 * tr_x_yy[i] + 4.0 * tr_xy_y[i] - 4.0 * tr_xy_yyy[i] * tke_0 - 10.0 * tr_xyy_yy[i] * tbe_0 - 4.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_yz[i] = 2.0 * tr_x_yz[i] + 2.0 * tr_xy_z[i] - 4.0 * tr_xy_yyz[i] * tke_0 - 10.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyy_zz[i] = 2.0 * tr_x_zz[i] - 4.0 * tr_xy_yzz[i] * tke_0 - 10.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 264-270 components of targeted buffer : FD

    auto tr_y_0_y_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 264);

    auto tr_y_0_y_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 265);

    auto tr_y_0_y_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 266);

    auto tr_y_0_y_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 267);

    auto tr_y_0_y_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 268);

    auto tr_y_0_y_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 269);

    #pragma omp simd aligned(tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyz_x, tr_xyyz_xxy, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_y, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xz_x, tr_xz_xxy, tr_xz_xyy, tr_xz_xyz, tr_xz_y, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_y_0_y_xyz_xx, tr_y_0_y_xyz_xy, tr_y_0_y_xyz_xz, tr_y_0_y_xyz_yy, tr_y_0_y_xyz_yz, tr_y_0_y_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xyz_xx[i] = -2.0 * tr_xz_xxy[i] * tke_0 - 6.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xy[i] = tr_xz_x[i] - 2.0 * tr_xz_xyy[i] * tke_0 - 6.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_xz[i] = -2.0 * tr_xz_xyz[i] * tke_0 - 6.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_yy[i] = 2.0 * tr_xz_y[i] - 2.0 * tr_xz_yyy[i] * tke_0 - 6.0 * tr_xyz_yy[i] * tbe_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_yz[i] = tr_xz_z[i] - 2.0 * tr_xz_yyz[i] * tke_0 - 6.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xyz_zz[i] = -2.0 * tr_xz_yzz[i] * tke_0 - 6.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 270-276 components of targeted buffer : FD

    auto tr_y_0_y_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 270);

    auto tr_y_0_y_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 271);

    auto tr_y_0_y_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 272);

    auto tr_y_0_y_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 273);

    auto tr_y_0_y_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 274);

    auto tr_y_0_y_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 275);

    #pragma omp simd aligned(tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyzz_x, tr_xyzz_xxy, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_y, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_y_0_y_xzz_xx, tr_y_0_y_xzz_xy, tr_y_0_y_xzz_xz, tr_y_0_y_xzz_yy, tr_y_0_y_xzz_yz, tr_y_0_y_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_xzz_xx[i] = -2.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xy[i] = -2.0 * tr_xzz_xy[i] * tbe_0 - 2.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_xz[i] = -2.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_yy[i] = -2.0 * tr_xzz_yy[i] * tbe_0 - 4.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_yz[i] = -2.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_xzz_zz[i] = -2.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 276-282 components of targeted buffer : FD

    auto tr_y_0_y_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 276);

    auto tr_y_0_y_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 277);

    auto tr_y_0_y_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 278);

    auto tr_y_0_y_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 279);

    auto tr_y_0_y_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 280);

    auto tr_y_0_y_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 281);

    #pragma omp simd aligned(tr_y_0_y_yyy_xx, tr_y_0_y_yyy_xy, tr_y_0_y_yyy_xz, tr_y_0_y_yyy_yy, tr_y_0_y_yyy_yz, tr_y_0_y_yyy_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_xxy, tr_yy_xyy, tr_yy_xyz, tr_yy_y, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyy_x, tr_yyyy_xxy, tr_yyyy_xyy, tr_yyyy_xyz, tr_yyyy_y, tr_yyyy_yyy, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_z, tr_yyyyy_xx, tr_yyyyy_xy, tr_yyyyy_xz, tr_yyyyy_yy, tr_yyyyy_yz, tr_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyy_xx[i] = 6.0 * tr_y_xx[i] - 6.0 * tr_yy_xxy[i] * tke_0 - 14.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_yyyy_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xy[i] = 6.0 * tr_y_xy[i] + 3.0 * tr_yy_x[i] - 6.0 * tr_yy_xyy[i] * tke_0 - 14.0 * tr_yyy_xy[i] * tbe_0 - 2.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyy_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_xz[i] = 6.0 * tr_y_xz[i] - 6.0 * tr_yy_xyz[i] * tke_0 - 14.0 * tr_yyy_xz[i] * tbe_0 + 4.0 * tr_yyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_yy[i] = 6.0 * tr_y_yy[i] + 6.0 * tr_yy_y[i] - 6.0 * tr_yy_yyy[i] * tke_0 - 14.0 * tr_yyy_yy[i] * tbe_0 - 4.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_yyyy_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_yz[i] = 6.0 * tr_y_yz[i] + 3.0 * tr_yy_z[i] - 6.0 * tr_yy_yyz[i] * tke_0 - 14.0 * tr_yyy_yz[i] * tbe_0 - 2.0 * tr_yyyy_z[i] * tbe_0 + 4.0 * tr_yyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyy_zz[i] = 6.0 * tr_y_zz[i] - 6.0 * tr_yy_yzz[i] * tke_0 - 14.0 * tr_yyy_zz[i] * tbe_0 + 4.0 * tr_yyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyy_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 282-288 components of targeted buffer : FD

    auto tr_y_0_y_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 282);

    auto tr_y_0_y_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 283);

    auto tr_y_0_y_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 284);

    auto tr_y_0_y_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 285);

    auto tr_y_0_y_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 286);

    auto tr_y_0_y_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 287);

    #pragma omp simd aligned(tr_y_0_y_yyz_xx, tr_y_0_y_yyz_xy, tr_y_0_y_yyz_xz, tr_y_0_y_yyz_yy, tr_y_0_y_yyz_yz, tr_y_0_y_yyz_zz, tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyyz_x, tr_yyyz_xxy, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_y, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_z, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yz_x, tr_yz_xxy, tr_yz_xyy, tr_yz_xyz, tr_yz_y, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yyz_xx[i] = 2.0 * tr_z_xx[i] - 4.0 * tr_yz_xxy[i] * tke_0 - 10.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xy[i] = 2.0 * tr_z_xy[i] + 2.0 * tr_yz_x[i] - 4.0 * tr_yz_xyy[i] * tke_0 - 10.0 * tr_yyz_xy[i] * tbe_0 - 2.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_xz[i] = 2.0 * tr_z_xz[i] - 4.0 * tr_yz_xyz[i] * tke_0 - 10.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_yy[i] = 2.0 * tr_z_yy[i] + 4.0 * tr_yz_y[i] - 4.0 * tr_yz_yyy[i] * tke_0 - 10.0 * tr_yyz_yy[i] * tbe_0 - 4.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_yz[i] = 2.0 * tr_z_yz[i] + 2.0 * tr_yz_z[i] - 4.0 * tr_yz_yyz[i] * tke_0 - 10.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yyz_zz[i] = 2.0 * tr_z_zz[i] - 4.0 * tr_yz_yzz[i] * tke_0 - 10.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 288-294 components of targeted buffer : FD

    auto tr_y_0_y_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 288);

    auto tr_y_0_y_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 289);

    auto tr_y_0_y_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 290);

    auto tr_y_0_y_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 291);

    auto tr_y_0_y_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 292);

    auto tr_y_0_y_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 293);

    #pragma omp simd aligned(tr_y_0_y_yzz_xx, tr_y_0_y_yzz_xy, tr_y_0_y_yzz_xz, tr_y_0_y_yzz_yy, tr_y_0_y_yzz_yz, tr_y_0_y_yzz_zz, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_yyzz_x, tr_yyzz_xxy, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_y, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_z, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_zz_x, tr_zz_xxy, tr_zz_xyy, tr_zz_xyz, tr_zz_y, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_yzz_xx[i] = -2.0 * tr_zz_xxy[i] * tke_0 - 6.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xy[i] = tr_zz_x[i] - 2.0 * tr_zz_xyy[i] * tke_0 - 6.0 * tr_yzz_xy[i] * tbe_0 - 2.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_xz[i] = -2.0 * tr_zz_xyz[i] * tke_0 - 6.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_yy[i] = 2.0 * tr_zz_y[i] - 2.0 * tr_zz_yyy[i] * tke_0 - 6.0 * tr_yzz_yy[i] * tbe_0 - 4.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_yz[i] = tr_zz_z[i] - 2.0 * tr_zz_yyz[i] * tke_0 - 6.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_yzz_zz[i] = -2.0 * tr_zz_yzz[i] * tke_0 - 6.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 294-300 components of targeted buffer : FD

    auto tr_y_0_y_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 294);

    auto tr_y_0_y_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 295);

    auto tr_y_0_y_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 296);

    auto tr_y_0_y_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 297);

    auto tr_y_0_y_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 298);

    auto tr_y_0_y_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 299);

    #pragma omp simd aligned(tr_y_0_y_zzz_xx, tr_y_0_y_zzz_xy, tr_y_0_y_zzz_xz, tr_y_0_y_zzz_yy, tr_y_0_y_zzz_yz, tr_y_0_y_zzz_zz, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_yzzz_x, tr_yzzz_xxy, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_y, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_z, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_y_zzz_xx[i] = -2.0 * tr_zzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xy[i] = -2.0 * tr_zzz_xy[i] * tbe_0 - 2.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_xz[i] = -2.0 * tr_zzz_xz[i] * tbe_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_yy[i] = -2.0 * tr_zzz_yy[i] * tbe_0 - 4.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_yz[i] = -2.0 * tr_zzz_yz[i] * tbe_0 - 2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_y_zzz_zz[i] = -2.0 * tr_zzz_zz[i] * tbe_0 + 4.0 * tr_yzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 300-306 components of targeted buffer : FD

    auto tr_y_0_z_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 300);

    auto tr_y_0_z_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 301);

    auto tr_y_0_z_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 302);

    auto tr_y_0_z_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 303);

    auto tr_y_0_z_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 304);

    auto tr_y_0_z_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 305);

    #pragma omp simd aligned(tr_xxxy_x, tr_xxxy_xxz, tr_xxxy_xyz, tr_xxxy_xzz, tr_xxxy_y, tr_xxxy_yyz, tr_xxxy_yzz, tr_xxxy_z, tr_xxxy_zzz, tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_y_0_z_xxx_xx, tr_y_0_z_xxx_xy, tr_y_0_z_xxx_xz, tr_y_0_z_xxx_yy, tr_y_0_z_xxx_yz, tr_y_0_z_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxx_xx[i] = 4.0 * tr_xxxy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xy[i] = 4.0 * tr_xxxy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_xz[i] = -2.0 * tr_xxxy_x[i] * tbe_0 + 4.0 * tr_xxxy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_yy[i] = 4.0 * tr_xxxy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_yz[i] = -2.0 * tr_xxxy_y[i] * tbe_0 + 4.0 * tr_xxxy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxx_zz[i] = -4.0 * tr_xxxy_z[i] * tbe_0 + 4.0 * tr_xxxy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 306-312 components of targeted buffer : FD

    auto tr_y_0_z_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 306);

    auto tr_y_0_z_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 307);

    auto tr_y_0_z_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 308);

    auto tr_y_0_z_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 309);

    auto tr_y_0_z_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 310);

    auto tr_y_0_z_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 311);

    #pragma omp simd aligned(tr_xx_x, tr_xx_xxz, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxyy_x, tr_xxyy_xxz, tr_xxyy_xyz, tr_xxyy_xzz, tr_xxyy_y, tr_xxyy_yyz, tr_xxyy_yzz, tr_xxyy_z, tr_xxyy_zzz, tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_y_0_z_xxy_xx, tr_y_0_z_xxy_xy, tr_y_0_z_xxy_xz, tr_y_0_z_xxy_yy, tr_y_0_z_xxy_yz, tr_y_0_z_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxy_xx[i] = -2.0 * tr_xx_xxz[i] * tke_0 - 2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xy[i] = -2.0 * tr_xx_xyz[i] * tke_0 - 2.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_xz[i] = tr_xx_x[i] - 2.0 * tr_xx_xzz[i] * tke_0 - 2.0 * tr_xxz_xz[i] * tbe_0 - 2.0 * tr_xxyy_x[i] * tbe_0 + 4.0 * tr_xxyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_yy[i] = -2.0 * tr_xx_yyz[i] * tke_0 - 2.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_yz[i] = tr_xx_y[i] - 2.0 * tr_xx_yzz[i] * tke_0 - 2.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxyy_y[i] * tbe_0 + 4.0 * tr_xxyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxy_zz[i] = 2.0 * tr_xx_z[i] - 2.0 * tr_xx_zzz[i] * tke_0 - 2.0 * tr_xxz_zz[i] * tbe_0 - 4.0 * tr_xxyy_z[i] * tbe_0 + 4.0 * tr_xxyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 312-318 components of targeted buffer : FD

    auto tr_y_0_z_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 312);

    auto tr_y_0_z_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 313);

    auto tr_y_0_z_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 314);

    auto tr_y_0_z_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 315);

    auto tr_y_0_z_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 316);

    auto tr_y_0_z_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 317);

    #pragma omp simd aligned(tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_xxz, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_y_0_z_xxz_xx, tr_y_0_z_xxz_xy, tr_y_0_z_xxz_xz, tr_y_0_z_xxz_yy, tr_y_0_z_xxz_yz, tr_y_0_z_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xxz_xx[i] = -2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xy[i] = -2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_xz[i] = -2.0 * tr_xxy_xz[i] * tbe_0 - 2.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_yy[i] = -2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_yz[i] = -2.0 * tr_xxy_yz[i] * tbe_0 - 2.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xxz_zz[i] = -2.0 * tr_xxy_zz[i] * tbe_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 318-324 components of targeted buffer : FD

    auto tr_y_0_z_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 318);

    auto tr_y_0_z_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 319);

    auto tr_y_0_z_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 320);

    auto tr_y_0_z_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 321);

    auto tr_y_0_z_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 322);

    auto tr_y_0_z_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 323);

    #pragma omp simd aligned(tr_xy_x, tr_xy_xxz, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_xyyy_x, tr_xyyy_xxz, tr_xyyy_xyz, tr_xyyy_xzz, tr_xyyy_y, tr_xyyy_yyz, tr_xyyy_yzz, tr_xyyy_z, tr_xyyy_zzz, tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_y_0_z_xyy_xx, tr_y_0_z_xyy_xy, tr_y_0_z_xyy_xz, tr_y_0_z_xyy_yy, tr_y_0_z_xyy_yz, tr_y_0_z_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyy_xx[i] = -4.0 * tr_xy_xxz[i] * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xy[i] = -4.0 * tr_xy_xyz[i] * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_xz[i] = 2.0 * tr_xy_x[i] - 4.0 * tr_xy_xzz[i] * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xyyy_x[i] * tbe_0 + 4.0 * tr_xyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_yy[i] = -4.0 * tr_xy_yyz[i] * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_yz[i] = 2.0 * tr_xy_y[i] - 4.0 * tr_xy_yzz[i] * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyyy_y[i] * tbe_0 + 4.0 * tr_xyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyy_zz[i] = 4.0 * tr_xy_z[i] - 4.0 * tr_xy_zzz[i] * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 - 4.0 * tr_xyyy_z[i] * tbe_0 + 4.0 * tr_xyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 324-330 components of targeted buffer : FD

    auto tr_y_0_z_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 324);

    auto tr_y_0_z_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 325);

    auto tr_y_0_z_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 326);

    auto tr_y_0_z_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 327);

    auto tr_y_0_z_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 328);

    auto tr_y_0_z_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 329);

    #pragma omp simd aligned(tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_xxz, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xz_x, tr_xz_xxz, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_y_0_z_xyz_xx, tr_y_0_z_xyz_xy, tr_y_0_z_xyz_xz, tr_y_0_z_xyz_yy, tr_y_0_z_xyz_yz, tr_y_0_z_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xyz_xx[i] = tr_x_xx[i] - 2.0 * tr_xz_xxz[i] * tke_0 - 2.0 * tr_xzz_xx[i] * tbe_0 - 2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xy[i] = tr_x_xy[i] - 2.0 * tr_xz_xyz[i] * tke_0 - 2.0 * tr_xzz_xy[i] * tbe_0 - 2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_xz[i] = tr_x_xz[i] + tr_xz_x[i] - 2.0 * tr_xz_xzz[i] * tke_0 - 2.0 * tr_xzz_xz[i] * tbe_0 - 2.0 * tr_xyy_xz[i] * tbe_0 - 2.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_yy[i] = tr_x_yy[i] - 2.0 * tr_xz_yyz[i] * tke_0 - 2.0 * tr_xzz_yy[i] * tbe_0 - 2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_yz[i] = tr_x_yz[i] + tr_xz_y[i] - 2.0 * tr_xz_yzz[i] * tke_0 - 2.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xyz_zz[i] = tr_x_zz[i] + 2.0 * tr_xz_z[i] - 2.0 * tr_xz_zzz[i] * tke_0 - 2.0 * tr_xzz_zz[i] * tbe_0 - 2.0 * tr_xyy_zz[i] * tbe_0 - 4.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 330-336 components of targeted buffer : FD

    auto tr_y_0_z_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 330);

    auto tr_y_0_z_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 331);

    auto tr_y_0_z_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 332);

    auto tr_y_0_z_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 333);

    auto tr_y_0_z_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 334);

    auto tr_y_0_z_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 335);

    #pragma omp simd aligned(tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_xxz, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_y_0_z_xzz_xx, tr_y_0_z_xzz_xy, tr_y_0_z_xzz_xz, tr_y_0_z_xzz_yy, tr_y_0_z_xzz_yz, tr_y_0_z_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_xzz_xx[i] = -4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xy[i] = -4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_xz[i] = -4.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_yy[i] = -4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_yz[i] = -4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_xzz_zz[i] = -4.0 * tr_xyz_zz[i] * tbe_0 - 4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 336-342 components of targeted buffer : FD

    auto tr_y_0_z_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 336);

    auto tr_y_0_z_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 337);

    auto tr_y_0_z_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 338);

    auto tr_y_0_z_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 339);

    auto tr_y_0_z_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 340);

    auto tr_y_0_z_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 341);

    #pragma omp simd aligned(tr_y_0_z_yyy_xx, tr_y_0_z_yyy_xy, tr_y_0_z_yyy_xz, tr_y_0_z_yyy_yy, tr_y_0_z_yyy_yz, tr_y_0_z_yyy_zz, tr_yy_x, tr_yy_xxz, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, tr_yyyy_x, tr_yyyy_xxz, tr_yyyy_xyz, tr_yyyy_xzz, tr_yyyy_y, tr_yyyy_yyz, tr_yyyy_yzz, tr_yyyy_z, tr_yyyy_zzz, tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyy_xx[i] = -6.0 * tr_yy_xxz[i] * tke_0 - 6.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyyy_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xy[i] = -6.0 * tr_yy_xyz[i] * tke_0 - 6.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyyy_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_xz[i] = 3.0 * tr_yy_x[i] - 6.0 * tr_yy_xzz[i] * tke_0 - 6.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_yyyy_x[i] * tbe_0 + 4.0 * tr_yyyy_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_yy[i] = -6.0 * tr_yy_yyz[i] * tke_0 - 6.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyyy_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_yz[i] = 3.0 * tr_yy_y[i] - 6.0 * tr_yy_yzz[i] * tke_0 - 6.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_yyyy_y[i] * tbe_0 + 4.0 * tr_yyyy_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyy_zz[i] = 6.0 * tr_yy_z[i] - 6.0 * tr_yy_zzz[i] * tke_0 - 6.0 * tr_yyz_zz[i] * tbe_0 - 4.0 * tr_yyyy_z[i] * tbe_0 + 4.0 * tr_yyyy_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 342-348 components of targeted buffer : FD

    auto tr_y_0_z_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 342);

    auto tr_y_0_z_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 343);

    auto tr_y_0_z_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 344);

    auto tr_y_0_z_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 345);

    auto tr_y_0_z_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 346);

    auto tr_y_0_z_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 347);

    #pragma omp simd aligned(tr_y_0_z_yyz_xx, tr_y_0_z_yyz_xy, tr_y_0_z_yyz_xz, tr_y_0_z_yyz_yy, tr_y_0_z_yyz_yz, tr_y_0_z_yyz_zz, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_xxz, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_y, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_z, tr_yyyz_zzz, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_yz_x, tr_yz_xxz, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yyz_xx[i] = 2.0 * tr_y_xx[i] - 4.0 * tr_yz_xxz[i] * tke_0 - 4.0 * tr_yzz_xx[i] * tbe_0 - 2.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xy[i] = 2.0 * tr_y_xy[i] - 4.0 * tr_yz_xyz[i] * tke_0 - 4.0 * tr_yzz_xy[i] * tbe_0 - 2.0 * tr_yyy_xy[i] * tbe_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_xz[i] = 2.0 * tr_y_xz[i] + 2.0 * tr_yz_x[i] - 4.0 * tr_yz_xzz[i] * tke_0 - 4.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_yyy_xz[i] * tbe_0 - 2.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_yy[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_yz_yyz[i] * tke_0 - 4.0 * tr_yzz_yy[i] * tbe_0 - 2.0 * tr_yyy_yy[i] * tbe_0 + 4.0 * tr_yyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_yz[i] = 2.0 * tr_y_yz[i] + 2.0 * tr_yz_y[i] - 4.0 * tr_yz_yzz[i] * tke_0 - 4.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_yyy_yz[i] * tbe_0 - 2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yyz_zz[i] = 2.0 * tr_y_zz[i] + 4.0 * tr_yz_z[i] - 4.0 * tr_yz_zzz[i] * tke_0 - 4.0 * tr_yzz_zz[i] * tbe_0 - 2.0 * tr_yyy_zz[i] * tbe_0 - 4.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 348-354 components of targeted buffer : FD

    auto tr_y_0_z_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 348);

    auto tr_y_0_z_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 349);

    auto tr_y_0_z_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 350);

    auto tr_y_0_z_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 351);

    auto tr_y_0_z_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 352);

    auto tr_y_0_z_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 353);

    #pragma omp simd aligned(tr_y_0_z_yzz_xx, tr_y_0_z_yzz_xy, tr_y_0_z_yzz_xz, tr_y_0_z_yzz_yy, tr_y_0_z_yzz_yz, tr_y_0_z_yzz_zz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_xxz, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_y, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_z, tr_yyzz_zzz, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_xxz, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_yzz_xx[i] = 2.0 * tr_z_xx[i] - 2.0 * tr_zz_xxz[i] * tke_0 - 2.0 * tr_zzz_xx[i] * tbe_0 - 4.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xy[i] = 2.0 * tr_z_xy[i] - 2.0 * tr_zz_xyz[i] * tke_0 - 2.0 * tr_zzz_xy[i] * tbe_0 - 4.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_xz[i] = 2.0 * tr_z_xz[i] + tr_zz_x[i] - 2.0 * tr_zz_xzz[i] * tke_0 - 2.0 * tr_zzz_xz[i] * tbe_0 - 4.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_yy[i] = 2.0 * tr_z_yy[i] - 2.0 * tr_zz_yyz[i] * tke_0 - 2.0 * tr_zzz_yy[i] * tbe_0 - 4.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_yz[i] = 2.0 * tr_z_yz[i] + tr_zz_y[i] - 2.0 * tr_zz_yzz[i] * tke_0 - 2.0 * tr_zzz_yz[i] * tbe_0 - 4.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_yzz_zz[i] = 2.0 * tr_z_zz[i] + 2.0 * tr_zz_z[i] - 2.0 * tr_zz_zzz[i] * tke_0 - 2.0 * tr_zzz_zz[i] * tbe_0 - 4.0 * tr_yyz_zz[i] * tbe_0 - 4.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 354-360 components of targeted buffer : FD

    auto tr_y_0_z_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 354);

    auto tr_y_0_z_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 355);

    auto tr_y_0_z_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 356);

    auto tr_y_0_z_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 357);

    auto tr_y_0_z_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 358);

    auto tr_y_0_z_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 359);

    #pragma omp simd aligned(tr_y_0_z_zzz_xx, tr_y_0_z_zzz_xy, tr_y_0_z_zzz_xz, tr_y_0_z_zzz_yy, tr_y_0_z_zzz_yz, tr_y_0_z_zzz_zz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_xxz, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_y, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_z, tr_yzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_y_0_z_zzz_xx[i] = -6.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xx[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xy[i] = -6.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_xz[i] = -6.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_yy[i] = -6.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yy[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_yz[i] = -6.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tbe_0;

        tr_y_0_z_zzz_zz[i] = -6.0 * tr_yzz_zz[i] * tbe_0 - 4.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 360-366 components of targeted buffer : FD

    auto tr_z_0_x_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 360);

    auto tr_z_0_x_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 361);

    auto tr_z_0_x_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 362);

    auto tr_z_0_x_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 363);

    auto tr_z_0_x_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 364);

    auto tr_z_0_x_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 365);

    #pragma omp simd aligned(tr_xxxxz_xx, tr_xxxxz_xy, tr_xxxxz_xz, tr_xxxxz_yy, tr_xxxxz_yz, tr_xxxxz_zz, tr_xxxz_x, tr_xxxz_xxx, tr_xxxz_xxy, tr_xxxz_xxz, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_y, tr_xxxz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_z_0_x_xxx_xx, tr_z_0_x_xxx_xy, tr_z_0_x_xxx_xz, tr_z_0_x_xxx_yy, tr_z_0_x_xxx_yz, tr_z_0_x_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxx_xx[i] = -6.0 * tr_xxz_xx[i] * tbe_0 - 4.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xy[i] = -6.0 * tr_xxz_xy[i] * tbe_0 - 2.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_xz[i] = -6.0 * tr_xxz_xz[i] * tbe_0 - 2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_yy[i] = -6.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_yz[i] = -6.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxx_zz[i] = -6.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxxz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 366-372 components of targeted buffer : FD

    auto tr_z_0_x_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 366);

    auto tr_z_0_x_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 367);

    auto tr_z_0_x_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 368);

    auto tr_z_0_x_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 369);

    auto tr_z_0_x_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 370);

    auto tr_z_0_x_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 371);

    #pragma omp simd aligned(tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxyz_x, tr_xxyz_xxx, tr_xxyz_xxy, tr_xxyz_xxz, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_z_0_x_xxy_xx, tr_z_0_x_xxy_xy, tr_z_0_x_xxy_xz, tr_z_0_x_xxy_yy, tr_z_0_x_xxy_yz, tr_z_0_x_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxy_xx[i] = -4.0 * tr_xyz_xx[i] * tbe_0 - 4.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xy[i] = -4.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_xz[i] = -4.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_yy[i] = -4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_yz[i] = -4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxy_zz[i] = -4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 372-378 components of targeted buffer : FD

    auto tr_z_0_x_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 372);

    auto tr_z_0_x_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 373);

    auto tr_z_0_x_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 374);

    auto tr_z_0_x_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 375);

    auto tr_z_0_x_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 376);

    auto tr_z_0_x_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 377);

    #pragma omp simd aligned(tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xx_x, tr_xx_xxx, tr_xx_xxy, tr_xx_xxz, tr_xx_xyy, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_z, tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_xxzz_x, tr_xxzz_xxx, tr_xxzz_xxy, tr_xxzz_xxz, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_y, tr_xxzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_z_0_x_xxz_xx, tr_z_0_x_xxz_xy, tr_z_0_x_xxz_xz, tr_z_0_x_xxz_yy, tr_z_0_x_xxz_yz, tr_z_0_x_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xxz_xx[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_xzz_xx[i] * tbe_0 + 2.0 * tr_xx_x[i] - 2.0 * tr_xx_xxx[i] * tke_0 - 4.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxxzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xy[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_xzz_xy[i] * tbe_0 + tr_xx_y[i] - 2.0 * tr_xx_xxy[i] * tke_0 - 2.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xy[i] * tbe_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_xz[i] = 2.0 * tr_x_xz[i] - 4.0 * tr_xzz_xz[i] * tbe_0 + tr_xx_z[i] - 2.0 * tr_xx_xxz[i] * tke_0 - 2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_xz[i] * tbe_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_yy[i] = 2.0 * tr_x_yy[i] - 4.0 * tr_xzz_yy[i] * tbe_0 - 2.0 * tr_xx_xyy[i] * tke_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxxzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_yz[i] = 2.0 * tr_x_yz[i] - 4.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_yz[i] * tbe_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xxz_zz[i] = 2.0 * tr_x_zz[i] - 4.0 * tr_xzz_zz[i] * tbe_0 - 2.0 * tr_xx_xzz[i] * tke_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxx_zz[i] * tbe_0 + 4.0 * tr_xxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 378-384 components of targeted buffer : FD

    auto tr_z_0_x_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 378);

    auto tr_z_0_x_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 379);

    auto tr_z_0_x_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 380);

    auto tr_z_0_x_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 381);

    auto tr_z_0_x_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 382);

    auto tr_z_0_x_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 383);

    #pragma omp simd aligned(tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xyyz_x, tr_xyyz_xxx, tr_xyyz_xxy, tr_xyyz_xxz, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_z, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_z_0_x_xyy_xx, tr_z_0_x_xyy_xy, tr_z_0_x_xyy_xz, tr_z_0_x_xyy_yy, tr_z_0_x_xyy_yz, tr_z_0_x_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyy_xx[i] = -2.0 * tr_yyz_xx[i] * tbe_0 - 4.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xy[i] = -2.0 * tr_yyz_xy[i] * tbe_0 - 2.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_xz[i] = -2.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_yy[i] = -2.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_yz[i] = -2.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyy_zz[i] = -2.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 384-390 components of targeted buffer : FD

    auto tr_z_0_x_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 384);

    auto tr_z_0_x_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 385);

    auto tr_z_0_x_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 386);

    auto tr_z_0_x_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 387);

    auto tr_z_0_x_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 388);

    auto tr_z_0_x_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 389);

    #pragma omp simd aligned(tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xy_x, tr_xy_xxx, tr_xy_xxy, tr_xy_xxz, tr_xy_xyy, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_z, tr_xyzz_x, tr_xyzz_xxx, tr_xyzz_xxy, tr_xyzz_xxz, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_z, tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_z_0_x_xyz_xx, tr_z_0_x_xyz_xy, tr_z_0_x_xyz_xz, tr_z_0_x_xyz_yy, tr_z_0_x_xyz_yz, tr_z_0_x_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xyz_xx[i] = tr_y_xx[i] - 2.0 * tr_yzz_xx[i] * tbe_0 + 2.0 * tr_xy_x[i] - 2.0 * tr_xy_xxx[i] * tke_0 - 4.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xy[i] = tr_y_xy[i] - 2.0 * tr_yzz_xy[i] * tbe_0 + tr_xy_y[i] - 2.0 * tr_xy_xxy[i] * tke_0 - 2.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_xz[i] = tr_y_xz[i] - 2.0 * tr_yzz_xz[i] * tbe_0 + tr_xy_z[i] - 2.0 * tr_xy_xxz[i] * tke_0 - 2.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_yy[i] = tr_y_yy[i] - 2.0 * tr_yzz_yy[i] * tbe_0 - 2.0 * tr_xy_xyy[i] * tke_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_yz[i] = tr_y_yz[i] - 2.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xyz_zz[i] = tr_y_zz[i] - 2.0 * tr_yzz_zz[i] * tbe_0 - 2.0 * tr_xy_xzz[i] * tke_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 390-396 components of targeted buffer : FD

    auto tr_z_0_x_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 390);

    auto tr_z_0_x_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 391);

    auto tr_z_0_x_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 392);

    auto tr_z_0_x_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 393);

    auto tr_z_0_x_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 394);

    auto tr_z_0_x_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 395);

    #pragma omp simd aligned(tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_xz_x, tr_xz_xxx, tr_xz_xxy, tr_xz_xxz, tr_xz_xyy, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_z, tr_xzzz_x, tr_xzzz_xxx, tr_xzzz_xxy, tr_xzzz_xxz, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_y, tr_xzzz_z, tr_z_0_x_xzz_xx, tr_z_0_x_xzz_xy, tr_z_0_x_xzz_xz, tr_z_0_x_xzz_yy, tr_z_0_x_xzz_yz, tr_z_0_x_xzz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_xzz_xx[i] = 2.0 * tr_z_xx[i] - 2.0 * tr_zzz_xx[i] * tbe_0 + 4.0 * tr_xz_x[i] - 4.0 * tr_xz_xxx[i] * tke_0 - 4.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xy[i] = 2.0 * tr_z_xy[i] - 2.0 * tr_zzz_xy[i] * tbe_0 + 2.0 * tr_xz_y[i] - 4.0 * tr_xz_xxy[i] * tke_0 - 2.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_xz[i] = 2.0 * tr_z_xz[i] - 2.0 * tr_zzz_xz[i] * tbe_0 + 2.0 * tr_xz_z[i] - 4.0 * tr_xz_xxz[i] * tke_0 - 2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_yy[i] = 2.0 * tr_z_yy[i] - 2.0 * tr_zzz_yy[i] * tbe_0 - 4.0 * tr_xz_xyy[i] * tke_0 + 4.0 * tr_xzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_yz[i] = 2.0 * tr_z_yz[i] - 2.0 * tr_zzz_yz[i] * tbe_0 - 4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_yz[i] * tbe_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_xzz_zz[i] = 2.0 * tr_z_zz[i] - 2.0 * tr_zzz_zz[i] * tbe_0 - 4.0 * tr_xz_xzz[i] * tke_0 + 4.0 * tr_xzzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 396-402 components of targeted buffer : FD

    auto tr_z_0_x_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 396);

    auto tr_z_0_x_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 397);

    auto tr_z_0_x_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 398);

    auto tr_z_0_x_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 399);

    auto tr_z_0_x_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 400);

    auto tr_z_0_x_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 401);

    #pragma omp simd aligned(tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_yyyz_x, tr_yyyz_xxx, tr_yyyz_xxy, tr_yyyz_xxz, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_y, tr_yyyz_z, tr_z_0_x_yyy_xx, tr_z_0_x_yyy_xy, tr_z_0_x_yyy_xz, tr_z_0_x_yyy_yy, tr_z_0_x_yyy_yz, tr_z_0_x_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyy_xx[i] = -4.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xxx[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xy[i] = -2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_xz[i] = -2.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_yy[i] = 4.0 * tr_yyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_yz[i] = 4.0 * tr_yyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyy_zz[i] = 4.0 * tr_yyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 402-408 components of targeted buffer : FD

    auto tr_z_0_x_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 402);

    auto tr_z_0_x_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 403);

    auto tr_z_0_x_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 404);

    auto tr_z_0_x_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 405);

    auto tr_z_0_x_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 406);

    auto tr_z_0_x_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 407);

    #pragma omp simd aligned(tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_yy_x, tr_yy_xxx, tr_yy_xxy, tr_yy_xxz, tr_yy_xyy, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_z, tr_yyzz_x, tr_yyzz_xxx, tr_yyzz_xxy, tr_yyzz_xxz, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_y, tr_yyzz_z, tr_z_0_x_yyz_xx, tr_z_0_x_yyz_xy, tr_z_0_x_yyz_xz, tr_z_0_x_yyz_yy, tr_z_0_x_yyz_yz, tr_z_0_x_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yyz_xx[i] = 2.0 * tr_yy_x[i] - 2.0 * tr_yy_xxx[i] * tke_0 - 4.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xxx[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xy[i] = tr_yy_y[i] - 2.0 * tr_yy_xxy[i] * tke_0 - 2.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_xz[i] = tr_yy_z[i] - 2.0 * tr_yy_xxz[i] * tke_0 - 2.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_yy[i] = -2.0 * tr_yy_xyy[i] * tke_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_yz[i] = -2.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yyz_zz[i] = -2.0 * tr_yy_xzz[i] * tke_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 408-414 components of targeted buffer : FD

    auto tr_z_0_x_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 408);

    auto tr_z_0_x_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 409);

    auto tr_z_0_x_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 410);

    auto tr_z_0_x_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 411);

    auto tr_z_0_x_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 412);

    auto tr_z_0_x_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 413);

    #pragma omp simd aligned(tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_yz_x, tr_yz_xxx, tr_yz_xxy, tr_yz_xxz, tr_yz_xyy, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_z, tr_yzzz_x, tr_yzzz_xxx, tr_yzzz_xxy, tr_yzzz_xxz, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_y, tr_yzzz_z, tr_z_0_x_yzz_xx, tr_z_0_x_yzz_xy, tr_z_0_x_yzz_xz, tr_z_0_x_yzz_yy, tr_z_0_x_yzz_yz, tr_z_0_x_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_yzz_xx[i] = 4.0 * tr_yz_x[i] - 4.0 * tr_yz_xxx[i] * tke_0 - 4.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xxx[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xy[i] = 2.0 * tr_yz_y[i] - 4.0 * tr_yz_xxy[i] * tke_0 - 2.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_xz[i] = 2.0 * tr_yz_z[i] - 4.0 * tr_yz_xxz[i] * tke_0 - 2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzz_xxz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_yy[i] = -4.0 * tr_yz_xyy[i] * tke_0 + 4.0 * tr_yzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_yz[i] = -4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_yzz_zz[i] = -4.0 * tr_yz_xzz[i] * tke_0 + 4.0 * tr_yzzz_xzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 414-420 components of targeted buffer : FD

    auto tr_z_0_x_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 414);

    auto tr_z_0_x_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 415);

    auto tr_z_0_x_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 416);

    auto tr_z_0_x_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 417);

    auto tr_z_0_x_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 418);

    auto tr_z_0_x_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 419);

    #pragma omp simd aligned(tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, tr_z_0_x_zzz_xx, tr_z_0_x_zzz_xy, tr_z_0_x_zzz_xz, tr_z_0_x_zzz_yy, tr_z_0_x_zzz_yz, tr_z_0_x_zzz_zz, tr_zz_x, tr_zz_xxx, tr_zz_xxy, tr_zz_xxz, tr_zz_xyy, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_z, tr_zzzz_x, tr_zzzz_xxx, tr_zzzz_xxy, tr_zzzz_xxz, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_y, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_x_zzz_xx[i] = 6.0 * tr_zz_x[i] - 6.0 * tr_zz_xxx[i] * tke_0 - 4.0 * tr_zzzz_x[i] * tbe_0 + 4.0 * tr_zzzz_xxx[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xy[i] = 3.0 * tr_zz_y[i] - 6.0 * tr_zz_xxy[i] * tke_0 - 2.0 * tr_zzzz_y[i] * tbe_0 + 4.0 * tr_zzzz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_xz[i] = 3.0 * tr_zz_z[i] - 6.0 * tr_zz_xxz[i] * tke_0 - 2.0 * tr_zzzz_z[i] * tbe_0 + 4.0 * tr_zzzz_xxz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_xz[i] * tbe_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_yy[i] = -6.0 * tr_zz_xyy[i] * tke_0 + 4.0 * tr_zzzz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_yz[i] = -6.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zzzz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_yz[i] * tbe_0 + 4.0 * tr_xzzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_x_zzz_zz[i] = -6.0 * tr_zz_xzz[i] * tke_0 + 4.0 * tr_zzzz_xzz[i] * tbe_0 * tke_0 - 6.0 * tr_xzz_zz[i] * tbe_0 + 4.0 * tr_xzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 420-426 components of targeted buffer : FD

    auto tr_z_0_y_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 420);

    auto tr_z_0_y_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 421);

    auto tr_z_0_y_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 422);

    auto tr_z_0_y_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 423);

    auto tr_z_0_y_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 424);

    auto tr_z_0_y_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 425);

    #pragma omp simd aligned(tr_xxxyz_xx, tr_xxxyz_xy, tr_xxxyz_xz, tr_xxxyz_yy, tr_xxxyz_yz, tr_xxxyz_zz, tr_xxxz_x, tr_xxxz_xxy, tr_xxxz_xyy, tr_xxxz_xyz, tr_xxxz_y, tr_xxxz_yyy, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_z, tr_z_0_y_xxx_xx, tr_z_0_y_xxx_xy, tr_z_0_y_xxx_xz, tr_z_0_y_xxx_yy, tr_z_0_y_xxx_yz, tr_z_0_y_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxx_xx[i] = 4.0 * tr_xxxz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xy[i] = -2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_xz[i] = 4.0 * tr_xxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_yy[i] = -4.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_yz[i] = -2.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxx_zz[i] = 4.0 * tr_xxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 426-432 components of targeted buffer : FD

    auto tr_z_0_y_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 426);

    auto tr_z_0_y_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 427);

    auto tr_z_0_y_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 428);

    auto tr_z_0_y_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 429);

    auto tr_z_0_y_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 430);

    auto tr_z_0_y_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 431);

    #pragma omp simd aligned(tr_xxyyz_xx, tr_xxyyz_xy, tr_xxyyz_xz, tr_xxyyz_yy, tr_xxyyz_yz, tr_xxyyz_zz, tr_xxyz_x, tr_xxyz_xxy, tr_xxyz_xyy, tr_xxyz_xyz, tr_xxyz_y, tr_xxyz_yyy, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_z_0_y_xxy_xx, tr_z_0_y_xxy_xy, tr_z_0_y_xxy_xz, tr_z_0_y_xxy_yy, tr_z_0_y_xxy_yz, tr_z_0_y_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxy_xx[i] = -2.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xy[i] = -2.0 * tr_xxz_xy[i] * tbe_0 - 2.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_xz[i] = -2.0 * tr_xxz_xz[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_yy[i] = -2.0 * tr_xxz_yy[i] * tbe_0 - 4.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_yz[i] = -2.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxy_zz[i] = -2.0 * tr_xxz_zz[i] * tbe_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 432-438 components of targeted buffer : FD

    auto tr_z_0_y_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 432);

    auto tr_z_0_y_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 433);

    auto tr_z_0_y_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 434);

    auto tr_z_0_y_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 435);

    auto tr_z_0_y_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 436);

    auto tr_z_0_y_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 437);

    #pragma omp simd aligned(tr_xx_x, tr_xx_xxy, tr_xx_xyy, tr_xx_xyz, tr_xx_y, tr_xx_yyy, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_xxzz_x, tr_xxzz_xxy, tr_xxzz_xyy, tr_xxzz_xyz, tr_xxzz_y, tr_xxzz_yyy, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_z, tr_z_0_y_xxz_xx, tr_z_0_y_xxz_xy, tr_z_0_y_xxz_xz, tr_z_0_y_xxz_yy, tr_z_0_y_xxz_yz, tr_z_0_y_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xxz_xx[i] = -2.0 * tr_xx_xxy[i] * tke_0 + 4.0 * tr_xxzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xy[i] = tr_xx_x[i] - 2.0 * tr_xx_xyy[i] * tke_0 - 2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_xz[i] = -2.0 * tr_xx_xyz[i] * tke_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_xz[i] * tbe_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_yy[i] = 2.0 * tr_xx_y[i] - 2.0 * tr_xx_yyy[i] * tke_0 - 4.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_yz[i] = tr_xx_z[i] - 2.0 * tr_xx_yyz[i] * tke_0 - 2.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_yz[i] * tbe_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xxz_zz[i] = -2.0 * tr_xx_yzz[i] * tke_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_xxy_zz[i] * tbe_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 438-444 components of targeted buffer : FD

    auto tr_z_0_y_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 438);

    auto tr_z_0_y_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 439);

    auto tr_z_0_y_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 440);

    auto tr_z_0_y_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 441);

    auto tr_z_0_y_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 442);

    auto tr_z_0_y_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 443);

    #pragma omp simd aligned(tr_xyyyz_xx, tr_xyyyz_xy, tr_xyyyz_xz, tr_xyyyz_yy, tr_xyyyz_yz, tr_xyyyz_zz, tr_xyyz_x, tr_xyyz_xxy, tr_xyyz_xyy, tr_xyyz_xyz, tr_xyyz_y, tr_xyyz_yyy, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_z_0_y_xyy_xx, tr_z_0_y_xyy_xy, tr_z_0_y_xyy_xz, tr_z_0_y_xyy_yy, tr_z_0_y_xyy_yz, tr_z_0_y_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyy_xx[i] = -4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xy[i] = -4.0 * tr_xyz_xy[i] * tbe_0 - 2.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_xz[i] = -4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_yy[i] = -4.0 * tr_xyz_yy[i] * tbe_0 - 4.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_yz[i] = -4.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyy_zz[i] = -4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 444-450 components of targeted buffer : FD

    auto tr_z_0_y_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 444);

    auto tr_z_0_y_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 445);

    auto tr_z_0_y_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 446);

    auto tr_z_0_y_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 447);

    auto tr_z_0_y_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 448);

    auto tr_z_0_y_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 449);

    #pragma omp simd aligned(tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xy_x, tr_xy_xxy, tr_xy_xyy, tr_xy_xyz, tr_xy_y, tr_xy_yyy, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_xyzz_x, tr_xyzz_xxy, tr_xyzz_xyy, tr_xyzz_xyz, tr_xyzz_y, tr_xyzz_yyy, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_z_0_y_xyz_xx, tr_z_0_y_xyz_xy, tr_z_0_y_xyz_xz, tr_z_0_y_xyz_yy, tr_z_0_y_xyz_yz, tr_z_0_y_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xyz_xx[i] = tr_x_xx[i] - 2.0 * tr_xzz_xx[i] * tbe_0 - 2.0 * tr_xy_xxy[i] * tke_0 + 4.0 * tr_xyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xy[i] = tr_x_xy[i] - 2.0 * tr_xzz_xy[i] * tbe_0 + tr_xy_x[i] - 2.0 * tr_xy_xyy[i] * tke_0 - 2.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_xz[i] = tr_x_xz[i] - 2.0 * tr_xzz_xz[i] * tbe_0 - 2.0 * tr_xy_xyz[i] * tke_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_xz[i] * tbe_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_yy[i] = tr_x_yy[i] - 2.0 * tr_xzz_yy[i] * tbe_0 + 2.0 * tr_xy_y[i] - 2.0 * tr_xy_yyy[i] * tke_0 - 4.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_yz[i] = tr_x_yz[i] - 2.0 * tr_xzz_yz[i] * tbe_0 + tr_xy_z[i] - 2.0 * tr_xy_yyz[i] * tke_0 - 2.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_yz[i] * tbe_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xyz_zz[i] = tr_x_zz[i] - 2.0 * tr_xzz_zz[i] * tbe_0 - 2.0 * tr_xy_yzz[i] * tke_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_xyy_zz[i] * tbe_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 450-456 components of targeted buffer : FD

    auto tr_z_0_y_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 450);

    auto tr_z_0_y_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 451);

    auto tr_z_0_y_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 452);

    auto tr_z_0_y_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 453);

    auto tr_z_0_y_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 454);

    auto tr_z_0_y_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 455);

    #pragma omp simd aligned(tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_xz_x, tr_xz_xxy, tr_xz_xyy, tr_xz_xyz, tr_xz_y, tr_xz_yyy, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xzzz_x, tr_xzzz_xxy, tr_xzzz_xyy, tr_xzzz_xyz, tr_xzzz_y, tr_xzzz_yyy, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_z, tr_z_0_y_xzz_xx, tr_z_0_y_xzz_xy, tr_z_0_y_xzz_xz, tr_z_0_y_xzz_yy, tr_z_0_y_xzz_yz, tr_z_0_y_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_xzz_xx[i] = -4.0 * tr_xz_xxy[i] * tke_0 + 4.0 * tr_xzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xy[i] = 2.0 * tr_xz_x[i] - 4.0 * tr_xz_xyy[i] * tke_0 - 2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_xz[i] = -4.0 * tr_xz_xyz[i] * tke_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_xz[i] * tbe_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_yy[i] = 4.0 * tr_xz_y[i] - 4.0 * tr_xz_yyy[i] * tke_0 - 4.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yyy[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_yz[i] = 2.0 * tr_xz_z[i] - 4.0 * tr_xz_yyz[i] * tke_0 - 2.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzz_yyz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_yz[i] * tbe_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_xzz_zz[i] = -4.0 * tr_xz_yzz[i] * tke_0 + 4.0 * tr_xzzz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_xyz_zz[i] * tbe_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 456-462 components of targeted buffer : FD

    auto tr_z_0_y_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 456);

    auto tr_z_0_y_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 457);

    auto tr_z_0_y_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 458);

    auto tr_z_0_y_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 459);

    auto tr_z_0_y_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 460);

    auto tr_z_0_y_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 461);

    #pragma omp simd aligned(tr_yyyyz_xx, tr_yyyyz_xy, tr_yyyyz_xz, tr_yyyyz_yy, tr_yyyyz_yz, tr_yyyyz_zz, tr_yyyz_x, tr_yyyz_xxy, tr_yyyz_xyy, tr_yyyz_xyz, tr_yyyz_y, tr_yyyz_yyy, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_z, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_z_0_y_yyy_xx, tr_z_0_y_yyy_xy, tr_z_0_y_yyy_xz, tr_z_0_y_yyy_yy, tr_z_0_y_yyy_yz, tr_z_0_y_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyy_xx[i] = -6.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xy[i] = -6.0 * tr_yyz_xy[i] * tbe_0 - 2.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_xz[i] = -6.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_yy[i] = -6.0 * tr_yyz_yy[i] * tbe_0 - 4.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyz_yyy[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_yz[i] = -6.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyy_zz[i] = -6.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyyz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 462-468 components of targeted buffer : FD

    auto tr_z_0_y_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 462);

    auto tr_z_0_y_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 463);

    auto tr_z_0_y_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 464);

    auto tr_z_0_y_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 465);

    auto tr_z_0_y_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 466);

    auto tr_z_0_y_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 467);

    #pragma omp simd aligned(tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yy_x, tr_yy_xxy, tr_yy_xyy, tr_yy_xyz, tr_yy_y, tr_yy_yyy, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_yyzz_x, tr_yyzz_xxy, tr_yyzz_xyy, tr_yyzz_xyz, tr_yyzz_y, tr_yyzz_yyy, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_z, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_z_0_y_yyz_xx, tr_z_0_y_yyz_xy, tr_z_0_y_yyz_xz, tr_z_0_y_yyz_yy, tr_z_0_y_yyz_yz, tr_z_0_y_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yyz_xx[i] = 2.0 * tr_y_xx[i] - 4.0 * tr_yzz_xx[i] * tbe_0 - 2.0 * tr_yy_xxy[i] * tke_0 + 4.0 * tr_yyzz_xxy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_yyyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xy[i] = 2.0 * tr_y_xy[i] - 4.0 * tr_yzz_xy[i] * tbe_0 + tr_yy_x[i] - 2.0 * tr_yy_xyy[i] * tke_0 - 2.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xy[i] * tbe_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_xz[i] = 2.0 * tr_y_xz[i] - 4.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_yy_xyz[i] * tke_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_xz[i] * tbe_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_yy[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_yzz_yy[i] * tbe_0 + 2.0 * tr_yy_y[i] - 2.0 * tr_yy_yyy[i] * tke_0 - 4.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzz_yyy[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yy[i] * tbe_0 + 4.0 * tr_yyyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_yz[i] = 2.0 * tr_y_yz[i] - 4.0 * tr_yzz_yz[i] * tbe_0 + tr_yy_z[i] - 2.0 * tr_yy_yyz[i] * tke_0 - 2.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_yz[i] * tbe_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yyz_zz[i] = 2.0 * tr_y_zz[i] - 4.0 * tr_yzz_zz[i] * tbe_0 - 2.0 * tr_yy_yzz[i] * tke_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tke_0 - 2.0 * tr_yyy_zz[i] * tbe_0 + 4.0 * tr_yyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 468-474 components of targeted buffer : FD

    auto tr_z_0_y_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 468);

    auto tr_z_0_y_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 469);

    auto tr_z_0_y_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 470);

    auto tr_z_0_y_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 471);

    auto tr_z_0_y_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 472);

    auto tr_z_0_y_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 473);

    #pragma omp simd aligned(tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_yz_x, tr_yz_xxy, tr_yz_xyy, tr_yz_xyz, tr_yz_y, tr_yz_yyy, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yzzz_x, tr_yzzz_xxy, tr_yzzz_xyy, tr_yzzz_xyz, tr_yzzz_y, tr_yzzz_yyy, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_z, tr_z_0_y_yzz_xx, tr_z_0_y_yzz_xy, tr_z_0_y_yzz_xz, tr_z_0_y_yzz_yy, tr_z_0_y_yzz_yz, tr_z_0_y_yzz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_yzz_xx[i] = 2.0 * tr_z_xx[i] - 2.0 * tr_zzz_xx[i] * tbe_0 - 4.0 * tr_yz_xxy[i] * tke_0 + 4.0 * tr_yzzz_xxy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xy[i] = 2.0 * tr_z_xy[i] - 2.0 * tr_zzz_xy[i] * tbe_0 + 2.0 * tr_yz_x[i] - 4.0 * tr_yz_xyy[i] * tke_0 - 2.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_xz[i] = 2.0 * tr_z_xz[i] - 2.0 * tr_zzz_xz[i] * tbe_0 - 4.0 * tr_yz_xyz[i] * tke_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_xz[i] * tbe_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_yy[i] = 2.0 * tr_z_yy[i] - 2.0 * tr_zzz_yy[i] * tbe_0 + 4.0 * tr_yz_y[i] - 4.0 * tr_yz_yyy[i] * tke_0 - 4.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yyy[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_yz[i] = 2.0 * tr_z_yz[i] - 2.0 * tr_zzz_yz[i] * tbe_0 + 2.0 * tr_yz_z[i] - 4.0 * tr_yz_yyz[i] * tke_0 - 2.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzz_yyz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_yz[i] * tbe_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_yzz_zz[i] = 2.0 * tr_z_zz[i] - 2.0 * tr_zzz_zz[i] * tbe_0 - 4.0 * tr_yz_yzz[i] * tke_0 + 4.0 * tr_yzzz_yzz[i] * tbe_0 * tke_0 - 4.0 * tr_yyz_zz[i] * tbe_0 + 4.0 * tr_yyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 474-480 components of targeted buffer : FD

    auto tr_z_0_y_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 474);

    auto tr_z_0_y_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 475);

    auto tr_z_0_y_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 476);

    auto tr_z_0_y_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 477);

    auto tr_z_0_y_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 478);

    auto tr_z_0_y_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 479);

    #pragma omp simd aligned(tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, tr_z_0_y_zzz_xx, tr_z_0_y_zzz_xy, tr_z_0_y_zzz_xz, tr_z_0_y_zzz_yy, tr_z_0_y_zzz_yz, tr_z_0_y_zzz_zz, tr_zz_x, tr_zz_xxy, tr_zz_xyy, tr_zz_xyz, tr_zz_y, tr_zz_yyy, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zzzz_x, tr_zzzz_xxy, tr_zzzz_xyy, tr_zzzz_xyz, tr_zzzz_y, tr_zzzz_yyy, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_y_zzz_xx[i] = -6.0 * tr_zz_xxy[i] * tke_0 + 4.0 * tr_zzzz_xxy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xy[i] = 3.0 * tr_zz_x[i] - 6.0 * tr_zz_xyy[i] * tke_0 - 2.0 * tr_zzzz_x[i] * tbe_0 + 4.0 * tr_zzzz_xyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_xz[i] = -6.0 * tr_zz_xyz[i] * tke_0 + 4.0 * tr_zzzz_xyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_xz[i] * tbe_0 + 4.0 * tr_yzzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_yy[i] = 6.0 * tr_zz_y[i] - 6.0 * tr_zz_yyy[i] * tke_0 - 4.0 * tr_zzzz_y[i] * tbe_0 + 4.0 * tr_zzzz_yyy[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_yz[i] = 3.0 * tr_zz_z[i] - 6.0 * tr_zz_yyz[i] * tke_0 - 2.0 * tr_zzzz_z[i] * tbe_0 + 4.0 * tr_zzzz_yyz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_yz[i] * tbe_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_y_zzz_zz[i] = -6.0 * tr_zz_yzz[i] * tke_0 + 4.0 * tr_zzzz_yzz[i] * tbe_0 * tke_0 - 6.0 * tr_yzz_zz[i] * tbe_0 + 4.0 * tr_yzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 480-486 components of targeted buffer : FD

    auto tr_z_0_z_xxx_xx = pbuffer.data(idx_op_geom_110_fd + 480);

    auto tr_z_0_z_xxx_xy = pbuffer.data(idx_op_geom_110_fd + 481);

    auto tr_z_0_z_xxx_xz = pbuffer.data(idx_op_geom_110_fd + 482);

    auto tr_z_0_z_xxx_yy = pbuffer.data(idx_op_geom_110_fd + 483);

    auto tr_z_0_z_xxx_yz = pbuffer.data(idx_op_geom_110_fd + 484);

    auto tr_z_0_z_xxx_zz = pbuffer.data(idx_op_geom_110_fd + 485);

    #pragma omp simd aligned(tr_xxx_xx, tr_xxx_xy, tr_xxx_xz, tr_xxx_yy, tr_xxx_yz, tr_xxx_zz, tr_xxxz_x, tr_xxxz_xxz, tr_xxxz_xyz, tr_xxxz_xzz, tr_xxxz_y, tr_xxxz_yyz, tr_xxxz_yzz, tr_xxxz_z, tr_xxxz_zzz, tr_xxxzz_xx, tr_xxxzz_xy, tr_xxxzz_xz, tr_xxxzz_yy, tr_xxxzz_yz, tr_xxxzz_zz, tr_z_0_z_xxx_xx, tr_z_0_z_xxx_xy, tr_z_0_z_xxx_xz, tr_z_0_z_xxx_yy, tr_z_0_z_xxx_yz, tr_z_0_z_xxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxx_xx[i] = -2.0 * tr_xxx_xx[i] * tbe_0 + 4.0 * tr_xxxz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xy[i] = -2.0 * tr_xxx_xy[i] * tbe_0 + 4.0 * tr_xxxz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_xz[i] = -2.0 * tr_xxx_xz[i] * tbe_0 - 2.0 * tr_xxxz_x[i] * tbe_0 + 4.0 * tr_xxxz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_yy[i] = -2.0 * tr_xxx_yy[i] * tbe_0 + 4.0 * tr_xxxz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_yz[i] = -2.0 * tr_xxx_yz[i] * tbe_0 - 2.0 * tr_xxxz_y[i] * tbe_0 + 4.0 * tr_xxxz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxx_zz[i] = -2.0 * tr_xxx_zz[i] * tbe_0 - 4.0 * tr_xxxz_z[i] * tbe_0 + 4.0 * tr_xxxz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxxzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 486-492 components of targeted buffer : FD

    auto tr_z_0_z_xxy_xx = pbuffer.data(idx_op_geom_110_fd + 486);

    auto tr_z_0_z_xxy_xy = pbuffer.data(idx_op_geom_110_fd + 487);

    auto tr_z_0_z_xxy_xz = pbuffer.data(idx_op_geom_110_fd + 488);

    auto tr_z_0_z_xxy_yy = pbuffer.data(idx_op_geom_110_fd + 489);

    auto tr_z_0_z_xxy_yz = pbuffer.data(idx_op_geom_110_fd + 490);

    auto tr_z_0_z_xxy_zz = pbuffer.data(idx_op_geom_110_fd + 491);

    #pragma omp simd aligned(tr_xxy_xx, tr_xxy_xy, tr_xxy_xz, tr_xxy_yy, tr_xxy_yz, tr_xxy_zz, tr_xxyz_x, tr_xxyz_xxz, tr_xxyz_xyz, tr_xxyz_xzz, tr_xxyz_y, tr_xxyz_yyz, tr_xxyz_yzz, tr_xxyz_z, tr_xxyz_zzz, tr_xxyzz_xx, tr_xxyzz_xy, tr_xxyzz_xz, tr_xxyzz_yy, tr_xxyzz_yz, tr_xxyzz_zz, tr_z_0_z_xxy_xx, tr_z_0_z_xxy_xy, tr_z_0_z_xxy_xz, tr_z_0_z_xxy_yy, tr_z_0_z_xxy_yz, tr_z_0_z_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxy_xx[i] = -2.0 * tr_xxy_xx[i] * tbe_0 + 4.0 * tr_xxyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xy[i] = -2.0 * tr_xxy_xy[i] * tbe_0 + 4.0 * tr_xxyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_xz[i] = -2.0 * tr_xxy_xz[i] * tbe_0 - 2.0 * tr_xxyz_x[i] * tbe_0 + 4.0 * tr_xxyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_yy[i] = -2.0 * tr_xxy_yy[i] * tbe_0 + 4.0 * tr_xxyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_yz[i] = -2.0 * tr_xxy_yz[i] * tbe_0 - 2.0 * tr_xxyz_y[i] * tbe_0 + 4.0 * tr_xxyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxy_zz[i] = -2.0 * tr_xxy_zz[i] * tbe_0 - 4.0 * tr_xxyz_z[i] * tbe_0 + 4.0 * tr_xxyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 492-498 components of targeted buffer : FD

    auto tr_z_0_z_xxz_xx = pbuffer.data(idx_op_geom_110_fd + 492);

    auto tr_z_0_z_xxz_xy = pbuffer.data(idx_op_geom_110_fd + 493);

    auto tr_z_0_z_xxz_xz = pbuffer.data(idx_op_geom_110_fd + 494);

    auto tr_z_0_z_xxz_yy = pbuffer.data(idx_op_geom_110_fd + 495);

    auto tr_z_0_z_xxz_yz = pbuffer.data(idx_op_geom_110_fd + 496);

    auto tr_z_0_z_xxz_zz = pbuffer.data(idx_op_geom_110_fd + 497);

    #pragma omp simd aligned(tr_xx_x, tr_xx_xxz, tr_xx_xyz, tr_xx_xzz, tr_xx_y, tr_xx_yyz, tr_xx_yzz, tr_xx_z, tr_xx_zzz, tr_xxz_xx, tr_xxz_xy, tr_xxz_xz, tr_xxz_yy, tr_xxz_yz, tr_xxz_zz, tr_xxzz_x, tr_xxzz_xxz, tr_xxzz_xyz, tr_xxzz_xzz, tr_xxzz_y, tr_xxzz_yyz, tr_xxzz_yzz, tr_xxzz_z, tr_xxzz_zzz, tr_xxzzz_xx, tr_xxzzz_xy, tr_xxzzz_xz, tr_xxzzz_yy, tr_xxzzz_yz, tr_xxzzz_zz, tr_z_0_z_xxz_xx, tr_z_0_z_xxz_xy, tr_z_0_z_xxz_xz, tr_z_0_z_xxz_yy, tr_z_0_z_xxz_yz, tr_z_0_z_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xxz_xx[i] = -2.0 * tr_xx_xxz[i] * tke_0 - 6.0 * tr_xxz_xx[i] * tbe_0 + 4.0 * tr_xxzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xy[i] = -2.0 * tr_xx_xyz[i] * tke_0 - 6.0 * tr_xxz_xy[i] * tbe_0 + 4.0 * tr_xxzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_xz[i] = tr_xx_x[i] - 2.0 * tr_xx_xzz[i] * tke_0 - 6.0 * tr_xxz_xz[i] * tbe_0 - 2.0 * tr_xxzz_x[i] * tbe_0 + 4.0 * tr_xxzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_yy[i] = -2.0 * tr_xx_yyz[i] * tke_0 - 6.0 * tr_xxz_yy[i] * tbe_0 + 4.0 * tr_xxzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_yz[i] = tr_xx_y[i] - 2.0 * tr_xx_yzz[i] * tke_0 - 6.0 * tr_xxz_yz[i] * tbe_0 - 2.0 * tr_xxzz_y[i] * tbe_0 + 4.0 * tr_xxzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xxz_zz[i] = 2.0 * tr_xx_z[i] - 2.0 * tr_xx_zzz[i] * tke_0 - 6.0 * tr_xxz_zz[i] * tbe_0 - 4.0 * tr_xxzz_z[i] * tbe_0 + 4.0 * tr_xxzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xxzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 498-504 components of targeted buffer : FD

    auto tr_z_0_z_xyy_xx = pbuffer.data(idx_op_geom_110_fd + 498);

    auto tr_z_0_z_xyy_xy = pbuffer.data(idx_op_geom_110_fd + 499);

    auto tr_z_0_z_xyy_xz = pbuffer.data(idx_op_geom_110_fd + 500);

    auto tr_z_0_z_xyy_yy = pbuffer.data(idx_op_geom_110_fd + 501);

    auto tr_z_0_z_xyy_yz = pbuffer.data(idx_op_geom_110_fd + 502);

    auto tr_z_0_z_xyy_zz = pbuffer.data(idx_op_geom_110_fd + 503);

    #pragma omp simd aligned(tr_xyy_xx, tr_xyy_xy, tr_xyy_xz, tr_xyy_yy, tr_xyy_yz, tr_xyy_zz, tr_xyyz_x, tr_xyyz_xxz, tr_xyyz_xyz, tr_xyyz_xzz, tr_xyyz_y, tr_xyyz_yyz, tr_xyyz_yzz, tr_xyyz_z, tr_xyyz_zzz, tr_xyyzz_xx, tr_xyyzz_xy, tr_xyyzz_xz, tr_xyyzz_yy, tr_xyyzz_yz, tr_xyyzz_zz, tr_z_0_z_xyy_xx, tr_z_0_z_xyy_xy, tr_z_0_z_xyy_xz, tr_z_0_z_xyy_yy, tr_z_0_z_xyy_yz, tr_z_0_z_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyy_xx[i] = -2.0 * tr_xyy_xx[i] * tbe_0 + 4.0 * tr_xyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xy[i] = -2.0 * tr_xyy_xy[i] * tbe_0 + 4.0 * tr_xyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_xz[i] = -2.0 * tr_xyy_xz[i] * tbe_0 - 2.0 * tr_xyyz_x[i] * tbe_0 + 4.0 * tr_xyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_yy[i] = -2.0 * tr_xyy_yy[i] * tbe_0 + 4.0 * tr_xyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_yz[i] = -2.0 * tr_xyy_yz[i] * tbe_0 - 2.0 * tr_xyyz_y[i] * tbe_0 + 4.0 * tr_xyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyy_zz[i] = -2.0 * tr_xyy_zz[i] * tbe_0 - 4.0 * tr_xyyz_z[i] * tbe_0 + 4.0 * tr_xyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 504-510 components of targeted buffer : FD

    auto tr_z_0_z_xyz_xx = pbuffer.data(idx_op_geom_110_fd + 504);

    auto tr_z_0_z_xyz_xy = pbuffer.data(idx_op_geom_110_fd + 505);

    auto tr_z_0_z_xyz_xz = pbuffer.data(idx_op_geom_110_fd + 506);

    auto tr_z_0_z_xyz_yy = pbuffer.data(idx_op_geom_110_fd + 507);

    auto tr_z_0_z_xyz_yz = pbuffer.data(idx_op_geom_110_fd + 508);

    auto tr_z_0_z_xyz_zz = pbuffer.data(idx_op_geom_110_fd + 509);

    #pragma omp simd aligned(tr_xy_x, tr_xy_xxz, tr_xy_xyz, tr_xy_xzz, tr_xy_y, tr_xy_yyz, tr_xy_yzz, tr_xy_z, tr_xy_zzz, tr_xyz_xx, tr_xyz_xy, tr_xyz_xz, tr_xyz_yy, tr_xyz_yz, tr_xyz_zz, tr_xyzz_x, tr_xyzz_xxz, tr_xyzz_xyz, tr_xyzz_xzz, tr_xyzz_y, tr_xyzz_yyz, tr_xyzz_yzz, tr_xyzz_z, tr_xyzz_zzz, tr_xyzzz_xx, tr_xyzzz_xy, tr_xyzzz_xz, tr_xyzzz_yy, tr_xyzzz_yz, tr_xyzzz_zz, tr_z_0_z_xyz_xx, tr_z_0_z_xyz_xy, tr_z_0_z_xyz_xz, tr_z_0_z_xyz_yy, tr_z_0_z_xyz_yz, tr_z_0_z_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xyz_xx[i] = -2.0 * tr_xy_xxz[i] * tke_0 - 6.0 * tr_xyz_xx[i] * tbe_0 + 4.0 * tr_xyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xy[i] = -2.0 * tr_xy_xyz[i] * tke_0 - 6.0 * tr_xyz_xy[i] * tbe_0 + 4.0 * tr_xyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_xz[i] = tr_xy_x[i] - 2.0 * tr_xy_xzz[i] * tke_0 - 6.0 * tr_xyz_xz[i] * tbe_0 - 2.0 * tr_xyzz_x[i] * tbe_0 + 4.0 * tr_xyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_yy[i] = -2.0 * tr_xy_yyz[i] * tke_0 - 6.0 * tr_xyz_yy[i] * tbe_0 + 4.0 * tr_xyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_yz[i] = tr_xy_y[i] - 2.0 * tr_xy_yzz[i] * tke_0 - 6.0 * tr_xyz_yz[i] * tbe_0 - 2.0 * tr_xyzz_y[i] * tbe_0 + 4.0 * tr_xyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xyz_zz[i] = 2.0 * tr_xy_z[i] - 2.0 * tr_xy_zzz[i] * tke_0 - 6.0 * tr_xyz_zz[i] * tbe_0 - 4.0 * tr_xyzz_z[i] * tbe_0 + 4.0 * tr_xyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 510-516 components of targeted buffer : FD

    auto tr_z_0_z_xzz_xx = pbuffer.data(idx_op_geom_110_fd + 510);

    auto tr_z_0_z_xzz_xy = pbuffer.data(idx_op_geom_110_fd + 511);

    auto tr_z_0_z_xzz_xz = pbuffer.data(idx_op_geom_110_fd + 512);

    auto tr_z_0_z_xzz_yy = pbuffer.data(idx_op_geom_110_fd + 513);

    auto tr_z_0_z_xzz_yz = pbuffer.data(idx_op_geom_110_fd + 514);

    auto tr_z_0_z_xzz_zz = pbuffer.data(idx_op_geom_110_fd + 515);

    #pragma omp simd aligned(tr_x_xx, tr_x_xy, tr_x_xz, tr_x_yy, tr_x_yz, tr_x_zz, tr_xz_x, tr_xz_xxz, tr_xz_xyz, tr_xz_xzz, tr_xz_y, tr_xz_yyz, tr_xz_yzz, tr_xz_z, tr_xz_zzz, tr_xzz_xx, tr_xzz_xy, tr_xzz_xz, tr_xzz_yy, tr_xzz_yz, tr_xzz_zz, tr_xzzz_x, tr_xzzz_xxz, tr_xzzz_xyz, tr_xzzz_xzz, tr_xzzz_y, tr_xzzz_yyz, tr_xzzz_yzz, tr_xzzz_z, tr_xzzz_zzz, tr_xzzzz_xx, tr_xzzzz_xy, tr_xzzzz_xz, tr_xzzzz_yy, tr_xzzzz_yz, tr_xzzzz_zz, tr_z_0_z_xzz_xx, tr_z_0_z_xzz_xy, tr_z_0_z_xzz_xz, tr_z_0_z_xzz_yy, tr_z_0_z_xzz_yz, tr_z_0_z_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_xzz_xx[i] = 2.0 * tr_x_xx[i] - 4.0 * tr_xz_xxz[i] * tke_0 - 10.0 * tr_xzz_xx[i] * tbe_0 + 4.0 * tr_xzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xy[i] = 2.0 * tr_x_xy[i] - 4.0 * tr_xz_xyz[i] * tke_0 - 10.0 * tr_xzz_xy[i] * tbe_0 + 4.0 * tr_xzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_xz[i] = 2.0 * tr_x_xz[i] + 2.0 * tr_xz_x[i] - 4.0 * tr_xz_xzz[i] * tke_0 - 10.0 * tr_xzz_xz[i] * tbe_0 - 2.0 * tr_xzzz_x[i] * tbe_0 + 4.0 * tr_xzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_yy[i] = 2.0 * tr_x_yy[i] - 4.0 * tr_xz_yyz[i] * tke_0 - 10.0 * tr_xzz_yy[i] * tbe_0 + 4.0 * tr_xzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_yz[i] = 2.0 * tr_x_yz[i] + 2.0 * tr_xz_y[i] - 4.0 * tr_xz_yzz[i] * tke_0 - 10.0 * tr_xzz_yz[i] * tbe_0 - 2.0 * tr_xzzz_y[i] * tbe_0 + 4.0 * tr_xzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_xzz_zz[i] = 2.0 * tr_x_zz[i] + 4.0 * tr_xz_z[i] - 4.0 * tr_xz_zzz[i] * tke_0 - 10.0 * tr_xzz_zz[i] * tbe_0 - 4.0 * tr_xzzz_z[i] * tbe_0 + 4.0 * tr_xzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_xzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 516-522 components of targeted buffer : FD

    auto tr_z_0_z_yyy_xx = pbuffer.data(idx_op_geom_110_fd + 516);

    auto tr_z_0_z_yyy_xy = pbuffer.data(idx_op_geom_110_fd + 517);

    auto tr_z_0_z_yyy_xz = pbuffer.data(idx_op_geom_110_fd + 518);

    auto tr_z_0_z_yyy_yy = pbuffer.data(idx_op_geom_110_fd + 519);

    auto tr_z_0_z_yyy_yz = pbuffer.data(idx_op_geom_110_fd + 520);

    auto tr_z_0_z_yyy_zz = pbuffer.data(idx_op_geom_110_fd + 521);

    #pragma omp simd aligned(tr_yyy_xx, tr_yyy_xy, tr_yyy_xz, tr_yyy_yy, tr_yyy_yz, tr_yyy_zz, tr_yyyz_x, tr_yyyz_xxz, tr_yyyz_xyz, tr_yyyz_xzz, tr_yyyz_y, tr_yyyz_yyz, tr_yyyz_yzz, tr_yyyz_z, tr_yyyz_zzz, tr_yyyzz_xx, tr_yyyzz_xy, tr_yyyzz_xz, tr_yyyzz_yy, tr_yyyzz_yz, tr_yyyzz_zz, tr_z_0_z_yyy_xx, tr_z_0_z_yyy_xy, tr_z_0_z_yyy_xz, tr_z_0_z_yyy_yy, tr_z_0_z_yyy_yz, tr_z_0_z_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyy_xx[i] = -2.0 * tr_yyy_xx[i] * tbe_0 + 4.0 * tr_yyyz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xy[i] = -2.0 * tr_yyy_xy[i] * tbe_0 + 4.0 * tr_yyyz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_xz[i] = -2.0 * tr_yyy_xz[i] * tbe_0 - 2.0 * tr_yyyz_x[i] * tbe_0 + 4.0 * tr_yyyz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_yy[i] = -2.0 * tr_yyy_yy[i] * tbe_0 + 4.0 * tr_yyyz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_yz[i] = -2.0 * tr_yyy_yz[i] * tbe_0 - 2.0 * tr_yyyz_y[i] * tbe_0 + 4.0 * tr_yyyz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyy_zz[i] = -2.0 * tr_yyy_zz[i] * tbe_0 - 4.0 * tr_yyyz_z[i] * tbe_0 + 4.0 * tr_yyyz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyyzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 522-528 components of targeted buffer : FD

    auto tr_z_0_z_yyz_xx = pbuffer.data(idx_op_geom_110_fd + 522);

    auto tr_z_0_z_yyz_xy = pbuffer.data(idx_op_geom_110_fd + 523);

    auto tr_z_0_z_yyz_xz = pbuffer.data(idx_op_geom_110_fd + 524);

    auto tr_z_0_z_yyz_yy = pbuffer.data(idx_op_geom_110_fd + 525);

    auto tr_z_0_z_yyz_yz = pbuffer.data(idx_op_geom_110_fd + 526);

    auto tr_z_0_z_yyz_zz = pbuffer.data(idx_op_geom_110_fd + 527);

    #pragma omp simd aligned(tr_yy_x, tr_yy_xxz, tr_yy_xyz, tr_yy_xzz, tr_yy_y, tr_yy_yyz, tr_yy_yzz, tr_yy_z, tr_yy_zzz, tr_yyz_xx, tr_yyz_xy, tr_yyz_xz, tr_yyz_yy, tr_yyz_yz, tr_yyz_zz, tr_yyzz_x, tr_yyzz_xxz, tr_yyzz_xyz, tr_yyzz_xzz, tr_yyzz_y, tr_yyzz_yyz, tr_yyzz_yzz, tr_yyzz_z, tr_yyzz_zzz, tr_yyzzz_xx, tr_yyzzz_xy, tr_yyzzz_xz, tr_yyzzz_yy, tr_yyzzz_yz, tr_yyzzz_zz, tr_z_0_z_yyz_xx, tr_z_0_z_yyz_xy, tr_z_0_z_yyz_xz, tr_z_0_z_yyz_yy, tr_z_0_z_yyz_yz, tr_z_0_z_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yyz_xx[i] = -2.0 * tr_yy_xxz[i] * tke_0 - 6.0 * tr_yyz_xx[i] * tbe_0 + 4.0 * tr_yyzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xy[i] = -2.0 * tr_yy_xyz[i] * tke_0 - 6.0 * tr_yyz_xy[i] * tbe_0 + 4.0 * tr_yyzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_xz[i] = tr_yy_x[i] - 2.0 * tr_yy_xzz[i] * tke_0 - 6.0 * tr_yyz_xz[i] * tbe_0 - 2.0 * tr_yyzz_x[i] * tbe_0 + 4.0 * tr_yyzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_yy[i] = -2.0 * tr_yy_yyz[i] * tke_0 - 6.0 * tr_yyz_yy[i] * tbe_0 + 4.0 * tr_yyzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_yz[i] = tr_yy_y[i] - 2.0 * tr_yy_yzz[i] * tke_0 - 6.0 * tr_yyz_yz[i] * tbe_0 - 2.0 * tr_yyzz_y[i] * tbe_0 + 4.0 * tr_yyzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yyz_zz[i] = 2.0 * tr_yy_z[i] - 2.0 * tr_yy_zzz[i] * tke_0 - 6.0 * tr_yyz_zz[i] * tbe_0 - 4.0 * tr_yyzz_z[i] * tbe_0 + 4.0 * tr_yyzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yyzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 528-534 components of targeted buffer : FD

    auto tr_z_0_z_yzz_xx = pbuffer.data(idx_op_geom_110_fd + 528);

    auto tr_z_0_z_yzz_xy = pbuffer.data(idx_op_geom_110_fd + 529);

    auto tr_z_0_z_yzz_xz = pbuffer.data(idx_op_geom_110_fd + 530);

    auto tr_z_0_z_yzz_yy = pbuffer.data(idx_op_geom_110_fd + 531);

    auto tr_z_0_z_yzz_yz = pbuffer.data(idx_op_geom_110_fd + 532);

    auto tr_z_0_z_yzz_zz = pbuffer.data(idx_op_geom_110_fd + 533);

    #pragma omp simd aligned(tr_y_xx, tr_y_xy, tr_y_xz, tr_y_yy, tr_y_yz, tr_y_zz, tr_yz_x, tr_yz_xxz, tr_yz_xyz, tr_yz_xzz, tr_yz_y, tr_yz_yyz, tr_yz_yzz, tr_yz_z, tr_yz_zzz, tr_yzz_xx, tr_yzz_xy, tr_yzz_xz, tr_yzz_yy, tr_yzz_yz, tr_yzz_zz, tr_yzzz_x, tr_yzzz_xxz, tr_yzzz_xyz, tr_yzzz_xzz, tr_yzzz_y, tr_yzzz_yyz, tr_yzzz_yzz, tr_yzzz_z, tr_yzzz_zzz, tr_yzzzz_xx, tr_yzzzz_xy, tr_yzzzz_xz, tr_yzzzz_yy, tr_yzzzz_yz, tr_yzzzz_zz, tr_z_0_z_yzz_xx, tr_z_0_z_yzz_xy, tr_z_0_z_yzz_xz, tr_z_0_z_yzz_yy, tr_z_0_z_yzz_yz, tr_z_0_z_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_yzz_xx[i] = 2.0 * tr_y_xx[i] - 4.0 * tr_yz_xxz[i] * tke_0 - 10.0 * tr_yzz_xx[i] * tbe_0 + 4.0 * tr_yzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xy[i] = 2.0 * tr_y_xy[i] - 4.0 * tr_yz_xyz[i] * tke_0 - 10.0 * tr_yzz_xy[i] * tbe_0 + 4.0 * tr_yzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_xz[i] = 2.0 * tr_y_xz[i] + 2.0 * tr_yz_x[i] - 4.0 * tr_yz_xzz[i] * tke_0 - 10.0 * tr_yzz_xz[i] * tbe_0 - 2.0 * tr_yzzz_x[i] * tbe_0 + 4.0 * tr_yzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_yy[i] = 2.0 * tr_y_yy[i] - 4.0 * tr_yz_yyz[i] * tke_0 - 10.0 * tr_yzz_yy[i] * tbe_0 + 4.0 * tr_yzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_yz[i] = 2.0 * tr_y_yz[i] + 2.0 * tr_yz_y[i] - 4.0 * tr_yz_yzz[i] * tke_0 - 10.0 * tr_yzz_yz[i] * tbe_0 - 2.0 * tr_yzzz_y[i] * tbe_0 + 4.0 * tr_yzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_yzz_zz[i] = 2.0 * tr_y_zz[i] + 4.0 * tr_yz_z[i] - 4.0 * tr_yz_zzz[i] * tke_0 - 10.0 * tr_yzz_zz[i] * tbe_0 - 4.0 * tr_yzzz_z[i] * tbe_0 + 4.0 * tr_yzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_yzzzz_zz[i] * tbe_0 * tbe_0;
    }

    // Set up 534-540 components of targeted buffer : FD

    auto tr_z_0_z_zzz_xx = pbuffer.data(idx_op_geom_110_fd + 534);

    auto tr_z_0_z_zzz_xy = pbuffer.data(idx_op_geom_110_fd + 535);

    auto tr_z_0_z_zzz_xz = pbuffer.data(idx_op_geom_110_fd + 536);

    auto tr_z_0_z_zzz_yy = pbuffer.data(idx_op_geom_110_fd + 537);

    auto tr_z_0_z_zzz_yz = pbuffer.data(idx_op_geom_110_fd + 538);

    auto tr_z_0_z_zzz_zz = pbuffer.data(idx_op_geom_110_fd + 539);

    #pragma omp simd aligned(tr_z_0_z_zzz_xx, tr_z_0_z_zzz_xy, tr_z_0_z_zzz_xz, tr_z_0_z_zzz_yy, tr_z_0_z_zzz_yz, tr_z_0_z_zzz_zz, tr_z_xx, tr_z_xy, tr_z_xz, tr_z_yy, tr_z_yz, tr_z_zz, tr_zz_x, tr_zz_xxz, tr_zz_xyz, tr_zz_xzz, tr_zz_y, tr_zz_yyz, tr_zz_yzz, tr_zz_z, tr_zz_zzz, tr_zzz_xx, tr_zzz_xy, tr_zzz_xz, tr_zzz_yy, tr_zzz_yz, tr_zzz_zz, tr_zzzz_x, tr_zzzz_xxz, tr_zzzz_xyz, tr_zzzz_xzz, tr_zzzz_y, tr_zzzz_yyz, tr_zzzz_yzz, tr_zzzz_z, tr_zzzz_zzz, tr_zzzzz_xx, tr_zzzzz_xy, tr_zzzzz_xz, tr_zzzzz_yy, tr_zzzzz_yz, tr_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tbe_0 = a_exp;

        const double tke_0 = b_exps[i];

        tr_z_0_z_zzz_xx[i] = 6.0 * tr_z_xx[i] - 6.0 * tr_zz_xxz[i] * tke_0 - 14.0 * tr_zzz_xx[i] * tbe_0 + 4.0 * tr_zzzz_xxz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xx[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xy[i] = 6.0 * tr_z_xy[i] - 6.0 * tr_zz_xyz[i] * tke_0 - 14.0 * tr_zzz_xy[i] * tbe_0 + 4.0 * tr_zzzz_xyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_xz[i] = 6.0 * tr_z_xz[i] + 3.0 * tr_zz_x[i] - 6.0 * tr_zz_xzz[i] * tke_0 - 14.0 * tr_zzz_xz[i] * tbe_0 - 2.0 * tr_zzzz_x[i] * tbe_0 + 4.0 * tr_zzzz_xzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_xz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_yy[i] = 6.0 * tr_z_yy[i] - 6.0 * tr_zz_yyz[i] * tke_0 - 14.0 * tr_zzz_yy[i] * tbe_0 + 4.0 * tr_zzzz_yyz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yy[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_yz[i] = 6.0 * tr_z_yz[i] + 3.0 * tr_zz_y[i] - 6.0 * tr_zz_yzz[i] * tke_0 - 14.0 * tr_zzz_yz[i] * tbe_0 - 2.0 * tr_zzzz_y[i] * tbe_0 + 4.0 * tr_zzzz_yzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_yz[i] * tbe_0 * tbe_0;

        tr_z_0_z_zzz_zz[i] = 6.0 * tr_z_zz[i] + 6.0 * tr_zz_z[i] - 6.0 * tr_zz_zzz[i] * tke_0 - 14.0 * tr_zzz_zz[i] * tbe_0 - 4.0 * tr_zzzz_z[i] * tbe_0 + 4.0 * tr_zzzz_zzz[i] * tbe_0 * tke_0 + 4.0 * tr_zzzzz_zz[i] * tbe_0 * tbe_0;
    }

}

} // t2cgeom namespace

