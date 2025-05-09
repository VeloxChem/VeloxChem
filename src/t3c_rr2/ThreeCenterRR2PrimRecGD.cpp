#include "ThreeCenterRR2PrimRecGD.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_gd(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_gd,
                  const size_t idx_fd,
                  const size_t idx_g_fd,
                  const size_t idx_gp,
                  const size_t idx_g_gp,
                  const size_t idx_gd,
                  const size_t idx_g_gd,
                  const CSimdArray<double>& factors,
                  const size_t idx_rgc,
                  const double a_exp,
                  const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_fd);

    auto ts_xxx_xy = pbuffer.data(idx_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_fd + 10);

    auto ts_xxy_zz = pbuffer.data(idx_fd + 11);

    auto ts_xxz_xx = pbuffer.data(idx_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_fd + 14);

    auto ts_xxz_yy = pbuffer.data(idx_fd + 15);

    auto ts_xxz_yz = pbuffer.data(idx_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_fd + 17);

    auto ts_xyy_xx = pbuffer.data(idx_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_fd + 19);

    auto ts_xyy_xz = pbuffer.data(idx_fd + 20);

    auto ts_xyy_yy = pbuffer.data(idx_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_fd + 23);

    auto ts_xyz_xx = pbuffer.data(idx_fd + 24);

    auto ts_xyz_xy = pbuffer.data(idx_fd + 25);

    auto ts_xyz_xz = pbuffer.data(idx_fd + 26);

    auto ts_xyz_yy = pbuffer.data(idx_fd + 27);

    auto ts_xyz_yz = pbuffer.data(idx_fd + 28);

    auto ts_xyz_zz = pbuffer.data(idx_fd + 29);

    auto ts_xzz_xx = pbuffer.data(idx_fd + 30);

    auto ts_xzz_xy = pbuffer.data(idx_fd + 31);

    auto ts_xzz_xz = pbuffer.data(idx_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_fd + 41);

    auto ts_yyz_xx = pbuffer.data(idx_fd + 42);

    auto ts_yyz_xy = pbuffer.data(idx_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_fd + 47);

    auto ts_yzz_xx = pbuffer.data(idx_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto gr_xxx_xx = pbuffer.data(idx_g_fd);

    auto gr_xxx_xy = pbuffer.data(idx_g_fd + 1);

    auto gr_xxx_xz = pbuffer.data(idx_g_fd + 2);

    auto gr_xxx_yy = pbuffer.data(idx_g_fd + 3);

    auto gr_xxx_yz = pbuffer.data(idx_g_fd + 4);

    auto gr_xxx_zz = pbuffer.data(idx_g_fd + 5);

    auto gr_xxy_xx = pbuffer.data(idx_g_fd + 6);

    auto gr_xxy_xy = pbuffer.data(idx_g_fd + 7);

    auto gr_xxy_xz = pbuffer.data(idx_g_fd + 8);

    auto gr_xxy_yy = pbuffer.data(idx_g_fd + 9);

    auto gr_xxy_yz = pbuffer.data(idx_g_fd + 10);

    auto gr_xxy_zz = pbuffer.data(idx_g_fd + 11);

    auto gr_xxz_xx = pbuffer.data(idx_g_fd + 12);

    auto gr_xxz_xy = pbuffer.data(idx_g_fd + 13);

    auto gr_xxz_xz = pbuffer.data(idx_g_fd + 14);

    auto gr_xxz_yy = pbuffer.data(idx_g_fd + 15);

    auto gr_xxz_yz = pbuffer.data(idx_g_fd + 16);

    auto gr_xxz_zz = pbuffer.data(idx_g_fd + 17);

    auto gr_xyy_xx = pbuffer.data(idx_g_fd + 18);

    auto gr_xyy_xy = pbuffer.data(idx_g_fd + 19);

    auto gr_xyy_xz = pbuffer.data(idx_g_fd + 20);

    auto gr_xyy_yy = pbuffer.data(idx_g_fd + 21);

    auto gr_xyy_yz = pbuffer.data(idx_g_fd + 22);

    auto gr_xyy_zz = pbuffer.data(idx_g_fd + 23);

    auto gr_xyz_xx = pbuffer.data(idx_g_fd + 24);

    auto gr_xyz_xy = pbuffer.data(idx_g_fd + 25);

    auto gr_xyz_xz = pbuffer.data(idx_g_fd + 26);

    auto gr_xyz_yy = pbuffer.data(idx_g_fd + 27);

    auto gr_xyz_yz = pbuffer.data(idx_g_fd + 28);

    auto gr_xyz_zz = pbuffer.data(idx_g_fd + 29);

    auto gr_xzz_xx = pbuffer.data(idx_g_fd + 30);

    auto gr_xzz_xy = pbuffer.data(idx_g_fd + 31);

    auto gr_xzz_xz = pbuffer.data(idx_g_fd + 32);

    auto gr_xzz_yy = pbuffer.data(idx_g_fd + 33);

    auto gr_xzz_yz = pbuffer.data(idx_g_fd + 34);

    auto gr_xzz_zz = pbuffer.data(idx_g_fd + 35);

    auto gr_yyy_xx = pbuffer.data(idx_g_fd + 36);

    auto gr_yyy_xy = pbuffer.data(idx_g_fd + 37);

    auto gr_yyy_xz = pbuffer.data(idx_g_fd + 38);

    auto gr_yyy_yy = pbuffer.data(idx_g_fd + 39);

    auto gr_yyy_yz = pbuffer.data(idx_g_fd + 40);

    auto gr_yyy_zz = pbuffer.data(idx_g_fd + 41);

    auto gr_yyz_xx = pbuffer.data(idx_g_fd + 42);

    auto gr_yyz_xy = pbuffer.data(idx_g_fd + 43);

    auto gr_yyz_xz = pbuffer.data(idx_g_fd + 44);

    auto gr_yyz_yy = pbuffer.data(idx_g_fd + 45);

    auto gr_yyz_yz = pbuffer.data(idx_g_fd + 46);

    auto gr_yyz_zz = pbuffer.data(idx_g_fd + 47);

    auto gr_yzz_xx = pbuffer.data(idx_g_fd + 48);

    auto gr_yzz_xy = pbuffer.data(idx_g_fd + 49);

    auto gr_yzz_xz = pbuffer.data(idx_g_fd + 50);

    auto gr_yzz_yy = pbuffer.data(idx_g_fd + 51);

    auto gr_yzz_yz = pbuffer.data(idx_g_fd + 52);

    auto gr_yzz_zz = pbuffer.data(idx_g_fd + 53);

    auto gr_zzz_xx = pbuffer.data(idx_g_fd + 54);

    auto gr_zzz_xy = pbuffer.data(idx_g_fd + 55);

    auto gr_zzz_xz = pbuffer.data(idx_g_fd + 56);

    auto gr_zzz_yy = pbuffer.data(idx_g_fd + 57);

    auto gr_zzz_yz = pbuffer.data(idx_g_fd + 58);

    auto gr_zzz_zz = pbuffer.data(idx_g_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_gp + 5);

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    auto ts_xxyz_x = pbuffer.data(idx_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_gp + 14);

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    auto ts_xyyz_x = pbuffer.data(idx_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    auto ts_xyzz_x = pbuffer.data(idx_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_gp + 26);

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

    auto ts_yyyz_x = pbuffer.data(idx_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto gr_xxxx_x = pbuffer.data(idx_g_gp);

    auto gr_xxxx_y = pbuffer.data(idx_g_gp + 1);

    auto gr_xxxx_z = pbuffer.data(idx_g_gp + 2);

    auto gr_xxxy_x = pbuffer.data(idx_g_gp + 3);

    auto gr_xxxy_y = pbuffer.data(idx_g_gp + 4);

    auto gr_xxxy_z = pbuffer.data(idx_g_gp + 5);

    auto gr_xxxz_x = pbuffer.data(idx_g_gp + 6);

    auto gr_xxxz_y = pbuffer.data(idx_g_gp + 7);

    auto gr_xxxz_z = pbuffer.data(idx_g_gp + 8);

    auto gr_xxyy_x = pbuffer.data(idx_g_gp + 9);

    auto gr_xxyy_y = pbuffer.data(idx_g_gp + 10);

    auto gr_xxyy_z = pbuffer.data(idx_g_gp + 11);

    auto gr_xxyz_x = pbuffer.data(idx_g_gp + 12);

    auto gr_xxyz_y = pbuffer.data(idx_g_gp + 13);

    auto gr_xxyz_z = pbuffer.data(idx_g_gp + 14);

    auto gr_xxzz_x = pbuffer.data(idx_g_gp + 15);

    auto gr_xxzz_y = pbuffer.data(idx_g_gp + 16);

    auto gr_xxzz_z = pbuffer.data(idx_g_gp + 17);

    auto gr_xyyy_x = pbuffer.data(idx_g_gp + 18);

    auto gr_xyyy_y = pbuffer.data(idx_g_gp + 19);

    auto gr_xyyy_z = pbuffer.data(idx_g_gp + 20);

    auto gr_xyyz_x = pbuffer.data(idx_g_gp + 21);

    auto gr_xyyz_y = pbuffer.data(idx_g_gp + 22);

    auto gr_xyyz_z = pbuffer.data(idx_g_gp + 23);

    auto gr_xyzz_x = pbuffer.data(idx_g_gp + 24);

    auto gr_xyzz_y = pbuffer.data(idx_g_gp + 25);

    auto gr_xyzz_z = pbuffer.data(idx_g_gp + 26);

    auto gr_xzzz_x = pbuffer.data(idx_g_gp + 27);

    auto gr_xzzz_y = pbuffer.data(idx_g_gp + 28);

    auto gr_xzzz_z = pbuffer.data(idx_g_gp + 29);

    auto gr_yyyy_x = pbuffer.data(idx_g_gp + 30);

    auto gr_yyyy_y = pbuffer.data(idx_g_gp + 31);

    auto gr_yyyy_z = pbuffer.data(idx_g_gp + 32);

    auto gr_yyyz_x = pbuffer.data(idx_g_gp + 33);

    auto gr_yyyz_y = pbuffer.data(idx_g_gp + 34);

    auto gr_yyyz_z = pbuffer.data(idx_g_gp + 35);

    auto gr_yyzz_x = pbuffer.data(idx_g_gp + 36);

    auto gr_yyzz_y = pbuffer.data(idx_g_gp + 37);

    auto gr_yyzz_z = pbuffer.data(idx_g_gp + 38);

    auto gr_yzzz_x = pbuffer.data(idx_g_gp + 39);

    auto gr_yzzz_y = pbuffer.data(idx_g_gp + 40);

    auto gr_yzzz_z = pbuffer.data(idx_g_gp + 41);

    auto gr_zzzz_x = pbuffer.data(idx_g_gp + 42);

    auto gr_zzzz_y = pbuffer.data(idx_g_gp + 43);

    auto gr_zzzz_z = pbuffer.data(idx_g_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_gd + 11);

    auto ts_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto ts_xxyz_xx = pbuffer.data(idx_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_gd + 29);

    auto ts_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto ts_xyyz_xx = pbuffer.data(idx_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto ts_xyzz_xx = pbuffer.data(idx_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_gd + 53);

    auto ts_xzzz_xx = pbuffer.data(idx_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto ts_yyyz_xx = pbuffer.data(idx_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto gr_xxxx_xx = pbuffer.data(idx_g_gd);

    auto gr_xxxx_xy = pbuffer.data(idx_g_gd + 1);

    auto gr_xxxx_xz = pbuffer.data(idx_g_gd + 2);

    auto gr_xxxx_yy = pbuffer.data(idx_g_gd + 3);

    auto gr_xxxx_yz = pbuffer.data(idx_g_gd + 4);

    auto gr_xxxx_zz = pbuffer.data(idx_g_gd + 5);

    auto gr_xxxy_xx = pbuffer.data(idx_g_gd + 6);

    auto gr_xxxy_xy = pbuffer.data(idx_g_gd + 7);

    auto gr_xxxy_xz = pbuffer.data(idx_g_gd + 8);

    auto gr_xxxy_yy = pbuffer.data(idx_g_gd + 9);

    auto gr_xxxy_yz = pbuffer.data(idx_g_gd + 10);

    auto gr_xxxy_zz = pbuffer.data(idx_g_gd + 11);

    auto gr_xxxz_xx = pbuffer.data(idx_g_gd + 12);

    auto gr_xxxz_xy = pbuffer.data(idx_g_gd + 13);

    auto gr_xxxz_xz = pbuffer.data(idx_g_gd + 14);

    auto gr_xxxz_yy = pbuffer.data(idx_g_gd + 15);

    auto gr_xxxz_yz = pbuffer.data(idx_g_gd + 16);

    auto gr_xxxz_zz = pbuffer.data(idx_g_gd + 17);

    auto gr_xxyy_xx = pbuffer.data(idx_g_gd + 18);

    auto gr_xxyy_xy = pbuffer.data(idx_g_gd + 19);

    auto gr_xxyy_xz = pbuffer.data(idx_g_gd + 20);

    auto gr_xxyy_yy = pbuffer.data(idx_g_gd + 21);

    auto gr_xxyy_yz = pbuffer.data(idx_g_gd + 22);

    auto gr_xxyy_zz = pbuffer.data(idx_g_gd + 23);

    auto gr_xxyz_xx = pbuffer.data(idx_g_gd + 24);

    auto gr_xxyz_xy = pbuffer.data(idx_g_gd + 25);

    auto gr_xxyz_xz = pbuffer.data(idx_g_gd + 26);

    auto gr_xxyz_yy = pbuffer.data(idx_g_gd + 27);

    auto gr_xxyz_yz = pbuffer.data(idx_g_gd + 28);

    auto gr_xxyz_zz = pbuffer.data(idx_g_gd + 29);

    auto gr_xxzz_xx = pbuffer.data(idx_g_gd + 30);

    auto gr_xxzz_xy = pbuffer.data(idx_g_gd + 31);

    auto gr_xxzz_xz = pbuffer.data(idx_g_gd + 32);

    auto gr_xxzz_yy = pbuffer.data(idx_g_gd + 33);

    auto gr_xxzz_yz = pbuffer.data(idx_g_gd + 34);

    auto gr_xxzz_zz = pbuffer.data(idx_g_gd + 35);

    auto gr_xyyy_xx = pbuffer.data(idx_g_gd + 36);

    auto gr_xyyy_xy = pbuffer.data(idx_g_gd + 37);

    auto gr_xyyy_xz = pbuffer.data(idx_g_gd + 38);

    auto gr_xyyy_yy = pbuffer.data(idx_g_gd + 39);

    auto gr_xyyy_yz = pbuffer.data(idx_g_gd + 40);

    auto gr_xyyy_zz = pbuffer.data(idx_g_gd + 41);

    auto gr_xyyz_xx = pbuffer.data(idx_g_gd + 42);

    auto gr_xyyz_xy = pbuffer.data(idx_g_gd + 43);

    auto gr_xyyz_xz = pbuffer.data(idx_g_gd + 44);

    auto gr_xyyz_yy = pbuffer.data(idx_g_gd + 45);

    auto gr_xyyz_yz = pbuffer.data(idx_g_gd + 46);

    auto gr_xyyz_zz = pbuffer.data(idx_g_gd + 47);

    auto gr_xyzz_xx = pbuffer.data(idx_g_gd + 48);

    auto gr_xyzz_xy = pbuffer.data(idx_g_gd + 49);

    auto gr_xyzz_xz = pbuffer.data(idx_g_gd + 50);

    auto gr_xyzz_yy = pbuffer.data(idx_g_gd + 51);

    auto gr_xyzz_yz = pbuffer.data(idx_g_gd + 52);

    auto gr_xyzz_zz = pbuffer.data(idx_g_gd + 53);

    auto gr_xzzz_xx = pbuffer.data(idx_g_gd + 54);

    auto gr_xzzz_xy = pbuffer.data(idx_g_gd + 55);

    auto gr_xzzz_xz = pbuffer.data(idx_g_gd + 56);

    auto gr_xzzz_yy = pbuffer.data(idx_g_gd + 57);

    auto gr_xzzz_yz = pbuffer.data(idx_g_gd + 58);

    auto gr_xzzz_zz = pbuffer.data(idx_g_gd + 59);

    auto gr_yyyy_xx = pbuffer.data(idx_g_gd + 60);

    auto gr_yyyy_xy = pbuffer.data(idx_g_gd + 61);

    auto gr_yyyy_xz = pbuffer.data(idx_g_gd + 62);

    auto gr_yyyy_yy = pbuffer.data(idx_g_gd + 63);

    auto gr_yyyy_yz = pbuffer.data(idx_g_gd + 64);

    auto gr_yyyy_zz = pbuffer.data(idx_g_gd + 65);

    auto gr_yyyz_xx = pbuffer.data(idx_g_gd + 66);

    auto gr_yyyz_xy = pbuffer.data(idx_g_gd + 67);

    auto gr_yyyz_xz = pbuffer.data(idx_g_gd + 68);

    auto gr_yyyz_yy = pbuffer.data(idx_g_gd + 69);

    auto gr_yyyz_yz = pbuffer.data(idx_g_gd + 70);

    auto gr_yyyz_zz = pbuffer.data(idx_g_gd + 71);

    auto gr_yyzz_xx = pbuffer.data(idx_g_gd + 72);

    auto gr_yyzz_xy = pbuffer.data(idx_g_gd + 73);

    auto gr_yyzz_xz = pbuffer.data(idx_g_gd + 74);

    auto gr_yyzz_yy = pbuffer.data(idx_g_gd + 75);

    auto gr_yyzz_yz = pbuffer.data(idx_g_gd + 76);

    auto gr_yyzz_zz = pbuffer.data(idx_g_gd + 77);

    auto gr_yzzz_xx = pbuffer.data(idx_g_gd + 78);

    auto gr_yzzz_xy = pbuffer.data(idx_g_gd + 79);

    auto gr_yzzz_xz = pbuffer.data(idx_g_gd + 80);

    auto gr_yzzz_yy = pbuffer.data(idx_g_gd + 81);

    auto gr_yzzz_yz = pbuffer.data(idx_g_gd + 82);

    auto gr_yzzz_zz = pbuffer.data(idx_g_gd + 83);

    auto gr_zzzz_xx = pbuffer.data(idx_g_gd + 84);

    auto gr_zzzz_xy = pbuffer.data(idx_g_gd + 85);

    auto gr_zzzz_xz = pbuffer.data(idx_g_gd + 86);

    auto gr_zzzz_yy = pbuffer.data(idx_g_gd + 87);

    auto gr_zzzz_yz = pbuffer.data(idx_g_gd + 88);

    auto gr_zzzz_zz = pbuffer.data(idx_g_gd + 89);

    // Set up 0-6 components of targeted buffer : GD

    auto grr_x_xxxx_xx = pbuffer.data(idx_gr_gd);

    auto grr_x_xxxx_xy = pbuffer.data(idx_gr_gd + 1);

    auto grr_x_xxxx_xz = pbuffer.data(idx_gr_gd + 2);

    auto grr_x_xxxx_yy = pbuffer.data(idx_gr_gd + 3);

    auto grr_x_xxxx_yz = pbuffer.data(idx_gr_gd + 4);

    auto grr_x_xxxx_zz = pbuffer.data(idx_gr_gd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xx, gr_xxx_xy, gr_xxx_xz, gr_xxx_yy, gr_xxx_yz, gr_xxx_zz, gr_xxxx_x, gr_xxxx_xx, gr_xxxx_xy, gr_xxxx_xz, gr_xxxx_y, gr_xxxx_yy, gr_xxxx_yz, gr_xxxx_z, gr_xxxx_zz, grr_x_xxxx_xx, grr_x_xxxx_xy, grr_x_xxxx_xz, grr_x_xxxx_yy, grr_x_xxxx_yz, grr_x_xxxx_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxx_xx[i] = 4.0 * ts_xxx_xx[i] * gfe2_0 + 4.0 * gr_xxx_xx[i] * gfe_0 + 2.0 * ts_xxxx_x[i] * gfe2_0 + 2.0 * gr_xxxx_x[i] * gfe_0 + ts_xxxx_xx[i] * gfe_0 * gc_x[i] + gr_xxxx_xx[i] * gc_x[i];

        grr_x_xxxx_xy[i] = 4.0 * ts_xxx_xy[i] * gfe2_0 + 4.0 * gr_xxx_xy[i] * gfe_0 + ts_xxxx_y[i] * gfe2_0 + gr_xxxx_y[i] * gfe_0 + ts_xxxx_xy[i] * gfe_0 * gc_x[i] + gr_xxxx_xy[i] * gc_x[i];

        grr_x_xxxx_xz[i] = 4.0 * ts_xxx_xz[i] * gfe2_0 + 4.0 * gr_xxx_xz[i] * gfe_0 + ts_xxxx_z[i] * gfe2_0 + gr_xxxx_z[i] * gfe_0 + ts_xxxx_xz[i] * gfe_0 * gc_x[i] + gr_xxxx_xz[i] * gc_x[i];

        grr_x_xxxx_yy[i] = 4.0 * ts_xxx_yy[i] * gfe2_0 + 4.0 * gr_xxx_yy[i] * gfe_0 + ts_xxxx_yy[i] * gfe_0 * gc_x[i] + gr_xxxx_yy[i] * gc_x[i];

        grr_x_xxxx_yz[i] = 4.0 * ts_xxx_yz[i] * gfe2_0 + 4.0 * gr_xxx_yz[i] * gfe_0 + ts_xxxx_yz[i] * gfe_0 * gc_x[i] + gr_xxxx_yz[i] * gc_x[i];

        grr_x_xxxx_zz[i] = 4.0 * ts_xxx_zz[i] * gfe2_0 + 4.0 * gr_xxx_zz[i] * gfe_0 + ts_xxxx_zz[i] * gfe_0 * gc_x[i] + gr_xxxx_zz[i] * gc_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto grr_x_xxxy_xx = pbuffer.data(idx_gr_gd + 6);

    auto grr_x_xxxy_xy = pbuffer.data(idx_gr_gd + 7);

    auto grr_x_xxxy_xz = pbuffer.data(idx_gr_gd + 8);

    auto grr_x_xxxy_yy = pbuffer.data(idx_gr_gd + 9);

    auto grr_x_xxxy_yz = pbuffer.data(idx_gr_gd + 10);

    auto grr_x_xxxy_zz = pbuffer.data(idx_gr_gd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_x, gr_xxxy_xx, gr_xxxy_xy, gr_xxxy_xz, gr_xxxy_y, gr_xxxy_yy, gr_xxxy_yz, gr_xxxy_z, gr_xxxy_zz, gr_xxy_xx, gr_xxy_xy, gr_xxy_xz, gr_xxy_yy, gr_xxy_yz, gr_xxy_zz, grr_x_xxxy_xx, grr_x_xxxy_xy, grr_x_xxxy_xz, grr_x_xxxy_yy, grr_x_xxxy_yz, grr_x_xxxy_zz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_yy, ts_xxy_yz, ts_xxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxy_xx[i] = 3.0 * ts_xxy_xx[i] * gfe2_0 + 3.0 * gr_xxy_xx[i] * gfe_0 + 2.0 * ts_xxxy_x[i] * gfe2_0 + 2.0 * gr_xxxy_x[i] * gfe_0 + ts_xxxy_xx[i] * gfe_0 * gc_x[i] + gr_xxxy_xx[i] * gc_x[i];

        grr_x_xxxy_xy[i] = 3.0 * ts_xxy_xy[i] * gfe2_0 + 3.0 * gr_xxy_xy[i] * gfe_0 + ts_xxxy_y[i] * gfe2_0 + gr_xxxy_y[i] * gfe_0 + ts_xxxy_xy[i] * gfe_0 * gc_x[i] + gr_xxxy_xy[i] * gc_x[i];

        grr_x_xxxy_xz[i] = 3.0 * ts_xxy_xz[i] * gfe2_0 + 3.0 * gr_xxy_xz[i] * gfe_0 + ts_xxxy_z[i] * gfe2_0 + gr_xxxy_z[i] * gfe_0 + ts_xxxy_xz[i] * gfe_0 * gc_x[i] + gr_xxxy_xz[i] * gc_x[i];

        grr_x_xxxy_yy[i] = 3.0 * ts_xxy_yy[i] * gfe2_0 + 3.0 * gr_xxy_yy[i] * gfe_0 + ts_xxxy_yy[i] * gfe_0 * gc_x[i] + gr_xxxy_yy[i] * gc_x[i];

        grr_x_xxxy_yz[i] = 3.0 * ts_xxy_yz[i] * gfe2_0 + 3.0 * gr_xxy_yz[i] * gfe_0 + ts_xxxy_yz[i] * gfe_0 * gc_x[i] + gr_xxxy_yz[i] * gc_x[i];

        grr_x_xxxy_zz[i] = 3.0 * ts_xxy_zz[i] * gfe2_0 + 3.0 * gr_xxy_zz[i] * gfe_0 + ts_xxxy_zz[i] * gfe_0 * gc_x[i] + gr_xxxy_zz[i] * gc_x[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto grr_x_xxxz_xx = pbuffer.data(idx_gr_gd + 12);

    auto grr_x_xxxz_xy = pbuffer.data(idx_gr_gd + 13);

    auto grr_x_xxxz_xz = pbuffer.data(idx_gr_gd + 14);

    auto grr_x_xxxz_yy = pbuffer.data(idx_gr_gd + 15);

    auto grr_x_xxxz_yz = pbuffer.data(idx_gr_gd + 16);

    auto grr_x_xxxz_zz = pbuffer.data(idx_gr_gd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_x, gr_xxxz_xx, gr_xxxz_xy, gr_xxxz_xz, gr_xxxz_y, gr_xxxz_yy, gr_xxxz_yz, gr_xxxz_z, gr_xxxz_zz, gr_xxz_xx, gr_xxz_xy, gr_xxz_xz, gr_xxz_yy, gr_xxz_yz, gr_xxz_zz, grr_x_xxxz_xx, grr_x_xxxz_xy, grr_x_xxxz_xz, grr_x_xxxz_yy, grr_x_xxxz_yz, grr_x_xxxz_zz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_yy, ts_xxz_yz, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxz_xx[i] = 3.0 * ts_xxz_xx[i] * gfe2_0 + 3.0 * gr_xxz_xx[i] * gfe_0 + 2.0 * ts_xxxz_x[i] * gfe2_0 + 2.0 * gr_xxxz_x[i] * gfe_0 + ts_xxxz_xx[i] * gfe_0 * gc_x[i] + gr_xxxz_xx[i] * gc_x[i];

        grr_x_xxxz_xy[i] = 3.0 * ts_xxz_xy[i] * gfe2_0 + 3.0 * gr_xxz_xy[i] * gfe_0 + ts_xxxz_y[i] * gfe2_0 + gr_xxxz_y[i] * gfe_0 + ts_xxxz_xy[i] * gfe_0 * gc_x[i] + gr_xxxz_xy[i] * gc_x[i];

        grr_x_xxxz_xz[i] = 3.0 * ts_xxz_xz[i] * gfe2_0 + 3.0 * gr_xxz_xz[i] * gfe_0 + ts_xxxz_z[i] * gfe2_0 + gr_xxxz_z[i] * gfe_0 + ts_xxxz_xz[i] * gfe_0 * gc_x[i] + gr_xxxz_xz[i] * gc_x[i];

        grr_x_xxxz_yy[i] = 3.0 * ts_xxz_yy[i] * gfe2_0 + 3.0 * gr_xxz_yy[i] * gfe_0 + ts_xxxz_yy[i] * gfe_0 * gc_x[i] + gr_xxxz_yy[i] * gc_x[i];

        grr_x_xxxz_yz[i] = 3.0 * ts_xxz_yz[i] * gfe2_0 + 3.0 * gr_xxz_yz[i] * gfe_0 + ts_xxxz_yz[i] * gfe_0 * gc_x[i] + gr_xxxz_yz[i] * gc_x[i];

        grr_x_xxxz_zz[i] = 3.0 * ts_xxz_zz[i] * gfe2_0 + 3.0 * gr_xxz_zz[i] * gfe_0 + ts_xxxz_zz[i] * gfe_0 * gc_x[i] + gr_xxxz_zz[i] * gc_x[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto grr_x_xxyy_xx = pbuffer.data(idx_gr_gd + 18);

    auto grr_x_xxyy_xy = pbuffer.data(idx_gr_gd + 19);

    auto grr_x_xxyy_xz = pbuffer.data(idx_gr_gd + 20);

    auto grr_x_xxyy_yy = pbuffer.data(idx_gr_gd + 21);

    auto grr_x_xxyy_yz = pbuffer.data(idx_gr_gd + 22);

    auto grr_x_xxyy_zz = pbuffer.data(idx_gr_gd + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_x, gr_xxyy_xx, gr_xxyy_xy, gr_xxyy_xz, gr_xxyy_y, gr_xxyy_yy, gr_xxyy_yz, gr_xxyy_z, gr_xxyy_zz, gr_xyy_xx, gr_xyy_xy, gr_xyy_xz, gr_xyy_yy, gr_xyy_yz, gr_xyy_zz, grr_x_xxyy_xx, grr_x_xxyy_xy, grr_x_xxyy_xz, grr_x_xxyy_yy, grr_x_xxyy_yz, grr_x_xxyy_zz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxyy_xx[i] = 2.0 * ts_xyy_xx[i] * gfe2_0 + 2.0 * gr_xyy_xx[i] * gfe_0 + 2.0 * ts_xxyy_x[i] * gfe2_0 + 2.0 * gr_xxyy_x[i] * gfe_0 + ts_xxyy_xx[i] * gfe_0 * gc_x[i] + gr_xxyy_xx[i] * gc_x[i];

        grr_x_xxyy_xy[i] = 2.0 * ts_xyy_xy[i] * gfe2_0 + 2.0 * gr_xyy_xy[i] * gfe_0 + ts_xxyy_y[i] * gfe2_0 + gr_xxyy_y[i] * gfe_0 + ts_xxyy_xy[i] * gfe_0 * gc_x[i] + gr_xxyy_xy[i] * gc_x[i];

        grr_x_xxyy_xz[i] = 2.0 * ts_xyy_xz[i] * gfe2_0 + 2.0 * gr_xyy_xz[i] * gfe_0 + ts_xxyy_z[i] * gfe2_0 + gr_xxyy_z[i] * gfe_0 + ts_xxyy_xz[i] * gfe_0 * gc_x[i] + gr_xxyy_xz[i] * gc_x[i];

        grr_x_xxyy_yy[i] = 2.0 * ts_xyy_yy[i] * gfe2_0 + 2.0 * gr_xyy_yy[i] * gfe_0 + ts_xxyy_yy[i] * gfe_0 * gc_x[i] + gr_xxyy_yy[i] * gc_x[i];

        grr_x_xxyy_yz[i] = 2.0 * ts_xyy_yz[i] * gfe2_0 + 2.0 * gr_xyy_yz[i] * gfe_0 + ts_xxyy_yz[i] * gfe_0 * gc_x[i] + gr_xxyy_yz[i] * gc_x[i];

        grr_x_xxyy_zz[i] = 2.0 * ts_xyy_zz[i] * gfe2_0 + 2.0 * gr_xyy_zz[i] * gfe_0 + ts_xxyy_zz[i] * gfe_0 * gc_x[i] + gr_xxyy_zz[i] * gc_x[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto grr_x_xxyz_xx = pbuffer.data(idx_gr_gd + 24);

    auto grr_x_xxyz_xy = pbuffer.data(idx_gr_gd + 25);

    auto grr_x_xxyz_xz = pbuffer.data(idx_gr_gd + 26);

    auto grr_x_xxyz_yy = pbuffer.data(idx_gr_gd + 27);

    auto grr_x_xxyz_yz = pbuffer.data(idx_gr_gd + 28);

    auto grr_x_xxyz_zz = pbuffer.data(idx_gr_gd + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_x, gr_xxyz_xx, gr_xxyz_xy, gr_xxyz_xz, gr_xxyz_y, gr_xxyz_yy, gr_xxyz_yz, gr_xxyz_z, gr_xxyz_zz, gr_xyz_xx, gr_xyz_xy, gr_xyz_xz, gr_xyz_yy, gr_xyz_yz, gr_xyz_zz, grr_x_xxyz_xx, grr_x_xxyz_xy, grr_x_xxyz_xz, grr_x_xxyz_yy, grr_x_xxyz_yz, grr_x_xxyz_zz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_yy, ts_xyz_yz, ts_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxyz_xx[i] = 2.0 * ts_xyz_xx[i] * gfe2_0 + 2.0 * gr_xyz_xx[i] * gfe_0 + 2.0 * ts_xxyz_x[i] * gfe2_0 + 2.0 * gr_xxyz_x[i] * gfe_0 + ts_xxyz_xx[i] * gfe_0 * gc_x[i] + gr_xxyz_xx[i] * gc_x[i];

        grr_x_xxyz_xy[i] = 2.0 * ts_xyz_xy[i] * gfe2_0 + 2.0 * gr_xyz_xy[i] * gfe_0 + ts_xxyz_y[i] * gfe2_0 + gr_xxyz_y[i] * gfe_0 + ts_xxyz_xy[i] * gfe_0 * gc_x[i] + gr_xxyz_xy[i] * gc_x[i];

        grr_x_xxyz_xz[i] = 2.0 * ts_xyz_xz[i] * gfe2_0 + 2.0 * gr_xyz_xz[i] * gfe_0 + ts_xxyz_z[i] * gfe2_0 + gr_xxyz_z[i] * gfe_0 + ts_xxyz_xz[i] * gfe_0 * gc_x[i] + gr_xxyz_xz[i] * gc_x[i];

        grr_x_xxyz_yy[i] = 2.0 * ts_xyz_yy[i] * gfe2_0 + 2.0 * gr_xyz_yy[i] * gfe_0 + ts_xxyz_yy[i] * gfe_0 * gc_x[i] + gr_xxyz_yy[i] * gc_x[i];

        grr_x_xxyz_yz[i] = 2.0 * ts_xyz_yz[i] * gfe2_0 + 2.0 * gr_xyz_yz[i] * gfe_0 + ts_xxyz_yz[i] * gfe_0 * gc_x[i] + gr_xxyz_yz[i] * gc_x[i];

        grr_x_xxyz_zz[i] = 2.0 * ts_xyz_zz[i] * gfe2_0 + 2.0 * gr_xyz_zz[i] * gfe_0 + ts_xxyz_zz[i] * gfe_0 * gc_x[i] + gr_xxyz_zz[i] * gc_x[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto grr_x_xxzz_xx = pbuffer.data(idx_gr_gd + 30);

    auto grr_x_xxzz_xy = pbuffer.data(idx_gr_gd + 31);

    auto grr_x_xxzz_xz = pbuffer.data(idx_gr_gd + 32);

    auto grr_x_xxzz_yy = pbuffer.data(idx_gr_gd + 33);

    auto grr_x_xxzz_yz = pbuffer.data(idx_gr_gd + 34);

    auto grr_x_xxzz_zz = pbuffer.data(idx_gr_gd + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_x, gr_xxzz_xx, gr_xxzz_xy, gr_xxzz_xz, gr_xxzz_y, gr_xxzz_yy, gr_xxzz_yz, gr_xxzz_z, gr_xxzz_zz, gr_xzz_xx, gr_xzz_xy, gr_xzz_xz, gr_xzz_yy, gr_xzz_yz, gr_xzz_zz, grr_x_xxzz_xx, grr_x_xxzz_xy, grr_x_xxzz_xz, grr_x_xxzz_yy, grr_x_xxzz_yz, grr_x_xxzz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxzz_xx[i] = 2.0 * ts_xzz_xx[i] * gfe2_0 + 2.0 * gr_xzz_xx[i] * gfe_0 + 2.0 * ts_xxzz_x[i] * gfe2_0 + 2.0 * gr_xxzz_x[i] * gfe_0 + ts_xxzz_xx[i] * gfe_0 * gc_x[i] + gr_xxzz_xx[i] * gc_x[i];

        grr_x_xxzz_xy[i] = 2.0 * ts_xzz_xy[i] * gfe2_0 + 2.0 * gr_xzz_xy[i] * gfe_0 + ts_xxzz_y[i] * gfe2_0 + gr_xxzz_y[i] * gfe_0 + ts_xxzz_xy[i] * gfe_0 * gc_x[i] + gr_xxzz_xy[i] * gc_x[i];

        grr_x_xxzz_xz[i] = 2.0 * ts_xzz_xz[i] * gfe2_0 + 2.0 * gr_xzz_xz[i] * gfe_0 + ts_xxzz_z[i] * gfe2_0 + gr_xxzz_z[i] * gfe_0 + ts_xxzz_xz[i] * gfe_0 * gc_x[i] + gr_xxzz_xz[i] * gc_x[i];

        grr_x_xxzz_yy[i] = 2.0 * ts_xzz_yy[i] * gfe2_0 + 2.0 * gr_xzz_yy[i] * gfe_0 + ts_xxzz_yy[i] * gfe_0 * gc_x[i] + gr_xxzz_yy[i] * gc_x[i];

        grr_x_xxzz_yz[i] = 2.0 * ts_xzz_yz[i] * gfe2_0 + 2.0 * gr_xzz_yz[i] * gfe_0 + ts_xxzz_yz[i] * gfe_0 * gc_x[i] + gr_xxzz_yz[i] * gc_x[i];

        grr_x_xxzz_zz[i] = 2.0 * ts_xzz_zz[i] * gfe2_0 + 2.0 * gr_xzz_zz[i] * gfe_0 + ts_xxzz_zz[i] * gfe_0 * gc_x[i] + gr_xxzz_zz[i] * gc_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto grr_x_xyyy_xx = pbuffer.data(idx_gr_gd + 36);

    auto grr_x_xyyy_xy = pbuffer.data(idx_gr_gd + 37);

    auto grr_x_xyyy_xz = pbuffer.data(idx_gr_gd + 38);

    auto grr_x_xyyy_yy = pbuffer.data(idx_gr_gd + 39);

    auto grr_x_xyyy_yz = pbuffer.data(idx_gr_gd + 40);

    auto grr_x_xyyy_zz = pbuffer.data(idx_gr_gd + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_x, gr_xyyy_xx, gr_xyyy_xy, gr_xyyy_xz, gr_xyyy_y, gr_xyyy_yy, gr_xyyy_yz, gr_xyyy_z, gr_xyyy_zz, gr_yyy_xx, gr_yyy_xy, gr_yyy_xz, gr_yyy_yy, gr_yyy_yz, gr_yyy_zz, grr_x_xyyy_xx, grr_x_xyyy_xy, grr_x_xyyy_xz, grr_x_xyyy_yy, grr_x_xyyy_yz, grr_x_xyyy_zz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyyy_xx[i] = ts_yyy_xx[i] * gfe2_0 + gr_yyy_xx[i] * gfe_0 + 2.0 * ts_xyyy_x[i] * gfe2_0 + 2.0 * gr_xyyy_x[i] * gfe_0 + ts_xyyy_xx[i] * gfe_0 * gc_x[i] + gr_xyyy_xx[i] * gc_x[i];

        grr_x_xyyy_xy[i] = ts_yyy_xy[i] * gfe2_0 + gr_yyy_xy[i] * gfe_0 + ts_xyyy_y[i] * gfe2_0 + gr_xyyy_y[i] * gfe_0 + ts_xyyy_xy[i] * gfe_0 * gc_x[i] + gr_xyyy_xy[i] * gc_x[i];

        grr_x_xyyy_xz[i] = ts_yyy_xz[i] * gfe2_0 + gr_yyy_xz[i] * gfe_0 + ts_xyyy_z[i] * gfe2_0 + gr_xyyy_z[i] * gfe_0 + ts_xyyy_xz[i] * gfe_0 * gc_x[i] + gr_xyyy_xz[i] * gc_x[i];

        grr_x_xyyy_yy[i] = ts_yyy_yy[i] * gfe2_0 + gr_yyy_yy[i] * gfe_0 + ts_xyyy_yy[i] * gfe_0 * gc_x[i] + gr_xyyy_yy[i] * gc_x[i];

        grr_x_xyyy_yz[i] = ts_yyy_yz[i] * gfe2_0 + gr_yyy_yz[i] * gfe_0 + ts_xyyy_yz[i] * gfe_0 * gc_x[i] + gr_xyyy_yz[i] * gc_x[i];

        grr_x_xyyy_zz[i] = ts_yyy_zz[i] * gfe2_0 + gr_yyy_zz[i] * gfe_0 + ts_xyyy_zz[i] * gfe_0 * gc_x[i] + gr_xyyy_zz[i] * gc_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto grr_x_xyyz_xx = pbuffer.data(idx_gr_gd + 42);

    auto grr_x_xyyz_xy = pbuffer.data(idx_gr_gd + 43);

    auto grr_x_xyyz_xz = pbuffer.data(idx_gr_gd + 44);

    auto grr_x_xyyz_yy = pbuffer.data(idx_gr_gd + 45);

    auto grr_x_xyyz_yz = pbuffer.data(idx_gr_gd + 46);

    auto grr_x_xyyz_zz = pbuffer.data(idx_gr_gd + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_x, gr_xyyz_xx, gr_xyyz_xy, gr_xyyz_xz, gr_xyyz_y, gr_xyyz_yy, gr_xyyz_yz, gr_xyyz_z, gr_xyyz_zz, gr_yyz_xx, gr_yyz_xy, gr_yyz_xz, gr_yyz_yy, gr_yyz_yz, gr_yyz_zz, grr_x_xyyz_xx, grr_x_xyyz_xy, grr_x_xyyz_xz, grr_x_xyyz_yy, grr_x_xyyz_yz, grr_x_xyyz_zz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_yy, ts_yyz_yz, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyyz_xx[i] = ts_yyz_xx[i] * gfe2_0 + gr_yyz_xx[i] * gfe_0 + 2.0 * ts_xyyz_x[i] * gfe2_0 + 2.0 * gr_xyyz_x[i] * gfe_0 + ts_xyyz_xx[i] * gfe_0 * gc_x[i] + gr_xyyz_xx[i] * gc_x[i];

        grr_x_xyyz_xy[i] = ts_yyz_xy[i] * gfe2_0 + gr_yyz_xy[i] * gfe_0 + ts_xyyz_y[i] * gfe2_0 + gr_xyyz_y[i] * gfe_0 + ts_xyyz_xy[i] * gfe_0 * gc_x[i] + gr_xyyz_xy[i] * gc_x[i];

        grr_x_xyyz_xz[i] = ts_yyz_xz[i] * gfe2_0 + gr_yyz_xz[i] * gfe_0 + ts_xyyz_z[i] * gfe2_0 + gr_xyyz_z[i] * gfe_0 + ts_xyyz_xz[i] * gfe_0 * gc_x[i] + gr_xyyz_xz[i] * gc_x[i];

        grr_x_xyyz_yy[i] = ts_yyz_yy[i] * gfe2_0 + gr_yyz_yy[i] * gfe_0 + ts_xyyz_yy[i] * gfe_0 * gc_x[i] + gr_xyyz_yy[i] * gc_x[i];

        grr_x_xyyz_yz[i] = ts_yyz_yz[i] * gfe2_0 + gr_yyz_yz[i] * gfe_0 + ts_xyyz_yz[i] * gfe_0 * gc_x[i] + gr_xyyz_yz[i] * gc_x[i];

        grr_x_xyyz_zz[i] = ts_yyz_zz[i] * gfe2_0 + gr_yyz_zz[i] * gfe_0 + ts_xyyz_zz[i] * gfe_0 * gc_x[i] + gr_xyyz_zz[i] * gc_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto grr_x_xyzz_xx = pbuffer.data(idx_gr_gd + 48);

    auto grr_x_xyzz_xy = pbuffer.data(idx_gr_gd + 49);

    auto grr_x_xyzz_xz = pbuffer.data(idx_gr_gd + 50);

    auto grr_x_xyzz_yy = pbuffer.data(idx_gr_gd + 51);

    auto grr_x_xyzz_yz = pbuffer.data(idx_gr_gd + 52);

    auto grr_x_xyzz_zz = pbuffer.data(idx_gr_gd + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_x, gr_xyzz_xx, gr_xyzz_xy, gr_xyzz_xz, gr_xyzz_y, gr_xyzz_yy, gr_xyzz_yz, gr_xyzz_z, gr_xyzz_zz, gr_yzz_xx, gr_yzz_xy, gr_yzz_xz, gr_yzz_yy, gr_yzz_yz, gr_yzz_zz, grr_x_xyzz_xx, grr_x_xyzz_xy, grr_x_xyzz_xz, grr_x_xyzz_yy, grr_x_xyzz_yz, grr_x_xyzz_zz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_yy, ts_yzz_yz, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyzz_xx[i] = ts_yzz_xx[i] * gfe2_0 + gr_yzz_xx[i] * gfe_0 + 2.0 * ts_xyzz_x[i] * gfe2_0 + 2.0 * gr_xyzz_x[i] * gfe_0 + ts_xyzz_xx[i] * gfe_0 * gc_x[i] + gr_xyzz_xx[i] * gc_x[i];

        grr_x_xyzz_xy[i] = ts_yzz_xy[i] * gfe2_0 + gr_yzz_xy[i] * gfe_0 + ts_xyzz_y[i] * gfe2_0 + gr_xyzz_y[i] * gfe_0 + ts_xyzz_xy[i] * gfe_0 * gc_x[i] + gr_xyzz_xy[i] * gc_x[i];

        grr_x_xyzz_xz[i] = ts_yzz_xz[i] * gfe2_0 + gr_yzz_xz[i] * gfe_0 + ts_xyzz_z[i] * gfe2_0 + gr_xyzz_z[i] * gfe_0 + ts_xyzz_xz[i] * gfe_0 * gc_x[i] + gr_xyzz_xz[i] * gc_x[i];

        grr_x_xyzz_yy[i] = ts_yzz_yy[i] * gfe2_0 + gr_yzz_yy[i] * gfe_0 + ts_xyzz_yy[i] * gfe_0 * gc_x[i] + gr_xyzz_yy[i] * gc_x[i];

        grr_x_xyzz_yz[i] = ts_yzz_yz[i] * gfe2_0 + gr_yzz_yz[i] * gfe_0 + ts_xyzz_yz[i] * gfe_0 * gc_x[i] + gr_xyzz_yz[i] * gc_x[i];

        grr_x_xyzz_zz[i] = ts_yzz_zz[i] * gfe2_0 + gr_yzz_zz[i] * gfe_0 + ts_xyzz_zz[i] * gfe_0 * gc_x[i] + gr_xyzz_zz[i] * gc_x[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto grr_x_xzzz_xx = pbuffer.data(idx_gr_gd + 54);

    auto grr_x_xzzz_xy = pbuffer.data(idx_gr_gd + 55);

    auto grr_x_xzzz_xz = pbuffer.data(idx_gr_gd + 56);

    auto grr_x_xzzz_yy = pbuffer.data(idx_gr_gd + 57);

    auto grr_x_xzzz_yz = pbuffer.data(idx_gr_gd + 58);

    auto grr_x_xzzz_zz = pbuffer.data(idx_gr_gd + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_x, gr_xzzz_xx, gr_xzzz_xy, gr_xzzz_xz, gr_xzzz_y, gr_xzzz_yy, gr_xzzz_yz, gr_xzzz_z, gr_xzzz_zz, gr_zzz_xx, gr_zzz_xy, gr_zzz_xz, gr_zzz_yy, gr_zzz_yz, gr_zzz_zz, grr_x_xzzz_xx, grr_x_xzzz_xy, grr_x_xzzz_xz, grr_x_xzzz_yy, grr_x_xzzz_yz, grr_x_xzzz_zz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xzzz_xx[i] = ts_zzz_xx[i] * gfe2_0 + gr_zzz_xx[i] * gfe_0 + 2.0 * ts_xzzz_x[i] * gfe2_0 + 2.0 * gr_xzzz_x[i] * gfe_0 + ts_xzzz_xx[i] * gfe_0 * gc_x[i] + gr_xzzz_xx[i] * gc_x[i];

        grr_x_xzzz_xy[i] = ts_zzz_xy[i] * gfe2_0 + gr_zzz_xy[i] * gfe_0 + ts_xzzz_y[i] * gfe2_0 + gr_xzzz_y[i] * gfe_0 + ts_xzzz_xy[i] * gfe_0 * gc_x[i] + gr_xzzz_xy[i] * gc_x[i];

        grr_x_xzzz_xz[i] = ts_zzz_xz[i] * gfe2_0 + gr_zzz_xz[i] * gfe_0 + ts_xzzz_z[i] * gfe2_0 + gr_xzzz_z[i] * gfe_0 + ts_xzzz_xz[i] * gfe_0 * gc_x[i] + gr_xzzz_xz[i] * gc_x[i];

        grr_x_xzzz_yy[i] = ts_zzz_yy[i] * gfe2_0 + gr_zzz_yy[i] * gfe_0 + ts_xzzz_yy[i] * gfe_0 * gc_x[i] + gr_xzzz_yy[i] * gc_x[i];

        grr_x_xzzz_yz[i] = ts_zzz_yz[i] * gfe2_0 + gr_zzz_yz[i] * gfe_0 + ts_xzzz_yz[i] * gfe_0 * gc_x[i] + gr_xzzz_yz[i] * gc_x[i];

        grr_x_xzzz_zz[i] = ts_zzz_zz[i] * gfe2_0 + gr_zzz_zz[i] * gfe_0 + ts_xzzz_zz[i] * gfe_0 * gc_x[i] + gr_xzzz_zz[i] * gc_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto grr_x_yyyy_xx = pbuffer.data(idx_gr_gd + 60);

    auto grr_x_yyyy_xy = pbuffer.data(idx_gr_gd + 61);

    auto grr_x_yyyy_xz = pbuffer.data(idx_gr_gd + 62);

    auto grr_x_yyyy_yy = pbuffer.data(idx_gr_gd + 63);

    auto grr_x_yyyy_yz = pbuffer.data(idx_gr_gd + 64);

    auto grr_x_yyyy_zz = pbuffer.data(idx_gr_gd + 65);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_x, gr_yyyy_xx, gr_yyyy_xy, gr_yyyy_xz, gr_yyyy_y, gr_yyyy_yy, gr_yyyy_yz, gr_yyyy_z, gr_yyyy_zz, grr_x_yyyy_xx, grr_x_yyyy_xy, grr_x_yyyy_xz, grr_x_yyyy_yy, grr_x_yyyy_yz, grr_x_yyyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyyy_xx[i] = 2.0 * ts_yyyy_x[i] * gfe2_0 + 2.0 * gr_yyyy_x[i] * gfe_0 + ts_yyyy_xx[i] * gfe_0 * gc_x[i] + gr_yyyy_xx[i] * gc_x[i];

        grr_x_yyyy_xy[i] = ts_yyyy_y[i] * gfe2_0 + gr_yyyy_y[i] * gfe_0 + ts_yyyy_xy[i] * gfe_0 * gc_x[i] + gr_yyyy_xy[i] * gc_x[i];

        grr_x_yyyy_xz[i] = ts_yyyy_z[i] * gfe2_0 + gr_yyyy_z[i] * gfe_0 + ts_yyyy_xz[i] * gfe_0 * gc_x[i] + gr_yyyy_xz[i] * gc_x[i];

        grr_x_yyyy_yy[i] = ts_yyyy_yy[i] * gfe_0 * gc_x[i] + gr_yyyy_yy[i] * gc_x[i];

        grr_x_yyyy_yz[i] = ts_yyyy_yz[i] * gfe_0 * gc_x[i] + gr_yyyy_yz[i] * gc_x[i];

        grr_x_yyyy_zz[i] = ts_yyyy_zz[i] * gfe_0 * gc_x[i] + gr_yyyy_zz[i] * gc_x[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto grr_x_yyyz_xx = pbuffer.data(idx_gr_gd + 66);

    auto grr_x_yyyz_xy = pbuffer.data(idx_gr_gd + 67);

    auto grr_x_yyyz_xz = pbuffer.data(idx_gr_gd + 68);

    auto grr_x_yyyz_yy = pbuffer.data(idx_gr_gd + 69);

    auto grr_x_yyyz_yz = pbuffer.data(idx_gr_gd + 70);

    auto grr_x_yyyz_zz = pbuffer.data(idx_gr_gd + 71);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_x, gr_yyyz_xx, gr_yyyz_xy, gr_yyyz_xz, gr_yyyz_y, gr_yyyz_yy, gr_yyyz_yz, gr_yyyz_z, gr_yyyz_zz, grr_x_yyyz_xx, grr_x_yyyz_xy, grr_x_yyyz_xz, grr_x_yyyz_yy, grr_x_yyyz_yz, grr_x_yyyz_zz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyyz_xx[i] = 2.0 * ts_yyyz_x[i] * gfe2_0 + 2.0 * gr_yyyz_x[i] * gfe_0 + ts_yyyz_xx[i] * gfe_0 * gc_x[i] + gr_yyyz_xx[i] * gc_x[i];

        grr_x_yyyz_xy[i] = ts_yyyz_y[i] * gfe2_0 + gr_yyyz_y[i] * gfe_0 + ts_yyyz_xy[i] * gfe_0 * gc_x[i] + gr_yyyz_xy[i] * gc_x[i];

        grr_x_yyyz_xz[i] = ts_yyyz_z[i] * gfe2_0 + gr_yyyz_z[i] * gfe_0 + ts_yyyz_xz[i] * gfe_0 * gc_x[i] + gr_yyyz_xz[i] * gc_x[i];

        grr_x_yyyz_yy[i] = ts_yyyz_yy[i] * gfe_0 * gc_x[i] + gr_yyyz_yy[i] * gc_x[i];

        grr_x_yyyz_yz[i] = ts_yyyz_yz[i] * gfe_0 * gc_x[i] + gr_yyyz_yz[i] * gc_x[i];

        grr_x_yyyz_zz[i] = ts_yyyz_zz[i] * gfe_0 * gc_x[i] + gr_yyyz_zz[i] * gc_x[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto grr_x_yyzz_xx = pbuffer.data(idx_gr_gd + 72);

    auto grr_x_yyzz_xy = pbuffer.data(idx_gr_gd + 73);

    auto grr_x_yyzz_xz = pbuffer.data(idx_gr_gd + 74);

    auto grr_x_yyzz_yy = pbuffer.data(idx_gr_gd + 75);

    auto grr_x_yyzz_yz = pbuffer.data(idx_gr_gd + 76);

    auto grr_x_yyzz_zz = pbuffer.data(idx_gr_gd + 77);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_x, gr_yyzz_xx, gr_yyzz_xy, gr_yyzz_xz, gr_yyzz_y, gr_yyzz_yy, gr_yyzz_yz, gr_yyzz_z, gr_yyzz_zz, grr_x_yyzz_xx, grr_x_yyzz_xy, grr_x_yyzz_xz, grr_x_yyzz_yy, grr_x_yyzz_yz, grr_x_yyzz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyzz_xx[i] = 2.0 * ts_yyzz_x[i] * gfe2_0 + 2.0 * gr_yyzz_x[i] * gfe_0 + ts_yyzz_xx[i] * gfe_0 * gc_x[i] + gr_yyzz_xx[i] * gc_x[i];

        grr_x_yyzz_xy[i] = ts_yyzz_y[i] * gfe2_0 + gr_yyzz_y[i] * gfe_0 + ts_yyzz_xy[i] * gfe_0 * gc_x[i] + gr_yyzz_xy[i] * gc_x[i];

        grr_x_yyzz_xz[i] = ts_yyzz_z[i] * gfe2_0 + gr_yyzz_z[i] * gfe_0 + ts_yyzz_xz[i] * gfe_0 * gc_x[i] + gr_yyzz_xz[i] * gc_x[i];

        grr_x_yyzz_yy[i] = ts_yyzz_yy[i] * gfe_0 * gc_x[i] + gr_yyzz_yy[i] * gc_x[i];

        grr_x_yyzz_yz[i] = ts_yyzz_yz[i] * gfe_0 * gc_x[i] + gr_yyzz_yz[i] * gc_x[i];

        grr_x_yyzz_zz[i] = ts_yyzz_zz[i] * gfe_0 * gc_x[i] + gr_yyzz_zz[i] * gc_x[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto grr_x_yzzz_xx = pbuffer.data(idx_gr_gd + 78);

    auto grr_x_yzzz_xy = pbuffer.data(idx_gr_gd + 79);

    auto grr_x_yzzz_xz = pbuffer.data(idx_gr_gd + 80);

    auto grr_x_yzzz_yy = pbuffer.data(idx_gr_gd + 81);

    auto grr_x_yzzz_yz = pbuffer.data(idx_gr_gd + 82);

    auto grr_x_yzzz_zz = pbuffer.data(idx_gr_gd + 83);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_x, gr_yzzz_xx, gr_yzzz_xy, gr_yzzz_xz, gr_yzzz_y, gr_yzzz_yy, gr_yzzz_yz, gr_yzzz_z, gr_yzzz_zz, grr_x_yzzz_xx, grr_x_yzzz_xy, grr_x_yzzz_xz, grr_x_yzzz_yy, grr_x_yzzz_yz, grr_x_yzzz_zz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yzzz_xx[i] = 2.0 * ts_yzzz_x[i] * gfe2_0 + 2.0 * gr_yzzz_x[i] * gfe_0 + ts_yzzz_xx[i] * gfe_0 * gc_x[i] + gr_yzzz_xx[i] * gc_x[i];

        grr_x_yzzz_xy[i] = ts_yzzz_y[i] * gfe2_0 + gr_yzzz_y[i] * gfe_0 + ts_yzzz_xy[i] * gfe_0 * gc_x[i] + gr_yzzz_xy[i] * gc_x[i];

        grr_x_yzzz_xz[i] = ts_yzzz_z[i] * gfe2_0 + gr_yzzz_z[i] * gfe_0 + ts_yzzz_xz[i] * gfe_0 * gc_x[i] + gr_yzzz_xz[i] * gc_x[i];

        grr_x_yzzz_yy[i] = ts_yzzz_yy[i] * gfe_0 * gc_x[i] + gr_yzzz_yy[i] * gc_x[i];

        grr_x_yzzz_yz[i] = ts_yzzz_yz[i] * gfe_0 * gc_x[i] + gr_yzzz_yz[i] * gc_x[i];

        grr_x_yzzz_zz[i] = ts_yzzz_zz[i] * gfe_0 * gc_x[i] + gr_yzzz_zz[i] * gc_x[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto grr_x_zzzz_xx = pbuffer.data(idx_gr_gd + 84);

    auto grr_x_zzzz_xy = pbuffer.data(idx_gr_gd + 85);

    auto grr_x_zzzz_xz = pbuffer.data(idx_gr_gd + 86);

    auto grr_x_zzzz_yy = pbuffer.data(idx_gr_gd + 87);

    auto grr_x_zzzz_yz = pbuffer.data(idx_gr_gd + 88);

    auto grr_x_zzzz_zz = pbuffer.data(idx_gr_gd + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_x, gr_zzzz_xx, gr_zzzz_xy, gr_zzzz_xz, gr_zzzz_y, gr_zzzz_yy, gr_zzzz_yz, gr_zzzz_z, gr_zzzz_zz, grr_x_zzzz_xx, grr_x_zzzz_xy, grr_x_zzzz_xz, grr_x_zzzz_yy, grr_x_zzzz_yz, grr_x_zzzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zzzz_xx[i] = 2.0 * ts_zzzz_x[i] * gfe2_0 + 2.0 * gr_zzzz_x[i] * gfe_0 + ts_zzzz_xx[i] * gfe_0 * gc_x[i] + gr_zzzz_xx[i] * gc_x[i];

        grr_x_zzzz_xy[i] = ts_zzzz_y[i] * gfe2_0 + gr_zzzz_y[i] * gfe_0 + ts_zzzz_xy[i] * gfe_0 * gc_x[i] + gr_zzzz_xy[i] * gc_x[i];

        grr_x_zzzz_xz[i] = ts_zzzz_z[i] * gfe2_0 + gr_zzzz_z[i] * gfe_0 + ts_zzzz_xz[i] * gfe_0 * gc_x[i] + gr_zzzz_xz[i] * gc_x[i];

        grr_x_zzzz_yy[i] = ts_zzzz_yy[i] * gfe_0 * gc_x[i] + gr_zzzz_yy[i] * gc_x[i];

        grr_x_zzzz_yz[i] = ts_zzzz_yz[i] * gfe_0 * gc_x[i] + gr_zzzz_yz[i] * gc_x[i];

        grr_x_zzzz_zz[i] = ts_zzzz_zz[i] * gfe_0 * gc_x[i] + gr_zzzz_zz[i] * gc_x[i];
    }

    // Set up 90-96 components of targeted buffer : GD

    auto grr_y_xxxx_xx = pbuffer.data(idx_gr_gd + 90);

    auto grr_y_xxxx_xy = pbuffer.data(idx_gr_gd + 91);

    auto grr_y_xxxx_xz = pbuffer.data(idx_gr_gd + 92);

    auto grr_y_xxxx_yy = pbuffer.data(idx_gr_gd + 93);

    auto grr_y_xxxx_yz = pbuffer.data(idx_gr_gd + 94);

    auto grr_y_xxxx_zz = pbuffer.data(idx_gr_gd + 95);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_x, gr_xxxx_xx, gr_xxxx_xy, gr_xxxx_xz, gr_xxxx_y, gr_xxxx_yy, gr_xxxx_yz, gr_xxxx_z, gr_xxxx_zz, grr_y_xxxx_xx, grr_y_xxxx_xy, grr_y_xxxx_xz, grr_y_xxxx_yy, grr_y_xxxx_yz, grr_y_xxxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxxx_xx[i] = ts_xxxx_xx[i] * gfe_0 * gc_y[i] + gr_xxxx_xx[i] * gc_y[i];

        grr_y_xxxx_xy[i] = ts_xxxx_x[i] * gfe2_0 + gr_xxxx_x[i] * gfe_0 + ts_xxxx_xy[i] * gfe_0 * gc_y[i] + gr_xxxx_xy[i] * gc_y[i];

        grr_y_xxxx_xz[i] = ts_xxxx_xz[i] * gfe_0 * gc_y[i] + gr_xxxx_xz[i] * gc_y[i];

        grr_y_xxxx_yy[i] = 2.0 * ts_xxxx_y[i] * gfe2_0 + 2.0 * gr_xxxx_y[i] * gfe_0 + ts_xxxx_yy[i] * gfe_0 * gc_y[i] + gr_xxxx_yy[i] * gc_y[i];

        grr_y_xxxx_yz[i] = ts_xxxx_z[i] * gfe2_0 + gr_xxxx_z[i] * gfe_0 + ts_xxxx_yz[i] * gfe_0 * gc_y[i] + gr_xxxx_yz[i] * gc_y[i];

        grr_y_xxxx_zz[i] = ts_xxxx_zz[i] * gfe_0 * gc_y[i] + gr_xxxx_zz[i] * gc_y[i];
    }

    // Set up 96-102 components of targeted buffer : GD

    auto grr_y_xxxy_xx = pbuffer.data(idx_gr_gd + 96);

    auto grr_y_xxxy_xy = pbuffer.data(idx_gr_gd + 97);

    auto grr_y_xxxy_xz = pbuffer.data(idx_gr_gd + 98);

    auto grr_y_xxxy_yy = pbuffer.data(idx_gr_gd + 99);

    auto grr_y_xxxy_yz = pbuffer.data(idx_gr_gd + 100);

    auto grr_y_xxxy_zz = pbuffer.data(idx_gr_gd + 101);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xx, gr_xxx_xy, gr_xxx_xz, gr_xxx_yy, gr_xxx_yz, gr_xxx_zz, gr_xxxy_x, gr_xxxy_xx, gr_xxxy_xy, gr_xxxy_xz, gr_xxxy_y, gr_xxxy_yy, gr_xxxy_yz, gr_xxxy_z, gr_xxxy_zz, grr_y_xxxy_xx, grr_y_xxxy_xy, grr_y_xxxy_xz, grr_y_xxxy_yy, grr_y_xxxy_yz, grr_y_xxxy_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxxy_xx[i] = ts_xxx_xx[i] * gfe2_0 + gr_xxx_xx[i] * gfe_0 + ts_xxxy_xx[i] * gfe_0 * gc_y[i] + gr_xxxy_xx[i] * gc_y[i];

        grr_y_xxxy_xy[i] = ts_xxx_xy[i] * gfe2_0 + gr_xxx_xy[i] * gfe_0 + ts_xxxy_x[i] * gfe2_0 + gr_xxxy_x[i] * gfe_0 + ts_xxxy_xy[i] * gfe_0 * gc_y[i] + gr_xxxy_xy[i] * gc_y[i];

        grr_y_xxxy_xz[i] = ts_xxx_xz[i] * gfe2_0 + gr_xxx_xz[i] * gfe_0 + ts_xxxy_xz[i] * gfe_0 * gc_y[i] + gr_xxxy_xz[i] * gc_y[i];

        grr_y_xxxy_yy[i] = ts_xxx_yy[i] * gfe2_0 + gr_xxx_yy[i] * gfe_0 + 2.0 * ts_xxxy_y[i] * gfe2_0 + 2.0 * gr_xxxy_y[i] * gfe_0 + ts_xxxy_yy[i] * gfe_0 * gc_y[i] + gr_xxxy_yy[i] * gc_y[i];

        grr_y_xxxy_yz[i] = ts_xxx_yz[i] * gfe2_0 + gr_xxx_yz[i] * gfe_0 + ts_xxxy_z[i] * gfe2_0 + gr_xxxy_z[i] * gfe_0 + ts_xxxy_yz[i] * gfe_0 * gc_y[i] + gr_xxxy_yz[i] * gc_y[i];

        grr_y_xxxy_zz[i] = ts_xxx_zz[i] * gfe2_0 + gr_xxx_zz[i] * gfe_0 + ts_xxxy_zz[i] * gfe_0 * gc_y[i] + gr_xxxy_zz[i] * gc_y[i];
    }

    // Set up 102-108 components of targeted buffer : GD

    auto grr_y_xxxz_xx = pbuffer.data(idx_gr_gd + 102);

    auto grr_y_xxxz_xy = pbuffer.data(idx_gr_gd + 103);

    auto grr_y_xxxz_xz = pbuffer.data(idx_gr_gd + 104);

    auto grr_y_xxxz_yy = pbuffer.data(idx_gr_gd + 105);

    auto grr_y_xxxz_yz = pbuffer.data(idx_gr_gd + 106);

    auto grr_y_xxxz_zz = pbuffer.data(idx_gr_gd + 107);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_x, gr_xxxz_xx, gr_xxxz_xy, gr_xxxz_xz, gr_xxxz_y, gr_xxxz_yy, gr_xxxz_yz, gr_xxxz_z, gr_xxxz_zz, grr_y_xxxz_xx, grr_y_xxxz_xy, grr_y_xxxz_xz, grr_y_xxxz_yy, grr_y_xxxz_yz, grr_y_xxxz_zz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxxz_xx[i] = ts_xxxz_xx[i] * gfe_0 * gc_y[i] + gr_xxxz_xx[i] * gc_y[i];

        grr_y_xxxz_xy[i] = ts_xxxz_x[i] * gfe2_0 + gr_xxxz_x[i] * gfe_0 + ts_xxxz_xy[i] * gfe_0 * gc_y[i] + gr_xxxz_xy[i] * gc_y[i];

        grr_y_xxxz_xz[i] = ts_xxxz_xz[i] * gfe_0 * gc_y[i] + gr_xxxz_xz[i] * gc_y[i];

        grr_y_xxxz_yy[i] = 2.0 * ts_xxxz_y[i] * gfe2_0 + 2.0 * gr_xxxz_y[i] * gfe_0 + ts_xxxz_yy[i] * gfe_0 * gc_y[i] + gr_xxxz_yy[i] * gc_y[i];

        grr_y_xxxz_yz[i] = ts_xxxz_z[i] * gfe2_0 + gr_xxxz_z[i] * gfe_0 + ts_xxxz_yz[i] * gfe_0 * gc_y[i] + gr_xxxz_yz[i] * gc_y[i];

        grr_y_xxxz_zz[i] = ts_xxxz_zz[i] * gfe_0 * gc_y[i] + gr_xxxz_zz[i] * gc_y[i];
    }

    // Set up 108-114 components of targeted buffer : GD

    auto grr_y_xxyy_xx = pbuffer.data(idx_gr_gd + 108);

    auto grr_y_xxyy_xy = pbuffer.data(idx_gr_gd + 109);

    auto grr_y_xxyy_xz = pbuffer.data(idx_gr_gd + 110);

    auto grr_y_xxyy_yy = pbuffer.data(idx_gr_gd + 111);

    auto grr_y_xxyy_yz = pbuffer.data(idx_gr_gd + 112);

    auto grr_y_xxyy_zz = pbuffer.data(idx_gr_gd + 113);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xx, gr_xxy_xy, gr_xxy_xz, gr_xxy_yy, gr_xxy_yz, gr_xxy_zz, gr_xxyy_x, gr_xxyy_xx, gr_xxyy_xy, gr_xxyy_xz, gr_xxyy_y, gr_xxyy_yy, gr_xxyy_yz, gr_xxyy_z, gr_xxyy_zz, grr_y_xxyy_xx, grr_y_xxyy_xy, grr_y_xxyy_xz, grr_y_xxyy_yy, grr_y_xxyy_yz, grr_y_xxyy_zz, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_yy, ts_xxy_yz, ts_xxy_zz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxyy_xx[i] = 2.0 * ts_xxy_xx[i] * gfe2_0 + 2.0 * gr_xxy_xx[i] * gfe_0 + ts_xxyy_xx[i] * gfe_0 * gc_y[i] + gr_xxyy_xx[i] * gc_y[i];

        grr_y_xxyy_xy[i] = 2.0 * ts_xxy_xy[i] * gfe2_0 + 2.0 * gr_xxy_xy[i] * gfe_0 + ts_xxyy_x[i] * gfe2_0 + gr_xxyy_x[i] * gfe_0 + ts_xxyy_xy[i] * gfe_0 * gc_y[i] + gr_xxyy_xy[i] * gc_y[i];

        grr_y_xxyy_xz[i] = 2.0 * ts_xxy_xz[i] * gfe2_0 + 2.0 * gr_xxy_xz[i] * gfe_0 + ts_xxyy_xz[i] * gfe_0 * gc_y[i] + gr_xxyy_xz[i] * gc_y[i];

        grr_y_xxyy_yy[i] = 2.0 * ts_xxy_yy[i] * gfe2_0 + 2.0 * gr_xxy_yy[i] * gfe_0 + 2.0 * ts_xxyy_y[i] * gfe2_0 + 2.0 * gr_xxyy_y[i] * gfe_0 + ts_xxyy_yy[i] * gfe_0 * gc_y[i] + gr_xxyy_yy[i] * gc_y[i];

        grr_y_xxyy_yz[i] = 2.0 * ts_xxy_yz[i] * gfe2_0 + 2.0 * gr_xxy_yz[i] * gfe_0 + ts_xxyy_z[i] * gfe2_0 + gr_xxyy_z[i] * gfe_0 + ts_xxyy_yz[i] * gfe_0 * gc_y[i] + gr_xxyy_yz[i] * gc_y[i];

        grr_y_xxyy_zz[i] = 2.0 * ts_xxy_zz[i] * gfe2_0 + 2.0 * gr_xxy_zz[i] * gfe_0 + ts_xxyy_zz[i] * gfe_0 * gc_y[i] + gr_xxyy_zz[i] * gc_y[i];
    }

    // Set up 114-120 components of targeted buffer : GD

    auto grr_y_xxyz_xx = pbuffer.data(idx_gr_gd + 114);

    auto grr_y_xxyz_xy = pbuffer.data(idx_gr_gd + 115);

    auto grr_y_xxyz_xz = pbuffer.data(idx_gr_gd + 116);

    auto grr_y_xxyz_yy = pbuffer.data(idx_gr_gd + 117);

    auto grr_y_xxyz_yz = pbuffer.data(idx_gr_gd + 118);

    auto grr_y_xxyz_zz = pbuffer.data(idx_gr_gd + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_x, gr_xxyz_xx, gr_xxyz_xy, gr_xxyz_xz, gr_xxyz_y, gr_xxyz_yy, gr_xxyz_yz, gr_xxyz_z, gr_xxyz_zz, gr_xxz_xx, gr_xxz_xy, gr_xxz_xz, gr_xxz_yy, gr_xxz_yz, gr_xxz_zz, grr_y_xxyz_xx, grr_y_xxyz_xy, grr_y_xxyz_xz, grr_y_xxyz_yy, grr_y_xxyz_yz, grr_y_xxyz_zz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_yy, ts_xxz_yz, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxyz_xx[i] = ts_xxz_xx[i] * gfe2_0 + gr_xxz_xx[i] * gfe_0 + ts_xxyz_xx[i] * gfe_0 * gc_y[i] + gr_xxyz_xx[i] * gc_y[i];

        grr_y_xxyz_xy[i] = ts_xxz_xy[i] * gfe2_0 + gr_xxz_xy[i] * gfe_0 + ts_xxyz_x[i] * gfe2_0 + gr_xxyz_x[i] * gfe_0 + ts_xxyz_xy[i] * gfe_0 * gc_y[i] + gr_xxyz_xy[i] * gc_y[i];

        grr_y_xxyz_xz[i] = ts_xxz_xz[i] * gfe2_0 + gr_xxz_xz[i] * gfe_0 + ts_xxyz_xz[i] * gfe_0 * gc_y[i] + gr_xxyz_xz[i] * gc_y[i];

        grr_y_xxyz_yy[i] = ts_xxz_yy[i] * gfe2_0 + gr_xxz_yy[i] * gfe_0 + 2.0 * ts_xxyz_y[i] * gfe2_0 + 2.0 * gr_xxyz_y[i] * gfe_0 + ts_xxyz_yy[i] * gfe_0 * gc_y[i] + gr_xxyz_yy[i] * gc_y[i];

        grr_y_xxyz_yz[i] = ts_xxz_yz[i] * gfe2_0 + gr_xxz_yz[i] * gfe_0 + ts_xxyz_z[i] * gfe2_0 + gr_xxyz_z[i] * gfe_0 + ts_xxyz_yz[i] * gfe_0 * gc_y[i] + gr_xxyz_yz[i] * gc_y[i];

        grr_y_xxyz_zz[i] = ts_xxz_zz[i] * gfe2_0 + gr_xxz_zz[i] * gfe_0 + ts_xxyz_zz[i] * gfe_0 * gc_y[i] + gr_xxyz_zz[i] * gc_y[i];
    }

    // Set up 120-126 components of targeted buffer : GD

    auto grr_y_xxzz_xx = pbuffer.data(idx_gr_gd + 120);

    auto grr_y_xxzz_xy = pbuffer.data(idx_gr_gd + 121);

    auto grr_y_xxzz_xz = pbuffer.data(idx_gr_gd + 122);

    auto grr_y_xxzz_yy = pbuffer.data(idx_gr_gd + 123);

    auto grr_y_xxzz_yz = pbuffer.data(idx_gr_gd + 124);

    auto grr_y_xxzz_zz = pbuffer.data(idx_gr_gd + 125);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_x, gr_xxzz_xx, gr_xxzz_xy, gr_xxzz_xz, gr_xxzz_y, gr_xxzz_yy, gr_xxzz_yz, gr_xxzz_z, gr_xxzz_zz, grr_y_xxzz_xx, grr_y_xxzz_xy, grr_y_xxzz_xz, grr_y_xxzz_yy, grr_y_xxzz_yz, grr_y_xxzz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxzz_xx[i] = ts_xxzz_xx[i] * gfe_0 * gc_y[i] + gr_xxzz_xx[i] * gc_y[i];

        grr_y_xxzz_xy[i] = ts_xxzz_x[i] * gfe2_0 + gr_xxzz_x[i] * gfe_0 + ts_xxzz_xy[i] * gfe_0 * gc_y[i] + gr_xxzz_xy[i] * gc_y[i];

        grr_y_xxzz_xz[i] = ts_xxzz_xz[i] * gfe_0 * gc_y[i] + gr_xxzz_xz[i] * gc_y[i];

        grr_y_xxzz_yy[i] = 2.0 * ts_xxzz_y[i] * gfe2_0 + 2.0 * gr_xxzz_y[i] * gfe_0 + ts_xxzz_yy[i] * gfe_0 * gc_y[i] + gr_xxzz_yy[i] * gc_y[i];

        grr_y_xxzz_yz[i] = ts_xxzz_z[i] * gfe2_0 + gr_xxzz_z[i] * gfe_0 + ts_xxzz_yz[i] * gfe_0 * gc_y[i] + gr_xxzz_yz[i] * gc_y[i];

        grr_y_xxzz_zz[i] = ts_xxzz_zz[i] * gfe_0 * gc_y[i] + gr_xxzz_zz[i] * gc_y[i];
    }

    // Set up 126-132 components of targeted buffer : GD

    auto grr_y_xyyy_xx = pbuffer.data(idx_gr_gd + 126);

    auto grr_y_xyyy_xy = pbuffer.data(idx_gr_gd + 127);

    auto grr_y_xyyy_xz = pbuffer.data(idx_gr_gd + 128);

    auto grr_y_xyyy_yy = pbuffer.data(idx_gr_gd + 129);

    auto grr_y_xyyy_yz = pbuffer.data(idx_gr_gd + 130);

    auto grr_y_xyyy_zz = pbuffer.data(idx_gr_gd + 131);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xx, gr_xyy_xy, gr_xyy_xz, gr_xyy_yy, gr_xyy_yz, gr_xyy_zz, gr_xyyy_x, gr_xyyy_xx, gr_xyyy_xy, gr_xyyy_xz, gr_xyyy_y, gr_xyyy_yy, gr_xyyy_yz, gr_xyyy_z, gr_xyyy_zz, grr_y_xyyy_xx, grr_y_xyyy_xy, grr_y_xyyy_xz, grr_y_xyyy_yy, grr_y_xyyy_yz, grr_y_xyyy_zz, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyyy_xx[i] = 3.0 * ts_xyy_xx[i] * gfe2_0 + 3.0 * gr_xyy_xx[i] * gfe_0 + ts_xyyy_xx[i] * gfe_0 * gc_y[i] + gr_xyyy_xx[i] * gc_y[i];

        grr_y_xyyy_xy[i] = 3.0 * ts_xyy_xy[i] * gfe2_0 + 3.0 * gr_xyy_xy[i] * gfe_0 + ts_xyyy_x[i] * gfe2_0 + gr_xyyy_x[i] * gfe_0 + ts_xyyy_xy[i] * gfe_0 * gc_y[i] + gr_xyyy_xy[i] * gc_y[i];

        grr_y_xyyy_xz[i] = 3.0 * ts_xyy_xz[i] * gfe2_0 + 3.0 * gr_xyy_xz[i] * gfe_0 + ts_xyyy_xz[i] * gfe_0 * gc_y[i] + gr_xyyy_xz[i] * gc_y[i];

        grr_y_xyyy_yy[i] = 3.0 * ts_xyy_yy[i] * gfe2_0 + 3.0 * gr_xyy_yy[i] * gfe_0 + 2.0 * ts_xyyy_y[i] * gfe2_0 + 2.0 * gr_xyyy_y[i] * gfe_0 + ts_xyyy_yy[i] * gfe_0 * gc_y[i] + gr_xyyy_yy[i] * gc_y[i];

        grr_y_xyyy_yz[i] = 3.0 * ts_xyy_yz[i] * gfe2_0 + 3.0 * gr_xyy_yz[i] * gfe_0 + ts_xyyy_z[i] * gfe2_0 + gr_xyyy_z[i] * gfe_0 + ts_xyyy_yz[i] * gfe_0 * gc_y[i] + gr_xyyy_yz[i] * gc_y[i];

        grr_y_xyyy_zz[i] = 3.0 * ts_xyy_zz[i] * gfe2_0 + 3.0 * gr_xyy_zz[i] * gfe_0 + ts_xyyy_zz[i] * gfe_0 * gc_y[i] + gr_xyyy_zz[i] * gc_y[i];
    }

    // Set up 132-138 components of targeted buffer : GD

    auto grr_y_xyyz_xx = pbuffer.data(idx_gr_gd + 132);

    auto grr_y_xyyz_xy = pbuffer.data(idx_gr_gd + 133);

    auto grr_y_xyyz_xz = pbuffer.data(idx_gr_gd + 134);

    auto grr_y_xyyz_yy = pbuffer.data(idx_gr_gd + 135);

    auto grr_y_xyyz_yz = pbuffer.data(idx_gr_gd + 136);

    auto grr_y_xyyz_zz = pbuffer.data(idx_gr_gd + 137);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_x, gr_xyyz_xx, gr_xyyz_xy, gr_xyyz_xz, gr_xyyz_y, gr_xyyz_yy, gr_xyyz_yz, gr_xyyz_z, gr_xyyz_zz, gr_xyz_xx, gr_xyz_xy, gr_xyz_xz, gr_xyz_yy, gr_xyz_yz, gr_xyz_zz, grr_y_xyyz_xx, grr_y_xyyz_xy, grr_y_xyyz_xz, grr_y_xyyz_yy, grr_y_xyyz_yz, grr_y_xyyz_zz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_yy, ts_xyz_yz, ts_xyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyyz_xx[i] = 2.0 * ts_xyz_xx[i] * gfe2_0 + 2.0 * gr_xyz_xx[i] * gfe_0 + ts_xyyz_xx[i] * gfe_0 * gc_y[i] + gr_xyyz_xx[i] * gc_y[i];

        grr_y_xyyz_xy[i] = 2.0 * ts_xyz_xy[i] * gfe2_0 + 2.0 * gr_xyz_xy[i] * gfe_0 + ts_xyyz_x[i] * gfe2_0 + gr_xyyz_x[i] * gfe_0 + ts_xyyz_xy[i] * gfe_0 * gc_y[i] + gr_xyyz_xy[i] * gc_y[i];

        grr_y_xyyz_xz[i] = 2.0 * ts_xyz_xz[i] * gfe2_0 + 2.0 * gr_xyz_xz[i] * gfe_0 + ts_xyyz_xz[i] * gfe_0 * gc_y[i] + gr_xyyz_xz[i] * gc_y[i];

        grr_y_xyyz_yy[i] = 2.0 * ts_xyz_yy[i] * gfe2_0 + 2.0 * gr_xyz_yy[i] * gfe_0 + 2.0 * ts_xyyz_y[i] * gfe2_0 + 2.0 * gr_xyyz_y[i] * gfe_0 + ts_xyyz_yy[i] * gfe_0 * gc_y[i] + gr_xyyz_yy[i] * gc_y[i];

        grr_y_xyyz_yz[i] = 2.0 * ts_xyz_yz[i] * gfe2_0 + 2.0 * gr_xyz_yz[i] * gfe_0 + ts_xyyz_z[i] * gfe2_0 + gr_xyyz_z[i] * gfe_0 + ts_xyyz_yz[i] * gfe_0 * gc_y[i] + gr_xyyz_yz[i] * gc_y[i];

        grr_y_xyyz_zz[i] = 2.0 * ts_xyz_zz[i] * gfe2_0 + 2.0 * gr_xyz_zz[i] * gfe_0 + ts_xyyz_zz[i] * gfe_0 * gc_y[i] + gr_xyyz_zz[i] * gc_y[i];
    }

    // Set up 138-144 components of targeted buffer : GD

    auto grr_y_xyzz_xx = pbuffer.data(idx_gr_gd + 138);

    auto grr_y_xyzz_xy = pbuffer.data(idx_gr_gd + 139);

    auto grr_y_xyzz_xz = pbuffer.data(idx_gr_gd + 140);

    auto grr_y_xyzz_yy = pbuffer.data(idx_gr_gd + 141);

    auto grr_y_xyzz_yz = pbuffer.data(idx_gr_gd + 142);

    auto grr_y_xyzz_zz = pbuffer.data(idx_gr_gd + 143);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_x, gr_xyzz_xx, gr_xyzz_xy, gr_xyzz_xz, gr_xyzz_y, gr_xyzz_yy, gr_xyzz_yz, gr_xyzz_z, gr_xyzz_zz, gr_xzz_xx, gr_xzz_xy, gr_xzz_xz, gr_xzz_yy, gr_xzz_yz, gr_xzz_zz, grr_y_xyzz_xx, grr_y_xyzz_xy, grr_y_xyzz_xz, grr_y_xyzz_yy, grr_y_xyzz_yz, grr_y_xyzz_zz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyzz_xx[i] = ts_xzz_xx[i] * gfe2_0 + gr_xzz_xx[i] * gfe_0 + ts_xyzz_xx[i] * gfe_0 * gc_y[i] + gr_xyzz_xx[i] * gc_y[i];

        grr_y_xyzz_xy[i] = ts_xzz_xy[i] * gfe2_0 + gr_xzz_xy[i] * gfe_0 + ts_xyzz_x[i] * gfe2_0 + gr_xyzz_x[i] * gfe_0 + ts_xyzz_xy[i] * gfe_0 * gc_y[i] + gr_xyzz_xy[i] * gc_y[i];

        grr_y_xyzz_xz[i] = ts_xzz_xz[i] * gfe2_0 + gr_xzz_xz[i] * gfe_0 + ts_xyzz_xz[i] * gfe_0 * gc_y[i] + gr_xyzz_xz[i] * gc_y[i];

        grr_y_xyzz_yy[i] = ts_xzz_yy[i] * gfe2_0 + gr_xzz_yy[i] * gfe_0 + 2.0 * ts_xyzz_y[i] * gfe2_0 + 2.0 * gr_xyzz_y[i] * gfe_0 + ts_xyzz_yy[i] * gfe_0 * gc_y[i] + gr_xyzz_yy[i] * gc_y[i];

        grr_y_xyzz_yz[i] = ts_xzz_yz[i] * gfe2_0 + gr_xzz_yz[i] * gfe_0 + ts_xyzz_z[i] * gfe2_0 + gr_xyzz_z[i] * gfe_0 + ts_xyzz_yz[i] * gfe_0 * gc_y[i] + gr_xyzz_yz[i] * gc_y[i];

        grr_y_xyzz_zz[i] = ts_xzz_zz[i] * gfe2_0 + gr_xzz_zz[i] * gfe_0 + ts_xyzz_zz[i] * gfe_0 * gc_y[i] + gr_xyzz_zz[i] * gc_y[i];
    }

    // Set up 144-150 components of targeted buffer : GD

    auto grr_y_xzzz_xx = pbuffer.data(idx_gr_gd + 144);

    auto grr_y_xzzz_xy = pbuffer.data(idx_gr_gd + 145);

    auto grr_y_xzzz_xz = pbuffer.data(idx_gr_gd + 146);

    auto grr_y_xzzz_yy = pbuffer.data(idx_gr_gd + 147);

    auto grr_y_xzzz_yz = pbuffer.data(idx_gr_gd + 148);

    auto grr_y_xzzz_zz = pbuffer.data(idx_gr_gd + 149);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_x, gr_xzzz_xx, gr_xzzz_xy, gr_xzzz_xz, gr_xzzz_y, gr_xzzz_yy, gr_xzzz_yz, gr_xzzz_z, gr_xzzz_zz, grr_y_xzzz_xx, grr_y_xzzz_xy, grr_y_xzzz_xz, grr_y_xzzz_yy, grr_y_xzzz_yz, grr_y_xzzz_zz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xzzz_xx[i] = ts_xzzz_xx[i] * gfe_0 * gc_y[i] + gr_xzzz_xx[i] * gc_y[i];

        grr_y_xzzz_xy[i] = ts_xzzz_x[i] * gfe2_0 + gr_xzzz_x[i] * gfe_0 + ts_xzzz_xy[i] * gfe_0 * gc_y[i] + gr_xzzz_xy[i] * gc_y[i];

        grr_y_xzzz_xz[i] = ts_xzzz_xz[i] * gfe_0 * gc_y[i] + gr_xzzz_xz[i] * gc_y[i];

        grr_y_xzzz_yy[i] = 2.0 * ts_xzzz_y[i] * gfe2_0 + 2.0 * gr_xzzz_y[i] * gfe_0 + ts_xzzz_yy[i] * gfe_0 * gc_y[i] + gr_xzzz_yy[i] * gc_y[i];

        grr_y_xzzz_yz[i] = ts_xzzz_z[i] * gfe2_0 + gr_xzzz_z[i] * gfe_0 + ts_xzzz_yz[i] * gfe_0 * gc_y[i] + gr_xzzz_yz[i] * gc_y[i];

        grr_y_xzzz_zz[i] = ts_xzzz_zz[i] * gfe_0 * gc_y[i] + gr_xzzz_zz[i] * gc_y[i];
    }

    // Set up 150-156 components of targeted buffer : GD

    auto grr_y_yyyy_xx = pbuffer.data(idx_gr_gd + 150);

    auto grr_y_yyyy_xy = pbuffer.data(idx_gr_gd + 151);

    auto grr_y_yyyy_xz = pbuffer.data(idx_gr_gd + 152);

    auto grr_y_yyyy_yy = pbuffer.data(idx_gr_gd + 153);

    auto grr_y_yyyy_yz = pbuffer.data(idx_gr_gd + 154);

    auto grr_y_yyyy_zz = pbuffer.data(idx_gr_gd + 155);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xx, gr_yyy_xy, gr_yyy_xz, gr_yyy_yy, gr_yyy_yz, gr_yyy_zz, gr_yyyy_x, gr_yyyy_xx, gr_yyyy_xy, gr_yyyy_xz, gr_yyyy_y, gr_yyyy_yy, gr_yyyy_yz, gr_yyyy_z, gr_yyyy_zz, grr_y_yyyy_xx, grr_y_yyyy_xy, grr_y_yyyy_xz, grr_y_yyyy_yy, grr_y_yyyy_yz, grr_y_yyyy_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyyy_xx[i] = 4.0 * ts_yyy_xx[i] * gfe2_0 + 4.0 * gr_yyy_xx[i] * gfe_0 + ts_yyyy_xx[i] * gfe_0 * gc_y[i] + gr_yyyy_xx[i] * gc_y[i];

        grr_y_yyyy_xy[i] = 4.0 * ts_yyy_xy[i] * gfe2_0 + 4.0 * gr_yyy_xy[i] * gfe_0 + ts_yyyy_x[i] * gfe2_0 + gr_yyyy_x[i] * gfe_0 + ts_yyyy_xy[i] * gfe_0 * gc_y[i] + gr_yyyy_xy[i] * gc_y[i];

        grr_y_yyyy_xz[i] = 4.0 * ts_yyy_xz[i] * gfe2_0 + 4.0 * gr_yyy_xz[i] * gfe_0 + ts_yyyy_xz[i] * gfe_0 * gc_y[i] + gr_yyyy_xz[i] * gc_y[i];

        grr_y_yyyy_yy[i] = 4.0 * ts_yyy_yy[i] * gfe2_0 + 4.0 * gr_yyy_yy[i] * gfe_0 + 2.0 * ts_yyyy_y[i] * gfe2_0 + 2.0 * gr_yyyy_y[i] * gfe_0 + ts_yyyy_yy[i] * gfe_0 * gc_y[i] + gr_yyyy_yy[i] * gc_y[i];

        grr_y_yyyy_yz[i] = 4.0 * ts_yyy_yz[i] * gfe2_0 + 4.0 * gr_yyy_yz[i] * gfe_0 + ts_yyyy_z[i] * gfe2_0 + gr_yyyy_z[i] * gfe_0 + ts_yyyy_yz[i] * gfe_0 * gc_y[i] + gr_yyyy_yz[i] * gc_y[i];

        grr_y_yyyy_zz[i] = 4.0 * ts_yyy_zz[i] * gfe2_0 + 4.0 * gr_yyy_zz[i] * gfe_0 + ts_yyyy_zz[i] * gfe_0 * gc_y[i] + gr_yyyy_zz[i] * gc_y[i];
    }

    // Set up 156-162 components of targeted buffer : GD

    auto grr_y_yyyz_xx = pbuffer.data(idx_gr_gd + 156);

    auto grr_y_yyyz_xy = pbuffer.data(idx_gr_gd + 157);

    auto grr_y_yyyz_xz = pbuffer.data(idx_gr_gd + 158);

    auto grr_y_yyyz_yy = pbuffer.data(idx_gr_gd + 159);

    auto grr_y_yyyz_yz = pbuffer.data(idx_gr_gd + 160);

    auto grr_y_yyyz_zz = pbuffer.data(idx_gr_gd + 161);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_x, gr_yyyz_xx, gr_yyyz_xy, gr_yyyz_xz, gr_yyyz_y, gr_yyyz_yy, gr_yyyz_yz, gr_yyyz_z, gr_yyyz_zz, gr_yyz_xx, gr_yyz_xy, gr_yyz_xz, gr_yyz_yy, gr_yyz_yz, gr_yyz_zz, grr_y_yyyz_xx, grr_y_yyyz_xy, grr_y_yyyz_xz, grr_y_yyyz_yy, grr_y_yyyz_yz, grr_y_yyyz_zz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_yy, ts_yyz_yz, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyyz_xx[i] = 3.0 * ts_yyz_xx[i] * gfe2_0 + 3.0 * gr_yyz_xx[i] * gfe_0 + ts_yyyz_xx[i] * gfe_0 * gc_y[i] + gr_yyyz_xx[i] * gc_y[i];

        grr_y_yyyz_xy[i] = 3.0 * ts_yyz_xy[i] * gfe2_0 + 3.0 * gr_yyz_xy[i] * gfe_0 + ts_yyyz_x[i] * gfe2_0 + gr_yyyz_x[i] * gfe_0 + ts_yyyz_xy[i] * gfe_0 * gc_y[i] + gr_yyyz_xy[i] * gc_y[i];

        grr_y_yyyz_xz[i] = 3.0 * ts_yyz_xz[i] * gfe2_0 + 3.0 * gr_yyz_xz[i] * gfe_0 + ts_yyyz_xz[i] * gfe_0 * gc_y[i] + gr_yyyz_xz[i] * gc_y[i];

        grr_y_yyyz_yy[i] = 3.0 * ts_yyz_yy[i] * gfe2_0 + 3.0 * gr_yyz_yy[i] * gfe_0 + 2.0 * ts_yyyz_y[i] * gfe2_0 + 2.0 * gr_yyyz_y[i] * gfe_0 + ts_yyyz_yy[i] * gfe_0 * gc_y[i] + gr_yyyz_yy[i] * gc_y[i];

        grr_y_yyyz_yz[i] = 3.0 * ts_yyz_yz[i] * gfe2_0 + 3.0 * gr_yyz_yz[i] * gfe_0 + ts_yyyz_z[i] * gfe2_0 + gr_yyyz_z[i] * gfe_0 + ts_yyyz_yz[i] * gfe_0 * gc_y[i] + gr_yyyz_yz[i] * gc_y[i];

        grr_y_yyyz_zz[i] = 3.0 * ts_yyz_zz[i] * gfe2_0 + 3.0 * gr_yyz_zz[i] * gfe_0 + ts_yyyz_zz[i] * gfe_0 * gc_y[i] + gr_yyyz_zz[i] * gc_y[i];
    }

    // Set up 162-168 components of targeted buffer : GD

    auto grr_y_yyzz_xx = pbuffer.data(idx_gr_gd + 162);

    auto grr_y_yyzz_xy = pbuffer.data(idx_gr_gd + 163);

    auto grr_y_yyzz_xz = pbuffer.data(idx_gr_gd + 164);

    auto grr_y_yyzz_yy = pbuffer.data(idx_gr_gd + 165);

    auto grr_y_yyzz_yz = pbuffer.data(idx_gr_gd + 166);

    auto grr_y_yyzz_zz = pbuffer.data(idx_gr_gd + 167);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_x, gr_yyzz_xx, gr_yyzz_xy, gr_yyzz_xz, gr_yyzz_y, gr_yyzz_yy, gr_yyzz_yz, gr_yyzz_z, gr_yyzz_zz, gr_yzz_xx, gr_yzz_xy, gr_yzz_xz, gr_yzz_yy, gr_yzz_yz, gr_yzz_zz, grr_y_yyzz_xx, grr_y_yyzz_xy, grr_y_yyzz_xz, grr_y_yyzz_yy, grr_y_yyzz_yz, grr_y_yyzz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_yy, ts_yzz_yz, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyzz_xx[i] = 2.0 * ts_yzz_xx[i] * gfe2_0 + 2.0 * gr_yzz_xx[i] * gfe_0 + ts_yyzz_xx[i] * gfe_0 * gc_y[i] + gr_yyzz_xx[i] * gc_y[i];

        grr_y_yyzz_xy[i] = 2.0 * ts_yzz_xy[i] * gfe2_0 + 2.0 * gr_yzz_xy[i] * gfe_0 + ts_yyzz_x[i] * gfe2_0 + gr_yyzz_x[i] * gfe_0 + ts_yyzz_xy[i] * gfe_0 * gc_y[i] + gr_yyzz_xy[i] * gc_y[i];

        grr_y_yyzz_xz[i] = 2.0 * ts_yzz_xz[i] * gfe2_0 + 2.0 * gr_yzz_xz[i] * gfe_0 + ts_yyzz_xz[i] * gfe_0 * gc_y[i] + gr_yyzz_xz[i] * gc_y[i];

        grr_y_yyzz_yy[i] = 2.0 * ts_yzz_yy[i] * gfe2_0 + 2.0 * gr_yzz_yy[i] * gfe_0 + 2.0 * ts_yyzz_y[i] * gfe2_0 + 2.0 * gr_yyzz_y[i] * gfe_0 + ts_yyzz_yy[i] * gfe_0 * gc_y[i] + gr_yyzz_yy[i] * gc_y[i];

        grr_y_yyzz_yz[i] = 2.0 * ts_yzz_yz[i] * gfe2_0 + 2.0 * gr_yzz_yz[i] * gfe_0 + ts_yyzz_z[i] * gfe2_0 + gr_yyzz_z[i] * gfe_0 + ts_yyzz_yz[i] * gfe_0 * gc_y[i] + gr_yyzz_yz[i] * gc_y[i];

        grr_y_yyzz_zz[i] = 2.0 * ts_yzz_zz[i] * gfe2_0 + 2.0 * gr_yzz_zz[i] * gfe_0 + ts_yyzz_zz[i] * gfe_0 * gc_y[i] + gr_yyzz_zz[i] * gc_y[i];
    }

    // Set up 168-174 components of targeted buffer : GD

    auto grr_y_yzzz_xx = pbuffer.data(idx_gr_gd + 168);

    auto grr_y_yzzz_xy = pbuffer.data(idx_gr_gd + 169);

    auto grr_y_yzzz_xz = pbuffer.data(idx_gr_gd + 170);

    auto grr_y_yzzz_yy = pbuffer.data(idx_gr_gd + 171);

    auto grr_y_yzzz_yz = pbuffer.data(idx_gr_gd + 172);

    auto grr_y_yzzz_zz = pbuffer.data(idx_gr_gd + 173);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_x, gr_yzzz_xx, gr_yzzz_xy, gr_yzzz_xz, gr_yzzz_y, gr_yzzz_yy, gr_yzzz_yz, gr_yzzz_z, gr_yzzz_zz, gr_zzz_xx, gr_zzz_xy, gr_zzz_xz, gr_zzz_yy, gr_zzz_yz, gr_zzz_zz, grr_y_yzzz_xx, grr_y_yzzz_xy, grr_y_yzzz_xz, grr_y_yzzz_yy, grr_y_yzzz_yz, grr_y_yzzz_zz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yzzz_xx[i] = ts_zzz_xx[i] * gfe2_0 + gr_zzz_xx[i] * gfe_0 + ts_yzzz_xx[i] * gfe_0 * gc_y[i] + gr_yzzz_xx[i] * gc_y[i];

        grr_y_yzzz_xy[i] = ts_zzz_xy[i] * gfe2_0 + gr_zzz_xy[i] * gfe_0 + ts_yzzz_x[i] * gfe2_0 + gr_yzzz_x[i] * gfe_0 + ts_yzzz_xy[i] * gfe_0 * gc_y[i] + gr_yzzz_xy[i] * gc_y[i];

        grr_y_yzzz_xz[i] = ts_zzz_xz[i] * gfe2_0 + gr_zzz_xz[i] * gfe_0 + ts_yzzz_xz[i] * gfe_0 * gc_y[i] + gr_yzzz_xz[i] * gc_y[i];

        grr_y_yzzz_yy[i] = ts_zzz_yy[i] * gfe2_0 + gr_zzz_yy[i] * gfe_0 + 2.0 * ts_yzzz_y[i] * gfe2_0 + 2.0 * gr_yzzz_y[i] * gfe_0 + ts_yzzz_yy[i] * gfe_0 * gc_y[i] + gr_yzzz_yy[i] * gc_y[i];

        grr_y_yzzz_yz[i] = ts_zzz_yz[i] * gfe2_0 + gr_zzz_yz[i] * gfe_0 + ts_yzzz_z[i] * gfe2_0 + gr_yzzz_z[i] * gfe_0 + ts_yzzz_yz[i] * gfe_0 * gc_y[i] + gr_yzzz_yz[i] * gc_y[i];

        grr_y_yzzz_zz[i] = ts_zzz_zz[i] * gfe2_0 + gr_zzz_zz[i] * gfe_0 + ts_yzzz_zz[i] * gfe_0 * gc_y[i] + gr_yzzz_zz[i] * gc_y[i];
    }

    // Set up 174-180 components of targeted buffer : GD

    auto grr_y_zzzz_xx = pbuffer.data(idx_gr_gd + 174);

    auto grr_y_zzzz_xy = pbuffer.data(idx_gr_gd + 175);

    auto grr_y_zzzz_xz = pbuffer.data(idx_gr_gd + 176);

    auto grr_y_zzzz_yy = pbuffer.data(idx_gr_gd + 177);

    auto grr_y_zzzz_yz = pbuffer.data(idx_gr_gd + 178);

    auto grr_y_zzzz_zz = pbuffer.data(idx_gr_gd + 179);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_x, gr_zzzz_xx, gr_zzzz_xy, gr_zzzz_xz, gr_zzzz_y, gr_zzzz_yy, gr_zzzz_yz, gr_zzzz_z, gr_zzzz_zz, grr_y_zzzz_xx, grr_y_zzzz_xy, grr_y_zzzz_xz, grr_y_zzzz_yy, grr_y_zzzz_yz, grr_y_zzzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zzzz_xx[i] = ts_zzzz_xx[i] * gfe_0 * gc_y[i] + gr_zzzz_xx[i] * gc_y[i];

        grr_y_zzzz_xy[i] = ts_zzzz_x[i] * gfe2_0 + gr_zzzz_x[i] * gfe_0 + ts_zzzz_xy[i] * gfe_0 * gc_y[i] + gr_zzzz_xy[i] * gc_y[i];

        grr_y_zzzz_xz[i] = ts_zzzz_xz[i] * gfe_0 * gc_y[i] + gr_zzzz_xz[i] * gc_y[i];

        grr_y_zzzz_yy[i] = 2.0 * ts_zzzz_y[i] * gfe2_0 + 2.0 * gr_zzzz_y[i] * gfe_0 + ts_zzzz_yy[i] * gfe_0 * gc_y[i] + gr_zzzz_yy[i] * gc_y[i];

        grr_y_zzzz_yz[i] = ts_zzzz_z[i] * gfe2_0 + gr_zzzz_z[i] * gfe_0 + ts_zzzz_yz[i] * gfe_0 * gc_y[i] + gr_zzzz_yz[i] * gc_y[i];

        grr_y_zzzz_zz[i] = ts_zzzz_zz[i] * gfe_0 * gc_y[i] + gr_zzzz_zz[i] * gc_y[i];
    }

    // Set up 180-186 components of targeted buffer : GD

    auto grr_z_xxxx_xx = pbuffer.data(idx_gr_gd + 180);

    auto grr_z_xxxx_xy = pbuffer.data(idx_gr_gd + 181);

    auto grr_z_xxxx_xz = pbuffer.data(idx_gr_gd + 182);

    auto grr_z_xxxx_yy = pbuffer.data(idx_gr_gd + 183);

    auto grr_z_xxxx_yz = pbuffer.data(idx_gr_gd + 184);

    auto grr_z_xxxx_zz = pbuffer.data(idx_gr_gd + 185);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_x, gr_xxxx_xx, gr_xxxx_xy, gr_xxxx_xz, gr_xxxx_y, gr_xxxx_yy, gr_xxxx_yz, gr_xxxx_z, gr_xxxx_zz, grr_z_xxxx_xx, grr_z_xxxx_xy, grr_z_xxxx_xz, grr_z_xxxx_yy, grr_z_xxxx_yz, grr_z_xxxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxxx_xx[i] = ts_xxxx_xx[i] * gfe_0 * gc_z[i] + gr_xxxx_xx[i] * gc_z[i];

        grr_z_xxxx_xy[i] = ts_xxxx_xy[i] * gfe_0 * gc_z[i] + gr_xxxx_xy[i] * gc_z[i];

        grr_z_xxxx_xz[i] = ts_xxxx_x[i] * gfe2_0 + gr_xxxx_x[i] * gfe_0 + ts_xxxx_xz[i] * gfe_0 * gc_z[i] + gr_xxxx_xz[i] * gc_z[i];

        grr_z_xxxx_yy[i] = ts_xxxx_yy[i] * gfe_0 * gc_z[i] + gr_xxxx_yy[i] * gc_z[i];

        grr_z_xxxx_yz[i] = ts_xxxx_y[i] * gfe2_0 + gr_xxxx_y[i] * gfe_0 + ts_xxxx_yz[i] * gfe_0 * gc_z[i] + gr_xxxx_yz[i] * gc_z[i];

        grr_z_xxxx_zz[i] = 2.0 * ts_xxxx_z[i] * gfe2_0 + 2.0 * gr_xxxx_z[i] * gfe_0 + ts_xxxx_zz[i] * gfe_0 * gc_z[i] + gr_xxxx_zz[i] * gc_z[i];
    }

    // Set up 186-192 components of targeted buffer : GD

    auto grr_z_xxxy_xx = pbuffer.data(idx_gr_gd + 186);

    auto grr_z_xxxy_xy = pbuffer.data(idx_gr_gd + 187);

    auto grr_z_xxxy_xz = pbuffer.data(idx_gr_gd + 188);

    auto grr_z_xxxy_yy = pbuffer.data(idx_gr_gd + 189);

    auto grr_z_xxxy_yz = pbuffer.data(idx_gr_gd + 190);

    auto grr_z_xxxy_zz = pbuffer.data(idx_gr_gd + 191);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_x, gr_xxxy_xx, gr_xxxy_xy, gr_xxxy_xz, gr_xxxy_y, gr_xxxy_yy, gr_xxxy_yz, gr_xxxy_z, gr_xxxy_zz, grr_z_xxxy_xx, grr_z_xxxy_xy, grr_z_xxxy_xz, grr_z_xxxy_yy, grr_z_xxxy_yz, grr_z_xxxy_zz, ts_xxxy_x, ts_xxxy_xx, ts_xxxy_xy, ts_xxxy_xz, ts_xxxy_y, ts_xxxy_yy, ts_xxxy_yz, ts_xxxy_z, ts_xxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxxy_xx[i] = ts_xxxy_xx[i] * gfe_0 * gc_z[i] + gr_xxxy_xx[i] * gc_z[i];

        grr_z_xxxy_xy[i] = ts_xxxy_xy[i] * gfe_0 * gc_z[i] + gr_xxxy_xy[i] * gc_z[i];

        grr_z_xxxy_xz[i] = ts_xxxy_x[i] * gfe2_0 + gr_xxxy_x[i] * gfe_0 + ts_xxxy_xz[i] * gfe_0 * gc_z[i] + gr_xxxy_xz[i] * gc_z[i];

        grr_z_xxxy_yy[i] = ts_xxxy_yy[i] * gfe_0 * gc_z[i] + gr_xxxy_yy[i] * gc_z[i];

        grr_z_xxxy_yz[i] = ts_xxxy_y[i] * gfe2_0 + gr_xxxy_y[i] * gfe_0 + ts_xxxy_yz[i] * gfe_0 * gc_z[i] + gr_xxxy_yz[i] * gc_z[i];

        grr_z_xxxy_zz[i] = 2.0 * ts_xxxy_z[i] * gfe2_0 + 2.0 * gr_xxxy_z[i] * gfe_0 + ts_xxxy_zz[i] * gfe_0 * gc_z[i] + gr_xxxy_zz[i] * gc_z[i];
    }

    // Set up 192-198 components of targeted buffer : GD

    auto grr_z_xxxz_xx = pbuffer.data(idx_gr_gd + 192);

    auto grr_z_xxxz_xy = pbuffer.data(idx_gr_gd + 193);

    auto grr_z_xxxz_xz = pbuffer.data(idx_gr_gd + 194);

    auto grr_z_xxxz_yy = pbuffer.data(idx_gr_gd + 195);

    auto grr_z_xxxz_yz = pbuffer.data(idx_gr_gd + 196);

    auto grr_z_xxxz_zz = pbuffer.data(idx_gr_gd + 197);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_xx, gr_xxx_xy, gr_xxx_xz, gr_xxx_yy, gr_xxx_yz, gr_xxx_zz, gr_xxxz_x, gr_xxxz_xx, gr_xxxz_xy, gr_xxxz_xz, gr_xxxz_y, gr_xxxz_yy, gr_xxxz_yz, gr_xxxz_z, gr_xxxz_zz, grr_z_xxxz_xx, grr_z_xxxz_xy, grr_z_xxxz_xz, grr_z_xxxz_yy, grr_z_xxxz_yz, grr_z_xxxz_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxz_x, ts_xxxz_xx, ts_xxxz_xy, ts_xxxz_xz, ts_xxxz_y, ts_xxxz_yy, ts_xxxz_yz, ts_xxxz_z, ts_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxxz_xx[i] = ts_xxx_xx[i] * gfe2_0 + gr_xxx_xx[i] * gfe_0 + ts_xxxz_xx[i] * gfe_0 * gc_z[i] + gr_xxxz_xx[i] * gc_z[i];

        grr_z_xxxz_xy[i] = ts_xxx_xy[i] * gfe2_0 + gr_xxx_xy[i] * gfe_0 + ts_xxxz_xy[i] * gfe_0 * gc_z[i] + gr_xxxz_xy[i] * gc_z[i];

        grr_z_xxxz_xz[i] = ts_xxx_xz[i] * gfe2_0 + gr_xxx_xz[i] * gfe_0 + ts_xxxz_x[i] * gfe2_0 + gr_xxxz_x[i] * gfe_0 + ts_xxxz_xz[i] * gfe_0 * gc_z[i] + gr_xxxz_xz[i] * gc_z[i];

        grr_z_xxxz_yy[i] = ts_xxx_yy[i] * gfe2_0 + gr_xxx_yy[i] * gfe_0 + ts_xxxz_yy[i] * gfe_0 * gc_z[i] + gr_xxxz_yy[i] * gc_z[i];

        grr_z_xxxz_yz[i] = ts_xxx_yz[i] * gfe2_0 + gr_xxx_yz[i] * gfe_0 + ts_xxxz_y[i] * gfe2_0 + gr_xxxz_y[i] * gfe_0 + ts_xxxz_yz[i] * gfe_0 * gc_z[i] + gr_xxxz_yz[i] * gc_z[i];

        grr_z_xxxz_zz[i] = ts_xxx_zz[i] * gfe2_0 + gr_xxx_zz[i] * gfe_0 + 2.0 * ts_xxxz_z[i] * gfe2_0 + 2.0 * gr_xxxz_z[i] * gfe_0 + ts_xxxz_zz[i] * gfe_0 * gc_z[i] + gr_xxxz_zz[i] * gc_z[i];
    }

    // Set up 198-204 components of targeted buffer : GD

    auto grr_z_xxyy_xx = pbuffer.data(idx_gr_gd + 198);

    auto grr_z_xxyy_xy = pbuffer.data(idx_gr_gd + 199);

    auto grr_z_xxyy_xz = pbuffer.data(idx_gr_gd + 200);

    auto grr_z_xxyy_yy = pbuffer.data(idx_gr_gd + 201);

    auto grr_z_xxyy_yz = pbuffer.data(idx_gr_gd + 202);

    auto grr_z_xxyy_zz = pbuffer.data(idx_gr_gd + 203);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_x, gr_xxyy_xx, gr_xxyy_xy, gr_xxyy_xz, gr_xxyy_y, gr_xxyy_yy, gr_xxyy_yz, gr_xxyy_z, gr_xxyy_zz, grr_z_xxyy_xx, grr_z_xxyy_xy, grr_z_xxyy_xz, grr_z_xxyy_yy, grr_z_xxyy_yz, grr_z_xxyy_zz, ts_xxyy_x, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_xz, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_z, ts_xxyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxyy_xx[i] = ts_xxyy_xx[i] * gfe_0 * gc_z[i] + gr_xxyy_xx[i] * gc_z[i];

        grr_z_xxyy_xy[i] = ts_xxyy_xy[i] * gfe_0 * gc_z[i] + gr_xxyy_xy[i] * gc_z[i];

        grr_z_xxyy_xz[i] = ts_xxyy_x[i] * gfe2_0 + gr_xxyy_x[i] * gfe_0 + ts_xxyy_xz[i] * gfe_0 * gc_z[i] + gr_xxyy_xz[i] * gc_z[i];

        grr_z_xxyy_yy[i] = ts_xxyy_yy[i] * gfe_0 * gc_z[i] + gr_xxyy_yy[i] * gc_z[i];

        grr_z_xxyy_yz[i] = ts_xxyy_y[i] * gfe2_0 + gr_xxyy_y[i] * gfe_0 + ts_xxyy_yz[i] * gfe_0 * gc_z[i] + gr_xxyy_yz[i] * gc_z[i];

        grr_z_xxyy_zz[i] = 2.0 * ts_xxyy_z[i] * gfe2_0 + 2.0 * gr_xxyy_z[i] * gfe_0 + ts_xxyy_zz[i] * gfe_0 * gc_z[i] + gr_xxyy_zz[i] * gc_z[i];
    }

    // Set up 204-210 components of targeted buffer : GD

    auto grr_z_xxyz_xx = pbuffer.data(idx_gr_gd + 204);

    auto grr_z_xxyz_xy = pbuffer.data(idx_gr_gd + 205);

    auto grr_z_xxyz_xz = pbuffer.data(idx_gr_gd + 206);

    auto grr_z_xxyz_yy = pbuffer.data(idx_gr_gd + 207);

    auto grr_z_xxyz_yz = pbuffer.data(idx_gr_gd + 208);

    auto grr_z_xxyz_zz = pbuffer.data(idx_gr_gd + 209);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_xx, gr_xxy_xy, gr_xxy_xz, gr_xxy_yy, gr_xxy_yz, gr_xxy_zz, gr_xxyz_x, gr_xxyz_xx, gr_xxyz_xy, gr_xxyz_xz, gr_xxyz_y, gr_xxyz_yy, gr_xxyz_yz, gr_xxyz_z, gr_xxyz_zz, grr_z_xxyz_xx, grr_z_xxyz_xy, grr_z_xxyz_xz, grr_z_xxyz_yy, grr_z_xxyz_yz, grr_z_xxyz_zz, ts_xxy_xx, ts_xxy_xy, ts_xxy_xz, ts_xxy_yy, ts_xxy_yz, ts_xxy_zz, ts_xxyz_x, ts_xxyz_xx, ts_xxyz_xy, ts_xxyz_xz, ts_xxyz_y, ts_xxyz_yy, ts_xxyz_yz, ts_xxyz_z, ts_xxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxyz_xx[i] = ts_xxy_xx[i] * gfe2_0 + gr_xxy_xx[i] * gfe_0 + ts_xxyz_xx[i] * gfe_0 * gc_z[i] + gr_xxyz_xx[i] * gc_z[i];

        grr_z_xxyz_xy[i] = ts_xxy_xy[i] * gfe2_0 + gr_xxy_xy[i] * gfe_0 + ts_xxyz_xy[i] * gfe_0 * gc_z[i] + gr_xxyz_xy[i] * gc_z[i];

        grr_z_xxyz_xz[i] = ts_xxy_xz[i] * gfe2_0 + gr_xxy_xz[i] * gfe_0 + ts_xxyz_x[i] * gfe2_0 + gr_xxyz_x[i] * gfe_0 + ts_xxyz_xz[i] * gfe_0 * gc_z[i] + gr_xxyz_xz[i] * gc_z[i];

        grr_z_xxyz_yy[i] = ts_xxy_yy[i] * gfe2_0 + gr_xxy_yy[i] * gfe_0 + ts_xxyz_yy[i] * gfe_0 * gc_z[i] + gr_xxyz_yy[i] * gc_z[i];

        grr_z_xxyz_yz[i] = ts_xxy_yz[i] * gfe2_0 + gr_xxy_yz[i] * gfe_0 + ts_xxyz_y[i] * gfe2_0 + gr_xxyz_y[i] * gfe_0 + ts_xxyz_yz[i] * gfe_0 * gc_z[i] + gr_xxyz_yz[i] * gc_z[i];

        grr_z_xxyz_zz[i] = ts_xxy_zz[i] * gfe2_0 + gr_xxy_zz[i] * gfe_0 + 2.0 * ts_xxyz_z[i] * gfe2_0 + 2.0 * gr_xxyz_z[i] * gfe_0 + ts_xxyz_zz[i] * gfe_0 * gc_z[i] + gr_xxyz_zz[i] * gc_z[i];
    }

    // Set up 210-216 components of targeted buffer : GD

    auto grr_z_xxzz_xx = pbuffer.data(idx_gr_gd + 210);

    auto grr_z_xxzz_xy = pbuffer.data(idx_gr_gd + 211);

    auto grr_z_xxzz_xz = pbuffer.data(idx_gr_gd + 212);

    auto grr_z_xxzz_yy = pbuffer.data(idx_gr_gd + 213);

    auto grr_z_xxzz_yz = pbuffer.data(idx_gr_gd + 214);

    auto grr_z_xxzz_zz = pbuffer.data(idx_gr_gd + 215);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_xx, gr_xxz_xy, gr_xxz_xz, gr_xxz_yy, gr_xxz_yz, gr_xxz_zz, gr_xxzz_x, gr_xxzz_xx, gr_xxzz_xy, gr_xxzz_xz, gr_xxzz_y, gr_xxzz_yy, gr_xxzz_yz, gr_xxzz_z, gr_xxzz_zz, grr_z_xxzz_xx, grr_z_xxzz_xy, grr_z_xxzz_xz, grr_z_xxzz_yy, grr_z_xxzz_yz, grr_z_xxzz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxz_xz, ts_xxz_yy, ts_xxz_yz, ts_xxz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_y, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxzz_xx[i] = 2.0 * ts_xxz_xx[i] * gfe2_0 + 2.0 * gr_xxz_xx[i] * gfe_0 + ts_xxzz_xx[i] * gfe_0 * gc_z[i] + gr_xxzz_xx[i] * gc_z[i];

        grr_z_xxzz_xy[i] = 2.0 * ts_xxz_xy[i] * gfe2_0 + 2.0 * gr_xxz_xy[i] * gfe_0 + ts_xxzz_xy[i] * gfe_0 * gc_z[i] + gr_xxzz_xy[i] * gc_z[i];

        grr_z_xxzz_xz[i] = 2.0 * ts_xxz_xz[i] * gfe2_0 + 2.0 * gr_xxz_xz[i] * gfe_0 + ts_xxzz_x[i] * gfe2_0 + gr_xxzz_x[i] * gfe_0 + ts_xxzz_xz[i] * gfe_0 * gc_z[i] + gr_xxzz_xz[i] * gc_z[i];

        grr_z_xxzz_yy[i] = 2.0 * ts_xxz_yy[i] * gfe2_0 + 2.0 * gr_xxz_yy[i] * gfe_0 + ts_xxzz_yy[i] * gfe_0 * gc_z[i] + gr_xxzz_yy[i] * gc_z[i];

        grr_z_xxzz_yz[i] = 2.0 * ts_xxz_yz[i] * gfe2_0 + 2.0 * gr_xxz_yz[i] * gfe_0 + ts_xxzz_y[i] * gfe2_0 + gr_xxzz_y[i] * gfe_0 + ts_xxzz_yz[i] * gfe_0 * gc_z[i] + gr_xxzz_yz[i] * gc_z[i];

        grr_z_xxzz_zz[i] = 2.0 * ts_xxz_zz[i] * gfe2_0 + 2.0 * gr_xxz_zz[i] * gfe_0 + 2.0 * ts_xxzz_z[i] * gfe2_0 + 2.0 * gr_xxzz_z[i] * gfe_0 + ts_xxzz_zz[i] * gfe_0 * gc_z[i] + gr_xxzz_zz[i] * gc_z[i];
    }

    // Set up 216-222 components of targeted buffer : GD

    auto grr_z_xyyy_xx = pbuffer.data(idx_gr_gd + 216);

    auto grr_z_xyyy_xy = pbuffer.data(idx_gr_gd + 217);

    auto grr_z_xyyy_xz = pbuffer.data(idx_gr_gd + 218);

    auto grr_z_xyyy_yy = pbuffer.data(idx_gr_gd + 219);

    auto grr_z_xyyy_yz = pbuffer.data(idx_gr_gd + 220);

    auto grr_z_xyyy_zz = pbuffer.data(idx_gr_gd + 221);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_x, gr_xyyy_xx, gr_xyyy_xy, gr_xyyy_xz, gr_xyyy_y, gr_xyyy_yy, gr_xyyy_yz, gr_xyyy_z, gr_xyyy_zz, grr_z_xyyy_xx, grr_z_xyyy_xy, grr_z_xyyy_xz, grr_z_xyyy_yy, grr_z_xyyy_yz, grr_z_xyyy_zz, ts_xyyy_x, ts_xyyy_xx, ts_xyyy_xy, ts_xyyy_xz, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_z, ts_xyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyyy_xx[i] = ts_xyyy_xx[i] * gfe_0 * gc_z[i] + gr_xyyy_xx[i] * gc_z[i];

        grr_z_xyyy_xy[i] = ts_xyyy_xy[i] * gfe_0 * gc_z[i] + gr_xyyy_xy[i] * gc_z[i];

        grr_z_xyyy_xz[i] = ts_xyyy_x[i] * gfe2_0 + gr_xyyy_x[i] * gfe_0 + ts_xyyy_xz[i] * gfe_0 * gc_z[i] + gr_xyyy_xz[i] * gc_z[i];

        grr_z_xyyy_yy[i] = ts_xyyy_yy[i] * gfe_0 * gc_z[i] + gr_xyyy_yy[i] * gc_z[i];

        grr_z_xyyy_yz[i] = ts_xyyy_y[i] * gfe2_0 + gr_xyyy_y[i] * gfe_0 + ts_xyyy_yz[i] * gfe_0 * gc_z[i] + gr_xyyy_yz[i] * gc_z[i];

        grr_z_xyyy_zz[i] = 2.0 * ts_xyyy_z[i] * gfe2_0 + 2.0 * gr_xyyy_z[i] * gfe_0 + ts_xyyy_zz[i] * gfe_0 * gc_z[i] + gr_xyyy_zz[i] * gc_z[i];
    }

    // Set up 222-228 components of targeted buffer : GD

    auto grr_z_xyyz_xx = pbuffer.data(idx_gr_gd + 222);

    auto grr_z_xyyz_xy = pbuffer.data(idx_gr_gd + 223);

    auto grr_z_xyyz_xz = pbuffer.data(idx_gr_gd + 224);

    auto grr_z_xyyz_yy = pbuffer.data(idx_gr_gd + 225);

    auto grr_z_xyyz_yz = pbuffer.data(idx_gr_gd + 226);

    auto grr_z_xyyz_zz = pbuffer.data(idx_gr_gd + 227);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_xx, gr_xyy_xy, gr_xyy_xz, gr_xyy_yy, gr_xyy_yz, gr_xyy_zz, gr_xyyz_x, gr_xyyz_xx, gr_xyyz_xy, gr_xyyz_xz, gr_xyyz_y, gr_xyyz_yy, gr_xyyz_yz, gr_xyyz_z, gr_xyyz_zz, grr_z_xyyz_xx, grr_z_xyyz_xy, grr_z_xyyz_xz, grr_z_xyyz_yy, grr_z_xyyz_yz, grr_z_xyyz_zz, ts_xyy_xx, ts_xyy_xy, ts_xyy_xz, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, ts_xyyz_x, ts_xyyz_xx, ts_xyyz_xy, ts_xyyz_xz, ts_xyyz_y, ts_xyyz_yy, ts_xyyz_yz, ts_xyyz_z, ts_xyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyyz_xx[i] = ts_xyy_xx[i] * gfe2_0 + gr_xyy_xx[i] * gfe_0 + ts_xyyz_xx[i] * gfe_0 * gc_z[i] + gr_xyyz_xx[i] * gc_z[i];

        grr_z_xyyz_xy[i] = ts_xyy_xy[i] * gfe2_0 + gr_xyy_xy[i] * gfe_0 + ts_xyyz_xy[i] * gfe_0 * gc_z[i] + gr_xyyz_xy[i] * gc_z[i];

        grr_z_xyyz_xz[i] = ts_xyy_xz[i] * gfe2_0 + gr_xyy_xz[i] * gfe_0 + ts_xyyz_x[i] * gfe2_0 + gr_xyyz_x[i] * gfe_0 + ts_xyyz_xz[i] * gfe_0 * gc_z[i] + gr_xyyz_xz[i] * gc_z[i];

        grr_z_xyyz_yy[i] = ts_xyy_yy[i] * gfe2_0 + gr_xyy_yy[i] * gfe_0 + ts_xyyz_yy[i] * gfe_0 * gc_z[i] + gr_xyyz_yy[i] * gc_z[i];

        grr_z_xyyz_yz[i] = ts_xyy_yz[i] * gfe2_0 + gr_xyy_yz[i] * gfe_0 + ts_xyyz_y[i] * gfe2_0 + gr_xyyz_y[i] * gfe_0 + ts_xyyz_yz[i] * gfe_0 * gc_z[i] + gr_xyyz_yz[i] * gc_z[i];

        grr_z_xyyz_zz[i] = ts_xyy_zz[i] * gfe2_0 + gr_xyy_zz[i] * gfe_0 + 2.0 * ts_xyyz_z[i] * gfe2_0 + 2.0 * gr_xyyz_z[i] * gfe_0 + ts_xyyz_zz[i] * gfe_0 * gc_z[i] + gr_xyyz_zz[i] * gc_z[i];
    }

    // Set up 228-234 components of targeted buffer : GD

    auto grr_z_xyzz_xx = pbuffer.data(idx_gr_gd + 228);

    auto grr_z_xyzz_xy = pbuffer.data(idx_gr_gd + 229);

    auto grr_z_xyzz_xz = pbuffer.data(idx_gr_gd + 230);

    auto grr_z_xyzz_yy = pbuffer.data(idx_gr_gd + 231);

    auto grr_z_xyzz_yz = pbuffer.data(idx_gr_gd + 232);

    auto grr_z_xyzz_zz = pbuffer.data(idx_gr_gd + 233);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_xx, gr_xyz_xy, gr_xyz_xz, gr_xyz_yy, gr_xyz_yz, gr_xyz_zz, gr_xyzz_x, gr_xyzz_xx, gr_xyzz_xy, gr_xyzz_xz, gr_xyzz_y, gr_xyzz_yy, gr_xyzz_yz, gr_xyzz_z, gr_xyzz_zz, grr_z_xyzz_xx, grr_z_xyzz_xy, grr_z_xyzz_xz, grr_z_xyzz_yy, grr_z_xyzz_yz, grr_z_xyzz_zz, ts_xyz_xx, ts_xyz_xy, ts_xyz_xz, ts_xyz_yy, ts_xyz_yz, ts_xyz_zz, ts_xyzz_x, ts_xyzz_xx, ts_xyzz_xy, ts_xyzz_xz, ts_xyzz_y, ts_xyzz_yy, ts_xyzz_yz, ts_xyzz_z, ts_xyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyzz_xx[i] = 2.0 * ts_xyz_xx[i] * gfe2_0 + 2.0 * gr_xyz_xx[i] * gfe_0 + ts_xyzz_xx[i] * gfe_0 * gc_z[i] + gr_xyzz_xx[i] * gc_z[i];

        grr_z_xyzz_xy[i] = 2.0 * ts_xyz_xy[i] * gfe2_0 + 2.0 * gr_xyz_xy[i] * gfe_0 + ts_xyzz_xy[i] * gfe_0 * gc_z[i] + gr_xyzz_xy[i] * gc_z[i];

        grr_z_xyzz_xz[i] = 2.0 * ts_xyz_xz[i] * gfe2_0 + 2.0 * gr_xyz_xz[i] * gfe_0 + ts_xyzz_x[i] * gfe2_0 + gr_xyzz_x[i] * gfe_0 + ts_xyzz_xz[i] * gfe_0 * gc_z[i] + gr_xyzz_xz[i] * gc_z[i];

        grr_z_xyzz_yy[i] = 2.0 * ts_xyz_yy[i] * gfe2_0 + 2.0 * gr_xyz_yy[i] * gfe_0 + ts_xyzz_yy[i] * gfe_0 * gc_z[i] + gr_xyzz_yy[i] * gc_z[i];

        grr_z_xyzz_yz[i] = 2.0 * ts_xyz_yz[i] * gfe2_0 + 2.0 * gr_xyz_yz[i] * gfe_0 + ts_xyzz_y[i] * gfe2_0 + gr_xyzz_y[i] * gfe_0 + ts_xyzz_yz[i] * gfe_0 * gc_z[i] + gr_xyzz_yz[i] * gc_z[i];

        grr_z_xyzz_zz[i] = 2.0 * ts_xyz_zz[i] * gfe2_0 + 2.0 * gr_xyz_zz[i] * gfe_0 + 2.0 * ts_xyzz_z[i] * gfe2_0 + 2.0 * gr_xyzz_z[i] * gfe_0 + ts_xyzz_zz[i] * gfe_0 * gc_z[i] + gr_xyzz_zz[i] * gc_z[i];
    }

    // Set up 234-240 components of targeted buffer : GD

    auto grr_z_xzzz_xx = pbuffer.data(idx_gr_gd + 234);

    auto grr_z_xzzz_xy = pbuffer.data(idx_gr_gd + 235);

    auto grr_z_xzzz_xz = pbuffer.data(idx_gr_gd + 236);

    auto grr_z_xzzz_yy = pbuffer.data(idx_gr_gd + 237);

    auto grr_z_xzzz_yz = pbuffer.data(idx_gr_gd + 238);

    auto grr_z_xzzz_zz = pbuffer.data(idx_gr_gd + 239);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_xx, gr_xzz_xy, gr_xzz_xz, gr_xzz_yy, gr_xzz_yz, gr_xzz_zz, gr_xzzz_x, gr_xzzz_xx, gr_xzzz_xy, gr_xzzz_xz, gr_xzzz_y, gr_xzzz_yy, gr_xzzz_yz, gr_xzzz_z, gr_xzzz_zz, grr_z_xzzz_xx, grr_z_xzzz_xy, grr_z_xzzz_xz, grr_z_xzzz_yy, grr_z_xzzz_yz, grr_z_xzzz_zz, ts_xzz_xx, ts_xzz_xy, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, ts_xzzz_x, ts_xzzz_xx, ts_xzzz_xy, ts_xzzz_xz, ts_xzzz_y, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xzzz_xx[i] = 3.0 * ts_xzz_xx[i] * gfe2_0 + 3.0 * gr_xzz_xx[i] * gfe_0 + ts_xzzz_xx[i] * gfe_0 * gc_z[i] + gr_xzzz_xx[i] * gc_z[i];

        grr_z_xzzz_xy[i] = 3.0 * ts_xzz_xy[i] * gfe2_0 + 3.0 * gr_xzz_xy[i] * gfe_0 + ts_xzzz_xy[i] * gfe_0 * gc_z[i] + gr_xzzz_xy[i] * gc_z[i];

        grr_z_xzzz_xz[i] = 3.0 * ts_xzz_xz[i] * gfe2_0 + 3.0 * gr_xzz_xz[i] * gfe_0 + ts_xzzz_x[i] * gfe2_0 + gr_xzzz_x[i] * gfe_0 + ts_xzzz_xz[i] * gfe_0 * gc_z[i] + gr_xzzz_xz[i] * gc_z[i];

        grr_z_xzzz_yy[i] = 3.0 * ts_xzz_yy[i] * gfe2_0 + 3.0 * gr_xzz_yy[i] * gfe_0 + ts_xzzz_yy[i] * gfe_0 * gc_z[i] + gr_xzzz_yy[i] * gc_z[i];

        grr_z_xzzz_yz[i] = 3.0 * ts_xzz_yz[i] * gfe2_0 + 3.0 * gr_xzz_yz[i] * gfe_0 + ts_xzzz_y[i] * gfe2_0 + gr_xzzz_y[i] * gfe_0 + ts_xzzz_yz[i] * gfe_0 * gc_z[i] + gr_xzzz_yz[i] * gc_z[i];

        grr_z_xzzz_zz[i] = 3.0 * ts_xzz_zz[i] * gfe2_0 + 3.0 * gr_xzz_zz[i] * gfe_0 + 2.0 * ts_xzzz_z[i] * gfe2_0 + 2.0 * gr_xzzz_z[i] * gfe_0 + ts_xzzz_zz[i] * gfe_0 * gc_z[i] + gr_xzzz_zz[i] * gc_z[i];
    }

    // Set up 240-246 components of targeted buffer : GD

    auto grr_z_yyyy_xx = pbuffer.data(idx_gr_gd + 240);

    auto grr_z_yyyy_xy = pbuffer.data(idx_gr_gd + 241);

    auto grr_z_yyyy_xz = pbuffer.data(idx_gr_gd + 242);

    auto grr_z_yyyy_yy = pbuffer.data(idx_gr_gd + 243);

    auto grr_z_yyyy_yz = pbuffer.data(idx_gr_gd + 244);

    auto grr_z_yyyy_zz = pbuffer.data(idx_gr_gd + 245);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_x, gr_yyyy_xx, gr_yyyy_xy, gr_yyyy_xz, gr_yyyy_y, gr_yyyy_yy, gr_yyyy_yz, gr_yyyy_z, gr_yyyy_zz, grr_z_yyyy_xx, grr_z_yyyy_xy, grr_z_yyyy_xz, grr_z_yyyy_yy, grr_z_yyyy_yz, grr_z_yyyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyyy_xx[i] = ts_yyyy_xx[i] * gfe_0 * gc_z[i] + gr_yyyy_xx[i] * gc_z[i];

        grr_z_yyyy_xy[i] = ts_yyyy_xy[i] * gfe_0 * gc_z[i] + gr_yyyy_xy[i] * gc_z[i];

        grr_z_yyyy_xz[i] = ts_yyyy_x[i] * gfe2_0 + gr_yyyy_x[i] * gfe_0 + ts_yyyy_xz[i] * gfe_0 * gc_z[i] + gr_yyyy_xz[i] * gc_z[i];

        grr_z_yyyy_yy[i] = ts_yyyy_yy[i] * gfe_0 * gc_z[i] + gr_yyyy_yy[i] * gc_z[i];

        grr_z_yyyy_yz[i] = ts_yyyy_y[i] * gfe2_0 + gr_yyyy_y[i] * gfe_0 + ts_yyyy_yz[i] * gfe_0 * gc_z[i] + gr_yyyy_yz[i] * gc_z[i];

        grr_z_yyyy_zz[i] = 2.0 * ts_yyyy_z[i] * gfe2_0 + 2.0 * gr_yyyy_z[i] * gfe_0 + ts_yyyy_zz[i] * gfe_0 * gc_z[i] + gr_yyyy_zz[i] * gc_z[i];
    }

    // Set up 246-252 components of targeted buffer : GD

    auto grr_z_yyyz_xx = pbuffer.data(idx_gr_gd + 246);

    auto grr_z_yyyz_xy = pbuffer.data(idx_gr_gd + 247);

    auto grr_z_yyyz_xz = pbuffer.data(idx_gr_gd + 248);

    auto grr_z_yyyz_yy = pbuffer.data(idx_gr_gd + 249);

    auto grr_z_yyyz_yz = pbuffer.data(idx_gr_gd + 250);

    auto grr_z_yyyz_zz = pbuffer.data(idx_gr_gd + 251);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_xx, gr_yyy_xy, gr_yyy_xz, gr_yyy_yy, gr_yyy_yz, gr_yyy_zz, gr_yyyz_x, gr_yyyz_xx, gr_yyyz_xy, gr_yyyz_xz, gr_yyyz_y, gr_yyyz_yy, gr_yyyz_yz, gr_yyyz_z, gr_yyyz_zz, grr_z_yyyz_xx, grr_z_yyyz_xy, grr_z_yyyz_xz, grr_z_yyyz_yy, grr_z_yyyz_yz, grr_z_yyyz_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, ts_yyyz_x, ts_yyyz_xx, ts_yyyz_xy, ts_yyyz_xz, ts_yyyz_y, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyyz_xx[i] = ts_yyy_xx[i] * gfe2_0 + gr_yyy_xx[i] * gfe_0 + ts_yyyz_xx[i] * gfe_0 * gc_z[i] + gr_yyyz_xx[i] * gc_z[i];

        grr_z_yyyz_xy[i] = ts_yyy_xy[i] * gfe2_0 + gr_yyy_xy[i] * gfe_0 + ts_yyyz_xy[i] * gfe_0 * gc_z[i] + gr_yyyz_xy[i] * gc_z[i];

        grr_z_yyyz_xz[i] = ts_yyy_xz[i] * gfe2_0 + gr_yyy_xz[i] * gfe_0 + ts_yyyz_x[i] * gfe2_0 + gr_yyyz_x[i] * gfe_0 + ts_yyyz_xz[i] * gfe_0 * gc_z[i] + gr_yyyz_xz[i] * gc_z[i];

        grr_z_yyyz_yy[i] = ts_yyy_yy[i] * gfe2_0 + gr_yyy_yy[i] * gfe_0 + ts_yyyz_yy[i] * gfe_0 * gc_z[i] + gr_yyyz_yy[i] * gc_z[i];

        grr_z_yyyz_yz[i] = ts_yyy_yz[i] * gfe2_0 + gr_yyy_yz[i] * gfe_0 + ts_yyyz_y[i] * gfe2_0 + gr_yyyz_y[i] * gfe_0 + ts_yyyz_yz[i] * gfe_0 * gc_z[i] + gr_yyyz_yz[i] * gc_z[i];

        grr_z_yyyz_zz[i] = ts_yyy_zz[i] * gfe2_0 + gr_yyy_zz[i] * gfe_0 + 2.0 * ts_yyyz_z[i] * gfe2_0 + 2.0 * gr_yyyz_z[i] * gfe_0 + ts_yyyz_zz[i] * gfe_0 * gc_z[i] + gr_yyyz_zz[i] * gc_z[i];
    }

    // Set up 252-258 components of targeted buffer : GD

    auto grr_z_yyzz_xx = pbuffer.data(idx_gr_gd + 252);

    auto grr_z_yyzz_xy = pbuffer.data(idx_gr_gd + 253);

    auto grr_z_yyzz_xz = pbuffer.data(idx_gr_gd + 254);

    auto grr_z_yyzz_yy = pbuffer.data(idx_gr_gd + 255);

    auto grr_z_yyzz_yz = pbuffer.data(idx_gr_gd + 256);

    auto grr_z_yyzz_zz = pbuffer.data(idx_gr_gd + 257);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_xx, gr_yyz_xy, gr_yyz_xz, gr_yyz_yy, gr_yyz_yz, gr_yyz_zz, gr_yyzz_x, gr_yyzz_xx, gr_yyzz_xy, gr_yyzz_xz, gr_yyzz_y, gr_yyzz_yy, gr_yyzz_yz, gr_yyzz_z, gr_yyzz_zz, grr_z_yyzz_xx, grr_z_yyzz_xy, grr_z_yyzz_xz, grr_z_yyzz_yy, grr_z_yyzz_yz, grr_z_yyzz_zz, ts_yyz_xx, ts_yyz_xy, ts_yyz_xz, ts_yyz_yy, ts_yyz_yz, ts_yyz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyzz_xx[i] = 2.0 * ts_yyz_xx[i] * gfe2_0 + 2.0 * gr_yyz_xx[i] * gfe_0 + ts_yyzz_xx[i] * gfe_0 * gc_z[i] + gr_yyzz_xx[i] * gc_z[i];

        grr_z_yyzz_xy[i] = 2.0 * ts_yyz_xy[i] * gfe2_0 + 2.0 * gr_yyz_xy[i] * gfe_0 + ts_yyzz_xy[i] * gfe_0 * gc_z[i] + gr_yyzz_xy[i] * gc_z[i];

        grr_z_yyzz_xz[i] = 2.0 * ts_yyz_xz[i] * gfe2_0 + 2.0 * gr_yyz_xz[i] * gfe_0 + ts_yyzz_x[i] * gfe2_0 + gr_yyzz_x[i] * gfe_0 + ts_yyzz_xz[i] * gfe_0 * gc_z[i] + gr_yyzz_xz[i] * gc_z[i];

        grr_z_yyzz_yy[i] = 2.0 * ts_yyz_yy[i] * gfe2_0 + 2.0 * gr_yyz_yy[i] * gfe_0 + ts_yyzz_yy[i] * gfe_0 * gc_z[i] + gr_yyzz_yy[i] * gc_z[i];

        grr_z_yyzz_yz[i] = 2.0 * ts_yyz_yz[i] * gfe2_0 + 2.0 * gr_yyz_yz[i] * gfe_0 + ts_yyzz_y[i] * gfe2_0 + gr_yyzz_y[i] * gfe_0 + ts_yyzz_yz[i] * gfe_0 * gc_z[i] + gr_yyzz_yz[i] * gc_z[i];

        grr_z_yyzz_zz[i] = 2.0 * ts_yyz_zz[i] * gfe2_0 + 2.0 * gr_yyz_zz[i] * gfe_0 + 2.0 * ts_yyzz_z[i] * gfe2_0 + 2.0 * gr_yyzz_z[i] * gfe_0 + ts_yyzz_zz[i] * gfe_0 * gc_z[i] + gr_yyzz_zz[i] * gc_z[i];
    }

    // Set up 258-264 components of targeted buffer : GD

    auto grr_z_yzzz_xx = pbuffer.data(idx_gr_gd + 258);

    auto grr_z_yzzz_xy = pbuffer.data(idx_gr_gd + 259);

    auto grr_z_yzzz_xz = pbuffer.data(idx_gr_gd + 260);

    auto grr_z_yzzz_yy = pbuffer.data(idx_gr_gd + 261);

    auto grr_z_yzzz_yz = pbuffer.data(idx_gr_gd + 262);

    auto grr_z_yzzz_zz = pbuffer.data(idx_gr_gd + 263);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_xx, gr_yzz_xy, gr_yzz_xz, gr_yzz_yy, gr_yzz_yz, gr_yzz_zz, gr_yzzz_x, gr_yzzz_xx, gr_yzzz_xy, gr_yzzz_xz, gr_yzzz_y, gr_yzzz_yy, gr_yzzz_yz, gr_yzzz_z, gr_yzzz_zz, grr_z_yzzz_xx, grr_z_yzzz_xy, grr_z_yzzz_xz, grr_z_yzzz_yy, grr_z_yzzz_yz, grr_z_yzzz_zz, ts_yzz_xx, ts_yzz_xy, ts_yzz_xz, ts_yzz_yy, ts_yzz_yz, ts_yzz_zz, ts_yzzz_x, ts_yzzz_xx, ts_yzzz_xy, ts_yzzz_xz, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yzzz_xx[i] = 3.0 * ts_yzz_xx[i] * gfe2_0 + 3.0 * gr_yzz_xx[i] * gfe_0 + ts_yzzz_xx[i] * gfe_0 * gc_z[i] + gr_yzzz_xx[i] * gc_z[i];

        grr_z_yzzz_xy[i] = 3.0 * ts_yzz_xy[i] * gfe2_0 + 3.0 * gr_yzz_xy[i] * gfe_0 + ts_yzzz_xy[i] * gfe_0 * gc_z[i] + gr_yzzz_xy[i] * gc_z[i];

        grr_z_yzzz_xz[i] = 3.0 * ts_yzz_xz[i] * gfe2_0 + 3.0 * gr_yzz_xz[i] * gfe_0 + ts_yzzz_x[i] * gfe2_0 + gr_yzzz_x[i] * gfe_0 + ts_yzzz_xz[i] * gfe_0 * gc_z[i] + gr_yzzz_xz[i] * gc_z[i];

        grr_z_yzzz_yy[i] = 3.0 * ts_yzz_yy[i] * gfe2_0 + 3.0 * gr_yzz_yy[i] * gfe_0 + ts_yzzz_yy[i] * gfe_0 * gc_z[i] + gr_yzzz_yy[i] * gc_z[i];

        grr_z_yzzz_yz[i] = 3.0 * ts_yzz_yz[i] * gfe2_0 + 3.0 * gr_yzz_yz[i] * gfe_0 + ts_yzzz_y[i] * gfe2_0 + gr_yzzz_y[i] * gfe_0 + ts_yzzz_yz[i] * gfe_0 * gc_z[i] + gr_yzzz_yz[i] * gc_z[i];

        grr_z_yzzz_zz[i] = 3.0 * ts_yzz_zz[i] * gfe2_0 + 3.0 * gr_yzz_zz[i] * gfe_0 + 2.0 * ts_yzzz_z[i] * gfe2_0 + 2.0 * gr_yzzz_z[i] * gfe_0 + ts_yzzz_zz[i] * gfe_0 * gc_z[i] + gr_yzzz_zz[i] * gc_z[i];
    }

    // Set up 264-270 components of targeted buffer : GD

    auto grr_z_zzzz_xx = pbuffer.data(idx_gr_gd + 264);

    auto grr_z_zzzz_xy = pbuffer.data(idx_gr_gd + 265);

    auto grr_z_zzzz_xz = pbuffer.data(idx_gr_gd + 266);

    auto grr_z_zzzz_yy = pbuffer.data(idx_gr_gd + 267);

    auto grr_z_zzzz_yz = pbuffer.data(idx_gr_gd + 268);

    auto grr_z_zzzz_zz = pbuffer.data(idx_gr_gd + 269);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_xx, gr_zzz_xy, gr_zzz_xz, gr_zzz_yy, gr_zzz_yz, gr_zzz_zz, gr_zzzz_x, gr_zzzz_xx, gr_zzzz_xy, gr_zzzz_xz, gr_zzzz_y, gr_zzzz_yy, gr_zzzz_yz, gr_zzzz_z, gr_zzzz_zz, grr_z_zzzz_xx, grr_z_zzzz_xy, grr_z_zzzz_xz, grr_z_zzzz_yy, grr_z_zzzz_yz, grr_z_zzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zzzz_xx[i] = 4.0 * ts_zzz_xx[i] * gfe2_0 + 4.0 * gr_zzz_xx[i] * gfe_0 + ts_zzzz_xx[i] * gfe_0 * gc_z[i] + gr_zzzz_xx[i] * gc_z[i];

        grr_z_zzzz_xy[i] = 4.0 * ts_zzz_xy[i] * gfe2_0 + 4.0 * gr_zzz_xy[i] * gfe_0 + ts_zzzz_xy[i] * gfe_0 * gc_z[i] + gr_zzzz_xy[i] * gc_z[i];

        grr_z_zzzz_xz[i] = 4.0 * ts_zzz_xz[i] * gfe2_0 + 4.0 * gr_zzz_xz[i] * gfe_0 + ts_zzzz_x[i] * gfe2_0 + gr_zzzz_x[i] * gfe_0 + ts_zzzz_xz[i] * gfe_0 * gc_z[i] + gr_zzzz_xz[i] * gc_z[i];

        grr_z_zzzz_yy[i] = 4.0 * ts_zzz_yy[i] * gfe2_0 + 4.0 * gr_zzz_yy[i] * gfe_0 + ts_zzzz_yy[i] * gfe_0 * gc_z[i] + gr_zzzz_yy[i] * gc_z[i];

        grr_z_zzzz_yz[i] = 4.0 * ts_zzz_yz[i] * gfe2_0 + 4.0 * gr_zzz_yz[i] * gfe_0 + ts_zzzz_y[i] * gfe2_0 + gr_zzzz_y[i] * gfe_0 + ts_zzzz_yz[i] * gfe_0 * gc_z[i] + gr_zzzz_yz[i] * gc_z[i];

        grr_z_zzzz_zz[i] = 4.0 * ts_zzz_zz[i] * gfe2_0 + 4.0 * gr_zzz_zz[i] * gfe_0 + 2.0 * ts_zzzz_z[i] * gfe2_0 + 2.0 * gr_zzzz_z[i] * gfe_0 + ts_zzzz_zz[i] * gfe_0 * gc_z[i] + gr_zzzz_zz[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

