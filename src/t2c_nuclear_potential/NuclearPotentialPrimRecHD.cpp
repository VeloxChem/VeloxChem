#include "NuclearPotentialPrimRecHD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hd(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hd,
                               const size_t              idx_npot_0_fd,
                               const size_t              idx_npot_1_fd,
                               const size_t              idx_npot_0_gp,
                               const size_t              idx_npot_1_gp,
                               const size_t              idx_npot_0_gd,
                               const size_t              idx_npot_1_gd,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpa,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_0 = pbuffer.data(idx_npot_0_fd);

    auto ta_xxx_xy_0 = pbuffer.data(idx_npot_0_fd + 1);

    auto ta_xxx_xz_0 = pbuffer.data(idx_npot_0_fd + 2);

    auto ta_xxx_yy_0 = pbuffer.data(idx_npot_0_fd + 3);

    auto ta_xxx_yz_0 = pbuffer.data(idx_npot_0_fd + 4);

    auto ta_xxx_zz_0 = pbuffer.data(idx_npot_0_fd + 5);

    auto ta_xxy_xx_0 = pbuffer.data(idx_npot_0_fd + 6);

    auto ta_xxy_xz_0 = pbuffer.data(idx_npot_0_fd + 8);

    auto ta_xxy_yy_0 = pbuffer.data(idx_npot_0_fd + 9);

    auto ta_xxy_yz_0 = pbuffer.data(idx_npot_0_fd + 10);

    auto ta_xxz_xx_0 = pbuffer.data(idx_npot_0_fd + 12);

    auto ta_xxz_xy_0 = pbuffer.data(idx_npot_0_fd + 13);

    auto ta_xxz_xz_0 = pbuffer.data(idx_npot_0_fd + 14);

    auto ta_xxz_yz_0 = pbuffer.data(idx_npot_0_fd + 16);

    auto ta_xxz_zz_0 = pbuffer.data(idx_npot_0_fd + 17);

    auto ta_xyy_xy_0 = pbuffer.data(idx_npot_0_fd + 19);

    auto ta_xyy_yy_0 = pbuffer.data(idx_npot_0_fd + 21);

    auto ta_xyy_yz_0 = pbuffer.data(idx_npot_0_fd + 22);

    auto ta_xyy_zz_0 = pbuffer.data(idx_npot_0_fd + 23);

    auto ta_xyz_yz_0 = pbuffer.data(idx_npot_0_fd + 28);

    auto ta_xzz_xz_0 = pbuffer.data(idx_npot_0_fd + 32);

    auto ta_xzz_yy_0 = pbuffer.data(idx_npot_0_fd + 33);

    auto ta_xzz_yz_0 = pbuffer.data(idx_npot_0_fd + 34);

    auto ta_xzz_zz_0 = pbuffer.data(idx_npot_0_fd + 35);

    auto ta_yyy_xx_0 = pbuffer.data(idx_npot_0_fd + 36);

    auto ta_yyy_xy_0 = pbuffer.data(idx_npot_0_fd + 37);

    auto ta_yyy_xz_0 = pbuffer.data(idx_npot_0_fd + 38);

    auto ta_yyy_yy_0 = pbuffer.data(idx_npot_0_fd + 39);

    auto ta_yyy_yz_0 = pbuffer.data(idx_npot_0_fd + 40);

    auto ta_yyy_zz_0 = pbuffer.data(idx_npot_0_fd + 41);

    auto ta_yyz_xy_0 = pbuffer.data(idx_npot_0_fd + 43);

    auto ta_yyz_xz_0 = pbuffer.data(idx_npot_0_fd + 44);

    auto ta_yyz_yy_0 = pbuffer.data(idx_npot_0_fd + 45);

    auto ta_yyz_yz_0 = pbuffer.data(idx_npot_0_fd + 46);

    auto ta_yyz_zz_0 = pbuffer.data(idx_npot_0_fd + 47);

    auto ta_yzz_xx_0 = pbuffer.data(idx_npot_0_fd + 48);

    auto ta_yzz_xz_0 = pbuffer.data(idx_npot_0_fd + 50);

    auto ta_yzz_yy_0 = pbuffer.data(idx_npot_0_fd + 51);

    auto ta_yzz_yz_0 = pbuffer.data(idx_npot_0_fd + 52);

    auto ta_yzz_zz_0 = pbuffer.data(idx_npot_0_fd + 53);

    auto ta_zzz_xx_0 = pbuffer.data(idx_npot_0_fd + 54);

    auto ta_zzz_xy_0 = pbuffer.data(idx_npot_0_fd + 55);

    auto ta_zzz_xz_0 = pbuffer.data(idx_npot_0_fd + 56);

    auto ta_zzz_yy_0 = pbuffer.data(idx_npot_0_fd + 57);

    auto ta_zzz_yz_0 = pbuffer.data(idx_npot_0_fd + 58);

    auto ta_zzz_zz_0 = pbuffer.data(idx_npot_0_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_1 = pbuffer.data(idx_npot_1_fd);

    auto ta_xxx_xy_1 = pbuffer.data(idx_npot_1_fd + 1);

    auto ta_xxx_xz_1 = pbuffer.data(idx_npot_1_fd + 2);

    auto ta_xxx_yy_1 = pbuffer.data(idx_npot_1_fd + 3);

    auto ta_xxx_yz_1 = pbuffer.data(idx_npot_1_fd + 4);

    auto ta_xxx_zz_1 = pbuffer.data(idx_npot_1_fd + 5);

    auto ta_xxy_xx_1 = pbuffer.data(idx_npot_1_fd + 6);

    auto ta_xxy_xz_1 = pbuffer.data(idx_npot_1_fd + 8);

    auto ta_xxy_yy_1 = pbuffer.data(idx_npot_1_fd + 9);

    auto ta_xxy_yz_1 = pbuffer.data(idx_npot_1_fd + 10);

    auto ta_xxz_xx_1 = pbuffer.data(idx_npot_1_fd + 12);

    auto ta_xxz_xy_1 = pbuffer.data(idx_npot_1_fd + 13);

    auto ta_xxz_xz_1 = pbuffer.data(idx_npot_1_fd + 14);

    auto ta_xxz_yz_1 = pbuffer.data(idx_npot_1_fd + 16);

    auto ta_xxz_zz_1 = pbuffer.data(idx_npot_1_fd + 17);

    auto ta_xyy_xy_1 = pbuffer.data(idx_npot_1_fd + 19);

    auto ta_xyy_yy_1 = pbuffer.data(idx_npot_1_fd + 21);

    auto ta_xyy_yz_1 = pbuffer.data(idx_npot_1_fd + 22);

    auto ta_xyy_zz_1 = pbuffer.data(idx_npot_1_fd + 23);

    auto ta_xyz_yz_1 = pbuffer.data(idx_npot_1_fd + 28);

    auto ta_xzz_xz_1 = pbuffer.data(idx_npot_1_fd + 32);

    auto ta_xzz_yy_1 = pbuffer.data(idx_npot_1_fd + 33);

    auto ta_xzz_yz_1 = pbuffer.data(idx_npot_1_fd + 34);

    auto ta_xzz_zz_1 = pbuffer.data(idx_npot_1_fd + 35);

    auto ta_yyy_xx_1 = pbuffer.data(idx_npot_1_fd + 36);

    auto ta_yyy_xy_1 = pbuffer.data(idx_npot_1_fd + 37);

    auto ta_yyy_xz_1 = pbuffer.data(idx_npot_1_fd + 38);

    auto ta_yyy_yy_1 = pbuffer.data(idx_npot_1_fd + 39);

    auto ta_yyy_yz_1 = pbuffer.data(idx_npot_1_fd + 40);

    auto ta_yyy_zz_1 = pbuffer.data(idx_npot_1_fd + 41);

    auto ta_yyz_xy_1 = pbuffer.data(idx_npot_1_fd + 43);

    auto ta_yyz_xz_1 = pbuffer.data(idx_npot_1_fd + 44);

    auto ta_yyz_yy_1 = pbuffer.data(idx_npot_1_fd + 45);

    auto ta_yyz_yz_1 = pbuffer.data(idx_npot_1_fd + 46);

    auto ta_yyz_zz_1 = pbuffer.data(idx_npot_1_fd + 47);

    auto ta_yzz_xx_1 = pbuffer.data(idx_npot_1_fd + 48);

    auto ta_yzz_xz_1 = pbuffer.data(idx_npot_1_fd + 50);

    auto ta_yzz_yy_1 = pbuffer.data(idx_npot_1_fd + 51);

    auto ta_yzz_yz_1 = pbuffer.data(idx_npot_1_fd + 52);

    auto ta_yzz_zz_1 = pbuffer.data(idx_npot_1_fd + 53);

    auto ta_zzz_xx_1 = pbuffer.data(idx_npot_1_fd + 54);

    auto ta_zzz_xy_1 = pbuffer.data(idx_npot_1_fd + 55);

    auto ta_zzz_xz_1 = pbuffer.data(idx_npot_1_fd + 56);

    auto ta_zzz_yy_1 = pbuffer.data(idx_npot_1_fd + 57);

    auto ta_zzz_yz_1 = pbuffer.data(idx_npot_1_fd + 58);

    auto ta_zzz_zz_1 = pbuffer.data(idx_npot_1_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_0 = pbuffer.data(idx_npot_0_gp);

    auto ta_xxxx_y_0 = pbuffer.data(idx_npot_0_gp + 1);

    auto ta_xxxx_z_0 = pbuffer.data(idx_npot_0_gp + 2);

    auto ta_xxyy_y_0 = pbuffer.data(idx_npot_0_gp + 10);

    auto ta_xxzz_x_0 = pbuffer.data(idx_npot_0_gp + 15);

    auto ta_xxzz_z_0 = pbuffer.data(idx_npot_0_gp + 17);

    auto ta_xyyy_y_0 = pbuffer.data(idx_npot_0_gp + 19);

    auto ta_xzzz_z_0 = pbuffer.data(idx_npot_0_gp + 29);

    auto ta_yyyy_x_0 = pbuffer.data(idx_npot_0_gp + 30);

    auto ta_yyyy_y_0 = pbuffer.data(idx_npot_0_gp + 31);

    auto ta_yyyy_z_0 = pbuffer.data(idx_npot_0_gp + 32);

    auto ta_yyyz_z_0 = pbuffer.data(idx_npot_0_gp + 35);

    auto ta_yyzz_x_0 = pbuffer.data(idx_npot_0_gp + 36);

    auto ta_yyzz_y_0 = pbuffer.data(idx_npot_0_gp + 37);

    auto ta_yyzz_z_0 = pbuffer.data(idx_npot_0_gp + 38);

    auto ta_yzzz_y_0 = pbuffer.data(idx_npot_0_gp + 40);

    auto ta_yzzz_z_0 = pbuffer.data(idx_npot_0_gp + 41);

    auto ta_zzzz_x_0 = pbuffer.data(idx_npot_0_gp + 42);

    auto ta_zzzz_y_0 = pbuffer.data(idx_npot_0_gp + 43);

    auto ta_zzzz_z_0 = pbuffer.data(idx_npot_0_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto ta_xxxx_x_1 = pbuffer.data(idx_npot_1_gp);

    auto ta_xxxx_y_1 = pbuffer.data(idx_npot_1_gp + 1);

    auto ta_xxxx_z_1 = pbuffer.data(idx_npot_1_gp + 2);

    auto ta_xxyy_y_1 = pbuffer.data(idx_npot_1_gp + 10);

    auto ta_xxzz_x_1 = pbuffer.data(idx_npot_1_gp + 15);

    auto ta_xxzz_z_1 = pbuffer.data(idx_npot_1_gp + 17);

    auto ta_xyyy_y_1 = pbuffer.data(idx_npot_1_gp + 19);

    auto ta_xzzz_z_1 = pbuffer.data(idx_npot_1_gp + 29);

    auto ta_yyyy_x_1 = pbuffer.data(idx_npot_1_gp + 30);

    auto ta_yyyy_y_1 = pbuffer.data(idx_npot_1_gp + 31);

    auto ta_yyyy_z_1 = pbuffer.data(idx_npot_1_gp + 32);

    auto ta_yyyz_z_1 = pbuffer.data(idx_npot_1_gp + 35);

    auto ta_yyzz_x_1 = pbuffer.data(idx_npot_1_gp + 36);

    auto ta_yyzz_y_1 = pbuffer.data(idx_npot_1_gp + 37);

    auto ta_yyzz_z_1 = pbuffer.data(idx_npot_1_gp + 38);

    auto ta_yzzz_y_1 = pbuffer.data(idx_npot_1_gp + 40);

    auto ta_yzzz_z_1 = pbuffer.data(idx_npot_1_gp + 41);

    auto ta_zzzz_x_1 = pbuffer.data(idx_npot_1_gp + 42);

    auto ta_zzzz_y_1 = pbuffer.data(idx_npot_1_gp + 43);

    auto ta_zzzz_z_1 = pbuffer.data(idx_npot_1_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto ta_xxxx_xx_0 = pbuffer.data(idx_npot_0_gd);

    auto ta_xxxx_xy_0 = pbuffer.data(idx_npot_0_gd + 1);

    auto ta_xxxx_xz_0 = pbuffer.data(idx_npot_0_gd + 2);

    auto ta_xxxx_yy_0 = pbuffer.data(idx_npot_0_gd + 3);

    auto ta_xxxx_yz_0 = pbuffer.data(idx_npot_0_gd + 4);

    auto ta_xxxx_zz_0 = pbuffer.data(idx_npot_0_gd + 5);

    auto ta_xxxy_xx_0 = pbuffer.data(idx_npot_0_gd + 6);

    auto ta_xxxy_xy_0 = pbuffer.data(idx_npot_0_gd + 7);

    auto ta_xxxy_xz_0 = pbuffer.data(idx_npot_0_gd + 8);

    auto ta_xxxy_yy_0 = pbuffer.data(idx_npot_0_gd + 9);

    auto ta_xxxy_yz_0 = pbuffer.data(idx_npot_0_gd + 10);

    auto ta_xxxz_xx_0 = pbuffer.data(idx_npot_0_gd + 12);

    auto ta_xxxz_xy_0 = pbuffer.data(idx_npot_0_gd + 13);

    auto ta_xxxz_xz_0 = pbuffer.data(idx_npot_0_gd + 14);

    auto ta_xxxz_yz_0 = pbuffer.data(idx_npot_0_gd + 16);

    auto ta_xxxz_zz_0 = pbuffer.data(idx_npot_0_gd + 17);

    auto ta_xxyy_xx_0 = pbuffer.data(idx_npot_0_gd + 18);

    auto ta_xxyy_xy_0 = pbuffer.data(idx_npot_0_gd + 19);

    auto ta_xxyy_xz_0 = pbuffer.data(idx_npot_0_gd + 20);

    auto ta_xxyy_yy_0 = pbuffer.data(idx_npot_0_gd + 21);

    auto ta_xxyy_yz_0 = pbuffer.data(idx_npot_0_gd + 22);

    auto ta_xxyy_zz_0 = pbuffer.data(idx_npot_0_gd + 23);

    auto ta_xxyz_xz_0 = pbuffer.data(idx_npot_0_gd + 26);

    auto ta_xxyz_yz_0 = pbuffer.data(idx_npot_0_gd + 28);

    auto ta_xxzz_xx_0 = pbuffer.data(idx_npot_0_gd + 30);

    auto ta_xxzz_xy_0 = pbuffer.data(idx_npot_0_gd + 31);

    auto ta_xxzz_xz_0 = pbuffer.data(idx_npot_0_gd + 32);

    auto ta_xxzz_yy_0 = pbuffer.data(idx_npot_0_gd + 33);

    auto ta_xxzz_yz_0 = pbuffer.data(idx_npot_0_gd + 34);

    auto ta_xxzz_zz_0 = pbuffer.data(idx_npot_0_gd + 35);

    auto ta_xyyy_xx_0 = pbuffer.data(idx_npot_0_gd + 36);

    auto ta_xyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 37);

    auto ta_xyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 39);

    auto ta_xyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 40);

    auto ta_xyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 41);

    auto ta_xyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 46);

    auto ta_xyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 47);

    auto ta_xyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 51);

    auto ta_xyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 52);

    auto ta_xzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 54);

    auto ta_xzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 56);

    auto ta_xzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 57);

    auto ta_xzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 58);

    auto ta_xzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 59);

    auto ta_yyyy_xx_0 = pbuffer.data(idx_npot_0_gd + 60);

    auto ta_yyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 61);

    auto ta_yyyy_xz_0 = pbuffer.data(idx_npot_0_gd + 62);

    auto ta_yyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 63);

    auto ta_yyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 64);

    auto ta_yyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 65);

    auto ta_yyyz_xy_0 = pbuffer.data(idx_npot_0_gd + 67);

    auto ta_yyyz_xz_0 = pbuffer.data(idx_npot_0_gd + 68);

    auto ta_yyyz_yy_0 = pbuffer.data(idx_npot_0_gd + 69);

    auto ta_yyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 70);

    auto ta_yyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 71);

    auto ta_yyzz_xx_0 = pbuffer.data(idx_npot_0_gd + 72);

    auto ta_yyzz_xy_0 = pbuffer.data(idx_npot_0_gd + 73);

    auto ta_yyzz_xz_0 = pbuffer.data(idx_npot_0_gd + 74);

    auto ta_yyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 75);

    auto ta_yyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 76);

    auto ta_yyzz_zz_0 = pbuffer.data(idx_npot_0_gd + 77);

    auto ta_yzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 78);

    auto ta_yzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 79);

    auto ta_yzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 80);

    auto ta_yzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 81);

    auto ta_yzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 82);

    auto ta_yzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 83);

    auto ta_zzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 84);

    auto ta_zzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 85);

    auto ta_zzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 86);

    auto ta_zzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 87);

    auto ta_zzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 88);

    auto ta_zzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto ta_xxxx_xx_1 = pbuffer.data(idx_npot_1_gd);

    auto ta_xxxx_xy_1 = pbuffer.data(idx_npot_1_gd + 1);

    auto ta_xxxx_xz_1 = pbuffer.data(idx_npot_1_gd + 2);

    auto ta_xxxx_yy_1 = pbuffer.data(idx_npot_1_gd + 3);

    auto ta_xxxx_yz_1 = pbuffer.data(idx_npot_1_gd + 4);

    auto ta_xxxx_zz_1 = pbuffer.data(idx_npot_1_gd + 5);

    auto ta_xxxy_xx_1 = pbuffer.data(idx_npot_1_gd + 6);

    auto ta_xxxy_xy_1 = pbuffer.data(idx_npot_1_gd + 7);

    auto ta_xxxy_xz_1 = pbuffer.data(idx_npot_1_gd + 8);

    auto ta_xxxy_yy_1 = pbuffer.data(idx_npot_1_gd + 9);

    auto ta_xxxy_yz_1 = pbuffer.data(idx_npot_1_gd + 10);

    auto ta_xxxz_xx_1 = pbuffer.data(idx_npot_1_gd + 12);

    auto ta_xxxz_xy_1 = pbuffer.data(idx_npot_1_gd + 13);

    auto ta_xxxz_xz_1 = pbuffer.data(idx_npot_1_gd + 14);

    auto ta_xxxz_yz_1 = pbuffer.data(idx_npot_1_gd + 16);

    auto ta_xxxz_zz_1 = pbuffer.data(idx_npot_1_gd + 17);

    auto ta_xxyy_xx_1 = pbuffer.data(idx_npot_1_gd + 18);

    auto ta_xxyy_xy_1 = pbuffer.data(idx_npot_1_gd + 19);

    auto ta_xxyy_xz_1 = pbuffer.data(idx_npot_1_gd + 20);

    auto ta_xxyy_yy_1 = pbuffer.data(idx_npot_1_gd + 21);

    auto ta_xxyy_yz_1 = pbuffer.data(idx_npot_1_gd + 22);

    auto ta_xxyy_zz_1 = pbuffer.data(idx_npot_1_gd + 23);

    auto ta_xxyz_xz_1 = pbuffer.data(idx_npot_1_gd + 26);

    auto ta_xxyz_yz_1 = pbuffer.data(idx_npot_1_gd + 28);

    auto ta_xxzz_xx_1 = pbuffer.data(idx_npot_1_gd + 30);

    auto ta_xxzz_xy_1 = pbuffer.data(idx_npot_1_gd + 31);

    auto ta_xxzz_xz_1 = pbuffer.data(idx_npot_1_gd + 32);

    auto ta_xxzz_yy_1 = pbuffer.data(idx_npot_1_gd + 33);

    auto ta_xxzz_yz_1 = pbuffer.data(idx_npot_1_gd + 34);

    auto ta_xxzz_zz_1 = pbuffer.data(idx_npot_1_gd + 35);

    auto ta_xyyy_xx_1 = pbuffer.data(idx_npot_1_gd + 36);

    auto ta_xyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 37);

    auto ta_xyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 39);

    auto ta_xyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 40);

    auto ta_xyyy_zz_1 = pbuffer.data(idx_npot_1_gd + 41);

    auto ta_xyyz_yz_1 = pbuffer.data(idx_npot_1_gd + 46);

    auto ta_xyyz_zz_1 = pbuffer.data(idx_npot_1_gd + 47);

    auto ta_xyzz_yy_1 = pbuffer.data(idx_npot_1_gd + 51);

    auto ta_xyzz_yz_1 = pbuffer.data(idx_npot_1_gd + 52);

    auto ta_xzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 54);

    auto ta_xzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 56);

    auto ta_xzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 57);

    auto ta_xzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 58);

    auto ta_xzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 59);

    auto ta_yyyy_xx_1 = pbuffer.data(idx_npot_1_gd + 60);

    auto ta_yyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 61);

    auto ta_yyyy_xz_1 = pbuffer.data(idx_npot_1_gd + 62);

    auto ta_yyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 63);

    auto ta_yyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 64);

    auto ta_yyyy_zz_1 = pbuffer.data(idx_npot_1_gd + 65);

    auto ta_yyyz_xy_1 = pbuffer.data(idx_npot_1_gd + 67);

    auto ta_yyyz_xz_1 = pbuffer.data(idx_npot_1_gd + 68);

    auto ta_yyyz_yy_1 = pbuffer.data(idx_npot_1_gd + 69);

    auto ta_yyyz_yz_1 = pbuffer.data(idx_npot_1_gd + 70);

    auto ta_yyyz_zz_1 = pbuffer.data(idx_npot_1_gd + 71);

    auto ta_yyzz_xx_1 = pbuffer.data(idx_npot_1_gd + 72);

    auto ta_yyzz_xy_1 = pbuffer.data(idx_npot_1_gd + 73);

    auto ta_yyzz_xz_1 = pbuffer.data(idx_npot_1_gd + 74);

    auto ta_yyzz_yy_1 = pbuffer.data(idx_npot_1_gd + 75);

    auto ta_yyzz_yz_1 = pbuffer.data(idx_npot_1_gd + 76);

    auto ta_yyzz_zz_1 = pbuffer.data(idx_npot_1_gd + 77);

    auto ta_yzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 78);

    auto ta_yzzz_xy_1 = pbuffer.data(idx_npot_1_gd + 79);

    auto ta_yzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 80);

    auto ta_yzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 81);

    auto ta_yzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 82);

    auto ta_yzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 83);

    auto ta_zzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 84);

    auto ta_zzzz_xy_1 = pbuffer.data(idx_npot_1_gd + 85);

    auto ta_zzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 86);

    auto ta_zzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 87);

    auto ta_zzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 88);

    auto ta_zzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 89);

    // Set up 0-6 components of targeted buffer : HD

    auto ta_xxxxx_xx_0 = pbuffer.data(idx_npot_0_hd);

    auto ta_xxxxx_xy_0 = pbuffer.data(idx_npot_0_hd + 1);

    auto ta_xxxxx_xz_0 = pbuffer.data(idx_npot_0_hd + 2);

    auto ta_xxxxx_yy_0 = pbuffer.data(idx_npot_0_hd + 3);

    auto ta_xxxxx_yz_0 = pbuffer.data(idx_npot_0_hd + 4);

    auto ta_xxxxx_zz_0 = pbuffer.data(idx_npot_0_hd + 5);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xxx_xx_0,   \
                             ta_xxx_xx_1,   \
                             ta_xxx_xy_0,   \
                             ta_xxx_xy_1,   \
                             ta_xxx_xz_0,   \
                             ta_xxx_xz_1,   \
                             ta_xxx_yy_0,   \
                             ta_xxx_yy_1,   \
                             ta_xxx_yz_0,   \
                             ta_xxx_yz_1,   \
                             ta_xxx_zz_0,   \
                             ta_xxx_zz_1,   \
                             ta_xxxx_x_0,   \
                             ta_xxxx_x_1,   \
                             ta_xxxx_xx_0,  \
                             ta_xxxx_xx_1,  \
                             ta_xxxx_xy_0,  \
                             ta_xxxx_xy_1,  \
                             ta_xxxx_xz_0,  \
                             ta_xxxx_xz_1,  \
                             ta_xxxx_y_0,   \
                             ta_xxxx_y_1,   \
                             ta_xxxx_yy_0,  \
                             ta_xxxx_yy_1,  \
                             ta_xxxx_yz_0,  \
                             ta_xxxx_yz_1,  \
                             ta_xxxx_z_0,   \
                             ta_xxxx_z_1,   \
                             ta_xxxx_zz_0,  \
                             ta_xxxx_zz_1,  \
                             ta_xxxxx_xx_0, \
                             ta_xxxxx_xy_0, \
                             ta_xxxxx_xz_0, \
                             ta_xxxxx_yy_0, \
                             ta_xxxxx_yz_0, \
                             ta_xxxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_xx_0[i] = 4.0 * ta_xxx_xx_0[i] * fe_0 - 4.0 * ta_xxx_xx_1[i] * fe_0 + 2.0 * ta_xxxx_x_0[i] * fe_0 - 2.0 * ta_xxxx_x_1[i] * fe_0 +
                           ta_xxxx_xx_0[i] * pa_x[i] - ta_xxxx_xx_1[i] * pc_x[i];

        ta_xxxxx_xy_0[i] = 4.0 * ta_xxx_xy_0[i] * fe_0 - 4.0 * ta_xxx_xy_1[i] * fe_0 + ta_xxxx_y_0[i] * fe_0 - ta_xxxx_y_1[i] * fe_0 +
                           ta_xxxx_xy_0[i] * pa_x[i] - ta_xxxx_xy_1[i] * pc_x[i];

        ta_xxxxx_xz_0[i] = 4.0 * ta_xxx_xz_0[i] * fe_0 - 4.0 * ta_xxx_xz_1[i] * fe_0 + ta_xxxx_z_0[i] * fe_0 - ta_xxxx_z_1[i] * fe_0 +
                           ta_xxxx_xz_0[i] * pa_x[i] - ta_xxxx_xz_1[i] * pc_x[i];

        ta_xxxxx_yy_0[i] = 4.0 * ta_xxx_yy_0[i] * fe_0 - 4.0 * ta_xxx_yy_1[i] * fe_0 + ta_xxxx_yy_0[i] * pa_x[i] - ta_xxxx_yy_1[i] * pc_x[i];

        ta_xxxxx_yz_0[i] = 4.0 * ta_xxx_yz_0[i] * fe_0 - 4.0 * ta_xxx_yz_1[i] * fe_0 + ta_xxxx_yz_0[i] * pa_x[i] - ta_xxxx_yz_1[i] * pc_x[i];

        ta_xxxxx_zz_0[i] = 4.0 * ta_xxx_zz_0[i] * fe_0 - 4.0 * ta_xxx_zz_1[i] * fe_0 + ta_xxxx_zz_0[i] * pa_x[i] - ta_xxxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : HD

    auto ta_xxxxy_xx_0 = pbuffer.data(idx_npot_0_hd + 6);

    auto ta_xxxxy_xy_0 = pbuffer.data(idx_npot_0_hd + 7);

    auto ta_xxxxy_xz_0 = pbuffer.data(idx_npot_0_hd + 8);

    auto ta_xxxxy_yy_0 = pbuffer.data(idx_npot_0_hd + 9);

    auto ta_xxxxy_yz_0 = pbuffer.data(idx_npot_0_hd + 10);

    auto ta_xxxxy_zz_0 = pbuffer.data(idx_npot_0_hd + 11);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxxx_x_0,   \
                             ta_xxxx_x_1,   \
                             ta_xxxx_xx_0,  \
                             ta_xxxx_xx_1,  \
                             ta_xxxx_xy_0,  \
                             ta_xxxx_xy_1,  \
                             ta_xxxx_xz_0,  \
                             ta_xxxx_xz_1,  \
                             ta_xxxx_zz_0,  \
                             ta_xxxx_zz_1,  \
                             ta_xxxxy_xx_0, \
                             ta_xxxxy_xy_0, \
                             ta_xxxxy_xz_0, \
                             ta_xxxxy_yy_0, \
                             ta_xxxxy_yz_0, \
                             ta_xxxxy_zz_0, \
                             ta_xxxy_yy_0,  \
                             ta_xxxy_yy_1,  \
                             ta_xxxy_yz_0,  \
                             ta_xxxy_yz_1,  \
                             ta_xxy_yy_0,   \
                             ta_xxy_yy_1,   \
                             ta_xxy_yz_0,   \
                             ta_xxy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxy_xx_0[i] = ta_xxxx_xx_0[i] * pa_y[i] - ta_xxxx_xx_1[i] * pc_y[i];

        ta_xxxxy_xy_0[i] = ta_xxxx_x_0[i] * fe_0 - ta_xxxx_x_1[i] * fe_0 + ta_xxxx_xy_0[i] * pa_y[i] - ta_xxxx_xy_1[i] * pc_y[i];

        ta_xxxxy_xz_0[i] = ta_xxxx_xz_0[i] * pa_y[i] - ta_xxxx_xz_1[i] * pc_y[i];

        ta_xxxxy_yy_0[i] = 3.0 * ta_xxy_yy_0[i] * fe_0 - 3.0 * ta_xxy_yy_1[i] * fe_0 + ta_xxxy_yy_0[i] * pa_x[i] - ta_xxxy_yy_1[i] * pc_x[i];

        ta_xxxxy_yz_0[i] = 3.0 * ta_xxy_yz_0[i] * fe_0 - 3.0 * ta_xxy_yz_1[i] * fe_0 + ta_xxxy_yz_0[i] * pa_x[i] - ta_xxxy_yz_1[i] * pc_x[i];

        ta_xxxxy_zz_0[i] = ta_xxxx_zz_0[i] * pa_y[i] - ta_xxxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : HD

    auto ta_xxxxz_xx_0 = pbuffer.data(idx_npot_0_hd + 12);

    auto ta_xxxxz_xy_0 = pbuffer.data(idx_npot_0_hd + 13);

    auto ta_xxxxz_xz_0 = pbuffer.data(idx_npot_0_hd + 14);

    auto ta_xxxxz_yy_0 = pbuffer.data(idx_npot_0_hd + 15);

    auto ta_xxxxz_yz_0 = pbuffer.data(idx_npot_0_hd + 16);

    auto ta_xxxxz_zz_0 = pbuffer.data(idx_npot_0_hd + 17);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxxx_x_0,   \
                             ta_xxxx_x_1,   \
                             ta_xxxx_xx_0,  \
                             ta_xxxx_xx_1,  \
                             ta_xxxx_xy_0,  \
                             ta_xxxx_xy_1,  \
                             ta_xxxx_xz_0,  \
                             ta_xxxx_xz_1,  \
                             ta_xxxx_yy_0,  \
                             ta_xxxx_yy_1,  \
                             ta_xxxxz_xx_0, \
                             ta_xxxxz_xy_0, \
                             ta_xxxxz_xz_0, \
                             ta_xxxxz_yy_0, \
                             ta_xxxxz_yz_0, \
                             ta_xxxxz_zz_0, \
                             ta_xxxz_yz_0,  \
                             ta_xxxz_yz_1,  \
                             ta_xxxz_zz_0,  \
                             ta_xxxz_zz_1,  \
                             ta_xxz_yz_0,   \
                             ta_xxz_yz_1,   \
                             ta_xxz_zz_0,   \
                             ta_xxz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxz_xx_0[i] = ta_xxxx_xx_0[i] * pa_z[i] - ta_xxxx_xx_1[i] * pc_z[i];

        ta_xxxxz_xy_0[i] = ta_xxxx_xy_0[i] * pa_z[i] - ta_xxxx_xy_1[i] * pc_z[i];

        ta_xxxxz_xz_0[i] = ta_xxxx_x_0[i] * fe_0 - ta_xxxx_x_1[i] * fe_0 + ta_xxxx_xz_0[i] * pa_z[i] - ta_xxxx_xz_1[i] * pc_z[i];

        ta_xxxxz_yy_0[i] = ta_xxxx_yy_0[i] * pa_z[i] - ta_xxxx_yy_1[i] * pc_z[i];

        ta_xxxxz_yz_0[i] = 3.0 * ta_xxz_yz_0[i] * fe_0 - 3.0 * ta_xxz_yz_1[i] * fe_0 + ta_xxxz_yz_0[i] * pa_x[i] - ta_xxxz_yz_1[i] * pc_x[i];

        ta_xxxxz_zz_0[i] = 3.0 * ta_xxz_zz_0[i] * fe_0 - 3.0 * ta_xxz_zz_1[i] * fe_0 + ta_xxxz_zz_0[i] * pa_x[i] - ta_xxxz_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : HD

    auto ta_xxxyy_xx_0 = pbuffer.data(idx_npot_0_hd + 18);

    auto ta_xxxyy_xy_0 = pbuffer.data(idx_npot_0_hd + 19);

    auto ta_xxxyy_xz_0 = pbuffer.data(idx_npot_0_hd + 20);

    auto ta_xxxyy_yy_0 = pbuffer.data(idx_npot_0_hd + 21);

    auto ta_xxxyy_yz_0 = pbuffer.data(idx_npot_0_hd + 22);

    auto ta_xxxyy_zz_0 = pbuffer.data(idx_npot_0_hd + 23);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxx_xx_0,   \
                             ta_xxx_xx_1,   \
                             ta_xxx_xz_0,   \
                             ta_xxx_xz_1,   \
                             ta_xxxy_xx_0,  \
                             ta_xxxy_xx_1,  \
                             ta_xxxy_xz_0,  \
                             ta_xxxy_xz_1,  \
                             ta_xxxyy_xx_0, \
                             ta_xxxyy_xy_0, \
                             ta_xxxyy_xz_0, \
                             ta_xxxyy_yy_0, \
                             ta_xxxyy_yz_0, \
                             ta_xxxyy_zz_0, \
                             ta_xxyy_xy_0,  \
                             ta_xxyy_xy_1,  \
                             ta_xxyy_y_0,   \
                             ta_xxyy_y_1,   \
                             ta_xxyy_yy_0,  \
                             ta_xxyy_yy_1,  \
                             ta_xxyy_yz_0,  \
                             ta_xxyy_yz_1,  \
                             ta_xxyy_zz_0,  \
                             ta_xxyy_zz_1,  \
                             ta_xyy_xy_0,   \
                             ta_xyy_xy_1,   \
                             ta_xyy_yy_0,   \
                             ta_xyy_yy_1,   \
                             ta_xyy_yz_0,   \
                             ta_xyy_yz_1,   \
                             ta_xyy_zz_0,   \
                             ta_xyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyy_xx_0[i] = ta_xxx_xx_0[i] * fe_0 - ta_xxx_xx_1[i] * fe_0 + ta_xxxy_xx_0[i] * pa_y[i] - ta_xxxy_xx_1[i] * pc_y[i];

        ta_xxxyy_xy_0[i] = 2.0 * ta_xyy_xy_0[i] * fe_0 - 2.0 * ta_xyy_xy_1[i] * fe_0 + ta_xxyy_y_0[i] * fe_0 - ta_xxyy_y_1[i] * fe_0 +
                           ta_xxyy_xy_0[i] * pa_x[i] - ta_xxyy_xy_1[i] * pc_x[i];

        ta_xxxyy_xz_0[i] = ta_xxx_xz_0[i] * fe_0 - ta_xxx_xz_1[i] * fe_0 + ta_xxxy_xz_0[i] * pa_y[i] - ta_xxxy_xz_1[i] * pc_y[i];

        ta_xxxyy_yy_0[i] = 2.0 * ta_xyy_yy_0[i] * fe_0 - 2.0 * ta_xyy_yy_1[i] * fe_0 + ta_xxyy_yy_0[i] * pa_x[i] - ta_xxyy_yy_1[i] * pc_x[i];

        ta_xxxyy_yz_0[i] = 2.0 * ta_xyy_yz_0[i] * fe_0 - 2.0 * ta_xyy_yz_1[i] * fe_0 + ta_xxyy_yz_0[i] * pa_x[i] - ta_xxyy_yz_1[i] * pc_x[i];

        ta_xxxyy_zz_0[i] = 2.0 * ta_xyy_zz_0[i] * fe_0 - 2.0 * ta_xyy_zz_1[i] * fe_0 + ta_xxyy_zz_0[i] * pa_x[i] - ta_xxyy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : HD

    auto ta_xxxyz_xx_0 = pbuffer.data(idx_npot_0_hd + 24);

    auto ta_xxxyz_xy_0 = pbuffer.data(idx_npot_0_hd + 25);

    auto ta_xxxyz_xz_0 = pbuffer.data(idx_npot_0_hd + 26);

    auto ta_xxxyz_yy_0 = pbuffer.data(idx_npot_0_hd + 27);

    auto ta_xxxyz_yz_0 = pbuffer.data(idx_npot_0_hd + 28);

    auto ta_xxxyz_zz_0 = pbuffer.data(idx_npot_0_hd + 29);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta_xxxy_xy_0,  \
                             ta_xxxy_xy_1,  \
                             ta_xxxy_yy_0,  \
                             ta_xxxy_yy_1,  \
                             ta_xxxyz_xx_0, \
                             ta_xxxyz_xy_0, \
                             ta_xxxyz_xz_0, \
                             ta_xxxyz_yy_0, \
                             ta_xxxyz_yz_0, \
                             ta_xxxyz_zz_0, \
                             ta_xxxz_xx_0,  \
                             ta_xxxz_xx_1,  \
                             ta_xxxz_xz_0,  \
                             ta_xxxz_xz_1,  \
                             ta_xxxz_zz_0,  \
                             ta_xxxz_zz_1,  \
                             ta_xxyz_yz_0,  \
                             ta_xxyz_yz_1,  \
                             ta_xyz_yz_0,   \
                             ta_xyz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyz_xx_0[i] = ta_xxxz_xx_0[i] * pa_y[i] - ta_xxxz_xx_1[i] * pc_y[i];

        ta_xxxyz_xy_0[i] = ta_xxxy_xy_0[i] * pa_z[i] - ta_xxxy_xy_1[i] * pc_z[i];

        ta_xxxyz_xz_0[i] = ta_xxxz_xz_0[i] * pa_y[i] - ta_xxxz_xz_1[i] * pc_y[i];

        ta_xxxyz_yy_0[i] = ta_xxxy_yy_0[i] * pa_z[i] - ta_xxxy_yy_1[i] * pc_z[i];

        ta_xxxyz_yz_0[i] = 2.0 * ta_xyz_yz_0[i] * fe_0 - 2.0 * ta_xyz_yz_1[i] * fe_0 + ta_xxyz_yz_0[i] * pa_x[i] - ta_xxyz_yz_1[i] * pc_x[i];

        ta_xxxyz_zz_0[i] = ta_xxxz_zz_0[i] * pa_y[i] - ta_xxxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : HD

    auto ta_xxxzz_xx_0 = pbuffer.data(idx_npot_0_hd + 30);

    auto ta_xxxzz_xy_0 = pbuffer.data(idx_npot_0_hd + 31);

    auto ta_xxxzz_xz_0 = pbuffer.data(idx_npot_0_hd + 32);

    auto ta_xxxzz_yy_0 = pbuffer.data(idx_npot_0_hd + 33);

    auto ta_xxxzz_yz_0 = pbuffer.data(idx_npot_0_hd + 34);

    auto ta_xxxzz_zz_0 = pbuffer.data(idx_npot_0_hd + 35);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxx_xx_0,   \
                             ta_xxx_xx_1,   \
                             ta_xxx_xy_0,   \
                             ta_xxx_xy_1,   \
                             ta_xxxz_xx_0,  \
                             ta_xxxz_xx_1,  \
                             ta_xxxz_xy_0,  \
                             ta_xxxz_xy_1,  \
                             ta_xxxzz_xx_0, \
                             ta_xxxzz_xy_0, \
                             ta_xxxzz_xz_0, \
                             ta_xxxzz_yy_0, \
                             ta_xxxzz_yz_0, \
                             ta_xxxzz_zz_0, \
                             ta_xxzz_xz_0,  \
                             ta_xxzz_xz_1,  \
                             ta_xxzz_yy_0,  \
                             ta_xxzz_yy_1,  \
                             ta_xxzz_yz_0,  \
                             ta_xxzz_yz_1,  \
                             ta_xxzz_z_0,   \
                             ta_xxzz_z_1,   \
                             ta_xxzz_zz_0,  \
                             ta_xxzz_zz_1,  \
                             ta_xzz_xz_0,   \
                             ta_xzz_xz_1,   \
                             ta_xzz_yy_0,   \
                             ta_xzz_yy_1,   \
                             ta_xzz_yz_0,   \
                             ta_xzz_yz_1,   \
                             ta_xzz_zz_0,   \
                             ta_xzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzz_xx_0[i] = ta_xxx_xx_0[i] * fe_0 - ta_xxx_xx_1[i] * fe_0 + ta_xxxz_xx_0[i] * pa_z[i] - ta_xxxz_xx_1[i] * pc_z[i];

        ta_xxxzz_xy_0[i] = ta_xxx_xy_0[i] * fe_0 - ta_xxx_xy_1[i] * fe_0 + ta_xxxz_xy_0[i] * pa_z[i] - ta_xxxz_xy_1[i] * pc_z[i];

        ta_xxxzz_xz_0[i] = 2.0 * ta_xzz_xz_0[i] * fe_0 - 2.0 * ta_xzz_xz_1[i] * fe_0 + ta_xxzz_z_0[i] * fe_0 - ta_xxzz_z_1[i] * fe_0 +
                           ta_xxzz_xz_0[i] * pa_x[i] - ta_xxzz_xz_1[i] * pc_x[i];

        ta_xxxzz_yy_0[i] = 2.0 * ta_xzz_yy_0[i] * fe_0 - 2.0 * ta_xzz_yy_1[i] * fe_0 + ta_xxzz_yy_0[i] * pa_x[i] - ta_xxzz_yy_1[i] * pc_x[i];

        ta_xxxzz_yz_0[i] = 2.0 * ta_xzz_yz_0[i] * fe_0 - 2.0 * ta_xzz_yz_1[i] * fe_0 + ta_xxzz_yz_0[i] * pa_x[i] - ta_xxzz_yz_1[i] * pc_x[i];

        ta_xxxzz_zz_0[i] = 2.0 * ta_xzz_zz_0[i] * fe_0 - 2.0 * ta_xzz_zz_1[i] * fe_0 + ta_xxzz_zz_0[i] * pa_x[i] - ta_xxzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : HD

    auto ta_xxyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 36);

    auto ta_xxyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 37);

    auto ta_xxyyy_xz_0 = pbuffer.data(idx_npot_0_hd + 38);

    auto ta_xxyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 39);

    auto ta_xxyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 40);

    auto ta_xxyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 41);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxy_xx_0,   \
                             ta_xxy_xx_1,   \
                             ta_xxy_xz_0,   \
                             ta_xxy_xz_1,   \
                             ta_xxyy_xx_0,  \
                             ta_xxyy_xx_1,  \
                             ta_xxyy_xz_0,  \
                             ta_xxyy_xz_1,  \
                             ta_xxyyy_xx_0, \
                             ta_xxyyy_xy_0, \
                             ta_xxyyy_xz_0, \
                             ta_xxyyy_yy_0, \
                             ta_xxyyy_yz_0, \
                             ta_xxyyy_zz_0, \
                             ta_xyyy_xy_0,  \
                             ta_xyyy_xy_1,  \
                             ta_xyyy_y_0,   \
                             ta_xyyy_y_1,   \
                             ta_xyyy_yy_0,  \
                             ta_xyyy_yy_1,  \
                             ta_xyyy_yz_0,  \
                             ta_xyyy_yz_1,  \
                             ta_xyyy_zz_0,  \
                             ta_xyyy_zz_1,  \
                             ta_yyy_xy_0,   \
                             ta_yyy_xy_1,   \
                             ta_yyy_yy_0,   \
                             ta_yyy_yy_1,   \
                             ta_yyy_yz_0,   \
                             ta_yyy_yz_1,   \
                             ta_yyy_zz_0,   \
                             ta_yyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyy_xx_0[i] = 2.0 * ta_xxy_xx_0[i] * fe_0 - 2.0 * ta_xxy_xx_1[i] * fe_0 + ta_xxyy_xx_0[i] * pa_y[i] - ta_xxyy_xx_1[i] * pc_y[i];

        ta_xxyyy_xy_0[i] = ta_yyy_xy_0[i] * fe_0 - ta_yyy_xy_1[i] * fe_0 + ta_xyyy_y_0[i] * fe_0 - ta_xyyy_y_1[i] * fe_0 + ta_xyyy_xy_0[i] * pa_x[i] -
                           ta_xyyy_xy_1[i] * pc_x[i];

        ta_xxyyy_xz_0[i] = 2.0 * ta_xxy_xz_0[i] * fe_0 - 2.0 * ta_xxy_xz_1[i] * fe_0 + ta_xxyy_xz_0[i] * pa_y[i] - ta_xxyy_xz_1[i] * pc_y[i];

        ta_xxyyy_yy_0[i] = ta_yyy_yy_0[i] * fe_0 - ta_yyy_yy_1[i] * fe_0 + ta_xyyy_yy_0[i] * pa_x[i] - ta_xyyy_yy_1[i] * pc_x[i];

        ta_xxyyy_yz_0[i] = ta_yyy_yz_0[i] * fe_0 - ta_yyy_yz_1[i] * fe_0 + ta_xyyy_yz_0[i] * pa_x[i] - ta_xyyy_yz_1[i] * pc_x[i];

        ta_xxyyy_zz_0[i] = ta_yyy_zz_0[i] * fe_0 - ta_yyy_zz_1[i] * fe_0 + ta_xyyy_zz_0[i] * pa_x[i] - ta_xyyy_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : HD

    auto ta_xxyyz_xx_0 = pbuffer.data(idx_npot_0_hd + 42);

    auto ta_xxyyz_xy_0 = pbuffer.data(idx_npot_0_hd + 43);

    auto ta_xxyyz_xz_0 = pbuffer.data(idx_npot_0_hd + 44);

    auto ta_xxyyz_yy_0 = pbuffer.data(idx_npot_0_hd + 45);

    auto ta_xxyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 46);

    auto ta_xxyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 47);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta_xxyy_xx_0,  \
                             ta_xxyy_xx_1,  \
                             ta_xxyy_xy_0,  \
                             ta_xxyy_xy_1,  \
                             ta_xxyy_yy_0,  \
                             ta_xxyy_yy_1,  \
                             ta_xxyyz_xx_0, \
                             ta_xxyyz_xy_0, \
                             ta_xxyyz_xz_0, \
                             ta_xxyyz_yy_0, \
                             ta_xxyyz_yz_0, \
                             ta_xxyyz_zz_0, \
                             ta_xxyz_xz_0,  \
                             ta_xxyz_xz_1,  \
                             ta_xxz_xz_0,   \
                             ta_xxz_xz_1,   \
                             ta_xyyz_yz_0,  \
                             ta_xyyz_yz_1,  \
                             ta_xyyz_zz_0,  \
                             ta_xyyz_zz_1,  \
                             ta_yyz_yz_0,   \
                             ta_yyz_yz_1,   \
                             ta_yyz_zz_0,   \
                             ta_yyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyz_xx_0[i] = ta_xxyy_xx_0[i] * pa_z[i] - ta_xxyy_xx_1[i] * pc_z[i];

        ta_xxyyz_xy_0[i] = ta_xxyy_xy_0[i] * pa_z[i] - ta_xxyy_xy_1[i] * pc_z[i];

        ta_xxyyz_xz_0[i] = ta_xxz_xz_0[i] * fe_0 - ta_xxz_xz_1[i] * fe_0 + ta_xxyz_xz_0[i] * pa_y[i] - ta_xxyz_xz_1[i] * pc_y[i];

        ta_xxyyz_yy_0[i] = ta_xxyy_yy_0[i] * pa_z[i] - ta_xxyy_yy_1[i] * pc_z[i];

        ta_xxyyz_yz_0[i] = ta_yyz_yz_0[i] * fe_0 - ta_yyz_yz_1[i] * fe_0 + ta_xyyz_yz_0[i] * pa_x[i] - ta_xyyz_yz_1[i] * pc_x[i];

        ta_xxyyz_zz_0[i] = ta_yyz_zz_0[i] * fe_0 - ta_yyz_zz_1[i] * fe_0 + ta_xyyz_zz_0[i] * pa_x[i] - ta_xyyz_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : HD

    auto ta_xxyzz_xx_0 = pbuffer.data(idx_npot_0_hd + 48);

    auto ta_xxyzz_xy_0 = pbuffer.data(idx_npot_0_hd + 49);

    auto ta_xxyzz_xz_0 = pbuffer.data(idx_npot_0_hd + 50);

    auto ta_xxyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 51);

    auto ta_xxyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 52);

    auto ta_xxyzz_zz_0 = pbuffer.data(idx_npot_0_hd + 53);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxyzz_xx_0, \
                             ta_xxyzz_xy_0, \
                             ta_xxyzz_xz_0, \
                             ta_xxyzz_yy_0, \
                             ta_xxyzz_yz_0, \
                             ta_xxyzz_zz_0, \
                             ta_xxzz_x_0,   \
                             ta_xxzz_x_1,   \
                             ta_xxzz_xx_0,  \
                             ta_xxzz_xx_1,  \
                             ta_xxzz_xy_0,  \
                             ta_xxzz_xy_1,  \
                             ta_xxzz_xz_0,  \
                             ta_xxzz_xz_1,  \
                             ta_xxzz_zz_0,  \
                             ta_xxzz_zz_1,  \
                             ta_xyzz_yy_0,  \
                             ta_xyzz_yy_1,  \
                             ta_xyzz_yz_0,  \
                             ta_xyzz_yz_1,  \
                             ta_yzz_yy_0,   \
                             ta_yzz_yy_1,   \
                             ta_yzz_yz_0,   \
                             ta_yzz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzz_xx_0[i] = ta_xxzz_xx_0[i] * pa_y[i] - ta_xxzz_xx_1[i] * pc_y[i];

        ta_xxyzz_xy_0[i] = ta_xxzz_x_0[i] * fe_0 - ta_xxzz_x_1[i] * fe_0 + ta_xxzz_xy_0[i] * pa_y[i] - ta_xxzz_xy_1[i] * pc_y[i];

        ta_xxyzz_xz_0[i] = ta_xxzz_xz_0[i] * pa_y[i] - ta_xxzz_xz_1[i] * pc_y[i];

        ta_xxyzz_yy_0[i] = ta_yzz_yy_0[i] * fe_0 - ta_yzz_yy_1[i] * fe_0 + ta_xyzz_yy_0[i] * pa_x[i] - ta_xyzz_yy_1[i] * pc_x[i];

        ta_xxyzz_yz_0[i] = ta_yzz_yz_0[i] * fe_0 - ta_yzz_yz_1[i] * fe_0 + ta_xyzz_yz_0[i] * pa_x[i] - ta_xyzz_yz_1[i] * pc_x[i];

        ta_xxyzz_zz_0[i] = ta_xxzz_zz_0[i] * pa_y[i] - ta_xxzz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : HD

    auto ta_xxzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 54);

    auto ta_xxzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 55);

    auto ta_xxzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 56);

    auto ta_xxzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 57);

    auto ta_xxzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 58);

    auto ta_xxzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxz_xx_0,   \
                             ta_xxz_xx_1,   \
                             ta_xxz_xy_0,   \
                             ta_xxz_xy_1,   \
                             ta_xxzz_xx_0,  \
                             ta_xxzz_xx_1,  \
                             ta_xxzz_xy_0,  \
                             ta_xxzz_xy_1,  \
                             ta_xxzzz_xx_0, \
                             ta_xxzzz_xy_0, \
                             ta_xxzzz_xz_0, \
                             ta_xxzzz_yy_0, \
                             ta_xxzzz_yz_0, \
                             ta_xxzzz_zz_0, \
                             ta_xzzz_xz_0,  \
                             ta_xzzz_xz_1,  \
                             ta_xzzz_yy_0,  \
                             ta_xzzz_yy_1,  \
                             ta_xzzz_yz_0,  \
                             ta_xzzz_yz_1,  \
                             ta_xzzz_z_0,   \
                             ta_xzzz_z_1,   \
                             ta_xzzz_zz_0,  \
                             ta_xzzz_zz_1,  \
                             ta_zzz_xz_0,   \
                             ta_zzz_xz_1,   \
                             ta_zzz_yy_0,   \
                             ta_zzz_yy_1,   \
                             ta_zzz_yz_0,   \
                             ta_zzz_yz_1,   \
                             ta_zzz_zz_0,   \
                             ta_zzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzz_xx_0[i] = 2.0 * ta_xxz_xx_0[i] * fe_0 - 2.0 * ta_xxz_xx_1[i] * fe_0 + ta_xxzz_xx_0[i] * pa_z[i] - ta_xxzz_xx_1[i] * pc_z[i];

        ta_xxzzz_xy_0[i] = 2.0 * ta_xxz_xy_0[i] * fe_0 - 2.0 * ta_xxz_xy_1[i] * fe_0 + ta_xxzz_xy_0[i] * pa_z[i] - ta_xxzz_xy_1[i] * pc_z[i];

        ta_xxzzz_xz_0[i] = ta_zzz_xz_0[i] * fe_0 - ta_zzz_xz_1[i] * fe_0 + ta_xzzz_z_0[i] * fe_0 - ta_xzzz_z_1[i] * fe_0 + ta_xzzz_xz_0[i] * pa_x[i] -
                           ta_xzzz_xz_1[i] * pc_x[i];

        ta_xxzzz_yy_0[i] = ta_zzz_yy_0[i] * fe_0 - ta_zzz_yy_1[i] * fe_0 + ta_xzzz_yy_0[i] * pa_x[i] - ta_xzzz_yy_1[i] * pc_x[i];

        ta_xxzzz_yz_0[i] = ta_zzz_yz_0[i] * fe_0 - ta_zzz_yz_1[i] * fe_0 + ta_xzzz_yz_0[i] * pa_x[i] - ta_xzzz_yz_1[i] * pc_x[i];

        ta_xxzzz_zz_0[i] = ta_zzz_zz_0[i] * fe_0 - ta_zzz_zz_1[i] * fe_0 + ta_xzzz_zz_0[i] * pa_x[i] - ta_xzzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : HD

    auto ta_xyyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 60);

    auto ta_xyyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 61);

    auto ta_xyyyy_xz_0 = pbuffer.data(idx_npot_0_hd + 62);

    auto ta_xyyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 63);

    auto ta_xyyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 64);

    auto ta_xyyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 65);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyyyy_xx_0, \
                             ta_xyyyy_xy_0, \
                             ta_xyyyy_xz_0, \
                             ta_xyyyy_yy_0, \
                             ta_xyyyy_yz_0, \
                             ta_xyyyy_zz_0, \
                             ta_yyyy_x_0,   \
                             ta_yyyy_x_1,   \
                             ta_yyyy_xx_0,  \
                             ta_yyyy_xx_1,  \
                             ta_yyyy_xy_0,  \
                             ta_yyyy_xy_1,  \
                             ta_yyyy_xz_0,  \
                             ta_yyyy_xz_1,  \
                             ta_yyyy_y_0,   \
                             ta_yyyy_y_1,   \
                             ta_yyyy_yy_0,  \
                             ta_yyyy_yy_1,  \
                             ta_yyyy_yz_0,  \
                             ta_yyyy_yz_1,  \
                             ta_yyyy_z_0,   \
                             ta_yyyy_z_1,   \
                             ta_yyyy_zz_0,  \
                             ta_yyyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyy_xx_0[i] = 2.0 * ta_yyyy_x_0[i] * fe_0 - 2.0 * ta_yyyy_x_1[i] * fe_0 + ta_yyyy_xx_0[i] * pa_x[i] - ta_yyyy_xx_1[i] * pc_x[i];

        ta_xyyyy_xy_0[i] = ta_yyyy_y_0[i] * fe_0 - ta_yyyy_y_1[i] * fe_0 + ta_yyyy_xy_0[i] * pa_x[i] - ta_yyyy_xy_1[i] * pc_x[i];

        ta_xyyyy_xz_0[i] = ta_yyyy_z_0[i] * fe_0 - ta_yyyy_z_1[i] * fe_0 + ta_yyyy_xz_0[i] * pa_x[i] - ta_yyyy_xz_1[i] * pc_x[i];

        ta_xyyyy_yy_0[i] = ta_yyyy_yy_0[i] * pa_x[i] - ta_yyyy_yy_1[i] * pc_x[i];

        ta_xyyyy_yz_0[i] = ta_yyyy_yz_0[i] * pa_x[i] - ta_yyyy_yz_1[i] * pc_x[i];

        ta_xyyyy_zz_0[i] = ta_yyyy_zz_0[i] * pa_x[i] - ta_yyyy_zz_1[i] * pc_x[i];
    }

    // Set up 66-72 components of targeted buffer : HD

    auto ta_xyyyz_xx_0 = pbuffer.data(idx_npot_0_hd + 66);

    auto ta_xyyyz_xy_0 = pbuffer.data(idx_npot_0_hd + 67);

    auto ta_xyyyz_xz_0 = pbuffer.data(idx_npot_0_hd + 68);

    auto ta_xyyyz_yy_0 = pbuffer.data(idx_npot_0_hd + 69);

    auto ta_xyyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 70);

    auto ta_xyyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 71);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xyyy_xx_0,  \
                             ta_xyyy_xx_1,  \
                             ta_xyyy_xy_0,  \
                             ta_xyyy_xy_1,  \
                             ta_xyyyz_xx_0, \
                             ta_xyyyz_xy_0, \
                             ta_xyyyz_xz_0, \
                             ta_xyyyz_yy_0, \
                             ta_xyyyz_yz_0, \
                             ta_xyyyz_zz_0, \
                             ta_yyyz_xz_0,  \
                             ta_yyyz_xz_1,  \
                             ta_yyyz_yy_0,  \
                             ta_yyyz_yy_1,  \
                             ta_yyyz_yz_0,  \
                             ta_yyyz_yz_1,  \
                             ta_yyyz_z_0,   \
                             ta_yyyz_z_1,   \
                             ta_yyyz_zz_0,  \
                             ta_yyyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyz_xx_0[i] = ta_xyyy_xx_0[i] * pa_z[i] - ta_xyyy_xx_1[i] * pc_z[i];

        ta_xyyyz_xy_0[i] = ta_xyyy_xy_0[i] * pa_z[i] - ta_xyyy_xy_1[i] * pc_z[i];

        ta_xyyyz_xz_0[i] = ta_yyyz_z_0[i] * fe_0 - ta_yyyz_z_1[i] * fe_0 + ta_yyyz_xz_0[i] * pa_x[i] - ta_yyyz_xz_1[i] * pc_x[i];

        ta_xyyyz_yy_0[i] = ta_yyyz_yy_0[i] * pa_x[i] - ta_yyyz_yy_1[i] * pc_x[i];

        ta_xyyyz_yz_0[i] = ta_yyyz_yz_0[i] * pa_x[i] - ta_yyyz_yz_1[i] * pc_x[i];

        ta_xyyyz_zz_0[i] = ta_yyyz_zz_0[i] * pa_x[i] - ta_yyyz_zz_1[i] * pc_x[i];
    }

    // Set up 72-78 components of targeted buffer : HD

    auto ta_xyyzz_xx_0 = pbuffer.data(idx_npot_0_hd + 72);

    auto ta_xyyzz_xy_0 = pbuffer.data(idx_npot_0_hd + 73);

    auto ta_xyyzz_xz_0 = pbuffer.data(idx_npot_0_hd + 74);

    auto ta_xyyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 75);

    auto ta_xyyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 76);

    auto ta_xyyzz_zz_0 = pbuffer.data(idx_npot_0_hd + 77);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyyzz_xx_0, \
                             ta_xyyzz_xy_0, \
                             ta_xyyzz_xz_0, \
                             ta_xyyzz_yy_0, \
                             ta_xyyzz_yz_0, \
                             ta_xyyzz_zz_0, \
                             ta_yyzz_x_0,   \
                             ta_yyzz_x_1,   \
                             ta_yyzz_xx_0,  \
                             ta_yyzz_xx_1,  \
                             ta_yyzz_xy_0,  \
                             ta_yyzz_xy_1,  \
                             ta_yyzz_xz_0,  \
                             ta_yyzz_xz_1,  \
                             ta_yyzz_y_0,   \
                             ta_yyzz_y_1,   \
                             ta_yyzz_yy_0,  \
                             ta_yyzz_yy_1,  \
                             ta_yyzz_yz_0,  \
                             ta_yyzz_yz_1,  \
                             ta_yyzz_z_0,   \
                             ta_yyzz_z_1,   \
                             ta_yyzz_zz_0,  \
                             ta_yyzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzz_xx_0[i] = 2.0 * ta_yyzz_x_0[i] * fe_0 - 2.0 * ta_yyzz_x_1[i] * fe_0 + ta_yyzz_xx_0[i] * pa_x[i] - ta_yyzz_xx_1[i] * pc_x[i];

        ta_xyyzz_xy_0[i] = ta_yyzz_y_0[i] * fe_0 - ta_yyzz_y_1[i] * fe_0 + ta_yyzz_xy_0[i] * pa_x[i] - ta_yyzz_xy_1[i] * pc_x[i];

        ta_xyyzz_xz_0[i] = ta_yyzz_z_0[i] * fe_0 - ta_yyzz_z_1[i] * fe_0 + ta_yyzz_xz_0[i] * pa_x[i] - ta_yyzz_xz_1[i] * pc_x[i];

        ta_xyyzz_yy_0[i] = ta_yyzz_yy_0[i] * pa_x[i] - ta_yyzz_yy_1[i] * pc_x[i];

        ta_xyyzz_yz_0[i] = ta_yyzz_yz_0[i] * pa_x[i] - ta_yyzz_yz_1[i] * pc_x[i];

        ta_xyyzz_zz_0[i] = ta_yyzz_zz_0[i] * pa_x[i] - ta_yyzz_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : HD

    auto ta_xyzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 78);

    auto ta_xyzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 79);

    auto ta_xyzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 80);

    auto ta_xyzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 81);

    auto ta_xyzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 82);

    auto ta_xyzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 83);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xyzzz_xx_0, \
                             ta_xyzzz_xy_0, \
                             ta_xyzzz_xz_0, \
                             ta_xyzzz_yy_0, \
                             ta_xyzzz_yz_0, \
                             ta_xyzzz_zz_0, \
                             ta_xzzz_xx_0,  \
                             ta_xzzz_xx_1,  \
                             ta_xzzz_xz_0,  \
                             ta_xzzz_xz_1,  \
                             ta_yzzz_xy_0,  \
                             ta_yzzz_xy_1,  \
                             ta_yzzz_y_0,   \
                             ta_yzzz_y_1,   \
                             ta_yzzz_yy_0,  \
                             ta_yzzz_yy_1,  \
                             ta_yzzz_yz_0,  \
                             ta_yzzz_yz_1,  \
                             ta_yzzz_zz_0,  \
                             ta_yzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzz_xx_0[i] = ta_xzzz_xx_0[i] * pa_y[i] - ta_xzzz_xx_1[i] * pc_y[i];

        ta_xyzzz_xy_0[i] = ta_yzzz_y_0[i] * fe_0 - ta_yzzz_y_1[i] * fe_0 + ta_yzzz_xy_0[i] * pa_x[i] - ta_yzzz_xy_1[i] * pc_x[i];

        ta_xyzzz_xz_0[i] = ta_xzzz_xz_0[i] * pa_y[i] - ta_xzzz_xz_1[i] * pc_y[i];

        ta_xyzzz_yy_0[i] = ta_yzzz_yy_0[i] * pa_x[i] - ta_yzzz_yy_1[i] * pc_x[i];

        ta_xyzzz_yz_0[i] = ta_yzzz_yz_0[i] * pa_x[i] - ta_yzzz_yz_1[i] * pc_x[i];

        ta_xyzzz_zz_0[i] = ta_yzzz_zz_0[i] * pa_x[i] - ta_yzzz_zz_1[i] * pc_x[i];
    }

    // Set up 84-90 components of targeted buffer : HD

    auto ta_xzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 84);

    auto ta_xzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 85);

    auto ta_xzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 86);

    auto ta_xzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 87);

    auto ta_xzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 88);

    auto ta_xzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 89);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xzzzz_xx_0, \
                             ta_xzzzz_xy_0, \
                             ta_xzzzz_xz_0, \
                             ta_xzzzz_yy_0, \
                             ta_xzzzz_yz_0, \
                             ta_xzzzz_zz_0, \
                             ta_zzzz_x_0,   \
                             ta_zzzz_x_1,   \
                             ta_zzzz_xx_0,  \
                             ta_zzzz_xx_1,  \
                             ta_zzzz_xy_0,  \
                             ta_zzzz_xy_1,  \
                             ta_zzzz_xz_0,  \
                             ta_zzzz_xz_1,  \
                             ta_zzzz_y_0,   \
                             ta_zzzz_y_1,   \
                             ta_zzzz_yy_0,  \
                             ta_zzzz_yy_1,  \
                             ta_zzzz_yz_0,  \
                             ta_zzzz_yz_1,  \
                             ta_zzzz_z_0,   \
                             ta_zzzz_z_1,   \
                             ta_zzzz_zz_0,  \
                             ta_zzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzz_xx_0[i] = 2.0 * ta_zzzz_x_0[i] * fe_0 - 2.0 * ta_zzzz_x_1[i] * fe_0 + ta_zzzz_xx_0[i] * pa_x[i] - ta_zzzz_xx_1[i] * pc_x[i];

        ta_xzzzz_xy_0[i] = ta_zzzz_y_0[i] * fe_0 - ta_zzzz_y_1[i] * fe_0 + ta_zzzz_xy_0[i] * pa_x[i] - ta_zzzz_xy_1[i] * pc_x[i];

        ta_xzzzz_xz_0[i] = ta_zzzz_z_0[i] * fe_0 - ta_zzzz_z_1[i] * fe_0 + ta_zzzz_xz_0[i] * pa_x[i] - ta_zzzz_xz_1[i] * pc_x[i];

        ta_xzzzz_yy_0[i] = ta_zzzz_yy_0[i] * pa_x[i] - ta_zzzz_yy_1[i] * pc_x[i];

        ta_xzzzz_yz_0[i] = ta_zzzz_yz_0[i] * pa_x[i] - ta_zzzz_yz_1[i] * pc_x[i];

        ta_xzzzz_zz_0[i] = ta_zzzz_zz_0[i] * pa_x[i] - ta_zzzz_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : HD

    auto ta_yyyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 90);

    auto ta_yyyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 91);

    auto ta_yyyyy_xz_0 = pbuffer.data(idx_npot_0_hd + 92);

    auto ta_yyyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 93);

    auto ta_yyyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 94);

    auto ta_yyyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 95);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yyy_xx_0,   \
                             ta_yyy_xx_1,   \
                             ta_yyy_xy_0,   \
                             ta_yyy_xy_1,   \
                             ta_yyy_xz_0,   \
                             ta_yyy_xz_1,   \
                             ta_yyy_yy_0,   \
                             ta_yyy_yy_1,   \
                             ta_yyy_yz_0,   \
                             ta_yyy_yz_1,   \
                             ta_yyy_zz_0,   \
                             ta_yyy_zz_1,   \
                             ta_yyyy_x_0,   \
                             ta_yyyy_x_1,   \
                             ta_yyyy_xx_0,  \
                             ta_yyyy_xx_1,  \
                             ta_yyyy_xy_0,  \
                             ta_yyyy_xy_1,  \
                             ta_yyyy_xz_0,  \
                             ta_yyyy_xz_1,  \
                             ta_yyyy_y_0,   \
                             ta_yyyy_y_1,   \
                             ta_yyyy_yy_0,  \
                             ta_yyyy_yy_1,  \
                             ta_yyyy_yz_0,  \
                             ta_yyyy_yz_1,  \
                             ta_yyyy_z_0,   \
                             ta_yyyy_z_1,   \
                             ta_yyyy_zz_0,  \
                             ta_yyyy_zz_1,  \
                             ta_yyyyy_xx_0, \
                             ta_yyyyy_xy_0, \
                             ta_yyyyy_xz_0, \
                             ta_yyyyy_yy_0, \
                             ta_yyyyy_yz_0, \
                             ta_yyyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyy_xx_0[i] = 4.0 * ta_yyy_xx_0[i] * fe_0 - 4.0 * ta_yyy_xx_1[i] * fe_0 + ta_yyyy_xx_0[i] * pa_y[i] - ta_yyyy_xx_1[i] * pc_y[i];

        ta_yyyyy_xy_0[i] = 4.0 * ta_yyy_xy_0[i] * fe_0 - 4.0 * ta_yyy_xy_1[i] * fe_0 + ta_yyyy_x_0[i] * fe_0 - ta_yyyy_x_1[i] * fe_0 +
                           ta_yyyy_xy_0[i] * pa_y[i] - ta_yyyy_xy_1[i] * pc_y[i];

        ta_yyyyy_xz_0[i] = 4.0 * ta_yyy_xz_0[i] * fe_0 - 4.0 * ta_yyy_xz_1[i] * fe_0 + ta_yyyy_xz_0[i] * pa_y[i] - ta_yyyy_xz_1[i] * pc_y[i];

        ta_yyyyy_yy_0[i] = 4.0 * ta_yyy_yy_0[i] * fe_0 - 4.0 * ta_yyy_yy_1[i] * fe_0 + 2.0 * ta_yyyy_y_0[i] * fe_0 - 2.0 * ta_yyyy_y_1[i] * fe_0 +
                           ta_yyyy_yy_0[i] * pa_y[i] - ta_yyyy_yy_1[i] * pc_y[i];

        ta_yyyyy_yz_0[i] = 4.0 * ta_yyy_yz_0[i] * fe_0 - 4.0 * ta_yyy_yz_1[i] * fe_0 + ta_yyyy_z_0[i] * fe_0 - ta_yyyy_z_1[i] * fe_0 +
                           ta_yyyy_yz_0[i] * pa_y[i] - ta_yyyy_yz_1[i] * pc_y[i];

        ta_yyyyy_zz_0[i] = 4.0 * ta_yyy_zz_0[i] * fe_0 - 4.0 * ta_yyy_zz_1[i] * fe_0 + ta_yyyy_zz_0[i] * pa_y[i] - ta_yyyy_zz_1[i] * pc_y[i];
    }

    // Set up 96-102 components of targeted buffer : HD

    auto ta_yyyyz_xx_0 = pbuffer.data(idx_npot_0_hd + 96);

    auto ta_yyyyz_xy_0 = pbuffer.data(idx_npot_0_hd + 97);

    auto ta_yyyyz_xz_0 = pbuffer.data(idx_npot_0_hd + 98);

    auto ta_yyyyz_yy_0 = pbuffer.data(idx_npot_0_hd + 99);

    auto ta_yyyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 100);

    auto ta_yyyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 101);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyyy_xx_0,  \
                             ta_yyyy_xx_1,  \
                             ta_yyyy_xy_0,  \
                             ta_yyyy_xy_1,  \
                             ta_yyyy_y_0,   \
                             ta_yyyy_y_1,   \
                             ta_yyyy_yy_0,  \
                             ta_yyyy_yy_1,  \
                             ta_yyyy_yz_0,  \
                             ta_yyyy_yz_1,  \
                             ta_yyyyz_xx_0, \
                             ta_yyyyz_xy_0, \
                             ta_yyyyz_xz_0, \
                             ta_yyyyz_yy_0, \
                             ta_yyyyz_yz_0, \
                             ta_yyyyz_zz_0, \
                             ta_yyyz_xz_0,  \
                             ta_yyyz_xz_1,  \
                             ta_yyyz_zz_0,  \
                             ta_yyyz_zz_1,  \
                             ta_yyz_xz_0,   \
                             ta_yyz_xz_1,   \
                             ta_yyz_zz_0,   \
                             ta_yyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyz_xx_0[i] = ta_yyyy_xx_0[i] * pa_z[i] - ta_yyyy_xx_1[i] * pc_z[i];

        ta_yyyyz_xy_0[i] = ta_yyyy_xy_0[i] * pa_z[i] - ta_yyyy_xy_1[i] * pc_z[i];

        ta_yyyyz_xz_0[i] = 3.0 * ta_yyz_xz_0[i] * fe_0 - 3.0 * ta_yyz_xz_1[i] * fe_0 + ta_yyyz_xz_0[i] * pa_y[i] - ta_yyyz_xz_1[i] * pc_y[i];

        ta_yyyyz_yy_0[i] = ta_yyyy_yy_0[i] * pa_z[i] - ta_yyyy_yy_1[i] * pc_z[i];

        ta_yyyyz_yz_0[i] = ta_yyyy_y_0[i] * fe_0 - ta_yyyy_y_1[i] * fe_0 + ta_yyyy_yz_0[i] * pa_z[i] - ta_yyyy_yz_1[i] * pc_z[i];

        ta_yyyyz_zz_0[i] = 3.0 * ta_yyz_zz_0[i] * fe_0 - 3.0 * ta_yyz_zz_1[i] * fe_0 + ta_yyyz_zz_0[i] * pa_y[i] - ta_yyyz_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : HD

    auto ta_yyyzz_xx_0 = pbuffer.data(idx_npot_0_hd + 102);

    auto ta_yyyzz_xy_0 = pbuffer.data(idx_npot_0_hd + 103);

    auto ta_yyyzz_xz_0 = pbuffer.data(idx_npot_0_hd + 104);

    auto ta_yyyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 105);

    auto ta_yyyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 106);

    auto ta_yyyzz_zz_0 = pbuffer.data(idx_npot_0_hd + 107);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyy_xy_0,   \
                             ta_yyy_xy_1,   \
                             ta_yyy_yy_0,   \
                             ta_yyy_yy_1,   \
                             ta_yyyz_xy_0,  \
                             ta_yyyz_xy_1,  \
                             ta_yyyz_yy_0,  \
                             ta_yyyz_yy_1,  \
                             ta_yyyzz_xx_0, \
                             ta_yyyzz_xy_0, \
                             ta_yyyzz_xz_0, \
                             ta_yyyzz_yy_0, \
                             ta_yyyzz_yz_0, \
                             ta_yyyzz_zz_0, \
                             ta_yyzz_xx_0,  \
                             ta_yyzz_xx_1,  \
                             ta_yyzz_xz_0,  \
                             ta_yyzz_xz_1,  \
                             ta_yyzz_yz_0,  \
                             ta_yyzz_yz_1,  \
                             ta_yyzz_z_0,   \
                             ta_yyzz_z_1,   \
                             ta_yyzz_zz_0,  \
                             ta_yyzz_zz_1,  \
                             ta_yzz_xx_0,   \
                             ta_yzz_xx_1,   \
                             ta_yzz_xz_0,   \
                             ta_yzz_xz_1,   \
                             ta_yzz_yz_0,   \
                             ta_yzz_yz_1,   \
                             ta_yzz_zz_0,   \
                             ta_yzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzz_xx_0[i] = 2.0 * ta_yzz_xx_0[i] * fe_0 - 2.0 * ta_yzz_xx_1[i] * fe_0 + ta_yyzz_xx_0[i] * pa_y[i] - ta_yyzz_xx_1[i] * pc_y[i];

        ta_yyyzz_xy_0[i] = ta_yyy_xy_0[i] * fe_0 - ta_yyy_xy_1[i] * fe_0 + ta_yyyz_xy_0[i] * pa_z[i] - ta_yyyz_xy_1[i] * pc_z[i];

        ta_yyyzz_xz_0[i] = 2.0 * ta_yzz_xz_0[i] * fe_0 - 2.0 * ta_yzz_xz_1[i] * fe_0 + ta_yyzz_xz_0[i] * pa_y[i] - ta_yyzz_xz_1[i] * pc_y[i];

        ta_yyyzz_yy_0[i] = ta_yyy_yy_0[i] * fe_0 - ta_yyy_yy_1[i] * fe_0 + ta_yyyz_yy_0[i] * pa_z[i] - ta_yyyz_yy_1[i] * pc_z[i];

        ta_yyyzz_yz_0[i] = 2.0 * ta_yzz_yz_0[i] * fe_0 - 2.0 * ta_yzz_yz_1[i] * fe_0 + ta_yyzz_z_0[i] * fe_0 - ta_yyzz_z_1[i] * fe_0 +
                           ta_yyzz_yz_0[i] * pa_y[i] - ta_yyzz_yz_1[i] * pc_y[i];

        ta_yyyzz_zz_0[i] = 2.0 * ta_yzz_zz_0[i] * fe_0 - 2.0 * ta_yzz_zz_1[i] * fe_0 + ta_yyzz_zz_0[i] * pa_y[i] - ta_yyzz_zz_1[i] * pc_y[i];
    }

    // Set up 108-114 components of targeted buffer : HD

    auto ta_yyzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 108);

    auto ta_yyzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 109);

    auto ta_yyzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 110);

    auto ta_yyzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 111);

    auto ta_yyzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 112);

    auto ta_yyzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 113);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyz_xy_0,   \
                             ta_yyz_xy_1,   \
                             ta_yyz_yy_0,   \
                             ta_yyz_yy_1,   \
                             ta_yyzz_xy_0,  \
                             ta_yyzz_xy_1,  \
                             ta_yyzz_yy_0,  \
                             ta_yyzz_yy_1,  \
                             ta_yyzzz_xx_0, \
                             ta_yyzzz_xy_0, \
                             ta_yyzzz_xz_0, \
                             ta_yyzzz_yy_0, \
                             ta_yyzzz_yz_0, \
                             ta_yyzzz_zz_0, \
                             ta_yzzz_xx_0,  \
                             ta_yzzz_xx_1,  \
                             ta_yzzz_xz_0,  \
                             ta_yzzz_xz_1,  \
                             ta_yzzz_yz_0,  \
                             ta_yzzz_yz_1,  \
                             ta_yzzz_z_0,   \
                             ta_yzzz_z_1,   \
                             ta_yzzz_zz_0,  \
                             ta_yzzz_zz_1,  \
                             ta_zzz_xx_0,   \
                             ta_zzz_xx_1,   \
                             ta_zzz_xz_0,   \
                             ta_zzz_xz_1,   \
                             ta_zzz_yz_0,   \
                             ta_zzz_yz_1,   \
                             ta_zzz_zz_0,   \
                             ta_zzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzz_xx_0[i] = ta_zzz_xx_0[i] * fe_0 - ta_zzz_xx_1[i] * fe_0 + ta_yzzz_xx_0[i] * pa_y[i] - ta_yzzz_xx_1[i] * pc_y[i];

        ta_yyzzz_xy_0[i] = 2.0 * ta_yyz_xy_0[i] * fe_0 - 2.0 * ta_yyz_xy_1[i] * fe_0 + ta_yyzz_xy_0[i] * pa_z[i] - ta_yyzz_xy_1[i] * pc_z[i];

        ta_yyzzz_xz_0[i] = ta_zzz_xz_0[i] * fe_0 - ta_zzz_xz_1[i] * fe_0 + ta_yzzz_xz_0[i] * pa_y[i] - ta_yzzz_xz_1[i] * pc_y[i];

        ta_yyzzz_yy_0[i] = 2.0 * ta_yyz_yy_0[i] * fe_0 - 2.0 * ta_yyz_yy_1[i] * fe_0 + ta_yyzz_yy_0[i] * pa_z[i] - ta_yyzz_yy_1[i] * pc_z[i];

        ta_yyzzz_yz_0[i] = ta_zzz_yz_0[i] * fe_0 - ta_zzz_yz_1[i] * fe_0 + ta_yzzz_z_0[i] * fe_0 - ta_yzzz_z_1[i] * fe_0 + ta_yzzz_yz_0[i] * pa_y[i] -
                           ta_yzzz_yz_1[i] * pc_y[i];

        ta_yyzzz_zz_0[i] = ta_zzz_zz_0[i] * fe_0 - ta_zzz_zz_1[i] * fe_0 + ta_yzzz_zz_0[i] * pa_y[i] - ta_yzzz_zz_1[i] * pc_y[i];
    }

    // Set up 114-120 components of targeted buffer : HD

    auto ta_yzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 114);

    auto ta_yzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 115);

    auto ta_yzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 116);

    auto ta_yzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 117);

    auto ta_yzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 118);

    auto ta_yzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 119);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yzzzz_xx_0, \
                             ta_yzzzz_xy_0, \
                             ta_yzzzz_xz_0, \
                             ta_yzzzz_yy_0, \
                             ta_yzzzz_yz_0, \
                             ta_yzzzz_zz_0, \
                             ta_zzzz_x_0,   \
                             ta_zzzz_x_1,   \
                             ta_zzzz_xx_0,  \
                             ta_zzzz_xx_1,  \
                             ta_zzzz_xy_0,  \
                             ta_zzzz_xy_1,  \
                             ta_zzzz_xz_0,  \
                             ta_zzzz_xz_1,  \
                             ta_zzzz_y_0,   \
                             ta_zzzz_y_1,   \
                             ta_zzzz_yy_0,  \
                             ta_zzzz_yy_1,  \
                             ta_zzzz_yz_0,  \
                             ta_zzzz_yz_1,  \
                             ta_zzzz_z_0,   \
                             ta_zzzz_z_1,   \
                             ta_zzzz_zz_0,  \
                             ta_zzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzz_xx_0[i] = ta_zzzz_xx_0[i] * pa_y[i] - ta_zzzz_xx_1[i] * pc_y[i];

        ta_yzzzz_xy_0[i] = ta_zzzz_x_0[i] * fe_0 - ta_zzzz_x_1[i] * fe_0 + ta_zzzz_xy_0[i] * pa_y[i] - ta_zzzz_xy_1[i] * pc_y[i];

        ta_yzzzz_xz_0[i] = ta_zzzz_xz_0[i] * pa_y[i] - ta_zzzz_xz_1[i] * pc_y[i];

        ta_yzzzz_yy_0[i] = 2.0 * ta_zzzz_y_0[i] * fe_0 - 2.0 * ta_zzzz_y_1[i] * fe_0 + ta_zzzz_yy_0[i] * pa_y[i] - ta_zzzz_yy_1[i] * pc_y[i];

        ta_yzzzz_yz_0[i] = ta_zzzz_z_0[i] * fe_0 - ta_zzzz_z_1[i] * fe_0 + ta_zzzz_yz_0[i] * pa_y[i] - ta_zzzz_yz_1[i] * pc_y[i];

        ta_yzzzz_zz_0[i] = ta_zzzz_zz_0[i] * pa_y[i] - ta_zzzz_zz_1[i] * pc_y[i];
    }

    // Set up 120-126 components of targeted buffer : HD

    auto ta_zzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 120);

    auto ta_zzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 121);

    auto ta_zzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 122);

    auto ta_zzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 123);

    auto ta_zzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 124);

    auto ta_zzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 125);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta_zzz_xx_0,   \
                             ta_zzz_xx_1,   \
                             ta_zzz_xy_0,   \
                             ta_zzz_xy_1,   \
                             ta_zzz_xz_0,   \
                             ta_zzz_xz_1,   \
                             ta_zzz_yy_0,   \
                             ta_zzz_yy_1,   \
                             ta_zzz_yz_0,   \
                             ta_zzz_yz_1,   \
                             ta_zzz_zz_0,   \
                             ta_zzz_zz_1,   \
                             ta_zzzz_x_0,   \
                             ta_zzzz_x_1,   \
                             ta_zzzz_xx_0,  \
                             ta_zzzz_xx_1,  \
                             ta_zzzz_xy_0,  \
                             ta_zzzz_xy_1,  \
                             ta_zzzz_xz_0,  \
                             ta_zzzz_xz_1,  \
                             ta_zzzz_y_0,   \
                             ta_zzzz_y_1,   \
                             ta_zzzz_yy_0,  \
                             ta_zzzz_yy_1,  \
                             ta_zzzz_yz_0,  \
                             ta_zzzz_yz_1,  \
                             ta_zzzz_z_0,   \
                             ta_zzzz_z_1,   \
                             ta_zzzz_zz_0,  \
                             ta_zzzz_zz_1,  \
                             ta_zzzzz_xx_0, \
                             ta_zzzzz_xy_0, \
                             ta_zzzzz_xz_0, \
                             ta_zzzzz_yy_0, \
                             ta_zzzzz_yz_0, \
                             ta_zzzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzz_xx_0[i] = 4.0 * ta_zzz_xx_0[i] * fe_0 - 4.0 * ta_zzz_xx_1[i] * fe_0 + ta_zzzz_xx_0[i] * pa_z[i] - ta_zzzz_xx_1[i] * pc_z[i];

        ta_zzzzz_xy_0[i] = 4.0 * ta_zzz_xy_0[i] * fe_0 - 4.0 * ta_zzz_xy_1[i] * fe_0 + ta_zzzz_xy_0[i] * pa_z[i] - ta_zzzz_xy_1[i] * pc_z[i];

        ta_zzzzz_xz_0[i] = 4.0 * ta_zzz_xz_0[i] * fe_0 - 4.0 * ta_zzz_xz_1[i] * fe_0 + ta_zzzz_x_0[i] * fe_0 - ta_zzzz_x_1[i] * fe_0 +
                           ta_zzzz_xz_0[i] * pa_z[i] - ta_zzzz_xz_1[i] * pc_z[i];

        ta_zzzzz_yy_0[i] = 4.0 * ta_zzz_yy_0[i] * fe_0 - 4.0 * ta_zzz_yy_1[i] * fe_0 + ta_zzzz_yy_0[i] * pa_z[i] - ta_zzzz_yy_1[i] * pc_z[i];

        ta_zzzzz_yz_0[i] = 4.0 * ta_zzz_yz_0[i] * fe_0 - 4.0 * ta_zzz_yz_1[i] * fe_0 + ta_zzzz_y_0[i] * fe_0 - ta_zzzz_y_1[i] * fe_0 +
                           ta_zzzz_yz_0[i] * pa_z[i] - ta_zzzz_yz_1[i] * pc_z[i];

        ta_zzzzz_zz_0[i] = 4.0 * ta_zzz_zz_0[i] * fe_0 - 4.0 * ta_zzz_zz_1[i] * fe_0 + 2.0 * ta_zzzz_z_0[i] * fe_0 - 2.0 * ta_zzzz_z_1[i] * fe_0 +
                           ta_zzzz_zz_0[i] * pa_z[i] - ta_zzzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
