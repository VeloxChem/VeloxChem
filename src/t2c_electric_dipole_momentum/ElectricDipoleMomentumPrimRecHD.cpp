#include "ElectricDipoleMomentumPrimRecHD.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_hd(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_hd,
                                      const size_t idx_dip_fd,
                                      const size_t idx_dip_gp,
                                      const size_t idx_ovl_gd,
                                      const size_t idx_dip_gd,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FD

    auto tr_x_xxx_xx = pbuffer.data(idx_dip_fd);

    auto tr_x_xxx_xy = pbuffer.data(idx_dip_fd + 1);

    auto tr_x_xxx_xz = pbuffer.data(idx_dip_fd + 2);

    auto tr_x_xxx_yy = pbuffer.data(idx_dip_fd + 3);

    auto tr_x_xxx_yz = pbuffer.data(idx_dip_fd + 4);

    auto tr_x_xxx_zz = pbuffer.data(idx_dip_fd + 5);

    auto tr_x_xxy_xx = pbuffer.data(idx_dip_fd + 6);

    auto tr_x_xxy_xy = pbuffer.data(idx_dip_fd + 7);

    auto tr_x_xxy_xz = pbuffer.data(idx_dip_fd + 8);

    auto tr_x_xxy_zz = pbuffer.data(idx_dip_fd + 11);

    auto tr_x_xxz_xx = pbuffer.data(idx_dip_fd + 12);

    auto tr_x_xxz_xy = pbuffer.data(idx_dip_fd + 13);

    auto tr_x_xxz_xz = pbuffer.data(idx_dip_fd + 14);

    auto tr_x_xxz_yy = pbuffer.data(idx_dip_fd + 15);

    auto tr_x_xxz_zz = pbuffer.data(idx_dip_fd + 17);

    auto tr_x_xyy_xx = pbuffer.data(idx_dip_fd + 18);

    auto tr_x_xyy_xy = pbuffer.data(idx_dip_fd + 19);

    auto tr_x_xyy_xz = pbuffer.data(idx_dip_fd + 20);

    auto tr_x_xyy_yy = pbuffer.data(idx_dip_fd + 21);

    auto tr_x_xyy_yz = pbuffer.data(idx_dip_fd + 22);

    auto tr_x_xyz_xz = pbuffer.data(idx_dip_fd + 26);

    auto tr_x_xzz_xx = pbuffer.data(idx_dip_fd + 30);

    auto tr_x_xzz_xy = pbuffer.data(idx_dip_fd + 31);

    auto tr_x_xzz_xz = pbuffer.data(idx_dip_fd + 32);

    auto tr_x_xzz_yz = pbuffer.data(idx_dip_fd + 34);

    auto tr_x_xzz_zz = pbuffer.data(idx_dip_fd + 35);

    auto tr_x_yyy_xx = pbuffer.data(idx_dip_fd + 36);

    auto tr_x_yyy_xy = pbuffer.data(idx_dip_fd + 37);

    auto tr_x_yyy_xz = pbuffer.data(idx_dip_fd + 38);

    auto tr_x_yyy_yy = pbuffer.data(idx_dip_fd + 39);

    auto tr_x_yyy_yz = pbuffer.data(idx_dip_fd + 40);

    auto tr_x_yyy_zz = pbuffer.data(idx_dip_fd + 41);

    auto tr_x_yyz_xy = pbuffer.data(idx_dip_fd + 43);

    auto tr_x_yyz_xz = pbuffer.data(idx_dip_fd + 44);

    auto tr_x_yyz_yy = pbuffer.data(idx_dip_fd + 45);

    auto tr_x_yyz_zz = pbuffer.data(idx_dip_fd + 47);

    auto tr_x_yzz_xx = pbuffer.data(idx_dip_fd + 48);

    auto tr_x_yzz_xz = pbuffer.data(idx_dip_fd + 50);

    auto tr_x_yzz_yz = pbuffer.data(idx_dip_fd + 52);

    auto tr_x_yzz_zz = pbuffer.data(idx_dip_fd + 53);

    auto tr_x_zzz_xx = pbuffer.data(idx_dip_fd + 54);

    auto tr_x_zzz_xy = pbuffer.data(idx_dip_fd + 55);

    auto tr_x_zzz_xz = pbuffer.data(idx_dip_fd + 56);

    auto tr_x_zzz_yy = pbuffer.data(idx_dip_fd + 57);

    auto tr_x_zzz_yz = pbuffer.data(idx_dip_fd + 58);

    auto tr_x_zzz_zz = pbuffer.data(idx_dip_fd + 59);

    auto tr_y_xxx_xx = pbuffer.data(idx_dip_fd + 60);

    auto tr_y_xxx_xy = pbuffer.data(idx_dip_fd + 61);

    auto tr_y_xxx_xz = pbuffer.data(idx_dip_fd + 62);

    auto tr_y_xxx_yy = pbuffer.data(idx_dip_fd + 63);

    auto tr_y_xxx_yz = pbuffer.data(idx_dip_fd + 64);

    auto tr_y_xxx_zz = pbuffer.data(idx_dip_fd + 65);

    auto tr_y_xxy_xy = pbuffer.data(idx_dip_fd + 67);

    auto tr_y_xxy_yy = pbuffer.data(idx_dip_fd + 69);

    auto tr_y_xxy_yz = pbuffer.data(idx_dip_fd + 70);

    auto tr_y_xxy_zz = pbuffer.data(idx_dip_fd + 71);

    auto tr_y_xxz_xx = pbuffer.data(idx_dip_fd + 72);

    auto tr_y_xxz_xy = pbuffer.data(idx_dip_fd + 73);

    auto tr_y_xxz_yz = pbuffer.data(idx_dip_fd + 76);

    auto tr_y_xxz_zz = pbuffer.data(idx_dip_fd + 77);

    auto tr_y_xyy_xx = pbuffer.data(idx_dip_fd + 78);

    auto tr_y_xyy_xy = pbuffer.data(idx_dip_fd + 79);

    auto tr_y_xyy_xz = pbuffer.data(idx_dip_fd + 80);

    auto tr_y_xyy_yy = pbuffer.data(idx_dip_fd + 81);

    auto tr_y_xyy_yz = pbuffer.data(idx_dip_fd + 82);

    auto tr_y_xyy_zz = pbuffer.data(idx_dip_fd + 83);

    auto tr_y_xyz_yz = pbuffer.data(idx_dip_fd + 88);

    auto tr_y_xyz_zz = pbuffer.data(idx_dip_fd + 89);

    auto tr_y_xzz_xz = pbuffer.data(idx_dip_fd + 92);

    auto tr_y_xzz_yy = pbuffer.data(idx_dip_fd + 93);

    auto tr_y_xzz_yz = pbuffer.data(idx_dip_fd + 94);

    auto tr_y_xzz_zz = pbuffer.data(idx_dip_fd + 95);

    auto tr_y_yyy_xx = pbuffer.data(idx_dip_fd + 96);

    auto tr_y_yyy_xy = pbuffer.data(idx_dip_fd + 97);

    auto tr_y_yyy_xz = pbuffer.data(idx_dip_fd + 98);

    auto tr_y_yyy_yy = pbuffer.data(idx_dip_fd + 99);

    auto tr_y_yyy_yz = pbuffer.data(idx_dip_fd + 100);

    auto tr_y_yyy_zz = pbuffer.data(idx_dip_fd + 101);

    auto tr_y_yyz_xx = pbuffer.data(idx_dip_fd + 102);

    auto tr_y_yyz_xy = pbuffer.data(idx_dip_fd + 103);

    auto tr_y_yyz_yy = pbuffer.data(idx_dip_fd + 105);

    auto tr_y_yyz_yz = pbuffer.data(idx_dip_fd + 106);

    auto tr_y_yyz_zz = pbuffer.data(idx_dip_fd + 107);

    auto tr_y_yzz_xy = pbuffer.data(idx_dip_fd + 109);

    auto tr_y_yzz_xz = pbuffer.data(idx_dip_fd + 110);

    auto tr_y_yzz_yy = pbuffer.data(idx_dip_fd + 111);

    auto tr_y_yzz_yz = pbuffer.data(idx_dip_fd + 112);

    auto tr_y_yzz_zz = pbuffer.data(idx_dip_fd + 113);

    auto tr_y_zzz_xx = pbuffer.data(idx_dip_fd + 114);

    auto tr_y_zzz_xy = pbuffer.data(idx_dip_fd + 115);

    auto tr_y_zzz_xz = pbuffer.data(idx_dip_fd + 116);

    auto tr_y_zzz_yy = pbuffer.data(idx_dip_fd + 117);

    auto tr_y_zzz_yz = pbuffer.data(idx_dip_fd + 118);

    auto tr_y_zzz_zz = pbuffer.data(idx_dip_fd + 119);

    auto tr_z_xxx_xx = pbuffer.data(idx_dip_fd + 120);

    auto tr_z_xxx_xy = pbuffer.data(idx_dip_fd + 121);

    auto tr_z_xxx_xz = pbuffer.data(idx_dip_fd + 122);

    auto tr_z_xxx_yy = pbuffer.data(idx_dip_fd + 123);

    auto tr_z_xxx_yz = pbuffer.data(idx_dip_fd + 124);

    auto tr_z_xxx_zz = pbuffer.data(idx_dip_fd + 125);

    auto tr_z_xxy_xx = pbuffer.data(idx_dip_fd + 126);

    auto tr_z_xxy_xz = pbuffer.data(idx_dip_fd + 128);

    auto tr_z_xxy_yy = pbuffer.data(idx_dip_fd + 129);

    auto tr_z_xxy_yz = pbuffer.data(idx_dip_fd + 130);

    auto tr_z_xxz_xx = pbuffer.data(idx_dip_fd + 132);

    auto tr_z_xxz_xz = pbuffer.data(idx_dip_fd + 134);

    auto tr_z_xxz_yy = pbuffer.data(idx_dip_fd + 135);

    auto tr_z_xxz_yz = pbuffer.data(idx_dip_fd + 136);

    auto tr_z_xxz_zz = pbuffer.data(idx_dip_fd + 137);

    auto tr_z_xyy_xy = pbuffer.data(idx_dip_fd + 139);

    auto tr_z_xyy_yy = pbuffer.data(idx_dip_fd + 141);

    auto tr_z_xyy_yz = pbuffer.data(idx_dip_fd + 142);

    auto tr_z_xyy_zz = pbuffer.data(idx_dip_fd + 143);

    auto tr_z_xyz_yy = pbuffer.data(idx_dip_fd + 147);

    auto tr_z_xyz_yz = pbuffer.data(idx_dip_fd + 148);

    auto tr_z_xzz_xx = pbuffer.data(idx_dip_fd + 150);

    auto tr_z_xzz_xy = pbuffer.data(idx_dip_fd + 151);

    auto tr_z_xzz_xz = pbuffer.data(idx_dip_fd + 152);

    auto tr_z_xzz_yy = pbuffer.data(idx_dip_fd + 153);

    auto tr_z_xzz_yz = pbuffer.data(idx_dip_fd + 154);

    auto tr_z_xzz_zz = pbuffer.data(idx_dip_fd + 155);

    auto tr_z_yyy_xx = pbuffer.data(idx_dip_fd + 156);

    auto tr_z_yyy_xy = pbuffer.data(idx_dip_fd + 157);

    auto tr_z_yyy_xz = pbuffer.data(idx_dip_fd + 158);

    auto tr_z_yyy_yy = pbuffer.data(idx_dip_fd + 159);

    auto tr_z_yyy_yz = pbuffer.data(idx_dip_fd + 160);

    auto tr_z_yyy_zz = pbuffer.data(idx_dip_fd + 161);

    auto tr_z_yyz_xx = pbuffer.data(idx_dip_fd + 162);

    auto tr_z_yyz_xz = pbuffer.data(idx_dip_fd + 164);

    auto tr_z_yyz_yy = pbuffer.data(idx_dip_fd + 165);

    auto tr_z_yyz_yz = pbuffer.data(idx_dip_fd + 166);

    auto tr_z_yyz_zz = pbuffer.data(idx_dip_fd + 167);

    auto tr_z_yzz_xx = pbuffer.data(idx_dip_fd + 168);

    auto tr_z_yzz_xy = pbuffer.data(idx_dip_fd + 169);

    auto tr_z_yzz_xz = pbuffer.data(idx_dip_fd + 170);

    auto tr_z_yzz_yy = pbuffer.data(idx_dip_fd + 171);

    auto tr_z_yzz_yz = pbuffer.data(idx_dip_fd + 172);

    auto tr_z_yzz_zz = pbuffer.data(idx_dip_fd + 173);

    auto tr_z_zzz_xx = pbuffer.data(idx_dip_fd + 174);

    auto tr_z_zzz_xy = pbuffer.data(idx_dip_fd + 175);

    auto tr_z_zzz_xz = pbuffer.data(idx_dip_fd + 176);

    auto tr_z_zzz_yy = pbuffer.data(idx_dip_fd + 177);

    auto tr_z_zzz_yz = pbuffer.data(idx_dip_fd + 178);

    auto tr_z_zzz_zz = pbuffer.data(idx_dip_fd + 179);

    // Set up components of auxiliary buffer : GP

    auto tr_x_xxxx_x = pbuffer.data(idx_dip_gp);

    auto tr_x_xxxx_y = pbuffer.data(idx_dip_gp + 1);

    auto tr_x_xxxx_z = pbuffer.data(idx_dip_gp + 2);

    auto tr_x_xxxy_x = pbuffer.data(idx_dip_gp + 3);

    auto tr_x_xxxz_x = pbuffer.data(idx_dip_gp + 6);

    auto tr_x_xxxz_z = pbuffer.data(idx_dip_gp + 8);

    auto tr_x_xxyy_x = pbuffer.data(idx_dip_gp + 9);

    auto tr_x_xxyy_y = pbuffer.data(idx_dip_gp + 10);

    auto tr_x_xxzz_x = pbuffer.data(idx_dip_gp + 15);

    auto tr_x_xxzz_y = pbuffer.data(idx_dip_gp + 16);

    auto tr_x_xxzz_z = pbuffer.data(idx_dip_gp + 17);

    auto tr_x_xzzz_x = pbuffer.data(idx_dip_gp + 27);

    auto tr_x_yyyy_x = pbuffer.data(idx_dip_gp + 30);

    auto tr_x_yyyy_y = pbuffer.data(idx_dip_gp + 31);

    auto tr_x_yyyy_z = pbuffer.data(idx_dip_gp + 32);

    auto tr_x_yyzz_z = pbuffer.data(idx_dip_gp + 38);

    auto tr_x_yzzz_z = pbuffer.data(idx_dip_gp + 41);

    auto tr_x_zzzz_x = pbuffer.data(idx_dip_gp + 42);

    auto tr_x_zzzz_y = pbuffer.data(idx_dip_gp + 43);

    auto tr_x_zzzz_z = pbuffer.data(idx_dip_gp + 44);

    auto tr_y_xxxx_x = pbuffer.data(idx_dip_gp + 45);

    auto tr_y_xxxx_y = pbuffer.data(idx_dip_gp + 46);

    auto tr_y_xxxx_z = pbuffer.data(idx_dip_gp + 47);

    auto tr_y_xxxy_y = pbuffer.data(idx_dip_gp + 49);

    auto tr_y_xxyy_x = pbuffer.data(idx_dip_gp + 54);

    auto tr_y_xxyy_y = pbuffer.data(idx_dip_gp + 55);

    auto tr_y_xxyy_z = pbuffer.data(idx_dip_gp + 56);

    auto tr_y_xxzz_z = pbuffer.data(idx_dip_gp + 62);

    auto tr_y_xyyy_x = pbuffer.data(idx_dip_gp + 63);

    auto tr_y_xyyy_y = pbuffer.data(idx_dip_gp + 64);

    auto tr_y_xyyy_z = pbuffer.data(idx_dip_gp + 65);

    auto tr_y_xzzz_z = pbuffer.data(idx_dip_gp + 74);

    auto tr_y_yyyy_x = pbuffer.data(idx_dip_gp + 75);

    auto tr_y_yyyy_y = pbuffer.data(idx_dip_gp + 76);

    auto tr_y_yyyy_z = pbuffer.data(idx_dip_gp + 77);

    auto tr_y_yyyz_y = pbuffer.data(idx_dip_gp + 79);

    auto tr_y_yyyz_z = pbuffer.data(idx_dip_gp + 80);

    auto tr_y_yyzz_x = pbuffer.data(idx_dip_gp + 81);

    auto tr_y_yyzz_y = pbuffer.data(idx_dip_gp + 82);

    auto tr_y_yyzz_z = pbuffer.data(idx_dip_gp + 83);

    auto tr_y_yzzz_x = pbuffer.data(idx_dip_gp + 84);

    auto tr_y_yzzz_y = pbuffer.data(idx_dip_gp + 85);

    auto tr_y_yzzz_z = pbuffer.data(idx_dip_gp + 86);

    auto tr_y_zzzz_x = pbuffer.data(idx_dip_gp + 87);

    auto tr_y_zzzz_y = pbuffer.data(idx_dip_gp + 88);

    auto tr_y_zzzz_z = pbuffer.data(idx_dip_gp + 89);

    auto tr_z_xxxx_x = pbuffer.data(idx_dip_gp + 90);

    auto tr_z_xxxx_y = pbuffer.data(idx_dip_gp + 91);

    auto tr_z_xxxx_z = pbuffer.data(idx_dip_gp + 92);

    auto tr_z_xxxz_x = pbuffer.data(idx_dip_gp + 96);

    auto tr_z_xxxz_z = pbuffer.data(idx_dip_gp + 98);

    auto tr_z_xxyy_y = pbuffer.data(idx_dip_gp + 100);

    auto tr_z_xxzz_x = pbuffer.data(idx_dip_gp + 105);

    auto tr_z_xxzz_y = pbuffer.data(idx_dip_gp + 106);

    auto tr_z_xxzz_z = pbuffer.data(idx_dip_gp + 107);

    auto tr_z_xyyy_y = pbuffer.data(idx_dip_gp + 109);

    auto tr_z_xzzz_x = pbuffer.data(idx_dip_gp + 117);

    auto tr_z_xzzz_y = pbuffer.data(idx_dip_gp + 118);

    auto tr_z_xzzz_z = pbuffer.data(idx_dip_gp + 119);

    auto tr_z_yyyy_x = pbuffer.data(idx_dip_gp + 120);

    auto tr_z_yyyy_y = pbuffer.data(idx_dip_gp + 121);

    auto tr_z_yyyy_z = pbuffer.data(idx_dip_gp + 122);

    auto tr_z_yyyz_x = pbuffer.data(idx_dip_gp + 123);

    auto tr_z_yyyz_y = pbuffer.data(idx_dip_gp + 124);

    auto tr_z_yyyz_z = pbuffer.data(idx_dip_gp + 125);

    auto tr_z_yyzz_x = pbuffer.data(idx_dip_gp + 126);

    auto tr_z_yyzz_y = pbuffer.data(idx_dip_gp + 127);

    auto tr_z_yyzz_z = pbuffer.data(idx_dip_gp + 128);

    auto tr_z_yzzz_x = pbuffer.data(idx_dip_gp + 129);

    auto tr_z_yzzz_y = pbuffer.data(idx_dip_gp + 130);

    auto tr_z_yzzz_z = pbuffer.data(idx_dip_gp + 131);

    auto tr_z_zzzz_x = pbuffer.data(idx_dip_gp + 132);

    auto tr_z_zzzz_y = pbuffer.data(idx_dip_gp + 133);

    auto tr_z_zzzz_z = pbuffer.data(idx_dip_gp + 134);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

    auto ts_xxxz_xz = pbuffer.data(idx_ovl_gd + 14);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

    auto ts_yyyz_yz = pbuffer.data(idx_ovl_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_ovl_gd + 71);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

    auto ts_yzzz_xz = pbuffer.data(idx_ovl_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_ovl_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_ovl_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_ovl_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_ovl_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_ovl_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_ovl_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_ovl_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_ovl_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_ovl_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto tr_x_xxxx_xx = pbuffer.data(idx_dip_gd);

    auto tr_x_xxxx_xy = pbuffer.data(idx_dip_gd + 1);

    auto tr_x_xxxx_xz = pbuffer.data(idx_dip_gd + 2);

    auto tr_x_xxxx_yy = pbuffer.data(idx_dip_gd + 3);

    auto tr_x_xxxx_yz = pbuffer.data(idx_dip_gd + 4);

    auto tr_x_xxxx_zz = pbuffer.data(idx_dip_gd + 5);

    auto tr_x_xxxy_xx = pbuffer.data(idx_dip_gd + 6);

    auto tr_x_xxxy_xy = pbuffer.data(idx_dip_gd + 7);

    auto tr_x_xxxy_xz = pbuffer.data(idx_dip_gd + 8);

    auto tr_x_xxxy_yy = pbuffer.data(idx_dip_gd + 9);

    auto tr_x_xxxy_zz = pbuffer.data(idx_dip_gd + 11);

    auto tr_x_xxxz_xx = pbuffer.data(idx_dip_gd + 12);

    auto tr_x_xxxz_xy = pbuffer.data(idx_dip_gd + 13);

    auto tr_x_xxxz_xz = pbuffer.data(idx_dip_gd + 14);

    auto tr_x_xxxz_yy = pbuffer.data(idx_dip_gd + 15);

    auto tr_x_xxxz_yz = pbuffer.data(idx_dip_gd + 16);

    auto tr_x_xxxz_zz = pbuffer.data(idx_dip_gd + 17);

    auto tr_x_xxyy_xx = pbuffer.data(idx_dip_gd + 18);

    auto tr_x_xxyy_xy = pbuffer.data(idx_dip_gd + 19);

    auto tr_x_xxyy_xz = pbuffer.data(idx_dip_gd + 20);

    auto tr_x_xxyy_yy = pbuffer.data(idx_dip_gd + 21);

    auto tr_x_xxyy_yz = pbuffer.data(idx_dip_gd + 22);

    auto tr_x_xxyy_zz = pbuffer.data(idx_dip_gd + 23);

    auto tr_x_xxyz_xz = pbuffer.data(idx_dip_gd + 26);

    auto tr_x_xxyz_zz = pbuffer.data(idx_dip_gd + 29);

    auto tr_x_xxzz_xx = pbuffer.data(idx_dip_gd + 30);

    auto tr_x_xxzz_xy = pbuffer.data(idx_dip_gd + 31);

    auto tr_x_xxzz_xz = pbuffer.data(idx_dip_gd + 32);

    auto tr_x_xxzz_yy = pbuffer.data(idx_dip_gd + 33);

    auto tr_x_xxzz_yz = pbuffer.data(idx_dip_gd + 34);

    auto tr_x_xxzz_zz = pbuffer.data(idx_dip_gd + 35);

    auto tr_x_xyyy_xx = pbuffer.data(idx_dip_gd + 36);

    auto tr_x_xyyy_xy = pbuffer.data(idx_dip_gd + 37);

    auto tr_x_xyyy_xz = pbuffer.data(idx_dip_gd + 38);

    auto tr_x_xyyy_yy = pbuffer.data(idx_dip_gd + 39);

    auto tr_x_xyyy_yz = pbuffer.data(idx_dip_gd + 40);

    auto tr_x_xyyz_xy = pbuffer.data(idx_dip_gd + 43);

    auto tr_x_xyyz_xz = pbuffer.data(idx_dip_gd + 44);

    auto tr_x_xyzz_xx = pbuffer.data(idx_dip_gd + 48);

    auto tr_x_xyzz_xz = pbuffer.data(idx_dip_gd + 50);

    auto tr_x_xzzz_xx = pbuffer.data(idx_dip_gd + 54);

    auto tr_x_xzzz_xy = pbuffer.data(idx_dip_gd + 55);

    auto tr_x_xzzz_xz = pbuffer.data(idx_dip_gd + 56);

    auto tr_x_xzzz_yz = pbuffer.data(idx_dip_gd + 58);

    auto tr_x_xzzz_zz = pbuffer.data(idx_dip_gd + 59);

    auto tr_x_yyyy_xx = pbuffer.data(idx_dip_gd + 60);

    auto tr_x_yyyy_xy = pbuffer.data(idx_dip_gd + 61);

    auto tr_x_yyyy_xz = pbuffer.data(idx_dip_gd + 62);

    auto tr_x_yyyy_yy = pbuffer.data(idx_dip_gd + 63);

    auto tr_x_yyyy_yz = pbuffer.data(idx_dip_gd + 64);

    auto tr_x_yyyy_zz = pbuffer.data(idx_dip_gd + 65);

    auto tr_x_yyyz_xy = pbuffer.data(idx_dip_gd + 67);

    auto tr_x_yyyz_xz = pbuffer.data(idx_dip_gd + 68);

    auto tr_x_yyyz_yy = pbuffer.data(idx_dip_gd + 69);

    auto tr_x_yyyz_yz = pbuffer.data(idx_dip_gd + 70);

    auto tr_x_yyyz_zz = pbuffer.data(idx_dip_gd + 71);

    auto tr_x_yyzz_xx = pbuffer.data(idx_dip_gd + 72);

    auto tr_x_yyzz_xy = pbuffer.data(idx_dip_gd + 73);

    auto tr_x_yyzz_xz = pbuffer.data(idx_dip_gd + 74);

    auto tr_x_yyzz_yy = pbuffer.data(idx_dip_gd + 75);

    auto tr_x_yyzz_yz = pbuffer.data(idx_dip_gd + 76);

    auto tr_x_yyzz_zz = pbuffer.data(idx_dip_gd + 77);

    auto tr_x_yzzz_xx = pbuffer.data(idx_dip_gd + 78);

    auto tr_x_yzzz_xz = pbuffer.data(idx_dip_gd + 80);

    auto tr_x_yzzz_yy = pbuffer.data(idx_dip_gd + 81);

    auto tr_x_yzzz_yz = pbuffer.data(idx_dip_gd + 82);

    auto tr_x_yzzz_zz = pbuffer.data(idx_dip_gd + 83);

    auto tr_x_zzzz_xx = pbuffer.data(idx_dip_gd + 84);

    auto tr_x_zzzz_xy = pbuffer.data(idx_dip_gd + 85);

    auto tr_x_zzzz_xz = pbuffer.data(idx_dip_gd + 86);

    auto tr_x_zzzz_yy = pbuffer.data(idx_dip_gd + 87);

    auto tr_x_zzzz_yz = pbuffer.data(idx_dip_gd + 88);

    auto tr_x_zzzz_zz = pbuffer.data(idx_dip_gd + 89);

    auto tr_y_xxxx_xx = pbuffer.data(idx_dip_gd + 90);

    auto tr_y_xxxx_xy = pbuffer.data(idx_dip_gd + 91);

    auto tr_y_xxxx_xz = pbuffer.data(idx_dip_gd + 92);

    auto tr_y_xxxx_yy = pbuffer.data(idx_dip_gd + 93);

    auto tr_y_xxxx_yz = pbuffer.data(idx_dip_gd + 94);

    auto tr_y_xxxx_zz = pbuffer.data(idx_dip_gd + 95);

    auto tr_y_xxxy_xx = pbuffer.data(idx_dip_gd + 96);

    auto tr_y_xxxy_xy = pbuffer.data(idx_dip_gd + 97);

    auto tr_y_xxxy_yy = pbuffer.data(idx_dip_gd + 99);

    auto tr_y_xxxy_yz = pbuffer.data(idx_dip_gd + 100);

    auto tr_y_xxxy_zz = pbuffer.data(idx_dip_gd + 101);

    auto tr_y_xxxz_xx = pbuffer.data(idx_dip_gd + 102);

    auto tr_y_xxxz_xy = pbuffer.data(idx_dip_gd + 103);

    auto tr_y_xxxz_xz = pbuffer.data(idx_dip_gd + 104);

    auto tr_y_xxxz_yz = pbuffer.data(idx_dip_gd + 106);

    auto tr_y_xxxz_zz = pbuffer.data(idx_dip_gd + 107);

    auto tr_y_xxyy_xx = pbuffer.data(idx_dip_gd + 108);

    auto tr_y_xxyy_xy = pbuffer.data(idx_dip_gd + 109);

    auto tr_y_xxyy_xz = pbuffer.data(idx_dip_gd + 110);

    auto tr_y_xxyy_yy = pbuffer.data(idx_dip_gd + 111);

    auto tr_y_xxyy_yz = pbuffer.data(idx_dip_gd + 112);

    auto tr_y_xxyy_zz = pbuffer.data(idx_dip_gd + 113);

    auto tr_y_xxyz_xy = pbuffer.data(idx_dip_gd + 115);

    auto tr_y_xxyz_yz = pbuffer.data(idx_dip_gd + 118);

    auto tr_y_xxyz_zz = pbuffer.data(idx_dip_gd + 119);

    auto tr_y_xxzz_xx = pbuffer.data(idx_dip_gd + 120);

    auto tr_y_xxzz_xy = pbuffer.data(idx_dip_gd + 121);

    auto tr_y_xxzz_xz = pbuffer.data(idx_dip_gd + 122);

    auto tr_y_xxzz_yy = pbuffer.data(idx_dip_gd + 123);

    auto tr_y_xxzz_yz = pbuffer.data(idx_dip_gd + 124);

    auto tr_y_xxzz_zz = pbuffer.data(idx_dip_gd + 125);

    auto tr_y_xyyy_xx = pbuffer.data(idx_dip_gd + 126);

    auto tr_y_xyyy_xy = pbuffer.data(idx_dip_gd + 127);

    auto tr_y_xyyy_xz = pbuffer.data(idx_dip_gd + 128);

    auto tr_y_xyyy_yy = pbuffer.data(idx_dip_gd + 129);

    auto tr_y_xyyy_yz = pbuffer.data(idx_dip_gd + 130);

    auto tr_y_xyyy_zz = pbuffer.data(idx_dip_gd + 131);

    auto tr_y_xyyz_yz = pbuffer.data(idx_dip_gd + 136);

    auto tr_y_xyyz_zz = pbuffer.data(idx_dip_gd + 137);

    auto tr_y_xyzz_yy = pbuffer.data(idx_dip_gd + 141);

    auto tr_y_xyzz_yz = pbuffer.data(idx_dip_gd + 142);

    auto tr_y_xyzz_zz = pbuffer.data(idx_dip_gd + 143);

    auto tr_y_xzzz_xz = pbuffer.data(idx_dip_gd + 146);

    auto tr_y_xzzz_yy = pbuffer.data(idx_dip_gd + 147);

    auto tr_y_xzzz_yz = pbuffer.data(idx_dip_gd + 148);

    auto tr_y_xzzz_zz = pbuffer.data(idx_dip_gd + 149);

    auto tr_y_yyyy_xx = pbuffer.data(idx_dip_gd + 150);

    auto tr_y_yyyy_xy = pbuffer.data(idx_dip_gd + 151);

    auto tr_y_yyyy_xz = pbuffer.data(idx_dip_gd + 152);

    auto tr_y_yyyy_yy = pbuffer.data(idx_dip_gd + 153);

    auto tr_y_yyyy_yz = pbuffer.data(idx_dip_gd + 154);

    auto tr_y_yyyy_zz = pbuffer.data(idx_dip_gd + 155);

    auto tr_y_yyyz_xx = pbuffer.data(idx_dip_gd + 156);

    auto tr_y_yyyz_xy = pbuffer.data(idx_dip_gd + 157);

    auto tr_y_yyyz_xz = pbuffer.data(idx_dip_gd + 158);

    auto tr_y_yyyz_yy = pbuffer.data(idx_dip_gd + 159);

    auto tr_y_yyyz_yz = pbuffer.data(idx_dip_gd + 160);

    auto tr_y_yyyz_zz = pbuffer.data(idx_dip_gd + 161);

    auto tr_y_yyzz_xx = pbuffer.data(idx_dip_gd + 162);

    auto tr_y_yyzz_xy = pbuffer.data(idx_dip_gd + 163);

    auto tr_y_yyzz_xz = pbuffer.data(idx_dip_gd + 164);

    auto tr_y_yyzz_yy = pbuffer.data(idx_dip_gd + 165);

    auto tr_y_yyzz_yz = pbuffer.data(idx_dip_gd + 166);

    auto tr_y_yyzz_zz = pbuffer.data(idx_dip_gd + 167);

    auto tr_y_yzzz_xx = pbuffer.data(idx_dip_gd + 168);

    auto tr_y_yzzz_xy = pbuffer.data(idx_dip_gd + 169);

    auto tr_y_yzzz_xz = pbuffer.data(idx_dip_gd + 170);

    auto tr_y_yzzz_yy = pbuffer.data(idx_dip_gd + 171);

    auto tr_y_yzzz_yz = pbuffer.data(idx_dip_gd + 172);

    auto tr_y_yzzz_zz = pbuffer.data(idx_dip_gd + 173);

    auto tr_y_zzzz_xx = pbuffer.data(idx_dip_gd + 174);

    auto tr_y_zzzz_xy = pbuffer.data(idx_dip_gd + 175);

    auto tr_y_zzzz_xz = pbuffer.data(idx_dip_gd + 176);

    auto tr_y_zzzz_yy = pbuffer.data(idx_dip_gd + 177);

    auto tr_y_zzzz_yz = pbuffer.data(idx_dip_gd + 178);

    auto tr_y_zzzz_zz = pbuffer.data(idx_dip_gd + 179);

    auto tr_z_xxxx_xx = pbuffer.data(idx_dip_gd + 180);

    auto tr_z_xxxx_xy = pbuffer.data(idx_dip_gd + 181);

    auto tr_z_xxxx_xz = pbuffer.data(idx_dip_gd + 182);

    auto tr_z_xxxx_yy = pbuffer.data(idx_dip_gd + 183);

    auto tr_z_xxxx_yz = pbuffer.data(idx_dip_gd + 184);

    auto tr_z_xxxx_zz = pbuffer.data(idx_dip_gd + 185);

    auto tr_z_xxxy_xx = pbuffer.data(idx_dip_gd + 186);

    auto tr_z_xxxy_xz = pbuffer.data(idx_dip_gd + 188);

    auto tr_z_xxxy_yy = pbuffer.data(idx_dip_gd + 189);

    auto tr_z_xxxy_yz = pbuffer.data(idx_dip_gd + 190);

    auto tr_z_xxxz_xx = pbuffer.data(idx_dip_gd + 192);

    auto tr_z_xxxz_xy = pbuffer.data(idx_dip_gd + 193);

    auto tr_z_xxxz_xz = pbuffer.data(idx_dip_gd + 194);

    auto tr_z_xxxz_yy = pbuffer.data(idx_dip_gd + 195);

    auto tr_z_xxxz_yz = pbuffer.data(idx_dip_gd + 196);

    auto tr_z_xxxz_zz = pbuffer.data(idx_dip_gd + 197);

    auto tr_z_xxyy_xx = pbuffer.data(idx_dip_gd + 198);

    auto tr_z_xxyy_xy = pbuffer.data(idx_dip_gd + 199);

    auto tr_z_xxyy_xz = pbuffer.data(idx_dip_gd + 200);

    auto tr_z_xxyy_yy = pbuffer.data(idx_dip_gd + 201);

    auto tr_z_xxyy_yz = pbuffer.data(idx_dip_gd + 202);

    auto tr_z_xxyy_zz = pbuffer.data(idx_dip_gd + 203);

    auto tr_z_xxyz_xx = pbuffer.data(idx_dip_gd + 204);

    auto tr_z_xxyz_xz = pbuffer.data(idx_dip_gd + 206);

    auto tr_z_xxyz_yy = pbuffer.data(idx_dip_gd + 207);

    auto tr_z_xxyz_yz = pbuffer.data(idx_dip_gd + 208);

    auto tr_z_xxzz_xx = pbuffer.data(idx_dip_gd + 210);

    auto tr_z_xxzz_xy = pbuffer.data(idx_dip_gd + 211);

    auto tr_z_xxzz_xz = pbuffer.data(idx_dip_gd + 212);

    auto tr_z_xxzz_yy = pbuffer.data(idx_dip_gd + 213);

    auto tr_z_xxzz_yz = pbuffer.data(idx_dip_gd + 214);

    auto tr_z_xxzz_zz = pbuffer.data(idx_dip_gd + 215);

    auto tr_z_xyyy_xy = pbuffer.data(idx_dip_gd + 217);

    auto tr_z_xyyy_yy = pbuffer.data(idx_dip_gd + 219);

    auto tr_z_xyyy_yz = pbuffer.data(idx_dip_gd + 220);

    auto tr_z_xyyy_zz = pbuffer.data(idx_dip_gd + 221);

    auto tr_z_xyyz_yy = pbuffer.data(idx_dip_gd + 225);

    auto tr_z_xyyz_yz = pbuffer.data(idx_dip_gd + 226);

    auto tr_z_xyyz_zz = pbuffer.data(idx_dip_gd + 227);

    auto tr_z_xyzz_yy = pbuffer.data(idx_dip_gd + 231);

    auto tr_z_xyzz_yz = pbuffer.data(idx_dip_gd + 232);

    auto tr_z_xzzz_xx = pbuffer.data(idx_dip_gd + 234);

    auto tr_z_xzzz_xy = pbuffer.data(idx_dip_gd + 235);

    auto tr_z_xzzz_xz = pbuffer.data(idx_dip_gd + 236);

    auto tr_z_xzzz_yy = pbuffer.data(idx_dip_gd + 237);

    auto tr_z_xzzz_yz = pbuffer.data(idx_dip_gd + 238);

    auto tr_z_xzzz_zz = pbuffer.data(idx_dip_gd + 239);

    auto tr_z_yyyy_xx = pbuffer.data(idx_dip_gd + 240);

    auto tr_z_yyyy_xy = pbuffer.data(idx_dip_gd + 241);

    auto tr_z_yyyy_xz = pbuffer.data(idx_dip_gd + 242);

    auto tr_z_yyyy_yy = pbuffer.data(idx_dip_gd + 243);

    auto tr_z_yyyy_yz = pbuffer.data(idx_dip_gd + 244);

    auto tr_z_yyyy_zz = pbuffer.data(idx_dip_gd + 245);

    auto tr_z_yyyz_xx = pbuffer.data(idx_dip_gd + 246);

    auto tr_z_yyyz_xy = pbuffer.data(idx_dip_gd + 247);

    auto tr_z_yyyz_xz = pbuffer.data(idx_dip_gd + 248);

    auto tr_z_yyyz_yy = pbuffer.data(idx_dip_gd + 249);

    auto tr_z_yyyz_yz = pbuffer.data(idx_dip_gd + 250);

    auto tr_z_yyyz_zz = pbuffer.data(idx_dip_gd + 251);

    auto tr_z_yyzz_xx = pbuffer.data(idx_dip_gd + 252);

    auto tr_z_yyzz_xy = pbuffer.data(idx_dip_gd + 253);

    auto tr_z_yyzz_xz = pbuffer.data(idx_dip_gd + 254);

    auto tr_z_yyzz_yy = pbuffer.data(idx_dip_gd + 255);

    auto tr_z_yyzz_yz = pbuffer.data(idx_dip_gd + 256);

    auto tr_z_yyzz_zz = pbuffer.data(idx_dip_gd + 257);

    auto tr_z_yzzz_xx = pbuffer.data(idx_dip_gd + 258);

    auto tr_z_yzzz_xy = pbuffer.data(idx_dip_gd + 259);

    auto tr_z_yzzz_xz = pbuffer.data(idx_dip_gd + 260);

    auto tr_z_yzzz_yy = pbuffer.data(idx_dip_gd + 261);

    auto tr_z_yzzz_yz = pbuffer.data(idx_dip_gd + 262);

    auto tr_z_yzzz_zz = pbuffer.data(idx_dip_gd + 263);

    auto tr_z_zzzz_xx = pbuffer.data(idx_dip_gd + 264);

    auto tr_z_zzzz_xy = pbuffer.data(idx_dip_gd + 265);

    auto tr_z_zzzz_xz = pbuffer.data(idx_dip_gd + 266);

    auto tr_z_zzzz_yy = pbuffer.data(idx_dip_gd + 267);

    auto tr_z_zzzz_yz = pbuffer.data(idx_dip_gd + 268);

    auto tr_z_zzzz_zz = pbuffer.data(idx_dip_gd + 269);

    // Set up 0-6 components of targeted buffer : HD

    auto tr_x_xxxxx_xx = pbuffer.data(idx_dip_hd);

    auto tr_x_xxxxx_xy = pbuffer.data(idx_dip_hd + 1);

    auto tr_x_xxxxx_xz = pbuffer.data(idx_dip_hd + 2);

    auto tr_x_xxxxx_yy = pbuffer.data(idx_dip_hd + 3);

    auto tr_x_xxxxx_yz = pbuffer.data(idx_dip_hd + 4);

    auto tr_x_xxxxx_zz = pbuffer.data(idx_dip_hd + 5);

    #pragma omp simd aligned(pa_x, tr_x_xxx_xx, tr_x_xxx_xy, tr_x_xxx_xz, tr_x_xxx_yy, tr_x_xxx_yz, tr_x_xxx_zz, tr_x_xxxx_x, tr_x_xxxx_xx, tr_x_xxxx_xy, tr_x_xxxx_xz, tr_x_xxxx_y, tr_x_xxxx_yy, tr_x_xxxx_yz, tr_x_xxxx_z, tr_x_xxxx_zz, tr_x_xxxxx_xx, tr_x_xxxxx_xy, tr_x_xxxxx_xz, tr_x_xxxxx_yy, tr_x_xxxxx_yz, tr_x_xxxxx_zz, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_xx[i] = 4.0 * tr_x_xxx_xx[i] * fe_0 + 2.0 * tr_x_xxxx_x[i] * fe_0 + ts_xxxx_xx[i] * fe_0 + tr_x_xxxx_xx[i] * pa_x[i];

        tr_x_xxxxx_xy[i] = 4.0 * tr_x_xxx_xy[i] * fe_0 + tr_x_xxxx_y[i] * fe_0 + ts_xxxx_xy[i] * fe_0 + tr_x_xxxx_xy[i] * pa_x[i];

        tr_x_xxxxx_xz[i] = 4.0 * tr_x_xxx_xz[i] * fe_0 + tr_x_xxxx_z[i] * fe_0 + ts_xxxx_xz[i] * fe_0 + tr_x_xxxx_xz[i] * pa_x[i];

        tr_x_xxxxx_yy[i] = 4.0 * tr_x_xxx_yy[i] * fe_0 + ts_xxxx_yy[i] * fe_0 + tr_x_xxxx_yy[i] * pa_x[i];

        tr_x_xxxxx_yz[i] = 4.0 * tr_x_xxx_yz[i] * fe_0 + ts_xxxx_yz[i] * fe_0 + tr_x_xxxx_yz[i] * pa_x[i];

        tr_x_xxxxx_zz[i] = 4.0 * tr_x_xxx_zz[i] * fe_0 + ts_xxxx_zz[i] * fe_0 + tr_x_xxxx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : HD

    auto tr_x_xxxxy_xx = pbuffer.data(idx_dip_hd + 6);

    auto tr_x_xxxxy_xy = pbuffer.data(idx_dip_hd + 7);

    auto tr_x_xxxxy_xz = pbuffer.data(idx_dip_hd + 8);

    auto tr_x_xxxxy_yy = pbuffer.data(idx_dip_hd + 9);

    auto tr_x_xxxxy_yz = pbuffer.data(idx_dip_hd + 10);

    auto tr_x_xxxxy_zz = pbuffer.data(idx_dip_hd + 11);

    #pragma omp simd aligned(pa_y, tr_x_xxxx_x, tr_x_xxxx_xx, tr_x_xxxx_xy, tr_x_xxxx_xz, tr_x_xxxx_y, tr_x_xxxx_yy, tr_x_xxxx_yz, tr_x_xxxx_z, tr_x_xxxx_zz, tr_x_xxxxy_xx, tr_x_xxxxy_xy, tr_x_xxxxy_xz, tr_x_xxxxy_yy, tr_x_xxxxy_yz, tr_x_xxxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxy_xx[i] = tr_x_xxxx_xx[i] * pa_y[i];

        tr_x_xxxxy_xy[i] = tr_x_xxxx_x[i] * fe_0 + tr_x_xxxx_xy[i] * pa_y[i];

        tr_x_xxxxy_xz[i] = tr_x_xxxx_xz[i] * pa_y[i];

        tr_x_xxxxy_yy[i] = 2.0 * tr_x_xxxx_y[i] * fe_0 + tr_x_xxxx_yy[i] * pa_y[i];

        tr_x_xxxxy_yz[i] = tr_x_xxxx_z[i] * fe_0 + tr_x_xxxx_yz[i] * pa_y[i];

        tr_x_xxxxy_zz[i] = tr_x_xxxx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : HD

    auto tr_x_xxxxz_xx = pbuffer.data(idx_dip_hd + 12);

    auto tr_x_xxxxz_xy = pbuffer.data(idx_dip_hd + 13);

    auto tr_x_xxxxz_xz = pbuffer.data(idx_dip_hd + 14);

    auto tr_x_xxxxz_yy = pbuffer.data(idx_dip_hd + 15);

    auto tr_x_xxxxz_yz = pbuffer.data(idx_dip_hd + 16);

    auto tr_x_xxxxz_zz = pbuffer.data(idx_dip_hd + 17);

    #pragma omp simd aligned(pa_z, tr_x_xxxx_x, tr_x_xxxx_xx, tr_x_xxxx_xy, tr_x_xxxx_xz, tr_x_xxxx_y, tr_x_xxxx_yy, tr_x_xxxx_yz, tr_x_xxxx_z, tr_x_xxxx_zz, tr_x_xxxxz_xx, tr_x_xxxxz_xy, tr_x_xxxxz_xz, tr_x_xxxxz_yy, tr_x_xxxxz_yz, tr_x_xxxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxz_xx[i] = tr_x_xxxx_xx[i] * pa_z[i];

        tr_x_xxxxz_xy[i] = tr_x_xxxx_xy[i] * pa_z[i];

        tr_x_xxxxz_xz[i] = tr_x_xxxx_x[i] * fe_0 + tr_x_xxxx_xz[i] * pa_z[i];

        tr_x_xxxxz_yy[i] = tr_x_xxxx_yy[i] * pa_z[i];

        tr_x_xxxxz_yz[i] = tr_x_xxxx_y[i] * fe_0 + tr_x_xxxx_yz[i] * pa_z[i];

        tr_x_xxxxz_zz[i] = 2.0 * tr_x_xxxx_z[i] * fe_0 + tr_x_xxxx_zz[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : HD

    auto tr_x_xxxyy_xx = pbuffer.data(idx_dip_hd + 18);

    auto tr_x_xxxyy_xy = pbuffer.data(idx_dip_hd + 19);

    auto tr_x_xxxyy_xz = pbuffer.data(idx_dip_hd + 20);

    auto tr_x_xxxyy_yy = pbuffer.data(idx_dip_hd + 21);

    auto tr_x_xxxyy_yz = pbuffer.data(idx_dip_hd + 22);

    auto tr_x_xxxyy_zz = pbuffer.data(idx_dip_hd + 23);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxx_xx, tr_x_xxx_xy, tr_x_xxx_xz, tr_x_xxx_zz, tr_x_xxxy_x, tr_x_xxxy_xx, tr_x_xxxy_xy, tr_x_xxxy_xz, tr_x_xxxy_zz, tr_x_xxxyy_xx, tr_x_xxxyy_xy, tr_x_xxxyy_xz, tr_x_xxxyy_yy, tr_x_xxxyy_yz, tr_x_xxxyy_zz, tr_x_xxyy_yy, tr_x_xxyy_yz, tr_x_xyy_yy, tr_x_xyy_yz, ts_xxyy_yy, ts_xxyy_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyy_xx[i] = tr_x_xxx_xx[i] * fe_0 + tr_x_xxxy_xx[i] * pa_y[i];

        tr_x_xxxyy_xy[i] = tr_x_xxx_xy[i] * fe_0 + tr_x_xxxy_x[i] * fe_0 + tr_x_xxxy_xy[i] * pa_y[i];

        tr_x_xxxyy_xz[i] = tr_x_xxx_xz[i] * fe_0 + tr_x_xxxy_xz[i] * pa_y[i];

        tr_x_xxxyy_yy[i] = 2.0 * tr_x_xyy_yy[i] * fe_0 + ts_xxyy_yy[i] * fe_0 + tr_x_xxyy_yy[i] * pa_x[i];

        tr_x_xxxyy_yz[i] = 2.0 * tr_x_xyy_yz[i] * fe_0 + ts_xxyy_yz[i] * fe_0 + tr_x_xxyy_yz[i] * pa_x[i];

        tr_x_xxxyy_zz[i] = tr_x_xxx_zz[i] * fe_0 + tr_x_xxxy_zz[i] * pa_y[i];
    }

    // Set up 24-30 components of targeted buffer : HD

    auto tr_x_xxxyz_xx = pbuffer.data(idx_dip_hd + 24);

    auto tr_x_xxxyz_xy = pbuffer.data(idx_dip_hd + 25);

    auto tr_x_xxxyz_xz = pbuffer.data(idx_dip_hd + 26);

    auto tr_x_xxxyz_yy = pbuffer.data(idx_dip_hd + 27);

    auto tr_x_xxxyz_yz = pbuffer.data(idx_dip_hd + 28);

    auto tr_x_xxxyz_zz = pbuffer.data(idx_dip_hd + 29);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxxy_xy, tr_x_xxxy_yy, tr_x_xxxyz_xx, tr_x_xxxyz_xy, tr_x_xxxyz_xz, tr_x_xxxyz_yy, tr_x_xxxyz_yz, tr_x_xxxyz_zz, tr_x_xxxz_xx, tr_x_xxxz_xz, tr_x_xxxz_yz, tr_x_xxxz_z, tr_x_xxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyz_xx[i] = tr_x_xxxz_xx[i] * pa_y[i];

        tr_x_xxxyz_xy[i] = tr_x_xxxy_xy[i] * pa_z[i];

        tr_x_xxxyz_xz[i] = tr_x_xxxz_xz[i] * pa_y[i];

        tr_x_xxxyz_yy[i] = tr_x_xxxy_yy[i] * pa_z[i];

        tr_x_xxxyz_yz[i] = tr_x_xxxz_z[i] * fe_0 + tr_x_xxxz_yz[i] * pa_y[i];

        tr_x_xxxyz_zz[i] = tr_x_xxxz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : HD

    auto tr_x_xxxzz_xx = pbuffer.data(idx_dip_hd + 30);

    auto tr_x_xxxzz_xy = pbuffer.data(idx_dip_hd + 31);

    auto tr_x_xxxzz_xz = pbuffer.data(idx_dip_hd + 32);

    auto tr_x_xxxzz_yy = pbuffer.data(idx_dip_hd + 33);

    auto tr_x_xxxzz_yz = pbuffer.data(idx_dip_hd + 34);

    auto tr_x_xxxzz_zz = pbuffer.data(idx_dip_hd + 35);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxx_xx, tr_x_xxx_xy, tr_x_xxx_xz, tr_x_xxx_yy, tr_x_xxxz_x, tr_x_xxxz_xx, tr_x_xxxz_xy, tr_x_xxxz_xz, tr_x_xxxz_yy, tr_x_xxxzz_xx, tr_x_xxxzz_xy, tr_x_xxxzz_xz, tr_x_xxxzz_yy, tr_x_xxxzz_yz, tr_x_xxxzz_zz, tr_x_xxzz_yz, tr_x_xxzz_zz, tr_x_xzz_yz, tr_x_xzz_zz, ts_xxzz_yz, ts_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzz_xx[i] = tr_x_xxx_xx[i] * fe_0 + tr_x_xxxz_xx[i] * pa_z[i];

        tr_x_xxxzz_xy[i] = tr_x_xxx_xy[i] * fe_0 + tr_x_xxxz_xy[i] * pa_z[i];

        tr_x_xxxzz_xz[i] = tr_x_xxx_xz[i] * fe_0 + tr_x_xxxz_x[i] * fe_0 + tr_x_xxxz_xz[i] * pa_z[i];

        tr_x_xxxzz_yy[i] = tr_x_xxx_yy[i] * fe_0 + tr_x_xxxz_yy[i] * pa_z[i];

        tr_x_xxxzz_yz[i] = 2.0 * tr_x_xzz_yz[i] * fe_0 + ts_xxzz_yz[i] * fe_0 + tr_x_xxzz_yz[i] * pa_x[i];

        tr_x_xxxzz_zz[i] = 2.0 * tr_x_xzz_zz[i] * fe_0 + ts_xxzz_zz[i] * fe_0 + tr_x_xxzz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : HD

    auto tr_x_xxyyy_xx = pbuffer.data(idx_dip_hd + 36);

    auto tr_x_xxyyy_xy = pbuffer.data(idx_dip_hd + 37);

    auto tr_x_xxyyy_xz = pbuffer.data(idx_dip_hd + 38);

    auto tr_x_xxyyy_yy = pbuffer.data(idx_dip_hd + 39);

    auto tr_x_xxyyy_yz = pbuffer.data(idx_dip_hd + 40);

    auto tr_x_xxyyy_zz = pbuffer.data(idx_dip_hd + 41);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xxy_xx, tr_x_xxy_xy, tr_x_xxy_xz, tr_x_xxy_zz, tr_x_xxyy_x, tr_x_xxyy_xx, tr_x_xxyy_xy, tr_x_xxyy_xz, tr_x_xxyy_zz, tr_x_xxyyy_xx, tr_x_xxyyy_xy, tr_x_xxyyy_xz, tr_x_xxyyy_yy, tr_x_xxyyy_yz, tr_x_xxyyy_zz, tr_x_xyyy_yy, tr_x_xyyy_yz, tr_x_yyy_yy, tr_x_yyy_yz, ts_xyyy_yy, ts_xyyy_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyy_xx[i] = 2.0 * tr_x_xxy_xx[i] * fe_0 + tr_x_xxyy_xx[i] * pa_y[i];

        tr_x_xxyyy_xy[i] = 2.0 * tr_x_xxy_xy[i] * fe_0 + tr_x_xxyy_x[i] * fe_0 + tr_x_xxyy_xy[i] * pa_y[i];

        tr_x_xxyyy_xz[i] = 2.0 * tr_x_xxy_xz[i] * fe_0 + tr_x_xxyy_xz[i] * pa_y[i];

        tr_x_xxyyy_yy[i] = tr_x_yyy_yy[i] * fe_0 + ts_xyyy_yy[i] * fe_0 + tr_x_xyyy_yy[i] * pa_x[i];

        tr_x_xxyyy_yz[i] = tr_x_yyy_yz[i] * fe_0 + ts_xyyy_yz[i] * fe_0 + tr_x_xyyy_yz[i] * pa_x[i];

        tr_x_xxyyy_zz[i] = 2.0 * tr_x_xxy_zz[i] * fe_0 + tr_x_xxyy_zz[i] * pa_y[i];
    }

    // Set up 42-48 components of targeted buffer : HD

    auto tr_x_xxyyz_xx = pbuffer.data(idx_dip_hd + 42);

    auto tr_x_xxyyz_xy = pbuffer.data(idx_dip_hd + 43);

    auto tr_x_xxyyz_xz = pbuffer.data(idx_dip_hd + 44);

    auto tr_x_xxyyz_yy = pbuffer.data(idx_dip_hd + 45);

    auto tr_x_xxyyz_yz = pbuffer.data(idx_dip_hd + 46);

    auto tr_x_xxyyz_zz = pbuffer.data(idx_dip_hd + 47);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_xxyy_xx, tr_x_xxyy_xy, tr_x_xxyy_y, tr_x_xxyy_yy, tr_x_xxyy_yz, tr_x_xxyyz_xx, tr_x_xxyyz_xy, tr_x_xxyyz_xz, tr_x_xxyyz_yy, tr_x_xxyyz_yz, tr_x_xxyyz_zz, tr_x_xxyz_xz, tr_x_xxyz_zz, tr_x_xxz_xz, tr_x_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyz_xx[i] = tr_x_xxyy_xx[i] * pa_z[i];

        tr_x_xxyyz_xy[i] = tr_x_xxyy_xy[i] * pa_z[i];

        tr_x_xxyyz_xz[i] = tr_x_xxz_xz[i] * fe_0 + tr_x_xxyz_xz[i] * pa_y[i];

        tr_x_xxyyz_yy[i] = tr_x_xxyy_yy[i] * pa_z[i];

        tr_x_xxyyz_yz[i] = tr_x_xxyy_y[i] * fe_0 + tr_x_xxyy_yz[i] * pa_z[i];

        tr_x_xxyyz_zz[i] = tr_x_xxz_zz[i] * fe_0 + tr_x_xxyz_zz[i] * pa_y[i];
    }

    // Set up 48-54 components of targeted buffer : HD

    auto tr_x_xxyzz_xx = pbuffer.data(idx_dip_hd + 48);

    auto tr_x_xxyzz_xy = pbuffer.data(idx_dip_hd + 49);

    auto tr_x_xxyzz_xz = pbuffer.data(idx_dip_hd + 50);

    auto tr_x_xxyzz_yy = pbuffer.data(idx_dip_hd + 51);

    auto tr_x_xxyzz_yz = pbuffer.data(idx_dip_hd + 52);

    auto tr_x_xxyzz_zz = pbuffer.data(idx_dip_hd + 53);

    #pragma omp simd aligned(pa_y, tr_x_xxyzz_xx, tr_x_xxyzz_xy, tr_x_xxyzz_xz, tr_x_xxyzz_yy, tr_x_xxyzz_yz, tr_x_xxyzz_zz, tr_x_xxzz_x, tr_x_xxzz_xx, tr_x_xxzz_xy, tr_x_xxzz_xz, tr_x_xxzz_y, tr_x_xxzz_yy, tr_x_xxzz_yz, tr_x_xxzz_z, tr_x_xxzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzz_xx[i] = tr_x_xxzz_xx[i] * pa_y[i];

        tr_x_xxyzz_xy[i] = tr_x_xxzz_x[i] * fe_0 + tr_x_xxzz_xy[i] * pa_y[i];

        tr_x_xxyzz_xz[i] = tr_x_xxzz_xz[i] * pa_y[i];

        tr_x_xxyzz_yy[i] = 2.0 * tr_x_xxzz_y[i] * fe_0 + tr_x_xxzz_yy[i] * pa_y[i];

        tr_x_xxyzz_yz[i] = tr_x_xxzz_z[i] * fe_0 + tr_x_xxzz_yz[i] * pa_y[i];

        tr_x_xxyzz_zz[i] = tr_x_xxzz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : HD

    auto tr_x_xxzzz_xx = pbuffer.data(idx_dip_hd + 54);

    auto tr_x_xxzzz_xy = pbuffer.data(idx_dip_hd + 55);

    auto tr_x_xxzzz_xz = pbuffer.data(idx_dip_hd + 56);

    auto tr_x_xxzzz_yy = pbuffer.data(idx_dip_hd + 57);

    auto tr_x_xxzzz_yz = pbuffer.data(idx_dip_hd + 58);

    auto tr_x_xxzzz_zz = pbuffer.data(idx_dip_hd + 59);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xxz_xx, tr_x_xxz_xy, tr_x_xxz_xz, tr_x_xxz_yy, tr_x_xxzz_x, tr_x_xxzz_xx, tr_x_xxzz_xy, tr_x_xxzz_xz, tr_x_xxzz_yy, tr_x_xxzzz_xx, tr_x_xxzzz_xy, tr_x_xxzzz_xz, tr_x_xxzzz_yy, tr_x_xxzzz_yz, tr_x_xxzzz_zz, tr_x_xzzz_yz, tr_x_xzzz_zz, tr_x_zzz_yz, tr_x_zzz_zz, ts_xzzz_yz, ts_xzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzz_xx[i] = 2.0 * tr_x_xxz_xx[i] * fe_0 + tr_x_xxzz_xx[i] * pa_z[i];

        tr_x_xxzzz_xy[i] = 2.0 * tr_x_xxz_xy[i] * fe_0 + tr_x_xxzz_xy[i] * pa_z[i];

        tr_x_xxzzz_xz[i] = 2.0 * tr_x_xxz_xz[i] * fe_0 + tr_x_xxzz_x[i] * fe_0 + tr_x_xxzz_xz[i] * pa_z[i];

        tr_x_xxzzz_yy[i] = 2.0 * tr_x_xxz_yy[i] * fe_0 + tr_x_xxzz_yy[i] * pa_z[i];

        tr_x_xxzzz_yz[i] = tr_x_zzz_yz[i] * fe_0 + ts_xzzz_yz[i] * fe_0 + tr_x_xzzz_yz[i] * pa_x[i];

        tr_x_xxzzz_zz[i] = tr_x_zzz_zz[i] * fe_0 + ts_xzzz_zz[i] * fe_0 + tr_x_xzzz_zz[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : HD

    auto tr_x_xyyyy_xx = pbuffer.data(idx_dip_hd + 60);

    auto tr_x_xyyyy_xy = pbuffer.data(idx_dip_hd + 61);

    auto tr_x_xyyyy_xz = pbuffer.data(idx_dip_hd + 62);

    auto tr_x_xyyyy_yy = pbuffer.data(idx_dip_hd + 63);

    auto tr_x_xyyyy_yz = pbuffer.data(idx_dip_hd + 64);

    auto tr_x_xyyyy_zz = pbuffer.data(idx_dip_hd + 65);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyy_xx, tr_x_xyy_xz, tr_x_xyyy_xx, tr_x_xyyy_xz, tr_x_xyyyy_xx, tr_x_xyyyy_xy, tr_x_xyyyy_xz, tr_x_xyyyy_yy, tr_x_xyyyy_yz, tr_x_xyyyy_zz, tr_x_yyyy_xy, tr_x_yyyy_y, tr_x_yyyy_yy, tr_x_yyyy_yz, tr_x_yyyy_zz, ts_yyyy_xy, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyy_xx[i] = 3.0 * tr_x_xyy_xx[i] * fe_0 + tr_x_xyyy_xx[i] * pa_y[i];

        tr_x_xyyyy_xy[i] = tr_x_yyyy_y[i] * fe_0 + ts_yyyy_xy[i] * fe_0 + tr_x_yyyy_xy[i] * pa_x[i];

        tr_x_xyyyy_xz[i] = 3.0 * tr_x_xyy_xz[i] * fe_0 + tr_x_xyyy_xz[i] * pa_y[i];

        tr_x_xyyyy_yy[i] = ts_yyyy_yy[i] * fe_0 + tr_x_yyyy_yy[i] * pa_x[i];

        tr_x_xyyyy_yz[i] = ts_yyyy_yz[i] * fe_0 + tr_x_yyyy_yz[i] * pa_x[i];

        tr_x_xyyyy_zz[i] = ts_yyyy_zz[i] * fe_0 + tr_x_yyyy_zz[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : HD

    auto tr_x_xyyyz_xx = pbuffer.data(idx_dip_hd + 66);

    auto tr_x_xyyyz_xy = pbuffer.data(idx_dip_hd + 67);

    auto tr_x_xyyyz_xz = pbuffer.data(idx_dip_hd + 68);

    auto tr_x_xyyyz_yy = pbuffer.data(idx_dip_hd + 69);

    auto tr_x_xyyyz_yz = pbuffer.data(idx_dip_hd + 70);

    auto tr_x_xyyyz_zz = pbuffer.data(idx_dip_hd + 71);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyyy_xx, tr_x_xyyy_xy, tr_x_xyyy_yy, tr_x_xyyyz_xx, tr_x_xyyyz_xy, tr_x_xyyyz_xz, tr_x_xyyyz_yy, tr_x_xyyyz_yz, tr_x_xyyyz_zz, tr_x_xyyz_xz, tr_x_xyz_xz, tr_x_yyyz_yz, tr_x_yyyz_zz, ts_yyyz_yz, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyz_xx[i] = tr_x_xyyy_xx[i] * pa_z[i];

        tr_x_xyyyz_xy[i] = tr_x_xyyy_xy[i] * pa_z[i];

        tr_x_xyyyz_xz[i] = 2.0 * tr_x_xyz_xz[i] * fe_0 + tr_x_xyyz_xz[i] * pa_y[i];

        tr_x_xyyyz_yy[i] = tr_x_xyyy_yy[i] * pa_z[i];

        tr_x_xyyyz_yz[i] = ts_yyyz_yz[i] * fe_0 + tr_x_yyyz_yz[i] * pa_x[i];

        tr_x_xyyyz_zz[i] = ts_yyyz_zz[i] * fe_0 + tr_x_yyyz_zz[i] * pa_x[i];
    }

    // Set up 72-78 components of targeted buffer : HD

    auto tr_x_xyyzz_xx = pbuffer.data(idx_dip_hd + 72);

    auto tr_x_xyyzz_xy = pbuffer.data(idx_dip_hd + 73);

    auto tr_x_xyyzz_xz = pbuffer.data(idx_dip_hd + 74);

    auto tr_x_xyyzz_yy = pbuffer.data(idx_dip_hd + 75);

    auto tr_x_xyyzz_yz = pbuffer.data(idx_dip_hd + 76);

    auto tr_x_xyyzz_zz = pbuffer.data(idx_dip_hd + 77);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xyy_xy, tr_x_xyyz_xy, tr_x_xyyzz_xx, tr_x_xyyzz_xy, tr_x_xyyzz_xz, tr_x_xyyzz_yy, tr_x_xyyzz_yz, tr_x_xyyzz_zz, tr_x_xyzz_xx, tr_x_xyzz_xz, tr_x_xzz_xx, tr_x_xzz_xz, tr_x_yyzz_yy, tr_x_yyzz_yz, tr_x_yyzz_zz, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzz_xx[i] = tr_x_xzz_xx[i] * fe_0 + tr_x_xyzz_xx[i] * pa_y[i];

        tr_x_xyyzz_xy[i] = tr_x_xyy_xy[i] * fe_0 + tr_x_xyyz_xy[i] * pa_z[i];

        tr_x_xyyzz_xz[i] = tr_x_xzz_xz[i] * fe_0 + tr_x_xyzz_xz[i] * pa_y[i];

        tr_x_xyyzz_yy[i] = ts_yyzz_yy[i] * fe_0 + tr_x_yyzz_yy[i] * pa_x[i];

        tr_x_xyyzz_yz[i] = ts_yyzz_yz[i] * fe_0 + tr_x_yyzz_yz[i] * pa_x[i];

        tr_x_xyyzz_zz[i] = ts_yyzz_zz[i] * fe_0 + tr_x_yyzz_zz[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : HD

    auto tr_x_xyzzz_xx = pbuffer.data(idx_dip_hd + 78);

    auto tr_x_xyzzz_xy = pbuffer.data(idx_dip_hd + 79);

    auto tr_x_xyzzz_xz = pbuffer.data(idx_dip_hd + 80);

    auto tr_x_xyzzz_yy = pbuffer.data(idx_dip_hd + 81);

    auto tr_x_xyzzz_yz = pbuffer.data(idx_dip_hd + 82);

    auto tr_x_xyzzz_zz = pbuffer.data(idx_dip_hd + 83);

    #pragma omp simd aligned(pa_x, pa_y, tr_x_xyzzz_xx, tr_x_xyzzz_xy, tr_x_xyzzz_xz, tr_x_xyzzz_yy, tr_x_xyzzz_yz, tr_x_xyzzz_zz, tr_x_xzzz_x, tr_x_xzzz_xx, tr_x_xzzz_xy, tr_x_xzzz_xz, tr_x_xzzz_zz, tr_x_yzzz_yy, tr_x_yzzz_yz, ts_yzzz_yy, ts_yzzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzz_xx[i] = tr_x_xzzz_xx[i] * pa_y[i];

        tr_x_xyzzz_xy[i] = tr_x_xzzz_x[i] * fe_0 + tr_x_xzzz_xy[i] * pa_y[i];

        tr_x_xyzzz_xz[i] = tr_x_xzzz_xz[i] * pa_y[i];

        tr_x_xyzzz_yy[i] = ts_yzzz_yy[i] * fe_0 + tr_x_yzzz_yy[i] * pa_x[i];

        tr_x_xyzzz_yz[i] = ts_yzzz_yz[i] * fe_0 + tr_x_yzzz_yz[i] * pa_x[i];

        tr_x_xyzzz_zz[i] = tr_x_xzzz_zz[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : HD

    auto tr_x_xzzzz_xx = pbuffer.data(idx_dip_hd + 84);

    auto tr_x_xzzzz_xy = pbuffer.data(idx_dip_hd + 85);

    auto tr_x_xzzzz_xz = pbuffer.data(idx_dip_hd + 86);

    auto tr_x_xzzzz_yy = pbuffer.data(idx_dip_hd + 87);

    auto tr_x_xzzzz_yz = pbuffer.data(idx_dip_hd + 88);

    auto tr_x_xzzzz_zz = pbuffer.data(idx_dip_hd + 89);

    #pragma omp simd aligned(pa_x, pa_z, tr_x_xzz_xx, tr_x_xzz_xy, tr_x_xzzz_xx, tr_x_xzzz_xy, tr_x_xzzzz_xx, tr_x_xzzzz_xy, tr_x_xzzzz_xz, tr_x_xzzzz_yy, tr_x_xzzzz_yz, tr_x_xzzzz_zz, tr_x_zzzz_xz, tr_x_zzzz_yy, tr_x_zzzz_yz, tr_x_zzzz_z, tr_x_zzzz_zz, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzz_xx[i] = 3.0 * tr_x_xzz_xx[i] * fe_0 + tr_x_xzzz_xx[i] * pa_z[i];

        tr_x_xzzzz_xy[i] = 3.0 * tr_x_xzz_xy[i] * fe_0 + tr_x_xzzz_xy[i] * pa_z[i];

        tr_x_xzzzz_xz[i] = tr_x_zzzz_z[i] * fe_0 + ts_zzzz_xz[i] * fe_0 + tr_x_zzzz_xz[i] * pa_x[i];

        tr_x_xzzzz_yy[i] = ts_zzzz_yy[i] * fe_0 + tr_x_zzzz_yy[i] * pa_x[i];

        tr_x_xzzzz_yz[i] = ts_zzzz_yz[i] * fe_0 + tr_x_zzzz_yz[i] * pa_x[i];

        tr_x_xzzzz_zz[i] = ts_zzzz_zz[i] * fe_0 + tr_x_zzzz_zz[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : HD

    auto tr_x_yyyyy_xx = pbuffer.data(idx_dip_hd + 90);

    auto tr_x_yyyyy_xy = pbuffer.data(idx_dip_hd + 91);

    auto tr_x_yyyyy_xz = pbuffer.data(idx_dip_hd + 92);

    auto tr_x_yyyyy_yy = pbuffer.data(idx_dip_hd + 93);

    auto tr_x_yyyyy_yz = pbuffer.data(idx_dip_hd + 94);

    auto tr_x_yyyyy_zz = pbuffer.data(idx_dip_hd + 95);

    #pragma omp simd aligned(pa_y, tr_x_yyy_xx, tr_x_yyy_xy, tr_x_yyy_xz, tr_x_yyy_yy, tr_x_yyy_yz, tr_x_yyy_zz, tr_x_yyyy_x, tr_x_yyyy_xx, tr_x_yyyy_xy, tr_x_yyyy_xz, tr_x_yyyy_y, tr_x_yyyy_yy, tr_x_yyyy_yz, tr_x_yyyy_z, tr_x_yyyy_zz, tr_x_yyyyy_xx, tr_x_yyyyy_xy, tr_x_yyyyy_xz, tr_x_yyyyy_yy, tr_x_yyyyy_yz, tr_x_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyy_xx[i] = 4.0 * tr_x_yyy_xx[i] * fe_0 + tr_x_yyyy_xx[i] * pa_y[i];

        tr_x_yyyyy_xy[i] = 4.0 * tr_x_yyy_xy[i] * fe_0 + tr_x_yyyy_x[i] * fe_0 + tr_x_yyyy_xy[i] * pa_y[i];

        tr_x_yyyyy_xz[i] = 4.0 * tr_x_yyy_xz[i] * fe_0 + tr_x_yyyy_xz[i] * pa_y[i];

        tr_x_yyyyy_yy[i] = 4.0 * tr_x_yyy_yy[i] * fe_0 + 2.0 * tr_x_yyyy_y[i] * fe_0 + tr_x_yyyy_yy[i] * pa_y[i];

        tr_x_yyyyy_yz[i] = 4.0 * tr_x_yyy_yz[i] * fe_0 + tr_x_yyyy_z[i] * fe_0 + tr_x_yyyy_yz[i] * pa_y[i];

        tr_x_yyyyy_zz[i] = 4.0 * tr_x_yyy_zz[i] * fe_0 + tr_x_yyyy_zz[i] * pa_y[i];
    }

    // Set up 96-102 components of targeted buffer : HD

    auto tr_x_yyyyz_xx = pbuffer.data(idx_dip_hd + 96);

    auto tr_x_yyyyz_xy = pbuffer.data(idx_dip_hd + 97);

    auto tr_x_yyyyz_xz = pbuffer.data(idx_dip_hd + 98);

    auto tr_x_yyyyz_yy = pbuffer.data(idx_dip_hd + 99);

    auto tr_x_yyyyz_yz = pbuffer.data(idx_dip_hd + 100);

    auto tr_x_yyyyz_zz = pbuffer.data(idx_dip_hd + 101);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyyy_xx, tr_x_yyyy_xy, tr_x_yyyy_y, tr_x_yyyy_yy, tr_x_yyyy_yz, tr_x_yyyyz_xx, tr_x_yyyyz_xy, tr_x_yyyyz_xz, tr_x_yyyyz_yy, tr_x_yyyyz_yz, tr_x_yyyyz_zz, tr_x_yyyz_xz, tr_x_yyyz_zz, tr_x_yyz_xz, tr_x_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyz_xx[i] = tr_x_yyyy_xx[i] * pa_z[i];

        tr_x_yyyyz_xy[i] = tr_x_yyyy_xy[i] * pa_z[i];

        tr_x_yyyyz_xz[i] = 3.0 * tr_x_yyz_xz[i] * fe_0 + tr_x_yyyz_xz[i] * pa_y[i];

        tr_x_yyyyz_yy[i] = tr_x_yyyy_yy[i] * pa_z[i];

        tr_x_yyyyz_yz[i] = tr_x_yyyy_y[i] * fe_0 + tr_x_yyyy_yz[i] * pa_z[i];

        tr_x_yyyyz_zz[i] = 3.0 * tr_x_yyz_zz[i] * fe_0 + tr_x_yyyz_zz[i] * pa_y[i];
    }

    // Set up 102-108 components of targeted buffer : HD

    auto tr_x_yyyzz_xx = pbuffer.data(idx_dip_hd + 102);

    auto tr_x_yyyzz_xy = pbuffer.data(idx_dip_hd + 103);

    auto tr_x_yyyzz_xz = pbuffer.data(idx_dip_hd + 104);

    auto tr_x_yyyzz_yy = pbuffer.data(idx_dip_hd + 105);

    auto tr_x_yyyzz_yz = pbuffer.data(idx_dip_hd + 106);

    auto tr_x_yyyzz_zz = pbuffer.data(idx_dip_hd + 107);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyy_xy, tr_x_yyy_yy, tr_x_yyyz_xy, tr_x_yyyz_yy, tr_x_yyyzz_xx, tr_x_yyyzz_xy, tr_x_yyyzz_xz, tr_x_yyyzz_yy, tr_x_yyyzz_yz, tr_x_yyyzz_zz, tr_x_yyzz_xx, tr_x_yyzz_xz, tr_x_yyzz_yz, tr_x_yyzz_z, tr_x_yyzz_zz, tr_x_yzz_xx, tr_x_yzz_xz, tr_x_yzz_yz, tr_x_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzz_xx[i] = 2.0 * tr_x_yzz_xx[i] * fe_0 + tr_x_yyzz_xx[i] * pa_y[i];

        tr_x_yyyzz_xy[i] = tr_x_yyy_xy[i] * fe_0 + tr_x_yyyz_xy[i] * pa_z[i];

        tr_x_yyyzz_xz[i] = 2.0 * tr_x_yzz_xz[i] * fe_0 + tr_x_yyzz_xz[i] * pa_y[i];

        tr_x_yyyzz_yy[i] = tr_x_yyy_yy[i] * fe_0 + tr_x_yyyz_yy[i] * pa_z[i];

        tr_x_yyyzz_yz[i] = 2.0 * tr_x_yzz_yz[i] * fe_0 + tr_x_yyzz_z[i] * fe_0 + tr_x_yyzz_yz[i] * pa_y[i];

        tr_x_yyyzz_zz[i] = 2.0 * tr_x_yzz_zz[i] * fe_0 + tr_x_yyzz_zz[i] * pa_y[i];
    }

    // Set up 108-114 components of targeted buffer : HD

    auto tr_x_yyzzz_xx = pbuffer.data(idx_dip_hd + 108);

    auto tr_x_yyzzz_xy = pbuffer.data(idx_dip_hd + 109);

    auto tr_x_yyzzz_xz = pbuffer.data(idx_dip_hd + 110);

    auto tr_x_yyzzz_yy = pbuffer.data(idx_dip_hd + 111);

    auto tr_x_yyzzz_yz = pbuffer.data(idx_dip_hd + 112);

    auto tr_x_yyzzz_zz = pbuffer.data(idx_dip_hd + 113);

    #pragma omp simd aligned(pa_y, pa_z, tr_x_yyz_xy, tr_x_yyz_yy, tr_x_yyzz_xy, tr_x_yyzz_yy, tr_x_yyzzz_xx, tr_x_yyzzz_xy, tr_x_yyzzz_xz, tr_x_yyzzz_yy, tr_x_yyzzz_yz, tr_x_yyzzz_zz, tr_x_yzzz_xx, tr_x_yzzz_xz, tr_x_yzzz_yz, tr_x_yzzz_z, tr_x_yzzz_zz, tr_x_zzz_xx, tr_x_zzz_xz, tr_x_zzz_yz, tr_x_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzz_xx[i] = tr_x_zzz_xx[i] * fe_0 + tr_x_yzzz_xx[i] * pa_y[i];

        tr_x_yyzzz_xy[i] = 2.0 * tr_x_yyz_xy[i] * fe_0 + tr_x_yyzz_xy[i] * pa_z[i];

        tr_x_yyzzz_xz[i] = tr_x_zzz_xz[i] * fe_0 + tr_x_yzzz_xz[i] * pa_y[i];

        tr_x_yyzzz_yy[i] = 2.0 * tr_x_yyz_yy[i] * fe_0 + tr_x_yyzz_yy[i] * pa_z[i];

        tr_x_yyzzz_yz[i] = tr_x_zzz_yz[i] * fe_0 + tr_x_yzzz_z[i] * fe_0 + tr_x_yzzz_yz[i] * pa_y[i];

        tr_x_yyzzz_zz[i] = tr_x_zzz_zz[i] * fe_0 + tr_x_yzzz_zz[i] * pa_y[i];
    }

    // Set up 114-120 components of targeted buffer : HD

    auto tr_x_yzzzz_xx = pbuffer.data(idx_dip_hd + 114);

    auto tr_x_yzzzz_xy = pbuffer.data(idx_dip_hd + 115);

    auto tr_x_yzzzz_xz = pbuffer.data(idx_dip_hd + 116);

    auto tr_x_yzzzz_yy = pbuffer.data(idx_dip_hd + 117);

    auto tr_x_yzzzz_yz = pbuffer.data(idx_dip_hd + 118);

    auto tr_x_yzzzz_zz = pbuffer.data(idx_dip_hd + 119);

    #pragma omp simd aligned(pa_y, tr_x_yzzzz_xx, tr_x_yzzzz_xy, tr_x_yzzzz_xz, tr_x_yzzzz_yy, tr_x_yzzzz_yz, tr_x_yzzzz_zz, tr_x_zzzz_x, tr_x_zzzz_xx, tr_x_zzzz_xy, tr_x_zzzz_xz, tr_x_zzzz_y, tr_x_zzzz_yy, tr_x_zzzz_yz, tr_x_zzzz_z, tr_x_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzz_xx[i] = tr_x_zzzz_xx[i] * pa_y[i];

        tr_x_yzzzz_xy[i] = tr_x_zzzz_x[i] * fe_0 + tr_x_zzzz_xy[i] * pa_y[i];

        tr_x_yzzzz_xz[i] = tr_x_zzzz_xz[i] * pa_y[i];

        tr_x_yzzzz_yy[i] = 2.0 * tr_x_zzzz_y[i] * fe_0 + tr_x_zzzz_yy[i] * pa_y[i];

        tr_x_yzzzz_yz[i] = tr_x_zzzz_z[i] * fe_0 + tr_x_zzzz_yz[i] * pa_y[i];

        tr_x_yzzzz_zz[i] = tr_x_zzzz_zz[i] * pa_y[i];
    }

    // Set up 120-126 components of targeted buffer : HD

    auto tr_x_zzzzz_xx = pbuffer.data(idx_dip_hd + 120);

    auto tr_x_zzzzz_xy = pbuffer.data(idx_dip_hd + 121);

    auto tr_x_zzzzz_xz = pbuffer.data(idx_dip_hd + 122);

    auto tr_x_zzzzz_yy = pbuffer.data(idx_dip_hd + 123);

    auto tr_x_zzzzz_yz = pbuffer.data(idx_dip_hd + 124);

    auto tr_x_zzzzz_zz = pbuffer.data(idx_dip_hd + 125);

    #pragma omp simd aligned(pa_z, tr_x_zzz_xx, tr_x_zzz_xy, tr_x_zzz_xz, tr_x_zzz_yy, tr_x_zzz_yz, tr_x_zzz_zz, tr_x_zzzz_x, tr_x_zzzz_xx, tr_x_zzzz_xy, tr_x_zzzz_xz, tr_x_zzzz_y, tr_x_zzzz_yy, tr_x_zzzz_yz, tr_x_zzzz_z, tr_x_zzzz_zz, tr_x_zzzzz_xx, tr_x_zzzzz_xy, tr_x_zzzzz_xz, tr_x_zzzzz_yy, tr_x_zzzzz_yz, tr_x_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzz_xx[i] = 4.0 * tr_x_zzz_xx[i] * fe_0 + tr_x_zzzz_xx[i] * pa_z[i];

        tr_x_zzzzz_xy[i] = 4.0 * tr_x_zzz_xy[i] * fe_0 + tr_x_zzzz_xy[i] * pa_z[i];

        tr_x_zzzzz_xz[i] = 4.0 * tr_x_zzz_xz[i] * fe_0 + tr_x_zzzz_x[i] * fe_0 + tr_x_zzzz_xz[i] * pa_z[i];

        tr_x_zzzzz_yy[i] = 4.0 * tr_x_zzz_yy[i] * fe_0 + tr_x_zzzz_yy[i] * pa_z[i];

        tr_x_zzzzz_yz[i] = 4.0 * tr_x_zzz_yz[i] * fe_0 + tr_x_zzzz_y[i] * fe_0 + tr_x_zzzz_yz[i] * pa_z[i];

        tr_x_zzzzz_zz[i] = 4.0 * tr_x_zzz_zz[i] * fe_0 + 2.0 * tr_x_zzzz_z[i] * fe_0 + tr_x_zzzz_zz[i] * pa_z[i];
    }

    // Set up 126-132 components of targeted buffer : HD

    auto tr_y_xxxxx_xx = pbuffer.data(idx_dip_hd + 126);

    auto tr_y_xxxxx_xy = pbuffer.data(idx_dip_hd + 127);

    auto tr_y_xxxxx_xz = pbuffer.data(idx_dip_hd + 128);

    auto tr_y_xxxxx_yy = pbuffer.data(idx_dip_hd + 129);

    auto tr_y_xxxxx_yz = pbuffer.data(idx_dip_hd + 130);

    auto tr_y_xxxxx_zz = pbuffer.data(idx_dip_hd + 131);

    #pragma omp simd aligned(pa_x, tr_y_xxx_xx, tr_y_xxx_xy, tr_y_xxx_xz, tr_y_xxx_yy, tr_y_xxx_yz, tr_y_xxx_zz, tr_y_xxxx_x, tr_y_xxxx_xx, tr_y_xxxx_xy, tr_y_xxxx_xz, tr_y_xxxx_y, tr_y_xxxx_yy, tr_y_xxxx_yz, tr_y_xxxx_z, tr_y_xxxx_zz, tr_y_xxxxx_xx, tr_y_xxxxx_xy, tr_y_xxxxx_xz, tr_y_xxxxx_yy, tr_y_xxxxx_yz, tr_y_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxx_xx[i] = 4.0 * tr_y_xxx_xx[i] * fe_0 + 2.0 * tr_y_xxxx_x[i] * fe_0 + tr_y_xxxx_xx[i] * pa_x[i];

        tr_y_xxxxx_xy[i] = 4.0 * tr_y_xxx_xy[i] * fe_0 + tr_y_xxxx_y[i] * fe_0 + tr_y_xxxx_xy[i] * pa_x[i];

        tr_y_xxxxx_xz[i] = 4.0 * tr_y_xxx_xz[i] * fe_0 + tr_y_xxxx_z[i] * fe_0 + tr_y_xxxx_xz[i] * pa_x[i];

        tr_y_xxxxx_yy[i] = 4.0 * tr_y_xxx_yy[i] * fe_0 + tr_y_xxxx_yy[i] * pa_x[i];

        tr_y_xxxxx_yz[i] = 4.0 * tr_y_xxx_yz[i] * fe_0 + tr_y_xxxx_yz[i] * pa_x[i];

        tr_y_xxxxx_zz[i] = 4.0 * tr_y_xxx_zz[i] * fe_0 + tr_y_xxxx_zz[i] * pa_x[i];
    }

    // Set up 132-138 components of targeted buffer : HD

    auto tr_y_xxxxy_xx = pbuffer.data(idx_dip_hd + 132);

    auto tr_y_xxxxy_xy = pbuffer.data(idx_dip_hd + 133);

    auto tr_y_xxxxy_xz = pbuffer.data(idx_dip_hd + 134);

    auto tr_y_xxxxy_yy = pbuffer.data(idx_dip_hd + 135);

    auto tr_y_xxxxy_yz = pbuffer.data(idx_dip_hd + 136);

    auto tr_y_xxxxy_zz = pbuffer.data(idx_dip_hd + 137);

    #pragma omp simd aligned(pa_x, pa_y, tr_y_xxxx_xx, tr_y_xxxx_xz, tr_y_xxxxy_xx, tr_y_xxxxy_xy, tr_y_xxxxy_xz, tr_y_xxxxy_yy, tr_y_xxxxy_yz, tr_y_xxxxy_zz, tr_y_xxxy_xy, tr_y_xxxy_y, tr_y_xxxy_yy, tr_y_xxxy_yz, tr_y_xxxy_zz, tr_y_xxy_xy, tr_y_xxy_yy, tr_y_xxy_yz, tr_y_xxy_zz, ts_xxxx_xx, ts_xxxx_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxy_xx[i] = ts_xxxx_xx[i] * fe_0 + tr_y_xxxx_xx[i] * pa_y[i];

        tr_y_xxxxy_xy[i] = 3.0 * tr_y_xxy_xy[i] * fe_0 + tr_y_xxxy_y[i] * fe_0 + tr_y_xxxy_xy[i] * pa_x[i];

        tr_y_xxxxy_xz[i] = ts_xxxx_xz[i] * fe_0 + tr_y_xxxx_xz[i] * pa_y[i];

        tr_y_xxxxy_yy[i] = 3.0 * tr_y_xxy_yy[i] * fe_0 + tr_y_xxxy_yy[i] * pa_x[i];

        tr_y_xxxxy_yz[i] = 3.0 * tr_y_xxy_yz[i] * fe_0 + tr_y_xxxy_yz[i] * pa_x[i];

        tr_y_xxxxy_zz[i] = 3.0 * tr_y_xxy_zz[i] * fe_0 + tr_y_xxxy_zz[i] * pa_x[i];
    }

    // Set up 138-144 components of targeted buffer : HD

    auto tr_y_xxxxz_xx = pbuffer.data(idx_dip_hd + 138);

    auto tr_y_xxxxz_xy = pbuffer.data(idx_dip_hd + 139);

    auto tr_y_xxxxz_xz = pbuffer.data(idx_dip_hd + 140);

    auto tr_y_xxxxz_yy = pbuffer.data(idx_dip_hd + 141);

    auto tr_y_xxxxz_yz = pbuffer.data(idx_dip_hd + 142);

    auto tr_y_xxxxz_zz = pbuffer.data(idx_dip_hd + 143);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxxx_x, tr_y_xxxx_xx, tr_y_xxxx_xy, tr_y_xxxx_xz, tr_y_xxxx_yy, tr_y_xxxxz_xx, tr_y_xxxxz_xy, tr_y_xxxxz_xz, tr_y_xxxxz_yy, tr_y_xxxxz_yz, tr_y_xxxxz_zz, tr_y_xxxz_yz, tr_y_xxxz_zz, tr_y_xxz_yz, tr_y_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxz_xx[i] = tr_y_xxxx_xx[i] * pa_z[i];

        tr_y_xxxxz_xy[i] = tr_y_xxxx_xy[i] * pa_z[i];

        tr_y_xxxxz_xz[i] = tr_y_xxxx_x[i] * fe_0 + tr_y_xxxx_xz[i] * pa_z[i];

        tr_y_xxxxz_yy[i] = tr_y_xxxx_yy[i] * pa_z[i];

        tr_y_xxxxz_yz[i] = 3.0 * tr_y_xxz_yz[i] * fe_0 + tr_y_xxxz_yz[i] * pa_x[i];

        tr_y_xxxxz_zz[i] = 3.0 * tr_y_xxz_zz[i] * fe_0 + tr_y_xxxz_zz[i] * pa_x[i];
    }

    // Set up 144-150 components of targeted buffer : HD

    auto tr_y_xxxyy_xx = pbuffer.data(idx_dip_hd + 144);

    auto tr_y_xxxyy_xy = pbuffer.data(idx_dip_hd + 145);

    auto tr_y_xxxyy_xz = pbuffer.data(idx_dip_hd + 146);

    auto tr_y_xxxyy_yy = pbuffer.data(idx_dip_hd + 147);

    auto tr_y_xxxyy_yz = pbuffer.data(idx_dip_hd + 148);

    auto tr_y_xxxyy_zz = pbuffer.data(idx_dip_hd + 149);

    #pragma omp simd aligned(pa_x, tr_y_xxxyy_xx, tr_y_xxxyy_xy, tr_y_xxxyy_xz, tr_y_xxxyy_yy, tr_y_xxxyy_yz, tr_y_xxxyy_zz, tr_y_xxyy_x, tr_y_xxyy_xx, tr_y_xxyy_xy, tr_y_xxyy_xz, tr_y_xxyy_y, tr_y_xxyy_yy, tr_y_xxyy_yz, tr_y_xxyy_z, tr_y_xxyy_zz, tr_y_xyy_xx, tr_y_xyy_xy, tr_y_xyy_xz, tr_y_xyy_yy, tr_y_xyy_yz, tr_y_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyy_xx[i] = 2.0 * tr_y_xyy_xx[i] * fe_0 + 2.0 * tr_y_xxyy_x[i] * fe_0 + tr_y_xxyy_xx[i] * pa_x[i];

        tr_y_xxxyy_xy[i] = 2.0 * tr_y_xyy_xy[i] * fe_0 + tr_y_xxyy_y[i] * fe_0 + tr_y_xxyy_xy[i] * pa_x[i];

        tr_y_xxxyy_xz[i] = 2.0 * tr_y_xyy_xz[i] * fe_0 + tr_y_xxyy_z[i] * fe_0 + tr_y_xxyy_xz[i] * pa_x[i];

        tr_y_xxxyy_yy[i] = 2.0 * tr_y_xyy_yy[i] * fe_0 + tr_y_xxyy_yy[i] * pa_x[i];

        tr_y_xxxyy_yz[i] = 2.0 * tr_y_xyy_yz[i] * fe_0 + tr_y_xxyy_yz[i] * pa_x[i];

        tr_y_xxxyy_zz[i] = 2.0 * tr_y_xyy_zz[i] * fe_0 + tr_y_xxyy_zz[i] * pa_x[i];
    }

    // Set up 150-156 components of targeted buffer : HD

    auto tr_y_xxxyz_xx = pbuffer.data(idx_dip_hd + 150);

    auto tr_y_xxxyz_xy = pbuffer.data(idx_dip_hd + 151);

    auto tr_y_xxxyz_xz = pbuffer.data(idx_dip_hd + 152);

    auto tr_y_xxxyz_yy = pbuffer.data(idx_dip_hd + 153);

    auto tr_y_xxxyz_yz = pbuffer.data(idx_dip_hd + 154);

    auto tr_y_xxxyz_zz = pbuffer.data(idx_dip_hd + 155);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxxy_xx, tr_y_xxxy_xy, tr_y_xxxy_yy, tr_y_xxxyz_xx, tr_y_xxxyz_xy, tr_y_xxxyz_xz, tr_y_xxxyz_yy, tr_y_xxxyz_yz, tr_y_xxxyz_zz, tr_y_xxxz_xz, tr_y_xxyz_yz, tr_y_xxyz_zz, tr_y_xyz_yz, tr_y_xyz_zz, ts_xxxz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyz_xx[i] = tr_y_xxxy_xx[i] * pa_z[i];

        tr_y_xxxyz_xy[i] = tr_y_xxxy_xy[i] * pa_z[i];

        tr_y_xxxyz_xz[i] = ts_xxxz_xz[i] * fe_0 + tr_y_xxxz_xz[i] * pa_y[i];

        tr_y_xxxyz_yy[i] = tr_y_xxxy_yy[i] * pa_z[i];

        tr_y_xxxyz_yz[i] = 2.0 * tr_y_xyz_yz[i] * fe_0 + tr_y_xxyz_yz[i] * pa_x[i];

        tr_y_xxxyz_zz[i] = 2.0 * tr_y_xyz_zz[i] * fe_0 + tr_y_xxyz_zz[i] * pa_x[i];
    }

    // Set up 156-162 components of targeted buffer : HD

    auto tr_y_xxxzz_xx = pbuffer.data(idx_dip_hd + 156);

    auto tr_y_xxxzz_xy = pbuffer.data(idx_dip_hd + 157);

    auto tr_y_xxxzz_xz = pbuffer.data(idx_dip_hd + 158);

    auto tr_y_xxxzz_yy = pbuffer.data(idx_dip_hd + 159);

    auto tr_y_xxxzz_yz = pbuffer.data(idx_dip_hd + 160);

    auto tr_y_xxxzz_zz = pbuffer.data(idx_dip_hd + 161);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxx_xx, tr_y_xxx_xy, tr_y_xxxz_xx, tr_y_xxxz_xy, tr_y_xxxzz_xx, tr_y_xxxzz_xy, tr_y_xxxzz_xz, tr_y_xxxzz_yy, tr_y_xxxzz_yz, tr_y_xxxzz_zz, tr_y_xxzz_xz, tr_y_xxzz_yy, tr_y_xxzz_yz, tr_y_xxzz_z, tr_y_xxzz_zz, tr_y_xzz_xz, tr_y_xzz_yy, tr_y_xzz_yz, tr_y_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzz_xx[i] = tr_y_xxx_xx[i] * fe_0 + tr_y_xxxz_xx[i] * pa_z[i];

        tr_y_xxxzz_xy[i] = tr_y_xxx_xy[i] * fe_0 + tr_y_xxxz_xy[i] * pa_z[i];

        tr_y_xxxzz_xz[i] = 2.0 * tr_y_xzz_xz[i] * fe_0 + tr_y_xxzz_z[i] * fe_0 + tr_y_xxzz_xz[i] * pa_x[i];

        tr_y_xxxzz_yy[i] = 2.0 * tr_y_xzz_yy[i] * fe_0 + tr_y_xxzz_yy[i] * pa_x[i];

        tr_y_xxxzz_yz[i] = 2.0 * tr_y_xzz_yz[i] * fe_0 + tr_y_xxzz_yz[i] * pa_x[i];

        tr_y_xxxzz_zz[i] = 2.0 * tr_y_xzz_zz[i] * fe_0 + tr_y_xxzz_zz[i] * pa_x[i];
    }

    // Set up 162-168 components of targeted buffer : HD

    auto tr_y_xxyyy_xx = pbuffer.data(idx_dip_hd + 162);

    auto tr_y_xxyyy_xy = pbuffer.data(idx_dip_hd + 163);

    auto tr_y_xxyyy_xz = pbuffer.data(idx_dip_hd + 164);

    auto tr_y_xxyyy_yy = pbuffer.data(idx_dip_hd + 165);

    auto tr_y_xxyyy_yz = pbuffer.data(idx_dip_hd + 166);

    auto tr_y_xxyyy_zz = pbuffer.data(idx_dip_hd + 167);

    #pragma omp simd aligned(pa_x, tr_y_xxyyy_xx, tr_y_xxyyy_xy, tr_y_xxyyy_xz, tr_y_xxyyy_yy, tr_y_xxyyy_yz, tr_y_xxyyy_zz, tr_y_xyyy_x, tr_y_xyyy_xx, tr_y_xyyy_xy, tr_y_xyyy_xz, tr_y_xyyy_y, tr_y_xyyy_yy, tr_y_xyyy_yz, tr_y_xyyy_z, tr_y_xyyy_zz, tr_y_yyy_xx, tr_y_yyy_xy, tr_y_yyy_xz, tr_y_yyy_yy, tr_y_yyy_yz, tr_y_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyy_xx[i] = tr_y_yyy_xx[i] * fe_0 + 2.0 * tr_y_xyyy_x[i] * fe_0 + tr_y_xyyy_xx[i] * pa_x[i];

        tr_y_xxyyy_xy[i] = tr_y_yyy_xy[i] * fe_0 + tr_y_xyyy_y[i] * fe_0 + tr_y_xyyy_xy[i] * pa_x[i];

        tr_y_xxyyy_xz[i] = tr_y_yyy_xz[i] * fe_0 + tr_y_xyyy_z[i] * fe_0 + tr_y_xyyy_xz[i] * pa_x[i];

        tr_y_xxyyy_yy[i] = tr_y_yyy_yy[i] * fe_0 + tr_y_xyyy_yy[i] * pa_x[i];

        tr_y_xxyyy_yz[i] = tr_y_yyy_yz[i] * fe_0 + tr_y_xyyy_yz[i] * pa_x[i];

        tr_y_xxyyy_zz[i] = tr_y_yyy_zz[i] * fe_0 + tr_y_xyyy_zz[i] * pa_x[i];
    }

    // Set up 168-174 components of targeted buffer : HD

    auto tr_y_xxyyz_xx = pbuffer.data(idx_dip_hd + 168);

    auto tr_y_xxyyz_xy = pbuffer.data(idx_dip_hd + 169);

    auto tr_y_xxyyz_xz = pbuffer.data(idx_dip_hd + 170);

    auto tr_y_xxyyz_yy = pbuffer.data(idx_dip_hd + 171);

    auto tr_y_xxyyz_yz = pbuffer.data(idx_dip_hd + 172);

    auto tr_y_xxyyz_zz = pbuffer.data(idx_dip_hd + 173);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxyy_x, tr_y_xxyy_xx, tr_y_xxyy_xy, tr_y_xxyy_xz, tr_y_xxyy_yy, tr_y_xxyyz_xx, tr_y_xxyyz_xy, tr_y_xxyyz_xz, tr_y_xxyyz_yy, tr_y_xxyyz_yz, tr_y_xxyyz_zz, tr_y_xyyz_yz, tr_y_xyyz_zz, tr_y_yyz_yz, tr_y_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyz_xx[i] = tr_y_xxyy_xx[i] * pa_z[i];

        tr_y_xxyyz_xy[i] = tr_y_xxyy_xy[i] * pa_z[i];

        tr_y_xxyyz_xz[i] = tr_y_xxyy_x[i] * fe_0 + tr_y_xxyy_xz[i] * pa_z[i];

        tr_y_xxyyz_yy[i] = tr_y_xxyy_yy[i] * pa_z[i];

        tr_y_xxyyz_yz[i] = tr_y_yyz_yz[i] * fe_0 + tr_y_xyyz_yz[i] * pa_x[i];

        tr_y_xxyyz_zz[i] = tr_y_yyz_zz[i] * fe_0 + tr_y_xyyz_zz[i] * pa_x[i];
    }

    // Set up 174-180 components of targeted buffer : HD

    auto tr_y_xxyzz_xx = pbuffer.data(idx_dip_hd + 174);

    auto tr_y_xxyzz_xy = pbuffer.data(idx_dip_hd + 175);

    auto tr_y_xxyzz_xz = pbuffer.data(idx_dip_hd + 176);

    auto tr_y_xxyzz_yy = pbuffer.data(idx_dip_hd + 177);

    auto tr_y_xxyzz_yz = pbuffer.data(idx_dip_hd + 178);

    auto tr_y_xxyzz_zz = pbuffer.data(idx_dip_hd + 179);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_y_xxy_xy, tr_y_xxyz_xy, tr_y_xxyzz_xx, tr_y_xxyzz_xy, tr_y_xxyzz_xz, tr_y_xxyzz_yy, tr_y_xxyzz_yz, tr_y_xxyzz_zz, tr_y_xxzz_xx, tr_y_xxzz_xz, tr_y_xyzz_yy, tr_y_xyzz_yz, tr_y_xyzz_zz, tr_y_yzz_yy, tr_y_yzz_yz, tr_y_yzz_zz, ts_xxzz_xx, ts_xxzz_xz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzz_xx[i] = ts_xxzz_xx[i] * fe_0 + tr_y_xxzz_xx[i] * pa_y[i];

        tr_y_xxyzz_xy[i] = tr_y_xxy_xy[i] * fe_0 + tr_y_xxyz_xy[i] * pa_z[i];

        tr_y_xxyzz_xz[i] = ts_xxzz_xz[i] * fe_0 + tr_y_xxzz_xz[i] * pa_y[i];

        tr_y_xxyzz_yy[i] = tr_y_yzz_yy[i] * fe_0 + tr_y_xyzz_yy[i] * pa_x[i];

        tr_y_xxyzz_yz[i] = tr_y_yzz_yz[i] * fe_0 + tr_y_xyzz_yz[i] * pa_x[i];

        tr_y_xxyzz_zz[i] = tr_y_yzz_zz[i] * fe_0 + tr_y_xyzz_zz[i] * pa_x[i];
    }

    // Set up 180-186 components of targeted buffer : HD

    auto tr_y_xxzzz_xx = pbuffer.data(idx_dip_hd + 180);

    auto tr_y_xxzzz_xy = pbuffer.data(idx_dip_hd + 181);

    auto tr_y_xxzzz_xz = pbuffer.data(idx_dip_hd + 182);

    auto tr_y_xxzzz_yy = pbuffer.data(idx_dip_hd + 183);

    auto tr_y_xxzzz_yz = pbuffer.data(idx_dip_hd + 184);

    auto tr_y_xxzzz_zz = pbuffer.data(idx_dip_hd + 185);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xxz_xx, tr_y_xxz_xy, tr_y_xxzz_xx, tr_y_xxzz_xy, tr_y_xxzzz_xx, tr_y_xxzzz_xy, tr_y_xxzzz_xz, tr_y_xxzzz_yy, tr_y_xxzzz_yz, tr_y_xxzzz_zz, tr_y_xzzz_xz, tr_y_xzzz_yy, tr_y_xzzz_yz, tr_y_xzzz_z, tr_y_xzzz_zz, tr_y_zzz_xz, tr_y_zzz_yy, tr_y_zzz_yz, tr_y_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzz_xx[i] = 2.0 * tr_y_xxz_xx[i] * fe_0 + tr_y_xxzz_xx[i] * pa_z[i];

        tr_y_xxzzz_xy[i] = 2.0 * tr_y_xxz_xy[i] * fe_0 + tr_y_xxzz_xy[i] * pa_z[i];

        tr_y_xxzzz_xz[i] = tr_y_zzz_xz[i] * fe_0 + tr_y_xzzz_z[i] * fe_0 + tr_y_xzzz_xz[i] * pa_x[i];

        tr_y_xxzzz_yy[i] = tr_y_zzz_yy[i] * fe_0 + tr_y_xzzz_yy[i] * pa_x[i];

        tr_y_xxzzz_yz[i] = tr_y_zzz_yz[i] * fe_0 + tr_y_xzzz_yz[i] * pa_x[i];

        tr_y_xxzzz_zz[i] = tr_y_zzz_zz[i] * fe_0 + tr_y_xzzz_zz[i] * pa_x[i];
    }

    // Set up 186-192 components of targeted buffer : HD

    auto tr_y_xyyyy_xx = pbuffer.data(idx_dip_hd + 186);

    auto tr_y_xyyyy_xy = pbuffer.data(idx_dip_hd + 187);

    auto tr_y_xyyyy_xz = pbuffer.data(idx_dip_hd + 188);

    auto tr_y_xyyyy_yy = pbuffer.data(idx_dip_hd + 189);

    auto tr_y_xyyyy_yz = pbuffer.data(idx_dip_hd + 190);

    auto tr_y_xyyyy_zz = pbuffer.data(idx_dip_hd + 191);

    #pragma omp simd aligned(pa_x, tr_y_xyyyy_xx, tr_y_xyyyy_xy, tr_y_xyyyy_xz, tr_y_xyyyy_yy, tr_y_xyyyy_yz, tr_y_xyyyy_zz, tr_y_yyyy_x, tr_y_yyyy_xx, tr_y_yyyy_xy, tr_y_yyyy_xz, tr_y_yyyy_y, tr_y_yyyy_yy, tr_y_yyyy_yz, tr_y_yyyy_z, tr_y_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyy_xx[i] = 2.0 * tr_y_yyyy_x[i] * fe_0 + tr_y_yyyy_xx[i] * pa_x[i];

        tr_y_xyyyy_xy[i] = tr_y_yyyy_y[i] * fe_0 + tr_y_yyyy_xy[i] * pa_x[i];

        tr_y_xyyyy_xz[i] = tr_y_yyyy_z[i] * fe_0 + tr_y_yyyy_xz[i] * pa_x[i];

        tr_y_xyyyy_yy[i] = tr_y_yyyy_yy[i] * pa_x[i];

        tr_y_xyyyy_yz[i] = tr_y_yyyy_yz[i] * pa_x[i];

        tr_y_xyyyy_zz[i] = tr_y_yyyy_zz[i] * pa_x[i];
    }

    // Set up 192-198 components of targeted buffer : HD

    auto tr_y_xyyyz_xx = pbuffer.data(idx_dip_hd + 192);

    auto tr_y_xyyyz_xy = pbuffer.data(idx_dip_hd + 193);

    auto tr_y_xyyyz_xz = pbuffer.data(idx_dip_hd + 194);

    auto tr_y_xyyyz_yy = pbuffer.data(idx_dip_hd + 195);

    auto tr_y_xyyyz_yz = pbuffer.data(idx_dip_hd + 196);

    auto tr_y_xyyyz_zz = pbuffer.data(idx_dip_hd + 197);

    #pragma omp simd aligned(pa_x, pa_z, tr_y_xyyy_xx, tr_y_xyyy_xy, tr_y_xyyyz_xx, tr_y_xyyyz_xy, tr_y_xyyyz_xz, tr_y_xyyyz_yy, tr_y_xyyyz_yz, tr_y_xyyyz_zz, tr_y_yyyz_xz, tr_y_yyyz_yy, tr_y_yyyz_yz, tr_y_yyyz_z, tr_y_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyz_xx[i] = tr_y_xyyy_xx[i] * pa_z[i];

        tr_y_xyyyz_xy[i] = tr_y_xyyy_xy[i] * pa_z[i];

        tr_y_xyyyz_xz[i] = tr_y_yyyz_z[i] * fe_0 + tr_y_yyyz_xz[i] * pa_x[i];

        tr_y_xyyyz_yy[i] = tr_y_yyyz_yy[i] * pa_x[i];

        tr_y_xyyyz_yz[i] = tr_y_yyyz_yz[i] * pa_x[i];

        tr_y_xyyyz_zz[i] = tr_y_yyyz_zz[i] * pa_x[i];
    }

    // Set up 198-204 components of targeted buffer : HD

    auto tr_y_xyyzz_xx = pbuffer.data(idx_dip_hd + 198);

    auto tr_y_xyyzz_xy = pbuffer.data(idx_dip_hd + 199);

    auto tr_y_xyyzz_xz = pbuffer.data(idx_dip_hd + 200);

    auto tr_y_xyyzz_yy = pbuffer.data(idx_dip_hd + 201);

    auto tr_y_xyyzz_yz = pbuffer.data(idx_dip_hd + 202);

    auto tr_y_xyyzz_zz = pbuffer.data(idx_dip_hd + 203);

    #pragma omp simd aligned(pa_x, tr_y_xyyzz_xx, tr_y_xyyzz_xy, tr_y_xyyzz_xz, tr_y_xyyzz_yy, tr_y_xyyzz_yz, tr_y_xyyzz_zz, tr_y_yyzz_x, tr_y_yyzz_xx, tr_y_yyzz_xy, tr_y_yyzz_xz, tr_y_yyzz_y, tr_y_yyzz_yy, tr_y_yyzz_yz, tr_y_yyzz_z, tr_y_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzz_xx[i] = 2.0 * tr_y_yyzz_x[i] * fe_0 + tr_y_yyzz_xx[i] * pa_x[i];

        tr_y_xyyzz_xy[i] = tr_y_yyzz_y[i] * fe_0 + tr_y_yyzz_xy[i] * pa_x[i];

        tr_y_xyyzz_xz[i] = tr_y_yyzz_z[i] * fe_0 + tr_y_yyzz_xz[i] * pa_x[i];

        tr_y_xyyzz_yy[i] = tr_y_yyzz_yy[i] * pa_x[i];

        tr_y_xyyzz_yz[i] = tr_y_yyzz_yz[i] * pa_x[i];

        tr_y_xyyzz_zz[i] = tr_y_yyzz_zz[i] * pa_x[i];
    }

    // Set up 204-210 components of targeted buffer : HD

    auto tr_y_xyzzz_xx = pbuffer.data(idx_dip_hd + 204);

    auto tr_y_xyzzz_xy = pbuffer.data(idx_dip_hd + 205);

    auto tr_y_xyzzz_xz = pbuffer.data(idx_dip_hd + 206);

    auto tr_y_xyzzz_yy = pbuffer.data(idx_dip_hd + 207);

    auto tr_y_xyzzz_yz = pbuffer.data(idx_dip_hd + 208);

    auto tr_y_xyzzz_zz = pbuffer.data(idx_dip_hd + 209);

    #pragma omp simd aligned(pa_x, tr_y_xyzzz_xx, tr_y_xyzzz_xy, tr_y_xyzzz_xz, tr_y_xyzzz_yy, tr_y_xyzzz_yz, tr_y_xyzzz_zz, tr_y_yzzz_x, tr_y_yzzz_xx, tr_y_yzzz_xy, tr_y_yzzz_xz, tr_y_yzzz_y, tr_y_yzzz_yy, tr_y_yzzz_yz, tr_y_yzzz_z, tr_y_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzz_xx[i] = 2.0 * tr_y_yzzz_x[i] * fe_0 + tr_y_yzzz_xx[i] * pa_x[i];

        tr_y_xyzzz_xy[i] = tr_y_yzzz_y[i] * fe_0 + tr_y_yzzz_xy[i] * pa_x[i];

        tr_y_xyzzz_xz[i] = tr_y_yzzz_z[i] * fe_0 + tr_y_yzzz_xz[i] * pa_x[i];

        tr_y_xyzzz_yy[i] = tr_y_yzzz_yy[i] * pa_x[i];

        tr_y_xyzzz_yz[i] = tr_y_yzzz_yz[i] * pa_x[i];

        tr_y_xyzzz_zz[i] = tr_y_yzzz_zz[i] * pa_x[i];
    }

    // Set up 210-216 components of targeted buffer : HD

    auto tr_y_xzzzz_xx = pbuffer.data(idx_dip_hd + 210);

    auto tr_y_xzzzz_xy = pbuffer.data(idx_dip_hd + 211);

    auto tr_y_xzzzz_xz = pbuffer.data(idx_dip_hd + 212);

    auto tr_y_xzzzz_yy = pbuffer.data(idx_dip_hd + 213);

    auto tr_y_xzzzz_yz = pbuffer.data(idx_dip_hd + 214);

    auto tr_y_xzzzz_zz = pbuffer.data(idx_dip_hd + 215);

    #pragma omp simd aligned(pa_x, tr_y_xzzzz_xx, tr_y_xzzzz_xy, tr_y_xzzzz_xz, tr_y_xzzzz_yy, tr_y_xzzzz_yz, tr_y_xzzzz_zz, tr_y_zzzz_x, tr_y_zzzz_xx, tr_y_zzzz_xy, tr_y_zzzz_xz, tr_y_zzzz_y, tr_y_zzzz_yy, tr_y_zzzz_yz, tr_y_zzzz_z, tr_y_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzz_xx[i] = 2.0 * tr_y_zzzz_x[i] * fe_0 + tr_y_zzzz_xx[i] * pa_x[i];

        tr_y_xzzzz_xy[i] = tr_y_zzzz_y[i] * fe_0 + tr_y_zzzz_xy[i] * pa_x[i];

        tr_y_xzzzz_xz[i] = tr_y_zzzz_z[i] * fe_0 + tr_y_zzzz_xz[i] * pa_x[i];

        tr_y_xzzzz_yy[i] = tr_y_zzzz_yy[i] * pa_x[i];

        tr_y_xzzzz_yz[i] = tr_y_zzzz_yz[i] * pa_x[i];

        tr_y_xzzzz_zz[i] = tr_y_zzzz_zz[i] * pa_x[i];
    }

    // Set up 216-222 components of targeted buffer : HD

    auto tr_y_yyyyy_xx = pbuffer.data(idx_dip_hd + 216);

    auto tr_y_yyyyy_xy = pbuffer.data(idx_dip_hd + 217);

    auto tr_y_yyyyy_xz = pbuffer.data(idx_dip_hd + 218);

    auto tr_y_yyyyy_yy = pbuffer.data(idx_dip_hd + 219);

    auto tr_y_yyyyy_yz = pbuffer.data(idx_dip_hd + 220);

    auto tr_y_yyyyy_zz = pbuffer.data(idx_dip_hd + 221);

    #pragma omp simd aligned(pa_y, tr_y_yyy_xx, tr_y_yyy_xy, tr_y_yyy_xz, tr_y_yyy_yy, tr_y_yyy_yz, tr_y_yyy_zz, tr_y_yyyy_x, tr_y_yyyy_xx, tr_y_yyyy_xy, tr_y_yyyy_xz, tr_y_yyyy_y, tr_y_yyyy_yy, tr_y_yyyy_yz, tr_y_yyyy_z, tr_y_yyyy_zz, tr_y_yyyyy_xx, tr_y_yyyyy_xy, tr_y_yyyyy_xz, tr_y_yyyyy_yy, tr_y_yyyyy_yz, tr_y_yyyyy_zz, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyy_xx[i] = 4.0 * tr_y_yyy_xx[i] * fe_0 + ts_yyyy_xx[i] * fe_0 + tr_y_yyyy_xx[i] * pa_y[i];

        tr_y_yyyyy_xy[i] = 4.0 * tr_y_yyy_xy[i] * fe_0 + tr_y_yyyy_x[i] * fe_0 + ts_yyyy_xy[i] * fe_0 + tr_y_yyyy_xy[i] * pa_y[i];

        tr_y_yyyyy_xz[i] = 4.0 * tr_y_yyy_xz[i] * fe_0 + ts_yyyy_xz[i] * fe_0 + tr_y_yyyy_xz[i] * pa_y[i];

        tr_y_yyyyy_yy[i] = 4.0 * tr_y_yyy_yy[i] * fe_0 + 2.0 * tr_y_yyyy_y[i] * fe_0 + ts_yyyy_yy[i] * fe_0 + tr_y_yyyy_yy[i] * pa_y[i];

        tr_y_yyyyy_yz[i] = 4.0 * tr_y_yyy_yz[i] * fe_0 + tr_y_yyyy_z[i] * fe_0 + ts_yyyy_yz[i] * fe_0 + tr_y_yyyy_yz[i] * pa_y[i];

        tr_y_yyyyy_zz[i] = 4.0 * tr_y_yyy_zz[i] * fe_0 + ts_yyyy_zz[i] * fe_0 + tr_y_yyyy_zz[i] * pa_y[i];
    }

    // Set up 222-228 components of targeted buffer : HD

    auto tr_y_yyyyz_xx = pbuffer.data(idx_dip_hd + 222);

    auto tr_y_yyyyz_xy = pbuffer.data(idx_dip_hd + 223);

    auto tr_y_yyyyz_xz = pbuffer.data(idx_dip_hd + 224);

    auto tr_y_yyyyz_yy = pbuffer.data(idx_dip_hd + 225);

    auto tr_y_yyyyz_yz = pbuffer.data(idx_dip_hd + 226);

    auto tr_y_yyyyz_zz = pbuffer.data(idx_dip_hd + 227);

    #pragma omp simd aligned(pa_z, tr_y_yyyy_x, tr_y_yyyy_xx, tr_y_yyyy_xy, tr_y_yyyy_xz, tr_y_yyyy_y, tr_y_yyyy_yy, tr_y_yyyy_yz, tr_y_yyyy_z, tr_y_yyyy_zz, tr_y_yyyyz_xx, tr_y_yyyyz_xy, tr_y_yyyyz_xz, tr_y_yyyyz_yy, tr_y_yyyyz_yz, tr_y_yyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyz_xx[i] = tr_y_yyyy_xx[i] * pa_z[i];

        tr_y_yyyyz_xy[i] = tr_y_yyyy_xy[i] * pa_z[i];

        tr_y_yyyyz_xz[i] = tr_y_yyyy_x[i] * fe_0 + tr_y_yyyy_xz[i] * pa_z[i];

        tr_y_yyyyz_yy[i] = tr_y_yyyy_yy[i] * pa_z[i];

        tr_y_yyyyz_yz[i] = tr_y_yyyy_y[i] * fe_0 + tr_y_yyyy_yz[i] * pa_z[i];

        tr_y_yyyyz_zz[i] = 2.0 * tr_y_yyyy_z[i] * fe_0 + tr_y_yyyy_zz[i] * pa_z[i];
    }

    // Set up 228-234 components of targeted buffer : HD

    auto tr_y_yyyzz_xx = pbuffer.data(idx_dip_hd + 228);

    auto tr_y_yyyzz_xy = pbuffer.data(idx_dip_hd + 229);

    auto tr_y_yyyzz_xz = pbuffer.data(idx_dip_hd + 230);

    auto tr_y_yyyzz_yy = pbuffer.data(idx_dip_hd + 231);

    auto tr_y_yyyzz_yz = pbuffer.data(idx_dip_hd + 232);

    auto tr_y_yyyzz_zz = pbuffer.data(idx_dip_hd + 233);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyy_xx, tr_y_yyy_xy, tr_y_yyy_yy, tr_y_yyy_yz, tr_y_yyyz_xx, tr_y_yyyz_xy, tr_y_yyyz_y, tr_y_yyyz_yy, tr_y_yyyz_yz, tr_y_yyyzz_xx, tr_y_yyyzz_xy, tr_y_yyyzz_xz, tr_y_yyyzz_yy, tr_y_yyyzz_yz, tr_y_yyyzz_zz, tr_y_yyzz_xz, tr_y_yyzz_zz, tr_y_yzz_xz, tr_y_yzz_zz, ts_yyzz_xz, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzz_xx[i] = tr_y_yyy_xx[i] * fe_0 + tr_y_yyyz_xx[i] * pa_z[i];

        tr_y_yyyzz_xy[i] = tr_y_yyy_xy[i] * fe_0 + tr_y_yyyz_xy[i] * pa_z[i];

        tr_y_yyyzz_xz[i] = 2.0 * tr_y_yzz_xz[i] * fe_0 + ts_yyzz_xz[i] * fe_0 + tr_y_yyzz_xz[i] * pa_y[i];

        tr_y_yyyzz_yy[i] = tr_y_yyy_yy[i] * fe_0 + tr_y_yyyz_yy[i] * pa_z[i];

        tr_y_yyyzz_yz[i] = tr_y_yyy_yz[i] * fe_0 + tr_y_yyyz_y[i] * fe_0 + tr_y_yyyz_yz[i] * pa_z[i];

        tr_y_yyyzz_zz[i] = 2.0 * tr_y_yzz_zz[i] * fe_0 + ts_yyzz_zz[i] * fe_0 + tr_y_yyzz_zz[i] * pa_y[i];
    }

    // Set up 234-240 components of targeted buffer : HD

    auto tr_y_yyzzz_xx = pbuffer.data(idx_dip_hd + 234);

    auto tr_y_yyzzz_xy = pbuffer.data(idx_dip_hd + 235);

    auto tr_y_yyzzz_xz = pbuffer.data(idx_dip_hd + 236);

    auto tr_y_yyzzz_yy = pbuffer.data(idx_dip_hd + 237);

    auto tr_y_yyzzz_yz = pbuffer.data(idx_dip_hd + 238);

    auto tr_y_yyzzz_zz = pbuffer.data(idx_dip_hd + 239);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yyz_xx, tr_y_yyz_xy, tr_y_yyz_yy, tr_y_yyz_yz, tr_y_yyzz_xx, tr_y_yyzz_xy, tr_y_yyzz_y, tr_y_yyzz_yy, tr_y_yyzz_yz, tr_y_yyzzz_xx, tr_y_yyzzz_xy, tr_y_yyzzz_xz, tr_y_yyzzz_yy, tr_y_yyzzz_yz, tr_y_yyzzz_zz, tr_y_yzzz_xz, tr_y_yzzz_zz, tr_y_zzz_xz, tr_y_zzz_zz, ts_yzzz_xz, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzz_xx[i] = 2.0 * tr_y_yyz_xx[i] * fe_0 + tr_y_yyzz_xx[i] * pa_z[i];

        tr_y_yyzzz_xy[i] = 2.0 * tr_y_yyz_xy[i] * fe_0 + tr_y_yyzz_xy[i] * pa_z[i];

        tr_y_yyzzz_xz[i] = tr_y_zzz_xz[i] * fe_0 + ts_yzzz_xz[i] * fe_0 + tr_y_yzzz_xz[i] * pa_y[i];

        tr_y_yyzzz_yy[i] = 2.0 * tr_y_yyz_yy[i] * fe_0 + tr_y_yyzz_yy[i] * pa_z[i];

        tr_y_yyzzz_yz[i] = 2.0 * tr_y_yyz_yz[i] * fe_0 + tr_y_yyzz_y[i] * fe_0 + tr_y_yyzz_yz[i] * pa_z[i];

        tr_y_yyzzz_zz[i] = tr_y_zzz_zz[i] * fe_0 + ts_yzzz_zz[i] * fe_0 + tr_y_yzzz_zz[i] * pa_y[i];
    }

    // Set up 240-246 components of targeted buffer : HD

    auto tr_y_yzzzz_xx = pbuffer.data(idx_dip_hd + 240);

    auto tr_y_yzzzz_xy = pbuffer.data(idx_dip_hd + 241);

    auto tr_y_yzzzz_xz = pbuffer.data(idx_dip_hd + 242);

    auto tr_y_yzzzz_yy = pbuffer.data(idx_dip_hd + 243);

    auto tr_y_yzzzz_yz = pbuffer.data(idx_dip_hd + 244);

    auto tr_y_yzzzz_zz = pbuffer.data(idx_dip_hd + 245);

    #pragma omp simd aligned(pa_y, pa_z, tr_y_yzz_xy, tr_y_yzz_yy, tr_y_yzzz_xy, tr_y_yzzz_yy, tr_y_yzzzz_xx, tr_y_yzzzz_xy, tr_y_yzzzz_xz, tr_y_yzzzz_yy, tr_y_yzzzz_yz, tr_y_yzzzz_zz, tr_y_zzzz_xx, tr_y_zzzz_xz, tr_y_zzzz_yz, tr_y_zzzz_z, tr_y_zzzz_zz, ts_zzzz_xx, ts_zzzz_xz, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzz_xx[i] = ts_zzzz_xx[i] * fe_0 + tr_y_zzzz_xx[i] * pa_y[i];

        tr_y_yzzzz_xy[i] = 3.0 * tr_y_yzz_xy[i] * fe_0 + tr_y_yzzz_xy[i] * pa_z[i];

        tr_y_yzzzz_xz[i] = ts_zzzz_xz[i] * fe_0 + tr_y_zzzz_xz[i] * pa_y[i];

        tr_y_yzzzz_yy[i] = 3.0 * tr_y_yzz_yy[i] * fe_0 + tr_y_yzzz_yy[i] * pa_z[i];

        tr_y_yzzzz_yz[i] = tr_y_zzzz_z[i] * fe_0 + ts_zzzz_yz[i] * fe_0 + tr_y_zzzz_yz[i] * pa_y[i];

        tr_y_yzzzz_zz[i] = ts_zzzz_zz[i] * fe_0 + tr_y_zzzz_zz[i] * pa_y[i];
    }

    // Set up 246-252 components of targeted buffer : HD

    auto tr_y_zzzzz_xx = pbuffer.data(idx_dip_hd + 246);

    auto tr_y_zzzzz_xy = pbuffer.data(idx_dip_hd + 247);

    auto tr_y_zzzzz_xz = pbuffer.data(idx_dip_hd + 248);

    auto tr_y_zzzzz_yy = pbuffer.data(idx_dip_hd + 249);

    auto tr_y_zzzzz_yz = pbuffer.data(idx_dip_hd + 250);

    auto tr_y_zzzzz_zz = pbuffer.data(idx_dip_hd + 251);

    #pragma omp simd aligned(pa_z, tr_y_zzz_xx, tr_y_zzz_xy, tr_y_zzz_xz, tr_y_zzz_yy, tr_y_zzz_yz, tr_y_zzz_zz, tr_y_zzzz_x, tr_y_zzzz_xx, tr_y_zzzz_xy, tr_y_zzzz_xz, tr_y_zzzz_y, tr_y_zzzz_yy, tr_y_zzzz_yz, tr_y_zzzz_z, tr_y_zzzz_zz, tr_y_zzzzz_xx, tr_y_zzzzz_xy, tr_y_zzzzz_xz, tr_y_zzzzz_yy, tr_y_zzzzz_yz, tr_y_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzz_xx[i] = 4.0 * tr_y_zzz_xx[i] * fe_0 + tr_y_zzzz_xx[i] * pa_z[i];

        tr_y_zzzzz_xy[i] = 4.0 * tr_y_zzz_xy[i] * fe_0 + tr_y_zzzz_xy[i] * pa_z[i];

        tr_y_zzzzz_xz[i] = 4.0 * tr_y_zzz_xz[i] * fe_0 + tr_y_zzzz_x[i] * fe_0 + tr_y_zzzz_xz[i] * pa_z[i];

        tr_y_zzzzz_yy[i] = 4.0 * tr_y_zzz_yy[i] * fe_0 + tr_y_zzzz_yy[i] * pa_z[i];

        tr_y_zzzzz_yz[i] = 4.0 * tr_y_zzz_yz[i] * fe_0 + tr_y_zzzz_y[i] * fe_0 + tr_y_zzzz_yz[i] * pa_z[i];

        tr_y_zzzzz_zz[i] = 4.0 * tr_y_zzz_zz[i] * fe_0 + 2.0 * tr_y_zzzz_z[i] * fe_0 + tr_y_zzzz_zz[i] * pa_z[i];
    }

    // Set up 252-258 components of targeted buffer : HD

    auto tr_z_xxxxx_xx = pbuffer.data(idx_dip_hd + 252);

    auto tr_z_xxxxx_xy = pbuffer.data(idx_dip_hd + 253);

    auto tr_z_xxxxx_xz = pbuffer.data(idx_dip_hd + 254);

    auto tr_z_xxxxx_yy = pbuffer.data(idx_dip_hd + 255);

    auto tr_z_xxxxx_yz = pbuffer.data(idx_dip_hd + 256);

    auto tr_z_xxxxx_zz = pbuffer.data(idx_dip_hd + 257);

    #pragma omp simd aligned(pa_x, tr_z_xxx_xx, tr_z_xxx_xy, tr_z_xxx_xz, tr_z_xxx_yy, tr_z_xxx_yz, tr_z_xxx_zz, tr_z_xxxx_x, tr_z_xxxx_xx, tr_z_xxxx_xy, tr_z_xxxx_xz, tr_z_xxxx_y, tr_z_xxxx_yy, tr_z_xxxx_yz, tr_z_xxxx_z, tr_z_xxxx_zz, tr_z_xxxxx_xx, tr_z_xxxxx_xy, tr_z_xxxxx_xz, tr_z_xxxxx_yy, tr_z_xxxxx_yz, tr_z_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxx_xx[i] = 4.0 * tr_z_xxx_xx[i] * fe_0 + 2.0 * tr_z_xxxx_x[i] * fe_0 + tr_z_xxxx_xx[i] * pa_x[i];

        tr_z_xxxxx_xy[i] = 4.0 * tr_z_xxx_xy[i] * fe_0 + tr_z_xxxx_y[i] * fe_0 + tr_z_xxxx_xy[i] * pa_x[i];

        tr_z_xxxxx_xz[i] = 4.0 * tr_z_xxx_xz[i] * fe_0 + tr_z_xxxx_z[i] * fe_0 + tr_z_xxxx_xz[i] * pa_x[i];

        tr_z_xxxxx_yy[i] = 4.0 * tr_z_xxx_yy[i] * fe_0 + tr_z_xxxx_yy[i] * pa_x[i];

        tr_z_xxxxx_yz[i] = 4.0 * tr_z_xxx_yz[i] * fe_0 + tr_z_xxxx_yz[i] * pa_x[i];

        tr_z_xxxxx_zz[i] = 4.0 * tr_z_xxx_zz[i] * fe_0 + tr_z_xxxx_zz[i] * pa_x[i];
    }

    // Set up 258-264 components of targeted buffer : HD

    auto tr_z_xxxxy_xx = pbuffer.data(idx_dip_hd + 258);

    auto tr_z_xxxxy_xy = pbuffer.data(idx_dip_hd + 259);

    auto tr_z_xxxxy_xz = pbuffer.data(idx_dip_hd + 260);

    auto tr_z_xxxxy_yy = pbuffer.data(idx_dip_hd + 261);

    auto tr_z_xxxxy_yz = pbuffer.data(idx_dip_hd + 262);

    auto tr_z_xxxxy_zz = pbuffer.data(idx_dip_hd + 263);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxx_x, tr_z_xxxx_xx, tr_z_xxxx_xy, tr_z_xxxx_xz, tr_z_xxxx_zz, tr_z_xxxxy_xx, tr_z_xxxxy_xy, tr_z_xxxxy_xz, tr_z_xxxxy_yy, tr_z_xxxxy_yz, tr_z_xxxxy_zz, tr_z_xxxy_yy, tr_z_xxxy_yz, tr_z_xxy_yy, tr_z_xxy_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxy_xx[i] = tr_z_xxxx_xx[i] * pa_y[i];

        tr_z_xxxxy_xy[i] = tr_z_xxxx_x[i] * fe_0 + tr_z_xxxx_xy[i] * pa_y[i];

        tr_z_xxxxy_xz[i] = tr_z_xxxx_xz[i] * pa_y[i];

        tr_z_xxxxy_yy[i] = 3.0 * tr_z_xxy_yy[i] * fe_0 + tr_z_xxxy_yy[i] * pa_x[i];

        tr_z_xxxxy_yz[i] = 3.0 * tr_z_xxy_yz[i] * fe_0 + tr_z_xxxy_yz[i] * pa_x[i];

        tr_z_xxxxy_zz[i] = tr_z_xxxx_zz[i] * pa_y[i];
    }

    // Set up 264-270 components of targeted buffer : HD

    auto tr_z_xxxxz_xx = pbuffer.data(idx_dip_hd + 264);

    auto tr_z_xxxxz_xy = pbuffer.data(idx_dip_hd + 265);

    auto tr_z_xxxxz_xz = pbuffer.data(idx_dip_hd + 266);

    auto tr_z_xxxxz_yy = pbuffer.data(idx_dip_hd + 267);

    auto tr_z_xxxxz_yz = pbuffer.data(idx_dip_hd + 268);

    auto tr_z_xxxxz_zz = pbuffer.data(idx_dip_hd + 269);

    #pragma omp simd aligned(pa_x, pa_z, tr_z_xxxx_xx, tr_z_xxxx_xy, tr_z_xxxxz_xx, tr_z_xxxxz_xy, tr_z_xxxxz_xz, tr_z_xxxxz_yy, tr_z_xxxxz_yz, tr_z_xxxxz_zz, tr_z_xxxz_xz, tr_z_xxxz_yy, tr_z_xxxz_yz, tr_z_xxxz_z, tr_z_xxxz_zz, tr_z_xxz_xz, tr_z_xxz_yy, tr_z_xxz_yz, tr_z_xxz_zz, ts_xxxx_xx, ts_xxxx_xy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxz_xx[i] = ts_xxxx_xx[i] * fe_0 + tr_z_xxxx_xx[i] * pa_z[i];

        tr_z_xxxxz_xy[i] = ts_xxxx_xy[i] * fe_0 + tr_z_xxxx_xy[i] * pa_z[i];

        tr_z_xxxxz_xz[i] = 3.0 * tr_z_xxz_xz[i] * fe_0 + tr_z_xxxz_z[i] * fe_0 + tr_z_xxxz_xz[i] * pa_x[i];

        tr_z_xxxxz_yy[i] = 3.0 * tr_z_xxz_yy[i] * fe_0 + tr_z_xxxz_yy[i] * pa_x[i];

        tr_z_xxxxz_yz[i] = 3.0 * tr_z_xxz_yz[i] * fe_0 + tr_z_xxxz_yz[i] * pa_x[i];

        tr_z_xxxxz_zz[i] = 3.0 * tr_z_xxz_zz[i] * fe_0 + tr_z_xxxz_zz[i] * pa_x[i];
    }

    // Set up 270-276 components of targeted buffer : HD

    auto tr_z_xxxyy_xx = pbuffer.data(idx_dip_hd + 270);

    auto tr_z_xxxyy_xy = pbuffer.data(idx_dip_hd + 271);

    auto tr_z_xxxyy_xz = pbuffer.data(idx_dip_hd + 272);

    auto tr_z_xxxyy_yy = pbuffer.data(idx_dip_hd + 273);

    auto tr_z_xxxyy_yz = pbuffer.data(idx_dip_hd + 274);

    auto tr_z_xxxyy_zz = pbuffer.data(idx_dip_hd + 275);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxx_xx, tr_z_xxx_xz, tr_z_xxxy_xx, tr_z_xxxy_xz, tr_z_xxxyy_xx, tr_z_xxxyy_xy, tr_z_xxxyy_xz, tr_z_xxxyy_yy, tr_z_xxxyy_yz, tr_z_xxxyy_zz, tr_z_xxyy_xy, tr_z_xxyy_y, tr_z_xxyy_yy, tr_z_xxyy_yz, tr_z_xxyy_zz, tr_z_xyy_xy, tr_z_xyy_yy, tr_z_xyy_yz, tr_z_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyy_xx[i] = tr_z_xxx_xx[i] * fe_0 + tr_z_xxxy_xx[i] * pa_y[i];

        tr_z_xxxyy_xy[i] = 2.0 * tr_z_xyy_xy[i] * fe_0 + tr_z_xxyy_y[i] * fe_0 + tr_z_xxyy_xy[i] * pa_x[i];

        tr_z_xxxyy_xz[i] = tr_z_xxx_xz[i] * fe_0 + tr_z_xxxy_xz[i] * pa_y[i];

        tr_z_xxxyy_yy[i] = 2.0 * tr_z_xyy_yy[i] * fe_0 + tr_z_xxyy_yy[i] * pa_x[i];

        tr_z_xxxyy_yz[i] = 2.0 * tr_z_xyy_yz[i] * fe_0 + tr_z_xxyy_yz[i] * pa_x[i];

        tr_z_xxxyy_zz[i] = 2.0 * tr_z_xyy_zz[i] * fe_0 + tr_z_xxyy_zz[i] * pa_x[i];
    }

    // Set up 276-282 components of targeted buffer : HD

    auto tr_z_xxxyz_xx = pbuffer.data(idx_dip_hd + 276);

    auto tr_z_xxxyz_xy = pbuffer.data(idx_dip_hd + 277);

    auto tr_z_xxxyz_xz = pbuffer.data(idx_dip_hd + 278);

    auto tr_z_xxxyz_yy = pbuffer.data(idx_dip_hd + 279);

    auto tr_z_xxxyz_yz = pbuffer.data(idx_dip_hd + 280);

    auto tr_z_xxxyz_zz = pbuffer.data(idx_dip_hd + 281);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxxyz_xx, tr_z_xxxyz_xy, tr_z_xxxyz_xz, tr_z_xxxyz_yy, tr_z_xxxyz_yz, tr_z_xxxyz_zz, tr_z_xxxz_x, tr_z_xxxz_xx, tr_z_xxxz_xy, tr_z_xxxz_xz, tr_z_xxxz_zz, tr_z_xxyz_yy, tr_z_xxyz_yz, tr_z_xyz_yy, tr_z_xyz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyz_xx[i] = tr_z_xxxz_xx[i] * pa_y[i];

        tr_z_xxxyz_xy[i] = tr_z_xxxz_x[i] * fe_0 + tr_z_xxxz_xy[i] * pa_y[i];

        tr_z_xxxyz_xz[i] = tr_z_xxxz_xz[i] * pa_y[i];

        tr_z_xxxyz_yy[i] = 2.0 * tr_z_xyz_yy[i] * fe_0 + tr_z_xxyz_yy[i] * pa_x[i];

        tr_z_xxxyz_yz[i] = 2.0 * tr_z_xyz_yz[i] * fe_0 + tr_z_xxyz_yz[i] * pa_x[i];

        tr_z_xxxyz_zz[i] = tr_z_xxxz_zz[i] * pa_y[i];
    }

    // Set up 282-288 components of targeted buffer : HD

    auto tr_z_xxxzz_xx = pbuffer.data(idx_dip_hd + 282);

    auto tr_z_xxxzz_xy = pbuffer.data(idx_dip_hd + 283);

    auto tr_z_xxxzz_xz = pbuffer.data(idx_dip_hd + 284);

    auto tr_z_xxxzz_yy = pbuffer.data(idx_dip_hd + 285);

    auto tr_z_xxxzz_yz = pbuffer.data(idx_dip_hd + 286);

    auto tr_z_xxxzz_zz = pbuffer.data(idx_dip_hd + 287);

    #pragma omp simd aligned(pa_x, tr_z_xxxzz_xx, tr_z_xxxzz_xy, tr_z_xxxzz_xz, tr_z_xxxzz_yy, tr_z_xxxzz_yz, tr_z_xxxzz_zz, tr_z_xxzz_x, tr_z_xxzz_xx, tr_z_xxzz_xy, tr_z_xxzz_xz, tr_z_xxzz_y, tr_z_xxzz_yy, tr_z_xxzz_yz, tr_z_xxzz_z, tr_z_xxzz_zz, tr_z_xzz_xx, tr_z_xzz_xy, tr_z_xzz_xz, tr_z_xzz_yy, tr_z_xzz_yz, tr_z_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzz_xx[i] = 2.0 * tr_z_xzz_xx[i] * fe_0 + 2.0 * tr_z_xxzz_x[i] * fe_0 + tr_z_xxzz_xx[i] * pa_x[i];

        tr_z_xxxzz_xy[i] = 2.0 * tr_z_xzz_xy[i] * fe_0 + tr_z_xxzz_y[i] * fe_0 + tr_z_xxzz_xy[i] * pa_x[i];

        tr_z_xxxzz_xz[i] = 2.0 * tr_z_xzz_xz[i] * fe_0 + tr_z_xxzz_z[i] * fe_0 + tr_z_xxzz_xz[i] * pa_x[i];

        tr_z_xxxzz_yy[i] = 2.0 * tr_z_xzz_yy[i] * fe_0 + tr_z_xxzz_yy[i] * pa_x[i];

        tr_z_xxxzz_yz[i] = 2.0 * tr_z_xzz_yz[i] * fe_0 + tr_z_xxzz_yz[i] * pa_x[i];

        tr_z_xxxzz_zz[i] = 2.0 * tr_z_xzz_zz[i] * fe_0 + tr_z_xxzz_zz[i] * pa_x[i];
    }

    // Set up 288-294 components of targeted buffer : HD

    auto tr_z_xxyyy_xx = pbuffer.data(idx_dip_hd + 288);

    auto tr_z_xxyyy_xy = pbuffer.data(idx_dip_hd + 289);

    auto tr_z_xxyyy_xz = pbuffer.data(idx_dip_hd + 290);

    auto tr_z_xxyyy_yy = pbuffer.data(idx_dip_hd + 291);

    auto tr_z_xxyyy_yz = pbuffer.data(idx_dip_hd + 292);

    auto tr_z_xxyyy_zz = pbuffer.data(idx_dip_hd + 293);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxy_xx, tr_z_xxy_xz, tr_z_xxyy_xx, tr_z_xxyy_xz, tr_z_xxyyy_xx, tr_z_xxyyy_xy, tr_z_xxyyy_xz, tr_z_xxyyy_yy, tr_z_xxyyy_yz, tr_z_xxyyy_zz, tr_z_xyyy_xy, tr_z_xyyy_y, tr_z_xyyy_yy, tr_z_xyyy_yz, tr_z_xyyy_zz, tr_z_yyy_xy, tr_z_yyy_yy, tr_z_yyy_yz, tr_z_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyy_xx[i] = 2.0 * tr_z_xxy_xx[i] * fe_0 + tr_z_xxyy_xx[i] * pa_y[i];

        tr_z_xxyyy_xy[i] = tr_z_yyy_xy[i] * fe_0 + tr_z_xyyy_y[i] * fe_0 + tr_z_xyyy_xy[i] * pa_x[i];

        tr_z_xxyyy_xz[i] = 2.0 * tr_z_xxy_xz[i] * fe_0 + tr_z_xxyy_xz[i] * pa_y[i];

        tr_z_xxyyy_yy[i] = tr_z_yyy_yy[i] * fe_0 + tr_z_xyyy_yy[i] * pa_x[i];

        tr_z_xxyyy_yz[i] = tr_z_yyy_yz[i] * fe_0 + tr_z_xyyy_yz[i] * pa_x[i];

        tr_z_xxyyy_zz[i] = tr_z_yyy_zz[i] * fe_0 + tr_z_xyyy_zz[i] * pa_x[i];
    }

    // Set up 294-300 components of targeted buffer : HD

    auto tr_z_xxyyz_xx = pbuffer.data(idx_dip_hd + 294);

    auto tr_z_xxyyz_xy = pbuffer.data(idx_dip_hd + 295);

    auto tr_z_xxyyz_xz = pbuffer.data(idx_dip_hd + 296);

    auto tr_z_xxyyz_yy = pbuffer.data(idx_dip_hd + 297);

    auto tr_z_xxyyz_yz = pbuffer.data(idx_dip_hd + 298);

    auto tr_z_xxyyz_zz = pbuffer.data(idx_dip_hd + 299);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_z_xxyy_xy, tr_z_xxyyz_xx, tr_z_xxyyz_xy, tr_z_xxyyz_xz, tr_z_xxyyz_yy, tr_z_xxyyz_yz, tr_z_xxyyz_zz, tr_z_xxyz_xx, tr_z_xxyz_xz, tr_z_xxz_xx, tr_z_xxz_xz, tr_z_xyyz_yy, tr_z_xyyz_yz, tr_z_xyyz_zz, tr_z_yyz_yy, tr_z_yyz_yz, tr_z_yyz_zz, ts_xxyy_xy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyz_xx[i] = tr_z_xxz_xx[i] * fe_0 + tr_z_xxyz_xx[i] * pa_y[i];

        tr_z_xxyyz_xy[i] = ts_xxyy_xy[i] * fe_0 + tr_z_xxyy_xy[i] * pa_z[i];

        tr_z_xxyyz_xz[i] = tr_z_xxz_xz[i] * fe_0 + tr_z_xxyz_xz[i] * pa_y[i];

        tr_z_xxyyz_yy[i] = tr_z_yyz_yy[i] * fe_0 + tr_z_xyyz_yy[i] * pa_x[i];

        tr_z_xxyyz_yz[i] = tr_z_yyz_yz[i] * fe_0 + tr_z_xyyz_yz[i] * pa_x[i];

        tr_z_xxyyz_zz[i] = tr_z_yyz_zz[i] * fe_0 + tr_z_xyyz_zz[i] * pa_x[i];
    }

    // Set up 300-306 components of targeted buffer : HD

    auto tr_z_xxyzz_xx = pbuffer.data(idx_dip_hd + 300);

    auto tr_z_xxyzz_xy = pbuffer.data(idx_dip_hd + 301);

    auto tr_z_xxyzz_xz = pbuffer.data(idx_dip_hd + 302);

    auto tr_z_xxyzz_yy = pbuffer.data(idx_dip_hd + 303);

    auto tr_z_xxyzz_yz = pbuffer.data(idx_dip_hd + 304);

    auto tr_z_xxyzz_zz = pbuffer.data(idx_dip_hd + 305);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xxyzz_xx, tr_z_xxyzz_xy, tr_z_xxyzz_xz, tr_z_xxyzz_yy, tr_z_xxyzz_yz, tr_z_xxyzz_zz, tr_z_xxzz_x, tr_z_xxzz_xx, tr_z_xxzz_xy, tr_z_xxzz_xz, tr_z_xxzz_zz, tr_z_xyzz_yy, tr_z_xyzz_yz, tr_z_yzz_yy, tr_z_yzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzz_xx[i] = tr_z_xxzz_xx[i] * pa_y[i];

        tr_z_xxyzz_xy[i] = tr_z_xxzz_x[i] * fe_0 + tr_z_xxzz_xy[i] * pa_y[i];

        tr_z_xxyzz_xz[i] = tr_z_xxzz_xz[i] * pa_y[i];

        tr_z_xxyzz_yy[i] = tr_z_yzz_yy[i] * fe_0 + tr_z_xyzz_yy[i] * pa_x[i];

        tr_z_xxyzz_yz[i] = tr_z_yzz_yz[i] * fe_0 + tr_z_xyzz_yz[i] * pa_x[i];

        tr_z_xxyzz_zz[i] = tr_z_xxzz_zz[i] * pa_y[i];
    }

    // Set up 306-312 components of targeted buffer : HD

    auto tr_z_xxzzz_xx = pbuffer.data(idx_dip_hd + 306);

    auto tr_z_xxzzz_xy = pbuffer.data(idx_dip_hd + 307);

    auto tr_z_xxzzz_xz = pbuffer.data(idx_dip_hd + 308);

    auto tr_z_xxzzz_yy = pbuffer.data(idx_dip_hd + 309);

    auto tr_z_xxzzz_yz = pbuffer.data(idx_dip_hd + 310);

    auto tr_z_xxzzz_zz = pbuffer.data(idx_dip_hd + 311);

    #pragma omp simd aligned(pa_x, tr_z_xxzzz_xx, tr_z_xxzzz_xy, tr_z_xxzzz_xz, tr_z_xxzzz_yy, tr_z_xxzzz_yz, tr_z_xxzzz_zz, tr_z_xzzz_x, tr_z_xzzz_xx, tr_z_xzzz_xy, tr_z_xzzz_xz, tr_z_xzzz_y, tr_z_xzzz_yy, tr_z_xzzz_yz, tr_z_xzzz_z, tr_z_xzzz_zz, tr_z_zzz_xx, tr_z_zzz_xy, tr_z_zzz_xz, tr_z_zzz_yy, tr_z_zzz_yz, tr_z_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzz_xx[i] = tr_z_zzz_xx[i] * fe_0 + 2.0 * tr_z_xzzz_x[i] * fe_0 + tr_z_xzzz_xx[i] * pa_x[i];

        tr_z_xxzzz_xy[i] = tr_z_zzz_xy[i] * fe_0 + tr_z_xzzz_y[i] * fe_0 + tr_z_xzzz_xy[i] * pa_x[i];

        tr_z_xxzzz_xz[i] = tr_z_zzz_xz[i] * fe_0 + tr_z_xzzz_z[i] * fe_0 + tr_z_xzzz_xz[i] * pa_x[i];

        tr_z_xxzzz_yy[i] = tr_z_zzz_yy[i] * fe_0 + tr_z_xzzz_yy[i] * pa_x[i];

        tr_z_xxzzz_yz[i] = tr_z_zzz_yz[i] * fe_0 + tr_z_xzzz_yz[i] * pa_x[i];

        tr_z_xxzzz_zz[i] = tr_z_zzz_zz[i] * fe_0 + tr_z_xzzz_zz[i] * pa_x[i];
    }

    // Set up 312-318 components of targeted buffer : HD

    auto tr_z_xyyyy_xx = pbuffer.data(idx_dip_hd + 312);

    auto tr_z_xyyyy_xy = pbuffer.data(idx_dip_hd + 313);

    auto tr_z_xyyyy_xz = pbuffer.data(idx_dip_hd + 314);

    auto tr_z_xyyyy_yy = pbuffer.data(idx_dip_hd + 315);

    auto tr_z_xyyyy_yz = pbuffer.data(idx_dip_hd + 316);

    auto tr_z_xyyyy_zz = pbuffer.data(idx_dip_hd + 317);

    #pragma omp simd aligned(pa_x, tr_z_xyyyy_xx, tr_z_xyyyy_xy, tr_z_xyyyy_xz, tr_z_xyyyy_yy, tr_z_xyyyy_yz, tr_z_xyyyy_zz, tr_z_yyyy_x, tr_z_yyyy_xx, tr_z_yyyy_xy, tr_z_yyyy_xz, tr_z_yyyy_y, tr_z_yyyy_yy, tr_z_yyyy_yz, tr_z_yyyy_z, tr_z_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyy_xx[i] = 2.0 * tr_z_yyyy_x[i] * fe_0 + tr_z_yyyy_xx[i] * pa_x[i];

        tr_z_xyyyy_xy[i] = tr_z_yyyy_y[i] * fe_0 + tr_z_yyyy_xy[i] * pa_x[i];

        tr_z_xyyyy_xz[i] = tr_z_yyyy_z[i] * fe_0 + tr_z_yyyy_xz[i] * pa_x[i];

        tr_z_xyyyy_yy[i] = tr_z_yyyy_yy[i] * pa_x[i];

        tr_z_xyyyy_yz[i] = tr_z_yyyy_yz[i] * pa_x[i];

        tr_z_xyyyy_zz[i] = tr_z_yyyy_zz[i] * pa_x[i];
    }

    // Set up 318-324 components of targeted buffer : HD

    auto tr_z_xyyyz_xx = pbuffer.data(idx_dip_hd + 318);

    auto tr_z_xyyyz_xy = pbuffer.data(idx_dip_hd + 319);

    auto tr_z_xyyyz_xz = pbuffer.data(idx_dip_hd + 320);

    auto tr_z_xyyyz_yy = pbuffer.data(idx_dip_hd + 321);

    auto tr_z_xyyyz_yz = pbuffer.data(idx_dip_hd + 322);

    auto tr_z_xyyyz_zz = pbuffer.data(idx_dip_hd + 323);

    #pragma omp simd aligned(pa_x, tr_z_xyyyz_xx, tr_z_xyyyz_xy, tr_z_xyyyz_xz, tr_z_xyyyz_yy, tr_z_xyyyz_yz, tr_z_xyyyz_zz, tr_z_yyyz_x, tr_z_yyyz_xx, tr_z_yyyz_xy, tr_z_yyyz_xz, tr_z_yyyz_y, tr_z_yyyz_yy, tr_z_yyyz_yz, tr_z_yyyz_z, tr_z_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyz_xx[i] = 2.0 * tr_z_yyyz_x[i] * fe_0 + tr_z_yyyz_xx[i] * pa_x[i];

        tr_z_xyyyz_xy[i] = tr_z_yyyz_y[i] * fe_0 + tr_z_yyyz_xy[i] * pa_x[i];

        tr_z_xyyyz_xz[i] = tr_z_yyyz_z[i] * fe_0 + tr_z_yyyz_xz[i] * pa_x[i];

        tr_z_xyyyz_yy[i] = tr_z_yyyz_yy[i] * pa_x[i];

        tr_z_xyyyz_yz[i] = tr_z_yyyz_yz[i] * pa_x[i];

        tr_z_xyyyz_zz[i] = tr_z_yyyz_zz[i] * pa_x[i];
    }

    // Set up 324-330 components of targeted buffer : HD

    auto tr_z_xyyzz_xx = pbuffer.data(idx_dip_hd + 324);

    auto tr_z_xyyzz_xy = pbuffer.data(idx_dip_hd + 325);

    auto tr_z_xyyzz_xz = pbuffer.data(idx_dip_hd + 326);

    auto tr_z_xyyzz_yy = pbuffer.data(idx_dip_hd + 327);

    auto tr_z_xyyzz_yz = pbuffer.data(idx_dip_hd + 328);

    auto tr_z_xyyzz_zz = pbuffer.data(idx_dip_hd + 329);

    #pragma omp simd aligned(pa_x, tr_z_xyyzz_xx, tr_z_xyyzz_xy, tr_z_xyyzz_xz, tr_z_xyyzz_yy, tr_z_xyyzz_yz, tr_z_xyyzz_zz, tr_z_yyzz_x, tr_z_yyzz_xx, tr_z_yyzz_xy, tr_z_yyzz_xz, tr_z_yyzz_y, tr_z_yyzz_yy, tr_z_yyzz_yz, tr_z_yyzz_z, tr_z_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzz_xx[i] = 2.0 * tr_z_yyzz_x[i] * fe_0 + tr_z_yyzz_xx[i] * pa_x[i];

        tr_z_xyyzz_xy[i] = tr_z_yyzz_y[i] * fe_0 + tr_z_yyzz_xy[i] * pa_x[i];

        tr_z_xyyzz_xz[i] = tr_z_yyzz_z[i] * fe_0 + tr_z_yyzz_xz[i] * pa_x[i];

        tr_z_xyyzz_yy[i] = tr_z_yyzz_yy[i] * pa_x[i];

        tr_z_xyyzz_yz[i] = tr_z_yyzz_yz[i] * pa_x[i];

        tr_z_xyyzz_zz[i] = tr_z_yyzz_zz[i] * pa_x[i];
    }

    // Set up 330-336 components of targeted buffer : HD

    auto tr_z_xyzzz_xx = pbuffer.data(idx_dip_hd + 330);

    auto tr_z_xyzzz_xy = pbuffer.data(idx_dip_hd + 331);

    auto tr_z_xyzzz_xz = pbuffer.data(idx_dip_hd + 332);

    auto tr_z_xyzzz_yy = pbuffer.data(idx_dip_hd + 333);

    auto tr_z_xyzzz_yz = pbuffer.data(idx_dip_hd + 334);

    auto tr_z_xyzzz_zz = pbuffer.data(idx_dip_hd + 335);

    #pragma omp simd aligned(pa_x, pa_y, tr_z_xyzzz_xx, tr_z_xyzzz_xy, tr_z_xyzzz_xz, tr_z_xyzzz_yy, tr_z_xyzzz_yz, tr_z_xyzzz_zz, tr_z_xzzz_xx, tr_z_xzzz_xz, tr_z_yzzz_xy, tr_z_yzzz_y, tr_z_yzzz_yy, tr_z_yzzz_yz, tr_z_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzz_xx[i] = tr_z_xzzz_xx[i] * pa_y[i];

        tr_z_xyzzz_xy[i] = tr_z_yzzz_y[i] * fe_0 + tr_z_yzzz_xy[i] * pa_x[i];

        tr_z_xyzzz_xz[i] = tr_z_xzzz_xz[i] * pa_y[i];

        tr_z_xyzzz_yy[i] = tr_z_yzzz_yy[i] * pa_x[i];

        tr_z_xyzzz_yz[i] = tr_z_yzzz_yz[i] * pa_x[i];

        tr_z_xyzzz_zz[i] = tr_z_yzzz_zz[i] * pa_x[i];
    }

    // Set up 336-342 components of targeted buffer : HD

    auto tr_z_xzzzz_xx = pbuffer.data(idx_dip_hd + 336);

    auto tr_z_xzzzz_xy = pbuffer.data(idx_dip_hd + 337);

    auto tr_z_xzzzz_xz = pbuffer.data(idx_dip_hd + 338);

    auto tr_z_xzzzz_yy = pbuffer.data(idx_dip_hd + 339);

    auto tr_z_xzzzz_yz = pbuffer.data(idx_dip_hd + 340);

    auto tr_z_xzzzz_zz = pbuffer.data(idx_dip_hd + 341);

    #pragma omp simd aligned(pa_x, tr_z_xzzzz_xx, tr_z_xzzzz_xy, tr_z_xzzzz_xz, tr_z_xzzzz_yy, tr_z_xzzzz_yz, tr_z_xzzzz_zz, tr_z_zzzz_x, tr_z_zzzz_xx, tr_z_zzzz_xy, tr_z_zzzz_xz, tr_z_zzzz_y, tr_z_zzzz_yy, tr_z_zzzz_yz, tr_z_zzzz_z, tr_z_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzz_xx[i] = 2.0 * tr_z_zzzz_x[i] * fe_0 + tr_z_zzzz_xx[i] * pa_x[i];

        tr_z_xzzzz_xy[i] = tr_z_zzzz_y[i] * fe_0 + tr_z_zzzz_xy[i] * pa_x[i];

        tr_z_xzzzz_xz[i] = tr_z_zzzz_z[i] * fe_0 + tr_z_zzzz_xz[i] * pa_x[i];

        tr_z_xzzzz_yy[i] = tr_z_zzzz_yy[i] * pa_x[i];

        tr_z_xzzzz_yz[i] = tr_z_zzzz_yz[i] * pa_x[i];

        tr_z_xzzzz_zz[i] = tr_z_zzzz_zz[i] * pa_x[i];
    }

    // Set up 342-348 components of targeted buffer : HD

    auto tr_z_yyyyy_xx = pbuffer.data(idx_dip_hd + 342);

    auto tr_z_yyyyy_xy = pbuffer.data(idx_dip_hd + 343);

    auto tr_z_yyyyy_xz = pbuffer.data(idx_dip_hd + 344);

    auto tr_z_yyyyy_yy = pbuffer.data(idx_dip_hd + 345);

    auto tr_z_yyyyy_yz = pbuffer.data(idx_dip_hd + 346);

    auto tr_z_yyyyy_zz = pbuffer.data(idx_dip_hd + 347);

    #pragma omp simd aligned(pa_y, tr_z_yyy_xx, tr_z_yyy_xy, tr_z_yyy_xz, tr_z_yyy_yy, tr_z_yyy_yz, tr_z_yyy_zz, tr_z_yyyy_x, tr_z_yyyy_xx, tr_z_yyyy_xy, tr_z_yyyy_xz, tr_z_yyyy_y, tr_z_yyyy_yy, tr_z_yyyy_yz, tr_z_yyyy_z, tr_z_yyyy_zz, tr_z_yyyyy_xx, tr_z_yyyyy_xy, tr_z_yyyyy_xz, tr_z_yyyyy_yy, tr_z_yyyyy_yz, tr_z_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyy_xx[i] = 4.0 * tr_z_yyy_xx[i] * fe_0 + tr_z_yyyy_xx[i] * pa_y[i];

        tr_z_yyyyy_xy[i] = 4.0 * tr_z_yyy_xy[i] * fe_0 + tr_z_yyyy_x[i] * fe_0 + tr_z_yyyy_xy[i] * pa_y[i];

        tr_z_yyyyy_xz[i] = 4.0 * tr_z_yyy_xz[i] * fe_0 + tr_z_yyyy_xz[i] * pa_y[i];

        tr_z_yyyyy_yy[i] = 4.0 * tr_z_yyy_yy[i] * fe_0 + 2.0 * tr_z_yyyy_y[i] * fe_0 + tr_z_yyyy_yy[i] * pa_y[i];

        tr_z_yyyyy_yz[i] = 4.0 * tr_z_yyy_yz[i] * fe_0 + tr_z_yyyy_z[i] * fe_0 + tr_z_yyyy_yz[i] * pa_y[i];

        tr_z_yyyyy_zz[i] = 4.0 * tr_z_yyy_zz[i] * fe_0 + tr_z_yyyy_zz[i] * pa_y[i];
    }

    // Set up 348-354 components of targeted buffer : HD

    auto tr_z_yyyyz_xx = pbuffer.data(idx_dip_hd + 348);

    auto tr_z_yyyyz_xy = pbuffer.data(idx_dip_hd + 349);

    auto tr_z_yyyyz_xz = pbuffer.data(idx_dip_hd + 350);

    auto tr_z_yyyyz_yy = pbuffer.data(idx_dip_hd + 351);

    auto tr_z_yyyyz_yz = pbuffer.data(idx_dip_hd + 352);

    auto tr_z_yyyyz_zz = pbuffer.data(idx_dip_hd + 353);

    #pragma omp simd aligned(pa_y, pa_z, tr_z_yyyy_xy, tr_z_yyyy_yy, tr_z_yyyyz_xx, tr_z_yyyyz_xy, tr_z_yyyyz_xz, tr_z_yyyyz_yy, tr_z_yyyyz_yz, tr_z_yyyyz_zz, tr_z_yyyz_xx, tr_z_yyyz_xz, tr_z_yyyz_yz, tr_z_yyyz_z, tr_z_yyyz_zz, tr_z_yyz_xx, tr_z_yyz_xz, tr_z_yyz_yz, tr_z_yyz_zz, ts_yyyy_xy, ts_yyyy_yy, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyz_xx[i] = 3.0 * tr_z_yyz_xx[i] * fe_0 + tr_z_yyyz_xx[i] * pa_y[i];

        tr_z_yyyyz_xy[i] = ts_yyyy_xy[i] * fe_0 + tr_z_yyyy_xy[i] * pa_z[i];

        tr_z_yyyyz_xz[i] = 3.0 * tr_z_yyz_xz[i] * fe_0 + tr_z_yyyz_xz[i] * pa_y[i];

        tr_z_yyyyz_yy[i] = ts_yyyy_yy[i] * fe_0 + tr_z_yyyy_yy[i] * pa_z[i];

        tr_z_yyyyz_yz[i] = 3.0 * tr_z_yyz_yz[i] * fe_0 + tr_z_yyyz_z[i] * fe_0 + tr_z_yyyz_yz[i] * pa_y[i];

        tr_z_yyyyz_zz[i] = 3.0 * tr_z_yyz_zz[i] * fe_0 + tr_z_yyyz_zz[i] * pa_y[i];
    }

    // Set up 354-360 components of targeted buffer : HD

    auto tr_z_yyyzz_xx = pbuffer.data(idx_dip_hd + 354);

    auto tr_z_yyyzz_xy = pbuffer.data(idx_dip_hd + 355);

    auto tr_z_yyyzz_xz = pbuffer.data(idx_dip_hd + 356);

    auto tr_z_yyyzz_yy = pbuffer.data(idx_dip_hd + 357);

    auto tr_z_yyyzz_yz = pbuffer.data(idx_dip_hd + 358);

    auto tr_z_yyyzz_zz = pbuffer.data(idx_dip_hd + 359);

    #pragma omp simd aligned(pa_y, tr_z_yyyzz_xx, tr_z_yyyzz_xy, tr_z_yyyzz_xz, tr_z_yyyzz_yy, tr_z_yyyzz_yz, tr_z_yyyzz_zz, tr_z_yyzz_x, tr_z_yyzz_xx, tr_z_yyzz_xy, tr_z_yyzz_xz, tr_z_yyzz_y, tr_z_yyzz_yy, tr_z_yyzz_yz, tr_z_yyzz_z, tr_z_yyzz_zz, tr_z_yzz_xx, tr_z_yzz_xy, tr_z_yzz_xz, tr_z_yzz_yy, tr_z_yzz_yz, tr_z_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzz_xx[i] = 2.0 * tr_z_yzz_xx[i] * fe_0 + tr_z_yyzz_xx[i] * pa_y[i];

        tr_z_yyyzz_xy[i] = 2.0 * tr_z_yzz_xy[i] * fe_0 + tr_z_yyzz_x[i] * fe_0 + tr_z_yyzz_xy[i] * pa_y[i];

        tr_z_yyyzz_xz[i] = 2.0 * tr_z_yzz_xz[i] * fe_0 + tr_z_yyzz_xz[i] * pa_y[i];

        tr_z_yyyzz_yy[i] = 2.0 * tr_z_yzz_yy[i] * fe_0 + 2.0 * tr_z_yyzz_y[i] * fe_0 + tr_z_yyzz_yy[i] * pa_y[i];

        tr_z_yyyzz_yz[i] = 2.0 * tr_z_yzz_yz[i] * fe_0 + tr_z_yyzz_z[i] * fe_0 + tr_z_yyzz_yz[i] * pa_y[i];

        tr_z_yyyzz_zz[i] = 2.0 * tr_z_yzz_zz[i] * fe_0 + tr_z_yyzz_zz[i] * pa_y[i];
    }

    // Set up 360-366 components of targeted buffer : HD

    auto tr_z_yyzzz_xx = pbuffer.data(idx_dip_hd + 360);

    auto tr_z_yyzzz_xy = pbuffer.data(idx_dip_hd + 361);

    auto tr_z_yyzzz_xz = pbuffer.data(idx_dip_hd + 362);

    auto tr_z_yyzzz_yy = pbuffer.data(idx_dip_hd + 363);

    auto tr_z_yyzzz_yz = pbuffer.data(idx_dip_hd + 364);

    auto tr_z_yyzzz_zz = pbuffer.data(idx_dip_hd + 365);

    #pragma omp simd aligned(pa_y, tr_z_yyzzz_xx, tr_z_yyzzz_xy, tr_z_yyzzz_xz, tr_z_yyzzz_yy, tr_z_yyzzz_yz, tr_z_yyzzz_zz, tr_z_yzzz_x, tr_z_yzzz_xx, tr_z_yzzz_xy, tr_z_yzzz_xz, tr_z_yzzz_y, tr_z_yzzz_yy, tr_z_yzzz_yz, tr_z_yzzz_z, tr_z_yzzz_zz, tr_z_zzz_xx, tr_z_zzz_xy, tr_z_zzz_xz, tr_z_zzz_yy, tr_z_zzz_yz, tr_z_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzz_xx[i] = tr_z_zzz_xx[i] * fe_0 + tr_z_yzzz_xx[i] * pa_y[i];

        tr_z_yyzzz_xy[i] = tr_z_zzz_xy[i] * fe_0 + tr_z_yzzz_x[i] * fe_0 + tr_z_yzzz_xy[i] * pa_y[i];

        tr_z_yyzzz_xz[i] = tr_z_zzz_xz[i] * fe_0 + tr_z_yzzz_xz[i] * pa_y[i];

        tr_z_yyzzz_yy[i] = tr_z_zzz_yy[i] * fe_0 + 2.0 * tr_z_yzzz_y[i] * fe_0 + tr_z_yzzz_yy[i] * pa_y[i];

        tr_z_yyzzz_yz[i] = tr_z_zzz_yz[i] * fe_0 + tr_z_yzzz_z[i] * fe_0 + tr_z_yzzz_yz[i] * pa_y[i];

        tr_z_yyzzz_zz[i] = tr_z_zzz_zz[i] * fe_0 + tr_z_yzzz_zz[i] * pa_y[i];
    }

    // Set up 366-372 components of targeted buffer : HD

    auto tr_z_yzzzz_xx = pbuffer.data(idx_dip_hd + 366);

    auto tr_z_yzzzz_xy = pbuffer.data(idx_dip_hd + 367);

    auto tr_z_yzzzz_xz = pbuffer.data(idx_dip_hd + 368);

    auto tr_z_yzzzz_yy = pbuffer.data(idx_dip_hd + 369);

    auto tr_z_yzzzz_yz = pbuffer.data(idx_dip_hd + 370);

    auto tr_z_yzzzz_zz = pbuffer.data(idx_dip_hd + 371);

    #pragma omp simd aligned(pa_y, tr_z_yzzzz_xx, tr_z_yzzzz_xy, tr_z_yzzzz_xz, tr_z_yzzzz_yy, tr_z_yzzzz_yz, tr_z_yzzzz_zz, tr_z_zzzz_x, tr_z_zzzz_xx, tr_z_zzzz_xy, tr_z_zzzz_xz, tr_z_zzzz_y, tr_z_zzzz_yy, tr_z_zzzz_yz, tr_z_zzzz_z, tr_z_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzz_xx[i] = tr_z_zzzz_xx[i] * pa_y[i];

        tr_z_yzzzz_xy[i] = tr_z_zzzz_x[i] * fe_0 + tr_z_zzzz_xy[i] * pa_y[i];

        tr_z_yzzzz_xz[i] = tr_z_zzzz_xz[i] * pa_y[i];

        tr_z_yzzzz_yy[i] = 2.0 * tr_z_zzzz_y[i] * fe_0 + tr_z_zzzz_yy[i] * pa_y[i];

        tr_z_yzzzz_yz[i] = tr_z_zzzz_z[i] * fe_0 + tr_z_zzzz_yz[i] * pa_y[i];

        tr_z_yzzzz_zz[i] = tr_z_zzzz_zz[i] * pa_y[i];
    }

    // Set up 372-378 components of targeted buffer : HD

    auto tr_z_zzzzz_xx = pbuffer.data(idx_dip_hd + 372);

    auto tr_z_zzzzz_xy = pbuffer.data(idx_dip_hd + 373);

    auto tr_z_zzzzz_xz = pbuffer.data(idx_dip_hd + 374);

    auto tr_z_zzzzz_yy = pbuffer.data(idx_dip_hd + 375);

    auto tr_z_zzzzz_yz = pbuffer.data(idx_dip_hd + 376);

    auto tr_z_zzzzz_zz = pbuffer.data(idx_dip_hd + 377);

    #pragma omp simd aligned(pa_z, tr_z_zzz_xx, tr_z_zzz_xy, tr_z_zzz_xz, tr_z_zzz_yy, tr_z_zzz_yz, tr_z_zzz_zz, tr_z_zzzz_x, tr_z_zzzz_xx, tr_z_zzzz_xy, tr_z_zzzz_xz, tr_z_zzzz_y, tr_z_zzzz_yy, tr_z_zzzz_yz, tr_z_zzzz_z, tr_z_zzzz_zz, tr_z_zzzzz_xx, tr_z_zzzzz_xy, tr_z_zzzzz_xz, tr_z_zzzzz_yy, tr_z_zzzzz_yz, tr_z_zzzzz_zz, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzz_xx[i] = 4.0 * tr_z_zzz_xx[i] * fe_0 + ts_zzzz_xx[i] * fe_0 + tr_z_zzzz_xx[i] * pa_z[i];

        tr_z_zzzzz_xy[i] = 4.0 * tr_z_zzz_xy[i] * fe_0 + ts_zzzz_xy[i] * fe_0 + tr_z_zzzz_xy[i] * pa_z[i];

        tr_z_zzzzz_xz[i] = 4.0 * tr_z_zzz_xz[i] * fe_0 + tr_z_zzzz_x[i] * fe_0 + ts_zzzz_xz[i] * fe_0 + tr_z_zzzz_xz[i] * pa_z[i];

        tr_z_zzzzz_yy[i] = 4.0 * tr_z_zzz_yy[i] * fe_0 + ts_zzzz_yy[i] * fe_0 + tr_z_zzzz_yy[i] * pa_z[i];

        tr_z_zzzzz_yz[i] = 4.0 * tr_z_zzz_yz[i] * fe_0 + tr_z_zzzz_y[i] * fe_0 + ts_zzzz_yz[i] * fe_0 + tr_z_zzzz_yz[i] * pa_z[i];

        tr_z_zzzzz_zz[i] = 4.0 * tr_z_zzz_zz[i] * fe_0 + 2.0 * tr_z_zzzz_z[i] * fe_0 + ts_zzzz_zz[i] * fe_0 + tr_z_zzzz_zz[i] * pa_z[i];
    }

}

} // diprec namespace

