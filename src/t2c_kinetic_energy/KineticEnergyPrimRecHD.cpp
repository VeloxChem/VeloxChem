#include "KineticEnergyPrimRecHD.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_hd(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_hd,
                            const size_t idx_ovl_fd,
                            const size_t idx_kin_fd,
                            const size_t idx_kin_gp,
                            const size_t idx_kin_gd,
                            const size_t idx_ovl_hd,
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

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_ovl_fd + 6);

    auto ts_xxy_xz = pbuffer.data(idx_ovl_fd + 8);

    auto ts_xxz_xx = pbuffer.data(idx_ovl_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_ovl_fd + 13);

    auto ts_xyy_xy = pbuffer.data(idx_ovl_fd + 19);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_ovl_fd + 23);

    auto ts_xzz_xz = pbuffer.data(idx_ovl_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_ovl_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_ovl_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_ovl_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_ovl_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_ovl_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_ovl_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_ovl_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_ovl_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_ovl_fd + 41);

    auto ts_yyz_xy = pbuffer.data(idx_ovl_fd + 43);

    auto ts_yyz_yy = pbuffer.data(idx_ovl_fd + 45);

    auto ts_yzz_xx = pbuffer.data(idx_ovl_fd + 48);

    auto ts_yzz_xz = pbuffer.data(idx_ovl_fd + 50);

    auto ts_yzz_yz = pbuffer.data(idx_ovl_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_ovl_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_ovl_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_ovl_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_ovl_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_ovl_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_ovl_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_ovl_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto tk_xxx_xx = pbuffer.data(idx_kin_fd);

    auto tk_xxx_xy = pbuffer.data(idx_kin_fd + 1);

    auto tk_xxx_xz = pbuffer.data(idx_kin_fd + 2);

    auto tk_xxx_yy = pbuffer.data(idx_kin_fd + 3);

    auto tk_xxx_yz = pbuffer.data(idx_kin_fd + 4);

    auto tk_xxx_zz = pbuffer.data(idx_kin_fd + 5);

    auto tk_xxy_xx = pbuffer.data(idx_kin_fd + 6);

    auto tk_xxy_xz = pbuffer.data(idx_kin_fd + 8);

    auto tk_xxz_xx = pbuffer.data(idx_kin_fd + 12);

    auto tk_xxz_xy = pbuffer.data(idx_kin_fd + 13);

    auto tk_xyy_xy = pbuffer.data(idx_kin_fd + 19);

    auto tk_xyy_yy = pbuffer.data(idx_kin_fd + 21);

    auto tk_xyy_yz = pbuffer.data(idx_kin_fd + 22);

    auto tk_xyy_zz = pbuffer.data(idx_kin_fd + 23);

    auto tk_xzz_xz = pbuffer.data(idx_kin_fd + 32);

    auto tk_xzz_yy = pbuffer.data(idx_kin_fd + 33);

    auto tk_xzz_yz = pbuffer.data(idx_kin_fd + 34);

    auto tk_xzz_zz = pbuffer.data(idx_kin_fd + 35);

    auto tk_yyy_xx = pbuffer.data(idx_kin_fd + 36);

    auto tk_yyy_xy = pbuffer.data(idx_kin_fd + 37);

    auto tk_yyy_xz = pbuffer.data(idx_kin_fd + 38);

    auto tk_yyy_yy = pbuffer.data(idx_kin_fd + 39);

    auto tk_yyy_yz = pbuffer.data(idx_kin_fd + 40);

    auto tk_yyy_zz = pbuffer.data(idx_kin_fd + 41);

    auto tk_yyz_xy = pbuffer.data(idx_kin_fd + 43);

    auto tk_yyz_yy = pbuffer.data(idx_kin_fd + 45);

    auto tk_yzz_xx = pbuffer.data(idx_kin_fd + 48);

    auto tk_yzz_xz = pbuffer.data(idx_kin_fd + 50);

    auto tk_yzz_yz = pbuffer.data(idx_kin_fd + 52);

    auto tk_yzz_zz = pbuffer.data(idx_kin_fd + 53);

    auto tk_zzz_xx = pbuffer.data(idx_kin_fd + 54);

    auto tk_zzz_xy = pbuffer.data(idx_kin_fd + 55);

    auto tk_zzz_xz = pbuffer.data(idx_kin_fd + 56);

    auto tk_zzz_yy = pbuffer.data(idx_kin_fd + 57);

    auto tk_zzz_yz = pbuffer.data(idx_kin_fd + 58);

    auto tk_zzz_zz = pbuffer.data(idx_kin_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto tk_xxxx_x = pbuffer.data(idx_kin_gp);

    auto tk_xxxx_y = pbuffer.data(idx_kin_gp + 1);

    auto tk_xxxx_z = pbuffer.data(idx_kin_gp + 2);

    auto tk_xxxz_z = pbuffer.data(idx_kin_gp + 8);

    auto tk_xxyy_x = pbuffer.data(idx_kin_gp + 9);

    auto tk_xxyy_y = pbuffer.data(idx_kin_gp + 10);

    auto tk_xxyy_z = pbuffer.data(idx_kin_gp + 11);

    auto tk_xxzz_x = pbuffer.data(idx_kin_gp + 15);

    auto tk_xxzz_y = pbuffer.data(idx_kin_gp + 16);

    auto tk_xxzz_z = pbuffer.data(idx_kin_gp + 17);

    auto tk_xyyy_y = pbuffer.data(idx_kin_gp + 19);

    auto tk_xzzz_z = pbuffer.data(idx_kin_gp + 29);

    auto tk_yyyy_x = pbuffer.data(idx_kin_gp + 30);

    auto tk_yyyy_y = pbuffer.data(idx_kin_gp + 31);

    auto tk_yyyy_z = pbuffer.data(idx_kin_gp + 32);

    auto tk_yyyz_z = pbuffer.data(idx_kin_gp + 35);

    auto tk_yyzz_x = pbuffer.data(idx_kin_gp + 36);

    auto tk_yyzz_y = pbuffer.data(idx_kin_gp + 37);

    auto tk_yyzz_z = pbuffer.data(idx_kin_gp + 38);

    auto tk_yzzz_y = pbuffer.data(idx_kin_gp + 40);

    auto tk_yzzz_z = pbuffer.data(idx_kin_gp + 41);

    auto tk_zzzz_x = pbuffer.data(idx_kin_gp + 42);

    auto tk_zzzz_y = pbuffer.data(idx_kin_gp + 43);

    auto tk_zzzz_z = pbuffer.data(idx_kin_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto tk_xxxx_xx = pbuffer.data(idx_kin_gd);

    auto tk_xxxx_xy = pbuffer.data(idx_kin_gd + 1);

    auto tk_xxxx_xz = pbuffer.data(idx_kin_gd + 2);

    auto tk_xxxx_yy = pbuffer.data(idx_kin_gd + 3);

    auto tk_xxxx_yz = pbuffer.data(idx_kin_gd + 4);

    auto tk_xxxx_zz = pbuffer.data(idx_kin_gd + 5);

    auto tk_xxxy_xx = pbuffer.data(idx_kin_gd + 6);

    auto tk_xxxy_xy = pbuffer.data(idx_kin_gd + 7);

    auto tk_xxxy_xz = pbuffer.data(idx_kin_gd + 8);

    auto tk_xxxy_yy = pbuffer.data(idx_kin_gd + 9);

    auto tk_xxxz_xx = pbuffer.data(idx_kin_gd + 12);

    auto tk_xxxz_xy = pbuffer.data(idx_kin_gd + 13);

    auto tk_xxxz_xz = pbuffer.data(idx_kin_gd + 14);

    auto tk_xxxz_yz = pbuffer.data(idx_kin_gd + 16);

    auto tk_xxxz_zz = pbuffer.data(idx_kin_gd + 17);

    auto tk_xxyy_xx = pbuffer.data(idx_kin_gd + 18);

    auto tk_xxyy_xy = pbuffer.data(idx_kin_gd + 19);

    auto tk_xxyy_xz = pbuffer.data(idx_kin_gd + 20);

    auto tk_xxyy_yy = pbuffer.data(idx_kin_gd + 21);

    auto tk_xxyy_yz = pbuffer.data(idx_kin_gd + 22);

    auto tk_xxyy_zz = pbuffer.data(idx_kin_gd + 23);

    auto tk_xxzz_xx = pbuffer.data(idx_kin_gd + 30);

    auto tk_xxzz_xy = pbuffer.data(idx_kin_gd + 31);

    auto tk_xxzz_xz = pbuffer.data(idx_kin_gd + 32);

    auto tk_xxzz_yy = pbuffer.data(idx_kin_gd + 33);

    auto tk_xxzz_yz = pbuffer.data(idx_kin_gd + 34);

    auto tk_xxzz_zz = pbuffer.data(idx_kin_gd + 35);

    auto tk_xyyy_xx = pbuffer.data(idx_kin_gd + 36);

    auto tk_xyyy_xy = pbuffer.data(idx_kin_gd + 37);

    auto tk_xyyy_yy = pbuffer.data(idx_kin_gd + 39);

    auto tk_xyyy_yz = pbuffer.data(idx_kin_gd + 40);

    auto tk_xyyy_zz = pbuffer.data(idx_kin_gd + 41);

    auto tk_xzzz_xx = pbuffer.data(idx_kin_gd + 54);

    auto tk_xzzz_xz = pbuffer.data(idx_kin_gd + 56);

    auto tk_xzzz_yy = pbuffer.data(idx_kin_gd + 57);

    auto tk_xzzz_yz = pbuffer.data(idx_kin_gd + 58);

    auto tk_xzzz_zz = pbuffer.data(idx_kin_gd + 59);

    auto tk_yyyy_xx = pbuffer.data(idx_kin_gd + 60);

    auto tk_yyyy_xy = pbuffer.data(idx_kin_gd + 61);

    auto tk_yyyy_xz = pbuffer.data(idx_kin_gd + 62);

    auto tk_yyyy_yy = pbuffer.data(idx_kin_gd + 63);

    auto tk_yyyy_yz = pbuffer.data(idx_kin_gd + 64);

    auto tk_yyyy_zz = pbuffer.data(idx_kin_gd + 65);

    auto tk_yyyz_xy = pbuffer.data(idx_kin_gd + 67);

    auto tk_yyyz_xz = pbuffer.data(idx_kin_gd + 68);

    auto tk_yyyz_yy = pbuffer.data(idx_kin_gd + 69);

    auto tk_yyyz_yz = pbuffer.data(idx_kin_gd + 70);

    auto tk_yyyz_zz = pbuffer.data(idx_kin_gd + 71);

    auto tk_yyzz_xx = pbuffer.data(idx_kin_gd + 72);

    auto tk_yyzz_xy = pbuffer.data(idx_kin_gd + 73);

    auto tk_yyzz_xz = pbuffer.data(idx_kin_gd + 74);

    auto tk_yyzz_yy = pbuffer.data(idx_kin_gd + 75);

    auto tk_yyzz_yz = pbuffer.data(idx_kin_gd + 76);

    auto tk_yyzz_zz = pbuffer.data(idx_kin_gd + 77);

    auto tk_yzzz_xx = pbuffer.data(idx_kin_gd + 78);

    auto tk_yzzz_xy = pbuffer.data(idx_kin_gd + 79);

    auto tk_yzzz_xz = pbuffer.data(idx_kin_gd + 80);

    auto tk_yzzz_yy = pbuffer.data(idx_kin_gd + 81);

    auto tk_yzzz_yz = pbuffer.data(idx_kin_gd + 82);

    auto tk_yzzz_zz = pbuffer.data(idx_kin_gd + 83);

    auto tk_zzzz_xx = pbuffer.data(idx_kin_gd + 84);

    auto tk_zzzz_xy = pbuffer.data(idx_kin_gd + 85);

    auto tk_zzzz_xz = pbuffer.data(idx_kin_gd + 86);

    auto tk_zzzz_yy = pbuffer.data(idx_kin_gd + 87);

    auto tk_zzzz_yz = pbuffer.data(idx_kin_gd + 88);

    auto tk_zzzz_zz = pbuffer.data(idx_kin_gd + 89);

    // Set up components of auxiliary buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_ovl_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_ovl_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_ovl_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_ovl_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_ovl_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_ovl_hd + 5);

    auto ts_xxxxy_xx = pbuffer.data(idx_ovl_hd + 6);

    auto ts_xxxxy_xy = pbuffer.data(idx_ovl_hd + 7);

    auto ts_xxxxy_xz = pbuffer.data(idx_ovl_hd + 8);

    auto ts_xxxxy_yy = pbuffer.data(idx_ovl_hd + 9);

    auto ts_xxxxy_yz = pbuffer.data(idx_ovl_hd + 10);

    auto ts_xxxxy_zz = pbuffer.data(idx_ovl_hd + 11);

    auto ts_xxxxz_xx = pbuffer.data(idx_ovl_hd + 12);

    auto ts_xxxxz_xy = pbuffer.data(idx_ovl_hd + 13);

    auto ts_xxxxz_xz = pbuffer.data(idx_ovl_hd + 14);

    auto ts_xxxxz_yy = pbuffer.data(idx_ovl_hd + 15);

    auto ts_xxxxz_yz = pbuffer.data(idx_ovl_hd + 16);

    auto ts_xxxxz_zz = pbuffer.data(idx_ovl_hd + 17);

    auto ts_xxxyy_xx = pbuffer.data(idx_ovl_hd + 18);

    auto ts_xxxyy_xy = pbuffer.data(idx_ovl_hd + 19);

    auto ts_xxxyy_xz = pbuffer.data(idx_ovl_hd + 20);

    auto ts_xxxyy_yy = pbuffer.data(idx_ovl_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_ovl_hd + 22);

    auto ts_xxxyy_zz = pbuffer.data(idx_ovl_hd + 23);

    auto ts_xxxyz_xx = pbuffer.data(idx_ovl_hd + 24);

    auto ts_xxxyz_xy = pbuffer.data(idx_ovl_hd + 25);

    auto ts_xxxyz_xz = pbuffer.data(idx_ovl_hd + 26);

    auto ts_xxxyz_yy = pbuffer.data(idx_ovl_hd + 27);

    auto ts_xxxyz_yz = pbuffer.data(idx_ovl_hd + 28);

    auto ts_xxxyz_zz = pbuffer.data(idx_ovl_hd + 29);

    auto ts_xxxzz_xx = pbuffer.data(idx_ovl_hd + 30);

    auto ts_xxxzz_xy = pbuffer.data(idx_ovl_hd + 31);

    auto ts_xxxzz_xz = pbuffer.data(idx_ovl_hd + 32);

    auto ts_xxxzz_yy = pbuffer.data(idx_ovl_hd + 33);

    auto ts_xxxzz_yz = pbuffer.data(idx_ovl_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_ovl_hd + 35);

    auto ts_xxyyy_xx = pbuffer.data(idx_ovl_hd + 36);

    auto ts_xxyyy_xy = pbuffer.data(idx_ovl_hd + 37);

    auto ts_xxyyy_xz = pbuffer.data(idx_ovl_hd + 38);

    auto ts_xxyyy_yy = pbuffer.data(idx_ovl_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_ovl_hd + 40);

    auto ts_xxyyy_zz = pbuffer.data(idx_ovl_hd + 41);

    auto ts_xxyyz_xx = pbuffer.data(idx_ovl_hd + 42);

    auto ts_xxyyz_xy = pbuffer.data(idx_ovl_hd + 43);

    auto ts_xxyyz_xz = pbuffer.data(idx_ovl_hd + 44);

    auto ts_xxyyz_yy = pbuffer.data(idx_ovl_hd + 45);

    auto ts_xxyyz_yz = pbuffer.data(idx_ovl_hd + 46);

    auto ts_xxyyz_zz = pbuffer.data(idx_ovl_hd + 47);

    auto ts_xxyzz_xx = pbuffer.data(idx_ovl_hd + 48);

    auto ts_xxyzz_xy = pbuffer.data(idx_ovl_hd + 49);

    auto ts_xxyzz_xz = pbuffer.data(idx_ovl_hd + 50);

    auto ts_xxyzz_yy = pbuffer.data(idx_ovl_hd + 51);

    auto ts_xxyzz_yz = pbuffer.data(idx_ovl_hd + 52);

    auto ts_xxyzz_zz = pbuffer.data(idx_ovl_hd + 53);

    auto ts_xxzzz_xx = pbuffer.data(idx_ovl_hd + 54);

    auto ts_xxzzz_xy = pbuffer.data(idx_ovl_hd + 55);

    auto ts_xxzzz_xz = pbuffer.data(idx_ovl_hd + 56);

    auto ts_xxzzz_yy = pbuffer.data(idx_ovl_hd + 57);

    auto ts_xxzzz_yz = pbuffer.data(idx_ovl_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_ovl_hd + 59);

    auto ts_xyyyy_xx = pbuffer.data(idx_ovl_hd + 60);

    auto ts_xyyyy_xy = pbuffer.data(idx_ovl_hd + 61);

    auto ts_xyyyy_xz = pbuffer.data(idx_ovl_hd + 62);

    auto ts_xyyyy_yy = pbuffer.data(idx_ovl_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_ovl_hd + 64);

    auto ts_xyyyy_zz = pbuffer.data(idx_ovl_hd + 65);

    auto ts_xyyyz_xx = pbuffer.data(idx_ovl_hd + 66);

    auto ts_xyyyz_xy = pbuffer.data(idx_ovl_hd + 67);

    auto ts_xyyyz_xz = pbuffer.data(idx_ovl_hd + 68);

    auto ts_xyyyz_yy = pbuffer.data(idx_ovl_hd + 69);

    auto ts_xyyyz_yz = pbuffer.data(idx_ovl_hd + 70);

    auto ts_xyyyz_zz = pbuffer.data(idx_ovl_hd + 71);

    auto ts_xyyzz_xx = pbuffer.data(idx_ovl_hd + 72);

    auto ts_xyyzz_xy = pbuffer.data(idx_ovl_hd + 73);

    auto ts_xyyzz_xz = pbuffer.data(idx_ovl_hd + 74);

    auto ts_xyyzz_yy = pbuffer.data(idx_ovl_hd + 75);

    auto ts_xyyzz_yz = pbuffer.data(idx_ovl_hd + 76);

    auto ts_xyyzz_zz = pbuffer.data(idx_ovl_hd + 77);

    auto ts_xyzzz_xx = pbuffer.data(idx_ovl_hd + 78);

    auto ts_xyzzz_xy = pbuffer.data(idx_ovl_hd + 79);

    auto ts_xyzzz_xz = pbuffer.data(idx_ovl_hd + 80);

    auto ts_xyzzz_yy = pbuffer.data(idx_ovl_hd + 81);

    auto ts_xyzzz_yz = pbuffer.data(idx_ovl_hd + 82);

    auto ts_xyzzz_zz = pbuffer.data(idx_ovl_hd + 83);

    auto ts_xzzzz_xx = pbuffer.data(idx_ovl_hd + 84);

    auto ts_xzzzz_xy = pbuffer.data(idx_ovl_hd + 85);

    auto ts_xzzzz_xz = pbuffer.data(idx_ovl_hd + 86);

    auto ts_xzzzz_yy = pbuffer.data(idx_ovl_hd + 87);

    auto ts_xzzzz_yz = pbuffer.data(idx_ovl_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_ovl_hd + 89);

    auto ts_yyyyy_xx = pbuffer.data(idx_ovl_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_ovl_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_ovl_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_ovl_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_ovl_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_ovl_hd + 95);

    auto ts_yyyyz_xx = pbuffer.data(idx_ovl_hd + 96);

    auto ts_yyyyz_xy = pbuffer.data(idx_ovl_hd + 97);

    auto ts_yyyyz_xz = pbuffer.data(idx_ovl_hd + 98);

    auto ts_yyyyz_yy = pbuffer.data(idx_ovl_hd + 99);

    auto ts_yyyyz_yz = pbuffer.data(idx_ovl_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_ovl_hd + 101);

    auto ts_yyyzz_xx = pbuffer.data(idx_ovl_hd + 102);

    auto ts_yyyzz_xy = pbuffer.data(idx_ovl_hd + 103);

    auto ts_yyyzz_xz = pbuffer.data(idx_ovl_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_ovl_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_ovl_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_ovl_hd + 107);

    auto ts_yyzzz_xx = pbuffer.data(idx_ovl_hd + 108);

    auto ts_yyzzz_xy = pbuffer.data(idx_ovl_hd + 109);

    auto ts_yyzzz_xz = pbuffer.data(idx_ovl_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_ovl_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_ovl_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_ovl_hd + 113);

    auto ts_yzzzz_xx = pbuffer.data(idx_ovl_hd + 114);

    auto ts_yzzzz_xy = pbuffer.data(idx_ovl_hd + 115);

    auto ts_yzzzz_xz = pbuffer.data(idx_ovl_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_ovl_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_ovl_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_ovl_hd + 119);

    auto ts_zzzzz_xx = pbuffer.data(idx_ovl_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_ovl_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_ovl_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_ovl_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_ovl_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_ovl_hd + 125);

    // Set up 0-6 components of targeted buffer : HD

    auto tk_xxxxx_xx = pbuffer.data(idx_kin_hd);

    auto tk_xxxxx_xy = pbuffer.data(idx_kin_hd + 1);

    auto tk_xxxxx_xz = pbuffer.data(idx_kin_hd + 2);

    auto tk_xxxxx_yy = pbuffer.data(idx_kin_hd + 3);

    auto tk_xxxxx_yz = pbuffer.data(idx_kin_hd + 4);

    auto tk_xxxxx_zz = pbuffer.data(idx_kin_hd + 5);

    #pragma omp simd aligned(pa_x, tk_xxx_xx, tk_xxx_xy, tk_xxx_xz, tk_xxx_yy, tk_xxx_yz, tk_xxx_zz, tk_xxxx_x, tk_xxxx_xx, tk_xxxx_xy, tk_xxxx_xz, tk_xxxx_y, tk_xxxx_yy, tk_xxxx_yz, tk_xxxx_z, tk_xxxx_zz, tk_xxxxx_xx, tk_xxxxx_xy, tk_xxxxx_xz, tk_xxxxx_yy, tk_xxxxx_yz, tk_xxxxx_zz, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_yy, ts_xxxxx_yz, ts_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxx_xx[i] = -8.0 * ts_xxx_xx[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xx[i] * fe_0 + 2.0 * tk_xxxx_x[i] * fe_0 + tk_xxxx_xx[i] * pa_x[i] + 2.0 * ts_xxxxx_xx[i] * fz_0;

        tk_xxxxx_xy[i] = -8.0 * ts_xxx_xy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xy[i] * fe_0 + tk_xxxx_y[i] * fe_0 + tk_xxxx_xy[i] * pa_x[i] + 2.0 * ts_xxxxx_xy[i] * fz_0;

        tk_xxxxx_xz[i] = -8.0 * ts_xxx_xz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_xz[i] * fe_0 + tk_xxxx_z[i] * fe_0 + tk_xxxx_xz[i] * pa_x[i] + 2.0 * ts_xxxxx_xz[i] * fz_0;

        tk_xxxxx_yy[i] = -8.0 * ts_xxx_yy[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yy[i] * fe_0 + tk_xxxx_yy[i] * pa_x[i] + 2.0 * ts_xxxxx_yy[i] * fz_0;

        tk_xxxxx_yz[i] = -8.0 * ts_xxx_yz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_yz[i] * fe_0 + tk_xxxx_yz[i] * pa_x[i] + 2.0 * ts_xxxxx_yz[i] * fz_0;

        tk_xxxxx_zz[i] = -8.0 * ts_xxx_zz[i] * fbe_0 * fz_0 + 4.0 * tk_xxx_zz[i] * fe_0 + tk_xxxx_zz[i] * pa_x[i] + 2.0 * ts_xxxxx_zz[i] * fz_0;
    }

    // Set up 6-12 components of targeted buffer : HD

    auto tk_xxxxy_xx = pbuffer.data(idx_kin_hd + 6);

    auto tk_xxxxy_xy = pbuffer.data(idx_kin_hd + 7);

    auto tk_xxxxy_xz = pbuffer.data(idx_kin_hd + 8);

    auto tk_xxxxy_yy = pbuffer.data(idx_kin_hd + 9);

    auto tk_xxxxy_yz = pbuffer.data(idx_kin_hd + 10);

    auto tk_xxxxy_zz = pbuffer.data(idx_kin_hd + 11);

    #pragma omp simd aligned(pa_y, tk_xxxx_x, tk_xxxx_xx, tk_xxxx_xy, tk_xxxx_xz, tk_xxxx_y, tk_xxxx_yy, tk_xxxx_yz, tk_xxxx_z, tk_xxxx_zz, tk_xxxxy_xx, tk_xxxxy_xy, tk_xxxxy_xz, tk_xxxxy_yy, tk_xxxxy_yz, tk_xxxxy_zz, ts_xxxxy_xx, ts_xxxxy_xy, ts_xxxxy_xz, ts_xxxxy_yy, ts_xxxxy_yz, ts_xxxxy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxy_xx[i] = tk_xxxx_xx[i] * pa_y[i] + 2.0 * ts_xxxxy_xx[i] * fz_0;

        tk_xxxxy_xy[i] = tk_xxxx_x[i] * fe_0 + tk_xxxx_xy[i] * pa_y[i] + 2.0 * ts_xxxxy_xy[i] * fz_0;

        tk_xxxxy_xz[i] = tk_xxxx_xz[i] * pa_y[i] + 2.0 * ts_xxxxy_xz[i] * fz_0;

        tk_xxxxy_yy[i] = 2.0 * tk_xxxx_y[i] * fe_0 + tk_xxxx_yy[i] * pa_y[i] + 2.0 * ts_xxxxy_yy[i] * fz_0;

        tk_xxxxy_yz[i] = tk_xxxx_z[i] * fe_0 + tk_xxxx_yz[i] * pa_y[i] + 2.0 * ts_xxxxy_yz[i] * fz_0;

        tk_xxxxy_zz[i] = tk_xxxx_zz[i] * pa_y[i] + 2.0 * ts_xxxxy_zz[i] * fz_0;
    }

    // Set up 12-18 components of targeted buffer : HD

    auto tk_xxxxz_xx = pbuffer.data(idx_kin_hd + 12);

    auto tk_xxxxz_xy = pbuffer.data(idx_kin_hd + 13);

    auto tk_xxxxz_xz = pbuffer.data(idx_kin_hd + 14);

    auto tk_xxxxz_yy = pbuffer.data(idx_kin_hd + 15);

    auto tk_xxxxz_yz = pbuffer.data(idx_kin_hd + 16);

    auto tk_xxxxz_zz = pbuffer.data(idx_kin_hd + 17);

    #pragma omp simd aligned(pa_z, tk_xxxx_x, tk_xxxx_xx, tk_xxxx_xy, tk_xxxx_xz, tk_xxxx_y, tk_xxxx_yy, tk_xxxx_yz, tk_xxxx_z, tk_xxxx_zz, tk_xxxxz_xx, tk_xxxxz_xy, tk_xxxxz_xz, tk_xxxxz_yy, tk_xxxxz_yz, tk_xxxxz_zz, ts_xxxxz_xx, ts_xxxxz_xy, ts_xxxxz_xz, ts_xxxxz_yy, ts_xxxxz_yz, ts_xxxxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxz_xx[i] = tk_xxxx_xx[i] * pa_z[i] + 2.0 * ts_xxxxz_xx[i] * fz_0;

        tk_xxxxz_xy[i] = tk_xxxx_xy[i] * pa_z[i] + 2.0 * ts_xxxxz_xy[i] * fz_0;

        tk_xxxxz_xz[i] = tk_xxxx_x[i] * fe_0 + tk_xxxx_xz[i] * pa_z[i] + 2.0 * ts_xxxxz_xz[i] * fz_0;

        tk_xxxxz_yy[i] = tk_xxxx_yy[i] * pa_z[i] + 2.0 * ts_xxxxz_yy[i] * fz_0;

        tk_xxxxz_yz[i] = tk_xxxx_y[i] * fe_0 + tk_xxxx_yz[i] * pa_z[i] + 2.0 * ts_xxxxz_yz[i] * fz_0;

        tk_xxxxz_zz[i] = 2.0 * tk_xxxx_z[i] * fe_0 + tk_xxxx_zz[i] * pa_z[i] + 2.0 * ts_xxxxz_zz[i] * fz_0;
    }

    // Set up 18-24 components of targeted buffer : HD

    auto tk_xxxyy_xx = pbuffer.data(idx_kin_hd + 18);

    auto tk_xxxyy_xy = pbuffer.data(idx_kin_hd + 19);

    auto tk_xxxyy_xz = pbuffer.data(idx_kin_hd + 20);

    auto tk_xxxyy_yy = pbuffer.data(idx_kin_hd + 21);

    auto tk_xxxyy_yz = pbuffer.data(idx_kin_hd + 22);

    auto tk_xxxyy_zz = pbuffer.data(idx_kin_hd + 23);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxx_xx, tk_xxx_xz, tk_xxxy_xx, tk_xxxy_xz, tk_xxxyy_xx, tk_xxxyy_xy, tk_xxxyy_xz, tk_xxxyy_yy, tk_xxxyy_yz, tk_xxxyy_zz, tk_xxyy_xy, tk_xxyy_y, tk_xxyy_yy, tk_xxyy_yz, tk_xxyy_zz, tk_xyy_xy, tk_xyy_yy, tk_xyy_yz, tk_xyy_zz, ts_xxx_xx, ts_xxx_xz, ts_xxxyy_xx, ts_xxxyy_xy, ts_xxxyy_xz, ts_xxxyy_yy, ts_xxxyy_yz, ts_xxxyy_zz, ts_xyy_xy, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyy_xx[i] = -2.0 * ts_xxx_xx[i] * fbe_0 * fz_0 + tk_xxx_xx[i] * fe_0 + tk_xxxy_xx[i] * pa_y[i] + 2.0 * ts_xxxyy_xx[i] * fz_0;

        tk_xxxyy_xy[i] = -4.0 * ts_xyy_xy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_xy[i] * fe_0 + tk_xxyy_y[i] * fe_0 + tk_xxyy_xy[i] * pa_x[i] + 2.0 * ts_xxxyy_xy[i] * fz_0;

        tk_xxxyy_xz[i] = -2.0 * ts_xxx_xz[i] * fbe_0 * fz_0 + tk_xxx_xz[i] * fe_0 + tk_xxxy_xz[i] * pa_y[i] + 2.0 * ts_xxxyy_xz[i] * fz_0;

        tk_xxxyy_yy[i] = -4.0 * ts_xyy_yy[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yy[i] * fe_0 + tk_xxyy_yy[i] * pa_x[i] + 2.0 * ts_xxxyy_yy[i] * fz_0;

        tk_xxxyy_yz[i] = -4.0 * ts_xyy_yz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_yz[i] * fe_0 + tk_xxyy_yz[i] * pa_x[i] + 2.0 * ts_xxxyy_yz[i] * fz_0;

        tk_xxxyy_zz[i] = -4.0 * ts_xyy_zz[i] * fbe_0 * fz_0 + 2.0 * tk_xyy_zz[i] * fe_0 + tk_xxyy_zz[i] * pa_x[i] + 2.0 * ts_xxxyy_zz[i] * fz_0;
    }

    // Set up 24-30 components of targeted buffer : HD

    auto tk_xxxyz_xx = pbuffer.data(idx_kin_hd + 24);

    auto tk_xxxyz_xy = pbuffer.data(idx_kin_hd + 25);

    auto tk_xxxyz_xz = pbuffer.data(idx_kin_hd + 26);

    auto tk_xxxyz_yy = pbuffer.data(idx_kin_hd + 27);

    auto tk_xxxyz_yz = pbuffer.data(idx_kin_hd + 28);

    auto tk_xxxyz_zz = pbuffer.data(idx_kin_hd + 29);

    #pragma omp simd aligned(pa_y, pa_z, tk_xxxy_xy, tk_xxxy_yy, tk_xxxyz_xx, tk_xxxyz_xy, tk_xxxyz_xz, tk_xxxyz_yy, tk_xxxyz_yz, tk_xxxyz_zz, tk_xxxz_xx, tk_xxxz_xz, tk_xxxz_yz, tk_xxxz_z, tk_xxxz_zz, ts_xxxyz_xx, ts_xxxyz_xy, ts_xxxyz_xz, ts_xxxyz_yy, ts_xxxyz_yz, ts_xxxyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyz_xx[i] = tk_xxxz_xx[i] * pa_y[i] + 2.0 * ts_xxxyz_xx[i] * fz_0;

        tk_xxxyz_xy[i] = tk_xxxy_xy[i] * pa_z[i] + 2.0 * ts_xxxyz_xy[i] * fz_0;

        tk_xxxyz_xz[i] = tk_xxxz_xz[i] * pa_y[i] + 2.0 * ts_xxxyz_xz[i] * fz_0;

        tk_xxxyz_yy[i] = tk_xxxy_yy[i] * pa_z[i] + 2.0 * ts_xxxyz_yy[i] * fz_0;

        tk_xxxyz_yz[i] = tk_xxxz_z[i] * fe_0 + tk_xxxz_yz[i] * pa_y[i] + 2.0 * ts_xxxyz_yz[i] * fz_0;

        tk_xxxyz_zz[i] = tk_xxxz_zz[i] * pa_y[i] + 2.0 * ts_xxxyz_zz[i] * fz_0;
    }

    // Set up 30-36 components of targeted buffer : HD

    auto tk_xxxzz_xx = pbuffer.data(idx_kin_hd + 30);

    auto tk_xxxzz_xy = pbuffer.data(idx_kin_hd + 31);

    auto tk_xxxzz_xz = pbuffer.data(idx_kin_hd + 32);

    auto tk_xxxzz_yy = pbuffer.data(idx_kin_hd + 33);

    auto tk_xxxzz_yz = pbuffer.data(idx_kin_hd + 34);

    auto tk_xxxzz_zz = pbuffer.data(idx_kin_hd + 35);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxx_xx, tk_xxx_xy, tk_xxxz_xx, tk_xxxz_xy, tk_xxxzz_xx, tk_xxxzz_xy, tk_xxxzz_xz, tk_xxxzz_yy, tk_xxxzz_yz, tk_xxxzz_zz, tk_xxzz_xz, tk_xxzz_yy, tk_xxzz_yz, tk_xxzz_z, tk_xxzz_zz, tk_xzz_xz, tk_xzz_yy, tk_xzz_yz, tk_xzz_zz, ts_xxx_xx, ts_xxx_xy, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzz_xz, ts_xxxzz_yy, ts_xxxzz_yz, ts_xxxzz_zz, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzz_xx[i] = -2.0 * ts_xxx_xx[i] * fbe_0 * fz_0 + tk_xxx_xx[i] * fe_0 + tk_xxxz_xx[i] * pa_z[i] + 2.0 * ts_xxxzz_xx[i] * fz_0;

        tk_xxxzz_xy[i] = -2.0 * ts_xxx_xy[i] * fbe_0 * fz_0 + tk_xxx_xy[i] * fe_0 + tk_xxxz_xy[i] * pa_z[i] + 2.0 * ts_xxxzz_xy[i] * fz_0;

        tk_xxxzz_xz[i] = -4.0 * ts_xzz_xz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_xz[i] * fe_0 + tk_xxzz_z[i] * fe_0 + tk_xxzz_xz[i] * pa_x[i] + 2.0 * ts_xxxzz_xz[i] * fz_0;

        tk_xxxzz_yy[i] = -4.0 * ts_xzz_yy[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yy[i] * fe_0 + tk_xxzz_yy[i] * pa_x[i] + 2.0 * ts_xxxzz_yy[i] * fz_0;

        tk_xxxzz_yz[i] = -4.0 * ts_xzz_yz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_yz[i] * fe_0 + tk_xxzz_yz[i] * pa_x[i] + 2.0 * ts_xxxzz_yz[i] * fz_0;

        tk_xxxzz_zz[i] = -4.0 * ts_xzz_zz[i] * fbe_0 * fz_0 + 2.0 * tk_xzz_zz[i] * fe_0 + tk_xxzz_zz[i] * pa_x[i] + 2.0 * ts_xxxzz_zz[i] * fz_0;
    }

    // Set up 36-42 components of targeted buffer : HD

    auto tk_xxyyy_xx = pbuffer.data(idx_kin_hd + 36);

    auto tk_xxyyy_xy = pbuffer.data(idx_kin_hd + 37);

    auto tk_xxyyy_xz = pbuffer.data(idx_kin_hd + 38);

    auto tk_xxyyy_yy = pbuffer.data(idx_kin_hd + 39);

    auto tk_xxyyy_yz = pbuffer.data(idx_kin_hd + 40);

    auto tk_xxyyy_zz = pbuffer.data(idx_kin_hd + 41);

    #pragma omp simd aligned(pa_x, pa_y, tk_xxy_xx, tk_xxy_xz, tk_xxyy_xx, tk_xxyy_xz, tk_xxyyy_xx, tk_xxyyy_xy, tk_xxyyy_xz, tk_xxyyy_yy, tk_xxyyy_yz, tk_xxyyy_zz, tk_xyyy_xy, tk_xyyy_y, tk_xyyy_yy, tk_xyyy_yz, tk_xyyy_zz, tk_yyy_xy, tk_yyy_yy, tk_yyy_yz, tk_yyy_zz, ts_xxy_xx, ts_xxy_xz, ts_xxyyy_xx, ts_xxyyy_xy, ts_xxyyy_xz, ts_xxyyy_yy, ts_xxyyy_yz, ts_xxyyy_zz, ts_yyy_xy, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyy_xx[i] = -4.0 * ts_xxy_xx[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xx[i] * fe_0 + tk_xxyy_xx[i] * pa_y[i] + 2.0 * ts_xxyyy_xx[i] * fz_0;

        tk_xxyyy_xy[i] = -2.0 * ts_yyy_xy[i] * fbe_0 * fz_0 + tk_yyy_xy[i] * fe_0 + tk_xyyy_y[i] * fe_0 + tk_xyyy_xy[i] * pa_x[i] + 2.0 * ts_xxyyy_xy[i] * fz_0;

        tk_xxyyy_xz[i] = -4.0 * ts_xxy_xz[i] * fbe_0 * fz_0 + 2.0 * tk_xxy_xz[i] * fe_0 + tk_xxyy_xz[i] * pa_y[i] + 2.0 * ts_xxyyy_xz[i] * fz_0;

        tk_xxyyy_yy[i] = -2.0 * ts_yyy_yy[i] * fbe_0 * fz_0 + tk_yyy_yy[i] * fe_0 + tk_xyyy_yy[i] * pa_x[i] + 2.0 * ts_xxyyy_yy[i] * fz_0;

        tk_xxyyy_yz[i] = -2.0 * ts_yyy_yz[i] * fbe_0 * fz_0 + tk_yyy_yz[i] * fe_0 + tk_xyyy_yz[i] * pa_x[i] + 2.0 * ts_xxyyy_yz[i] * fz_0;

        tk_xxyyy_zz[i] = -2.0 * ts_yyy_zz[i] * fbe_0 * fz_0 + tk_yyy_zz[i] * fe_0 + tk_xyyy_zz[i] * pa_x[i] + 2.0 * ts_xxyyy_zz[i] * fz_0;
    }

    // Set up 42-48 components of targeted buffer : HD

    auto tk_xxyyz_xx = pbuffer.data(idx_kin_hd + 42);

    auto tk_xxyyz_xy = pbuffer.data(idx_kin_hd + 43);

    auto tk_xxyyz_xz = pbuffer.data(idx_kin_hd + 44);

    auto tk_xxyyz_yy = pbuffer.data(idx_kin_hd + 45);

    auto tk_xxyyz_yz = pbuffer.data(idx_kin_hd + 46);

    auto tk_xxyyz_zz = pbuffer.data(idx_kin_hd + 47);

    #pragma omp simd aligned(pa_z, tk_xxyy_x, tk_xxyy_xx, tk_xxyy_xy, tk_xxyy_xz, tk_xxyy_y, tk_xxyy_yy, tk_xxyy_yz, tk_xxyy_z, tk_xxyy_zz, tk_xxyyz_xx, tk_xxyyz_xy, tk_xxyyz_xz, tk_xxyyz_yy, tk_xxyyz_yz, tk_xxyyz_zz, ts_xxyyz_xx, ts_xxyyz_xy, ts_xxyyz_xz, ts_xxyyz_yy, ts_xxyyz_yz, ts_xxyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyz_xx[i] = tk_xxyy_xx[i] * pa_z[i] + 2.0 * ts_xxyyz_xx[i] * fz_0;

        tk_xxyyz_xy[i] = tk_xxyy_xy[i] * pa_z[i] + 2.0 * ts_xxyyz_xy[i] * fz_0;

        tk_xxyyz_xz[i] = tk_xxyy_x[i] * fe_0 + tk_xxyy_xz[i] * pa_z[i] + 2.0 * ts_xxyyz_xz[i] * fz_0;

        tk_xxyyz_yy[i] = tk_xxyy_yy[i] * pa_z[i] + 2.0 * ts_xxyyz_yy[i] * fz_0;

        tk_xxyyz_yz[i] = tk_xxyy_y[i] * fe_0 + tk_xxyy_yz[i] * pa_z[i] + 2.0 * ts_xxyyz_yz[i] * fz_0;

        tk_xxyyz_zz[i] = 2.0 * tk_xxyy_z[i] * fe_0 + tk_xxyy_zz[i] * pa_z[i] + 2.0 * ts_xxyyz_zz[i] * fz_0;
    }

    // Set up 48-54 components of targeted buffer : HD

    auto tk_xxyzz_xx = pbuffer.data(idx_kin_hd + 48);

    auto tk_xxyzz_xy = pbuffer.data(idx_kin_hd + 49);

    auto tk_xxyzz_xz = pbuffer.data(idx_kin_hd + 50);

    auto tk_xxyzz_yy = pbuffer.data(idx_kin_hd + 51);

    auto tk_xxyzz_yz = pbuffer.data(idx_kin_hd + 52);

    auto tk_xxyzz_zz = pbuffer.data(idx_kin_hd + 53);

    #pragma omp simd aligned(pa_y, tk_xxyzz_xx, tk_xxyzz_xy, tk_xxyzz_xz, tk_xxyzz_yy, tk_xxyzz_yz, tk_xxyzz_zz, tk_xxzz_x, tk_xxzz_xx, tk_xxzz_xy, tk_xxzz_xz, tk_xxzz_y, tk_xxzz_yy, tk_xxzz_yz, tk_xxzz_z, tk_xxzz_zz, ts_xxyzz_xx, ts_xxyzz_xy, ts_xxyzz_xz, ts_xxyzz_yy, ts_xxyzz_yz, ts_xxyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzz_xx[i] = tk_xxzz_xx[i] * pa_y[i] + 2.0 * ts_xxyzz_xx[i] * fz_0;

        tk_xxyzz_xy[i] = tk_xxzz_x[i] * fe_0 + tk_xxzz_xy[i] * pa_y[i] + 2.0 * ts_xxyzz_xy[i] * fz_0;

        tk_xxyzz_xz[i] = tk_xxzz_xz[i] * pa_y[i] + 2.0 * ts_xxyzz_xz[i] * fz_0;

        tk_xxyzz_yy[i] = 2.0 * tk_xxzz_y[i] * fe_0 + tk_xxzz_yy[i] * pa_y[i] + 2.0 * ts_xxyzz_yy[i] * fz_0;

        tk_xxyzz_yz[i] = tk_xxzz_z[i] * fe_0 + tk_xxzz_yz[i] * pa_y[i] + 2.0 * ts_xxyzz_yz[i] * fz_0;

        tk_xxyzz_zz[i] = tk_xxzz_zz[i] * pa_y[i] + 2.0 * ts_xxyzz_zz[i] * fz_0;
    }

    // Set up 54-60 components of targeted buffer : HD

    auto tk_xxzzz_xx = pbuffer.data(idx_kin_hd + 54);

    auto tk_xxzzz_xy = pbuffer.data(idx_kin_hd + 55);

    auto tk_xxzzz_xz = pbuffer.data(idx_kin_hd + 56);

    auto tk_xxzzz_yy = pbuffer.data(idx_kin_hd + 57);

    auto tk_xxzzz_yz = pbuffer.data(idx_kin_hd + 58);

    auto tk_xxzzz_zz = pbuffer.data(idx_kin_hd + 59);

    #pragma omp simd aligned(pa_x, pa_z, tk_xxz_xx, tk_xxz_xy, tk_xxzz_xx, tk_xxzz_xy, tk_xxzzz_xx, tk_xxzzz_xy, tk_xxzzz_xz, tk_xxzzz_yy, tk_xxzzz_yz, tk_xxzzz_zz, tk_xzzz_xz, tk_xzzz_yy, tk_xzzz_yz, tk_xzzz_z, tk_xzzz_zz, tk_zzz_xz, tk_zzz_yy, tk_zzz_yz, tk_zzz_zz, ts_xxz_xx, ts_xxz_xy, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzz_xz, ts_xxzzz_yy, ts_xxzzz_yz, ts_xxzzz_zz, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzz_xx[i] = -4.0 * ts_xxz_xx[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xx[i] * fe_0 + tk_xxzz_xx[i] * pa_z[i] + 2.0 * ts_xxzzz_xx[i] * fz_0;

        tk_xxzzz_xy[i] = -4.0 * ts_xxz_xy[i] * fbe_0 * fz_0 + 2.0 * tk_xxz_xy[i] * fe_0 + tk_xxzz_xy[i] * pa_z[i] + 2.0 * ts_xxzzz_xy[i] * fz_0;

        tk_xxzzz_xz[i] = -2.0 * ts_zzz_xz[i] * fbe_0 * fz_0 + tk_zzz_xz[i] * fe_0 + tk_xzzz_z[i] * fe_0 + tk_xzzz_xz[i] * pa_x[i] + 2.0 * ts_xxzzz_xz[i] * fz_0;

        tk_xxzzz_yy[i] = -2.0 * ts_zzz_yy[i] * fbe_0 * fz_0 + tk_zzz_yy[i] * fe_0 + tk_xzzz_yy[i] * pa_x[i] + 2.0 * ts_xxzzz_yy[i] * fz_0;

        tk_xxzzz_yz[i] = -2.0 * ts_zzz_yz[i] * fbe_0 * fz_0 + tk_zzz_yz[i] * fe_0 + tk_xzzz_yz[i] * pa_x[i] + 2.0 * ts_xxzzz_yz[i] * fz_0;

        tk_xxzzz_zz[i] = -2.0 * ts_zzz_zz[i] * fbe_0 * fz_0 + tk_zzz_zz[i] * fe_0 + tk_xzzz_zz[i] * pa_x[i] + 2.0 * ts_xxzzz_zz[i] * fz_0;
    }

    // Set up 60-66 components of targeted buffer : HD

    auto tk_xyyyy_xx = pbuffer.data(idx_kin_hd + 60);

    auto tk_xyyyy_xy = pbuffer.data(idx_kin_hd + 61);

    auto tk_xyyyy_xz = pbuffer.data(idx_kin_hd + 62);

    auto tk_xyyyy_yy = pbuffer.data(idx_kin_hd + 63);

    auto tk_xyyyy_yz = pbuffer.data(idx_kin_hd + 64);

    auto tk_xyyyy_zz = pbuffer.data(idx_kin_hd + 65);

    #pragma omp simd aligned(pa_x, tk_xyyyy_xx, tk_xyyyy_xy, tk_xyyyy_xz, tk_xyyyy_yy, tk_xyyyy_yz, tk_xyyyy_zz, tk_yyyy_x, tk_yyyy_xx, tk_yyyy_xy, tk_yyyy_xz, tk_yyyy_y, tk_yyyy_yy, tk_yyyy_yz, tk_yyyy_z, tk_yyyy_zz, ts_xyyyy_xx, ts_xyyyy_xy, ts_xyyyy_xz, ts_xyyyy_yy, ts_xyyyy_yz, ts_xyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyy_xx[i] = 2.0 * tk_yyyy_x[i] * fe_0 + tk_yyyy_xx[i] * pa_x[i] + 2.0 * ts_xyyyy_xx[i] * fz_0;

        tk_xyyyy_xy[i] = tk_yyyy_y[i] * fe_0 + tk_yyyy_xy[i] * pa_x[i] + 2.0 * ts_xyyyy_xy[i] * fz_0;

        tk_xyyyy_xz[i] = tk_yyyy_z[i] * fe_0 + tk_yyyy_xz[i] * pa_x[i] + 2.0 * ts_xyyyy_xz[i] * fz_0;

        tk_xyyyy_yy[i] = tk_yyyy_yy[i] * pa_x[i] + 2.0 * ts_xyyyy_yy[i] * fz_0;

        tk_xyyyy_yz[i] = tk_yyyy_yz[i] * pa_x[i] + 2.0 * ts_xyyyy_yz[i] * fz_0;

        tk_xyyyy_zz[i] = tk_yyyy_zz[i] * pa_x[i] + 2.0 * ts_xyyyy_zz[i] * fz_0;
    }

    // Set up 66-72 components of targeted buffer : HD

    auto tk_xyyyz_xx = pbuffer.data(idx_kin_hd + 66);

    auto tk_xyyyz_xy = pbuffer.data(idx_kin_hd + 67);

    auto tk_xyyyz_xz = pbuffer.data(idx_kin_hd + 68);

    auto tk_xyyyz_yy = pbuffer.data(idx_kin_hd + 69);

    auto tk_xyyyz_yz = pbuffer.data(idx_kin_hd + 70);

    auto tk_xyyyz_zz = pbuffer.data(idx_kin_hd + 71);

    #pragma omp simd aligned(pa_x, pa_z, tk_xyyy_xx, tk_xyyy_xy, tk_xyyyz_xx, tk_xyyyz_xy, tk_xyyyz_xz, tk_xyyyz_yy, tk_xyyyz_yz, tk_xyyyz_zz, tk_yyyz_xz, tk_yyyz_yy, tk_yyyz_yz, tk_yyyz_z, tk_yyyz_zz, ts_xyyyz_xx, ts_xyyyz_xy, ts_xyyyz_xz, ts_xyyyz_yy, ts_xyyyz_yz, ts_xyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyz_xx[i] = tk_xyyy_xx[i] * pa_z[i] + 2.0 * ts_xyyyz_xx[i] * fz_0;

        tk_xyyyz_xy[i] = tk_xyyy_xy[i] * pa_z[i] + 2.0 * ts_xyyyz_xy[i] * fz_0;

        tk_xyyyz_xz[i] = tk_yyyz_z[i] * fe_0 + tk_yyyz_xz[i] * pa_x[i] + 2.0 * ts_xyyyz_xz[i] * fz_0;

        tk_xyyyz_yy[i] = tk_yyyz_yy[i] * pa_x[i] + 2.0 * ts_xyyyz_yy[i] * fz_0;

        tk_xyyyz_yz[i] = tk_yyyz_yz[i] * pa_x[i] + 2.0 * ts_xyyyz_yz[i] * fz_0;

        tk_xyyyz_zz[i] = tk_yyyz_zz[i] * pa_x[i] + 2.0 * ts_xyyyz_zz[i] * fz_0;
    }

    // Set up 72-78 components of targeted buffer : HD

    auto tk_xyyzz_xx = pbuffer.data(idx_kin_hd + 72);

    auto tk_xyyzz_xy = pbuffer.data(idx_kin_hd + 73);

    auto tk_xyyzz_xz = pbuffer.data(idx_kin_hd + 74);

    auto tk_xyyzz_yy = pbuffer.data(idx_kin_hd + 75);

    auto tk_xyyzz_yz = pbuffer.data(idx_kin_hd + 76);

    auto tk_xyyzz_zz = pbuffer.data(idx_kin_hd + 77);

    #pragma omp simd aligned(pa_x, tk_xyyzz_xx, tk_xyyzz_xy, tk_xyyzz_xz, tk_xyyzz_yy, tk_xyyzz_yz, tk_xyyzz_zz, tk_yyzz_x, tk_yyzz_xx, tk_yyzz_xy, tk_yyzz_xz, tk_yyzz_y, tk_yyzz_yy, tk_yyzz_yz, tk_yyzz_z, tk_yyzz_zz, ts_xyyzz_xx, ts_xyyzz_xy, ts_xyyzz_xz, ts_xyyzz_yy, ts_xyyzz_yz, ts_xyyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzz_xx[i] = 2.0 * tk_yyzz_x[i] * fe_0 + tk_yyzz_xx[i] * pa_x[i] + 2.0 * ts_xyyzz_xx[i] * fz_0;

        tk_xyyzz_xy[i] = tk_yyzz_y[i] * fe_0 + tk_yyzz_xy[i] * pa_x[i] + 2.0 * ts_xyyzz_xy[i] * fz_0;

        tk_xyyzz_xz[i] = tk_yyzz_z[i] * fe_0 + tk_yyzz_xz[i] * pa_x[i] + 2.0 * ts_xyyzz_xz[i] * fz_0;

        tk_xyyzz_yy[i] = tk_yyzz_yy[i] * pa_x[i] + 2.0 * ts_xyyzz_yy[i] * fz_0;

        tk_xyyzz_yz[i] = tk_yyzz_yz[i] * pa_x[i] + 2.0 * ts_xyyzz_yz[i] * fz_0;

        tk_xyyzz_zz[i] = tk_yyzz_zz[i] * pa_x[i] + 2.0 * ts_xyyzz_zz[i] * fz_0;
    }

    // Set up 78-84 components of targeted buffer : HD

    auto tk_xyzzz_xx = pbuffer.data(idx_kin_hd + 78);

    auto tk_xyzzz_xy = pbuffer.data(idx_kin_hd + 79);

    auto tk_xyzzz_xz = pbuffer.data(idx_kin_hd + 80);

    auto tk_xyzzz_yy = pbuffer.data(idx_kin_hd + 81);

    auto tk_xyzzz_yz = pbuffer.data(idx_kin_hd + 82);

    auto tk_xyzzz_zz = pbuffer.data(idx_kin_hd + 83);

    #pragma omp simd aligned(pa_x, pa_y, tk_xyzzz_xx, tk_xyzzz_xy, tk_xyzzz_xz, tk_xyzzz_yy, tk_xyzzz_yz, tk_xyzzz_zz, tk_xzzz_xx, tk_xzzz_xz, tk_yzzz_xy, tk_yzzz_y, tk_yzzz_yy, tk_yzzz_yz, tk_yzzz_zz, ts_xyzzz_xx, ts_xyzzz_xy, ts_xyzzz_xz, ts_xyzzz_yy, ts_xyzzz_yz, ts_xyzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzz_xx[i] = tk_xzzz_xx[i] * pa_y[i] + 2.0 * ts_xyzzz_xx[i] * fz_0;

        tk_xyzzz_xy[i] = tk_yzzz_y[i] * fe_0 + tk_yzzz_xy[i] * pa_x[i] + 2.0 * ts_xyzzz_xy[i] * fz_0;

        tk_xyzzz_xz[i] = tk_xzzz_xz[i] * pa_y[i] + 2.0 * ts_xyzzz_xz[i] * fz_0;

        tk_xyzzz_yy[i] = tk_yzzz_yy[i] * pa_x[i] + 2.0 * ts_xyzzz_yy[i] * fz_0;

        tk_xyzzz_yz[i] = tk_yzzz_yz[i] * pa_x[i] + 2.0 * ts_xyzzz_yz[i] * fz_0;

        tk_xyzzz_zz[i] = tk_yzzz_zz[i] * pa_x[i] + 2.0 * ts_xyzzz_zz[i] * fz_0;
    }

    // Set up 84-90 components of targeted buffer : HD

    auto tk_xzzzz_xx = pbuffer.data(idx_kin_hd + 84);

    auto tk_xzzzz_xy = pbuffer.data(idx_kin_hd + 85);

    auto tk_xzzzz_xz = pbuffer.data(idx_kin_hd + 86);

    auto tk_xzzzz_yy = pbuffer.data(idx_kin_hd + 87);

    auto tk_xzzzz_yz = pbuffer.data(idx_kin_hd + 88);

    auto tk_xzzzz_zz = pbuffer.data(idx_kin_hd + 89);

    #pragma omp simd aligned(pa_x, tk_xzzzz_xx, tk_xzzzz_xy, tk_xzzzz_xz, tk_xzzzz_yy, tk_xzzzz_yz, tk_xzzzz_zz, tk_zzzz_x, tk_zzzz_xx, tk_zzzz_xy, tk_zzzz_xz, tk_zzzz_y, tk_zzzz_yy, tk_zzzz_yz, tk_zzzz_z, tk_zzzz_zz, ts_xzzzz_xx, ts_xzzzz_xy, ts_xzzzz_xz, ts_xzzzz_yy, ts_xzzzz_yz, ts_xzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzz_xx[i] = 2.0 * tk_zzzz_x[i] * fe_0 + tk_zzzz_xx[i] * pa_x[i] + 2.0 * ts_xzzzz_xx[i] * fz_0;

        tk_xzzzz_xy[i] = tk_zzzz_y[i] * fe_0 + tk_zzzz_xy[i] * pa_x[i] + 2.0 * ts_xzzzz_xy[i] * fz_0;

        tk_xzzzz_xz[i] = tk_zzzz_z[i] * fe_0 + tk_zzzz_xz[i] * pa_x[i] + 2.0 * ts_xzzzz_xz[i] * fz_0;

        tk_xzzzz_yy[i] = tk_zzzz_yy[i] * pa_x[i] + 2.0 * ts_xzzzz_yy[i] * fz_0;

        tk_xzzzz_yz[i] = tk_zzzz_yz[i] * pa_x[i] + 2.0 * ts_xzzzz_yz[i] * fz_0;

        tk_xzzzz_zz[i] = tk_zzzz_zz[i] * pa_x[i] + 2.0 * ts_xzzzz_zz[i] * fz_0;
    }

    // Set up 90-96 components of targeted buffer : HD

    auto tk_yyyyy_xx = pbuffer.data(idx_kin_hd + 90);

    auto tk_yyyyy_xy = pbuffer.data(idx_kin_hd + 91);

    auto tk_yyyyy_xz = pbuffer.data(idx_kin_hd + 92);

    auto tk_yyyyy_yy = pbuffer.data(idx_kin_hd + 93);

    auto tk_yyyyy_yz = pbuffer.data(idx_kin_hd + 94);

    auto tk_yyyyy_zz = pbuffer.data(idx_kin_hd + 95);

    #pragma omp simd aligned(pa_y, tk_yyy_xx, tk_yyy_xy, tk_yyy_xz, tk_yyy_yy, tk_yyy_yz, tk_yyy_zz, tk_yyyy_x, tk_yyyy_xx, tk_yyyy_xy, tk_yyyy_xz, tk_yyyy_y, tk_yyyy_yy, tk_yyyy_yz, tk_yyyy_z, tk_yyyy_zz, tk_yyyyy_xx, tk_yyyyy_xy, tk_yyyyy_xz, tk_yyyyy_yy, tk_yyyyy_yz, tk_yyyyy_zz, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyy_xx[i] = -8.0 * ts_yyy_xx[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xx[i] * fe_0 + tk_yyyy_xx[i] * pa_y[i] + 2.0 * ts_yyyyy_xx[i] * fz_0;

        tk_yyyyy_xy[i] = -8.0 * ts_yyy_xy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xy[i] * fe_0 + tk_yyyy_x[i] * fe_0 + tk_yyyy_xy[i] * pa_y[i] + 2.0 * ts_yyyyy_xy[i] * fz_0;

        tk_yyyyy_xz[i] = -8.0 * ts_yyy_xz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_xz[i] * fe_0 + tk_yyyy_xz[i] * pa_y[i] + 2.0 * ts_yyyyy_xz[i] * fz_0;

        tk_yyyyy_yy[i] = -8.0 * ts_yyy_yy[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yy[i] * fe_0 + 2.0 * tk_yyyy_y[i] * fe_0 + tk_yyyy_yy[i] * pa_y[i] + 2.0 * ts_yyyyy_yy[i] * fz_0;

        tk_yyyyy_yz[i] = -8.0 * ts_yyy_yz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_yz[i] * fe_0 + tk_yyyy_z[i] * fe_0 + tk_yyyy_yz[i] * pa_y[i] + 2.0 * ts_yyyyy_yz[i] * fz_0;

        tk_yyyyy_zz[i] = -8.0 * ts_yyy_zz[i] * fbe_0 * fz_0 + 4.0 * tk_yyy_zz[i] * fe_0 + tk_yyyy_zz[i] * pa_y[i] + 2.0 * ts_yyyyy_zz[i] * fz_0;
    }

    // Set up 96-102 components of targeted buffer : HD

    auto tk_yyyyz_xx = pbuffer.data(idx_kin_hd + 96);

    auto tk_yyyyz_xy = pbuffer.data(idx_kin_hd + 97);

    auto tk_yyyyz_xz = pbuffer.data(idx_kin_hd + 98);

    auto tk_yyyyz_yy = pbuffer.data(idx_kin_hd + 99);

    auto tk_yyyyz_yz = pbuffer.data(idx_kin_hd + 100);

    auto tk_yyyyz_zz = pbuffer.data(idx_kin_hd + 101);

    #pragma omp simd aligned(pa_z, tk_yyyy_x, tk_yyyy_xx, tk_yyyy_xy, tk_yyyy_xz, tk_yyyy_y, tk_yyyy_yy, tk_yyyy_yz, tk_yyyy_z, tk_yyyy_zz, tk_yyyyz_xx, tk_yyyyz_xy, tk_yyyyz_xz, tk_yyyyz_yy, tk_yyyyz_yz, tk_yyyyz_zz, ts_yyyyz_xx, ts_yyyyz_xy, ts_yyyyz_xz, ts_yyyyz_yy, ts_yyyyz_yz, ts_yyyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyz_xx[i] = tk_yyyy_xx[i] * pa_z[i] + 2.0 * ts_yyyyz_xx[i] * fz_0;

        tk_yyyyz_xy[i] = tk_yyyy_xy[i] * pa_z[i] + 2.0 * ts_yyyyz_xy[i] * fz_0;

        tk_yyyyz_xz[i] = tk_yyyy_x[i] * fe_0 + tk_yyyy_xz[i] * pa_z[i] + 2.0 * ts_yyyyz_xz[i] * fz_0;

        tk_yyyyz_yy[i] = tk_yyyy_yy[i] * pa_z[i] + 2.0 * ts_yyyyz_yy[i] * fz_0;

        tk_yyyyz_yz[i] = tk_yyyy_y[i] * fe_0 + tk_yyyy_yz[i] * pa_z[i] + 2.0 * ts_yyyyz_yz[i] * fz_0;

        tk_yyyyz_zz[i] = 2.0 * tk_yyyy_z[i] * fe_0 + tk_yyyy_zz[i] * pa_z[i] + 2.0 * ts_yyyyz_zz[i] * fz_0;
    }

    // Set up 102-108 components of targeted buffer : HD

    auto tk_yyyzz_xx = pbuffer.data(idx_kin_hd + 102);

    auto tk_yyyzz_xy = pbuffer.data(idx_kin_hd + 103);

    auto tk_yyyzz_xz = pbuffer.data(idx_kin_hd + 104);

    auto tk_yyyzz_yy = pbuffer.data(idx_kin_hd + 105);

    auto tk_yyyzz_yz = pbuffer.data(idx_kin_hd + 106);

    auto tk_yyyzz_zz = pbuffer.data(idx_kin_hd + 107);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyy_xy, tk_yyy_yy, tk_yyyz_xy, tk_yyyz_yy, tk_yyyzz_xx, tk_yyyzz_xy, tk_yyyzz_xz, tk_yyyzz_yy, tk_yyyzz_yz, tk_yyyzz_zz, tk_yyzz_xx, tk_yyzz_xz, tk_yyzz_yz, tk_yyzz_z, tk_yyzz_zz, tk_yzz_xx, tk_yzz_xz, tk_yzz_yz, tk_yzz_zz, ts_yyy_xy, ts_yyy_yy, ts_yyyzz_xx, ts_yyyzz_xy, ts_yyyzz_xz, ts_yyyzz_yy, ts_yyyzz_yz, ts_yyyzz_zz, ts_yzz_xx, ts_yzz_xz, ts_yzz_yz, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzz_xx[i] = -4.0 * ts_yzz_xx[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xx[i] * fe_0 + tk_yyzz_xx[i] * pa_y[i] + 2.0 * ts_yyyzz_xx[i] * fz_0;

        tk_yyyzz_xy[i] = -2.0 * ts_yyy_xy[i] * fbe_0 * fz_0 + tk_yyy_xy[i] * fe_0 + tk_yyyz_xy[i] * pa_z[i] + 2.0 * ts_yyyzz_xy[i] * fz_0;

        tk_yyyzz_xz[i] = -4.0 * ts_yzz_xz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_xz[i] * fe_0 + tk_yyzz_xz[i] * pa_y[i] + 2.0 * ts_yyyzz_xz[i] * fz_0;

        tk_yyyzz_yy[i] = -2.0 * ts_yyy_yy[i] * fbe_0 * fz_0 + tk_yyy_yy[i] * fe_0 + tk_yyyz_yy[i] * pa_z[i] + 2.0 * ts_yyyzz_yy[i] * fz_0;

        tk_yyyzz_yz[i] = -4.0 * ts_yzz_yz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_yz[i] * fe_0 + tk_yyzz_z[i] * fe_0 + tk_yyzz_yz[i] * pa_y[i] + 2.0 * ts_yyyzz_yz[i] * fz_0;

        tk_yyyzz_zz[i] = -4.0 * ts_yzz_zz[i] * fbe_0 * fz_0 + 2.0 * tk_yzz_zz[i] * fe_0 + tk_yyzz_zz[i] * pa_y[i] + 2.0 * ts_yyyzz_zz[i] * fz_0;
    }

    // Set up 108-114 components of targeted buffer : HD

    auto tk_yyzzz_xx = pbuffer.data(idx_kin_hd + 108);

    auto tk_yyzzz_xy = pbuffer.data(idx_kin_hd + 109);

    auto tk_yyzzz_xz = pbuffer.data(idx_kin_hd + 110);

    auto tk_yyzzz_yy = pbuffer.data(idx_kin_hd + 111);

    auto tk_yyzzz_yz = pbuffer.data(idx_kin_hd + 112);

    auto tk_yyzzz_zz = pbuffer.data(idx_kin_hd + 113);

    #pragma omp simd aligned(pa_y, pa_z, tk_yyz_xy, tk_yyz_yy, tk_yyzz_xy, tk_yyzz_yy, tk_yyzzz_xx, tk_yyzzz_xy, tk_yyzzz_xz, tk_yyzzz_yy, tk_yyzzz_yz, tk_yyzzz_zz, tk_yzzz_xx, tk_yzzz_xz, tk_yzzz_yz, tk_yzzz_z, tk_yzzz_zz, tk_zzz_xx, tk_zzz_xz, tk_zzz_yz, tk_zzz_zz, ts_yyz_xy, ts_yyz_yy, ts_yyzzz_xx, ts_yyzzz_xy, ts_yyzzz_xz, ts_yyzzz_yy, ts_yyzzz_yz, ts_yyzzz_zz, ts_zzz_xx, ts_zzz_xz, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzz_xx[i] = -2.0 * ts_zzz_xx[i] * fbe_0 * fz_0 + tk_zzz_xx[i] * fe_0 + tk_yzzz_xx[i] * pa_y[i] + 2.0 * ts_yyzzz_xx[i] * fz_0;

        tk_yyzzz_xy[i] = -4.0 * ts_yyz_xy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_xy[i] * fe_0 + tk_yyzz_xy[i] * pa_z[i] + 2.0 * ts_yyzzz_xy[i] * fz_0;

        tk_yyzzz_xz[i] = -2.0 * ts_zzz_xz[i] * fbe_0 * fz_0 + tk_zzz_xz[i] * fe_0 + tk_yzzz_xz[i] * pa_y[i] + 2.0 * ts_yyzzz_xz[i] * fz_0;

        tk_yyzzz_yy[i] = -4.0 * ts_yyz_yy[i] * fbe_0 * fz_0 + 2.0 * tk_yyz_yy[i] * fe_0 + tk_yyzz_yy[i] * pa_z[i] + 2.0 * ts_yyzzz_yy[i] * fz_0;

        tk_yyzzz_yz[i] = -2.0 * ts_zzz_yz[i] * fbe_0 * fz_0 + tk_zzz_yz[i] * fe_0 + tk_yzzz_z[i] * fe_0 + tk_yzzz_yz[i] * pa_y[i] + 2.0 * ts_yyzzz_yz[i] * fz_0;

        tk_yyzzz_zz[i] = -2.0 * ts_zzz_zz[i] * fbe_0 * fz_0 + tk_zzz_zz[i] * fe_0 + tk_yzzz_zz[i] * pa_y[i] + 2.0 * ts_yyzzz_zz[i] * fz_0;
    }

    // Set up 114-120 components of targeted buffer : HD

    auto tk_yzzzz_xx = pbuffer.data(idx_kin_hd + 114);

    auto tk_yzzzz_xy = pbuffer.data(idx_kin_hd + 115);

    auto tk_yzzzz_xz = pbuffer.data(idx_kin_hd + 116);

    auto tk_yzzzz_yy = pbuffer.data(idx_kin_hd + 117);

    auto tk_yzzzz_yz = pbuffer.data(idx_kin_hd + 118);

    auto tk_yzzzz_zz = pbuffer.data(idx_kin_hd + 119);

    #pragma omp simd aligned(pa_y, tk_yzzzz_xx, tk_yzzzz_xy, tk_yzzzz_xz, tk_yzzzz_yy, tk_yzzzz_yz, tk_yzzzz_zz, tk_zzzz_x, tk_zzzz_xx, tk_zzzz_xy, tk_zzzz_xz, tk_zzzz_y, tk_zzzz_yy, tk_zzzz_yz, tk_zzzz_z, tk_zzzz_zz, ts_yzzzz_xx, ts_yzzzz_xy, ts_yzzzz_xz, ts_yzzzz_yy, ts_yzzzz_yz, ts_yzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzz_xx[i] = tk_zzzz_xx[i] * pa_y[i] + 2.0 * ts_yzzzz_xx[i] * fz_0;

        tk_yzzzz_xy[i] = tk_zzzz_x[i] * fe_0 + tk_zzzz_xy[i] * pa_y[i] + 2.0 * ts_yzzzz_xy[i] * fz_0;

        tk_yzzzz_xz[i] = tk_zzzz_xz[i] * pa_y[i] + 2.0 * ts_yzzzz_xz[i] * fz_0;

        tk_yzzzz_yy[i] = 2.0 * tk_zzzz_y[i] * fe_0 + tk_zzzz_yy[i] * pa_y[i] + 2.0 * ts_yzzzz_yy[i] * fz_0;

        tk_yzzzz_yz[i] = tk_zzzz_z[i] * fe_0 + tk_zzzz_yz[i] * pa_y[i] + 2.0 * ts_yzzzz_yz[i] * fz_0;

        tk_yzzzz_zz[i] = tk_zzzz_zz[i] * pa_y[i] + 2.0 * ts_yzzzz_zz[i] * fz_0;
    }

    // Set up 120-126 components of targeted buffer : HD

    auto tk_zzzzz_xx = pbuffer.data(idx_kin_hd + 120);

    auto tk_zzzzz_xy = pbuffer.data(idx_kin_hd + 121);

    auto tk_zzzzz_xz = pbuffer.data(idx_kin_hd + 122);

    auto tk_zzzzz_yy = pbuffer.data(idx_kin_hd + 123);

    auto tk_zzzzz_yz = pbuffer.data(idx_kin_hd + 124);

    auto tk_zzzzz_zz = pbuffer.data(idx_kin_hd + 125);

    #pragma omp simd aligned(pa_z, tk_zzz_xx, tk_zzz_xy, tk_zzz_xz, tk_zzz_yy, tk_zzz_yz, tk_zzz_zz, tk_zzzz_x, tk_zzzz_xx, tk_zzzz_xy, tk_zzzz_xz, tk_zzzz_y, tk_zzzz_yy, tk_zzzz_yz, tk_zzzz_z, tk_zzzz_zz, tk_zzzzz_xx, tk_zzzzz_xy, tk_zzzzz_xz, tk_zzzzz_yy, tk_zzzzz_yz, tk_zzzzz_zz, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzz_xx[i] = -8.0 * ts_zzz_xx[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xx[i] * fe_0 + tk_zzzz_xx[i] * pa_z[i] + 2.0 * ts_zzzzz_xx[i] * fz_0;

        tk_zzzzz_xy[i] = -8.0 * ts_zzz_xy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xy[i] * fe_0 + tk_zzzz_xy[i] * pa_z[i] + 2.0 * ts_zzzzz_xy[i] * fz_0;

        tk_zzzzz_xz[i] = -8.0 * ts_zzz_xz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_xz[i] * fe_0 + tk_zzzz_x[i] * fe_0 + tk_zzzz_xz[i] * pa_z[i] + 2.0 * ts_zzzzz_xz[i] * fz_0;

        tk_zzzzz_yy[i] = -8.0 * ts_zzz_yy[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yy[i] * fe_0 + tk_zzzz_yy[i] * pa_z[i] + 2.0 * ts_zzzzz_yy[i] * fz_0;

        tk_zzzzz_yz[i] = -8.0 * ts_zzz_yz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_yz[i] * fe_0 + tk_zzzz_y[i] * fe_0 + tk_zzzz_yz[i] * pa_z[i] + 2.0 * ts_zzzzz_yz[i] * fz_0;

        tk_zzzzz_zz[i] = -8.0 * ts_zzz_zz[i] * fbe_0 * fz_0 + 4.0 * tk_zzz_zz[i] * fe_0 + 2.0 * tk_zzzz_z[i] * fe_0 + tk_zzzz_zz[i] * pa_z[i] + 2.0 * ts_zzzzz_zz[i] * fz_0;
    }

}

} // kinrec namespace

