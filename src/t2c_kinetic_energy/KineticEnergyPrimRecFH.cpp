#include "KineticEnergyPrimRecFH.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_fh(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_fh,
                            const size_t idx_ovl_ph,
                            const size_t idx_kin_ph,
                            const size_t idx_kin_dg,
                            const size_t idx_kin_dh,
                            const size_t idx_ovl_fh,
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

    // Set up components of auxiliary buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ovl_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ovl_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ovl_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ovl_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ovl_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ovl_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ovl_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ovl_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ovl_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ovl_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ovl_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ovl_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ovl_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ovl_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ovl_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ovl_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ovl_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ovl_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ovl_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ovl_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ovl_ph + 20);

    auto ts_y_xxxxx = pbuffer.data(idx_ovl_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ovl_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ovl_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ovl_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ovl_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ovl_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ovl_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ovl_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ovl_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ovl_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ovl_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ovl_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ovl_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ovl_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ovl_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ovl_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ovl_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ovl_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ovl_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ovl_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ovl_ph + 41);

    auto ts_z_xxxxx = pbuffer.data(idx_ovl_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ovl_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ovl_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ovl_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ovl_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ovl_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ovl_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ovl_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ovl_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ovl_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ovl_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ovl_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ovl_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ovl_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ovl_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ovl_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ovl_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ovl_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ovl_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ovl_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ovl_ph + 62);

    // Set up components of auxiliary buffer : PH

    auto tk_x_xxxxx = pbuffer.data(idx_kin_ph);

    auto tk_x_xxxxy = pbuffer.data(idx_kin_ph + 1);

    auto tk_x_xxxxz = pbuffer.data(idx_kin_ph + 2);

    auto tk_x_xxxyy = pbuffer.data(idx_kin_ph + 3);

    auto tk_x_xxxyz = pbuffer.data(idx_kin_ph + 4);

    auto tk_x_xxxzz = pbuffer.data(idx_kin_ph + 5);

    auto tk_x_xxyyy = pbuffer.data(idx_kin_ph + 6);

    auto tk_x_xxyyz = pbuffer.data(idx_kin_ph + 7);

    auto tk_x_xxyzz = pbuffer.data(idx_kin_ph + 8);

    auto tk_x_xxzzz = pbuffer.data(idx_kin_ph + 9);

    auto tk_x_xyyyy = pbuffer.data(idx_kin_ph + 10);

    auto tk_x_xyyyz = pbuffer.data(idx_kin_ph + 11);

    auto tk_x_xyyzz = pbuffer.data(idx_kin_ph + 12);

    auto tk_x_xyzzz = pbuffer.data(idx_kin_ph + 13);

    auto tk_x_xzzzz = pbuffer.data(idx_kin_ph + 14);

    auto tk_x_yyyyy = pbuffer.data(idx_kin_ph + 15);

    auto tk_x_yyyyz = pbuffer.data(idx_kin_ph + 16);

    auto tk_x_yyyzz = pbuffer.data(idx_kin_ph + 17);

    auto tk_x_yyzzz = pbuffer.data(idx_kin_ph + 18);

    auto tk_x_yzzzz = pbuffer.data(idx_kin_ph + 19);

    auto tk_x_zzzzz = pbuffer.data(idx_kin_ph + 20);

    auto tk_y_xxxxx = pbuffer.data(idx_kin_ph + 21);

    auto tk_y_xxxxy = pbuffer.data(idx_kin_ph + 22);

    auto tk_y_xxxxz = pbuffer.data(idx_kin_ph + 23);

    auto tk_y_xxxyy = pbuffer.data(idx_kin_ph + 24);

    auto tk_y_xxxyz = pbuffer.data(idx_kin_ph + 25);

    auto tk_y_xxxzz = pbuffer.data(idx_kin_ph + 26);

    auto tk_y_xxyyy = pbuffer.data(idx_kin_ph + 27);

    auto tk_y_xxyyz = pbuffer.data(idx_kin_ph + 28);

    auto tk_y_xxyzz = pbuffer.data(idx_kin_ph + 29);

    auto tk_y_xxzzz = pbuffer.data(idx_kin_ph + 30);

    auto tk_y_xyyyy = pbuffer.data(idx_kin_ph + 31);

    auto tk_y_xyyyz = pbuffer.data(idx_kin_ph + 32);

    auto tk_y_xyyzz = pbuffer.data(idx_kin_ph + 33);

    auto tk_y_xyzzz = pbuffer.data(idx_kin_ph + 34);

    auto tk_y_xzzzz = pbuffer.data(idx_kin_ph + 35);

    auto tk_y_yyyyy = pbuffer.data(idx_kin_ph + 36);

    auto tk_y_yyyyz = pbuffer.data(idx_kin_ph + 37);

    auto tk_y_yyyzz = pbuffer.data(idx_kin_ph + 38);

    auto tk_y_yyzzz = pbuffer.data(idx_kin_ph + 39);

    auto tk_y_yzzzz = pbuffer.data(idx_kin_ph + 40);

    auto tk_y_zzzzz = pbuffer.data(idx_kin_ph + 41);

    auto tk_z_xxxxx = pbuffer.data(idx_kin_ph + 42);

    auto tk_z_xxxxy = pbuffer.data(idx_kin_ph + 43);

    auto tk_z_xxxxz = pbuffer.data(idx_kin_ph + 44);

    auto tk_z_xxxyy = pbuffer.data(idx_kin_ph + 45);

    auto tk_z_xxxyz = pbuffer.data(idx_kin_ph + 46);

    auto tk_z_xxxzz = pbuffer.data(idx_kin_ph + 47);

    auto tk_z_xxyyy = pbuffer.data(idx_kin_ph + 48);

    auto tk_z_xxyyz = pbuffer.data(idx_kin_ph + 49);

    auto tk_z_xxyzz = pbuffer.data(idx_kin_ph + 50);

    auto tk_z_xxzzz = pbuffer.data(idx_kin_ph + 51);

    auto tk_z_xyyyy = pbuffer.data(idx_kin_ph + 52);

    auto tk_z_xyyyz = pbuffer.data(idx_kin_ph + 53);

    auto tk_z_xyyzz = pbuffer.data(idx_kin_ph + 54);

    auto tk_z_xyzzz = pbuffer.data(idx_kin_ph + 55);

    auto tk_z_xzzzz = pbuffer.data(idx_kin_ph + 56);

    auto tk_z_yyyyy = pbuffer.data(idx_kin_ph + 57);

    auto tk_z_yyyyz = pbuffer.data(idx_kin_ph + 58);

    auto tk_z_yyyzz = pbuffer.data(idx_kin_ph + 59);

    auto tk_z_yyzzz = pbuffer.data(idx_kin_ph + 60);

    auto tk_z_yzzzz = pbuffer.data(idx_kin_ph + 61);

    auto tk_z_zzzzz = pbuffer.data(idx_kin_ph + 62);

    // Set up components of auxiliary buffer : DG

    auto tk_xx_xxxx = pbuffer.data(idx_kin_dg);

    auto tk_xx_xxxy = pbuffer.data(idx_kin_dg + 1);

    auto tk_xx_xxxz = pbuffer.data(idx_kin_dg + 2);

    auto tk_xx_xxyy = pbuffer.data(idx_kin_dg + 3);

    auto tk_xx_xxyz = pbuffer.data(idx_kin_dg + 4);

    auto tk_xx_xxzz = pbuffer.data(idx_kin_dg + 5);

    auto tk_xx_xyyy = pbuffer.data(idx_kin_dg + 6);

    auto tk_xx_xyyz = pbuffer.data(idx_kin_dg + 7);

    auto tk_xx_xyzz = pbuffer.data(idx_kin_dg + 8);

    auto tk_xx_xzzz = pbuffer.data(idx_kin_dg + 9);

    auto tk_xx_yyyy = pbuffer.data(idx_kin_dg + 10);

    auto tk_xx_yyyz = pbuffer.data(idx_kin_dg + 11);

    auto tk_xx_yyzz = pbuffer.data(idx_kin_dg + 12);

    auto tk_xx_yzzz = pbuffer.data(idx_kin_dg + 13);

    auto tk_xx_zzzz = pbuffer.data(idx_kin_dg + 14);

    auto tk_yy_xxxx = pbuffer.data(idx_kin_dg + 45);

    auto tk_yy_xxxy = pbuffer.data(idx_kin_dg + 46);

    auto tk_yy_xxxz = pbuffer.data(idx_kin_dg + 47);

    auto tk_yy_xxyy = pbuffer.data(idx_kin_dg + 48);

    auto tk_yy_xxyz = pbuffer.data(idx_kin_dg + 49);

    auto tk_yy_xxzz = pbuffer.data(idx_kin_dg + 50);

    auto tk_yy_xyyy = pbuffer.data(idx_kin_dg + 51);

    auto tk_yy_xyyz = pbuffer.data(idx_kin_dg + 52);

    auto tk_yy_xyzz = pbuffer.data(idx_kin_dg + 53);

    auto tk_yy_xzzz = pbuffer.data(idx_kin_dg + 54);

    auto tk_yy_yyyy = pbuffer.data(idx_kin_dg + 55);

    auto tk_yy_yyyz = pbuffer.data(idx_kin_dg + 56);

    auto tk_yy_yyzz = pbuffer.data(idx_kin_dg + 57);

    auto tk_yy_yzzz = pbuffer.data(idx_kin_dg + 58);

    auto tk_yy_zzzz = pbuffer.data(idx_kin_dg + 59);

    auto tk_yz_xxyz = pbuffer.data(idx_kin_dg + 64);

    auto tk_yz_xyyz = pbuffer.data(idx_kin_dg + 67);

    auto tk_yz_xyzz = pbuffer.data(idx_kin_dg + 68);

    auto tk_yz_yyyz = pbuffer.data(idx_kin_dg + 71);

    auto tk_yz_yyzz = pbuffer.data(idx_kin_dg + 72);

    auto tk_yz_yzzz = pbuffer.data(idx_kin_dg + 73);

    auto tk_zz_xxxx = pbuffer.data(idx_kin_dg + 75);

    auto tk_zz_xxxy = pbuffer.data(idx_kin_dg + 76);

    auto tk_zz_xxxz = pbuffer.data(idx_kin_dg + 77);

    auto tk_zz_xxyy = pbuffer.data(idx_kin_dg + 78);

    auto tk_zz_xxyz = pbuffer.data(idx_kin_dg + 79);

    auto tk_zz_xxzz = pbuffer.data(idx_kin_dg + 80);

    auto tk_zz_xyyy = pbuffer.data(idx_kin_dg + 81);

    auto tk_zz_xyyz = pbuffer.data(idx_kin_dg + 82);

    auto tk_zz_xyzz = pbuffer.data(idx_kin_dg + 83);

    auto tk_zz_xzzz = pbuffer.data(idx_kin_dg + 84);

    auto tk_zz_yyyy = pbuffer.data(idx_kin_dg + 85);

    auto tk_zz_yyyz = pbuffer.data(idx_kin_dg + 86);

    auto tk_zz_yyzz = pbuffer.data(idx_kin_dg + 87);

    auto tk_zz_yzzz = pbuffer.data(idx_kin_dg + 88);

    auto tk_zz_zzzz = pbuffer.data(idx_kin_dg + 89);

    // Set up components of auxiliary buffer : DH

    auto tk_xx_xxxxx = pbuffer.data(idx_kin_dh);

    auto tk_xx_xxxxy = pbuffer.data(idx_kin_dh + 1);

    auto tk_xx_xxxxz = pbuffer.data(idx_kin_dh + 2);

    auto tk_xx_xxxyy = pbuffer.data(idx_kin_dh + 3);

    auto tk_xx_xxxyz = pbuffer.data(idx_kin_dh + 4);

    auto tk_xx_xxxzz = pbuffer.data(idx_kin_dh + 5);

    auto tk_xx_xxyyy = pbuffer.data(idx_kin_dh + 6);

    auto tk_xx_xxyyz = pbuffer.data(idx_kin_dh + 7);

    auto tk_xx_xxyzz = pbuffer.data(idx_kin_dh + 8);

    auto tk_xx_xxzzz = pbuffer.data(idx_kin_dh + 9);

    auto tk_xx_xyyyy = pbuffer.data(idx_kin_dh + 10);

    auto tk_xx_xyyyz = pbuffer.data(idx_kin_dh + 11);

    auto tk_xx_xyyzz = pbuffer.data(idx_kin_dh + 12);

    auto tk_xx_xyzzz = pbuffer.data(idx_kin_dh + 13);

    auto tk_xx_xzzzz = pbuffer.data(idx_kin_dh + 14);

    auto tk_xx_yyyyy = pbuffer.data(idx_kin_dh + 15);

    auto tk_xx_yyyyz = pbuffer.data(idx_kin_dh + 16);

    auto tk_xx_yyyzz = pbuffer.data(idx_kin_dh + 17);

    auto tk_xx_yyzzz = pbuffer.data(idx_kin_dh + 18);

    auto tk_xx_yzzzz = pbuffer.data(idx_kin_dh + 19);

    auto tk_xx_zzzzz = pbuffer.data(idx_kin_dh + 20);

    auto tk_xy_xxxxy = pbuffer.data(idx_kin_dh + 22);

    auto tk_xy_xxxyy = pbuffer.data(idx_kin_dh + 24);

    auto tk_xy_xxyyy = pbuffer.data(idx_kin_dh + 27);

    auto tk_xy_xyyyy = pbuffer.data(idx_kin_dh + 31);

    auto tk_xz_xxxxx = pbuffer.data(idx_kin_dh + 42);

    auto tk_xz_xxxxz = pbuffer.data(idx_kin_dh + 44);

    auto tk_xz_xxxzz = pbuffer.data(idx_kin_dh + 47);

    auto tk_xz_xxzzz = pbuffer.data(idx_kin_dh + 51);

    auto tk_xz_xzzzz = pbuffer.data(idx_kin_dh + 56);

    auto tk_yy_xxxxx = pbuffer.data(idx_kin_dh + 63);

    auto tk_yy_xxxxy = pbuffer.data(idx_kin_dh + 64);

    auto tk_yy_xxxxz = pbuffer.data(idx_kin_dh + 65);

    auto tk_yy_xxxyy = pbuffer.data(idx_kin_dh + 66);

    auto tk_yy_xxxyz = pbuffer.data(idx_kin_dh + 67);

    auto tk_yy_xxxzz = pbuffer.data(idx_kin_dh + 68);

    auto tk_yy_xxyyy = pbuffer.data(idx_kin_dh + 69);

    auto tk_yy_xxyyz = pbuffer.data(idx_kin_dh + 70);

    auto tk_yy_xxyzz = pbuffer.data(idx_kin_dh + 71);

    auto tk_yy_xxzzz = pbuffer.data(idx_kin_dh + 72);

    auto tk_yy_xyyyy = pbuffer.data(idx_kin_dh + 73);

    auto tk_yy_xyyyz = pbuffer.data(idx_kin_dh + 74);

    auto tk_yy_xyyzz = pbuffer.data(idx_kin_dh + 75);

    auto tk_yy_xyzzz = pbuffer.data(idx_kin_dh + 76);

    auto tk_yy_xzzzz = pbuffer.data(idx_kin_dh + 77);

    auto tk_yy_yyyyy = pbuffer.data(idx_kin_dh + 78);

    auto tk_yy_yyyyz = pbuffer.data(idx_kin_dh + 79);

    auto tk_yy_yyyzz = pbuffer.data(idx_kin_dh + 80);

    auto tk_yy_yyzzz = pbuffer.data(idx_kin_dh + 81);

    auto tk_yy_yzzzz = pbuffer.data(idx_kin_dh + 82);

    auto tk_yy_zzzzz = pbuffer.data(idx_kin_dh + 83);

    auto tk_yz_xxxyz = pbuffer.data(idx_kin_dh + 88);

    auto tk_yz_xxyyz = pbuffer.data(idx_kin_dh + 91);

    auto tk_yz_xxyzz = pbuffer.data(idx_kin_dh + 92);

    auto tk_yz_xyyyz = pbuffer.data(idx_kin_dh + 95);

    auto tk_yz_xyyzz = pbuffer.data(idx_kin_dh + 96);

    auto tk_yz_xyzzz = pbuffer.data(idx_kin_dh + 97);

    auto tk_yz_yyyyy = pbuffer.data(idx_kin_dh + 99);

    auto tk_yz_yyyyz = pbuffer.data(idx_kin_dh + 100);

    auto tk_yz_yyyzz = pbuffer.data(idx_kin_dh + 101);

    auto tk_yz_yyzzz = pbuffer.data(idx_kin_dh + 102);

    auto tk_yz_yzzzz = pbuffer.data(idx_kin_dh + 103);

    auto tk_yz_zzzzz = pbuffer.data(idx_kin_dh + 104);

    auto tk_zz_xxxxx = pbuffer.data(idx_kin_dh + 105);

    auto tk_zz_xxxxy = pbuffer.data(idx_kin_dh + 106);

    auto tk_zz_xxxxz = pbuffer.data(idx_kin_dh + 107);

    auto tk_zz_xxxyy = pbuffer.data(idx_kin_dh + 108);

    auto tk_zz_xxxyz = pbuffer.data(idx_kin_dh + 109);

    auto tk_zz_xxxzz = pbuffer.data(idx_kin_dh + 110);

    auto tk_zz_xxyyy = pbuffer.data(idx_kin_dh + 111);

    auto tk_zz_xxyyz = pbuffer.data(idx_kin_dh + 112);

    auto tk_zz_xxyzz = pbuffer.data(idx_kin_dh + 113);

    auto tk_zz_xxzzz = pbuffer.data(idx_kin_dh + 114);

    auto tk_zz_xyyyy = pbuffer.data(idx_kin_dh + 115);

    auto tk_zz_xyyyz = pbuffer.data(idx_kin_dh + 116);

    auto tk_zz_xyyzz = pbuffer.data(idx_kin_dh + 117);

    auto tk_zz_xyzzz = pbuffer.data(idx_kin_dh + 118);

    auto tk_zz_xzzzz = pbuffer.data(idx_kin_dh + 119);

    auto tk_zz_yyyyy = pbuffer.data(idx_kin_dh + 120);

    auto tk_zz_yyyyz = pbuffer.data(idx_kin_dh + 121);

    auto tk_zz_yyyzz = pbuffer.data(idx_kin_dh + 122);

    auto tk_zz_yyzzz = pbuffer.data(idx_kin_dh + 123);

    auto tk_zz_yzzzz = pbuffer.data(idx_kin_dh + 124);

    auto tk_zz_zzzzz = pbuffer.data(idx_kin_dh + 125);

    // Set up components of auxiliary buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_ovl_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_ovl_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_ovl_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_ovl_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_ovl_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_ovl_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_ovl_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_ovl_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_ovl_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_ovl_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_ovl_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_ovl_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_ovl_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_ovl_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_ovl_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_ovl_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_ovl_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_ovl_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_ovl_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_ovl_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_ovl_fh + 20);

    auto ts_xxy_xxxxx = pbuffer.data(idx_ovl_fh + 21);

    auto ts_xxy_xxxxy = pbuffer.data(idx_ovl_fh + 22);

    auto ts_xxy_xxxxz = pbuffer.data(idx_ovl_fh + 23);

    auto ts_xxy_xxxyy = pbuffer.data(idx_ovl_fh + 24);

    auto ts_xxy_xxxyz = pbuffer.data(idx_ovl_fh + 25);

    auto ts_xxy_xxxzz = pbuffer.data(idx_ovl_fh + 26);

    auto ts_xxy_xxyyy = pbuffer.data(idx_ovl_fh + 27);

    auto ts_xxy_xxyyz = pbuffer.data(idx_ovl_fh + 28);

    auto ts_xxy_xxyzz = pbuffer.data(idx_ovl_fh + 29);

    auto ts_xxy_xxzzz = pbuffer.data(idx_ovl_fh + 30);

    auto ts_xxy_xyyyy = pbuffer.data(idx_ovl_fh + 31);

    auto ts_xxy_xyyyz = pbuffer.data(idx_ovl_fh + 32);

    auto ts_xxy_xyyzz = pbuffer.data(idx_ovl_fh + 33);

    auto ts_xxy_xyzzz = pbuffer.data(idx_ovl_fh + 34);

    auto ts_xxy_xzzzz = pbuffer.data(idx_ovl_fh + 35);

    auto ts_xxy_yyyyy = pbuffer.data(idx_ovl_fh + 36);

    auto ts_xxy_yyyyz = pbuffer.data(idx_ovl_fh + 37);

    auto ts_xxy_yyyzz = pbuffer.data(idx_ovl_fh + 38);

    auto ts_xxy_yyzzz = pbuffer.data(idx_ovl_fh + 39);

    auto ts_xxy_yzzzz = pbuffer.data(idx_ovl_fh + 40);

    auto ts_xxy_zzzzz = pbuffer.data(idx_ovl_fh + 41);

    auto ts_xxz_xxxxx = pbuffer.data(idx_ovl_fh + 42);

    auto ts_xxz_xxxxy = pbuffer.data(idx_ovl_fh + 43);

    auto ts_xxz_xxxxz = pbuffer.data(idx_ovl_fh + 44);

    auto ts_xxz_xxxyy = pbuffer.data(idx_ovl_fh + 45);

    auto ts_xxz_xxxyz = pbuffer.data(idx_ovl_fh + 46);

    auto ts_xxz_xxxzz = pbuffer.data(idx_ovl_fh + 47);

    auto ts_xxz_xxyyy = pbuffer.data(idx_ovl_fh + 48);

    auto ts_xxz_xxyyz = pbuffer.data(idx_ovl_fh + 49);

    auto ts_xxz_xxyzz = pbuffer.data(idx_ovl_fh + 50);

    auto ts_xxz_xxzzz = pbuffer.data(idx_ovl_fh + 51);

    auto ts_xxz_xyyyy = pbuffer.data(idx_ovl_fh + 52);

    auto ts_xxz_xyyyz = pbuffer.data(idx_ovl_fh + 53);

    auto ts_xxz_xyyzz = pbuffer.data(idx_ovl_fh + 54);

    auto ts_xxz_xyzzz = pbuffer.data(idx_ovl_fh + 55);

    auto ts_xxz_xzzzz = pbuffer.data(idx_ovl_fh + 56);

    auto ts_xxz_yyyyy = pbuffer.data(idx_ovl_fh + 57);

    auto ts_xxz_yyyyz = pbuffer.data(idx_ovl_fh + 58);

    auto ts_xxz_yyyzz = pbuffer.data(idx_ovl_fh + 59);

    auto ts_xxz_yyzzz = pbuffer.data(idx_ovl_fh + 60);

    auto ts_xxz_yzzzz = pbuffer.data(idx_ovl_fh + 61);

    auto ts_xxz_zzzzz = pbuffer.data(idx_ovl_fh + 62);

    auto ts_xyy_xxxxx = pbuffer.data(idx_ovl_fh + 63);

    auto ts_xyy_xxxxy = pbuffer.data(idx_ovl_fh + 64);

    auto ts_xyy_xxxxz = pbuffer.data(idx_ovl_fh + 65);

    auto ts_xyy_xxxyy = pbuffer.data(idx_ovl_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_ovl_fh + 67);

    auto ts_xyy_xxxzz = pbuffer.data(idx_ovl_fh + 68);

    auto ts_xyy_xxyyy = pbuffer.data(idx_ovl_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_ovl_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_ovl_fh + 71);

    auto ts_xyy_xxzzz = pbuffer.data(idx_ovl_fh + 72);

    auto ts_xyy_xyyyy = pbuffer.data(idx_ovl_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_ovl_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_ovl_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_ovl_fh + 76);

    auto ts_xyy_xzzzz = pbuffer.data(idx_ovl_fh + 77);

    auto ts_xyy_yyyyy = pbuffer.data(idx_ovl_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_ovl_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_ovl_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_ovl_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_ovl_fh + 82);

    auto ts_xyy_zzzzz = pbuffer.data(idx_ovl_fh + 83);

    auto ts_xyz_xxxxx = pbuffer.data(idx_ovl_fh + 84);

    auto ts_xyz_xxxxy = pbuffer.data(idx_ovl_fh + 85);

    auto ts_xyz_xxxxz = pbuffer.data(idx_ovl_fh + 86);

    auto ts_xyz_xxxyy = pbuffer.data(idx_ovl_fh + 87);

    auto ts_xyz_xxxyz = pbuffer.data(idx_ovl_fh + 88);

    auto ts_xyz_xxxzz = pbuffer.data(idx_ovl_fh + 89);

    auto ts_xyz_xxyyy = pbuffer.data(idx_ovl_fh + 90);

    auto ts_xyz_xxyyz = pbuffer.data(idx_ovl_fh + 91);

    auto ts_xyz_xxyzz = pbuffer.data(idx_ovl_fh + 92);

    auto ts_xyz_xxzzz = pbuffer.data(idx_ovl_fh + 93);

    auto ts_xyz_xyyyy = pbuffer.data(idx_ovl_fh + 94);

    auto ts_xyz_xyyyz = pbuffer.data(idx_ovl_fh + 95);

    auto ts_xyz_xyyzz = pbuffer.data(idx_ovl_fh + 96);

    auto ts_xyz_xyzzz = pbuffer.data(idx_ovl_fh + 97);

    auto ts_xyz_xzzzz = pbuffer.data(idx_ovl_fh + 98);

    auto ts_xyz_yyyyy = pbuffer.data(idx_ovl_fh + 99);

    auto ts_xyz_yyyyz = pbuffer.data(idx_ovl_fh + 100);

    auto ts_xyz_yyyzz = pbuffer.data(idx_ovl_fh + 101);

    auto ts_xyz_yyzzz = pbuffer.data(idx_ovl_fh + 102);

    auto ts_xyz_yzzzz = pbuffer.data(idx_ovl_fh + 103);

    auto ts_xyz_zzzzz = pbuffer.data(idx_ovl_fh + 104);

    auto ts_xzz_xxxxx = pbuffer.data(idx_ovl_fh + 105);

    auto ts_xzz_xxxxy = pbuffer.data(idx_ovl_fh + 106);

    auto ts_xzz_xxxxz = pbuffer.data(idx_ovl_fh + 107);

    auto ts_xzz_xxxyy = pbuffer.data(idx_ovl_fh + 108);

    auto ts_xzz_xxxyz = pbuffer.data(idx_ovl_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_ovl_fh + 110);

    auto ts_xzz_xxyyy = pbuffer.data(idx_ovl_fh + 111);

    auto ts_xzz_xxyyz = pbuffer.data(idx_ovl_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_ovl_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_ovl_fh + 114);

    auto ts_xzz_xyyyy = pbuffer.data(idx_ovl_fh + 115);

    auto ts_xzz_xyyyz = pbuffer.data(idx_ovl_fh + 116);

    auto ts_xzz_xyyzz = pbuffer.data(idx_ovl_fh + 117);

    auto ts_xzz_xyzzz = pbuffer.data(idx_ovl_fh + 118);

    auto ts_xzz_xzzzz = pbuffer.data(idx_ovl_fh + 119);

    auto ts_xzz_yyyyy = pbuffer.data(idx_ovl_fh + 120);

    auto ts_xzz_yyyyz = pbuffer.data(idx_ovl_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_ovl_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_ovl_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_ovl_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_ovl_fh + 125);

    auto ts_yyy_xxxxx = pbuffer.data(idx_ovl_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_ovl_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_ovl_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_ovl_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_ovl_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_ovl_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_ovl_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_ovl_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_ovl_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_ovl_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_ovl_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_ovl_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_ovl_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_ovl_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_ovl_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_ovl_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_ovl_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_ovl_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_ovl_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_ovl_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_ovl_fh + 146);

    auto ts_yyz_xxxxx = pbuffer.data(idx_ovl_fh + 147);

    auto ts_yyz_xxxxy = pbuffer.data(idx_ovl_fh + 148);

    auto ts_yyz_xxxxz = pbuffer.data(idx_ovl_fh + 149);

    auto ts_yyz_xxxyy = pbuffer.data(idx_ovl_fh + 150);

    auto ts_yyz_xxxyz = pbuffer.data(idx_ovl_fh + 151);

    auto ts_yyz_xxxzz = pbuffer.data(idx_ovl_fh + 152);

    auto ts_yyz_xxyyy = pbuffer.data(idx_ovl_fh + 153);

    auto ts_yyz_xxyyz = pbuffer.data(idx_ovl_fh + 154);

    auto ts_yyz_xxyzz = pbuffer.data(idx_ovl_fh + 155);

    auto ts_yyz_xxzzz = pbuffer.data(idx_ovl_fh + 156);

    auto ts_yyz_xyyyy = pbuffer.data(idx_ovl_fh + 157);

    auto ts_yyz_xyyyz = pbuffer.data(idx_ovl_fh + 158);

    auto ts_yyz_xyyzz = pbuffer.data(idx_ovl_fh + 159);

    auto ts_yyz_xyzzz = pbuffer.data(idx_ovl_fh + 160);

    auto ts_yyz_xzzzz = pbuffer.data(idx_ovl_fh + 161);

    auto ts_yyz_yyyyy = pbuffer.data(idx_ovl_fh + 162);

    auto ts_yyz_yyyyz = pbuffer.data(idx_ovl_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_ovl_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_ovl_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_ovl_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_ovl_fh + 167);

    auto ts_yzz_xxxxx = pbuffer.data(idx_ovl_fh + 168);

    auto ts_yzz_xxxxy = pbuffer.data(idx_ovl_fh + 169);

    auto ts_yzz_xxxxz = pbuffer.data(idx_ovl_fh + 170);

    auto ts_yzz_xxxyy = pbuffer.data(idx_ovl_fh + 171);

    auto ts_yzz_xxxyz = pbuffer.data(idx_ovl_fh + 172);

    auto ts_yzz_xxxzz = pbuffer.data(idx_ovl_fh + 173);

    auto ts_yzz_xxyyy = pbuffer.data(idx_ovl_fh + 174);

    auto ts_yzz_xxyyz = pbuffer.data(idx_ovl_fh + 175);

    auto ts_yzz_xxyzz = pbuffer.data(idx_ovl_fh + 176);

    auto ts_yzz_xxzzz = pbuffer.data(idx_ovl_fh + 177);

    auto ts_yzz_xyyyy = pbuffer.data(idx_ovl_fh + 178);

    auto ts_yzz_xyyyz = pbuffer.data(idx_ovl_fh + 179);

    auto ts_yzz_xyyzz = pbuffer.data(idx_ovl_fh + 180);

    auto ts_yzz_xyzzz = pbuffer.data(idx_ovl_fh + 181);

    auto ts_yzz_xzzzz = pbuffer.data(idx_ovl_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_ovl_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_ovl_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_ovl_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_ovl_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_ovl_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_ovl_fh + 188);

    auto ts_zzz_xxxxx = pbuffer.data(idx_ovl_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_ovl_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_ovl_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_ovl_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_ovl_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_ovl_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_ovl_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_ovl_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_ovl_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_ovl_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_ovl_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_ovl_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_ovl_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_ovl_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_ovl_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_ovl_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_ovl_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_ovl_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_ovl_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_ovl_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_ovl_fh + 209);

    // Set up 0-21 components of targeted buffer : FH

    auto tk_xxx_xxxxx = pbuffer.data(idx_kin_fh);

    auto tk_xxx_xxxxy = pbuffer.data(idx_kin_fh + 1);

    auto tk_xxx_xxxxz = pbuffer.data(idx_kin_fh + 2);

    auto tk_xxx_xxxyy = pbuffer.data(idx_kin_fh + 3);

    auto tk_xxx_xxxyz = pbuffer.data(idx_kin_fh + 4);

    auto tk_xxx_xxxzz = pbuffer.data(idx_kin_fh + 5);

    auto tk_xxx_xxyyy = pbuffer.data(idx_kin_fh + 6);

    auto tk_xxx_xxyyz = pbuffer.data(idx_kin_fh + 7);

    auto tk_xxx_xxyzz = pbuffer.data(idx_kin_fh + 8);

    auto tk_xxx_xxzzz = pbuffer.data(idx_kin_fh + 9);

    auto tk_xxx_xyyyy = pbuffer.data(idx_kin_fh + 10);

    auto tk_xxx_xyyyz = pbuffer.data(idx_kin_fh + 11);

    auto tk_xxx_xyyzz = pbuffer.data(idx_kin_fh + 12);

    auto tk_xxx_xyzzz = pbuffer.data(idx_kin_fh + 13);

    auto tk_xxx_xzzzz = pbuffer.data(idx_kin_fh + 14);

    auto tk_xxx_yyyyy = pbuffer.data(idx_kin_fh + 15);

    auto tk_xxx_yyyyz = pbuffer.data(idx_kin_fh + 16);

    auto tk_xxx_yyyzz = pbuffer.data(idx_kin_fh + 17);

    auto tk_xxx_yyzzz = pbuffer.data(idx_kin_fh + 18);

    auto tk_xxx_yzzzz = pbuffer.data(idx_kin_fh + 19);

    auto tk_xxx_zzzzz = pbuffer.data(idx_kin_fh + 20);

    #pragma omp simd aligned(pa_x, tk_x_xxxxx, tk_x_xxxxy, tk_x_xxxxz, tk_x_xxxyy, tk_x_xxxyz, tk_x_xxxzz, tk_x_xxyyy, tk_x_xxyyz, tk_x_xxyzz, tk_x_xxzzz, tk_x_xyyyy, tk_x_xyyyz, tk_x_xyyzz, tk_x_xyzzz, tk_x_xzzzz, tk_x_yyyyy, tk_x_yyyyz, tk_x_yyyzz, tk_x_yyzzz, tk_x_yzzzz, tk_x_zzzzz, tk_xx_xxxx, tk_xx_xxxxx, tk_xx_xxxxy, tk_xx_xxxxz, tk_xx_xxxy, tk_xx_xxxyy, tk_xx_xxxyz, tk_xx_xxxz, tk_xx_xxxzz, tk_xx_xxyy, tk_xx_xxyyy, tk_xx_xxyyz, tk_xx_xxyz, tk_xx_xxyzz, tk_xx_xxzz, tk_xx_xxzzz, tk_xx_xyyy, tk_xx_xyyyy, tk_xx_xyyyz, tk_xx_xyyz, tk_xx_xyyzz, tk_xx_xyzz, tk_xx_xyzzz, tk_xx_xzzz, tk_xx_xzzzz, tk_xx_yyyy, tk_xx_yyyyy, tk_xx_yyyyz, tk_xx_yyyz, tk_xx_yyyzz, tk_xx_yyzz, tk_xx_yyzzz, tk_xx_yzzz, tk_xx_yzzzz, tk_xx_zzzz, tk_xx_zzzzz, tk_xxx_xxxxx, tk_xxx_xxxxy, tk_xxx_xxxxz, tk_xxx_xxxyy, tk_xxx_xxxyz, tk_xxx_xxxzz, tk_xxx_xxyyy, tk_xxx_xxyyz, tk_xxx_xxyzz, tk_xxx_xxzzz, tk_xxx_xyyyy, tk_xxx_xyyyz, tk_xxx_xyyzz, tk_xxx_xyzzz, tk_xxx_xzzzz, tk_xxx_yyyyy, tk_xxx_yyyyz, tk_xxx_yyyzz, tk_xxx_yyzzz, tk_xxx_yzzzz, tk_xxx_zzzzz, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxzz, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyzz, ts_x_xxzzz, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyzz, ts_x_xyzzz, ts_x_xzzzz, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyzz, ts_x_yyzzz, ts_x_yzzzz, ts_x_zzzzz, ts_xxx_xxxxx, ts_xxx_xxxxy, ts_xxx_xxxxz, ts_xxx_xxxyy, ts_xxx_xxxyz, ts_xxx_xxxzz, ts_xxx_xxyyy, ts_xxx_xxyyz, ts_xxx_xxyzz, ts_xxx_xxzzz, ts_xxx_xyyyy, ts_xxx_xyyyz, ts_xxx_xyyzz, ts_xxx_xyzzz, ts_xxx_xzzzz, ts_xxx_yyyyy, ts_xxx_yyyyz, ts_xxx_yyyzz, ts_xxx_yyzzz, ts_xxx_yzzzz, ts_xxx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_xxxxx[i] = -4.0 * ts_x_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxx[i] * fe_0 + 5.0 * tk_xx_xxxx[i] * fe_0 + tk_xx_xxxxx[i] * pa_x[i] + 2.0 * ts_xxx_xxxxx[i] * fz_0;

        tk_xxx_xxxxy[i] = -4.0 * ts_x_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxy[i] * fe_0 + 4.0 * tk_xx_xxxy[i] * fe_0 + tk_xx_xxxxy[i] * pa_x[i] + 2.0 * ts_xxx_xxxxy[i] * fz_0;

        tk_xxx_xxxxz[i] = -4.0 * ts_x_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxz[i] * fe_0 + 4.0 * tk_xx_xxxz[i] * fe_0 + tk_xx_xxxxz[i] * pa_x[i] + 2.0 * ts_xxx_xxxxz[i] * fz_0;

        tk_xxx_xxxyy[i] = -4.0 * ts_x_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxyy[i] * fe_0 + 3.0 * tk_xx_xxyy[i] * fe_0 + tk_xx_xxxyy[i] * pa_x[i] + 2.0 * ts_xxx_xxxyy[i] * fz_0;

        tk_xxx_xxxyz[i] = -4.0 * ts_x_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxyz[i] * fe_0 + 3.0 * tk_xx_xxyz[i] * fe_0 + tk_xx_xxxyz[i] * pa_x[i] + 2.0 * ts_xxx_xxxyz[i] * fz_0;

        tk_xxx_xxxzz[i] = -4.0 * ts_x_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxzz[i] * fe_0 + 3.0 * tk_xx_xxzz[i] * fe_0 + tk_xx_xxxzz[i] * pa_x[i] + 2.0 * ts_xxx_xxxzz[i] * fz_0;

        tk_xxx_xxyyy[i] = -4.0 * ts_x_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyyy[i] * fe_0 + 2.0 * tk_xx_xyyy[i] * fe_0 + tk_xx_xxyyy[i] * pa_x[i] + 2.0 * ts_xxx_xxyyy[i] * fz_0;

        tk_xxx_xxyyz[i] = -4.0 * ts_x_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyyz[i] * fe_0 + 2.0 * tk_xx_xyyz[i] * fe_0 + tk_xx_xxyyz[i] * pa_x[i] + 2.0 * ts_xxx_xxyyz[i] * fz_0;

        tk_xxx_xxyzz[i] = -4.0 * ts_x_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyzz[i] * fe_0 + 2.0 * tk_xx_xyzz[i] * fe_0 + tk_xx_xxyzz[i] * pa_x[i] + 2.0 * ts_xxx_xxyzz[i] * fz_0;

        tk_xxx_xxzzz[i] = -4.0 * ts_x_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxzzz[i] * fe_0 + 2.0 * tk_xx_xzzz[i] * fe_0 + tk_xx_xxzzz[i] * pa_x[i] + 2.0 * ts_xxx_xxzzz[i] * fz_0;

        tk_xxx_xyyyy[i] = -4.0 * ts_x_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyyy[i] * fe_0 + tk_xx_yyyy[i] * fe_0 + tk_xx_xyyyy[i] * pa_x[i] + 2.0 * ts_xxx_xyyyy[i] * fz_0;

        tk_xxx_xyyyz[i] = -4.0 * ts_x_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyyz[i] * fe_0 + tk_xx_yyyz[i] * fe_0 + tk_xx_xyyyz[i] * pa_x[i] + 2.0 * ts_xxx_xyyyz[i] * fz_0;

        tk_xxx_xyyzz[i] = -4.0 * ts_x_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyzz[i] * fe_0 + tk_xx_yyzz[i] * fe_0 + tk_xx_xyyzz[i] * pa_x[i] + 2.0 * ts_xxx_xyyzz[i] * fz_0;

        tk_xxx_xyzzz[i] = -4.0 * ts_x_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyzzz[i] * fe_0 + tk_xx_yzzz[i] * fe_0 + tk_xx_xyzzz[i] * pa_x[i] + 2.0 * ts_xxx_xyzzz[i] * fz_0;

        tk_xxx_xzzzz[i] = -4.0 * ts_x_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xzzzz[i] * fe_0 + tk_xx_zzzz[i] * fe_0 + tk_xx_xzzzz[i] * pa_x[i] + 2.0 * ts_xxx_xzzzz[i] * fz_0;

        tk_xxx_yyyyy[i] = -4.0 * ts_x_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyyy[i] * fe_0 + tk_xx_yyyyy[i] * pa_x[i] + 2.0 * ts_xxx_yyyyy[i] * fz_0;

        tk_xxx_yyyyz[i] = -4.0 * ts_x_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyyz[i] * fe_0 + tk_xx_yyyyz[i] * pa_x[i] + 2.0 * ts_xxx_yyyyz[i] * fz_0;

        tk_xxx_yyyzz[i] = -4.0 * ts_x_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyzz[i] * fe_0 + tk_xx_yyyzz[i] * pa_x[i] + 2.0 * ts_xxx_yyyzz[i] * fz_0;

        tk_xxx_yyzzz[i] = -4.0 * ts_x_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyzzz[i] * fe_0 + tk_xx_yyzzz[i] * pa_x[i] + 2.0 * ts_xxx_yyzzz[i] * fz_0;

        tk_xxx_yzzzz[i] = -4.0 * ts_x_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yzzzz[i] * fe_0 + tk_xx_yzzzz[i] * pa_x[i] + 2.0 * ts_xxx_yzzzz[i] * fz_0;

        tk_xxx_zzzzz[i] = -4.0 * ts_x_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_zzzzz[i] * fe_0 + tk_xx_zzzzz[i] * pa_x[i] + 2.0 * ts_xxx_zzzzz[i] * fz_0;
    }

    // Set up 21-42 components of targeted buffer : FH

    auto tk_xxy_xxxxx = pbuffer.data(idx_kin_fh + 21);

    auto tk_xxy_xxxxy = pbuffer.data(idx_kin_fh + 22);

    auto tk_xxy_xxxxz = pbuffer.data(idx_kin_fh + 23);

    auto tk_xxy_xxxyy = pbuffer.data(idx_kin_fh + 24);

    auto tk_xxy_xxxyz = pbuffer.data(idx_kin_fh + 25);

    auto tk_xxy_xxxzz = pbuffer.data(idx_kin_fh + 26);

    auto tk_xxy_xxyyy = pbuffer.data(idx_kin_fh + 27);

    auto tk_xxy_xxyyz = pbuffer.data(idx_kin_fh + 28);

    auto tk_xxy_xxyzz = pbuffer.data(idx_kin_fh + 29);

    auto tk_xxy_xxzzz = pbuffer.data(idx_kin_fh + 30);

    auto tk_xxy_xyyyy = pbuffer.data(idx_kin_fh + 31);

    auto tk_xxy_xyyyz = pbuffer.data(idx_kin_fh + 32);

    auto tk_xxy_xyyzz = pbuffer.data(idx_kin_fh + 33);

    auto tk_xxy_xyzzz = pbuffer.data(idx_kin_fh + 34);

    auto tk_xxy_xzzzz = pbuffer.data(idx_kin_fh + 35);

    auto tk_xxy_yyyyy = pbuffer.data(idx_kin_fh + 36);

    auto tk_xxy_yyyyz = pbuffer.data(idx_kin_fh + 37);

    auto tk_xxy_yyyzz = pbuffer.data(idx_kin_fh + 38);

    auto tk_xxy_yyzzz = pbuffer.data(idx_kin_fh + 39);

    auto tk_xxy_yzzzz = pbuffer.data(idx_kin_fh + 40);

    auto tk_xxy_zzzzz = pbuffer.data(idx_kin_fh + 41);

    #pragma omp simd aligned(pa_y, tk_xx_xxxx, tk_xx_xxxxx, tk_xx_xxxxy, tk_xx_xxxxz, tk_xx_xxxy, tk_xx_xxxyy, tk_xx_xxxyz, tk_xx_xxxz, tk_xx_xxxzz, tk_xx_xxyy, tk_xx_xxyyy, tk_xx_xxyyz, tk_xx_xxyz, tk_xx_xxyzz, tk_xx_xxzz, tk_xx_xxzzz, tk_xx_xyyy, tk_xx_xyyyy, tk_xx_xyyyz, tk_xx_xyyz, tk_xx_xyyzz, tk_xx_xyzz, tk_xx_xyzzz, tk_xx_xzzz, tk_xx_xzzzz, tk_xx_yyyy, tk_xx_yyyyy, tk_xx_yyyyz, tk_xx_yyyz, tk_xx_yyyzz, tk_xx_yyzz, tk_xx_yyzzz, tk_xx_yzzz, tk_xx_yzzzz, tk_xx_zzzz, tk_xx_zzzzz, tk_xxy_xxxxx, tk_xxy_xxxxy, tk_xxy_xxxxz, tk_xxy_xxxyy, tk_xxy_xxxyz, tk_xxy_xxxzz, tk_xxy_xxyyy, tk_xxy_xxyyz, tk_xxy_xxyzz, tk_xxy_xxzzz, tk_xxy_xyyyy, tk_xxy_xyyyz, tk_xxy_xyyzz, tk_xxy_xyzzz, tk_xxy_xzzzz, tk_xxy_yyyyy, tk_xxy_yyyyz, tk_xxy_yyyzz, tk_xxy_yyzzz, tk_xxy_yzzzz, tk_xxy_zzzzz, ts_xxy_xxxxx, ts_xxy_xxxxy, ts_xxy_xxxxz, ts_xxy_xxxyy, ts_xxy_xxxyz, ts_xxy_xxxzz, ts_xxy_xxyyy, ts_xxy_xxyyz, ts_xxy_xxyzz, ts_xxy_xxzzz, ts_xxy_xyyyy, ts_xxy_xyyyz, ts_xxy_xyyzz, ts_xxy_xyzzz, ts_xxy_xzzzz, ts_xxy_yyyyy, ts_xxy_yyyyz, ts_xxy_yyyzz, ts_xxy_yyzzz, ts_xxy_yzzzz, ts_xxy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxy_xxxxx[i] = tk_xx_xxxxx[i] * pa_y[i] + 2.0 * ts_xxy_xxxxx[i] * fz_0;

        tk_xxy_xxxxy[i] = tk_xx_xxxx[i] * fe_0 + tk_xx_xxxxy[i] * pa_y[i] + 2.0 * ts_xxy_xxxxy[i] * fz_0;

        tk_xxy_xxxxz[i] = tk_xx_xxxxz[i] * pa_y[i] + 2.0 * ts_xxy_xxxxz[i] * fz_0;

        tk_xxy_xxxyy[i] = 2.0 * tk_xx_xxxy[i] * fe_0 + tk_xx_xxxyy[i] * pa_y[i] + 2.0 * ts_xxy_xxxyy[i] * fz_0;

        tk_xxy_xxxyz[i] = tk_xx_xxxz[i] * fe_0 + tk_xx_xxxyz[i] * pa_y[i] + 2.0 * ts_xxy_xxxyz[i] * fz_0;

        tk_xxy_xxxzz[i] = tk_xx_xxxzz[i] * pa_y[i] + 2.0 * ts_xxy_xxxzz[i] * fz_0;

        tk_xxy_xxyyy[i] = 3.0 * tk_xx_xxyy[i] * fe_0 + tk_xx_xxyyy[i] * pa_y[i] + 2.0 * ts_xxy_xxyyy[i] * fz_0;

        tk_xxy_xxyyz[i] = 2.0 * tk_xx_xxyz[i] * fe_0 + tk_xx_xxyyz[i] * pa_y[i] + 2.0 * ts_xxy_xxyyz[i] * fz_0;

        tk_xxy_xxyzz[i] = tk_xx_xxzz[i] * fe_0 + tk_xx_xxyzz[i] * pa_y[i] + 2.0 * ts_xxy_xxyzz[i] * fz_0;

        tk_xxy_xxzzz[i] = tk_xx_xxzzz[i] * pa_y[i] + 2.0 * ts_xxy_xxzzz[i] * fz_0;

        tk_xxy_xyyyy[i] = 4.0 * tk_xx_xyyy[i] * fe_0 + tk_xx_xyyyy[i] * pa_y[i] + 2.0 * ts_xxy_xyyyy[i] * fz_0;

        tk_xxy_xyyyz[i] = 3.0 * tk_xx_xyyz[i] * fe_0 + tk_xx_xyyyz[i] * pa_y[i] + 2.0 * ts_xxy_xyyyz[i] * fz_0;

        tk_xxy_xyyzz[i] = 2.0 * tk_xx_xyzz[i] * fe_0 + tk_xx_xyyzz[i] * pa_y[i] + 2.0 * ts_xxy_xyyzz[i] * fz_0;

        tk_xxy_xyzzz[i] = tk_xx_xzzz[i] * fe_0 + tk_xx_xyzzz[i] * pa_y[i] + 2.0 * ts_xxy_xyzzz[i] * fz_0;

        tk_xxy_xzzzz[i] = tk_xx_xzzzz[i] * pa_y[i] + 2.0 * ts_xxy_xzzzz[i] * fz_0;

        tk_xxy_yyyyy[i] = 5.0 * tk_xx_yyyy[i] * fe_0 + tk_xx_yyyyy[i] * pa_y[i] + 2.0 * ts_xxy_yyyyy[i] * fz_0;

        tk_xxy_yyyyz[i] = 4.0 * tk_xx_yyyz[i] * fe_0 + tk_xx_yyyyz[i] * pa_y[i] + 2.0 * ts_xxy_yyyyz[i] * fz_0;

        tk_xxy_yyyzz[i] = 3.0 * tk_xx_yyzz[i] * fe_0 + tk_xx_yyyzz[i] * pa_y[i] + 2.0 * ts_xxy_yyyzz[i] * fz_0;

        tk_xxy_yyzzz[i] = 2.0 * tk_xx_yzzz[i] * fe_0 + tk_xx_yyzzz[i] * pa_y[i] + 2.0 * ts_xxy_yyzzz[i] * fz_0;

        tk_xxy_yzzzz[i] = tk_xx_zzzz[i] * fe_0 + tk_xx_yzzzz[i] * pa_y[i] + 2.0 * ts_xxy_yzzzz[i] * fz_0;

        tk_xxy_zzzzz[i] = tk_xx_zzzzz[i] * pa_y[i] + 2.0 * ts_xxy_zzzzz[i] * fz_0;
    }

    // Set up 42-63 components of targeted buffer : FH

    auto tk_xxz_xxxxx = pbuffer.data(idx_kin_fh + 42);

    auto tk_xxz_xxxxy = pbuffer.data(idx_kin_fh + 43);

    auto tk_xxz_xxxxz = pbuffer.data(idx_kin_fh + 44);

    auto tk_xxz_xxxyy = pbuffer.data(idx_kin_fh + 45);

    auto tk_xxz_xxxyz = pbuffer.data(idx_kin_fh + 46);

    auto tk_xxz_xxxzz = pbuffer.data(idx_kin_fh + 47);

    auto tk_xxz_xxyyy = pbuffer.data(idx_kin_fh + 48);

    auto tk_xxz_xxyyz = pbuffer.data(idx_kin_fh + 49);

    auto tk_xxz_xxyzz = pbuffer.data(idx_kin_fh + 50);

    auto tk_xxz_xxzzz = pbuffer.data(idx_kin_fh + 51);

    auto tk_xxz_xyyyy = pbuffer.data(idx_kin_fh + 52);

    auto tk_xxz_xyyyz = pbuffer.data(idx_kin_fh + 53);

    auto tk_xxz_xyyzz = pbuffer.data(idx_kin_fh + 54);

    auto tk_xxz_xyzzz = pbuffer.data(idx_kin_fh + 55);

    auto tk_xxz_xzzzz = pbuffer.data(idx_kin_fh + 56);

    auto tk_xxz_yyyyy = pbuffer.data(idx_kin_fh + 57);

    auto tk_xxz_yyyyz = pbuffer.data(idx_kin_fh + 58);

    auto tk_xxz_yyyzz = pbuffer.data(idx_kin_fh + 59);

    auto tk_xxz_yyzzz = pbuffer.data(idx_kin_fh + 60);

    auto tk_xxz_yzzzz = pbuffer.data(idx_kin_fh + 61);

    auto tk_xxz_zzzzz = pbuffer.data(idx_kin_fh + 62);

    #pragma omp simd aligned(pa_z, tk_xx_xxxx, tk_xx_xxxxx, tk_xx_xxxxy, tk_xx_xxxxz, tk_xx_xxxy, tk_xx_xxxyy, tk_xx_xxxyz, tk_xx_xxxz, tk_xx_xxxzz, tk_xx_xxyy, tk_xx_xxyyy, tk_xx_xxyyz, tk_xx_xxyz, tk_xx_xxyzz, tk_xx_xxzz, tk_xx_xxzzz, tk_xx_xyyy, tk_xx_xyyyy, tk_xx_xyyyz, tk_xx_xyyz, tk_xx_xyyzz, tk_xx_xyzz, tk_xx_xyzzz, tk_xx_xzzz, tk_xx_xzzzz, tk_xx_yyyy, tk_xx_yyyyy, tk_xx_yyyyz, tk_xx_yyyz, tk_xx_yyyzz, tk_xx_yyzz, tk_xx_yyzzz, tk_xx_yzzz, tk_xx_yzzzz, tk_xx_zzzz, tk_xx_zzzzz, tk_xxz_xxxxx, tk_xxz_xxxxy, tk_xxz_xxxxz, tk_xxz_xxxyy, tk_xxz_xxxyz, tk_xxz_xxxzz, tk_xxz_xxyyy, tk_xxz_xxyyz, tk_xxz_xxyzz, tk_xxz_xxzzz, tk_xxz_xyyyy, tk_xxz_xyyyz, tk_xxz_xyyzz, tk_xxz_xyzzz, tk_xxz_xzzzz, tk_xxz_yyyyy, tk_xxz_yyyyz, tk_xxz_yyyzz, tk_xxz_yyzzz, tk_xxz_yzzzz, tk_xxz_zzzzz, ts_xxz_xxxxx, ts_xxz_xxxxy, ts_xxz_xxxxz, ts_xxz_xxxyy, ts_xxz_xxxyz, ts_xxz_xxxzz, ts_xxz_xxyyy, ts_xxz_xxyyz, ts_xxz_xxyzz, ts_xxz_xxzzz, ts_xxz_xyyyy, ts_xxz_xyyyz, ts_xxz_xyyzz, ts_xxz_xyzzz, ts_xxz_xzzzz, ts_xxz_yyyyy, ts_xxz_yyyyz, ts_xxz_yyyzz, ts_xxz_yyzzz, ts_xxz_yzzzz, ts_xxz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxz_xxxxx[i] = tk_xx_xxxxx[i] * pa_z[i] + 2.0 * ts_xxz_xxxxx[i] * fz_0;

        tk_xxz_xxxxy[i] = tk_xx_xxxxy[i] * pa_z[i] + 2.0 * ts_xxz_xxxxy[i] * fz_0;

        tk_xxz_xxxxz[i] = tk_xx_xxxx[i] * fe_0 + tk_xx_xxxxz[i] * pa_z[i] + 2.0 * ts_xxz_xxxxz[i] * fz_0;

        tk_xxz_xxxyy[i] = tk_xx_xxxyy[i] * pa_z[i] + 2.0 * ts_xxz_xxxyy[i] * fz_0;

        tk_xxz_xxxyz[i] = tk_xx_xxxy[i] * fe_0 + tk_xx_xxxyz[i] * pa_z[i] + 2.0 * ts_xxz_xxxyz[i] * fz_0;

        tk_xxz_xxxzz[i] = 2.0 * tk_xx_xxxz[i] * fe_0 + tk_xx_xxxzz[i] * pa_z[i] + 2.0 * ts_xxz_xxxzz[i] * fz_0;

        tk_xxz_xxyyy[i] = tk_xx_xxyyy[i] * pa_z[i] + 2.0 * ts_xxz_xxyyy[i] * fz_0;

        tk_xxz_xxyyz[i] = tk_xx_xxyy[i] * fe_0 + tk_xx_xxyyz[i] * pa_z[i] + 2.0 * ts_xxz_xxyyz[i] * fz_0;

        tk_xxz_xxyzz[i] = 2.0 * tk_xx_xxyz[i] * fe_0 + tk_xx_xxyzz[i] * pa_z[i] + 2.0 * ts_xxz_xxyzz[i] * fz_0;

        tk_xxz_xxzzz[i] = 3.0 * tk_xx_xxzz[i] * fe_0 + tk_xx_xxzzz[i] * pa_z[i] + 2.0 * ts_xxz_xxzzz[i] * fz_0;

        tk_xxz_xyyyy[i] = tk_xx_xyyyy[i] * pa_z[i] + 2.0 * ts_xxz_xyyyy[i] * fz_0;

        tk_xxz_xyyyz[i] = tk_xx_xyyy[i] * fe_0 + tk_xx_xyyyz[i] * pa_z[i] + 2.0 * ts_xxz_xyyyz[i] * fz_0;

        tk_xxz_xyyzz[i] = 2.0 * tk_xx_xyyz[i] * fe_0 + tk_xx_xyyzz[i] * pa_z[i] + 2.0 * ts_xxz_xyyzz[i] * fz_0;

        tk_xxz_xyzzz[i] = 3.0 * tk_xx_xyzz[i] * fe_0 + tk_xx_xyzzz[i] * pa_z[i] + 2.0 * ts_xxz_xyzzz[i] * fz_0;

        tk_xxz_xzzzz[i] = 4.0 * tk_xx_xzzz[i] * fe_0 + tk_xx_xzzzz[i] * pa_z[i] + 2.0 * ts_xxz_xzzzz[i] * fz_0;

        tk_xxz_yyyyy[i] = tk_xx_yyyyy[i] * pa_z[i] + 2.0 * ts_xxz_yyyyy[i] * fz_0;

        tk_xxz_yyyyz[i] = tk_xx_yyyy[i] * fe_0 + tk_xx_yyyyz[i] * pa_z[i] + 2.0 * ts_xxz_yyyyz[i] * fz_0;

        tk_xxz_yyyzz[i] = 2.0 * tk_xx_yyyz[i] * fe_0 + tk_xx_yyyzz[i] * pa_z[i] + 2.0 * ts_xxz_yyyzz[i] * fz_0;

        tk_xxz_yyzzz[i] = 3.0 * tk_xx_yyzz[i] * fe_0 + tk_xx_yyzzz[i] * pa_z[i] + 2.0 * ts_xxz_yyzzz[i] * fz_0;

        tk_xxz_yzzzz[i] = 4.0 * tk_xx_yzzz[i] * fe_0 + tk_xx_yzzzz[i] * pa_z[i] + 2.0 * ts_xxz_yzzzz[i] * fz_0;

        tk_xxz_zzzzz[i] = 5.0 * tk_xx_zzzz[i] * fe_0 + tk_xx_zzzzz[i] * pa_z[i] + 2.0 * ts_xxz_zzzzz[i] * fz_0;
    }

    // Set up 63-84 components of targeted buffer : FH

    auto tk_xyy_xxxxx = pbuffer.data(idx_kin_fh + 63);

    auto tk_xyy_xxxxy = pbuffer.data(idx_kin_fh + 64);

    auto tk_xyy_xxxxz = pbuffer.data(idx_kin_fh + 65);

    auto tk_xyy_xxxyy = pbuffer.data(idx_kin_fh + 66);

    auto tk_xyy_xxxyz = pbuffer.data(idx_kin_fh + 67);

    auto tk_xyy_xxxzz = pbuffer.data(idx_kin_fh + 68);

    auto tk_xyy_xxyyy = pbuffer.data(idx_kin_fh + 69);

    auto tk_xyy_xxyyz = pbuffer.data(idx_kin_fh + 70);

    auto tk_xyy_xxyzz = pbuffer.data(idx_kin_fh + 71);

    auto tk_xyy_xxzzz = pbuffer.data(idx_kin_fh + 72);

    auto tk_xyy_xyyyy = pbuffer.data(idx_kin_fh + 73);

    auto tk_xyy_xyyyz = pbuffer.data(idx_kin_fh + 74);

    auto tk_xyy_xyyzz = pbuffer.data(idx_kin_fh + 75);

    auto tk_xyy_xyzzz = pbuffer.data(idx_kin_fh + 76);

    auto tk_xyy_xzzzz = pbuffer.data(idx_kin_fh + 77);

    auto tk_xyy_yyyyy = pbuffer.data(idx_kin_fh + 78);

    auto tk_xyy_yyyyz = pbuffer.data(idx_kin_fh + 79);

    auto tk_xyy_yyyzz = pbuffer.data(idx_kin_fh + 80);

    auto tk_xyy_yyzzz = pbuffer.data(idx_kin_fh + 81);

    auto tk_xyy_yzzzz = pbuffer.data(idx_kin_fh + 82);

    auto tk_xyy_zzzzz = pbuffer.data(idx_kin_fh + 83);

    #pragma omp simd aligned(pa_x, tk_xyy_xxxxx, tk_xyy_xxxxy, tk_xyy_xxxxz, tk_xyy_xxxyy, tk_xyy_xxxyz, tk_xyy_xxxzz, tk_xyy_xxyyy, tk_xyy_xxyyz, tk_xyy_xxyzz, tk_xyy_xxzzz, tk_xyy_xyyyy, tk_xyy_xyyyz, tk_xyy_xyyzz, tk_xyy_xyzzz, tk_xyy_xzzzz, tk_xyy_yyyyy, tk_xyy_yyyyz, tk_xyy_yyyzz, tk_xyy_yyzzz, tk_xyy_yzzzz, tk_xyy_zzzzz, tk_yy_xxxx, tk_yy_xxxxx, tk_yy_xxxxy, tk_yy_xxxxz, tk_yy_xxxy, tk_yy_xxxyy, tk_yy_xxxyz, tk_yy_xxxz, tk_yy_xxxzz, tk_yy_xxyy, tk_yy_xxyyy, tk_yy_xxyyz, tk_yy_xxyz, tk_yy_xxyzz, tk_yy_xxzz, tk_yy_xxzzz, tk_yy_xyyy, tk_yy_xyyyy, tk_yy_xyyyz, tk_yy_xyyz, tk_yy_xyyzz, tk_yy_xyzz, tk_yy_xyzzz, tk_yy_xzzz, tk_yy_xzzzz, tk_yy_yyyy, tk_yy_yyyyy, tk_yy_yyyyz, tk_yy_yyyz, tk_yy_yyyzz, tk_yy_yyzz, tk_yy_yyzzz, tk_yy_yzzz, tk_yy_yzzzz, tk_yy_zzzz, tk_yy_zzzzz, ts_xyy_xxxxx, ts_xyy_xxxxy, ts_xyy_xxxxz, ts_xyy_xxxyy, ts_xyy_xxxyz, ts_xyy_xxxzz, ts_xyy_xxyyy, ts_xyy_xxyyz, ts_xyy_xxyzz, ts_xyy_xxzzz, ts_xyy_xyyyy, ts_xyy_xyyyz, ts_xyy_xyyzz, ts_xyy_xyzzz, ts_xyy_xzzzz, ts_xyy_yyyyy, ts_xyy_yyyyz, ts_xyy_yyyzz, ts_xyy_yyzzz, ts_xyy_yzzzz, ts_xyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyy_xxxxx[i] = 5.0 * tk_yy_xxxx[i] * fe_0 + tk_yy_xxxxx[i] * pa_x[i] + 2.0 * ts_xyy_xxxxx[i] * fz_0;

        tk_xyy_xxxxy[i] = 4.0 * tk_yy_xxxy[i] * fe_0 + tk_yy_xxxxy[i] * pa_x[i] + 2.0 * ts_xyy_xxxxy[i] * fz_0;

        tk_xyy_xxxxz[i] = 4.0 * tk_yy_xxxz[i] * fe_0 + tk_yy_xxxxz[i] * pa_x[i] + 2.0 * ts_xyy_xxxxz[i] * fz_0;

        tk_xyy_xxxyy[i] = 3.0 * tk_yy_xxyy[i] * fe_0 + tk_yy_xxxyy[i] * pa_x[i] + 2.0 * ts_xyy_xxxyy[i] * fz_0;

        tk_xyy_xxxyz[i] = 3.0 * tk_yy_xxyz[i] * fe_0 + tk_yy_xxxyz[i] * pa_x[i] + 2.0 * ts_xyy_xxxyz[i] * fz_0;

        tk_xyy_xxxzz[i] = 3.0 * tk_yy_xxzz[i] * fe_0 + tk_yy_xxxzz[i] * pa_x[i] + 2.0 * ts_xyy_xxxzz[i] * fz_0;

        tk_xyy_xxyyy[i] = 2.0 * tk_yy_xyyy[i] * fe_0 + tk_yy_xxyyy[i] * pa_x[i] + 2.0 * ts_xyy_xxyyy[i] * fz_0;

        tk_xyy_xxyyz[i] = 2.0 * tk_yy_xyyz[i] * fe_0 + tk_yy_xxyyz[i] * pa_x[i] + 2.0 * ts_xyy_xxyyz[i] * fz_0;

        tk_xyy_xxyzz[i] = 2.0 * tk_yy_xyzz[i] * fe_0 + tk_yy_xxyzz[i] * pa_x[i] + 2.0 * ts_xyy_xxyzz[i] * fz_0;

        tk_xyy_xxzzz[i] = 2.0 * tk_yy_xzzz[i] * fe_0 + tk_yy_xxzzz[i] * pa_x[i] + 2.0 * ts_xyy_xxzzz[i] * fz_0;

        tk_xyy_xyyyy[i] = tk_yy_yyyy[i] * fe_0 + tk_yy_xyyyy[i] * pa_x[i] + 2.0 * ts_xyy_xyyyy[i] * fz_0;

        tk_xyy_xyyyz[i] = tk_yy_yyyz[i] * fe_0 + tk_yy_xyyyz[i] * pa_x[i] + 2.0 * ts_xyy_xyyyz[i] * fz_0;

        tk_xyy_xyyzz[i] = tk_yy_yyzz[i] * fe_0 + tk_yy_xyyzz[i] * pa_x[i] + 2.0 * ts_xyy_xyyzz[i] * fz_0;

        tk_xyy_xyzzz[i] = tk_yy_yzzz[i] * fe_0 + tk_yy_xyzzz[i] * pa_x[i] + 2.0 * ts_xyy_xyzzz[i] * fz_0;

        tk_xyy_xzzzz[i] = tk_yy_zzzz[i] * fe_0 + tk_yy_xzzzz[i] * pa_x[i] + 2.0 * ts_xyy_xzzzz[i] * fz_0;

        tk_xyy_yyyyy[i] = tk_yy_yyyyy[i] * pa_x[i] + 2.0 * ts_xyy_yyyyy[i] * fz_0;

        tk_xyy_yyyyz[i] = tk_yy_yyyyz[i] * pa_x[i] + 2.0 * ts_xyy_yyyyz[i] * fz_0;

        tk_xyy_yyyzz[i] = tk_yy_yyyzz[i] * pa_x[i] + 2.0 * ts_xyy_yyyzz[i] * fz_0;

        tk_xyy_yyzzz[i] = tk_yy_yyzzz[i] * pa_x[i] + 2.0 * ts_xyy_yyzzz[i] * fz_0;

        tk_xyy_yzzzz[i] = tk_yy_yzzzz[i] * pa_x[i] + 2.0 * ts_xyy_yzzzz[i] * fz_0;

        tk_xyy_zzzzz[i] = tk_yy_zzzzz[i] * pa_x[i] + 2.0 * ts_xyy_zzzzz[i] * fz_0;
    }

    // Set up 84-105 components of targeted buffer : FH

    auto tk_xyz_xxxxx = pbuffer.data(idx_kin_fh + 84);

    auto tk_xyz_xxxxy = pbuffer.data(idx_kin_fh + 85);

    auto tk_xyz_xxxxz = pbuffer.data(idx_kin_fh + 86);

    auto tk_xyz_xxxyy = pbuffer.data(idx_kin_fh + 87);

    auto tk_xyz_xxxyz = pbuffer.data(idx_kin_fh + 88);

    auto tk_xyz_xxxzz = pbuffer.data(idx_kin_fh + 89);

    auto tk_xyz_xxyyy = pbuffer.data(idx_kin_fh + 90);

    auto tk_xyz_xxyyz = pbuffer.data(idx_kin_fh + 91);

    auto tk_xyz_xxyzz = pbuffer.data(idx_kin_fh + 92);

    auto tk_xyz_xxzzz = pbuffer.data(idx_kin_fh + 93);

    auto tk_xyz_xyyyy = pbuffer.data(idx_kin_fh + 94);

    auto tk_xyz_xyyyz = pbuffer.data(idx_kin_fh + 95);

    auto tk_xyz_xyyzz = pbuffer.data(idx_kin_fh + 96);

    auto tk_xyz_xyzzz = pbuffer.data(idx_kin_fh + 97);

    auto tk_xyz_xzzzz = pbuffer.data(idx_kin_fh + 98);

    auto tk_xyz_yyyyy = pbuffer.data(idx_kin_fh + 99);

    auto tk_xyz_yyyyz = pbuffer.data(idx_kin_fh + 100);

    auto tk_xyz_yyyzz = pbuffer.data(idx_kin_fh + 101);

    auto tk_xyz_yyzzz = pbuffer.data(idx_kin_fh + 102);

    auto tk_xyz_yzzzz = pbuffer.data(idx_kin_fh + 103);

    auto tk_xyz_zzzzz = pbuffer.data(idx_kin_fh + 104);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tk_xy_xxxxy, tk_xy_xxxyy, tk_xy_xxyyy, tk_xy_xyyyy, tk_xyz_xxxxx, tk_xyz_xxxxy, tk_xyz_xxxxz, tk_xyz_xxxyy, tk_xyz_xxxyz, tk_xyz_xxxzz, tk_xyz_xxyyy, tk_xyz_xxyyz, tk_xyz_xxyzz, tk_xyz_xxzzz, tk_xyz_xyyyy, tk_xyz_xyyyz, tk_xyz_xyyzz, tk_xyz_xyzzz, tk_xyz_xzzzz, tk_xyz_yyyyy, tk_xyz_yyyyz, tk_xyz_yyyzz, tk_xyz_yyzzz, tk_xyz_yzzzz, tk_xyz_zzzzz, tk_xz_xxxxx, tk_xz_xxxxz, tk_xz_xxxzz, tk_xz_xxzzz, tk_xz_xzzzz, tk_yz_xxxyz, tk_yz_xxyyz, tk_yz_xxyz, tk_yz_xxyzz, tk_yz_xyyyz, tk_yz_xyyz, tk_yz_xyyzz, tk_yz_xyzz, tk_yz_xyzzz, tk_yz_yyyyy, tk_yz_yyyyz, tk_yz_yyyz, tk_yz_yyyzz, tk_yz_yyzz, tk_yz_yyzzz, tk_yz_yzzz, tk_yz_yzzzz, tk_yz_zzzzz, ts_xyz_xxxxx, ts_xyz_xxxxy, ts_xyz_xxxxz, ts_xyz_xxxyy, ts_xyz_xxxyz, ts_xyz_xxxzz, ts_xyz_xxyyy, ts_xyz_xxyyz, ts_xyz_xxyzz, ts_xyz_xxzzz, ts_xyz_xyyyy, ts_xyz_xyyyz, ts_xyz_xyyzz, ts_xyz_xyzzz, ts_xyz_xzzzz, ts_xyz_yyyyy, ts_xyz_yyyyz, ts_xyz_yyyzz, ts_xyz_yyzzz, ts_xyz_yzzzz, ts_xyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyz_xxxxx[i] = tk_xz_xxxxx[i] * pa_y[i] + 2.0 * ts_xyz_xxxxx[i] * fz_0;

        tk_xyz_xxxxy[i] = tk_xy_xxxxy[i] * pa_z[i] + 2.0 * ts_xyz_xxxxy[i] * fz_0;

        tk_xyz_xxxxz[i] = tk_xz_xxxxz[i] * pa_y[i] + 2.0 * ts_xyz_xxxxz[i] * fz_0;

        tk_xyz_xxxyy[i] = tk_xy_xxxyy[i] * pa_z[i] + 2.0 * ts_xyz_xxxyy[i] * fz_0;

        tk_xyz_xxxyz[i] = 3.0 * tk_yz_xxyz[i] * fe_0 + tk_yz_xxxyz[i] * pa_x[i] + 2.0 * ts_xyz_xxxyz[i] * fz_0;

        tk_xyz_xxxzz[i] = tk_xz_xxxzz[i] * pa_y[i] + 2.0 * ts_xyz_xxxzz[i] * fz_0;

        tk_xyz_xxyyy[i] = tk_xy_xxyyy[i] * pa_z[i] + 2.0 * ts_xyz_xxyyy[i] * fz_0;

        tk_xyz_xxyyz[i] = 2.0 * tk_yz_xyyz[i] * fe_0 + tk_yz_xxyyz[i] * pa_x[i] + 2.0 * ts_xyz_xxyyz[i] * fz_0;

        tk_xyz_xxyzz[i] = 2.0 * tk_yz_xyzz[i] * fe_0 + tk_yz_xxyzz[i] * pa_x[i] + 2.0 * ts_xyz_xxyzz[i] * fz_0;

        tk_xyz_xxzzz[i] = tk_xz_xxzzz[i] * pa_y[i] + 2.0 * ts_xyz_xxzzz[i] * fz_0;

        tk_xyz_xyyyy[i] = tk_xy_xyyyy[i] * pa_z[i] + 2.0 * ts_xyz_xyyyy[i] * fz_0;

        tk_xyz_xyyyz[i] = tk_yz_yyyz[i] * fe_0 + tk_yz_xyyyz[i] * pa_x[i] + 2.0 * ts_xyz_xyyyz[i] * fz_0;

        tk_xyz_xyyzz[i] = tk_yz_yyzz[i] * fe_0 + tk_yz_xyyzz[i] * pa_x[i] + 2.0 * ts_xyz_xyyzz[i] * fz_0;

        tk_xyz_xyzzz[i] = tk_yz_yzzz[i] * fe_0 + tk_yz_xyzzz[i] * pa_x[i] + 2.0 * ts_xyz_xyzzz[i] * fz_0;

        tk_xyz_xzzzz[i] = tk_xz_xzzzz[i] * pa_y[i] + 2.0 * ts_xyz_xzzzz[i] * fz_0;

        tk_xyz_yyyyy[i] = tk_yz_yyyyy[i] * pa_x[i] + 2.0 * ts_xyz_yyyyy[i] * fz_0;

        tk_xyz_yyyyz[i] = tk_yz_yyyyz[i] * pa_x[i] + 2.0 * ts_xyz_yyyyz[i] * fz_0;

        tk_xyz_yyyzz[i] = tk_yz_yyyzz[i] * pa_x[i] + 2.0 * ts_xyz_yyyzz[i] * fz_0;

        tk_xyz_yyzzz[i] = tk_yz_yyzzz[i] * pa_x[i] + 2.0 * ts_xyz_yyzzz[i] * fz_0;

        tk_xyz_yzzzz[i] = tk_yz_yzzzz[i] * pa_x[i] + 2.0 * ts_xyz_yzzzz[i] * fz_0;

        tk_xyz_zzzzz[i] = tk_yz_zzzzz[i] * pa_x[i] + 2.0 * ts_xyz_zzzzz[i] * fz_0;
    }

    // Set up 105-126 components of targeted buffer : FH

    auto tk_xzz_xxxxx = pbuffer.data(idx_kin_fh + 105);

    auto tk_xzz_xxxxy = pbuffer.data(idx_kin_fh + 106);

    auto tk_xzz_xxxxz = pbuffer.data(idx_kin_fh + 107);

    auto tk_xzz_xxxyy = pbuffer.data(idx_kin_fh + 108);

    auto tk_xzz_xxxyz = pbuffer.data(idx_kin_fh + 109);

    auto tk_xzz_xxxzz = pbuffer.data(idx_kin_fh + 110);

    auto tk_xzz_xxyyy = pbuffer.data(idx_kin_fh + 111);

    auto tk_xzz_xxyyz = pbuffer.data(idx_kin_fh + 112);

    auto tk_xzz_xxyzz = pbuffer.data(idx_kin_fh + 113);

    auto tk_xzz_xxzzz = pbuffer.data(idx_kin_fh + 114);

    auto tk_xzz_xyyyy = pbuffer.data(idx_kin_fh + 115);

    auto tk_xzz_xyyyz = pbuffer.data(idx_kin_fh + 116);

    auto tk_xzz_xyyzz = pbuffer.data(idx_kin_fh + 117);

    auto tk_xzz_xyzzz = pbuffer.data(idx_kin_fh + 118);

    auto tk_xzz_xzzzz = pbuffer.data(idx_kin_fh + 119);

    auto tk_xzz_yyyyy = pbuffer.data(idx_kin_fh + 120);

    auto tk_xzz_yyyyz = pbuffer.data(idx_kin_fh + 121);

    auto tk_xzz_yyyzz = pbuffer.data(idx_kin_fh + 122);

    auto tk_xzz_yyzzz = pbuffer.data(idx_kin_fh + 123);

    auto tk_xzz_yzzzz = pbuffer.data(idx_kin_fh + 124);

    auto tk_xzz_zzzzz = pbuffer.data(idx_kin_fh + 125);

    #pragma omp simd aligned(pa_x, tk_xzz_xxxxx, tk_xzz_xxxxy, tk_xzz_xxxxz, tk_xzz_xxxyy, tk_xzz_xxxyz, tk_xzz_xxxzz, tk_xzz_xxyyy, tk_xzz_xxyyz, tk_xzz_xxyzz, tk_xzz_xxzzz, tk_xzz_xyyyy, tk_xzz_xyyyz, tk_xzz_xyyzz, tk_xzz_xyzzz, tk_xzz_xzzzz, tk_xzz_yyyyy, tk_xzz_yyyyz, tk_xzz_yyyzz, tk_xzz_yyzzz, tk_xzz_yzzzz, tk_xzz_zzzzz, tk_zz_xxxx, tk_zz_xxxxx, tk_zz_xxxxy, tk_zz_xxxxz, tk_zz_xxxy, tk_zz_xxxyy, tk_zz_xxxyz, tk_zz_xxxz, tk_zz_xxxzz, tk_zz_xxyy, tk_zz_xxyyy, tk_zz_xxyyz, tk_zz_xxyz, tk_zz_xxyzz, tk_zz_xxzz, tk_zz_xxzzz, tk_zz_xyyy, tk_zz_xyyyy, tk_zz_xyyyz, tk_zz_xyyz, tk_zz_xyyzz, tk_zz_xyzz, tk_zz_xyzzz, tk_zz_xzzz, tk_zz_xzzzz, tk_zz_yyyy, tk_zz_yyyyy, tk_zz_yyyyz, tk_zz_yyyz, tk_zz_yyyzz, tk_zz_yyzz, tk_zz_yyzzz, tk_zz_yzzz, tk_zz_yzzzz, tk_zz_zzzz, tk_zz_zzzzz, ts_xzz_xxxxx, ts_xzz_xxxxy, ts_xzz_xxxxz, ts_xzz_xxxyy, ts_xzz_xxxyz, ts_xzz_xxxzz, ts_xzz_xxyyy, ts_xzz_xxyyz, ts_xzz_xxyzz, ts_xzz_xxzzz, ts_xzz_xyyyy, ts_xzz_xyyyz, ts_xzz_xyyzz, ts_xzz_xyzzz, ts_xzz_xzzzz, ts_xzz_yyyyy, ts_xzz_yyyyz, ts_xzz_yyyzz, ts_xzz_yyzzz, ts_xzz_yzzzz, ts_xzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzz_xxxxx[i] = 5.0 * tk_zz_xxxx[i] * fe_0 + tk_zz_xxxxx[i] * pa_x[i] + 2.0 * ts_xzz_xxxxx[i] * fz_0;

        tk_xzz_xxxxy[i] = 4.0 * tk_zz_xxxy[i] * fe_0 + tk_zz_xxxxy[i] * pa_x[i] + 2.0 * ts_xzz_xxxxy[i] * fz_0;

        tk_xzz_xxxxz[i] = 4.0 * tk_zz_xxxz[i] * fe_0 + tk_zz_xxxxz[i] * pa_x[i] + 2.0 * ts_xzz_xxxxz[i] * fz_0;

        tk_xzz_xxxyy[i] = 3.0 * tk_zz_xxyy[i] * fe_0 + tk_zz_xxxyy[i] * pa_x[i] + 2.0 * ts_xzz_xxxyy[i] * fz_0;

        tk_xzz_xxxyz[i] = 3.0 * tk_zz_xxyz[i] * fe_0 + tk_zz_xxxyz[i] * pa_x[i] + 2.0 * ts_xzz_xxxyz[i] * fz_0;

        tk_xzz_xxxzz[i] = 3.0 * tk_zz_xxzz[i] * fe_0 + tk_zz_xxxzz[i] * pa_x[i] + 2.0 * ts_xzz_xxxzz[i] * fz_0;

        tk_xzz_xxyyy[i] = 2.0 * tk_zz_xyyy[i] * fe_0 + tk_zz_xxyyy[i] * pa_x[i] + 2.0 * ts_xzz_xxyyy[i] * fz_0;

        tk_xzz_xxyyz[i] = 2.0 * tk_zz_xyyz[i] * fe_0 + tk_zz_xxyyz[i] * pa_x[i] + 2.0 * ts_xzz_xxyyz[i] * fz_0;

        tk_xzz_xxyzz[i] = 2.0 * tk_zz_xyzz[i] * fe_0 + tk_zz_xxyzz[i] * pa_x[i] + 2.0 * ts_xzz_xxyzz[i] * fz_0;

        tk_xzz_xxzzz[i] = 2.0 * tk_zz_xzzz[i] * fe_0 + tk_zz_xxzzz[i] * pa_x[i] + 2.0 * ts_xzz_xxzzz[i] * fz_0;

        tk_xzz_xyyyy[i] = tk_zz_yyyy[i] * fe_0 + tk_zz_xyyyy[i] * pa_x[i] + 2.0 * ts_xzz_xyyyy[i] * fz_0;

        tk_xzz_xyyyz[i] = tk_zz_yyyz[i] * fe_0 + tk_zz_xyyyz[i] * pa_x[i] + 2.0 * ts_xzz_xyyyz[i] * fz_0;

        tk_xzz_xyyzz[i] = tk_zz_yyzz[i] * fe_0 + tk_zz_xyyzz[i] * pa_x[i] + 2.0 * ts_xzz_xyyzz[i] * fz_0;

        tk_xzz_xyzzz[i] = tk_zz_yzzz[i] * fe_0 + tk_zz_xyzzz[i] * pa_x[i] + 2.0 * ts_xzz_xyzzz[i] * fz_0;

        tk_xzz_xzzzz[i] = tk_zz_zzzz[i] * fe_0 + tk_zz_xzzzz[i] * pa_x[i] + 2.0 * ts_xzz_xzzzz[i] * fz_0;

        tk_xzz_yyyyy[i] = tk_zz_yyyyy[i] * pa_x[i] + 2.0 * ts_xzz_yyyyy[i] * fz_0;

        tk_xzz_yyyyz[i] = tk_zz_yyyyz[i] * pa_x[i] + 2.0 * ts_xzz_yyyyz[i] * fz_0;

        tk_xzz_yyyzz[i] = tk_zz_yyyzz[i] * pa_x[i] + 2.0 * ts_xzz_yyyzz[i] * fz_0;

        tk_xzz_yyzzz[i] = tk_zz_yyzzz[i] * pa_x[i] + 2.0 * ts_xzz_yyzzz[i] * fz_0;

        tk_xzz_yzzzz[i] = tk_zz_yzzzz[i] * pa_x[i] + 2.0 * ts_xzz_yzzzz[i] * fz_0;

        tk_xzz_zzzzz[i] = tk_zz_zzzzz[i] * pa_x[i] + 2.0 * ts_xzz_zzzzz[i] * fz_0;
    }

    // Set up 126-147 components of targeted buffer : FH

    auto tk_yyy_xxxxx = pbuffer.data(idx_kin_fh + 126);

    auto tk_yyy_xxxxy = pbuffer.data(idx_kin_fh + 127);

    auto tk_yyy_xxxxz = pbuffer.data(idx_kin_fh + 128);

    auto tk_yyy_xxxyy = pbuffer.data(idx_kin_fh + 129);

    auto tk_yyy_xxxyz = pbuffer.data(idx_kin_fh + 130);

    auto tk_yyy_xxxzz = pbuffer.data(idx_kin_fh + 131);

    auto tk_yyy_xxyyy = pbuffer.data(idx_kin_fh + 132);

    auto tk_yyy_xxyyz = pbuffer.data(idx_kin_fh + 133);

    auto tk_yyy_xxyzz = pbuffer.data(idx_kin_fh + 134);

    auto tk_yyy_xxzzz = pbuffer.data(idx_kin_fh + 135);

    auto tk_yyy_xyyyy = pbuffer.data(idx_kin_fh + 136);

    auto tk_yyy_xyyyz = pbuffer.data(idx_kin_fh + 137);

    auto tk_yyy_xyyzz = pbuffer.data(idx_kin_fh + 138);

    auto tk_yyy_xyzzz = pbuffer.data(idx_kin_fh + 139);

    auto tk_yyy_xzzzz = pbuffer.data(idx_kin_fh + 140);

    auto tk_yyy_yyyyy = pbuffer.data(idx_kin_fh + 141);

    auto tk_yyy_yyyyz = pbuffer.data(idx_kin_fh + 142);

    auto tk_yyy_yyyzz = pbuffer.data(idx_kin_fh + 143);

    auto tk_yyy_yyzzz = pbuffer.data(idx_kin_fh + 144);

    auto tk_yyy_yzzzz = pbuffer.data(idx_kin_fh + 145);

    auto tk_yyy_zzzzz = pbuffer.data(idx_kin_fh + 146);

    #pragma omp simd aligned(pa_y, tk_y_xxxxx, tk_y_xxxxy, tk_y_xxxxz, tk_y_xxxyy, tk_y_xxxyz, tk_y_xxxzz, tk_y_xxyyy, tk_y_xxyyz, tk_y_xxyzz, tk_y_xxzzz, tk_y_xyyyy, tk_y_xyyyz, tk_y_xyyzz, tk_y_xyzzz, tk_y_xzzzz, tk_y_yyyyy, tk_y_yyyyz, tk_y_yyyzz, tk_y_yyzzz, tk_y_yzzzz, tk_y_zzzzz, tk_yy_xxxx, tk_yy_xxxxx, tk_yy_xxxxy, tk_yy_xxxxz, tk_yy_xxxy, tk_yy_xxxyy, tk_yy_xxxyz, tk_yy_xxxz, tk_yy_xxxzz, tk_yy_xxyy, tk_yy_xxyyy, tk_yy_xxyyz, tk_yy_xxyz, tk_yy_xxyzz, tk_yy_xxzz, tk_yy_xxzzz, tk_yy_xyyy, tk_yy_xyyyy, tk_yy_xyyyz, tk_yy_xyyz, tk_yy_xyyzz, tk_yy_xyzz, tk_yy_xyzzz, tk_yy_xzzz, tk_yy_xzzzz, tk_yy_yyyy, tk_yy_yyyyy, tk_yy_yyyyz, tk_yy_yyyz, tk_yy_yyyzz, tk_yy_yyzz, tk_yy_yyzzz, tk_yy_yzzz, tk_yy_yzzzz, tk_yy_zzzz, tk_yy_zzzzz, tk_yyy_xxxxx, tk_yyy_xxxxy, tk_yyy_xxxxz, tk_yyy_xxxyy, tk_yyy_xxxyz, tk_yyy_xxxzz, tk_yyy_xxyyy, tk_yyy_xxyyz, tk_yyy_xxyzz, tk_yyy_xxzzz, tk_yyy_xyyyy, tk_yyy_xyyyz, tk_yyy_xyyzz, tk_yyy_xyzzz, tk_yyy_xzzzz, tk_yyy_yyyyy, tk_yyy_yyyyz, tk_yyy_yyyzz, tk_yyy_yyzzz, tk_yyy_yzzzz, tk_yyy_zzzzz, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxzz, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyzz, ts_y_xxzzz, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyzz, ts_y_xyzzz, ts_y_xzzzz, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyzz, ts_y_yyzzz, ts_y_yzzzz, ts_y_zzzzz, ts_yyy_xxxxx, ts_yyy_xxxxy, ts_yyy_xxxxz, ts_yyy_xxxyy, ts_yyy_xxxyz, ts_yyy_xxxzz, ts_yyy_xxyyy, ts_yyy_xxyyz, ts_yyy_xxyzz, ts_yyy_xxzzz, ts_yyy_xyyyy, ts_yyy_xyyyz, ts_yyy_xyyzz, ts_yyy_xyzzz, ts_yyy_xzzzz, ts_yyy_yyyyy, ts_yyy_yyyyz, ts_yyy_yyyzz, ts_yyy_yyzzz, ts_yyy_yzzzz, ts_yyy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyy_xxxxx[i] = -4.0 * ts_y_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxx[i] * fe_0 + tk_yy_xxxxx[i] * pa_y[i] + 2.0 * ts_yyy_xxxxx[i] * fz_0;

        tk_yyy_xxxxy[i] = -4.0 * ts_y_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxy[i] * fe_0 + tk_yy_xxxx[i] * fe_0 + tk_yy_xxxxy[i] * pa_y[i] + 2.0 * ts_yyy_xxxxy[i] * fz_0;

        tk_yyy_xxxxz[i] = -4.0 * ts_y_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxz[i] * fe_0 + tk_yy_xxxxz[i] * pa_y[i] + 2.0 * ts_yyy_xxxxz[i] * fz_0;

        tk_yyy_xxxyy[i] = -4.0 * ts_y_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxyy[i] * fe_0 + 2.0 * tk_yy_xxxy[i] * fe_0 + tk_yy_xxxyy[i] * pa_y[i] + 2.0 * ts_yyy_xxxyy[i] * fz_0;

        tk_yyy_xxxyz[i] = -4.0 * ts_y_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxyz[i] * fe_0 + tk_yy_xxxz[i] * fe_0 + tk_yy_xxxyz[i] * pa_y[i] + 2.0 * ts_yyy_xxxyz[i] * fz_0;

        tk_yyy_xxxzz[i] = -4.0 * ts_y_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxzz[i] * fe_0 + tk_yy_xxxzz[i] * pa_y[i] + 2.0 * ts_yyy_xxxzz[i] * fz_0;

        tk_yyy_xxyyy[i] = -4.0 * ts_y_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyyy[i] * fe_0 + 3.0 * tk_yy_xxyy[i] * fe_0 + tk_yy_xxyyy[i] * pa_y[i] + 2.0 * ts_yyy_xxyyy[i] * fz_0;

        tk_yyy_xxyyz[i] = -4.0 * ts_y_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyyz[i] * fe_0 + 2.0 * tk_yy_xxyz[i] * fe_0 + tk_yy_xxyyz[i] * pa_y[i] + 2.0 * ts_yyy_xxyyz[i] * fz_0;

        tk_yyy_xxyzz[i] = -4.0 * ts_y_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyzz[i] * fe_0 + tk_yy_xxzz[i] * fe_0 + tk_yy_xxyzz[i] * pa_y[i] + 2.0 * ts_yyy_xxyzz[i] * fz_0;

        tk_yyy_xxzzz[i] = -4.0 * ts_y_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxzzz[i] * fe_0 + tk_yy_xxzzz[i] * pa_y[i] + 2.0 * ts_yyy_xxzzz[i] * fz_0;

        tk_yyy_xyyyy[i] = -4.0 * ts_y_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyyy[i] * fe_0 + 4.0 * tk_yy_xyyy[i] * fe_0 + tk_yy_xyyyy[i] * pa_y[i] + 2.0 * ts_yyy_xyyyy[i] * fz_0;

        tk_yyy_xyyyz[i] = -4.0 * ts_y_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyyz[i] * fe_0 + 3.0 * tk_yy_xyyz[i] * fe_0 + tk_yy_xyyyz[i] * pa_y[i] + 2.0 * ts_yyy_xyyyz[i] * fz_0;

        tk_yyy_xyyzz[i] = -4.0 * ts_y_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyzz[i] * fe_0 + 2.0 * tk_yy_xyzz[i] * fe_0 + tk_yy_xyyzz[i] * pa_y[i] + 2.0 * ts_yyy_xyyzz[i] * fz_0;

        tk_yyy_xyzzz[i] = -4.0 * ts_y_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyzzz[i] * fe_0 + tk_yy_xzzz[i] * fe_0 + tk_yy_xyzzz[i] * pa_y[i] + 2.0 * ts_yyy_xyzzz[i] * fz_0;

        tk_yyy_xzzzz[i] = -4.0 * ts_y_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xzzzz[i] * fe_0 + tk_yy_xzzzz[i] * pa_y[i] + 2.0 * ts_yyy_xzzzz[i] * fz_0;

        tk_yyy_yyyyy[i] = -4.0 * ts_y_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyyy[i] * fe_0 + 5.0 * tk_yy_yyyy[i] * fe_0 + tk_yy_yyyyy[i] * pa_y[i] + 2.0 * ts_yyy_yyyyy[i] * fz_0;

        tk_yyy_yyyyz[i] = -4.0 * ts_y_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyyz[i] * fe_0 + 4.0 * tk_yy_yyyz[i] * fe_0 + tk_yy_yyyyz[i] * pa_y[i] + 2.0 * ts_yyy_yyyyz[i] * fz_0;

        tk_yyy_yyyzz[i] = -4.0 * ts_y_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyzz[i] * fe_0 + 3.0 * tk_yy_yyzz[i] * fe_0 + tk_yy_yyyzz[i] * pa_y[i] + 2.0 * ts_yyy_yyyzz[i] * fz_0;

        tk_yyy_yyzzz[i] = -4.0 * ts_y_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyzzz[i] * fe_0 + 2.0 * tk_yy_yzzz[i] * fe_0 + tk_yy_yyzzz[i] * pa_y[i] + 2.0 * ts_yyy_yyzzz[i] * fz_0;

        tk_yyy_yzzzz[i] = -4.0 * ts_y_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yzzzz[i] * fe_0 + tk_yy_zzzz[i] * fe_0 + tk_yy_yzzzz[i] * pa_y[i] + 2.0 * ts_yyy_yzzzz[i] * fz_0;

        tk_yyy_zzzzz[i] = -4.0 * ts_y_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_zzzzz[i] * fe_0 + tk_yy_zzzzz[i] * pa_y[i] + 2.0 * ts_yyy_zzzzz[i] * fz_0;
    }

    // Set up 147-168 components of targeted buffer : FH

    auto tk_yyz_xxxxx = pbuffer.data(idx_kin_fh + 147);

    auto tk_yyz_xxxxy = pbuffer.data(idx_kin_fh + 148);

    auto tk_yyz_xxxxz = pbuffer.data(idx_kin_fh + 149);

    auto tk_yyz_xxxyy = pbuffer.data(idx_kin_fh + 150);

    auto tk_yyz_xxxyz = pbuffer.data(idx_kin_fh + 151);

    auto tk_yyz_xxxzz = pbuffer.data(idx_kin_fh + 152);

    auto tk_yyz_xxyyy = pbuffer.data(idx_kin_fh + 153);

    auto tk_yyz_xxyyz = pbuffer.data(idx_kin_fh + 154);

    auto tk_yyz_xxyzz = pbuffer.data(idx_kin_fh + 155);

    auto tk_yyz_xxzzz = pbuffer.data(idx_kin_fh + 156);

    auto tk_yyz_xyyyy = pbuffer.data(idx_kin_fh + 157);

    auto tk_yyz_xyyyz = pbuffer.data(idx_kin_fh + 158);

    auto tk_yyz_xyyzz = pbuffer.data(idx_kin_fh + 159);

    auto tk_yyz_xyzzz = pbuffer.data(idx_kin_fh + 160);

    auto tk_yyz_xzzzz = pbuffer.data(idx_kin_fh + 161);

    auto tk_yyz_yyyyy = pbuffer.data(idx_kin_fh + 162);

    auto tk_yyz_yyyyz = pbuffer.data(idx_kin_fh + 163);

    auto tk_yyz_yyyzz = pbuffer.data(idx_kin_fh + 164);

    auto tk_yyz_yyzzz = pbuffer.data(idx_kin_fh + 165);

    auto tk_yyz_yzzzz = pbuffer.data(idx_kin_fh + 166);

    auto tk_yyz_zzzzz = pbuffer.data(idx_kin_fh + 167);

    #pragma omp simd aligned(pa_z, tk_yy_xxxx, tk_yy_xxxxx, tk_yy_xxxxy, tk_yy_xxxxz, tk_yy_xxxy, tk_yy_xxxyy, tk_yy_xxxyz, tk_yy_xxxz, tk_yy_xxxzz, tk_yy_xxyy, tk_yy_xxyyy, tk_yy_xxyyz, tk_yy_xxyz, tk_yy_xxyzz, tk_yy_xxzz, tk_yy_xxzzz, tk_yy_xyyy, tk_yy_xyyyy, tk_yy_xyyyz, tk_yy_xyyz, tk_yy_xyyzz, tk_yy_xyzz, tk_yy_xyzzz, tk_yy_xzzz, tk_yy_xzzzz, tk_yy_yyyy, tk_yy_yyyyy, tk_yy_yyyyz, tk_yy_yyyz, tk_yy_yyyzz, tk_yy_yyzz, tk_yy_yyzzz, tk_yy_yzzz, tk_yy_yzzzz, tk_yy_zzzz, tk_yy_zzzzz, tk_yyz_xxxxx, tk_yyz_xxxxy, tk_yyz_xxxxz, tk_yyz_xxxyy, tk_yyz_xxxyz, tk_yyz_xxxzz, tk_yyz_xxyyy, tk_yyz_xxyyz, tk_yyz_xxyzz, tk_yyz_xxzzz, tk_yyz_xyyyy, tk_yyz_xyyyz, tk_yyz_xyyzz, tk_yyz_xyzzz, tk_yyz_xzzzz, tk_yyz_yyyyy, tk_yyz_yyyyz, tk_yyz_yyyzz, tk_yyz_yyzzz, tk_yyz_yzzzz, tk_yyz_zzzzz, ts_yyz_xxxxx, ts_yyz_xxxxy, ts_yyz_xxxxz, ts_yyz_xxxyy, ts_yyz_xxxyz, ts_yyz_xxxzz, ts_yyz_xxyyy, ts_yyz_xxyyz, ts_yyz_xxyzz, ts_yyz_xxzzz, ts_yyz_xyyyy, ts_yyz_xyyyz, ts_yyz_xyyzz, ts_yyz_xyzzz, ts_yyz_xzzzz, ts_yyz_yyyyy, ts_yyz_yyyyz, ts_yyz_yyyzz, ts_yyz_yyzzz, ts_yyz_yzzzz, ts_yyz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyz_xxxxx[i] = tk_yy_xxxxx[i] * pa_z[i] + 2.0 * ts_yyz_xxxxx[i] * fz_0;

        tk_yyz_xxxxy[i] = tk_yy_xxxxy[i] * pa_z[i] + 2.0 * ts_yyz_xxxxy[i] * fz_0;

        tk_yyz_xxxxz[i] = tk_yy_xxxx[i] * fe_0 + tk_yy_xxxxz[i] * pa_z[i] + 2.0 * ts_yyz_xxxxz[i] * fz_0;

        tk_yyz_xxxyy[i] = tk_yy_xxxyy[i] * pa_z[i] + 2.0 * ts_yyz_xxxyy[i] * fz_0;

        tk_yyz_xxxyz[i] = tk_yy_xxxy[i] * fe_0 + tk_yy_xxxyz[i] * pa_z[i] + 2.0 * ts_yyz_xxxyz[i] * fz_0;

        tk_yyz_xxxzz[i] = 2.0 * tk_yy_xxxz[i] * fe_0 + tk_yy_xxxzz[i] * pa_z[i] + 2.0 * ts_yyz_xxxzz[i] * fz_0;

        tk_yyz_xxyyy[i] = tk_yy_xxyyy[i] * pa_z[i] + 2.0 * ts_yyz_xxyyy[i] * fz_0;

        tk_yyz_xxyyz[i] = tk_yy_xxyy[i] * fe_0 + tk_yy_xxyyz[i] * pa_z[i] + 2.0 * ts_yyz_xxyyz[i] * fz_0;

        tk_yyz_xxyzz[i] = 2.0 * tk_yy_xxyz[i] * fe_0 + tk_yy_xxyzz[i] * pa_z[i] + 2.0 * ts_yyz_xxyzz[i] * fz_0;

        tk_yyz_xxzzz[i] = 3.0 * tk_yy_xxzz[i] * fe_0 + tk_yy_xxzzz[i] * pa_z[i] + 2.0 * ts_yyz_xxzzz[i] * fz_0;

        tk_yyz_xyyyy[i] = tk_yy_xyyyy[i] * pa_z[i] + 2.0 * ts_yyz_xyyyy[i] * fz_0;

        tk_yyz_xyyyz[i] = tk_yy_xyyy[i] * fe_0 + tk_yy_xyyyz[i] * pa_z[i] + 2.0 * ts_yyz_xyyyz[i] * fz_0;

        tk_yyz_xyyzz[i] = 2.0 * tk_yy_xyyz[i] * fe_0 + tk_yy_xyyzz[i] * pa_z[i] + 2.0 * ts_yyz_xyyzz[i] * fz_0;

        tk_yyz_xyzzz[i] = 3.0 * tk_yy_xyzz[i] * fe_0 + tk_yy_xyzzz[i] * pa_z[i] + 2.0 * ts_yyz_xyzzz[i] * fz_0;

        tk_yyz_xzzzz[i] = 4.0 * tk_yy_xzzz[i] * fe_0 + tk_yy_xzzzz[i] * pa_z[i] + 2.0 * ts_yyz_xzzzz[i] * fz_0;

        tk_yyz_yyyyy[i] = tk_yy_yyyyy[i] * pa_z[i] + 2.0 * ts_yyz_yyyyy[i] * fz_0;

        tk_yyz_yyyyz[i] = tk_yy_yyyy[i] * fe_0 + tk_yy_yyyyz[i] * pa_z[i] + 2.0 * ts_yyz_yyyyz[i] * fz_0;

        tk_yyz_yyyzz[i] = 2.0 * tk_yy_yyyz[i] * fe_0 + tk_yy_yyyzz[i] * pa_z[i] + 2.0 * ts_yyz_yyyzz[i] * fz_0;

        tk_yyz_yyzzz[i] = 3.0 * tk_yy_yyzz[i] * fe_0 + tk_yy_yyzzz[i] * pa_z[i] + 2.0 * ts_yyz_yyzzz[i] * fz_0;

        tk_yyz_yzzzz[i] = 4.0 * tk_yy_yzzz[i] * fe_0 + tk_yy_yzzzz[i] * pa_z[i] + 2.0 * ts_yyz_yzzzz[i] * fz_0;

        tk_yyz_zzzzz[i] = 5.0 * tk_yy_zzzz[i] * fe_0 + tk_yy_zzzzz[i] * pa_z[i] + 2.0 * ts_yyz_zzzzz[i] * fz_0;
    }

    // Set up 168-189 components of targeted buffer : FH

    auto tk_yzz_xxxxx = pbuffer.data(idx_kin_fh + 168);

    auto tk_yzz_xxxxy = pbuffer.data(idx_kin_fh + 169);

    auto tk_yzz_xxxxz = pbuffer.data(idx_kin_fh + 170);

    auto tk_yzz_xxxyy = pbuffer.data(idx_kin_fh + 171);

    auto tk_yzz_xxxyz = pbuffer.data(idx_kin_fh + 172);

    auto tk_yzz_xxxzz = pbuffer.data(idx_kin_fh + 173);

    auto tk_yzz_xxyyy = pbuffer.data(idx_kin_fh + 174);

    auto tk_yzz_xxyyz = pbuffer.data(idx_kin_fh + 175);

    auto tk_yzz_xxyzz = pbuffer.data(idx_kin_fh + 176);

    auto tk_yzz_xxzzz = pbuffer.data(idx_kin_fh + 177);

    auto tk_yzz_xyyyy = pbuffer.data(idx_kin_fh + 178);

    auto tk_yzz_xyyyz = pbuffer.data(idx_kin_fh + 179);

    auto tk_yzz_xyyzz = pbuffer.data(idx_kin_fh + 180);

    auto tk_yzz_xyzzz = pbuffer.data(idx_kin_fh + 181);

    auto tk_yzz_xzzzz = pbuffer.data(idx_kin_fh + 182);

    auto tk_yzz_yyyyy = pbuffer.data(idx_kin_fh + 183);

    auto tk_yzz_yyyyz = pbuffer.data(idx_kin_fh + 184);

    auto tk_yzz_yyyzz = pbuffer.data(idx_kin_fh + 185);

    auto tk_yzz_yyzzz = pbuffer.data(idx_kin_fh + 186);

    auto tk_yzz_yzzzz = pbuffer.data(idx_kin_fh + 187);

    auto tk_yzz_zzzzz = pbuffer.data(idx_kin_fh + 188);

    #pragma omp simd aligned(pa_y, tk_yzz_xxxxx, tk_yzz_xxxxy, tk_yzz_xxxxz, tk_yzz_xxxyy, tk_yzz_xxxyz, tk_yzz_xxxzz, tk_yzz_xxyyy, tk_yzz_xxyyz, tk_yzz_xxyzz, tk_yzz_xxzzz, tk_yzz_xyyyy, tk_yzz_xyyyz, tk_yzz_xyyzz, tk_yzz_xyzzz, tk_yzz_xzzzz, tk_yzz_yyyyy, tk_yzz_yyyyz, tk_yzz_yyyzz, tk_yzz_yyzzz, tk_yzz_yzzzz, tk_yzz_zzzzz, tk_zz_xxxx, tk_zz_xxxxx, tk_zz_xxxxy, tk_zz_xxxxz, tk_zz_xxxy, tk_zz_xxxyy, tk_zz_xxxyz, tk_zz_xxxz, tk_zz_xxxzz, tk_zz_xxyy, tk_zz_xxyyy, tk_zz_xxyyz, tk_zz_xxyz, tk_zz_xxyzz, tk_zz_xxzz, tk_zz_xxzzz, tk_zz_xyyy, tk_zz_xyyyy, tk_zz_xyyyz, tk_zz_xyyz, tk_zz_xyyzz, tk_zz_xyzz, tk_zz_xyzzz, tk_zz_xzzz, tk_zz_xzzzz, tk_zz_yyyy, tk_zz_yyyyy, tk_zz_yyyyz, tk_zz_yyyz, tk_zz_yyyzz, tk_zz_yyzz, tk_zz_yyzzz, tk_zz_yzzz, tk_zz_yzzzz, tk_zz_zzzz, tk_zz_zzzzz, ts_yzz_xxxxx, ts_yzz_xxxxy, ts_yzz_xxxxz, ts_yzz_xxxyy, ts_yzz_xxxyz, ts_yzz_xxxzz, ts_yzz_xxyyy, ts_yzz_xxyyz, ts_yzz_xxyzz, ts_yzz_xxzzz, ts_yzz_xyyyy, ts_yzz_xyyyz, ts_yzz_xyyzz, ts_yzz_xyzzz, ts_yzz_xzzzz, ts_yzz_yyyyy, ts_yzz_yyyyz, ts_yzz_yyyzz, ts_yzz_yyzzz, ts_yzz_yzzzz, ts_yzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzz_xxxxx[i] = tk_zz_xxxxx[i] * pa_y[i] + 2.0 * ts_yzz_xxxxx[i] * fz_0;

        tk_yzz_xxxxy[i] = tk_zz_xxxx[i] * fe_0 + tk_zz_xxxxy[i] * pa_y[i] + 2.0 * ts_yzz_xxxxy[i] * fz_0;

        tk_yzz_xxxxz[i] = tk_zz_xxxxz[i] * pa_y[i] + 2.0 * ts_yzz_xxxxz[i] * fz_0;

        tk_yzz_xxxyy[i] = 2.0 * tk_zz_xxxy[i] * fe_0 + tk_zz_xxxyy[i] * pa_y[i] + 2.0 * ts_yzz_xxxyy[i] * fz_0;

        tk_yzz_xxxyz[i] = tk_zz_xxxz[i] * fe_0 + tk_zz_xxxyz[i] * pa_y[i] + 2.0 * ts_yzz_xxxyz[i] * fz_0;

        tk_yzz_xxxzz[i] = tk_zz_xxxzz[i] * pa_y[i] + 2.0 * ts_yzz_xxxzz[i] * fz_0;

        tk_yzz_xxyyy[i] = 3.0 * tk_zz_xxyy[i] * fe_0 + tk_zz_xxyyy[i] * pa_y[i] + 2.0 * ts_yzz_xxyyy[i] * fz_0;

        tk_yzz_xxyyz[i] = 2.0 * tk_zz_xxyz[i] * fe_0 + tk_zz_xxyyz[i] * pa_y[i] + 2.0 * ts_yzz_xxyyz[i] * fz_0;

        tk_yzz_xxyzz[i] = tk_zz_xxzz[i] * fe_0 + tk_zz_xxyzz[i] * pa_y[i] + 2.0 * ts_yzz_xxyzz[i] * fz_0;

        tk_yzz_xxzzz[i] = tk_zz_xxzzz[i] * pa_y[i] + 2.0 * ts_yzz_xxzzz[i] * fz_0;

        tk_yzz_xyyyy[i] = 4.0 * tk_zz_xyyy[i] * fe_0 + tk_zz_xyyyy[i] * pa_y[i] + 2.0 * ts_yzz_xyyyy[i] * fz_0;

        tk_yzz_xyyyz[i] = 3.0 * tk_zz_xyyz[i] * fe_0 + tk_zz_xyyyz[i] * pa_y[i] + 2.0 * ts_yzz_xyyyz[i] * fz_0;

        tk_yzz_xyyzz[i] = 2.0 * tk_zz_xyzz[i] * fe_0 + tk_zz_xyyzz[i] * pa_y[i] + 2.0 * ts_yzz_xyyzz[i] * fz_0;

        tk_yzz_xyzzz[i] = tk_zz_xzzz[i] * fe_0 + tk_zz_xyzzz[i] * pa_y[i] + 2.0 * ts_yzz_xyzzz[i] * fz_0;

        tk_yzz_xzzzz[i] = tk_zz_xzzzz[i] * pa_y[i] + 2.0 * ts_yzz_xzzzz[i] * fz_0;

        tk_yzz_yyyyy[i] = 5.0 * tk_zz_yyyy[i] * fe_0 + tk_zz_yyyyy[i] * pa_y[i] + 2.0 * ts_yzz_yyyyy[i] * fz_0;

        tk_yzz_yyyyz[i] = 4.0 * tk_zz_yyyz[i] * fe_0 + tk_zz_yyyyz[i] * pa_y[i] + 2.0 * ts_yzz_yyyyz[i] * fz_0;

        tk_yzz_yyyzz[i] = 3.0 * tk_zz_yyzz[i] * fe_0 + tk_zz_yyyzz[i] * pa_y[i] + 2.0 * ts_yzz_yyyzz[i] * fz_0;

        tk_yzz_yyzzz[i] = 2.0 * tk_zz_yzzz[i] * fe_0 + tk_zz_yyzzz[i] * pa_y[i] + 2.0 * ts_yzz_yyzzz[i] * fz_0;

        tk_yzz_yzzzz[i] = tk_zz_zzzz[i] * fe_0 + tk_zz_yzzzz[i] * pa_y[i] + 2.0 * ts_yzz_yzzzz[i] * fz_0;

        tk_yzz_zzzzz[i] = tk_zz_zzzzz[i] * pa_y[i] + 2.0 * ts_yzz_zzzzz[i] * fz_0;
    }

    // Set up 189-210 components of targeted buffer : FH

    auto tk_zzz_xxxxx = pbuffer.data(idx_kin_fh + 189);

    auto tk_zzz_xxxxy = pbuffer.data(idx_kin_fh + 190);

    auto tk_zzz_xxxxz = pbuffer.data(idx_kin_fh + 191);

    auto tk_zzz_xxxyy = pbuffer.data(idx_kin_fh + 192);

    auto tk_zzz_xxxyz = pbuffer.data(idx_kin_fh + 193);

    auto tk_zzz_xxxzz = pbuffer.data(idx_kin_fh + 194);

    auto tk_zzz_xxyyy = pbuffer.data(idx_kin_fh + 195);

    auto tk_zzz_xxyyz = pbuffer.data(idx_kin_fh + 196);

    auto tk_zzz_xxyzz = pbuffer.data(idx_kin_fh + 197);

    auto tk_zzz_xxzzz = pbuffer.data(idx_kin_fh + 198);

    auto tk_zzz_xyyyy = pbuffer.data(idx_kin_fh + 199);

    auto tk_zzz_xyyyz = pbuffer.data(idx_kin_fh + 200);

    auto tk_zzz_xyyzz = pbuffer.data(idx_kin_fh + 201);

    auto tk_zzz_xyzzz = pbuffer.data(idx_kin_fh + 202);

    auto tk_zzz_xzzzz = pbuffer.data(idx_kin_fh + 203);

    auto tk_zzz_yyyyy = pbuffer.data(idx_kin_fh + 204);

    auto tk_zzz_yyyyz = pbuffer.data(idx_kin_fh + 205);

    auto tk_zzz_yyyzz = pbuffer.data(idx_kin_fh + 206);

    auto tk_zzz_yyzzz = pbuffer.data(idx_kin_fh + 207);

    auto tk_zzz_yzzzz = pbuffer.data(idx_kin_fh + 208);

    auto tk_zzz_zzzzz = pbuffer.data(idx_kin_fh + 209);

    #pragma omp simd aligned(pa_z, tk_z_xxxxx, tk_z_xxxxy, tk_z_xxxxz, tk_z_xxxyy, tk_z_xxxyz, tk_z_xxxzz, tk_z_xxyyy, tk_z_xxyyz, tk_z_xxyzz, tk_z_xxzzz, tk_z_xyyyy, tk_z_xyyyz, tk_z_xyyzz, tk_z_xyzzz, tk_z_xzzzz, tk_z_yyyyy, tk_z_yyyyz, tk_z_yyyzz, tk_z_yyzzz, tk_z_yzzzz, tk_z_zzzzz, tk_zz_xxxx, tk_zz_xxxxx, tk_zz_xxxxy, tk_zz_xxxxz, tk_zz_xxxy, tk_zz_xxxyy, tk_zz_xxxyz, tk_zz_xxxz, tk_zz_xxxzz, tk_zz_xxyy, tk_zz_xxyyy, tk_zz_xxyyz, tk_zz_xxyz, tk_zz_xxyzz, tk_zz_xxzz, tk_zz_xxzzz, tk_zz_xyyy, tk_zz_xyyyy, tk_zz_xyyyz, tk_zz_xyyz, tk_zz_xyyzz, tk_zz_xyzz, tk_zz_xyzzz, tk_zz_xzzz, tk_zz_xzzzz, tk_zz_yyyy, tk_zz_yyyyy, tk_zz_yyyyz, tk_zz_yyyz, tk_zz_yyyzz, tk_zz_yyzz, tk_zz_yyzzz, tk_zz_yzzz, tk_zz_yzzzz, tk_zz_zzzz, tk_zz_zzzzz, tk_zzz_xxxxx, tk_zzz_xxxxy, tk_zzz_xxxxz, tk_zzz_xxxyy, tk_zzz_xxxyz, tk_zzz_xxxzz, tk_zzz_xxyyy, tk_zzz_xxyyz, tk_zzz_xxyzz, tk_zzz_xxzzz, tk_zzz_xyyyy, tk_zzz_xyyyz, tk_zzz_xyyzz, tk_zzz_xyzzz, tk_zzz_xzzzz, tk_zzz_yyyyy, tk_zzz_yyyyz, tk_zzz_yyyzz, tk_zzz_yyzzz, tk_zzz_yzzzz, tk_zzz_zzzzz, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxzz, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyzz, ts_z_xxzzz, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyzz, ts_z_xyzzz, ts_z_xzzzz, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyzz, ts_z_yyzzz, ts_z_yzzzz, ts_z_zzzzz, ts_zzz_xxxxx, ts_zzz_xxxxy, ts_zzz_xxxxz, ts_zzz_xxxyy, ts_zzz_xxxyz, ts_zzz_xxxzz, ts_zzz_xxyyy, ts_zzz_xxyyz, ts_zzz_xxyzz, ts_zzz_xxzzz, ts_zzz_xyyyy, ts_zzz_xyyyz, ts_zzz_xyyzz, ts_zzz_xyzzz, ts_zzz_xzzzz, ts_zzz_yyyyy, ts_zzz_yyyyz, ts_zzz_yyyzz, ts_zzz_yyzzz, ts_zzz_yzzzz, ts_zzz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzz_xxxxx[i] = -4.0 * ts_z_xxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxx[i] * fe_0 + tk_zz_xxxxx[i] * pa_z[i] + 2.0 * ts_zzz_xxxxx[i] * fz_0;

        tk_zzz_xxxxy[i] = -4.0 * ts_z_xxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxy[i] * fe_0 + tk_zz_xxxxy[i] * pa_z[i] + 2.0 * ts_zzz_xxxxy[i] * fz_0;

        tk_zzz_xxxxz[i] = -4.0 * ts_z_xxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxz[i] * fe_0 + tk_zz_xxxx[i] * fe_0 + tk_zz_xxxxz[i] * pa_z[i] + 2.0 * ts_zzz_xxxxz[i] * fz_0;

        tk_zzz_xxxyy[i] = -4.0 * ts_z_xxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxyy[i] * fe_0 + tk_zz_xxxyy[i] * pa_z[i] + 2.0 * ts_zzz_xxxyy[i] * fz_0;

        tk_zzz_xxxyz[i] = -4.0 * ts_z_xxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxyz[i] * fe_0 + tk_zz_xxxy[i] * fe_0 + tk_zz_xxxyz[i] * pa_z[i] + 2.0 * ts_zzz_xxxyz[i] * fz_0;

        tk_zzz_xxxzz[i] = -4.0 * ts_z_xxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxzz[i] * fe_0 + 2.0 * tk_zz_xxxz[i] * fe_0 + tk_zz_xxxzz[i] * pa_z[i] + 2.0 * ts_zzz_xxxzz[i] * fz_0;

        tk_zzz_xxyyy[i] = -4.0 * ts_z_xxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyyy[i] * fe_0 + tk_zz_xxyyy[i] * pa_z[i] + 2.0 * ts_zzz_xxyyy[i] * fz_0;

        tk_zzz_xxyyz[i] = -4.0 * ts_z_xxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyyz[i] * fe_0 + tk_zz_xxyy[i] * fe_0 + tk_zz_xxyyz[i] * pa_z[i] + 2.0 * ts_zzz_xxyyz[i] * fz_0;

        tk_zzz_xxyzz[i] = -4.0 * ts_z_xxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyzz[i] * fe_0 + 2.0 * tk_zz_xxyz[i] * fe_0 + tk_zz_xxyzz[i] * pa_z[i] + 2.0 * ts_zzz_xxyzz[i] * fz_0;

        tk_zzz_xxzzz[i] = -4.0 * ts_z_xxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxzzz[i] * fe_0 + 3.0 * tk_zz_xxzz[i] * fe_0 + tk_zz_xxzzz[i] * pa_z[i] + 2.0 * ts_zzz_xxzzz[i] * fz_0;

        tk_zzz_xyyyy[i] = -4.0 * ts_z_xyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyyy[i] * fe_0 + tk_zz_xyyyy[i] * pa_z[i] + 2.0 * ts_zzz_xyyyy[i] * fz_0;

        tk_zzz_xyyyz[i] = -4.0 * ts_z_xyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyyz[i] * fe_0 + tk_zz_xyyy[i] * fe_0 + tk_zz_xyyyz[i] * pa_z[i] + 2.0 * ts_zzz_xyyyz[i] * fz_0;

        tk_zzz_xyyzz[i] = -4.0 * ts_z_xyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyzz[i] * fe_0 + 2.0 * tk_zz_xyyz[i] * fe_0 + tk_zz_xyyzz[i] * pa_z[i] + 2.0 * ts_zzz_xyyzz[i] * fz_0;

        tk_zzz_xyzzz[i] = -4.0 * ts_z_xyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyzzz[i] * fe_0 + 3.0 * tk_zz_xyzz[i] * fe_0 + tk_zz_xyzzz[i] * pa_z[i] + 2.0 * ts_zzz_xyzzz[i] * fz_0;

        tk_zzz_xzzzz[i] = -4.0 * ts_z_xzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xzzzz[i] * fe_0 + 4.0 * tk_zz_xzzz[i] * fe_0 + tk_zz_xzzzz[i] * pa_z[i] + 2.0 * ts_zzz_xzzzz[i] * fz_0;

        tk_zzz_yyyyy[i] = -4.0 * ts_z_yyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyyy[i] * fe_0 + tk_zz_yyyyy[i] * pa_z[i] + 2.0 * ts_zzz_yyyyy[i] * fz_0;

        tk_zzz_yyyyz[i] = -4.0 * ts_z_yyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyyz[i] * fe_0 + tk_zz_yyyy[i] * fe_0 + tk_zz_yyyyz[i] * pa_z[i] + 2.0 * ts_zzz_yyyyz[i] * fz_0;

        tk_zzz_yyyzz[i] = -4.0 * ts_z_yyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyzz[i] * fe_0 + 2.0 * tk_zz_yyyz[i] * fe_0 + tk_zz_yyyzz[i] * pa_z[i] + 2.0 * ts_zzz_yyyzz[i] * fz_0;

        tk_zzz_yyzzz[i] = -4.0 * ts_z_yyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyzzz[i] * fe_0 + 3.0 * tk_zz_yyzz[i] * fe_0 + tk_zz_yyzzz[i] * pa_z[i] + 2.0 * ts_zzz_yyzzz[i] * fz_0;

        tk_zzz_yzzzz[i] = -4.0 * ts_z_yzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yzzzz[i] * fe_0 + 4.0 * tk_zz_yzzz[i] * fe_0 + tk_zz_yzzzz[i] * pa_z[i] + 2.0 * ts_zzz_yzzzz[i] * fz_0;

        tk_zzz_zzzzz[i] = -4.0 * ts_z_zzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_zzzzz[i] * fe_0 + 5.0 * tk_zz_zzzz[i] * fe_0 + tk_zz_zzzzz[i] * pa_z[i] + 2.0 * ts_zzz_zzzzz[i] * fz_0;
    }

}

} // kinrec namespace

