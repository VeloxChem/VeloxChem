#include "KineticEnergyPrimRecDH.hpp"

namespace kinrec { // kinrec namespace

auto
comp_prim_kinetic_energy_dh(CSimdArray<double>& pbuffer, 
                            const size_t idx_kin_dh,
                            const size_t idx_ovl_sh,
                            const size_t idx_kin_sh,
                            const size_t idx_kin_pg,
                            const size_t idx_kin_ph,
                            const size_t idx_ovl_dh,
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

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_ovl_sh);

    auto ts_0_xxxxy = pbuffer.data(idx_ovl_sh + 1);

    auto ts_0_xxxxz = pbuffer.data(idx_ovl_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_ovl_sh + 3);

    auto ts_0_xxxyz = pbuffer.data(idx_ovl_sh + 4);

    auto ts_0_xxxzz = pbuffer.data(idx_ovl_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_ovl_sh + 6);

    auto ts_0_xxyyz = pbuffer.data(idx_ovl_sh + 7);

    auto ts_0_xxyzz = pbuffer.data(idx_ovl_sh + 8);

    auto ts_0_xxzzz = pbuffer.data(idx_ovl_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_ovl_sh + 10);

    auto ts_0_xyyyz = pbuffer.data(idx_ovl_sh + 11);

    auto ts_0_xyyzz = pbuffer.data(idx_ovl_sh + 12);

    auto ts_0_xyzzz = pbuffer.data(idx_ovl_sh + 13);

    auto ts_0_xzzzz = pbuffer.data(idx_ovl_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_ovl_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_ovl_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_ovl_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_ovl_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_ovl_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_ovl_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto tk_0_xxxxx = pbuffer.data(idx_kin_sh);

    auto tk_0_xxxxy = pbuffer.data(idx_kin_sh + 1);

    auto tk_0_xxxxz = pbuffer.data(idx_kin_sh + 2);

    auto tk_0_xxxyy = pbuffer.data(idx_kin_sh + 3);

    auto tk_0_xxxyz = pbuffer.data(idx_kin_sh + 4);

    auto tk_0_xxxzz = pbuffer.data(idx_kin_sh + 5);

    auto tk_0_xxyyy = pbuffer.data(idx_kin_sh + 6);

    auto tk_0_xxyyz = pbuffer.data(idx_kin_sh + 7);

    auto tk_0_xxyzz = pbuffer.data(idx_kin_sh + 8);

    auto tk_0_xxzzz = pbuffer.data(idx_kin_sh + 9);

    auto tk_0_xyyyy = pbuffer.data(idx_kin_sh + 10);

    auto tk_0_xyyyz = pbuffer.data(idx_kin_sh + 11);

    auto tk_0_xyyzz = pbuffer.data(idx_kin_sh + 12);

    auto tk_0_xyzzz = pbuffer.data(idx_kin_sh + 13);

    auto tk_0_xzzzz = pbuffer.data(idx_kin_sh + 14);

    auto tk_0_yyyyy = pbuffer.data(idx_kin_sh + 15);

    auto tk_0_yyyyz = pbuffer.data(idx_kin_sh + 16);

    auto tk_0_yyyzz = pbuffer.data(idx_kin_sh + 17);

    auto tk_0_yyzzz = pbuffer.data(idx_kin_sh + 18);

    auto tk_0_yzzzz = pbuffer.data(idx_kin_sh + 19);

    auto tk_0_zzzzz = pbuffer.data(idx_kin_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto tk_x_xxxx = pbuffer.data(idx_kin_pg);

    auto tk_x_xxxy = pbuffer.data(idx_kin_pg + 1);

    auto tk_x_xxxz = pbuffer.data(idx_kin_pg + 2);

    auto tk_x_xxyy = pbuffer.data(idx_kin_pg + 3);

    auto tk_x_xxyz = pbuffer.data(idx_kin_pg + 4);

    auto tk_x_xxzz = pbuffer.data(idx_kin_pg + 5);

    auto tk_x_xyyy = pbuffer.data(idx_kin_pg + 6);

    auto tk_x_xyyz = pbuffer.data(idx_kin_pg + 7);

    auto tk_x_xyzz = pbuffer.data(idx_kin_pg + 8);

    auto tk_x_xzzz = pbuffer.data(idx_kin_pg + 9);

    auto tk_x_yyyy = pbuffer.data(idx_kin_pg + 10);

    auto tk_x_yyyz = pbuffer.data(idx_kin_pg + 11);

    auto tk_x_yyzz = pbuffer.data(idx_kin_pg + 12);

    auto tk_x_yzzz = pbuffer.data(idx_kin_pg + 13);

    auto tk_x_zzzz = pbuffer.data(idx_kin_pg + 14);

    auto tk_y_xxxx = pbuffer.data(idx_kin_pg + 15);

    auto tk_y_xxxy = pbuffer.data(idx_kin_pg + 16);

    auto tk_y_xxxz = pbuffer.data(idx_kin_pg + 17);

    auto tk_y_xxyy = pbuffer.data(idx_kin_pg + 18);

    auto tk_y_xxyz = pbuffer.data(idx_kin_pg + 19);

    auto tk_y_xxzz = pbuffer.data(idx_kin_pg + 20);

    auto tk_y_xyyy = pbuffer.data(idx_kin_pg + 21);

    auto tk_y_xyyz = pbuffer.data(idx_kin_pg + 22);

    auto tk_y_xyzz = pbuffer.data(idx_kin_pg + 23);

    auto tk_y_xzzz = pbuffer.data(idx_kin_pg + 24);

    auto tk_y_yyyy = pbuffer.data(idx_kin_pg + 25);

    auto tk_y_yyyz = pbuffer.data(idx_kin_pg + 26);

    auto tk_y_yyzz = pbuffer.data(idx_kin_pg + 27);

    auto tk_y_yzzz = pbuffer.data(idx_kin_pg + 28);

    auto tk_y_zzzz = pbuffer.data(idx_kin_pg + 29);

    auto tk_z_xxxx = pbuffer.data(idx_kin_pg + 30);

    auto tk_z_xxxy = pbuffer.data(idx_kin_pg + 31);

    auto tk_z_xxxz = pbuffer.data(idx_kin_pg + 32);

    auto tk_z_xxyy = pbuffer.data(idx_kin_pg + 33);

    auto tk_z_xxyz = pbuffer.data(idx_kin_pg + 34);

    auto tk_z_xxzz = pbuffer.data(idx_kin_pg + 35);

    auto tk_z_xyyy = pbuffer.data(idx_kin_pg + 36);

    auto tk_z_xyyz = pbuffer.data(idx_kin_pg + 37);

    auto tk_z_xyzz = pbuffer.data(idx_kin_pg + 38);

    auto tk_z_xzzz = pbuffer.data(idx_kin_pg + 39);

    auto tk_z_yyyy = pbuffer.data(idx_kin_pg + 40);

    auto tk_z_yyyz = pbuffer.data(idx_kin_pg + 41);

    auto tk_z_yyzz = pbuffer.data(idx_kin_pg + 42);

    auto tk_z_yzzz = pbuffer.data(idx_kin_pg + 43);

    auto tk_z_zzzz = pbuffer.data(idx_kin_pg + 44);

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

    // Set up components of auxiliary buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_ovl_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_ovl_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_ovl_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_ovl_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_ovl_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_ovl_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_ovl_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_ovl_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_ovl_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_ovl_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_ovl_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_ovl_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_ovl_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_ovl_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_ovl_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_ovl_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_ovl_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_ovl_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_ovl_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_ovl_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_ovl_dh + 20);

    auto ts_xy_xxxxx = pbuffer.data(idx_ovl_dh + 21);

    auto ts_xy_xxxxy = pbuffer.data(idx_ovl_dh + 22);

    auto ts_xy_xxxxz = pbuffer.data(idx_ovl_dh + 23);

    auto ts_xy_xxxyy = pbuffer.data(idx_ovl_dh + 24);

    auto ts_xy_xxxyz = pbuffer.data(idx_ovl_dh + 25);

    auto ts_xy_xxxzz = pbuffer.data(idx_ovl_dh + 26);

    auto ts_xy_xxyyy = pbuffer.data(idx_ovl_dh + 27);

    auto ts_xy_xxyyz = pbuffer.data(idx_ovl_dh + 28);

    auto ts_xy_xxyzz = pbuffer.data(idx_ovl_dh + 29);

    auto ts_xy_xxzzz = pbuffer.data(idx_ovl_dh + 30);

    auto ts_xy_xyyyy = pbuffer.data(idx_ovl_dh + 31);

    auto ts_xy_xyyyz = pbuffer.data(idx_ovl_dh + 32);

    auto ts_xy_xyyzz = pbuffer.data(idx_ovl_dh + 33);

    auto ts_xy_xyzzz = pbuffer.data(idx_ovl_dh + 34);

    auto ts_xy_xzzzz = pbuffer.data(idx_ovl_dh + 35);

    auto ts_xy_yyyyy = pbuffer.data(idx_ovl_dh + 36);

    auto ts_xy_yyyyz = pbuffer.data(idx_ovl_dh + 37);

    auto ts_xy_yyyzz = pbuffer.data(idx_ovl_dh + 38);

    auto ts_xy_yyzzz = pbuffer.data(idx_ovl_dh + 39);

    auto ts_xy_yzzzz = pbuffer.data(idx_ovl_dh + 40);

    auto ts_xy_zzzzz = pbuffer.data(idx_ovl_dh + 41);

    auto ts_xz_xxxxx = pbuffer.data(idx_ovl_dh + 42);

    auto ts_xz_xxxxy = pbuffer.data(idx_ovl_dh + 43);

    auto ts_xz_xxxxz = pbuffer.data(idx_ovl_dh + 44);

    auto ts_xz_xxxyy = pbuffer.data(idx_ovl_dh + 45);

    auto ts_xz_xxxyz = pbuffer.data(idx_ovl_dh + 46);

    auto ts_xz_xxxzz = pbuffer.data(idx_ovl_dh + 47);

    auto ts_xz_xxyyy = pbuffer.data(idx_ovl_dh + 48);

    auto ts_xz_xxyyz = pbuffer.data(idx_ovl_dh + 49);

    auto ts_xz_xxyzz = pbuffer.data(idx_ovl_dh + 50);

    auto ts_xz_xxzzz = pbuffer.data(idx_ovl_dh + 51);

    auto ts_xz_xyyyy = pbuffer.data(idx_ovl_dh + 52);

    auto ts_xz_xyyyz = pbuffer.data(idx_ovl_dh + 53);

    auto ts_xz_xyyzz = pbuffer.data(idx_ovl_dh + 54);

    auto ts_xz_xyzzz = pbuffer.data(idx_ovl_dh + 55);

    auto ts_xz_xzzzz = pbuffer.data(idx_ovl_dh + 56);

    auto ts_xz_yyyyy = pbuffer.data(idx_ovl_dh + 57);

    auto ts_xz_yyyyz = pbuffer.data(idx_ovl_dh + 58);

    auto ts_xz_yyyzz = pbuffer.data(idx_ovl_dh + 59);

    auto ts_xz_yyzzz = pbuffer.data(idx_ovl_dh + 60);

    auto ts_xz_yzzzz = pbuffer.data(idx_ovl_dh + 61);

    auto ts_xz_zzzzz = pbuffer.data(idx_ovl_dh + 62);

    auto ts_yy_xxxxx = pbuffer.data(idx_ovl_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_ovl_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_ovl_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_ovl_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_ovl_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_ovl_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_ovl_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_ovl_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_ovl_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_ovl_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_ovl_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_ovl_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_ovl_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_ovl_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_ovl_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_ovl_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_ovl_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_ovl_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_ovl_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_ovl_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_ovl_dh + 83);

    auto ts_yz_xxxxx = pbuffer.data(idx_ovl_dh + 84);

    auto ts_yz_xxxxy = pbuffer.data(idx_ovl_dh + 85);

    auto ts_yz_xxxxz = pbuffer.data(idx_ovl_dh + 86);

    auto ts_yz_xxxyy = pbuffer.data(idx_ovl_dh + 87);

    auto ts_yz_xxxyz = pbuffer.data(idx_ovl_dh + 88);

    auto ts_yz_xxxzz = pbuffer.data(idx_ovl_dh + 89);

    auto ts_yz_xxyyy = pbuffer.data(idx_ovl_dh + 90);

    auto ts_yz_xxyyz = pbuffer.data(idx_ovl_dh + 91);

    auto ts_yz_xxyzz = pbuffer.data(idx_ovl_dh + 92);

    auto ts_yz_xxzzz = pbuffer.data(idx_ovl_dh + 93);

    auto ts_yz_xyyyy = pbuffer.data(idx_ovl_dh + 94);

    auto ts_yz_xyyyz = pbuffer.data(idx_ovl_dh + 95);

    auto ts_yz_xyyzz = pbuffer.data(idx_ovl_dh + 96);

    auto ts_yz_xyzzz = pbuffer.data(idx_ovl_dh + 97);

    auto ts_yz_xzzzz = pbuffer.data(idx_ovl_dh + 98);

    auto ts_yz_yyyyy = pbuffer.data(idx_ovl_dh + 99);

    auto ts_yz_yyyyz = pbuffer.data(idx_ovl_dh + 100);

    auto ts_yz_yyyzz = pbuffer.data(idx_ovl_dh + 101);

    auto ts_yz_yyzzz = pbuffer.data(idx_ovl_dh + 102);

    auto ts_yz_yzzzz = pbuffer.data(idx_ovl_dh + 103);

    auto ts_yz_zzzzz = pbuffer.data(idx_ovl_dh + 104);

    auto ts_zz_xxxxx = pbuffer.data(idx_ovl_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_ovl_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_ovl_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_ovl_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_ovl_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_ovl_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_ovl_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_ovl_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_ovl_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_ovl_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_ovl_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_ovl_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_ovl_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_ovl_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_ovl_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_ovl_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_ovl_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_ovl_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_ovl_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_ovl_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_ovl_dh + 125);

    // Set up 0-21 components of targeted buffer : DH

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

    #pragma omp simd aligned(pa_x, tk_0_xxxxx, tk_0_xxxxy, tk_0_xxxxz, tk_0_xxxyy, tk_0_xxxyz, tk_0_xxxzz, tk_0_xxyyy, tk_0_xxyyz, tk_0_xxyzz, tk_0_xxzzz, tk_0_xyyyy, tk_0_xyyyz, tk_0_xyyzz, tk_0_xyzzz, tk_0_xzzzz, tk_0_yyyyy, tk_0_yyyyz, tk_0_yyyzz, tk_0_yyzzz, tk_0_yzzzz, tk_0_zzzzz, tk_x_xxxx, tk_x_xxxxx, tk_x_xxxxy, tk_x_xxxxz, tk_x_xxxy, tk_x_xxxyy, tk_x_xxxyz, tk_x_xxxz, tk_x_xxxzz, tk_x_xxyy, tk_x_xxyyy, tk_x_xxyyz, tk_x_xxyz, tk_x_xxyzz, tk_x_xxzz, tk_x_xxzzz, tk_x_xyyy, tk_x_xyyyy, tk_x_xyyyz, tk_x_xyyz, tk_x_xyyzz, tk_x_xyzz, tk_x_xyzzz, tk_x_xzzz, tk_x_xzzzz, tk_x_yyyy, tk_x_yyyyy, tk_x_yyyyz, tk_x_yyyz, tk_x_yyyzz, tk_x_yyzz, tk_x_yyzzz, tk_x_yzzz, tk_x_yzzzz, tk_x_zzzz, tk_x_zzzzz, tk_xx_xxxxx, tk_xx_xxxxy, tk_xx_xxxxz, tk_xx_xxxyy, tk_xx_xxxyz, tk_xx_xxxzz, tk_xx_xxyyy, tk_xx_xxyyz, tk_xx_xxyzz, tk_xx_xxzzz, tk_xx_xyyyy, tk_xx_xyyyz, tk_xx_xyyzz, tk_xx_xyzzz, tk_xx_xzzzz, tk_xx_yyyyy, tk_xx_yyyyz, tk_xx_yyyzz, tk_xx_yyzzz, tk_xx_yzzzz, tk_xx_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_xx_xxxxx, ts_xx_xxxxy, ts_xx_xxxxz, ts_xx_xxxyy, ts_xx_xxxyz, ts_xx_xxxzz, ts_xx_xxyyy, ts_xx_xxyyz, ts_xx_xxyzz, ts_xx_xxzzz, ts_xx_xyyyy, ts_xx_xyyyz, ts_xx_xyyzz, ts_xx_xyzzz, ts_xx_xzzzz, ts_xx_yyyyy, ts_xx_yyyyz, ts_xx_yyyzz, ts_xx_yyzzz, ts_xx_yzzzz, ts_xx_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_xxxxx[i] = -2.0 * ts_0_xxxxx[i] * fbe_0 * fz_0 + tk_0_xxxxx[i] * fe_0 + 5.0 * tk_x_xxxx[i] * fe_0 + tk_x_xxxxx[i] * pa_x[i] + 2.0 * ts_xx_xxxxx[i] * fz_0;

        tk_xx_xxxxy[i] = -2.0 * ts_0_xxxxy[i] * fbe_0 * fz_0 + tk_0_xxxxy[i] * fe_0 + 4.0 * tk_x_xxxy[i] * fe_0 + tk_x_xxxxy[i] * pa_x[i] + 2.0 * ts_xx_xxxxy[i] * fz_0;

        tk_xx_xxxxz[i] = -2.0 * ts_0_xxxxz[i] * fbe_0 * fz_0 + tk_0_xxxxz[i] * fe_0 + 4.0 * tk_x_xxxz[i] * fe_0 + tk_x_xxxxz[i] * pa_x[i] + 2.0 * ts_xx_xxxxz[i] * fz_0;

        tk_xx_xxxyy[i] = -2.0 * ts_0_xxxyy[i] * fbe_0 * fz_0 + tk_0_xxxyy[i] * fe_0 + 3.0 * tk_x_xxyy[i] * fe_0 + tk_x_xxxyy[i] * pa_x[i] + 2.0 * ts_xx_xxxyy[i] * fz_0;

        tk_xx_xxxyz[i] = -2.0 * ts_0_xxxyz[i] * fbe_0 * fz_0 + tk_0_xxxyz[i] * fe_0 + 3.0 * tk_x_xxyz[i] * fe_0 + tk_x_xxxyz[i] * pa_x[i] + 2.0 * ts_xx_xxxyz[i] * fz_0;

        tk_xx_xxxzz[i] = -2.0 * ts_0_xxxzz[i] * fbe_0 * fz_0 + tk_0_xxxzz[i] * fe_0 + 3.0 * tk_x_xxzz[i] * fe_0 + tk_x_xxxzz[i] * pa_x[i] + 2.0 * ts_xx_xxxzz[i] * fz_0;

        tk_xx_xxyyy[i] = -2.0 * ts_0_xxyyy[i] * fbe_0 * fz_0 + tk_0_xxyyy[i] * fe_0 + 2.0 * tk_x_xyyy[i] * fe_0 + tk_x_xxyyy[i] * pa_x[i] + 2.0 * ts_xx_xxyyy[i] * fz_0;

        tk_xx_xxyyz[i] = -2.0 * ts_0_xxyyz[i] * fbe_0 * fz_0 + tk_0_xxyyz[i] * fe_0 + 2.0 * tk_x_xyyz[i] * fe_0 + tk_x_xxyyz[i] * pa_x[i] + 2.0 * ts_xx_xxyyz[i] * fz_0;

        tk_xx_xxyzz[i] = -2.0 * ts_0_xxyzz[i] * fbe_0 * fz_0 + tk_0_xxyzz[i] * fe_0 + 2.0 * tk_x_xyzz[i] * fe_0 + tk_x_xxyzz[i] * pa_x[i] + 2.0 * ts_xx_xxyzz[i] * fz_0;

        tk_xx_xxzzz[i] = -2.0 * ts_0_xxzzz[i] * fbe_0 * fz_0 + tk_0_xxzzz[i] * fe_0 + 2.0 * tk_x_xzzz[i] * fe_0 + tk_x_xxzzz[i] * pa_x[i] + 2.0 * ts_xx_xxzzz[i] * fz_0;

        tk_xx_xyyyy[i] = -2.0 * ts_0_xyyyy[i] * fbe_0 * fz_0 + tk_0_xyyyy[i] * fe_0 + tk_x_yyyy[i] * fe_0 + tk_x_xyyyy[i] * pa_x[i] + 2.0 * ts_xx_xyyyy[i] * fz_0;

        tk_xx_xyyyz[i] = -2.0 * ts_0_xyyyz[i] * fbe_0 * fz_0 + tk_0_xyyyz[i] * fe_0 + tk_x_yyyz[i] * fe_0 + tk_x_xyyyz[i] * pa_x[i] + 2.0 * ts_xx_xyyyz[i] * fz_0;

        tk_xx_xyyzz[i] = -2.0 * ts_0_xyyzz[i] * fbe_0 * fz_0 + tk_0_xyyzz[i] * fe_0 + tk_x_yyzz[i] * fe_0 + tk_x_xyyzz[i] * pa_x[i] + 2.0 * ts_xx_xyyzz[i] * fz_0;

        tk_xx_xyzzz[i] = -2.0 * ts_0_xyzzz[i] * fbe_0 * fz_0 + tk_0_xyzzz[i] * fe_0 + tk_x_yzzz[i] * fe_0 + tk_x_xyzzz[i] * pa_x[i] + 2.0 * ts_xx_xyzzz[i] * fz_0;

        tk_xx_xzzzz[i] = -2.0 * ts_0_xzzzz[i] * fbe_0 * fz_0 + tk_0_xzzzz[i] * fe_0 + tk_x_zzzz[i] * fe_0 + tk_x_xzzzz[i] * pa_x[i] + 2.0 * ts_xx_xzzzz[i] * fz_0;

        tk_xx_yyyyy[i] = -2.0 * ts_0_yyyyy[i] * fbe_0 * fz_0 + tk_0_yyyyy[i] * fe_0 + tk_x_yyyyy[i] * pa_x[i] + 2.0 * ts_xx_yyyyy[i] * fz_0;

        tk_xx_yyyyz[i] = -2.0 * ts_0_yyyyz[i] * fbe_0 * fz_0 + tk_0_yyyyz[i] * fe_0 + tk_x_yyyyz[i] * pa_x[i] + 2.0 * ts_xx_yyyyz[i] * fz_0;

        tk_xx_yyyzz[i] = -2.0 * ts_0_yyyzz[i] * fbe_0 * fz_0 + tk_0_yyyzz[i] * fe_0 + tk_x_yyyzz[i] * pa_x[i] + 2.0 * ts_xx_yyyzz[i] * fz_0;

        tk_xx_yyzzz[i] = -2.0 * ts_0_yyzzz[i] * fbe_0 * fz_0 + tk_0_yyzzz[i] * fe_0 + tk_x_yyzzz[i] * pa_x[i] + 2.0 * ts_xx_yyzzz[i] * fz_0;

        tk_xx_yzzzz[i] = -2.0 * ts_0_yzzzz[i] * fbe_0 * fz_0 + tk_0_yzzzz[i] * fe_0 + tk_x_yzzzz[i] * pa_x[i] + 2.0 * ts_xx_yzzzz[i] * fz_0;

        tk_xx_zzzzz[i] = -2.0 * ts_0_zzzzz[i] * fbe_0 * fz_0 + tk_0_zzzzz[i] * fe_0 + tk_x_zzzzz[i] * pa_x[i] + 2.0 * ts_xx_zzzzz[i] * fz_0;
    }

    // Set up 21-42 components of targeted buffer : DH

    auto tk_xy_xxxxx = pbuffer.data(idx_kin_dh + 21);

    auto tk_xy_xxxxy = pbuffer.data(idx_kin_dh + 22);

    auto tk_xy_xxxxz = pbuffer.data(idx_kin_dh + 23);

    auto tk_xy_xxxyy = pbuffer.data(idx_kin_dh + 24);

    auto tk_xy_xxxyz = pbuffer.data(idx_kin_dh + 25);

    auto tk_xy_xxxzz = pbuffer.data(idx_kin_dh + 26);

    auto tk_xy_xxyyy = pbuffer.data(idx_kin_dh + 27);

    auto tk_xy_xxyyz = pbuffer.data(idx_kin_dh + 28);

    auto tk_xy_xxyzz = pbuffer.data(idx_kin_dh + 29);

    auto tk_xy_xxzzz = pbuffer.data(idx_kin_dh + 30);

    auto tk_xy_xyyyy = pbuffer.data(idx_kin_dh + 31);

    auto tk_xy_xyyyz = pbuffer.data(idx_kin_dh + 32);

    auto tk_xy_xyyzz = pbuffer.data(idx_kin_dh + 33);

    auto tk_xy_xyzzz = pbuffer.data(idx_kin_dh + 34);

    auto tk_xy_xzzzz = pbuffer.data(idx_kin_dh + 35);

    auto tk_xy_yyyyy = pbuffer.data(idx_kin_dh + 36);

    auto tk_xy_yyyyz = pbuffer.data(idx_kin_dh + 37);

    auto tk_xy_yyyzz = pbuffer.data(idx_kin_dh + 38);

    auto tk_xy_yyzzz = pbuffer.data(idx_kin_dh + 39);

    auto tk_xy_yzzzz = pbuffer.data(idx_kin_dh + 40);

    auto tk_xy_zzzzz = pbuffer.data(idx_kin_dh + 41);

    #pragma omp simd aligned(pa_x, pa_y, tk_x_xxxxx, tk_x_xxxxz, tk_x_xxxzz, tk_x_xxzzz, tk_x_xzzzz, tk_xy_xxxxx, tk_xy_xxxxy, tk_xy_xxxxz, tk_xy_xxxyy, tk_xy_xxxyz, tk_xy_xxxzz, tk_xy_xxyyy, tk_xy_xxyyz, tk_xy_xxyzz, tk_xy_xxzzz, tk_xy_xyyyy, tk_xy_xyyyz, tk_xy_xyyzz, tk_xy_xyzzz, tk_xy_xzzzz, tk_xy_yyyyy, tk_xy_yyyyz, tk_xy_yyyzz, tk_xy_yyzzz, tk_xy_yzzzz, tk_xy_zzzzz, tk_y_xxxxy, tk_y_xxxy, tk_y_xxxyy, tk_y_xxxyz, tk_y_xxyy, tk_y_xxyyy, tk_y_xxyyz, tk_y_xxyz, tk_y_xxyzz, tk_y_xyyy, tk_y_xyyyy, tk_y_xyyyz, tk_y_xyyz, tk_y_xyyzz, tk_y_xyzz, tk_y_xyzzz, tk_y_yyyy, tk_y_yyyyy, tk_y_yyyyz, tk_y_yyyz, tk_y_yyyzz, tk_y_yyzz, tk_y_yyzzz, tk_y_yzzz, tk_y_yzzzz, tk_y_zzzzz, ts_xy_xxxxx, ts_xy_xxxxy, ts_xy_xxxxz, ts_xy_xxxyy, ts_xy_xxxyz, ts_xy_xxxzz, ts_xy_xxyyy, ts_xy_xxyyz, ts_xy_xxyzz, ts_xy_xxzzz, ts_xy_xyyyy, ts_xy_xyyyz, ts_xy_xyyzz, ts_xy_xyzzz, ts_xy_xzzzz, ts_xy_yyyyy, ts_xy_yyyyz, ts_xy_yyyzz, ts_xy_yyzzz, ts_xy_yzzzz, ts_xy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xy_xxxxx[i] = tk_x_xxxxx[i] * pa_y[i] + 2.0 * ts_xy_xxxxx[i] * fz_0;

        tk_xy_xxxxy[i] = 4.0 * tk_y_xxxy[i] * fe_0 + tk_y_xxxxy[i] * pa_x[i] + 2.0 * ts_xy_xxxxy[i] * fz_0;

        tk_xy_xxxxz[i] = tk_x_xxxxz[i] * pa_y[i] + 2.0 * ts_xy_xxxxz[i] * fz_0;

        tk_xy_xxxyy[i] = 3.0 * tk_y_xxyy[i] * fe_0 + tk_y_xxxyy[i] * pa_x[i] + 2.0 * ts_xy_xxxyy[i] * fz_0;

        tk_xy_xxxyz[i] = 3.0 * tk_y_xxyz[i] * fe_0 + tk_y_xxxyz[i] * pa_x[i] + 2.0 * ts_xy_xxxyz[i] * fz_0;

        tk_xy_xxxzz[i] = tk_x_xxxzz[i] * pa_y[i] + 2.0 * ts_xy_xxxzz[i] * fz_0;

        tk_xy_xxyyy[i] = 2.0 * tk_y_xyyy[i] * fe_0 + tk_y_xxyyy[i] * pa_x[i] + 2.0 * ts_xy_xxyyy[i] * fz_0;

        tk_xy_xxyyz[i] = 2.0 * tk_y_xyyz[i] * fe_0 + tk_y_xxyyz[i] * pa_x[i] + 2.0 * ts_xy_xxyyz[i] * fz_0;

        tk_xy_xxyzz[i] = 2.0 * tk_y_xyzz[i] * fe_0 + tk_y_xxyzz[i] * pa_x[i] + 2.0 * ts_xy_xxyzz[i] * fz_0;

        tk_xy_xxzzz[i] = tk_x_xxzzz[i] * pa_y[i] + 2.0 * ts_xy_xxzzz[i] * fz_0;

        tk_xy_xyyyy[i] = tk_y_yyyy[i] * fe_0 + tk_y_xyyyy[i] * pa_x[i] + 2.0 * ts_xy_xyyyy[i] * fz_0;

        tk_xy_xyyyz[i] = tk_y_yyyz[i] * fe_0 + tk_y_xyyyz[i] * pa_x[i] + 2.0 * ts_xy_xyyyz[i] * fz_0;

        tk_xy_xyyzz[i] = tk_y_yyzz[i] * fe_0 + tk_y_xyyzz[i] * pa_x[i] + 2.0 * ts_xy_xyyzz[i] * fz_0;

        tk_xy_xyzzz[i] = tk_y_yzzz[i] * fe_0 + tk_y_xyzzz[i] * pa_x[i] + 2.0 * ts_xy_xyzzz[i] * fz_0;

        tk_xy_xzzzz[i] = tk_x_xzzzz[i] * pa_y[i] + 2.0 * ts_xy_xzzzz[i] * fz_0;

        tk_xy_yyyyy[i] = tk_y_yyyyy[i] * pa_x[i] + 2.0 * ts_xy_yyyyy[i] * fz_0;

        tk_xy_yyyyz[i] = tk_y_yyyyz[i] * pa_x[i] + 2.0 * ts_xy_yyyyz[i] * fz_0;

        tk_xy_yyyzz[i] = tk_y_yyyzz[i] * pa_x[i] + 2.0 * ts_xy_yyyzz[i] * fz_0;

        tk_xy_yyzzz[i] = tk_y_yyzzz[i] * pa_x[i] + 2.0 * ts_xy_yyzzz[i] * fz_0;

        tk_xy_yzzzz[i] = tk_y_yzzzz[i] * pa_x[i] + 2.0 * ts_xy_yzzzz[i] * fz_0;

        tk_xy_zzzzz[i] = tk_y_zzzzz[i] * pa_x[i] + 2.0 * ts_xy_zzzzz[i] * fz_0;
    }

    // Set up 42-63 components of targeted buffer : DH

    auto tk_xz_xxxxx = pbuffer.data(idx_kin_dh + 42);

    auto tk_xz_xxxxy = pbuffer.data(idx_kin_dh + 43);

    auto tk_xz_xxxxz = pbuffer.data(idx_kin_dh + 44);

    auto tk_xz_xxxyy = pbuffer.data(idx_kin_dh + 45);

    auto tk_xz_xxxyz = pbuffer.data(idx_kin_dh + 46);

    auto tk_xz_xxxzz = pbuffer.data(idx_kin_dh + 47);

    auto tk_xz_xxyyy = pbuffer.data(idx_kin_dh + 48);

    auto tk_xz_xxyyz = pbuffer.data(idx_kin_dh + 49);

    auto tk_xz_xxyzz = pbuffer.data(idx_kin_dh + 50);

    auto tk_xz_xxzzz = pbuffer.data(idx_kin_dh + 51);

    auto tk_xz_xyyyy = pbuffer.data(idx_kin_dh + 52);

    auto tk_xz_xyyyz = pbuffer.data(idx_kin_dh + 53);

    auto tk_xz_xyyzz = pbuffer.data(idx_kin_dh + 54);

    auto tk_xz_xyzzz = pbuffer.data(idx_kin_dh + 55);

    auto tk_xz_xzzzz = pbuffer.data(idx_kin_dh + 56);

    auto tk_xz_yyyyy = pbuffer.data(idx_kin_dh + 57);

    auto tk_xz_yyyyz = pbuffer.data(idx_kin_dh + 58);

    auto tk_xz_yyyzz = pbuffer.data(idx_kin_dh + 59);

    auto tk_xz_yyzzz = pbuffer.data(idx_kin_dh + 60);

    auto tk_xz_yzzzz = pbuffer.data(idx_kin_dh + 61);

    auto tk_xz_zzzzz = pbuffer.data(idx_kin_dh + 62);

    #pragma omp simd aligned(pa_x, pa_z, tk_x_xxxxx, tk_x_xxxxy, tk_x_xxxyy, tk_x_xxyyy, tk_x_xyyyy, tk_xz_xxxxx, tk_xz_xxxxy, tk_xz_xxxxz, tk_xz_xxxyy, tk_xz_xxxyz, tk_xz_xxxzz, tk_xz_xxyyy, tk_xz_xxyyz, tk_xz_xxyzz, tk_xz_xxzzz, tk_xz_xyyyy, tk_xz_xyyyz, tk_xz_xyyzz, tk_xz_xyzzz, tk_xz_xzzzz, tk_xz_yyyyy, tk_xz_yyyyz, tk_xz_yyyzz, tk_xz_yyzzz, tk_xz_yzzzz, tk_xz_zzzzz, tk_z_xxxxz, tk_z_xxxyz, tk_z_xxxz, tk_z_xxxzz, tk_z_xxyyz, tk_z_xxyz, tk_z_xxyzz, tk_z_xxzz, tk_z_xxzzz, tk_z_xyyyz, tk_z_xyyz, tk_z_xyyzz, tk_z_xyzz, tk_z_xyzzz, tk_z_xzzz, tk_z_xzzzz, tk_z_yyyyy, tk_z_yyyyz, tk_z_yyyz, tk_z_yyyzz, tk_z_yyzz, tk_z_yyzzz, tk_z_yzzz, tk_z_yzzzz, tk_z_zzzz, tk_z_zzzzz, ts_xz_xxxxx, ts_xz_xxxxy, ts_xz_xxxxz, ts_xz_xxxyy, ts_xz_xxxyz, ts_xz_xxxzz, ts_xz_xxyyy, ts_xz_xxyyz, ts_xz_xxyzz, ts_xz_xxzzz, ts_xz_xyyyy, ts_xz_xyyyz, ts_xz_xyyzz, ts_xz_xyzzz, ts_xz_xzzzz, ts_xz_yyyyy, ts_xz_yyyyz, ts_xz_yyyzz, ts_xz_yyzzz, ts_xz_yzzzz, ts_xz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xz_xxxxx[i] = tk_x_xxxxx[i] * pa_z[i] + 2.0 * ts_xz_xxxxx[i] * fz_0;

        tk_xz_xxxxy[i] = tk_x_xxxxy[i] * pa_z[i] + 2.0 * ts_xz_xxxxy[i] * fz_0;

        tk_xz_xxxxz[i] = 4.0 * tk_z_xxxz[i] * fe_0 + tk_z_xxxxz[i] * pa_x[i] + 2.0 * ts_xz_xxxxz[i] * fz_0;

        tk_xz_xxxyy[i] = tk_x_xxxyy[i] * pa_z[i] + 2.0 * ts_xz_xxxyy[i] * fz_0;

        tk_xz_xxxyz[i] = 3.0 * tk_z_xxyz[i] * fe_0 + tk_z_xxxyz[i] * pa_x[i] + 2.0 * ts_xz_xxxyz[i] * fz_0;

        tk_xz_xxxzz[i] = 3.0 * tk_z_xxzz[i] * fe_0 + tk_z_xxxzz[i] * pa_x[i] + 2.0 * ts_xz_xxxzz[i] * fz_0;

        tk_xz_xxyyy[i] = tk_x_xxyyy[i] * pa_z[i] + 2.0 * ts_xz_xxyyy[i] * fz_0;

        tk_xz_xxyyz[i] = 2.0 * tk_z_xyyz[i] * fe_0 + tk_z_xxyyz[i] * pa_x[i] + 2.0 * ts_xz_xxyyz[i] * fz_0;

        tk_xz_xxyzz[i] = 2.0 * tk_z_xyzz[i] * fe_0 + tk_z_xxyzz[i] * pa_x[i] + 2.0 * ts_xz_xxyzz[i] * fz_0;

        tk_xz_xxzzz[i] = 2.0 * tk_z_xzzz[i] * fe_0 + tk_z_xxzzz[i] * pa_x[i] + 2.0 * ts_xz_xxzzz[i] * fz_0;

        tk_xz_xyyyy[i] = tk_x_xyyyy[i] * pa_z[i] + 2.0 * ts_xz_xyyyy[i] * fz_0;

        tk_xz_xyyyz[i] = tk_z_yyyz[i] * fe_0 + tk_z_xyyyz[i] * pa_x[i] + 2.0 * ts_xz_xyyyz[i] * fz_0;

        tk_xz_xyyzz[i] = tk_z_yyzz[i] * fe_0 + tk_z_xyyzz[i] * pa_x[i] + 2.0 * ts_xz_xyyzz[i] * fz_0;

        tk_xz_xyzzz[i] = tk_z_yzzz[i] * fe_0 + tk_z_xyzzz[i] * pa_x[i] + 2.0 * ts_xz_xyzzz[i] * fz_0;

        tk_xz_xzzzz[i] = tk_z_zzzz[i] * fe_0 + tk_z_xzzzz[i] * pa_x[i] + 2.0 * ts_xz_xzzzz[i] * fz_0;

        tk_xz_yyyyy[i] = tk_z_yyyyy[i] * pa_x[i] + 2.0 * ts_xz_yyyyy[i] * fz_0;

        tk_xz_yyyyz[i] = tk_z_yyyyz[i] * pa_x[i] + 2.0 * ts_xz_yyyyz[i] * fz_0;

        tk_xz_yyyzz[i] = tk_z_yyyzz[i] * pa_x[i] + 2.0 * ts_xz_yyyzz[i] * fz_0;

        tk_xz_yyzzz[i] = tk_z_yyzzz[i] * pa_x[i] + 2.0 * ts_xz_yyzzz[i] * fz_0;

        tk_xz_yzzzz[i] = tk_z_yzzzz[i] * pa_x[i] + 2.0 * ts_xz_yzzzz[i] * fz_0;

        tk_xz_zzzzz[i] = tk_z_zzzzz[i] * pa_x[i] + 2.0 * ts_xz_zzzzz[i] * fz_0;
    }

    // Set up 63-84 components of targeted buffer : DH

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

    #pragma omp simd aligned(pa_y, tk_0_xxxxx, tk_0_xxxxy, tk_0_xxxxz, tk_0_xxxyy, tk_0_xxxyz, tk_0_xxxzz, tk_0_xxyyy, tk_0_xxyyz, tk_0_xxyzz, tk_0_xxzzz, tk_0_xyyyy, tk_0_xyyyz, tk_0_xyyzz, tk_0_xyzzz, tk_0_xzzzz, tk_0_yyyyy, tk_0_yyyyz, tk_0_yyyzz, tk_0_yyzzz, tk_0_yzzzz, tk_0_zzzzz, tk_y_xxxx, tk_y_xxxxx, tk_y_xxxxy, tk_y_xxxxz, tk_y_xxxy, tk_y_xxxyy, tk_y_xxxyz, tk_y_xxxz, tk_y_xxxzz, tk_y_xxyy, tk_y_xxyyy, tk_y_xxyyz, tk_y_xxyz, tk_y_xxyzz, tk_y_xxzz, tk_y_xxzzz, tk_y_xyyy, tk_y_xyyyy, tk_y_xyyyz, tk_y_xyyz, tk_y_xyyzz, tk_y_xyzz, tk_y_xyzzz, tk_y_xzzz, tk_y_xzzzz, tk_y_yyyy, tk_y_yyyyy, tk_y_yyyyz, tk_y_yyyz, tk_y_yyyzz, tk_y_yyzz, tk_y_yyzzz, tk_y_yzzz, tk_y_yzzzz, tk_y_zzzz, tk_y_zzzzz, tk_yy_xxxxx, tk_yy_xxxxy, tk_yy_xxxxz, tk_yy_xxxyy, tk_yy_xxxyz, tk_yy_xxxzz, tk_yy_xxyyy, tk_yy_xxyyz, tk_yy_xxyzz, tk_yy_xxzzz, tk_yy_xyyyy, tk_yy_xyyyz, tk_yy_xyyzz, tk_yy_xyzzz, tk_yy_xzzzz, tk_yy_yyyyy, tk_yy_yyyyz, tk_yy_yyyzz, tk_yy_yyzzz, tk_yy_yzzzz, tk_yy_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_yy_xxxxx, ts_yy_xxxxy, ts_yy_xxxxz, ts_yy_xxxyy, ts_yy_xxxyz, ts_yy_xxxzz, ts_yy_xxyyy, ts_yy_xxyyz, ts_yy_xxyzz, ts_yy_xxzzz, ts_yy_xyyyy, ts_yy_xyyyz, ts_yy_xyyzz, ts_yy_xyzzz, ts_yy_xzzzz, ts_yy_yyyyy, ts_yy_yyyyz, ts_yy_yyyzz, ts_yy_yyzzz, ts_yy_yzzzz, ts_yy_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yy_xxxxx[i] = -2.0 * ts_0_xxxxx[i] * fbe_0 * fz_0 + tk_0_xxxxx[i] * fe_0 + tk_y_xxxxx[i] * pa_y[i] + 2.0 * ts_yy_xxxxx[i] * fz_0;

        tk_yy_xxxxy[i] = -2.0 * ts_0_xxxxy[i] * fbe_0 * fz_0 + tk_0_xxxxy[i] * fe_0 + tk_y_xxxx[i] * fe_0 + tk_y_xxxxy[i] * pa_y[i] + 2.0 * ts_yy_xxxxy[i] * fz_0;

        tk_yy_xxxxz[i] = -2.0 * ts_0_xxxxz[i] * fbe_0 * fz_0 + tk_0_xxxxz[i] * fe_0 + tk_y_xxxxz[i] * pa_y[i] + 2.0 * ts_yy_xxxxz[i] * fz_0;

        tk_yy_xxxyy[i] = -2.0 * ts_0_xxxyy[i] * fbe_0 * fz_0 + tk_0_xxxyy[i] * fe_0 + 2.0 * tk_y_xxxy[i] * fe_0 + tk_y_xxxyy[i] * pa_y[i] + 2.0 * ts_yy_xxxyy[i] * fz_0;

        tk_yy_xxxyz[i] = -2.0 * ts_0_xxxyz[i] * fbe_0 * fz_0 + tk_0_xxxyz[i] * fe_0 + tk_y_xxxz[i] * fe_0 + tk_y_xxxyz[i] * pa_y[i] + 2.0 * ts_yy_xxxyz[i] * fz_0;

        tk_yy_xxxzz[i] = -2.0 * ts_0_xxxzz[i] * fbe_0 * fz_0 + tk_0_xxxzz[i] * fe_0 + tk_y_xxxzz[i] * pa_y[i] + 2.0 * ts_yy_xxxzz[i] * fz_0;

        tk_yy_xxyyy[i] = -2.0 * ts_0_xxyyy[i] * fbe_0 * fz_0 + tk_0_xxyyy[i] * fe_0 + 3.0 * tk_y_xxyy[i] * fe_0 + tk_y_xxyyy[i] * pa_y[i] + 2.0 * ts_yy_xxyyy[i] * fz_0;

        tk_yy_xxyyz[i] = -2.0 * ts_0_xxyyz[i] * fbe_0 * fz_0 + tk_0_xxyyz[i] * fe_0 + 2.0 * tk_y_xxyz[i] * fe_0 + tk_y_xxyyz[i] * pa_y[i] + 2.0 * ts_yy_xxyyz[i] * fz_0;

        tk_yy_xxyzz[i] = -2.0 * ts_0_xxyzz[i] * fbe_0 * fz_0 + tk_0_xxyzz[i] * fe_0 + tk_y_xxzz[i] * fe_0 + tk_y_xxyzz[i] * pa_y[i] + 2.0 * ts_yy_xxyzz[i] * fz_0;

        tk_yy_xxzzz[i] = -2.0 * ts_0_xxzzz[i] * fbe_0 * fz_0 + tk_0_xxzzz[i] * fe_0 + tk_y_xxzzz[i] * pa_y[i] + 2.0 * ts_yy_xxzzz[i] * fz_0;

        tk_yy_xyyyy[i] = -2.0 * ts_0_xyyyy[i] * fbe_0 * fz_0 + tk_0_xyyyy[i] * fe_0 + 4.0 * tk_y_xyyy[i] * fe_0 + tk_y_xyyyy[i] * pa_y[i] + 2.0 * ts_yy_xyyyy[i] * fz_0;

        tk_yy_xyyyz[i] = -2.0 * ts_0_xyyyz[i] * fbe_0 * fz_0 + tk_0_xyyyz[i] * fe_0 + 3.0 * tk_y_xyyz[i] * fe_0 + tk_y_xyyyz[i] * pa_y[i] + 2.0 * ts_yy_xyyyz[i] * fz_0;

        tk_yy_xyyzz[i] = -2.0 * ts_0_xyyzz[i] * fbe_0 * fz_0 + tk_0_xyyzz[i] * fe_0 + 2.0 * tk_y_xyzz[i] * fe_0 + tk_y_xyyzz[i] * pa_y[i] + 2.0 * ts_yy_xyyzz[i] * fz_0;

        tk_yy_xyzzz[i] = -2.0 * ts_0_xyzzz[i] * fbe_0 * fz_0 + tk_0_xyzzz[i] * fe_0 + tk_y_xzzz[i] * fe_0 + tk_y_xyzzz[i] * pa_y[i] + 2.0 * ts_yy_xyzzz[i] * fz_0;

        tk_yy_xzzzz[i] = -2.0 * ts_0_xzzzz[i] * fbe_0 * fz_0 + tk_0_xzzzz[i] * fe_0 + tk_y_xzzzz[i] * pa_y[i] + 2.0 * ts_yy_xzzzz[i] * fz_0;

        tk_yy_yyyyy[i] = -2.0 * ts_0_yyyyy[i] * fbe_0 * fz_0 + tk_0_yyyyy[i] * fe_0 + 5.0 * tk_y_yyyy[i] * fe_0 + tk_y_yyyyy[i] * pa_y[i] + 2.0 * ts_yy_yyyyy[i] * fz_0;

        tk_yy_yyyyz[i] = -2.0 * ts_0_yyyyz[i] * fbe_0 * fz_0 + tk_0_yyyyz[i] * fe_0 + 4.0 * tk_y_yyyz[i] * fe_0 + tk_y_yyyyz[i] * pa_y[i] + 2.0 * ts_yy_yyyyz[i] * fz_0;

        tk_yy_yyyzz[i] = -2.0 * ts_0_yyyzz[i] * fbe_0 * fz_0 + tk_0_yyyzz[i] * fe_0 + 3.0 * tk_y_yyzz[i] * fe_0 + tk_y_yyyzz[i] * pa_y[i] + 2.0 * ts_yy_yyyzz[i] * fz_0;

        tk_yy_yyzzz[i] = -2.0 * ts_0_yyzzz[i] * fbe_0 * fz_0 + tk_0_yyzzz[i] * fe_0 + 2.0 * tk_y_yzzz[i] * fe_0 + tk_y_yyzzz[i] * pa_y[i] + 2.0 * ts_yy_yyzzz[i] * fz_0;

        tk_yy_yzzzz[i] = -2.0 * ts_0_yzzzz[i] * fbe_0 * fz_0 + tk_0_yzzzz[i] * fe_0 + tk_y_zzzz[i] * fe_0 + tk_y_yzzzz[i] * pa_y[i] + 2.0 * ts_yy_yzzzz[i] * fz_0;

        tk_yy_zzzzz[i] = -2.0 * ts_0_zzzzz[i] * fbe_0 * fz_0 + tk_0_zzzzz[i] * fe_0 + tk_y_zzzzz[i] * pa_y[i] + 2.0 * ts_yy_zzzzz[i] * fz_0;
    }

    // Set up 84-105 components of targeted buffer : DH

    auto tk_yz_xxxxx = pbuffer.data(idx_kin_dh + 84);

    auto tk_yz_xxxxy = pbuffer.data(idx_kin_dh + 85);

    auto tk_yz_xxxxz = pbuffer.data(idx_kin_dh + 86);

    auto tk_yz_xxxyy = pbuffer.data(idx_kin_dh + 87);

    auto tk_yz_xxxyz = pbuffer.data(idx_kin_dh + 88);

    auto tk_yz_xxxzz = pbuffer.data(idx_kin_dh + 89);

    auto tk_yz_xxyyy = pbuffer.data(idx_kin_dh + 90);

    auto tk_yz_xxyyz = pbuffer.data(idx_kin_dh + 91);

    auto tk_yz_xxyzz = pbuffer.data(idx_kin_dh + 92);

    auto tk_yz_xxzzz = pbuffer.data(idx_kin_dh + 93);

    auto tk_yz_xyyyy = pbuffer.data(idx_kin_dh + 94);

    auto tk_yz_xyyyz = pbuffer.data(idx_kin_dh + 95);

    auto tk_yz_xyyzz = pbuffer.data(idx_kin_dh + 96);

    auto tk_yz_xyzzz = pbuffer.data(idx_kin_dh + 97);

    auto tk_yz_xzzzz = pbuffer.data(idx_kin_dh + 98);

    auto tk_yz_yyyyy = pbuffer.data(idx_kin_dh + 99);

    auto tk_yz_yyyyz = pbuffer.data(idx_kin_dh + 100);

    auto tk_yz_yyyzz = pbuffer.data(idx_kin_dh + 101);

    auto tk_yz_yyzzz = pbuffer.data(idx_kin_dh + 102);

    auto tk_yz_yzzzz = pbuffer.data(idx_kin_dh + 103);

    auto tk_yz_zzzzz = pbuffer.data(idx_kin_dh + 104);

    #pragma omp simd aligned(pa_y, pa_z, tk_y_xxxxy, tk_y_xxxyy, tk_y_xxyyy, tk_y_xyyyy, tk_y_yyyyy, tk_yz_xxxxx, tk_yz_xxxxy, tk_yz_xxxxz, tk_yz_xxxyy, tk_yz_xxxyz, tk_yz_xxxzz, tk_yz_xxyyy, tk_yz_xxyyz, tk_yz_xxyzz, tk_yz_xxzzz, tk_yz_xyyyy, tk_yz_xyyyz, tk_yz_xyyzz, tk_yz_xyzzz, tk_yz_xzzzz, tk_yz_yyyyy, tk_yz_yyyyz, tk_yz_yyyzz, tk_yz_yyzzz, tk_yz_yzzzz, tk_yz_zzzzz, tk_z_xxxxx, tk_z_xxxxz, tk_z_xxxyz, tk_z_xxxz, tk_z_xxxzz, tk_z_xxyyz, tk_z_xxyz, tk_z_xxyzz, tk_z_xxzz, tk_z_xxzzz, tk_z_xyyyz, tk_z_xyyz, tk_z_xyyzz, tk_z_xyzz, tk_z_xyzzz, tk_z_xzzz, tk_z_xzzzz, tk_z_yyyyz, tk_z_yyyz, tk_z_yyyzz, tk_z_yyzz, tk_z_yyzzz, tk_z_yzzz, tk_z_yzzzz, tk_z_zzzz, tk_z_zzzzz, ts_yz_xxxxx, ts_yz_xxxxy, ts_yz_xxxxz, ts_yz_xxxyy, ts_yz_xxxyz, ts_yz_xxxzz, ts_yz_xxyyy, ts_yz_xxyyz, ts_yz_xxyzz, ts_yz_xxzzz, ts_yz_xyyyy, ts_yz_xyyyz, ts_yz_xyyzz, ts_yz_xyzzz, ts_yz_xzzzz, ts_yz_yyyyy, ts_yz_yyyyz, ts_yz_yyyzz, ts_yz_yyzzz, ts_yz_yzzzz, ts_yz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yz_xxxxx[i] = tk_z_xxxxx[i] * pa_y[i] + 2.0 * ts_yz_xxxxx[i] * fz_0;

        tk_yz_xxxxy[i] = tk_y_xxxxy[i] * pa_z[i] + 2.0 * ts_yz_xxxxy[i] * fz_0;

        tk_yz_xxxxz[i] = tk_z_xxxxz[i] * pa_y[i] + 2.0 * ts_yz_xxxxz[i] * fz_0;

        tk_yz_xxxyy[i] = tk_y_xxxyy[i] * pa_z[i] + 2.0 * ts_yz_xxxyy[i] * fz_0;

        tk_yz_xxxyz[i] = tk_z_xxxz[i] * fe_0 + tk_z_xxxyz[i] * pa_y[i] + 2.0 * ts_yz_xxxyz[i] * fz_0;

        tk_yz_xxxzz[i] = tk_z_xxxzz[i] * pa_y[i] + 2.0 * ts_yz_xxxzz[i] * fz_0;

        tk_yz_xxyyy[i] = tk_y_xxyyy[i] * pa_z[i] + 2.0 * ts_yz_xxyyy[i] * fz_0;

        tk_yz_xxyyz[i] = 2.0 * tk_z_xxyz[i] * fe_0 + tk_z_xxyyz[i] * pa_y[i] + 2.0 * ts_yz_xxyyz[i] * fz_0;

        tk_yz_xxyzz[i] = tk_z_xxzz[i] * fe_0 + tk_z_xxyzz[i] * pa_y[i] + 2.0 * ts_yz_xxyzz[i] * fz_0;

        tk_yz_xxzzz[i] = tk_z_xxzzz[i] * pa_y[i] + 2.0 * ts_yz_xxzzz[i] * fz_0;

        tk_yz_xyyyy[i] = tk_y_xyyyy[i] * pa_z[i] + 2.0 * ts_yz_xyyyy[i] * fz_0;

        tk_yz_xyyyz[i] = 3.0 * tk_z_xyyz[i] * fe_0 + tk_z_xyyyz[i] * pa_y[i] + 2.0 * ts_yz_xyyyz[i] * fz_0;

        tk_yz_xyyzz[i] = 2.0 * tk_z_xyzz[i] * fe_0 + tk_z_xyyzz[i] * pa_y[i] + 2.0 * ts_yz_xyyzz[i] * fz_0;

        tk_yz_xyzzz[i] = tk_z_xzzz[i] * fe_0 + tk_z_xyzzz[i] * pa_y[i] + 2.0 * ts_yz_xyzzz[i] * fz_0;

        tk_yz_xzzzz[i] = tk_z_xzzzz[i] * pa_y[i] + 2.0 * ts_yz_xzzzz[i] * fz_0;

        tk_yz_yyyyy[i] = tk_y_yyyyy[i] * pa_z[i] + 2.0 * ts_yz_yyyyy[i] * fz_0;

        tk_yz_yyyyz[i] = 4.0 * tk_z_yyyz[i] * fe_0 + tk_z_yyyyz[i] * pa_y[i] + 2.0 * ts_yz_yyyyz[i] * fz_0;

        tk_yz_yyyzz[i] = 3.0 * tk_z_yyzz[i] * fe_0 + tk_z_yyyzz[i] * pa_y[i] + 2.0 * ts_yz_yyyzz[i] * fz_0;

        tk_yz_yyzzz[i] = 2.0 * tk_z_yzzz[i] * fe_0 + tk_z_yyzzz[i] * pa_y[i] + 2.0 * ts_yz_yyzzz[i] * fz_0;

        tk_yz_yzzzz[i] = tk_z_zzzz[i] * fe_0 + tk_z_yzzzz[i] * pa_y[i] + 2.0 * ts_yz_yzzzz[i] * fz_0;

        tk_yz_zzzzz[i] = tk_z_zzzzz[i] * pa_y[i] + 2.0 * ts_yz_zzzzz[i] * fz_0;
    }

    // Set up 105-126 components of targeted buffer : DH

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

    #pragma omp simd aligned(pa_z, tk_0_xxxxx, tk_0_xxxxy, tk_0_xxxxz, tk_0_xxxyy, tk_0_xxxyz, tk_0_xxxzz, tk_0_xxyyy, tk_0_xxyyz, tk_0_xxyzz, tk_0_xxzzz, tk_0_xyyyy, tk_0_xyyyz, tk_0_xyyzz, tk_0_xyzzz, tk_0_xzzzz, tk_0_yyyyy, tk_0_yyyyz, tk_0_yyyzz, tk_0_yyzzz, tk_0_yzzzz, tk_0_zzzzz, tk_z_xxxx, tk_z_xxxxx, tk_z_xxxxy, tk_z_xxxxz, tk_z_xxxy, tk_z_xxxyy, tk_z_xxxyz, tk_z_xxxz, tk_z_xxxzz, tk_z_xxyy, tk_z_xxyyy, tk_z_xxyyz, tk_z_xxyz, tk_z_xxyzz, tk_z_xxzz, tk_z_xxzzz, tk_z_xyyy, tk_z_xyyyy, tk_z_xyyyz, tk_z_xyyz, tk_z_xyyzz, tk_z_xyzz, tk_z_xyzzz, tk_z_xzzz, tk_z_xzzzz, tk_z_yyyy, tk_z_yyyyy, tk_z_yyyyz, tk_z_yyyz, tk_z_yyyzz, tk_z_yyzz, tk_z_yyzzz, tk_z_yzzz, tk_z_yzzzz, tk_z_zzzz, tk_z_zzzzz, tk_zz_xxxxx, tk_zz_xxxxy, tk_zz_xxxxz, tk_zz_xxxyy, tk_zz_xxxyz, tk_zz_xxxzz, tk_zz_xxyyy, tk_zz_xxyyz, tk_zz_xxyzz, tk_zz_xxzzz, tk_zz_xyyyy, tk_zz_xyyyz, tk_zz_xyyzz, tk_zz_xyzzz, tk_zz_xzzzz, tk_zz_yyyyy, tk_zz_yyyyz, tk_zz_yyyzz, tk_zz_yyzzz, tk_zz_yzzzz, tk_zz_zzzzz, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxzz, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzzz, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzzzz, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyzz, ts_0_yyzzz, ts_0_yzzzz, ts_0_zzzzz, ts_zz_xxxxx, ts_zz_xxxxy, ts_zz_xxxxz, ts_zz_xxxyy, ts_zz_xxxyz, ts_zz_xxxzz, ts_zz_xxyyy, ts_zz_xxyyz, ts_zz_xxyzz, ts_zz_xxzzz, ts_zz_xyyyy, ts_zz_xyyyz, ts_zz_xyyzz, ts_zz_xyzzz, ts_zz_xzzzz, ts_zz_yyyyy, ts_zz_yyyyz, ts_zz_yyyzz, ts_zz_yyzzz, ts_zz_yzzzz, ts_zz_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zz_xxxxx[i] = -2.0 * ts_0_xxxxx[i] * fbe_0 * fz_0 + tk_0_xxxxx[i] * fe_0 + tk_z_xxxxx[i] * pa_z[i] + 2.0 * ts_zz_xxxxx[i] * fz_0;

        tk_zz_xxxxy[i] = -2.0 * ts_0_xxxxy[i] * fbe_0 * fz_0 + tk_0_xxxxy[i] * fe_0 + tk_z_xxxxy[i] * pa_z[i] + 2.0 * ts_zz_xxxxy[i] * fz_0;

        tk_zz_xxxxz[i] = -2.0 * ts_0_xxxxz[i] * fbe_0 * fz_0 + tk_0_xxxxz[i] * fe_0 + tk_z_xxxx[i] * fe_0 + tk_z_xxxxz[i] * pa_z[i] + 2.0 * ts_zz_xxxxz[i] * fz_0;

        tk_zz_xxxyy[i] = -2.0 * ts_0_xxxyy[i] * fbe_0 * fz_0 + tk_0_xxxyy[i] * fe_0 + tk_z_xxxyy[i] * pa_z[i] + 2.0 * ts_zz_xxxyy[i] * fz_0;

        tk_zz_xxxyz[i] = -2.0 * ts_0_xxxyz[i] * fbe_0 * fz_0 + tk_0_xxxyz[i] * fe_0 + tk_z_xxxy[i] * fe_0 + tk_z_xxxyz[i] * pa_z[i] + 2.0 * ts_zz_xxxyz[i] * fz_0;

        tk_zz_xxxzz[i] = -2.0 * ts_0_xxxzz[i] * fbe_0 * fz_0 + tk_0_xxxzz[i] * fe_0 + 2.0 * tk_z_xxxz[i] * fe_0 + tk_z_xxxzz[i] * pa_z[i] + 2.0 * ts_zz_xxxzz[i] * fz_0;

        tk_zz_xxyyy[i] = -2.0 * ts_0_xxyyy[i] * fbe_0 * fz_0 + tk_0_xxyyy[i] * fe_0 + tk_z_xxyyy[i] * pa_z[i] + 2.0 * ts_zz_xxyyy[i] * fz_0;

        tk_zz_xxyyz[i] = -2.0 * ts_0_xxyyz[i] * fbe_0 * fz_0 + tk_0_xxyyz[i] * fe_0 + tk_z_xxyy[i] * fe_0 + tk_z_xxyyz[i] * pa_z[i] + 2.0 * ts_zz_xxyyz[i] * fz_0;

        tk_zz_xxyzz[i] = -2.0 * ts_0_xxyzz[i] * fbe_0 * fz_0 + tk_0_xxyzz[i] * fe_0 + 2.0 * tk_z_xxyz[i] * fe_0 + tk_z_xxyzz[i] * pa_z[i] + 2.0 * ts_zz_xxyzz[i] * fz_0;

        tk_zz_xxzzz[i] = -2.0 * ts_0_xxzzz[i] * fbe_0 * fz_0 + tk_0_xxzzz[i] * fe_0 + 3.0 * tk_z_xxzz[i] * fe_0 + tk_z_xxzzz[i] * pa_z[i] + 2.0 * ts_zz_xxzzz[i] * fz_0;

        tk_zz_xyyyy[i] = -2.0 * ts_0_xyyyy[i] * fbe_0 * fz_0 + tk_0_xyyyy[i] * fe_0 + tk_z_xyyyy[i] * pa_z[i] + 2.0 * ts_zz_xyyyy[i] * fz_0;

        tk_zz_xyyyz[i] = -2.0 * ts_0_xyyyz[i] * fbe_0 * fz_0 + tk_0_xyyyz[i] * fe_0 + tk_z_xyyy[i] * fe_0 + tk_z_xyyyz[i] * pa_z[i] + 2.0 * ts_zz_xyyyz[i] * fz_0;

        tk_zz_xyyzz[i] = -2.0 * ts_0_xyyzz[i] * fbe_0 * fz_0 + tk_0_xyyzz[i] * fe_0 + 2.0 * tk_z_xyyz[i] * fe_0 + tk_z_xyyzz[i] * pa_z[i] + 2.0 * ts_zz_xyyzz[i] * fz_0;

        tk_zz_xyzzz[i] = -2.0 * ts_0_xyzzz[i] * fbe_0 * fz_0 + tk_0_xyzzz[i] * fe_0 + 3.0 * tk_z_xyzz[i] * fe_0 + tk_z_xyzzz[i] * pa_z[i] + 2.0 * ts_zz_xyzzz[i] * fz_0;

        tk_zz_xzzzz[i] = -2.0 * ts_0_xzzzz[i] * fbe_0 * fz_0 + tk_0_xzzzz[i] * fe_0 + 4.0 * tk_z_xzzz[i] * fe_0 + tk_z_xzzzz[i] * pa_z[i] + 2.0 * ts_zz_xzzzz[i] * fz_0;

        tk_zz_yyyyy[i] = -2.0 * ts_0_yyyyy[i] * fbe_0 * fz_0 + tk_0_yyyyy[i] * fe_0 + tk_z_yyyyy[i] * pa_z[i] + 2.0 * ts_zz_yyyyy[i] * fz_0;

        tk_zz_yyyyz[i] = -2.0 * ts_0_yyyyz[i] * fbe_0 * fz_0 + tk_0_yyyyz[i] * fe_0 + tk_z_yyyy[i] * fe_0 + tk_z_yyyyz[i] * pa_z[i] + 2.0 * ts_zz_yyyyz[i] * fz_0;

        tk_zz_yyyzz[i] = -2.0 * ts_0_yyyzz[i] * fbe_0 * fz_0 + tk_0_yyyzz[i] * fe_0 + 2.0 * tk_z_yyyz[i] * fe_0 + tk_z_yyyzz[i] * pa_z[i] + 2.0 * ts_zz_yyyzz[i] * fz_0;

        tk_zz_yyzzz[i] = -2.0 * ts_0_yyzzz[i] * fbe_0 * fz_0 + tk_0_yyzzz[i] * fe_0 + 3.0 * tk_z_yyzz[i] * fe_0 + tk_z_yyzzz[i] * pa_z[i] + 2.0 * ts_zz_yyzzz[i] * fz_0;

        tk_zz_yzzzz[i] = -2.0 * ts_0_yzzzz[i] * fbe_0 * fz_0 + tk_0_yzzzz[i] * fe_0 + 4.0 * tk_z_yzzz[i] * fe_0 + tk_z_yzzzz[i] * pa_z[i] + 2.0 * ts_zz_yzzzz[i] * fz_0;

        tk_zz_zzzzz[i] = -2.0 * ts_0_zzzzz[i] * fbe_0 * fz_0 + tk_0_zzzzz[i] * fe_0 + 5.0 * tk_z_zzzz[i] * fe_0 + tk_z_zzzzz[i] * pa_z[i] + 2.0 * ts_zz_zzzzz[i] * fz_0;
    }

}

} // kinrec namespace

