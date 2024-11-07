#include "KineticEnergyPrimRecFI.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_fi(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_fi,
                            const size_t              idx_ovl_pi,
                            const size_t              idx_kin_pi,
                            const size_t              idx_kin_dh,
                            const size_t              idx_kin_di,
                            const size_t              idx_ovl_fi,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : PI

    auto ts_x_xxxxxx = pbuffer.data(idx_ovl_pi);

    auto ts_x_xxxxxy = pbuffer.data(idx_ovl_pi + 1);

    auto ts_x_xxxxxz = pbuffer.data(idx_ovl_pi + 2);

    auto ts_x_xxxxyy = pbuffer.data(idx_ovl_pi + 3);

    auto ts_x_xxxxyz = pbuffer.data(idx_ovl_pi + 4);

    auto ts_x_xxxxzz = pbuffer.data(idx_ovl_pi + 5);

    auto ts_x_xxxyyy = pbuffer.data(idx_ovl_pi + 6);

    auto ts_x_xxxyyz = pbuffer.data(idx_ovl_pi + 7);

    auto ts_x_xxxyzz = pbuffer.data(idx_ovl_pi + 8);

    auto ts_x_xxxzzz = pbuffer.data(idx_ovl_pi + 9);

    auto ts_x_xxyyyy = pbuffer.data(idx_ovl_pi + 10);

    auto ts_x_xxyyyz = pbuffer.data(idx_ovl_pi + 11);

    auto ts_x_xxyyzz = pbuffer.data(idx_ovl_pi + 12);

    auto ts_x_xxyzzz = pbuffer.data(idx_ovl_pi + 13);

    auto ts_x_xxzzzz = pbuffer.data(idx_ovl_pi + 14);

    auto ts_x_xyyyyy = pbuffer.data(idx_ovl_pi + 15);

    auto ts_x_xyyyyz = pbuffer.data(idx_ovl_pi + 16);

    auto ts_x_xyyyzz = pbuffer.data(idx_ovl_pi + 17);

    auto ts_x_xyyzzz = pbuffer.data(idx_ovl_pi + 18);

    auto ts_x_xyzzzz = pbuffer.data(idx_ovl_pi + 19);

    auto ts_x_xzzzzz = pbuffer.data(idx_ovl_pi + 20);

    auto ts_x_yyyyyy = pbuffer.data(idx_ovl_pi + 21);

    auto ts_x_yyyyyz = pbuffer.data(idx_ovl_pi + 22);

    auto ts_x_yyyyzz = pbuffer.data(idx_ovl_pi + 23);

    auto ts_x_yyyzzz = pbuffer.data(idx_ovl_pi + 24);

    auto ts_x_yyzzzz = pbuffer.data(idx_ovl_pi + 25);

    auto ts_x_yzzzzz = pbuffer.data(idx_ovl_pi + 26);

    auto ts_x_zzzzzz = pbuffer.data(idx_ovl_pi + 27);

    auto ts_y_xxxxxx = pbuffer.data(idx_ovl_pi + 28);

    auto ts_y_xxxxxy = pbuffer.data(idx_ovl_pi + 29);

    auto ts_y_xxxxxz = pbuffer.data(idx_ovl_pi + 30);

    auto ts_y_xxxxyy = pbuffer.data(idx_ovl_pi + 31);

    auto ts_y_xxxxyz = pbuffer.data(idx_ovl_pi + 32);

    auto ts_y_xxxxzz = pbuffer.data(idx_ovl_pi + 33);

    auto ts_y_xxxyyy = pbuffer.data(idx_ovl_pi + 34);

    auto ts_y_xxxyyz = pbuffer.data(idx_ovl_pi + 35);

    auto ts_y_xxxyzz = pbuffer.data(idx_ovl_pi + 36);

    auto ts_y_xxxzzz = pbuffer.data(idx_ovl_pi + 37);

    auto ts_y_xxyyyy = pbuffer.data(idx_ovl_pi + 38);

    auto ts_y_xxyyyz = pbuffer.data(idx_ovl_pi + 39);

    auto ts_y_xxyyzz = pbuffer.data(idx_ovl_pi + 40);

    auto ts_y_xxyzzz = pbuffer.data(idx_ovl_pi + 41);

    auto ts_y_xxzzzz = pbuffer.data(idx_ovl_pi + 42);

    auto ts_y_xyyyyy = pbuffer.data(idx_ovl_pi + 43);

    auto ts_y_xyyyyz = pbuffer.data(idx_ovl_pi + 44);

    auto ts_y_xyyyzz = pbuffer.data(idx_ovl_pi + 45);

    auto ts_y_xyyzzz = pbuffer.data(idx_ovl_pi + 46);

    auto ts_y_xyzzzz = pbuffer.data(idx_ovl_pi + 47);

    auto ts_y_xzzzzz = pbuffer.data(idx_ovl_pi + 48);

    auto ts_y_yyyyyy = pbuffer.data(idx_ovl_pi + 49);

    auto ts_y_yyyyyz = pbuffer.data(idx_ovl_pi + 50);

    auto ts_y_yyyyzz = pbuffer.data(idx_ovl_pi + 51);

    auto ts_y_yyyzzz = pbuffer.data(idx_ovl_pi + 52);

    auto ts_y_yyzzzz = pbuffer.data(idx_ovl_pi + 53);

    auto ts_y_yzzzzz = pbuffer.data(idx_ovl_pi + 54);

    auto ts_y_zzzzzz = pbuffer.data(idx_ovl_pi + 55);

    auto ts_z_xxxxxx = pbuffer.data(idx_ovl_pi + 56);

    auto ts_z_xxxxxy = pbuffer.data(idx_ovl_pi + 57);

    auto ts_z_xxxxxz = pbuffer.data(idx_ovl_pi + 58);

    auto ts_z_xxxxyy = pbuffer.data(idx_ovl_pi + 59);

    auto ts_z_xxxxyz = pbuffer.data(idx_ovl_pi + 60);

    auto ts_z_xxxxzz = pbuffer.data(idx_ovl_pi + 61);

    auto ts_z_xxxyyy = pbuffer.data(idx_ovl_pi + 62);

    auto ts_z_xxxyyz = pbuffer.data(idx_ovl_pi + 63);

    auto ts_z_xxxyzz = pbuffer.data(idx_ovl_pi + 64);

    auto ts_z_xxxzzz = pbuffer.data(idx_ovl_pi + 65);

    auto ts_z_xxyyyy = pbuffer.data(idx_ovl_pi + 66);

    auto ts_z_xxyyyz = pbuffer.data(idx_ovl_pi + 67);

    auto ts_z_xxyyzz = pbuffer.data(idx_ovl_pi + 68);

    auto ts_z_xxyzzz = pbuffer.data(idx_ovl_pi + 69);

    auto ts_z_xxzzzz = pbuffer.data(idx_ovl_pi + 70);

    auto ts_z_xyyyyy = pbuffer.data(idx_ovl_pi + 71);

    auto ts_z_xyyyyz = pbuffer.data(idx_ovl_pi + 72);

    auto ts_z_xyyyzz = pbuffer.data(idx_ovl_pi + 73);

    auto ts_z_xyyzzz = pbuffer.data(idx_ovl_pi + 74);

    auto ts_z_xyzzzz = pbuffer.data(idx_ovl_pi + 75);

    auto ts_z_xzzzzz = pbuffer.data(idx_ovl_pi + 76);

    auto ts_z_yyyyyy = pbuffer.data(idx_ovl_pi + 77);

    auto ts_z_yyyyyz = pbuffer.data(idx_ovl_pi + 78);

    auto ts_z_yyyyzz = pbuffer.data(idx_ovl_pi + 79);

    auto ts_z_yyyzzz = pbuffer.data(idx_ovl_pi + 80);

    auto ts_z_yyzzzz = pbuffer.data(idx_ovl_pi + 81);

    auto ts_z_yzzzzz = pbuffer.data(idx_ovl_pi + 82);

    auto ts_z_zzzzzz = pbuffer.data(idx_ovl_pi + 83);

    // Set up components of auxiliary buffer : PI

    auto tk_x_xxxxxx = pbuffer.data(idx_kin_pi);

    auto tk_x_xxxxxy = pbuffer.data(idx_kin_pi + 1);

    auto tk_x_xxxxxz = pbuffer.data(idx_kin_pi + 2);

    auto tk_x_xxxxyy = pbuffer.data(idx_kin_pi + 3);

    auto tk_x_xxxxyz = pbuffer.data(idx_kin_pi + 4);

    auto tk_x_xxxxzz = pbuffer.data(idx_kin_pi + 5);

    auto tk_x_xxxyyy = pbuffer.data(idx_kin_pi + 6);

    auto tk_x_xxxyyz = pbuffer.data(idx_kin_pi + 7);

    auto tk_x_xxxyzz = pbuffer.data(idx_kin_pi + 8);

    auto tk_x_xxxzzz = pbuffer.data(idx_kin_pi + 9);

    auto tk_x_xxyyyy = pbuffer.data(idx_kin_pi + 10);

    auto tk_x_xxyyyz = pbuffer.data(idx_kin_pi + 11);

    auto tk_x_xxyyzz = pbuffer.data(idx_kin_pi + 12);

    auto tk_x_xxyzzz = pbuffer.data(idx_kin_pi + 13);

    auto tk_x_xxzzzz = pbuffer.data(idx_kin_pi + 14);

    auto tk_x_xyyyyy = pbuffer.data(idx_kin_pi + 15);

    auto tk_x_xyyyyz = pbuffer.data(idx_kin_pi + 16);

    auto tk_x_xyyyzz = pbuffer.data(idx_kin_pi + 17);

    auto tk_x_xyyzzz = pbuffer.data(idx_kin_pi + 18);

    auto tk_x_xyzzzz = pbuffer.data(idx_kin_pi + 19);

    auto tk_x_xzzzzz = pbuffer.data(idx_kin_pi + 20);

    auto tk_x_yyyyyy = pbuffer.data(idx_kin_pi + 21);

    auto tk_x_yyyyyz = pbuffer.data(idx_kin_pi + 22);

    auto tk_x_yyyyzz = pbuffer.data(idx_kin_pi + 23);

    auto tk_x_yyyzzz = pbuffer.data(idx_kin_pi + 24);

    auto tk_x_yyzzzz = pbuffer.data(idx_kin_pi + 25);

    auto tk_x_yzzzzz = pbuffer.data(idx_kin_pi + 26);

    auto tk_x_zzzzzz = pbuffer.data(idx_kin_pi + 27);

    auto tk_y_xxxxxx = pbuffer.data(idx_kin_pi + 28);

    auto tk_y_xxxxxy = pbuffer.data(idx_kin_pi + 29);

    auto tk_y_xxxxxz = pbuffer.data(idx_kin_pi + 30);

    auto tk_y_xxxxyy = pbuffer.data(idx_kin_pi + 31);

    auto tk_y_xxxxyz = pbuffer.data(idx_kin_pi + 32);

    auto tk_y_xxxxzz = pbuffer.data(idx_kin_pi + 33);

    auto tk_y_xxxyyy = pbuffer.data(idx_kin_pi + 34);

    auto tk_y_xxxyyz = pbuffer.data(idx_kin_pi + 35);

    auto tk_y_xxxyzz = pbuffer.data(idx_kin_pi + 36);

    auto tk_y_xxxzzz = pbuffer.data(idx_kin_pi + 37);

    auto tk_y_xxyyyy = pbuffer.data(idx_kin_pi + 38);

    auto tk_y_xxyyyz = pbuffer.data(idx_kin_pi + 39);

    auto tk_y_xxyyzz = pbuffer.data(idx_kin_pi + 40);

    auto tk_y_xxyzzz = pbuffer.data(idx_kin_pi + 41);

    auto tk_y_xxzzzz = pbuffer.data(idx_kin_pi + 42);

    auto tk_y_xyyyyy = pbuffer.data(idx_kin_pi + 43);

    auto tk_y_xyyyyz = pbuffer.data(idx_kin_pi + 44);

    auto tk_y_xyyyzz = pbuffer.data(idx_kin_pi + 45);

    auto tk_y_xyyzzz = pbuffer.data(idx_kin_pi + 46);

    auto tk_y_xyzzzz = pbuffer.data(idx_kin_pi + 47);

    auto tk_y_xzzzzz = pbuffer.data(idx_kin_pi + 48);

    auto tk_y_yyyyyy = pbuffer.data(idx_kin_pi + 49);

    auto tk_y_yyyyyz = pbuffer.data(idx_kin_pi + 50);

    auto tk_y_yyyyzz = pbuffer.data(idx_kin_pi + 51);

    auto tk_y_yyyzzz = pbuffer.data(idx_kin_pi + 52);

    auto tk_y_yyzzzz = pbuffer.data(idx_kin_pi + 53);

    auto tk_y_yzzzzz = pbuffer.data(idx_kin_pi + 54);

    auto tk_y_zzzzzz = pbuffer.data(idx_kin_pi + 55);

    auto tk_z_xxxxxx = pbuffer.data(idx_kin_pi + 56);

    auto tk_z_xxxxxy = pbuffer.data(idx_kin_pi + 57);

    auto tk_z_xxxxxz = pbuffer.data(idx_kin_pi + 58);

    auto tk_z_xxxxyy = pbuffer.data(idx_kin_pi + 59);

    auto tk_z_xxxxyz = pbuffer.data(idx_kin_pi + 60);

    auto tk_z_xxxxzz = pbuffer.data(idx_kin_pi + 61);

    auto tk_z_xxxyyy = pbuffer.data(idx_kin_pi + 62);

    auto tk_z_xxxyyz = pbuffer.data(idx_kin_pi + 63);

    auto tk_z_xxxyzz = pbuffer.data(idx_kin_pi + 64);

    auto tk_z_xxxzzz = pbuffer.data(idx_kin_pi + 65);

    auto tk_z_xxyyyy = pbuffer.data(idx_kin_pi + 66);

    auto tk_z_xxyyyz = pbuffer.data(idx_kin_pi + 67);

    auto tk_z_xxyyzz = pbuffer.data(idx_kin_pi + 68);

    auto tk_z_xxyzzz = pbuffer.data(idx_kin_pi + 69);

    auto tk_z_xxzzzz = pbuffer.data(idx_kin_pi + 70);

    auto tk_z_xyyyyy = pbuffer.data(idx_kin_pi + 71);

    auto tk_z_xyyyyz = pbuffer.data(idx_kin_pi + 72);

    auto tk_z_xyyyzz = pbuffer.data(idx_kin_pi + 73);

    auto tk_z_xyyzzz = pbuffer.data(idx_kin_pi + 74);

    auto tk_z_xyzzzz = pbuffer.data(idx_kin_pi + 75);

    auto tk_z_xzzzzz = pbuffer.data(idx_kin_pi + 76);

    auto tk_z_yyyyyy = pbuffer.data(idx_kin_pi + 77);

    auto tk_z_yyyyyz = pbuffer.data(idx_kin_pi + 78);

    auto tk_z_yyyyzz = pbuffer.data(idx_kin_pi + 79);

    auto tk_z_yyyzzz = pbuffer.data(idx_kin_pi + 80);

    auto tk_z_yyzzzz = pbuffer.data(idx_kin_pi + 81);

    auto tk_z_yzzzzz = pbuffer.data(idx_kin_pi + 82);

    auto tk_z_zzzzzz = pbuffer.data(idx_kin_pi + 83);

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

    auto tk_yz_yyyyz = pbuffer.data(idx_kin_dh + 100);

    auto tk_yz_yyyzz = pbuffer.data(idx_kin_dh + 101);

    auto tk_yz_yyzzz = pbuffer.data(idx_kin_dh + 102);

    auto tk_yz_yzzzz = pbuffer.data(idx_kin_dh + 103);

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

    // Set up components of auxiliary buffer : DI

    auto tk_xx_xxxxxx = pbuffer.data(idx_kin_di);

    auto tk_xx_xxxxxy = pbuffer.data(idx_kin_di + 1);

    auto tk_xx_xxxxxz = pbuffer.data(idx_kin_di + 2);

    auto tk_xx_xxxxyy = pbuffer.data(idx_kin_di + 3);

    auto tk_xx_xxxxyz = pbuffer.data(idx_kin_di + 4);

    auto tk_xx_xxxxzz = pbuffer.data(idx_kin_di + 5);

    auto tk_xx_xxxyyy = pbuffer.data(idx_kin_di + 6);

    auto tk_xx_xxxyyz = pbuffer.data(idx_kin_di + 7);

    auto tk_xx_xxxyzz = pbuffer.data(idx_kin_di + 8);

    auto tk_xx_xxxzzz = pbuffer.data(idx_kin_di + 9);

    auto tk_xx_xxyyyy = pbuffer.data(idx_kin_di + 10);

    auto tk_xx_xxyyyz = pbuffer.data(idx_kin_di + 11);

    auto tk_xx_xxyyzz = pbuffer.data(idx_kin_di + 12);

    auto tk_xx_xxyzzz = pbuffer.data(idx_kin_di + 13);

    auto tk_xx_xxzzzz = pbuffer.data(idx_kin_di + 14);

    auto tk_xx_xyyyyy = pbuffer.data(idx_kin_di + 15);

    auto tk_xx_xyyyyz = pbuffer.data(idx_kin_di + 16);

    auto tk_xx_xyyyzz = pbuffer.data(idx_kin_di + 17);

    auto tk_xx_xyyzzz = pbuffer.data(idx_kin_di + 18);

    auto tk_xx_xyzzzz = pbuffer.data(idx_kin_di + 19);

    auto tk_xx_xzzzzz = pbuffer.data(idx_kin_di + 20);

    auto tk_xx_yyyyyy = pbuffer.data(idx_kin_di + 21);

    auto tk_xx_yyyyyz = pbuffer.data(idx_kin_di + 22);

    auto tk_xx_yyyyzz = pbuffer.data(idx_kin_di + 23);

    auto tk_xx_yyyzzz = pbuffer.data(idx_kin_di + 24);

    auto tk_xx_yyzzzz = pbuffer.data(idx_kin_di + 25);

    auto tk_xx_yzzzzz = pbuffer.data(idx_kin_di + 26);

    auto tk_xx_zzzzzz = pbuffer.data(idx_kin_di + 27);

    auto tk_xy_xxxxxy = pbuffer.data(idx_kin_di + 29);

    auto tk_xy_xxxxyy = pbuffer.data(idx_kin_di + 31);

    auto tk_xy_xxxyyy = pbuffer.data(idx_kin_di + 34);

    auto tk_xy_xxyyyy = pbuffer.data(idx_kin_di + 38);

    auto tk_xy_xyyyyy = pbuffer.data(idx_kin_di + 43);

    auto tk_xz_xxxxxx = pbuffer.data(idx_kin_di + 56);

    auto tk_xz_xxxxxz = pbuffer.data(idx_kin_di + 58);

    auto tk_xz_xxxxzz = pbuffer.data(idx_kin_di + 61);

    auto tk_xz_xxxzzz = pbuffer.data(idx_kin_di + 65);

    auto tk_xz_xxzzzz = pbuffer.data(idx_kin_di + 70);

    auto tk_xz_xzzzzz = pbuffer.data(idx_kin_di + 76);

    auto tk_yy_xxxxxx = pbuffer.data(idx_kin_di + 84);

    auto tk_yy_xxxxxy = pbuffer.data(idx_kin_di + 85);

    auto tk_yy_xxxxxz = pbuffer.data(idx_kin_di + 86);

    auto tk_yy_xxxxyy = pbuffer.data(idx_kin_di + 87);

    auto tk_yy_xxxxyz = pbuffer.data(idx_kin_di + 88);

    auto tk_yy_xxxxzz = pbuffer.data(idx_kin_di + 89);

    auto tk_yy_xxxyyy = pbuffer.data(idx_kin_di + 90);

    auto tk_yy_xxxyyz = pbuffer.data(idx_kin_di + 91);

    auto tk_yy_xxxyzz = pbuffer.data(idx_kin_di + 92);

    auto tk_yy_xxxzzz = pbuffer.data(idx_kin_di + 93);

    auto tk_yy_xxyyyy = pbuffer.data(idx_kin_di + 94);

    auto tk_yy_xxyyyz = pbuffer.data(idx_kin_di + 95);

    auto tk_yy_xxyyzz = pbuffer.data(idx_kin_di + 96);

    auto tk_yy_xxyzzz = pbuffer.data(idx_kin_di + 97);

    auto tk_yy_xxzzzz = pbuffer.data(idx_kin_di + 98);

    auto tk_yy_xyyyyy = pbuffer.data(idx_kin_di + 99);

    auto tk_yy_xyyyyz = pbuffer.data(idx_kin_di + 100);

    auto tk_yy_xyyyzz = pbuffer.data(idx_kin_di + 101);

    auto tk_yy_xyyzzz = pbuffer.data(idx_kin_di + 102);

    auto tk_yy_xyzzzz = pbuffer.data(idx_kin_di + 103);

    auto tk_yy_xzzzzz = pbuffer.data(idx_kin_di + 104);

    auto tk_yy_yyyyyy = pbuffer.data(idx_kin_di + 105);

    auto tk_yy_yyyyyz = pbuffer.data(idx_kin_di + 106);

    auto tk_yy_yyyyzz = pbuffer.data(idx_kin_di + 107);

    auto tk_yy_yyyzzz = pbuffer.data(idx_kin_di + 108);

    auto tk_yy_yyzzzz = pbuffer.data(idx_kin_di + 109);

    auto tk_yy_yzzzzz = pbuffer.data(idx_kin_di + 110);

    auto tk_yy_zzzzzz = pbuffer.data(idx_kin_di + 111);

    auto tk_yz_xxxxyz = pbuffer.data(idx_kin_di + 116);

    auto tk_yz_xxxyyz = pbuffer.data(idx_kin_di + 119);

    auto tk_yz_xxxyzz = pbuffer.data(idx_kin_di + 120);

    auto tk_yz_xxyyyz = pbuffer.data(idx_kin_di + 123);

    auto tk_yz_xxyyzz = pbuffer.data(idx_kin_di + 124);

    auto tk_yz_xxyzzz = pbuffer.data(idx_kin_di + 125);

    auto tk_yz_xyyyyz = pbuffer.data(idx_kin_di + 128);

    auto tk_yz_xyyyzz = pbuffer.data(idx_kin_di + 129);

    auto tk_yz_xyyzzz = pbuffer.data(idx_kin_di + 130);

    auto tk_yz_xyzzzz = pbuffer.data(idx_kin_di + 131);

    auto tk_yz_yyyyyy = pbuffer.data(idx_kin_di + 133);

    auto tk_yz_yyyyyz = pbuffer.data(idx_kin_di + 134);

    auto tk_yz_yyyyzz = pbuffer.data(idx_kin_di + 135);

    auto tk_yz_yyyzzz = pbuffer.data(idx_kin_di + 136);

    auto tk_yz_yyzzzz = pbuffer.data(idx_kin_di + 137);

    auto tk_yz_yzzzzz = pbuffer.data(idx_kin_di + 138);

    auto tk_yz_zzzzzz = pbuffer.data(idx_kin_di + 139);

    auto tk_zz_xxxxxx = pbuffer.data(idx_kin_di + 140);

    auto tk_zz_xxxxxy = pbuffer.data(idx_kin_di + 141);

    auto tk_zz_xxxxxz = pbuffer.data(idx_kin_di + 142);

    auto tk_zz_xxxxyy = pbuffer.data(idx_kin_di + 143);

    auto tk_zz_xxxxyz = pbuffer.data(idx_kin_di + 144);

    auto tk_zz_xxxxzz = pbuffer.data(idx_kin_di + 145);

    auto tk_zz_xxxyyy = pbuffer.data(idx_kin_di + 146);

    auto tk_zz_xxxyyz = pbuffer.data(idx_kin_di + 147);

    auto tk_zz_xxxyzz = pbuffer.data(idx_kin_di + 148);

    auto tk_zz_xxxzzz = pbuffer.data(idx_kin_di + 149);

    auto tk_zz_xxyyyy = pbuffer.data(idx_kin_di + 150);

    auto tk_zz_xxyyyz = pbuffer.data(idx_kin_di + 151);

    auto tk_zz_xxyyzz = pbuffer.data(idx_kin_di + 152);

    auto tk_zz_xxyzzz = pbuffer.data(idx_kin_di + 153);

    auto tk_zz_xxzzzz = pbuffer.data(idx_kin_di + 154);

    auto tk_zz_xyyyyy = pbuffer.data(idx_kin_di + 155);

    auto tk_zz_xyyyyz = pbuffer.data(idx_kin_di + 156);

    auto tk_zz_xyyyzz = pbuffer.data(idx_kin_di + 157);

    auto tk_zz_xyyzzz = pbuffer.data(idx_kin_di + 158);

    auto tk_zz_xyzzzz = pbuffer.data(idx_kin_di + 159);

    auto tk_zz_xzzzzz = pbuffer.data(idx_kin_di + 160);

    auto tk_zz_yyyyyy = pbuffer.data(idx_kin_di + 161);

    auto tk_zz_yyyyyz = pbuffer.data(idx_kin_di + 162);

    auto tk_zz_yyyyzz = pbuffer.data(idx_kin_di + 163);

    auto tk_zz_yyyzzz = pbuffer.data(idx_kin_di + 164);

    auto tk_zz_yyzzzz = pbuffer.data(idx_kin_di + 165);

    auto tk_zz_yzzzzz = pbuffer.data(idx_kin_di + 166);

    auto tk_zz_zzzzzz = pbuffer.data(idx_kin_di + 167);

    // Set up components of auxiliary buffer : FI

    auto ts_xxx_xxxxxx = pbuffer.data(idx_ovl_fi);

    auto ts_xxx_xxxxxy = pbuffer.data(idx_ovl_fi + 1);

    auto ts_xxx_xxxxxz = pbuffer.data(idx_ovl_fi + 2);

    auto ts_xxx_xxxxyy = pbuffer.data(idx_ovl_fi + 3);

    auto ts_xxx_xxxxyz = pbuffer.data(idx_ovl_fi + 4);

    auto ts_xxx_xxxxzz = pbuffer.data(idx_ovl_fi + 5);

    auto ts_xxx_xxxyyy = pbuffer.data(idx_ovl_fi + 6);

    auto ts_xxx_xxxyyz = pbuffer.data(idx_ovl_fi + 7);

    auto ts_xxx_xxxyzz = pbuffer.data(idx_ovl_fi + 8);

    auto ts_xxx_xxxzzz = pbuffer.data(idx_ovl_fi + 9);

    auto ts_xxx_xxyyyy = pbuffer.data(idx_ovl_fi + 10);

    auto ts_xxx_xxyyyz = pbuffer.data(idx_ovl_fi + 11);

    auto ts_xxx_xxyyzz = pbuffer.data(idx_ovl_fi + 12);

    auto ts_xxx_xxyzzz = pbuffer.data(idx_ovl_fi + 13);

    auto ts_xxx_xxzzzz = pbuffer.data(idx_ovl_fi + 14);

    auto ts_xxx_xyyyyy = pbuffer.data(idx_ovl_fi + 15);

    auto ts_xxx_xyyyyz = pbuffer.data(idx_ovl_fi + 16);

    auto ts_xxx_xyyyzz = pbuffer.data(idx_ovl_fi + 17);

    auto ts_xxx_xyyzzz = pbuffer.data(idx_ovl_fi + 18);

    auto ts_xxx_xyzzzz = pbuffer.data(idx_ovl_fi + 19);

    auto ts_xxx_xzzzzz = pbuffer.data(idx_ovl_fi + 20);

    auto ts_xxx_yyyyyy = pbuffer.data(idx_ovl_fi + 21);

    auto ts_xxx_yyyyyz = pbuffer.data(idx_ovl_fi + 22);

    auto ts_xxx_yyyyzz = pbuffer.data(idx_ovl_fi + 23);

    auto ts_xxx_yyyzzz = pbuffer.data(idx_ovl_fi + 24);

    auto ts_xxx_yyzzzz = pbuffer.data(idx_ovl_fi + 25);

    auto ts_xxx_yzzzzz = pbuffer.data(idx_ovl_fi + 26);

    auto ts_xxx_zzzzzz = pbuffer.data(idx_ovl_fi + 27);

    auto ts_xxy_xxxxxx = pbuffer.data(idx_ovl_fi + 28);

    auto ts_xxy_xxxxxy = pbuffer.data(idx_ovl_fi + 29);

    auto ts_xxy_xxxxxz = pbuffer.data(idx_ovl_fi + 30);

    auto ts_xxy_xxxxyy = pbuffer.data(idx_ovl_fi + 31);

    auto ts_xxy_xxxxyz = pbuffer.data(idx_ovl_fi + 32);

    auto ts_xxy_xxxxzz = pbuffer.data(idx_ovl_fi + 33);

    auto ts_xxy_xxxyyy = pbuffer.data(idx_ovl_fi + 34);

    auto ts_xxy_xxxyyz = pbuffer.data(idx_ovl_fi + 35);

    auto ts_xxy_xxxyzz = pbuffer.data(idx_ovl_fi + 36);

    auto ts_xxy_xxxzzz = pbuffer.data(idx_ovl_fi + 37);

    auto ts_xxy_xxyyyy = pbuffer.data(idx_ovl_fi + 38);

    auto ts_xxy_xxyyyz = pbuffer.data(idx_ovl_fi + 39);

    auto ts_xxy_xxyyzz = pbuffer.data(idx_ovl_fi + 40);

    auto ts_xxy_xxyzzz = pbuffer.data(idx_ovl_fi + 41);

    auto ts_xxy_xxzzzz = pbuffer.data(idx_ovl_fi + 42);

    auto ts_xxy_xyyyyy = pbuffer.data(idx_ovl_fi + 43);

    auto ts_xxy_xyyyyz = pbuffer.data(idx_ovl_fi + 44);

    auto ts_xxy_xyyyzz = pbuffer.data(idx_ovl_fi + 45);

    auto ts_xxy_xyyzzz = pbuffer.data(idx_ovl_fi + 46);

    auto ts_xxy_xyzzzz = pbuffer.data(idx_ovl_fi + 47);

    auto ts_xxy_xzzzzz = pbuffer.data(idx_ovl_fi + 48);

    auto ts_xxy_yyyyyy = pbuffer.data(idx_ovl_fi + 49);

    auto ts_xxy_yyyyyz = pbuffer.data(idx_ovl_fi + 50);

    auto ts_xxy_yyyyzz = pbuffer.data(idx_ovl_fi + 51);

    auto ts_xxy_yyyzzz = pbuffer.data(idx_ovl_fi + 52);

    auto ts_xxy_yyzzzz = pbuffer.data(idx_ovl_fi + 53);

    auto ts_xxy_yzzzzz = pbuffer.data(idx_ovl_fi + 54);

    auto ts_xxy_zzzzzz = pbuffer.data(idx_ovl_fi + 55);

    auto ts_xxz_xxxxxx = pbuffer.data(idx_ovl_fi + 56);

    auto ts_xxz_xxxxxy = pbuffer.data(idx_ovl_fi + 57);

    auto ts_xxz_xxxxxz = pbuffer.data(idx_ovl_fi + 58);

    auto ts_xxz_xxxxyy = pbuffer.data(idx_ovl_fi + 59);

    auto ts_xxz_xxxxyz = pbuffer.data(idx_ovl_fi + 60);

    auto ts_xxz_xxxxzz = pbuffer.data(idx_ovl_fi + 61);

    auto ts_xxz_xxxyyy = pbuffer.data(idx_ovl_fi + 62);

    auto ts_xxz_xxxyyz = pbuffer.data(idx_ovl_fi + 63);

    auto ts_xxz_xxxyzz = pbuffer.data(idx_ovl_fi + 64);

    auto ts_xxz_xxxzzz = pbuffer.data(idx_ovl_fi + 65);

    auto ts_xxz_xxyyyy = pbuffer.data(idx_ovl_fi + 66);

    auto ts_xxz_xxyyyz = pbuffer.data(idx_ovl_fi + 67);

    auto ts_xxz_xxyyzz = pbuffer.data(idx_ovl_fi + 68);

    auto ts_xxz_xxyzzz = pbuffer.data(idx_ovl_fi + 69);

    auto ts_xxz_xxzzzz = pbuffer.data(idx_ovl_fi + 70);

    auto ts_xxz_xyyyyy = pbuffer.data(idx_ovl_fi + 71);

    auto ts_xxz_xyyyyz = pbuffer.data(idx_ovl_fi + 72);

    auto ts_xxz_xyyyzz = pbuffer.data(idx_ovl_fi + 73);

    auto ts_xxz_xyyzzz = pbuffer.data(idx_ovl_fi + 74);

    auto ts_xxz_xyzzzz = pbuffer.data(idx_ovl_fi + 75);

    auto ts_xxz_xzzzzz = pbuffer.data(idx_ovl_fi + 76);

    auto ts_xxz_yyyyyy = pbuffer.data(idx_ovl_fi + 77);

    auto ts_xxz_yyyyyz = pbuffer.data(idx_ovl_fi + 78);

    auto ts_xxz_yyyyzz = pbuffer.data(idx_ovl_fi + 79);

    auto ts_xxz_yyyzzz = pbuffer.data(idx_ovl_fi + 80);

    auto ts_xxz_yyzzzz = pbuffer.data(idx_ovl_fi + 81);

    auto ts_xxz_yzzzzz = pbuffer.data(idx_ovl_fi + 82);

    auto ts_xxz_zzzzzz = pbuffer.data(idx_ovl_fi + 83);

    auto ts_xyy_xxxxxx = pbuffer.data(idx_ovl_fi + 84);

    auto ts_xyy_xxxxxy = pbuffer.data(idx_ovl_fi + 85);

    auto ts_xyy_xxxxxz = pbuffer.data(idx_ovl_fi + 86);

    auto ts_xyy_xxxxyy = pbuffer.data(idx_ovl_fi + 87);

    auto ts_xyy_xxxxyz = pbuffer.data(idx_ovl_fi + 88);

    auto ts_xyy_xxxxzz = pbuffer.data(idx_ovl_fi + 89);

    auto ts_xyy_xxxyyy = pbuffer.data(idx_ovl_fi + 90);

    auto ts_xyy_xxxyyz = pbuffer.data(idx_ovl_fi + 91);

    auto ts_xyy_xxxyzz = pbuffer.data(idx_ovl_fi + 92);

    auto ts_xyy_xxxzzz = pbuffer.data(idx_ovl_fi + 93);

    auto ts_xyy_xxyyyy = pbuffer.data(idx_ovl_fi + 94);

    auto ts_xyy_xxyyyz = pbuffer.data(idx_ovl_fi + 95);

    auto ts_xyy_xxyyzz = pbuffer.data(idx_ovl_fi + 96);

    auto ts_xyy_xxyzzz = pbuffer.data(idx_ovl_fi + 97);

    auto ts_xyy_xxzzzz = pbuffer.data(idx_ovl_fi + 98);

    auto ts_xyy_xyyyyy = pbuffer.data(idx_ovl_fi + 99);

    auto ts_xyy_xyyyyz = pbuffer.data(idx_ovl_fi + 100);

    auto ts_xyy_xyyyzz = pbuffer.data(idx_ovl_fi + 101);

    auto ts_xyy_xyyzzz = pbuffer.data(idx_ovl_fi + 102);

    auto ts_xyy_xyzzzz = pbuffer.data(idx_ovl_fi + 103);

    auto ts_xyy_xzzzzz = pbuffer.data(idx_ovl_fi + 104);

    auto ts_xyy_yyyyyy = pbuffer.data(idx_ovl_fi + 105);

    auto ts_xyy_yyyyyz = pbuffer.data(idx_ovl_fi + 106);

    auto ts_xyy_yyyyzz = pbuffer.data(idx_ovl_fi + 107);

    auto ts_xyy_yyyzzz = pbuffer.data(idx_ovl_fi + 108);

    auto ts_xyy_yyzzzz = pbuffer.data(idx_ovl_fi + 109);

    auto ts_xyy_yzzzzz = pbuffer.data(idx_ovl_fi + 110);

    auto ts_xyy_zzzzzz = pbuffer.data(idx_ovl_fi + 111);

    auto ts_xyz_xxxxxx = pbuffer.data(idx_ovl_fi + 112);

    auto ts_xyz_xxxxxy = pbuffer.data(idx_ovl_fi + 113);

    auto ts_xyz_xxxxxz = pbuffer.data(idx_ovl_fi + 114);

    auto ts_xyz_xxxxyy = pbuffer.data(idx_ovl_fi + 115);

    auto ts_xyz_xxxxyz = pbuffer.data(idx_ovl_fi + 116);

    auto ts_xyz_xxxxzz = pbuffer.data(idx_ovl_fi + 117);

    auto ts_xyz_xxxyyy = pbuffer.data(idx_ovl_fi + 118);

    auto ts_xyz_xxxyyz = pbuffer.data(idx_ovl_fi + 119);

    auto ts_xyz_xxxyzz = pbuffer.data(idx_ovl_fi + 120);

    auto ts_xyz_xxxzzz = pbuffer.data(idx_ovl_fi + 121);

    auto ts_xyz_xxyyyy = pbuffer.data(idx_ovl_fi + 122);

    auto ts_xyz_xxyyyz = pbuffer.data(idx_ovl_fi + 123);

    auto ts_xyz_xxyyzz = pbuffer.data(idx_ovl_fi + 124);

    auto ts_xyz_xxyzzz = pbuffer.data(idx_ovl_fi + 125);

    auto ts_xyz_xxzzzz = pbuffer.data(idx_ovl_fi + 126);

    auto ts_xyz_xyyyyy = pbuffer.data(idx_ovl_fi + 127);

    auto ts_xyz_xyyyyz = pbuffer.data(idx_ovl_fi + 128);

    auto ts_xyz_xyyyzz = pbuffer.data(idx_ovl_fi + 129);

    auto ts_xyz_xyyzzz = pbuffer.data(idx_ovl_fi + 130);

    auto ts_xyz_xyzzzz = pbuffer.data(idx_ovl_fi + 131);

    auto ts_xyz_xzzzzz = pbuffer.data(idx_ovl_fi + 132);

    auto ts_xyz_yyyyyy = pbuffer.data(idx_ovl_fi + 133);

    auto ts_xyz_yyyyyz = pbuffer.data(idx_ovl_fi + 134);

    auto ts_xyz_yyyyzz = pbuffer.data(idx_ovl_fi + 135);

    auto ts_xyz_yyyzzz = pbuffer.data(idx_ovl_fi + 136);

    auto ts_xyz_yyzzzz = pbuffer.data(idx_ovl_fi + 137);

    auto ts_xyz_yzzzzz = pbuffer.data(idx_ovl_fi + 138);

    auto ts_xyz_zzzzzz = pbuffer.data(idx_ovl_fi + 139);

    auto ts_xzz_xxxxxx = pbuffer.data(idx_ovl_fi + 140);

    auto ts_xzz_xxxxxy = pbuffer.data(idx_ovl_fi + 141);

    auto ts_xzz_xxxxxz = pbuffer.data(idx_ovl_fi + 142);

    auto ts_xzz_xxxxyy = pbuffer.data(idx_ovl_fi + 143);

    auto ts_xzz_xxxxyz = pbuffer.data(idx_ovl_fi + 144);

    auto ts_xzz_xxxxzz = pbuffer.data(idx_ovl_fi + 145);

    auto ts_xzz_xxxyyy = pbuffer.data(idx_ovl_fi + 146);

    auto ts_xzz_xxxyyz = pbuffer.data(idx_ovl_fi + 147);

    auto ts_xzz_xxxyzz = pbuffer.data(idx_ovl_fi + 148);

    auto ts_xzz_xxxzzz = pbuffer.data(idx_ovl_fi + 149);

    auto ts_xzz_xxyyyy = pbuffer.data(idx_ovl_fi + 150);

    auto ts_xzz_xxyyyz = pbuffer.data(idx_ovl_fi + 151);

    auto ts_xzz_xxyyzz = pbuffer.data(idx_ovl_fi + 152);

    auto ts_xzz_xxyzzz = pbuffer.data(idx_ovl_fi + 153);

    auto ts_xzz_xxzzzz = pbuffer.data(idx_ovl_fi + 154);

    auto ts_xzz_xyyyyy = pbuffer.data(idx_ovl_fi + 155);

    auto ts_xzz_xyyyyz = pbuffer.data(idx_ovl_fi + 156);

    auto ts_xzz_xyyyzz = pbuffer.data(idx_ovl_fi + 157);

    auto ts_xzz_xyyzzz = pbuffer.data(idx_ovl_fi + 158);

    auto ts_xzz_xyzzzz = pbuffer.data(idx_ovl_fi + 159);

    auto ts_xzz_xzzzzz = pbuffer.data(idx_ovl_fi + 160);

    auto ts_xzz_yyyyyy = pbuffer.data(idx_ovl_fi + 161);

    auto ts_xzz_yyyyyz = pbuffer.data(idx_ovl_fi + 162);

    auto ts_xzz_yyyyzz = pbuffer.data(idx_ovl_fi + 163);

    auto ts_xzz_yyyzzz = pbuffer.data(idx_ovl_fi + 164);

    auto ts_xzz_yyzzzz = pbuffer.data(idx_ovl_fi + 165);

    auto ts_xzz_yzzzzz = pbuffer.data(idx_ovl_fi + 166);

    auto ts_xzz_zzzzzz = pbuffer.data(idx_ovl_fi + 167);

    auto ts_yyy_xxxxxx = pbuffer.data(idx_ovl_fi + 168);

    auto ts_yyy_xxxxxy = pbuffer.data(idx_ovl_fi + 169);

    auto ts_yyy_xxxxxz = pbuffer.data(idx_ovl_fi + 170);

    auto ts_yyy_xxxxyy = pbuffer.data(idx_ovl_fi + 171);

    auto ts_yyy_xxxxyz = pbuffer.data(idx_ovl_fi + 172);

    auto ts_yyy_xxxxzz = pbuffer.data(idx_ovl_fi + 173);

    auto ts_yyy_xxxyyy = pbuffer.data(idx_ovl_fi + 174);

    auto ts_yyy_xxxyyz = pbuffer.data(idx_ovl_fi + 175);

    auto ts_yyy_xxxyzz = pbuffer.data(idx_ovl_fi + 176);

    auto ts_yyy_xxxzzz = pbuffer.data(idx_ovl_fi + 177);

    auto ts_yyy_xxyyyy = pbuffer.data(idx_ovl_fi + 178);

    auto ts_yyy_xxyyyz = pbuffer.data(idx_ovl_fi + 179);

    auto ts_yyy_xxyyzz = pbuffer.data(idx_ovl_fi + 180);

    auto ts_yyy_xxyzzz = pbuffer.data(idx_ovl_fi + 181);

    auto ts_yyy_xxzzzz = pbuffer.data(idx_ovl_fi + 182);

    auto ts_yyy_xyyyyy = pbuffer.data(idx_ovl_fi + 183);

    auto ts_yyy_xyyyyz = pbuffer.data(idx_ovl_fi + 184);

    auto ts_yyy_xyyyzz = pbuffer.data(idx_ovl_fi + 185);

    auto ts_yyy_xyyzzz = pbuffer.data(idx_ovl_fi + 186);

    auto ts_yyy_xyzzzz = pbuffer.data(idx_ovl_fi + 187);

    auto ts_yyy_xzzzzz = pbuffer.data(idx_ovl_fi + 188);

    auto ts_yyy_yyyyyy = pbuffer.data(idx_ovl_fi + 189);

    auto ts_yyy_yyyyyz = pbuffer.data(idx_ovl_fi + 190);

    auto ts_yyy_yyyyzz = pbuffer.data(idx_ovl_fi + 191);

    auto ts_yyy_yyyzzz = pbuffer.data(idx_ovl_fi + 192);

    auto ts_yyy_yyzzzz = pbuffer.data(idx_ovl_fi + 193);

    auto ts_yyy_yzzzzz = pbuffer.data(idx_ovl_fi + 194);

    auto ts_yyy_zzzzzz = pbuffer.data(idx_ovl_fi + 195);

    auto ts_yyz_xxxxxx = pbuffer.data(idx_ovl_fi + 196);

    auto ts_yyz_xxxxxy = pbuffer.data(idx_ovl_fi + 197);

    auto ts_yyz_xxxxxz = pbuffer.data(idx_ovl_fi + 198);

    auto ts_yyz_xxxxyy = pbuffer.data(idx_ovl_fi + 199);

    auto ts_yyz_xxxxyz = pbuffer.data(idx_ovl_fi + 200);

    auto ts_yyz_xxxxzz = pbuffer.data(idx_ovl_fi + 201);

    auto ts_yyz_xxxyyy = pbuffer.data(idx_ovl_fi + 202);

    auto ts_yyz_xxxyyz = pbuffer.data(idx_ovl_fi + 203);

    auto ts_yyz_xxxyzz = pbuffer.data(idx_ovl_fi + 204);

    auto ts_yyz_xxxzzz = pbuffer.data(idx_ovl_fi + 205);

    auto ts_yyz_xxyyyy = pbuffer.data(idx_ovl_fi + 206);

    auto ts_yyz_xxyyyz = pbuffer.data(idx_ovl_fi + 207);

    auto ts_yyz_xxyyzz = pbuffer.data(idx_ovl_fi + 208);

    auto ts_yyz_xxyzzz = pbuffer.data(idx_ovl_fi + 209);

    auto ts_yyz_xxzzzz = pbuffer.data(idx_ovl_fi + 210);

    auto ts_yyz_xyyyyy = pbuffer.data(idx_ovl_fi + 211);

    auto ts_yyz_xyyyyz = pbuffer.data(idx_ovl_fi + 212);

    auto ts_yyz_xyyyzz = pbuffer.data(idx_ovl_fi + 213);

    auto ts_yyz_xyyzzz = pbuffer.data(idx_ovl_fi + 214);

    auto ts_yyz_xyzzzz = pbuffer.data(idx_ovl_fi + 215);

    auto ts_yyz_xzzzzz = pbuffer.data(idx_ovl_fi + 216);

    auto ts_yyz_yyyyyy = pbuffer.data(idx_ovl_fi + 217);

    auto ts_yyz_yyyyyz = pbuffer.data(idx_ovl_fi + 218);

    auto ts_yyz_yyyyzz = pbuffer.data(idx_ovl_fi + 219);

    auto ts_yyz_yyyzzz = pbuffer.data(idx_ovl_fi + 220);

    auto ts_yyz_yyzzzz = pbuffer.data(idx_ovl_fi + 221);

    auto ts_yyz_yzzzzz = pbuffer.data(idx_ovl_fi + 222);

    auto ts_yyz_zzzzzz = pbuffer.data(idx_ovl_fi + 223);

    auto ts_yzz_xxxxxx = pbuffer.data(idx_ovl_fi + 224);

    auto ts_yzz_xxxxxy = pbuffer.data(idx_ovl_fi + 225);

    auto ts_yzz_xxxxxz = pbuffer.data(idx_ovl_fi + 226);

    auto ts_yzz_xxxxyy = pbuffer.data(idx_ovl_fi + 227);

    auto ts_yzz_xxxxyz = pbuffer.data(idx_ovl_fi + 228);

    auto ts_yzz_xxxxzz = pbuffer.data(idx_ovl_fi + 229);

    auto ts_yzz_xxxyyy = pbuffer.data(idx_ovl_fi + 230);

    auto ts_yzz_xxxyyz = pbuffer.data(idx_ovl_fi + 231);

    auto ts_yzz_xxxyzz = pbuffer.data(idx_ovl_fi + 232);

    auto ts_yzz_xxxzzz = pbuffer.data(idx_ovl_fi + 233);

    auto ts_yzz_xxyyyy = pbuffer.data(idx_ovl_fi + 234);

    auto ts_yzz_xxyyyz = pbuffer.data(idx_ovl_fi + 235);

    auto ts_yzz_xxyyzz = pbuffer.data(idx_ovl_fi + 236);

    auto ts_yzz_xxyzzz = pbuffer.data(idx_ovl_fi + 237);

    auto ts_yzz_xxzzzz = pbuffer.data(idx_ovl_fi + 238);

    auto ts_yzz_xyyyyy = pbuffer.data(idx_ovl_fi + 239);

    auto ts_yzz_xyyyyz = pbuffer.data(idx_ovl_fi + 240);

    auto ts_yzz_xyyyzz = pbuffer.data(idx_ovl_fi + 241);

    auto ts_yzz_xyyzzz = pbuffer.data(idx_ovl_fi + 242);

    auto ts_yzz_xyzzzz = pbuffer.data(idx_ovl_fi + 243);

    auto ts_yzz_xzzzzz = pbuffer.data(idx_ovl_fi + 244);

    auto ts_yzz_yyyyyy = pbuffer.data(idx_ovl_fi + 245);

    auto ts_yzz_yyyyyz = pbuffer.data(idx_ovl_fi + 246);

    auto ts_yzz_yyyyzz = pbuffer.data(idx_ovl_fi + 247);

    auto ts_yzz_yyyzzz = pbuffer.data(idx_ovl_fi + 248);

    auto ts_yzz_yyzzzz = pbuffer.data(idx_ovl_fi + 249);

    auto ts_yzz_yzzzzz = pbuffer.data(idx_ovl_fi + 250);

    auto ts_yzz_zzzzzz = pbuffer.data(idx_ovl_fi + 251);

    auto ts_zzz_xxxxxx = pbuffer.data(idx_ovl_fi + 252);

    auto ts_zzz_xxxxxy = pbuffer.data(idx_ovl_fi + 253);

    auto ts_zzz_xxxxxz = pbuffer.data(idx_ovl_fi + 254);

    auto ts_zzz_xxxxyy = pbuffer.data(idx_ovl_fi + 255);

    auto ts_zzz_xxxxyz = pbuffer.data(idx_ovl_fi + 256);

    auto ts_zzz_xxxxzz = pbuffer.data(idx_ovl_fi + 257);

    auto ts_zzz_xxxyyy = pbuffer.data(idx_ovl_fi + 258);

    auto ts_zzz_xxxyyz = pbuffer.data(idx_ovl_fi + 259);

    auto ts_zzz_xxxyzz = pbuffer.data(idx_ovl_fi + 260);

    auto ts_zzz_xxxzzz = pbuffer.data(idx_ovl_fi + 261);

    auto ts_zzz_xxyyyy = pbuffer.data(idx_ovl_fi + 262);

    auto ts_zzz_xxyyyz = pbuffer.data(idx_ovl_fi + 263);

    auto ts_zzz_xxyyzz = pbuffer.data(idx_ovl_fi + 264);

    auto ts_zzz_xxyzzz = pbuffer.data(idx_ovl_fi + 265);

    auto ts_zzz_xxzzzz = pbuffer.data(idx_ovl_fi + 266);

    auto ts_zzz_xyyyyy = pbuffer.data(idx_ovl_fi + 267);

    auto ts_zzz_xyyyyz = pbuffer.data(idx_ovl_fi + 268);

    auto ts_zzz_xyyyzz = pbuffer.data(idx_ovl_fi + 269);

    auto ts_zzz_xyyzzz = pbuffer.data(idx_ovl_fi + 270);

    auto ts_zzz_xyzzzz = pbuffer.data(idx_ovl_fi + 271);

    auto ts_zzz_xzzzzz = pbuffer.data(idx_ovl_fi + 272);

    auto ts_zzz_yyyyyy = pbuffer.data(idx_ovl_fi + 273);

    auto ts_zzz_yyyyyz = pbuffer.data(idx_ovl_fi + 274);

    auto ts_zzz_yyyyzz = pbuffer.data(idx_ovl_fi + 275);

    auto ts_zzz_yyyzzz = pbuffer.data(idx_ovl_fi + 276);

    auto ts_zzz_yyzzzz = pbuffer.data(idx_ovl_fi + 277);

    auto ts_zzz_yzzzzz = pbuffer.data(idx_ovl_fi + 278);

    auto ts_zzz_zzzzzz = pbuffer.data(idx_ovl_fi + 279);

    // Set up 0-28 components of targeted buffer : FI

    auto tk_xxx_xxxxxx = pbuffer.data(idx_kin_fi);

    auto tk_xxx_xxxxxy = pbuffer.data(idx_kin_fi + 1);

    auto tk_xxx_xxxxxz = pbuffer.data(idx_kin_fi + 2);

    auto tk_xxx_xxxxyy = pbuffer.data(idx_kin_fi + 3);

    auto tk_xxx_xxxxyz = pbuffer.data(idx_kin_fi + 4);

    auto tk_xxx_xxxxzz = pbuffer.data(idx_kin_fi + 5);

    auto tk_xxx_xxxyyy = pbuffer.data(idx_kin_fi + 6);

    auto tk_xxx_xxxyyz = pbuffer.data(idx_kin_fi + 7);

    auto tk_xxx_xxxyzz = pbuffer.data(idx_kin_fi + 8);

    auto tk_xxx_xxxzzz = pbuffer.data(idx_kin_fi + 9);

    auto tk_xxx_xxyyyy = pbuffer.data(idx_kin_fi + 10);

    auto tk_xxx_xxyyyz = pbuffer.data(idx_kin_fi + 11);

    auto tk_xxx_xxyyzz = pbuffer.data(idx_kin_fi + 12);

    auto tk_xxx_xxyzzz = pbuffer.data(idx_kin_fi + 13);

    auto tk_xxx_xxzzzz = pbuffer.data(idx_kin_fi + 14);

    auto tk_xxx_xyyyyy = pbuffer.data(idx_kin_fi + 15);

    auto tk_xxx_xyyyyz = pbuffer.data(idx_kin_fi + 16);

    auto tk_xxx_xyyyzz = pbuffer.data(idx_kin_fi + 17);

    auto tk_xxx_xyyzzz = pbuffer.data(idx_kin_fi + 18);

    auto tk_xxx_xyzzzz = pbuffer.data(idx_kin_fi + 19);

    auto tk_xxx_xzzzzz = pbuffer.data(idx_kin_fi + 20);

    auto tk_xxx_yyyyyy = pbuffer.data(idx_kin_fi + 21);

    auto tk_xxx_yyyyyz = pbuffer.data(idx_kin_fi + 22);

    auto tk_xxx_yyyyzz = pbuffer.data(idx_kin_fi + 23);

    auto tk_xxx_yyyzzz = pbuffer.data(idx_kin_fi + 24);

    auto tk_xxx_yyzzzz = pbuffer.data(idx_kin_fi + 25);

    auto tk_xxx_yzzzzz = pbuffer.data(idx_kin_fi + 26);

    auto tk_xxx_zzzzzz = pbuffer.data(idx_kin_fi + 27);

#pragma omp simd aligned(pa_x,              \
                             tk_x_xxxxxx,   \
                             tk_x_xxxxxy,   \
                             tk_x_xxxxxz,   \
                             tk_x_xxxxyy,   \
                             tk_x_xxxxyz,   \
                             tk_x_xxxxzz,   \
                             tk_x_xxxyyy,   \
                             tk_x_xxxyyz,   \
                             tk_x_xxxyzz,   \
                             tk_x_xxxzzz,   \
                             tk_x_xxyyyy,   \
                             tk_x_xxyyyz,   \
                             tk_x_xxyyzz,   \
                             tk_x_xxyzzz,   \
                             tk_x_xxzzzz,   \
                             tk_x_xyyyyy,   \
                             tk_x_xyyyyz,   \
                             tk_x_xyyyzz,   \
                             tk_x_xyyzzz,   \
                             tk_x_xyzzzz,   \
                             tk_x_xzzzzz,   \
                             tk_x_yyyyyy,   \
                             tk_x_yyyyyz,   \
                             tk_x_yyyyzz,   \
                             tk_x_yyyzzz,   \
                             tk_x_yyzzzz,   \
                             tk_x_yzzzzz,   \
                             tk_x_zzzzzz,   \
                             tk_xx_xxxxx,   \
                             tk_xx_xxxxxx,  \
                             tk_xx_xxxxxy,  \
                             tk_xx_xxxxxz,  \
                             tk_xx_xxxxy,   \
                             tk_xx_xxxxyy,  \
                             tk_xx_xxxxyz,  \
                             tk_xx_xxxxz,   \
                             tk_xx_xxxxzz,  \
                             tk_xx_xxxyy,   \
                             tk_xx_xxxyyy,  \
                             tk_xx_xxxyyz,  \
                             tk_xx_xxxyz,   \
                             tk_xx_xxxyzz,  \
                             tk_xx_xxxzz,   \
                             tk_xx_xxxzzz,  \
                             tk_xx_xxyyy,   \
                             tk_xx_xxyyyy,  \
                             tk_xx_xxyyyz,  \
                             tk_xx_xxyyz,   \
                             tk_xx_xxyyzz,  \
                             tk_xx_xxyzz,   \
                             tk_xx_xxyzzz,  \
                             tk_xx_xxzzz,   \
                             tk_xx_xxzzzz,  \
                             tk_xx_xyyyy,   \
                             tk_xx_xyyyyy,  \
                             tk_xx_xyyyyz,  \
                             tk_xx_xyyyz,   \
                             tk_xx_xyyyzz,  \
                             tk_xx_xyyzz,   \
                             tk_xx_xyyzzz,  \
                             tk_xx_xyzzz,   \
                             tk_xx_xyzzzz,  \
                             tk_xx_xzzzz,   \
                             tk_xx_xzzzzz,  \
                             tk_xx_yyyyy,   \
                             tk_xx_yyyyyy,  \
                             tk_xx_yyyyyz,  \
                             tk_xx_yyyyz,   \
                             tk_xx_yyyyzz,  \
                             tk_xx_yyyzz,   \
                             tk_xx_yyyzzz,  \
                             tk_xx_yyzzz,   \
                             tk_xx_yyzzzz,  \
                             tk_xx_yzzzz,   \
                             tk_xx_yzzzzz,  \
                             tk_xx_zzzzz,   \
                             tk_xx_zzzzzz,  \
                             tk_xxx_xxxxxx, \
                             tk_xxx_xxxxxy, \
                             tk_xxx_xxxxxz, \
                             tk_xxx_xxxxyy, \
                             tk_xxx_xxxxyz, \
                             tk_xxx_xxxxzz, \
                             tk_xxx_xxxyyy, \
                             tk_xxx_xxxyyz, \
                             tk_xxx_xxxyzz, \
                             tk_xxx_xxxzzz, \
                             tk_xxx_xxyyyy, \
                             tk_xxx_xxyyyz, \
                             tk_xxx_xxyyzz, \
                             tk_xxx_xxyzzz, \
                             tk_xxx_xxzzzz, \
                             tk_xxx_xyyyyy, \
                             tk_xxx_xyyyyz, \
                             tk_xxx_xyyyzz, \
                             tk_xxx_xyyzzz, \
                             tk_xxx_xyzzzz, \
                             tk_xxx_xzzzzz, \
                             tk_xxx_yyyyyy, \
                             tk_xxx_yyyyyz, \
                             tk_xxx_yyyyzz, \
                             tk_xxx_yyyzzz, \
                             tk_xxx_yyzzzz, \
                             tk_xxx_yzzzzz, \
                             tk_xxx_zzzzzz, \
                             ts_x_xxxxxx,   \
                             ts_x_xxxxxy,   \
                             ts_x_xxxxxz,   \
                             ts_x_xxxxyy,   \
                             ts_x_xxxxyz,   \
                             ts_x_xxxxzz,   \
                             ts_x_xxxyyy,   \
                             ts_x_xxxyyz,   \
                             ts_x_xxxyzz,   \
                             ts_x_xxxzzz,   \
                             ts_x_xxyyyy,   \
                             ts_x_xxyyyz,   \
                             ts_x_xxyyzz,   \
                             ts_x_xxyzzz,   \
                             ts_x_xxzzzz,   \
                             ts_x_xyyyyy,   \
                             ts_x_xyyyyz,   \
                             ts_x_xyyyzz,   \
                             ts_x_xyyzzz,   \
                             ts_x_xyzzzz,   \
                             ts_x_xzzzzz,   \
                             ts_x_yyyyyy,   \
                             ts_x_yyyyyz,   \
                             ts_x_yyyyzz,   \
                             ts_x_yyyzzz,   \
                             ts_x_yyzzzz,   \
                             ts_x_yzzzzz,   \
                             ts_x_zzzzzz,   \
                             ts_xxx_xxxxxx, \
                             ts_xxx_xxxxxy, \
                             ts_xxx_xxxxxz, \
                             ts_xxx_xxxxyy, \
                             ts_xxx_xxxxyz, \
                             ts_xxx_xxxxzz, \
                             ts_xxx_xxxyyy, \
                             ts_xxx_xxxyyz, \
                             ts_xxx_xxxyzz, \
                             ts_xxx_xxxzzz, \
                             ts_xxx_xxyyyy, \
                             ts_xxx_xxyyyz, \
                             ts_xxx_xxyyzz, \
                             ts_xxx_xxyzzz, \
                             ts_xxx_xxzzzz, \
                             ts_xxx_xyyyyy, \
                             ts_xxx_xyyyyz, \
                             ts_xxx_xyyyzz, \
                             ts_xxx_xyyzzz, \
                             ts_xxx_xyzzzz, \
                             ts_xxx_xzzzzz, \
                             ts_xxx_yyyyyy, \
                             ts_xxx_yyyyyz, \
                             ts_xxx_yyyyzz, \
                             ts_xxx_yyyzzz, \
                             ts_xxx_yyzzzz, \
                             ts_xxx_yzzzzz, \
                             ts_xxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxx_xxxxxx[i] = -4.0 * ts_x_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxxx[i] * fe_0 + 6.0 * tk_xx_xxxxx[i] * fe_0 +
                           tk_xx_xxxxxx[i] * pa_x[i] + 2.0 * ts_xxx_xxxxxx[i] * fz_0;

        tk_xxx_xxxxxy[i] = -4.0 * ts_x_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxxy[i] * fe_0 + 5.0 * tk_xx_xxxxy[i] * fe_0 +
                           tk_xx_xxxxxy[i] * pa_x[i] + 2.0 * ts_xxx_xxxxxy[i] * fz_0;

        tk_xxx_xxxxxz[i] = -4.0 * ts_x_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxxz[i] * fe_0 + 5.0 * tk_xx_xxxxz[i] * fe_0 +
                           tk_xx_xxxxxz[i] * pa_x[i] + 2.0 * ts_xxx_xxxxxz[i] * fz_0;

        tk_xxx_xxxxyy[i] = -4.0 * ts_x_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxyy[i] * fe_0 + 4.0 * tk_xx_xxxyy[i] * fe_0 +
                           tk_xx_xxxxyy[i] * pa_x[i] + 2.0 * ts_xxx_xxxxyy[i] * fz_0;

        tk_xxx_xxxxyz[i] = -4.0 * ts_x_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxyz[i] * fe_0 + 4.0 * tk_xx_xxxyz[i] * fe_0 +
                           tk_xx_xxxxyz[i] * pa_x[i] + 2.0 * ts_xxx_xxxxyz[i] * fz_0;

        tk_xxx_xxxxzz[i] = -4.0 * ts_x_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxxzz[i] * fe_0 + 4.0 * tk_xx_xxxzz[i] * fe_0 +
                           tk_xx_xxxxzz[i] * pa_x[i] + 2.0 * ts_xxx_xxxxzz[i] * fz_0;

        tk_xxx_xxxyyy[i] = -4.0 * ts_x_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxyyy[i] * fe_0 + 3.0 * tk_xx_xxyyy[i] * fe_0 +
                           tk_xx_xxxyyy[i] * pa_x[i] + 2.0 * ts_xxx_xxxyyy[i] * fz_0;

        tk_xxx_xxxyyz[i] = -4.0 * ts_x_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxyyz[i] * fe_0 + 3.0 * tk_xx_xxyyz[i] * fe_0 +
                           tk_xx_xxxyyz[i] * pa_x[i] + 2.0 * ts_xxx_xxxyyz[i] * fz_0;

        tk_xxx_xxxyzz[i] = -4.0 * ts_x_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxyzz[i] * fe_0 + 3.0 * tk_xx_xxyzz[i] * fe_0 +
                           tk_xx_xxxyzz[i] * pa_x[i] + 2.0 * ts_xxx_xxxyzz[i] * fz_0;

        tk_xxx_xxxzzz[i] = -4.0 * ts_x_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxxzzz[i] * fe_0 + 3.0 * tk_xx_xxzzz[i] * fe_0 +
                           tk_xx_xxxzzz[i] * pa_x[i] + 2.0 * ts_xxx_xxxzzz[i] * fz_0;

        tk_xxx_xxyyyy[i] = -4.0 * ts_x_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyyyy[i] * fe_0 + 2.0 * tk_xx_xyyyy[i] * fe_0 +
                           tk_xx_xxyyyy[i] * pa_x[i] + 2.0 * ts_xxx_xxyyyy[i] * fz_0;

        tk_xxx_xxyyyz[i] = -4.0 * ts_x_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyyyz[i] * fe_0 + 2.0 * tk_xx_xyyyz[i] * fe_0 +
                           tk_xx_xxyyyz[i] * pa_x[i] + 2.0 * ts_xxx_xxyyyz[i] * fz_0;

        tk_xxx_xxyyzz[i] = -4.0 * ts_x_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyyzz[i] * fe_0 + 2.0 * tk_xx_xyyzz[i] * fe_0 +
                           tk_xx_xxyyzz[i] * pa_x[i] + 2.0 * ts_xxx_xxyyzz[i] * fz_0;

        tk_xxx_xxyzzz[i] = -4.0 * ts_x_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxyzzz[i] * fe_0 + 2.0 * tk_xx_xyzzz[i] * fe_0 +
                           tk_xx_xxyzzz[i] * pa_x[i] + 2.0 * ts_xxx_xxyzzz[i] * fz_0;

        tk_xxx_xxzzzz[i] = -4.0 * ts_x_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xxzzzz[i] * fe_0 + 2.0 * tk_xx_xzzzz[i] * fe_0 +
                           tk_xx_xxzzzz[i] * pa_x[i] + 2.0 * ts_xxx_xxzzzz[i] * fz_0;

        tk_xxx_xyyyyy[i] = -4.0 * ts_x_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyyyy[i] * fe_0 + tk_xx_yyyyy[i] * fe_0 + tk_xx_xyyyyy[i] * pa_x[i] +
                           2.0 * ts_xxx_xyyyyy[i] * fz_0;

        tk_xxx_xyyyyz[i] = -4.0 * ts_x_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyyyz[i] * fe_0 + tk_xx_yyyyz[i] * fe_0 + tk_xx_xyyyyz[i] * pa_x[i] +
                           2.0 * ts_xxx_xyyyyz[i] * fz_0;

        tk_xxx_xyyyzz[i] = -4.0 * ts_x_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyyzz[i] * fe_0 + tk_xx_yyyzz[i] * fe_0 + tk_xx_xyyyzz[i] * pa_x[i] +
                           2.0 * ts_xxx_xyyyzz[i] * fz_0;

        tk_xxx_xyyzzz[i] = -4.0 * ts_x_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyyzzz[i] * fe_0 + tk_xx_yyzzz[i] * fe_0 + tk_xx_xyyzzz[i] * pa_x[i] +
                           2.0 * ts_xxx_xyyzzz[i] * fz_0;

        tk_xxx_xyzzzz[i] = -4.0 * ts_x_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xyzzzz[i] * fe_0 + tk_xx_yzzzz[i] * fe_0 + tk_xx_xyzzzz[i] * pa_x[i] +
                           2.0 * ts_xxx_xyzzzz[i] * fz_0;

        tk_xxx_xzzzzz[i] = -4.0 * ts_x_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_xzzzzz[i] * fe_0 + tk_xx_zzzzz[i] * fe_0 + tk_xx_xzzzzz[i] * pa_x[i] +
                           2.0 * ts_xxx_xzzzzz[i] * fz_0;

        tk_xxx_yyyyyy[i] =
            -4.0 * ts_x_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyyyy[i] * fe_0 + tk_xx_yyyyyy[i] * pa_x[i] + 2.0 * ts_xxx_yyyyyy[i] * fz_0;

        tk_xxx_yyyyyz[i] =
            -4.0 * ts_x_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyyyz[i] * fe_0 + tk_xx_yyyyyz[i] * pa_x[i] + 2.0 * ts_xxx_yyyyyz[i] * fz_0;

        tk_xxx_yyyyzz[i] =
            -4.0 * ts_x_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyyzz[i] * fe_0 + tk_xx_yyyyzz[i] * pa_x[i] + 2.0 * ts_xxx_yyyyzz[i] * fz_0;

        tk_xxx_yyyzzz[i] =
            -4.0 * ts_x_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyyzzz[i] * fe_0 + tk_xx_yyyzzz[i] * pa_x[i] + 2.0 * ts_xxx_yyyzzz[i] * fz_0;

        tk_xxx_yyzzzz[i] =
            -4.0 * ts_x_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yyzzzz[i] * fe_0 + tk_xx_yyzzzz[i] * pa_x[i] + 2.0 * ts_xxx_yyzzzz[i] * fz_0;

        tk_xxx_yzzzzz[i] =
            -4.0 * ts_x_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_yzzzzz[i] * fe_0 + tk_xx_yzzzzz[i] * pa_x[i] + 2.0 * ts_xxx_yzzzzz[i] * fz_0;

        tk_xxx_zzzzzz[i] =
            -4.0 * ts_x_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_x_zzzzzz[i] * fe_0 + tk_xx_zzzzzz[i] * pa_x[i] + 2.0 * ts_xxx_zzzzzz[i] * fz_0;
    }

    // Set up 28-56 components of targeted buffer : FI

    auto tk_xxy_xxxxxx = pbuffer.data(idx_kin_fi + 28);

    auto tk_xxy_xxxxxy = pbuffer.data(idx_kin_fi + 29);

    auto tk_xxy_xxxxxz = pbuffer.data(idx_kin_fi + 30);

    auto tk_xxy_xxxxyy = pbuffer.data(idx_kin_fi + 31);

    auto tk_xxy_xxxxyz = pbuffer.data(idx_kin_fi + 32);

    auto tk_xxy_xxxxzz = pbuffer.data(idx_kin_fi + 33);

    auto tk_xxy_xxxyyy = pbuffer.data(idx_kin_fi + 34);

    auto tk_xxy_xxxyyz = pbuffer.data(idx_kin_fi + 35);

    auto tk_xxy_xxxyzz = pbuffer.data(idx_kin_fi + 36);

    auto tk_xxy_xxxzzz = pbuffer.data(idx_kin_fi + 37);

    auto tk_xxy_xxyyyy = pbuffer.data(idx_kin_fi + 38);

    auto tk_xxy_xxyyyz = pbuffer.data(idx_kin_fi + 39);

    auto tk_xxy_xxyyzz = pbuffer.data(idx_kin_fi + 40);

    auto tk_xxy_xxyzzz = pbuffer.data(idx_kin_fi + 41);

    auto tk_xxy_xxzzzz = pbuffer.data(idx_kin_fi + 42);

    auto tk_xxy_xyyyyy = pbuffer.data(idx_kin_fi + 43);

    auto tk_xxy_xyyyyz = pbuffer.data(idx_kin_fi + 44);

    auto tk_xxy_xyyyzz = pbuffer.data(idx_kin_fi + 45);

    auto tk_xxy_xyyzzz = pbuffer.data(idx_kin_fi + 46);

    auto tk_xxy_xyzzzz = pbuffer.data(idx_kin_fi + 47);

    auto tk_xxy_xzzzzz = pbuffer.data(idx_kin_fi + 48);

    auto tk_xxy_yyyyyy = pbuffer.data(idx_kin_fi + 49);

    auto tk_xxy_yyyyyz = pbuffer.data(idx_kin_fi + 50);

    auto tk_xxy_yyyyzz = pbuffer.data(idx_kin_fi + 51);

    auto tk_xxy_yyyzzz = pbuffer.data(idx_kin_fi + 52);

    auto tk_xxy_yyzzzz = pbuffer.data(idx_kin_fi + 53);

    auto tk_xxy_yzzzzz = pbuffer.data(idx_kin_fi + 54);

    auto tk_xxy_zzzzzz = pbuffer.data(idx_kin_fi + 55);

#pragma omp simd aligned(pa_y,              \
                             tk_xx_xxxxx,   \
                             tk_xx_xxxxxx,  \
                             tk_xx_xxxxxy,  \
                             tk_xx_xxxxxz,  \
                             tk_xx_xxxxy,   \
                             tk_xx_xxxxyy,  \
                             tk_xx_xxxxyz,  \
                             tk_xx_xxxxz,   \
                             tk_xx_xxxxzz,  \
                             tk_xx_xxxyy,   \
                             tk_xx_xxxyyy,  \
                             tk_xx_xxxyyz,  \
                             tk_xx_xxxyz,   \
                             tk_xx_xxxyzz,  \
                             tk_xx_xxxzz,   \
                             tk_xx_xxxzzz,  \
                             tk_xx_xxyyy,   \
                             tk_xx_xxyyyy,  \
                             tk_xx_xxyyyz,  \
                             tk_xx_xxyyz,   \
                             tk_xx_xxyyzz,  \
                             tk_xx_xxyzz,   \
                             tk_xx_xxyzzz,  \
                             tk_xx_xxzzz,   \
                             tk_xx_xxzzzz,  \
                             tk_xx_xyyyy,   \
                             tk_xx_xyyyyy,  \
                             tk_xx_xyyyyz,  \
                             tk_xx_xyyyz,   \
                             tk_xx_xyyyzz,  \
                             tk_xx_xyyzz,   \
                             tk_xx_xyyzzz,  \
                             tk_xx_xyzzz,   \
                             tk_xx_xyzzzz,  \
                             tk_xx_xzzzz,   \
                             tk_xx_xzzzzz,  \
                             tk_xx_yyyyy,   \
                             tk_xx_yyyyyy,  \
                             tk_xx_yyyyyz,  \
                             tk_xx_yyyyz,   \
                             tk_xx_yyyyzz,  \
                             tk_xx_yyyzz,   \
                             tk_xx_yyyzzz,  \
                             tk_xx_yyzzz,   \
                             tk_xx_yyzzzz,  \
                             tk_xx_yzzzz,   \
                             tk_xx_yzzzzz,  \
                             tk_xx_zzzzz,   \
                             tk_xx_zzzzzz,  \
                             tk_xxy_xxxxxx, \
                             tk_xxy_xxxxxy, \
                             tk_xxy_xxxxxz, \
                             tk_xxy_xxxxyy, \
                             tk_xxy_xxxxyz, \
                             tk_xxy_xxxxzz, \
                             tk_xxy_xxxyyy, \
                             tk_xxy_xxxyyz, \
                             tk_xxy_xxxyzz, \
                             tk_xxy_xxxzzz, \
                             tk_xxy_xxyyyy, \
                             tk_xxy_xxyyyz, \
                             tk_xxy_xxyyzz, \
                             tk_xxy_xxyzzz, \
                             tk_xxy_xxzzzz, \
                             tk_xxy_xyyyyy, \
                             tk_xxy_xyyyyz, \
                             tk_xxy_xyyyzz, \
                             tk_xxy_xyyzzz, \
                             tk_xxy_xyzzzz, \
                             tk_xxy_xzzzzz, \
                             tk_xxy_yyyyyy, \
                             tk_xxy_yyyyyz, \
                             tk_xxy_yyyyzz, \
                             tk_xxy_yyyzzz, \
                             tk_xxy_yyzzzz, \
                             tk_xxy_yzzzzz, \
                             tk_xxy_zzzzzz, \
                             ts_xxy_xxxxxx, \
                             ts_xxy_xxxxxy, \
                             ts_xxy_xxxxxz, \
                             ts_xxy_xxxxyy, \
                             ts_xxy_xxxxyz, \
                             ts_xxy_xxxxzz, \
                             ts_xxy_xxxyyy, \
                             ts_xxy_xxxyyz, \
                             ts_xxy_xxxyzz, \
                             ts_xxy_xxxzzz, \
                             ts_xxy_xxyyyy, \
                             ts_xxy_xxyyyz, \
                             ts_xxy_xxyyzz, \
                             ts_xxy_xxyzzz, \
                             ts_xxy_xxzzzz, \
                             ts_xxy_xyyyyy, \
                             ts_xxy_xyyyyz, \
                             ts_xxy_xyyyzz, \
                             ts_xxy_xyyzzz, \
                             ts_xxy_xyzzzz, \
                             ts_xxy_xzzzzz, \
                             ts_xxy_yyyyyy, \
                             ts_xxy_yyyyyz, \
                             ts_xxy_yyyyzz, \
                             ts_xxy_yyyzzz, \
                             ts_xxy_yyzzzz, \
                             ts_xxy_yzzzzz, \
                             ts_xxy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxy_xxxxxx[i] = tk_xx_xxxxxx[i] * pa_y[i] + 2.0 * ts_xxy_xxxxxx[i] * fz_0;

        tk_xxy_xxxxxy[i] = tk_xx_xxxxx[i] * fe_0 + tk_xx_xxxxxy[i] * pa_y[i] + 2.0 * ts_xxy_xxxxxy[i] * fz_0;

        tk_xxy_xxxxxz[i] = tk_xx_xxxxxz[i] * pa_y[i] + 2.0 * ts_xxy_xxxxxz[i] * fz_0;

        tk_xxy_xxxxyy[i] = 2.0 * tk_xx_xxxxy[i] * fe_0 + tk_xx_xxxxyy[i] * pa_y[i] + 2.0 * ts_xxy_xxxxyy[i] * fz_0;

        tk_xxy_xxxxyz[i] = tk_xx_xxxxz[i] * fe_0 + tk_xx_xxxxyz[i] * pa_y[i] + 2.0 * ts_xxy_xxxxyz[i] * fz_0;

        tk_xxy_xxxxzz[i] = tk_xx_xxxxzz[i] * pa_y[i] + 2.0 * ts_xxy_xxxxzz[i] * fz_0;

        tk_xxy_xxxyyy[i] = 3.0 * tk_xx_xxxyy[i] * fe_0 + tk_xx_xxxyyy[i] * pa_y[i] + 2.0 * ts_xxy_xxxyyy[i] * fz_0;

        tk_xxy_xxxyyz[i] = 2.0 * tk_xx_xxxyz[i] * fe_0 + tk_xx_xxxyyz[i] * pa_y[i] + 2.0 * ts_xxy_xxxyyz[i] * fz_0;

        tk_xxy_xxxyzz[i] = tk_xx_xxxzz[i] * fe_0 + tk_xx_xxxyzz[i] * pa_y[i] + 2.0 * ts_xxy_xxxyzz[i] * fz_0;

        tk_xxy_xxxzzz[i] = tk_xx_xxxzzz[i] * pa_y[i] + 2.0 * ts_xxy_xxxzzz[i] * fz_0;

        tk_xxy_xxyyyy[i] = 4.0 * tk_xx_xxyyy[i] * fe_0 + tk_xx_xxyyyy[i] * pa_y[i] + 2.0 * ts_xxy_xxyyyy[i] * fz_0;

        tk_xxy_xxyyyz[i] = 3.0 * tk_xx_xxyyz[i] * fe_0 + tk_xx_xxyyyz[i] * pa_y[i] + 2.0 * ts_xxy_xxyyyz[i] * fz_0;

        tk_xxy_xxyyzz[i] = 2.0 * tk_xx_xxyzz[i] * fe_0 + tk_xx_xxyyzz[i] * pa_y[i] + 2.0 * ts_xxy_xxyyzz[i] * fz_0;

        tk_xxy_xxyzzz[i] = tk_xx_xxzzz[i] * fe_0 + tk_xx_xxyzzz[i] * pa_y[i] + 2.0 * ts_xxy_xxyzzz[i] * fz_0;

        tk_xxy_xxzzzz[i] = tk_xx_xxzzzz[i] * pa_y[i] + 2.0 * ts_xxy_xxzzzz[i] * fz_0;

        tk_xxy_xyyyyy[i] = 5.0 * tk_xx_xyyyy[i] * fe_0 + tk_xx_xyyyyy[i] * pa_y[i] + 2.0 * ts_xxy_xyyyyy[i] * fz_0;

        tk_xxy_xyyyyz[i] = 4.0 * tk_xx_xyyyz[i] * fe_0 + tk_xx_xyyyyz[i] * pa_y[i] + 2.0 * ts_xxy_xyyyyz[i] * fz_0;

        tk_xxy_xyyyzz[i] = 3.0 * tk_xx_xyyzz[i] * fe_0 + tk_xx_xyyyzz[i] * pa_y[i] + 2.0 * ts_xxy_xyyyzz[i] * fz_0;

        tk_xxy_xyyzzz[i] = 2.0 * tk_xx_xyzzz[i] * fe_0 + tk_xx_xyyzzz[i] * pa_y[i] + 2.0 * ts_xxy_xyyzzz[i] * fz_0;

        tk_xxy_xyzzzz[i] = tk_xx_xzzzz[i] * fe_0 + tk_xx_xyzzzz[i] * pa_y[i] + 2.0 * ts_xxy_xyzzzz[i] * fz_0;

        tk_xxy_xzzzzz[i] = tk_xx_xzzzzz[i] * pa_y[i] + 2.0 * ts_xxy_xzzzzz[i] * fz_0;

        tk_xxy_yyyyyy[i] = 6.0 * tk_xx_yyyyy[i] * fe_0 + tk_xx_yyyyyy[i] * pa_y[i] + 2.0 * ts_xxy_yyyyyy[i] * fz_0;

        tk_xxy_yyyyyz[i] = 5.0 * tk_xx_yyyyz[i] * fe_0 + tk_xx_yyyyyz[i] * pa_y[i] + 2.0 * ts_xxy_yyyyyz[i] * fz_0;

        tk_xxy_yyyyzz[i] = 4.0 * tk_xx_yyyzz[i] * fe_0 + tk_xx_yyyyzz[i] * pa_y[i] + 2.0 * ts_xxy_yyyyzz[i] * fz_0;

        tk_xxy_yyyzzz[i] = 3.0 * tk_xx_yyzzz[i] * fe_0 + tk_xx_yyyzzz[i] * pa_y[i] + 2.0 * ts_xxy_yyyzzz[i] * fz_0;

        tk_xxy_yyzzzz[i] = 2.0 * tk_xx_yzzzz[i] * fe_0 + tk_xx_yyzzzz[i] * pa_y[i] + 2.0 * ts_xxy_yyzzzz[i] * fz_0;

        tk_xxy_yzzzzz[i] = tk_xx_zzzzz[i] * fe_0 + tk_xx_yzzzzz[i] * pa_y[i] + 2.0 * ts_xxy_yzzzzz[i] * fz_0;

        tk_xxy_zzzzzz[i] = tk_xx_zzzzzz[i] * pa_y[i] + 2.0 * ts_xxy_zzzzzz[i] * fz_0;
    }

    // Set up 56-84 components of targeted buffer : FI

    auto tk_xxz_xxxxxx = pbuffer.data(idx_kin_fi + 56);

    auto tk_xxz_xxxxxy = pbuffer.data(idx_kin_fi + 57);

    auto tk_xxz_xxxxxz = pbuffer.data(idx_kin_fi + 58);

    auto tk_xxz_xxxxyy = pbuffer.data(idx_kin_fi + 59);

    auto tk_xxz_xxxxyz = pbuffer.data(idx_kin_fi + 60);

    auto tk_xxz_xxxxzz = pbuffer.data(idx_kin_fi + 61);

    auto tk_xxz_xxxyyy = pbuffer.data(idx_kin_fi + 62);

    auto tk_xxz_xxxyyz = pbuffer.data(idx_kin_fi + 63);

    auto tk_xxz_xxxyzz = pbuffer.data(idx_kin_fi + 64);

    auto tk_xxz_xxxzzz = pbuffer.data(idx_kin_fi + 65);

    auto tk_xxz_xxyyyy = pbuffer.data(idx_kin_fi + 66);

    auto tk_xxz_xxyyyz = pbuffer.data(idx_kin_fi + 67);

    auto tk_xxz_xxyyzz = pbuffer.data(idx_kin_fi + 68);

    auto tk_xxz_xxyzzz = pbuffer.data(idx_kin_fi + 69);

    auto tk_xxz_xxzzzz = pbuffer.data(idx_kin_fi + 70);

    auto tk_xxz_xyyyyy = pbuffer.data(idx_kin_fi + 71);

    auto tk_xxz_xyyyyz = pbuffer.data(idx_kin_fi + 72);

    auto tk_xxz_xyyyzz = pbuffer.data(idx_kin_fi + 73);

    auto tk_xxz_xyyzzz = pbuffer.data(idx_kin_fi + 74);

    auto tk_xxz_xyzzzz = pbuffer.data(idx_kin_fi + 75);

    auto tk_xxz_xzzzzz = pbuffer.data(idx_kin_fi + 76);

    auto tk_xxz_yyyyyy = pbuffer.data(idx_kin_fi + 77);

    auto tk_xxz_yyyyyz = pbuffer.data(idx_kin_fi + 78);

    auto tk_xxz_yyyyzz = pbuffer.data(idx_kin_fi + 79);

    auto tk_xxz_yyyzzz = pbuffer.data(idx_kin_fi + 80);

    auto tk_xxz_yyzzzz = pbuffer.data(idx_kin_fi + 81);

    auto tk_xxz_yzzzzz = pbuffer.data(idx_kin_fi + 82);

    auto tk_xxz_zzzzzz = pbuffer.data(idx_kin_fi + 83);

#pragma omp simd aligned(pa_z,              \
                             tk_xx_xxxxx,   \
                             tk_xx_xxxxxx,  \
                             tk_xx_xxxxxy,  \
                             tk_xx_xxxxxz,  \
                             tk_xx_xxxxy,   \
                             tk_xx_xxxxyy,  \
                             tk_xx_xxxxyz,  \
                             tk_xx_xxxxz,   \
                             tk_xx_xxxxzz,  \
                             tk_xx_xxxyy,   \
                             tk_xx_xxxyyy,  \
                             tk_xx_xxxyyz,  \
                             tk_xx_xxxyz,   \
                             tk_xx_xxxyzz,  \
                             tk_xx_xxxzz,   \
                             tk_xx_xxxzzz,  \
                             tk_xx_xxyyy,   \
                             tk_xx_xxyyyy,  \
                             tk_xx_xxyyyz,  \
                             tk_xx_xxyyz,   \
                             tk_xx_xxyyzz,  \
                             tk_xx_xxyzz,   \
                             tk_xx_xxyzzz,  \
                             tk_xx_xxzzz,   \
                             tk_xx_xxzzzz,  \
                             tk_xx_xyyyy,   \
                             tk_xx_xyyyyy,  \
                             tk_xx_xyyyyz,  \
                             tk_xx_xyyyz,   \
                             tk_xx_xyyyzz,  \
                             tk_xx_xyyzz,   \
                             tk_xx_xyyzzz,  \
                             tk_xx_xyzzz,   \
                             tk_xx_xyzzzz,  \
                             tk_xx_xzzzz,   \
                             tk_xx_xzzzzz,  \
                             tk_xx_yyyyy,   \
                             tk_xx_yyyyyy,  \
                             tk_xx_yyyyyz,  \
                             tk_xx_yyyyz,   \
                             tk_xx_yyyyzz,  \
                             tk_xx_yyyzz,   \
                             tk_xx_yyyzzz,  \
                             tk_xx_yyzzz,   \
                             tk_xx_yyzzzz,  \
                             tk_xx_yzzzz,   \
                             tk_xx_yzzzzz,  \
                             tk_xx_zzzzz,   \
                             tk_xx_zzzzzz,  \
                             tk_xxz_xxxxxx, \
                             tk_xxz_xxxxxy, \
                             tk_xxz_xxxxxz, \
                             tk_xxz_xxxxyy, \
                             tk_xxz_xxxxyz, \
                             tk_xxz_xxxxzz, \
                             tk_xxz_xxxyyy, \
                             tk_xxz_xxxyyz, \
                             tk_xxz_xxxyzz, \
                             tk_xxz_xxxzzz, \
                             tk_xxz_xxyyyy, \
                             tk_xxz_xxyyyz, \
                             tk_xxz_xxyyzz, \
                             tk_xxz_xxyzzz, \
                             tk_xxz_xxzzzz, \
                             tk_xxz_xyyyyy, \
                             tk_xxz_xyyyyz, \
                             tk_xxz_xyyyzz, \
                             tk_xxz_xyyzzz, \
                             tk_xxz_xyzzzz, \
                             tk_xxz_xzzzzz, \
                             tk_xxz_yyyyyy, \
                             tk_xxz_yyyyyz, \
                             tk_xxz_yyyyzz, \
                             tk_xxz_yyyzzz, \
                             tk_xxz_yyzzzz, \
                             tk_xxz_yzzzzz, \
                             tk_xxz_zzzzzz, \
                             ts_xxz_xxxxxx, \
                             ts_xxz_xxxxxy, \
                             ts_xxz_xxxxxz, \
                             ts_xxz_xxxxyy, \
                             ts_xxz_xxxxyz, \
                             ts_xxz_xxxxzz, \
                             ts_xxz_xxxyyy, \
                             ts_xxz_xxxyyz, \
                             ts_xxz_xxxyzz, \
                             ts_xxz_xxxzzz, \
                             ts_xxz_xxyyyy, \
                             ts_xxz_xxyyyz, \
                             ts_xxz_xxyyzz, \
                             ts_xxz_xxyzzz, \
                             ts_xxz_xxzzzz, \
                             ts_xxz_xyyyyy, \
                             ts_xxz_xyyyyz, \
                             ts_xxz_xyyyzz, \
                             ts_xxz_xyyzzz, \
                             ts_xxz_xyzzzz, \
                             ts_xxz_xzzzzz, \
                             ts_xxz_yyyyyy, \
                             ts_xxz_yyyyyz, \
                             ts_xxz_yyyyzz, \
                             ts_xxz_yyyzzz, \
                             ts_xxz_yyzzzz, \
                             ts_xxz_yzzzzz, \
                             ts_xxz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxz_xxxxxx[i] = tk_xx_xxxxxx[i] * pa_z[i] + 2.0 * ts_xxz_xxxxxx[i] * fz_0;

        tk_xxz_xxxxxy[i] = tk_xx_xxxxxy[i] * pa_z[i] + 2.0 * ts_xxz_xxxxxy[i] * fz_0;

        tk_xxz_xxxxxz[i] = tk_xx_xxxxx[i] * fe_0 + tk_xx_xxxxxz[i] * pa_z[i] + 2.0 * ts_xxz_xxxxxz[i] * fz_0;

        tk_xxz_xxxxyy[i] = tk_xx_xxxxyy[i] * pa_z[i] + 2.0 * ts_xxz_xxxxyy[i] * fz_0;

        tk_xxz_xxxxyz[i] = tk_xx_xxxxy[i] * fe_0 + tk_xx_xxxxyz[i] * pa_z[i] + 2.0 * ts_xxz_xxxxyz[i] * fz_0;

        tk_xxz_xxxxzz[i] = 2.0 * tk_xx_xxxxz[i] * fe_0 + tk_xx_xxxxzz[i] * pa_z[i] + 2.0 * ts_xxz_xxxxzz[i] * fz_0;

        tk_xxz_xxxyyy[i] = tk_xx_xxxyyy[i] * pa_z[i] + 2.0 * ts_xxz_xxxyyy[i] * fz_0;

        tk_xxz_xxxyyz[i] = tk_xx_xxxyy[i] * fe_0 + tk_xx_xxxyyz[i] * pa_z[i] + 2.0 * ts_xxz_xxxyyz[i] * fz_0;

        tk_xxz_xxxyzz[i] = 2.0 * tk_xx_xxxyz[i] * fe_0 + tk_xx_xxxyzz[i] * pa_z[i] + 2.0 * ts_xxz_xxxyzz[i] * fz_0;

        tk_xxz_xxxzzz[i] = 3.0 * tk_xx_xxxzz[i] * fe_0 + tk_xx_xxxzzz[i] * pa_z[i] + 2.0 * ts_xxz_xxxzzz[i] * fz_0;

        tk_xxz_xxyyyy[i] = tk_xx_xxyyyy[i] * pa_z[i] + 2.0 * ts_xxz_xxyyyy[i] * fz_0;

        tk_xxz_xxyyyz[i] = tk_xx_xxyyy[i] * fe_0 + tk_xx_xxyyyz[i] * pa_z[i] + 2.0 * ts_xxz_xxyyyz[i] * fz_0;

        tk_xxz_xxyyzz[i] = 2.0 * tk_xx_xxyyz[i] * fe_0 + tk_xx_xxyyzz[i] * pa_z[i] + 2.0 * ts_xxz_xxyyzz[i] * fz_0;

        tk_xxz_xxyzzz[i] = 3.0 * tk_xx_xxyzz[i] * fe_0 + tk_xx_xxyzzz[i] * pa_z[i] + 2.0 * ts_xxz_xxyzzz[i] * fz_0;

        tk_xxz_xxzzzz[i] = 4.0 * tk_xx_xxzzz[i] * fe_0 + tk_xx_xxzzzz[i] * pa_z[i] + 2.0 * ts_xxz_xxzzzz[i] * fz_0;

        tk_xxz_xyyyyy[i] = tk_xx_xyyyyy[i] * pa_z[i] + 2.0 * ts_xxz_xyyyyy[i] * fz_0;

        tk_xxz_xyyyyz[i] = tk_xx_xyyyy[i] * fe_0 + tk_xx_xyyyyz[i] * pa_z[i] + 2.0 * ts_xxz_xyyyyz[i] * fz_0;

        tk_xxz_xyyyzz[i] = 2.0 * tk_xx_xyyyz[i] * fe_0 + tk_xx_xyyyzz[i] * pa_z[i] + 2.0 * ts_xxz_xyyyzz[i] * fz_0;

        tk_xxz_xyyzzz[i] = 3.0 * tk_xx_xyyzz[i] * fe_0 + tk_xx_xyyzzz[i] * pa_z[i] + 2.0 * ts_xxz_xyyzzz[i] * fz_0;

        tk_xxz_xyzzzz[i] = 4.0 * tk_xx_xyzzz[i] * fe_0 + tk_xx_xyzzzz[i] * pa_z[i] + 2.0 * ts_xxz_xyzzzz[i] * fz_0;

        tk_xxz_xzzzzz[i] = 5.0 * tk_xx_xzzzz[i] * fe_0 + tk_xx_xzzzzz[i] * pa_z[i] + 2.0 * ts_xxz_xzzzzz[i] * fz_0;

        tk_xxz_yyyyyy[i] = tk_xx_yyyyyy[i] * pa_z[i] + 2.0 * ts_xxz_yyyyyy[i] * fz_0;

        tk_xxz_yyyyyz[i] = tk_xx_yyyyy[i] * fe_0 + tk_xx_yyyyyz[i] * pa_z[i] + 2.0 * ts_xxz_yyyyyz[i] * fz_0;

        tk_xxz_yyyyzz[i] = 2.0 * tk_xx_yyyyz[i] * fe_0 + tk_xx_yyyyzz[i] * pa_z[i] + 2.0 * ts_xxz_yyyyzz[i] * fz_0;

        tk_xxz_yyyzzz[i] = 3.0 * tk_xx_yyyzz[i] * fe_0 + tk_xx_yyyzzz[i] * pa_z[i] + 2.0 * ts_xxz_yyyzzz[i] * fz_0;

        tk_xxz_yyzzzz[i] = 4.0 * tk_xx_yyzzz[i] * fe_0 + tk_xx_yyzzzz[i] * pa_z[i] + 2.0 * ts_xxz_yyzzzz[i] * fz_0;

        tk_xxz_yzzzzz[i] = 5.0 * tk_xx_yzzzz[i] * fe_0 + tk_xx_yzzzzz[i] * pa_z[i] + 2.0 * ts_xxz_yzzzzz[i] * fz_0;

        tk_xxz_zzzzzz[i] = 6.0 * tk_xx_zzzzz[i] * fe_0 + tk_xx_zzzzzz[i] * pa_z[i] + 2.0 * ts_xxz_zzzzzz[i] * fz_0;
    }

    // Set up 84-112 components of targeted buffer : FI

    auto tk_xyy_xxxxxx = pbuffer.data(idx_kin_fi + 84);

    auto tk_xyy_xxxxxy = pbuffer.data(idx_kin_fi + 85);

    auto tk_xyy_xxxxxz = pbuffer.data(idx_kin_fi + 86);

    auto tk_xyy_xxxxyy = pbuffer.data(idx_kin_fi + 87);

    auto tk_xyy_xxxxyz = pbuffer.data(idx_kin_fi + 88);

    auto tk_xyy_xxxxzz = pbuffer.data(idx_kin_fi + 89);

    auto tk_xyy_xxxyyy = pbuffer.data(idx_kin_fi + 90);

    auto tk_xyy_xxxyyz = pbuffer.data(idx_kin_fi + 91);

    auto tk_xyy_xxxyzz = pbuffer.data(idx_kin_fi + 92);

    auto tk_xyy_xxxzzz = pbuffer.data(idx_kin_fi + 93);

    auto tk_xyy_xxyyyy = pbuffer.data(idx_kin_fi + 94);

    auto tk_xyy_xxyyyz = pbuffer.data(idx_kin_fi + 95);

    auto tk_xyy_xxyyzz = pbuffer.data(idx_kin_fi + 96);

    auto tk_xyy_xxyzzz = pbuffer.data(idx_kin_fi + 97);

    auto tk_xyy_xxzzzz = pbuffer.data(idx_kin_fi + 98);

    auto tk_xyy_xyyyyy = pbuffer.data(idx_kin_fi + 99);

    auto tk_xyy_xyyyyz = pbuffer.data(idx_kin_fi + 100);

    auto tk_xyy_xyyyzz = pbuffer.data(idx_kin_fi + 101);

    auto tk_xyy_xyyzzz = pbuffer.data(idx_kin_fi + 102);

    auto tk_xyy_xyzzzz = pbuffer.data(idx_kin_fi + 103);

    auto tk_xyy_xzzzzz = pbuffer.data(idx_kin_fi + 104);

    auto tk_xyy_yyyyyy = pbuffer.data(idx_kin_fi + 105);

    auto tk_xyy_yyyyyz = pbuffer.data(idx_kin_fi + 106);

    auto tk_xyy_yyyyzz = pbuffer.data(idx_kin_fi + 107);

    auto tk_xyy_yyyzzz = pbuffer.data(idx_kin_fi + 108);

    auto tk_xyy_yyzzzz = pbuffer.data(idx_kin_fi + 109);

    auto tk_xyy_yzzzzz = pbuffer.data(idx_kin_fi + 110);

    auto tk_xyy_zzzzzz = pbuffer.data(idx_kin_fi + 111);

#pragma omp simd aligned(pa_x,              \
                             tk_xyy_xxxxxx, \
                             tk_xyy_xxxxxy, \
                             tk_xyy_xxxxxz, \
                             tk_xyy_xxxxyy, \
                             tk_xyy_xxxxyz, \
                             tk_xyy_xxxxzz, \
                             tk_xyy_xxxyyy, \
                             tk_xyy_xxxyyz, \
                             tk_xyy_xxxyzz, \
                             tk_xyy_xxxzzz, \
                             tk_xyy_xxyyyy, \
                             tk_xyy_xxyyyz, \
                             tk_xyy_xxyyzz, \
                             tk_xyy_xxyzzz, \
                             tk_xyy_xxzzzz, \
                             tk_xyy_xyyyyy, \
                             tk_xyy_xyyyyz, \
                             tk_xyy_xyyyzz, \
                             tk_xyy_xyyzzz, \
                             tk_xyy_xyzzzz, \
                             tk_xyy_xzzzzz, \
                             tk_xyy_yyyyyy, \
                             tk_xyy_yyyyyz, \
                             tk_xyy_yyyyzz, \
                             tk_xyy_yyyzzz, \
                             tk_xyy_yyzzzz, \
                             tk_xyy_yzzzzz, \
                             tk_xyy_zzzzzz, \
                             tk_yy_xxxxx,   \
                             tk_yy_xxxxxx,  \
                             tk_yy_xxxxxy,  \
                             tk_yy_xxxxxz,  \
                             tk_yy_xxxxy,   \
                             tk_yy_xxxxyy,  \
                             tk_yy_xxxxyz,  \
                             tk_yy_xxxxz,   \
                             tk_yy_xxxxzz,  \
                             tk_yy_xxxyy,   \
                             tk_yy_xxxyyy,  \
                             tk_yy_xxxyyz,  \
                             tk_yy_xxxyz,   \
                             tk_yy_xxxyzz,  \
                             tk_yy_xxxzz,   \
                             tk_yy_xxxzzz,  \
                             tk_yy_xxyyy,   \
                             tk_yy_xxyyyy,  \
                             tk_yy_xxyyyz,  \
                             tk_yy_xxyyz,   \
                             tk_yy_xxyyzz,  \
                             tk_yy_xxyzz,   \
                             tk_yy_xxyzzz,  \
                             tk_yy_xxzzz,   \
                             tk_yy_xxzzzz,  \
                             tk_yy_xyyyy,   \
                             tk_yy_xyyyyy,  \
                             tk_yy_xyyyyz,  \
                             tk_yy_xyyyz,   \
                             tk_yy_xyyyzz,  \
                             tk_yy_xyyzz,   \
                             tk_yy_xyyzzz,  \
                             tk_yy_xyzzz,   \
                             tk_yy_xyzzzz,  \
                             tk_yy_xzzzz,   \
                             tk_yy_xzzzzz,  \
                             tk_yy_yyyyy,   \
                             tk_yy_yyyyyy,  \
                             tk_yy_yyyyyz,  \
                             tk_yy_yyyyz,   \
                             tk_yy_yyyyzz,  \
                             tk_yy_yyyzz,   \
                             tk_yy_yyyzzz,  \
                             tk_yy_yyzzz,   \
                             tk_yy_yyzzzz,  \
                             tk_yy_yzzzz,   \
                             tk_yy_yzzzzz,  \
                             tk_yy_zzzzz,   \
                             tk_yy_zzzzzz,  \
                             ts_xyy_xxxxxx, \
                             ts_xyy_xxxxxy, \
                             ts_xyy_xxxxxz, \
                             ts_xyy_xxxxyy, \
                             ts_xyy_xxxxyz, \
                             ts_xyy_xxxxzz, \
                             ts_xyy_xxxyyy, \
                             ts_xyy_xxxyyz, \
                             ts_xyy_xxxyzz, \
                             ts_xyy_xxxzzz, \
                             ts_xyy_xxyyyy, \
                             ts_xyy_xxyyyz, \
                             ts_xyy_xxyyzz, \
                             ts_xyy_xxyzzz, \
                             ts_xyy_xxzzzz, \
                             ts_xyy_xyyyyy, \
                             ts_xyy_xyyyyz, \
                             ts_xyy_xyyyzz, \
                             ts_xyy_xyyzzz, \
                             ts_xyy_xyzzzz, \
                             ts_xyy_xzzzzz, \
                             ts_xyy_yyyyyy, \
                             ts_xyy_yyyyyz, \
                             ts_xyy_yyyyzz, \
                             ts_xyy_yyyzzz, \
                             ts_xyy_yyzzzz, \
                             ts_xyy_yzzzzz, \
                             ts_xyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyy_xxxxxx[i] = 6.0 * tk_yy_xxxxx[i] * fe_0 + tk_yy_xxxxxx[i] * pa_x[i] + 2.0 * ts_xyy_xxxxxx[i] * fz_0;

        tk_xyy_xxxxxy[i] = 5.0 * tk_yy_xxxxy[i] * fe_0 + tk_yy_xxxxxy[i] * pa_x[i] + 2.0 * ts_xyy_xxxxxy[i] * fz_0;

        tk_xyy_xxxxxz[i] = 5.0 * tk_yy_xxxxz[i] * fe_0 + tk_yy_xxxxxz[i] * pa_x[i] + 2.0 * ts_xyy_xxxxxz[i] * fz_0;

        tk_xyy_xxxxyy[i] = 4.0 * tk_yy_xxxyy[i] * fe_0 + tk_yy_xxxxyy[i] * pa_x[i] + 2.0 * ts_xyy_xxxxyy[i] * fz_0;

        tk_xyy_xxxxyz[i] = 4.0 * tk_yy_xxxyz[i] * fe_0 + tk_yy_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyy_xxxxyz[i] * fz_0;

        tk_xyy_xxxxzz[i] = 4.0 * tk_yy_xxxzz[i] * fe_0 + tk_yy_xxxxzz[i] * pa_x[i] + 2.0 * ts_xyy_xxxxzz[i] * fz_0;

        tk_xyy_xxxyyy[i] = 3.0 * tk_yy_xxyyy[i] * fe_0 + tk_yy_xxxyyy[i] * pa_x[i] + 2.0 * ts_xyy_xxxyyy[i] * fz_0;

        tk_xyy_xxxyyz[i] = 3.0 * tk_yy_xxyyz[i] * fe_0 + tk_yy_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyy_xxxyyz[i] * fz_0;

        tk_xyy_xxxyzz[i] = 3.0 * tk_yy_xxyzz[i] * fe_0 + tk_yy_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyy_xxxyzz[i] * fz_0;

        tk_xyy_xxxzzz[i] = 3.0 * tk_yy_xxzzz[i] * fe_0 + tk_yy_xxxzzz[i] * pa_x[i] + 2.0 * ts_xyy_xxxzzz[i] * fz_0;

        tk_xyy_xxyyyy[i] = 2.0 * tk_yy_xyyyy[i] * fe_0 + tk_yy_xxyyyy[i] * pa_x[i] + 2.0 * ts_xyy_xxyyyy[i] * fz_0;

        tk_xyy_xxyyyz[i] = 2.0 * tk_yy_xyyyz[i] * fe_0 + tk_yy_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyy_xxyyyz[i] * fz_0;

        tk_xyy_xxyyzz[i] = 2.0 * tk_yy_xyyzz[i] * fe_0 + tk_yy_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyy_xxyyzz[i] * fz_0;

        tk_xyy_xxyzzz[i] = 2.0 * tk_yy_xyzzz[i] * fe_0 + tk_yy_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyy_xxyzzz[i] * fz_0;

        tk_xyy_xxzzzz[i] = 2.0 * tk_yy_xzzzz[i] * fe_0 + tk_yy_xxzzzz[i] * pa_x[i] + 2.0 * ts_xyy_xxzzzz[i] * fz_0;

        tk_xyy_xyyyyy[i] = tk_yy_yyyyy[i] * fe_0 + tk_yy_xyyyyy[i] * pa_x[i] + 2.0 * ts_xyy_xyyyyy[i] * fz_0;

        tk_xyy_xyyyyz[i] = tk_yy_yyyyz[i] * fe_0 + tk_yy_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyy_xyyyyz[i] * fz_0;

        tk_xyy_xyyyzz[i] = tk_yy_yyyzz[i] * fe_0 + tk_yy_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyy_xyyyzz[i] * fz_0;

        tk_xyy_xyyzzz[i] = tk_yy_yyzzz[i] * fe_0 + tk_yy_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyy_xyyzzz[i] * fz_0;

        tk_xyy_xyzzzz[i] = tk_yy_yzzzz[i] * fe_0 + tk_yy_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyy_xyzzzz[i] * fz_0;

        tk_xyy_xzzzzz[i] = tk_yy_zzzzz[i] * fe_0 + tk_yy_xzzzzz[i] * pa_x[i] + 2.0 * ts_xyy_xzzzzz[i] * fz_0;

        tk_xyy_yyyyyy[i] = tk_yy_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyy_yyyyyy[i] * fz_0;

        tk_xyy_yyyyyz[i] = tk_yy_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyy_yyyyyz[i] * fz_0;

        tk_xyy_yyyyzz[i] = tk_yy_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyy_yyyyzz[i] * fz_0;

        tk_xyy_yyyzzz[i] = tk_yy_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyy_yyyzzz[i] * fz_0;

        tk_xyy_yyzzzz[i] = tk_yy_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyy_yyzzzz[i] * fz_0;

        tk_xyy_yzzzzz[i] = tk_yy_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyy_yzzzzz[i] * fz_0;

        tk_xyy_zzzzzz[i] = tk_yy_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyy_zzzzzz[i] * fz_0;
    }

    // Set up 112-140 components of targeted buffer : FI

    auto tk_xyz_xxxxxx = pbuffer.data(idx_kin_fi + 112);

    auto tk_xyz_xxxxxy = pbuffer.data(idx_kin_fi + 113);

    auto tk_xyz_xxxxxz = pbuffer.data(idx_kin_fi + 114);

    auto tk_xyz_xxxxyy = pbuffer.data(idx_kin_fi + 115);

    auto tk_xyz_xxxxyz = pbuffer.data(idx_kin_fi + 116);

    auto tk_xyz_xxxxzz = pbuffer.data(idx_kin_fi + 117);

    auto tk_xyz_xxxyyy = pbuffer.data(idx_kin_fi + 118);

    auto tk_xyz_xxxyyz = pbuffer.data(idx_kin_fi + 119);

    auto tk_xyz_xxxyzz = pbuffer.data(idx_kin_fi + 120);

    auto tk_xyz_xxxzzz = pbuffer.data(idx_kin_fi + 121);

    auto tk_xyz_xxyyyy = pbuffer.data(idx_kin_fi + 122);

    auto tk_xyz_xxyyyz = pbuffer.data(idx_kin_fi + 123);

    auto tk_xyz_xxyyzz = pbuffer.data(idx_kin_fi + 124);

    auto tk_xyz_xxyzzz = pbuffer.data(idx_kin_fi + 125);

    auto tk_xyz_xxzzzz = pbuffer.data(idx_kin_fi + 126);

    auto tk_xyz_xyyyyy = pbuffer.data(idx_kin_fi + 127);

    auto tk_xyz_xyyyyz = pbuffer.data(idx_kin_fi + 128);

    auto tk_xyz_xyyyzz = pbuffer.data(idx_kin_fi + 129);

    auto tk_xyz_xyyzzz = pbuffer.data(idx_kin_fi + 130);

    auto tk_xyz_xyzzzz = pbuffer.data(idx_kin_fi + 131);

    auto tk_xyz_xzzzzz = pbuffer.data(idx_kin_fi + 132);

    auto tk_xyz_yyyyyy = pbuffer.data(idx_kin_fi + 133);

    auto tk_xyz_yyyyyz = pbuffer.data(idx_kin_fi + 134);

    auto tk_xyz_yyyyzz = pbuffer.data(idx_kin_fi + 135);

    auto tk_xyz_yyyzzz = pbuffer.data(idx_kin_fi + 136);

    auto tk_xyz_yyzzzz = pbuffer.data(idx_kin_fi + 137);

    auto tk_xyz_yzzzzz = pbuffer.data(idx_kin_fi + 138);

    auto tk_xyz_zzzzzz = pbuffer.data(idx_kin_fi + 139);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             tk_xy_xxxxxy,  \
                             tk_xy_xxxxyy,  \
                             tk_xy_xxxyyy,  \
                             tk_xy_xxyyyy,  \
                             tk_xy_xyyyyy,  \
                             tk_xyz_xxxxxx, \
                             tk_xyz_xxxxxy, \
                             tk_xyz_xxxxxz, \
                             tk_xyz_xxxxyy, \
                             tk_xyz_xxxxyz, \
                             tk_xyz_xxxxzz, \
                             tk_xyz_xxxyyy, \
                             tk_xyz_xxxyyz, \
                             tk_xyz_xxxyzz, \
                             tk_xyz_xxxzzz, \
                             tk_xyz_xxyyyy, \
                             tk_xyz_xxyyyz, \
                             tk_xyz_xxyyzz, \
                             tk_xyz_xxyzzz, \
                             tk_xyz_xxzzzz, \
                             tk_xyz_xyyyyy, \
                             tk_xyz_xyyyyz, \
                             tk_xyz_xyyyzz, \
                             tk_xyz_xyyzzz, \
                             tk_xyz_xyzzzz, \
                             tk_xyz_xzzzzz, \
                             tk_xyz_yyyyyy, \
                             tk_xyz_yyyyyz, \
                             tk_xyz_yyyyzz, \
                             tk_xyz_yyyzzz, \
                             tk_xyz_yyzzzz, \
                             tk_xyz_yzzzzz, \
                             tk_xyz_zzzzzz, \
                             tk_xz_xxxxxx,  \
                             tk_xz_xxxxxz,  \
                             tk_xz_xxxxzz,  \
                             tk_xz_xxxzzz,  \
                             tk_xz_xxzzzz,  \
                             tk_xz_xzzzzz,  \
                             tk_yz_xxxxyz,  \
                             tk_yz_xxxyyz,  \
                             tk_yz_xxxyz,   \
                             tk_yz_xxxyzz,  \
                             tk_yz_xxyyyz,  \
                             tk_yz_xxyyz,   \
                             tk_yz_xxyyzz,  \
                             tk_yz_xxyzz,   \
                             tk_yz_xxyzzz,  \
                             tk_yz_xyyyyz,  \
                             tk_yz_xyyyz,   \
                             tk_yz_xyyyzz,  \
                             tk_yz_xyyzz,   \
                             tk_yz_xyyzzz,  \
                             tk_yz_xyzzz,   \
                             tk_yz_xyzzzz,  \
                             tk_yz_yyyyyy,  \
                             tk_yz_yyyyyz,  \
                             tk_yz_yyyyz,   \
                             tk_yz_yyyyzz,  \
                             tk_yz_yyyzz,   \
                             tk_yz_yyyzzz,  \
                             tk_yz_yyzzz,   \
                             tk_yz_yyzzzz,  \
                             tk_yz_yzzzz,   \
                             tk_yz_yzzzzz,  \
                             tk_yz_zzzzzz,  \
                             ts_xyz_xxxxxx, \
                             ts_xyz_xxxxxy, \
                             ts_xyz_xxxxxz, \
                             ts_xyz_xxxxyy, \
                             ts_xyz_xxxxyz, \
                             ts_xyz_xxxxzz, \
                             ts_xyz_xxxyyy, \
                             ts_xyz_xxxyyz, \
                             ts_xyz_xxxyzz, \
                             ts_xyz_xxxzzz, \
                             ts_xyz_xxyyyy, \
                             ts_xyz_xxyyyz, \
                             ts_xyz_xxyyzz, \
                             ts_xyz_xxyzzz, \
                             ts_xyz_xxzzzz, \
                             ts_xyz_xyyyyy, \
                             ts_xyz_xyyyyz, \
                             ts_xyz_xyyyzz, \
                             ts_xyz_xyyzzz, \
                             ts_xyz_xyzzzz, \
                             ts_xyz_xzzzzz, \
                             ts_xyz_yyyyyy, \
                             ts_xyz_yyyyyz, \
                             ts_xyz_yyyyzz, \
                             ts_xyz_yyyzzz, \
                             ts_xyz_yyzzzz, \
                             ts_xyz_yzzzzz, \
                             ts_xyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyz_xxxxxx[i] = tk_xz_xxxxxx[i] * pa_y[i] + 2.0 * ts_xyz_xxxxxx[i] * fz_0;

        tk_xyz_xxxxxy[i] = tk_xy_xxxxxy[i] * pa_z[i] + 2.0 * ts_xyz_xxxxxy[i] * fz_0;

        tk_xyz_xxxxxz[i] = tk_xz_xxxxxz[i] * pa_y[i] + 2.0 * ts_xyz_xxxxxz[i] * fz_0;

        tk_xyz_xxxxyy[i] = tk_xy_xxxxyy[i] * pa_z[i] + 2.0 * ts_xyz_xxxxyy[i] * fz_0;

        tk_xyz_xxxxyz[i] = 4.0 * tk_yz_xxxyz[i] * fe_0 + tk_yz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xyz_xxxxyz[i] * fz_0;

        tk_xyz_xxxxzz[i] = tk_xz_xxxxzz[i] * pa_y[i] + 2.0 * ts_xyz_xxxxzz[i] * fz_0;

        tk_xyz_xxxyyy[i] = tk_xy_xxxyyy[i] * pa_z[i] + 2.0 * ts_xyz_xxxyyy[i] * fz_0;

        tk_xyz_xxxyyz[i] = 3.0 * tk_yz_xxyyz[i] * fe_0 + tk_yz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xyz_xxxyyz[i] * fz_0;

        tk_xyz_xxxyzz[i] = 3.0 * tk_yz_xxyzz[i] * fe_0 + tk_yz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xyz_xxxyzz[i] * fz_0;

        tk_xyz_xxxzzz[i] = tk_xz_xxxzzz[i] * pa_y[i] + 2.0 * ts_xyz_xxxzzz[i] * fz_0;

        tk_xyz_xxyyyy[i] = tk_xy_xxyyyy[i] * pa_z[i] + 2.0 * ts_xyz_xxyyyy[i] * fz_0;

        tk_xyz_xxyyyz[i] = 2.0 * tk_yz_xyyyz[i] * fe_0 + tk_yz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xyz_xxyyyz[i] * fz_0;

        tk_xyz_xxyyzz[i] = 2.0 * tk_yz_xyyzz[i] * fe_0 + tk_yz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xyz_xxyyzz[i] * fz_0;

        tk_xyz_xxyzzz[i] = 2.0 * tk_yz_xyzzz[i] * fe_0 + tk_yz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xyz_xxyzzz[i] * fz_0;

        tk_xyz_xxzzzz[i] = tk_xz_xxzzzz[i] * pa_y[i] + 2.0 * ts_xyz_xxzzzz[i] * fz_0;

        tk_xyz_xyyyyy[i] = tk_xy_xyyyyy[i] * pa_z[i] + 2.0 * ts_xyz_xyyyyy[i] * fz_0;

        tk_xyz_xyyyyz[i] = tk_yz_yyyyz[i] * fe_0 + tk_yz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xyz_xyyyyz[i] * fz_0;

        tk_xyz_xyyyzz[i] = tk_yz_yyyzz[i] * fe_0 + tk_yz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xyz_xyyyzz[i] * fz_0;

        tk_xyz_xyyzzz[i] = tk_yz_yyzzz[i] * fe_0 + tk_yz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xyz_xyyzzz[i] * fz_0;

        tk_xyz_xyzzzz[i] = tk_yz_yzzzz[i] * fe_0 + tk_yz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xyz_xyzzzz[i] * fz_0;

        tk_xyz_xzzzzz[i] = tk_xz_xzzzzz[i] * pa_y[i] + 2.0 * ts_xyz_xzzzzz[i] * fz_0;

        tk_xyz_yyyyyy[i] = tk_yz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xyz_yyyyyy[i] * fz_0;

        tk_xyz_yyyyyz[i] = tk_yz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xyz_yyyyyz[i] * fz_0;

        tk_xyz_yyyyzz[i] = tk_yz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xyz_yyyyzz[i] * fz_0;

        tk_xyz_yyyzzz[i] = tk_yz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xyz_yyyzzz[i] * fz_0;

        tk_xyz_yyzzzz[i] = tk_yz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xyz_yyzzzz[i] * fz_0;

        tk_xyz_yzzzzz[i] = tk_yz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xyz_yzzzzz[i] * fz_0;

        tk_xyz_zzzzzz[i] = tk_yz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xyz_zzzzzz[i] * fz_0;
    }

    // Set up 140-168 components of targeted buffer : FI

    auto tk_xzz_xxxxxx = pbuffer.data(idx_kin_fi + 140);

    auto tk_xzz_xxxxxy = pbuffer.data(idx_kin_fi + 141);

    auto tk_xzz_xxxxxz = pbuffer.data(idx_kin_fi + 142);

    auto tk_xzz_xxxxyy = pbuffer.data(idx_kin_fi + 143);

    auto tk_xzz_xxxxyz = pbuffer.data(idx_kin_fi + 144);

    auto tk_xzz_xxxxzz = pbuffer.data(idx_kin_fi + 145);

    auto tk_xzz_xxxyyy = pbuffer.data(idx_kin_fi + 146);

    auto tk_xzz_xxxyyz = pbuffer.data(idx_kin_fi + 147);

    auto tk_xzz_xxxyzz = pbuffer.data(idx_kin_fi + 148);

    auto tk_xzz_xxxzzz = pbuffer.data(idx_kin_fi + 149);

    auto tk_xzz_xxyyyy = pbuffer.data(idx_kin_fi + 150);

    auto tk_xzz_xxyyyz = pbuffer.data(idx_kin_fi + 151);

    auto tk_xzz_xxyyzz = pbuffer.data(idx_kin_fi + 152);

    auto tk_xzz_xxyzzz = pbuffer.data(idx_kin_fi + 153);

    auto tk_xzz_xxzzzz = pbuffer.data(idx_kin_fi + 154);

    auto tk_xzz_xyyyyy = pbuffer.data(idx_kin_fi + 155);

    auto tk_xzz_xyyyyz = pbuffer.data(idx_kin_fi + 156);

    auto tk_xzz_xyyyzz = pbuffer.data(idx_kin_fi + 157);

    auto tk_xzz_xyyzzz = pbuffer.data(idx_kin_fi + 158);

    auto tk_xzz_xyzzzz = pbuffer.data(idx_kin_fi + 159);

    auto tk_xzz_xzzzzz = pbuffer.data(idx_kin_fi + 160);

    auto tk_xzz_yyyyyy = pbuffer.data(idx_kin_fi + 161);

    auto tk_xzz_yyyyyz = pbuffer.data(idx_kin_fi + 162);

    auto tk_xzz_yyyyzz = pbuffer.data(idx_kin_fi + 163);

    auto tk_xzz_yyyzzz = pbuffer.data(idx_kin_fi + 164);

    auto tk_xzz_yyzzzz = pbuffer.data(idx_kin_fi + 165);

    auto tk_xzz_yzzzzz = pbuffer.data(idx_kin_fi + 166);

    auto tk_xzz_zzzzzz = pbuffer.data(idx_kin_fi + 167);

#pragma omp simd aligned(pa_x,              \
                             tk_xzz_xxxxxx, \
                             tk_xzz_xxxxxy, \
                             tk_xzz_xxxxxz, \
                             tk_xzz_xxxxyy, \
                             tk_xzz_xxxxyz, \
                             tk_xzz_xxxxzz, \
                             tk_xzz_xxxyyy, \
                             tk_xzz_xxxyyz, \
                             tk_xzz_xxxyzz, \
                             tk_xzz_xxxzzz, \
                             tk_xzz_xxyyyy, \
                             tk_xzz_xxyyyz, \
                             tk_xzz_xxyyzz, \
                             tk_xzz_xxyzzz, \
                             tk_xzz_xxzzzz, \
                             tk_xzz_xyyyyy, \
                             tk_xzz_xyyyyz, \
                             tk_xzz_xyyyzz, \
                             tk_xzz_xyyzzz, \
                             tk_xzz_xyzzzz, \
                             tk_xzz_xzzzzz, \
                             tk_xzz_yyyyyy, \
                             tk_xzz_yyyyyz, \
                             tk_xzz_yyyyzz, \
                             tk_xzz_yyyzzz, \
                             tk_xzz_yyzzzz, \
                             tk_xzz_yzzzzz, \
                             tk_xzz_zzzzzz, \
                             tk_zz_xxxxx,   \
                             tk_zz_xxxxxx,  \
                             tk_zz_xxxxxy,  \
                             tk_zz_xxxxxz,  \
                             tk_zz_xxxxy,   \
                             tk_zz_xxxxyy,  \
                             tk_zz_xxxxyz,  \
                             tk_zz_xxxxz,   \
                             tk_zz_xxxxzz,  \
                             tk_zz_xxxyy,   \
                             tk_zz_xxxyyy,  \
                             tk_zz_xxxyyz,  \
                             tk_zz_xxxyz,   \
                             tk_zz_xxxyzz,  \
                             tk_zz_xxxzz,   \
                             tk_zz_xxxzzz,  \
                             tk_zz_xxyyy,   \
                             tk_zz_xxyyyy,  \
                             tk_zz_xxyyyz,  \
                             tk_zz_xxyyz,   \
                             tk_zz_xxyyzz,  \
                             tk_zz_xxyzz,   \
                             tk_zz_xxyzzz,  \
                             tk_zz_xxzzz,   \
                             tk_zz_xxzzzz,  \
                             tk_zz_xyyyy,   \
                             tk_zz_xyyyyy,  \
                             tk_zz_xyyyyz,  \
                             tk_zz_xyyyz,   \
                             tk_zz_xyyyzz,  \
                             tk_zz_xyyzz,   \
                             tk_zz_xyyzzz,  \
                             tk_zz_xyzzz,   \
                             tk_zz_xyzzzz,  \
                             tk_zz_xzzzz,   \
                             tk_zz_xzzzzz,  \
                             tk_zz_yyyyy,   \
                             tk_zz_yyyyyy,  \
                             tk_zz_yyyyyz,  \
                             tk_zz_yyyyz,   \
                             tk_zz_yyyyzz,  \
                             tk_zz_yyyzz,   \
                             tk_zz_yyyzzz,  \
                             tk_zz_yyzzz,   \
                             tk_zz_yyzzzz,  \
                             tk_zz_yzzzz,   \
                             tk_zz_yzzzzz,  \
                             tk_zz_zzzzz,   \
                             tk_zz_zzzzzz,  \
                             ts_xzz_xxxxxx, \
                             ts_xzz_xxxxxy, \
                             ts_xzz_xxxxxz, \
                             ts_xzz_xxxxyy, \
                             ts_xzz_xxxxyz, \
                             ts_xzz_xxxxzz, \
                             ts_xzz_xxxyyy, \
                             ts_xzz_xxxyyz, \
                             ts_xzz_xxxyzz, \
                             ts_xzz_xxxzzz, \
                             ts_xzz_xxyyyy, \
                             ts_xzz_xxyyyz, \
                             ts_xzz_xxyyzz, \
                             ts_xzz_xxyzzz, \
                             ts_xzz_xxzzzz, \
                             ts_xzz_xyyyyy, \
                             ts_xzz_xyyyyz, \
                             ts_xzz_xyyyzz, \
                             ts_xzz_xyyzzz, \
                             ts_xzz_xyzzzz, \
                             ts_xzz_xzzzzz, \
                             ts_xzz_yyyyyy, \
                             ts_xzz_yyyyyz, \
                             ts_xzz_yyyyzz, \
                             ts_xzz_yyyzzz, \
                             ts_xzz_yyzzzz, \
                             ts_xzz_yzzzzz, \
                             ts_xzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzz_xxxxxx[i] = 6.0 * tk_zz_xxxxx[i] * fe_0 + tk_zz_xxxxxx[i] * pa_x[i] + 2.0 * ts_xzz_xxxxxx[i] * fz_0;

        tk_xzz_xxxxxy[i] = 5.0 * tk_zz_xxxxy[i] * fe_0 + tk_zz_xxxxxy[i] * pa_x[i] + 2.0 * ts_xzz_xxxxxy[i] * fz_0;

        tk_xzz_xxxxxz[i] = 5.0 * tk_zz_xxxxz[i] * fe_0 + tk_zz_xxxxxz[i] * pa_x[i] + 2.0 * ts_xzz_xxxxxz[i] * fz_0;

        tk_xzz_xxxxyy[i] = 4.0 * tk_zz_xxxyy[i] * fe_0 + tk_zz_xxxxyy[i] * pa_x[i] + 2.0 * ts_xzz_xxxxyy[i] * fz_0;

        tk_xzz_xxxxyz[i] = 4.0 * tk_zz_xxxyz[i] * fe_0 + tk_zz_xxxxyz[i] * pa_x[i] + 2.0 * ts_xzz_xxxxyz[i] * fz_0;

        tk_xzz_xxxxzz[i] = 4.0 * tk_zz_xxxzz[i] * fe_0 + tk_zz_xxxxzz[i] * pa_x[i] + 2.0 * ts_xzz_xxxxzz[i] * fz_0;

        tk_xzz_xxxyyy[i] = 3.0 * tk_zz_xxyyy[i] * fe_0 + tk_zz_xxxyyy[i] * pa_x[i] + 2.0 * ts_xzz_xxxyyy[i] * fz_0;

        tk_xzz_xxxyyz[i] = 3.0 * tk_zz_xxyyz[i] * fe_0 + tk_zz_xxxyyz[i] * pa_x[i] + 2.0 * ts_xzz_xxxyyz[i] * fz_0;

        tk_xzz_xxxyzz[i] = 3.0 * tk_zz_xxyzz[i] * fe_0 + tk_zz_xxxyzz[i] * pa_x[i] + 2.0 * ts_xzz_xxxyzz[i] * fz_0;

        tk_xzz_xxxzzz[i] = 3.0 * tk_zz_xxzzz[i] * fe_0 + tk_zz_xxxzzz[i] * pa_x[i] + 2.0 * ts_xzz_xxxzzz[i] * fz_0;

        tk_xzz_xxyyyy[i] = 2.0 * tk_zz_xyyyy[i] * fe_0 + tk_zz_xxyyyy[i] * pa_x[i] + 2.0 * ts_xzz_xxyyyy[i] * fz_0;

        tk_xzz_xxyyyz[i] = 2.0 * tk_zz_xyyyz[i] * fe_0 + tk_zz_xxyyyz[i] * pa_x[i] + 2.0 * ts_xzz_xxyyyz[i] * fz_0;

        tk_xzz_xxyyzz[i] = 2.0 * tk_zz_xyyzz[i] * fe_0 + tk_zz_xxyyzz[i] * pa_x[i] + 2.0 * ts_xzz_xxyyzz[i] * fz_0;

        tk_xzz_xxyzzz[i] = 2.0 * tk_zz_xyzzz[i] * fe_0 + tk_zz_xxyzzz[i] * pa_x[i] + 2.0 * ts_xzz_xxyzzz[i] * fz_0;

        tk_xzz_xxzzzz[i] = 2.0 * tk_zz_xzzzz[i] * fe_0 + tk_zz_xxzzzz[i] * pa_x[i] + 2.0 * ts_xzz_xxzzzz[i] * fz_0;

        tk_xzz_xyyyyy[i] = tk_zz_yyyyy[i] * fe_0 + tk_zz_xyyyyy[i] * pa_x[i] + 2.0 * ts_xzz_xyyyyy[i] * fz_0;

        tk_xzz_xyyyyz[i] = tk_zz_yyyyz[i] * fe_0 + tk_zz_xyyyyz[i] * pa_x[i] + 2.0 * ts_xzz_xyyyyz[i] * fz_0;

        tk_xzz_xyyyzz[i] = tk_zz_yyyzz[i] * fe_0 + tk_zz_xyyyzz[i] * pa_x[i] + 2.0 * ts_xzz_xyyyzz[i] * fz_0;

        tk_xzz_xyyzzz[i] = tk_zz_yyzzz[i] * fe_0 + tk_zz_xyyzzz[i] * pa_x[i] + 2.0 * ts_xzz_xyyzzz[i] * fz_0;

        tk_xzz_xyzzzz[i] = tk_zz_yzzzz[i] * fe_0 + tk_zz_xyzzzz[i] * pa_x[i] + 2.0 * ts_xzz_xyzzzz[i] * fz_0;

        tk_xzz_xzzzzz[i] = tk_zz_zzzzz[i] * fe_0 + tk_zz_xzzzzz[i] * pa_x[i] + 2.0 * ts_xzz_xzzzzz[i] * fz_0;

        tk_xzz_yyyyyy[i] = tk_zz_yyyyyy[i] * pa_x[i] + 2.0 * ts_xzz_yyyyyy[i] * fz_0;

        tk_xzz_yyyyyz[i] = tk_zz_yyyyyz[i] * pa_x[i] + 2.0 * ts_xzz_yyyyyz[i] * fz_0;

        tk_xzz_yyyyzz[i] = tk_zz_yyyyzz[i] * pa_x[i] + 2.0 * ts_xzz_yyyyzz[i] * fz_0;

        tk_xzz_yyyzzz[i] = tk_zz_yyyzzz[i] * pa_x[i] + 2.0 * ts_xzz_yyyzzz[i] * fz_0;

        tk_xzz_yyzzzz[i] = tk_zz_yyzzzz[i] * pa_x[i] + 2.0 * ts_xzz_yyzzzz[i] * fz_0;

        tk_xzz_yzzzzz[i] = tk_zz_yzzzzz[i] * pa_x[i] + 2.0 * ts_xzz_yzzzzz[i] * fz_0;

        tk_xzz_zzzzzz[i] = tk_zz_zzzzzz[i] * pa_x[i] + 2.0 * ts_xzz_zzzzzz[i] * fz_0;
    }

    // Set up 168-196 components of targeted buffer : FI

    auto tk_yyy_xxxxxx = pbuffer.data(idx_kin_fi + 168);

    auto tk_yyy_xxxxxy = pbuffer.data(idx_kin_fi + 169);

    auto tk_yyy_xxxxxz = pbuffer.data(idx_kin_fi + 170);

    auto tk_yyy_xxxxyy = pbuffer.data(idx_kin_fi + 171);

    auto tk_yyy_xxxxyz = pbuffer.data(idx_kin_fi + 172);

    auto tk_yyy_xxxxzz = pbuffer.data(idx_kin_fi + 173);

    auto tk_yyy_xxxyyy = pbuffer.data(idx_kin_fi + 174);

    auto tk_yyy_xxxyyz = pbuffer.data(idx_kin_fi + 175);

    auto tk_yyy_xxxyzz = pbuffer.data(idx_kin_fi + 176);

    auto tk_yyy_xxxzzz = pbuffer.data(idx_kin_fi + 177);

    auto tk_yyy_xxyyyy = pbuffer.data(idx_kin_fi + 178);

    auto tk_yyy_xxyyyz = pbuffer.data(idx_kin_fi + 179);

    auto tk_yyy_xxyyzz = pbuffer.data(idx_kin_fi + 180);

    auto tk_yyy_xxyzzz = pbuffer.data(idx_kin_fi + 181);

    auto tk_yyy_xxzzzz = pbuffer.data(idx_kin_fi + 182);

    auto tk_yyy_xyyyyy = pbuffer.data(idx_kin_fi + 183);

    auto tk_yyy_xyyyyz = pbuffer.data(idx_kin_fi + 184);

    auto tk_yyy_xyyyzz = pbuffer.data(idx_kin_fi + 185);

    auto tk_yyy_xyyzzz = pbuffer.data(idx_kin_fi + 186);

    auto tk_yyy_xyzzzz = pbuffer.data(idx_kin_fi + 187);

    auto tk_yyy_xzzzzz = pbuffer.data(idx_kin_fi + 188);

    auto tk_yyy_yyyyyy = pbuffer.data(idx_kin_fi + 189);

    auto tk_yyy_yyyyyz = pbuffer.data(idx_kin_fi + 190);

    auto tk_yyy_yyyyzz = pbuffer.data(idx_kin_fi + 191);

    auto tk_yyy_yyyzzz = pbuffer.data(idx_kin_fi + 192);

    auto tk_yyy_yyzzzz = pbuffer.data(idx_kin_fi + 193);

    auto tk_yyy_yzzzzz = pbuffer.data(idx_kin_fi + 194);

    auto tk_yyy_zzzzzz = pbuffer.data(idx_kin_fi + 195);

#pragma omp simd aligned(pa_y,              \
                             tk_y_xxxxxx,   \
                             tk_y_xxxxxy,   \
                             tk_y_xxxxxz,   \
                             tk_y_xxxxyy,   \
                             tk_y_xxxxyz,   \
                             tk_y_xxxxzz,   \
                             tk_y_xxxyyy,   \
                             tk_y_xxxyyz,   \
                             tk_y_xxxyzz,   \
                             tk_y_xxxzzz,   \
                             tk_y_xxyyyy,   \
                             tk_y_xxyyyz,   \
                             tk_y_xxyyzz,   \
                             tk_y_xxyzzz,   \
                             tk_y_xxzzzz,   \
                             tk_y_xyyyyy,   \
                             tk_y_xyyyyz,   \
                             tk_y_xyyyzz,   \
                             tk_y_xyyzzz,   \
                             tk_y_xyzzzz,   \
                             tk_y_xzzzzz,   \
                             tk_y_yyyyyy,   \
                             tk_y_yyyyyz,   \
                             tk_y_yyyyzz,   \
                             tk_y_yyyzzz,   \
                             tk_y_yyzzzz,   \
                             tk_y_yzzzzz,   \
                             tk_y_zzzzzz,   \
                             tk_yy_xxxxx,   \
                             tk_yy_xxxxxx,  \
                             tk_yy_xxxxxy,  \
                             tk_yy_xxxxxz,  \
                             tk_yy_xxxxy,   \
                             tk_yy_xxxxyy,  \
                             tk_yy_xxxxyz,  \
                             tk_yy_xxxxz,   \
                             tk_yy_xxxxzz,  \
                             tk_yy_xxxyy,   \
                             tk_yy_xxxyyy,  \
                             tk_yy_xxxyyz,  \
                             tk_yy_xxxyz,   \
                             tk_yy_xxxyzz,  \
                             tk_yy_xxxzz,   \
                             tk_yy_xxxzzz,  \
                             tk_yy_xxyyy,   \
                             tk_yy_xxyyyy,  \
                             tk_yy_xxyyyz,  \
                             tk_yy_xxyyz,   \
                             tk_yy_xxyyzz,  \
                             tk_yy_xxyzz,   \
                             tk_yy_xxyzzz,  \
                             tk_yy_xxzzz,   \
                             tk_yy_xxzzzz,  \
                             tk_yy_xyyyy,   \
                             tk_yy_xyyyyy,  \
                             tk_yy_xyyyyz,  \
                             tk_yy_xyyyz,   \
                             tk_yy_xyyyzz,  \
                             tk_yy_xyyzz,   \
                             tk_yy_xyyzzz,  \
                             tk_yy_xyzzz,   \
                             tk_yy_xyzzzz,  \
                             tk_yy_xzzzz,   \
                             tk_yy_xzzzzz,  \
                             tk_yy_yyyyy,   \
                             tk_yy_yyyyyy,  \
                             tk_yy_yyyyyz,  \
                             tk_yy_yyyyz,   \
                             tk_yy_yyyyzz,  \
                             tk_yy_yyyzz,   \
                             tk_yy_yyyzzz,  \
                             tk_yy_yyzzz,   \
                             tk_yy_yyzzzz,  \
                             tk_yy_yzzzz,   \
                             tk_yy_yzzzzz,  \
                             tk_yy_zzzzz,   \
                             tk_yy_zzzzzz,  \
                             tk_yyy_xxxxxx, \
                             tk_yyy_xxxxxy, \
                             tk_yyy_xxxxxz, \
                             tk_yyy_xxxxyy, \
                             tk_yyy_xxxxyz, \
                             tk_yyy_xxxxzz, \
                             tk_yyy_xxxyyy, \
                             tk_yyy_xxxyyz, \
                             tk_yyy_xxxyzz, \
                             tk_yyy_xxxzzz, \
                             tk_yyy_xxyyyy, \
                             tk_yyy_xxyyyz, \
                             tk_yyy_xxyyzz, \
                             tk_yyy_xxyzzz, \
                             tk_yyy_xxzzzz, \
                             tk_yyy_xyyyyy, \
                             tk_yyy_xyyyyz, \
                             tk_yyy_xyyyzz, \
                             tk_yyy_xyyzzz, \
                             tk_yyy_xyzzzz, \
                             tk_yyy_xzzzzz, \
                             tk_yyy_yyyyyy, \
                             tk_yyy_yyyyyz, \
                             tk_yyy_yyyyzz, \
                             tk_yyy_yyyzzz, \
                             tk_yyy_yyzzzz, \
                             tk_yyy_yzzzzz, \
                             tk_yyy_zzzzzz, \
                             ts_y_xxxxxx,   \
                             ts_y_xxxxxy,   \
                             ts_y_xxxxxz,   \
                             ts_y_xxxxyy,   \
                             ts_y_xxxxyz,   \
                             ts_y_xxxxzz,   \
                             ts_y_xxxyyy,   \
                             ts_y_xxxyyz,   \
                             ts_y_xxxyzz,   \
                             ts_y_xxxzzz,   \
                             ts_y_xxyyyy,   \
                             ts_y_xxyyyz,   \
                             ts_y_xxyyzz,   \
                             ts_y_xxyzzz,   \
                             ts_y_xxzzzz,   \
                             ts_y_xyyyyy,   \
                             ts_y_xyyyyz,   \
                             ts_y_xyyyzz,   \
                             ts_y_xyyzzz,   \
                             ts_y_xyzzzz,   \
                             ts_y_xzzzzz,   \
                             ts_y_yyyyyy,   \
                             ts_y_yyyyyz,   \
                             ts_y_yyyyzz,   \
                             ts_y_yyyzzz,   \
                             ts_y_yyzzzz,   \
                             ts_y_yzzzzz,   \
                             ts_y_zzzzzz,   \
                             ts_yyy_xxxxxx, \
                             ts_yyy_xxxxxy, \
                             ts_yyy_xxxxxz, \
                             ts_yyy_xxxxyy, \
                             ts_yyy_xxxxyz, \
                             ts_yyy_xxxxzz, \
                             ts_yyy_xxxyyy, \
                             ts_yyy_xxxyyz, \
                             ts_yyy_xxxyzz, \
                             ts_yyy_xxxzzz, \
                             ts_yyy_xxyyyy, \
                             ts_yyy_xxyyyz, \
                             ts_yyy_xxyyzz, \
                             ts_yyy_xxyzzz, \
                             ts_yyy_xxzzzz, \
                             ts_yyy_xyyyyy, \
                             ts_yyy_xyyyyz, \
                             ts_yyy_xyyyzz, \
                             ts_yyy_xyyzzz, \
                             ts_yyy_xyzzzz, \
                             ts_yyy_xzzzzz, \
                             ts_yyy_yyyyyy, \
                             ts_yyy_yyyyyz, \
                             ts_yyy_yyyyzz, \
                             ts_yyy_yyyzzz, \
                             ts_yyy_yyzzzz, \
                             ts_yyy_yzzzzz, \
                             ts_yyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyy_xxxxxx[i] =
            -4.0 * ts_y_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxxx[i] * fe_0 + tk_yy_xxxxxx[i] * pa_y[i] + 2.0 * ts_yyy_xxxxxx[i] * fz_0;

        tk_yyy_xxxxxy[i] = -4.0 * ts_y_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxxy[i] * fe_0 + tk_yy_xxxxx[i] * fe_0 + tk_yy_xxxxxy[i] * pa_y[i] +
                           2.0 * ts_yyy_xxxxxy[i] * fz_0;

        tk_yyy_xxxxxz[i] =
            -4.0 * ts_y_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxxz[i] * fe_0 + tk_yy_xxxxxz[i] * pa_y[i] + 2.0 * ts_yyy_xxxxxz[i] * fz_0;

        tk_yyy_xxxxyy[i] = -4.0 * ts_y_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxyy[i] * fe_0 + 2.0 * tk_yy_xxxxy[i] * fe_0 +
                           tk_yy_xxxxyy[i] * pa_y[i] + 2.0 * ts_yyy_xxxxyy[i] * fz_0;

        tk_yyy_xxxxyz[i] = -4.0 * ts_y_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxyz[i] * fe_0 + tk_yy_xxxxz[i] * fe_0 + tk_yy_xxxxyz[i] * pa_y[i] +
                           2.0 * ts_yyy_xxxxyz[i] * fz_0;

        tk_yyy_xxxxzz[i] =
            -4.0 * ts_y_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxxzz[i] * fe_0 + tk_yy_xxxxzz[i] * pa_y[i] + 2.0 * ts_yyy_xxxxzz[i] * fz_0;

        tk_yyy_xxxyyy[i] = -4.0 * ts_y_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxyyy[i] * fe_0 + 3.0 * tk_yy_xxxyy[i] * fe_0 +
                           tk_yy_xxxyyy[i] * pa_y[i] + 2.0 * ts_yyy_xxxyyy[i] * fz_0;

        tk_yyy_xxxyyz[i] = -4.0 * ts_y_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxyyz[i] * fe_0 + 2.0 * tk_yy_xxxyz[i] * fe_0 +
                           tk_yy_xxxyyz[i] * pa_y[i] + 2.0 * ts_yyy_xxxyyz[i] * fz_0;

        tk_yyy_xxxyzz[i] = -4.0 * ts_y_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxyzz[i] * fe_0 + tk_yy_xxxzz[i] * fe_0 + tk_yy_xxxyzz[i] * pa_y[i] +
                           2.0 * ts_yyy_xxxyzz[i] * fz_0;

        tk_yyy_xxxzzz[i] =
            -4.0 * ts_y_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxxzzz[i] * fe_0 + tk_yy_xxxzzz[i] * pa_y[i] + 2.0 * ts_yyy_xxxzzz[i] * fz_0;

        tk_yyy_xxyyyy[i] = -4.0 * ts_y_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyyyy[i] * fe_0 + 4.0 * tk_yy_xxyyy[i] * fe_0 +
                           tk_yy_xxyyyy[i] * pa_y[i] + 2.0 * ts_yyy_xxyyyy[i] * fz_0;

        tk_yyy_xxyyyz[i] = -4.0 * ts_y_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyyyz[i] * fe_0 + 3.0 * tk_yy_xxyyz[i] * fe_0 +
                           tk_yy_xxyyyz[i] * pa_y[i] + 2.0 * ts_yyy_xxyyyz[i] * fz_0;

        tk_yyy_xxyyzz[i] = -4.0 * ts_y_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyyzz[i] * fe_0 + 2.0 * tk_yy_xxyzz[i] * fe_0 +
                           tk_yy_xxyyzz[i] * pa_y[i] + 2.0 * ts_yyy_xxyyzz[i] * fz_0;

        tk_yyy_xxyzzz[i] = -4.0 * ts_y_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxyzzz[i] * fe_0 + tk_yy_xxzzz[i] * fe_0 + tk_yy_xxyzzz[i] * pa_y[i] +
                           2.0 * ts_yyy_xxyzzz[i] * fz_0;

        tk_yyy_xxzzzz[i] =
            -4.0 * ts_y_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xxzzzz[i] * fe_0 + tk_yy_xxzzzz[i] * pa_y[i] + 2.0 * ts_yyy_xxzzzz[i] * fz_0;

        tk_yyy_xyyyyy[i] = -4.0 * ts_y_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyyyy[i] * fe_0 + 5.0 * tk_yy_xyyyy[i] * fe_0 +
                           tk_yy_xyyyyy[i] * pa_y[i] + 2.0 * ts_yyy_xyyyyy[i] * fz_0;

        tk_yyy_xyyyyz[i] = -4.0 * ts_y_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyyyz[i] * fe_0 + 4.0 * tk_yy_xyyyz[i] * fe_0 +
                           tk_yy_xyyyyz[i] * pa_y[i] + 2.0 * ts_yyy_xyyyyz[i] * fz_0;

        tk_yyy_xyyyzz[i] = -4.0 * ts_y_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyyzz[i] * fe_0 + 3.0 * tk_yy_xyyzz[i] * fe_0 +
                           tk_yy_xyyyzz[i] * pa_y[i] + 2.0 * ts_yyy_xyyyzz[i] * fz_0;

        tk_yyy_xyyzzz[i] = -4.0 * ts_y_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyyzzz[i] * fe_0 + 2.0 * tk_yy_xyzzz[i] * fe_0 +
                           tk_yy_xyyzzz[i] * pa_y[i] + 2.0 * ts_yyy_xyyzzz[i] * fz_0;

        tk_yyy_xyzzzz[i] = -4.0 * ts_y_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xyzzzz[i] * fe_0 + tk_yy_xzzzz[i] * fe_0 + tk_yy_xyzzzz[i] * pa_y[i] +
                           2.0 * ts_yyy_xyzzzz[i] * fz_0;

        tk_yyy_xzzzzz[i] =
            -4.0 * ts_y_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_xzzzzz[i] * fe_0 + tk_yy_xzzzzz[i] * pa_y[i] + 2.0 * ts_yyy_xzzzzz[i] * fz_0;

        tk_yyy_yyyyyy[i] = -4.0 * ts_y_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyyyy[i] * fe_0 + 6.0 * tk_yy_yyyyy[i] * fe_0 +
                           tk_yy_yyyyyy[i] * pa_y[i] + 2.0 * ts_yyy_yyyyyy[i] * fz_0;

        tk_yyy_yyyyyz[i] = -4.0 * ts_y_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyyyz[i] * fe_0 + 5.0 * tk_yy_yyyyz[i] * fe_0 +
                           tk_yy_yyyyyz[i] * pa_y[i] + 2.0 * ts_yyy_yyyyyz[i] * fz_0;

        tk_yyy_yyyyzz[i] = -4.0 * ts_y_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyyzz[i] * fe_0 + 4.0 * tk_yy_yyyzz[i] * fe_0 +
                           tk_yy_yyyyzz[i] * pa_y[i] + 2.0 * ts_yyy_yyyyzz[i] * fz_0;

        tk_yyy_yyyzzz[i] = -4.0 * ts_y_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyyzzz[i] * fe_0 + 3.0 * tk_yy_yyzzz[i] * fe_0 +
                           tk_yy_yyyzzz[i] * pa_y[i] + 2.0 * ts_yyy_yyyzzz[i] * fz_0;

        tk_yyy_yyzzzz[i] = -4.0 * ts_y_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yyzzzz[i] * fe_0 + 2.0 * tk_yy_yzzzz[i] * fe_0 +
                           tk_yy_yyzzzz[i] * pa_y[i] + 2.0 * ts_yyy_yyzzzz[i] * fz_0;

        tk_yyy_yzzzzz[i] = -4.0 * ts_y_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_yzzzzz[i] * fe_0 + tk_yy_zzzzz[i] * fe_0 + tk_yy_yzzzzz[i] * pa_y[i] +
                           2.0 * ts_yyy_yzzzzz[i] * fz_0;

        tk_yyy_zzzzzz[i] =
            -4.0 * ts_y_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_y_zzzzzz[i] * fe_0 + tk_yy_zzzzzz[i] * pa_y[i] + 2.0 * ts_yyy_zzzzzz[i] * fz_0;
    }

    // Set up 196-224 components of targeted buffer : FI

    auto tk_yyz_xxxxxx = pbuffer.data(idx_kin_fi + 196);

    auto tk_yyz_xxxxxy = pbuffer.data(idx_kin_fi + 197);

    auto tk_yyz_xxxxxz = pbuffer.data(idx_kin_fi + 198);

    auto tk_yyz_xxxxyy = pbuffer.data(idx_kin_fi + 199);

    auto tk_yyz_xxxxyz = pbuffer.data(idx_kin_fi + 200);

    auto tk_yyz_xxxxzz = pbuffer.data(idx_kin_fi + 201);

    auto tk_yyz_xxxyyy = pbuffer.data(idx_kin_fi + 202);

    auto tk_yyz_xxxyyz = pbuffer.data(idx_kin_fi + 203);

    auto tk_yyz_xxxyzz = pbuffer.data(idx_kin_fi + 204);

    auto tk_yyz_xxxzzz = pbuffer.data(idx_kin_fi + 205);

    auto tk_yyz_xxyyyy = pbuffer.data(idx_kin_fi + 206);

    auto tk_yyz_xxyyyz = pbuffer.data(idx_kin_fi + 207);

    auto tk_yyz_xxyyzz = pbuffer.data(idx_kin_fi + 208);

    auto tk_yyz_xxyzzz = pbuffer.data(idx_kin_fi + 209);

    auto tk_yyz_xxzzzz = pbuffer.data(idx_kin_fi + 210);

    auto tk_yyz_xyyyyy = pbuffer.data(idx_kin_fi + 211);

    auto tk_yyz_xyyyyz = pbuffer.data(idx_kin_fi + 212);

    auto tk_yyz_xyyyzz = pbuffer.data(idx_kin_fi + 213);

    auto tk_yyz_xyyzzz = pbuffer.data(idx_kin_fi + 214);

    auto tk_yyz_xyzzzz = pbuffer.data(idx_kin_fi + 215);

    auto tk_yyz_xzzzzz = pbuffer.data(idx_kin_fi + 216);

    auto tk_yyz_yyyyyy = pbuffer.data(idx_kin_fi + 217);

    auto tk_yyz_yyyyyz = pbuffer.data(idx_kin_fi + 218);

    auto tk_yyz_yyyyzz = pbuffer.data(idx_kin_fi + 219);

    auto tk_yyz_yyyzzz = pbuffer.data(idx_kin_fi + 220);

    auto tk_yyz_yyzzzz = pbuffer.data(idx_kin_fi + 221);

    auto tk_yyz_yzzzzz = pbuffer.data(idx_kin_fi + 222);

    auto tk_yyz_zzzzzz = pbuffer.data(idx_kin_fi + 223);

#pragma omp simd aligned(pa_z,              \
                             tk_yy_xxxxx,   \
                             tk_yy_xxxxxx,  \
                             tk_yy_xxxxxy,  \
                             tk_yy_xxxxxz,  \
                             tk_yy_xxxxy,   \
                             tk_yy_xxxxyy,  \
                             tk_yy_xxxxyz,  \
                             tk_yy_xxxxz,   \
                             tk_yy_xxxxzz,  \
                             tk_yy_xxxyy,   \
                             tk_yy_xxxyyy,  \
                             tk_yy_xxxyyz,  \
                             tk_yy_xxxyz,   \
                             tk_yy_xxxyzz,  \
                             tk_yy_xxxzz,   \
                             tk_yy_xxxzzz,  \
                             tk_yy_xxyyy,   \
                             tk_yy_xxyyyy,  \
                             tk_yy_xxyyyz,  \
                             tk_yy_xxyyz,   \
                             tk_yy_xxyyzz,  \
                             tk_yy_xxyzz,   \
                             tk_yy_xxyzzz,  \
                             tk_yy_xxzzz,   \
                             tk_yy_xxzzzz,  \
                             tk_yy_xyyyy,   \
                             tk_yy_xyyyyy,  \
                             tk_yy_xyyyyz,  \
                             tk_yy_xyyyz,   \
                             tk_yy_xyyyzz,  \
                             tk_yy_xyyzz,   \
                             tk_yy_xyyzzz,  \
                             tk_yy_xyzzz,   \
                             tk_yy_xyzzzz,  \
                             tk_yy_xzzzz,   \
                             tk_yy_xzzzzz,  \
                             tk_yy_yyyyy,   \
                             tk_yy_yyyyyy,  \
                             tk_yy_yyyyyz,  \
                             tk_yy_yyyyz,   \
                             tk_yy_yyyyzz,  \
                             tk_yy_yyyzz,   \
                             tk_yy_yyyzzz,  \
                             tk_yy_yyzzz,   \
                             tk_yy_yyzzzz,  \
                             tk_yy_yzzzz,   \
                             tk_yy_yzzzzz,  \
                             tk_yy_zzzzz,   \
                             tk_yy_zzzzzz,  \
                             tk_yyz_xxxxxx, \
                             tk_yyz_xxxxxy, \
                             tk_yyz_xxxxxz, \
                             tk_yyz_xxxxyy, \
                             tk_yyz_xxxxyz, \
                             tk_yyz_xxxxzz, \
                             tk_yyz_xxxyyy, \
                             tk_yyz_xxxyyz, \
                             tk_yyz_xxxyzz, \
                             tk_yyz_xxxzzz, \
                             tk_yyz_xxyyyy, \
                             tk_yyz_xxyyyz, \
                             tk_yyz_xxyyzz, \
                             tk_yyz_xxyzzz, \
                             tk_yyz_xxzzzz, \
                             tk_yyz_xyyyyy, \
                             tk_yyz_xyyyyz, \
                             tk_yyz_xyyyzz, \
                             tk_yyz_xyyzzz, \
                             tk_yyz_xyzzzz, \
                             tk_yyz_xzzzzz, \
                             tk_yyz_yyyyyy, \
                             tk_yyz_yyyyyz, \
                             tk_yyz_yyyyzz, \
                             tk_yyz_yyyzzz, \
                             tk_yyz_yyzzzz, \
                             tk_yyz_yzzzzz, \
                             tk_yyz_zzzzzz, \
                             ts_yyz_xxxxxx, \
                             ts_yyz_xxxxxy, \
                             ts_yyz_xxxxxz, \
                             ts_yyz_xxxxyy, \
                             ts_yyz_xxxxyz, \
                             ts_yyz_xxxxzz, \
                             ts_yyz_xxxyyy, \
                             ts_yyz_xxxyyz, \
                             ts_yyz_xxxyzz, \
                             ts_yyz_xxxzzz, \
                             ts_yyz_xxyyyy, \
                             ts_yyz_xxyyyz, \
                             ts_yyz_xxyyzz, \
                             ts_yyz_xxyzzz, \
                             ts_yyz_xxzzzz, \
                             ts_yyz_xyyyyy, \
                             ts_yyz_xyyyyz, \
                             ts_yyz_xyyyzz, \
                             ts_yyz_xyyzzz, \
                             ts_yyz_xyzzzz, \
                             ts_yyz_xzzzzz, \
                             ts_yyz_yyyyyy, \
                             ts_yyz_yyyyyz, \
                             ts_yyz_yyyyzz, \
                             ts_yyz_yyyzzz, \
                             ts_yyz_yyzzzz, \
                             ts_yyz_yzzzzz, \
                             ts_yyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyz_xxxxxx[i] = tk_yy_xxxxxx[i] * pa_z[i] + 2.0 * ts_yyz_xxxxxx[i] * fz_0;

        tk_yyz_xxxxxy[i] = tk_yy_xxxxxy[i] * pa_z[i] + 2.0 * ts_yyz_xxxxxy[i] * fz_0;

        tk_yyz_xxxxxz[i] = tk_yy_xxxxx[i] * fe_0 + tk_yy_xxxxxz[i] * pa_z[i] + 2.0 * ts_yyz_xxxxxz[i] * fz_0;

        tk_yyz_xxxxyy[i] = tk_yy_xxxxyy[i] * pa_z[i] + 2.0 * ts_yyz_xxxxyy[i] * fz_0;

        tk_yyz_xxxxyz[i] = tk_yy_xxxxy[i] * fe_0 + tk_yy_xxxxyz[i] * pa_z[i] + 2.0 * ts_yyz_xxxxyz[i] * fz_0;

        tk_yyz_xxxxzz[i] = 2.0 * tk_yy_xxxxz[i] * fe_0 + tk_yy_xxxxzz[i] * pa_z[i] + 2.0 * ts_yyz_xxxxzz[i] * fz_0;

        tk_yyz_xxxyyy[i] = tk_yy_xxxyyy[i] * pa_z[i] + 2.0 * ts_yyz_xxxyyy[i] * fz_0;

        tk_yyz_xxxyyz[i] = tk_yy_xxxyy[i] * fe_0 + tk_yy_xxxyyz[i] * pa_z[i] + 2.0 * ts_yyz_xxxyyz[i] * fz_0;

        tk_yyz_xxxyzz[i] = 2.0 * tk_yy_xxxyz[i] * fe_0 + tk_yy_xxxyzz[i] * pa_z[i] + 2.0 * ts_yyz_xxxyzz[i] * fz_0;

        tk_yyz_xxxzzz[i] = 3.0 * tk_yy_xxxzz[i] * fe_0 + tk_yy_xxxzzz[i] * pa_z[i] + 2.0 * ts_yyz_xxxzzz[i] * fz_0;

        tk_yyz_xxyyyy[i] = tk_yy_xxyyyy[i] * pa_z[i] + 2.0 * ts_yyz_xxyyyy[i] * fz_0;

        tk_yyz_xxyyyz[i] = tk_yy_xxyyy[i] * fe_0 + tk_yy_xxyyyz[i] * pa_z[i] + 2.0 * ts_yyz_xxyyyz[i] * fz_0;

        tk_yyz_xxyyzz[i] = 2.0 * tk_yy_xxyyz[i] * fe_0 + tk_yy_xxyyzz[i] * pa_z[i] + 2.0 * ts_yyz_xxyyzz[i] * fz_0;

        tk_yyz_xxyzzz[i] = 3.0 * tk_yy_xxyzz[i] * fe_0 + tk_yy_xxyzzz[i] * pa_z[i] + 2.0 * ts_yyz_xxyzzz[i] * fz_0;

        tk_yyz_xxzzzz[i] = 4.0 * tk_yy_xxzzz[i] * fe_0 + tk_yy_xxzzzz[i] * pa_z[i] + 2.0 * ts_yyz_xxzzzz[i] * fz_0;

        tk_yyz_xyyyyy[i] = tk_yy_xyyyyy[i] * pa_z[i] + 2.0 * ts_yyz_xyyyyy[i] * fz_0;

        tk_yyz_xyyyyz[i] = tk_yy_xyyyy[i] * fe_0 + tk_yy_xyyyyz[i] * pa_z[i] + 2.0 * ts_yyz_xyyyyz[i] * fz_0;

        tk_yyz_xyyyzz[i] = 2.0 * tk_yy_xyyyz[i] * fe_0 + tk_yy_xyyyzz[i] * pa_z[i] + 2.0 * ts_yyz_xyyyzz[i] * fz_0;

        tk_yyz_xyyzzz[i] = 3.0 * tk_yy_xyyzz[i] * fe_0 + tk_yy_xyyzzz[i] * pa_z[i] + 2.0 * ts_yyz_xyyzzz[i] * fz_0;

        tk_yyz_xyzzzz[i] = 4.0 * tk_yy_xyzzz[i] * fe_0 + tk_yy_xyzzzz[i] * pa_z[i] + 2.0 * ts_yyz_xyzzzz[i] * fz_0;

        tk_yyz_xzzzzz[i] = 5.0 * tk_yy_xzzzz[i] * fe_0 + tk_yy_xzzzzz[i] * pa_z[i] + 2.0 * ts_yyz_xzzzzz[i] * fz_0;

        tk_yyz_yyyyyy[i] = tk_yy_yyyyyy[i] * pa_z[i] + 2.0 * ts_yyz_yyyyyy[i] * fz_0;

        tk_yyz_yyyyyz[i] = tk_yy_yyyyy[i] * fe_0 + tk_yy_yyyyyz[i] * pa_z[i] + 2.0 * ts_yyz_yyyyyz[i] * fz_0;

        tk_yyz_yyyyzz[i] = 2.0 * tk_yy_yyyyz[i] * fe_0 + tk_yy_yyyyzz[i] * pa_z[i] + 2.0 * ts_yyz_yyyyzz[i] * fz_0;

        tk_yyz_yyyzzz[i] = 3.0 * tk_yy_yyyzz[i] * fe_0 + tk_yy_yyyzzz[i] * pa_z[i] + 2.0 * ts_yyz_yyyzzz[i] * fz_0;

        tk_yyz_yyzzzz[i] = 4.0 * tk_yy_yyzzz[i] * fe_0 + tk_yy_yyzzzz[i] * pa_z[i] + 2.0 * ts_yyz_yyzzzz[i] * fz_0;

        tk_yyz_yzzzzz[i] = 5.0 * tk_yy_yzzzz[i] * fe_0 + tk_yy_yzzzzz[i] * pa_z[i] + 2.0 * ts_yyz_yzzzzz[i] * fz_0;

        tk_yyz_zzzzzz[i] = 6.0 * tk_yy_zzzzz[i] * fe_0 + tk_yy_zzzzzz[i] * pa_z[i] + 2.0 * ts_yyz_zzzzzz[i] * fz_0;
    }

    // Set up 224-252 components of targeted buffer : FI

    auto tk_yzz_xxxxxx = pbuffer.data(idx_kin_fi + 224);

    auto tk_yzz_xxxxxy = pbuffer.data(idx_kin_fi + 225);

    auto tk_yzz_xxxxxz = pbuffer.data(idx_kin_fi + 226);

    auto tk_yzz_xxxxyy = pbuffer.data(idx_kin_fi + 227);

    auto tk_yzz_xxxxyz = pbuffer.data(idx_kin_fi + 228);

    auto tk_yzz_xxxxzz = pbuffer.data(idx_kin_fi + 229);

    auto tk_yzz_xxxyyy = pbuffer.data(idx_kin_fi + 230);

    auto tk_yzz_xxxyyz = pbuffer.data(idx_kin_fi + 231);

    auto tk_yzz_xxxyzz = pbuffer.data(idx_kin_fi + 232);

    auto tk_yzz_xxxzzz = pbuffer.data(idx_kin_fi + 233);

    auto tk_yzz_xxyyyy = pbuffer.data(idx_kin_fi + 234);

    auto tk_yzz_xxyyyz = pbuffer.data(idx_kin_fi + 235);

    auto tk_yzz_xxyyzz = pbuffer.data(idx_kin_fi + 236);

    auto tk_yzz_xxyzzz = pbuffer.data(idx_kin_fi + 237);

    auto tk_yzz_xxzzzz = pbuffer.data(idx_kin_fi + 238);

    auto tk_yzz_xyyyyy = pbuffer.data(idx_kin_fi + 239);

    auto tk_yzz_xyyyyz = pbuffer.data(idx_kin_fi + 240);

    auto tk_yzz_xyyyzz = pbuffer.data(idx_kin_fi + 241);

    auto tk_yzz_xyyzzz = pbuffer.data(idx_kin_fi + 242);

    auto tk_yzz_xyzzzz = pbuffer.data(idx_kin_fi + 243);

    auto tk_yzz_xzzzzz = pbuffer.data(idx_kin_fi + 244);

    auto tk_yzz_yyyyyy = pbuffer.data(idx_kin_fi + 245);

    auto tk_yzz_yyyyyz = pbuffer.data(idx_kin_fi + 246);

    auto tk_yzz_yyyyzz = pbuffer.data(idx_kin_fi + 247);

    auto tk_yzz_yyyzzz = pbuffer.data(idx_kin_fi + 248);

    auto tk_yzz_yyzzzz = pbuffer.data(idx_kin_fi + 249);

    auto tk_yzz_yzzzzz = pbuffer.data(idx_kin_fi + 250);

    auto tk_yzz_zzzzzz = pbuffer.data(idx_kin_fi + 251);

#pragma omp simd aligned(pa_y,              \
                             tk_yzz_xxxxxx, \
                             tk_yzz_xxxxxy, \
                             tk_yzz_xxxxxz, \
                             tk_yzz_xxxxyy, \
                             tk_yzz_xxxxyz, \
                             tk_yzz_xxxxzz, \
                             tk_yzz_xxxyyy, \
                             tk_yzz_xxxyyz, \
                             tk_yzz_xxxyzz, \
                             tk_yzz_xxxzzz, \
                             tk_yzz_xxyyyy, \
                             tk_yzz_xxyyyz, \
                             tk_yzz_xxyyzz, \
                             tk_yzz_xxyzzz, \
                             tk_yzz_xxzzzz, \
                             tk_yzz_xyyyyy, \
                             tk_yzz_xyyyyz, \
                             tk_yzz_xyyyzz, \
                             tk_yzz_xyyzzz, \
                             tk_yzz_xyzzzz, \
                             tk_yzz_xzzzzz, \
                             tk_yzz_yyyyyy, \
                             tk_yzz_yyyyyz, \
                             tk_yzz_yyyyzz, \
                             tk_yzz_yyyzzz, \
                             tk_yzz_yyzzzz, \
                             tk_yzz_yzzzzz, \
                             tk_yzz_zzzzzz, \
                             tk_zz_xxxxx,   \
                             tk_zz_xxxxxx,  \
                             tk_zz_xxxxxy,  \
                             tk_zz_xxxxxz,  \
                             tk_zz_xxxxy,   \
                             tk_zz_xxxxyy,  \
                             tk_zz_xxxxyz,  \
                             tk_zz_xxxxz,   \
                             tk_zz_xxxxzz,  \
                             tk_zz_xxxyy,   \
                             tk_zz_xxxyyy,  \
                             tk_zz_xxxyyz,  \
                             tk_zz_xxxyz,   \
                             tk_zz_xxxyzz,  \
                             tk_zz_xxxzz,   \
                             tk_zz_xxxzzz,  \
                             tk_zz_xxyyy,   \
                             tk_zz_xxyyyy,  \
                             tk_zz_xxyyyz,  \
                             tk_zz_xxyyz,   \
                             tk_zz_xxyyzz,  \
                             tk_zz_xxyzz,   \
                             tk_zz_xxyzzz,  \
                             tk_zz_xxzzz,   \
                             tk_zz_xxzzzz,  \
                             tk_zz_xyyyy,   \
                             tk_zz_xyyyyy,  \
                             tk_zz_xyyyyz,  \
                             tk_zz_xyyyz,   \
                             tk_zz_xyyyzz,  \
                             tk_zz_xyyzz,   \
                             tk_zz_xyyzzz,  \
                             tk_zz_xyzzz,   \
                             tk_zz_xyzzzz,  \
                             tk_zz_xzzzz,   \
                             tk_zz_xzzzzz,  \
                             tk_zz_yyyyy,   \
                             tk_zz_yyyyyy,  \
                             tk_zz_yyyyyz,  \
                             tk_zz_yyyyz,   \
                             tk_zz_yyyyzz,  \
                             tk_zz_yyyzz,   \
                             tk_zz_yyyzzz,  \
                             tk_zz_yyzzz,   \
                             tk_zz_yyzzzz,  \
                             tk_zz_yzzzz,   \
                             tk_zz_yzzzzz,  \
                             tk_zz_zzzzz,   \
                             tk_zz_zzzzzz,  \
                             ts_yzz_xxxxxx, \
                             ts_yzz_xxxxxy, \
                             ts_yzz_xxxxxz, \
                             ts_yzz_xxxxyy, \
                             ts_yzz_xxxxyz, \
                             ts_yzz_xxxxzz, \
                             ts_yzz_xxxyyy, \
                             ts_yzz_xxxyyz, \
                             ts_yzz_xxxyzz, \
                             ts_yzz_xxxzzz, \
                             ts_yzz_xxyyyy, \
                             ts_yzz_xxyyyz, \
                             ts_yzz_xxyyzz, \
                             ts_yzz_xxyzzz, \
                             ts_yzz_xxzzzz, \
                             ts_yzz_xyyyyy, \
                             ts_yzz_xyyyyz, \
                             ts_yzz_xyyyzz, \
                             ts_yzz_xyyzzz, \
                             ts_yzz_xyzzzz, \
                             ts_yzz_xzzzzz, \
                             ts_yzz_yyyyyy, \
                             ts_yzz_yyyyyz, \
                             ts_yzz_yyyyzz, \
                             ts_yzz_yyyzzz, \
                             ts_yzz_yyzzzz, \
                             ts_yzz_yzzzzz, \
                             ts_yzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzz_xxxxxx[i] = tk_zz_xxxxxx[i] * pa_y[i] + 2.0 * ts_yzz_xxxxxx[i] * fz_0;

        tk_yzz_xxxxxy[i] = tk_zz_xxxxx[i] * fe_0 + tk_zz_xxxxxy[i] * pa_y[i] + 2.0 * ts_yzz_xxxxxy[i] * fz_0;

        tk_yzz_xxxxxz[i] = tk_zz_xxxxxz[i] * pa_y[i] + 2.0 * ts_yzz_xxxxxz[i] * fz_0;

        tk_yzz_xxxxyy[i] = 2.0 * tk_zz_xxxxy[i] * fe_0 + tk_zz_xxxxyy[i] * pa_y[i] + 2.0 * ts_yzz_xxxxyy[i] * fz_0;

        tk_yzz_xxxxyz[i] = tk_zz_xxxxz[i] * fe_0 + tk_zz_xxxxyz[i] * pa_y[i] + 2.0 * ts_yzz_xxxxyz[i] * fz_0;

        tk_yzz_xxxxzz[i] = tk_zz_xxxxzz[i] * pa_y[i] + 2.0 * ts_yzz_xxxxzz[i] * fz_0;

        tk_yzz_xxxyyy[i] = 3.0 * tk_zz_xxxyy[i] * fe_0 + tk_zz_xxxyyy[i] * pa_y[i] + 2.0 * ts_yzz_xxxyyy[i] * fz_0;

        tk_yzz_xxxyyz[i] = 2.0 * tk_zz_xxxyz[i] * fe_0 + tk_zz_xxxyyz[i] * pa_y[i] + 2.0 * ts_yzz_xxxyyz[i] * fz_0;

        tk_yzz_xxxyzz[i] = tk_zz_xxxzz[i] * fe_0 + tk_zz_xxxyzz[i] * pa_y[i] + 2.0 * ts_yzz_xxxyzz[i] * fz_0;

        tk_yzz_xxxzzz[i] = tk_zz_xxxzzz[i] * pa_y[i] + 2.0 * ts_yzz_xxxzzz[i] * fz_0;

        tk_yzz_xxyyyy[i] = 4.0 * tk_zz_xxyyy[i] * fe_0 + tk_zz_xxyyyy[i] * pa_y[i] + 2.0 * ts_yzz_xxyyyy[i] * fz_0;

        tk_yzz_xxyyyz[i] = 3.0 * tk_zz_xxyyz[i] * fe_0 + tk_zz_xxyyyz[i] * pa_y[i] + 2.0 * ts_yzz_xxyyyz[i] * fz_0;

        tk_yzz_xxyyzz[i] = 2.0 * tk_zz_xxyzz[i] * fe_0 + tk_zz_xxyyzz[i] * pa_y[i] + 2.0 * ts_yzz_xxyyzz[i] * fz_0;

        tk_yzz_xxyzzz[i] = tk_zz_xxzzz[i] * fe_0 + tk_zz_xxyzzz[i] * pa_y[i] + 2.0 * ts_yzz_xxyzzz[i] * fz_0;

        tk_yzz_xxzzzz[i] = tk_zz_xxzzzz[i] * pa_y[i] + 2.0 * ts_yzz_xxzzzz[i] * fz_0;

        tk_yzz_xyyyyy[i] = 5.0 * tk_zz_xyyyy[i] * fe_0 + tk_zz_xyyyyy[i] * pa_y[i] + 2.0 * ts_yzz_xyyyyy[i] * fz_0;

        tk_yzz_xyyyyz[i] = 4.0 * tk_zz_xyyyz[i] * fe_0 + tk_zz_xyyyyz[i] * pa_y[i] + 2.0 * ts_yzz_xyyyyz[i] * fz_0;

        tk_yzz_xyyyzz[i] = 3.0 * tk_zz_xyyzz[i] * fe_0 + tk_zz_xyyyzz[i] * pa_y[i] + 2.0 * ts_yzz_xyyyzz[i] * fz_0;

        tk_yzz_xyyzzz[i] = 2.0 * tk_zz_xyzzz[i] * fe_0 + tk_zz_xyyzzz[i] * pa_y[i] + 2.0 * ts_yzz_xyyzzz[i] * fz_0;

        tk_yzz_xyzzzz[i] = tk_zz_xzzzz[i] * fe_0 + tk_zz_xyzzzz[i] * pa_y[i] + 2.0 * ts_yzz_xyzzzz[i] * fz_0;

        tk_yzz_xzzzzz[i] = tk_zz_xzzzzz[i] * pa_y[i] + 2.0 * ts_yzz_xzzzzz[i] * fz_0;

        tk_yzz_yyyyyy[i] = 6.0 * tk_zz_yyyyy[i] * fe_0 + tk_zz_yyyyyy[i] * pa_y[i] + 2.0 * ts_yzz_yyyyyy[i] * fz_0;

        tk_yzz_yyyyyz[i] = 5.0 * tk_zz_yyyyz[i] * fe_0 + tk_zz_yyyyyz[i] * pa_y[i] + 2.0 * ts_yzz_yyyyyz[i] * fz_0;

        tk_yzz_yyyyzz[i] = 4.0 * tk_zz_yyyzz[i] * fe_0 + tk_zz_yyyyzz[i] * pa_y[i] + 2.0 * ts_yzz_yyyyzz[i] * fz_0;

        tk_yzz_yyyzzz[i] = 3.0 * tk_zz_yyzzz[i] * fe_0 + tk_zz_yyyzzz[i] * pa_y[i] + 2.0 * ts_yzz_yyyzzz[i] * fz_0;

        tk_yzz_yyzzzz[i] = 2.0 * tk_zz_yzzzz[i] * fe_0 + tk_zz_yyzzzz[i] * pa_y[i] + 2.0 * ts_yzz_yyzzzz[i] * fz_0;

        tk_yzz_yzzzzz[i] = tk_zz_zzzzz[i] * fe_0 + tk_zz_yzzzzz[i] * pa_y[i] + 2.0 * ts_yzz_yzzzzz[i] * fz_0;

        tk_yzz_zzzzzz[i] = tk_zz_zzzzzz[i] * pa_y[i] + 2.0 * ts_yzz_zzzzzz[i] * fz_0;
    }

    // Set up 252-280 components of targeted buffer : FI

    auto tk_zzz_xxxxxx = pbuffer.data(idx_kin_fi + 252);

    auto tk_zzz_xxxxxy = pbuffer.data(idx_kin_fi + 253);

    auto tk_zzz_xxxxxz = pbuffer.data(idx_kin_fi + 254);

    auto tk_zzz_xxxxyy = pbuffer.data(idx_kin_fi + 255);

    auto tk_zzz_xxxxyz = pbuffer.data(idx_kin_fi + 256);

    auto tk_zzz_xxxxzz = pbuffer.data(idx_kin_fi + 257);

    auto tk_zzz_xxxyyy = pbuffer.data(idx_kin_fi + 258);

    auto tk_zzz_xxxyyz = pbuffer.data(idx_kin_fi + 259);

    auto tk_zzz_xxxyzz = pbuffer.data(idx_kin_fi + 260);

    auto tk_zzz_xxxzzz = pbuffer.data(idx_kin_fi + 261);

    auto tk_zzz_xxyyyy = pbuffer.data(idx_kin_fi + 262);

    auto tk_zzz_xxyyyz = pbuffer.data(idx_kin_fi + 263);

    auto tk_zzz_xxyyzz = pbuffer.data(idx_kin_fi + 264);

    auto tk_zzz_xxyzzz = pbuffer.data(idx_kin_fi + 265);

    auto tk_zzz_xxzzzz = pbuffer.data(idx_kin_fi + 266);

    auto tk_zzz_xyyyyy = pbuffer.data(idx_kin_fi + 267);

    auto tk_zzz_xyyyyz = pbuffer.data(idx_kin_fi + 268);

    auto tk_zzz_xyyyzz = pbuffer.data(idx_kin_fi + 269);

    auto tk_zzz_xyyzzz = pbuffer.data(idx_kin_fi + 270);

    auto tk_zzz_xyzzzz = pbuffer.data(idx_kin_fi + 271);

    auto tk_zzz_xzzzzz = pbuffer.data(idx_kin_fi + 272);

    auto tk_zzz_yyyyyy = pbuffer.data(idx_kin_fi + 273);

    auto tk_zzz_yyyyyz = pbuffer.data(idx_kin_fi + 274);

    auto tk_zzz_yyyyzz = pbuffer.data(idx_kin_fi + 275);

    auto tk_zzz_yyyzzz = pbuffer.data(idx_kin_fi + 276);

    auto tk_zzz_yyzzzz = pbuffer.data(idx_kin_fi + 277);

    auto tk_zzz_yzzzzz = pbuffer.data(idx_kin_fi + 278);

    auto tk_zzz_zzzzzz = pbuffer.data(idx_kin_fi + 279);

#pragma omp simd aligned(pa_z,              \
                             tk_z_xxxxxx,   \
                             tk_z_xxxxxy,   \
                             tk_z_xxxxxz,   \
                             tk_z_xxxxyy,   \
                             tk_z_xxxxyz,   \
                             tk_z_xxxxzz,   \
                             tk_z_xxxyyy,   \
                             tk_z_xxxyyz,   \
                             tk_z_xxxyzz,   \
                             tk_z_xxxzzz,   \
                             tk_z_xxyyyy,   \
                             tk_z_xxyyyz,   \
                             tk_z_xxyyzz,   \
                             tk_z_xxyzzz,   \
                             tk_z_xxzzzz,   \
                             tk_z_xyyyyy,   \
                             tk_z_xyyyyz,   \
                             tk_z_xyyyzz,   \
                             tk_z_xyyzzz,   \
                             tk_z_xyzzzz,   \
                             tk_z_xzzzzz,   \
                             tk_z_yyyyyy,   \
                             tk_z_yyyyyz,   \
                             tk_z_yyyyzz,   \
                             tk_z_yyyzzz,   \
                             tk_z_yyzzzz,   \
                             tk_z_yzzzzz,   \
                             tk_z_zzzzzz,   \
                             tk_zz_xxxxx,   \
                             tk_zz_xxxxxx,  \
                             tk_zz_xxxxxy,  \
                             tk_zz_xxxxxz,  \
                             tk_zz_xxxxy,   \
                             tk_zz_xxxxyy,  \
                             tk_zz_xxxxyz,  \
                             tk_zz_xxxxz,   \
                             tk_zz_xxxxzz,  \
                             tk_zz_xxxyy,   \
                             tk_zz_xxxyyy,  \
                             tk_zz_xxxyyz,  \
                             tk_zz_xxxyz,   \
                             tk_zz_xxxyzz,  \
                             tk_zz_xxxzz,   \
                             tk_zz_xxxzzz,  \
                             tk_zz_xxyyy,   \
                             tk_zz_xxyyyy,  \
                             tk_zz_xxyyyz,  \
                             tk_zz_xxyyz,   \
                             tk_zz_xxyyzz,  \
                             tk_zz_xxyzz,   \
                             tk_zz_xxyzzz,  \
                             tk_zz_xxzzz,   \
                             tk_zz_xxzzzz,  \
                             tk_zz_xyyyy,   \
                             tk_zz_xyyyyy,  \
                             tk_zz_xyyyyz,  \
                             tk_zz_xyyyz,   \
                             tk_zz_xyyyzz,  \
                             tk_zz_xyyzz,   \
                             tk_zz_xyyzzz,  \
                             tk_zz_xyzzz,   \
                             tk_zz_xyzzzz,  \
                             tk_zz_xzzzz,   \
                             tk_zz_xzzzzz,  \
                             tk_zz_yyyyy,   \
                             tk_zz_yyyyyy,  \
                             tk_zz_yyyyyz,  \
                             tk_zz_yyyyz,   \
                             tk_zz_yyyyzz,  \
                             tk_zz_yyyzz,   \
                             tk_zz_yyyzzz,  \
                             tk_zz_yyzzz,   \
                             tk_zz_yyzzzz,  \
                             tk_zz_yzzzz,   \
                             tk_zz_yzzzzz,  \
                             tk_zz_zzzzz,   \
                             tk_zz_zzzzzz,  \
                             tk_zzz_xxxxxx, \
                             tk_zzz_xxxxxy, \
                             tk_zzz_xxxxxz, \
                             tk_zzz_xxxxyy, \
                             tk_zzz_xxxxyz, \
                             tk_zzz_xxxxzz, \
                             tk_zzz_xxxyyy, \
                             tk_zzz_xxxyyz, \
                             tk_zzz_xxxyzz, \
                             tk_zzz_xxxzzz, \
                             tk_zzz_xxyyyy, \
                             tk_zzz_xxyyyz, \
                             tk_zzz_xxyyzz, \
                             tk_zzz_xxyzzz, \
                             tk_zzz_xxzzzz, \
                             tk_zzz_xyyyyy, \
                             tk_zzz_xyyyyz, \
                             tk_zzz_xyyyzz, \
                             tk_zzz_xyyzzz, \
                             tk_zzz_xyzzzz, \
                             tk_zzz_xzzzzz, \
                             tk_zzz_yyyyyy, \
                             tk_zzz_yyyyyz, \
                             tk_zzz_yyyyzz, \
                             tk_zzz_yyyzzz, \
                             tk_zzz_yyzzzz, \
                             tk_zzz_yzzzzz, \
                             tk_zzz_zzzzzz, \
                             ts_z_xxxxxx,   \
                             ts_z_xxxxxy,   \
                             ts_z_xxxxxz,   \
                             ts_z_xxxxyy,   \
                             ts_z_xxxxyz,   \
                             ts_z_xxxxzz,   \
                             ts_z_xxxyyy,   \
                             ts_z_xxxyyz,   \
                             ts_z_xxxyzz,   \
                             ts_z_xxxzzz,   \
                             ts_z_xxyyyy,   \
                             ts_z_xxyyyz,   \
                             ts_z_xxyyzz,   \
                             ts_z_xxyzzz,   \
                             ts_z_xxzzzz,   \
                             ts_z_xyyyyy,   \
                             ts_z_xyyyyz,   \
                             ts_z_xyyyzz,   \
                             ts_z_xyyzzz,   \
                             ts_z_xyzzzz,   \
                             ts_z_xzzzzz,   \
                             ts_z_yyyyyy,   \
                             ts_z_yyyyyz,   \
                             ts_z_yyyyzz,   \
                             ts_z_yyyzzz,   \
                             ts_z_yyzzzz,   \
                             ts_z_yzzzzz,   \
                             ts_z_zzzzzz,   \
                             ts_zzz_xxxxxx, \
                             ts_zzz_xxxxxy, \
                             ts_zzz_xxxxxz, \
                             ts_zzz_xxxxyy, \
                             ts_zzz_xxxxyz, \
                             ts_zzz_xxxxzz, \
                             ts_zzz_xxxyyy, \
                             ts_zzz_xxxyyz, \
                             ts_zzz_xxxyzz, \
                             ts_zzz_xxxzzz, \
                             ts_zzz_xxyyyy, \
                             ts_zzz_xxyyyz, \
                             ts_zzz_xxyyzz, \
                             ts_zzz_xxyzzz, \
                             ts_zzz_xxzzzz, \
                             ts_zzz_xyyyyy, \
                             ts_zzz_xyyyyz, \
                             ts_zzz_xyyyzz, \
                             ts_zzz_xyyzzz, \
                             ts_zzz_xyzzzz, \
                             ts_zzz_xzzzzz, \
                             ts_zzz_yyyyyy, \
                             ts_zzz_yyyyyz, \
                             ts_zzz_yyyyzz, \
                             ts_zzz_yyyzzz, \
                             ts_zzz_yyzzzz, \
                             ts_zzz_yzzzzz, \
                             ts_zzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzz_xxxxxx[i] =
            -4.0 * ts_z_xxxxxx[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxxx[i] * fe_0 + tk_zz_xxxxxx[i] * pa_z[i] + 2.0 * ts_zzz_xxxxxx[i] * fz_0;

        tk_zzz_xxxxxy[i] =
            -4.0 * ts_z_xxxxxy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxxy[i] * fe_0 + tk_zz_xxxxxy[i] * pa_z[i] + 2.0 * ts_zzz_xxxxxy[i] * fz_0;

        tk_zzz_xxxxxz[i] = -4.0 * ts_z_xxxxxz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxxz[i] * fe_0 + tk_zz_xxxxx[i] * fe_0 + tk_zz_xxxxxz[i] * pa_z[i] +
                           2.0 * ts_zzz_xxxxxz[i] * fz_0;

        tk_zzz_xxxxyy[i] =
            -4.0 * ts_z_xxxxyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxyy[i] * fe_0 + tk_zz_xxxxyy[i] * pa_z[i] + 2.0 * ts_zzz_xxxxyy[i] * fz_0;

        tk_zzz_xxxxyz[i] = -4.0 * ts_z_xxxxyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxyz[i] * fe_0 + tk_zz_xxxxy[i] * fe_0 + tk_zz_xxxxyz[i] * pa_z[i] +
                           2.0 * ts_zzz_xxxxyz[i] * fz_0;

        tk_zzz_xxxxzz[i] = -4.0 * ts_z_xxxxzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxxzz[i] * fe_0 + 2.0 * tk_zz_xxxxz[i] * fe_0 +
                           tk_zz_xxxxzz[i] * pa_z[i] + 2.0 * ts_zzz_xxxxzz[i] * fz_0;

        tk_zzz_xxxyyy[i] =
            -4.0 * ts_z_xxxyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxyyy[i] * fe_0 + tk_zz_xxxyyy[i] * pa_z[i] + 2.0 * ts_zzz_xxxyyy[i] * fz_0;

        tk_zzz_xxxyyz[i] = -4.0 * ts_z_xxxyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxyyz[i] * fe_0 + tk_zz_xxxyy[i] * fe_0 + tk_zz_xxxyyz[i] * pa_z[i] +
                           2.0 * ts_zzz_xxxyyz[i] * fz_0;

        tk_zzz_xxxyzz[i] = -4.0 * ts_z_xxxyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxyzz[i] * fe_0 + 2.0 * tk_zz_xxxyz[i] * fe_0 +
                           tk_zz_xxxyzz[i] * pa_z[i] + 2.0 * ts_zzz_xxxyzz[i] * fz_0;

        tk_zzz_xxxzzz[i] = -4.0 * ts_z_xxxzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxxzzz[i] * fe_0 + 3.0 * tk_zz_xxxzz[i] * fe_0 +
                           tk_zz_xxxzzz[i] * pa_z[i] + 2.0 * ts_zzz_xxxzzz[i] * fz_0;

        tk_zzz_xxyyyy[i] =
            -4.0 * ts_z_xxyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyyyy[i] * fe_0 + tk_zz_xxyyyy[i] * pa_z[i] + 2.0 * ts_zzz_xxyyyy[i] * fz_0;

        tk_zzz_xxyyyz[i] = -4.0 * ts_z_xxyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyyyz[i] * fe_0 + tk_zz_xxyyy[i] * fe_0 + tk_zz_xxyyyz[i] * pa_z[i] +
                           2.0 * ts_zzz_xxyyyz[i] * fz_0;

        tk_zzz_xxyyzz[i] = -4.0 * ts_z_xxyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyyzz[i] * fe_0 + 2.0 * tk_zz_xxyyz[i] * fe_0 +
                           tk_zz_xxyyzz[i] * pa_z[i] + 2.0 * ts_zzz_xxyyzz[i] * fz_0;

        tk_zzz_xxyzzz[i] = -4.0 * ts_z_xxyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxyzzz[i] * fe_0 + 3.0 * tk_zz_xxyzz[i] * fe_0 +
                           tk_zz_xxyzzz[i] * pa_z[i] + 2.0 * ts_zzz_xxyzzz[i] * fz_0;

        tk_zzz_xxzzzz[i] = -4.0 * ts_z_xxzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xxzzzz[i] * fe_0 + 4.0 * tk_zz_xxzzz[i] * fe_0 +
                           tk_zz_xxzzzz[i] * pa_z[i] + 2.0 * ts_zzz_xxzzzz[i] * fz_0;

        tk_zzz_xyyyyy[i] =
            -4.0 * ts_z_xyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyyyy[i] * fe_0 + tk_zz_xyyyyy[i] * pa_z[i] + 2.0 * ts_zzz_xyyyyy[i] * fz_0;

        tk_zzz_xyyyyz[i] = -4.0 * ts_z_xyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyyyz[i] * fe_0 + tk_zz_xyyyy[i] * fe_0 + tk_zz_xyyyyz[i] * pa_z[i] +
                           2.0 * ts_zzz_xyyyyz[i] * fz_0;

        tk_zzz_xyyyzz[i] = -4.0 * ts_z_xyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyyzz[i] * fe_0 + 2.0 * tk_zz_xyyyz[i] * fe_0 +
                           tk_zz_xyyyzz[i] * pa_z[i] + 2.0 * ts_zzz_xyyyzz[i] * fz_0;

        tk_zzz_xyyzzz[i] = -4.0 * ts_z_xyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyyzzz[i] * fe_0 + 3.0 * tk_zz_xyyzz[i] * fe_0 +
                           tk_zz_xyyzzz[i] * pa_z[i] + 2.0 * ts_zzz_xyyzzz[i] * fz_0;

        tk_zzz_xyzzzz[i] = -4.0 * ts_z_xyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xyzzzz[i] * fe_0 + 4.0 * tk_zz_xyzzz[i] * fe_0 +
                           tk_zz_xyzzzz[i] * pa_z[i] + 2.0 * ts_zzz_xyzzzz[i] * fz_0;

        tk_zzz_xzzzzz[i] = -4.0 * ts_z_xzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_xzzzzz[i] * fe_0 + 5.0 * tk_zz_xzzzz[i] * fe_0 +
                           tk_zz_xzzzzz[i] * pa_z[i] + 2.0 * ts_zzz_xzzzzz[i] * fz_0;

        tk_zzz_yyyyyy[i] =
            -4.0 * ts_z_yyyyyy[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyyyy[i] * fe_0 + tk_zz_yyyyyy[i] * pa_z[i] + 2.0 * ts_zzz_yyyyyy[i] * fz_0;

        tk_zzz_yyyyyz[i] = -4.0 * ts_z_yyyyyz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyyyz[i] * fe_0 + tk_zz_yyyyy[i] * fe_0 + tk_zz_yyyyyz[i] * pa_z[i] +
                           2.0 * ts_zzz_yyyyyz[i] * fz_0;

        tk_zzz_yyyyzz[i] = -4.0 * ts_z_yyyyzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyyzz[i] * fe_0 + 2.0 * tk_zz_yyyyz[i] * fe_0 +
                           tk_zz_yyyyzz[i] * pa_z[i] + 2.0 * ts_zzz_yyyyzz[i] * fz_0;

        tk_zzz_yyyzzz[i] = -4.0 * ts_z_yyyzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyyzzz[i] * fe_0 + 3.0 * tk_zz_yyyzz[i] * fe_0 +
                           tk_zz_yyyzzz[i] * pa_z[i] + 2.0 * ts_zzz_yyyzzz[i] * fz_0;

        tk_zzz_yyzzzz[i] = -4.0 * ts_z_yyzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yyzzzz[i] * fe_0 + 4.0 * tk_zz_yyzzz[i] * fe_0 +
                           tk_zz_yyzzzz[i] * pa_z[i] + 2.0 * ts_zzz_yyzzzz[i] * fz_0;

        tk_zzz_yzzzzz[i] = -4.0 * ts_z_yzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_yzzzzz[i] * fe_0 + 5.0 * tk_zz_yzzzz[i] * fe_0 +
                           tk_zz_yzzzzz[i] * pa_z[i] + 2.0 * ts_zzz_yzzzzz[i] * fz_0;

        tk_zzz_zzzzzz[i] = -4.0 * ts_z_zzzzzz[i] * fbe_0 * fz_0 + 2.0 * tk_z_zzzzzz[i] * fe_0 + 6.0 * tk_zz_zzzzz[i] * fe_0 +
                           tk_zz_zzzzzz[i] * pa_z[i] + 2.0 * ts_zzz_zzzzzz[i] * fz_0;
    }
}

}  // namespace kinrec
