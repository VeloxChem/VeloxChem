#include "TwoCenterElectronRepulsionPrimRecFI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_fi(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fi,
                                const size_t idx_eri_0_pi,
                                const size_t idx_eri_1_pi,
                                const size_t idx_eri_1_dh,
                                const size_t idx_eri_1_di,
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

    // Set up components of auxiliary buffer : PI

    auto g_x_xxxxxx_0 = pbuffer.data(idx_eri_0_pi);

    auto g_x_xxxxxy_0 = pbuffer.data(idx_eri_0_pi + 1);

    auto g_x_xxxxxz_0 = pbuffer.data(idx_eri_0_pi + 2);

    auto g_x_xxxxyy_0 = pbuffer.data(idx_eri_0_pi + 3);

    auto g_x_xxxxyz_0 = pbuffer.data(idx_eri_0_pi + 4);

    auto g_x_xxxxzz_0 = pbuffer.data(idx_eri_0_pi + 5);

    auto g_x_xxxyyy_0 = pbuffer.data(idx_eri_0_pi + 6);

    auto g_x_xxxyyz_0 = pbuffer.data(idx_eri_0_pi + 7);

    auto g_x_xxxyzz_0 = pbuffer.data(idx_eri_0_pi + 8);

    auto g_x_xxxzzz_0 = pbuffer.data(idx_eri_0_pi + 9);

    auto g_x_xxyyyy_0 = pbuffer.data(idx_eri_0_pi + 10);

    auto g_x_xxyyyz_0 = pbuffer.data(idx_eri_0_pi + 11);

    auto g_x_xxyyzz_0 = pbuffer.data(idx_eri_0_pi + 12);

    auto g_x_xxyzzz_0 = pbuffer.data(idx_eri_0_pi + 13);

    auto g_x_xxzzzz_0 = pbuffer.data(idx_eri_0_pi + 14);

    auto g_x_xyyyyy_0 = pbuffer.data(idx_eri_0_pi + 15);

    auto g_x_xyyyyz_0 = pbuffer.data(idx_eri_0_pi + 16);

    auto g_x_xyyyzz_0 = pbuffer.data(idx_eri_0_pi + 17);

    auto g_x_xyyzzz_0 = pbuffer.data(idx_eri_0_pi + 18);

    auto g_x_xyzzzz_0 = pbuffer.data(idx_eri_0_pi + 19);

    auto g_x_xzzzzz_0 = pbuffer.data(idx_eri_0_pi + 20);

    auto g_x_yyyyyy_0 = pbuffer.data(idx_eri_0_pi + 21);

    auto g_x_yyyyyz_0 = pbuffer.data(idx_eri_0_pi + 22);

    auto g_x_yyyyzz_0 = pbuffer.data(idx_eri_0_pi + 23);

    auto g_x_yyyzzz_0 = pbuffer.data(idx_eri_0_pi + 24);

    auto g_x_yyzzzz_0 = pbuffer.data(idx_eri_0_pi + 25);

    auto g_x_yzzzzz_0 = pbuffer.data(idx_eri_0_pi + 26);

    auto g_x_zzzzzz_0 = pbuffer.data(idx_eri_0_pi + 27);

    auto g_y_xxxxxx_0 = pbuffer.data(idx_eri_0_pi + 28);

    auto g_y_xxxxxy_0 = pbuffer.data(idx_eri_0_pi + 29);

    auto g_y_xxxxxz_0 = pbuffer.data(idx_eri_0_pi + 30);

    auto g_y_xxxxyy_0 = pbuffer.data(idx_eri_0_pi + 31);

    auto g_y_xxxxyz_0 = pbuffer.data(idx_eri_0_pi + 32);

    auto g_y_xxxxzz_0 = pbuffer.data(idx_eri_0_pi + 33);

    auto g_y_xxxyyy_0 = pbuffer.data(idx_eri_0_pi + 34);

    auto g_y_xxxyyz_0 = pbuffer.data(idx_eri_0_pi + 35);

    auto g_y_xxxyzz_0 = pbuffer.data(idx_eri_0_pi + 36);

    auto g_y_xxxzzz_0 = pbuffer.data(idx_eri_0_pi + 37);

    auto g_y_xxyyyy_0 = pbuffer.data(idx_eri_0_pi + 38);

    auto g_y_xxyyyz_0 = pbuffer.data(idx_eri_0_pi + 39);

    auto g_y_xxyyzz_0 = pbuffer.data(idx_eri_0_pi + 40);

    auto g_y_xxyzzz_0 = pbuffer.data(idx_eri_0_pi + 41);

    auto g_y_xxzzzz_0 = pbuffer.data(idx_eri_0_pi + 42);

    auto g_y_xyyyyy_0 = pbuffer.data(idx_eri_0_pi + 43);

    auto g_y_xyyyyz_0 = pbuffer.data(idx_eri_0_pi + 44);

    auto g_y_xyyyzz_0 = pbuffer.data(idx_eri_0_pi + 45);

    auto g_y_xyyzzz_0 = pbuffer.data(idx_eri_0_pi + 46);

    auto g_y_xyzzzz_0 = pbuffer.data(idx_eri_0_pi + 47);

    auto g_y_xzzzzz_0 = pbuffer.data(idx_eri_0_pi + 48);

    auto g_y_yyyyyy_0 = pbuffer.data(idx_eri_0_pi + 49);

    auto g_y_yyyyyz_0 = pbuffer.data(idx_eri_0_pi + 50);

    auto g_y_yyyyzz_0 = pbuffer.data(idx_eri_0_pi + 51);

    auto g_y_yyyzzz_0 = pbuffer.data(idx_eri_0_pi + 52);

    auto g_y_yyzzzz_0 = pbuffer.data(idx_eri_0_pi + 53);

    auto g_y_yzzzzz_0 = pbuffer.data(idx_eri_0_pi + 54);

    auto g_y_zzzzzz_0 = pbuffer.data(idx_eri_0_pi + 55);

    auto g_z_xxxxxx_0 = pbuffer.data(idx_eri_0_pi + 56);

    auto g_z_xxxxxy_0 = pbuffer.data(idx_eri_0_pi + 57);

    auto g_z_xxxxxz_0 = pbuffer.data(idx_eri_0_pi + 58);

    auto g_z_xxxxyy_0 = pbuffer.data(idx_eri_0_pi + 59);

    auto g_z_xxxxyz_0 = pbuffer.data(idx_eri_0_pi + 60);

    auto g_z_xxxxzz_0 = pbuffer.data(idx_eri_0_pi + 61);

    auto g_z_xxxyyy_0 = pbuffer.data(idx_eri_0_pi + 62);

    auto g_z_xxxyyz_0 = pbuffer.data(idx_eri_0_pi + 63);

    auto g_z_xxxyzz_0 = pbuffer.data(idx_eri_0_pi + 64);

    auto g_z_xxxzzz_0 = pbuffer.data(idx_eri_0_pi + 65);

    auto g_z_xxyyyy_0 = pbuffer.data(idx_eri_0_pi + 66);

    auto g_z_xxyyyz_0 = pbuffer.data(idx_eri_0_pi + 67);

    auto g_z_xxyyzz_0 = pbuffer.data(idx_eri_0_pi + 68);

    auto g_z_xxyzzz_0 = pbuffer.data(idx_eri_0_pi + 69);

    auto g_z_xxzzzz_0 = pbuffer.data(idx_eri_0_pi + 70);

    auto g_z_xyyyyy_0 = pbuffer.data(idx_eri_0_pi + 71);

    auto g_z_xyyyyz_0 = pbuffer.data(idx_eri_0_pi + 72);

    auto g_z_xyyyzz_0 = pbuffer.data(idx_eri_0_pi + 73);

    auto g_z_xyyzzz_0 = pbuffer.data(idx_eri_0_pi + 74);

    auto g_z_xyzzzz_0 = pbuffer.data(idx_eri_0_pi + 75);

    auto g_z_xzzzzz_0 = pbuffer.data(idx_eri_0_pi + 76);

    auto g_z_yyyyyy_0 = pbuffer.data(idx_eri_0_pi + 77);

    auto g_z_yyyyyz_0 = pbuffer.data(idx_eri_0_pi + 78);

    auto g_z_yyyyzz_0 = pbuffer.data(idx_eri_0_pi + 79);

    auto g_z_yyyzzz_0 = pbuffer.data(idx_eri_0_pi + 80);

    auto g_z_yyzzzz_0 = pbuffer.data(idx_eri_0_pi + 81);

    auto g_z_yzzzzz_0 = pbuffer.data(idx_eri_0_pi + 82);

    auto g_z_zzzzzz_0 = pbuffer.data(idx_eri_0_pi + 83);

    // Set up components of auxiliary buffer : PI

    auto g_x_xxxxxx_1 = pbuffer.data(idx_eri_1_pi);

    auto g_x_xxxxxy_1 = pbuffer.data(idx_eri_1_pi + 1);

    auto g_x_xxxxxz_1 = pbuffer.data(idx_eri_1_pi + 2);

    auto g_x_xxxxyy_1 = pbuffer.data(idx_eri_1_pi + 3);

    auto g_x_xxxxyz_1 = pbuffer.data(idx_eri_1_pi + 4);

    auto g_x_xxxxzz_1 = pbuffer.data(idx_eri_1_pi + 5);

    auto g_x_xxxyyy_1 = pbuffer.data(idx_eri_1_pi + 6);

    auto g_x_xxxyyz_1 = pbuffer.data(idx_eri_1_pi + 7);

    auto g_x_xxxyzz_1 = pbuffer.data(idx_eri_1_pi + 8);

    auto g_x_xxxzzz_1 = pbuffer.data(idx_eri_1_pi + 9);

    auto g_x_xxyyyy_1 = pbuffer.data(idx_eri_1_pi + 10);

    auto g_x_xxyyyz_1 = pbuffer.data(idx_eri_1_pi + 11);

    auto g_x_xxyyzz_1 = pbuffer.data(idx_eri_1_pi + 12);

    auto g_x_xxyzzz_1 = pbuffer.data(idx_eri_1_pi + 13);

    auto g_x_xxzzzz_1 = pbuffer.data(idx_eri_1_pi + 14);

    auto g_x_xyyyyy_1 = pbuffer.data(idx_eri_1_pi + 15);

    auto g_x_xyyyyz_1 = pbuffer.data(idx_eri_1_pi + 16);

    auto g_x_xyyyzz_1 = pbuffer.data(idx_eri_1_pi + 17);

    auto g_x_xyyzzz_1 = pbuffer.data(idx_eri_1_pi + 18);

    auto g_x_xyzzzz_1 = pbuffer.data(idx_eri_1_pi + 19);

    auto g_x_xzzzzz_1 = pbuffer.data(idx_eri_1_pi + 20);

    auto g_x_yyyyyy_1 = pbuffer.data(idx_eri_1_pi + 21);

    auto g_x_yyyyyz_1 = pbuffer.data(idx_eri_1_pi + 22);

    auto g_x_yyyyzz_1 = pbuffer.data(idx_eri_1_pi + 23);

    auto g_x_yyyzzz_1 = pbuffer.data(idx_eri_1_pi + 24);

    auto g_x_yyzzzz_1 = pbuffer.data(idx_eri_1_pi + 25);

    auto g_x_yzzzzz_1 = pbuffer.data(idx_eri_1_pi + 26);

    auto g_x_zzzzzz_1 = pbuffer.data(idx_eri_1_pi + 27);

    auto g_y_xxxxxx_1 = pbuffer.data(idx_eri_1_pi + 28);

    auto g_y_xxxxxy_1 = pbuffer.data(idx_eri_1_pi + 29);

    auto g_y_xxxxxz_1 = pbuffer.data(idx_eri_1_pi + 30);

    auto g_y_xxxxyy_1 = pbuffer.data(idx_eri_1_pi + 31);

    auto g_y_xxxxyz_1 = pbuffer.data(idx_eri_1_pi + 32);

    auto g_y_xxxxzz_1 = pbuffer.data(idx_eri_1_pi + 33);

    auto g_y_xxxyyy_1 = pbuffer.data(idx_eri_1_pi + 34);

    auto g_y_xxxyyz_1 = pbuffer.data(idx_eri_1_pi + 35);

    auto g_y_xxxyzz_1 = pbuffer.data(idx_eri_1_pi + 36);

    auto g_y_xxxzzz_1 = pbuffer.data(idx_eri_1_pi + 37);

    auto g_y_xxyyyy_1 = pbuffer.data(idx_eri_1_pi + 38);

    auto g_y_xxyyyz_1 = pbuffer.data(idx_eri_1_pi + 39);

    auto g_y_xxyyzz_1 = pbuffer.data(idx_eri_1_pi + 40);

    auto g_y_xxyzzz_1 = pbuffer.data(idx_eri_1_pi + 41);

    auto g_y_xxzzzz_1 = pbuffer.data(idx_eri_1_pi + 42);

    auto g_y_xyyyyy_1 = pbuffer.data(idx_eri_1_pi + 43);

    auto g_y_xyyyyz_1 = pbuffer.data(idx_eri_1_pi + 44);

    auto g_y_xyyyzz_1 = pbuffer.data(idx_eri_1_pi + 45);

    auto g_y_xyyzzz_1 = pbuffer.data(idx_eri_1_pi + 46);

    auto g_y_xyzzzz_1 = pbuffer.data(idx_eri_1_pi + 47);

    auto g_y_xzzzzz_1 = pbuffer.data(idx_eri_1_pi + 48);

    auto g_y_yyyyyy_1 = pbuffer.data(idx_eri_1_pi + 49);

    auto g_y_yyyyyz_1 = pbuffer.data(idx_eri_1_pi + 50);

    auto g_y_yyyyzz_1 = pbuffer.data(idx_eri_1_pi + 51);

    auto g_y_yyyzzz_1 = pbuffer.data(idx_eri_1_pi + 52);

    auto g_y_yyzzzz_1 = pbuffer.data(idx_eri_1_pi + 53);

    auto g_y_yzzzzz_1 = pbuffer.data(idx_eri_1_pi + 54);

    auto g_y_zzzzzz_1 = pbuffer.data(idx_eri_1_pi + 55);

    auto g_z_xxxxxx_1 = pbuffer.data(idx_eri_1_pi + 56);

    auto g_z_xxxxxy_1 = pbuffer.data(idx_eri_1_pi + 57);

    auto g_z_xxxxxz_1 = pbuffer.data(idx_eri_1_pi + 58);

    auto g_z_xxxxyy_1 = pbuffer.data(idx_eri_1_pi + 59);

    auto g_z_xxxxyz_1 = pbuffer.data(idx_eri_1_pi + 60);

    auto g_z_xxxxzz_1 = pbuffer.data(idx_eri_1_pi + 61);

    auto g_z_xxxyyy_1 = pbuffer.data(idx_eri_1_pi + 62);

    auto g_z_xxxyyz_1 = pbuffer.data(idx_eri_1_pi + 63);

    auto g_z_xxxyzz_1 = pbuffer.data(idx_eri_1_pi + 64);

    auto g_z_xxxzzz_1 = pbuffer.data(idx_eri_1_pi + 65);

    auto g_z_xxyyyy_1 = pbuffer.data(idx_eri_1_pi + 66);

    auto g_z_xxyyyz_1 = pbuffer.data(idx_eri_1_pi + 67);

    auto g_z_xxyyzz_1 = pbuffer.data(idx_eri_1_pi + 68);

    auto g_z_xxyzzz_1 = pbuffer.data(idx_eri_1_pi + 69);

    auto g_z_xxzzzz_1 = pbuffer.data(idx_eri_1_pi + 70);

    auto g_z_xyyyyy_1 = pbuffer.data(idx_eri_1_pi + 71);

    auto g_z_xyyyyz_1 = pbuffer.data(idx_eri_1_pi + 72);

    auto g_z_xyyyzz_1 = pbuffer.data(idx_eri_1_pi + 73);

    auto g_z_xyyzzz_1 = pbuffer.data(idx_eri_1_pi + 74);

    auto g_z_xyzzzz_1 = pbuffer.data(idx_eri_1_pi + 75);

    auto g_z_xzzzzz_1 = pbuffer.data(idx_eri_1_pi + 76);

    auto g_z_yyyyyy_1 = pbuffer.data(idx_eri_1_pi + 77);

    auto g_z_yyyyyz_1 = pbuffer.data(idx_eri_1_pi + 78);

    auto g_z_yyyyzz_1 = pbuffer.data(idx_eri_1_pi + 79);

    auto g_z_yyyzzz_1 = pbuffer.data(idx_eri_1_pi + 80);

    auto g_z_yyzzzz_1 = pbuffer.data(idx_eri_1_pi + 81);

    auto g_z_yzzzzz_1 = pbuffer.data(idx_eri_1_pi + 82);

    auto g_z_zzzzzz_1 = pbuffer.data(idx_eri_1_pi + 83);

    // Set up components of auxiliary buffer : DH

    auto g_xx_xxxxx_1 = pbuffer.data(idx_eri_1_dh);

    auto g_xx_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 1);

    auto g_xx_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 2);

    auto g_xx_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 3);

    auto g_xx_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 4);

    auto g_xx_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 5);

    auto g_xx_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 6);

    auto g_xx_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 7);

    auto g_xx_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 8);

    auto g_xx_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 9);

    auto g_xx_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 10);

    auto g_xx_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 11);

    auto g_xx_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 12);

    auto g_xx_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 13);

    auto g_xx_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 14);

    auto g_xx_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 15);

    auto g_xx_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 16);

    auto g_xx_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 17);

    auto g_xx_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 18);

    auto g_xx_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 19);

    auto g_xx_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 20);

    auto g_yy_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 63);

    auto g_yy_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 64);

    auto g_yy_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 65);

    auto g_yy_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 66);

    auto g_yy_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 67);

    auto g_yy_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 68);

    auto g_yy_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 69);

    auto g_yy_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 70);

    auto g_yy_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 71);

    auto g_yy_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 72);

    auto g_yy_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 73);

    auto g_yy_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 74);

    auto g_yy_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 75);

    auto g_yy_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 76);

    auto g_yy_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 77);

    auto g_yy_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 78);

    auto g_yy_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 79);

    auto g_yy_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 80);

    auto g_yy_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 81);

    auto g_yy_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 82);

    auto g_yy_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 83);

    auto g_yz_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 88);

    auto g_yz_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 91);

    auto g_yz_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 92);

    auto g_yz_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 95);

    auto g_yz_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 96);

    auto g_yz_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 97);

    auto g_yz_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 100);

    auto g_yz_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 101);

    auto g_yz_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 102);

    auto g_yz_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 103);

    auto g_zz_xxxxx_1 = pbuffer.data(idx_eri_1_dh + 105);

    auto g_zz_xxxxy_1 = pbuffer.data(idx_eri_1_dh + 106);

    auto g_zz_xxxxz_1 = pbuffer.data(idx_eri_1_dh + 107);

    auto g_zz_xxxyy_1 = pbuffer.data(idx_eri_1_dh + 108);

    auto g_zz_xxxyz_1 = pbuffer.data(idx_eri_1_dh + 109);

    auto g_zz_xxxzz_1 = pbuffer.data(idx_eri_1_dh + 110);

    auto g_zz_xxyyy_1 = pbuffer.data(idx_eri_1_dh + 111);

    auto g_zz_xxyyz_1 = pbuffer.data(idx_eri_1_dh + 112);

    auto g_zz_xxyzz_1 = pbuffer.data(idx_eri_1_dh + 113);

    auto g_zz_xxzzz_1 = pbuffer.data(idx_eri_1_dh + 114);

    auto g_zz_xyyyy_1 = pbuffer.data(idx_eri_1_dh + 115);

    auto g_zz_xyyyz_1 = pbuffer.data(idx_eri_1_dh + 116);

    auto g_zz_xyyzz_1 = pbuffer.data(idx_eri_1_dh + 117);

    auto g_zz_xyzzz_1 = pbuffer.data(idx_eri_1_dh + 118);

    auto g_zz_xzzzz_1 = pbuffer.data(idx_eri_1_dh + 119);

    auto g_zz_yyyyy_1 = pbuffer.data(idx_eri_1_dh + 120);

    auto g_zz_yyyyz_1 = pbuffer.data(idx_eri_1_dh + 121);

    auto g_zz_yyyzz_1 = pbuffer.data(idx_eri_1_dh + 122);

    auto g_zz_yyzzz_1 = pbuffer.data(idx_eri_1_dh + 123);

    auto g_zz_yzzzz_1 = pbuffer.data(idx_eri_1_dh + 124);

    auto g_zz_zzzzz_1 = pbuffer.data(idx_eri_1_dh + 125);

    // Set up components of auxiliary buffer : DI

    auto g_xx_xxxxxx_1 = pbuffer.data(idx_eri_1_di);

    auto g_xx_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 1);

    auto g_xx_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 2);

    auto g_xx_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 3);

    auto g_xx_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 4);

    auto g_xx_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 5);

    auto g_xx_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 6);

    auto g_xx_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 7);

    auto g_xx_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 8);

    auto g_xx_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 9);

    auto g_xx_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 10);

    auto g_xx_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 11);

    auto g_xx_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 12);

    auto g_xx_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 13);

    auto g_xx_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 14);

    auto g_xx_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 15);

    auto g_xx_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 16);

    auto g_xx_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 17);

    auto g_xx_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 18);

    auto g_xx_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 19);

    auto g_xx_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 20);

    auto g_xx_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 21);

    auto g_xx_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 22);

    auto g_xx_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 23);

    auto g_xx_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 24);

    auto g_xx_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 25);

    auto g_xx_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 26);

    auto g_xx_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 27);

    auto g_xy_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 29);

    auto g_xy_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 31);

    auto g_xy_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 34);

    auto g_xy_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 38);

    auto g_xy_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 43);

    auto g_xz_xxxxxx_1 = pbuffer.data(idx_eri_1_di + 56);

    auto g_xz_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 58);

    auto g_xz_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 61);

    auto g_xz_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 65);

    auto g_xz_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 70);

    auto g_xz_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 76);

    auto g_yy_xxxxxx_1 = pbuffer.data(idx_eri_1_di + 84);

    auto g_yy_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 85);

    auto g_yy_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 86);

    auto g_yy_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 87);

    auto g_yy_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 88);

    auto g_yy_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 89);

    auto g_yy_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 90);

    auto g_yy_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 91);

    auto g_yy_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 92);

    auto g_yy_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 93);

    auto g_yy_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 94);

    auto g_yy_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 95);

    auto g_yy_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 96);

    auto g_yy_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 97);

    auto g_yy_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 98);

    auto g_yy_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 99);

    auto g_yy_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 100);

    auto g_yy_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 101);

    auto g_yy_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 102);

    auto g_yy_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 103);

    auto g_yy_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 104);

    auto g_yy_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 105);

    auto g_yy_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 106);

    auto g_yy_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 107);

    auto g_yy_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 108);

    auto g_yy_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 109);

    auto g_yy_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 110);

    auto g_yy_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 111);

    auto g_yz_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 116);

    auto g_yz_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 119);

    auto g_yz_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 120);

    auto g_yz_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 123);

    auto g_yz_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 124);

    auto g_yz_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 125);

    auto g_yz_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 128);

    auto g_yz_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 129);

    auto g_yz_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 130);

    auto g_yz_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 131);

    auto g_yz_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 133);

    auto g_yz_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 134);

    auto g_yz_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 135);

    auto g_yz_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 136);

    auto g_yz_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 137);

    auto g_yz_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 138);

    auto g_yz_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 139);

    auto g_zz_xxxxxx_1 = pbuffer.data(idx_eri_1_di + 140);

    auto g_zz_xxxxxy_1 = pbuffer.data(idx_eri_1_di + 141);

    auto g_zz_xxxxxz_1 = pbuffer.data(idx_eri_1_di + 142);

    auto g_zz_xxxxyy_1 = pbuffer.data(idx_eri_1_di + 143);

    auto g_zz_xxxxyz_1 = pbuffer.data(idx_eri_1_di + 144);

    auto g_zz_xxxxzz_1 = pbuffer.data(idx_eri_1_di + 145);

    auto g_zz_xxxyyy_1 = pbuffer.data(idx_eri_1_di + 146);

    auto g_zz_xxxyyz_1 = pbuffer.data(idx_eri_1_di + 147);

    auto g_zz_xxxyzz_1 = pbuffer.data(idx_eri_1_di + 148);

    auto g_zz_xxxzzz_1 = pbuffer.data(idx_eri_1_di + 149);

    auto g_zz_xxyyyy_1 = pbuffer.data(idx_eri_1_di + 150);

    auto g_zz_xxyyyz_1 = pbuffer.data(idx_eri_1_di + 151);

    auto g_zz_xxyyzz_1 = pbuffer.data(idx_eri_1_di + 152);

    auto g_zz_xxyzzz_1 = pbuffer.data(idx_eri_1_di + 153);

    auto g_zz_xxzzzz_1 = pbuffer.data(idx_eri_1_di + 154);

    auto g_zz_xyyyyy_1 = pbuffer.data(idx_eri_1_di + 155);

    auto g_zz_xyyyyz_1 = pbuffer.data(idx_eri_1_di + 156);

    auto g_zz_xyyyzz_1 = pbuffer.data(idx_eri_1_di + 157);

    auto g_zz_xyyzzz_1 = pbuffer.data(idx_eri_1_di + 158);

    auto g_zz_xyzzzz_1 = pbuffer.data(idx_eri_1_di + 159);

    auto g_zz_xzzzzz_1 = pbuffer.data(idx_eri_1_di + 160);

    auto g_zz_yyyyyy_1 = pbuffer.data(idx_eri_1_di + 161);

    auto g_zz_yyyyyz_1 = pbuffer.data(idx_eri_1_di + 162);

    auto g_zz_yyyyzz_1 = pbuffer.data(idx_eri_1_di + 163);

    auto g_zz_yyyzzz_1 = pbuffer.data(idx_eri_1_di + 164);

    auto g_zz_yyzzzz_1 = pbuffer.data(idx_eri_1_di + 165);

    auto g_zz_yzzzzz_1 = pbuffer.data(idx_eri_1_di + 166);

    auto g_zz_zzzzzz_1 = pbuffer.data(idx_eri_1_di + 167);

    // Set up 0-28 components of targeted buffer : FI

    auto g_xxx_xxxxxx_0 = pbuffer.data(idx_eri_0_fi);

    auto g_xxx_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 1);

    auto g_xxx_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 2);

    auto g_xxx_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 3);

    auto g_xxx_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 4);

    auto g_xxx_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 5);

    auto g_xxx_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 6);

    auto g_xxx_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 7);

    auto g_xxx_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 8);

    auto g_xxx_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 9);

    auto g_xxx_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 10);

    auto g_xxx_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 11);

    auto g_xxx_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 12);

    auto g_xxx_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 13);

    auto g_xxx_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 14);

    auto g_xxx_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 15);

    auto g_xxx_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 16);

    auto g_xxx_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 17);

    auto g_xxx_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 18);

    auto g_xxx_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 19);

    auto g_xxx_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 20);

    auto g_xxx_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 21);

    auto g_xxx_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 22);

    auto g_xxx_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 23);

    auto g_xxx_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 24);

    auto g_xxx_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 25);

    auto g_xxx_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 26);

    auto g_xxx_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 27);

    #pragma omp simd aligned(g_x_xxxxxx_0, g_x_xxxxxx_1, g_x_xxxxxy_0, g_x_xxxxxy_1, g_x_xxxxxz_0, g_x_xxxxxz_1, g_x_xxxxyy_0, g_x_xxxxyy_1, g_x_xxxxyz_0, g_x_xxxxyz_1, g_x_xxxxzz_0, g_x_xxxxzz_1, g_x_xxxyyy_0, g_x_xxxyyy_1, g_x_xxxyyz_0, g_x_xxxyyz_1, g_x_xxxyzz_0, g_x_xxxyzz_1, g_x_xxxzzz_0, g_x_xxxzzz_1, g_x_xxyyyy_0, g_x_xxyyyy_1, g_x_xxyyyz_0, g_x_xxyyyz_1, g_x_xxyyzz_0, g_x_xxyyzz_1, g_x_xxyzzz_0, g_x_xxyzzz_1, g_x_xxzzzz_0, g_x_xxzzzz_1, g_x_xyyyyy_0, g_x_xyyyyy_1, g_x_xyyyyz_0, g_x_xyyyyz_1, g_x_xyyyzz_0, g_x_xyyyzz_1, g_x_xyyzzz_0, g_x_xyyzzz_1, g_x_xyzzzz_0, g_x_xyzzzz_1, g_x_xzzzzz_0, g_x_xzzzzz_1, g_x_yyyyyy_0, g_x_yyyyyy_1, g_x_yyyyyz_0, g_x_yyyyyz_1, g_x_yyyyzz_0, g_x_yyyyzz_1, g_x_yyyzzz_0, g_x_yyyzzz_1, g_x_yyzzzz_0, g_x_yyzzzz_1, g_x_yzzzzz_0, g_x_yzzzzz_1, g_x_zzzzzz_0, g_x_zzzzzz_1, g_xx_xxxxx_1, g_xx_xxxxxx_1, g_xx_xxxxxy_1, g_xx_xxxxxz_1, g_xx_xxxxy_1, g_xx_xxxxyy_1, g_xx_xxxxyz_1, g_xx_xxxxz_1, g_xx_xxxxzz_1, g_xx_xxxyy_1, g_xx_xxxyyy_1, g_xx_xxxyyz_1, g_xx_xxxyz_1, g_xx_xxxyzz_1, g_xx_xxxzz_1, g_xx_xxxzzz_1, g_xx_xxyyy_1, g_xx_xxyyyy_1, g_xx_xxyyyz_1, g_xx_xxyyz_1, g_xx_xxyyzz_1, g_xx_xxyzz_1, g_xx_xxyzzz_1, g_xx_xxzzz_1, g_xx_xxzzzz_1, g_xx_xyyyy_1, g_xx_xyyyyy_1, g_xx_xyyyyz_1, g_xx_xyyyz_1, g_xx_xyyyzz_1, g_xx_xyyzz_1, g_xx_xyyzzz_1, g_xx_xyzzz_1, g_xx_xyzzzz_1, g_xx_xzzzz_1, g_xx_xzzzzz_1, g_xx_yyyyy_1, g_xx_yyyyyy_1, g_xx_yyyyyz_1, g_xx_yyyyz_1, g_xx_yyyyzz_1, g_xx_yyyzz_1, g_xx_yyyzzz_1, g_xx_yyzzz_1, g_xx_yyzzzz_1, g_xx_yzzzz_1, g_xx_yzzzzz_1, g_xx_zzzzz_1, g_xx_zzzzzz_1, g_xxx_xxxxxx_0, g_xxx_xxxxxy_0, g_xxx_xxxxxz_0, g_xxx_xxxxyy_0, g_xxx_xxxxyz_0, g_xxx_xxxxzz_0, g_xxx_xxxyyy_0, g_xxx_xxxyyz_0, g_xxx_xxxyzz_0, g_xxx_xxxzzz_0, g_xxx_xxyyyy_0, g_xxx_xxyyyz_0, g_xxx_xxyyzz_0, g_xxx_xxyzzz_0, g_xxx_xxzzzz_0, g_xxx_xyyyyy_0, g_xxx_xyyyyz_0, g_xxx_xyyyzz_0, g_xxx_xyyzzz_0, g_xxx_xyzzzz_0, g_xxx_xzzzzz_0, g_xxx_yyyyyy_0, g_xxx_yyyyyz_0, g_xxx_yyyyzz_0, g_xxx_yyyzzz_0, g_xxx_yyzzzz_0, g_xxx_yzzzzz_0, g_xxx_zzzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxx_xxxxxx_0[i] = 2.0 * g_x_xxxxxx_0[i] * fbe_0 - 2.0 * g_x_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xx_xxxxx_1[i] * fe_0 + g_xx_xxxxxx_1[i] * pa_x[i];

        g_xxx_xxxxxy_0[i] = 2.0 * g_x_xxxxxy_0[i] * fbe_0 - 2.0 * g_x_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xx_xxxxy_1[i] * fe_0 + g_xx_xxxxxy_1[i] * pa_x[i];

        g_xxx_xxxxxz_0[i] = 2.0 * g_x_xxxxxz_0[i] * fbe_0 - 2.0 * g_x_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xx_xxxxz_1[i] * fe_0 + g_xx_xxxxxz_1[i] * pa_x[i];

        g_xxx_xxxxyy_0[i] = 2.0 * g_x_xxxxyy_0[i] * fbe_0 - 2.0 * g_x_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xx_xxxyy_1[i] * fe_0 + g_xx_xxxxyy_1[i] * pa_x[i];

        g_xxx_xxxxyz_0[i] = 2.0 * g_x_xxxxyz_0[i] * fbe_0 - 2.0 * g_x_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xx_xxxyz_1[i] * fe_0 + g_xx_xxxxyz_1[i] * pa_x[i];

        g_xxx_xxxxzz_0[i] = 2.0 * g_x_xxxxzz_0[i] * fbe_0 - 2.0 * g_x_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xx_xxxzz_1[i] * fe_0 + g_xx_xxxxzz_1[i] * pa_x[i];

        g_xxx_xxxyyy_0[i] = 2.0 * g_x_xxxyyy_0[i] * fbe_0 - 2.0 * g_x_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xx_xxyyy_1[i] * fe_0 + g_xx_xxxyyy_1[i] * pa_x[i];

        g_xxx_xxxyyz_0[i] = 2.0 * g_x_xxxyyz_0[i] * fbe_0 - 2.0 * g_x_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xx_xxyyz_1[i] * fe_0 + g_xx_xxxyyz_1[i] * pa_x[i];

        g_xxx_xxxyzz_0[i] = 2.0 * g_x_xxxyzz_0[i] * fbe_0 - 2.0 * g_x_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xx_xxyzz_1[i] * fe_0 + g_xx_xxxyzz_1[i] * pa_x[i];

        g_xxx_xxxzzz_0[i] = 2.0 * g_x_xxxzzz_0[i] * fbe_0 - 2.0 * g_x_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xx_xxzzz_1[i] * fe_0 + g_xx_xxxzzz_1[i] * pa_x[i];

        g_xxx_xxyyyy_0[i] = 2.0 * g_x_xxyyyy_0[i] * fbe_0 - 2.0 * g_x_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xx_xyyyy_1[i] * fe_0 + g_xx_xxyyyy_1[i] * pa_x[i];

        g_xxx_xxyyyz_0[i] = 2.0 * g_x_xxyyyz_0[i] * fbe_0 - 2.0 * g_x_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xx_xyyyz_1[i] * fe_0 + g_xx_xxyyyz_1[i] * pa_x[i];

        g_xxx_xxyyzz_0[i] = 2.0 * g_x_xxyyzz_0[i] * fbe_0 - 2.0 * g_x_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xx_xyyzz_1[i] * fe_0 + g_xx_xxyyzz_1[i] * pa_x[i];

        g_xxx_xxyzzz_0[i] = 2.0 * g_x_xxyzzz_0[i] * fbe_0 - 2.0 * g_x_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xx_xyzzz_1[i] * fe_0 + g_xx_xxyzzz_1[i] * pa_x[i];

        g_xxx_xxzzzz_0[i] = 2.0 * g_x_xxzzzz_0[i] * fbe_0 - 2.0 * g_x_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xx_xzzzz_1[i] * fe_0 + g_xx_xxzzzz_1[i] * pa_x[i];

        g_xxx_xyyyyy_0[i] = 2.0 * g_x_xyyyyy_0[i] * fbe_0 - 2.0 * g_x_xyyyyy_1[i] * fz_be_0 + g_xx_yyyyy_1[i] * fe_0 + g_xx_xyyyyy_1[i] * pa_x[i];

        g_xxx_xyyyyz_0[i] = 2.0 * g_x_xyyyyz_0[i] * fbe_0 - 2.0 * g_x_xyyyyz_1[i] * fz_be_0 + g_xx_yyyyz_1[i] * fe_0 + g_xx_xyyyyz_1[i] * pa_x[i];

        g_xxx_xyyyzz_0[i] = 2.0 * g_x_xyyyzz_0[i] * fbe_0 - 2.0 * g_x_xyyyzz_1[i] * fz_be_0 + g_xx_yyyzz_1[i] * fe_0 + g_xx_xyyyzz_1[i] * pa_x[i];

        g_xxx_xyyzzz_0[i] = 2.0 * g_x_xyyzzz_0[i] * fbe_0 - 2.0 * g_x_xyyzzz_1[i] * fz_be_0 + g_xx_yyzzz_1[i] * fe_0 + g_xx_xyyzzz_1[i] * pa_x[i];

        g_xxx_xyzzzz_0[i] = 2.0 * g_x_xyzzzz_0[i] * fbe_0 - 2.0 * g_x_xyzzzz_1[i] * fz_be_0 + g_xx_yzzzz_1[i] * fe_0 + g_xx_xyzzzz_1[i] * pa_x[i];

        g_xxx_xzzzzz_0[i] = 2.0 * g_x_xzzzzz_0[i] * fbe_0 - 2.0 * g_x_xzzzzz_1[i] * fz_be_0 + g_xx_zzzzz_1[i] * fe_0 + g_xx_xzzzzz_1[i] * pa_x[i];

        g_xxx_yyyyyy_0[i] = 2.0 * g_x_yyyyyy_0[i] * fbe_0 - 2.0 * g_x_yyyyyy_1[i] * fz_be_0 + g_xx_yyyyyy_1[i] * pa_x[i];

        g_xxx_yyyyyz_0[i] = 2.0 * g_x_yyyyyz_0[i] * fbe_0 - 2.0 * g_x_yyyyyz_1[i] * fz_be_0 + g_xx_yyyyyz_1[i] * pa_x[i];

        g_xxx_yyyyzz_0[i] = 2.0 * g_x_yyyyzz_0[i] * fbe_0 - 2.0 * g_x_yyyyzz_1[i] * fz_be_0 + g_xx_yyyyzz_1[i] * pa_x[i];

        g_xxx_yyyzzz_0[i] = 2.0 * g_x_yyyzzz_0[i] * fbe_0 - 2.0 * g_x_yyyzzz_1[i] * fz_be_0 + g_xx_yyyzzz_1[i] * pa_x[i];

        g_xxx_yyzzzz_0[i] = 2.0 * g_x_yyzzzz_0[i] * fbe_0 - 2.0 * g_x_yyzzzz_1[i] * fz_be_0 + g_xx_yyzzzz_1[i] * pa_x[i];

        g_xxx_yzzzzz_0[i] = 2.0 * g_x_yzzzzz_0[i] * fbe_0 - 2.0 * g_x_yzzzzz_1[i] * fz_be_0 + g_xx_yzzzzz_1[i] * pa_x[i];

        g_xxx_zzzzzz_0[i] = 2.0 * g_x_zzzzzz_0[i] * fbe_0 - 2.0 * g_x_zzzzzz_1[i] * fz_be_0 + g_xx_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : FI

    auto g_xxy_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 28);

    auto g_xxy_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 29);

    auto g_xxy_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 30);

    auto g_xxy_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 31);

    auto g_xxy_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 32);

    auto g_xxy_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 33);

    auto g_xxy_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 34);

    auto g_xxy_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 35);

    auto g_xxy_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 36);

    auto g_xxy_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 37);

    auto g_xxy_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 38);

    auto g_xxy_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 39);

    auto g_xxy_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 40);

    auto g_xxy_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 41);

    auto g_xxy_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 42);

    auto g_xxy_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 43);

    auto g_xxy_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 44);

    auto g_xxy_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 45);

    auto g_xxy_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 46);

    auto g_xxy_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 47);

    auto g_xxy_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 48);

    auto g_xxy_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 49);

    auto g_xxy_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 50);

    auto g_xxy_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 51);

    auto g_xxy_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 52);

    auto g_xxy_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 53);

    auto g_xxy_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 54);

    auto g_xxy_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 55);

    #pragma omp simd aligned(g_xx_xxxxx_1, g_xx_xxxxxx_1, g_xx_xxxxxy_1, g_xx_xxxxxz_1, g_xx_xxxxy_1, g_xx_xxxxyy_1, g_xx_xxxxyz_1, g_xx_xxxxz_1, g_xx_xxxxzz_1, g_xx_xxxyy_1, g_xx_xxxyyy_1, g_xx_xxxyyz_1, g_xx_xxxyz_1, g_xx_xxxyzz_1, g_xx_xxxzz_1, g_xx_xxxzzz_1, g_xx_xxyyy_1, g_xx_xxyyyy_1, g_xx_xxyyyz_1, g_xx_xxyyz_1, g_xx_xxyyzz_1, g_xx_xxyzz_1, g_xx_xxyzzz_1, g_xx_xxzzz_1, g_xx_xxzzzz_1, g_xx_xyyyy_1, g_xx_xyyyyy_1, g_xx_xyyyyz_1, g_xx_xyyyz_1, g_xx_xyyyzz_1, g_xx_xyyzz_1, g_xx_xyyzzz_1, g_xx_xyzzz_1, g_xx_xyzzzz_1, g_xx_xzzzz_1, g_xx_xzzzzz_1, g_xx_yyyyy_1, g_xx_yyyyyy_1, g_xx_yyyyyz_1, g_xx_yyyyz_1, g_xx_yyyyzz_1, g_xx_yyyzz_1, g_xx_yyyzzz_1, g_xx_yyzzz_1, g_xx_yyzzzz_1, g_xx_yzzzz_1, g_xx_yzzzzz_1, g_xx_zzzzz_1, g_xx_zzzzzz_1, g_xxy_xxxxxx_0, g_xxy_xxxxxy_0, g_xxy_xxxxxz_0, g_xxy_xxxxyy_0, g_xxy_xxxxyz_0, g_xxy_xxxxzz_0, g_xxy_xxxyyy_0, g_xxy_xxxyyz_0, g_xxy_xxxyzz_0, g_xxy_xxxzzz_0, g_xxy_xxyyyy_0, g_xxy_xxyyyz_0, g_xxy_xxyyzz_0, g_xxy_xxyzzz_0, g_xxy_xxzzzz_0, g_xxy_xyyyyy_0, g_xxy_xyyyyz_0, g_xxy_xyyyzz_0, g_xxy_xyyzzz_0, g_xxy_xyzzzz_0, g_xxy_xzzzzz_0, g_xxy_yyyyyy_0, g_xxy_yyyyyz_0, g_xxy_yyyyzz_0, g_xxy_yyyzzz_0, g_xxy_yyzzzz_0, g_xxy_yzzzzz_0, g_xxy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxy_xxxxxx_0[i] = g_xx_xxxxxx_1[i] * pa_y[i];

        g_xxy_xxxxxy_0[i] = g_xx_xxxxx_1[i] * fe_0 + g_xx_xxxxxy_1[i] * pa_y[i];

        g_xxy_xxxxxz_0[i] = g_xx_xxxxxz_1[i] * pa_y[i];

        g_xxy_xxxxyy_0[i] = 2.0 * g_xx_xxxxy_1[i] * fe_0 + g_xx_xxxxyy_1[i] * pa_y[i];

        g_xxy_xxxxyz_0[i] = g_xx_xxxxz_1[i] * fe_0 + g_xx_xxxxyz_1[i] * pa_y[i];

        g_xxy_xxxxzz_0[i] = g_xx_xxxxzz_1[i] * pa_y[i];

        g_xxy_xxxyyy_0[i] = 3.0 * g_xx_xxxyy_1[i] * fe_0 + g_xx_xxxyyy_1[i] * pa_y[i];

        g_xxy_xxxyyz_0[i] = 2.0 * g_xx_xxxyz_1[i] * fe_0 + g_xx_xxxyyz_1[i] * pa_y[i];

        g_xxy_xxxyzz_0[i] = g_xx_xxxzz_1[i] * fe_0 + g_xx_xxxyzz_1[i] * pa_y[i];

        g_xxy_xxxzzz_0[i] = g_xx_xxxzzz_1[i] * pa_y[i];

        g_xxy_xxyyyy_0[i] = 4.0 * g_xx_xxyyy_1[i] * fe_0 + g_xx_xxyyyy_1[i] * pa_y[i];

        g_xxy_xxyyyz_0[i] = 3.0 * g_xx_xxyyz_1[i] * fe_0 + g_xx_xxyyyz_1[i] * pa_y[i];

        g_xxy_xxyyzz_0[i] = 2.0 * g_xx_xxyzz_1[i] * fe_0 + g_xx_xxyyzz_1[i] * pa_y[i];

        g_xxy_xxyzzz_0[i] = g_xx_xxzzz_1[i] * fe_0 + g_xx_xxyzzz_1[i] * pa_y[i];

        g_xxy_xxzzzz_0[i] = g_xx_xxzzzz_1[i] * pa_y[i];

        g_xxy_xyyyyy_0[i] = 5.0 * g_xx_xyyyy_1[i] * fe_0 + g_xx_xyyyyy_1[i] * pa_y[i];

        g_xxy_xyyyyz_0[i] = 4.0 * g_xx_xyyyz_1[i] * fe_0 + g_xx_xyyyyz_1[i] * pa_y[i];

        g_xxy_xyyyzz_0[i] = 3.0 * g_xx_xyyzz_1[i] * fe_0 + g_xx_xyyyzz_1[i] * pa_y[i];

        g_xxy_xyyzzz_0[i] = 2.0 * g_xx_xyzzz_1[i] * fe_0 + g_xx_xyyzzz_1[i] * pa_y[i];

        g_xxy_xyzzzz_0[i] = g_xx_xzzzz_1[i] * fe_0 + g_xx_xyzzzz_1[i] * pa_y[i];

        g_xxy_xzzzzz_0[i] = g_xx_xzzzzz_1[i] * pa_y[i];

        g_xxy_yyyyyy_0[i] = 6.0 * g_xx_yyyyy_1[i] * fe_0 + g_xx_yyyyyy_1[i] * pa_y[i];

        g_xxy_yyyyyz_0[i] = 5.0 * g_xx_yyyyz_1[i] * fe_0 + g_xx_yyyyyz_1[i] * pa_y[i];

        g_xxy_yyyyzz_0[i] = 4.0 * g_xx_yyyzz_1[i] * fe_0 + g_xx_yyyyzz_1[i] * pa_y[i];

        g_xxy_yyyzzz_0[i] = 3.0 * g_xx_yyzzz_1[i] * fe_0 + g_xx_yyyzzz_1[i] * pa_y[i];

        g_xxy_yyzzzz_0[i] = 2.0 * g_xx_yzzzz_1[i] * fe_0 + g_xx_yyzzzz_1[i] * pa_y[i];

        g_xxy_yzzzzz_0[i] = g_xx_zzzzz_1[i] * fe_0 + g_xx_yzzzzz_1[i] * pa_y[i];

        g_xxy_zzzzzz_0[i] = g_xx_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : FI

    auto g_xxz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 56);

    auto g_xxz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 57);

    auto g_xxz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 58);

    auto g_xxz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 59);

    auto g_xxz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 60);

    auto g_xxz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 61);

    auto g_xxz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 62);

    auto g_xxz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 63);

    auto g_xxz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 64);

    auto g_xxz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 65);

    auto g_xxz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 66);

    auto g_xxz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 67);

    auto g_xxz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 68);

    auto g_xxz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 69);

    auto g_xxz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 70);

    auto g_xxz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 71);

    auto g_xxz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 72);

    auto g_xxz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 73);

    auto g_xxz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 74);

    auto g_xxz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 75);

    auto g_xxz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 76);

    auto g_xxz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 77);

    auto g_xxz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 78);

    auto g_xxz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 79);

    auto g_xxz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 80);

    auto g_xxz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 81);

    auto g_xxz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 82);

    auto g_xxz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 83);

    #pragma omp simd aligned(g_xx_xxxxx_1, g_xx_xxxxxx_1, g_xx_xxxxxy_1, g_xx_xxxxxz_1, g_xx_xxxxy_1, g_xx_xxxxyy_1, g_xx_xxxxyz_1, g_xx_xxxxz_1, g_xx_xxxxzz_1, g_xx_xxxyy_1, g_xx_xxxyyy_1, g_xx_xxxyyz_1, g_xx_xxxyz_1, g_xx_xxxyzz_1, g_xx_xxxzz_1, g_xx_xxxzzz_1, g_xx_xxyyy_1, g_xx_xxyyyy_1, g_xx_xxyyyz_1, g_xx_xxyyz_1, g_xx_xxyyzz_1, g_xx_xxyzz_1, g_xx_xxyzzz_1, g_xx_xxzzz_1, g_xx_xxzzzz_1, g_xx_xyyyy_1, g_xx_xyyyyy_1, g_xx_xyyyyz_1, g_xx_xyyyz_1, g_xx_xyyyzz_1, g_xx_xyyzz_1, g_xx_xyyzzz_1, g_xx_xyzzz_1, g_xx_xyzzzz_1, g_xx_xzzzz_1, g_xx_xzzzzz_1, g_xx_yyyyy_1, g_xx_yyyyyy_1, g_xx_yyyyyz_1, g_xx_yyyyz_1, g_xx_yyyyzz_1, g_xx_yyyzz_1, g_xx_yyyzzz_1, g_xx_yyzzz_1, g_xx_yyzzzz_1, g_xx_yzzzz_1, g_xx_yzzzzz_1, g_xx_zzzzz_1, g_xx_zzzzzz_1, g_xxz_xxxxxx_0, g_xxz_xxxxxy_0, g_xxz_xxxxxz_0, g_xxz_xxxxyy_0, g_xxz_xxxxyz_0, g_xxz_xxxxzz_0, g_xxz_xxxyyy_0, g_xxz_xxxyyz_0, g_xxz_xxxyzz_0, g_xxz_xxxzzz_0, g_xxz_xxyyyy_0, g_xxz_xxyyyz_0, g_xxz_xxyyzz_0, g_xxz_xxyzzz_0, g_xxz_xxzzzz_0, g_xxz_xyyyyy_0, g_xxz_xyyyyz_0, g_xxz_xyyyzz_0, g_xxz_xyyzzz_0, g_xxz_xyzzzz_0, g_xxz_xzzzzz_0, g_xxz_yyyyyy_0, g_xxz_yyyyyz_0, g_xxz_yyyyzz_0, g_xxz_yyyzzz_0, g_xxz_yyzzzz_0, g_xxz_yzzzzz_0, g_xxz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxz_xxxxxx_0[i] = g_xx_xxxxxx_1[i] * pa_z[i];

        g_xxz_xxxxxy_0[i] = g_xx_xxxxxy_1[i] * pa_z[i];

        g_xxz_xxxxxz_0[i] = g_xx_xxxxx_1[i] * fe_0 + g_xx_xxxxxz_1[i] * pa_z[i];

        g_xxz_xxxxyy_0[i] = g_xx_xxxxyy_1[i] * pa_z[i];

        g_xxz_xxxxyz_0[i] = g_xx_xxxxy_1[i] * fe_0 + g_xx_xxxxyz_1[i] * pa_z[i];

        g_xxz_xxxxzz_0[i] = 2.0 * g_xx_xxxxz_1[i] * fe_0 + g_xx_xxxxzz_1[i] * pa_z[i];

        g_xxz_xxxyyy_0[i] = g_xx_xxxyyy_1[i] * pa_z[i];

        g_xxz_xxxyyz_0[i] = g_xx_xxxyy_1[i] * fe_0 + g_xx_xxxyyz_1[i] * pa_z[i];

        g_xxz_xxxyzz_0[i] = 2.0 * g_xx_xxxyz_1[i] * fe_0 + g_xx_xxxyzz_1[i] * pa_z[i];

        g_xxz_xxxzzz_0[i] = 3.0 * g_xx_xxxzz_1[i] * fe_0 + g_xx_xxxzzz_1[i] * pa_z[i];

        g_xxz_xxyyyy_0[i] = g_xx_xxyyyy_1[i] * pa_z[i];

        g_xxz_xxyyyz_0[i] = g_xx_xxyyy_1[i] * fe_0 + g_xx_xxyyyz_1[i] * pa_z[i];

        g_xxz_xxyyzz_0[i] = 2.0 * g_xx_xxyyz_1[i] * fe_0 + g_xx_xxyyzz_1[i] * pa_z[i];

        g_xxz_xxyzzz_0[i] = 3.0 * g_xx_xxyzz_1[i] * fe_0 + g_xx_xxyzzz_1[i] * pa_z[i];

        g_xxz_xxzzzz_0[i] = 4.0 * g_xx_xxzzz_1[i] * fe_0 + g_xx_xxzzzz_1[i] * pa_z[i];

        g_xxz_xyyyyy_0[i] = g_xx_xyyyyy_1[i] * pa_z[i];

        g_xxz_xyyyyz_0[i] = g_xx_xyyyy_1[i] * fe_0 + g_xx_xyyyyz_1[i] * pa_z[i];

        g_xxz_xyyyzz_0[i] = 2.0 * g_xx_xyyyz_1[i] * fe_0 + g_xx_xyyyzz_1[i] * pa_z[i];

        g_xxz_xyyzzz_0[i] = 3.0 * g_xx_xyyzz_1[i] * fe_0 + g_xx_xyyzzz_1[i] * pa_z[i];

        g_xxz_xyzzzz_0[i] = 4.0 * g_xx_xyzzz_1[i] * fe_0 + g_xx_xyzzzz_1[i] * pa_z[i];

        g_xxz_xzzzzz_0[i] = 5.0 * g_xx_xzzzz_1[i] * fe_0 + g_xx_xzzzzz_1[i] * pa_z[i];

        g_xxz_yyyyyy_0[i] = g_xx_yyyyyy_1[i] * pa_z[i];

        g_xxz_yyyyyz_0[i] = g_xx_yyyyy_1[i] * fe_0 + g_xx_yyyyyz_1[i] * pa_z[i];

        g_xxz_yyyyzz_0[i] = 2.0 * g_xx_yyyyz_1[i] * fe_0 + g_xx_yyyyzz_1[i] * pa_z[i];

        g_xxz_yyyzzz_0[i] = 3.0 * g_xx_yyyzz_1[i] * fe_0 + g_xx_yyyzzz_1[i] * pa_z[i];

        g_xxz_yyzzzz_0[i] = 4.0 * g_xx_yyzzz_1[i] * fe_0 + g_xx_yyzzzz_1[i] * pa_z[i];

        g_xxz_yzzzzz_0[i] = 5.0 * g_xx_yzzzz_1[i] * fe_0 + g_xx_yzzzzz_1[i] * pa_z[i];

        g_xxz_zzzzzz_0[i] = 6.0 * g_xx_zzzzz_1[i] * fe_0 + g_xx_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : FI

    auto g_xyy_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 84);

    auto g_xyy_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 85);

    auto g_xyy_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 86);

    auto g_xyy_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 87);

    auto g_xyy_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 88);

    auto g_xyy_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 89);

    auto g_xyy_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 90);

    auto g_xyy_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 91);

    auto g_xyy_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 92);

    auto g_xyy_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 93);

    auto g_xyy_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 94);

    auto g_xyy_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 95);

    auto g_xyy_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 96);

    auto g_xyy_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 97);

    auto g_xyy_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 98);

    auto g_xyy_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 99);

    auto g_xyy_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 100);

    auto g_xyy_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 101);

    auto g_xyy_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 102);

    auto g_xyy_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 103);

    auto g_xyy_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 104);

    auto g_xyy_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 105);

    auto g_xyy_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 106);

    auto g_xyy_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 107);

    auto g_xyy_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 108);

    auto g_xyy_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 109);

    auto g_xyy_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 110);

    auto g_xyy_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 111);

    #pragma omp simd aligned(g_xyy_xxxxxx_0, g_xyy_xxxxxy_0, g_xyy_xxxxxz_0, g_xyy_xxxxyy_0, g_xyy_xxxxyz_0, g_xyy_xxxxzz_0, g_xyy_xxxyyy_0, g_xyy_xxxyyz_0, g_xyy_xxxyzz_0, g_xyy_xxxzzz_0, g_xyy_xxyyyy_0, g_xyy_xxyyyz_0, g_xyy_xxyyzz_0, g_xyy_xxyzzz_0, g_xyy_xxzzzz_0, g_xyy_xyyyyy_0, g_xyy_xyyyyz_0, g_xyy_xyyyzz_0, g_xyy_xyyzzz_0, g_xyy_xyzzzz_0, g_xyy_xzzzzz_0, g_xyy_yyyyyy_0, g_xyy_yyyyyz_0, g_xyy_yyyyzz_0, g_xyy_yyyzzz_0, g_xyy_yyzzzz_0, g_xyy_yzzzzz_0, g_xyy_zzzzzz_0, g_yy_xxxxx_1, g_yy_xxxxxx_1, g_yy_xxxxxy_1, g_yy_xxxxxz_1, g_yy_xxxxy_1, g_yy_xxxxyy_1, g_yy_xxxxyz_1, g_yy_xxxxz_1, g_yy_xxxxzz_1, g_yy_xxxyy_1, g_yy_xxxyyy_1, g_yy_xxxyyz_1, g_yy_xxxyz_1, g_yy_xxxyzz_1, g_yy_xxxzz_1, g_yy_xxxzzz_1, g_yy_xxyyy_1, g_yy_xxyyyy_1, g_yy_xxyyyz_1, g_yy_xxyyz_1, g_yy_xxyyzz_1, g_yy_xxyzz_1, g_yy_xxyzzz_1, g_yy_xxzzz_1, g_yy_xxzzzz_1, g_yy_xyyyy_1, g_yy_xyyyyy_1, g_yy_xyyyyz_1, g_yy_xyyyz_1, g_yy_xyyyzz_1, g_yy_xyyzz_1, g_yy_xyyzzz_1, g_yy_xyzzz_1, g_yy_xyzzzz_1, g_yy_xzzzz_1, g_yy_xzzzzz_1, g_yy_yyyyy_1, g_yy_yyyyyy_1, g_yy_yyyyyz_1, g_yy_yyyyz_1, g_yy_yyyyzz_1, g_yy_yyyzz_1, g_yy_yyyzzz_1, g_yy_yyzzz_1, g_yy_yyzzzz_1, g_yy_yzzzz_1, g_yy_yzzzzz_1, g_yy_zzzzz_1, g_yy_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyy_xxxxxx_0[i] = 6.0 * g_yy_xxxxx_1[i] * fe_0 + g_yy_xxxxxx_1[i] * pa_x[i];

        g_xyy_xxxxxy_0[i] = 5.0 * g_yy_xxxxy_1[i] * fe_0 + g_yy_xxxxxy_1[i] * pa_x[i];

        g_xyy_xxxxxz_0[i] = 5.0 * g_yy_xxxxz_1[i] * fe_0 + g_yy_xxxxxz_1[i] * pa_x[i];

        g_xyy_xxxxyy_0[i] = 4.0 * g_yy_xxxyy_1[i] * fe_0 + g_yy_xxxxyy_1[i] * pa_x[i];

        g_xyy_xxxxyz_0[i] = 4.0 * g_yy_xxxyz_1[i] * fe_0 + g_yy_xxxxyz_1[i] * pa_x[i];

        g_xyy_xxxxzz_0[i] = 4.0 * g_yy_xxxzz_1[i] * fe_0 + g_yy_xxxxzz_1[i] * pa_x[i];

        g_xyy_xxxyyy_0[i] = 3.0 * g_yy_xxyyy_1[i] * fe_0 + g_yy_xxxyyy_1[i] * pa_x[i];

        g_xyy_xxxyyz_0[i] = 3.0 * g_yy_xxyyz_1[i] * fe_0 + g_yy_xxxyyz_1[i] * pa_x[i];

        g_xyy_xxxyzz_0[i] = 3.0 * g_yy_xxyzz_1[i] * fe_0 + g_yy_xxxyzz_1[i] * pa_x[i];

        g_xyy_xxxzzz_0[i] = 3.0 * g_yy_xxzzz_1[i] * fe_0 + g_yy_xxxzzz_1[i] * pa_x[i];

        g_xyy_xxyyyy_0[i] = 2.0 * g_yy_xyyyy_1[i] * fe_0 + g_yy_xxyyyy_1[i] * pa_x[i];

        g_xyy_xxyyyz_0[i] = 2.0 * g_yy_xyyyz_1[i] * fe_0 + g_yy_xxyyyz_1[i] * pa_x[i];

        g_xyy_xxyyzz_0[i] = 2.0 * g_yy_xyyzz_1[i] * fe_0 + g_yy_xxyyzz_1[i] * pa_x[i];

        g_xyy_xxyzzz_0[i] = 2.0 * g_yy_xyzzz_1[i] * fe_0 + g_yy_xxyzzz_1[i] * pa_x[i];

        g_xyy_xxzzzz_0[i] = 2.0 * g_yy_xzzzz_1[i] * fe_0 + g_yy_xxzzzz_1[i] * pa_x[i];

        g_xyy_xyyyyy_0[i] = g_yy_yyyyy_1[i] * fe_0 + g_yy_xyyyyy_1[i] * pa_x[i];

        g_xyy_xyyyyz_0[i] = g_yy_yyyyz_1[i] * fe_0 + g_yy_xyyyyz_1[i] * pa_x[i];

        g_xyy_xyyyzz_0[i] = g_yy_yyyzz_1[i] * fe_0 + g_yy_xyyyzz_1[i] * pa_x[i];

        g_xyy_xyyzzz_0[i] = g_yy_yyzzz_1[i] * fe_0 + g_yy_xyyzzz_1[i] * pa_x[i];

        g_xyy_xyzzzz_0[i] = g_yy_yzzzz_1[i] * fe_0 + g_yy_xyzzzz_1[i] * pa_x[i];

        g_xyy_xzzzzz_0[i] = g_yy_zzzzz_1[i] * fe_0 + g_yy_xzzzzz_1[i] * pa_x[i];

        g_xyy_yyyyyy_0[i] = g_yy_yyyyyy_1[i] * pa_x[i];

        g_xyy_yyyyyz_0[i] = g_yy_yyyyyz_1[i] * pa_x[i];

        g_xyy_yyyyzz_0[i] = g_yy_yyyyzz_1[i] * pa_x[i];

        g_xyy_yyyzzz_0[i] = g_yy_yyyzzz_1[i] * pa_x[i];

        g_xyy_yyzzzz_0[i] = g_yy_yyzzzz_1[i] * pa_x[i];

        g_xyy_yzzzzz_0[i] = g_yy_yzzzzz_1[i] * pa_x[i];

        g_xyy_zzzzzz_0[i] = g_yy_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : FI

    auto g_xyz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 112);

    auto g_xyz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 113);

    auto g_xyz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 114);

    auto g_xyz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 115);

    auto g_xyz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 116);

    auto g_xyz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 117);

    auto g_xyz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 118);

    auto g_xyz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 119);

    auto g_xyz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 120);

    auto g_xyz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 121);

    auto g_xyz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 122);

    auto g_xyz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 123);

    auto g_xyz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 124);

    auto g_xyz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 125);

    auto g_xyz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 126);

    auto g_xyz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 127);

    auto g_xyz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 128);

    auto g_xyz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 129);

    auto g_xyz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 130);

    auto g_xyz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 131);

    auto g_xyz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 132);

    auto g_xyz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 133);

    auto g_xyz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 134);

    auto g_xyz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 135);

    auto g_xyz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 136);

    auto g_xyz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 137);

    auto g_xyz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 138);

    auto g_xyz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 139);

    #pragma omp simd aligned(g_xy_xxxxxy_1, g_xy_xxxxyy_1, g_xy_xxxyyy_1, g_xy_xxyyyy_1, g_xy_xyyyyy_1, g_xyz_xxxxxx_0, g_xyz_xxxxxy_0, g_xyz_xxxxxz_0, g_xyz_xxxxyy_0, g_xyz_xxxxyz_0, g_xyz_xxxxzz_0, g_xyz_xxxyyy_0, g_xyz_xxxyyz_0, g_xyz_xxxyzz_0, g_xyz_xxxzzz_0, g_xyz_xxyyyy_0, g_xyz_xxyyyz_0, g_xyz_xxyyzz_0, g_xyz_xxyzzz_0, g_xyz_xxzzzz_0, g_xyz_xyyyyy_0, g_xyz_xyyyyz_0, g_xyz_xyyyzz_0, g_xyz_xyyzzz_0, g_xyz_xyzzzz_0, g_xyz_xzzzzz_0, g_xyz_yyyyyy_0, g_xyz_yyyyyz_0, g_xyz_yyyyzz_0, g_xyz_yyyzzz_0, g_xyz_yyzzzz_0, g_xyz_yzzzzz_0, g_xyz_zzzzzz_0, g_xz_xxxxxx_1, g_xz_xxxxxz_1, g_xz_xxxxzz_1, g_xz_xxxzzz_1, g_xz_xxzzzz_1, g_xz_xzzzzz_1, g_yz_xxxxyz_1, g_yz_xxxyyz_1, g_yz_xxxyz_1, g_yz_xxxyzz_1, g_yz_xxyyyz_1, g_yz_xxyyz_1, g_yz_xxyyzz_1, g_yz_xxyzz_1, g_yz_xxyzzz_1, g_yz_xyyyyz_1, g_yz_xyyyz_1, g_yz_xyyyzz_1, g_yz_xyyzz_1, g_yz_xyyzzz_1, g_yz_xyzzz_1, g_yz_xyzzzz_1, g_yz_yyyyyy_1, g_yz_yyyyyz_1, g_yz_yyyyz_1, g_yz_yyyyzz_1, g_yz_yyyzz_1, g_yz_yyyzzz_1, g_yz_yyzzz_1, g_yz_yyzzzz_1, g_yz_yzzzz_1, g_yz_yzzzzz_1, g_yz_zzzzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyz_xxxxxx_0[i] = g_xz_xxxxxx_1[i] * pa_y[i];

        g_xyz_xxxxxy_0[i] = g_xy_xxxxxy_1[i] * pa_z[i];

        g_xyz_xxxxxz_0[i] = g_xz_xxxxxz_1[i] * pa_y[i];

        g_xyz_xxxxyy_0[i] = g_xy_xxxxyy_1[i] * pa_z[i];

        g_xyz_xxxxyz_0[i] = 4.0 * g_yz_xxxyz_1[i] * fe_0 + g_yz_xxxxyz_1[i] * pa_x[i];

        g_xyz_xxxxzz_0[i] = g_xz_xxxxzz_1[i] * pa_y[i];

        g_xyz_xxxyyy_0[i] = g_xy_xxxyyy_1[i] * pa_z[i];

        g_xyz_xxxyyz_0[i] = 3.0 * g_yz_xxyyz_1[i] * fe_0 + g_yz_xxxyyz_1[i] * pa_x[i];

        g_xyz_xxxyzz_0[i] = 3.0 * g_yz_xxyzz_1[i] * fe_0 + g_yz_xxxyzz_1[i] * pa_x[i];

        g_xyz_xxxzzz_0[i] = g_xz_xxxzzz_1[i] * pa_y[i];

        g_xyz_xxyyyy_0[i] = g_xy_xxyyyy_1[i] * pa_z[i];

        g_xyz_xxyyyz_0[i] = 2.0 * g_yz_xyyyz_1[i] * fe_0 + g_yz_xxyyyz_1[i] * pa_x[i];

        g_xyz_xxyyzz_0[i] = 2.0 * g_yz_xyyzz_1[i] * fe_0 + g_yz_xxyyzz_1[i] * pa_x[i];

        g_xyz_xxyzzz_0[i] = 2.0 * g_yz_xyzzz_1[i] * fe_0 + g_yz_xxyzzz_1[i] * pa_x[i];

        g_xyz_xxzzzz_0[i] = g_xz_xxzzzz_1[i] * pa_y[i];

        g_xyz_xyyyyy_0[i] = g_xy_xyyyyy_1[i] * pa_z[i];

        g_xyz_xyyyyz_0[i] = g_yz_yyyyz_1[i] * fe_0 + g_yz_xyyyyz_1[i] * pa_x[i];

        g_xyz_xyyyzz_0[i] = g_yz_yyyzz_1[i] * fe_0 + g_yz_xyyyzz_1[i] * pa_x[i];

        g_xyz_xyyzzz_0[i] = g_yz_yyzzz_1[i] * fe_0 + g_yz_xyyzzz_1[i] * pa_x[i];

        g_xyz_xyzzzz_0[i] = g_yz_yzzzz_1[i] * fe_0 + g_yz_xyzzzz_1[i] * pa_x[i];

        g_xyz_xzzzzz_0[i] = g_xz_xzzzzz_1[i] * pa_y[i];

        g_xyz_yyyyyy_0[i] = g_yz_yyyyyy_1[i] * pa_x[i];

        g_xyz_yyyyyz_0[i] = g_yz_yyyyyz_1[i] * pa_x[i];

        g_xyz_yyyyzz_0[i] = g_yz_yyyyzz_1[i] * pa_x[i];

        g_xyz_yyyzzz_0[i] = g_yz_yyyzzz_1[i] * pa_x[i];

        g_xyz_yyzzzz_0[i] = g_yz_yyzzzz_1[i] * pa_x[i];

        g_xyz_yzzzzz_0[i] = g_yz_yzzzzz_1[i] * pa_x[i];

        g_xyz_zzzzzz_0[i] = g_yz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 140-168 components of targeted buffer : FI

    auto g_xzz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 140);

    auto g_xzz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 141);

    auto g_xzz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 142);

    auto g_xzz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 143);

    auto g_xzz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 144);

    auto g_xzz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 145);

    auto g_xzz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 146);

    auto g_xzz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 147);

    auto g_xzz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 148);

    auto g_xzz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 149);

    auto g_xzz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 150);

    auto g_xzz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 151);

    auto g_xzz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 152);

    auto g_xzz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 153);

    auto g_xzz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 154);

    auto g_xzz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 155);

    auto g_xzz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 156);

    auto g_xzz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 157);

    auto g_xzz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 158);

    auto g_xzz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 159);

    auto g_xzz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 160);

    auto g_xzz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 161);

    auto g_xzz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 162);

    auto g_xzz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 163);

    auto g_xzz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 164);

    auto g_xzz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 165);

    auto g_xzz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 166);

    auto g_xzz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 167);

    #pragma omp simd aligned(g_xzz_xxxxxx_0, g_xzz_xxxxxy_0, g_xzz_xxxxxz_0, g_xzz_xxxxyy_0, g_xzz_xxxxyz_0, g_xzz_xxxxzz_0, g_xzz_xxxyyy_0, g_xzz_xxxyyz_0, g_xzz_xxxyzz_0, g_xzz_xxxzzz_0, g_xzz_xxyyyy_0, g_xzz_xxyyyz_0, g_xzz_xxyyzz_0, g_xzz_xxyzzz_0, g_xzz_xxzzzz_0, g_xzz_xyyyyy_0, g_xzz_xyyyyz_0, g_xzz_xyyyzz_0, g_xzz_xyyzzz_0, g_xzz_xyzzzz_0, g_xzz_xzzzzz_0, g_xzz_yyyyyy_0, g_xzz_yyyyyz_0, g_xzz_yyyyzz_0, g_xzz_yyyzzz_0, g_xzz_yyzzzz_0, g_xzz_yzzzzz_0, g_xzz_zzzzzz_0, g_zz_xxxxx_1, g_zz_xxxxxx_1, g_zz_xxxxxy_1, g_zz_xxxxxz_1, g_zz_xxxxy_1, g_zz_xxxxyy_1, g_zz_xxxxyz_1, g_zz_xxxxz_1, g_zz_xxxxzz_1, g_zz_xxxyy_1, g_zz_xxxyyy_1, g_zz_xxxyyz_1, g_zz_xxxyz_1, g_zz_xxxyzz_1, g_zz_xxxzz_1, g_zz_xxxzzz_1, g_zz_xxyyy_1, g_zz_xxyyyy_1, g_zz_xxyyyz_1, g_zz_xxyyz_1, g_zz_xxyyzz_1, g_zz_xxyzz_1, g_zz_xxyzzz_1, g_zz_xxzzz_1, g_zz_xxzzzz_1, g_zz_xyyyy_1, g_zz_xyyyyy_1, g_zz_xyyyyz_1, g_zz_xyyyz_1, g_zz_xyyyzz_1, g_zz_xyyzz_1, g_zz_xyyzzz_1, g_zz_xyzzz_1, g_zz_xyzzzz_1, g_zz_xzzzz_1, g_zz_xzzzzz_1, g_zz_yyyyy_1, g_zz_yyyyyy_1, g_zz_yyyyyz_1, g_zz_yyyyz_1, g_zz_yyyyzz_1, g_zz_yyyzz_1, g_zz_yyyzzz_1, g_zz_yyzzz_1, g_zz_yyzzzz_1, g_zz_yzzzz_1, g_zz_yzzzzz_1, g_zz_zzzzz_1, g_zz_zzzzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzz_xxxxxx_0[i] = 6.0 * g_zz_xxxxx_1[i] * fe_0 + g_zz_xxxxxx_1[i] * pa_x[i];

        g_xzz_xxxxxy_0[i] = 5.0 * g_zz_xxxxy_1[i] * fe_0 + g_zz_xxxxxy_1[i] * pa_x[i];

        g_xzz_xxxxxz_0[i] = 5.0 * g_zz_xxxxz_1[i] * fe_0 + g_zz_xxxxxz_1[i] * pa_x[i];

        g_xzz_xxxxyy_0[i] = 4.0 * g_zz_xxxyy_1[i] * fe_0 + g_zz_xxxxyy_1[i] * pa_x[i];

        g_xzz_xxxxyz_0[i] = 4.0 * g_zz_xxxyz_1[i] * fe_0 + g_zz_xxxxyz_1[i] * pa_x[i];

        g_xzz_xxxxzz_0[i] = 4.0 * g_zz_xxxzz_1[i] * fe_0 + g_zz_xxxxzz_1[i] * pa_x[i];

        g_xzz_xxxyyy_0[i] = 3.0 * g_zz_xxyyy_1[i] * fe_0 + g_zz_xxxyyy_1[i] * pa_x[i];

        g_xzz_xxxyyz_0[i] = 3.0 * g_zz_xxyyz_1[i] * fe_0 + g_zz_xxxyyz_1[i] * pa_x[i];

        g_xzz_xxxyzz_0[i] = 3.0 * g_zz_xxyzz_1[i] * fe_0 + g_zz_xxxyzz_1[i] * pa_x[i];

        g_xzz_xxxzzz_0[i] = 3.0 * g_zz_xxzzz_1[i] * fe_0 + g_zz_xxxzzz_1[i] * pa_x[i];

        g_xzz_xxyyyy_0[i] = 2.0 * g_zz_xyyyy_1[i] * fe_0 + g_zz_xxyyyy_1[i] * pa_x[i];

        g_xzz_xxyyyz_0[i] = 2.0 * g_zz_xyyyz_1[i] * fe_0 + g_zz_xxyyyz_1[i] * pa_x[i];

        g_xzz_xxyyzz_0[i] = 2.0 * g_zz_xyyzz_1[i] * fe_0 + g_zz_xxyyzz_1[i] * pa_x[i];

        g_xzz_xxyzzz_0[i] = 2.0 * g_zz_xyzzz_1[i] * fe_0 + g_zz_xxyzzz_1[i] * pa_x[i];

        g_xzz_xxzzzz_0[i] = 2.0 * g_zz_xzzzz_1[i] * fe_0 + g_zz_xxzzzz_1[i] * pa_x[i];

        g_xzz_xyyyyy_0[i] = g_zz_yyyyy_1[i] * fe_0 + g_zz_xyyyyy_1[i] * pa_x[i];

        g_xzz_xyyyyz_0[i] = g_zz_yyyyz_1[i] * fe_0 + g_zz_xyyyyz_1[i] * pa_x[i];

        g_xzz_xyyyzz_0[i] = g_zz_yyyzz_1[i] * fe_0 + g_zz_xyyyzz_1[i] * pa_x[i];

        g_xzz_xyyzzz_0[i] = g_zz_yyzzz_1[i] * fe_0 + g_zz_xyyzzz_1[i] * pa_x[i];

        g_xzz_xyzzzz_0[i] = g_zz_yzzzz_1[i] * fe_0 + g_zz_xyzzzz_1[i] * pa_x[i];

        g_xzz_xzzzzz_0[i] = g_zz_zzzzz_1[i] * fe_0 + g_zz_xzzzzz_1[i] * pa_x[i];

        g_xzz_yyyyyy_0[i] = g_zz_yyyyyy_1[i] * pa_x[i];

        g_xzz_yyyyyz_0[i] = g_zz_yyyyyz_1[i] * pa_x[i];

        g_xzz_yyyyzz_0[i] = g_zz_yyyyzz_1[i] * pa_x[i];

        g_xzz_yyyzzz_0[i] = g_zz_yyyzzz_1[i] * pa_x[i];

        g_xzz_yyzzzz_0[i] = g_zz_yyzzzz_1[i] * pa_x[i];

        g_xzz_yzzzzz_0[i] = g_zz_yzzzzz_1[i] * pa_x[i];

        g_xzz_zzzzzz_0[i] = g_zz_zzzzzz_1[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : FI

    auto g_yyy_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 168);

    auto g_yyy_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 169);

    auto g_yyy_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 170);

    auto g_yyy_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 171);

    auto g_yyy_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 172);

    auto g_yyy_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 173);

    auto g_yyy_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 174);

    auto g_yyy_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 175);

    auto g_yyy_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 176);

    auto g_yyy_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 177);

    auto g_yyy_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 178);

    auto g_yyy_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 179);

    auto g_yyy_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 180);

    auto g_yyy_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 181);

    auto g_yyy_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 182);

    auto g_yyy_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 183);

    auto g_yyy_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 184);

    auto g_yyy_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 185);

    auto g_yyy_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 186);

    auto g_yyy_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 187);

    auto g_yyy_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 188);

    auto g_yyy_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 189);

    auto g_yyy_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 190);

    auto g_yyy_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 191);

    auto g_yyy_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 192);

    auto g_yyy_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 193);

    auto g_yyy_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 194);

    auto g_yyy_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 195);

    #pragma omp simd aligned(g_y_xxxxxx_0, g_y_xxxxxx_1, g_y_xxxxxy_0, g_y_xxxxxy_1, g_y_xxxxxz_0, g_y_xxxxxz_1, g_y_xxxxyy_0, g_y_xxxxyy_1, g_y_xxxxyz_0, g_y_xxxxyz_1, g_y_xxxxzz_0, g_y_xxxxzz_1, g_y_xxxyyy_0, g_y_xxxyyy_1, g_y_xxxyyz_0, g_y_xxxyyz_1, g_y_xxxyzz_0, g_y_xxxyzz_1, g_y_xxxzzz_0, g_y_xxxzzz_1, g_y_xxyyyy_0, g_y_xxyyyy_1, g_y_xxyyyz_0, g_y_xxyyyz_1, g_y_xxyyzz_0, g_y_xxyyzz_1, g_y_xxyzzz_0, g_y_xxyzzz_1, g_y_xxzzzz_0, g_y_xxzzzz_1, g_y_xyyyyy_0, g_y_xyyyyy_1, g_y_xyyyyz_0, g_y_xyyyyz_1, g_y_xyyyzz_0, g_y_xyyyzz_1, g_y_xyyzzz_0, g_y_xyyzzz_1, g_y_xyzzzz_0, g_y_xyzzzz_1, g_y_xzzzzz_0, g_y_xzzzzz_1, g_y_yyyyyy_0, g_y_yyyyyy_1, g_y_yyyyyz_0, g_y_yyyyyz_1, g_y_yyyyzz_0, g_y_yyyyzz_1, g_y_yyyzzz_0, g_y_yyyzzz_1, g_y_yyzzzz_0, g_y_yyzzzz_1, g_y_yzzzzz_0, g_y_yzzzzz_1, g_y_zzzzzz_0, g_y_zzzzzz_1, g_yy_xxxxx_1, g_yy_xxxxxx_1, g_yy_xxxxxy_1, g_yy_xxxxxz_1, g_yy_xxxxy_1, g_yy_xxxxyy_1, g_yy_xxxxyz_1, g_yy_xxxxz_1, g_yy_xxxxzz_1, g_yy_xxxyy_1, g_yy_xxxyyy_1, g_yy_xxxyyz_1, g_yy_xxxyz_1, g_yy_xxxyzz_1, g_yy_xxxzz_1, g_yy_xxxzzz_1, g_yy_xxyyy_1, g_yy_xxyyyy_1, g_yy_xxyyyz_1, g_yy_xxyyz_1, g_yy_xxyyzz_1, g_yy_xxyzz_1, g_yy_xxyzzz_1, g_yy_xxzzz_1, g_yy_xxzzzz_1, g_yy_xyyyy_1, g_yy_xyyyyy_1, g_yy_xyyyyz_1, g_yy_xyyyz_1, g_yy_xyyyzz_1, g_yy_xyyzz_1, g_yy_xyyzzz_1, g_yy_xyzzz_1, g_yy_xyzzzz_1, g_yy_xzzzz_1, g_yy_xzzzzz_1, g_yy_yyyyy_1, g_yy_yyyyyy_1, g_yy_yyyyyz_1, g_yy_yyyyz_1, g_yy_yyyyzz_1, g_yy_yyyzz_1, g_yy_yyyzzz_1, g_yy_yyzzz_1, g_yy_yyzzzz_1, g_yy_yzzzz_1, g_yy_yzzzzz_1, g_yy_zzzzz_1, g_yy_zzzzzz_1, g_yyy_xxxxxx_0, g_yyy_xxxxxy_0, g_yyy_xxxxxz_0, g_yyy_xxxxyy_0, g_yyy_xxxxyz_0, g_yyy_xxxxzz_0, g_yyy_xxxyyy_0, g_yyy_xxxyyz_0, g_yyy_xxxyzz_0, g_yyy_xxxzzz_0, g_yyy_xxyyyy_0, g_yyy_xxyyyz_0, g_yyy_xxyyzz_0, g_yyy_xxyzzz_0, g_yyy_xxzzzz_0, g_yyy_xyyyyy_0, g_yyy_xyyyyz_0, g_yyy_xyyyzz_0, g_yyy_xyyzzz_0, g_yyy_xyzzzz_0, g_yyy_xzzzzz_0, g_yyy_yyyyyy_0, g_yyy_yyyyyz_0, g_yyy_yyyyzz_0, g_yyy_yyyzzz_0, g_yyy_yyzzzz_0, g_yyy_yzzzzz_0, g_yyy_zzzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyy_xxxxxx_0[i] = 2.0 * g_y_xxxxxx_0[i] * fbe_0 - 2.0 * g_y_xxxxxx_1[i] * fz_be_0 + g_yy_xxxxxx_1[i] * pa_y[i];

        g_yyy_xxxxxy_0[i] = 2.0 * g_y_xxxxxy_0[i] * fbe_0 - 2.0 * g_y_xxxxxy_1[i] * fz_be_0 + g_yy_xxxxx_1[i] * fe_0 + g_yy_xxxxxy_1[i] * pa_y[i];

        g_yyy_xxxxxz_0[i] = 2.0 * g_y_xxxxxz_0[i] * fbe_0 - 2.0 * g_y_xxxxxz_1[i] * fz_be_0 + g_yy_xxxxxz_1[i] * pa_y[i];

        g_yyy_xxxxyy_0[i] = 2.0 * g_y_xxxxyy_0[i] * fbe_0 - 2.0 * g_y_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yy_xxxxy_1[i] * fe_0 + g_yy_xxxxyy_1[i] * pa_y[i];

        g_yyy_xxxxyz_0[i] = 2.0 * g_y_xxxxyz_0[i] * fbe_0 - 2.0 * g_y_xxxxyz_1[i] * fz_be_0 + g_yy_xxxxz_1[i] * fe_0 + g_yy_xxxxyz_1[i] * pa_y[i];

        g_yyy_xxxxzz_0[i] = 2.0 * g_y_xxxxzz_0[i] * fbe_0 - 2.0 * g_y_xxxxzz_1[i] * fz_be_0 + g_yy_xxxxzz_1[i] * pa_y[i];

        g_yyy_xxxyyy_0[i] = 2.0 * g_y_xxxyyy_0[i] * fbe_0 - 2.0 * g_y_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yy_xxxyy_1[i] * fe_0 + g_yy_xxxyyy_1[i] * pa_y[i];

        g_yyy_xxxyyz_0[i] = 2.0 * g_y_xxxyyz_0[i] * fbe_0 - 2.0 * g_y_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yy_xxxyz_1[i] * fe_0 + g_yy_xxxyyz_1[i] * pa_y[i];

        g_yyy_xxxyzz_0[i] = 2.0 * g_y_xxxyzz_0[i] * fbe_0 - 2.0 * g_y_xxxyzz_1[i] * fz_be_0 + g_yy_xxxzz_1[i] * fe_0 + g_yy_xxxyzz_1[i] * pa_y[i];

        g_yyy_xxxzzz_0[i] = 2.0 * g_y_xxxzzz_0[i] * fbe_0 - 2.0 * g_y_xxxzzz_1[i] * fz_be_0 + g_yy_xxxzzz_1[i] * pa_y[i];

        g_yyy_xxyyyy_0[i] = 2.0 * g_y_xxyyyy_0[i] * fbe_0 - 2.0 * g_y_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yy_xxyyy_1[i] * fe_0 + g_yy_xxyyyy_1[i] * pa_y[i];

        g_yyy_xxyyyz_0[i] = 2.0 * g_y_xxyyyz_0[i] * fbe_0 - 2.0 * g_y_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yy_xxyyz_1[i] * fe_0 + g_yy_xxyyyz_1[i] * pa_y[i];

        g_yyy_xxyyzz_0[i] = 2.0 * g_y_xxyyzz_0[i] * fbe_0 - 2.0 * g_y_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yy_xxyzz_1[i] * fe_0 + g_yy_xxyyzz_1[i] * pa_y[i];

        g_yyy_xxyzzz_0[i] = 2.0 * g_y_xxyzzz_0[i] * fbe_0 - 2.0 * g_y_xxyzzz_1[i] * fz_be_0 + g_yy_xxzzz_1[i] * fe_0 + g_yy_xxyzzz_1[i] * pa_y[i];

        g_yyy_xxzzzz_0[i] = 2.0 * g_y_xxzzzz_0[i] * fbe_0 - 2.0 * g_y_xxzzzz_1[i] * fz_be_0 + g_yy_xxzzzz_1[i] * pa_y[i];

        g_yyy_xyyyyy_0[i] = 2.0 * g_y_xyyyyy_0[i] * fbe_0 - 2.0 * g_y_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yy_xyyyy_1[i] * fe_0 + g_yy_xyyyyy_1[i] * pa_y[i];

        g_yyy_xyyyyz_0[i] = 2.0 * g_y_xyyyyz_0[i] * fbe_0 - 2.0 * g_y_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yy_xyyyz_1[i] * fe_0 + g_yy_xyyyyz_1[i] * pa_y[i];

        g_yyy_xyyyzz_0[i] = 2.0 * g_y_xyyyzz_0[i] * fbe_0 - 2.0 * g_y_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yy_xyyzz_1[i] * fe_0 + g_yy_xyyyzz_1[i] * pa_y[i];

        g_yyy_xyyzzz_0[i] = 2.0 * g_y_xyyzzz_0[i] * fbe_0 - 2.0 * g_y_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yy_xyzzz_1[i] * fe_0 + g_yy_xyyzzz_1[i] * pa_y[i];

        g_yyy_xyzzzz_0[i] = 2.0 * g_y_xyzzzz_0[i] * fbe_0 - 2.0 * g_y_xyzzzz_1[i] * fz_be_0 + g_yy_xzzzz_1[i] * fe_0 + g_yy_xyzzzz_1[i] * pa_y[i];

        g_yyy_xzzzzz_0[i] = 2.0 * g_y_xzzzzz_0[i] * fbe_0 - 2.0 * g_y_xzzzzz_1[i] * fz_be_0 + g_yy_xzzzzz_1[i] * pa_y[i];

        g_yyy_yyyyyy_0[i] = 2.0 * g_y_yyyyyy_0[i] * fbe_0 - 2.0 * g_y_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yy_yyyyy_1[i] * fe_0 + g_yy_yyyyyy_1[i] * pa_y[i];

        g_yyy_yyyyyz_0[i] = 2.0 * g_y_yyyyyz_0[i] * fbe_0 - 2.0 * g_y_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yy_yyyyz_1[i] * fe_0 + g_yy_yyyyyz_1[i] * pa_y[i];

        g_yyy_yyyyzz_0[i] = 2.0 * g_y_yyyyzz_0[i] * fbe_0 - 2.0 * g_y_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yy_yyyzz_1[i] * fe_0 + g_yy_yyyyzz_1[i] * pa_y[i];

        g_yyy_yyyzzz_0[i] = 2.0 * g_y_yyyzzz_0[i] * fbe_0 - 2.0 * g_y_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yy_yyzzz_1[i] * fe_0 + g_yy_yyyzzz_1[i] * pa_y[i];

        g_yyy_yyzzzz_0[i] = 2.0 * g_y_yyzzzz_0[i] * fbe_0 - 2.0 * g_y_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yy_yzzzz_1[i] * fe_0 + g_yy_yyzzzz_1[i] * pa_y[i];

        g_yyy_yzzzzz_0[i] = 2.0 * g_y_yzzzzz_0[i] * fbe_0 - 2.0 * g_y_yzzzzz_1[i] * fz_be_0 + g_yy_zzzzz_1[i] * fe_0 + g_yy_yzzzzz_1[i] * pa_y[i];

        g_yyy_zzzzzz_0[i] = 2.0 * g_y_zzzzzz_0[i] * fbe_0 - 2.0 * g_y_zzzzzz_1[i] * fz_be_0 + g_yy_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 196-224 components of targeted buffer : FI

    auto g_yyz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 196);

    auto g_yyz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 197);

    auto g_yyz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 198);

    auto g_yyz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 199);

    auto g_yyz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 200);

    auto g_yyz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 201);

    auto g_yyz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 202);

    auto g_yyz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 203);

    auto g_yyz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 204);

    auto g_yyz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 205);

    auto g_yyz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 206);

    auto g_yyz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 207);

    auto g_yyz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 208);

    auto g_yyz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 209);

    auto g_yyz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 210);

    auto g_yyz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 211);

    auto g_yyz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 212);

    auto g_yyz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 213);

    auto g_yyz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 214);

    auto g_yyz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 215);

    auto g_yyz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 216);

    auto g_yyz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 217);

    auto g_yyz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 218);

    auto g_yyz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 219);

    auto g_yyz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 220);

    auto g_yyz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 221);

    auto g_yyz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 222);

    auto g_yyz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 223);

    #pragma omp simd aligned(g_yy_xxxxx_1, g_yy_xxxxxx_1, g_yy_xxxxxy_1, g_yy_xxxxxz_1, g_yy_xxxxy_1, g_yy_xxxxyy_1, g_yy_xxxxyz_1, g_yy_xxxxz_1, g_yy_xxxxzz_1, g_yy_xxxyy_1, g_yy_xxxyyy_1, g_yy_xxxyyz_1, g_yy_xxxyz_1, g_yy_xxxyzz_1, g_yy_xxxzz_1, g_yy_xxxzzz_1, g_yy_xxyyy_1, g_yy_xxyyyy_1, g_yy_xxyyyz_1, g_yy_xxyyz_1, g_yy_xxyyzz_1, g_yy_xxyzz_1, g_yy_xxyzzz_1, g_yy_xxzzz_1, g_yy_xxzzzz_1, g_yy_xyyyy_1, g_yy_xyyyyy_1, g_yy_xyyyyz_1, g_yy_xyyyz_1, g_yy_xyyyzz_1, g_yy_xyyzz_1, g_yy_xyyzzz_1, g_yy_xyzzz_1, g_yy_xyzzzz_1, g_yy_xzzzz_1, g_yy_xzzzzz_1, g_yy_yyyyy_1, g_yy_yyyyyy_1, g_yy_yyyyyz_1, g_yy_yyyyz_1, g_yy_yyyyzz_1, g_yy_yyyzz_1, g_yy_yyyzzz_1, g_yy_yyzzz_1, g_yy_yyzzzz_1, g_yy_yzzzz_1, g_yy_yzzzzz_1, g_yy_zzzzz_1, g_yy_zzzzzz_1, g_yyz_xxxxxx_0, g_yyz_xxxxxy_0, g_yyz_xxxxxz_0, g_yyz_xxxxyy_0, g_yyz_xxxxyz_0, g_yyz_xxxxzz_0, g_yyz_xxxyyy_0, g_yyz_xxxyyz_0, g_yyz_xxxyzz_0, g_yyz_xxxzzz_0, g_yyz_xxyyyy_0, g_yyz_xxyyyz_0, g_yyz_xxyyzz_0, g_yyz_xxyzzz_0, g_yyz_xxzzzz_0, g_yyz_xyyyyy_0, g_yyz_xyyyyz_0, g_yyz_xyyyzz_0, g_yyz_xyyzzz_0, g_yyz_xyzzzz_0, g_yyz_xzzzzz_0, g_yyz_yyyyyy_0, g_yyz_yyyyyz_0, g_yyz_yyyyzz_0, g_yyz_yyyzzz_0, g_yyz_yyzzzz_0, g_yyz_yzzzzz_0, g_yyz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyz_xxxxxx_0[i] = g_yy_xxxxxx_1[i] * pa_z[i];

        g_yyz_xxxxxy_0[i] = g_yy_xxxxxy_1[i] * pa_z[i];

        g_yyz_xxxxxz_0[i] = g_yy_xxxxx_1[i] * fe_0 + g_yy_xxxxxz_1[i] * pa_z[i];

        g_yyz_xxxxyy_0[i] = g_yy_xxxxyy_1[i] * pa_z[i];

        g_yyz_xxxxyz_0[i] = g_yy_xxxxy_1[i] * fe_0 + g_yy_xxxxyz_1[i] * pa_z[i];

        g_yyz_xxxxzz_0[i] = 2.0 * g_yy_xxxxz_1[i] * fe_0 + g_yy_xxxxzz_1[i] * pa_z[i];

        g_yyz_xxxyyy_0[i] = g_yy_xxxyyy_1[i] * pa_z[i];

        g_yyz_xxxyyz_0[i] = g_yy_xxxyy_1[i] * fe_0 + g_yy_xxxyyz_1[i] * pa_z[i];

        g_yyz_xxxyzz_0[i] = 2.0 * g_yy_xxxyz_1[i] * fe_0 + g_yy_xxxyzz_1[i] * pa_z[i];

        g_yyz_xxxzzz_0[i] = 3.0 * g_yy_xxxzz_1[i] * fe_0 + g_yy_xxxzzz_1[i] * pa_z[i];

        g_yyz_xxyyyy_0[i] = g_yy_xxyyyy_1[i] * pa_z[i];

        g_yyz_xxyyyz_0[i] = g_yy_xxyyy_1[i] * fe_0 + g_yy_xxyyyz_1[i] * pa_z[i];

        g_yyz_xxyyzz_0[i] = 2.0 * g_yy_xxyyz_1[i] * fe_0 + g_yy_xxyyzz_1[i] * pa_z[i];

        g_yyz_xxyzzz_0[i] = 3.0 * g_yy_xxyzz_1[i] * fe_0 + g_yy_xxyzzz_1[i] * pa_z[i];

        g_yyz_xxzzzz_0[i] = 4.0 * g_yy_xxzzz_1[i] * fe_0 + g_yy_xxzzzz_1[i] * pa_z[i];

        g_yyz_xyyyyy_0[i] = g_yy_xyyyyy_1[i] * pa_z[i];

        g_yyz_xyyyyz_0[i] = g_yy_xyyyy_1[i] * fe_0 + g_yy_xyyyyz_1[i] * pa_z[i];

        g_yyz_xyyyzz_0[i] = 2.0 * g_yy_xyyyz_1[i] * fe_0 + g_yy_xyyyzz_1[i] * pa_z[i];

        g_yyz_xyyzzz_0[i] = 3.0 * g_yy_xyyzz_1[i] * fe_0 + g_yy_xyyzzz_1[i] * pa_z[i];

        g_yyz_xyzzzz_0[i] = 4.0 * g_yy_xyzzz_1[i] * fe_0 + g_yy_xyzzzz_1[i] * pa_z[i];

        g_yyz_xzzzzz_0[i] = 5.0 * g_yy_xzzzz_1[i] * fe_0 + g_yy_xzzzzz_1[i] * pa_z[i];

        g_yyz_yyyyyy_0[i] = g_yy_yyyyyy_1[i] * pa_z[i];

        g_yyz_yyyyyz_0[i] = g_yy_yyyyy_1[i] * fe_0 + g_yy_yyyyyz_1[i] * pa_z[i];

        g_yyz_yyyyzz_0[i] = 2.0 * g_yy_yyyyz_1[i] * fe_0 + g_yy_yyyyzz_1[i] * pa_z[i];

        g_yyz_yyyzzz_0[i] = 3.0 * g_yy_yyyzz_1[i] * fe_0 + g_yy_yyyzzz_1[i] * pa_z[i];

        g_yyz_yyzzzz_0[i] = 4.0 * g_yy_yyzzz_1[i] * fe_0 + g_yy_yyzzzz_1[i] * pa_z[i];

        g_yyz_yzzzzz_0[i] = 5.0 * g_yy_yzzzz_1[i] * fe_0 + g_yy_yzzzzz_1[i] * pa_z[i];

        g_yyz_zzzzzz_0[i] = 6.0 * g_yy_zzzzz_1[i] * fe_0 + g_yy_zzzzzz_1[i] * pa_z[i];
    }

    // Set up 224-252 components of targeted buffer : FI

    auto g_yzz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 224);

    auto g_yzz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 225);

    auto g_yzz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 226);

    auto g_yzz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 227);

    auto g_yzz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 228);

    auto g_yzz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 229);

    auto g_yzz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 230);

    auto g_yzz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 231);

    auto g_yzz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 232);

    auto g_yzz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 233);

    auto g_yzz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 234);

    auto g_yzz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 235);

    auto g_yzz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 236);

    auto g_yzz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 237);

    auto g_yzz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 238);

    auto g_yzz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 239);

    auto g_yzz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 240);

    auto g_yzz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 241);

    auto g_yzz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 242);

    auto g_yzz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 243);

    auto g_yzz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 244);

    auto g_yzz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 245);

    auto g_yzz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 246);

    auto g_yzz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 247);

    auto g_yzz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 248);

    auto g_yzz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 249);

    auto g_yzz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 250);

    auto g_yzz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 251);

    #pragma omp simd aligned(g_yzz_xxxxxx_0, g_yzz_xxxxxy_0, g_yzz_xxxxxz_0, g_yzz_xxxxyy_0, g_yzz_xxxxyz_0, g_yzz_xxxxzz_0, g_yzz_xxxyyy_0, g_yzz_xxxyyz_0, g_yzz_xxxyzz_0, g_yzz_xxxzzz_0, g_yzz_xxyyyy_0, g_yzz_xxyyyz_0, g_yzz_xxyyzz_0, g_yzz_xxyzzz_0, g_yzz_xxzzzz_0, g_yzz_xyyyyy_0, g_yzz_xyyyyz_0, g_yzz_xyyyzz_0, g_yzz_xyyzzz_0, g_yzz_xyzzzz_0, g_yzz_xzzzzz_0, g_yzz_yyyyyy_0, g_yzz_yyyyyz_0, g_yzz_yyyyzz_0, g_yzz_yyyzzz_0, g_yzz_yyzzzz_0, g_yzz_yzzzzz_0, g_yzz_zzzzzz_0, g_zz_xxxxx_1, g_zz_xxxxxx_1, g_zz_xxxxxy_1, g_zz_xxxxxz_1, g_zz_xxxxy_1, g_zz_xxxxyy_1, g_zz_xxxxyz_1, g_zz_xxxxz_1, g_zz_xxxxzz_1, g_zz_xxxyy_1, g_zz_xxxyyy_1, g_zz_xxxyyz_1, g_zz_xxxyz_1, g_zz_xxxyzz_1, g_zz_xxxzz_1, g_zz_xxxzzz_1, g_zz_xxyyy_1, g_zz_xxyyyy_1, g_zz_xxyyyz_1, g_zz_xxyyz_1, g_zz_xxyyzz_1, g_zz_xxyzz_1, g_zz_xxyzzz_1, g_zz_xxzzz_1, g_zz_xxzzzz_1, g_zz_xyyyy_1, g_zz_xyyyyy_1, g_zz_xyyyyz_1, g_zz_xyyyz_1, g_zz_xyyyzz_1, g_zz_xyyzz_1, g_zz_xyyzzz_1, g_zz_xyzzz_1, g_zz_xyzzzz_1, g_zz_xzzzz_1, g_zz_xzzzzz_1, g_zz_yyyyy_1, g_zz_yyyyyy_1, g_zz_yyyyyz_1, g_zz_yyyyz_1, g_zz_yyyyzz_1, g_zz_yyyzz_1, g_zz_yyyzzz_1, g_zz_yyzzz_1, g_zz_yyzzzz_1, g_zz_yzzzz_1, g_zz_yzzzzz_1, g_zz_zzzzz_1, g_zz_zzzzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzz_xxxxxx_0[i] = g_zz_xxxxxx_1[i] * pa_y[i];

        g_yzz_xxxxxy_0[i] = g_zz_xxxxx_1[i] * fe_0 + g_zz_xxxxxy_1[i] * pa_y[i];

        g_yzz_xxxxxz_0[i] = g_zz_xxxxxz_1[i] * pa_y[i];

        g_yzz_xxxxyy_0[i] = 2.0 * g_zz_xxxxy_1[i] * fe_0 + g_zz_xxxxyy_1[i] * pa_y[i];

        g_yzz_xxxxyz_0[i] = g_zz_xxxxz_1[i] * fe_0 + g_zz_xxxxyz_1[i] * pa_y[i];

        g_yzz_xxxxzz_0[i] = g_zz_xxxxzz_1[i] * pa_y[i];

        g_yzz_xxxyyy_0[i] = 3.0 * g_zz_xxxyy_1[i] * fe_0 + g_zz_xxxyyy_1[i] * pa_y[i];

        g_yzz_xxxyyz_0[i] = 2.0 * g_zz_xxxyz_1[i] * fe_0 + g_zz_xxxyyz_1[i] * pa_y[i];

        g_yzz_xxxyzz_0[i] = g_zz_xxxzz_1[i] * fe_0 + g_zz_xxxyzz_1[i] * pa_y[i];

        g_yzz_xxxzzz_0[i] = g_zz_xxxzzz_1[i] * pa_y[i];

        g_yzz_xxyyyy_0[i] = 4.0 * g_zz_xxyyy_1[i] * fe_0 + g_zz_xxyyyy_1[i] * pa_y[i];

        g_yzz_xxyyyz_0[i] = 3.0 * g_zz_xxyyz_1[i] * fe_0 + g_zz_xxyyyz_1[i] * pa_y[i];

        g_yzz_xxyyzz_0[i] = 2.0 * g_zz_xxyzz_1[i] * fe_0 + g_zz_xxyyzz_1[i] * pa_y[i];

        g_yzz_xxyzzz_0[i] = g_zz_xxzzz_1[i] * fe_0 + g_zz_xxyzzz_1[i] * pa_y[i];

        g_yzz_xxzzzz_0[i] = g_zz_xxzzzz_1[i] * pa_y[i];

        g_yzz_xyyyyy_0[i] = 5.0 * g_zz_xyyyy_1[i] * fe_0 + g_zz_xyyyyy_1[i] * pa_y[i];

        g_yzz_xyyyyz_0[i] = 4.0 * g_zz_xyyyz_1[i] * fe_0 + g_zz_xyyyyz_1[i] * pa_y[i];

        g_yzz_xyyyzz_0[i] = 3.0 * g_zz_xyyzz_1[i] * fe_0 + g_zz_xyyyzz_1[i] * pa_y[i];

        g_yzz_xyyzzz_0[i] = 2.0 * g_zz_xyzzz_1[i] * fe_0 + g_zz_xyyzzz_1[i] * pa_y[i];

        g_yzz_xyzzzz_0[i] = g_zz_xzzzz_1[i] * fe_0 + g_zz_xyzzzz_1[i] * pa_y[i];

        g_yzz_xzzzzz_0[i] = g_zz_xzzzzz_1[i] * pa_y[i];

        g_yzz_yyyyyy_0[i] = 6.0 * g_zz_yyyyy_1[i] * fe_0 + g_zz_yyyyyy_1[i] * pa_y[i];

        g_yzz_yyyyyz_0[i] = 5.0 * g_zz_yyyyz_1[i] * fe_0 + g_zz_yyyyyz_1[i] * pa_y[i];

        g_yzz_yyyyzz_0[i] = 4.0 * g_zz_yyyzz_1[i] * fe_0 + g_zz_yyyyzz_1[i] * pa_y[i];

        g_yzz_yyyzzz_0[i] = 3.0 * g_zz_yyzzz_1[i] * fe_0 + g_zz_yyyzzz_1[i] * pa_y[i];

        g_yzz_yyzzzz_0[i] = 2.0 * g_zz_yzzzz_1[i] * fe_0 + g_zz_yyzzzz_1[i] * pa_y[i];

        g_yzz_yzzzzz_0[i] = g_zz_zzzzz_1[i] * fe_0 + g_zz_yzzzzz_1[i] * pa_y[i];

        g_yzz_zzzzzz_0[i] = g_zz_zzzzzz_1[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : FI

    auto g_zzz_xxxxxx_0 = pbuffer.data(idx_eri_0_fi + 252);

    auto g_zzz_xxxxxy_0 = pbuffer.data(idx_eri_0_fi + 253);

    auto g_zzz_xxxxxz_0 = pbuffer.data(idx_eri_0_fi + 254);

    auto g_zzz_xxxxyy_0 = pbuffer.data(idx_eri_0_fi + 255);

    auto g_zzz_xxxxyz_0 = pbuffer.data(idx_eri_0_fi + 256);

    auto g_zzz_xxxxzz_0 = pbuffer.data(idx_eri_0_fi + 257);

    auto g_zzz_xxxyyy_0 = pbuffer.data(idx_eri_0_fi + 258);

    auto g_zzz_xxxyyz_0 = pbuffer.data(idx_eri_0_fi + 259);

    auto g_zzz_xxxyzz_0 = pbuffer.data(idx_eri_0_fi + 260);

    auto g_zzz_xxxzzz_0 = pbuffer.data(idx_eri_0_fi + 261);

    auto g_zzz_xxyyyy_0 = pbuffer.data(idx_eri_0_fi + 262);

    auto g_zzz_xxyyyz_0 = pbuffer.data(idx_eri_0_fi + 263);

    auto g_zzz_xxyyzz_0 = pbuffer.data(idx_eri_0_fi + 264);

    auto g_zzz_xxyzzz_0 = pbuffer.data(idx_eri_0_fi + 265);

    auto g_zzz_xxzzzz_0 = pbuffer.data(idx_eri_0_fi + 266);

    auto g_zzz_xyyyyy_0 = pbuffer.data(idx_eri_0_fi + 267);

    auto g_zzz_xyyyyz_0 = pbuffer.data(idx_eri_0_fi + 268);

    auto g_zzz_xyyyzz_0 = pbuffer.data(idx_eri_0_fi + 269);

    auto g_zzz_xyyzzz_0 = pbuffer.data(idx_eri_0_fi + 270);

    auto g_zzz_xyzzzz_0 = pbuffer.data(idx_eri_0_fi + 271);

    auto g_zzz_xzzzzz_0 = pbuffer.data(idx_eri_0_fi + 272);

    auto g_zzz_yyyyyy_0 = pbuffer.data(idx_eri_0_fi + 273);

    auto g_zzz_yyyyyz_0 = pbuffer.data(idx_eri_0_fi + 274);

    auto g_zzz_yyyyzz_0 = pbuffer.data(idx_eri_0_fi + 275);

    auto g_zzz_yyyzzz_0 = pbuffer.data(idx_eri_0_fi + 276);

    auto g_zzz_yyzzzz_0 = pbuffer.data(idx_eri_0_fi + 277);

    auto g_zzz_yzzzzz_0 = pbuffer.data(idx_eri_0_fi + 278);

    auto g_zzz_zzzzzz_0 = pbuffer.data(idx_eri_0_fi + 279);

    #pragma omp simd aligned(g_z_xxxxxx_0, g_z_xxxxxx_1, g_z_xxxxxy_0, g_z_xxxxxy_1, g_z_xxxxxz_0, g_z_xxxxxz_1, g_z_xxxxyy_0, g_z_xxxxyy_1, g_z_xxxxyz_0, g_z_xxxxyz_1, g_z_xxxxzz_0, g_z_xxxxzz_1, g_z_xxxyyy_0, g_z_xxxyyy_1, g_z_xxxyyz_0, g_z_xxxyyz_1, g_z_xxxyzz_0, g_z_xxxyzz_1, g_z_xxxzzz_0, g_z_xxxzzz_1, g_z_xxyyyy_0, g_z_xxyyyy_1, g_z_xxyyyz_0, g_z_xxyyyz_1, g_z_xxyyzz_0, g_z_xxyyzz_1, g_z_xxyzzz_0, g_z_xxyzzz_1, g_z_xxzzzz_0, g_z_xxzzzz_1, g_z_xyyyyy_0, g_z_xyyyyy_1, g_z_xyyyyz_0, g_z_xyyyyz_1, g_z_xyyyzz_0, g_z_xyyyzz_1, g_z_xyyzzz_0, g_z_xyyzzz_1, g_z_xyzzzz_0, g_z_xyzzzz_1, g_z_xzzzzz_0, g_z_xzzzzz_1, g_z_yyyyyy_0, g_z_yyyyyy_1, g_z_yyyyyz_0, g_z_yyyyyz_1, g_z_yyyyzz_0, g_z_yyyyzz_1, g_z_yyyzzz_0, g_z_yyyzzz_1, g_z_yyzzzz_0, g_z_yyzzzz_1, g_z_yzzzzz_0, g_z_yzzzzz_1, g_z_zzzzzz_0, g_z_zzzzzz_1, g_zz_xxxxx_1, g_zz_xxxxxx_1, g_zz_xxxxxy_1, g_zz_xxxxxz_1, g_zz_xxxxy_1, g_zz_xxxxyy_1, g_zz_xxxxyz_1, g_zz_xxxxz_1, g_zz_xxxxzz_1, g_zz_xxxyy_1, g_zz_xxxyyy_1, g_zz_xxxyyz_1, g_zz_xxxyz_1, g_zz_xxxyzz_1, g_zz_xxxzz_1, g_zz_xxxzzz_1, g_zz_xxyyy_1, g_zz_xxyyyy_1, g_zz_xxyyyz_1, g_zz_xxyyz_1, g_zz_xxyyzz_1, g_zz_xxyzz_1, g_zz_xxyzzz_1, g_zz_xxzzz_1, g_zz_xxzzzz_1, g_zz_xyyyy_1, g_zz_xyyyyy_1, g_zz_xyyyyz_1, g_zz_xyyyz_1, g_zz_xyyyzz_1, g_zz_xyyzz_1, g_zz_xyyzzz_1, g_zz_xyzzz_1, g_zz_xyzzzz_1, g_zz_xzzzz_1, g_zz_xzzzzz_1, g_zz_yyyyy_1, g_zz_yyyyyy_1, g_zz_yyyyyz_1, g_zz_yyyyz_1, g_zz_yyyyzz_1, g_zz_yyyzz_1, g_zz_yyyzzz_1, g_zz_yyzzz_1, g_zz_yyzzzz_1, g_zz_yzzzz_1, g_zz_yzzzzz_1, g_zz_zzzzz_1, g_zz_zzzzzz_1, g_zzz_xxxxxx_0, g_zzz_xxxxxy_0, g_zzz_xxxxxz_0, g_zzz_xxxxyy_0, g_zzz_xxxxyz_0, g_zzz_xxxxzz_0, g_zzz_xxxyyy_0, g_zzz_xxxyyz_0, g_zzz_xxxyzz_0, g_zzz_xxxzzz_0, g_zzz_xxyyyy_0, g_zzz_xxyyyz_0, g_zzz_xxyyzz_0, g_zzz_xxyzzz_0, g_zzz_xxzzzz_0, g_zzz_xyyyyy_0, g_zzz_xyyyyz_0, g_zzz_xyyyzz_0, g_zzz_xyyzzz_0, g_zzz_xyzzzz_0, g_zzz_xzzzzz_0, g_zzz_yyyyyy_0, g_zzz_yyyyyz_0, g_zzz_yyyyzz_0, g_zzz_yyyzzz_0, g_zzz_yyzzzz_0, g_zzz_yzzzzz_0, g_zzz_zzzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzz_xxxxxx_0[i] = 2.0 * g_z_xxxxxx_0[i] * fbe_0 - 2.0 * g_z_xxxxxx_1[i] * fz_be_0 + g_zz_xxxxxx_1[i] * pa_z[i];

        g_zzz_xxxxxy_0[i] = 2.0 * g_z_xxxxxy_0[i] * fbe_0 - 2.0 * g_z_xxxxxy_1[i] * fz_be_0 + g_zz_xxxxxy_1[i] * pa_z[i];

        g_zzz_xxxxxz_0[i] = 2.0 * g_z_xxxxxz_0[i] * fbe_0 - 2.0 * g_z_xxxxxz_1[i] * fz_be_0 + g_zz_xxxxx_1[i] * fe_0 + g_zz_xxxxxz_1[i] * pa_z[i];

        g_zzz_xxxxyy_0[i] = 2.0 * g_z_xxxxyy_0[i] * fbe_0 - 2.0 * g_z_xxxxyy_1[i] * fz_be_0 + g_zz_xxxxyy_1[i] * pa_z[i];

        g_zzz_xxxxyz_0[i] = 2.0 * g_z_xxxxyz_0[i] * fbe_0 - 2.0 * g_z_xxxxyz_1[i] * fz_be_0 + g_zz_xxxxy_1[i] * fe_0 + g_zz_xxxxyz_1[i] * pa_z[i];

        g_zzz_xxxxzz_0[i] = 2.0 * g_z_xxxxzz_0[i] * fbe_0 - 2.0 * g_z_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zz_xxxxz_1[i] * fe_0 + g_zz_xxxxzz_1[i] * pa_z[i];

        g_zzz_xxxyyy_0[i] = 2.0 * g_z_xxxyyy_0[i] * fbe_0 - 2.0 * g_z_xxxyyy_1[i] * fz_be_0 + g_zz_xxxyyy_1[i] * pa_z[i];

        g_zzz_xxxyyz_0[i] = 2.0 * g_z_xxxyyz_0[i] * fbe_0 - 2.0 * g_z_xxxyyz_1[i] * fz_be_0 + g_zz_xxxyy_1[i] * fe_0 + g_zz_xxxyyz_1[i] * pa_z[i];

        g_zzz_xxxyzz_0[i] = 2.0 * g_z_xxxyzz_0[i] * fbe_0 - 2.0 * g_z_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zz_xxxyz_1[i] * fe_0 + g_zz_xxxyzz_1[i] * pa_z[i];

        g_zzz_xxxzzz_0[i] = 2.0 * g_z_xxxzzz_0[i] * fbe_0 - 2.0 * g_z_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zz_xxxzz_1[i] * fe_0 + g_zz_xxxzzz_1[i] * pa_z[i];

        g_zzz_xxyyyy_0[i] = 2.0 * g_z_xxyyyy_0[i] * fbe_0 - 2.0 * g_z_xxyyyy_1[i] * fz_be_0 + g_zz_xxyyyy_1[i] * pa_z[i];

        g_zzz_xxyyyz_0[i] = 2.0 * g_z_xxyyyz_0[i] * fbe_0 - 2.0 * g_z_xxyyyz_1[i] * fz_be_0 + g_zz_xxyyy_1[i] * fe_0 + g_zz_xxyyyz_1[i] * pa_z[i];

        g_zzz_xxyyzz_0[i] = 2.0 * g_z_xxyyzz_0[i] * fbe_0 - 2.0 * g_z_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zz_xxyyz_1[i] * fe_0 + g_zz_xxyyzz_1[i] * pa_z[i];

        g_zzz_xxyzzz_0[i] = 2.0 * g_z_xxyzzz_0[i] * fbe_0 - 2.0 * g_z_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zz_xxyzz_1[i] * fe_0 + g_zz_xxyzzz_1[i] * pa_z[i];

        g_zzz_xxzzzz_0[i] = 2.0 * g_z_xxzzzz_0[i] * fbe_0 - 2.0 * g_z_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zz_xxzzz_1[i] * fe_0 + g_zz_xxzzzz_1[i] * pa_z[i];

        g_zzz_xyyyyy_0[i] = 2.0 * g_z_xyyyyy_0[i] * fbe_0 - 2.0 * g_z_xyyyyy_1[i] * fz_be_0 + g_zz_xyyyyy_1[i] * pa_z[i];

        g_zzz_xyyyyz_0[i] = 2.0 * g_z_xyyyyz_0[i] * fbe_0 - 2.0 * g_z_xyyyyz_1[i] * fz_be_0 + g_zz_xyyyy_1[i] * fe_0 + g_zz_xyyyyz_1[i] * pa_z[i];

        g_zzz_xyyyzz_0[i] = 2.0 * g_z_xyyyzz_0[i] * fbe_0 - 2.0 * g_z_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_xyyyz_1[i] * fe_0 + g_zz_xyyyzz_1[i] * pa_z[i];

        g_zzz_xyyzzz_0[i] = 2.0 * g_z_xyyzzz_0[i] * fbe_0 - 2.0 * g_z_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_xyyzz_1[i] * fe_0 + g_zz_xyyzzz_1[i] * pa_z[i];

        g_zzz_xyzzzz_0[i] = 2.0 * g_z_xyzzzz_0[i] * fbe_0 - 2.0 * g_z_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_xyzzz_1[i] * fe_0 + g_zz_xyzzzz_1[i] * pa_z[i];

        g_zzz_xzzzzz_0[i] = 2.0 * g_z_xzzzzz_0[i] * fbe_0 - 2.0 * g_z_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_xzzzz_1[i] * fe_0 + g_zz_xzzzzz_1[i] * pa_z[i];

        g_zzz_yyyyyy_0[i] = 2.0 * g_z_yyyyyy_0[i] * fbe_0 - 2.0 * g_z_yyyyyy_1[i] * fz_be_0 + g_zz_yyyyyy_1[i] * pa_z[i];

        g_zzz_yyyyyz_0[i] = 2.0 * g_z_yyyyyz_0[i] * fbe_0 - 2.0 * g_z_yyyyyz_1[i] * fz_be_0 + g_zz_yyyyy_1[i] * fe_0 + g_zz_yyyyyz_1[i] * pa_z[i];

        g_zzz_yyyyzz_0[i] = 2.0 * g_z_yyyyzz_0[i] * fbe_0 - 2.0 * g_z_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zz_yyyyz_1[i] * fe_0 + g_zz_yyyyzz_1[i] * pa_z[i];

        g_zzz_yyyzzz_0[i] = 2.0 * g_z_yyyzzz_0[i] * fbe_0 - 2.0 * g_z_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zz_yyyzz_1[i] * fe_0 + g_zz_yyyzzz_1[i] * pa_z[i];

        g_zzz_yyzzzz_0[i] = 2.0 * g_z_yyzzzz_0[i] * fbe_0 - 2.0 * g_z_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zz_yyzzz_1[i] * fe_0 + g_zz_yyzzzz_1[i] * pa_z[i];

        g_zzz_yzzzzz_0[i] = 2.0 * g_z_yzzzzz_0[i] * fbe_0 - 2.0 * g_z_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zz_yzzzz_1[i] * fe_0 + g_zz_yzzzzz_1[i] * pa_z[i];

        g_zzz_zzzzzz_0[i] = 2.0 * g_z_zzzzzz_0[i] * fbe_0 - 2.0 * g_z_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zz_zzzzz_1[i] * fe_0 + g_zz_zzzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

