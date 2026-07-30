#include "LocalCorePotentialPrimRecDI.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_di(CSimdArray<double>& pbuffer, 
                                  const size_t idx_di,
                                  const size_t idx_si,
                                  const size_t idx_ph,
                                  const size_t idx_pi,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SI

    auto tg_0_xxxxxx = pbuffer.data(idx_si);

    auto tg_0_xxxxxy = pbuffer.data(idx_si + 1);

    auto tg_0_xxxxxz = pbuffer.data(idx_si + 2);

    auto tg_0_xxxxyy = pbuffer.data(idx_si + 3);

    auto tg_0_xxxxyz = pbuffer.data(idx_si + 4);

    auto tg_0_xxxxzz = pbuffer.data(idx_si + 5);

    auto tg_0_xxxyyy = pbuffer.data(idx_si + 6);

    auto tg_0_xxxyyz = pbuffer.data(idx_si + 7);

    auto tg_0_xxxyzz = pbuffer.data(idx_si + 8);

    auto tg_0_xxxzzz = pbuffer.data(idx_si + 9);

    auto tg_0_xxyyyy = pbuffer.data(idx_si + 10);

    auto tg_0_xxyyyz = pbuffer.data(idx_si + 11);

    auto tg_0_xxyyzz = pbuffer.data(idx_si + 12);

    auto tg_0_xxyzzz = pbuffer.data(idx_si + 13);

    auto tg_0_xxzzzz = pbuffer.data(idx_si + 14);

    auto tg_0_xyyyyy = pbuffer.data(idx_si + 15);

    auto tg_0_xyyyyz = pbuffer.data(idx_si + 16);

    auto tg_0_xyyyzz = pbuffer.data(idx_si + 17);

    auto tg_0_xyyzzz = pbuffer.data(idx_si + 18);

    auto tg_0_xyzzzz = pbuffer.data(idx_si + 19);

    auto tg_0_xzzzzz = pbuffer.data(idx_si + 20);

    auto tg_0_yyyyyy = pbuffer.data(idx_si + 21);

    auto tg_0_yyyyyz = pbuffer.data(idx_si + 22);

    auto tg_0_yyyyzz = pbuffer.data(idx_si + 23);

    auto tg_0_yyyzzz = pbuffer.data(idx_si + 24);

    auto tg_0_yyzzzz = pbuffer.data(idx_si + 25);

    auto tg_0_yzzzzz = pbuffer.data(idx_si + 26);

    auto tg_0_zzzzzz = pbuffer.data(idx_si + 27);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx = pbuffer.data(idx_ph);

    auto tg_x_xxxxy = pbuffer.data(idx_ph + 1);

    auto tg_x_xxxxz = pbuffer.data(idx_ph + 2);

    auto tg_x_xxxyy = pbuffer.data(idx_ph + 3);

    auto tg_x_xxxyz = pbuffer.data(idx_ph + 4);

    auto tg_x_xxxzz = pbuffer.data(idx_ph + 5);

    auto tg_x_xxyyy = pbuffer.data(idx_ph + 6);

    auto tg_x_xxyyz = pbuffer.data(idx_ph + 7);

    auto tg_x_xxyzz = pbuffer.data(idx_ph + 8);

    auto tg_x_xxzzz = pbuffer.data(idx_ph + 9);

    auto tg_x_xyyyy = pbuffer.data(idx_ph + 10);

    auto tg_x_xyyyz = pbuffer.data(idx_ph + 11);

    auto tg_x_xyyzz = pbuffer.data(idx_ph + 12);

    auto tg_x_xyzzz = pbuffer.data(idx_ph + 13);

    auto tg_x_xzzzz = pbuffer.data(idx_ph + 14);

    auto tg_x_yyyyy = pbuffer.data(idx_ph + 15);

    auto tg_x_yyyyz = pbuffer.data(idx_ph + 16);

    auto tg_x_yyyzz = pbuffer.data(idx_ph + 17);

    auto tg_x_yyzzz = pbuffer.data(idx_ph + 18);

    auto tg_x_yzzzz = pbuffer.data(idx_ph + 19);

    auto tg_x_zzzzz = pbuffer.data(idx_ph + 20);

    auto tg_y_xxxxx = pbuffer.data(idx_ph + 21);

    auto tg_y_xxxxy = pbuffer.data(idx_ph + 22);

    auto tg_y_xxxxz = pbuffer.data(idx_ph + 23);

    auto tg_y_xxxyy = pbuffer.data(idx_ph + 24);

    auto tg_y_xxxyz = pbuffer.data(idx_ph + 25);

    auto tg_y_xxxzz = pbuffer.data(idx_ph + 26);

    auto tg_y_xxyyy = pbuffer.data(idx_ph + 27);

    auto tg_y_xxyyz = pbuffer.data(idx_ph + 28);

    auto tg_y_xxyzz = pbuffer.data(idx_ph + 29);

    auto tg_y_xxzzz = pbuffer.data(idx_ph + 30);

    auto tg_y_xyyyy = pbuffer.data(idx_ph + 31);

    auto tg_y_xyyyz = pbuffer.data(idx_ph + 32);

    auto tg_y_xyyzz = pbuffer.data(idx_ph + 33);

    auto tg_y_xyzzz = pbuffer.data(idx_ph + 34);

    auto tg_y_xzzzz = pbuffer.data(idx_ph + 35);

    auto tg_y_yyyyy = pbuffer.data(idx_ph + 36);

    auto tg_y_yyyyz = pbuffer.data(idx_ph + 37);

    auto tg_y_yyyzz = pbuffer.data(idx_ph + 38);

    auto tg_y_yyzzz = pbuffer.data(idx_ph + 39);

    auto tg_y_yzzzz = pbuffer.data(idx_ph + 40);

    auto tg_y_zzzzz = pbuffer.data(idx_ph + 41);

    auto tg_z_xxxxx = pbuffer.data(idx_ph + 42);

    auto tg_z_xxxxy = pbuffer.data(idx_ph + 43);

    auto tg_z_xxxxz = pbuffer.data(idx_ph + 44);

    auto tg_z_xxxyy = pbuffer.data(idx_ph + 45);

    auto tg_z_xxxyz = pbuffer.data(idx_ph + 46);

    auto tg_z_xxxzz = pbuffer.data(idx_ph + 47);

    auto tg_z_xxyyy = pbuffer.data(idx_ph + 48);

    auto tg_z_xxyyz = pbuffer.data(idx_ph + 49);

    auto tg_z_xxyzz = pbuffer.data(idx_ph + 50);

    auto tg_z_xxzzz = pbuffer.data(idx_ph + 51);

    auto tg_z_xyyyy = pbuffer.data(idx_ph + 52);

    auto tg_z_xyyyz = pbuffer.data(idx_ph + 53);

    auto tg_z_xyyzz = pbuffer.data(idx_ph + 54);

    auto tg_z_xyzzz = pbuffer.data(idx_ph + 55);

    auto tg_z_xzzzz = pbuffer.data(idx_ph + 56);

    auto tg_z_yyyyy = pbuffer.data(idx_ph + 57);

    auto tg_z_yyyyz = pbuffer.data(idx_ph + 58);

    auto tg_z_yyyzz = pbuffer.data(idx_ph + 59);

    auto tg_z_yyzzz = pbuffer.data(idx_ph + 60);

    auto tg_z_yzzzz = pbuffer.data(idx_ph + 61);

    auto tg_z_zzzzz = pbuffer.data(idx_ph + 62);

    // Set up components of auxiliary buffer : PI

    auto tg_x_xxxxxx = pbuffer.data(idx_pi);

    auto tg_x_xxxxxy = pbuffer.data(idx_pi + 1);

    auto tg_x_xxxxxz = pbuffer.data(idx_pi + 2);

    auto tg_x_xxxxyy = pbuffer.data(idx_pi + 3);

    auto tg_x_xxxxyz = pbuffer.data(idx_pi + 4);

    auto tg_x_xxxxzz = pbuffer.data(idx_pi + 5);

    auto tg_x_xxxyyy = pbuffer.data(idx_pi + 6);

    auto tg_x_xxxyyz = pbuffer.data(idx_pi + 7);

    auto tg_x_xxxyzz = pbuffer.data(idx_pi + 8);

    auto tg_x_xxxzzz = pbuffer.data(idx_pi + 9);

    auto tg_x_xxyyyy = pbuffer.data(idx_pi + 10);

    auto tg_x_xxyyyz = pbuffer.data(idx_pi + 11);

    auto tg_x_xxyyzz = pbuffer.data(idx_pi + 12);

    auto tg_x_xxyzzz = pbuffer.data(idx_pi + 13);

    auto tg_x_xxzzzz = pbuffer.data(idx_pi + 14);

    auto tg_x_xyyyyy = pbuffer.data(idx_pi + 15);

    auto tg_x_xyyyyz = pbuffer.data(idx_pi + 16);

    auto tg_x_xyyyzz = pbuffer.data(idx_pi + 17);

    auto tg_x_xyyzzz = pbuffer.data(idx_pi + 18);

    auto tg_x_xyzzzz = pbuffer.data(idx_pi + 19);

    auto tg_x_xzzzzz = pbuffer.data(idx_pi + 20);

    auto tg_x_yyyyyy = pbuffer.data(idx_pi + 21);

    auto tg_x_yyyyyz = pbuffer.data(idx_pi + 22);

    auto tg_x_yyyyzz = pbuffer.data(idx_pi + 23);

    auto tg_x_yyyzzz = pbuffer.data(idx_pi + 24);

    auto tg_x_yyzzzz = pbuffer.data(idx_pi + 25);

    auto tg_x_yzzzzz = pbuffer.data(idx_pi + 26);

    auto tg_x_zzzzzz = pbuffer.data(idx_pi + 27);

    auto tg_y_xxxxxx = pbuffer.data(idx_pi + 28);

    auto tg_y_xxxxxy = pbuffer.data(idx_pi + 29);

    auto tg_y_xxxxxz = pbuffer.data(idx_pi + 30);

    auto tg_y_xxxxyy = pbuffer.data(idx_pi + 31);

    auto tg_y_xxxxyz = pbuffer.data(idx_pi + 32);

    auto tg_y_xxxxzz = pbuffer.data(idx_pi + 33);

    auto tg_y_xxxyyy = pbuffer.data(idx_pi + 34);

    auto tg_y_xxxyyz = pbuffer.data(idx_pi + 35);

    auto tg_y_xxxyzz = pbuffer.data(idx_pi + 36);

    auto tg_y_xxxzzz = pbuffer.data(idx_pi + 37);

    auto tg_y_xxyyyy = pbuffer.data(idx_pi + 38);

    auto tg_y_xxyyyz = pbuffer.data(idx_pi + 39);

    auto tg_y_xxyyzz = pbuffer.data(idx_pi + 40);

    auto tg_y_xxyzzz = pbuffer.data(idx_pi + 41);

    auto tg_y_xxzzzz = pbuffer.data(idx_pi + 42);

    auto tg_y_xyyyyy = pbuffer.data(idx_pi + 43);

    auto tg_y_xyyyyz = pbuffer.data(idx_pi + 44);

    auto tg_y_xyyyzz = pbuffer.data(idx_pi + 45);

    auto tg_y_xyyzzz = pbuffer.data(idx_pi + 46);

    auto tg_y_xyzzzz = pbuffer.data(idx_pi + 47);

    auto tg_y_xzzzzz = pbuffer.data(idx_pi + 48);

    auto tg_y_yyyyyy = pbuffer.data(idx_pi + 49);

    auto tg_y_yyyyyz = pbuffer.data(idx_pi + 50);

    auto tg_y_yyyyzz = pbuffer.data(idx_pi + 51);

    auto tg_y_yyyzzz = pbuffer.data(idx_pi + 52);

    auto tg_y_yyzzzz = pbuffer.data(idx_pi + 53);

    auto tg_y_yzzzzz = pbuffer.data(idx_pi + 54);

    auto tg_y_zzzzzz = pbuffer.data(idx_pi + 55);

    auto tg_z_xxxxxx = pbuffer.data(idx_pi + 56);

    auto tg_z_xxxxxy = pbuffer.data(idx_pi + 57);

    auto tg_z_xxxxxz = pbuffer.data(idx_pi + 58);

    auto tg_z_xxxxyy = pbuffer.data(idx_pi + 59);

    auto tg_z_xxxxyz = pbuffer.data(idx_pi + 60);

    auto tg_z_xxxxzz = pbuffer.data(idx_pi + 61);

    auto tg_z_xxxyyy = pbuffer.data(idx_pi + 62);

    auto tg_z_xxxyyz = pbuffer.data(idx_pi + 63);

    auto tg_z_xxxyzz = pbuffer.data(idx_pi + 64);

    auto tg_z_xxxzzz = pbuffer.data(idx_pi + 65);

    auto tg_z_xxyyyy = pbuffer.data(idx_pi + 66);

    auto tg_z_xxyyyz = pbuffer.data(idx_pi + 67);

    auto tg_z_xxyyzz = pbuffer.data(idx_pi + 68);

    auto tg_z_xxyzzz = pbuffer.data(idx_pi + 69);

    auto tg_z_xxzzzz = pbuffer.data(idx_pi + 70);

    auto tg_z_xyyyyy = pbuffer.data(idx_pi + 71);

    auto tg_z_xyyyyz = pbuffer.data(idx_pi + 72);

    auto tg_z_xyyyzz = pbuffer.data(idx_pi + 73);

    auto tg_z_xyyzzz = pbuffer.data(idx_pi + 74);

    auto tg_z_xyzzzz = pbuffer.data(idx_pi + 75);

    auto tg_z_xzzzzz = pbuffer.data(idx_pi + 76);

    auto tg_z_yyyyyy = pbuffer.data(idx_pi + 77);

    auto tg_z_yyyyyz = pbuffer.data(idx_pi + 78);

    auto tg_z_yyyyzz = pbuffer.data(idx_pi + 79);

    auto tg_z_yyyzzz = pbuffer.data(idx_pi + 80);

    auto tg_z_yyzzzz = pbuffer.data(idx_pi + 81);

    auto tg_z_yzzzzz = pbuffer.data(idx_pi + 82);

    auto tg_z_zzzzzz = pbuffer.data(idx_pi + 83);

    // Set up components of targeted buffer : DI

    auto tg_xx_xxxxxx = pbuffer.data(idx_di);

    auto tg_xx_xxxxxy = pbuffer.data(idx_di + 1);

    auto tg_xx_xxxxxz = pbuffer.data(idx_di + 2);

    auto tg_xx_xxxxyy = pbuffer.data(idx_di + 3);

    auto tg_xx_xxxxyz = pbuffer.data(idx_di + 4);

    auto tg_xx_xxxxzz = pbuffer.data(idx_di + 5);

    auto tg_xx_xxxyyy = pbuffer.data(idx_di + 6);

    auto tg_xx_xxxyyz = pbuffer.data(idx_di + 7);

    auto tg_xx_xxxyzz = pbuffer.data(idx_di + 8);

    auto tg_xx_xxxzzz = pbuffer.data(idx_di + 9);

    auto tg_xx_xxyyyy = pbuffer.data(idx_di + 10);

    auto tg_xx_xxyyyz = pbuffer.data(idx_di + 11);

    auto tg_xx_xxyyzz = pbuffer.data(idx_di + 12);

    auto tg_xx_xxyzzz = pbuffer.data(idx_di + 13);

    auto tg_xx_xxzzzz = pbuffer.data(idx_di + 14);

    auto tg_xx_xyyyyy = pbuffer.data(idx_di + 15);

    auto tg_xx_xyyyyz = pbuffer.data(idx_di + 16);

    auto tg_xx_xyyyzz = pbuffer.data(idx_di + 17);

    auto tg_xx_xyyzzz = pbuffer.data(idx_di + 18);

    auto tg_xx_xyzzzz = pbuffer.data(idx_di + 19);

    auto tg_xx_xzzzzz = pbuffer.data(idx_di + 20);

    auto tg_xx_yyyyyy = pbuffer.data(idx_di + 21);

    auto tg_xx_yyyyyz = pbuffer.data(idx_di + 22);

    auto tg_xx_yyyyzz = pbuffer.data(idx_di + 23);

    auto tg_xx_yyyzzz = pbuffer.data(idx_di + 24);

    auto tg_xx_yyzzzz = pbuffer.data(idx_di + 25);

    auto tg_xx_yzzzzz = pbuffer.data(idx_di + 26);

    auto tg_xx_zzzzzz = pbuffer.data(idx_di + 27);

    auto tg_xy_xxxxxx = pbuffer.data(idx_di + 28);

    auto tg_xy_xxxxxy = pbuffer.data(idx_di + 29);

    auto tg_xy_xxxxxz = pbuffer.data(idx_di + 30);

    auto tg_xy_xxxxyy = pbuffer.data(idx_di + 31);

    auto tg_xy_xxxxyz = pbuffer.data(idx_di + 32);

    auto tg_xy_xxxxzz = pbuffer.data(idx_di + 33);

    auto tg_xy_xxxyyy = pbuffer.data(idx_di + 34);

    auto tg_xy_xxxyyz = pbuffer.data(idx_di + 35);

    auto tg_xy_xxxyzz = pbuffer.data(idx_di + 36);

    auto tg_xy_xxxzzz = pbuffer.data(idx_di + 37);

    auto tg_xy_xxyyyy = pbuffer.data(idx_di + 38);

    auto tg_xy_xxyyyz = pbuffer.data(idx_di + 39);

    auto tg_xy_xxyyzz = pbuffer.data(idx_di + 40);

    auto tg_xy_xxyzzz = pbuffer.data(idx_di + 41);

    auto tg_xy_xxzzzz = pbuffer.data(idx_di + 42);

    auto tg_xy_xyyyyy = pbuffer.data(idx_di + 43);

    auto tg_xy_xyyyyz = pbuffer.data(idx_di + 44);

    auto tg_xy_xyyyzz = pbuffer.data(idx_di + 45);

    auto tg_xy_xyyzzz = pbuffer.data(idx_di + 46);

    auto tg_xy_xyzzzz = pbuffer.data(idx_di + 47);

    auto tg_xy_xzzzzz = pbuffer.data(idx_di + 48);

    auto tg_xy_yyyyyy = pbuffer.data(idx_di + 49);

    auto tg_xy_yyyyyz = pbuffer.data(idx_di + 50);

    auto tg_xy_yyyyzz = pbuffer.data(idx_di + 51);

    auto tg_xy_yyyzzz = pbuffer.data(idx_di + 52);

    auto tg_xy_yyzzzz = pbuffer.data(idx_di + 53);

    auto tg_xy_yzzzzz = pbuffer.data(idx_di + 54);

    auto tg_xy_zzzzzz = pbuffer.data(idx_di + 55);

    auto tg_xz_xxxxxx = pbuffer.data(idx_di + 56);

    auto tg_xz_xxxxxy = pbuffer.data(idx_di + 57);

    auto tg_xz_xxxxxz = pbuffer.data(idx_di + 58);

    auto tg_xz_xxxxyy = pbuffer.data(idx_di + 59);

    auto tg_xz_xxxxyz = pbuffer.data(idx_di + 60);

    auto tg_xz_xxxxzz = pbuffer.data(idx_di + 61);

    auto tg_xz_xxxyyy = pbuffer.data(idx_di + 62);

    auto tg_xz_xxxyyz = pbuffer.data(idx_di + 63);

    auto tg_xz_xxxyzz = pbuffer.data(idx_di + 64);

    auto tg_xz_xxxzzz = pbuffer.data(idx_di + 65);

    auto tg_xz_xxyyyy = pbuffer.data(idx_di + 66);

    auto tg_xz_xxyyyz = pbuffer.data(idx_di + 67);

    auto tg_xz_xxyyzz = pbuffer.data(idx_di + 68);

    auto tg_xz_xxyzzz = pbuffer.data(idx_di + 69);

    auto tg_xz_xxzzzz = pbuffer.data(idx_di + 70);

    auto tg_xz_xyyyyy = pbuffer.data(idx_di + 71);

    auto tg_xz_xyyyyz = pbuffer.data(idx_di + 72);

    auto tg_xz_xyyyzz = pbuffer.data(idx_di + 73);

    auto tg_xz_xyyzzz = pbuffer.data(idx_di + 74);

    auto tg_xz_xyzzzz = pbuffer.data(idx_di + 75);

    auto tg_xz_xzzzzz = pbuffer.data(idx_di + 76);

    auto tg_xz_yyyyyy = pbuffer.data(idx_di + 77);

    auto tg_xz_yyyyyz = pbuffer.data(idx_di + 78);

    auto tg_xz_yyyyzz = pbuffer.data(idx_di + 79);

    auto tg_xz_yyyzzz = pbuffer.data(idx_di + 80);

    auto tg_xz_yyzzzz = pbuffer.data(idx_di + 81);

    auto tg_xz_yzzzzz = pbuffer.data(idx_di + 82);

    auto tg_xz_zzzzzz = pbuffer.data(idx_di + 83);

    auto tg_yy_xxxxxx = pbuffer.data(idx_di + 84);

    auto tg_yy_xxxxxy = pbuffer.data(idx_di + 85);

    auto tg_yy_xxxxxz = pbuffer.data(idx_di + 86);

    auto tg_yy_xxxxyy = pbuffer.data(idx_di + 87);

    auto tg_yy_xxxxyz = pbuffer.data(idx_di + 88);

    auto tg_yy_xxxxzz = pbuffer.data(idx_di + 89);

    auto tg_yy_xxxyyy = pbuffer.data(idx_di + 90);

    auto tg_yy_xxxyyz = pbuffer.data(idx_di + 91);

    auto tg_yy_xxxyzz = pbuffer.data(idx_di + 92);

    auto tg_yy_xxxzzz = pbuffer.data(idx_di + 93);

    auto tg_yy_xxyyyy = pbuffer.data(idx_di + 94);

    auto tg_yy_xxyyyz = pbuffer.data(idx_di + 95);

    auto tg_yy_xxyyzz = pbuffer.data(idx_di + 96);

    auto tg_yy_xxyzzz = pbuffer.data(idx_di + 97);

    auto tg_yy_xxzzzz = pbuffer.data(idx_di + 98);

    auto tg_yy_xyyyyy = pbuffer.data(idx_di + 99);

    auto tg_yy_xyyyyz = pbuffer.data(idx_di + 100);

    auto tg_yy_xyyyzz = pbuffer.data(idx_di + 101);

    auto tg_yy_xyyzzz = pbuffer.data(idx_di + 102);

    auto tg_yy_xyzzzz = pbuffer.data(idx_di + 103);

    auto tg_yy_xzzzzz = pbuffer.data(idx_di + 104);

    auto tg_yy_yyyyyy = pbuffer.data(idx_di + 105);

    auto tg_yy_yyyyyz = pbuffer.data(idx_di + 106);

    auto tg_yy_yyyyzz = pbuffer.data(idx_di + 107);

    auto tg_yy_yyyzzz = pbuffer.data(idx_di + 108);

    auto tg_yy_yyzzzz = pbuffer.data(idx_di + 109);

    auto tg_yy_yzzzzz = pbuffer.data(idx_di + 110);

    auto tg_yy_zzzzzz = pbuffer.data(idx_di + 111);

    auto tg_yz_xxxxxx = pbuffer.data(idx_di + 112);

    auto tg_yz_xxxxxy = pbuffer.data(idx_di + 113);

    auto tg_yz_xxxxxz = pbuffer.data(idx_di + 114);

    auto tg_yz_xxxxyy = pbuffer.data(idx_di + 115);

    auto tg_yz_xxxxyz = pbuffer.data(idx_di + 116);

    auto tg_yz_xxxxzz = pbuffer.data(idx_di + 117);

    auto tg_yz_xxxyyy = pbuffer.data(idx_di + 118);

    auto tg_yz_xxxyyz = pbuffer.data(idx_di + 119);

    auto tg_yz_xxxyzz = pbuffer.data(idx_di + 120);

    auto tg_yz_xxxzzz = pbuffer.data(idx_di + 121);

    auto tg_yz_xxyyyy = pbuffer.data(idx_di + 122);

    auto tg_yz_xxyyyz = pbuffer.data(idx_di + 123);

    auto tg_yz_xxyyzz = pbuffer.data(idx_di + 124);

    auto tg_yz_xxyzzz = pbuffer.data(idx_di + 125);

    auto tg_yz_xxzzzz = pbuffer.data(idx_di + 126);

    auto tg_yz_xyyyyy = pbuffer.data(idx_di + 127);

    auto tg_yz_xyyyyz = pbuffer.data(idx_di + 128);

    auto tg_yz_xyyyzz = pbuffer.data(idx_di + 129);

    auto tg_yz_xyyzzz = pbuffer.data(idx_di + 130);

    auto tg_yz_xyzzzz = pbuffer.data(idx_di + 131);

    auto tg_yz_xzzzzz = pbuffer.data(idx_di + 132);

    auto tg_yz_yyyyyy = pbuffer.data(idx_di + 133);

    auto tg_yz_yyyyyz = pbuffer.data(idx_di + 134);

    auto tg_yz_yyyyzz = pbuffer.data(idx_di + 135);

    auto tg_yz_yyyzzz = pbuffer.data(idx_di + 136);

    auto tg_yz_yyzzzz = pbuffer.data(idx_di + 137);

    auto tg_yz_yzzzzz = pbuffer.data(idx_di + 138);

    auto tg_yz_zzzzzz = pbuffer.data(idx_di + 139);

    auto tg_zz_xxxxxx = pbuffer.data(idx_di + 140);

    auto tg_zz_xxxxxy = pbuffer.data(idx_di + 141);

    auto tg_zz_xxxxxz = pbuffer.data(idx_di + 142);

    auto tg_zz_xxxxyy = pbuffer.data(idx_di + 143);

    auto tg_zz_xxxxyz = pbuffer.data(idx_di + 144);

    auto tg_zz_xxxxzz = pbuffer.data(idx_di + 145);

    auto tg_zz_xxxyyy = pbuffer.data(idx_di + 146);

    auto tg_zz_xxxyyz = pbuffer.data(idx_di + 147);

    auto tg_zz_xxxyzz = pbuffer.data(idx_di + 148);

    auto tg_zz_xxxzzz = pbuffer.data(idx_di + 149);

    auto tg_zz_xxyyyy = pbuffer.data(idx_di + 150);

    auto tg_zz_xxyyyz = pbuffer.data(idx_di + 151);

    auto tg_zz_xxyyzz = pbuffer.data(idx_di + 152);

    auto tg_zz_xxyzzz = pbuffer.data(idx_di + 153);

    auto tg_zz_xxzzzz = pbuffer.data(idx_di + 154);

    auto tg_zz_xyyyyy = pbuffer.data(idx_di + 155);

    auto tg_zz_xyyyyz = pbuffer.data(idx_di + 156);

    auto tg_zz_xyyyzz = pbuffer.data(idx_di + 157);

    auto tg_zz_xyyzzz = pbuffer.data(idx_di + 158);

    auto tg_zz_xyzzzz = pbuffer.data(idx_di + 159);

    auto tg_zz_xzzzzz = pbuffer.data(idx_di + 160);

    auto tg_zz_yyyyyy = pbuffer.data(idx_di + 161);

    auto tg_zz_yyyyyz = pbuffer.data(idx_di + 162);

    auto tg_zz_yyyyzz = pbuffer.data(idx_di + 163);

    auto tg_zz_yyyzzz = pbuffer.data(idx_di + 164);

    auto tg_zz_yyzzzz = pbuffer.data(idx_di + 165);

    auto tg_zz_yzzzzz = pbuffer.data(idx_di + 166);

    auto tg_zz_zzzzzz = pbuffer.data(idx_di + 167);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxxxxx, tg_0_xxxxxy, tg_0_xxxxxz, tg_0_xxxxyy, tg_0_xxxxyz, tg_0_xxxxzz, tg_0_xxxyyy, tg_0_xxxyyz, tg_0_xxxyzz, tg_0_xxxzzz, tg_0_xxyyyy, tg_0_xxyyyz, tg_0_xxyyzz, tg_0_xxyzzz, tg_0_xxzzzz, tg_0_xyyyyy, tg_0_xyyyyz, tg_0_xyyyzz, tg_0_xyyzzz, tg_0_xyzzzz, tg_0_xzzzzz, tg_0_yyyyyy, tg_0_yyyyyz, tg_0_yyyyzz, tg_0_yyyzzz, tg_0_yyzzzz, tg_0_yzzzzz, tg_0_zzzzzz, tg_x_xxxxx, tg_x_xxxxxx, tg_x_xxxxxy, tg_x_xxxxxz, tg_x_xxxxy, tg_x_xxxxyy, tg_x_xxxxyz, tg_x_xxxxz, tg_x_xxxxzz, tg_x_xxxyy, tg_x_xxxyyy, tg_x_xxxyyz, tg_x_xxxyz, tg_x_xxxyzz, tg_x_xxxzz, tg_x_xxxzzz, tg_x_xxyyy, tg_x_xxyyyy, tg_x_xxyyyz, tg_x_xxyyz, tg_x_xxyyzz, tg_x_xxyzz, tg_x_xxyzzz, tg_x_xxzzz, tg_x_xxzzzz, tg_x_xyyyy, tg_x_xyyyyy, tg_x_xyyyyz, tg_x_xyyyz, tg_x_xyyyzz, tg_x_xyyzz, tg_x_xyyzzz, tg_x_xyzzz, tg_x_xyzzzz, tg_x_xzzzz, tg_x_xzzzzz, tg_x_yyyyy, tg_x_yyyyyy, tg_x_yyyyyz, tg_x_yyyyz, tg_x_yyyyzz, tg_x_yyyzz, tg_x_yyyzzz, tg_x_yyzzz, tg_x_yyzzzz, tg_x_yzzzz, tg_x_yzzzzz, tg_x_zzzzz, tg_x_zzzzzz, tg_xx_xxxxxx, tg_xx_xxxxxy, tg_xx_xxxxxz, tg_xx_xxxxyy, tg_xx_xxxxyz, tg_xx_xxxxzz, tg_xx_xxxyyy, tg_xx_xxxyyz, tg_xx_xxxyzz, tg_xx_xxxzzz, tg_xx_xxyyyy, tg_xx_xxyyyz, tg_xx_xxyyzz, tg_xx_xxyzzz, tg_xx_xxzzzz, tg_xx_xyyyyy, tg_xx_xyyyyz, tg_xx_xyyyzz, tg_xx_xyyzzz, tg_xx_xyzzzz, tg_xx_xzzzzz, tg_xx_yyyyyy, tg_xx_yyyyyz, tg_xx_yyyyzz, tg_xx_yyyzzz, tg_xx_yyzzzz, tg_xx_yzzzzz, tg_xx_zzzzzz, tg_xy_xxxxxx, tg_xy_xxxxxy, tg_xy_xxxxxz, tg_xy_xxxxyy, tg_xy_xxxxyz, tg_xy_xxxxzz, tg_xy_xxxyyy, tg_xy_xxxyyz, tg_xy_xxxyzz, tg_xy_xxxzzz, tg_xy_xxyyyy, tg_xy_xxyyyz, tg_xy_xxyyzz, tg_xy_xxyzzz, tg_xy_xxzzzz, tg_xy_xyyyyy, tg_xy_xyyyyz, tg_xy_xyyyzz, tg_xy_xyyzzz, tg_xy_xyzzzz, tg_xy_xzzzzz, tg_xy_yyyyyy, tg_xy_yyyyyz, tg_xy_yyyyzz, tg_xy_yyyzzz, tg_xy_yyzzzz, tg_xy_yzzzzz, tg_xy_zzzzzz, tg_xz_xxxxxx, tg_xz_xxxxxy, tg_xz_xxxxxz, tg_xz_xxxxyy, tg_xz_xxxxyz, tg_xz_xxxxzz, tg_xz_xxxyyy, tg_xz_xxxyyz, tg_xz_xxxyzz, tg_xz_xxxzzz, tg_xz_xxyyyy, tg_xz_xxyyyz, tg_xz_xxyyzz, tg_xz_xxyzzz, tg_xz_xxzzzz, tg_xz_xyyyyy, tg_xz_xyyyyz, tg_xz_xyyyzz, tg_xz_xyyzzz, tg_xz_xyzzzz, tg_xz_xzzzzz, tg_xz_yyyyyy, tg_xz_yyyyyz, tg_xz_yyyyzz, tg_xz_yyyzzz, tg_xz_yyzzzz, tg_xz_yzzzzz, tg_xz_zzzzzz, tg_y_xxxxx, tg_y_xxxxxx, tg_y_xxxxxy, tg_y_xxxxxz, tg_y_xxxxy, tg_y_xxxxyy, tg_y_xxxxyz, tg_y_xxxxz, tg_y_xxxxzz, tg_y_xxxyy, tg_y_xxxyyy, tg_y_xxxyyz, tg_y_xxxyz, tg_y_xxxyzz, tg_y_xxxzz, tg_y_xxxzzz, tg_y_xxyyy, tg_y_xxyyyy, tg_y_xxyyyz, tg_y_xxyyz, tg_y_xxyyzz, tg_y_xxyzz, tg_y_xxyzzz, tg_y_xxzzz, tg_y_xxzzzz, tg_y_xyyyy, tg_y_xyyyyy, tg_y_xyyyyz, tg_y_xyyyz, tg_y_xyyyzz, tg_y_xyyzz, tg_y_xyyzzz, tg_y_xyzzz, tg_y_xyzzzz, tg_y_xzzzz, tg_y_xzzzzz, tg_y_yyyyy, tg_y_yyyyyy, tg_y_yyyyyz, tg_y_yyyyz, tg_y_yyyyzz, tg_y_yyyzz, tg_y_yyyzzz, tg_y_yyzzz, tg_y_yyzzzz, tg_y_yzzzz, tg_y_yzzzzz, tg_y_zzzzz, tg_y_zzzzzz, tg_yy_xxxxxx, tg_yy_xxxxxy, tg_yy_xxxxxz, tg_yy_xxxxyy, tg_yy_xxxxyz, tg_yy_xxxxzz, tg_yy_xxxyyy, tg_yy_xxxyyz, tg_yy_xxxyzz, tg_yy_xxxzzz, tg_yy_xxyyyy, tg_yy_xxyyyz, tg_yy_xxyyzz, tg_yy_xxyzzz, tg_yy_xxzzzz, tg_yy_xyyyyy, tg_yy_xyyyyz, tg_yy_xyyyzz, tg_yy_xyyzzz, tg_yy_xyzzzz, tg_yy_xzzzzz, tg_yy_yyyyyy, tg_yy_yyyyyz, tg_yy_yyyyzz, tg_yy_yyyzzz, tg_yy_yyzzzz, tg_yy_yzzzzz, tg_yy_zzzzzz, tg_yz_xxxxxx, tg_yz_xxxxxy, tg_yz_xxxxxz, tg_yz_xxxxyy, tg_yz_xxxxyz, tg_yz_xxxxzz, tg_yz_xxxyyy, tg_yz_xxxyyz, tg_yz_xxxyzz, tg_yz_xxxzzz, tg_yz_xxyyyy, tg_yz_xxyyyz, tg_yz_xxyyzz, tg_yz_xxyzzz, tg_yz_xxzzzz, tg_yz_xyyyyy, tg_yz_xyyyyz, tg_yz_xyyyzz, tg_yz_xyyzzz, tg_yz_xyzzzz, tg_yz_xzzzzz, tg_yz_yyyyyy, tg_yz_yyyyyz, tg_yz_yyyyzz, tg_yz_yyyzzz, tg_yz_yyzzzz, tg_yz_yzzzzz, tg_yz_zzzzzz, tg_z_xxxxx, tg_z_xxxxxx, tg_z_xxxxxy, tg_z_xxxxxz, tg_z_xxxxy, tg_z_xxxxyy, tg_z_xxxxyz, tg_z_xxxxz, tg_z_xxxxzz, tg_z_xxxyy, tg_z_xxxyyy, tg_z_xxxyyz, tg_z_xxxyz, tg_z_xxxyzz, tg_z_xxxzz, tg_z_xxxzzz, tg_z_xxyyy, tg_z_xxyyyy, tg_z_xxyyyz, tg_z_xxyyz, tg_z_xxyyzz, tg_z_xxyzz, tg_z_xxyzzz, tg_z_xxzzz, tg_z_xxzzzz, tg_z_xyyyy, tg_z_xyyyyy, tg_z_xyyyyz, tg_z_xyyyz, tg_z_xyyyzz, tg_z_xyyzz, tg_z_xyyzzz, tg_z_xyzzz, tg_z_xyzzzz, tg_z_xzzzz, tg_z_xzzzzz, tg_z_yyyyy, tg_z_yyyyyy, tg_z_yyyyyz, tg_z_yyyyz, tg_z_yyyyzz, tg_z_yyyzz, tg_z_yyyzzz, tg_z_yyzzz, tg_z_yyzzzz, tg_z_yzzzz, tg_z_yzzzzz, tg_z_zzzzz, tg_z_zzzzzz, tg_zz_xxxxxx, tg_zz_xxxxxy, tg_zz_xxxxxz, tg_zz_xxxxyy, tg_zz_xxxxyz, tg_zz_xxxxzz, tg_zz_xxxyyy, tg_zz_xxxyyz, tg_zz_xxxyzz, tg_zz_xxxzzz, tg_zz_xxyyyy, tg_zz_xxyyyz, tg_zz_xxyyzz, tg_zz_xxyzzz, tg_zz_xxzzzz, tg_zz_xyyyyy, tg_zz_xyyyyz, tg_zz_xyyyzz, tg_zz_xyyzzz, tg_zz_xyzzzz, tg_zz_xzzzzz, tg_zz_yyyyyy, tg_zz_yyyyyz, tg_zz_yyyyzz, tg_zz_yyyzzz, tg_zz_yyzzzz, tg_zz_yzzzzz, tg_zz_zzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xx_xxxxxx[i] = tg_0_xxxxxx[i] * fxi[i] + 6.0 * tg_x_xxxxx[i] * fxi[i] + tg_x_xxxxxx[i] * ra_x[i];

        tg_xx_xxxxxy[i] = tg_0_xxxxxy[i] * fxi[i] + 5.0 * tg_x_xxxxy[i] * fxi[i] + tg_x_xxxxxy[i] * ra_x[i];

        tg_xx_xxxxxz[i] = tg_0_xxxxxz[i] * fxi[i] + 5.0 * tg_x_xxxxz[i] * fxi[i] + tg_x_xxxxxz[i] * ra_x[i];

        tg_xx_xxxxyy[i] = tg_0_xxxxyy[i] * fxi[i] + 4.0 * tg_x_xxxyy[i] * fxi[i] + tg_x_xxxxyy[i] * ra_x[i];

        tg_xx_xxxxyz[i] = tg_0_xxxxyz[i] * fxi[i] + 4.0 * tg_x_xxxyz[i] * fxi[i] + tg_x_xxxxyz[i] * ra_x[i];

        tg_xx_xxxxzz[i] = tg_0_xxxxzz[i] * fxi[i] + 4.0 * tg_x_xxxzz[i] * fxi[i] + tg_x_xxxxzz[i] * ra_x[i];

        tg_xx_xxxyyy[i] = tg_0_xxxyyy[i] * fxi[i] + 3.0 * tg_x_xxyyy[i] * fxi[i] + tg_x_xxxyyy[i] * ra_x[i];

        tg_xx_xxxyyz[i] = tg_0_xxxyyz[i] * fxi[i] + 3.0 * tg_x_xxyyz[i] * fxi[i] + tg_x_xxxyyz[i] * ra_x[i];

        tg_xx_xxxyzz[i] = tg_0_xxxyzz[i] * fxi[i] + 3.0 * tg_x_xxyzz[i] * fxi[i] + tg_x_xxxyzz[i] * ra_x[i];

        tg_xx_xxxzzz[i] = tg_0_xxxzzz[i] * fxi[i] + 3.0 * tg_x_xxzzz[i] * fxi[i] + tg_x_xxxzzz[i] * ra_x[i];

        tg_xx_xxyyyy[i] = tg_0_xxyyyy[i] * fxi[i] + 2.0 * tg_x_xyyyy[i] * fxi[i] + tg_x_xxyyyy[i] * ra_x[i];

        tg_xx_xxyyyz[i] = tg_0_xxyyyz[i] * fxi[i] + 2.0 * tg_x_xyyyz[i] * fxi[i] + tg_x_xxyyyz[i] * ra_x[i];

        tg_xx_xxyyzz[i] = tg_0_xxyyzz[i] * fxi[i] + 2.0 * tg_x_xyyzz[i] * fxi[i] + tg_x_xxyyzz[i] * ra_x[i];

        tg_xx_xxyzzz[i] = tg_0_xxyzzz[i] * fxi[i] + 2.0 * tg_x_xyzzz[i] * fxi[i] + tg_x_xxyzzz[i] * ra_x[i];

        tg_xx_xxzzzz[i] = tg_0_xxzzzz[i] * fxi[i] + 2.0 * tg_x_xzzzz[i] * fxi[i] + tg_x_xxzzzz[i] * ra_x[i];

        tg_xx_xyyyyy[i] = tg_0_xyyyyy[i] * fxi[i] + tg_x_yyyyy[i] * fxi[i] + tg_x_xyyyyy[i] * ra_x[i];

        tg_xx_xyyyyz[i] = tg_0_xyyyyz[i] * fxi[i] + tg_x_yyyyz[i] * fxi[i] + tg_x_xyyyyz[i] * ra_x[i];

        tg_xx_xyyyzz[i] = tg_0_xyyyzz[i] * fxi[i] + tg_x_yyyzz[i] * fxi[i] + tg_x_xyyyzz[i] * ra_x[i];

        tg_xx_xyyzzz[i] = tg_0_xyyzzz[i] * fxi[i] + tg_x_yyzzz[i] * fxi[i] + tg_x_xyyzzz[i] * ra_x[i];

        tg_xx_xyzzzz[i] = tg_0_xyzzzz[i] * fxi[i] + tg_x_yzzzz[i] * fxi[i] + tg_x_xyzzzz[i] * ra_x[i];

        tg_xx_xzzzzz[i] = tg_0_xzzzzz[i] * fxi[i] + tg_x_zzzzz[i] * fxi[i] + tg_x_xzzzzz[i] * ra_x[i];

        tg_xx_yyyyyy[i] = tg_0_yyyyyy[i] * fxi[i] + tg_x_yyyyyy[i] * ra_x[i];

        tg_xx_yyyyyz[i] = tg_0_yyyyyz[i] * fxi[i] + tg_x_yyyyyz[i] * ra_x[i];

        tg_xx_yyyyzz[i] = tg_0_yyyyzz[i] * fxi[i] + tg_x_yyyyzz[i] * ra_x[i];

        tg_xx_yyyzzz[i] = tg_0_yyyzzz[i] * fxi[i] + tg_x_yyyzzz[i] * ra_x[i];

        tg_xx_yyzzzz[i] = tg_0_yyzzzz[i] * fxi[i] + tg_x_yyzzzz[i] * ra_x[i];

        tg_xx_yzzzzz[i] = tg_0_yzzzzz[i] * fxi[i] + tg_x_yzzzzz[i] * ra_x[i];

        tg_xx_zzzzzz[i] = tg_0_zzzzzz[i] * fxi[i] + tg_x_zzzzzz[i] * ra_x[i];

        tg_xy_xxxxxx[i] = tg_x_xxxxxx[i] * ra_y[i];

        tg_xy_xxxxxy[i] = 5.0 * tg_y_xxxxy[i] * fxi[i] + tg_y_xxxxxy[i] * ra_x[i];

        tg_xy_xxxxxz[i] = tg_x_xxxxxz[i] * ra_y[i];

        tg_xy_xxxxyy[i] = 4.0 * tg_y_xxxyy[i] * fxi[i] + tg_y_xxxxyy[i] * ra_x[i];

        tg_xy_xxxxyz[i] = 4.0 * tg_y_xxxyz[i] * fxi[i] + tg_y_xxxxyz[i] * ra_x[i];

        tg_xy_xxxxzz[i] = tg_x_xxxxzz[i] * ra_y[i];

        tg_xy_xxxyyy[i] = 3.0 * tg_y_xxyyy[i] * fxi[i] + tg_y_xxxyyy[i] * ra_x[i];

        tg_xy_xxxyyz[i] = 3.0 * tg_y_xxyyz[i] * fxi[i] + tg_y_xxxyyz[i] * ra_x[i];

        tg_xy_xxxyzz[i] = 3.0 * tg_y_xxyzz[i] * fxi[i] + tg_y_xxxyzz[i] * ra_x[i];

        tg_xy_xxxzzz[i] = tg_x_xxxzzz[i] * ra_y[i];

        tg_xy_xxyyyy[i] = 2.0 * tg_y_xyyyy[i] * fxi[i] + tg_y_xxyyyy[i] * ra_x[i];

        tg_xy_xxyyyz[i] = 2.0 * tg_y_xyyyz[i] * fxi[i] + tg_y_xxyyyz[i] * ra_x[i];

        tg_xy_xxyyzz[i] = 2.0 * tg_y_xyyzz[i] * fxi[i] + tg_y_xxyyzz[i] * ra_x[i];

        tg_xy_xxyzzz[i] = 2.0 * tg_y_xyzzz[i] * fxi[i] + tg_y_xxyzzz[i] * ra_x[i];

        tg_xy_xxzzzz[i] = tg_x_xxzzzz[i] * ra_y[i];

        tg_xy_xyyyyy[i] = tg_y_yyyyy[i] * fxi[i] + tg_y_xyyyyy[i] * ra_x[i];

        tg_xy_xyyyyz[i] = tg_y_yyyyz[i] * fxi[i] + tg_y_xyyyyz[i] * ra_x[i];

        tg_xy_xyyyzz[i] = tg_y_yyyzz[i] * fxi[i] + tg_y_xyyyzz[i] * ra_x[i];

        tg_xy_xyyzzz[i] = tg_y_yyzzz[i] * fxi[i] + tg_y_xyyzzz[i] * ra_x[i];

        tg_xy_xyzzzz[i] = tg_y_yzzzz[i] * fxi[i] + tg_y_xyzzzz[i] * ra_x[i];

        tg_xy_xzzzzz[i] = tg_x_xzzzzz[i] * ra_y[i];

        tg_xy_yyyyyy[i] = tg_y_yyyyyy[i] * ra_x[i];

        tg_xy_yyyyyz[i] = tg_y_yyyyyz[i] * ra_x[i];

        tg_xy_yyyyzz[i] = tg_y_yyyyzz[i] * ra_x[i];

        tg_xy_yyyzzz[i] = tg_y_yyyzzz[i] * ra_x[i];

        tg_xy_yyzzzz[i] = tg_y_yyzzzz[i] * ra_x[i];

        tg_xy_yzzzzz[i] = tg_y_yzzzzz[i] * ra_x[i];

        tg_xy_zzzzzz[i] = tg_y_zzzzzz[i] * ra_x[i];

        tg_xz_xxxxxx[i] = tg_x_xxxxxx[i] * ra_z[i];

        tg_xz_xxxxxy[i] = tg_x_xxxxxy[i] * ra_z[i];

        tg_xz_xxxxxz[i] = 5.0 * tg_z_xxxxz[i] * fxi[i] + tg_z_xxxxxz[i] * ra_x[i];

        tg_xz_xxxxyy[i] = tg_x_xxxxyy[i] * ra_z[i];

        tg_xz_xxxxyz[i] = 4.0 * tg_z_xxxyz[i] * fxi[i] + tg_z_xxxxyz[i] * ra_x[i];

        tg_xz_xxxxzz[i] = 4.0 * tg_z_xxxzz[i] * fxi[i] + tg_z_xxxxzz[i] * ra_x[i];

        tg_xz_xxxyyy[i] = tg_x_xxxyyy[i] * ra_z[i];

        tg_xz_xxxyyz[i] = 3.0 * tg_z_xxyyz[i] * fxi[i] + tg_z_xxxyyz[i] * ra_x[i];

        tg_xz_xxxyzz[i] = 3.0 * tg_z_xxyzz[i] * fxi[i] + tg_z_xxxyzz[i] * ra_x[i];

        tg_xz_xxxzzz[i] = 3.0 * tg_z_xxzzz[i] * fxi[i] + tg_z_xxxzzz[i] * ra_x[i];

        tg_xz_xxyyyy[i] = tg_x_xxyyyy[i] * ra_z[i];

        tg_xz_xxyyyz[i] = 2.0 * tg_z_xyyyz[i] * fxi[i] + tg_z_xxyyyz[i] * ra_x[i];

        tg_xz_xxyyzz[i] = 2.0 * tg_z_xyyzz[i] * fxi[i] + tg_z_xxyyzz[i] * ra_x[i];

        tg_xz_xxyzzz[i] = 2.0 * tg_z_xyzzz[i] * fxi[i] + tg_z_xxyzzz[i] * ra_x[i];

        tg_xz_xxzzzz[i] = 2.0 * tg_z_xzzzz[i] * fxi[i] + tg_z_xxzzzz[i] * ra_x[i];

        tg_xz_xyyyyy[i] = tg_x_xyyyyy[i] * ra_z[i];

        tg_xz_xyyyyz[i] = tg_z_yyyyz[i] * fxi[i] + tg_z_xyyyyz[i] * ra_x[i];

        tg_xz_xyyyzz[i] = tg_z_yyyzz[i] * fxi[i] + tg_z_xyyyzz[i] * ra_x[i];

        tg_xz_xyyzzz[i] = tg_z_yyzzz[i] * fxi[i] + tg_z_xyyzzz[i] * ra_x[i];

        tg_xz_xyzzzz[i] = tg_z_yzzzz[i] * fxi[i] + tg_z_xyzzzz[i] * ra_x[i];

        tg_xz_xzzzzz[i] = tg_z_zzzzz[i] * fxi[i] + tg_z_xzzzzz[i] * ra_x[i];

        tg_xz_yyyyyy[i] = tg_z_yyyyyy[i] * ra_x[i];

        tg_xz_yyyyyz[i] = tg_z_yyyyyz[i] * ra_x[i];

        tg_xz_yyyyzz[i] = tg_z_yyyyzz[i] * ra_x[i];

        tg_xz_yyyzzz[i] = tg_z_yyyzzz[i] * ra_x[i];

        tg_xz_yyzzzz[i] = tg_z_yyzzzz[i] * ra_x[i];

        tg_xz_yzzzzz[i] = tg_z_yzzzzz[i] * ra_x[i];

        tg_xz_zzzzzz[i] = tg_z_zzzzzz[i] * ra_x[i];

        tg_yy_xxxxxx[i] = tg_0_xxxxxx[i] * fxi[i] + tg_y_xxxxxx[i] * ra_y[i];

        tg_yy_xxxxxy[i] = tg_0_xxxxxy[i] * fxi[i] + tg_y_xxxxx[i] * fxi[i] + tg_y_xxxxxy[i] * ra_y[i];

        tg_yy_xxxxxz[i] = tg_0_xxxxxz[i] * fxi[i] + tg_y_xxxxxz[i] * ra_y[i];

        tg_yy_xxxxyy[i] = tg_0_xxxxyy[i] * fxi[i] + 2.0 * tg_y_xxxxy[i] * fxi[i] + tg_y_xxxxyy[i] * ra_y[i];

        tg_yy_xxxxyz[i] = tg_0_xxxxyz[i] * fxi[i] + tg_y_xxxxz[i] * fxi[i] + tg_y_xxxxyz[i] * ra_y[i];

        tg_yy_xxxxzz[i] = tg_0_xxxxzz[i] * fxi[i] + tg_y_xxxxzz[i] * ra_y[i];

        tg_yy_xxxyyy[i] = tg_0_xxxyyy[i] * fxi[i] + 3.0 * tg_y_xxxyy[i] * fxi[i] + tg_y_xxxyyy[i] * ra_y[i];

        tg_yy_xxxyyz[i] = tg_0_xxxyyz[i] * fxi[i] + 2.0 * tg_y_xxxyz[i] * fxi[i] + tg_y_xxxyyz[i] * ra_y[i];

        tg_yy_xxxyzz[i] = tg_0_xxxyzz[i] * fxi[i] + tg_y_xxxzz[i] * fxi[i] + tg_y_xxxyzz[i] * ra_y[i];

        tg_yy_xxxzzz[i] = tg_0_xxxzzz[i] * fxi[i] + tg_y_xxxzzz[i] * ra_y[i];

        tg_yy_xxyyyy[i] = tg_0_xxyyyy[i] * fxi[i] + 4.0 * tg_y_xxyyy[i] * fxi[i] + tg_y_xxyyyy[i] * ra_y[i];

        tg_yy_xxyyyz[i] = tg_0_xxyyyz[i] * fxi[i] + 3.0 * tg_y_xxyyz[i] * fxi[i] + tg_y_xxyyyz[i] * ra_y[i];

        tg_yy_xxyyzz[i] = tg_0_xxyyzz[i] * fxi[i] + 2.0 * tg_y_xxyzz[i] * fxi[i] + tg_y_xxyyzz[i] * ra_y[i];

        tg_yy_xxyzzz[i] = tg_0_xxyzzz[i] * fxi[i] + tg_y_xxzzz[i] * fxi[i] + tg_y_xxyzzz[i] * ra_y[i];

        tg_yy_xxzzzz[i] = tg_0_xxzzzz[i] * fxi[i] + tg_y_xxzzzz[i] * ra_y[i];

        tg_yy_xyyyyy[i] = tg_0_xyyyyy[i] * fxi[i] + 5.0 * tg_y_xyyyy[i] * fxi[i] + tg_y_xyyyyy[i] * ra_y[i];

        tg_yy_xyyyyz[i] = tg_0_xyyyyz[i] * fxi[i] + 4.0 * tg_y_xyyyz[i] * fxi[i] + tg_y_xyyyyz[i] * ra_y[i];

        tg_yy_xyyyzz[i] = tg_0_xyyyzz[i] * fxi[i] + 3.0 * tg_y_xyyzz[i] * fxi[i] + tg_y_xyyyzz[i] * ra_y[i];

        tg_yy_xyyzzz[i] = tg_0_xyyzzz[i] * fxi[i] + 2.0 * tg_y_xyzzz[i] * fxi[i] + tg_y_xyyzzz[i] * ra_y[i];

        tg_yy_xyzzzz[i] = tg_0_xyzzzz[i] * fxi[i] + tg_y_xzzzz[i] * fxi[i] + tg_y_xyzzzz[i] * ra_y[i];

        tg_yy_xzzzzz[i] = tg_0_xzzzzz[i] * fxi[i] + tg_y_xzzzzz[i] * ra_y[i];

        tg_yy_yyyyyy[i] = tg_0_yyyyyy[i] * fxi[i] + 6.0 * tg_y_yyyyy[i] * fxi[i] + tg_y_yyyyyy[i] * ra_y[i];

        tg_yy_yyyyyz[i] = tg_0_yyyyyz[i] * fxi[i] + 5.0 * tg_y_yyyyz[i] * fxi[i] + tg_y_yyyyyz[i] * ra_y[i];

        tg_yy_yyyyzz[i] = tg_0_yyyyzz[i] * fxi[i] + 4.0 * tg_y_yyyzz[i] * fxi[i] + tg_y_yyyyzz[i] * ra_y[i];

        tg_yy_yyyzzz[i] = tg_0_yyyzzz[i] * fxi[i] + 3.0 * tg_y_yyzzz[i] * fxi[i] + tg_y_yyyzzz[i] * ra_y[i];

        tg_yy_yyzzzz[i] = tg_0_yyzzzz[i] * fxi[i] + 2.0 * tg_y_yzzzz[i] * fxi[i] + tg_y_yyzzzz[i] * ra_y[i];

        tg_yy_yzzzzz[i] = tg_0_yzzzzz[i] * fxi[i] + tg_y_zzzzz[i] * fxi[i] + tg_y_yzzzzz[i] * ra_y[i];

        tg_yy_zzzzzz[i] = tg_0_zzzzzz[i] * fxi[i] + tg_y_zzzzzz[i] * ra_y[i];

        tg_yz_xxxxxx[i] = tg_z_xxxxxx[i] * ra_y[i];

        tg_yz_xxxxxy[i] = tg_y_xxxxxy[i] * ra_z[i];

        tg_yz_xxxxxz[i] = tg_z_xxxxxz[i] * ra_y[i];

        tg_yz_xxxxyy[i] = tg_y_xxxxyy[i] * ra_z[i];

        tg_yz_xxxxyz[i] = tg_z_xxxxz[i] * fxi[i] + tg_z_xxxxyz[i] * ra_y[i];

        tg_yz_xxxxzz[i] = tg_z_xxxxzz[i] * ra_y[i];

        tg_yz_xxxyyy[i] = tg_y_xxxyyy[i] * ra_z[i];

        tg_yz_xxxyyz[i] = 2.0 * tg_z_xxxyz[i] * fxi[i] + tg_z_xxxyyz[i] * ra_y[i];

        tg_yz_xxxyzz[i] = tg_z_xxxzz[i] * fxi[i] + tg_z_xxxyzz[i] * ra_y[i];

        tg_yz_xxxzzz[i] = tg_z_xxxzzz[i] * ra_y[i];

        tg_yz_xxyyyy[i] = tg_y_xxyyyy[i] * ra_z[i];

        tg_yz_xxyyyz[i] = 3.0 * tg_z_xxyyz[i] * fxi[i] + tg_z_xxyyyz[i] * ra_y[i];

        tg_yz_xxyyzz[i] = 2.0 * tg_z_xxyzz[i] * fxi[i] + tg_z_xxyyzz[i] * ra_y[i];

        tg_yz_xxyzzz[i] = tg_z_xxzzz[i] * fxi[i] + tg_z_xxyzzz[i] * ra_y[i];

        tg_yz_xxzzzz[i] = tg_z_xxzzzz[i] * ra_y[i];

        tg_yz_xyyyyy[i] = tg_y_xyyyyy[i] * ra_z[i];

        tg_yz_xyyyyz[i] = 4.0 * tg_z_xyyyz[i] * fxi[i] + tg_z_xyyyyz[i] * ra_y[i];

        tg_yz_xyyyzz[i] = 3.0 * tg_z_xyyzz[i] * fxi[i] + tg_z_xyyyzz[i] * ra_y[i];

        tg_yz_xyyzzz[i] = 2.0 * tg_z_xyzzz[i] * fxi[i] + tg_z_xyyzzz[i] * ra_y[i];

        tg_yz_xyzzzz[i] = tg_z_xzzzz[i] * fxi[i] + tg_z_xyzzzz[i] * ra_y[i];

        tg_yz_xzzzzz[i] = tg_z_xzzzzz[i] * ra_y[i];

        tg_yz_yyyyyy[i] = tg_y_yyyyyy[i] * ra_z[i];

        tg_yz_yyyyyz[i] = 5.0 * tg_z_yyyyz[i] * fxi[i] + tg_z_yyyyyz[i] * ra_y[i];

        tg_yz_yyyyzz[i] = 4.0 * tg_z_yyyzz[i] * fxi[i] + tg_z_yyyyzz[i] * ra_y[i];

        tg_yz_yyyzzz[i] = 3.0 * tg_z_yyzzz[i] * fxi[i] + tg_z_yyyzzz[i] * ra_y[i];

        tg_yz_yyzzzz[i] = 2.0 * tg_z_yzzzz[i] * fxi[i] + tg_z_yyzzzz[i] * ra_y[i];

        tg_yz_yzzzzz[i] = tg_z_zzzzz[i] * fxi[i] + tg_z_yzzzzz[i] * ra_y[i];

        tg_yz_zzzzzz[i] = tg_z_zzzzzz[i] * ra_y[i];

        tg_zz_xxxxxx[i] = tg_0_xxxxxx[i] * fxi[i] + tg_z_xxxxxx[i] * ra_z[i];

        tg_zz_xxxxxy[i] = tg_0_xxxxxy[i] * fxi[i] + tg_z_xxxxxy[i] * ra_z[i];

        tg_zz_xxxxxz[i] = tg_0_xxxxxz[i] * fxi[i] + tg_z_xxxxx[i] * fxi[i] + tg_z_xxxxxz[i] * ra_z[i];

        tg_zz_xxxxyy[i] = tg_0_xxxxyy[i] * fxi[i] + tg_z_xxxxyy[i] * ra_z[i];

        tg_zz_xxxxyz[i] = tg_0_xxxxyz[i] * fxi[i] + tg_z_xxxxy[i] * fxi[i] + tg_z_xxxxyz[i] * ra_z[i];

        tg_zz_xxxxzz[i] = tg_0_xxxxzz[i] * fxi[i] + 2.0 * tg_z_xxxxz[i] * fxi[i] + tg_z_xxxxzz[i] * ra_z[i];

        tg_zz_xxxyyy[i] = tg_0_xxxyyy[i] * fxi[i] + tg_z_xxxyyy[i] * ra_z[i];

        tg_zz_xxxyyz[i] = tg_0_xxxyyz[i] * fxi[i] + tg_z_xxxyy[i] * fxi[i] + tg_z_xxxyyz[i] * ra_z[i];

        tg_zz_xxxyzz[i] = tg_0_xxxyzz[i] * fxi[i] + 2.0 * tg_z_xxxyz[i] * fxi[i] + tg_z_xxxyzz[i] * ra_z[i];

        tg_zz_xxxzzz[i] = tg_0_xxxzzz[i] * fxi[i] + 3.0 * tg_z_xxxzz[i] * fxi[i] + tg_z_xxxzzz[i] * ra_z[i];

        tg_zz_xxyyyy[i] = tg_0_xxyyyy[i] * fxi[i] + tg_z_xxyyyy[i] * ra_z[i];

        tg_zz_xxyyyz[i] = tg_0_xxyyyz[i] * fxi[i] + tg_z_xxyyy[i] * fxi[i] + tg_z_xxyyyz[i] * ra_z[i];

        tg_zz_xxyyzz[i] = tg_0_xxyyzz[i] * fxi[i] + 2.0 * tg_z_xxyyz[i] * fxi[i] + tg_z_xxyyzz[i] * ra_z[i];

        tg_zz_xxyzzz[i] = tg_0_xxyzzz[i] * fxi[i] + 3.0 * tg_z_xxyzz[i] * fxi[i] + tg_z_xxyzzz[i] * ra_z[i];

        tg_zz_xxzzzz[i] = tg_0_xxzzzz[i] * fxi[i] + 4.0 * tg_z_xxzzz[i] * fxi[i] + tg_z_xxzzzz[i] * ra_z[i];

        tg_zz_xyyyyy[i] = tg_0_xyyyyy[i] * fxi[i] + tg_z_xyyyyy[i] * ra_z[i];

        tg_zz_xyyyyz[i] = tg_0_xyyyyz[i] * fxi[i] + tg_z_xyyyy[i] * fxi[i] + tg_z_xyyyyz[i] * ra_z[i];

        tg_zz_xyyyzz[i] = tg_0_xyyyzz[i] * fxi[i] + 2.0 * tg_z_xyyyz[i] * fxi[i] + tg_z_xyyyzz[i] * ra_z[i];

        tg_zz_xyyzzz[i] = tg_0_xyyzzz[i] * fxi[i] + 3.0 * tg_z_xyyzz[i] * fxi[i] + tg_z_xyyzzz[i] * ra_z[i];

        tg_zz_xyzzzz[i] = tg_0_xyzzzz[i] * fxi[i] + 4.0 * tg_z_xyzzz[i] * fxi[i] + tg_z_xyzzzz[i] * ra_z[i];

        tg_zz_xzzzzz[i] = tg_0_xzzzzz[i] * fxi[i] + 5.0 * tg_z_xzzzz[i] * fxi[i] + tg_z_xzzzzz[i] * ra_z[i];

        tg_zz_yyyyyy[i] = tg_0_yyyyyy[i] * fxi[i] + tg_z_yyyyyy[i] * ra_z[i];

        tg_zz_yyyyyz[i] = tg_0_yyyyyz[i] * fxi[i] + tg_z_yyyyy[i] * fxi[i] + tg_z_yyyyyz[i] * ra_z[i];

        tg_zz_yyyyzz[i] = tg_0_yyyyzz[i] * fxi[i] + 2.0 * tg_z_yyyyz[i] * fxi[i] + tg_z_yyyyzz[i] * ra_z[i];

        tg_zz_yyyzzz[i] = tg_0_yyyzzz[i] * fxi[i] + 3.0 * tg_z_yyyzz[i] * fxi[i] + tg_z_yyyzzz[i] * ra_z[i];

        tg_zz_yyzzzz[i] = tg_0_yyzzzz[i] * fxi[i] + 4.0 * tg_z_yyzzz[i] * fxi[i] + tg_z_yyzzzz[i] * ra_z[i];

        tg_zz_yzzzzz[i] = tg_0_yzzzzz[i] * fxi[i] + 5.0 * tg_z_yzzzz[i] * fxi[i] + tg_z_yzzzzz[i] * ra_z[i];

        tg_zz_zzzzzz[i] = tg_0_zzzzzz[i] * fxi[i] + 6.0 * tg_z_zzzzz[i] * fxi[i] + tg_z_zzzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

