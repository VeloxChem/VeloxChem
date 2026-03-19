#include "LocalCorePotentialPrimRecFI.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_fi(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fi,
                                  const size_t idx_pi,
                                  const size_t idx_dh,
                                  const size_t idx_di,
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

    // Set up components of auxiliary buffer : DH

    auto tg_xx_xxxxx = pbuffer.data(idx_dh);

    auto tg_xx_xxxxy = pbuffer.data(idx_dh + 1);

    auto tg_xx_xxxxz = pbuffer.data(idx_dh + 2);

    auto tg_xx_xxxyy = pbuffer.data(idx_dh + 3);

    auto tg_xx_xxxyz = pbuffer.data(idx_dh + 4);

    auto tg_xx_xxxzz = pbuffer.data(idx_dh + 5);

    auto tg_xx_xxyyy = pbuffer.data(idx_dh + 6);

    auto tg_xx_xxyyz = pbuffer.data(idx_dh + 7);

    auto tg_xx_xxyzz = pbuffer.data(idx_dh + 8);

    auto tg_xx_xxzzz = pbuffer.data(idx_dh + 9);

    auto tg_xx_xyyyy = pbuffer.data(idx_dh + 10);

    auto tg_xx_xyyyz = pbuffer.data(idx_dh + 11);

    auto tg_xx_xyyzz = pbuffer.data(idx_dh + 12);

    auto tg_xx_xyzzz = pbuffer.data(idx_dh + 13);

    auto tg_xx_xzzzz = pbuffer.data(idx_dh + 14);

    auto tg_xx_yyyyy = pbuffer.data(idx_dh + 15);

    auto tg_xx_yyyyz = pbuffer.data(idx_dh + 16);

    auto tg_xx_yyyzz = pbuffer.data(idx_dh + 17);

    auto tg_xx_yyzzz = pbuffer.data(idx_dh + 18);

    auto tg_xx_yzzzz = pbuffer.data(idx_dh + 19);

    auto tg_xx_zzzzz = pbuffer.data(idx_dh + 20);

    auto tg_yy_xxxxx = pbuffer.data(idx_dh + 63);

    auto tg_yy_xxxxy = pbuffer.data(idx_dh + 64);

    auto tg_yy_xxxxz = pbuffer.data(idx_dh + 65);

    auto tg_yy_xxxyy = pbuffer.data(idx_dh + 66);

    auto tg_yy_xxxyz = pbuffer.data(idx_dh + 67);

    auto tg_yy_xxxzz = pbuffer.data(idx_dh + 68);

    auto tg_yy_xxyyy = pbuffer.data(idx_dh + 69);

    auto tg_yy_xxyyz = pbuffer.data(idx_dh + 70);

    auto tg_yy_xxyzz = pbuffer.data(idx_dh + 71);

    auto tg_yy_xxzzz = pbuffer.data(idx_dh + 72);

    auto tg_yy_xyyyy = pbuffer.data(idx_dh + 73);

    auto tg_yy_xyyyz = pbuffer.data(idx_dh + 74);

    auto tg_yy_xyyzz = pbuffer.data(idx_dh + 75);

    auto tg_yy_xyzzz = pbuffer.data(idx_dh + 76);

    auto tg_yy_xzzzz = pbuffer.data(idx_dh + 77);

    auto tg_yy_yyyyy = pbuffer.data(idx_dh + 78);

    auto tg_yy_yyyyz = pbuffer.data(idx_dh + 79);

    auto tg_yy_yyyzz = pbuffer.data(idx_dh + 80);

    auto tg_yy_yyzzz = pbuffer.data(idx_dh + 81);

    auto tg_yy_yzzzz = pbuffer.data(idx_dh + 82);

    auto tg_yy_zzzzz = pbuffer.data(idx_dh + 83);

    auto tg_yz_xxxyz = pbuffer.data(idx_dh + 88);

    auto tg_yz_xxyyz = pbuffer.data(idx_dh + 91);

    auto tg_yz_xxyzz = pbuffer.data(idx_dh + 92);

    auto tg_yz_xyyyz = pbuffer.data(idx_dh + 95);

    auto tg_yz_xyyzz = pbuffer.data(idx_dh + 96);

    auto tg_yz_xyzzz = pbuffer.data(idx_dh + 97);

    auto tg_yz_yyyyz = pbuffer.data(idx_dh + 100);

    auto tg_yz_yyyzz = pbuffer.data(idx_dh + 101);

    auto tg_yz_yyzzz = pbuffer.data(idx_dh + 102);

    auto tg_yz_yzzzz = pbuffer.data(idx_dh + 103);

    auto tg_zz_xxxxx = pbuffer.data(idx_dh + 105);

    auto tg_zz_xxxxy = pbuffer.data(idx_dh + 106);

    auto tg_zz_xxxxz = pbuffer.data(idx_dh + 107);

    auto tg_zz_xxxyy = pbuffer.data(idx_dh + 108);

    auto tg_zz_xxxyz = pbuffer.data(idx_dh + 109);

    auto tg_zz_xxxzz = pbuffer.data(idx_dh + 110);

    auto tg_zz_xxyyy = pbuffer.data(idx_dh + 111);

    auto tg_zz_xxyyz = pbuffer.data(idx_dh + 112);

    auto tg_zz_xxyzz = pbuffer.data(idx_dh + 113);

    auto tg_zz_xxzzz = pbuffer.data(idx_dh + 114);

    auto tg_zz_xyyyy = pbuffer.data(idx_dh + 115);

    auto tg_zz_xyyyz = pbuffer.data(idx_dh + 116);

    auto tg_zz_xyyzz = pbuffer.data(idx_dh + 117);

    auto tg_zz_xyzzz = pbuffer.data(idx_dh + 118);

    auto tg_zz_xzzzz = pbuffer.data(idx_dh + 119);

    auto tg_zz_yyyyy = pbuffer.data(idx_dh + 120);

    auto tg_zz_yyyyz = pbuffer.data(idx_dh + 121);

    auto tg_zz_yyyzz = pbuffer.data(idx_dh + 122);

    auto tg_zz_yyzzz = pbuffer.data(idx_dh + 123);

    auto tg_zz_yzzzz = pbuffer.data(idx_dh + 124);

    auto tg_zz_zzzzz = pbuffer.data(idx_dh + 125);

    // Set up components of auxiliary buffer : DI

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

    auto tg_xy_xxxxxy = pbuffer.data(idx_di + 29);

    auto tg_xy_xxxxyy = pbuffer.data(idx_di + 31);

    auto tg_xy_xxxyyy = pbuffer.data(idx_di + 34);

    auto tg_xy_xxyyyy = pbuffer.data(idx_di + 38);

    auto tg_xy_xyyyyy = pbuffer.data(idx_di + 43);

    auto tg_xy_yyyyyy = pbuffer.data(idx_di + 49);

    auto tg_xy_yyyyyz = pbuffer.data(idx_di + 50);

    auto tg_xy_yyyyzz = pbuffer.data(idx_di + 51);

    auto tg_xy_yyyzzz = pbuffer.data(idx_di + 52);

    auto tg_xy_yyzzzz = pbuffer.data(idx_di + 53);

    auto tg_xy_yzzzzz = pbuffer.data(idx_di + 54);

    auto tg_xz_xxxxxx = pbuffer.data(idx_di + 56);

    auto tg_xz_xxxxxz = pbuffer.data(idx_di + 58);

    auto tg_xz_xxxxzz = pbuffer.data(idx_di + 61);

    auto tg_xz_xxxzzz = pbuffer.data(idx_di + 65);

    auto tg_xz_xxzzzz = pbuffer.data(idx_di + 70);

    auto tg_xz_xzzzzz = pbuffer.data(idx_di + 76);

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

    auto tg_yz_xxxxxz = pbuffer.data(idx_di + 114);

    auto tg_yz_xxxxyz = pbuffer.data(idx_di + 116);

    auto tg_yz_xxxxzz = pbuffer.data(idx_di + 117);

    auto tg_yz_xxxyyz = pbuffer.data(idx_di + 119);

    auto tg_yz_xxxyzz = pbuffer.data(idx_di + 120);

    auto tg_yz_xxxzzz = pbuffer.data(idx_di + 121);

    auto tg_yz_xxyyyz = pbuffer.data(idx_di + 123);

    auto tg_yz_xxyyzz = pbuffer.data(idx_di + 124);

    auto tg_yz_xxyzzz = pbuffer.data(idx_di + 125);

    auto tg_yz_xxzzzz = pbuffer.data(idx_di + 126);

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

    // Set up components of targeted buffer : FI

    auto tg_xxx_xxxxxx = pbuffer.data(idx_fi);

    auto tg_xxx_xxxxxy = pbuffer.data(idx_fi + 1);

    auto tg_xxx_xxxxxz = pbuffer.data(idx_fi + 2);

    auto tg_xxx_xxxxyy = pbuffer.data(idx_fi + 3);

    auto tg_xxx_xxxxyz = pbuffer.data(idx_fi + 4);

    auto tg_xxx_xxxxzz = pbuffer.data(idx_fi + 5);

    auto tg_xxx_xxxyyy = pbuffer.data(idx_fi + 6);

    auto tg_xxx_xxxyyz = pbuffer.data(idx_fi + 7);

    auto tg_xxx_xxxyzz = pbuffer.data(idx_fi + 8);

    auto tg_xxx_xxxzzz = pbuffer.data(idx_fi + 9);

    auto tg_xxx_xxyyyy = pbuffer.data(idx_fi + 10);

    auto tg_xxx_xxyyyz = pbuffer.data(idx_fi + 11);

    auto tg_xxx_xxyyzz = pbuffer.data(idx_fi + 12);

    auto tg_xxx_xxyzzz = pbuffer.data(idx_fi + 13);

    auto tg_xxx_xxzzzz = pbuffer.data(idx_fi + 14);

    auto tg_xxx_xyyyyy = pbuffer.data(idx_fi + 15);

    auto tg_xxx_xyyyyz = pbuffer.data(idx_fi + 16);

    auto tg_xxx_xyyyzz = pbuffer.data(idx_fi + 17);

    auto tg_xxx_xyyzzz = pbuffer.data(idx_fi + 18);

    auto tg_xxx_xyzzzz = pbuffer.data(idx_fi + 19);

    auto tg_xxx_xzzzzz = pbuffer.data(idx_fi + 20);

    auto tg_xxx_yyyyyy = pbuffer.data(idx_fi + 21);

    auto tg_xxx_yyyyyz = pbuffer.data(idx_fi + 22);

    auto tg_xxx_yyyyzz = pbuffer.data(idx_fi + 23);

    auto tg_xxx_yyyzzz = pbuffer.data(idx_fi + 24);

    auto tg_xxx_yyzzzz = pbuffer.data(idx_fi + 25);

    auto tg_xxx_yzzzzz = pbuffer.data(idx_fi + 26);

    auto tg_xxx_zzzzzz = pbuffer.data(idx_fi + 27);

    auto tg_xxy_xxxxxx = pbuffer.data(idx_fi + 28);

    auto tg_xxy_xxxxxy = pbuffer.data(idx_fi + 29);

    auto tg_xxy_xxxxxz = pbuffer.data(idx_fi + 30);

    auto tg_xxy_xxxxyy = pbuffer.data(idx_fi + 31);

    auto tg_xxy_xxxxyz = pbuffer.data(idx_fi + 32);

    auto tg_xxy_xxxxzz = pbuffer.data(idx_fi + 33);

    auto tg_xxy_xxxyyy = pbuffer.data(idx_fi + 34);

    auto tg_xxy_xxxyyz = pbuffer.data(idx_fi + 35);

    auto tg_xxy_xxxyzz = pbuffer.data(idx_fi + 36);

    auto tg_xxy_xxxzzz = pbuffer.data(idx_fi + 37);

    auto tg_xxy_xxyyyy = pbuffer.data(idx_fi + 38);

    auto tg_xxy_xxyyyz = pbuffer.data(idx_fi + 39);

    auto tg_xxy_xxyyzz = pbuffer.data(idx_fi + 40);

    auto tg_xxy_xxyzzz = pbuffer.data(idx_fi + 41);

    auto tg_xxy_xxzzzz = pbuffer.data(idx_fi + 42);

    auto tg_xxy_xyyyyy = pbuffer.data(idx_fi + 43);

    auto tg_xxy_xyyyyz = pbuffer.data(idx_fi + 44);

    auto tg_xxy_xyyyzz = pbuffer.data(idx_fi + 45);

    auto tg_xxy_xyyzzz = pbuffer.data(idx_fi + 46);

    auto tg_xxy_xyzzzz = pbuffer.data(idx_fi + 47);

    auto tg_xxy_xzzzzz = pbuffer.data(idx_fi + 48);

    auto tg_xxy_yyyyyy = pbuffer.data(idx_fi + 49);

    auto tg_xxy_yyyyyz = pbuffer.data(idx_fi + 50);

    auto tg_xxy_yyyyzz = pbuffer.data(idx_fi + 51);

    auto tg_xxy_yyyzzz = pbuffer.data(idx_fi + 52);

    auto tg_xxy_yyzzzz = pbuffer.data(idx_fi + 53);

    auto tg_xxy_yzzzzz = pbuffer.data(idx_fi + 54);

    auto tg_xxy_zzzzzz = pbuffer.data(idx_fi + 55);

    auto tg_xxz_xxxxxx = pbuffer.data(idx_fi + 56);

    auto tg_xxz_xxxxxy = pbuffer.data(idx_fi + 57);

    auto tg_xxz_xxxxxz = pbuffer.data(idx_fi + 58);

    auto tg_xxz_xxxxyy = pbuffer.data(idx_fi + 59);

    auto tg_xxz_xxxxyz = pbuffer.data(idx_fi + 60);

    auto tg_xxz_xxxxzz = pbuffer.data(idx_fi + 61);

    auto tg_xxz_xxxyyy = pbuffer.data(idx_fi + 62);

    auto tg_xxz_xxxyyz = pbuffer.data(idx_fi + 63);

    auto tg_xxz_xxxyzz = pbuffer.data(idx_fi + 64);

    auto tg_xxz_xxxzzz = pbuffer.data(idx_fi + 65);

    auto tg_xxz_xxyyyy = pbuffer.data(idx_fi + 66);

    auto tg_xxz_xxyyyz = pbuffer.data(idx_fi + 67);

    auto tg_xxz_xxyyzz = pbuffer.data(idx_fi + 68);

    auto tg_xxz_xxyzzz = pbuffer.data(idx_fi + 69);

    auto tg_xxz_xxzzzz = pbuffer.data(idx_fi + 70);

    auto tg_xxz_xyyyyy = pbuffer.data(idx_fi + 71);

    auto tg_xxz_xyyyyz = pbuffer.data(idx_fi + 72);

    auto tg_xxz_xyyyzz = pbuffer.data(idx_fi + 73);

    auto tg_xxz_xyyzzz = pbuffer.data(idx_fi + 74);

    auto tg_xxz_xyzzzz = pbuffer.data(idx_fi + 75);

    auto tg_xxz_xzzzzz = pbuffer.data(idx_fi + 76);

    auto tg_xxz_yyyyyy = pbuffer.data(idx_fi + 77);

    auto tg_xxz_yyyyyz = pbuffer.data(idx_fi + 78);

    auto tg_xxz_yyyyzz = pbuffer.data(idx_fi + 79);

    auto tg_xxz_yyyzzz = pbuffer.data(idx_fi + 80);

    auto tg_xxz_yyzzzz = pbuffer.data(idx_fi + 81);

    auto tg_xxz_yzzzzz = pbuffer.data(idx_fi + 82);

    auto tg_xxz_zzzzzz = pbuffer.data(idx_fi + 83);

    auto tg_xyy_xxxxxx = pbuffer.data(idx_fi + 84);

    auto tg_xyy_xxxxxy = pbuffer.data(idx_fi + 85);

    auto tg_xyy_xxxxxz = pbuffer.data(idx_fi + 86);

    auto tg_xyy_xxxxyy = pbuffer.data(idx_fi + 87);

    auto tg_xyy_xxxxyz = pbuffer.data(idx_fi + 88);

    auto tg_xyy_xxxxzz = pbuffer.data(idx_fi + 89);

    auto tg_xyy_xxxyyy = pbuffer.data(idx_fi + 90);

    auto tg_xyy_xxxyyz = pbuffer.data(idx_fi + 91);

    auto tg_xyy_xxxyzz = pbuffer.data(idx_fi + 92);

    auto tg_xyy_xxxzzz = pbuffer.data(idx_fi + 93);

    auto tg_xyy_xxyyyy = pbuffer.data(idx_fi + 94);

    auto tg_xyy_xxyyyz = pbuffer.data(idx_fi + 95);

    auto tg_xyy_xxyyzz = pbuffer.data(idx_fi + 96);

    auto tg_xyy_xxyzzz = pbuffer.data(idx_fi + 97);

    auto tg_xyy_xxzzzz = pbuffer.data(idx_fi + 98);

    auto tg_xyy_xyyyyy = pbuffer.data(idx_fi + 99);

    auto tg_xyy_xyyyyz = pbuffer.data(idx_fi + 100);

    auto tg_xyy_xyyyzz = pbuffer.data(idx_fi + 101);

    auto tg_xyy_xyyzzz = pbuffer.data(idx_fi + 102);

    auto tg_xyy_xyzzzz = pbuffer.data(idx_fi + 103);

    auto tg_xyy_xzzzzz = pbuffer.data(idx_fi + 104);

    auto tg_xyy_yyyyyy = pbuffer.data(idx_fi + 105);

    auto tg_xyy_yyyyyz = pbuffer.data(idx_fi + 106);

    auto tg_xyy_yyyyzz = pbuffer.data(idx_fi + 107);

    auto tg_xyy_yyyzzz = pbuffer.data(idx_fi + 108);

    auto tg_xyy_yyzzzz = pbuffer.data(idx_fi + 109);

    auto tg_xyy_yzzzzz = pbuffer.data(idx_fi + 110);

    auto tg_xyy_zzzzzz = pbuffer.data(idx_fi + 111);

    auto tg_xyz_xxxxxx = pbuffer.data(idx_fi + 112);

    auto tg_xyz_xxxxxy = pbuffer.data(idx_fi + 113);

    auto tg_xyz_xxxxxz = pbuffer.data(idx_fi + 114);

    auto tg_xyz_xxxxyy = pbuffer.data(idx_fi + 115);

    auto tg_xyz_xxxxyz = pbuffer.data(idx_fi + 116);

    auto tg_xyz_xxxxzz = pbuffer.data(idx_fi + 117);

    auto tg_xyz_xxxyyy = pbuffer.data(idx_fi + 118);

    auto tg_xyz_xxxyyz = pbuffer.data(idx_fi + 119);

    auto tg_xyz_xxxyzz = pbuffer.data(idx_fi + 120);

    auto tg_xyz_xxxzzz = pbuffer.data(idx_fi + 121);

    auto tg_xyz_xxyyyy = pbuffer.data(idx_fi + 122);

    auto tg_xyz_xxyyyz = pbuffer.data(idx_fi + 123);

    auto tg_xyz_xxyyzz = pbuffer.data(idx_fi + 124);

    auto tg_xyz_xxyzzz = pbuffer.data(idx_fi + 125);

    auto tg_xyz_xxzzzz = pbuffer.data(idx_fi + 126);

    auto tg_xyz_xyyyyy = pbuffer.data(idx_fi + 127);

    auto tg_xyz_xyyyyz = pbuffer.data(idx_fi + 128);

    auto tg_xyz_xyyyzz = pbuffer.data(idx_fi + 129);

    auto tg_xyz_xyyzzz = pbuffer.data(idx_fi + 130);

    auto tg_xyz_xyzzzz = pbuffer.data(idx_fi + 131);

    auto tg_xyz_xzzzzz = pbuffer.data(idx_fi + 132);

    auto tg_xyz_yyyyyy = pbuffer.data(idx_fi + 133);

    auto tg_xyz_yyyyyz = pbuffer.data(idx_fi + 134);

    auto tg_xyz_yyyyzz = pbuffer.data(idx_fi + 135);

    auto tg_xyz_yyyzzz = pbuffer.data(idx_fi + 136);

    auto tg_xyz_yyzzzz = pbuffer.data(idx_fi + 137);

    auto tg_xyz_yzzzzz = pbuffer.data(idx_fi + 138);

    auto tg_xyz_zzzzzz = pbuffer.data(idx_fi + 139);

    auto tg_xzz_xxxxxx = pbuffer.data(idx_fi + 140);

    auto tg_xzz_xxxxxy = pbuffer.data(idx_fi + 141);

    auto tg_xzz_xxxxxz = pbuffer.data(idx_fi + 142);

    auto tg_xzz_xxxxyy = pbuffer.data(idx_fi + 143);

    auto tg_xzz_xxxxyz = pbuffer.data(idx_fi + 144);

    auto tg_xzz_xxxxzz = pbuffer.data(idx_fi + 145);

    auto tg_xzz_xxxyyy = pbuffer.data(idx_fi + 146);

    auto tg_xzz_xxxyyz = pbuffer.data(idx_fi + 147);

    auto tg_xzz_xxxyzz = pbuffer.data(idx_fi + 148);

    auto tg_xzz_xxxzzz = pbuffer.data(idx_fi + 149);

    auto tg_xzz_xxyyyy = pbuffer.data(idx_fi + 150);

    auto tg_xzz_xxyyyz = pbuffer.data(idx_fi + 151);

    auto tg_xzz_xxyyzz = pbuffer.data(idx_fi + 152);

    auto tg_xzz_xxyzzz = pbuffer.data(idx_fi + 153);

    auto tg_xzz_xxzzzz = pbuffer.data(idx_fi + 154);

    auto tg_xzz_xyyyyy = pbuffer.data(idx_fi + 155);

    auto tg_xzz_xyyyyz = pbuffer.data(idx_fi + 156);

    auto tg_xzz_xyyyzz = pbuffer.data(idx_fi + 157);

    auto tg_xzz_xyyzzz = pbuffer.data(idx_fi + 158);

    auto tg_xzz_xyzzzz = pbuffer.data(idx_fi + 159);

    auto tg_xzz_xzzzzz = pbuffer.data(idx_fi + 160);

    auto tg_xzz_yyyyyy = pbuffer.data(idx_fi + 161);

    auto tg_xzz_yyyyyz = pbuffer.data(idx_fi + 162);

    auto tg_xzz_yyyyzz = pbuffer.data(idx_fi + 163);

    auto tg_xzz_yyyzzz = pbuffer.data(idx_fi + 164);

    auto tg_xzz_yyzzzz = pbuffer.data(idx_fi + 165);

    auto tg_xzz_yzzzzz = pbuffer.data(idx_fi + 166);

    auto tg_xzz_zzzzzz = pbuffer.data(idx_fi + 167);

    auto tg_yyy_xxxxxx = pbuffer.data(idx_fi + 168);

    auto tg_yyy_xxxxxy = pbuffer.data(idx_fi + 169);

    auto tg_yyy_xxxxxz = pbuffer.data(idx_fi + 170);

    auto tg_yyy_xxxxyy = pbuffer.data(idx_fi + 171);

    auto tg_yyy_xxxxyz = pbuffer.data(idx_fi + 172);

    auto tg_yyy_xxxxzz = pbuffer.data(idx_fi + 173);

    auto tg_yyy_xxxyyy = pbuffer.data(idx_fi + 174);

    auto tg_yyy_xxxyyz = pbuffer.data(idx_fi + 175);

    auto tg_yyy_xxxyzz = pbuffer.data(idx_fi + 176);

    auto tg_yyy_xxxzzz = pbuffer.data(idx_fi + 177);

    auto tg_yyy_xxyyyy = pbuffer.data(idx_fi + 178);

    auto tg_yyy_xxyyyz = pbuffer.data(idx_fi + 179);

    auto tg_yyy_xxyyzz = pbuffer.data(idx_fi + 180);

    auto tg_yyy_xxyzzz = pbuffer.data(idx_fi + 181);

    auto tg_yyy_xxzzzz = pbuffer.data(idx_fi + 182);

    auto tg_yyy_xyyyyy = pbuffer.data(idx_fi + 183);

    auto tg_yyy_xyyyyz = pbuffer.data(idx_fi + 184);

    auto tg_yyy_xyyyzz = pbuffer.data(idx_fi + 185);

    auto tg_yyy_xyyzzz = pbuffer.data(idx_fi + 186);

    auto tg_yyy_xyzzzz = pbuffer.data(idx_fi + 187);

    auto tg_yyy_xzzzzz = pbuffer.data(idx_fi + 188);

    auto tg_yyy_yyyyyy = pbuffer.data(idx_fi + 189);

    auto tg_yyy_yyyyyz = pbuffer.data(idx_fi + 190);

    auto tg_yyy_yyyyzz = pbuffer.data(idx_fi + 191);

    auto tg_yyy_yyyzzz = pbuffer.data(idx_fi + 192);

    auto tg_yyy_yyzzzz = pbuffer.data(idx_fi + 193);

    auto tg_yyy_yzzzzz = pbuffer.data(idx_fi + 194);

    auto tg_yyy_zzzzzz = pbuffer.data(idx_fi + 195);

    auto tg_yyz_xxxxxx = pbuffer.data(idx_fi + 196);

    auto tg_yyz_xxxxxy = pbuffer.data(idx_fi + 197);

    auto tg_yyz_xxxxxz = pbuffer.data(idx_fi + 198);

    auto tg_yyz_xxxxyy = pbuffer.data(idx_fi + 199);

    auto tg_yyz_xxxxyz = pbuffer.data(idx_fi + 200);

    auto tg_yyz_xxxxzz = pbuffer.data(idx_fi + 201);

    auto tg_yyz_xxxyyy = pbuffer.data(idx_fi + 202);

    auto tg_yyz_xxxyyz = pbuffer.data(idx_fi + 203);

    auto tg_yyz_xxxyzz = pbuffer.data(idx_fi + 204);

    auto tg_yyz_xxxzzz = pbuffer.data(idx_fi + 205);

    auto tg_yyz_xxyyyy = pbuffer.data(idx_fi + 206);

    auto tg_yyz_xxyyyz = pbuffer.data(idx_fi + 207);

    auto tg_yyz_xxyyzz = pbuffer.data(idx_fi + 208);

    auto tg_yyz_xxyzzz = pbuffer.data(idx_fi + 209);

    auto tg_yyz_xxzzzz = pbuffer.data(idx_fi + 210);

    auto tg_yyz_xyyyyy = pbuffer.data(idx_fi + 211);

    auto tg_yyz_xyyyyz = pbuffer.data(idx_fi + 212);

    auto tg_yyz_xyyyzz = pbuffer.data(idx_fi + 213);

    auto tg_yyz_xyyzzz = pbuffer.data(idx_fi + 214);

    auto tg_yyz_xyzzzz = pbuffer.data(idx_fi + 215);

    auto tg_yyz_xzzzzz = pbuffer.data(idx_fi + 216);

    auto tg_yyz_yyyyyy = pbuffer.data(idx_fi + 217);

    auto tg_yyz_yyyyyz = pbuffer.data(idx_fi + 218);

    auto tg_yyz_yyyyzz = pbuffer.data(idx_fi + 219);

    auto tg_yyz_yyyzzz = pbuffer.data(idx_fi + 220);

    auto tg_yyz_yyzzzz = pbuffer.data(idx_fi + 221);

    auto tg_yyz_yzzzzz = pbuffer.data(idx_fi + 222);

    auto tg_yyz_zzzzzz = pbuffer.data(idx_fi + 223);

    auto tg_yzz_xxxxxx = pbuffer.data(idx_fi + 224);

    auto tg_yzz_xxxxxy = pbuffer.data(idx_fi + 225);

    auto tg_yzz_xxxxxz = pbuffer.data(idx_fi + 226);

    auto tg_yzz_xxxxyy = pbuffer.data(idx_fi + 227);

    auto tg_yzz_xxxxyz = pbuffer.data(idx_fi + 228);

    auto tg_yzz_xxxxzz = pbuffer.data(idx_fi + 229);

    auto tg_yzz_xxxyyy = pbuffer.data(idx_fi + 230);

    auto tg_yzz_xxxyyz = pbuffer.data(idx_fi + 231);

    auto tg_yzz_xxxyzz = pbuffer.data(idx_fi + 232);

    auto tg_yzz_xxxzzz = pbuffer.data(idx_fi + 233);

    auto tg_yzz_xxyyyy = pbuffer.data(idx_fi + 234);

    auto tg_yzz_xxyyyz = pbuffer.data(idx_fi + 235);

    auto tg_yzz_xxyyzz = pbuffer.data(idx_fi + 236);

    auto tg_yzz_xxyzzz = pbuffer.data(idx_fi + 237);

    auto tg_yzz_xxzzzz = pbuffer.data(idx_fi + 238);

    auto tg_yzz_xyyyyy = pbuffer.data(idx_fi + 239);

    auto tg_yzz_xyyyyz = pbuffer.data(idx_fi + 240);

    auto tg_yzz_xyyyzz = pbuffer.data(idx_fi + 241);

    auto tg_yzz_xyyzzz = pbuffer.data(idx_fi + 242);

    auto tg_yzz_xyzzzz = pbuffer.data(idx_fi + 243);

    auto tg_yzz_xzzzzz = pbuffer.data(idx_fi + 244);

    auto tg_yzz_yyyyyy = pbuffer.data(idx_fi + 245);

    auto tg_yzz_yyyyyz = pbuffer.data(idx_fi + 246);

    auto tg_yzz_yyyyzz = pbuffer.data(idx_fi + 247);

    auto tg_yzz_yyyzzz = pbuffer.data(idx_fi + 248);

    auto tg_yzz_yyzzzz = pbuffer.data(idx_fi + 249);

    auto tg_yzz_yzzzzz = pbuffer.data(idx_fi + 250);

    auto tg_yzz_zzzzzz = pbuffer.data(idx_fi + 251);

    auto tg_zzz_xxxxxx = pbuffer.data(idx_fi + 252);

    auto tg_zzz_xxxxxy = pbuffer.data(idx_fi + 253);

    auto tg_zzz_xxxxxz = pbuffer.data(idx_fi + 254);

    auto tg_zzz_xxxxyy = pbuffer.data(idx_fi + 255);

    auto tg_zzz_xxxxyz = pbuffer.data(idx_fi + 256);

    auto tg_zzz_xxxxzz = pbuffer.data(idx_fi + 257);

    auto tg_zzz_xxxyyy = pbuffer.data(idx_fi + 258);

    auto tg_zzz_xxxyyz = pbuffer.data(idx_fi + 259);

    auto tg_zzz_xxxyzz = pbuffer.data(idx_fi + 260);

    auto tg_zzz_xxxzzz = pbuffer.data(idx_fi + 261);

    auto tg_zzz_xxyyyy = pbuffer.data(idx_fi + 262);

    auto tg_zzz_xxyyyz = pbuffer.data(idx_fi + 263);

    auto tg_zzz_xxyyzz = pbuffer.data(idx_fi + 264);

    auto tg_zzz_xxyzzz = pbuffer.data(idx_fi + 265);

    auto tg_zzz_xxzzzz = pbuffer.data(idx_fi + 266);

    auto tg_zzz_xyyyyy = pbuffer.data(idx_fi + 267);

    auto tg_zzz_xyyyyz = pbuffer.data(idx_fi + 268);

    auto tg_zzz_xyyyzz = pbuffer.data(idx_fi + 269);

    auto tg_zzz_xyyzzz = pbuffer.data(idx_fi + 270);

    auto tg_zzz_xyzzzz = pbuffer.data(idx_fi + 271);

    auto tg_zzz_xzzzzz = pbuffer.data(idx_fi + 272);

    auto tg_zzz_yyyyyy = pbuffer.data(idx_fi + 273);

    auto tg_zzz_yyyyyz = pbuffer.data(idx_fi + 274);

    auto tg_zzz_yyyyzz = pbuffer.data(idx_fi + 275);

    auto tg_zzz_yyyzzz = pbuffer.data(idx_fi + 276);

    auto tg_zzz_yyzzzz = pbuffer.data(idx_fi + 277);

    auto tg_zzz_yzzzzz = pbuffer.data(idx_fi + 278);

    auto tg_zzz_zzzzzz = pbuffer.data(idx_fi + 279);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_xxxxxx, tg_x_xxxxxy, tg_x_xxxxxz, tg_x_xxxxyy, tg_x_xxxxyz, tg_x_xxxxzz, tg_x_xxxyyy, tg_x_xxxyyz, tg_x_xxxyzz, tg_x_xxxzzz, tg_x_xxyyyy, tg_x_xxyyyz, tg_x_xxyyzz, tg_x_xxyzzz, tg_x_xxzzzz, tg_x_xyyyyy, tg_x_xyyyyz, tg_x_xyyyzz, tg_x_xyyzzz, tg_x_xyzzzz, tg_x_xzzzzz, tg_x_yyyyyy, tg_x_yyyyyz, tg_x_yyyyzz, tg_x_yyyzzz, tg_x_yyzzzz, tg_x_yzzzzz, tg_x_zzzzzz, tg_xx_xxxxx, tg_xx_xxxxxx, tg_xx_xxxxxy, tg_xx_xxxxxz, tg_xx_xxxxy, tg_xx_xxxxyy, tg_xx_xxxxyz, tg_xx_xxxxz, tg_xx_xxxxzz, tg_xx_xxxyy, tg_xx_xxxyyy, tg_xx_xxxyyz, tg_xx_xxxyz, tg_xx_xxxyzz, tg_xx_xxxzz, tg_xx_xxxzzz, tg_xx_xxyyy, tg_xx_xxyyyy, tg_xx_xxyyyz, tg_xx_xxyyz, tg_xx_xxyyzz, tg_xx_xxyzz, tg_xx_xxyzzz, tg_xx_xxzzz, tg_xx_xxzzzz, tg_xx_xyyyy, tg_xx_xyyyyy, tg_xx_xyyyyz, tg_xx_xyyyz, tg_xx_xyyyzz, tg_xx_xyyzz, tg_xx_xyyzzz, tg_xx_xyzzz, tg_xx_xyzzzz, tg_xx_xzzzz, tg_xx_xzzzzz, tg_xx_yyyyy, tg_xx_yyyyyy, tg_xx_yyyyyz, tg_xx_yyyyz, tg_xx_yyyyzz, tg_xx_yyyzz, tg_xx_yyyzzz, tg_xx_yyzzz, tg_xx_yyzzzz, tg_xx_yzzzz, tg_xx_yzzzzz, tg_xx_zzzzz, tg_xx_zzzzzz, tg_xxx_xxxxxx, tg_xxx_xxxxxy, tg_xxx_xxxxxz, tg_xxx_xxxxyy, tg_xxx_xxxxyz, tg_xxx_xxxxzz, tg_xxx_xxxyyy, tg_xxx_xxxyyz, tg_xxx_xxxyzz, tg_xxx_xxxzzz, tg_xxx_xxyyyy, tg_xxx_xxyyyz, tg_xxx_xxyyzz, tg_xxx_xxyzzz, tg_xxx_xxzzzz, tg_xxx_xyyyyy, tg_xxx_xyyyyz, tg_xxx_xyyyzz, tg_xxx_xyyzzz, tg_xxx_xyzzzz, tg_xxx_xzzzzz, tg_xxx_yyyyyy, tg_xxx_yyyyyz, tg_xxx_yyyyzz, tg_xxx_yyyzzz, tg_xxx_yyzzzz, tg_xxx_yzzzzz, tg_xxx_zzzzzz, tg_xxy_xxxxxx, tg_xxy_xxxxxy, tg_xxy_xxxxxz, tg_xxy_xxxxyy, tg_xxy_xxxxyz, tg_xxy_xxxxzz, tg_xxy_xxxyyy, tg_xxy_xxxyyz, tg_xxy_xxxyzz, tg_xxy_xxxzzz, tg_xxy_xxyyyy, tg_xxy_xxyyyz, tg_xxy_xxyyzz, tg_xxy_xxyzzz, tg_xxy_xxzzzz, tg_xxy_xyyyyy, tg_xxy_xyyyyz, tg_xxy_xyyyzz, tg_xxy_xyyzzz, tg_xxy_xyzzzz, tg_xxy_xzzzzz, tg_xxy_yyyyyy, tg_xxy_yyyyyz, tg_xxy_yyyyzz, tg_xxy_yyyzzz, tg_xxy_yyzzzz, tg_xxy_yzzzzz, tg_xxy_zzzzzz, tg_xxz_xxxxxx, tg_xxz_xxxxxy, tg_xxz_xxxxxz, tg_xxz_xxxxyy, tg_xxz_xxxxyz, tg_xxz_xxxxzz, tg_xxz_xxxyyy, tg_xxz_xxxyyz, tg_xxz_xxxyzz, tg_xxz_xxxzzz, tg_xxz_xxyyyy, tg_xxz_xxyyyz, tg_xxz_xxyyzz, tg_xxz_xxyzzz, tg_xxz_xxzzzz, tg_xxz_xyyyyy, tg_xxz_xyyyyz, tg_xxz_xyyyzz, tg_xxz_xyyzzz, tg_xxz_xyzzzz, tg_xxz_xzzzzz, tg_xxz_yyyyyy, tg_xxz_yyyyyz, tg_xxz_yyyyzz, tg_xxz_yyyzzz, tg_xxz_yyzzzz, tg_xxz_yzzzzz, tg_xxz_zzzzzz, tg_xy_xxxxxy, tg_xy_xxxxyy, tg_xy_xxxyyy, tg_xy_xxyyyy, tg_xy_xyyyyy, tg_xy_yyyyyy, tg_xy_yyyyyz, tg_xy_yyyyzz, tg_xy_yyyzzz, tg_xy_yyzzzz, tg_xy_yzzzzz, tg_xyy_xxxxxx, tg_xyy_xxxxxy, tg_xyy_xxxxxz, tg_xyy_xxxxyy, tg_xyy_xxxxyz, tg_xyy_xxxxzz, tg_xyy_xxxyyy, tg_xyy_xxxyyz, tg_xyy_xxxyzz, tg_xyy_xxxzzz, tg_xyy_xxyyyy, tg_xyy_xxyyyz, tg_xyy_xxyyzz, tg_xyy_xxyzzz, tg_xyy_xxzzzz, tg_xyy_xyyyyy, tg_xyy_xyyyyz, tg_xyy_xyyyzz, tg_xyy_xyyzzz, tg_xyy_xyzzzz, tg_xyy_xzzzzz, tg_xyy_yyyyyy, tg_xyy_yyyyyz, tg_xyy_yyyyzz, tg_xyy_yyyzzz, tg_xyy_yyzzzz, tg_xyy_yzzzzz, tg_xyy_zzzzzz, tg_xyz_xxxxxx, tg_xyz_xxxxxy, tg_xyz_xxxxxz, tg_xyz_xxxxyy, tg_xyz_xxxxyz, tg_xyz_xxxxzz, tg_xyz_xxxyyy, tg_xyz_xxxyyz, tg_xyz_xxxyzz, tg_xyz_xxxzzz, tg_xyz_xxyyyy, tg_xyz_xxyyyz, tg_xyz_xxyyzz, tg_xyz_xxyzzz, tg_xyz_xxzzzz, tg_xyz_xyyyyy, tg_xyz_xyyyyz, tg_xyz_xyyyzz, tg_xyz_xyyzzz, tg_xyz_xyzzzz, tg_xyz_xzzzzz, tg_xyz_yyyyyy, tg_xyz_yyyyyz, tg_xyz_yyyyzz, tg_xyz_yyyzzz, tg_xyz_yyzzzz, tg_xyz_yzzzzz, tg_xyz_zzzzzz, tg_xz_xxxxxx, tg_xz_xxxxxz, tg_xz_xxxxzz, tg_xz_xxxzzz, tg_xz_xxzzzz, tg_xz_xzzzzz, tg_xz_yyyyyz, tg_xz_yyyyzz, tg_xz_yyyzzz, tg_xz_yyzzzz, tg_xz_yzzzzz, tg_xz_zzzzzz, tg_xzz_xxxxxx, tg_xzz_xxxxxy, tg_xzz_xxxxxz, tg_xzz_xxxxyy, tg_xzz_xxxxyz, tg_xzz_xxxxzz, tg_xzz_xxxyyy, tg_xzz_xxxyyz, tg_xzz_xxxyzz, tg_xzz_xxxzzz, tg_xzz_xxyyyy, tg_xzz_xxyyyz, tg_xzz_xxyyzz, tg_xzz_xxyzzz, tg_xzz_xxzzzz, tg_xzz_xyyyyy, tg_xzz_xyyyyz, tg_xzz_xyyyzz, tg_xzz_xyyzzz, tg_xzz_xyzzzz, tg_xzz_xzzzzz, tg_xzz_yyyyyy, tg_xzz_yyyyyz, tg_xzz_yyyyzz, tg_xzz_yyyzzz, tg_xzz_yyzzzz, tg_xzz_yzzzzz, tg_xzz_zzzzzz, tg_y_xxxxxx, tg_y_xxxxxy, tg_y_xxxxxz, tg_y_xxxxyy, tg_y_xxxxyz, tg_y_xxxxzz, tg_y_xxxyyy, tg_y_xxxyyz, tg_y_xxxyzz, tg_y_xxxzzz, tg_y_xxyyyy, tg_y_xxyyyz, tg_y_xxyyzz, tg_y_xxyzzz, tg_y_xxzzzz, tg_y_xyyyyy, tg_y_xyyyyz, tg_y_xyyyzz, tg_y_xyyzzz, tg_y_xyzzzz, tg_y_xzzzzz, tg_y_yyyyyy, tg_y_yyyyyz, tg_y_yyyyzz, tg_y_yyyzzz, tg_y_yyzzzz, tg_y_yzzzzz, tg_y_zzzzzz, tg_yy_xxxxx, tg_yy_xxxxxx, tg_yy_xxxxxy, tg_yy_xxxxxz, tg_yy_xxxxy, tg_yy_xxxxyy, tg_yy_xxxxyz, tg_yy_xxxxz, tg_yy_xxxxzz, tg_yy_xxxyy, tg_yy_xxxyyy, tg_yy_xxxyyz, tg_yy_xxxyz, tg_yy_xxxyzz, tg_yy_xxxzz, tg_yy_xxxzzz, tg_yy_xxyyy, tg_yy_xxyyyy, tg_yy_xxyyyz, tg_yy_xxyyz, tg_yy_xxyyzz, tg_yy_xxyzz, tg_yy_xxyzzz, tg_yy_xxzzz, tg_yy_xxzzzz, tg_yy_xyyyy, tg_yy_xyyyyy, tg_yy_xyyyyz, tg_yy_xyyyz, tg_yy_xyyyzz, tg_yy_xyyzz, tg_yy_xyyzzz, tg_yy_xyzzz, tg_yy_xyzzzz, tg_yy_xzzzz, tg_yy_xzzzzz, tg_yy_yyyyy, tg_yy_yyyyyy, tg_yy_yyyyyz, tg_yy_yyyyz, tg_yy_yyyyzz, tg_yy_yyyzz, tg_yy_yyyzzz, tg_yy_yyzzz, tg_yy_yyzzzz, tg_yy_yzzzz, tg_yy_yzzzzz, tg_yy_zzzzz, tg_yy_zzzzzz, tg_yyy_xxxxxx, tg_yyy_xxxxxy, tg_yyy_xxxxxz, tg_yyy_xxxxyy, tg_yyy_xxxxyz, tg_yyy_xxxxzz, tg_yyy_xxxyyy, tg_yyy_xxxyyz, tg_yyy_xxxyzz, tg_yyy_xxxzzz, tg_yyy_xxyyyy, tg_yyy_xxyyyz, tg_yyy_xxyyzz, tg_yyy_xxyzzz, tg_yyy_xxzzzz, tg_yyy_xyyyyy, tg_yyy_xyyyyz, tg_yyy_xyyyzz, tg_yyy_xyyzzz, tg_yyy_xyzzzz, tg_yyy_xzzzzz, tg_yyy_yyyyyy, tg_yyy_yyyyyz, tg_yyy_yyyyzz, tg_yyy_yyyzzz, tg_yyy_yyzzzz, tg_yyy_yzzzzz, tg_yyy_zzzzzz, tg_yyz_xxxxxx, tg_yyz_xxxxxy, tg_yyz_xxxxxz, tg_yyz_xxxxyy, tg_yyz_xxxxyz, tg_yyz_xxxxzz, tg_yyz_xxxyyy, tg_yyz_xxxyyz, tg_yyz_xxxyzz, tg_yyz_xxxzzz, tg_yyz_xxyyyy, tg_yyz_xxyyyz, tg_yyz_xxyyzz, tg_yyz_xxyzzz, tg_yyz_xxzzzz, tg_yyz_xyyyyy, tg_yyz_xyyyyz, tg_yyz_xyyyzz, tg_yyz_xyyzzz, tg_yyz_xyzzzz, tg_yyz_xzzzzz, tg_yyz_yyyyyy, tg_yyz_yyyyyz, tg_yyz_yyyyzz, tg_yyz_yyyzzz, tg_yyz_yyzzzz, tg_yyz_yzzzzz, tg_yyz_zzzzzz, tg_yz_xxxxxz, tg_yz_xxxxyz, tg_yz_xxxxzz, tg_yz_xxxyyz, tg_yz_xxxyz, tg_yz_xxxyzz, tg_yz_xxxzzz, tg_yz_xxyyyz, tg_yz_xxyyz, tg_yz_xxyyzz, tg_yz_xxyzz, tg_yz_xxyzzz, tg_yz_xxzzzz, tg_yz_xyyyyz, tg_yz_xyyyz, tg_yz_xyyyzz, tg_yz_xyyzz, tg_yz_xyyzzz, tg_yz_xyzzz, tg_yz_xyzzzz, tg_yz_xzzzzz, tg_yz_yyyyyy, tg_yz_yyyyyz, tg_yz_yyyyz, tg_yz_yyyyzz, tg_yz_yyyzz, tg_yz_yyyzzz, tg_yz_yyzzz, tg_yz_yyzzzz, tg_yz_yzzzz, tg_yz_yzzzzz, tg_yz_zzzzzz, tg_yzz_xxxxxx, tg_yzz_xxxxxy, tg_yzz_xxxxxz, tg_yzz_xxxxyy, tg_yzz_xxxxyz, tg_yzz_xxxxzz, tg_yzz_xxxyyy, tg_yzz_xxxyyz, tg_yzz_xxxyzz, tg_yzz_xxxzzz, tg_yzz_xxyyyy, tg_yzz_xxyyyz, tg_yzz_xxyyzz, tg_yzz_xxyzzz, tg_yzz_xxzzzz, tg_yzz_xyyyyy, tg_yzz_xyyyyz, tg_yzz_xyyyzz, tg_yzz_xyyzzz, tg_yzz_xyzzzz, tg_yzz_xzzzzz, tg_yzz_yyyyyy, tg_yzz_yyyyyz, tg_yzz_yyyyzz, tg_yzz_yyyzzz, tg_yzz_yyzzzz, tg_yzz_yzzzzz, tg_yzz_zzzzzz, tg_z_xxxxxx, tg_z_xxxxxy, tg_z_xxxxxz, tg_z_xxxxyy, tg_z_xxxxyz, tg_z_xxxxzz, tg_z_xxxyyy, tg_z_xxxyyz, tg_z_xxxyzz, tg_z_xxxzzz, tg_z_xxyyyy, tg_z_xxyyyz, tg_z_xxyyzz, tg_z_xxyzzz, tg_z_xxzzzz, tg_z_xyyyyy, tg_z_xyyyyz, tg_z_xyyyzz, tg_z_xyyzzz, tg_z_xyzzzz, tg_z_xzzzzz, tg_z_yyyyyy, tg_z_yyyyyz, tg_z_yyyyzz, tg_z_yyyzzz, tg_z_yyzzzz, tg_z_yzzzzz, tg_z_zzzzzz, tg_zz_xxxxx, tg_zz_xxxxxx, tg_zz_xxxxxy, tg_zz_xxxxxz, tg_zz_xxxxy, tg_zz_xxxxyy, tg_zz_xxxxyz, tg_zz_xxxxz, tg_zz_xxxxzz, tg_zz_xxxyy, tg_zz_xxxyyy, tg_zz_xxxyyz, tg_zz_xxxyz, tg_zz_xxxyzz, tg_zz_xxxzz, tg_zz_xxxzzz, tg_zz_xxyyy, tg_zz_xxyyyy, tg_zz_xxyyyz, tg_zz_xxyyz, tg_zz_xxyyzz, tg_zz_xxyzz, tg_zz_xxyzzz, tg_zz_xxzzz, tg_zz_xxzzzz, tg_zz_xyyyy, tg_zz_xyyyyy, tg_zz_xyyyyz, tg_zz_xyyyz, tg_zz_xyyyzz, tg_zz_xyyzz, tg_zz_xyyzzz, tg_zz_xyzzz, tg_zz_xyzzzz, tg_zz_xzzzz, tg_zz_xzzzzz, tg_zz_yyyyy, tg_zz_yyyyyy, tg_zz_yyyyyz, tg_zz_yyyyz, tg_zz_yyyyzz, tg_zz_yyyzz, tg_zz_yyyzzz, tg_zz_yyzzz, tg_zz_yyzzzz, tg_zz_yzzzz, tg_zz_yzzzzz, tg_zz_zzzzz, tg_zz_zzzzzz, tg_zzz_xxxxxx, tg_zzz_xxxxxy, tg_zzz_xxxxxz, tg_zzz_xxxxyy, tg_zzz_xxxxyz, tg_zzz_xxxxzz, tg_zzz_xxxyyy, tg_zzz_xxxyyz, tg_zzz_xxxyzz, tg_zzz_xxxzzz, tg_zzz_xxyyyy, tg_zzz_xxyyyz, tg_zzz_xxyyzz, tg_zzz_xxyzzz, tg_zzz_xxzzzz, tg_zzz_xyyyyy, tg_zzz_xyyyyz, tg_zzz_xyyyzz, tg_zzz_xyyzzz, tg_zzz_xyzzzz, tg_zzz_xzzzzz, tg_zzz_yyyyyy, tg_zzz_yyyyyz, tg_zzz_yyyyzz, tg_zzz_yyyzzz, tg_zzz_yyzzzz, tg_zzz_yzzzzz, tg_zzz_zzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_xxxxxx[i] = 2.0 * tg_x_xxxxxx[i] * fxi[i] + 6.0 * tg_xx_xxxxx[i] * fxi[i] + tg_xx_xxxxxx[i] * ra_x[i];

        tg_xxx_xxxxxy[i] = 2.0 * tg_x_xxxxxy[i] * fxi[i] + 5.0 * tg_xx_xxxxy[i] * fxi[i] + tg_xx_xxxxxy[i] * ra_x[i];

        tg_xxx_xxxxxz[i] = 2.0 * tg_x_xxxxxz[i] * fxi[i] + 5.0 * tg_xx_xxxxz[i] * fxi[i] + tg_xx_xxxxxz[i] * ra_x[i];

        tg_xxx_xxxxyy[i] = 2.0 * tg_x_xxxxyy[i] * fxi[i] + 4.0 * tg_xx_xxxyy[i] * fxi[i] + tg_xx_xxxxyy[i] * ra_x[i];

        tg_xxx_xxxxyz[i] = 2.0 * tg_x_xxxxyz[i] * fxi[i] + 4.0 * tg_xx_xxxyz[i] * fxi[i] + tg_xx_xxxxyz[i] * ra_x[i];

        tg_xxx_xxxxzz[i] = 2.0 * tg_x_xxxxzz[i] * fxi[i] + 4.0 * tg_xx_xxxzz[i] * fxi[i] + tg_xx_xxxxzz[i] * ra_x[i];

        tg_xxx_xxxyyy[i] = 2.0 * tg_x_xxxyyy[i] * fxi[i] + 3.0 * tg_xx_xxyyy[i] * fxi[i] + tg_xx_xxxyyy[i] * ra_x[i];

        tg_xxx_xxxyyz[i] = 2.0 * tg_x_xxxyyz[i] * fxi[i] + 3.0 * tg_xx_xxyyz[i] * fxi[i] + tg_xx_xxxyyz[i] * ra_x[i];

        tg_xxx_xxxyzz[i] = 2.0 * tg_x_xxxyzz[i] * fxi[i] + 3.0 * tg_xx_xxyzz[i] * fxi[i] + tg_xx_xxxyzz[i] * ra_x[i];

        tg_xxx_xxxzzz[i] = 2.0 * tg_x_xxxzzz[i] * fxi[i] + 3.0 * tg_xx_xxzzz[i] * fxi[i] + tg_xx_xxxzzz[i] * ra_x[i];

        tg_xxx_xxyyyy[i] = 2.0 * tg_x_xxyyyy[i] * fxi[i] + 2.0 * tg_xx_xyyyy[i] * fxi[i] + tg_xx_xxyyyy[i] * ra_x[i];

        tg_xxx_xxyyyz[i] = 2.0 * tg_x_xxyyyz[i] * fxi[i] + 2.0 * tg_xx_xyyyz[i] * fxi[i] + tg_xx_xxyyyz[i] * ra_x[i];

        tg_xxx_xxyyzz[i] = 2.0 * tg_x_xxyyzz[i] * fxi[i] + 2.0 * tg_xx_xyyzz[i] * fxi[i] + tg_xx_xxyyzz[i] * ra_x[i];

        tg_xxx_xxyzzz[i] = 2.0 * tg_x_xxyzzz[i] * fxi[i] + 2.0 * tg_xx_xyzzz[i] * fxi[i] + tg_xx_xxyzzz[i] * ra_x[i];

        tg_xxx_xxzzzz[i] = 2.0 * tg_x_xxzzzz[i] * fxi[i] + 2.0 * tg_xx_xzzzz[i] * fxi[i] + tg_xx_xxzzzz[i] * ra_x[i];

        tg_xxx_xyyyyy[i] = 2.0 * tg_x_xyyyyy[i] * fxi[i] + tg_xx_yyyyy[i] * fxi[i] + tg_xx_xyyyyy[i] * ra_x[i];

        tg_xxx_xyyyyz[i] = 2.0 * tg_x_xyyyyz[i] * fxi[i] + tg_xx_yyyyz[i] * fxi[i] + tg_xx_xyyyyz[i] * ra_x[i];

        tg_xxx_xyyyzz[i] = 2.0 * tg_x_xyyyzz[i] * fxi[i] + tg_xx_yyyzz[i] * fxi[i] + tg_xx_xyyyzz[i] * ra_x[i];

        tg_xxx_xyyzzz[i] = 2.0 * tg_x_xyyzzz[i] * fxi[i] + tg_xx_yyzzz[i] * fxi[i] + tg_xx_xyyzzz[i] * ra_x[i];

        tg_xxx_xyzzzz[i] = 2.0 * tg_x_xyzzzz[i] * fxi[i] + tg_xx_yzzzz[i] * fxi[i] + tg_xx_xyzzzz[i] * ra_x[i];

        tg_xxx_xzzzzz[i] = 2.0 * tg_x_xzzzzz[i] * fxi[i] + tg_xx_zzzzz[i] * fxi[i] + tg_xx_xzzzzz[i] * ra_x[i];

        tg_xxx_yyyyyy[i] = 2.0 * tg_x_yyyyyy[i] * fxi[i] + tg_xx_yyyyyy[i] * ra_x[i];

        tg_xxx_yyyyyz[i] = 2.0 * tg_x_yyyyyz[i] * fxi[i] + tg_xx_yyyyyz[i] * ra_x[i];

        tg_xxx_yyyyzz[i] = 2.0 * tg_x_yyyyzz[i] * fxi[i] + tg_xx_yyyyzz[i] * ra_x[i];

        tg_xxx_yyyzzz[i] = 2.0 * tg_x_yyyzzz[i] * fxi[i] + tg_xx_yyyzzz[i] * ra_x[i];

        tg_xxx_yyzzzz[i] = 2.0 * tg_x_yyzzzz[i] * fxi[i] + tg_xx_yyzzzz[i] * ra_x[i];

        tg_xxx_yzzzzz[i] = 2.0 * tg_x_yzzzzz[i] * fxi[i] + tg_xx_yzzzzz[i] * ra_x[i];

        tg_xxx_zzzzzz[i] = 2.0 * tg_x_zzzzzz[i] * fxi[i] + tg_xx_zzzzzz[i] * ra_x[i];

        tg_xxy_xxxxxx[i] = tg_xx_xxxxxx[i] * ra_y[i];

        tg_xxy_xxxxxy[i] = tg_xx_xxxxx[i] * fxi[i] + tg_xx_xxxxxy[i] * ra_y[i];

        tg_xxy_xxxxxz[i] = tg_xx_xxxxxz[i] * ra_y[i];

        tg_xxy_xxxxyy[i] = 2.0 * tg_xx_xxxxy[i] * fxi[i] + tg_xx_xxxxyy[i] * ra_y[i];

        tg_xxy_xxxxyz[i] = tg_xx_xxxxz[i] * fxi[i] + tg_xx_xxxxyz[i] * ra_y[i];

        tg_xxy_xxxxzz[i] = tg_xx_xxxxzz[i] * ra_y[i];

        tg_xxy_xxxyyy[i] = 3.0 * tg_xx_xxxyy[i] * fxi[i] + tg_xx_xxxyyy[i] * ra_y[i];

        tg_xxy_xxxyyz[i] = 2.0 * tg_xx_xxxyz[i] * fxi[i] + tg_xx_xxxyyz[i] * ra_y[i];

        tg_xxy_xxxyzz[i] = tg_xx_xxxzz[i] * fxi[i] + tg_xx_xxxyzz[i] * ra_y[i];

        tg_xxy_xxxzzz[i] = tg_xx_xxxzzz[i] * ra_y[i];

        tg_xxy_xxyyyy[i] = 4.0 * tg_xx_xxyyy[i] * fxi[i] + tg_xx_xxyyyy[i] * ra_y[i];

        tg_xxy_xxyyyz[i] = 3.0 * tg_xx_xxyyz[i] * fxi[i] + tg_xx_xxyyyz[i] * ra_y[i];

        tg_xxy_xxyyzz[i] = 2.0 * tg_xx_xxyzz[i] * fxi[i] + tg_xx_xxyyzz[i] * ra_y[i];

        tg_xxy_xxyzzz[i] = tg_xx_xxzzz[i] * fxi[i] + tg_xx_xxyzzz[i] * ra_y[i];

        tg_xxy_xxzzzz[i] = tg_xx_xxzzzz[i] * ra_y[i];

        tg_xxy_xyyyyy[i] = 5.0 * tg_xx_xyyyy[i] * fxi[i] + tg_xx_xyyyyy[i] * ra_y[i];

        tg_xxy_xyyyyz[i] = 4.0 * tg_xx_xyyyz[i] * fxi[i] + tg_xx_xyyyyz[i] * ra_y[i];

        tg_xxy_xyyyzz[i] = 3.0 * tg_xx_xyyzz[i] * fxi[i] + tg_xx_xyyyzz[i] * ra_y[i];

        tg_xxy_xyyzzz[i] = 2.0 * tg_xx_xyzzz[i] * fxi[i] + tg_xx_xyyzzz[i] * ra_y[i];

        tg_xxy_xyzzzz[i] = tg_xx_xzzzz[i] * fxi[i] + tg_xx_xyzzzz[i] * ra_y[i];

        tg_xxy_xzzzzz[i] = tg_xx_xzzzzz[i] * ra_y[i];

        tg_xxy_yyyyyy[i] = tg_y_yyyyyy[i] * fxi[i] + tg_xy_yyyyyy[i] * ra_x[i];

        tg_xxy_yyyyyz[i] = tg_y_yyyyyz[i] * fxi[i] + tg_xy_yyyyyz[i] * ra_x[i];

        tg_xxy_yyyyzz[i] = tg_y_yyyyzz[i] * fxi[i] + tg_xy_yyyyzz[i] * ra_x[i];

        tg_xxy_yyyzzz[i] = tg_y_yyyzzz[i] * fxi[i] + tg_xy_yyyzzz[i] * ra_x[i];

        tg_xxy_yyzzzz[i] = tg_y_yyzzzz[i] * fxi[i] + tg_xy_yyzzzz[i] * ra_x[i];

        tg_xxy_yzzzzz[i] = tg_y_yzzzzz[i] * fxi[i] + tg_xy_yzzzzz[i] * ra_x[i];

        tg_xxy_zzzzzz[i] = tg_xx_zzzzzz[i] * ra_y[i];

        tg_xxz_xxxxxx[i] = tg_xx_xxxxxx[i] * ra_z[i];

        tg_xxz_xxxxxy[i] = tg_xx_xxxxxy[i] * ra_z[i];

        tg_xxz_xxxxxz[i] = tg_xx_xxxxx[i] * fxi[i] + tg_xx_xxxxxz[i] * ra_z[i];

        tg_xxz_xxxxyy[i] = tg_xx_xxxxyy[i] * ra_z[i];

        tg_xxz_xxxxyz[i] = tg_xx_xxxxy[i] * fxi[i] + tg_xx_xxxxyz[i] * ra_z[i];

        tg_xxz_xxxxzz[i] = 2.0 * tg_xx_xxxxz[i] * fxi[i] + tg_xx_xxxxzz[i] * ra_z[i];

        tg_xxz_xxxyyy[i] = tg_xx_xxxyyy[i] * ra_z[i];

        tg_xxz_xxxyyz[i] = tg_xx_xxxyy[i] * fxi[i] + tg_xx_xxxyyz[i] * ra_z[i];

        tg_xxz_xxxyzz[i] = 2.0 * tg_xx_xxxyz[i] * fxi[i] + tg_xx_xxxyzz[i] * ra_z[i];

        tg_xxz_xxxzzz[i] = 3.0 * tg_xx_xxxzz[i] * fxi[i] + tg_xx_xxxzzz[i] * ra_z[i];

        tg_xxz_xxyyyy[i] = tg_xx_xxyyyy[i] * ra_z[i];

        tg_xxz_xxyyyz[i] = tg_xx_xxyyy[i] * fxi[i] + tg_xx_xxyyyz[i] * ra_z[i];

        tg_xxz_xxyyzz[i] = 2.0 * tg_xx_xxyyz[i] * fxi[i] + tg_xx_xxyyzz[i] * ra_z[i];

        tg_xxz_xxyzzz[i] = 3.0 * tg_xx_xxyzz[i] * fxi[i] + tg_xx_xxyzzz[i] * ra_z[i];

        tg_xxz_xxzzzz[i] = 4.0 * tg_xx_xxzzz[i] * fxi[i] + tg_xx_xxzzzz[i] * ra_z[i];

        tg_xxz_xyyyyy[i] = tg_xx_xyyyyy[i] * ra_z[i];

        tg_xxz_xyyyyz[i] = tg_xx_xyyyy[i] * fxi[i] + tg_xx_xyyyyz[i] * ra_z[i];

        tg_xxz_xyyyzz[i] = 2.0 * tg_xx_xyyyz[i] * fxi[i] + tg_xx_xyyyzz[i] * ra_z[i];

        tg_xxz_xyyzzz[i] = 3.0 * tg_xx_xyyzz[i] * fxi[i] + tg_xx_xyyzzz[i] * ra_z[i];

        tg_xxz_xyzzzz[i] = 4.0 * tg_xx_xyzzz[i] * fxi[i] + tg_xx_xyzzzz[i] * ra_z[i];

        tg_xxz_xzzzzz[i] = 5.0 * tg_xx_xzzzz[i] * fxi[i] + tg_xx_xzzzzz[i] * ra_z[i];

        tg_xxz_yyyyyy[i] = tg_xx_yyyyyy[i] * ra_z[i];

        tg_xxz_yyyyyz[i] = tg_z_yyyyyz[i] * fxi[i] + tg_xz_yyyyyz[i] * ra_x[i];

        tg_xxz_yyyyzz[i] = tg_z_yyyyzz[i] * fxi[i] + tg_xz_yyyyzz[i] * ra_x[i];

        tg_xxz_yyyzzz[i] = tg_z_yyyzzz[i] * fxi[i] + tg_xz_yyyzzz[i] * ra_x[i];

        tg_xxz_yyzzzz[i] = tg_z_yyzzzz[i] * fxi[i] + tg_xz_yyzzzz[i] * ra_x[i];

        tg_xxz_yzzzzz[i] = tg_z_yzzzzz[i] * fxi[i] + tg_xz_yzzzzz[i] * ra_x[i];

        tg_xxz_zzzzzz[i] = tg_z_zzzzzz[i] * fxi[i] + tg_xz_zzzzzz[i] * ra_x[i];

        tg_xyy_xxxxxx[i] = 6.0 * tg_yy_xxxxx[i] * fxi[i] + tg_yy_xxxxxx[i] * ra_x[i];

        tg_xyy_xxxxxy[i] = 5.0 * tg_yy_xxxxy[i] * fxi[i] + tg_yy_xxxxxy[i] * ra_x[i];

        tg_xyy_xxxxxz[i] = 5.0 * tg_yy_xxxxz[i] * fxi[i] + tg_yy_xxxxxz[i] * ra_x[i];

        tg_xyy_xxxxyy[i] = 4.0 * tg_yy_xxxyy[i] * fxi[i] + tg_yy_xxxxyy[i] * ra_x[i];

        tg_xyy_xxxxyz[i] = 4.0 * tg_yy_xxxyz[i] * fxi[i] + tg_yy_xxxxyz[i] * ra_x[i];

        tg_xyy_xxxxzz[i] = 4.0 * tg_yy_xxxzz[i] * fxi[i] + tg_yy_xxxxzz[i] * ra_x[i];

        tg_xyy_xxxyyy[i] = 3.0 * tg_yy_xxyyy[i] * fxi[i] + tg_yy_xxxyyy[i] * ra_x[i];

        tg_xyy_xxxyyz[i] = 3.0 * tg_yy_xxyyz[i] * fxi[i] + tg_yy_xxxyyz[i] * ra_x[i];

        tg_xyy_xxxyzz[i] = 3.0 * tg_yy_xxyzz[i] * fxi[i] + tg_yy_xxxyzz[i] * ra_x[i];

        tg_xyy_xxxzzz[i] = 3.0 * tg_yy_xxzzz[i] * fxi[i] + tg_yy_xxxzzz[i] * ra_x[i];

        tg_xyy_xxyyyy[i] = 2.0 * tg_yy_xyyyy[i] * fxi[i] + tg_yy_xxyyyy[i] * ra_x[i];

        tg_xyy_xxyyyz[i] = 2.0 * tg_yy_xyyyz[i] * fxi[i] + tg_yy_xxyyyz[i] * ra_x[i];

        tg_xyy_xxyyzz[i] = 2.0 * tg_yy_xyyzz[i] * fxi[i] + tg_yy_xxyyzz[i] * ra_x[i];

        tg_xyy_xxyzzz[i] = 2.0 * tg_yy_xyzzz[i] * fxi[i] + tg_yy_xxyzzz[i] * ra_x[i];

        tg_xyy_xxzzzz[i] = 2.0 * tg_yy_xzzzz[i] * fxi[i] + tg_yy_xxzzzz[i] * ra_x[i];

        tg_xyy_xyyyyy[i] = tg_yy_yyyyy[i] * fxi[i] + tg_yy_xyyyyy[i] * ra_x[i];

        tg_xyy_xyyyyz[i] = tg_yy_yyyyz[i] * fxi[i] + tg_yy_xyyyyz[i] * ra_x[i];

        tg_xyy_xyyyzz[i] = tg_yy_yyyzz[i] * fxi[i] + tg_yy_xyyyzz[i] * ra_x[i];

        tg_xyy_xyyzzz[i] = tg_yy_yyzzz[i] * fxi[i] + tg_yy_xyyzzz[i] * ra_x[i];

        tg_xyy_xyzzzz[i] = tg_yy_yzzzz[i] * fxi[i] + tg_yy_xyzzzz[i] * ra_x[i];

        tg_xyy_xzzzzz[i] = tg_yy_zzzzz[i] * fxi[i] + tg_yy_xzzzzz[i] * ra_x[i];

        tg_xyy_yyyyyy[i] = tg_yy_yyyyyy[i] * ra_x[i];

        tg_xyy_yyyyyz[i] = tg_yy_yyyyyz[i] * ra_x[i];

        tg_xyy_yyyyzz[i] = tg_yy_yyyyzz[i] * ra_x[i];

        tg_xyy_yyyzzz[i] = tg_yy_yyyzzz[i] * ra_x[i];

        tg_xyy_yyzzzz[i] = tg_yy_yyzzzz[i] * ra_x[i];

        tg_xyy_yzzzzz[i] = tg_yy_yzzzzz[i] * ra_x[i];

        tg_xyy_zzzzzz[i] = tg_yy_zzzzzz[i] * ra_x[i];

        tg_xyz_xxxxxx[i] = tg_xz_xxxxxx[i] * ra_y[i];

        tg_xyz_xxxxxy[i] = tg_xy_xxxxxy[i] * ra_z[i];

        tg_xyz_xxxxxz[i] = tg_xz_xxxxxz[i] * ra_y[i];

        tg_xyz_xxxxyy[i] = tg_xy_xxxxyy[i] * ra_z[i];

        tg_xyz_xxxxyz[i] = 4.0 * tg_yz_xxxyz[i] * fxi[i] + tg_yz_xxxxyz[i] * ra_x[i];

        tg_xyz_xxxxzz[i] = tg_xz_xxxxzz[i] * ra_y[i];

        tg_xyz_xxxyyy[i] = tg_xy_xxxyyy[i] * ra_z[i];

        tg_xyz_xxxyyz[i] = 3.0 * tg_yz_xxyyz[i] * fxi[i] + tg_yz_xxxyyz[i] * ra_x[i];

        tg_xyz_xxxyzz[i] = 3.0 * tg_yz_xxyzz[i] * fxi[i] + tg_yz_xxxyzz[i] * ra_x[i];

        tg_xyz_xxxzzz[i] = tg_xz_xxxzzz[i] * ra_y[i];

        tg_xyz_xxyyyy[i] = tg_xy_xxyyyy[i] * ra_z[i];

        tg_xyz_xxyyyz[i] = 2.0 * tg_yz_xyyyz[i] * fxi[i] + tg_yz_xxyyyz[i] * ra_x[i];

        tg_xyz_xxyyzz[i] = 2.0 * tg_yz_xyyzz[i] * fxi[i] + tg_yz_xxyyzz[i] * ra_x[i];

        tg_xyz_xxyzzz[i] = 2.0 * tg_yz_xyzzz[i] * fxi[i] + tg_yz_xxyzzz[i] * ra_x[i];

        tg_xyz_xxzzzz[i] = tg_xz_xxzzzz[i] * ra_y[i];

        tg_xyz_xyyyyy[i] = tg_xy_xyyyyy[i] * ra_z[i];

        tg_xyz_xyyyyz[i] = tg_yz_yyyyz[i] * fxi[i] + tg_yz_xyyyyz[i] * ra_x[i];

        tg_xyz_xyyyzz[i] = tg_yz_yyyzz[i] * fxi[i] + tg_yz_xyyyzz[i] * ra_x[i];

        tg_xyz_xyyzzz[i] = tg_yz_yyzzz[i] * fxi[i] + tg_yz_xyyzzz[i] * ra_x[i];

        tg_xyz_xyzzzz[i] = tg_yz_yzzzz[i] * fxi[i] + tg_yz_xyzzzz[i] * ra_x[i];

        tg_xyz_xzzzzz[i] = tg_xz_xzzzzz[i] * ra_y[i];

        tg_xyz_yyyyyy[i] = tg_yz_yyyyyy[i] * ra_x[i];

        tg_xyz_yyyyyz[i] = tg_yz_yyyyyz[i] * ra_x[i];

        tg_xyz_yyyyzz[i] = tg_yz_yyyyzz[i] * ra_x[i];

        tg_xyz_yyyzzz[i] = tg_yz_yyyzzz[i] * ra_x[i];

        tg_xyz_yyzzzz[i] = tg_yz_yyzzzz[i] * ra_x[i];

        tg_xyz_yzzzzz[i] = tg_yz_yzzzzz[i] * ra_x[i];

        tg_xyz_zzzzzz[i] = tg_yz_zzzzzz[i] * ra_x[i];

        tg_xzz_xxxxxx[i] = 6.0 * tg_zz_xxxxx[i] * fxi[i] + tg_zz_xxxxxx[i] * ra_x[i];

        tg_xzz_xxxxxy[i] = 5.0 * tg_zz_xxxxy[i] * fxi[i] + tg_zz_xxxxxy[i] * ra_x[i];

        tg_xzz_xxxxxz[i] = 5.0 * tg_zz_xxxxz[i] * fxi[i] + tg_zz_xxxxxz[i] * ra_x[i];

        tg_xzz_xxxxyy[i] = 4.0 * tg_zz_xxxyy[i] * fxi[i] + tg_zz_xxxxyy[i] * ra_x[i];

        tg_xzz_xxxxyz[i] = 4.0 * tg_zz_xxxyz[i] * fxi[i] + tg_zz_xxxxyz[i] * ra_x[i];

        tg_xzz_xxxxzz[i] = 4.0 * tg_zz_xxxzz[i] * fxi[i] + tg_zz_xxxxzz[i] * ra_x[i];

        tg_xzz_xxxyyy[i] = 3.0 * tg_zz_xxyyy[i] * fxi[i] + tg_zz_xxxyyy[i] * ra_x[i];

        tg_xzz_xxxyyz[i] = 3.0 * tg_zz_xxyyz[i] * fxi[i] + tg_zz_xxxyyz[i] * ra_x[i];

        tg_xzz_xxxyzz[i] = 3.0 * tg_zz_xxyzz[i] * fxi[i] + tg_zz_xxxyzz[i] * ra_x[i];

        tg_xzz_xxxzzz[i] = 3.0 * tg_zz_xxzzz[i] * fxi[i] + tg_zz_xxxzzz[i] * ra_x[i];

        tg_xzz_xxyyyy[i] = 2.0 * tg_zz_xyyyy[i] * fxi[i] + tg_zz_xxyyyy[i] * ra_x[i];

        tg_xzz_xxyyyz[i] = 2.0 * tg_zz_xyyyz[i] * fxi[i] + tg_zz_xxyyyz[i] * ra_x[i];

        tg_xzz_xxyyzz[i] = 2.0 * tg_zz_xyyzz[i] * fxi[i] + tg_zz_xxyyzz[i] * ra_x[i];

        tg_xzz_xxyzzz[i] = 2.0 * tg_zz_xyzzz[i] * fxi[i] + tg_zz_xxyzzz[i] * ra_x[i];

        tg_xzz_xxzzzz[i] = 2.0 * tg_zz_xzzzz[i] * fxi[i] + tg_zz_xxzzzz[i] * ra_x[i];

        tg_xzz_xyyyyy[i] = tg_zz_yyyyy[i] * fxi[i] + tg_zz_xyyyyy[i] * ra_x[i];

        tg_xzz_xyyyyz[i] = tg_zz_yyyyz[i] * fxi[i] + tg_zz_xyyyyz[i] * ra_x[i];

        tg_xzz_xyyyzz[i] = tg_zz_yyyzz[i] * fxi[i] + tg_zz_xyyyzz[i] * ra_x[i];

        tg_xzz_xyyzzz[i] = tg_zz_yyzzz[i] * fxi[i] + tg_zz_xyyzzz[i] * ra_x[i];

        tg_xzz_xyzzzz[i] = tg_zz_yzzzz[i] * fxi[i] + tg_zz_xyzzzz[i] * ra_x[i];

        tg_xzz_xzzzzz[i] = tg_zz_zzzzz[i] * fxi[i] + tg_zz_xzzzzz[i] * ra_x[i];

        tg_xzz_yyyyyy[i] = tg_zz_yyyyyy[i] * ra_x[i];

        tg_xzz_yyyyyz[i] = tg_zz_yyyyyz[i] * ra_x[i];

        tg_xzz_yyyyzz[i] = tg_zz_yyyyzz[i] * ra_x[i];

        tg_xzz_yyyzzz[i] = tg_zz_yyyzzz[i] * ra_x[i];

        tg_xzz_yyzzzz[i] = tg_zz_yyzzzz[i] * ra_x[i];

        tg_xzz_yzzzzz[i] = tg_zz_yzzzzz[i] * ra_x[i];

        tg_xzz_zzzzzz[i] = tg_zz_zzzzzz[i] * ra_x[i];

        tg_yyy_xxxxxx[i] = 2.0 * tg_y_xxxxxx[i] * fxi[i] + tg_yy_xxxxxx[i] * ra_y[i];

        tg_yyy_xxxxxy[i] = 2.0 * tg_y_xxxxxy[i] * fxi[i] + tg_yy_xxxxx[i] * fxi[i] + tg_yy_xxxxxy[i] * ra_y[i];

        tg_yyy_xxxxxz[i] = 2.0 * tg_y_xxxxxz[i] * fxi[i] + tg_yy_xxxxxz[i] * ra_y[i];

        tg_yyy_xxxxyy[i] = 2.0 * tg_y_xxxxyy[i] * fxi[i] + 2.0 * tg_yy_xxxxy[i] * fxi[i] + tg_yy_xxxxyy[i] * ra_y[i];

        tg_yyy_xxxxyz[i] = 2.0 * tg_y_xxxxyz[i] * fxi[i] + tg_yy_xxxxz[i] * fxi[i] + tg_yy_xxxxyz[i] * ra_y[i];

        tg_yyy_xxxxzz[i] = 2.0 * tg_y_xxxxzz[i] * fxi[i] + tg_yy_xxxxzz[i] * ra_y[i];

        tg_yyy_xxxyyy[i] = 2.0 * tg_y_xxxyyy[i] * fxi[i] + 3.0 * tg_yy_xxxyy[i] * fxi[i] + tg_yy_xxxyyy[i] * ra_y[i];

        tg_yyy_xxxyyz[i] = 2.0 * tg_y_xxxyyz[i] * fxi[i] + 2.0 * tg_yy_xxxyz[i] * fxi[i] + tg_yy_xxxyyz[i] * ra_y[i];

        tg_yyy_xxxyzz[i] = 2.0 * tg_y_xxxyzz[i] * fxi[i] + tg_yy_xxxzz[i] * fxi[i] + tg_yy_xxxyzz[i] * ra_y[i];

        tg_yyy_xxxzzz[i] = 2.0 * tg_y_xxxzzz[i] * fxi[i] + tg_yy_xxxzzz[i] * ra_y[i];

        tg_yyy_xxyyyy[i] = 2.0 * tg_y_xxyyyy[i] * fxi[i] + 4.0 * tg_yy_xxyyy[i] * fxi[i] + tg_yy_xxyyyy[i] * ra_y[i];

        tg_yyy_xxyyyz[i] = 2.0 * tg_y_xxyyyz[i] * fxi[i] + 3.0 * tg_yy_xxyyz[i] * fxi[i] + tg_yy_xxyyyz[i] * ra_y[i];

        tg_yyy_xxyyzz[i] = 2.0 * tg_y_xxyyzz[i] * fxi[i] + 2.0 * tg_yy_xxyzz[i] * fxi[i] + tg_yy_xxyyzz[i] * ra_y[i];

        tg_yyy_xxyzzz[i] = 2.0 * tg_y_xxyzzz[i] * fxi[i] + tg_yy_xxzzz[i] * fxi[i] + tg_yy_xxyzzz[i] * ra_y[i];

        tg_yyy_xxzzzz[i] = 2.0 * tg_y_xxzzzz[i] * fxi[i] + tg_yy_xxzzzz[i] * ra_y[i];

        tg_yyy_xyyyyy[i] = 2.0 * tg_y_xyyyyy[i] * fxi[i] + 5.0 * tg_yy_xyyyy[i] * fxi[i] + tg_yy_xyyyyy[i] * ra_y[i];

        tg_yyy_xyyyyz[i] = 2.0 * tg_y_xyyyyz[i] * fxi[i] + 4.0 * tg_yy_xyyyz[i] * fxi[i] + tg_yy_xyyyyz[i] * ra_y[i];

        tg_yyy_xyyyzz[i] = 2.0 * tg_y_xyyyzz[i] * fxi[i] + 3.0 * tg_yy_xyyzz[i] * fxi[i] + tg_yy_xyyyzz[i] * ra_y[i];

        tg_yyy_xyyzzz[i] = 2.0 * tg_y_xyyzzz[i] * fxi[i] + 2.0 * tg_yy_xyzzz[i] * fxi[i] + tg_yy_xyyzzz[i] * ra_y[i];

        tg_yyy_xyzzzz[i] = 2.0 * tg_y_xyzzzz[i] * fxi[i] + tg_yy_xzzzz[i] * fxi[i] + tg_yy_xyzzzz[i] * ra_y[i];

        tg_yyy_xzzzzz[i] = 2.0 * tg_y_xzzzzz[i] * fxi[i] + tg_yy_xzzzzz[i] * ra_y[i];

        tg_yyy_yyyyyy[i] = 2.0 * tg_y_yyyyyy[i] * fxi[i] + 6.0 * tg_yy_yyyyy[i] * fxi[i] + tg_yy_yyyyyy[i] * ra_y[i];

        tg_yyy_yyyyyz[i] = 2.0 * tg_y_yyyyyz[i] * fxi[i] + 5.0 * tg_yy_yyyyz[i] * fxi[i] + tg_yy_yyyyyz[i] * ra_y[i];

        tg_yyy_yyyyzz[i] = 2.0 * tg_y_yyyyzz[i] * fxi[i] + 4.0 * tg_yy_yyyzz[i] * fxi[i] + tg_yy_yyyyzz[i] * ra_y[i];

        tg_yyy_yyyzzz[i] = 2.0 * tg_y_yyyzzz[i] * fxi[i] + 3.0 * tg_yy_yyzzz[i] * fxi[i] + tg_yy_yyyzzz[i] * ra_y[i];

        tg_yyy_yyzzzz[i] = 2.0 * tg_y_yyzzzz[i] * fxi[i] + 2.0 * tg_yy_yzzzz[i] * fxi[i] + tg_yy_yyzzzz[i] * ra_y[i];

        tg_yyy_yzzzzz[i] = 2.0 * tg_y_yzzzzz[i] * fxi[i] + tg_yy_zzzzz[i] * fxi[i] + tg_yy_yzzzzz[i] * ra_y[i];

        tg_yyy_zzzzzz[i] = 2.0 * tg_y_zzzzzz[i] * fxi[i] + tg_yy_zzzzzz[i] * ra_y[i];

        tg_yyz_xxxxxx[i] = tg_yy_xxxxxx[i] * ra_z[i];

        tg_yyz_xxxxxy[i] = tg_yy_xxxxxy[i] * ra_z[i];

        tg_yyz_xxxxxz[i] = tg_z_xxxxxz[i] * fxi[i] + tg_yz_xxxxxz[i] * ra_y[i];

        tg_yyz_xxxxyy[i] = tg_yy_xxxxyy[i] * ra_z[i];

        tg_yyz_xxxxyz[i] = tg_yy_xxxxy[i] * fxi[i] + tg_yy_xxxxyz[i] * ra_z[i];

        tg_yyz_xxxxzz[i] = tg_z_xxxxzz[i] * fxi[i] + tg_yz_xxxxzz[i] * ra_y[i];

        tg_yyz_xxxyyy[i] = tg_yy_xxxyyy[i] * ra_z[i];

        tg_yyz_xxxyyz[i] = tg_yy_xxxyy[i] * fxi[i] + tg_yy_xxxyyz[i] * ra_z[i];

        tg_yyz_xxxyzz[i] = 2.0 * tg_yy_xxxyz[i] * fxi[i] + tg_yy_xxxyzz[i] * ra_z[i];

        tg_yyz_xxxzzz[i] = tg_z_xxxzzz[i] * fxi[i] + tg_yz_xxxzzz[i] * ra_y[i];

        tg_yyz_xxyyyy[i] = tg_yy_xxyyyy[i] * ra_z[i];

        tg_yyz_xxyyyz[i] = tg_yy_xxyyy[i] * fxi[i] + tg_yy_xxyyyz[i] * ra_z[i];

        tg_yyz_xxyyzz[i] = 2.0 * tg_yy_xxyyz[i] * fxi[i] + tg_yy_xxyyzz[i] * ra_z[i];

        tg_yyz_xxyzzz[i] = 3.0 * tg_yy_xxyzz[i] * fxi[i] + tg_yy_xxyzzz[i] * ra_z[i];

        tg_yyz_xxzzzz[i] = tg_z_xxzzzz[i] * fxi[i] + tg_yz_xxzzzz[i] * ra_y[i];

        tg_yyz_xyyyyy[i] = tg_yy_xyyyyy[i] * ra_z[i];

        tg_yyz_xyyyyz[i] = tg_yy_xyyyy[i] * fxi[i] + tg_yy_xyyyyz[i] * ra_z[i];

        tg_yyz_xyyyzz[i] = 2.0 * tg_yy_xyyyz[i] * fxi[i] + tg_yy_xyyyzz[i] * ra_z[i];

        tg_yyz_xyyzzz[i] = 3.0 * tg_yy_xyyzz[i] * fxi[i] + tg_yy_xyyzzz[i] * ra_z[i];

        tg_yyz_xyzzzz[i] = 4.0 * tg_yy_xyzzz[i] * fxi[i] + tg_yy_xyzzzz[i] * ra_z[i];

        tg_yyz_xzzzzz[i] = tg_z_xzzzzz[i] * fxi[i] + tg_yz_xzzzzz[i] * ra_y[i];

        tg_yyz_yyyyyy[i] = tg_yy_yyyyyy[i] * ra_z[i];

        tg_yyz_yyyyyz[i] = tg_yy_yyyyy[i] * fxi[i] + tg_yy_yyyyyz[i] * ra_z[i];

        tg_yyz_yyyyzz[i] = 2.0 * tg_yy_yyyyz[i] * fxi[i] + tg_yy_yyyyzz[i] * ra_z[i];

        tg_yyz_yyyzzz[i] = 3.0 * tg_yy_yyyzz[i] * fxi[i] + tg_yy_yyyzzz[i] * ra_z[i];

        tg_yyz_yyzzzz[i] = 4.0 * tg_yy_yyzzz[i] * fxi[i] + tg_yy_yyzzzz[i] * ra_z[i];

        tg_yyz_yzzzzz[i] = 5.0 * tg_yy_yzzzz[i] * fxi[i] + tg_yy_yzzzzz[i] * ra_z[i];

        tg_yyz_zzzzzz[i] = tg_z_zzzzzz[i] * fxi[i] + tg_yz_zzzzzz[i] * ra_y[i];

        tg_yzz_xxxxxx[i] = tg_zz_xxxxxx[i] * ra_y[i];

        tg_yzz_xxxxxy[i] = tg_zz_xxxxx[i] * fxi[i] + tg_zz_xxxxxy[i] * ra_y[i];

        tg_yzz_xxxxxz[i] = tg_zz_xxxxxz[i] * ra_y[i];

        tg_yzz_xxxxyy[i] = 2.0 * tg_zz_xxxxy[i] * fxi[i] + tg_zz_xxxxyy[i] * ra_y[i];

        tg_yzz_xxxxyz[i] = tg_zz_xxxxz[i] * fxi[i] + tg_zz_xxxxyz[i] * ra_y[i];

        tg_yzz_xxxxzz[i] = tg_zz_xxxxzz[i] * ra_y[i];

        tg_yzz_xxxyyy[i] = 3.0 * tg_zz_xxxyy[i] * fxi[i] + tg_zz_xxxyyy[i] * ra_y[i];

        tg_yzz_xxxyyz[i] = 2.0 * tg_zz_xxxyz[i] * fxi[i] + tg_zz_xxxyyz[i] * ra_y[i];

        tg_yzz_xxxyzz[i] = tg_zz_xxxzz[i] * fxi[i] + tg_zz_xxxyzz[i] * ra_y[i];

        tg_yzz_xxxzzz[i] = tg_zz_xxxzzz[i] * ra_y[i];

        tg_yzz_xxyyyy[i] = 4.0 * tg_zz_xxyyy[i] * fxi[i] + tg_zz_xxyyyy[i] * ra_y[i];

        tg_yzz_xxyyyz[i] = 3.0 * tg_zz_xxyyz[i] * fxi[i] + tg_zz_xxyyyz[i] * ra_y[i];

        tg_yzz_xxyyzz[i] = 2.0 * tg_zz_xxyzz[i] * fxi[i] + tg_zz_xxyyzz[i] * ra_y[i];

        tg_yzz_xxyzzz[i] = tg_zz_xxzzz[i] * fxi[i] + tg_zz_xxyzzz[i] * ra_y[i];

        tg_yzz_xxzzzz[i] = tg_zz_xxzzzz[i] * ra_y[i];

        tg_yzz_xyyyyy[i] = 5.0 * tg_zz_xyyyy[i] * fxi[i] + tg_zz_xyyyyy[i] * ra_y[i];

        tg_yzz_xyyyyz[i] = 4.0 * tg_zz_xyyyz[i] * fxi[i] + tg_zz_xyyyyz[i] * ra_y[i];

        tg_yzz_xyyyzz[i] = 3.0 * tg_zz_xyyzz[i] * fxi[i] + tg_zz_xyyyzz[i] * ra_y[i];

        tg_yzz_xyyzzz[i] = 2.0 * tg_zz_xyzzz[i] * fxi[i] + tg_zz_xyyzzz[i] * ra_y[i];

        tg_yzz_xyzzzz[i] = tg_zz_xzzzz[i] * fxi[i] + tg_zz_xyzzzz[i] * ra_y[i];

        tg_yzz_xzzzzz[i] = tg_zz_xzzzzz[i] * ra_y[i];

        tg_yzz_yyyyyy[i] = 6.0 * tg_zz_yyyyy[i] * fxi[i] + tg_zz_yyyyyy[i] * ra_y[i];

        tg_yzz_yyyyyz[i] = 5.0 * tg_zz_yyyyz[i] * fxi[i] + tg_zz_yyyyyz[i] * ra_y[i];

        tg_yzz_yyyyzz[i] = 4.0 * tg_zz_yyyzz[i] * fxi[i] + tg_zz_yyyyzz[i] * ra_y[i];

        tg_yzz_yyyzzz[i] = 3.0 * tg_zz_yyzzz[i] * fxi[i] + tg_zz_yyyzzz[i] * ra_y[i];

        tg_yzz_yyzzzz[i] = 2.0 * tg_zz_yzzzz[i] * fxi[i] + tg_zz_yyzzzz[i] * ra_y[i];

        tg_yzz_yzzzzz[i] = tg_zz_zzzzz[i] * fxi[i] + tg_zz_yzzzzz[i] * ra_y[i];

        tg_yzz_zzzzzz[i] = tg_zz_zzzzzz[i] * ra_y[i];

        tg_zzz_xxxxxx[i] = 2.0 * tg_z_xxxxxx[i] * fxi[i] + tg_zz_xxxxxx[i] * ra_z[i];

        tg_zzz_xxxxxy[i] = 2.0 * tg_z_xxxxxy[i] * fxi[i] + tg_zz_xxxxxy[i] * ra_z[i];

        tg_zzz_xxxxxz[i] = 2.0 * tg_z_xxxxxz[i] * fxi[i] + tg_zz_xxxxx[i] * fxi[i] + tg_zz_xxxxxz[i] * ra_z[i];

        tg_zzz_xxxxyy[i] = 2.0 * tg_z_xxxxyy[i] * fxi[i] + tg_zz_xxxxyy[i] * ra_z[i];

        tg_zzz_xxxxyz[i] = 2.0 * tg_z_xxxxyz[i] * fxi[i] + tg_zz_xxxxy[i] * fxi[i] + tg_zz_xxxxyz[i] * ra_z[i];

        tg_zzz_xxxxzz[i] = 2.0 * tg_z_xxxxzz[i] * fxi[i] + 2.0 * tg_zz_xxxxz[i] * fxi[i] + tg_zz_xxxxzz[i] * ra_z[i];

        tg_zzz_xxxyyy[i] = 2.0 * tg_z_xxxyyy[i] * fxi[i] + tg_zz_xxxyyy[i] * ra_z[i];

        tg_zzz_xxxyyz[i] = 2.0 * tg_z_xxxyyz[i] * fxi[i] + tg_zz_xxxyy[i] * fxi[i] + tg_zz_xxxyyz[i] * ra_z[i];

        tg_zzz_xxxyzz[i] = 2.0 * tg_z_xxxyzz[i] * fxi[i] + 2.0 * tg_zz_xxxyz[i] * fxi[i] + tg_zz_xxxyzz[i] * ra_z[i];

        tg_zzz_xxxzzz[i] = 2.0 * tg_z_xxxzzz[i] * fxi[i] + 3.0 * tg_zz_xxxzz[i] * fxi[i] + tg_zz_xxxzzz[i] * ra_z[i];

        tg_zzz_xxyyyy[i] = 2.0 * tg_z_xxyyyy[i] * fxi[i] + tg_zz_xxyyyy[i] * ra_z[i];

        tg_zzz_xxyyyz[i] = 2.0 * tg_z_xxyyyz[i] * fxi[i] + tg_zz_xxyyy[i] * fxi[i] + tg_zz_xxyyyz[i] * ra_z[i];

        tg_zzz_xxyyzz[i] = 2.0 * tg_z_xxyyzz[i] * fxi[i] + 2.0 * tg_zz_xxyyz[i] * fxi[i] + tg_zz_xxyyzz[i] * ra_z[i];

        tg_zzz_xxyzzz[i] = 2.0 * tg_z_xxyzzz[i] * fxi[i] + 3.0 * tg_zz_xxyzz[i] * fxi[i] + tg_zz_xxyzzz[i] * ra_z[i];

        tg_zzz_xxzzzz[i] = 2.0 * tg_z_xxzzzz[i] * fxi[i] + 4.0 * tg_zz_xxzzz[i] * fxi[i] + tg_zz_xxzzzz[i] * ra_z[i];

        tg_zzz_xyyyyy[i] = 2.0 * tg_z_xyyyyy[i] * fxi[i] + tg_zz_xyyyyy[i] * ra_z[i];

        tg_zzz_xyyyyz[i] = 2.0 * tg_z_xyyyyz[i] * fxi[i] + tg_zz_xyyyy[i] * fxi[i] + tg_zz_xyyyyz[i] * ra_z[i];

        tg_zzz_xyyyzz[i] = 2.0 * tg_z_xyyyzz[i] * fxi[i] + 2.0 * tg_zz_xyyyz[i] * fxi[i] + tg_zz_xyyyzz[i] * ra_z[i];

        tg_zzz_xyyzzz[i] = 2.0 * tg_z_xyyzzz[i] * fxi[i] + 3.0 * tg_zz_xyyzz[i] * fxi[i] + tg_zz_xyyzzz[i] * ra_z[i];

        tg_zzz_xyzzzz[i] = 2.0 * tg_z_xyzzzz[i] * fxi[i] + 4.0 * tg_zz_xyzzz[i] * fxi[i] + tg_zz_xyzzzz[i] * ra_z[i];

        tg_zzz_xzzzzz[i] = 2.0 * tg_z_xzzzzz[i] * fxi[i] + 5.0 * tg_zz_xzzzz[i] * fxi[i] + tg_zz_xzzzzz[i] * ra_z[i];

        tg_zzz_yyyyyy[i] = 2.0 * tg_z_yyyyyy[i] * fxi[i] + tg_zz_yyyyyy[i] * ra_z[i];

        tg_zzz_yyyyyz[i] = 2.0 * tg_z_yyyyyz[i] * fxi[i] + tg_zz_yyyyy[i] * fxi[i] + tg_zz_yyyyyz[i] * ra_z[i];

        tg_zzz_yyyyzz[i] = 2.0 * tg_z_yyyyzz[i] * fxi[i] + 2.0 * tg_zz_yyyyz[i] * fxi[i] + tg_zz_yyyyzz[i] * ra_z[i];

        tg_zzz_yyyzzz[i] = 2.0 * tg_z_yyyzzz[i] * fxi[i] + 3.0 * tg_zz_yyyzz[i] * fxi[i] + tg_zz_yyyzzz[i] * ra_z[i];

        tg_zzz_yyzzzz[i] = 2.0 * tg_z_yyzzzz[i] * fxi[i] + 4.0 * tg_zz_yyzzz[i] * fxi[i] + tg_zz_yyzzzz[i] * ra_z[i];

        tg_zzz_yzzzzz[i] = 2.0 * tg_z_yzzzzz[i] * fxi[i] + 5.0 * tg_zz_yzzzz[i] * fxi[i] + tg_zz_yzzzzz[i] * ra_z[i];

        tg_zzz_zzzzzz[i] = 2.0 * tg_z_zzzzzz[i] * fxi[i] + 6.0 * tg_zz_zzzzz[i] * fxi[i] + tg_zz_zzzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

