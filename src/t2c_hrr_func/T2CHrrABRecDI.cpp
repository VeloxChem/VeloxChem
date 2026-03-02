#include "T2CHrrABRecDI.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_di(CSimdArray<double>& cbuffer, 
            const size_t idx_di,
            const size_t idx_pi,
            const size_t idx_pk,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : PI

    auto t_x_xxxxxx = cbuffer.data(idx_pi);

    auto t_x_xxxxxy = cbuffer.data(idx_pi + 1);

    auto t_x_xxxxxz = cbuffer.data(idx_pi + 2);

    auto t_x_xxxxyy = cbuffer.data(idx_pi + 3);

    auto t_x_xxxxyz = cbuffer.data(idx_pi + 4);

    auto t_x_xxxxzz = cbuffer.data(idx_pi + 5);

    auto t_x_xxxyyy = cbuffer.data(idx_pi + 6);

    auto t_x_xxxyyz = cbuffer.data(idx_pi + 7);

    auto t_x_xxxyzz = cbuffer.data(idx_pi + 8);

    auto t_x_xxxzzz = cbuffer.data(idx_pi + 9);

    auto t_x_xxyyyy = cbuffer.data(idx_pi + 10);

    auto t_x_xxyyyz = cbuffer.data(idx_pi + 11);

    auto t_x_xxyyzz = cbuffer.data(idx_pi + 12);

    auto t_x_xxyzzz = cbuffer.data(idx_pi + 13);

    auto t_x_xxzzzz = cbuffer.data(idx_pi + 14);

    auto t_x_xyyyyy = cbuffer.data(idx_pi + 15);

    auto t_x_xyyyyz = cbuffer.data(idx_pi + 16);

    auto t_x_xyyyzz = cbuffer.data(idx_pi + 17);

    auto t_x_xyyzzz = cbuffer.data(idx_pi + 18);

    auto t_x_xyzzzz = cbuffer.data(idx_pi + 19);

    auto t_x_xzzzzz = cbuffer.data(idx_pi + 20);

    auto t_x_yyyyyy = cbuffer.data(idx_pi + 21);

    auto t_x_yyyyyz = cbuffer.data(idx_pi + 22);

    auto t_x_yyyyzz = cbuffer.data(idx_pi + 23);

    auto t_x_yyyzzz = cbuffer.data(idx_pi + 24);

    auto t_x_yyzzzz = cbuffer.data(idx_pi + 25);

    auto t_x_yzzzzz = cbuffer.data(idx_pi + 26);

    auto t_x_zzzzzz = cbuffer.data(idx_pi + 27);

    auto t_y_xxxxxx = cbuffer.data(idx_pi + 28);

    auto t_y_xxxxxy = cbuffer.data(idx_pi + 29);

    auto t_y_xxxxxz = cbuffer.data(idx_pi + 30);

    auto t_y_xxxxyy = cbuffer.data(idx_pi + 31);

    auto t_y_xxxxyz = cbuffer.data(idx_pi + 32);

    auto t_y_xxxxzz = cbuffer.data(idx_pi + 33);

    auto t_y_xxxyyy = cbuffer.data(idx_pi + 34);

    auto t_y_xxxyyz = cbuffer.data(idx_pi + 35);

    auto t_y_xxxyzz = cbuffer.data(idx_pi + 36);

    auto t_y_xxxzzz = cbuffer.data(idx_pi + 37);

    auto t_y_xxyyyy = cbuffer.data(idx_pi + 38);

    auto t_y_xxyyyz = cbuffer.data(idx_pi + 39);

    auto t_y_xxyyzz = cbuffer.data(idx_pi + 40);

    auto t_y_xxyzzz = cbuffer.data(idx_pi + 41);

    auto t_y_xxzzzz = cbuffer.data(idx_pi + 42);

    auto t_y_xyyyyy = cbuffer.data(idx_pi + 43);

    auto t_y_xyyyyz = cbuffer.data(idx_pi + 44);

    auto t_y_xyyyzz = cbuffer.data(idx_pi + 45);

    auto t_y_xyyzzz = cbuffer.data(idx_pi + 46);

    auto t_y_xyzzzz = cbuffer.data(idx_pi + 47);

    auto t_y_xzzzzz = cbuffer.data(idx_pi + 48);

    auto t_y_yyyyyy = cbuffer.data(idx_pi + 49);

    auto t_y_yyyyyz = cbuffer.data(idx_pi + 50);

    auto t_y_yyyyzz = cbuffer.data(idx_pi + 51);

    auto t_y_yyyzzz = cbuffer.data(idx_pi + 52);

    auto t_y_yyzzzz = cbuffer.data(idx_pi + 53);

    auto t_y_yzzzzz = cbuffer.data(idx_pi + 54);

    auto t_y_zzzzzz = cbuffer.data(idx_pi + 55);

    auto t_z_xxxxxx = cbuffer.data(idx_pi + 56);

    auto t_z_xxxxxy = cbuffer.data(idx_pi + 57);

    auto t_z_xxxxxz = cbuffer.data(idx_pi + 58);

    auto t_z_xxxxyy = cbuffer.data(idx_pi + 59);

    auto t_z_xxxxyz = cbuffer.data(idx_pi + 60);

    auto t_z_xxxxzz = cbuffer.data(idx_pi + 61);

    auto t_z_xxxyyy = cbuffer.data(idx_pi + 62);

    auto t_z_xxxyyz = cbuffer.data(idx_pi + 63);

    auto t_z_xxxyzz = cbuffer.data(idx_pi + 64);

    auto t_z_xxxzzz = cbuffer.data(idx_pi + 65);

    auto t_z_xxyyyy = cbuffer.data(idx_pi + 66);

    auto t_z_xxyyyz = cbuffer.data(idx_pi + 67);

    auto t_z_xxyyzz = cbuffer.data(idx_pi + 68);

    auto t_z_xxyzzz = cbuffer.data(idx_pi + 69);

    auto t_z_xxzzzz = cbuffer.data(idx_pi + 70);

    auto t_z_xyyyyy = cbuffer.data(idx_pi + 71);

    auto t_z_xyyyyz = cbuffer.data(idx_pi + 72);

    auto t_z_xyyyzz = cbuffer.data(idx_pi + 73);

    auto t_z_xyyzzz = cbuffer.data(idx_pi + 74);

    auto t_z_xyzzzz = cbuffer.data(idx_pi + 75);

    auto t_z_xzzzzz = cbuffer.data(idx_pi + 76);

    auto t_z_yyyyyy = cbuffer.data(idx_pi + 77);

    auto t_z_yyyyyz = cbuffer.data(idx_pi + 78);

    auto t_z_yyyyzz = cbuffer.data(idx_pi + 79);

    auto t_z_yyyzzz = cbuffer.data(idx_pi + 80);

    auto t_z_yyzzzz = cbuffer.data(idx_pi + 81);

    auto t_z_yzzzzz = cbuffer.data(idx_pi + 82);

    auto t_z_zzzzzz = cbuffer.data(idx_pi + 83);

    // Set up components of auxiliary buffer : PK

    auto t_x_xxxxxxx = cbuffer.data(idx_pk);

    auto t_x_xxxxxxy = cbuffer.data(idx_pk + 1);

    auto t_x_xxxxxxz = cbuffer.data(idx_pk + 2);

    auto t_x_xxxxxyy = cbuffer.data(idx_pk + 3);

    auto t_x_xxxxxyz = cbuffer.data(idx_pk + 4);

    auto t_x_xxxxxzz = cbuffer.data(idx_pk + 5);

    auto t_x_xxxxyyy = cbuffer.data(idx_pk + 6);

    auto t_x_xxxxyyz = cbuffer.data(idx_pk + 7);

    auto t_x_xxxxyzz = cbuffer.data(idx_pk + 8);

    auto t_x_xxxxzzz = cbuffer.data(idx_pk + 9);

    auto t_x_xxxyyyy = cbuffer.data(idx_pk + 10);

    auto t_x_xxxyyyz = cbuffer.data(idx_pk + 11);

    auto t_x_xxxyyzz = cbuffer.data(idx_pk + 12);

    auto t_x_xxxyzzz = cbuffer.data(idx_pk + 13);

    auto t_x_xxxzzzz = cbuffer.data(idx_pk + 14);

    auto t_x_xxyyyyy = cbuffer.data(idx_pk + 15);

    auto t_x_xxyyyyz = cbuffer.data(idx_pk + 16);

    auto t_x_xxyyyzz = cbuffer.data(idx_pk + 17);

    auto t_x_xxyyzzz = cbuffer.data(idx_pk + 18);

    auto t_x_xxyzzzz = cbuffer.data(idx_pk + 19);

    auto t_x_xxzzzzz = cbuffer.data(idx_pk + 20);

    auto t_x_xyyyyyy = cbuffer.data(idx_pk + 21);

    auto t_x_xyyyyyz = cbuffer.data(idx_pk + 22);

    auto t_x_xyyyyzz = cbuffer.data(idx_pk + 23);

    auto t_x_xyyyzzz = cbuffer.data(idx_pk + 24);

    auto t_x_xyyzzzz = cbuffer.data(idx_pk + 25);

    auto t_x_xyzzzzz = cbuffer.data(idx_pk + 26);

    auto t_x_xzzzzzz = cbuffer.data(idx_pk + 27);

    auto t_y_xxxxxxx = cbuffer.data(idx_pk + 36);

    auto t_y_xxxxxxy = cbuffer.data(idx_pk + 37);

    auto t_y_xxxxxxz = cbuffer.data(idx_pk + 38);

    auto t_y_xxxxxyy = cbuffer.data(idx_pk + 39);

    auto t_y_xxxxxyz = cbuffer.data(idx_pk + 40);

    auto t_y_xxxxxzz = cbuffer.data(idx_pk + 41);

    auto t_y_xxxxyyy = cbuffer.data(idx_pk + 42);

    auto t_y_xxxxyyz = cbuffer.data(idx_pk + 43);

    auto t_y_xxxxyzz = cbuffer.data(idx_pk + 44);

    auto t_y_xxxxzzz = cbuffer.data(idx_pk + 45);

    auto t_y_xxxyyyy = cbuffer.data(idx_pk + 46);

    auto t_y_xxxyyyz = cbuffer.data(idx_pk + 47);

    auto t_y_xxxyyzz = cbuffer.data(idx_pk + 48);

    auto t_y_xxxyzzz = cbuffer.data(idx_pk + 49);

    auto t_y_xxxzzzz = cbuffer.data(idx_pk + 50);

    auto t_y_xxyyyyy = cbuffer.data(idx_pk + 51);

    auto t_y_xxyyyyz = cbuffer.data(idx_pk + 52);

    auto t_y_xxyyyzz = cbuffer.data(idx_pk + 53);

    auto t_y_xxyyzzz = cbuffer.data(idx_pk + 54);

    auto t_y_xxyzzzz = cbuffer.data(idx_pk + 55);

    auto t_y_xxzzzzz = cbuffer.data(idx_pk + 56);

    auto t_y_xyyyyyy = cbuffer.data(idx_pk + 57);

    auto t_y_xyyyyyz = cbuffer.data(idx_pk + 58);

    auto t_y_xyyyyzz = cbuffer.data(idx_pk + 59);

    auto t_y_xyyyzzz = cbuffer.data(idx_pk + 60);

    auto t_y_xyyzzzz = cbuffer.data(idx_pk + 61);

    auto t_y_xyzzzzz = cbuffer.data(idx_pk + 62);

    auto t_y_xzzzzzz = cbuffer.data(idx_pk + 63);

    auto t_y_yyyyyyy = cbuffer.data(idx_pk + 64);

    auto t_y_yyyyyyz = cbuffer.data(idx_pk + 65);

    auto t_y_yyyyyzz = cbuffer.data(idx_pk + 66);

    auto t_y_yyyyzzz = cbuffer.data(idx_pk + 67);

    auto t_y_yyyzzzz = cbuffer.data(idx_pk + 68);

    auto t_y_yyzzzzz = cbuffer.data(idx_pk + 69);

    auto t_y_yzzzzzz = cbuffer.data(idx_pk + 70);

    auto t_z_xxxxxxx = cbuffer.data(idx_pk + 72);

    auto t_z_xxxxxxy = cbuffer.data(idx_pk + 73);

    auto t_z_xxxxxxz = cbuffer.data(idx_pk + 74);

    auto t_z_xxxxxyy = cbuffer.data(idx_pk + 75);

    auto t_z_xxxxxyz = cbuffer.data(idx_pk + 76);

    auto t_z_xxxxxzz = cbuffer.data(idx_pk + 77);

    auto t_z_xxxxyyy = cbuffer.data(idx_pk + 78);

    auto t_z_xxxxyyz = cbuffer.data(idx_pk + 79);

    auto t_z_xxxxyzz = cbuffer.data(idx_pk + 80);

    auto t_z_xxxxzzz = cbuffer.data(idx_pk + 81);

    auto t_z_xxxyyyy = cbuffer.data(idx_pk + 82);

    auto t_z_xxxyyyz = cbuffer.data(idx_pk + 83);

    auto t_z_xxxyyzz = cbuffer.data(idx_pk + 84);

    auto t_z_xxxyzzz = cbuffer.data(idx_pk + 85);

    auto t_z_xxxzzzz = cbuffer.data(idx_pk + 86);

    auto t_z_xxyyyyy = cbuffer.data(idx_pk + 87);

    auto t_z_xxyyyyz = cbuffer.data(idx_pk + 88);

    auto t_z_xxyyyzz = cbuffer.data(idx_pk + 89);

    auto t_z_xxyyzzz = cbuffer.data(idx_pk + 90);

    auto t_z_xxyzzzz = cbuffer.data(idx_pk + 91);

    auto t_z_xxzzzzz = cbuffer.data(idx_pk + 92);

    auto t_z_xyyyyyy = cbuffer.data(idx_pk + 93);

    auto t_z_xyyyyyz = cbuffer.data(idx_pk + 94);

    auto t_z_xyyyyzz = cbuffer.data(idx_pk + 95);

    auto t_z_xyyyzzz = cbuffer.data(idx_pk + 96);

    auto t_z_xyyzzzz = cbuffer.data(idx_pk + 97);

    auto t_z_xyzzzzz = cbuffer.data(idx_pk + 98);

    auto t_z_xzzzzzz = cbuffer.data(idx_pk + 99);

    auto t_z_yyyyyyy = cbuffer.data(idx_pk + 100);

    auto t_z_yyyyyyz = cbuffer.data(idx_pk + 101);

    auto t_z_yyyyyzz = cbuffer.data(idx_pk + 102);

    auto t_z_yyyyzzz = cbuffer.data(idx_pk + 103);

    auto t_z_yyyzzzz = cbuffer.data(idx_pk + 104);

    auto t_z_yyzzzzz = cbuffer.data(idx_pk + 105);

    auto t_z_yzzzzzz = cbuffer.data(idx_pk + 106);

    auto t_z_zzzzzzz = cbuffer.data(idx_pk + 107);

    // Set up components of targeted buffer : DI

    auto t_xx_xxxxxx = cbuffer.data(idx_di);

    auto t_xx_xxxxxy = cbuffer.data(idx_di + 1);

    auto t_xx_xxxxxz = cbuffer.data(idx_di + 2);

    auto t_xx_xxxxyy = cbuffer.data(idx_di + 3);

    auto t_xx_xxxxyz = cbuffer.data(idx_di + 4);

    auto t_xx_xxxxzz = cbuffer.data(idx_di + 5);

    auto t_xx_xxxyyy = cbuffer.data(idx_di + 6);

    auto t_xx_xxxyyz = cbuffer.data(idx_di + 7);

    auto t_xx_xxxyzz = cbuffer.data(idx_di + 8);

    auto t_xx_xxxzzz = cbuffer.data(idx_di + 9);

    auto t_xx_xxyyyy = cbuffer.data(idx_di + 10);

    auto t_xx_xxyyyz = cbuffer.data(idx_di + 11);

    auto t_xx_xxyyzz = cbuffer.data(idx_di + 12);

    auto t_xx_xxyzzz = cbuffer.data(idx_di + 13);

    auto t_xx_xxzzzz = cbuffer.data(idx_di + 14);

    auto t_xx_xyyyyy = cbuffer.data(idx_di + 15);

    auto t_xx_xyyyyz = cbuffer.data(idx_di + 16);

    auto t_xx_xyyyzz = cbuffer.data(idx_di + 17);

    auto t_xx_xyyzzz = cbuffer.data(idx_di + 18);

    auto t_xx_xyzzzz = cbuffer.data(idx_di + 19);

    auto t_xx_xzzzzz = cbuffer.data(idx_di + 20);

    auto t_xx_yyyyyy = cbuffer.data(idx_di + 21);

    auto t_xx_yyyyyz = cbuffer.data(idx_di + 22);

    auto t_xx_yyyyzz = cbuffer.data(idx_di + 23);

    auto t_xx_yyyzzz = cbuffer.data(idx_di + 24);

    auto t_xx_yyzzzz = cbuffer.data(idx_di + 25);

    auto t_xx_yzzzzz = cbuffer.data(idx_di + 26);

    auto t_xx_zzzzzz = cbuffer.data(idx_di + 27);

    auto t_xy_xxxxxx = cbuffer.data(idx_di + 28);

    auto t_xy_xxxxxy = cbuffer.data(idx_di + 29);

    auto t_xy_xxxxxz = cbuffer.data(idx_di + 30);

    auto t_xy_xxxxyy = cbuffer.data(idx_di + 31);

    auto t_xy_xxxxyz = cbuffer.data(idx_di + 32);

    auto t_xy_xxxxzz = cbuffer.data(idx_di + 33);

    auto t_xy_xxxyyy = cbuffer.data(idx_di + 34);

    auto t_xy_xxxyyz = cbuffer.data(idx_di + 35);

    auto t_xy_xxxyzz = cbuffer.data(idx_di + 36);

    auto t_xy_xxxzzz = cbuffer.data(idx_di + 37);

    auto t_xy_xxyyyy = cbuffer.data(idx_di + 38);

    auto t_xy_xxyyyz = cbuffer.data(idx_di + 39);

    auto t_xy_xxyyzz = cbuffer.data(idx_di + 40);

    auto t_xy_xxyzzz = cbuffer.data(idx_di + 41);

    auto t_xy_xxzzzz = cbuffer.data(idx_di + 42);

    auto t_xy_xyyyyy = cbuffer.data(idx_di + 43);

    auto t_xy_xyyyyz = cbuffer.data(idx_di + 44);

    auto t_xy_xyyyzz = cbuffer.data(idx_di + 45);

    auto t_xy_xyyzzz = cbuffer.data(idx_di + 46);

    auto t_xy_xyzzzz = cbuffer.data(idx_di + 47);

    auto t_xy_xzzzzz = cbuffer.data(idx_di + 48);

    auto t_xy_yyyyyy = cbuffer.data(idx_di + 49);

    auto t_xy_yyyyyz = cbuffer.data(idx_di + 50);

    auto t_xy_yyyyzz = cbuffer.data(idx_di + 51);

    auto t_xy_yyyzzz = cbuffer.data(idx_di + 52);

    auto t_xy_yyzzzz = cbuffer.data(idx_di + 53);

    auto t_xy_yzzzzz = cbuffer.data(idx_di + 54);

    auto t_xy_zzzzzz = cbuffer.data(idx_di + 55);

    auto t_xz_xxxxxx = cbuffer.data(idx_di + 56);

    auto t_xz_xxxxxy = cbuffer.data(idx_di + 57);

    auto t_xz_xxxxxz = cbuffer.data(idx_di + 58);

    auto t_xz_xxxxyy = cbuffer.data(idx_di + 59);

    auto t_xz_xxxxyz = cbuffer.data(idx_di + 60);

    auto t_xz_xxxxzz = cbuffer.data(idx_di + 61);

    auto t_xz_xxxyyy = cbuffer.data(idx_di + 62);

    auto t_xz_xxxyyz = cbuffer.data(idx_di + 63);

    auto t_xz_xxxyzz = cbuffer.data(idx_di + 64);

    auto t_xz_xxxzzz = cbuffer.data(idx_di + 65);

    auto t_xz_xxyyyy = cbuffer.data(idx_di + 66);

    auto t_xz_xxyyyz = cbuffer.data(idx_di + 67);

    auto t_xz_xxyyzz = cbuffer.data(idx_di + 68);

    auto t_xz_xxyzzz = cbuffer.data(idx_di + 69);

    auto t_xz_xxzzzz = cbuffer.data(idx_di + 70);

    auto t_xz_xyyyyy = cbuffer.data(idx_di + 71);

    auto t_xz_xyyyyz = cbuffer.data(idx_di + 72);

    auto t_xz_xyyyzz = cbuffer.data(idx_di + 73);

    auto t_xz_xyyzzz = cbuffer.data(idx_di + 74);

    auto t_xz_xyzzzz = cbuffer.data(idx_di + 75);

    auto t_xz_xzzzzz = cbuffer.data(idx_di + 76);

    auto t_xz_yyyyyy = cbuffer.data(idx_di + 77);

    auto t_xz_yyyyyz = cbuffer.data(idx_di + 78);

    auto t_xz_yyyyzz = cbuffer.data(idx_di + 79);

    auto t_xz_yyyzzz = cbuffer.data(idx_di + 80);

    auto t_xz_yyzzzz = cbuffer.data(idx_di + 81);

    auto t_xz_yzzzzz = cbuffer.data(idx_di + 82);

    auto t_xz_zzzzzz = cbuffer.data(idx_di + 83);

    auto t_yy_xxxxxx = cbuffer.data(idx_di + 84);

    auto t_yy_xxxxxy = cbuffer.data(idx_di + 85);

    auto t_yy_xxxxxz = cbuffer.data(idx_di + 86);

    auto t_yy_xxxxyy = cbuffer.data(idx_di + 87);

    auto t_yy_xxxxyz = cbuffer.data(idx_di + 88);

    auto t_yy_xxxxzz = cbuffer.data(idx_di + 89);

    auto t_yy_xxxyyy = cbuffer.data(idx_di + 90);

    auto t_yy_xxxyyz = cbuffer.data(idx_di + 91);

    auto t_yy_xxxyzz = cbuffer.data(idx_di + 92);

    auto t_yy_xxxzzz = cbuffer.data(idx_di + 93);

    auto t_yy_xxyyyy = cbuffer.data(idx_di + 94);

    auto t_yy_xxyyyz = cbuffer.data(idx_di + 95);

    auto t_yy_xxyyzz = cbuffer.data(idx_di + 96);

    auto t_yy_xxyzzz = cbuffer.data(idx_di + 97);

    auto t_yy_xxzzzz = cbuffer.data(idx_di + 98);

    auto t_yy_xyyyyy = cbuffer.data(idx_di + 99);

    auto t_yy_xyyyyz = cbuffer.data(idx_di + 100);

    auto t_yy_xyyyzz = cbuffer.data(idx_di + 101);

    auto t_yy_xyyzzz = cbuffer.data(idx_di + 102);

    auto t_yy_xyzzzz = cbuffer.data(idx_di + 103);

    auto t_yy_xzzzzz = cbuffer.data(idx_di + 104);

    auto t_yy_yyyyyy = cbuffer.data(idx_di + 105);

    auto t_yy_yyyyyz = cbuffer.data(idx_di + 106);

    auto t_yy_yyyyzz = cbuffer.data(idx_di + 107);

    auto t_yy_yyyzzz = cbuffer.data(idx_di + 108);

    auto t_yy_yyzzzz = cbuffer.data(idx_di + 109);

    auto t_yy_yzzzzz = cbuffer.data(idx_di + 110);

    auto t_yy_zzzzzz = cbuffer.data(idx_di + 111);

    auto t_yz_xxxxxx = cbuffer.data(idx_di + 112);

    auto t_yz_xxxxxy = cbuffer.data(idx_di + 113);

    auto t_yz_xxxxxz = cbuffer.data(idx_di + 114);

    auto t_yz_xxxxyy = cbuffer.data(idx_di + 115);

    auto t_yz_xxxxyz = cbuffer.data(idx_di + 116);

    auto t_yz_xxxxzz = cbuffer.data(idx_di + 117);

    auto t_yz_xxxyyy = cbuffer.data(idx_di + 118);

    auto t_yz_xxxyyz = cbuffer.data(idx_di + 119);

    auto t_yz_xxxyzz = cbuffer.data(idx_di + 120);

    auto t_yz_xxxzzz = cbuffer.data(idx_di + 121);

    auto t_yz_xxyyyy = cbuffer.data(idx_di + 122);

    auto t_yz_xxyyyz = cbuffer.data(idx_di + 123);

    auto t_yz_xxyyzz = cbuffer.data(idx_di + 124);

    auto t_yz_xxyzzz = cbuffer.data(idx_di + 125);

    auto t_yz_xxzzzz = cbuffer.data(idx_di + 126);

    auto t_yz_xyyyyy = cbuffer.data(idx_di + 127);

    auto t_yz_xyyyyz = cbuffer.data(idx_di + 128);

    auto t_yz_xyyyzz = cbuffer.data(idx_di + 129);

    auto t_yz_xyyzzz = cbuffer.data(idx_di + 130);

    auto t_yz_xyzzzz = cbuffer.data(idx_di + 131);

    auto t_yz_xzzzzz = cbuffer.data(idx_di + 132);

    auto t_yz_yyyyyy = cbuffer.data(idx_di + 133);

    auto t_yz_yyyyyz = cbuffer.data(idx_di + 134);

    auto t_yz_yyyyzz = cbuffer.data(idx_di + 135);

    auto t_yz_yyyzzz = cbuffer.data(idx_di + 136);

    auto t_yz_yyzzzz = cbuffer.data(idx_di + 137);

    auto t_yz_yzzzzz = cbuffer.data(idx_di + 138);

    auto t_yz_zzzzzz = cbuffer.data(idx_di + 139);

    auto t_zz_xxxxxx = cbuffer.data(idx_di + 140);

    auto t_zz_xxxxxy = cbuffer.data(idx_di + 141);

    auto t_zz_xxxxxz = cbuffer.data(idx_di + 142);

    auto t_zz_xxxxyy = cbuffer.data(idx_di + 143);

    auto t_zz_xxxxyz = cbuffer.data(idx_di + 144);

    auto t_zz_xxxxzz = cbuffer.data(idx_di + 145);

    auto t_zz_xxxyyy = cbuffer.data(idx_di + 146);

    auto t_zz_xxxyyz = cbuffer.data(idx_di + 147);

    auto t_zz_xxxyzz = cbuffer.data(idx_di + 148);

    auto t_zz_xxxzzz = cbuffer.data(idx_di + 149);

    auto t_zz_xxyyyy = cbuffer.data(idx_di + 150);

    auto t_zz_xxyyyz = cbuffer.data(idx_di + 151);

    auto t_zz_xxyyzz = cbuffer.data(idx_di + 152);

    auto t_zz_xxyzzz = cbuffer.data(idx_di + 153);

    auto t_zz_xxzzzz = cbuffer.data(idx_di + 154);

    auto t_zz_xyyyyy = cbuffer.data(idx_di + 155);

    auto t_zz_xyyyyz = cbuffer.data(idx_di + 156);

    auto t_zz_xyyyzz = cbuffer.data(idx_di + 157);

    auto t_zz_xyyzzz = cbuffer.data(idx_di + 158);

    auto t_zz_xyzzzz = cbuffer.data(idx_di + 159);

    auto t_zz_xzzzzz = cbuffer.data(idx_di + 160);

    auto t_zz_yyyyyy = cbuffer.data(idx_di + 161);

    auto t_zz_yyyyyz = cbuffer.data(idx_di + 162);

    auto t_zz_yyyyzz = cbuffer.data(idx_di + 163);

    auto t_zz_yyyzzz = cbuffer.data(idx_di + 164);

    auto t_zz_yyzzzz = cbuffer.data(idx_di + 165);

    auto t_zz_yzzzzz = cbuffer.data(idx_di + 166);

    auto t_zz_zzzzzz = cbuffer.data(idx_di + 167);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_x_xxxxxx, t_x_xxxxxxx, t_x_xxxxxxy, t_x_xxxxxxz, t_x_xxxxxy, t_x_xxxxxyy, t_x_xxxxxyz, t_x_xxxxxz, t_x_xxxxxzz, t_x_xxxxyy, t_x_xxxxyyy, t_x_xxxxyyz, t_x_xxxxyz, t_x_xxxxyzz, t_x_xxxxzz, t_x_xxxxzzz, t_x_xxxyyy, t_x_xxxyyyy, t_x_xxxyyyz, t_x_xxxyyz, t_x_xxxyyzz, t_x_xxxyzz, t_x_xxxyzzz, t_x_xxxzzz, t_x_xxxzzzz, t_x_xxyyyy, t_x_xxyyyyy, t_x_xxyyyyz, t_x_xxyyyz, t_x_xxyyyzz, t_x_xxyyzz, t_x_xxyyzzz, t_x_xxyzzz, t_x_xxyzzzz, t_x_xxzzzz, t_x_xxzzzzz, t_x_xyyyyy, t_x_xyyyyyy, t_x_xyyyyyz, t_x_xyyyyz, t_x_xyyyyzz, t_x_xyyyzz, t_x_xyyyzzz, t_x_xyyzzz, t_x_xyyzzzz, t_x_xyzzzz, t_x_xyzzzzz, t_x_xzzzzz, t_x_xzzzzzz, t_x_yyyyyy, t_x_yyyyyz, t_x_yyyyzz, t_x_yyyzzz, t_x_yyzzzz, t_x_yzzzzz, t_x_zzzzzz, t_xx_xxxxxx, t_xx_xxxxxy, t_xx_xxxxxz, t_xx_xxxxyy, t_xx_xxxxyz, t_xx_xxxxzz, t_xx_xxxyyy, t_xx_xxxyyz, t_xx_xxxyzz, t_xx_xxxzzz, t_xx_xxyyyy, t_xx_xxyyyz, t_xx_xxyyzz, t_xx_xxyzzz, t_xx_xxzzzz, t_xx_xyyyyy, t_xx_xyyyyz, t_xx_xyyyzz, t_xx_xyyzzz, t_xx_xyzzzz, t_xx_xzzzzz, t_xx_yyyyyy, t_xx_yyyyyz, t_xx_yyyyzz, t_xx_yyyzzz, t_xx_yyzzzz, t_xx_yzzzzz, t_xx_zzzzzz, t_xy_xxxxxx, t_xy_xxxxxy, t_xy_xxxxxz, t_xy_xxxxyy, t_xy_xxxxyz, t_xy_xxxxzz, t_xy_xxxyyy, t_xy_xxxyyz, t_xy_xxxyzz, t_xy_xxxzzz, t_xy_xxyyyy, t_xy_xxyyyz, t_xy_xxyyzz, t_xy_xxyzzz, t_xy_xxzzzz, t_xy_xyyyyy, t_xy_xyyyyz, t_xy_xyyyzz, t_xy_xyyzzz, t_xy_xyzzzz, t_xy_xzzzzz, t_xy_yyyyyy, t_xy_yyyyyz, t_xy_yyyyzz, t_xy_yyyzzz, t_xy_yyzzzz, t_xy_yzzzzz, t_xy_zzzzzz, t_xz_xxxxxx, t_xz_xxxxxy, t_xz_xxxxxz, t_xz_xxxxyy, t_xz_xxxxyz, t_xz_xxxxzz, t_xz_xxxyyy, t_xz_xxxyyz, t_xz_xxxyzz, t_xz_xxxzzz, t_xz_xxyyyy, t_xz_xxyyyz, t_xz_xxyyzz, t_xz_xxyzzz, t_xz_xxzzzz, t_xz_xyyyyy, t_xz_xyyyyz, t_xz_xyyyzz, t_xz_xyyzzz, t_xz_xyzzzz, t_xz_xzzzzz, t_xz_yyyyyy, t_xz_yyyyyz, t_xz_yyyyzz, t_xz_yyyzzz, t_xz_yyzzzz, t_xz_yzzzzz, t_xz_zzzzzz, t_y_xxxxxx, t_y_xxxxxxx, t_y_xxxxxxy, t_y_xxxxxxz, t_y_xxxxxy, t_y_xxxxxyy, t_y_xxxxxyz, t_y_xxxxxz, t_y_xxxxxzz, t_y_xxxxyy, t_y_xxxxyyy, t_y_xxxxyyz, t_y_xxxxyz, t_y_xxxxyzz, t_y_xxxxzz, t_y_xxxxzzz, t_y_xxxyyy, t_y_xxxyyyy, t_y_xxxyyyz, t_y_xxxyyz, t_y_xxxyyzz, t_y_xxxyzz, t_y_xxxyzzz, t_y_xxxzzz, t_y_xxxzzzz, t_y_xxyyyy, t_y_xxyyyyy, t_y_xxyyyyz, t_y_xxyyyz, t_y_xxyyyzz, t_y_xxyyzz, t_y_xxyyzzz, t_y_xxyzzz, t_y_xxyzzzz, t_y_xxzzzz, t_y_xxzzzzz, t_y_xyyyyy, t_y_xyyyyyy, t_y_xyyyyyz, t_y_xyyyyz, t_y_xyyyyzz, t_y_xyyyzz, t_y_xyyyzzz, t_y_xyyzzz, t_y_xyyzzzz, t_y_xyzzzz, t_y_xyzzzzz, t_y_xzzzzz, t_y_xzzzzzz, t_y_yyyyyy, t_y_yyyyyyy, t_y_yyyyyyz, t_y_yyyyyz, t_y_yyyyyzz, t_y_yyyyzz, t_y_yyyyzzz, t_y_yyyzzz, t_y_yyyzzzz, t_y_yyzzzz, t_y_yyzzzzz, t_y_yzzzzz, t_y_yzzzzzz, t_y_zzzzzz, t_yy_xxxxxx, t_yy_xxxxxy, t_yy_xxxxxz, t_yy_xxxxyy, t_yy_xxxxyz, t_yy_xxxxzz, t_yy_xxxyyy, t_yy_xxxyyz, t_yy_xxxyzz, t_yy_xxxzzz, t_yy_xxyyyy, t_yy_xxyyyz, t_yy_xxyyzz, t_yy_xxyzzz, t_yy_xxzzzz, t_yy_xyyyyy, t_yy_xyyyyz, t_yy_xyyyzz, t_yy_xyyzzz, t_yy_xyzzzz, t_yy_xzzzzz, t_yy_yyyyyy, t_yy_yyyyyz, t_yy_yyyyzz, t_yy_yyyzzz, t_yy_yyzzzz, t_yy_yzzzzz, t_yy_zzzzzz, t_yz_xxxxxx, t_yz_xxxxxy, t_yz_xxxxxz, t_yz_xxxxyy, t_yz_xxxxyz, t_yz_xxxxzz, t_yz_xxxyyy, t_yz_xxxyyz, t_yz_xxxyzz, t_yz_xxxzzz, t_yz_xxyyyy, t_yz_xxyyyz, t_yz_xxyyzz, t_yz_xxyzzz, t_yz_xxzzzz, t_yz_xyyyyy, t_yz_xyyyyz, t_yz_xyyyzz, t_yz_xyyzzz, t_yz_xyzzzz, t_yz_xzzzzz, t_yz_yyyyyy, t_yz_yyyyyz, t_yz_yyyyzz, t_yz_yyyzzz, t_yz_yyzzzz, t_yz_yzzzzz, t_yz_zzzzzz, t_z_xxxxxx, t_z_xxxxxxx, t_z_xxxxxxy, t_z_xxxxxxz, t_z_xxxxxy, t_z_xxxxxyy, t_z_xxxxxyz, t_z_xxxxxz, t_z_xxxxxzz, t_z_xxxxyy, t_z_xxxxyyy, t_z_xxxxyyz, t_z_xxxxyz, t_z_xxxxyzz, t_z_xxxxzz, t_z_xxxxzzz, t_z_xxxyyy, t_z_xxxyyyy, t_z_xxxyyyz, t_z_xxxyyz, t_z_xxxyyzz, t_z_xxxyzz, t_z_xxxyzzz, t_z_xxxzzz, t_z_xxxzzzz, t_z_xxyyyy, t_z_xxyyyyy, t_z_xxyyyyz, t_z_xxyyyz, t_z_xxyyyzz, t_z_xxyyzz, t_z_xxyyzzz, t_z_xxyzzz, t_z_xxyzzzz, t_z_xxzzzz, t_z_xxzzzzz, t_z_xyyyyy, t_z_xyyyyyy, t_z_xyyyyyz, t_z_xyyyyz, t_z_xyyyyzz, t_z_xyyyzz, t_z_xyyyzzz, t_z_xyyzzz, t_z_xyyzzzz, t_z_xyzzzz, t_z_xyzzzzz, t_z_xzzzzz, t_z_xzzzzzz, t_z_yyyyyy, t_z_yyyyyyy, t_z_yyyyyyz, t_z_yyyyyz, t_z_yyyyyzz, t_z_yyyyzz, t_z_yyyyzzz, t_z_yyyzzz, t_z_yyyzzzz, t_z_yyzzzz, t_z_yyzzzzz, t_z_yzzzzz, t_z_yzzzzzz, t_z_zzzzzz, t_z_zzzzzzz, t_zz_xxxxxx, t_zz_xxxxxy, t_zz_xxxxxz, t_zz_xxxxyy, t_zz_xxxxyz, t_zz_xxxxzz, t_zz_xxxyyy, t_zz_xxxyyz, t_zz_xxxyzz, t_zz_xxxzzz, t_zz_xxyyyy, t_zz_xxyyyz, t_zz_xxyyzz, t_zz_xxyzzz, t_zz_xxzzzz, t_zz_xyyyyy, t_zz_xyyyyz, t_zz_xyyyzz, t_zz_xyyzzz, t_zz_xyzzzz, t_zz_xzzzzz, t_zz_yyyyyy, t_zz_yyyyyz, t_zz_yyyyzz, t_zz_yyyzzz, t_zz_yyzzzz, t_zz_yzzzzz, t_zz_zzzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xx_xxxxxx[i] = -t_x_xxxxxx[i] * ab_x[i] + t_x_xxxxxxx[i];

        t_xx_xxxxxy[i] = -t_x_xxxxxy[i] * ab_x[i] + t_x_xxxxxxy[i];

        t_xx_xxxxxz[i] = -t_x_xxxxxz[i] * ab_x[i] + t_x_xxxxxxz[i];

        t_xx_xxxxyy[i] = -t_x_xxxxyy[i] * ab_x[i] + t_x_xxxxxyy[i];

        t_xx_xxxxyz[i] = -t_x_xxxxyz[i] * ab_x[i] + t_x_xxxxxyz[i];

        t_xx_xxxxzz[i] = -t_x_xxxxzz[i] * ab_x[i] + t_x_xxxxxzz[i];

        t_xx_xxxyyy[i] = -t_x_xxxyyy[i] * ab_x[i] + t_x_xxxxyyy[i];

        t_xx_xxxyyz[i] = -t_x_xxxyyz[i] * ab_x[i] + t_x_xxxxyyz[i];

        t_xx_xxxyzz[i] = -t_x_xxxyzz[i] * ab_x[i] + t_x_xxxxyzz[i];

        t_xx_xxxzzz[i] = -t_x_xxxzzz[i] * ab_x[i] + t_x_xxxxzzz[i];

        t_xx_xxyyyy[i] = -t_x_xxyyyy[i] * ab_x[i] + t_x_xxxyyyy[i];

        t_xx_xxyyyz[i] = -t_x_xxyyyz[i] * ab_x[i] + t_x_xxxyyyz[i];

        t_xx_xxyyzz[i] = -t_x_xxyyzz[i] * ab_x[i] + t_x_xxxyyzz[i];

        t_xx_xxyzzz[i] = -t_x_xxyzzz[i] * ab_x[i] + t_x_xxxyzzz[i];

        t_xx_xxzzzz[i] = -t_x_xxzzzz[i] * ab_x[i] + t_x_xxxzzzz[i];

        t_xx_xyyyyy[i] = -t_x_xyyyyy[i] * ab_x[i] + t_x_xxyyyyy[i];

        t_xx_xyyyyz[i] = -t_x_xyyyyz[i] * ab_x[i] + t_x_xxyyyyz[i];

        t_xx_xyyyzz[i] = -t_x_xyyyzz[i] * ab_x[i] + t_x_xxyyyzz[i];

        t_xx_xyyzzz[i] = -t_x_xyyzzz[i] * ab_x[i] + t_x_xxyyzzz[i];

        t_xx_xyzzzz[i] = -t_x_xyzzzz[i] * ab_x[i] + t_x_xxyzzzz[i];

        t_xx_xzzzzz[i] = -t_x_xzzzzz[i] * ab_x[i] + t_x_xxzzzzz[i];

        t_xx_yyyyyy[i] = -t_x_yyyyyy[i] * ab_x[i] + t_x_xyyyyyy[i];

        t_xx_yyyyyz[i] = -t_x_yyyyyz[i] * ab_x[i] + t_x_xyyyyyz[i];

        t_xx_yyyyzz[i] = -t_x_yyyyzz[i] * ab_x[i] + t_x_xyyyyzz[i];

        t_xx_yyyzzz[i] = -t_x_yyyzzz[i] * ab_x[i] + t_x_xyyyzzz[i];

        t_xx_yyzzzz[i] = -t_x_yyzzzz[i] * ab_x[i] + t_x_xyyzzzz[i];

        t_xx_yzzzzz[i] = -t_x_yzzzzz[i] * ab_x[i] + t_x_xyzzzzz[i];

        t_xx_zzzzzz[i] = -t_x_zzzzzz[i] * ab_x[i] + t_x_xzzzzzz[i];

        t_xy_xxxxxx[i] = -t_y_xxxxxx[i] * ab_x[i] + t_y_xxxxxxx[i];

        t_xy_xxxxxy[i] = -t_y_xxxxxy[i] * ab_x[i] + t_y_xxxxxxy[i];

        t_xy_xxxxxz[i] = -t_y_xxxxxz[i] * ab_x[i] + t_y_xxxxxxz[i];

        t_xy_xxxxyy[i] = -t_y_xxxxyy[i] * ab_x[i] + t_y_xxxxxyy[i];

        t_xy_xxxxyz[i] = -t_y_xxxxyz[i] * ab_x[i] + t_y_xxxxxyz[i];

        t_xy_xxxxzz[i] = -t_y_xxxxzz[i] * ab_x[i] + t_y_xxxxxzz[i];

        t_xy_xxxyyy[i] = -t_y_xxxyyy[i] * ab_x[i] + t_y_xxxxyyy[i];

        t_xy_xxxyyz[i] = -t_y_xxxyyz[i] * ab_x[i] + t_y_xxxxyyz[i];

        t_xy_xxxyzz[i] = -t_y_xxxyzz[i] * ab_x[i] + t_y_xxxxyzz[i];

        t_xy_xxxzzz[i] = -t_y_xxxzzz[i] * ab_x[i] + t_y_xxxxzzz[i];

        t_xy_xxyyyy[i] = -t_y_xxyyyy[i] * ab_x[i] + t_y_xxxyyyy[i];

        t_xy_xxyyyz[i] = -t_y_xxyyyz[i] * ab_x[i] + t_y_xxxyyyz[i];

        t_xy_xxyyzz[i] = -t_y_xxyyzz[i] * ab_x[i] + t_y_xxxyyzz[i];

        t_xy_xxyzzz[i] = -t_y_xxyzzz[i] * ab_x[i] + t_y_xxxyzzz[i];

        t_xy_xxzzzz[i] = -t_y_xxzzzz[i] * ab_x[i] + t_y_xxxzzzz[i];

        t_xy_xyyyyy[i] = -t_y_xyyyyy[i] * ab_x[i] + t_y_xxyyyyy[i];

        t_xy_xyyyyz[i] = -t_y_xyyyyz[i] * ab_x[i] + t_y_xxyyyyz[i];

        t_xy_xyyyzz[i] = -t_y_xyyyzz[i] * ab_x[i] + t_y_xxyyyzz[i];

        t_xy_xyyzzz[i] = -t_y_xyyzzz[i] * ab_x[i] + t_y_xxyyzzz[i];

        t_xy_xyzzzz[i] = -t_y_xyzzzz[i] * ab_x[i] + t_y_xxyzzzz[i];

        t_xy_xzzzzz[i] = -t_y_xzzzzz[i] * ab_x[i] + t_y_xxzzzzz[i];

        t_xy_yyyyyy[i] = -t_y_yyyyyy[i] * ab_x[i] + t_y_xyyyyyy[i];

        t_xy_yyyyyz[i] = -t_y_yyyyyz[i] * ab_x[i] + t_y_xyyyyyz[i];

        t_xy_yyyyzz[i] = -t_y_yyyyzz[i] * ab_x[i] + t_y_xyyyyzz[i];

        t_xy_yyyzzz[i] = -t_y_yyyzzz[i] * ab_x[i] + t_y_xyyyzzz[i];

        t_xy_yyzzzz[i] = -t_y_yyzzzz[i] * ab_x[i] + t_y_xyyzzzz[i];

        t_xy_yzzzzz[i] = -t_y_yzzzzz[i] * ab_x[i] + t_y_xyzzzzz[i];

        t_xy_zzzzzz[i] = -t_y_zzzzzz[i] * ab_x[i] + t_y_xzzzzzz[i];

        t_xz_xxxxxx[i] = -t_z_xxxxxx[i] * ab_x[i] + t_z_xxxxxxx[i];

        t_xz_xxxxxy[i] = -t_z_xxxxxy[i] * ab_x[i] + t_z_xxxxxxy[i];

        t_xz_xxxxxz[i] = -t_z_xxxxxz[i] * ab_x[i] + t_z_xxxxxxz[i];

        t_xz_xxxxyy[i] = -t_z_xxxxyy[i] * ab_x[i] + t_z_xxxxxyy[i];

        t_xz_xxxxyz[i] = -t_z_xxxxyz[i] * ab_x[i] + t_z_xxxxxyz[i];

        t_xz_xxxxzz[i] = -t_z_xxxxzz[i] * ab_x[i] + t_z_xxxxxzz[i];

        t_xz_xxxyyy[i] = -t_z_xxxyyy[i] * ab_x[i] + t_z_xxxxyyy[i];

        t_xz_xxxyyz[i] = -t_z_xxxyyz[i] * ab_x[i] + t_z_xxxxyyz[i];

        t_xz_xxxyzz[i] = -t_z_xxxyzz[i] * ab_x[i] + t_z_xxxxyzz[i];

        t_xz_xxxzzz[i] = -t_z_xxxzzz[i] * ab_x[i] + t_z_xxxxzzz[i];

        t_xz_xxyyyy[i] = -t_z_xxyyyy[i] * ab_x[i] + t_z_xxxyyyy[i];

        t_xz_xxyyyz[i] = -t_z_xxyyyz[i] * ab_x[i] + t_z_xxxyyyz[i];

        t_xz_xxyyzz[i] = -t_z_xxyyzz[i] * ab_x[i] + t_z_xxxyyzz[i];

        t_xz_xxyzzz[i] = -t_z_xxyzzz[i] * ab_x[i] + t_z_xxxyzzz[i];

        t_xz_xxzzzz[i] = -t_z_xxzzzz[i] * ab_x[i] + t_z_xxxzzzz[i];

        t_xz_xyyyyy[i] = -t_z_xyyyyy[i] * ab_x[i] + t_z_xxyyyyy[i];

        t_xz_xyyyyz[i] = -t_z_xyyyyz[i] * ab_x[i] + t_z_xxyyyyz[i];

        t_xz_xyyyzz[i] = -t_z_xyyyzz[i] * ab_x[i] + t_z_xxyyyzz[i];

        t_xz_xyyzzz[i] = -t_z_xyyzzz[i] * ab_x[i] + t_z_xxyyzzz[i];

        t_xz_xyzzzz[i] = -t_z_xyzzzz[i] * ab_x[i] + t_z_xxyzzzz[i];

        t_xz_xzzzzz[i] = -t_z_xzzzzz[i] * ab_x[i] + t_z_xxzzzzz[i];

        t_xz_yyyyyy[i] = -t_z_yyyyyy[i] * ab_x[i] + t_z_xyyyyyy[i];

        t_xz_yyyyyz[i] = -t_z_yyyyyz[i] * ab_x[i] + t_z_xyyyyyz[i];

        t_xz_yyyyzz[i] = -t_z_yyyyzz[i] * ab_x[i] + t_z_xyyyyzz[i];

        t_xz_yyyzzz[i] = -t_z_yyyzzz[i] * ab_x[i] + t_z_xyyyzzz[i];

        t_xz_yyzzzz[i] = -t_z_yyzzzz[i] * ab_x[i] + t_z_xyyzzzz[i];

        t_xz_yzzzzz[i] = -t_z_yzzzzz[i] * ab_x[i] + t_z_xyzzzzz[i];

        t_xz_zzzzzz[i] = -t_z_zzzzzz[i] * ab_x[i] + t_z_xzzzzzz[i];

        t_yy_xxxxxx[i] = -t_y_xxxxxx[i] * ab_y[i] + t_y_xxxxxxy[i];

        t_yy_xxxxxy[i] = -t_y_xxxxxy[i] * ab_y[i] + t_y_xxxxxyy[i];

        t_yy_xxxxxz[i] = -t_y_xxxxxz[i] * ab_y[i] + t_y_xxxxxyz[i];

        t_yy_xxxxyy[i] = -t_y_xxxxyy[i] * ab_y[i] + t_y_xxxxyyy[i];

        t_yy_xxxxyz[i] = -t_y_xxxxyz[i] * ab_y[i] + t_y_xxxxyyz[i];

        t_yy_xxxxzz[i] = -t_y_xxxxzz[i] * ab_y[i] + t_y_xxxxyzz[i];

        t_yy_xxxyyy[i] = -t_y_xxxyyy[i] * ab_y[i] + t_y_xxxyyyy[i];

        t_yy_xxxyyz[i] = -t_y_xxxyyz[i] * ab_y[i] + t_y_xxxyyyz[i];

        t_yy_xxxyzz[i] = -t_y_xxxyzz[i] * ab_y[i] + t_y_xxxyyzz[i];

        t_yy_xxxzzz[i] = -t_y_xxxzzz[i] * ab_y[i] + t_y_xxxyzzz[i];

        t_yy_xxyyyy[i] = -t_y_xxyyyy[i] * ab_y[i] + t_y_xxyyyyy[i];

        t_yy_xxyyyz[i] = -t_y_xxyyyz[i] * ab_y[i] + t_y_xxyyyyz[i];

        t_yy_xxyyzz[i] = -t_y_xxyyzz[i] * ab_y[i] + t_y_xxyyyzz[i];

        t_yy_xxyzzz[i] = -t_y_xxyzzz[i] * ab_y[i] + t_y_xxyyzzz[i];

        t_yy_xxzzzz[i] = -t_y_xxzzzz[i] * ab_y[i] + t_y_xxyzzzz[i];

        t_yy_xyyyyy[i] = -t_y_xyyyyy[i] * ab_y[i] + t_y_xyyyyyy[i];

        t_yy_xyyyyz[i] = -t_y_xyyyyz[i] * ab_y[i] + t_y_xyyyyyz[i];

        t_yy_xyyyzz[i] = -t_y_xyyyzz[i] * ab_y[i] + t_y_xyyyyzz[i];

        t_yy_xyyzzz[i] = -t_y_xyyzzz[i] * ab_y[i] + t_y_xyyyzzz[i];

        t_yy_xyzzzz[i] = -t_y_xyzzzz[i] * ab_y[i] + t_y_xyyzzzz[i];

        t_yy_xzzzzz[i] = -t_y_xzzzzz[i] * ab_y[i] + t_y_xyzzzzz[i];

        t_yy_yyyyyy[i] = -t_y_yyyyyy[i] * ab_y[i] + t_y_yyyyyyy[i];

        t_yy_yyyyyz[i] = -t_y_yyyyyz[i] * ab_y[i] + t_y_yyyyyyz[i];

        t_yy_yyyyzz[i] = -t_y_yyyyzz[i] * ab_y[i] + t_y_yyyyyzz[i];

        t_yy_yyyzzz[i] = -t_y_yyyzzz[i] * ab_y[i] + t_y_yyyyzzz[i];

        t_yy_yyzzzz[i] = -t_y_yyzzzz[i] * ab_y[i] + t_y_yyyzzzz[i];

        t_yy_yzzzzz[i] = -t_y_yzzzzz[i] * ab_y[i] + t_y_yyzzzzz[i];

        t_yy_zzzzzz[i] = -t_y_zzzzzz[i] * ab_y[i] + t_y_yzzzzzz[i];

        t_yz_xxxxxx[i] = -t_z_xxxxxx[i] * ab_y[i] + t_z_xxxxxxy[i];

        t_yz_xxxxxy[i] = -t_z_xxxxxy[i] * ab_y[i] + t_z_xxxxxyy[i];

        t_yz_xxxxxz[i] = -t_z_xxxxxz[i] * ab_y[i] + t_z_xxxxxyz[i];

        t_yz_xxxxyy[i] = -t_z_xxxxyy[i] * ab_y[i] + t_z_xxxxyyy[i];

        t_yz_xxxxyz[i] = -t_z_xxxxyz[i] * ab_y[i] + t_z_xxxxyyz[i];

        t_yz_xxxxzz[i] = -t_z_xxxxzz[i] * ab_y[i] + t_z_xxxxyzz[i];

        t_yz_xxxyyy[i] = -t_z_xxxyyy[i] * ab_y[i] + t_z_xxxyyyy[i];

        t_yz_xxxyyz[i] = -t_z_xxxyyz[i] * ab_y[i] + t_z_xxxyyyz[i];

        t_yz_xxxyzz[i] = -t_z_xxxyzz[i] * ab_y[i] + t_z_xxxyyzz[i];

        t_yz_xxxzzz[i] = -t_z_xxxzzz[i] * ab_y[i] + t_z_xxxyzzz[i];

        t_yz_xxyyyy[i] = -t_z_xxyyyy[i] * ab_y[i] + t_z_xxyyyyy[i];

        t_yz_xxyyyz[i] = -t_z_xxyyyz[i] * ab_y[i] + t_z_xxyyyyz[i];

        t_yz_xxyyzz[i] = -t_z_xxyyzz[i] * ab_y[i] + t_z_xxyyyzz[i];

        t_yz_xxyzzz[i] = -t_z_xxyzzz[i] * ab_y[i] + t_z_xxyyzzz[i];

        t_yz_xxzzzz[i] = -t_z_xxzzzz[i] * ab_y[i] + t_z_xxyzzzz[i];

        t_yz_xyyyyy[i] = -t_z_xyyyyy[i] * ab_y[i] + t_z_xyyyyyy[i];

        t_yz_xyyyyz[i] = -t_z_xyyyyz[i] * ab_y[i] + t_z_xyyyyyz[i];

        t_yz_xyyyzz[i] = -t_z_xyyyzz[i] * ab_y[i] + t_z_xyyyyzz[i];

        t_yz_xyyzzz[i] = -t_z_xyyzzz[i] * ab_y[i] + t_z_xyyyzzz[i];

        t_yz_xyzzzz[i] = -t_z_xyzzzz[i] * ab_y[i] + t_z_xyyzzzz[i];

        t_yz_xzzzzz[i] = -t_z_xzzzzz[i] * ab_y[i] + t_z_xyzzzzz[i];

        t_yz_yyyyyy[i] = -t_z_yyyyyy[i] * ab_y[i] + t_z_yyyyyyy[i];

        t_yz_yyyyyz[i] = -t_z_yyyyyz[i] * ab_y[i] + t_z_yyyyyyz[i];

        t_yz_yyyyzz[i] = -t_z_yyyyzz[i] * ab_y[i] + t_z_yyyyyzz[i];

        t_yz_yyyzzz[i] = -t_z_yyyzzz[i] * ab_y[i] + t_z_yyyyzzz[i];

        t_yz_yyzzzz[i] = -t_z_yyzzzz[i] * ab_y[i] + t_z_yyyzzzz[i];

        t_yz_yzzzzz[i] = -t_z_yzzzzz[i] * ab_y[i] + t_z_yyzzzzz[i];

        t_yz_zzzzzz[i] = -t_z_zzzzzz[i] * ab_y[i] + t_z_yzzzzzz[i];

        t_zz_xxxxxx[i] = -t_z_xxxxxx[i] * ab_z[i] + t_z_xxxxxxz[i];

        t_zz_xxxxxy[i] = -t_z_xxxxxy[i] * ab_z[i] + t_z_xxxxxyz[i];

        t_zz_xxxxxz[i] = -t_z_xxxxxz[i] * ab_z[i] + t_z_xxxxxzz[i];

        t_zz_xxxxyy[i] = -t_z_xxxxyy[i] * ab_z[i] + t_z_xxxxyyz[i];

        t_zz_xxxxyz[i] = -t_z_xxxxyz[i] * ab_z[i] + t_z_xxxxyzz[i];

        t_zz_xxxxzz[i] = -t_z_xxxxzz[i] * ab_z[i] + t_z_xxxxzzz[i];

        t_zz_xxxyyy[i] = -t_z_xxxyyy[i] * ab_z[i] + t_z_xxxyyyz[i];

        t_zz_xxxyyz[i] = -t_z_xxxyyz[i] * ab_z[i] + t_z_xxxyyzz[i];

        t_zz_xxxyzz[i] = -t_z_xxxyzz[i] * ab_z[i] + t_z_xxxyzzz[i];

        t_zz_xxxzzz[i] = -t_z_xxxzzz[i] * ab_z[i] + t_z_xxxzzzz[i];

        t_zz_xxyyyy[i] = -t_z_xxyyyy[i] * ab_z[i] + t_z_xxyyyyz[i];

        t_zz_xxyyyz[i] = -t_z_xxyyyz[i] * ab_z[i] + t_z_xxyyyzz[i];

        t_zz_xxyyzz[i] = -t_z_xxyyzz[i] * ab_z[i] + t_z_xxyyzzz[i];

        t_zz_xxyzzz[i] = -t_z_xxyzzz[i] * ab_z[i] + t_z_xxyzzzz[i];

        t_zz_xxzzzz[i] = -t_z_xxzzzz[i] * ab_z[i] + t_z_xxzzzzz[i];

        t_zz_xyyyyy[i] = -t_z_xyyyyy[i] * ab_z[i] + t_z_xyyyyyz[i];

        t_zz_xyyyyz[i] = -t_z_xyyyyz[i] * ab_z[i] + t_z_xyyyyzz[i];

        t_zz_xyyyzz[i] = -t_z_xyyyzz[i] * ab_z[i] + t_z_xyyyzzz[i];

        t_zz_xyyzzz[i] = -t_z_xyyzzz[i] * ab_z[i] + t_z_xyyzzzz[i];

        t_zz_xyzzzz[i] = -t_z_xyzzzz[i] * ab_z[i] + t_z_xyzzzzz[i];

        t_zz_xzzzzz[i] = -t_z_xzzzzz[i] * ab_z[i] + t_z_xzzzzzz[i];

        t_zz_yyyyyy[i] = -t_z_yyyyyy[i] * ab_z[i] + t_z_yyyyyyz[i];

        t_zz_yyyyyz[i] = -t_z_yyyyyz[i] * ab_z[i] + t_z_yyyyyzz[i];

        t_zz_yyyyzz[i] = -t_z_yyyyzz[i] * ab_z[i] + t_z_yyyyzzz[i];

        t_zz_yyyzzz[i] = -t_z_yyyzzz[i] * ab_z[i] + t_z_yyyzzzz[i];

        t_zz_yyzzzz[i] = -t_z_yyzzzz[i] * ab_z[i] + t_z_yyzzzzz[i];

        t_zz_yzzzzz[i] = -t_z_yzzzzz[i] * ab_z[i] + t_z_yzzzzzz[i];

        t_zz_zzzzzz[i] = -t_z_zzzzzz[i] * ab_z[i] + t_z_zzzzzzz[i];
    }
}

} // t2chrr namespace

