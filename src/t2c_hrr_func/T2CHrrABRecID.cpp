#include "T2CHrrABRecID.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_id(CSimdArray<double>& cbuffer, 
            const size_t idx_id,
            const size_t idx_ip,
            const size_t idx_kp,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : IP

    auto t_xxxxxx_x = cbuffer.data(idx_ip);

    auto t_xxxxxx_y = cbuffer.data(idx_ip + 1);

    auto t_xxxxxx_z = cbuffer.data(idx_ip + 2);

    auto t_xxxxxy_x = cbuffer.data(idx_ip + 3);

    auto t_xxxxxy_y = cbuffer.data(idx_ip + 4);

    auto t_xxxxxy_z = cbuffer.data(idx_ip + 5);

    auto t_xxxxxz_x = cbuffer.data(idx_ip + 6);

    auto t_xxxxxz_y = cbuffer.data(idx_ip + 7);

    auto t_xxxxxz_z = cbuffer.data(idx_ip + 8);

    auto t_xxxxyy_x = cbuffer.data(idx_ip + 9);

    auto t_xxxxyy_y = cbuffer.data(idx_ip + 10);

    auto t_xxxxyy_z = cbuffer.data(idx_ip + 11);

    auto t_xxxxyz_x = cbuffer.data(idx_ip + 12);

    auto t_xxxxyz_y = cbuffer.data(idx_ip + 13);

    auto t_xxxxyz_z = cbuffer.data(idx_ip + 14);

    auto t_xxxxzz_x = cbuffer.data(idx_ip + 15);

    auto t_xxxxzz_y = cbuffer.data(idx_ip + 16);

    auto t_xxxxzz_z = cbuffer.data(idx_ip + 17);

    auto t_xxxyyy_x = cbuffer.data(idx_ip + 18);

    auto t_xxxyyy_y = cbuffer.data(idx_ip + 19);

    auto t_xxxyyy_z = cbuffer.data(idx_ip + 20);

    auto t_xxxyyz_x = cbuffer.data(idx_ip + 21);

    auto t_xxxyyz_y = cbuffer.data(idx_ip + 22);

    auto t_xxxyyz_z = cbuffer.data(idx_ip + 23);

    auto t_xxxyzz_x = cbuffer.data(idx_ip + 24);

    auto t_xxxyzz_y = cbuffer.data(idx_ip + 25);

    auto t_xxxyzz_z = cbuffer.data(idx_ip + 26);

    auto t_xxxzzz_x = cbuffer.data(idx_ip + 27);

    auto t_xxxzzz_y = cbuffer.data(idx_ip + 28);

    auto t_xxxzzz_z = cbuffer.data(idx_ip + 29);

    auto t_xxyyyy_x = cbuffer.data(idx_ip + 30);

    auto t_xxyyyy_y = cbuffer.data(idx_ip + 31);

    auto t_xxyyyy_z = cbuffer.data(idx_ip + 32);

    auto t_xxyyyz_x = cbuffer.data(idx_ip + 33);

    auto t_xxyyyz_y = cbuffer.data(idx_ip + 34);

    auto t_xxyyyz_z = cbuffer.data(idx_ip + 35);

    auto t_xxyyzz_x = cbuffer.data(idx_ip + 36);

    auto t_xxyyzz_y = cbuffer.data(idx_ip + 37);

    auto t_xxyyzz_z = cbuffer.data(idx_ip + 38);

    auto t_xxyzzz_x = cbuffer.data(idx_ip + 39);

    auto t_xxyzzz_y = cbuffer.data(idx_ip + 40);

    auto t_xxyzzz_z = cbuffer.data(idx_ip + 41);

    auto t_xxzzzz_x = cbuffer.data(idx_ip + 42);

    auto t_xxzzzz_y = cbuffer.data(idx_ip + 43);

    auto t_xxzzzz_z = cbuffer.data(idx_ip + 44);

    auto t_xyyyyy_x = cbuffer.data(idx_ip + 45);

    auto t_xyyyyy_y = cbuffer.data(idx_ip + 46);

    auto t_xyyyyy_z = cbuffer.data(idx_ip + 47);

    auto t_xyyyyz_x = cbuffer.data(idx_ip + 48);

    auto t_xyyyyz_y = cbuffer.data(idx_ip + 49);

    auto t_xyyyyz_z = cbuffer.data(idx_ip + 50);

    auto t_xyyyzz_x = cbuffer.data(idx_ip + 51);

    auto t_xyyyzz_y = cbuffer.data(idx_ip + 52);

    auto t_xyyyzz_z = cbuffer.data(idx_ip + 53);

    auto t_xyyzzz_x = cbuffer.data(idx_ip + 54);

    auto t_xyyzzz_y = cbuffer.data(idx_ip + 55);

    auto t_xyyzzz_z = cbuffer.data(idx_ip + 56);

    auto t_xyzzzz_x = cbuffer.data(idx_ip + 57);

    auto t_xyzzzz_y = cbuffer.data(idx_ip + 58);

    auto t_xyzzzz_z = cbuffer.data(idx_ip + 59);

    auto t_xzzzzz_x = cbuffer.data(idx_ip + 60);

    auto t_xzzzzz_y = cbuffer.data(idx_ip + 61);

    auto t_xzzzzz_z = cbuffer.data(idx_ip + 62);

    auto t_yyyyyy_x = cbuffer.data(idx_ip + 63);

    auto t_yyyyyy_y = cbuffer.data(idx_ip + 64);

    auto t_yyyyyy_z = cbuffer.data(idx_ip + 65);

    auto t_yyyyyz_x = cbuffer.data(idx_ip + 66);

    auto t_yyyyyz_y = cbuffer.data(idx_ip + 67);

    auto t_yyyyyz_z = cbuffer.data(idx_ip + 68);

    auto t_yyyyzz_x = cbuffer.data(idx_ip + 69);

    auto t_yyyyzz_y = cbuffer.data(idx_ip + 70);

    auto t_yyyyzz_z = cbuffer.data(idx_ip + 71);

    auto t_yyyzzz_x = cbuffer.data(idx_ip + 72);

    auto t_yyyzzz_y = cbuffer.data(idx_ip + 73);

    auto t_yyyzzz_z = cbuffer.data(idx_ip + 74);

    auto t_yyzzzz_x = cbuffer.data(idx_ip + 75);

    auto t_yyzzzz_y = cbuffer.data(idx_ip + 76);

    auto t_yyzzzz_z = cbuffer.data(idx_ip + 77);

    auto t_yzzzzz_x = cbuffer.data(idx_ip + 78);

    auto t_yzzzzz_y = cbuffer.data(idx_ip + 79);

    auto t_yzzzzz_z = cbuffer.data(idx_ip + 80);

    auto t_zzzzzz_x = cbuffer.data(idx_ip + 81);

    auto t_zzzzzz_y = cbuffer.data(idx_ip + 82);

    auto t_zzzzzz_z = cbuffer.data(idx_ip + 83);

    // Set up components of auxiliary buffer : KP

    auto t_xxxxxxx_x = cbuffer.data(idx_kp);

    auto t_xxxxxxx_y = cbuffer.data(idx_kp + 1);

    auto t_xxxxxxx_z = cbuffer.data(idx_kp + 2);

    auto t_xxxxxxy_x = cbuffer.data(idx_kp + 3);

    auto t_xxxxxxy_y = cbuffer.data(idx_kp + 4);

    auto t_xxxxxxy_z = cbuffer.data(idx_kp + 5);

    auto t_xxxxxxz_x = cbuffer.data(idx_kp + 6);

    auto t_xxxxxxz_y = cbuffer.data(idx_kp + 7);

    auto t_xxxxxxz_z = cbuffer.data(idx_kp + 8);

    auto t_xxxxxyy_x = cbuffer.data(idx_kp + 9);

    auto t_xxxxxyy_y = cbuffer.data(idx_kp + 10);

    auto t_xxxxxyy_z = cbuffer.data(idx_kp + 11);

    auto t_xxxxxyz_x = cbuffer.data(idx_kp + 12);

    auto t_xxxxxyz_y = cbuffer.data(idx_kp + 13);

    auto t_xxxxxyz_z = cbuffer.data(idx_kp + 14);

    auto t_xxxxxzz_x = cbuffer.data(idx_kp + 15);

    auto t_xxxxxzz_y = cbuffer.data(idx_kp + 16);

    auto t_xxxxxzz_z = cbuffer.data(idx_kp + 17);

    auto t_xxxxyyy_x = cbuffer.data(idx_kp + 18);

    auto t_xxxxyyy_y = cbuffer.data(idx_kp + 19);

    auto t_xxxxyyy_z = cbuffer.data(idx_kp + 20);

    auto t_xxxxyyz_x = cbuffer.data(idx_kp + 21);

    auto t_xxxxyyz_y = cbuffer.data(idx_kp + 22);

    auto t_xxxxyyz_z = cbuffer.data(idx_kp + 23);

    auto t_xxxxyzz_x = cbuffer.data(idx_kp + 24);

    auto t_xxxxyzz_y = cbuffer.data(idx_kp + 25);

    auto t_xxxxyzz_z = cbuffer.data(idx_kp + 26);

    auto t_xxxxzzz_x = cbuffer.data(idx_kp + 27);

    auto t_xxxxzzz_y = cbuffer.data(idx_kp + 28);

    auto t_xxxxzzz_z = cbuffer.data(idx_kp + 29);

    auto t_xxxyyyy_x = cbuffer.data(idx_kp + 30);

    auto t_xxxyyyy_y = cbuffer.data(idx_kp + 31);

    auto t_xxxyyyy_z = cbuffer.data(idx_kp + 32);

    auto t_xxxyyyz_x = cbuffer.data(idx_kp + 33);

    auto t_xxxyyyz_y = cbuffer.data(idx_kp + 34);

    auto t_xxxyyyz_z = cbuffer.data(idx_kp + 35);

    auto t_xxxyyzz_x = cbuffer.data(idx_kp + 36);

    auto t_xxxyyzz_y = cbuffer.data(idx_kp + 37);

    auto t_xxxyyzz_z = cbuffer.data(idx_kp + 38);

    auto t_xxxyzzz_x = cbuffer.data(idx_kp + 39);

    auto t_xxxyzzz_y = cbuffer.data(idx_kp + 40);

    auto t_xxxyzzz_z = cbuffer.data(idx_kp + 41);

    auto t_xxxzzzz_x = cbuffer.data(idx_kp + 42);

    auto t_xxxzzzz_y = cbuffer.data(idx_kp + 43);

    auto t_xxxzzzz_z = cbuffer.data(idx_kp + 44);

    auto t_xxyyyyy_x = cbuffer.data(idx_kp + 45);

    auto t_xxyyyyy_y = cbuffer.data(idx_kp + 46);

    auto t_xxyyyyy_z = cbuffer.data(idx_kp + 47);

    auto t_xxyyyyz_x = cbuffer.data(idx_kp + 48);

    auto t_xxyyyyz_y = cbuffer.data(idx_kp + 49);

    auto t_xxyyyyz_z = cbuffer.data(idx_kp + 50);

    auto t_xxyyyzz_x = cbuffer.data(idx_kp + 51);

    auto t_xxyyyzz_y = cbuffer.data(idx_kp + 52);

    auto t_xxyyyzz_z = cbuffer.data(idx_kp + 53);

    auto t_xxyyzzz_x = cbuffer.data(idx_kp + 54);

    auto t_xxyyzzz_y = cbuffer.data(idx_kp + 55);

    auto t_xxyyzzz_z = cbuffer.data(idx_kp + 56);

    auto t_xxyzzzz_x = cbuffer.data(idx_kp + 57);

    auto t_xxyzzzz_y = cbuffer.data(idx_kp + 58);

    auto t_xxyzzzz_z = cbuffer.data(idx_kp + 59);

    auto t_xxzzzzz_x = cbuffer.data(idx_kp + 60);

    auto t_xxzzzzz_y = cbuffer.data(idx_kp + 61);

    auto t_xxzzzzz_z = cbuffer.data(idx_kp + 62);

    auto t_xyyyyyy_x = cbuffer.data(idx_kp + 63);

    auto t_xyyyyyy_y = cbuffer.data(idx_kp + 64);

    auto t_xyyyyyy_z = cbuffer.data(idx_kp + 65);

    auto t_xyyyyyz_x = cbuffer.data(idx_kp + 66);

    auto t_xyyyyyz_y = cbuffer.data(idx_kp + 67);

    auto t_xyyyyyz_z = cbuffer.data(idx_kp + 68);

    auto t_xyyyyzz_x = cbuffer.data(idx_kp + 69);

    auto t_xyyyyzz_y = cbuffer.data(idx_kp + 70);

    auto t_xyyyyzz_z = cbuffer.data(idx_kp + 71);

    auto t_xyyyzzz_x = cbuffer.data(idx_kp + 72);

    auto t_xyyyzzz_y = cbuffer.data(idx_kp + 73);

    auto t_xyyyzzz_z = cbuffer.data(idx_kp + 74);

    auto t_xyyzzzz_x = cbuffer.data(idx_kp + 75);

    auto t_xyyzzzz_y = cbuffer.data(idx_kp + 76);

    auto t_xyyzzzz_z = cbuffer.data(idx_kp + 77);

    auto t_xyzzzzz_x = cbuffer.data(idx_kp + 78);

    auto t_xyzzzzz_y = cbuffer.data(idx_kp + 79);

    auto t_xyzzzzz_z = cbuffer.data(idx_kp + 80);

    auto t_xzzzzzz_x = cbuffer.data(idx_kp + 81);

    auto t_xzzzzzz_y = cbuffer.data(idx_kp + 82);

    auto t_xzzzzzz_z = cbuffer.data(idx_kp + 83);

    auto t_yyyyyyy_y = cbuffer.data(idx_kp + 85);

    auto t_yyyyyyy_z = cbuffer.data(idx_kp + 86);

    auto t_yyyyyyz_y = cbuffer.data(idx_kp + 88);

    auto t_yyyyyyz_z = cbuffer.data(idx_kp + 89);

    auto t_yyyyyzz_y = cbuffer.data(idx_kp + 91);

    auto t_yyyyyzz_z = cbuffer.data(idx_kp + 92);

    auto t_yyyyzzz_y = cbuffer.data(idx_kp + 94);

    auto t_yyyyzzz_z = cbuffer.data(idx_kp + 95);

    auto t_yyyzzzz_y = cbuffer.data(idx_kp + 97);

    auto t_yyyzzzz_z = cbuffer.data(idx_kp + 98);

    auto t_yyzzzzz_y = cbuffer.data(idx_kp + 100);

    auto t_yyzzzzz_z = cbuffer.data(idx_kp + 101);

    auto t_yzzzzzz_y = cbuffer.data(idx_kp + 103);

    auto t_yzzzzzz_z = cbuffer.data(idx_kp + 104);

    auto t_zzzzzzz_z = cbuffer.data(idx_kp + 107);

    // Set up components of targeted buffer : ID

    auto t_xxxxxx_xx = cbuffer.data(idx_id);

    auto t_xxxxxx_xy = cbuffer.data(idx_id + 1);

    auto t_xxxxxx_xz = cbuffer.data(idx_id + 2);

    auto t_xxxxxx_yy = cbuffer.data(idx_id + 3);

    auto t_xxxxxx_yz = cbuffer.data(idx_id + 4);

    auto t_xxxxxx_zz = cbuffer.data(idx_id + 5);

    auto t_xxxxxy_xx = cbuffer.data(idx_id + 6);

    auto t_xxxxxy_xy = cbuffer.data(idx_id + 7);

    auto t_xxxxxy_xz = cbuffer.data(idx_id + 8);

    auto t_xxxxxy_yy = cbuffer.data(idx_id + 9);

    auto t_xxxxxy_yz = cbuffer.data(idx_id + 10);

    auto t_xxxxxy_zz = cbuffer.data(idx_id + 11);

    auto t_xxxxxz_xx = cbuffer.data(idx_id + 12);

    auto t_xxxxxz_xy = cbuffer.data(idx_id + 13);

    auto t_xxxxxz_xz = cbuffer.data(idx_id + 14);

    auto t_xxxxxz_yy = cbuffer.data(idx_id + 15);

    auto t_xxxxxz_yz = cbuffer.data(idx_id + 16);

    auto t_xxxxxz_zz = cbuffer.data(idx_id + 17);

    auto t_xxxxyy_xx = cbuffer.data(idx_id + 18);

    auto t_xxxxyy_xy = cbuffer.data(idx_id + 19);

    auto t_xxxxyy_xz = cbuffer.data(idx_id + 20);

    auto t_xxxxyy_yy = cbuffer.data(idx_id + 21);

    auto t_xxxxyy_yz = cbuffer.data(idx_id + 22);

    auto t_xxxxyy_zz = cbuffer.data(idx_id + 23);

    auto t_xxxxyz_xx = cbuffer.data(idx_id + 24);

    auto t_xxxxyz_xy = cbuffer.data(idx_id + 25);

    auto t_xxxxyz_xz = cbuffer.data(idx_id + 26);

    auto t_xxxxyz_yy = cbuffer.data(idx_id + 27);

    auto t_xxxxyz_yz = cbuffer.data(idx_id + 28);

    auto t_xxxxyz_zz = cbuffer.data(idx_id + 29);

    auto t_xxxxzz_xx = cbuffer.data(idx_id + 30);

    auto t_xxxxzz_xy = cbuffer.data(idx_id + 31);

    auto t_xxxxzz_xz = cbuffer.data(idx_id + 32);

    auto t_xxxxzz_yy = cbuffer.data(idx_id + 33);

    auto t_xxxxzz_yz = cbuffer.data(idx_id + 34);

    auto t_xxxxzz_zz = cbuffer.data(idx_id + 35);

    auto t_xxxyyy_xx = cbuffer.data(idx_id + 36);

    auto t_xxxyyy_xy = cbuffer.data(idx_id + 37);

    auto t_xxxyyy_xz = cbuffer.data(idx_id + 38);

    auto t_xxxyyy_yy = cbuffer.data(idx_id + 39);

    auto t_xxxyyy_yz = cbuffer.data(idx_id + 40);

    auto t_xxxyyy_zz = cbuffer.data(idx_id + 41);

    auto t_xxxyyz_xx = cbuffer.data(idx_id + 42);

    auto t_xxxyyz_xy = cbuffer.data(idx_id + 43);

    auto t_xxxyyz_xz = cbuffer.data(idx_id + 44);

    auto t_xxxyyz_yy = cbuffer.data(idx_id + 45);

    auto t_xxxyyz_yz = cbuffer.data(idx_id + 46);

    auto t_xxxyyz_zz = cbuffer.data(idx_id + 47);

    auto t_xxxyzz_xx = cbuffer.data(idx_id + 48);

    auto t_xxxyzz_xy = cbuffer.data(idx_id + 49);

    auto t_xxxyzz_xz = cbuffer.data(idx_id + 50);

    auto t_xxxyzz_yy = cbuffer.data(idx_id + 51);

    auto t_xxxyzz_yz = cbuffer.data(idx_id + 52);

    auto t_xxxyzz_zz = cbuffer.data(idx_id + 53);

    auto t_xxxzzz_xx = cbuffer.data(idx_id + 54);

    auto t_xxxzzz_xy = cbuffer.data(idx_id + 55);

    auto t_xxxzzz_xz = cbuffer.data(idx_id + 56);

    auto t_xxxzzz_yy = cbuffer.data(idx_id + 57);

    auto t_xxxzzz_yz = cbuffer.data(idx_id + 58);

    auto t_xxxzzz_zz = cbuffer.data(idx_id + 59);

    auto t_xxyyyy_xx = cbuffer.data(idx_id + 60);

    auto t_xxyyyy_xy = cbuffer.data(idx_id + 61);

    auto t_xxyyyy_xz = cbuffer.data(idx_id + 62);

    auto t_xxyyyy_yy = cbuffer.data(idx_id + 63);

    auto t_xxyyyy_yz = cbuffer.data(idx_id + 64);

    auto t_xxyyyy_zz = cbuffer.data(idx_id + 65);

    auto t_xxyyyz_xx = cbuffer.data(idx_id + 66);

    auto t_xxyyyz_xy = cbuffer.data(idx_id + 67);

    auto t_xxyyyz_xz = cbuffer.data(idx_id + 68);

    auto t_xxyyyz_yy = cbuffer.data(idx_id + 69);

    auto t_xxyyyz_yz = cbuffer.data(idx_id + 70);

    auto t_xxyyyz_zz = cbuffer.data(idx_id + 71);

    auto t_xxyyzz_xx = cbuffer.data(idx_id + 72);

    auto t_xxyyzz_xy = cbuffer.data(idx_id + 73);

    auto t_xxyyzz_xz = cbuffer.data(idx_id + 74);

    auto t_xxyyzz_yy = cbuffer.data(idx_id + 75);

    auto t_xxyyzz_yz = cbuffer.data(idx_id + 76);

    auto t_xxyyzz_zz = cbuffer.data(idx_id + 77);

    auto t_xxyzzz_xx = cbuffer.data(idx_id + 78);

    auto t_xxyzzz_xy = cbuffer.data(idx_id + 79);

    auto t_xxyzzz_xz = cbuffer.data(idx_id + 80);

    auto t_xxyzzz_yy = cbuffer.data(idx_id + 81);

    auto t_xxyzzz_yz = cbuffer.data(idx_id + 82);

    auto t_xxyzzz_zz = cbuffer.data(idx_id + 83);

    auto t_xxzzzz_xx = cbuffer.data(idx_id + 84);

    auto t_xxzzzz_xy = cbuffer.data(idx_id + 85);

    auto t_xxzzzz_xz = cbuffer.data(idx_id + 86);

    auto t_xxzzzz_yy = cbuffer.data(idx_id + 87);

    auto t_xxzzzz_yz = cbuffer.data(idx_id + 88);

    auto t_xxzzzz_zz = cbuffer.data(idx_id + 89);

    auto t_xyyyyy_xx = cbuffer.data(idx_id + 90);

    auto t_xyyyyy_xy = cbuffer.data(idx_id + 91);

    auto t_xyyyyy_xz = cbuffer.data(idx_id + 92);

    auto t_xyyyyy_yy = cbuffer.data(idx_id + 93);

    auto t_xyyyyy_yz = cbuffer.data(idx_id + 94);

    auto t_xyyyyy_zz = cbuffer.data(idx_id + 95);

    auto t_xyyyyz_xx = cbuffer.data(idx_id + 96);

    auto t_xyyyyz_xy = cbuffer.data(idx_id + 97);

    auto t_xyyyyz_xz = cbuffer.data(idx_id + 98);

    auto t_xyyyyz_yy = cbuffer.data(idx_id + 99);

    auto t_xyyyyz_yz = cbuffer.data(idx_id + 100);

    auto t_xyyyyz_zz = cbuffer.data(idx_id + 101);

    auto t_xyyyzz_xx = cbuffer.data(idx_id + 102);

    auto t_xyyyzz_xy = cbuffer.data(idx_id + 103);

    auto t_xyyyzz_xz = cbuffer.data(idx_id + 104);

    auto t_xyyyzz_yy = cbuffer.data(idx_id + 105);

    auto t_xyyyzz_yz = cbuffer.data(idx_id + 106);

    auto t_xyyyzz_zz = cbuffer.data(idx_id + 107);

    auto t_xyyzzz_xx = cbuffer.data(idx_id + 108);

    auto t_xyyzzz_xy = cbuffer.data(idx_id + 109);

    auto t_xyyzzz_xz = cbuffer.data(idx_id + 110);

    auto t_xyyzzz_yy = cbuffer.data(idx_id + 111);

    auto t_xyyzzz_yz = cbuffer.data(idx_id + 112);

    auto t_xyyzzz_zz = cbuffer.data(idx_id + 113);

    auto t_xyzzzz_xx = cbuffer.data(idx_id + 114);

    auto t_xyzzzz_xy = cbuffer.data(idx_id + 115);

    auto t_xyzzzz_xz = cbuffer.data(idx_id + 116);

    auto t_xyzzzz_yy = cbuffer.data(idx_id + 117);

    auto t_xyzzzz_yz = cbuffer.data(idx_id + 118);

    auto t_xyzzzz_zz = cbuffer.data(idx_id + 119);

    auto t_xzzzzz_xx = cbuffer.data(idx_id + 120);

    auto t_xzzzzz_xy = cbuffer.data(idx_id + 121);

    auto t_xzzzzz_xz = cbuffer.data(idx_id + 122);

    auto t_xzzzzz_yy = cbuffer.data(idx_id + 123);

    auto t_xzzzzz_yz = cbuffer.data(idx_id + 124);

    auto t_xzzzzz_zz = cbuffer.data(idx_id + 125);

    auto t_yyyyyy_xx = cbuffer.data(idx_id + 126);

    auto t_yyyyyy_xy = cbuffer.data(idx_id + 127);

    auto t_yyyyyy_xz = cbuffer.data(idx_id + 128);

    auto t_yyyyyy_yy = cbuffer.data(idx_id + 129);

    auto t_yyyyyy_yz = cbuffer.data(idx_id + 130);

    auto t_yyyyyy_zz = cbuffer.data(idx_id + 131);

    auto t_yyyyyz_xx = cbuffer.data(idx_id + 132);

    auto t_yyyyyz_xy = cbuffer.data(idx_id + 133);

    auto t_yyyyyz_xz = cbuffer.data(idx_id + 134);

    auto t_yyyyyz_yy = cbuffer.data(idx_id + 135);

    auto t_yyyyyz_yz = cbuffer.data(idx_id + 136);

    auto t_yyyyyz_zz = cbuffer.data(idx_id + 137);

    auto t_yyyyzz_xx = cbuffer.data(idx_id + 138);

    auto t_yyyyzz_xy = cbuffer.data(idx_id + 139);

    auto t_yyyyzz_xz = cbuffer.data(idx_id + 140);

    auto t_yyyyzz_yy = cbuffer.data(idx_id + 141);

    auto t_yyyyzz_yz = cbuffer.data(idx_id + 142);

    auto t_yyyyzz_zz = cbuffer.data(idx_id + 143);

    auto t_yyyzzz_xx = cbuffer.data(idx_id + 144);

    auto t_yyyzzz_xy = cbuffer.data(idx_id + 145);

    auto t_yyyzzz_xz = cbuffer.data(idx_id + 146);

    auto t_yyyzzz_yy = cbuffer.data(idx_id + 147);

    auto t_yyyzzz_yz = cbuffer.data(idx_id + 148);

    auto t_yyyzzz_zz = cbuffer.data(idx_id + 149);

    auto t_yyzzzz_xx = cbuffer.data(idx_id + 150);

    auto t_yyzzzz_xy = cbuffer.data(idx_id + 151);

    auto t_yyzzzz_xz = cbuffer.data(idx_id + 152);

    auto t_yyzzzz_yy = cbuffer.data(idx_id + 153);

    auto t_yyzzzz_yz = cbuffer.data(idx_id + 154);

    auto t_yyzzzz_zz = cbuffer.data(idx_id + 155);

    auto t_yzzzzz_xx = cbuffer.data(idx_id + 156);

    auto t_yzzzzz_xy = cbuffer.data(idx_id + 157);

    auto t_yzzzzz_xz = cbuffer.data(idx_id + 158);

    auto t_yzzzzz_yy = cbuffer.data(idx_id + 159);

    auto t_yzzzzz_yz = cbuffer.data(idx_id + 160);

    auto t_yzzzzz_zz = cbuffer.data(idx_id + 161);

    auto t_zzzzzz_xx = cbuffer.data(idx_id + 162);

    auto t_zzzzzz_xy = cbuffer.data(idx_id + 163);

    auto t_zzzzzz_xz = cbuffer.data(idx_id + 164);

    auto t_zzzzzz_yy = cbuffer.data(idx_id + 165);

    auto t_zzzzzz_yz = cbuffer.data(idx_id + 166);

    auto t_zzzzzz_zz = cbuffer.data(idx_id + 167);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxxxx_x, t_xxxxxx_xx, t_xxxxxx_xy, t_xxxxxx_xz, t_xxxxxx_y, t_xxxxxx_yy, t_xxxxxx_yz, t_xxxxxx_z, t_xxxxxx_zz, t_xxxxxxx_x, t_xxxxxxx_y, t_xxxxxxx_z, t_xxxxxxy_x, t_xxxxxxy_y, t_xxxxxxy_z, t_xxxxxxz_x, t_xxxxxxz_y, t_xxxxxxz_z, t_xxxxxy_x, t_xxxxxy_xx, t_xxxxxy_xy, t_xxxxxy_xz, t_xxxxxy_y, t_xxxxxy_yy, t_xxxxxy_yz, t_xxxxxy_z, t_xxxxxy_zz, t_xxxxxyy_x, t_xxxxxyy_y, t_xxxxxyy_z, t_xxxxxyz_x, t_xxxxxyz_y, t_xxxxxyz_z, t_xxxxxz_x, t_xxxxxz_xx, t_xxxxxz_xy, t_xxxxxz_xz, t_xxxxxz_y, t_xxxxxz_yy, t_xxxxxz_yz, t_xxxxxz_z, t_xxxxxz_zz, t_xxxxxzz_x, t_xxxxxzz_y, t_xxxxxzz_z, t_xxxxyy_x, t_xxxxyy_xx, t_xxxxyy_xy, t_xxxxyy_xz, t_xxxxyy_y, t_xxxxyy_yy, t_xxxxyy_yz, t_xxxxyy_z, t_xxxxyy_zz, t_xxxxyyy_x, t_xxxxyyy_y, t_xxxxyyy_z, t_xxxxyyz_x, t_xxxxyyz_y, t_xxxxyyz_z, t_xxxxyz_x, t_xxxxyz_xx, t_xxxxyz_xy, t_xxxxyz_xz, t_xxxxyz_y, t_xxxxyz_yy, t_xxxxyz_yz, t_xxxxyz_z, t_xxxxyz_zz, t_xxxxyzz_x, t_xxxxyzz_y, t_xxxxyzz_z, t_xxxxzz_x, t_xxxxzz_xx, t_xxxxzz_xy, t_xxxxzz_xz, t_xxxxzz_y, t_xxxxzz_yy, t_xxxxzz_yz, t_xxxxzz_z, t_xxxxzz_zz, t_xxxxzzz_x, t_xxxxzzz_y, t_xxxxzzz_z, t_xxxyyy_x, t_xxxyyy_xx, t_xxxyyy_xy, t_xxxyyy_xz, t_xxxyyy_y, t_xxxyyy_yy, t_xxxyyy_yz, t_xxxyyy_z, t_xxxyyy_zz, t_xxxyyyy_x, t_xxxyyyy_y, t_xxxyyyy_z, t_xxxyyyz_x, t_xxxyyyz_y, t_xxxyyyz_z, t_xxxyyz_x, t_xxxyyz_xx, t_xxxyyz_xy, t_xxxyyz_xz, t_xxxyyz_y, t_xxxyyz_yy, t_xxxyyz_yz, t_xxxyyz_z, t_xxxyyz_zz, t_xxxyyzz_x, t_xxxyyzz_y, t_xxxyyzz_z, t_xxxyzz_x, t_xxxyzz_xx, t_xxxyzz_xy, t_xxxyzz_xz, t_xxxyzz_y, t_xxxyzz_yy, t_xxxyzz_yz, t_xxxyzz_z, t_xxxyzz_zz, t_xxxyzzz_x, t_xxxyzzz_y, t_xxxyzzz_z, t_xxxzzz_x, t_xxxzzz_xx, t_xxxzzz_xy, t_xxxzzz_xz, t_xxxzzz_y, t_xxxzzz_yy, t_xxxzzz_yz, t_xxxzzz_z, t_xxxzzz_zz, t_xxxzzzz_x, t_xxxzzzz_y, t_xxxzzzz_z, t_xxyyyy_x, t_xxyyyy_xx, t_xxyyyy_xy, t_xxyyyy_xz, t_xxyyyy_y, t_xxyyyy_yy, t_xxyyyy_yz, t_xxyyyy_z, t_xxyyyy_zz, t_xxyyyyy_x, t_xxyyyyy_y, t_xxyyyyy_z, t_xxyyyyz_x, t_xxyyyyz_y, t_xxyyyyz_z, t_xxyyyz_x, t_xxyyyz_xx, t_xxyyyz_xy, t_xxyyyz_xz, t_xxyyyz_y, t_xxyyyz_yy, t_xxyyyz_yz, t_xxyyyz_z, t_xxyyyz_zz, t_xxyyyzz_x, t_xxyyyzz_y, t_xxyyyzz_z, t_xxyyzz_x, t_xxyyzz_xx, t_xxyyzz_xy, t_xxyyzz_xz, t_xxyyzz_y, t_xxyyzz_yy, t_xxyyzz_yz, t_xxyyzz_z, t_xxyyzz_zz, t_xxyyzzz_x, t_xxyyzzz_y, t_xxyyzzz_z, t_xxyzzz_x, t_xxyzzz_xx, t_xxyzzz_xy, t_xxyzzz_xz, t_xxyzzz_y, t_xxyzzz_yy, t_xxyzzz_yz, t_xxyzzz_z, t_xxyzzz_zz, t_xxyzzzz_x, t_xxyzzzz_y, t_xxyzzzz_z, t_xxzzzz_x, t_xxzzzz_xx, t_xxzzzz_xy, t_xxzzzz_xz, t_xxzzzz_y, t_xxzzzz_yy, t_xxzzzz_yz, t_xxzzzz_z, t_xxzzzz_zz, t_xxzzzzz_x, t_xxzzzzz_y, t_xxzzzzz_z, t_xyyyyy_x, t_xyyyyy_xx, t_xyyyyy_xy, t_xyyyyy_xz, t_xyyyyy_y, t_xyyyyy_yy, t_xyyyyy_yz, t_xyyyyy_z, t_xyyyyy_zz, t_xyyyyyy_x, t_xyyyyyy_y, t_xyyyyyy_z, t_xyyyyyz_x, t_xyyyyyz_y, t_xyyyyyz_z, t_xyyyyz_x, t_xyyyyz_xx, t_xyyyyz_xy, t_xyyyyz_xz, t_xyyyyz_y, t_xyyyyz_yy, t_xyyyyz_yz, t_xyyyyz_z, t_xyyyyz_zz, t_xyyyyzz_x, t_xyyyyzz_y, t_xyyyyzz_z, t_xyyyzz_x, t_xyyyzz_xx, t_xyyyzz_xy, t_xyyyzz_xz, t_xyyyzz_y, t_xyyyzz_yy, t_xyyyzz_yz, t_xyyyzz_z, t_xyyyzz_zz, t_xyyyzzz_x, t_xyyyzzz_y, t_xyyyzzz_z, t_xyyzzz_x, t_xyyzzz_xx, t_xyyzzz_xy, t_xyyzzz_xz, t_xyyzzz_y, t_xyyzzz_yy, t_xyyzzz_yz, t_xyyzzz_z, t_xyyzzz_zz, t_xyyzzzz_x, t_xyyzzzz_y, t_xyyzzzz_z, t_xyzzzz_x, t_xyzzzz_xx, t_xyzzzz_xy, t_xyzzzz_xz, t_xyzzzz_y, t_xyzzzz_yy, t_xyzzzz_yz, t_xyzzzz_z, t_xyzzzz_zz, t_xyzzzzz_x, t_xyzzzzz_y, t_xyzzzzz_z, t_xzzzzz_x, t_xzzzzz_xx, t_xzzzzz_xy, t_xzzzzz_xz, t_xzzzzz_y, t_xzzzzz_yy, t_xzzzzz_yz, t_xzzzzz_z, t_xzzzzz_zz, t_xzzzzzz_x, t_xzzzzzz_y, t_xzzzzzz_z, t_yyyyyy_x, t_yyyyyy_xx, t_yyyyyy_xy, t_yyyyyy_xz, t_yyyyyy_y, t_yyyyyy_yy, t_yyyyyy_yz, t_yyyyyy_z, t_yyyyyy_zz, t_yyyyyyy_y, t_yyyyyyy_z, t_yyyyyyz_y, t_yyyyyyz_z, t_yyyyyz_x, t_yyyyyz_xx, t_yyyyyz_xy, t_yyyyyz_xz, t_yyyyyz_y, t_yyyyyz_yy, t_yyyyyz_yz, t_yyyyyz_z, t_yyyyyz_zz, t_yyyyyzz_y, t_yyyyyzz_z, t_yyyyzz_x, t_yyyyzz_xx, t_yyyyzz_xy, t_yyyyzz_xz, t_yyyyzz_y, t_yyyyzz_yy, t_yyyyzz_yz, t_yyyyzz_z, t_yyyyzz_zz, t_yyyyzzz_y, t_yyyyzzz_z, t_yyyzzz_x, t_yyyzzz_xx, t_yyyzzz_xy, t_yyyzzz_xz, t_yyyzzz_y, t_yyyzzz_yy, t_yyyzzz_yz, t_yyyzzz_z, t_yyyzzz_zz, t_yyyzzzz_y, t_yyyzzzz_z, t_yyzzzz_x, t_yyzzzz_xx, t_yyzzzz_xy, t_yyzzzz_xz, t_yyzzzz_y, t_yyzzzz_yy, t_yyzzzz_yz, t_yyzzzz_z, t_yyzzzz_zz, t_yyzzzzz_y, t_yyzzzzz_z, t_yzzzzz_x, t_yzzzzz_xx, t_yzzzzz_xy, t_yzzzzz_xz, t_yzzzzz_y, t_yzzzzz_yy, t_yzzzzz_yz, t_yzzzzz_z, t_yzzzzz_zz, t_yzzzzzz_y, t_yzzzzzz_z, t_zzzzzz_x, t_zzzzzz_xx, t_zzzzzz_xy, t_zzzzzz_xz, t_zzzzzz_y, t_zzzzzz_yy, t_zzzzzz_yz, t_zzzzzz_z, t_zzzzzz_zz, t_zzzzzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxxxx_xx[i] = t_xxxxxx_x[i] * ab_x[i] + t_xxxxxxx_x[i];

        t_xxxxxx_xy[i] = t_xxxxxx_y[i] * ab_x[i] + t_xxxxxxx_y[i];

        t_xxxxxx_xz[i] = t_xxxxxx_z[i] * ab_x[i] + t_xxxxxxx_z[i];

        t_xxxxxx_yy[i] = t_xxxxxx_y[i] * ab_y[i] + t_xxxxxxy_y[i];

        t_xxxxxx_yz[i] = t_xxxxxx_z[i] * ab_y[i] + t_xxxxxxy_z[i];

        t_xxxxxx_zz[i] = t_xxxxxx_z[i] * ab_z[i] + t_xxxxxxz_z[i];

        t_xxxxxy_xx[i] = t_xxxxxy_x[i] * ab_x[i] + t_xxxxxxy_x[i];

        t_xxxxxy_xy[i] = t_xxxxxy_y[i] * ab_x[i] + t_xxxxxxy_y[i];

        t_xxxxxy_xz[i] = t_xxxxxy_z[i] * ab_x[i] + t_xxxxxxy_z[i];

        t_xxxxxy_yy[i] = t_xxxxxy_y[i] * ab_y[i] + t_xxxxxyy_y[i];

        t_xxxxxy_yz[i] = t_xxxxxy_z[i] * ab_y[i] + t_xxxxxyy_z[i];

        t_xxxxxy_zz[i] = t_xxxxxy_z[i] * ab_z[i] + t_xxxxxyz_z[i];

        t_xxxxxz_xx[i] = t_xxxxxz_x[i] * ab_x[i] + t_xxxxxxz_x[i];

        t_xxxxxz_xy[i] = t_xxxxxz_y[i] * ab_x[i] + t_xxxxxxz_y[i];

        t_xxxxxz_xz[i] = t_xxxxxz_z[i] * ab_x[i] + t_xxxxxxz_z[i];

        t_xxxxxz_yy[i] = t_xxxxxz_y[i] * ab_y[i] + t_xxxxxyz_y[i];

        t_xxxxxz_yz[i] = t_xxxxxz_z[i] * ab_y[i] + t_xxxxxyz_z[i];

        t_xxxxxz_zz[i] = t_xxxxxz_z[i] * ab_z[i] + t_xxxxxzz_z[i];

        t_xxxxyy_xx[i] = t_xxxxyy_x[i] * ab_x[i] + t_xxxxxyy_x[i];

        t_xxxxyy_xy[i] = t_xxxxyy_y[i] * ab_x[i] + t_xxxxxyy_y[i];

        t_xxxxyy_xz[i] = t_xxxxyy_z[i] * ab_x[i] + t_xxxxxyy_z[i];

        t_xxxxyy_yy[i] = t_xxxxyy_y[i] * ab_y[i] + t_xxxxyyy_y[i];

        t_xxxxyy_yz[i] = t_xxxxyy_z[i] * ab_y[i] + t_xxxxyyy_z[i];

        t_xxxxyy_zz[i] = t_xxxxyy_z[i] * ab_z[i] + t_xxxxyyz_z[i];

        t_xxxxyz_xx[i] = t_xxxxyz_x[i] * ab_x[i] + t_xxxxxyz_x[i];

        t_xxxxyz_xy[i] = t_xxxxyz_y[i] * ab_x[i] + t_xxxxxyz_y[i];

        t_xxxxyz_xz[i] = t_xxxxyz_z[i] * ab_x[i] + t_xxxxxyz_z[i];

        t_xxxxyz_yy[i] = t_xxxxyz_y[i] * ab_y[i] + t_xxxxyyz_y[i];

        t_xxxxyz_yz[i] = t_xxxxyz_z[i] * ab_y[i] + t_xxxxyyz_z[i];

        t_xxxxyz_zz[i] = t_xxxxyz_z[i] * ab_z[i] + t_xxxxyzz_z[i];

        t_xxxxzz_xx[i] = t_xxxxzz_x[i] * ab_x[i] + t_xxxxxzz_x[i];

        t_xxxxzz_xy[i] = t_xxxxzz_y[i] * ab_x[i] + t_xxxxxzz_y[i];

        t_xxxxzz_xz[i] = t_xxxxzz_z[i] * ab_x[i] + t_xxxxxzz_z[i];

        t_xxxxzz_yy[i] = t_xxxxzz_y[i] * ab_y[i] + t_xxxxyzz_y[i];

        t_xxxxzz_yz[i] = t_xxxxzz_z[i] * ab_y[i] + t_xxxxyzz_z[i];

        t_xxxxzz_zz[i] = t_xxxxzz_z[i] * ab_z[i] + t_xxxxzzz_z[i];

        t_xxxyyy_xx[i] = t_xxxyyy_x[i] * ab_x[i] + t_xxxxyyy_x[i];

        t_xxxyyy_xy[i] = t_xxxyyy_y[i] * ab_x[i] + t_xxxxyyy_y[i];

        t_xxxyyy_xz[i] = t_xxxyyy_z[i] * ab_x[i] + t_xxxxyyy_z[i];

        t_xxxyyy_yy[i] = t_xxxyyy_y[i] * ab_y[i] + t_xxxyyyy_y[i];

        t_xxxyyy_yz[i] = t_xxxyyy_z[i] * ab_y[i] + t_xxxyyyy_z[i];

        t_xxxyyy_zz[i] = t_xxxyyy_z[i] * ab_z[i] + t_xxxyyyz_z[i];

        t_xxxyyz_xx[i] = t_xxxyyz_x[i] * ab_x[i] + t_xxxxyyz_x[i];

        t_xxxyyz_xy[i] = t_xxxyyz_y[i] * ab_x[i] + t_xxxxyyz_y[i];

        t_xxxyyz_xz[i] = t_xxxyyz_z[i] * ab_x[i] + t_xxxxyyz_z[i];

        t_xxxyyz_yy[i] = t_xxxyyz_y[i] * ab_y[i] + t_xxxyyyz_y[i];

        t_xxxyyz_yz[i] = t_xxxyyz_z[i] * ab_y[i] + t_xxxyyyz_z[i];

        t_xxxyyz_zz[i] = t_xxxyyz_z[i] * ab_z[i] + t_xxxyyzz_z[i];

        t_xxxyzz_xx[i] = t_xxxyzz_x[i] * ab_x[i] + t_xxxxyzz_x[i];

        t_xxxyzz_xy[i] = t_xxxyzz_y[i] * ab_x[i] + t_xxxxyzz_y[i];

        t_xxxyzz_xz[i] = t_xxxyzz_z[i] * ab_x[i] + t_xxxxyzz_z[i];

        t_xxxyzz_yy[i] = t_xxxyzz_y[i] * ab_y[i] + t_xxxyyzz_y[i];

        t_xxxyzz_yz[i] = t_xxxyzz_z[i] * ab_y[i] + t_xxxyyzz_z[i];

        t_xxxyzz_zz[i] = t_xxxyzz_z[i] * ab_z[i] + t_xxxyzzz_z[i];

        t_xxxzzz_xx[i] = t_xxxzzz_x[i] * ab_x[i] + t_xxxxzzz_x[i];

        t_xxxzzz_xy[i] = t_xxxzzz_y[i] * ab_x[i] + t_xxxxzzz_y[i];

        t_xxxzzz_xz[i] = t_xxxzzz_z[i] * ab_x[i] + t_xxxxzzz_z[i];

        t_xxxzzz_yy[i] = t_xxxzzz_y[i] * ab_y[i] + t_xxxyzzz_y[i];

        t_xxxzzz_yz[i] = t_xxxzzz_z[i] * ab_y[i] + t_xxxyzzz_z[i];

        t_xxxzzz_zz[i] = t_xxxzzz_z[i] * ab_z[i] + t_xxxzzzz_z[i];

        t_xxyyyy_xx[i] = t_xxyyyy_x[i] * ab_x[i] + t_xxxyyyy_x[i];

        t_xxyyyy_xy[i] = t_xxyyyy_y[i] * ab_x[i] + t_xxxyyyy_y[i];

        t_xxyyyy_xz[i] = t_xxyyyy_z[i] * ab_x[i] + t_xxxyyyy_z[i];

        t_xxyyyy_yy[i] = t_xxyyyy_y[i] * ab_y[i] + t_xxyyyyy_y[i];

        t_xxyyyy_yz[i] = t_xxyyyy_z[i] * ab_y[i] + t_xxyyyyy_z[i];

        t_xxyyyy_zz[i] = t_xxyyyy_z[i] * ab_z[i] + t_xxyyyyz_z[i];

        t_xxyyyz_xx[i] = t_xxyyyz_x[i] * ab_x[i] + t_xxxyyyz_x[i];

        t_xxyyyz_xy[i] = t_xxyyyz_y[i] * ab_x[i] + t_xxxyyyz_y[i];

        t_xxyyyz_xz[i] = t_xxyyyz_z[i] * ab_x[i] + t_xxxyyyz_z[i];

        t_xxyyyz_yy[i] = t_xxyyyz_y[i] * ab_y[i] + t_xxyyyyz_y[i];

        t_xxyyyz_yz[i] = t_xxyyyz_z[i] * ab_y[i] + t_xxyyyyz_z[i];

        t_xxyyyz_zz[i] = t_xxyyyz_z[i] * ab_z[i] + t_xxyyyzz_z[i];

        t_xxyyzz_xx[i] = t_xxyyzz_x[i] * ab_x[i] + t_xxxyyzz_x[i];

        t_xxyyzz_xy[i] = t_xxyyzz_y[i] * ab_x[i] + t_xxxyyzz_y[i];

        t_xxyyzz_xz[i] = t_xxyyzz_z[i] * ab_x[i] + t_xxxyyzz_z[i];

        t_xxyyzz_yy[i] = t_xxyyzz_y[i] * ab_y[i] + t_xxyyyzz_y[i];

        t_xxyyzz_yz[i] = t_xxyyzz_z[i] * ab_y[i] + t_xxyyyzz_z[i];

        t_xxyyzz_zz[i] = t_xxyyzz_z[i] * ab_z[i] + t_xxyyzzz_z[i];

        t_xxyzzz_xx[i] = t_xxyzzz_x[i] * ab_x[i] + t_xxxyzzz_x[i];

        t_xxyzzz_xy[i] = t_xxyzzz_y[i] * ab_x[i] + t_xxxyzzz_y[i];

        t_xxyzzz_xz[i] = t_xxyzzz_z[i] * ab_x[i] + t_xxxyzzz_z[i];

        t_xxyzzz_yy[i] = t_xxyzzz_y[i] * ab_y[i] + t_xxyyzzz_y[i];

        t_xxyzzz_yz[i] = t_xxyzzz_z[i] * ab_y[i] + t_xxyyzzz_z[i];

        t_xxyzzz_zz[i] = t_xxyzzz_z[i] * ab_z[i] + t_xxyzzzz_z[i];

        t_xxzzzz_xx[i] = t_xxzzzz_x[i] * ab_x[i] + t_xxxzzzz_x[i];

        t_xxzzzz_xy[i] = t_xxzzzz_y[i] * ab_x[i] + t_xxxzzzz_y[i];

        t_xxzzzz_xz[i] = t_xxzzzz_z[i] * ab_x[i] + t_xxxzzzz_z[i];

        t_xxzzzz_yy[i] = t_xxzzzz_y[i] * ab_y[i] + t_xxyzzzz_y[i];

        t_xxzzzz_yz[i] = t_xxzzzz_z[i] * ab_y[i] + t_xxyzzzz_z[i];

        t_xxzzzz_zz[i] = t_xxzzzz_z[i] * ab_z[i] + t_xxzzzzz_z[i];

        t_xyyyyy_xx[i] = t_xyyyyy_x[i] * ab_x[i] + t_xxyyyyy_x[i];

        t_xyyyyy_xy[i] = t_xyyyyy_y[i] * ab_x[i] + t_xxyyyyy_y[i];

        t_xyyyyy_xz[i] = t_xyyyyy_z[i] * ab_x[i] + t_xxyyyyy_z[i];

        t_xyyyyy_yy[i] = t_xyyyyy_y[i] * ab_y[i] + t_xyyyyyy_y[i];

        t_xyyyyy_yz[i] = t_xyyyyy_z[i] * ab_y[i] + t_xyyyyyy_z[i];

        t_xyyyyy_zz[i] = t_xyyyyy_z[i] * ab_z[i] + t_xyyyyyz_z[i];

        t_xyyyyz_xx[i] = t_xyyyyz_x[i] * ab_x[i] + t_xxyyyyz_x[i];

        t_xyyyyz_xy[i] = t_xyyyyz_y[i] * ab_x[i] + t_xxyyyyz_y[i];

        t_xyyyyz_xz[i] = t_xyyyyz_z[i] * ab_x[i] + t_xxyyyyz_z[i];

        t_xyyyyz_yy[i] = t_xyyyyz_y[i] * ab_y[i] + t_xyyyyyz_y[i];

        t_xyyyyz_yz[i] = t_xyyyyz_z[i] * ab_y[i] + t_xyyyyyz_z[i];

        t_xyyyyz_zz[i] = t_xyyyyz_z[i] * ab_z[i] + t_xyyyyzz_z[i];

        t_xyyyzz_xx[i] = t_xyyyzz_x[i] * ab_x[i] + t_xxyyyzz_x[i];

        t_xyyyzz_xy[i] = t_xyyyzz_y[i] * ab_x[i] + t_xxyyyzz_y[i];

        t_xyyyzz_xz[i] = t_xyyyzz_z[i] * ab_x[i] + t_xxyyyzz_z[i];

        t_xyyyzz_yy[i] = t_xyyyzz_y[i] * ab_y[i] + t_xyyyyzz_y[i];

        t_xyyyzz_yz[i] = t_xyyyzz_z[i] * ab_y[i] + t_xyyyyzz_z[i];

        t_xyyyzz_zz[i] = t_xyyyzz_z[i] * ab_z[i] + t_xyyyzzz_z[i];

        t_xyyzzz_xx[i] = t_xyyzzz_x[i] * ab_x[i] + t_xxyyzzz_x[i];

        t_xyyzzz_xy[i] = t_xyyzzz_y[i] * ab_x[i] + t_xxyyzzz_y[i];

        t_xyyzzz_xz[i] = t_xyyzzz_z[i] * ab_x[i] + t_xxyyzzz_z[i];

        t_xyyzzz_yy[i] = t_xyyzzz_y[i] * ab_y[i] + t_xyyyzzz_y[i];

        t_xyyzzz_yz[i] = t_xyyzzz_z[i] * ab_y[i] + t_xyyyzzz_z[i];

        t_xyyzzz_zz[i] = t_xyyzzz_z[i] * ab_z[i] + t_xyyzzzz_z[i];

        t_xyzzzz_xx[i] = t_xyzzzz_x[i] * ab_x[i] + t_xxyzzzz_x[i];

        t_xyzzzz_xy[i] = t_xyzzzz_y[i] * ab_x[i] + t_xxyzzzz_y[i];

        t_xyzzzz_xz[i] = t_xyzzzz_z[i] * ab_x[i] + t_xxyzzzz_z[i];

        t_xyzzzz_yy[i] = t_xyzzzz_y[i] * ab_y[i] + t_xyyzzzz_y[i];

        t_xyzzzz_yz[i] = t_xyzzzz_z[i] * ab_y[i] + t_xyyzzzz_z[i];

        t_xyzzzz_zz[i] = t_xyzzzz_z[i] * ab_z[i] + t_xyzzzzz_z[i];

        t_xzzzzz_xx[i] = t_xzzzzz_x[i] * ab_x[i] + t_xxzzzzz_x[i];

        t_xzzzzz_xy[i] = t_xzzzzz_y[i] * ab_x[i] + t_xxzzzzz_y[i];

        t_xzzzzz_xz[i] = t_xzzzzz_z[i] * ab_x[i] + t_xxzzzzz_z[i];

        t_xzzzzz_yy[i] = t_xzzzzz_y[i] * ab_y[i] + t_xyzzzzz_y[i];

        t_xzzzzz_yz[i] = t_xzzzzz_z[i] * ab_y[i] + t_xyzzzzz_z[i];

        t_xzzzzz_zz[i] = t_xzzzzz_z[i] * ab_z[i] + t_xzzzzzz_z[i];

        t_yyyyyy_xx[i] = t_yyyyyy_x[i] * ab_x[i] + t_xyyyyyy_x[i];

        t_yyyyyy_xy[i] = t_yyyyyy_y[i] * ab_x[i] + t_xyyyyyy_y[i];

        t_yyyyyy_xz[i] = t_yyyyyy_z[i] * ab_x[i] + t_xyyyyyy_z[i];

        t_yyyyyy_yy[i] = t_yyyyyy_y[i] * ab_y[i] + t_yyyyyyy_y[i];

        t_yyyyyy_yz[i] = t_yyyyyy_z[i] * ab_y[i] + t_yyyyyyy_z[i];

        t_yyyyyy_zz[i] = t_yyyyyy_z[i] * ab_z[i] + t_yyyyyyz_z[i];

        t_yyyyyz_xx[i] = t_yyyyyz_x[i] * ab_x[i] + t_xyyyyyz_x[i];

        t_yyyyyz_xy[i] = t_yyyyyz_y[i] * ab_x[i] + t_xyyyyyz_y[i];

        t_yyyyyz_xz[i] = t_yyyyyz_z[i] * ab_x[i] + t_xyyyyyz_z[i];

        t_yyyyyz_yy[i] = t_yyyyyz_y[i] * ab_y[i] + t_yyyyyyz_y[i];

        t_yyyyyz_yz[i] = t_yyyyyz_z[i] * ab_y[i] + t_yyyyyyz_z[i];

        t_yyyyyz_zz[i] = t_yyyyyz_z[i] * ab_z[i] + t_yyyyyzz_z[i];

        t_yyyyzz_xx[i] = t_yyyyzz_x[i] * ab_x[i] + t_xyyyyzz_x[i];

        t_yyyyzz_xy[i] = t_yyyyzz_y[i] * ab_x[i] + t_xyyyyzz_y[i];

        t_yyyyzz_xz[i] = t_yyyyzz_z[i] * ab_x[i] + t_xyyyyzz_z[i];

        t_yyyyzz_yy[i] = t_yyyyzz_y[i] * ab_y[i] + t_yyyyyzz_y[i];

        t_yyyyzz_yz[i] = t_yyyyzz_z[i] * ab_y[i] + t_yyyyyzz_z[i];

        t_yyyyzz_zz[i] = t_yyyyzz_z[i] * ab_z[i] + t_yyyyzzz_z[i];

        t_yyyzzz_xx[i] = t_yyyzzz_x[i] * ab_x[i] + t_xyyyzzz_x[i];

        t_yyyzzz_xy[i] = t_yyyzzz_y[i] * ab_x[i] + t_xyyyzzz_y[i];

        t_yyyzzz_xz[i] = t_yyyzzz_z[i] * ab_x[i] + t_xyyyzzz_z[i];

        t_yyyzzz_yy[i] = t_yyyzzz_y[i] * ab_y[i] + t_yyyyzzz_y[i];

        t_yyyzzz_yz[i] = t_yyyzzz_z[i] * ab_y[i] + t_yyyyzzz_z[i];

        t_yyyzzz_zz[i] = t_yyyzzz_z[i] * ab_z[i] + t_yyyzzzz_z[i];

        t_yyzzzz_xx[i] = t_yyzzzz_x[i] * ab_x[i] + t_xyyzzzz_x[i];

        t_yyzzzz_xy[i] = t_yyzzzz_y[i] * ab_x[i] + t_xyyzzzz_y[i];

        t_yyzzzz_xz[i] = t_yyzzzz_z[i] * ab_x[i] + t_xyyzzzz_z[i];

        t_yyzzzz_yy[i] = t_yyzzzz_y[i] * ab_y[i] + t_yyyzzzz_y[i];

        t_yyzzzz_yz[i] = t_yyzzzz_z[i] * ab_y[i] + t_yyyzzzz_z[i];

        t_yyzzzz_zz[i] = t_yyzzzz_z[i] * ab_z[i] + t_yyzzzzz_z[i];

        t_yzzzzz_xx[i] = t_yzzzzz_x[i] * ab_x[i] + t_xyzzzzz_x[i];

        t_yzzzzz_xy[i] = t_yzzzzz_y[i] * ab_x[i] + t_xyzzzzz_y[i];

        t_yzzzzz_xz[i] = t_yzzzzz_z[i] * ab_x[i] + t_xyzzzzz_z[i];

        t_yzzzzz_yy[i] = t_yzzzzz_y[i] * ab_y[i] + t_yyzzzzz_y[i];

        t_yzzzzz_yz[i] = t_yzzzzz_z[i] * ab_y[i] + t_yyzzzzz_z[i];

        t_yzzzzz_zz[i] = t_yzzzzz_z[i] * ab_z[i] + t_yzzzzzz_z[i];

        t_zzzzzz_xx[i] = t_zzzzzz_x[i] * ab_x[i] + t_xzzzzzz_x[i];

        t_zzzzzz_xy[i] = t_zzzzzz_y[i] * ab_x[i] + t_xzzzzzz_y[i];

        t_zzzzzz_xz[i] = t_zzzzzz_z[i] * ab_x[i] + t_xzzzzzz_z[i];

        t_zzzzzz_yy[i] = t_zzzzzz_y[i] * ab_y[i] + t_yzzzzzz_y[i];

        t_zzzzzz_yz[i] = t_zzzzzz_z[i] * ab_y[i] + t_yzzzzzz_z[i];

        t_zzzzzz_zz[i] = t_zzzzzz_z[i] * ab_z[i] + t_zzzzzzz_z[i];
    }
}

} // t2chrr namespace

