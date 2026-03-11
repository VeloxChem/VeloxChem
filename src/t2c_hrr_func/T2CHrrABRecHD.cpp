#include "T2CHrrABRecHD.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_hd(CSimdArray<double>& cbuffer, 
            const size_t idx_hd,
            const size_t idx_hp,
            const size_t idx_ip,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : HP

    auto t_xxxxx_x = cbuffer.data(idx_hp);

    auto t_xxxxx_y = cbuffer.data(idx_hp + 1);

    auto t_xxxxx_z = cbuffer.data(idx_hp + 2);

    auto t_xxxxy_x = cbuffer.data(idx_hp + 3);

    auto t_xxxxy_y = cbuffer.data(idx_hp + 4);

    auto t_xxxxy_z = cbuffer.data(idx_hp + 5);

    auto t_xxxxz_x = cbuffer.data(idx_hp + 6);

    auto t_xxxxz_y = cbuffer.data(idx_hp + 7);

    auto t_xxxxz_z = cbuffer.data(idx_hp + 8);

    auto t_xxxyy_x = cbuffer.data(idx_hp + 9);

    auto t_xxxyy_y = cbuffer.data(idx_hp + 10);

    auto t_xxxyy_z = cbuffer.data(idx_hp + 11);

    auto t_xxxyz_x = cbuffer.data(idx_hp + 12);

    auto t_xxxyz_y = cbuffer.data(idx_hp + 13);

    auto t_xxxyz_z = cbuffer.data(idx_hp + 14);

    auto t_xxxzz_x = cbuffer.data(idx_hp + 15);

    auto t_xxxzz_y = cbuffer.data(idx_hp + 16);

    auto t_xxxzz_z = cbuffer.data(idx_hp + 17);

    auto t_xxyyy_x = cbuffer.data(idx_hp + 18);

    auto t_xxyyy_y = cbuffer.data(idx_hp + 19);

    auto t_xxyyy_z = cbuffer.data(idx_hp + 20);

    auto t_xxyyz_x = cbuffer.data(idx_hp + 21);

    auto t_xxyyz_y = cbuffer.data(idx_hp + 22);

    auto t_xxyyz_z = cbuffer.data(idx_hp + 23);

    auto t_xxyzz_x = cbuffer.data(idx_hp + 24);

    auto t_xxyzz_y = cbuffer.data(idx_hp + 25);

    auto t_xxyzz_z = cbuffer.data(idx_hp + 26);

    auto t_xxzzz_x = cbuffer.data(idx_hp + 27);

    auto t_xxzzz_y = cbuffer.data(idx_hp + 28);

    auto t_xxzzz_z = cbuffer.data(idx_hp + 29);

    auto t_xyyyy_x = cbuffer.data(idx_hp + 30);

    auto t_xyyyy_y = cbuffer.data(idx_hp + 31);

    auto t_xyyyy_z = cbuffer.data(idx_hp + 32);

    auto t_xyyyz_x = cbuffer.data(idx_hp + 33);

    auto t_xyyyz_y = cbuffer.data(idx_hp + 34);

    auto t_xyyyz_z = cbuffer.data(idx_hp + 35);

    auto t_xyyzz_x = cbuffer.data(idx_hp + 36);

    auto t_xyyzz_y = cbuffer.data(idx_hp + 37);

    auto t_xyyzz_z = cbuffer.data(idx_hp + 38);

    auto t_xyzzz_x = cbuffer.data(idx_hp + 39);

    auto t_xyzzz_y = cbuffer.data(idx_hp + 40);

    auto t_xyzzz_z = cbuffer.data(idx_hp + 41);

    auto t_xzzzz_x = cbuffer.data(idx_hp + 42);

    auto t_xzzzz_y = cbuffer.data(idx_hp + 43);

    auto t_xzzzz_z = cbuffer.data(idx_hp + 44);

    auto t_yyyyy_x = cbuffer.data(idx_hp + 45);

    auto t_yyyyy_y = cbuffer.data(idx_hp + 46);

    auto t_yyyyy_z = cbuffer.data(idx_hp + 47);

    auto t_yyyyz_x = cbuffer.data(idx_hp + 48);

    auto t_yyyyz_y = cbuffer.data(idx_hp + 49);

    auto t_yyyyz_z = cbuffer.data(idx_hp + 50);

    auto t_yyyzz_x = cbuffer.data(idx_hp + 51);

    auto t_yyyzz_y = cbuffer.data(idx_hp + 52);

    auto t_yyyzz_z = cbuffer.data(idx_hp + 53);

    auto t_yyzzz_x = cbuffer.data(idx_hp + 54);

    auto t_yyzzz_y = cbuffer.data(idx_hp + 55);

    auto t_yyzzz_z = cbuffer.data(idx_hp + 56);

    auto t_yzzzz_x = cbuffer.data(idx_hp + 57);

    auto t_yzzzz_y = cbuffer.data(idx_hp + 58);

    auto t_yzzzz_z = cbuffer.data(idx_hp + 59);

    auto t_zzzzz_x = cbuffer.data(idx_hp + 60);

    auto t_zzzzz_y = cbuffer.data(idx_hp + 61);

    auto t_zzzzz_z = cbuffer.data(idx_hp + 62);

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

    auto t_yyyyyy_y = cbuffer.data(idx_ip + 64);

    auto t_yyyyyy_z = cbuffer.data(idx_ip + 65);

    auto t_yyyyyz_y = cbuffer.data(idx_ip + 67);

    auto t_yyyyyz_z = cbuffer.data(idx_ip + 68);

    auto t_yyyyzz_y = cbuffer.data(idx_ip + 70);

    auto t_yyyyzz_z = cbuffer.data(idx_ip + 71);

    auto t_yyyzzz_y = cbuffer.data(idx_ip + 73);

    auto t_yyyzzz_z = cbuffer.data(idx_ip + 74);

    auto t_yyzzzz_y = cbuffer.data(idx_ip + 76);

    auto t_yyzzzz_z = cbuffer.data(idx_ip + 77);

    auto t_yzzzzz_y = cbuffer.data(idx_ip + 79);

    auto t_yzzzzz_z = cbuffer.data(idx_ip + 80);

    auto t_zzzzzz_z = cbuffer.data(idx_ip + 83);

    // Set up components of targeted buffer : HD

    auto t_xxxxx_xx = cbuffer.data(idx_hd);

    auto t_xxxxx_xy = cbuffer.data(idx_hd + 1);

    auto t_xxxxx_xz = cbuffer.data(idx_hd + 2);

    auto t_xxxxx_yy = cbuffer.data(idx_hd + 3);

    auto t_xxxxx_yz = cbuffer.data(idx_hd + 4);

    auto t_xxxxx_zz = cbuffer.data(idx_hd + 5);

    auto t_xxxxy_xx = cbuffer.data(idx_hd + 6);

    auto t_xxxxy_xy = cbuffer.data(idx_hd + 7);

    auto t_xxxxy_xz = cbuffer.data(idx_hd + 8);

    auto t_xxxxy_yy = cbuffer.data(idx_hd + 9);

    auto t_xxxxy_yz = cbuffer.data(idx_hd + 10);

    auto t_xxxxy_zz = cbuffer.data(idx_hd + 11);

    auto t_xxxxz_xx = cbuffer.data(idx_hd + 12);

    auto t_xxxxz_xy = cbuffer.data(idx_hd + 13);

    auto t_xxxxz_xz = cbuffer.data(idx_hd + 14);

    auto t_xxxxz_yy = cbuffer.data(idx_hd + 15);

    auto t_xxxxz_yz = cbuffer.data(idx_hd + 16);

    auto t_xxxxz_zz = cbuffer.data(idx_hd + 17);

    auto t_xxxyy_xx = cbuffer.data(idx_hd + 18);

    auto t_xxxyy_xy = cbuffer.data(idx_hd + 19);

    auto t_xxxyy_xz = cbuffer.data(idx_hd + 20);

    auto t_xxxyy_yy = cbuffer.data(idx_hd + 21);

    auto t_xxxyy_yz = cbuffer.data(idx_hd + 22);

    auto t_xxxyy_zz = cbuffer.data(idx_hd + 23);

    auto t_xxxyz_xx = cbuffer.data(idx_hd + 24);

    auto t_xxxyz_xy = cbuffer.data(idx_hd + 25);

    auto t_xxxyz_xz = cbuffer.data(idx_hd + 26);

    auto t_xxxyz_yy = cbuffer.data(idx_hd + 27);

    auto t_xxxyz_yz = cbuffer.data(idx_hd + 28);

    auto t_xxxyz_zz = cbuffer.data(idx_hd + 29);

    auto t_xxxzz_xx = cbuffer.data(idx_hd + 30);

    auto t_xxxzz_xy = cbuffer.data(idx_hd + 31);

    auto t_xxxzz_xz = cbuffer.data(idx_hd + 32);

    auto t_xxxzz_yy = cbuffer.data(idx_hd + 33);

    auto t_xxxzz_yz = cbuffer.data(idx_hd + 34);

    auto t_xxxzz_zz = cbuffer.data(idx_hd + 35);

    auto t_xxyyy_xx = cbuffer.data(idx_hd + 36);

    auto t_xxyyy_xy = cbuffer.data(idx_hd + 37);

    auto t_xxyyy_xz = cbuffer.data(idx_hd + 38);

    auto t_xxyyy_yy = cbuffer.data(idx_hd + 39);

    auto t_xxyyy_yz = cbuffer.data(idx_hd + 40);

    auto t_xxyyy_zz = cbuffer.data(idx_hd + 41);

    auto t_xxyyz_xx = cbuffer.data(idx_hd + 42);

    auto t_xxyyz_xy = cbuffer.data(idx_hd + 43);

    auto t_xxyyz_xz = cbuffer.data(idx_hd + 44);

    auto t_xxyyz_yy = cbuffer.data(idx_hd + 45);

    auto t_xxyyz_yz = cbuffer.data(idx_hd + 46);

    auto t_xxyyz_zz = cbuffer.data(idx_hd + 47);

    auto t_xxyzz_xx = cbuffer.data(idx_hd + 48);

    auto t_xxyzz_xy = cbuffer.data(idx_hd + 49);

    auto t_xxyzz_xz = cbuffer.data(idx_hd + 50);

    auto t_xxyzz_yy = cbuffer.data(idx_hd + 51);

    auto t_xxyzz_yz = cbuffer.data(idx_hd + 52);

    auto t_xxyzz_zz = cbuffer.data(idx_hd + 53);

    auto t_xxzzz_xx = cbuffer.data(idx_hd + 54);

    auto t_xxzzz_xy = cbuffer.data(idx_hd + 55);

    auto t_xxzzz_xz = cbuffer.data(idx_hd + 56);

    auto t_xxzzz_yy = cbuffer.data(idx_hd + 57);

    auto t_xxzzz_yz = cbuffer.data(idx_hd + 58);

    auto t_xxzzz_zz = cbuffer.data(idx_hd + 59);

    auto t_xyyyy_xx = cbuffer.data(idx_hd + 60);

    auto t_xyyyy_xy = cbuffer.data(idx_hd + 61);

    auto t_xyyyy_xz = cbuffer.data(idx_hd + 62);

    auto t_xyyyy_yy = cbuffer.data(idx_hd + 63);

    auto t_xyyyy_yz = cbuffer.data(idx_hd + 64);

    auto t_xyyyy_zz = cbuffer.data(idx_hd + 65);

    auto t_xyyyz_xx = cbuffer.data(idx_hd + 66);

    auto t_xyyyz_xy = cbuffer.data(idx_hd + 67);

    auto t_xyyyz_xz = cbuffer.data(idx_hd + 68);

    auto t_xyyyz_yy = cbuffer.data(idx_hd + 69);

    auto t_xyyyz_yz = cbuffer.data(idx_hd + 70);

    auto t_xyyyz_zz = cbuffer.data(idx_hd + 71);

    auto t_xyyzz_xx = cbuffer.data(idx_hd + 72);

    auto t_xyyzz_xy = cbuffer.data(idx_hd + 73);

    auto t_xyyzz_xz = cbuffer.data(idx_hd + 74);

    auto t_xyyzz_yy = cbuffer.data(idx_hd + 75);

    auto t_xyyzz_yz = cbuffer.data(idx_hd + 76);

    auto t_xyyzz_zz = cbuffer.data(idx_hd + 77);

    auto t_xyzzz_xx = cbuffer.data(idx_hd + 78);

    auto t_xyzzz_xy = cbuffer.data(idx_hd + 79);

    auto t_xyzzz_xz = cbuffer.data(idx_hd + 80);

    auto t_xyzzz_yy = cbuffer.data(idx_hd + 81);

    auto t_xyzzz_yz = cbuffer.data(idx_hd + 82);

    auto t_xyzzz_zz = cbuffer.data(idx_hd + 83);

    auto t_xzzzz_xx = cbuffer.data(idx_hd + 84);

    auto t_xzzzz_xy = cbuffer.data(idx_hd + 85);

    auto t_xzzzz_xz = cbuffer.data(idx_hd + 86);

    auto t_xzzzz_yy = cbuffer.data(idx_hd + 87);

    auto t_xzzzz_yz = cbuffer.data(idx_hd + 88);

    auto t_xzzzz_zz = cbuffer.data(idx_hd + 89);

    auto t_yyyyy_xx = cbuffer.data(idx_hd + 90);

    auto t_yyyyy_xy = cbuffer.data(idx_hd + 91);

    auto t_yyyyy_xz = cbuffer.data(idx_hd + 92);

    auto t_yyyyy_yy = cbuffer.data(idx_hd + 93);

    auto t_yyyyy_yz = cbuffer.data(idx_hd + 94);

    auto t_yyyyy_zz = cbuffer.data(idx_hd + 95);

    auto t_yyyyz_xx = cbuffer.data(idx_hd + 96);

    auto t_yyyyz_xy = cbuffer.data(idx_hd + 97);

    auto t_yyyyz_xz = cbuffer.data(idx_hd + 98);

    auto t_yyyyz_yy = cbuffer.data(idx_hd + 99);

    auto t_yyyyz_yz = cbuffer.data(idx_hd + 100);

    auto t_yyyyz_zz = cbuffer.data(idx_hd + 101);

    auto t_yyyzz_xx = cbuffer.data(idx_hd + 102);

    auto t_yyyzz_xy = cbuffer.data(idx_hd + 103);

    auto t_yyyzz_xz = cbuffer.data(idx_hd + 104);

    auto t_yyyzz_yy = cbuffer.data(idx_hd + 105);

    auto t_yyyzz_yz = cbuffer.data(idx_hd + 106);

    auto t_yyyzz_zz = cbuffer.data(idx_hd + 107);

    auto t_yyzzz_xx = cbuffer.data(idx_hd + 108);

    auto t_yyzzz_xy = cbuffer.data(idx_hd + 109);

    auto t_yyzzz_xz = cbuffer.data(idx_hd + 110);

    auto t_yyzzz_yy = cbuffer.data(idx_hd + 111);

    auto t_yyzzz_yz = cbuffer.data(idx_hd + 112);

    auto t_yyzzz_zz = cbuffer.data(idx_hd + 113);

    auto t_yzzzz_xx = cbuffer.data(idx_hd + 114);

    auto t_yzzzz_xy = cbuffer.data(idx_hd + 115);

    auto t_yzzzz_xz = cbuffer.data(idx_hd + 116);

    auto t_yzzzz_yy = cbuffer.data(idx_hd + 117);

    auto t_yzzzz_yz = cbuffer.data(idx_hd + 118);

    auto t_yzzzz_zz = cbuffer.data(idx_hd + 119);

    auto t_zzzzz_xx = cbuffer.data(idx_hd + 120);

    auto t_zzzzz_xy = cbuffer.data(idx_hd + 121);

    auto t_zzzzz_xz = cbuffer.data(idx_hd + 122);

    auto t_zzzzz_yy = cbuffer.data(idx_hd + 123);

    auto t_zzzzz_yz = cbuffer.data(idx_hd + 124);

    auto t_zzzzz_zz = cbuffer.data(idx_hd + 125);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxxx_x, t_xxxxx_xx, t_xxxxx_xy, t_xxxxx_xz, t_xxxxx_y, t_xxxxx_yy, t_xxxxx_yz, t_xxxxx_z, t_xxxxx_zz, t_xxxxxx_x, t_xxxxxx_y, t_xxxxxx_z, t_xxxxxy_x, t_xxxxxy_y, t_xxxxxy_z, t_xxxxxz_x, t_xxxxxz_y, t_xxxxxz_z, t_xxxxy_x, t_xxxxy_xx, t_xxxxy_xy, t_xxxxy_xz, t_xxxxy_y, t_xxxxy_yy, t_xxxxy_yz, t_xxxxy_z, t_xxxxy_zz, t_xxxxyy_x, t_xxxxyy_y, t_xxxxyy_z, t_xxxxyz_x, t_xxxxyz_y, t_xxxxyz_z, t_xxxxz_x, t_xxxxz_xx, t_xxxxz_xy, t_xxxxz_xz, t_xxxxz_y, t_xxxxz_yy, t_xxxxz_yz, t_xxxxz_z, t_xxxxz_zz, t_xxxxzz_x, t_xxxxzz_y, t_xxxxzz_z, t_xxxyy_x, t_xxxyy_xx, t_xxxyy_xy, t_xxxyy_xz, t_xxxyy_y, t_xxxyy_yy, t_xxxyy_yz, t_xxxyy_z, t_xxxyy_zz, t_xxxyyy_x, t_xxxyyy_y, t_xxxyyy_z, t_xxxyyz_x, t_xxxyyz_y, t_xxxyyz_z, t_xxxyz_x, t_xxxyz_xx, t_xxxyz_xy, t_xxxyz_xz, t_xxxyz_y, t_xxxyz_yy, t_xxxyz_yz, t_xxxyz_z, t_xxxyz_zz, t_xxxyzz_x, t_xxxyzz_y, t_xxxyzz_z, t_xxxzz_x, t_xxxzz_xx, t_xxxzz_xy, t_xxxzz_xz, t_xxxzz_y, t_xxxzz_yy, t_xxxzz_yz, t_xxxzz_z, t_xxxzz_zz, t_xxxzzz_x, t_xxxzzz_y, t_xxxzzz_z, t_xxyyy_x, t_xxyyy_xx, t_xxyyy_xy, t_xxyyy_xz, t_xxyyy_y, t_xxyyy_yy, t_xxyyy_yz, t_xxyyy_z, t_xxyyy_zz, t_xxyyyy_x, t_xxyyyy_y, t_xxyyyy_z, t_xxyyyz_x, t_xxyyyz_y, t_xxyyyz_z, t_xxyyz_x, t_xxyyz_xx, t_xxyyz_xy, t_xxyyz_xz, t_xxyyz_y, t_xxyyz_yy, t_xxyyz_yz, t_xxyyz_z, t_xxyyz_zz, t_xxyyzz_x, t_xxyyzz_y, t_xxyyzz_z, t_xxyzz_x, t_xxyzz_xx, t_xxyzz_xy, t_xxyzz_xz, t_xxyzz_y, t_xxyzz_yy, t_xxyzz_yz, t_xxyzz_z, t_xxyzz_zz, t_xxyzzz_x, t_xxyzzz_y, t_xxyzzz_z, t_xxzzz_x, t_xxzzz_xx, t_xxzzz_xy, t_xxzzz_xz, t_xxzzz_y, t_xxzzz_yy, t_xxzzz_yz, t_xxzzz_z, t_xxzzz_zz, t_xxzzzz_x, t_xxzzzz_y, t_xxzzzz_z, t_xyyyy_x, t_xyyyy_xx, t_xyyyy_xy, t_xyyyy_xz, t_xyyyy_y, t_xyyyy_yy, t_xyyyy_yz, t_xyyyy_z, t_xyyyy_zz, t_xyyyyy_x, t_xyyyyy_y, t_xyyyyy_z, t_xyyyyz_x, t_xyyyyz_y, t_xyyyyz_z, t_xyyyz_x, t_xyyyz_xx, t_xyyyz_xy, t_xyyyz_xz, t_xyyyz_y, t_xyyyz_yy, t_xyyyz_yz, t_xyyyz_z, t_xyyyz_zz, t_xyyyzz_x, t_xyyyzz_y, t_xyyyzz_z, t_xyyzz_x, t_xyyzz_xx, t_xyyzz_xy, t_xyyzz_xz, t_xyyzz_y, t_xyyzz_yy, t_xyyzz_yz, t_xyyzz_z, t_xyyzz_zz, t_xyyzzz_x, t_xyyzzz_y, t_xyyzzz_z, t_xyzzz_x, t_xyzzz_xx, t_xyzzz_xy, t_xyzzz_xz, t_xyzzz_y, t_xyzzz_yy, t_xyzzz_yz, t_xyzzz_z, t_xyzzz_zz, t_xyzzzz_x, t_xyzzzz_y, t_xyzzzz_z, t_xzzzz_x, t_xzzzz_xx, t_xzzzz_xy, t_xzzzz_xz, t_xzzzz_y, t_xzzzz_yy, t_xzzzz_yz, t_xzzzz_z, t_xzzzz_zz, t_xzzzzz_x, t_xzzzzz_y, t_xzzzzz_z, t_yyyyy_x, t_yyyyy_xx, t_yyyyy_xy, t_yyyyy_xz, t_yyyyy_y, t_yyyyy_yy, t_yyyyy_yz, t_yyyyy_z, t_yyyyy_zz, t_yyyyyy_y, t_yyyyyy_z, t_yyyyyz_y, t_yyyyyz_z, t_yyyyz_x, t_yyyyz_xx, t_yyyyz_xy, t_yyyyz_xz, t_yyyyz_y, t_yyyyz_yy, t_yyyyz_yz, t_yyyyz_z, t_yyyyz_zz, t_yyyyzz_y, t_yyyyzz_z, t_yyyzz_x, t_yyyzz_xx, t_yyyzz_xy, t_yyyzz_xz, t_yyyzz_y, t_yyyzz_yy, t_yyyzz_yz, t_yyyzz_z, t_yyyzz_zz, t_yyyzzz_y, t_yyyzzz_z, t_yyzzz_x, t_yyzzz_xx, t_yyzzz_xy, t_yyzzz_xz, t_yyzzz_y, t_yyzzz_yy, t_yyzzz_yz, t_yyzzz_z, t_yyzzz_zz, t_yyzzzz_y, t_yyzzzz_z, t_yzzzz_x, t_yzzzz_xx, t_yzzzz_xy, t_yzzzz_xz, t_yzzzz_y, t_yzzzz_yy, t_yzzzz_yz, t_yzzzz_z, t_yzzzz_zz, t_yzzzzz_y, t_yzzzzz_z, t_zzzzz_x, t_zzzzz_xx, t_zzzzz_xy, t_zzzzz_xz, t_zzzzz_y, t_zzzzz_yy, t_zzzzz_yz, t_zzzzz_z, t_zzzzz_zz, t_zzzzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxxx_xx[i] = t_xxxxx_x[i] * ab_x[i] + t_xxxxxx_x[i];

        t_xxxxx_xy[i] = t_xxxxx_y[i] * ab_x[i] + t_xxxxxx_y[i];

        t_xxxxx_xz[i] = t_xxxxx_z[i] * ab_x[i] + t_xxxxxx_z[i];

        t_xxxxx_yy[i] = t_xxxxx_y[i] * ab_y[i] + t_xxxxxy_y[i];

        t_xxxxx_yz[i] = t_xxxxx_z[i] * ab_y[i] + t_xxxxxy_z[i];

        t_xxxxx_zz[i] = t_xxxxx_z[i] * ab_z[i] + t_xxxxxz_z[i];

        t_xxxxy_xx[i] = t_xxxxy_x[i] * ab_x[i] + t_xxxxxy_x[i];

        t_xxxxy_xy[i] = t_xxxxy_y[i] * ab_x[i] + t_xxxxxy_y[i];

        t_xxxxy_xz[i] = t_xxxxy_z[i] * ab_x[i] + t_xxxxxy_z[i];

        t_xxxxy_yy[i] = t_xxxxy_y[i] * ab_y[i] + t_xxxxyy_y[i];

        t_xxxxy_yz[i] = t_xxxxy_z[i] * ab_y[i] + t_xxxxyy_z[i];

        t_xxxxy_zz[i] = t_xxxxy_z[i] * ab_z[i] + t_xxxxyz_z[i];

        t_xxxxz_xx[i] = t_xxxxz_x[i] * ab_x[i] + t_xxxxxz_x[i];

        t_xxxxz_xy[i] = t_xxxxz_y[i] * ab_x[i] + t_xxxxxz_y[i];

        t_xxxxz_xz[i] = t_xxxxz_z[i] * ab_x[i] + t_xxxxxz_z[i];

        t_xxxxz_yy[i] = t_xxxxz_y[i] * ab_y[i] + t_xxxxyz_y[i];

        t_xxxxz_yz[i] = t_xxxxz_z[i] * ab_y[i] + t_xxxxyz_z[i];

        t_xxxxz_zz[i] = t_xxxxz_z[i] * ab_z[i] + t_xxxxzz_z[i];

        t_xxxyy_xx[i] = t_xxxyy_x[i] * ab_x[i] + t_xxxxyy_x[i];

        t_xxxyy_xy[i] = t_xxxyy_y[i] * ab_x[i] + t_xxxxyy_y[i];

        t_xxxyy_xz[i] = t_xxxyy_z[i] * ab_x[i] + t_xxxxyy_z[i];

        t_xxxyy_yy[i] = t_xxxyy_y[i] * ab_y[i] + t_xxxyyy_y[i];

        t_xxxyy_yz[i] = t_xxxyy_z[i] * ab_y[i] + t_xxxyyy_z[i];

        t_xxxyy_zz[i] = t_xxxyy_z[i] * ab_z[i] + t_xxxyyz_z[i];

        t_xxxyz_xx[i] = t_xxxyz_x[i] * ab_x[i] + t_xxxxyz_x[i];

        t_xxxyz_xy[i] = t_xxxyz_y[i] * ab_x[i] + t_xxxxyz_y[i];

        t_xxxyz_xz[i] = t_xxxyz_z[i] * ab_x[i] + t_xxxxyz_z[i];

        t_xxxyz_yy[i] = t_xxxyz_y[i] * ab_y[i] + t_xxxyyz_y[i];

        t_xxxyz_yz[i] = t_xxxyz_z[i] * ab_y[i] + t_xxxyyz_z[i];

        t_xxxyz_zz[i] = t_xxxyz_z[i] * ab_z[i] + t_xxxyzz_z[i];

        t_xxxzz_xx[i] = t_xxxzz_x[i] * ab_x[i] + t_xxxxzz_x[i];

        t_xxxzz_xy[i] = t_xxxzz_y[i] * ab_x[i] + t_xxxxzz_y[i];

        t_xxxzz_xz[i] = t_xxxzz_z[i] * ab_x[i] + t_xxxxzz_z[i];

        t_xxxzz_yy[i] = t_xxxzz_y[i] * ab_y[i] + t_xxxyzz_y[i];

        t_xxxzz_yz[i] = t_xxxzz_z[i] * ab_y[i] + t_xxxyzz_z[i];

        t_xxxzz_zz[i] = t_xxxzz_z[i] * ab_z[i] + t_xxxzzz_z[i];

        t_xxyyy_xx[i] = t_xxyyy_x[i] * ab_x[i] + t_xxxyyy_x[i];

        t_xxyyy_xy[i] = t_xxyyy_y[i] * ab_x[i] + t_xxxyyy_y[i];

        t_xxyyy_xz[i] = t_xxyyy_z[i] * ab_x[i] + t_xxxyyy_z[i];

        t_xxyyy_yy[i] = t_xxyyy_y[i] * ab_y[i] + t_xxyyyy_y[i];

        t_xxyyy_yz[i] = t_xxyyy_z[i] * ab_y[i] + t_xxyyyy_z[i];

        t_xxyyy_zz[i] = t_xxyyy_z[i] * ab_z[i] + t_xxyyyz_z[i];

        t_xxyyz_xx[i] = t_xxyyz_x[i] * ab_x[i] + t_xxxyyz_x[i];

        t_xxyyz_xy[i] = t_xxyyz_y[i] * ab_x[i] + t_xxxyyz_y[i];

        t_xxyyz_xz[i] = t_xxyyz_z[i] * ab_x[i] + t_xxxyyz_z[i];

        t_xxyyz_yy[i] = t_xxyyz_y[i] * ab_y[i] + t_xxyyyz_y[i];

        t_xxyyz_yz[i] = t_xxyyz_z[i] * ab_y[i] + t_xxyyyz_z[i];

        t_xxyyz_zz[i] = t_xxyyz_z[i] * ab_z[i] + t_xxyyzz_z[i];

        t_xxyzz_xx[i] = t_xxyzz_x[i] * ab_x[i] + t_xxxyzz_x[i];

        t_xxyzz_xy[i] = t_xxyzz_y[i] * ab_x[i] + t_xxxyzz_y[i];

        t_xxyzz_xz[i] = t_xxyzz_z[i] * ab_x[i] + t_xxxyzz_z[i];

        t_xxyzz_yy[i] = t_xxyzz_y[i] * ab_y[i] + t_xxyyzz_y[i];

        t_xxyzz_yz[i] = t_xxyzz_z[i] * ab_y[i] + t_xxyyzz_z[i];

        t_xxyzz_zz[i] = t_xxyzz_z[i] * ab_z[i] + t_xxyzzz_z[i];

        t_xxzzz_xx[i] = t_xxzzz_x[i] * ab_x[i] + t_xxxzzz_x[i];

        t_xxzzz_xy[i] = t_xxzzz_y[i] * ab_x[i] + t_xxxzzz_y[i];

        t_xxzzz_xz[i] = t_xxzzz_z[i] * ab_x[i] + t_xxxzzz_z[i];

        t_xxzzz_yy[i] = t_xxzzz_y[i] * ab_y[i] + t_xxyzzz_y[i];

        t_xxzzz_yz[i] = t_xxzzz_z[i] * ab_y[i] + t_xxyzzz_z[i];

        t_xxzzz_zz[i] = t_xxzzz_z[i] * ab_z[i] + t_xxzzzz_z[i];

        t_xyyyy_xx[i] = t_xyyyy_x[i] * ab_x[i] + t_xxyyyy_x[i];

        t_xyyyy_xy[i] = t_xyyyy_y[i] * ab_x[i] + t_xxyyyy_y[i];

        t_xyyyy_xz[i] = t_xyyyy_z[i] * ab_x[i] + t_xxyyyy_z[i];

        t_xyyyy_yy[i] = t_xyyyy_y[i] * ab_y[i] + t_xyyyyy_y[i];

        t_xyyyy_yz[i] = t_xyyyy_z[i] * ab_y[i] + t_xyyyyy_z[i];

        t_xyyyy_zz[i] = t_xyyyy_z[i] * ab_z[i] + t_xyyyyz_z[i];

        t_xyyyz_xx[i] = t_xyyyz_x[i] * ab_x[i] + t_xxyyyz_x[i];

        t_xyyyz_xy[i] = t_xyyyz_y[i] * ab_x[i] + t_xxyyyz_y[i];

        t_xyyyz_xz[i] = t_xyyyz_z[i] * ab_x[i] + t_xxyyyz_z[i];

        t_xyyyz_yy[i] = t_xyyyz_y[i] * ab_y[i] + t_xyyyyz_y[i];

        t_xyyyz_yz[i] = t_xyyyz_z[i] * ab_y[i] + t_xyyyyz_z[i];

        t_xyyyz_zz[i] = t_xyyyz_z[i] * ab_z[i] + t_xyyyzz_z[i];

        t_xyyzz_xx[i] = t_xyyzz_x[i] * ab_x[i] + t_xxyyzz_x[i];

        t_xyyzz_xy[i] = t_xyyzz_y[i] * ab_x[i] + t_xxyyzz_y[i];

        t_xyyzz_xz[i] = t_xyyzz_z[i] * ab_x[i] + t_xxyyzz_z[i];

        t_xyyzz_yy[i] = t_xyyzz_y[i] * ab_y[i] + t_xyyyzz_y[i];

        t_xyyzz_yz[i] = t_xyyzz_z[i] * ab_y[i] + t_xyyyzz_z[i];

        t_xyyzz_zz[i] = t_xyyzz_z[i] * ab_z[i] + t_xyyzzz_z[i];

        t_xyzzz_xx[i] = t_xyzzz_x[i] * ab_x[i] + t_xxyzzz_x[i];

        t_xyzzz_xy[i] = t_xyzzz_y[i] * ab_x[i] + t_xxyzzz_y[i];

        t_xyzzz_xz[i] = t_xyzzz_z[i] * ab_x[i] + t_xxyzzz_z[i];

        t_xyzzz_yy[i] = t_xyzzz_y[i] * ab_y[i] + t_xyyzzz_y[i];

        t_xyzzz_yz[i] = t_xyzzz_z[i] * ab_y[i] + t_xyyzzz_z[i];

        t_xyzzz_zz[i] = t_xyzzz_z[i] * ab_z[i] + t_xyzzzz_z[i];

        t_xzzzz_xx[i] = t_xzzzz_x[i] * ab_x[i] + t_xxzzzz_x[i];

        t_xzzzz_xy[i] = t_xzzzz_y[i] * ab_x[i] + t_xxzzzz_y[i];

        t_xzzzz_xz[i] = t_xzzzz_z[i] * ab_x[i] + t_xxzzzz_z[i];

        t_xzzzz_yy[i] = t_xzzzz_y[i] * ab_y[i] + t_xyzzzz_y[i];

        t_xzzzz_yz[i] = t_xzzzz_z[i] * ab_y[i] + t_xyzzzz_z[i];

        t_xzzzz_zz[i] = t_xzzzz_z[i] * ab_z[i] + t_xzzzzz_z[i];

        t_yyyyy_xx[i] = t_yyyyy_x[i] * ab_x[i] + t_xyyyyy_x[i];

        t_yyyyy_xy[i] = t_yyyyy_y[i] * ab_x[i] + t_xyyyyy_y[i];

        t_yyyyy_xz[i] = t_yyyyy_z[i] * ab_x[i] + t_xyyyyy_z[i];

        t_yyyyy_yy[i] = t_yyyyy_y[i] * ab_y[i] + t_yyyyyy_y[i];

        t_yyyyy_yz[i] = t_yyyyy_z[i] * ab_y[i] + t_yyyyyy_z[i];

        t_yyyyy_zz[i] = t_yyyyy_z[i] * ab_z[i] + t_yyyyyz_z[i];

        t_yyyyz_xx[i] = t_yyyyz_x[i] * ab_x[i] + t_xyyyyz_x[i];

        t_yyyyz_xy[i] = t_yyyyz_y[i] * ab_x[i] + t_xyyyyz_y[i];

        t_yyyyz_xz[i] = t_yyyyz_z[i] * ab_x[i] + t_xyyyyz_z[i];

        t_yyyyz_yy[i] = t_yyyyz_y[i] * ab_y[i] + t_yyyyyz_y[i];

        t_yyyyz_yz[i] = t_yyyyz_z[i] * ab_y[i] + t_yyyyyz_z[i];

        t_yyyyz_zz[i] = t_yyyyz_z[i] * ab_z[i] + t_yyyyzz_z[i];

        t_yyyzz_xx[i] = t_yyyzz_x[i] * ab_x[i] + t_xyyyzz_x[i];

        t_yyyzz_xy[i] = t_yyyzz_y[i] * ab_x[i] + t_xyyyzz_y[i];

        t_yyyzz_xz[i] = t_yyyzz_z[i] * ab_x[i] + t_xyyyzz_z[i];

        t_yyyzz_yy[i] = t_yyyzz_y[i] * ab_y[i] + t_yyyyzz_y[i];

        t_yyyzz_yz[i] = t_yyyzz_z[i] * ab_y[i] + t_yyyyzz_z[i];

        t_yyyzz_zz[i] = t_yyyzz_z[i] * ab_z[i] + t_yyyzzz_z[i];

        t_yyzzz_xx[i] = t_yyzzz_x[i] * ab_x[i] + t_xyyzzz_x[i];

        t_yyzzz_xy[i] = t_yyzzz_y[i] * ab_x[i] + t_xyyzzz_y[i];

        t_yyzzz_xz[i] = t_yyzzz_z[i] * ab_x[i] + t_xyyzzz_z[i];

        t_yyzzz_yy[i] = t_yyzzz_y[i] * ab_y[i] + t_yyyzzz_y[i];

        t_yyzzz_yz[i] = t_yyzzz_z[i] * ab_y[i] + t_yyyzzz_z[i];

        t_yyzzz_zz[i] = t_yyzzz_z[i] * ab_z[i] + t_yyzzzz_z[i];

        t_yzzzz_xx[i] = t_yzzzz_x[i] * ab_x[i] + t_xyzzzz_x[i];

        t_yzzzz_xy[i] = t_yzzzz_y[i] * ab_x[i] + t_xyzzzz_y[i];

        t_yzzzz_xz[i] = t_yzzzz_z[i] * ab_x[i] + t_xyzzzz_z[i];

        t_yzzzz_yy[i] = t_yzzzz_y[i] * ab_y[i] + t_yyzzzz_y[i];

        t_yzzzz_yz[i] = t_yzzzz_z[i] * ab_y[i] + t_yyzzzz_z[i];

        t_yzzzz_zz[i] = t_yzzzz_z[i] * ab_z[i] + t_yzzzzz_z[i];

        t_zzzzz_xx[i] = t_zzzzz_x[i] * ab_x[i] + t_xzzzzz_x[i];

        t_zzzzz_xy[i] = t_zzzzz_y[i] * ab_x[i] + t_xzzzzz_y[i];

        t_zzzzz_xz[i] = t_zzzzz_z[i] * ab_x[i] + t_xzzzzz_z[i];

        t_zzzzz_yy[i] = t_zzzzz_y[i] * ab_y[i] + t_yzzzzz_y[i];

        t_zzzzz_yz[i] = t_zzzzz_z[i] * ab_y[i] + t_yzzzzz_z[i];

        t_zzzzz_zz[i] = t_zzzzz_z[i] * ab_z[i] + t_zzzzzz_z[i];
    }
}

} // t2chrr namespace

