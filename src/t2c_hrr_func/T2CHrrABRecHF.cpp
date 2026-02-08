#include "T2CHrrABRecHF.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_hf(CSimdArray<double>& cbuffer, 
            const size_t idx_hf,
            const size_t idx_hd,
            const size_t idx_id,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : HD

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

    // Set up components of auxiliary buffer : ID

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

    auto t_yyyyyy_yy = cbuffer.data(idx_id + 129);

    auto t_yyyyyy_yz = cbuffer.data(idx_id + 130);

    auto t_yyyyyy_zz = cbuffer.data(idx_id + 131);

    auto t_yyyyyz_yy = cbuffer.data(idx_id + 135);

    auto t_yyyyyz_yz = cbuffer.data(idx_id + 136);

    auto t_yyyyyz_zz = cbuffer.data(idx_id + 137);

    auto t_yyyyzz_yy = cbuffer.data(idx_id + 141);

    auto t_yyyyzz_yz = cbuffer.data(idx_id + 142);

    auto t_yyyyzz_zz = cbuffer.data(idx_id + 143);

    auto t_yyyzzz_yy = cbuffer.data(idx_id + 147);

    auto t_yyyzzz_yz = cbuffer.data(idx_id + 148);

    auto t_yyyzzz_zz = cbuffer.data(idx_id + 149);

    auto t_yyzzzz_yy = cbuffer.data(idx_id + 153);

    auto t_yyzzzz_yz = cbuffer.data(idx_id + 154);

    auto t_yyzzzz_zz = cbuffer.data(idx_id + 155);

    auto t_yzzzzz_yy = cbuffer.data(idx_id + 159);

    auto t_yzzzzz_yz = cbuffer.data(idx_id + 160);

    auto t_yzzzzz_zz = cbuffer.data(idx_id + 161);

    auto t_zzzzzz_zz = cbuffer.data(idx_id + 167);

    // Set up components of targeted buffer : HF

    auto t_xxxxx_xxx = cbuffer.data(idx_hf);

    auto t_xxxxx_xxy = cbuffer.data(idx_hf + 1);

    auto t_xxxxx_xxz = cbuffer.data(idx_hf + 2);

    auto t_xxxxx_xyy = cbuffer.data(idx_hf + 3);

    auto t_xxxxx_xyz = cbuffer.data(idx_hf + 4);

    auto t_xxxxx_xzz = cbuffer.data(idx_hf + 5);

    auto t_xxxxx_yyy = cbuffer.data(idx_hf + 6);

    auto t_xxxxx_yyz = cbuffer.data(idx_hf + 7);

    auto t_xxxxx_yzz = cbuffer.data(idx_hf + 8);

    auto t_xxxxx_zzz = cbuffer.data(idx_hf + 9);

    auto t_xxxxy_xxx = cbuffer.data(idx_hf + 10);

    auto t_xxxxy_xxy = cbuffer.data(idx_hf + 11);

    auto t_xxxxy_xxz = cbuffer.data(idx_hf + 12);

    auto t_xxxxy_xyy = cbuffer.data(idx_hf + 13);

    auto t_xxxxy_xyz = cbuffer.data(idx_hf + 14);

    auto t_xxxxy_xzz = cbuffer.data(idx_hf + 15);

    auto t_xxxxy_yyy = cbuffer.data(idx_hf + 16);

    auto t_xxxxy_yyz = cbuffer.data(idx_hf + 17);

    auto t_xxxxy_yzz = cbuffer.data(idx_hf + 18);

    auto t_xxxxy_zzz = cbuffer.data(idx_hf + 19);

    auto t_xxxxz_xxx = cbuffer.data(idx_hf + 20);

    auto t_xxxxz_xxy = cbuffer.data(idx_hf + 21);

    auto t_xxxxz_xxz = cbuffer.data(idx_hf + 22);

    auto t_xxxxz_xyy = cbuffer.data(idx_hf + 23);

    auto t_xxxxz_xyz = cbuffer.data(idx_hf + 24);

    auto t_xxxxz_xzz = cbuffer.data(idx_hf + 25);

    auto t_xxxxz_yyy = cbuffer.data(idx_hf + 26);

    auto t_xxxxz_yyz = cbuffer.data(idx_hf + 27);

    auto t_xxxxz_yzz = cbuffer.data(idx_hf + 28);

    auto t_xxxxz_zzz = cbuffer.data(idx_hf + 29);

    auto t_xxxyy_xxx = cbuffer.data(idx_hf + 30);

    auto t_xxxyy_xxy = cbuffer.data(idx_hf + 31);

    auto t_xxxyy_xxz = cbuffer.data(idx_hf + 32);

    auto t_xxxyy_xyy = cbuffer.data(idx_hf + 33);

    auto t_xxxyy_xyz = cbuffer.data(idx_hf + 34);

    auto t_xxxyy_xzz = cbuffer.data(idx_hf + 35);

    auto t_xxxyy_yyy = cbuffer.data(idx_hf + 36);

    auto t_xxxyy_yyz = cbuffer.data(idx_hf + 37);

    auto t_xxxyy_yzz = cbuffer.data(idx_hf + 38);

    auto t_xxxyy_zzz = cbuffer.data(idx_hf + 39);

    auto t_xxxyz_xxx = cbuffer.data(idx_hf + 40);

    auto t_xxxyz_xxy = cbuffer.data(idx_hf + 41);

    auto t_xxxyz_xxz = cbuffer.data(idx_hf + 42);

    auto t_xxxyz_xyy = cbuffer.data(idx_hf + 43);

    auto t_xxxyz_xyz = cbuffer.data(idx_hf + 44);

    auto t_xxxyz_xzz = cbuffer.data(idx_hf + 45);

    auto t_xxxyz_yyy = cbuffer.data(idx_hf + 46);

    auto t_xxxyz_yyz = cbuffer.data(idx_hf + 47);

    auto t_xxxyz_yzz = cbuffer.data(idx_hf + 48);

    auto t_xxxyz_zzz = cbuffer.data(idx_hf + 49);

    auto t_xxxzz_xxx = cbuffer.data(idx_hf + 50);

    auto t_xxxzz_xxy = cbuffer.data(idx_hf + 51);

    auto t_xxxzz_xxz = cbuffer.data(idx_hf + 52);

    auto t_xxxzz_xyy = cbuffer.data(idx_hf + 53);

    auto t_xxxzz_xyz = cbuffer.data(idx_hf + 54);

    auto t_xxxzz_xzz = cbuffer.data(idx_hf + 55);

    auto t_xxxzz_yyy = cbuffer.data(idx_hf + 56);

    auto t_xxxzz_yyz = cbuffer.data(idx_hf + 57);

    auto t_xxxzz_yzz = cbuffer.data(idx_hf + 58);

    auto t_xxxzz_zzz = cbuffer.data(idx_hf + 59);

    auto t_xxyyy_xxx = cbuffer.data(idx_hf + 60);

    auto t_xxyyy_xxy = cbuffer.data(idx_hf + 61);

    auto t_xxyyy_xxz = cbuffer.data(idx_hf + 62);

    auto t_xxyyy_xyy = cbuffer.data(idx_hf + 63);

    auto t_xxyyy_xyz = cbuffer.data(idx_hf + 64);

    auto t_xxyyy_xzz = cbuffer.data(idx_hf + 65);

    auto t_xxyyy_yyy = cbuffer.data(idx_hf + 66);

    auto t_xxyyy_yyz = cbuffer.data(idx_hf + 67);

    auto t_xxyyy_yzz = cbuffer.data(idx_hf + 68);

    auto t_xxyyy_zzz = cbuffer.data(idx_hf + 69);

    auto t_xxyyz_xxx = cbuffer.data(idx_hf + 70);

    auto t_xxyyz_xxy = cbuffer.data(idx_hf + 71);

    auto t_xxyyz_xxz = cbuffer.data(idx_hf + 72);

    auto t_xxyyz_xyy = cbuffer.data(idx_hf + 73);

    auto t_xxyyz_xyz = cbuffer.data(idx_hf + 74);

    auto t_xxyyz_xzz = cbuffer.data(idx_hf + 75);

    auto t_xxyyz_yyy = cbuffer.data(idx_hf + 76);

    auto t_xxyyz_yyz = cbuffer.data(idx_hf + 77);

    auto t_xxyyz_yzz = cbuffer.data(idx_hf + 78);

    auto t_xxyyz_zzz = cbuffer.data(idx_hf + 79);

    auto t_xxyzz_xxx = cbuffer.data(idx_hf + 80);

    auto t_xxyzz_xxy = cbuffer.data(idx_hf + 81);

    auto t_xxyzz_xxz = cbuffer.data(idx_hf + 82);

    auto t_xxyzz_xyy = cbuffer.data(idx_hf + 83);

    auto t_xxyzz_xyz = cbuffer.data(idx_hf + 84);

    auto t_xxyzz_xzz = cbuffer.data(idx_hf + 85);

    auto t_xxyzz_yyy = cbuffer.data(idx_hf + 86);

    auto t_xxyzz_yyz = cbuffer.data(idx_hf + 87);

    auto t_xxyzz_yzz = cbuffer.data(idx_hf + 88);

    auto t_xxyzz_zzz = cbuffer.data(idx_hf + 89);

    auto t_xxzzz_xxx = cbuffer.data(idx_hf + 90);

    auto t_xxzzz_xxy = cbuffer.data(idx_hf + 91);

    auto t_xxzzz_xxz = cbuffer.data(idx_hf + 92);

    auto t_xxzzz_xyy = cbuffer.data(idx_hf + 93);

    auto t_xxzzz_xyz = cbuffer.data(idx_hf + 94);

    auto t_xxzzz_xzz = cbuffer.data(idx_hf + 95);

    auto t_xxzzz_yyy = cbuffer.data(idx_hf + 96);

    auto t_xxzzz_yyz = cbuffer.data(idx_hf + 97);

    auto t_xxzzz_yzz = cbuffer.data(idx_hf + 98);

    auto t_xxzzz_zzz = cbuffer.data(idx_hf + 99);

    auto t_xyyyy_xxx = cbuffer.data(idx_hf + 100);

    auto t_xyyyy_xxy = cbuffer.data(idx_hf + 101);

    auto t_xyyyy_xxz = cbuffer.data(idx_hf + 102);

    auto t_xyyyy_xyy = cbuffer.data(idx_hf + 103);

    auto t_xyyyy_xyz = cbuffer.data(idx_hf + 104);

    auto t_xyyyy_xzz = cbuffer.data(idx_hf + 105);

    auto t_xyyyy_yyy = cbuffer.data(idx_hf + 106);

    auto t_xyyyy_yyz = cbuffer.data(idx_hf + 107);

    auto t_xyyyy_yzz = cbuffer.data(idx_hf + 108);

    auto t_xyyyy_zzz = cbuffer.data(idx_hf + 109);

    auto t_xyyyz_xxx = cbuffer.data(idx_hf + 110);

    auto t_xyyyz_xxy = cbuffer.data(idx_hf + 111);

    auto t_xyyyz_xxz = cbuffer.data(idx_hf + 112);

    auto t_xyyyz_xyy = cbuffer.data(idx_hf + 113);

    auto t_xyyyz_xyz = cbuffer.data(idx_hf + 114);

    auto t_xyyyz_xzz = cbuffer.data(idx_hf + 115);

    auto t_xyyyz_yyy = cbuffer.data(idx_hf + 116);

    auto t_xyyyz_yyz = cbuffer.data(idx_hf + 117);

    auto t_xyyyz_yzz = cbuffer.data(idx_hf + 118);

    auto t_xyyyz_zzz = cbuffer.data(idx_hf + 119);

    auto t_xyyzz_xxx = cbuffer.data(idx_hf + 120);

    auto t_xyyzz_xxy = cbuffer.data(idx_hf + 121);

    auto t_xyyzz_xxz = cbuffer.data(idx_hf + 122);

    auto t_xyyzz_xyy = cbuffer.data(idx_hf + 123);

    auto t_xyyzz_xyz = cbuffer.data(idx_hf + 124);

    auto t_xyyzz_xzz = cbuffer.data(idx_hf + 125);

    auto t_xyyzz_yyy = cbuffer.data(idx_hf + 126);

    auto t_xyyzz_yyz = cbuffer.data(idx_hf + 127);

    auto t_xyyzz_yzz = cbuffer.data(idx_hf + 128);

    auto t_xyyzz_zzz = cbuffer.data(idx_hf + 129);

    auto t_xyzzz_xxx = cbuffer.data(idx_hf + 130);

    auto t_xyzzz_xxy = cbuffer.data(idx_hf + 131);

    auto t_xyzzz_xxz = cbuffer.data(idx_hf + 132);

    auto t_xyzzz_xyy = cbuffer.data(idx_hf + 133);

    auto t_xyzzz_xyz = cbuffer.data(idx_hf + 134);

    auto t_xyzzz_xzz = cbuffer.data(idx_hf + 135);

    auto t_xyzzz_yyy = cbuffer.data(idx_hf + 136);

    auto t_xyzzz_yyz = cbuffer.data(idx_hf + 137);

    auto t_xyzzz_yzz = cbuffer.data(idx_hf + 138);

    auto t_xyzzz_zzz = cbuffer.data(idx_hf + 139);

    auto t_xzzzz_xxx = cbuffer.data(idx_hf + 140);

    auto t_xzzzz_xxy = cbuffer.data(idx_hf + 141);

    auto t_xzzzz_xxz = cbuffer.data(idx_hf + 142);

    auto t_xzzzz_xyy = cbuffer.data(idx_hf + 143);

    auto t_xzzzz_xyz = cbuffer.data(idx_hf + 144);

    auto t_xzzzz_xzz = cbuffer.data(idx_hf + 145);

    auto t_xzzzz_yyy = cbuffer.data(idx_hf + 146);

    auto t_xzzzz_yyz = cbuffer.data(idx_hf + 147);

    auto t_xzzzz_yzz = cbuffer.data(idx_hf + 148);

    auto t_xzzzz_zzz = cbuffer.data(idx_hf + 149);

    auto t_yyyyy_xxx = cbuffer.data(idx_hf + 150);

    auto t_yyyyy_xxy = cbuffer.data(idx_hf + 151);

    auto t_yyyyy_xxz = cbuffer.data(idx_hf + 152);

    auto t_yyyyy_xyy = cbuffer.data(idx_hf + 153);

    auto t_yyyyy_xyz = cbuffer.data(idx_hf + 154);

    auto t_yyyyy_xzz = cbuffer.data(idx_hf + 155);

    auto t_yyyyy_yyy = cbuffer.data(idx_hf + 156);

    auto t_yyyyy_yyz = cbuffer.data(idx_hf + 157);

    auto t_yyyyy_yzz = cbuffer.data(idx_hf + 158);

    auto t_yyyyy_zzz = cbuffer.data(idx_hf + 159);

    auto t_yyyyz_xxx = cbuffer.data(idx_hf + 160);

    auto t_yyyyz_xxy = cbuffer.data(idx_hf + 161);

    auto t_yyyyz_xxz = cbuffer.data(idx_hf + 162);

    auto t_yyyyz_xyy = cbuffer.data(idx_hf + 163);

    auto t_yyyyz_xyz = cbuffer.data(idx_hf + 164);

    auto t_yyyyz_xzz = cbuffer.data(idx_hf + 165);

    auto t_yyyyz_yyy = cbuffer.data(idx_hf + 166);

    auto t_yyyyz_yyz = cbuffer.data(idx_hf + 167);

    auto t_yyyyz_yzz = cbuffer.data(idx_hf + 168);

    auto t_yyyyz_zzz = cbuffer.data(idx_hf + 169);

    auto t_yyyzz_xxx = cbuffer.data(idx_hf + 170);

    auto t_yyyzz_xxy = cbuffer.data(idx_hf + 171);

    auto t_yyyzz_xxz = cbuffer.data(idx_hf + 172);

    auto t_yyyzz_xyy = cbuffer.data(idx_hf + 173);

    auto t_yyyzz_xyz = cbuffer.data(idx_hf + 174);

    auto t_yyyzz_xzz = cbuffer.data(idx_hf + 175);

    auto t_yyyzz_yyy = cbuffer.data(idx_hf + 176);

    auto t_yyyzz_yyz = cbuffer.data(idx_hf + 177);

    auto t_yyyzz_yzz = cbuffer.data(idx_hf + 178);

    auto t_yyyzz_zzz = cbuffer.data(idx_hf + 179);

    auto t_yyzzz_xxx = cbuffer.data(idx_hf + 180);

    auto t_yyzzz_xxy = cbuffer.data(idx_hf + 181);

    auto t_yyzzz_xxz = cbuffer.data(idx_hf + 182);

    auto t_yyzzz_xyy = cbuffer.data(idx_hf + 183);

    auto t_yyzzz_xyz = cbuffer.data(idx_hf + 184);

    auto t_yyzzz_xzz = cbuffer.data(idx_hf + 185);

    auto t_yyzzz_yyy = cbuffer.data(idx_hf + 186);

    auto t_yyzzz_yyz = cbuffer.data(idx_hf + 187);

    auto t_yyzzz_yzz = cbuffer.data(idx_hf + 188);

    auto t_yyzzz_zzz = cbuffer.data(idx_hf + 189);

    auto t_yzzzz_xxx = cbuffer.data(idx_hf + 190);

    auto t_yzzzz_xxy = cbuffer.data(idx_hf + 191);

    auto t_yzzzz_xxz = cbuffer.data(idx_hf + 192);

    auto t_yzzzz_xyy = cbuffer.data(idx_hf + 193);

    auto t_yzzzz_xyz = cbuffer.data(idx_hf + 194);

    auto t_yzzzz_xzz = cbuffer.data(idx_hf + 195);

    auto t_yzzzz_yyy = cbuffer.data(idx_hf + 196);

    auto t_yzzzz_yyz = cbuffer.data(idx_hf + 197);

    auto t_yzzzz_yzz = cbuffer.data(idx_hf + 198);

    auto t_yzzzz_zzz = cbuffer.data(idx_hf + 199);

    auto t_zzzzz_xxx = cbuffer.data(idx_hf + 200);

    auto t_zzzzz_xxy = cbuffer.data(idx_hf + 201);

    auto t_zzzzz_xxz = cbuffer.data(idx_hf + 202);

    auto t_zzzzz_xyy = cbuffer.data(idx_hf + 203);

    auto t_zzzzz_xyz = cbuffer.data(idx_hf + 204);

    auto t_zzzzz_xzz = cbuffer.data(idx_hf + 205);

    auto t_zzzzz_yyy = cbuffer.data(idx_hf + 206);

    auto t_zzzzz_yyz = cbuffer.data(idx_hf + 207);

    auto t_zzzzz_yzz = cbuffer.data(idx_hf + 208);

    auto t_zzzzz_zzz = cbuffer.data(idx_hf + 209);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxxx_xx, t_xxxxx_xxx, t_xxxxx_xxy, t_xxxxx_xxz, t_xxxxx_xy, t_xxxxx_xyy, t_xxxxx_xyz, t_xxxxx_xz, t_xxxxx_xzz, t_xxxxx_yy, t_xxxxx_yyy, t_xxxxx_yyz, t_xxxxx_yz, t_xxxxx_yzz, t_xxxxx_zz, t_xxxxx_zzz, t_xxxxxx_xx, t_xxxxxx_xy, t_xxxxxx_xz, t_xxxxxx_yy, t_xxxxxx_yz, t_xxxxxx_zz, t_xxxxxy_xx, t_xxxxxy_xy, t_xxxxxy_xz, t_xxxxxy_yy, t_xxxxxy_yz, t_xxxxxy_zz, t_xxxxxz_xx, t_xxxxxz_xy, t_xxxxxz_xz, t_xxxxxz_yy, t_xxxxxz_yz, t_xxxxxz_zz, t_xxxxy_xx, t_xxxxy_xxx, t_xxxxy_xxy, t_xxxxy_xxz, t_xxxxy_xy, t_xxxxy_xyy, t_xxxxy_xyz, t_xxxxy_xz, t_xxxxy_xzz, t_xxxxy_yy, t_xxxxy_yyy, t_xxxxy_yyz, t_xxxxy_yz, t_xxxxy_yzz, t_xxxxy_zz, t_xxxxy_zzz, t_xxxxyy_xx, t_xxxxyy_xy, t_xxxxyy_xz, t_xxxxyy_yy, t_xxxxyy_yz, t_xxxxyy_zz, t_xxxxyz_xx, t_xxxxyz_xy, t_xxxxyz_xz, t_xxxxyz_yy, t_xxxxyz_yz, t_xxxxyz_zz, t_xxxxz_xx, t_xxxxz_xxx, t_xxxxz_xxy, t_xxxxz_xxz, t_xxxxz_xy, t_xxxxz_xyy, t_xxxxz_xyz, t_xxxxz_xz, t_xxxxz_xzz, t_xxxxz_yy, t_xxxxz_yyy, t_xxxxz_yyz, t_xxxxz_yz, t_xxxxz_yzz, t_xxxxz_zz, t_xxxxz_zzz, t_xxxxzz_xx, t_xxxxzz_xy, t_xxxxzz_xz, t_xxxxzz_yy, t_xxxxzz_yz, t_xxxxzz_zz, t_xxxyy_xx, t_xxxyy_xxx, t_xxxyy_xxy, t_xxxyy_xxz, t_xxxyy_xy, t_xxxyy_xyy, t_xxxyy_xyz, t_xxxyy_xz, t_xxxyy_xzz, t_xxxyy_yy, t_xxxyy_yyy, t_xxxyy_yyz, t_xxxyy_yz, t_xxxyy_yzz, t_xxxyy_zz, t_xxxyy_zzz, t_xxxyyy_xx, t_xxxyyy_xy, t_xxxyyy_xz, t_xxxyyy_yy, t_xxxyyy_yz, t_xxxyyy_zz, t_xxxyyz_xx, t_xxxyyz_xy, t_xxxyyz_xz, t_xxxyyz_yy, t_xxxyyz_yz, t_xxxyyz_zz, t_xxxyz_xx, t_xxxyz_xxx, t_xxxyz_xxy, t_xxxyz_xxz, t_xxxyz_xy, t_xxxyz_xyy, t_xxxyz_xyz, t_xxxyz_xz, t_xxxyz_xzz, t_xxxyz_yy, t_xxxyz_yyy, t_xxxyz_yyz, t_xxxyz_yz, t_xxxyz_yzz, t_xxxyz_zz, t_xxxyz_zzz, t_xxxyzz_xx, t_xxxyzz_xy, t_xxxyzz_xz, t_xxxyzz_yy, t_xxxyzz_yz, t_xxxyzz_zz, t_xxxzz_xx, t_xxxzz_xxx, t_xxxzz_xxy, t_xxxzz_xxz, t_xxxzz_xy, t_xxxzz_xyy, t_xxxzz_xyz, t_xxxzz_xz, t_xxxzz_xzz, t_xxxzz_yy, t_xxxzz_yyy, t_xxxzz_yyz, t_xxxzz_yz, t_xxxzz_yzz, t_xxxzz_zz, t_xxxzz_zzz, t_xxxzzz_xx, t_xxxzzz_xy, t_xxxzzz_xz, t_xxxzzz_yy, t_xxxzzz_yz, t_xxxzzz_zz, t_xxyyy_xx, t_xxyyy_xxx, t_xxyyy_xxy, t_xxyyy_xxz, t_xxyyy_xy, t_xxyyy_xyy, t_xxyyy_xyz, t_xxyyy_xz, t_xxyyy_xzz, t_xxyyy_yy, t_xxyyy_yyy, t_xxyyy_yyz, t_xxyyy_yz, t_xxyyy_yzz, t_xxyyy_zz, t_xxyyy_zzz, t_xxyyyy_xx, t_xxyyyy_xy, t_xxyyyy_xz, t_xxyyyy_yy, t_xxyyyy_yz, t_xxyyyy_zz, t_xxyyyz_xx, t_xxyyyz_xy, t_xxyyyz_xz, t_xxyyyz_yy, t_xxyyyz_yz, t_xxyyyz_zz, t_xxyyz_xx, t_xxyyz_xxx, t_xxyyz_xxy, t_xxyyz_xxz, t_xxyyz_xy, t_xxyyz_xyy, t_xxyyz_xyz, t_xxyyz_xz, t_xxyyz_xzz, t_xxyyz_yy, t_xxyyz_yyy, t_xxyyz_yyz, t_xxyyz_yz, t_xxyyz_yzz, t_xxyyz_zz, t_xxyyz_zzz, t_xxyyzz_xx, t_xxyyzz_xy, t_xxyyzz_xz, t_xxyyzz_yy, t_xxyyzz_yz, t_xxyyzz_zz, t_xxyzz_xx, t_xxyzz_xxx, t_xxyzz_xxy, t_xxyzz_xxz, t_xxyzz_xy, t_xxyzz_xyy, t_xxyzz_xyz, t_xxyzz_xz, t_xxyzz_xzz, t_xxyzz_yy, t_xxyzz_yyy, t_xxyzz_yyz, t_xxyzz_yz, t_xxyzz_yzz, t_xxyzz_zz, t_xxyzz_zzz, t_xxyzzz_xx, t_xxyzzz_xy, t_xxyzzz_xz, t_xxyzzz_yy, t_xxyzzz_yz, t_xxyzzz_zz, t_xxzzz_xx, t_xxzzz_xxx, t_xxzzz_xxy, t_xxzzz_xxz, t_xxzzz_xy, t_xxzzz_xyy, t_xxzzz_xyz, t_xxzzz_xz, t_xxzzz_xzz, t_xxzzz_yy, t_xxzzz_yyy, t_xxzzz_yyz, t_xxzzz_yz, t_xxzzz_yzz, t_xxzzz_zz, t_xxzzz_zzz, t_xxzzzz_xx, t_xxzzzz_xy, t_xxzzzz_xz, t_xxzzzz_yy, t_xxzzzz_yz, t_xxzzzz_zz, t_xyyyy_xx, t_xyyyy_xxx, t_xyyyy_xxy, t_xyyyy_xxz, t_xyyyy_xy, t_xyyyy_xyy, t_xyyyy_xyz, t_xyyyy_xz, t_xyyyy_xzz, t_xyyyy_yy, t_xyyyy_yyy, t_xyyyy_yyz, t_xyyyy_yz, t_xyyyy_yzz, t_xyyyy_zz, t_xyyyy_zzz, t_xyyyyy_xx, t_xyyyyy_xy, t_xyyyyy_xz, t_xyyyyy_yy, t_xyyyyy_yz, t_xyyyyy_zz, t_xyyyyz_xx, t_xyyyyz_xy, t_xyyyyz_xz, t_xyyyyz_yy, t_xyyyyz_yz, t_xyyyyz_zz, t_xyyyz_xx, t_xyyyz_xxx, t_xyyyz_xxy, t_xyyyz_xxz, t_xyyyz_xy, t_xyyyz_xyy, t_xyyyz_xyz, t_xyyyz_xz, t_xyyyz_xzz, t_xyyyz_yy, t_xyyyz_yyy, t_xyyyz_yyz, t_xyyyz_yz, t_xyyyz_yzz, t_xyyyz_zz, t_xyyyz_zzz, t_xyyyzz_xx, t_xyyyzz_xy, t_xyyyzz_xz, t_xyyyzz_yy, t_xyyyzz_yz, t_xyyyzz_zz, t_xyyzz_xx, t_xyyzz_xxx, t_xyyzz_xxy, t_xyyzz_xxz, t_xyyzz_xy, t_xyyzz_xyy, t_xyyzz_xyz, t_xyyzz_xz, t_xyyzz_xzz, t_xyyzz_yy, t_xyyzz_yyy, t_xyyzz_yyz, t_xyyzz_yz, t_xyyzz_yzz, t_xyyzz_zz, t_xyyzz_zzz, t_xyyzzz_xx, t_xyyzzz_xy, t_xyyzzz_xz, t_xyyzzz_yy, t_xyyzzz_yz, t_xyyzzz_zz, t_xyzzz_xx, t_xyzzz_xxx, t_xyzzz_xxy, t_xyzzz_xxz, t_xyzzz_xy, t_xyzzz_xyy, t_xyzzz_xyz, t_xyzzz_xz, t_xyzzz_xzz, t_xyzzz_yy, t_xyzzz_yyy, t_xyzzz_yyz, t_xyzzz_yz, t_xyzzz_yzz, t_xyzzz_zz, t_xyzzz_zzz, t_xyzzzz_xx, t_xyzzzz_xy, t_xyzzzz_xz, t_xyzzzz_yy, t_xyzzzz_yz, t_xyzzzz_zz, t_xzzzz_xx, t_xzzzz_xxx, t_xzzzz_xxy, t_xzzzz_xxz, t_xzzzz_xy, t_xzzzz_xyy, t_xzzzz_xyz, t_xzzzz_xz, t_xzzzz_xzz, t_xzzzz_yy, t_xzzzz_yyy, t_xzzzz_yyz, t_xzzzz_yz, t_xzzzz_yzz, t_xzzzz_zz, t_xzzzz_zzz, t_xzzzzz_xx, t_xzzzzz_xy, t_xzzzzz_xz, t_xzzzzz_yy, t_xzzzzz_yz, t_xzzzzz_zz, t_yyyyy_xx, t_yyyyy_xxx, t_yyyyy_xxy, t_yyyyy_xxz, t_yyyyy_xy, t_yyyyy_xyy, t_yyyyy_xyz, t_yyyyy_xz, t_yyyyy_xzz, t_yyyyy_yy, t_yyyyy_yyy, t_yyyyy_yyz, t_yyyyy_yz, t_yyyyy_yzz, t_yyyyy_zz, t_yyyyy_zzz, t_yyyyyy_yy, t_yyyyyy_yz, t_yyyyyy_zz, t_yyyyyz_yy, t_yyyyyz_yz, t_yyyyyz_zz, t_yyyyz_xx, t_yyyyz_xxx, t_yyyyz_xxy, t_yyyyz_xxz, t_yyyyz_xy, t_yyyyz_xyy, t_yyyyz_xyz, t_yyyyz_xz, t_yyyyz_xzz, t_yyyyz_yy, t_yyyyz_yyy, t_yyyyz_yyz, t_yyyyz_yz, t_yyyyz_yzz, t_yyyyz_zz, t_yyyyz_zzz, t_yyyyzz_yy, t_yyyyzz_yz, t_yyyyzz_zz, t_yyyzz_xx, t_yyyzz_xxx, t_yyyzz_xxy, t_yyyzz_xxz, t_yyyzz_xy, t_yyyzz_xyy, t_yyyzz_xyz, t_yyyzz_xz, t_yyyzz_xzz, t_yyyzz_yy, t_yyyzz_yyy, t_yyyzz_yyz, t_yyyzz_yz, t_yyyzz_yzz, t_yyyzz_zz, t_yyyzz_zzz, t_yyyzzz_yy, t_yyyzzz_yz, t_yyyzzz_zz, t_yyzzz_xx, t_yyzzz_xxx, t_yyzzz_xxy, t_yyzzz_xxz, t_yyzzz_xy, t_yyzzz_xyy, t_yyzzz_xyz, t_yyzzz_xz, t_yyzzz_xzz, t_yyzzz_yy, t_yyzzz_yyy, t_yyzzz_yyz, t_yyzzz_yz, t_yyzzz_yzz, t_yyzzz_zz, t_yyzzz_zzz, t_yyzzzz_yy, t_yyzzzz_yz, t_yyzzzz_zz, t_yzzzz_xx, t_yzzzz_xxx, t_yzzzz_xxy, t_yzzzz_xxz, t_yzzzz_xy, t_yzzzz_xyy, t_yzzzz_xyz, t_yzzzz_xz, t_yzzzz_xzz, t_yzzzz_yy, t_yzzzz_yyy, t_yzzzz_yyz, t_yzzzz_yz, t_yzzzz_yzz, t_yzzzz_zz, t_yzzzz_zzz, t_yzzzzz_yy, t_yzzzzz_yz, t_yzzzzz_zz, t_zzzzz_xx, t_zzzzz_xxx, t_zzzzz_xxy, t_zzzzz_xxz, t_zzzzz_xy, t_zzzzz_xyy, t_zzzzz_xyz, t_zzzzz_xz, t_zzzzz_xzz, t_zzzzz_yy, t_zzzzz_yyy, t_zzzzz_yyz, t_zzzzz_yz, t_zzzzz_yzz, t_zzzzz_zz, t_zzzzz_zzz, t_zzzzzz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxxx_xxx[i] = t_xxxxx_xx[i] * ab_x[i] + t_xxxxxx_xx[i];

        t_xxxxx_xxy[i] = t_xxxxx_xy[i] * ab_x[i] + t_xxxxxx_xy[i];

        t_xxxxx_xxz[i] = t_xxxxx_xz[i] * ab_x[i] + t_xxxxxx_xz[i];

        t_xxxxx_xyy[i] = t_xxxxx_yy[i] * ab_x[i] + t_xxxxxx_yy[i];

        t_xxxxx_xyz[i] = t_xxxxx_yz[i] * ab_x[i] + t_xxxxxx_yz[i];

        t_xxxxx_xzz[i] = t_xxxxx_zz[i] * ab_x[i] + t_xxxxxx_zz[i];

        t_xxxxx_yyy[i] = t_xxxxx_yy[i] * ab_y[i] + t_xxxxxy_yy[i];

        t_xxxxx_yyz[i] = t_xxxxx_yz[i] * ab_y[i] + t_xxxxxy_yz[i];

        t_xxxxx_yzz[i] = t_xxxxx_zz[i] * ab_y[i] + t_xxxxxy_zz[i];

        t_xxxxx_zzz[i] = t_xxxxx_zz[i] * ab_z[i] + t_xxxxxz_zz[i];

        t_xxxxy_xxx[i] = t_xxxxy_xx[i] * ab_x[i] + t_xxxxxy_xx[i];

        t_xxxxy_xxy[i] = t_xxxxy_xy[i] * ab_x[i] + t_xxxxxy_xy[i];

        t_xxxxy_xxz[i] = t_xxxxy_xz[i] * ab_x[i] + t_xxxxxy_xz[i];

        t_xxxxy_xyy[i] = t_xxxxy_yy[i] * ab_x[i] + t_xxxxxy_yy[i];

        t_xxxxy_xyz[i] = t_xxxxy_yz[i] * ab_x[i] + t_xxxxxy_yz[i];

        t_xxxxy_xzz[i] = t_xxxxy_zz[i] * ab_x[i] + t_xxxxxy_zz[i];

        t_xxxxy_yyy[i] = t_xxxxy_yy[i] * ab_y[i] + t_xxxxyy_yy[i];

        t_xxxxy_yyz[i] = t_xxxxy_yz[i] * ab_y[i] + t_xxxxyy_yz[i];

        t_xxxxy_yzz[i] = t_xxxxy_zz[i] * ab_y[i] + t_xxxxyy_zz[i];

        t_xxxxy_zzz[i] = t_xxxxy_zz[i] * ab_z[i] + t_xxxxyz_zz[i];

        t_xxxxz_xxx[i] = t_xxxxz_xx[i] * ab_x[i] + t_xxxxxz_xx[i];

        t_xxxxz_xxy[i] = t_xxxxz_xy[i] * ab_x[i] + t_xxxxxz_xy[i];

        t_xxxxz_xxz[i] = t_xxxxz_xz[i] * ab_x[i] + t_xxxxxz_xz[i];

        t_xxxxz_xyy[i] = t_xxxxz_yy[i] * ab_x[i] + t_xxxxxz_yy[i];

        t_xxxxz_xyz[i] = t_xxxxz_yz[i] * ab_x[i] + t_xxxxxz_yz[i];

        t_xxxxz_xzz[i] = t_xxxxz_zz[i] * ab_x[i] + t_xxxxxz_zz[i];

        t_xxxxz_yyy[i] = t_xxxxz_yy[i] * ab_y[i] + t_xxxxyz_yy[i];

        t_xxxxz_yyz[i] = t_xxxxz_yz[i] * ab_y[i] + t_xxxxyz_yz[i];

        t_xxxxz_yzz[i] = t_xxxxz_zz[i] * ab_y[i] + t_xxxxyz_zz[i];

        t_xxxxz_zzz[i] = t_xxxxz_zz[i] * ab_z[i] + t_xxxxzz_zz[i];

        t_xxxyy_xxx[i] = t_xxxyy_xx[i] * ab_x[i] + t_xxxxyy_xx[i];

        t_xxxyy_xxy[i] = t_xxxyy_xy[i] * ab_x[i] + t_xxxxyy_xy[i];

        t_xxxyy_xxz[i] = t_xxxyy_xz[i] * ab_x[i] + t_xxxxyy_xz[i];

        t_xxxyy_xyy[i] = t_xxxyy_yy[i] * ab_x[i] + t_xxxxyy_yy[i];

        t_xxxyy_xyz[i] = t_xxxyy_yz[i] * ab_x[i] + t_xxxxyy_yz[i];

        t_xxxyy_xzz[i] = t_xxxyy_zz[i] * ab_x[i] + t_xxxxyy_zz[i];

        t_xxxyy_yyy[i] = t_xxxyy_yy[i] * ab_y[i] + t_xxxyyy_yy[i];

        t_xxxyy_yyz[i] = t_xxxyy_yz[i] * ab_y[i] + t_xxxyyy_yz[i];

        t_xxxyy_yzz[i] = t_xxxyy_zz[i] * ab_y[i] + t_xxxyyy_zz[i];

        t_xxxyy_zzz[i] = t_xxxyy_zz[i] * ab_z[i] + t_xxxyyz_zz[i];

        t_xxxyz_xxx[i] = t_xxxyz_xx[i] * ab_x[i] + t_xxxxyz_xx[i];

        t_xxxyz_xxy[i] = t_xxxyz_xy[i] * ab_x[i] + t_xxxxyz_xy[i];

        t_xxxyz_xxz[i] = t_xxxyz_xz[i] * ab_x[i] + t_xxxxyz_xz[i];

        t_xxxyz_xyy[i] = t_xxxyz_yy[i] * ab_x[i] + t_xxxxyz_yy[i];

        t_xxxyz_xyz[i] = t_xxxyz_yz[i] * ab_x[i] + t_xxxxyz_yz[i];

        t_xxxyz_xzz[i] = t_xxxyz_zz[i] * ab_x[i] + t_xxxxyz_zz[i];

        t_xxxyz_yyy[i] = t_xxxyz_yy[i] * ab_y[i] + t_xxxyyz_yy[i];

        t_xxxyz_yyz[i] = t_xxxyz_yz[i] * ab_y[i] + t_xxxyyz_yz[i];

        t_xxxyz_yzz[i] = t_xxxyz_zz[i] * ab_y[i] + t_xxxyyz_zz[i];

        t_xxxyz_zzz[i] = t_xxxyz_zz[i] * ab_z[i] + t_xxxyzz_zz[i];

        t_xxxzz_xxx[i] = t_xxxzz_xx[i] * ab_x[i] + t_xxxxzz_xx[i];

        t_xxxzz_xxy[i] = t_xxxzz_xy[i] * ab_x[i] + t_xxxxzz_xy[i];

        t_xxxzz_xxz[i] = t_xxxzz_xz[i] * ab_x[i] + t_xxxxzz_xz[i];

        t_xxxzz_xyy[i] = t_xxxzz_yy[i] * ab_x[i] + t_xxxxzz_yy[i];

        t_xxxzz_xyz[i] = t_xxxzz_yz[i] * ab_x[i] + t_xxxxzz_yz[i];

        t_xxxzz_xzz[i] = t_xxxzz_zz[i] * ab_x[i] + t_xxxxzz_zz[i];

        t_xxxzz_yyy[i] = t_xxxzz_yy[i] * ab_y[i] + t_xxxyzz_yy[i];

        t_xxxzz_yyz[i] = t_xxxzz_yz[i] * ab_y[i] + t_xxxyzz_yz[i];

        t_xxxzz_yzz[i] = t_xxxzz_zz[i] * ab_y[i] + t_xxxyzz_zz[i];

        t_xxxzz_zzz[i] = t_xxxzz_zz[i] * ab_z[i] + t_xxxzzz_zz[i];

        t_xxyyy_xxx[i] = t_xxyyy_xx[i] * ab_x[i] + t_xxxyyy_xx[i];

        t_xxyyy_xxy[i] = t_xxyyy_xy[i] * ab_x[i] + t_xxxyyy_xy[i];

        t_xxyyy_xxz[i] = t_xxyyy_xz[i] * ab_x[i] + t_xxxyyy_xz[i];

        t_xxyyy_xyy[i] = t_xxyyy_yy[i] * ab_x[i] + t_xxxyyy_yy[i];

        t_xxyyy_xyz[i] = t_xxyyy_yz[i] * ab_x[i] + t_xxxyyy_yz[i];

        t_xxyyy_xzz[i] = t_xxyyy_zz[i] * ab_x[i] + t_xxxyyy_zz[i];

        t_xxyyy_yyy[i] = t_xxyyy_yy[i] * ab_y[i] + t_xxyyyy_yy[i];

        t_xxyyy_yyz[i] = t_xxyyy_yz[i] * ab_y[i] + t_xxyyyy_yz[i];

        t_xxyyy_yzz[i] = t_xxyyy_zz[i] * ab_y[i] + t_xxyyyy_zz[i];

        t_xxyyy_zzz[i] = t_xxyyy_zz[i] * ab_z[i] + t_xxyyyz_zz[i];

        t_xxyyz_xxx[i] = t_xxyyz_xx[i] * ab_x[i] + t_xxxyyz_xx[i];

        t_xxyyz_xxy[i] = t_xxyyz_xy[i] * ab_x[i] + t_xxxyyz_xy[i];

        t_xxyyz_xxz[i] = t_xxyyz_xz[i] * ab_x[i] + t_xxxyyz_xz[i];

        t_xxyyz_xyy[i] = t_xxyyz_yy[i] * ab_x[i] + t_xxxyyz_yy[i];

        t_xxyyz_xyz[i] = t_xxyyz_yz[i] * ab_x[i] + t_xxxyyz_yz[i];

        t_xxyyz_xzz[i] = t_xxyyz_zz[i] * ab_x[i] + t_xxxyyz_zz[i];

        t_xxyyz_yyy[i] = t_xxyyz_yy[i] * ab_y[i] + t_xxyyyz_yy[i];

        t_xxyyz_yyz[i] = t_xxyyz_yz[i] * ab_y[i] + t_xxyyyz_yz[i];

        t_xxyyz_yzz[i] = t_xxyyz_zz[i] * ab_y[i] + t_xxyyyz_zz[i];

        t_xxyyz_zzz[i] = t_xxyyz_zz[i] * ab_z[i] + t_xxyyzz_zz[i];

        t_xxyzz_xxx[i] = t_xxyzz_xx[i] * ab_x[i] + t_xxxyzz_xx[i];

        t_xxyzz_xxy[i] = t_xxyzz_xy[i] * ab_x[i] + t_xxxyzz_xy[i];

        t_xxyzz_xxz[i] = t_xxyzz_xz[i] * ab_x[i] + t_xxxyzz_xz[i];

        t_xxyzz_xyy[i] = t_xxyzz_yy[i] * ab_x[i] + t_xxxyzz_yy[i];

        t_xxyzz_xyz[i] = t_xxyzz_yz[i] * ab_x[i] + t_xxxyzz_yz[i];

        t_xxyzz_xzz[i] = t_xxyzz_zz[i] * ab_x[i] + t_xxxyzz_zz[i];

        t_xxyzz_yyy[i] = t_xxyzz_yy[i] * ab_y[i] + t_xxyyzz_yy[i];

        t_xxyzz_yyz[i] = t_xxyzz_yz[i] * ab_y[i] + t_xxyyzz_yz[i];

        t_xxyzz_yzz[i] = t_xxyzz_zz[i] * ab_y[i] + t_xxyyzz_zz[i];

        t_xxyzz_zzz[i] = t_xxyzz_zz[i] * ab_z[i] + t_xxyzzz_zz[i];

        t_xxzzz_xxx[i] = t_xxzzz_xx[i] * ab_x[i] + t_xxxzzz_xx[i];

        t_xxzzz_xxy[i] = t_xxzzz_xy[i] * ab_x[i] + t_xxxzzz_xy[i];

        t_xxzzz_xxz[i] = t_xxzzz_xz[i] * ab_x[i] + t_xxxzzz_xz[i];

        t_xxzzz_xyy[i] = t_xxzzz_yy[i] * ab_x[i] + t_xxxzzz_yy[i];

        t_xxzzz_xyz[i] = t_xxzzz_yz[i] * ab_x[i] + t_xxxzzz_yz[i];

        t_xxzzz_xzz[i] = t_xxzzz_zz[i] * ab_x[i] + t_xxxzzz_zz[i];

        t_xxzzz_yyy[i] = t_xxzzz_yy[i] * ab_y[i] + t_xxyzzz_yy[i];

        t_xxzzz_yyz[i] = t_xxzzz_yz[i] * ab_y[i] + t_xxyzzz_yz[i];

        t_xxzzz_yzz[i] = t_xxzzz_zz[i] * ab_y[i] + t_xxyzzz_zz[i];

        t_xxzzz_zzz[i] = t_xxzzz_zz[i] * ab_z[i] + t_xxzzzz_zz[i];

        t_xyyyy_xxx[i] = t_xyyyy_xx[i] * ab_x[i] + t_xxyyyy_xx[i];

        t_xyyyy_xxy[i] = t_xyyyy_xy[i] * ab_x[i] + t_xxyyyy_xy[i];

        t_xyyyy_xxz[i] = t_xyyyy_xz[i] * ab_x[i] + t_xxyyyy_xz[i];

        t_xyyyy_xyy[i] = t_xyyyy_yy[i] * ab_x[i] + t_xxyyyy_yy[i];

        t_xyyyy_xyz[i] = t_xyyyy_yz[i] * ab_x[i] + t_xxyyyy_yz[i];

        t_xyyyy_xzz[i] = t_xyyyy_zz[i] * ab_x[i] + t_xxyyyy_zz[i];

        t_xyyyy_yyy[i] = t_xyyyy_yy[i] * ab_y[i] + t_xyyyyy_yy[i];

        t_xyyyy_yyz[i] = t_xyyyy_yz[i] * ab_y[i] + t_xyyyyy_yz[i];

        t_xyyyy_yzz[i] = t_xyyyy_zz[i] * ab_y[i] + t_xyyyyy_zz[i];

        t_xyyyy_zzz[i] = t_xyyyy_zz[i] * ab_z[i] + t_xyyyyz_zz[i];

        t_xyyyz_xxx[i] = t_xyyyz_xx[i] * ab_x[i] + t_xxyyyz_xx[i];

        t_xyyyz_xxy[i] = t_xyyyz_xy[i] * ab_x[i] + t_xxyyyz_xy[i];

        t_xyyyz_xxz[i] = t_xyyyz_xz[i] * ab_x[i] + t_xxyyyz_xz[i];

        t_xyyyz_xyy[i] = t_xyyyz_yy[i] * ab_x[i] + t_xxyyyz_yy[i];

        t_xyyyz_xyz[i] = t_xyyyz_yz[i] * ab_x[i] + t_xxyyyz_yz[i];

        t_xyyyz_xzz[i] = t_xyyyz_zz[i] * ab_x[i] + t_xxyyyz_zz[i];

        t_xyyyz_yyy[i] = t_xyyyz_yy[i] * ab_y[i] + t_xyyyyz_yy[i];

        t_xyyyz_yyz[i] = t_xyyyz_yz[i] * ab_y[i] + t_xyyyyz_yz[i];

        t_xyyyz_yzz[i] = t_xyyyz_zz[i] * ab_y[i] + t_xyyyyz_zz[i];

        t_xyyyz_zzz[i] = t_xyyyz_zz[i] * ab_z[i] + t_xyyyzz_zz[i];

        t_xyyzz_xxx[i] = t_xyyzz_xx[i] * ab_x[i] + t_xxyyzz_xx[i];

        t_xyyzz_xxy[i] = t_xyyzz_xy[i] * ab_x[i] + t_xxyyzz_xy[i];

        t_xyyzz_xxz[i] = t_xyyzz_xz[i] * ab_x[i] + t_xxyyzz_xz[i];

        t_xyyzz_xyy[i] = t_xyyzz_yy[i] * ab_x[i] + t_xxyyzz_yy[i];

        t_xyyzz_xyz[i] = t_xyyzz_yz[i] * ab_x[i] + t_xxyyzz_yz[i];

        t_xyyzz_xzz[i] = t_xyyzz_zz[i] * ab_x[i] + t_xxyyzz_zz[i];

        t_xyyzz_yyy[i] = t_xyyzz_yy[i] * ab_y[i] + t_xyyyzz_yy[i];

        t_xyyzz_yyz[i] = t_xyyzz_yz[i] * ab_y[i] + t_xyyyzz_yz[i];

        t_xyyzz_yzz[i] = t_xyyzz_zz[i] * ab_y[i] + t_xyyyzz_zz[i];

        t_xyyzz_zzz[i] = t_xyyzz_zz[i] * ab_z[i] + t_xyyzzz_zz[i];

        t_xyzzz_xxx[i] = t_xyzzz_xx[i] * ab_x[i] + t_xxyzzz_xx[i];

        t_xyzzz_xxy[i] = t_xyzzz_xy[i] * ab_x[i] + t_xxyzzz_xy[i];

        t_xyzzz_xxz[i] = t_xyzzz_xz[i] * ab_x[i] + t_xxyzzz_xz[i];

        t_xyzzz_xyy[i] = t_xyzzz_yy[i] * ab_x[i] + t_xxyzzz_yy[i];

        t_xyzzz_xyz[i] = t_xyzzz_yz[i] * ab_x[i] + t_xxyzzz_yz[i];

        t_xyzzz_xzz[i] = t_xyzzz_zz[i] * ab_x[i] + t_xxyzzz_zz[i];

        t_xyzzz_yyy[i] = t_xyzzz_yy[i] * ab_y[i] + t_xyyzzz_yy[i];

        t_xyzzz_yyz[i] = t_xyzzz_yz[i] * ab_y[i] + t_xyyzzz_yz[i];

        t_xyzzz_yzz[i] = t_xyzzz_zz[i] * ab_y[i] + t_xyyzzz_zz[i];

        t_xyzzz_zzz[i] = t_xyzzz_zz[i] * ab_z[i] + t_xyzzzz_zz[i];

        t_xzzzz_xxx[i] = t_xzzzz_xx[i] * ab_x[i] + t_xxzzzz_xx[i];

        t_xzzzz_xxy[i] = t_xzzzz_xy[i] * ab_x[i] + t_xxzzzz_xy[i];

        t_xzzzz_xxz[i] = t_xzzzz_xz[i] * ab_x[i] + t_xxzzzz_xz[i];

        t_xzzzz_xyy[i] = t_xzzzz_yy[i] * ab_x[i] + t_xxzzzz_yy[i];

        t_xzzzz_xyz[i] = t_xzzzz_yz[i] * ab_x[i] + t_xxzzzz_yz[i];

        t_xzzzz_xzz[i] = t_xzzzz_zz[i] * ab_x[i] + t_xxzzzz_zz[i];

        t_xzzzz_yyy[i] = t_xzzzz_yy[i] * ab_y[i] + t_xyzzzz_yy[i];

        t_xzzzz_yyz[i] = t_xzzzz_yz[i] * ab_y[i] + t_xyzzzz_yz[i];

        t_xzzzz_yzz[i] = t_xzzzz_zz[i] * ab_y[i] + t_xyzzzz_zz[i];

        t_xzzzz_zzz[i] = t_xzzzz_zz[i] * ab_z[i] + t_xzzzzz_zz[i];

        t_yyyyy_xxx[i] = t_yyyyy_xx[i] * ab_x[i] + t_xyyyyy_xx[i];

        t_yyyyy_xxy[i] = t_yyyyy_xy[i] * ab_x[i] + t_xyyyyy_xy[i];

        t_yyyyy_xxz[i] = t_yyyyy_xz[i] * ab_x[i] + t_xyyyyy_xz[i];

        t_yyyyy_xyy[i] = t_yyyyy_yy[i] * ab_x[i] + t_xyyyyy_yy[i];

        t_yyyyy_xyz[i] = t_yyyyy_yz[i] * ab_x[i] + t_xyyyyy_yz[i];

        t_yyyyy_xzz[i] = t_yyyyy_zz[i] * ab_x[i] + t_xyyyyy_zz[i];

        t_yyyyy_yyy[i] = t_yyyyy_yy[i] * ab_y[i] + t_yyyyyy_yy[i];

        t_yyyyy_yyz[i] = t_yyyyy_yz[i] * ab_y[i] + t_yyyyyy_yz[i];

        t_yyyyy_yzz[i] = t_yyyyy_zz[i] * ab_y[i] + t_yyyyyy_zz[i];

        t_yyyyy_zzz[i] = t_yyyyy_zz[i] * ab_z[i] + t_yyyyyz_zz[i];

        t_yyyyz_xxx[i] = t_yyyyz_xx[i] * ab_x[i] + t_xyyyyz_xx[i];

        t_yyyyz_xxy[i] = t_yyyyz_xy[i] * ab_x[i] + t_xyyyyz_xy[i];

        t_yyyyz_xxz[i] = t_yyyyz_xz[i] * ab_x[i] + t_xyyyyz_xz[i];

        t_yyyyz_xyy[i] = t_yyyyz_yy[i] * ab_x[i] + t_xyyyyz_yy[i];

        t_yyyyz_xyz[i] = t_yyyyz_yz[i] * ab_x[i] + t_xyyyyz_yz[i];

        t_yyyyz_xzz[i] = t_yyyyz_zz[i] * ab_x[i] + t_xyyyyz_zz[i];

        t_yyyyz_yyy[i] = t_yyyyz_yy[i] * ab_y[i] + t_yyyyyz_yy[i];

        t_yyyyz_yyz[i] = t_yyyyz_yz[i] * ab_y[i] + t_yyyyyz_yz[i];

        t_yyyyz_yzz[i] = t_yyyyz_zz[i] * ab_y[i] + t_yyyyyz_zz[i];

        t_yyyyz_zzz[i] = t_yyyyz_zz[i] * ab_z[i] + t_yyyyzz_zz[i];

        t_yyyzz_xxx[i] = t_yyyzz_xx[i] * ab_x[i] + t_xyyyzz_xx[i];

        t_yyyzz_xxy[i] = t_yyyzz_xy[i] * ab_x[i] + t_xyyyzz_xy[i];

        t_yyyzz_xxz[i] = t_yyyzz_xz[i] * ab_x[i] + t_xyyyzz_xz[i];

        t_yyyzz_xyy[i] = t_yyyzz_yy[i] * ab_x[i] + t_xyyyzz_yy[i];

        t_yyyzz_xyz[i] = t_yyyzz_yz[i] * ab_x[i] + t_xyyyzz_yz[i];

        t_yyyzz_xzz[i] = t_yyyzz_zz[i] * ab_x[i] + t_xyyyzz_zz[i];

        t_yyyzz_yyy[i] = t_yyyzz_yy[i] * ab_y[i] + t_yyyyzz_yy[i];

        t_yyyzz_yyz[i] = t_yyyzz_yz[i] * ab_y[i] + t_yyyyzz_yz[i];

        t_yyyzz_yzz[i] = t_yyyzz_zz[i] * ab_y[i] + t_yyyyzz_zz[i];

        t_yyyzz_zzz[i] = t_yyyzz_zz[i] * ab_z[i] + t_yyyzzz_zz[i];

        t_yyzzz_xxx[i] = t_yyzzz_xx[i] * ab_x[i] + t_xyyzzz_xx[i];

        t_yyzzz_xxy[i] = t_yyzzz_xy[i] * ab_x[i] + t_xyyzzz_xy[i];

        t_yyzzz_xxz[i] = t_yyzzz_xz[i] * ab_x[i] + t_xyyzzz_xz[i];

        t_yyzzz_xyy[i] = t_yyzzz_yy[i] * ab_x[i] + t_xyyzzz_yy[i];

        t_yyzzz_xyz[i] = t_yyzzz_yz[i] * ab_x[i] + t_xyyzzz_yz[i];

        t_yyzzz_xzz[i] = t_yyzzz_zz[i] * ab_x[i] + t_xyyzzz_zz[i];

        t_yyzzz_yyy[i] = t_yyzzz_yy[i] * ab_y[i] + t_yyyzzz_yy[i];

        t_yyzzz_yyz[i] = t_yyzzz_yz[i] * ab_y[i] + t_yyyzzz_yz[i];

        t_yyzzz_yzz[i] = t_yyzzz_zz[i] * ab_y[i] + t_yyyzzz_zz[i];

        t_yyzzz_zzz[i] = t_yyzzz_zz[i] * ab_z[i] + t_yyzzzz_zz[i];

        t_yzzzz_xxx[i] = t_yzzzz_xx[i] * ab_x[i] + t_xyzzzz_xx[i];

        t_yzzzz_xxy[i] = t_yzzzz_xy[i] * ab_x[i] + t_xyzzzz_xy[i];

        t_yzzzz_xxz[i] = t_yzzzz_xz[i] * ab_x[i] + t_xyzzzz_xz[i];

        t_yzzzz_xyy[i] = t_yzzzz_yy[i] * ab_x[i] + t_xyzzzz_yy[i];

        t_yzzzz_xyz[i] = t_yzzzz_yz[i] * ab_x[i] + t_xyzzzz_yz[i];

        t_yzzzz_xzz[i] = t_yzzzz_zz[i] * ab_x[i] + t_xyzzzz_zz[i];

        t_yzzzz_yyy[i] = t_yzzzz_yy[i] * ab_y[i] + t_yyzzzz_yy[i];

        t_yzzzz_yyz[i] = t_yzzzz_yz[i] * ab_y[i] + t_yyzzzz_yz[i];

        t_yzzzz_yzz[i] = t_yzzzz_zz[i] * ab_y[i] + t_yyzzzz_zz[i];

        t_yzzzz_zzz[i] = t_yzzzz_zz[i] * ab_z[i] + t_yzzzzz_zz[i];

        t_zzzzz_xxx[i] = t_zzzzz_xx[i] * ab_x[i] + t_xzzzzz_xx[i];

        t_zzzzz_xxy[i] = t_zzzzz_xy[i] * ab_x[i] + t_xzzzzz_xy[i];

        t_zzzzz_xxz[i] = t_zzzzz_xz[i] * ab_x[i] + t_xzzzzz_xz[i];

        t_zzzzz_xyy[i] = t_zzzzz_yy[i] * ab_x[i] + t_xzzzzz_yy[i];

        t_zzzzz_xyz[i] = t_zzzzz_yz[i] * ab_x[i] + t_xzzzzz_yz[i];

        t_zzzzz_xzz[i] = t_zzzzz_zz[i] * ab_x[i] + t_xzzzzz_zz[i];

        t_zzzzz_yyy[i] = t_zzzzz_yy[i] * ab_y[i] + t_yzzzzz_yy[i];

        t_zzzzz_yyz[i] = t_zzzzz_yz[i] * ab_y[i] + t_yzzzzz_yz[i];

        t_zzzzz_yzz[i] = t_zzzzz_zz[i] * ab_y[i] + t_yzzzzz_zz[i];

        t_zzzzz_zzz[i] = t_zzzzz_zz[i] * ab_z[i] + t_zzzzzz_zz[i];
    }
}

} // t2chrr namespace

