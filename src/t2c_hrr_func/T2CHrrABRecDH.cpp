#include "T2CHrrABRecDH.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_dh(CSimdArray<double>& cbuffer, 
            const size_t idx_dh,
            const size_t idx_ph,
            const size_t idx_pi,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : PH

    auto t_x_xxxxx = cbuffer.data(idx_ph);

    auto t_x_xxxxy = cbuffer.data(idx_ph + 1);

    auto t_x_xxxxz = cbuffer.data(idx_ph + 2);

    auto t_x_xxxyy = cbuffer.data(idx_ph + 3);

    auto t_x_xxxyz = cbuffer.data(idx_ph + 4);

    auto t_x_xxxzz = cbuffer.data(idx_ph + 5);

    auto t_x_xxyyy = cbuffer.data(idx_ph + 6);

    auto t_x_xxyyz = cbuffer.data(idx_ph + 7);

    auto t_x_xxyzz = cbuffer.data(idx_ph + 8);

    auto t_x_xxzzz = cbuffer.data(idx_ph + 9);

    auto t_x_xyyyy = cbuffer.data(idx_ph + 10);

    auto t_x_xyyyz = cbuffer.data(idx_ph + 11);

    auto t_x_xyyzz = cbuffer.data(idx_ph + 12);

    auto t_x_xyzzz = cbuffer.data(idx_ph + 13);

    auto t_x_xzzzz = cbuffer.data(idx_ph + 14);

    auto t_x_yyyyy = cbuffer.data(idx_ph + 15);

    auto t_x_yyyyz = cbuffer.data(idx_ph + 16);

    auto t_x_yyyzz = cbuffer.data(idx_ph + 17);

    auto t_x_yyzzz = cbuffer.data(idx_ph + 18);

    auto t_x_yzzzz = cbuffer.data(idx_ph + 19);

    auto t_x_zzzzz = cbuffer.data(idx_ph + 20);

    auto t_y_xxxxx = cbuffer.data(idx_ph + 21);

    auto t_y_xxxxy = cbuffer.data(idx_ph + 22);

    auto t_y_xxxxz = cbuffer.data(idx_ph + 23);

    auto t_y_xxxyy = cbuffer.data(idx_ph + 24);

    auto t_y_xxxyz = cbuffer.data(idx_ph + 25);

    auto t_y_xxxzz = cbuffer.data(idx_ph + 26);

    auto t_y_xxyyy = cbuffer.data(idx_ph + 27);

    auto t_y_xxyyz = cbuffer.data(idx_ph + 28);

    auto t_y_xxyzz = cbuffer.data(idx_ph + 29);

    auto t_y_xxzzz = cbuffer.data(idx_ph + 30);

    auto t_y_xyyyy = cbuffer.data(idx_ph + 31);

    auto t_y_xyyyz = cbuffer.data(idx_ph + 32);

    auto t_y_xyyzz = cbuffer.data(idx_ph + 33);

    auto t_y_xyzzz = cbuffer.data(idx_ph + 34);

    auto t_y_xzzzz = cbuffer.data(idx_ph + 35);

    auto t_y_yyyyy = cbuffer.data(idx_ph + 36);

    auto t_y_yyyyz = cbuffer.data(idx_ph + 37);

    auto t_y_yyyzz = cbuffer.data(idx_ph + 38);

    auto t_y_yyzzz = cbuffer.data(idx_ph + 39);

    auto t_y_yzzzz = cbuffer.data(idx_ph + 40);

    auto t_y_zzzzz = cbuffer.data(idx_ph + 41);

    auto t_z_xxxxx = cbuffer.data(idx_ph + 42);

    auto t_z_xxxxy = cbuffer.data(idx_ph + 43);

    auto t_z_xxxxz = cbuffer.data(idx_ph + 44);

    auto t_z_xxxyy = cbuffer.data(idx_ph + 45);

    auto t_z_xxxyz = cbuffer.data(idx_ph + 46);

    auto t_z_xxxzz = cbuffer.data(idx_ph + 47);

    auto t_z_xxyyy = cbuffer.data(idx_ph + 48);

    auto t_z_xxyyz = cbuffer.data(idx_ph + 49);

    auto t_z_xxyzz = cbuffer.data(idx_ph + 50);

    auto t_z_xxzzz = cbuffer.data(idx_ph + 51);

    auto t_z_xyyyy = cbuffer.data(idx_ph + 52);

    auto t_z_xyyyz = cbuffer.data(idx_ph + 53);

    auto t_z_xyyzz = cbuffer.data(idx_ph + 54);

    auto t_z_xyzzz = cbuffer.data(idx_ph + 55);

    auto t_z_xzzzz = cbuffer.data(idx_ph + 56);

    auto t_z_yyyyy = cbuffer.data(idx_ph + 57);

    auto t_z_yyyyz = cbuffer.data(idx_ph + 58);

    auto t_z_yyyzz = cbuffer.data(idx_ph + 59);

    auto t_z_yyzzz = cbuffer.data(idx_ph + 60);

    auto t_z_yzzzz = cbuffer.data(idx_ph + 61);

    auto t_z_zzzzz = cbuffer.data(idx_ph + 62);

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

    // Set up components of targeted buffer : DH

    auto t_xx_xxxxx = cbuffer.data(idx_dh);

    auto t_xx_xxxxy = cbuffer.data(idx_dh + 1);

    auto t_xx_xxxxz = cbuffer.data(idx_dh + 2);

    auto t_xx_xxxyy = cbuffer.data(idx_dh + 3);

    auto t_xx_xxxyz = cbuffer.data(idx_dh + 4);

    auto t_xx_xxxzz = cbuffer.data(idx_dh + 5);

    auto t_xx_xxyyy = cbuffer.data(idx_dh + 6);

    auto t_xx_xxyyz = cbuffer.data(idx_dh + 7);

    auto t_xx_xxyzz = cbuffer.data(idx_dh + 8);

    auto t_xx_xxzzz = cbuffer.data(idx_dh + 9);

    auto t_xx_xyyyy = cbuffer.data(idx_dh + 10);

    auto t_xx_xyyyz = cbuffer.data(idx_dh + 11);

    auto t_xx_xyyzz = cbuffer.data(idx_dh + 12);

    auto t_xx_xyzzz = cbuffer.data(idx_dh + 13);

    auto t_xx_xzzzz = cbuffer.data(idx_dh + 14);

    auto t_xx_yyyyy = cbuffer.data(idx_dh + 15);

    auto t_xx_yyyyz = cbuffer.data(idx_dh + 16);

    auto t_xx_yyyzz = cbuffer.data(idx_dh + 17);

    auto t_xx_yyzzz = cbuffer.data(idx_dh + 18);

    auto t_xx_yzzzz = cbuffer.data(idx_dh + 19);

    auto t_xx_zzzzz = cbuffer.data(idx_dh + 20);

    auto t_xy_xxxxx = cbuffer.data(idx_dh + 21);

    auto t_xy_xxxxy = cbuffer.data(idx_dh + 22);

    auto t_xy_xxxxz = cbuffer.data(idx_dh + 23);

    auto t_xy_xxxyy = cbuffer.data(idx_dh + 24);

    auto t_xy_xxxyz = cbuffer.data(idx_dh + 25);

    auto t_xy_xxxzz = cbuffer.data(idx_dh + 26);

    auto t_xy_xxyyy = cbuffer.data(idx_dh + 27);

    auto t_xy_xxyyz = cbuffer.data(idx_dh + 28);

    auto t_xy_xxyzz = cbuffer.data(idx_dh + 29);

    auto t_xy_xxzzz = cbuffer.data(idx_dh + 30);

    auto t_xy_xyyyy = cbuffer.data(idx_dh + 31);

    auto t_xy_xyyyz = cbuffer.data(idx_dh + 32);

    auto t_xy_xyyzz = cbuffer.data(idx_dh + 33);

    auto t_xy_xyzzz = cbuffer.data(idx_dh + 34);

    auto t_xy_xzzzz = cbuffer.data(idx_dh + 35);

    auto t_xy_yyyyy = cbuffer.data(idx_dh + 36);

    auto t_xy_yyyyz = cbuffer.data(idx_dh + 37);

    auto t_xy_yyyzz = cbuffer.data(idx_dh + 38);

    auto t_xy_yyzzz = cbuffer.data(idx_dh + 39);

    auto t_xy_yzzzz = cbuffer.data(idx_dh + 40);

    auto t_xy_zzzzz = cbuffer.data(idx_dh + 41);

    auto t_xz_xxxxx = cbuffer.data(idx_dh + 42);

    auto t_xz_xxxxy = cbuffer.data(idx_dh + 43);

    auto t_xz_xxxxz = cbuffer.data(idx_dh + 44);

    auto t_xz_xxxyy = cbuffer.data(idx_dh + 45);

    auto t_xz_xxxyz = cbuffer.data(idx_dh + 46);

    auto t_xz_xxxzz = cbuffer.data(idx_dh + 47);

    auto t_xz_xxyyy = cbuffer.data(idx_dh + 48);

    auto t_xz_xxyyz = cbuffer.data(idx_dh + 49);

    auto t_xz_xxyzz = cbuffer.data(idx_dh + 50);

    auto t_xz_xxzzz = cbuffer.data(idx_dh + 51);

    auto t_xz_xyyyy = cbuffer.data(idx_dh + 52);

    auto t_xz_xyyyz = cbuffer.data(idx_dh + 53);

    auto t_xz_xyyzz = cbuffer.data(idx_dh + 54);

    auto t_xz_xyzzz = cbuffer.data(idx_dh + 55);

    auto t_xz_xzzzz = cbuffer.data(idx_dh + 56);

    auto t_xz_yyyyy = cbuffer.data(idx_dh + 57);

    auto t_xz_yyyyz = cbuffer.data(idx_dh + 58);

    auto t_xz_yyyzz = cbuffer.data(idx_dh + 59);

    auto t_xz_yyzzz = cbuffer.data(idx_dh + 60);

    auto t_xz_yzzzz = cbuffer.data(idx_dh + 61);

    auto t_xz_zzzzz = cbuffer.data(idx_dh + 62);

    auto t_yy_xxxxx = cbuffer.data(idx_dh + 63);

    auto t_yy_xxxxy = cbuffer.data(idx_dh + 64);

    auto t_yy_xxxxz = cbuffer.data(idx_dh + 65);

    auto t_yy_xxxyy = cbuffer.data(idx_dh + 66);

    auto t_yy_xxxyz = cbuffer.data(idx_dh + 67);

    auto t_yy_xxxzz = cbuffer.data(idx_dh + 68);

    auto t_yy_xxyyy = cbuffer.data(idx_dh + 69);

    auto t_yy_xxyyz = cbuffer.data(idx_dh + 70);

    auto t_yy_xxyzz = cbuffer.data(idx_dh + 71);

    auto t_yy_xxzzz = cbuffer.data(idx_dh + 72);

    auto t_yy_xyyyy = cbuffer.data(idx_dh + 73);

    auto t_yy_xyyyz = cbuffer.data(idx_dh + 74);

    auto t_yy_xyyzz = cbuffer.data(idx_dh + 75);

    auto t_yy_xyzzz = cbuffer.data(idx_dh + 76);

    auto t_yy_xzzzz = cbuffer.data(idx_dh + 77);

    auto t_yy_yyyyy = cbuffer.data(idx_dh + 78);

    auto t_yy_yyyyz = cbuffer.data(idx_dh + 79);

    auto t_yy_yyyzz = cbuffer.data(idx_dh + 80);

    auto t_yy_yyzzz = cbuffer.data(idx_dh + 81);

    auto t_yy_yzzzz = cbuffer.data(idx_dh + 82);

    auto t_yy_zzzzz = cbuffer.data(idx_dh + 83);

    auto t_yz_xxxxx = cbuffer.data(idx_dh + 84);

    auto t_yz_xxxxy = cbuffer.data(idx_dh + 85);

    auto t_yz_xxxxz = cbuffer.data(idx_dh + 86);

    auto t_yz_xxxyy = cbuffer.data(idx_dh + 87);

    auto t_yz_xxxyz = cbuffer.data(idx_dh + 88);

    auto t_yz_xxxzz = cbuffer.data(idx_dh + 89);

    auto t_yz_xxyyy = cbuffer.data(idx_dh + 90);

    auto t_yz_xxyyz = cbuffer.data(idx_dh + 91);

    auto t_yz_xxyzz = cbuffer.data(idx_dh + 92);

    auto t_yz_xxzzz = cbuffer.data(idx_dh + 93);

    auto t_yz_xyyyy = cbuffer.data(idx_dh + 94);

    auto t_yz_xyyyz = cbuffer.data(idx_dh + 95);

    auto t_yz_xyyzz = cbuffer.data(idx_dh + 96);

    auto t_yz_xyzzz = cbuffer.data(idx_dh + 97);

    auto t_yz_xzzzz = cbuffer.data(idx_dh + 98);

    auto t_yz_yyyyy = cbuffer.data(idx_dh + 99);

    auto t_yz_yyyyz = cbuffer.data(idx_dh + 100);

    auto t_yz_yyyzz = cbuffer.data(idx_dh + 101);

    auto t_yz_yyzzz = cbuffer.data(idx_dh + 102);

    auto t_yz_yzzzz = cbuffer.data(idx_dh + 103);

    auto t_yz_zzzzz = cbuffer.data(idx_dh + 104);

    auto t_zz_xxxxx = cbuffer.data(idx_dh + 105);

    auto t_zz_xxxxy = cbuffer.data(idx_dh + 106);

    auto t_zz_xxxxz = cbuffer.data(idx_dh + 107);

    auto t_zz_xxxyy = cbuffer.data(idx_dh + 108);

    auto t_zz_xxxyz = cbuffer.data(idx_dh + 109);

    auto t_zz_xxxzz = cbuffer.data(idx_dh + 110);

    auto t_zz_xxyyy = cbuffer.data(idx_dh + 111);

    auto t_zz_xxyyz = cbuffer.data(idx_dh + 112);

    auto t_zz_xxyzz = cbuffer.data(idx_dh + 113);

    auto t_zz_xxzzz = cbuffer.data(idx_dh + 114);

    auto t_zz_xyyyy = cbuffer.data(idx_dh + 115);

    auto t_zz_xyyyz = cbuffer.data(idx_dh + 116);

    auto t_zz_xyyzz = cbuffer.data(idx_dh + 117);

    auto t_zz_xyzzz = cbuffer.data(idx_dh + 118);

    auto t_zz_xzzzz = cbuffer.data(idx_dh + 119);

    auto t_zz_yyyyy = cbuffer.data(idx_dh + 120);

    auto t_zz_yyyyz = cbuffer.data(idx_dh + 121);

    auto t_zz_yyyzz = cbuffer.data(idx_dh + 122);

    auto t_zz_yyzzz = cbuffer.data(idx_dh + 123);

    auto t_zz_yzzzz = cbuffer.data(idx_dh + 124);

    auto t_zz_zzzzz = cbuffer.data(idx_dh + 125);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_x_xxxxx, t_x_xxxxxx, t_x_xxxxxy, t_x_xxxxxz, t_x_xxxxy, t_x_xxxxyy, t_x_xxxxyz, t_x_xxxxz, t_x_xxxxzz, t_x_xxxyy, t_x_xxxyyy, t_x_xxxyyz, t_x_xxxyz, t_x_xxxyzz, t_x_xxxzz, t_x_xxxzzz, t_x_xxyyy, t_x_xxyyyy, t_x_xxyyyz, t_x_xxyyz, t_x_xxyyzz, t_x_xxyzz, t_x_xxyzzz, t_x_xxzzz, t_x_xxzzzz, t_x_xyyyy, t_x_xyyyyy, t_x_xyyyyz, t_x_xyyyz, t_x_xyyyzz, t_x_xyyzz, t_x_xyyzzz, t_x_xyzzz, t_x_xyzzzz, t_x_xzzzz, t_x_xzzzzz, t_x_yyyyy, t_x_yyyyz, t_x_yyyzz, t_x_yyzzz, t_x_yzzzz, t_x_zzzzz, t_xx_xxxxx, t_xx_xxxxy, t_xx_xxxxz, t_xx_xxxyy, t_xx_xxxyz, t_xx_xxxzz, t_xx_xxyyy, t_xx_xxyyz, t_xx_xxyzz, t_xx_xxzzz, t_xx_xyyyy, t_xx_xyyyz, t_xx_xyyzz, t_xx_xyzzz, t_xx_xzzzz, t_xx_yyyyy, t_xx_yyyyz, t_xx_yyyzz, t_xx_yyzzz, t_xx_yzzzz, t_xx_zzzzz, t_xy_xxxxx, t_xy_xxxxy, t_xy_xxxxz, t_xy_xxxyy, t_xy_xxxyz, t_xy_xxxzz, t_xy_xxyyy, t_xy_xxyyz, t_xy_xxyzz, t_xy_xxzzz, t_xy_xyyyy, t_xy_xyyyz, t_xy_xyyzz, t_xy_xyzzz, t_xy_xzzzz, t_xy_yyyyy, t_xy_yyyyz, t_xy_yyyzz, t_xy_yyzzz, t_xy_yzzzz, t_xy_zzzzz, t_xz_xxxxx, t_xz_xxxxy, t_xz_xxxxz, t_xz_xxxyy, t_xz_xxxyz, t_xz_xxxzz, t_xz_xxyyy, t_xz_xxyyz, t_xz_xxyzz, t_xz_xxzzz, t_xz_xyyyy, t_xz_xyyyz, t_xz_xyyzz, t_xz_xyzzz, t_xz_xzzzz, t_xz_yyyyy, t_xz_yyyyz, t_xz_yyyzz, t_xz_yyzzz, t_xz_yzzzz, t_xz_zzzzz, t_y_xxxxx, t_y_xxxxxx, t_y_xxxxxy, t_y_xxxxxz, t_y_xxxxy, t_y_xxxxyy, t_y_xxxxyz, t_y_xxxxz, t_y_xxxxzz, t_y_xxxyy, t_y_xxxyyy, t_y_xxxyyz, t_y_xxxyz, t_y_xxxyzz, t_y_xxxzz, t_y_xxxzzz, t_y_xxyyy, t_y_xxyyyy, t_y_xxyyyz, t_y_xxyyz, t_y_xxyyzz, t_y_xxyzz, t_y_xxyzzz, t_y_xxzzz, t_y_xxzzzz, t_y_xyyyy, t_y_xyyyyy, t_y_xyyyyz, t_y_xyyyz, t_y_xyyyzz, t_y_xyyzz, t_y_xyyzzz, t_y_xyzzz, t_y_xyzzzz, t_y_xzzzz, t_y_xzzzzz, t_y_yyyyy, t_y_yyyyyy, t_y_yyyyyz, t_y_yyyyz, t_y_yyyyzz, t_y_yyyzz, t_y_yyyzzz, t_y_yyzzz, t_y_yyzzzz, t_y_yzzzz, t_y_yzzzzz, t_y_zzzzz, t_yy_xxxxx, t_yy_xxxxy, t_yy_xxxxz, t_yy_xxxyy, t_yy_xxxyz, t_yy_xxxzz, t_yy_xxyyy, t_yy_xxyyz, t_yy_xxyzz, t_yy_xxzzz, t_yy_xyyyy, t_yy_xyyyz, t_yy_xyyzz, t_yy_xyzzz, t_yy_xzzzz, t_yy_yyyyy, t_yy_yyyyz, t_yy_yyyzz, t_yy_yyzzz, t_yy_yzzzz, t_yy_zzzzz, t_yz_xxxxx, t_yz_xxxxy, t_yz_xxxxz, t_yz_xxxyy, t_yz_xxxyz, t_yz_xxxzz, t_yz_xxyyy, t_yz_xxyyz, t_yz_xxyzz, t_yz_xxzzz, t_yz_xyyyy, t_yz_xyyyz, t_yz_xyyzz, t_yz_xyzzz, t_yz_xzzzz, t_yz_yyyyy, t_yz_yyyyz, t_yz_yyyzz, t_yz_yyzzz, t_yz_yzzzz, t_yz_zzzzz, t_z_xxxxx, t_z_xxxxxx, t_z_xxxxxy, t_z_xxxxxz, t_z_xxxxy, t_z_xxxxyy, t_z_xxxxyz, t_z_xxxxz, t_z_xxxxzz, t_z_xxxyy, t_z_xxxyyy, t_z_xxxyyz, t_z_xxxyz, t_z_xxxyzz, t_z_xxxzz, t_z_xxxzzz, t_z_xxyyy, t_z_xxyyyy, t_z_xxyyyz, t_z_xxyyz, t_z_xxyyzz, t_z_xxyzz, t_z_xxyzzz, t_z_xxzzz, t_z_xxzzzz, t_z_xyyyy, t_z_xyyyyy, t_z_xyyyyz, t_z_xyyyz, t_z_xyyyzz, t_z_xyyzz, t_z_xyyzzz, t_z_xyzzz, t_z_xyzzzz, t_z_xzzzz, t_z_xzzzzz, t_z_yyyyy, t_z_yyyyyy, t_z_yyyyyz, t_z_yyyyz, t_z_yyyyzz, t_z_yyyzz, t_z_yyyzzz, t_z_yyzzz, t_z_yyzzzz, t_z_yzzzz, t_z_yzzzzz, t_z_zzzzz, t_z_zzzzzz, t_zz_xxxxx, t_zz_xxxxy, t_zz_xxxxz, t_zz_xxxyy, t_zz_xxxyz, t_zz_xxxzz, t_zz_xxyyy, t_zz_xxyyz, t_zz_xxyzz, t_zz_xxzzz, t_zz_xyyyy, t_zz_xyyyz, t_zz_xyyzz, t_zz_xyzzz, t_zz_xzzzz, t_zz_yyyyy, t_zz_yyyyz, t_zz_yyyzz, t_zz_yyzzz, t_zz_yzzzz, t_zz_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xx_xxxxx[i] = -t_x_xxxxx[i] * ab_x[i] + t_x_xxxxxx[i];

        t_xx_xxxxy[i] = -t_x_xxxxy[i] * ab_x[i] + t_x_xxxxxy[i];

        t_xx_xxxxz[i] = -t_x_xxxxz[i] * ab_x[i] + t_x_xxxxxz[i];

        t_xx_xxxyy[i] = -t_x_xxxyy[i] * ab_x[i] + t_x_xxxxyy[i];

        t_xx_xxxyz[i] = -t_x_xxxyz[i] * ab_x[i] + t_x_xxxxyz[i];

        t_xx_xxxzz[i] = -t_x_xxxzz[i] * ab_x[i] + t_x_xxxxzz[i];

        t_xx_xxyyy[i] = -t_x_xxyyy[i] * ab_x[i] + t_x_xxxyyy[i];

        t_xx_xxyyz[i] = -t_x_xxyyz[i] * ab_x[i] + t_x_xxxyyz[i];

        t_xx_xxyzz[i] = -t_x_xxyzz[i] * ab_x[i] + t_x_xxxyzz[i];

        t_xx_xxzzz[i] = -t_x_xxzzz[i] * ab_x[i] + t_x_xxxzzz[i];

        t_xx_xyyyy[i] = -t_x_xyyyy[i] * ab_x[i] + t_x_xxyyyy[i];

        t_xx_xyyyz[i] = -t_x_xyyyz[i] * ab_x[i] + t_x_xxyyyz[i];

        t_xx_xyyzz[i] = -t_x_xyyzz[i] * ab_x[i] + t_x_xxyyzz[i];

        t_xx_xyzzz[i] = -t_x_xyzzz[i] * ab_x[i] + t_x_xxyzzz[i];

        t_xx_xzzzz[i] = -t_x_xzzzz[i] * ab_x[i] + t_x_xxzzzz[i];

        t_xx_yyyyy[i] = -t_x_yyyyy[i] * ab_x[i] + t_x_xyyyyy[i];

        t_xx_yyyyz[i] = -t_x_yyyyz[i] * ab_x[i] + t_x_xyyyyz[i];

        t_xx_yyyzz[i] = -t_x_yyyzz[i] * ab_x[i] + t_x_xyyyzz[i];

        t_xx_yyzzz[i] = -t_x_yyzzz[i] * ab_x[i] + t_x_xyyzzz[i];

        t_xx_yzzzz[i] = -t_x_yzzzz[i] * ab_x[i] + t_x_xyzzzz[i];

        t_xx_zzzzz[i] = -t_x_zzzzz[i] * ab_x[i] + t_x_xzzzzz[i];

        t_xy_xxxxx[i] = -t_y_xxxxx[i] * ab_x[i] + t_y_xxxxxx[i];

        t_xy_xxxxy[i] = -t_y_xxxxy[i] * ab_x[i] + t_y_xxxxxy[i];

        t_xy_xxxxz[i] = -t_y_xxxxz[i] * ab_x[i] + t_y_xxxxxz[i];

        t_xy_xxxyy[i] = -t_y_xxxyy[i] * ab_x[i] + t_y_xxxxyy[i];

        t_xy_xxxyz[i] = -t_y_xxxyz[i] * ab_x[i] + t_y_xxxxyz[i];

        t_xy_xxxzz[i] = -t_y_xxxzz[i] * ab_x[i] + t_y_xxxxzz[i];

        t_xy_xxyyy[i] = -t_y_xxyyy[i] * ab_x[i] + t_y_xxxyyy[i];

        t_xy_xxyyz[i] = -t_y_xxyyz[i] * ab_x[i] + t_y_xxxyyz[i];

        t_xy_xxyzz[i] = -t_y_xxyzz[i] * ab_x[i] + t_y_xxxyzz[i];

        t_xy_xxzzz[i] = -t_y_xxzzz[i] * ab_x[i] + t_y_xxxzzz[i];

        t_xy_xyyyy[i] = -t_y_xyyyy[i] * ab_x[i] + t_y_xxyyyy[i];

        t_xy_xyyyz[i] = -t_y_xyyyz[i] * ab_x[i] + t_y_xxyyyz[i];

        t_xy_xyyzz[i] = -t_y_xyyzz[i] * ab_x[i] + t_y_xxyyzz[i];

        t_xy_xyzzz[i] = -t_y_xyzzz[i] * ab_x[i] + t_y_xxyzzz[i];

        t_xy_xzzzz[i] = -t_y_xzzzz[i] * ab_x[i] + t_y_xxzzzz[i];

        t_xy_yyyyy[i] = -t_y_yyyyy[i] * ab_x[i] + t_y_xyyyyy[i];

        t_xy_yyyyz[i] = -t_y_yyyyz[i] * ab_x[i] + t_y_xyyyyz[i];

        t_xy_yyyzz[i] = -t_y_yyyzz[i] * ab_x[i] + t_y_xyyyzz[i];

        t_xy_yyzzz[i] = -t_y_yyzzz[i] * ab_x[i] + t_y_xyyzzz[i];

        t_xy_yzzzz[i] = -t_y_yzzzz[i] * ab_x[i] + t_y_xyzzzz[i];

        t_xy_zzzzz[i] = -t_y_zzzzz[i] * ab_x[i] + t_y_xzzzzz[i];

        t_xz_xxxxx[i] = -t_z_xxxxx[i] * ab_x[i] + t_z_xxxxxx[i];

        t_xz_xxxxy[i] = -t_z_xxxxy[i] * ab_x[i] + t_z_xxxxxy[i];

        t_xz_xxxxz[i] = -t_z_xxxxz[i] * ab_x[i] + t_z_xxxxxz[i];

        t_xz_xxxyy[i] = -t_z_xxxyy[i] * ab_x[i] + t_z_xxxxyy[i];

        t_xz_xxxyz[i] = -t_z_xxxyz[i] * ab_x[i] + t_z_xxxxyz[i];

        t_xz_xxxzz[i] = -t_z_xxxzz[i] * ab_x[i] + t_z_xxxxzz[i];

        t_xz_xxyyy[i] = -t_z_xxyyy[i] * ab_x[i] + t_z_xxxyyy[i];

        t_xz_xxyyz[i] = -t_z_xxyyz[i] * ab_x[i] + t_z_xxxyyz[i];

        t_xz_xxyzz[i] = -t_z_xxyzz[i] * ab_x[i] + t_z_xxxyzz[i];

        t_xz_xxzzz[i] = -t_z_xxzzz[i] * ab_x[i] + t_z_xxxzzz[i];

        t_xz_xyyyy[i] = -t_z_xyyyy[i] * ab_x[i] + t_z_xxyyyy[i];

        t_xz_xyyyz[i] = -t_z_xyyyz[i] * ab_x[i] + t_z_xxyyyz[i];

        t_xz_xyyzz[i] = -t_z_xyyzz[i] * ab_x[i] + t_z_xxyyzz[i];

        t_xz_xyzzz[i] = -t_z_xyzzz[i] * ab_x[i] + t_z_xxyzzz[i];

        t_xz_xzzzz[i] = -t_z_xzzzz[i] * ab_x[i] + t_z_xxzzzz[i];

        t_xz_yyyyy[i] = -t_z_yyyyy[i] * ab_x[i] + t_z_xyyyyy[i];

        t_xz_yyyyz[i] = -t_z_yyyyz[i] * ab_x[i] + t_z_xyyyyz[i];

        t_xz_yyyzz[i] = -t_z_yyyzz[i] * ab_x[i] + t_z_xyyyzz[i];

        t_xz_yyzzz[i] = -t_z_yyzzz[i] * ab_x[i] + t_z_xyyzzz[i];

        t_xz_yzzzz[i] = -t_z_yzzzz[i] * ab_x[i] + t_z_xyzzzz[i];

        t_xz_zzzzz[i] = -t_z_zzzzz[i] * ab_x[i] + t_z_xzzzzz[i];

        t_yy_xxxxx[i] = -t_y_xxxxx[i] * ab_y[i] + t_y_xxxxxy[i];

        t_yy_xxxxy[i] = -t_y_xxxxy[i] * ab_y[i] + t_y_xxxxyy[i];

        t_yy_xxxxz[i] = -t_y_xxxxz[i] * ab_y[i] + t_y_xxxxyz[i];

        t_yy_xxxyy[i] = -t_y_xxxyy[i] * ab_y[i] + t_y_xxxyyy[i];

        t_yy_xxxyz[i] = -t_y_xxxyz[i] * ab_y[i] + t_y_xxxyyz[i];

        t_yy_xxxzz[i] = -t_y_xxxzz[i] * ab_y[i] + t_y_xxxyzz[i];

        t_yy_xxyyy[i] = -t_y_xxyyy[i] * ab_y[i] + t_y_xxyyyy[i];

        t_yy_xxyyz[i] = -t_y_xxyyz[i] * ab_y[i] + t_y_xxyyyz[i];

        t_yy_xxyzz[i] = -t_y_xxyzz[i] * ab_y[i] + t_y_xxyyzz[i];

        t_yy_xxzzz[i] = -t_y_xxzzz[i] * ab_y[i] + t_y_xxyzzz[i];

        t_yy_xyyyy[i] = -t_y_xyyyy[i] * ab_y[i] + t_y_xyyyyy[i];

        t_yy_xyyyz[i] = -t_y_xyyyz[i] * ab_y[i] + t_y_xyyyyz[i];

        t_yy_xyyzz[i] = -t_y_xyyzz[i] * ab_y[i] + t_y_xyyyzz[i];

        t_yy_xyzzz[i] = -t_y_xyzzz[i] * ab_y[i] + t_y_xyyzzz[i];

        t_yy_xzzzz[i] = -t_y_xzzzz[i] * ab_y[i] + t_y_xyzzzz[i];

        t_yy_yyyyy[i] = -t_y_yyyyy[i] * ab_y[i] + t_y_yyyyyy[i];

        t_yy_yyyyz[i] = -t_y_yyyyz[i] * ab_y[i] + t_y_yyyyyz[i];

        t_yy_yyyzz[i] = -t_y_yyyzz[i] * ab_y[i] + t_y_yyyyzz[i];

        t_yy_yyzzz[i] = -t_y_yyzzz[i] * ab_y[i] + t_y_yyyzzz[i];

        t_yy_yzzzz[i] = -t_y_yzzzz[i] * ab_y[i] + t_y_yyzzzz[i];

        t_yy_zzzzz[i] = -t_y_zzzzz[i] * ab_y[i] + t_y_yzzzzz[i];

        t_yz_xxxxx[i] = -t_z_xxxxx[i] * ab_y[i] + t_z_xxxxxy[i];

        t_yz_xxxxy[i] = -t_z_xxxxy[i] * ab_y[i] + t_z_xxxxyy[i];

        t_yz_xxxxz[i] = -t_z_xxxxz[i] * ab_y[i] + t_z_xxxxyz[i];

        t_yz_xxxyy[i] = -t_z_xxxyy[i] * ab_y[i] + t_z_xxxyyy[i];

        t_yz_xxxyz[i] = -t_z_xxxyz[i] * ab_y[i] + t_z_xxxyyz[i];

        t_yz_xxxzz[i] = -t_z_xxxzz[i] * ab_y[i] + t_z_xxxyzz[i];

        t_yz_xxyyy[i] = -t_z_xxyyy[i] * ab_y[i] + t_z_xxyyyy[i];

        t_yz_xxyyz[i] = -t_z_xxyyz[i] * ab_y[i] + t_z_xxyyyz[i];

        t_yz_xxyzz[i] = -t_z_xxyzz[i] * ab_y[i] + t_z_xxyyzz[i];

        t_yz_xxzzz[i] = -t_z_xxzzz[i] * ab_y[i] + t_z_xxyzzz[i];

        t_yz_xyyyy[i] = -t_z_xyyyy[i] * ab_y[i] + t_z_xyyyyy[i];

        t_yz_xyyyz[i] = -t_z_xyyyz[i] * ab_y[i] + t_z_xyyyyz[i];

        t_yz_xyyzz[i] = -t_z_xyyzz[i] * ab_y[i] + t_z_xyyyzz[i];

        t_yz_xyzzz[i] = -t_z_xyzzz[i] * ab_y[i] + t_z_xyyzzz[i];

        t_yz_xzzzz[i] = -t_z_xzzzz[i] * ab_y[i] + t_z_xyzzzz[i];

        t_yz_yyyyy[i] = -t_z_yyyyy[i] * ab_y[i] + t_z_yyyyyy[i];

        t_yz_yyyyz[i] = -t_z_yyyyz[i] * ab_y[i] + t_z_yyyyyz[i];

        t_yz_yyyzz[i] = -t_z_yyyzz[i] * ab_y[i] + t_z_yyyyzz[i];

        t_yz_yyzzz[i] = -t_z_yyzzz[i] * ab_y[i] + t_z_yyyzzz[i];

        t_yz_yzzzz[i] = -t_z_yzzzz[i] * ab_y[i] + t_z_yyzzzz[i];

        t_yz_zzzzz[i] = -t_z_zzzzz[i] * ab_y[i] + t_z_yzzzzz[i];

        t_zz_xxxxx[i] = -t_z_xxxxx[i] * ab_z[i] + t_z_xxxxxz[i];

        t_zz_xxxxy[i] = -t_z_xxxxy[i] * ab_z[i] + t_z_xxxxyz[i];

        t_zz_xxxxz[i] = -t_z_xxxxz[i] * ab_z[i] + t_z_xxxxzz[i];

        t_zz_xxxyy[i] = -t_z_xxxyy[i] * ab_z[i] + t_z_xxxyyz[i];

        t_zz_xxxyz[i] = -t_z_xxxyz[i] * ab_z[i] + t_z_xxxyzz[i];

        t_zz_xxxzz[i] = -t_z_xxxzz[i] * ab_z[i] + t_z_xxxzzz[i];

        t_zz_xxyyy[i] = -t_z_xxyyy[i] * ab_z[i] + t_z_xxyyyz[i];

        t_zz_xxyyz[i] = -t_z_xxyyz[i] * ab_z[i] + t_z_xxyyzz[i];

        t_zz_xxyzz[i] = -t_z_xxyzz[i] * ab_z[i] + t_z_xxyzzz[i];

        t_zz_xxzzz[i] = -t_z_xxzzz[i] * ab_z[i] + t_z_xxzzzz[i];

        t_zz_xyyyy[i] = -t_z_xyyyy[i] * ab_z[i] + t_z_xyyyyz[i];

        t_zz_xyyyz[i] = -t_z_xyyyz[i] * ab_z[i] + t_z_xyyyzz[i];

        t_zz_xyyzz[i] = -t_z_xyyzz[i] * ab_z[i] + t_z_xyyzzz[i];

        t_zz_xyzzz[i] = -t_z_xyzzz[i] * ab_z[i] + t_z_xyzzzz[i];

        t_zz_xzzzz[i] = -t_z_xzzzz[i] * ab_z[i] + t_z_xzzzzz[i];

        t_zz_yyyyy[i] = -t_z_yyyyy[i] * ab_z[i] + t_z_yyyyyz[i];

        t_zz_yyyyz[i] = -t_z_yyyyz[i] * ab_z[i] + t_z_yyyyzz[i];

        t_zz_yyyzz[i] = -t_z_yyyzz[i] * ab_z[i] + t_z_yyyzzz[i];

        t_zz_yyzzz[i] = -t_z_yyzzz[i] * ab_z[i] + t_z_yyzzzz[i];

        t_zz_yzzzz[i] = -t_z_yzzzz[i] * ab_z[i] + t_z_yzzzzz[i];

        t_zz_zzzzz[i] = -t_z_zzzzz[i] * ab_z[i] + t_z_zzzzzz[i];
    }
}

} // t2chrr namespace

