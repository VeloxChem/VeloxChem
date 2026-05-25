#include "LocalCorePotentialPrimRecID.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_id(CSimdArray<double>& pbuffer, 
                                  const size_t idx_id,
                                  const size_t idx_gd,
                                  const size_t idx_hp,
                                  const size_t idx_hd,
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

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx = pbuffer.data(idx_gd);

    auto tg_xxxx_xy = pbuffer.data(idx_gd + 1);

    auto tg_xxxx_xz = pbuffer.data(idx_gd + 2);

    auto tg_xxxx_yy = pbuffer.data(idx_gd + 3);

    auto tg_xxxx_yz = pbuffer.data(idx_gd + 4);

    auto tg_xxxx_zz = pbuffer.data(idx_gd + 5);

    auto tg_xxxy_xx = pbuffer.data(idx_gd + 6);

    auto tg_xxxy_xz = pbuffer.data(idx_gd + 8);

    auto tg_xxxy_yy = pbuffer.data(idx_gd + 9);

    auto tg_xxxy_yz = pbuffer.data(idx_gd + 10);

    auto tg_xxxz_xx = pbuffer.data(idx_gd + 12);

    auto tg_xxxz_xy = pbuffer.data(idx_gd + 13);

    auto tg_xxxz_xz = pbuffer.data(idx_gd + 14);

    auto tg_xxxz_yz = pbuffer.data(idx_gd + 16);

    auto tg_xxxz_zz = pbuffer.data(idx_gd + 17);

    auto tg_xxyy_xx = pbuffer.data(idx_gd + 18);

    auto tg_xxyy_xy = pbuffer.data(idx_gd + 19);

    auto tg_xxyy_xz = pbuffer.data(idx_gd + 20);

    auto tg_xxyy_yy = pbuffer.data(idx_gd + 21);

    auto tg_xxyy_yz = pbuffer.data(idx_gd + 22);

    auto tg_xxyy_zz = pbuffer.data(idx_gd + 23);

    auto tg_xxyz_xz = pbuffer.data(idx_gd + 26);

    auto tg_xxyz_yz = pbuffer.data(idx_gd + 28);

    auto tg_xxzz_xx = pbuffer.data(idx_gd + 30);

    auto tg_xxzz_xy = pbuffer.data(idx_gd + 31);

    auto tg_xxzz_xz = pbuffer.data(idx_gd + 32);

    auto tg_xxzz_yy = pbuffer.data(idx_gd + 33);

    auto tg_xxzz_yz = pbuffer.data(idx_gd + 34);

    auto tg_xxzz_zz = pbuffer.data(idx_gd + 35);

    auto tg_xyyy_xy = pbuffer.data(idx_gd + 37);

    auto tg_xyyy_yy = pbuffer.data(idx_gd + 39);

    auto tg_xyyy_yz = pbuffer.data(idx_gd + 40);

    auto tg_xyyy_zz = pbuffer.data(idx_gd + 41);

    auto tg_xyyz_yz = pbuffer.data(idx_gd + 46);

    auto tg_xyyz_zz = pbuffer.data(idx_gd + 47);

    auto tg_xyzz_yy = pbuffer.data(idx_gd + 51);

    auto tg_xyzz_yz = pbuffer.data(idx_gd + 52);

    auto tg_xzzz_xz = pbuffer.data(idx_gd + 56);

    auto tg_xzzz_yy = pbuffer.data(idx_gd + 57);

    auto tg_xzzz_yz = pbuffer.data(idx_gd + 58);

    auto tg_xzzz_zz = pbuffer.data(idx_gd + 59);

    auto tg_yyyy_xx = pbuffer.data(idx_gd + 60);

    auto tg_yyyy_xy = pbuffer.data(idx_gd + 61);

    auto tg_yyyy_xz = pbuffer.data(idx_gd + 62);

    auto tg_yyyy_yy = pbuffer.data(idx_gd + 63);

    auto tg_yyyy_yz = pbuffer.data(idx_gd + 64);

    auto tg_yyyy_zz = pbuffer.data(idx_gd + 65);

    auto tg_yyyz_xy = pbuffer.data(idx_gd + 67);

    auto tg_yyyz_xz = pbuffer.data(idx_gd + 68);

    auto tg_yyyz_yy = pbuffer.data(idx_gd + 69);

    auto tg_yyyz_yz = pbuffer.data(idx_gd + 70);

    auto tg_yyyz_zz = pbuffer.data(idx_gd + 71);

    auto tg_yyzz_xx = pbuffer.data(idx_gd + 72);

    auto tg_yyzz_xy = pbuffer.data(idx_gd + 73);

    auto tg_yyzz_xz = pbuffer.data(idx_gd + 74);

    auto tg_yyzz_yy = pbuffer.data(idx_gd + 75);

    auto tg_yyzz_yz = pbuffer.data(idx_gd + 76);

    auto tg_yyzz_zz = pbuffer.data(idx_gd + 77);

    auto tg_yzzz_xx = pbuffer.data(idx_gd + 78);

    auto tg_yzzz_xz = pbuffer.data(idx_gd + 80);

    auto tg_yzzz_yy = pbuffer.data(idx_gd + 81);

    auto tg_yzzz_yz = pbuffer.data(idx_gd + 82);

    auto tg_yzzz_zz = pbuffer.data(idx_gd + 83);

    auto tg_zzzz_xx = pbuffer.data(idx_gd + 84);

    auto tg_zzzz_xy = pbuffer.data(idx_gd + 85);

    auto tg_zzzz_xz = pbuffer.data(idx_gd + 86);

    auto tg_zzzz_yy = pbuffer.data(idx_gd + 87);

    auto tg_zzzz_yz = pbuffer.data(idx_gd + 88);

    auto tg_zzzz_zz = pbuffer.data(idx_gd + 89);

    // Set up components of auxiliary buffer : HP

    auto tg_xxxxx_x = pbuffer.data(idx_hp);

    auto tg_xxxxx_y = pbuffer.data(idx_hp + 1);

    auto tg_xxxxx_z = pbuffer.data(idx_hp + 2);

    auto tg_xxxyy_y = pbuffer.data(idx_hp + 10);

    auto tg_xxxzz_x = pbuffer.data(idx_hp + 15);

    auto tg_xxxzz_z = pbuffer.data(idx_hp + 17);

    auto tg_xxyyy_y = pbuffer.data(idx_hp + 19);

    auto tg_xxzzz_x = pbuffer.data(idx_hp + 27);

    auto tg_xxzzz_z = pbuffer.data(idx_hp + 29);

    auto tg_xyyyy_y = pbuffer.data(idx_hp + 31);

    auto tg_xzzzz_z = pbuffer.data(idx_hp + 44);

    auto tg_yyyyy_x = pbuffer.data(idx_hp + 45);

    auto tg_yyyyy_y = pbuffer.data(idx_hp + 46);

    auto tg_yyyyy_z = pbuffer.data(idx_hp + 47);

    auto tg_yyyyz_z = pbuffer.data(idx_hp + 50);

    auto tg_yyyzz_x = pbuffer.data(idx_hp + 51);

    auto tg_yyyzz_y = pbuffer.data(idx_hp + 52);

    auto tg_yyyzz_z = pbuffer.data(idx_hp + 53);

    auto tg_yyzzz_x = pbuffer.data(idx_hp + 54);

    auto tg_yyzzz_y = pbuffer.data(idx_hp + 55);

    auto tg_yyzzz_z = pbuffer.data(idx_hp + 56);

    auto tg_yzzzz_y = pbuffer.data(idx_hp + 58);

    auto tg_yzzzz_z = pbuffer.data(idx_hp + 59);

    auto tg_zzzzz_x = pbuffer.data(idx_hp + 60);

    auto tg_zzzzz_y = pbuffer.data(idx_hp + 61);

    auto tg_zzzzz_z = pbuffer.data(idx_hp + 62);

    // Set up components of auxiliary buffer : HD

    auto tg_xxxxx_xx = pbuffer.data(idx_hd);

    auto tg_xxxxx_xy = pbuffer.data(idx_hd + 1);

    auto tg_xxxxx_xz = pbuffer.data(idx_hd + 2);

    auto tg_xxxxx_yy = pbuffer.data(idx_hd + 3);

    auto tg_xxxxx_yz = pbuffer.data(idx_hd + 4);

    auto tg_xxxxx_zz = pbuffer.data(idx_hd + 5);

    auto tg_xxxxy_xx = pbuffer.data(idx_hd + 6);

    auto tg_xxxxy_xy = pbuffer.data(idx_hd + 7);

    auto tg_xxxxy_xz = pbuffer.data(idx_hd + 8);

    auto tg_xxxxy_yy = pbuffer.data(idx_hd + 9);

    auto tg_xxxxy_yz = pbuffer.data(idx_hd + 10);

    auto tg_xxxxz_xx = pbuffer.data(idx_hd + 12);

    auto tg_xxxxz_xy = pbuffer.data(idx_hd + 13);

    auto tg_xxxxz_xz = pbuffer.data(idx_hd + 14);

    auto tg_xxxxz_yz = pbuffer.data(idx_hd + 16);

    auto tg_xxxxz_zz = pbuffer.data(idx_hd + 17);

    auto tg_xxxyy_xx = pbuffer.data(idx_hd + 18);

    auto tg_xxxyy_xy = pbuffer.data(idx_hd + 19);

    auto tg_xxxyy_xz = pbuffer.data(idx_hd + 20);

    auto tg_xxxyy_yy = pbuffer.data(idx_hd + 21);

    auto tg_xxxyy_yz = pbuffer.data(idx_hd + 22);

    auto tg_xxxyy_zz = pbuffer.data(idx_hd + 23);

    auto tg_xxxyz_xz = pbuffer.data(idx_hd + 26);

    auto tg_xxxyz_yz = pbuffer.data(idx_hd + 28);

    auto tg_xxxzz_xx = pbuffer.data(idx_hd + 30);

    auto tg_xxxzz_xy = pbuffer.data(idx_hd + 31);

    auto tg_xxxzz_xz = pbuffer.data(idx_hd + 32);

    auto tg_xxxzz_yy = pbuffer.data(idx_hd + 33);

    auto tg_xxxzz_yz = pbuffer.data(idx_hd + 34);

    auto tg_xxxzz_zz = pbuffer.data(idx_hd + 35);

    auto tg_xxyyy_xx = pbuffer.data(idx_hd + 36);

    auto tg_xxyyy_xy = pbuffer.data(idx_hd + 37);

    auto tg_xxyyy_xz = pbuffer.data(idx_hd + 38);

    auto tg_xxyyy_yy = pbuffer.data(idx_hd + 39);

    auto tg_xxyyy_yz = pbuffer.data(idx_hd + 40);

    auto tg_xxyyy_zz = pbuffer.data(idx_hd + 41);

    auto tg_xxyyz_xy = pbuffer.data(idx_hd + 43);

    auto tg_xxyyz_xz = pbuffer.data(idx_hd + 44);

    auto tg_xxyyz_yz = pbuffer.data(idx_hd + 46);

    auto tg_xxyyz_zz = pbuffer.data(idx_hd + 47);

    auto tg_xxyzz_xx = pbuffer.data(idx_hd + 48);

    auto tg_xxyzz_xz = pbuffer.data(idx_hd + 50);

    auto tg_xxyzz_yy = pbuffer.data(idx_hd + 51);

    auto tg_xxyzz_yz = pbuffer.data(idx_hd + 52);

    auto tg_xxzzz_xx = pbuffer.data(idx_hd + 54);

    auto tg_xxzzz_xy = pbuffer.data(idx_hd + 55);

    auto tg_xxzzz_xz = pbuffer.data(idx_hd + 56);

    auto tg_xxzzz_yy = pbuffer.data(idx_hd + 57);

    auto tg_xxzzz_yz = pbuffer.data(idx_hd + 58);

    auto tg_xxzzz_zz = pbuffer.data(idx_hd + 59);

    auto tg_xyyyy_xx = pbuffer.data(idx_hd + 60);

    auto tg_xyyyy_xy = pbuffer.data(idx_hd + 61);

    auto tg_xyyyy_yy = pbuffer.data(idx_hd + 63);

    auto tg_xyyyy_yz = pbuffer.data(idx_hd + 64);

    auto tg_xyyyy_zz = pbuffer.data(idx_hd + 65);

    auto tg_xyyyz_yz = pbuffer.data(idx_hd + 70);

    auto tg_xyyyz_zz = pbuffer.data(idx_hd + 71);

    auto tg_xyyzz_yy = pbuffer.data(idx_hd + 75);

    auto tg_xyyzz_yz = pbuffer.data(idx_hd + 76);

    auto tg_xyyzz_zz = pbuffer.data(idx_hd + 77);

    auto tg_xyzzz_yy = pbuffer.data(idx_hd + 81);

    auto tg_xyzzz_yz = pbuffer.data(idx_hd + 82);

    auto tg_xzzzz_xx = pbuffer.data(idx_hd + 84);

    auto tg_xzzzz_xz = pbuffer.data(idx_hd + 86);

    auto tg_xzzzz_yy = pbuffer.data(idx_hd + 87);

    auto tg_xzzzz_yz = pbuffer.data(idx_hd + 88);

    auto tg_xzzzz_zz = pbuffer.data(idx_hd + 89);

    auto tg_yyyyy_xx = pbuffer.data(idx_hd + 90);

    auto tg_yyyyy_xy = pbuffer.data(idx_hd + 91);

    auto tg_yyyyy_xz = pbuffer.data(idx_hd + 92);

    auto tg_yyyyy_yy = pbuffer.data(idx_hd + 93);

    auto tg_yyyyy_yz = pbuffer.data(idx_hd + 94);

    auto tg_yyyyy_zz = pbuffer.data(idx_hd + 95);

    auto tg_yyyyz_xy = pbuffer.data(idx_hd + 97);

    auto tg_yyyyz_xz = pbuffer.data(idx_hd + 98);

    auto tg_yyyyz_yy = pbuffer.data(idx_hd + 99);

    auto tg_yyyyz_yz = pbuffer.data(idx_hd + 100);

    auto tg_yyyyz_zz = pbuffer.data(idx_hd + 101);

    auto tg_yyyzz_xx = pbuffer.data(idx_hd + 102);

    auto tg_yyyzz_xy = pbuffer.data(idx_hd + 103);

    auto tg_yyyzz_xz = pbuffer.data(idx_hd + 104);

    auto tg_yyyzz_yy = pbuffer.data(idx_hd + 105);

    auto tg_yyyzz_yz = pbuffer.data(idx_hd + 106);

    auto tg_yyyzz_zz = pbuffer.data(idx_hd + 107);

    auto tg_yyzzz_xx = pbuffer.data(idx_hd + 108);

    auto tg_yyzzz_xy = pbuffer.data(idx_hd + 109);

    auto tg_yyzzz_xz = pbuffer.data(idx_hd + 110);

    auto tg_yyzzz_yy = pbuffer.data(idx_hd + 111);

    auto tg_yyzzz_yz = pbuffer.data(idx_hd + 112);

    auto tg_yyzzz_zz = pbuffer.data(idx_hd + 113);

    auto tg_yzzzz_xx = pbuffer.data(idx_hd + 114);

    auto tg_yzzzz_xy = pbuffer.data(idx_hd + 115);

    auto tg_yzzzz_xz = pbuffer.data(idx_hd + 116);

    auto tg_yzzzz_yy = pbuffer.data(idx_hd + 117);

    auto tg_yzzzz_yz = pbuffer.data(idx_hd + 118);

    auto tg_yzzzz_zz = pbuffer.data(idx_hd + 119);

    auto tg_zzzzz_xx = pbuffer.data(idx_hd + 120);

    auto tg_zzzzz_xy = pbuffer.data(idx_hd + 121);

    auto tg_zzzzz_xz = pbuffer.data(idx_hd + 122);

    auto tg_zzzzz_yy = pbuffer.data(idx_hd + 123);

    auto tg_zzzzz_yz = pbuffer.data(idx_hd + 124);

    auto tg_zzzzz_zz = pbuffer.data(idx_hd + 125);

    // Set up components of targeted buffer : ID

    auto tg_xxxxxx_xx = pbuffer.data(idx_id);

    auto tg_xxxxxx_xy = pbuffer.data(idx_id + 1);

    auto tg_xxxxxx_xz = pbuffer.data(idx_id + 2);

    auto tg_xxxxxx_yy = pbuffer.data(idx_id + 3);

    auto tg_xxxxxx_yz = pbuffer.data(idx_id + 4);

    auto tg_xxxxxx_zz = pbuffer.data(idx_id + 5);

    auto tg_xxxxxy_xx = pbuffer.data(idx_id + 6);

    auto tg_xxxxxy_xy = pbuffer.data(idx_id + 7);

    auto tg_xxxxxy_xz = pbuffer.data(idx_id + 8);

    auto tg_xxxxxy_yy = pbuffer.data(idx_id + 9);

    auto tg_xxxxxy_yz = pbuffer.data(idx_id + 10);

    auto tg_xxxxxy_zz = pbuffer.data(idx_id + 11);

    auto tg_xxxxxz_xx = pbuffer.data(idx_id + 12);

    auto tg_xxxxxz_xy = pbuffer.data(idx_id + 13);

    auto tg_xxxxxz_xz = pbuffer.data(idx_id + 14);

    auto tg_xxxxxz_yy = pbuffer.data(idx_id + 15);

    auto tg_xxxxxz_yz = pbuffer.data(idx_id + 16);

    auto tg_xxxxxz_zz = pbuffer.data(idx_id + 17);

    auto tg_xxxxyy_xx = pbuffer.data(idx_id + 18);

    auto tg_xxxxyy_xy = pbuffer.data(idx_id + 19);

    auto tg_xxxxyy_xz = pbuffer.data(idx_id + 20);

    auto tg_xxxxyy_yy = pbuffer.data(idx_id + 21);

    auto tg_xxxxyy_yz = pbuffer.data(idx_id + 22);

    auto tg_xxxxyy_zz = pbuffer.data(idx_id + 23);

    auto tg_xxxxyz_xx = pbuffer.data(idx_id + 24);

    auto tg_xxxxyz_xy = pbuffer.data(idx_id + 25);

    auto tg_xxxxyz_xz = pbuffer.data(idx_id + 26);

    auto tg_xxxxyz_yy = pbuffer.data(idx_id + 27);

    auto tg_xxxxyz_yz = pbuffer.data(idx_id + 28);

    auto tg_xxxxyz_zz = pbuffer.data(idx_id + 29);

    auto tg_xxxxzz_xx = pbuffer.data(idx_id + 30);

    auto tg_xxxxzz_xy = pbuffer.data(idx_id + 31);

    auto tg_xxxxzz_xz = pbuffer.data(idx_id + 32);

    auto tg_xxxxzz_yy = pbuffer.data(idx_id + 33);

    auto tg_xxxxzz_yz = pbuffer.data(idx_id + 34);

    auto tg_xxxxzz_zz = pbuffer.data(idx_id + 35);

    auto tg_xxxyyy_xx = pbuffer.data(idx_id + 36);

    auto tg_xxxyyy_xy = pbuffer.data(idx_id + 37);

    auto tg_xxxyyy_xz = pbuffer.data(idx_id + 38);

    auto tg_xxxyyy_yy = pbuffer.data(idx_id + 39);

    auto tg_xxxyyy_yz = pbuffer.data(idx_id + 40);

    auto tg_xxxyyy_zz = pbuffer.data(idx_id + 41);

    auto tg_xxxyyz_xx = pbuffer.data(idx_id + 42);

    auto tg_xxxyyz_xy = pbuffer.data(idx_id + 43);

    auto tg_xxxyyz_xz = pbuffer.data(idx_id + 44);

    auto tg_xxxyyz_yy = pbuffer.data(idx_id + 45);

    auto tg_xxxyyz_yz = pbuffer.data(idx_id + 46);

    auto tg_xxxyyz_zz = pbuffer.data(idx_id + 47);

    auto tg_xxxyzz_xx = pbuffer.data(idx_id + 48);

    auto tg_xxxyzz_xy = pbuffer.data(idx_id + 49);

    auto tg_xxxyzz_xz = pbuffer.data(idx_id + 50);

    auto tg_xxxyzz_yy = pbuffer.data(idx_id + 51);

    auto tg_xxxyzz_yz = pbuffer.data(idx_id + 52);

    auto tg_xxxyzz_zz = pbuffer.data(idx_id + 53);

    auto tg_xxxzzz_xx = pbuffer.data(idx_id + 54);

    auto tg_xxxzzz_xy = pbuffer.data(idx_id + 55);

    auto tg_xxxzzz_xz = pbuffer.data(idx_id + 56);

    auto tg_xxxzzz_yy = pbuffer.data(idx_id + 57);

    auto tg_xxxzzz_yz = pbuffer.data(idx_id + 58);

    auto tg_xxxzzz_zz = pbuffer.data(idx_id + 59);

    auto tg_xxyyyy_xx = pbuffer.data(idx_id + 60);

    auto tg_xxyyyy_xy = pbuffer.data(idx_id + 61);

    auto tg_xxyyyy_xz = pbuffer.data(idx_id + 62);

    auto tg_xxyyyy_yy = pbuffer.data(idx_id + 63);

    auto tg_xxyyyy_yz = pbuffer.data(idx_id + 64);

    auto tg_xxyyyy_zz = pbuffer.data(idx_id + 65);

    auto tg_xxyyyz_xx = pbuffer.data(idx_id + 66);

    auto tg_xxyyyz_xy = pbuffer.data(idx_id + 67);

    auto tg_xxyyyz_xz = pbuffer.data(idx_id + 68);

    auto tg_xxyyyz_yy = pbuffer.data(idx_id + 69);

    auto tg_xxyyyz_yz = pbuffer.data(idx_id + 70);

    auto tg_xxyyyz_zz = pbuffer.data(idx_id + 71);

    auto tg_xxyyzz_xx = pbuffer.data(idx_id + 72);

    auto tg_xxyyzz_xy = pbuffer.data(idx_id + 73);

    auto tg_xxyyzz_xz = pbuffer.data(idx_id + 74);

    auto tg_xxyyzz_yy = pbuffer.data(idx_id + 75);

    auto tg_xxyyzz_yz = pbuffer.data(idx_id + 76);

    auto tg_xxyyzz_zz = pbuffer.data(idx_id + 77);

    auto tg_xxyzzz_xx = pbuffer.data(idx_id + 78);

    auto tg_xxyzzz_xy = pbuffer.data(idx_id + 79);

    auto tg_xxyzzz_xz = pbuffer.data(idx_id + 80);

    auto tg_xxyzzz_yy = pbuffer.data(idx_id + 81);

    auto tg_xxyzzz_yz = pbuffer.data(idx_id + 82);

    auto tg_xxyzzz_zz = pbuffer.data(idx_id + 83);

    auto tg_xxzzzz_xx = pbuffer.data(idx_id + 84);

    auto tg_xxzzzz_xy = pbuffer.data(idx_id + 85);

    auto tg_xxzzzz_xz = pbuffer.data(idx_id + 86);

    auto tg_xxzzzz_yy = pbuffer.data(idx_id + 87);

    auto tg_xxzzzz_yz = pbuffer.data(idx_id + 88);

    auto tg_xxzzzz_zz = pbuffer.data(idx_id + 89);

    auto tg_xyyyyy_xx = pbuffer.data(idx_id + 90);

    auto tg_xyyyyy_xy = pbuffer.data(idx_id + 91);

    auto tg_xyyyyy_xz = pbuffer.data(idx_id + 92);

    auto tg_xyyyyy_yy = pbuffer.data(idx_id + 93);

    auto tg_xyyyyy_yz = pbuffer.data(idx_id + 94);

    auto tg_xyyyyy_zz = pbuffer.data(idx_id + 95);

    auto tg_xyyyyz_xx = pbuffer.data(idx_id + 96);

    auto tg_xyyyyz_xy = pbuffer.data(idx_id + 97);

    auto tg_xyyyyz_xz = pbuffer.data(idx_id + 98);

    auto tg_xyyyyz_yy = pbuffer.data(idx_id + 99);

    auto tg_xyyyyz_yz = pbuffer.data(idx_id + 100);

    auto tg_xyyyyz_zz = pbuffer.data(idx_id + 101);

    auto tg_xyyyzz_xx = pbuffer.data(idx_id + 102);

    auto tg_xyyyzz_xy = pbuffer.data(idx_id + 103);

    auto tg_xyyyzz_xz = pbuffer.data(idx_id + 104);

    auto tg_xyyyzz_yy = pbuffer.data(idx_id + 105);

    auto tg_xyyyzz_yz = pbuffer.data(idx_id + 106);

    auto tg_xyyyzz_zz = pbuffer.data(idx_id + 107);

    auto tg_xyyzzz_xx = pbuffer.data(idx_id + 108);

    auto tg_xyyzzz_xy = pbuffer.data(idx_id + 109);

    auto tg_xyyzzz_xz = pbuffer.data(idx_id + 110);

    auto tg_xyyzzz_yy = pbuffer.data(idx_id + 111);

    auto tg_xyyzzz_yz = pbuffer.data(idx_id + 112);

    auto tg_xyyzzz_zz = pbuffer.data(idx_id + 113);

    auto tg_xyzzzz_xx = pbuffer.data(idx_id + 114);

    auto tg_xyzzzz_xy = pbuffer.data(idx_id + 115);

    auto tg_xyzzzz_xz = pbuffer.data(idx_id + 116);

    auto tg_xyzzzz_yy = pbuffer.data(idx_id + 117);

    auto tg_xyzzzz_yz = pbuffer.data(idx_id + 118);

    auto tg_xyzzzz_zz = pbuffer.data(idx_id + 119);

    auto tg_xzzzzz_xx = pbuffer.data(idx_id + 120);

    auto tg_xzzzzz_xy = pbuffer.data(idx_id + 121);

    auto tg_xzzzzz_xz = pbuffer.data(idx_id + 122);

    auto tg_xzzzzz_yy = pbuffer.data(idx_id + 123);

    auto tg_xzzzzz_yz = pbuffer.data(idx_id + 124);

    auto tg_xzzzzz_zz = pbuffer.data(idx_id + 125);

    auto tg_yyyyyy_xx = pbuffer.data(idx_id + 126);

    auto tg_yyyyyy_xy = pbuffer.data(idx_id + 127);

    auto tg_yyyyyy_xz = pbuffer.data(idx_id + 128);

    auto tg_yyyyyy_yy = pbuffer.data(idx_id + 129);

    auto tg_yyyyyy_yz = pbuffer.data(idx_id + 130);

    auto tg_yyyyyy_zz = pbuffer.data(idx_id + 131);

    auto tg_yyyyyz_xx = pbuffer.data(idx_id + 132);

    auto tg_yyyyyz_xy = pbuffer.data(idx_id + 133);

    auto tg_yyyyyz_xz = pbuffer.data(idx_id + 134);

    auto tg_yyyyyz_yy = pbuffer.data(idx_id + 135);

    auto tg_yyyyyz_yz = pbuffer.data(idx_id + 136);

    auto tg_yyyyyz_zz = pbuffer.data(idx_id + 137);

    auto tg_yyyyzz_xx = pbuffer.data(idx_id + 138);

    auto tg_yyyyzz_xy = pbuffer.data(idx_id + 139);

    auto tg_yyyyzz_xz = pbuffer.data(idx_id + 140);

    auto tg_yyyyzz_yy = pbuffer.data(idx_id + 141);

    auto tg_yyyyzz_yz = pbuffer.data(idx_id + 142);

    auto tg_yyyyzz_zz = pbuffer.data(idx_id + 143);

    auto tg_yyyzzz_xx = pbuffer.data(idx_id + 144);

    auto tg_yyyzzz_xy = pbuffer.data(idx_id + 145);

    auto tg_yyyzzz_xz = pbuffer.data(idx_id + 146);

    auto tg_yyyzzz_yy = pbuffer.data(idx_id + 147);

    auto tg_yyyzzz_yz = pbuffer.data(idx_id + 148);

    auto tg_yyyzzz_zz = pbuffer.data(idx_id + 149);

    auto tg_yyzzzz_xx = pbuffer.data(idx_id + 150);

    auto tg_yyzzzz_xy = pbuffer.data(idx_id + 151);

    auto tg_yyzzzz_xz = pbuffer.data(idx_id + 152);

    auto tg_yyzzzz_yy = pbuffer.data(idx_id + 153);

    auto tg_yyzzzz_yz = pbuffer.data(idx_id + 154);

    auto tg_yyzzzz_zz = pbuffer.data(idx_id + 155);

    auto tg_yzzzzz_xx = pbuffer.data(idx_id + 156);

    auto tg_yzzzzz_xy = pbuffer.data(idx_id + 157);

    auto tg_yzzzzz_xz = pbuffer.data(idx_id + 158);

    auto tg_yzzzzz_yy = pbuffer.data(idx_id + 159);

    auto tg_yzzzzz_yz = pbuffer.data(idx_id + 160);

    auto tg_yzzzzz_zz = pbuffer.data(idx_id + 161);

    auto tg_zzzzzz_xx = pbuffer.data(idx_id + 162);

    auto tg_zzzzzz_xy = pbuffer.data(idx_id + 163);

    auto tg_zzzzzz_xz = pbuffer.data(idx_id + 164);

    auto tg_zzzzzz_yy = pbuffer.data(idx_id + 165);

    auto tg_zzzzzz_yz = pbuffer.data(idx_id + 166);

    auto tg_zzzzzz_zz = pbuffer.data(idx_id + 167);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_xxxx_xx, tg_xxxx_xy, tg_xxxx_xz, tg_xxxx_yy, tg_xxxx_yz, tg_xxxx_zz, tg_xxxxx_x, tg_xxxxx_xx, tg_xxxxx_xy, tg_xxxxx_xz, tg_xxxxx_y, tg_xxxxx_yy, tg_xxxxx_yz, tg_xxxxx_z, tg_xxxxx_zz, tg_xxxxxx_xx, tg_xxxxxx_xy, tg_xxxxxx_xz, tg_xxxxxx_yy, tg_xxxxxx_yz, tg_xxxxxx_zz, tg_xxxxxy_xx, tg_xxxxxy_xy, tg_xxxxxy_xz, tg_xxxxxy_yy, tg_xxxxxy_yz, tg_xxxxxy_zz, tg_xxxxxz_xx, tg_xxxxxz_xy, tg_xxxxxz_xz, tg_xxxxxz_yy, tg_xxxxxz_yz, tg_xxxxxz_zz, tg_xxxxy_xx, tg_xxxxy_xy, tg_xxxxy_xz, tg_xxxxy_yy, tg_xxxxy_yz, tg_xxxxyy_xx, tg_xxxxyy_xy, tg_xxxxyy_xz, tg_xxxxyy_yy, tg_xxxxyy_yz, tg_xxxxyy_zz, tg_xxxxyz_xx, tg_xxxxyz_xy, tg_xxxxyz_xz, tg_xxxxyz_yy, tg_xxxxyz_yz, tg_xxxxyz_zz, tg_xxxxz_xx, tg_xxxxz_xy, tg_xxxxz_xz, tg_xxxxz_yz, tg_xxxxz_zz, tg_xxxxzz_xx, tg_xxxxzz_xy, tg_xxxxzz_xz, tg_xxxxzz_yy, tg_xxxxzz_yz, tg_xxxxzz_zz, tg_xxxy_xx, tg_xxxy_xz, tg_xxxy_yy, tg_xxxy_yz, tg_xxxyy_xx, tg_xxxyy_xy, tg_xxxyy_xz, tg_xxxyy_y, tg_xxxyy_yy, tg_xxxyy_yz, tg_xxxyy_zz, tg_xxxyyy_xx, tg_xxxyyy_xy, tg_xxxyyy_xz, tg_xxxyyy_yy, tg_xxxyyy_yz, tg_xxxyyy_zz, tg_xxxyyz_xx, tg_xxxyyz_xy, tg_xxxyyz_xz, tg_xxxyyz_yy, tg_xxxyyz_yz, tg_xxxyyz_zz, tg_xxxyz_xz, tg_xxxyz_yz, tg_xxxyzz_xx, tg_xxxyzz_xy, tg_xxxyzz_xz, tg_xxxyzz_yy, tg_xxxyzz_yz, tg_xxxyzz_zz, tg_xxxz_xx, tg_xxxz_xy, tg_xxxz_xz, tg_xxxz_yz, tg_xxxz_zz, tg_xxxzz_x, tg_xxxzz_xx, tg_xxxzz_xy, tg_xxxzz_xz, tg_xxxzz_yy, tg_xxxzz_yz, tg_xxxzz_z, tg_xxxzz_zz, tg_xxxzzz_xx, tg_xxxzzz_xy, tg_xxxzzz_xz, tg_xxxzzz_yy, tg_xxxzzz_yz, tg_xxxzzz_zz, tg_xxyy_xx, tg_xxyy_xy, tg_xxyy_xz, tg_xxyy_yy, tg_xxyy_yz, tg_xxyy_zz, tg_xxyyy_xx, tg_xxyyy_xy, tg_xxyyy_xz, tg_xxyyy_y, tg_xxyyy_yy, tg_xxyyy_yz, tg_xxyyy_zz, tg_xxyyyy_xx, tg_xxyyyy_xy, tg_xxyyyy_xz, tg_xxyyyy_yy, tg_xxyyyy_yz, tg_xxyyyy_zz, tg_xxyyyz_xx, tg_xxyyyz_xy, tg_xxyyyz_xz, tg_xxyyyz_yy, tg_xxyyyz_yz, tg_xxyyyz_zz, tg_xxyyz_xy, tg_xxyyz_xz, tg_xxyyz_yz, tg_xxyyz_zz, tg_xxyyzz_xx, tg_xxyyzz_xy, tg_xxyyzz_xz, tg_xxyyzz_yy, tg_xxyyzz_yz, tg_xxyyzz_zz, tg_xxyz_xz, tg_xxyz_yz, tg_xxyzz_xx, tg_xxyzz_xz, tg_xxyzz_yy, tg_xxyzz_yz, tg_xxyzzz_xx, tg_xxyzzz_xy, tg_xxyzzz_xz, tg_xxyzzz_yy, tg_xxyzzz_yz, tg_xxyzzz_zz, tg_xxzz_xx, tg_xxzz_xy, tg_xxzz_xz, tg_xxzz_yy, tg_xxzz_yz, tg_xxzz_zz, tg_xxzzz_x, tg_xxzzz_xx, tg_xxzzz_xy, tg_xxzzz_xz, tg_xxzzz_yy, tg_xxzzz_yz, tg_xxzzz_z, tg_xxzzz_zz, tg_xxzzzz_xx, tg_xxzzzz_xy, tg_xxzzzz_xz, tg_xxzzzz_yy, tg_xxzzzz_yz, tg_xxzzzz_zz, tg_xyyy_xy, tg_xyyy_yy, tg_xyyy_yz, tg_xyyy_zz, tg_xyyyy_xx, tg_xyyyy_xy, tg_xyyyy_y, tg_xyyyy_yy, tg_xyyyy_yz, tg_xyyyy_zz, tg_xyyyyy_xx, tg_xyyyyy_xy, tg_xyyyyy_xz, tg_xyyyyy_yy, tg_xyyyyy_yz, tg_xyyyyy_zz, tg_xyyyyz_xx, tg_xyyyyz_xy, tg_xyyyyz_xz, tg_xyyyyz_yy, tg_xyyyyz_yz, tg_xyyyyz_zz, tg_xyyyz_yz, tg_xyyyz_zz, tg_xyyyzz_xx, tg_xyyyzz_xy, tg_xyyyzz_xz, tg_xyyyzz_yy, tg_xyyyzz_yz, tg_xyyyzz_zz, tg_xyyz_yz, tg_xyyz_zz, tg_xyyzz_yy, tg_xyyzz_yz, tg_xyyzz_zz, tg_xyyzzz_xx, tg_xyyzzz_xy, tg_xyyzzz_xz, tg_xyyzzz_yy, tg_xyyzzz_yz, tg_xyyzzz_zz, tg_xyzz_yy, tg_xyzz_yz, tg_xyzzz_yy, tg_xyzzz_yz, tg_xyzzzz_xx, tg_xyzzzz_xy, tg_xyzzzz_xz, tg_xyzzzz_yy, tg_xyzzzz_yz, tg_xyzzzz_zz, tg_xzzz_xz, tg_xzzz_yy, tg_xzzz_yz, tg_xzzz_zz, tg_xzzzz_xx, tg_xzzzz_xz, tg_xzzzz_yy, tg_xzzzz_yz, tg_xzzzz_z, tg_xzzzz_zz, tg_xzzzzz_xx, tg_xzzzzz_xy, tg_xzzzzz_xz, tg_xzzzzz_yy, tg_xzzzzz_yz, tg_xzzzzz_zz, tg_yyyy_xx, tg_yyyy_xy, tg_yyyy_xz, tg_yyyy_yy, tg_yyyy_yz, tg_yyyy_zz, tg_yyyyy_x, tg_yyyyy_xx, tg_yyyyy_xy, tg_yyyyy_xz, tg_yyyyy_y, tg_yyyyy_yy, tg_yyyyy_yz, tg_yyyyy_z, tg_yyyyy_zz, tg_yyyyyy_xx, tg_yyyyyy_xy, tg_yyyyyy_xz, tg_yyyyyy_yy, tg_yyyyyy_yz, tg_yyyyyy_zz, tg_yyyyyz_xx, tg_yyyyyz_xy, tg_yyyyyz_xz, tg_yyyyyz_yy, tg_yyyyyz_yz, tg_yyyyyz_zz, tg_yyyyz_xy, tg_yyyyz_xz, tg_yyyyz_yy, tg_yyyyz_yz, tg_yyyyz_z, tg_yyyyz_zz, tg_yyyyzz_xx, tg_yyyyzz_xy, tg_yyyyzz_xz, tg_yyyyzz_yy, tg_yyyyzz_yz, tg_yyyyzz_zz, tg_yyyz_xy, tg_yyyz_xz, tg_yyyz_yy, tg_yyyz_yz, tg_yyyz_zz, tg_yyyzz_x, tg_yyyzz_xx, tg_yyyzz_xy, tg_yyyzz_xz, tg_yyyzz_y, tg_yyyzz_yy, tg_yyyzz_yz, tg_yyyzz_z, tg_yyyzz_zz, tg_yyyzzz_xx, tg_yyyzzz_xy, tg_yyyzzz_xz, tg_yyyzzz_yy, tg_yyyzzz_yz, tg_yyyzzz_zz, tg_yyzz_xx, tg_yyzz_xy, tg_yyzz_xz, tg_yyzz_yy, tg_yyzz_yz, tg_yyzz_zz, tg_yyzzz_x, tg_yyzzz_xx, tg_yyzzz_xy, tg_yyzzz_xz, tg_yyzzz_y, tg_yyzzz_yy, tg_yyzzz_yz, tg_yyzzz_z, tg_yyzzz_zz, tg_yyzzzz_xx, tg_yyzzzz_xy, tg_yyzzzz_xz, tg_yyzzzz_yy, tg_yyzzzz_yz, tg_yyzzzz_zz, tg_yzzz_xx, tg_yzzz_xz, tg_yzzz_yy, tg_yzzz_yz, tg_yzzz_zz, tg_yzzzz_xx, tg_yzzzz_xy, tg_yzzzz_xz, tg_yzzzz_y, tg_yzzzz_yy, tg_yzzzz_yz, tg_yzzzz_z, tg_yzzzz_zz, tg_yzzzzz_xx, tg_yzzzzz_xy, tg_yzzzzz_xz, tg_yzzzzz_yy, tg_yzzzzz_yz, tg_yzzzzz_zz, tg_zzzz_xx, tg_zzzz_xy, tg_zzzz_xz, tg_zzzz_yy, tg_zzzz_yz, tg_zzzz_zz, tg_zzzzz_x, tg_zzzzz_xx, tg_zzzzz_xy, tg_zzzzz_xz, tg_zzzzz_y, tg_zzzzz_yy, tg_zzzzz_yz, tg_zzzzz_z, tg_zzzzz_zz, tg_zzzzzz_xx, tg_zzzzzz_xy, tg_zzzzzz_xz, tg_zzzzzz_yy, tg_zzzzzz_yz, tg_zzzzzz_zz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxxxxx_xx[i] = 5.0 * tg_xxxx_xx[i] * fxi[i] + 2.0 * tg_xxxxx_x[i] * fxi[i] + tg_xxxxx_xx[i] * ra_x[i];

        tg_xxxxxx_xy[i] = 5.0 * tg_xxxx_xy[i] * fxi[i] + tg_xxxxx_y[i] * fxi[i] + tg_xxxxx_xy[i] * ra_x[i];

        tg_xxxxxx_xz[i] = 5.0 * tg_xxxx_xz[i] * fxi[i] + tg_xxxxx_z[i] * fxi[i] + tg_xxxxx_xz[i] * ra_x[i];

        tg_xxxxxx_yy[i] = 5.0 * tg_xxxx_yy[i] * fxi[i] + tg_xxxxx_yy[i] * ra_x[i];

        tg_xxxxxx_yz[i] = 5.0 * tg_xxxx_yz[i] * fxi[i] + tg_xxxxx_yz[i] * ra_x[i];

        tg_xxxxxx_zz[i] = 5.0 * tg_xxxx_zz[i] * fxi[i] + tg_xxxxx_zz[i] * ra_x[i];

        tg_xxxxxy_xx[i] = tg_xxxxx_xx[i] * ra_y[i];

        tg_xxxxxy_xy[i] = tg_xxxxx_x[i] * fxi[i] + tg_xxxxx_xy[i] * ra_y[i];

        tg_xxxxxy_xz[i] = tg_xxxxx_xz[i] * ra_y[i];

        tg_xxxxxy_yy[i] = 4.0 * tg_xxxy_yy[i] * fxi[i] + tg_xxxxy_yy[i] * ra_x[i];

        tg_xxxxxy_yz[i] = 4.0 * tg_xxxy_yz[i] * fxi[i] + tg_xxxxy_yz[i] * ra_x[i];

        tg_xxxxxy_zz[i] = tg_xxxxx_zz[i] * ra_y[i];

        tg_xxxxxz_xx[i] = tg_xxxxx_xx[i] * ra_z[i];

        tg_xxxxxz_xy[i] = tg_xxxxx_xy[i] * ra_z[i];

        tg_xxxxxz_xz[i] = tg_xxxxx_x[i] * fxi[i] + tg_xxxxx_xz[i] * ra_z[i];

        tg_xxxxxz_yy[i] = tg_xxxxx_yy[i] * ra_z[i];

        tg_xxxxxz_yz[i] = 4.0 * tg_xxxz_yz[i] * fxi[i] + tg_xxxxz_yz[i] * ra_x[i];

        tg_xxxxxz_zz[i] = 4.0 * tg_xxxz_zz[i] * fxi[i] + tg_xxxxz_zz[i] * ra_x[i];

        tg_xxxxyy_xx[i] = tg_xxxx_xx[i] * fxi[i] + tg_xxxxy_xx[i] * ra_y[i];

        tg_xxxxyy_xy[i] = 3.0 * tg_xxyy_xy[i] * fxi[i] + tg_xxxyy_y[i] * fxi[i] + tg_xxxyy_xy[i] * ra_x[i];

        tg_xxxxyy_xz[i] = tg_xxxx_xz[i] * fxi[i] + tg_xxxxy_xz[i] * ra_y[i];

        tg_xxxxyy_yy[i] = 3.0 * tg_xxyy_yy[i] * fxi[i] + tg_xxxyy_yy[i] * ra_x[i];

        tg_xxxxyy_yz[i] = 3.0 * tg_xxyy_yz[i] * fxi[i] + tg_xxxyy_yz[i] * ra_x[i];

        tg_xxxxyy_zz[i] = 3.0 * tg_xxyy_zz[i] * fxi[i] + tg_xxxyy_zz[i] * ra_x[i];

        tg_xxxxyz_xx[i] = tg_xxxxz_xx[i] * ra_y[i];

        tg_xxxxyz_xy[i] = tg_xxxxy_xy[i] * ra_z[i];

        tg_xxxxyz_xz[i] = tg_xxxxz_xz[i] * ra_y[i];

        tg_xxxxyz_yy[i] = tg_xxxxy_yy[i] * ra_z[i];

        tg_xxxxyz_yz[i] = 3.0 * tg_xxyz_yz[i] * fxi[i] + tg_xxxyz_yz[i] * ra_x[i];

        tg_xxxxyz_zz[i] = tg_xxxxz_zz[i] * ra_y[i];

        tg_xxxxzz_xx[i] = tg_xxxx_xx[i] * fxi[i] + tg_xxxxz_xx[i] * ra_z[i];

        tg_xxxxzz_xy[i] = tg_xxxx_xy[i] * fxi[i] + tg_xxxxz_xy[i] * ra_z[i];

        tg_xxxxzz_xz[i] = 3.0 * tg_xxzz_xz[i] * fxi[i] + tg_xxxzz_z[i] * fxi[i] + tg_xxxzz_xz[i] * ra_x[i];

        tg_xxxxzz_yy[i] = 3.0 * tg_xxzz_yy[i] * fxi[i] + tg_xxxzz_yy[i] * ra_x[i];

        tg_xxxxzz_yz[i] = 3.0 * tg_xxzz_yz[i] * fxi[i] + tg_xxxzz_yz[i] * ra_x[i];

        tg_xxxxzz_zz[i] = 3.0 * tg_xxzz_zz[i] * fxi[i] + tg_xxxzz_zz[i] * ra_x[i];

        tg_xxxyyy_xx[i] = 2.0 * tg_xxxy_xx[i] * fxi[i] + tg_xxxyy_xx[i] * ra_y[i];

        tg_xxxyyy_xy[i] = 2.0 * tg_xyyy_xy[i] * fxi[i] + tg_xxyyy_y[i] * fxi[i] + tg_xxyyy_xy[i] * ra_x[i];

        tg_xxxyyy_xz[i] = 2.0 * tg_xxxy_xz[i] * fxi[i] + tg_xxxyy_xz[i] * ra_y[i];

        tg_xxxyyy_yy[i] = 2.0 * tg_xyyy_yy[i] * fxi[i] + tg_xxyyy_yy[i] * ra_x[i];

        tg_xxxyyy_yz[i] = 2.0 * tg_xyyy_yz[i] * fxi[i] + tg_xxyyy_yz[i] * ra_x[i];

        tg_xxxyyy_zz[i] = 2.0 * tg_xyyy_zz[i] * fxi[i] + tg_xxyyy_zz[i] * ra_x[i];

        tg_xxxyyz_xx[i] = tg_xxxyy_xx[i] * ra_z[i];

        tg_xxxyyz_xy[i] = tg_xxxyy_xy[i] * ra_z[i];

        tg_xxxyyz_xz[i] = tg_xxxz_xz[i] * fxi[i] + tg_xxxyz_xz[i] * ra_y[i];

        tg_xxxyyz_yy[i] = tg_xxxyy_yy[i] * ra_z[i];

        tg_xxxyyz_yz[i] = 2.0 * tg_xyyz_yz[i] * fxi[i] + tg_xxyyz_yz[i] * ra_x[i];

        tg_xxxyyz_zz[i] = 2.0 * tg_xyyz_zz[i] * fxi[i] + tg_xxyyz_zz[i] * ra_x[i];

        tg_xxxyzz_xx[i] = tg_xxxzz_xx[i] * ra_y[i];

        tg_xxxyzz_xy[i] = tg_xxxzz_x[i] * fxi[i] + tg_xxxzz_xy[i] * ra_y[i];

        tg_xxxyzz_xz[i] = tg_xxxzz_xz[i] * ra_y[i];

        tg_xxxyzz_yy[i] = 2.0 * tg_xyzz_yy[i] * fxi[i] + tg_xxyzz_yy[i] * ra_x[i];

        tg_xxxyzz_yz[i] = 2.0 * tg_xyzz_yz[i] * fxi[i] + tg_xxyzz_yz[i] * ra_x[i];

        tg_xxxyzz_zz[i] = tg_xxxzz_zz[i] * ra_y[i];

        tg_xxxzzz_xx[i] = 2.0 * tg_xxxz_xx[i] * fxi[i] + tg_xxxzz_xx[i] * ra_z[i];

        tg_xxxzzz_xy[i] = 2.0 * tg_xxxz_xy[i] * fxi[i] + tg_xxxzz_xy[i] * ra_z[i];

        tg_xxxzzz_xz[i] = 2.0 * tg_xzzz_xz[i] * fxi[i] + tg_xxzzz_z[i] * fxi[i] + tg_xxzzz_xz[i] * ra_x[i];

        tg_xxxzzz_yy[i] = 2.0 * tg_xzzz_yy[i] * fxi[i] + tg_xxzzz_yy[i] * ra_x[i];

        tg_xxxzzz_yz[i] = 2.0 * tg_xzzz_yz[i] * fxi[i] + tg_xxzzz_yz[i] * ra_x[i];

        tg_xxxzzz_zz[i] = 2.0 * tg_xzzz_zz[i] * fxi[i] + tg_xxzzz_zz[i] * ra_x[i];

        tg_xxyyyy_xx[i] = 3.0 * tg_xxyy_xx[i] * fxi[i] + tg_xxyyy_xx[i] * ra_y[i];

        tg_xxyyyy_xy[i] = tg_yyyy_xy[i] * fxi[i] + tg_xyyyy_y[i] * fxi[i] + tg_xyyyy_xy[i] * ra_x[i];

        tg_xxyyyy_xz[i] = 3.0 * tg_xxyy_xz[i] * fxi[i] + tg_xxyyy_xz[i] * ra_y[i];

        tg_xxyyyy_yy[i] = tg_yyyy_yy[i] * fxi[i] + tg_xyyyy_yy[i] * ra_x[i];

        tg_xxyyyy_yz[i] = tg_yyyy_yz[i] * fxi[i] + tg_xyyyy_yz[i] * ra_x[i];

        tg_xxyyyy_zz[i] = tg_yyyy_zz[i] * fxi[i] + tg_xyyyy_zz[i] * ra_x[i];

        tg_xxyyyz_xx[i] = tg_xxyyy_xx[i] * ra_z[i];

        tg_xxyyyz_xy[i] = tg_xxyyy_xy[i] * ra_z[i];

        tg_xxyyyz_xz[i] = 2.0 * tg_xxyz_xz[i] * fxi[i] + tg_xxyyz_xz[i] * ra_y[i];

        tg_xxyyyz_yy[i] = tg_xxyyy_yy[i] * ra_z[i];

        tg_xxyyyz_yz[i] = tg_yyyz_yz[i] * fxi[i] + tg_xyyyz_yz[i] * ra_x[i];

        tg_xxyyyz_zz[i] = tg_yyyz_zz[i] * fxi[i] + tg_xyyyz_zz[i] * ra_x[i];

        tg_xxyyzz_xx[i] = tg_xxzz_xx[i] * fxi[i] + tg_xxyzz_xx[i] * ra_y[i];

        tg_xxyyzz_xy[i] = tg_xxyy_xy[i] * fxi[i] + tg_xxyyz_xy[i] * ra_z[i];

        tg_xxyyzz_xz[i] = tg_xxzz_xz[i] * fxi[i] + tg_xxyzz_xz[i] * ra_y[i];

        tg_xxyyzz_yy[i] = tg_yyzz_yy[i] * fxi[i] + tg_xyyzz_yy[i] * ra_x[i];

        tg_xxyyzz_yz[i] = tg_yyzz_yz[i] * fxi[i] + tg_xyyzz_yz[i] * ra_x[i];

        tg_xxyyzz_zz[i] = tg_yyzz_zz[i] * fxi[i] + tg_xyyzz_zz[i] * ra_x[i];

        tg_xxyzzz_xx[i] = tg_xxzzz_xx[i] * ra_y[i];

        tg_xxyzzz_xy[i] = tg_xxzzz_x[i] * fxi[i] + tg_xxzzz_xy[i] * ra_y[i];

        tg_xxyzzz_xz[i] = tg_xxzzz_xz[i] * ra_y[i];

        tg_xxyzzz_yy[i] = tg_yzzz_yy[i] * fxi[i] + tg_xyzzz_yy[i] * ra_x[i];

        tg_xxyzzz_yz[i] = tg_yzzz_yz[i] * fxi[i] + tg_xyzzz_yz[i] * ra_x[i];

        tg_xxyzzz_zz[i] = tg_xxzzz_zz[i] * ra_y[i];

        tg_xxzzzz_xx[i] = 3.0 * tg_xxzz_xx[i] * fxi[i] + tg_xxzzz_xx[i] * ra_z[i];

        tg_xxzzzz_xy[i] = 3.0 * tg_xxzz_xy[i] * fxi[i] + tg_xxzzz_xy[i] * ra_z[i];

        tg_xxzzzz_xz[i] = tg_zzzz_xz[i] * fxi[i] + tg_xzzzz_z[i] * fxi[i] + tg_xzzzz_xz[i] * ra_x[i];

        tg_xxzzzz_yy[i] = tg_zzzz_yy[i] * fxi[i] + tg_xzzzz_yy[i] * ra_x[i];

        tg_xxzzzz_yz[i] = tg_zzzz_yz[i] * fxi[i] + tg_xzzzz_yz[i] * ra_x[i];

        tg_xxzzzz_zz[i] = tg_zzzz_zz[i] * fxi[i] + tg_xzzzz_zz[i] * ra_x[i];

        tg_xyyyyy_xx[i] = 2.0 * tg_yyyyy_x[i] * fxi[i] + tg_yyyyy_xx[i] * ra_x[i];

        tg_xyyyyy_xy[i] = tg_yyyyy_y[i] * fxi[i] + tg_yyyyy_xy[i] * ra_x[i];

        tg_xyyyyy_xz[i] = tg_yyyyy_z[i] * fxi[i] + tg_yyyyy_xz[i] * ra_x[i];

        tg_xyyyyy_yy[i] = tg_yyyyy_yy[i] * ra_x[i];

        tg_xyyyyy_yz[i] = tg_yyyyy_yz[i] * ra_x[i];

        tg_xyyyyy_zz[i] = tg_yyyyy_zz[i] * ra_x[i];

        tg_xyyyyz_xx[i] = tg_xyyyy_xx[i] * ra_z[i];

        tg_xyyyyz_xy[i] = tg_xyyyy_xy[i] * ra_z[i];

        tg_xyyyyz_xz[i] = tg_yyyyz_z[i] * fxi[i] + tg_yyyyz_xz[i] * ra_x[i];

        tg_xyyyyz_yy[i] = tg_yyyyz_yy[i] * ra_x[i];

        tg_xyyyyz_yz[i] = tg_yyyyz_yz[i] * ra_x[i];

        tg_xyyyyz_zz[i] = tg_yyyyz_zz[i] * ra_x[i];

        tg_xyyyzz_xx[i] = 2.0 * tg_yyyzz_x[i] * fxi[i] + tg_yyyzz_xx[i] * ra_x[i];

        tg_xyyyzz_xy[i] = tg_yyyzz_y[i] * fxi[i] + tg_yyyzz_xy[i] * ra_x[i];

        tg_xyyyzz_xz[i] = tg_yyyzz_z[i] * fxi[i] + tg_yyyzz_xz[i] * ra_x[i];

        tg_xyyyzz_yy[i] = tg_yyyzz_yy[i] * ra_x[i];

        tg_xyyyzz_yz[i] = tg_yyyzz_yz[i] * ra_x[i];

        tg_xyyyzz_zz[i] = tg_yyyzz_zz[i] * ra_x[i];

        tg_xyyzzz_xx[i] = 2.0 * tg_yyzzz_x[i] * fxi[i] + tg_yyzzz_xx[i] * ra_x[i];

        tg_xyyzzz_xy[i] = tg_yyzzz_y[i] * fxi[i] + tg_yyzzz_xy[i] * ra_x[i];

        tg_xyyzzz_xz[i] = tg_yyzzz_z[i] * fxi[i] + tg_yyzzz_xz[i] * ra_x[i];

        tg_xyyzzz_yy[i] = tg_yyzzz_yy[i] * ra_x[i];

        tg_xyyzzz_yz[i] = tg_yyzzz_yz[i] * ra_x[i];

        tg_xyyzzz_zz[i] = tg_yyzzz_zz[i] * ra_x[i];

        tg_xyzzzz_xx[i] = tg_xzzzz_xx[i] * ra_y[i];

        tg_xyzzzz_xy[i] = tg_yzzzz_y[i] * fxi[i] + tg_yzzzz_xy[i] * ra_x[i];

        tg_xyzzzz_xz[i] = tg_xzzzz_xz[i] * ra_y[i];

        tg_xyzzzz_yy[i] = tg_yzzzz_yy[i] * ra_x[i];

        tg_xyzzzz_yz[i] = tg_yzzzz_yz[i] * ra_x[i];

        tg_xyzzzz_zz[i] = tg_yzzzz_zz[i] * ra_x[i];

        tg_xzzzzz_xx[i] = 2.0 * tg_zzzzz_x[i] * fxi[i] + tg_zzzzz_xx[i] * ra_x[i];

        tg_xzzzzz_xy[i] = tg_zzzzz_y[i] * fxi[i] + tg_zzzzz_xy[i] * ra_x[i];

        tg_xzzzzz_xz[i] = tg_zzzzz_z[i] * fxi[i] + tg_zzzzz_xz[i] * ra_x[i];

        tg_xzzzzz_yy[i] = tg_zzzzz_yy[i] * ra_x[i];

        tg_xzzzzz_yz[i] = tg_zzzzz_yz[i] * ra_x[i];

        tg_xzzzzz_zz[i] = tg_zzzzz_zz[i] * ra_x[i];

        tg_yyyyyy_xx[i] = 5.0 * tg_yyyy_xx[i] * fxi[i] + tg_yyyyy_xx[i] * ra_y[i];

        tg_yyyyyy_xy[i] = 5.0 * tg_yyyy_xy[i] * fxi[i] + tg_yyyyy_x[i] * fxi[i] + tg_yyyyy_xy[i] * ra_y[i];

        tg_yyyyyy_xz[i] = 5.0 * tg_yyyy_xz[i] * fxi[i] + tg_yyyyy_xz[i] * ra_y[i];

        tg_yyyyyy_yy[i] = 5.0 * tg_yyyy_yy[i] * fxi[i] + 2.0 * tg_yyyyy_y[i] * fxi[i] + tg_yyyyy_yy[i] * ra_y[i];

        tg_yyyyyy_yz[i] = 5.0 * tg_yyyy_yz[i] * fxi[i] + tg_yyyyy_z[i] * fxi[i] + tg_yyyyy_yz[i] * ra_y[i];

        tg_yyyyyy_zz[i] = 5.0 * tg_yyyy_zz[i] * fxi[i] + tg_yyyyy_zz[i] * ra_y[i];

        tg_yyyyyz_xx[i] = tg_yyyyy_xx[i] * ra_z[i];

        tg_yyyyyz_xy[i] = tg_yyyyy_xy[i] * ra_z[i];

        tg_yyyyyz_xz[i] = 4.0 * tg_yyyz_xz[i] * fxi[i] + tg_yyyyz_xz[i] * ra_y[i];

        tg_yyyyyz_yy[i] = tg_yyyyy_yy[i] * ra_z[i];

        tg_yyyyyz_yz[i] = tg_yyyyy_y[i] * fxi[i] + tg_yyyyy_yz[i] * ra_z[i];

        tg_yyyyyz_zz[i] = 4.0 * tg_yyyz_zz[i] * fxi[i] + tg_yyyyz_zz[i] * ra_y[i];

        tg_yyyyzz_xx[i] = 3.0 * tg_yyzz_xx[i] * fxi[i] + tg_yyyzz_xx[i] * ra_y[i];

        tg_yyyyzz_xy[i] = tg_yyyy_xy[i] * fxi[i] + tg_yyyyz_xy[i] * ra_z[i];

        tg_yyyyzz_xz[i] = 3.0 * tg_yyzz_xz[i] * fxi[i] + tg_yyyzz_xz[i] * ra_y[i];

        tg_yyyyzz_yy[i] = tg_yyyy_yy[i] * fxi[i] + tg_yyyyz_yy[i] * ra_z[i];

        tg_yyyyzz_yz[i] = 3.0 * tg_yyzz_yz[i] * fxi[i] + tg_yyyzz_z[i] * fxi[i] + tg_yyyzz_yz[i] * ra_y[i];

        tg_yyyyzz_zz[i] = 3.0 * tg_yyzz_zz[i] * fxi[i] + tg_yyyzz_zz[i] * ra_y[i];

        tg_yyyzzz_xx[i] = 2.0 * tg_yzzz_xx[i] * fxi[i] + tg_yyzzz_xx[i] * ra_y[i];

        tg_yyyzzz_xy[i] = 2.0 * tg_yyyz_xy[i] * fxi[i] + tg_yyyzz_xy[i] * ra_z[i];

        tg_yyyzzz_xz[i] = 2.0 * tg_yzzz_xz[i] * fxi[i] + tg_yyzzz_xz[i] * ra_y[i];

        tg_yyyzzz_yy[i] = 2.0 * tg_yyyz_yy[i] * fxi[i] + tg_yyyzz_yy[i] * ra_z[i];

        tg_yyyzzz_yz[i] = 2.0 * tg_yzzz_yz[i] * fxi[i] + tg_yyzzz_z[i] * fxi[i] + tg_yyzzz_yz[i] * ra_y[i];

        tg_yyyzzz_zz[i] = 2.0 * tg_yzzz_zz[i] * fxi[i] + tg_yyzzz_zz[i] * ra_y[i];

        tg_yyzzzz_xx[i] = tg_zzzz_xx[i] * fxi[i] + tg_yzzzz_xx[i] * ra_y[i];

        tg_yyzzzz_xy[i] = 3.0 * tg_yyzz_xy[i] * fxi[i] + tg_yyzzz_xy[i] * ra_z[i];

        tg_yyzzzz_xz[i] = tg_zzzz_xz[i] * fxi[i] + tg_yzzzz_xz[i] * ra_y[i];

        tg_yyzzzz_yy[i] = 3.0 * tg_yyzz_yy[i] * fxi[i] + tg_yyzzz_yy[i] * ra_z[i];

        tg_yyzzzz_yz[i] = tg_zzzz_yz[i] * fxi[i] + tg_yzzzz_z[i] * fxi[i] + tg_yzzzz_yz[i] * ra_y[i];

        tg_yyzzzz_zz[i] = tg_zzzz_zz[i] * fxi[i] + tg_yzzzz_zz[i] * ra_y[i];

        tg_yzzzzz_xx[i] = tg_zzzzz_xx[i] * ra_y[i];

        tg_yzzzzz_xy[i] = tg_zzzzz_x[i] * fxi[i] + tg_zzzzz_xy[i] * ra_y[i];

        tg_yzzzzz_xz[i] = tg_zzzzz_xz[i] * ra_y[i];

        tg_yzzzzz_yy[i] = 2.0 * tg_zzzzz_y[i] * fxi[i] + tg_zzzzz_yy[i] * ra_y[i];

        tg_yzzzzz_yz[i] = tg_zzzzz_z[i] * fxi[i] + tg_zzzzz_yz[i] * ra_y[i];

        tg_yzzzzz_zz[i] = tg_zzzzz_zz[i] * ra_y[i];

        tg_zzzzzz_xx[i] = 5.0 * tg_zzzz_xx[i] * fxi[i] + tg_zzzzz_xx[i] * ra_z[i];

        tg_zzzzzz_xy[i] = 5.0 * tg_zzzz_xy[i] * fxi[i] + tg_zzzzz_xy[i] * ra_z[i];

        tg_zzzzzz_xz[i] = 5.0 * tg_zzzz_xz[i] * fxi[i] + tg_zzzzz_x[i] * fxi[i] + tg_zzzzz_xz[i] * ra_z[i];

        tg_zzzzzz_yy[i] = 5.0 * tg_zzzz_yy[i] * fxi[i] + tg_zzzzz_yy[i] * ra_z[i];

        tg_zzzzzz_yz[i] = 5.0 * tg_zzzz_yz[i] * fxi[i] + tg_zzzzz_y[i] * fxi[i] + tg_zzzzz_yz[i] * ra_z[i];

        tg_zzzzzz_zz[i] = 5.0 * tg_zzzz_zz[i] * fxi[i] + 2.0 * tg_zzzzz_z[i] * fxi[i] + tg_zzzzz_zz[i] * ra_z[i];
    }
}

} // t2lecp namespace

