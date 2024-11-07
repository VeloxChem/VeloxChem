#include "OverlapPrimRecFI.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_fi(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_fi,
                     const size_t              idx_ovl_pi,
                     const size_t              idx_ovl_dh,
                     const size_t              idx_ovl_di,
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

    // Set up components of auxiliary buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_ovl_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_ovl_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_ovl_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_ovl_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_ovl_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_ovl_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_ovl_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_ovl_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_ovl_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_ovl_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_ovl_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_ovl_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_ovl_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_ovl_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_ovl_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_ovl_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_ovl_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_ovl_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_ovl_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_ovl_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_ovl_dh + 20);

    auto ts_yy_xxxxx = pbuffer.data(idx_ovl_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_ovl_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_ovl_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_ovl_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_ovl_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_ovl_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_ovl_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_ovl_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_ovl_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_ovl_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_ovl_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_ovl_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_ovl_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_ovl_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_ovl_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_ovl_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_ovl_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_ovl_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_ovl_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_ovl_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_ovl_dh + 83);

    auto ts_yz_xxxyz = pbuffer.data(idx_ovl_dh + 88);

    auto ts_yz_xxyyz = pbuffer.data(idx_ovl_dh + 91);

    auto ts_yz_xxyzz = pbuffer.data(idx_ovl_dh + 92);

    auto ts_yz_xyyyz = pbuffer.data(idx_ovl_dh + 95);

    auto ts_yz_xyyzz = pbuffer.data(idx_ovl_dh + 96);

    auto ts_yz_xyzzz = pbuffer.data(idx_ovl_dh + 97);

    auto ts_yz_yyyyz = pbuffer.data(idx_ovl_dh + 100);

    auto ts_yz_yyyzz = pbuffer.data(idx_ovl_dh + 101);

    auto ts_yz_yyzzz = pbuffer.data(idx_ovl_dh + 102);

    auto ts_yz_yzzzz = pbuffer.data(idx_ovl_dh + 103);

    auto ts_zz_xxxxx = pbuffer.data(idx_ovl_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_ovl_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_ovl_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_ovl_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_ovl_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_ovl_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_ovl_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_ovl_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_ovl_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_ovl_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_ovl_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_ovl_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_ovl_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_ovl_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_ovl_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_ovl_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_ovl_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_ovl_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_ovl_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_ovl_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_ovl_dh + 125);

    // Set up components of auxiliary buffer : DI

    auto ts_xx_xxxxxx = pbuffer.data(idx_ovl_di);

    auto ts_xx_xxxxxy = pbuffer.data(idx_ovl_di + 1);

    auto ts_xx_xxxxxz = pbuffer.data(idx_ovl_di + 2);

    auto ts_xx_xxxxyy = pbuffer.data(idx_ovl_di + 3);

    auto ts_xx_xxxxyz = pbuffer.data(idx_ovl_di + 4);

    auto ts_xx_xxxxzz = pbuffer.data(idx_ovl_di + 5);

    auto ts_xx_xxxyyy = pbuffer.data(idx_ovl_di + 6);

    auto ts_xx_xxxyyz = pbuffer.data(idx_ovl_di + 7);

    auto ts_xx_xxxyzz = pbuffer.data(idx_ovl_di + 8);

    auto ts_xx_xxxzzz = pbuffer.data(idx_ovl_di + 9);

    auto ts_xx_xxyyyy = pbuffer.data(idx_ovl_di + 10);

    auto ts_xx_xxyyyz = pbuffer.data(idx_ovl_di + 11);

    auto ts_xx_xxyyzz = pbuffer.data(idx_ovl_di + 12);

    auto ts_xx_xxyzzz = pbuffer.data(idx_ovl_di + 13);

    auto ts_xx_xxzzzz = pbuffer.data(idx_ovl_di + 14);

    auto ts_xx_xyyyyy = pbuffer.data(idx_ovl_di + 15);

    auto ts_xx_xyyyyz = pbuffer.data(idx_ovl_di + 16);

    auto ts_xx_xyyyzz = pbuffer.data(idx_ovl_di + 17);

    auto ts_xx_xyyzzz = pbuffer.data(idx_ovl_di + 18);

    auto ts_xx_xyzzzz = pbuffer.data(idx_ovl_di + 19);

    auto ts_xx_xzzzzz = pbuffer.data(idx_ovl_di + 20);

    auto ts_xx_yyyyyy = pbuffer.data(idx_ovl_di + 21);

    auto ts_xx_yyyyyz = pbuffer.data(idx_ovl_di + 22);

    auto ts_xx_yyyyzz = pbuffer.data(idx_ovl_di + 23);

    auto ts_xx_yyyzzz = pbuffer.data(idx_ovl_di + 24);

    auto ts_xx_yyzzzz = pbuffer.data(idx_ovl_di + 25);

    auto ts_xx_yzzzzz = pbuffer.data(idx_ovl_di + 26);

    auto ts_xx_zzzzzz = pbuffer.data(idx_ovl_di + 27);

    auto ts_xy_xxxxxy = pbuffer.data(idx_ovl_di + 29);

    auto ts_xy_xxxxyy = pbuffer.data(idx_ovl_di + 31);

    auto ts_xy_xxxyyy = pbuffer.data(idx_ovl_di + 34);

    auto ts_xy_xxyyyy = pbuffer.data(idx_ovl_di + 38);

    auto ts_xy_xyyyyy = pbuffer.data(idx_ovl_di + 43);

    auto ts_xy_yyyyyy = pbuffer.data(idx_ovl_di + 49);

    auto ts_xy_yyyyyz = pbuffer.data(idx_ovl_di + 50);

    auto ts_xy_yyyyzz = pbuffer.data(idx_ovl_di + 51);

    auto ts_xy_yyyzzz = pbuffer.data(idx_ovl_di + 52);

    auto ts_xy_yyzzzz = pbuffer.data(idx_ovl_di + 53);

    auto ts_xy_yzzzzz = pbuffer.data(idx_ovl_di + 54);

    auto ts_xz_xxxxxx = pbuffer.data(idx_ovl_di + 56);

    auto ts_xz_xxxxxz = pbuffer.data(idx_ovl_di + 58);

    auto ts_xz_xxxxzz = pbuffer.data(idx_ovl_di + 61);

    auto ts_xz_xxxzzz = pbuffer.data(idx_ovl_di + 65);

    auto ts_xz_xxzzzz = pbuffer.data(idx_ovl_di + 70);

    auto ts_xz_xzzzzz = pbuffer.data(idx_ovl_di + 76);

    auto ts_xz_yyyyyz = pbuffer.data(idx_ovl_di + 78);

    auto ts_xz_yyyyzz = pbuffer.data(idx_ovl_di + 79);

    auto ts_xz_yyyzzz = pbuffer.data(idx_ovl_di + 80);

    auto ts_xz_yyzzzz = pbuffer.data(idx_ovl_di + 81);

    auto ts_xz_yzzzzz = pbuffer.data(idx_ovl_di + 82);

    auto ts_xz_zzzzzz = pbuffer.data(idx_ovl_di + 83);

    auto ts_yy_xxxxxx = pbuffer.data(idx_ovl_di + 84);

    auto ts_yy_xxxxxy = pbuffer.data(idx_ovl_di + 85);

    auto ts_yy_xxxxxz = pbuffer.data(idx_ovl_di + 86);

    auto ts_yy_xxxxyy = pbuffer.data(idx_ovl_di + 87);

    auto ts_yy_xxxxyz = pbuffer.data(idx_ovl_di + 88);

    auto ts_yy_xxxxzz = pbuffer.data(idx_ovl_di + 89);

    auto ts_yy_xxxyyy = pbuffer.data(idx_ovl_di + 90);

    auto ts_yy_xxxyyz = pbuffer.data(idx_ovl_di + 91);

    auto ts_yy_xxxyzz = pbuffer.data(idx_ovl_di + 92);

    auto ts_yy_xxxzzz = pbuffer.data(idx_ovl_di + 93);

    auto ts_yy_xxyyyy = pbuffer.data(idx_ovl_di + 94);

    auto ts_yy_xxyyyz = pbuffer.data(idx_ovl_di + 95);

    auto ts_yy_xxyyzz = pbuffer.data(idx_ovl_di + 96);

    auto ts_yy_xxyzzz = pbuffer.data(idx_ovl_di + 97);

    auto ts_yy_xxzzzz = pbuffer.data(idx_ovl_di + 98);

    auto ts_yy_xyyyyy = pbuffer.data(idx_ovl_di + 99);

    auto ts_yy_xyyyyz = pbuffer.data(idx_ovl_di + 100);

    auto ts_yy_xyyyzz = pbuffer.data(idx_ovl_di + 101);

    auto ts_yy_xyyzzz = pbuffer.data(idx_ovl_di + 102);

    auto ts_yy_xyzzzz = pbuffer.data(idx_ovl_di + 103);

    auto ts_yy_xzzzzz = pbuffer.data(idx_ovl_di + 104);

    auto ts_yy_yyyyyy = pbuffer.data(idx_ovl_di + 105);

    auto ts_yy_yyyyyz = pbuffer.data(idx_ovl_di + 106);

    auto ts_yy_yyyyzz = pbuffer.data(idx_ovl_di + 107);

    auto ts_yy_yyyzzz = pbuffer.data(idx_ovl_di + 108);

    auto ts_yy_yyzzzz = pbuffer.data(idx_ovl_di + 109);

    auto ts_yy_yzzzzz = pbuffer.data(idx_ovl_di + 110);

    auto ts_yy_zzzzzz = pbuffer.data(idx_ovl_di + 111);

    auto ts_yz_xxxxxz = pbuffer.data(idx_ovl_di + 114);

    auto ts_yz_xxxxyz = pbuffer.data(idx_ovl_di + 116);

    auto ts_yz_xxxxzz = pbuffer.data(idx_ovl_di + 117);

    auto ts_yz_xxxyyz = pbuffer.data(idx_ovl_di + 119);

    auto ts_yz_xxxyzz = pbuffer.data(idx_ovl_di + 120);

    auto ts_yz_xxxzzz = pbuffer.data(idx_ovl_di + 121);

    auto ts_yz_xxyyyz = pbuffer.data(idx_ovl_di + 123);

    auto ts_yz_xxyyzz = pbuffer.data(idx_ovl_di + 124);

    auto ts_yz_xxyzzz = pbuffer.data(idx_ovl_di + 125);

    auto ts_yz_xxzzzz = pbuffer.data(idx_ovl_di + 126);

    auto ts_yz_xyyyyz = pbuffer.data(idx_ovl_di + 128);

    auto ts_yz_xyyyzz = pbuffer.data(idx_ovl_di + 129);

    auto ts_yz_xyyzzz = pbuffer.data(idx_ovl_di + 130);

    auto ts_yz_xyzzzz = pbuffer.data(idx_ovl_di + 131);

    auto ts_yz_xzzzzz = pbuffer.data(idx_ovl_di + 132);

    auto ts_yz_yyyyyy = pbuffer.data(idx_ovl_di + 133);

    auto ts_yz_yyyyyz = pbuffer.data(idx_ovl_di + 134);

    auto ts_yz_yyyyzz = pbuffer.data(idx_ovl_di + 135);

    auto ts_yz_yyyzzz = pbuffer.data(idx_ovl_di + 136);

    auto ts_yz_yyzzzz = pbuffer.data(idx_ovl_di + 137);

    auto ts_yz_yzzzzz = pbuffer.data(idx_ovl_di + 138);

    auto ts_yz_zzzzzz = pbuffer.data(idx_ovl_di + 139);

    auto ts_zz_xxxxxx = pbuffer.data(idx_ovl_di + 140);

    auto ts_zz_xxxxxy = pbuffer.data(idx_ovl_di + 141);

    auto ts_zz_xxxxxz = pbuffer.data(idx_ovl_di + 142);

    auto ts_zz_xxxxyy = pbuffer.data(idx_ovl_di + 143);

    auto ts_zz_xxxxyz = pbuffer.data(idx_ovl_di + 144);

    auto ts_zz_xxxxzz = pbuffer.data(idx_ovl_di + 145);

    auto ts_zz_xxxyyy = pbuffer.data(idx_ovl_di + 146);

    auto ts_zz_xxxyyz = pbuffer.data(idx_ovl_di + 147);

    auto ts_zz_xxxyzz = pbuffer.data(idx_ovl_di + 148);

    auto ts_zz_xxxzzz = pbuffer.data(idx_ovl_di + 149);

    auto ts_zz_xxyyyy = pbuffer.data(idx_ovl_di + 150);

    auto ts_zz_xxyyyz = pbuffer.data(idx_ovl_di + 151);

    auto ts_zz_xxyyzz = pbuffer.data(idx_ovl_di + 152);

    auto ts_zz_xxyzzz = pbuffer.data(idx_ovl_di + 153);

    auto ts_zz_xxzzzz = pbuffer.data(idx_ovl_di + 154);

    auto ts_zz_xyyyyy = pbuffer.data(idx_ovl_di + 155);

    auto ts_zz_xyyyyz = pbuffer.data(idx_ovl_di + 156);

    auto ts_zz_xyyyzz = pbuffer.data(idx_ovl_di + 157);

    auto ts_zz_xyyzzz = pbuffer.data(idx_ovl_di + 158);

    auto ts_zz_xyzzzz = pbuffer.data(idx_ovl_di + 159);

    auto ts_zz_xzzzzz = pbuffer.data(idx_ovl_di + 160);

    auto ts_zz_yyyyyy = pbuffer.data(idx_ovl_di + 161);

    auto ts_zz_yyyyyz = pbuffer.data(idx_ovl_di + 162);

    auto ts_zz_yyyyzz = pbuffer.data(idx_ovl_di + 163);

    auto ts_zz_yyyzzz = pbuffer.data(idx_ovl_di + 164);

    auto ts_zz_yyzzzz = pbuffer.data(idx_ovl_di + 165);

    auto ts_zz_yzzzzz = pbuffer.data(idx_ovl_di + 166);

    auto ts_zz_zzzzzz = pbuffer.data(idx_ovl_di + 167);

    // Set up 0-28 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_x,              \
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
                             ts_xx_xxxxx,   \
                             ts_xx_xxxxxx,  \
                             ts_xx_xxxxxy,  \
                             ts_xx_xxxxxz,  \
                             ts_xx_xxxxy,   \
                             ts_xx_xxxxyy,  \
                             ts_xx_xxxxyz,  \
                             ts_xx_xxxxz,   \
                             ts_xx_xxxxzz,  \
                             ts_xx_xxxyy,   \
                             ts_xx_xxxyyy,  \
                             ts_xx_xxxyyz,  \
                             ts_xx_xxxyz,   \
                             ts_xx_xxxyzz,  \
                             ts_xx_xxxzz,   \
                             ts_xx_xxxzzz,  \
                             ts_xx_xxyyy,   \
                             ts_xx_xxyyyy,  \
                             ts_xx_xxyyyz,  \
                             ts_xx_xxyyz,   \
                             ts_xx_xxyyzz,  \
                             ts_xx_xxyzz,   \
                             ts_xx_xxyzzz,  \
                             ts_xx_xxzzz,   \
                             ts_xx_xxzzzz,  \
                             ts_xx_xyyyy,   \
                             ts_xx_xyyyyy,  \
                             ts_xx_xyyyyz,  \
                             ts_xx_xyyyz,   \
                             ts_xx_xyyyzz,  \
                             ts_xx_xyyzz,   \
                             ts_xx_xyyzzz,  \
                             ts_xx_xyzzz,   \
                             ts_xx_xyzzzz,  \
                             ts_xx_xzzzz,   \
                             ts_xx_xzzzzz,  \
                             ts_xx_yyyyy,   \
                             ts_xx_yyyyyy,  \
                             ts_xx_yyyyyz,  \
                             ts_xx_yyyyz,   \
                             ts_xx_yyyyzz,  \
                             ts_xx_yyyzz,   \
                             ts_xx_yyyzzz,  \
                             ts_xx_yyzzz,   \
                             ts_xx_yyzzzz,  \
                             ts_xx_yzzzz,   \
                             ts_xx_yzzzzz,  \
                             ts_xx_zzzzz,   \
                             ts_xx_zzzzzz,  \
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

        ts_xxx_xxxxxx[i] = 2.0 * ts_x_xxxxxx[i] * fe_0 + 6.0 * ts_xx_xxxxx[i] * fe_0 + ts_xx_xxxxxx[i] * pa_x[i];

        ts_xxx_xxxxxy[i] = 2.0 * ts_x_xxxxxy[i] * fe_0 + 5.0 * ts_xx_xxxxy[i] * fe_0 + ts_xx_xxxxxy[i] * pa_x[i];

        ts_xxx_xxxxxz[i] = 2.0 * ts_x_xxxxxz[i] * fe_0 + 5.0 * ts_xx_xxxxz[i] * fe_0 + ts_xx_xxxxxz[i] * pa_x[i];

        ts_xxx_xxxxyy[i] = 2.0 * ts_x_xxxxyy[i] * fe_0 + 4.0 * ts_xx_xxxyy[i] * fe_0 + ts_xx_xxxxyy[i] * pa_x[i];

        ts_xxx_xxxxyz[i] = 2.0 * ts_x_xxxxyz[i] * fe_0 + 4.0 * ts_xx_xxxyz[i] * fe_0 + ts_xx_xxxxyz[i] * pa_x[i];

        ts_xxx_xxxxzz[i] = 2.0 * ts_x_xxxxzz[i] * fe_0 + 4.0 * ts_xx_xxxzz[i] * fe_0 + ts_xx_xxxxzz[i] * pa_x[i];

        ts_xxx_xxxyyy[i] = 2.0 * ts_x_xxxyyy[i] * fe_0 + 3.0 * ts_xx_xxyyy[i] * fe_0 + ts_xx_xxxyyy[i] * pa_x[i];

        ts_xxx_xxxyyz[i] = 2.0 * ts_x_xxxyyz[i] * fe_0 + 3.0 * ts_xx_xxyyz[i] * fe_0 + ts_xx_xxxyyz[i] * pa_x[i];

        ts_xxx_xxxyzz[i] = 2.0 * ts_x_xxxyzz[i] * fe_0 + 3.0 * ts_xx_xxyzz[i] * fe_0 + ts_xx_xxxyzz[i] * pa_x[i];

        ts_xxx_xxxzzz[i] = 2.0 * ts_x_xxxzzz[i] * fe_0 + 3.0 * ts_xx_xxzzz[i] * fe_0 + ts_xx_xxxzzz[i] * pa_x[i];

        ts_xxx_xxyyyy[i] = 2.0 * ts_x_xxyyyy[i] * fe_0 + 2.0 * ts_xx_xyyyy[i] * fe_0 + ts_xx_xxyyyy[i] * pa_x[i];

        ts_xxx_xxyyyz[i] = 2.0 * ts_x_xxyyyz[i] * fe_0 + 2.0 * ts_xx_xyyyz[i] * fe_0 + ts_xx_xxyyyz[i] * pa_x[i];

        ts_xxx_xxyyzz[i] = 2.0 * ts_x_xxyyzz[i] * fe_0 + 2.0 * ts_xx_xyyzz[i] * fe_0 + ts_xx_xxyyzz[i] * pa_x[i];

        ts_xxx_xxyzzz[i] = 2.0 * ts_x_xxyzzz[i] * fe_0 + 2.0 * ts_xx_xyzzz[i] * fe_0 + ts_xx_xxyzzz[i] * pa_x[i];

        ts_xxx_xxzzzz[i] = 2.0 * ts_x_xxzzzz[i] * fe_0 + 2.0 * ts_xx_xzzzz[i] * fe_0 + ts_xx_xxzzzz[i] * pa_x[i];

        ts_xxx_xyyyyy[i] = 2.0 * ts_x_xyyyyy[i] * fe_0 + ts_xx_yyyyy[i] * fe_0 + ts_xx_xyyyyy[i] * pa_x[i];

        ts_xxx_xyyyyz[i] = 2.0 * ts_x_xyyyyz[i] * fe_0 + ts_xx_yyyyz[i] * fe_0 + ts_xx_xyyyyz[i] * pa_x[i];

        ts_xxx_xyyyzz[i] = 2.0 * ts_x_xyyyzz[i] * fe_0 + ts_xx_yyyzz[i] * fe_0 + ts_xx_xyyyzz[i] * pa_x[i];

        ts_xxx_xyyzzz[i] = 2.0 * ts_x_xyyzzz[i] * fe_0 + ts_xx_yyzzz[i] * fe_0 + ts_xx_xyyzzz[i] * pa_x[i];

        ts_xxx_xyzzzz[i] = 2.0 * ts_x_xyzzzz[i] * fe_0 + ts_xx_yzzzz[i] * fe_0 + ts_xx_xyzzzz[i] * pa_x[i];

        ts_xxx_xzzzzz[i] = 2.0 * ts_x_xzzzzz[i] * fe_0 + ts_xx_zzzzz[i] * fe_0 + ts_xx_xzzzzz[i] * pa_x[i];

        ts_xxx_yyyyyy[i] = 2.0 * ts_x_yyyyyy[i] * fe_0 + ts_xx_yyyyyy[i] * pa_x[i];

        ts_xxx_yyyyyz[i] = 2.0 * ts_x_yyyyyz[i] * fe_0 + ts_xx_yyyyyz[i] * pa_x[i];

        ts_xxx_yyyyzz[i] = 2.0 * ts_x_yyyyzz[i] * fe_0 + ts_xx_yyyyzz[i] * pa_x[i];

        ts_xxx_yyyzzz[i] = 2.0 * ts_x_yyyzzz[i] * fe_0 + ts_xx_yyyzzz[i] * pa_x[i];

        ts_xxx_yyzzzz[i] = 2.0 * ts_x_yyzzzz[i] * fe_0 + ts_xx_yyzzzz[i] * pa_x[i];

        ts_xxx_yzzzzz[i] = 2.0 * ts_x_yzzzzz[i] * fe_0 + ts_xx_yzzzzz[i] * pa_x[i];

        ts_xxx_zzzzzz[i] = 2.0 * ts_x_zzzzzz[i] * fe_0 + ts_xx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             ts_xx_xxxxx,   \
                             ts_xx_xxxxxx,  \
                             ts_xx_xxxxxy,  \
                             ts_xx_xxxxxz,  \
                             ts_xx_xxxxy,   \
                             ts_xx_xxxxyy,  \
                             ts_xx_xxxxyz,  \
                             ts_xx_xxxxz,   \
                             ts_xx_xxxxzz,  \
                             ts_xx_xxxyy,   \
                             ts_xx_xxxyyy,  \
                             ts_xx_xxxyyz,  \
                             ts_xx_xxxyz,   \
                             ts_xx_xxxyzz,  \
                             ts_xx_xxxzz,   \
                             ts_xx_xxxzzz,  \
                             ts_xx_xxyyy,   \
                             ts_xx_xxyyyy,  \
                             ts_xx_xxyyyz,  \
                             ts_xx_xxyyz,   \
                             ts_xx_xxyyzz,  \
                             ts_xx_xxyzz,   \
                             ts_xx_xxyzzz,  \
                             ts_xx_xxzzz,   \
                             ts_xx_xxzzzz,  \
                             ts_xx_xyyyy,   \
                             ts_xx_xyyyyy,  \
                             ts_xx_xyyyyz,  \
                             ts_xx_xyyyz,   \
                             ts_xx_xyyyzz,  \
                             ts_xx_xyyzz,   \
                             ts_xx_xyyzzz,  \
                             ts_xx_xyzzz,   \
                             ts_xx_xyzzzz,  \
                             ts_xx_xzzzz,   \
                             ts_xx_xzzzzz,  \
                             ts_xx_zzzzzz,  \
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
                             ts_xy_yyyyyy,  \
                             ts_xy_yyyyyz,  \
                             ts_xy_yyyyzz,  \
                             ts_xy_yyyzzz,  \
                             ts_xy_yyzzzz,  \
                             ts_xy_yzzzzz,  \
                             ts_y_yyyyyy,   \
                             ts_y_yyyyyz,   \
                             ts_y_yyyyzz,   \
                             ts_y_yyyzzz,   \
                             ts_y_yyzzzz,   \
                             ts_y_yzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxy_xxxxxx[i] = ts_xx_xxxxxx[i] * pa_y[i];

        ts_xxy_xxxxxy[i] = ts_xx_xxxxx[i] * fe_0 + ts_xx_xxxxxy[i] * pa_y[i];

        ts_xxy_xxxxxz[i] = ts_xx_xxxxxz[i] * pa_y[i];

        ts_xxy_xxxxyy[i] = 2.0 * ts_xx_xxxxy[i] * fe_0 + ts_xx_xxxxyy[i] * pa_y[i];

        ts_xxy_xxxxyz[i] = ts_xx_xxxxz[i] * fe_0 + ts_xx_xxxxyz[i] * pa_y[i];

        ts_xxy_xxxxzz[i] = ts_xx_xxxxzz[i] * pa_y[i];

        ts_xxy_xxxyyy[i] = 3.0 * ts_xx_xxxyy[i] * fe_0 + ts_xx_xxxyyy[i] * pa_y[i];

        ts_xxy_xxxyyz[i] = 2.0 * ts_xx_xxxyz[i] * fe_0 + ts_xx_xxxyyz[i] * pa_y[i];

        ts_xxy_xxxyzz[i] = ts_xx_xxxzz[i] * fe_0 + ts_xx_xxxyzz[i] * pa_y[i];

        ts_xxy_xxxzzz[i] = ts_xx_xxxzzz[i] * pa_y[i];

        ts_xxy_xxyyyy[i] = 4.0 * ts_xx_xxyyy[i] * fe_0 + ts_xx_xxyyyy[i] * pa_y[i];

        ts_xxy_xxyyyz[i] = 3.0 * ts_xx_xxyyz[i] * fe_0 + ts_xx_xxyyyz[i] * pa_y[i];

        ts_xxy_xxyyzz[i] = 2.0 * ts_xx_xxyzz[i] * fe_0 + ts_xx_xxyyzz[i] * pa_y[i];

        ts_xxy_xxyzzz[i] = ts_xx_xxzzz[i] * fe_0 + ts_xx_xxyzzz[i] * pa_y[i];

        ts_xxy_xxzzzz[i] = ts_xx_xxzzzz[i] * pa_y[i];

        ts_xxy_xyyyyy[i] = 5.0 * ts_xx_xyyyy[i] * fe_0 + ts_xx_xyyyyy[i] * pa_y[i];

        ts_xxy_xyyyyz[i] = 4.0 * ts_xx_xyyyz[i] * fe_0 + ts_xx_xyyyyz[i] * pa_y[i];

        ts_xxy_xyyyzz[i] = 3.0 * ts_xx_xyyzz[i] * fe_0 + ts_xx_xyyyzz[i] * pa_y[i];

        ts_xxy_xyyzzz[i] = 2.0 * ts_xx_xyzzz[i] * fe_0 + ts_xx_xyyzzz[i] * pa_y[i];

        ts_xxy_xyzzzz[i] = ts_xx_xzzzz[i] * fe_0 + ts_xx_xyzzzz[i] * pa_y[i];

        ts_xxy_xzzzzz[i] = ts_xx_xzzzzz[i] * pa_y[i];

        ts_xxy_yyyyyy[i] = ts_y_yyyyyy[i] * fe_0 + ts_xy_yyyyyy[i] * pa_x[i];

        ts_xxy_yyyyyz[i] = ts_y_yyyyyz[i] * fe_0 + ts_xy_yyyyyz[i] * pa_x[i];

        ts_xxy_yyyyzz[i] = ts_y_yyyyzz[i] * fe_0 + ts_xy_yyyyzz[i] * pa_x[i];

        ts_xxy_yyyzzz[i] = ts_y_yyyzzz[i] * fe_0 + ts_xy_yyyzzz[i] * pa_x[i];

        ts_xxy_yyzzzz[i] = ts_y_yyzzzz[i] * fe_0 + ts_xy_yyzzzz[i] * pa_x[i];

        ts_xxy_yzzzzz[i] = ts_y_yzzzzz[i] * fe_0 + ts_xy_yzzzzz[i] * pa_x[i];

        ts_xxy_zzzzzz[i] = ts_xx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             ts_xx_xxxxx,   \
                             ts_xx_xxxxxx,  \
                             ts_xx_xxxxxy,  \
                             ts_xx_xxxxxz,  \
                             ts_xx_xxxxy,   \
                             ts_xx_xxxxyy,  \
                             ts_xx_xxxxyz,  \
                             ts_xx_xxxxz,   \
                             ts_xx_xxxxzz,  \
                             ts_xx_xxxyy,   \
                             ts_xx_xxxyyy,  \
                             ts_xx_xxxyyz,  \
                             ts_xx_xxxyz,   \
                             ts_xx_xxxyzz,  \
                             ts_xx_xxxzz,   \
                             ts_xx_xxxzzz,  \
                             ts_xx_xxyyy,   \
                             ts_xx_xxyyyy,  \
                             ts_xx_xxyyyz,  \
                             ts_xx_xxyyz,   \
                             ts_xx_xxyyzz,  \
                             ts_xx_xxyzz,   \
                             ts_xx_xxyzzz,  \
                             ts_xx_xxzzz,   \
                             ts_xx_xxzzzz,  \
                             ts_xx_xyyyy,   \
                             ts_xx_xyyyyy,  \
                             ts_xx_xyyyyz,  \
                             ts_xx_xyyyz,   \
                             ts_xx_xyyyzz,  \
                             ts_xx_xyyzz,   \
                             ts_xx_xyyzzz,  \
                             ts_xx_xyzzz,   \
                             ts_xx_xyzzzz,  \
                             ts_xx_xzzzz,   \
                             ts_xx_xzzzzz,  \
                             ts_xx_yyyyyy,  \
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
                             ts_xz_yyyyyz,  \
                             ts_xz_yyyyzz,  \
                             ts_xz_yyyzzz,  \
                             ts_xz_yyzzzz,  \
                             ts_xz_yzzzzz,  \
                             ts_xz_zzzzzz,  \
                             ts_z_yyyyyz,   \
                             ts_z_yyyyzz,   \
                             ts_z_yyyzzz,   \
                             ts_z_yyzzzz,   \
                             ts_z_yzzzzz,   \
                             ts_z_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxz_xxxxxx[i] = ts_xx_xxxxxx[i] * pa_z[i];

        ts_xxz_xxxxxy[i] = ts_xx_xxxxxy[i] * pa_z[i];

        ts_xxz_xxxxxz[i] = ts_xx_xxxxx[i] * fe_0 + ts_xx_xxxxxz[i] * pa_z[i];

        ts_xxz_xxxxyy[i] = ts_xx_xxxxyy[i] * pa_z[i];

        ts_xxz_xxxxyz[i] = ts_xx_xxxxy[i] * fe_0 + ts_xx_xxxxyz[i] * pa_z[i];

        ts_xxz_xxxxzz[i] = 2.0 * ts_xx_xxxxz[i] * fe_0 + ts_xx_xxxxzz[i] * pa_z[i];

        ts_xxz_xxxyyy[i] = ts_xx_xxxyyy[i] * pa_z[i];

        ts_xxz_xxxyyz[i] = ts_xx_xxxyy[i] * fe_0 + ts_xx_xxxyyz[i] * pa_z[i];

        ts_xxz_xxxyzz[i] = 2.0 * ts_xx_xxxyz[i] * fe_0 + ts_xx_xxxyzz[i] * pa_z[i];

        ts_xxz_xxxzzz[i] = 3.0 * ts_xx_xxxzz[i] * fe_0 + ts_xx_xxxzzz[i] * pa_z[i];

        ts_xxz_xxyyyy[i] = ts_xx_xxyyyy[i] * pa_z[i];

        ts_xxz_xxyyyz[i] = ts_xx_xxyyy[i] * fe_0 + ts_xx_xxyyyz[i] * pa_z[i];

        ts_xxz_xxyyzz[i] = 2.0 * ts_xx_xxyyz[i] * fe_0 + ts_xx_xxyyzz[i] * pa_z[i];

        ts_xxz_xxyzzz[i] = 3.0 * ts_xx_xxyzz[i] * fe_0 + ts_xx_xxyzzz[i] * pa_z[i];

        ts_xxz_xxzzzz[i] = 4.0 * ts_xx_xxzzz[i] * fe_0 + ts_xx_xxzzzz[i] * pa_z[i];

        ts_xxz_xyyyyy[i] = ts_xx_xyyyyy[i] * pa_z[i];

        ts_xxz_xyyyyz[i] = ts_xx_xyyyy[i] * fe_0 + ts_xx_xyyyyz[i] * pa_z[i];

        ts_xxz_xyyyzz[i] = 2.0 * ts_xx_xyyyz[i] * fe_0 + ts_xx_xyyyzz[i] * pa_z[i];

        ts_xxz_xyyzzz[i] = 3.0 * ts_xx_xyyzz[i] * fe_0 + ts_xx_xyyzzz[i] * pa_z[i];

        ts_xxz_xyzzzz[i] = 4.0 * ts_xx_xyzzz[i] * fe_0 + ts_xx_xyzzzz[i] * pa_z[i];

        ts_xxz_xzzzzz[i] = 5.0 * ts_xx_xzzzz[i] * fe_0 + ts_xx_xzzzzz[i] * pa_z[i];

        ts_xxz_yyyyyy[i] = ts_xx_yyyyyy[i] * pa_z[i];

        ts_xxz_yyyyyz[i] = ts_z_yyyyyz[i] * fe_0 + ts_xz_yyyyyz[i] * pa_x[i];

        ts_xxz_yyyyzz[i] = ts_z_yyyyzz[i] * fe_0 + ts_xz_yyyyzz[i] * pa_x[i];

        ts_xxz_yyyzzz[i] = ts_z_yyyzzz[i] * fe_0 + ts_xz_yyyzzz[i] * pa_x[i];

        ts_xxz_yyzzzz[i] = ts_z_yyzzzz[i] * fe_0 + ts_xz_yyzzzz[i] * pa_x[i];

        ts_xxz_yzzzzz[i] = ts_z_yzzzzz[i] * fe_0 + ts_xz_yzzzzz[i] * pa_x[i];

        ts_xxz_zzzzzz[i] = ts_z_zzzzzz[i] * fe_0 + ts_xz_zzzzzz[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_x,              \
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
                             ts_yy_xxxxx,   \
                             ts_yy_xxxxxx,  \
                             ts_yy_xxxxxy,  \
                             ts_yy_xxxxxz,  \
                             ts_yy_xxxxy,   \
                             ts_yy_xxxxyy,  \
                             ts_yy_xxxxyz,  \
                             ts_yy_xxxxz,   \
                             ts_yy_xxxxzz,  \
                             ts_yy_xxxyy,   \
                             ts_yy_xxxyyy,  \
                             ts_yy_xxxyyz,  \
                             ts_yy_xxxyz,   \
                             ts_yy_xxxyzz,  \
                             ts_yy_xxxzz,   \
                             ts_yy_xxxzzz,  \
                             ts_yy_xxyyy,   \
                             ts_yy_xxyyyy,  \
                             ts_yy_xxyyyz,  \
                             ts_yy_xxyyz,   \
                             ts_yy_xxyyzz,  \
                             ts_yy_xxyzz,   \
                             ts_yy_xxyzzz,  \
                             ts_yy_xxzzz,   \
                             ts_yy_xxzzzz,  \
                             ts_yy_xyyyy,   \
                             ts_yy_xyyyyy,  \
                             ts_yy_xyyyyz,  \
                             ts_yy_xyyyz,   \
                             ts_yy_xyyyzz,  \
                             ts_yy_xyyzz,   \
                             ts_yy_xyyzzz,  \
                             ts_yy_xyzzz,   \
                             ts_yy_xyzzzz,  \
                             ts_yy_xzzzz,   \
                             ts_yy_xzzzzz,  \
                             ts_yy_yyyyy,   \
                             ts_yy_yyyyyy,  \
                             ts_yy_yyyyyz,  \
                             ts_yy_yyyyz,   \
                             ts_yy_yyyyzz,  \
                             ts_yy_yyyzz,   \
                             ts_yy_yyyzzz,  \
                             ts_yy_yyzzz,   \
                             ts_yy_yyzzzz,  \
                             ts_yy_yzzzz,   \
                             ts_yy_yzzzzz,  \
                             ts_yy_zzzzz,   \
                             ts_yy_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyy_xxxxxx[i] = 6.0 * ts_yy_xxxxx[i] * fe_0 + ts_yy_xxxxxx[i] * pa_x[i];

        ts_xyy_xxxxxy[i] = 5.0 * ts_yy_xxxxy[i] * fe_0 + ts_yy_xxxxxy[i] * pa_x[i];

        ts_xyy_xxxxxz[i] = 5.0 * ts_yy_xxxxz[i] * fe_0 + ts_yy_xxxxxz[i] * pa_x[i];

        ts_xyy_xxxxyy[i] = 4.0 * ts_yy_xxxyy[i] * fe_0 + ts_yy_xxxxyy[i] * pa_x[i];

        ts_xyy_xxxxyz[i] = 4.0 * ts_yy_xxxyz[i] * fe_0 + ts_yy_xxxxyz[i] * pa_x[i];

        ts_xyy_xxxxzz[i] = 4.0 * ts_yy_xxxzz[i] * fe_0 + ts_yy_xxxxzz[i] * pa_x[i];

        ts_xyy_xxxyyy[i] = 3.0 * ts_yy_xxyyy[i] * fe_0 + ts_yy_xxxyyy[i] * pa_x[i];

        ts_xyy_xxxyyz[i] = 3.0 * ts_yy_xxyyz[i] * fe_0 + ts_yy_xxxyyz[i] * pa_x[i];

        ts_xyy_xxxyzz[i] = 3.0 * ts_yy_xxyzz[i] * fe_0 + ts_yy_xxxyzz[i] * pa_x[i];

        ts_xyy_xxxzzz[i] = 3.0 * ts_yy_xxzzz[i] * fe_0 + ts_yy_xxxzzz[i] * pa_x[i];

        ts_xyy_xxyyyy[i] = 2.0 * ts_yy_xyyyy[i] * fe_0 + ts_yy_xxyyyy[i] * pa_x[i];

        ts_xyy_xxyyyz[i] = 2.0 * ts_yy_xyyyz[i] * fe_0 + ts_yy_xxyyyz[i] * pa_x[i];

        ts_xyy_xxyyzz[i] = 2.0 * ts_yy_xyyzz[i] * fe_0 + ts_yy_xxyyzz[i] * pa_x[i];

        ts_xyy_xxyzzz[i] = 2.0 * ts_yy_xyzzz[i] * fe_0 + ts_yy_xxyzzz[i] * pa_x[i];

        ts_xyy_xxzzzz[i] = 2.0 * ts_yy_xzzzz[i] * fe_0 + ts_yy_xxzzzz[i] * pa_x[i];

        ts_xyy_xyyyyy[i] = ts_yy_yyyyy[i] * fe_0 + ts_yy_xyyyyy[i] * pa_x[i];

        ts_xyy_xyyyyz[i] = ts_yy_yyyyz[i] * fe_0 + ts_yy_xyyyyz[i] * pa_x[i];

        ts_xyy_xyyyzz[i] = ts_yy_yyyzz[i] * fe_0 + ts_yy_xyyyzz[i] * pa_x[i];

        ts_xyy_xyyzzz[i] = ts_yy_yyzzz[i] * fe_0 + ts_yy_xyyzzz[i] * pa_x[i];

        ts_xyy_xyzzzz[i] = ts_yy_yzzzz[i] * fe_0 + ts_yy_xyzzzz[i] * pa_x[i];

        ts_xyy_xzzzzz[i] = ts_yy_zzzzz[i] * fe_0 + ts_yy_xzzzzz[i] * pa_x[i];

        ts_xyy_yyyyyy[i] = ts_yy_yyyyyy[i] * pa_x[i];

        ts_xyy_yyyyyz[i] = ts_yy_yyyyyz[i] * pa_x[i];

        ts_xyy_yyyyzz[i] = ts_yy_yyyyzz[i] * pa_x[i];

        ts_xyy_yyyzzz[i] = ts_yy_yyyzzz[i] * pa_x[i];

        ts_xyy_yyzzzz[i] = ts_yy_yyzzzz[i] * pa_x[i];

        ts_xyy_yzzzzz[i] = ts_yy_yzzzzz[i] * pa_x[i];

        ts_xyy_zzzzzz[i] = ts_yy_zzzzzz[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             ts_xy_xxxxxy,  \
                             ts_xy_xxxxyy,  \
                             ts_xy_xxxyyy,  \
                             ts_xy_xxyyyy,  \
                             ts_xy_xyyyyy,  \
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
                             ts_xz_xxxxxx,  \
                             ts_xz_xxxxxz,  \
                             ts_xz_xxxxzz,  \
                             ts_xz_xxxzzz,  \
                             ts_xz_xxzzzz,  \
                             ts_xz_xzzzzz,  \
                             ts_yz_xxxxyz,  \
                             ts_yz_xxxyyz,  \
                             ts_yz_xxxyz,   \
                             ts_yz_xxxyzz,  \
                             ts_yz_xxyyyz,  \
                             ts_yz_xxyyz,   \
                             ts_yz_xxyyzz,  \
                             ts_yz_xxyzz,   \
                             ts_yz_xxyzzz,  \
                             ts_yz_xyyyyz,  \
                             ts_yz_xyyyz,   \
                             ts_yz_xyyyzz,  \
                             ts_yz_xyyzz,   \
                             ts_yz_xyyzzz,  \
                             ts_yz_xyzzz,   \
                             ts_yz_xyzzzz,  \
                             ts_yz_yyyyyy,  \
                             ts_yz_yyyyyz,  \
                             ts_yz_yyyyz,   \
                             ts_yz_yyyyzz,  \
                             ts_yz_yyyzz,   \
                             ts_yz_yyyzzz,  \
                             ts_yz_yyzzz,   \
                             ts_yz_yyzzzz,  \
                             ts_yz_yzzzz,   \
                             ts_yz_yzzzzz,  \
                             ts_yz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyz_xxxxxx[i] = ts_xz_xxxxxx[i] * pa_y[i];

        ts_xyz_xxxxxy[i] = ts_xy_xxxxxy[i] * pa_z[i];

        ts_xyz_xxxxxz[i] = ts_xz_xxxxxz[i] * pa_y[i];

        ts_xyz_xxxxyy[i] = ts_xy_xxxxyy[i] * pa_z[i];

        ts_xyz_xxxxyz[i] = 4.0 * ts_yz_xxxyz[i] * fe_0 + ts_yz_xxxxyz[i] * pa_x[i];

        ts_xyz_xxxxzz[i] = ts_xz_xxxxzz[i] * pa_y[i];

        ts_xyz_xxxyyy[i] = ts_xy_xxxyyy[i] * pa_z[i];

        ts_xyz_xxxyyz[i] = 3.0 * ts_yz_xxyyz[i] * fe_0 + ts_yz_xxxyyz[i] * pa_x[i];

        ts_xyz_xxxyzz[i] = 3.0 * ts_yz_xxyzz[i] * fe_0 + ts_yz_xxxyzz[i] * pa_x[i];

        ts_xyz_xxxzzz[i] = ts_xz_xxxzzz[i] * pa_y[i];

        ts_xyz_xxyyyy[i] = ts_xy_xxyyyy[i] * pa_z[i];

        ts_xyz_xxyyyz[i] = 2.0 * ts_yz_xyyyz[i] * fe_0 + ts_yz_xxyyyz[i] * pa_x[i];

        ts_xyz_xxyyzz[i] = 2.0 * ts_yz_xyyzz[i] * fe_0 + ts_yz_xxyyzz[i] * pa_x[i];

        ts_xyz_xxyzzz[i] = 2.0 * ts_yz_xyzzz[i] * fe_0 + ts_yz_xxyzzz[i] * pa_x[i];

        ts_xyz_xxzzzz[i] = ts_xz_xxzzzz[i] * pa_y[i];

        ts_xyz_xyyyyy[i] = ts_xy_xyyyyy[i] * pa_z[i];

        ts_xyz_xyyyyz[i] = ts_yz_yyyyz[i] * fe_0 + ts_yz_xyyyyz[i] * pa_x[i];

        ts_xyz_xyyyzz[i] = ts_yz_yyyzz[i] * fe_0 + ts_yz_xyyyzz[i] * pa_x[i];

        ts_xyz_xyyzzz[i] = ts_yz_yyzzz[i] * fe_0 + ts_yz_xyyzzz[i] * pa_x[i];

        ts_xyz_xyzzzz[i] = ts_yz_yzzzz[i] * fe_0 + ts_yz_xyzzzz[i] * pa_x[i];

        ts_xyz_xzzzzz[i] = ts_xz_xzzzzz[i] * pa_y[i];

        ts_xyz_yyyyyy[i] = ts_yz_yyyyyy[i] * pa_x[i];

        ts_xyz_yyyyyz[i] = ts_yz_yyyyyz[i] * pa_x[i];

        ts_xyz_yyyyzz[i] = ts_yz_yyyyzz[i] * pa_x[i];

        ts_xyz_yyyzzz[i] = ts_yz_yyyzzz[i] * pa_x[i];

        ts_xyz_yyzzzz[i] = ts_yz_yyzzzz[i] * pa_x[i];

        ts_xyz_yzzzzz[i] = ts_yz_yzzzzz[i] * pa_x[i];

        ts_xyz_zzzzzz[i] = ts_yz_zzzzzz[i] * pa_x[i];
    }

    // Set up 140-168 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_x,              \
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
                             ts_zz_xxxxx,   \
                             ts_zz_xxxxxx,  \
                             ts_zz_xxxxxy,  \
                             ts_zz_xxxxxz,  \
                             ts_zz_xxxxy,   \
                             ts_zz_xxxxyy,  \
                             ts_zz_xxxxyz,  \
                             ts_zz_xxxxz,   \
                             ts_zz_xxxxzz,  \
                             ts_zz_xxxyy,   \
                             ts_zz_xxxyyy,  \
                             ts_zz_xxxyyz,  \
                             ts_zz_xxxyz,   \
                             ts_zz_xxxyzz,  \
                             ts_zz_xxxzz,   \
                             ts_zz_xxxzzz,  \
                             ts_zz_xxyyy,   \
                             ts_zz_xxyyyy,  \
                             ts_zz_xxyyyz,  \
                             ts_zz_xxyyz,   \
                             ts_zz_xxyyzz,  \
                             ts_zz_xxyzz,   \
                             ts_zz_xxyzzz,  \
                             ts_zz_xxzzz,   \
                             ts_zz_xxzzzz,  \
                             ts_zz_xyyyy,   \
                             ts_zz_xyyyyy,  \
                             ts_zz_xyyyyz,  \
                             ts_zz_xyyyz,   \
                             ts_zz_xyyyzz,  \
                             ts_zz_xyyzz,   \
                             ts_zz_xyyzzz,  \
                             ts_zz_xyzzz,   \
                             ts_zz_xyzzzz,  \
                             ts_zz_xzzzz,   \
                             ts_zz_xzzzzz,  \
                             ts_zz_yyyyy,   \
                             ts_zz_yyyyyy,  \
                             ts_zz_yyyyyz,  \
                             ts_zz_yyyyz,   \
                             ts_zz_yyyyzz,  \
                             ts_zz_yyyzz,   \
                             ts_zz_yyyzzz,  \
                             ts_zz_yyzzz,   \
                             ts_zz_yyzzzz,  \
                             ts_zz_yzzzz,   \
                             ts_zz_yzzzzz,  \
                             ts_zz_zzzzz,   \
                             ts_zz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzz_xxxxxx[i] = 6.0 * ts_zz_xxxxx[i] * fe_0 + ts_zz_xxxxxx[i] * pa_x[i];

        ts_xzz_xxxxxy[i] = 5.0 * ts_zz_xxxxy[i] * fe_0 + ts_zz_xxxxxy[i] * pa_x[i];

        ts_xzz_xxxxxz[i] = 5.0 * ts_zz_xxxxz[i] * fe_0 + ts_zz_xxxxxz[i] * pa_x[i];

        ts_xzz_xxxxyy[i] = 4.0 * ts_zz_xxxyy[i] * fe_0 + ts_zz_xxxxyy[i] * pa_x[i];

        ts_xzz_xxxxyz[i] = 4.0 * ts_zz_xxxyz[i] * fe_0 + ts_zz_xxxxyz[i] * pa_x[i];

        ts_xzz_xxxxzz[i] = 4.0 * ts_zz_xxxzz[i] * fe_0 + ts_zz_xxxxzz[i] * pa_x[i];

        ts_xzz_xxxyyy[i] = 3.0 * ts_zz_xxyyy[i] * fe_0 + ts_zz_xxxyyy[i] * pa_x[i];

        ts_xzz_xxxyyz[i] = 3.0 * ts_zz_xxyyz[i] * fe_0 + ts_zz_xxxyyz[i] * pa_x[i];

        ts_xzz_xxxyzz[i] = 3.0 * ts_zz_xxyzz[i] * fe_0 + ts_zz_xxxyzz[i] * pa_x[i];

        ts_xzz_xxxzzz[i] = 3.0 * ts_zz_xxzzz[i] * fe_0 + ts_zz_xxxzzz[i] * pa_x[i];

        ts_xzz_xxyyyy[i] = 2.0 * ts_zz_xyyyy[i] * fe_0 + ts_zz_xxyyyy[i] * pa_x[i];

        ts_xzz_xxyyyz[i] = 2.0 * ts_zz_xyyyz[i] * fe_0 + ts_zz_xxyyyz[i] * pa_x[i];

        ts_xzz_xxyyzz[i] = 2.0 * ts_zz_xyyzz[i] * fe_0 + ts_zz_xxyyzz[i] * pa_x[i];

        ts_xzz_xxyzzz[i] = 2.0 * ts_zz_xyzzz[i] * fe_0 + ts_zz_xxyzzz[i] * pa_x[i];

        ts_xzz_xxzzzz[i] = 2.0 * ts_zz_xzzzz[i] * fe_0 + ts_zz_xxzzzz[i] * pa_x[i];

        ts_xzz_xyyyyy[i] = ts_zz_yyyyy[i] * fe_0 + ts_zz_xyyyyy[i] * pa_x[i];

        ts_xzz_xyyyyz[i] = ts_zz_yyyyz[i] * fe_0 + ts_zz_xyyyyz[i] * pa_x[i];

        ts_xzz_xyyyzz[i] = ts_zz_yyyzz[i] * fe_0 + ts_zz_xyyyzz[i] * pa_x[i];

        ts_xzz_xyyzzz[i] = ts_zz_yyzzz[i] * fe_0 + ts_zz_xyyzzz[i] * pa_x[i];

        ts_xzz_xyzzzz[i] = ts_zz_yzzzz[i] * fe_0 + ts_zz_xyzzzz[i] * pa_x[i];

        ts_xzz_xzzzzz[i] = ts_zz_zzzzz[i] * fe_0 + ts_zz_xzzzzz[i] * pa_x[i];

        ts_xzz_yyyyyy[i] = ts_zz_yyyyyy[i] * pa_x[i];

        ts_xzz_yyyyyz[i] = ts_zz_yyyyyz[i] * pa_x[i];

        ts_xzz_yyyyzz[i] = ts_zz_yyyyzz[i] * pa_x[i];

        ts_xzz_yyyzzz[i] = ts_zz_yyyzzz[i] * pa_x[i];

        ts_xzz_yyzzzz[i] = ts_zz_yyzzzz[i] * pa_x[i];

        ts_xzz_yzzzzz[i] = ts_zz_yzzzzz[i] * pa_x[i];

        ts_xzz_zzzzzz[i] = ts_zz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_y,              \
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
                             ts_yy_xxxxx,   \
                             ts_yy_xxxxxx,  \
                             ts_yy_xxxxxy,  \
                             ts_yy_xxxxxz,  \
                             ts_yy_xxxxy,   \
                             ts_yy_xxxxyy,  \
                             ts_yy_xxxxyz,  \
                             ts_yy_xxxxz,   \
                             ts_yy_xxxxzz,  \
                             ts_yy_xxxyy,   \
                             ts_yy_xxxyyy,  \
                             ts_yy_xxxyyz,  \
                             ts_yy_xxxyz,   \
                             ts_yy_xxxyzz,  \
                             ts_yy_xxxzz,   \
                             ts_yy_xxxzzz,  \
                             ts_yy_xxyyy,   \
                             ts_yy_xxyyyy,  \
                             ts_yy_xxyyyz,  \
                             ts_yy_xxyyz,   \
                             ts_yy_xxyyzz,  \
                             ts_yy_xxyzz,   \
                             ts_yy_xxyzzz,  \
                             ts_yy_xxzzz,   \
                             ts_yy_xxzzzz,  \
                             ts_yy_xyyyy,   \
                             ts_yy_xyyyyy,  \
                             ts_yy_xyyyyz,  \
                             ts_yy_xyyyz,   \
                             ts_yy_xyyyzz,  \
                             ts_yy_xyyzz,   \
                             ts_yy_xyyzzz,  \
                             ts_yy_xyzzz,   \
                             ts_yy_xyzzzz,  \
                             ts_yy_xzzzz,   \
                             ts_yy_xzzzzz,  \
                             ts_yy_yyyyy,   \
                             ts_yy_yyyyyy,  \
                             ts_yy_yyyyyz,  \
                             ts_yy_yyyyz,   \
                             ts_yy_yyyyzz,  \
                             ts_yy_yyyzz,   \
                             ts_yy_yyyzzz,  \
                             ts_yy_yyzzz,   \
                             ts_yy_yyzzzz,  \
                             ts_yy_yzzzz,   \
                             ts_yy_yzzzzz,  \
                             ts_yy_zzzzz,   \
                             ts_yy_zzzzzz,  \
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

        ts_yyy_xxxxxx[i] = 2.0 * ts_y_xxxxxx[i] * fe_0 + ts_yy_xxxxxx[i] * pa_y[i];

        ts_yyy_xxxxxy[i] = 2.0 * ts_y_xxxxxy[i] * fe_0 + ts_yy_xxxxx[i] * fe_0 + ts_yy_xxxxxy[i] * pa_y[i];

        ts_yyy_xxxxxz[i] = 2.0 * ts_y_xxxxxz[i] * fe_0 + ts_yy_xxxxxz[i] * pa_y[i];

        ts_yyy_xxxxyy[i] = 2.0 * ts_y_xxxxyy[i] * fe_0 + 2.0 * ts_yy_xxxxy[i] * fe_0 + ts_yy_xxxxyy[i] * pa_y[i];

        ts_yyy_xxxxyz[i] = 2.0 * ts_y_xxxxyz[i] * fe_0 + ts_yy_xxxxz[i] * fe_0 + ts_yy_xxxxyz[i] * pa_y[i];

        ts_yyy_xxxxzz[i] = 2.0 * ts_y_xxxxzz[i] * fe_0 + ts_yy_xxxxzz[i] * pa_y[i];

        ts_yyy_xxxyyy[i] = 2.0 * ts_y_xxxyyy[i] * fe_0 + 3.0 * ts_yy_xxxyy[i] * fe_0 + ts_yy_xxxyyy[i] * pa_y[i];

        ts_yyy_xxxyyz[i] = 2.0 * ts_y_xxxyyz[i] * fe_0 + 2.0 * ts_yy_xxxyz[i] * fe_0 + ts_yy_xxxyyz[i] * pa_y[i];

        ts_yyy_xxxyzz[i] = 2.0 * ts_y_xxxyzz[i] * fe_0 + ts_yy_xxxzz[i] * fe_0 + ts_yy_xxxyzz[i] * pa_y[i];

        ts_yyy_xxxzzz[i] = 2.0 * ts_y_xxxzzz[i] * fe_0 + ts_yy_xxxzzz[i] * pa_y[i];

        ts_yyy_xxyyyy[i] = 2.0 * ts_y_xxyyyy[i] * fe_0 + 4.0 * ts_yy_xxyyy[i] * fe_0 + ts_yy_xxyyyy[i] * pa_y[i];

        ts_yyy_xxyyyz[i] = 2.0 * ts_y_xxyyyz[i] * fe_0 + 3.0 * ts_yy_xxyyz[i] * fe_0 + ts_yy_xxyyyz[i] * pa_y[i];

        ts_yyy_xxyyzz[i] = 2.0 * ts_y_xxyyzz[i] * fe_0 + 2.0 * ts_yy_xxyzz[i] * fe_0 + ts_yy_xxyyzz[i] * pa_y[i];

        ts_yyy_xxyzzz[i] = 2.0 * ts_y_xxyzzz[i] * fe_0 + ts_yy_xxzzz[i] * fe_0 + ts_yy_xxyzzz[i] * pa_y[i];

        ts_yyy_xxzzzz[i] = 2.0 * ts_y_xxzzzz[i] * fe_0 + ts_yy_xxzzzz[i] * pa_y[i];

        ts_yyy_xyyyyy[i] = 2.0 * ts_y_xyyyyy[i] * fe_0 + 5.0 * ts_yy_xyyyy[i] * fe_0 + ts_yy_xyyyyy[i] * pa_y[i];

        ts_yyy_xyyyyz[i] = 2.0 * ts_y_xyyyyz[i] * fe_0 + 4.0 * ts_yy_xyyyz[i] * fe_0 + ts_yy_xyyyyz[i] * pa_y[i];

        ts_yyy_xyyyzz[i] = 2.0 * ts_y_xyyyzz[i] * fe_0 + 3.0 * ts_yy_xyyzz[i] * fe_0 + ts_yy_xyyyzz[i] * pa_y[i];

        ts_yyy_xyyzzz[i] = 2.0 * ts_y_xyyzzz[i] * fe_0 + 2.0 * ts_yy_xyzzz[i] * fe_0 + ts_yy_xyyzzz[i] * pa_y[i];

        ts_yyy_xyzzzz[i] = 2.0 * ts_y_xyzzzz[i] * fe_0 + ts_yy_xzzzz[i] * fe_0 + ts_yy_xyzzzz[i] * pa_y[i];

        ts_yyy_xzzzzz[i] = 2.0 * ts_y_xzzzzz[i] * fe_0 + ts_yy_xzzzzz[i] * pa_y[i];

        ts_yyy_yyyyyy[i] = 2.0 * ts_y_yyyyyy[i] * fe_0 + 6.0 * ts_yy_yyyyy[i] * fe_0 + ts_yy_yyyyyy[i] * pa_y[i];

        ts_yyy_yyyyyz[i] = 2.0 * ts_y_yyyyyz[i] * fe_0 + 5.0 * ts_yy_yyyyz[i] * fe_0 + ts_yy_yyyyyz[i] * pa_y[i];

        ts_yyy_yyyyzz[i] = 2.0 * ts_y_yyyyzz[i] * fe_0 + 4.0 * ts_yy_yyyzz[i] * fe_0 + ts_yy_yyyyzz[i] * pa_y[i];

        ts_yyy_yyyzzz[i] = 2.0 * ts_y_yyyzzz[i] * fe_0 + 3.0 * ts_yy_yyzzz[i] * fe_0 + ts_yy_yyyzzz[i] * pa_y[i];

        ts_yyy_yyzzzz[i] = 2.0 * ts_y_yyzzzz[i] * fe_0 + 2.0 * ts_yy_yzzzz[i] * fe_0 + ts_yy_yyzzzz[i] * pa_y[i];

        ts_yyy_yzzzzz[i] = 2.0 * ts_y_yzzzzz[i] * fe_0 + ts_yy_zzzzz[i] * fe_0 + ts_yy_yzzzzz[i] * pa_y[i];

        ts_yyy_zzzzzz[i] = 2.0 * ts_y_zzzzzz[i] * fe_0 + ts_yy_zzzzzz[i] * pa_y[i];
    }

    // Set up 196-224 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             ts_yy_xxxxxx,  \
                             ts_yy_xxxxxy,  \
                             ts_yy_xxxxy,   \
                             ts_yy_xxxxyy,  \
                             ts_yy_xxxxyz,  \
                             ts_yy_xxxyy,   \
                             ts_yy_xxxyyy,  \
                             ts_yy_xxxyyz,  \
                             ts_yy_xxxyz,   \
                             ts_yy_xxxyzz,  \
                             ts_yy_xxyyy,   \
                             ts_yy_xxyyyy,  \
                             ts_yy_xxyyyz,  \
                             ts_yy_xxyyz,   \
                             ts_yy_xxyyzz,  \
                             ts_yy_xxyzz,   \
                             ts_yy_xxyzzz,  \
                             ts_yy_xyyyy,   \
                             ts_yy_xyyyyy,  \
                             ts_yy_xyyyyz,  \
                             ts_yy_xyyyz,   \
                             ts_yy_xyyyzz,  \
                             ts_yy_xyyzz,   \
                             ts_yy_xyyzzz,  \
                             ts_yy_xyzzz,   \
                             ts_yy_xyzzzz,  \
                             ts_yy_yyyyy,   \
                             ts_yy_yyyyyy,  \
                             ts_yy_yyyyyz,  \
                             ts_yy_yyyyz,   \
                             ts_yy_yyyyzz,  \
                             ts_yy_yyyzz,   \
                             ts_yy_yyyzzz,  \
                             ts_yy_yyzzz,   \
                             ts_yy_yyzzzz,  \
                             ts_yy_yzzzz,   \
                             ts_yy_yzzzzz,  \
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
                             ts_yz_xxxxxz,  \
                             ts_yz_xxxxzz,  \
                             ts_yz_xxxzzz,  \
                             ts_yz_xxzzzz,  \
                             ts_yz_xzzzzz,  \
                             ts_yz_zzzzzz,  \
                             ts_z_xxxxxz,   \
                             ts_z_xxxxzz,   \
                             ts_z_xxxzzz,   \
                             ts_z_xxzzzz,   \
                             ts_z_xzzzzz,   \
                             ts_z_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyz_xxxxxx[i] = ts_yy_xxxxxx[i] * pa_z[i];

        ts_yyz_xxxxxy[i] = ts_yy_xxxxxy[i] * pa_z[i];

        ts_yyz_xxxxxz[i] = ts_z_xxxxxz[i] * fe_0 + ts_yz_xxxxxz[i] * pa_y[i];

        ts_yyz_xxxxyy[i] = ts_yy_xxxxyy[i] * pa_z[i];

        ts_yyz_xxxxyz[i] = ts_yy_xxxxy[i] * fe_0 + ts_yy_xxxxyz[i] * pa_z[i];

        ts_yyz_xxxxzz[i] = ts_z_xxxxzz[i] * fe_0 + ts_yz_xxxxzz[i] * pa_y[i];

        ts_yyz_xxxyyy[i] = ts_yy_xxxyyy[i] * pa_z[i];

        ts_yyz_xxxyyz[i] = ts_yy_xxxyy[i] * fe_0 + ts_yy_xxxyyz[i] * pa_z[i];

        ts_yyz_xxxyzz[i] = 2.0 * ts_yy_xxxyz[i] * fe_0 + ts_yy_xxxyzz[i] * pa_z[i];

        ts_yyz_xxxzzz[i] = ts_z_xxxzzz[i] * fe_0 + ts_yz_xxxzzz[i] * pa_y[i];

        ts_yyz_xxyyyy[i] = ts_yy_xxyyyy[i] * pa_z[i];

        ts_yyz_xxyyyz[i] = ts_yy_xxyyy[i] * fe_0 + ts_yy_xxyyyz[i] * pa_z[i];

        ts_yyz_xxyyzz[i] = 2.0 * ts_yy_xxyyz[i] * fe_0 + ts_yy_xxyyzz[i] * pa_z[i];

        ts_yyz_xxyzzz[i] = 3.0 * ts_yy_xxyzz[i] * fe_0 + ts_yy_xxyzzz[i] * pa_z[i];

        ts_yyz_xxzzzz[i] = ts_z_xxzzzz[i] * fe_0 + ts_yz_xxzzzz[i] * pa_y[i];

        ts_yyz_xyyyyy[i] = ts_yy_xyyyyy[i] * pa_z[i];

        ts_yyz_xyyyyz[i] = ts_yy_xyyyy[i] * fe_0 + ts_yy_xyyyyz[i] * pa_z[i];

        ts_yyz_xyyyzz[i] = 2.0 * ts_yy_xyyyz[i] * fe_0 + ts_yy_xyyyzz[i] * pa_z[i];

        ts_yyz_xyyzzz[i] = 3.0 * ts_yy_xyyzz[i] * fe_0 + ts_yy_xyyzzz[i] * pa_z[i];

        ts_yyz_xyzzzz[i] = 4.0 * ts_yy_xyzzz[i] * fe_0 + ts_yy_xyzzzz[i] * pa_z[i];

        ts_yyz_xzzzzz[i] = ts_z_xzzzzz[i] * fe_0 + ts_yz_xzzzzz[i] * pa_y[i];

        ts_yyz_yyyyyy[i] = ts_yy_yyyyyy[i] * pa_z[i];

        ts_yyz_yyyyyz[i] = ts_yy_yyyyy[i] * fe_0 + ts_yy_yyyyyz[i] * pa_z[i];

        ts_yyz_yyyyzz[i] = 2.0 * ts_yy_yyyyz[i] * fe_0 + ts_yy_yyyyzz[i] * pa_z[i];

        ts_yyz_yyyzzz[i] = 3.0 * ts_yy_yyyzz[i] * fe_0 + ts_yy_yyyzzz[i] * pa_z[i];

        ts_yyz_yyzzzz[i] = 4.0 * ts_yy_yyzzz[i] * fe_0 + ts_yy_yyzzzz[i] * pa_z[i];

        ts_yyz_yzzzzz[i] = 5.0 * ts_yy_yzzzz[i] * fe_0 + ts_yy_yzzzzz[i] * pa_z[i];

        ts_yyz_zzzzzz[i] = ts_z_zzzzzz[i] * fe_0 + ts_yz_zzzzzz[i] * pa_y[i];
    }

    // Set up 224-252 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_y,              \
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
                             ts_zz_xxxxx,   \
                             ts_zz_xxxxxx,  \
                             ts_zz_xxxxxy,  \
                             ts_zz_xxxxxz,  \
                             ts_zz_xxxxy,   \
                             ts_zz_xxxxyy,  \
                             ts_zz_xxxxyz,  \
                             ts_zz_xxxxz,   \
                             ts_zz_xxxxzz,  \
                             ts_zz_xxxyy,   \
                             ts_zz_xxxyyy,  \
                             ts_zz_xxxyyz,  \
                             ts_zz_xxxyz,   \
                             ts_zz_xxxyzz,  \
                             ts_zz_xxxzz,   \
                             ts_zz_xxxzzz,  \
                             ts_zz_xxyyy,   \
                             ts_zz_xxyyyy,  \
                             ts_zz_xxyyyz,  \
                             ts_zz_xxyyz,   \
                             ts_zz_xxyyzz,  \
                             ts_zz_xxyzz,   \
                             ts_zz_xxyzzz,  \
                             ts_zz_xxzzz,   \
                             ts_zz_xxzzzz,  \
                             ts_zz_xyyyy,   \
                             ts_zz_xyyyyy,  \
                             ts_zz_xyyyyz,  \
                             ts_zz_xyyyz,   \
                             ts_zz_xyyyzz,  \
                             ts_zz_xyyzz,   \
                             ts_zz_xyyzzz,  \
                             ts_zz_xyzzz,   \
                             ts_zz_xyzzzz,  \
                             ts_zz_xzzzz,   \
                             ts_zz_xzzzzz,  \
                             ts_zz_yyyyy,   \
                             ts_zz_yyyyyy,  \
                             ts_zz_yyyyyz,  \
                             ts_zz_yyyyz,   \
                             ts_zz_yyyyzz,  \
                             ts_zz_yyyzz,   \
                             ts_zz_yyyzzz,  \
                             ts_zz_yyzzz,   \
                             ts_zz_yyzzzz,  \
                             ts_zz_yzzzz,   \
                             ts_zz_yzzzzz,  \
                             ts_zz_zzzzz,   \
                             ts_zz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzz_xxxxxx[i] = ts_zz_xxxxxx[i] * pa_y[i];

        ts_yzz_xxxxxy[i] = ts_zz_xxxxx[i] * fe_0 + ts_zz_xxxxxy[i] * pa_y[i];

        ts_yzz_xxxxxz[i] = ts_zz_xxxxxz[i] * pa_y[i];

        ts_yzz_xxxxyy[i] = 2.0 * ts_zz_xxxxy[i] * fe_0 + ts_zz_xxxxyy[i] * pa_y[i];

        ts_yzz_xxxxyz[i] = ts_zz_xxxxz[i] * fe_0 + ts_zz_xxxxyz[i] * pa_y[i];

        ts_yzz_xxxxzz[i] = ts_zz_xxxxzz[i] * pa_y[i];

        ts_yzz_xxxyyy[i] = 3.0 * ts_zz_xxxyy[i] * fe_0 + ts_zz_xxxyyy[i] * pa_y[i];

        ts_yzz_xxxyyz[i] = 2.0 * ts_zz_xxxyz[i] * fe_0 + ts_zz_xxxyyz[i] * pa_y[i];

        ts_yzz_xxxyzz[i] = ts_zz_xxxzz[i] * fe_0 + ts_zz_xxxyzz[i] * pa_y[i];

        ts_yzz_xxxzzz[i] = ts_zz_xxxzzz[i] * pa_y[i];

        ts_yzz_xxyyyy[i] = 4.0 * ts_zz_xxyyy[i] * fe_0 + ts_zz_xxyyyy[i] * pa_y[i];

        ts_yzz_xxyyyz[i] = 3.0 * ts_zz_xxyyz[i] * fe_0 + ts_zz_xxyyyz[i] * pa_y[i];

        ts_yzz_xxyyzz[i] = 2.0 * ts_zz_xxyzz[i] * fe_0 + ts_zz_xxyyzz[i] * pa_y[i];

        ts_yzz_xxyzzz[i] = ts_zz_xxzzz[i] * fe_0 + ts_zz_xxyzzz[i] * pa_y[i];

        ts_yzz_xxzzzz[i] = ts_zz_xxzzzz[i] * pa_y[i];

        ts_yzz_xyyyyy[i] = 5.0 * ts_zz_xyyyy[i] * fe_0 + ts_zz_xyyyyy[i] * pa_y[i];

        ts_yzz_xyyyyz[i] = 4.0 * ts_zz_xyyyz[i] * fe_0 + ts_zz_xyyyyz[i] * pa_y[i];

        ts_yzz_xyyyzz[i] = 3.0 * ts_zz_xyyzz[i] * fe_0 + ts_zz_xyyyzz[i] * pa_y[i];

        ts_yzz_xyyzzz[i] = 2.0 * ts_zz_xyzzz[i] * fe_0 + ts_zz_xyyzzz[i] * pa_y[i];

        ts_yzz_xyzzzz[i] = ts_zz_xzzzz[i] * fe_0 + ts_zz_xyzzzz[i] * pa_y[i];

        ts_yzz_xzzzzz[i] = ts_zz_xzzzzz[i] * pa_y[i];

        ts_yzz_yyyyyy[i] = 6.0 * ts_zz_yyyyy[i] * fe_0 + ts_zz_yyyyyy[i] * pa_y[i];

        ts_yzz_yyyyyz[i] = 5.0 * ts_zz_yyyyz[i] * fe_0 + ts_zz_yyyyyz[i] * pa_y[i];

        ts_yzz_yyyyzz[i] = 4.0 * ts_zz_yyyzz[i] * fe_0 + ts_zz_yyyyzz[i] * pa_y[i];

        ts_yzz_yyyzzz[i] = 3.0 * ts_zz_yyzzz[i] * fe_0 + ts_zz_yyyzzz[i] * pa_y[i];

        ts_yzz_yyzzzz[i] = 2.0 * ts_zz_yzzzz[i] * fe_0 + ts_zz_yyzzzz[i] * pa_y[i];

        ts_yzz_yzzzzz[i] = ts_zz_zzzzz[i] * fe_0 + ts_zz_yzzzzz[i] * pa_y[i];

        ts_yzz_zzzzzz[i] = ts_zz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : FI

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

#pragma omp simd aligned(pa_z,              \
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
                             ts_zz_xxxxx,   \
                             ts_zz_xxxxxx,  \
                             ts_zz_xxxxxy,  \
                             ts_zz_xxxxxz,  \
                             ts_zz_xxxxy,   \
                             ts_zz_xxxxyy,  \
                             ts_zz_xxxxyz,  \
                             ts_zz_xxxxz,   \
                             ts_zz_xxxxzz,  \
                             ts_zz_xxxyy,   \
                             ts_zz_xxxyyy,  \
                             ts_zz_xxxyyz,  \
                             ts_zz_xxxyz,   \
                             ts_zz_xxxyzz,  \
                             ts_zz_xxxzz,   \
                             ts_zz_xxxzzz,  \
                             ts_zz_xxyyy,   \
                             ts_zz_xxyyyy,  \
                             ts_zz_xxyyyz,  \
                             ts_zz_xxyyz,   \
                             ts_zz_xxyyzz,  \
                             ts_zz_xxyzz,   \
                             ts_zz_xxyzzz,  \
                             ts_zz_xxzzz,   \
                             ts_zz_xxzzzz,  \
                             ts_zz_xyyyy,   \
                             ts_zz_xyyyyy,  \
                             ts_zz_xyyyyz,  \
                             ts_zz_xyyyz,   \
                             ts_zz_xyyyzz,  \
                             ts_zz_xyyzz,   \
                             ts_zz_xyyzzz,  \
                             ts_zz_xyzzz,   \
                             ts_zz_xyzzzz,  \
                             ts_zz_xzzzz,   \
                             ts_zz_xzzzzz,  \
                             ts_zz_yyyyy,   \
                             ts_zz_yyyyyy,  \
                             ts_zz_yyyyyz,  \
                             ts_zz_yyyyz,   \
                             ts_zz_yyyyzz,  \
                             ts_zz_yyyzz,   \
                             ts_zz_yyyzzz,  \
                             ts_zz_yyzzz,   \
                             ts_zz_yyzzzz,  \
                             ts_zz_yzzzz,   \
                             ts_zz_yzzzzz,  \
                             ts_zz_zzzzz,   \
                             ts_zz_zzzzzz,  \
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

        ts_zzz_xxxxxx[i] = 2.0 * ts_z_xxxxxx[i] * fe_0 + ts_zz_xxxxxx[i] * pa_z[i];

        ts_zzz_xxxxxy[i] = 2.0 * ts_z_xxxxxy[i] * fe_0 + ts_zz_xxxxxy[i] * pa_z[i];

        ts_zzz_xxxxxz[i] = 2.0 * ts_z_xxxxxz[i] * fe_0 + ts_zz_xxxxx[i] * fe_0 + ts_zz_xxxxxz[i] * pa_z[i];

        ts_zzz_xxxxyy[i] = 2.0 * ts_z_xxxxyy[i] * fe_0 + ts_zz_xxxxyy[i] * pa_z[i];

        ts_zzz_xxxxyz[i] = 2.0 * ts_z_xxxxyz[i] * fe_0 + ts_zz_xxxxy[i] * fe_0 + ts_zz_xxxxyz[i] * pa_z[i];

        ts_zzz_xxxxzz[i] = 2.0 * ts_z_xxxxzz[i] * fe_0 + 2.0 * ts_zz_xxxxz[i] * fe_0 + ts_zz_xxxxzz[i] * pa_z[i];

        ts_zzz_xxxyyy[i] = 2.0 * ts_z_xxxyyy[i] * fe_0 + ts_zz_xxxyyy[i] * pa_z[i];

        ts_zzz_xxxyyz[i] = 2.0 * ts_z_xxxyyz[i] * fe_0 + ts_zz_xxxyy[i] * fe_0 + ts_zz_xxxyyz[i] * pa_z[i];

        ts_zzz_xxxyzz[i] = 2.0 * ts_z_xxxyzz[i] * fe_0 + 2.0 * ts_zz_xxxyz[i] * fe_0 + ts_zz_xxxyzz[i] * pa_z[i];

        ts_zzz_xxxzzz[i] = 2.0 * ts_z_xxxzzz[i] * fe_0 + 3.0 * ts_zz_xxxzz[i] * fe_0 + ts_zz_xxxzzz[i] * pa_z[i];

        ts_zzz_xxyyyy[i] = 2.0 * ts_z_xxyyyy[i] * fe_0 + ts_zz_xxyyyy[i] * pa_z[i];

        ts_zzz_xxyyyz[i] = 2.0 * ts_z_xxyyyz[i] * fe_0 + ts_zz_xxyyy[i] * fe_0 + ts_zz_xxyyyz[i] * pa_z[i];

        ts_zzz_xxyyzz[i] = 2.0 * ts_z_xxyyzz[i] * fe_0 + 2.0 * ts_zz_xxyyz[i] * fe_0 + ts_zz_xxyyzz[i] * pa_z[i];

        ts_zzz_xxyzzz[i] = 2.0 * ts_z_xxyzzz[i] * fe_0 + 3.0 * ts_zz_xxyzz[i] * fe_0 + ts_zz_xxyzzz[i] * pa_z[i];

        ts_zzz_xxzzzz[i] = 2.0 * ts_z_xxzzzz[i] * fe_0 + 4.0 * ts_zz_xxzzz[i] * fe_0 + ts_zz_xxzzzz[i] * pa_z[i];

        ts_zzz_xyyyyy[i] = 2.0 * ts_z_xyyyyy[i] * fe_0 + ts_zz_xyyyyy[i] * pa_z[i];

        ts_zzz_xyyyyz[i] = 2.0 * ts_z_xyyyyz[i] * fe_0 + ts_zz_xyyyy[i] * fe_0 + ts_zz_xyyyyz[i] * pa_z[i];

        ts_zzz_xyyyzz[i] = 2.0 * ts_z_xyyyzz[i] * fe_0 + 2.0 * ts_zz_xyyyz[i] * fe_0 + ts_zz_xyyyzz[i] * pa_z[i];

        ts_zzz_xyyzzz[i] = 2.0 * ts_z_xyyzzz[i] * fe_0 + 3.0 * ts_zz_xyyzz[i] * fe_0 + ts_zz_xyyzzz[i] * pa_z[i];

        ts_zzz_xyzzzz[i] = 2.0 * ts_z_xyzzzz[i] * fe_0 + 4.0 * ts_zz_xyzzz[i] * fe_0 + ts_zz_xyzzzz[i] * pa_z[i];

        ts_zzz_xzzzzz[i] = 2.0 * ts_z_xzzzzz[i] * fe_0 + 5.0 * ts_zz_xzzzz[i] * fe_0 + ts_zz_xzzzzz[i] * pa_z[i];

        ts_zzz_yyyyyy[i] = 2.0 * ts_z_yyyyyy[i] * fe_0 + ts_zz_yyyyyy[i] * pa_z[i];

        ts_zzz_yyyyyz[i] = 2.0 * ts_z_yyyyyz[i] * fe_0 + ts_zz_yyyyy[i] * fe_0 + ts_zz_yyyyyz[i] * pa_z[i];

        ts_zzz_yyyyzz[i] = 2.0 * ts_z_yyyyzz[i] * fe_0 + 2.0 * ts_zz_yyyyz[i] * fe_0 + ts_zz_yyyyzz[i] * pa_z[i];

        ts_zzz_yyyzzz[i] = 2.0 * ts_z_yyyzzz[i] * fe_0 + 3.0 * ts_zz_yyyzz[i] * fe_0 + ts_zz_yyyzzz[i] * pa_z[i];

        ts_zzz_yyzzzz[i] = 2.0 * ts_z_yyzzzz[i] * fe_0 + 4.0 * ts_zz_yyzzz[i] * fe_0 + ts_zz_yyzzzz[i] * pa_z[i];

        ts_zzz_yzzzzz[i] = 2.0 * ts_z_yzzzzz[i] * fe_0 + 5.0 * ts_zz_yzzzz[i] * fe_0 + ts_zz_yzzzzz[i] * pa_z[i];

        ts_zzz_zzzzzz[i] = 2.0 * ts_z_zzzzzz[i] * fe_0 + 6.0 * ts_zz_zzzzz[i] * fe_0 + ts_zz_zzzzzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
