#include "OverlapPrimRecGI.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_gi(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_gi,
                     const size_t idx_ovl_di,
                     const size_t idx_ovl_fh,
                     const size_t idx_ovl_fi,
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

    auto ts_xy_yyyyyy = pbuffer.data(idx_ovl_di + 49);

    auto ts_xy_yyyyyz = pbuffer.data(idx_ovl_di + 50);

    auto ts_xy_yyyyzz = pbuffer.data(idx_ovl_di + 51);

    auto ts_xy_yyyzzz = pbuffer.data(idx_ovl_di + 52);

    auto ts_xy_yyzzzz = pbuffer.data(idx_ovl_di + 53);

    auto ts_xy_yzzzzz = pbuffer.data(idx_ovl_di + 54);

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

    auto ts_yz_xxxxzz = pbuffer.data(idx_ovl_di + 117);

    auto ts_yz_xxxzzz = pbuffer.data(idx_ovl_di + 121);

    auto ts_yz_xxzzzz = pbuffer.data(idx_ovl_di + 126);

    auto ts_yz_xzzzzz = pbuffer.data(idx_ovl_di + 132);

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

    // Set up components of auxiliary buffer : FH

    auto ts_xxx_xxxxx = pbuffer.data(idx_ovl_fh);

    auto ts_xxx_xxxxy = pbuffer.data(idx_ovl_fh + 1);

    auto ts_xxx_xxxxz = pbuffer.data(idx_ovl_fh + 2);

    auto ts_xxx_xxxyy = pbuffer.data(idx_ovl_fh + 3);

    auto ts_xxx_xxxyz = pbuffer.data(idx_ovl_fh + 4);

    auto ts_xxx_xxxzz = pbuffer.data(idx_ovl_fh + 5);

    auto ts_xxx_xxyyy = pbuffer.data(idx_ovl_fh + 6);

    auto ts_xxx_xxyyz = pbuffer.data(idx_ovl_fh + 7);

    auto ts_xxx_xxyzz = pbuffer.data(idx_ovl_fh + 8);

    auto ts_xxx_xxzzz = pbuffer.data(idx_ovl_fh + 9);

    auto ts_xxx_xyyyy = pbuffer.data(idx_ovl_fh + 10);

    auto ts_xxx_xyyyz = pbuffer.data(idx_ovl_fh + 11);

    auto ts_xxx_xyyzz = pbuffer.data(idx_ovl_fh + 12);

    auto ts_xxx_xyzzz = pbuffer.data(idx_ovl_fh + 13);

    auto ts_xxx_xzzzz = pbuffer.data(idx_ovl_fh + 14);

    auto ts_xxx_yyyyy = pbuffer.data(idx_ovl_fh + 15);

    auto ts_xxx_yyyyz = pbuffer.data(idx_ovl_fh + 16);

    auto ts_xxx_yyyzz = pbuffer.data(idx_ovl_fh + 17);

    auto ts_xxx_yyzzz = pbuffer.data(idx_ovl_fh + 18);

    auto ts_xxx_yzzzz = pbuffer.data(idx_ovl_fh + 19);

    auto ts_xxx_zzzzz = pbuffer.data(idx_ovl_fh + 20);

    auto ts_xxz_xxxxz = pbuffer.data(idx_ovl_fh + 44);

    auto ts_xxz_xxxyz = pbuffer.data(idx_ovl_fh + 46);

    auto ts_xxz_xxxzz = pbuffer.data(idx_ovl_fh + 47);

    auto ts_xxz_xxyyz = pbuffer.data(idx_ovl_fh + 49);

    auto ts_xxz_xxyzz = pbuffer.data(idx_ovl_fh + 50);

    auto ts_xxz_xxzzz = pbuffer.data(idx_ovl_fh + 51);

    auto ts_xxz_xyyyz = pbuffer.data(idx_ovl_fh + 53);

    auto ts_xxz_xyyzz = pbuffer.data(idx_ovl_fh + 54);

    auto ts_xxz_xyzzz = pbuffer.data(idx_ovl_fh + 55);

    auto ts_xxz_xzzzz = pbuffer.data(idx_ovl_fh + 56);

    auto ts_xyy_xxxxy = pbuffer.data(idx_ovl_fh + 64);

    auto ts_xyy_xxxyy = pbuffer.data(idx_ovl_fh + 66);

    auto ts_xyy_xxxyz = pbuffer.data(idx_ovl_fh + 67);

    auto ts_xyy_xxyyy = pbuffer.data(idx_ovl_fh + 69);

    auto ts_xyy_xxyyz = pbuffer.data(idx_ovl_fh + 70);

    auto ts_xyy_xxyzz = pbuffer.data(idx_ovl_fh + 71);

    auto ts_xyy_xyyyy = pbuffer.data(idx_ovl_fh + 73);

    auto ts_xyy_xyyyz = pbuffer.data(idx_ovl_fh + 74);

    auto ts_xyy_xyyzz = pbuffer.data(idx_ovl_fh + 75);

    auto ts_xyy_xyzzz = pbuffer.data(idx_ovl_fh + 76);

    auto ts_xyy_yyyyy = pbuffer.data(idx_ovl_fh + 78);

    auto ts_xyy_yyyyz = pbuffer.data(idx_ovl_fh + 79);

    auto ts_xyy_yyyzz = pbuffer.data(idx_ovl_fh + 80);

    auto ts_xyy_yyzzz = pbuffer.data(idx_ovl_fh + 81);

    auto ts_xyy_yzzzz = pbuffer.data(idx_ovl_fh + 82);

    auto ts_xzz_xxxxz = pbuffer.data(idx_ovl_fh + 107);

    auto ts_xzz_xxxyz = pbuffer.data(idx_ovl_fh + 109);

    auto ts_xzz_xxxzz = pbuffer.data(idx_ovl_fh + 110);

    auto ts_xzz_xxyyz = pbuffer.data(idx_ovl_fh + 112);

    auto ts_xzz_xxyzz = pbuffer.data(idx_ovl_fh + 113);

    auto ts_xzz_xxzzz = pbuffer.data(idx_ovl_fh + 114);

    auto ts_xzz_xyyyz = pbuffer.data(idx_ovl_fh + 116);

    auto ts_xzz_xyyzz = pbuffer.data(idx_ovl_fh + 117);

    auto ts_xzz_xyzzz = pbuffer.data(idx_ovl_fh + 118);

    auto ts_xzz_xzzzz = pbuffer.data(idx_ovl_fh + 119);

    auto ts_xzz_yyyyz = pbuffer.data(idx_ovl_fh + 121);

    auto ts_xzz_yyyzz = pbuffer.data(idx_ovl_fh + 122);

    auto ts_xzz_yyzzz = pbuffer.data(idx_ovl_fh + 123);

    auto ts_xzz_yzzzz = pbuffer.data(idx_ovl_fh + 124);

    auto ts_xzz_zzzzz = pbuffer.data(idx_ovl_fh + 125);

    auto ts_yyy_xxxxx = pbuffer.data(idx_ovl_fh + 126);

    auto ts_yyy_xxxxy = pbuffer.data(idx_ovl_fh + 127);

    auto ts_yyy_xxxxz = pbuffer.data(idx_ovl_fh + 128);

    auto ts_yyy_xxxyy = pbuffer.data(idx_ovl_fh + 129);

    auto ts_yyy_xxxyz = pbuffer.data(idx_ovl_fh + 130);

    auto ts_yyy_xxxzz = pbuffer.data(idx_ovl_fh + 131);

    auto ts_yyy_xxyyy = pbuffer.data(idx_ovl_fh + 132);

    auto ts_yyy_xxyyz = pbuffer.data(idx_ovl_fh + 133);

    auto ts_yyy_xxyzz = pbuffer.data(idx_ovl_fh + 134);

    auto ts_yyy_xxzzz = pbuffer.data(idx_ovl_fh + 135);

    auto ts_yyy_xyyyy = pbuffer.data(idx_ovl_fh + 136);

    auto ts_yyy_xyyyz = pbuffer.data(idx_ovl_fh + 137);

    auto ts_yyy_xyyzz = pbuffer.data(idx_ovl_fh + 138);

    auto ts_yyy_xyzzz = pbuffer.data(idx_ovl_fh + 139);

    auto ts_yyy_xzzzz = pbuffer.data(idx_ovl_fh + 140);

    auto ts_yyy_yyyyy = pbuffer.data(idx_ovl_fh + 141);

    auto ts_yyy_yyyyz = pbuffer.data(idx_ovl_fh + 142);

    auto ts_yyy_yyyzz = pbuffer.data(idx_ovl_fh + 143);

    auto ts_yyy_yyzzz = pbuffer.data(idx_ovl_fh + 144);

    auto ts_yyy_yzzzz = pbuffer.data(idx_ovl_fh + 145);

    auto ts_yyy_zzzzz = pbuffer.data(idx_ovl_fh + 146);

    auto ts_yyz_xxxxz = pbuffer.data(idx_ovl_fh + 149);

    auto ts_yyz_xxxyz = pbuffer.data(idx_ovl_fh + 151);

    auto ts_yyz_xxxzz = pbuffer.data(idx_ovl_fh + 152);

    auto ts_yyz_xxyyz = pbuffer.data(idx_ovl_fh + 154);

    auto ts_yyz_xxyzz = pbuffer.data(idx_ovl_fh + 155);

    auto ts_yyz_xxzzz = pbuffer.data(idx_ovl_fh + 156);

    auto ts_yyz_xyyyz = pbuffer.data(idx_ovl_fh + 158);

    auto ts_yyz_xyyzz = pbuffer.data(idx_ovl_fh + 159);

    auto ts_yyz_xyzzz = pbuffer.data(idx_ovl_fh + 160);

    auto ts_yyz_xzzzz = pbuffer.data(idx_ovl_fh + 161);

    auto ts_yyz_yyyyz = pbuffer.data(idx_ovl_fh + 163);

    auto ts_yyz_yyyzz = pbuffer.data(idx_ovl_fh + 164);

    auto ts_yyz_yyzzz = pbuffer.data(idx_ovl_fh + 165);

    auto ts_yyz_yzzzz = pbuffer.data(idx_ovl_fh + 166);

    auto ts_yyz_zzzzz = pbuffer.data(idx_ovl_fh + 167);

    auto ts_yzz_xxxxy = pbuffer.data(idx_ovl_fh + 169);

    auto ts_yzz_xxxxz = pbuffer.data(idx_ovl_fh + 170);

    auto ts_yzz_xxxyy = pbuffer.data(idx_ovl_fh + 171);

    auto ts_yzz_xxxyz = pbuffer.data(idx_ovl_fh + 172);

    auto ts_yzz_xxxzz = pbuffer.data(idx_ovl_fh + 173);

    auto ts_yzz_xxyyy = pbuffer.data(idx_ovl_fh + 174);

    auto ts_yzz_xxyyz = pbuffer.data(idx_ovl_fh + 175);

    auto ts_yzz_xxyzz = pbuffer.data(idx_ovl_fh + 176);

    auto ts_yzz_xxzzz = pbuffer.data(idx_ovl_fh + 177);

    auto ts_yzz_xyyyy = pbuffer.data(idx_ovl_fh + 178);

    auto ts_yzz_xyyyz = pbuffer.data(idx_ovl_fh + 179);

    auto ts_yzz_xyyzz = pbuffer.data(idx_ovl_fh + 180);

    auto ts_yzz_xyzzz = pbuffer.data(idx_ovl_fh + 181);

    auto ts_yzz_xzzzz = pbuffer.data(idx_ovl_fh + 182);

    auto ts_yzz_yyyyy = pbuffer.data(idx_ovl_fh + 183);

    auto ts_yzz_yyyyz = pbuffer.data(idx_ovl_fh + 184);

    auto ts_yzz_yyyzz = pbuffer.data(idx_ovl_fh + 185);

    auto ts_yzz_yyzzz = pbuffer.data(idx_ovl_fh + 186);

    auto ts_yzz_yzzzz = pbuffer.data(idx_ovl_fh + 187);

    auto ts_yzz_zzzzz = pbuffer.data(idx_ovl_fh + 188);

    auto ts_zzz_xxxxx = pbuffer.data(idx_ovl_fh + 189);

    auto ts_zzz_xxxxy = pbuffer.data(idx_ovl_fh + 190);

    auto ts_zzz_xxxxz = pbuffer.data(idx_ovl_fh + 191);

    auto ts_zzz_xxxyy = pbuffer.data(idx_ovl_fh + 192);

    auto ts_zzz_xxxyz = pbuffer.data(idx_ovl_fh + 193);

    auto ts_zzz_xxxzz = pbuffer.data(idx_ovl_fh + 194);

    auto ts_zzz_xxyyy = pbuffer.data(idx_ovl_fh + 195);

    auto ts_zzz_xxyyz = pbuffer.data(idx_ovl_fh + 196);

    auto ts_zzz_xxyzz = pbuffer.data(idx_ovl_fh + 197);

    auto ts_zzz_xxzzz = pbuffer.data(idx_ovl_fh + 198);

    auto ts_zzz_xyyyy = pbuffer.data(idx_ovl_fh + 199);

    auto ts_zzz_xyyyz = pbuffer.data(idx_ovl_fh + 200);

    auto ts_zzz_xyyzz = pbuffer.data(idx_ovl_fh + 201);

    auto ts_zzz_xyzzz = pbuffer.data(idx_ovl_fh + 202);

    auto ts_zzz_xzzzz = pbuffer.data(idx_ovl_fh + 203);

    auto ts_zzz_yyyyy = pbuffer.data(idx_ovl_fh + 204);

    auto ts_zzz_yyyyz = pbuffer.data(idx_ovl_fh + 205);

    auto ts_zzz_yyyzz = pbuffer.data(idx_ovl_fh + 206);

    auto ts_zzz_yyzzz = pbuffer.data(idx_ovl_fh + 207);

    auto ts_zzz_yzzzz = pbuffer.data(idx_ovl_fh + 208);

    auto ts_zzz_zzzzz = pbuffer.data(idx_ovl_fh + 209);

    // Set up components of auxiliary buffer : FI

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

    auto ts_xxy_xxxxxx = pbuffer.data(idx_ovl_fi + 28);

    auto ts_xxy_xxxxxy = pbuffer.data(idx_ovl_fi + 29);

    auto ts_xxy_xxxxxz = pbuffer.data(idx_ovl_fi + 30);

    auto ts_xxy_xxxxyy = pbuffer.data(idx_ovl_fi + 31);

    auto ts_xxy_xxxxzz = pbuffer.data(idx_ovl_fi + 33);

    auto ts_xxy_xxxyyy = pbuffer.data(idx_ovl_fi + 34);

    auto ts_xxy_xxxzzz = pbuffer.data(idx_ovl_fi + 37);

    auto ts_xxy_xxyyyy = pbuffer.data(idx_ovl_fi + 38);

    auto ts_xxy_xxzzzz = pbuffer.data(idx_ovl_fi + 42);

    auto ts_xxy_xyyyyy = pbuffer.data(idx_ovl_fi + 43);

    auto ts_xxy_xzzzzz = pbuffer.data(idx_ovl_fi + 48);

    auto ts_xxy_yyyyyy = pbuffer.data(idx_ovl_fi + 49);

    auto ts_xxy_yyyyyz = pbuffer.data(idx_ovl_fi + 50);

    auto ts_xxy_yyyyzz = pbuffer.data(idx_ovl_fi + 51);

    auto ts_xxy_yyyzzz = pbuffer.data(idx_ovl_fi + 52);

    auto ts_xxy_yyzzzz = pbuffer.data(idx_ovl_fi + 53);

    auto ts_xxy_yzzzzz = pbuffer.data(idx_ovl_fi + 54);

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

    auto ts_xxz_yyyyyz = pbuffer.data(idx_ovl_fi + 78);

    auto ts_xxz_yyyyzz = pbuffer.data(idx_ovl_fi + 79);

    auto ts_xxz_yyyzzz = pbuffer.data(idx_ovl_fi + 80);

    auto ts_xxz_yyzzzz = pbuffer.data(idx_ovl_fi + 81);

    auto ts_xxz_yzzzzz = pbuffer.data(idx_ovl_fi + 82);

    auto ts_xxz_zzzzzz = pbuffer.data(idx_ovl_fi + 83);

    auto ts_xyy_xxxxxx = pbuffer.data(idx_ovl_fi + 84);

    auto ts_xyy_xxxxxy = pbuffer.data(idx_ovl_fi + 85);

    auto ts_xyy_xxxxyy = pbuffer.data(idx_ovl_fi + 87);

    auto ts_xyy_xxxxyz = pbuffer.data(idx_ovl_fi + 88);

    auto ts_xyy_xxxyyy = pbuffer.data(idx_ovl_fi + 90);

    auto ts_xyy_xxxyyz = pbuffer.data(idx_ovl_fi + 91);

    auto ts_xyy_xxxyzz = pbuffer.data(idx_ovl_fi + 92);

    auto ts_xyy_xxyyyy = pbuffer.data(idx_ovl_fi + 94);

    auto ts_xyy_xxyyyz = pbuffer.data(idx_ovl_fi + 95);

    auto ts_xyy_xxyyzz = pbuffer.data(idx_ovl_fi + 96);

    auto ts_xyy_xxyzzz = pbuffer.data(idx_ovl_fi + 97);

    auto ts_xyy_xyyyyy = pbuffer.data(idx_ovl_fi + 99);

    auto ts_xyy_xyyyyz = pbuffer.data(idx_ovl_fi + 100);

    auto ts_xyy_xyyyzz = pbuffer.data(idx_ovl_fi + 101);

    auto ts_xyy_xyyzzz = pbuffer.data(idx_ovl_fi + 102);

    auto ts_xyy_xyzzzz = pbuffer.data(idx_ovl_fi + 103);

    auto ts_xyy_yyyyyy = pbuffer.data(idx_ovl_fi + 105);

    auto ts_xyy_yyyyyz = pbuffer.data(idx_ovl_fi + 106);

    auto ts_xyy_yyyyzz = pbuffer.data(idx_ovl_fi + 107);

    auto ts_xyy_yyyzzz = pbuffer.data(idx_ovl_fi + 108);

    auto ts_xyy_yyzzzz = pbuffer.data(idx_ovl_fi + 109);

    auto ts_xyy_yzzzzz = pbuffer.data(idx_ovl_fi + 110);

    auto ts_xyy_zzzzzz = pbuffer.data(idx_ovl_fi + 111);

    auto ts_xyz_yyyyyz = pbuffer.data(idx_ovl_fi + 134);

    auto ts_xyz_yyyyzz = pbuffer.data(idx_ovl_fi + 135);

    auto ts_xyz_yyyzzz = pbuffer.data(idx_ovl_fi + 136);

    auto ts_xyz_yyzzzz = pbuffer.data(idx_ovl_fi + 137);

    auto ts_xyz_yzzzzz = pbuffer.data(idx_ovl_fi + 138);

    auto ts_xzz_xxxxxx = pbuffer.data(idx_ovl_fi + 140);

    auto ts_xzz_xxxxxz = pbuffer.data(idx_ovl_fi + 142);

    auto ts_xzz_xxxxyz = pbuffer.data(idx_ovl_fi + 144);

    auto ts_xzz_xxxxzz = pbuffer.data(idx_ovl_fi + 145);

    auto ts_xzz_xxxyyz = pbuffer.data(idx_ovl_fi + 147);

    auto ts_xzz_xxxyzz = pbuffer.data(idx_ovl_fi + 148);

    auto ts_xzz_xxxzzz = pbuffer.data(idx_ovl_fi + 149);

    auto ts_xzz_xxyyyz = pbuffer.data(idx_ovl_fi + 151);

    auto ts_xzz_xxyyzz = pbuffer.data(idx_ovl_fi + 152);

    auto ts_xzz_xxyzzz = pbuffer.data(idx_ovl_fi + 153);

    auto ts_xzz_xxzzzz = pbuffer.data(idx_ovl_fi + 154);

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

    // Set up 0-28 components of targeted buffer : GI

    auto ts_xxxx_xxxxxx = pbuffer.data(idx_ovl_gi);

    auto ts_xxxx_xxxxxy = pbuffer.data(idx_ovl_gi + 1);

    auto ts_xxxx_xxxxxz = pbuffer.data(idx_ovl_gi + 2);

    auto ts_xxxx_xxxxyy = pbuffer.data(idx_ovl_gi + 3);

    auto ts_xxxx_xxxxyz = pbuffer.data(idx_ovl_gi + 4);

    auto ts_xxxx_xxxxzz = pbuffer.data(idx_ovl_gi + 5);

    auto ts_xxxx_xxxyyy = pbuffer.data(idx_ovl_gi + 6);

    auto ts_xxxx_xxxyyz = pbuffer.data(idx_ovl_gi + 7);

    auto ts_xxxx_xxxyzz = pbuffer.data(idx_ovl_gi + 8);

    auto ts_xxxx_xxxzzz = pbuffer.data(idx_ovl_gi + 9);

    auto ts_xxxx_xxyyyy = pbuffer.data(idx_ovl_gi + 10);

    auto ts_xxxx_xxyyyz = pbuffer.data(idx_ovl_gi + 11);

    auto ts_xxxx_xxyyzz = pbuffer.data(idx_ovl_gi + 12);

    auto ts_xxxx_xxyzzz = pbuffer.data(idx_ovl_gi + 13);

    auto ts_xxxx_xxzzzz = pbuffer.data(idx_ovl_gi + 14);

    auto ts_xxxx_xyyyyy = pbuffer.data(idx_ovl_gi + 15);

    auto ts_xxxx_xyyyyz = pbuffer.data(idx_ovl_gi + 16);

    auto ts_xxxx_xyyyzz = pbuffer.data(idx_ovl_gi + 17);

    auto ts_xxxx_xyyzzz = pbuffer.data(idx_ovl_gi + 18);

    auto ts_xxxx_xyzzzz = pbuffer.data(idx_ovl_gi + 19);

    auto ts_xxxx_xzzzzz = pbuffer.data(idx_ovl_gi + 20);

    auto ts_xxxx_yyyyyy = pbuffer.data(idx_ovl_gi + 21);

    auto ts_xxxx_yyyyyz = pbuffer.data(idx_ovl_gi + 22);

    auto ts_xxxx_yyyyzz = pbuffer.data(idx_ovl_gi + 23);

    auto ts_xxxx_yyyzzz = pbuffer.data(idx_ovl_gi + 24);

    auto ts_xxxx_yyzzzz = pbuffer.data(idx_ovl_gi + 25);

    auto ts_xxxx_yzzzzz = pbuffer.data(idx_ovl_gi + 26);

    auto ts_xxxx_zzzzzz = pbuffer.data(idx_ovl_gi + 27);

    #pragma omp simd aligned(pa_x, ts_xx_xxxxxx, ts_xx_xxxxxy, ts_xx_xxxxxz, ts_xx_xxxxyy, ts_xx_xxxxyz, ts_xx_xxxxzz, ts_xx_xxxyyy, ts_xx_xxxyyz, ts_xx_xxxyzz, ts_xx_xxxzzz, ts_xx_xxyyyy, ts_xx_xxyyyz, ts_xx_xxyyzz, ts_xx_xxyzzz, ts_xx_xxzzzz, ts_xx_xyyyyy, ts_xx_xyyyyz, ts_xx_xyyyzz, ts_xx_xyyzzz, ts_xx_xyzzzz, ts_xx_xzzzzz, ts_xx_yyyyyy, ts_xx_yyyyyz, ts_xx_yyyyzz, ts_xx_yyyzzz, ts_xx_yyzzzz, ts_xx_yzzzzz, ts_xx_zzzzzz, ts_xxx_xxxxx, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxxz, ts_xxx_xxxxy, ts_xxx_xxxxyy, ts_xxx_xxxxyz, ts_xxx_xxxxz, ts_xxx_xxxxzz, ts_xxx_xxxyy, ts_xxx_xxxyyy, ts_xxx_xxxyyz, ts_xxx_xxxyz, ts_xxx_xxxyzz, ts_xxx_xxxzz, ts_xxx_xxxzzz, ts_xxx_xxyyy, ts_xxx_xxyyyy, ts_xxx_xxyyyz, ts_xxx_xxyyz, ts_xxx_xxyyzz, ts_xxx_xxyzz, ts_xxx_xxyzzz, ts_xxx_xxzzz, ts_xxx_xxzzzz, ts_xxx_xyyyy, ts_xxx_xyyyyy, ts_xxx_xyyyyz, ts_xxx_xyyyz, ts_xxx_xyyyzz, ts_xxx_xyyzz, ts_xxx_xyyzzz, ts_xxx_xyzzz, ts_xxx_xyzzzz, ts_xxx_xzzzz, ts_xxx_xzzzzz, ts_xxx_yyyyy, ts_xxx_yyyyyy, ts_xxx_yyyyyz, ts_xxx_yyyyz, ts_xxx_yyyyzz, ts_xxx_yyyzz, ts_xxx_yyyzzz, ts_xxx_yyzzz, ts_xxx_yyzzzz, ts_xxx_yzzzz, ts_xxx_yzzzzz, ts_xxx_zzzzz, ts_xxx_zzzzzz, ts_xxxx_xxxxxx, ts_xxxx_xxxxxy, ts_xxxx_xxxxxz, ts_xxxx_xxxxyy, ts_xxxx_xxxxyz, ts_xxxx_xxxxzz, ts_xxxx_xxxyyy, ts_xxxx_xxxyyz, ts_xxxx_xxxyzz, ts_xxxx_xxxzzz, ts_xxxx_xxyyyy, ts_xxxx_xxyyyz, ts_xxxx_xxyyzz, ts_xxxx_xxyzzz, ts_xxxx_xxzzzz, ts_xxxx_xyyyyy, ts_xxxx_xyyyyz, ts_xxxx_xyyyzz, ts_xxxx_xyyzzz, ts_xxxx_xyzzzz, ts_xxxx_xzzzzz, ts_xxxx_yyyyyy, ts_xxxx_yyyyyz, ts_xxxx_yyyyzz, ts_xxxx_yyyzzz, ts_xxxx_yyzzzz, ts_xxxx_yzzzzz, ts_xxxx_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxx_xxxxxx[i] = 3.0 * ts_xx_xxxxxx[i] * fe_0 + 6.0 * ts_xxx_xxxxx[i] * fe_0 + ts_xxx_xxxxxx[i] * pa_x[i];

        ts_xxxx_xxxxxy[i] = 3.0 * ts_xx_xxxxxy[i] * fe_0 + 5.0 * ts_xxx_xxxxy[i] * fe_0 + ts_xxx_xxxxxy[i] * pa_x[i];

        ts_xxxx_xxxxxz[i] = 3.0 * ts_xx_xxxxxz[i] * fe_0 + 5.0 * ts_xxx_xxxxz[i] * fe_0 + ts_xxx_xxxxxz[i] * pa_x[i];

        ts_xxxx_xxxxyy[i] = 3.0 * ts_xx_xxxxyy[i] * fe_0 + 4.0 * ts_xxx_xxxyy[i] * fe_0 + ts_xxx_xxxxyy[i] * pa_x[i];

        ts_xxxx_xxxxyz[i] = 3.0 * ts_xx_xxxxyz[i] * fe_0 + 4.0 * ts_xxx_xxxyz[i] * fe_0 + ts_xxx_xxxxyz[i] * pa_x[i];

        ts_xxxx_xxxxzz[i] = 3.0 * ts_xx_xxxxzz[i] * fe_0 + 4.0 * ts_xxx_xxxzz[i] * fe_0 + ts_xxx_xxxxzz[i] * pa_x[i];

        ts_xxxx_xxxyyy[i] = 3.0 * ts_xx_xxxyyy[i] * fe_0 + 3.0 * ts_xxx_xxyyy[i] * fe_0 + ts_xxx_xxxyyy[i] * pa_x[i];

        ts_xxxx_xxxyyz[i] = 3.0 * ts_xx_xxxyyz[i] * fe_0 + 3.0 * ts_xxx_xxyyz[i] * fe_0 + ts_xxx_xxxyyz[i] * pa_x[i];

        ts_xxxx_xxxyzz[i] = 3.0 * ts_xx_xxxyzz[i] * fe_0 + 3.0 * ts_xxx_xxyzz[i] * fe_0 + ts_xxx_xxxyzz[i] * pa_x[i];

        ts_xxxx_xxxzzz[i] = 3.0 * ts_xx_xxxzzz[i] * fe_0 + 3.0 * ts_xxx_xxzzz[i] * fe_0 + ts_xxx_xxxzzz[i] * pa_x[i];

        ts_xxxx_xxyyyy[i] = 3.0 * ts_xx_xxyyyy[i] * fe_0 + 2.0 * ts_xxx_xyyyy[i] * fe_0 + ts_xxx_xxyyyy[i] * pa_x[i];

        ts_xxxx_xxyyyz[i] = 3.0 * ts_xx_xxyyyz[i] * fe_0 + 2.0 * ts_xxx_xyyyz[i] * fe_0 + ts_xxx_xxyyyz[i] * pa_x[i];

        ts_xxxx_xxyyzz[i] = 3.0 * ts_xx_xxyyzz[i] * fe_0 + 2.0 * ts_xxx_xyyzz[i] * fe_0 + ts_xxx_xxyyzz[i] * pa_x[i];

        ts_xxxx_xxyzzz[i] = 3.0 * ts_xx_xxyzzz[i] * fe_0 + 2.0 * ts_xxx_xyzzz[i] * fe_0 + ts_xxx_xxyzzz[i] * pa_x[i];

        ts_xxxx_xxzzzz[i] = 3.0 * ts_xx_xxzzzz[i] * fe_0 + 2.0 * ts_xxx_xzzzz[i] * fe_0 + ts_xxx_xxzzzz[i] * pa_x[i];

        ts_xxxx_xyyyyy[i] = 3.0 * ts_xx_xyyyyy[i] * fe_0 + ts_xxx_yyyyy[i] * fe_0 + ts_xxx_xyyyyy[i] * pa_x[i];

        ts_xxxx_xyyyyz[i] = 3.0 * ts_xx_xyyyyz[i] * fe_0 + ts_xxx_yyyyz[i] * fe_0 + ts_xxx_xyyyyz[i] * pa_x[i];

        ts_xxxx_xyyyzz[i] = 3.0 * ts_xx_xyyyzz[i] * fe_0 + ts_xxx_yyyzz[i] * fe_0 + ts_xxx_xyyyzz[i] * pa_x[i];

        ts_xxxx_xyyzzz[i] = 3.0 * ts_xx_xyyzzz[i] * fe_0 + ts_xxx_yyzzz[i] * fe_0 + ts_xxx_xyyzzz[i] * pa_x[i];

        ts_xxxx_xyzzzz[i] = 3.0 * ts_xx_xyzzzz[i] * fe_0 + ts_xxx_yzzzz[i] * fe_0 + ts_xxx_xyzzzz[i] * pa_x[i];

        ts_xxxx_xzzzzz[i] = 3.0 * ts_xx_xzzzzz[i] * fe_0 + ts_xxx_zzzzz[i] * fe_0 + ts_xxx_xzzzzz[i] * pa_x[i];

        ts_xxxx_yyyyyy[i] = 3.0 * ts_xx_yyyyyy[i] * fe_0 + ts_xxx_yyyyyy[i] * pa_x[i];

        ts_xxxx_yyyyyz[i] = 3.0 * ts_xx_yyyyyz[i] * fe_0 + ts_xxx_yyyyyz[i] * pa_x[i];

        ts_xxxx_yyyyzz[i] = 3.0 * ts_xx_yyyyzz[i] * fe_0 + ts_xxx_yyyyzz[i] * pa_x[i];

        ts_xxxx_yyyzzz[i] = 3.0 * ts_xx_yyyzzz[i] * fe_0 + ts_xxx_yyyzzz[i] * pa_x[i];

        ts_xxxx_yyzzzz[i] = 3.0 * ts_xx_yyzzzz[i] * fe_0 + ts_xxx_yyzzzz[i] * pa_x[i];

        ts_xxxx_yzzzzz[i] = 3.0 * ts_xx_yzzzzz[i] * fe_0 + ts_xxx_yzzzzz[i] * pa_x[i];

        ts_xxxx_zzzzzz[i] = 3.0 * ts_xx_zzzzzz[i] * fe_0 + ts_xxx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : GI

    auto ts_xxxy_xxxxxx = pbuffer.data(idx_ovl_gi + 28);

    auto ts_xxxy_xxxxxy = pbuffer.data(idx_ovl_gi + 29);

    auto ts_xxxy_xxxxxz = pbuffer.data(idx_ovl_gi + 30);

    auto ts_xxxy_xxxxyy = pbuffer.data(idx_ovl_gi + 31);

    auto ts_xxxy_xxxxyz = pbuffer.data(idx_ovl_gi + 32);

    auto ts_xxxy_xxxxzz = pbuffer.data(idx_ovl_gi + 33);

    auto ts_xxxy_xxxyyy = pbuffer.data(idx_ovl_gi + 34);

    auto ts_xxxy_xxxyyz = pbuffer.data(idx_ovl_gi + 35);

    auto ts_xxxy_xxxyzz = pbuffer.data(idx_ovl_gi + 36);

    auto ts_xxxy_xxxzzz = pbuffer.data(idx_ovl_gi + 37);

    auto ts_xxxy_xxyyyy = pbuffer.data(idx_ovl_gi + 38);

    auto ts_xxxy_xxyyyz = pbuffer.data(idx_ovl_gi + 39);

    auto ts_xxxy_xxyyzz = pbuffer.data(idx_ovl_gi + 40);

    auto ts_xxxy_xxyzzz = pbuffer.data(idx_ovl_gi + 41);

    auto ts_xxxy_xxzzzz = pbuffer.data(idx_ovl_gi + 42);

    auto ts_xxxy_xyyyyy = pbuffer.data(idx_ovl_gi + 43);

    auto ts_xxxy_xyyyyz = pbuffer.data(idx_ovl_gi + 44);

    auto ts_xxxy_xyyyzz = pbuffer.data(idx_ovl_gi + 45);

    auto ts_xxxy_xyyzzz = pbuffer.data(idx_ovl_gi + 46);

    auto ts_xxxy_xyzzzz = pbuffer.data(idx_ovl_gi + 47);

    auto ts_xxxy_xzzzzz = pbuffer.data(idx_ovl_gi + 48);

    auto ts_xxxy_yyyyyy = pbuffer.data(idx_ovl_gi + 49);

    auto ts_xxxy_yyyyyz = pbuffer.data(idx_ovl_gi + 50);

    auto ts_xxxy_yyyyzz = pbuffer.data(idx_ovl_gi + 51);

    auto ts_xxxy_yyyzzz = pbuffer.data(idx_ovl_gi + 52);

    auto ts_xxxy_yyzzzz = pbuffer.data(idx_ovl_gi + 53);

    auto ts_xxxy_yzzzzz = pbuffer.data(idx_ovl_gi + 54);

    auto ts_xxxy_zzzzzz = pbuffer.data(idx_ovl_gi + 55);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxx_xxxxx, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxxz, ts_xxx_xxxxy, ts_xxx_xxxxyy, ts_xxx_xxxxyz, ts_xxx_xxxxz, ts_xxx_xxxxzz, ts_xxx_xxxyy, ts_xxx_xxxyyy, ts_xxx_xxxyyz, ts_xxx_xxxyz, ts_xxx_xxxyzz, ts_xxx_xxxzz, ts_xxx_xxxzzz, ts_xxx_xxyyy, ts_xxx_xxyyyy, ts_xxx_xxyyyz, ts_xxx_xxyyz, ts_xxx_xxyyzz, ts_xxx_xxyzz, ts_xxx_xxyzzz, ts_xxx_xxzzz, ts_xxx_xxzzzz, ts_xxx_xyyyy, ts_xxx_xyyyyy, ts_xxx_xyyyyz, ts_xxx_xyyyz, ts_xxx_xyyyzz, ts_xxx_xyyzz, ts_xxx_xyyzzz, ts_xxx_xyzzz, ts_xxx_xyzzzz, ts_xxx_xzzzz, ts_xxx_xzzzzz, ts_xxx_zzzzzz, ts_xxxy_xxxxxx, ts_xxxy_xxxxxy, ts_xxxy_xxxxxz, ts_xxxy_xxxxyy, ts_xxxy_xxxxyz, ts_xxxy_xxxxzz, ts_xxxy_xxxyyy, ts_xxxy_xxxyyz, ts_xxxy_xxxyzz, ts_xxxy_xxxzzz, ts_xxxy_xxyyyy, ts_xxxy_xxyyyz, ts_xxxy_xxyyzz, ts_xxxy_xxyzzz, ts_xxxy_xxzzzz, ts_xxxy_xyyyyy, ts_xxxy_xyyyyz, ts_xxxy_xyyyzz, ts_xxxy_xyyzzz, ts_xxxy_xyzzzz, ts_xxxy_xzzzzz, ts_xxxy_yyyyyy, ts_xxxy_yyyyyz, ts_xxxy_yyyyzz, ts_xxxy_yyyzzz, ts_xxxy_yyzzzz, ts_xxxy_yzzzzz, ts_xxxy_zzzzzz, ts_xxy_yyyyyy, ts_xxy_yyyyyz, ts_xxy_yyyyzz, ts_xxy_yyyzzz, ts_xxy_yyzzzz, ts_xxy_yzzzzz, ts_xy_yyyyyy, ts_xy_yyyyyz, ts_xy_yyyyzz, ts_xy_yyyzzz, ts_xy_yyzzzz, ts_xy_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxy_xxxxxx[i] = ts_xxx_xxxxxx[i] * pa_y[i];

        ts_xxxy_xxxxxy[i] = ts_xxx_xxxxx[i] * fe_0 + ts_xxx_xxxxxy[i] * pa_y[i];

        ts_xxxy_xxxxxz[i] = ts_xxx_xxxxxz[i] * pa_y[i];

        ts_xxxy_xxxxyy[i] = 2.0 * ts_xxx_xxxxy[i] * fe_0 + ts_xxx_xxxxyy[i] * pa_y[i];

        ts_xxxy_xxxxyz[i] = ts_xxx_xxxxz[i] * fe_0 + ts_xxx_xxxxyz[i] * pa_y[i];

        ts_xxxy_xxxxzz[i] = ts_xxx_xxxxzz[i] * pa_y[i];

        ts_xxxy_xxxyyy[i] = 3.0 * ts_xxx_xxxyy[i] * fe_0 + ts_xxx_xxxyyy[i] * pa_y[i];

        ts_xxxy_xxxyyz[i] = 2.0 * ts_xxx_xxxyz[i] * fe_0 + ts_xxx_xxxyyz[i] * pa_y[i];

        ts_xxxy_xxxyzz[i] = ts_xxx_xxxzz[i] * fe_0 + ts_xxx_xxxyzz[i] * pa_y[i];

        ts_xxxy_xxxzzz[i] = ts_xxx_xxxzzz[i] * pa_y[i];

        ts_xxxy_xxyyyy[i] = 4.0 * ts_xxx_xxyyy[i] * fe_0 + ts_xxx_xxyyyy[i] * pa_y[i];

        ts_xxxy_xxyyyz[i] = 3.0 * ts_xxx_xxyyz[i] * fe_0 + ts_xxx_xxyyyz[i] * pa_y[i];

        ts_xxxy_xxyyzz[i] = 2.0 * ts_xxx_xxyzz[i] * fe_0 + ts_xxx_xxyyzz[i] * pa_y[i];

        ts_xxxy_xxyzzz[i] = ts_xxx_xxzzz[i] * fe_0 + ts_xxx_xxyzzz[i] * pa_y[i];

        ts_xxxy_xxzzzz[i] = ts_xxx_xxzzzz[i] * pa_y[i];

        ts_xxxy_xyyyyy[i] = 5.0 * ts_xxx_xyyyy[i] * fe_0 + ts_xxx_xyyyyy[i] * pa_y[i];

        ts_xxxy_xyyyyz[i] = 4.0 * ts_xxx_xyyyz[i] * fe_0 + ts_xxx_xyyyyz[i] * pa_y[i];

        ts_xxxy_xyyyzz[i] = 3.0 * ts_xxx_xyyzz[i] * fe_0 + ts_xxx_xyyyzz[i] * pa_y[i];

        ts_xxxy_xyyzzz[i] = 2.0 * ts_xxx_xyzzz[i] * fe_0 + ts_xxx_xyyzzz[i] * pa_y[i];

        ts_xxxy_xyzzzz[i] = ts_xxx_xzzzz[i] * fe_0 + ts_xxx_xyzzzz[i] * pa_y[i];

        ts_xxxy_xzzzzz[i] = ts_xxx_xzzzzz[i] * pa_y[i];

        ts_xxxy_yyyyyy[i] = 2.0 * ts_xy_yyyyyy[i] * fe_0 + ts_xxy_yyyyyy[i] * pa_x[i];

        ts_xxxy_yyyyyz[i] = 2.0 * ts_xy_yyyyyz[i] * fe_0 + ts_xxy_yyyyyz[i] * pa_x[i];

        ts_xxxy_yyyyzz[i] = 2.0 * ts_xy_yyyyzz[i] * fe_0 + ts_xxy_yyyyzz[i] * pa_x[i];

        ts_xxxy_yyyzzz[i] = 2.0 * ts_xy_yyyzzz[i] * fe_0 + ts_xxy_yyyzzz[i] * pa_x[i];

        ts_xxxy_yyzzzz[i] = 2.0 * ts_xy_yyzzzz[i] * fe_0 + ts_xxy_yyzzzz[i] * pa_x[i];

        ts_xxxy_yzzzzz[i] = 2.0 * ts_xy_yzzzzz[i] * fe_0 + ts_xxy_yzzzzz[i] * pa_x[i];

        ts_xxxy_zzzzzz[i] = ts_xxx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : GI

    auto ts_xxxz_xxxxxx = pbuffer.data(idx_ovl_gi + 56);

    auto ts_xxxz_xxxxxy = pbuffer.data(idx_ovl_gi + 57);

    auto ts_xxxz_xxxxxz = pbuffer.data(idx_ovl_gi + 58);

    auto ts_xxxz_xxxxyy = pbuffer.data(idx_ovl_gi + 59);

    auto ts_xxxz_xxxxyz = pbuffer.data(idx_ovl_gi + 60);

    auto ts_xxxz_xxxxzz = pbuffer.data(idx_ovl_gi + 61);

    auto ts_xxxz_xxxyyy = pbuffer.data(idx_ovl_gi + 62);

    auto ts_xxxz_xxxyyz = pbuffer.data(idx_ovl_gi + 63);

    auto ts_xxxz_xxxyzz = pbuffer.data(idx_ovl_gi + 64);

    auto ts_xxxz_xxxzzz = pbuffer.data(idx_ovl_gi + 65);

    auto ts_xxxz_xxyyyy = pbuffer.data(idx_ovl_gi + 66);

    auto ts_xxxz_xxyyyz = pbuffer.data(idx_ovl_gi + 67);

    auto ts_xxxz_xxyyzz = pbuffer.data(idx_ovl_gi + 68);

    auto ts_xxxz_xxyzzz = pbuffer.data(idx_ovl_gi + 69);

    auto ts_xxxz_xxzzzz = pbuffer.data(idx_ovl_gi + 70);

    auto ts_xxxz_xyyyyy = pbuffer.data(idx_ovl_gi + 71);

    auto ts_xxxz_xyyyyz = pbuffer.data(idx_ovl_gi + 72);

    auto ts_xxxz_xyyyzz = pbuffer.data(idx_ovl_gi + 73);

    auto ts_xxxz_xyyzzz = pbuffer.data(idx_ovl_gi + 74);

    auto ts_xxxz_xyzzzz = pbuffer.data(idx_ovl_gi + 75);

    auto ts_xxxz_xzzzzz = pbuffer.data(idx_ovl_gi + 76);

    auto ts_xxxz_yyyyyy = pbuffer.data(idx_ovl_gi + 77);

    auto ts_xxxz_yyyyyz = pbuffer.data(idx_ovl_gi + 78);

    auto ts_xxxz_yyyyzz = pbuffer.data(idx_ovl_gi + 79);

    auto ts_xxxz_yyyzzz = pbuffer.data(idx_ovl_gi + 80);

    auto ts_xxxz_yyzzzz = pbuffer.data(idx_ovl_gi + 81);

    auto ts_xxxz_yzzzzz = pbuffer.data(idx_ovl_gi + 82);

    auto ts_xxxz_zzzzzz = pbuffer.data(idx_ovl_gi + 83);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxx_xxxxx, ts_xxx_xxxxxx, ts_xxx_xxxxxy, ts_xxx_xxxxxz, ts_xxx_xxxxy, ts_xxx_xxxxyy, ts_xxx_xxxxyz, ts_xxx_xxxxz, ts_xxx_xxxxzz, ts_xxx_xxxyy, ts_xxx_xxxyyy, ts_xxx_xxxyyz, ts_xxx_xxxyz, ts_xxx_xxxyzz, ts_xxx_xxxzz, ts_xxx_xxxzzz, ts_xxx_xxyyy, ts_xxx_xxyyyy, ts_xxx_xxyyyz, ts_xxx_xxyyz, ts_xxx_xxyyzz, ts_xxx_xxyzz, ts_xxx_xxyzzz, ts_xxx_xxzzz, ts_xxx_xxzzzz, ts_xxx_xyyyy, ts_xxx_xyyyyy, ts_xxx_xyyyyz, ts_xxx_xyyyz, ts_xxx_xyyyzz, ts_xxx_xyyzz, ts_xxx_xyyzzz, ts_xxx_xyzzz, ts_xxx_xyzzzz, ts_xxx_xzzzz, ts_xxx_xzzzzz, ts_xxx_yyyyyy, ts_xxxz_xxxxxx, ts_xxxz_xxxxxy, ts_xxxz_xxxxxz, ts_xxxz_xxxxyy, ts_xxxz_xxxxyz, ts_xxxz_xxxxzz, ts_xxxz_xxxyyy, ts_xxxz_xxxyyz, ts_xxxz_xxxyzz, ts_xxxz_xxxzzz, ts_xxxz_xxyyyy, ts_xxxz_xxyyyz, ts_xxxz_xxyyzz, ts_xxxz_xxyzzz, ts_xxxz_xxzzzz, ts_xxxz_xyyyyy, ts_xxxz_xyyyyz, ts_xxxz_xyyyzz, ts_xxxz_xyyzzz, ts_xxxz_xyzzzz, ts_xxxz_xzzzzz, ts_xxxz_yyyyyy, ts_xxxz_yyyyyz, ts_xxxz_yyyyzz, ts_xxxz_yyyzzz, ts_xxxz_yyzzzz, ts_xxxz_yzzzzz, ts_xxxz_zzzzzz, ts_xxz_yyyyyz, ts_xxz_yyyyzz, ts_xxz_yyyzzz, ts_xxz_yyzzzz, ts_xxz_yzzzzz, ts_xxz_zzzzzz, ts_xz_yyyyyz, ts_xz_yyyyzz, ts_xz_yyyzzz, ts_xz_yyzzzz, ts_xz_yzzzzz, ts_xz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxz_xxxxxx[i] = ts_xxx_xxxxxx[i] * pa_z[i];

        ts_xxxz_xxxxxy[i] = ts_xxx_xxxxxy[i] * pa_z[i];

        ts_xxxz_xxxxxz[i] = ts_xxx_xxxxx[i] * fe_0 + ts_xxx_xxxxxz[i] * pa_z[i];

        ts_xxxz_xxxxyy[i] = ts_xxx_xxxxyy[i] * pa_z[i];

        ts_xxxz_xxxxyz[i] = ts_xxx_xxxxy[i] * fe_0 + ts_xxx_xxxxyz[i] * pa_z[i];

        ts_xxxz_xxxxzz[i] = 2.0 * ts_xxx_xxxxz[i] * fe_0 + ts_xxx_xxxxzz[i] * pa_z[i];

        ts_xxxz_xxxyyy[i] = ts_xxx_xxxyyy[i] * pa_z[i];

        ts_xxxz_xxxyyz[i] = ts_xxx_xxxyy[i] * fe_0 + ts_xxx_xxxyyz[i] * pa_z[i];

        ts_xxxz_xxxyzz[i] = 2.0 * ts_xxx_xxxyz[i] * fe_0 + ts_xxx_xxxyzz[i] * pa_z[i];

        ts_xxxz_xxxzzz[i] = 3.0 * ts_xxx_xxxzz[i] * fe_0 + ts_xxx_xxxzzz[i] * pa_z[i];

        ts_xxxz_xxyyyy[i] = ts_xxx_xxyyyy[i] * pa_z[i];

        ts_xxxz_xxyyyz[i] = ts_xxx_xxyyy[i] * fe_0 + ts_xxx_xxyyyz[i] * pa_z[i];

        ts_xxxz_xxyyzz[i] = 2.0 * ts_xxx_xxyyz[i] * fe_0 + ts_xxx_xxyyzz[i] * pa_z[i];

        ts_xxxz_xxyzzz[i] = 3.0 * ts_xxx_xxyzz[i] * fe_0 + ts_xxx_xxyzzz[i] * pa_z[i];

        ts_xxxz_xxzzzz[i] = 4.0 * ts_xxx_xxzzz[i] * fe_0 + ts_xxx_xxzzzz[i] * pa_z[i];

        ts_xxxz_xyyyyy[i] = ts_xxx_xyyyyy[i] * pa_z[i];

        ts_xxxz_xyyyyz[i] = ts_xxx_xyyyy[i] * fe_0 + ts_xxx_xyyyyz[i] * pa_z[i];

        ts_xxxz_xyyyzz[i] = 2.0 * ts_xxx_xyyyz[i] * fe_0 + ts_xxx_xyyyzz[i] * pa_z[i];

        ts_xxxz_xyyzzz[i] = 3.0 * ts_xxx_xyyzz[i] * fe_0 + ts_xxx_xyyzzz[i] * pa_z[i];

        ts_xxxz_xyzzzz[i] = 4.0 * ts_xxx_xyzzz[i] * fe_0 + ts_xxx_xyzzzz[i] * pa_z[i];

        ts_xxxz_xzzzzz[i] = 5.0 * ts_xxx_xzzzz[i] * fe_0 + ts_xxx_xzzzzz[i] * pa_z[i];

        ts_xxxz_yyyyyy[i] = ts_xxx_yyyyyy[i] * pa_z[i];

        ts_xxxz_yyyyyz[i] = 2.0 * ts_xz_yyyyyz[i] * fe_0 + ts_xxz_yyyyyz[i] * pa_x[i];

        ts_xxxz_yyyyzz[i] = 2.0 * ts_xz_yyyyzz[i] * fe_0 + ts_xxz_yyyyzz[i] * pa_x[i];

        ts_xxxz_yyyzzz[i] = 2.0 * ts_xz_yyyzzz[i] * fe_0 + ts_xxz_yyyzzz[i] * pa_x[i];

        ts_xxxz_yyzzzz[i] = 2.0 * ts_xz_yyzzzz[i] * fe_0 + ts_xxz_yyzzzz[i] * pa_x[i];

        ts_xxxz_yzzzzz[i] = 2.0 * ts_xz_yzzzzz[i] * fe_0 + ts_xxz_yzzzzz[i] * pa_x[i];

        ts_xxxz_zzzzzz[i] = 2.0 * ts_xz_zzzzzz[i] * fe_0 + ts_xxz_zzzzzz[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : GI

    auto ts_xxyy_xxxxxx = pbuffer.data(idx_ovl_gi + 84);

    auto ts_xxyy_xxxxxy = pbuffer.data(idx_ovl_gi + 85);

    auto ts_xxyy_xxxxxz = pbuffer.data(idx_ovl_gi + 86);

    auto ts_xxyy_xxxxyy = pbuffer.data(idx_ovl_gi + 87);

    auto ts_xxyy_xxxxyz = pbuffer.data(idx_ovl_gi + 88);

    auto ts_xxyy_xxxxzz = pbuffer.data(idx_ovl_gi + 89);

    auto ts_xxyy_xxxyyy = pbuffer.data(idx_ovl_gi + 90);

    auto ts_xxyy_xxxyyz = pbuffer.data(idx_ovl_gi + 91);

    auto ts_xxyy_xxxyzz = pbuffer.data(idx_ovl_gi + 92);

    auto ts_xxyy_xxxzzz = pbuffer.data(idx_ovl_gi + 93);

    auto ts_xxyy_xxyyyy = pbuffer.data(idx_ovl_gi + 94);

    auto ts_xxyy_xxyyyz = pbuffer.data(idx_ovl_gi + 95);

    auto ts_xxyy_xxyyzz = pbuffer.data(idx_ovl_gi + 96);

    auto ts_xxyy_xxyzzz = pbuffer.data(idx_ovl_gi + 97);

    auto ts_xxyy_xxzzzz = pbuffer.data(idx_ovl_gi + 98);

    auto ts_xxyy_xyyyyy = pbuffer.data(idx_ovl_gi + 99);

    auto ts_xxyy_xyyyyz = pbuffer.data(idx_ovl_gi + 100);

    auto ts_xxyy_xyyyzz = pbuffer.data(idx_ovl_gi + 101);

    auto ts_xxyy_xyyzzz = pbuffer.data(idx_ovl_gi + 102);

    auto ts_xxyy_xyzzzz = pbuffer.data(idx_ovl_gi + 103);

    auto ts_xxyy_xzzzzz = pbuffer.data(idx_ovl_gi + 104);

    auto ts_xxyy_yyyyyy = pbuffer.data(idx_ovl_gi + 105);

    auto ts_xxyy_yyyyyz = pbuffer.data(idx_ovl_gi + 106);

    auto ts_xxyy_yyyyzz = pbuffer.data(idx_ovl_gi + 107);

    auto ts_xxyy_yyyzzz = pbuffer.data(idx_ovl_gi + 108);

    auto ts_xxyy_yyzzzz = pbuffer.data(idx_ovl_gi + 109);

    auto ts_xxyy_yzzzzz = pbuffer.data(idx_ovl_gi + 110);

    auto ts_xxyy_zzzzzz = pbuffer.data(idx_ovl_gi + 111);

    #pragma omp simd aligned(pa_x, pa_y, ts_xx_xxxxxx, ts_xx_xxxxxz, ts_xx_xxxxzz, ts_xx_xxxzzz, ts_xx_xxzzzz, ts_xx_xzzzzz, ts_xxy_xxxxxx, ts_xxy_xxxxxz, ts_xxy_xxxxzz, ts_xxy_xxxzzz, ts_xxy_xxzzzz, ts_xxy_xzzzzz, ts_xxyy_xxxxxx, ts_xxyy_xxxxxy, ts_xxyy_xxxxxz, ts_xxyy_xxxxyy, ts_xxyy_xxxxyz, ts_xxyy_xxxxzz, ts_xxyy_xxxyyy, ts_xxyy_xxxyyz, ts_xxyy_xxxyzz, ts_xxyy_xxxzzz, ts_xxyy_xxyyyy, ts_xxyy_xxyyyz, ts_xxyy_xxyyzz, ts_xxyy_xxyzzz, ts_xxyy_xxzzzz, ts_xxyy_xyyyyy, ts_xxyy_xyyyyz, ts_xxyy_xyyyzz, ts_xxyy_xyyzzz, ts_xxyy_xyzzzz, ts_xxyy_xzzzzz, ts_xxyy_yyyyyy, ts_xxyy_yyyyyz, ts_xxyy_yyyyzz, ts_xxyy_yyyzzz, ts_xxyy_yyzzzz, ts_xxyy_yzzzzz, ts_xxyy_zzzzzz, ts_xyy_xxxxxy, ts_xyy_xxxxy, ts_xyy_xxxxyy, ts_xyy_xxxxyz, ts_xyy_xxxyy, ts_xyy_xxxyyy, ts_xyy_xxxyyz, ts_xyy_xxxyz, ts_xyy_xxxyzz, ts_xyy_xxyyy, ts_xyy_xxyyyy, ts_xyy_xxyyyz, ts_xyy_xxyyz, ts_xyy_xxyyzz, ts_xyy_xxyzz, ts_xyy_xxyzzz, ts_xyy_xyyyy, ts_xyy_xyyyyy, ts_xyy_xyyyyz, ts_xyy_xyyyz, ts_xyy_xyyyzz, ts_xyy_xyyzz, ts_xyy_xyyzzz, ts_xyy_xyzzz, ts_xyy_xyzzzz, ts_xyy_yyyyy, ts_xyy_yyyyyy, ts_xyy_yyyyyz, ts_xyy_yyyyz, ts_xyy_yyyyzz, ts_xyy_yyyzz, ts_xyy_yyyzzz, ts_xyy_yyzzz, ts_xyy_yyzzzz, ts_xyy_yzzzz, ts_xyy_yzzzzz, ts_xyy_zzzzzz, ts_yy_xxxxxy, ts_yy_xxxxyy, ts_yy_xxxxyz, ts_yy_xxxyyy, ts_yy_xxxyyz, ts_yy_xxxyzz, ts_yy_xxyyyy, ts_yy_xxyyyz, ts_yy_xxyyzz, ts_yy_xxyzzz, ts_yy_xyyyyy, ts_yy_xyyyyz, ts_yy_xyyyzz, ts_yy_xyyzzz, ts_yy_xyzzzz, ts_yy_yyyyyy, ts_yy_yyyyyz, ts_yy_yyyyzz, ts_yy_yyyzzz, ts_yy_yyzzzz, ts_yy_yzzzzz, ts_yy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyy_xxxxxx[i] = ts_xx_xxxxxx[i] * fe_0 + ts_xxy_xxxxxx[i] * pa_y[i];

        ts_xxyy_xxxxxy[i] = ts_yy_xxxxxy[i] * fe_0 + 5.0 * ts_xyy_xxxxy[i] * fe_0 + ts_xyy_xxxxxy[i] * pa_x[i];

        ts_xxyy_xxxxxz[i] = ts_xx_xxxxxz[i] * fe_0 + ts_xxy_xxxxxz[i] * pa_y[i];

        ts_xxyy_xxxxyy[i] = ts_yy_xxxxyy[i] * fe_0 + 4.0 * ts_xyy_xxxyy[i] * fe_0 + ts_xyy_xxxxyy[i] * pa_x[i];

        ts_xxyy_xxxxyz[i] = ts_yy_xxxxyz[i] * fe_0 + 4.0 * ts_xyy_xxxyz[i] * fe_0 + ts_xyy_xxxxyz[i] * pa_x[i];

        ts_xxyy_xxxxzz[i] = ts_xx_xxxxzz[i] * fe_0 + ts_xxy_xxxxzz[i] * pa_y[i];

        ts_xxyy_xxxyyy[i] = ts_yy_xxxyyy[i] * fe_0 + 3.0 * ts_xyy_xxyyy[i] * fe_0 + ts_xyy_xxxyyy[i] * pa_x[i];

        ts_xxyy_xxxyyz[i] = ts_yy_xxxyyz[i] * fe_0 + 3.0 * ts_xyy_xxyyz[i] * fe_0 + ts_xyy_xxxyyz[i] * pa_x[i];

        ts_xxyy_xxxyzz[i] = ts_yy_xxxyzz[i] * fe_0 + 3.0 * ts_xyy_xxyzz[i] * fe_0 + ts_xyy_xxxyzz[i] * pa_x[i];

        ts_xxyy_xxxzzz[i] = ts_xx_xxxzzz[i] * fe_0 + ts_xxy_xxxzzz[i] * pa_y[i];

        ts_xxyy_xxyyyy[i] = ts_yy_xxyyyy[i] * fe_0 + 2.0 * ts_xyy_xyyyy[i] * fe_0 + ts_xyy_xxyyyy[i] * pa_x[i];

        ts_xxyy_xxyyyz[i] = ts_yy_xxyyyz[i] * fe_0 + 2.0 * ts_xyy_xyyyz[i] * fe_0 + ts_xyy_xxyyyz[i] * pa_x[i];

        ts_xxyy_xxyyzz[i] = ts_yy_xxyyzz[i] * fe_0 + 2.0 * ts_xyy_xyyzz[i] * fe_0 + ts_xyy_xxyyzz[i] * pa_x[i];

        ts_xxyy_xxyzzz[i] = ts_yy_xxyzzz[i] * fe_0 + 2.0 * ts_xyy_xyzzz[i] * fe_0 + ts_xyy_xxyzzz[i] * pa_x[i];

        ts_xxyy_xxzzzz[i] = ts_xx_xxzzzz[i] * fe_0 + ts_xxy_xxzzzz[i] * pa_y[i];

        ts_xxyy_xyyyyy[i] = ts_yy_xyyyyy[i] * fe_0 + ts_xyy_yyyyy[i] * fe_0 + ts_xyy_xyyyyy[i] * pa_x[i];

        ts_xxyy_xyyyyz[i] = ts_yy_xyyyyz[i] * fe_0 + ts_xyy_yyyyz[i] * fe_0 + ts_xyy_xyyyyz[i] * pa_x[i];

        ts_xxyy_xyyyzz[i] = ts_yy_xyyyzz[i] * fe_0 + ts_xyy_yyyzz[i] * fe_0 + ts_xyy_xyyyzz[i] * pa_x[i];

        ts_xxyy_xyyzzz[i] = ts_yy_xyyzzz[i] * fe_0 + ts_xyy_yyzzz[i] * fe_0 + ts_xyy_xyyzzz[i] * pa_x[i];

        ts_xxyy_xyzzzz[i] = ts_yy_xyzzzz[i] * fe_0 + ts_xyy_yzzzz[i] * fe_0 + ts_xyy_xyzzzz[i] * pa_x[i];

        ts_xxyy_xzzzzz[i] = ts_xx_xzzzzz[i] * fe_0 + ts_xxy_xzzzzz[i] * pa_y[i];

        ts_xxyy_yyyyyy[i] = ts_yy_yyyyyy[i] * fe_0 + ts_xyy_yyyyyy[i] * pa_x[i];

        ts_xxyy_yyyyyz[i] = ts_yy_yyyyyz[i] * fe_0 + ts_xyy_yyyyyz[i] * pa_x[i];

        ts_xxyy_yyyyzz[i] = ts_yy_yyyyzz[i] * fe_0 + ts_xyy_yyyyzz[i] * pa_x[i];

        ts_xxyy_yyyzzz[i] = ts_yy_yyyzzz[i] * fe_0 + ts_xyy_yyyzzz[i] * pa_x[i];

        ts_xxyy_yyzzzz[i] = ts_yy_yyzzzz[i] * fe_0 + ts_xyy_yyzzzz[i] * pa_x[i];

        ts_xxyy_yzzzzz[i] = ts_yy_yzzzzz[i] * fe_0 + ts_xyy_yzzzzz[i] * pa_x[i];

        ts_xxyy_zzzzzz[i] = ts_yy_zzzzzz[i] * fe_0 + ts_xyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : GI

    auto ts_xxyz_xxxxxx = pbuffer.data(idx_ovl_gi + 112);

    auto ts_xxyz_xxxxxy = pbuffer.data(idx_ovl_gi + 113);

    auto ts_xxyz_xxxxxz = pbuffer.data(idx_ovl_gi + 114);

    auto ts_xxyz_xxxxyy = pbuffer.data(idx_ovl_gi + 115);

    auto ts_xxyz_xxxxyz = pbuffer.data(idx_ovl_gi + 116);

    auto ts_xxyz_xxxxzz = pbuffer.data(idx_ovl_gi + 117);

    auto ts_xxyz_xxxyyy = pbuffer.data(idx_ovl_gi + 118);

    auto ts_xxyz_xxxyyz = pbuffer.data(idx_ovl_gi + 119);

    auto ts_xxyz_xxxyzz = pbuffer.data(idx_ovl_gi + 120);

    auto ts_xxyz_xxxzzz = pbuffer.data(idx_ovl_gi + 121);

    auto ts_xxyz_xxyyyy = pbuffer.data(idx_ovl_gi + 122);

    auto ts_xxyz_xxyyyz = pbuffer.data(idx_ovl_gi + 123);

    auto ts_xxyz_xxyyzz = pbuffer.data(idx_ovl_gi + 124);

    auto ts_xxyz_xxyzzz = pbuffer.data(idx_ovl_gi + 125);

    auto ts_xxyz_xxzzzz = pbuffer.data(idx_ovl_gi + 126);

    auto ts_xxyz_xyyyyy = pbuffer.data(idx_ovl_gi + 127);

    auto ts_xxyz_xyyyyz = pbuffer.data(idx_ovl_gi + 128);

    auto ts_xxyz_xyyyzz = pbuffer.data(idx_ovl_gi + 129);

    auto ts_xxyz_xyyzzz = pbuffer.data(idx_ovl_gi + 130);

    auto ts_xxyz_xyzzzz = pbuffer.data(idx_ovl_gi + 131);

    auto ts_xxyz_xzzzzz = pbuffer.data(idx_ovl_gi + 132);

    auto ts_xxyz_yyyyyy = pbuffer.data(idx_ovl_gi + 133);

    auto ts_xxyz_yyyyyz = pbuffer.data(idx_ovl_gi + 134);

    auto ts_xxyz_yyyyzz = pbuffer.data(idx_ovl_gi + 135);

    auto ts_xxyz_yyyzzz = pbuffer.data(idx_ovl_gi + 136);

    auto ts_xxyz_yyzzzz = pbuffer.data(idx_ovl_gi + 137);

    auto ts_xxyz_yzzzzz = pbuffer.data(idx_ovl_gi + 138);

    auto ts_xxyz_zzzzzz = pbuffer.data(idx_ovl_gi + 139);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxy_xxxxxy, ts_xxy_xxxxyy, ts_xxy_xxxyyy, ts_xxy_xxyyyy, ts_xxy_xyyyyy, ts_xxy_yyyyyy, ts_xxyz_xxxxxx, ts_xxyz_xxxxxy, ts_xxyz_xxxxxz, ts_xxyz_xxxxyy, ts_xxyz_xxxxyz, ts_xxyz_xxxxzz, ts_xxyz_xxxyyy, ts_xxyz_xxxyyz, ts_xxyz_xxxyzz, ts_xxyz_xxxzzz, ts_xxyz_xxyyyy, ts_xxyz_xxyyyz, ts_xxyz_xxyyzz, ts_xxyz_xxyzzz, ts_xxyz_xxzzzz, ts_xxyz_xyyyyy, ts_xxyz_xyyyyz, ts_xxyz_xyyyzz, ts_xxyz_xyyzzz, ts_xxyz_xyzzzz, ts_xxyz_xzzzzz, ts_xxyz_yyyyyy, ts_xxyz_yyyyyz, ts_xxyz_yyyyzz, ts_xxyz_yyyzzz, ts_xxyz_yyzzzz, ts_xxyz_yzzzzz, ts_xxyz_zzzzzz, ts_xxz_xxxxxx, ts_xxz_xxxxxz, ts_xxz_xxxxyz, ts_xxz_xxxxz, ts_xxz_xxxxzz, ts_xxz_xxxyyz, ts_xxz_xxxyz, ts_xxz_xxxyzz, ts_xxz_xxxzz, ts_xxz_xxxzzz, ts_xxz_xxyyyz, ts_xxz_xxyyz, ts_xxz_xxyyzz, ts_xxz_xxyzz, ts_xxz_xxyzzz, ts_xxz_xxzzz, ts_xxz_xxzzzz, ts_xxz_xyyyyz, ts_xxz_xyyyz, ts_xxz_xyyyzz, ts_xxz_xyyzz, ts_xxz_xyyzzz, ts_xxz_xyzzz, ts_xxz_xyzzzz, ts_xxz_xzzzz, ts_xxz_xzzzzz, ts_xxz_zzzzzz, ts_xyz_yyyyyz, ts_xyz_yyyyzz, ts_xyz_yyyzzz, ts_xyz_yyzzzz, ts_xyz_yzzzzz, ts_yz_yyyyyz, ts_yz_yyyyzz, ts_yz_yyyzzz, ts_yz_yyzzzz, ts_yz_yzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyz_xxxxxx[i] = ts_xxz_xxxxxx[i] * pa_y[i];

        ts_xxyz_xxxxxy[i] = ts_xxy_xxxxxy[i] * pa_z[i];

        ts_xxyz_xxxxxz[i] = ts_xxz_xxxxxz[i] * pa_y[i];

        ts_xxyz_xxxxyy[i] = ts_xxy_xxxxyy[i] * pa_z[i];

        ts_xxyz_xxxxyz[i] = ts_xxz_xxxxz[i] * fe_0 + ts_xxz_xxxxyz[i] * pa_y[i];

        ts_xxyz_xxxxzz[i] = ts_xxz_xxxxzz[i] * pa_y[i];

        ts_xxyz_xxxyyy[i] = ts_xxy_xxxyyy[i] * pa_z[i];

        ts_xxyz_xxxyyz[i] = 2.0 * ts_xxz_xxxyz[i] * fe_0 + ts_xxz_xxxyyz[i] * pa_y[i];

        ts_xxyz_xxxyzz[i] = ts_xxz_xxxzz[i] * fe_0 + ts_xxz_xxxyzz[i] * pa_y[i];

        ts_xxyz_xxxzzz[i] = ts_xxz_xxxzzz[i] * pa_y[i];

        ts_xxyz_xxyyyy[i] = ts_xxy_xxyyyy[i] * pa_z[i];

        ts_xxyz_xxyyyz[i] = 3.0 * ts_xxz_xxyyz[i] * fe_0 + ts_xxz_xxyyyz[i] * pa_y[i];

        ts_xxyz_xxyyzz[i] = 2.0 * ts_xxz_xxyzz[i] * fe_0 + ts_xxz_xxyyzz[i] * pa_y[i];

        ts_xxyz_xxyzzz[i] = ts_xxz_xxzzz[i] * fe_0 + ts_xxz_xxyzzz[i] * pa_y[i];

        ts_xxyz_xxzzzz[i] = ts_xxz_xxzzzz[i] * pa_y[i];

        ts_xxyz_xyyyyy[i] = ts_xxy_xyyyyy[i] * pa_z[i];

        ts_xxyz_xyyyyz[i] = 4.0 * ts_xxz_xyyyz[i] * fe_0 + ts_xxz_xyyyyz[i] * pa_y[i];

        ts_xxyz_xyyyzz[i] = 3.0 * ts_xxz_xyyzz[i] * fe_0 + ts_xxz_xyyyzz[i] * pa_y[i];

        ts_xxyz_xyyzzz[i] = 2.0 * ts_xxz_xyzzz[i] * fe_0 + ts_xxz_xyyzzz[i] * pa_y[i];

        ts_xxyz_xyzzzz[i] = ts_xxz_xzzzz[i] * fe_0 + ts_xxz_xyzzzz[i] * pa_y[i];

        ts_xxyz_xzzzzz[i] = ts_xxz_xzzzzz[i] * pa_y[i];

        ts_xxyz_yyyyyy[i] = ts_xxy_yyyyyy[i] * pa_z[i];

        ts_xxyz_yyyyyz[i] = ts_yz_yyyyyz[i] * fe_0 + ts_xyz_yyyyyz[i] * pa_x[i];

        ts_xxyz_yyyyzz[i] = ts_yz_yyyyzz[i] * fe_0 + ts_xyz_yyyyzz[i] * pa_x[i];

        ts_xxyz_yyyzzz[i] = ts_yz_yyyzzz[i] * fe_0 + ts_xyz_yyyzzz[i] * pa_x[i];

        ts_xxyz_yyzzzz[i] = ts_yz_yyzzzz[i] * fe_0 + ts_xyz_yyzzzz[i] * pa_x[i];

        ts_xxyz_yzzzzz[i] = ts_yz_yzzzzz[i] * fe_0 + ts_xyz_yzzzzz[i] * pa_x[i];

        ts_xxyz_zzzzzz[i] = ts_xxz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : GI

    auto ts_xxzz_xxxxxx = pbuffer.data(idx_ovl_gi + 140);

    auto ts_xxzz_xxxxxy = pbuffer.data(idx_ovl_gi + 141);

    auto ts_xxzz_xxxxxz = pbuffer.data(idx_ovl_gi + 142);

    auto ts_xxzz_xxxxyy = pbuffer.data(idx_ovl_gi + 143);

    auto ts_xxzz_xxxxyz = pbuffer.data(idx_ovl_gi + 144);

    auto ts_xxzz_xxxxzz = pbuffer.data(idx_ovl_gi + 145);

    auto ts_xxzz_xxxyyy = pbuffer.data(idx_ovl_gi + 146);

    auto ts_xxzz_xxxyyz = pbuffer.data(idx_ovl_gi + 147);

    auto ts_xxzz_xxxyzz = pbuffer.data(idx_ovl_gi + 148);

    auto ts_xxzz_xxxzzz = pbuffer.data(idx_ovl_gi + 149);

    auto ts_xxzz_xxyyyy = pbuffer.data(idx_ovl_gi + 150);

    auto ts_xxzz_xxyyyz = pbuffer.data(idx_ovl_gi + 151);

    auto ts_xxzz_xxyyzz = pbuffer.data(idx_ovl_gi + 152);

    auto ts_xxzz_xxyzzz = pbuffer.data(idx_ovl_gi + 153);

    auto ts_xxzz_xxzzzz = pbuffer.data(idx_ovl_gi + 154);

    auto ts_xxzz_xyyyyy = pbuffer.data(idx_ovl_gi + 155);

    auto ts_xxzz_xyyyyz = pbuffer.data(idx_ovl_gi + 156);

    auto ts_xxzz_xyyyzz = pbuffer.data(idx_ovl_gi + 157);

    auto ts_xxzz_xyyzzz = pbuffer.data(idx_ovl_gi + 158);

    auto ts_xxzz_xyzzzz = pbuffer.data(idx_ovl_gi + 159);

    auto ts_xxzz_xzzzzz = pbuffer.data(idx_ovl_gi + 160);

    auto ts_xxzz_yyyyyy = pbuffer.data(idx_ovl_gi + 161);

    auto ts_xxzz_yyyyyz = pbuffer.data(idx_ovl_gi + 162);

    auto ts_xxzz_yyyyzz = pbuffer.data(idx_ovl_gi + 163);

    auto ts_xxzz_yyyzzz = pbuffer.data(idx_ovl_gi + 164);

    auto ts_xxzz_yyzzzz = pbuffer.data(idx_ovl_gi + 165);

    auto ts_xxzz_yzzzzz = pbuffer.data(idx_ovl_gi + 166);

    auto ts_xxzz_zzzzzz = pbuffer.data(idx_ovl_gi + 167);

    #pragma omp simd aligned(pa_x, pa_z, ts_xx_xxxxxx, ts_xx_xxxxxy, ts_xx_xxxxyy, ts_xx_xxxyyy, ts_xx_xxyyyy, ts_xx_xyyyyy, ts_xxz_xxxxxx, ts_xxz_xxxxxy, ts_xxz_xxxxyy, ts_xxz_xxxyyy, ts_xxz_xxyyyy, ts_xxz_xyyyyy, ts_xxzz_xxxxxx, ts_xxzz_xxxxxy, ts_xxzz_xxxxxz, ts_xxzz_xxxxyy, ts_xxzz_xxxxyz, ts_xxzz_xxxxzz, ts_xxzz_xxxyyy, ts_xxzz_xxxyyz, ts_xxzz_xxxyzz, ts_xxzz_xxxzzz, ts_xxzz_xxyyyy, ts_xxzz_xxyyyz, ts_xxzz_xxyyzz, ts_xxzz_xxyzzz, ts_xxzz_xxzzzz, ts_xxzz_xyyyyy, ts_xxzz_xyyyyz, ts_xxzz_xyyyzz, ts_xxzz_xyyzzz, ts_xxzz_xyzzzz, ts_xxzz_xzzzzz, ts_xxzz_yyyyyy, ts_xxzz_yyyyyz, ts_xxzz_yyyyzz, ts_xxzz_yyyzzz, ts_xxzz_yyzzzz, ts_xxzz_yzzzzz, ts_xxzz_zzzzzz, ts_xzz_xxxxxz, ts_xzz_xxxxyz, ts_xzz_xxxxz, ts_xzz_xxxxzz, ts_xzz_xxxyyz, ts_xzz_xxxyz, ts_xzz_xxxyzz, ts_xzz_xxxzz, ts_xzz_xxxzzz, ts_xzz_xxyyyz, ts_xzz_xxyyz, ts_xzz_xxyyzz, ts_xzz_xxyzz, ts_xzz_xxyzzz, ts_xzz_xxzzz, ts_xzz_xxzzzz, ts_xzz_xyyyyz, ts_xzz_xyyyz, ts_xzz_xyyyzz, ts_xzz_xyyzz, ts_xzz_xyyzzz, ts_xzz_xyzzz, ts_xzz_xyzzzz, ts_xzz_xzzzz, ts_xzz_xzzzzz, ts_xzz_yyyyyy, ts_xzz_yyyyyz, ts_xzz_yyyyz, ts_xzz_yyyyzz, ts_xzz_yyyzz, ts_xzz_yyyzzz, ts_xzz_yyzzz, ts_xzz_yyzzzz, ts_xzz_yzzzz, ts_xzz_yzzzzz, ts_xzz_zzzzz, ts_xzz_zzzzzz, ts_zz_xxxxxz, ts_zz_xxxxyz, ts_zz_xxxxzz, ts_zz_xxxyyz, ts_zz_xxxyzz, ts_zz_xxxzzz, ts_zz_xxyyyz, ts_zz_xxyyzz, ts_zz_xxyzzz, ts_zz_xxzzzz, ts_zz_xyyyyz, ts_zz_xyyyzz, ts_zz_xyyzzz, ts_zz_xyzzzz, ts_zz_xzzzzz, ts_zz_yyyyyy, ts_zz_yyyyyz, ts_zz_yyyyzz, ts_zz_yyyzzz, ts_zz_yyzzzz, ts_zz_yzzzzz, ts_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzz_xxxxxx[i] = ts_xx_xxxxxx[i] * fe_0 + ts_xxz_xxxxxx[i] * pa_z[i];

        ts_xxzz_xxxxxy[i] = ts_xx_xxxxxy[i] * fe_0 + ts_xxz_xxxxxy[i] * pa_z[i];

        ts_xxzz_xxxxxz[i] = ts_zz_xxxxxz[i] * fe_0 + 5.0 * ts_xzz_xxxxz[i] * fe_0 + ts_xzz_xxxxxz[i] * pa_x[i];

        ts_xxzz_xxxxyy[i] = ts_xx_xxxxyy[i] * fe_0 + ts_xxz_xxxxyy[i] * pa_z[i];

        ts_xxzz_xxxxyz[i] = ts_zz_xxxxyz[i] * fe_0 + 4.0 * ts_xzz_xxxyz[i] * fe_0 + ts_xzz_xxxxyz[i] * pa_x[i];

        ts_xxzz_xxxxzz[i] = ts_zz_xxxxzz[i] * fe_0 + 4.0 * ts_xzz_xxxzz[i] * fe_0 + ts_xzz_xxxxzz[i] * pa_x[i];

        ts_xxzz_xxxyyy[i] = ts_xx_xxxyyy[i] * fe_0 + ts_xxz_xxxyyy[i] * pa_z[i];

        ts_xxzz_xxxyyz[i] = ts_zz_xxxyyz[i] * fe_0 + 3.0 * ts_xzz_xxyyz[i] * fe_0 + ts_xzz_xxxyyz[i] * pa_x[i];

        ts_xxzz_xxxyzz[i] = ts_zz_xxxyzz[i] * fe_0 + 3.0 * ts_xzz_xxyzz[i] * fe_0 + ts_xzz_xxxyzz[i] * pa_x[i];

        ts_xxzz_xxxzzz[i] = ts_zz_xxxzzz[i] * fe_0 + 3.0 * ts_xzz_xxzzz[i] * fe_0 + ts_xzz_xxxzzz[i] * pa_x[i];

        ts_xxzz_xxyyyy[i] = ts_xx_xxyyyy[i] * fe_0 + ts_xxz_xxyyyy[i] * pa_z[i];

        ts_xxzz_xxyyyz[i] = ts_zz_xxyyyz[i] * fe_0 + 2.0 * ts_xzz_xyyyz[i] * fe_0 + ts_xzz_xxyyyz[i] * pa_x[i];

        ts_xxzz_xxyyzz[i] = ts_zz_xxyyzz[i] * fe_0 + 2.0 * ts_xzz_xyyzz[i] * fe_0 + ts_xzz_xxyyzz[i] * pa_x[i];

        ts_xxzz_xxyzzz[i] = ts_zz_xxyzzz[i] * fe_0 + 2.0 * ts_xzz_xyzzz[i] * fe_0 + ts_xzz_xxyzzz[i] * pa_x[i];

        ts_xxzz_xxzzzz[i] = ts_zz_xxzzzz[i] * fe_0 + 2.0 * ts_xzz_xzzzz[i] * fe_0 + ts_xzz_xxzzzz[i] * pa_x[i];

        ts_xxzz_xyyyyy[i] = ts_xx_xyyyyy[i] * fe_0 + ts_xxz_xyyyyy[i] * pa_z[i];

        ts_xxzz_xyyyyz[i] = ts_zz_xyyyyz[i] * fe_0 + ts_xzz_yyyyz[i] * fe_0 + ts_xzz_xyyyyz[i] * pa_x[i];

        ts_xxzz_xyyyzz[i] = ts_zz_xyyyzz[i] * fe_0 + ts_xzz_yyyzz[i] * fe_0 + ts_xzz_xyyyzz[i] * pa_x[i];

        ts_xxzz_xyyzzz[i] = ts_zz_xyyzzz[i] * fe_0 + ts_xzz_yyzzz[i] * fe_0 + ts_xzz_xyyzzz[i] * pa_x[i];

        ts_xxzz_xyzzzz[i] = ts_zz_xyzzzz[i] * fe_0 + ts_xzz_yzzzz[i] * fe_0 + ts_xzz_xyzzzz[i] * pa_x[i];

        ts_xxzz_xzzzzz[i] = ts_zz_xzzzzz[i] * fe_0 + ts_xzz_zzzzz[i] * fe_0 + ts_xzz_xzzzzz[i] * pa_x[i];

        ts_xxzz_yyyyyy[i] = ts_zz_yyyyyy[i] * fe_0 + ts_xzz_yyyyyy[i] * pa_x[i];

        ts_xxzz_yyyyyz[i] = ts_zz_yyyyyz[i] * fe_0 + ts_xzz_yyyyyz[i] * pa_x[i];

        ts_xxzz_yyyyzz[i] = ts_zz_yyyyzz[i] * fe_0 + ts_xzz_yyyyzz[i] * pa_x[i];

        ts_xxzz_yyyzzz[i] = ts_zz_yyyzzz[i] * fe_0 + ts_xzz_yyyzzz[i] * pa_x[i];

        ts_xxzz_yyzzzz[i] = ts_zz_yyzzzz[i] * fe_0 + ts_xzz_yyzzzz[i] * pa_x[i];

        ts_xxzz_yzzzzz[i] = ts_zz_yzzzzz[i] * fe_0 + ts_xzz_yzzzzz[i] * pa_x[i];

        ts_xxzz_zzzzzz[i] = ts_zz_zzzzzz[i] * fe_0 + ts_xzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : GI

    auto ts_xyyy_xxxxxx = pbuffer.data(idx_ovl_gi + 168);

    auto ts_xyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 169);

    auto ts_xyyy_xxxxxz = pbuffer.data(idx_ovl_gi + 170);

    auto ts_xyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 171);

    auto ts_xyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 172);

    auto ts_xyyy_xxxxzz = pbuffer.data(idx_ovl_gi + 173);

    auto ts_xyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 174);

    auto ts_xyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 175);

    auto ts_xyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 176);

    auto ts_xyyy_xxxzzz = pbuffer.data(idx_ovl_gi + 177);

    auto ts_xyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 178);

    auto ts_xyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 179);

    auto ts_xyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 180);

    auto ts_xyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 181);

    auto ts_xyyy_xxzzzz = pbuffer.data(idx_ovl_gi + 182);

    auto ts_xyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 183);

    auto ts_xyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 184);

    auto ts_xyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 185);

    auto ts_xyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 186);

    auto ts_xyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 187);

    auto ts_xyyy_xzzzzz = pbuffer.data(idx_ovl_gi + 188);

    auto ts_xyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 189);

    auto ts_xyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 190);

    auto ts_xyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 191);

    auto ts_xyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 192);

    auto ts_xyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 193);

    auto ts_xyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 194);

    auto ts_xyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 195);

    #pragma omp simd aligned(pa_x, ts_xyyy_xxxxxx, ts_xyyy_xxxxxy, ts_xyyy_xxxxxz, ts_xyyy_xxxxyy, ts_xyyy_xxxxyz, ts_xyyy_xxxxzz, ts_xyyy_xxxyyy, ts_xyyy_xxxyyz, ts_xyyy_xxxyzz, ts_xyyy_xxxzzz, ts_xyyy_xxyyyy, ts_xyyy_xxyyyz, ts_xyyy_xxyyzz, ts_xyyy_xxyzzz, ts_xyyy_xxzzzz, ts_xyyy_xyyyyy, ts_xyyy_xyyyyz, ts_xyyy_xyyyzz, ts_xyyy_xyyzzz, ts_xyyy_xyzzzz, ts_xyyy_xzzzzz, ts_xyyy_yyyyyy, ts_xyyy_yyyyyz, ts_xyyy_yyyyzz, ts_xyyy_yyyzzz, ts_xyyy_yyzzzz, ts_xyyy_yzzzzz, ts_xyyy_zzzzzz, ts_yyy_xxxxx, ts_yyy_xxxxxx, ts_yyy_xxxxxy, ts_yyy_xxxxxz, ts_yyy_xxxxy, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxxz, ts_yyy_xxxxzz, ts_yyy_xxxyy, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyz, ts_yyy_xxxyzz, ts_yyy_xxxzz, ts_yyy_xxxzzz, ts_yyy_xxyyy, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyz, ts_yyy_xxyyzz, ts_yyy_xxyzz, ts_yyy_xxyzzz, ts_yyy_xxzzz, ts_yyy_xxzzzz, ts_yyy_xyyyy, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzz, ts_yyy_xyzzzz, ts_yyy_xzzzz, ts_yyy_xzzzzz, ts_yyy_yyyyy, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzz, ts_yyy_yzzzzz, ts_yyy_zzzzz, ts_yyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyy_xxxxxx[i] = 6.0 * ts_yyy_xxxxx[i] * fe_0 + ts_yyy_xxxxxx[i] * pa_x[i];

        ts_xyyy_xxxxxy[i] = 5.0 * ts_yyy_xxxxy[i] * fe_0 + ts_yyy_xxxxxy[i] * pa_x[i];

        ts_xyyy_xxxxxz[i] = 5.0 * ts_yyy_xxxxz[i] * fe_0 + ts_yyy_xxxxxz[i] * pa_x[i];

        ts_xyyy_xxxxyy[i] = 4.0 * ts_yyy_xxxyy[i] * fe_0 + ts_yyy_xxxxyy[i] * pa_x[i];

        ts_xyyy_xxxxyz[i] = 4.0 * ts_yyy_xxxyz[i] * fe_0 + ts_yyy_xxxxyz[i] * pa_x[i];

        ts_xyyy_xxxxzz[i] = 4.0 * ts_yyy_xxxzz[i] * fe_0 + ts_yyy_xxxxzz[i] * pa_x[i];

        ts_xyyy_xxxyyy[i] = 3.0 * ts_yyy_xxyyy[i] * fe_0 + ts_yyy_xxxyyy[i] * pa_x[i];

        ts_xyyy_xxxyyz[i] = 3.0 * ts_yyy_xxyyz[i] * fe_0 + ts_yyy_xxxyyz[i] * pa_x[i];

        ts_xyyy_xxxyzz[i] = 3.0 * ts_yyy_xxyzz[i] * fe_0 + ts_yyy_xxxyzz[i] * pa_x[i];

        ts_xyyy_xxxzzz[i] = 3.0 * ts_yyy_xxzzz[i] * fe_0 + ts_yyy_xxxzzz[i] * pa_x[i];

        ts_xyyy_xxyyyy[i] = 2.0 * ts_yyy_xyyyy[i] * fe_0 + ts_yyy_xxyyyy[i] * pa_x[i];

        ts_xyyy_xxyyyz[i] = 2.0 * ts_yyy_xyyyz[i] * fe_0 + ts_yyy_xxyyyz[i] * pa_x[i];

        ts_xyyy_xxyyzz[i] = 2.0 * ts_yyy_xyyzz[i] * fe_0 + ts_yyy_xxyyzz[i] * pa_x[i];

        ts_xyyy_xxyzzz[i] = 2.0 * ts_yyy_xyzzz[i] * fe_0 + ts_yyy_xxyzzz[i] * pa_x[i];

        ts_xyyy_xxzzzz[i] = 2.0 * ts_yyy_xzzzz[i] * fe_0 + ts_yyy_xxzzzz[i] * pa_x[i];

        ts_xyyy_xyyyyy[i] = ts_yyy_yyyyy[i] * fe_0 + ts_yyy_xyyyyy[i] * pa_x[i];

        ts_xyyy_xyyyyz[i] = ts_yyy_yyyyz[i] * fe_0 + ts_yyy_xyyyyz[i] * pa_x[i];

        ts_xyyy_xyyyzz[i] = ts_yyy_yyyzz[i] * fe_0 + ts_yyy_xyyyzz[i] * pa_x[i];

        ts_xyyy_xyyzzz[i] = ts_yyy_yyzzz[i] * fe_0 + ts_yyy_xyyzzz[i] * pa_x[i];

        ts_xyyy_xyzzzz[i] = ts_yyy_yzzzz[i] * fe_0 + ts_yyy_xyzzzz[i] * pa_x[i];

        ts_xyyy_xzzzzz[i] = ts_yyy_zzzzz[i] * fe_0 + ts_yyy_xzzzzz[i] * pa_x[i];

        ts_xyyy_yyyyyy[i] = ts_yyy_yyyyyy[i] * pa_x[i];

        ts_xyyy_yyyyyz[i] = ts_yyy_yyyyyz[i] * pa_x[i];

        ts_xyyy_yyyyzz[i] = ts_yyy_yyyyzz[i] * pa_x[i];

        ts_xyyy_yyyzzz[i] = ts_yyy_yyyzzz[i] * pa_x[i];

        ts_xyyy_yyzzzz[i] = ts_yyy_yyzzzz[i] * pa_x[i];

        ts_xyyy_yzzzzz[i] = ts_yyy_yzzzzz[i] * pa_x[i];

        ts_xyyy_zzzzzz[i] = ts_yyy_zzzzzz[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : GI

    auto ts_xyyz_xxxxxx = pbuffer.data(idx_ovl_gi + 196);

    auto ts_xyyz_xxxxxy = pbuffer.data(idx_ovl_gi + 197);

    auto ts_xyyz_xxxxxz = pbuffer.data(idx_ovl_gi + 198);

    auto ts_xyyz_xxxxyy = pbuffer.data(idx_ovl_gi + 199);

    auto ts_xyyz_xxxxyz = pbuffer.data(idx_ovl_gi + 200);

    auto ts_xyyz_xxxxzz = pbuffer.data(idx_ovl_gi + 201);

    auto ts_xyyz_xxxyyy = pbuffer.data(idx_ovl_gi + 202);

    auto ts_xyyz_xxxyyz = pbuffer.data(idx_ovl_gi + 203);

    auto ts_xyyz_xxxyzz = pbuffer.data(idx_ovl_gi + 204);

    auto ts_xyyz_xxxzzz = pbuffer.data(idx_ovl_gi + 205);

    auto ts_xyyz_xxyyyy = pbuffer.data(idx_ovl_gi + 206);

    auto ts_xyyz_xxyyyz = pbuffer.data(idx_ovl_gi + 207);

    auto ts_xyyz_xxyyzz = pbuffer.data(idx_ovl_gi + 208);

    auto ts_xyyz_xxyzzz = pbuffer.data(idx_ovl_gi + 209);

    auto ts_xyyz_xxzzzz = pbuffer.data(idx_ovl_gi + 210);

    auto ts_xyyz_xyyyyy = pbuffer.data(idx_ovl_gi + 211);

    auto ts_xyyz_xyyyyz = pbuffer.data(idx_ovl_gi + 212);

    auto ts_xyyz_xyyyzz = pbuffer.data(idx_ovl_gi + 213);

    auto ts_xyyz_xyyzzz = pbuffer.data(idx_ovl_gi + 214);

    auto ts_xyyz_xyzzzz = pbuffer.data(idx_ovl_gi + 215);

    auto ts_xyyz_xzzzzz = pbuffer.data(idx_ovl_gi + 216);

    auto ts_xyyz_yyyyyy = pbuffer.data(idx_ovl_gi + 217);

    auto ts_xyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 218);

    auto ts_xyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 219);

    auto ts_xyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 220);

    auto ts_xyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 221);

    auto ts_xyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 222);

    auto ts_xyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 223);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyy_xxxxxx, ts_xyy_xxxxxy, ts_xyy_xxxxyy, ts_xyy_xxxyyy, ts_xyy_xxyyyy, ts_xyy_xyyyyy, ts_xyyz_xxxxxx, ts_xyyz_xxxxxy, ts_xyyz_xxxxxz, ts_xyyz_xxxxyy, ts_xyyz_xxxxyz, ts_xyyz_xxxxzz, ts_xyyz_xxxyyy, ts_xyyz_xxxyyz, ts_xyyz_xxxyzz, ts_xyyz_xxxzzz, ts_xyyz_xxyyyy, ts_xyyz_xxyyyz, ts_xyyz_xxyyzz, ts_xyyz_xxyzzz, ts_xyyz_xxzzzz, ts_xyyz_xyyyyy, ts_xyyz_xyyyyz, ts_xyyz_xyyyzz, ts_xyyz_xyyzzz, ts_xyyz_xyzzzz, ts_xyyz_xzzzzz, ts_xyyz_yyyyyy, ts_xyyz_yyyyyz, ts_xyyz_yyyyzz, ts_xyyz_yyyzzz, ts_xyyz_yyzzzz, ts_xyyz_yzzzzz, ts_xyyz_zzzzzz, ts_yyz_xxxxxz, ts_yyz_xxxxyz, ts_yyz_xxxxz, ts_yyz_xxxxzz, ts_yyz_xxxyyz, ts_yyz_xxxyz, ts_yyz_xxxyzz, ts_yyz_xxxzz, ts_yyz_xxxzzz, ts_yyz_xxyyyz, ts_yyz_xxyyz, ts_yyz_xxyyzz, ts_yyz_xxyzz, ts_yyz_xxyzzz, ts_yyz_xxzzz, ts_yyz_xxzzzz, ts_yyz_xyyyyz, ts_yyz_xyyyz, ts_yyz_xyyyzz, ts_yyz_xyyzz, ts_yyz_xyyzzz, ts_yyz_xyzzz, ts_yyz_xyzzzz, ts_yyz_xzzzz, ts_yyz_xzzzzz, ts_yyz_yyyyyy, ts_yyz_yyyyyz, ts_yyz_yyyyz, ts_yyz_yyyyzz, ts_yyz_yyyzz, ts_yyz_yyyzzz, ts_yyz_yyzzz, ts_yyz_yyzzzz, ts_yyz_yzzzz, ts_yyz_yzzzzz, ts_yyz_zzzzz, ts_yyz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyz_xxxxxx[i] = ts_xyy_xxxxxx[i] * pa_z[i];

        ts_xyyz_xxxxxy[i] = ts_xyy_xxxxxy[i] * pa_z[i];

        ts_xyyz_xxxxxz[i] = 5.0 * ts_yyz_xxxxz[i] * fe_0 + ts_yyz_xxxxxz[i] * pa_x[i];

        ts_xyyz_xxxxyy[i] = ts_xyy_xxxxyy[i] * pa_z[i];

        ts_xyyz_xxxxyz[i] = 4.0 * ts_yyz_xxxyz[i] * fe_0 + ts_yyz_xxxxyz[i] * pa_x[i];

        ts_xyyz_xxxxzz[i] = 4.0 * ts_yyz_xxxzz[i] * fe_0 + ts_yyz_xxxxzz[i] * pa_x[i];

        ts_xyyz_xxxyyy[i] = ts_xyy_xxxyyy[i] * pa_z[i];

        ts_xyyz_xxxyyz[i] = 3.0 * ts_yyz_xxyyz[i] * fe_0 + ts_yyz_xxxyyz[i] * pa_x[i];

        ts_xyyz_xxxyzz[i] = 3.0 * ts_yyz_xxyzz[i] * fe_0 + ts_yyz_xxxyzz[i] * pa_x[i];

        ts_xyyz_xxxzzz[i] = 3.0 * ts_yyz_xxzzz[i] * fe_0 + ts_yyz_xxxzzz[i] * pa_x[i];

        ts_xyyz_xxyyyy[i] = ts_xyy_xxyyyy[i] * pa_z[i];

        ts_xyyz_xxyyyz[i] = 2.0 * ts_yyz_xyyyz[i] * fe_0 + ts_yyz_xxyyyz[i] * pa_x[i];

        ts_xyyz_xxyyzz[i] = 2.0 * ts_yyz_xyyzz[i] * fe_0 + ts_yyz_xxyyzz[i] * pa_x[i];

        ts_xyyz_xxyzzz[i] = 2.0 * ts_yyz_xyzzz[i] * fe_0 + ts_yyz_xxyzzz[i] * pa_x[i];

        ts_xyyz_xxzzzz[i] = 2.0 * ts_yyz_xzzzz[i] * fe_0 + ts_yyz_xxzzzz[i] * pa_x[i];

        ts_xyyz_xyyyyy[i] = ts_xyy_xyyyyy[i] * pa_z[i];

        ts_xyyz_xyyyyz[i] = ts_yyz_yyyyz[i] * fe_0 + ts_yyz_xyyyyz[i] * pa_x[i];

        ts_xyyz_xyyyzz[i] = ts_yyz_yyyzz[i] * fe_0 + ts_yyz_xyyyzz[i] * pa_x[i];

        ts_xyyz_xyyzzz[i] = ts_yyz_yyzzz[i] * fe_0 + ts_yyz_xyyzzz[i] * pa_x[i];

        ts_xyyz_xyzzzz[i] = ts_yyz_yzzzz[i] * fe_0 + ts_yyz_xyzzzz[i] * pa_x[i];

        ts_xyyz_xzzzzz[i] = ts_yyz_zzzzz[i] * fe_0 + ts_yyz_xzzzzz[i] * pa_x[i];

        ts_xyyz_yyyyyy[i] = ts_yyz_yyyyyy[i] * pa_x[i];

        ts_xyyz_yyyyyz[i] = ts_yyz_yyyyyz[i] * pa_x[i];

        ts_xyyz_yyyyzz[i] = ts_yyz_yyyyzz[i] * pa_x[i];

        ts_xyyz_yyyzzz[i] = ts_yyz_yyyzzz[i] * pa_x[i];

        ts_xyyz_yyzzzz[i] = ts_yyz_yyzzzz[i] * pa_x[i];

        ts_xyyz_yzzzzz[i] = ts_yyz_yzzzzz[i] * pa_x[i];

        ts_xyyz_zzzzzz[i] = ts_yyz_zzzzzz[i] * pa_x[i];
    }

    // Set up 224-252 components of targeted buffer : GI

    auto ts_xyzz_xxxxxx = pbuffer.data(idx_ovl_gi + 224);

    auto ts_xyzz_xxxxxy = pbuffer.data(idx_ovl_gi + 225);

    auto ts_xyzz_xxxxxz = pbuffer.data(idx_ovl_gi + 226);

    auto ts_xyzz_xxxxyy = pbuffer.data(idx_ovl_gi + 227);

    auto ts_xyzz_xxxxyz = pbuffer.data(idx_ovl_gi + 228);

    auto ts_xyzz_xxxxzz = pbuffer.data(idx_ovl_gi + 229);

    auto ts_xyzz_xxxyyy = pbuffer.data(idx_ovl_gi + 230);

    auto ts_xyzz_xxxyyz = pbuffer.data(idx_ovl_gi + 231);

    auto ts_xyzz_xxxyzz = pbuffer.data(idx_ovl_gi + 232);

    auto ts_xyzz_xxxzzz = pbuffer.data(idx_ovl_gi + 233);

    auto ts_xyzz_xxyyyy = pbuffer.data(idx_ovl_gi + 234);

    auto ts_xyzz_xxyyyz = pbuffer.data(idx_ovl_gi + 235);

    auto ts_xyzz_xxyyzz = pbuffer.data(idx_ovl_gi + 236);

    auto ts_xyzz_xxyzzz = pbuffer.data(idx_ovl_gi + 237);

    auto ts_xyzz_xxzzzz = pbuffer.data(idx_ovl_gi + 238);

    auto ts_xyzz_xyyyyy = pbuffer.data(idx_ovl_gi + 239);

    auto ts_xyzz_xyyyyz = pbuffer.data(idx_ovl_gi + 240);

    auto ts_xyzz_xyyyzz = pbuffer.data(idx_ovl_gi + 241);

    auto ts_xyzz_xyyzzz = pbuffer.data(idx_ovl_gi + 242);

    auto ts_xyzz_xyzzzz = pbuffer.data(idx_ovl_gi + 243);

    auto ts_xyzz_xzzzzz = pbuffer.data(idx_ovl_gi + 244);

    auto ts_xyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 245);

    auto ts_xyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 246);

    auto ts_xyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 247);

    auto ts_xyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 248);

    auto ts_xyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 249);

    auto ts_xyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 250);

    auto ts_xyzz_zzzzzz = pbuffer.data(idx_ovl_gi + 251);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzz_xxxxxx, ts_xyzz_xxxxxy, ts_xyzz_xxxxxz, ts_xyzz_xxxxyy, ts_xyzz_xxxxyz, ts_xyzz_xxxxzz, ts_xyzz_xxxyyy, ts_xyzz_xxxyyz, ts_xyzz_xxxyzz, ts_xyzz_xxxzzz, ts_xyzz_xxyyyy, ts_xyzz_xxyyyz, ts_xyzz_xxyyzz, ts_xyzz_xxyzzz, ts_xyzz_xxzzzz, ts_xyzz_xyyyyy, ts_xyzz_xyyyyz, ts_xyzz_xyyyzz, ts_xyzz_xyyzzz, ts_xyzz_xyzzzz, ts_xyzz_xzzzzz, ts_xyzz_yyyyyy, ts_xyzz_yyyyyz, ts_xyzz_yyyyzz, ts_xyzz_yyyzzz, ts_xyzz_yyzzzz, ts_xyzz_yzzzzz, ts_xyzz_zzzzzz, ts_xzz_xxxxxx, ts_xzz_xxxxxz, ts_xzz_xxxxzz, ts_xzz_xxxzzz, ts_xzz_xxzzzz, ts_xzz_xzzzzz, ts_yzz_xxxxxy, ts_yzz_xxxxy, ts_yzz_xxxxyy, ts_yzz_xxxxyz, ts_yzz_xxxyy, ts_yzz_xxxyyy, ts_yzz_xxxyyz, ts_yzz_xxxyz, ts_yzz_xxxyzz, ts_yzz_xxyyy, ts_yzz_xxyyyy, ts_yzz_xxyyyz, ts_yzz_xxyyz, ts_yzz_xxyyzz, ts_yzz_xxyzz, ts_yzz_xxyzzz, ts_yzz_xyyyy, ts_yzz_xyyyyy, ts_yzz_xyyyyz, ts_yzz_xyyyz, ts_yzz_xyyyzz, ts_yzz_xyyzz, ts_yzz_xyyzzz, ts_yzz_xyzzz, ts_yzz_xyzzzz, ts_yzz_yyyyy, ts_yzz_yyyyyy, ts_yzz_yyyyyz, ts_yzz_yyyyz, ts_yzz_yyyyzz, ts_yzz_yyyzz, ts_yzz_yyyzzz, ts_yzz_yyzzz, ts_yzz_yyzzzz, ts_yzz_yzzzz, ts_yzz_yzzzzz, ts_yzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzz_xxxxxx[i] = ts_xzz_xxxxxx[i] * pa_y[i];

        ts_xyzz_xxxxxy[i] = 5.0 * ts_yzz_xxxxy[i] * fe_0 + ts_yzz_xxxxxy[i] * pa_x[i];

        ts_xyzz_xxxxxz[i] = ts_xzz_xxxxxz[i] * pa_y[i];

        ts_xyzz_xxxxyy[i] = 4.0 * ts_yzz_xxxyy[i] * fe_0 + ts_yzz_xxxxyy[i] * pa_x[i];

        ts_xyzz_xxxxyz[i] = 4.0 * ts_yzz_xxxyz[i] * fe_0 + ts_yzz_xxxxyz[i] * pa_x[i];

        ts_xyzz_xxxxzz[i] = ts_xzz_xxxxzz[i] * pa_y[i];

        ts_xyzz_xxxyyy[i] = 3.0 * ts_yzz_xxyyy[i] * fe_0 + ts_yzz_xxxyyy[i] * pa_x[i];

        ts_xyzz_xxxyyz[i] = 3.0 * ts_yzz_xxyyz[i] * fe_0 + ts_yzz_xxxyyz[i] * pa_x[i];

        ts_xyzz_xxxyzz[i] = 3.0 * ts_yzz_xxyzz[i] * fe_0 + ts_yzz_xxxyzz[i] * pa_x[i];

        ts_xyzz_xxxzzz[i] = ts_xzz_xxxzzz[i] * pa_y[i];

        ts_xyzz_xxyyyy[i] = 2.0 * ts_yzz_xyyyy[i] * fe_0 + ts_yzz_xxyyyy[i] * pa_x[i];

        ts_xyzz_xxyyyz[i] = 2.0 * ts_yzz_xyyyz[i] * fe_0 + ts_yzz_xxyyyz[i] * pa_x[i];

        ts_xyzz_xxyyzz[i] = 2.0 * ts_yzz_xyyzz[i] * fe_0 + ts_yzz_xxyyzz[i] * pa_x[i];

        ts_xyzz_xxyzzz[i] = 2.0 * ts_yzz_xyzzz[i] * fe_0 + ts_yzz_xxyzzz[i] * pa_x[i];

        ts_xyzz_xxzzzz[i] = ts_xzz_xxzzzz[i] * pa_y[i];

        ts_xyzz_xyyyyy[i] = ts_yzz_yyyyy[i] * fe_0 + ts_yzz_xyyyyy[i] * pa_x[i];

        ts_xyzz_xyyyyz[i] = ts_yzz_yyyyz[i] * fe_0 + ts_yzz_xyyyyz[i] * pa_x[i];

        ts_xyzz_xyyyzz[i] = ts_yzz_yyyzz[i] * fe_0 + ts_yzz_xyyyzz[i] * pa_x[i];

        ts_xyzz_xyyzzz[i] = ts_yzz_yyzzz[i] * fe_0 + ts_yzz_xyyzzz[i] * pa_x[i];

        ts_xyzz_xyzzzz[i] = ts_yzz_yzzzz[i] * fe_0 + ts_yzz_xyzzzz[i] * pa_x[i];

        ts_xyzz_xzzzzz[i] = ts_xzz_xzzzzz[i] * pa_y[i];

        ts_xyzz_yyyyyy[i] = ts_yzz_yyyyyy[i] * pa_x[i];

        ts_xyzz_yyyyyz[i] = ts_yzz_yyyyyz[i] * pa_x[i];

        ts_xyzz_yyyyzz[i] = ts_yzz_yyyyzz[i] * pa_x[i];

        ts_xyzz_yyyzzz[i] = ts_yzz_yyyzzz[i] * pa_x[i];

        ts_xyzz_yyzzzz[i] = ts_yzz_yyzzzz[i] * pa_x[i];

        ts_xyzz_yzzzzz[i] = ts_yzz_yzzzzz[i] * pa_x[i];

        ts_xyzz_zzzzzz[i] = ts_yzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 252-280 components of targeted buffer : GI

    auto ts_xzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 252);

    auto ts_xzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 253);

    auto ts_xzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 254);

    auto ts_xzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 255);

    auto ts_xzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 256);

    auto ts_xzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 257);

    auto ts_xzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 258);

    auto ts_xzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 259);

    auto ts_xzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 260);

    auto ts_xzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 261);

    auto ts_xzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 262);

    auto ts_xzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 263);

    auto ts_xzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 264);

    auto ts_xzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 265);

    auto ts_xzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 266);

    auto ts_xzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 267);

    auto ts_xzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 268);

    auto ts_xzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 269);

    auto ts_xzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 270);

    auto ts_xzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 271);

    auto ts_xzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 272);

    auto ts_xzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 273);

    auto ts_xzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 274);

    auto ts_xzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 275);

    auto ts_xzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 276);

    auto ts_xzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 277);

    auto ts_xzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 278);

    auto ts_xzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 279);

    #pragma omp simd aligned(pa_x, ts_xzzz_xxxxxx, ts_xzzz_xxxxxy, ts_xzzz_xxxxxz, ts_xzzz_xxxxyy, ts_xzzz_xxxxyz, ts_xzzz_xxxxzz, ts_xzzz_xxxyyy, ts_xzzz_xxxyyz, ts_xzzz_xxxyzz, ts_xzzz_xxxzzz, ts_xzzz_xxyyyy, ts_xzzz_xxyyyz, ts_xzzz_xxyyzz, ts_xzzz_xxyzzz, ts_xzzz_xxzzzz, ts_xzzz_xyyyyy, ts_xzzz_xyyyyz, ts_xzzz_xyyyzz, ts_xzzz_xyyzzz, ts_xzzz_xyzzzz, ts_xzzz_xzzzzz, ts_xzzz_yyyyyy, ts_xzzz_yyyyyz, ts_xzzz_yyyyzz, ts_xzzz_yyyzzz, ts_xzzz_yyzzzz, ts_xzzz_yzzzzz, ts_xzzz_zzzzzz, ts_zzz_xxxxx, ts_zzz_xxxxxx, ts_zzz_xxxxxy, ts_zzz_xxxxxz, ts_zzz_xxxxy, ts_zzz_xxxxyy, ts_zzz_xxxxyz, ts_zzz_xxxxz, ts_zzz_xxxxzz, ts_zzz_xxxyy, ts_zzz_xxxyyy, ts_zzz_xxxyyz, ts_zzz_xxxyz, ts_zzz_xxxyzz, ts_zzz_xxxzz, ts_zzz_xxxzzz, ts_zzz_xxyyy, ts_zzz_xxyyyy, ts_zzz_xxyyyz, ts_zzz_xxyyz, ts_zzz_xxyyzz, ts_zzz_xxyzz, ts_zzz_xxyzzz, ts_zzz_xxzzz, ts_zzz_xxzzzz, ts_zzz_xyyyy, ts_zzz_xyyyyy, ts_zzz_xyyyyz, ts_zzz_xyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyy, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzz_xxxxxx[i] = 6.0 * ts_zzz_xxxxx[i] * fe_0 + ts_zzz_xxxxxx[i] * pa_x[i];

        ts_xzzz_xxxxxy[i] = 5.0 * ts_zzz_xxxxy[i] * fe_0 + ts_zzz_xxxxxy[i] * pa_x[i];

        ts_xzzz_xxxxxz[i] = 5.0 * ts_zzz_xxxxz[i] * fe_0 + ts_zzz_xxxxxz[i] * pa_x[i];

        ts_xzzz_xxxxyy[i] = 4.0 * ts_zzz_xxxyy[i] * fe_0 + ts_zzz_xxxxyy[i] * pa_x[i];

        ts_xzzz_xxxxyz[i] = 4.0 * ts_zzz_xxxyz[i] * fe_0 + ts_zzz_xxxxyz[i] * pa_x[i];

        ts_xzzz_xxxxzz[i] = 4.0 * ts_zzz_xxxzz[i] * fe_0 + ts_zzz_xxxxzz[i] * pa_x[i];

        ts_xzzz_xxxyyy[i] = 3.0 * ts_zzz_xxyyy[i] * fe_0 + ts_zzz_xxxyyy[i] * pa_x[i];

        ts_xzzz_xxxyyz[i] = 3.0 * ts_zzz_xxyyz[i] * fe_0 + ts_zzz_xxxyyz[i] * pa_x[i];

        ts_xzzz_xxxyzz[i] = 3.0 * ts_zzz_xxyzz[i] * fe_0 + ts_zzz_xxxyzz[i] * pa_x[i];

        ts_xzzz_xxxzzz[i] = 3.0 * ts_zzz_xxzzz[i] * fe_0 + ts_zzz_xxxzzz[i] * pa_x[i];

        ts_xzzz_xxyyyy[i] = 2.0 * ts_zzz_xyyyy[i] * fe_0 + ts_zzz_xxyyyy[i] * pa_x[i];

        ts_xzzz_xxyyyz[i] = 2.0 * ts_zzz_xyyyz[i] * fe_0 + ts_zzz_xxyyyz[i] * pa_x[i];

        ts_xzzz_xxyyzz[i] = 2.0 * ts_zzz_xyyzz[i] * fe_0 + ts_zzz_xxyyzz[i] * pa_x[i];

        ts_xzzz_xxyzzz[i] = 2.0 * ts_zzz_xyzzz[i] * fe_0 + ts_zzz_xxyzzz[i] * pa_x[i];

        ts_xzzz_xxzzzz[i] = 2.0 * ts_zzz_xzzzz[i] * fe_0 + ts_zzz_xxzzzz[i] * pa_x[i];

        ts_xzzz_xyyyyy[i] = ts_zzz_yyyyy[i] * fe_0 + ts_zzz_xyyyyy[i] * pa_x[i];

        ts_xzzz_xyyyyz[i] = ts_zzz_yyyyz[i] * fe_0 + ts_zzz_xyyyyz[i] * pa_x[i];

        ts_xzzz_xyyyzz[i] = ts_zzz_yyyzz[i] * fe_0 + ts_zzz_xyyyzz[i] * pa_x[i];

        ts_xzzz_xyyzzz[i] = ts_zzz_yyzzz[i] * fe_0 + ts_zzz_xyyzzz[i] * pa_x[i];

        ts_xzzz_xyzzzz[i] = ts_zzz_yzzzz[i] * fe_0 + ts_zzz_xyzzzz[i] * pa_x[i];

        ts_xzzz_xzzzzz[i] = ts_zzz_zzzzz[i] * fe_0 + ts_zzz_xzzzzz[i] * pa_x[i];

        ts_xzzz_yyyyyy[i] = ts_zzz_yyyyyy[i] * pa_x[i];

        ts_xzzz_yyyyyz[i] = ts_zzz_yyyyyz[i] * pa_x[i];

        ts_xzzz_yyyyzz[i] = ts_zzz_yyyyzz[i] * pa_x[i];

        ts_xzzz_yyyzzz[i] = ts_zzz_yyyzzz[i] * pa_x[i];

        ts_xzzz_yyzzzz[i] = ts_zzz_yyzzzz[i] * pa_x[i];

        ts_xzzz_yzzzzz[i] = ts_zzz_yzzzzz[i] * pa_x[i];

        ts_xzzz_zzzzzz[i] = ts_zzz_zzzzzz[i] * pa_x[i];
    }

    // Set up 280-308 components of targeted buffer : GI

    auto ts_yyyy_xxxxxx = pbuffer.data(idx_ovl_gi + 280);

    auto ts_yyyy_xxxxxy = pbuffer.data(idx_ovl_gi + 281);

    auto ts_yyyy_xxxxxz = pbuffer.data(idx_ovl_gi + 282);

    auto ts_yyyy_xxxxyy = pbuffer.data(idx_ovl_gi + 283);

    auto ts_yyyy_xxxxyz = pbuffer.data(idx_ovl_gi + 284);

    auto ts_yyyy_xxxxzz = pbuffer.data(idx_ovl_gi + 285);

    auto ts_yyyy_xxxyyy = pbuffer.data(idx_ovl_gi + 286);

    auto ts_yyyy_xxxyyz = pbuffer.data(idx_ovl_gi + 287);

    auto ts_yyyy_xxxyzz = pbuffer.data(idx_ovl_gi + 288);

    auto ts_yyyy_xxxzzz = pbuffer.data(idx_ovl_gi + 289);

    auto ts_yyyy_xxyyyy = pbuffer.data(idx_ovl_gi + 290);

    auto ts_yyyy_xxyyyz = pbuffer.data(idx_ovl_gi + 291);

    auto ts_yyyy_xxyyzz = pbuffer.data(idx_ovl_gi + 292);

    auto ts_yyyy_xxyzzz = pbuffer.data(idx_ovl_gi + 293);

    auto ts_yyyy_xxzzzz = pbuffer.data(idx_ovl_gi + 294);

    auto ts_yyyy_xyyyyy = pbuffer.data(idx_ovl_gi + 295);

    auto ts_yyyy_xyyyyz = pbuffer.data(idx_ovl_gi + 296);

    auto ts_yyyy_xyyyzz = pbuffer.data(idx_ovl_gi + 297);

    auto ts_yyyy_xyyzzz = pbuffer.data(idx_ovl_gi + 298);

    auto ts_yyyy_xyzzzz = pbuffer.data(idx_ovl_gi + 299);

    auto ts_yyyy_xzzzzz = pbuffer.data(idx_ovl_gi + 300);

    auto ts_yyyy_yyyyyy = pbuffer.data(idx_ovl_gi + 301);

    auto ts_yyyy_yyyyyz = pbuffer.data(idx_ovl_gi + 302);

    auto ts_yyyy_yyyyzz = pbuffer.data(idx_ovl_gi + 303);

    auto ts_yyyy_yyyzzz = pbuffer.data(idx_ovl_gi + 304);

    auto ts_yyyy_yyzzzz = pbuffer.data(idx_ovl_gi + 305);

    auto ts_yyyy_yzzzzz = pbuffer.data(idx_ovl_gi + 306);

    auto ts_yyyy_zzzzzz = pbuffer.data(idx_ovl_gi + 307);

    #pragma omp simd aligned(pa_y, ts_yy_xxxxxx, ts_yy_xxxxxy, ts_yy_xxxxxz, ts_yy_xxxxyy, ts_yy_xxxxyz, ts_yy_xxxxzz, ts_yy_xxxyyy, ts_yy_xxxyyz, ts_yy_xxxyzz, ts_yy_xxxzzz, ts_yy_xxyyyy, ts_yy_xxyyyz, ts_yy_xxyyzz, ts_yy_xxyzzz, ts_yy_xxzzzz, ts_yy_xyyyyy, ts_yy_xyyyyz, ts_yy_xyyyzz, ts_yy_xyyzzz, ts_yy_xyzzzz, ts_yy_xzzzzz, ts_yy_yyyyyy, ts_yy_yyyyyz, ts_yy_yyyyzz, ts_yy_yyyzzz, ts_yy_yyzzzz, ts_yy_yzzzzz, ts_yy_zzzzzz, ts_yyy_xxxxx, ts_yyy_xxxxxx, ts_yyy_xxxxxy, ts_yyy_xxxxxz, ts_yyy_xxxxy, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxxz, ts_yyy_xxxxzz, ts_yyy_xxxyy, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyz, ts_yyy_xxxyzz, ts_yyy_xxxzz, ts_yyy_xxxzzz, ts_yyy_xxyyy, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyz, ts_yyy_xxyyzz, ts_yyy_xxyzz, ts_yyy_xxyzzz, ts_yyy_xxzzz, ts_yyy_xxzzzz, ts_yyy_xyyyy, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzz, ts_yyy_xyzzzz, ts_yyy_xzzzz, ts_yyy_xzzzzz, ts_yyy_yyyyy, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzz, ts_yyy_yzzzzz, ts_yyy_zzzzz, ts_yyy_zzzzzz, ts_yyyy_xxxxxx, ts_yyyy_xxxxxy, ts_yyyy_xxxxxz, ts_yyyy_xxxxyy, ts_yyyy_xxxxyz, ts_yyyy_xxxxzz, ts_yyyy_xxxyyy, ts_yyyy_xxxyyz, ts_yyyy_xxxyzz, ts_yyyy_xxxzzz, ts_yyyy_xxyyyy, ts_yyyy_xxyyyz, ts_yyyy_xxyyzz, ts_yyyy_xxyzzz, ts_yyyy_xxzzzz, ts_yyyy_xyyyyy, ts_yyyy_xyyyyz, ts_yyyy_xyyyzz, ts_yyyy_xyyzzz, ts_yyyy_xyzzzz, ts_yyyy_xzzzzz, ts_yyyy_yyyyyy, ts_yyyy_yyyyyz, ts_yyyy_yyyyzz, ts_yyyy_yyyzzz, ts_yyyy_yyzzzz, ts_yyyy_yzzzzz, ts_yyyy_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyy_xxxxxx[i] = 3.0 * ts_yy_xxxxxx[i] * fe_0 + ts_yyy_xxxxxx[i] * pa_y[i];

        ts_yyyy_xxxxxy[i] = 3.0 * ts_yy_xxxxxy[i] * fe_0 + ts_yyy_xxxxx[i] * fe_0 + ts_yyy_xxxxxy[i] * pa_y[i];

        ts_yyyy_xxxxxz[i] = 3.0 * ts_yy_xxxxxz[i] * fe_0 + ts_yyy_xxxxxz[i] * pa_y[i];

        ts_yyyy_xxxxyy[i] = 3.0 * ts_yy_xxxxyy[i] * fe_0 + 2.0 * ts_yyy_xxxxy[i] * fe_0 + ts_yyy_xxxxyy[i] * pa_y[i];

        ts_yyyy_xxxxyz[i] = 3.0 * ts_yy_xxxxyz[i] * fe_0 + ts_yyy_xxxxz[i] * fe_0 + ts_yyy_xxxxyz[i] * pa_y[i];

        ts_yyyy_xxxxzz[i] = 3.0 * ts_yy_xxxxzz[i] * fe_0 + ts_yyy_xxxxzz[i] * pa_y[i];

        ts_yyyy_xxxyyy[i] = 3.0 * ts_yy_xxxyyy[i] * fe_0 + 3.0 * ts_yyy_xxxyy[i] * fe_0 + ts_yyy_xxxyyy[i] * pa_y[i];

        ts_yyyy_xxxyyz[i] = 3.0 * ts_yy_xxxyyz[i] * fe_0 + 2.0 * ts_yyy_xxxyz[i] * fe_0 + ts_yyy_xxxyyz[i] * pa_y[i];

        ts_yyyy_xxxyzz[i] = 3.0 * ts_yy_xxxyzz[i] * fe_0 + ts_yyy_xxxzz[i] * fe_0 + ts_yyy_xxxyzz[i] * pa_y[i];

        ts_yyyy_xxxzzz[i] = 3.0 * ts_yy_xxxzzz[i] * fe_0 + ts_yyy_xxxzzz[i] * pa_y[i];

        ts_yyyy_xxyyyy[i] = 3.0 * ts_yy_xxyyyy[i] * fe_0 + 4.0 * ts_yyy_xxyyy[i] * fe_0 + ts_yyy_xxyyyy[i] * pa_y[i];

        ts_yyyy_xxyyyz[i] = 3.0 * ts_yy_xxyyyz[i] * fe_0 + 3.0 * ts_yyy_xxyyz[i] * fe_0 + ts_yyy_xxyyyz[i] * pa_y[i];

        ts_yyyy_xxyyzz[i] = 3.0 * ts_yy_xxyyzz[i] * fe_0 + 2.0 * ts_yyy_xxyzz[i] * fe_0 + ts_yyy_xxyyzz[i] * pa_y[i];

        ts_yyyy_xxyzzz[i] = 3.0 * ts_yy_xxyzzz[i] * fe_0 + ts_yyy_xxzzz[i] * fe_0 + ts_yyy_xxyzzz[i] * pa_y[i];

        ts_yyyy_xxzzzz[i] = 3.0 * ts_yy_xxzzzz[i] * fe_0 + ts_yyy_xxzzzz[i] * pa_y[i];

        ts_yyyy_xyyyyy[i] = 3.0 * ts_yy_xyyyyy[i] * fe_0 + 5.0 * ts_yyy_xyyyy[i] * fe_0 + ts_yyy_xyyyyy[i] * pa_y[i];

        ts_yyyy_xyyyyz[i] = 3.0 * ts_yy_xyyyyz[i] * fe_0 + 4.0 * ts_yyy_xyyyz[i] * fe_0 + ts_yyy_xyyyyz[i] * pa_y[i];

        ts_yyyy_xyyyzz[i] = 3.0 * ts_yy_xyyyzz[i] * fe_0 + 3.0 * ts_yyy_xyyzz[i] * fe_0 + ts_yyy_xyyyzz[i] * pa_y[i];

        ts_yyyy_xyyzzz[i] = 3.0 * ts_yy_xyyzzz[i] * fe_0 + 2.0 * ts_yyy_xyzzz[i] * fe_0 + ts_yyy_xyyzzz[i] * pa_y[i];

        ts_yyyy_xyzzzz[i] = 3.0 * ts_yy_xyzzzz[i] * fe_0 + ts_yyy_xzzzz[i] * fe_0 + ts_yyy_xyzzzz[i] * pa_y[i];

        ts_yyyy_xzzzzz[i] = 3.0 * ts_yy_xzzzzz[i] * fe_0 + ts_yyy_xzzzzz[i] * pa_y[i];

        ts_yyyy_yyyyyy[i] = 3.0 * ts_yy_yyyyyy[i] * fe_0 + 6.0 * ts_yyy_yyyyy[i] * fe_0 + ts_yyy_yyyyyy[i] * pa_y[i];

        ts_yyyy_yyyyyz[i] = 3.0 * ts_yy_yyyyyz[i] * fe_0 + 5.0 * ts_yyy_yyyyz[i] * fe_0 + ts_yyy_yyyyyz[i] * pa_y[i];

        ts_yyyy_yyyyzz[i] = 3.0 * ts_yy_yyyyzz[i] * fe_0 + 4.0 * ts_yyy_yyyzz[i] * fe_0 + ts_yyy_yyyyzz[i] * pa_y[i];

        ts_yyyy_yyyzzz[i] = 3.0 * ts_yy_yyyzzz[i] * fe_0 + 3.0 * ts_yyy_yyzzz[i] * fe_0 + ts_yyy_yyyzzz[i] * pa_y[i];

        ts_yyyy_yyzzzz[i] = 3.0 * ts_yy_yyzzzz[i] * fe_0 + 2.0 * ts_yyy_yzzzz[i] * fe_0 + ts_yyy_yyzzzz[i] * pa_y[i];

        ts_yyyy_yzzzzz[i] = 3.0 * ts_yy_yzzzzz[i] * fe_0 + ts_yyy_zzzzz[i] * fe_0 + ts_yyy_yzzzzz[i] * pa_y[i];

        ts_yyyy_zzzzzz[i] = 3.0 * ts_yy_zzzzzz[i] * fe_0 + ts_yyy_zzzzzz[i] * pa_y[i];
    }

    // Set up 308-336 components of targeted buffer : GI

    auto ts_yyyz_xxxxxx = pbuffer.data(idx_ovl_gi + 308);

    auto ts_yyyz_xxxxxy = pbuffer.data(idx_ovl_gi + 309);

    auto ts_yyyz_xxxxxz = pbuffer.data(idx_ovl_gi + 310);

    auto ts_yyyz_xxxxyy = pbuffer.data(idx_ovl_gi + 311);

    auto ts_yyyz_xxxxyz = pbuffer.data(idx_ovl_gi + 312);

    auto ts_yyyz_xxxxzz = pbuffer.data(idx_ovl_gi + 313);

    auto ts_yyyz_xxxyyy = pbuffer.data(idx_ovl_gi + 314);

    auto ts_yyyz_xxxyyz = pbuffer.data(idx_ovl_gi + 315);

    auto ts_yyyz_xxxyzz = pbuffer.data(idx_ovl_gi + 316);

    auto ts_yyyz_xxxzzz = pbuffer.data(idx_ovl_gi + 317);

    auto ts_yyyz_xxyyyy = pbuffer.data(idx_ovl_gi + 318);

    auto ts_yyyz_xxyyyz = pbuffer.data(idx_ovl_gi + 319);

    auto ts_yyyz_xxyyzz = pbuffer.data(idx_ovl_gi + 320);

    auto ts_yyyz_xxyzzz = pbuffer.data(idx_ovl_gi + 321);

    auto ts_yyyz_xxzzzz = pbuffer.data(idx_ovl_gi + 322);

    auto ts_yyyz_xyyyyy = pbuffer.data(idx_ovl_gi + 323);

    auto ts_yyyz_xyyyyz = pbuffer.data(idx_ovl_gi + 324);

    auto ts_yyyz_xyyyzz = pbuffer.data(idx_ovl_gi + 325);

    auto ts_yyyz_xyyzzz = pbuffer.data(idx_ovl_gi + 326);

    auto ts_yyyz_xyzzzz = pbuffer.data(idx_ovl_gi + 327);

    auto ts_yyyz_xzzzzz = pbuffer.data(idx_ovl_gi + 328);

    auto ts_yyyz_yyyyyy = pbuffer.data(idx_ovl_gi + 329);

    auto ts_yyyz_yyyyyz = pbuffer.data(idx_ovl_gi + 330);

    auto ts_yyyz_yyyyzz = pbuffer.data(idx_ovl_gi + 331);

    auto ts_yyyz_yyyzzz = pbuffer.data(idx_ovl_gi + 332);

    auto ts_yyyz_yyzzzz = pbuffer.data(idx_ovl_gi + 333);

    auto ts_yyyz_yzzzzz = pbuffer.data(idx_ovl_gi + 334);

    auto ts_yyyz_zzzzzz = pbuffer.data(idx_ovl_gi + 335);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyy_xxxxxx, ts_yyy_xxxxxy, ts_yyy_xxxxy, ts_yyy_xxxxyy, ts_yyy_xxxxyz, ts_yyy_xxxyy, ts_yyy_xxxyyy, ts_yyy_xxxyyz, ts_yyy_xxxyz, ts_yyy_xxxyzz, ts_yyy_xxyyy, ts_yyy_xxyyyy, ts_yyy_xxyyyz, ts_yyy_xxyyz, ts_yyy_xxyyzz, ts_yyy_xxyzz, ts_yyy_xxyzzz, ts_yyy_xyyyy, ts_yyy_xyyyyy, ts_yyy_xyyyyz, ts_yyy_xyyyz, ts_yyy_xyyyzz, ts_yyy_xyyzz, ts_yyy_xyyzzz, ts_yyy_xyzzz, ts_yyy_xyzzzz, ts_yyy_yyyyy, ts_yyy_yyyyyy, ts_yyy_yyyyyz, ts_yyy_yyyyz, ts_yyy_yyyyzz, ts_yyy_yyyzz, ts_yyy_yyyzzz, ts_yyy_yyzzz, ts_yyy_yyzzzz, ts_yyy_yzzzz, ts_yyy_yzzzzz, ts_yyyz_xxxxxx, ts_yyyz_xxxxxy, ts_yyyz_xxxxxz, ts_yyyz_xxxxyy, ts_yyyz_xxxxyz, ts_yyyz_xxxxzz, ts_yyyz_xxxyyy, ts_yyyz_xxxyyz, ts_yyyz_xxxyzz, ts_yyyz_xxxzzz, ts_yyyz_xxyyyy, ts_yyyz_xxyyyz, ts_yyyz_xxyyzz, ts_yyyz_xxyzzz, ts_yyyz_xxzzzz, ts_yyyz_xyyyyy, ts_yyyz_xyyyyz, ts_yyyz_xyyyzz, ts_yyyz_xyyzzz, ts_yyyz_xyzzzz, ts_yyyz_xzzzzz, ts_yyyz_yyyyyy, ts_yyyz_yyyyyz, ts_yyyz_yyyyzz, ts_yyyz_yyyzzz, ts_yyyz_yyzzzz, ts_yyyz_yzzzzz, ts_yyyz_zzzzzz, ts_yyz_xxxxxz, ts_yyz_xxxxzz, ts_yyz_xxxzzz, ts_yyz_xxzzzz, ts_yyz_xzzzzz, ts_yyz_zzzzzz, ts_yz_xxxxxz, ts_yz_xxxxzz, ts_yz_xxxzzz, ts_yz_xxzzzz, ts_yz_xzzzzz, ts_yz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyz_xxxxxx[i] = ts_yyy_xxxxxx[i] * pa_z[i];

        ts_yyyz_xxxxxy[i] = ts_yyy_xxxxxy[i] * pa_z[i];

        ts_yyyz_xxxxxz[i] = 2.0 * ts_yz_xxxxxz[i] * fe_0 + ts_yyz_xxxxxz[i] * pa_y[i];

        ts_yyyz_xxxxyy[i] = ts_yyy_xxxxyy[i] * pa_z[i];

        ts_yyyz_xxxxyz[i] = ts_yyy_xxxxy[i] * fe_0 + ts_yyy_xxxxyz[i] * pa_z[i];

        ts_yyyz_xxxxzz[i] = 2.0 * ts_yz_xxxxzz[i] * fe_0 + ts_yyz_xxxxzz[i] * pa_y[i];

        ts_yyyz_xxxyyy[i] = ts_yyy_xxxyyy[i] * pa_z[i];

        ts_yyyz_xxxyyz[i] = ts_yyy_xxxyy[i] * fe_0 + ts_yyy_xxxyyz[i] * pa_z[i];

        ts_yyyz_xxxyzz[i] = 2.0 * ts_yyy_xxxyz[i] * fe_0 + ts_yyy_xxxyzz[i] * pa_z[i];

        ts_yyyz_xxxzzz[i] = 2.0 * ts_yz_xxxzzz[i] * fe_0 + ts_yyz_xxxzzz[i] * pa_y[i];

        ts_yyyz_xxyyyy[i] = ts_yyy_xxyyyy[i] * pa_z[i];

        ts_yyyz_xxyyyz[i] = ts_yyy_xxyyy[i] * fe_0 + ts_yyy_xxyyyz[i] * pa_z[i];

        ts_yyyz_xxyyzz[i] = 2.0 * ts_yyy_xxyyz[i] * fe_0 + ts_yyy_xxyyzz[i] * pa_z[i];

        ts_yyyz_xxyzzz[i] = 3.0 * ts_yyy_xxyzz[i] * fe_0 + ts_yyy_xxyzzz[i] * pa_z[i];

        ts_yyyz_xxzzzz[i] = 2.0 * ts_yz_xxzzzz[i] * fe_0 + ts_yyz_xxzzzz[i] * pa_y[i];

        ts_yyyz_xyyyyy[i] = ts_yyy_xyyyyy[i] * pa_z[i];

        ts_yyyz_xyyyyz[i] = ts_yyy_xyyyy[i] * fe_0 + ts_yyy_xyyyyz[i] * pa_z[i];

        ts_yyyz_xyyyzz[i] = 2.0 * ts_yyy_xyyyz[i] * fe_0 + ts_yyy_xyyyzz[i] * pa_z[i];

        ts_yyyz_xyyzzz[i] = 3.0 * ts_yyy_xyyzz[i] * fe_0 + ts_yyy_xyyzzz[i] * pa_z[i];

        ts_yyyz_xyzzzz[i] = 4.0 * ts_yyy_xyzzz[i] * fe_0 + ts_yyy_xyzzzz[i] * pa_z[i];

        ts_yyyz_xzzzzz[i] = 2.0 * ts_yz_xzzzzz[i] * fe_0 + ts_yyz_xzzzzz[i] * pa_y[i];

        ts_yyyz_yyyyyy[i] = ts_yyy_yyyyyy[i] * pa_z[i];

        ts_yyyz_yyyyyz[i] = ts_yyy_yyyyy[i] * fe_0 + ts_yyy_yyyyyz[i] * pa_z[i];

        ts_yyyz_yyyyzz[i] = 2.0 * ts_yyy_yyyyz[i] * fe_0 + ts_yyy_yyyyzz[i] * pa_z[i];

        ts_yyyz_yyyzzz[i] = 3.0 * ts_yyy_yyyzz[i] * fe_0 + ts_yyy_yyyzzz[i] * pa_z[i];

        ts_yyyz_yyzzzz[i] = 4.0 * ts_yyy_yyzzz[i] * fe_0 + ts_yyy_yyzzzz[i] * pa_z[i];

        ts_yyyz_yzzzzz[i] = 5.0 * ts_yyy_yzzzz[i] * fe_0 + ts_yyy_yzzzzz[i] * pa_z[i];

        ts_yyyz_zzzzzz[i] = 2.0 * ts_yz_zzzzzz[i] * fe_0 + ts_yyz_zzzzzz[i] * pa_y[i];
    }

    // Set up 336-364 components of targeted buffer : GI

    auto ts_yyzz_xxxxxx = pbuffer.data(idx_ovl_gi + 336);

    auto ts_yyzz_xxxxxy = pbuffer.data(idx_ovl_gi + 337);

    auto ts_yyzz_xxxxxz = pbuffer.data(idx_ovl_gi + 338);

    auto ts_yyzz_xxxxyy = pbuffer.data(idx_ovl_gi + 339);

    auto ts_yyzz_xxxxyz = pbuffer.data(idx_ovl_gi + 340);

    auto ts_yyzz_xxxxzz = pbuffer.data(idx_ovl_gi + 341);

    auto ts_yyzz_xxxyyy = pbuffer.data(idx_ovl_gi + 342);

    auto ts_yyzz_xxxyyz = pbuffer.data(idx_ovl_gi + 343);

    auto ts_yyzz_xxxyzz = pbuffer.data(idx_ovl_gi + 344);

    auto ts_yyzz_xxxzzz = pbuffer.data(idx_ovl_gi + 345);

    auto ts_yyzz_xxyyyy = pbuffer.data(idx_ovl_gi + 346);

    auto ts_yyzz_xxyyyz = pbuffer.data(idx_ovl_gi + 347);

    auto ts_yyzz_xxyyzz = pbuffer.data(idx_ovl_gi + 348);

    auto ts_yyzz_xxyzzz = pbuffer.data(idx_ovl_gi + 349);

    auto ts_yyzz_xxzzzz = pbuffer.data(idx_ovl_gi + 350);

    auto ts_yyzz_xyyyyy = pbuffer.data(idx_ovl_gi + 351);

    auto ts_yyzz_xyyyyz = pbuffer.data(idx_ovl_gi + 352);

    auto ts_yyzz_xyyyzz = pbuffer.data(idx_ovl_gi + 353);

    auto ts_yyzz_xyyzzz = pbuffer.data(idx_ovl_gi + 354);

    auto ts_yyzz_xyzzzz = pbuffer.data(idx_ovl_gi + 355);

    auto ts_yyzz_xzzzzz = pbuffer.data(idx_ovl_gi + 356);

    auto ts_yyzz_yyyyyy = pbuffer.data(idx_ovl_gi + 357);

    auto ts_yyzz_yyyyyz = pbuffer.data(idx_ovl_gi + 358);

    auto ts_yyzz_yyyyzz = pbuffer.data(idx_ovl_gi + 359);

    auto ts_yyzz_yyyzzz = pbuffer.data(idx_ovl_gi + 360);

    auto ts_yyzz_yyzzzz = pbuffer.data(idx_ovl_gi + 361);

    auto ts_yyzz_yzzzzz = pbuffer.data(idx_ovl_gi + 362);

    auto ts_yyzz_zzzzzz = pbuffer.data(idx_ovl_gi + 363);

    #pragma omp simd aligned(pa_y, pa_z, ts_yy_xxxxxy, ts_yy_xxxxyy, ts_yy_xxxyyy, ts_yy_xxyyyy, ts_yy_xyyyyy, ts_yy_yyyyyy, ts_yyz_xxxxxy, ts_yyz_xxxxyy, ts_yyz_xxxyyy, ts_yyz_xxyyyy, ts_yyz_xyyyyy, ts_yyz_yyyyyy, ts_yyzz_xxxxxx, ts_yyzz_xxxxxy, ts_yyzz_xxxxxz, ts_yyzz_xxxxyy, ts_yyzz_xxxxyz, ts_yyzz_xxxxzz, ts_yyzz_xxxyyy, ts_yyzz_xxxyyz, ts_yyzz_xxxyzz, ts_yyzz_xxxzzz, ts_yyzz_xxyyyy, ts_yyzz_xxyyyz, ts_yyzz_xxyyzz, ts_yyzz_xxyzzz, ts_yyzz_xxzzzz, ts_yyzz_xyyyyy, ts_yyzz_xyyyyz, ts_yyzz_xyyyzz, ts_yyzz_xyyzzz, ts_yyzz_xyzzzz, ts_yyzz_xzzzzz, ts_yyzz_yyyyyy, ts_yyzz_yyyyyz, ts_yyzz_yyyyzz, ts_yyzz_yyyzzz, ts_yyzz_yyzzzz, ts_yyzz_yzzzzz, ts_yyzz_zzzzzz, ts_yzz_xxxxxx, ts_yzz_xxxxxz, ts_yzz_xxxxyz, ts_yzz_xxxxz, ts_yzz_xxxxzz, ts_yzz_xxxyyz, ts_yzz_xxxyz, ts_yzz_xxxyzz, ts_yzz_xxxzz, ts_yzz_xxxzzz, ts_yzz_xxyyyz, ts_yzz_xxyyz, ts_yzz_xxyyzz, ts_yzz_xxyzz, ts_yzz_xxyzzz, ts_yzz_xxzzz, ts_yzz_xxzzzz, ts_yzz_xyyyyz, ts_yzz_xyyyz, ts_yzz_xyyyzz, ts_yzz_xyyzz, ts_yzz_xyyzzz, ts_yzz_xyzzz, ts_yzz_xyzzzz, ts_yzz_xzzzz, ts_yzz_xzzzzz, ts_yzz_yyyyyz, ts_yzz_yyyyz, ts_yzz_yyyyzz, ts_yzz_yyyzz, ts_yzz_yyyzzz, ts_yzz_yyzzz, ts_yzz_yyzzzz, ts_yzz_yzzzz, ts_yzz_yzzzzz, ts_yzz_zzzzz, ts_yzz_zzzzzz, ts_zz_xxxxxx, ts_zz_xxxxxz, ts_zz_xxxxyz, ts_zz_xxxxzz, ts_zz_xxxyyz, ts_zz_xxxyzz, ts_zz_xxxzzz, ts_zz_xxyyyz, ts_zz_xxyyzz, ts_zz_xxyzzz, ts_zz_xxzzzz, ts_zz_xyyyyz, ts_zz_xyyyzz, ts_zz_xyyzzz, ts_zz_xyzzzz, ts_zz_xzzzzz, ts_zz_yyyyyz, ts_zz_yyyyzz, ts_zz_yyyzzz, ts_zz_yyzzzz, ts_zz_yzzzzz, ts_zz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzz_xxxxxx[i] = ts_zz_xxxxxx[i] * fe_0 + ts_yzz_xxxxxx[i] * pa_y[i];

        ts_yyzz_xxxxxy[i] = ts_yy_xxxxxy[i] * fe_0 + ts_yyz_xxxxxy[i] * pa_z[i];

        ts_yyzz_xxxxxz[i] = ts_zz_xxxxxz[i] * fe_0 + ts_yzz_xxxxxz[i] * pa_y[i];

        ts_yyzz_xxxxyy[i] = ts_yy_xxxxyy[i] * fe_0 + ts_yyz_xxxxyy[i] * pa_z[i];

        ts_yyzz_xxxxyz[i] = ts_zz_xxxxyz[i] * fe_0 + ts_yzz_xxxxz[i] * fe_0 + ts_yzz_xxxxyz[i] * pa_y[i];

        ts_yyzz_xxxxzz[i] = ts_zz_xxxxzz[i] * fe_0 + ts_yzz_xxxxzz[i] * pa_y[i];

        ts_yyzz_xxxyyy[i] = ts_yy_xxxyyy[i] * fe_0 + ts_yyz_xxxyyy[i] * pa_z[i];

        ts_yyzz_xxxyyz[i] = ts_zz_xxxyyz[i] * fe_0 + 2.0 * ts_yzz_xxxyz[i] * fe_0 + ts_yzz_xxxyyz[i] * pa_y[i];

        ts_yyzz_xxxyzz[i] = ts_zz_xxxyzz[i] * fe_0 + ts_yzz_xxxzz[i] * fe_0 + ts_yzz_xxxyzz[i] * pa_y[i];

        ts_yyzz_xxxzzz[i] = ts_zz_xxxzzz[i] * fe_0 + ts_yzz_xxxzzz[i] * pa_y[i];

        ts_yyzz_xxyyyy[i] = ts_yy_xxyyyy[i] * fe_0 + ts_yyz_xxyyyy[i] * pa_z[i];

        ts_yyzz_xxyyyz[i] = ts_zz_xxyyyz[i] * fe_0 + 3.0 * ts_yzz_xxyyz[i] * fe_0 + ts_yzz_xxyyyz[i] * pa_y[i];

        ts_yyzz_xxyyzz[i] = ts_zz_xxyyzz[i] * fe_0 + 2.0 * ts_yzz_xxyzz[i] * fe_0 + ts_yzz_xxyyzz[i] * pa_y[i];

        ts_yyzz_xxyzzz[i] = ts_zz_xxyzzz[i] * fe_0 + ts_yzz_xxzzz[i] * fe_0 + ts_yzz_xxyzzz[i] * pa_y[i];

        ts_yyzz_xxzzzz[i] = ts_zz_xxzzzz[i] * fe_0 + ts_yzz_xxzzzz[i] * pa_y[i];

        ts_yyzz_xyyyyy[i] = ts_yy_xyyyyy[i] * fe_0 + ts_yyz_xyyyyy[i] * pa_z[i];

        ts_yyzz_xyyyyz[i] = ts_zz_xyyyyz[i] * fe_0 + 4.0 * ts_yzz_xyyyz[i] * fe_0 + ts_yzz_xyyyyz[i] * pa_y[i];

        ts_yyzz_xyyyzz[i] = ts_zz_xyyyzz[i] * fe_0 + 3.0 * ts_yzz_xyyzz[i] * fe_0 + ts_yzz_xyyyzz[i] * pa_y[i];

        ts_yyzz_xyyzzz[i] = ts_zz_xyyzzz[i] * fe_0 + 2.0 * ts_yzz_xyzzz[i] * fe_0 + ts_yzz_xyyzzz[i] * pa_y[i];

        ts_yyzz_xyzzzz[i] = ts_zz_xyzzzz[i] * fe_0 + ts_yzz_xzzzz[i] * fe_0 + ts_yzz_xyzzzz[i] * pa_y[i];

        ts_yyzz_xzzzzz[i] = ts_zz_xzzzzz[i] * fe_0 + ts_yzz_xzzzzz[i] * pa_y[i];

        ts_yyzz_yyyyyy[i] = ts_yy_yyyyyy[i] * fe_0 + ts_yyz_yyyyyy[i] * pa_z[i];

        ts_yyzz_yyyyyz[i] = ts_zz_yyyyyz[i] * fe_0 + 5.0 * ts_yzz_yyyyz[i] * fe_0 + ts_yzz_yyyyyz[i] * pa_y[i];

        ts_yyzz_yyyyzz[i] = ts_zz_yyyyzz[i] * fe_0 + 4.0 * ts_yzz_yyyzz[i] * fe_0 + ts_yzz_yyyyzz[i] * pa_y[i];

        ts_yyzz_yyyzzz[i] = ts_zz_yyyzzz[i] * fe_0 + 3.0 * ts_yzz_yyzzz[i] * fe_0 + ts_yzz_yyyzzz[i] * pa_y[i];

        ts_yyzz_yyzzzz[i] = ts_zz_yyzzzz[i] * fe_0 + 2.0 * ts_yzz_yzzzz[i] * fe_0 + ts_yzz_yyzzzz[i] * pa_y[i];

        ts_yyzz_yzzzzz[i] = ts_zz_yzzzzz[i] * fe_0 + ts_yzz_zzzzz[i] * fe_0 + ts_yzz_yzzzzz[i] * pa_y[i];

        ts_yyzz_zzzzzz[i] = ts_zz_zzzzzz[i] * fe_0 + ts_yzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 364-392 components of targeted buffer : GI

    auto ts_yzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 364);

    auto ts_yzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 365);

    auto ts_yzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 366);

    auto ts_yzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 367);

    auto ts_yzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 368);

    auto ts_yzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 369);

    auto ts_yzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 370);

    auto ts_yzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 371);

    auto ts_yzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 372);

    auto ts_yzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 373);

    auto ts_yzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 374);

    auto ts_yzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 375);

    auto ts_yzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 376);

    auto ts_yzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 377);

    auto ts_yzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 378);

    auto ts_yzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 379);

    auto ts_yzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 380);

    auto ts_yzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 381);

    auto ts_yzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 382);

    auto ts_yzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 383);

    auto ts_yzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 384);

    auto ts_yzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 385);

    auto ts_yzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 386);

    auto ts_yzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 387);

    auto ts_yzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 388);

    auto ts_yzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 389);

    auto ts_yzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 390);

    auto ts_yzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 391);

    #pragma omp simd aligned(pa_y, ts_yzzz_xxxxxx, ts_yzzz_xxxxxy, ts_yzzz_xxxxxz, ts_yzzz_xxxxyy, ts_yzzz_xxxxyz, ts_yzzz_xxxxzz, ts_yzzz_xxxyyy, ts_yzzz_xxxyyz, ts_yzzz_xxxyzz, ts_yzzz_xxxzzz, ts_yzzz_xxyyyy, ts_yzzz_xxyyyz, ts_yzzz_xxyyzz, ts_yzzz_xxyzzz, ts_yzzz_xxzzzz, ts_yzzz_xyyyyy, ts_yzzz_xyyyyz, ts_yzzz_xyyyzz, ts_yzzz_xyyzzz, ts_yzzz_xyzzzz, ts_yzzz_xzzzzz, ts_yzzz_yyyyyy, ts_yzzz_yyyyyz, ts_yzzz_yyyyzz, ts_yzzz_yyyzzz, ts_yzzz_yyzzzz, ts_yzzz_yzzzzz, ts_yzzz_zzzzzz, ts_zzz_xxxxx, ts_zzz_xxxxxx, ts_zzz_xxxxxy, ts_zzz_xxxxxz, ts_zzz_xxxxy, ts_zzz_xxxxyy, ts_zzz_xxxxyz, ts_zzz_xxxxz, ts_zzz_xxxxzz, ts_zzz_xxxyy, ts_zzz_xxxyyy, ts_zzz_xxxyyz, ts_zzz_xxxyz, ts_zzz_xxxyzz, ts_zzz_xxxzz, ts_zzz_xxxzzz, ts_zzz_xxyyy, ts_zzz_xxyyyy, ts_zzz_xxyyyz, ts_zzz_xxyyz, ts_zzz_xxyyzz, ts_zzz_xxyzz, ts_zzz_xxyzzz, ts_zzz_xxzzz, ts_zzz_xxzzzz, ts_zzz_xyyyy, ts_zzz_xyyyyy, ts_zzz_xyyyyz, ts_zzz_xyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyy, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzz, ts_zzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzz_xxxxxx[i] = ts_zzz_xxxxxx[i] * pa_y[i];

        ts_yzzz_xxxxxy[i] = ts_zzz_xxxxx[i] * fe_0 + ts_zzz_xxxxxy[i] * pa_y[i];

        ts_yzzz_xxxxxz[i] = ts_zzz_xxxxxz[i] * pa_y[i];

        ts_yzzz_xxxxyy[i] = 2.0 * ts_zzz_xxxxy[i] * fe_0 + ts_zzz_xxxxyy[i] * pa_y[i];

        ts_yzzz_xxxxyz[i] = ts_zzz_xxxxz[i] * fe_0 + ts_zzz_xxxxyz[i] * pa_y[i];

        ts_yzzz_xxxxzz[i] = ts_zzz_xxxxzz[i] * pa_y[i];

        ts_yzzz_xxxyyy[i] = 3.0 * ts_zzz_xxxyy[i] * fe_0 + ts_zzz_xxxyyy[i] * pa_y[i];

        ts_yzzz_xxxyyz[i] = 2.0 * ts_zzz_xxxyz[i] * fe_0 + ts_zzz_xxxyyz[i] * pa_y[i];

        ts_yzzz_xxxyzz[i] = ts_zzz_xxxzz[i] * fe_0 + ts_zzz_xxxyzz[i] * pa_y[i];

        ts_yzzz_xxxzzz[i] = ts_zzz_xxxzzz[i] * pa_y[i];

        ts_yzzz_xxyyyy[i] = 4.0 * ts_zzz_xxyyy[i] * fe_0 + ts_zzz_xxyyyy[i] * pa_y[i];

        ts_yzzz_xxyyyz[i] = 3.0 * ts_zzz_xxyyz[i] * fe_0 + ts_zzz_xxyyyz[i] * pa_y[i];

        ts_yzzz_xxyyzz[i] = 2.0 * ts_zzz_xxyzz[i] * fe_0 + ts_zzz_xxyyzz[i] * pa_y[i];

        ts_yzzz_xxyzzz[i] = ts_zzz_xxzzz[i] * fe_0 + ts_zzz_xxyzzz[i] * pa_y[i];

        ts_yzzz_xxzzzz[i] = ts_zzz_xxzzzz[i] * pa_y[i];

        ts_yzzz_xyyyyy[i] = 5.0 * ts_zzz_xyyyy[i] * fe_0 + ts_zzz_xyyyyy[i] * pa_y[i];

        ts_yzzz_xyyyyz[i] = 4.0 * ts_zzz_xyyyz[i] * fe_0 + ts_zzz_xyyyyz[i] * pa_y[i];

        ts_yzzz_xyyyzz[i] = 3.0 * ts_zzz_xyyzz[i] * fe_0 + ts_zzz_xyyyzz[i] * pa_y[i];

        ts_yzzz_xyyzzz[i] = 2.0 * ts_zzz_xyzzz[i] * fe_0 + ts_zzz_xyyzzz[i] * pa_y[i];

        ts_yzzz_xyzzzz[i] = ts_zzz_xzzzz[i] * fe_0 + ts_zzz_xyzzzz[i] * pa_y[i];

        ts_yzzz_xzzzzz[i] = ts_zzz_xzzzzz[i] * pa_y[i];

        ts_yzzz_yyyyyy[i] = 6.0 * ts_zzz_yyyyy[i] * fe_0 + ts_zzz_yyyyyy[i] * pa_y[i];

        ts_yzzz_yyyyyz[i] = 5.0 * ts_zzz_yyyyz[i] * fe_0 + ts_zzz_yyyyyz[i] * pa_y[i];

        ts_yzzz_yyyyzz[i] = 4.0 * ts_zzz_yyyzz[i] * fe_0 + ts_zzz_yyyyzz[i] * pa_y[i];

        ts_yzzz_yyyzzz[i] = 3.0 * ts_zzz_yyzzz[i] * fe_0 + ts_zzz_yyyzzz[i] * pa_y[i];

        ts_yzzz_yyzzzz[i] = 2.0 * ts_zzz_yzzzz[i] * fe_0 + ts_zzz_yyzzzz[i] * pa_y[i];

        ts_yzzz_yzzzzz[i] = ts_zzz_zzzzz[i] * fe_0 + ts_zzz_yzzzzz[i] * pa_y[i];

        ts_yzzz_zzzzzz[i] = ts_zzz_zzzzzz[i] * pa_y[i];
    }

    // Set up 392-420 components of targeted buffer : GI

    auto ts_zzzz_xxxxxx = pbuffer.data(idx_ovl_gi + 392);

    auto ts_zzzz_xxxxxy = pbuffer.data(idx_ovl_gi + 393);

    auto ts_zzzz_xxxxxz = pbuffer.data(idx_ovl_gi + 394);

    auto ts_zzzz_xxxxyy = pbuffer.data(idx_ovl_gi + 395);

    auto ts_zzzz_xxxxyz = pbuffer.data(idx_ovl_gi + 396);

    auto ts_zzzz_xxxxzz = pbuffer.data(idx_ovl_gi + 397);

    auto ts_zzzz_xxxyyy = pbuffer.data(idx_ovl_gi + 398);

    auto ts_zzzz_xxxyyz = pbuffer.data(idx_ovl_gi + 399);

    auto ts_zzzz_xxxyzz = pbuffer.data(idx_ovl_gi + 400);

    auto ts_zzzz_xxxzzz = pbuffer.data(idx_ovl_gi + 401);

    auto ts_zzzz_xxyyyy = pbuffer.data(idx_ovl_gi + 402);

    auto ts_zzzz_xxyyyz = pbuffer.data(idx_ovl_gi + 403);

    auto ts_zzzz_xxyyzz = pbuffer.data(idx_ovl_gi + 404);

    auto ts_zzzz_xxyzzz = pbuffer.data(idx_ovl_gi + 405);

    auto ts_zzzz_xxzzzz = pbuffer.data(idx_ovl_gi + 406);

    auto ts_zzzz_xyyyyy = pbuffer.data(idx_ovl_gi + 407);

    auto ts_zzzz_xyyyyz = pbuffer.data(idx_ovl_gi + 408);

    auto ts_zzzz_xyyyzz = pbuffer.data(idx_ovl_gi + 409);

    auto ts_zzzz_xyyzzz = pbuffer.data(idx_ovl_gi + 410);

    auto ts_zzzz_xyzzzz = pbuffer.data(idx_ovl_gi + 411);

    auto ts_zzzz_xzzzzz = pbuffer.data(idx_ovl_gi + 412);

    auto ts_zzzz_yyyyyy = pbuffer.data(idx_ovl_gi + 413);

    auto ts_zzzz_yyyyyz = pbuffer.data(idx_ovl_gi + 414);

    auto ts_zzzz_yyyyzz = pbuffer.data(idx_ovl_gi + 415);

    auto ts_zzzz_yyyzzz = pbuffer.data(idx_ovl_gi + 416);

    auto ts_zzzz_yyzzzz = pbuffer.data(idx_ovl_gi + 417);

    auto ts_zzzz_yzzzzz = pbuffer.data(idx_ovl_gi + 418);

    auto ts_zzzz_zzzzzz = pbuffer.data(idx_ovl_gi + 419);

    #pragma omp simd aligned(pa_z, ts_zz_xxxxxx, ts_zz_xxxxxy, ts_zz_xxxxxz, ts_zz_xxxxyy, ts_zz_xxxxyz, ts_zz_xxxxzz, ts_zz_xxxyyy, ts_zz_xxxyyz, ts_zz_xxxyzz, ts_zz_xxxzzz, ts_zz_xxyyyy, ts_zz_xxyyyz, ts_zz_xxyyzz, ts_zz_xxyzzz, ts_zz_xxzzzz, ts_zz_xyyyyy, ts_zz_xyyyyz, ts_zz_xyyyzz, ts_zz_xyyzzz, ts_zz_xyzzzz, ts_zz_xzzzzz, ts_zz_yyyyyy, ts_zz_yyyyyz, ts_zz_yyyyzz, ts_zz_yyyzzz, ts_zz_yyzzzz, ts_zz_yzzzzz, ts_zz_zzzzzz, ts_zzz_xxxxx, ts_zzz_xxxxxx, ts_zzz_xxxxxy, ts_zzz_xxxxxz, ts_zzz_xxxxy, ts_zzz_xxxxyy, ts_zzz_xxxxyz, ts_zzz_xxxxz, ts_zzz_xxxxzz, ts_zzz_xxxyy, ts_zzz_xxxyyy, ts_zzz_xxxyyz, ts_zzz_xxxyz, ts_zzz_xxxyzz, ts_zzz_xxxzz, ts_zzz_xxxzzz, ts_zzz_xxyyy, ts_zzz_xxyyyy, ts_zzz_xxyyyz, ts_zzz_xxyyz, ts_zzz_xxyyzz, ts_zzz_xxyzz, ts_zzz_xxyzzz, ts_zzz_xxzzz, ts_zzz_xxzzzz, ts_zzz_xyyyy, ts_zzz_xyyyyy, ts_zzz_xyyyyz, ts_zzz_xyyyz, ts_zzz_xyyyzz, ts_zzz_xyyzz, ts_zzz_xyyzzz, ts_zzz_xyzzz, ts_zzz_xyzzzz, ts_zzz_xzzzz, ts_zzz_xzzzzz, ts_zzz_yyyyy, ts_zzz_yyyyyy, ts_zzz_yyyyyz, ts_zzz_yyyyz, ts_zzz_yyyyzz, ts_zzz_yyyzz, ts_zzz_yyyzzz, ts_zzz_yyzzz, ts_zzz_yyzzzz, ts_zzz_yzzzz, ts_zzz_yzzzzz, ts_zzz_zzzzz, ts_zzz_zzzzzz, ts_zzzz_xxxxxx, ts_zzzz_xxxxxy, ts_zzzz_xxxxxz, ts_zzzz_xxxxyy, ts_zzzz_xxxxyz, ts_zzzz_xxxxzz, ts_zzzz_xxxyyy, ts_zzzz_xxxyyz, ts_zzzz_xxxyzz, ts_zzzz_xxxzzz, ts_zzzz_xxyyyy, ts_zzzz_xxyyyz, ts_zzzz_xxyyzz, ts_zzzz_xxyzzz, ts_zzzz_xxzzzz, ts_zzzz_xyyyyy, ts_zzzz_xyyyyz, ts_zzzz_xyyyzz, ts_zzzz_xyyzzz, ts_zzzz_xyzzzz, ts_zzzz_xzzzzz, ts_zzzz_yyyyyy, ts_zzzz_yyyyyz, ts_zzzz_yyyyzz, ts_zzzz_yyyzzz, ts_zzzz_yyzzzz, ts_zzzz_yzzzzz, ts_zzzz_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzz_xxxxxx[i] = 3.0 * ts_zz_xxxxxx[i] * fe_0 + ts_zzz_xxxxxx[i] * pa_z[i];

        ts_zzzz_xxxxxy[i] = 3.0 * ts_zz_xxxxxy[i] * fe_0 + ts_zzz_xxxxxy[i] * pa_z[i];

        ts_zzzz_xxxxxz[i] = 3.0 * ts_zz_xxxxxz[i] * fe_0 + ts_zzz_xxxxx[i] * fe_0 + ts_zzz_xxxxxz[i] * pa_z[i];

        ts_zzzz_xxxxyy[i] = 3.0 * ts_zz_xxxxyy[i] * fe_0 + ts_zzz_xxxxyy[i] * pa_z[i];

        ts_zzzz_xxxxyz[i] = 3.0 * ts_zz_xxxxyz[i] * fe_0 + ts_zzz_xxxxy[i] * fe_0 + ts_zzz_xxxxyz[i] * pa_z[i];

        ts_zzzz_xxxxzz[i] = 3.0 * ts_zz_xxxxzz[i] * fe_0 + 2.0 * ts_zzz_xxxxz[i] * fe_0 + ts_zzz_xxxxzz[i] * pa_z[i];

        ts_zzzz_xxxyyy[i] = 3.0 * ts_zz_xxxyyy[i] * fe_0 + ts_zzz_xxxyyy[i] * pa_z[i];

        ts_zzzz_xxxyyz[i] = 3.0 * ts_zz_xxxyyz[i] * fe_0 + ts_zzz_xxxyy[i] * fe_0 + ts_zzz_xxxyyz[i] * pa_z[i];

        ts_zzzz_xxxyzz[i] = 3.0 * ts_zz_xxxyzz[i] * fe_0 + 2.0 * ts_zzz_xxxyz[i] * fe_0 + ts_zzz_xxxyzz[i] * pa_z[i];

        ts_zzzz_xxxzzz[i] = 3.0 * ts_zz_xxxzzz[i] * fe_0 + 3.0 * ts_zzz_xxxzz[i] * fe_0 + ts_zzz_xxxzzz[i] * pa_z[i];

        ts_zzzz_xxyyyy[i] = 3.0 * ts_zz_xxyyyy[i] * fe_0 + ts_zzz_xxyyyy[i] * pa_z[i];

        ts_zzzz_xxyyyz[i] = 3.0 * ts_zz_xxyyyz[i] * fe_0 + ts_zzz_xxyyy[i] * fe_0 + ts_zzz_xxyyyz[i] * pa_z[i];

        ts_zzzz_xxyyzz[i] = 3.0 * ts_zz_xxyyzz[i] * fe_0 + 2.0 * ts_zzz_xxyyz[i] * fe_0 + ts_zzz_xxyyzz[i] * pa_z[i];

        ts_zzzz_xxyzzz[i] = 3.0 * ts_zz_xxyzzz[i] * fe_0 + 3.0 * ts_zzz_xxyzz[i] * fe_0 + ts_zzz_xxyzzz[i] * pa_z[i];

        ts_zzzz_xxzzzz[i] = 3.0 * ts_zz_xxzzzz[i] * fe_0 + 4.0 * ts_zzz_xxzzz[i] * fe_0 + ts_zzz_xxzzzz[i] * pa_z[i];

        ts_zzzz_xyyyyy[i] = 3.0 * ts_zz_xyyyyy[i] * fe_0 + ts_zzz_xyyyyy[i] * pa_z[i];

        ts_zzzz_xyyyyz[i] = 3.0 * ts_zz_xyyyyz[i] * fe_0 + ts_zzz_xyyyy[i] * fe_0 + ts_zzz_xyyyyz[i] * pa_z[i];

        ts_zzzz_xyyyzz[i] = 3.0 * ts_zz_xyyyzz[i] * fe_0 + 2.0 * ts_zzz_xyyyz[i] * fe_0 + ts_zzz_xyyyzz[i] * pa_z[i];

        ts_zzzz_xyyzzz[i] = 3.0 * ts_zz_xyyzzz[i] * fe_0 + 3.0 * ts_zzz_xyyzz[i] * fe_0 + ts_zzz_xyyzzz[i] * pa_z[i];

        ts_zzzz_xyzzzz[i] = 3.0 * ts_zz_xyzzzz[i] * fe_0 + 4.0 * ts_zzz_xyzzz[i] * fe_0 + ts_zzz_xyzzzz[i] * pa_z[i];

        ts_zzzz_xzzzzz[i] = 3.0 * ts_zz_xzzzzz[i] * fe_0 + 5.0 * ts_zzz_xzzzz[i] * fe_0 + ts_zzz_xzzzzz[i] * pa_z[i];

        ts_zzzz_yyyyyy[i] = 3.0 * ts_zz_yyyyyy[i] * fe_0 + ts_zzz_yyyyyy[i] * pa_z[i];

        ts_zzzz_yyyyyz[i] = 3.0 * ts_zz_yyyyyz[i] * fe_0 + ts_zzz_yyyyy[i] * fe_0 + ts_zzz_yyyyyz[i] * pa_z[i];

        ts_zzzz_yyyyzz[i] = 3.0 * ts_zz_yyyyzz[i] * fe_0 + 2.0 * ts_zzz_yyyyz[i] * fe_0 + ts_zzz_yyyyzz[i] * pa_z[i];

        ts_zzzz_yyyzzz[i] = 3.0 * ts_zz_yyyzzz[i] * fe_0 + 3.0 * ts_zzz_yyyzz[i] * fe_0 + ts_zzz_yyyzzz[i] * pa_z[i];

        ts_zzzz_yyzzzz[i] = 3.0 * ts_zz_yyzzzz[i] * fe_0 + 4.0 * ts_zzz_yyzzz[i] * fe_0 + ts_zzz_yyzzzz[i] * pa_z[i];

        ts_zzzz_yzzzzz[i] = 3.0 * ts_zz_yzzzzz[i] * fe_0 + 5.0 * ts_zzz_yzzzz[i] * fe_0 + ts_zzz_yzzzzz[i] * pa_z[i];

        ts_zzzz_zzzzzz[i] = 3.0 * ts_zz_zzzzzz[i] * fe_0 + 6.0 * ts_zzz_zzzzz[i] * fe_0 + ts_zzz_zzzzzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

