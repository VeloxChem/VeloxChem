#include "NuclearPotentialPrimRecGI.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_gi(CSimdArray<double>& pbuffer, 
                               const size_t idx_npot_0_gi,
                               const size_t idx_npot_0_di,
                               const size_t idx_npot_1_di,
                               const size_t idx_npot_0_fh,
                               const size_t idx_npot_1_fh,
                               const size_t idx_npot_0_fi,
                               const size_t idx_npot_1_fi,
                               const CSimdArray<double>& factors,
                               const size_t idx_rpa,
                               const size_t idx_rpc,
                               const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : DI

    auto ta_xx_xxxxxx_0 = pbuffer.data(idx_npot_0_di);

    auto ta_xx_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 1);

    auto ta_xx_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 2);

    auto ta_xx_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 3);

    auto ta_xx_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 4);

    auto ta_xx_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 5);

    auto ta_xx_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 6);

    auto ta_xx_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 7);

    auto ta_xx_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 8);

    auto ta_xx_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 9);

    auto ta_xx_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 10);

    auto ta_xx_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 11);

    auto ta_xx_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 12);

    auto ta_xx_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 13);

    auto ta_xx_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 14);

    auto ta_xx_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 15);

    auto ta_xx_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 16);

    auto ta_xx_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 17);

    auto ta_xx_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 18);

    auto ta_xx_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 19);

    auto ta_xx_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 20);

    auto ta_xx_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 21);

    auto ta_xx_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 22);

    auto ta_xx_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 23);

    auto ta_xx_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 24);

    auto ta_xx_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 25);

    auto ta_xx_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 26);

    auto ta_xx_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 27);

    auto ta_xy_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 49);

    auto ta_xy_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 50);

    auto ta_xy_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 51);

    auto ta_xy_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 52);

    auto ta_xy_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 53);

    auto ta_xy_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 54);

    auto ta_xz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 78);

    auto ta_xz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 79);

    auto ta_xz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 80);

    auto ta_xz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 81);

    auto ta_xz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 82);

    auto ta_xz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 83);

    auto ta_yy_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 84);

    auto ta_yy_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 85);

    auto ta_yy_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 86);

    auto ta_yy_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 87);

    auto ta_yy_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 88);

    auto ta_yy_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 89);

    auto ta_yy_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 90);

    auto ta_yy_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 91);

    auto ta_yy_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 92);

    auto ta_yy_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 93);

    auto ta_yy_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 94);

    auto ta_yy_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 95);

    auto ta_yy_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 96);

    auto ta_yy_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 97);

    auto ta_yy_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 98);

    auto ta_yy_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 99);

    auto ta_yy_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 100);

    auto ta_yy_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 101);

    auto ta_yy_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 102);

    auto ta_yy_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 103);

    auto ta_yy_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 104);

    auto ta_yy_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 105);

    auto ta_yy_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 106);

    auto ta_yy_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 107);

    auto ta_yy_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 108);

    auto ta_yy_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 109);

    auto ta_yy_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 110);

    auto ta_yy_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 111);

    auto ta_yz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 114);

    auto ta_yz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 117);

    auto ta_yz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 121);

    auto ta_yz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 126);

    auto ta_yz_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 132);

    auto ta_yz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 134);

    auto ta_yz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 135);

    auto ta_yz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 136);

    auto ta_yz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 137);

    auto ta_yz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 138);

    auto ta_yz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 139);

    auto ta_zz_xxxxxx_0 = pbuffer.data(idx_npot_0_di + 140);

    auto ta_zz_xxxxxy_0 = pbuffer.data(idx_npot_0_di + 141);

    auto ta_zz_xxxxxz_0 = pbuffer.data(idx_npot_0_di + 142);

    auto ta_zz_xxxxyy_0 = pbuffer.data(idx_npot_0_di + 143);

    auto ta_zz_xxxxyz_0 = pbuffer.data(idx_npot_0_di + 144);

    auto ta_zz_xxxxzz_0 = pbuffer.data(idx_npot_0_di + 145);

    auto ta_zz_xxxyyy_0 = pbuffer.data(idx_npot_0_di + 146);

    auto ta_zz_xxxyyz_0 = pbuffer.data(idx_npot_0_di + 147);

    auto ta_zz_xxxyzz_0 = pbuffer.data(idx_npot_0_di + 148);

    auto ta_zz_xxxzzz_0 = pbuffer.data(idx_npot_0_di + 149);

    auto ta_zz_xxyyyy_0 = pbuffer.data(idx_npot_0_di + 150);

    auto ta_zz_xxyyyz_0 = pbuffer.data(idx_npot_0_di + 151);

    auto ta_zz_xxyyzz_0 = pbuffer.data(idx_npot_0_di + 152);

    auto ta_zz_xxyzzz_0 = pbuffer.data(idx_npot_0_di + 153);

    auto ta_zz_xxzzzz_0 = pbuffer.data(idx_npot_0_di + 154);

    auto ta_zz_xyyyyy_0 = pbuffer.data(idx_npot_0_di + 155);

    auto ta_zz_xyyyyz_0 = pbuffer.data(idx_npot_0_di + 156);

    auto ta_zz_xyyyzz_0 = pbuffer.data(idx_npot_0_di + 157);

    auto ta_zz_xyyzzz_0 = pbuffer.data(idx_npot_0_di + 158);

    auto ta_zz_xyzzzz_0 = pbuffer.data(idx_npot_0_di + 159);

    auto ta_zz_xzzzzz_0 = pbuffer.data(idx_npot_0_di + 160);

    auto ta_zz_yyyyyy_0 = pbuffer.data(idx_npot_0_di + 161);

    auto ta_zz_yyyyyz_0 = pbuffer.data(idx_npot_0_di + 162);

    auto ta_zz_yyyyzz_0 = pbuffer.data(idx_npot_0_di + 163);

    auto ta_zz_yyyzzz_0 = pbuffer.data(idx_npot_0_di + 164);

    auto ta_zz_yyzzzz_0 = pbuffer.data(idx_npot_0_di + 165);

    auto ta_zz_yzzzzz_0 = pbuffer.data(idx_npot_0_di + 166);

    auto ta_zz_zzzzzz_0 = pbuffer.data(idx_npot_0_di + 167);

    // Set up components of auxiliary buffer : DI

    auto ta_xx_xxxxxx_1 = pbuffer.data(idx_npot_1_di);

    auto ta_xx_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 1);

    auto ta_xx_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 2);

    auto ta_xx_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 3);

    auto ta_xx_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 4);

    auto ta_xx_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 5);

    auto ta_xx_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 6);

    auto ta_xx_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 7);

    auto ta_xx_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 8);

    auto ta_xx_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 9);

    auto ta_xx_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 10);

    auto ta_xx_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 11);

    auto ta_xx_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 12);

    auto ta_xx_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 13);

    auto ta_xx_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 14);

    auto ta_xx_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 15);

    auto ta_xx_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 16);

    auto ta_xx_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 17);

    auto ta_xx_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 18);

    auto ta_xx_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 19);

    auto ta_xx_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 20);

    auto ta_xx_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 21);

    auto ta_xx_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 22);

    auto ta_xx_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 23);

    auto ta_xx_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 24);

    auto ta_xx_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 25);

    auto ta_xx_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 26);

    auto ta_xx_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 27);

    auto ta_xy_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 49);

    auto ta_xy_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 50);

    auto ta_xy_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 51);

    auto ta_xy_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 52);

    auto ta_xy_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 53);

    auto ta_xy_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 54);

    auto ta_xz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 78);

    auto ta_xz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 79);

    auto ta_xz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 80);

    auto ta_xz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 81);

    auto ta_xz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 82);

    auto ta_xz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 83);

    auto ta_yy_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 84);

    auto ta_yy_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 85);

    auto ta_yy_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 86);

    auto ta_yy_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 87);

    auto ta_yy_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 88);

    auto ta_yy_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 89);

    auto ta_yy_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 90);

    auto ta_yy_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 91);

    auto ta_yy_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 92);

    auto ta_yy_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 93);

    auto ta_yy_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 94);

    auto ta_yy_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 95);

    auto ta_yy_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 96);

    auto ta_yy_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 97);

    auto ta_yy_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 98);

    auto ta_yy_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 99);

    auto ta_yy_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 100);

    auto ta_yy_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 101);

    auto ta_yy_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 102);

    auto ta_yy_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 103);

    auto ta_yy_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 104);

    auto ta_yy_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 105);

    auto ta_yy_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 106);

    auto ta_yy_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 107);

    auto ta_yy_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 108);

    auto ta_yy_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 109);

    auto ta_yy_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 110);

    auto ta_yy_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 111);

    auto ta_yz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 114);

    auto ta_yz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 117);

    auto ta_yz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 121);

    auto ta_yz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 126);

    auto ta_yz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 132);

    auto ta_yz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 134);

    auto ta_yz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 135);

    auto ta_yz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 136);

    auto ta_yz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 137);

    auto ta_yz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 138);

    auto ta_yz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 139);

    auto ta_zz_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 140);

    auto ta_zz_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 141);

    auto ta_zz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 142);

    auto ta_zz_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 143);

    auto ta_zz_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 144);

    auto ta_zz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 145);

    auto ta_zz_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 146);

    auto ta_zz_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 147);

    auto ta_zz_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 148);

    auto ta_zz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 149);

    auto ta_zz_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 150);

    auto ta_zz_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 151);

    auto ta_zz_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 152);

    auto ta_zz_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 153);

    auto ta_zz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 154);

    auto ta_zz_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 155);

    auto ta_zz_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 156);

    auto ta_zz_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 157);

    auto ta_zz_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 158);

    auto ta_zz_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 159);

    auto ta_zz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 160);

    auto ta_zz_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 161);

    auto ta_zz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 162);

    auto ta_zz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 163);

    auto ta_zz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 164);

    auto ta_zz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 165);

    auto ta_zz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 166);

    auto ta_zz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 167);

    // Set up components of auxiliary buffer : FH

    auto ta_xxx_xxxxx_0 = pbuffer.data(idx_npot_0_fh);

    auto ta_xxx_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 1);

    auto ta_xxx_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 2);

    auto ta_xxx_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 3);

    auto ta_xxx_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 4);

    auto ta_xxx_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 5);

    auto ta_xxx_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 6);

    auto ta_xxx_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 7);

    auto ta_xxx_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 8);

    auto ta_xxx_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 9);

    auto ta_xxx_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 10);

    auto ta_xxx_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 11);

    auto ta_xxx_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 12);

    auto ta_xxx_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 13);

    auto ta_xxx_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 14);

    auto ta_xxx_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 15);

    auto ta_xxx_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 16);

    auto ta_xxx_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 17);

    auto ta_xxx_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 18);

    auto ta_xxx_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 19);

    auto ta_xxx_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 20);

    auto ta_xxz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 44);

    auto ta_xxz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 46);

    auto ta_xxz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 47);

    auto ta_xxz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 49);

    auto ta_xxz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 50);

    auto ta_xxz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 51);

    auto ta_xxz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 53);

    auto ta_xxz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 54);

    auto ta_xxz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 55);

    auto ta_xxz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 56);

    auto ta_xyy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 64);

    auto ta_xyy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 66);

    auto ta_xyy_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 67);

    auto ta_xyy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 69);

    auto ta_xyy_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 70);

    auto ta_xyy_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 71);

    auto ta_xyy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 73);

    auto ta_xyy_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 74);

    auto ta_xyy_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 75);

    auto ta_xyy_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 76);

    auto ta_xyy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 78);

    auto ta_xyy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 79);

    auto ta_xyy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 80);

    auto ta_xyy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 81);

    auto ta_xyy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 82);

    auto ta_xzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 107);

    auto ta_xzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 109);

    auto ta_xzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 110);

    auto ta_xzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 112);

    auto ta_xzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 113);

    auto ta_xzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 114);

    auto ta_xzz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 116);

    auto ta_xzz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 117);

    auto ta_xzz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 118);

    auto ta_xzz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 119);

    auto ta_xzz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 121);

    auto ta_xzz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 122);

    auto ta_xzz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 123);

    auto ta_xzz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 124);

    auto ta_xzz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 125);

    auto ta_yyy_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 126);

    auto ta_yyy_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 127);

    auto ta_yyy_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 128);

    auto ta_yyy_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 129);

    auto ta_yyy_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 130);

    auto ta_yyy_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 131);

    auto ta_yyy_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 132);

    auto ta_yyy_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 133);

    auto ta_yyy_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 134);

    auto ta_yyy_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 135);

    auto ta_yyy_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 136);

    auto ta_yyy_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 137);

    auto ta_yyy_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 138);

    auto ta_yyy_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 139);

    auto ta_yyy_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 140);

    auto ta_yyy_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 141);

    auto ta_yyy_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 142);

    auto ta_yyy_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 143);

    auto ta_yyy_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 144);

    auto ta_yyy_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 145);

    auto ta_yyy_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 146);

    auto ta_yyz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 149);

    auto ta_yyz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 151);

    auto ta_yyz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 152);

    auto ta_yyz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 154);

    auto ta_yyz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 155);

    auto ta_yyz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 156);

    auto ta_yyz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 158);

    auto ta_yyz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 159);

    auto ta_yyz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 160);

    auto ta_yyz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 161);

    auto ta_yyz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 163);

    auto ta_yyz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 164);

    auto ta_yyz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 165);

    auto ta_yyz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 166);

    auto ta_yyz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 167);

    auto ta_yzz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 169);

    auto ta_yzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 170);

    auto ta_yzz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 171);

    auto ta_yzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 172);

    auto ta_yzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 173);

    auto ta_yzz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 174);

    auto ta_yzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 175);

    auto ta_yzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 176);

    auto ta_yzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 177);

    auto ta_yzz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 178);

    auto ta_yzz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 179);

    auto ta_yzz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 180);

    auto ta_yzz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 181);

    auto ta_yzz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 182);

    auto ta_yzz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 183);

    auto ta_yzz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 184);

    auto ta_yzz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 185);

    auto ta_yzz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 186);

    auto ta_yzz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 187);

    auto ta_yzz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 188);

    auto ta_zzz_xxxxx_0 = pbuffer.data(idx_npot_0_fh + 189);

    auto ta_zzz_xxxxy_0 = pbuffer.data(idx_npot_0_fh + 190);

    auto ta_zzz_xxxxz_0 = pbuffer.data(idx_npot_0_fh + 191);

    auto ta_zzz_xxxyy_0 = pbuffer.data(idx_npot_0_fh + 192);

    auto ta_zzz_xxxyz_0 = pbuffer.data(idx_npot_0_fh + 193);

    auto ta_zzz_xxxzz_0 = pbuffer.data(idx_npot_0_fh + 194);

    auto ta_zzz_xxyyy_0 = pbuffer.data(idx_npot_0_fh + 195);

    auto ta_zzz_xxyyz_0 = pbuffer.data(idx_npot_0_fh + 196);

    auto ta_zzz_xxyzz_0 = pbuffer.data(idx_npot_0_fh + 197);

    auto ta_zzz_xxzzz_0 = pbuffer.data(idx_npot_0_fh + 198);

    auto ta_zzz_xyyyy_0 = pbuffer.data(idx_npot_0_fh + 199);

    auto ta_zzz_xyyyz_0 = pbuffer.data(idx_npot_0_fh + 200);

    auto ta_zzz_xyyzz_0 = pbuffer.data(idx_npot_0_fh + 201);

    auto ta_zzz_xyzzz_0 = pbuffer.data(idx_npot_0_fh + 202);

    auto ta_zzz_xzzzz_0 = pbuffer.data(idx_npot_0_fh + 203);

    auto ta_zzz_yyyyy_0 = pbuffer.data(idx_npot_0_fh + 204);

    auto ta_zzz_yyyyz_0 = pbuffer.data(idx_npot_0_fh + 205);

    auto ta_zzz_yyyzz_0 = pbuffer.data(idx_npot_0_fh + 206);

    auto ta_zzz_yyzzz_0 = pbuffer.data(idx_npot_0_fh + 207);

    auto ta_zzz_yzzzz_0 = pbuffer.data(idx_npot_0_fh + 208);

    auto ta_zzz_zzzzz_0 = pbuffer.data(idx_npot_0_fh + 209);

    // Set up components of auxiliary buffer : FH

    auto ta_xxx_xxxxx_1 = pbuffer.data(idx_npot_1_fh);

    auto ta_xxx_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 1);

    auto ta_xxx_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 2);

    auto ta_xxx_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 3);

    auto ta_xxx_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 4);

    auto ta_xxx_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 5);

    auto ta_xxx_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 6);

    auto ta_xxx_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 7);

    auto ta_xxx_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 8);

    auto ta_xxx_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 9);

    auto ta_xxx_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 10);

    auto ta_xxx_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 11);

    auto ta_xxx_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 12);

    auto ta_xxx_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 13);

    auto ta_xxx_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 14);

    auto ta_xxx_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 15);

    auto ta_xxx_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 16);

    auto ta_xxx_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 17);

    auto ta_xxx_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 18);

    auto ta_xxx_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 19);

    auto ta_xxx_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 20);

    auto ta_xxz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 44);

    auto ta_xxz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 46);

    auto ta_xxz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 47);

    auto ta_xxz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 49);

    auto ta_xxz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 50);

    auto ta_xxz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 51);

    auto ta_xxz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 53);

    auto ta_xxz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 54);

    auto ta_xxz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 55);

    auto ta_xxz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 56);

    auto ta_xyy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 64);

    auto ta_xyy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 66);

    auto ta_xyy_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 67);

    auto ta_xyy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 69);

    auto ta_xyy_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 70);

    auto ta_xyy_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 71);

    auto ta_xyy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 73);

    auto ta_xyy_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 74);

    auto ta_xyy_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 75);

    auto ta_xyy_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 76);

    auto ta_xyy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 78);

    auto ta_xyy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 79);

    auto ta_xyy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 80);

    auto ta_xyy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 81);

    auto ta_xyy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 82);

    auto ta_xzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 107);

    auto ta_xzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 109);

    auto ta_xzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 110);

    auto ta_xzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 112);

    auto ta_xzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 113);

    auto ta_xzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 114);

    auto ta_xzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 116);

    auto ta_xzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 117);

    auto ta_xzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 118);

    auto ta_xzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 119);

    auto ta_xzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 121);

    auto ta_xzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 122);

    auto ta_xzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 123);

    auto ta_xzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 124);

    auto ta_xzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 125);

    auto ta_yyy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 126);

    auto ta_yyy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 127);

    auto ta_yyy_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 128);

    auto ta_yyy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 129);

    auto ta_yyy_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 130);

    auto ta_yyy_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 131);

    auto ta_yyy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 132);

    auto ta_yyy_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 133);

    auto ta_yyy_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 134);

    auto ta_yyy_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 135);

    auto ta_yyy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 136);

    auto ta_yyy_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 137);

    auto ta_yyy_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 138);

    auto ta_yyy_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 139);

    auto ta_yyy_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 140);

    auto ta_yyy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 141);

    auto ta_yyy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 142);

    auto ta_yyy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 143);

    auto ta_yyy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 144);

    auto ta_yyy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 145);

    auto ta_yyy_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 146);

    auto ta_yyz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 149);

    auto ta_yyz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 151);

    auto ta_yyz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 152);

    auto ta_yyz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 154);

    auto ta_yyz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 155);

    auto ta_yyz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 156);

    auto ta_yyz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 158);

    auto ta_yyz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 159);

    auto ta_yyz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 160);

    auto ta_yyz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 161);

    auto ta_yyz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 163);

    auto ta_yyz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 164);

    auto ta_yyz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 165);

    auto ta_yyz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 166);

    auto ta_yyz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 167);

    auto ta_yzz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 169);

    auto ta_yzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 170);

    auto ta_yzz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 171);

    auto ta_yzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 172);

    auto ta_yzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 173);

    auto ta_yzz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 174);

    auto ta_yzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 175);

    auto ta_yzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 176);

    auto ta_yzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 177);

    auto ta_yzz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 178);

    auto ta_yzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 179);

    auto ta_yzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 180);

    auto ta_yzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 181);

    auto ta_yzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 182);

    auto ta_yzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 183);

    auto ta_yzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 184);

    auto ta_yzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 185);

    auto ta_yzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 186);

    auto ta_yzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 187);

    auto ta_yzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 188);

    auto ta_zzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 189);

    auto ta_zzz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 190);

    auto ta_zzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 191);

    auto ta_zzz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 192);

    auto ta_zzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 193);

    auto ta_zzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 194);

    auto ta_zzz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 195);

    auto ta_zzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 196);

    auto ta_zzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 197);

    auto ta_zzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 198);

    auto ta_zzz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 199);

    auto ta_zzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 200);

    auto ta_zzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 201);

    auto ta_zzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 202);

    auto ta_zzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 203);

    auto ta_zzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 204);

    auto ta_zzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 205);

    auto ta_zzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 206);

    auto ta_zzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 207);

    auto ta_zzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 208);

    auto ta_zzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 209);

    // Set up components of auxiliary buffer : FI

    auto ta_xxx_xxxxxx_0 = pbuffer.data(idx_npot_0_fi);

    auto ta_xxx_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 1);

    auto ta_xxx_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 2);

    auto ta_xxx_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 3);

    auto ta_xxx_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 4);

    auto ta_xxx_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 5);

    auto ta_xxx_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 6);

    auto ta_xxx_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 7);

    auto ta_xxx_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 8);

    auto ta_xxx_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 9);

    auto ta_xxx_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 10);

    auto ta_xxx_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 11);

    auto ta_xxx_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 12);

    auto ta_xxx_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 13);

    auto ta_xxx_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 14);

    auto ta_xxx_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 15);

    auto ta_xxx_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 16);

    auto ta_xxx_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 17);

    auto ta_xxx_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 18);

    auto ta_xxx_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 19);

    auto ta_xxx_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 20);

    auto ta_xxx_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 21);

    auto ta_xxx_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 22);

    auto ta_xxx_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 23);

    auto ta_xxx_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 24);

    auto ta_xxx_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 25);

    auto ta_xxx_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 26);

    auto ta_xxx_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 27);

    auto ta_xxy_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 28);

    auto ta_xxy_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 29);

    auto ta_xxy_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 30);

    auto ta_xxy_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 31);

    auto ta_xxy_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 33);

    auto ta_xxy_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 34);

    auto ta_xxy_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 37);

    auto ta_xxy_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 38);

    auto ta_xxy_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 42);

    auto ta_xxy_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 43);

    auto ta_xxy_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 48);

    auto ta_xxy_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 49);

    auto ta_xxy_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 50);

    auto ta_xxy_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 51);

    auto ta_xxy_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 52);

    auto ta_xxy_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 53);

    auto ta_xxy_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 54);

    auto ta_xxz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 56);

    auto ta_xxz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 57);

    auto ta_xxz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 58);

    auto ta_xxz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 59);

    auto ta_xxz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 60);

    auto ta_xxz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 61);

    auto ta_xxz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 62);

    auto ta_xxz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 63);

    auto ta_xxz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 64);

    auto ta_xxz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 65);

    auto ta_xxz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 66);

    auto ta_xxz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 67);

    auto ta_xxz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 68);

    auto ta_xxz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 69);

    auto ta_xxz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 70);

    auto ta_xxz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 71);

    auto ta_xxz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 72);

    auto ta_xxz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 73);

    auto ta_xxz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 74);

    auto ta_xxz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 75);

    auto ta_xxz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 76);

    auto ta_xxz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 78);

    auto ta_xxz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 79);

    auto ta_xxz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 80);

    auto ta_xxz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 81);

    auto ta_xxz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 82);

    auto ta_xxz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 83);

    auto ta_xyy_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 84);

    auto ta_xyy_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 85);

    auto ta_xyy_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 87);

    auto ta_xyy_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 88);

    auto ta_xyy_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 90);

    auto ta_xyy_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 91);

    auto ta_xyy_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 92);

    auto ta_xyy_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 94);

    auto ta_xyy_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 95);

    auto ta_xyy_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 96);

    auto ta_xyy_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 97);

    auto ta_xyy_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 99);

    auto ta_xyy_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 100);

    auto ta_xyy_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 101);

    auto ta_xyy_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 102);

    auto ta_xyy_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 103);

    auto ta_xyy_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 105);

    auto ta_xyy_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 106);

    auto ta_xyy_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 107);

    auto ta_xyy_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 108);

    auto ta_xyy_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 109);

    auto ta_xyy_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 110);

    auto ta_xyy_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 111);

    auto ta_xyz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 134);

    auto ta_xyz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 135);

    auto ta_xyz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 136);

    auto ta_xyz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 137);

    auto ta_xyz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 138);

    auto ta_xzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 140);

    auto ta_xzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 142);

    auto ta_xzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 144);

    auto ta_xzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 145);

    auto ta_xzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 147);

    auto ta_xzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 148);

    auto ta_xzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 149);

    auto ta_xzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 151);

    auto ta_xzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 152);

    auto ta_xzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 153);

    auto ta_xzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 154);

    auto ta_xzz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 156);

    auto ta_xzz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 157);

    auto ta_xzz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 158);

    auto ta_xzz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 159);

    auto ta_xzz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 160);

    auto ta_xzz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 161);

    auto ta_xzz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 162);

    auto ta_xzz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 163);

    auto ta_xzz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 164);

    auto ta_xzz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 165);

    auto ta_xzz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 166);

    auto ta_xzz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 167);

    auto ta_yyy_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 168);

    auto ta_yyy_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 169);

    auto ta_yyy_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 170);

    auto ta_yyy_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 171);

    auto ta_yyy_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 172);

    auto ta_yyy_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 173);

    auto ta_yyy_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 174);

    auto ta_yyy_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 175);

    auto ta_yyy_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 176);

    auto ta_yyy_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 177);

    auto ta_yyy_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 178);

    auto ta_yyy_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 179);

    auto ta_yyy_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 180);

    auto ta_yyy_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 181);

    auto ta_yyy_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 182);

    auto ta_yyy_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 183);

    auto ta_yyy_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 184);

    auto ta_yyy_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 185);

    auto ta_yyy_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 186);

    auto ta_yyy_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 187);

    auto ta_yyy_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 188);

    auto ta_yyy_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 189);

    auto ta_yyy_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 190);

    auto ta_yyy_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 191);

    auto ta_yyy_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 192);

    auto ta_yyy_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 193);

    auto ta_yyy_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 194);

    auto ta_yyy_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 195);

    auto ta_yyz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 197);

    auto ta_yyz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 198);

    auto ta_yyz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 199);

    auto ta_yyz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 200);

    auto ta_yyz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 201);

    auto ta_yyz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 202);

    auto ta_yyz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 203);

    auto ta_yyz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 204);

    auto ta_yyz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 205);

    auto ta_yyz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 206);

    auto ta_yyz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 207);

    auto ta_yyz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 208);

    auto ta_yyz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 209);

    auto ta_yyz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 210);

    auto ta_yyz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 211);

    auto ta_yyz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 212);

    auto ta_yyz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 213);

    auto ta_yyz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 214);

    auto ta_yyz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 215);

    auto ta_yyz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 216);

    auto ta_yyz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 217);

    auto ta_yyz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 218);

    auto ta_yyz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 219);

    auto ta_yyz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 220);

    auto ta_yyz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 221);

    auto ta_yyz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 222);

    auto ta_yyz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 223);

    auto ta_yzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 224);

    auto ta_yzz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 225);

    auto ta_yzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 226);

    auto ta_yzz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 227);

    auto ta_yzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 228);

    auto ta_yzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 229);

    auto ta_yzz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 230);

    auto ta_yzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 231);

    auto ta_yzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 232);

    auto ta_yzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 233);

    auto ta_yzz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 234);

    auto ta_yzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 235);

    auto ta_yzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 236);

    auto ta_yzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 237);

    auto ta_yzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 238);

    auto ta_yzz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 239);

    auto ta_yzz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 240);

    auto ta_yzz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 241);

    auto ta_yzz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 242);

    auto ta_yzz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 243);

    auto ta_yzz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 244);

    auto ta_yzz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 245);

    auto ta_yzz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 246);

    auto ta_yzz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 247);

    auto ta_yzz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 248);

    auto ta_yzz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 249);

    auto ta_yzz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 250);

    auto ta_yzz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 251);

    auto ta_zzz_xxxxxx_0 = pbuffer.data(idx_npot_0_fi + 252);

    auto ta_zzz_xxxxxy_0 = pbuffer.data(idx_npot_0_fi + 253);

    auto ta_zzz_xxxxxz_0 = pbuffer.data(idx_npot_0_fi + 254);

    auto ta_zzz_xxxxyy_0 = pbuffer.data(idx_npot_0_fi + 255);

    auto ta_zzz_xxxxyz_0 = pbuffer.data(idx_npot_0_fi + 256);

    auto ta_zzz_xxxxzz_0 = pbuffer.data(idx_npot_0_fi + 257);

    auto ta_zzz_xxxyyy_0 = pbuffer.data(idx_npot_0_fi + 258);

    auto ta_zzz_xxxyyz_0 = pbuffer.data(idx_npot_0_fi + 259);

    auto ta_zzz_xxxyzz_0 = pbuffer.data(idx_npot_0_fi + 260);

    auto ta_zzz_xxxzzz_0 = pbuffer.data(idx_npot_0_fi + 261);

    auto ta_zzz_xxyyyy_0 = pbuffer.data(idx_npot_0_fi + 262);

    auto ta_zzz_xxyyyz_0 = pbuffer.data(idx_npot_0_fi + 263);

    auto ta_zzz_xxyyzz_0 = pbuffer.data(idx_npot_0_fi + 264);

    auto ta_zzz_xxyzzz_0 = pbuffer.data(idx_npot_0_fi + 265);

    auto ta_zzz_xxzzzz_0 = pbuffer.data(idx_npot_0_fi + 266);

    auto ta_zzz_xyyyyy_0 = pbuffer.data(idx_npot_0_fi + 267);

    auto ta_zzz_xyyyyz_0 = pbuffer.data(idx_npot_0_fi + 268);

    auto ta_zzz_xyyyzz_0 = pbuffer.data(idx_npot_0_fi + 269);

    auto ta_zzz_xyyzzz_0 = pbuffer.data(idx_npot_0_fi + 270);

    auto ta_zzz_xyzzzz_0 = pbuffer.data(idx_npot_0_fi + 271);

    auto ta_zzz_xzzzzz_0 = pbuffer.data(idx_npot_0_fi + 272);

    auto ta_zzz_yyyyyy_0 = pbuffer.data(idx_npot_0_fi + 273);

    auto ta_zzz_yyyyyz_0 = pbuffer.data(idx_npot_0_fi + 274);

    auto ta_zzz_yyyyzz_0 = pbuffer.data(idx_npot_0_fi + 275);

    auto ta_zzz_yyyzzz_0 = pbuffer.data(idx_npot_0_fi + 276);

    auto ta_zzz_yyzzzz_0 = pbuffer.data(idx_npot_0_fi + 277);

    auto ta_zzz_yzzzzz_0 = pbuffer.data(idx_npot_0_fi + 278);

    auto ta_zzz_zzzzzz_0 = pbuffer.data(idx_npot_0_fi + 279);

    // Set up components of auxiliary buffer : FI

    auto ta_xxx_xxxxxx_1 = pbuffer.data(idx_npot_1_fi);

    auto ta_xxx_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 1);

    auto ta_xxx_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 2);

    auto ta_xxx_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 3);

    auto ta_xxx_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 4);

    auto ta_xxx_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 5);

    auto ta_xxx_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 6);

    auto ta_xxx_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 7);

    auto ta_xxx_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 8);

    auto ta_xxx_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 9);

    auto ta_xxx_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 10);

    auto ta_xxx_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 11);

    auto ta_xxx_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 12);

    auto ta_xxx_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 13);

    auto ta_xxx_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 14);

    auto ta_xxx_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 15);

    auto ta_xxx_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 16);

    auto ta_xxx_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 17);

    auto ta_xxx_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 18);

    auto ta_xxx_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 19);

    auto ta_xxx_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 20);

    auto ta_xxx_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 21);

    auto ta_xxx_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 22);

    auto ta_xxx_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 23);

    auto ta_xxx_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 24);

    auto ta_xxx_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 25);

    auto ta_xxx_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 26);

    auto ta_xxx_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 27);

    auto ta_xxy_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 28);

    auto ta_xxy_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 29);

    auto ta_xxy_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 30);

    auto ta_xxy_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 31);

    auto ta_xxy_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 33);

    auto ta_xxy_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 34);

    auto ta_xxy_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 37);

    auto ta_xxy_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 38);

    auto ta_xxy_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 42);

    auto ta_xxy_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 43);

    auto ta_xxy_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 48);

    auto ta_xxy_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 49);

    auto ta_xxy_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 50);

    auto ta_xxy_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 51);

    auto ta_xxy_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 52);

    auto ta_xxy_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 53);

    auto ta_xxy_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 54);

    auto ta_xxz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 56);

    auto ta_xxz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 57);

    auto ta_xxz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 58);

    auto ta_xxz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 59);

    auto ta_xxz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 60);

    auto ta_xxz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 61);

    auto ta_xxz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 62);

    auto ta_xxz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 63);

    auto ta_xxz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 64);

    auto ta_xxz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 65);

    auto ta_xxz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 66);

    auto ta_xxz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 67);

    auto ta_xxz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 68);

    auto ta_xxz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 69);

    auto ta_xxz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 70);

    auto ta_xxz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 71);

    auto ta_xxz_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 72);

    auto ta_xxz_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 73);

    auto ta_xxz_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 74);

    auto ta_xxz_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 75);

    auto ta_xxz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 76);

    auto ta_xxz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 78);

    auto ta_xxz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 79);

    auto ta_xxz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 80);

    auto ta_xxz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 81);

    auto ta_xxz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 82);

    auto ta_xxz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 83);

    auto ta_xyy_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 84);

    auto ta_xyy_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 85);

    auto ta_xyy_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 87);

    auto ta_xyy_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 88);

    auto ta_xyy_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 90);

    auto ta_xyy_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 91);

    auto ta_xyy_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 92);

    auto ta_xyy_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 94);

    auto ta_xyy_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 95);

    auto ta_xyy_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 96);

    auto ta_xyy_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 97);

    auto ta_xyy_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 99);

    auto ta_xyy_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 100);

    auto ta_xyy_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 101);

    auto ta_xyy_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 102);

    auto ta_xyy_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 103);

    auto ta_xyy_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 105);

    auto ta_xyy_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 106);

    auto ta_xyy_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 107);

    auto ta_xyy_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 108);

    auto ta_xyy_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 109);

    auto ta_xyy_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 110);

    auto ta_xyy_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 111);

    auto ta_xyz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 134);

    auto ta_xyz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 135);

    auto ta_xyz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 136);

    auto ta_xyz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 137);

    auto ta_xyz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 138);

    auto ta_xzz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 140);

    auto ta_xzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 142);

    auto ta_xzz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 144);

    auto ta_xzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 145);

    auto ta_xzz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 147);

    auto ta_xzz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 148);

    auto ta_xzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 149);

    auto ta_xzz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 151);

    auto ta_xzz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 152);

    auto ta_xzz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 153);

    auto ta_xzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 154);

    auto ta_xzz_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 156);

    auto ta_xzz_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 157);

    auto ta_xzz_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 158);

    auto ta_xzz_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 159);

    auto ta_xzz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 160);

    auto ta_xzz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 161);

    auto ta_xzz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 162);

    auto ta_xzz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 163);

    auto ta_xzz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 164);

    auto ta_xzz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 165);

    auto ta_xzz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 166);

    auto ta_xzz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 167);

    auto ta_yyy_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 168);

    auto ta_yyy_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 169);

    auto ta_yyy_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 170);

    auto ta_yyy_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 171);

    auto ta_yyy_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 172);

    auto ta_yyy_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 173);

    auto ta_yyy_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 174);

    auto ta_yyy_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 175);

    auto ta_yyy_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 176);

    auto ta_yyy_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 177);

    auto ta_yyy_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 178);

    auto ta_yyy_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 179);

    auto ta_yyy_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 180);

    auto ta_yyy_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 181);

    auto ta_yyy_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 182);

    auto ta_yyy_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 183);

    auto ta_yyy_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 184);

    auto ta_yyy_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 185);

    auto ta_yyy_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 186);

    auto ta_yyy_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 187);

    auto ta_yyy_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 188);

    auto ta_yyy_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 189);

    auto ta_yyy_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 190);

    auto ta_yyy_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 191);

    auto ta_yyy_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 192);

    auto ta_yyy_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 193);

    auto ta_yyy_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 194);

    auto ta_yyy_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 195);

    auto ta_yyz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 197);

    auto ta_yyz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 198);

    auto ta_yyz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 199);

    auto ta_yyz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 200);

    auto ta_yyz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 201);

    auto ta_yyz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 202);

    auto ta_yyz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 203);

    auto ta_yyz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 204);

    auto ta_yyz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 205);

    auto ta_yyz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 206);

    auto ta_yyz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 207);

    auto ta_yyz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 208);

    auto ta_yyz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 209);

    auto ta_yyz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 210);

    auto ta_yyz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 211);

    auto ta_yyz_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 212);

    auto ta_yyz_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 213);

    auto ta_yyz_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 214);

    auto ta_yyz_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 215);

    auto ta_yyz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 216);

    auto ta_yyz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 217);

    auto ta_yyz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 218);

    auto ta_yyz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 219);

    auto ta_yyz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 220);

    auto ta_yyz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 221);

    auto ta_yyz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 222);

    auto ta_yyz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 223);

    auto ta_yzz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 224);

    auto ta_yzz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 225);

    auto ta_yzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 226);

    auto ta_yzz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 227);

    auto ta_yzz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 228);

    auto ta_yzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 229);

    auto ta_yzz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 230);

    auto ta_yzz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 231);

    auto ta_yzz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 232);

    auto ta_yzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 233);

    auto ta_yzz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 234);

    auto ta_yzz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 235);

    auto ta_yzz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 236);

    auto ta_yzz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 237);

    auto ta_yzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 238);

    auto ta_yzz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 239);

    auto ta_yzz_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 240);

    auto ta_yzz_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 241);

    auto ta_yzz_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 242);

    auto ta_yzz_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 243);

    auto ta_yzz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 244);

    auto ta_yzz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 245);

    auto ta_yzz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 246);

    auto ta_yzz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 247);

    auto ta_yzz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 248);

    auto ta_yzz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 249);

    auto ta_yzz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 250);

    auto ta_yzz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 251);

    auto ta_zzz_xxxxxx_1 = pbuffer.data(idx_npot_1_fi + 252);

    auto ta_zzz_xxxxxy_1 = pbuffer.data(idx_npot_1_fi + 253);

    auto ta_zzz_xxxxxz_1 = pbuffer.data(idx_npot_1_fi + 254);

    auto ta_zzz_xxxxyy_1 = pbuffer.data(idx_npot_1_fi + 255);

    auto ta_zzz_xxxxyz_1 = pbuffer.data(idx_npot_1_fi + 256);

    auto ta_zzz_xxxxzz_1 = pbuffer.data(idx_npot_1_fi + 257);

    auto ta_zzz_xxxyyy_1 = pbuffer.data(idx_npot_1_fi + 258);

    auto ta_zzz_xxxyyz_1 = pbuffer.data(idx_npot_1_fi + 259);

    auto ta_zzz_xxxyzz_1 = pbuffer.data(idx_npot_1_fi + 260);

    auto ta_zzz_xxxzzz_1 = pbuffer.data(idx_npot_1_fi + 261);

    auto ta_zzz_xxyyyy_1 = pbuffer.data(idx_npot_1_fi + 262);

    auto ta_zzz_xxyyyz_1 = pbuffer.data(idx_npot_1_fi + 263);

    auto ta_zzz_xxyyzz_1 = pbuffer.data(idx_npot_1_fi + 264);

    auto ta_zzz_xxyzzz_1 = pbuffer.data(idx_npot_1_fi + 265);

    auto ta_zzz_xxzzzz_1 = pbuffer.data(idx_npot_1_fi + 266);

    auto ta_zzz_xyyyyy_1 = pbuffer.data(idx_npot_1_fi + 267);

    auto ta_zzz_xyyyyz_1 = pbuffer.data(idx_npot_1_fi + 268);

    auto ta_zzz_xyyyzz_1 = pbuffer.data(idx_npot_1_fi + 269);

    auto ta_zzz_xyyzzz_1 = pbuffer.data(idx_npot_1_fi + 270);

    auto ta_zzz_xyzzzz_1 = pbuffer.data(idx_npot_1_fi + 271);

    auto ta_zzz_xzzzzz_1 = pbuffer.data(idx_npot_1_fi + 272);

    auto ta_zzz_yyyyyy_1 = pbuffer.data(idx_npot_1_fi + 273);

    auto ta_zzz_yyyyyz_1 = pbuffer.data(idx_npot_1_fi + 274);

    auto ta_zzz_yyyyzz_1 = pbuffer.data(idx_npot_1_fi + 275);

    auto ta_zzz_yyyzzz_1 = pbuffer.data(idx_npot_1_fi + 276);

    auto ta_zzz_yyzzzz_1 = pbuffer.data(idx_npot_1_fi + 277);

    auto ta_zzz_yzzzzz_1 = pbuffer.data(idx_npot_1_fi + 278);

    auto ta_zzz_zzzzzz_1 = pbuffer.data(idx_npot_1_fi + 279);

    // Set up 0-28 components of targeted buffer : GI

    auto ta_xxxx_xxxxxx_0 = pbuffer.data(idx_npot_0_gi);

    auto ta_xxxx_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 1);

    auto ta_xxxx_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 2);

    auto ta_xxxx_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 3);

    auto ta_xxxx_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 4);

    auto ta_xxxx_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 5);

    auto ta_xxxx_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 6);

    auto ta_xxxx_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 7);

    auto ta_xxxx_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 8);

    auto ta_xxxx_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 9);

    auto ta_xxxx_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 10);

    auto ta_xxxx_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 11);

    auto ta_xxxx_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 12);

    auto ta_xxxx_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 13);

    auto ta_xxxx_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 14);

    auto ta_xxxx_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 15);

    auto ta_xxxx_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 16);

    auto ta_xxxx_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 17);

    auto ta_xxxx_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 18);

    auto ta_xxxx_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 19);

    auto ta_xxxx_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 20);

    auto ta_xxxx_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 21);

    auto ta_xxxx_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 22);

    auto ta_xxxx_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 23);

    auto ta_xxxx_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 24);

    auto ta_xxxx_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 25);

    auto ta_xxxx_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 26);

    auto ta_xxxx_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 27);

    #pragma omp simd aligned(pa_x, pc_x, ta_xx_xxxxxx_0, ta_xx_xxxxxx_1, ta_xx_xxxxxy_0, ta_xx_xxxxxy_1, ta_xx_xxxxxz_0, ta_xx_xxxxxz_1, ta_xx_xxxxyy_0, ta_xx_xxxxyy_1, ta_xx_xxxxyz_0, ta_xx_xxxxyz_1, ta_xx_xxxxzz_0, ta_xx_xxxxzz_1, ta_xx_xxxyyy_0, ta_xx_xxxyyy_1, ta_xx_xxxyyz_0, ta_xx_xxxyyz_1, ta_xx_xxxyzz_0, ta_xx_xxxyzz_1, ta_xx_xxxzzz_0, ta_xx_xxxzzz_1, ta_xx_xxyyyy_0, ta_xx_xxyyyy_1, ta_xx_xxyyyz_0, ta_xx_xxyyyz_1, ta_xx_xxyyzz_0, ta_xx_xxyyzz_1, ta_xx_xxyzzz_0, ta_xx_xxyzzz_1, ta_xx_xxzzzz_0, ta_xx_xxzzzz_1, ta_xx_xyyyyy_0, ta_xx_xyyyyy_1, ta_xx_xyyyyz_0, ta_xx_xyyyyz_1, ta_xx_xyyyzz_0, ta_xx_xyyyzz_1, ta_xx_xyyzzz_0, ta_xx_xyyzzz_1, ta_xx_xyzzzz_0, ta_xx_xyzzzz_1, ta_xx_xzzzzz_0, ta_xx_xzzzzz_1, ta_xx_yyyyyy_0, ta_xx_yyyyyy_1, ta_xx_yyyyyz_0, ta_xx_yyyyyz_1, ta_xx_yyyyzz_0, ta_xx_yyyyzz_1, ta_xx_yyyzzz_0, ta_xx_yyyzzz_1, ta_xx_yyzzzz_0, ta_xx_yyzzzz_1, ta_xx_yzzzzz_0, ta_xx_yzzzzz_1, ta_xx_zzzzzz_0, ta_xx_zzzzzz_1, ta_xxx_xxxxx_0, ta_xxx_xxxxx_1, ta_xxx_xxxxxx_0, ta_xxx_xxxxxx_1, ta_xxx_xxxxxy_0, ta_xxx_xxxxxy_1, ta_xxx_xxxxxz_0, ta_xxx_xxxxxz_1, ta_xxx_xxxxy_0, ta_xxx_xxxxy_1, ta_xxx_xxxxyy_0, ta_xxx_xxxxyy_1, ta_xxx_xxxxyz_0, ta_xxx_xxxxyz_1, ta_xxx_xxxxz_0, ta_xxx_xxxxz_1, ta_xxx_xxxxzz_0, ta_xxx_xxxxzz_1, ta_xxx_xxxyy_0, ta_xxx_xxxyy_1, ta_xxx_xxxyyy_0, ta_xxx_xxxyyy_1, ta_xxx_xxxyyz_0, ta_xxx_xxxyyz_1, ta_xxx_xxxyz_0, ta_xxx_xxxyz_1, ta_xxx_xxxyzz_0, ta_xxx_xxxyzz_1, ta_xxx_xxxzz_0, ta_xxx_xxxzz_1, ta_xxx_xxxzzz_0, ta_xxx_xxxzzz_1, ta_xxx_xxyyy_0, ta_xxx_xxyyy_1, ta_xxx_xxyyyy_0, ta_xxx_xxyyyy_1, ta_xxx_xxyyyz_0, ta_xxx_xxyyyz_1, ta_xxx_xxyyz_0, ta_xxx_xxyyz_1, ta_xxx_xxyyzz_0, ta_xxx_xxyyzz_1, ta_xxx_xxyzz_0, ta_xxx_xxyzz_1, ta_xxx_xxyzzz_0, ta_xxx_xxyzzz_1, ta_xxx_xxzzz_0, ta_xxx_xxzzz_1, ta_xxx_xxzzzz_0, ta_xxx_xxzzzz_1, ta_xxx_xyyyy_0, ta_xxx_xyyyy_1, ta_xxx_xyyyyy_0, ta_xxx_xyyyyy_1, ta_xxx_xyyyyz_0, ta_xxx_xyyyyz_1, ta_xxx_xyyyz_0, ta_xxx_xyyyz_1, ta_xxx_xyyyzz_0, ta_xxx_xyyyzz_1, ta_xxx_xyyzz_0, ta_xxx_xyyzz_1, ta_xxx_xyyzzz_0, ta_xxx_xyyzzz_1, ta_xxx_xyzzz_0, ta_xxx_xyzzz_1, ta_xxx_xyzzzz_0, ta_xxx_xyzzzz_1, ta_xxx_xzzzz_0, ta_xxx_xzzzz_1, ta_xxx_xzzzzz_0, ta_xxx_xzzzzz_1, ta_xxx_yyyyy_0, ta_xxx_yyyyy_1, ta_xxx_yyyyyy_0, ta_xxx_yyyyyy_1, ta_xxx_yyyyyz_0, ta_xxx_yyyyyz_1, ta_xxx_yyyyz_0, ta_xxx_yyyyz_1, ta_xxx_yyyyzz_0, ta_xxx_yyyyzz_1, ta_xxx_yyyzz_0, ta_xxx_yyyzz_1, ta_xxx_yyyzzz_0, ta_xxx_yyyzzz_1, ta_xxx_yyzzz_0, ta_xxx_yyzzz_1, ta_xxx_yyzzzz_0, ta_xxx_yyzzzz_1, ta_xxx_yzzzz_0, ta_xxx_yzzzz_1, ta_xxx_yzzzzz_0, ta_xxx_yzzzzz_1, ta_xxx_zzzzz_0, ta_xxx_zzzzz_1, ta_xxx_zzzzzz_0, ta_xxx_zzzzzz_1, ta_xxxx_xxxxxx_0, ta_xxxx_xxxxxy_0, ta_xxxx_xxxxxz_0, ta_xxxx_xxxxyy_0, ta_xxxx_xxxxyz_0, ta_xxxx_xxxxzz_0, ta_xxxx_xxxyyy_0, ta_xxxx_xxxyyz_0, ta_xxxx_xxxyzz_0, ta_xxxx_xxxzzz_0, ta_xxxx_xxyyyy_0, ta_xxxx_xxyyyz_0, ta_xxxx_xxyyzz_0, ta_xxxx_xxyzzz_0, ta_xxxx_xxzzzz_0, ta_xxxx_xyyyyy_0, ta_xxxx_xyyyyz_0, ta_xxxx_xyyyzz_0, ta_xxxx_xyyzzz_0, ta_xxxx_xyzzzz_0, ta_xxxx_xzzzzz_0, ta_xxxx_yyyyyy_0, ta_xxxx_yyyyyz_0, ta_xxxx_yyyyzz_0, ta_xxxx_yyyzzz_0, ta_xxxx_yyzzzz_0, ta_xxxx_yzzzzz_0, ta_xxxx_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_xxxxxx_0[i] = 3.0 * ta_xx_xxxxxx_0[i] * fe_0 - 3.0 * ta_xx_xxxxxx_1[i] * fe_0 + 6.0 * ta_xxx_xxxxx_0[i] * fe_0 - 6.0 * ta_xxx_xxxxx_1[i] * fe_0 + ta_xxx_xxxxxx_0[i] * pa_x[i] - ta_xxx_xxxxxx_1[i] * pc_x[i];

        ta_xxxx_xxxxxy_0[i] = 3.0 * ta_xx_xxxxxy_0[i] * fe_0 - 3.0 * ta_xx_xxxxxy_1[i] * fe_0 + 5.0 * ta_xxx_xxxxy_0[i] * fe_0 - 5.0 * ta_xxx_xxxxy_1[i] * fe_0 + ta_xxx_xxxxxy_0[i] * pa_x[i] - ta_xxx_xxxxxy_1[i] * pc_x[i];

        ta_xxxx_xxxxxz_0[i] = 3.0 * ta_xx_xxxxxz_0[i] * fe_0 - 3.0 * ta_xx_xxxxxz_1[i] * fe_0 + 5.0 * ta_xxx_xxxxz_0[i] * fe_0 - 5.0 * ta_xxx_xxxxz_1[i] * fe_0 + ta_xxx_xxxxxz_0[i] * pa_x[i] - ta_xxx_xxxxxz_1[i] * pc_x[i];

        ta_xxxx_xxxxyy_0[i] = 3.0 * ta_xx_xxxxyy_0[i] * fe_0 - 3.0 * ta_xx_xxxxyy_1[i] * fe_0 + 4.0 * ta_xxx_xxxyy_0[i] * fe_0 - 4.0 * ta_xxx_xxxyy_1[i] * fe_0 + ta_xxx_xxxxyy_0[i] * pa_x[i] - ta_xxx_xxxxyy_1[i] * pc_x[i];

        ta_xxxx_xxxxyz_0[i] = 3.0 * ta_xx_xxxxyz_0[i] * fe_0 - 3.0 * ta_xx_xxxxyz_1[i] * fe_0 + 4.0 * ta_xxx_xxxyz_0[i] * fe_0 - 4.0 * ta_xxx_xxxyz_1[i] * fe_0 + ta_xxx_xxxxyz_0[i] * pa_x[i] - ta_xxx_xxxxyz_1[i] * pc_x[i];

        ta_xxxx_xxxxzz_0[i] = 3.0 * ta_xx_xxxxzz_0[i] * fe_0 - 3.0 * ta_xx_xxxxzz_1[i] * fe_0 + 4.0 * ta_xxx_xxxzz_0[i] * fe_0 - 4.0 * ta_xxx_xxxzz_1[i] * fe_0 + ta_xxx_xxxxzz_0[i] * pa_x[i] - ta_xxx_xxxxzz_1[i] * pc_x[i];

        ta_xxxx_xxxyyy_0[i] = 3.0 * ta_xx_xxxyyy_0[i] * fe_0 - 3.0 * ta_xx_xxxyyy_1[i] * fe_0 + 3.0 * ta_xxx_xxyyy_0[i] * fe_0 - 3.0 * ta_xxx_xxyyy_1[i] * fe_0 + ta_xxx_xxxyyy_0[i] * pa_x[i] - ta_xxx_xxxyyy_1[i] * pc_x[i];

        ta_xxxx_xxxyyz_0[i] = 3.0 * ta_xx_xxxyyz_0[i] * fe_0 - 3.0 * ta_xx_xxxyyz_1[i] * fe_0 + 3.0 * ta_xxx_xxyyz_0[i] * fe_0 - 3.0 * ta_xxx_xxyyz_1[i] * fe_0 + ta_xxx_xxxyyz_0[i] * pa_x[i] - ta_xxx_xxxyyz_1[i] * pc_x[i];

        ta_xxxx_xxxyzz_0[i] = 3.0 * ta_xx_xxxyzz_0[i] * fe_0 - 3.0 * ta_xx_xxxyzz_1[i] * fe_0 + 3.0 * ta_xxx_xxyzz_0[i] * fe_0 - 3.0 * ta_xxx_xxyzz_1[i] * fe_0 + ta_xxx_xxxyzz_0[i] * pa_x[i] - ta_xxx_xxxyzz_1[i] * pc_x[i];

        ta_xxxx_xxxzzz_0[i] = 3.0 * ta_xx_xxxzzz_0[i] * fe_0 - 3.0 * ta_xx_xxxzzz_1[i] * fe_0 + 3.0 * ta_xxx_xxzzz_0[i] * fe_0 - 3.0 * ta_xxx_xxzzz_1[i] * fe_0 + ta_xxx_xxxzzz_0[i] * pa_x[i] - ta_xxx_xxxzzz_1[i] * pc_x[i];

        ta_xxxx_xxyyyy_0[i] = 3.0 * ta_xx_xxyyyy_0[i] * fe_0 - 3.0 * ta_xx_xxyyyy_1[i] * fe_0 + 2.0 * ta_xxx_xyyyy_0[i] * fe_0 - 2.0 * ta_xxx_xyyyy_1[i] * fe_0 + ta_xxx_xxyyyy_0[i] * pa_x[i] - ta_xxx_xxyyyy_1[i] * pc_x[i];

        ta_xxxx_xxyyyz_0[i] = 3.0 * ta_xx_xxyyyz_0[i] * fe_0 - 3.0 * ta_xx_xxyyyz_1[i] * fe_0 + 2.0 * ta_xxx_xyyyz_0[i] * fe_0 - 2.0 * ta_xxx_xyyyz_1[i] * fe_0 + ta_xxx_xxyyyz_0[i] * pa_x[i] - ta_xxx_xxyyyz_1[i] * pc_x[i];

        ta_xxxx_xxyyzz_0[i] = 3.0 * ta_xx_xxyyzz_0[i] * fe_0 - 3.0 * ta_xx_xxyyzz_1[i] * fe_0 + 2.0 * ta_xxx_xyyzz_0[i] * fe_0 - 2.0 * ta_xxx_xyyzz_1[i] * fe_0 + ta_xxx_xxyyzz_0[i] * pa_x[i] - ta_xxx_xxyyzz_1[i] * pc_x[i];

        ta_xxxx_xxyzzz_0[i] = 3.0 * ta_xx_xxyzzz_0[i] * fe_0 - 3.0 * ta_xx_xxyzzz_1[i] * fe_0 + 2.0 * ta_xxx_xyzzz_0[i] * fe_0 - 2.0 * ta_xxx_xyzzz_1[i] * fe_0 + ta_xxx_xxyzzz_0[i] * pa_x[i] - ta_xxx_xxyzzz_1[i] * pc_x[i];

        ta_xxxx_xxzzzz_0[i] = 3.0 * ta_xx_xxzzzz_0[i] * fe_0 - 3.0 * ta_xx_xxzzzz_1[i] * fe_0 + 2.0 * ta_xxx_xzzzz_0[i] * fe_0 - 2.0 * ta_xxx_xzzzz_1[i] * fe_0 + ta_xxx_xxzzzz_0[i] * pa_x[i] - ta_xxx_xxzzzz_1[i] * pc_x[i];

        ta_xxxx_xyyyyy_0[i] = 3.0 * ta_xx_xyyyyy_0[i] * fe_0 - 3.0 * ta_xx_xyyyyy_1[i] * fe_0 + ta_xxx_yyyyy_0[i] * fe_0 - ta_xxx_yyyyy_1[i] * fe_0 + ta_xxx_xyyyyy_0[i] * pa_x[i] - ta_xxx_xyyyyy_1[i] * pc_x[i];

        ta_xxxx_xyyyyz_0[i] = 3.0 * ta_xx_xyyyyz_0[i] * fe_0 - 3.0 * ta_xx_xyyyyz_1[i] * fe_0 + ta_xxx_yyyyz_0[i] * fe_0 - ta_xxx_yyyyz_1[i] * fe_0 + ta_xxx_xyyyyz_0[i] * pa_x[i] - ta_xxx_xyyyyz_1[i] * pc_x[i];

        ta_xxxx_xyyyzz_0[i] = 3.0 * ta_xx_xyyyzz_0[i] * fe_0 - 3.0 * ta_xx_xyyyzz_1[i] * fe_0 + ta_xxx_yyyzz_0[i] * fe_0 - ta_xxx_yyyzz_1[i] * fe_0 + ta_xxx_xyyyzz_0[i] * pa_x[i] - ta_xxx_xyyyzz_1[i] * pc_x[i];

        ta_xxxx_xyyzzz_0[i] = 3.0 * ta_xx_xyyzzz_0[i] * fe_0 - 3.0 * ta_xx_xyyzzz_1[i] * fe_0 + ta_xxx_yyzzz_0[i] * fe_0 - ta_xxx_yyzzz_1[i] * fe_0 + ta_xxx_xyyzzz_0[i] * pa_x[i] - ta_xxx_xyyzzz_1[i] * pc_x[i];

        ta_xxxx_xyzzzz_0[i] = 3.0 * ta_xx_xyzzzz_0[i] * fe_0 - 3.0 * ta_xx_xyzzzz_1[i] * fe_0 + ta_xxx_yzzzz_0[i] * fe_0 - ta_xxx_yzzzz_1[i] * fe_0 + ta_xxx_xyzzzz_0[i] * pa_x[i] - ta_xxx_xyzzzz_1[i] * pc_x[i];

        ta_xxxx_xzzzzz_0[i] = 3.0 * ta_xx_xzzzzz_0[i] * fe_0 - 3.0 * ta_xx_xzzzzz_1[i] * fe_0 + ta_xxx_zzzzz_0[i] * fe_0 - ta_xxx_zzzzz_1[i] * fe_0 + ta_xxx_xzzzzz_0[i] * pa_x[i] - ta_xxx_xzzzzz_1[i] * pc_x[i];

        ta_xxxx_yyyyyy_0[i] = 3.0 * ta_xx_yyyyyy_0[i] * fe_0 - 3.0 * ta_xx_yyyyyy_1[i] * fe_0 + ta_xxx_yyyyyy_0[i] * pa_x[i] - ta_xxx_yyyyyy_1[i] * pc_x[i];

        ta_xxxx_yyyyyz_0[i] = 3.0 * ta_xx_yyyyyz_0[i] * fe_0 - 3.0 * ta_xx_yyyyyz_1[i] * fe_0 + ta_xxx_yyyyyz_0[i] * pa_x[i] - ta_xxx_yyyyyz_1[i] * pc_x[i];

        ta_xxxx_yyyyzz_0[i] = 3.0 * ta_xx_yyyyzz_0[i] * fe_0 - 3.0 * ta_xx_yyyyzz_1[i] * fe_0 + ta_xxx_yyyyzz_0[i] * pa_x[i] - ta_xxx_yyyyzz_1[i] * pc_x[i];

        ta_xxxx_yyyzzz_0[i] = 3.0 * ta_xx_yyyzzz_0[i] * fe_0 - 3.0 * ta_xx_yyyzzz_1[i] * fe_0 + ta_xxx_yyyzzz_0[i] * pa_x[i] - ta_xxx_yyyzzz_1[i] * pc_x[i];

        ta_xxxx_yyzzzz_0[i] = 3.0 * ta_xx_yyzzzz_0[i] * fe_0 - 3.0 * ta_xx_yyzzzz_1[i] * fe_0 + ta_xxx_yyzzzz_0[i] * pa_x[i] - ta_xxx_yyzzzz_1[i] * pc_x[i];

        ta_xxxx_yzzzzz_0[i] = 3.0 * ta_xx_yzzzzz_0[i] * fe_0 - 3.0 * ta_xx_yzzzzz_1[i] * fe_0 + ta_xxx_yzzzzz_0[i] * pa_x[i] - ta_xxx_yzzzzz_1[i] * pc_x[i];

        ta_xxxx_zzzzzz_0[i] = 3.0 * ta_xx_zzzzzz_0[i] * fe_0 - 3.0 * ta_xx_zzzzzz_1[i] * fe_0 + ta_xxx_zzzzzz_0[i] * pa_x[i] - ta_xxx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : GI

    auto ta_xxxy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 28);

    auto ta_xxxy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 29);

    auto ta_xxxy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 30);

    auto ta_xxxy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 31);

    auto ta_xxxy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 32);

    auto ta_xxxy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 33);

    auto ta_xxxy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 34);

    auto ta_xxxy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 35);

    auto ta_xxxy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 36);

    auto ta_xxxy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 37);

    auto ta_xxxy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 38);

    auto ta_xxxy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 39);

    auto ta_xxxy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 40);

    auto ta_xxxy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 41);

    auto ta_xxxy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 42);

    auto ta_xxxy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 43);

    auto ta_xxxy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 44);

    auto ta_xxxy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 45);

    auto ta_xxxy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 46);

    auto ta_xxxy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 47);

    auto ta_xxxy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 48);

    auto ta_xxxy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 49);

    auto ta_xxxy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 50);

    auto ta_xxxy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 51);

    auto ta_xxxy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 52);

    auto ta_xxxy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 53);

    auto ta_xxxy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 54);

    auto ta_xxxy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 55);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xxx_xxxxx_0, ta_xxx_xxxxx_1, ta_xxx_xxxxxx_0, ta_xxx_xxxxxx_1, ta_xxx_xxxxxy_0, ta_xxx_xxxxxy_1, ta_xxx_xxxxxz_0, ta_xxx_xxxxxz_1, ta_xxx_xxxxy_0, ta_xxx_xxxxy_1, ta_xxx_xxxxyy_0, ta_xxx_xxxxyy_1, ta_xxx_xxxxyz_0, ta_xxx_xxxxyz_1, ta_xxx_xxxxz_0, ta_xxx_xxxxz_1, ta_xxx_xxxxzz_0, ta_xxx_xxxxzz_1, ta_xxx_xxxyy_0, ta_xxx_xxxyy_1, ta_xxx_xxxyyy_0, ta_xxx_xxxyyy_1, ta_xxx_xxxyyz_0, ta_xxx_xxxyyz_1, ta_xxx_xxxyz_0, ta_xxx_xxxyz_1, ta_xxx_xxxyzz_0, ta_xxx_xxxyzz_1, ta_xxx_xxxzz_0, ta_xxx_xxxzz_1, ta_xxx_xxxzzz_0, ta_xxx_xxxzzz_1, ta_xxx_xxyyy_0, ta_xxx_xxyyy_1, ta_xxx_xxyyyy_0, ta_xxx_xxyyyy_1, ta_xxx_xxyyyz_0, ta_xxx_xxyyyz_1, ta_xxx_xxyyz_0, ta_xxx_xxyyz_1, ta_xxx_xxyyzz_0, ta_xxx_xxyyzz_1, ta_xxx_xxyzz_0, ta_xxx_xxyzz_1, ta_xxx_xxyzzz_0, ta_xxx_xxyzzz_1, ta_xxx_xxzzz_0, ta_xxx_xxzzz_1, ta_xxx_xxzzzz_0, ta_xxx_xxzzzz_1, ta_xxx_xyyyy_0, ta_xxx_xyyyy_1, ta_xxx_xyyyyy_0, ta_xxx_xyyyyy_1, ta_xxx_xyyyyz_0, ta_xxx_xyyyyz_1, ta_xxx_xyyyz_0, ta_xxx_xyyyz_1, ta_xxx_xyyyzz_0, ta_xxx_xyyyzz_1, ta_xxx_xyyzz_0, ta_xxx_xyyzz_1, ta_xxx_xyyzzz_0, ta_xxx_xyyzzz_1, ta_xxx_xyzzz_0, ta_xxx_xyzzz_1, ta_xxx_xyzzzz_0, ta_xxx_xyzzzz_1, ta_xxx_xzzzz_0, ta_xxx_xzzzz_1, ta_xxx_xzzzzz_0, ta_xxx_xzzzzz_1, ta_xxx_zzzzzz_0, ta_xxx_zzzzzz_1, ta_xxxy_xxxxxx_0, ta_xxxy_xxxxxy_0, ta_xxxy_xxxxxz_0, ta_xxxy_xxxxyy_0, ta_xxxy_xxxxyz_0, ta_xxxy_xxxxzz_0, ta_xxxy_xxxyyy_0, ta_xxxy_xxxyyz_0, ta_xxxy_xxxyzz_0, ta_xxxy_xxxzzz_0, ta_xxxy_xxyyyy_0, ta_xxxy_xxyyyz_0, ta_xxxy_xxyyzz_0, ta_xxxy_xxyzzz_0, ta_xxxy_xxzzzz_0, ta_xxxy_xyyyyy_0, ta_xxxy_xyyyyz_0, ta_xxxy_xyyyzz_0, ta_xxxy_xyyzzz_0, ta_xxxy_xyzzzz_0, ta_xxxy_xzzzzz_0, ta_xxxy_yyyyyy_0, ta_xxxy_yyyyyz_0, ta_xxxy_yyyyzz_0, ta_xxxy_yyyzzz_0, ta_xxxy_yyzzzz_0, ta_xxxy_yzzzzz_0, ta_xxxy_zzzzzz_0, ta_xxy_yyyyyy_0, ta_xxy_yyyyyy_1, ta_xxy_yyyyyz_0, ta_xxy_yyyyyz_1, ta_xxy_yyyyzz_0, ta_xxy_yyyyzz_1, ta_xxy_yyyzzz_0, ta_xxy_yyyzzz_1, ta_xxy_yyzzzz_0, ta_xxy_yyzzzz_1, ta_xxy_yzzzzz_0, ta_xxy_yzzzzz_1, ta_xy_yyyyyy_0, ta_xy_yyyyyy_1, ta_xy_yyyyyz_0, ta_xy_yyyyyz_1, ta_xy_yyyyzz_0, ta_xy_yyyyzz_1, ta_xy_yyyzzz_0, ta_xy_yyyzzz_1, ta_xy_yyzzzz_0, ta_xy_yyzzzz_1, ta_xy_yzzzzz_0, ta_xy_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxy_xxxxxx_0[i] = ta_xxx_xxxxxx_0[i] * pa_y[i] - ta_xxx_xxxxxx_1[i] * pc_y[i];

        ta_xxxy_xxxxxy_0[i] = ta_xxx_xxxxx_0[i] * fe_0 - ta_xxx_xxxxx_1[i] * fe_0 + ta_xxx_xxxxxy_0[i] * pa_y[i] - ta_xxx_xxxxxy_1[i] * pc_y[i];

        ta_xxxy_xxxxxz_0[i] = ta_xxx_xxxxxz_0[i] * pa_y[i] - ta_xxx_xxxxxz_1[i] * pc_y[i];

        ta_xxxy_xxxxyy_0[i] = 2.0 * ta_xxx_xxxxy_0[i] * fe_0 - 2.0 * ta_xxx_xxxxy_1[i] * fe_0 + ta_xxx_xxxxyy_0[i] * pa_y[i] - ta_xxx_xxxxyy_1[i] * pc_y[i];

        ta_xxxy_xxxxyz_0[i] = ta_xxx_xxxxz_0[i] * fe_0 - ta_xxx_xxxxz_1[i] * fe_0 + ta_xxx_xxxxyz_0[i] * pa_y[i] - ta_xxx_xxxxyz_1[i] * pc_y[i];

        ta_xxxy_xxxxzz_0[i] = ta_xxx_xxxxzz_0[i] * pa_y[i] - ta_xxx_xxxxzz_1[i] * pc_y[i];

        ta_xxxy_xxxyyy_0[i] = 3.0 * ta_xxx_xxxyy_0[i] * fe_0 - 3.0 * ta_xxx_xxxyy_1[i] * fe_0 + ta_xxx_xxxyyy_0[i] * pa_y[i] - ta_xxx_xxxyyy_1[i] * pc_y[i];

        ta_xxxy_xxxyyz_0[i] = 2.0 * ta_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta_xxx_xxxyz_1[i] * fe_0 + ta_xxx_xxxyyz_0[i] * pa_y[i] - ta_xxx_xxxyyz_1[i] * pc_y[i];

        ta_xxxy_xxxyzz_0[i] = ta_xxx_xxxzz_0[i] * fe_0 - ta_xxx_xxxzz_1[i] * fe_0 + ta_xxx_xxxyzz_0[i] * pa_y[i] - ta_xxx_xxxyzz_1[i] * pc_y[i];

        ta_xxxy_xxxzzz_0[i] = ta_xxx_xxxzzz_0[i] * pa_y[i] - ta_xxx_xxxzzz_1[i] * pc_y[i];

        ta_xxxy_xxyyyy_0[i] = 4.0 * ta_xxx_xxyyy_0[i] * fe_0 - 4.0 * ta_xxx_xxyyy_1[i] * fe_0 + ta_xxx_xxyyyy_0[i] * pa_y[i] - ta_xxx_xxyyyy_1[i] * pc_y[i];

        ta_xxxy_xxyyyz_0[i] = 3.0 * ta_xxx_xxyyz_0[i] * fe_0 - 3.0 * ta_xxx_xxyyz_1[i] * fe_0 + ta_xxx_xxyyyz_0[i] * pa_y[i] - ta_xxx_xxyyyz_1[i] * pc_y[i];

        ta_xxxy_xxyyzz_0[i] = 2.0 * ta_xxx_xxyzz_0[i] * fe_0 - 2.0 * ta_xxx_xxyzz_1[i] * fe_0 + ta_xxx_xxyyzz_0[i] * pa_y[i] - ta_xxx_xxyyzz_1[i] * pc_y[i];

        ta_xxxy_xxyzzz_0[i] = ta_xxx_xxzzz_0[i] * fe_0 - ta_xxx_xxzzz_1[i] * fe_0 + ta_xxx_xxyzzz_0[i] * pa_y[i] - ta_xxx_xxyzzz_1[i] * pc_y[i];

        ta_xxxy_xxzzzz_0[i] = ta_xxx_xxzzzz_0[i] * pa_y[i] - ta_xxx_xxzzzz_1[i] * pc_y[i];

        ta_xxxy_xyyyyy_0[i] = 5.0 * ta_xxx_xyyyy_0[i] * fe_0 - 5.0 * ta_xxx_xyyyy_1[i] * fe_0 + ta_xxx_xyyyyy_0[i] * pa_y[i] - ta_xxx_xyyyyy_1[i] * pc_y[i];

        ta_xxxy_xyyyyz_0[i] = 4.0 * ta_xxx_xyyyz_0[i] * fe_0 - 4.0 * ta_xxx_xyyyz_1[i] * fe_0 + ta_xxx_xyyyyz_0[i] * pa_y[i] - ta_xxx_xyyyyz_1[i] * pc_y[i];

        ta_xxxy_xyyyzz_0[i] = 3.0 * ta_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta_xxx_xyyzz_1[i] * fe_0 + ta_xxx_xyyyzz_0[i] * pa_y[i] - ta_xxx_xyyyzz_1[i] * pc_y[i];

        ta_xxxy_xyyzzz_0[i] = 2.0 * ta_xxx_xyzzz_0[i] * fe_0 - 2.0 * ta_xxx_xyzzz_1[i] * fe_0 + ta_xxx_xyyzzz_0[i] * pa_y[i] - ta_xxx_xyyzzz_1[i] * pc_y[i];

        ta_xxxy_xyzzzz_0[i] = ta_xxx_xzzzz_0[i] * fe_0 - ta_xxx_xzzzz_1[i] * fe_0 + ta_xxx_xyzzzz_0[i] * pa_y[i] - ta_xxx_xyzzzz_1[i] * pc_y[i];

        ta_xxxy_xzzzzz_0[i] = ta_xxx_xzzzzz_0[i] * pa_y[i] - ta_xxx_xzzzzz_1[i] * pc_y[i];

        ta_xxxy_yyyyyy_0[i] = 2.0 * ta_xy_yyyyyy_0[i] * fe_0 - 2.0 * ta_xy_yyyyyy_1[i] * fe_0 + ta_xxy_yyyyyy_0[i] * pa_x[i] - ta_xxy_yyyyyy_1[i] * pc_x[i];

        ta_xxxy_yyyyyz_0[i] = 2.0 * ta_xy_yyyyyz_0[i] * fe_0 - 2.0 * ta_xy_yyyyyz_1[i] * fe_0 + ta_xxy_yyyyyz_0[i] * pa_x[i] - ta_xxy_yyyyyz_1[i] * pc_x[i];

        ta_xxxy_yyyyzz_0[i] = 2.0 * ta_xy_yyyyzz_0[i] * fe_0 - 2.0 * ta_xy_yyyyzz_1[i] * fe_0 + ta_xxy_yyyyzz_0[i] * pa_x[i] - ta_xxy_yyyyzz_1[i] * pc_x[i];

        ta_xxxy_yyyzzz_0[i] = 2.0 * ta_xy_yyyzzz_0[i] * fe_0 - 2.0 * ta_xy_yyyzzz_1[i] * fe_0 + ta_xxy_yyyzzz_0[i] * pa_x[i] - ta_xxy_yyyzzz_1[i] * pc_x[i];

        ta_xxxy_yyzzzz_0[i] = 2.0 * ta_xy_yyzzzz_0[i] * fe_0 - 2.0 * ta_xy_yyzzzz_1[i] * fe_0 + ta_xxy_yyzzzz_0[i] * pa_x[i] - ta_xxy_yyzzzz_1[i] * pc_x[i];

        ta_xxxy_yzzzzz_0[i] = 2.0 * ta_xy_yzzzzz_0[i] * fe_0 - 2.0 * ta_xy_yzzzzz_1[i] * fe_0 + ta_xxy_yzzzzz_0[i] * pa_x[i] - ta_xxy_yzzzzz_1[i] * pc_x[i];

        ta_xxxy_zzzzzz_0[i] = ta_xxx_zzzzzz_0[i] * pa_y[i] - ta_xxx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : GI

    auto ta_xxxz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 56);

    auto ta_xxxz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 57);

    auto ta_xxxz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 58);

    auto ta_xxxz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 59);

    auto ta_xxxz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 60);

    auto ta_xxxz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 61);

    auto ta_xxxz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 62);

    auto ta_xxxz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 63);

    auto ta_xxxz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 64);

    auto ta_xxxz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 65);

    auto ta_xxxz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 66);

    auto ta_xxxz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 67);

    auto ta_xxxz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 68);

    auto ta_xxxz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 69);

    auto ta_xxxz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 70);

    auto ta_xxxz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 71);

    auto ta_xxxz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 72);

    auto ta_xxxz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 73);

    auto ta_xxxz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 74);

    auto ta_xxxz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 75);

    auto ta_xxxz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 76);

    auto ta_xxxz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 77);

    auto ta_xxxz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 78);

    auto ta_xxxz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 79);

    auto ta_xxxz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 80);

    auto ta_xxxz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 81);

    auto ta_xxxz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 82);

    auto ta_xxxz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 83);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xxx_xxxxx_0, ta_xxx_xxxxx_1, ta_xxx_xxxxxx_0, ta_xxx_xxxxxx_1, ta_xxx_xxxxxy_0, ta_xxx_xxxxxy_1, ta_xxx_xxxxxz_0, ta_xxx_xxxxxz_1, ta_xxx_xxxxy_0, ta_xxx_xxxxy_1, ta_xxx_xxxxyy_0, ta_xxx_xxxxyy_1, ta_xxx_xxxxyz_0, ta_xxx_xxxxyz_1, ta_xxx_xxxxz_0, ta_xxx_xxxxz_1, ta_xxx_xxxxzz_0, ta_xxx_xxxxzz_1, ta_xxx_xxxyy_0, ta_xxx_xxxyy_1, ta_xxx_xxxyyy_0, ta_xxx_xxxyyy_1, ta_xxx_xxxyyz_0, ta_xxx_xxxyyz_1, ta_xxx_xxxyz_0, ta_xxx_xxxyz_1, ta_xxx_xxxyzz_0, ta_xxx_xxxyzz_1, ta_xxx_xxxzz_0, ta_xxx_xxxzz_1, ta_xxx_xxxzzz_0, ta_xxx_xxxzzz_1, ta_xxx_xxyyy_0, ta_xxx_xxyyy_1, ta_xxx_xxyyyy_0, ta_xxx_xxyyyy_1, ta_xxx_xxyyyz_0, ta_xxx_xxyyyz_1, ta_xxx_xxyyz_0, ta_xxx_xxyyz_1, ta_xxx_xxyyzz_0, ta_xxx_xxyyzz_1, ta_xxx_xxyzz_0, ta_xxx_xxyzz_1, ta_xxx_xxyzzz_0, ta_xxx_xxyzzz_1, ta_xxx_xxzzz_0, ta_xxx_xxzzz_1, ta_xxx_xxzzzz_0, ta_xxx_xxzzzz_1, ta_xxx_xyyyy_0, ta_xxx_xyyyy_1, ta_xxx_xyyyyy_0, ta_xxx_xyyyyy_1, ta_xxx_xyyyyz_0, ta_xxx_xyyyyz_1, ta_xxx_xyyyz_0, ta_xxx_xyyyz_1, ta_xxx_xyyyzz_0, ta_xxx_xyyyzz_1, ta_xxx_xyyzz_0, ta_xxx_xyyzz_1, ta_xxx_xyyzzz_0, ta_xxx_xyyzzz_1, ta_xxx_xyzzz_0, ta_xxx_xyzzz_1, ta_xxx_xyzzzz_0, ta_xxx_xyzzzz_1, ta_xxx_xzzzz_0, ta_xxx_xzzzz_1, ta_xxx_xzzzzz_0, ta_xxx_xzzzzz_1, ta_xxx_yyyyyy_0, ta_xxx_yyyyyy_1, ta_xxxz_xxxxxx_0, ta_xxxz_xxxxxy_0, ta_xxxz_xxxxxz_0, ta_xxxz_xxxxyy_0, ta_xxxz_xxxxyz_0, ta_xxxz_xxxxzz_0, ta_xxxz_xxxyyy_0, ta_xxxz_xxxyyz_0, ta_xxxz_xxxyzz_0, ta_xxxz_xxxzzz_0, ta_xxxz_xxyyyy_0, ta_xxxz_xxyyyz_0, ta_xxxz_xxyyzz_0, ta_xxxz_xxyzzz_0, ta_xxxz_xxzzzz_0, ta_xxxz_xyyyyy_0, ta_xxxz_xyyyyz_0, ta_xxxz_xyyyzz_0, ta_xxxz_xyyzzz_0, ta_xxxz_xyzzzz_0, ta_xxxz_xzzzzz_0, ta_xxxz_yyyyyy_0, ta_xxxz_yyyyyz_0, ta_xxxz_yyyyzz_0, ta_xxxz_yyyzzz_0, ta_xxxz_yyzzzz_0, ta_xxxz_yzzzzz_0, ta_xxxz_zzzzzz_0, ta_xxz_yyyyyz_0, ta_xxz_yyyyyz_1, ta_xxz_yyyyzz_0, ta_xxz_yyyyzz_1, ta_xxz_yyyzzz_0, ta_xxz_yyyzzz_1, ta_xxz_yyzzzz_0, ta_xxz_yyzzzz_1, ta_xxz_yzzzzz_0, ta_xxz_yzzzzz_1, ta_xxz_zzzzzz_0, ta_xxz_zzzzzz_1, ta_xz_yyyyyz_0, ta_xz_yyyyyz_1, ta_xz_yyyyzz_0, ta_xz_yyyyzz_1, ta_xz_yyyzzz_0, ta_xz_yyyzzz_1, ta_xz_yyzzzz_0, ta_xz_yyzzzz_1, ta_xz_yzzzzz_0, ta_xz_yzzzzz_1, ta_xz_zzzzzz_0, ta_xz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxz_xxxxxx_0[i] = ta_xxx_xxxxxx_0[i] * pa_z[i] - ta_xxx_xxxxxx_1[i] * pc_z[i];

        ta_xxxz_xxxxxy_0[i] = ta_xxx_xxxxxy_0[i] * pa_z[i] - ta_xxx_xxxxxy_1[i] * pc_z[i];

        ta_xxxz_xxxxxz_0[i] = ta_xxx_xxxxx_0[i] * fe_0 - ta_xxx_xxxxx_1[i] * fe_0 + ta_xxx_xxxxxz_0[i] * pa_z[i] - ta_xxx_xxxxxz_1[i] * pc_z[i];

        ta_xxxz_xxxxyy_0[i] = ta_xxx_xxxxyy_0[i] * pa_z[i] - ta_xxx_xxxxyy_1[i] * pc_z[i];

        ta_xxxz_xxxxyz_0[i] = ta_xxx_xxxxy_0[i] * fe_0 - ta_xxx_xxxxy_1[i] * fe_0 + ta_xxx_xxxxyz_0[i] * pa_z[i] - ta_xxx_xxxxyz_1[i] * pc_z[i];

        ta_xxxz_xxxxzz_0[i] = 2.0 * ta_xxx_xxxxz_0[i] * fe_0 - 2.0 * ta_xxx_xxxxz_1[i] * fe_0 + ta_xxx_xxxxzz_0[i] * pa_z[i] - ta_xxx_xxxxzz_1[i] * pc_z[i];

        ta_xxxz_xxxyyy_0[i] = ta_xxx_xxxyyy_0[i] * pa_z[i] - ta_xxx_xxxyyy_1[i] * pc_z[i];

        ta_xxxz_xxxyyz_0[i] = ta_xxx_xxxyy_0[i] * fe_0 - ta_xxx_xxxyy_1[i] * fe_0 + ta_xxx_xxxyyz_0[i] * pa_z[i] - ta_xxx_xxxyyz_1[i] * pc_z[i];

        ta_xxxz_xxxyzz_0[i] = 2.0 * ta_xxx_xxxyz_0[i] * fe_0 - 2.0 * ta_xxx_xxxyz_1[i] * fe_0 + ta_xxx_xxxyzz_0[i] * pa_z[i] - ta_xxx_xxxyzz_1[i] * pc_z[i];

        ta_xxxz_xxxzzz_0[i] = 3.0 * ta_xxx_xxxzz_0[i] * fe_0 - 3.0 * ta_xxx_xxxzz_1[i] * fe_0 + ta_xxx_xxxzzz_0[i] * pa_z[i] - ta_xxx_xxxzzz_1[i] * pc_z[i];

        ta_xxxz_xxyyyy_0[i] = ta_xxx_xxyyyy_0[i] * pa_z[i] - ta_xxx_xxyyyy_1[i] * pc_z[i];

        ta_xxxz_xxyyyz_0[i] = ta_xxx_xxyyy_0[i] * fe_0 - ta_xxx_xxyyy_1[i] * fe_0 + ta_xxx_xxyyyz_0[i] * pa_z[i] - ta_xxx_xxyyyz_1[i] * pc_z[i];

        ta_xxxz_xxyyzz_0[i] = 2.0 * ta_xxx_xxyyz_0[i] * fe_0 - 2.0 * ta_xxx_xxyyz_1[i] * fe_0 + ta_xxx_xxyyzz_0[i] * pa_z[i] - ta_xxx_xxyyzz_1[i] * pc_z[i];

        ta_xxxz_xxyzzz_0[i] = 3.0 * ta_xxx_xxyzz_0[i] * fe_0 - 3.0 * ta_xxx_xxyzz_1[i] * fe_0 + ta_xxx_xxyzzz_0[i] * pa_z[i] - ta_xxx_xxyzzz_1[i] * pc_z[i];

        ta_xxxz_xxzzzz_0[i] = 4.0 * ta_xxx_xxzzz_0[i] * fe_0 - 4.0 * ta_xxx_xxzzz_1[i] * fe_0 + ta_xxx_xxzzzz_0[i] * pa_z[i] - ta_xxx_xxzzzz_1[i] * pc_z[i];

        ta_xxxz_xyyyyy_0[i] = ta_xxx_xyyyyy_0[i] * pa_z[i] - ta_xxx_xyyyyy_1[i] * pc_z[i];

        ta_xxxz_xyyyyz_0[i] = ta_xxx_xyyyy_0[i] * fe_0 - ta_xxx_xyyyy_1[i] * fe_0 + ta_xxx_xyyyyz_0[i] * pa_z[i] - ta_xxx_xyyyyz_1[i] * pc_z[i];

        ta_xxxz_xyyyzz_0[i] = 2.0 * ta_xxx_xyyyz_0[i] * fe_0 - 2.0 * ta_xxx_xyyyz_1[i] * fe_0 + ta_xxx_xyyyzz_0[i] * pa_z[i] - ta_xxx_xyyyzz_1[i] * pc_z[i];

        ta_xxxz_xyyzzz_0[i] = 3.0 * ta_xxx_xyyzz_0[i] * fe_0 - 3.0 * ta_xxx_xyyzz_1[i] * fe_0 + ta_xxx_xyyzzz_0[i] * pa_z[i] - ta_xxx_xyyzzz_1[i] * pc_z[i];

        ta_xxxz_xyzzzz_0[i] = 4.0 * ta_xxx_xyzzz_0[i] * fe_0 - 4.0 * ta_xxx_xyzzz_1[i] * fe_0 + ta_xxx_xyzzzz_0[i] * pa_z[i] - ta_xxx_xyzzzz_1[i] * pc_z[i];

        ta_xxxz_xzzzzz_0[i] = 5.0 * ta_xxx_xzzzz_0[i] * fe_0 - 5.0 * ta_xxx_xzzzz_1[i] * fe_0 + ta_xxx_xzzzzz_0[i] * pa_z[i] - ta_xxx_xzzzzz_1[i] * pc_z[i];

        ta_xxxz_yyyyyy_0[i] = ta_xxx_yyyyyy_0[i] * pa_z[i] - ta_xxx_yyyyyy_1[i] * pc_z[i];

        ta_xxxz_yyyyyz_0[i] = 2.0 * ta_xz_yyyyyz_0[i] * fe_0 - 2.0 * ta_xz_yyyyyz_1[i] * fe_0 + ta_xxz_yyyyyz_0[i] * pa_x[i] - ta_xxz_yyyyyz_1[i] * pc_x[i];

        ta_xxxz_yyyyzz_0[i] = 2.0 * ta_xz_yyyyzz_0[i] * fe_0 - 2.0 * ta_xz_yyyyzz_1[i] * fe_0 + ta_xxz_yyyyzz_0[i] * pa_x[i] - ta_xxz_yyyyzz_1[i] * pc_x[i];

        ta_xxxz_yyyzzz_0[i] = 2.0 * ta_xz_yyyzzz_0[i] * fe_0 - 2.0 * ta_xz_yyyzzz_1[i] * fe_0 + ta_xxz_yyyzzz_0[i] * pa_x[i] - ta_xxz_yyyzzz_1[i] * pc_x[i];

        ta_xxxz_yyzzzz_0[i] = 2.0 * ta_xz_yyzzzz_0[i] * fe_0 - 2.0 * ta_xz_yyzzzz_1[i] * fe_0 + ta_xxz_yyzzzz_0[i] * pa_x[i] - ta_xxz_yyzzzz_1[i] * pc_x[i];

        ta_xxxz_yzzzzz_0[i] = 2.0 * ta_xz_yzzzzz_0[i] * fe_0 - 2.0 * ta_xz_yzzzzz_1[i] * fe_0 + ta_xxz_yzzzzz_0[i] * pa_x[i] - ta_xxz_yzzzzz_1[i] * pc_x[i];

        ta_xxxz_zzzzzz_0[i] = 2.0 * ta_xz_zzzzzz_0[i] * fe_0 - 2.0 * ta_xz_zzzzzz_1[i] * fe_0 + ta_xxz_zzzzzz_0[i] * pa_x[i] - ta_xxz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 84-112 components of targeted buffer : GI

    auto ta_xxyy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 84);

    auto ta_xxyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 85);

    auto ta_xxyy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 86);

    auto ta_xxyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 87);

    auto ta_xxyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 88);

    auto ta_xxyy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 89);

    auto ta_xxyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 90);

    auto ta_xxyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 91);

    auto ta_xxyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 92);

    auto ta_xxyy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 93);

    auto ta_xxyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 94);

    auto ta_xxyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 95);

    auto ta_xxyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 96);

    auto ta_xxyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 97);

    auto ta_xxyy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 98);

    auto ta_xxyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 99);

    auto ta_xxyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 100);

    auto ta_xxyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 101);

    auto ta_xxyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 102);

    auto ta_xxyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 103);

    auto ta_xxyy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 104);

    auto ta_xxyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 105);

    auto ta_xxyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 106);

    auto ta_xxyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 107);

    auto ta_xxyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 108);

    auto ta_xxyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 109);

    auto ta_xxyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 110);

    auto ta_xxyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 111);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xx_xxxxxx_0, ta_xx_xxxxxx_1, ta_xx_xxxxxz_0, ta_xx_xxxxxz_1, ta_xx_xxxxzz_0, ta_xx_xxxxzz_1, ta_xx_xxxzzz_0, ta_xx_xxxzzz_1, ta_xx_xxzzzz_0, ta_xx_xxzzzz_1, ta_xx_xzzzzz_0, ta_xx_xzzzzz_1, ta_xxy_xxxxxx_0, ta_xxy_xxxxxx_1, ta_xxy_xxxxxz_0, ta_xxy_xxxxxz_1, ta_xxy_xxxxzz_0, ta_xxy_xxxxzz_1, ta_xxy_xxxzzz_0, ta_xxy_xxxzzz_1, ta_xxy_xxzzzz_0, ta_xxy_xxzzzz_1, ta_xxy_xzzzzz_0, ta_xxy_xzzzzz_1, ta_xxyy_xxxxxx_0, ta_xxyy_xxxxxy_0, ta_xxyy_xxxxxz_0, ta_xxyy_xxxxyy_0, ta_xxyy_xxxxyz_0, ta_xxyy_xxxxzz_0, ta_xxyy_xxxyyy_0, ta_xxyy_xxxyyz_0, ta_xxyy_xxxyzz_0, ta_xxyy_xxxzzz_0, ta_xxyy_xxyyyy_0, ta_xxyy_xxyyyz_0, ta_xxyy_xxyyzz_0, ta_xxyy_xxyzzz_0, ta_xxyy_xxzzzz_0, ta_xxyy_xyyyyy_0, ta_xxyy_xyyyyz_0, ta_xxyy_xyyyzz_0, ta_xxyy_xyyzzz_0, ta_xxyy_xyzzzz_0, ta_xxyy_xzzzzz_0, ta_xxyy_yyyyyy_0, ta_xxyy_yyyyyz_0, ta_xxyy_yyyyzz_0, ta_xxyy_yyyzzz_0, ta_xxyy_yyzzzz_0, ta_xxyy_yzzzzz_0, ta_xxyy_zzzzzz_0, ta_xyy_xxxxxy_0, ta_xyy_xxxxxy_1, ta_xyy_xxxxy_0, ta_xyy_xxxxy_1, ta_xyy_xxxxyy_0, ta_xyy_xxxxyy_1, ta_xyy_xxxxyz_0, ta_xyy_xxxxyz_1, ta_xyy_xxxyy_0, ta_xyy_xxxyy_1, ta_xyy_xxxyyy_0, ta_xyy_xxxyyy_1, ta_xyy_xxxyyz_0, ta_xyy_xxxyyz_1, ta_xyy_xxxyz_0, ta_xyy_xxxyz_1, ta_xyy_xxxyzz_0, ta_xyy_xxxyzz_1, ta_xyy_xxyyy_0, ta_xyy_xxyyy_1, ta_xyy_xxyyyy_0, ta_xyy_xxyyyy_1, ta_xyy_xxyyyz_0, ta_xyy_xxyyyz_1, ta_xyy_xxyyz_0, ta_xyy_xxyyz_1, ta_xyy_xxyyzz_0, ta_xyy_xxyyzz_1, ta_xyy_xxyzz_0, ta_xyy_xxyzz_1, ta_xyy_xxyzzz_0, ta_xyy_xxyzzz_1, ta_xyy_xyyyy_0, ta_xyy_xyyyy_1, ta_xyy_xyyyyy_0, ta_xyy_xyyyyy_1, ta_xyy_xyyyyz_0, ta_xyy_xyyyyz_1, ta_xyy_xyyyz_0, ta_xyy_xyyyz_1, ta_xyy_xyyyzz_0, ta_xyy_xyyyzz_1, ta_xyy_xyyzz_0, ta_xyy_xyyzz_1, ta_xyy_xyyzzz_0, ta_xyy_xyyzzz_1, ta_xyy_xyzzz_0, ta_xyy_xyzzz_1, ta_xyy_xyzzzz_0, ta_xyy_xyzzzz_1, ta_xyy_yyyyy_0, ta_xyy_yyyyy_1, ta_xyy_yyyyyy_0, ta_xyy_yyyyyy_1, ta_xyy_yyyyyz_0, ta_xyy_yyyyyz_1, ta_xyy_yyyyz_0, ta_xyy_yyyyz_1, ta_xyy_yyyyzz_0, ta_xyy_yyyyzz_1, ta_xyy_yyyzz_0, ta_xyy_yyyzz_1, ta_xyy_yyyzzz_0, ta_xyy_yyyzzz_1, ta_xyy_yyzzz_0, ta_xyy_yyzzz_1, ta_xyy_yyzzzz_0, ta_xyy_yyzzzz_1, ta_xyy_yzzzz_0, ta_xyy_yzzzz_1, ta_xyy_yzzzzz_0, ta_xyy_yzzzzz_1, ta_xyy_zzzzzz_0, ta_xyy_zzzzzz_1, ta_yy_xxxxxy_0, ta_yy_xxxxxy_1, ta_yy_xxxxyy_0, ta_yy_xxxxyy_1, ta_yy_xxxxyz_0, ta_yy_xxxxyz_1, ta_yy_xxxyyy_0, ta_yy_xxxyyy_1, ta_yy_xxxyyz_0, ta_yy_xxxyyz_1, ta_yy_xxxyzz_0, ta_yy_xxxyzz_1, ta_yy_xxyyyy_0, ta_yy_xxyyyy_1, ta_yy_xxyyyz_0, ta_yy_xxyyyz_1, ta_yy_xxyyzz_0, ta_yy_xxyyzz_1, ta_yy_xxyzzz_0, ta_yy_xxyzzz_1, ta_yy_xyyyyy_0, ta_yy_xyyyyy_1, ta_yy_xyyyyz_0, ta_yy_xyyyyz_1, ta_yy_xyyyzz_0, ta_yy_xyyyzz_1, ta_yy_xyyzzz_0, ta_yy_xyyzzz_1, ta_yy_xyzzzz_0, ta_yy_xyzzzz_1, ta_yy_yyyyyy_0, ta_yy_yyyyyy_1, ta_yy_yyyyyz_0, ta_yy_yyyyyz_1, ta_yy_yyyyzz_0, ta_yy_yyyyzz_1, ta_yy_yyyzzz_0, ta_yy_yyyzzz_1, ta_yy_yyzzzz_0, ta_yy_yyzzzz_1, ta_yy_yzzzzz_0, ta_yy_yzzzzz_1, ta_yy_zzzzzz_0, ta_yy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyy_xxxxxx_0[i] = ta_xx_xxxxxx_0[i] * fe_0 - ta_xx_xxxxxx_1[i] * fe_0 + ta_xxy_xxxxxx_0[i] * pa_y[i] - ta_xxy_xxxxxx_1[i] * pc_y[i];

        ta_xxyy_xxxxxy_0[i] = ta_yy_xxxxxy_0[i] * fe_0 - ta_yy_xxxxxy_1[i] * fe_0 + 5.0 * ta_xyy_xxxxy_0[i] * fe_0 - 5.0 * ta_xyy_xxxxy_1[i] * fe_0 + ta_xyy_xxxxxy_0[i] * pa_x[i] - ta_xyy_xxxxxy_1[i] * pc_x[i];

        ta_xxyy_xxxxxz_0[i] = ta_xx_xxxxxz_0[i] * fe_0 - ta_xx_xxxxxz_1[i] * fe_0 + ta_xxy_xxxxxz_0[i] * pa_y[i] - ta_xxy_xxxxxz_1[i] * pc_y[i];

        ta_xxyy_xxxxyy_0[i] = ta_yy_xxxxyy_0[i] * fe_0 - ta_yy_xxxxyy_1[i] * fe_0 + 4.0 * ta_xyy_xxxyy_0[i] * fe_0 - 4.0 * ta_xyy_xxxyy_1[i] * fe_0 + ta_xyy_xxxxyy_0[i] * pa_x[i] - ta_xyy_xxxxyy_1[i] * pc_x[i];

        ta_xxyy_xxxxyz_0[i] = ta_yy_xxxxyz_0[i] * fe_0 - ta_yy_xxxxyz_1[i] * fe_0 + 4.0 * ta_xyy_xxxyz_0[i] * fe_0 - 4.0 * ta_xyy_xxxyz_1[i] * fe_0 + ta_xyy_xxxxyz_0[i] * pa_x[i] - ta_xyy_xxxxyz_1[i] * pc_x[i];

        ta_xxyy_xxxxzz_0[i] = ta_xx_xxxxzz_0[i] * fe_0 - ta_xx_xxxxzz_1[i] * fe_0 + ta_xxy_xxxxzz_0[i] * pa_y[i] - ta_xxy_xxxxzz_1[i] * pc_y[i];

        ta_xxyy_xxxyyy_0[i] = ta_yy_xxxyyy_0[i] * fe_0 - ta_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta_xyy_xxyyy_0[i] * fe_0 - 3.0 * ta_xyy_xxyyy_1[i] * fe_0 + ta_xyy_xxxyyy_0[i] * pa_x[i] - ta_xyy_xxxyyy_1[i] * pc_x[i];

        ta_xxyy_xxxyyz_0[i] = ta_yy_xxxyyz_0[i] * fe_0 - ta_yy_xxxyyz_1[i] * fe_0 + 3.0 * ta_xyy_xxyyz_0[i] * fe_0 - 3.0 * ta_xyy_xxyyz_1[i] * fe_0 + ta_xyy_xxxyyz_0[i] * pa_x[i] - ta_xyy_xxxyyz_1[i] * pc_x[i];

        ta_xxyy_xxxyzz_0[i] = ta_yy_xxxyzz_0[i] * fe_0 - ta_yy_xxxyzz_1[i] * fe_0 + 3.0 * ta_xyy_xxyzz_0[i] * fe_0 - 3.0 * ta_xyy_xxyzz_1[i] * fe_0 + ta_xyy_xxxyzz_0[i] * pa_x[i] - ta_xyy_xxxyzz_1[i] * pc_x[i];

        ta_xxyy_xxxzzz_0[i] = ta_xx_xxxzzz_0[i] * fe_0 - ta_xx_xxxzzz_1[i] * fe_0 + ta_xxy_xxxzzz_0[i] * pa_y[i] - ta_xxy_xxxzzz_1[i] * pc_y[i];

        ta_xxyy_xxyyyy_0[i] = ta_yy_xxyyyy_0[i] * fe_0 - ta_yy_xxyyyy_1[i] * fe_0 + 2.0 * ta_xyy_xyyyy_0[i] * fe_0 - 2.0 * ta_xyy_xyyyy_1[i] * fe_0 + ta_xyy_xxyyyy_0[i] * pa_x[i] - ta_xyy_xxyyyy_1[i] * pc_x[i];

        ta_xxyy_xxyyyz_0[i] = ta_yy_xxyyyz_0[i] * fe_0 - ta_yy_xxyyyz_1[i] * fe_0 + 2.0 * ta_xyy_xyyyz_0[i] * fe_0 - 2.0 * ta_xyy_xyyyz_1[i] * fe_0 + ta_xyy_xxyyyz_0[i] * pa_x[i] - ta_xyy_xxyyyz_1[i] * pc_x[i];

        ta_xxyy_xxyyzz_0[i] = ta_yy_xxyyzz_0[i] * fe_0 - ta_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta_xyy_xyyzz_0[i] * fe_0 - 2.0 * ta_xyy_xyyzz_1[i] * fe_0 + ta_xyy_xxyyzz_0[i] * pa_x[i] - ta_xyy_xxyyzz_1[i] * pc_x[i];

        ta_xxyy_xxyzzz_0[i] = ta_yy_xxyzzz_0[i] * fe_0 - ta_yy_xxyzzz_1[i] * fe_0 + 2.0 * ta_xyy_xyzzz_0[i] * fe_0 - 2.0 * ta_xyy_xyzzz_1[i] * fe_0 + ta_xyy_xxyzzz_0[i] * pa_x[i] - ta_xyy_xxyzzz_1[i] * pc_x[i];

        ta_xxyy_xxzzzz_0[i] = ta_xx_xxzzzz_0[i] * fe_0 - ta_xx_xxzzzz_1[i] * fe_0 + ta_xxy_xxzzzz_0[i] * pa_y[i] - ta_xxy_xxzzzz_1[i] * pc_y[i];

        ta_xxyy_xyyyyy_0[i] = ta_yy_xyyyyy_0[i] * fe_0 - ta_yy_xyyyyy_1[i] * fe_0 + ta_xyy_yyyyy_0[i] * fe_0 - ta_xyy_yyyyy_1[i] * fe_0 + ta_xyy_xyyyyy_0[i] * pa_x[i] - ta_xyy_xyyyyy_1[i] * pc_x[i];

        ta_xxyy_xyyyyz_0[i] = ta_yy_xyyyyz_0[i] * fe_0 - ta_yy_xyyyyz_1[i] * fe_0 + ta_xyy_yyyyz_0[i] * fe_0 - ta_xyy_yyyyz_1[i] * fe_0 + ta_xyy_xyyyyz_0[i] * pa_x[i] - ta_xyy_xyyyyz_1[i] * pc_x[i];

        ta_xxyy_xyyyzz_0[i] = ta_yy_xyyyzz_0[i] * fe_0 - ta_yy_xyyyzz_1[i] * fe_0 + ta_xyy_yyyzz_0[i] * fe_0 - ta_xyy_yyyzz_1[i] * fe_0 + ta_xyy_xyyyzz_0[i] * pa_x[i] - ta_xyy_xyyyzz_1[i] * pc_x[i];

        ta_xxyy_xyyzzz_0[i] = ta_yy_xyyzzz_0[i] * fe_0 - ta_yy_xyyzzz_1[i] * fe_0 + ta_xyy_yyzzz_0[i] * fe_0 - ta_xyy_yyzzz_1[i] * fe_0 + ta_xyy_xyyzzz_0[i] * pa_x[i] - ta_xyy_xyyzzz_1[i] * pc_x[i];

        ta_xxyy_xyzzzz_0[i] = ta_yy_xyzzzz_0[i] * fe_0 - ta_yy_xyzzzz_1[i] * fe_0 + ta_xyy_yzzzz_0[i] * fe_0 - ta_xyy_yzzzz_1[i] * fe_0 + ta_xyy_xyzzzz_0[i] * pa_x[i] - ta_xyy_xyzzzz_1[i] * pc_x[i];

        ta_xxyy_xzzzzz_0[i] = ta_xx_xzzzzz_0[i] * fe_0 - ta_xx_xzzzzz_1[i] * fe_0 + ta_xxy_xzzzzz_0[i] * pa_y[i] - ta_xxy_xzzzzz_1[i] * pc_y[i];

        ta_xxyy_yyyyyy_0[i] = ta_yy_yyyyyy_0[i] * fe_0 - ta_yy_yyyyyy_1[i] * fe_0 + ta_xyy_yyyyyy_0[i] * pa_x[i] - ta_xyy_yyyyyy_1[i] * pc_x[i];

        ta_xxyy_yyyyyz_0[i] = ta_yy_yyyyyz_0[i] * fe_0 - ta_yy_yyyyyz_1[i] * fe_0 + ta_xyy_yyyyyz_0[i] * pa_x[i] - ta_xyy_yyyyyz_1[i] * pc_x[i];

        ta_xxyy_yyyyzz_0[i] = ta_yy_yyyyzz_0[i] * fe_0 - ta_yy_yyyyzz_1[i] * fe_0 + ta_xyy_yyyyzz_0[i] * pa_x[i] - ta_xyy_yyyyzz_1[i] * pc_x[i];

        ta_xxyy_yyyzzz_0[i] = ta_yy_yyyzzz_0[i] * fe_0 - ta_yy_yyyzzz_1[i] * fe_0 + ta_xyy_yyyzzz_0[i] * pa_x[i] - ta_xyy_yyyzzz_1[i] * pc_x[i];

        ta_xxyy_yyzzzz_0[i] = ta_yy_yyzzzz_0[i] * fe_0 - ta_yy_yyzzzz_1[i] * fe_0 + ta_xyy_yyzzzz_0[i] * pa_x[i] - ta_xyy_yyzzzz_1[i] * pc_x[i];

        ta_xxyy_yzzzzz_0[i] = ta_yy_yzzzzz_0[i] * fe_0 - ta_yy_yzzzzz_1[i] * fe_0 + ta_xyy_yzzzzz_0[i] * pa_x[i] - ta_xyy_yzzzzz_1[i] * pc_x[i];

        ta_xxyy_zzzzzz_0[i] = ta_yy_zzzzzz_0[i] * fe_0 - ta_yy_zzzzzz_1[i] * fe_0 + ta_xyy_zzzzzz_0[i] * pa_x[i] - ta_xyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 112-140 components of targeted buffer : GI

    auto ta_xxyz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 112);

    auto ta_xxyz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 113);

    auto ta_xxyz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 114);

    auto ta_xxyz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 115);

    auto ta_xxyz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 116);

    auto ta_xxyz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 117);

    auto ta_xxyz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 118);

    auto ta_xxyz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 119);

    auto ta_xxyz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 120);

    auto ta_xxyz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 121);

    auto ta_xxyz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 122);

    auto ta_xxyz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 123);

    auto ta_xxyz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 124);

    auto ta_xxyz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 125);

    auto ta_xxyz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 126);

    auto ta_xxyz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 127);

    auto ta_xxyz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 128);

    auto ta_xxyz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 129);

    auto ta_xxyz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 130);

    auto ta_xxyz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 131);

    auto ta_xxyz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 132);

    auto ta_xxyz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 133);

    auto ta_xxyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 134);

    auto ta_xxyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 135);

    auto ta_xxyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 136);

    auto ta_xxyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 137);

    auto ta_xxyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 138);

    auto ta_xxyz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 139);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxy_xxxxxy_0, ta_xxy_xxxxxy_1, ta_xxy_xxxxyy_0, ta_xxy_xxxxyy_1, ta_xxy_xxxyyy_0, ta_xxy_xxxyyy_1, ta_xxy_xxyyyy_0, ta_xxy_xxyyyy_1, ta_xxy_xyyyyy_0, ta_xxy_xyyyyy_1, ta_xxy_yyyyyy_0, ta_xxy_yyyyyy_1, ta_xxyz_xxxxxx_0, ta_xxyz_xxxxxy_0, ta_xxyz_xxxxxz_0, ta_xxyz_xxxxyy_0, ta_xxyz_xxxxyz_0, ta_xxyz_xxxxzz_0, ta_xxyz_xxxyyy_0, ta_xxyz_xxxyyz_0, ta_xxyz_xxxyzz_0, ta_xxyz_xxxzzz_0, ta_xxyz_xxyyyy_0, ta_xxyz_xxyyyz_0, ta_xxyz_xxyyzz_0, ta_xxyz_xxyzzz_0, ta_xxyz_xxzzzz_0, ta_xxyz_xyyyyy_0, ta_xxyz_xyyyyz_0, ta_xxyz_xyyyzz_0, ta_xxyz_xyyzzz_0, ta_xxyz_xyzzzz_0, ta_xxyz_xzzzzz_0, ta_xxyz_yyyyyy_0, ta_xxyz_yyyyyz_0, ta_xxyz_yyyyzz_0, ta_xxyz_yyyzzz_0, ta_xxyz_yyzzzz_0, ta_xxyz_yzzzzz_0, ta_xxyz_zzzzzz_0, ta_xxz_xxxxxx_0, ta_xxz_xxxxxx_1, ta_xxz_xxxxxz_0, ta_xxz_xxxxxz_1, ta_xxz_xxxxyz_0, ta_xxz_xxxxyz_1, ta_xxz_xxxxz_0, ta_xxz_xxxxz_1, ta_xxz_xxxxzz_0, ta_xxz_xxxxzz_1, ta_xxz_xxxyyz_0, ta_xxz_xxxyyz_1, ta_xxz_xxxyz_0, ta_xxz_xxxyz_1, ta_xxz_xxxyzz_0, ta_xxz_xxxyzz_1, ta_xxz_xxxzz_0, ta_xxz_xxxzz_1, ta_xxz_xxxzzz_0, ta_xxz_xxxzzz_1, ta_xxz_xxyyyz_0, ta_xxz_xxyyyz_1, ta_xxz_xxyyz_0, ta_xxz_xxyyz_1, ta_xxz_xxyyzz_0, ta_xxz_xxyyzz_1, ta_xxz_xxyzz_0, ta_xxz_xxyzz_1, ta_xxz_xxyzzz_0, ta_xxz_xxyzzz_1, ta_xxz_xxzzz_0, ta_xxz_xxzzz_1, ta_xxz_xxzzzz_0, ta_xxz_xxzzzz_1, ta_xxz_xyyyyz_0, ta_xxz_xyyyyz_1, ta_xxz_xyyyz_0, ta_xxz_xyyyz_1, ta_xxz_xyyyzz_0, ta_xxz_xyyyzz_1, ta_xxz_xyyzz_0, ta_xxz_xyyzz_1, ta_xxz_xyyzzz_0, ta_xxz_xyyzzz_1, ta_xxz_xyzzz_0, ta_xxz_xyzzz_1, ta_xxz_xyzzzz_0, ta_xxz_xyzzzz_1, ta_xxz_xzzzz_0, ta_xxz_xzzzz_1, ta_xxz_xzzzzz_0, ta_xxz_xzzzzz_1, ta_xxz_zzzzzz_0, ta_xxz_zzzzzz_1, ta_xyz_yyyyyz_0, ta_xyz_yyyyyz_1, ta_xyz_yyyyzz_0, ta_xyz_yyyyzz_1, ta_xyz_yyyzzz_0, ta_xyz_yyyzzz_1, ta_xyz_yyzzzz_0, ta_xyz_yyzzzz_1, ta_xyz_yzzzzz_0, ta_xyz_yzzzzz_1, ta_yz_yyyyyz_0, ta_yz_yyyyyz_1, ta_yz_yyyyzz_0, ta_yz_yyyyzz_1, ta_yz_yyyzzz_0, ta_yz_yyyzzz_1, ta_yz_yyzzzz_0, ta_yz_yyzzzz_1, ta_yz_yzzzzz_0, ta_yz_yzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyz_xxxxxx_0[i] = ta_xxz_xxxxxx_0[i] * pa_y[i] - ta_xxz_xxxxxx_1[i] * pc_y[i];

        ta_xxyz_xxxxxy_0[i] = ta_xxy_xxxxxy_0[i] * pa_z[i] - ta_xxy_xxxxxy_1[i] * pc_z[i];

        ta_xxyz_xxxxxz_0[i] = ta_xxz_xxxxxz_0[i] * pa_y[i] - ta_xxz_xxxxxz_1[i] * pc_y[i];

        ta_xxyz_xxxxyy_0[i] = ta_xxy_xxxxyy_0[i] * pa_z[i] - ta_xxy_xxxxyy_1[i] * pc_z[i];

        ta_xxyz_xxxxyz_0[i] = ta_xxz_xxxxz_0[i] * fe_0 - ta_xxz_xxxxz_1[i] * fe_0 + ta_xxz_xxxxyz_0[i] * pa_y[i] - ta_xxz_xxxxyz_1[i] * pc_y[i];

        ta_xxyz_xxxxzz_0[i] = ta_xxz_xxxxzz_0[i] * pa_y[i] - ta_xxz_xxxxzz_1[i] * pc_y[i];

        ta_xxyz_xxxyyy_0[i] = ta_xxy_xxxyyy_0[i] * pa_z[i] - ta_xxy_xxxyyy_1[i] * pc_z[i];

        ta_xxyz_xxxyyz_0[i] = 2.0 * ta_xxz_xxxyz_0[i] * fe_0 - 2.0 * ta_xxz_xxxyz_1[i] * fe_0 + ta_xxz_xxxyyz_0[i] * pa_y[i] - ta_xxz_xxxyyz_1[i] * pc_y[i];

        ta_xxyz_xxxyzz_0[i] = ta_xxz_xxxzz_0[i] * fe_0 - ta_xxz_xxxzz_1[i] * fe_0 + ta_xxz_xxxyzz_0[i] * pa_y[i] - ta_xxz_xxxyzz_1[i] * pc_y[i];

        ta_xxyz_xxxzzz_0[i] = ta_xxz_xxxzzz_0[i] * pa_y[i] - ta_xxz_xxxzzz_1[i] * pc_y[i];

        ta_xxyz_xxyyyy_0[i] = ta_xxy_xxyyyy_0[i] * pa_z[i] - ta_xxy_xxyyyy_1[i] * pc_z[i];

        ta_xxyz_xxyyyz_0[i] = 3.0 * ta_xxz_xxyyz_0[i] * fe_0 - 3.0 * ta_xxz_xxyyz_1[i] * fe_0 + ta_xxz_xxyyyz_0[i] * pa_y[i] - ta_xxz_xxyyyz_1[i] * pc_y[i];

        ta_xxyz_xxyyzz_0[i] = 2.0 * ta_xxz_xxyzz_0[i] * fe_0 - 2.0 * ta_xxz_xxyzz_1[i] * fe_0 + ta_xxz_xxyyzz_0[i] * pa_y[i] - ta_xxz_xxyyzz_1[i] * pc_y[i];

        ta_xxyz_xxyzzz_0[i] = ta_xxz_xxzzz_0[i] * fe_0 - ta_xxz_xxzzz_1[i] * fe_0 + ta_xxz_xxyzzz_0[i] * pa_y[i] - ta_xxz_xxyzzz_1[i] * pc_y[i];

        ta_xxyz_xxzzzz_0[i] = ta_xxz_xxzzzz_0[i] * pa_y[i] - ta_xxz_xxzzzz_1[i] * pc_y[i];

        ta_xxyz_xyyyyy_0[i] = ta_xxy_xyyyyy_0[i] * pa_z[i] - ta_xxy_xyyyyy_1[i] * pc_z[i];

        ta_xxyz_xyyyyz_0[i] = 4.0 * ta_xxz_xyyyz_0[i] * fe_0 - 4.0 * ta_xxz_xyyyz_1[i] * fe_0 + ta_xxz_xyyyyz_0[i] * pa_y[i] - ta_xxz_xyyyyz_1[i] * pc_y[i];

        ta_xxyz_xyyyzz_0[i] = 3.0 * ta_xxz_xyyzz_0[i] * fe_0 - 3.0 * ta_xxz_xyyzz_1[i] * fe_0 + ta_xxz_xyyyzz_0[i] * pa_y[i] - ta_xxz_xyyyzz_1[i] * pc_y[i];

        ta_xxyz_xyyzzz_0[i] = 2.0 * ta_xxz_xyzzz_0[i] * fe_0 - 2.0 * ta_xxz_xyzzz_1[i] * fe_0 + ta_xxz_xyyzzz_0[i] * pa_y[i] - ta_xxz_xyyzzz_1[i] * pc_y[i];

        ta_xxyz_xyzzzz_0[i] = ta_xxz_xzzzz_0[i] * fe_0 - ta_xxz_xzzzz_1[i] * fe_0 + ta_xxz_xyzzzz_0[i] * pa_y[i] - ta_xxz_xyzzzz_1[i] * pc_y[i];

        ta_xxyz_xzzzzz_0[i] = ta_xxz_xzzzzz_0[i] * pa_y[i] - ta_xxz_xzzzzz_1[i] * pc_y[i];

        ta_xxyz_yyyyyy_0[i] = ta_xxy_yyyyyy_0[i] * pa_z[i] - ta_xxy_yyyyyy_1[i] * pc_z[i];

        ta_xxyz_yyyyyz_0[i] = ta_yz_yyyyyz_0[i] * fe_0 - ta_yz_yyyyyz_1[i] * fe_0 + ta_xyz_yyyyyz_0[i] * pa_x[i] - ta_xyz_yyyyyz_1[i] * pc_x[i];

        ta_xxyz_yyyyzz_0[i] = ta_yz_yyyyzz_0[i] * fe_0 - ta_yz_yyyyzz_1[i] * fe_0 + ta_xyz_yyyyzz_0[i] * pa_x[i] - ta_xyz_yyyyzz_1[i] * pc_x[i];

        ta_xxyz_yyyzzz_0[i] = ta_yz_yyyzzz_0[i] * fe_0 - ta_yz_yyyzzz_1[i] * fe_0 + ta_xyz_yyyzzz_0[i] * pa_x[i] - ta_xyz_yyyzzz_1[i] * pc_x[i];

        ta_xxyz_yyzzzz_0[i] = ta_yz_yyzzzz_0[i] * fe_0 - ta_yz_yyzzzz_1[i] * fe_0 + ta_xyz_yyzzzz_0[i] * pa_x[i] - ta_xyz_yyzzzz_1[i] * pc_x[i];

        ta_xxyz_yzzzzz_0[i] = ta_yz_yzzzzz_0[i] * fe_0 - ta_yz_yzzzzz_1[i] * fe_0 + ta_xyz_yzzzzz_0[i] * pa_x[i] - ta_xyz_yzzzzz_1[i] * pc_x[i];

        ta_xxyz_zzzzzz_0[i] = ta_xxz_zzzzzz_0[i] * pa_y[i] - ta_xxz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : GI

    auto ta_xxzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 140);

    auto ta_xxzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 141);

    auto ta_xxzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 142);

    auto ta_xxzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 143);

    auto ta_xxzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 144);

    auto ta_xxzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 145);

    auto ta_xxzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 146);

    auto ta_xxzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 147);

    auto ta_xxzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 148);

    auto ta_xxzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 149);

    auto ta_xxzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 150);

    auto ta_xxzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 151);

    auto ta_xxzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 152);

    auto ta_xxzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 153);

    auto ta_xxzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 154);

    auto ta_xxzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 155);

    auto ta_xxzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 156);

    auto ta_xxzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 157);

    auto ta_xxzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 158);

    auto ta_xxzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 159);

    auto ta_xxzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 160);

    auto ta_xxzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 161);

    auto ta_xxzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 162);

    auto ta_xxzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 163);

    auto ta_xxzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 164);

    auto ta_xxzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 165);

    auto ta_xxzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 166);

    auto ta_xxzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 167);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xx_xxxxxx_0, ta_xx_xxxxxx_1, ta_xx_xxxxxy_0, ta_xx_xxxxxy_1, ta_xx_xxxxyy_0, ta_xx_xxxxyy_1, ta_xx_xxxyyy_0, ta_xx_xxxyyy_1, ta_xx_xxyyyy_0, ta_xx_xxyyyy_1, ta_xx_xyyyyy_0, ta_xx_xyyyyy_1, ta_xxz_xxxxxx_0, ta_xxz_xxxxxx_1, ta_xxz_xxxxxy_0, ta_xxz_xxxxxy_1, ta_xxz_xxxxyy_0, ta_xxz_xxxxyy_1, ta_xxz_xxxyyy_0, ta_xxz_xxxyyy_1, ta_xxz_xxyyyy_0, ta_xxz_xxyyyy_1, ta_xxz_xyyyyy_0, ta_xxz_xyyyyy_1, ta_xxzz_xxxxxx_0, ta_xxzz_xxxxxy_0, ta_xxzz_xxxxxz_0, ta_xxzz_xxxxyy_0, ta_xxzz_xxxxyz_0, ta_xxzz_xxxxzz_0, ta_xxzz_xxxyyy_0, ta_xxzz_xxxyyz_0, ta_xxzz_xxxyzz_0, ta_xxzz_xxxzzz_0, ta_xxzz_xxyyyy_0, ta_xxzz_xxyyyz_0, ta_xxzz_xxyyzz_0, ta_xxzz_xxyzzz_0, ta_xxzz_xxzzzz_0, ta_xxzz_xyyyyy_0, ta_xxzz_xyyyyz_0, ta_xxzz_xyyyzz_0, ta_xxzz_xyyzzz_0, ta_xxzz_xyzzzz_0, ta_xxzz_xzzzzz_0, ta_xxzz_yyyyyy_0, ta_xxzz_yyyyyz_0, ta_xxzz_yyyyzz_0, ta_xxzz_yyyzzz_0, ta_xxzz_yyzzzz_0, ta_xxzz_yzzzzz_0, ta_xxzz_zzzzzz_0, ta_xzz_xxxxxz_0, ta_xzz_xxxxxz_1, ta_xzz_xxxxyz_0, ta_xzz_xxxxyz_1, ta_xzz_xxxxz_0, ta_xzz_xxxxz_1, ta_xzz_xxxxzz_0, ta_xzz_xxxxzz_1, ta_xzz_xxxyyz_0, ta_xzz_xxxyyz_1, ta_xzz_xxxyz_0, ta_xzz_xxxyz_1, ta_xzz_xxxyzz_0, ta_xzz_xxxyzz_1, ta_xzz_xxxzz_0, ta_xzz_xxxzz_1, ta_xzz_xxxzzz_0, ta_xzz_xxxzzz_1, ta_xzz_xxyyyz_0, ta_xzz_xxyyyz_1, ta_xzz_xxyyz_0, ta_xzz_xxyyz_1, ta_xzz_xxyyzz_0, ta_xzz_xxyyzz_1, ta_xzz_xxyzz_0, ta_xzz_xxyzz_1, ta_xzz_xxyzzz_0, ta_xzz_xxyzzz_1, ta_xzz_xxzzz_0, ta_xzz_xxzzz_1, ta_xzz_xxzzzz_0, ta_xzz_xxzzzz_1, ta_xzz_xyyyyz_0, ta_xzz_xyyyyz_1, ta_xzz_xyyyz_0, ta_xzz_xyyyz_1, ta_xzz_xyyyzz_0, ta_xzz_xyyyzz_1, ta_xzz_xyyzz_0, ta_xzz_xyyzz_1, ta_xzz_xyyzzz_0, ta_xzz_xyyzzz_1, ta_xzz_xyzzz_0, ta_xzz_xyzzz_1, ta_xzz_xyzzzz_0, ta_xzz_xyzzzz_1, ta_xzz_xzzzz_0, ta_xzz_xzzzz_1, ta_xzz_xzzzzz_0, ta_xzz_xzzzzz_1, ta_xzz_yyyyyy_0, ta_xzz_yyyyyy_1, ta_xzz_yyyyyz_0, ta_xzz_yyyyyz_1, ta_xzz_yyyyz_0, ta_xzz_yyyyz_1, ta_xzz_yyyyzz_0, ta_xzz_yyyyzz_1, ta_xzz_yyyzz_0, ta_xzz_yyyzz_1, ta_xzz_yyyzzz_0, ta_xzz_yyyzzz_1, ta_xzz_yyzzz_0, ta_xzz_yyzzz_1, ta_xzz_yyzzzz_0, ta_xzz_yyzzzz_1, ta_xzz_yzzzz_0, ta_xzz_yzzzz_1, ta_xzz_yzzzzz_0, ta_xzz_yzzzzz_1, ta_xzz_zzzzz_0, ta_xzz_zzzzz_1, ta_xzz_zzzzzz_0, ta_xzz_zzzzzz_1, ta_zz_xxxxxz_0, ta_zz_xxxxxz_1, ta_zz_xxxxyz_0, ta_zz_xxxxyz_1, ta_zz_xxxxzz_0, ta_zz_xxxxzz_1, ta_zz_xxxyyz_0, ta_zz_xxxyyz_1, ta_zz_xxxyzz_0, ta_zz_xxxyzz_1, ta_zz_xxxzzz_0, ta_zz_xxxzzz_1, ta_zz_xxyyyz_0, ta_zz_xxyyyz_1, ta_zz_xxyyzz_0, ta_zz_xxyyzz_1, ta_zz_xxyzzz_0, ta_zz_xxyzzz_1, ta_zz_xxzzzz_0, ta_zz_xxzzzz_1, ta_zz_xyyyyz_0, ta_zz_xyyyyz_1, ta_zz_xyyyzz_0, ta_zz_xyyyzz_1, ta_zz_xyyzzz_0, ta_zz_xyyzzz_1, ta_zz_xyzzzz_0, ta_zz_xyzzzz_1, ta_zz_xzzzzz_0, ta_zz_xzzzzz_1, ta_zz_yyyyyy_0, ta_zz_yyyyyy_1, ta_zz_yyyyyz_0, ta_zz_yyyyyz_1, ta_zz_yyyyzz_0, ta_zz_yyyyzz_1, ta_zz_yyyzzz_0, ta_zz_yyyzzz_1, ta_zz_yyzzzz_0, ta_zz_yyzzzz_1, ta_zz_yzzzzz_0, ta_zz_yzzzzz_1, ta_zz_zzzzzz_0, ta_zz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzz_xxxxxx_0[i] = ta_xx_xxxxxx_0[i] * fe_0 - ta_xx_xxxxxx_1[i] * fe_0 + ta_xxz_xxxxxx_0[i] * pa_z[i] - ta_xxz_xxxxxx_1[i] * pc_z[i];

        ta_xxzz_xxxxxy_0[i] = ta_xx_xxxxxy_0[i] * fe_0 - ta_xx_xxxxxy_1[i] * fe_0 + ta_xxz_xxxxxy_0[i] * pa_z[i] - ta_xxz_xxxxxy_1[i] * pc_z[i];

        ta_xxzz_xxxxxz_0[i] = ta_zz_xxxxxz_0[i] * fe_0 - ta_zz_xxxxxz_1[i] * fe_0 + 5.0 * ta_xzz_xxxxz_0[i] * fe_0 - 5.0 * ta_xzz_xxxxz_1[i] * fe_0 + ta_xzz_xxxxxz_0[i] * pa_x[i] - ta_xzz_xxxxxz_1[i] * pc_x[i];

        ta_xxzz_xxxxyy_0[i] = ta_xx_xxxxyy_0[i] * fe_0 - ta_xx_xxxxyy_1[i] * fe_0 + ta_xxz_xxxxyy_0[i] * pa_z[i] - ta_xxz_xxxxyy_1[i] * pc_z[i];

        ta_xxzz_xxxxyz_0[i] = ta_zz_xxxxyz_0[i] * fe_0 - ta_zz_xxxxyz_1[i] * fe_0 + 4.0 * ta_xzz_xxxyz_0[i] * fe_0 - 4.0 * ta_xzz_xxxyz_1[i] * fe_0 + ta_xzz_xxxxyz_0[i] * pa_x[i] - ta_xzz_xxxxyz_1[i] * pc_x[i];

        ta_xxzz_xxxxzz_0[i] = ta_zz_xxxxzz_0[i] * fe_0 - ta_zz_xxxxzz_1[i] * fe_0 + 4.0 * ta_xzz_xxxzz_0[i] * fe_0 - 4.0 * ta_xzz_xxxzz_1[i] * fe_0 + ta_xzz_xxxxzz_0[i] * pa_x[i] - ta_xzz_xxxxzz_1[i] * pc_x[i];

        ta_xxzz_xxxyyy_0[i] = ta_xx_xxxyyy_0[i] * fe_0 - ta_xx_xxxyyy_1[i] * fe_0 + ta_xxz_xxxyyy_0[i] * pa_z[i] - ta_xxz_xxxyyy_1[i] * pc_z[i];

        ta_xxzz_xxxyyz_0[i] = ta_zz_xxxyyz_0[i] * fe_0 - ta_zz_xxxyyz_1[i] * fe_0 + 3.0 * ta_xzz_xxyyz_0[i] * fe_0 - 3.0 * ta_xzz_xxyyz_1[i] * fe_0 + ta_xzz_xxxyyz_0[i] * pa_x[i] - ta_xzz_xxxyyz_1[i] * pc_x[i];

        ta_xxzz_xxxyzz_0[i] = ta_zz_xxxyzz_0[i] * fe_0 - ta_zz_xxxyzz_1[i] * fe_0 + 3.0 * ta_xzz_xxyzz_0[i] * fe_0 - 3.0 * ta_xzz_xxyzz_1[i] * fe_0 + ta_xzz_xxxyzz_0[i] * pa_x[i] - ta_xzz_xxxyzz_1[i] * pc_x[i];

        ta_xxzz_xxxzzz_0[i] = ta_zz_xxxzzz_0[i] * fe_0 - ta_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta_xzz_xxzzz_0[i] * fe_0 - 3.0 * ta_xzz_xxzzz_1[i] * fe_0 + ta_xzz_xxxzzz_0[i] * pa_x[i] - ta_xzz_xxxzzz_1[i] * pc_x[i];

        ta_xxzz_xxyyyy_0[i] = ta_xx_xxyyyy_0[i] * fe_0 - ta_xx_xxyyyy_1[i] * fe_0 + ta_xxz_xxyyyy_0[i] * pa_z[i] - ta_xxz_xxyyyy_1[i] * pc_z[i];

        ta_xxzz_xxyyyz_0[i] = ta_zz_xxyyyz_0[i] * fe_0 - ta_zz_xxyyyz_1[i] * fe_0 + 2.0 * ta_xzz_xyyyz_0[i] * fe_0 - 2.0 * ta_xzz_xyyyz_1[i] * fe_0 + ta_xzz_xxyyyz_0[i] * pa_x[i] - ta_xzz_xxyyyz_1[i] * pc_x[i];

        ta_xxzz_xxyyzz_0[i] = ta_zz_xxyyzz_0[i] * fe_0 - ta_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta_xzz_xyyzz_0[i] * fe_0 - 2.0 * ta_xzz_xyyzz_1[i] * fe_0 + ta_xzz_xxyyzz_0[i] * pa_x[i] - ta_xzz_xxyyzz_1[i] * pc_x[i];

        ta_xxzz_xxyzzz_0[i] = ta_zz_xxyzzz_0[i] * fe_0 - ta_zz_xxyzzz_1[i] * fe_0 + 2.0 * ta_xzz_xyzzz_0[i] * fe_0 - 2.0 * ta_xzz_xyzzz_1[i] * fe_0 + ta_xzz_xxyzzz_0[i] * pa_x[i] - ta_xzz_xxyzzz_1[i] * pc_x[i];

        ta_xxzz_xxzzzz_0[i] = ta_zz_xxzzzz_0[i] * fe_0 - ta_zz_xxzzzz_1[i] * fe_0 + 2.0 * ta_xzz_xzzzz_0[i] * fe_0 - 2.0 * ta_xzz_xzzzz_1[i] * fe_0 + ta_xzz_xxzzzz_0[i] * pa_x[i] - ta_xzz_xxzzzz_1[i] * pc_x[i];

        ta_xxzz_xyyyyy_0[i] = ta_xx_xyyyyy_0[i] * fe_0 - ta_xx_xyyyyy_1[i] * fe_0 + ta_xxz_xyyyyy_0[i] * pa_z[i] - ta_xxz_xyyyyy_1[i] * pc_z[i];

        ta_xxzz_xyyyyz_0[i] = ta_zz_xyyyyz_0[i] * fe_0 - ta_zz_xyyyyz_1[i] * fe_0 + ta_xzz_yyyyz_0[i] * fe_0 - ta_xzz_yyyyz_1[i] * fe_0 + ta_xzz_xyyyyz_0[i] * pa_x[i] - ta_xzz_xyyyyz_1[i] * pc_x[i];

        ta_xxzz_xyyyzz_0[i] = ta_zz_xyyyzz_0[i] * fe_0 - ta_zz_xyyyzz_1[i] * fe_0 + ta_xzz_yyyzz_0[i] * fe_0 - ta_xzz_yyyzz_1[i] * fe_0 + ta_xzz_xyyyzz_0[i] * pa_x[i] - ta_xzz_xyyyzz_1[i] * pc_x[i];

        ta_xxzz_xyyzzz_0[i] = ta_zz_xyyzzz_0[i] * fe_0 - ta_zz_xyyzzz_1[i] * fe_0 + ta_xzz_yyzzz_0[i] * fe_0 - ta_xzz_yyzzz_1[i] * fe_0 + ta_xzz_xyyzzz_0[i] * pa_x[i] - ta_xzz_xyyzzz_1[i] * pc_x[i];

        ta_xxzz_xyzzzz_0[i] = ta_zz_xyzzzz_0[i] * fe_0 - ta_zz_xyzzzz_1[i] * fe_0 + ta_xzz_yzzzz_0[i] * fe_0 - ta_xzz_yzzzz_1[i] * fe_0 + ta_xzz_xyzzzz_0[i] * pa_x[i] - ta_xzz_xyzzzz_1[i] * pc_x[i];

        ta_xxzz_xzzzzz_0[i] = ta_zz_xzzzzz_0[i] * fe_0 - ta_zz_xzzzzz_1[i] * fe_0 + ta_xzz_zzzzz_0[i] * fe_0 - ta_xzz_zzzzz_1[i] * fe_0 + ta_xzz_xzzzzz_0[i] * pa_x[i] - ta_xzz_xzzzzz_1[i] * pc_x[i];

        ta_xxzz_yyyyyy_0[i] = ta_zz_yyyyyy_0[i] * fe_0 - ta_zz_yyyyyy_1[i] * fe_0 + ta_xzz_yyyyyy_0[i] * pa_x[i] - ta_xzz_yyyyyy_1[i] * pc_x[i];

        ta_xxzz_yyyyyz_0[i] = ta_zz_yyyyyz_0[i] * fe_0 - ta_zz_yyyyyz_1[i] * fe_0 + ta_xzz_yyyyyz_0[i] * pa_x[i] - ta_xzz_yyyyyz_1[i] * pc_x[i];

        ta_xxzz_yyyyzz_0[i] = ta_zz_yyyyzz_0[i] * fe_0 - ta_zz_yyyyzz_1[i] * fe_0 + ta_xzz_yyyyzz_0[i] * pa_x[i] - ta_xzz_yyyyzz_1[i] * pc_x[i];

        ta_xxzz_yyyzzz_0[i] = ta_zz_yyyzzz_0[i] * fe_0 - ta_zz_yyyzzz_1[i] * fe_0 + ta_xzz_yyyzzz_0[i] * pa_x[i] - ta_xzz_yyyzzz_1[i] * pc_x[i];

        ta_xxzz_yyzzzz_0[i] = ta_zz_yyzzzz_0[i] * fe_0 - ta_zz_yyzzzz_1[i] * fe_0 + ta_xzz_yyzzzz_0[i] * pa_x[i] - ta_xzz_yyzzzz_1[i] * pc_x[i];

        ta_xxzz_yzzzzz_0[i] = ta_zz_yzzzzz_0[i] * fe_0 - ta_zz_yzzzzz_1[i] * fe_0 + ta_xzz_yzzzzz_0[i] * pa_x[i] - ta_xzz_yzzzzz_1[i] * pc_x[i];

        ta_xxzz_zzzzzz_0[i] = ta_zz_zzzzzz_0[i] * fe_0 - ta_zz_zzzzzz_1[i] * fe_0 + ta_xzz_zzzzzz_0[i] * pa_x[i] - ta_xzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 168-196 components of targeted buffer : GI

    auto ta_xyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 168);

    auto ta_xyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 169);

    auto ta_xyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 170);

    auto ta_xyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 171);

    auto ta_xyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 172);

    auto ta_xyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 173);

    auto ta_xyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 174);

    auto ta_xyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 175);

    auto ta_xyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 176);

    auto ta_xyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 177);

    auto ta_xyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 178);

    auto ta_xyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 179);

    auto ta_xyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 180);

    auto ta_xyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 181);

    auto ta_xyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 182);

    auto ta_xyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 183);

    auto ta_xyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 184);

    auto ta_xyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 185);

    auto ta_xyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 186);

    auto ta_xyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 187);

    auto ta_xyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 188);

    auto ta_xyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 189);

    auto ta_xyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 190);

    auto ta_xyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 191);

    auto ta_xyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 192);

    auto ta_xyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 193);

    auto ta_xyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 194);

    auto ta_xyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 195);

    #pragma omp simd aligned(pa_x, pc_x, ta_xyyy_xxxxxx_0, ta_xyyy_xxxxxy_0, ta_xyyy_xxxxxz_0, ta_xyyy_xxxxyy_0, ta_xyyy_xxxxyz_0, ta_xyyy_xxxxzz_0, ta_xyyy_xxxyyy_0, ta_xyyy_xxxyyz_0, ta_xyyy_xxxyzz_0, ta_xyyy_xxxzzz_0, ta_xyyy_xxyyyy_0, ta_xyyy_xxyyyz_0, ta_xyyy_xxyyzz_0, ta_xyyy_xxyzzz_0, ta_xyyy_xxzzzz_0, ta_xyyy_xyyyyy_0, ta_xyyy_xyyyyz_0, ta_xyyy_xyyyzz_0, ta_xyyy_xyyzzz_0, ta_xyyy_xyzzzz_0, ta_xyyy_xzzzzz_0, ta_xyyy_yyyyyy_0, ta_xyyy_yyyyyz_0, ta_xyyy_yyyyzz_0, ta_xyyy_yyyzzz_0, ta_xyyy_yyzzzz_0, ta_xyyy_yzzzzz_0, ta_xyyy_zzzzzz_0, ta_yyy_xxxxx_0, ta_yyy_xxxxx_1, ta_yyy_xxxxxx_0, ta_yyy_xxxxxx_1, ta_yyy_xxxxxy_0, ta_yyy_xxxxxy_1, ta_yyy_xxxxxz_0, ta_yyy_xxxxxz_1, ta_yyy_xxxxy_0, ta_yyy_xxxxy_1, ta_yyy_xxxxyy_0, ta_yyy_xxxxyy_1, ta_yyy_xxxxyz_0, ta_yyy_xxxxyz_1, ta_yyy_xxxxz_0, ta_yyy_xxxxz_1, ta_yyy_xxxxzz_0, ta_yyy_xxxxzz_1, ta_yyy_xxxyy_0, ta_yyy_xxxyy_1, ta_yyy_xxxyyy_0, ta_yyy_xxxyyy_1, ta_yyy_xxxyyz_0, ta_yyy_xxxyyz_1, ta_yyy_xxxyz_0, ta_yyy_xxxyz_1, ta_yyy_xxxyzz_0, ta_yyy_xxxyzz_1, ta_yyy_xxxzz_0, ta_yyy_xxxzz_1, ta_yyy_xxxzzz_0, ta_yyy_xxxzzz_1, ta_yyy_xxyyy_0, ta_yyy_xxyyy_1, ta_yyy_xxyyyy_0, ta_yyy_xxyyyy_1, ta_yyy_xxyyyz_0, ta_yyy_xxyyyz_1, ta_yyy_xxyyz_0, ta_yyy_xxyyz_1, ta_yyy_xxyyzz_0, ta_yyy_xxyyzz_1, ta_yyy_xxyzz_0, ta_yyy_xxyzz_1, ta_yyy_xxyzzz_0, ta_yyy_xxyzzz_1, ta_yyy_xxzzz_0, ta_yyy_xxzzz_1, ta_yyy_xxzzzz_0, ta_yyy_xxzzzz_1, ta_yyy_xyyyy_0, ta_yyy_xyyyy_1, ta_yyy_xyyyyy_0, ta_yyy_xyyyyy_1, ta_yyy_xyyyyz_0, ta_yyy_xyyyyz_1, ta_yyy_xyyyz_0, ta_yyy_xyyyz_1, ta_yyy_xyyyzz_0, ta_yyy_xyyyzz_1, ta_yyy_xyyzz_0, ta_yyy_xyyzz_1, ta_yyy_xyyzzz_0, ta_yyy_xyyzzz_1, ta_yyy_xyzzz_0, ta_yyy_xyzzz_1, ta_yyy_xyzzzz_0, ta_yyy_xyzzzz_1, ta_yyy_xzzzz_0, ta_yyy_xzzzz_1, ta_yyy_xzzzzz_0, ta_yyy_xzzzzz_1, ta_yyy_yyyyy_0, ta_yyy_yyyyy_1, ta_yyy_yyyyyy_0, ta_yyy_yyyyyy_1, ta_yyy_yyyyyz_0, ta_yyy_yyyyyz_1, ta_yyy_yyyyz_0, ta_yyy_yyyyz_1, ta_yyy_yyyyzz_0, ta_yyy_yyyyzz_1, ta_yyy_yyyzz_0, ta_yyy_yyyzz_1, ta_yyy_yyyzzz_0, ta_yyy_yyyzzz_1, ta_yyy_yyzzz_0, ta_yyy_yyzzz_1, ta_yyy_yyzzzz_0, ta_yyy_yyzzzz_1, ta_yyy_yzzzz_0, ta_yyy_yzzzz_1, ta_yyy_yzzzzz_0, ta_yyy_yzzzzz_1, ta_yyy_zzzzz_0, ta_yyy_zzzzz_1, ta_yyy_zzzzzz_0, ta_yyy_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyy_xxxxxx_0[i] = 6.0 * ta_yyy_xxxxx_0[i] * fe_0 - 6.0 * ta_yyy_xxxxx_1[i] * fe_0 + ta_yyy_xxxxxx_0[i] * pa_x[i] - ta_yyy_xxxxxx_1[i] * pc_x[i];

        ta_xyyy_xxxxxy_0[i] = 5.0 * ta_yyy_xxxxy_0[i] * fe_0 - 5.0 * ta_yyy_xxxxy_1[i] * fe_0 + ta_yyy_xxxxxy_0[i] * pa_x[i] - ta_yyy_xxxxxy_1[i] * pc_x[i];

        ta_xyyy_xxxxxz_0[i] = 5.0 * ta_yyy_xxxxz_0[i] * fe_0 - 5.0 * ta_yyy_xxxxz_1[i] * fe_0 + ta_yyy_xxxxxz_0[i] * pa_x[i] - ta_yyy_xxxxxz_1[i] * pc_x[i];

        ta_xyyy_xxxxyy_0[i] = 4.0 * ta_yyy_xxxyy_0[i] * fe_0 - 4.0 * ta_yyy_xxxyy_1[i] * fe_0 + ta_yyy_xxxxyy_0[i] * pa_x[i] - ta_yyy_xxxxyy_1[i] * pc_x[i];

        ta_xyyy_xxxxyz_0[i] = 4.0 * ta_yyy_xxxyz_0[i] * fe_0 - 4.0 * ta_yyy_xxxyz_1[i] * fe_0 + ta_yyy_xxxxyz_0[i] * pa_x[i] - ta_yyy_xxxxyz_1[i] * pc_x[i];

        ta_xyyy_xxxxzz_0[i] = 4.0 * ta_yyy_xxxzz_0[i] * fe_0 - 4.0 * ta_yyy_xxxzz_1[i] * fe_0 + ta_yyy_xxxxzz_0[i] * pa_x[i] - ta_yyy_xxxxzz_1[i] * pc_x[i];

        ta_xyyy_xxxyyy_0[i] = 3.0 * ta_yyy_xxyyy_0[i] * fe_0 - 3.0 * ta_yyy_xxyyy_1[i] * fe_0 + ta_yyy_xxxyyy_0[i] * pa_x[i] - ta_yyy_xxxyyy_1[i] * pc_x[i];

        ta_xyyy_xxxyyz_0[i] = 3.0 * ta_yyy_xxyyz_0[i] * fe_0 - 3.0 * ta_yyy_xxyyz_1[i] * fe_0 + ta_yyy_xxxyyz_0[i] * pa_x[i] - ta_yyy_xxxyyz_1[i] * pc_x[i];

        ta_xyyy_xxxyzz_0[i] = 3.0 * ta_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta_yyy_xxyzz_1[i] * fe_0 + ta_yyy_xxxyzz_0[i] * pa_x[i] - ta_yyy_xxxyzz_1[i] * pc_x[i];

        ta_xyyy_xxxzzz_0[i] = 3.0 * ta_yyy_xxzzz_0[i] * fe_0 - 3.0 * ta_yyy_xxzzz_1[i] * fe_0 + ta_yyy_xxxzzz_0[i] * pa_x[i] - ta_yyy_xxxzzz_1[i] * pc_x[i];

        ta_xyyy_xxyyyy_0[i] = 2.0 * ta_yyy_xyyyy_0[i] * fe_0 - 2.0 * ta_yyy_xyyyy_1[i] * fe_0 + ta_yyy_xxyyyy_0[i] * pa_x[i] - ta_yyy_xxyyyy_1[i] * pc_x[i];

        ta_xyyy_xxyyyz_0[i] = 2.0 * ta_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta_yyy_xyyyz_1[i] * fe_0 + ta_yyy_xxyyyz_0[i] * pa_x[i] - ta_yyy_xxyyyz_1[i] * pc_x[i];

        ta_xyyy_xxyyzz_0[i] = 2.0 * ta_yyy_xyyzz_0[i] * fe_0 - 2.0 * ta_yyy_xyyzz_1[i] * fe_0 + ta_yyy_xxyyzz_0[i] * pa_x[i] - ta_yyy_xxyyzz_1[i] * pc_x[i];

        ta_xyyy_xxyzzz_0[i] = 2.0 * ta_yyy_xyzzz_0[i] * fe_0 - 2.0 * ta_yyy_xyzzz_1[i] * fe_0 + ta_yyy_xxyzzz_0[i] * pa_x[i] - ta_yyy_xxyzzz_1[i] * pc_x[i];

        ta_xyyy_xxzzzz_0[i] = 2.0 * ta_yyy_xzzzz_0[i] * fe_0 - 2.0 * ta_yyy_xzzzz_1[i] * fe_0 + ta_yyy_xxzzzz_0[i] * pa_x[i] - ta_yyy_xxzzzz_1[i] * pc_x[i];

        ta_xyyy_xyyyyy_0[i] = ta_yyy_yyyyy_0[i] * fe_0 - ta_yyy_yyyyy_1[i] * fe_0 + ta_yyy_xyyyyy_0[i] * pa_x[i] - ta_yyy_xyyyyy_1[i] * pc_x[i];

        ta_xyyy_xyyyyz_0[i] = ta_yyy_yyyyz_0[i] * fe_0 - ta_yyy_yyyyz_1[i] * fe_0 + ta_yyy_xyyyyz_0[i] * pa_x[i] - ta_yyy_xyyyyz_1[i] * pc_x[i];

        ta_xyyy_xyyyzz_0[i] = ta_yyy_yyyzz_0[i] * fe_0 - ta_yyy_yyyzz_1[i] * fe_0 + ta_yyy_xyyyzz_0[i] * pa_x[i] - ta_yyy_xyyyzz_1[i] * pc_x[i];

        ta_xyyy_xyyzzz_0[i] = ta_yyy_yyzzz_0[i] * fe_0 - ta_yyy_yyzzz_1[i] * fe_0 + ta_yyy_xyyzzz_0[i] * pa_x[i] - ta_yyy_xyyzzz_1[i] * pc_x[i];

        ta_xyyy_xyzzzz_0[i] = ta_yyy_yzzzz_0[i] * fe_0 - ta_yyy_yzzzz_1[i] * fe_0 + ta_yyy_xyzzzz_0[i] * pa_x[i] - ta_yyy_xyzzzz_1[i] * pc_x[i];

        ta_xyyy_xzzzzz_0[i] = ta_yyy_zzzzz_0[i] * fe_0 - ta_yyy_zzzzz_1[i] * fe_0 + ta_yyy_xzzzzz_0[i] * pa_x[i] - ta_yyy_xzzzzz_1[i] * pc_x[i];

        ta_xyyy_yyyyyy_0[i] = ta_yyy_yyyyyy_0[i] * pa_x[i] - ta_yyy_yyyyyy_1[i] * pc_x[i];

        ta_xyyy_yyyyyz_0[i] = ta_yyy_yyyyyz_0[i] * pa_x[i] - ta_yyy_yyyyyz_1[i] * pc_x[i];

        ta_xyyy_yyyyzz_0[i] = ta_yyy_yyyyzz_0[i] * pa_x[i] - ta_yyy_yyyyzz_1[i] * pc_x[i];

        ta_xyyy_yyyzzz_0[i] = ta_yyy_yyyzzz_0[i] * pa_x[i] - ta_yyy_yyyzzz_1[i] * pc_x[i];

        ta_xyyy_yyzzzz_0[i] = ta_yyy_yyzzzz_0[i] * pa_x[i] - ta_yyy_yyzzzz_1[i] * pc_x[i];

        ta_xyyy_yzzzzz_0[i] = ta_yyy_yzzzzz_0[i] * pa_x[i] - ta_yyy_yzzzzz_1[i] * pc_x[i];

        ta_xyyy_zzzzzz_0[i] = ta_yyy_zzzzzz_0[i] * pa_x[i] - ta_yyy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 196-224 components of targeted buffer : GI

    auto ta_xyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 196);

    auto ta_xyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 197);

    auto ta_xyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 198);

    auto ta_xyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 199);

    auto ta_xyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 200);

    auto ta_xyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 201);

    auto ta_xyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 202);

    auto ta_xyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 203);

    auto ta_xyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 204);

    auto ta_xyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 205);

    auto ta_xyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 206);

    auto ta_xyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 207);

    auto ta_xyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 208);

    auto ta_xyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 209);

    auto ta_xyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 210);

    auto ta_xyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 211);

    auto ta_xyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 212);

    auto ta_xyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 213);

    auto ta_xyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 214);

    auto ta_xyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 215);

    auto ta_xyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 216);

    auto ta_xyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 217);

    auto ta_xyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 218);

    auto ta_xyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 219);

    auto ta_xyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 220);

    auto ta_xyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 221);

    auto ta_xyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 222);

    auto ta_xyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 223);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta_xyy_xxxxxx_0, ta_xyy_xxxxxx_1, ta_xyy_xxxxxy_0, ta_xyy_xxxxxy_1, ta_xyy_xxxxyy_0, ta_xyy_xxxxyy_1, ta_xyy_xxxyyy_0, ta_xyy_xxxyyy_1, ta_xyy_xxyyyy_0, ta_xyy_xxyyyy_1, ta_xyy_xyyyyy_0, ta_xyy_xyyyyy_1, ta_xyyz_xxxxxx_0, ta_xyyz_xxxxxy_0, ta_xyyz_xxxxxz_0, ta_xyyz_xxxxyy_0, ta_xyyz_xxxxyz_0, ta_xyyz_xxxxzz_0, ta_xyyz_xxxyyy_0, ta_xyyz_xxxyyz_0, ta_xyyz_xxxyzz_0, ta_xyyz_xxxzzz_0, ta_xyyz_xxyyyy_0, ta_xyyz_xxyyyz_0, ta_xyyz_xxyyzz_0, ta_xyyz_xxyzzz_0, ta_xyyz_xxzzzz_0, ta_xyyz_xyyyyy_0, ta_xyyz_xyyyyz_0, ta_xyyz_xyyyzz_0, ta_xyyz_xyyzzz_0, ta_xyyz_xyzzzz_0, ta_xyyz_xzzzzz_0, ta_xyyz_yyyyyy_0, ta_xyyz_yyyyyz_0, ta_xyyz_yyyyzz_0, ta_xyyz_yyyzzz_0, ta_xyyz_yyzzzz_0, ta_xyyz_yzzzzz_0, ta_xyyz_zzzzzz_0, ta_yyz_xxxxxz_0, ta_yyz_xxxxxz_1, ta_yyz_xxxxyz_0, ta_yyz_xxxxyz_1, ta_yyz_xxxxz_0, ta_yyz_xxxxz_1, ta_yyz_xxxxzz_0, ta_yyz_xxxxzz_1, ta_yyz_xxxyyz_0, ta_yyz_xxxyyz_1, ta_yyz_xxxyz_0, ta_yyz_xxxyz_1, ta_yyz_xxxyzz_0, ta_yyz_xxxyzz_1, ta_yyz_xxxzz_0, ta_yyz_xxxzz_1, ta_yyz_xxxzzz_0, ta_yyz_xxxzzz_1, ta_yyz_xxyyyz_0, ta_yyz_xxyyyz_1, ta_yyz_xxyyz_0, ta_yyz_xxyyz_1, ta_yyz_xxyyzz_0, ta_yyz_xxyyzz_1, ta_yyz_xxyzz_0, ta_yyz_xxyzz_1, ta_yyz_xxyzzz_0, ta_yyz_xxyzzz_1, ta_yyz_xxzzz_0, ta_yyz_xxzzz_1, ta_yyz_xxzzzz_0, ta_yyz_xxzzzz_1, ta_yyz_xyyyyz_0, ta_yyz_xyyyyz_1, ta_yyz_xyyyz_0, ta_yyz_xyyyz_1, ta_yyz_xyyyzz_0, ta_yyz_xyyyzz_1, ta_yyz_xyyzz_0, ta_yyz_xyyzz_1, ta_yyz_xyyzzz_0, ta_yyz_xyyzzz_1, ta_yyz_xyzzz_0, ta_yyz_xyzzz_1, ta_yyz_xyzzzz_0, ta_yyz_xyzzzz_1, ta_yyz_xzzzz_0, ta_yyz_xzzzz_1, ta_yyz_xzzzzz_0, ta_yyz_xzzzzz_1, ta_yyz_yyyyyy_0, ta_yyz_yyyyyy_1, ta_yyz_yyyyyz_0, ta_yyz_yyyyyz_1, ta_yyz_yyyyz_0, ta_yyz_yyyyz_1, ta_yyz_yyyyzz_0, ta_yyz_yyyyzz_1, ta_yyz_yyyzz_0, ta_yyz_yyyzz_1, ta_yyz_yyyzzz_0, ta_yyz_yyyzzz_1, ta_yyz_yyzzz_0, ta_yyz_yyzzz_1, ta_yyz_yyzzzz_0, ta_yyz_yyzzzz_1, ta_yyz_yzzzz_0, ta_yyz_yzzzz_1, ta_yyz_yzzzzz_0, ta_yyz_yzzzzz_1, ta_yyz_zzzzz_0, ta_yyz_zzzzz_1, ta_yyz_zzzzzz_0, ta_yyz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyz_xxxxxx_0[i] = ta_xyy_xxxxxx_0[i] * pa_z[i] - ta_xyy_xxxxxx_1[i] * pc_z[i];

        ta_xyyz_xxxxxy_0[i] = ta_xyy_xxxxxy_0[i] * pa_z[i] - ta_xyy_xxxxxy_1[i] * pc_z[i];

        ta_xyyz_xxxxxz_0[i] = 5.0 * ta_yyz_xxxxz_0[i] * fe_0 - 5.0 * ta_yyz_xxxxz_1[i] * fe_0 + ta_yyz_xxxxxz_0[i] * pa_x[i] - ta_yyz_xxxxxz_1[i] * pc_x[i];

        ta_xyyz_xxxxyy_0[i] = ta_xyy_xxxxyy_0[i] * pa_z[i] - ta_xyy_xxxxyy_1[i] * pc_z[i];

        ta_xyyz_xxxxyz_0[i] = 4.0 * ta_yyz_xxxyz_0[i] * fe_0 - 4.0 * ta_yyz_xxxyz_1[i] * fe_0 + ta_yyz_xxxxyz_0[i] * pa_x[i] - ta_yyz_xxxxyz_1[i] * pc_x[i];

        ta_xyyz_xxxxzz_0[i] = 4.0 * ta_yyz_xxxzz_0[i] * fe_0 - 4.0 * ta_yyz_xxxzz_1[i] * fe_0 + ta_yyz_xxxxzz_0[i] * pa_x[i] - ta_yyz_xxxxzz_1[i] * pc_x[i];

        ta_xyyz_xxxyyy_0[i] = ta_xyy_xxxyyy_0[i] * pa_z[i] - ta_xyy_xxxyyy_1[i] * pc_z[i];

        ta_xyyz_xxxyyz_0[i] = 3.0 * ta_yyz_xxyyz_0[i] * fe_0 - 3.0 * ta_yyz_xxyyz_1[i] * fe_0 + ta_yyz_xxxyyz_0[i] * pa_x[i] - ta_yyz_xxxyyz_1[i] * pc_x[i];

        ta_xyyz_xxxyzz_0[i] = 3.0 * ta_yyz_xxyzz_0[i] * fe_0 - 3.0 * ta_yyz_xxyzz_1[i] * fe_0 + ta_yyz_xxxyzz_0[i] * pa_x[i] - ta_yyz_xxxyzz_1[i] * pc_x[i];

        ta_xyyz_xxxzzz_0[i] = 3.0 * ta_yyz_xxzzz_0[i] * fe_0 - 3.0 * ta_yyz_xxzzz_1[i] * fe_0 + ta_yyz_xxxzzz_0[i] * pa_x[i] - ta_yyz_xxxzzz_1[i] * pc_x[i];

        ta_xyyz_xxyyyy_0[i] = ta_xyy_xxyyyy_0[i] * pa_z[i] - ta_xyy_xxyyyy_1[i] * pc_z[i];

        ta_xyyz_xxyyyz_0[i] = 2.0 * ta_yyz_xyyyz_0[i] * fe_0 - 2.0 * ta_yyz_xyyyz_1[i] * fe_0 + ta_yyz_xxyyyz_0[i] * pa_x[i] - ta_yyz_xxyyyz_1[i] * pc_x[i];

        ta_xyyz_xxyyzz_0[i] = 2.0 * ta_yyz_xyyzz_0[i] * fe_0 - 2.0 * ta_yyz_xyyzz_1[i] * fe_0 + ta_yyz_xxyyzz_0[i] * pa_x[i] - ta_yyz_xxyyzz_1[i] * pc_x[i];

        ta_xyyz_xxyzzz_0[i] = 2.0 * ta_yyz_xyzzz_0[i] * fe_0 - 2.0 * ta_yyz_xyzzz_1[i] * fe_0 + ta_yyz_xxyzzz_0[i] * pa_x[i] - ta_yyz_xxyzzz_1[i] * pc_x[i];

        ta_xyyz_xxzzzz_0[i] = 2.0 * ta_yyz_xzzzz_0[i] * fe_0 - 2.0 * ta_yyz_xzzzz_1[i] * fe_0 + ta_yyz_xxzzzz_0[i] * pa_x[i] - ta_yyz_xxzzzz_1[i] * pc_x[i];

        ta_xyyz_xyyyyy_0[i] = ta_xyy_xyyyyy_0[i] * pa_z[i] - ta_xyy_xyyyyy_1[i] * pc_z[i];

        ta_xyyz_xyyyyz_0[i] = ta_yyz_yyyyz_0[i] * fe_0 - ta_yyz_yyyyz_1[i] * fe_0 + ta_yyz_xyyyyz_0[i] * pa_x[i] - ta_yyz_xyyyyz_1[i] * pc_x[i];

        ta_xyyz_xyyyzz_0[i] = ta_yyz_yyyzz_0[i] * fe_0 - ta_yyz_yyyzz_1[i] * fe_0 + ta_yyz_xyyyzz_0[i] * pa_x[i] - ta_yyz_xyyyzz_1[i] * pc_x[i];

        ta_xyyz_xyyzzz_0[i] = ta_yyz_yyzzz_0[i] * fe_0 - ta_yyz_yyzzz_1[i] * fe_0 + ta_yyz_xyyzzz_0[i] * pa_x[i] - ta_yyz_xyyzzz_1[i] * pc_x[i];

        ta_xyyz_xyzzzz_0[i] = ta_yyz_yzzzz_0[i] * fe_0 - ta_yyz_yzzzz_1[i] * fe_0 + ta_yyz_xyzzzz_0[i] * pa_x[i] - ta_yyz_xyzzzz_1[i] * pc_x[i];

        ta_xyyz_xzzzzz_0[i] = ta_yyz_zzzzz_0[i] * fe_0 - ta_yyz_zzzzz_1[i] * fe_0 + ta_yyz_xzzzzz_0[i] * pa_x[i] - ta_yyz_xzzzzz_1[i] * pc_x[i];

        ta_xyyz_yyyyyy_0[i] = ta_yyz_yyyyyy_0[i] * pa_x[i] - ta_yyz_yyyyyy_1[i] * pc_x[i];

        ta_xyyz_yyyyyz_0[i] = ta_yyz_yyyyyz_0[i] * pa_x[i] - ta_yyz_yyyyyz_1[i] * pc_x[i];

        ta_xyyz_yyyyzz_0[i] = ta_yyz_yyyyzz_0[i] * pa_x[i] - ta_yyz_yyyyzz_1[i] * pc_x[i];

        ta_xyyz_yyyzzz_0[i] = ta_yyz_yyyzzz_0[i] * pa_x[i] - ta_yyz_yyyzzz_1[i] * pc_x[i];

        ta_xyyz_yyzzzz_0[i] = ta_yyz_yyzzzz_0[i] * pa_x[i] - ta_yyz_yyzzzz_1[i] * pc_x[i];

        ta_xyyz_yzzzzz_0[i] = ta_yyz_yzzzzz_0[i] * pa_x[i] - ta_yyz_yzzzzz_1[i] * pc_x[i];

        ta_xyyz_zzzzzz_0[i] = ta_yyz_zzzzzz_0[i] * pa_x[i] - ta_yyz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 224-252 components of targeted buffer : GI

    auto ta_xyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 224);

    auto ta_xyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 225);

    auto ta_xyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 226);

    auto ta_xyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 227);

    auto ta_xyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 228);

    auto ta_xyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 229);

    auto ta_xyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 230);

    auto ta_xyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 231);

    auto ta_xyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 232);

    auto ta_xyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 233);

    auto ta_xyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 234);

    auto ta_xyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 235);

    auto ta_xyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 236);

    auto ta_xyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 237);

    auto ta_xyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 238);

    auto ta_xyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 239);

    auto ta_xyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 240);

    auto ta_xyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 241);

    auto ta_xyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 242);

    auto ta_xyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 243);

    auto ta_xyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 244);

    auto ta_xyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 245);

    auto ta_xyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 246);

    auto ta_xyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 247);

    auto ta_xyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 248);

    auto ta_xyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 249);

    auto ta_xyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 250);

    auto ta_xyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 251);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta_xyzz_xxxxxx_0, ta_xyzz_xxxxxy_0, ta_xyzz_xxxxxz_0, ta_xyzz_xxxxyy_0, ta_xyzz_xxxxyz_0, ta_xyzz_xxxxzz_0, ta_xyzz_xxxyyy_0, ta_xyzz_xxxyyz_0, ta_xyzz_xxxyzz_0, ta_xyzz_xxxzzz_0, ta_xyzz_xxyyyy_0, ta_xyzz_xxyyyz_0, ta_xyzz_xxyyzz_0, ta_xyzz_xxyzzz_0, ta_xyzz_xxzzzz_0, ta_xyzz_xyyyyy_0, ta_xyzz_xyyyyz_0, ta_xyzz_xyyyzz_0, ta_xyzz_xyyzzz_0, ta_xyzz_xyzzzz_0, ta_xyzz_xzzzzz_0, ta_xyzz_yyyyyy_0, ta_xyzz_yyyyyz_0, ta_xyzz_yyyyzz_0, ta_xyzz_yyyzzz_0, ta_xyzz_yyzzzz_0, ta_xyzz_yzzzzz_0, ta_xyzz_zzzzzz_0, ta_xzz_xxxxxx_0, ta_xzz_xxxxxx_1, ta_xzz_xxxxxz_0, ta_xzz_xxxxxz_1, ta_xzz_xxxxzz_0, ta_xzz_xxxxzz_1, ta_xzz_xxxzzz_0, ta_xzz_xxxzzz_1, ta_xzz_xxzzzz_0, ta_xzz_xxzzzz_1, ta_xzz_xzzzzz_0, ta_xzz_xzzzzz_1, ta_yzz_xxxxxy_0, ta_yzz_xxxxxy_1, ta_yzz_xxxxy_0, ta_yzz_xxxxy_1, ta_yzz_xxxxyy_0, ta_yzz_xxxxyy_1, ta_yzz_xxxxyz_0, ta_yzz_xxxxyz_1, ta_yzz_xxxyy_0, ta_yzz_xxxyy_1, ta_yzz_xxxyyy_0, ta_yzz_xxxyyy_1, ta_yzz_xxxyyz_0, ta_yzz_xxxyyz_1, ta_yzz_xxxyz_0, ta_yzz_xxxyz_1, ta_yzz_xxxyzz_0, ta_yzz_xxxyzz_1, ta_yzz_xxyyy_0, ta_yzz_xxyyy_1, ta_yzz_xxyyyy_0, ta_yzz_xxyyyy_1, ta_yzz_xxyyyz_0, ta_yzz_xxyyyz_1, ta_yzz_xxyyz_0, ta_yzz_xxyyz_1, ta_yzz_xxyyzz_0, ta_yzz_xxyyzz_1, ta_yzz_xxyzz_0, ta_yzz_xxyzz_1, ta_yzz_xxyzzz_0, ta_yzz_xxyzzz_1, ta_yzz_xyyyy_0, ta_yzz_xyyyy_1, ta_yzz_xyyyyy_0, ta_yzz_xyyyyy_1, ta_yzz_xyyyyz_0, ta_yzz_xyyyyz_1, ta_yzz_xyyyz_0, ta_yzz_xyyyz_1, ta_yzz_xyyyzz_0, ta_yzz_xyyyzz_1, ta_yzz_xyyzz_0, ta_yzz_xyyzz_1, ta_yzz_xyyzzz_0, ta_yzz_xyyzzz_1, ta_yzz_xyzzz_0, ta_yzz_xyzzz_1, ta_yzz_xyzzzz_0, ta_yzz_xyzzzz_1, ta_yzz_yyyyy_0, ta_yzz_yyyyy_1, ta_yzz_yyyyyy_0, ta_yzz_yyyyyy_1, ta_yzz_yyyyyz_0, ta_yzz_yyyyyz_1, ta_yzz_yyyyz_0, ta_yzz_yyyyz_1, ta_yzz_yyyyzz_0, ta_yzz_yyyyzz_1, ta_yzz_yyyzz_0, ta_yzz_yyyzz_1, ta_yzz_yyyzzz_0, ta_yzz_yyyzzz_1, ta_yzz_yyzzz_0, ta_yzz_yyzzz_1, ta_yzz_yyzzzz_0, ta_yzz_yyzzzz_1, ta_yzz_yzzzz_0, ta_yzz_yzzzz_1, ta_yzz_yzzzzz_0, ta_yzz_yzzzzz_1, ta_yzz_zzzzzz_0, ta_yzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzz_xxxxxx_0[i] = ta_xzz_xxxxxx_0[i] * pa_y[i] - ta_xzz_xxxxxx_1[i] * pc_y[i];

        ta_xyzz_xxxxxy_0[i] = 5.0 * ta_yzz_xxxxy_0[i] * fe_0 - 5.0 * ta_yzz_xxxxy_1[i] * fe_0 + ta_yzz_xxxxxy_0[i] * pa_x[i] - ta_yzz_xxxxxy_1[i] * pc_x[i];

        ta_xyzz_xxxxxz_0[i] = ta_xzz_xxxxxz_0[i] * pa_y[i] - ta_xzz_xxxxxz_1[i] * pc_y[i];

        ta_xyzz_xxxxyy_0[i] = 4.0 * ta_yzz_xxxyy_0[i] * fe_0 - 4.0 * ta_yzz_xxxyy_1[i] * fe_0 + ta_yzz_xxxxyy_0[i] * pa_x[i] - ta_yzz_xxxxyy_1[i] * pc_x[i];

        ta_xyzz_xxxxyz_0[i] = 4.0 * ta_yzz_xxxyz_0[i] * fe_0 - 4.0 * ta_yzz_xxxyz_1[i] * fe_0 + ta_yzz_xxxxyz_0[i] * pa_x[i] - ta_yzz_xxxxyz_1[i] * pc_x[i];

        ta_xyzz_xxxxzz_0[i] = ta_xzz_xxxxzz_0[i] * pa_y[i] - ta_xzz_xxxxzz_1[i] * pc_y[i];

        ta_xyzz_xxxyyy_0[i] = 3.0 * ta_yzz_xxyyy_0[i] * fe_0 - 3.0 * ta_yzz_xxyyy_1[i] * fe_0 + ta_yzz_xxxyyy_0[i] * pa_x[i] - ta_yzz_xxxyyy_1[i] * pc_x[i];

        ta_xyzz_xxxyyz_0[i] = 3.0 * ta_yzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yzz_xxyyz_1[i] * fe_0 + ta_yzz_xxxyyz_0[i] * pa_x[i] - ta_yzz_xxxyyz_1[i] * pc_x[i];

        ta_xyzz_xxxyzz_0[i] = 3.0 * ta_yzz_xxyzz_0[i] * fe_0 - 3.0 * ta_yzz_xxyzz_1[i] * fe_0 + ta_yzz_xxxyzz_0[i] * pa_x[i] - ta_yzz_xxxyzz_1[i] * pc_x[i];

        ta_xyzz_xxxzzz_0[i] = ta_xzz_xxxzzz_0[i] * pa_y[i] - ta_xzz_xxxzzz_1[i] * pc_y[i];

        ta_xyzz_xxyyyy_0[i] = 2.0 * ta_yzz_xyyyy_0[i] * fe_0 - 2.0 * ta_yzz_xyyyy_1[i] * fe_0 + ta_yzz_xxyyyy_0[i] * pa_x[i] - ta_yzz_xxyyyy_1[i] * pc_x[i];

        ta_xyzz_xxyyyz_0[i] = 2.0 * ta_yzz_xyyyz_0[i] * fe_0 - 2.0 * ta_yzz_xyyyz_1[i] * fe_0 + ta_yzz_xxyyyz_0[i] * pa_x[i] - ta_yzz_xxyyyz_1[i] * pc_x[i];

        ta_xyzz_xxyyzz_0[i] = 2.0 * ta_yzz_xyyzz_0[i] * fe_0 - 2.0 * ta_yzz_xyyzz_1[i] * fe_0 + ta_yzz_xxyyzz_0[i] * pa_x[i] - ta_yzz_xxyyzz_1[i] * pc_x[i];

        ta_xyzz_xxyzzz_0[i] = 2.0 * ta_yzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzzz_1[i] * fe_0 + ta_yzz_xxyzzz_0[i] * pa_x[i] - ta_yzz_xxyzzz_1[i] * pc_x[i];

        ta_xyzz_xxzzzz_0[i] = ta_xzz_xxzzzz_0[i] * pa_y[i] - ta_xzz_xxzzzz_1[i] * pc_y[i];

        ta_xyzz_xyyyyy_0[i] = ta_yzz_yyyyy_0[i] * fe_0 - ta_yzz_yyyyy_1[i] * fe_0 + ta_yzz_xyyyyy_0[i] * pa_x[i] - ta_yzz_xyyyyy_1[i] * pc_x[i];

        ta_xyzz_xyyyyz_0[i] = ta_yzz_yyyyz_0[i] * fe_0 - ta_yzz_yyyyz_1[i] * fe_0 + ta_yzz_xyyyyz_0[i] * pa_x[i] - ta_yzz_xyyyyz_1[i] * pc_x[i];

        ta_xyzz_xyyyzz_0[i] = ta_yzz_yyyzz_0[i] * fe_0 - ta_yzz_yyyzz_1[i] * fe_0 + ta_yzz_xyyyzz_0[i] * pa_x[i] - ta_yzz_xyyyzz_1[i] * pc_x[i];

        ta_xyzz_xyyzzz_0[i] = ta_yzz_yyzzz_0[i] * fe_0 - ta_yzz_yyzzz_1[i] * fe_0 + ta_yzz_xyyzzz_0[i] * pa_x[i] - ta_yzz_xyyzzz_1[i] * pc_x[i];

        ta_xyzz_xyzzzz_0[i] = ta_yzz_yzzzz_0[i] * fe_0 - ta_yzz_yzzzz_1[i] * fe_0 + ta_yzz_xyzzzz_0[i] * pa_x[i] - ta_yzz_xyzzzz_1[i] * pc_x[i];

        ta_xyzz_xzzzzz_0[i] = ta_xzz_xzzzzz_0[i] * pa_y[i] - ta_xzz_xzzzzz_1[i] * pc_y[i];

        ta_xyzz_yyyyyy_0[i] = ta_yzz_yyyyyy_0[i] * pa_x[i] - ta_yzz_yyyyyy_1[i] * pc_x[i];

        ta_xyzz_yyyyyz_0[i] = ta_yzz_yyyyyz_0[i] * pa_x[i] - ta_yzz_yyyyyz_1[i] * pc_x[i];

        ta_xyzz_yyyyzz_0[i] = ta_yzz_yyyyzz_0[i] * pa_x[i] - ta_yzz_yyyyzz_1[i] * pc_x[i];

        ta_xyzz_yyyzzz_0[i] = ta_yzz_yyyzzz_0[i] * pa_x[i] - ta_yzz_yyyzzz_1[i] * pc_x[i];

        ta_xyzz_yyzzzz_0[i] = ta_yzz_yyzzzz_0[i] * pa_x[i] - ta_yzz_yyzzzz_1[i] * pc_x[i];

        ta_xyzz_yzzzzz_0[i] = ta_yzz_yzzzzz_0[i] * pa_x[i] - ta_yzz_yzzzzz_1[i] * pc_x[i];

        ta_xyzz_zzzzzz_0[i] = ta_yzz_zzzzzz_0[i] * pa_x[i] - ta_yzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 252-280 components of targeted buffer : GI

    auto ta_xzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 252);

    auto ta_xzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 253);

    auto ta_xzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 254);

    auto ta_xzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 255);

    auto ta_xzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 256);

    auto ta_xzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 257);

    auto ta_xzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 258);

    auto ta_xzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 259);

    auto ta_xzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 260);

    auto ta_xzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 261);

    auto ta_xzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 262);

    auto ta_xzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 263);

    auto ta_xzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 264);

    auto ta_xzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 265);

    auto ta_xzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 266);

    auto ta_xzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 267);

    auto ta_xzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 268);

    auto ta_xzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 269);

    auto ta_xzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 270);

    auto ta_xzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 271);

    auto ta_xzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 272);

    auto ta_xzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 273);

    auto ta_xzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 274);

    auto ta_xzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 275);

    auto ta_xzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 276);

    auto ta_xzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 277);

    auto ta_xzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 278);

    auto ta_xzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 279);

    #pragma omp simd aligned(pa_x, pc_x, ta_xzzz_xxxxxx_0, ta_xzzz_xxxxxy_0, ta_xzzz_xxxxxz_0, ta_xzzz_xxxxyy_0, ta_xzzz_xxxxyz_0, ta_xzzz_xxxxzz_0, ta_xzzz_xxxyyy_0, ta_xzzz_xxxyyz_0, ta_xzzz_xxxyzz_0, ta_xzzz_xxxzzz_0, ta_xzzz_xxyyyy_0, ta_xzzz_xxyyyz_0, ta_xzzz_xxyyzz_0, ta_xzzz_xxyzzz_0, ta_xzzz_xxzzzz_0, ta_xzzz_xyyyyy_0, ta_xzzz_xyyyyz_0, ta_xzzz_xyyyzz_0, ta_xzzz_xyyzzz_0, ta_xzzz_xyzzzz_0, ta_xzzz_xzzzzz_0, ta_xzzz_yyyyyy_0, ta_xzzz_yyyyyz_0, ta_xzzz_yyyyzz_0, ta_xzzz_yyyzzz_0, ta_xzzz_yyzzzz_0, ta_xzzz_yzzzzz_0, ta_xzzz_zzzzzz_0, ta_zzz_xxxxx_0, ta_zzz_xxxxx_1, ta_zzz_xxxxxx_0, ta_zzz_xxxxxx_1, ta_zzz_xxxxxy_0, ta_zzz_xxxxxy_1, ta_zzz_xxxxxz_0, ta_zzz_xxxxxz_1, ta_zzz_xxxxy_0, ta_zzz_xxxxy_1, ta_zzz_xxxxyy_0, ta_zzz_xxxxyy_1, ta_zzz_xxxxyz_0, ta_zzz_xxxxyz_1, ta_zzz_xxxxz_0, ta_zzz_xxxxz_1, ta_zzz_xxxxzz_0, ta_zzz_xxxxzz_1, ta_zzz_xxxyy_0, ta_zzz_xxxyy_1, ta_zzz_xxxyyy_0, ta_zzz_xxxyyy_1, ta_zzz_xxxyyz_0, ta_zzz_xxxyyz_1, ta_zzz_xxxyz_0, ta_zzz_xxxyz_1, ta_zzz_xxxyzz_0, ta_zzz_xxxyzz_1, ta_zzz_xxxzz_0, ta_zzz_xxxzz_1, ta_zzz_xxxzzz_0, ta_zzz_xxxzzz_1, ta_zzz_xxyyy_0, ta_zzz_xxyyy_1, ta_zzz_xxyyyy_0, ta_zzz_xxyyyy_1, ta_zzz_xxyyyz_0, ta_zzz_xxyyyz_1, ta_zzz_xxyyz_0, ta_zzz_xxyyz_1, ta_zzz_xxyyzz_0, ta_zzz_xxyyzz_1, ta_zzz_xxyzz_0, ta_zzz_xxyzz_1, ta_zzz_xxyzzz_0, ta_zzz_xxyzzz_1, ta_zzz_xxzzz_0, ta_zzz_xxzzz_1, ta_zzz_xxzzzz_0, ta_zzz_xxzzzz_1, ta_zzz_xyyyy_0, ta_zzz_xyyyy_1, ta_zzz_xyyyyy_0, ta_zzz_xyyyyy_1, ta_zzz_xyyyyz_0, ta_zzz_xyyyyz_1, ta_zzz_xyyyz_0, ta_zzz_xyyyz_1, ta_zzz_xyyyzz_0, ta_zzz_xyyyzz_1, ta_zzz_xyyzz_0, ta_zzz_xyyzz_1, ta_zzz_xyyzzz_0, ta_zzz_xyyzzz_1, ta_zzz_xyzzz_0, ta_zzz_xyzzz_1, ta_zzz_xyzzzz_0, ta_zzz_xyzzzz_1, ta_zzz_xzzzz_0, ta_zzz_xzzzz_1, ta_zzz_xzzzzz_0, ta_zzz_xzzzzz_1, ta_zzz_yyyyy_0, ta_zzz_yyyyy_1, ta_zzz_yyyyyy_0, ta_zzz_yyyyyy_1, ta_zzz_yyyyyz_0, ta_zzz_yyyyyz_1, ta_zzz_yyyyz_0, ta_zzz_yyyyz_1, ta_zzz_yyyyzz_0, ta_zzz_yyyyzz_1, ta_zzz_yyyzz_0, ta_zzz_yyyzz_1, ta_zzz_yyyzzz_0, ta_zzz_yyyzzz_1, ta_zzz_yyzzz_0, ta_zzz_yyzzz_1, ta_zzz_yyzzzz_0, ta_zzz_yyzzzz_1, ta_zzz_yzzzz_0, ta_zzz_yzzzz_1, ta_zzz_yzzzzz_0, ta_zzz_yzzzzz_1, ta_zzz_zzzzz_0, ta_zzz_zzzzz_1, ta_zzz_zzzzzz_0, ta_zzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzz_xxxxxx_0[i] = 6.0 * ta_zzz_xxxxx_0[i] * fe_0 - 6.0 * ta_zzz_xxxxx_1[i] * fe_0 + ta_zzz_xxxxxx_0[i] * pa_x[i] - ta_zzz_xxxxxx_1[i] * pc_x[i];

        ta_xzzz_xxxxxy_0[i] = 5.0 * ta_zzz_xxxxy_0[i] * fe_0 - 5.0 * ta_zzz_xxxxy_1[i] * fe_0 + ta_zzz_xxxxxy_0[i] * pa_x[i] - ta_zzz_xxxxxy_1[i] * pc_x[i];

        ta_xzzz_xxxxxz_0[i] = 5.0 * ta_zzz_xxxxz_0[i] * fe_0 - 5.0 * ta_zzz_xxxxz_1[i] * fe_0 + ta_zzz_xxxxxz_0[i] * pa_x[i] - ta_zzz_xxxxxz_1[i] * pc_x[i];

        ta_xzzz_xxxxyy_0[i] = 4.0 * ta_zzz_xxxyy_0[i] * fe_0 - 4.0 * ta_zzz_xxxyy_1[i] * fe_0 + ta_zzz_xxxxyy_0[i] * pa_x[i] - ta_zzz_xxxxyy_1[i] * pc_x[i];

        ta_xzzz_xxxxyz_0[i] = 4.0 * ta_zzz_xxxyz_0[i] * fe_0 - 4.0 * ta_zzz_xxxyz_1[i] * fe_0 + ta_zzz_xxxxyz_0[i] * pa_x[i] - ta_zzz_xxxxyz_1[i] * pc_x[i];

        ta_xzzz_xxxxzz_0[i] = 4.0 * ta_zzz_xxxzz_0[i] * fe_0 - 4.0 * ta_zzz_xxxzz_1[i] * fe_0 + ta_zzz_xxxxzz_0[i] * pa_x[i] - ta_zzz_xxxxzz_1[i] * pc_x[i];

        ta_xzzz_xxxyyy_0[i] = 3.0 * ta_zzz_xxyyy_0[i] * fe_0 - 3.0 * ta_zzz_xxyyy_1[i] * fe_0 + ta_zzz_xxxyyy_0[i] * pa_x[i] - ta_zzz_xxxyyy_1[i] * pc_x[i];

        ta_xzzz_xxxyyz_0[i] = 3.0 * ta_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta_zzz_xxyyz_1[i] * fe_0 + ta_zzz_xxxyyz_0[i] * pa_x[i] - ta_zzz_xxxyyz_1[i] * pc_x[i];

        ta_xzzz_xxxyzz_0[i] = 3.0 * ta_zzz_xxyzz_0[i] * fe_0 - 3.0 * ta_zzz_xxyzz_1[i] * fe_0 + ta_zzz_xxxyzz_0[i] * pa_x[i] - ta_zzz_xxxyzz_1[i] * pc_x[i];

        ta_xzzz_xxxzzz_0[i] = 3.0 * ta_zzz_xxzzz_0[i] * fe_0 - 3.0 * ta_zzz_xxzzz_1[i] * fe_0 + ta_zzz_xxxzzz_0[i] * pa_x[i] - ta_zzz_xxxzzz_1[i] * pc_x[i];

        ta_xzzz_xxyyyy_0[i] = 2.0 * ta_zzz_xyyyy_0[i] * fe_0 - 2.0 * ta_zzz_xyyyy_1[i] * fe_0 + ta_zzz_xxyyyy_0[i] * pa_x[i] - ta_zzz_xxyyyy_1[i] * pc_x[i];

        ta_xzzz_xxyyyz_0[i] = 2.0 * ta_zzz_xyyyz_0[i] * fe_0 - 2.0 * ta_zzz_xyyyz_1[i] * fe_0 + ta_zzz_xxyyyz_0[i] * pa_x[i] - ta_zzz_xxyyyz_1[i] * pc_x[i];

        ta_xzzz_xxyyzz_0[i] = 2.0 * ta_zzz_xyyzz_0[i] * fe_0 - 2.0 * ta_zzz_xyyzz_1[i] * fe_0 + ta_zzz_xxyyzz_0[i] * pa_x[i] - ta_zzz_xxyyzz_1[i] * pc_x[i];

        ta_xzzz_xxyzzz_0[i] = 2.0 * ta_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta_zzz_xyzzz_1[i] * fe_0 + ta_zzz_xxyzzz_0[i] * pa_x[i] - ta_zzz_xxyzzz_1[i] * pc_x[i];

        ta_xzzz_xxzzzz_0[i] = 2.0 * ta_zzz_xzzzz_0[i] * fe_0 - 2.0 * ta_zzz_xzzzz_1[i] * fe_0 + ta_zzz_xxzzzz_0[i] * pa_x[i] - ta_zzz_xxzzzz_1[i] * pc_x[i];

        ta_xzzz_xyyyyy_0[i] = ta_zzz_yyyyy_0[i] * fe_0 - ta_zzz_yyyyy_1[i] * fe_0 + ta_zzz_xyyyyy_0[i] * pa_x[i] - ta_zzz_xyyyyy_1[i] * pc_x[i];

        ta_xzzz_xyyyyz_0[i] = ta_zzz_yyyyz_0[i] * fe_0 - ta_zzz_yyyyz_1[i] * fe_0 + ta_zzz_xyyyyz_0[i] * pa_x[i] - ta_zzz_xyyyyz_1[i] * pc_x[i];

        ta_xzzz_xyyyzz_0[i] = ta_zzz_yyyzz_0[i] * fe_0 - ta_zzz_yyyzz_1[i] * fe_0 + ta_zzz_xyyyzz_0[i] * pa_x[i] - ta_zzz_xyyyzz_1[i] * pc_x[i];

        ta_xzzz_xyyzzz_0[i] = ta_zzz_yyzzz_0[i] * fe_0 - ta_zzz_yyzzz_1[i] * fe_0 + ta_zzz_xyyzzz_0[i] * pa_x[i] - ta_zzz_xyyzzz_1[i] * pc_x[i];

        ta_xzzz_xyzzzz_0[i] = ta_zzz_yzzzz_0[i] * fe_0 - ta_zzz_yzzzz_1[i] * fe_0 + ta_zzz_xyzzzz_0[i] * pa_x[i] - ta_zzz_xyzzzz_1[i] * pc_x[i];

        ta_xzzz_xzzzzz_0[i] = ta_zzz_zzzzz_0[i] * fe_0 - ta_zzz_zzzzz_1[i] * fe_0 + ta_zzz_xzzzzz_0[i] * pa_x[i] - ta_zzz_xzzzzz_1[i] * pc_x[i];

        ta_xzzz_yyyyyy_0[i] = ta_zzz_yyyyyy_0[i] * pa_x[i] - ta_zzz_yyyyyy_1[i] * pc_x[i];

        ta_xzzz_yyyyyz_0[i] = ta_zzz_yyyyyz_0[i] * pa_x[i] - ta_zzz_yyyyyz_1[i] * pc_x[i];

        ta_xzzz_yyyyzz_0[i] = ta_zzz_yyyyzz_0[i] * pa_x[i] - ta_zzz_yyyyzz_1[i] * pc_x[i];

        ta_xzzz_yyyzzz_0[i] = ta_zzz_yyyzzz_0[i] * pa_x[i] - ta_zzz_yyyzzz_1[i] * pc_x[i];

        ta_xzzz_yyzzzz_0[i] = ta_zzz_yyzzzz_0[i] * pa_x[i] - ta_zzz_yyzzzz_1[i] * pc_x[i];

        ta_xzzz_yzzzzz_0[i] = ta_zzz_yzzzzz_0[i] * pa_x[i] - ta_zzz_yzzzzz_1[i] * pc_x[i];

        ta_xzzz_zzzzzz_0[i] = ta_zzz_zzzzzz_0[i] * pa_x[i] - ta_zzz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 280-308 components of targeted buffer : GI

    auto ta_yyyy_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 280);

    auto ta_yyyy_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 281);

    auto ta_yyyy_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 282);

    auto ta_yyyy_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 283);

    auto ta_yyyy_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 284);

    auto ta_yyyy_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 285);

    auto ta_yyyy_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 286);

    auto ta_yyyy_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 287);

    auto ta_yyyy_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 288);

    auto ta_yyyy_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 289);

    auto ta_yyyy_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 290);

    auto ta_yyyy_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 291);

    auto ta_yyyy_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 292);

    auto ta_yyyy_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 293);

    auto ta_yyyy_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 294);

    auto ta_yyyy_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 295);

    auto ta_yyyy_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 296);

    auto ta_yyyy_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 297);

    auto ta_yyyy_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 298);

    auto ta_yyyy_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 299);

    auto ta_yyyy_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 300);

    auto ta_yyyy_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 301);

    auto ta_yyyy_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 302);

    auto ta_yyyy_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 303);

    auto ta_yyyy_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 304);

    auto ta_yyyy_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 305);

    auto ta_yyyy_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 306);

    auto ta_yyyy_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 307);

    #pragma omp simd aligned(pa_y, pc_y, ta_yy_xxxxxx_0, ta_yy_xxxxxx_1, ta_yy_xxxxxy_0, ta_yy_xxxxxy_1, ta_yy_xxxxxz_0, ta_yy_xxxxxz_1, ta_yy_xxxxyy_0, ta_yy_xxxxyy_1, ta_yy_xxxxyz_0, ta_yy_xxxxyz_1, ta_yy_xxxxzz_0, ta_yy_xxxxzz_1, ta_yy_xxxyyy_0, ta_yy_xxxyyy_1, ta_yy_xxxyyz_0, ta_yy_xxxyyz_1, ta_yy_xxxyzz_0, ta_yy_xxxyzz_1, ta_yy_xxxzzz_0, ta_yy_xxxzzz_1, ta_yy_xxyyyy_0, ta_yy_xxyyyy_1, ta_yy_xxyyyz_0, ta_yy_xxyyyz_1, ta_yy_xxyyzz_0, ta_yy_xxyyzz_1, ta_yy_xxyzzz_0, ta_yy_xxyzzz_1, ta_yy_xxzzzz_0, ta_yy_xxzzzz_1, ta_yy_xyyyyy_0, ta_yy_xyyyyy_1, ta_yy_xyyyyz_0, ta_yy_xyyyyz_1, ta_yy_xyyyzz_0, ta_yy_xyyyzz_1, ta_yy_xyyzzz_0, ta_yy_xyyzzz_1, ta_yy_xyzzzz_0, ta_yy_xyzzzz_1, ta_yy_xzzzzz_0, ta_yy_xzzzzz_1, ta_yy_yyyyyy_0, ta_yy_yyyyyy_1, ta_yy_yyyyyz_0, ta_yy_yyyyyz_1, ta_yy_yyyyzz_0, ta_yy_yyyyzz_1, ta_yy_yyyzzz_0, ta_yy_yyyzzz_1, ta_yy_yyzzzz_0, ta_yy_yyzzzz_1, ta_yy_yzzzzz_0, ta_yy_yzzzzz_1, ta_yy_zzzzzz_0, ta_yy_zzzzzz_1, ta_yyy_xxxxx_0, ta_yyy_xxxxx_1, ta_yyy_xxxxxx_0, ta_yyy_xxxxxx_1, ta_yyy_xxxxxy_0, ta_yyy_xxxxxy_1, ta_yyy_xxxxxz_0, ta_yyy_xxxxxz_1, ta_yyy_xxxxy_0, ta_yyy_xxxxy_1, ta_yyy_xxxxyy_0, ta_yyy_xxxxyy_1, ta_yyy_xxxxyz_0, ta_yyy_xxxxyz_1, ta_yyy_xxxxz_0, ta_yyy_xxxxz_1, ta_yyy_xxxxzz_0, ta_yyy_xxxxzz_1, ta_yyy_xxxyy_0, ta_yyy_xxxyy_1, ta_yyy_xxxyyy_0, ta_yyy_xxxyyy_1, ta_yyy_xxxyyz_0, ta_yyy_xxxyyz_1, ta_yyy_xxxyz_0, ta_yyy_xxxyz_1, ta_yyy_xxxyzz_0, ta_yyy_xxxyzz_1, ta_yyy_xxxzz_0, ta_yyy_xxxzz_1, ta_yyy_xxxzzz_0, ta_yyy_xxxzzz_1, ta_yyy_xxyyy_0, ta_yyy_xxyyy_1, ta_yyy_xxyyyy_0, ta_yyy_xxyyyy_1, ta_yyy_xxyyyz_0, ta_yyy_xxyyyz_1, ta_yyy_xxyyz_0, ta_yyy_xxyyz_1, ta_yyy_xxyyzz_0, ta_yyy_xxyyzz_1, ta_yyy_xxyzz_0, ta_yyy_xxyzz_1, ta_yyy_xxyzzz_0, ta_yyy_xxyzzz_1, ta_yyy_xxzzz_0, ta_yyy_xxzzz_1, ta_yyy_xxzzzz_0, ta_yyy_xxzzzz_1, ta_yyy_xyyyy_0, ta_yyy_xyyyy_1, ta_yyy_xyyyyy_0, ta_yyy_xyyyyy_1, ta_yyy_xyyyyz_0, ta_yyy_xyyyyz_1, ta_yyy_xyyyz_0, ta_yyy_xyyyz_1, ta_yyy_xyyyzz_0, ta_yyy_xyyyzz_1, ta_yyy_xyyzz_0, ta_yyy_xyyzz_1, ta_yyy_xyyzzz_0, ta_yyy_xyyzzz_1, ta_yyy_xyzzz_0, ta_yyy_xyzzz_1, ta_yyy_xyzzzz_0, ta_yyy_xyzzzz_1, ta_yyy_xzzzz_0, ta_yyy_xzzzz_1, ta_yyy_xzzzzz_0, ta_yyy_xzzzzz_1, ta_yyy_yyyyy_0, ta_yyy_yyyyy_1, ta_yyy_yyyyyy_0, ta_yyy_yyyyyy_1, ta_yyy_yyyyyz_0, ta_yyy_yyyyyz_1, ta_yyy_yyyyz_0, ta_yyy_yyyyz_1, ta_yyy_yyyyzz_0, ta_yyy_yyyyzz_1, ta_yyy_yyyzz_0, ta_yyy_yyyzz_1, ta_yyy_yyyzzz_0, ta_yyy_yyyzzz_1, ta_yyy_yyzzz_0, ta_yyy_yyzzz_1, ta_yyy_yyzzzz_0, ta_yyy_yyzzzz_1, ta_yyy_yzzzz_0, ta_yyy_yzzzz_1, ta_yyy_yzzzzz_0, ta_yyy_yzzzzz_1, ta_yyy_zzzzz_0, ta_yyy_zzzzz_1, ta_yyy_zzzzzz_0, ta_yyy_zzzzzz_1, ta_yyyy_xxxxxx_0, ta_yyyy_xxxxxy_0, ta_yyyy_xxxxxz_0, ta_yyyy_xxxxyy_0, ta_yyyy_xxxxyz_0, ta_yyyy_xxxxzz_0, ta_yyyy_xxxyyy_0, ta_yyyy_xxxyyz_0, ta_yyyy_xxxyzz_0, ta_yyyy_xxxzzz_0, ta_yyyy_xxyyyy_0, ta_yyyy_xxyyyz_0, ta_yyyy_xxyyzz_0, ta_yyyy_xxyzzz_0, ta_yyyy_xxzzzz_0, ta_yyyy_xyyyyy_0, ta_yyyy_xyyyyz_0, ta_yyyy_xyyyzz_0, ta_yyyy_xyyzzz_0, ta_yyyy_xyzzzz_0, ta_yyyy_xzzzzz_0, ta_yyyy_yyyyyy_0, ta_yyyy_yyyyyz_0, ta_yyyy_yyyyzz_0, ta_yyyy_yyyzzz_0, ta_yyyy_yyzzzz_0, ta_yyyy_yzzzzz_0, ta_yyyy_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyy_xxxxxx_0[i] = 3.0 * ta_yy_xxxxxx_0[i] * fe_0 - 3.0 * ta_yy_xxxxxx_1[i] * fe_0 + ta_yyy_xxxxxx_0[i] * pa_y[i] - ta_yyy_xxxxxx_1[i] * pc_y[i];

        ta_yyyy_xxxxxy_0[i] = 3.0 * ta_yy_xxxxxy_0[i] * fe_0 - 3.0 * ta_yy_xxxxxy_1[i] * fe_0 + ta_yyy_xxxxx_0[i] * fe_0 - ta_yyy_xxxxx_1[i] * fe_0 + ta_yyy_xxxxxy_0[i] * pa_y[i] - ta_yyy_xxxxxy_1[i] * pc_y[i];

        ta_yyyy_xxxxxz_0[i] = 3.0 * ta_yy_xxxxxz_0[i] * fe_0 - 3.0 * ta_yy_xxxxxz_1[i] * fe_0 + ta_yyy_xxxxxz_0[i] * pa_y[i] - ta_yyy_xxxxxz_1[i] * pc_y[i];

        ta_yyyy_xxxxyy_0[i] = 3.0 * ta_yy_xxxxyy_0[i] * fe_0 - 3.0 * ta_yy_xxxxyy_1[i] * fe_0 + 2.0 * ta_yyy_xxxxy_0[i] * fe_0 - 2.0 * ta_yyy_xxxxy_1[i] * fe_0 + ta_yyy_xxxxyy_0[i] * pa_y[i] - ta_yyy_xxxxyy_1[i] * pc_y[i];

        ta_yyyy_xxxxyz_0[i] = 3.0 * ta_yy_xxxxyz_0[i] * fe_0 - 3.0 * ta_yy_xxxxyz_1[i] * fe_0 + ta_yyy_xxxxz_0[i] * fe_0 - ta_yyy_xxxxz_1[i] * fe_0 + ta_yyy_xxxxyz_0[i] * pa_y[i] - ta_yyy_xxxxyz_1[i] * pc_y[i];

        ta_yyyy_xxxxzz_0[i] = 3.0 * ta_yy_xxxxzz_0[i] * fe_0 - 3.0 * ta_yy_xxxxzz_1[i] * fe_0 + ta_yyy_xxxxzz_0[i] * pa_y[i] - ta_yyy_xxxxzz_1[i] * pc_y[i];

        ta_yyyy_xxxyyy_0[i] = 3.0 * ta_yy_xxxyyy_0[i] * fe_0 - 3.0 * ta_yy_xxxyyy_1[i] * fe_0 + 3.0 * ta_yyy_xxxyy_0[i] * fe_0 - 3.0 * ta_yyy_xxxyy_1[i] * fe_0 + ta_yyy_xxxyyy_0[i] * pa_y[i] - ta_yyy_xxxyyy_1[i] * pc_y[i];

        ta_yyyy_xxxyyz_0[i] = 3.0 * ta_yy_xxxyyz_0[i] * fe_0 - 3.0 * ta_yy_xxxyyz_1[i] * fe_0 + 2.0 * ta_yyy_xxxyz_0[i] * fe_0 - 2.0 * ta_yyy_xxxyz_1[i] * fe_0 + ta_yyy_xxxyyz_0[i] * pa_y[i] - ta_yyy_xxxyyz_1[i] * pc_y[i];

        ta_yyyy_xxxyzz_0[i] = 3.0 * ta_yy_xxxyzz_0[i] * fe_0 - 3.0 * ta_yy_xxxyzz_1[i] * fe_0 + ta_yyy_xxxzz_0[i] * fe_0 - ta_yyy_xxxzz_1[i] * fe_0 + ta_yyy_xxxyzz_0[i] * pa_y[i] - ta_yyy_xxxyzz_1[i] * pc_y[i];

        ta_yyyy_xxxzzz_0[i] = 3.0 * ta_yy_xxxzzz_0[i] * fe_0 - 3.0 * ta_yy_xxxzzz_1[i] * fe_0 + ta_yyy_xxxzzz_0[i] * pa_y[i] - ta_yyy_xxxzzz_1[i] * pc_y[i];

        ta_yyyy_xxyyyy_0[i] = 3.0 * ta_yy_xxyyyy_0[i] * fe_0 - 3.0 * ta_yy_xxyyyy_1[i] * fe_0 + 4.0 * ta_yyy_xxyyy_0[i] * fe_0 - 4.0 * ta_yyy_xxyyy_1[i] * fe_0 + ta_yyy_xxyyyy_0[i] * pa_y[i] - ta_yyy_xxyyyy_1[i] * pc_y[i];

        ta_yyyy_xxyyyz_0[i] = 3.0 * ta_yy_xxyyyz_0[i] * fe_0 - 3.0 * ta_yy_xxyyyz_1[i] * fe_0 + 3.0 * ta_yyy_xxyyz_0[i] * fe_0 - 3.0 * ta_yyy_xxyyz_1[i] * fe_0 + ta_yyy_xxyyyz_0[i] * pa_y[i] - ta_yyy_xxyyyz_1[i] * pc_y[i];

        ta_yyyy_xxyyzz_0[i] = 3.0 * ta_yy_xxyyzz_0[i] * fe_0 - 3.0 * ta_yy_xxyyzz_1[i] * fe_0 + 2.0 * ta_yyy_xxyzz_0[i] * fe_0 - 2.0 * ta_yyy_xxyzz_1[i] * fe_0 + ta_yyy_xxyyzz_0[i] * pa_y[i] - ta_yyy_xxyyzz_1[i] * pc_y[i];

        ta_yyyy_xxyzzz_0[i] = 3.0 * ta_yy_xxyzzz_0[i] * fe_0 - 3.0 * ta_yy_xxyzzz_1[i] * fe_0 + ta_yyy_xxzzz_0[i] * fe_0 - ta_yyy_xxzzz_1[i] * fe_0 + ta_yyy_xxyzzz_0[i] * pa_y[i] - ta_yyy_xxyzzz_1[i] * pc_y[i];

        ta_yyyy_xxzzzz_0[i] = 3.0 * ta_yy_xxzzzz_0[i] * fe_0 - 3.0 * ta_yy_xxzzzz_1[i] * fe_0 + ta_yyy_xxzzzz_0[i] * pa_y[i] - ta_yyy_xxzzzz_1[i] * pc_y[i];

        ta_yyyy_xyyyyy_0[i] = 3.0 * ta_yy_xyyyyy_0[i] * fe_0 - 3.0 * ta_yy_xyyyyy_1[i] * fe_0 + 5.0 * ta_yyy_xyyyy_0[i] * fe_0 - 5.0 * ta_yyy_xyyyy_1[i] * fe_0 + ta_yyy_xyyyyy_0[i] * pa_y[i] - ta_yyy_xyyyyy_1[i] * pc_y[i];

        ta_yyyy_xyyyyz_0[i] = 3.0 * ta_yy_xyyyyz_0[i] * fe_0 - 3.0 * ta_yy_xyyyyz_1[i] * fe_0 + 4.0 * ta_yyy_xyyyz_0[i] * fe_0 - 4.0 * ta_yyy_xyyyz_1[i] * fe_0 + ta_yyy_xyyyyz_0[i] * pa_y[i] - ta_yyy_xyyyyz_1[i] * pc_y[i];

        ta_yyyy_xyyyzz_0[i] = 3.0 * ta_yy_xyyyzz_0[i] * fe_0 - 3.0 * ta_yy_xyyyzz_1[i] * fe_0 + 3.0 * ta_yyy_xyyzz_0[i] * fe_0 - 3.0 * ta_yyy_xyyzz_1[i] * fe_0 + ta_yyy_xyyyzz_0[i] * pa_y[i] - ta_yyy_xyyyzz_1[i] * pc_y[i];

        ta_yyyy_xyyzzz_0[i] = 3.0 * ta_yy_xyyzzz_0[i] * fe_0 - 3.0 * ta_yy_xyyzzz_1[i] * fe_0 + 2.0 * ta_yyy_xyzzz_0[i] * fe_0 - 2.0 * ta_yyy_xyzzz_1[i] * fe_0 + ta_yyy_xyyzzz_0[i] * pa_y[i] - ta_yyy_xyyzzz_1[i] * pc_y[i];

        ta_yyyy_xyzzzz_0[i] = 3.0 * ta_yy_xyzzzz_0[i] * fe_0 - 3.0 * ta_yy_xyzzzz_1[i] * fe_0 + ta_yyy_xzzzz_0[i] * fe_0 - ta_yyy_xzzzz_1[i] * fe_0 + ta_yyy_xyzzzz_0[i] * pa_y[i] - ta_yyy_xyzzzz_1[i] * pc_y[i];

        ta_yyyy_xzzzzz_0[i] = 3.0 * ta_yy_xzzzzz_0[i] * fe_0 - 3.0 * ta_yy_xzzzzz_1[i] * fe_0 + ta_yyy_xzzzzz_0[i] * pa_y[i] - ta_yyy_xzzzzz_1[i] * pc_y[i];

        ta_yyyy_yyyyyy_0[i] = 3.0 * ta_yy_yyyyyy_0[i] * fe_0 - 3.0 * ta_yy_yyyyyy_1[i] * fe_0 + 6.0 * ta_yyy_yyyyy_0[i] * fe_0 - 6.0 * ta_yyy_yyyyy_1[i] * fe_0 + ta_yyy_yyyyyy_0[i] * pa_y[i] - ta_yyy_yyyyyy_1[i] * pc_y[i];

        ta_yyyy_yyyyyz_0[i] = 3.0 * ta_yy_yyyyyz_0[i] * fe_0 - 3.0 * ta_yy_yyyyyz_1[i] * fe_0 + 5.0 * ta_yyy_yyyyz_0[i] * fe_0 - 5.0 * ta_yyy_yyyyz_1[i] * fe_0 + ta_yyy_yyyyyz_0[i] * pa_y[i] - ta_yyy_yyyyyz_1[i] * pc_y[i];

        ta_yyyy_yyyyzz_0[i] = 3.0 * ta_yy_yyyyzz_0[i] * fe_0 - 3.0 * ta_yy_yyyyzz_1[i] * fe_0 + 4.0 * ta_yyy_yyyzz_0[i] * fe_0 - 4.0 * ta_yyy_yyyzz_1[i] * fe_0 + ta_yyy_yyyyzz_0[i] * pa_y[i] - ta_yyy_yyyyzz_1[i] * pc_y[i];

        ta_yyyy_yyyzzz_0[i] = 3.0 * ta_yy_yyyzzz_0[i] * fe_0 - 3.0 * ta_yy_yyyzzz_1[i] * fe_0 + 3.0 * ta_yyy_yyzzz_0[i] * fe_0 - 3.0 * ta_yyy_yyzzz_1[i] * fe_0 + ta_yyy_yyyzzz_0[i] * pa_y[i] - ta_yyy_yyyzzz_1[i] * pc_y[i];

        ta_yyyy_yyzzzz_0[i] = 3.0 * ta_yy_yyzzzz_0[i] * fe_0 - 3.0 * ta_yy_yyzzzz_1[i] * fe_0 + 2.0 * ta_yyy_yzzzz_0[i] * fe_0 - 2.0 * ta_yyy_yzzzz_1[i] * fe_0 + ta_yyy_yyzzzz_0[i] * pa_y[i] - ta_yyy_yyzzzz_1[i] * pc_y[i];

        ta_yyyy_yzzzzz_0[i] = 3.0 * ta_yy_yzzzzz_0[i] * fe_0 - 3.0 * ta_yy_yzzzzz_1[i] * fe_0 + ta_yyy_zzzzz_0[i] * fe_0 - ta_yyy_zzzzz_1[i] * fe_0 + ta_yyy_yzzzzz_0[i] * pa_y[i] - ta_yyy_yzzzzz_1[i] * pc_y[i];

        ta_yyyy_zzzzzz_0[i] = 3.0 * ta_yy_zzzzzz_0[i] * fe_0 - 3.0 * ta_yy_zzzzzz_1[i] * fe_0 + ta_yyy_zzzzzz_0[i] * pa_y[i] - ta_yyy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 308-336 components of targeted buffer : GI

    auto ta_yyyz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 308);

    auto ta_yyyz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 309);

    auto ta_yyyz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 310);

    auto ta_yyyz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 311);

    auto ta_yyyz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 312);

    auto ta_yyyz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 313);

    auto ta_yyyz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 314);

    auto ta_yyyz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 315);

    auto ta_yyyz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 316);

    auto ta_yyyz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 317);

    auto ta_yyyz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 318);

    auto ta_yyyz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 319);

    auto ta_yyyz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 320);

    auto ta_yyyz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 321);

    auto ta_yyyz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 322);

    auto ta_yyyz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 323);

    auto ta_yyyz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 324);

    auto ta_yyyz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 325);

    auto ta_yyyz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 326);

    auto ta_yyyz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 327);

    auto ta_yyyz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 328);

    auto ta_yyyz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 329);

    auto ta_yyyz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 330);

    auto ta_yyyz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 331);

    auto ta_yyyz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 332);

    auto ta_yyyz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 333);

    auto ta_yyyz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 334);

    auto ta_yyyz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 335);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yyy_xxxxxx_0, ta_yyy_xxxxxx_1, ta_yyy_xxxxxy_0, ta_yyy_xxxxxy_1, ta_yyy_xxxxy_0, ta_yyy_xxxxy_1, ta_yyy_xxxxyy_0, ta_yyy_xxxxyy_1, ta_yyy_xxxxyz_0, ta_yyy_xxxxyz_1, ta_yyy_xxxyy_0, ta_yyy_xxxyy_1, ta_yyy_xxxyyy_0, ta_yyy_xxxyyy_1, ta_yyy_xxxyyz_0, ta_yyy_xxxyyz_1, ta_yyy_xxxyz_0, ta_yyy_xxxyz_1, ta_yyy_xxxyzz_0, ta_yyy_xxxyzz_1, ta_yyy_xxyyy_0, ta_yyy_xxyyy_1, ta_yyy_xxyyyy_0, ta_yyy_xxyyyy_1, ta_yyy_xxyyyz_0, ta_yyy_xxyyyz_1, ta_yyy_xxyyz_0, ta_yyy_xxyyz_1, ta_yyy_xxyyzz_0, ta_yyy_xxyyzz_1, ta_yyy_xxyzz_0, ta_yyy_xxyzz_1, ta_yyy_xxyzzz_0, ta_yyy_xxyzzz_1, ta_yyy_xyyyy_0, ta_yyy_xyyyy_1, ta_yyy_xyyyyy_0, ta_yyy_xyyyyy_1, ta_yyy_xyyyyz_0, ta_yyy_xyyyyz_1, ta_yyy_xyyyz_0, ta_yyy_xyyyz_1, ta_yyy_xyyyzz_0, ta_yyy_xyyyzz_1, ta_yyy_xyyzz_0, ta_yyy_xyyzz_1, ta_yyy_xyyzzz_0, ta_yyy_xyyzzz_1, ta_yyy_xyzzz_0, ta_yyy_xyzzz_1, ta_yyy_xyzzzz_0, ta_yyy_xyzzzz_1, ta_yyy_yyyyy_0, ta_yyy_yyyyy_1, ta_yyy_yyyyyy_0, ta_yyy_yyyyyy_1, ta_yyy_yyyyyz_0, ta_yyy_yyyyyz_1, ta_yyy_yyyyz_0, ta_yyy_yyyyz_1, ta_yyy_yyyyzz_0, ta_yyy_yyyyzz_1, ta_yyy_yyyzz_0, ta_yyy_yyyzz_1, ta_yyy_yyyzzz_0, ta_yyy_yyyzzz_1, ta_yyy_yyzzz_0, ta_yyy_yyzzz_1, ta_yyy_yyzzzz_0, ta_yyy_yyzzzz_1, ta_yyy_yzzzz_0, ta_yyy_yzzzz_1, ta_yyy_yzzzzz_0, ta_yyy_yzzzzz_1, ta_yyyz_xxxxxx_0, ta_yyyz_xxxxxy_0, ta_yyyz_xxxxxz_0, ta_yyyz_xxxxyy_0, ta_yyyz_xxxxyz_0, ta_yyyz_xxxxzz_0, ta_yyyz_xxxyyy_0, ta_yyyz_xxxyyz_0, ta_yyyz_xxxyzz_0, ta_yyyz_xxxzzz_0, ta_yyyz_xxyyyy_0, ta_yyyz_xxyyyz_0, ta_yyyz_xxyyzz_0, ta_yyyz_xxyzzz_0, ta_yyyz_xxzzzz_0, ta_yyyz_xyyyyy_0, ta_yyyz_xyyyyz_0, ta_yyyz_xyyyzz_0, ta_yyyz_xyyzzz_0, ta_yyyz_xyzzzz_0, ta_yyyz_xzzzzz_0, ta_yyyz_yyyyyy_0, ta_yyyz_yyyyyz_0, ta_yyyz_yyyyzz_0, ta_yyyz_yyyzzz_0, ta_yyyz_yyzzzz_0, ta_yyyz_yzzzzz_0, ta_yyyz_zzzzzz_0, ta_yyz_xxxxxz_0, ta_yyz_xxxxxz_1, ta_yyz_xxxxzz_0, ta_yyz_xxxxzz_1, ta_yyz_xxxzzz_0, ta_yyz_xxxzzz_1, ta_yyz_xxzzzz_0, ta_yyz_xxzzzz_1, ta_yyz_xzzzzz_0, ta_yyz_xzzzzz_1, ta_yyz_zzzzzz_0, ta_yyz_zzzzzz_1, ta_yz_xxxxxz_0, ta_yz_xxxxxz_1, ta_yz_xxxxzz_0, ta_yz_xxxxzz_1, ta_yz_xxxzzz_0, ta_yz_xxxzzz_1, ta_yz_xxzzzz_0, ta_yz_xxzzzz_1, ta_yz_xzzzzz_0, ta_yz_xzzzzz_1, ta_yz_zzzzzz_0, ta_yz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyz_xxxxxx_0[i] = ta_yyy_xxxxxx_0[i] * pa_z[i] - ta_yyy_xxxxxx_1[i] * pc_z[i];

        ta_yyyz_xxxxxy_0[i] = ta_yyy_xxxxxy_0[i] * pa_z[i] - ta_yyy_xxxxxy_1[i] * pc_z[i];

        ta_yyyz_xxxxxz_0[i] = 2.0 * ta_yz_xxxxxz_0[i] * fe_0 - 2.0 * ta_yz_xxxxxz_1[i] * fe_0 + ta_yyz_xxxxxz_0[i] * pa_y[i] - ta_yyz_xxxxxz_1[i] * pc_y[i];

        ta_yyyz_xxxxyy_0[i] = ta_yyy_xxxxyy_0[i] * pa_z[i] - ta_yyy_xxxxyy_1[i] * pc_z[i];

        ta_yyyz_xxxxyz_0[i] = ta_yyy_xxxxy_0[i] * fe_0 - ta_yyy_xxxxy_1[i] * fe_0 + ta_yyy_xxxxyz_0[i] * pa_z[i] - ta_yyy_xxxxyz_1[i] * pc_z[i];

        ta_yyyz_xxxxzz_0[i] = 2.0 * ta_yz_xxxxzz_0[i] * fe_0 - 2.0 * ta_yz_xxxxzz_1[i] * fe_0 + ta_yyz_xxxxzz_0[i] * pa_y[i] - ta_yyz_xxxxzz_1[i] * pc_y[i];

        ta_yyyz_xxxyyy_0[i] = ta_yyy_xxxyyy_0[i] * pa_z[i] - ta_yyy_xxxyyy_1[i] * pc_z[i];

        ta_yyyz_xxxyyz_0[i] = ta_yyy_xxxyy_0[i] * fe_0 - ta_yyy_xxxyy_1[i] * fe_0 + ta_yyy_xxxyyz_0[i] * pa_z[i] - ta_yyy_xxxyyz_1[i] * pc_z[i];

        ta_yyyz_xxxyzz_0[i] = 2.0 * ta_yyy_xxxyz_0[i] * fe_0 - 2.0 * ta_yyy_xxxyz_1[i] * fe_0 + ta_yyy_xxxyzz_0[i] * pa_z[i] - ta_yyy_xxxyzz_1[i] * pc_z[i];

        ta_yyyz_xxxzzz_0[i] = 2.0 * ta_yz_xxxzzz_0[i] * fe_0 - 2.0 * ta_yz_xxxzzz_1[i] * fe_0 + ta_yyz_xxxzzz_0[i] * pa_y[i] - ta_yyz_xxxzzz_1[i] * pc_y[i];

        ta_yyyz_xxyyyy_0[i] = ta_yyy_xxyyyy_0[i] * pa_z[i] - ta_yyy_xxyyyy_1[i] * pc_z[i];

        ta_yyyz_xxyyyz_0[i] = ta_yyy_xxyyy_0[i] * fe_0 - ta_yyy_xxyyy_1[i] * fe_0 + ta_yyy_xxyyyz_0[i] * pa_z[i] - ta_yyy_xxyyyz_1[i] * pc_z[i];

        ta_yyyz_xxyyzz_0[i] = 2.0 * ta_yyy_xxyyz_0[i] * fe_0 - 2.0 * ta_yyy_xxyyz_1[i] * fe_0 + ta_yyy_xxyyzz_0[i] * pa_z[i] - ta_yyy_xxyyzz_1[i] * pc_z[i];

        ta_yyyz_xxyzzz_0[i] = 3.0 * ta_yyy_xxyzz_0[i] * fe_0 - 3.0 * ta_yyy_xxyzz_1[i] * fe_0 + ta_yyy_xxyzzz_0[i] * pa_z[i] - ta_yyy_xxyzzz_1[i] * pc_z[i];

        ta_yyyz_xxzzzz_0[i] = 2.0 * ta_yz_xxzzzz_0[i] * fe_0 - 2.0 * ta_yz_xxzzzz_1[i] * fe_0 + ta_yyz_xxzzzz_0[i] * pa_y[i] - ta_yyz_xxzzzz_1[i] * pc_y[i];

        ta_yyyz_xyyyyy_0[i] = ta_yyy_xyyyyy_0[i] * pa_z[i] - ta_yyy_xyyyyy_1[i] * pc_z[i];

        ta_yyyz_xyyyyz_0[i] = ta_yyy_xyyyy_0[i] * fe_0 - ta_yyy_xyyyy_1[i] * fe_0 + ta_yyy_xyyyyz_0[i] * pa_z[i] - ta_yyy_xyyyyz_1[i] * pc_z[i];

        ta_yyyz_xyyyzz_0[i] = 2.0 * ta_yyy_xyyyz_0[i] * fe_0 - 2.0 * ta_yyy_xyyyz_1[i] * fe_0 + ta_yyy_xyyyzz_0[i] * pa_z[i] - ta_yyy_xyyyzz_1[i] * pc_z[i];

        ta_yyyz_xyyzzz_0[i] = 3.0 * ta_yyy_xyyzz_0[i] * fe_0 - 3.0 * ta_yyy_xyyzz_1[i] * fe_0 + ta_yyy_xyyzzz_0[i] * pa_z[i] - ta_yyy_xyyzzz_1[i] * pc_z[i];

        ta_yyyz_xyzzzz_0[i] = 4.0 * ta_yyy_xyzzz_0[i] * fe_0 - 4.0 * ta_yyy_xyzzz_1[i] * fe_0 + ta_yyy_xyzzzz_0[i] * pa_z[i] - ta_yyy_xyzzzz_1[i] * pc_z[i];

        ta_yyyz_xzzzzz_0[i] = 2.0 * ta_yz_xzzzzz_0[i] * fe_0 - 2.0 * ta_yz_xzzzzz_1[i] * fe_0 + ta_yyz_xzzzzz_0[i] * pa_y[i] - ta_yyz_xzzzzz_1[i] * pc_y[i];

        ta_yyyz_yyyyyy_0[i] = ta_yyy_yyyyyy_0[i] * pa_z[i] - ta_yyy_yyyyyy_1[i] * pc_z[i];

        ta_yyyz_yyyyyz_0[i] = ta_yyy_yyyyy_0[i] * fe_0 - ta_yyy_yyyyy_1[i] * fe_0 + ta_yyy_yyyyyz_0[i] * pa_z[i] - ta_yyy_yyyyyz_1[i] * pc_z[i];

        ta_yyyz_yyyyzz_0[i] = 2.0 * ta_yyy_yyyyz_0[i] * fe_0 - 2.0 * ta_yyy_yyyyz_1[i] * fe_0 + ta_yyy_yyyyzz_0[i] * pa_z[i] - ta_yyy_yyyyzz_1[i] * pc_z[i];

        ta_yyyz_yyyzzz_0[i] = 3.0 * ta_yyy_yyyzz_0[i] * fe_0 - 3.0 * ta_yyy_yyyzz_1[i] * fe_0 + ta_yyy_yyyzzz_0[i] * pa_z[i] - ta_yyy_yyyzzz_1[i] * pc_z[i];

        ta_yyyz_yyzzzz_0[i] = 4.0 * ta_yyy_yyzzz_0[i] * fe_0 - 4.0 * ta_yyy_yyzzz_1[i] * fe_0 + ta_yyy_yyzzzz_0[i] * pa_z[i] - ta_yyy_yyzzzz_1[i] * pc_z[i];

        ta_yyyz_yzzzzz_0[i] = 5.0 * ta_yyy_yzzzz_0[i] * fe_0 - 5.0 * ta_yyy_yzzzz_1[i] * fe_0 + ta_yyy_yzzzzz_0[i] * pa_z[i] - ta_yyy_yzzzzz_1[i] * pc_z[i];

        ta_yyyz_zzzzzz_0[i] = 2.0 * ta_yz_zzzzzz_0[i] * fe_0 - 2.0 * ta_yz_zzzzzz_1[i] * fe_0 + ta_yyz_zzzzzz_0[i] * pa_y[i] - ta_yyz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 336-364 components of targeted buffer : GI

    auto ta_yyzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 336);

    auto ta_yyzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 337);

    auto ta_yyzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 338);

    auto ta_yyzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 339);

    auto ta_yyzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 340);

    auto ta_yyzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 341);

    auto ta_yyzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 342);

    auto ta_yyzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 343);

    auto ta_yyzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 344);

    auto ta_yyzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 345);

    auto ta_yyzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 346);

    auto ta_yyzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 347);

    auto ta_yyzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 348);

    auto ta_yyzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 349);

    auto ta_yyzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 350);

    auto ta_yyzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 351);

    auto ta_yyzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 352);

    auto ta_yyzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 353);

    auto ta_yyzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 354);

    auto ta_yyzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 355);

    auto ta_yyzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 356);

    auto ta_yyzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 357);

    auto ta_yyzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 358);

    auto ta_yyzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 359);

    auto ta_yyzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 360);

    auto ta_yyzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 361);

    auto ta_yyzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 362);

    auto ta_yyzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 363);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta_yy_xxxxxy_0, ta_yy_xxxxxy_1, ta_yy_xxxxyy_0, ta_yy_xxxxyy_1, ta_yy_xxxyyy_0, ta_yy_xxxyyy_1, ta_yy_xxyyyy_0, ta_yy_xxyyyy_1, ta_yy_xyyyyy_0, ta_yy_xyyyyy_1, ta_yy_yyyyyy_0, ta_yy_yyyyyy_1, ta_yyz_xxxxxy_0, ta_yyz_xxxxxy_1, ta_yyz_xxxxyy_0, ta_yyz_xxxxyy_1, ta_yyz_xxxyyy_0, ta_yyz_xxxyyy_1, ta_yyz_xxyyyy_0, ta_yyz_xxyyyy_1, ta_yyz_xyyyyy_0, ta_yyz_xyyyyy_1, ta_yyz_yyyyyy_0, ta_yyz_yyyyyy_1, ta_yyzz_xxxxxx_0, ta_yyzz_xxxxxy_0, ta_yyzz_xxxxxz_0, ta_yyzz_xxxxyy_0, ta_yyzz_xxxxyz_0, ta_yyzz_xxxxzz_0, ta_yyzz_xxxyyy_0, ta_yyzz_xxxyyz_0, ta_yyzz_xxxyzz_0, ta_yyzz_xxxzzz_0, ta_yyzz_xxyyyy_0, ta_yyzz_xxyyyz_0, ta_yyzz_xxyyzz_0, ta_yyzz_xxyzzz_0, ta_yyzz_xxzzzz_0, ta_yyzz_xyyyyy_0, ta_yyzz_xyyyyz_0, ta_yyzz_xyyyzz_0, ta_yyzz_xyyzzz_0, ta_yyzz_xyzzzz_0, ta_yyzz_xzzzzz_0, ta_yyzz_yyyyyy_0, ta_yyzz_yyyyyz_0, ta_yyzz_yyyyzz_0, ta_yyzz_yyyzzz_0, ta_yyzz_yyzzzz_0, ta_yyzz_yzzzzz_0, ta_yyzz_zzzzzz_0, ta_yzz_xxxxxx_0, ta_yzz_xxxxxx_1, ta_yzz_xxxxxz_0, ta_yzz_xxxxxz_1, ta_yzz_xxxxyz_0, ta_yzz_xxxxyz_1, ta_yzz_xxxxz_0, ta_yzz_xxxxz_1, ta_yzz_xxxxzz_0, ta_yzz_xxxxzz_1, ta_yzz_xxxyyz_0, ta_yzz_xxxyyz_1, ta_yzz_xxxyz_0, ta_yzz_xxxyz_1, ta_yzz_xxxyzz_0, ta_yzz_xxxyzz_1, ta_yzz_xxxzz_0, ta_yzz_xxxzz_1, ta_yzz_xxxzzz_0, ta_yzz_xxxzzz_1, ta_yzz_xxyyyz_0, ta_yzz_xxyyyz_1, ta_yzz_xxyyz_0, ta_yzz_xxyyz_1, ta_yzz_xxyyzz_0, ta_yzz_xxyyzz_1, ta_yzz_xxyzz_0, ta_yzz_xxyzz_1, ta_yzz_xxyzzz_0, ta_yzz_xxyzzz_1, ta_yzz_xxzzz_0, ta_yzz_xxzzz_1, ta_yzz_xxzzzz_0, ta_yzz_xxzzzz_1, ta_yzz_xyyyyz_0, ta_yzz_xyyyyz_1, ta_yzz_xyyyz_0, ta_yzz_xyyyz_1, ta_yzz_xyyyzz_0, ta_yzz_xyyyzz_1, ta_yzz_xyyzz_0, ta_yzz_xyyzz_1, ta_yzz_xyyzzz_0, ta_yzz_xyyzzz_1, ta_yzz_xyzzz_0, ta_yzz_xyzzz_1, ta_yzz_xyzzzz_0, ta_yzz_xyzzzz_1, ta_yzz_xzzzz_0, ta_yzz_xzzzz_1, ta_yzz_xzzzzz_0, ta_yzz_xzzzzz_1, ta_yzz_yyyyyz_0, ta_yzz_yyyyyz_1, ta_yzz_yyyyz_0, ta_yzz_yyyyz_1, ta_yzz_yyyyzz_0, ta_yzz_yyyyzz_1, ta_yzz_yyyzz_0, ta_yzz_yyyzz_1, ta_yzz_yyyzzz_0, ta_yzz_yyyzzz_1, ta_yzz_yyzzz_0, ta_yzz_yyzzz_1, ta_yzz_yyzzzz_0, ta_yzz_yyzzzz_1, ta_yzz_yzzzz_0, ta_yzz_yzzzz_1, ta_yzz_yzzzzz_0, ta_yzz_yzzzzz_1, ta_yzz_zzzzz_0, ta_yzz_zzzzz_1, ta_yzz_zzzzzz_0, ta_yzz_zzzzzz_1, ta_zz_xxxxxx_0, ta_zz_xxxxxx_1, ta_zz_xxxxxz_0, ta_zz_xxxxxz_1, ta_zz_xxxxyz_0, ta_zz_xxxxyz_1, ta_zz_xxxxzz_0, ta_zz_xxxxzz_1, ta_zz_xxxyyz_0, ta_zz_xxxyyz_1, ta_zz_xxxyzz_0, ta_zz_xxxyzz_1, ta_zz_xxxzzz_0, ta_zz_xxxzzz_1, ta_zz_xxyyyz_0, ta_zz_xxyyyz_1, ta_zz_xxyyzz_0, ta_zz_xxyyzz_1, ta_zz_xxyzzz_0, ta_zz_xxyzzz_1, ta_zz_xxzzzz_0, ta_zz_xxzzzz_1, ta_zz_xyyyyz_0, ta_zz_xyyyyz_1, ta_zz_xyyyzz_0, ta_zz_xyyyzz_1, ta_zz_xyyzzz_0, ta_zz_xyyzzz_1, ta_zz_xyzzzz_0, ta_zz_xyzzzz_1, ta_zz_xzzzzz_0, ta_zz_xzzzzz_1, ta_zz_yyyyyz_0, ta_zz_yyyyyz_1, ta_zz_yyyyzz_0, ta_zz_yyyyzz_1, ta_zz_yyyzzz_0, ta_zz_yyyzzz_1, ta_zz_yyzzzz_0, ta_zz_yyzzzz_1, ta_zz_yzzzzz_0, ta_zz_yzzzzz_1, ta_zz_zzzzzz_0, ta_zz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzz_xxxxxx_0[i] = ta_zz_xxxxxx_0[i] * fe_0 - ta_zz_xxxxxx_1[i] * fe_0 + ta_yzz_xxxxxx_0[i] * pa_y[i] - ta_yzz_xxxxxx_1[i] * pc_y[i];

        ta_yyzz_xxxxxy_0[i] = ta_yy_xxxxxy_0[i] * fe_0 - ta_yy_xxxxxy_1[i] * fe_0 + ta_yyz_xxxxxy_0[i] * pa_z[i] - ta_yyz_xxxxxy_1[i] * pc_z[i];

        ta_yyzz_xxxxxz_0[i] = ta_zz_xxxxxz_0[i] * fe_0 - ta_zz_xxxxxz_1[i] * fe_0 + ta_yzz_xxxxxz_0[i] * pa_y[i] - ta_yzz_xxxxxz_1[i] * pc_y[i];

        ta_yyzz_xxxxyy_0[i] = ta_yy_xxxxyy_0[i] * fe_0 - ta_yy_xxxxyy_1[i] * fe_0 + ta_yyz_xxxxyy_0[i] * pa_z[i] - ta_yyz_xxxxyy_1[i] * pc_z[i];

        ta_yyzz_xxxxyz_0[i] = ta_zz_xxxxyz_0[i] * fe_0 - ta_zz_xxxxyz_1[i] * fe_0 + ta_yzz_xxxxz_0[i] * fe_0 - ta_yzz_xxxxz_1[i] * fe_0 + ta_yzz_xxxxyz_0[i] * pa_y[i] - ta_yzz_xxxxyz_1[i] * pc_y[i];

        ta_yyzz_xxxxzz_0[i] = ta_zz_xxxxzz_0[i] * fe_0 - ta_zz_xxxxzz_1[i] * fe_0 + ta_yzz_xxxxzz_0[i] * pa_y[i] - ta_yzz_xxxxzz_1[i] * pc_y[i];

        ta_yyzz_xxxyyy_0[i] = ta_yy_xxxyyy_0[i] * fe_0 - ta_yy_xxxyyy_1[i] * fe_0 + ta_yyz_xxxyyy_0[i] * pa_z[i] - ta_yyz_xxxyyy_1[i] * pc_z[i];

        ta_yyzz_xxxyyz_0[i] = ta_zz_xxxyyz_0[i] * fe_0 - ta_zz_xxxyyz_1[i] * fe_0 + 2.0 * ta_yzz_xxxyz_0[i] * fe_0 - 2.0 * ta_yzz_xxxyz_1[i] * fe_0 + ta_yzz_xxxyyz_0[i] * pa_y[i] - ta_yzz_xxxyyz_1[i] * pc_y[i];

        ta_yyzz_xxxyzz_0[i] = ta_zz_xxxyzz_0[i] * fe_0 - ta_zz_xxxyzz_1[i] * fe_0 + ta_yzz_xxxzz_0[i] * fe_0 - ta_yzz_xxxzz_1[i] * fe_0 + ta_yzz_xxxyzz_0[i] * pa_y[i] - ta_yzz_xxxyzz_1[i] * pc_y[i];

        ta_yyzz_xxxzzz_0[i] = ta_zz_xxxzzz_0[i] * fe_0 - ta_zz_xxxzzz_1[i] * fe_0 + ta_yzz_xxxzzz_0[i] * pa_y[i] - ta_yzz_xxxzzz_1[i] * pc_y[i];

        ta_yyzz_xxyyyy_0[i] = ta_yy_xxyyyy_0[i] * fe_0 - ta_yy_xxyyyy_1[i] * fe_0 + ta_yyz_xxyyyy_0[i] * pa_z[i] - ta_yyz_xxyyyy_1[i] * pc_z[i];

        ta_yyzz_xxyyyz_0[i] = ta_zz_xxyyyz_0[i] * fe_0 - ta_zz_xxyyyz_1[i] * fe_0 + 3.0 * ta_yzz_xxyyz_0[i] * fe_0 - 3.0 * ta_yzz_xxyyz_1[i] * fe_0 + ta_yzz_xxyyyz_0[i] * pa_y[i] - ta_yzz_xxyyyz_1[i] * pc_y[i];

        ta_yyzz_xxyyzz_0[i] = ta_zz_xxyyzz_0[i] * fe_0 - ta_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta_yzz_xxyzz_0[i] * fe_0 - 2.0 * ta_yzz_xxyzz_1[i] * fe_0 + ta_yzz_xxyyzz_0[i] * pa_y[i] - ta_yzz_xxyyzz_1[i] * pc_y[i];

        ta_yyzz_xxyzzz_0[i] = ta_zz_xxyzzz_0[i] * fe_0 - ta_zz_xxyzzz_1[i] * fe_0 + ta_yzz_xxzzz_0[i] * fe_0 - ta_yzz_xxzzz_1[i] * fe_0 + ta_yzz_xxyzzz_0[i] * pa_y[i] - ta_yzz_xxyzzz_1[i] * pc_y[i];

        ta_yyzz_xxzzzz_0[i] = ta_zz_xxzzzz_0[i] * fe_0 - ta_zz_xxzzzz_1[i] * fe_0 + ta_yzz_xxzzzz_0[i] * pa_y[i] - ta_yzz_xxzzzz_1[i] * pc_y[i];

        ta_yyzz_xyyyyy_0[i] = ta_yy_xyyyyy_0[i] * fe_0 - ta_yy_xyyyyy_1[i] * fe_0 + ta_yyz_xyyyyy_0[i] * pa_z[i] - ta_yyz_xyyyyy_1[i] * pc_z[i];

        ta_yyzz_xyyyyz_0[i] = ta_zz_xyyyyz_0[i] * fe_0 - ta_zz_xyyyyz_1[i] * fe_0 + 4.0 * ta_yzz_xyyyz_0[i] * fe_0 - 4.0 * ta_yzz_xyyyz_1[i] * fe_0 + ta_yzz_xyyyyz_0[i] * pa_y[i] - ta_yzz_xyyyyz_1[i] * pc_y[i];

        ta_yyzz_xyyyzz_0[i] = ta_zz_xyyyzz_0[i] * fe_0 - ta_zz_xyyyzz_1[i] * fe_0 + 3.0 * ta_yzz_xyyzz_0[i] * fe_0 - 3.0 * ta_yzz_xyyzz_1[i] * fe_0 + ta_yzz_xyyyzz_0[i] * pa_y[i] - ta_yzz_xyyyzz_1[i] * pc_y[i];

        ta_yyzz_xyyzzz_0[i] = ta_zz_xyyzzz_0[i] * fe_0 - ta_zz_xyyzzz_1[i] * fe_0 + 2.0 * ta_yzz_xyzzz_0[i] * fe_0 - 2.0 * ta_yzz_xyzzz_1[i] * fe_0 + ta_yzz_xyyzzz_0[i] * pa_y[i] - ta_yzz_xyyzzz_1[i] * pc_y[i];

        ta_yyzz_xyzzzz_0[i] = ta_zz_xyzzzz_0[i] * fe_0 - ta_zz_xyzzzz_1[i] * fe_0 + ta_yzz_xzzzz_0[i] * fe_0 - ta_yzz_xzzzz_1[i] * fe_0 + ta_yzz_xyzzzz_0[i] * pa_y[i] - ta_yzz_xyzzzz_1[i] * pc_y[i];

        ta_yyzz_xzzzzz_0[i] = ta_zz_xzzzzz_0[i] * fe_0 - ta_zz_xzzzzz_1[i] * fe_0 + ta_yzz_xzzzzz_0[i] * pa_y[i] - ta_yzz_xzzzzz_1[i] * pc_y[i];

        ta_yyzz_yyyyyy_0[i] = ta_yy_yyyyyy_0[i] * fe_0 - ta_yy_yyyyyy_1[i] * fe_0 + ta_yyz_yyyyyy_0[i] * pa_z[i] - ta_yyz_yyyyyy_1[i] * pc_z[i];

        ta_yyzz_yyyyyz_0[i] = ta_zz_yyyyyz_0[i] * fe_0 - ta_zz_yyyyyz_1[i] * fe_0 + 5.0 * ta_yzz_yyyyz_0[i] * fe_0 - 5.0 * ta_yzz_yyyyz_1[i] * fe_0 + ta_yzz_yyyyyz_0[i] * pa_y[i] - ta_yzz_yyyyyz_1[i] * pc_y[i];

        ta_yyzz_yyyyzz_0[i] = ta_zz_yyyyzz_0[i] * fe_0 - ta_zz_yyyyzz_1[i] * fe_0 + 4.0 * ta_yzz_yyyzz_0[i] * fe_0 - 4.0 * ta_yzz_yyyzz_1[i] * fe_0 + ta_yzz_yyyyzz_0[i] * pa_y[i] - ta_yzz_yyyyzz_1[i] * pc_y[i];

        ta_yyzz_yyyzzz_0[i] = ta_zz_yyyzzz_0[i] * fe_0 - ta_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta_yzz_yyzzz_0[i] * fe_0 - 3.0 * ta_yzz_yyzzz_1[i] * fe_0 + ta_yzz_yyyzzz_0[i] * pa_y[i] - ta_yzz_yyyzzz_1[i] * pc_y[i];

        ta_yyzz_yyzzzz_0[i] = ta_zz_yyzzzz_0[i] * fe_0 - ta_zz_yyzzzz_1[i] * fe_0 + 2.0 * ta_yzz_yzzzz_0[i] * fe_0 - 2.0 * ta_yzz_yzzzz_1[i] * fe_0 + ta_yzz_yyzzzz_0[i] * pa_y[i] - ta_yzz_yyzzzz_1[i] * pc_y[i];

        ta_yyzz_yzzzzz_0[i] = ta_zz_yzzzzz_0[i] * fe_0 - ta_zz_yzzzzz_1[i] * fe_0 + ta_yzz_zzzzz_0[i] * fe_0 - ta_yzz_zzzzz_1[i] * fe_0 + ta_yzz_yzzzzz_0[i] * pa_y[i] - ta_yzz_yzzzzz_1[i] * pc_y[i];

        ta_yyzz_zzzzzz_0[i] = ta_zz_zzzzzz_0[i] * fe_0 - ta_zz_zzzzzz_1[i] * fe_0 + ta_yzz_zzzzzz_0[i] * pa_y[i] - ta_yzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 364-392 components of targeted buffer : GI

    auto ta_yzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 364);

    auto ta_yzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 365);

    auto ta_yzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 366);

    auto ta_yzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 367);

    auto ta_yzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 368);

    auto ta_yzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 369);

    auto ta_yzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 370);

    auto ta_yzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 371);

    auto ta_yzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 372);

    auto ta_yzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 373);

    auto ta_yzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 374);

    auto ta_yzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 375);

    auto ta_yzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 376);

    auto ta_yzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 377);

    auto ta_yzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 378);

    auto ta_yzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 379);

    auto ta_yzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 380);

    auto ta_yzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 381);

    auto ta_yzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 382);

    auto ta_yzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 383);

    auto ta_yzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 384);

    auto ta_yzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 385);

    auto ta_yzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 386);

    auto ta_yzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 387);

    auto ta_yzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 388);

    auto ta_yzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 389);

    auto ta_yzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 390);

    auto ta_yzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 391);

    #pragma omp simd aligned(pa_y, pc_y, ta_yzzz_xxxxxx_0, ta_yzzz_xxxxxy_0, ta_yzzz_xxxxxz_0, ta_yzzz_xxxxyy_0, ta_yzzz_xxxxyz_0, ta_yzzz_xxxxzz_0, ta_yzzz_xxxyyy_0, ta_yzzz_xxxyyz_0, ta_yzzz_xxxyzz_0, ta_yzzz_xxxzzz_0, ta_yzzz_xxyyyy_0, ta_yzzz_xxyyyz_0, ta_yzzz_xxyyzz_0, ta_yzzz_xxyzzz_0, ta_yzzz_xxzzzz_0, ta_yzzz_xyyyyy_0, ta_yzzz_xyyyyz_0, ta_yzzz_xyyyzz_0, ta_yzzz_xyyzzz_0, ta_yzzz_xyzzzz_0, ta_yzzz_xzzzzz_0, ta_yzzz_yyyyyy_0, ta_yzzz_yyyyyz_0, ta_yzzz_yyyyzz_0, ta_yzzz_yyyzzz_0, ta_yzzz_yyzzzz_0, ta_yzzz_yzzzzz_0, ta_yzzz_zzzzzz_0, ta_zzz_xxxxx_0, ta_zzz_xxxxx_1, ta_zzz_xxxxxx_0, ta_zzz_xxxxxx_1, ta_zzz_xxxxxy_0, ta_zzz_xxxxxy_1, ta_zzz_xxxxxz_0, ta_zzz_xxxxxz_1, ta_zzz_xxxxy_0, ta_zzz_xxxxy_1, ta_zzz_xxxxyy_0, ta_zzz_xxxxyy_1, ta_zzz_xxxxyz_0, ta_zzz_xxxxyz_1, ta_zzz_xxxxz_0, ta_zzz_xxxxz_1, ta_zzz_xxxxzz_0, ta_zzz_xxxxzz_1, ta_zzz_xxxyy_0, ta_zzz_xxxyy_1, ta_zzz_xxxyyy_0, ta_zzz_xxxyyy_1, ta_zzz_xxxyyz_0, ta_zzz_xxxyyz_1, ta_zzz_xxxyz_0, ta_zzz_xxxyz_1, ta_zzz_xxxyzz_0, ta_zzz_xxxyzz_1, ta_zzz_xxxzz_0, ta_zzz_xxxzz_1, ta_zzz_xxxzzz_0, ta_zzz_xxxzzz_1, ta_zzz_xxyyy_0, ta_zzz_xxyyy_1, ta_zzz_xxyyyy_0, ta_zzz_xxyyyy_1, ta_zzz_xxyyyz_0, ta_zzz_xxyyyz_1, ta_zzz_xxyyz_0, ta_zzz_xxyyz_1, ta_zzz_xxyyzz_0, ta_zzz_xxyyzz_1, ta_zzz_xxyzz_0, ta_zzz_xxyzz_1, ta_zzz_xxyzzz_0, ta_zzz_xxyzzz_1, ta_zzz_xxzzz_0, ta_zzz_xxzzz_1, ta_zzz_xxzzzz_0, ta_zzz_xxzzzz_1, ta_zzz_xyyyy_0, ta_zzz_xyyyy_1, ta_zzz_xyyyyy_0, ta_zzz_xyyyyy_1, ta_zzz_xyyyyz_0, ta_zzz_xyyyyz_1, ta_zzz_xyyyz_0, ta_zzz_xyyyz_1, ta_zzz_xyyyzz_0, ta_zzz_xyyyzz_1, ta_zzz_xyyzz_0, ta_zzz_xyyzz_1, ta_zzz_xyyzzz_0, ta_zzz_xyyzzz_1, ta_zzz_xyzzz_0, ta_zzz_xyzzz_1, ta_zzz_xyzzzz_0, ta_zzz_xyzzzz_1, ta_zzz_xzzzz_0, ta_zzz_xzzzz_1, ta_zzz_xzzzzz_0, ta_zzz_xzzzzz_1, ta_zzz_yyyyy_0, ta_zzz_yyyyy_1, ta_zzz_yyyyyy_0, ta_zzz_yyyyyy_1, ta_zzz_yyyyyz_0, ta_zzz_yyyyyz_1, ta_zzz_yyyyz_0, ta_zzz_yyyyz_1, ta_zzz_yyyyzz_0, ta_zzz_yyyyzz_1, ta_zzz_yyyzz_0, ta_zzz_yyyzz_1, ta_zzz_yyyzzz_0, ta_zzz_yyyzzz_1, ta_zzz_yyzzz_0, ta_zzz_yyzzz_1, ta_zzz_yyzzzz_0, ta_zzz_yyzzzz_1, ta_zzz_yzzzz_0, ta_zzz_yzzzz_1, ta_zzz_yzzzzz_0, ta_zzz_yzzzzz_1, ta_zzz_zzzzz_0, ta_zzz_zzzzz_1, ta_zzz_zzzzzz_0, ta_zzz_zzzzzz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzz_xxxxxx_0[i] = ta_zzz_xxxxxx_0[i] * pa_y[i] - ta_zzz_xxxxxx_1[i] * pc_y[i];

        ta_yzzz_xxxxxy_0[i] = ta_zzz_xxxxx_0[i] * fe_0 - ta_zzz_xxxxx_1[i] * fe_0 + ta_zzz_xxxxxy_0[i] * pa_y[i] - ta_zzz_xxxxxy_1[i] * pc_y[i];

        ta_yzzz_xxxxxz_0[i] = ta_zzz_xxxxxz_0[i] * pa_y[i] - ta_zzz_xxxxxz_1[i] * pc_y[i];

        ta_yzzz_xxxxyy_0[i] = 2.0 * ta_zzz_xxxxy_0[i] * fe_0 - 2.0 * ta_zzz_xxxxy_1[i] * fe_0 + ta_zzz_xxxxyy_0[i] * pa_y[i] - ta_zzz_xxxxyy_1[i] * pc_y[i];

        ta_yzzz_xxxxyz_0[i] = ta_zzz_xxxxz_0[i] * fe_0 - ta_zzz_xxxxz_1[i] * fe_0 + ta_zzz_xxxxyz_0[i] * pa_y[i] - ta_zzz_xxxxyz_1[i] * pc_y[i];

        ta_yzzz_xxxxzz_0[i] = ta_zzz_xxxxzz_0[i] * pa_y[i] - ta_zzz_xxxxzz_1[i] * pc_y[i];

        ta_yzzz_xxxyyy_0[i] = 3.0 * ta_zzz_xxxyy_0[i] * fe_0 - 3.0 * ta_zzz_xxxyy_1[i] * fe_0 + ta_zzz_xxxyyy_0[i] * pa_y[i] - ta_zzz_xxxyyy_1[i] * pc_y[i];

        ta_yzzz_xxxyyz_0[i] = 2.0 * ta_zzz_xxxyz_0[i] * fe_0 - 2.0 * ta_zzz_xxxyz_1[i] * fe_0 + ta_zzz_xxxyyz_0[i] * pa_y[i] - ta_zzz_xxxyyz_1[i] * pc_y[i];

        ta_yzzz_xxxyzz_0[i] = ta_zzz_xxxzz_0[i] * fe_0 - ta_zzz_xxxzz_1[i] * fe_0 + ta_zzz_xxxyzz_0[i] * pa_y[i] - ta_zzz_xxxyzz_1[i] * pc_y[i];

        ta_yzzz_xxxzzz_0[i] = ta_zzz_xxxzzz_0[i] * pa_y[i] - ta_zzz_xxxzzz_1[i] * pc_y[i];

        ta_yzzz_xxyyyy_0[i] = 4.0 * ta_zzz_xxyyy_0[i] * fe_0 - 4.0 * ta_zzz_xxyyy_1[i] * fe_0 + ta_zzz_xxyyyy_0[i] * pa_y[i] - ta_zzz_xxyyyy_1[i] * pc_y[i];

        ta_yzzz_xxyyyz_0[i] = 3.0 * ta_zzz_xxyyz_0[i] * fe_0 - 3.0 * ta_zzz_xxyyz_1[i] * fe_0 + ta_zzz_xxyyyz_0[i] * pa_y[i] - ta_zzz_xxyyyz_1[i] * pc_y[i];

        ta_yzzz_xxyyzz_0[i] = 2.0 * ta_zzz_xxyzz_0[i] * fe_0 - 2.0 * ta_zzz_xxyzz_1[i] * fe_0 + ta_zzz_xxyyzz_0[i] * pa_y[i] - ta_zzz_xxyyzz_1[i] * pc_y[i];

        ta_yzzz_xxyzzz_0[i] = ta_zzz_xxzzz_0[i] * fe_0 - ta_zzz_xxzzz_1[i] * fe_0 + ta_zzz_xxyzzz_0[i] * pa_y[i] - ta_zzz_xxyzzz_1[i] * pc_y[i];

        ta_yzzz_xxzzzz_0[i] = ta_zzz_xxzzzz_0[i] * pa_y[i] - ta_zzz_xxzzzz_1[i] * pc_y[i];

        ta_yzzz_xyyyyy_0[i] = 5.0 * ta_zzz_xyyyy_0[i] * fe_0 - 5.0 * ta_zzz_xyyyy_1[i] * fe_0 + ta_zzz_xyyyyy_0[i] * pa_y[i] - ta_zzz_xyyyyy_1[i] * pc_y[i];

        ta_yzzz_xyyyyz_0[i] = 4.0 * ta_zzz_xyyyz_0[i] * fe_0 - 4.0 * ta_zzz_xyyyz_1[i] * fe_0 + ta_zzz_xyyyyz_0[i] * pa_y[i] - ta_zzz_xyyyyz_1[i] * pc_y[i];

        ta_yzzz_xyyyzz_0[i] = 3.0 * ta_zzz_xyyzz_0[i] * fe_0 - 3.0 * ta_zzz_xyyzz_1[i] * fe_0 + ta_zzz_xyyyzz_0[i] * pa_y[i] - ta_zzz_xyyyzz_1[i] * pc_y[i];

        ta_yzzz_xyyzzz_0[i] = 2.0 * ta_zzz_xyzzz_0[i] * fe_0 - 2.0 * ta_zzz_xyzzz_1[i] * fe_0 + ta_zzz_xyyzzz_0[i] * pa_y[i] - ta_zzz_xyyzzz_1[i] * pc_y[i];

        ta_yzzz_xyzzzz_0[i] = ta_zzz_xzzzz_0[i] * fe_0 - ta_zzz_xzzzz_1[i] * fe_0 + ta_zzz_xyzzzz_0[i] * pa_y[i] - ta_zzz_xyzzzz_1[i] * pc_y[i];

        ta_yzzz_xzzzzz_0[i] = ta_zzz_xzzzzz_0[i] * pa_y[i] - ta_zzz_xzzzzz_1[i] * pc_y[i];

        ta_yzzz_yyyyyy_0[i] = 6.0 * ta_zzz_yyyyy_0[i] * fe_0 - 6.0 * ta_zzz_yyyyy_1[i] * fe_0 + ta_zzz_yyyyyy_0[i] * pa_y[i] - ta_zzz_yyyyyy_1[i] * pc_y[i];

        ta_yzzz_yyyyyz_0[i] = 5.0 * ta_zzz_yyyyz_0[i] * fe_0 - 5.0 * ta_zzz_yyyyz_1[i] * fe_0 + ta_zzz_yyyyyz_0[i] * pa_y[i] - ta_zzz_yyyyyz_1[i] * pc_y[i];

        ta_yzzz_yyyyzz_0[i] = 4.0 * ta_zzz_yyyzz_0[i] * fe_0 - 4.0 * ta_zzz_yyyzz_1[i] * fe_0 + ta_zzz_yyyyzz_0[i] * pa_y[i] - ta_zzz_yyyyzz_1[i] * pc_y[i];

        ta_yzzz_yyyzzz_0[i] = 3.0 * ta_zzz_yyzzz_0[i] * fe_0 - 3.0 * ta_zzz_yyzzz_1[i] * fe_0 + ta_zzz_yyyzzz_0[i] * pa_y[i] - ta_zzz_yyyzzz_1[i] * pc_y[i];

        ta_yzzz_yyzzzz_0[i] = 2.0 * ta_zzz_yzzzz_0[i] * fe_0 - 2.0 * ta_zzz_yzzzz_1[i] * fe_0 + ta_zzz_yyzzzz_0[i] * pa_y[i] - ta_zzz_yyzzzz_1[i] * pc_y[i];

        ta_yzzz_yzzzzz_0[i] = ta_zzz_zzzzz_0[i] * fe_0 - ta_zzz_zzzzz_1[i] * fe_0 + ta_zzz_yzzzzz_0[i] * pa_y[i] - ta_zzz_yzzzzz_1[i] * pc_y[i];

        ta_yzzz_zzzzzz_0[i] = ta_zzz_zzzzzz_0[i] * pa_y[i] - ta_zzz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 392-420 components of targeted buffer : GI

    auto ta_zzzz_xxxxxx_0 = pbuffer.data(idx_npot_0_gi + 392);

    auto ta_zzzz_xxxxxy_0 = pbuffer.data(idx_npot_0_gi + 393);

    auto ta_zzzz_xxxxxz_0 = pbuffer.data(idx_npot_0_gi + 394);

    auto ta_zzzz_xxxxyy_0 = pbuffer.data(idx_npot_0_gi + 395);

    auto ta_zzzz_xxxxyz_0 = pbuffer.data(idx_npot_0_gi + 396);

    auto ta_zzzz_xxxxzz_0 = pbuffer.data(idx_npot_0_gi + 397);

    auto ta_zzzz_xxxyyy_0 = pbuffer.data(idx_npot_0_gi + 398);

    auto ta_zzzz_xxxyyz_0 = pbuffer.data(idx_npot_0_gi + 399);

    auto ta_zzzz_xxxyzz_0 = pbuffer.data(idx_npot_0_gi + 400);

    auto ta_zzzz_xxxzzz_0 = pbuffer.data(idx_npot_0_gi + 401);

    auto ta_zzzz_xxyyyy_0 = pbuffer.data(idx_npot_0_gi + 402);

    auto ta_zzzz_xxyyyz_0 = pbuffer.data(idx_npot_0_gi + 403);

    auto ta_zzzz_xxyyzz_0 = pbuffer.data(idx_npot_0_gi + 404);

    auto ta_zzzz_xxyzzz_0 = pbuffer.data(idx_npot_0_gi + 405);

    auto ta_zzzz_xxzzzz_0 = pbuffer.data(idx_npot_0_gi + 406);

    auto ta_zzzz_xyyyyy_0 = pbuffer.data(idx_npot_0_gi + 407);

    auto ta_zzzz_xyyyyz_0 = pbuffer.data(idx_npot_0_gi + 408);

    auto ta_zzzz_xyyyzz_0 = pbuffer.data(idx_npot_0_gi + 409);

    auto ta_zzzz_xyyzzz_0 = pbuffer.data(idx_npot_0_gi + 410);

    auto ta_zzzz_xyzzzz_0 = pbuffer.data(idx_npot_0_gi + 411);

    auto ta_zzzz_xzzzzz_0 = pbuffer.data(idx_npot_0_gi + 412);

    auto ta_zzzz_yyyyyy_0 = pbuffer.data(idx_npot_0_gi + 413);

    auto ta_zzzz_yyyyyz_0 = pbuffer.data(idx_npot_0_gi + 414);

    auto ta_zzzz_yyyyzz_0 = pbuffer.data(idx_npot_0_gi + 415);

    auto ta_zzzz_yyyzzz_0 = pbuffer.data(idx_npot_0_gi + 416);

    auto ta_zzzz_yyzzzz_0 = pbuffer.data(idx_npot_0_gi + 417);

    auto ta_zzzz_yzzzzz_0 = pbuffer.data(idx_npot_0_gi + 418);

    auto ta_zzzz_zzzzzz_0 = pbuffer.data(idx_npot_0_gi + 419);

    #pragma omp simd aligned(pa_z, pc_z, ta_zz_xxxxxx_0, ta_zz_xxxxxx_1, ta_zz_xxxxxy_0, ta_zz_xxxxxy_1, ta_zz_xxxxxz_0, ta_zz_xxxxxz_1, ta_zz_xxxxyy_0, ta_zz_xxxxyy_1, ta_zz_xxxxyz_0, ta_zz_xxxxyz_1, ta_zz_xxxxzz_0, ta_zz_xxxxzz_1, ta_zz_xxxyyy_0, ta_zz_xxxyyy_1, ta_zz_xxxyyz_0, ta_zz_xxxyyz_1, ta_zz_xxxyzz_0, ta_zz_xxxyzz_1, ta_zz_xxxzzz_0, ta_zz_xxxzzz_1, ta_zz_xxyyyy_0, ta_zz_xxyyyy_1, ta_zz_xxyyyz_0, ta_zz_xxyyyz_1, ta_zz_xxyyzz_0, ta_zz_xxyyzz_1, ta_zz_xxyzzz_0, ta_zz_xxyzzz_1, ta_zz_xxzzzz_0, ta_zz_xxzzzz_1, ta_zz_xyyyyy_0, ta_zz_xyyyyy_1, ta_zz_xyyyyz_0, ta_zz_xyyyyz_1, ta_zz_xyyyzz_0, ta_zz_xyyyzz_1, ta_zz_xyyzzz_0, ta_zz_xyyzzz_1, ta_zz_xyzzzz_0, ta_zz_xyzzzz_1, ta_zz_xzzzzz_0, ta_zz_xzzzzz_1, ta_zz_yyyyyy_0, ta_zz_yyyyyy_1, ta_zz_yyyyyz_0, ta_zz_yyyyyz_1, ta_zz_yyyyzz_0, ta_zz_yyyyzz_1, ta_zz_yyyzzz_0, ta_zz_yyyzzz_1, ta_zz_yyzzzz_0, ta_zz_yyzzzz_1, ta_zz_yzzzzz_0, ta_zz_yzzzzz_1, ta_zz_zzzzzz_0, ta_zz_zzzzzz_1, ta_zzz_xxxxx_0, ta_zzz_xxxxx_1, ta_zzz_xxxxxx_0, ta_zzz_xxxxxx_1, ta_zzz_xxxxxy_0, ta_zzz_xxxxxy_1, ta_zzz_xxxxxz_0, ta_zzz_xxxxxz_1, ta_zzz_xxxxy_0, ta_zzz_xxxxy_1, ta_zzz_xxxxyy_0, ta_zzz_xxxxyy_1, ta_zzz_xxxxyz_0, ta_zzz_xxxxyz_1, ta_zzz_xxxxz_0, ta_zzz_xxxxz_1, ta_zzz_xxxxzz_0, ta_zzz_xxxxzz_1, ta_zzz_xxxyy_0, ta_zzz_xxxyy_1, ta_zzz_xxxyyy_0, ta_zzz_xxxyyy_1, ta_zzz_xxxyyz_0, ta_zzz_xxxyyz_1, ta_zzz_xxxyz_0, ta_zzz_xxxyz_1, ta_zzz_xxxyzz_0, ta_zzz_xxxyzz_1, ta_zzz_xxxzz_0, ta_zzz_xxxzz_1, ta_zzz_xxxzzz_0, ta_zzz_xxxzzz_1, ta_zzz_xxyyy_0, ta_zzz_xxyyy_1, ta_zzz_xxyyyy_0, ta_zzz_xxyyyy_1, ta_zzz_xxyyyz_0, ta_zzz_xxyyyz_1, ta_zzz_xxyyz_0, ta_zzz_xxyyz_1, ta_zzz_xxyyzz_0, ta_zzz_xxyyzz_1, ta_zzz_xxyzz_0, ta_zzz_xxyzz_1, ta_zzz_xxyzzz_0, ta_zzz_xxyzzz_1, ta_zzz_xxzzz_0, ta_zzz_xxzzz_1, ta_zzz_xxzzzz_0, ta_zzz_xxzzzz_1, ta_zzz_xyyyy_0, ta_zzz_xyyyy_1, ta_zzz_xyyyyy_0, ta_zzz_xyyyyy_1, ta_zzz_xyyyyz_0, ta_zzz_xyyyyz_1, ta_zzz_xyyyz_0, ta_zzz_xyyyz_1, ta_zzz_xyyyzz_0, ta_zzz_xyyyzz_1, ta_zzz_xyyzz_0, ta_zzz_xyyzz_1, ta_zzz_xyyzzz_0, ta_zzz_xyyzzz_1, ta_zzz_xyzzz_0, ta_zzz_xyzzz_1, ta_zzz_xyzzzz_0, ta_zzz_xyzzzz_1, ta_zzz_xzzzz_0, ta_zzz_xzzzz_1, ta_zzz_xzzzzz_0, ta_zzz_xzzzzz_1, ta_zzz_yyyyy_0, ta_zzz_yyyyy_1, ta_zzz_yyyyyy_0, ta_zzz_yyyyyy_1, ta_zzz_yyyyyz_0, ta_zzz_yyyyyz_1, ta_zzz_yyyyz_0, ta_zzz_yyyyz_1, ta_zzz_yyyyzz_0, ta_zzz_yyyyzz_1, ta_zzz_yyyzz_0, ta_zzz_yyyzz_1, ta_zzz_yyyzzz_0, ta_zzz_yyyzzz_1, ta_zzz_yyzzz_0, ta_zzz_yyzzz_1, ta_zzz_yyzzzz_0, ta_zzz_yyzzzz_1, ta_zzz_yzzzz_0, ta_zzz_yzzzz_1, ta_zzz_yzzzzz_0, ta_zzz_yzzzzz_1, ta_zzz_zzzzz_0, ta_zzz_zzzzz_1, ta_zzz_zzzzzz_0, ta_zzz_zzzzzz_1, ta_zzzz_xxxxxx_0, ta_zzzz_xxxxxy_0, ta_zzzz_xxxxxz_0, ta_zzzz_xxxxyy_0, ta_zzzz_xxxxyz_0, ta_zzzz_xxxxzz_0, ta_zzzz_xxxyyy_0, ta_zzzz_xxxyyz_0, ta_zzzz_xxxyzz_0, ta_zzzz_xxxzzz_0, ta_zzzz_xxyyyy_0, ta_zzzz_xxyyyz_0, ta_zzzz_xxyyzz_0, ta_zzzz_xxyzzz_0, ta_zzzz_xxzzzz_0, ta_zzzz_xyyyyy_0, ta_zzzz_xyyyyz_0, ta_zzzz_xyyyzz_0, ta_zzzz_xyyzzz_0, ta_zzzz_xyzzzz_0, ta_zzzz_xzzzzz_0, ta_zzzz_yyyyyy_0, ta_zzzz_yyyyyz_0, ta_zzzz_yyyyzz_0, ta_zzzz_yyyzzz_0, ta_zzzz_yyzzzz_0, ta_zzzz_yzzzzz_0, ta_zzzz_zzzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzz_xxxxxx_0[i] = 3.0 * ta_zz_xxxxxx_0[i] * fe_0 - 3.0 * ta_zz_xxxxxx_1[i] * fe_0 + ta_zzz_xxxxxx_0[i] * pa_z[i] - ta_zzz_xxxxxx_1[i] * pc_z[i];

        ta_zzzz_xxxxxy_0[i] = 3.0 * ta_zz_xxxxxy_0[i] * fe_0 - 3.0 * ta_zz_xxxxxy_1[i] * fe_0 + ta_zzz_xxxxxy_0[i] * pa_z[i] - ta_zzz_xxxxxy_1[i] * pc_z[i];

        ta_zzzz_xxxxxz_0[i] = 3.0 * ta_zz_xxxxxz_0[i] * fe_0 - 3.0 * ta_zz_xxxxxz_1[i] * fe_0 + ta_zzz_xxxxx_0[i] * fe_0 - ta_zzz_xxxxx_1[i] * fe_0 + ta_zzz_xxxxxz_0[i] * pa_z[i] - ta_zzz_xxxxxz_1[i] * pc_z[i];

        ta_zzzz_xxxxyy_0[i] = 3.0 * ta_zz_xxxxyy_0[i] * fe_0 - 3.0 * ta_zz_xxxxyy_1[i] * fe_0 + ta_zzz_xxxxyy_0[i] * pa_z[i] - ta_zzz_xxxxyy_1[i] * pc_z[i];

        ta_zzzz_xxxxyz_0[i] = 3.0 * ta_zz_xxxxyz_0[i] * fe_0 - 3.0 * ta_zz_xxxxyz_1[i] * fe_0 + ta_zzz_xxxxy_0[i] * fe_0 - ta_zzz_xxxxy_1[i] * fe_0 + ta_zzz_xxxxyz_0[i] * pa_z[i] - ta_zzz_xxxxyz_1[i] * pc_z[i];

        ta_zzzz_xxxxzz_0[i] = 3.0 * ta_zz_xxxxzz_0[i] * fe_0 - 3.0 * ta_zz_xxxxzz_1[i] * fe_0 + 2.0 * ta_zzz_xxxxz_0[i] * fe_0 - 2.0 * ta_zzz_xxxxz_1[i] * fe_0 + ta_zzz_xxxxzz_0[i] * pa_z[i] - ta_zzz_xxxxzz_1[i] * pc_z[i];

        ta_zzzz_xxxyyy_0[i] = 3.0 * ta_zz_xxxyyy_0[i] * fe_0 - 3.0 * ta_zz_xxxyyy_1[i] * fe_0 + ta_zzz_xxxyyy_0[i] * pa_z[i] - ta_zzz_xxxyyy_1[i] * pc_z[i];

        ta_zzzz_xxxyyz_0[i] = 3.0 * ta_zz_xxxyyz_0[i] * fe_0 - 3.0 * ta_zz_xxxyyz_1[i] * fe_0 + ta_zzz_xxxyy_0[i] * fe_0 - ta_zzz_xxxyy_1[i] * fe_0 + ta_zzz_xxxyyz_0[i] * pa_z[i] - ta_zzz_xxxyyz_1[i] * pc_z[i];

        ta_zzzz_xxxyzz_0[i] = 3.0 * ta_zz_xxxyzz_0[i] * fe_0 - 3.0 * ta_zz_xxxyzz_1[i] * fe_0 + 2.0 * ta_zzz_xxxyz_0[i] * fe_0 - 2.0 * ta_zzz_xxxyz_1[i] * fe_0 + ta_zzz_xxxyzz_0[i] * pa_z[i] - ta_zzz_xxxyzz_1[i] * pc_z[i];

        ta_zzzz_xxxzzz_0[i] = 3.0 * ta_zz_xxxzzz_0[i] * fe_0 - 3.0 * ta_zz_xxxzzz_1[i] * fe_0 + 3.0 * ta_zzz_xxxzz_0[i] * fe_0 - 3.0 * ta_zzz_xxxzz_1[i] * fe_0 + ta_zzz_xxxzzz_0[i] * pa_z[i] - ta_zzz_xxxzzz_1[i] * pc_z[i];

        ta_zzzz_xxyyyy_0[i] = 3.0 * ta_zz_xxyyyy_0[i] * fe_0 - 3.0 * ta_zz_xxyyyy_1[i] * fe_0 + ta_zzz_xxyyyy_0[i] * pa_z[i] - ta_zzz_xxyyyy_1[i] * pc_z[i];

        ta_zzzz_xxyyyz_0[i] = 3.0 * ta_zz_xxyyyz_0[i] * fe_0 - 3.0 * ta_zz_xxyyyz_1[i] * fe_0 + ta_zzz_xxyyy_0[i] * fe_0 - ta_zzz_xxyyy_1[i] * fe_0 + ta_zzz_xxyyyz_0[i] * pa_z[i] - ta_zzz_xxyyyz_1[i] * pc_z[i];

        ta_zzzz_xxyyzz_0[i] = 3.0 * ta_zz_xxyyzz_0[i] * fe_0 - 3.0 * ta_zz_xxyyzz_1[i] * fe_0 + 2.0 * ta_zzz_xxyyz_0[i] * fe_0 - 2.0 * ta_zzz_xxyyz_1[i] * fe_0 + ta_zzz_xxyyzz_0[i] * pa_z[i] - ta_zzz_xxyyzz_1[i] * pc_z[i];

        ta_zzzz_xxyzzz_0[i] = 3.0 * ta_zz_xxyzzz_0[i] * fe_0 - 3.0 * ta_zz_xxyzzz_1[i] * fe_0 + 3.0 * ta_zzz_xxyzz_0[i] * fe_0 - 3.0 * ta_zzz_xxyzz_1[i] * fe_0 + ta_zzz_xxyzzz_0[i] * pa_z[i] - ta_zzz_xxyzzz_1[i] * pc_z[i];

        ta_zzzz_xxzzzz_0[i] = 3.0 * ta_zz_xxzzzz_0[i] * fe_0 - 3.0 * ta_zz_xxzzzz_1[i] * fe_0 + 4.0 * ta_zzz_xxzzz_0[i] * fe_0 - 4.0 * ta_zzz_xxzzz_1[i] * fe_0 + ta_zzz_xxzzzz_0[i] * pa_z[i] - ta_zzz_xxzzzz_1[i] * pc_z[i];

        ta_zzzz_xyyyyy_0[i] = 3.0 * ta_zz_xyyyyy_0[i] * fe_0 - 3.0 * ta_zz_xyyyyy_1[i] * fe_0 + ta_zzz_xyyyyy_0[i] * pa_z[i] - ta_zzz_xyyyyy_1[i] * pc_z[i];

        ta_zzzz_xyyyyz_0[i] = 3.0 * ta_zz_xyyyyz_0[i] * fe_0 - 3.0 * ta_zz_xyyyyz_1[i] * fe_0 + ta_zzz_xyyyy_0[i] * fe_0 - ta_zzz_xyyyy_1[i] * fe_0 + ta_zzz_xyyyyz_0[i] * pa_z[i] - ta_zzz_xyyyyz_1[i] * pc_z[i];

        ta_zzzz_xyyyzz_0[i] = 3.0 * ta_zz_xyyyzz_0[i] * fe_0 - 3.0 * ta_zz_xyyyzz_1[i] * fe_0 + 2.0 * ta_zzz_xyyyz_0[i] * fe_0 - 2.0 * ta_zzz_xyyyz_1[i] * fe_0 + ta_zzz_xyyyzz_0[i] * pa_z[i] - ta_zzz_xyyyzz_1[i] * pc_z[i];

        ta_zzzz_xyyzzz_0[i] = 3.0 * ta_zz_xyyzzz_0[i] * fe_0 - 3.0 * ta_zz_xyyzzz_1[i] * fe_0 + 3.0 * ta_zzz_xyyzz_0[i] * fe_0 - 3.0 * ta_zzz_xyyzz_1[i] * fe_0 + ta_zzz_xyyzzz_0[i] * pa_z[i] - ta_zzz_xyyzzz_1[i] * pc_z[i];

        ta_zzzz_xyzzzz_0[i] = 3.0 * ta_zz_xyzzzz_0[i] * fe_0 - 3.0 * ta_zz_xyzzzz_1[i] * fe_0 + 4.0 * ta_zzz_xyzzz_0[i] * fe_0 - 4.0 * ta_zzz_xyzzz_1[i] * fe_0 + ta_zzz_xyzzzz_0[i] * pa_z[i] - ta_zzz_xyzzzz_1[i] * pc_z[i];

        ta_zzzz_xzzzzz_0[i] = 3.0 * ta_zz_xzzzzz_0[i] * fe_0 - 3.0 * ta_zz_xzzzzz_1[i] * fe_0 + 5.0 * ta_zzz_xzzzz_0[i] * fe_0 - 5.0 * ta_zzz_xzzzz_1[i] * fe_0 + ta_zzz_xzzzzz_0[i] * pa_z[i] - ta_zzz_xzzzzz_1[i] * pc_z[i];

        ta_zzzz_yyyyyy_0[i] = 3.0 * ta_zz_yyyyyy_0[i] * fe_0 - 3.0 * ta_zz_yyyyyy_1[i] * fe_0 + ta_zzz_yyyyyy_0[i] * pa_z[i] - ta_zzz_yyyyyy_1[i] * pc_z[i];

        ta_zzzz_yyyyyz_0[i] = 3.0 * ta_zz_yyyyyz_0[i] * fe_0 - 3.0 * ta_zz_yyyyyz_1[i] * fe_0 + ta_zzz_yyyyy_0[i] * fe_0 - ta_zzz_yyyyy_1[i] * fe_0 + ta_zzz_yyyyyz_0[i] * pa_z[i] - ta_zzz_yyyyyz_1[i] * pc_z[i];

        ta_zzzz_yyyyzz_0[i] = 3.0 * ta_zz_yyyyzz_0[i] * fe_0 - 3.0 * ta_zz_yyyyzz_1[i] * fe_0 + 2.0 * ta_zzz_yyyyz_0[i] * fe_0 - 2.0 * ta_zzz_yyyyz_1[i] * fe_0 + ta_zzz_yyyyzz_0[i] * pa_z[i] - ta_zzz_yyyyzz_1[i] * pc_z[i];

        ta_zzzz_yyyzzz_0[i] = 3.0 * ta_zz_yyyzzz_0[i] * fe_0 - 3.0 * ta_zz_yyyzzz_1[i] * fe_0 + 3.0 * ta_zzz_yyyzz_0[i] * fe_0 - 3.0 * ta_zzz_yyyzz_1[i] * fe_0 + ta_zzz_yyyzzz_0[i] * pa_z[i] - ta_zzz_yyyzzz_1[i] * pc_z[i];

        ta_zzzz_yyzzzz_0[i] = 3.0 * ta_zz_yyzzzz_0[i] * fe_0 - 3.0 * ta_zz_yyzzzz_1[i] * fe_0 + 4.0 * ta_zzz_yyzzz_0[i] * fe_0 - 4.0 * ta_zzz_yyzzz_1[i] * fe_0 + ta_zzz_yyzzzz_0[i] * pa_z[i] - ta_zzz_yyzzzz_1[i] * pc_z[i];

        ta_zzzz_yzzzzz_0[i] = 3.0 * ta_zz_yzzzzz_0[i] * fe_0 - 3.0 * ta_zz_yzzzzz_1[i] * fe_0 + 5.0 * ta_zzz_yzzzz_0[i] * fe_0 - 5.0 * ta_zzz_yzzzz_1[i] * fe_0 + ta_zzz_yzzzzz_0[i] * pa_z[i] - ta_zzz_yzzzzz_1[i] * pc_z[i];

        ta_zzzz_zzzzzz_0[i] = 3.0 * ta_zz_zzzzzz_0[i] * fe_0 - 3.0 * ta_zz_zzzzzz_1[i] * fe_0 + 6.0 * ta_zzz_zzzzz_0[i] * fe_0 - 6.0 * ta_zzz_zzzzz_1[i] * fe_0 + ta_zzz_zzzzzz_0[i] * pa_z[i] - ta_zzz_zzzzzz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

