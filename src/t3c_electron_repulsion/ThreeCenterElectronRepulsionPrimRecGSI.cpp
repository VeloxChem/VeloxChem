#include "ThreeCenterElectronRepulsionPrimRecGSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsi,
                                 size_t idx_eri_0_dsi,
                                 size_t idx_eri_1_dsi,
                                 size_t idx_eri_1_fsh,
                                 size_t idx_eri_1_fsi,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : DSI

    auto g_xx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi);

    auto g_xx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 1);

    auto g_xx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 2);

    auto g_xx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 3);

    auto g_xx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 4);

    auto g_xx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 5);

    auto g_xx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 6);

    auto g_xx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 7);

    auto g_xx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 8);

    auto g_xx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 9);

    auto g_xx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 10);

    auto g_xx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 11);

    auto g_xx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 12);

    auto g_xx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 13);

    auto g_xx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 14);

    auto g_xx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 15);

    auto g_xx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 16);

    auto g_xx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 17);

    auto g_xx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 18);

    auto g_xx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 19);

    auto g_xx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 20);

    auto g_xx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 21);

    auto g_xx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 22);

    auto g_xx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 23);

    auto g_xx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 24);

    auto g_xx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 25);

    auto g_xx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 26);

    auto g_xx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 27);

    auto g_yy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 84);

    auto g_yy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 85);

    auto g_yy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 86);

    auto g_yy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 87);

    auto g_yy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 88);

    auto g_yy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 89);

    auto g_yy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 90);

    auto g_yy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 91);

    auto g_yy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 92);

    auto g_yy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 93);

    auto g_yy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 94);

    auto g_yy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 95);

    auto g_yy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 96);

    auto g_yy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 97);

    auto g_yy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 98);

    auto g_yy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 99);

    auto g_yy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 100);

    auto g_yy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 101);

    auto g_yy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 102);

    auto g_yy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 103);

    auto g_yy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 104);

    auto g_yy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 105);

    auto g_yy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 106);

    auto g_yy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 107);

    auto g_yy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 108);

    auto g_yy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 109);

    auto g_yy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 110);

    auto g_yy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 111);

    auto g_zz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_dsi + 140);

    auto g_zz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_dsi + 141);

    auto g_zz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_dsi + 142);

    auto g_zz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_dsi + 143);

    auto g_zz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_dsi + 144);

    auto g_zz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_dsi + 145);

    auto g_zz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_dsi + 146);

    auto g_zz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_dsi + 147);

    auto g_zz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_dsi + 148);

    auto g_zz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_dsi + 149);

    auto g_zz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_dsi + 150);

    auto g_zz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_dsi + 151);

    auto g_zz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_dsi + 152);

    auto g_zz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_dsi + 153);

    auto g_zz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_dsi + 154);

    auto g_zz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 155);

    auto g_zz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 156);

    auto g_zz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 157);

    auto g_zz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 158);

    auto g_zz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 159);

    auto g_zz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 160);

    auto g_zz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_dsi + 161);

    auto g_zz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_dsi + 162);

    auto g_zz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_dsi + 163);

    auto g_zz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_dsi + 164);

    auto g_zz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_dsi + 165);

    auto g_zz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 166);

    auto g_zz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_dsi + 167);

    /// Set up components of auxilary buffer : DSI

    auto g_xx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi);

    auto g_xx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 1);

    auto g_xx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 2);

    auto g_xx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 3);

    auto g_xx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 4);

    auto g_xx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 5);

    auto g_xx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 6);

    auto g_xx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 7);

    auto g_xx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 8);

    auto g_xx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 9);

    auto g_xx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 10);

    auto g_xx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 11);

    auto g_xx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 12);

    auto g_xx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 13);

    auto g_xx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 14);

    auto g_xx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 15);

    auto g_xx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 16);

    auto g_xx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 17);

    auto g_xx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 18);

    auto g_xx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 19);

    auto g_xx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 20);

    auto g_xx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 21);

    auto g_xx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 22);

    auto g_xx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 23);

    auto g_xx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 24);

    auto g_xx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 25);

    auto g_xx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 26);

    auto g_xx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 27);

    auto g_yy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi + 84);

    auto g_yy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 85);

    auto g_yy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 86);

    auto g_yy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 87);

    auto g_yy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 88);

    auto g_yy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 89);

    auto g_yy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 90);

    auto g_yy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 91);

    auto g_yy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 92);

    auto g_yy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 93);

    auto g_yy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 94);

    auto g_yy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 95);

    auto g_yy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 96);

    auto g_yy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 97);

    auto g_yy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 98);

    auto g_yy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 99);

    auto g_yy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 100);

    auto g_yy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 101);

    auto g_yy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 102);

    auto g_yy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 103);

    auto g_yy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 104);

    auto g_yy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 105);

    auto g_yy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 106);

    auto g_yy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 107);

    auto g_yy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 108);

    auto g_yy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 109);

    auto g_yy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 110);

    auto g_yy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 111);

    auto g_zz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_dsi + 140);

    auto g_zz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_dsi + 141);

    auto g_zz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_dsi + 142);

    auto g_zz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_dsi + 143);

    auto g_zz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_dsi + 144);

    auto g_zz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_dsi + 145);

    auto g_zz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_dsi + 146);

    auto g_zz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_dsi + 147);

    auto g_zz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_dsi + 148);

    auto g_zz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_dsi + 149);

    auto g_zz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_dsi + 150);

    auto g_zz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_dsi + 151);

    auto g_zz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_dsi + 152);

    auto g_zz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_dsi + 153);

    auto g_zz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_dsi + 154);

    auto g_zz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 155);

    auto g_zz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 156);

    auto g_zz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 157);

    auto g_zz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 158);

    auto g_zz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 159);

    auto g_zz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 160);

    auto g_zz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_dsi + 161);

    auto g_zz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_dsi + 162);

    auto g_zz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_dsi + 163);

    auto g_zz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_dsi + 164);

    auto g_zz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_dsi + 165);

    auto g_zz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 166);

    auto g_zz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_dsi + 167);

    /// Set up components of auxilary buffer : FSH

    auto g_xxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh);

    auto g_xxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 1);

    auto g_xxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 2);

    auto g_xxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 3);

    auto g_xxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 4);

    auto g_xxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 5);

    auto g_xxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 6);

    auto g_xxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 7);

    auto g_xxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 8);

    auto g_xxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 9);

    auto g_xxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 10);

    auto g_xxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 11);

    auto g_xxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 12);

    auto g_xxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 13);

    auto g_xxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 14);

    auto g_xxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 15);

    auto g_xxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 16);

    auto g_xxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 17);

    auto g_xxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 18);

    auto g_xxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 19);

    auto g_xxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 20);

    auto g_xxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 44);

    auto g_xxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 46);

    auto g_xxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 47);

    auto g_xxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 49);

    auto g_xxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 50);

    auto g_xxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 51);

    auto g_xxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 53);

    auto g_xxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 54);

    auto g_xxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 55);

    auto g_xxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 56);

    auto g_xxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 58);

    auto g_xxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 59);

    auto g_xxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 60);

    auto g_xxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 61);

    auto g_xxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 62);

    auto g_xyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 64);

    auto g_xyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 66);

    auto g_xyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 67);

    auto g_xyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 69);

    auto g_xyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 70);

    auto g_xyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 71);

    auto g_xyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 73);

    auto g_xyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 74);

    auto g_xyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 75);

    auto g_xyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 76);

    auto g_xyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 78);

    auto g_xyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 79);

    auto g_xyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 80);

    auto g_xyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 81);

    auto g_xyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 82);

    auto g_xzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 107);

    auto g_xzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 109);

    auto g_xzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 110);

    auto g_xzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 112);

    auto g_xzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 113);

    auto g_xzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 114);

    auto g_xzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 116);

    auto g_xzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 117);

    auto g_xzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 118);

    auto g_xzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 119);

    auto g_xzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 121);

    auto g_xzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 122);

    auto g_xzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 123);

    auto g_xzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 124);

    auto g_xzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 125);

    auto g_yyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 126);

    auto g_yyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 127);

    auto g_yyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 128);

    auto g_yyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 129);

    auto g_yyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 130);

    auto g_yyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 131);

    auto g_yyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 132);

    auto g_yyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 133);

    auto g_yyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 134);

    auto g_yyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 135);

    auto g_yyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 136);

    auto g_yyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 137);

    auto g_yyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 138);

    auto g_yyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 139);

    auto g_yyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 140);

    auto g_yyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 141);

    auto g_yyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 142);

    auto g_yyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 143);

    auto g_yyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 144);

    auto g_yyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 145);

    auto g_yyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 146);

    auto g_yyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 149);

    auto g_yyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 151);

    auto g_yyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 152);

    auto g_yyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 154);

    auto g_yyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 155);

    auto g_yyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 156);

    auto g_yyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 158);

    auto g_yyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 159);

    auto g_yyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 160);

    auto g_yyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 161);

    auto g_yyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 163);

    auto g_yyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 164);

    auto g_yyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 165);

    auto g_yyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 166);

    auto g_yyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 167);

    auto g_yzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 169);

    auto g_yzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 170);

    auto g_yzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 171);

    auto g_yzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 172);

    auto g_yzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 173);

    auto g_yzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 174);

    auto g_yzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 175);

    auto g_yzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 176);

    auto g_yzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 177);

    auto g_yzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 178);

    auto g_yzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 179);

    auto g_yzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 180);

    auto g_yzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 181);

    auto g_yzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 182);

    auto g_yzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 183);

    auto g_yzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 184);

    auto g_yzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 185);

    auto g_yzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 186);

    auto g_yzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 187);

    auto g_yzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 188);

    auto g_zzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_fsh + 189);

    auto g_zzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_fsh + 190);

    auto g_zzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_fsh + 191);

    auto g_zzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_fsh + 192);

    auto g_zzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_fsh + 193);

    auto g_zzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_fsh + 194);

    auto g_zzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_fsh + 195);

    auto g_zzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_fsh + 196);

    auto g_zzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_fsh + 197);

    auto g_zzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_fsh + 198);

    auto g_zzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_fsh + 199);

    auto g_zzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_fsh + 200);

    auto g_zzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_fsh + 201);

    auto g_zzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_fsh + 202);

    auto g_zzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_fsh + 203);

    auto g_zzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_fsh + 204);

    auto g_zzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_fsh + 205);

    auto g_zzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_fsh + 206);

    auto g_zzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_fsh + 207);

    auto g_zzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_fsh + 208);

    auto g_zzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_fsh + 209);

    /// Set up components of auxilary buffer : FSI

    auto g_xxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi);

    auto g_xxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 1);

    auto g_xxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 2);

    auto g_xxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 3);

    auto g_xxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 4);

    auto g_xxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 5);

    auto g_xxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 6);

    auto g_xxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 7);

    auto g_xxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 8);

    auto g_xxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 9);

    auto g_xxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 10);

    auto g_xxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 11);

    auto g_xxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 12);

    auto g_xxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 13);

    auto g_xxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 14);

    auto g_xxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 15);

    auto g_xxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 16);

    auto g_xxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 17);

    auto g_xxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 18);

    auto g_xxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 19);

    auto g_xxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 20);

    auto g_xxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 21);

    auto g_xxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 22);

    auto g_xxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 23);

    auto g_xxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 24);

    auto g_xxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 25);

    auto g_xxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 26);

    auto g_xxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 27);

    auto g_xxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 28);

    auto g_xxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 29);

    auto g_xxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 30);

    auto g_xxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 31);

    auto g_xxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 33);

    auto g_xxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 34);

    auto g_xxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 37);

    auto g_xxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 38);

    auto g_xxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 42);

    auto g_xxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 43);

    auto g_xxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 48);

    auto g_xxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 49);

    auto g_xxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 56);

    auto g_xxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 57);

    auto g_xxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 58);

    auto g_xxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 59);

    auto g_xxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 60);

    auto g_xxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 61);

    auto g_xxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 62);

    auto g_xxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 63);

    auto g_xxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 64);

    auto g_xxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 65);

    auto g_xxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 66);

    auto g_xxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 67);

    auto g_xxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 68);

    auto g_xxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 69);

    auto g_xxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 70);

    auto g_xxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 71);

    auto g_xxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 72);

    auto g_xxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 73);

    auto g_xxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 74);

    auto g_xxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 75);

    auto g_xxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 76);

    auto g_xxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 78);

    auto g_xxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 79);

    auto g_xxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 80);

    auto g_xxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 81);

    auto g_xxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 82);

    auto g_xxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 83);

    auto g_xyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 84);

    auto g_xyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 85);

    auto g_xyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 87);

    auto g_xyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 88);

    auto g_xyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 90);

    auto g_xyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 91);

    auto g_xyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 92);

    auto g_xyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 94);

    auto g_xyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 95);

    auto g_xyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 96);

    auto g_xyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 97);

    auto g_xyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 99);

    auto g_xyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 100);

    auto g_xyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 101);

    auto g_xyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 102);

    auto g_xyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 103);

    auto g_xyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 105);

    auto g_xyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 106);

    auto g_xyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 107);

    auto g_xyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 108);

    auto g_xyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 109);

    auto g_xyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 110);

    auto g_xyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 111);

    auto g_xzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 140);

    auto g_xzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 142);

    auto g_xzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 144);

    auto g_xzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 145);

    auto g_xzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 147);

    auto g_xzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 148);

    auto g_xzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 149);

    auto g_xzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 151);

    auto g_xzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 152);

    auto g_xzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 153);

    auto g_xzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 154);

    auto g_xzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 156);

    auto g_xzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 157);

    auto g_xzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 158);

    auto g_xzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 159);

    auto g_xzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 160);

    auto g_xzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 161);

    auto g_xzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 162);

    auto g_xzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 163);

    auto g_xzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 164);

    auto g_xzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 165);

    auto g_xzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 166);

    auto g_xzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 167);

    auto g_yyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 168);

    auto g_yyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 169);

    auto g_yyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 170);

    auto g_yyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 171);

    auto g_yyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 172);

    auto g_yyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 173);

    auto g_yyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 174);

    auto g_yyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 175);

    auto g_yyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 176);

    auto g_yyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 177);

    auto g_yyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 178);

    auto g_yyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 179);

    auto g_yyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 180);

    auto g_yyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 181);

    auto g_yyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 182);

    auto g_yyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 183);

    auto g_yyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 184);

    auto g_yyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 185);

    auto g_yyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 186);

    auto g_yyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 187);

    auto g_yyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 188);

    auto g_yyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 189);

    auto g_yyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 190);

    auto g_yyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 191);

    auto g_yyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 192);

    auto g_yyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 193);

    auto g_yyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 194);

    auto g_yyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 195);

    auto g_yyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 197);

    auto g_yyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 198);

    auto g_yyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 199);

    auto g_yyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 200);

    auto g_yyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 201);

    auto g_yyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 202);

    auto g_yyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 203);

    auto g_yyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 204);

    auto g_yyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 205);

    auto g_yyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 206);

    auto g_yyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 207);

    auto g_yyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 208);

    auto g_yyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 209);

    auto g_yyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 210);

    auto g_yyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 211);

    auto g_yyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 212);

    auto g_yyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 213);

    auto g_yyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 214);

    auto g_yyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 215);

    auto g_yyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 216);

    auto g_yyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 217);

    auto g_yyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 218);

    auto g_yyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 219);

    auto g_yyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 220);

    auto g_yyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 221);

    auto g_yyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 222);

    auto g_yyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 223);

    auto g_yzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 224);

    auto g_yzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 225);

    auto g_yzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 226);

    auto g_yzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 227);

    auto g_yzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 228);

    auto g_yzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 229);

    auto g_yzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 230);

    auto g_yzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 231);

    auto g_yzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 232);

    auto g_yzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 233);

    auto g_yzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 234);

    auto g_yzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 235);

    auto g_yzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 236);

    auto g_yzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 237);

    auto g_yzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 238);

    auto g_yzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 239);

    auto g_yzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 240);

    auto g_yzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 241);

    auto g_yzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 242);

    auto g_yzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 243);

    auto g_yzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 244);

    auto g_yzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 245);

    auto g_yzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 246);

    auto g_yzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 247);

    auto g_yzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 248);

    auto g_yzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 249);

    auto g_yzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 250);

    auto g_yzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 251);

    auto g_zzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_fsi + 252);

    auto g_zzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_fsi + 253);

    auto g_zzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_fsi + 254);

    auto g_zzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_fsi + 255);

    auto g_zzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_fsi + 256);

    auto g_zzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_fsi + 257);

    auto g_zzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_fsi + 258);

    auto g_zzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_fsi + 259);

    auto g_zzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_fsi + 260);

    auto g_zzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_fsi + 261);

    auto g_zzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_fsi + 262);

    auto g_zzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_fsi + 263);

    auto g_zzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_fsi + 264);

    auto g_zzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_fsi + 265);

    auto g_zzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_fsi + 266);

    auto g_zzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 267);

    auto g_zzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 268);

    auto g_zzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 269);

    auto g_zzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 270);

    auto g_zzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 271);

    auto g_zzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 272);

    auto g_zzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_fsi + 273);

    auto g_zzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_fsi + 274);

    auto g_zzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_fsi + 275);

    auto g_zzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_fsi + 276);

    auto g_zzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_fsi + 277);

    auto g_zzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 278);

    auto g_zzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_fsi + 279);

    /// Set up 0-28 components of targeted buffer : GSI

    auto g_xxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi);

    auto g_xxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 1);

    auto g_xxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 2);

    auto g_xxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 3);

    auto g_xxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 4);

    auto g_xxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 5);

    auto g_xxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 6);

    auto g_xxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 7);

    auto g_xxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 8);

    auto g_xxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 9);

    auto g_xxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 10);

    auto g_xxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 11);

    auto g_xxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 12);

    auto g_xxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 13);

    auto g_xxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 14);

    auto g_xxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 15);

    auto g_xxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 16);

    auto g_xxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 17);

    auto g_xxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 18);

    auto g_xxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 19);

    auto g_xxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 20);

    auto g_xxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 21);

    auto g_xxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 22);

    auto g_xxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 23);

    auto g_xxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 24);

    auto g_xxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 25);

    auto g_xxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 26);

    auto g_xxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 27);

    #pragma omp simd aligned(g_xx_0_xxxxxx_0, g_xx_0_xxxxxx_1, g_xx_0_xxxxxy_0, g_xx_0_xxxxxy_1, g_xx_0_xxxxxz_0, g_xx_0_xxxxxz_1, g_xx_0_xxxxyy_0, g_xx_0_xxxxyy_1, g_xx_0_xxxxyz_0, g_xx_0_xxxxyz_1, g_xx_0_xxxxzz_0, g_xx_0_xxxxzz_1, g_xx_0_xxxyyy_0, g_xx_0_xxxyyy_1, g_xx_0_xxxyyz_0, g_xx_0_xxxyyz_1, g_xx_0_xxxyzz_0, g_xx_0_xxxyzz_1, g_xx_0_xxxzzz_0, g_xx_0_xxxzzz_1, g_xx_0_xxyyyy_0, g_xx_0_xxyyyy_1, g_xx_0_xxyyyz_0, g_xx_0_xxyyyz_1, g_xx_0_xxyyzz_0, g_xx_0_xxyyzz_1, g_xx_0_xxyzzz_0, g_xx_0_xxyzzz_1, g_xx_0_xxzzzz_0, g_xx_0_xxzzzz_1, g_xx_0_xyyyyy_0, g_xx_0_xyyyyy_1, g_xx_0_xyyyyz_0, g_xx_0_xyyyyz_1, g_xx_0_xyyyzz_0, g_xx_0_xyyyzz_1, g_xx_0_xyyzzz_0, g_xx_0_xyyzzz_1, g_xx_0_xyzzzz_0, g_xx_0_xyzzzz_1, g_xx_0_xzzzzz_0, g_xx_0_xzzzzz_1, g_xx_0_yyyyyy_0, g_xx_0_yyyyyy_1, g_xx_0_yyyyyz_0, g_xx_0_yyyyyz_1, g_xx_0_yyyyzz_0, g_xx_0_yyyyzz_1, g_xx_0_yyyzzz_0, g_xx_0_yyyzzz_1, g_xx_0_yyzzzz_0, g_xx_0_yyzzzz_1, g_xx_0_yzzzzz_0, g_xx_0_yzzzzz_1, g_xx_0_zzzzzz_0, g_xx_0_zzzzzz_1, g_xxx_0_xxxxx_1, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxy_1, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxz_1, g_xxx_0_xxxxzz_1, g_xxx_0_xxxyy_1, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyz_1, g_xxx_0_xxxyzz_1, g_xxx_0_xxxzz_1, g_xxx_0_xxxzzz_1, g_xxx_0_xxyyy_1, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyz_1, g_xxx_0_xxyyzz_1, g_xxx_0_xxyzz_1, g_xxx_0_xxyzzz_1, g_xxx_0_xxzzz_1, g_xxx_0_xxzzzz_1, g_xxx_0_xyyyy_1, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyz_1, g_xxx_0_xyyyzz_1, g_xxx_0_xyyzz_1, g_xxx_0_xyyzzz_1, g_xxx_0_xyzzz_1, g_xxx_0_xyzzzz_1, g_xxx_0_xzzzz_1, g_xxx_0_xzzzzz_1, g_xxx_0_yyyyy_1, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyz_1, g_xxx_0_yyyyzz_1, g_xxx_0_yyyzz_1, g_xxx_0_yyyzzz_1, g_xxx_0_yyzzz_1, g_xxx_0_yyzzzz_1, g_xxx_0_yzzzz_1, g_xxx_0_yzzzzz_1, g_xxx_0_zzzzz_1, g_xxx_0_zzzzzz_1, g_xxxx_0_xxxxxx_0, g_xxxx_0_xxxxxy_0, g_xxxx_0_xxxxxz_0, g_xxxx_0_xxxxyy_0, g_xxxx_0_xxxxyz_0, g_xxxx_0_xxxxzz_0, g_xxxx_0_xxxyyy_0, g_xxxx_0_xxxyyz_0, g_xxxx_0_xxxyzz_0, g_xxxx_0_xxxzzz_0, g_xxxx_0_xxyyyy_0, g_xxxx_0_xxyyyz_0, g_xxxx_0_xxyyzz_0, g_xxxx_0_xxyzzz_0, g_xxxx_0_xxzzzz_0, g_xxxx_0_xyyyyy_0, g_xxxx_0_xyyyyz_0, g_xxxx_0_xyyyzz_0, g_xxxx_0_xyyzzz_0, g_xxxx_0_xyzzzz_0, g_xxxx_0_xzzzzz_0, g_xxxx_0_yyyyyy_0, g_xxxx_0_yyyyyz_0, g_xxxx_0_yyyyzz_0, g_xxxx_0_yyyzzz_0, g_xxxx_0_yyzzzz_0, g_xxxx_0_yzzzzz_0, g_xxxx_0_zzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxxxxx_0[i] = 3.0 * g_xx_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxx_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxx_1[i] * wa_x[i];

        g_xxxx_0_xxxxxy_0[i] = 3.0 * g_xx_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxz_0[i] = 3.0 * g_xx_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyy_0[i] = 3.0 * g_xx_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxyz_0[i] = 3.0 * g_xx_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxzz_0[i] = 3.0 * g_xx_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyy_0[i] = 3.0 * g_xx_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxyyz_0[i] = 3.0 * g_xx_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxyzz_0[i] = 3.0 * g_xx_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxzzz_0[i] = 3.0 * g_xx_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyy_0[i] = 3.0 * g_xx_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxyyyz_0[i] = 3.0 * g_xx_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxyyzz_0[i] = 3.0 * g_xx_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxyzzz_0[i] = 3.0 * g_xx_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxzzzz_0[i] = 3.0 * g_xx_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyy_0[i] = 3.0 * g_xx_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xyyyyz_0[i] = 3.0 * g_xx_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xyyyzz_0[i] = 3.0 * g_xx_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xyyzzz_0[i] = 3.0 * g_xx_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyzzz_1[i] * fz_be_0 + g_xxx_0_yyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xyzzzz_0[i] = 3.0 * g_xx_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xzzzzz_0[i] = 3.0 * g_xx_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyy_0[i] = 3.0 * g_xx_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyy_1[i] * wa_x[i];

        g_xxxx_0_yyyyyz_0[i] = 3.0 * g_xx_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyz_1[i] * wa_x[i];

        g_xxxx_0_yyyyzz_0[i] = 3.0 * g_xx_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyzz_1[i] * wa_x[i];

        g_xxxx_0_yyyzzz_0[i] = 3.0 * g_xx_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzz_1[i] * wa_x[i];

        g_xxxx_0_yyzzzz_0[i] = 3.0 * g_xx_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzz_1[i] * wa_x[i];

        g_xxxx_0_yzzzzz_0[i] = 3.0 * g_xx_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzz_1[i] * wa_x[i];

        g_xxxx_0_zzzzzz_0[i] = 3.0 * g_xx_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 28-56 components of targeted buffer : GSI

    auto g_xxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 28);

    auto g_xxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 29);

    auto g_xxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 30);

    auto g_xxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 31);

    auto g_xxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 32);

    auto g_xxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 33);

    auto g_xxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 34);

    auto g_xxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 35);

    auto g_xxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 36);

    auto g_xxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 37);

    auto g_xxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 38);

    auto g_xxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 39);

    auto g_xxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 40);

    auto g_xxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 41);

    auto g_xxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 42);

    auto g_xxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 43);

    auto g_xxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 44);

    auto g_xxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 45);

    auto g_xxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 46);

    auto g_xxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 47);

    auto g_xxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 48);

    auto g_xxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 49);

    auto g_xxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 50);

    auto g_xxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 51);

    auto g_xxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 52);

    auto g_xxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 53);

    auto g_xxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 54);

    auto g_xxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 55);

    #pragma omp simd aligned(g_xxx_0_xxxxx_1, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxy_1, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxz_1, g_xxx_0_xxxxzz_1, g_xxx_0_xxxyy_1, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyz_1, g_xxx_0_xxxyzz_1, g_xxx_0_xxxzz_1, g_xxx_0_xxxzzz_1, g_xxx_0_xxyyy_1, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyz_1, g_xxx_0_xxyyzz_1, g_xxx_0_xxyzz_1, g_xxx_0_xxyzzz_1, g_xxx_0_xxzzz_1, g_xxx_0_xxzzzz_1, g_xxx_0_xyyyy_1, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyz_1, g_xxx_0_xyyyzz_1, g_xxx_0_xyyzz_1, g_xxx_0_xyyzzz_1, g_xxx_0_xyzzz_1, g_xxx_0_xyzzzz_1, g_xxx_0_xzzzz_1, g_xxx_0_xzzzzz_1, g_xxx_0_yyyyy_1, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyz_1, g_xxx_0_yyyyzz_1, g_xxx_0_yyyzz_1, g_xxx_0_yyyzzz_1, g_xxx_0_yyzzz_1, g_xxx_0_yyzzzz_1, g_xxx_0_yzzzz_1, g_xxx_0_yzzzzz_1, g_xxx_0_zzzzz_1, g_xxx_0_zzzzzz_1, g_xxxy_0_xxxxxx_0, g_xxxy_0_xxxxxy_0, g_xxxy_0_xxxxxz_0, g_xxxy_0_xxxxyy_0, g_xxxy_0_xxxxyz_0, g_xxxy_0_xxxxzz_0, g_xxxy_0_xxxyyy_0, g_xxxy_0_xxxyyz_0, g_xxxy_0_xxxyzz_0, g_xxxy_0_xxxzzz_0, g_xxxy_0_xxyyyy_0, g_xxxy_0_xxyyyz_0, g_xxxy_0_xxyyzz_0, g_xxxy_0_xxyzzz_0, g_xxxy_0_xxzzzz_0, g_xxxy_0_xyyyyy_0, g_xxxy_0_xyyyyz_0, g_xxxy_0_xyyyzz_0, g_xxxy_0_xyyzzz_0, g_xxxy_0_xyzzzz_0, g_xxxy_0_xzzzzz_0, g_xxxy_0_yyyyyy_0, g_xxxy_0_yyyyyz_0, g_xxxy_0_yyyyzz_0, g_xxxy_0_yyyzzz_0, g_xxxy_0_yyzzzz_0, g_xxxy_0_yzzzzz_0, g_xxxy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxxxxx_0[i] = g_xxx_0_xxxxxx_1[i] * wa_y[i];

        g_xxxy_0_xxxxxy_0[i] = g_xxx_0_xxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxz_0[i] = g_xxx_0_xxxxxz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyy_0[i] = 2.0 * g_xxx_0_xxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxyz_0[i] = g_xxx_0_xxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxzz_0[i] = g_xxx_0_xxxxzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyy_0[i] = 3.0 * g_xxx_0_xxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxyyz_0[i] = 2.0 * g_xxx_0_xxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxyzz_0[i] = g_xxx_0_xxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxzzz_0[i] = g_xxx_0_xxxzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyy_0[i] = 4.0 * g_xxx_0_xxyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxyyyz_0[i] = 3.0 * g_xxx_0_xxyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxyyzz_0[i] = 2.0 * g_xxx_0_xxyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxyzzz_0[i] = g_xxx_0_xxzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxzzzz_0[i] = g_xxx_0_xxzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyy_0[i] = 5.0 * g_xxx_0_xyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xyyyyz_0[i] = 4.0 * g_xxx_0_xyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xyyyzz_0[i] = 3.0 * g_xxx_0_xyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xyyzzz_0[i] = 2.0 * g_xxx_0_xyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xyzzzz_0[i] = g_xxx_0_xzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xzzzzz_0[i] = g_xxx_0_xzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyy_0[i] = 6.0 * g_xxx_0_yyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyy_1[i] * wa_y[i];

        g_xxxy_0_yyyyyz_0[i] = 5.0 * g_xxx_0_yyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyz_1[i] * wa_y[i];

        g_xxxy_0_yyyyzz_0[i] = 4.0 * g_xxx_0_yyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzz_1[i] * wa_y[i];

        g_xxxy_0_yyyzzz_0[i] = 3.0 * g_xxx_0_yyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzz_1[i] * wa_y[i];

        g_xxxy_0_yyzzzz_0[i] = 2.0 * g_xxx_0_yzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzz_1[i] * wa_y[i];

        g_xxxy_0_yzzzzz_0[i] = g_xxx_0_zzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzz_1[i] * wa_y[i];

        g_xxxy_0_zzzzzz_0[i] = g_xxx_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 56-84 components of targeted buffer : GSI

    auto g_xxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 56);

    auto g_xxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 57);

    auto g_xxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 58);

    auto g_xxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 59);

    auto g_xxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 60);

    auto g_xxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 61);

    auto g_xxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 62);

    auto g_xxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 63);

    auto g_xxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 64);

    auto g_xxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 65);

    auto g_xxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 66);

    auto g_xxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 67);

    auto g_xxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 68);

    auto g_xxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 69);

    auto g_xxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 70);

    auto g_xxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 71);

    auto g_xxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 72);

    auto g_xxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 73);

    auto g_xxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 74);

    auto g_xxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 75);

    auto g_xxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 76);

    auto g_xxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 77);

    auto g_xxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 78);

    auto g_xxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 79);

    auto g_xxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 80);

    auto g_xxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 81);

    auto g_xxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 82);

    auto g_xxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 83);

    #pragma omp simd aligned(g_xxx_0_xxxxx_1, g_xxx_0_xxxxxx_1, g_xxx_0_xxxxxy_1, g_xxx_0_xxxxxz_1, g_xxx_0_xxxxy_1, g_xxx_0_xxxxyy_1, g_xxx_0_xxxxyz_1, g_xxx_0_xxxxz_1, g_xxx_0_xxxxzz_1, g_xxx_0_xxxyy_1, g_xxx_0_xxxyyy_1, g_xxx_0_xxxyyz_1, g_xxx_0_xxxyz_1, g_xxx_0_xxxyzz_1, g_xxx_0_xxxzz_1, g_xxx_0_xxxzzz_1, g_xxx_0_xxyyy_1, g_xxx_0_xxyyyy_1, g_xxx_0_xxyyyz_1, g_xxx_0_xxyyz_1, g_xxx_0_xxyyzz_1, g_xxx_0_xxyzz_1, g_xxx_0_xxyzzz_1, g_xxx_0_xxzzz_1, g_xxx_0_xxzzzz_1, g_xxx_0_xyyyy_1, g_xxx_0_xyyyyy_1, g_xxx_0_xyyyyz_1, g_xxx_0_xyyyz_1, g_xxx_0_xyyyzz_1, g_xxx_0_xyyzz_1, g_xxx_0_xyyzzz_1, g_xxx_0_xyzzz_1, g_xxx_0_xyzzzz_1, g_xxx_0_xzzzz_1, g_xxx_0_xzzzzz_1, g_xxx_0_yyyyy_1, g_xxx_0_yyyyyy_1, g_xxx_0_yyyyyz_1, g_xxx_0_yyyyz_1, g_xxx_0_yyyyzz_1, g_xxx_0_yyyzz_1, g_xxx_0_yyyzzz_1, g_xxx_0_yyzzz_1, g_xxx_0_yyzzzz_1, g_xxx_0_yzzzz_1, g_xxx_0_yzzzzz_1, g_xxx_0_zzzzz_1, g_xxx_0_zzzzzz_1, g_xxxz_0_xxxxxx_0, g_xxxz_0_xxxxxy_0, g_xxxz_0_xxxxxz_0, g_xxxz_0_xxxxyy_0, g_xxxz_0_xxxxyz_0, g_xxxz_0_xxxxzz_0, g_xxxz_0_xxxyyy_0, g_xxxz_0_xxxyyz_0, g_xxxz_0_xxxyzz_0, g_xxxz_0_xxxzzz_0, g_xxxz_0_xxyyyy_0, g_xxxz_0_xxyyyz_0, g_xxxz_0_xxyyzz_0, g_xxxz_0_xxyzzz_0, g_xxxz_0_xxzzzz_0, g_xxxz_0_xyyyyy_0, g_xxxz_0_xyyyyz_0, g_xxxz_0_xyyyzz_0, g_xxxz_0_xyyzzz_0, g_xxxz_0_xyzzzz_0, g_xxxz_0_xzzzzz_0, g_xxxz_0_yyyyyy_0, g_xxxz_0_yyyyyz_0, g_xxxz_0_yyyyzz_0, g_xxxz_0_yyyzzz_0, g_xxxz_0_yyzzzz_0, g_xxxz_0_yzzzzz_0, g_xxxz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxxxxx_0[i] = g_xxx_0_xxxxxx_1[i] * wa_z[i];

        g_xxxz_0_xxxxxy_0[i] = g_xxx_0_xxxxxy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxz_0[i] = g_xxx_0_xxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyy_0[i] = g_xxx_0_xxxxyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxyz_0[i] = g_xxx_0_xxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxzz_0[i] = 2.0 * g_xxx_0_xxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyy_0[i] = g_xxx_0_xxxyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxyyz_0[i] = g_xxx_0_xxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxyzz_0[i] = 2.0 * g_xxx_0_xxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxzzz_0[i] = 3.0 * g_xxx_0_xxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyy_0[i] = g_xxx_0_xxyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxyyyz_0[i] = g_xxx_0_xxyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxyyzz_0[i] = 2.0 * g_xxx_0_xxyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxyzzz_0[i] = 3.0 * g_xxx_0_xxyzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxzzzz_0[i] = 4.0 * g_xxx_0_xxzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyy_0[i] = g_xxx_0_xyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xyyyyz_0[i] = g_xxx_0_xyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xyyyzz_0[i] = 2.0 * g_xxx_0_xyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xyyzzz_0[i] = 3.0 * g_xxx_0_xyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xyzzzz_0[i] = 4.0 * g_xxx_0_xyzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xzzzzz_0[i] = 5.0 * g_xxx_0_xzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyy_0[i] = g_xxx_0_yyyyyy_1[i] * wa_z[i];

        g_xxxz_0_yyyyyz_0[i] = g_xxx_0_yyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyz_1[i] * wa_z[i];

        g_xxxz_0_yyyyzz_0[i] = 2.0 * g_xxx_0_yyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyzz_1[i] * wa_z[i];

        g_xxxz_0_yyyzzz_0[i] = 3.0 * g_xxx_0_yyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzz_1[i] * wa_z[i];

        g_xxxz_0_yyzzzz_0[i] = 4.0 * g_xxx_0_yyzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzz_1[i] * wa_z[i];

        g_xxxz_0_yzzzzz_0[i] = 5.0 * g_xxx_0_yzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzz_1[i] * wa_z[i];

        g_xxxz_0_zzzzzz_0[i] = 6.0 * g_xxx_0_zzzzz_1[i] * fi_acd_0 + g_xxx_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 84-112 components of targeted buffer : GSI

    auto g_xxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 84);

    auto g_xxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 85);

    auto g_xxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 86);

    auto g_xxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 87);

    auto g_xxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 88);

    auto g_xxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 89);

    auto g_xxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 90);

    auto g_xxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 91);

    auto g_xxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 92);

    auto g_xxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 93);

    auto g_xxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 94);

    auto g_xxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 95);

    auto g_xxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 96);

    auto g_xxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 97);

    auto g_xxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 98);

    auto g_xxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 99);

    auto g_xxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 100);

    auto g_xxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 101);

    auto g_xxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 102);

    auto g_xxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 103);

    auto g_xxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 104);

    auto g_xxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 105);

    auto g_xxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 106);

    auto g_xxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 107);

    auto g_xxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 108);

    auto g_xxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 109);

    auto g_xxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 110);

    auto g_xxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 111);

    #pragma omp simd aligned(g_xx_0_xxxxxx_0, g_xx_0_xxxxxx_1, g_xx_0_xxxxxz_0, g_xx_0_xxxxxz_1, g_xx_0_xxxxzz_0, g_xx_0_xxxxzz_1, g_xx_0_xxxzzz_0, g_xx_0_xxxzzz_1, g_xx_0_xxzzzz_0, g_xx_0_xxzzzz_1, g_xx_0_xzzzzz_0, g_xx_0_xzzzzz_1, g_xxy_0_xxxxxx_1, g_xxy_0_xxxxxz_1, g_xxy_0_xxxxzz_1, g_xxy_0_xxxzzz_1, g_xxy_0_xxzzzz_1, g_xxy_0_xzzzzz_1, g_xxyy_0_xxxxxx_0, g_xxyy_0_xxxxxy_0, g_xxyy_0_xxxxxz_0, g_xxyy_0_xxxxyy_0, g_xxyy_0_xxxxyz_0, g_xxyy_0_xxxxzz_0, g_xxyy_0_xxxyyy_0, g_xxyy_0_xxxyyz_0, g_xxyy_0_xxxyzz_0, g_xxyy_0_xxxzzz_0, g_xxyy_0_xxyyyy_0, g_xxyy_0_xxyyyz_0, g_xxyy_0_xxyyzz_0, g_xxyy_0_xxyzzz_0, g_xxyy_0_xxzzzz_0, g_xxyy_0_xyyyyy_0, g_xxyy_0_xyyyyz_0, g_xxyy_0_xyyyzz_0, g_xxyy_0_xyyzzz_0, g_xxyy_0_xyzzzz_0, g_xxyy_0_xzzzzz_0, g_xxyy_0_yyyyyy_0, g_xxyy_0_yyyyyz_0, g_xxyy_0_yyyyzz_0, g_xxyy_0_yyyzzz_0, g_xxyy_0_yyzzzz_0, g_xxyy_0_yzzzzz_0, g_xxyy_0_zzzzzz_0, g_xyy_0_xxxxxy_1, g_xyy_0_xxxxy_1, g_xyy_0_xxxxyy_1, g_xyy_0_xxxxyz_1, g_xyy_0_xxxyy_1, g_xyy_0_xxxyyy_1, g_xyy_0_xxxyyz_1, g_xyy_0_xxxyz_1, g_xyy_0_xxxyzz_1, g_xyy_0_xxyyy_1, g_xyy_0_xxyyyy_1, g_xyy_0_xxyyyz_1, g_xyy_0_xxyyz_1, g_xyy_0_xxyyzz_1, g_xyy_0_xxyzz_1, g_xyy_0_xxyzzz_1, g_xyy_0_xyyyy_1, g_xyy_0_xyyyyy_1, g_xyy_0_xyyyyz_1, g_xyy_0_xyyyz_1, g_xyy_0_xyyyzz_1, g_xyy_0_xyyzz_1, g_xyy_0_xyyzzz_1, g_xyy_0_xyzzz_1, g_xyy_0_xyzzzz_1, g_xyy_0_yyyyy_1, g_xyy_0_yyyyyy_1, g_xyy_0_yyyyyz_1, g_xyy_0_yyyyz_1, g_xyy_0_yyyyzz_1, g_xyy_0_yyyzz_1, g_xyy_0_yyyzzz_1, g_xyy_0_yyzzz_1, g_xyy_0_yyzzzz_1, g_xyy_0_yzzzz_1, g_xyy_0_yzzzzz_1, g_xyy_0_zzzzzz_1, g_yy_0_xxxxxy_0, g_yy_0_xxxxxy_1, g_yy_0_xxxxyy_0, g_yy_0_xxxxyy_1, g_yy_0_xxxxyz_0, g_yy_0_xxxxyz_1, g_yy_0_xxxyyy_0, g_yy_0_xxxyyy_1, g_yy_0_xxxyyz_0, g_yy_0_xxxyyz_1, g_yy_0_xxxyzz_0, g_yy_0_xxxyzz_1, g_yy_0_xxyyyy_0, g_yy_0_xxyyyy_1, g_yy_0_xxyyyz_0, g_yy_0_xxyyyz_1, g_yy_0_xxyyzz_0, g_yy_0_xxyyzz_1, g_yy_0_xxyzzz_0, g_yy_0_xxyzzz_1, g_yy_0_xyyyyy_0, g_yy_0_xyyyyy_1, g_yy_0_xyyyyz_0, g_yy_0_xyyyyz_1, g_yy_0_xyyyzz_0, g_yy_0_xyyyzz_1, g_yy_0_xyyzzz_0, g_yy_0_xyyzzz_1, g_yy_0_xyzzzz_0, g_yy_0_xyzzzz_1, g_yy_0_yyyyyy_0, g_yy_0_yyyyyy_1, g_yy_0_yyyyyz_0, g_yy_0_yyyyyz_1, g_yy_0_yyyyzz_0, g_yy_0_yyyyzz_1, g_yy_0_yyyzzz_0, g_yy_0_yyyzzz_1, g_yy_0_yyzzzz_0, g_yy_0_yyzzzz_1, g_yy_0_yzzzzz_0, g_yy_0_yzzzzz_1, g_yy_0_zzzzzz_0, g_yy_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxxxxx_0[i] = g_xx_0_xxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxx_1[i] * fz_be_0 + g_xxy_0_xxxxxx_1[i] * wa_y[i];

        g_xxyy_0_xxxxxy_0[i] = g_yy_0_xxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxy_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxy_1[i] * fi_acd_0 + g_xyy_0_xxxxxy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxz_0[i] = g_xx_0_xxxxxz_0[i] * fbe_0 - g_xx_0_xxxxxz_1[i] * fz_be_0 + g_xxy_0_xxxxxz_1[i] * wa_y[i];

        g_xxyy_0_xxxxyy_0[i] = g_yy_0_xxxxyy_0[i] * fbe_0 - g_yy_0_xxxxyy_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyy_1[i] * fi_acd_0 + g_xyy_0_xxxxyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxyz_0[i] = g_yy_0_xxxxyz_0[i] * fbe_0 - g_yy_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyz_1[i] * fi_acd_0 + g_xyy_0_xxxxyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxzz_0[i] = g_xx_0_xxxxzz_0[i] * fbe_0 - g_xx_0_xxxxzz_1[i] * fz_be_0 + g_xxy_0_xxxxzz_1[i] * wa_y[i];

        g_xxyy_0_xxxyyy_0[i] = g_yy_0_xxxyyy_0[i] * fbe_0 - g_yy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyy_1[i] * fi_acd_0 + g_xyy_0_xxxyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxyyz_0[i] = g_yy_0_xxxyyz_0[i] * fbe_0 - g_yy_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyz_1[i] * fi_acd_0 + g_xyy_0_xxxyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxyzz_0[i] = g_yy_0_xxxyzz_0[i] * fbe_0 - g_yy_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyzz_1[i] * fi_acd_0 + g_xyy_0_xxxyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxzzz_0[i] = g_xx_0_xxxzzz_0[i] * fbe_0 - g_xx_0_xxxzzz_1[i] * fz_be_0 + g_xxy_0_xxxzzz_1[i] * wa_y[i];

        g_xxyy_0_xxyyyy_0[i] = g_yy_0_xxyyyy_0[i] * fbe_0 - g_yy_0_xxyyyy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyy_1[i] * fi_acd_0 + g_xyy_0_xxyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxyyyz_0[i] = g_yy_0_xxyyyz_0[i] * fbe_0 - g_yy_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyz_1[i] * fi_acd_0 + g_xyy_0_xxyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxyyzz_0[i] = g_yy_0_xxyyzz_0[i] * fbe_0 - g_yy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyzz_1[i] * fi_acd_0 + g_xyy_0_xxyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxyzzz_0[i] = g_yy_0_xxyzzz_0[i] * fbe_0 - g_yy_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyzzz_1[i] * fi_acd_0 + g_xyy_0_xxyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxzzzz_0[i] = g_xx_0_xxzzzz_0[i] * fbe_0 - g_xx_0_xxzzzz_1[i] * fz_be_0 + g_xxy_0_xxzzzz_1[i] * wa_y[i];

        g_xxyy_0_xyyyyy_0[i] = g_yy_0_xyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyy_1[i] * fi_acd_0 + g_xyy_0_xyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xyyyyz_0[i] = g_yy_0_xyyyyz_0[i] * fbe_0 - g_yy_0_xyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyz_1[i] * fi_acd_0 + g_xyy_0_xyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xyyyzz_0[i] = g_yy_0_xyyyzz_0[i] * fbe_0 - g_yy_0_xyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyzz_1[i] * fi_acd_0 + g_xyy_0_xyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xyyzzz_0[i] = g_yy_0_xyyzzz_0[i] * fbe_0 - g_yy_0_xyyzzz_1[i] * fz_be_0 + g_xyy_0_yyzzz_1[i] * fi_acd_0 + g_xyy_0_xyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xyzzzz_0[i] = g_yy_0_xyzzzz_0[i] * fbe_0 - g_yy_0_xyzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzz_1[i] * fi_acd_0 + g_xyy_0_xyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xzzzzz_0[i] = g_xx_0_xzzzzz_0[i] * fbe_0 - g_xx_0_xzzzzz_1[i] * fz_be_0 + g_xxy_0_xzzzzz_1[i] * wa_y[i];

        g_xxyy_0_yyyyyy_0[i] = g_yy_0_yyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyy_1[i] * wa_x[i];

        g_xxyy_0_yyyyyz_0[i] = g_yy_0_yyyyyz_0[i] * fbe_0 - g_yy_0_yyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyz_1[i] * wa_x[i];

        g_xxyy_0_yyyyzz_0[i] = g_yy_0_yyyyzz_0[i] * fbe_0 - g_yy_0_yyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyzz_1[i] * wa_x[i];

        g_xxyy_0_yyyzzz_0[i] = g_yy_0_yyyzzz_0[i] * fbe_0 - g_yy_0_yyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzz_1[i] * wa_x[i];

        g_xxyy_0_yyzzzz_0[i] = g_yy_0_yyzzzz_0[i] * fbe_0 - g_yy_0_yyzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzz_1[i] * wa_x[i];

        g_xxyy_0_yzzzzz_0[i] = g_yy_0_yzzzzz_0[i] * fbe_0 - g_yy_0_yzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzz_1[i] * wa_x[i];

        g_xxyy_0_zzzzzz_0[i] = g_yy_0_zzzzzz_0[i] * fbe_0 - g_yy_0_zzzzzz_1[i] * fz_be_0 + g_xyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 112-140 components of targeted buffer : GSI

    auto g_xxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 112);

    auto g_xxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 113);

    auto g_xxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 114);

    auto g_xxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 115);

    auto g_xxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 116);

    auto g_xxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 117);

    auto g_xxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 118);

    auto g_xxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 119);

    auto g_xxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 120);

    auto g_xxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 121);

    auto g_xxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 122);

    auto g_xxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 123);

    auto g_xxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 124);

    auto g_xxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 125);

    auto g_xxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 126);

    auto g_xxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 127);

    auto g_xxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 128);

    auto g_xxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 129);

    auto g_xxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 130);

    auto g_xxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 131);

    auto g_xxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 132);

    auto g_xxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 133);

    auto g_xxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 134);

    auto g_xxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 135);

    auto g_xxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 136);

    auto g_xxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 137);

    auto g_xxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 138);

    auto g_xxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 139);

    #pragma omp simd aligned(g_xxy_0_xxxxxy_1, g_xxy_0_xxxxyy_1, g_xxy_0_xxxyyy_1, g_xxy_0_xxyyyy_1, g_xxy_0_xyyyyy_1, g_xxy_0_yyyyyy_1, g_xxyz_0_xxxxxx_0, g_xxyz_0_xxxxxy_0, g_xxyz_0_xxxxxz_0, g_xxyz_0_xxxxyy_0, g_xxyz_0_xxxxyz_0, g_xxyz_0_xxxxzz_0, g_xxyz_0_xxxyyy_0, g_xxyz_0_xxxyyz_0, g_xxyz_0_xxxyzz_0, g_xxyz_0_xxxzzz_0, g_xxyz_0_xxyyyy_0, g_xxyz_0_xxyyyz_0, g_xxyz_0_xxyyzz_0, g_xxyz_0_xxyzzz_0, g_xxyz_0_xxzzzz_0, g_xxyz_0_xyyyyy_0, g_xxyz_0_xyyyyz_0, g_xxyz_0_xyyyzz_0, g_xxyz_0_xyyzzz_0, g_xxyz_0_xyzzzz_0, g_xxyz_0_xzzzzz_0, g_xxyz_0_yyyyyy_0, g_xxyz_0_yyyyyz_0, g_xxyz_0_yyyyzz_0, g_xxyz_0_yyyzzz_0, g_xxyz_0_yyzzzz_0, g_xxyz_0_yzzzzz_0, g_xxyz_0_zzzzzz_0, g_xxz_0_xxxxxx_1, g_xxz_0_xxxxxz_1, g_xxz_0_xxxxyz_1, g_xxz_0_xxxxz_1, g_xxz_0_xxxxzz_1, g_xxz_0_xxxyyz_1, g_xxz_0_xxxyz_1, g_xxz_0_xxxyzz_1, g_xxz_0_xxxzz_1, g_xxz_0_xxxzzz_1, g_xxz_0_xxyyyz_1, g_xxz_0_xxyyz_1, g_xxz_0_xxyyzz_1, g_xxz_0_xxyzz_1, g_xxz_0_xxyzzz_1, g_xxz_0_xxzzz_1, g_xxz_0_xxzzzz_1, g_xxz_0_xyyyyz_1, g_xxz_0_xyyyz_1, g_xxz_0_xyyyzz_1, g_xxz_0_xyyzz_1, g_xxz_0_xyyzzz_1, g_xxz_0_xyzzz_1, g_xxz_0_xyzzzz_1, g_xxz_0_xzzzz_1, g_xxz_0_xzzzzz_1, g_xxz_0_yyyyyz_1, g_xxz_0_yyyyz_1, g_xxz_0_yyyyzz_1, g_xxz_0_yyyzz_1, g_xxz_0_yyyzzz_1, g_xxz_0_yyzzz_1, g_xxz_0_yyzzzz_1, g_xxz_0_yzzzz_1, g_xxz_0_yzzzzz_1, g_xxz_0_zzzzz_1, g_xxz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxxxxx_0[i] = g_xxz_0_xxxxxx_1[i] * wa_y[i];

        g_xxyz_0_xxxxxy_0[i] = g_xxy_0_xxxxxy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxz_0[i] = g_xxz_0_xxxxxz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyy_0[i] = g_xxy_0_xxxxyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxyz_0[i] = g_xxz_0_xxxxz_1[i] * fi_acd_0 + g_xxz_0_xxxxyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxzz_0[i] = g_xxz_0_xxxxzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyy_0[i] = g_xxy_0_xxxyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxyyz_0[i] = 2.0 * g_xxz_0_xxxyz_1[i] * fi_acd_0 + g_xxz_0_xxxyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxyzz_0[i] = g_xxz_0_xxxzz_1[i] * fi_acd_0 + g_xxz_0_xxxyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxzzz_0[i] = g_xxz_0_xxxzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyy_0[i] = g_xxy_0_xxyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxyyyz_0[i] = 3.0 * g_xxz_0_xxyyz_1[i] * fi_acd_0 + g_xxz_0_xxyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxyyzz_0[i] = 2.0 * g_xxz_0_xxyzz_1[i] * fi_acd_0 + g_xxz_0_xxyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxyzzz_0[i] = g_xxz_0_xxzzz_1[i] * fi_acd_0 + g_xxz_0_xxyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxzzzz_0[i] = g_xxz_0_xxzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyy_0[i] = g_xxy_0_xyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xyyyyz_0[i] = 4.0 * g_xxz_0_xyyyz_1[i] * fi_acd_0 + g_xxz_0_xyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xyyyzz_0[i] = 3.0 * g_xxz_0_xyyzz_1[i] * fi_acd_0 + g_xxz_0_xyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xyyzzz_0[i] = 2.0 * g_xxz_0_xyzzz_1[i] * fi_acd_0 + g_xxz_0_xyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xyzzzz_0[i] = g_xxz_0_xzzzz_1[i] * fi_acd_0 + g_xxz_0_xyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xzzzzz_0[i] = g_xxz_0_xzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyy_0[i] = g_xxy_0_yyyyyy_1[i] * wa_z[i];

        g_xxyz_0_yyyyyz_0[i] = 5.0 * g_xxz_0_yyyyz_1[i] * fi_acd_0 + g_xxz_0_yyyyyz_1[i] * wa_y[i];

        g_xxyz_0_yyyyzz_0[i] = 4.0 * g_xxz_0_yyyzz_1[i] * fi_acd_0 + g_xxz_0_yyyyzz_1[i] * wa_y[i];

        g_xxyz_0_yyyzzz_0[i] = 3.0 * g_xxz_0_yyzzz_1[i] * fi_acd_0 + g_xxz_0_yyyzzz_1[i] * wa_y[i];

        g_xxyz_0_yyzzzz_0[i] = 2.0 * g_xxz_0_yzzzz_1[i] * fi_acd_0 + g_xxz_0_yyzzzz_1[i] * wa_y[i];

        g_xxyz_0_yzzzzz_0[i] = g_xxz_0_zzzzz_1[i] * fi_acd_0 + g_xxz_0_yzzzzz_1[i] * wa_y[i];

        g_xxyz_0_zzzzzz_0[i] = g_xxz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 140-168 components of targeted buffer : GSI

    auto g_xxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 140);

    auto g_xxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 141);

    auto g_xxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 142);

    auto g_xxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 143);

    auto g_xxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 144);

    auto g_xxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 145);

    auto g_xxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 146);

    auto g_xxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 147);

    auto g_xxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 148);

    auto g_xxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 149);

    auto g_xxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 150);

    auto g_xxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 151);

    auto g_xxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 152);

    auto g_xxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 153);

    auto g_xxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 154);

    auto g_xxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 155);

    auto g_xxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 156);

    auto g_xxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 157);

    auto g_xxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 158);

    auto g_xxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 159);

    auto g_xxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 160);

    auto g_xxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 161);

    auto g_xxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 162);

    auto g_xxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 163);

    auto g_xxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 164);

    auto g_xxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 165);

    auto g_xxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 166);

    auto g_xxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 167);

    #pragma omp simd aligned(g_xx_0_xxxxxx_0, g_xx_0_xxxxxx_1, g_xx_0_xxxxxy_0, g_xx_0_xxxxxy_1, g_xx_0_xxxxyy_0, g_xx_0_xxxxyy_1, g_xx_0_xxxyyy_0, g_xx_0_xxxyyy_1, g_xx_0_xxyyyy_0, g_xx_0_xxyyyy_1, g_xx_0_xyyyyy_0, g_xx_0_xyyyyy_1, g_xxz_0_xxxxxx_1, g_xxz_0_xxxxxy_1, g_xxz_0_xxxxyy_1, g_xxz_0_xxxyyy_1, g_xxz_0_xxyyyy_1, g_xxz_0_xyyyyy_1, g_xxzz_0_xxxxxx_0, g_xxzz_0_xxxxxy_0, g_xxzz_0_xxxxxz_0, g_xxzz_0_xxxxyy_0, g_xxzz_0_xxxxyz_0, g_xxzz_0_xxxxzz_0, g_xxzz_0_xxxyyy_0, g_xxzz_0_xxxyyz_0, g_xxzz_0_xxxyzz_0, g_xxzz_0_xxxzzz_0, g_xxzz_0_xxyyyy_0, g_xxzz_0_xxyyyz_0, g_xxzz_0_xxyyzz_0, g_xxzz_0_xxyzzz_0, g_xxzz_0_xxzzzz_0, g_xxzz_0_xyyyyy_0, g_xxzz_0_xyyyyz_0, g_xxzz_0_xyyyzz_0, g_xxzz_0_xyyzzz_0, g_xxzz_0_xyzzzz_0, g_xxzz_0_xzzzzz_0, g_xxzz_0_yyyyyy_0, g_xxzz_0_yyyyyz_0, g_xxzz_0_yyyyzz_0, g_xxzz_0_yyyzzz_0, g_xxzz_0_yyzzzz_0, g_xxzz_0_yzzzzz_0, g_xxzz_0_zzzzzz_0, g_xzz_0_xxxxxz_1, g_xzz_0_xxxxyz_1, g_xzz_0_xxxxz_1, g_xzz_0_xxxxzz_1, g_xzz_0_xxxyyz_1, g_xzz_0_xxxyz_1, g_xzz_0_xxxyzz_1, g_xzz_0_xxxzz_1, g_xzz_0_xxxzzz_1, g_xzz_0_xxyyyz_1, g_xzz_0_xxyyz_1, g_xzz_0_xxyyzz_1, g_xzz_0_xxyzz_1, g_xzz_0_xxyzzz_1, g_xzz_0_xxzzz_1, g_xzz_0_xxzzzz_1, g_xzz_0_xyyyyz_1, g_xzz_0_xyyyz_1, g_xzz_0_xyyyzz_1, g_xzz_0_xyyzz_1, g_xzz_0_xyyzzz_1, g_xzz_0_xyzzz_1, g_xzz_0_xyzzzz_1, g_xzz_0_xzzzz_1, g_xzz_0_xzzzzz_1, g_xzz_0_yyyyyy_1, g_xzz_0_yyyyyz_1, g_xzz_0_yyyyz_1, g_xzz_0_yyyyzz_1, g_xzz_0_yyyzz_1, g_xzz_0_yyyzzz_1, g_xzz_0_yyzzz_1, g_xzz_0_yyzzzz_1, g_xzz_0_yzzzz_1, g_xzz_0_yzzzzz_1, g_xzz_0_zzzzz_1, g_xzz_0_zzzzzz_1, g_zz_0_xxxxxz_0, g_zz_0_xxxxxz_1, g_zz_0_xxxxyz_0, g_zz_0_xxxxyz_1, g_zz_0_xxxxzz_0, g_zz_0_xxxxzz_1, g_zz_0_xxxyyz_0, g_zz_0_xxxyyz_1, g_zz_0_xxxyzz_0, g_zz_0_xxxyzz_1, g_zz_0_xxxzzz_0, g_zz_0_xxxzzz_1, g_zz_0_xxyyyz_0, g_zz_0_xxyyyz_1, g_zz_0_xxyyzz_0, g_zz_0_xxyyzz_1, g_zz_0_xxyzzz_0, g_zz_0_xxyzzz_1, g_zz_0_xxzzzz_0, g_zz_0_xxzzzz_1, g_zz_0_xyyyyz_0, g_zz_0_xyyyyz_1, g_zz_0_xyyyzz_0, g_zz_0_xyyyzz_1, g_zz_0_xyyzzz_0, g_zz_0_xyyzzz_1, g_zz_0_xyzzzz_0, g_zz_0_xyzzzz_1, g_zz_0_xzzzzz_0, g_zz_0_xzzzzz_1, g_zz_0_yyyyyy_0, g_zz_0_yyyyyy_1, g_zz_0_yyyyyz_0, g_zz_0_yyyyyz_1, g_zz_0_yyyyzz_0, g_zz_0_yyyyzz_1, g_zz_0_yyyzzz_0, g_zz_0_yyyzzz_1, g_zz_0_yyzzzz_0, g_zz_0_yyzzzz_1, g_zz_0_yzzzzz_0, g_zz_0_yzzzzz_1, g_zz_0_zzzzzz_0, g_zz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxxxxx_0[i] = g_xx_0_xxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxx_1[i] * fz_be_0 + g_xxz_0_xxxxxx_1[i] * wa_z[i];

        g_xxzz_0_xxxxxy_0[i] = g_xx_0_xxxxxy_0[i] * fbe_0 - g_xx_0_xxxxxy_1[i] * fz_be_0 + g_xxz_0_xxxxxy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxz_0[i] = g_zz_0_xxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxz_1[i] * fi_acd_0 + g_xzz_0_xxxxxz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyy_0[i] = g_xx_0_xxxxyy_0[i] * fbe_0 - g_xx_0_xxxxyy_1[i] * fz_be_0 + g_xxz_0_xxxxyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxyz_0[i] = g_zz_0_xxxxyz_0[i] * fbe_0 - g_zz_0_xxxxyz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyz_1[i] * fi_acd_0 + g_xzz_0_xxxxyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxzz_0[i] = g_zz_0_xxxxzz_0[i] * fbe_0 - g_zz_0_xxxxzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxzz_1[i] * fi_acd_0 + g_xzz_0_xxxxzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyy_0[i] = g_xx_0_xxxyyy_0[i] * fbe_0 - g_xx_0_xxxyyy_1[i] * fz_be_0 + g_xxz_0_xxxyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxyyz_0[i] = g_zz_0_xxxyyz_0[i] * fbe_0 - g_zz_0_xxxyyz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyz_1[i] * fi_acd_0 + g_xzz_0_xxxyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxyzz_0[i] = g_zz_0_xxxyzz_0[i] * fbe_0 - g_zz_0_xxxyzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyzz_1[i] * fi_acd_0 + g_xzz_0_xxxyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxzzz_0[i] = g_zz_0_xxxzzz_0[i] * fbe_0 - g_zz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxzzz_1[i] * fi_acd_0 + g_xzz_0_xxxzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyy_0[i] = g_xx_0_xxyyyy_0[i] * fbe_0 - g_xx_0_xxyyyy_1[i] * fz_be_0 + g_xxz_0_xxyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxyyyz_0[i] = g_zz_0_xxyyyz_0[i] * fbe_0 - g_zz_0_xxyyyz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyz_1[i] * fi_acd_0 + g_xzz_0_xxyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxyyzz_0[i] = g_zz_0_xxyyzz_0[i] * fbe_0 - g_zz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyzz_1[i] * fi_acd_0 + g_xzz_0_xxyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxyzzz_0[i] = g_zz_0_xxyzzz_0[i] * fbe_0 - g_zz_0_xxyzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyzzz_1[i] * fi_acd_0 + g_xzz_0_xxyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxzzzz_0[i] = g_zz_0_xxzzzz_0[i] * fbe_0 - g_zz_0_xxzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xzzzz_1[i] * fi_acd_0 + g_xzz_0_xxzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyy_0[i] = g_xx_0_xyyyyy_0[i] * fbe_0 - g_xx_0_xyyyyy_1[i] * fz_be_0 + g_xxz_0_xyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xyyyyz_0[i] = g_zz_0_xyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyz_1[i] * fi_acd_0 + g_xzz_0_xyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xyyyzz_0[i] = g_zz_0_xyyyzz_0[i] * fbe_0 - g_zz_0_xyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyzz_1[i] * fi_acd_0 + g_xzz_0_xyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xyyzzz_0[i] = g_zz_0_xyyzzz_0[i] * fbe_0 - g_zz_0_xyyzzz_1[i] * fz_be_0 + g_xzz_0_yyzzz_1[i] * fi_acd_0 + g_xzz_0_xyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xyzzzz_0[i] = g_zz_0_xyzzzz_0[i] * fbe_0 - g_zz_0_xyzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzz_1[i] * fi_acd_0 + g_xzz_0_xyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xzzzzz_0[i] = g_zz_0_xzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzz_1[i] * fi_acd_0 + g_xzz_0_xzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyy_0[i] = g_zz_0_yyyyyy_0[i] * fbe_0 - g_zz_0_yyyyyy_1[i] * fz_be_0 + g_xzz_0_yyyyyy_1[i] * wa_x[i];

        g_xxzz_0_yyyyyz_0[i] = g_zz_0_yyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyz_1[i] * wa_x[i];

        g_xxzz_0_yyyyzz_0[i] = g_zz_0_yyyyzz_0[i] * fbe_0 - g_zz_0_yyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyzz_1[i] * wa_x[i];

        g_xxzz_0_yyyzzz_0[i] = g_zz_0_yyyzzz_0[i] * fbe_0 - g_zz_0_yyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzz_1[i] * wa_x[i];

        g_xxzz_0_yyzzzz_0[i] = g_zz_0_yyzzzz_0[i] * fbe_0 - g_zz_0_yyzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzz_1[i] * wa_x[i];

        g_xxzz_0_yzzzzz_0[i] = g_zz_0_yzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzz_1[i] * wa_x[i];

        g_xxzz_0_zzzzzz_0[i] = g_zz_0_zzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 168-196 components of targeted buffer : GSI

    auto g_xyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 168);

    auto g_xyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 169);

    auto g_xyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 170);

    auto g_xyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 171);

    auto g_xyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 172);

    auto g_xyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 173);

    auto g_xyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 174);

    auto g_xyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 175);

    auto g_xyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 176);

    auto g_xyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 177);

    auto g_xyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 178);

    auto g_xyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 179);

    auto g_xyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 180);

    auto g_xyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 181);

    auto g_xyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 182);

    auto g_xyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 183);

    auto g_xyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 184);

    auto g_xyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 185);

    auto g_xyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 186);

    auto g_xyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 187);

    auto g_xyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 188);

    auto g_xyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 189);

    auto g_xyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 190);

    auto g_xyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 191);

    auto g_xyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 192);

    auto g_xyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 193);

    auto g_xyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 194);

    auto g_xyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 195);

    #pragma omp simd aligned(g_xyyy_0_xxxxxx_0, g_xyyy_0_xxxxxy_0, g_xyyy_0_xxxxxz_0, g_xyyy_0_xxxxyy_0, g_xyyy_0_xxxxyz_0, g_xyyy_0_xxxxzz_0, g_xyyy_0_xxxyyy_0, g_xyyy_0_xxxyyz_0, g_xyyy_0_xxxyzz_0, g_xyyy_0_xxxzzz_0, g_xyyy_0_xxyyyy_0, g_xyyy_0_xxyyyz_0, g_xyyy_0_xxyyzz_0, g_xyyy_0_xxyzzz_0, g_xyyy_0_xxzzzz_0, g_xyyy_0_xyyyyy_0, g_xyyy_0_xyyyyz_0, g_xyyy_0_xyyyzz_0, g_xyyy_0_xyyzzz_0, g_xyyy_0_xyzzzz_0, g_xyyy_0_xzzzzz_0, g_xyyy_0_yyyyyy_0, g_xyyy_0_yyyyyz_0, g_xyyy_0_yyyyzz_0, g_xyyy_0_yyyzzz_0, g_xyyy_0_yyzzzz_0, g_xyyy_0_yzzzzz_0, g_xyyy_0_zzzzzz_0, g_yyy_0_xxxxx_1, g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxy_1, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxz_1, g_yyy_0_xxxxzz_1, g_yyy_0_xxxyy_1, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyz_1, g_yyy_0_xxxyzz_1, g_yyy_0_xxxzz_1, g_yyy_0_xxxzzz_1, g_yyy_0_xxyyy_1, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyz_1, g_yyy_0_xxyyzz_1, g_yyy_0_xxyzz_1, g_yyy_0_xxyzzz_1, g_yyy_0_xxzzz_1, g_yyy_0_xxzzzz_1, g_yyy_0_xyyyy_1, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyz_1, g_yyy_0_xyyyzz_1, g_yyy_0_xyyzz_1, g_yyy_0_xyyzzz_1, g_yyy_0_xyzzz_1, g_yyy_0_xyzzzz_1, g_yyy_0_xzzzz_1, g_yyy_0_xzzzzz_1, g_yyy_0_yyyyy_1, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyz_1, g_yyy_0_yyyyzz_1, g_yyy_0_yyyzz_1, g_yyy_0_yyyzzz_1, g_yyy_0_yyzzz_1, g_yyy_0_yyzzzz_1, g_yyy_0_yzzzz_1, g_yyy_0_yzzzzz_1, g_yyy_0_zzzzz_1, g_yyy_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxxxxx_0[i] = 6.0 * g_yyy_0_xxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxx_1[i] * wa_x[i];

        g_xyyy_0_xxxxxy_0[i] = 5.0 * g_yyy_0_xxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxz_0[i] = 5.0 * g_yyy_0_xxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyy_0[i] = 4.0 * g_yyy_0_xxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxyz_0[i] = 4.0 * g_yyy_0_xxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxzz_0[i] = 4.0 * g_yyy_0_xxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyy_0[i] = 3.0 * g_yyy_0_xxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxyyz_0[i] = 3.0 * g_yyy_0_xxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxyzz_0[i] = 3.0 * g_yyy_0_xxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxzzz_0[i] = 3.0 * g_yyy_0_xxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyy_0[i] = 2.0 * g_yyy_0_xyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxyyyz_0[i] = 2.0 * g_yyy_0_xyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxyyzz_0[i] = 2.0 * g_yyy_0_xyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxyzzz_0[i] = 2.0 * g_yyy_0_xyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxzzzz_0[i] = 2.0 * g_yyy_0_xzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyy_0[i] = g_yyy_0_yyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xyyyyz_0[i] = g_yyy_0_yyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xyyyzz_0[i] = g_yyy_0_yyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xyyzzz_0[i] = g_yyy_0_yyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xyzzzz_0[i] = g_yyy_0_yzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xzzzzz_0[i] = g_yyy_0_zzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyy_0[i] = g_yyy_0_yyyyyy_1[i] * wa_x[i];

        g_xyyy_0_yyyyyz_0[i] = g_yyy_0_yyyyyz_1[i] * wa_x[i];

        g_xyyy_0_yyyyzz_0[i] = g_yyy_0_yyyyzz_1[i] * wa_x[i];

        g_xyyy_0_yyyzzz_0[i] = g_yyy_0_yyyzzz_1[i] * wa_x[i];

        g_xyyy_0_yyzzzz_0[i] = g_yyy_0_yyzzzz_1[i] * wa_x[i];

        g_xyyy_0_yzzzzz_0[i] = g_yyy_0_yzzzzz_1[i] * wa_x[i];

        g_xyyy_0_zzzzzz_0[i] = g_yyy_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 196-224 components of targeted buffer : GSI

    auto g_xyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 196);

    auto g_xyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 197);

    auto g_xyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 198);

    auto g_xyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 199);

    auto g_xyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 200);

    auto g_xyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 201);

    auto g_xyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 202);

    auto g_xyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 203);

    auto g_xyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 204);

    auto g_xyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 205);

    auto g_xyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 206);

    auto g_xyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 207);

    auto g_xyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 208);

    auto g_xyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 209);

    auto g_xyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 210);

    auto g_xyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 211);

    auto g_xyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 212);

    auto g_xyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 213);

    auto g_xyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 214);

    auto g_xyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 215);

    auto g_xyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 216);

    auto g_xyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 217);

    auto g_xyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 218);

    auto g_xyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 219);

    auto g_xyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 220);

    auto g_xyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 221);

    auto g_xyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 222);

    auto g_xyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 223);

    #pragma omp simd aligned(g_xyy_0_xxxxxx_1, g_xyy_0_xxxxxy_1, g_xyy_0_xxxxyy_1, g_xyy_0_xxxyyy_1, g_xyy_0_xxyyyy_1, g_xyy_0_xyyyyy_1, g_xyyz_0_xxxxxx_0, g_xyyz_0_xxxxxy_0, g_xyyz_0_xxxxxz_0, g_xyyz_0_xxxxyy_0, g_xyyz_0_xxxxyz_0, g_xyyz_0_xxxxzz_0, g_xyyz_0_xxxyyy_0, g_xyyz_0_xxxyyz_0, g_xyyz_0_xxxyzz_0, g_xyyz_0_xxxzzz_0, g_xyyz_0_xxyyyy_0, g_xyyz_0_xxyyyz_0, g_xyyz_0_xxyyzz_0, g_xyyz_0_xxyzzz_0, g_xyyz_0_xxzzzz_0, g_xyyz_0_xyyyyy_0, g_xyyz_0_xyyyyz_0, g_xyyz_0_xyyyzz_0, g_xyyz_0_xyyzzz_0, g_xyyz_0_xyzzzz_0, g_xyyz_0_xzzzzz_0, g_xyyz_0_yyyyyy_0, g_xyyz_0_yyyyyz_0, g_xyyz_0_yyyyzz_0, g_xyyz_0_yyyzzz_0, g_xyyz_0_yyzzzz_0, g_xyyz_0_yzzzzz_0, g_xyyz_0_zzzzzz_0, g_yyz_0_xxxxxz_1, g_yyz_0_xxxxyz_1, g_yyz_0_xxxxz_1, g_yyz_0_xxxxzz_1, g_yyz_0_xxxyyz_1, g_yyz_0_xxxyz_1, g_yyz_0_xxxyzz_1, g_yyz_0_xxxzz_1, g_yyz_0_xxxzzz_1, g_yyz_0_xxyyyz_1, g_yyz_0_xxyyz_1, g_yyz_0_xxyyzz_1, g_yyz_0_xxyzz_1, g_yyz_0_xxyzzz_1, g_yyz_0_xxzzz_1, g_yyz_0_xxzzzz_1, g_yyz_0_xyyyyz_1, g_yyz_0_xyyyz_1, g_yyz_0_xyyyzz_1, g_yyz_0_xyyzz_1, g_yyz_0_xyyzzz_1, g_yyz_0_xyzzz_1, g_yyz_0_xyzzzz_1, g_yyz_0_xzzzz_1, g_yyz_0_xzzzzz_1, g_yyz_0_yyyyyy_1, g_yyz_0_yyyyyz_1, g_yyz_0_yyyyz_1, g_yyz_0_yyyyzz_1, g_yyz_0_yyyzz_1, g_yyz_0_yyyzzz_1, g_yyz_0_yyzzz_1, g_yyz_0_yyzzzz_1, g_yyz_0_yzzzz_1, g_yyz_0_yzzzzz_1, g_yyz_0_zzzzz_1, g_yyz_0_zzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxxxxx_0[i] = g_xyy_0_xxxxxx_1[i] * wa_z[i];

        g_xyyz_0_xxxxxy_0[i] = g_xyy_0_xxxxxy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxz_0[i] = 5.0 * g_yyz_0_xxxxz_1[i] * fi_acd_0 + g_yyz_0_xxxxxz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyy_0[i] = g_xyy_0_xxxxyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxyz_0[i] = 4.0 * g_yyz_0_xxxyz_1[i] * fi_acd_0 + g_yyz_0_xxxxyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxzz_0[i] = 4.0 * g_yyz_0_xxxzz_1[i] * fi_acd_0 + g_yyz_0_xxxxzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyy_0[i] = g_xyy_0_xxxyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxyyz_0[i] = 3.0 * g_yyz_0_xxyyz_1[i] * fi_acd_0 + g_yyz_0_xxxyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxyzz_0[i] = 3.0 * g_yyz_0_xxyzz_1[i] * fi_acd_0 + g_yyz_0_xxxyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxzzz_0[i] = 3.0 * g_yyz_0_xxzzz_1[i] * fi_acd_0 + g_yyz_0_xxxzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyy_0[i] = g_xyy_0_xxyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxyyyz_0[i] = 2.0 * g_yyz_0_xyyyz_1[i] * fi_acd_0 + g_yyz_0_xxyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxyyzz_0[i] = 2.0 * g_yyz_0_xyyzz_1[i] * fi_acd_0 + g_yyz_0_xxyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxyzzz_0[i] = 2.0 * g_yyz_0_xyzzz_1[i] * fi_acd_0 + g_yyz_0_xxyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxzzzz_0[i] = 2.0 * g_yyz_0_xzzzz_1[i] * fi_acd_0 + g_yyz_0_xxzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyy_0[i] = g_xyy_0_xyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xyyyyz_0[i] = g_yyz_0_yyyyz_1[i] * fi_acd_0 + g_yyz_0_xyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xyyyzz_0[i] = g_yyz_0_yyyzz_1[i] * fi_acd_0 + g_yyz_0_xyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xyyzzz_0[i] = g_yyz_0_yyzzz_1[i] * fi_acd_0 + g_yyz_0_xyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xyzzzz_0[i] = g_yyz_0_yzzzz_1[i] * fi_acd_0 + g_yyz_0_xyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xzzzzz_0[i] = g_yyz_0_zzzzz_1[i] * fi_acd_0 + g_yyz_0_xzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyy_0[i] = g_yyz_0_yyyyyy_1[i] * wa_x[i];

        g_xyyz_0_yyyyyz_0[i] = g_yyz_0_yyyyyz_1[i] * wa_x[i];

        g_xyyz_0_yyyyzz_0[i] = g_yyz_0_yyyyzz_1[i] * wa_x[i];

        g_xyyz_0_yyyzzz_0[i] = g_yyz_0_yyyzzz_1[i] * wa_x[i];

        g_xyyz_0_yyzzzz_0[i] = g_yyz_0_yyzzzz_1[i] * wa_x[i];

        g_xyyz_0_yzzzzz_0[i] = g_yyz_0_yzzzzz_1[i] * wa_x[i];

        g_xyyz_0_zzzzzz_0[i] = g_yyz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 224-252 components of targeted buffer : GSI

    auto g_xyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 224);

    auto g_xyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 225);

    auto g_xyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 226);

    auto g_xyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 227);

    auto g_xyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 228);

    auto g_xyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 229);

    auto g_xyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 230);

    auto g_xyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 231);

    auto g_xyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 232);

    auto g_xyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 233);

    auto g_xyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 234);

    auto g_xyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 235);

    auto g_xyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 236);

    auto g_xyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 237);

    auto g_xyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 238);

    auto g_xyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 239);

    auto g_xyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 240);

    auto g_xyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 241);

    auto g_xyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 242);

    auto g_xyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 243);

    auto g_xyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 244);

    auto g_xyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 245);

    auto g_xyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 246);

    auto g_xyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 247);

    auto g_xyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 248);

    auto g_xyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 249);

    auto g_xyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 250);

    auto g_xyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 251);

    #pragma omp simd aligned(g_xyzz_0_xxxxxx_0, g_xyzz_0_xxxxxy_0, g_xyzz_0_xxxxxz_0, g_xyzz_0_xxxxyy_0, g_xyzz_0_xxxxyz_0, g_xyzz_0_xxxxzz_0, g_xyzz_0_xxxyyy_0, g_xyzz_0_xxxyyz_0, g_xyzz_0_xxxyzz_0, g_xyzz_0_xxxzzz_0, g_xyzz_0_xxyyyy_0, g_xyzz_0_xxyyyz_0, g_xyzz_0_xxyyzz_0, g_xyzz_0_xxyzzz_0, g_xyzz_0_xxzzzz_0, g_xyzz_0_xyyyyy_0, g_xyzz_0_xyyyyz_0, g_xyzz_0_xyyyzz_0, g_xyzz_0_xyyzzz_0, g_xyzz_0_xyzzzz_0, g_xyzz_0_xzzzzz_0, g_xyzz_0_yyyyyy_0, g_xyzz_0_yyyyyz_0, g_xyzz_0_yyyyzz_0, g_xyzz_0_yyyzzz_0, g_xyzz_0_yyzzzz_0, g_xyzz_0_yzzzzz_0, g_xyzz_0_zzzzzz_0, g_xzz_0_xxxxxx_1, g_xzz_0_xxxxxz_1, g_xzz_0_xxxxzz_1, g_xzz_0_xxxzzz_1, g_xzz_0_xxzzzz_1, g_xzz_0_xzzzzz_1, g_yzz_0_xxxxxy_1, g_yzz_0_xxxxy_1, g_yzz_0_xxxxyy_1, g_yzz_0_xxxxyz_1, g_yzz_0_xxxyy_1, g_yzz_0_xxxyyy_1, g_yzz_0_xxxyyz_1, g_yzz_0_xxxyz_1, g_yzz_0_xxxyzz_1, g_yzz_0_xxyyy_1, g_yzz_0_xxyyyy_1, g_yzz_0_xxyyyz_1, g_yzz_0_xxyyz_1, g_yzz_0_xxyyzz_1, g_yzz_0_xxyzz_1, g_yzz_0_xxyzzz_1, g_yzz_0_xyyyy_1, g_yzz_0_xyyyyy_1, g_yzz_0_xyyyyz_1, g_yzz_0_xyyyz_1, g_yzz_0_xyyyzz_1, g_yzz_0_xyyzz_1, g_yzz_0_xyyzzz_1, g_yzz_0_xyzzz_1, g_yzz_0_xyzzzz_1, g_yzz_0_yyyyy_1, g_yzz_0_yyyyyy_1, g_yzz_0_yyyyyz_1, g_yzz_0_yyyyz_1, g_yzz_0_yyyyzz_1, g_yzz_0_yyyzz_1, g_yzz_0_yyyzzz_1, g_yzz_0_yyzzz_1, g_yzz_0_yyzzzz_1, g_yzz_0_yzzzz_1, g_yzz_0_yzzzzz_1, g_yzz_0_zzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxxxxx_0[i] = g_xzz_0_xxxxxx_1[i] * wa_y[i];

        g_xyzz_0_xxxxxy_0[i] = 5.0 * g_yzz_0_xxxxy_1[i] * fi_acd_0 + g_yzz_0_xxxxxy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxz_0[i] = g_xzz_0_xxxxxz_1[i] * wa_y[i];

        g_xyzz_0_xxxxyy_0[i] = 4.0 * g_yzz_0_xxxyy_1[i] * fi_acd_0 + g_yzz_0_xxxxyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxyz_0[i] = 4.0 * g_yzz_0_xxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxzz_0[i] = g_xzz_0_xxxxzz_1[i] * wa_y[i];

        g_xyzz_0_xxxyyy_0[i] = 3.0 * g_yzz_0_xxyyy_1[i] * fi_acd_0 + g_yzz_0_xxxyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxyyz_0[i] = 3.0 * g_yzz_0_xxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxyzz_0[i] = 3.0 * g_yzz_0_xxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxzzz_0[i] = g_xzz_0_xxxzzz_1[i] * wa_y[i];

        g_xyzz_0_xxyyyy_0[i] = 2.0 * g_yzz_0_xyyyy_1[i] * fi_acd_0 + g_yzz_0_xxyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxyyyz_0[i] = 2.0 * g_yzz_0_xyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxyyzz_0[i] = 2.0 * g_yzz_0_xyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxyzzz_0[i] = 2.0 * g_yzz_0_xyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxzzzz_0[i] = g_xzz_0_xxzzzz_1[i] * wa_y[i];

        g_xyzz_0_xyyyyy_0[i] = g_yzz_0_yyyyy_1[i] * fi_acd_0 + g_yzz_0_xyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xyyyyz_0[i] = g_yzz_0_yyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xyyyzz_0[i] = g_yzz_0_yyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xyyzzz_0[i] = g_yzz_0_yyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xyzzzz_0[i] = g_yzz_0_yzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xzzzzz_0[i] = g_xzz_0_xzzzzz_1[i] * wa_y[i];

        g_xyzz_0_yyyyyy_0[i] = g_yzz_0_yyyyyy_1[i] * wa_x[i];

        g_xyzz_0_yyyyyz_0[i] = g_yzz_0_yyyyyz_1[i] * wa_x[i];

        g_xyzz_0_yyyyzz_0[i] = g_yzz_0_yyyyzz_1[i] * wa_x[i];

        g_xyzz_0_yyyzzz_0[i] = g_yzz_0_yyyzzz_1[i] * wa_x[i];

        g_xyzz_0_yyzzzz_0[i] = g_yzz_0_yyzzzz_1[i] * wa_x[i];

        g_xyzz_0_yzzzzz_0[i] = g_yzz_0_yzzzzz_1[i] * wa_x[i];

        g_xyzz_0_zzzzzz_0[i] = g_yzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 252-280 components of targeted buffer : GSI

    auto g_xzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 252);

    auto g_xzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 253);

    auto g_xzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 254);

    auto g_xzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 255);

    auto g_xzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 256);

    auto g_xzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 257);

    auto g_xzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 258);

    auto g_xzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 259);

    auto g_xzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 260);

    auto g_xzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 261);

    auto g_xzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 262);

    auto g_xzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 263);

    auto g_xzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 264);

    auto g_xzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 265);

    auto g_xzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 266);

    auto g_xzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 267);

    auto g_xzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 268);

    auto g_xzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 269);

    auto g_xzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 270);

    auto g_xzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 271);

    auto g_xzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 272);

    auto g_xzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 273);

    auto g_xzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 274);

    auto g_xzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 275);

    auto g_xzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 276);

    auto g_xzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 277);

    auto g_xzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 278);

    auto g_xzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 279);

    #pragma omp simd aligned(g_xzzz_0_xxxxxx_0, g_xzzz_0_xxxxxy_0, g_xzzz_0_xxxxxz_0, g_xzzz_0_xxxxyy_0, g_xzzz_0_xxxxyz_0, g_xzzz_0_xxxxzz_0, g_xzzz_0_xxxyyy_0, g_xzzz_0_xxxyyz_0, g_xzzz_0_xxxyzz_0, g_xzzz_0_xxxzzz_0, g_xzzz_0_xxyyyy_0, g_xzzz_0_xxyyyz_0, g_xzzz_0_xxyyzz_0, g_xzzz_0_xxyzzz_0, g_xzzz_0_xxzzzz_0, g_xzzz_0_xyyyyy_0, g_xzzz_0_xyyyyz_0, g_xzzz_0_xyyyzz_0, g_xzzz_0_xyyzzz_0, g_xzzz_0_xyzzzz_0, g_xzzz_0_xzzzzz_0, g_xzzz_0_yyyyyy_0, g_xzzz_0_yyyyyz_0, g_xzzz_0_yyyyzz_0, g_xzzz_0_yyyzzz_0, g_xzzz_0_yyzzzz_0, g_xzzz_0_yzzzzz_0, g_xzzz_0_zzzzzz_0, g_zzz_0_xxxxx_1, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxy_1, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxz_1, g_zzz_0_xxxxzz_1, g_zzz_0_xxxyy_1, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyz_1, g_zzz_0_xxxyzz_1, g_zzz_0_xxxzz_1, g_zzz_0_xxxzzz_1, g_zzz_0_xxyyy_1, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyz_1, g_zzz_0_xxyyzz_1, g_zzz_0_xxyzz_1, g_zzz_0_xxyzzz_1, g_zzz_0_xxzzz_1, g_zzz_0_xxzzzz_1, g_zzz_0_xyyyy_1, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyz_1, g_zzz_0_xyyyzz_1, g_zzz_0_xyyzz_1, g_zzz_0_xyyzzz_1, g_zzz_0_xyzzz_1, g_zzz_0_xyzzzz_1, g_zzz_0_xzzzz_1, g_zzz_0_xzzzzz_1, g_zzz_0_yyyyy_1, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyz_1, g_zzz_0_yyyyzz_1, g_zzz_0_yyyzz_1, g_zzz_0_yyyzzz_1, g_zzz_0_yyzzz_1, g_zzz_0_yyzzzz_1, g_zzz_0_yzzzz_1, g_zzz_0_yzzzzz_1, g_zzz_0_zzzzz_1, g_zzz_0_zzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxxxxx_0[i] = 6.0 * g_zzz_0_xxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxx_1[i] * wa_x[i];

        g_xzzz_0_xxxxxy_0[i] = 5.0 * g_zzz_0_xxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxz_0[i] = 5.0 * g_zzz_0_xxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyy_0[i] = 4.0 * g_zzz_0_xxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxyz_0[i] = 4.0 * g_zzz_0_xxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxzz_0[i] = 4.0 * g_zzz_0_xxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyy_0[i] = 3.0 * g_zzz_0_xxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxyyz_0[i] = 3.0 * g_zzz_0_xxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxyzz_0[i] = 3.0 * g_zzz_0_xxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxzzz_0[i] = 3.0 * g_zzz_0_xxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyy_0[i] = 2.0 * g_zzz_0_xyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxyyyz_0[i] = 2.0 * g_zzz_0_xyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxyyzz_0[i] = 2.0 * g_zzz_0_xyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxyzzz_0[i] = 2.0 * g_zzz_0_xyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxzzzz_0[i] = 2.0 * g_zzz_0_xzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyy_0[i] = g_zzz_0_yyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xyyyyz_0[i] = g_zzz_0_yyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xyyyzz_0[i] = g_zzz_0_yyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xyyzzz_0[i] = g_zzz_0_yyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xyzzzz_0[i] = g_zzz_0_yzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xzzzzz_0[i] = g_zzz_0_zzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyy_0[i] = g_zzz_0_yyyyyy_1[i] * wa_x[i];

        g_xzzz_0_yyyyyz_0[i] = g_zzz_0_yyyyyz_1[i] * wa_x[i];

        g_xzzz_0_yyyyzz_0[i] = g_zzz_0_yyyyzz_1[i] * wa_x[i];

        g_xzzz_0_yyyzzz_0[i] = g_zzz_0_yyyzzz_1[i] * wa_x[i];

        g_xzzz_0_yyzzzz_0[i] = g_zzz_0_yyzzzz_1[i] * wa_x[i];

        g_xzzz_0_yzzzzz_0[i] = g_zzz_0_yzzzzz_1[i] * wa_x[i];

        g_xzzz_0_zzzzzz_0[i] = g_zzz_0_zzzzzz_1[i] * wa_x[i];
    }

    /// Set up 280-308 components of targeted buffer : GSI

    auto g_yyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 280);

    auto g_yyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 281);

    auto g_yyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 282);

    auto g_yyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 283);

    auto g_yyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 284);

    auto g_yyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 285);

    auto g_yyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 286);

    auto g_yyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 287);

    auto g_yyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 288);

    auto g_yyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 289);

    auto g_yyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 290);

    auto g_yyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 291);

    auto g_yyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 292);

    auto g_yyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 293);

    auto g_yyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 294);

    auto g_yyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 295);

    auto g_yyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 296);

    auto g_yyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 297);

    auto g_yyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 298);

    auto g_yyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 299);

    auto g_yyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 300);

    auto g_yyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 301);

    auto g_yyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 302);

    auto g_yyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 303);

    auto g_yyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 304);

    auto g_yyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 305);

    auto g_yyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 306);

    auto g_yyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 307);

    #pragma omp simd aligned(g_yy_0_xxxxxx_0, g_yy_0_xxxxxx_1, g_yy_0_xxxxxy_0, g_yy_0_xxxxxy_1, g_yy_0_xxxxxz_0, g_yy_0_xxxxxz_1, g_yy_0_xxxxyy_0, g_yy_0_xxxxyy_1, g_yy_0_xxxxyz_0, g_yy_0_xxxxyz_1, g_yy_0_xxxxzz_0, g_yy_0_xxxxzz_1, g_yy_0_xxxyyy_0, g_yy_0_xxxyyy_1, g_yy_0_xxxyyz_0, g_yy_0_xxxyyz_1, g_yy_0_xxxyzz_0, g_yy_0_xxxyzz_1, g_yy_0_xxxzzz_0, g_yy_0_xxxzzz_1, g_yy_0_xxyyyy_0, g_yy_0_xxyyyy_1, g_yy_0_xxyyyz_0, g_yy_0_xxyyyz_1, g_yy_0_xxyyzz_0, g_yy_0_xxyyzz_1, g_yy_0_xxyzzz_0, g_yy_0_xxyzzz_1, g_yy_0_xxzzzz_0, g_yy_0_xxzzzz_1, g_yy_0_xyyyyy_0, g_yy_0_xyyyyy_1, g_yy_0_xyyyyz_0, g_yy_0_xyyyyz_1, g_yy_0_xyyyzz_0, g_yy_0_xyyyzz_1, g_yy_0_xyyzzz_0, g_yy_0_xyyzzz_1, g_yy_0_xyzzzz_0, g_yy_0_xyzzzz_1, g_yy_0_xzzzzz_0, g_yy_0_xzzzzz_1, g_yy_0_yyyyyy_0, g_yy_0_yyyyyy_1, g_yy_0_yyyyyz_0, g_yy_0_yyyyyz_1, g_yy_0_yyyyzz_0, g_yy_0_yyyyzz_1, g_yy_0_yyyzzz_0, g_yy_0_yyyzzz_1, g_yy_0_yyzzzz_0, g_yy_0_yyzzzz_1, g_yy_0_yzzzzz_0, g_yy_0_yzzzzz_1, g_yy_0_zzzzzz_0, g_yy_0_zzzzzz_1, g_yyy_0_xxxxx_1, g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxy_1, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxz_1, g_yyy_0_xxxxzz_1, g_yyy_0_xxxyy_1, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyz_1, g_yyy_0_xxxyzz_1, g_yyy_0_xxxzz_1, g_yyy_0_xxxzzz_1, g_yyy_0_xxyyy_1, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyz_1, g_yyy_0_xxyyzz_1, g_yyy_0_xxyzz_1, g_yyy_0_xxyzzz_1, g_yyy_0_xxzzz_1, g_yyy_0_xxzzzz_1, g_yyy_0_xyyyy_1, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyz_1, g_yyy_0_xyyyzz_1, g_yyy_0_xyyzz_1, g_yyy_0_xyyzzz_1, g_yyy_0_xyzzz_1, g_yyy_0_xyzzzz_1, g_yyy_0_xzzzz_1, g_yyy_0_xzzzzz_1, g_yyy_0_yyyyy_1, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyz_1, g_yyy_0_yyyyzz_1, g_yyy_0_yyyzz_1, g_yyy_0_yyyzzz_1, g_yyy_0_yyzzz_1, g_yyy_0_yyzzzz_1, g_yyy_0_yzzzz_1, g_yyy_0_yzzzzz_1, g_yyy_0_zzzzz_1, g_yyy_0_zzzzzz_1, g_yyyy_0_xxxxxx_0, g_yyyy_0_xxxxxy_0, g_yyyy_0_xxxxxz_0, g_yyyy_0_xxxxyy_0, g_yyyy_0_xxxxyz_0, g_yyyy_0_xxxxzz_0, g_yyyy_0_xxxyyy_0, g_yyyy_0_xxxyyz_0, g_yyyy_0_xxxyzz_0, g_yyyy_0_xxxzzz_0, g_yyyy_0_xxyyyy_0, g_yyyy_0_xxyyyz_0, g_yyyy_0_xxyyzz_0, g_yyyy_0_xxyzzz_0, g_yyyy_0_xxzzzz_0, g_yyyy_0_xyyyyy_0, g_yyyy_0_xyyyyz_0, g_yyyy_0_xyyyzz_0, g_yyyy_0_xyyzzz_0, g_yyyy_0_xyzzzz_0, g_yyyy_0_xzzzzz_0, g_yyyy_0_yyyyyy_0, g_yyyy_0_yyyyyz_0, g_yyyy_0_yyyyzz_0, g_yyyy_0_yyyzzz_0, g_yyyy_0_yyzzzz_0, g_yyyy_0_yzzzzz_0, g_yyyy_0_zzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxxxxx_0[i] = 3.0 * g_yy_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxx_1[i] * fz_be_0 + g_yyy_0_xxxxxx_1[i] * wa_y[i];

        g_yyyy_0_xxxxxy_0[i] = 3.0 * g_yy_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxy_1[i] * fz_be_0 + g_yyy_0_xxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxz_0[i] = 3.0 * g_yy_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxz_1[i] * fz_be_0 + g_yyy_0_xxxxxz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyy_0[i] = 3.0 * g_yy_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxyz_0[i] = 3.0 * g_yy_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyz_1[i] * fz_be_0 + g_yyy_0_xxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxzz_0[i] = 3.0 * g_yy_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxzz_1[i] * fz_be_0 + g_yyy_0_xxxxzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyy_0[i] = 3.0 * g_yy_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxyyz_0[i] = 3.0 * g_yy_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxyzz_0[i] = 3.0 * g_yy_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyzz_1[i] * fz_be_0 + g_yyy_0_xxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxzzz_0[i] = 3.0 * g_yy_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyy_0[i] = 3.0 * g_yy_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxyyyz_0[i] = 3.0 * g_yy_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxyyzz_0[i] = 3.0 * g_yy_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxyzzz_0[i] = 3.0 * g_yy_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyzzz_1[i] * fz_be_0 + g_yyy_0_xxzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxzzzz_0[i] = 3.0 * g_yy_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyy_0[i] = 3.0 * g_yy_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_0_xyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xyyyyz_0[i] = 3.0 * g_yy_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xyyyzz_0[i] = 3.0 * g_yy_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xyyzzz_0[i] = 3.0 * g_yy_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xyzzzz_0[i] = 3.0 * g_yy_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xzzzzz_0[i] = 3.0 * g_yy_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyy_0[i] = 3.0 * g_yy_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyy_1[i] * fz_be_0 + 6.0 * g_yyy_0_yyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyy_1[i] * wa_y[i];

        g_yyyy_0_yyyyyz_0[i] = 3.0 * g_yy_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yyy_0_yyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyz_1[i] * wa_y[i];

        g_yyyy_0_yyyyzz_0[i] = 3.0 * g_yy_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_yyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzz_1[i] * wa_y[i];

        g_yyyy_0_yyyzzz_0[i] = 3.0 * g_yy_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_yyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzz_1[i] * wa_y[i];

        g_yyyy_0_yyzzzz_0[i] = 3.0 * g_yy_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzz_1[i] * wa_y[i];

        g_yyyy_0_yzzzzz_0[i] = 3.0 * g_yy_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzz_1[i] * wa_y[i];

        g_yyyy_0_zzzzzz_0[i] = 3.0 * g_yy_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 308-336 components of targeted buffer : GSI

    auto g_yyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 308);

    auto g_yyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 309);

    auto g_yyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 310);

    auto g_yyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 311);

    auto g_yyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 312);

    auto g_yyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 313);

    auto g_yyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 314);

    auto g_yyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 315);

    auto g_yyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 316);

    auto g_yyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 317);

    auto g_yyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 318);

    auto g_yyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 319);

    auto g_yyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 320);

    auto g_yyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 321);

    auto g_yyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 322);

    auto g_yyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 323);

    auto g_yyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 324);

    auto g_yyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 325);

    auto g_yyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 326);

    auto g_yyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 327);

    auto g_yyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 328);

    auto g_yyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 329);

    auto g_yyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 330);

    auto g_yyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 331);

    auto g_yyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 332);

    auto g_yyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 333);

    auto g_yyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 334);

    auto g_yyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 335);

    #pragma omp simd aligned(g_yyy_0_xxxxx_1, g_yyy_0_xxxxxx_1, g_yyy_0_xxxxxy_1, g_yyy_0_xxxxxz_1, g_yyy_0_xxxxy_1, g_yyy_0_xxxxyy_1, g_yyy_0_xxxxyz_1, g_yyy_0_xxxxz_1, g_yyy_0_xxxxzz_1, g_yyy_0_xxxyy_1, g_yyy_0_xxxyyy_1, g_yyy_0_xxxyyz_1, g_yyy_0_xxxyz_1, g_yyy_0_xxxyzz_1, g_yyy_0_xxxzz_1, g_yyy_0_xxxzzz_1, g_yyy_0_xxyyy_1, g_yyy_0_xxyyyy_1, g_yyy_0_xxyyyz_1, g_yyy_0_xxyyz_1, g_yyy_0_xxyyzz_1, g_yyy_0_xxyzz_1, g_yyy_0_xxyzzz_1, g_yyy_0_xxzzz_1, g_yyy_0_xxzzzz_1, g_yyy_0_xyyyy_1, g_yyy_0_xyyyyy_1, g_yyy_0_xyyyyz_1, g_yyy_0_xyyyz_1, g_yyy_0_xyyyzz_1, g_yyy_0_xyyzz_1, g_yyy_0_xyyzzz_1, g_yyy_0_xyzzz_1, g_yyy_0_xyzzzz_1, g_yyy_0_xzzzz_1, g_yyy_0_xzzzzz_1, g_yyy_0_yyyyy_1, g_yyy_0_yyyyyy_1, g_yyy_0_yyyyyz_1, g_yyy_0_yyyyz_1, g_yyy_0_yyyyzz_1, g_yyy_0_yyyzz_1, g_yyy_0_yyyzzz_1, g_yyy_0_yyzzz_1, g_yyy_0_yyzzzz_1, g_yyy_0_yzzzz_1, g_yyy_0_yzzzzz_1, g_yyy_0_zzzzz_1, g_yyy_0_zzzzzz_1, g_yyyz_0_xxxxxx_0, g_yyyz_0_xxxxxy_0, g_yyyz_0_xxxxxz_0, g_yyyz_0_xxxxyy_0, g_yyyz_0_xxxxyz_0, g_yyyz_0_xxxxzz_0, g_yyyz_0_xxxyyy_0, g_yyyz_0_xxxyyz_0, g_yyyz_0_xxxyzz_0, g_yyyz_0_xxxzzz_0, g_yyyz_0_xxyyyy_0, g_yyyz_0_xxyyyz_0, g_yyyz_0_xxyyzz_0, g_yyyz_0_xxyzzz_0, g_yyyz_0_xxzzzz_0, g_yyyz_0_xyyyyy_0, g_yyyz_0_xyyyyz_0, g_yyyz_0_xyyyzz_0, g_yyyz_0_xyyzzz_0, g_yyyz_0_xyzzzz_0, g_yyyz_0_xzzzzz_0, g_yyyz_0_yyyyyy_0, g_yyyz_0_yyyyyz_0, g_yyyz_0_yyyyzz_0, g_yyyz_0_yyyzzz_0, g_yyyz_0_yyzzzz_0, g_yyyz_0_yzzzzz_0, g_yyyz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxxxxx_0[i] = g_yyy_0_xxxxxx_1[i] * wa_z[i];

        g_yyyz_0_xxxxxy_0[i] = g_yyy_0_xxxxxy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxz_0[i] = g_yyy_0_xxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyy_0[i] = g_yyy_0_xxxxyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxyz_0[i] = g_yyy_0_xxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxzz_0[i] = 2.0 * g_yyy_0_xxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyy_0[i] = g_yyy_0_xxxyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxyyz_0[i] = g_yyy_0_xxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxyzz_0[i] = 2.0 * g_yyy_0_xxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxzzz_0[i] = 3.0 * g_yyy_0_xxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyy_0[i] = g_yyy_0_xxyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxyyyz_0[i] = g_yyy_0_xxyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxyyzz_0[i] = 2.0 * g_yyy_0_xxyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxyzzz_0[i] = 3.0 * g_yyy_0_xxyzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxzzzz_0[i] = 4.0 * g_yyy_0_xxzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyy_0[i] = g_yyy_0_xyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xyyyyz_0[i] = g_yyy_0_xyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xyyyzz_0[i] = 2.0 * g_yyy_0_xyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xyyzzz_0[i] = 3.0 * g_yyy_0_xyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xyzzzz_0[i] = 4.0 * g_yyy_0_xyzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xzzzzz_0[i] = 5.0 * g_yyy_0_xzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyy_0[i] = g_yyy_0_yyyyyy_1[i] * wa_z[i];

        g_yyyz_0_yyyyyz_0[i] = g_yyy_0_yyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyz_1[i] * wa_z[i];

        g_yyyz_0_yyyyzz_0[i] = 2.0 * g_yyy_0_yyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyzz_1[i] * wa_z[i];

        g_yyyz_0_yyyzzz_0[i] = 3.0 * g_yyy_0_yyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzz_1[i] * wa_z[i];

        g_yyyz_0_yyzzzz_0[i] = 4.0 * g_yyy_0_yyzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzz_1[i] * wa_z[i];

        g_yyyz_0_yzzzzz_0[i] = 5.0 * g_yyy_0_yzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzz_1[i] * wa_z[i];

        g_yyyz_0_zzzzzz_0[i] = 6.0 * g_yyy_0_zzzzz_1[i] * fi_acd_0 + g_yyy_0_zzzzzz_1[i] * wa_z[i];
    }

    /// Set up 336-364 components of targeted buffer : GSI

    auto g_yyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 336);

    auto g_yyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 337);

    auto g_yyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 338);

    auto g_yyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 339);

    auto g_yyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 340);

    auto g_yyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 341);

    auto g_yyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 342);

    auto g_yyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 343);

    auto g_yyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 344);

    auto g_yyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 345);

    auto g_yyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 346);

    auto g_yyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 347);

    auto g_yyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 348);

    auto g_yyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 349);

    auto g_yyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 350);

    auto g_yyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 351);

    auto g_yyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 352);

    auto g_yyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 353);

    auto g_yyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 354);

    auto g_yyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 355);

    auto g_yyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 356);

    auto g_yyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 357);

    auto g_yyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 358);

    auto g_yyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 359);

    auto g_yyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 360);

    auto g_yyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 361);

    auto g_yyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 362);

    auto g_yyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 363);

    #pragma omp simd aligned(g_yy_0_xxxxxy_0, g_yy_0_xxxxxy_1, g_yy_0_xxxxyy_0, g_yy_0_xxxxyy_1, g_yy_0_xxxyyy_0, g_yy_0_xxxyyy_1, g_yy_0_xxyyyy_0, g_yy_0_xxyyyy_1, g_yy_0_xyyyyy_0, g_yy_0_xyyyyy_1, g_yy_0_yyyyyy_0, g_yy_0_yyyyyy_1, g_yyz_0_xxxxxy_1, g_yyz_0_xxxxyy_1, g_yyz_0_xxxyyy_1, g_yyz_0_xxyyyy_1, g_yyz_0_xyyyyy_1, g_yyz_0_yyyyyy_1, g_yyzz_0_xxxxxx_0, g_yyzz_0_xxxxxy_0, g_yyzz_0_xxxxxz_0, g_yyzz_0_xxxxyy_0, g_yyzz_0_xxxxyz_0, g_yyzz_0_xxxxzz_0, g_yyzz_0_xxxyyy_0, g_yyzz_0_xxxyyz_0, g_yyzz_0_xxxyzz_0, g_yyzz_0_xxxzzz_0, g_yyzz_0_xxyyyy_0, g_yyzz_0_xxyyyz_0, g_yyzz_0_xxyyzz_0, g_yyzz_0_xxyzzz_0, g_yyzz_0_xxzzzz_0, g_yyzz_0_xyyyyy_0, g_yyzz_0_xyyyyz_0, g_yyzz_0_xyyyzz_0, g_yyzz_0_xyyzzz_0, g_yyzz_0_xyzzzz_0, g_yyzz_0_xzzzzz_0, g_yyzz_0_yyyyyy_0, g_yyzz_0_yyyyyz_0, g_yyzz_0_yyyyzz_0, g_yyzz_0_yyyzzz_0, g_yyzz_0_yyzzzz_0, g_yyzz_0_yzzzzz_0, g_yyzz_0_zzzzzz_0, g_yzz_0_xxxxxx_1, g_yzz_0_xxxxxz_1, g_yzz_0_xxxxyz_1, g_yzz_0_xxxxz_1, g_yzz_0_xxxxzz_1, g_yzz_0_xxxyyz_1, g_yzz_0_xxxyz_1, g_yzz_0_xxxyzz_1, g_yzz_0_xxxzz_1, g_yzz_0_xxxzzz_1, g_yzz_0_xxyyyz_1, g_yzz_0_xxyyz_1, g_yzz_0_xxyyzz_1, g_yzz_0_xxyzz_1, g_yzz_0_xxyzzz_1, g_yzz_0_xxzzz_1, g_yzz_0_xxzzzz_1, g_yzz_0_xyyyyz_1, g_yzz_0_xyyyz_1, g_yzz_0_xyyyzz_1, g_yzz_0_xyyzz_1, g_yzz_0_xyyzzz_1, g_yzz_0_xyzzz_1, g_yzz_0_xyzzzz_1, g_yzz_0_xzzzz_1, g_yzz_0_xzzzzz_1, g_yzz_0_yyyyyz_1, g_yzz_0_yyyyz_1, g_yzz_0_yyyyzz_1, g_yzz_0_yyyzz_1, g_yzz_0_yyyzzz_1, g_yzz_0_yyzzz_1, g_yzz_0_yyzzzz_1, g_yzz_0_yzzzz_1, g_yzz_0_yzzzzz_1, g_yzz_0_zzzzz_1, g_yzz_0_zzzzzz_1, g_zz_0_xxxxxx_0, g_zz_0_xxxxxx_1, g_zz_0_xxxxxz_0, g_zz_0_xxxxxz_1, g_zz_0_xxxxyz_0, g_zz_0_xxxxyz_1, g_zz_0_xxxxzz_0, g_zz_0_xxxxzz_1, g_zz_0_xxxyyz_0, g_zz_0_xxxyyz_1, g_zz_0_xxxyzz_0, g_zz_0_xxxyzz_1, g_zz_0_xxxzzz_0, g_zz_0_xxxzzz_1, g_zz_0_xxyyyz_0, g_zz_0_xxyyyz_1, g_zz_0_xxyyzz_0, g_zz_0_xxyyzz_1, g_zz_0_xxyzzz_0, g_zz_0_xxyzzz_1, g_zz_0_xxzzzz_0, g_zz_0_xxzzzz_1, g_zz_0_xyyyyz_0, g_zz_0_xyyyyz_1, g_zz_0_xyyyzz_0, g_zz_0_xyyyzz_1, g_zz_0_xyyzzz_0, g_zz_0_xyyzzz_1, g_zz_0_xyzzzz_0, g_zz_0_xyzzzz_1, g_zz_0_xzzzzz_0, g_zz_0_xzzzzz_1, g_zz_0_yyyyyz_0, g_zz_0_yyyyyz_1, g_zz_0_yyyyzz_0, g_zz_0_yyyyzz_1, g_zz_0_yyyzzz_0, g_zz_0_yyyzzz_1, g_zz_0_yyzzzz_0, g_zz_0_yyzzzz_1, g_zz_0_yzzzzz_0, g_zz_0_yzzzzz_1, g_zz_0_zzzzzz_0, g_zz_0_zzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxxxxx_0[i] = g_zz_0_xxxxxx_0[i] * fbe_0 - g_zz_0_xxxxxx_1[i] * fz_be_0 + g_yzz_0_xxxxxx_1[i] * wa_y[i];

        g_yyzz_0_xxxxxy_0[i] = g_yy_0_xxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxy_1[i] * fz_be_0 + g_yyz_0_xxxxxy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxz_0[i] = g_zz_0_xxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxz_1[i] * fz_be_0 + g_yzz_0_xxxxxz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyy_0[i] = g_yy_0_xxxxyy_0[i] * fbe_0 - g_yy_0_xxxxyy_1[i] * fz_be_0 + g_yyz_0_xxxxyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxyz_0[i] = g_zz_0_xxxxyz_0[i] * fbe_0 - g_zz_0_xxxxyz_1[i] * fz_be_0 + g_yzz_0_xxxxz_1[i] * fi_acd_0 + g_yzz_0_xxxxyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxzz_0[i] = g_zz_0_xxxxzz_0[i] * fbe_0 - g_zz_0_xxxxzz_1[i] * fz_be_0 + g_yzz_0_xxxxzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyy_0[i] = g_yy_0_xxxyyy_0[i] * fbe_0 - g_yy_0_xxxyyy_1[i] * fz_be_0 + g_yyz_0_xxxyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxyyz_0[i] = g_zz_0_xxxyyz_0[i] * fbe_0 - g_zz_0_xxxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxyzz_0[i] = g_zz_0_xxxyzz_0[i] * fbe_0 - g_zz_0_xxxyzz_1[i] * fz_be_0 + g_yzz_0_xxxzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxzzz_0[i] = g_zz_0_xxxzzz_0[i] * fbe_0 - g_zz_0_xxxzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyy_0[i] = g_yy_0_xxyyyy_0[i] * fbe_0 - g_yy_0_xxyyyy_1[i] * fz_be_0 + g_yyz_0_xxyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxyyyz_0[i] = g_zz_0_xxyyyz_0[i] * fbe_0 - g_zz_0_xxyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxyyzz_0[i] = g_zz_0_xxyyzz_0[i] * fbe_0 - g_zz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxyzzz_0[i] = g_zz_0_xxyzzz_0[i] * fbe_0 - g_zz_0_xxyzzz_1[i] * fz_be_0 + g_yzz_0_xxzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxzzzz_0[i] = g_zz_0_xxzzzz_0[i] * fbe_0 - g_zz_0_xxzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyy_0[i] = g_yy_0_xyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyy_1[i] * fz_be_0 + g_yyz_0_xyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xyyyyz_0[i] = g_zz_0_xyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xyyyzz_0[i] = g_zz_0_xyyyzz_0[i] * fbe_0 - g_zz_0_xyyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xyyzzz_0[i] = g_zz_0_xyyzzz_0[i] * fbe_0 - g_zz_0_xyyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xyzzzz_0[i] = g_zz_0_xyzzzz_0[i] * fbe_0 - g_zz_0_xyzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xzzzzz_0[i] = g_zz_0_xzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyy_0[i] = g_yy_0_yyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyy_1[i] * fz_be_0 + g_yyz_0_yyyyyy_1[i] * wa_z[i];

        g_yyzz_0_yyyyyz_0[i] = g_zz_0_yyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyz_1[i] * fz_be_0 + 5.0 * g_yzz_0_yyyyz_1[i] * fi_acd_0 + g_yzz_0_yyyyyz_1[i] * wa_y[i];

        g_yyzz_0_yyyyzz_0[i] = g_zz_0_yyyyzz_0[i] * fbe_0 - g_zz_0_yyyyzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_yyyzz_1[i] * fi_acd_0 + g_yzz_0_yyyyzz_1[i] * wa_y[i];

        g_yyzz_0_yyyzzz_0[i] = g_zz_0_yyyzzz_0[i] * fbe_0 - g_zz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_yyzzz_1[i] * fi_acd_0 + g_yzz_0_yyyzzz_1[i] * wa_y[i];

        g_yyzz_0_yyzzzz_0[i] = g_zz_0_yyzzzz_0[i] * fbe_0 - g_zz_0_yyzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yzzzz_1[i] * fi_acd_0 + g_yzz_0_yyzzzz_1[i] * wa_y[i];

        g_yyzz_0_yzzzzz_0[i] = g_zz_0_yzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzz_1[i] * fi_acd_0 + g_yzz_0_yzzzzz_1[i] * wa_y[i];

        g_yyzz_0_zzzzzz_0[i] = g_zz_0_zzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 364-392 components of targeted buffer : GSI

    auto g_yzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 364);

    auto g_yzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 365);

    auto g_yzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 366);

    auto g_yzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 367);

    auto g_yzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 368);

    auto g_yzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 369);

    auto g_yzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 370);

    auto g_yzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 371);

    auto g_yzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 372);

    auto g_yzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 373);

    auto g_yzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 374);

    auto g_yzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 375);

    auto g_yzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 376);

    auto g_yzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 377);

    auto g_yzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 378);

    auto g_yzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 379);

    auto g_yzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 380);

    auto g_yzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 381);

    auto g_yzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 382);

    auto g_yzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 383);

    auto g_yzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 384);

    auto g_yzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 385);

    auto g_yzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 386);

    auto g_yzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 387);

    auto g_yzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 388);

    auto g_yzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 389);

    auto g_yzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 390);

    auto g_yzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 391);

    #pragma omp simd aligned(g_yzzz_0_xxxxxx_0, g_yzzz_0_xxxxxy_0, g_yzzz_0_xxxxxz_0, g_yzzz_0_xxxxyy_0, g_yzzz_0_xxxxyz_0, g_yzzz_0_xxxxzz_0, g_yzzz_0_xxxyyy_0, g_yzzz_0_xxxyyz_0, g_yzzz_0_xxxyzz_0, g_yzzz_0_xxxzzz_0, g_yzzz_0_xxyyyy_0, g_yzzz_0_xxyyyz_0, g_yzzz_0_xxyyzz_0, g_yzzz_0_xxyzzz_0, g_yzzz_0_xxzzzz_0, g_yzzz_0_xyyyyy_0, g_yzzz_0_xyyyyz_0, g_yzzz_0_xyyyzz_0, g_yzzz_0_xyyzzz_0, g_yzzz_0_xyzzzz_0, g_yzzz_0_xzzzzz_0, g_yzzz_0_yyyyyy_0, g_yzzz_0_yyyyyz_0, g_yzzz_0_yyyyzz_0, g_yzzz_0_yyyzzz_0, g_yzzz_0_yyzzzz_0, g_yzzz_0_yzzzzz_0, g_yzzz_0_zzzzzz_0, g_zzz_0_xxxxx_1, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxy_1, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxz_1, g_zzz_0_xxxxzz_1, g_zzz_0_xxxyy_1, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyz_1, g_zzz_0_xxxyzz_1, g_zzz_0_xxxzz_1, g_zzz_0_xxxzzz_1, g_zzz_0_xxyyy_1, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyz_1, g_zzz_0_xxyyzz_1, g_zzz_0_xxyzz_1, g_zzz_0_xxyzzz_1, g_zzz_0_xxzzz_1, g_zzz_0_xxzzzz_1, g_zzz_0_xyyyy_1, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyz_1, g_zzz_0_xyyyzz_1, g_zzz_0_xyyzz_1, g_zzz_0_xyyzzz_1, g_zzz_0_xyzzz_1, g_zzz_0_xyzzzz_1, g_zzz_0_xzzzz_1, g_zzz_0_xzzzzz_1, g_zzz_0_yyyyy_1, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyz_1, g_zzz_0_yyyyzz_1, g_zzz_0_yyyzz_1, g_zzz_0_yyyzzz_1, g_zzz_0_yyzzz_1, g_zzz_0_yyzzzz_1, g_zzz_0_yzzzz_1, g_zzz_0_yzzzzz_1, g_zzz_0_zzzzz_1, g_zzz_0_zzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxxxxx_0[i] = g_zzz_0_xxxxxx_1[i] * wa_y[i];

        g_yzzz_0_xxxxxy_0[i] = g_zzz_0_xxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxz_0[i] = g_zzz_0_xxxxxz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyy_0[i] = 2.0 * g_zzz_0_xxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxyz_0[i] = g_zzz_0_xxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxzz_0[i] = g_zzz_0_xxxxzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyy_0[i] = 3.0 * g_zzz_0_xxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxyyz_0[i] = 2.0 * g_zzz_0_xxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxyzz_0[i] = g_zzz_0_xxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxzzz_0[i] = g_zzz_0_xxxzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyy_0[i] = 4.0 * g_zzz_0_xxyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxyyyz_0[i] = 3.0 * g_zzz_0_xxyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxyyzz_0[i] = 2.0 * g_zzz_0_xxyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxyzzz_0[i] = g_zzz_0_xxzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxzzzz_0[i] = g_zzz_0_xxzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyy_0[i] = 5.0 * g_zzz_0_xyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xyyyyz_0[i] = 4.0 * g_zzz_0_xyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xyyyzz_0[i] = 3.0 * g_zzz_0_xyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xyyzzz_0[i] = 2.0 * g_zzz_0_xyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xyzzzz_0[i] = g_zzz_0_xzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xzzzzz_0[i] = g_zzz_0_xzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyy_0[i] = 6.0 * g_zzz_0_yyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyy_1[i] * wa_y[i];

        g_yzzz_0_yyyyyz_0[i] = 5.0 * g_zzz_0_yyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyz_1[i] * wa_y[i];

        g_yzzz_0_yyyyzz_0[i] = 4.0 * g_zzz_0_yyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzz_1[i] * wa_y[i];

        g_yzzz_0_yyyzzz_0[i] = 3.0 * g_zzz_0_yyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzz_1[i] * wa_y[i];

        g_yzzz_0_yyzzzz_0[i] = 2.0 * g_zzz_0_yzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzz_1[i] * wa_y[i];

        g_yzzz_0_yzzzzz_0[i] = g_zzz_0_zzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzz_1[i] * wa_y[i];

        g_yzzz_0_zzzzzz_0[i] = g_zzz_0_zzzzzz_1[i] * wa_y[i];
    }

    /// Set up 392-420 components of targeted buffer : GSI

    auto g_zzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_gsi + 392);

    auto g_zzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_gsi + 393);

    auto g_zzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_gsi + 394);

    auto g_zzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_gsi + 395);

    auto g_zzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_gsi + 396);

    auto g_zzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_gsi + 397);

    auto g_zzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_gsi + 398);

    auto g_zzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_gsi + 399);

    auto g_zzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_gsi + 400);

    auto g_zzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_gsi + 401);

    auto g_zzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_gsi + 402);

    auto g_zzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_gsi + 403);

    auto g_zzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_gsi + 404);

    auto g_zzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_gsi + 405);

    auto g_zzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_gsi + 406);

    auto g_zzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 407);

    auto g_zzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 408);

    auto g_zzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 409);

    auto g_zzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 410);

    auto g_zzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 411);

    auto g_zzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 412);

    auto g_zzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_gsi + 413);

    auto g_zzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_gsi + 414);

    auto g_zzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_gsi + 415);

    auto g_zzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_gsi + 416);

    auto g_zzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_gsi + 417);

    auto g_zzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 418);

    auto g_zzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_gsi + 419);

    #pragma omp simd aligned(g_zz_0_xxxxxx_0, g_zz_0_xxxxxx_1, g_zz_0_xxxxxy_0, g_zz_0_xxxxxy_1, g_zz_0_xxxxxz_0, g_zz_0_xxxxxz_1, g_zz_0_xxxxyy_0, g_zz_0_xxxxyy_1, g_zz_0_xxxxyz_0, g_zz_0_xxxxyz_1, g_zz_0_xxxxzz_0, g_zz_0_xxxxzz_1, g_zz_0_xxxyyy_0, g_zz_0_xxxyyy_1, g_zz_0_xxxyyz_0, g_zz_0_xxxyyz_1, g_zz_0_xxxyzz_0, g_zz_0_xxxyzz_1, g_zz_0_xxxzzz_0, g_zz_0_xxxzzz_1, g_zz_0_xxyyyy_0, g_zz_0_xxyyyy_1, g_zz_0_xxyyyz_0, g_zz_0_xxyyyz_1, g_zz_0_xxyyzz_0, g_zz_0_xxyyzz_1, g_zz_0_xxyzzz_0, g_zz_0_xxyzzz_1, g_zz_0_xxzzzz_0, g_zz_0_xxzzzz_1, g_zz_0_xyyyyy_0, g_zz_0_xyyyyy_1, g_zz_0_xyyyyz_0, g_zz_0_xyyyyz_1, g_zz_0_xyyyzz_0, g_zz_0_xyyyzz_1, g_zz_0_xyyzzz_0, g_zz_0_xyyzzz_1, g_zz_0_xyzzzz_0, g_zz_0_xyzzzz_1, g_zz_0_xzzzzz_0, g_zz_0_xzzzzz_1, g_zz_0_yyyyyy_0, g_zz_0_yyyyyy_1, g_zz_0_yyyyyz_0, g_zz_0_yyyyyz_1, g_zz_0_yyyyzz_0, g_zz_0_yyyyzz_1, g_zz_0_yyyzzz_0, g_zz_0_yyyzzz_1, g_zz_0_yyzzzz_0, g_zz_0_yyzzzz_1, g_zz_0_yzzzzz_0, g_zz_0_yzzzzz_1, g_zz_0_zzzzzz_0, g_zz_0_zzzzzz_1, g_zzz_0_xxxxx_1, g_zzz_0_xxxxxx_1, g_zzz_0_xxxxxy_1, g_zzz_0_xxxxxz_1, g_zzz_0_xxxxy_1, g_zzz_0_xxxxyy_1, g_zzz_0_xxxxyz_1, g_zzz_0_xxxxz_1, g_zzz_0_xxxxzz_1, g_zzz_0_xxxyy_1, g_zzz_0_xxxyyy_1, g_zzz_0_xxxyyz_1, g_zzz_0_xxxyz_1, g_zzz_0_xxxyzz_1, g_zzz_0_xxxzz_1, g_zzz_0_xxxzzz_1, g_zzz_0_xxyyy_1, g_zzz_0_xxyyyy_1, g_zzz_0_xxyyyz_1, g_zzz_0_xxyyz_1, g_zzz_0_xxyyzz_1, g_zzz_0_xxyzz_1, g_zzz_0_xxyzzz_1, g_zzz_0_xxzzz_1, g_zzz_0_xxzzzz_1, g_zzz_0_xyyyy_1, g_zzz_0_xyyyyy_1, g_zzz_0_xyyyyz_1, g_zzz_0_xyyyz_1, g_zzz_0_xyyyzz_1, g_zzz_0_xyyzz_1, g_zzz_0_xyyzzz_1, g_zzz_0_xyzzz_1, g_zzz_0_xyzzzz_1, g_zzz_0_xzzzz_1, g_zzz_0_xzzzzz_1, g_zzz_0_yyyyy_1, g_zzz_0_yyyyyy_1, g_zzz_0_yyyyyz_1, g_zzz_0_yyyyz_1, g_zzz_0_yyyyzz_1, g_zzz_0_yyyzz_1, g_zzz_0_yyyzzz_1, g_zzz_0_yyzzz_1, g_zzz_0_yyzzzz_1, g_zzz_0_yzzzz_1, g_zzz_0_yzzzzz_1, g_zzz_0_zzzzz_1, g_zzz_0_zzzzzz_1, g_zzzz_0_xxxxxx_0, g_zzzz_0_xxxxxy_0, g_zzzz_0_xxxxxz_0, g_zzzz_0_xxxxyy_0, g_zzzz_0_xxxxyz_0, g_zzzz_0_xxxxzz_0, g_zzzz_0_xxxyyy_0, g_zzzz_0_xxxyyz_0, g_zzzz_0_xxxyzz_0, g_zzzz_0_xxxzzz_0, g_zzzz_0_xxyyyy_0, g_zzzz_0_xxyyyz_0, g_zzzz_0_xxyyzz_0, g_zzzz_0_xxyzzz_0, g_zzzz_0_xxzzzz_0, g_zzzz_0_xyyyyy_0, g_zzzz_0_xyyyyz_0, g_zzzz_0_xyyyzz_0, g_zzzz_0_xyyzzz_0, g_zzzz_0_xyzzzz_0, g_zzzz_0_xzzzzz_0, g_zzzz_0_yyyyyy_0, g_zzzz_0_yyyyyz_0, g_zzzz_0_yyyyzz_0, g_zzzz_0_yyyzzz_0, g_zzzz_0_yyzzzz_0, g_zzzz_0_yzzzzz_0, g_zzzz_0_zzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxxxxx_0[i] = 3.0 * g_zz_0_xxxxxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxx_1[i] * fz_be_0 + g_zzz_0_xxxxxx_1[i] * wa_z[i];

        g_zzzz_0_xxxxxy_0[i] = 3.0 * g_zz_0_xxxxxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxy_1[i] * fz_be_0 + g_zzz_0_xxxxxy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxz_0[i] = 3.0 * g_zz_0_xxxxxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxz_1[i] * fz_be_0 + g_zzz_0_xxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyy_0[i] = 3.0 * g_zz_0_xxxxyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyy_1[i] * fz_be_0 + g_zzz_0_xxxxyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxyz_0[i] = 3.0 * g_zz_0_xxxxyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyz_1[i] * fz_be_0 + g_zzz_0_xxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxzz_0[i] = 3.0 * g_zz_0_xxxxzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyy_0[i] = 3.0 * g_zz_0_xxxyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyy_1[i] * fz_be_0 + g_zzz_0_xxxyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxyyz_0[i] = 3.0 * g_zz_0_xxxyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyz_1[i] * fz_be_0 + g_zzz_0_xxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxyzz_0[i] = 3.0 * g_zz_0_xxxyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxzzz_0[i] = 3.0 * g_zz_0_xxxzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyy_0[i] = 3.0 * g_zz_0_xxyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyy_1[i] * fz_be_0 + g_zzz_0_xxyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxyyyz_0[i] = 3.0 * g_zz_0_xxyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyz_1[i] * fz_be_0 + g_zzz_0_xxyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxyyzz_0[i] = 3.0 * g_zz_0_xxyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxyzzz_0[i] = 3.0 * g_zz_0_xxyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxyzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxzzzz_0[i] = 3.0 * g_zz_0_xxzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyy_0[i] = 3.0 * g_zz_0_xyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyy_1[i] * fz_be_0 + g_zzz_0_xyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xyyyyz_0[i] = 3.0 * g_zz_0_xyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyz_1[i] * fz_be_0 + g_zzz_0_xyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xyyyzz_0[i] = 3.0 * g_zz_0_xyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xyyzzz_0[i] = 3.0 * g_zz_0_xyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xyzzzz_0[i] = 3.0 * g_zz_0_xyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xyzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xzzzzz_0[i] = 3.0 * g_zz_0_xzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyy_0[i] = 3.0 * g_zz_0_yyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyy_1[i] * fz_be_0 + g_zzz_0_yyyyyy_1[i] * wa_z[i];

        g_zzzz_0_yyyyyz_0[i] = 3.0 * g_zz_0_yyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyz_1[i] * fz_be_0 + g_zzz_0_yyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyz_1[i] * wa_z[i];

        g_zzzz_0_yyyyzz_0[i] = 3.0 * g_zz_0_yyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyzz_1[i] * wa_z[i];

        g_zzzz_0_yyyzzz_0[i] = 3.0 * g_zz_0_yyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_yyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzz_1[i] * wa_z[i];

        g_zzzz_0_yyzzzz_0[i] = 3.0 * g_zz_0_yyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_yyzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzz_1[i] * wa_z[i];

        g_zzzz_0_yzzzzz_0[i] = 3.0 * g_zz_0_yzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_yzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzz_1[i] * wa_z[i];

        g_zzzz_0_zzzzzz_0[i] = 3.0 * g_zz_0_zzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_zzzzz_1[i] * fi_acd_0 + g_zzz_0_zzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

