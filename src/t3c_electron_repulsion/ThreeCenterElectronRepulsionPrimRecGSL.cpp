#include "ThreeCenterElectronRepulsionPrimRecGSL.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsl(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsl,
                                 size_t idx_eri_0_dsl,
                                 size_t idx_eri_1_dsl,
                                 size_t idx_eri_1_fsk,
                                 size_t idx_eri_1_fsl,
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

    /// Set up components of auxilary buffer : DSL

    auto g_xx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl);

    auto g_xx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 1);

    auto g_xx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 2);

    auto g_xx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 3);

    auto g_xx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 4);

    auto g_xx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 5);

    auto g_xx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 6);

    auto g_xx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 7);

    auto g_xx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 8);

    auto g_xx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 9);

    auto g_xx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 10);

    auto g_xx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 11);

    auto g_xx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 12);

    auto g_xx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 13);

    auto g_xx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 14);

    auto g_xx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 15);

    auto g_xx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 16);

    auto g_xx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 17);

    auto g_xx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 18);

    auto g_xx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 19);

    auto g_xx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 20);

    auto g_xx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 21);

    auto g_xx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 22);

    auto g_xx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 23);

    auto g_xx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 24);

    auto g_xx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 25);

    auto g_xx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 26);

    auto g_xx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 27);

    auto g_xx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 28);

    auto g_xx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 29);

    auto g_xx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 30);

    auto g_xx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 31);

    auto g_xx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 32);

    auto g_xx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 33);

    auto g_xx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 34);

    auto g_xx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 35);

    auto g_xx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 36);

    auto g_xx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 37);

    auto g_xx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 38);

    auto g_xx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 39);

    auto g_xx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 40);

    auto g_xx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 41);

    auto g_xx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 42);

    auto g_xx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 43);

    auto g_xx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 44);

    auto g_yy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 135);

    auto g_yy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 136);

    auto g_yy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 137);

    auto g_yy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 138);

    auto g_yy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 139);

    auto g_yy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 140);

    auto g_yy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 141);

    auto g_yy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 142);

    auto g_yy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 143);

    auto g_yy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 144);

    auto g_yy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 145);

    auto g_yy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 146);

    auto g_yy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 147);

    auto g_yy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 148);

    auto g_yy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 149);

    auto g_yy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 150);

    auto g_yy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 151);

    auto g_yy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 152);

    auto g_yy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 153);

    auto g_yy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 154);

    auto g_yy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 155);

    auto g_yy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 156);

    auto g_yy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 157);

    auto g_yy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 158);

    auto g_yy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 159);

    auto g_yy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 160);

    auto g_yy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 161);

    auto g_yy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 162);

    auto g_yy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 163);

    auto g_yy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 164);

    auto g_yy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 165);

    auto g_yy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 166);

    auto g_yy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 167);

    auto g_yy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 168);

    auto g_yy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 169);

    auto g_yy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 170);

    auto g_yy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 171);

    auto g_yy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 172);

    auto g_yy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 173);

    auto g_yy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 174);

    auto g_yy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 175);

    auto g_yy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 176);

    auto g_yy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 177);

    auto g_yy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 178);

    auto g_yy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 179);

    auto g_zz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_dsl + 225);

    auto g_zz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_dsl + 226);

    auto g_zz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_dsl + 227);

    auto g_zz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_dsl + 228);

    auto g_zz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_dsl + 229);

    auto g_zz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_dsl + 230);

    auto g_zz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_dsl + 231);

    auto g_zz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_dsl + 232);

    auto g_zz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_dsl + 233);

    auto g_zz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_dsl + 234);

    auto g_zz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_dsl + 235);

    auto g_zz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_dsl + 236);

    auto g_zz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_dsl + 237);

    auto g_zz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_dsl + 238);

    auto g_zz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_dsl + 239);

    auto g_zz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 240);

    auto g_zz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 241);

    auto g_zz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 242);

    auto g_zz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 243);

    auto g_zz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 244);

    auto g_zz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 245);

    auto g_zz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 246);

    auto g_zz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 247);

    auto g_zz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 248);

    auto g_zz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 249);

    auto g_zz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 250);

    auto g_zz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 251);

    auto g_zz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 252);

    auto g_zz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 253);

    auto g_zz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 254);

    auto g_zz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 255);

    auto g_zz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 256);

    auto g_zz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 257);

    auto g_zz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 258);

    auto g_zz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 259);

    auto g_zz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 260);

    auto g_zz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_dsl + 261);

    auto g_zz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_dsl + 262);

    auto g_zz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_dsl + 263);

    auto g_zz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_dsl + 264);

    auto g_zz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_dsl + 265);

    auto g_zz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 266);

    auto g_zz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 267);

    auto g_zz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 268);

    auto g_zz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_dsl + 269);

    /// Set up components of auxilary buffer : DSL

    auto g_xx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl);

    auto g_xx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 1);

    auto g_xx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 2);

    auto g_xx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 3);

    auto g_xx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 4);

    auto g_xx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 5);

    auto g_xx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 6);

    auto g_xx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 7);

    auto g_xx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 8);

    auto g_xx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 9);

    auto g_xx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 10);

    auto g_xx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 11);

    auto g_xx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 12);

    auto g_xx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 13);

    auto g_xx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 14);

    auto g_xx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 15);

    auto g_xx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 16);

    auto g_xx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 17);

    auto g_xx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 18);

    auto g_xx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 19);

    auto g_xx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 20);

    auto g_xx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 21);

    auto g_xx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 22);

    auto g_xx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 23);

    auto g_xx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 24);

    auto g_xx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 25);

    auto g_xx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 26);

    auto g_xx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 27);

    auto g_xx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 28);

    auto g_xx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 29);

    auto g_xx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 30);

    auto g_xx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 31);

    auto g_xx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 32);

    auto g_xx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 33);

    auto g_xx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 34);

    auto g_xx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 35);

    auto g_xx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 36);

    auto g_xx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 37);

    auto g_xx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 38);

    auto g_xx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 39);

    auto g_xx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 40);

    auto g_xx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 41);

    auto g_xx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 42);

    auto g_xx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 43);

    auto g_xx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 44);

    auto g_yy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl + 135);

    auto g_yy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 136);

    auto g_yy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 137);

    auto g_yy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 138);

    auto g_yy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 139);

    auto g_yy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 140);

    auto g_yy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 141);

    auto g_yy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 142);

    auto g_yy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 143);

    auto g_yy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 144);

    auto g_yy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 145);

    auto g_yy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 146);

    auto g_yy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 147);

    auto g_yy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 148);

    auto g_yy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 149);

    auto g_yy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 150);

    auto g_yy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 151);

    auto g_yy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 152);

    auto g_yy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 153);

    auto g_yy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 154);

    auto g_yy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 155);

    auto g_yy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 156);

    auto g_yy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 157);

    auto g_yy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 158);

    auto g_yy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 159);

    auto g_yy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 160);

    auto g_yy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 161);

    auto g_yy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 162);

    auto g_yy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 163);

    auto g_yy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 164);

    auto g_yy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 165);

    auto g_yy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 166);

    auto g_yy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 167);

    auto g_yy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 168);

    auto g_yy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 169);

    auto g_yy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 170);

    auto g_yy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 171);

    auto g_yy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 172);

    auto g_yy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 173);

    auto g_yy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 174);

    auto g_yy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 175);

    auto g_yy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 176);

    auto g_yy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 177);

    auto g_yy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 178);

    auto g_yy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 179);

    auto g_zz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_dsl + 225);

    auto g_zz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_dsl + 226);

    auto g_zz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_dsl + 227);

    auto g_zz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_dsl + 228);

    auto g_zz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_dsl + 229);

    auto g_zz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_dsl + 230);

    auto g_zz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_dsl + 231);

    auto g_zz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_dsl + 232);

    auto g_zz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_dsl + 233);

    auto g_zz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_dsl + 234);

    auto g_zz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_dsl + 235);

    auto g_zz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_dsl + 236);

    auto g_zz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_dsl + 237);

    auto g_zz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_dsl + 238);

    auto g_zz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_dsl + 239);

    auto g_zz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 240);

    auto g_zz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 241);

    auto g_zz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 242);

    auto g_zz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 243);

    auto g_zz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 244);

    auto g_zz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 245);

    auto g_zz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 246);

    auto g_zz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 247);

    auto g_zz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 248);

    auto g_zz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 249);

    auto g_zz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 250);

    auto g_zz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 251);

    auto g_zz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 252);

    auto g_zz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 253);

    auto g_zz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 254);

    auto g_zz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 255);

    auto g_zz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 256);

    auto g_zz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 257);

    auto g_zz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 258);

    auto g_zz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 259);

    auto g_zz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 260);

    auto g_zz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_dsl + 261);

    auto g_zz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_dsl + 262);

    auto g_zz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_dsl + 263);

    auto g_zz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_dsl + 264);

    auto g_zz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_dsl + 265);

    auto g_zz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 266);

    auto g_zz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 267);

    auto g_zz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 268);

    auto g_zz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_dsl + 269);

    /// Set up components of auxilary buffer : FSK

    auto g_xxx_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk);

    auto g_xxx_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 1);

    auto g_xxx_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 2);

    auto g_xxx_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 3);

    auto g_xxx_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 4);

    auto g_xxx_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 5);

    auto g_xxx_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 6);

    auto g_xxx_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 7);

    auto g_xxx_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 8);

    auto g_xxx_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 9);

    auto g_xxx_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 10);

    auto g_xxx_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 11);

    auto g_xxx_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 12);

    auto g_xxx_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 13);

    auto g_xxx_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 14);

    auto g_xxx_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 15);

    auto g_xxx_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 16);

    auto g_xxx_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 17);

    auto g_xxx_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 18);

    auto g_xxx_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 19);

    auto g_xxx_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 20);

    auto g_xxx_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 21);

    auto g_xxx_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 22);

    auto g_xxx_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 23);

    auto g_xxx_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 24);

    auto g_xxx_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 25);

    auto g_xxx_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 26);

    auto g_xxx_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 27);

    auto g_xxx_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 28);

    auto g_xxx_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 29);

    auto g_xxx_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 30);

    auto g_xxx_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 31);

    auto g_xxx_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 32);

    auto g_xxx_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 33);

    auto g_xxx_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 34);

    auto g_xxx_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 35);

    auto g_xxz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 74);

    auto g_xxz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 76);

    auto g_xxz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 77);

    auto g_xxz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 79);

    auto g_xxz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 80);

    auto g_xxz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 81);

    auto g_xxz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 83);

    auto g_xxz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 84);

    auto g_xxz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 85);

    auto g_xxz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 86);

    auto g_xxz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 88);

    auto g_xxz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 89);

    auto g_xxz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 90);

    auto g_xxz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 91);

    auto g_xxz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 92);

    auto g_xxz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 94);

    auto g_xxz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 95);

    auto g_xxz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 96);

    auto g_xxz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 97);

    auto g_xxz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 98);

    auto g_xxz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 99);

    auto g_xxz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 101);

    auto g_xxz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 102);

    auto g_xxz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 103);

    auto g_xxz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 104);

    auto g_xxz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 105);

    auto g_xxz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 106);

    auto g_xxz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 107);

    auto g_xyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 109);

    auto g_xyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 111);

    auto g_xyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 112);

    auto g_xyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 114);

    auto g_xyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 115);

    auto g_xyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 116);

    auto g_xyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 118);

    auto g_xyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 119);

    auto g_xyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 120);

    auto g_xyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 121);

    auto g_xyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 123);

    auto g_xyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 124);

    auto g_xyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 125);

    auto g_xyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 126);

    auto g_xyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 127);

    auto g_xyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 129);

    auto g_xyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 130);

    auto g_xyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 131);

    auto g_xyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 132);

    auto g_xyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 133);

    auto g_xyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 134);

    auto g_xyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 136);

    auto g_xyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 137);

    auto g_xyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 138);

    auto g_xyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 139);

    auto g_xyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 140);

    auto g_xyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 141);

    auto g_xyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 142);

    auto g_xzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 182);

    auto g_xzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 184);

    auto g_xzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 185);

    auto g_xzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 187);

    auto g_xzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 188);

    auto g_xzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 189);

    auto g_xzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 191);

    auto g_xzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 192);

    auto g_xzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 193);

    auto g_xzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 194);

    auto g_xzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 196);

    auto g_xzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 197);

    auto g_xzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 198);

    auto g_xzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 199);

    auto g_xzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 200);

    auto g_xzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 202);

    auto g_xzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 203);

    auto g_xzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 204);

    auto g_xzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 205);

    auto g_xzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 206);

    auto g_xzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 207);

    auto g_xzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 209);

    auto g_xzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 210);

    auto g_xzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 211);

    auto g_xzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 212);

    auto g_xzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 213);

    auto g_xzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 214);

    auto g_xzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 215);

    auto g_yyy_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 216);

    auto g_yyy_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 217);

    auto g_yyy_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 218);

    auto g_yyy_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 219);

    auto g_yyy_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 220);

    auto g_yyy_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 221);

    auto g_yyy_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 222);

    auto g_yyy_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 223);

    auto g_yyy_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 224);

    auto g_yyy_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 225);

    auto g_yyy_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 226);

    auto g_yyy_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 227);

    auto g_yyy_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 228);

    auto g_yyy_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 229);

    auto g_yyy_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 230);

    auto g_yyy_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 231);

    auto g_yyy_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 232);

    auto g_yyy_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 233);

    auto g_yyy_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 234);

    auto g_yyy_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 235);

    auto g_yyy_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 236);

    auto g_yyy_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 237);

    auto g_yyy_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 238);

    auto g_yyy_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 239);

    auto g_yyy_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 240);

    auto g_yyy_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 241);

    auto g_yyy_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 242);

    auto g_yyy_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 243);

    auto g_yyy_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 244);

    auto g_yyy_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 245);

    auto g_yyy_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 246);

    auto g_yyy_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 247);

    auto g_yyy_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 248);

    auto g_yyy_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 249);

    auto g_yyy_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 250);

    auto g_yyy_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 251);

    auto g_yyz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 254);

    auto g_yyz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 256);

    auto g_yyz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 257);

    auto g_yyz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 259);

    auto g_yyz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 260);

    auto g_yyz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 261);

    auto g_yyz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 263);

    auto g_yyz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 264);

    auto g_yyz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 265);

    auto g_yyz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 266);

    auto g_yyz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 268);

    auto g_yyz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 269);

    auto g_yyz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 270);

    auto g_yyz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 271);

    auto g_yyz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 272);

    auto g_yyz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 274);

    auto g_yyz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 275);

    auto g_yyz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 276);

    auto g_yyz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 277);

    auto g_yyz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 278);

    auto g_yyz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 279);

    auto g_yyz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 281);

    auto g_yyz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 282);

    auto g_yyz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 283);

    auto g_yyz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 284);

    auto g_yyz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 285);

    auto g_yyz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 286);

    auto g_yyz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 287);

    auto g_yzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 289);

    auto g_yzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 290);

    auto g_yzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 291);

    auto g_yzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 292);

    auto g_yzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 293);

    auto g_yzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 294);

    auto g_yzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 295);

    auto g_yzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 296);

    auto g_yzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 297);

    auto g_yzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 298);

    auto g_yzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 299);

    auto g_yzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 300);

    auto g_yzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 301);

    auto g_yzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 302);

    auto g_yzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 303);

    auto g_yzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 304);

    auto g_yzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 305);

    auto g_yzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 306);

    auto g_yzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 307);

    auto g_yzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 308);

    auto g_yzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 309);

    auto g_yzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 310);

    auto g_yzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 311);

    auto g_yzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 312);

    auto g_yzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 313);

    auto g_yzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 314);

    auto g_yzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 315);

    auto g_yzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 316);

    auto g_yzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 317);

    auto g_yzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 318);

    auto g_yzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 319);

    auto g_yzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 320);

    auto g_yzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 321);

    auto g_yzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 322);

    auto g_yzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 323);

    auto g_zzz_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_fsk + 324);

    auto g_zzz_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_fsk + 325);

    auto g_zzz_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_fsk + 326);

    auto g_zzz_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_fsk + 327);

    auto g_zzz_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_fsk + 328);

    auto g_zzz_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_fsk + 329);

    auto g_zzz_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_fsk + 330);

    auto g_zzz_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_fsk + 331);

    auto g_zzz_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_fsk + 332);

    auto g_zzz_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_fsk + 333);

    auto g_zzz_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_fsk + 334);

    auto g_zzz_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_fsk + 335);

    auto g_zzz_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_fsk + 336);

    auto g_zzz_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_fsk + 337);

    auto g_zzz_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_fsk + 338);

    auto g_zzz_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 339);

    auto g_zzz_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 340);

    auto g_zzz_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 341);

    auto g_zzz_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 342);

    auto g_zzz_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 343);

    auto g_zzz_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 344);

    auto g_zzz_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 345);

    auto g_zzz_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 346);

    auto g_zzz_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 347);

    auto g_zzz_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 348);

    auto g_zzz_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 349);

    auto g_zzz_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 350);

    auto g_zzz_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 351);

    auto g_zzz_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_fsk + 352);

    auto g_zzz_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_fsk + 353);

    auto g_zzz_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_fsk + 354);

    auto g_zzz_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_fsk + 355);

    auto g_zzz_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_fsk + 356);

    auto g_zzz_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 357);

    auto g_zzz_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 358);

    auto g_zzz_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_fsk + 359);

    /// Set up components of auxilary buffer : FSL

    auto g_xxx_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl);

    auto g_xxx_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 1);

    auto g_xxx_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 2);

    auto g_xxx_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 3);

    auto g_xxx_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 4);

    auto g_xxx_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 5);

    auto g_xxx_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 6);

    auto g_xxx_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 7);

    auto g_xxx_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 8);

    auto g_xxx_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 9);

    auto g_xxx_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 10);

    auto g_xxx_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 11);

    auto g_xxx_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 12);

    auto g_xxx_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 13);

    auto g_xxx_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 14);

    auto g_xxx_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 15);

    auto g_xxx_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 16);

    auto g_xxx_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 17);

    auto g_xxx_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 18);

    auto g_xxx_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 19);

    auto g_xxx_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 20);

    auto g_xxx_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 21);

    auto g_xxx_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 22);

    auto g_xxx_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 23);

    auto g_xxx_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 24);

    auto g_xxx_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 25);

    auto g_xxx_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 26);

    auto g_xxx_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 27);

    auto g_xxx_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 28);

    auto g_xxx_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 29);

    auto g_xxx_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 30);

    auto g_xxx_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 31);

    auto g_xxx_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 32);

    auto g_xxx_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 33);

    auto g_xxx_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 34);

    auto g_xxx_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 35);

    auto g_xxx_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 36);

    auto g_xxx_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 37);

    auto g_xxx_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 38);

    auto g_xxx_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 39);

    auto g_xxx_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 40);

    auto g_xxx_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 41);

    auto g_xxx_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 42);

    auto g_xxx_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 43);

    auto g_xxx_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 44);

    auto g_xxy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 45);

    auto g_xxy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 46);

    auto g_xxy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 47);

    auto g_xxy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 48);

    auto g_xxy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 50);

    auto g_xxy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 51);

    auto g_xxy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 54);

    auto g_xxy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 55);

    auto g_xxy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 59);

    auto g_xxy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 60);

    auto g_xxy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 65);

    auto g_xxy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 66);

    auto g_xxy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 72);

    auto g_xxy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 73);

    auto g_xxy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 80);

    auto g_xxy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 81);

    auto g_xxz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 90);

    auto g_xxz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 91);

    auto g_xxz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 92);

    auto g_xxz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 93);

    auto g_xxz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 94);

    auto g_xxz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 95);

    auto g_xxz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 96);

    auto g_xxz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 97);

    auto g_xxz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 98);

    auto g_xxz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 99);

    auto g_xxz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 100);

    auto g_xxz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 101);

    auto g_xxz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 102);

    auto g_xxz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 103);

    auto g_xxz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 104);

    auto g_xxz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 105);

    auto g_xxz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 106);

    auto g_xxz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 107);

    auto g_xxz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 108);

    auto g_xxz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 109);

    auto g_xxz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 110);

    auto g_xxz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 111);

    auto g_xxz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 112);

    auto g_xxz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 113);

    auto g_xxz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 114);

    auto g_xxz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 115);

    auto g_xxz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 116);

    auto g_xxz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 117);

    auto g_xxz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 118);

    auto g_xxz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 119);

    auto g_xxz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 120);

    auto g_xxz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 121);

    auto g_xxz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 122);

    auto g_xxz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 123);

    auto g_xxz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 124);

    auto g_xxz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 125);

    auto g_xxz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 127);

    auto g_xxz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 128);

    auto g_xxz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 129);

    auto g_xxz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 130);

    auto g_xxz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 131);

    auto g_xxz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 132);

    auto g_xxz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 133);

    auto g_xxz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 134);

    auto g_xyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 135);

    auto g_xyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 136);

    auto g_xyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 138);

    auto g_xyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 139);

    auto g_xyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 141);

    auto g_xyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 142);

    auto g_xyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 143);

    auto g_xyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 145);

    auto g_xyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 146);

    auto g_xyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 147);

    auto g_xyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 148);

    auto g_xyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 150);

    auto g_xyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 151);

    auto g_xyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 152);

    auto g_xyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 153);

    auto g_xyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 154);

    auto g_xyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 156);

    auto g_xyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 157);

    auto g_xyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 158);

    auto g_xyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 159);

    auto g_xyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 160);

    auto g_xyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 161);

    auto g_xyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 163);

    auto g_xyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 164);

    auto g_xyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 165);

    auto g_xyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 166);

    auto g_xyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 167);

    auto g_xyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 168);

    auto g_xyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 169);

    auto g_xyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 171);

    auto g_xyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 172);

    auto g_xyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 173);

    auto g_xyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 174);

    auto g_xyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 175);

    auto g_xyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 176);

    auto g_xyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 177);

    auto g_xyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 178);

    auto g_xyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 179);

    auto g_xzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 225);

    auto g_xzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 227);

    auto g_xzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 229);

    auto g_xzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 230);

    auto g_xzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 232);

    auto g_xzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 233);

    auto g_xzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 234);

    auto g_xzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 236);

    auto g_xzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 237);

    auto g_xzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 238);

    auto g_xzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 239);

    auto g_xzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 241);

    auto g_xzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 242);

    auto g_xzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 243);

    auto g_xzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 244);

    auto g_xzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 245);

    auto g_xzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 247);

    auto g_xzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 248);

    auto g_xzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 249);

    auto g_xzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 250);

    auto g_xzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 251);

    auto g_xzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 252);

    auto g_xzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 254);

    auto g_xzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 255);

    auto g_xzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 256);

    auto g_xzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 257);

    auto g_xzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 258);

    auto g_xzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 259);

    auto g_xzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 260);

    auto g_xzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 261);

    auto g_xzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 262);

    auto g_xzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 263);

    auto g_xzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 264);

    auto g_xzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 265);

    auto g_xzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 266);

    auto g_xzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 267);

    auto g_xzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 268);

    auto g_xzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 269);

    auto g_yyy_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 270);

    auto g_yyy_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 271);

    auto g_yyy_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 272);

    auto g_yyy_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 273);

    auto g_yyy_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 274);

    auto g_yyy_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 275);

    auto g_yyy_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 276);

    auto g_yyy_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 277);

    auto g_yyy_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 278);

    auto g_yyy_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 279);

    auto g_yyy_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 280);

    auto g_yyy_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 281);

    auto g_yyy_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 282);

    auto g_yyy_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 283);

    auto g_yyy_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 284);

    auto g_yyy_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 285);

    auto g_yyy_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 286);

    auto g_yyy_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 287);

    auto g_yyy_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 288);

    auto g_yyy_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 289);

    auto g_yyy_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 290);

    auto g_yyy_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 291);

    auto g_yyy_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 292);

    auto g_yyy_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 293);

    auto g_yyy_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 294);

    auto g_yyy_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 295);

    auto g_yyy_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 296);

    auto g_yyy_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 297);

    auto g_yyy_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 298);

    auto g_yyy_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 299);

    auto g_yyy_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 300);

    auto g_yyy_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 301);

    auto g_yyy_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 302);

    auto g_yyy_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 303);

    auto g_yyy_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 304);

    auto g_yyy_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 305);

    auto g_yyy_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 306);

    auto g_yyy_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 307);

    auto g_yyy_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 308);

    auto g_yyy_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 309);

    auto g_yyy_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 310);

    auto g_yyy_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 311);

    auto g_yyy_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 312);

    auto g_yyy_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 313);

    auto g_yyy_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 314);

    auto g_yyz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 316);

    auto g_yyz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 317);

    auto g_yyz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 318);

    auto g_yyz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 319);

    auto g_yyz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 320);

    auto g_yyz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 321);

    auto g_yyz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 322);

    auto g_yyz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 323);

    auto g_yyz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 324);

    auto g_yyz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 325);

    auto g_yyz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 326);

    auto g_yyz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 327);

    auto g_yyz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 328);

    auto g_yyz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 329);

    auto g_yyz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 330);

    auto g_yyz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 331);

    auto g_yyz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 332);

    auto g_yyz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 333);

    auto g_yyz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 334);

    auto g_yyz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 335);

    auto g_yyz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 336);

    auto g_yyz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 337);

    auto g_yyz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 338);

    auto g_yyz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 339);

    auto g_yyz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 340);

    auto g_yyz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 341);

    auto g_yyz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 342);

    auto g_yyz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 343);

    auto g_yyz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 344);

    auto g_yyz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 345);

    auto g_yyz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 346);

    auto g_yyz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 347);

    auto g_yyz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 348);

    auto g_yyz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 349);

    auto g_yyz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 350);

    auto g_yyz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 351);

    auto g_yyz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 352);

    auto g_yyz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 353);

    auto g_yyz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 354);

    auto g_yyz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 355);

    auto g_yyz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 356);

    auto g_yyz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 357);

    auto g_yyz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 358);

    auto g_yyz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 359);

    auto g_yzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 360);

    auto g_yzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 361);

    auto g_yzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 362);

    auto g_yzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 363);

    auto g_yzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 364);

    auto g_yzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 365);

    auto g_yzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 366);

    auto g_yzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 367);

    auto g_yzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 368);

    auto g_yzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 369);

    auto g_yzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 370);

    auto g_yzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 371);

    auto g_yzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 372);

    auto g_yzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 373);

    auto g_yzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 374);

    auto g_yzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 375);

    auto g_yzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 376);

    auto g_yzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 377);

    auto g_yzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 378);

    auto g_yzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 379);

    auto g_yzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 380);

    auto g_yzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 381);

    auto g_yzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 382);

    auto g_yzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 383);

    auto g_yzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 384);

    auto g_yzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 385);

    auto g_yzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 386);

    auto g_yzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 387);

    auto g_yzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 388);

    auto g_yzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 389);

    auto g_yzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 390);

    auto g_yzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 391);

    auto g_yzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 392);

    auto g_yzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 393);

    auto g_yzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 394);

    auto g_yzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 395);

    auto g_yzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 396);

    auto g_yzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 397);

    auto g_yzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 398);

    auto g_yzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 399);

    auto g_yzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 400);

    auto g_yzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 401);

    auto g_yzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 402);

    auto g_yzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 403);

    auto g_yzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 404);

    auto g_zzz_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_fsl + 405);

    auto g_zzz_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_fsl + 406);

    auto g_zzz_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_fsl + 407);

    auto g_zzz_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_fsl + 408);

    auto g_zzz_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_fsl + 409);

    auto g_zzz_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_fsl + 410);

    auto g_zzz_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_fsl + 411);

    auto g_zzz_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_fsl + 412);

    auto g_zzz_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_fsl + 413);

    auto g_zzz_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_fsl + 414);

    auto g_zzz_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_fsl + 415);

    auto g_zzz_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_fsl + 416);

    auto g_zzz_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_fsl + 417);

    auto g_zzz_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_fsl + 418);

    auto g_zzz_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_fsl + 419);

    auto g_zzz_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 420);

    auto g_zzz_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 421);

    auto g_zzz_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 422);

    auto g_zzz_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 423);

    auto g_zzz_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 424);

    auto g_zzz_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 425);

    auto g_zzz_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 426);

    auto g_zzz_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 427);

    auto g_zzz_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 428);

    auto g_zzz_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 429);

    auto g_zzz_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 430);

    auto g_zzz_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 431);

    auto g_zzz_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 432);

    auto g_zzz_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 433);

    auto g_zzz_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 434);

    auto g_zzz_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 435);

    auto g_zzz_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 436);

    auto g_zzz_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 437);

    auto g_zzz_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 438);

    auto g_zzz_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 439);

    auto g_zzz_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 440);

    auto g_zzz_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_fsl + 441);

    auto g_zzz_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_fsl + 442);

    auto g_zzz_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_fsl + 443);

    auto g_zzz_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_fsl + 444);

    auto g_zzz_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_fsl + 445);

    auto g_zzz_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 446);

    auto g_zzz_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 447);

    auto g_zzz_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 448);

    auto g_zzz_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_fsl + 449);

    /// Set up 0-45 components of targeted buffer : GSL

    auto g_xxxx_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl);

    auto g_xxxx_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 1);

    auto g_xxxx_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 2);

    auto g_xxxx_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 3);

    auto g_xxxx_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 4);

    auto g_xxxx_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 5);

    auto g_xxxx_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 6);

    auto g_xxxx_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 7);

    auto g_xxxx_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 8);

    auto g_xxxx_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 9);

    auto g_xxxx_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 10);

    auto g_xxxx_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 11);

    auto g_xxxx_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 12);

    auto g_xxxx_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 13);

    auto g_xxxx_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 14);

    auto g_xxxx_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 15);

    auto g_xxxx_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 16);

    auto g_xxxx_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 17);

    auto g_xxxx_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 18);

    auto g_xxxx_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 19);

    auto g_xxxx_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 20);

    auto g_xxxx_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 21);

    auto g_xxxx_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 22);

    auto g_xxxx_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 23);

    auto g_xxxx_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 24);

    auto g_xxxx_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 25);

    auto g_xxxx_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 26);

    auto g_xxxx_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 27);

    auto g_xxxx_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 28);

    auto g_xxxx_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 29);

    auto g_xxxx_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 30);

    auto g_xxxx_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 31);

    auto g_xxxx_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 32);

    auto g_xxxx_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 33);

    auto g_xxxx_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 34);

    auto g_xxxx_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 35);

    auto g_xxxx_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 36);

    auto g_xxxx_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 37);

    auto g_xxxx_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 38);

    auto g_xxxx_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 39);

    auto g_xxxx_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 40);

    auto g_xxxx_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 41);

    auto g_xxxx_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 42);

    auto g_xxxx_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 43);

    auto g_xxxx_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 44);

    #pragma omp simd aligned(g_xx_0_xxxxxxxx_0, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxy_0, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxxz_0, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxyy_0, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxxyz_0, g_xx_0_xxxxxxyz_1, g_xx_0_xxxxxxzz_0, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxyyy_0, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxxyyz_0, g_xx_0_xxxxxyyz_1, g_xx_0_xxxxxyzz_0, g_xx_0_xxxxxyzz_1, g_xx_0_xxxxxzzz_0, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxyyyy_0, g_xx_0_xxxxyyyy_1, g_xx_0_xxxxyyyz_0, g_xx_0_xxxxyyyz_1, g_xx_0_xxxxyyzz_0, g_xx_0_xxxxyyzz_1, g_xx_0_xxxxyzzz_0, g_xx_0_xxxxyzzz_1, g_xx_0_xxxxzzzz_0, g_xx_0_xxxxzzzz_1, g_xx_0_xxxyyyyy_0, g_xx_0_xxxyyyyy_1, g_xx_0_xxxyyyyz_0, g_xx_0_xxxyyyyz_1, g_xx_0_xxxyyyzz_0, g_xx_0_xxxyyyzz_1, g_xx_0_xxxyyzzz_0, g_xx_0_xxxyyzzz_1, g_xx_0_xxxyzzzz_0, g_xx_0_xxxyzzzz_1, g_xx_0_xxxzzzzz_0, g_xx_0_xxxzzzzz_1, g_xx_0_xxyyyyyy_0, g_xx_0_xxyyyyyy_1, g_xx_0_xxyyyyyz_0, g_xx_0_xxyyyyyz_1, g_xx_0_xxyyyyzz_0, g_xx_0_xxyyyyzz_1, g_xx_0_xxyyyzzz_0, g_xx_0_xxyyyzzz_1, g_xx_0_xxyyzzzz_0, g_xx_0_xxyyzzzz_1, g_xx_0_xxyzzzzz_0, g_xx_0_xxyzzzzz_1, g_xx_0_xxzzzzzz_0, g_xx_0_xxzzzzzz_1, g_xx_0_xyyyyyyy_0, g_xx_0_xyyyyyyy_1, g_xx_0_xyyyyyyz_0, g_xx_0_xyyyyyyz_1, g_xx_0_xyyyyyzz_0, g_xx_0_xyyyyyzz_1, g_xx_0_xyyyyzzz_0, g_xx_0_xyyyyzzz_1, g_xx_0_xyyyzzzz_0, g_xx_0_xyyyzzzz_1, g_xx_0_xyyzzzzz_0, g_xx_0_xyyzzzzz_1, g_xx_0_xyzzzzzz_0, g_xx_0_xyzzzzzz_1, g_xx_0_xzzzzzzz_0, g_xx_0_xzzzzzzz_1, g_xx_0_yyyyyyyy_0, g_xx_0_yyyyyyyy_1, g_xx_0_yyyyyyyz_0, g_xx_0_yyyyyyyz_1, g_xx_0_yyyyyyzz_0, g_xx_0_yyyyyyzz_1, g_xx_0_yyyyyzzz_0, g_xx_0_yyyyyzzz_1, g_xx_0_yyyyzzzz_0, g_xx_0_yyyyzzzz_1, g_xx_0_yyyzzzzz_0, g_xx_0_yyyzzzzz_1, g_xx_0_yyzzzzzz_0, g_xx_0_yyzzzzzz_1, g_xx_0_yzzzzzzz_0, g_xx_0_yzzzzzzz_1, g_xx_0_zzzzzzzz_0, g_xx_0_zzzzzzzz_1, g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxzzzz_1, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyzzzz_1, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxzzzzz_1, g_xxx_0_xxzzzzzz_1, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyzzzz_1, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyzzzzz_1, g_xxx_0_xyzzzzzz_1, g_xxx_0_xzzzzzz_1, g_xxx_0_xzzzzzzz_1, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyzzzz_1, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyzzzzz_1, g_xxx_0_yyzzzzzz_1, g_xxx_0_yzzzzzz_1, g_xxx_0_yzzzzzzz_1, g_xxx_0_zzzzzzz_1, g_xxx_0_zzzzzzzz_1, g_xxxx_0_xxxxxxxx_0, g_xxxx_0_xxxxxxxy_0, g_xxxx_0_xxxxxxxz_0, g_xxxx_0_xxxxxxyy_0, g_xxxx_0_xxxxxxyz_0, g_xxxx_0_xxxxxxzz_0, g_xxxx_0_xxxxxyyy_0, g_xxxx_0_xxxxxyyz_0, g_xxxx_0_xxxxxyzz_0, g_xxxx_0_xxxxxzzz_0, g_xxxx_0_xxxxyyyy_0, g_xxxx_0_xxxxyyyz_0, g_xxxx_0_xxxxyyzz_0, g_xxxx_0_xxxxyzzz_0, g_xxxx_0_xxxxzzzz_0, g_xxxx_0_xxxyyyyy_0, g_xxxx_0_xxxyyyyz_0, g_xxxx_0_xxxyyyzz_0, g_xxxx_0_xxxyyzzz_0, g_xxxx_0_xxxyzzzz_0, g_xxxx_0_xxxzzzzz_0, g_xxxx_0_xxyyyyyy_0, g_xxxx_0_xxyyyyyz_0, g_xxxx_0_xxyyyyzz_0, g_xxxx_0_xxyyyzzz_0, g_xxxx_0_xxyyzzzz_0, g_xxxx_0_xxyzzzzz_0, g_xxxx_0_xxzzzzzz_0, g_xxxx_0_xyyyyyyy_0, g_xxxx_0_xyyyyyyz_0, g_xxxx_0_xyyyyyzz_0, g_xxxx_0_xyyyyzzz_0, g_xxxx_0_xyyyzzzz_0, g_xxxx_0_xyyzzzzz_0, g_xxxx_0_xyzzzzzz_0, g_xxxx_0_xzzzzzzz_0, g_xxxx_0_yyyyyyyy_0, g_xxxx_0_yyyyyyyz_0, g_xxxx_0_yyyyyyzz_0, g_xxxx_0_yyyyyzzz_0, g_xxxx_0_yyyyzzzz_0, g_xxxx_0_yyyzzzzz_0, g_xxxx_0_yyzzzzzz_0, g_xxxx_0_yzzzzzzz_0, g_xxxx_0_zzzzzzzz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xxxxxxxx_0[i] = 3.0 * g_xx_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxx_1[i] * fz_be_0 + 8.0 * g_xxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxx_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxy_0[i] = 3.0 * g_xx_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxxz_0[i] = 3.0 * g_xx_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxyy_0[i] = 3.0 * g_xx_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxyz_0[i] = 3.0 * g_xx_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxxzz_0[i] = 3.0 * g_xx_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyyy_0[i] = 3.0 * g_xx_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyyz_0[i] = 3.0 * g_xx_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxyzz_0[i] = 3.0 * g_xx_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxxzzz_0[i] = 3.0 * g_xx_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyyy_0[i] = 3.0 * g_xx_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyyz_0[i] = 3.0 * g_xx_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyyzz_0[i] = 3.0 * g_xx_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxyzzz_0[i] = 3.0 * g_xx_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxxzzzz_0[i] = 3.0 * g_xx_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyyy_0[i] = 3.0 * g_xx_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyyz_0[i] = 3.0 * g_xx_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyyzz_0[i] = 3.0 * g_xx_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyyzzz_0[i] = 3.0 * g_xx_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxyzzzz_0[i] = 3.0 * g_xx_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxxzzzzz_0[i] = 3.0 * g_xx_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyyy_0[i] = 3.0 * g_xx_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyyz_0[i] = 3.0 * g_xx_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyyzz_0[i] = 3.0 * g_xx_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyyzzz_0[i] = 3.0 * g_xx_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyyzzzz_0[i] = 3.0 * g_xx_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxyzzzzz_0[i] = 3.0 * g_xx_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xxzzzzzz_0[i] = 3.0 * g_xx_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyyy_0[i] = 3.0 * g_xx_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyyz_0[i] = 3.0 * g_xx_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyyzz_0[i] = 3.0 * g_xx_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyyzzz_0[i] = 3.0 * g_xx_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyyzzzz_0[i] = 3.0 * g_xx_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyyzzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyyzzzzz_0[i] = 3.0 * g_xx_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyyzzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xyzzzzzz_0[i] = 3.0 * g_xx_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xyzzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_xzzzzzzz_0[i] = 3.0 * g_xx_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyyy_0[i] = 3.0 * g_xx_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyyy_1[i] * fz_be_0 + g_xxx_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyyz_0[i] = 3.0 * g_xx_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyyz_1[i] * fz_be_0 + g_xxx_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyyzz_0[i] = 3.0 * g_xx_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyyzz_1[i] * fz_be_0 + g_xxx_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyyzzz_0[i] = 3.0 * g_xx_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyyzzz_1[i] * fz_be_0 + g_xxx_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyyzzzz_0[i] = 3.0 * g_xx_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyyzzzz_1[i] * fz_be_0 + g_xxx_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyyzzzzz_0[i] = 3.0 * g_xx_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyyzzzzz_1[i] * fz_be_0 + g_xxx_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yyzzzzzz_0[i] = 3.0 * g_xx_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yyzzzzzz_1[i] * fz_be_0 + g_xxx_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_yzzzzzzz_0[i] = 3.0 * g_xx_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_yzzzzzzz_1[i] * fz_be_0 + g_xxx_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxxx_0_zzzzzzzz_0[i] = 3.0 * g_xx_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_xx_0_zzzzzzzz_1[i] * fz_be_0 + g_xxx_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 45-90 components of targeted buffer : GSL

    auto g_xxxy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 45);

    auto g_xxxy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 46);

    auto g_xxxy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 47);

    auto g_xxxy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 48);

    auto g_xxxy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 49);

    auto g_xxxy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 50);

    auto g_xxxy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 51);

    auto g_xxxy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 52);

    auto g_xxxy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 53);

    auto g_xxxy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 54);

    auto g_xxxy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 55);

    auto g_xxxy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 56);

    auto g_xxxy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 57);

    auto g_xxxy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 58);

    auto g_xxxy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 59);

    auto g_xxxy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 60);

    auto g_xxxy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 61);

    auto g_xxxy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 62);

    auto g_xxxy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 63);

    auto g_xxxy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 64);

    auto g_xxxy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 65);

    auto g_xxxy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 66);

    auto g_xxxy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 67);

    auto g_xxxy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 68);

    auto g_xxxy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 69);

    auto g_xxxy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 70);

    auto g_xxxy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 71);

    auto g_xxxy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 72);

    auto g_xxxy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 73);

    auto g_xxxy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 74);

    auto g_xxxy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 75);

    auto g_xxxy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 76);

    auto g_xxxy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 77);

    auto g_xxxy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 78);

    auto g_xxxy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 79);

    auto g_xxxy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 80);

    auto g_xxxy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 81);

    auto g_xxxy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 82);

    auto g_xxxy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 83);

    auto g_xxxy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 84);

    auto g_xxxy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 85);

    auto g_xxxy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 86);

    auto g_xxxy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 87);

    auto g_xxxy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 88);

    auto g_xxxy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 89);

    #pragma omp simd aligned(g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxzzzz_1, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyzzzz_1, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxzzzzz_1, g_xxx_0_xxzzzzzz_1, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyzzzz_1, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyzzzzz_1, g_xxx_0_xyzzzzzz_1, g_xxx_0_xzzzzzz_1, g_xxx_0_xzzzzzzz_1, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyzzzz_1, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyzzzzz_1, g_xxx_0_yyzzzzzz_1, g_xxx_0_yzzzzzz_1, g_xxx_0_yzzzzzzz_1, g_xxx_0_zzzzzzz_1, g_xxx_0_zzzzzzzz_1, g_xxxy_0_xxxxxxxx_0, g_xxxy_0_xxxxxxxy_0, g_xxxy_0_xxxxxxxz_0, g_xxxy_0_xxxxxxyy_0, g_xxxy_0_xxxxxxyz_0, g_xxxy_0_xxxxxxzz_0, g_xxxy_0_xxxxxyyy_0, g_xxxy_0_xxxxxyyz_0, g_xxxy_0_xxxxxyzz_0, g_xxxy_0_xxxxxzzz_0, g_xxxy_0_xxxxyyyy_0, g_xxxy_0_xxxxyyyz_0, g_xxxy_0_xxxxyyzz_0, g_xxxy_0_xxxxyzzz_0, g_xxxy_0_xxxxzzzz_0, g_xxxy_0_xxxyyyyy_0, g_xxxy_0_xxxyyyyz_0, g_xxxy_0_xxxyyyzz_0, g_xxxy_0_xxxyyzzz_0, g_xxxy_0_xxxyzzzz_0, g_xxxy_0_xxxzzzzz_0, g_xxxy_0_xxyyyyyy_0, g_xxxy_0_xxyyyyyz_0, g_xxxy_0_xxyyyyzz_0, g_xxxy_0_xxyyyzzz_0, g_xxxy_0_xxyyzzzz_0, g_xxxy_0_xxyzzzzz_0, g_xxxy_0_xxzzzzzz_0, g_xxxy_0_xyyyyyyy_0, g_xxxy_0_xyyyyyyz_0, g_xxxy_0_xyyyyyzz_0, g_xxxy_0_xyyyyzzz_0, g_xxxy_0_xyyyzzzz_0, g_xxxy_0_xyyzzzzz_0, g_xxxy_0_xyzzzzzz_0, g_xxxy_0_xzzzzzzz_0, g_xxxy_0_yyyyyyyy_0, g_xxxy_0_yyyyyyyz_0, g_xxxy_0_yyyyyyzz_0, g_xxxy_0_yyyyyzzz_0, g_xxxy_0_yyyyzzzz_0, g_xxxy_0_yyyzzzzz_0, g_xxxy_0_yyzzzzzz_0, g_xxxy_0_yzzzzzzz_0, g_xxxy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xxxxxxxx_0[i] = g_xxx_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxy_0[i] = g_xxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxxz_0[i] = g_xxx_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxyy_0[i] = 2.0 * g_xxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxyz_0[i] = g_xxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxxzz_0[i] = g_xxx_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyyy_0[i] = 3.0 * g_xxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyyz_0[i] = 2.0 * g_xxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxyzz_0[i] = g_xxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxxzzz_0[i] = g_xxx_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyyy_0[i] = 4.0 * g_xxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyyz_0[i] = 3.0 * g_xxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyyzz_0[i] = 2.0 * g_xxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxyzzz_0[i] = g_xxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxxzzzz_0[i] = g_xxx_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyyy_0[i] = 5.0 * g_xxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyyz_0[i] = 4.0 * g_xxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyyzz_0[i] = 3.0 * g_xxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyyzzz_0[i] = 2.0 * g_xxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxyzzzz_0[i] = g_xxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxxzzzzz_0[i] = g_xxx_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyyy_0[i] = 6.0 * g_xxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyyz_0[i] = 5.0 * g_xxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyyzz_0[i] = 4.0 * g_xxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyyzzz_0[i] = 3.0 * g_xxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyyzzzz_0[i] = 2.0 * g_xxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxyzzzzz_0[i] = g_xxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xxzzzzzz_0[i] = g_xxx_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyyy_0[i] = 7.0 * g_xxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyyz_0[i] = 6.0 * g_xxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyyzz_0[i] = 5.0 * g_xxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyyzzz_0[i] = 4.0 * g_xxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyyzzzz_0[i] = 3.0 * g_xxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyyzzzzz_0[i] = 2.0 * g_xxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xyzzzzzz_0[i] = g_xxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_xzzzzzzz_0[i] = g_xxx_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyyy_0[i] = 8.0 * g_xxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyy_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyyz_0[i] = 7.0 * g_xxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyyzz_0[i] = 6.0 * g_xxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyyzzz_0[i] = 5.0 * g_xxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyyzzzz_0[i] = 4.0 * g_xxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyyzzzzz_0[i] = 3.0 * g_xxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yyzzzzzz_0[i] = 2.0 * g_xxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_yzzzzzzz_0[i] = g_xxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxxy_0_zzzzzzzz_0[i] = g_xxx_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 90-135 components of targeted buffer : GSL

    auto g_xxxz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 90);

    auto g_xxxz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 91);

    auto g_xxxz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 92);

    auto g_xxxz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 93);

    auto g_xxxz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 94);

    auto g_xxxz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 95);

    auto g_xxxz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 96);

    auto g_xxxz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 97);

    auto g_xxxz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 98);

    auto g_xxxz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 99);

    auto g_xxxz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 100);

    auto g_xxxz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 101);

    auto g_xxxz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 102);

    auto g_xxxz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 103);

    auto g_xxxz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 104);

    auto g_xxxz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 105);

    auto g_xxxz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 106);

    auto g_xxxz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 107);

    auto g_xxxz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 108);

    auto g_xxxz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 109);

    auto g_xxxz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 110);

    auto g_xxxz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 111);

    auto g_xxxz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 112);

    auto g_xxxz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 113);

    auto g_xxxz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 114);

    auto g_xxxz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 115);

    auto g_xxxz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 116);

    auto g_xxxz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 117);

    auto g_xxxz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 118);

    auto g_xxxz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 119);

    auto g_xxxz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 120);

    auto g_xxxz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 121);

    auto g_xxxz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 122);

    auto g_xxxz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 123);

    auto g_xxxz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 124);

    auto g_xxxz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 125);

    auto g_xxxz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 126);

    auto g_xxxz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 127);

    auto g_xxxz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 128);

    auto g_xxxz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 129);

    auto g_xxxz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 130);

    auto g_xxxz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 131);

    auto g_xxxz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 132);

    auto g_xxxz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 133);

    auto g_xxxz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 134);

    #pragma omp simd aligned(g_xxx_0_xxxxxxx_1, g_xxx_0_xxxxxxxx_1, g_xxx_0_xxxxxxxy_1, g_xxx_0_xxxxxxxz_1, g_xxx_0_xxxxxxy_1, g_xxx_0_xxxxxxyy_1, g_xxx_0_xxxxxxyz_1, g_xxx_0_xxxxxxz_1, g_xxx_0_xxxxxxzz_1, g_xxx_0_xxxxxyy_1, g_xxx_0_xxxxxyyy_1, g_xxx_0_xxxxxyyz_1, g_xxx_0_xxxxxyz_1, g_xxx_0_xxxxxyzz_1, g_xxx_0_xxxxxzz_1, g_xxx_0_xxxxxzzz_1, g_xxx_0_xxxxyyy_1, g_xxx_0_xxxxyyyy_1, g_xxx_0_xxxxyyyz_1, g_xxx_0_xxxxyyz_1, g_xxx_0_xxxxyyzz_1, g_xxx_0_xxxxyzz_1, g_xxx_0_xxxxyzzz_1, g_xxx_0_xxxxzzz_1, g_xxx_0_xxxxzzzz_1, g_xxx_0_xxxyyyy_1, g_xxx_0_xxxyyyyy_1, g_xxx_0_xxxyyyyz_1, g_xxx_0_xxxyyyz_1, g_xxx_0_xxxyyyzz_1, g_xxx_0_xxxyyzz_1, g_xxx_0_xxxyyzzz_1, g_xxx_0_xxxyzzz_1, g_xxx_0_xxxyzzzz_1, g_xxx_0_xxxzzzz_1, g_xxx_0_xxxzzzzz_1, g_xxx_0_xxyyyyy_1, g_xxx_0_xxyyyyyy_1, g_xxx_0_xxyyyyyz_1, g_xxx_0_xxyyyyz_1, g_xxx_0_xxyyyyzz_1, g_xxx_0_xxyyyzz_1, g_xxx_0_xxyyyzzz_1, g_xxx_0_xxyyzzz_1, g_xxx_0_xxyyzzzz_1, g_xxx_0_xxyzzzz_1, g_xxx_0_xxyzzzzz_1, g_xxx_0_xxzzzzz_1, g_xxx_0_xxzzzzzz_1, g_xxx_0_xyyyyyy_1, g_xxx_0_xyyyyyyy_1, g_xxx_0_xyyyyyyz_1, g_xxx_0_xyyyyyz_1, g_xxx_0_xyyyyyzz_1, g_xxx_0_xyyyyzz_1, g_xxx_0_xyyyyzzz_1, g_xxx_0_xyyyzzz_1, g_xxx_0_xyyyzzzz_1, g_xxx_0_xyyzzzz_1, g_xxx_0_xyyzzzzz_1, g_xxx_0_xyzzzzz_1, g_xxx_0_xyzzzzzz_1, g_xxx_0_xzzzzzz_1, g_xxx_0_xzzzzzzz_1, g_xxx_0_yyyyyyy_1, g_xxx_0_yyyyyyyy_1, g_xxx_0_yyyyyyyz_1, g_xxx_0_yyyyyyz_1, g_xxx_0_yyyyyyzz_1, g_xxx_0_yyyyyzz_1, g_xxx_0_yyyyyzzz_1, g_xxx_0_yyyyzzz_1, g_xxx_0_yyyyzzzz_1, g_xxx_0_yyyzzzz_1, g_xxx_0_yyyzzzzz_1, g_xxx_0_yyzzzzz_1, g_xxx_0_yyzzzzzz_1, g_xxx_0_yzzzzzz_1, g_xxx_0_yzzzzzzz_1, g_xxx_0_zzzzzzz_1, g_xxx_0_zzzzzzzz_1, g_xxxz_0_xxxxxxxx_0, g_xxxz_0_xxxxxxxy_0, g_xxxz_0_xxxxxxxz_0, g_xxxz_0_xxxxxxyy_0, g_xxxz_0_xxxxxxyz_0, g_xxxz_0_xxxxxxzz_0, g_xxxz_0_xxxxxyyy_0, g_xxxz_0_xxxxxyyz_0, g_xxxz_0_xxxxxyzz_0, g_xxxz_0_xxxxxzzz_0, g_xxxz_0_xxxxyyyy_0, g_xxxz_0_xxxxyyyz_0, g_xxxz_0_xxxxyyzz_0, g_xxxz_0_xxxxyzzz_0, g_xxxz_0_xxxxzzzz_0, g_xxxz_0_xxxyyyyy_0, g_xxxz_0_xxxyyyyz_0, g_xxxz_0_xxxyyyzz_0, g_xxxz_0_xxxyyzzz_0, g_xxxz_0_xxxyzzzz_0, g_xxxz_0_xxxzzzzz_0, g_xxxz_0_xxyyyyyy_0, g_xxxz_0_xxyyyyyz_0, g_xxxz_0_xxyyyyzz_0, g_xxxz_0_xxyyyzzz_0, g_xxxz_0_xxyyzzzz_0, g_xxxz_0_xxyzzzzz_0, g_xxxz_0_xxzzzzzz_0, g_xxxz_0_xyyyyyyy_0, g_xxxz_0_xyyyyyyz_0, g_xxxz_0_xyyyyyzz_0, g_xxxz_0_xyyyyzzz_0, g_xxxz_0_xyyyzzzz_0, g_xxxz_0_xyyzzzzz_0, g_xxxz_0_xyzzzzzz_0, g_xxxz_0_xzzzzzzz_0, g_xxxz_0_yyyyyyyy_0, g_xxxz_0_yyyyyyyz_0, g_xxxz_0_yyyyyyzz_0, g_xxxz_0_yyyyyzzz_0, g_xxxz_0_yyyyzzzz_0, g_xxxz_0_yyyzzzzz_0, g_xxxz_0_yyzzzzzz_0, g_xxxz_0_yzzzzzzz_0, g_xxxz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xxxxxxxx_0[i] = g_xxx_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxy_0[i] = g_xxx_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxxz_0[i] = g_xxx_0_xxxxxxx_1[i] * fi_acd_0 + g_xxx_0_xxxxxxxz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxyy_0[i] = g_xxx_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxyz_0[i] = g_xxx_0_xxxxxxy_1[i] * fi_acd_0 + g_xxx_0_xxxxxxyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxxzz_0[i] = 2.0 * g_xxx_0_xxxxxxz_1[i] * fi_acd_0 + g_xxx_0_xxxxxxzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyyy_0[i] = g_xxx_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyyz_0[i] = g_xxx_0_xxxxxyy_1[i] * fi_acd_0 + g_xxx_0_xxxxxyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxyzz_0[i] = 2.0 * g_xxx_0_xxxxxyz_1[i] * fi_acd_0 + g_xxx_0_xxxxxyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxxzzz_0[i] = 3.0 * g_xxx_0_xxxxxzz_1[i] * fi_acd_0 + g_xxx_0_xxxxxzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyyy_0[i] = g_xxx_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyyz_0[i] = g_xxx_0_xxxxyyy_1[i] * fi_acd_0 + g_xxx_0_xxxxyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyyzz_0[i] = 2.0 * g_xxx_0_xxxxyyz_1[i] * fi_acd_0 + g_xxx_0_xxxxyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxyzzz_0[i] = 3.0 * g_xxx_0_xxxxyzz_1[i] * fi_acd_0 + g_xxx_0_xxxxyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxxzzzz_0[i] = 4.0 * g_xxx_0_xxxxzzz_1[i] * fi_acd_0 + g_xxx_0_xxxxzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyyy_0[i] = g_xxx_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyyz_0[i] = g_xxx_0_xxxyyyy_1[i] * fi_acd_0 + g_xxx_0_xxxyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyyzz_0[i] = 2.0 * g_xxx_0_xxxyyyz_1[i] * fi_acd_0 + g_xxx_0_xxxyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyyzzz_0[i] = 3.0 * g_xxx_0_xxxyyzz_1[i] * fi_acd_0 + g_xxx_0_xxxyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxyzzzz_0[i] = 4.0 * g_xxx_0_xxxyzzz_1[i] * fi_acd_0 + g_xxx_0_xxxyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxxzzzzz_0[i] = 5.0 * g_xxx_0_xxxzzzz_1[i] * fi_acd_0 + g_xxx_0_xxxzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyyy_0[i] = g_xxx_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyyz_0[i] = g_xxx_0_xxyyyyy_1[i] * fi_acd_0 + g_xxx_0_xxyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyyzz_0[i] = 2.0 * g_xxx_0_xxyyyyz_1[i] * fi_acd_0 + g_xxx_0_xxyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyyzzz_0[i] = 3.0 * g_xxx_0_xxyyyzz_1[i] * fi_acd_0 + g_xxx_0_xxyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyyzzzz_0[i] = 4.0 * g_xxx_0_xxyyzzz_1[i] * fi_acd_0 + g_xxx_0_xxyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxyzzzzz_0[i] = 5.0 * g_xxx_0_xxyzzzz_1[i] * fi_acd_0 + g_xxx_0_xxyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xxzzzzzz_0[i] = 6.0 * g_xxx_0_xxzzzzz_1[i] * fi_acd_0 + g_xxx_0_xxzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyyy_0[i] = g_xxx_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyyz_0[i] = g_xxx_0_xyyyyyy_1[i] * fi_acd_0 + g_xxx_0_xyyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyyzz_0[i] = 2.0 * g_xxx_0_xyyyyyz_1[i] * fi_acd_0 + g_xxx_0_xyyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyyzzz_0[i] = 3.0 * g_xxx_0_xyyyyzz_1[i] * fi_acd_0 + g_xxx_0_xyyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyyzzzz_0[i] = 4.0 * g_xxx_0_xyyyzzz_1[i] * fi_acd_0 + g_xxx_0_xyyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyyzzzzz_0[i] = 5.0 * g_xxx_0_xyyzzzz_1[i] * fi_acd_0 + g_xxx_0_xyyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xyzzzzzz_0[i] = 6.0 * g_xxx_0_xyzzzzz_1[i] * fi_acd_0 + g_xxx_0_xyzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_xzzzzzzz_0[i] = 7.0 * g_xxx_0_xzzzzzz_1[i] * fi_acd_0 + g_xxx_0_xzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyyy_0[i] = g_xxx_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyyz_0[i] = g_xxx_0_yyyyyyy_1[i] * fi_acd_0 + g_xxx_0_yyyyyyyz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyyzz_0[i] = 2.0 * g_xxx_0_yyyyyyz_1[i] * fi_acd_0 + g_xxx_0_yyyyyyzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyyzzz_0[i] = 3.0 * g_xxx_0_yyyyyzz_1[i] * fi_acd_0 + g_xxx_0_yyyyyzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyyzzzz_0[i] = 4.0 * g_xxx_0_yyyyzzz_1[i] * fi_acd_0 + g_xxx_0_yyyyzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyyzzzzz_0[i] = 5.0 * g_xxx_0_yyyzzzz_1[i] * fi_acd_0 + g_xxx_0_yyyzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yyzzzzzz_0[i] = 6.0 * g_xxx_0_yyzzzzz_1[i] * fi_acd_0 + g_xxx_0_yyzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_yzzzzzzz_0[i] = 7.0 * g_xxx_0_yzzzzzz_1[i] * fi_acd_0 + g_xxx_0_yzzzzzzz_1[i] * wa_z[i];

        g_xxxz_0_zzzzzzzz_0[i] = 8.0 * g_xxx_0_zzzzzzz_1[i] * fi_acd_0 + g_xxx_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 135-180 components of targeted buffer : GSL

    auto g_xxyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 135);

    auto g_xxyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 136);

    auto g_xxyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 137);

    auto g_xxyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 138);

    auto g_xxyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 139);

    auto g_xxyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 140);

    auto g_xxyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 141);

    auto g_xxyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 142);

    auto g_xxyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 143);

    auto g_xxyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 144);

    auto g_xxyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 145);

    auto g_xxyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 146);

    auto g_xxyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 147);

    auto g_xxyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 148);

    auto g_xxyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 149);

    auto g_xxyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 150);

    auto g_xxyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 151);

    auto g_xxyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 152);

    auto g_xxyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 153);

    auto g_xxyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 154);

    auto g_xxyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 155);

    auto g_xxyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 156);

    auto g_xxyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 157);

    auto g_xxyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 158);

    auto g_xxyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 159);

    auto g_xxyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 160);

    auto g_xxyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 161);

    auto g_xxyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 162);

    auto g_xxyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 163);

    auto g_xxyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 164);

    auto g_xxyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 165);

    auto g_xxyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 166);

    auto g_xxyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 167);

    auto g_xxyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 168);

    auto g_xxyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 169);

    auto g_xxyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 170);

    auto g_xxyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 171);

    auto g_xxyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 172);

    auto g_xxyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 173);

    auto g_xxyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 174);

    auto g_xxyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 175);

    auto g_xxyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 176);

    auto g_xxyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 177);

    auto g_xxyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 178);

    auto g_xxyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 179);

    #pragma omp simd aligned(g_xx_0_xxxxxxxx_0, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxz_0, g_xx_0_xxxxxxxz_1, g_xx_0_xxxxxxzz_0, g_xx_0_xxxxxxzz_1, g_xx_0_xxxxxzzz_0, g_xx_0_xxxxxzzz_1, g_xx_0_xxxxzzzz_0, g_xx_0_xxxxzzzz_1, g_xx_0_xxxzzzzz_0, g_xx_0_xxxzzzzz_1, g_xx_0_xxzzzzzz_0, g_xx_0_xxzzzzzz_1, g_xx_0_xzzzzzzz_0, g_xx_0_xzzzzzzz_1, g_xxy_0_xxxxxxxx_1, g_xxy_0_xxxxxxxz_1, g_xxy_0_xxxxxxzz_1, g_xxy_0_xxxxxzzz_1, g_xxy_0_xxxxzzzz_1, g_xxy_0_xxxzzzzz_1, g_xxy_0_xxzzzzzz_1, g_xxy_0_xzzzzzzz_1, g_xxyy_0_xxxxxxxx_0, g_xxyy_0_xxxxxxxy_0, g_xxyy_0_xxxxxxxz_0, g_xxyy_0_xxxxxxyy_0, g_xxyy_0_xxxxxxyz_0, g_xxyy_0_xxxxxxzz_0, g_xxyy_0_xxxxxyyy_0, g_xxyy_0_xxxxxyyz_0, g_xxyy_0_xxxxxyzz_0, g_xxyy_0_xxxxxzzz_0, g_xxyy_0_xxxxyyyy_0, g_xxyy_0_xxxxyyyz_0, g_xxyy_0_xxxxyyzz_0, g_xxyy_0_xxxxyzzz_0, g_xxyy_0_xxxxzzzz_0, g_xxyy_0_xxxyyyyy_0, g_xxyy_0_xxxyyyyz_0, g_xxyy_0_xxxyyyzz_0, g_xxyy_0_xxxyyzzz_0, g_xxyy_0_xxxyzzzz_0, g_xxyy_0_xxxzzzzz_0, g_xxyy_0_xxyyyyyy_0, g_xxyy_0_xxyyyyyz_0, g_xxyy_0_xxyyyyzz_0, g_xxyy_0_xxyyyzzz_0, g_xxyy_0_xxyyzzzz_0, g_xxyy_0_xxyzzzzz_0, g_xxyy_0_xxzzzzzz_0, g_xxyy_0_xyyyyyyy_0, g_xxyy_0_xyyyyyyz_0, g_xxyy_0_xyyyyyzz_0, g_xxyy_0_xyyyyzzz_0, g_xxyy_0_xyyyzzzz_0, g_xxyy_0_xyyzzzzz_0, g_xxyy_0_xyzzzzzz_0, g_xxyy_0_xzzzzzzz_0, g_xxyy_0_yyyyyyyy_0, g_xxyy_0_yyyyyyyz_0, g_xxyy_0_yyyyyyzz_0, g_xxyy_0_yyyyyzzz_0, g_xxyy_0_yyyyzzzz_0, g_xxyy_0_yyyzzzzz_0, g_xxyy_0_yyzzzzzz_0, g_xxyy_0_yzzzzzzz_0, g_xxyy_0_zzzzzzzz_0, g_xyy_0_xxxxxxxy_1, g_xyy_0_xxxxxxy_1, g_xyy_0_xxxxxxyy_1, g_xyy_0_xxxxxxyz_1, g_xyy_0_xxxxxyy_1, g_xyy_0_xxxxxyyy_1, g_xyy_0_xxxxxyyz_1, g_xyy_0_xxxxxyz_1, g_xyy_0_xxxxxyzz_1, g_xyy_0_xxxxyyy_1, g_xyy_0_xxxxyyyy_1, g_xyy_0_xxxxyyyz_1, g_xyy_0_xxxxyyz_1, g_xyy_0_xxxxyyzz_1, g_xyy_0_xxxxyzz_1, g_xyy_0_xxxxyzzz_1, g_xyy_0_xxxyyyy_1, g_xyy_0_xxxyyyyy_1, g_xyy_0_xxxyyyyz_1, g_xyy_0_xxxyyyz_1, g_xyy_0_xxxyyyzz_1, g_xyy_0_xxxyyzz_1, g_xyy_0_xxxyyzzz_1, g_xyy_0_xxxyzzz_1, g_xyy_0_xxxyzzzz_1, g_xyy_0_xxyyyyy_1, g_xyy_0_xxyyyyyy_1, g_xyy_0_xxyyyyyz_1, g_xyy_0_xxyyyyz_1, g_xyy_0_xxyyyyzz_1, g_xyy_0_xxyyyzz_1, g_xyy_0_xxyyyzzz_1, g_xyy_0_xxyyzzz_1, g_xyy_0_xxyyzzzz_1, g_xyy_0_xxyzzzz_1, g_xyy_0_xxyzzzzz_1, g_xyy_0_xyyyyyy_1, g_xyy_0_xyyyyyyy_1, g_xyy_0_xyyyyyyz_1, g_xyy_0_xyyyyyz_1, g_xyy_0_xyyyyyzz_1, g_xyy_0_xyyyyzz_1, g_xyy_0_xyyyyzzz_1, g_xyy_0_xyyyzzz_1, g_xyy_0_xyyyzzzz_1, g_xyy_0_xyyzzzz_1, g_xyy_0_xyyzzzzz_1, g_xyy_0_xyzzzzz_1, g_xyy_0_xyzzzzzz_1, g_xyy_0_yyyyyyy_1, g_xyy_0_yyyyyyyy_1, g_xyy_0_yyyyyyyz_1, g_xyy_0_yyyyyyz_1, g_xyy_0_yyyyyyzz_1, g_xyy_0_yyyyyzz_1, g_xyy_0_yyyyyzzz_1, g_xyy_0_yyyyzzz_1, g_xyy_0_yyyyzzzz_1, g_xyy_0_yyyzzzz_1, g_xyy_0_yyyzzzzz_1, g_xyy_0_yyzzzzz_1, g_xyy_0_yyzzzzzz_1, g_xyy_0_yzzzzzz_1, g_xyy_0_yzzzzzzz_1, g_xyy_0_zzzzzzzz_1, g_yy_0_xxxxxxxy_0, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxyy_0, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyz_0, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxyyy_0, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyz_0, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyzz_0, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxyyyy_0, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyz_0, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyzz_0, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyzzz_0, g_yy_0_xxxxyzzz_1, g_yy_0_xxxyyyyy_0, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyz_0, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyzz_0, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyzzz_0, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyzzzz_0, g_yy_0_xxxyzzzz_1, g_yy_0_xxyyyyyy_0, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyz_0, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyzz_0, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyzzz_0, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyzzzz_0, g_yy_0_xxyyzzzz_1, g_yy_0_xxyzzzzz_0, g_yy_0_xxyzzzzz_1, g_yy_0_xyyyyyyy_0, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyz_0, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyzz_0, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyzzz_0, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyzzzz_0, g_yy_0_xyyyzzzz_1, g_yy_0_xyyzzzzz_0, g_yy_0_xyyzzzzz_1, g_yy_0_xyzzzzzz_0, g_yy_0_xyzzzzzz_1, g_yy_0_yyyyyyyy_0, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyz_0, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyzz_0, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyzzz_0, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyzzzz_0, g_yy_0_yyyyzzzz_1, g_yy_0_yyyzzzzz_0, g_yy_0_yyyzzzzz_1, g_yy_0_yyzzzzzz_0, g_yy_0_yyzzzzzz_1, g_yy_0_yzzzzzzz_0, g_yy_0_yzzzzzzz_1, g_yy_0_zzzzzzzz_0, g_yy_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xxxxxxxx_0[i] = g_xx_0_xxxxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxy_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyy_0_xxxxxxxy_0[i] = g_yy_0_xxxxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxxxy_1[i] * fz_be_0 + 7.0 * g_xyy_0_xxxxxxy_1[i] * fi_acd_0 + g_xyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxxz_0[i] = g_xx_0_xxxxxxxz_0[i] * fbe_0 - g_xx_0_xxxxxxxz_1[i] * fz_be_0 + g_xxy_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyy_0_xxxxxxyy_0[i] = g_yy_0_xxxxxxyy_0[i] * fbe_0 - g_yy_0_xxxxxxyy_1[i] * fz_be_0 + 6.0 * g_xyy_0_xxxxxyy_1[i] * fi_acd_0 + g_xyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxyz_0[i] = g_yy_0_xxxxxxyz_0[i] * fbe_0 - g_yy_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xyy_0_xxxxxyz_1[i] * fi_acd_0 + g_xyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxxzz_0[i] = g_xx_0_xxxxxxzz_0[i] * fbe_0 - g_xx_0_xxxxxxzz_1[i] * fz_be_0 + g_xxy_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyy_0_xxxxxyyy_0[i] = g_yy_0_xxxxxyyy_0[i] * fbe_0 - g_yy_0_xxxxxyyy_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyyy_1[i] * fi_acd_0 + g_xyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxxyyz_0[i] = g_yy_0_xxxxxyyz_0[i] * fbe_0 - g_yy_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyyz_1[i] * fi_acd_0 + g_xyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxyzz_0[i] = g_yy_0_xxxxxyzz_0[i] * fbe_0 - g_yy_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xyy_0_xxxxyzz_1[i] * fi_acd_0 + g_xyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxxzzz_0[i] = g_xx_0_xxxxxzzz_0[i] * fbe_0 - g_xx_0_xxxxxzzz_1[i] * fz_be_0 + g_xxy_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyy_0_xxxxyyyy_0[i] = g_yy_0_xxxxyyyy_0[i] * fbe_0 - g_yy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyyy_1[i] * fi_acd_0 + g_xyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxxyyyz_0[i] = g_yy_0_xxxxyyyz_0[i] * fbe_0 - g_yy_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyyz_1[i] * fi_acd_0 + g_xyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxxyyzz_0[i] = g_yy_0_xxxxyyzz_0[i] * fbe_0 - g_yy_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyyzz_1[i] * fi_acd_0 + g_xyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxyzzz_0[i] = g_yy_0_xxxxyzzz_0[i] * fbe_0 - g_yy_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xyy_0_xxxyzzz_1[i] * fi_acd_0 + g_xyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxxzzzz_0[i] = g_xx_0_xxxxzzzz_0[i] * fbe_0 - g_xx_0_xxxxzzzz_1[i] * fz_be_0 + g_xxy_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyy_0_xxxyyyyy_0[i] = g_yy_0_xxxyyyyy_0[i] * fbe_0 - g_yy_0_xxxyyyyy_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyyy_1[i] * fi_acd_0 + g_xyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxxyyyyz_0[i] = g_yy_0_xxxyyyyz_0[i] * fbe_0 - g_yy_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyyz_1[i] * fi_acd_0 + g_xyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxxyyyzz_0[i] = g_yy_0_xxxyyyzz_0[i] * fbe_0 - g_yy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyyzz_1[i] * fi_acd_0 + g_xyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxxyyzzz_0[i] = g_yy_0_xxxyyzzz_0[i] * fbe_0 - g_yy_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyyzzz_1[i] * fi_acd_0 + g_xyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxyzzzz_0[i] = g_yy_0_xxxyzzzz_0[i] * fbe_0 - g_yy_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xyy_0_xxyzzzz_1[i] * fi_acd_0 + g_xyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxxzzzzz_0[i] = g_xx_0_xxxzzzzz_0[i] * fbe_0 - g_xx_0_xxxzzzzz_1[i] * fz_be_0 + g_xxy_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyy_0_xxyyyyyy_0[i] = g_yy_0_xxyyyyyy_0[i] * fbe_0 - g_yy_0_xxyyyyyy_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyyy_1[i] * fi_acd_0 + g_xyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xxyyyyyz_0[i] = g_yy_0_xxyyyyyz_0[i] * fbe_0 - g_yy_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyyz_1[i] * fi_acd_0 + g_xyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xxyyyyzz_0[i] = g_yy_0_xxyyyyzz_0[i] * fbe_0 - g_yy_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyyzz_1[i] * fi_acd_0 + g_xyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xxyyyzzz_0[i] = g_yy_0_xxyyyzzz_0[i] * fbe_0 - g_yy_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyyzzz_1[i] * fi_acd_0 + g_xyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xxyyzzzz_0[i] = g_yy_0_xxyyzzzz_0[i] * fbe_0 - g_yy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyyzzzz_1[i] * fi_acd_0 + g_xyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxyzzzzz_0[i] = g_yy_0_xxyzzzzz_0[i] * fbe_0 - g_yy_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xyy_0_xyzzzzz_1[i] * fi_acd_0 + g_xyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xxzzzzzz_0[i] = g_xx_0_xxzzzzzz_0[i] * fbe_0 - g_xx_0_xxzzzzzz_1[i] * fz_be_0 + g_xxy_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyy_0_xyyyyyyy_0[i] = g_yy_0_xyyyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyyy_1[i] * fi_acd_0 + g_xyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_xyyyyyyz_0[i] = g_yy_0_xyyyyyyz_0[i] * fbe_0 - g_yy_0_xyyyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyyz_1[i] * fi_acd_0 + g_xyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_xyyyyyzz_0[i] = g_yy_0_xyyyyyzz_0[i] * fbe_0 - g_yy_0_xyyyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyyzz_1[i] * fi_acd_0 + g_xyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_xyyyyzzz_0[i] = g_yy_0_xyyyyzzz_0[i] * fbe_0 - g_yy_0_xyyyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyyzzz_1[i] * fi_acd_0 + g_xyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_xyyyzzzz_0[i] = g_yy_0_xyyyzzzz_0[i] * fbe_0 - g_yy_0_xyyyzzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzzz_1[i] * fi_acd_0 + g_xyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_xyyzzzzz_0[i] = g_yy_0_xyyzzzzz_0[i] * fbe_0 - g_yy_0_xyyzzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzzz_1[i] * fi_acd_0 + g_xyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xyzzzzzz_0[i] = g_yy_0_xyzzzzzz_0[i] * fbe_0 - g_yy_0_xyzzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzzz_1[i] * fi_acd_0 + g_xyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_xzzzzzzz_0[i] = g_xx_0_xzzzzzzz_0[i] * fbe_0 - g_xx_0_xzzzzzzz_1[i] * fz_be_0 + g_xxy_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyy_0_yyyyyyyy_0[i] = g_yy_0_yyyyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyyyy_1[i] * fz_be_0 + g_xyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxyy_0_yyyyyyyz_0[i] = g_yy_0_yyyyyyyz_0[i] * fbe_0 - g_yy_0_yyyyyyyz_1[i] * fz_be_0 + g_xyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxyy_0_yyyyyyzz_0[i] = g_yy_0_yyyyyyzz_0[i] * fbe_0 - g_yy_0_yyyyyyzz_1[i] * fz_be_0 + g_xyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxyy_0_yyyyyzzz_0[i] = g_yy_0_yyyyyzzz_0[i] * fbe_0 - g_yy_0_yyyyyzzz_1[i] * fz_be_0 + g_xyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxyy_0_yyyyzzzz_0[i] = g_yy_0_yyyyzzzz_0[i] * fbe_0 - g_yy_0_yyyyzzzz_1[i] * fz_be_0 + g_xyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxyy_0_yyyzzzzz_0[i] = g_yy_0_yyyzzzzz_0[i] * fbe_0 - g_yy_0_yyyzzzzz_1[i] * fz_be_0 + g_xyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxyy_0_yyzzzzzz_0[i] = g_yy_0_yyzzzzzz_0[i] * fbe_0 - g_yy_0_yyzzzzzz_1[i] * fz_be_0 + g_xyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_yzzzzzzz_0[i] = g_yy_0_yzzzzzzz_0[i] * fbe_0 - g_yy_0_yzzzzzzz_1[i] * fz_be_0 + g_xyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxyy_0_zzzzzzzz_0[i] = g_yy_0_zzzzzzzz_0[i] * fbe_0 - g_yy_0_zzzzzzzz_1[i] * fz_be_0 + g_xyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 180-225 components of targeted buffer : GSL

    auto g_xxyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 180);

    auto g_xxyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 181);

    auto g_xxyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 182);

    auto g_xxyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 183);

    auto g_xxyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 184);

    auto g_xxyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 185);

    auto g_xxyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 186);

    auto g_xxyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 187);

    auto g_xxyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 188);

    auto g_xxyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 189);

    auto g_xxyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 190);

    auto g_xxyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 191);

    auto g_xxyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 192);

    auto g_xxyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 193);

    auto g_xxyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 194);

    auto g_xxyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 195);

    auto g_xxyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 196);

    auto g_xxyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 197);

    auto g_xxyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 198);

    auto g_xxyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 199);

    auto g_xxyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 200);

    auto g_xxyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 201);

    auto g_xxyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 202);

    auto g_xxyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 203);

    auto g_xxyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 204);

    auto g_xxyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 205);

    auto g_xxyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 206);

    auto g_xxyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 207);

    auto g_xxyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 208);

    auto g_xxyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 209);

    auto g_xxyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 210);

    auto g_xxyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 211);

    auto g_xxyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 212);

    auto g_xxyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 213);

    auto g_xxyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 214);

    auto g_xxyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 215);

    auto g_xxyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 216);

    auto g_xxyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 217);

    auto g_xxyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 218);

    auto g_xxyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 219);

    auto g_xxyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 220);

    auto g_xxyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 221);

    auto g_xxyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 222);

    auto g_xxyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 223);

    auto g_xxyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 224);

    #pragma omp simd aligned(g_xxy_0_xxxxxxxy_1, g_xxy_0_xxxxxxyy_1, g_xxy_0_xxxxxyyy_1, g_xxy_0_xxxxyyyy_1, g_xxy_0_xxxyyyyy_1, g_xxy_0_xxyyyyyy_1, g_xxy_0_xyyyyyyy_1, g_xxy_0_yyyyyyyy_1, g_xxyz_0_xxxxxxxx_0, g_xxyz_0_xxxxxxxy_0, g_xxyz_0_xxxxxxxz_0, g_xxyz_0_xxxxxxyy_0, g_xxyz_0_xxxxxxyz_0, g_xxyz_0_xxxxxxzz_0, g_xxyz_0_xxxxxyyy_0, g_xxyz_0_xxxxxyyz_0, g_xxyz_0_xxxxxyzz_0, g_xxyz_0_xxxxxzzz_0, g_xxyz_0_xxxxyyyy_0, g_xxyz_0_xxxxyyyz_0, g_xxyz_0_xxxxyyzz_0, g_xxyz_0_xxxxyzzz_0, g_xxyz_0_xxxxzzzz_0, g_xxyz_0_xxxyyyyy_0, g_xxyz_0_xxxyyyyz_0, g_xxyz_0_xxxyyyzz_0, g_xxyz_0_xxxyyzzz_0, g_xxyz_0_xxxyzzzz_0, g_xxyz_0_xxxzzzzz_0, g_xxyz_0_xxyyyyyy_0, g_xxyz_0_xxyyyyyz_0, g_xxyz_0_xxyyyyzz_0, g_xxyz_0_xxyyyzzz_0, g_xxyz_0_xxyyzzzz_0, g_xxyz_0_xxyzzzzz_0, g_xxyz_0_xxzzzzzz_0, g_xxyz_0_xyyyyyyy_0, g_xxyz_0_xyyyyyyz_0, g_xxyz_0_xyyyyyzz_0, g_xxyz_0_xyyyyzzz_0, g_xxyz_0_xyyyzzzz_0, g_xxyz_0_xyyzzzzz_0, g_xxyz_0_xyzzzzzz_0, g_xxyz_0_xzzzzzzz_0, g_xxyz_0_yyyyyyyy_0, g_xxyz_0_yyyyyyyz_0, g_xxyz_0_yyyyyyzz_0, g_xxyz_0_yyyyyzzz_0, g_xxyz_0_yyyyzzzz_0, g_xxyz_0_yyyzzzzz_0, g_xxyz_0_yyzzzzzz_0, g_xxyz_0_yzzzzzzz_0, g_xxyz_0_zzzzzzzz_0, g_xxz_0_xxxxxxxx_1, g_xxz_0_xxxxxxxz_1, g_xxz_0_xxxxxxyz_1, g_xxz_0_xxxxxxz_1, g_xxz_0_xxxxxxzz_1, g_xxz_0_xxxxxyyz_1, g_xxz_0_xxxxxyz_1, g_xxz_0_xxxxxyzz_1, g_xxz_0_xxxxxzz_1, g_xxz_0_xxxxxzzz_1, g_xxz_0_xxxxyyyz_1, g_xxz_0_xxxxyyz_1, g_xxz_0_xxxxyyzz_1, g_xxz_0_xxxxyzz_1, g_xxz_0_xxxxyzzz_1, g_xxz_0_xxxxzzz_1, g_xxz_0_xxxxzzzz_1, g_xxz_0_xxxyyyyz_1, g_xxz_0_xxxyyyz_1, g_xxz_0_xxxyyyzz_1, g_xxz_0_xxxyyzz_1, g_xxz_0_xxxyyzzz_1, g_xxz_0_xxxyzzz_1, g_xxz_0_xxxyzzzz_1, g_xxz_0_xxxzzzz_1, g_xxz_0_xxxzzzzz_1, g_xxz_0_xxyyyyyz_1, g_xxz_0_xxyyyyz_1, g_xxz_0_xxyyyyzz_1, g_xxz_0_xxyyyzz_1, g_xxz_0_xxyyyzzz_1, g_xxz_0_xxyyzzz_1, g_xxz_0_xxyyzzzz_1, g_xxz_0_xxyzzzz_1, g_xxz_0_xxyzzzzz_1, g_xxz_0_xxzzzzz_1, g_xxz_0_xxzzzzzz_1, g_xxz_0_xyyyyyyz_1, g_xxz_0_xyyyyyz_1, g_xxz_0_xyyyyyzz_1, g_xxz_0_xyyyyzz_1, g_xxz_0_xyyyyzzz_1, g_xxz_0_xyyyzzz_1, g_xxz_0_xyyyzzzz_1, g_xxz_0_xyyzzzz_1, g_xxz_0_xyyzzzzz_1, g_xxz_0_xyzzzzz_1, g_xxz_0_xyzzzzzz_1, g_xxz_0_xzzzzzz_1, g_xxz_0_xzzzzzzz_1, g_xxz_0_yyyyyyyz_1, g_xxz_0_yyyyyyz_1, g_xxz_0_yyyyyyzz_1, g_xxz_0_yyyyyzz_1, g_xxz_0_yyyyyzzz_1, g_xxz_0_yyyyzzz_1, g_xxz_0_yyyyzzzz_1, g_xxz_0_yyyzzzz_1, g_xxz_0_yyyzzzzz_1, g_xxz_0_yyzzzzz_1, g_xxz_0_yyzzzzzz_1, g_xxz_0_yzzzzzz_1, g_xxz_0_yzzzzzzz_1, g_xxz_0_zzzzzzz_1, g_xxz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xxxxxxxx_0[i] = g_xxz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxxy_0[i] = g_xxy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxxxz_0[i] = g_xxz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxyy_0[i] = g_xxy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxxyz_0[i] = g_xxz_0_xxxxxxz_1[i] * fi_acd_0 + g_xxz_0_xxxxxxyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxxzz_0[i] = g_xxz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxyyy_0[i] = g_xxy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxxyyz_0[i] = 2.0 * g_xxz_0_xxxxxyz_1[i] * fi_acd_0 + g_xxz_0_xxxxxyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxyzz_0[i] = g_xxz_0_xxxxxzz_1[i] * fi_acd_0 + g_xxz_0_xxxxxyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxxzzz_0[i] = g_xxz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyyyy_0[i] = g_xxy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxxyyyz_0[i] = 3.0 * g_xxz_0_xxxxyyz_1[i] * fi_acd_0 + g_xxz_0_xxxxyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyyzz_0[i] = 2.0 * g_xxz_0_xxxxyzz_1[i] * fi_acd_0 + g_xxz_0_xxxxyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxyzzz_0[i] = g_xxz_0_xxxxzzz_1[i] * fi_acd_0 + g_xxz_0_xxxxyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxxzzzz_0[i] = g_xxz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyyyy_0[i] = g_xxy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxxyyyyz_0[i] = 4.0 * g_xxz_0_xxxyyyz_1[i] * fi_acd_0 + g_xxz_0_xxxyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyyzz_0[i] = 3.0 * g_xxz_0_xxxyyzz_1[i] * fi_acd_0 + g_xxz_0_xxxyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyyzzz_0[i] = 2.0 * g_xxz_0_xxxyzzz_1[i] * fi_acd_0 + g_xxz_0_xxxyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxyzzzz_0[i] = g_xxz_0_xxxzzzz_1[i] * fi_acd_0 + g_xxz_0_xxxyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxxzzzzz_0[i] = g_xxz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyyyy_0[i] = g_xxy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xxyyyyyz_0[i] = 5.0 * g_xxz_0_xxyyyyz_1[i] * fi_acd_0 + g_xxz_0_xxyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyyzz_0[i] = 4.0 * g_xxz_0_xxyyyzz_1[i] * fi_acd_0 + g_xxz_0_xxyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyyzzz_0[i] = 3.0 * g_xxz_0_xxyyzzz_1[i] * fi_acd_0 + g_xxz_0_xxyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyyzzzz_0[i] = 2.0 * g_xxz_0_xxyzzzz_1[i] * fi_acd_0 + g_xxz_0_xxyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxyzzzzz_0[i] = g_xxz_0_xxzzzzz_1[i] * fi_acd_0 + g_xxz_0_xxyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xxzzzzzz_0[i] = g_xxz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyyyy_0[i] = g_xxy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_xyyyyyyz_0[i] = 6.0 * g_xxz_0_xyyyyyz_1[i] * fi_acd_0 + g_xxz_0_xyyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyyzz_0[i] = 5.0 * g_xxz_0_xyyyyzz_1[i] * fi_acd_0 + g_xxz_0_xyyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyyzzz_0[i] = 4.0 * g_xxz_0_xyyyzzz_1[i] * fi_acd_0 + g_xxz_0_xyyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyyzzzz_0[i] = 3.0 * g_xxz_0_xyyzzzz_1[i] * fi_acd_0 + g_xxz_0_xyyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyyzzzzz_0[i] = 2.0 * g_xxz_0_xyzzzzz_1[i] * fi_acd_0 + g_xxz_0_xyyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xyzzzzzz_0[i] = g_xxz_0_xzzzzzz_1[i] * fi_acd_0 + g_xxz_0_xyzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_xzzzzzzz_0[i] = g_xxz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyyyy_0[i] = g_xxy_0_yyyyyyyy_1[i] * wa_z[i];

        g_xxyz_0_yyyyyyyz_0[i] = 7.0 * g_xxz_0_yyyyyyz_1[i] * fi_acd_0 + g_xxz_0_yyyyyyyz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyyzz_0[i] = 6.0 * g_xxz_0_yyyyyzz_1[i] * fi_acd_0 + g_xxz_0_yyyyyyzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyyzzz_0[i] = 5.0 * g_xxz_0_yyyyzzz_1[i] * fi_acd_0 + g_xxz_0_yyyyyzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyyzzzz_0[i] = 4.0 * g_xxz_0_yyyzzzz_1[i] * fi_acd_0 + g_xxz_0_yyyyzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyyzzzzz_0[i] = 3.0 * g_xxz_0_yyzzzzz_1[i] * fi_acd_0 + g_xxz_0_yyyzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yyzzzzzz_0[i] = 2.0 * g_xxz_0_yzzzzzz_1[i] * fi_acd_0 + g_xxz_0_yyzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_yzzzzzzz_0[i] = g_xxz_0_zzzzzzz_1[i] * fi_acd_0 + g_xxz_0_yzzzzzzz_1[i] * wa_y[i];

        g_xxyz_0_zzzzzzzz_0[i] = g_xxz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 225-270 components of targeted buffer : GSL

    auto g_xxzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 225);

    auto g_xxzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 226);

    auto g_xxzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 227);

    auto g_xxzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 228);

    auto g_xxzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 229);

    auto g_xxzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 230);

    auto g_xxzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 231);

    auto g_xxzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 232);

    auto g_xxzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 233);

    auto g_xxzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 234);

    auto g_xxzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 235);

    auto g_xxzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 236);

    auto g_xxzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 237);

    auto g_xxzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 238);

    auto g_xxzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 239);

    auto g_xxzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 240);

    auto g_xxzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 241);

    auto g_xxzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 242);

    auto g_xxzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 243);

    auto g_xxzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 244);

    auto g_xxzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 245);

    auto g_xxzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 246);

    auto g_xxzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 247);

    auto g_xxzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 248);

    auto g_xxzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 249);

    auto g_xxzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 250);

    auto g_xxzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 251);

    auto g_xxzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 252);

    auto g_xxzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 253);

    auto g_xxzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 254);

    auto g_xxzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 255);

    auto g_xxzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 256);

    auto g_xxzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 257);

    auto g_xxzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 258);

    auto g_xxzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 259);

    auto g_xxzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 260);

    auto g_xxzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 261);

    auto g_xxzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 262);

    auto g_xxzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 263);

    auto g_xxzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 264);

    auto g_xxzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 265);

    auto g_xxzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 266);

    auto g_xxzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 267);

    auto g_xxzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 268);

    auto g_xxzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 269);

    #pragma omp simd aligned(g_xx_0_xxxxxxxx_0, g_xx_0_xxxxxxxx_1, g_xx_0_xxxxxxxy_0, g_xx_0_xxxxxxxy_1, g_xx_0_xxxxxxyy_0, g_xx_0_xxxxxxyy_1, g_xx_0_xxxxxyyy_0, g_xx_0_xxxxxyyy_1, g_xx_0_xxxxyyyy_0, g_xx_0_xxxxyyyy_1, g_xx_0_xxxyyyyy_0, g_xx_0_xxxyyyyy_1, g_xx_0_xxyyyyyy_0, g_xx_0_xxyyyyyy_1, g_xx_0_xyyyyyyy_0, g_xx_0_xyyyyyyy_1, g_xxz_0_xxxxxxxx_1, g_xxz_0_xxxxxxxy_1, g_xxz_0_xxxxxxyy_1, g_xxz_0_xxxxxyyy_1, g_xxz_0_xxxxyyyy_1, g_xxz_0_xxxyyyyy_1, g_xxz_0_xxyyyyyy_1, g_xxz_0_xyyyyyyy_1, g_xxzz_0_xxxxxxxx_0, g_xxzz_0_xxxxxxxy_0, g_xxzz_0_xxxxxxxz_0, g_xxzz_0_xxxxxxyy_0, g_xxzz_0_xxxxxxyz_0, g_xxzz_0_xxxxxxzz_0, g_xxzz_0_xxxxxyyy_0, g_xxzz_0_xxxxxyyz_0, g_xxzz_0_xxxxxyzz_0, g_xxzz_0_xxxxxzzz_0, g_xxzz_0_xxxxyyyy_0, g_xxzz_0_xxxxyyyz_0, g_xxzz_0_xxxxyyzz_0, g_xxzz_0_xxxxyzzz_0, g_xxzz_0_xxxxzzzz_0, g_xxzz_0_xxxyyyyy_0, g_xxzz_0_xxxyyyyz_0, g_xxzz_0_xxxyyyzz_0, g_xxzz_0_xxxyyzzz_0, g_xxzz_0_xxxyzzzz_0, g_xxzz_0_xxxzzzzz_0, g_xxzz_0_xxyyyyyy_0, g_xxzz_0_xxyyyyyz_0, g_xxzz_0_xxyyyyzz_0, g_xxzz_0_xxyyyzzz_0, g_xxzz_0_xxyyzzzz_0, g_xxzz_0_xxyzzzzz_0, g_xxzz_0_xxzzzzzz_0, g_xxzz_0_xyyyyyyy_0, g_xxzz_0_xyyyyyyz_0, g_xxzz_0_xyyyyyzz_0, g_xxzz_0_xyyyyzzz_0, g_xxzz_0_xyyyzzzz_0, g_xxzz_0_xyyzzzzz_0, g_xxzz_0_xyzzzzzz_0, g_xxzz_0_xzzzzzzz_0, g_xxzz_0_yyyyyyyy_0, g_xxzz_0_yyyyyyyz_0, g_xxzz_0_yyyyyyzz_0, g_xxzz_0_yyyyyzzz_0, g_xxzz_0_yyyyzzzz_0, g_xxzz_0_yyyzzzzz_0, g_xxzz_0_yyzzzzzz_0, g_xxzz_0_yzzzzzzz_0, g_xxzz_0_zzzzzzzz_0, g_xzz_0_xxxxxxxz_1, g_xzz_0_xxxxxxyz_1, g_xzz_0_xxxxxxz_1, g_xzz_0_xxxxxxzz_1, g_xzz_0_xxxxxyyz_1, g_xzz_0_xxxxxyz_1, g_xzz_0_xxxxxyzz_1, g_xzz_0_xxxxxzz_1, g_xzz_0_xxxxxzzz_1, g_xzz_0_xxxxyyyz_1, g_xzz_0_xxxxyyz_1, g_xzz_0_xxxxyyzz_1, g_xzz_0_xxxxyzz_1, g_xzz_0_xxxxyzzz_1, g_xzz_0_xxxxzzz_1, g_xzz_0_xxxxzzzz_1, g_xzz_0_xxxyyyyz_1, g_xzz_0_xxxyyyz_1, g_xzz_0_xxxyyyzz_1, g_xzz_0_xxxyyzz_1, g_xzz_0_xxxyyzzz_1, g_xzz_0_xxxyzzz_1, g_xzz_0_xxxyzzzz_1, g_xzz_0_xxxzzzz_1, g_xzz_0_xxxzzzzz_1, g_xzz_0_xxyyyyyz_1, g_xzz_0_xxyyyyz_1, g_xzz_0_xxyyyyzz_1, g_xzz_0_xxyyyzz_1, g_xzz_0_xxyyyzzz_1, g_xzz_0_xxyyzzz_1, g_xzz_0_xxyyzzzz_1, g_xzz_0_xxyzzzz_1, g_xzz_0_xxyzzzzz_1, g_xzz_0_xxzzzzz_1, g_xzz_0_xxzzzzzz_1, g_xzz_0_xyyyyyyz_1, g_xzz_0_xyyyyyz_1, g_xzz_0_xyyyyyzz_1, g_xzz_0_xyyyyzz_1, g_xzz_0_xyyyyzzz_1, g_xzz_0_xyyyzzz_1, g_xzz_0_xyyyzzzz_1, g_xzz_0_xyyzzzz_1, g_xzz_0_xyyzzzzz_1, g_xzz_0_xyzzzzz_1, g_xzz_0_xyzzzzzz_1, g_xzz_0_xzzzzzz_1, g_xzz_0_xzzzzzzz_1, g_xzz_0_yyyyyyyy_1, g_xzz_0_yyyyyyyz_1, g_xzz_0_yyyyyyz_1, g_xzz_0_yyyyyyzz_1, g_xzz_0_yyyyyzz_1, g_xzz_0_yyyyyzzz_1, g_xzz_0_yyyyzzz_1, g_xzz_0_yyyyzzzz_1, g_xzz_0_yyyzzzz_1, g_xzz_0_yyyzzzzz_1, g_xzz_0_yyzzzzz_1, g_xzz_0_yyzzzzzz_1, g_xzz_0_yzzzzzz_1, g_xzz_0_yzzzzzzz_1, g_xzz_0_zzzzzzz_1, g_xzz_0_zzzzzzzz_1, g_zz_0_xxxxxxxz_0, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxyz_0, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxzz_0, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxyyz_0, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyzz_0, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxzzz_0, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxyyyz_0, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyzz_0, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyzzz_0, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxzzzz_0, g_zz_0_xxxxzzzz_1, g_zz_0_xxxyyyyz_0, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyzz_0, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyzzz_0, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyzzzz_0, g_zz_0_xxxyzzzz_1, g_zz_0_xxxzzzzz_0, g_zz_0_xxxzzzzz_1, g_zz_0_xxyyyyyz_0, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyzz_0, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyzzz_0, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyzzzz_0, g_zz_0_xxyyzzzz_1, g_zz_0_xxyzzzzz_0, g_zz_0_xxyzzzzz_1, g_zz_0_xxzzzzzz_0, g_zz_0_xxzzzzzz_1, g_zz_0_xyyyyyyz_0, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyzz_0, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyzzz_0, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyzzzz_0, g_zz_0_xyyyzzzz_1, g_zz_0_xyyzzzzz_0, g_zz_0_xyyzzzzz_1, g_zz_0_xyzzzzzz_0, g_zz_0_xyzzzzzz_1, g_zz_0_xzzzzzzz_0, g_zz_0_xzzzzzzz_1, g_zz_0_yyyyyyyy_0, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyz_0, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyzz_0, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyzzz_0, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyzzzz_0, g_zz_0_yyyyzzzz_1, g_zz_0_yyyzzzzz_0, g_zz_0_yyyzzzzz_1, g_zz_0_yyzzzzzz_0, g_zz_0_yyzzzzzz_1, g_zz_0_yzzzzzzz_0, g_zz_0_yzzzzzzz_1, g_zz_0_zzzzzzzz_0, g_zz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xxxxxxxx_0[i] = g_xx_0_xxxxxxxx_0[i] * fbe_0 - g_xx_0_xxxxxxxx_1[i] * fz_be_0 + g_xxz_0_xxxxxxxx_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxxy_0[i] = g_xx_0_xxxxxxxy_0[i] * fbe_0 - g_xx_0_xxxxxxxy_1[i] * fz_be_0 + g_xxz_0_xxxxxxxy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxxz_0[i] = g_zz_0_xxxxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxxxz_1[i] * fz_be_0 + 7.0 * g_xzz_0_xxxxxxz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxyy_0[i] = g_xx_0_xxxxxxyy_0[i] * fbe_0 - g_xx_0_xxxxxxyy_1[i] * fz_be_0 + g_xxz_0_xxxxxxyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxxyz_0[i] = g_zz_0_xxxxxxyz_0[i] * fbe_0 - g_zz_0_xxxxxxyz_1[i] * fz_be_0 + 6.0 * g_xzz_0_xxxxxyz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxxzz_0[i] = g_zz_0_xxxxxxzz_0[i] * fbe_0 - g_zz_0_xxxxxxzz_1[i] * fz_be_0 + 6.0 * g_xzz_0_xxxxxzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxyyy_0[i] = g_xx_0_xxxxxyyy_0[i] * fbe_0 - g_xx_0_xxxxxyyy_1[i] * fz_be_0 + g_xxz_0_xxxxxyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxxyyz_0[i] = g_zz_0_xxxxxyyz_0[i] * fbe_0 - g_zz_0_xxxxxyyz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxyyz_1[i] * fi_acd_0 + g_xzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxyzz_0[i] = g_zz_0_xxxxxyzz_0[i] * fbe_0 - g_zz_0_xxxxxyzz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxyzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxxzzz_0[i] = g_zz_0_xxxxxzzz_0[i] * fbe_0 - g_zz_0_xxxxxzzz_1[i] * fz_be_0 + 5.0 * g_xzz_0_xxxxzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyyyy_0[i] = g_xx_0_xxxxyyyy_0[i] * fbe_0 - g_xx_0_xxxxyyyy_1[i] * fz_be_0 + g_xxz_0_xxxxyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxxyyyz_0[i] = g_zz_0_xxxxyyyz_0[i] * fbe_0 - g_zz_0_xxxxyyyz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyyyz_1[i] * fi_acd_0 + g_xzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyyzz_0[i] = g_zz_0_xxxxyyzz_0[i] * fbe_0 - g_zz_0_xxxxyyzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyyzz_1[i] * fi_acd_0 + g_xzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxyzzz_0[i] = g_zz_0_xxxxyzzz_0[i] * fbe_0 - g_zz_0_xxxxyzzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxyzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxxzzzz_0[i] = g_zz_0_xxxxzzzz_0[i] * fbe_0 - g_zz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_xzz_0_xxxzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyyyy_0[i] = g_xx_0_xxxyyyyy_0[i] * fbe_0 - g_xx_0_xxxyyyyy_1[i] * fz_be_0 + g_xxz_0_xxxyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxxyyyyz_0[i] = g_zz_0_xxxyyyyz_0[i] * fbe_0 - g_zz_0_xxxyyyyz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyyyz_1[i] * fi_acd_0 + g_xzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyyzz_0[i] = g_zz_0_xxxyyyzz_0[i] * fbe_0 - g_zz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyyzz_1[i] * fi_acd_0 + g_xzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyyzzz_0[i] = g_zz_0_xxxyyzzz_0[i] * fbe_0 - g_zz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyyzzz_1[i] * fi_acd_0 + g_xzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxyzzzz_0[i] = g_zz_0_xxxyzzzz_0[i] * fbe_0 - g_zz_0_xxxyzzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxyzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxxzzzzz_0[i] = g_zz_0_xxxzzzzz_0[i] * fbe_0 - g_zz_0_xxxzzzzz_1[i] * fz_be_0 + 3.0 * g_xzz_0_xxzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyyyy_0[i] = g_xx_0_xxyyyyyy_0[i] * fbe_0 - g_xx_0_xxyyyyyy_1[i] * fz_be_0 + g_xxz_0_xxyyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xxyyyyyz_0[i] = g_zz_0_xxyyyyyz_0[i] * fbe_0 - g_zz_0_xxyyyyyz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyyyz_1[i] * fi_acd_0 + g_xzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyyzz_0[i] = g_zz_0_xxyyyyzz_0[i] * fbe_0 - g_zz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyyzz_1[i] * fi_acd_0 + g_xzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyyzzz_0[i] = g_zz_0_xxyyyzzz_0[i] * fbe_0 - g_zz_0_xxyyyzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyyzzz_1[i] * fi_acd_0 + g_xzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyyzzzz_0[i] = g_zz_0_xxyyzzzz_0[i] * fbe_0 - g_zz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyyzzzz_1[i] * fi_acd_0 + g_xzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxyzzzzz_0[i] = g_zz_0_xxyzzzzz_0[i] * fbe_0 - g_zz_0_xxyzzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xyzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xxzzzzzz_0[i] = g_zz_0_xxzzzzzz_0[i] * fbe_0 - g_zz_0_xxzzzzzz_1[i] * fz_be_0 + 2.0 * g_xzz_0_xzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyyyy_0[i] = g_xx_0_xyyyyyyy_0[i] * fbe_0 - g_xx_0_xyyyyyyy_1[i] * fz_be_0 + g_xxz_0_xyyyyyyy_1[i] * wa_z[i];

        g_xxzz_0_xyyyyyyz_0[i] = g_zz_0_xyyyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyyz_1[i] * fi_acd_0 + g_xzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyyzz_0[i] = g_zz_0_xyyyyyzz_0[i] * fbe_0 - g_zz_0_xyyyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyyzz_1[i] * fi_acd_0 + g_xzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyyzzz_0[i] = g_zz_0_xyyyyzzz_0[i] * fbe_0 - g_zz_0_xyyyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyyzzz_1[i] * fi_acd_0 + g_xzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyyzzzz_0[i] = g_zz_0_xyyyzzzz_0[i] * fbe_0 - g_zz_0_xyyyzzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzzz_1[i] * fi_acd_0 + g_xzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyyzzzzz_0[i] = g_zz_0_xyyzzzzz_0[i] * fbe_0 - g_zz_0_xyyzzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzzz_1[i] * fi_acd_0 + g_xzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xyzzzzzz_0[i] = g_zz_0_xyzzzzzz_0[i] * fbe_0 - g_zz_0_xyzzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_xzzzzzzz_0[i] = g_zz_0_xzzzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzzz_1[i] * fi_acd_0 + g_xzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyyy_0[i] = g_zz_0_yyyyyyyy_0[i] * fbe_0 - g_zz_0_yyyyyyyy_1[i] * fz_be_0 + g_xzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyyz_0[i] = g_zz_0_yyyyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyyyz_1[i] * fz_be_0 + g_xzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyyzz_0[i] = g_zz_0_yyyyyyzz_0[i] * fbe_0 - g_zz_0_yyyyyyzz_1[i] * fz_be_0 + g_xzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyyzzz_0[i] = g_zz_0_yyyyyzzz_0[i] * fbe_0 - g_zz_0_yyyyyzzz_1[i] * fz_be_0 + g_xzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyyzzzz_0[i] = g_zz_0_yyyyzzzz_0[i] * fbe_0 - g_zz_0_yyyyzzzz_1[i] * fz_be_0 + g_xzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyyzzzzz_0[i] = g_zz_0_yyyzzzzz_0[i] * fbe_0 - g_zz_0_yyyzzzzz_1[i] * fz_be_0 + g_xzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yyzzzzzz_0[i] = g_zz_0_yyzzzzzz_0[i] * fbe_0 - g_zz_0_yyzzzzzz_1[i] * fz_be_0 + g_xzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_yzzzzzzz_0[i] = g_zz_0_yzzzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzzzz_1[i] * fz_be_0 + g_xzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xxzz_0_zzzzzzzz_0[i] = g_zz_0_zzzzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzzzz_1[i] * fz_be_0 + g_xzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 270-315 components of targeted buffer : GSL

    auto g_xyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 270);

    auto g_xyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 271);

    auto g_xyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 272);

    auto g_xyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 273);

    auto g_xyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 274);

    auto g_xyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 275);

    auto g_xyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 276);

    auto g_xyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 277);

    auto g_xyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 278);

    auto g_xyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 279);

    auto g_xyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 280);

    auto g_xyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 281);

    auto g_xyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 282);

    auto g_xyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 283);

    auto g_xyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 284);

    auto g_xyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 285);

    auto g_xyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 286);

    auto g_xyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 287);

    auto g_xyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 288);

    auto g_xyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 289);

    auto g_xyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 290);

    auto g_xyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 291);

    auto g_xyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 292);

    auto g_xyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 293);

    auto g_xyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 294);

    auto g_xyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 295);

    auto g_xyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 296);

    auto g_xyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 297);

    auto g_xyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 298);

    auto g_xyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 299);

    auto g_xyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 300);

    auto g_xyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 301);

    auto g_xyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 302);

    auto g_xyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 303);

    auto g_xyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 304);

    auto g_xyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 305);

    auto g_xyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 306);

    auto g_xyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 307);

    auto g_xyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 308);

    auto g_xyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 309);

    auto g_xyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 310);

    auto g_xyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 311);

    auto g_xyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 312);

    auto g_xyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 313);

    auto g_xyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 314);

    #pragma omp simd aligned(g_xyyy_0_xxxxxxxx_0, g_xyyy_0_xxxxxxxy_0, g_xyyy_0_xxxxxxxz_0, g_xyyy_0_xxxxxxyy_0, g_xyyy_0_xxxxxxyz_0, g_xyyy_0_xxxxxxzz_0, g_xyyy_0_xxxxxyyy_0, g_xyyy_0_xxxxxyyz_0, g_xyyy_0_xxxxxyzz_0, g_xyyy_0_xxxxxzzz_0, g_xyyy_0_xxxxyyyy_0, g_xyyy_0_xxxxyyyz_0, g_xyyy_0_xxxxyyzz_0, g_xyyy_0_xxxxyzzz_0, g_xyyy_0_xxxxzzzz_0, g_xyyy_0_xxxyyyyy_0, g_xyyy_0_xxxyyyyz_0, g_xyyy_0_xxxyyyzz_0, g_xyyy_0_xxxyyzzz_0, g_xyyy_0_xxxyzzzz_0, g_xyyy_0_xxxzzzzz_0, g_xyyy_0_xxyyyyyy_0, g_xyyy_0_xxyyyyyz_0, g_xyyy_0_xxyyyyzz_0, g_xyyy_0_xxyyyzzz_0, g_xyyy_0_xxyyzzzz_0, g_xyyy_0_xxyzzzzz_0, g_xyyy_0_xxzzzzzz_0, g_xyyy_0_xyyyyyyy_0, g_xyyy_0_xyyyyyyz_0, g_xyyy_0_xyyyyyzz_0, g_xyyy_0_xyyyyzzz_0, g_xyyy_0_xyyyzzzz_0, g_xyyy_0_xyyzzzzz_0, g_xyyy_0_xyzzzzzz_0, g_xyyy_0_xzzzzzzz_0, g_xyyy_0_yyyyyyyy_0, g_xyyy_0_yyyyyyyz_0, g_xyyy_0_yyyyyyzz_0, g_xyyy_0_yyyyyzzz_0, g_xyyy_0_yyyyzzzz_0, g_xyyy_0_yyyzzzzz_0, g_xyyy_0_yyzzzzzz_0, g_xyyy_0_yzzzzzzz_0, g_xyyy_0_zzzzzzzz_0, g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxzzzz_1, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyzzzz_1, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxzzzzz_1, g_yyy_0_xxzzzzzz_1, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyzzzz_1, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyzzzzz_1, g_yyy_0_xyzzzzzz_1, g_yyy_0_xzzzzzz_1, g_yyy_0_xzzzzzzz_1, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyzzzz_1, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyzzzzz_1, g_yyy_0_yyzzzzzz_1, g_yyy_0_yzzzzzz_1, g_yyy_0_yzzzzzzz_1, g_yyy_0_zzzzzzz_1, g_yyy_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xxxxxxxx_0[i] = 8.0 * g_yyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxx_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxy_0[i] = 7.0 * g_yyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxxz_0[i] = 7.0 * g_yyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxyy_0[i] = 6.0 * g_yyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxyz_0[i] = 6.0 * g_yyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxxzz_0[i] = 6.0 * g_yyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyyy_0[i] = 5.0 * g_yyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyyz_0[i] = 5.0 * g_yyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxyzz_0[i] = 5.0 * g_yyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxxzzz_0[i] = 5.0 * g_yyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyyy_0[i] = 4.0 * g_yyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyyz_0[i] = 4.0 * g_yyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyyzz_0[i] = 4.0 * g_yyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxyzzz_0[i] = 4.0 * g_yyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxxzzzz_0[i] = 4.0 * g_yyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyyy_0[i] = 3.0 * g_yyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyyz_0[i] = 3.0 * g_yyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyyzz_0[i] = 3.0 * g_yyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyyzzz_0[i] = 3.0 * g_yyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxyzzzz_0[i] = 3.0 * g_yyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxxzzzzz_0[i] = 3.0 * g_yyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyyy_0[i] = 2.0 * g_yyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyyz_0[i] = 2.0 * g_yyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyyzz_0[i] = 2.0 * g_yyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyyzzz_0[i] = 2.0 * g_yyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyyzzzz_0[i] = 2.0 * g_yyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxyzzzzz_0[i] = 2.0 * g_yyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xxzzzzzz_0[i] = 2.0 * g_yyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyyy_0[i] = g_yyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyyz_0[i] = g_yyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyyzz_0[i] = g_yyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyyzzz_0[i] = g_yyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyyzzzz_0[i] = g_yyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyyzzzzz_0[i] = g_yyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xyzzzzzz_0[i] = g_yyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_xzzzzzzz_0[i] = g_yyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyyy_0[i] = g_yyy_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyyz_0[i] = g_yyy_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyyzz_0[i] = g_yyy_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyyzzz_0[i] = g_yyy_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyyzzzz_0[i] = g_yyy_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyyzzzzz_0[i] = g_yyy_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yyzzzzzz_0[i] = g_yyy_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_yzzzzzzz_0[i] = g_yyy_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyy_0_zzzzzzzz_0[i] = g_yyy_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 315-360 components of targeted buffer : GSL

    auto g_xyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 315);

    auto g_xyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 316);

    auto g_xyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 317);

    auto g_xyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 318);

    auto g_xyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 319);

    auto g_xyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 320);

    auto g_xyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 321);

    auto g_xyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 322);

    auto g_xyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 323);

    auto g_xyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 324);

    auto g_xyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 325);

    auto g_xyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 326);

    auto g_xyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 327);

    auto g_xyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 328);

    auto g_xyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 329);

    auto g_xyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 330);

    auto g_xyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 331);

    auto g_xyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 332);

    auto g_xyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 333);

    auto g_xyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 334);

    auto g_xyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 335);

    auto g_xyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 336);

    auto g_xyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 337);

    auto g_xyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 338);

    auto g_xyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 339);

    auto g_xyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 340);

    auto g_xyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 341);

    auto g_xyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 342);

    auto g_xyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 343);

    auto g_xyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 344);

    auto g_xyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 345);

    auto g_xyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 346);

    auto g_xyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 347);

    auto g_xyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 348);

    auto g_xyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 349);

    auto g_xyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 350);

    auto g_xyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 351);

    auto g_xyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 352);

    auto g_xyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 353);

    auto g_xyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 354);

    auto g_xyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 355);

    auto g_xyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 356);

    auto g_xyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 357);

    auto g_xyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 358);

    auto g_xyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 359);

    #pragma omp simd aligned(g_xyy_0_xxxxxxxx_1, g_xyy_0_xxxxxxxy_1, g_xyy_0_xxxxxxyy_1, g_xyy_0_xxxxxyyy_1, g_xyy_0_xxxxyyyy_1, g_xyy_0_xxxyyyyy_1, g_xyy_0_xxyyyyyy_1, g_xyy_0_xyyyyyyy_1, g_xyyz_0_xxxxxxxx_0, g_xyyz_0_xxxxxxxy_0, g_xyyz_0_xxxxxxxz_0, g_xyyz_0_xxxxxxyy_0, g_xyyz_0_xxxxxxyz_0, g_xyyz_0_xxxxxxzz_0, g_xyyz_0_xxxxxyyy_0, g_xyyz_0_xxxxxyyz_0, g_xyyz_0_xxxxxyzz_0, g_xyyz_0_xxxxxzzz_0, g_xyyz_0_xxxxyyyy_0, g_xyyz_0_xxxxyyyz_0, g_xyyz_0_xxxxyyzz_0, g_xyyz_0_xxxxyzzz_0, g_xyyz_0_xxxxzzzz_0, g_xyyz_0_xxxyyyyy_0, g_xyyz_0_xxxyyyyz_0, g_xyyz_0_xxxyyyzz_0, g_xyyz_0_xxxyyzzz_0, g_xyyz_0_xxxyzzzz_0, g_xyyz_0_xxxzzzzz_0, g_xyyz_0_xxyyyyyy_0, g_xyyz_0_xxyyyyyz_0, g_xyyz_0_xxyyyyzz_0, g_xyyz_0_xxyyyzzz_0, g_xyyz_0_xxyyzzzz_0, g_xyyz_0_xxyzzzzz_0, g_xyyz_0_xxzzzzzz_0, g_xyyz_0_xyyyyyyy_0, g_xyyz_0_xyyyyyyz_0, g_xyyz_0_xyyyyyzz_0, g_xyyz_0_xyyyyzzz_0, g_xyyz_0_xyyyzzzz_0, g_xyyz_0_xyyzzzzz_0, g_xyyz_0_xyzzzzzz_0, g_xyyz_0_xzzzzzzz_0, g_xyyz_0_yyyyyyyy_0, g_xyyz_0_yyyyyyyz_0, g_xyyz_0_yyyyyyzz_0, g_xyyz_0_yyyyyzzz_0, g_xyyz_0_yyyyzzzz_0, g_xyyz_0_yyyzzzzz_0, g_xyyz_0_yyzzzzzz_0, g_xyyz_0_yzzzzzzz_0, g_xyyz_0_zzzzzzzz_0, g_yyz_0_xxxxxxxz_1, g_yyz_0_xxxxxxyz_1, g_yyz_0_xxxxxxz_1, g_yyz_0_xxxxxxzz_1, g_yyz_0_xxxxxyyz_1, g_yyz_0_xxxxxyz_1, g_yyz_0_xxxxxyzz_1, g_yyz_0_xxxxxzz_1, g_yyz_0_xxxxxzzz_1, g_yyz_0_xxxxyyyz_1, g_yyz_0_xxxxyyz_1, g_yyz_0_xxxxyyzz_1, g_yyz_0_xxxxyzz_1, g_yyz_0_xxxxyzzz_1, g_yyz_0_xxxxzzz_1, g_yyz_0_xxxxzzzz_1, g_yyz_0_xxxyyyyz_1, g_yyz_0_xxxyyyz_1, g_yyz_0_xxxyyyzz_1, g_yyz_0_xxxyyzz_1, g_yyz_0_xxxyyzzz_1, g_yyz_0_xxxyzzz_1, g_yyz_0_xxxyzzzz_1, g_yyz_0_xxxzzzz_1, g_yyz_0_xxxzzzzz_1, g_yyz_0_xxyyyyyz_1, g_yyz_0_xxyyyyz_1, g_yyz_0_xxyyyyzz_1, g_yyz_0_xxyyyzz_1, g_yyz_0_xxyyyzzz_1, g_yyz_0_xxyyzzz_1, g_yyz_0_xxyyzzzz_1, g_yyz_0_xxyzzzz_1, g_yyz_0_xxyzzzzz_1, g_yyz_0_xxzzzzz_1, g_yyz_0_xxzzzzzz_1, g_yyz_0_xyyyyyyz_1, g_yyz_0_xyyyyyz_1, g_yyz_0_xyyyyyzz_1, g_yyz_0_xyyyyzz_1, g_yyz_0_xyyyyzzz_1, g_yyz_0_xyyyzzz_1, g_yyz_0_xyyyzzzz_1, g_yyz_0_xyyzzzz_1, g_yyz_0_xyyzzzzz_1, g_yyz_0_xyzzzzz_1, g_yyz_0_xyzzzzzz_1, g_yyz_0_xzzzzzz_1, g_yyz_0_xzzzzzzz_1, g_yyz_0_yyyyyyyy_1, g_yyz_0_yyyyyyyz_1, g_yyz_0_yyyyyyz_1, g_yyz_0_yyyyyyzz_1, g_yyz_0_yyyyyzz_1, g_yyz_0_yyyyyzzz_1, g_yyz_0_yyyyzzz_1, g_yyz_0_yyyyzzzz_1, g_yyz_0_yyyzzzz_1, g_yyz_0_yyyzzzzz_1, g_yyz_0_yyzzzzz_1, g_yyz_0_yyzzzzzz_1, g_yyz_0_yzzzzzz_1, g_yyz_0_yzzzzzzz_1, g_yyz_0_zzzzzzz_1, g_yyz_0_zzzzzzzz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xxxxxxxx_0[i] = g_xyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxxy_0[i] = g_xyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxxz_0[i] = 7.0 * g_yyz_0_xxxxxxz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxyy_0[i] = g_xyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxxyz_0[i] = 6.0 * g_yyz_0_xxxxxyz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxxzz_0[i] = 6.0 * g_yyz_0_xxxxxzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxyyy_0[i] = g_xyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxxyyz_0[i] = 5.0 * g_yyz_0_xxxxyyz_1[i] * fi_acd_0 + g_yyz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxyzz_0[i] = 5.0 * g_yyz_0_xxxxyzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxxzzz_0[i] = 5.0 * g_yyz_0_xxxxzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyyyy_0[i] = g_xyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxxyyyz_0[i] = 4.0 * g_yyz_0_xxxyyyz_1[i] * fi_acd_0 + g_yyz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyyzz_0[i] = 4.0 * g_yyz_0_xxxyyzz_1[i] * fi_acd_0 + g_yyz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxyzzz_0[i] = 4.0 * g_yyz_0_xxxyzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyz_0_xxxzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyyyy_0[i] = g_xyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxxyyyyz_0[i] = 3.0 * g_yyz_0_xxyyyyz_1[i] * fi_acd_0 + g_yyz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyyzz_0[i] = 3.0 * g_yyz_0_xxyyyzz_1[i] * fi_acd_0 + g_yyz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyz_0_xxyyzzz_1[i] * fi_acd_0 + g_yyz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxyzzzz_0[i] = 3.0 * g_yyz_0_xxyzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxxzzzzz_0[i] = 3.0 * g_yyz_0_xxzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyyyy_0[i] = g_xyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xxyyyyyz_0[i] = 2.0 * g_yyz_0_xyyyyyz_1[i] * fi_acd_0 + g_yyz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyz_0_xyyyyzz_1[i] * fi_acd_0 + g_yyz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyyzzz_0[i] = 2.0 * g_yyz_0_xyyyzzz_1[i] * fi_acd_0 + g_yyz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyyzzzz_0[i] = 2.0 * g_yyz_0_xyyzzzz_1[i] * fi_acd_0 + g_yyz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxyzzzzz_0[i] = 2.0 * g_yyz_0_xyzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xxzzzzzz_0[i] = 2.0 * g_yyz_0_xzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyyyy_0[i] = g_xyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_xyyz_0_xyyyyyyz_0[i] = g_yyz_0_yyyyyyz_1[i] * fi_acd_0 + g_yyz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyyzz_0[i] = g_yyz_0_yyyyyzz_1[i] * fi_acd_0 + g_yyz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyyzzz_0[i] = g_yyz_0_yyyyzzz_1[i] * fi_acd_0 + g_yyz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyyzzzz_0[i] = g_yyz_0_yyyzzzz_1[i] * fi_acd_0 + g_yyz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyyzzzzz_0[i] = g_yyz_0_yyzzzzz_1[i] * fi_acd_0 + g_yyz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xyzzzzzz_0[i] = g_yyz_0_yzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_xzzzzzzz_0[i] = g_yyz_0_zzzzzzz_1[i] * fi_acd_0 + g_yyz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyyy_0[i] = g_yyz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyyz_0[i] = g_yyz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyyzz_0[i] = g_yyz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyyzzz_0[i] = g_yyz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyyzzzz_0[i] = g_yyz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyyzzzzz_0[i] = g_yyz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yyzzzzzz_0[i] = g_yyz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_yzzzzzzz_0[i] = g_yyz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyyz_0_zzzzzzzz_0[i] = g_yyz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 360-405 components of targeted buffer : GSL

    auto g_xyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 360);

    auto g_xyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 361);

    auto g_xyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 362);

    auto g_xyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 363);

    auto g_xyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 364);

    auto g_xyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 365);

    auto g_xyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 366);

    auto g_xyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 367);

    auto g_xyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 368);

    auto g_xyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 369);

    auto g_xyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 370);

    auto g_xyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 371);

    auto g_xyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 372);

    auto g_xyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 373);

    auto g_xyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 374);

    auto g_xyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 375);

    auto g_xyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 376);

    auto g_xyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 377);

    auto g_xyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 378);

    auto g_xyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 379);

    auto g_xyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 380);

    auto g_xyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 381);

    auto g_xyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 382);

    auto g_xyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 383);

    auto g_xyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 384);

    auto g_xyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 385);

    auto g_xyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 386);

    auto g_xyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 387);

    auto g_xyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 388);

    auto g_xyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 389);

    auto g_xyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 390);

    auto g_xyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 391);

    auto g_xyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 392);

    auto g_xyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 393);

    auto g_xyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 394);

    auto g_xyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 395);

    auto g_xyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 396);

    auto g_xyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 397);

    auto g_xyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 398);

    auto g_xyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 399);

    auto g_xyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 400);

    auto g_xyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 401);

    auto g_xyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 402);

    auto g_xyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 403);

    auto g_xyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 404);

    #pragma omp simd aligned(g_xyzz_0_xxxxxxxx_0, g_xyzz_0_xxxxxxxy_0, g_xyzz_0_xxxxxxxz_0, g_xyzz_0_xxxxxxyy_0, g_xyzz_0_xxxxxxyz_0, g_xyzz_0_xxxxxxzz_0, g_xyzz_0_xxxxxyyy_0, g_xyzz_0_xxxxxyyz_0, g_xyzz_0_xxxxxyzz_0, g_xyzz_0_xxxxxzzz_0, g_xyzz_0_xxxxyyyy_0, g_xyzz_0_xxxxyyyz_0, g_xyzz_0_xxxxyyzz_0, g_xyzz_0_xxxxyzzz_0, g_xyzz_0_xxxxzzzz_0, g_xyzz_0_xxxyyyyy_0, g_xyzz_0_xxxyyyyz_0, g_xyzz_0_xxxyyyzz_0, g_xyzz_0_xxxyyzzz_0, g_xyzz_0_xxxyzzzz_0, g_xyzz_0_xxxzzzzz_0, g_xyzz_0_xxyyyyyy_0, g_xyzz_0_xxyyyyyz_0, g_xyzz_0_xxyyyyzz_0, g_xyzz_0_xxyyyzzz_0, g_xyzz_0_xxyyzzzz_0, g_xyzz_0_xxyzzzzz_0, g_xyzz_0_xxzzzzzz_0, g_xyzz_0_xyyyyyyy_0, g_xyzz_0_xyyyyyyz_0, g_xyzz_0_xyyyyyzz_0, g_xyzz_0_xyyyyzzz_0, g_xyzz_0_xyyyzzzz_0, g_xyzz_0_xyyzzzzz_0, g_xyzz_0_xyzzzzzz_0, g_xyzz_0_xzzzzzzz_0, g_xyzz_0_yyyyyyyy_0, g_xyzz_0_yyyyyyyz_0, g_xyzz_0_yyyyyyzz_0, g_xyzz_0_yyyyyzzz_0, g_xyzz_0_yyyyzzzz_0, g_xyzz_0_yyyzzzzz_0, g_xyzz_0_yyzzzzzz_0, g_xyzz_0_yzzzzzzz_0, g_xyzz_0_zzzzzzzz_0, g_xzz_0_xxxxxxxx_1, g_xzz_0_xxxxxxxz_1, g_xzz_0_xxxxxxzz_1, g_xzz_0_xxxxxzzz_1, g_xzz_0_xxxxzzzz_1, g_xzz_0_xxxzzzzz_1, g_xzz_0_xxzzzzzz_1, g_xzz_0_xzzzzzzz_1, g_yzz_0_xxxxxxxy_1, g_yzz_0_xxxxxxy_1, g_yzz_0_xxxxxxyy_1, g_yzz_0_xxxxxxyz_1, g_yzz_0_xxxxxyy_1, g_yzz_0_xxxxxyyy_1, g_yzz_0_xxxxxyyz_1, g_yzz_0_xxxxxyz_1, g_yzz_0_xxxxxyzz_1, g_yzz_0_xxxxyyy_1, g_yzz_0_xxxxyyyy_1, g_yzz_0_xxxxyyyz_1, g_yzz_0_xxxxyyz_1, g_yzz_0_xxxxyyzz_1, g_yzz_0_xxxxyzz_1, g_yzz_0_xxxxyzzz_1, g_yzz_0_xxxyyyy_1, g_yzz_0_xxxyyyyy_1, g_yzz_0_xxxyyyyz_1, g_yzz_0_xxxyyyz_1, g_yzz_0_xxxyyyzz_1, g_yzz_0_xxxyyzz_1, g_yzz_0_xxxyyzzz_1, g_yzz_0_xxxyzzz_1, g_yzz_0_xxxyzzzz_1, g_yzz_0_xxyyyyy_1, g_yzz_0_xxyyyyyy_1, g_yzz_0_xxyyyyyz_1, g_yzz_0_xxyyyyz_1, g_yzz_0_xxyyyyzz_1, g_yzz_0_xxyyyzz_1, g_yzz_0_xxyyyzzz_1, g_yzz_0_xxyyzzz_1, g_yzz_0_xxyyzzzz_1, g_yzz_0_xxyzzzz_1, g_yzz_0_xxyzzzzz_1, g_yzz_0_xyyyyyy_1, g_yzz_0_xyyyyyyy_1, g_yzz_0_xyyyyyyz_1, g_yzz_0_xyyyyyz_1, g_yzz_0_xyyyyyzz_1, g_yzz_0_xyyyyzz_1, g_yzz_0_xyyyyzzz_1, g_yzz_0_xyyyzzz_1, g_yzz_0_xyyyzzzz_1, g_yzz_0_xyyzzzz_1, g_yzz_0_xyyzzzzz_1, g_yzz_0_xyzzzzz_1, g_yzz_0_xyzzzzzz_1, g_yzz_0_yyyyyyy_1, g_yzz_0_yyyyyyyy_1, g_yzz_0_yyyyyyyz_1, g_yzz_0_yyyyyyz_1, g_yzz_0_yyyyyyzz_1, g_yzz_0_yyyyyzz_1, g_yzz_0_yyyyyzzz_1, g_yzz_0_yyyyzzz_1, g_yzz_0_yyyyzzzz_1, g_yzz_0_yyyzzzz_1, g_yzz_0_yyyzzzzz_1, g_yzz_0_yyzzzzz_1, g_yzz_0_yyzzzzzz_1, g_yzz_0_yzzzzzz_1, g_yzz_0_yzzzzzzz_1, g_yzz_0_zzzzzzzz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xxxxxxxx_0[i] = g_xzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_xyzz_0_xxxxxxxy_0[i] = 7.0 * g_yzz_0_xxxxxxy_1[i] * fi_acd_0 + g_yzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxxz_0[i] = g_xzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_xyzz_0_xxxxxxyy_0[i] = 6.0 * g_yzz_0_xxxxxyy_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxyz_0[i] = 6.0 * g_yzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxxzz_0[i] = g_xzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_xyzz_0_xxxxxyyy_0[i] = 5.0 * g_yzz_0_xxxxyyy_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxxyyz_0[i] = 5.0 * g_yzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxyzz_0[i] = 5.0 * g_yzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxxzzz_0[i] = g_xzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_xyzz_0_xxxxyyyy_0[i] = 4.0 * g_yzz_0_xxxyyyy_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxxyyyz_0[i] = 4.0 * g_yzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxxyyzz_0[i] = 4.0 * g_yzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxyzzz_0[i] = 4.0 * g_yzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxxzzzz_0[i] = g_xzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_xyzz_0_xxxyyyyy_0[i] = 3.0 * g_yzz_0_xxyyyyy_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxxyyyyz_0[i] = 3.0 * g_yzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxxyyyzz_0[i] = 3.0 * g_yzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxxyyzzz_0[i] = 3.0 * g_yzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxyzzzz_0[i] = 3.0 * g_yzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxxzzzzz_0[i] = g_xzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_xyzz_0_xxyyyyyy_0[i] = 2.0 * g_yzz_0_xyyyyyy_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xxyyyyyz_0[i] = 2.0 * g_yzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xxyyyyzz_0[i] = 2.0 * g_yzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xxyyyzzz_0[i] = 2.0 * g_yzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xxyyzzzz_0[i] = 2.0 * g_yzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxyzzzzz_0[i] = 2.0 * g_yzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xxzzzzzz_0[i] = g_xzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_xyzz_0_xyyyyyyy_0[i] = g_yzz_0_yyyyyyy_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_xyyyyyyz_0[i] = g_yzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_xyyyyyzz_0[i] = g_yzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_xyyyyzzz_0[i] = g_yzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_xyyyzzzz_0[i] = g_yzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_xyyzzzzz_0[i] = g_yzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xyzzzzzz_0[i] = g_yzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_xzzzzzzz_0[i] = g_xzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_xyzz_0_yyyyyyyy_0[i] = g_yzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xyzz_0_yyyyyyyz_0[i] = g_yzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xyzz_0_yyyyyyzz_0[i] = g_yzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xyzz_0_yyyyyzzz_0[i] = g_yzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xyzz_0_yyyyzzzz_0[i] = g_yzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xyzz_0_yyyzzzzz_0[i] = g_yzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xyzz_0_yyzzzzzz_0[i] = g_yzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_yzzzzzzz_0[i] = g_yzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xyzz_0_zzzzzzzz_0[i] = g_yzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 405-450 components of targeted buffer : GSL

    auto g_xzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 405);

    auto g_xzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 406);

    auto g_xzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 407);

    auto g_xzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 408);

    auto g_xzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 409);

    auto g_xzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 410);

    auto g_xzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 411);

    auto g_xzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 412);

    auto g_xzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 413);

    auto g_xzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 414);

    auto g_xzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 415);

    auto g_xzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 416);

    auto g_xzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 417);

    auto g_xzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 418);

    auto g_xzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 419);

    auto g_xzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 420);

    auto g_xzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 421);

    auto g_xzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 422);

    auto g_xzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 423);

    auto g_xzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 424);

    auto g_xzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 425);

    auto g_xzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 426);

    auto g_xzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 427);

    auto g_xzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 428);

    auto g_xzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 429);

    auto g_xzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 430);

    auto g_xzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 431);

    auto g_xzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 432);

    auto g_xzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 433);

    auto g_xzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 434);

    auto g_xzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 435);

    auto g_xzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 436);

    auto g_xzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 437);

    auto g_xzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 438);

    auto g_xzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 439);

    auto g_xzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 440);

    auto g_xzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 441);

    auto g_xzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 442);

    auto g_xzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 443);

    auto g_xzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 444);

    auto g_xzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 445);

    auto g_xzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 446);

    auto g_xzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 447);

    auto g_xzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 448);

    auto g_xzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 449);

    #pragma omp simd aligned(g_xzzz_0_xxxxxxxx_0, g_xzzz_0_xxxxxxxy_0, g_xzzz_0_xxxxxxxz_0, g_xzzz_0_xxxxxxyy_0, g_xzzz_0_xxxxxxyz_0, g_xzzz_0_xxxxxxzz_0, g_xzzz_0_xxxxxyyy_0, g_xzzz_0_xxxxxyyz_0, g_xzzz_0_xxxxxyzz_0, g_xzzz_0_xxxxxzzz_0, g_xzzz_0_xxxxyyyy_0, g_xzzz_0_xxxxyyyz_0, g_xzzz_0_xxxxyyzz_0, g_xzzz_0_xxxxyzzz_0, g_xzzz_0_xxxxzzzz_0, g_xzzz_0_xxxyyyyy_0, g_xzzz_0_xxxyyyyz_0, g_xzzz_0_xxxyyyzz_0, g_xzzz_0_xxxyyzzz_0, g_xzzz_0_xxxyzzzz_0, g_xzzz_0_xxxzzzzz_0, g_xzzz_0_xxyyyyyy_0, g_xzzz_0_xxyyyyyz_0, g_xzzz_0_xxyyyyzz_0, g_xzzz_0_xxyyyzzz_0, g_xzzz_0_xxyyzzzz_0, g_xzzz_0_xxyzzzzz_0, g_xzzz_0_xxzzzzzz_0, g_xzzz_0_xyyyyyyy_0, g_xzzz_0_xyyyyyyz_0, g_xzzz_0_xyyyyyzz_0, g_xzzz_0_xyyyyzzz_0, g_xzzz_0_xyyyzzzz_0, g_xzzz_0_xyyzzzzz_0, g_xzzz_0_xyzzzzzz_0, g_xzzz_0_xzzzzzzz_0, g_xzzz_0_yyyyyyyy_0, g_xzzz_0_yyyyyyyz_0, g_xzzz_0_yyyyyyzz_0, g_xzzz_0_yyyyyzzz_0, g_xzzz_0_yyyyzzzz_0, g_xzzz_0_yyyzzzzz_0, g_xzzz_0_yyzzzzzz_0, g_xzzz_0_yzzzzzzz_0, g_xzzz_0_zzzzzzzz_0, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxzzzz_1, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyzzzz_1, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxzzzzz_1, g_zzz_0_xxzzzzzz_1, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyzzzz_1, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyzzzzz_1, g_zzz_0_xyzzzzzz_1, g_zzz_0_xzzzzzz_1, g_zzz_0_xzzzzzzz_1, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyzzzz_1, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyzzzzz_1, g_zzz_0_yyzzzzzz_1, g_zzz_0_yzzzzzz_1, g_zzz_0_yzzzzzzz_1, g_zzz_0_zzzzzzz_1, g_zzz_0_zzzzzzzz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xxxxxxxx_0[i] = 8.0 * g_zzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxx_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxy_0[i] = 7.0 * g_zzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxxz_0[i] = 7.0 * g_zzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxyy_0[i] = 6.0 * g_zzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxyz_0[i] = 6.0 * g_zzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxxzz_0[i] = 6.0 * g_zzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyyy_0[i] = 5.0 * g_zzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyyz_0[i] = 5.0 * g_zzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxyzz_0[i] = 5.0 * g_zzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxxzzz_0[i] = 5.0 * g_zzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyyz_0[i] = 4.0 * g_zzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyyzz_0[i] = 4.0 * g_zzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxyzzz_0[i] = 4.0 * g_zzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxxzzzz_0[i] = 4.0 * g_zzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyyy_0[i] = 3.0 * g_zzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyyz_0[i] = 3.0 * g_zzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyyzzz_0[i] = 3.0 * g_zzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxyzzzz_0[i] = 3.0 * g_zzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxxzzzzz_0[i] = 3.0 * g_zzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyyy_0[i] = 2.0 * g_zzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyyz_0[i] = 2.0 * g_zzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyyzz_0[i] = 2.0 * g_zzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyyzzz_0[i] = 2.0 * g_zzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxyzzzzz_0[i] = 2.0 * g_zzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xxzzzzzz_0[i] = 2.0 * g_zzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyyy_0[i] = g_zzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyyz_0[i] = g_zzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyyzz_0[i] = g_zzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyyzzz_0[i] = g_zzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyyzzzz_0[i] = g_zzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyyzzzzz_0[i] = g_zzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xyzzzzzz_0[i] = g_zzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_xzzzzzzz_0[i] = g_zzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyyy_0[i] = g_zzz_0_yyyyyyyy_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyyz_0[i] = g_zzz_0_yyyyyyyz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyyzz_0[i] = g_zzz_0_yyyyyyzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyyzzz_0[i] = g_zzz_0_yyyyyzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyyzzzz_0[i] = g_zzz_0_yyyyzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyyzzzzz_0[i] = g_zzz_0_yyyzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yyzzzzzz_0[i] = g_zzz_0_yyzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_yzzzzzzz_0[i] = g_zzz_0_yzzzzzzz_1[i] * wa_x[i];

        g_xzzz_0_zzzzzzzz_0[i] = g_zzz_0_zzzzzzzz_1[i] * wa_x[i];
    }

    /// Set up 450-495 components of targeted buffer : GSL

    auto g_yyyy_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 450);

    auto g_yyyy_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 451);

    auto g_yyyy_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 452);

    auto g_yyyy_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 453);

    auto g_yyyy_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 454);

    auto g_yyyy_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 455);

    auto g_yyyy_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 456);

    auto g_yyyy_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 457);

    auto g_yyyy_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 458);

    auto g_yyyy_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 459);

    auto g_yyyy_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 460);

    auto g_yyyy_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 461);

    auto g_yyyy_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 462);

    auto g_yyyy_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 463);

    auto g_yyyy_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 464);

    auto g_yyyy_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 465);

    auto g_yyyy_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 466);

    auto g_yyyy_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 467);

    auto g_yyyy_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 468);

    auto g_yyyy_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 469);

    auto g_yyyy_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 470);

    auto g_yyyy_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 471);

    auto g_yyyy_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 472);

    auto g_yyyy_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 473);

    auto g_yyyy_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 474);

    auto g_yyyy_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 475);

    auto g_yyyy_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 476);

    auto g_yyyy_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 477);

    auto g_yyyy_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 478);

    auto g_yyyy_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 479);

    auto g_yyyy_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 480);

    auto g_yyyy_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 481);

    auto g_yyyy_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 482);

    auto g_yyyy_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 483);

    auto g_yyyy_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 484);

    auto g_yyyy_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 485);

    auto g_yyyy_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 486);

    auto g_yyyy_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 487);

    auto g_yyyy_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 488);

    auto g_yyyy_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 489);

    auto g_yyyy_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 490);

    auto g_yyyy_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 491);

    auto g_yyyy_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 492);

    auto g_yyyy_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 493);

    auto g_yyyy_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 494);

    #pragma omp simd aligned(g_yy_0_xxxxxxxx_0, g_yy_0_xxxxxxxx_1, g_yy_0_xxxxxxxy_0, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxxz_0, g_yy_0_xxxxxxxz_1, g_yy_0_xxxxxxyy_0, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxxyz_0, g_yy_0_xxxxxxyz_1, g_yy_0_xxxxxxzz_0, g_yy_0_xxxxxxzz_1, g_yy_0_xxxxxyyy_0, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxxyyz_0, g_yy_0_xxxxxyyz_1, g_yy_0_xxxxxyzz_0, g_yy_0_xxxxxyzz_1, g_yy_0_xxxxxzzz_0, g_yy_0_xxxxxzzz_1, g_yy_0_xxxxyyyy_0, g_yy_0_xxxxyyyy_1, g_yy_0_xxxxyyyz_0, g_yy_0_xxxxyyyz_1, g_yy_0_xxxxyyzz_0, g_yy_0_xxxxyyzz_1, g_yy_0_xxxxyzzz_0, g_yy_0_xxxxyzzz_1, g_yy_0_xxxxzzzz_0, g_yy_0_xxxxzzzz_1, g_yy_0_xxxyyyyy_0, g_yy_0_xxxyyyyy_1, g_yy_0_xxxyyyyz_0, g_yy_0_xxxyyyyz_1, g_yy_0_xxxyyyzz_0, g_yy_0_xxxyyyzz_1, g_yy_0_xxxyyzzz_0, g_yy_0_xxxyyzzz_1, g_yy_0_xxxyzzzz_0, g_yy_0_xxxyzzzz_1, g_yy_0_xxxzzzzz_0, g_yy_0_xxxzzzzz_1, g_yy_0_xxyyyyyy_0, g_yy_0_xxyyyyyy_1, g_yy_0_xxyyyyyz_0, g_yy_0_xxyyyyyz_1, g_yy_0_xxyyyyzz_0, g_yy_0_xxyyyyzz_1, g_yy_0_xxyyyzzz_0, g_yy_0_xxyyyzzz_1, g_yy_0_xxyyzzzz_0, g_yy_0_xxyyzzzz_1, g_yy_0_xxyzzzzz_0, g_yy_0_xxyzzzzz_1, g_yy_0_xxzzzzzz_0, g_yy_0_xxzzzzzz_1, g_yy_0_xyyyyyyy_0, g_yy_0_xyyyyyyy_1, g_yy_0_xyyyyyyz_0, g_yy_0_xyyyyyyz_1, g_yy_0_xyyyyyzz_0, g_yy_0_xyyyyyzz_1, g_yy_0_xyyyyzzz_0, g_yy_0_xyyyyzzz_1, g_yy_0_xyyyzzzz_0, g_yy_0_xyyyzzzz_1, g_yy_0_xyyzzzzz_0, g_yy_0_xyyzzzzz_1, g_yy_0_xyzzzzzz_0, g_yy_0_xyzzzzzz_1, g_yy_0_xzzzzzzz_0, g_yy_0_xzzzzzzz_1, g_yy_0_yyyyyyyy_0, g_yy_0_yyyyyyyy_1, g_yy_0_yyyyyyyz_0, g_yy_0_yyyyyyyz_1, g_yy_0_yyyyyyzz_0, g_yy_0_yyyyyyzz_1, g_yy_0_yyyyyzzz_0, g_yy_0_yyyyyzzz_1, g_yy_0_yyyyzzzz_0, g_yy_0_yyyyzzzz_1, g_yy_0_yyyzzzzz_0, g_yy_0_yyyzzzzz_1, g_yy_0_yyzzzzzz_0, g_yy_0_yyzzzzzz_1, g_yy_0_yzzzzzzz_0, g_yy_0_yzzzzzzz_1, g_yy_0_zzzzzzzz_0, g_yy_0_zzzzzzzz_1, g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxzzzz_1, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyzzzz_1, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxzzzzz_1, g_yyy_0_xxzzzzzz_1, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyzzzz_1, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyzzzzz_1, g_yyy_0_xyzzzzzz_1, g_yyy_0_xzzzzzz_1, g_yyy_0_xzzzzzzz_1, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyzzzz_1, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyzzzzz_1, g_yyy_0_yyzzzzzz_1, g_yyy_0_yzzzzzz_1, g_yyy_0_yzzzzzzz_1, g_yyy_0_zzzzzzz_1, g_yyy_0_zzzzzzzz_1, g_yyyy_0_xxxxxxxx_0, g_yyyy_0_xxxxxxxy_0, g_yyyy_0_xxxxxxxz_0, g_yyyy_0_xxxxxxyy_0, g_yyyy_0_xxxxxxyz_0, g_yyyy_0_xxxxxxzz_0, g_yyyy_0_xxxxxyyy_0, g_yyyy_0_xxxxxyyz_0, g_yyyy_0_xxxxxyzz_0, g_yyyy_0_xxxxxzzz_0, g_yyyy_0_xxxxyyyy_0, g_yyyy_0_xxxxyyyz_0, g_yyyy_0_xxxxyyzz_0, g_yyyy_0_xxxxyzzz_0, g_yyyy_0_xxxxzzzz_0, g_yyyy_0_xxxyyyyy_0, g_yyyy_0_xxxyyyyz_0, g_yyyy_0_xxxyyyzz_0, g_yyyy_0_xxxyyzzz_0, g_yyyy_0_xxxyzzzz_0, g_yyyy_0_xxxzzzzz_0, g_yyyy_0_xxyyyyyy_0, g_yyyy_0_xxyyyyyz_0, g_yyyy_0_xxyyyyzz_0, g_yyyy_0_xxyyyzzz_0, g_yyyy_0_xxyyzzzz_0, g_yyyy_0_xxyzzzzz_0, g_yyyy_0_xxzzzzzz_0, g_yyyy_0_xyyyyyyy_0, g_yyyy_0_xyyyyyyz_0, g_yyyy_0_xyyyyyzz_0, g_yyyy_0_xyyyyzzz_0, g_yyyy_0_xyyyzzzz_0, g_yyyy_0_xyyzzzzz_0, g_yyyy_0_xyzzzzzz_0, g_yyyy_0_xzzzzzzz_0, g_yyyy_0_yyyyyyyy_0, g_yyyy_0_yyyyyyyz_0, g_yyyy_0_yyyyyyzz_0, g_yyyy_0_yyyyyzzz_0, g_yyyy_0_yyyyzzzz_0, g_yyyy_0_yyyzzzzz_0, g_yyyy_0_yyzzzzzz_0, g_yyyy_0_yzzzzzzz_0, g_yyyy_0_zzzzzzzz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xxxxxxxx_0[i] = 3.0 * g_yy_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxx_1[i] * fz_be_0 + g_yyy_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxy_0[i] = 3.0 * g_yy_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxxz_0[i] = 3.0 * g_yy_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxxz_1[i] * fz_be_0 + g_yyy_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxyy_0[i] = 3.0 * g_yy_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxyy_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxyz_0[i] = 3.0 * g_yy_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxyz_1[i] * fz_be_0 + g_yyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxxzz_0[i] = 3.0 * g_yy_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxxzz_1[i] * fz_be_0 + g_yyy_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyyy_0[i] = 3.0 * g_yy_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyyy_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyyz_0[i] = 3.0 * g_yy_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxyzz_0[i] = 3.0 * g_yy_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxyzz_1[i] * fz_be_0 + g_yyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxxzzz_0[i] = 3.0 * g_yy_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxxzzz_1[i] * fz_be_0 + g_yyy_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyyy_0[i] = 3.0 * g_yy_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyyy_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyyz_0[i] = 3.0 * g_yy_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyyzz_0[i] = 3.0 * g_yy_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxyzzz_0[i] = 3.0 * g_yy_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxyzzz_1[i] * fz_be_0 + g_yyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxxzzzz_0[i] = 3.0 * g_yy_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxxzzzz_1[i] * fz_be_0 + g_yyy_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyyy_0[i] = 3.0 * g_yy_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyyy_1[i] * fz_be_0 + 5.0 * g_yyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyyz_0[i] = 3.0 * g_yy_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyyzz_0[i] = 3.0 * g_yy_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyyzzz_0[i] = 3.0 * g_yy_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxyzzzz_0[i] = 3.0 * g_yy_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxyzzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxxzzzzz_0[i] = 3.0 * g_yy_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxxzzzzz_1[i] * fz_be_0 + g_yyy_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyyy_0[i] = 3.0 * g_yy_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyyy_1[i] * fz_be_0 + 6.0 * g_yyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyyz_0[i] = 3.0 * g_yy_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyyzz_0[i] = 3.0 * g_yy_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyyzzz_0[i] = 3.0 * g_yy_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyyzzzz_0[i] = 3.0 * g_yy_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxyzzzzz_0[i] = 3.0 * g_yy_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxyzzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xxzzzzzz_0[i] = 3.0 * g_yy_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xxzzzzzz_1[i] * fz_be_0 + g_yyy_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyyy_0[i] = 3.0 * g_yy_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyyy_1[i] * fz_be_0 + 7.0 * g_yyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyyz_0[i] = 3.0 * g_yy_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyyzz_0[i] = 3.0 * g_yy_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyyzzz_0[i] = 3.0 * g_yy_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyyzzzz_0[i] = 3.0 * g_yy_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyyzzzzz_0[i] = 3.0 * g_yy_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xyzzzzzz_0[i] = 3.0 * g_yy_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xyzzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_xzzzzzzz_0[i] = 3.0 * g_yy_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_xzzzzzzz_1[i] * fz_be_0 + g_yyy_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyyy_0[i] = 3.0 * g_yy_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyyy_1[i] * fz_be_0 + 8.0 * g_yyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyy_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyyz_0[i] = 3.0 * g_yy_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyyzz_0[i] = 3.0 * g_yy_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyyzzz_0[i] = 3.0 * g_yy_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyyzzzz_0[i] = 3.0 * g_yy_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyyzzzzz_0[i] = 3.0 * g_yy_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yyzzzzzz_0[i] = 3.0 * g_yy_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_yzzzzzzz_0[i] = 3.0 * g_yy_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_yzzzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyyy_0_zzzzzzzz_0[i] = 3.0 * g_yy_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_yy_0_zzzzzzzz_1[i] * fz_be_0 + g_yyy_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 495-540 components of targeted buffer : GSL

    auto g_yyyz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 495);

    auto g_yyyz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 496);

    auto g_yyyz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 497);

    auto g_yyyz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 498);

    auto g_yyyz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 499);

    auto g_yyyz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 500);

    auto g_yyyz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 501);

    auto g_yyyz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 502);

    auto g_yyyz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 503);

    auto g_yyyz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 504);

    auto g_yyyz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 505);

    auto g_yyyz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 506);

    auto g_yyyz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 507);

    auto g_yyyz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 508);

    auto g_yyyz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 509);

    auto g_yyyz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 510);

    auto g_yyyz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 511);

    auto g_yyyz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 512);

    auto g_yyyz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 513);

    auto g_yyyz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 514);

    auto g_yyyz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 515);

    auto g_yyyz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 516);

    auto g_yyyz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 517);

    auto g_yyyz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 518);

    auto g_yyyz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 519);

    auto g_yyyz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 520);

    auto g_yyyz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 521);

    auto g_yyyz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 522);

    auto g_yyyz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 523);

    auto g_yyyz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 524);

    auto g_yyyz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 525);

    auto g_yyyz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 526);

    auto g_yyyz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 527);

    auto g_yyyz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 528);

    auto g_yyyz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 529);

    auto g_yyyz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 530);

    auto g_yyyz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 531);

    auto g_yyyz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 532);

    auto g_yyyz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 533);

    auto g_yyyz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 534);

    auto g_yyyz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 535);

    auto g_yyyz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 536);

    auto g_yyyz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 537);

    auto g_yyyz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 538);

    auto g_yyyz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 539);

    #pragma omp simd aligned(g_yyy_0_xxxxxxx_1, g_yyy_0_xxxxxxxx_1, g_yyy_0_xxxxxxxy_1, g_yyy_0_xxxxxxxz_1, g_yyy_0_xxxxxxy_1, g_yyy_0_xxxxxxyy_1, g_yyy_0_xxxxxxyz_1, g_yyy_0_xxxxxxz_1, g_yyy_0_xxxxxxzz_1, g_yyy_0_xxxxxyy_1, g_yyy_0_xxxxxyyy_1, g_yyy_0_xxxxxyyz_1, g_yyy_0_xxxxxyz_1, g_yyy_0_xxxxxyzz_1, g_yyy_0_xxxxxzz_1, g_yyy_0_xxxxxzzz_1, g_yyy_0_xxxxyyy_1, g_yyy_0_xxxxyyyy_1, g_yyy_0_xxxxyyyz_1, g_yyy_0_xxxxyyz_1, g_yyy_0_xxxxyyzz_1, g_yyy_0_xxxxyzz_1, g_yyy_0_xxxxyzzz_1, g_yyy_0_xxxxzzz_1, g_yyy_0_xxxxzzzz_1, g_yyy_0_xxxyyyy_1, g_yyy_0_xxxyyyyy_1, g_yyy_0_xxxyyyyz_1, g_yyy_0_xxxyyyz_1, g_yyy_0_xxxyyyzz_1, g_yyy_0_xxxyyzz_1, g_yyy_0_xxxyyzzz_1, g_yyy_0_xxxyzzz_1, g_yyy_0_xxxyzzzz_1, g_yyy_0_xxxzzzz_1, g_yyy_0_xxxzzzzz_1, g_yyy_0_xxyyyyy_1, g_yyy_0_xxyyyyyy_1, g_yyy_0_xxyyyyyz_1, g_yyy_0_xxyyyyz_1, g_yyy_0_xxyyyyzz_1, g_yyy_0_xxyyyzz_1, g_yyy_0_xxyyyzzz_1, g_yyy_0_xxyyzzz_1, g_yyy_0_xxyyzzzz_1, g_yyy_0_xxyzzzz_1, g_yyy_0_xxyzzzzz_1, g_yyy_0_xxzzzzz_1, g_yyy_0_xxzzzzzz_1, g_yyy_0_xyyyyyy_1, g_yyy_0_xyyyyyyy_1, g_yyy_0_xyyyyyyz_1, g_yyy_0_xyyyyyz_1, g_yyy_0_xyyyyyzz_1, g_yyy_0_xyyyyzz_1, g_yyy_0_xyyyyzzz_1, g_yyy_0_xyyyzzz_1, g_yyy_0_xyyyzzzz_1, g_yyy_0_xyyzzzz_1, g_yyy_0_xyyzzzzz_1, g_yyy_0_xyzzzzz_1, g_yyy_0_xyzzzzzz_1, g_yyy_0_xzzzzzz_1, g_yyy_0_xzzzzzzz_1, g_yyy_0_yyyyyyy_1, g_yyy_0_yyyyyyyy_1, g_yyy_0_yyyyyyyz_1, g_yyy_0_yyyyyyz_1, g_yyy_0_yyyyyyzz_1, g_yyy_0_yyyyyzz_1, g_yyy_0_yyyyyzzz_1, g_yyy_0_yyyyzzz_1, g_yyy_0_yyyyzzzz_1, g_yyy_0_yyyzzzz_1, g_yyy_0_yyyzzzzz_1, g_yyy_0_yyzzzzz_1, g_yyy_0_yyzzzzzz_1, g_yyy_0_yzzzzzz_1, g_yyy_0_yzzzzzzz_1, g_yyy_0_zzzzzzz_1, g_yyy_0_zzzzzzzz_1, g_yyyz_0_xxxxxxxx_0, g_yyyz_0_xxxxxxxy_0, g_yyyz_0_xxxxxxxz_0, g_yyyz_0_xxxxxxyy_0, g_yyyz_0_xxxxxxyz_0, g_yyyz_0_xxxxxxzz_0, g_yyyz_0_xxxxxyyy_0, g_yyyz_0_xxxxxyyz_0, g_yyyz_0_xxxxxyzz_0, g_yyyz_0_xxxxxzzz_0, g_yyyz_0_xxxxyyyy_0, g_yyyz_0_xxxxyyyz_0, g_yyyz_0_xxxxyyzz_0, g_yyyz_0_xxxxyzzz_0, g_yyyz_0_xxxxzzzz_0, g_yyyz_0_xxxyyyyy_0, g_yyyz_0_xxxyyyyz_0, g_yyyz_0_xxxyyyzz_0, g_yyyz_0_xxxyyzzz_0, g_yyyz_0_xxxyzzzz_0, g_yyyz_0_xxxzzzzz_0, g_yyyz_0_xxyyyyyy_0, g_yyyz_0_xxyyyyyz_0, g_yyyz_0_xxyyyyzz_0, g_yyyz_0_xxyyyzzz_0, g_yyyz_0_xxyyzzzz_0, g_yyyz_0_xxyzzzzz_0, g_yyyz_0_xxzzzzzz_0, g_yyyz_0_xyyyyyyy_0, g_yyyz_0_xyyyyyyz_0, g_yyyz_0_xyyyyyzz_0, g_yyyz_0_xyyyyzzz_0, g_yyyz_0_xyyyzzzz_0, g_yyyz_0_xyyzzzzz_0, g_yyyz_0_xyzzzzzz_0, g_yyyz_0_xzzzzzzz_0, g_yyyz_0_yyyyyyyy_0, g_yyyz_0_yyyyyyyz_0, g_yyyz_0_yyyyyyzz_0, g_yyyz_0_yyyyyzzz_0, g_yyyz_0_yyyyzzzz_0, g_yyyz_0_yyyzzzzz_0, g_yyyz_0_yyzzzzzz_0, g_yyyz_0_yzzzzzzz_0, g_yyyz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xxxxxxxx_0[i] = g_yyy_0_xxxxxxxx_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxy_0[i] = g_yyy_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxxz_0[i] = g_yyy_0_xxxxxxx_1[i] * fi_acd_0 + g_yyy_0_xxxxxxxz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxyy_0[i] = g_yyy_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxyz_0[i] = g_yyy_0_xxxxxxy_1[i] * fi_acd_0 + g_yyy_0_xxxxxxyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxxzz_0[i] = 2.0 * g_yyy_0_xxxxxxz_1[i] * fi_acd_0 + g_yyy_0_xxxxxxzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyyy_0[i] = g_yyy_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyyz_0[i] = g_yyy_0_xxxxxyy_1[i] * fi_acd_0 + g_yyy_0_xxxxxyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxyzz_0[i] = 2.0 * g_yyy_0_xxxxxyz_1[i] * fi_acd_0 + g_yyy_0_xxxxxyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxxzzz_0[i] = 3.0 * g_yyy_0_xxxxxzz_1[i] * fi_acd_0 + g_yyy_0_xxxxxzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyyy_0[i] = g_yyy_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyyz_0[i] = g_yyy_0_xxxxyyy_1[i] * fi_acd_0 + g_yyy_0_xxxxyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyyzz_0[i] = 2.0 * g_yyy_0_xxxxyyz_1[i] * fi_acd_0 + g_yyy_0_xxxxyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxyzzz_0[i] = 3.0 * g_yyy_0_xxxxyzz_1[i] * fi_acd_0 + g_yyy_0_xxxxyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxxzzzz_0[i] = 4.0 * g_yyy_0_xxxxzzz_1[i] * fi_acd_0 + g_yyy_0_xxxxzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyyy_0[i] = g_yyy_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyyz_0[i] = g_yyy_0_xxxyyyy_1[i] * fi_acd_0 + g_yyy_0_xxxyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyyzz_0[i] = 2.0 * g_yyy_0_xxxyyyz_1[i] * fi_acd_0 + g_yyy_0_xxxyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyyzzz_0[i] = 3.0 * g_yyy_0_xxxyyzz_1[i] * fi_acd_0 + g_yyy_0_xxxyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxyzzzz_0[i] = 4.0 * g_yyy_0_xxxyzzz_1[i] * fi_acd_0 + g_yyy_0_xxxyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxxzzzzz_0[i] = 5.0 * g_yyy_0_xxxzzzz_1[i] * fi_acd_0 + g_yyy_0_xxxzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyyy_0[i] = g_yyy_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyyz_0[i] = g_yyy_0_xxyyyyy_1[i] * fi_acd_0 + g_yyy_0_xxyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyyzz_0[i] = 2.0 * g_yyy_0_xxyyyyz_1[i] * fi_acd_0 + g_yyy_0_xxyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyyzzz_0[i] = 3.0 * g_yyy_0_xxyyyzz_1[i] * fi_acd_0 + g_yyy_0_xxyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyyzzzz_0[i] = 4.0 * g_yyy_0_xxyyzzz_1[i] * fi_acd_0 + g_yyy_0_xxyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxyzzzzz_0[i] = 5.0 * g_yyy_0_xxyzzzz_1[i] * fi_acd_0 + g_yyy_0_xxyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xxzzzzzz_0[i] = 6.0 * g_yyy_0_xxzzzzz_1[i] * fi_acd_0 + g_yyy_0_xxzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyyy_0[i] = g_yyy_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyyz_0[i] = g_yyy_0_xyyyyyy_1[i] * fi_acd_0 + g_yyy_0_xyyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyyzz_0[i] = 2.0 * g_yyy_0_xyyyyyz_1[i] * fi_acd_0 + g_yyy_0_xyyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyyzzz_0[i] = 3.0 * g_yyy_0_xyyyyzz_1[i] * fi_acd_0 + g_yyy_0_xyyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyyzzzz_0[i] = 4.0 * g_yyy_0_xyyyzzz_1[i] * fi_acd_0 + g_yyy_0_xyyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyyzzzzz_0[i] = 5.0 * g_yyy_0_xyyzzzz_1[i] * fi_acd_0 + g_yyy_0_xyyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xyzzzzzz_0[i] = 6.0 * g_yyy_0_xyzzzzz_1[i] * fi_acd_0 + g_yyy_0_xyzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_xzzzzzzz_0[i] = 7.0 * g_yyy_0_xzzzzzz_1[i] * fi_acd_0 + g_yyy_0_xzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyyy_0[i] = g_yyy_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyyz_0[i] = g_yyy_0_yyyyyyy_1[i] * fi_acd_0 + g_yyy_0_yyyyyyyz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyyzz_0[i] = 2.0 * g_yyy_0_yyyyyyz_1[i] * fi_acd_0 + g_yyy_0_yyyyyyzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyyzzz_0[i] = 3.0 * g_yyy_0_yyyyyzz_1[i] * fi_acd_0 + g_yyy_0_yyyyyzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyyzzzz_0[i] = 4.0 * g_yyy_0_yyyyzzz_1[i] * fi_acd_0 + g_yyy_0_yyyyzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyyzzzzz_0[i] = 5.0 * g_yyy_0_yyyzzzz_1[i] * fi_acd_0 + g_yyy_0_yyyzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yyzzzzzz_0[i] = 6.0 * g_yyy_0_yyzzzzz_1[i] * fi_acd_0 + g_yyy_0_yyzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_yzzzzzzz_0[i] = 7.0 * g_yyy_0_yzzzzzz_1[i] * fi_acd_0 + g_yyy_0_yzzzzzzz_1[i] * wa_z[i];

        g_yyyz_0_zzzzzzzz_0[i] = 8.0 * g_yyy_0_zzzzzzz_1[i] * fi_acd_0 + g_yyy_0_zzzzzzzz_1[i] * wa_z[i];
    }

    /// Set up 540-585 components of targeted buffer : GSL

    auto g_yyzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 540);

    auto g_yyzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 541);

    auto g_yyzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 542);

    auto g_yyzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 543);

    auto g_yyzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 544);

    auto g_yyzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 545);

    auto g_yyzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 546);

    auto g_yyzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 547);

    auto g_yyzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 548);

    auto g_yyzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 549);

    auto g_yyzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 550);

    auto g_yyzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 551);

    auto g_yyzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 552);

    auto g_yyzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 553);

    auto g_yyzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 554);

    auto g_yyzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 555);

    auto g_yyzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 556);

    auto g_yyzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 557);

    auto g_yyzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 558);

    auto g_yyzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 559);

    auto g_yyzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 560);

    auto g_yyzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 561);

    auto g_yyzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 562);

    auto g_yyzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 563);

    auto g_yyzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 564);

    auto g_yyzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 565);

    auto g_yyzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 566);

    auto g_yyzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 567);

    auto g_yyzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 568);

    auto g_yyzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 569);

    auto g_yyzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 570);

    auto g_yyzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 571);

    auto g_yyzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 572);

    auto g_yyzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 573);

    auto g_yyzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 574);

    auto g_yyzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 575);

    auto g_yyzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 576);

    auto g_yyzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 577);

    auto g_yyzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 578);

    auto g_yyzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 579);

    auto g_yyzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 580);

    auto g_yyzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 581);

    auto g_yyzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 582);

    auto g_yyzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 583);

    auto g_yyzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 584);

    #pragma omp simd aligned(g_yy_0_xxxxxxxy_0, g_yy_0_xxxxxxxy_1, g_yy_0_xxxxxxyy_0, g_yy_0_xxxxxxyy_1, g_yy_0_xxxxxyyy_0, g_yy_0_xxxxxyyy_1, g_yy_0_xxxxyyyy_0, g_yy_0_xxxxyyyy_1, g_yy_0_xxxyyyyy_0, g_yy_0_xxxyyyyy_1, g_yy_0_xxyyyyyy_0, g_yy_0_xxyyyyyy_1, g_yy_0_xyyyyyyy_0, g_yy_0_xyyyyyyy_1, g_yy_0_yyyyyyyy_0, g_yy_0_yyyyyyyy_1, g_yyz_0_xxxxxxxy_1, g_yyz_0_xxxxxxyy_1, g_yyz_0_xxxxxyyy_1, g_yyz_0_xxxxyyyy_1, g_yyz_0_xxxyyyyy_1, g_yyz_0_xxyyyyyy_1, g_yyz_0_xyyyyyyy_1, g_yyz_0_yyyyyyyy_1, g_yyzz_0_xxxxxxxx_0, g_yyzz_0_xxxxxxxy_0, g_yyzz_0_xxxxxxxz_0, g_yyzz_0_xxxxxxyy_0, g_yyzz_0_xxxxxxyz_0, g_yyzz_0_xxxxxxzz_0, g_yyzz_0_xxxxxyyy_0, g_yyzz_0_xxxxxyyz_0, g_yyzz_0_xxxxxyzz_0, g_yyzz_0_xxxxxzzz_0, g_yyzz_0_xxxxyyyy_0, g_yyzz_0_xxxxyyyz_0, g_yyzz_0_xxxxyyzz_0, g_yyzz_0_xxxxyzzz_0, g_yyzz_0_xxxxzzzz_0, g_yyzz_0_xxxyyyyy_0, g_yyzz_0_xxxyyyyz_0, g_yyzz_0_xxxyyyzz_0, g_yyzz_0_xxxyyzzz_0, g_yyzz_0_xxxyzzzz_0, g_yyzz_0_xxxzzzzz_0, g_yyzz_0_xxyyyyyy_0, g_yyzz_0_xxyyyyyz_0, g_yyzz_0_xxyyyyzz_0, g_yyzz_0_xxyyyzzz_0, g_yyzz_0_xxyyzzzz_0, g_yyzz_0_xxyzzzzz_0, g_yyzz_0_xxzzzzzz_0, g_yyzz_0_xyyyyyyy_0, g_yyzz_0_xyyyyyyz_0, g_yyzz_0_xyyyyyzz_0, g_yyzz_0_xyyyyzzz_0, g_yyzz_0_xyyyzzzz_0, g_yyzz_0_xyyzzzzz_0, g_yyzz_0_xyzzzzzz_0, g_yyzz_0_xzzzzzzz_0, g_yyzz_0_yyyyyyyy_0, g_yyzz_0_yyyyyyyz_0, g_yyzz_0_yyyyyyzz_0, g_yyzz_0_yyyyyzzz_0, g_yyzz_0_yyyyzzzz_0, g_yyzz_0_yyyzzzzz_0, g_yyzz_0_yyzzzzzz_0, g_yyzz_0_yzzzzzzz_0, g_yyzz_0_zzzzzzzz_0, g_yzz_0_xxxxxxxx_1, g_yzz_0_xxxxxxxz_1, g_yzz_0_xxxxxxyz_1, g_yzz_0_xxxxxxz_1, g_yzz_0_xxxxxxzz_1, g_yzz_0_xxxxxyyz_1, g_yzz_0_xxxxxyz_1, g_yzz_0_xxxxxyzz_1, g_yzz_0_xxxxxzz_1, g_yzz_0_xxxxxzzz_1, g_yzz_0_xxxxyyyz_1, g_yzz_0_xxxxyyz_1, g_yzz_0_xxxxyyzz_1, g_yzz_0_xxxxyzz_1, g_yzz_0_xxxxyzzz_1, g_yzz_0_xxxxzzz_1, g_yzz_0_xxxxzzzz_1, g_yzz_0_xxxyyyyz_1, g_yzz_0_xxxyyyz_1, g_yzz_0_xxxyyyzz_1, g_yzz_0_xxxyyzz_1, g_yzz_0_xxxyyzzz_1, g_yzz_0_xxxyzzz_1, g_yzz_0_xxxyzzzz_1, g_yzz_0_xxxzzzz_1, g_yzz_0_xxxzzzzz_1, g_yzz_0_xxyyyyyz_1, g_yzz_0_xxyyyyz_1, g_yzz_0_xxyyyyzz_1, g_yzz_0_xxyyyzz_1, g_yzz_0_xxyyyzzz_1, g_yzz_0_xxyyzzz_1, g_yzz_0_xxyyzzzz_1, g_yzz_0_xxyzzzz_1, g_yzz_0_xxyzzzzz_1, g_yzz_0_xxzzzzz_1, g_yzz_0_xxzzzzzz_1, g_yzz_0_xyyyyyyz_1, g_yzz_0_xyyyyyz_1, g_yzz_0_xyyyyyzz_1, g_yzz_0_xyyyyzz_1, g_yzz_0_xyyyyzzz_1, g_yzz_0_xyyyzzz_1, g_yzz_0_xyyyzzzz_1, g_yzz_0_xyyzzzz_1, g_yzz_0_xyyzzzzz_1, g_yzz_0_xyzzzzz_1, g_yzz_0_xyzzzzzz_1, g_yzz_0_xzzzzzz_1, g_yzz_0_xzzzzzzz_1, g_yzz_0_yyyyyyyz_1, g_yzz_0_yyyyyyz_1, g_yzz_0_yyyyyyzz_1, g_yzz_0_yyyyyzz_1, g_yzz_0_yyyyyzzz_1, g_yzz_0_yyyyzzz_1, g_yzz_0_yyyyzzzz_1, g_yzz_0_yyyzzzz_1, g_yzz_0_yyyzzzzz_1, g_yzz_0_yyzzzzz_1, g_yzz_0_yyzzzzzz_1, g_yzz_0_yzzzzzz_1, g_yzz_0_yzzzzzzz_1, g_yzz_0_zzzzzzz_1, g_yzz_0_zzzzzzzz_1, g_zz_0_xxxxxxxx_0, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxz_0, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxyz_0, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxzz_0, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxyyz_0, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyzz_0, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxzzz_0, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxyyyz_0, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyzz_0, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyzzz_0, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxzzzz_0, g_zz_0_xxxxzzzz_1, g_zz_0_xxxyyyyz_0, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyzz_0, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyzzz_0, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyzzzz_0, g_zz_0_xxxyzzzz_1, g_zz_0_xxxzzzzz_0, g_zz_0_xxxzzzzz_1, g_zz_0_xxyyyyyz_0, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyzz_0, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyzzz_0, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyzzzz_0, g_zz_0_xxyyzzzz_1, g_zz_0_xxyzzzzz_0, g_zz_0_xxyzzzzz_1, g_zz_0_xxzzzzzz_0, g_zz_0_xxzzzzzz_1, g_zz_0_xyyyyyyz_0, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyzz_0, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyzzz_0, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyzzzz_0, g_zz_0_xyyyzzzz_1, g_zz_0_xyyzzzzz_0, g_zz_0_xyyzzzzz_1, g_zz_0_xyzzzzzz_0, g_zz_0_xyzzzzzz_1, g_zz_0_xzzzzzzz_0, g_zz_0_xzzzzzzz_1, g_zz_0_yyyyyyyz_0, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyzz_0, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyzzz_0, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyzzzz_0, g_zz_0_yyyyzzzz_1, g_zz_0_yyyzzzzz_0, g_zz_0_yyyzzzzz_1, g_zz_0_yyzzzzzz_0, g_zz_0_yyzzzzzz_1, g_zz_0_yzzzzzzz_0, g_zz_0_yzzzzzzz_1, g_zz_0_zzzzzzzz_0, g_zz_0_zzzzzzzz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xxxxxxxx_0[i] = g_zz_0_xxxxxxxx_0[i] * fbe_0 - g_zz_0_xxxxxxxx_1[i] * fz_be_0 + g_yzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxxy_0[i] = g_yy_0_xxxxxxxy_0[i] * fbe_0 - g_yy_0_xxxxxxxy_1[i] * fz_be_0 + g_yyz_0_xxxxxxxy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxxxz_0[i] = g_zz_0_xxxxxxxz_0[i] * fbe_0 - g_zz_0_xxxxxxxz_1[i] * fz_be_0 + g_yzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxyy_0[i] = g_yy_0_xxxxxxyy_0[i] * fbe_0 - g_yy_0_xxxxxxyy_1[i] * fz_be_0 + g_yyz_0_xxxxxxyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxxyz_0[i] = g_zz_0_xxxxxxyz_0[i] * fbe_0 - g_zz_0_xxxxxxyz_1[i] * fz_be_0 + g_yzz_0_xxxxxxz_1[i] * fi_acd_0 + g_yzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxxzz_0[i] = g_zz_0_xxxxxxzz_0[i] * fbe_0 - g_zz_0_xxxxxxzz_1[i] * fz_be_0 + g_yzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxyyy_0[i] = g_yy_0_xxxxxyyy_0[i] * fbe_0 - g_yy_0_xxxxxyyy_1[i] * fz_be_0 + g_yyz_0_xxxxxyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxxyyz_0[i] = g_zz_0_xxxxxyyz_0[i] * fbe_0 - g_zz_0_xxxxxyyz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxxxyz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxyzz_0[i] = g_zz_0_xxxxxyzz_0[i] * fbe_0 - g_zz_0_xxxxxyzz_1[i] * fz_be_0 + g_yzz_0_xxxxxzz_1[i] * fi_acd_0 + g_yzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxxzzz_0[i] = g_zz_0_xxxxxzzz_0[i] * fbe_0 - g_zz_0_xxxxxzzz_1[i] * fz_be_0 + g_yzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyyyy_0[i] = g_yy_0_xxxxyyyy_0[i] * fbe_0 - g_yy_0_xxxxyyyy_1[i] * fz_be_0 + g_yyz_0_xxxxyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxxyyyz_0[i] = g_zz_0_xxxxyyyz_0[i] * fbe_0 - g_zz_0_xxxxyyyz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxxxyyz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyyzz_0[i] = g_zz_0_xxxxyyzz_0[i] * fbe_0 - g_zz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxxyzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxyzzz_0[i] = g_zz_0_xxxxyzzz_0[i] * fbe_0 - g_zz_0_xxxxyzzz_1[i] * fz_be_0 + g_yzz_0_xxxxzzz_1[i] * fi_acd_0 + g_yzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxxzzzz_0[i] = g_zz_0_xxxxzzzz_0[i] * fbe_0 - g_zz_0_xxxxzzzz_1[i] * fz_be_0 + g_yzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyyyy_0[i] = g_yy_0_xxxyyyyy_0[i] * fbe_0 - g_yy_0_xxxyyyyy_1[i] * fz_be_0 + g_yyz_0_xxxyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxxyyyyz_0[i] = g_zz_0_xxxyyyyz_0[i] * fbe_0 - g_zz_0_xxxyyyyz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xxxyyyz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyyzz_0[i] = g_zz_0_xxxyyyzz_0[i] * fbe_0 - g_zz_0_xxxyyyzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxxyyzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyyzzz_0[i] = g_zz_0_xxxyyzzz_0[i] * fbe_0 - g_zz_0_xxxyyzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxxyzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxyzzzz_0[i] = g_zz_0_xxxyzzzz_0[i] * fbe_0 - g_zz_0_xxxyzzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzzz_1[i] * fi_acd_0 + g_yzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxxzzzzz_0[i] = g_zz_0_xxxzzzzz_0[i] * fbe_0 - g_zz_0_xxxzzzzz_1[i] * fz_be_0 + g_yzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyyyy_0[i] = g_yy_0_xxyyyyyy_0[i] * fbe_0 - g_yy_0_xxyyyyyy_1[i] * fz_be_0 + g_yyz_0_xxyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xxyyyyyz_0[i] = g_zz_0_xxyyyyyz_0[i] * fbe_0 - g_zz_0_xxyyyyyz_1[i] * fz_be_0 + 5.0 * g_yzz_0_xxyyyyz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyyzz_0[i] = g_zz_0_xxyyyyzz_0[i] * fbe_0 - g_zz_0_xxyyyyzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xxyyyzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyyzzz_0[i] = g_zz_0_xxyyyzzz_0[i] * fbe_0 - g_zz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xxyyzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyyzzzz_0[i] = g_zz_0_xxyyzzzz_0[i] * fbe_0 - g_zz_0_xxyyzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xxyzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxyzzzzz_0[i] = g_zz_0_xxyzzzzz_0[i] * fbe_0 - g_zz_0_xxyzzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzzz_1[i] * fi_acd_0 + g_yzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xxzzzzzz_0[i] = g_zz_0_xxzzzzzz_0[i] * fbe_0 - g_zz_0_xxzzzzzz_1[i] * fz_be_0 + g_yzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyyyy_0[i] = g_yy_0_xyyyyyyy_0[i] * fbe_0 - g_yy_0_xyyyyyyy_1[i] * fz_be_0 + g_yyz_0_xyyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_xyyyyyyz_0[i] = g_zz_0_xyyyyyyz_0[i] * fbe_0 - g_zz_0_xyyyyyyz_1[i] * fz_be_0 + 6.0 * g_yzz_0_xyyyyyz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyyzz_0[i] = g_zz_0_xyyyyyzz_0[i] * fbe_0 - g_zz_0_xyyyyyzz_1[i] * fz_be_0 + 5.0 * g_yzz_0_xyyyyzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyyzzz_0[i] = g_zz_0_xyyyyzzz_0[i] * fbe_0 - g_zz_0_xyyyyzzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_xyyyzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyyzzzz_0[i] = g_zz_0_xyyyzzzz_0[i] * fbe_0 - g_zz_0_xyyyzzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_xyyzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyyzzzzz_0[i] = g_zz_0_xyyzzzzz_0[i] * fbe_0 - g_zz_0_xyyzzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_xyzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xyzzzzzz_0[i] = g_zz_0_xyzzzzzz_0[i] * fbe_0 - g_zz_0_xyzzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzzz_1[i] * fi_acd_0 + g_yzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_xzzzzzzz_0[i] = g_zz_0_xzzzzzzz_0[i] * fbe_0 - g_zz_0_xzzzzzzz_1[i] * fz_be_0 + g_yzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyyyy_0[i] = g_yy_0_yyyyyyyy_0[i] * fbe_0 - g_yy_0_yyyyyyyy_1[i] * fz_be_0 + g_yyz_0_yyyyyyyy_1[i] * wa_z[i];

        g_yyzz_0_yyyyyyyz_0[i] = g_zz_0_yyyyyyyz_0[i] * fbe_0 - g_zz_0_yyyyyyyz_1[i] * fz_be_0 + 7.0 * g_yzz_0_yyyyyyz_1[i] * fi_acd_0 + g_yzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyyzz_0[i] = g_zz_0_yyyyyyzz_0[i] * fbe_0 - g_zz_0_yyyyyyzz_1[i] * fz_be_0 + 6.0 * g_yzz_0_yyyyyzz_1[i] * fi_acd_0 + g_yzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyyzzz_0[i] = g_zz_0_yyyyyzzz_0[i] * fbe_0 - g_zz_0_yyyyyzzz_1[i] * fz_be_0 + 5.0 * g_yzz_0_yyyyzzz_1[i] * fi_acd_0 + g_yzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyyzzzz_0[i] = g_zz_0_yyyyzzzz_0[i] * fbe_0 - g_zz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_yzz_0_yyyzzzz_1[i] * fi_acd_0 + g_yzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyyzzzzz_0[i] = g_zz_0_yyyzzzzz_0[i] * fbe_0 - g_zz_0_yyyzzzzz_1[i] * fz_be_0 + 3.0 * g_yzz_0_yyzzzzz_1[i] * fi_acd_0 + g_yzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yyzzzzzz_0[i] = g_zz_0_yyzzzzzz_0[i] * fbe_0 - g_zz_0_yyzzzzzz_1[i] * fz_be_0 + 2.0 * g_yzz_0_yzzzzzz_1[i] * fi_acd_0 + g_yzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_yzzzzzzz_0[i] = g_zz_0_yzzzzzzz_0[i] * fbe_0 - g_zz_0_yzzzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzzz_1[i] * fi_acd_0 + g_yzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yyzz_0_zzzzzzzz_0[i] = g_zz_0_zzzzzzzz_0[i] * fbe_0 - g_zz_0_zzzzzzzz_1[i] * fz_be_0 + g_yzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 585-630 components of targeted buffer : GSL

    auto g_yzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 585);

    auto g_yzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 586);

    auto g_yzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 587);

    auto g_yzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 588);

    auto g_yzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 589);

    auto g_yzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 590);

    auto g_yzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 591);

    auto g_yzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 592);

    auto g_yzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 593);

    auto g_yzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 594);

    auto g_yzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 595);

    auto g_yzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 596);

    auto g_yzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 597);

    auto g_yzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 598);

    auto g_yzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 599);

    auto g_yzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 600);

    auto g_yzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 601);

    auto g_yzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 602);

    auto g_yzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 603);

    auto g_yzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 604);

    auto g_yzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 605);

    auto g_yzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 606);

    auto g_yzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 607);

    auto g_yzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 608);

    auto g_yzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 609);

    auto g_yzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 610);

    auto g_yzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 611);

    auto g_yzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 612);

    auto g_yzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 613);

    auto g_yzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 614);

    auto g_yzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 615);

    auto g_yzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 616);

    auto g_yzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 617);

    auto g_yzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 618);

    auto g_yzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 619);

    auto g_yzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 620);

    auto g_yzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 621);

    auto g_yzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 622);

    auto g_yzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 623);

    auto g_yzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 624);

    auto g_yzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 625);

    auto g_yzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 626);

    auto g_yzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 627);

    auto g_yzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 628);

    auto g_yzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 629);

    #pragma omp simd aligned(g_yzzz_0_xxxxxxxx_0, g_yzzz_0_xxxxxxxy_0, g_yzzz_0_xxxxxxxz_0, g_yzzz_0_xxxxxxyy_0, g_yzzz_0_xxxxxxyz_0, g_yzzz_0_xxxxxxzz_0, g_yzzz_0_xxxxxyyy_0, g_yzzz_0_xxxxxyyz_0, g_yzzz_0_xxxxxyzz_0, g_yzzz_0_xxxxxzzz_0, g_yzzz_0_xxxxyyyy_0, g_yzzz_0_xxxxyyyz_0, g_yzzz_0_xxxxyyzz_0, g_yzzz_0_xxxxyzzz_0, g_yzzz_0_xxxxzzzz_0, g_yzzz_0_xxxyyyyy_0, g_yzzz_0_xxxyyyyz_0, g_yzzz_0_xxxyyyzz_0, g_yzzz_0_xxxyyzzz_0, g_yzzz_0_xxxyzzzz_0, g_yzzz_0_xxxzzzzz_0, g_yzzz_0_xxyyyyyy_0, g_yzzz_0_xxyyyyyz_0, g_yzzz_0_xxyyyyzz_0, g_yzzz_0_xxyyyzzz_0, g_yzzz_0_xxyyzzzz_0, g_yzzz_0_xxyzzzzz_0, g_yzzz_0_xxzzzzzz_0, g_yzzz_0_xyyyyyyy_0, g_yzzz_0_xyyyyyyz_0, g_yzzz_0_xyyyyyzz_0, g_yzzz_0_xyyyyzzz_0, g_yzzz_0_xyyyzzzz_0, g_yzzz_0_xyyzzzzz_0, g_yzzz_0_xyzzzzzz_0, g_yzzz_0_xzzzzzzz_0, g_yzzz_0_yyyyyyyy_0, g_yzzz_0_yyyyyyyz_0, g_yzzz_0_yyyyyyzz_0, g_yzzz_0_yyyyyzzz_0, g_yzzz_0_yyyyzzzz_0, g_yzzz_0_yyyzzzzz_0, g_yzzz_0_yyzzzzzz_0, g_yzzz_0_yzzzzzzz_0, g_yzzz_0_zzzzzzzz_0, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxzzzz_1, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyzzzz_1, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxzzzzz_1, g_zzz_0_xxzzzzzz_1, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyzzzz_1, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyzzzzz_1, g_zzz_0_xyzzzzzz_1, g_zzz_0_xzzzzzz_1, g_zzz_0_xzzzzzzz_1, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyzzzz_1, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyzzzzz_1, g_zzz_0_yyzzzzzz_1, g_zzz_0_yzzzzzz_1, g_zzz_0_yzzzzzzz_1, g_zzz_0_zzzzzzz_1, g_zzz_0_zzzzzzzz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xxxxxxxx_0[i] = g_zzz_0_xxxxxxxx_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxy_0[i] = g_zzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxxz_0[i] = g_zzz_0_xxxxxxxz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxyy_0[i] = 2.0 * g_zzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxyz_0[i] = g_zzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxxzz_0[i] = g_zzz_0_xxxxxxzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyyy_0[i] = 3.0 * g_zzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyyz_0[i] = 2.0 * g_zzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxyzz_0[i] = g_zzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxxzzz_0[i] = g_zzz_0_xxxxxzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyyy_0[i] = 4.0 * g_zzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyyz_0[i] = 3.0 * g_zzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyyzz_0[i] = 2.0 * g_zzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxyzzz_0[i] = g_zzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxxzzzz_0[i] = g_zzz_0_xxxxzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyyy_0[i] = 5.0 * g_zzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyyz_0[i] = 4.0 * g_zzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyyzz_0[i] = 3.0 * g_zzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyyzzz_0[i] = 2.0 * g_zzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxyzzzz_0[i] = g_zzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxxzzzzz_0[i] = g_zzz_0_xxxzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyyy_0[i] = 6.0 * g_zzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyyz_0[i] = 5.0 * g_zzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyyzz_0[i] = 4.0 * g_zzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyyzzz_0[i] = 3.0 * g_zzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyyzzzz_0[i] = 2.0 * g_zzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxyzzzzz_0[i] = g_zzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xxzzzzzz_0[i] = g_zzz_0_xxzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyyy_0[i] = 7.0 * g_zzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyyz_0[i] = 6.0 * g_zzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyyzz_0[i] = 5.0 * g_zzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyyzzz_0[i] = 4.0 * g_zzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyyzzzz_0[i] = 3.0 * g_zzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyyzzzzz_0[i] = 2.0 * g_zzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xyzzzzzz_0[i] = g_zzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_xzzzzzzz_0[i] = g_zzz_0_xzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyyy_0[i] = 8.0 * g_zzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyy_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyyz_0[i] = 7.0 * g_zzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyyzz_0[i] = 6.0 * g_zzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyyzzz_0[i] = 5.0 * g_zzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyyzzzz_0[i] = 4.0 * g_zzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyyzzzzz_0[i] = 3.0 * g_zzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yyzzzzzz_0[i] = 2.0 * g_zzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_yzzzzzzz_0[i] = g_zzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzzzz_1[i] * wa_y[i];

        g_yzzz_0_zzzzzzzz_0[i] = g_zzz_0_zzzzzzzz_1[i] * wa_y[i];
    }

    /// Set up 630-675 components of targeted buffer : GSL

    auto g_zzzz_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_gsl + 630);

    auto g_zzzz_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_gsl + 631);

    auto g_zzzz_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_gsl + 632);

    auto g_zzzz_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_gsl + 633);

    auto g_zzzz_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_gsl + 634);

    auto g_zzzz_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_gsl + 635);

    auto g_zzzz_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_gsl + 636);

    auto g_zzzz_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_gsl + 637);

    auto g_zzzz_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_gsl + 638);

    auto g_zzzz_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_gsl + 639);

    auto g_zzzz_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_gsl + 640);

    auto g_zzzz_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_gsl + 641);

    auto g_zzzz_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_gsl + 642);

    auto g_zzzz_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_gsl + 643);

    auto g_zzzz_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_gsl + 644);

    auto g_zzzz_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 645);

    auto g_zzzz_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 646);

    auto g_zzzz_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 647);

    auto g_zzzz_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 648);

    auto g_zzzz_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 649);

    auto g_zzzz_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 650);

    auto g_zzzz_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 651);

    auto g_zzzz_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 652);

    auto g_zzzz_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 653);

    auto g_zzzz_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 654);

    auto g_zzzz_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 655);

    auto g_zzzz_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 656);

    auto g_zzzz_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 657);

    auto g_zzzz_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 658);

    auto g_zzzz_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 659);

    auto g_zzzz_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 660);

    auto g_zzzz_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 661);

    auto g_zzzz_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 662);

    auto g_zzzz_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 663);

    auto g_zzzz_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 664);

    auto g_zzzz_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 665);

    auto g_zzzz_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_gsl + 666);

    auto g_zzzz_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_gsl + 667);

    auto g_zzzz_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_gsl + 668);

    auto g_zzzz_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_gsl + 669);

    auto g_zzzz_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_gsl + 670);

    auto g_zzzz_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 671);

    auto g_zzzz_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 672);

    auto g_zzzz_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 673);

    auto g_zzzz_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_gsl + 674);

    #pragma omp simd aligned(g_zz_0_xxxxxxxx_0, g_zz_0_xxxxxxxx_1, g_zz_0_xxxxxxxy_0, g_zz_0_xxxxxxxy_1, g_zz_0_xxxxxxxz_0, g_zz_0_xxxxxxxz_1, g_zz_0_xxxxxxyy_0, g_zz_0_xxxxxxyy_1, g_zz_0_xxxxxxyz_0, g_zz_0_xxxxxxyz_1, g_zz_0_xxxxxxzz_0, g_zz_0_xxxxxxzz_1, g_zz_0_xxxxxyyy_0, g_zz_0_xxxxxyyy_1, g_zz_0_xxxxxyyz_0, g_zz_0_xxxxxyyz_1, g_zz_0_xxxxxyzz_0, g_zz_0_xxxxxyzz_1, g_zz_0_xxxxxzzz_0, g_zz_0_xxxxxzzz_1, g_zz_0_xxxxyyyy_0, g_zz_0_xxxxyyyy_1, g_zz_0_xxxxyyyz_0, g_zz_0_xxxxyyyz_1, g_zz_0_xxxxyyzz_0, g_zz_0_xxxxyyzz_1, g_zz_0_xxxxyzzz_0, g_zz_0_xxxxyzzz_1, g_zz_0_xxxxzzzz_0, g_zz_0_xxxxzzzz_1, g_zz_0_xxxyyyyy_0, g_zz_0_xxxyyyyy_1, g_zz_0_xxxyyyyz_0, g_zz_0_xxxyyyyz_1, g_zz_0_xxxyyyzz_0, g_zz_0_xxxyyyzz_1, g_zz_0_xxxyyzzz_0, g_zz_0_xxxyyzzz_1, g_zz_0_xxxyzzzz_0, g_zz_0_xxxyzzzz_1, g_zz_0_xxxzzzzz_0, g_zz_0_xxxzzzzz_1, g_zz_0_xxyyyyyy_0, g_zz_0_xxyyyyyy_1, g_zz_0_xxyyyyyz_0, g_zz_0_xxyyyyyz_1, g_zz_0_xxyyyyzz_0, g_zz_0_xxyyyyzz_1, g_zz_0_xxyyyzzz_0, g_zz_0_xxyyyzzz_1, g_zz_0_xxyyzzzz_0, g_zz_0_xxyyzzzz_1, g_zz_0_xxyzzzzz_0, g_zz_0_xxyzzzzz_1, g_zz_0_xxzzzzzz_0, g_zz_0_xxzzzzzz_1, g_zz_0_xyyyyyyy_0, g_zz_0_xyyyyyyy_1, g_zz_0_xyyyyyyz_0, g_zz_0_xyyyyyyz_1, g_zz_0_xyyyyyzz_0, g_zz_0_xyyyyyzz_1, g_zz_0_xyyyyzzz_0, g_zz_0_xyyyyzzz_1, g_zz_0_xyyyzzzz_0, g_zz_0_xyyyzzzz_1, g_zz_0_xyyzzzzz_0, g_zz_0_xyyzzzzz_1, g_zz_0_xyzzzzzz_0, g_zz_0_xyzzzzzz_1, g_zz_0_xzzzzzzz_0, g_zz_0_xzzzzzzz_1, g_zz_0_yyyyyyyy_0, g_zz_0_yyyyyyyy_1, g_zz_0_yyyyyyyz_0, g_zz_0_yyyyyyyz_1, g_zz_0_yyyyyyzz_0, g_zz_0_yyyyyyzz_1, g_zz_0_yyyyyzzz_0, g_zz_0_yyyyyzzz_1, g_zz_0_yyyyzzzz_0, g_zz_0_yyyyzzzz_1, g_zz_0_yyyzzzzz_0, g_zz_0_yyyzzzzz_1, g_zz_0_yyzzzzzz_0, g_zz_0_yyzzzzzz_1, g_zz_0_yzzzzzzz_0, g_zz_0_yzzzzzzz_1, g_zz_0_zzzzzzzz_0, g_zz_0_zzzzzzzz_1, g_zzz_0_xxxxxxx_1, g_zzz_0_xxxxxxxx_1, g_zzz_0_xxxxxxxy_1, g_zzz_0_xxxxxxxz_1, g_zzz_0_xxxxxxy_1, g_zzz_0_xxxxxxyy_1, g_zzz_0_xxxxxxyz_1, g_zzz_0_xxxxxxz_1, g_zzz_0_xxxxxxzz_1, g_zzz_0_xxxxxyy_1, g_zzz_0_xxxxxyyy_1, g_zzz_0_xxxxxyyz_1, g_zzz_0_xxxxxyz_1, g_zzz_0_xxxxxyzz_1, g_zzz_0_xxxxxzz_1, g_zzz_0_xxxxxzzz_1, g_zzz_0_xxxxyyy_1, g_zzz_0_xxxxyyyy_1, g_zzz_0_xxxxyyyz_1, g_zzz_0_xxxxyyz_1, g_zzz_0_xxxxyyzz_1, g_zzz_0_xxxxyzz_1, g_zzz_0_xxxxyzzz_1, g_zzz_0_xxxxzzz_1, g_zzz_0_xxxxzzzz_1, g_zzz_0_xxxyyyy_1, g_zzz_0_xxxyyyyy_1, g_zzz_0_xxxyyyyz_1, g_zzz_0_xxxyyyz_1, g_zzz_0_xxxyyyzz_1, g_zzz_0_xxxyyzz_1, g_zzz_0_xxxyyzzz_1, g_zzz_0_xxxyzzz_1, g_zzz_0_xxxyzzzz_1, g_zzz_0_xxxzzzz_1, g_zzz_0_xxxzzzzz_1, g_zzz_0_xxyyyyy_1, g_zzz_0_xxyyyyyy_1, g_zzz_0_xxyyyyyz_1, g_zzz_0_xxyyyyz_1, g_zzz_0_xxyyyyzz_1, g_zzz_0_xxyyyzz_1, g_zzz_0_xxyyyzzz_1, g_zzz_0_xxyyzzz_1, g_zzz_0_xxyyzzzz_1, g_zzz_0_xxyzzzz_1, g_zzz_0_xxyzzzzz_1, g_zzz_0_xxzzzzz_1, g_zzz_0_xxzzzzzz_1, g_zzz_0_xyyyyyy_1, g_zzz_0_xyyyyyyy_1, g_zzz_0_xyyyyyyz_1, g_zzz_0_xyyyyyz_1, g_zzz_0_xyyyyyzz_1, g_zzz_0_xyyyyzz_1, g_zzz_0_xyyyyzzz_1, g_zzz_0_xyyyzzz_1, g_zzz_0_xyyyzzzz_1, g_zzz_0_xyyzzzz_1, g_zzz_0_xyyzzzzz_1, g_zzz_0_xyzzzzz_1, g_zzz_0_xyzzzzzz_1, g_zzz_0_xzzzzzz_1, g_zzz_0_xzzzzzzz_1, g_zzz_0_yyyyyyy_1, g_zzz_0_yyyyyyyy_1, g_zzz_0_yyyyyyyz_1, g_zzz_0_yyyyyyz_1, g_zzz_0_yyyyyyzz_1, g_zzz_0_yyyyyzz_1, g_zzz_0_yyyyyzzz_1, g_zzz_0_yyyyzzz_1, g_zzz_0_yyyyzzzz_1, g_zzz_0_yyyzzzz_1, g_zzz_0_yyyzzzzz_1, g_zzz_0_yyzzzzz_1, g_zzz_0_yyzzzzzz_1, g_zzz_0_yzzzzzz_1, g_zzz_0_yzzzzzzz_1, g_zzz_0_zzzzzzz_1, g_zzz_0_zzzzzzzz_1, g_zzzz_0_xxxxxxxx_0, g_zzzz_0_xxxxxxxy_0, g_zzzz_0_xxxxxxxz_0, g_zzzz_0_xxxxxxyy_0, g_zzzz_0_xxxxxxyz_0, g_zzzz_0_xxxxxxzz_0, g_zzzz_0_xxxxxyyy_0, g_zzzz_0_xxxxxyyz_0, g_zzzz_0_xxxxxyzz_0, g_zzzz_0_xxxxxzzz_0, g_zzzz_0_xxxxyyyy_0, g_zzzz_0_xxxxyyyz_0, g_zzzz_0_xxxxyyzz_0, g_zzzz_0_xxxxyzzz_0, g_zzzz_0_xxxxzzzz_0, g_zzzz_0_xxxyyyyy_0, g_zzzz_0_xxxyyyyz_0, g_zzzz_0_xxxyyyzz_0, g_zzzz_0_xxxyyzzz_0, g_zzzz_0_xxxyzzzz_0, g_zzzz_0_xxxzzzzz_0, g_zzzz_0_xxyyyyyy_0, g_zzzz_0_xxyyyyyz_0, g_zzzz_0_xxyyyyzz_0, g_zzzz_0_xxyyyzzz_0, g_zzzz_0_xxyyzzzz_0, g_zzzz_0_xxyzzzzz_0, g_zzzz_0_xxzzzzzz_0, g_zzzz_0_xyyyyyyy_0, g_zzzz_0_xyyyyyyz_0, g_zzzz_0_xyyyyyzz_0, g_zzzz_0_xyyyyzzz_0, g_zzzz_0_xyyyzzzz_0, g_zzzz_0_xyyzzzzz_0, g_zzzz_0_xyzzzzzz_0, g_zzzz_0_xzzzzzzz_0, g_zzzz_0_yyyyyyyy_0, g_zzzz_0_yyyyyyyz_0, g_zzzz_0_yyyyyyzz_0, g_zzzz_0_yyyyyzzz_0, g_zzzz_0_yyyyzzzz_0, g_zzzz_0_yyyzzzzz_0, g_zzzz_0_yyzzzzzz_0, g_zzzz_0_yzzzzzzz_0, g_zzzz_0_zzzzzzzz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xxxxxxxx_0[i] = 3.0 * g_zz_0_xxxxxxxx_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxx_1[i] * fz_be_0 + g_zzz_0_xxxxxxxx_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxy_0[i] = 3.0 * g_zz_0_xxxxxxxy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxy_1[i] * fz_be_0 + g_zzz_0_xxxxxxxy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxxz_0[i] = 3.0 * g_zz_0_xxxxxxxz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxxz_1[i] * fz_be_0 + g_zzz_0_xxxxxxx_1[i] * fi_acd_0 + g_zzz_0_xxxxxxxz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxyy_0[i] = 3.0 * g_zz_0_xxxxxxyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxyy_1[i] * fz_be_0 + g_zzz_0_xxxxxxyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxyz_0[i] = 3.0 * g_zz_0_xxxxxxyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxyz_1[i] * fz_be_0 + g_zzz_0_xxxxxxy_1[i] * fi_acd_0 + g_zzz_0_xxxxxxyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxxzz_0[i] = 3.0 * g_zz_0_xxxxxxzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxxzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxxxz_1[i] * fi_acd_0 + g_zzz_0_xxxxxxzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyyy_0[i] = 3.0 * g_zz_0_xxxxxyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyyy_1[i] * fz_be_0 + g_zzz_0_xxxxxyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyyz_0[i] = 3.0 * g_zz_0_xxxxxyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyyz_1[i] * fz_be_0 + g_zzz_0_xxxxxyy_1[i] * fi_acd_0 + g_zzz_0_xxxxxyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxyzz_0[i] = 3.0 * g_zz_0_xxxxxyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxxyz_1[i] * fi_acd_0 + g_zzz_0_xxxxxyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxxzzz_0[i] = 3.0 * g_zz_0_xxxxxzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxxzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxxxzz_1[i] * fi_acd_0 + g_zzz_0_xxxxxzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyyy_0[i] = 3.0 * g_zz_0_xxxxyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyyy_1[i] * fz_be_0 + g_zzz_0_xxxxyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyyz_0[i] = 3.0 * g_zz_0_xxxxyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyyz_1[i] * fz_be_0 + g_zzz_0_xxxxyyy_1[i] * fi_acd_0 + g_zzz_0_xxxxyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyyzz_0[i] = 3.0 * g_zz_0_xxxxyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxxyyz_1[i] * fi_acd_0 + g_zzz_0_xxxxyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxyzzz_0[i] = 3.0 * g_zz_0_xxxxyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxxyzz_1[i] * fi_acd_0 + g_zzz_0_xxxxyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxxzzzz_0[i] = 3.0 * g_zz_0_xxxxzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxxzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxxxzzz_1[i] * fi_acd_0 + g_zzz_0_xxxxzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyyy_0[i] = 3.0 * g_zz_0_xxxyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyyy_1[i] * fz_be_0 + g_zzz_0_xxxyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyyz_0[i] = 3.0 * g_zz_0_xxxyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyyz_1[i] * fz_be_0 + g_zzz_0_xxxyyyy_1[i] * fi_acd_0 + g_zzz_0_xxxyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyyzz_0[i] = 3.0 * g_zz_0_xxxyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxxyyyz_1[i] * fi_acd_0 + g_zzz_0_xxxyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyyzzz_0[i] = 3.0 * g_zz_0_xxxyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxxyyzz_1[i] * fi_acd_0 + g_zzz_0_xxxyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxyzzzz_0[i] = 3.0 * g_zz_0_xxxyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxxyzzz_1[i] * fi_acd_0 + g_zzz_0_xxxyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxxzzzzz_0[i] = 3.0 * g_zz_0_xxxzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxxzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xxxzzzz_1[i] * fi_acd_0 + g_zzz_0_xxxzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyyy_0[i] = 3.0 * g_zz_0_xxyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyyy_1[i] * fz_be_0 + g_zzz_0_xxyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyyz_0[i] = 3.0 * g_zz_0_xxyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyyz_1[i] * fz_be_0 + g_zzz_0_xxyyyyy_1[i] * fi_acd_0 + g_zzz_0_xxyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyyzz_0[i] = 3.0 * g_zz_0_xxyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xxyyyyz_1[i] * fi_acd_0 + g_zzz_0_xxyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyyzzz_0[i] = 3.0 * g_zz_0_xxyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xxyyyzz_1[i] * fi_acd_0 + g_zzz_0_xxyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyyzzzz_0[i] = 3.0 * g_zz_0_xxyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xxyyzzz_1[i] * fi_acd_0 + g_zzz_0_xxyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxyzzzzz_0[i] = 3.0 * g_zz_0_xxyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xxyzzzz_1[i] * fi_acd_0 + g_zzz_0_xxyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xxzzzzzz_0[i] = 3.0 * g_zz_0_xxzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xxzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_xxzzzzz_1[i] * fi_acd_0 + g_zzz_0_xxzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyyy_0[i] = 3.0 * g_zz_0_xyyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyyy_1[i] * fz_be_0 + g_zzz_0_xyyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyyz_0[i] = 3.0 * g_zz_0_xyyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyyz_1[i] * fz_be_0 + g_zzz_0_xyyyyyy_1[i] * fi_acd_0 + g_zzz_0_xyyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyyzz_0[i] = 3.0 * g_zz_0_xyyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_xyyyyyz_1[i] * fi_acd_0 + g_zzz_0_xyyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyyzzz_0[i] = 3.0 * g_zz_0_xyyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_xyyyyzz_1[i] * fi_acd_0 + g_zzz_0_xyyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyyzzzz_0[i] = 3.0 * g_zz_0_xyyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_xyyyzzz_1[i] * fi_acd_0 + g_zzz_0_xyyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyyzzzzz_0[i] = 3.0 * g_zz_0_xyyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_xyyzzzz_1[i] * fi_acd_0 + g_zzz_0_xyyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xyzzzzzz_0[i] = 3.0 * g_zz_0_xyzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_xyzzzzz_1[i] * fi_acd_0 + g_zzz_0_xyzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_xzzzzzzz_0[i] = 3.0 * g_zz_0_xzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_xzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzz_0_xzzzzzz_1[i] * fi_acd_0 + g_zzz_0_xzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyyy_0[i] = 3.0 * g_zz_0_yyyyyyyy_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyyy_1[i] * fz_be_0 + g_zzz_0_yyyyyyyy_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyyz_0[i] = 3.0 * g_zz_0_yyyyyyyz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyyz_1[i] * fz_be_0 + g_zzz_0_yyyyyyy_1[i] * fi_acd_0 + g_zzz_0_yyyyyyyz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyyzz_0[i] = 3.0 * g_zz_0_yyyyyyzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyyzz_1[i] * fz_be_0 + 2.0 * g_zzz_0_yyyyyyz_1[i] * fi_acd_0 + g_zzz_0_yyyyyyzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyyzzz_0[i] = 3.0 * g_zz_0_yyyyyzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyyzzz_1[i] * fz_be_0 + 3.0 * g_zzz_0_yyyyyzz_1[i] * fi_acd_0 + g_zzz_0_yyyyyzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyyzzzz_0[i] = 3.0 * g_zz_0_yyyyzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyyzzzz_1[i] * fz_be_0 + 4.0 * g_zzz_0_yyyyzzz_1[i] * fi_acd_0 + g_zzz_0_yyyyzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyyzzzzz_0[i] = 3.0 * g_zz_0_yyyzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyyzzzzz_1[i] * fz_be_0 + 5.0 * g_zzz_0_yyyzzzz_1[i] * fi_acd_0 + g_zzz_0_yyyzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yyzzzzzz_0[i] = 3.0 * g_zz_0_yyzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yyzzzzzz_1[i] * fz_be_0 + 6.0 * g_zzz_0_yyzzzzz_1[i] * fi_acd_0 + g_zzz_0_yyzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_yzzzzzzz_0[i] = 3.0 * g_zz_0_yzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_yzzzzzzz_1[i] * fz_be_0 + 7.0 * g_zzz_0_yzzzzzz_1[i] * fi_acd_0 + g_zzz_0_yzzzzzzz_1[i] * wa_z[i];

        g_zzzz_0_zzzzzzzz_0[i] = 3.0 * g_zz_0_zzzzzzzz_0[i] * fbe_0 - 3.0 * g_zz_0_zzzzzzzz_1[i] * fz_be_0 + 8.0 * g_zzz_0_zzzzzzz_1[i] * fi_acd_0 + g_zzz_0_zzzzzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

