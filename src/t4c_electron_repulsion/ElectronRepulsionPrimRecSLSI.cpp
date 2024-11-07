#include "ElectronRepulsionPrimRecSLSI.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_slsi(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_slsi,
                                  size_t                idx_eri_0_sisi,
                                  size_t                idx_eri_1_sisi,
                                  size_t                idx_eri_1_sksh,
                                  size_t                idx_eri_0_sksi,
                                  size_t                idx_eri_1_sksi,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SISI

    auto g_0_xxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi);

    auto g_0_xxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 1);

    auto g_0_xxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 2);

    auto g_0_xxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 3);

    auto g_0_xxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 4);

    auto g_0_xxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 5);

    auto g_0_xxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 6);

    auto g_0_xxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 7);

    auto g_0_xxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 8);

    auto g_0_xxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 9);

    auto g_0_xxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 10);

    auto g_0_xxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 11);

    auto g_0_xxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 12);

    auto g_0_xxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 13);

    auto g_0_xxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 14);

    auto g_0_xxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 15);

    auto g_0_xxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 16);

    auto g_0_xxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 17);

    auto g_0_xxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 18);

    auto g_0_xxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 19);

    auto g_0_xxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 20);

    auto g_0_xxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 21);

    auto g_0_xxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 22);

    auto g_0_xxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 23);

    auto g_0_xxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 24);

    auto g_0_xxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 25);

    auto g_0_xxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 26);

    auto g_0_xxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 27);

    auto g_0_xxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 28);

    auto g_0_xxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 30);

    auto g_0_xxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 33);

    auto g_0_xxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 37);

    auto g_0_xxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 42);

    auto g_0_xxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 48);

    auto g_0_xxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 56);

    auto g_0_xxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 57);

    auto g_0_xxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 59);

    auto g_0_xxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 62);

    auto g_0_xxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 66);

    auto g_0_xxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 71);

    auto g_0_xxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 84);

    auto g_0_xxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 85);

    auto g_0_xxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 86);

    auto g_0_xxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 87);

    auto g_0_xxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 88);

    auto g_0_xxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 89);

    auto g_0_xxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 90);

    auto g_0_xxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 91);

    auto g_0_xxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 92);

    auto g_0_xxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 93);

    auto g_0_xxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 94);

    auto g_0_xxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 95);

    auto g_0_xxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 96);

    auto g_0_xxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 97);

    auto g_0_xxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 98);

    auto g_0_xxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 99);

    auto g_0_xxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 100);

    auto g_0_xxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 101);

    auto g_0_xxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 102);

    auto g_0_xxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 103);

    auto g_0_xxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 104);

    auto g_0_xxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 105);

    auto g_0_xxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 106);

    auto g_0_xxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 107);

    auto g_0_xxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 108);

    auto g_0_xxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 109);

    auto g_0_xxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 110);

    auto g_0_xxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 111);

    auto g_0_xxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 140);

    auto g_0_xxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 141);

    auto g_0_xxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 142);

    auto g_0_xxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 143);

    auto g_0_xxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 144);

    auto g_0_xxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 145);

    auto g_0_xxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 146);

    auto g_0_xxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 147);

    auto g_0_xxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 148);

    auto g_0_xxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 149);

    auto g_0_xxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 150);

    auto g_0_xxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 151);

    auto g_0_xxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 152);

    auto g_0_xxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 153);

    auto g_0_xxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 154);

    auto g_0_xxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 155);

    auto g_0_xxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 156);

    auto g_0_xxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 157);

    auto g_0_xxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 158);

    auto g_0_xxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 159);

    auto g_0_xxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 160);

    auto g_0_xxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 161);

    auto g_0_xxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 162);

    auto g_0_xxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 163);

    auto g_0_xxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 164);

    auto g_0_xxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 165);

    auto g_0_xxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 166);

    auto g_0_xxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 167);

    auto g_0_xxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 168);

    auto g_0_xxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 169);

    auto g_0_xxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 170);

    auto g_0_xxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 171);

    auto g_0_xxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 172);

    auto g_0_xxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 173);

    auto g_0_xxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 174);

    auto g_0_xxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 175);

    auto g_0_xxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 176);

    auto g_0_xxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 177);

    auto g_0_xxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 178);

    auto g_0_xxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 179);

    auto g_0_xxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 180);

    auto g_0_xxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 181);

    auto g_0_xxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 182);

    auto g_0_xxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 183);

    auto g_0_xxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 184);

    auto g_0_xxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 185);

    auto g_0_xxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 186);

    auto g_0_xxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 187);

    auto g_0_xxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 188);

    auto g_0_xxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 189);

    auto g_0_xxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 190);

    auto g_0_xxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 191);

    auto g_0_xxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 192);

    auto g_0_xxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 193);

    auto g_0_xxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 194);

    auto g_0_xxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 195);

    auto g_0_xxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 197);

    auto g_0_xxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 199);

    auto g_0_xxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 202);

    auto g_0_xxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 206);

    auto g_0_xxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 211);

    auto g_0_xxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 224);

    auto g_0_xxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 226);

    auto g_0_xxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 229);

    auto g_0_xxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 233);

    auto g_0_xxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 238);

    auto g_0_xxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 244);

    auto g_0_xxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 252);

    auto g_0_xxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 253);

    auto g_0_xxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 254);

    auto g_0_xxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 255);

    auto g_0_xxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 256);

    auto g_0_xxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 257);

    auto g_0_xxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 258);

    auto g_0_xxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 259);

    auto g_0_xxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 260);

    auto g_0_xxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 261);

    auto g_0_xxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 262);

    auto g_0_xxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 263);

    auto g_0_xxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 264);

    auto g_0_xxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 265);

    auto g_0_xxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 266);

    auto g_0_xxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 267);

    auto g_0_xxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 268);

    auto g_0_xxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 269);

    auto g_0_xxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 270);

    auto g_0_xxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 271);

    auto g_0_xxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 272);

    auto g_0_xxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 273);

    auto g_0_xxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 274);

    auto g_0_xxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 275);

    auto g_0_xxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 276);

    auto g_0_xxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 277);

    auto g_0_xxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 278);

    auto g_0_xxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 279);

    auto g_0_xxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 280);

    auto g_0_xxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 281);

    auto g_0_xxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 282);

    auto g_0_xxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 283);

    auto g_0_xxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 284);

    auto g_0_xxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 285);

    auto g_0_xxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 286);

    auto g_0_xxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 287);

    auto g_0_xxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 288);

    auto g_0_xxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 289);

    auto g_0_xxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 290);

    auto g_0_xxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 291);

    auto g_0_xxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 292);

    auto g_0_xxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 293);

    auto g_0_xxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 294);

    auto g_0_xxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 295);

    auto g_0_xxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 296);

    auto g_0_xxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 297);

    auto g_0_xxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 298);

    auto g_0_xxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 299);

    auto g_0_xxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 300);

    auto g_0_xxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 301);

    auto g_0_xxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 302);

    auto g_0_xxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 303);

    auto g_0_xxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 304);

    auto g_0_xxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 305);

    auto g_0_xxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 306);

    auto g_0_xxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 307);

    auto g_0_xxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 309);

    auto g_0_xxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 311);

    auto g_0_xxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 314);

    auto g_0_xxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 318);

    auto g_0_xxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 323);

    auto g_0_xxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 336);

    auto g_0_xxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 337);

    auto g_0_xxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 338);

    auto g_0_xxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 339);

    auto g_0_xxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 340);

    auto g_0_xxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 341);

    auto g_0_xxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 342);

    auto g_0_xxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 343);

    auto g_0_xxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 344);

    auto g_0_xxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 345);

    auto g_0_xxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 346);

    auto g_0_xxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 347);

    auto g_0_xxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 348);

    auto g_0_xxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 349);

    auto g_0_xxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 350);

    auto g_0_xxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 351);

    auto g_0_xxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 352);

    auto g_0_xxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 353);

    auto g_0_xxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 354);

    auto g_0_xxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 355);

    auto g_0_xxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 356);

    auto g_0_xxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 357);

    auto g_0_xxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 358);

    auto g_0_xxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 359);

    auto g_0_xxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 360);

    auto g_0_xxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 361);

    auto g_0_xxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 362);

    auto g_0_xxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 363);

    auto g_0_xxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 364);

    auto g_0_xxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 366);

    auto g_0_xxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 369);

    auto g_0_xxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 373);

    auto g_0_xxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 378);

    auto g_0_xxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 384);

    auto g_0_xxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 392);

    auto g_0_xxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 393);

    auto g_0_xxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 394);

    auto g_0_xxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 395);

    auto g_0_xxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 396);

    auto g_0_xxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 397);

    auto g_0_xxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 398);

    auto g_0_xxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 399);

    auto g_0_xxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 400);

    auto g_0_xxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 401);

    auto g_0_xxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 402);

    auto g_0_xxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 403);

    auto g_0_xxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 404);

    auto g_0_xxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 405);

    auto g_0_xxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 406);

    auto g_0_xxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 407);

    auto g_0_xxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 408);

    auto g_0_xxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 409);

    auto g_0_xxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 410);

    auto g_0_xxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 411);

    auto g_0_xxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 412);

    auto g_0_xxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 413);

    auto g_0_xxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 414);

    auto g_0_xxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 415);

    auto g_0_xxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 416);

    auto g_0_xxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 417);

    auto g_0_xxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 418);

    auto g_0_xxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 419);

    auto g_0_xyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 421);

    auto g_0_xyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 423);

    auto g_0_xyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 424);

    auto g_0_xyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 426);

    auto g_0_xyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 427);

    auto g_0_xyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 428);

    auto g_0_xyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 430);

    auto g_0_xyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 431);

    auto g_0_xyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 432);

    auto g_0_xyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 433);

    auto g_0_xyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 435);

    auto g_0_xyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 436);

    auto g_0_xyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 437);

    auto g_0_xyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 438);

    auto g_0_xyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 439);

    auto g_0_xyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 441);

    auto g_0_xyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 442);

    auto g_0_xyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 443);

    auto g_0_xyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 444);

    auto g_0_xyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 445);

    auto g_0_xyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 446);

    auto g_0_xyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 447);

    auto g_0_xyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 480);

    auto g_0_xyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 483);

    auto g_0_xyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 484);

    auto g_0_xyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 487);

    auto g_0_xyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 488);

    auto g_0_xyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 489);

    auto g_0_xyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 492);

    auto g_0_xyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 493);

    auto g_0_xyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 494);

    auto g_0_xyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 495);

    auto g_0_xyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 497);

    auto g_0_xyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 498);

    auto g_0_xyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 499);

    auto g_0_xyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 500);

    auto g_0_xyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 501);

    auto g_0_xyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 502);

    auto g_0_xyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 503);

    auto g_0_xyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 508);

    auto g_0_xyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 511);

    auto g_0_xyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 512);

    auto g_0_xyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 515);

    auto g_0_xyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 516);

    auto g_0_xyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 517);

    auto g_0_xyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 520);

    auto g_0_xyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 521);

    auto g_0_xyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 522);

    auto g_0_xyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 523);

    auto g_0_xyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 525);

    auto g_0_xyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 526);

    auto g_0_xyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 527);

    auto g_0_xyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 528);

    auto g_0_xyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 529);

    auto g_0_xyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 530);

    auto g_0_xyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 531);

    auto g_0_xzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 562);

    auto g_0_xzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 564);

    auto g_0_xzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 565);

    auto g_0_xzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 567);

    auto g_0_xzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 568);

    auto g_0_xzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 569);

    auto g_0_xzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 571);

    auto g_0_xzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 572);

    auto g_0_xzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 573);

    auto g_0_xzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 574);

    auto g_0_xzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 576);

    auto g_0_xzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 577);

    auto g_0_xzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 578);

    auto g_0_xzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 579);

    auto g_0_xzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 580);

    auto g_0_xzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 581);

    auto g_0_xzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 582);

    auto g_0_xzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 583);

    auto g_0_xzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 584);

    auto g_0_xzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 585);

    auto g_0_xzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 586);

    auto g_0_xzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 587);

    auto g_0_yyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 588);

    auto g_0_yyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 589);

    auto g_0_yyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 590);

    auto g_0_yyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 591);

    auto g_0_yyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 592);

    auto g_0_yyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 593);

    auto g_0_yyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 594);

    auto g_0_yyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 595);

    auto g_0_yyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 596);

    auto g_0_yyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 597);

    auto g_0_yyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 598);

    auto g_0_yyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 599);

    auto g_0_yyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 600);

    auto g_0_yyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 601);

    auto g_0_yyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 602);

    auto g_0_yyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 603);

    auto g_0_yyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 604);

    auto g_0_yyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 605);

    auto g_0_yyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 606);

    auto g_0_yyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 607);

    auto g_0_yyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 608);

    auto g_0_yyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 609);

    auto g_0_yyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 610);

    auto g_0_yyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 611);

    auto g_0_yyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 612);

    auto g_0_yyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 613);

    auto g_0_yyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 614);

    auto g_0_yyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 615);

    auto g_0_yyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 617);

    auto g_0_yyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 619);

    auto g_0_yyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 622);

    auto g_0_yyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 626);

    auto g_0_yyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 631);

    auto g_0_yyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 637);

    auto g_0_yyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 644);

    auto g_0_yyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 645);

    auto g_0_yyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 646);

    auto g_0_yyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 647);

    auto g_0_yyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 648);

    auto g_0_yyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 649);

    auto g_0_yyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 650);

    auto g_0_yyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 651);

    auto g_0_yyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 652);

    auto g_0_yyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 653);

    auto g_0_yyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 654);

    auto g_0_yyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 655);

    auto g_0_yyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 656);

    auto g_0_yyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 657);

    auto g_0_yyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 658);

    auto g_0_yyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 659);

    auto g_0_yyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 660);

    auto g_0_yyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 661);

    auto g_0_yyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 662);

    auto g_0_yyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 663);

    auto g_0_yyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 664);

    auto g_0_yyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 665);

    auto g_0_yyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 666);

    auto g_0_yyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 667);

    auto g_0_yyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 668);

    auto g_0_yyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 669);

    auto g_0_yyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 670);

    auto g_0_yyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 671);

    auto g_0_yyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 672);

    auto g_0_yyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 673);

    auto g_0_yyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 674);

    auto g_0_yyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 675);

    auto g_0_yyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 676);

    auto g_0_yyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 677);

    auto g_0_yyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 678);

    auto g_0_yyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 679);

    auto g_0_yyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 680);

    auto g_0_yyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 681);

    auto g_0_yyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 682);

    auto g_0_yyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 683);

    auto g_0_yyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 684);

    auto g_0_yyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 685);

    auto g_0_yyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 686);

    auto g_0_yyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 687);

    auto g_0_yyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 688);

    auto g_0_yyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 689);

    auto g_0_yyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 690);

    auto g_0_yyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 691);

    auto g_0_yyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 692);

    auto g_0_yyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 693);

    auto g_0_yyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 694);

    auto g_0_yyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 695);

    auto g_0_yyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 696);

    auto g_0_yyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 697);

    auto g_0_yyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 698);

    auto g_0_yyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 699);

    auto g_0_yyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 700);

    auto g_0_yyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 701);

    auto g_0_yyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 702);

    auto g_0_yyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 703);

    auto g_0_yyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 704);

    auto g_0_yyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 705);

    auto g_0_yyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 706);

    auto g_0_yyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 707);

    auto g_0_yyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 708);

    auto g_0_yyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 709);

    auto g_0_yyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 710);

    auto g_0_yyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 711);

    auto g_0_yyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 712);

    auto g_0_yyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 713);

    auto g_0_yyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 714);

    auto g_0_yyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 715);

    auto g_0_yyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 716);

    auto g_0_yyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 717);

    auto g_0_yyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 718);

    auto g_0_yyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 719);

    auto g_0_yyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 720);

    auto g_0_yyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 721);

    auto g_0_yyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 722);

    auto g_0_yyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 723);

    auto g_0_yyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 724);

    auto g_0_yyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 725);

    auto g_0_yyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 726);

    auto g_0_yyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 727);

    auto g_0_yzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 728);

    auto g_0_yzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 730);

    auto g_0_yzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 732);

    auto g_0_yzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 733);

    auto g_0_yzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 735);

    auto g_0_yzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 736);

    auto g_0_yzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 737);

    auto g_0_yzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 739);

    auto g_0_yzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 740);

    auto g_0_yzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 741);

    auto g_0_yzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 742);

    auto g_0_yzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 744);

    auto g_0_yzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 745);

    auto g_0_yzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 746);

    auto g_0_yzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 747);

    auto g_0_yzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 748);

    auto g_0_yzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 750);

    auto g_0_yzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 751);

    auto g_0_yzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 752);

    auto g_0_yzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 753);

    auto g_0_yzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 754);

    auto g_0_yzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 755);

    auto g_0_zzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sisi + 756);

    auto g_0_zzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sisi + 757);

    auto g_0_zzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sisi + 758);

    auto g_0_zzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sisi + 759);

    auto g_0_zzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sisi + 760);

    auto g_0_zzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sisi + 761);

    auto g_0_zzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sisi + 762);

    auto g_0_zzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sisi + 763);

    auto g_0_zzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sisi + 764);

    auto g_0_zzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sisi + 765);

    auto g_0_zzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sisi + 766);

    auto g_0_zzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sisi + 767);

    auto g_0_zzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sisi + 768);

    auto g_0_zzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sisi + 769);

    auto g_0_zzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sisi + 770);

    auto g_0_zzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 771);

    auto g_0_zzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 772);

    auto g_0_zzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 773);

    auto g_0_zzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 774);

    auto g_0_zzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 775);

    auto g_0_zzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 776);

    auto g_0_zzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sisi + 777);

    auto g_0_zzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sisi + 778);

    auto g_0_zzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sisi + 779);

    auto g_0_zzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sisi + 780);

    auto g_0_zzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sisi + 781);

    auto g_0_zzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 782);

    auto g_0_zzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sisi + 783);

    /// Set up components of auxilary buffer : SISI

    auto g_0_xxxxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi);

    auto g_0_xxxxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 1);

    auto g_0_xxxxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 2);

    auto g_0_xxxxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 3);

    auto g_0_xxxxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 4);

    auto g_0_xxxxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 5);

    auto g_0_xxxxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 6);

    auto g_0_xxxxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 7);

    auto g_0_xxxxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 8);

    auto g_0_xxxxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 9);

    auto g_0_xxxxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 10);

    auto g_0_xxxxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 11);

    auto g_0_xxxxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 12);

    auto g_0_xxxxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 13);

    auto g_0_xxxxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 14);

    auto g_0_xxxxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 15);

    auto g_0_xxxxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 16);

    auto g_0_xxxxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 17);

    auto g_0_xxxxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 18);

    auto g_0_xxxxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 19);

    auto g_0_xxxxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 20);

    auto g_0_xxxxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 21);

    auto g_0_xxxxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 22);

    auto g_0_xxxxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 23);

    auto g_0_xxxxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 24);

    auto g_0_xxxxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 25);

    auto g_0_xxxxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 26);

    auto g_0_xxxxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 27);

    auto g_0_xxxxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 28);

    auto g_0_xxxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 30);

    auto g_0_xxxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 33);

    auto g_0_xxxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 37);

    auto g_0_xxxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 42);

    auto g_0_xxxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 48);

    auto g_0_xxxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 56);

    auto g_0_xxxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 57);

    auto g_0_xxxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 59);

    auto g_0_xxxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 62);

    auto g_0_xxxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 66);

    auto g_0_xxxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 71);

    auto g_0_xxxxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 84);

    auto g_0_xxxxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 85);

    auto g_0_xxxxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 86);

    auto g_0_xxxxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 87);

    auto g_0_xxxxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 88);

    auto g_0_xxxxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 89);

    auto g_0_xxxxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 90);

    auto g_0_xxxxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 91);

    auto g_0_xxxxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 92);

    auto g_0_xxxxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 93);

    auto g_0_xxxxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 94);

    auto g_0_xxxxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 95);

    auto g_0_xxxxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 96);

    auto g_0_xxxxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 97);

    auto g_0_xxxxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 98);

    auto g_0_xxxxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 99);

    auto g_0_xxxxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 100);

    auto g_0_xxxxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 101);

    auto g_0_xxxxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 102);

    auto g_0_xxxxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 103);

    auto g_0_xxxxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 104);

    auto g_0_xxxxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 105);

    auto g_0_xxxxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 106);

    auto g_0_xxxxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 107);

    auto g_0_xxxxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 108);

    auto g_0_xxxxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 109);

    auto g_0_xxxxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 110);

    auto g_0_xxxxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 111);

    auto g_0_xxxxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 140);

    auto g_0_xxxxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 141);

    auto g_0_xxxxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 142);

    auto g_0_xxxxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 143);

    auto g_0_xxxxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 144);

    auto g_0_xxxxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 145);

    auto g_0_xxxxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 146);

    auto g_0_xxxxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 147);

    auto g_0_xxxxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 148);

    auto g_0_xxxxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 149);

    auto g_0_xxxxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 150);

    auto g_0_xxxxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 151);

    auto g_0_xxxxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 152);

    auto g_0_xxxxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 153);

    auto g_0_xxxxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 154);

    auto g_0_xxxxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 155);

    auto g_0_xxxxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 156);

    auto g_0_xxxxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 157);

    auto g_0_xxxxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 158);

    auto g_0_xxxxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 159);

    auto g_0_xxxxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 160);

    auto g_0_xxxxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 161);

    auto g_0_xxxxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 162);

    auto g_0_xxxxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 163);

    auto g_0_xxxxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 164);

    auto g_0_xxxxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 165);

    auto g_0_xxxxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 166);

    auto g_0_xxxxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 167);

    auto g_0_xxxyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 168);

    auto g_0_xxxyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 169);

    auto g_0_xxxyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 170);

    auto g_0_xxxyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 171);

    auto g_0_xxxyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 172);

    auto g_0_xxxyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 173);

    auto g_0_xxxyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 174);

    auto g_0_xxxyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 175);

    auto g_0_xxxyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 176);

    auto g_0_xxxyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 177);

    auto g_0_xxxyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 178);

    auto g_0_xxxyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 179);

    auto g_0_xxxyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 180);

    auto g_0_xxxyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 181);

    auto g_0_xxxyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 182);

    auto g_0_xxxyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 183);

    auto g_0_xxxyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 184);

    auto g_0_xxxyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 185);

    auto g_0_xxxyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 186);

    auto g_0_xxxyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 187);

    auto g_0_xxxyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 188);

    auto g_0_xxxyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 189);

    auto g_0_xxxyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 190);

    auto g_0_xxxyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 191);

    auto g_0_xxxyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 192);

    auto g_0_xxxyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 193);

    auto g_0_xxxyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 194);

    auto g_0_xxxyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 195);

    auto g_0_xxxyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 197);

    auto g_0_xxxyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 199);

    auto g_0_xxxyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 202);

    auto g_0_xxxyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 206);

    auto g_0_xxxyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 211);

    auto g_0_xxxyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 224);

    auto g_0_xxxyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 226);

    auto g_0_xxxyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 229);

    auto g_0_xxxyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 233);

    auto g_0_xxxyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 238);

    auto g_0_xxxyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 244);

    auto g_0_xxxzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 252);

    auto g_0_xxxzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 253);

    auto g_0_xxxzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 254);

    auto g_0_xxxzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 255);

    auto g_0_xxxzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 256);

    auto g_0_xxxzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 257);

    auto g_0_xxxzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 258);

    auto g_0_xxxzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 259);

    auto g_0_xxxzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 260);

    auto g_0_xxxzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 261);

    auto g_0_xxxzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 262);

    auto g_0_xxxzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 263);

    auto g_0_xxxzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 264);

    auto g_0_xxxzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 265);

    auto g_0_xxxzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 266);

    auto g_0_xxxzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 267);

    auto g_0_xxxzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 268);

    auto g_0_xxxzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 269);

    auto g_0_xxxzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 270);

    auto g_0_xxxzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 271);

    auto g_0_xxxzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 272);

    auto g_0_xxxzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 273);

    auto g_0_xxxzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 274);

    auto g_0_xxxzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 275);

    auto g_0_xxxzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 276);

    auto g_0_xxxzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 277);

    auto g_0_xxxzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 278);

    auto g_0_xxxzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 279);

    auto g_0_xxyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 280);

    auto g_0_xxyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 281);

    auto g_0_xxyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 282);

    auto g_0_xxyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 283);

    auto g_0_xxyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 284);

    auto g_0_xxyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 285);

    auto g_0_xxyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 286);

    auto g_0_xxyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 287);

    auto g_0_xxyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 288);

    auto g_0_xxyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 289);

    auto g_0_xxyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 290);

    auto g_0_xxyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 291);

    auto g_0_xxyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 292);

    auto g_0_xxyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 293);

    auto g_0_xxyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 294);

    auto g_0_xxyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 295);

    auto g_0_xxyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 296);

    auto g_0_xxyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 297);

    auto g_0_xxyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 298);

    auto g_0_xxyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 299);

    auto g_0_xxyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 300);

    auto g_0_xxyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 301);

    auto g_0_xxyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 302);

    auto g_0_xxyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 303);

    auto g_0_xxyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 304);

    auto g_0_xxyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 305);

    auto g_0_xxyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 306);

    auto g_0_xxyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 307);

    auto g_0_xxyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 309);

    auto g_0_xxyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 311);

    auto g_0_xxyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 314);

    auto g_0_xxyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 318);

    auto g_0_xxyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 323);

    auto g_0_xxyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 336);

    auto g_0_xxyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 337);

    auto g_0_xxyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 338);

    auto g_0_xxyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 339);

    auto g_0_xxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 340);

    auto g_0_xxyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 341);

    auto g_0_xxyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 342);

    auto g_0_xxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 343);

    auto g_0_xxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 344);

    auto g_0_xxyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 345);

    auto g_0_xxyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 346);

    auto g_0_xxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 347);

    auto g_0_xxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 348);

    auto g_0_xxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 349);

    auto g_0_xxyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 350);

    auto g_0_xxyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 351);

    auto g_0_xxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 352);

    auto g_0_xxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 353);

    auto g_0_xxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 354);

    auto g_0_xxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 355);

    auto g_0_xxyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 356);

    auto g_0_xxyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 357);

    auto g_0_xxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 358);

    auto g_0_xxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 359);

    auto g_0_xxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 360);

    auto g_0_xxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 361);

    auto g_0_xxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 362);

    auto g_0_xxyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 363);

    auto g_0_xxyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 364);

    auto g_0_xxyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 366);

    auto g_0_xxyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 369);

    auto g_0_xxyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 373);

    auto g_0_xxyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 378);

    auto g_0_xxyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 384);

    auto g_0_xxzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 392);

    auto g_0_xxzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 393);

    auto g_0_xxzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 394);

    auto g_0_xxzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 395);

    auto g_0_xxzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 396);

    auto g_0_xxzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 397);

    auto g_0_xxzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 398);

    auto g_0_xxzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 399);

    auto g_0_xxzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 400);

    auto g_0_xxzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 401);

    auto g_0_xxzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 402);

    auto g_0_xxzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 403);

    auto g_0_xxzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 404);

    auto g_0_xxzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 405);

    auto g_0_xxzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 406);

    auto g_0_xxzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 407);

    auto g_0_xxzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 408);

    auto g_0_xxzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 409);

    auto g_0_xxzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 410);

    auto g_0_xxzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 411);

    auto g_0_xxzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 412);

    auto g_0_xxzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 413);

    auto g_0_xxzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 414);

    auto g_0_xxzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 415);

    auto g_0_xxzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 416);

    auto g_0_xxzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 417);

    auto g_0_xxzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 418);

    auto g_0_xxzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 419);

    auto g_0_xyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 421);

    auto g_0_xyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 423);

    auto g_0_xyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 424);

    auto g_0_xyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 426);

    auto g_0_xyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 427);

    auto g_0_xyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 428);

    auto g_0_xyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 430);

    auto g_0_xyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 431);

    auto g_0_xyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 432);

    auto g_0_xyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 433);

    auto g_0_xyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 435);

    auto g_0_xyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 436);

    auto g_0_xyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 437);

    auto g_0_xyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 438);

    auto g_0_xyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 439);

    auto g_0_xyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 441);

    auto g_0_xyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 442);

    auto g_0_xyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 443);

    auto g_0_xyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 444);

    auto g_0_xyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 445);

    auto g_0_xyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 446);

    auto g_0_xyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 447);

    auto g_0_xyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 480);

    auto g_0_xyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 483);

    auto g_0_xyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 484);

    auto g_0_xyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 487);

    auto g_0_xyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 488);

    auto g_0_xyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 489);

    auto g_0_xyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 492);

    auto g_0_xyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 493);

    auto g_0_xyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 494);

    auto g_0_xyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 495);

    auto g_0_xyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 497);

    auto g_0_xyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 498);

    auto g_0_xyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 499);

    auto g_0_xyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 500);

    auto g_0_xyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 501);

    auto g_0_xyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 502);

    auto g_0_xyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 503);

    auto g_0_xyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 508);

    auto g_0_xyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 511);

    auto g_0_xyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 512);

    auto g_0_xyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 515);

    auto g_0_xyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 516);

    auto g_0_xyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 517);

    auto g_0_xyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 520);

    auto g_0_xyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 521);

    auto g_0_xyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 522);

    auto g_0_xyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 523);

    auto g_0_xyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 525);

    auto g_0_xyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 526);

    auto g_0_xyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 527);

    auto g_0_xyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 528);

    auto g_0_xyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 529);

    auto g_0_xyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 530);

    auto g_0_xyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 531);

    auto g_0_xzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 562);

    auto g_0_xzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 564);

    auto g_0_xzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 565);

    auto g_0_xzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 567);

    auto g_0_xzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 568);

    auto g_0_xzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 569);

    auto g_0_xzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 571);

    auto g_0_xzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 572);

    auto g_0_xzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 573);

    auto g_0_xzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 574);

    auto g_0_xzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 576);

    auto g_0_xzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 577);

    auto g_0_xzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 578);

    auto g_0_xzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 579);

    auto g_0_xzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 580);

    auto g_0_xzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 581);

    auto g_0_xzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 582);

    auto g_0_xzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 583);

    auto g_0_xzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 584);

    auto g_0_xzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 585);

    auto g_0_xzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 586);

    auto g_0_xzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 587);

    auto g_0_yyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 588);

    auto g_0_yyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 589);

    auto g_0_yyyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 590);

    auto g_0_yyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 591);

    auto g_0_yyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 592);

    auto g_0_yyyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 593);

    auto g_0_yyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 594);

    auto g_0_yyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 595);

    auto g_0_yyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 596);

    auto g_0_yyyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 597);

    auto g_0_yyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 598);

    auto g_0_yyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 599);

    auto g_0_yyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 600);

    auto g_0_yyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 601);

    auto g_0_yyyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 602);

    auto g_0_yyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 603);

    auto g_0_yyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 604);

    auto g_0_yyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 605);

    auto g_0_yyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 606);

    auto g_0_yyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 607);

    auto g_0_yyyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 608);

    auto g_0_yyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 609);

    auto g_0_yyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 610);

    auto g_0_yyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 611);

    auto g_0_yyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 612);

    auto g_0_yyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 613);

    auto g_0_yyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 614);

    auto g_0_yyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 615);

    auto g_0_yyyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 617);

    auto g_0_yyyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 619);

    auto g_0_yyyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 622);

    auto g_0_yyyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 626);

    auto g_0_yyyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 631);

    auto g_0_yyyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 637);

    auto g_0_yyyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 644);

    auto g_0_yyyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 645);

    auto g_0_yyyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 646);

    auto g_0_yyyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 647);

    auto g_0_yyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 648);

    auto g_0_yyyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 649);

    auto g_0_yyyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 650);

    auto g_0_yyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 651);

    auto g_0_yyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 652);

    auto g_0_yyyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 653);

    auto g_0_yyyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 654);

    auto g_0_yyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 655);

    auto g_0_yyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 656);

    auto g_0_yyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 657);

    auto g_0_yyyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 658);

    auto g_0_yyyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 659);

    auto g_0_yyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 660);

    auto g_0_yyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 661);

    auto g_0_yyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 662);

    auto g_0_yyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 663);

    auto g_0_yyyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 664);

    auto g_0_yyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 665);

    auto g_0_yyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 666);

    auto g_0_yyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 667);

    auto g_0_yyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 668);

    auto g_0_yyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 669);

    auto g_0_yyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 670);

    auto g_0_yyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 671);

    auto g_0_yyyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 672);

    auto g_0_yyyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 673);

    auto g_0_yyyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 674);

    auto g_0_yyyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 675);

    auto g_0_yyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 676);

    auto g_0_yyyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 677);

    auto g_0_yyyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 678);

    auto g_0_yyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 679);

    auto g_0_yyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 680);

    auto g_0_yyyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 681);

    auto g_0_yyyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 682);

    auto g_0_yyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 683);

    auto g_0_yyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 684);

    auto g_0_yyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 685);

    auto g_0_yyyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 686);

    auto g_0_yyyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 687);

    auto g_0_yyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 688);

    auto g_0_yyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 689);

    auto g_0_yyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 690);

    auto g_0_yyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 691);

    auto g_0_yyyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 692);

    auto g_0_yyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 693);

    auto g_0_yyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 694);

    auto g_0_yyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 695);

    auto g_0_yyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 696);

    auto g_0_yyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 697);

    auto g_0_yyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 698);

    auto g_0_yyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 699);

    auto g_0_yyzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 700);

    auto g_0_yyzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 701);

    auto g_0_yyzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 702);

    auto g_0_yyzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 703);

    auto g_0_yyzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 704);

    auto g_0_yyzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 705);

    auto g_0_yyzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 706);

    auto g_0_yyzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 707);

    auto g_0_yyzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 708);

    auto g_0_yyzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 709);

    auto g_0_yyzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 710);

    auto g_0_yyzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 711);

    auto g_0_yyzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 712);

    auto g_0_yyzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 713);

    auto g_0_yyzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 714);

    auto g_0_yyzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 715);

    auto g_0_yyzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 716);

    auto g_0_yyzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 717);

    auto g_0_yyzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 718);

    auto g_0_yyzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 719);

    auto g_0_yyzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 720);

    auto g_0_yyzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 721);

    auto g_0_yyzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 722);

    auto g_0_yyzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 723);

    auto g_0_yyzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 724);

    auto g_0_yyzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 725);

    auto g_0_yyzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 726);

    auto g_0_yyzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 727);

    auto g_0_yzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 728);

    auto g_0_yzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 730);

    auto g_0_yzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 732);

    auto g_0_yzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 733);

    auto g_0_yzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 735);

    auto g_0_yzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 736);

    auto g_0_yzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 737);

    auto g_0_yzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 739);

    auto g_0_yzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 740);

    auto g_0_yzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 741);

    auto g_0_yzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 742);

    auto g_0_yzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 744);

    auto g_0_yzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 745);

    auto g_0_yzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 746);

    auto g_0_yzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 747);

    auto g_0_yzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 748);

    auto g_0_yzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 750);

    auto g_0_yzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 751);

    auto g_0_yzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 752);

    auto g_0_yzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 753);

    auto g_0_yzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 754);

    auto g_0_yzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 755);

    auto g_0_zzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sisi + 756);

    auto g_0_zzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sisi + 757);

    auto g_0_zzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sisi + 758);

    auto g_0_zzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sisi + 759);

    auto g_0_zzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sisi + 760);

    auto g_0_zzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sisi + 761);

    auto g_0_zzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sisi + 762);

    auto g_0_zzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sisi + 763);

    auto g_0_zzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sisi + 764);

    auto g_0_zzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sisi + 765);

    auto g_0_zzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sisi + 766);

    auto g_0_zzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sisi + 767);

    auto g_0_zzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sisi + 768);

    auto g_0_zzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sisi + 769);

    auto g_0_zzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sisi + 770);

    auto g_0_zzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 771);

    auto g_0_zzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 772);

    auto g_0_zzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 773);

    auto g_0_zzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 774);

    auto g_0_zzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 775);

    auto g_0_zzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 776);

    auto g_0_zzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sisi + 777);

    auto g_0_zzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sisi + 778);

    auto g_0_zzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sisi + 779);

    auto g_0_zzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sisi + 780);

    auto g_0_zzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sisi + 781);

    auto g_0_zzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 782);

    auto g_0_zzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sisi + 783);

    /// Set up components of auxilary buffer : SKSH

    auto g_0_xxxxxxx_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh);

    auto g_0_xxxxxxx_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 1);

    auto g_0_xxxxxxx_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 2);

    auto g_0_xxxxxxx_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 3);

    auto g_0_xxxxxxx_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 4);

    auto g_0_xxxxxxx_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 5);

    auto g_0_xxxxxxx_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 6);

    auto g_0_xxxxxxx_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 7);

    auto g_0_xxxxxxx_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 8);

    auto g_0_xxxxxxx_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 9);

    auto g_0_xxxxxxx_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 10);

    auto g_0_xxxxxxx_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 11);

    auto g_0_xxxxxxx_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 12);

    auto g_0_xxxxxxx_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 13);

    auto g_0_xxxxxxx_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 14);

    auto g_0_xxxxxxx_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 15);

    auto g_0_xxxxxxx_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 16);

    auto g_0_xxxxxxx_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 17);

    auto g_0_xxxxxxx_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 18);

    auto g_0_xxxxxxx_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 19);

    auto g_0_xxxxxxx_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 20);

    auto g_0_xxxxxxz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 44);

    auto g_0_xxxxxxz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 46);

    auto g_0_xxxxxxz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 47);

    auto g_0_xxxxxxz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 49);

    auto g_0_xxxxxxz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 50);

    auto g_0_xxxxxxz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 51);

    auto g_0_xxxxxxz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 53);

    auto g_0_xxxxxxz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 54);

    auto g_0_xxxxxxz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 55);

    auto g_0_xxxxxxz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 56);

    auto g_0_xxxxxxz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 58);

    auto g_0_xxxxxxz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 59);

    auto g_0_xxxxxxz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 60);

    auto g_0_xxxxxxz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 61);

    auto g_0_xxxxxxz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 62);

    auto g_0_xxxxxyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 63);

    auto g_0_xxxxxyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 64);

    auto g_0_xxxxxyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 65);

    auto g_0_xxxxxyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 66);

    auto g_0_xxxxxyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 67);

    auto g_0_xxxxxyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 68);

    auto g_0_xxxxxyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 69);

    auto g_0_xxxxxyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 70);

    auto g_0_xxxxxyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 71);

    auto g_0_xxxxxyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 72);

    auto g_0_xxxxxyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 73);

    auto g_0_xxxxxyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 74);

    auto g_0_xxxxxyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 75);

    auto g_0_xxxxxyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 76);

    auto g_0_xxxxxyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 77);

    auto g_0_xxxxxyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 78);

    auto g_0_xxxxxyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 79);

    auto g_0_xxxxxyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 80);

    auto g_0_xxxxxyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 81);

    auto g_0_xxxxxyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 82);

    auto g_0_xxxxxyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 83);

    auto g_0_xxxxxzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 105);

    auto g_0_xxxxxzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 106);

    auto g_0_xxxxxzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 107);

    auto g_0_xxxxxzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 108);

    auto g_0_xxxxxzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 109);

    auto g_0_xxxxxzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 110);

    auto g_0_xxxxxzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 111);

    auto g_0_xxxxxzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 112);

    auto g_0_xxxxxzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 113);

    auto g_0_xxxxxzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 114);

    auto g_0_xxxxxzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 115);

    auto g_0_xxxxxzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 116);

    auto g_0_xxxxxzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 117);

    auto g_0_xxxxxzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 118);

    auto g_0_xxxxxzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 119);

    auto g_0_xxxxxzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 120);

    auto g_0_xxxxxzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 121);

    auto g_0_xxxxxzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 122);

    auto g_0_xxxxxzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 123);

    auto g_0_xxxxxzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 124);

    auto g_0_xxxxxzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 125);

    auto g_0_xxxxyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 126);

    auto g_0_xxxxyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 127);

    auto g_0_xxxxyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 128);

    auto g_0_xxxxyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 129);

    auto g_0_xxxxyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 130);

    auto g_0_xxxxyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 131);

    auto g_0_xxxxyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 132);

    auto g_0_xxxxyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 133);

    auto g_0_xxxxyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 134);

    auto g_0_xxxxyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 135);

    auto g_0_xxxxyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 136);

    auto g_0_xxxxyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 137);

    auto g_0_xxxxyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 138);

    auto g_0_xxxxyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 139);

    auto g_0_xxxxyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 140);

    auto g_0_xxxxyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 141);

    auto g_0_xxxxyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 142);

    auto g_0_xxxxyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 143);

    auto g_0_xxxxyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 144);

    auto g_0_xxxxyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 145);

    auto g_0_xxxxyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 146);

    auto g_0_xxxxzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 189);

    auto g_0_xxxxzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 190);

    auto g_0_xxxxzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 191);

    auto g_0_xxxxzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 192);

    auto g_0_xxxxzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 193);

    auto g_0_xxxxzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 194);

    auto g_0_xxxxzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 195);

    auto g_0_xxxxzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 196);

    auto g_0_xxxxzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 197);

    auto g_0_xxxxzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 198);

    auto g_0_xxxxzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 199);

    auto g_0_xxxxzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 200);

    auto g_0_xxxxzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 201);

    auto g_0_xxxxzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 202);

    auto g_0_xxxxzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 203);

    auto g_0_xxxxzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 204);

    auto g_0_xxxxzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 205);

    auto g_0_xxxxzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 206);

    auto g_0_xxxxzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 207);

    auto g_0_xxxxzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 208);

    auto g_0_xxxxzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 209);

    auto g_0_xxxyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 210);

    auto g_0_xxxyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 211);

    auto g_0_xxxyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 212);

    auto g_0_xxxyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 213);

    auto g_0_xxxyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 214);

    auto g_0_xxxyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 215);

    auto g_0_xxxyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 216);

    auto g_0_xxxyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 217);

    auto g_0_xxxyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 218);

    auto g_0_xxxyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 219);

    auto g_0_xxxyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 220);

    auto g_0_xxxyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 221);

    auto g_0_xxxyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 222);

    auto g_0_xxxyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 223);

    auto g_0_xxxyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 224);

    auto g_0_xxxyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 225);

    auto g_0_xxxyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 226);

    auto g_0_xxxyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 227);

    auto g_0_xxxyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 228);

    auto g_0_xxxyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 229);

    auto g_0_xxxyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 230);

    auto g_0_xxxyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 256);

    auto g_0_xxxyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 259);

    auto g_0_xxxyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 260);

    auto g_0_xxxyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 263);

    auto g_0_xxxyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 264);

    auto g_0_xxxyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 265);

    auto g_0_xxxyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 268);

    auto g_0_xxxyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 269);

    auto g_0_xxxyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 270);

    auto g_0_xxxyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 271);

    auto g_0_xxxzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 294);

    auto g_0_xxxzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 295);

    auto g_0_xxxzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 296);

    auto g_0_xxxzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 297);

    auto g_0_xxxzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 298);

    auto g_0_xxxzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 299);

    auto g_0_xxxzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 300);

    auto g_0_xxxzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 301);

    auto g_0_xxxzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 302);

    auto g_0_xxxzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 303);

    auto g_0_xxxzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 304);

    auto g_0_xxxzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 305);

    auto g_0_xxxzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 306);

    auto g_0_xxxzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 307);

    auto g_0_xxxzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 308);

    auto g_0_xxxzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 309);

    auto g_0_xxxzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 310);

    auto g_0_xxxzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 311);

    auto g_0_xxxzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 312);

    auto g_0_xxxzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 313);

    auto g_0_xxxzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 314);

    auto g_0_xxyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 315);

    auto g_0_xxyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 316);

    auto g_0_xxyyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 317);

    auto g_0_xxyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 318);

    auto g_0_xxyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 319);

    auto g_0_xxyyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 320);

    auto g_0_xxyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 321);

    auto g_0_xxyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 322);

    auto g_0_xxyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 323);

    auto g_0_xxyyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 324);

    auto g_0_xxyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 325);

    auto g_0_xxyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 326);

    auto g_0_xxyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 327);

    auto g_0_xxyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 328);

    auto g_0_xxyyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 329);

    auto g_0_xxyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 330);

    auto g_0_xxyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 331);

    auto g_0_xxyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 332);

    auto g_0_xxyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 333);

    auto g_0_xxyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 334);

    auto g_0_xxyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 335);

    auto g_0_xxyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 361);

    auto g_0_xxyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 364);

    auto g_0_xxyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 365);

    auto g_0_xxyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 368);

    auto g_0_xxyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 369);

    auto g_0_xxyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 370);

    auto g_0_xxyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 373);

    auto g_0_xxyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 374);

    auto g_0_xxyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 375);

    auto g_0_xxyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 376);

    auto g_0_xxyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 382);

    auto g_0_xxyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 385);

    auto g_0_xxyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 386);

    auto g_0_xxyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 389);

    auto g_0_xxyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 390);

    auto g_0_xxyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 391);

    auto g_0_xxyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 394);

    auto g_0_xxyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 395);

    auto g_0_xxyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 396);

    auto g_0_xxyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 397);

    auto g_0_xxzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 420);

    auto g_0_xxzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 421);

    auto g_0_xxzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 422);

    auto g_0_xxzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 423);

    auto g_0_xxzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 424);

    auto g_0_xxzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 425);

    auto g_0_xxzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 426);

    auto g_0_xxzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 427);

    auto g_0_xxzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 428);

    auto g_0_xxzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 429);

    auto g_0_xxzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 430);

    auto g_0_xxzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 431);

    auto g_0_xxzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 432);

    auto g_0_xxzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 433);

    auto g_0_xxzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 434);

    auto g_0_xxzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 435);

    auto g_0_xxzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 436);

    auto g_0_xxzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 437);

    auto g_0_xxzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 438);

    auto g_0_xxzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 439);

    auto g_0_xxzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 440);

    auto g_0_xyyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 442);

    auto g_0_xyyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 444);

    auto g_0_xyyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 445);

    auto g_0_xyyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 447);

    auto g_0_xyyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 448);

    auto g_0_xyyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 449);

    auto g_0_xyyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 451);

    auto g_0_xyyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 452);

    auto g_0_xyyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 453);

    auto g_0_xyyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 454);

    auto g_0_xyyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 456);

    auto g_0_xyyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 457);

    auto g_0_xyyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 458);

    auto g_0_xyyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 459);

    auto g_0_xyyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 460);

    auto g_0_xyyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 487);

    auto g_0_xyyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 490);

    auto g_0_xyyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 491);

    auto g_0_xyyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 494);

    auto g_0_xyyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 495);

    auto g_0_xyyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 496);

    auto g_0_xyyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 499);

    auto g_0_xyyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 500);

    auto g_0_xyyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 501);

    auto g_0_xyyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 502);

    auto g_0_xyyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 508);

    auto g_0_xyyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 511);

    auto g_0_xyyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 512);

    auto g_0_xyyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 515);

    auto g_0_xyyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 516);

    auto g_0_xyyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 517);

    auto g_0_xyyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 520);

    auto g_0_xyyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 521);

    auto g_0_xyyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 522);

    auto g_0_xyyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 523);

    auto g_0_xyyzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 529);

    auto g_0_xyyzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 532);

    auto g_0_xyyzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 533);

    auto g_0_xyyzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 536);

    auto g_0_xyyzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 537);

    auto g_0_xyyzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 538);

    auto g_0_xyyzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 541);

    auto g_0_xyyzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 542);

    auto g_0_xyyzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 543);

    auto g_0_xyyzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 544);

    auto g_0_xzzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 569);

    auto g_0_xzzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 571);

    auto g_0_xzzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 572);

    auto g_0_xzzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 574);

    auto g_0_xzzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 575);

    auto g_0_xzzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 576);

    auto g_0_xzzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 578);

    auto g_0_xzzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 579);

    auto g_0_xzzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 580);

    auto g_0_xzzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 581);

    auto g_0_xzzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 583);

    auto g_0_xzzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 584);

    auto g_0_xzzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 585);

    auto g_0_xzzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 586);

    auto g_0_xzzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 587);

    auto g_0_yyyyyyy_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 588);

    auto g_0_yyyyyyy_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 589);

    auto g_0_yyyyyyy_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 590);

    auto g_0_yyyyyyy_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 591);

    auto g_0_yyyyyyy_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 592);

    auto g_0_yyyyyyy_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 593);

    auto g_0_yyyyyyy_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 594);

    auto g_0_yyyyyyy_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 595);

    auto g_0_yyyyyyy_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 596);

    auto g_0_yyyyyyy_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 597);

    auto g_0_yyyyyyy_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 598);

    auto g_0_yyyyyyy_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 599);

    auto g_0_yyyyyyy_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 600);

    auto g_0_yyyyyyy_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 601);

    auto g_0_yyyyyyy_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 602);

    auto g_0_yyyyyyy_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 603);

    auto g_0_yyyyyyy_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 604);

    auto g_0_yyyyyyy_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 605);

    auto g_0_yyyyyyy_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 606);

    auto g_0_yyyyyyy_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 607);

    auto g_0_yyyyyyy_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 608);

    auto g_0_yyyyyyz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 611);

    auto g_0_yyyyyyz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 613);

    auto g_0_yyyyyyz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 614);

    auto g_0_yyyyyyz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 616);

    auto g_0_yyyyyyz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 617);

    auto g_0_yyyyyyz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 618);

    auto g_0_yyyyyyz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 620);

    auto g_0_yyyyyyz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 621);

    auto g_0_yyyyyyz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 622);

    auto g_0_yyyyyyz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 623);

    auto g_0_yyyyyyz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 625);

    auto g_0_yyyyyyz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 626);

    auto g_0_yyyyyyz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 627);

    auto g_0_yyyyyyz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 628);

    auto g_0_yyyyyyz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 629);

    auto g_0_yyyyyzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 630);

    auto g_0_yyyyyzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 631);

    auto g_0_yyyyyzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 632);

    auto g_0_yyyyyzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 633);

    auto g_0_yyyyyzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 634);

    auto g_0_yyyyyzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 635);

    auto g_0_yyyyyzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 636);

    auto g_0_yyyyyzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 637);

    auto g_0_yyyyyzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 638);

    auto g_0_yyyyyzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 639);

    auto g_0_yyyyyzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 640);

    auto g_0_yyyyyzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 641);

    auto g_0_yyyyyzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 642);

    auto g_0_yyyyyzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 643);

    auto g_0_yyyyyzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 644);

    auto g_0_yyyyyzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 645);

    auto g_0_yyyyyzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 646);

    auto g_0_yyyyyzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 647);

    auto g_0_yyyyyzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 648);

    auto g_0_yyyyyzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 649);

    auto g_0_yyyyyzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 650);

    auto g_0_yyyyzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 651);

    auto g_0_yyyyzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 652);

    auto g_0_yyyyzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 653);

    auto g_0_yyyyzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 654);

    auto g_0_yyyyzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 655);

    auto g_0_yyyyzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 656);

    auto g_0_yyyyzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 657);

    auto g_0_yyyyzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 658);

    auto g_0_yyyyzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 659);

    auto g_0_yyyyzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 660);

    auto g_0_yyyyzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 661);

    auto g_0_yyyyzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 662);

    auto g_0_yyyyzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 663);

    auto g_0_yyyyzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 664);

    auto g_0_yyyyzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 665);

    auto g_0_yyyyzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 666);

    auto g_0_yyyyzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 667);

    auto g_0_yyyyzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 668);

    auto g_0_yyyyzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 669);

    auto g_0_yyyyzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 670);

    auto g_0_yyyyzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 671);

    auto g_0_yyyzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 672);

    auto g_0_yyyzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 673);

    auto g_0_yyyzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 674);

    auto g_0_yyyzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 675);

    auto g_0_yyyzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 676);

    auto g_0_yyyzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 677);

    auto g_0_yyyzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 678);

    auto g_0_yyyzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 679);

    auto g_0_yyyzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 680);

    auto g_0_yyyzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 681);

    auto g_0_yyyzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 682);

    auto g_0_yyyzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 683);

    auto g_0_yyyzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 684);

    auto g_0_yyyzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 685);

    auto g_0_yyyzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 686);

    auto g_0_yyyzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 687);

    auto g_0_yyyzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 688);

    auto g_0_yyyzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 689);

    auto g_0_yyyzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 690);

    auto g_0_yyyzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 691);

    auto g_0_yyyzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 692);

    auto g_0_yyzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 693);

    auto g_0_yyzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 694);

    auto g_0_yyzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 695);

    auto g_0_yyzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 696);

    auto g_0_yyzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 697);

    auto g_0_yyzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 698);

    auto g_0_yyzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 699);

    auto g_0_yyzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 700);

    auto g_0_yyzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 701);

    auto g_0_yyzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 702);

    auto g_0_yyzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 703);

    auto g_0_yyzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 704);

    auto g_0_yyzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 705);

    auto g_0_yyzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 706);

    auto g_0_yyzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 707);

    auto g_0_yyzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 708);

    auto g_0_yyzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 709);

    auto g_0_yyzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 710);

    auto g_0_yyzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 711);

    auto g_0_yyzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 712);

    auto g_0_yyzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 713);

    auto g_0_yzzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 715);

    auto g_0_yzzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 716);

    auto g_0_yzzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 717);

    auto g_0_yzzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 718);

    auto g_0_yzzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 719);

    auto g_0_yzzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 720);

    auto g_0_yzzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 721);

    auto g_0_yzzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 722);

    auto g_0_yzzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 723);

    auto g_0_yzzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 724);

    auto g_0_yzzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 725);

    auto g_0_yzzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 726);

    auto g_0_yzzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 727);

    auto g_0_yzzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 728);

    auto g_0_yzzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 729);

    auto g_0_yzzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 730);

    auto g_0_yzzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 731);

    auto g_0_yzzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 732);

    auto g_0_yzzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 733);

    auto g_0_yzzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 734);

    auto g_0_zzzzzzz_0_xxxxx_1 = pbuffer.data(idx_eri_1_sksh + 735);

    auto g_0_zzzzzzz_0_xxxxy_1 = pbuffer.data(idx_eri_1_sksh + 736);

    auto g_0_zzzzzzz_0_xxxxz_1 = pbuffer.data(idx_eri_1_sksh + 737);

    auto g_0_zzzzzzz_0_xxxyy_1 = pbuffer.data(idx_eri_1_sksh + 738);

    auto g_0_zzzzzzz_0_xxxyz_1 = pbuffer.data(idx_eri_1_sksh + 739);

    auto g_0_zzzzzzz_0_xxxzz_1 = pbuffer.data(idx_eri_1_sksh + 740);

    auto g_0_zzzzzzz_0_xxyyy_1 = pbuffer.data(idx_eri_1_sksh + 741);

    auto g_0_zzzzzzz_0_xxyyz_1 = pbuffer.data(idx_eri_1_sksh + 742);

    auto g_0_zzzzzzz_0_xxyzz_1 = pbuffer.data(idx_eri_1_sksh + 743);

    auto g_0_zzzzzzz_0_xxzzz_1 = pbuffer.data(idx_eri_1_sksh + 744);

    auto g_0_zzzzzzz_0_xyyyy_1 = pbuffer.data(idx_eri_1_sksh + 745);

    auto g_0_zzzzzzz_0_xyyyz_1 = pbuffer.data(idx_eri_1_sksh + 746);

    auto g_0_zzzzzzz_0_xyyzz_1 = pbuffer.data(idx_eri_1_sksh + 747);

    auto g_0_zzzzzzz_0_xyzzz_1 = pbuffer.data(idx_eri_1_sksh + 748);

    auto g_0_zzzzzzz_0_xzzzz_1 = pbuffer.data(idx_eri_1_sksh + 749);

    auto g_0_zzzzzzz_0_yyyyy_1 = pbuffer.data(idx_eri_1_sksh + 750);

    auto g_0_zzzzzzz_0_yyyyz_1 = pbuffer.data(idx_eri_1_sksh + 751);

    auto g_0_zzzzzzz_0_yyyzz_1 = pbuffer.data(idx_eri_1_sksh + 752);

    auto g_0_zzzzzzz_0_yyzzz_1 = pbuffer.data(idx_eri_1_sksh + 753);

    auto g_0_zzzzzzz_0_yzzzz_1 = pbuffer.data(idx_eri_1_sksh + 754);

    auto g_0_zzzzzzz_0_zzzzz_1 = pbuffer.data(idx_eri_1_sksh + 755);

    /// Set up components of auxilary buffer : SKSI

    auto g_0_xxxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi);

    auto g_0_xxxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 1);

    auto g_0_xxxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 2);

    auto g_0_xxxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 3);

    auto g_0_xxxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 4);

    auto g_0_xxxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 5);

    auto g_0_xxxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 6);

    auto g_0_xxxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 7);

    auto g_0_xxxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 8);

    auto g_0_xxxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 9);

    auto g_0_xxxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 10);

    auto g_0_xxxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 11);

    auto g_0_xxxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 12);

    auto g_0_xxxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 13);

    auto g_0_xxxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 14);

    auto g_0_xxxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 15);

    auto g_0_xxxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 16);

    auto g_0_xxxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 17);

    auto g_0_xxxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 18);

    auto g_0_xxxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 19);

    auto g_0_xxxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 20);

    auto g_0_xxxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 21);

    auto g_0_xxxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 22);

    auto g_0_xxxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 23);

    auto g_0_xxxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 24);

    auto g_0_xxxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 25);

    auto g_0_xxxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 26);

    auto g_0_xxxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 27);

    auto g_0_xxxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 28);

    auto g_0_xxxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 29);

    auto g_0_xxxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 30);

    auto g_0_xxxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 31);

    auto g_0_xxxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 33);

    auto g_0_xxxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 34);

    auto g_0_xxxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 37);

    auto g_0_xxxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 38);

    auto g_0_xxxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 42);

    auto g_0_xxxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 43);

    auto g_0_xxxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 48);

    auto g_0_xxxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 49);

    auto g_0_xxxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 56);

    auto g_0_xxxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 57);

    auto g_0_xxxxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 58);

    auto g_0_xxxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 59);

    auto g_0_xxxxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 60);

    auto g_0_xxxxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 61);

    auto g_0_xxxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 62);

    auto g_0_xxxxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 63);

    auto g_0_xxxxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 64);

    auto g_0_xxxxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 65);

    auto g_0_xxxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 66);

    auto g_0_xxxxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 67);

    auto g_0_xxxxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 68);

    auto g_0_xxxxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 69);

    auto g_0_xxxxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 70);

    auto g_0_xxxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 71);

    auto g_0_xxxxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 72);

    auto g_0_xxxxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 73);

    auto g_0_xxxxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 74);

    auto g_0_xxxxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 75);

    auto g_0_xxxxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 76);

    auto g_0_xxxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 78);

    auto g_0_xxxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 79);

    auto g_0_xxxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 80);

    auto g_0_xxxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 81);

    auto g_0_xxxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 82);

    auto g_0_xxxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 83);

    auto g_0_xxxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 84);

    auto g_0_xxxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 85);

    auto g_0_xxxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 86);

    auto g_0_xxxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 87);

    auto g_0_xxxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 88);

    auto g_0_xxxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 89);

    auto g_0_xxxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 90);

    auto g_0_xxxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 91);

    auto g_0_xxxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 92);

    auto g_0_xxxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 93);

    auto g_0_xxxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 94);

    auto g_0_xxxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 95);

    auto g_0_xxxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 96);

    auto g_0_xxxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 97);

    auto g_0_xxxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 98);

    auto g_0_xxxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 99);

    auto g_0_xxxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 100);

    auto g_0_xxxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 101);

    auto g_0_xxxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 102);

    auto g_0_xxxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 103);

    auto g_0_xxxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 104);

    auto g_0_xxxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 105);

    auto g_0_xxxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 106);

    auto g_0_xxxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 107);

    auto g_0_xxxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 108);

    auto g_0_xxxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 109);

    auto g_0_xxxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 110);

    auto g_0_xxxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 111);

    auto g_0_xxxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 140);

    auto g_0_xxxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 141);

    auto g_0_xxxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 142);

    auto g_0_xxxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 143);

    auto g_0_xxxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 144);

    auto g_0_xxxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 145);

    auto g_0_xxxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 146);

    auto g_0_xxxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 147);

    auto g_0_xxxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 148);

    auto g_0_xxxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 149);

    auto g_0_xxxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 150);

    auto g_0_xxxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 151);

    auto g_0_xxxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 152);

    auto g_0_xxxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 153);

    auto g_0_xxxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 154);

    auto g_0_xxxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 155);

    auto g_0_xxxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 156);

    auto g_0_xxxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 157);

    auto g_0_xxxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 158);

    auto g_0_xxxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 159);

    auto g_0_xxxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 160);

    auto g_0_xxxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 161);

    auto g_0_xxxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 162);

    auto g_0_xxxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 163);

    auto g_0_xxxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 164);

    auto g_0_xxxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 165);

    auto g_0_xxxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 166);

    auto g_0_xxxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 167);

    auto g_0_xxxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 168);

    auto g_0_xxxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 169);

    auto g_0_xxxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 170);

    auto g_0_xxxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 171);

    auto g_0_xxxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 172);

    auto g_0_xxxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 173);

    auto g_0_xxxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 174);

    auto g_0_xxxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 175);

    auto g_0_xxxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 176);

    auto g_0_xxxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 177);

    auto g_0_xxxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 178);

    auto g_0_xxxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 179);

    auto g_0_xxxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 180);

    auto g_0_xxxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 181);

    auto g_0_xxxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 182);

    auto g_0_xxxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 183);

    auto g_0_xxxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 184);

    auto g_0_xxxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 185);

    auto g_0_xxxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 186);

    auto g_0_xxxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 187);

    auto g_0_xxxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 188);

    auto g_0_xxxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 189);

    auto g_0_xxxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 190);

    auto g_0_xxxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 191);

    auto g_0_xxxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 192);

    auto g_0_xxxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 193);

    auto g_0_xxxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 194);

    auto g_0_xxxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 195);

    auto g_0_xxxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 197);

    auto g_0_xxxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 199);

    auto g_0_xxxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 202);

    auto g_0_xxxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 206);

    auto g_0_xxxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 211);

    auto g_0_xxxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 224);

    auto g_0_xxxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 226);

    auto g_0_xxxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 229);

    auto g_0_xxxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 233);

    auto g_0_xxxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 238);

    auto g_0_xxxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 244);

    auto g_0_xxxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 252);

    auto g_0_xxxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 253);

    auto g_0_xxxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 254);

    auto g_0_xxxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 255);

    auto g_0_xxxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 256);

    auto g_0_xxxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 257);

    auto g_0_xxxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 258);

    auto g_0_xxxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 259);

    auto g_0_xxxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 260);

    auto g_0_xxxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 261);

    auto g_0_xxxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 262);

    auto g_0_xxxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 263);

    auto g_0_xxxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 264);

    auto g_0_xxxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 265);

    auto g_0_xxxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 266);

    auto g_0_xxxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 267);

    auto g_0_xxxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 268);

    auto g_0_xxxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 269);

    auto g_0_xxxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 270);

    auto g_0_xxxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 271);

    auto g_0_xxxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 272);

    auto g_0_xxxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 273);

    auto g_0_xxxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 274);

    auto g_0_xxxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 275);

    auto g_0_xxxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 276);

    auto g_0_xxxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 277);

    auto g_0_xxxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 278);

    auto g_0_xxxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 279);

    auto g_0_xxxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 280);

    auto g_0_xxxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 281);

    auto g_0_xxxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 282);

    auto g_0_xxxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 283);

    auto g_0_xxxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 284);

    auto g_0_xxxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 285);

    auto g_0_xxxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 286);

    auto g_0_xxxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 287);

    auto g_0_xxxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 288);

    auto g_0_xxxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 289);

    auto g_0_xxxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 290);

    auto g_0_xxxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 291);

    auto g_0_xxxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 292);

    auto g_0_xxxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 293);

    auto g_0_xxxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 294);

    auto g_0_xxxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 295);

    auto g_0_xxxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 296);

    auto g_0_xxxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 297);

    auto g_0_xxxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 298);

    auto g_0_xxxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 299);

    auto g_0_xxxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 300);

    auto g_0_xxxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 301);

    auto g_0_xxxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 302);

    auto g_0_xxxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 303);

    auto g_0_xxxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 304);

    auto g_0_xxxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 305);

    auto g_0_xxxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 306);

    auto g_0_xxxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 307);

    auto g_0_xxxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 309);

    auto g_0_xxxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 311);

    auto g_0_xxxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 314);

    auto g_0_xxxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 318);

    auto g_0_xxxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 323);

    auto g_0_xxxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 336);

    auto g_0_xxxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 337);

    auto g_0_xxxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 338);

    auto g_0_xxxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 339);

    auto g_0_xxxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 340);

    auto g_0_xxxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 341);

    auto g_0_xxxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 342);

    auto g_0_xxxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 343);

    auto g_0_xxxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 344);

    auto g_0_xxxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 345);

    auto g_0_xxxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 346);

    auto g_0_xxxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 347);

    auto g_0_xxxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 348);

    auto g_0_xxxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 349);

    auto g_0_xxxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 350);

    auto g_0_xxxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 351);

    auto g_0_xxxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 352);

    auto g_0_xxxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 353);

    auto g_0_xxxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 354);

    auto g_0_xxxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 355);

    auto g_0_xxxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 356);

    auto g_0_xxxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 357);

    auto g_0_xxxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 358);

    auto g_0_xxxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 359);

    auto g_0_xxxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 360);

    auto g_0_xxxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 361);

    auto g_0_xxxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 362);

    auto g_0_xxxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 363);

    auto g_0_xxxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 364);

    auto g_0_xxxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 366);

    auto g_0_xxxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 369);

    auto g_0_xxxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 373);

    auto g_0_xxxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 378);

    auto g_0_xxxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 384);

    auto g_0_xxxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 392);

    auto g_0_xxxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 393);

    auto g_0_xxxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 394);

    auto g_0_xxxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 395);

    auto g_0_xxxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 396);

    auto g_0_xxxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 397);

    auto g_0_xxxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 398);

    auto g_0_xxxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 399);

    auto g_0_xxxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 400);

    auto g_0_xxxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 401);

    auto g_0_xxxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 402);

    auto g_0_xxxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 403);

    auto g_0_xxxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 404);

    auto g_0_xxxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 405);

    auto g_0_xxxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 406);

    auto g_0_xxxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 407);

    auto g_0_xxxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 408);

    auto g_0_xxxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 409);

    auto g_0_xxxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 410);

    auto g_0_xxxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 411);

    auto g_0_xxxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 412);

    auto g_0_xxxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 413);

    auto g_0_xxxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 414);

    auto g_0_xxxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 415);

    auto g_0_xxxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 416);

    auto g_0_xxxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 417);

    auto g_0_xxxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 418);

    auto g_0_xxxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 419);

    auto g_0_xxyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 420);

    auto g_0_xxyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 421);

    auto g_0_xxyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 422);

    auto g_0_xxyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 423);

    auto g_0_xxyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 424);

    auto g_0_xxyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 425);

    auto g_0_xxyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 426);

    auto g_0_xxyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 427);

    auto g_0_xxyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 428);

    auto g_0_xxyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 429);

    auto g_0_xxyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 430);

    auto g_0_xxyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 431);

    auto g_0_xxyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 432);

    auto g_0_xxyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 433);

    auto g_0_xxyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 434);

    auto g_0_xxyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 435);

    auto g_0_xxyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 436);

    auto g_0_xxyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 437);

    auto g_0_xxyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 438);

    auto g_0_xxyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 439);

    auto g_0_xxyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 440);

    auto g_0_xxyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 441);

    auto g_0_xxyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 442);

    auto g_0_xxyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 443);

    auto g_0_xxyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 444);

    auto g_0_xxyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 445);

    auto g_0_xxyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 446);

    auto g_0_xxyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 447);

    auto g_0_xxyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 449);

    auto g_0_xxyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 451);

    auto g_0_xxyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 454);

    auto g_0_xxyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 458);

    auto g_0_xxyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 463);

    auto g_0_xxyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 476);

    auto g_0_xxyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 477);

    auto g_0_xxyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 478);

    auto g_0_xxyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 479);

    auto g_0_xxyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 480);

    auto g_0_xxyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 481);

    auto g_0_xxyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 482);

    auto g_0_xxyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 483);

    auto g_0_xxyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 484);

    auto g_0_xxyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 485);

    auto g_0_xxyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 486);

    auto g_0_xxyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 487);

    auto g_0_xxyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 488);

    auto g_0_xxyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 489);

    auto g_0_xxyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 490);

    auto g_0_xxyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 491);

    auto g_0_xxyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 492);

    auto g_0_xxyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 493);

    auto g_0_xxyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 494);

    auto g_0_xxyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 495);

    auto g_0_xxyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 496);

    auto g_0_xxyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 497);

    auto g_0_xxyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 498);

    auto g_0_xxyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 499);

    auto g_0_xxyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 500);

    auto g_0_xxyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 501);

    auto g_0_xxyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 502);

    auto g_0_xxyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 503);

    auto g_0_xxyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 504);

    auto g_0_xxyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 505);

    auto g_0_xxyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 506);

    auto g_0_xxyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 507);

    auto g_0_xxyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 508);

    auto g_0_xxyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 509);

    auto g_0_xxyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 510);

    auto g_0_xxyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 511);

    auto g_0_xxyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 512);

    auto g_0_xxyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 513);

    auto g_0_xxyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 514);

    auto g_0_xxyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 515);

    auto g_0_xxyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 516);

    auto g_0_xxyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 517);

    auto g_0_xxyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 518);

    auto g_0_xxyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 519);

    auto g_0_xxyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 520);

    auto g_0_xxyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 521);

    auto g_0_xxyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 522);

    auto g_0_xxyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 523);

    auto g_0_xxyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 524);

    auto g_0_xxyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 525);

    auto g_0_xxyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 526);

    auto g_0_xxyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 527);

    auto g_0_xxyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 528);

    auto g_0_xxyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 529);

    auto g_0_xxyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 530);

    auto g_0_xxyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 531);

    auto g_0_xxyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 532);

    auto g_0_xxyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 534);

    auto g_0_xxyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 537);

    auto g_0_xxyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 541);

    auto g_0_xxyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 546);

    auto g_0_xxyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 552);

    auto g_0_xxzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 560);

    auto g_0_xxzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 561);

    auto g_0_xxzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 562);

    auto g_0_xxzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 563);

    auto g_0_xxzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 564);

    auto g_0_xxzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 565);

    auto g_0_xxzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 566);

    auto g_0_xxzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 567);

    auto g_0_xxzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 568);

    auto g_0_xxzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 569);

    auto g_0_xxzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 570);

    auto g_0_xxzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 571);

    auto g_0_xxzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 572);

    auto g_0_xxzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 573);

    auto g_0_xxzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 574);

    auto g_0_xxzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 575);

    auto g_0_xxzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 576);

    auto g_0_xxzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 577);

    auto g_0_xxzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 578);

    auto g_0_xxzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 579);

    auto g_0_xxzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 580);

    auto g_0_xxzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 581);

    auto g_0_xxzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 582);

    auto g_0_xxzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 583);

    auto g_0_xxzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 584);

    auto g_0_xxzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 585);

    auto g_0_xxzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 586);

    auto g_0_xxzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 587);

    auto g_0_xyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 588);

    auto g_0_xyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 589);

    auto g_0_xyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 591);

    auto g_0_xyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 592);

    auto g_0_xyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 594);

    auto g_0_xyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 595);

    auto g_0_xyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 596);

    auto g_0_xyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 598);

    auto g_0_xyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 599);

    auto g_0_xyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 600);

    auto g_0_xyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 601);

    auto g_0_xyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 603);

    auto g_0_xyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 604);

    auto g_0_xyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 605);

    auto g_0_xyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 606);

    auto g_0_xyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 607);

    auto g_0_xyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 609);

    auto g_0_xyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 610);

    auto g_0_xyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 611);

    auto g_0_xyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 612);

    auto g_0_xyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 613);

    auto g_0_xyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 614);

    auto g_0_xyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 615);

    auto g_0_xyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 648);

    auto g_0_xyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 651);

    auto g_0_xyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 652);

    auto g_0_xyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 655);

    auto g_0_xyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 656);

    auto g_0_xyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 657);

    auto g_0_xyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 660);

    auto g_0_xyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 661);

    auto g_0_xyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 662);

    auto g_0_xyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 663);

    auto g_0_xyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 665);

    auto g_0_xyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 666);

    auto g_0_xyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 667);

    auto g_0_xyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 668);

    auto g_0_xyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 669);

    auto g_0_xyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 670);

    auto g_0_xyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 671);

    auto g_0_xyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 676);

    auto g_0_xyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 679);

    auto g_0_xyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 680);

    auto g_0_xyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 683);

    auto g_0_xyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 684);

    auto g_0_xyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 685);

    auto g_0_xyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 688);

    auto g_0_xyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 689);

    auto g_0_xyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 690);

    auto g_0_xyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 691);

    auto g_0_xyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 693);

    auto g_0_xyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 694);

    auto g_0_xyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 695);

    auto g_0_xyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 696);

    auto g_0_xyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 697);

    auto g_0_xyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 698);

    auto g_0_xyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 699);

    auto g_0_xyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 704);

    auto g_0_xyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 707);

    auto g_0_xyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 708);

    auto g_0_xyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 711);

    auto g_0_xyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 712);

    auto g_0_xyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 713);

    auto g_0_xyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 716);

    auto g_0_xyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 717);

    auto g_0_xyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 718);

    auto g_0_xyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 719);

    auto g_0_xyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 721);

    auto g_0_xyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 722);

    auto g_0_xyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 723);

    auto g_0_xyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 724);

    auto g_0_xyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 725);

    auto g_0_xyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 726);

    auto g_0_xyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 727);

    auto g_0_xzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 756);

    auto g_0_xzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 758);

    auto g_0_xzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 760);

    auto g_0_xzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 761);

    auto g_0_xzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 763);

    auto g_0_xzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 764);

    auto g_0_xzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 765);

    auto g_0_xzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 767);

    auto g_0_xzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 768);

    auto g_0_xzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 769);

    auto g_0_xzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 770);

    auto g_0_xzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 772);

    auto g_0_xzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 773);

    auto g_0_xzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 774);

    auto g_0_xzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 775);

    auto g_0_xzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 776);

    auto g_0_xzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 777);

    auto g_0_xzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 778);

    auto g_0_xzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 779);

    auto g_0_xzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 780);

    auto g_0_xzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 781);

    auto g_0_xzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 782);

    auto g_0_xzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 783);

    auto g_0_yyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 784);

    auto g_0_yyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 785);

    auto g_0_yyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 786);

    auto g_0_yyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 787);

    auto g_0_yyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 788);

    auto g_0_yyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 789);

    auto g_0_yyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 790);

    auto g_0_yyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 791);

    auto g_0_yyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 792);

    auto g_0_yyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 793);

    auto g_0_yyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 794);

    auto g_0_yyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 795);

    auto g_0_yyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 796);

    auto g_0_yyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 797);

    auto g_0_yyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 798);

    auto g_0_yyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 799);

    auto g_0_yyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 800);

    auto g_0_yyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 801);

    auto g_0_yyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 802);

    auto g_0_yyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 803);

    auto g_0_yyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 804);

    auto g_0_yyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 805);

    auto g_0_yyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 806);

    auto g_0_yyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 807);

    auto g_0_yyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 808);

    auto g_0_yyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 809);

    auto g_0_yyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 810);

    auto g_0_yyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 811);

    auto g_0_yyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 813);

    auto g_0_yyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 814);

    auto g_0_yyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 815);

    auto g_0_yyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 816);

    auto g_0_yyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 817);

    auto g_0_yyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 818);

    auto g_0_yyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 819);

    auto g_0_yyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 820);

    auto g_0_yyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 821);

    auto g_0_yyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 822);

    auto g_0_yyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 823);

    auto g_0_yyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 824);

    auto g_0_yyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 825);

    auto g_0_yyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 826);

    auto g_0_yyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 827);

    auto g_0_yyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 828);

    auto g_0_yyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 829);

    auto g_0_yyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 830);

    auto g_0_yyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 831);

    auto g_0_yyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 832);

    auto g_0_yyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 833);

    auto g_0_yyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 834);

    auto g_0_yyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 835);

    auto g_0_yyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 836);

    auto g_0_yyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 837);

    auto g_0_yyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 838);

    auto g_0_yyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 839);

    auto g_0_yyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 840);

    auto g_0_yyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 841);

    auto g_0_yyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 842);

    auto g_0_yyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 843);

    auto g_0_yyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 844);

    auto g_0_yyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 845);

    auto g_0_yyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 846);

    auto g_0_yyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 847);

    auto g_0_yyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 848);

    auto g_0_yyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 849);

    auto g_0_yyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 850);

    auto g_0_yyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 851);

    auto g_0_yyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 852);

    auto g_0_yyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 853);

    auto g_0_yyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 854);

    auto g_0_yyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 855);

    auto g_0_yyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 856);

    auto g_0_yyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 857);

    auto g_0_yyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 858);

    auto g_0_yyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 859);

    auto g_0_yyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 860);

    auto g_0_yyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 861);

    auto g_0_yyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 862);

    auto g_0_yyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 863);

    auto g_0_yyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 864);

    auto g_0_yyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 865);

    auto g_0_yyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 866);

    auto g_0_yyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 867);

    auto g_0_yyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 868);

    auto g_0_yyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 869);

    auto g_0_yyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 870);

    auto g_0_yyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 871);

    auto g_0_yyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 872);

    auto g_0_yyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 873);

    auto g_0_yyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 874);

    auto g_0_yyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 875);

    auto g_0_yyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 876);

    auto g_0_yyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 877);

    auto g_0_yyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 878);

    auto g_0_yyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 879);

    auto g_0_yyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 880);

    auto g_0_yyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 881);

    auto g_0_yyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 882);

    auto g_0_yyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 883);

    auto g_0_yyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 884);

    auto g_0_yyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 885);

    auto g_0_yyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 886);

    auto g_0_yyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 887);

    auto g_0_yyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 888);

    auto g_0_yyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 889);

    auto g_0_yyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 890);

    auto g_0_yyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 891);

    auto g_0_yyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 892);

    auto g_0_yyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 893);

    auto g_0_yyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 894);

    auto g_0_yyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 895);

    auto g_0_yyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 896);

    auto g_0_yyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 897);

    auto g_0_yyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 898);

    auto g_0_yyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 899);

    auto g_0_yyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 900);

    auto g_0_yyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 901);

    auto g_0_yyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 902);

    auto g_0_yyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 903);

    auto g_0_yyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 904);

    auto g_0_yyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 905);

    auto g_0_yyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 906);

    auto g_0_yyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 907);

    auto g_0_yyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 908);

    auto g_0_yyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 909);

    auto g_0_yyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 910);

    auto g_0_yyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 911);

    auto g_0_yyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 912);

    auto g_0_yyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 913);

    auto g_0_yyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 914);

    auto g_0_yyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 915);

    auto g_0_yyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 916);

    auto g_0_yyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 917);

    auto g_0_yyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 918);

    auto g_0_yyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 919);

    auto g_0_yyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 920);

    auto g_0_yyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 921);

    auto g_0_yyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 922);

    auto g_0_yyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 923);

    auto g_0_yyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 924);

    auto g_0_yyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 925);

    auto g_0_yyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 926);

    auto g_0_yyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 927);

    auto g_0_yyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 928);

    auto g_0_yyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 929);

    auto g_0_yyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 930);

    auto g_0_yyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 931);

    auto g_0_yyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 932);

    auto g_0_yyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 933);

    auto g_0_yyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 934);

    auto g_0_yyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 935);

    auto g_0_yyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 936);

    auto g_0_yyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 937);

    auto g_0_yyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 938);

    auto g_0_yyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 939);

    auto g_0_yyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 940);

    auto g_0_yyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 941);

    auto g_0_yyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 942);

    auto g_0_yyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 943);

    auto g_0_yyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 944);

    auto g_0_yyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 945);

    auto g_0_yyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 946);

    auto g_0_yyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 947);

    auto g_0_yyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 948);

    auto g_0_yyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 949);

    auto g_0_yyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 950);

    auto g_0_yyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 951);

    auto g_0_yzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 952);

    auto g_0_yzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 953);

    auto g_0_yzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 954);

    auto g_0_yzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 955);

    auto g_0_yzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 956);

    auto g_0_yzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 957);

    auto g_0_yzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 958);

    auto g_0_yzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 959);

    auto g_0_yzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 960);

    auto g_0_yzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 961);

    auto g_0_yzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 962);

    auto g_0_yzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 963);

    auto g_0_yzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 964);

    auto g_0_yzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 965);

    auto g_0_yzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 966);

    auto g_0_yzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 967);

    auto g_0_yzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 968);

    auto g_0_yzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 969);

    auto g_0_yzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 970);

    auto g_0_yzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 971);

    auto g_0_yzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 972);

    auto g_0_yzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 973);

    auto g_0_yzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 974);

    auto g_0_yzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 975);

    auto g_0_yzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 976);

    auto g_0_yzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 977);

    auto g_0_yzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 978);

    auto g_0_yzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 979);

    auto g_0_zzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sksi + 980);

    auto g_0_zzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_sksi + 981);

    auto g_0_zzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sksi + 982);

    auto g_0_zzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sksi + 983);

    auto g_0_zzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_sksi + 984);

    auto g_0_zzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sksi + 985);

    auto g_0_zzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sksi + 986);

    auto g_0_zzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_sksi + 987);

    auto g_0_zzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_sksi + 988);

    auto g_0_zzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sksi + 989);

    auto g_0_zzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sksi + 990);

    auto g_0_zzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_sksi + 991);

    auto g_0_zzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sksi + 992);

    auto g_0_zzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_sksi + 993);

    auto g_0_zzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sksi + 994);

    auto g_0_zzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 995);

    auto g_0_zzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 996);

    auto g_0_zzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 997);

    auto g_0_zzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 998);

    auto g_0_zzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 999);

    auto g_0_zzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1000);

    auto g_0_zzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sksi + 1001);

    auto g_0_zzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sksi + 1002);

    auto g_0_zzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sksi + 1003);

    auto g_0_zzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sksi + 1004);

    auto g_0_zzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1005);

    auto g_0_zzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1006);

    auto g_0_zzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sksi + 1007);

    /// Set up components of auxilary buffer : SKSI

    auto g_0_xxxxxxx_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi);

    auto g_0_xxxxxxx_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 1);

    auto g_0_xxxxxxx_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 2);

    auto g_0_xxxxxxx_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 3);

    auto g_0_xxxxxxx_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 4);

    auto g_0_xxxxxxx_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 5);

    auto g_0_xxxxxxx_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 6);

    auto g_0_xxxxxxx_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 7);

    auto g_0_xxxxxxx_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 8);

    auto g_0_xxxxxxx_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 9);

    auto g_0_xxxxxxx_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 10);

    auto g_0_xxxxxxx_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 11);

    auto g_0_xxxxxxx_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 12);

    auto g_0_xxxxxxx_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 13);

    auto g_0_xxxxxxx_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 14);

    auto g_0_xxxxxxx_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 15);

    auto g_0_xxxxxxx_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 16);

    auto g_0_xxxxxxx_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 17);

    auto g_0_xxxxxxx_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 18);

    auto g_0_xxxxxxx_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 19);

    auto g_0_xxxxxxx_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 20);

    auto g_0_xxxxxxx_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 21);

    auto g_0_xxxxxxx_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 22);

    auto g_0_xxxxxxx_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 23);

    auto g_0_xxxxxxx_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 24);

    auto g_0_xxxxxxx_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 25);

    auto g_0_xxxxxxx_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 26);

    auto g_0_xxxxxxx_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 27);

    auto g_0_xxxxxxy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 28);

    auto g_0_xxxxxxy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 29);

    auto g_0_xxxxxxy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 30);

    auto g_0_xxxxxxy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 31);

    auto g_0_xxxxxxy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 33);

    auto g_0_xxxxxxy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 34);

    auto g_0_xxxxxxy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 37);

    auto g_0_xxxxxxy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 38);

    auto g_0_xxxxxxy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 42);

    auto g_0_xxxxxxy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 43);

    auto g_0_xxxxxxy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 48);

    auto g_0_xxxxxxy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 49);

    auto g_0_xxxxxxz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 56);

    auto g_0_xxxxxxz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 57);

    auto g_0_xxxxxxz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 58);

    auto g_0_xxxxxxz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 59);

    auto g_0_xxxxxxz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 60);

    auto g_0_xxxxxxz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 61);

    auto g_0_xxxxxxz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 62);

    auto g_0_xxxxxxz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 63);

    auto g_0_xxxxxxz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 64);

    auto g_0_xxxxxxz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 65);

    auto g_0_xxxxxxz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 66);

    auto g_0_xxxxxxz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 67);

    auto g_0_xxxxxxz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 68);

    auto g_0_xxxxxxz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 69);

    auto g_0_xxxxxxz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 70);

    auto g_0_xxxxxxz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 71);

    auto g_0_xxxxxxz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 72);

    auto g_0_xxxxxxz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 73);

    auto g_0_xxxxxxz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 74);

    auto g_0_xxxxxxz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 75);

    auto g_0_xxxxxxz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 76);

    auto g_0_xxxxxxz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 78);

    auto g_0_xxxxxxz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 79);

    auto g_0_xxxxxxz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 80);

    auto g_0_xxxxxxz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 81);

    auto g_0_xxxxxxz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 82);

    auto g_0_xxxxxxz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 83);

    auto g_0_xxxxxyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 84);

    auto g_0_xxxxxyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 85);

    auto g_0_xxxxxyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 86);

    auto g_0_xxxxxyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 87);

    auto g_0_xxxxxyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 88);

    auto g_0_xxxxxyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 89);

    auto g_0_xxxxxyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 90);

    auto g_0_xxxxxyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 91);

    auto g_0_xxxxxyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 92);

    auto g_0_xxxxxyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 93);

    auto g_0_xxxxxyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 94);

    auto g_0_xxxxxyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 95);

    auto g_0_xxxxxyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 96);

    auto g_0_xxxxxyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 97);

    auto g_0_xxxxxyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 98);

    auto g_0_xxxxxyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 99);

    auto g_0_xxxxxyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 100);

    auto g_0_xxxxxyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 101);

    auto g_0_xxxxxyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 102);

    auto g_0_xxxxxyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 103);

    auto g_0_xxxxxyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 104);

    auto g_0_xxxxxyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 105);

    auto g_0_xxxxxyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 106);

    auto g_0_xxxxxyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 107);

    auto g_0_xxxxxyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 108);

    auto g_0_xxxxxyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 109);

    auto g_0_xxxxxyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 110);

    auto g_0_xxxxxyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 111);

    auto g_0_xxxxxzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 140);

    auto g_0_xxxxxzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 141);

    auto g_0_xxxxxzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 142);

    auto g_0_xxxxxzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 143);

    auto g_0_xxxxxzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 144);

    auto g_0_xxxxxzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 145);

    auto g_0_xxxxxzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 146);

    auto g_0_xxxxxzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 147);

    auto g_0_xxxxxzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 148);

    auto g_0_xxxxxzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 149);

    auto g_0_xxxxxzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 150);

    auto g_0_xxxxxzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 151);

    auto g_0_xxxxxzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 152);

    auto g_0_xxxxxzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 153);

    auto g_0_xxxxxzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 154);

    auto g_0_xxxxxzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 155);

    auto g_0_xxxxxzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 156);

    auto g_0_xxxxxzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 157);

    auto g_0_xxxxxzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 158);

    auto g_0_xxxxxzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 159);

    auto g_0_xxxxxzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 160);

    auto g_0_xxxxxzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 161);

    auto g_0_xxxxxzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 162);

    auto g_0_xxxxxzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 163);

    auto g_0_xxxxxzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 164);

    auto g_0_xxxxxzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 165);

    auto g_0_xxxxxzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 166);

    auto g_0_xxxxxzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 167);

    auto g_0_xxxxyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 168);

    auto g_0_xxxxyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 169);

    auto g_0_xxxxyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 170);

    auto g_0_xxxxyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 171);

    auto g_0_xxxxyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 172);

    auto g_0_xxxxyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 173);

    auto g_0_xxxxyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 174);

    auto g_0_xxxxyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 175);

    auto g_0_xxxxyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 176);

    auto g_0_xxxxyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 177);

    auto g_0_xxxxyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 178);

    auto g_0_xxxxyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 179);

    auto g_0_xxxxyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 180);

    auto g_0_xxxxyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 181);

    auto g_0_xxxxyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 182);

    auto g_0_xxxxyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 183);

    auto g_0_xxxxyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 184);

    auto g_0_xxxxyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 185);

    auto g_0_xxxxyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 186);

    auto g_0_xxxxyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 187);

    auto g_0_xxxxyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 188);

    auto g_0_xxxxyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 189);

    auto g_0_xxxxyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 190);

    auto g_0_xxxxyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 191);

    auto g_0_xxxxyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 192);

    auto g_0_xxxxyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 193);

    auto g_0_xxxxyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 194);

    auto g_0_xxxxyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 195);

    auto g_0_xxxxyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 197);

    auto g_0_xxxxyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 199);

    auto g_0_xxxxyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 202);

    auto g_0_xxxxyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 206);

    auto g_0_xxxxyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 211);

    auto g_0_xxxxyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 224);

    auto g_0_xxxxyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 226);

    auto g_0_xxxxyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 229);

    auto g_0_xxxxyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 233);

    auto g_0_xxxxyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 238);

    auto g_0_xxxxyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 244);

    auto g_0_xxxxzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 252);

    auto g_0_xxxxzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 253);

    auto g_0_xxxxzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 254);

    auto g_0_xxxxzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 255);

    auto g_0_xxxxzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 256);

    auto g_0_xxxxzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 257);

    auto g_0_xxxxzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 258);

    auto g_0_xxxxzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 259);

    auto g_0_xxxxzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 260);

    auto g_0_xxxxzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 261);

    auto g_0_xxxxzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 262);

    auto g_0_xxxxzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 263);

    auto g_0_xxxxzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 264);

    auto g_0_xxxxzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 265);

    auto g_0_xxxxzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 266);

    auto g_0_xxxxzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 267);

    auto g_0_xxxxzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 268);

    auto g_0_xxxxzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 269);

    auto g_0_xxxxzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 270);

    auto g_0_xxxxzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 271);

    auto g_0_xxxxzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 272);

    auto g_0_xxxxzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 273);

    auto g_0_xxxxzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 274);

    auto g_0_xxxxzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 275);

    auto g_0_xxxxzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 276);

    auto g_0_xxxxzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 277);

    auto g_0_xxxxzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 278);

    auto g_0_xxxxzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 279);

    auto g_0_xxxyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 280);

    auto g_0_xxxyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 281);

    auto g_0_xxxyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 282);

    auto g_0_xxxyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 283);

    auto g_0_xxxyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 284);

    auto g_0_xxxyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 285);

    auto g_0_xxxyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 286);

    auto g_0_xxxyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 287);

    auto g_0_xxxyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 288);

    auto g_0_xxxyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 289);

    auto g_0_xxxyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 290);

    auto g_0_xxxyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 291);

    auto g_0_xxxyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 292);

    auto g_0_xxxyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 293);

    auto g_0_xxxyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 294);

    auto g_0_xxxyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 295);

    auto g_0_xxxyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 296);

    auto g_0_xxxyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 297);

    auto g_0_xxxyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 298);

    auto g_0_xxxyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 299);

    auto g_0_xxxyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 300);

    auto g_0_xxxyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 301);

    auto g_0_xxxyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 302);

    auto g_0_xxxyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 303);

    auto g_0_xxxyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 304);

    auto g_0_xxxyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 305);

    auto g_0_xxxyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 306);

    auto g_0_xxxyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 307);

    auto g_0_xxxyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 309);

    auto g_0_xxxyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 311);

    auto g_0_xxxyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 314);

    auto g_0_xxxyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 318);

    auto g_0_xxxyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 323);

    auto g_0_xxxyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 336);

    auto g_0_xxxyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 337);

    auto g_0_xxxyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 338);

    auto g_0_xxxyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 339);

    auto g_0_xxxyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 340);

    auto g_0_xxxyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 341);

    auto g_0_xxxyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 342);

    auto g_0_xxxyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 343);

    auto g_0_xxxyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 344);

    auto g_0_xxxyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 345);

    auto g_0_xxxyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 346);

    auto g_0_xxxyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 347);

    auto g_0_xxxyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 348);

    auto g_0_xxxyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 349);

    auto g_0_xxxyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 350);

    auto g_0_xxxyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 351);

    auto g_0_xxxyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 352);

    auto g_0_xxxyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 353);

    auto g_0_xxxyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 354);

    auto g_0_xxxyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 355);

    auto g_0_xxxyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 356);

    auto g_0_xxxyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 357);

    auto g_0_xxxyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 358);

    auto g_0_xxxyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 359);

    auto g_0_xxxyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 360);

    auto g_0_xxxyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 361);

    auto g_0_xxxyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 362);

    auto g_0_xxxyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 363);

    auto g_0_xxxyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 364);

    auto g_0_xxxyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 366);

    auto g_0_xxxyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 369);

    auto g_0_xxxyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 373);

    auto g_0_xxxyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 378);

    auto g_0_xxxyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 384);

    auto g_0_xxxzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 392);

    auto g_0_xxxzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 393);

    auto g_0_xxxzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 394);

    auto g_0_xxxzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 395);

    auto g_0_xxxzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 396);

    auto g_0_xxxzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 397);

    auto g_0_xxxzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 398);

    auto g_0_xxxzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 399);

    auto g_0_xxxzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 400);

    auto g_0_xxxzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 401);

    auto g_0_xxxzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 402);

    auto g_0_xxxzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 403);

    auto g_0_xxxzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 404);

    auto g_0_xxxzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 405);

    auto g_0_xxxzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 406);

    auto g_0_xxxzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 407);

    auto g_0_xxxzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 408);

    auto g_0_xxxzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 409);

    auto g_0_xxxzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 410);

    auto g_0_xxxzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 411);

    auto g_0_xxxzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 412);

    auto g_0_xxxzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 413);

    auto g_0_xxxzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 414);

    auto g_0_xxxzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 415);

    auto g_0_xxxzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 416);

    auto g_0_xxxzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 417);

    auto g_0_xxxzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 418);

    auto g_0_xxxzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 419);

    auto g_0_xxyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 420);

    auto g_0_xxyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 421);

    auto g_0_xxyyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 422);

    auto g_0_xxyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 423);

    auto g_0_xxyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 424);

    auto g_0_xxyyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 425);

    auto g_0_xxyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 426);

    auto g_0_xxyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 427);

    auto g_0_xxyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 428);

    auto g_0_xxyyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 429);

    auto g_0_xxyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 430);

    auto g_0_xxyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 431);

    auto g_0_xxyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 432);

    auto g_0_xxyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 433);

    auto g_0_xxyyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 434);

    auto g_0_xxyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 435);

    auto g_0_xxyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 436);

    auto g_0_xxyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 437);

    auto g_0_xxyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 438);

    auto g_0_xxyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 439);

    auto g_0_xxyyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 440);

    auto g_0_xxyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 441);

    auto g_0_xxyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 442);

    auto g_0_xxyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 443);

    auto g_0_xxyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 444);

    auto g_0_xxyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 445);

    auto g_0_xxyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 446);

    auto g_0_xxyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 447);

    auto g_0_xxyyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 449);

    auto g_0_xxyyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 451);

    auto g_0_xxyyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 454);

    auto g_0_xxyyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 458);

    auto g_0_xxyyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 463);

    auto g_0_xxyyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 476);

    auto g_0_xxyyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 477);

    auto g_0_xxyyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 478);

    auto g_0_xxyyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 479);

    auto g_0_xxyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 480);

    auto g_0_xxyyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 481);

    auto g_0_xxyyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 482);

    auto g_0_xxyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 483);

    auto g_0_xxyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 484);

    auto g_0_xxyyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 485);

    auto g_0_xxyyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 486);

    auto g_0_xxyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 487);

    auto g_0_xxyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 488);

    auto g_0_xxyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 489);

    auto g_0_xxyyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 490);

    auto g_0_xxyyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 491);

    auto g_0_xxyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 492);

    auto g_0_xxyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 493);

    auto g_0_xxyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 494);

    auto g_0_xxyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 495);

    auto g_0_xxyyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 496);

    auto g_0_xxyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 497);

    auto g_0_xxyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 498);

    auto g_0_xxyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 499);

    auto g_0_xxyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 500);

    auto g_0_xxyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 501);

    auto g_0_xxyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 502);

    auto g_0_xxyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 503);

    auto g_0_xxyyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 504);

    auto g_0_xxyyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 505);

    auto g_0_xxyyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 506);

    auto g_0_xxyyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 507);

    auto g_0_xxyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 508);

    auto g_0_xxyyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 509);

    auto g_0_xxyyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 510);

    auto g_0_xxyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 511);

    auto g_0_xxyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 512);

    auto g_0_xxyyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 513);

    auto g_0_xxyyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 514);

    auto g_0_xxyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 515);

    auto g_0_xxyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 516);

    auto g_0_xxyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 517);

    auto g_0_xxyyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 518);

    auto g_0_xxyyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 519);

    auto g_0_xxyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 520);

    auto g_0_xxyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 521);

    auto g_0_xxyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 522);

    auto g_0_xxyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 523);

    auto g_0_xxyyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 524);

    auto g_0_xxyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 525);

    auto g_0_xxyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 526);

    auto g_0_xxyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 527);

    auto g_0_xxyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 528);

    auto g_0_xxyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 529);

    auto g_0_xxyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 530);

    auto g_0_xxyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 531);

    auto g_0_xxyzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 532);

    auto g_0_xxyzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 534);

    auto g_0_xxyzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 537);

    auto g_0_xxyzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 541);

    auto g_0_xxyzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 546);

    auto g_0_xxyzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 552);

    auto g_0_xxzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 560);

    auto g_0_xxzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 561);

    auto g_0_xxzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 562);

    auto g_0_xxzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 563);

    auto g_0_xxzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 564);

    auto g_0_xxzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 565);

    auto g_0_xxzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 566);

    auto g_0_xxzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 567);

    auto g_0_xxzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 568);

    auto g_0_xxzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 569);

    auto g_0_xxzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 570);

    auto g_0_xxzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 571);

    auto g_0_xxzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 572);

    auto g_0_xxzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 573);

    auto g_0_xxzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 574);

    auto g_0_xxzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 575);

    auto g_0_xxzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 576);

    auto g_0_xxzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 577);

    auto g_0_xxzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 578);

    auto g_0_xxzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 579);

    auto g_0_xxzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 580);

    auto g_0_xxzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 581);

    auto g_0_xxzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 582);

    auto g_0_xxzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 583);

    auto g_0_xxzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 584);

    auto g_0_xxzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 585);

    auto g_0_xxzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 586);

    auto g_0_xxzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 587);

    auto g_0_xyyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 588);

    auto g_0_xyyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 589);

    auto g_0_xyyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 591);

    auto g_0_xyyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 592);

    auto g_0_xyyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 594);

    auto g_0_xyyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 595);

    auto g_0_xyyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 596);

    auto g_0_xyyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 598);

    auto g_0_xyyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 599);

    auto g_0_xyyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 600);

    auto g_0_xyyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 601);

    auto g_0_xyyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 603);

    auto g_0_xyyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 604);

    auto g_0_xyyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 605);

    auto g_0_xyyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 606);

    auto g_0_xyyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 607);

    auto g_0_xyyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 609);

    auto g_0_xyyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 610);

    auto g_0_xyyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 611);

    auto g_0_xyyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 612);

    auto g_0_xyyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 613);

    auto g_0_xyyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 614);

    auto g_0_xyyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 615);

    auto g_0_xyyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 648);

    auto g_0_xyyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 651);

    auto g_0_xyyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 652);

    auto g_0_xyyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 655);

    auto g_0_xyyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 656);

    auto g_0_xyyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 657);

    auto g_0_xyyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 660);

    auto g_0_xyyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 661);

    auto g_0_xyyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 662);

    auto g_0_xyyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 663);

    auto g_0_xyyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 665);

    auto g_0_xyyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 666);

    auto g_0_xyyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 667);

    auto g_0_xyyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 668);

    auto g_0_xyyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 669);

    auto g_0_xyyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 670);

    auto g_0_xyyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 671);

    auto g_0_xyyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 676);

    auto g_0_xyyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 679);

    auto g_0_xyyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 680);

    auto g_0_xyyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 683);

    auto g_0_xyyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 684);

    auto g_0_xyyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 685);

    auto g_0_xyyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 688);

    auto g_0_xyyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 689);

    auto g_0_xyyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 690);

    auto g_0_xyyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 691);

    auto g_0_xyyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 693);

    auto g_0_xyyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 694);

    auto g_0_xyyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 695);

    auto g_0_xyyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 696);

    auto g_0_xyyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 697);

    auto g_0_xyyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 698);

    auto g_0_xyyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 699);

    auto g_0_xyyzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 704);

    auto g_0_xyyzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 707);

    auto g_0_xyyzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 708);

    auto g_0_xyyzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 711);

    auto g_0_xyyzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 712);

    auto g_0_xyyzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 713);

    auto g_0_xyyzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 716);

    auto g_0_xyyzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 717);

    auto g_0_xyyzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 718);

    auto g_0_xyyzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 719);

    auto g_0_xyyzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 721);

    auto g_0_xyyzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 722);

    auto g_0_xyyzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 723);

    auto g_0_xyyzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 724);

    auto g_0_xyyzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 725);

    auto g_0_xyyzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 726);

    auto g_0_xyyzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 727);

    auto g_0_xzzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 756);

    auto g_0_xzzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 758);

    auto g_0_xzzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 760);

    auto g_0_xzzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 761);

    auto g_0_xzzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 763);

    auto g_0_xzzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 764);

    auto g_0_xzzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 765);

    auto g_0_xzzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 767);

    auto g_0_xzzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 768);

    auto g_0_xzzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 769);

    auto g_0_xzzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 770);

    auto g_0_xzzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 772);

    auto g_0_xzzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 773);

    auto g_0_xzzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 774);

    auto g_0_xzzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 775);

    auto g_0_xzzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 776);

    auto g_0_xzzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 777);

    auto g_0_xzzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 778);

    auto g_0_xzzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 779);

    auto g_0_xzzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 780);

    auto g_0_xzzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 781);

    auto g_0_xzzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 782);

    auto g_0_xzzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 783);

    auto g_0_yyyyyyy_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 784);

    auto g_0_yyyyyyy_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 785);

    auto g_0_yyyyyyy_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 786);

    auto g_0_yyyyyyy_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 787);

    auto g_0_yyyyyyy_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 788);

    auto g_0_yyyyyyy_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 789);

    auto g_0_yyyyyyy_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 790);

    auto g_0_yyyyyyy_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 791);

    auto g_0_yyyyyyy_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 792);

    auto g_0_yyyyyyy_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 793);

    auto g_0_yyyyyyy_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 794);

    auto g_0_yyyyyyy_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 795);

    auto g_0_yyyyyyy_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 796);

    auto g_0_yyyyyyy_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 797);

    auto g_0_yyyyyyy_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 798);

    auto g_0_yyyyyyy_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 799);

    auto g_0_yyyyyyy_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 800);

    auto g_0_yyyyyyy_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 801);

    auto g_0_yyyyyyy_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 802);

    auto g_0_yyyyyyy_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 803);

    auto g_0_yyyyyyy_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 804);

    auto g_0_yyyyyyy_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 805);

    auto g_0_yyyyyyy_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 806);

    auto g_0_yyyyyyy_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 807);

    auto g_0_yyyyyyy_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 808);

    auto g_0_yyyyyyy_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 809);

    auto g_0_yyyyyyy_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 810);

    auto g_0_yyyyyyy_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 811);

    auto g_0_yyyyyyz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 813);

    auto g_0_yyyyyyz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 814);

    auto g_0_yyyyyyz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 815);

    auto g_0_yyyyyyz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 816);

    auto g_0_yyyyyyz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 817);

    auto g_0_yyyyyyz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 818);

    auto g_0_yyyyyyz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 819);

    auto g_0_yyyyyyz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 820);

    auto g_0_yyyyyyz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 821);

    auto g_0_yyyyyyz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 822);

    auto g_0_yyyyyyz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 823);

    auto g_0_yyyyyyz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 824);

    auto g_0_yyyyyyz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 825);

    auto g_0_yyyyyyz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 826);

    auto g_0_yyyyyyz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 827);

    auto g_0_yyyyyyz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 828);

    auto g_0_yyyyyyz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 829);

    auto g_0_yyyyyyz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 830);

    auto g_0_yyyyyyz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 831);

    auto g_0_yyyyyyz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 832);

    auto g_0_yyyyyyz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 833);

    auto g_0_yyyyyyz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 834);

    auto g_0_yyyyyyz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 835);

    auto g_0_yyyyyyz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 836);

    auto g_0_yyyyyyz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 837);

    auto g_0_yyyyyyz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 838);

    auto g_0_yyyyyyz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 839);

    auto g_0_yyyyyzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 840);

    auto g_0_yyyyyzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 841);

    auto g_0_yyyyyzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 842);

    auto g_0_yyyyyzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 843);

    auto g_0_yyyyyzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 844);

    auto g_0_yyyyyzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 845);

    auto g_0_yyyyyzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 846);

    auto g_0_yyyyyzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 847);

    auto g_0_yyyyyzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 848);

    auto g_0_yyyyyzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 849);

    auto g_0_yyyyyzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 850);

    auto g_0_yyyyyzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 851);

    auto g_0_yyyyyzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 852);

    auto g_0_yyyyyzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 853);

    auto g_0_yyyyyzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 854);

    auto g_0_yyyyyzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 855);

    auto g_0_yyyyyzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 856);

    auto g_0_yyyyyzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 857);

    auto g_0_yyyyyzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 858);

    auto g_0_yyyyyzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 859);

    auto g_0_yyyyyzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 860);

    auto g_0_yyyyyzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 861);

    auto g_0_yyyyyzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 862);

    auto g_0_yyyyyzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 863);

    auto g_0_yyyyyzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 864);

    auto g_0_yyyyyzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 865);

    auto g_0_yyyyyzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 866);

    auto g_0_yyyyyzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 867);

    auto g_0_yyyyzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 868);

    auto g_0_yyyyzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 869);

    auto g_0_yyyyzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 870);

    auto g_0_yyyyzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 871);

    auto g_0_yyyyzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 872);

    auto g_0_yyyyzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 873);

    auto g_0_yyyyzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 874);

    auto g_0_yyyyzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 875);

    auto g_0_yyyyzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 876);

    auto g_0_yyyyzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 877);

    auto g_0_yyyyzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 878);

    auto g_0_yyyyzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 879);

    auto g_0_yyyyzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 880);

    auto g_0_yyyyzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 881);

    auto g_0_yyyyzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 882);

    auto g_0_yyyyzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 883);

    auto g_0_yyyyzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 884);

    auto g_0_yyyyzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 885);

    auto g_0_yyyyzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 886);

    auto g_0_yyyyzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 887);

    auto g_0_yyyyzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 888);

    auto g_0_yyyyzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 889);

    auto g_0_yyyyzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 890);

    auto g_0_yyyyzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 891);

    auto g_0_yyyyzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 892);

    auto g_0_yyyyzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 893);

    auto g_0_yyyyzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 894);

    auto g_0_yyyyzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 895);

    auto g_0_yyyzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 896);

    auto g_0_yyyzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 897);

    auto g_0_yyyzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 898);

    auto g_0_yyyzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 899);

    auto g_0_yyyzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 900);

    auto g_0_yyyzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 901);

    auto g_0_yyyzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 902);

    auto g_0_yyyzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 903);

    auto g_0_yyyzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 904);

    auto g_0_yyyzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 905);

    auto g_0_yyyzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 906);

    auto g_0_yyyzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 907);

    auto g_0_yyyzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 908);

    auto g_0_yyyzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 909);

    auto g_0_yyyzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 910);

    auto g_0_yyyzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 911);

    auto g_0_yyyzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 912);

    auto g_0_yyyzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 913);

    auto g_0_yyyzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 914);

    auto g_0_yyyzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 915);

    auto g_0_yyyzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 916);

    auto g_0_yyyzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 917);

    auto g_0_yyyzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 918);

    auto g_0_yyyzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 919);

    auto g_0_yyyzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 920);

    auto g_0_yyyzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 921);

    auto g_0_yyyzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 922);

    auto g_0_yyyzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 923);

    auto g_0_yyzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 924);

    auto g_0_yyzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 925);

    auto g_0_yyzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 926);

    auto g_0_yyzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 927);

    auto g_0_yyzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 928);

    auto g_0_yyzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 929);

    auto g_0_yyzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 930);

    auto g_0_yyzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 931);

    auto g_0_yyzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 932);

    auto g_0_yyzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 933);

    auto g_0_yyzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 934);

    auto g_0_yyzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 935);

    auto g_0_yyzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 936);

    auto g_0_yyzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 937);

    auto g_0_yyzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 938);

    auto g_0_yyzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 939);

    auto g_0_yyzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 940);

    auto g_0_yyzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 941);

    auto g_0_yyzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 942);

    auto g_0_yyzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 943);

    auto g_0_yyzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 944);

    auto g_0_yyzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 945);

    auto g_0_yyzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 946);

    auto g_0_yyzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 947);

    auto g_0_yyzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 948);

    auto g_0_yyzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 949);

    auto g_0_yyzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 950);

    auto g_0_yyzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 951);

    auto g_0_yzzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 952);

    auto g_0_yzzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 953);

    auto g_0_yzzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 954);

    auto g_0_yzzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 955);

    auto g_0_yzzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 956);

    auto g_0_yzzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 957);

    auto g_0_yzzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 958);

    auto g_0_yzzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 959);

    auto g_0_yzzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 960);

    auto g_0_yzzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 961);

    auto g_0_yzzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 962);

    auto g_0_yzzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 963);

    auto g_0_yzzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 964);

    auto g_0_yzzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 965);

    auto g_0_yzzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 966);

    auto g_0_yzzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 967);

    auto g_0_yzzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 968);

    auto g_0_yzzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 969);

    auto g_0_yzzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 970);

    auto g_0_yzzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 971);

    auto g_0_yzzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 972);

    auto g_0_yzzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 973);

    auto g_0_yzzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 974);

    auto g_0_yzzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 975);

    auto g_0_yzzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 976);

    auto g_0_yzzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 977);

    auto g_0_yzzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 978);

    auto g_0_yzzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 979);

    auto g_0_zzzzzzz_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sksi + 980);

    auto g_0_zzzzzzz_0_xxxxxy_1 = pbuffer.data(idx_eri_1_sksi + 981);

    auto g_0_zzzzzzz_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sksi + 982);

    auto g_0_zzzzzzz_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sksi + 983);

    auto g_0_zzzzzzz_0_xxxxyz_1 = pbuffer.data(idx_eri_1_sksi + 984);

    auto g_0_zzzzzzz_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sksi + 985);

    auto g_0_zzzzzzz_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sksi + 986);

    auto g_0_zzzzzzz_0_xxxyyz_1 = pbuffer.data(idx_eri_1_sksi + 987);

    auto g_0_zzzzzzz_0_xxxyzz_1 = pbuffer.data(idx_eri_1_sksi + 988);

    auto g_0_zzzzzzz_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sksi + 989);

    auto g_0_zzzzzzz_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sksi + 990);

    auto g_0_zzzzzzz_0_xxyyyz_1 = pbuffer.data(idx_eri_1_sksi + 991);

    auto g_0_zzzzzzz_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sksi + 992);

    auto g_0_zzzzzzz_0_xxyzzz_1 = pbuffer.data(idx_eri_1_sksi + 993);

    auto g_0_zzzzzzz_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sksi + 994);

    auto g_0_zzzzzzz_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 995);

    auto g_0_zzzzzzz_0_xyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 996);

    auto g_0_zzzzzzz_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 997);

    auto g_0_zzzzzzz_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 998);

    auto g_0_zzzzzzz_0_xyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 999);

    auto g_0_zzzzzzz_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 1000);

    auto g_0_zzzzzzz_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sksi + 1001);

    auto g_0_zzzzzzz_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sksi + 1002);

    auto g_0_zzzzzzz_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sksi + 1003);

    auto g_0_zzzzzzz_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sksi + 1004);

    auto g_0_zzzzzzz_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sksi + 1005);

    auto g_0_zzzzzzz_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 1006);

    auto g_0_zzzzzzz_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sksi + 1007);

    /// Set up 0-28 components of targeted buffer : SLSI

    auto g_0_xxxxxxxx_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi);

    auto g_0_xxxxxxxx_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1);

    auto g_0_xxxxxxxx_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 2);

    auto g_0_xxxxxxxx_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 3);

    auto g_0_xxxxxxxx_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 4);

    auto g_0_xxxxxxxx_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 5);

    auto g_0_xxxxxxxx_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 6);

    auto g_0_xxxxxxxx_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 7);

    auto g_0_xxxxxxxx_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 8);

    auto g_0_xxxxxxxx_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 9);

    auto g_0_xxxxxxxx_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 10);

    auto g_0_xxxxxxxx_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 11);

    auto g_0_xxxxxxxx_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 12);

    auto g_0_xxxxxxxx_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 13);

    auto g_0_xxxxxxxx_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 14);

    auto g_0_xxxxxxxx_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 15);

    auto g_0_xxxxxxxx_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 16);

    auto g_0_xxxxxxxx_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 17);

    auto g_0_xxxxxxxx_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 18);

    auto g_0_xxxxxxxx_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 19);

    auto g_0_xxxxxxxx_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 20);

    auto g_0_xxxxxxxx_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 21);

    auto g_0_xxxxxxxx_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 22);

    auto g_0_xxxxxxxx_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 23);

    auto g_0_xxxxxxxx_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 24);

    auto g_0_xxxxxxxx_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 25);

    auto g_0_xxxxxxxx_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 26);

    auto g_0_xxxxxxxx_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 27);

#pragma omp simd aligned(g_0_xxxxxx_0_xxxxxx_0,       \
                             g_0_xxxxxx_0_xxxxxx_1,   \
                             g_0_xxxxxx_0_xxxxxy_0,   \
                             g_0_xxxxxx_0_xxxxxy_1,   \
                             g_0_xxxxxx_0_xxxxxz_0,   \
                             g_0_xxxxxx_0_xxxxxz_1,   \
                             g_0_xxxxxx_0_xxxxyy_0,   \
                             g_0_xxxxxx_0_xxxxyy_1,   \
                             g_0_xxxxxx_0_xxxxyz_0,   \
                             g_0_xxxxxx_0_xxxxyz_1,   \
                             g_0_xxxxxx_0_xxxxzz_0,   \
                             g_0_xxxxxx_0_xxxxzz_1,   \
                             g_0_xxxxxx_0_xxxyyy_0,   \
                             g_0_xxxxxx_0_xxxyyy_1,   \
                             g_0_xxxxxx_0_xxxyyz_0,   \
                             g_0_xxxxxx_0_xxxyyz_1,   \
                             g_0_xxxxxx_0_xxxyzz_0,   \
                             g_0_xxxxxx_0_xxxyzz_1,   \
                             g_0_xxxxxx_0_xxxzzz_0,   \
                             g_0_xxxxxx_0_xxxzzz_1,   \
                             g_0_xxxxxx_0_xxyyyy_0,   \
                             g_0_xxxxxx_0_xxyyyy_1,   \
                             g_0_xxxxxx_0_xxyyyz_0,   \
                             g_0_xxxxxx_0_xxyyyz_1,   \
                             g_0_xxxxxx_0_xxyyzz_0,   \
                             g_0_xxxxxx_0_xxyyzz_1,   \
                             g_0_xxxxxx_0_xxyzzz_0,   \
                             g_0_xxxxxx_0_xxyzzz_1,   \
                             g_0_xxxxxx_0_xxzzzz_0,   \
                             g_0_xxxxxx_0_xxzzzz_1,   \
                             g_0_xxxxxx_0_xyyyyy_0,   \
                             g_0_xxxxxx_0_xyyyyy_1,   \
                             g_0_xxxxxx_0_xyyyyz_0,   \
                             g_0_xxxxxx_0_xyyyyz_1,   \
                             g_0_xxxxxx_0_xyyyzz_0,   \
                             g_0_xxxxxx_0_xyyyzz_1,   \
                             g_0_xxxxxx_0_xyyzzz_0,   \
                             g_0_xxxxxx_0_xyyzzz_1,   \
                             g_0_xxxxxx_0_xyzzzz_0,   \
                             g_0_xxxxxx_0_xyzzzz_1,   \
                             g_0_xxxxxx_0_xzzzzz_0,   \
                             g_0_xxxxxx_0_xzzzzz_1,   \
                             g_0_xxxxxx_0_yyyyyy_0,   \
                             g_0_xxxxxx_0_yyyyyy_1,   \
                             g_0_xxxxxx_0_yyyyyz_0,   \
                             g_0_xxxxxx_0_yyyyyz_1,   \
                             g_0_xxxxxx_0_yyyyzz_0,   \
                             g_0_xxxxxx_0_yyyyzz_1,   \
                             g_0_xxxxxx_0_yyyzzz_0,   \
                             g_0_xxxxxx_0_yyyzzz_1,   \
                             g_0_xxxxxx_0_yyzzzz_0,   \
                             g_0_xxxxxx_0_yyzzzz_1,   \
                             g_0_xxxxxx_0_yzzzzz_0,   \
                             g_0_xxxxxx_0_yzzzzz_1,   \
                             g_0_xxxxxx_0_zzzzzz_0,   \
                             g_0_xxxxxx_0_zzzzzz_1,   \
                             g_0_xxxxxxx_0_xxxxx_1,   \
                             g_0_xxxxxxx_0_xxxxxx_0,  \
                             g_0_xxxxxxx_0_xxxxxx_1,  \
                             g_0_xxxxxxx_0_xxxxxy_0,  \
                             g_0_xxxxxxx_0_xxxxxy_1,  \
                             g_0_xxxxxxx_0_xxxxxz_0,  \
                             g_0_xxxxxxx_0_xxxxxz_1,  \
                             g_0_xxxxxxx_0_xxxxy_1,   \
                             g_0_xxxxxxx_0_xxxxyy_0,  \
                             g_0_xxxxxxx_0_xxxxyy_1,  \
                             g_0_xxxxxxx_0_xxxxyz_0,  \
                             g_0_xxxxxxx_0_xxxxyz_1,  \
                             g_0_xxxxxxx_0_xxxxz_1,   \
                             g_0_xxxxxxx_0_xxxxzz_0,  \
                             g_0_xxxxxxx_0_xxxxzz_1,  \
                             g_0_xxxxxxx_0_xxxyy_1,   \
                             g_0_xxxxxxx_0_xxxyyy_0,  \
                             g_0_xxxxxxx_0_xxxyyy_1,  \
                             g_0_xxxxxxx_0_xxxyyz_0,  \
                             g_0_xxxxxxx_0_xxxyyz_1,  \
                             g_0_xxxxxxx_0_xxxyz_1,   \
                             g_0_xxxxxxx_0_xxxyzz_0,  \
                             g_0_xxxxxxx_0_xxxyzz_1,  \
                             g_0_xxxxxxx_0_xxxzz_1,   \
                             g_0_xxxxxxx_0_xxxzzz_0,  \
                             g_0_xxxxxxx_0_xxxzzz_1,  \
                             g_0_xxxxxxx_0_xxyyy_1,   \
                             g_0_xxxxxxx_0_xxyyyy_0,  \
                             g_0_xxxxxxx_0_xxyyyy_1,  \
                             g_0_xxxxxxx_0_xxyyyz_0,  \
                             g_0_xxxxxxx_0_xxyyyz_1,  \
                             g_0_xxxxxxx_0_xxyyz_1,   \
                             g_0_xxxxxxx_0_xxyyzz_0,  \
                             g_0_xxxxxxx_0_xxyyzz_1,  \
                             g_0_xxxxxxx_0_xxyzz_1,   \
                             g_0_xxxxxxx_0_xxyzzz_0,  \
                             g_0_xxxxxxx_0_xxyzzz_1,  \
                             g_0_xxxxxxx_0_xxzzz_1,   \
                             g_0_xxxxxxx_0_xxzzzz_0,  \
                             g_0_xxxxxxx_0_xxzzzz_1,  \
                             g_0_xxxxxxx_0_xyyyy_1,   \
                             g_0_xxxxxxx_0_xyyyyy_0,  \
                             g_0_xxxxxxx_0_xyyyyy_1,  \
                             g_0_xxxxxxx_0_xyyyyz_0,  \
                             g_0_xxxxxxx_0_xyyyyz_1,  \
                             g_0_xxxxxxx_0_xyyyz_1,   \
                             g_0_xxxxxxx_0_xyyyzz_0,  \
                             g_0_xxxxxxx_0_xyyyzz_1,  \
                             g_0_xxxxxxx_0_xyyzz_1,   \
                             g_0_xxxxxxx_0_xyyzzz_0,  \
                             g_0_xxxxxxx_0_xyyzzz_1,  \
                             g_0_xxxxxxx_0_xyzzz_1,   \
                             g_0_xxxxxxx_0_xyzzzz_0,  \
                             g_0_xxxxxxx_0_xyzzzz_1,  \
                             g_0_xxxxxxx_0_xzzzz_1,   \
                             g_0_xxxxxxx_0_xzzzzz_0,  \
                             g_0_xxxxxxx_0_xzzzzz_1,  \
                             g_0_xxxxxxx_0_yyyyy_1,   \
                             g_0_xxxxxxx_0_yyyyyy_0,  \
                             g_0_xxxxxxx_0_yyyyyy_1,  \
                             g_0_xxxxxxx_0_yyyyyz_0,  \
                             g_0_xxxxxxx_0_yyyyyz_1,  \
                             g_0_xxxxxxx_0_yyyyz_1,   \
                             g_0_xxxxxxx_0_yyyyzz_0,  \
                             g_0_xxxxxxx_0_yyyyzz_1,  \
                             g_0_xxxxxxx_0_yyyzz_1,   \
                             g_0_xxxxxxx_0_yyyzzz_0,  \
                             g_0_xxxxxxx_0_yyyzzz_1,  \
                             g_0_xxxxxxx_0_yyzzz_1,   \
                             g_0_xxxxxxx_0_yyzzzz_0,  \
                             g_0_xxxxxxx_0_yyzzzz_1,  \
                             g_0_xxxxxxx_0_yzzzz_1,   \
                             g_0_xxxxxxx_0_yzzzzz_0,  \
                             g_0_xxxxxxx_0_yzzzzz_1,  \
                             g_0_xxxxxxx_0_zzzzz_1,   \
                             g_0_xxxxxxx_0_zzzzzz_0,  \
                             g_0_xxxxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxxxx_0_xxxxxx_0, \
                             g_0_xxxxxxxx_0_xxxxxy_0, \
                             g_0_xxxxxxxx_0_xxxxxz_0, \
                             g_0_xxxxxxxx_0_xxxxyy_0, \
                             g_0_xxxxxxxx_0_xxxxyz_0, \
                             g_0_xxxxxxxx_0_xxxxzz_0, \
                             g_0_xxxxxxxx_0_xxxyyy_0, \
                             g_0_xxxxxxxx_0_xxxyyz_0, \
                             g_0_xxxxxxxx_0_xxxyzz_0, \
                             g_0_xxxxxxxx_0_xxxzzz_0, \
                             g_0_xxxxxxxx_0_xxyyyy_0, \
                             g_0_xxxxxxxx_0_xxyyyz_0, \
                             g_0_xxxxxxxx_0_xxyyzz_0, \
                             g_0_xxxxxxxx_0_xxyzzz_0, \
                             g_0_xxxxxxxx_0_xxzzzz_0, \
                             g_0_xxxxxxxx_0_xyyyyy_0, \
                             g_0_xxxxxxxx_0_xyyyyz_0, \
                             g_0_xxxxxxxx_0_xyyyzz_0, \
                             g_0_xxxxxxxx_0_xyyzzz_0, \
                             g_0_xxxxxxxx_0_xyzzzz_0, \
                             g_0_xxxxxxxx_0_xzzzzz_0, \
                             g_0_xxxxxxxx_0_yyyyyy_0, \
                             g_0_xxxxxxxx_0_yyyyyz_0, \
                             g_0_xxxxxxxx_0_yyyyzz_0, \
                             g_0_xxxxxxxx_0_yyyzzz_0, \
                             g_0_xxxxxxxx_0_yyzzzz_0, \
                             g_0_xxxxxxxx_0_yzzzzz_0, \
                             g_0_xxxxxxxx_0_zzzzzz_0, \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_xxxxxx_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxx_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxx_1[i] * fti_ab_0 +
                                     6.0 * g_0_xxxxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxx_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxxxx_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxy_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxy_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxxz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxxz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxxz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxyy_0[i] = 7.0 * g_0_xxxxxx_0_xxxxyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyy_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxyz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxzz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyyy_0[i] = 7.0 * g_0_xxxxxx_0_xxxyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyy_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyyz_0[i] = 7.0 * g_0_xxxxxx_0_xxxyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxzzz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyyy_0[i] = 7.0 * g_0_xxxxxx_0_xxyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyy_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyyz_0[i] = 7.0 * g_0_xxxxxx_0_xxyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyzz_0[i] = 7.0 * g_0_xxxxxx_0_xxyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzzzz_0[i] * pb_x +
                                     g_0_xxxxxxx_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyyy_0[i] = 7.0 * g_0_xxxxxx_0_xyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyy_0[i] * pb_x + g_0_xxxxxxx_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyyz_0[i] = 7.0 * g_0_xxxxxx_0_xyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyz_0[i] * pb_x + g_0_xxxxxxx_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyzz_0[i] = 7.0 * g_0_xxxxxx_0_xyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyzz_0[i] * pb_x + g_0_xxxxxxx_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyzzz_0[i] = 7.0 * g_0_xxxxxx_0_xyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzzz_0[i] * pb_x + g_0_xxxxxxx_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyyy_0[i] = 7.0 * g_0_xxxxxx_0_yyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyyyy_0[i] * pb_x + g_0_xxxxxxx_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyyz_0[i] = 7.0 * g_0_xxxxxx_0_yyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyyyz_0[i] * pb_x + g_0_xxxxxxx_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyzz_0[i] = 7.0 * g_0_xxxxxx_0_yyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyyzz_0[i] * pb_x + g_0_xxxxxxx_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyzzz_0[i] = 7.0 * g_0_xxxxxx_0_yyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyyzzz_0[i] * pb_x + g_0_xxxxxxx_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyzzzz_0[i] = 7.0 * g_0_xxxxxx_0_yyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yyzzzz_0[i] * pb_x + g_0_xxxxxxx_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_yzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_yzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_zzzzzz_0[i] = 7.0 * g_0_xxxxxx_0_zzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxxx_0_zzzzzz_0[i] * pb_x + g_0_xxxxxxx_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : SLSI

    auto g_0_xxxxxxxy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 28);

    auto g_0_xxxxxxxy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 29);

    auto g_0_xxxxxxxy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 30);

    auto g_0_xxxxxxxy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 31);

    auto g_0_xxxxxxxy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 32);

    auto g_0_xxxxxxxy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 33);

    auto g_0_xxxxxxxy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 34);

    auto g_0_xxxxxxxy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 35);

    auto g_0_xxxxxxxy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 36);

    auto g_0_xxxxxxxy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 37);

    auto g_0_xxxxxxxy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 38);

    auto g_0_xxxxxxxy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 39);

    auto g_0_xxxxxxxy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 40);

    auto g_0_xxxxxxxy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 41);

    auto g_0_xxxxxxxy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 42);

    auto g_0_xxxxxxxy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 43);

    auto g_0_xxxxxxxy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 44);

    auto g_0_xxxxxxxy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 45);

    auto g_0_xxxxxxxy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 46);

    auto g_0_xxxxxxxy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 47);

    auto g_0_xxxxxxxy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 48);

    auto g_0_xxxxxxxy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 49);

    auto g_0_xxxxxxxy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 50);

    auto g_0_xxxxxxxy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 51);

    auto g_0_xxxxxxxy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 52);

    auto g_0_xxxxxxxy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 53);

    auto g_0_xxxxxxxy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 54);

    auto g_0_xxxxxxxy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 55);

#pragma omp simd aligned(g_0_xxxxxxx_0_xxxxx_1,       \
                             g_0_xxxxxxx_0_xxxxxx_0,  \
                             g_0_xxxxxxx_0_xxxxxx_1,  \
                             g_0_xxxxxxx_0_xxxxxy_0,  \
                             g_0_xxxxxxx_0_xxxxxy_1,  \
                             g_0_xxxxxxx_0_xxxxxz_0,  \
                             g_0_xxxxxxx_0_xxxxxz_1,  \
                             g_0_xxxxxxx_0_xxxxy_1,   \
                             g_0_xxxxxxx_0_xxxxyy_0,  \
                             g_0_xxxxxxx_0_xxxxyy_1,  \
                             g_0_xxxxxxx_0_xxxxyz_0,  \
                             g_0_xxxxxxx_0_xxxxyz_1,  \
                             g_0_xxxxxxx_0_xxxxz_1,   \
                             g_0_xxxxxxx_0_xxxxzz_0,  \
                             g_0_xxxxxxx_0_xxxxzz_1,  \
                             g_0_xxxxxxx_0_xxxyy_1,   \
                             g_0_xxxxxxx_0_xxxyyy_0,  \
                             g_0_xxxxxxx_0_xxxyyy_1,  \
                             g_0_xxxxxxx_0_xxxyyz_0,  \
                             g_0_xxxxxxx_0_xxxyyz_1,  \
                             g_0_xxxxxxx_0_xxxyz_1,   \
                             g_0_xxxxxxx_0_xxxyzz_0,  \
                             g_0_xxxxxxx_0_xxxyzz_1,  \
                             g_0_xxxxxxx_0_xxxzz_1,   \
                             g_0_xxxxxxx_0_xxxzzz_0,  \
                             g_0_xxxxxxx_0_xxxzzz_1,  \
                             g_0_xxxxxxx_0_xxyyy_1,   \
                             g_0_xxxxxxx_0_xxyyyy_0,  \
                             g_0_xxxxxxx_0_xxyyyy_1,  \
                             g_0_xxxxxxx_0_xxyyyz_0,  \
                             g_0_xxxxxxx_0_xxyyyz_1,  \
                             g_0_xxxxxxx_0_xxyyz_1,   \
                             g_0_xxxxxxx_0_xxyyzz_0,  \
                             g_0_xxxxxxx_0_xxyyzz_1,  \
                             g_0_xxxxxxx_0_xxyzz_1,   \
                             g_0_xxxxxxx_0_xxyzzz_0,  \
                             g_0_xxxxxxx_0_xxyzzz_1,  \
                             g_0_xxxxxxx_0_xxzzz_1,   \
                             g_0_xxxxxxx_0_xxzzzz_0,  \
                             g_0_xxxxxxx_0_xxzzzz_1,  \
                             g_0_xxxxxxx_0_xyyyy_1,   \
                             g_0_xxxxxxx_0_xyyyyy_0,  \
                             g_0_xxxxxxx_0_xyyyyy_1,  \
                             g_0_xxxxxxx_0_xyyyyz_0,  \
                             g_0_xxxxxxx_0_xyyyyz_1,  \
                             g_0_xxxxxxx_0_xyyyz_1,   \
                             g_0_xxxxxxx_0_xyyyzz_0,  \
                             g_0_xxxxxxx_0_xyyyzz_1,  \
                             g_0_xxxxxxx_0_xyyzz_1,   \
                             g_0_xxxxxxx_0_xyyzzz_0,  \
                             g_0_xxxxxxx_0_xyyzzz_1,  \
                             g_0_xxxxxxx_0_xyzzz_1,   \
                             g_0_xxxxxxx_0_xyzzzz_0,  \
                             g_0_xxxxxxx_0_xyzzzz_1,  \
                             g_0_xxxxxxx_0_xzzzz_1,   \
                             g_0_xxxxxxx_0_xzzzzz_0,  \
                             g_0_xxxxxxx_0_xzzzzz_1,  \
                             g_0_xxxxxxx_0_yyyyy_1,   \
                             g_0_xxxxxxx_0_yyyyyy_0,  \
                             g_0_xxxxxxx_0_yyyyyy_1,  \
                             g_0_xxxxxxx_0_yyyyyz_0,  \
                             g_0_xxxxxxx_0_yyyyyz_1,  \
                             g_0_xxxxxxx_0_yyyyz_1,   \
                             g_0_xxxxxxx_0_yyyyzz_0,  \
                             g_0_xxxxxxx_0_yyyyzz_1,  \
                             g_0_xxxxxxx_0_yyyzz_1,   \
                             g_0_xxxxxxx_0_yyyzzz_0,  \
                             g_0_xxxxxxx_0_yyyzzz_1,  \
                             g_0_xxxxxxx_0_yyzzz_1,   \
                             g_0_xxxxxxx_0_yyzzzz_0,  \
                             g_0_xxxxxxx_0_yyzzzz_1,  \
                             g_0_xxxxxxx_0_yzzzz_1,   \
                             g_0_xxxxxxx_0_yzzzzz_0,  \
                             g_0_xxxxxxx_0_yzzzzz_1,  \
                             g_0_xxxxxxx_0_zzzzz_1,   \
                             g_0_xxxxxxx_0_zzzzzz_0,  \
                             g_0_xxxxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxxxy_0_xxxxxx_0, \
                             g_0_xxxxxxxy_0_xxxxxy_0, \
                             g_0_xxxxxxxy_0_xxxxxz_0, \
                             g_0_xxxxxxxy_0_xxxxyy_0, \
                             g_0_xxxxxxxy_0_xxxxyz_0, \
                             g_0_xxxxxxxy_0_xxxxzz_0, \
                             g_0_xxxxxxxy_0_xxxyyy_0, \
                             g_0_xxxxxxxy_0_xxxyyz_0, \
                             g_0_xxxxxxxy_0_xxxyzz_0, \
                             g_0_xxxxxxxy_0_xxxzzz_0, \
                             g_0_xxxxxxxy_0_xxyyyy_0, \
                             g_0_xxxxxxxy_0_xxyyyz_0, \
                             g_0_xxxxxxxy_0_xxyyzz_0, \
                             g_0_xxxxxxxy_0_xxyzzz_0, \
                             g_0_xxxxxxxy_0_xxzzzz_0, \
                             g_0_xxxxxxxy_0_xyyyyy_0, \
                             g_0_xxxxxxxy_0_xyyyyz_0, \
                             g_0_xxxxxxxy_0_xyyyzz_0, \
                             g_0_xxxxxxxy_0_xyyzzz_0, \
                             g_0_xxxxxxxy_0_xyzzzz_0, \
                             g_0_xxxxxxxy_0_xzzzzz_0, \
                             g_0_xxxxxxxy_0_yyyyyy_0, \
                             g_0_xxxxxxxy_0_yyyyyz_0, \
                             g_0_xxxxxxxy_0_yyyyzz_0, \
                             g_0_xxxxxxxy_0_yyyzzz_0, \
                             g_0_xxxxxxxy_0_yyzzzz_0, \
                             g_0_xxxxxxxy_0_yzzzzz_0, \
                             g_0_xxxxxxxy_0_zzzzzz_0, \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_xxxxxx_0[i] = g_0_xxxxxxx_0_xxxxxx_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxy_0[i] = g_0_xxxxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxy_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxxz_0[i] = g_0_xxxxxxx_0_xxxxxz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxyy_0[i] =
            2.0 * g_0_xxxxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyy_0[i] * pb_y + g_0_xxxxxxx_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxyz_0[i] = g_0_xxxxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxzz_0[i] = g_0_xxxxxxx_0_xxxxzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyyy_0[i] =
            3.0 * g_0_xxxxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyy_0[i] * pb_y + g_0_xxxxxxx_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyyz_0[i] =
            2.0 * g_0_xxxxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyz_0[i] * pb_y + g_0_xxxxxxx_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyzz_0[i] = g_0_xxxxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxzzz_0[i] = g_0_xxxxxxx_0_xxxzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyyy_0[i] =
            4.0 * g_0_xxxxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyy_0[i] * pb_y + g_0_xxxxxxx_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyyz_0[i] =
            3.0 * g_0_xxxxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyz_0[i] * pb_y + g_0_xxxxxxx_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyzz_0[i] * pb_y + g_0_xxxxxxx_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyzzz_0[i] = g_0_xxxxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxzzzz_0[i] = g_0_xxxxxxx_0_xxzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyyy_0[i] =
            5.0 * g_0_xxxxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyy_0[i] * pb_y + g_0_xxxxxxx_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyyz_0[i] =
            4.0 * g_0_xxxxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyz_0[i] * pb_y + g_0_xxxxxxx_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyzz_0[i] =
            3.0 * g_0_xxxxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyzz_0[i] * pb_y + g_0_xxxxxxx_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyzzz_0[i] =
            2.0 * g_0_xxxxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzzz_0[i] * pb_y + g_0_xxxxxxx_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyzzzz_0[i] = g_0_xxxxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xzzzzz_0[i] = g_0_xxxxxxx_0_xzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyyy_0[i] =
            6.0 * g_0_xxxxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyy_0[i] * pb_y + g_0_xxxxxxx_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyyz_0[i] =
            5.0 * g_0_xxxxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyz_0[i] * pb_y + g_0_xxxxxxx_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyzz_0[i] =
            4.0 * g_0_xxxxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyzz_0[i] * pb_y + g_0_xxxxxxx_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyzzz_0[i] * pb_y + g_0_xxxxxxx_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyzzzz_0[i] =
            2.0 * g_0_xxxxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzzzz_0[i] * pb_y + g_0_xxxxxxx_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yzzzzz_0[i] = g_0_xxxxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_zzzzzz_0[i] = g_0_xxxxxxx_0_zzzzzz_0[i] * pb_y + g_0_xxxxxxx_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : SLSI

    auto g_0_xxxxxxxz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 56);

    auto g_0_xxxxxxxz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 57);

    auto g_0_xxxxxxxz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 58);

    auto g_0_xxxxxxxz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 59);

    auto g_0_xxxxxxxz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 60);

    auto g_0_xxxxxxxz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 61);

    auto g_0_xxxxxxxz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 62);

    auto g_0_xxxxxxxz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 63);

    auto g_0_xxxxxxxz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 64);

    auto g_0_xxxxxxxz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 65);

    auto g_0_xxxxxxxz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 66);

    auto g_0_xxxxxxxz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 67);

    auto g_0_xxxxxxxz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 68);

    auto g_0_xxxxxxxz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 69);

    auto g_0_xxxxxxxz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 70);

    auto g_0_xxxxxxxz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 71);

    auto g_0_xxxxxxxz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 72);

    auto g_0_xxxxxxxz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 73);

    auto g_0_xxxxxxxz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 74);

    auto g_0_xxxxxxxz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 75);

    auto g_0_xxxxxxxz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 76);

    auto g_0_xxxxxxxz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 77);

    auto g_0_xxxxxxxz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 78);

    auto g_0_xxxxxxxz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 79);

    auto g_0_xxxxxxxz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 80);

    auto g_0_xxxxxxxz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 81);

    auto g_0_xxxxxxxz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 82);

    auto g_0_xxxxxxxz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 83);

#pragma omp simd aligned(g_0_xxxxxxx_0_xxxxx_1,       \
                             g_0_xxxxxxx_0_xxxxxx_0,  \
                             g_0_xxxxxxx_0_xxxxxx_1,  \
                             g_0_xxxxxxx_0_xxxxxy_0,  \
                             g_0_xxxxxxx_0_xxxxxy_1,  \
                             g_0_xxxxxxx_0_xxxxxz_0,  \
                             g_0_xxxxxxx_0_xxxxxz_1,  \
                             g_0_xxxxxxx_0_xxxxy_1,   \
                             g_0_xxxxxxx_0_xxxxyy_0,  \
                             g_0_xxxxxxx_0_xxxxyy_1,  \
                             g_0_xxxxxxx_0_xxxxyz_0,  \
                             g_0_xxxxxxx_0_xxxxyz_1,  \
                             g_0_xxxxxxx_0_xxxxz_1,   \
                             g_0_xxxxxxx_0_xxxxzz_0,  \
                             g_0_xxxxxxx_0_xxxxzz_1,  \
                             g_0_xxxxxxx_0_xxxyy_1,   \
                             g_0_xxxxxxx_0_xxxyyy_0,  \
                             g_0_xxxxxxx_0_xxxyyy_1,  \
                             g_0_xxxxxxx_0_xxxyyz_0,  \
                             g_0_xxxxxxx_0_xxxyyz_1,  \
                             g_0_xxxxxxx_0_xxxyz_1,   \
                             g_0_xxxxxxx_0_xxxyzz_0,  \
                             g_0_xxxxxxx_0_xxxyzz_1,  \
                             g_0_xxxxxxx_0_xxxzz_1,   \
                             g_0_xxxxxxx_0_xxxzzz_0,  \
                             g_0_xxxxxxx_0_xxxzzz_1,  \
                             g_0_xxxxxxx_0_xxyyy_1,   \
                             g_0_xxxxxxx_0_xxyyyy_0,  \
                             g_0_xxxxxxx_0_xxyyyy_1,  \
                             g_0_xxxxxxx_0_xxyyyz_0,  \
                             g_0_xxxxxxx_0_xxyyyz_1,  \
                             g_0_xxxxxxx_0_xxyyz_1,   \
                             g_0_xxxxxxx_0_xxyyzz_0,  \
                             g_0_xxxxxxx_0_xxyyzz_1,  \
                             g_0_xxxxxxx_0_xxyzz_1,   \
                             g_0_xxxxxxx_0_xxyzzz_0,  \
                             g_0_xxxxxxx_0_xxyzzz_1,  \
                             g_0_xxxxxxx_0_xxzzz_1,   \
                             g_0_xxxxxxx_0_xxzzzz_0,  \
                             g_0_xxxxxxx_0_xxzzzz_1,  \
                             g_0_xxxxxxx_0_xyyyy_1,   \
                             g_0_xxxxxxx_0_xyyyyy_0,  \
                             g_0_xxxxxxx_0_xyyyyy_1,  \
                             g_0_xxxxxxx_0_xyyyyz_0,  \
                             g_0_xxxxxxx_0_xyyyyz_1,  \
                             g_0_xxxxxxx_0_xyyyz_1,   \
                             g_0_xxxxxxx_0_xyyyzz_0,  \
                             g_0_xxxxxxx_0_xyyyzz_1,  \
                             g_0_xxxxxxx_0_xyyzz_1,   \
                             g_0_xxxxxxx_0_xyyzzz_0,  \
                             g_0_xxxxxxx_0_xyyzzz_1,  \
                             g_0_xxxxxxx_0_xyzzz_1,   \
                             g_0_xxxxxxx_0_xyzzzz_0,  \
                             g_0_xxxxxxx_0_xyzzzz_1,  \
                             g_0_xxxxxxx_0_xzzzz_1,   \
                             g_0_xxxxxxx_0_xzzzzz_0,  \
                             g_0_xxxxxxx_0_xzzzzz_1,  \
                             g_0_xxxxxxx_0_yyyyy_1,   \
                             g_0_xxxxxxx_0_yyyyyy_0,  \
                             g_0_xxxxxxx_0_yyyyyy_1,  \
                             g_0_xxxxxxx_0_yyyyyz_0,  \
                             g_0_xxxxxxx_0_yyyyyz_1,  \
                             g_0_xxxxxxx_0_yyyyz_1,   \
                             g_0_xxxxxxx_0_yyyyzz_0,  \
                             g_0_xxxxxxx_0_yyyyzz_1,  \
                             g_0_xxxxxxx_0_yyyzz_1,   \
                             g_0_xxxxxxx_0_yyyzzz_0,  \
                             g_0_xxxxxxx_0_yyyzzz_1,  \
                             g_0_xxxxxxx_0_yyzzz_1,   \
                             g_0_xxxxxxx_0_yyzzzz_0,  \
                             g_0_xxxxxxx_0_yyzzzz_1,  \
                             g_0_xxxxxxx_0_yzzzz_1,   \
                             g_0_xxxxxxx_0_yzzzzz_0,  \
                             g_0_xxxxxxx_0_yzzzzz_1,  \
                             g_0_xxxxxxx_0_zzzzz_1,   \
                             g_0_xxxxxxx_0_zzzzzz_0,  \
                             g_0_xxxxxxx_0_zzzzzz_1,  \
                             g_0_xxxxxxxz_0_xxxxxx_0, \
                             g_0_xxxxxxxz_0_xxxxxy_0, \
                             g_0_xxxxxxxz_0_xxxxxz_0, \
                             g_0_xxxxxxxz_0_xxxxyy_0, \
                             g_0_xxxxxxxz_0_xxxxyz_0, \
                             g_0_xxxxxxxz_0_xxxxzz_0, \
                             g_0_xxxxxxxz_0_xxxyyy_0, \
                             g_0_xxxxxxxz_0_xxxyyz_0, \
                             g_0_xxxxxxxz_0_xxxyzz_0, \
                             g_0_xxxxxxxz_0_xxxzzz_0, \
                             g_0_xxxxxxxz_0_xxyyyy_0, \
                             g_0_xxxxxxxz_0_xxyyyz_0, \
                             g_0_xxxxxxxz_0_xxyyzz_0, \
                             g_0_xxxxxxxz_0_xxyzzz_0, \
                             g_0_xxxxxxxz_0_xxzzzz_0, \
                             g_0_xxxxxxxz_0_xyyyyy_0, \
                             g_0_xxxxxxxz_0_xyyyyz_0, \
                             g_0_xxxxxxxz_0_xyyyzz_0, \
                             g_0_xxxxxxxz_0_xyyzzz_0, \
                             g_0_xxxxxxxz_0_xyzzzz_0, \
                             g_0_xxxxxxxz_0_xzzzzz_0, \
                             g_0_xxxxxxxz_0_yyyyyy_0, \
                             g_0_xxxxxxxz_0_yyyyyz_0, \
                             g_0_xxxxxxxz_0_yyyyzz_0, \
                             g_0_xxxxxxxz_0_yyyzzz_0, \
                             g_0_xxxxxxxz_0_yyzzzz_0, \
                             g_0_xxxxxxxz_0_yzzzzz_0, \
                             g_0_xxxxxxxz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_xxxxxx_0[i] = g_0_xxxxxxx_0_xxxxxx_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxy_0[i] = g_0_xxxxxxx_0_xxxxxy_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxxz_0[i] = g_0_xxxxxxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxxz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxyy_0[i] = g_0_xxxxxxx_0_xxxxyy_0[i] * pb_z + g_0_xxxxxxx_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxyz_0[i] = g_0_xxxxxxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxyz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxzz_0[i] =
            2.0 * g_0_xxxxxxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyyy_0[i] = g_0_xxxxxxx_0_xxxyyy_0[i] * pb_z + g_0_xxxxxxx_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyyz_0[i] = g_0_xxxxxxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyyz_0[i] * pb_z + g_0_xxxxxxx_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyzz_0[i] =
            2.0 * g_0_xxxxxxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxzzz_0[i] =
            3.0 * g_0_xxxxxxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyyy_0[i] = g_0_xxxxxxx_0_xxyyyy_0[i] * pb_z + g_0_xxxxxxx_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyyz_0[i] = g_0_xxxxxxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyyz_0[i] * pb_z + g_0_xxxxxxx_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxxxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyzz_0[i] * pb_z + g_0_xxxxxxx_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyzzz_0[i] =
            3.0 * g_0_xxxxxxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxzzzz_0[i] =
            4.0 * g_0_xxxxxxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyyy_0[i] = g_0_xxxxxxx_0_xyyyyy_0[i] * pb_z + g_0_xxxxxxx_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyyz_0[i] = g_0_xxxxxxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyyz_0[i] * pb_z + g_0_xxxxxxx_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyzz_0[i] =
            2.0 * g_0_xxxxxxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyzz_0[i] * pb_z + g_0_xxxxxxx_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyzzz_0[i] =
            3.0 * g_0_xxxxxxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzzz_0[i] * pb_z + g_0_xxxxxxx_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyzzzz_0[i] =
            4.0 * g_0_xxxxxxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xzzzzz_0[i] =
            5.0 * g_0_xxxxxxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyyy_0[i] = g_0_xxxxxxx_0_yyyyyy_0[i] * pb_z + g_0_xxxxxxx_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyyz_0[i] = g_0_xxxxxxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyyz_0[i] * pb_z + g_0_xxxxxxx_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyzz_0[i] =
            2.0 * g_0_xxxxxxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyzz_0[i] * pb_z + g_0_xxxxxxx_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxxxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyzzz_0[i] * pb_z + g_0_xxxxxxx_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyzzzz_0[i] =
            4.0 * g_0_xxxxxxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzzzz_0[i] * pb_z + g_0_xxxxxxx_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yzzzzz_0[i] =
            5.0 * g_0_xxxxxxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_zzzzzz_0[i] =
            6.0 * g_0_xxxxxxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_zzzzzz_0[i] * pb_z + g_0_xxxxxxx_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 84-112 components of targeted buffer : SLSI

    auto g_0_xxxxxxyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 84);

    auto g_0_xxxxxxyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 85);

    auto g_0_xxxxxxyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 86);

    auto g_0_xxxxxxyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 87);

    auto g_0_xxxxxxyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 88);

    auto g_0_xxxxxxyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 89);

    auto g_0_xxxxxxyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 90);

    auto g_0_xxxxxxyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 91);

    auto g_0_xxxxxxyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 92);

    auto g_0_xxxxxxyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 93);

    auto g_0_xxxxxxyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 94);

    auto g_0_xxxxxxyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 95);

    auto g_0_xxxxxxyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 96);

    auto g_0_xxxxxxyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 97);

    auto g_0_xxxxxxyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 98);

    auto g_0_xxxxxxyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 99);

    auto g_0_xxxxxxyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 100);

    auto g_0_xxxxxxyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 101);

    auto g_0_xxxxxxyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 102);

    auto g_0_xxxxxxyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 103);

    auto g_0_xxxxxxyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 104);

    auto g_0_xxxxxxyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 105);

    auto g_0_xxxxxxyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 106);

    auto g_0_xxxxxxyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 107);

    auto g_0_xxxxxxyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 108);

    auto g_0_xxxxxxyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 109);

    auto g_0_xxxxxxyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 110);

    auto g_0_xxxxxxyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 111);

#pragma omp simd aligned(g_0_xxxxxx_0_xxxxxx_0,       \
                             g_0_xxxxxx_0_xxxxxx_1,   \
                             g_0_xxxxxx_0_xxxxxz_0,   \
                             g_0_xxxxxx_0_xxxxxz_1,   \
                             g_0_xxxxxx_0_xxxxzz_0,   \
                             g_0_xxxxxx_0_xxxxzz_1,   \
                             g_0_xxxxxx_0_xxxzzz_0,   \
                             g_0_xxxxxx_0_xxxzzz_1,   \
                             g_0_xxxxxx_0_xxzzzz_0,   \
                             g_0_xxxxxx_0_xxzzzz_1,   \
                             g_0_xxxxxx_0_xzzzzz_0,   \
                             g_0_xxxxxx_0_xzzzzz_1,   \
                             g_0_xxxxxxy_0_xxxxxx_0,  \
                             g_0_xxxxxxy_0_xxxxxx_1,  \
                             g_0_xxxxxxy_0_xxxxxz_0,  \
                             g_0_xxxxxxy_0_xxxxxz_1,  \
                             g_0_xxxxxxy_0_xxxxzz_0,  \
                             g_0_xxxxxxy_0_xxxxzz_1,  \
                             g_0_xxxxxxy_0_xxxzzz_0,  \
                             g_0_xxxxxxy_0_xxxzzz_1,  \
                             g_0_xxxxxxy_0_xxzzzz_0,  \
                             g_0_xxxxxxy_0_xxzzzz_1,  \
                             g_0_xxxxxxy_0_xzzzzz_0,  \
                             g_0_xxxxxxy_0_xzzzzz_1,  \
                             g_0_xxxxxxyy_0_xxxxxx_0, \
                             g_0_xxxxxxyy_0_xxxxxy_0, \
                             g_0_xxxxxxyy_0_xxxxxz_0, \
                             g_0_xxxxxxyy_0_xxxxyy_0, \
                             g_0_xxxxxxyy_0_xxxxyz_0, \
                             g_0_xxxxxxyy_0_xxxxzz_0, \
                             g_0_xxxxxxyy_0_xxxyyy_0, \
                             g_0_xxxxxxyy_0_xxxyyz_0, \
                             g_0_xxxxxxyy_0_xxxyzz_0, \
                             g_0_xxxxxxyy_0_xxxzzz_0, \
                             g_0_xxxxxxyy_0_xxyyyy_0, \
                             g_0_xxxxxxyy_0_xxyyyz_0, \
                             g_0_xxxxxxyy_0_xxyyzz_0, \
                             g_0_xxxxxxyy_0_xxyzzz_0, \
                             g_0_xxxxxxyy_0_xxzzzz_0, \
                             g_0_xxxxxxyy_0_xyyyyy_0, \
                             g_0_xxxxxxyy_0_xyyyyz_0, \
                             g_0_xxxxxxyy_0_xyyyzz_0, \
                             g_0_xxxxxxyy_0_xyyzzz_0, \
                             g_0_xxxxxxyy_0_xyzzzz_0, \
                             g_0_xxxxxxyy_0_xzzzzz_0, \
                             g_0_xxxxxxyy_0_yyyyyy_0, \
                             g_0_xxxxxxyy_0_yyyyyz_0, \
                             g_0_xxxxxxyy_0_yyyyzz_0, \
                             g_0_xxxxxxyy_0_yyyzzz_0, \
                             g_0_xxxxxxyy_0_yyzzzz_0, \
                             g_0_xxxxxxyy_0_yzzzzz_0, \
                             g_0_xxxxxxyy_0_zzzzzz_0, \
                             g_0_xxxxxyy_0_xxxxxy_0,  \
                             g_0_xxxxxyy_0_xxxxxy_1,  \
                             g_0_xxxxxyy_0_xxxxy_1,   \
                             g_0_xxxxxyy_0_xxxxyy_0,  \
                             g_0_xxxxxyy_0_xxxxyy_1,  \
                             g_0_xxxxxyy_0_xxxxyz_0,  \
                             g_0_xxxxxyy_0_xxxxyz_1,  \
                             g_0_xxxxxyy_0_xxxyy_1,   \
                             g_0_xxxxxyy_0_xxxyyy_0,  \
                             g_0_xxxxxyy_0_xxxyyy_1,  \
                             g_0_xxxxxyy_0_xxxyyz_0,  \
                             g_0_xxxxxyy_0_xxxyyz_1,  \
                             g_0_xxxxxyy_0_xxxyz_1,   \
                             g_0_xxxxxyy_0_xxxyzz_0,  \
                             g_0_xxxxxyy_0_xxxyzz_1,  \
                             g_0_xxxxxyy_0_xxyyy_1,   \
                             g_0_xxxxxyy_0_xxyyyy_0,  \
                             g_0_xxxxxyy_0_xxyyyy_1,  \
                             g_0_xxxxxyy_0_xxyyyz_0,  \
                             g_0_xxxxxyy_0_xxyyyz_1,  \
                             g_0_xxxxxyy_0_xxyyz_1,   \
                             g_0_xxxxxyy_0_xxyyzz_0,  \
                             g_0_xxxxxyy_0_xxyyzz_1,  \
                             g_0_xxxxxyy_0_xxyzz_1,   \
                             g_0_xxxxxyy_0_xxyzzz_0,  \
                             g_0_xxxxxyy_0_xxyzzz_1,  \
                             g_0_xxxxxyy_0_xyyyy_1,   \
                             g_0_xxxxxyy_0_xyyyyy_0,  \
                             g_0_xxxxxyy_0_xyyyyy_1,  \
                             g_0_xxxxxyy_0_xyyyyz_0,  \
                             g_0_xxxxxyy_0_xyyyyz_1,  \
                             g_0_xxxxxyy_0_xyyyz_1,   \
                             g_0_xxxxxyy_0_xyyyzz_0,  \
                             g_0_xxxxxyy_0_xyyyzz_1,  \
                             g_0_xxxxxyy_0_xyyzz_1,   \
                             g_0_xxxxxyy_0_xyyzzz_0,  \
                             g_0_xxxxxyy_0_xyyzzz_1,  \
                             g_0_xxxxxyy_0_xyzzz_1,   \
                             g_0_xxxxxyy_0_xyzzzz_0,  \
                             g_0_xxxxxyy_0_xyzzzz_1,  \
                             g_0_xxxxxyy_0_yyyyy_1,   \
                             g_0_xxxxxyy_0_yyyyyy_0,  \
                             g_0_xxxxxyy_0_yyyyyy_1,  \
                             g_0_xxxxxyy_0_yyyyyz_0,  \
                             g_0_xxxxxyy_0_yyyyyz_1,  \
                             g_0_xxxxxyy_0_yyyyz_1,   \
                             g_0_xxxxxyy_0_yyyyzz_0,  \
                             g_0_xxxxxyy_0_yyyyzz_1,  \
                             g_0_xxxxxyy_0_yyyzz_1,   \
                             g_0_xxxxxyy_0_yyyzzz_0,  \
                             g_0_xxxxxyy_0_yyyzzz_1,  \
                             g_0_xxxxxyy_0_yyzzz_1,   \
                             g_0_xxxxxyy_0_yyzzzz_0,  \
                             g_0_xxxxxyy_0_yyzzzz_1,  \
                             g_0_xxxxxyy_0_yzzzz_1,   \
                             g_0_xxxxxyy_0_yzzzzz_0,  \
                             g_0_xxxxxyy_0_yzzzzz_1,  \
                             g_0_xxxxxyy_0_zzzzzz_0,  \
                             g_0_xxxxxyy_0_zzzzzz_1,  \
                             g_0_xxxxyy_0_xxxxxy_0,   \
                             g_0_xxxxyy_0_xxxxxy_1,   \
                             g_0_xxxxyy_0_xxxxyy_0,   \
                             g_0_xxxxyy_0_xxxxyy_1,   \
                             g_0_xxxxyy_0_xxxxyz_0,   \
                             g_0_xxxxyy_0_xxxxyz_1,   \
                             g_0_xxxxyy_0_xxxyyy_0,   \
                             g_0_xxxxyy_0_xxxyyy_1,   \
                             g_0_xxxxyy_0_xxxyyz_0,   \
                             g_0_xxxxyy_0_xxxyyz_1,   \
                             g_0_xxxxyy_0_xxxyzz_0,   \
                             g_0_xxxxyy_0_xxxyzz_1,   \
                             g_0_xxxxyy_0_xxyyyy_0,   \
                             g_0_xxxxyy_0_xxyyyy_1,   \
                             g_0_xxxxyy_0_xxyyyz_0,   \
                             g_0_xxxxyy_0_xxyyyz_1,   \
                             g_0_xxxxyy_0_xxyyzz_0,   \
                             g_0_xxxxyy_0_xxyyzz_1,   \
                             g_0_xxxxyy_0_xxyzzz_0,   \
                             g_0_xxxxyy_0_xxyzzz_1,   \
                             g_0_xxxxyy_0_xyyyyy_0,   \
                             g_0_xxxxyy_0_xyyyyy_1,   \
                             g_0_xxxxyy_0_xyyyyz_0,   \
                             g_0_xxxxyy_0_xyyyyz_1,   \
                             g_0_xxxxyy_0_xyyyzz_0,   \
                             g_0_xxxxyy_0_xyyyzz_1,   \
                             g_0_xxxxyy_0_xyyzzz_0,   \
                             g_0_xxxxyy_0_xyyzzz_1,   \
                             g_0_xxxxyy_0_xyzzzz_0,   \
                             g_0_xxxxyy_0_xyzzzz_1,   \
                             g_0_xxxxyy_0_yyyyyy_0,   \
                             g_0_xxxxyy_0_yyyyyy_1,   \
                             g_0_xxxxyy_0_yyyyyz_0,   \
                             g_0_xxxxyy_0_yyyyyz_1,   \
                             g_0_xxxxyy_0_yyyyzz_0,   \
                             g_0_xxxxyy_0_yyyyzz_1,   \
                             g_0_xxxxyy_0_yyyzzz_0,   \
                             g_0_xxxxyy_0_yyyzzz_1,   \
                             g_0_xxxxyy_0_yyzzzz_0,   \
                             g_0_xxxxyy_0_yyzzzz_1,   \
                             g_0_xxxxyy_0_yzzzzz_0,   \
                             g_0_xxxxyy_0_yzzzzz_1,   \
                             g_0_xxxxyy_0_zzzzzz_0,   \
                             g_0_xxxxyy_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_xxxxxx_0[i] = g_0_xxxxxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxxx_0[i] * pb_y +
                                     g_0_xxxxxxy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxxxy_0[i] = 5.0 * g_0_xxxxyy_0_xxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxxy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxy_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxxz_0[i] = g_0_xxxxxx_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxxz_0[i] * pb_y +
                                     g_0_xxxxxxy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxxyy_0[i] = 5.0 * g_0_xxxxyy_0_xxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyy_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxyz_0[i] = 5.0 * g_0_xxxxyy_0_xxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxzz_0[i] = g_0_xxxxxx_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxzz_0[i] * pb_y +
                                     g_0_xxxxxxy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxyyy_0[i] = 5.0 * g_0_xxxxyy_0_xxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyy_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxyyz_0[i] = 5.0 * g_0_xxxxyy_0_xxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxyzz_0[i] = 5.0 * g_0_xxxxyy_0_xxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxzzz_0[i] = g_0_xxxxxx_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxzzz_0[i] * pb_y +
                                     g_0_xxxxxxy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxyyyy_0[i] = 5.0 * g_0_xxxxyy_0_xxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyy_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyyyz_0[i] = 5.0 * g_0_xxxxyy_0_xxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyyzz_0[i] = 5.0 * g_0_xxxxyy_0_xxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyzzz_0[i] = 5.0 * g_0_xxxxyy_0_xxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxxxyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxzzzz_0[i] = g_0_xxxxxx_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxzzzz_0[i] * pb_y +
                                     g_0_xxxxxxy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xyyyyy_0[i] = 5.0 * g_0_xxxxyy_0_xyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyy_0[i] * pb_x + g_0_xxxxxyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyyyz_0[i] = 5.0 * g_0_xxxxyy_0_xyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyz_0[i] * pb_x + g_0_xxxxxyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyyzz_0[i] = 5.0 * g_0_xxxxyy_0_xyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyzz_0[i] * pb_x + g_0_xxxxxyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyzzz_0[i] * pb_x + g_0_xxxxxyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyzzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzzzz_0[i] * pb_x + g_0_xxxxxyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xzzzzz_0[i] = g_0_xxxxxx_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xzzzzz_0[i] * pb_y +
                                     g_0_xxxxxxy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_yyyyyy_0[i] = 5.0 * g_0_xxxxyy_0_yyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyyyy_0[i] * pb_x + g_0_xxxxxyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyyyz_0[i] = 5.0 * g_0_xxxxyy_0_yyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyyyz_0[i] * pb_x + g_0_xxxxxyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyyzz_0[i] = 5.0 * g_0_xxxxyy_0_yyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyyzz_0[i] * pb_x + g_0_xxxxxyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyyzzz_0[i] * pb_x + g_0_xxxxxyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yyzzzz_0[i] * pb_x + g_0_xxxxxyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_yzzzzz_0[i] * pb_x + g_0_xxxxxyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_zzzzzz_0[i] = 5.0 * g_0_xxxxyy_0_zzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_zzzzzz_0[i] * pb_x + g_0_xxxxxyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 112-140 components of targeted buffer : SLSI

    auto g_0_xxxxxxyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 112);

    auto g_0_xxxxxxyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 113);

    auto g_0_xxxxxxyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 114);

    auto g_0_xxxxxxyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 115);

    auto g_0_xxxxxxyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 116);

    auto g_0_xxxxxxyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 117);

    auto g_0_xxxxxxyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 118);

    auto g_0_xxxxxxyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 119);

    auto g_0_xxxxxxyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 120);

    auto g_0_xxxxxxyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 121);

    auto g_0_xxxxxxyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 122);

    auto g_0_xxxxxxyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 123);

    auto g_0_xxxxxxyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 124);

    auto g_0_xxxxxxyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 125);

    auto g_0_xxxxxxyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 126);

    auto g_0_xxxxxxyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 127);

    auto g_0_xxxxxxyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 128);

    auto g_0_xxxxxxyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 129);

    auto g_0_xxxxxxyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 130);

    auto g_0_xxxxxxyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 131);

    auto g_0_xxxxxxyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 132);

    auto g_0_xxxxxxyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 133);

    auto g_0_xxxxxxyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 134);

    auto g_0_xxxxxxyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 135);

    auto g_0_xxxxxxyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 136);

    auto g_0_xxxxxxyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 137);

    auto g_0_xxxxxxyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 138);

    auto g_0_xxxxxxyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 139);

#pragma omp simd aligned(g_0_xxxxxxy_0_xxxxxy_0,      \
                             g_0_xxxxxxy_0_xxxxxy_1,  \
                             g_0_xxxxxxy_0_xxxxyy_0,  \
                             g_0_xxxxxxy_0_xxxxyy_1,  \
                             g_0_xxxxxxy_0_xxxyyy_0,  \
                             g_0_xxxxxxy_0_xxxyyy_1,  \
                             g_0_xxxxxxy_0_xxyyyy_0,  \
                             g_0_xxxxxxy_0_xxyyyy_1,  \
                             g_0_xxxxxxy_0_xyyyyy_0,  \
                             g_0_xxxxxxy_0_xyyyyy_1,  \
                             g_0_xxxxxxy_0_yyyyyy_0,  \
                             g_0_xxxxxxy_0_yyyyyy_1,  \
                             g_0_xxxxxxyz_0_xxxxxx_0, \
                             g_0_xxxxxxyz_0_xxxxxy_0, \
                             g_0_xxxxxxyz_0_xxxxxz_0, \
                             g_0_xxxxxxyz_0_xxxxyy_0, \
                             g_0_xxxxxxyz_0_xxxxyz_0, \
                             g_0_xxxxxxyz_0_xxxxzz_0, \
                             g_0_xxxxxxyz_0_xxxyyy_0, \
                             g_0_xxxxxxyz_0_xxxyyz_0, \
                             g_0_xxxxxxyz_0_xxxyzz_0, \
                             g_0_xxxxxxyz_0_xxxzzz_0, \
                             g_0_xxxxxxyz_0_xxyyyy_0, \
                             g_0_xxxxxxyz_0_xxyyyz_0, \
                             g_0_xxxxxxyz_0_xxyyzz_0, \
                             g_0_xxxxxxyz_0_xxyzzz_0, \
                             g_0_xxxxxxyz_0_xxzzzz_0, \
                             g_0_xxxxxxyz_0_xyyyyy_0, \
                             g_0_xxxxxxyz_0_xyyyyz_0, \
                             g_0_xxxxxxyz_0_xyyyzz_0, \
                             g_0_xxxxxxyz_0_xyyzzz_0, \
                             g_0_xxxxxxyz_0_xyzzzz_0, \
                             g_0_xxxxxxyz_0_xzzzzz_0, \
                             g_0_xxxxxxyz_0_yyyyyy_0, \
                             g_0_xxxxxxyz_0_yyyyyz_0, \
                             g_0_xxxxxxyz_0_yyyyzz_0, \
                             g_0_xxxxxxyz_0_yyyzzz_0, \
                             g_0_xxxxxxyz_0_yyzzzz_0, \
                             g_0_xxxxxxyz_0_yzzzzz_0, \
                             g_0_xxxxxxyz_0_zzzzzz_0, \
                             g_0_xxxxxxz_0_xxxxxx_0,  \
                             g_0_xxxxxxz_0_xxxxxx_1,  \
                             g_0_xxxxxxz_0_xxxxxz_0,  \
                             g_0_xxxxxxz_0_xxxxxz_1,  \
                             g_0_xxxxxxz_0_xxxxyz_0,  \
                             g_0_xxxxxxz_0_xxxxyz_1,  \
                             g_0_xxxxxxz_0_xxxxz_1,   \
                             g_0_xxxxxxz_0_xxxxzz_0,  \
                             g_0_xxxxxxz_0_xxxxzz_1,  \
                             g_0_xxxxxxz_0_xxxyyz_0,  \
                             g_0_xxxxxxz_0_xxxyyz_1,  \
                             g_0_xxxxxxz_0_xxxyz_1,   \
                             g_0_xxxxxxz_0_xxxyzz_0,  \
                             g_0_xxxxxxz_0_xxxyzz_1,  \
                             g_0_xxxxxxz_0_xxxzz_1,   \
                             g_0_xxxxxxz_0_xxxzzz_0,  \
                             g_0_xxxxxxz_0_xxxzzz_1,  \
                             g_0_xxxxxxz_0_xxyyyz_0,  \
                             g_0_xxxxxxz_0_xxyyyz_1,  \
                             g_0_xxxxxxz_0_xxyyz_1,   \
                             g_0_xxxxxxz_0_xxyyzz_0,  \
                             g_0_xxxxxxz_0_xxyyzz_1,  \
                             g_0_xxxxxxz_0_xxyzz_1,   \
                             g_0_xxxxxxz_0_xxyzzz_0,  \
                             g_0_xxxxxxz_0_xxyzzz_1,  \
                             g_0_xxxxxxz_0_xxzzz_1,   \
                             g_0_xxxxxxz_0_xxzzzz_0,  \
                             g_0_xxxxxxz_0_xxzzzz_1,  \
                             g_0_xxxxxxz_0_xyyyyz_0,  \
                             g_0_xxxxxxz_0_xyyyyz_1,  \
                             g_0_xxxxxxz_0_xyyyz_1,   \
                             g_0_xxxxxxz_0_xyyyzz_0,  \
                             g_0_xxxxxxz_0_xyyyzz_1,  \
                             g_0_xxxxxxz_0_xyyzz_1,   \
                             g_0_xxxxxxz_0_xyyzzz_0,  \
                             g_0_xxxxxxz_0_xyyzzz_1,  \
                             g_0_xxxxxxz_0_xyzzz_1,   \
                             g_0_xxxxxxz_0_xyzzzz_0,  \
                             g_0_xxxxxxz_0_xyzzzz_1,  \
                             g_0_xxxxxxz_0_xzzzz_1,   \
                             g_0_xxxxxxz_0_xzzzzz_0,  \
                             g_0_xxxxxxz_0_xzzzzz_1,  \
                             g_0_xxxxxxz_0_yyyyyz_0,  \
                             g_0_xxxxxxz_0_yyyyyz_1,  \
                             g_0_xxxxxxz_0_yyyyz_1,   \
                             g_0_xxxxxxz_0_yyyyzz_0,  \
                             g_0_xxxxxxz_0_yyyyzz_1,  \
                             g_0_xxxxxxz_0_yyyzz_1,   \
                             g_0_xxxxxxz_0_yyyzzz_0,  \
                             g_0_xxxxxxz_0_yyyzzz_1,  \
                             g_0_xxxxxxz_0_yyzzz_1,   \
                             g_0_xxxxxxz_0_yyzzzz_0,  \
                             g_0_xxxxxxz_0_yyzzzz_1,  \
                             g_0_xxxxxxz_0_yzzzz_1,   \
                             g_0_xxxxxxz_0_yzzzzz_0,  \
                             g_0_xxxxxxz_0_yzzzzz_1,  \
                             g_0_xxxxxxz_0_zzzzz_1,   \
                             g_0_xxxxxxz_0_zzzzzz_0,  \
                             g_0_xxxxxxz_0_zzzzzz_1,  \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyz_0_xxxxxx_0[i] = g_0_xxxxxxz_0_xxxxxx_0[i] * pb_y + g_0_xxxxxxz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxxy_0[i] = g_0_xxxxxxy_0_xxxxxy_0[i] * pb_z + g_0_xxxxxxy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxxxz_0[i] = g_0_xxxxxxz_0_xxxxxz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxyy_0[i] = g_0_xxxxxxy_0_xxxxyy_0[i] * pb_z + g_0_xxxxxxy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxxyz_0[i] = g_0_xxxxxxz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxxyz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxzz_0[i] = g_0_xxxxxxz_0_xxxxzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxyyy_0[i] = g_0_xxxxxxy_0_xxxyyy_0[i] * pb_z + g_0_xxxxxxy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxyyz_0[i] =
            2.0 * g_0_xxxxxxz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxyyz_0[i] * pb_y + g_0_xxxxxxz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxyzz_0[i] = g_0_xxxxxxz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxyzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxzzz_0[i] = g_0_xxxxxxz_0_xxxzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyyyy_0[i] = g_0_xxxxxxy_0_xxyyyy_0[i] * pb_z + g_0_xxxxxxy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxyyyz_0[i] =
            3.0 * g_0_xxxxxxz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyyyz_0[i] * pb_y + g_0_xxxxxxz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxxxz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyyzz_0[i] * pb_y + g_0_xxxxxxz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyzzz_0[i] = g_0_xxxxxxz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxzzzz_0[i] = g_0_xxxxxxz_0_xxzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyyyy_0[i] = g_0_xxxxxxy_0_xyyyyy_0[i] * pb_z + g_0_xxxxxxy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xyyyyz_0[i] =
            4.0 * g_0_xxxxxxz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyyyz_0[i] * pb_y + g_0_xxxxxxz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyyzz_0[i] =
            3.0 * g_0_xxxxxxz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyyzz_0[i] * pb_y + g_0_xxxxxxz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyzzz_0[i] =
            2.0 * g_0_xxxxxxz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyzzz_0[i] * pb_y + g_0_xxxxxxz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyzzzz_0[i] = g_0_xxxxxxz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xzzzzz_0[i] = g_0_xxxxxxz_0_xzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyyyy_0[i] = g_0_xxxxxxy_0_yyyyyy_0[i] * pb_z + g_0_xxxxxxy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_yyyyyz_0[i] =
            5.0 * g_0_xxxxxxz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyyyz_0[i] * pb_y + g_0_xxxxxxz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyyzz_0[i] =
            4.0 * g_0_xxxxxxz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyyzz_0[i] * pb_y + g_0_xxxxxxz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxxxz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyzzz_0[i] * pb_y + g_0_xxxxxxz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyzzzz_0[i] =
            2.0 * g_0_xxxxxxz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyzzzz_0[i] * pb_y + g_0_xxxxxxz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yzzzzz_0[i] = g_0_xxxxxxz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_zzzzzz_0[i] = g_0_xxxxxxz_0_zzzzzz_0[i] * pb_y + g_0_xxxxxxz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 140-168 components of targeted buffer : SLSI

    auto g_0_xxxxxxzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 140);

    auto g_0_xxxxxxzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 141);

    auto g_0_xxxxxxzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 142);

    auto g_0_xxxxxxzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 143);

    auto g_0_xxxxxxzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 144);

    auto g_0_xxxxxxzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 145);

    auto g_0_xxxxxxzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 146);

    auto g_0_xxxxxxzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 147);

    auto g_0_xxxxxxzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 148);

    auto g_0_xxxxxxzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 149);

    auto g_0_xxxxxxzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 150);

    auto g_0_xxxxxxzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 151);

    auto g_0_xxxxxxzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 152);

    auto g_0_xxxxxxzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 153);

    auto g_0_xxxxxxzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 154);

    auto g_0_xxxxxxzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 155);

    auto g_0_xxxxxxzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 156);

    auto g_0_xxxxxxzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 157);

    auto g_0_xxxxxxzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 158);

    auto g_0_xxxxxxzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 159);

    auto g_0_xxxxxxzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 160);

    auto g_0_xxxxxxzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 161);

    auto g_0_xxxxxxzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 162);

    auto g_0_xxxxxxzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 163);

    auto g_0_xxxxxxzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 164);

    auto g_0_xxxxxxzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 165);

    auto g_0_xxxxxxzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 166);

    auto g_0_xxxxxxzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 167);

#pragma omp simd aligned(g_0_xxxxxx_0_xxxxxx_0,       \
                             g_0_xxxxxx_0_xxxxxx_1,   \
                             g_0_xxxxxx_0_xxxxxy_0,   \
                             g_0_xxxxxx_0_xxxxxy_1,   \
                             g_0_xxxxxx_0_xxxxyy_0,   \
                             g_0_xxxxxx_0_xxxxyy_1,   \
                             g_0_xxxxxx_0_xxxyyy_0,   \
                             g_0_xxxxxx_0_xxxyyy_1,   \
                             g_0_xxxxxx_0_xxyyyy_0,   \
                             g_0_xxxxxx_0_xxyyyy_1,   \
                             g_0_xxxxxx_0_xyyyyy_0,   \
                             g_0_xxxxxx_0_xyyyyy_1,   \
                             g_0_xxxxxxz_0_xxxxxx_0,  \
                             g_0_xxxxxxz_0_xxxxxx_1,  \
                             g_0_xxxxxxz_0_xxxxxy_0,  \
                             g_0_xxxxxxz_0_xxxxxy_1,  \
                             g_0_xxxxxxz_0_xxxxyy_0,  \
                             g_0_xxxxxxz_0_xxxxyy_1,  \
                             g_0_xxxxxxz_0_xxxyyy_0,  \
                             g_0_xxxxxxz_0_xxxyyy_1,  \
                             g_0_xxxxxxz_0_xxyyyy_0,  \
                             g_0_xxxxxxz_0_xxyyyy_1,  \
                             g_0_xxxxxxz_0_xyyyyy_0,  \
                             g_0_xxxxxxz_0_xyyyyy_1,  \
                             g_0_xxxxxxzz_0_xxxxxx_0, \
                             g_0_xxxxxxzz_0_xxxxxy_0, \
                             g_0_xxxxxxzz_0_xxxxxz_0, \
                             g_0_xxxxxxzz_0_xxxxyy_0, \
                             g_0_xxxxxxzz_0_xxxxyz_0, \
                             g_0_xxxxxxzz_0_xxxxzz_0, \
                             g_0_xxxxxxzz_0_xxxyyy_0, \
                             g_0_xxxxxxzz_0_xxxyyz_0, \
                             g_0_xxxxxxzz_0_xxxyzz_0, \
                             g_0_xxxxxxzz_0_xxxzzz_0, \
                             g_0_xxxxxxzz_0_xxyyyy_0, \
                             g_0_xxxxxxzz_0_xxyyyz_0, \
                             g_0_xxxxxxzz_0_xxyyzz_0, \
                             g_0_xxxxxxzz_0_xxyzzz_0, \
                             g_0_xxxxxxzz_0_xxzzzz_0, \
                             g_0_xxxxxxzz_0_xyyyyy_0, \
                             g_0_xxxxxxzz_0_xyyyyz_0, \
                             g_0_xxxxxxzz_0_xyyyzz_0, \
                             g_0_xxxxxxzz_0_xyyzzz_0, \
                             g_0_xxxxxxzz_0_xyzzzz_0, \
                             g_0_xxxxxxzz_0_xzzzzz_0, \
                             g_0_xxxxxxzz_0_yyyyyy_0, \
                             g_0_xxxxxxzz_0_yyyyyz_0, \
                             g_0_xxxxxxzz_0_yyyyzz_0, \
                             g_0_xxxxxxzz_0_yyyzzz_0, \
                             g_0_xxxxxxzz_0_yyzzzz_0, \
                             g_0_xxxxxxzz_0_yzzzzz_0, \
                             g_0_xxxxxxzz_0_zzzzzz_0, \
                             g_0_xxxxxzz_0_xxxxxz_0,  \
                             g_0_xxxxxzz_0_xxxxxz_1,  \
                             g_0_xxxxxzz_0_xxxxyz_0,  \
                             g_0_xxxxxzz_0_xxxxyz_1,  \
                             g_0_xxxxxzz_0_xxxxz_1,   \
                             g_0_xxxxxzz_0_xxxxzz_0,  \
                             g_0_xxxxxzz_0_xxxxzz_1,  \
                             g_0_xxxxxzz_0_xxxyyz_0,  \
                             g_0_xxxxxzz_0_xxxyyz_1,  \
                             g_0_xxxxxzz_0_xxxyz_1,   \
                             g_0_xxxxxzz_0_xxxyzz_0,  \
                             g_0_xxxxxzz_0_xxxyzz_1,  \
                             g_0_xxxxxzz_0_xxxzz_1,   \
                             g_0_xxxxxzz_0_xxxzzz_0,  \
                             g_0_xxxxxzz_0_xxxzzz_1,  \
                             g_0_xxxxxzz_0_xxyyyz_0,  \
                             g_0_xxxxxzz_0_xxyyyz_1,  \
                             g_0_xxxxxzz_0_xxyyz_1,   \
                             g_0_xxxxxzz_0_xxyyzz_0,  \
                             g_0_xxxxxzz_0_xxyyzz_1,  \
                             g_0_xxxxxzz_0_xxyzz_1,   \
                             g_0_xxxxxzz_0_xxyzzz_0,  \
                             g_0_xxxxxzz_0_xxyzzz_1,  \
                             g_0_xxxxxzz_0_xxzzz_1,   \
                             g_0_xxxxxzz_0_xxzzzz_0,  \
                             g_0_xxxxxzz_0_xxzzzz_1,  \
                             g_0_xxxxxzz_0_xyyyyz_0,  \
                             g_0_xxxxxzz_0_xyyyyz_1,  \
                             g_0_xxxxxzz_0_xyyyz_1,   \
                             g_0_xxxxxzz_0_xyyyzz_0,  \
                             g_0_xxxxxzz_0_xyyyzz_1,  \
                             g_0_xxxxxzz_0_xyyzz_1,   \
                             g_0_xxxxxzz_0_xyyzzz_0,  \
                             g_0_xxxxxzz_0_xyyzzz_1,  \
                             g_0_xxxxxzz_0_xyzzz_1,   \
                             g_0_xxxxxzz_0_xyzzzz_0,  \
                             g_0_xxxxxzz_0_xyzzzz_1,  \
                             g_0_xxxxxzz_0_xzzzz_1,   \
                             g_0_xxxxxzz_0_xzzzzz_0,  \
                             g_0_xxxxxzz_0_xzzzzz_1,  \
                             g_0_xxxxxzz_0_yyyyyy_0,  \
                             g_0_xxxxxzz_0_yyyyyy_1,  \
                             g_0_xxxxxzz_0_yyyyyz_0,  \
                             g_0_xxxxxzz_0_yyyyyz_1,  \
                             g_0_xxxxxzz_0_yyyyz_1,   \
                             g_0_xxxxxzz_0_yyyyzz_0,  \
                             g_0_xxxxxzz_0_yyyyzz_1,  \
                             g_0_xxxxxzz_0_yyyzz_1,   \
                             g_0_xxxxxzz_0_yyyzzz_0,  \
                             g_0_xxxxxzz_0_yyyzzz_1,  \
                             g_0_xxxxxzz_0_yyzzz_1,   \
                             g_0_xxxxxzz_0_yyzzzz_0,  \
                             g_0_xxxxxzz_0_yyzzzz_1,  \
                             g_0_xxxxxzz_0_yzzzz_1,   \
                             g_0_xxxxxzz_0_yzzzzz_0,  \
                             g_0_xxxxxzz_0_yzzzzz_1,  \
                             g_0_xxxxxzz_0_zzzzz_1,   \
                             g_0_xxxxxzz_0_zzzzzz_0,  \
                             g_0_xxxxxzz_0_zzzzzz_1,  \
                             g_0_xxxxzz_0_xxxxxz_0,   \
                             g_0_xxxxzz_0_xxxxxz_1,   \
                             g_0_xxxxzz_0_xxxxyz_0,   \
                             g_0_xxxxzz_0_xxxxyz_1,   \
                             g_0_xxxxzz_0_xxxxzz_0,   \
                             g_0_xxxxzz_0_xxxxzz_1,   \
                             g_0_xxxxzz_0_xxxyyz_0,   \
                             g_0_xxxxzz_0_xxxyyz_1,   \
                             g_0_xxxxzz_0_xxxyzz_0,   \
                             g_0_xxxxzz_0_xxxyzz_1,   \
                             g_0_xxxxzz_0_xxxzzz_0,   \
                             g_0_xxxxzz_0_xxxzzz_1,   \
                             g_0_xxxxzz_0_xxyyyz_0,   \
                             g_0_xxxxzz_0_xxyyyz_1,   \
                             g_0_xxxxzz_0_xxyyzz_0,   \
                             g_0_xxxxzz_0_xxyyzz_1,   \
                             g_0_xxxxzz_0_xxyzzz_0,   \
                             g_0_xxxxzz_0_xxyzzz_1,   \
                             g_0_xxxxzz_0_xxzzzz_0,   \
                             g_0_xxxxzz_0_xxzzzz_1,   \
                             g_0_xxxxzz_0_xyyyyz_0,   \
                             g_0_xxxxzz_0_xyyyyz_1,   \
                             g_0_xxxxzz_0_xyyyzz_0,   \
                             g_0_xxxxzz_0_xyyyzz_1,   \
                             g_0_xxxxzz_0_xyyzzz_0,   \
                             g_0_xxxxzz_0_xyyzzz_1,   \
                             g_0_xxxxzz_0_xyzzzz_0,   \
                             g_0_xxxxzz_0_xyzzzz_1,   \
                             g_0_xxxxzz_0_xzzzzz_0,   \
                             g_0_xxxxzz_0_xzzzzz_1,   \
                             g_0_xxxxzz_0_yyyyyy_0,   \
                             g_0_xxxxzz_0_yyyyyy_1,   \
                             g_0_xxxxzz_0_yyyyyz_0,   \
                             g_0_xxxxzz_0_yyyyyz_1,   \
                             g_0_xxxxzz_0_yyyyzz_0,   \
                             g_0_xxxxzz_0_yyyyzz_1,   \
                             g_0_xxxxzz_0_yyyzzz_0,   \
                             g_0_xxxxzz_0_yyyzzz_1,   \
                             g_0_xxxxzz_0_yyzzzz_0,   \
                             g_0_xxxxzz_0_yyzzzz_1,   \
                             g_0_xxxxzz_0_yzzzzz_0,   \
                             g_0_xxxxzz_0_yzzzzz_1,   \
                             g_0_xxxxzz_0_zzzzzz_0,   \
                             g_0_xxxxzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_xxxxxx_0[i] = g_0_xxxxxx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxxx_0[i] * pb_z +
                                     g_0_xxxxxxz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxxy_0[i] = g_0_xxxxxx_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxxy_0[i] * pb_z +
                                     g_0_xxxxxxz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxxz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxyy_0[i] = g_0_xxxxxx_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxyy_0[i] * pb_z +
                                     g_0_xxxxxxz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxyz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxxzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxzz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxyyy_0[i] = g_0_xxxxxx_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxyyy_0[i] * pb_z +
                                     g_0_xxxxxxz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxyyz_0[i] = 5.0 * g_0_xxxxzz_0_xxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxyzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxzzz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyyyy_0[i] = g_0_xxxxxx_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxyyyy_0[i] * pb_z +
                                     g_0_xxxxxxz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyyzz_0[i] = 5.0 * g_0_xxxxzz_0_xxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxzzzz_0[i] * pb_x +
                                     g_0_xxxxxzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyyyy_0[i] = g_0_xxxxxx_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xyyyyy_0[i] * pb_z +
                                     g_0_xxxxxxz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyz_0[i] * pb_x + g_0_xxxxxzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyyzz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyzz_0[i] * pb_x + g_0_xxxxxzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyzzz_0[i] = 5.0 * g_0_xxxxzz_0_xyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyzzz_0[i] * pb_x + g_0_xxxxxzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyyy_0[i] = 5.0 * g_0_xxxxzz_0_yyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyyyyy_0[i] * pb_x + g_0_xxxxxzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyyz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyyyyz_0[i] * pb_x + g_0_xxxxxzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyyyzz_0[i] * pb_x + g_0_xxxxxzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyzzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyyzzz_0[i] * pb_x + g_0_xxxxxzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyzzzz_0[i] = 5.0 * g_0_xxxxzz_0_yyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yyzzzz_0[i] * pb_x + g_0_xxxxxzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_yzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_yzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_zzzzzz_0[i] = 5.0 * g_0_xxxxzz_0_zzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_zzzzzz_0[i] * pb_x + g_0_xxxxxzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-196 components of targeted buffer : SLSI

    auto g_0_xxxxxyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 168);

    auto g_0_xxxxxyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 169);

    auto g_0_xxxxxyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 170);

    auto g_0_xxxxxyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 171);

    auto g_0_xxxxxyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 172);

    auto g_0_xxxxxyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 173);

    auto g_0_xxxxxyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 174);

    auto g_0_xxxxxyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 175);

    auto g_0_xxxxxyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 176);

    auto g_0_xxxxxyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 177);

    auto g_0_xxxxxyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 178);

    auto g_0_xxxxxyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 179);

    auto g_0_xxxxxyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 180);

    auto g_0_xxxxxyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 181);

    auto g_0_xxxxxyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 182);

    auto g_0_xxxxxyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 183);

    auto g_0_xxxxxyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 184);

    auto g_0_xxxxxyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 185);

    auto g_0_xxxxxyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 186);

    auto g_0_xxxxxyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 187);

    auto g_0_xxxxxyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 188);

    auto g_0_xxxxxyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 189);

    auto g_0_xxxxxyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 190);

    auto g_0_xxxxxyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 191);

    auto g_0_xxxxxyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 192);

    auto g_0_xxxxxyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 193);

    auto g_0_xxxxxyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 194);

    auto g_0_xxxxxyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 195);

#pragma omp simd aligned(g_0_xxxxxy_0_xxxxxx_0,       \
                             g_0_xxxxxy_0_xxxxxx_1,   \
                             g_0_xxxxxy_0_xxxxxz_0,   \
                             g_0_xxxxxy_0_xxxxxz_1,   \
                             g_0_xxxxxy_0_xxxxzz_0,   \
                             g_0_xxxxxy_0_xxxxzz_1,   \
                             g_0_xxxxxy_0_xxxzzz_0,   \
                             g_0_xxxxxy_0_xxxzzz_1,   \
                             g_0_xxxxxy_0_xxzzzz_0,   \
                             g_0_xxxxxy_0_xxzzzz_1,   \
                             g_0_xxxxxy_0_xzzzzz_0,   \
                             g_0_xxxxxy_0_xzzzzz_1,   \
                             g_0_xxxxxyy_0_xxxxxx_0,  \
                             g_0_xxxxxyy_0_xxxxxx_1,  \
                             g_0_xxxxxyy_0_xxxxxz_0,  \
                             g_0_xxxxxyy_0_xxxxxz_1,  \
                             g_0_xxxxxyy_0_xxxxzz_0,  \
                             g_0_xxxxxyy_0_xxxxzz_1,  \
                             g_0_xxxxxyy_0_xxxzzz_0,  \
                             g_0_xxxxxyy_0_xxxzzz_1,  \
                             g_0_xxxxxyy_0_xxzzzz_0,  \
                             g_0_xxxxxyy_0_xxzzzz_1,  \
                             g_0_xxxxxyy_0_xzzzzz_0,  \
                             g_0_xxxxxyy_0_xzzzzz_1,  \
                             g_0_xxxxxyyy_0_xxxxxx_0, \
                             g_0_xxxxxyyy_0_xxxxxy_0, \
                             g_0_xxxxxyyy_0_xxxxxz_0, \
                             g_0_xxxxxyyy_0_xxxxyy_0, \
                             g_0_xxxxxyyy_0_xxxxyz_0, \
                             g_0_xxxxxyyy_0_xxxxzz_0, \
                             g_0_xxxxxyyy_0_xxxyyy_0, \
                             g_0_xxxxxyyy_0_xxxyyz_0, \
                             g_0_xxxxxyyy_0_xxxyzz_0, \
                             g_0_xxxxxyyy_0_xxxzzz_0, \
                             g_0_xxxxxyyy_0_xxyyyy_0, \
                             g_0_xxxxxyyy_0_xxyyyz_0, \
                             g_0_xxxxxyyy_0_xxyyzz_0, \
                             g_0_xxxxxyyy_0_xxyzzz_0, \
                             g_0_xxxxxyyy_0_xxzzzz_0, \
                             g_0_xxxxxyyy_0_xyyyyy_0, \
                             g_0_xxxxxyyy_0_xyyyyz_0, \
                             g_0_xxxxxyyy_0_xyyyzz_0, \
                             g_0_xxxxxyyy_0_xyyzzz_0, \
                             g_0_xxxxxyyy_0_xyzzzz_0, \
                             g_0_xxxxxyyy_0_xzzzzz_0, \
                             g_0_xxxxxyyy_0_yyyyyy_0, \
                             g_0_xxxxxyyy_0_yyyyyz_0, \
                             g_0_xxxxxyyy_0_yyyyzz_0, \
                             g_0_xxxxxyyy_0_yyyzzz_0, \
                             g_0_xxxxxyyy_0_yyzzzz_0, \
                             g_0_xxxxxyyy_0_yzzzzz_0, \
                             g_0_xxxxxyyy_0_zzzzzz_0, \
                             g_0_xxxxyyy_0_xxxxxy_0,  \
                             g_0_xxxxyyy_0_xxxxxy_1,  \
                             g_0_xxxxyyy_0_xxxxy_1,   \
                             g_0_xxxxyyy_0_xxxxyy_0,  \
                             g_0_xxxxyyy_0_xxxxyy_1,  \
                             g_0_xxxxyyy_0_xxxxyz_0,  \
                             g_0_xxxxyyy_0_xxxxyz_1,  \
                             g_0_xxxxyyy_0_xxxyy_1,   \
                             g_0_xxxxyyy_0_xxxyyy_0,  \
                             g_0_xxxxyyy_0_xxxyyy_1,  \
                             g_0_xxxxyyy_0_xxxyyz_0,  \
                             g_0_xxxxyyy_0_xxxyyz_1,  \
                             g_0_xxxxyyy_0_xxxyz_1,   \
                             g_0_xxxxyyy_0_xxxyzz_0,  \
                             g_0_xxxxyyy_0_xxxyzz_1,  \
                             g_0_xxxxyyy_0_xxyyy_1,   \
                             g_0_xxxxyyy_0_xxyyyy_0,  \
                             g_0_xxxxyyy_0_xxyyyy_1,  \
                             g_0_xxxxyyy_0_xxyyyz_0,  \
                             g_0_xxxxyyy_0_xxyyyz_1,  \
                             g_0_xxxxyyy_0_xxyyz_1,   \
                             g_0_xxxxyyy_0_xxyyzz_0,  \
                             g_0_xxxxyyy_0_xxyyzz_1,  \
                             g_0_xxxxyyy_0_xxyzz_1,   \
                             g_0_xxxxyyy_0_xxyzzz_0,  \
                             g_0_xxxxyyy_0_xxyzzz_1,  \
                             g_0_xxxxyyy_0_xyyyy_1,   \
                             g_0_xxxxyyy_0_xyyyyy_0,  \
                             g_0_xxxxyyy_0_xyyyyy_1,  \
                             g_0_xxxxyyy_0_xyyyyz_0,  \
                             g_0_xxxxyyy_0_xyyyyz_1,  \
                             g_0_xxxxyyy_0_xyyyz_1,   \
                             g_0_xxxxyyy_0_xyyyzz_0,  \
                             g_0_xxxxyyy_0_xyyyzz_1,  \
                             g_0_xxxxyyy_0_xyyzz_1,   \
                             g_0_xxxxyyy_0_xyyzzz_0,  \
                             g_0_xxxxyyy_0_xyyzzz_1,  \
                             g_0_xxxxyyy_0_xyzzz_1,   \
                             g_0_xxxxyyy_0_xyzzzz_0,  \
                             g_0_xxxxyyy_0_xyzzzz_1,  \
                             g_0_xxxxyyy_0_yyyyy_1,   \
                             g_0_xxxxyyy_0_yyyyyy_0,  \
                             g_0_xxxxyyy_0_yyyyyy_1,  \
                             g_0_xxxxyyy_0_yyyyyz_0,  \
                             g_0_xxxxyyy_0_yyyyyz_1,  \
                             g_0_xxxxyyy_0_yyyyz_1,   \
                             g_0_xxxxyyy_0_yyyyzz_0,  \
                             g_0_xxxxyyy_0_yyyyzz_1,  \
                             g_0_xxxxyyy_0_yyyzz_1,   \
                             g_0_xxxxyyy_0_yyyzzz_0,  \
                             g_0_xxxxyyy_0_yyyzzz_1,  \
                             g_0_xxxxyyy_0_yyzzz_1,   \
                             g_0_xxxxyyy_0_yyzzzz_0,  \
                             g_0_xxxxyyy_0_yyzzzz_1,  \
                             g_0_xxxxyyy_0_yzzzz_1,   \
                             g_0_xxxxyyy_0_yzzzzz_0,  \
                             g_0_xxxxyyy_0_yzzzzz_1,  \
                             g_0_xxxxyyy_0_zzzzzz_0,  \
                             g_0_xxxxyyy_0_zzzzzz_1,  \
                             g_0_xxxyyy_0_xxxxxy_0,   \
                             g_0_xxxyyy_0_xxxxxy_1,   \
                             g_0_xxxyyy_0_xxxxyy_0,   \
                             g_0_xxxyyy_0_xxxxyy_1,   \
                             g_0_xxxyyy_0_xxxxyz_0,   \
                             g_0_xxxyyy_0_xxxxyz_1,   \
                             g_0_xxxyyy_0_xxxyyy_0,   \
                             g_0_xxxyyy_0_xxxyyy_1,   \
                             g_0_xxxyyy_0_xxxyyz_0,   \
                             g_0_xxxyyy_0_xxxyyz_1,   \
                             g_0_xxxyyy_0_xxxyzz_0,   \
                             g_0_xxxyyy_0_xxxyzz_1,   \
                             g_0_xxxyyy_0_xxyyyy_0,   \
                             g_0_xxxyyy_0_xxyyyy_1,   \
                             g_0_xxxyyy_0_xxyyyz_0,   \
                             g_0_xxxyyy_0_xxyyyz_1,   \
                             g_0_xxxyyy_0_xxyyzz_0,   \
                             g_0_xxxyyy_0_xxyyzz_1,   \
                             g_0_xxxyyy_0_xxyzzz_0,   \
                             g_0_xxxyyy_0_xxyzzz_1,   \
                             g_0_xxxyyy_0_xyyyyy_0,   \
                             g_0_xxxyyy_0_xyyyyy_1,   \
                             g_0_xxxyyy_0_xyyyyz_0,   \
                             g_0_xxxyyy_0_xyyyyz_1,   \
                             g_0_xxxyyy_0_xyyyzz_0,   \
                             g_0_xxxyyy_0_xyyyzz_1,   \
                             g_0_xxxyyy_0_xyyzzz_0,   \
                             g_0_xxxyyy_0_xyyzzz_1,   \
                             g_0_xxxyyy_0_xyzzzz_0,   \
                             g_0_xxxyyy_0_xyzzzz_1,   \
                             g_0_xxxyyy_0_yyyyyy_0,   \
                             g_0_xxxyyy_0_yyyyyy_1,   \
                             g_0_xxxyyy_0_yyyyyz_0,   \
                             g_0_xxxyyy_0_yyyyyz_1,   \
                             g_0_xxxyyy_0_yyyyzz_0,   \
                             g_0_xxxyyy_0_yyyyzz_1,   \
                             g_0_xxxyyy_0_yyyzzz_0,   \
                             g_0_xxxyyy_0_yyyzzz_1,   \
                             g_0_xxxyyy_0_yyzzzz_0,   \
                             g_0_xxxyyy_0_yyzzzz_1,   \
                             g_0_xxxyyy_0_yzzzzz_0,   \
                             g_0_xxxyyy_0_yzzzzz_1,   \
                             g_0_xxxyyy_0_zzzzzz_0,   \
                             g_0_xxxyyy_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_xxxxxx_0[i] = 2.0 * g_0_xxxxxy_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_xxxxxx_0[i] * pb_y + g_0_xxxxxyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxxxy_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxy_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxxz_0[i] = 2.0 * g_0_xxxxxy_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_xxxxxz_0[i] * pb_y + g_0_xxxxxyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxxyy_0[i] = 4.0 * g_0_xxxyyy_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyy_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxyz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxzz_0[i] = 2.0 * g_0_xxxxxy_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_xxxxzz_0[i] * pb_y + g_0_xxxxxyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxyyy_0[i] = 4.0 * g_0_xxxyyy_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyy_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxyyz_0[i] = 4.0 * g_0_xxxyyy_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxyzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxzzz_0[i] = 2.0 * g_0_xxxxxy_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_xxxzzz_0[i] * pb_y + g_0_xxxxxyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxyyyy_0[i] = 4.0 * g_0_xxxyyy_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyy_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyyyz_0[i] = 4.0 * g_0_xxxyyy_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyyzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxxyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxzzzz_0[i] = 2.0 * g_0_xxxxxy_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_xxzzzz_0[i] * pb_y + g_0_xxxxxyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xyyyyy_0[i] = 4.0 * g_0_xxxyyy_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyy_0[i] * pb_x + g_0_xxxxyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyyyz_0[i] = 4.0 * g_0_xxxyyy_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyz_0[i] * pb_x + g_0_xxxxyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyyzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyzz_0[i] * pb_x + g_0_xxxxyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyzzz_0[i] * pb_x + g_0_xxxxyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzzzz_0[i] * pb_x + g_0_xxxxyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xzzzzz_0[i] = 2.0 * g_0_xxxxxy_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxxyy_0_xzzzzz_0[i] * pb_y + g_0_xxxxxyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_yyyyyy_0[i] = 4.0 * g_0_xxxyyy_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyyyy_0[i] * pb_x + g_0_xxxxyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyyyz_0[i] = 4.0 * g_0_xxxyyy_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyyyz_0[i] * pb_x + g_0_xxxxyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyyzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyyzz_0[i] * pb_x + g_0_xxxxyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyyzzz_0[i] * pb_x + g_0_xxxxyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yyzzzz_0[i] * pb_x + g_0_xxxxyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_yzzzzz_0[i] * pb_x + g_0_xxxxyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_zzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_zzzzzz_0[i] * pb_x + g_0_xxxxyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 196-224 components of targeted buffer : SLSI

    auto g_0_xxxxxyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 196);

    auto g_0_xxxxxyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 197);

    auto g_0_xxxxxyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 198);

    auto g_0_xxxxxyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 199);

    auto g_0_xxxxxyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 200);

    auto g_0_xxxxxyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 201);

    auto g_0_xxxxxyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 202);

    auto g_0_xxxxxyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 203);

    auto g_0_xxxxxyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 204);

    auto g_0_xxxxxyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 205);

    auto g_0_xxxxxyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 206);

    auto g_0_xxxxxyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 207);

    auto g_0_xxxxxyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 208);

    auto g_0_xxxxxyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 209);

    auto g_0_xxxxxyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 210);

    auto g_0_xxxxxyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 211);

    auto g_0_xxxxxyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 212);

    auto g_0_xxxxxyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 213);

    auto g_0_xxxxxyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 214);

    auto g_0_xxxxxyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 215);

    auto g_0_xxxxxyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 216);

    auto g_0_xxxxxyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 217);

    auto g_0_xxxxxyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 218);

    auto g_0_xxxxxyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 219);

    auto g_0_xxxxxyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 220);

    auto g_0_xxxxxyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 221);

    auto g_0_xxxxxyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 222);

    auto g_0_xxxxxyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 223);

#pragma omp simd aligned(g_0_xxxxxyy_0_xxxxx_1,       \
                             g_0_xxxxxyy_0_xxxxxx_0,  \
                             g_0_xxxxxyy_0_xxxxxx_1,  \
                             g_0_xxxxxyy_0_xxxxxy_0,  \
                             g_0_xxxxxyy_0_xxxxxy_1,  \
                             g_0_xxxxxyy_0_xxxxxz_0,  \
                             g_0_xxxxxyy_0_xxxxxz_1,  \
                             g_0_xxxxxyy_0_xxxxy_1,   \
                             g_0_xxxxxyy_0_xxxxyy_0,  \
                             g_0_xxxxxyy_0_xxxxyy_1,  \
                             g_0_xxxxxyy_0_xxxxyz_0,  \
                             g_0_xxxxxyy_0_xxxxyz_1,  \
                             g_0_xxxxxyy_0_xxxxz_1,   \
                             g_0_xxxxxyy_0_xxxxzz_0,  \
                             g_0_xxxxxyy_0_xxxxzz_1,  \
                             g_0_xxxxxyy_0_xxxyy_1,   \
                             g_0_xxxxxyy_0_xxxyyy_0,  \
                             g_0_xxxxxyy_0_xxxyyy_1,  \
                             g_0_xxxxxyy_0_xxxyyz_0,  \
                             g_0_xxxxxyy_0_xxxyyz_1,  \
                             g_0_xxxxxyy_0_xxxyz_1,   \
                             g_0_xxxxxyy_0_xxxyzz_0,  \
                             g_0_xxxxxyy_0_xxxyzz_1,  \
                             g_0_xxxxxyy_0_xxxzz_1,   \
                             g_0_xxxxxyy_0_xxxzzz_0,  \
                             g_0_xxxxxyy_0_xxxzzz_1,  \
                             g_0_xxxxxyy_0_xxyyy_1,   \
                             g_0_xxxxxyy_0_xxyyyy_0,  \
                             g_0_xxxxxyy_0_xxyyyy_1,  \
                             g_0_xxxxxyy_0_xxyyyz_0,  \
                             g_0_xxxxxyy_0_xxyyyz_1,  \
                             g_0_xxxxxyy_0_xxyyz_1,   \
                             g_0_xxxxxyy_0_xxyyzz_0,  \
                             g_0_xxxxxyy_0_xxyyzz_1,  \
                             g_0_xxxxxyy_0_xxyzz_1,   \
                             g_0_xxxxxyy_0_xxyzzz_0,  \
                             g_0_xxxxxyy_0_xxyzzz_1,  \
                             g_0_xxxxxyy_0_xxzzz_1,   \
                             g_0_xxxxxyy_0_xxzzzz_0,  \
                             g_0_xxxxxyy_0_xxzzzz_1,  \
                             g_0_xxxxxyy_0_xyyyy_1,   \
                             g_0_xxxxxyy_0_xyyyyy_0,  \
                             g_0_xxxxxyy_0_xyyyyy_1,  \
                             g_0_xxxxxyy_0_xyyyyz_0,  \
                             g_0_xxxxxyy_0_xyyyyz_1,  \
                             g_0_xxxxxyy_0_xyyyz_1,   \
                             g_0_xxxxxyy_0_xyyyzz_0,  \
                             g_0_xxxxxyy_0_xyyyzz_1,  \
                             g_0_xxxxxyy_0_xyyzz_1,   \
                             g_0_xxxxxyy_0_xyyzzz_0,  \
                             g_0_xxxxxyy_0_xyyzzz_1,  \
                             g_0_xxxxxyy_0_xyzzz_1,   \
                             g_0_xxxxxyy_0_xyzzzz_0,  \
                             g_0_xxxxxyy_0_xyzzzz_1,  \
                             g_0_xxxxxyy_0_xzzzz_1,   \
                             g_0_xxxxxyy_0_xzzzzz_0,  \
                             g_0_xxxxxyy_0_xzzzzz_1,  \
                             g_0_xxxxxyy_0_yyyyy_1,   \
                             g_0_xxxxxyy_0_yyyyyy_0,  \
                             g_0_xxxxxyy_0_yyyyyy_1,  \
                             g_0_xxxxxyy_0_yyyyyz_0,  \
                             g_0_xxxxxyy_0_yyyyyz_1,  \
                             g_0_xxxxxyy_0_yyyyz_1,   \
                             g_0_xxxxxyy_0_yyyyzz_0,  \
                             g_0_xxxxxyy_0_yyyyzz_1,  \
                             g_0_xxxxxyy_0_yyyzz_1,   \
                             g_0_xxxxxyy_0_yyyzzz_0,  \
                             g_0_xxxxxyy_0_yyyzzz_1,  \
                             g_0_xxxxxyy_0_yyzzz_1,   \
                             g_0_xxxxxyy_0_yyzzzz_0,  \
                             g_0_xxxxxyy_0_yyzzzz_1,  \
                             g_0_xxxxxyy_0_yzzzz_1,   \
                             g_0_xxxxxyy_0_yzzzzz_0,  \
                             g_0_xxxxxyy_0_yzzzzz_1,  \
                             g_0_xxxxxyy_0_zzzzz_1,   \
                             g_0_xxxxxyy_0_zzzzzz_0,  \
                             g_0_xxxxxyy_0_zzzzzz_1,  \
                             g_0_xxxxxyyz_0_xxxxxx_0, \
                             g_0_xxxxxyyz_0_xxxxxy_0, \
                             g_0_xxxxxyyz_0_xxxxxz_0, \
                             g_0_xxxxxyyz_0_xxxxyy_0, \
                             g_0_xxxxxyyz_0_xxxxyz_0, \
                             g_0_xxxxxyyz_0_xxxxzz_0, \
                             g_0_xxxxxyyz_0_xxxyyy_0, \
                             g_0_xxxxxyyz_0_xxxyyz_0, \
                             g_0_xxxxxyyz_0_xxxyzz_0, \
                             g_0_xxxxxyyz_0_xxxzzz_0, \
                             g_0_xxxxxyyz_0_xxyyyy_0, \
                             g_0_xxxxxyyz_0_xxyyyz_0, \
                             g_0_xxxxxyyz_0_xxyyzz_0, \
                             g_0_xxxxxyyz_0_xxyzzz_0, \
                             g_0_xxxxxyyz_0_xxzzzz_0, \
                             g_0_xxxxxyyz_0_xyyyyy_0, \
                             g_0_xxxxxyyz_0_xyyyyz_0, \
                             g_0_xxxxxyyz_0_xyyyzz_0, \
                             g_0_xxxxxyyz_0_xyyzzz_0, \
                             g_0_xxxxxyyz_0_xyzzzz_0, \
                             g_0_xxxxxyyz_0_xzzzzz_0, \
                             g_0_xxxxxyyz_0_yyyyyy_0, \
                             g_0_xxxxxyyz_0_yyyyyz_0, \
                             g_0_xxxxxyyz_0_yyyyzz_0, \
                             g_0_xxxxxyyz_0_yyyzzz_0, \
                             g_0_xxxxxyyz_0_yyzzzz_0, \
                             g_0_xxxxxyyz_0_yzzzzz_0, \
                             g_0_xxxxxyyz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_xxxxxx_0[i] = g_0_xxxxxyy_0_xxxxxx_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxy_0[i] = g_0_xxxxxyy_0_xxxxxy_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxxz_0[i] = g_0_xxxxxyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxxz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxyy_0[i] = g_0_xxxxxyy_0_xxxxyy_0[i] * pb_z + g_0_xxxxxyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxyz_0[i] = g_0_xxxxxyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxyz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxzz_0[i] =
            2.0 * g_0_xxxxxyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyyy_0[i] = g_0_xxxxxyy_0_xxxyyy_0[i] * pb_z + g_0_xxxxxyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyyz_0[i] = g_0_xxxxxyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyyz_0[i] * pb_z + g_0_xxxxxyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyzz_0[i] =
            2.0 * g_0_xxxxxyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxzzz_0[i] =
            3.0 * g_0_xxxxxyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyyy_0[i] = g_0_xxxxxyy_0_xxyyyy_0[i] * pb_z + g_0_xxxxxyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyyz_0[i] = g_0_xxxxxyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyyz_0[i] * pb_z + g_0_xxxxxyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxxyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyzz_0[i] * pb_z + g_0_xxxxxyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyzzz_0[i] =
            3.0 * g_0_xxxxxyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxzzzz_0[i] =
            4.0 * g_0_xxxxxyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyyy_0[i] = g_0_xxxxxyy_0_xyyyyy_0[i] * pb_z + g_0_xxxxxyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyyz_0[i] = g_0_xxxxxyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyyz_0[i] * pb_z + g_0_xxxxxyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyzz_0[i] =
            2.0 * g_0_xxxxxyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyzz_0[i] * pb_z + g_0_xxxxxyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyzzz_0[i] =
            3.0 * g_0_xxxxxyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyzzz_0[i] * pb_z + g_0_xxxxxyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyzzzz_0[i] =
            4.0 * g_0_xxxxxyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xzzzzz_0[i] =
            5.0 * g_0_xxxxxyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyyy_0[i] = g_0_xxxxxyy_0_yyyyyy_0[i] * pb_z + g_0_xxxxxyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyyz_0[i] = g_0_xxxxxyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyyyz_0[i] * pb_z + g_0_xxxxxyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyzz_0[i] =
            2.0 * g_0_xxxxxyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyyzz_0[i] * pb_z + g_0_xxxxxyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxxyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyzzz_0[i] * pb_z + g_0_xxxxxyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyzzzz_0[i] =
            4.0 * g_0_xxxxxyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyzzzz_0[i] * pb_z + g_0_xxxxxyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yzzzzz_0[i] =
            5.0 * g_0_xxxxxyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_zzzzzz_0[i] =
            6.0 * g_0_xxxxxyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_zzzzzz_0[i] * pb_z + g_0_xxxxxyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 224-252 components of targeted buffer : SLSI

    auto g_0_xxxxxyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 224);

    auto g_0_xxxxxyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 225);

    auto g_0_xxxxxyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 226);

    auto g_0_xxxxxyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 227);

    auto g_0_xxxxxyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 228);

    auto g_0_xxxxxyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 229);

    auto g_0_xxxxxyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 230);

    auto g_0_xxxxxyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 231);

    auto g_0_xxxxxyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 232);

    auto g_0_xxxxxyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 233);

    auto g_0_xxxxxyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 234);

    auto g_0_xxxxxyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 235);

    auto g_0_xxxxxyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 236);

    auto g_0_xxxxxyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 237);

    auto g_0_xxxxxyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 238);

    auto g_0_xxxxxyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 239);

    auto g_0_xxxxxyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 240);

    auto g_0_xxxxxyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 241);

    auto g_0_xxxxxyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 242);

    auto g_0_xxxxxyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 243);

    auto g_0_xxxxxyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 244);

    auto g_0_xxxxxyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 245);

    auto g_0_xxxxxyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 246);

    auto g_0_xxxxxyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 247);

    auto g_0_xxxxxyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 248);

    auto g_0_xxxxxyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 249);

    auto g_0_xxxxxyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 250);

    auto g_0_xxxxxyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 251);

#pragma omp simd aligned(g_0_xxxxxyzz_0_xxxxxx_0,     \
                             g_0_xxxxxyzz_0_xxxxxy_0, \
                             g_0_xxxxxyzz_0_xxxxxz_0, \
                             g_0_xxxxxyzz_0_xxxxyy_0, \
                             g_0_xxxxxyzz_0_xxxxyz_0, \
                             g_0_xxxxxyzz_0_xxxxzz_0, \
                             g_0_xxxxxyzz_0_xxxyyy_0, \
                             g_0_xxxxxyzz_0_xxxyyz_0, \
                             g_0_xxxxxyzz_0_xxxyzz_0, \
                             g_0_xxxxxyzz_0_xxxzzz_0, \
                             g_0_xxxxxyzz_0_xxyyyy_0, \
                             g_0_xxxxxyzz_0_xxyyyz_0, \
                             g_0_xxxxxyzz_0_xxyyzz_0, \
                             g_0_xxxxxyzz_0_xxyzzz_0, \
                             g_0_xxxxxyzz_0_xxzzzz_0, \
                             g_0_xxxxxyzz_0_xyyyyy_0, \
                             g_0_xxxxxyzz_0_xyyyyz_0, \
                             g_0_xxxxxyzz_0_xyyyzz_0, \
                             g_0_xxxxxyzz_0_xyyzzz_0, \
                             g_0_xxxxxyzz_0_xyzzzz_0, \
                             g_0_xxxxxyzz_0_xzzzzz_0, \
                             g_0_xxxxxyzz_0_yyyyyy_0, \
                             g_0_xxxxxyzz_0_yyyyyz_0, \
                             g_0_xxxxxyzz_0_yyyyzz_0, \
                             g_0_xxxxxyzz_0_yyyzzz_0, \
                             g_0_xxxxxyzz_0_yyzzzz_0, \
                             g_0_xxxxxyzz_0_yzzzzz_0, \
                             g_0_xxxxxyzz_0_zzzzzz_0, \
                             g_0_xxxxxzz_0_xxxxx_1,   \
                             g_0_xxxxxzz_0_xxxxxx_0,  \
                             g_0_xxxxxzz_0_xxxxxx_1,  \
                             g_0_xxxxxzz_0_xxxxxy_0,  \
                             g_0_xxxxxzz_0_xxxxxy_1,  \
                             g_0_xxxxxzz_0_xxxxxz_0,  \
                             g_0_xxxxxzz_0_xxxxxz_1,  \
                             g_0_xxxxxzz_0_xxxxy_1,   \
                             g_0_xxxxxzz_0_xxxxyy_0,  \
                             g_0_xxxxxzz_0_xxxxyy_1,  \
                             g_0_xxxxxzz_0_xxxxyz_0,  \
                             g_0_xxxxxzz_0_xxxxyz_1,  \
                             g_0_xxxxxzz_0_xxxxz_1,   \
                             g_0_xxxxxzz_0_xxxxzz_0,  \
                             g_0_xxxxxzz_0_xxxxzz_1,  \
                             g_0_xxxxxzz_0_xxxyy_1,   \
                             g_0_xxxxxzz_0_xxxyyy_0,  \
                             g_0_xxxxxzz_0_xxxyyy_1,  \
                             g_0_xxxxxzz_0_xxxyyz_0,  \
                             g_0_xxxxxzz_0_xxxyyz_1,  \
                             g_0_xxxxxzz_0_xxxyz_1,   \
                             g_0_xxxxxzz_0_xxxyzz_0,  \
                             g_0_xxxxxzz_0_xxxyzz_1,  \
                             g_0_xxxxxzz_0_xxxzz_1,   \
                             g_0_xxxxxzz_0_xxxzzz_0,  \
                             g_0_xxxxxzz_0_xxxzzz_1,  \
                             g_0_xxxxxzz_0_xxyyy_1,   \
                             g_0_xxxxxzz_0_xxyyyy_0,  \
                             g_0_xxxxxzz_0_xxyyyy_1,  \
                             g_0_xxxxxzz_0_xxyyyz_0,  \
                             g_0_xxxxxzz_0_xxyyyz_1,  \
                             g_0_xxxxxzz_0_xxyyz_1,   \
                             g_0_xxxxxzz_0_xxyyzz_0,  \
                             g_0_xxxxxzz_0_xxyyzz_1,  \
                             g_0_xxxxxzz_0_xxyzz_1,   \
                             g_0_xxxxxzz_0_xxyzzz_0,  \
                             g_0_xxxxxzz_0_xxyzzz_1,  \
                             g_0_xxxxxzz_0_xxzzz_1,   \
                             g_0_xxxxxzz_0_xxzzzz_0,  \
                             g_0_xxxxxzz_0_xxzzzz_1,  \
                             g_0_xxxxxzz_0_xyyyy_1,   \
                             g_0_xxxxxzz_0_xyyyyy_0,  \
                             g_0_xxxxxzz_0_xyyyyy_1,  \
                             g_0_xxxxxzz_0_xyyyyz_0,  \
                             g_0_xxxxxzz_0_xyyyyz_1,  \
                             g_0_xxxxxzz_0_xyyyz_1,   \
                             g_0_xxxxxzz_0_xyyyzz_0,  \
                             g_0_xxxxxzz_0_xyyyzz_1,  \
                             g_0_xxxxxzz_0_xyyzz_1,   \
                             g_0_xxxxxzz_0_xyyzzz_0,  \
                             g_0_xxxxxzz_0_xyyzzz_1,  \
                             g_0_xxxxxzz_0_xyzzz_1,   \
                             g_0_xxxxxzz_0_xyzzzz_0,  \
                             g_0_xxxxxzz_0_xyzzzz_1,  \
                             g_0_xxxxxzz_0_xzzzz_1,   \
                             g_0_xxxxxzz_0_xzzzzz_0,  \
                             g_0_xxxxxzz_0_xzzzzz_1,  \
                             g_0_xxxxxzz_0_yyyyy_1,   \
                             g_0_xxxxxzz_0_yyyyyy_0,  \
                             g_0_xxxxxzz_0_yyyyyy_1,  \
                             g_0_xxxxxzz_0_yyyyyz_0,  \
                             g_0_xxxxxzz_0_yyyyyz_1,  \
                             g_0_xxxxxzz_0_yyyyz_1,   \
                             g_0_xxxxxzz_0_yyyyzz_0,  \
                             g_0_xxxxxzz_0_yyyyzz_1,  \
                             g_0_xxxxxzz_0_yyyzz_1,   \
                             g_0_xxxxxzz_0_yyyzzz_0,  \
                             g_0_xxxxxzz_0_yyyzzz_1,  \
                             g_0_xxxxxzz_0_yyzzz_1,   \
                             g_0_xxxxxzz_0_yyzzzz_0,  \
                             g_0_xxxxxzz_0_yyzzzz_1,  \
                             g_0_xxxxxzz_0_yzzzz_1,   \
                             g_0_xxxxxzz_0_yzzzzz_0,  \
                             g_0_xxxxxzz_0_yzzzzz_1,  \
                             g_0_xxxxxzz_0_zzzzz_1,   \
                             g_0_xxxxxzz_0_zzzzzz_0,  \
                             g_0_xxxxxzz_0_zzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_xxxxxx_0[i] = g_0_xxxxxzz_0_xxxxxx_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxy_0[i] = g_0_xxxxxzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxxy_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxxz_0[i] = g_0_xxxxxzz_0_xxxxxz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxyy_0[i] =
            2.0 * g_0_xxxxxzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyy_0[i] * pb_y + g_0_xxxxxzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxyz_0[i] = g_0_xxxxxzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxyz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxzz_0[i] = g_0_xxxxxzz_0_xxxxzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyyy_0[i] =
            3.0 * g_0_xxxxxzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyy_0[i] * pb_y + g_0_xxxxxzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyyz_0[i] =
            2.0 * g_0_xxxxxzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyyz_0[i] * pb_y + g_0_xxxxxzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyzz_0[i] = g_0_xxxxxzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxzzz_0[i] = g_0_xxxxxzz_0_xxxzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyyy_0[i] =
            4.0 * g_0_xxxxxzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyy_0[i] * pb_y + g_0_xxxxxzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyyz_0[i] =
            3.0 * g_0_xxxxxzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyyz_0[i] * pb_y + g_0_xxxxxzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxxzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyzz_0[i] * pb_y + g_0_xxxxxzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyzzz_0[i] = g_0_xxxxxzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxzzzz_0[i] = g_0_xxxxxzz_0_xxzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyyy_0[i] =
            5.0 * g_0_xxxxxzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyy_0[i] * pb_y + g_0_xxxxxzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyyz_0[i] =
            4.0 * g_0_xxxxxzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyyz_0[i] * pb_y + g_0_xxxxxzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyzz_0[i] =
            3.0 * g_0_xxxxxzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyzz_0[i] * pb_y + g_0_xxxxxzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyzzz_0[i] =
            2.0 * g_0_xxxxxzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyzzz_0[i] * pb_y + g_0_xxxxxzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyzzzz_0[i] = g_0_xxxxxzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xzzzzz_0[i] = g_0_xxxxxzz_0_xzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyyy_0[i] =
            6.0 * g_0_xxxxxzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyyy_0[i] * pb_y + g_0_xxxxxzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyyz_0[i] =
            5.0 * g_0_xxxxxzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyyz_0[i] * pb_y + g_0_xxxxxzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyzz_0[i] =
            4.0 * g_0_xxxxxzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyzz_0[i] * pb_y + g_0_xxxxxzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxxzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyzzz_0[i] * pb_y + g_0_xxxxxzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyzzzz_0[i] =
            2.0 * g_0_xxxxxzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyzzzz_0[i] * pb_y + g_0_xxxxxzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yzzzzz_0[i] = g_0_xxxxxzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_zzzzzz_0[i] = g_0_xxxxxzz_0_zzzzzz_0[i] * pb_y + g_0_xxxxxzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 252-280 components of targeted buffer : SLSI

    auto g_0_xxxxxzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 252);

    auto g_0_xxxxxzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 253);

    auto g_0_xxxxxzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 254);

    auto g_0_xxxxxzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 255);

    auto g_0_xxxxxzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 256);

    auto g_0_xxxxxzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 257);

    auto g_0_xxxxxzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 258);

    auto g_0_xxxxxzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 259);

    auto g_0_xxxxxzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 260);

    auto g_0_xxxxxzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 261);

    auto g_0_xxxxxzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 262);

    auto g_0_xxxxxzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 263);

    auto g_0_xxxxxzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 264);

    auto g_0_xxxxxzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 265);

    auto g_0_xxxxxzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 266);

    auto g_0_xxxxxzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 267);

    auto g_0_xxxxxzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 268);

    auto g_0_xxxxxzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 269);

    auto g_0_xxxxxzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 270);

    auto g_0_xxxxxzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 271);

    auto g_0_xxxxxzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 272);

    auto g_0_xxxxxzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 273);

    auto g_0_xxxxxzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 274);

    auto g_0_xxxxxzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 275);

    auto g_0_xxxxxzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 276);

    auto g_0_xxxxxzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 277);

    auto g_0_xxxxxzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 278);

    auto g_0_xxxxxzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 279);

#pragma omp simd aligned(g_0_xxxxxz_0_xxxxxx_0,       \
                             g_0_xxxxxz_0_xxxxxx_1,   \
                             g_0_xxxxxz_0_xxxxxy_0,   \
                             g_0_xxxxxz_0_xxxxxy_1,   \
                             g_0_xxxxxz_0_xxxxyy_0,   \
                             g_0_xxxxxz_0_xxxxyy_1,   \
                             g_0_xxxxxz_0_xxxyyy_0,   \
                             g_0_xxxxxz_0_xxxyyy_1,   \
                             g_0_xxxxxz_0_xxyyyy_0,   \
                             g_0_xxxxxz_0_xxyyyy_1,   \
                             g_0_xxxxxz_0_xyyyyy_0,   \
                             g_0_xxxxxz_0_xyyyyy_1,   \
                             g_0_xxxxxzz_0_xxxxxx_0,  \
                             g_0_xxxxxzz_0_xxxxxx_1,  \
                             g_0_xxxxxzz_0_xxxxxy_0,  \
                             g_0_xxxxxzz_0_xxxxxy_1,  \
                             g_0_xxxxxzz_0_xxxxyy_0,  \
                             g_0_xxxxxzz_0_xxxxyy_1,  \
                             g_0_xxxxxzz_0_xxxyyy_0,  \
                             g_0_xxxxxzz_0_xxxyyy_1,  \
                             g_0_xxxxxzz_0_xxyyyy_0,  \
                             g_0_xxxxxzz_0_xxyyyy_1,  \
                             g_0_xxxxxzz_0_xyyyyy_0,  \
                             g_0_xxxxxzz_0_xyyyyy_1,  \
                             g_0_xxxxxzzz_0_xxxxxx_0, \
                             g_0_xxxxxzzz_0_xxxxxy_0, \
                             g_0_xxxxxzzz_0_xxxxxz_0, \
                             g_0_xxxxxzzz_0_xxxxyy_0, \
                             g_0_xxxxxzzz_0_xxxxyz_0, \
                             g_0_xxxxxzzz_0_xxxxzz_0, \
                             g_0_xxxxxzzz_0_xxxyyy_0, \
                             g_0_xxxxxzzz_0_xxxyyz_0, \
                             g_0_xxxxxzzz_0_xxxyzz_0, \
                             g_0_xxxxxzzz_0_xxxzzz_0, \
                             g_0_xxxxxzzz_0_xxyyyy_0, \
                             g_0_xxxxxzzz_0_xxyyyz_0, \
                             g_0_xxxxxzzz_0_xxyyzz_0, \
                             g_0_xxxxxzzz_0_xxyzzz_0, \
                             g_0_xxxxxzzz_0_xxzzzz_0, \
                             g_0_xxxxxzzz_0_xyyyyy_0, \
                             g_0_xxxxxzzz_0_xyyyyz_0, \
                             g_0_xxxxxzzz_0_xyyyzz_0, \
                             g_0_xxxxxzzz_0_xyyzzz_0, \
                             g_0_xxxxxzzz_0_xyzzzz_0, \
                             g_0_xxxxxzzz_0_xzzzzz_0, \
                             g_0_xxxxxzzz_0_yyyyyy_0, \
                             g_0_xxxxxzzz_0_yyyyyz_0, \
                             g_0_xxxxxzzz_0_yyyyzz_0, \
                             g_0_xxxxxzzz_0_yyyzzz_0, \
                             g_0_xxxxxzzz_0_yyzzzz_0, \
                             g_0_xxxxxzzz_0_yzzzzz_0, \
                             g_0_xxxxxzzz_0_zzzzzz_0, \
                             g_0_xxxxzzz_0_xxxxxz_0,  \
                             g_0_xxxxzzz_0_xxxxxz_1,  \
                             g_0_xxxxzzz_0_xxxxyz_0,  \
                             g_0_xxxxzzz_0_xxxxyz_1,  \
                             g_0_xxxxzzz_0_xxxxz_1,   \
                             g_0_xxxxzzz_0_xxxxzz_0,  \
                             g_0_xxxxzzz_0_xxxxzz_1,  \
                             g_0_xxxxzzz_0_xxxyyz_0,  \
                             g_0_xxxxzzz_0_xxxyyz_1,  \
                             g_0_xxxxzzz_0_xxxyz_1,   \
                             g_0_xxxxzzz_0_xxxyzz_0,  \
                             g_0_xxxxzzz_0_xxxyzz_1,  \
                             g_0_xxxxzzz_0_xxxzz_1,   \
                             g_0_xxxxzzz_0_xxxzzz_0,  \
                             g_0_xxxxzzz_0_xxxzzz_1,  \
                             g_0_xxxxzzz_0_xxyyyz_0,  \
                             g_0_xxxxzzz_0_xxyyyz_1,  \
                             g_0_xxxxzzz_0_xxyyz_1,   \
                             g_0_xxxxzzz_0_xxyyzz_0,  \
                             g_0_xxxxzzz_0_xxyyzz_1,  \
                             g_0_xxxxzzz_0_xxyzz_1,   \
                             g_0_xxxxzzz_0_xxyzzz_0,  \
                             g_0_xxxxzzz_0_xxyzzz_1,  \
                             g_0_xxxxzzz_0_xxzzz_1,   \
                             g_0_xxxxzzz_0_xxzzzz_0,  \
                             g_0_xxxxzzz_0_xxzzzz_1,  \
                             g_0_xxxxzzz_0_xyyyyz_0,  \
                             g_0_xxxxzzz_0_xyyyyz_1,  \
                             g_0_xxxxzzz_0_xyyyz_1,   \
                             g_0_xxxxzzz_0_xyyyzz_0,  \
                             g_0_xxxxzzz_0_xyyyzz_1,  \
                             g_0_xxxxzzz_0_xyyzz_1,   \
                             g_0_xxxxzzz_0_xyyzzz_0,  \
                             g_0_xxxxzzz_0_xyyzzz_1,  \
                             g_0_xxxxzzz_0_xyzzz_1,   \
                             g_0_xxxxzzz_0_xyzzzz_0,  \
                             g_0_xxxxzzz_0_xyzzzz_1,  \
                             g_0_xxxxzzz_0_xzzzz_1,   \
                             g_0_xxxxzzz_0_xzzzzz_0,  \
                             g_0_xxxxzzz_0_xzzzzz_1,  \
                             g_0_xxxxzzz_0_yyyyyy_0,  \
                             g_0_xxxxzzz_0_yyyyyy_1,  \
                             g_0_xxxxzzz_0_yyyyyz_0,  \
                             g_0_xxxxzzz_0_yyyyyz_1,  \
                             g_0_xxxxzzz_0_yyyyz_1,   \
                             g_0_xxxxzzz_0_yyyyzz_0,  \
                             g_0_xxxxzzz_0_yyyyzz_1,  \
                             g_0_xxxxzzz_0_yyyzz_1,   \
                             g_0_xxxxzzz_0_yyyzzz_0,  \
                             g_0_xxxxzzz_0_yyyzzz_1,  \
                             g_0_xxxxzzz_0_yyzzz_1,   \
                             g_0_xxxxzzz_0_yyzzzz_0,  \
                             g_0_xxxxzzz_0_yyzzzz_1,  \
                             g_0_xxxxzzz_0_yzzzz_1,   \
                             g_0_xxxxzzz_0_yzzzzz_0,  \
                             g_0_xxxxzzz_0_yzzzzz_1,  \
                             g_0_xxxxzzz_0_zzzzz_1,   \
                             g_0_xxxxzzz_0_zzzzzz_0,  \
                             g_0_xxxxzzz_0_zzzzzz_1,  \
                             g_0_xxxzzz_0_xxxxxz_0,   \
                             g_0_xxxzzz_0_xxxxxz_1,   \
                             g_0_xxxzzz_0_xxxxyz_0,   \
                             g_0_xxxzzz_0_xxxxyz_1,   \
                             g_0_xxxzzz_0_xxxxzz_0,   \
                             g_0_xxxzzz_0_xxxxzz_1,   \
                             g_0_xxxzzz_0_xxxyyz_0,   \
                             g_0_xxxzzz_0_xxxyyz_1,   \
                             g_0_xxxzzz_0_xxxyzz_0,   \
                             g_0_xxxzzz_0_xxxyzz_1,   \
                             g_0_xxxzzz_0_xxxzzz_0,   \
                             g_0_xxxzzz_0_xxxzzz_1,   \
                             g_0_xxxzzz_0_xxyyyz_0,   \
                             g_0_xxxzzz_0_xxyyyz_1,   \
                             g_0_xxxzzz_0_xxyyzz_0,   \
                             g_0_xxxzzz_0_xxyyzz_1,   \
                             g_0_xxxzzz_0_xxyzzz_0,   \
                             g_0_xxxzzz_0_xxyzzz_1,   \
                             g_0_xxxzzz_0_xxzzzz_0,   \
                             g_0_xxxzzz_0_xxzzzz_1,   \
                             g_0_xxxzzz_0_xyyyyz_0,   \
                             g_0_xxxzzz_0_xyyyyz_1,   \
                             g_0_xxxzzz_0_xyyyzz_0,   \
                             g_0_xxxzzz_0_xyyyzz_1,   \
                             g_0_xxxzzz_0_xyyzzz_0,   \
                             g_0_xxxzzz_0_xyyzzz_1,   \
                             g_0_xxxzzz_0_xyzzzz_0,   \
                             g_0_xxxzzz_0_xyzzzz_1,   \
                             g_0_xxxzzz_0_xzzzzz_0,   \
                             g_0_xxxzzz_0_xzzzzz_1,   \
                             g_0_xxxzzz_0_yyyyyy_0,   \
                             g_0_xxxzzz_0_yyyyyy_1,   \
                             g_0_xxxzzz_0_yyyyyz_0,   \
                             g_0_xxxzzz_0_yyyyyz_1,   \
                             g_0_xxxzzz_0_yyyyzz_0,   \
                             g_0_xxxzzz_0_yyyyzz_1,   \
                             g_0_xxxzzz_0_yyyzzz_0,   \
                             g_0_xxxzzz_0_yyyzzz_1,   \
                             g_0_xxxzzz_0_yyzzzz_0,   \
                             g_0_xxxzzz_0_yyzzzz_1,   \
                             g_0_xxxzzz_0_yzzzzz_0,   \
                             g_0_xxxzzz_0_yzzzzz_1,   \
                             g_0_xxxzzz_0_zzzzzz_0,   \
                             g_0_xxxzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_xxxxxx_0[i] = 2.0 * g_0_xxxxxz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_xxxxxx_0[i] * pb_z + g_0_xxxxxzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxxxxz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_xxxxxy_0[i] * pb_z + g_0_xxxxxzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxxz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxxzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxxxz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_xxxxyy_0[i] * pb_z + g_0_xxxxxzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxyz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxxzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxxzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxzz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxxxxz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_xxxyyy_0[i] * pb_z + g_0_xxxxxzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxyzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxxzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxzzz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxxxxz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_xxyyyy_0[i] * pb_z + g_0_xxxxxzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxxzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxzzzz_0[i] * pb_x +
                                     g_0_xxxxzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxxxxz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxxzz_0_xyyyyy_0[i] * pb_z + g_0_xxxxxzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyz_0[i] * pb_x + g_0_xxxxzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyzz_0[i] * pb_x + g_0_xxxxzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyzzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyzzz_0[i] * pb_x + g_0_xxxxzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyyy_0[i] = 4.0 * g_0_xxxzzz_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyyyyy_0[i] * pb_x + g_0_xxxxzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyyz_0[i] = 4.0 * g_0_xxxzzz_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyyyyz_0[i] * pb_x + g_0_xxxxzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyyyzz_0[i] * pb_x + g_0_xxxxzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyyzzz_0[i] * pb_x + g_0_xxxxzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyzzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yyzzzz_0[i] * pb_x + g_0_xxxxzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_yzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_zzzzzz_0[i] = 4.0 * g_0_xxxzzz_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_zzzzzz_0[i] * pb_x + g_0_xxxxzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 280-308 components of targeted buffer : SLSI

    auto g_0_xxxxyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 280);

    auto g_0_xxxxyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 281);

    auto g_0_xxxxyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 282);

    auto g_0_xxxxyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 283);

    auto g_0_xxxxyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 284);

    auto g_0_xxxxyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 285);

    auto g_0_xxxxyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 286);

    auto g_0_xxxxyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 287);

    auto g_0_xxxxyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 288);

    auto g_0_xxxxyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 289);

    auto g_0_xxxxyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 290);

    auto g_0_xxxxyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 291);

    auto g_0_xxxxyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 292);

    auto g_0_xxxxyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 293);

    auto g_0_xxxxyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 294);

    auto g_0_xxxxyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 295);

    auto g_0_xxxxyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 296);

    auto g_0_xxxxyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 297);

    auto g_0_xxxxyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 298);

    auto g_0_xxxxyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 299);

    auto g_0_xxxxyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 300);

    auto g_0_xxxxyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 301);

    auto g_0_xxxxyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 302);

    auto g_0_xxxxyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 303);

    auto g_0_xxxxyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 304);

    auto g_0_xxxxyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 305);

    auto g_0_xxxxyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 306);

    auto g_0_xxxxyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 307);

#pragma omp simd aligned(g_0_xxxxyy_0_xxxxxx_0,       \
                             g_0_xxxxyy_0_xxxxxx_1,   \
                             g_0_xxxxyy_0_xxxxxz_0,   \
                             g_0_xxxxyy_0_xxxxxz_1,   \
                             g_0_xxxxyy_0_xxxxzz_0,   \
                             g_0_xxxxyy_0_xxxxzz_1,   \
                             g_0_xxxxyy_0_xxxzzz_0,   \
                             g_0_xxxxyy_0_xxxzzz_1,   \
                             g_0_xxxxyy_0_xxzzzz_0,   \
                             g_0_xxxxyy_0_xxzzzz_1,   \
                             g_0_xxxxyy_0_xzzzzz_0,   \
                             g_0_xxxxyy_0_xzzzzz_1,   \
                             g_0_xxxxyyy_0_xxxxxx_0,  \
                             g_0_xxxxyyy_0_xxxxxx_1,  \
                             g_0_xxxxyyy_0_xxxxxz_0,  \
                             g_0_xxxxyyy_0_xxxxxz_1,  \
                             g_0_xxxxyyy_0_xxxxzz_0,  \
                             g_0_xxxxyyy_0_xxxxzz_1,  \
                             g_0_xxxxyyy_0_xxxzzz_0,  \
                             g_0_xxxxyyy_0_xxxzzz_1,  \
                             g_0_xxxxyyy_0_xxzzzz_0,  \
                             g_0_xxxxyyy_0_xxzzzz_1,  \
                             g_0_xxxxyyy_0_xzzzzz_0,  \
                             g_0_xxxxyyy_0_xzzzzz_1,  \
                             g_0_xxxxyyyy_0_xxxxxx_0, \
                             g_0_xxxxyyyy_0_xxxxxy_0, \
                             g_0_xxxxyyyy_0_xxxxxz_0, \
                             g_0_xxxxyyyy_0_xxxxyy_0, \
                             g_0_xxxxyyyy_0_xxxxyz_0, \
                             g_0_xxxxyyyy_0_xxxxzz_0, \
                             g_0_xxxxyyyy_0_xxxyyy_0, \
                             g_0_xxxxyyyy_0_xxxyyz_0, \
                             g_0_xxxxyyyy_0_xxxyzz_0, \
                             g_0_xxxxyyyy_0_xxxzzz_0, \
                             g_0_xxxxyyyy_0_xxyyyy_0, \
                             g_0_xxxxyyyy_0_xxyyyz_0, \
                             g_0_xxxxyyyy_0_xxyyzz_0, \
                             g_0_xxxxyyyy_0_xxyzzz_0, \
                             g_0_xxxxyyyy_0_xxzzzz_0, \
                             g_0_xxxxyyyy_0_xyyyyy_0, \
                             g_0_xxxxyyyy_0_xyyyyz_0, \
                             g_0_xxxxyyyy_0_xyyyzz_0, \
                             g_0_xxxxyyyy_0_xyyzzz_0, \
                             g_0_xxxxyyyy_0_xyzzzz_0, \
                             g_0_xxxxyyyy_0_xzzzzz_0, \
                             g_0_xxxxyyyy_0_yyyyyy_0, \
                             g_0_xxxxyyyy_0_yyyyyz_0, \
                             g_0_xxxxyyyy_0_yyyyzz_0, \
                             g_0_xxxxyyyy_0_yyyzzz_0, \
                             g_0_xxxxyyyy_0_yyzzzz_0, \
                             g_0_xxxxyyyy_0_yzzzzz_0, \
                             g_0_xxxxyyyy_0_zzzzzz_0, \
                             g_0_xxxyyyy_0_xxxxxy_0,  \
                             g_0_xxxyyyy_0_xxxxxy_1,  \
                             g_0_xxxyyyy_0_xxxxy_1,   \
                             g_0_xxxyyyy_0_xxxxyy_0,  \
                             g_0_xxxyyyy_0_xxxxyy_1,  \
                             g_0_xxxyyyy_0_xxxxyz_0,  \
                             g_0_xxxyyyy_0_xxxxyz_1,  \
                             g_0_xxxyyyy_0_xxxyy_1,   \
                             g_0_xxxyyyy_0_xxxyyy_0,  \
                             g_0_xxxyyyy_0_xxxyyy_1,  \
                             g_0_xxxyyyy_0_xxxyyz_0,  \
                             g_0_xxxyyyy_0_xxxyyz_1,  \
                             g_0_xxxyyyy_0_xxxyz_1,   \
                             g_0_xxxyyyy_0_xxxyzz_0,  \
                             g_0_xxxyyyy_0_xxxyzz_1,  \
                             g_0_xxxyyyy_0_xxyyy_1,   \
                             g_0_xxxyyyy_0_xxyyyy_0,  \
                             g_0_xxxyyyy_0_xxyyyy_1,  \
                             g_0_xxxyyyy_0_xxyyyz_0,  \
                             g_0_xxxyyyy_0_xxyyyz_1,  \
                             g_0_xxxyyyy_0_xxyyz_1,   \
                             g_0_xxxyyyy_0_xxyyzz_0,  \
                             g_0_xxxyyyy_0_xxyyzz_1,  \
                             g_0_xxxyyyy_0_xxyzz_1,   \
                             g_0_xxxyyyy_0_xxyzzz_0,  \
                             g_0_xxxyyyy_0_xxyzzz_1,  \
                             g_0_xxxyyyy_0_xyyyy_1,   \
                             g_0_xxxyyyy_0_xyyyyy_0,  \
                             g_0_xxxyyyy_0_xyyyyy_1,  \
                             g_0_xxxyyyy_0_xyyyyz_0,  \
                             g_0_xxxyyyy_0_xyyyyz_1,  \
                             g_0_xxxyyyy_0_xyyyz_1,   \
                             g_0_xxxyyyy_0_xyyyzz_0,  \
                             g_0_xxxyyyy_0_xyyyzz_1,  \
                             g_0_xxxyyyy_0_xyyzz_1,   \
                             g_0_xxxyyyy_0_xyyzzz_0,  \
                             g_0_xxxyyyy_0_xyyzzz_1,  \
                             g_0_xxxyyyy_0_xyzzz_1,   \
                             g_0_xxxyyyy_0_xyzzzz_0,  \
                             g_0_xxxyyyy_0_xyzzzz_1,  \
                             g_0_xxxyyyy_0_yyyyy_1,   \
                             g_0_xxxyyyy_0_yyyyyy_0,  \
                             g_0_xxxyyyy_0_yyyyyy_1,  \
                             g_0_xxxyyyy_0_yyyyyz_0,  \
                             g_0_xxxyyyy_0_yyyyyz_1,  \
                             g_0_xxxyyyy_0_yyyyz_1,   \
                             g_0_xxxyyyy_0_yyyyzz_0,  \
                             g_0_xxxyyyy_0_yyyyzz_1,  \
                             g_0_xxxyyyy_0_yyyzz_1,   \
                             g_0_xxxyyyy_0_yyyzzz_0,  \
                             g_0_xxxyyyy_0_yyyzzz_1,  \
                             g_0_xxxyyyy_0_yyzzz_1,   \
                             g_0_xxxyyyy_0_yyzzzz_0,  \
                             g_0_xxxyyyy_0_yyzzzz_1,  \
                             g_0_xxxyyyy_0_yzzzz_1,   \
                             g_0_xxxyyyy_0_yzzzzz_0,  \
                             g_0_xxxyyyy_0_yzzzzz_1,  \
                             g_0_xxxyyyy_0_zzzzzz_0,  \
                             g_0_xxxyyyy_0_zzzzzz_1,  \
                             g_0_xxyyyy_0_xxxxxy_0,   \
                             g_0_xxyyyy_0_xxxxxy_1,   \
                             g_0_xxyyyy_0_xxxxyy_0,   \
                             g_0_xxyyyy_0_xxxxyy_1,   \
                             g_0_xxyyyy_0_xxxxyz_0,   \
                             g_0_xxyyyy_0_xxxxyz_1,   \
                             g_0_xxyyyy_0_xxxyyy_0,   \
                             g_0_xxyyyy_0_xxxyyy_1,   \
                             g_0_xxyyyy_0_xxxyyz_0,   \
                             g_0_xxyyyy_0_xxxyyz_1,   \
                             g_0_xxyyyy_0_xxxyzz_0,   \
                             g_0_xxyyyy_0_xxxyzz_1,   \
                             g_0_xxyyyy_0_xxyyyy_0,   \
                             g_0_xxyyyy_0_xxyyyy_1,   \
                             g_0_xxyyyy_0_xxyyyz_0,   \
                             g_0_xxyyyy_0_xxyyyz_1,   \
                             g_0_xxyyyy_0_xxyyzz_0,   \
                             g_0_xxyyyy_0_xxyyzz_1,   \
                             g_0_xxyyyy_0_xxyzzz_0,   \
                             g_0_xxyyyy_0_xxyzzz_1,   \
                             g_0_xxyyyy_0_xyyyyy_0,   \
                             g_0_xxyyyy_0_xyyyyy_1,   \
                             g_0_xxyyyy_0_xyyyyz_0,   \
                             g_0_xxyyyy_0_xyyyyz_1,   \
                             g_0_xxyyyy_0_xyyyzz_0,   \
                             g_0_xxyyyy_0_xyyyzz_1,   \
                             g_0_xxyyyy_0_xyyzzz_0,   \
                             g_0_xxyyyy_0_xyyzzz_1,   \
                             g_0_xxyyyy_0_xyzzzz_0,   \
                             g_0_xxyyyy_0_xyzzzz_1,   \
                             g_0_xxyyyy_0_yyyyyy_0,   \
                             g_0_xxyyyy_0_yyyyyy_1,   \
                             g_0_xxyyyy_0_yyyyyz_0,   \
                             g_0_xxyyyy_0_yyyyyz_1,   \
                             g_0_xxyyyy_0_yyyyzz_0,   \
                             g_0_xxyyyy_0_yyyyzz_1,   \
                             g_0_xxyyyy_0_yyyzzz_0,   \
                             g_0_xxyyyy_0_yyyzzz_1,   \
                             g_0_xxyyyy_0_yyzzzz_0,   \
                             g_0_xxyyyy_0_yyzzzz_1,   \
                             g_0_xxyyyy_0_yzzzzz_0,   \
                             g_0_xxyyyy_0_yzzzzz_1,   \
                             g_0_xxyyyy_0_zzzzzz_0,   \
                             g_0_xxyyyy_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_xxxxxx_0[i] = 3.0 * g_0_xxxxyy_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_xxxxxx_0[i] * pb_y + g_0_xxxxyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxxxy_0[i] = 3.0 * g_0_xxyyyy_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxy_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxxz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_xxxxxz_0[i] * pb_y + g_0_xxxxyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxxyy_0[i] = 3.0 * g_0_xxyyyy_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyy_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxyz_0[i] = 3.0 * g_0_xxyyyy_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_xxxxzz_0[i] * pb_y + g_0_xxxxyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxyyy_0[i] = 3.0 * g_0_xxyyyy_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyy_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxyyz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxyzz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_xxxzzz_0[i] * pb_y + g_0_xxxxyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxyyyy_0[i] = 3.0 * g_0_xxyyyy_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyy_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyyyz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyyzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxzzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_xxzzzz_0[i] * pb_y + g_0_xxxxyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xyyyyy_0[i] = 3.0 * g_0_xxyyyy_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyy_0[i] * pb_x + g_0_xxxyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyyyz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyz_0[i] * pb_x + g_0_xxxyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyyzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyzz_0[i] * pb_x + g_0_xxxyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyzzz_0[i] * pb_x + g_0_xxxyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyzzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzzzz_0[i] * pb_x + g_0_xxxyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xzzzzz_0[i] = 3.0 * g_0_xxxxyy_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxxyyy_0_xzzzzz_0[i] * pb_y + g_0_xxxxyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_yyyyyy_0[i] = 3.0 * g_0_xxyyyy_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyyyy_0[i] * pb_x + g_0_xxxyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyyyz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyyyz_0[i] * pb_x + g_0_xxxyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyyzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyyzz_0[i] * pb_x + g_0_xxxyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyyzzz_0[i] * pb_x + g_0_xxxyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyzzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yyzzzz_0[i] * pb_x + g_0_xxxyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yzzzzz_0[i] = 3.0 * g_0_xxyyyy_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_yzzzzz_0[i] * pb_x + g_0_xxxyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_zzzzzz_0[i] = 3.0 * g_0_xxyyyy_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_zzzzzz_0[i] * pb_x + g_0_xxxyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 308-336 components of targeted buffer : SLSI

    auto g_0_xxxxyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 308);

    auto g_0_xxxxyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 309);

    auto g_0_xxxxyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 310);

    auto g_0_xxxxyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 311);

    auto g_0_xxxxyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 312);

    auto g_0_xxxxyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 313);

    auto g_0_xxxxyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 314);

    auto g_0_xxxxyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 315);

    auto g_0_xxxxyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 316);

    auto g_0_xxxxyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 317);

    auto g_0_xxxxyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 318);

    auto g_0_xxxxyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 319);

    auto g_0_xxxxyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 320);

    auto g_0_xxxxyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 321);

    auto g_0_xxxxyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 322);

    auto g_0_xxxxyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 323);

    auto g_0_xxxxyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 324);

    auto g_0_xxxxyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 325);

    auto g_0_xxxxyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 326);

    auto g_0_xxxxyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 327);

    auto g_0_xxxxyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 328);

    auto g_0_xxxxyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 329);

    auto g_0_xxxxyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 330);

    auto g_0_xxxxyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 331);

    auto g_0_xxxxyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 332);

    auto g_0_xxxxyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 333);

    auto g_0_xxxxyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 334);

    auto g_0_xxxxyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 335);

#pragma omp simd aligned(g_0_xxxxyyy_0_xxxxx_1,       \
                             g_0_xxxxyyy_0_xxxxxx_0,  \
                             g_0_xxxxyyy_0_xxxxxx_1,  \
                             g_0_xxxxyyy_0_xxxxxy_0,  \
                             g_0_xxxxyyy_0_xxxxxy_1,  \
                             g_0_xxxxyyy_0_xxxxxz_0,  \
                             g_0_xxxxyyy_0_xxxxxz_1,  \
                             g_0_xxxxyyy_0_xxxxy_1,   \
                             g_0_xxxxyyy_0_xxxxyy_0,  \
                             g_0_xxxxyyy_0_xxxxyy_1,  \
                             g_0_xxxxyyy_0_xxxxyz_0,  \
                             g_0_xxxxyyy_0_xxxxyz_1,  \
                             g_0_xxxxyyy_0_xxxxz_1,   \
                             g_0_xxxxyyy_0_xxxxzz_0,  \
                             g_0_xxxxyyy_0_xxxxzz_1,  \
                             g_0_xxxxyyy_0_xxxyy_1,   \
                             g_0_xxxxyyy_0_xxxyyy_0,  \
                             g_0_xxxxyyy_0_xxxyyy_1,  \
                             g_0_xxxxyyy_0_xxxyyz_0,  \
                             g_0_xxxxyyy_0_xxxyyz_1,  \
                             g_0_xxxxyyy_0_xxxyz_1,   \
                             g_0_xxxxyyy_0_xxxyzz_0,  \
                             g_0_xxxxyyy_0_xxxyzz_1,  \
                             g_0_xxxxyyy_0_xxxzz_1,   \
                             g_0_xxxxyyy_0_xxxzzz_0,  \
                             g_0_xxxxyyy_0_xxxzzz_1,  \
                             g_0_xxxxyyy_0_xxyyy_1,   \
                             g_0_xxxxyyy_0_xxyyyy_0,  \
                             g_0_xxxxyyy_0_xxyyyy_1,  \
                             g_0_xxxxyyy_0_xxyyyz_0,  \
                             g_0_xxxxyyy_0_xxyyyz_1,  \
                             g_0_xxxxyyy_0_xxyyz_1,   \
                             g_0_xxxxyyy_0_xxyyzz_0,  \
                             g_0_xxxxyyy_0_xxyyzz_1,  \
                             g_0_xxxxyyy_0_xxyzz_1,   \
                             g_0_xxxxyyy_0_xxyzzz_0,  \
                             g_0_xxxxyyy_0_xxyzzz_1,  \
                             g_0_xxxxyyy_0_xxzzz_1,   \
                             g_0_xxxxyyy_0_xxzzzz_0,  \
                             g_0_xxxxyyy_0_xxzzzz_1,  \
                             g_0_xxxxyyy_0_xyyyy_1,   \
                             g_0_xxxxyyy_0_xyyyyy_0,  \
                             g_0_xxxxyyy_0_xyyyyy_1,  \
                             g_0_xxxxyyy_0_xyyyyz_0,  \
                             g_0_xxxxyyy_0_xyyyyz_1,  \
                             g_0_xxxxyyy_0_xyyyz_1,   \
                             g_0_xxxxyyy_0_xyyyzz_0,  \
                             g_0_xxxxyyy_0_xyyyzz_1,  \
                             g_0_xxxxyyy_0_xyyzz_1,   \
                             g_0_xxxxyyy_0_xyyzzz_0,  \
                             g_0_xxxxyyy_0_xyyzzz_1,  \
                             g_0_xxxxyyy_0_xyzzz_1,   \
                             g_0_xxxxyyy_0_xyzzzz_0,  \
                             g_0_xxxxyyy_0_xyzzzz_1,  \
                             g_0_xxxxyyy_0_xzzzz_1,   \
                             g_0_xxxxyyy_0_xzzzzz_0,  \
                             g_0_xxxxyyy_0_xzzzzz_1,  \
                             g_0_xxxxyyy_0_yyyyy_1,   \
                             g_0_xxxxyyy_0_yyyyyy_0,  \
                             g_0_xxxxyyy_0_yyyyyy_1,  \
                             g_0_xxxxyyy_0_yyyyyz_0,  \
                             g_0_xxxxyyy_0_yyyyyz_1,  \
                             g_0_xxxxyyy_0_yyyyz_1,   \
                             g_0_xxxxyyy_0_yyyyzz_0,  \
                             g_0_xxxxyyy_0_yyyyzz_1,  \
                             g_0_xxxxyyy_0_yyyzz_1,   \
                             g_0_xxxxyyy_0_yyyzzz_0,  \
                             g_0_xxxxyyy_0_yyyzzz_1,  \
                             g_0_xxxxyyy_0_yyzzz_1,   \
                             g_0_xxxxyyy_0_yyzzzz_0,  \
                             g_0_xxxxyyy_0_yyzzzz_1,  \
                             g_0_xxxxyyy_0_yzzzz_1,   \
                             g_0_xxxxyyy_0_yzzzzz_0,  \
                             g_0_xxxxyyy_0_yzzzzz_1,  \
                             g_0_xxxxyyy_0_zzzzz_1,   \
                             g_0_xxxxyyy_0_zzzzzz_0,  \
                             g_0_xxxxyyy_0_zzzzzz_1,  \
                             g_0_xxxxyyyz_0_xxxxxx_0, \
                             g_0_xxxxyyyz_0_xxxxxy_0, \
                             g_0_xxxxyyyz_0_xxxxxz_0, \
                             g_0_xxxxyyyz_0_xxxxyy_0, \
                             g_0_xxxxyyyz_0_xxxxyz_0, \
                             g_0_xxxxyyyz_0_xxxxzz_0, \
                             g_0_xxxxyyyz_0_xxxyyy_0, \
                             g_0_xxxxyyyz_0_xxxyyz_0, \
                             g_0_xxxxyyyz_0_xxxyzz_0, \
                             g_0_xxxxyyyz_0_xxxzzz_0, \
                             g_0_xxxxyyyz_0_xxyyyy_0, \
                             g_0_xxxxyyyz_0_xxyyyz_0, \
                             g_0_xxxxyyyz_0_xxyyzz_0, \
                             g_0_xxxxyyyz_0_xxyzzz_0, \
                             g_0_xxxxyyyz_0_xxzzzz_0, \
                             g_0_xxxxyyyz_0_xyyyyy_0, \
                             g_0_xxxxyyyz_0_xyyyyz_0, \
                             g_0_xxxxyyyz_0_xyyyzz_0, \
                             g_0_xxxxyyyz_0_xyyzzz_0, \
                             g_0_xxxxyyyz_0_xyzzzz_0, \
                             g_0_xxxxyyyz_0_xzzzzz_0, \
                             g_0_xxxxyyyz_0_yyyyyy_0, \
                             g_0_xxxxyyyz_0_yyyyyz_0, \
                             g_0_xxxxyyyz_0_yyyyzz_0, \
                             g_0_xxxxyyyz_0_yyyzzz_0, \
                             g_0_xxxxyyyz_0_yyzzzz_0, \
                             g_0_xxxxyyyz_0_yzzzzz_0, \
                             g_0_xxxxyyyz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_xxxxxx_0[i] = g_0_xxxxyyy_0_xxxxxx_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxy_0[i] = g_0_xxxxyyy_0_xxxxxy_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxxz_0[i] = g_0_xxxxyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxxz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxyy_0[i] = g_0_xxxxyyy_0_xxxxyy_0[i] * pb_z + g_0_xxxxyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxyz_0[i] = g_0_xxxxyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxyz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxzz_0[i] =
            2.0 * g_0_xxxxyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyyy_0[i] = g_0_xxxxyyy_0_xxxyyy_0[i] * pb_z + g_0_xxxxyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyyz_0[i] = g_0_xxxxyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyyz_0[i] * pb_z + g_0_xxxxyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyzz_0[i] =
            2.0 * g_0_xxxxyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxzzz_0[i] =
            3.0 * g_0_xxxxyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyyy_0[i] = g_0_xxxxyyy_0_xxyyyy_0[i] * pb_z + g_0_xxxxyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyyz_0[i] = g_0_xxxxyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyyz_0[i] * pb_z + g_0_xxxxyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyzz_0[i] * pb_z + g_0_xxxxyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyzzz_0[i] =
            3.0 * g_0_xxxxyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxzzzz_0[i] =
            4.0 * g_0_xxxxyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyyy_0[i] = g_0_xxxxyyy_0_xyyyyy_0[i] * pb_z + g_0_xxxxyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyyz_0[i] = g_0_xxxxyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyyz_0[i] * pb_z + g_0_xxxxyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyzz_0[i] =
            2.0 * g_0_xxxxyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyzz_0[i] * pb_z + g_0_xxxxyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyzzz_0[i] =
            3.0 * g_0_xxxxyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyzzz_0[i] * pb_z + g_0_xxxxyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyzzzz_0[i] =
            4.0 * g_0_xxxxyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xzzzzz_0[i] =
            5.0 * g_0_xxxxyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyyy_0[i] = g_0_xxxxyyy_0_yyyyyy_0[i] * pb_z + g_0_xxxxyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyyz_0[i] = g_0_xxxxyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyyyz_0[i] * pb_z + g_0_xxxxyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyzz_0[i] =
            2.0 * g_0_xxxxyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyyzz_0[i] * pb_z + g_0_xxxxyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyzzz_0[i] * pb_z + g_0_xxxxyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyzzzz_0[i] =
            4.0 * g_0_xxxxyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyzzzz_0[i] * pb_z + g_0_xxxxyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yzzzzz_0[i] =
            5.0 * g_0_xxxxyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_zzzzzz_0[i] =
            6.0 * g_0_xxxxyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_zzzzzz_0[i] * pb_z + g_0_xxxxyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 336-364 components of targeted buffer : SLSI

    auto g_0_xxxxyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 336);

    auto g_0_xxxxyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 337);

    auto g_0_xxxxyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 338);

    auto g_0_xxxxyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 339);

    auto g_0_xxxxyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 340);

    auto g_0_xxxxyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 341);

    auto g_0_xxxxyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 342);

    auto g_0_xxxxyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 343);

    auto g_0_xxxxyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 344);

    auto g_0_xxxxyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 345);

    auto g_0_xxxxyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 346);

    auto g_0_xxxxyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 347);

    auto g_0_xxxxyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 348);

    auto g_0_xxxxyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 349);

    auto g_0_xxxxyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 350);

    auto g_0_xxxxyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 351);

    auto g_0_xxxxyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 352);

    auto g_0_xxxxyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 353);

    auto g_0_xxxxyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 354);

    auto g_0_xxxxyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 355);

    auto g_0_xxxxyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 356);

    auto g_0_xxxxyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 357);

    auto g_0_xxxxyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 358);

    auto g_0_xxxxyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 359);

    auto g_0_xxxxyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 360);

    auto g_0_xxxxyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 361);

    auto g_0_xxxxyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 362);

    auto g_0_xxxxyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 363);

#pragma omp simd aligned(g_0_xxxxyy_0_xxxxxy_0,       \
                             g_0_xxxxyy_0_xxxxxy_1,   \
                             g_0_xxxxyy_0_xxxxyy_0,   \
                             g_0_xxxxyy_0_xxxxyy_1,   \
                             g_0_xxxxyy_0_xxxyyy_0,   \
                             g_0_xxxxyy_0_xxxyyy_1,   \
                             g_0_xxxxyy_0_xxyyyy_0,   \
                             g_0_xxxxyy_0_xxyyyy_1,   \
                             g_0_xxxxyy_0_xyyyyy_0,   \
                             g_0_xxxxyy_0_xyyyyy_1,   \
                             g_0_xxxxyyz_0_xxxxxy_0,  \
                             g_0_xxxxyyz_0_xxxxxy_1,  \
                             g_0_xxxxyyz_0_xxxxyy_0,  \
                             g_0_xxxxyyz_0_xxxxyy_1,  \
                             g_0_xxxxyyz_0_xxxyyy_0,  \
                             g_0_xxxxyyz_0_xxxyyy_1,  \
                             g_0_xxxxyyz_0_xxyyyy_0,  \
                             g_0_xxxxyyz_0_xxyyyy_1,  \
                             g_0_xxxxyyz_0_xyyyyy_0,  \
                             g_0_xxxxyyz_0_xyyyyy_1,  \
                             g_0_xxxxyyzz_0_xxxxxx_0, \
                             g_0_xxxxyyzz_0_xxxxxy_0, \
                             g_0_xxxxyyzz_0_xxxxxz_0, \
                             g_0_xxxxyyzz_0_xxxxyy_0, \
                             g_0_xxxxyyzz_0_xxxxyz_0, \
                             g_0_xxxxyyzz_0_xxxxzz_0, \
                             g_0_xxxxyyzz_0_xxxyyy_0, \
                             g_0_xxxxyyzz_0_xxxyyz_0, \
                             g_0_xxxxyyzz_0_xxxyzz_0, \
                             g_0_xxxxyyzz_0_xxxzzz_0, \
                             g_0_xxxxyyzz_0_xxyyyy_0, \
                             g_0_xxxxyyzz_0_xxyyyz_0, \
                             g_0_xxxxyyzz_0_xxyyzz_0, \
                             g_0_xxxxyyzz_0_xxyzzz_0, \
                             g_0_xxxxyyzz_0_xxzzzz_0, \
                             g_0_xxxxyyzz_0_xyyyyy_0, \
                             g_0_xxxxyyzz_0_xyyyyz_0, \
                             g_0_xxxxyyzz_0_xyyyzz_0, \
                             g_0_xxxxyyzz_0_xyyzzz_0, \
                             g_0_xxxxyyzz_0_xyzzzz_0, \
                             g_0_xxxxyyzz_0_xzzzzz_0, \
                             g_0_xxxxyyzz_0_yyyyyy_0, \
                             g_0_xxxxyyzz_0_yyyyyz_0, \
                             g_0_xxxxyyzz_0_yyyyzz_0, \
                             g_0_xxxxyyzz_0_yyyzzz_0, \
                             g_0_xxxxyyzz_0_yyzzzz_0, \
                             g_0_xxxxyyzz_0_yzzzzz_0, \
                             g_0_xxxxyyzz_0_zzzzzz_0, \
                             g_0_xxxxyzz_0_xxxxxx_0,  \
                             g_0_xxxxyzz_0_xxxxxx_1,  \
                             g_0_xxxxyzz_0_xxxxxz_0,  \
                             g_0_xxxxyzz_0_xxxxxz_1,  \
                             g_0_xxxxyzz_0_xxxxzz_0,  \
                             g_0_xxxxyzz_0_xxxxzz_1,  \
                             g_0_xxxxyzz_0_xxxzzz_0,  \
                             g_0_xxxxyzz_0_xxxzzz_1,  \
                             g_0_xxxxyzz_0_xxzzzz_0,  \
                             g_0_xxxxyzz_0_xxzzzz_1,  \
                             g_0_xxxxyzz_0_xzzzzz_0,  \
                             g_0_xxxxyzz_0_xzzzzz_1,  \
                             g_0_xxxxzz_0_xxxxxx_0,   \
                             g_0_xxxxzz_0_xxxxxx_1,   \
                             g_0_xxxxzz_0_xxxxxz_0,   \
                             g_0_xxxxzz_0_xxxxxz_1,   \
                             g_0_xxxxzz_0_xxxxzz_0,   \
                             g_0_xxxxzz_0_xxxxzz_1,   \
                             g_0_xxxxzz_0_xxxzzz_0,   \
                             g_0_xxxxzz_0_xxxzzz_1,   \
                             g_0_xxxxzz_0_xxzzzz_0,   \
                             g_0_xxxxzz_0_xxzzzz_1,   \
                             g_0_xxxxzz_0_xzzzzz_0,   \
                             g_0_xxxxzz_0_xzzzzz_1,   \
                             g_0_xxxyyzz_0_xxxxyz_0,  \
                             g_0_xxxyyzz_0_xxxxyz_1,  \
                             g_0_xxxyyzz_0_xxxyyz_0,  \
                             g_0_xxxyyzz_0_xxxyyz_1,  \
                             g_0_xxxyyzz_0_xxxyz_1,   \
                             g_0_xxxyyzz_0_xxxyzz_0,  \
                             g_0_xxxyyzz_0_xxxyzz_1,  \
                             g_0_xxxyyzz_0_xxyyyz_0,  \
                             g_0_xxxyyzz_0_xxyyyz_1,  \
                             g_0_xxxyyzz_0_xxyyz_1,   \
                             g_0_xxxyyzz_0_xxyyzz_0,  \
                             g_0_xxxyyzz_0_xxyyzz_1,  \
                             g_0_xxxyyzz_0_xxyzz_1,   \
                             g_0_xxxyyzz_0_xxyzzz_0,  \
                             g_0_xxxyyzz_0_xxyzzz_1,  \
                             g_0_xxxyyzz_0_xyyyyz_0,  \
                             g_0_xxxyyzz_0_xyyyyz_1,  \
                             g_0_xxxyyzz_0_xyyyz_1,   \
                             g_0_xxxyyzz_0_xyyyzz_0,  \
                             g_0_xxxyyzz_0_xyyyzz_1,  \
                             g_0_xxxyyzz_0_xyyzz_1,   \
                             g_0_xxxyyzz_0_xyyzzz_0,  \
                             g_0_xxxyyzz_0_xyyzzz_1,  \
                             g_0_xxxyyzz_0_xyzzz_1,   \
                             g_0_xxxyyzz_0_xyzzzz_0,  \
                             g_0_xxxyyzz_0_xyzzzz_1,  \
                             g_0_xxxyyzz_0_yyyyyy_0,  \
                             g_0_xxxyyzz_0_yyyyyy_1,  \
                             g_0_xxxyyzz_0_yyyyyz_0,  \
                             g_0_xxxyyzz_0_yyyyyz_1,  \
                             g_0_xxxyyzz_0_yyyyz_1,   \
                             g_0_xxxyyzz_0_yyyyzz_0,  \
                             g_0_xxxyyzz_0_yyyyzz_1,  \
                             g_0_xxxyyzz_0_yyyzz_1,   \
                             g_0_xxxyyzz_0_yyyzzz_0,  \
                             g_0_xxxyyzz_0_yyyzzz_1,  \
                             g_0_xxxyyzz_0_yyzzz_1,   \
                             g_0_xxxyyzz_0_yyzzzz_0,  \
                             g_0_xxxyyzz_0_yyzzzz_1,  \
                             g_0_xxxyyzz_0_yzzzz_1,   \
                             g_0_xxxyyzz_0_yzzzzz_0,  \
                             g_0_xxxyyzz_0_yzzzzz_1,  \
                             g_0_xxxyyzz_0_zzzzzz_0,  \
                             g_0_xxxyyzz_0_zzzzzz_1,  \
                             g_0_xxyyzz_0_xxxxyz_0,   \
                             g_0_xxyyzz_0_xxxxyz_1,   \
                             g_0_xxyyzz_0_xxxyyz_0,   \
                             g_0_xxyyzz_0_xxxyyz_1,   \
                             g_0_xxyyzz_0_xxxyzz_0,   \
                             g_0_xxyyzz_0_xxxyzz_1,   \
                             g_0_xxyyzz_0_xxyyyz_0,   \
                             g_0_xxyyzz_0_xxyyyz_1,   \
                             g_0_xxyyzz_0_xxyyzz_0,   \
                             g_0_xxyyzz_0_xxyyzz_1,   \
                             g_0_xxyyzz_0_xxyzzz_0,   \
                             g_0_xxyyzz_0_xxyzzz_1,   \
                             g_0_xxyyzz_0_xyyyyz_0,   \
                             g_0_xxyyzz_0_xyyyyz_1,   \
                             g_0_xxyyzz_0_xyyyzz_0,   \
                             g_0_xxyyzz_0_xyyyzz_1,   \
                             g_0_xxyyzz_0_xyyzzz_0,   \
                             g_0_xxyyzz_0_xyyzzz_1,   \
                             g_0_xxyyzz_0_xyzzzz_0,   \
                             g_0_xxyyzz_0_xyzzzz_1,   \
                             g_0_xxyyzz_0_yyyyyy_0,   \
                             g_0_xxyyzz_0_yyyyyy_1,   \
                             g_0_xxyyzz_0_yyyyyz_0,   \
                             g_0_xxyyzz_0_yyyyyz_1,   \
                             g_0_xxyyzz_0_yyyyzz_0,   \
                             g_0_xxyyzz_0_yyyyzz_1,   \
                             g_0_xxyyzz_0_yyyzzz_0,   \
                             g_0_xxyyzz_0_yyyzzz_1,   \
                             g_0_xxyyzz_0_yyzzzz_0,   \
                             g_0_xxyyzz_0_yyzzzz_1,   \
                             g_0_xxyyzz_0_yzzzzz_0,   \
                             g_0_xxyyzz_0_yzzzzz_1,   \
                             g_0_xxyyzz_0_zzzzzz_0,   \
                             g_0_xxyyzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_xxxxxx_0[i] = g_0_xxxxzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxxx_0[i] * pb_y +
                                     g_0_xxxxyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxxxy_0[i] = g_0_xxxxyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxxxy_0[i] * pb_z +
                                     g_0_xxxxyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxxxz_0[i] = g_0_xxxxzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxxz_0[i] * pb_y +
                                     g_0_xxxxyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxxyy_0[i] = g_0_xxxxyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxxyy_0[i] * pb_z +
                                     g_0_xxxxyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxxyz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxxzz_0[i] = g_0_xxxxzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxzz_0[i] * pb_y +
                                     g_0_xxxxyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxyyy_0[i] = g_0_xxxxyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxyyy_0[i] * pb_z +
                                     g_0_xxxxyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxyyz_0[i] = 3.0 * g_0_xxyyzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxyzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxzzz_0[i] = g_0_xxxxzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxzzz_0[i] * pb_y +
                                     g_0_xxxxyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxyyyy_0[i] = g_0_xxxxyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxyyyy_0[i] * pb_z +
                                     g_0_xxxxyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxyyyz_0[i] = 3.0 * g_0_xxyyzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxyyzz_0[i] = 3.0 * g_0_xxyyzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxyzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxzzzz_0[i] = g_0_xxxxzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxzzzz_0[i] * pb_y +
                                     g_0_xxxxyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xyyyyy_0[i] = g_0_xxxxyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xyyyyy_0[i] * pb_z +
                                     g_0_xxxxyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xyyyyz_0[i] = 3.0 * g_0_xxyyzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyyyz_0[i] * pb_x + g_0_xxxyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyyyzz_0[i] = 3.0 * g_0_xxyyzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyyzz_0[i] * pb_x + g_0_xxxyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyyzzz_0[i] = 3.0 * g_0_xxyyzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyzzz_0[i] * pb_x + g_0_xxxyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyzzzz_0[i] * pb_x + g_0_xxxyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xzzzzz_0[i] = g_0_xxxxzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xzzzzz_0[i] * pb_y +
                                     g_0_xxxxyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_yyyyyy_0[i] = 3.0 * g_0_xxyyzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyyyyy_0[i] * pb_x + g_0_xxxyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyyyz_0[i] = 3.0 * g_0_xxyyzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyyyyz_0[i] * pb_x + g_0_xxxyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyyzz_0[i] = 3.0 * g_0_xxyyzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyyyzz_0[i] * pb_x + g_0_xxxyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyzzz_0[i] = 3.0 * g_0_xxyyzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyyzzz_0[i] * pb_x + g_0_xxxyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyzzzz_0[i] = 3.0 * g_0_xxyyzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yyzzzz_0[i] * pb_x + g_0_xxxyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_yzzzzz_0[i] * pb_x + g_0_xxxyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_zzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_zzzzzz_0[i] * pb_x + g_0_xxxyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 364-392 components of targeted buffer : SLSI

    auto g_0_xxxxyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 364);

    auto g_0_xxxxyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 365);

    auto g_0_xxxxyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 366);

    auto g_0_xxxxyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 367);

    auto g_0_xxxxyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 368);

    auto g_0_xxxxyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 369);

    auto g_0_xxxxyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 370);

    auto g_0_xxxxyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 371);

    auto g_0_xxxxyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 372);

    auto g_0_xxxxyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 373);

    auto g_0_xxxxyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 374);

    auto g_0_xxxxyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 375);

    auto g_0_xxxxyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 376);

    auto g_0_xxxxyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 377);

    auto g_0_xxxxyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 378);

    auto g_0_xxxxyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 379);

    auto g_0_xxxxyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 380);

    auto g_0_xxxxyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 381);

    auto g_0_xxxxyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 382);

    auto g_0_xxxxyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 383);

    auto g_0_xxxxyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 384);

    auto g_0_xxxxyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 385);

    auto g_0_xxxxyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 386);

    auto g_0_xxxxyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 387);

    auto g_0_xxxxyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 388);

    auto g_0_xxxxyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 389);

    auto g_0_xxxxyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 390);

    auto g_0_xxxxyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 391);

#pragma omp simd aligned(g_0_xxxxyzzz_0_xxxxxx_0,     \
                             g_0_xxxxyzzz_0_xxxxxy_0, \
                             g_0_xxxxyzzz_0_xxxxxz_0, \
                             g_0_xxxxyzzz_0_xxxxyy_0, \
                             g_0_xxxxyzzz_0_xxxxyz_0, \
                             g_0_xxxxyzzz_0_xxxxzz_0, \
                             g_0_xxxxyzzz_0_xxxyyy_0, \
                             g_0_xxxxyzzz_0_xxxyyz_0, \
                             g_0_xxxxyzzz_0_xxxyzz_0, \
                             g_0_xxxxyzzz_0_xxxzzz_0, \
                             g_0_xxxxyzzz_0_xxyyyy_0, \
                             g_0_xxxxyzzz_0_xxyyyz_0, \
                             g_0_xxxxyzzz_0_xxyyzz_0, \
                             g_0_xxxxyzzz_0_xxyzzz_0, \
                             g_0_xxxxyzzz_0_xxzzzz_0, \
                             g_0_xxxxyzzz_0_xyyyyy_0, \
                             g_0_xxxxyzzz_0_xyyyyz_0, \
                             g_0_xxxxyzzz_0_xyyyzz_0, \
                             g_0_xxxxyzzz_0_xyyzzz_0, \
                             g_0_xxxxyzzz_0_xyzzzz_0, \
                             g_0_xxxxyzzz_0_xzzzzz_0, \
                             g_0_xxxxyzzz_0_yyyyyy_0, \
                             g_0_xxxxyzzz_0_yyyyyz_0, \
                             g_0_xxxxyzzz_0_yyyyzz_0, \
                             g_0_xxxxyzzz_0_yyyzzz_0, \
                             g_0_xxxxyzzz_0_yyzzzz_0, \
                             g_0_xxxxyzzz_0_yzzzzz_0, \
                             g_0_xxxxyzzz_0_zzzzzz_0, \
                             g_0_xxxxzzz_0_xxxxx_1,   \
                             g_0_xxxxzzz_0_xxxxxx_0,  \
                             g_0_xxxxzzz_0_xxxxxx_1,  \
                             g_0_xxxxzzz_0_xxxxxy_0,  \
                             g_0_xxxxzzz_0_xxxxxy_1,  \
                             g_0_xxxxzzz_0_xxxxxz_0,  \
                             g_0_xxxxzzz_0_xxxxxz_1,  \
                             g_0_xxxxzzz_0_xxxxy_1,   \
                             g_0_xxxxzzz_0_xxxxyy_0,  \
                             g_0_xxxxzzz_0_xxxxyy_1,  \
                             g_0_xxxxzzz_0_xxxxyz_0,  \
                             g_0_xxxxzzz_0_xxxxyz_1,  \
                             g_0_xxxxzzz_0_xxxxz_1,   \
                             g_0_xxxxzzz_0_xxxxzz_0,  \
                             g_0_xxxxzzz_0_xxxxzz_1,  \
                             g_0_xxxxzzz_0_xxxyy_1,   \
                             g_0_xxxxzzz_0_xxxyyy_0,  \
                             g_0_xxxxzzz_0_xxxyyy_1,  \
                             g_0_xxxxzzz_0_xxxyyz_0,  \
                             g_0_xxxxzzz_0_xxxyyz_1,  \
                             g_0_xxxxzzz_0_xxxyz_1,   \
                             g_0_xxxxzzz_0_xxxyzz_0,  \
                             g_0_xxxxzzz_0_xxxyzz_1,  \
                             g_0_xxxxzzz_0_xxxzz_1,   \
                             g_0_xxxxzzz_0_xxxzzz_0,  \
                             g_0_xxxxzzz_0_xxxzzz_1,  \
                             g_0_xxxxzzz_0_xxyyy_1,   \
                             g_0_xxxxzzz_0_xxyyyy_0,  \
                             g_0_xxxxzzz_0_xxyyyy_1,  \
                             g_0_xxxxzzz_0_xxyyyz_0,  \
                             g_0_xxxxzzz_0_xxyyyz_1,  \
                             g_0_xxxxzzz_0_xxyyz_1,   \
                             g_0_xxxxzzz_0_xxyyzz_0,  \
                             g_0_xxxxzzz_0_xxyyzz_1,  \
                             g_0_xxxxzzz_0_xxyzz_1,   \
                             g_0_xxxxzzz_0_xxyzzz_0,  \
                             g_0_xxxxzzz_0_xxyzzz_1,  \
                             g_0_xxxxzzz_0_xxzzz_1,   \
                             g_0_xxxxzzz_0_xxzzzz_0,  \
                             g_0_xxxxzzz_0_xxzzzz_1,  \
                             g_0_xxxxzzz_0_xyyyy_1,   \
                             g_0_xxxxzzz_0_xyyyyy_0,  \
                             g_0_xxxxzzz_0_xyyyyy_1,  \
                             g_0_xxxxzzz_0_xyyyyz_0,  \
                             g_0_xxxxzzz_0_xyyyyz_1,  \
                             g_0_xxxxzzz_0_xyyyz_1,   \
                             g_0_xxxxzzz_0_xyyyzz_0,  \
                             g_0_xxxxzzz_0_xyyyzz_1,  \
                             g_0_xxxxzzz_0_xyyzz_1,   \
                             g_0_xxxxzzz_0_xyyzzz_0,  \
                             g_0_xxxxzzz_0_xyyzzz_1,  \
                             g_0_xxxxzzz_0_xyzzz_1,   \
                             g_0_xxxxzzz_0_xyzzzz_0,  \
                             g_0_xxxxzzz_0_xyzzzz_1,  \
                             g_0_xxxxzzz_0_xzzzz_1,   \
                             g_0_xxxxzzz_0_xzzzzz_0,  \
                             g_0_xxxxzzz_0_xzzzzz_1,  \
                             g_0_xxxxzzz_0_yyyyy_1,   \
                             g_0_xxxxzzz_0_yyyyyy_0,  \
                             g_0_xxxxzzz_0_yyyyyy_1,  \
                             g_0_xxxxzzz_0_yyyyyz_0,  \
                             g_0_xxxxzzz_0_yyyyyz_1,  \
                             g_0_xxxxzzz_0_yyyyz_1,   \
                             g_0_xxxxzzz_0_yyyyzz_0,  \
                             g_0_xxxxzzz_0_yyyyzz_1,  \
                             g_0_xxxxzzz_0_yyyzz_1,   \
                             g_0_xxxxzzz_0_yyyzzz_0,  \
                             g_0_xxxxzzz_0_yyyzzz_1,  \
                             g_0_xxxxzzz_0_yyzzz_1,   \
                             g_0_xxxxzzz_0_yyzzzz_0,  \
                             g_0_xxxxzzz_0_yyzzzz_1,  \
                             g_0_xxxxzzz_0_yzzzz_1,   \
                             g_0_xxxxzzz_0_yzzzzz_0,  \
                             g_0_xxxxzzz_0_yzzzzz_1,  \
                             g_0_xxxxzzz_0_zzzzz_1,   \
                             g_0_xxxxzzz_0_zzzzzz_0,  \
                             g_0_xxxxzzz_0_zzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_xxxxxx_0[i] = g_0_xxxxzzz_0_xxxxxx_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxy_0[i] = g_0_xxxxzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxxy_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxxz_0[i] = g_0_xxxxzzz_0_xxxxxz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxyy_0[i] =
            2.0 * g_0_xxxxzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyy_0[i] * pb_y + g_0_xxxxzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxyz_0[i] = g_0_xxxxzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxyz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxzz_0[i] = g_0_xxxxzzz_0_xxxxzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyyy_0[i] =
            3.0 * g_0_xxxxzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyy_0[i] * pb_y + g_0_xxxxzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyyz_0[i] =
            2.0 * g_0_xxxxzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyyz_0[i] * pb_y + g_0_xxxxzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyzz_0[i] = g_0_xxxxzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxzzz_0[i] = g_0_xxxxzzz_0_xxxzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyyy_0[i] =
            4.0 * g_0_xxxxzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyy_0[i] * pb_y + g_0_xxxxzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyyz_0[i] =
            3.0 * g_0_xxxxzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyyz_0[i] * pb_y + g_0_xxxxzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxxzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyzz_0[i] * pb_y + g_0_xxxxzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyzzz_0[i] = g_0_xxxxzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxzzzz_0[i] = g_0_xxxxzzz_0_xxzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyyy_0[i] =
            5.0 * g_0_xxxxzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyy_0[i] * pb_y + g_0_xxxxzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyyz_0[i] =
            4.0 * g_0_xxxxzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyyz_0[i] * pb_y + g_0_xxxxzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyzz_0[i] =
            3.0 * g_0_xxxxzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyzz_0[i] * pb_y + g_0_xxxxzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyzzz_0[i] =
            2.0 * g_0_xxxxzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyzzz_0[i] * pb_y + g_0_xxxxzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyzzzz_0[i] = g_0_xxxxzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xzzzzz_0[i] = g_0_xxxxzzz_0_xzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyyy_0[i] =
            6.0 * g_0_xxxxzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyyy_0[i] * pb_y + g_0_xxxxzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyyz_0[i] =
            5.0 * g_0_xxxxzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyyz_0[i] * pb_y + g_0_xxxxzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyzz_0[i] =
            4.0 * g_0_xxxxzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyzz_0[i] * pb_y + g_0_xxxxzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxxzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyzzz_0[i] * pb_y + g_0_xxxxzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyzzzz_0[i] =
            2.0 * g_0_xxxxzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyzzzz_0[i] * pb_y + g_0_xxxxzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yzzzzz_0[i] = g_0_xxxxzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_zzzzzz_0[i] = g_0_xxxxzzz_0_zzzzzz_0[i] * pb_y + g_0_xxxxzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 392-420 components of targeted buffer : SLSI

    auto g_0_xxxxzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 392);

    auto g_0_xxxxzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 393);

    auto g_0_xxxxzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 394);

    auto g_0_xxxxzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 395);

    auto g_0_xxxxzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 396);

    auto g_0_xxxxzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 397);

    auto g_0_xxxxzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 398);

    auto g_0_xxxxzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 399);

    auto g_0_xxxxzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 400);

    auto g_0_xxxxzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 401);

    auto g_0_xxxxzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 402);

    auto g_0_xxxxzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 403);

    auto g_0_xxxxzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 404);

    auto g_0_xxxxzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 405);

    auto g_0_xxxxzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 406);

    auto g_0_xxxxzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 407);

    auto g_0_xxxxzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 408);

    auto g_0_xxxxzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 409);

    auto g_0_xxxxzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 410);

    auto g_0_xxxxzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 411);

    auto g_0_xxxxzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 412);

    auto g_0_xxxxzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 413);

    auto g_0_xxxxzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 414);

    auto g_0_xxxxzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 415);

    auto g_0_xxxxzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 416);

    auto g_0_xxxxzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 417);

    auto g_0_xxxxzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 418);

    auto g_0_xxxxzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 419);

#pragma omp simd aligned(g_0_xxxxzz_0_xxxxxx_0,       \
                             g_0_xxxxzz_0_xxxxxx_1,   \
                             g_0_xxxxzz_0_xxxxxy_0,   \
                             g_0_xxxxzz_0_xxxxxy_1,   \
                             g_0_xxxxzz_0_xxxxyy_0,   \
                             g_0_xxxxzz_0_xxxxyy_1,   \
                             g_0_xxxxzz_0_xxxyyy_0,   \
                             g_0_xxxxzz_0_xxxyyy_1,   \
                             g_0_xxxxzz_0_xxyyyy_0,   \
                             g_0_xxxxzz_0_xxyyyy_1,   \
                             g_0_xxxxzz_0_xyyyyy_0,   \
                             g_0_xxxxzz_0_xyyyyy_1,   \
                             g_0_xxxxzzz_0_xxxxxx_0,  \
                             g_0_xxxxzzz_0_xxxxxx_1,  \
                             g_0_xxxxzzz_0_xxxxxy_0,  \
                             g_0_xxxxzzz_0_xxxxxy_1,  \
                             g_0_xxxxzzz_0_xxxxyy_0,  \
                             g_0_xxxxzzz_0_xxxxyy_1,  \
                             g_0_xxxxzzz_0_xxxyyy_0,  \
                             g_0_xxxxzzz_0_xxxyyy_1,  \
                             g_0_xxxxzzz_0_xxyyyy_0,  \
                             g_0_xxxxzzz_0_xxyyyy_1,  \
                             g_0_xxxxzzz_0_xyyyyy_0,  \
                             g_0_xxxxzzz_0_xyyyyy_1,  \
                             g_0_xxxxzzzz_0_xxxxxx_0, \
                             g_0_xxxxzzzz_0_xxxxxy_0, \
                             g_0_xxxxzzzz_0_xxxxxz_0, \
                             g_0_xxxxzzzz_0_xxxxyy_0, \
                             g_0_xxxxzzzz_0_xxxxyz_0, \
                             g_0_xxxxzzzz_0_xxxxzz_0, \
                             g_0_xxxxzzzz_0_xxxyyy_0, \
                             g_0_xxxxzzzz_0_xxxyyz_0, \
                             g_0_xxxxzzzz_0_xxxyzz_0, \
                             g_0_xxxxzzzz_0_xxxzzz_0, \
                             g_0_xxxxzzzz_0_xxyyyy_0, \
                             g_0_xxxxzzzz_0_xxyyyz_0, \
                             g_0_xxxxzzzz_0_xxyyzz_0, \
                             g_0_xxxxzzzz_0_xxyzzz_0, \
                             g_0_xxxxzzzz_0_xxzzzz_0, \
                             g_0_xxxxzzzz_0_xyyyyy_0, \
                             g_0_xxxxzzzz_0_xyyyyz_0, \
                             g_0_xxxxzzzz_0_xyyyzz_0, \
                             g_0_xxxxzzzz_0_xyyzzz_0, \
                             g_0_xxxxzzzz_0_xyzzzz_0, \
                             g_0_xxxxzzzz_0_xzzzzz_0, \
                             g_0_xxxxzzzz_0_yyyyyy_0, \
                             g_0_xxxxzzzz_0_yyyyyz_0, \
                             g_0_xxxxzzzz_0_yyyyzz_0, \
                             g_0_xxxxzzzz_0_yyyzzz_0, \
                             g_0_xxxxzzzz_0_yyzzzz_0, \
                             g_0_xxxxzzzz_0_yzzzzz_0, \
                             g_0_xxxxzzzz_0_zzzzzz_0, \
                             g_0_xxxzzzz_0_xxxxxz_0,  \
                             g_0_xxxzzzz_0_xxxxxz_1,  \
                             g_0_xxxzzzz_0_xxxxyz_0,  \
                             g_0_xxxzzzz_0_xxxxyz_1,  \
                             g_0_xxxzzzz_0_xxxxz_1,   \
                             g_0_xxxzzzz_0_xxxxzz_0,  \
                             g_0_xxxzzzz_0_xxxxzz_1,  \
                             g_0_xxxzzzz_0_xxxyyz_0,  \
                             g_0_xxxzzzz_0_xxxyyz_1,  \
                             g_0_xxxzzzz_0_xxxyz_1,   \
                             g_0_xxxzzzz_0_xxxyzz_0,  \
                             g_0_xxxzzzz_0_xxxyzz_1,  \
                             g_0_xxxzzzz_0_xxxzz_1,   \
                             g_0_xxxzzzz_0_xxxzzz_0,  \
                             g_0_xxxzzzz_0_xxxzzz_1,  \
                             g_0_xxxzzzz_0_xxyyyz_0,  \
                             g_0_xxxzzzz_0_xxyyyz_1,  \
                             g_0_xxxzzzz_0_xxyyz_1,   \
                             g_0_xxxzzzz_0_xxyyzz_0,  \
                             g_0_xxxzzzz_0_xxyyzz_1,  \
                             g_0_xxxzzzz_0_xxyzz_1,   \
                             g_0_xxxzzzz_0_xxyzzz_0,  \
                             g_0_xxxzzzz_0_xxyzzz_1,  \
                             g_0_xxxzzzz_0_xxzzz_1,   \
                             g_0_xxxzzzz_0_xxzzzz_0,  \
                             g_0_xxxzzzz_0_xxzzzz_1,  \
                             g_0_xxxzzzz_0_xyyyyz_0,  \
                             g_0_xxxzzzz_0_xyyyyz_1,  \
                             g_0_xxxzzzz_0_xyyyz_1,   \
                             g_0_xxxzzzz_0_xyyyzz_0,  \
                             g_0_xxxzzzz_0_xyyyzz_1,  \
                             g_0_xxxzzzz_0_xyyzz_1,   \
                             g_0_xxxzzzz_0_xyyzzz_0,  \
                             g_0_xxxzzzz_0_xyyzzz_1,  \
                             g_0_xxxzzzz_0_xyzzz_1,   \
                             g_0_xxxzzzz_0_xyzzzz_0,  \
                             g_0_xxxzzzz_0_xyzzzz_1,  \
                             g_0_xxxzzzz_0_xzzzz_1,   \
                             g_0_xxxzzzz_0_xzzzzz_0,  \
                             g_0_xxxzzzz_0_xzzzzz_1,  \
                             g_0_xxxzzzz_0_yyyyyy_0,  \
                             g_0_xxxzzzz_0_yyyyyy_1,  \
                             g_0_xxxzzzz_0_yyyyyz_0,  \
                             g_0_xxxzzzz_0_yyyyyz_1,  \
                             g_0_xxxzzzz_0_yyyyz_1,   \
                             g_0_xxxzzzz_0_yyyyzz_0,  \
                             g_0_xxxzzzz_0_yyyyzz_1,  \
                             g_0_xxxzzzz_0_yyyzz_1,   \
                             g_0_xxxzzzz_0_yyyzzz_0,  \
                             g_0_xxxzzzz_0_yyyzzz_1,  \
                             g_0_xxxzzzz_0_yyzzz_1,   \
                             g_0_xxxzzzz_0_yyzzzz_0,  \
                             g_0_xxxzzzz_0_yyzzzz_1,  \
                             g_0_xxxzzzz_0_yzzzz_1,   \
                             g_0_xxxzzzz_0_yzzzzz_0,  \
                             g_0_xxxzzzz_0_yzzzzz_1,  \
                             g_0_xxxzzzz_0_zzzzz_1,   \
                             g_0_xxxzzzz_0_zzzzzz_0,  \
                             g_0_xxxzzzz_0_zzzzzz_1,  \
                             g_0_xxzzzz_0_xxxxxz_0,   \
                             g_0_xxzzzz_0_xxxxxz_1,   \
                             g_0_xxzzzz_0_xxxxyz_0,   \
                             g_0_xxzzzz_0_xxxxyz_1,   \
                             g_0_xxzzzz_0_xxxxzz_0,   \
                             g_0_xxzzzz_0_xxxxzz_1,   \
                             g_0_xxzzzz_0_xxxyyz_0,   \
                             g_0_xxzzzz_0_xxxyyz_1,   \
                             g_0_xxzzzz_0_xxxyzz_0,   \
                             g_0_xxzzzz_0_xxxyzz_1,   \
                             g_0_xxzzzz_0_xxxzzz_0,   \
                             g_0_xxzzzz_0_xxxzzz_1,   \
                             g_0_xxzzzz_0_xxyyyz_0,   \
                             g_0_xxzzzz_0_xxyyyz_1,   \
                             g_0_xxzzzz_0_xxyyzz_0,   \
                             g_0_xxzzzz_0_xxyyzz_1,   \
                             g_0_xxzzzz_0_xxyzzz_0,   \
                             g_0_xxzzzz_0_xxyzzz_1,   \
                             g_0_xxzzzz_0_xxzzzz_0,   \
                             g_0_xxzzzz_0_xxzzzz_1,   \
                             g_0_xxzzzz_0_xyyyyz_0,   \
                             g_0_xxzzzz_0_xyyyyz_1,   \
                             g_0_xxzzzz_0_xyyyzz_0,   \
                             g_0_xxzzzz_0_xyyyzz_1,   \
                             g_0_xxzzzz_0_xyyzzz_0,   \
                             g_0_xxzzzz_0_xyyzzz_1,   \
                             g_0_xxzzzz_0_xyzzzz_0,   \
                             g_0_xxzzzz_0_xyzzzz_1,   \
                             g_0_xxzzzz_0_xzzzzz_0,   \
                             g_0_xxzzzz_0_xzzzzz_1,   \
                             g_0_xxzzzz_0_yyyyyy_0,   \
                             g_0_xxzzzz_0_yyyyyy_1,   \
                             g_0_xxzzzz_0_yyyyyz_0,   \
                             g_0_xxzzzz_0_yyyyyz_1,   \
                             g_0_xxzzzz_0_yyyyzz_0,   \
                             g_0_xxzzzz_0_yyyyzz_1,   \
                             g_0_xxzzzz_0_yyyzzz_0,   \
                             g_0_xxzzzz_0_yyyzzz_1,   \
                             g_0_xxzzzz_0_yyzzzz_0,   \
                             g_0_xxzzzz_0_yyzzzz_1,   \
                             g_0_xxzzzz_0_yzzzzz_0,   \
                             g_0_xxzzzz_0_yzzzzz_1,   \
                             g_0_xxzzzz_0_zzzzzz_0,   \
                             g_0_xxzzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_xxxxxx_0[i] = 3.0 * g_0_xxxxzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_xxxxxx_0[i] * pb_z + g_0_xxxxzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxxy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_xxxxxy_0[i] * pb_z + g_0_xxxxzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxxz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxxzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_xxxxyy_0[i] * pb_z + g_0_xxxxzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxxzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxxzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxzz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_xxxyyy_0[i] * pb_z + g_0_xxxxzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxxzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxzzz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_xxyyyy_0[i] * pb_z + g_0_xxxxzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxxzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxzzzz_0[i] * pb_x +
                                     g_0_xxxzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyyyy_0[i] = 3.0 * g_0_xxxxzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxxzzz_0_xyyyyy_0[i] * pb_z + g_0_xxxxzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xyyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyz_0[i] * pb_x + g_0_xxxzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyzz_0[i] * pb_x + g_0_xxxzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyzzz_0[i] * pb_x + g_0_xxxzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyyy_0[i] = 3.0 * g_0_xxzzzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyyyyy_0[i] * pb_x + g_0_xxxzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyyz_0[i] = 3.0 * g_0_xxzzzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyyyyz_0[i] * pb_x + g_0_xxxzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyzz_0[i] = 3.0 * g_0_xxzzzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyyyzz_0[i] * pb_x + g_0_xxxzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyyzzz_0[i] * pb_x + g_0_xxxzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yyzzzz_0[i] * pb_x + g_0_xxxzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_yzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_zzzzzz_0[i] = 3.0 * g_0_xxzzzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_zzzzzz_0[i] * pb_x + g_0_xxxzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 420-448 components of targeted buffer : SLSI

    auto g_0_xxxyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 420);

    auto g_0_xxxyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 421);

    auto g_0_xxxyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 422);

    auto g_0_xxxyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 423);

    auto g_0_xxxyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 424);

    auto g_0_xxxyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 425);

    auto g_0_xxxyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 426);

    auto g_0_xxxyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 427);

    auto g_0_xxxyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 428);

    auto g_0_xxxyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 429);

    auto g_0_xxxyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 430);

    auto g_0_xxxyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 431);

    auto g_0_xxxyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 432);

    auto g_0_xxxyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 433);

    auto g_0_xxxyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 434);

    auto g_0_xxxyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 435);

    auto g_0_xxxyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 436);

    auto g_0_xxxyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 437);

    auto g_0_xxxyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 438);

    auto g_0_xxxyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 439);

    auto g_0_xxxyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 440);

    auto g_0_xxxyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 441);

    auto g_0_xxxyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 442);

    auto g_0_xxxyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 443);

    auto g_0_xxxyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 444);

    auto g_0_xxxyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 445);

    auto g_0_xxxyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 446);

    auto g_0_xxxyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 447);

#pragma omp simd aligned(g_0_xxxyyy_0_xxxxxx_0,       \
                             g_0_xxxyyy_0_xxxxxx_1,   \
                             g_0_xxxyyy_0_xxxxxz_0,   \
                             g_0_xxxyyy_0_xxxxxz_1,   \
                             g_0_xxxyyy_0_xxxxzz_0,   \
                             g_0_xxxyyy_0_xxxxzz_1,   \
                             g_0_xxxyyy_0_xxxzzz_0,   \
                             g_0_xxxyyy_0_xxxzzz_1,   \
                             g_0_xxxyyy_0_xxzzzz_0,   \
                             g_0_xxxyyy_0_xxzzzz_1,   \
                             g_0_xxxyyy_0_xzzzzz_0,   \
                             g_0_xxxyyy_0_xzzzzz_1,   \
                             g_0_xxxyyyy_0_xxxxxx_0,  \
                             g_0_xxxyyyy_0_xxxxxx_1,  \
                             g_0_xxxyyyy_0_xxxxxz_0,  \
                             g_0_xxxyyyy_0_xxxxxz_1,  \
                             g_0_xxxyyyy_0_xxxxzz_0,  \
                             g_0_xxxyyyy_0_xxxxzz_1,  \
                             g_0_xxxyyyy_0_xxxzzz_0,  \
                             g_0_xxxyyyy_0_xxxzzz_1,  \
                             g_0_xxxyyyy_0_xxzzzz_0,  \
                             g_0_xxxyyyy_0_xxzzzz_1,  \
                             g_0_xxxyyyy_0_xzzzzz_0,  \
                             g_0_xxxyyyy_0_xzzzzz_1,  \
                             g_0_xxxyyyyy_0_xxxxxx_0, \
                             g_0_xxxyyyyy_0_xxxxxy_0, \
                             g_0_xxxyyyyy_0_xxxxxz_0, \
                             g_0_xxxyyyyy_0_xxxxyy_0, \
                             g_0_xxxyyyyy_0_xxxxyz_0, \
                             g_0_xxxyyyyy_0_xxxxzz_0, \
                             g_0_xxxyyyyy_0_xxxyyy_0, \
                             g_0_xxxyyyyy_0_xxxyyz_0, \
                             g_0_xxxyyyyy_0_xxxyzz_0, \
                             g_0_xxxyyyyy_0_xxxzzz_0, \
                             g_0_xxxyyyyy_0_xxyyyy_0, \
                             g_0_xxxyyyyy_0_xxyyyz_0, \
                             g_0_xxxyyyyy_0_xxyyzz_0, \
                             g_0_xxxyyyyy_0_xxyzzz_0, \
                             g_0_xxxyyyyy_0_xxzzzz_0, \
                             g_0_xxxyyyyy_0_xyyyyy_0, \
                             g_0_xxxyyyyy_0_xyyyyz_0, \
                             g_0_xxxyyyyy_0_xyyyzz_0, \
                             g_0_xxxyyyyy_0_xyyzzz_0, \
                             g_0_xxxyyyyy_0_xyzzzz_0, \
                             g_0_xxxyyyyy_0_xzzzzz_0, \
                             g_0_xxxyyyyy_0_yyyyyy_0, \
                             g_0_xxxyyyyy_0_yyyyyz_0, \
                             g_0_xxxyyyyy_0_yyyyzz_0, \
                             g_0_xxxyyyyy_0_yyyzzz_0, \
                             g_0_xxxyyyyy_0_yyzzzz_0, \
                             g_0_xxxyyyyy_0_yzzzzz_0, \
                             g_0_xxxyyyyy_0_zzzzzz_0, \
                             g_0_xxyyyyy_0_xxxxxy_0,  \
                             g_0_xxyyyyy_0_xxxxxy_1,  \
                             g_0_xxyyyyy_0_xxxxy_1,   \
                             g_0_xxyyyyy_0_xxxxyy_0,  \
                             g_0_xxyyyyy_0_xxxxyy_1,  \
                             g_0_xxyyyyy_0_xxxxyz_0,  \
                             g_0_xxyyyyy_0_xxxxyz_1,  \
                             g_0_xxyyyyy_0_xxxyy_1,   \
                             g_0_xxyyyyy_0_xxxyyy_0,  \
                             g_0_xxyyyyy_0_xxxyyy_1,  \
                             g_0_xxyyyyy_0_xxxyyz_0,  \
                             g_0_xxyyyyy_0_xxxyyz_1,  \
                             g_0_xxyyyyy_0_xxxyz_1,   \
                             g_0_xxyyyyy_0_xxxyzz_0,  \
                             g_0_xxyyyyy_0_xxxyzz_1,  \
                             g_0_xxyyyyy_0_xxyyy_1,   \
                             g_0_xxyyyyy_0_xxyyyy_0,  \
                             g_0_xxyyyyy_0_xxyyyy_1,  \
                             g_0_xxyyyyy_0_xxyyyz_0,  \
                             g_0_xxyyyyy_0_xxyyyz_1,  \
                             g_0_xxyyyyy_0_xxyyz_1,   \
                             g_0_xxyyyyy_0_xxyyzz_0,  \
                             g_0_xxyyyyy_0_xxyyzz_1,  \
                             g_0_xxyyyyy_0_xxyzz_1,   \
                             g_0_xxyyyyy_0_xxyzzz_0,  \
                             g_0_xxyyyyy_0_xxyzzz_1,  \
                             g_0_xxyyyyy_0_xyyyy_1,   \
                             g_0_xxyyyyy_0_xyyyyy_0,  \
                             g_0_xxyyyyy_0_xyyyyy_1,  \
                             g_0_xxyyyyy_0_xyyyyz_0,  \
                             g_0_xxyyyyy_0_xyyyyz_1,  \
                             g_0_xxyyyyy_0_xyyyz_1,   \
                             g_0_xxyyyyy_0_xyyyzz_0,  \
                             g_0_xxyyyyy_0_xyyyzz_1,  \
                             g_0_xxyyyyy_0_xyyzz_1,   \
                             g_0_xxyyyyy_0_xyyzzz_0,  \
                             g_0_xxyyyyy_0_xyyzzz_1,  \
                             g_0_xxyyyyy_0_xyzzz_1,   \
                             g_0_xxyyyyy_0_xyzzzz_0,  \
                             g_0_xxyyyyy_0_xyzzzz_1,  \
                             g_0_xxyyyyy_0_yyyyy_1,   \
                             g_0_xxyyyyy_0_yyyyyy_0,  \
                             g_0_xxyyyyy_0_yyyyyy_1,  \
                             g_0_xxyyyyy_0_yyyyyz_0,  \
                             g_0_xxyyyyy_0_yyyyyz_1,  \
                             g_0_xxyyyyy_0_yyyyz_1,   \
                             g_0_xxyyyyy_0_yyyyzz_0,  \
                             g_0_xxyyyyy_0_yyyyzz_1,  \
                             g_0_xxyyyyy_0_yyyzz_1,   \
                             g_0_xxyyyyy_0_yyyzzz_0,  \
                             g_0_xxyyyyy_0_yyyzzz_1,  \
                             g_0_xxyyyyy_0_yyzzz_1,   \
                             g_0_xxyyyyy_0_yyzzzz_0,  \
                             g_0_xxyyyyy_0_yyzzzz_1,  \
                             g_0_xxyyyyy_0_yzzzz_1,   \
                             g_0_xxyyyyy_0_yzzzzz_0,  \
                             g_0_xxyyyyy_0_yzzzzz_1,  \
                             g_0_xxyyyyy_0_zzzzzz_0,  \
                             g_0_xxyyyyy_0_zzzzzz_1,  \
                             g_0_xyyyyy_0_xxxxxy_0,   \
                             g_0_xyyyyy_0_xxxxxy_1,   \
                             g_0_xyyyyy_0_xxxxyy_0,   \
                             g_0_xyyyyy_0_xxxxyy_1,   \
                             g_0_xyyyyy_0_xxxxyz_0,   \
                             g_0_xyyyyy_0_xxxxyz_1,   \
                             g_0_xyyyyy_0_xxxyyy_0,   \
                             g_0_xyyyyy_0_xxxyyy_1,   \
                             g_0_xyyyyy_0_xxxyyz_0,   \
                             g_0_xyyyyy_0_xxxyyz_1,   \
                             g_0_xyyyyy_0_xxxyzz_0,   \
                             g_0_xyyyyy_0_xxxyzz_1,   \
                             g_0_xyyyyy_0_xxyyyy_0,   \
                             g_0_xyyyyy_0_xxyyyy_1,   \
                             g_0_xyyyyy_0_xxyyyz_0,   \
                             g_0_xyyyyy_0_xxyyyz_1,   \
                             g_0_xyyyyy_0_xxyyzz_0,   \
                             g_0_xyyyyy_0_xxyyzz_1,   \
                             g_0_xyyyyy_0_xxyzzz_0,   \
                             g_0_xyyyyy_0_xxyzzz_1,   \
                             g_0_xyyyyy_0_xyyyyy_0,   \
                             g_0_xyyyyy_0_xyyyyy_1,   \
                             g_0_xyyyyy_0_xyyyyz_0,   \
                             g_0_xyyyyy_0_xyyyyz_1,   \
                             g_0_xyyyyy_0_xyyyzz_0,   \
                             g_0_xyyyyy_0_xyyyzz_1,   \
                             g_0_xyyyyy_0_xyyzzz_0,   \
                             g_0_xyyyyy_0_xyyzzz_1,   \
                             g_0_xyyyyy_0_xyzzzz_0,   \
                             g_0_xyyyyy_0_xyzzzz_1,   \
                             g_0_xyyyyy_0_yyyyyy_0,   \
                             g_0_xyyyyy_0_yyyyyy_1,   \
                             g_0_xyyyyy_0_yyyyyz_0,   \
                             g_0_xyyyyy_0_yyyyyz_1,   \
                             g_0_xyyyyy_0_yyyyzz_0,   \
                             g_0_xyyyyy_0_yyyyzz_1,   \
                             g_0_xyyyyy_0_yyyzzz_0,   \
                             g_0_xyyyyy_0_yyyzzz_1,   \
                             g_0_xyyyyy_0_yyzzzz_0,   \
                             g_0_xyyyyy_0_yyzzzz_1,   \
                             g_0_xyyyyy_0_yzzzzz_0,   \
                             g_0_xyyyyy_0_yzzzzz_1,   \
                             g_0_xyyyyy_0_zzzzzz_0,   \
                             g_0_xyyyyy_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_xxxxxx_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_xxxxxx_0[i] * pb_y + g_0_xxxyyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxxxy_0[i] = 2.0 * g_0_xyyyyy_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxy_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxxz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_xxxxxz_0[i] * pb_y + g_0_xxxyyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxxyy_0[i] = 2.0 * g_0_xyyyyy_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyy_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxyz_0[i] = 2.0 * g_0_xyyyyy_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_xxxxzz_0[i] * pb_y + g_0_xxxyyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxyyy_0[i] = 2.0 * g_0_xyyyyy_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyy_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxyyz_0[i] = 2.0 * g_0_xyyyyy_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxyzz_0[i] = 2.0 * g_0_xyyyyy_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_xxxzzz_0[i] * pb_y + g_0_xxxyyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxyyyy_0[i] = 2.0 * g_0_xyyyyy_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyy_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyyyz_0[i] = 2.0 * g_0_xyyyyy_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyyzz_0[i] = 2.0 * g_0_xyyyyy_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyzzz_0[i] = 2.0 * g_0_xyyyyy_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxyyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_xxzzzz_0[i] * pb_y + g_0_xxxyyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xyyyyy_0[i] = 2.0 * g_0_xyyyyy_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyy_0[i] * pb_x + g_0_xxyyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyyyz_0[i] = 2.0 * g_0_xyyyyy_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyz_0[i] * pb_x + g_0_xxyyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyyzz_0[i] = 2.0 * g_0_xyyyyy_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyzz_0[i] * pb_x + g_0_xxyyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyzzz_0[i] = 2.0 * g_0_xyyyyy_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyzzz_0[i] * pb_x + g_0_xxyyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyzzzz_0[i] = 2.0 * g_0_xyyyyy_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzzzz_0[i] * pb_x + g_0_xxyyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xzzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyyy_0_xzzzzz_0[i] * pb_y + g_0_xxxyyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_yyyyyy_0[i] = 2.0 * g_0_xyyyyy_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyyyy_0[i] * pb_x + g_0_xxyyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyyyz_0[i] = 2.0 * g_0_xyyyyy_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyyyz_0[i] * pb_x + g_0_xxyyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyyzz_0[i] = 2.0 * g_0_xyyyyy_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyyzz_0[i] * pb_x + g_0_xxyyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyzzz_0[i] = 2.0 * g_0_xyyyyy_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyyzzz_0[i] * pb_x + g_0_xxyyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyzzzz_0[i] = 2.0 * g_0_xyyyyy_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yyzzzz_0[i] * pb_x + g_0_xxyyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yzzzzz_0[i] = 2.0 * g_0_xyyyyy_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_yzzzzz_0[i] * pb_x + g_0_xxyyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_zzzzzz_0[i] = 2.0 * g_0_xyyyyy_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_zzzzzz_0[i] * pb_x + g_0_xxyyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 448-476 components of targeted buffer : SLSI

    auto g_0_xxxyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 448);

    auto g_0_xxxyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 449);

    auto g_0_xxxyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 450);

    auto g_0_xxxyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 451);

    auto g_0_xxxyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 452);

    auto g_0_xxxyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 453);

    auto g_0_xxxyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 454);

    auto g_0_xxxyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 455);

    auto g_0_xxxyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 456);

    auto g_0_xxxyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 457);

    auto g_0_xxxyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 458);

    auto g_0_xxxyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 459);

    auto g_0_xxxyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 460);

    auto g_0_xxxyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 461);

    auto g_0_xxxyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 462);

    auto g_0_xxxyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 463);

    auto g_0_xxxyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 464);

    auto g_0_xxxyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 465);

    auto g_0_xxxyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 466);

    auto g_0_xxxyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 467);

    auto g_0_xxxyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 468);

    auto g_0_xxxyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 469);

    auto g_0_xxxyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 470);

    auto g_0_xxxyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 471);

    auto g_0_xxxyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 472);

    auto g_0_xxxyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 473);

    auto g_0_xxxyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 474);

    auto g_0_xxxyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 475);

#pragma omp simd aligned(g_0_xxxyyyy_0_xxxxx_1,       \
                             g_0_xxxyyyy_0_xxxxxx_0,  \
                             g_0_xxxyyyy_0_xxxxxx_1,  \
                             g_0_xxxyyyy_0_xxxxxy_0,  \
                             g_0_xxxyyyy_0_xxxxxy_1,  \
                             g_0_xxxyyyy_0_xxxxxz_0,  \
                             g_0_xxxyyyy_0_xxxxxz_1,  \
                             g_0_xxxyyyy_0_xxxxy_1,   \
                             g_0_xxxyyyy_0_xxxxyy_0,  \
                             g_0_xxxyyyy_0_xxxxyy_1,  \
                             g_0_xxxyyyy_0_xxxxyz_0,  \
                             g_0_xxxyyyy_0_xxxxyz_1,  \
                             g_0_xxxyyyy_0_xxxxz_1,   \
                             g_0_xxxyyyy_0_xxxxzz_0,  \
                             g_0_xxxyyyy_0_xxxxzz_1,  \
                             g_0_xxxyyyy_0_xxxyy_1,   \
                             g_0_xxxyyyy_0_xxxyyy_0,  \
                             g_0_xxxyyyy_0_xxxyyy_1,  \
                             g_0_xxxyyyy_0_xxxyyz_0,  \
                             g_0_xxxyyyy_0_xxxyyz_1,  \
                             g_0_xxxyyyy_0_xxxyz_1,   \
                             g_0_xxxyyyy_0_xxxyzz_0,  \
                             g_0_xxxyyyy_0_xxxyzz_1,  \
                             g_0_xxxyyyy_0_xxxzz_1,   \
                             g_0_xxxyyyy_0_xxxzzz_0,  \
                             g_0_xxxyyyy_0_xxxzzz_1,  \
                             g_0_xxxyyyy_0_xxyyy_1,   \
                             g_0_xxxyyyy_0_xxyyyy_0,  \
                             g_0_xxxyyyy_0_xxyyyy_1,  \
                             g_0_xxxyyyy_0_xxyyyz_0,  \
                             g_0_xxxyyyy_0_xxyyyz_1,  \
                             g_0_xxxyyyy_0_xxyyz_1,   \
                             g_0_xxxyyyy_0_xxyyzz_0,  \
                             g_0_xxxyyyy_0_xxyyzz_1,  \
                             g_0_xxxyyyy_0_xxyzz_1,   \
                             g_0_xxxyyyy_0_xxyzzz_0,  \
                             g_0_xxxyyyy_0_xxyzzz_1,  \
                             g_0_xxxyyyy_0_xxzzz_1,   \
                             g_0_xxxyyyy_0_xxzzzz_0,  \
                             g_0_xxxyyyy_0_xxzzzz_1,  \
                             g_0_xxxyyyy_0_xyyyy_1,   \
                             g_0_xxxyyyy_0_xyyyyy_0,  \
                             g_0_xxxyyyy_0_xyyyyy_1,  \
                             g_0_xxxyyyy_0_xyyyyz_0,  \
                             g_0_xxxyyyy_0_xyyyyz_1,  \
                             g_0_xxxyyyy_0_xyyyz_1,   \
                             g_0_xxxyyyy_0_xyyyzz_0,  \
                             g_0_xxxyyyy_0_xyyyzz_1,  \
                             g_0_xxxyyyy_0_xyyzz_1,   \
                             g_0_xxxyyyy_0_xyyzzz_0,  \
                             g_0_xxxyyyy_0_xyyzzz_1,  \
                             g_0_xxxyyyy_0_xyzzz_1,   \
                             g_0_xxxyyyy_0_xyzzzz_0,  \
                             g_0_xxxyyyy_0_xyzzzz_1,  \
                             g_0_xxxyyyy_0_xzzzz_1,   \
                             g_0_xxxyyyy_0_xzzzzz_0,  \
                             g_0_xxxyyyy_0_xzzzzz_1,  \
                             g_0_xxxyyyy_0_yyyyy_1,   \
                             g_0_xxxyyyy_0_yyyyyy_0,  \
                             g_0_xxxyyyy_0_yyyyyy_1,  \
                             g_0_xxxyyyy_0_yyyyyz_0,  \
                             g_0_xxxyyyy_0_yyyyyz_1,  \
                             g_0_xxxyyyy_0_yyyyz_1,   \
                             g_0_xxxyyyy_0_yyyyzz_0,  \
                             g_0_xxxyyyy_0_yyyyzz_1,  \
                             g_0_xxxyyyy_0_yyyzz_1,   \
                             g_0_xxxyyyy_0_yyyzzz_0,  \
                             g_0_xxxyyyy_0_yyyzzz_1,  \
                             g_0_xxxyyyy_0_yyzzz_1,   \
                             g_0_xxxyyyy_0_yyzzzz_0,  \
                             g_0_xxxyyyy_0_yyzzzz_1,  \
                             g_0_xxxyyyy_0_yzzzz_1,   \
                             g_0_xxxyyyy_0_yzzzzz_0,  \
                             g_0_xxxyyyy_0_yzzzzz_1,  \
                             g_0_xxxyyyy_0_zzzzz_1,   \
                             g_0_xxxyyyy_0_zzzzzz_0,  \
                             g_0_xxxyyyy_0_zzzzzz_1,  \
                             g_0_xxxyyyyz_0_xxxxxx_0, \
                             g_0_xxxyyyyz_0_xxxxxy_0, \
                             g_0_xxxyyyyz_0_xxxxxz_0, \
                             g_0_xxxyyyyz_0_xxxxyy_0, \
                             g_0_xxxyyyyz_0_xxxxyz_0, \
                             g_0_xxxyyyyz_0_xxxxzz_0, \
                             g_0_xxxyyyyz_0_xxxyyy_0, \
                             g_0_xxxyyyyz_0_xxxyyz_0, \
                             g_0_xxxyyyyz_0_xxxyzz_0, \
                             g_0_xxxyyyyz_0_xxxzzz_0, \
                             g_0_xxxyyyyz_0_xxyyyy_0, \
                             g_0_xxxyyyyz_0_xxyyyz_0, \
                             g_0_xxxyyyyz_0_xxyyzz_0, \
                             g_0_xxxyyyyz_0_xxyzzz_0, \
                             g_0_xxxyyyyz_0_xxzzzz_0, \
                             g_0_xxxyyyyz_0_xyyyyy_0, \
                             g_0_xxxyyyyz_0_xyyyyz_0, \
                             g_0_xxxyyyyz_0_xyyyzz_0, \
                             g_0_xxxyyyyz_0_xyyzzz_0, \
                             g_0_xxxyyyyz_0_xyzzzz_0, \
                             g_0_xxxyyyyz_0_xzzzzz_0, \
                             g_0_xxxyyyyz_0_yyyyyy_0, \
                             g_0_xxxyyyyz_0_yyyyyz_0, \
                             g_0_xxxyyyyz_0_yyyyzz_0, \
                             g_0_xxxyyyyz_0_yyyzzz_0, \
                             g_0_xxxyyyyz_0_yyzzzz_0, \
                             g_0_xxxyyyyz_0_yzzzzz_0, \
                             g_0_xxxyyyyz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_xxxxxx_0[i] = g_0_xxxyyyy_0_xxxxxx_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxy_0[i] = g_0_xxxyyyy_0_xxxxxy_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxxz_0[i] = g_0_xxxyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxxz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxyy_0[i] = g_0_xxxyyyy_0_xxxxyy_0[i] * pb_z + g_0_xxxyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxyz_0[i] = g_0_xxxyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxyz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxzz_0[i] =
            2.0 * g_0_xxxyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyyy_0[i] = g_0_xxxyyyy_0_xxxyyy_0[i] * pb_z + g_0_xxxyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyyz_0[i] = g_0_xxxyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyyz_0[i] * pb_z + g_0_xxxyyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyzz_0[i] =
            2.0 * g_0_xxxyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxzzz_0[i] =
            3.0 * g_0_xxxyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyyy_0[i] = g_0_xxxyyyy_0_xxyyyy_0[i] * pb_z + g_0_xxxyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyyz_0[i] = g_0_xxxyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyyz_0[i] * pb_z + g_0_xxxyyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyzz_0[i] * pb_z + g_0_xxxyyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyzzz_0[i] =
            3.0 * g_0_xxxyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxzzzz_0[i] =
            4.0 * g_0_xxxyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyyy_0[i] = g_0_xxxyyyy_0_xyyyyy_0[i] * pb_z + g_0_xxxyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyyz_0[i] = g_0_xxxyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyyz_0[i] * pb_z + g_0_xxxyyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyzz_0[i] =
            2.0 * g_0_xxxyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyzz_0[i] * pb_z + g_0_xxxyyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyzzz_0[i] =
            3.0 * g_0_xxxyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyzzz_0[i] * pb_z + g_0_xxxyyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyzzzz_0[i] =
            4.0 * g_0_xxxyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xzzzzz_0[i] =
            5.0 * g_0_xxxyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyyy_0[i] = g_0_xxxyyyy_0_yyyyyy_0[i] * pb_z + g_0_xxxyyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyyz_0[i] = g_0_xxxyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyyyz_0[i] * pb_z + g_0_xxxyyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyzz_0[i] =
            2.0 * g_0_xxxyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyyzz_0[i] * pb_z + g_0_xxxyyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyzzz_0[i] * pb_z + g_0_xxxyyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyzzzz_0[i] =
            4.0 * g_0_xxxyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyzzzz_0[i] * pb_z + g_0_xxxyyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yzzzzz_0[i] =
            5.0 * g_0_xxxyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_zzzzzz_0[i] =
            6.0 * g_0_xxxyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_zzzzzz_0[i] * pb_z + g_0_xxxyyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 476-504 components of targeted buffer : SLSI

    auto g_0_xxxyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 476);

    auto g_0_xxxyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 477);

    auto g_0_xxxyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 478);

    auto g_0_xxxyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 479);

    auto g_0_xxxyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 480);

    auto g_0_xxxyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 481);

    auto g_0_xxxyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 482);

    auto g_0_xxxyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 483);

    auto g_0_xxxyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 484);

    auto g_0_xxxyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 485);

    auto g_0_xxxyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 486);

    auto g_0_xxxyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 487);

    auto g_0_xxxyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 488);

    auto g_0_xxxyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 489);

    auto g_0_xxxyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 490);

    auto g_0_xxxyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 491);

    auto g_0_xxxyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 492);

    auto g_0_xxxyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 493);

    auto g_0_xxxyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 494);

    auto g_0_xxxyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 495);

    auto g_0_xxxyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 496);

    auto g_0_xxxyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 497);

    auto g_0_xxxyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 498);

    auto g_0_xxxyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 499);

    auto g_0_xxxyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 500);

    auto g_0_xxxyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 501);

    auto g_0_xxxyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 502);

    auto g_0_xxxyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 503);

#pragma omp simd aligned(g_0_xxxyyy_0_xxxxxy_0,       \
                             g_0_xxxyyy_0_xxxxxy_1,   \
                             g_0_xxxyyy_0_xxxxyy_0,   \
                             g_0_xxxyyy_0_xxxxyy_1,   \
                             g_0_xxxyyy_0_xxxyyy_0,   \
                             g_0_xxxyyy_0_xxxyyy_1,   \
                             g_0_xxxyyy_0_xxyyyy_0,   \
                             g_0_xxxyyy_0_xxyyyy_1,   \
                             g_0_xxxyyy_0_xyyyyy_0,   \
                             g_0_xxxyyy_0_xyyyyy_1,   \
                             g_0_xxxyyyz_0_xxxxxy_0,  \
                             g_0_xxxyyyz_0_xxxxxy_1,  \
                             g_0_xxxyyyz_0_xxxxyy_0,  \
                             g_0_xxxyyyz_0_xxxxyy_1,  \
                             g_0_xxxyyyz_0_xxxyyy_0,  \
                             g_0_xxxyyyz_0_xxxyyy_1,  \
                             g_0_xxxyyyz_0_xxyyyy_0,  \
                             g_0_xxxyyyz_0_xxyyyy_1,  \
                             g_0_xxxyyyz_0_xyyyyy_0,  \
                             g_0_xxxyyyz_0_xyyyyy_1,  \
                             g_0_xxxyyyzz_0_xxxxxx_0, \
                             g_0_xxxyyyzz_0_xxxxxy_0, \
                             g_0_xxxyyyzz_0_xxxxxz_0, \
                             g_0_xxxyyyzz_0_xxxxyy_0, \
                             g_0_xxxyyyzz_0_xxxxyz_0, \
                             g_0_xxxyyyzz_0_xxxxzz_0, \
                             g_0_xxxyyyzz_0_xxxyyy_0, \
                             g_0_xxxyyyzz_0_xxxyyz_0, \
                             g_0_xxxyyyzz_0_xxxyzz_0, \
                             g_0_xxxyyyzz_0_xxxzzz_0, \
                             g_0_xxxyyyzz_0_xxyyyy_0, \
                             g_0_xxxyyyzz_0_xxyyyz_0, \
                             g_0_xxxyyyzz_0_xxyyzz_0, \
                             g_0_xxxyyyzz_0_xxyzzz_0, \
                             g_0_xxxyyyzz_0_xxzzzz_0, \
                             g_0_xxxyyyzz_0_xyyyyy_0, \
                             g_0_xxxyyyzz_0_xyyyyz_0, \
                             g_0_xxxyyyzz_0_xyyyzz_0, \
                             g_0_xxxyyyzz_0_xyyzzz_0, \
                             g_0_xxxyyyzz_0_xyzzzz_0, \
                             g_0_xxxyyyzz_0_xzzzzz_0, \
                             g_0_xxxyyyzz_0_yyyyyy_0, \
                             g_0_xxxyyyzz_0_yyyyyz_0, \
                             g_0_xxxyyyzz_0_yyyyzz_0, \
                             g_0_xxxyyyzz_0_yyyzzz_0, \
                             g_0_xxxyyyzz_0_yyzzzz_0, \
                             g_0_xxxyyyzz_0_yzzzzz_0, \
                             g_0_xxxyyyzz_0_zzzzzz_0, \
                             g_0_xxxyyzz_0_xxxxxx_0,  \
                             g_0_xxxyyzz_0_xxxxxx_1,  \
                             g_0_xxxyyzz_0_xxxxxz_0,  \
                             g_0_xxxyyzz_0_xxxxxz_1,  \
                             g_0_xxxyyzz_0_xxxxzz_0,  \
                             g_0_xxxyyzz_0_xxxxzz_1,  \
                             g_0_xxxyyzz_0_xxxzzz_0,  \
                             g_0_xxxyyzz_0_xxxzzz_1,  \
                             g_0_xxxyyzz_0_xxzzzz_0,  \
                             g_0_xxxyyzz_0_xxzzzz_1,  \
                             g_0_xxxyyzz_0_xzzzzz_0,  \
                             g_0_xxxyyzz_0_xzzzzz_1,  \
                             g_0_xxxyzz_0_xxxxxx_0,   \
                             g_0_xxxyzz_0_xxxxxx_1,   \
                             g_0_xxxyzz_0_xxxxxz_0,   \
                             g_0_xxxyzz_0_xxxxxz_1,   \
                             g_0_xxxyzz_0_xxxxzz_0,   \
                             g_0_xxxyzz_0_xxxxzz_1,   \
                             g_0_xxxyzz_0_xxxzzz_0,   \
                             g_0_xxxyzz_0_xxxzzz_1,   \
                             g_0_xxxyzz_0_xxzzzz_0,   \
                             g_0_xxxyzz_0_xxzzzz_1,   \
                             g_0_xxxyzz_0_xzzzzz_0,   \
                             g_0_xxxyzz_0_xzzzzz_1,   \
                             g_0_xxyyyzz_0_xxxxyz_0,  \
                             g_0_xxyyyzz_0_xxxxyz_1,  \
                             g_0_xxyyyzz_0_xxxyyz_0,  \
                             g_0_xxyyyzz_0_xxxyyz_1,  \
                             g_0_xxyyyzz_0_xxxyz_1,   \
                             g_0_xxyyyzz_0_xxxyzz_0,  \
                             g_0_xxyyyzz_0_xxxyzz_1,  \
                             g_0_xxyyyzz_0_xxyyyz_0,  \
                             g_0_xxyyyzz_0_xxyyyz_1,  \
                             g_0_xxyyyzz_0_xxyyz_1,   \
                             g_0_xxyyyzz_0_xxyyzz_0,  \
                             g_0_xxyyyzz_0_xxyyzz_1,  \
                             g_0_xxyyyzz_0_xxyzz_1,   \
                             g_0_xxyyyzz_0_xxyzzz_0,  \
                             g_0_xxyyyzz_0_xxyzzz_1,  \
                             g_0_xxyyyzz_0_xyyyyz_0,  \
                             g_0_xxyyyzz_0_xyyyyz_1,  \
                             g_0_xxyyyzz_0_xyyyz_1,   \
                             g_0_xxyyyzz_0_xyyyzz_0,  \
                             g_0_xxyyyzz_0_xyyyzz_1,  \
                             g_0_xxyyyzz_0_xyyzz_1,   \
                             g_0_xxyyyzz_0_xyyzzz_0,  \
                             g_0_xxyyyzz_0_xyyzzz_1,  \
                             g_0_xxyyyzz_0_xyzzz_1,   \
                             g_0_xxyyyzz_0_xyzzzz_0,  \
                             g_0_xxyyyzz_0_xyzzzz_1,  \
                             g_0_xxyyyzz_0_yyyyyy_0,  \
                             g_0_xxyyyzz_0_yyyyyy_1,  \
                             g_0_xxyyyzz_0_yyyyyz_0,  \
                             g_0_xxyyyzz_0_yyyyyz_1,  \
                             g_0_xxyyyzz_0_yyyyz_1,   \
                             g_0_xxyyyzz_0_yyyyzz_0,  \
                             g_0_xxyyyzz_0_yyyyzz_1,  \
                             g_0_xxyyyzz_0_yyyzz_1,   \
                             g_0_xxyyyzz_0_yyyzzz_0,  \
                             g_0_xxyyyzz_0_yyyzzz_1,  \
                             g_0_xxyyyzz_0_yyzzz_1,   \
                             g_0_xxyyyzz_0_yyzzzz_0,  \
                             g_0_xxyyyzz_0_yyzzzz_1,  \
                             g_0_xxyyyzz_0_yzzzz_1,   \
                             g_0_xxyyyzz_0_yzzzzz_0,  \
                             g_0_xxyyyzz_0_yzzzzz_1,  \
                             g_0_xxyyyzz_0_zzzzzz_0,  \
                             g_0_xxyyyzz_0_zzzzzz_1,  \
                             g_0_xyyyzz_0_xxxxyz_0,   \
                             g_0_xyyyzz_0_xxxxyz_1,   \
                             g_0_xyyyzz_0_xxxyyz_0,   \
                             g_0_xyyyzz_0_xxxyyz_1,   \
                             g_0_xyyyzz_0_xxxyzz_0,   \
                             g_0_xyyyzz_0_xxxyzz_1,   \
                             g_0_xyyyzz_0_xxyyyz_0,   \
                             g_0_xyyyzz_0_xxyyyz_1,   \
                             g_0_xyyyzz_0_xxyyzz_0,   \
                             g_0_xyyyzz_0_xxyyzz_1,   \
                             g_0_xyyyzz_0_xxyzzz_0,   \
                             g_0_xyyyzz_0_xxyzzz_1,   \
                             g_0_xyyyzz_0_xyyyyz_0,   \
                             g_0_xyyyzz_0_xyyyyz_1,   \
                             g_0_xyyyzz_0_xyyyzz_0,   \
                             g_0_xyyyzz_0_xyyyzz_1,   \
                             g_0_xyyyzz_0_xyyzzz_0,   \
                             g_0_xyyyzz_0_xyyzzz_1,   \
                             g_0_xyyyzz_0_xyzzzz_0,   \
                             g_0_xyyyzz_0_xyzzzz_1,   \
                             g_0_xyyyzz_0_yyyyyy_0,   \
                             g_0_xyyyzz_0_yyyyyy_1,   \
                             g_0_xyyyzz_0_yyyyyz_0,   \
                             g_0_xyyyzz_0_yyyyyz_1,   \
                             g_0_xyyyzz_0_yyyyzz_0,   \
                             g_0_xyyyzz_0_yyyyzz_1,   \
                             g_0_xyyyzz_0_yyyzzz_0,   \
                             g_0_xyyyzz_0_yyyzzz_1,   \
                             g_0_xyyyzz_0_yyzzzz_0,   \
                             g_0_xyyyzz_0_yyzzzz_1,   \
                             g_0_xyyyzz_0_yzzzzz_0,   \
                             g_0_xyyyzz_0_yzzzzz_1,   \
                             g_0_xyyyzz_0_zzzzzz_0,   \
                             g_0_xyyyzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_xxxxxx_0[i] = 2.0 * g_0_xxxyzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxxxx_0[i] * pb_y + g_0_xxxyyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxxxy_0[i] = g_0_xxxyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxxxy_0[i] * pb_z +
                                     g_0_xxxyyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxxxz_0[i] = 2.0 * g_0_xxxyzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxxxz_0[i] * pb_y + g_0_xxxyyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxxyy_0[i] = g_0_xxxyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxxyy_0[i] * pb_z +
                                     g_0_xxxyyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxxyz_0[i] = 2.0 * g_0_xyyyzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxyyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxxzz_0[i] = 2.0 * g_0_xxxyzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxxzz_0[i] * pb_y + g_0_xxxyyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxyyy_0[i] = g_0_xxxyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxyyy_0[i] * pb_z +
                                     g_0_xxxyyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxyyz_0[i] = 2.0 * g_0_xyyyzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxyyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxyzz_0[i] = 2.0 * g_0_xyyyzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxyyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxzzz_0[i] = 2.0 * g_0_xxxyzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxzzz_0[i] * pb_y + g_0_xxxyyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxyyyy_0[i] = g_0_xxxyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxyyyy_0[i] * pb_z +
                                     g_0_xxxyyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxyyyz_0[i] = 2.0 * g_0_xyyyzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxyyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxyyzz_0[i] = 2.0 * g_0_xyyyzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxyyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxyzzz_0[i] = 2.0 * g_0_xyyyzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxyyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxzzzz_0[i] = 2.0 * g_0_xxxyzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxzzzz_0[i] * pb_y + g_0_xxxyyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xyyyyy_0[i] = g_0_xxxyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xyyyyy_0[i] * pb_z +
                                     g_0_xxxyyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xyyyyz_0[i] = 2.0 * g_0_xyyyzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyyyz_0[i] * pb_x + g_0_xxyyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyyyzz_0[i] = 2.0 * g_0_xyyyzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyyzz_0[i] * pb_x + g_0_xxyyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyyzzz_0[i] = 2.0 * g_0_xyyyzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyzzz_0[i] * pb_x + g_0_xxyyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyzzzz_0[i] = 2.0 * g_0_xyyyzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyzzzz_0[i] * pb_x + g_0_xxyyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xzzzzz_0[i] = 2.0 * g_0_xxxyzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xzzzzz_0[i] * pb_y + g_0_xxxyyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_yyyyyy_0[i] = 2.0 * g_0_xyyyzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyyyyy_0[i] * pb_x + g_0_xxyyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyyyz_0[i] = 2.0 * g_0_xyyyzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyyyyz_0[i] * pb_x + g_0_xxyyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyyzz_0[i] = 2.0 * g_0_xyyyzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyyyzz_0[i] * pb_x + g_0_xxyyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyzzz_0[i] = 2.0 * g_0_xyyyzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyyzzz_0[i] * pb_x + g_0_xxyyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyzzzz_0[i] = 2.0 * g_0_xyyyzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yyzzzz_0[i] * pb_x + g_0_xxyyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yzzzzz_0[i] = 2.0 * g_0_xyyyzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_yzzzzz_0[i] * pb_x + g_0_xxyyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_zzzzzz_0[i] = 2.0 * g_0_xyyyzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_zzzzzz_0[i] * pb_x + g_0_xxyyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 504-532 components of targeted buffer : SLSI

    auto g_0_xxxyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 504);

    auto g_0_xxxyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 505);

    auto g_0_xxxyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 506);

    auto g_0_xxxyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 507);

    auto g_0_xxxyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 508);

    auto g_0_xxxyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 509);

    auto g_0_xxxyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 510);

    auto g_0_xxxyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 511);

    auto g_0_xxxyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 512);

    auto g_0_xxxyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 513);

    auto g_0_xxxyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 514);

    auto g_0_xxxyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 515);

    auto g_0_xxxyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 516);

    auto g_0_xxxyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 517);

    auto g_0_xxxyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 518);

    auto g_0_xxxyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 519);

    auto g_0_xxxyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 520);

    auto g_0_xxxyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 521);

    auto g_0_xxxyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 522);

    auto g_0_xxxyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 523);

    auto g_0_xxxyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 524);

    auto g_0_xxxyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 525);

    auto g_0_xxxyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 526);

    auto g_0_xxxyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 527);

    auto g_0_xxxyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 528);

    auto g_0_xxxyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 529);

    auto g_0_xxxyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 530);

    auto g_0_xxxyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 531);

#pragma omp simd aligned(g_0_xxxyyz_0_xxxxxy_0,       \
                             g_0_xxxyyz_0_xxxxxy_1,   \
                             g_0_xxxyyz_0_xxxxyy_0,   \
                             g_0_xxxyyz_0_xxxxyy_1,   \
                             g_0_xxxyyz_0_xxxyyy_0,   \
                             g_0_xxxyyz_0_xxxyyy_1,   \
                             g_0_xxxyyz_0_xxyyyy_0,   \
                             g_0_xxxyyz_0_xxyyyy_1,   \
                             g_0_xxxyyz_0_xyyyyy_0,   \
                             g_0_xxxyyz_0_xyyyyy_1,   \
                             g_0_xxxyyzz_0_xxxxxy_0,  \
                             g_0_xxxyyzz_0_xxxxxy_1,  \
                             g_0_xxxyyzz_0_xxxxyy_0,  \
                             g_0_xxxyyzz_0_xxxxyy_1,  \
                             g_0_xxxyyzz_0_xxxyyy_0,  \
                             g_0_xxxyyzz_0_xxxyyy_1,  \
                             g_0_xxxyyzz_0_xxyyyy_0,  \
                             g_0_xxxyyzz_0_xxyyyy_1,  \
                             g_0_xxxyyzz_0_xyyyyy_0,  \
                             g_0_xxxyyzz_0_xyyyyy_1,  \
                             g_0_xxxyyzzz_0_xxxxxx_0, \
                             g_0_xxxyyzzz_0_xxxxxy_0, \
                             g_0_xxxyyzzz_0_xxxxxz_0, \
                             g_0_xxxyyzzz_0_xxxxyy_0, \
                             g_0_xxxyyzzz_0_xxxxyz_0, \
                             g_0_xxxyyzzz_0_xxxxzz_0, \
                             g_0_xxxyyzzz_0_xxxyyy_0, \
                             g_0_xxxyyzzz_0_xxxyyz_0, \
                             g_0_xxxyyzzz_0_xxxyzz_0, \
                             g_0_xxxyyzzz_0_xxxzzz_0, \
                             g_0_xxxyyzzz_0_xxyyyy_0, \
                             g_0_xxxyyzzz_0_xxyyyz_0, \
                             g_0_xxxyyzzz_0_xxyyzz_0, \
                             g_0_xxxyyzzz_0_xxyzzz_0, \
                             g_0_xxxyyzzz_0_xxzzzz_0, \
                             g_0_xxxyyzzz_0_xyyyyy_0, \
                             g_0_xxxyyzzz_0_xyyyyz_0, \
                             g_0_xxxyyzzz_0_xyyyzz_0, \
                             g_0_xxxyyzzz_0_xyyzzz_0, \
                             g_0_xxxyyzzz_0_xyzzzz_0, \
                             g_0_xxxyyzzz_0_xzzzzz_0, \
                             g_0_xxxyyzzz_0_yyyyyy_0, \
                             g_0_xxxyyzzz_0_yyyyyz_0, \
                             g_0_xxxyyzzz_0_yyyyzz_0, \
                             g_0_xxxyyzzz_0_yyyzzz_0, \
                             g_0_xxxyyzzz_0_yyzzzz_0, \
                             g_0_xxxyyzzz_0_yzzzzz_0, \
                             g_0_xxxyyzzz_0_zzzzzz_0, \
                             g_0_xxxyzzz_0_xxxxxx_0,  \
                             g_0_xxxyzzz_0_xxxxxx_1,  \
                             g_0_xxxyzzz_0_xxxxxz_0,  \
                             g_0_xxxyzzz_0_xxxxxz_1,  \
                             g_0_xxxyzzz_0_xxxxzz_0,  \
                             g_0_xxxyzzz_0_xxxxzz_1,  \
                             g_0_xxxyzzz_0_xxxzzz_0,  \
                             g_0_xxxyzzz_0_xxxzzz_1,  \
                             g_0_xxxyzzz_0_xxzzzz_0,  \
                             g_0_xxxyzzz_0_xxzzzz_1,  \
                             g_0_xxxyzzz_0_xzzzzz_0,  \
                             g_0_xxxyzzz_0_xzzzzz_1,  \
                             g_0_xxxzzz_0_xxxxxx_0,   \
                             g_0_xxxzzz_0_xxxxxx_1,   \
                             g_0_xxxzzz_0_xxxxxz_0,   \
                             g_0_xxxzzz_0_xxxxxz_1,   \
                             g_0_xxxzzz_0_xxxxzz_0,   \
                             g_0_xxxzzz_0_xxxxzz_1,   \
                             g_0_xxxzzz_0_xxxzzz_0,   \
                             g_0_xxxzzz_0_xxxzzz_1,   \
                             g_0_xxxzzz_0_xxzzzz_0,   \
                             g_0_xxxzzz_0_xxzzzz_1,   \
                             g_0_xxxzzz_0_xzzzzz_0,   \
                             g_0_xxxzzz_0_xzzzzz_1,   \
                             g_0_xxyyzzz_0_xxxxyz_0,  \
                             g_0_xxyyzzz_0_xxxxyz_1,  \
                             g_0_xxyyzzz_0_xxxyyz_0,  \
                             g_0_xxyyzzz_0_xxxyyz_1,  \
                             g_0_xxyyzzz_0_xxxyz_1,   \
                             g_0_xxyyzzz_0_xxxyzz_0,  \
                             g_0_xxyyzzz_0_xxxyzz_1,  \
                             g_0_xxyyzzz_0_xxyyyz_0,  \
                             g_0_xxyyzzz_0_xxyyyz_1,  \
                             g_0_xxyyzzz_0_xxyyz_1,   \
                             g_0_xxyyzzz_0_xxyyzz_0,  \
                             g_0_xxyyzzz_0_xxyyzz_1,  \
                             g_0_xxyyzzz_0_xxyzz_1,   \
                             g_0_xxyyzzz_0_xxyzzz_0,  \
                             g_0_xxyyzzz_0_xxyzzz_1,  \
                             g_0_xxyyzzz_0_xyyyyz_0,  \
                             g_0_xxyyzzz_0_xyyyyz_1,  \
                             g_0_xxyyzzz_0_xyyyz_1,   \
                             g_0_xxyyzzz_0_xyyyzz_0,  \
                             g_0_xxyyzzz_0_xyyyzz_1,  \
                             g_0_xxyyzzz_0_xyyzz_1,   \
                             g_0_xxyyzzz_0_xyyzzz_0,  \
                             g_0_xxyyzzz_0_xyyzzz_1,  \
                             g_0_xxyyzzz_0_xyzzz_1,   \
                             g_0_xxyyzzz_0_xyzzzz_0,  \
                             g_0_xxyyzzz_0_xyzzzz_1,  \
                             g_0_xxyyzzz_0_yyyyyy_0,  \
                             g_0_xxyyzzz_0_yyyyyy_1,  \
                             g_0_xxyyzzz_0_yyyyyz_0,  \
                             g_0_xxyyzzz_0_yyyyyz_1,  \
                             g_0_xxyyzzz_0_yyyyz_1,   \
                             g_0_xxyyzzz_0_yyyyzz_0,  \
                             g_0_xxyyzzz_0_yyyyzz_1,  \
                             g_0_xxyyzzz_0_yyyzz_1,   \
                             g_0_xxyyzzz_0_yyyzzz_0,  \
                             g_0_xxyyzzz_0_yyyzzz_1,  \
                             g_0_xxyyzzz_0_yyzzz_1,   \
                             g_0_xxyyzzz_0_yyzzzz_0,  \
                             g_0_xxyyzzz_0_yyzzzz_1,  \
                             g_0_xxyyzzz_0_yzzzz_1,   \
                             g_0_xxyyzzz_0_yzzzzz_0,  \
                             g_0_xxyyzzz_0_yzzzzz_1,  \
                             g_0_xxyyzzz_0_zzzzzz_0,  \
                             g_0_xxyyzzz_0_zzzzzz_1,  \
                             g_0_xyyzzz_0_xxxxyz_0,   \
                             g_0_xyyzzz_0_xxxxyz_1,   \
                             g_0_xyyzzz_0_xxxyyz_0,   \
                             g_0_xyyzzz_0_xxxyyz_1,   \
                             g_0_xyyzzz_0_xxxyzz_0,   \
                             g_0_xyyzzz_0_xxxyzz_1,   \
                             g_0_xyyzzz_0_xxyyyz_0,   \
                             g_0_xyyzzz_0_xxyyyz_1,   \
                             g_0_xyyzzz_0_xxyyzz_0,   \
                             g_0_xyyzzz_0_xxyyzz_1,   \
                             g_0_xyyzzz_0_xxyzzz_0,   \
                             g_0_xyyzzz_0_xxyzzz_1,   \
                             g_0_xyyzzz_0_xyyyyz_0,   \
                             g_0_xyyzzz_0_xyyyyz_1,   \
                             g_0_xyyzzz_0_xyyyzz_0,   \
                             g_0_xyyzzz_0_xyyyzz_1,   \
                             g_0_xyyzzz_0_xyyzzz_0,   \
                             g_0_xyyzzz_0_xyyzzz_1,   \
                             g_0_xyyzzz_0_xyzzzz_0,   \
                             g_0_xyyzzz_0_xyzzzz_1,   \
                             g_0_xyyzzz_0_yyyyyy_0,   \
                             g_0_xyyzzz_0_yyyyyy_1,   \
                             g_0_xyyzzz_0_yyyyyz_0,   \
                             g_0_xyyzzz_0_yyyyyz_1,   \
                             g_0_xyyzzz_0_yyyyzz_0,   \
                             g_0_xyyzzz_0_yyyyzz_1,   \
                             g_0_xyyzzz_0_yyyzzz_0,   \
                             g_0_xyyzzz_0_yyyzzz_1,   \
                             g_0_xyyzzz_0_yyzzzz_0,   \
                             g_0_xyyzzz_0_yyzzzz_1,   \
                             g_0_xyyzzz_0_yzzzzz_0,   \
                             g_0_xyyzzz_0_yzzzzz_1,   \
                             g_0_xyyzzz_0_zzzzzz_0,   \
                             g_0_xyyzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_xxxxxx_0[i] = g_0_xxxzzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxxx_0[i] * pb_y +
                                     g_0_xxxyzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxxyyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxxxy_0[i] * pb_z + g_0_xxxyyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxxxz_0[i] = g_0_xxxzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxxz_0[i] * pb_y +
                                     g_0_xxxyzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxxyyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxxyy_0[i] * pb_z + g_0_xxxyyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxxyz_0[i] = 2.0 * g_0_xyyzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxyyzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxxzz_0[i] = g_0_xxxzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxzz_0[i] * pb_y +
                                     g_0_xxxyzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxxyyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxxyyy_0[i] * pb_z + g_0_xxxyyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxyyz_0[i] = 2.0 * g_0_xyyzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxyyzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxyzz_0[i] = 2.0 * g_0_xyyzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxyyzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxzzz_0[i] = g_0_xxxzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxzzz_0[i] * pb_y +
                                     g_0_xxxyzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxxyyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xxyyyy_0[i] * pb_z + g_0_xxxyyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxyyyz_0[i] = 2.0 * g_0_xyyzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxyyzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxyyzz_0[i] = 2.0 * g_0_xyyzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxyyzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxyzzz_0[i] = 2.0 * g_0_xyyzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxyyzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxzzzz_0[i] = g_0_xxxzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxzzzz_0[i] * pb_y +
                                     g_0_xxxyzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxxyyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxyyzz_0_xyyyyy_0[i] * pb_z + g_0_xxxyyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xyyyyz_0[i] = 2.0 * g_0_xyyzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyyyz_0[i] * pb_x + g_0_xxyyzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyyyzz_0[i] = 2.0 * g_0_xyyzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyyzz_0[i] * pb_x + g_0_xxyyzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyyzzz_0[i] = 2.0 * g_0_xyyzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyzzz_0[i] * pb_x + g_0_xxyyzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyzzzz_0[i] = 2.0 * g_0_xyyzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyzzzz_0[i] * pb_x + g_0_xxyyzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xzzzzz_0[i] = g_0_xxxzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xzzzzz_0[i] * pb_y +
                                     g_0_xxxyzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_yyyyyy_0[i] = 2.0 * g_0_xyyzzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyyyyy_0[i] * pb_x + g_0_xxyyzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyyyz_0[i] = 2.0 * g_0_xyyzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyyyyz_0[i] * pb_x + g_0_xxyyzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyyzz_0[i] = 2.0 * g_0_xyyzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyyyzz_0[i] * pb_x + g_0_xxyyzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyzzz_0[i] = 2.0 * g_0_xyyzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyyzzz_0[i] * pb_x + g_0_xxyyzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyzzzz_0[i] = 2.0 * g_0_xyyzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yyzzzz_0[i] * pb_x + g_0_xxyyzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yzzzzz_0[i] = 2.0 * g_0_xyyzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_yzzzzz_0[i] * pb_x + g_0_xxyyzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_zzzzzz_0[i] = 2.0 * g_0_xyyzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_zzzzzz_0[i] * pb_x + g_0_xxyyzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 532-560 components of targeted buffer : SLSI

    auto g_0_xxxyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 532);

    auto g_0_xxxyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 533);

    auto g_0_xxxyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 534);

    auto g_0_xxxyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 535);

    auto g_0_xxxyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 536);

    auto g_0_xxxyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 537);

    auto g_0_xxxyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 538);

    auto g_0_xxxyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 539);

    auto g_0_xxxyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 540);

    auto g_0_xxxyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 541);

    auto g_0_xxxyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 542);

    auto g_0_xxxyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 543);

    auto g_0_xxxyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 544);

    auto g_0_xxxyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 545);

    auto g_0_xxxyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 546);

    auto g_0_xxxyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 547);

    auto g_0_xxxyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 548);

    auto g_0_xxxyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 549);

    auto g_0_xxxyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 550);

    auto g_0_xxxyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 551);

    auto g_0_xxxyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 552);

    auto g_0_xxxyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 553);

    auto g_0_xxxyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 554);

    auto g_0_xxxyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 555);

    auto g_0_xxxyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 556);

    auto g_0_xxxyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 557);

    auto g_0_xxxyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 558);

    auto g_0_xxxyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 559);

#pragma omp simd aligned(g_0_xxxyzzzz_0_xxxxxx_0,     \
                             g_0_xxxyzzzz_0_xxxxxy_0, \
                             g_0_xxxyzzzz_0_xxxxxz_0, \
                             g_0_xxxyzzzz_0_xxxxyy_0, \
                             g_0_xxxyzzzz_0_xxxxyz_0, \
                             g_0_xxxyzzzz_0_xxxxzz_0, \
                             g_0_xxxyzzzz_0_xxxyyy_0, \
                             g_0_xxxyzzzz_0_xxxyyz_0, \
                             g_0_xxxyzzzz_0_xxxyzz_0, \
                             g_0_xxxyzzzz_0_xxxzzz_0, \
                             g_0_xxxyzzzz_0_xxyyyy_0, \
                             g_0_xxxyzzzz_0_xxyyyz_0, \
                             g_0_xxxyzzzz_0_xxyyzz_0, \
                             g_0_xxxyzzzz_0_xxyzzz_0, \
                             g_0_xxxyzzzz_0_xxzzzz_0, \
                             g_0_xxxyzzzz_0_xyyyyy_0, \
                             g_0_xxxyzzzz_0_xyyyyz_0, \
                             g_0_xxxyzzzz_0_xyyyzz_0, \
                             g_0_xxxyzzzz_0_xyyzzz_0, \
                             g_0_xxxyzzzz_0_xyzzzz_0, \
                             g_0_xxxyzzzz_0_xzzzzz_0, \
                             g_0_xxxyzzzz_0_yyyyyy_0, \
                             g_0_xxxyzzzz_0_yyyyyz_0, \
                             g_0_xxxyzzzz_0_yyyyzz_0, \
                             g_0_xxxyzzzz_0_yyyzzz_0, \
                             g_0_xxxyzzzz_0_yyzzzz_0, \
                             g_0_xxxyzzzz_0_yzzzzz_0, \
                             g_0_xxxyzzzz_0_zzzzzz_0, \
                             g_0_xxxzzzz_0_xxxxx_1,   \
                             g_0_xxxzzzz_0_xxxxxx_0,  \
                             g_0_xxxzzzz_0_xxxxxx_1,  \
                             g_0_xxxzzzz_0_xxxxxy_0,  \
                             g_0_xxxzzzz_0_xxxxxy_1,  \
                             g_0_xxxzzzz_0_xxxxxz_0,  \
                             g_0_xxxzzzz_0_xxxxxz_1,  \
                             g_0_xxxzzzz_0_xxxxy_1,   \
                             g_0_xxxzzzz_0_xxxxyy_0,  \
                             g_0_xxxzzzz_0_xxxxyy_1,  \
                             g_0_xxxzzzz_0_xxxxyz_0,  \
                             g_0_xxxzzzz_0_xxxxyz_1,  \
                             g_0_xxxzzzz_0_xxxxz_1,   \
                             g_0_xxxzzzz_0_xxxxzz_0,  \
                             g_0_xxxzzzz_0_xxxxzz_1,  \
                             g_0_xxxzzzz_0_xxxyy_1,   \
                             g_0_xxxzzzz_0_xxxyyy_0,  \
                             g_0_xxxzzzz_0_xxxyyy_1,  \
                             g_0_xxxzzzz_0_xxxyyz_0,  \
                             g_0_xxxzzzz_0_xxxyyz_1,  \
                             g_0_xxxzzzz_0_xxxyz_1,   \
                             g_0_xxxzzzz_0_xxxyzz_0,  \
                             g_0_xxxzzzz_0_xxxyzz_1,  \
                             g_0_xxxzzzz_0_xxxzz_1,   \
                             g_0_xxxzzzz_0_xxxzzz_0,  \
                             g_0_xxxzzzz_0_xxxzzz_1,  \
                             g_0_xxxzzzz_0_xxyyy_1,   \
                             g_0_xxxzzzz_0_xxyyyy_0,  \
                             g_0_xxxzzzz_0_xxyyyy_1,  \
                             g_0_xxxzzzz_0_xxyyyz_0,  \
                             g_0_xxxzzzz_0_xxyyyz_1,  \
                             g_0_xxxzzzz_0_xxyyz_1,   \
                             g_0_xxxzzzz_0_xxyyzz_0,  \
                             g_0_xxxzzzz_0_xxyyzz_1,  \
                             g_0_xxxzzzz_0_xxyzz_1,   \
                             g_0_xxxzzzz_0_xxyzzz_0,  \
                             g_0_xxxzzzz_0_xxyzzz_1,  \
                             g_0_xxxzzzz_0_xxzzz_1,   \
                             g_0_xxxzzzz_0_xxzzzz_0,  \
                             g_0_xxxzzzz_0_xxzzzz_1,  \
                             g_0_xxxzzzz_0_xyyyy_1,   \
                             g_0_xxxzzzz_0_xyyyyy_0,  \
                             g_0_xxxzzzz_0_xyyyyy_1,  \
                             g_0_xxxzzzz_0_xyyyyz_0,  \
                             g_0_xxxzzzz_0_xyyyyz_1,  \
                             g_0_xxxzzzz_0_xyyyz_1,   \
                             g_0_xxxzzzz_0_xyyyzz_0,  \
                             g_0_xxxzzzz_0_xyyyzz_1,  \
                             g_0_xxxzzzz_0_xyyzz_1,   \
                             g_0_xxxzzzz_0_xyyzzz_0,  \
                             g_0_xxxzzzz_0_xyyzzz_1,  \
                             g_0_xxxzzzz_0_xyzzz_1,   \
                             g_0_xxxzzzz_0_xyzzzz_0,  \
                             g_0_xxxzzzz_0_xyzzzz_1,  \
                             g_0_xxxzzzz_0_xzzzz_1,   \
                             g_0_xxxzzzz_0_xzzzzz_0,  \
                             g_0_xxxzzzz_0_xzzzzz_1,  \
                             g_0_xxxzzzz_0_yyyyy_1,   \
                             g_0_xxxzzzz_0_yyyyyy_0,  \
                             g_0_xxxzzzz_0_yyyyyy_1,  \
                             g_0_xxxzzzz_0_yyyyyz_0,  \
                             g_0_xxxzzzz_0_yyyyyz_1,  \
                             g_0_xxxzzzz_0_yyyyz_1,   \
                             g_0_xxxzzzz_0_yyyyzz_0,  \
                             g_0_xxxzzzz_0_yyyyzz_1,  \
                             g_0_xxxzzzz_0_yyyzz_1,   \
                             g_0_xxxzzzz_0_yyyzzz_0,  \
                             g_0_xxxzzzz_0_yyyzzz_1,  \
                             g_0_xxxzzzz_0_yyzzz_1,   \
                             g_0_xxxzzzz_0_yyzzzz_0,  \
                             g_0_xxxzzzz_0_yyzzzz_1,  \
                             g_0_xxxzzzz_0_yzzzz_1,   \
                             g_0_xxxzzzz_0_yzzzzz_0,  \
                             g_0_xxxzzzz_0_yzzzzz_1,  \
                             g_0_xxxzzzz_0_zzzzz_1,   \
                             g_0_xxxzzzz_0_zzzzzz_0,  \
                             g_0_xxxzzzz_0_zzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_xxxxxx_0[i] = g_0_xxxzzzz_0_xxxxxx_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxy_0[i] = g_0_xxxzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxxy_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxxz_0[i] = g_0_xxxzzzz_0_xxxxxz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxyy_0[i] =
            2.0 * g_0_xxxzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyy_0[i] * pb_y + g_0_xxxzzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxyz_0[i] = g_0_xxxzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxyz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxzz_0[i] = g_0_xxxzzzz_0_xxxxzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_xxxzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyy_0[i] * pb_y + g_0_xxxzzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyyz_0[i] =
            2.0 * g_0_xxxzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyyz_0[i] * pb_y + g_0_xxxzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyzz_0[i] = g_0_xxxzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxzzz_0[i] = g_0_xxxzzzz_0_xxxzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyyy_0[i] =
            4.0 * g_0_xxxzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyy_0[i] * pb_y + g_0_xxxzzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyyz_0[i] =
            3.0 * g_0_xxxzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyyz_0[i] * pb_y + g_0_xxxzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_xxxzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyzz_0[i] * pb_y + g_0_xxxzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyzzz_0[i] = g_0_xxxzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxzzzz_0[i] = g_0_xxxzzzz_0_xxzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyyy_0[i] =
            5.0 * g_0_xxxzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyy_0[i] * pb_y + g_0_xxxzzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyyz_0[i] =
            4.0 * g_0_xxxzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyyz_0[i] * pb_y + g_0_xxxzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyzz_0[i] =
            3.0 * g_0_xxxzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyzz_0[i] * pb_y + g_0_xxxzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyzzz_0[i] =
            2.0 * g_0_xxxzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyzzz_0[i] * pb_y + g_0_xxxzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyzzzz_0[i] = g_0_xxxzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xzzzzz_0[i] = g_0_xxxzzzz_0_xzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyyy_0[i] =
            6.0 * g_0_xxxzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyyy_0[i] * pb_y + g_0_xxxzzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyyz_0[i] =
            5.0 * g_0_xxxzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyyz_0[i] * pb_y + g_0_xxxzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyzz_0[i] =
            4.0 * g_0_xxxzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyzz_0[i] * pb_y + g_0_xxxzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyzzz_0[i] =
            3.0 * g_0_xxxzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyzzz_0[i] * pb_y + g_0_xxxzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyzzzz_0[i] =
            2.0 * g_0_xxxzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyzzzz_0[i] * pb_y + g_0_xxxzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yzzzzz_0[i] = g_0_xxxzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_zzzzzz_0[i] = g_0_xxxzzzz_0_zzzzzz_0[i] * pb_y + g_0_xxxzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 560-588 components of targeted buffer : SLSI

    auto g_0_xxxzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 560);

    auto g_0_xxxzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 561);

    auto g_0_xxxzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 562);

    auto g_0_xxxzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 563);

    auto g_0_xxxzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 564);

    auto g_0_xxxzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 565);

    auto g_0_xxxzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 566);

    auto g_0_xxxzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 567);

    auto g_0_xxxzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 568);

    auto g_0_xxxzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 569);

    auto g_0_xxxzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 570);

    auto g_0_xxxzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 571);

    auto g_0_xxxzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 572);

    auto g_0_xxxzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 573);

    auto g_0_xxxzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 574);

    auto g_0_xxxzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 575);

    auto g_0_xxxzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 576);

    auto g_0_xxxzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 577);

    auto g_0_xxxzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 578);

    auto g_0_xxxzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 579);

    auto g_0_xxxzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 580);

    auto g_0_xxxzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 581);

    auto g_0_xxxzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 582);

    auto g_0_xxxzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 583);

    auto g_0_xxxzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 584);

    auto g_0_xxxzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 585);

    auto g_0_xxxzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 586);

    auto g_0_xxxzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 587);

#pragma omp simd aligned(g_0_xxxzzz_0_xxxxxx_0,       \
                             g_0_xxxzzz_0_xxxxxx_1,   \
                             g_0_xxxzzz_0_xxxxxy_0,   \
                             g_0_xxxzzz_0_xxxxxy_1,   \
                             g_0_xxxzzz_0_xxxxyy_0,   \
                             g_0_xxxzzz_0_xxxxyy_1,   \
                             g_0_xxxzzz_0_xxxyyy_0,   \
                             g_0_xxxzzz_0_xxxyyy_1,   \
                             g_0_xxxzzz_0_xxyyyy_0,   \
                             g_0_xxxzzz_0_xxyyyy_1,   \
                             g_0_xxxzzz_0_xyyyyy_0,   \
                             g_0_xxxzzz_0_xyyyyy_1,   \
                             g_0_xxxzzzz_0_xxxxxx_0,  \
                             g_0_xxxzzzz_0_xxxxxx_1,  \
                             g_0_xxxzzzz_0_xxxxxy_0,  \
                             g_0_xxxzzzz_0_xxxxxy_1,  \
                             g_0_xxxzzzz_0_xxxxyy_0,  \
                             g_0_xxxzzzz_0_xxxxyy_1,  \
                             g_0_xxxzzzz_0_xxxyyy_0,  \
                             g_0_xxxzzzz_0_xxxyyy_1,  \
                             g_0_xxxzzzz_0_xxyyyy_0,  \
                             g_0_xxxzzzz_0_xxyyyy_1,  \
                             g_0_xxxzzzz_0_xyyyyy_0,  \
                             g_0_xxxzzzz_0_xyyyyy_1,  \
                             g_0_xxxzzzzz_0_xxxxxx_0, \
                             g_0_xxxzzzzz_0_xxxxxy_0, \
                             g_0_xxxzzzzz_0_xxxxxz_0, \
                             g_0_xxxzzzzz_0_xxxxyy_0, \
                             g_0_xxxzzzzz_0_xxxxyz_0, \
                             g_0_xxxzzzzz_0_xxxxzz_0, \
                             g_0_xxxzzzzz_0_xxxyyy_0, \
                             g_0_xxxzzzzz_0_xxxyyz_0, \
                             g_0_xxxzzzzz_0_xxxyzz_0, \
                             g_0_xxxzzzzz_0_xxxzzz_0, \
                             g_0_xxxzzzzz_0_xxyyyy_0, \
                             g_0_xxxzzzzz_0_xxyyyz_0, \
                             g_0_xxxzzzzz_0_xxyyzz_0, \
                             g_0_xxxzzzzz_0_xxyzzz_0, \
                             g_0_xxxzzzzz_0_xxzzzz_0, \
                             g_0_xxxzzzzz_0_xyyyyy_0, \
                             g_0_xxxzzzzz_0_xyyyyz_0, \
                             g_0_xxxzzzzz_0_xyyyzz_0, \
                             g_0_xxxzzzzz_0_xyyzzz_0, \
                             g_0_xxxzzzzz_0_xyzzzz_0, \
                             g_0_xxxzzzzz_0_xzzzzz_0, \
                             g_0_xxxzzzzz_0_yyyyyy_0, \
                             g_0_xxxzzzzz_0_yyyyyz_0, \
                             g_0_xxxzzzzz_0_yyyyzz_0, \
                             g_0_xxxzzzzz_0_yyyzzz_0, \
                             g_0_xxxzzzzz_0_yyzzzz_0, \
                             g_0_xxxzzzzz_0_yzzzzz_0, \
                             g_0_xxxzzzzz_0_zzzzzz_0, \
                             g_0_xxzzzzz_0_xxxxxz_0,  \
                             g_0_xxzzzzz_0_xxxxxz_1,  \
                             g_0_xxzzzzz_0_xxxxyz_0,  \
                             g_0_xxzzzzz_0_xxxxyz_1,  \
                             g_0_xxzzzzz_0_xxxxz_1,   \
                             g_0_xxzzzzz_0_xxxxzz_0,  \
                             g_0_xxzzzzz_0_xxxxzz_1,  \
                             g_0_xxzzzzz_0_xxxyyz_0,  \
                             g_0_xxzzzzz_0_xxxyyz_1,  \
                             g_0_xxzzzzz_0_xxxyz_1,   \
                             g_0_xxzzzzz_0_xxxyzz_0,  \
                             g_0_xxzzzzz_0_xxxyzz_1,  \
                             g_0_xxzzzzz_0_xxxzz_1,   \
                             g_0_xxzzzzz_0_xxxzzz_0,  \
                             g_0_xxzzzzz_0_xxxzzz_1,  \
                             g_0_xxzzzzz_0_xxyyyz_0,  \
                             g_0_xxzzzzz_0_xxyyyz_1,  \
                             g_0_xxzzzzz_0_xxyyz_1,   \
                             g_0_xxzzzzz_0_xxyyzz_0,  \
                             g_0_xxzzzzz_0_xxyyzz_1,  \
                             g_0_xxzzzzz_0_xxyzz_1,   \
                             g_0_xxzzzzz_0_xxyzzz_0,  \
                             g_0_xxzzzzz_0_xxyzzz_1,  \
                             g_0_xxzzzzz_0_xxzzz_1,   \
                             g_0_xxzzzzz_0_xxzzzz_0,  \
                             g_0_xxzzzzz_0_xxzzzz_1,  \
                             g_0_xxzzzzz_0_xyyyyz_0,  \
                             g_0_xxzzzzz_0_xyyyyz_1,  \
                             g_0_xxzzzzz_0_xyyyz_1,   \
                             g_0_xxzzzzz_0_xyyyzz_0,  \
                             g_0_xxzzzzz_0_xyyyzz_1,  \
                             g_0_xxzzzzz_0_xyyzz_1,   \
                             g_0_xxzzzzz_0_xyyzzz_0,  \
                             g_0_xxzzzzz_0_xyyzzz_1,  \
                             g_0_xxzzzzz_0_xyzzz_1,   \
                             g_0_xxzzzzz_0_xyzzzz_0,  \
                             g_0_xxzzzzz_0_xyzzzz_1,  \
                             g_0_xxzzzzz_0_xzzzz_1,   \
                             g_0_xxzzzzz_0_xzzzzz_0,  \
                             g_0_xxzzzzz_0_xzzzzz_1,  \
                             g_0_xxzzzzz_0_yyyyyy_0,  \
                             g_0_xxzzzzz_0_yyyyyy_1,  \
                             g_0_xxzzzzz_0_yyyyyz_0,  \
                             g_0_xxzzzzz_0_yyyyyz_1,  \
                             g_0_xxzzzzz_0_yyyyz_1,   \
                             g_0_xxzzzzz_0_yyyyzz_0,  \
                             g_0_xxzzzzz_0_yyyyzz_1,  \
                             g_0_xxzzzzz_0_yyyzz_1,   \
                             g_0_xxzzzzz_0_yyyzzz_0,  \
                             g_0_xxzzzzz_0_yyyzzz_1,  \
                             g_0_xxzzzzz_0_yyzzz_1,   \
                             g_0_xxzzzzz_0_yyzzzz_0,  \
                             g_0_xxzzzzz_0_yyzzzz_1,  \
                             g_0_xxzzzzz_0_yzzzz_1,   \
                             g_0_xxzzzzz_0_yzzzzz_0,  \
                             g_0_xxzzzzz_0_yzzzzz_1,  \
                             g_0_xxzzzzz_0_zzzzz_1,   \
                             g_0_xxzzzzz_0_zzzzzz_0,  \
                             g_0_xxzzzzz_0_zzzzzz_1,  \
                             g_0_xzzzzz_0_xxxxxz_0,   \
                             g_0_xzzzzz_0_xxxxxz_1,   \
                             g_0_xzzzzz_0_xxxxyz_0,   \
                             g_0_xzzzzz_0_xxxxyz_1,   \
                             g_0_xzzzzz_0_xxxxzz_0,   \
                             g_0_xzzzzz_0_xxxxzz_1,   \
                             g_0_xzzzzz_0_xxxyyz_0,   \
                             g_0_xzzzzz_0_xxxyyz_1,   \
                             g_0_xzzzzz_0_xxxyzz_0,   \
                             g_0_xzzzzz_0_xxxyzz_1,   \
                             g_0_xzzzzz_0_xxxzzz_0,   \
                             g_0_xzzzzz_0_xxxzzz_1,   \
                             g_0_xzzzzz_0_xxyyyz_0,   \
                             g_0_xzzzzz_0_xxyyyz_1,   \
                             g_0_xzzzzz_0_xxyyzz_0,   \
                             g_0_xzzzzz_0_xxyyzz_1,   \
                             g_0_xzzzzz_0_xxyzzz_0,   \
                             g_0_xzzzzz_0_xxyzzz_1,   \
                             g_0_xzzzzz_0_xxzzzz_0,   \
                             g_0_xzzzzz_0_xxzzzz_1,   \
                             g_0_xzzzzz_0_xyyyyz_0,   \
                             g_0_xzzzzz_0_xyyyyz_1,   \
                             g_0_xzzzzz_0_xyyyzz_0,   \
                             g_0_xzzzzz_0_xyyyzz_1,   \
                             g_0_xzzzzz_0_xyyzzz_0,   \
                             g_0_xzzzzz_0_xyyzzz_1,   \
                             g_0_xzzzzz_0_xyzzzz_0,   \
                             g_0_xzzzzz_0_xyzzzz_1,   \
                             g_0_xzzzzz_0_xzzzzz_0,   \
                             g_0_xzzzzz_0_xzzzzz_1,   \
                             g_0_xzzzzz_0_yyyyyy_0,   \
                             g_0_xzzzzz_0_yyyyyy_1,   \
                             g_0_xzzzzz_0_yyyyyz_0,   \
                             g_0_xzzzzz_0_yyyyyz_1,   \
                             g_0_xzzzzz_0_yyyyzz_0,   \
                             g_0_xzzzzz_0_yyyyzz_1,   \
                             g_0_xzzzzz_0_yyyzzz_0,   \
                             g_0_xzzzzz_0_yyyzzz_1,   \
                             g_0_xzzzzz_0_yyzzzz_0,   \
                             g_0_xzzzzz_0_yyzzzz_1,   \
                             g_0_xzzzzz_0_yzzzzz_0,   \
                             g_0_xzzzzz_0_yzzzzz_1,   \
                             g_0_xzzzzz_0_zzzzzz_0,   \
                             g_0_xzzzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_xxxxxx_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_xxxxxx_0[i] * pb_z + g_0_xxxzzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxxy_0[i] = 4.0 * g_0_xxxzzz_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_xxxxxy_0[i] * pb_z + g_0_xxxzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxxz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xxzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_xxxxyy_0[i] * pb_z + g_0_xxxzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxyz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxxzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xxzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxzz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_xxxyyy_0[i] * pb_z + g_0_xxxzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxyyz_0[i] = 2.0 * g_0_xzzzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxyzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xxzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxzzz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_xxyyyy_0[i] * pb_z + g_0_xxxzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxyyyz_0[i] = 2.0 * g_0_xzzzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_xzzzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xxzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxzzzz_0[i] * pb_x +
                                     g_0_xxzzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxxzzzz_0_xyyyyy_0[i] * pb_z + g_0_xxxzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xyyyyz_0[i] = 2.0 * g_0_xzzzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyz_0[i] * pb_x + g_0_xxzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyyzz_0[i] = 2.0 * g_0_xzzzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyzz_0[i] * pb_x + g_0_xxzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyzzz_0[i] = 2.0 * g_0_xzzzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyzzz_0[i] * pb_x + g_0_xxzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyyy_0[i] = 2.0 * g_0_xzzzzz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyyyyy_0[i] * pb_x + g_0_xxzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyyz_0[i] = 2.0 * g_0_xzzzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyyyyz_0[i] * pb_x + g_0_xxzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyzz_0[i] = 2.0 * g_0_xzzzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyyyzz_0[i] * pb_x + g_0_xxzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyzzz_0[i] = 2.0 * g_0_xzzzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyyzzz_0[i] * pb_x + g_0_xxzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyzzzz_0[i] = 2.0 * g_0_xzzzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yyzzzz_0[i] * pb_x + g_0_xxzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_yzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_zzzzzz_0[i] = 2.0 * g_0_xzzzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_zzzzzz_0[i] * pb_x + g_0_xxzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 588-616 components of targeted buffer : SLSI

    auto g_0_xxyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 588);

    auto g_0_xxyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 589);

    auto g_0_xxyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 590);

    auto g_0_xxyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 591);

    auto g_0_xxyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 592);

    auto g_0_xxyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 593);

    auto g_0_xxyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 594);

    auto g_0_xxyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 595);

    auto g_0_xxyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 596);

    auto g_0_xxyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 597);

    auto g_0_xxyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 598);

    auto g_0_xxyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 599);

    auto g_0_xxyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 600);

    auto g_0_xxyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 601);

    auto g_0_xxyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 602);

    auto g_0_xxyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 603);

    auto g_0_xxyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 604);

    auto g_0_xxyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 605);

    auto g_0_xxyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 606);

    auto g_0_xxyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 607);

    auto g_0_xxyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 608);

    auto g_0_xxyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 609);

    auto g_0_xxyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 610);

    auto g_0_xxyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 611);

    auto g_0_xxyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 612);

    auto g_0_xxyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 613);

    auto g_0_xxyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 614);

    auto g_0_xxyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 615);

#pragma omp simd aligned(g_0_xxyyyy_0_xxxxxx_0,       \
                             g_0_xxyyyy_0_xxxxxx_1,   \
                             g_0_xxyyyy_0_xxxxxz_0,   \
                             g_0_xxyyyy_0_xxxxxz_1,   \
                             g_0_xxyyyy_0_xxxxzz_0,   \
                             g_0_xxyyyy_0_xxxxzz_1,   \
                             g_0_xxyyyy_0_xxxzzz_0,   \
                             g_0_xxyyyy_0_xxxzzz_1,   \
                             g_0_xxyyyy_0_xxzzzz_0,   \
                             g_0_xxyyyy_0_xxzzzz_1,   \
                             g_0_xxyyyy_0_xzzzzz_0,   \
                             g_0_xxyyyy_0_xzzzzz_1,   \
                             g_0_xxyyyyy_0_xxxxxx_0,  \
                             g_0_xxyyyyy_0_xxxxxx_1,  \
                             g_0_xxyyyyy_0_xxxxxz_0,  \
                             g_0_xxyyyyy_0_xxxxxz_1,  \
                             g_0_xxyyyyy_0_xxxxzz_0,  \
                             g_0_xxyyyyy_0_xxxxzz_1,  \
                             g_0_xxyyyyy_0_xxxzzz_0,  \
                             g_0_xxyyyyy_0_xxxzzz_1,  \
                             g_0_xxyyyyy_0_xxzzzz_0,  \
                             g_0_xxyyyyy_0_xxzzzz_1,  \
                             g_0_xxyyyyy_0_xzzzzz_0,  \
                             g_0_xxyyyyy_0_xzzzzz_1,  \
                             g_0_xxyyyyyy_0_xxxxxx_0, \
                             g_0_xxyyyyyy_0_xxxxxy_0, \
                             g_0_xxyyyyyy_0_xxxxxz_0, \
                             g_0_xxyyyyyy_0_xxxxyy_0, \
                             g_0_xxyyyyyy_0_xxxxyz_0, \
                             g_0_xxyyyyyy_0_xxxxzz_0, \
                             g_0_xxyyyyyy_0_xxxyyy_0, \
                             g_0_xxyyyyyy_0_xxxyyz_0, \
                             g_0_xxyyyyyy_0_xxxyzz_0, \
                             g_0_xxyyyyyy_0_xxxzzz_0, \
                             g_0_xxyyyyyy_0_xxyyyy_0, \
                             g_0_xxyyyyyy_0_xxyyyz_0, \
                             g_0_xxyyyyyy_0_xxyyzz_0, \
                             g_0_xxyyyyyy_0_xxyzzz_0, \
                             g_0_xxyyyyyy_0_xxzzzz_0, \
                             g_0_xxyyyyyy_0_xyyyyy_0, \
                             g_0_xxyyyyyy_0_xyyyyz_0, \
                             g_0_xxyyyyyy_0_xyyyzz_0, \
                             g_0_xxyyyyyy_0_xyyzzz_0, \
                             g_0_xxyyyyyy_0_xyzzzz_0, \
                             g_0_xxyyyyyy_0_xzzzzz_0, \
                             g_0_xxyyyyyy_0_yyyyyy_0, \
                             g_0_xxyyyyyy_0_yyyyyz_0, \
                             g_0_xxyyyyyy_0_yyyyzz_0, \
                             g_0_xxyyyyyy_0_yyyzzz_0, \
                             g_0_xxyyyyyy_0_yyzzzz_0, \
                             g_0_xxyyyyyy_0_yzzzzz_0, \
                             g_0_xxyyyyyy_0_zzzzzz_0, \
                             g_0_xyyyyyy_0_xxxxxy_0,  \
                             g_0_xyyyyyy_0_xxxxxy_1,  \
                             g_0_xyyyyyy_0_xxxxy_1,   \
                             g_0_xyyyyyy_0_xxxxyy_0,  \
                             g_0_xyyyyyy_0_xxxxyy_1,  \
                             g_0_xyyyyyy_0_xxxxyz_0,  \
                             g_0_xyyyyyy_0_xxxxyz_1,  \
                             g_0_xyyyyyy_0_xxxyy_1,   \
                             g_0_xyyyyyy_0_xxxyyy_0,  \
                             g_0_xyyyyyy_0_xxxyyy_1,  \
                             g_0_xyyyyyy_0_xxxyyz_0,  \
                             g_0_xyyyyyy_0_xxxyyz_1,  \
                             g_0_xyyyyyy_0_xxxyz_1,   \
                             g_0_xyyyyyy_0_xxxyzz_0,  \
                             g_0_xyyyyyy_0_xxxyzz_1,  \
                             g_0_xyyyyyy_0_xxyyy_1,   \
                             g_0_xyyyyyy_0_xxyyyy_0,  \
                             g_0_xyyyyyy_0_xxyyyy_1,  \
                             g_0_xyyyyyy_0_xxyyyz_0,  \
                             g_0_xyyyyyy_0_xxyyyz_1,  \
                             g_0_xyyyyyy_0_xxyyz_1,   \
                             g_0_xyyyyyy_0_xxyyzz_0,  \
                             g_0_xyyyyyy_0_xxyyzz_1,  \
                             g_0_xyyyyyy_0_xxyzz_1,   \
                             g_0_xyyyyyy_0_xxyzzz_0,  \
                             g_0_xyyyyyy_0_xxyzzz_1,  \
                             g_0_xyyyyyy_0_xyyyy_1,   \
                             g_0_xyyyyyy_0_xyyyyy_0,  \
                             g_0_xyyyyyy_0_xyyyyy_1,  \
                             g_0_xyyyyyy_0_xyyyyz_0,  \
                             g_0_xyyyyyy_0_xyyyyz_1,  \
                             g_0_xyyyyyy_0_xyyyz_1,   \
                             g_0_xyyyyyy_0_xyyyzz_0,  \
                             g_0_xyyyyyy_0_xyyyzz_1,  \
                             g_0_xyyyyyy_0_xyyzz_1,   \
                             g_0_xyyyyyy_0_xyyzzz_0,  \
                             g_0_xyyyyyy_0_xyyzzz_1,  \
                             g_0_xyyyyyy_0_xyzzz_1,   \
                             g_0_xyyyyyy_0_xyzzzz_0,  \
                             g_0_xyyyyyy_0_xyzzzz_1,  \
                             g_0_xyyyyyy_0_yyyyy_1,   \
                             g_0_xyyyyyy_0_yyyyyy_0,  \
                             g_0_xyyyyyy_0_yyyyyy_1,  \
                             g_0_xyyyyyy_0_yyyyyz_0,  \
                             g_0_xyyyyyy_0_yyyyyz_1,  \
                             g_0_xyyyyyy_0_yyyyz_1,   \
                             g_0_xyyyyyy_0_yyyyzz_0,  \
                             g_0_xyyyyyy_0_yyyyzz_1,  \
                             g_0_xyyyyyy_0_yyyzz_1,   \
                             g_0_xyyyyyy_0_yyyzzz_0,  \
                             g_0_xyyyyyy_0_yyyzzz_1,  \
                             g_0_xyyyyyy_0_yyzzz_1,   \
                             g_0_xyyyyyy_0_yyzzzz_0,  \
                             g_0_xyyyyyy_0_yyzzzz_1,  \
                             g_0_xyyyyyy_0_yzzzz_1,   \
                             g_0_xyyyyyy_0_yzzzzz_0,  \
                             g_0_xyyyyyy_0_yzzzzz_1,  \
                             g_0_xyyyyyy_0_zzzzzz_0,  \
                             g_0_xyyyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyy_0_xxxxxy_0,   \
                             g_0_yyyyyy_0_xxxxxy_1,   \
                             g_0_yyyyyy_0_xxxxyy_0,   \
                             g_0_yyyyyy_0_xxxxyy_1,   \
                             g_0_yyyyyy_0_xxxxyz_0,   \
                             g_0_yyyyyy_0_xxxxyz_1,   \
                             g_0_yyyyyy_0_xxxyyy_0,   \
                             g_0_yyyyyy_0_xxxyyy_1,   \
                             g_0_yyyyyy_0_xxxyyz_0,   \
                             g_0_yyyyyy_0_xxxyyz_1,   \
                             g_0_yyyyyy_0_xxxyzz_0,   \
                             g_0_yyyyyy_0_xxxyzz_1,   \
                             g_0_yyyyyy_0_xxyyyy_0,   \
                             g_0_yyyyyy_0_xxyyyy_1,   \
                             g_0_yyyyyy_0_xxyyyz_0,   \
                             g_0_yyyyyy_0_xxyyyz_1,   \
                             g_0_yyyyyy_0_xxyyzz_0,   \
                             g_0_yyyyyy_0_xxyyzz_1,   \
                             g_0_yyyyyy_0_xxyzzz_0,   \
                             g_0_yyyyyy_0_xxyzzz_1,   \
                             g_0_yyyyyy_0_xyyyyy_0,   \
                             g_0_yyyyyy_0_xyyyyy_1,   \
                             g_0_yyyyyy_0_xyyyyz_0,   \
                             g_0_yyyyyy_0_xyyyyz_1,   \
                             g_0_yyyyyy_0_xyyyzz_0,   \
                             g_0_yyyyyy_0_xyyyzz_1,   \
                             g_0_yyyyyy_0_xyyzzz_0,   \
                             g_0_yyyyyy_0_xyyzzz_1,   \
                             g_0_yyyyyy_0_xyzzzz_0,   \
                             g_0_yyyyyy_0_xyzzzz_1,   \
                             g_0_yyyyyy_0_yyyyyy_0,   \
                             g_0_yyyyyy_0_yyyyyy_1,   \
                             g_0_yyyyyy_0_yyyyyz_0,   \
                             g_0_yyyyyy_0_yyyyyz_1,   \
                             g_0_yyyyyy_0_yyyyzz_0,   \
                             g_0_yyyyyy_0_yyyyzz_1,   \
                             g_0_yyyyyy_0_yyyzzz_0,   \
                             g_0_yyyyyy_0_yyyzzz_1,   \
                             g_0_yyyyyy_0_yyzzzz_0,   \
                             g_0_yyyyyy_0_yyzzzz_1,   \
                             g_0_yyyyyy_0_yzzzzz_0,   \
                             g_0_yyyyyy_0_yzzzzz_1,   \
                             g_0_yyyyyy_0_zzzzzz_0,   \
                             g_0_yyyyyy_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_xxxxxx_0[i] = 5.0 * g_0_xxyyyy_0_xxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_xxxxxx_0[i] * pb_y + g_0_xxyyyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxxxy_0[i] = g_0_yyyyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                     5.0 * g_0_xyyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxxy_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxxz_0[i] = 5.0 * g_0_xxyyyy_0_xxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_xxxxxz_0[i] * pb_y + g_0_xxyyyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxxyy_0[i] = g_0_yyyyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxyy_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxyz_0[i] = g_0_yyyyyy_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxyz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxzz_0[i] = 5.0 * g_0_xxyyyy_0_xxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_xxxxzz_0[i] * pb_y + g_0_xxyyyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxyyy_0[i] = g_0_yyyyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyyy_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxyyz_0[i] = g_0_yyyyyy_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyyz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxyzz_0[i] = g_0_yyyyyy_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_xxxzzz_0[i] * pb_y + g_0_xxyyyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxyyyy_0[i] = g_0_yyyyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyyy_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyyyz_0[i] = g_0_yyyyyy_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyyz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyyzz_0[i] = g_0_yyyyyy_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyzzz_0[i] = g_0_yyyyyy_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyzzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_xxzzzz_0[i] * pb_y + g_0_xxyyyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xyyyyy_0[i] = g_0_yyyyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyy_1[i] * fi_abcd_0 +
                                     g_0_xyyyyyy_0_xyyyyy_0[i] * pb_x + g_0_xyyyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyyyz_0[i] = g_0_yyyyyy_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyyy_0_xyyyyz_0[i] * pb_x + g_0_xyyyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyyzz_0[i] = g_0_yyyyyy_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyyy_0_xyyyzz_0[i] * pb_x + g_0_xyyyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyzzz_0[i] = g_0_yyyyyy_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyyy_0_xyyzzz_0[i] * pb_x + g_0_xyyyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyzzzz_0[i] = g_0_yyyyyy_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyyy_0_xyzzzz_0[i] * pb_x + g_0_xyyyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xzzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyyy_0_xzzzzz_0[i] * pb_y + g_0_xxyyyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_yyyyyy_0[i] = g_0_yyyyyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyy_0[i] * pb_x +
                                     g_0_xyyyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyyyz_0[i] = g_0_yyyyyy_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyyz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyyzz_0[i] = g_0_yyyyyy_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyzzz_0[i] = g_0_yyyyyy_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyzzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyzzzz_0[i] = g_0_yyyyyy_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzzzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yzzzzz_0[i] = g_0_yyyyyy_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzzzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_zzzzzz_0[i] = g_0_yyyyyy_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_zzzzzz_0[i] * pb_x +
                                     g_0_xyyyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 616-644 components of targeted buffer : SLSI

    auto g_0_xxyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 616);

    auto g_0_xxyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 617);

    auto g_0_xxyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 618);

    auto g_0_xxyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 619);

    auto g_0_xxyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 620);

    auto g_0_xxyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 621);

    auto g_0_xxyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 622);

    auto g_0_xxyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 623);

    auto g_0_xxyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 624);

    auto g_0_xxyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 625);

    auto g_0_xxyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 626);

    auto g_0_xxyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 627);

    auto g_0_xxyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 628);

    auto g_0_xxyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 629);

    auto g_0_xxyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 630);

    auto g_0_xxyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 631);

    auto g_0_xxyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 632);

    auto g_0_xxyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 633);

    auto g_0_xxyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 634);

    auto g_0_xxyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 635);

    auto g_0_xxyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 636);

    auto g_0_xxyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 637);

    auto g_0_xxyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 638);

    auto g_0_xxyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 639);

    auto g_0_xxyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 640);

    auto g_0_xxyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 641);

    auto g_0_xxyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 642);

    auto g_0_xxyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 643);

#pragma omp simd aligned(g_0_xxyyyyy_0_xxxxx_1,       \
                             g_0_xxyyyyy_0_xxxxxx_0,  \
                             g_0_xxyyyyy_0_xxxxxx_1,  \
                             g_0_xxyyyyy_0_xxxxxy_0,  \
                             g_0_xxyyyyy_0_xxxxxy_1,  \
                             g_0_xxyyyyy_0_xxxxxz_0,  \
                             g_0_xxyyyyy_0_xxxxxz_1,  \
                             g_0_xxyyyyy_0_xxxxy_1,   \
                             g_0_xxyyyyy_0_xxxxyy_0,  \
                             g_0_xxyyyyy_0_xxxxyy_1,  \
                             g_0_xxyyyyy_0_xxxxyz_0,  \
                             g_0_xxyyyyy_0_xxxxyz_1,  \
                             g_0_xxyyyyy_0_xxxxz_1,   \
                             g_0_xxyyyyy_0_xxxxzz_0,  \
                             g_0_xxyyyyy_0_xxxxzz_1,  \
                             g_0_xxyyyyy_0_xxxyy_1,   \
                             g_0_xxyyyyy_0_xxxyyy_0,  \
                             g_0_xxyyyyy_0_xxxyyy_1,  \
                             g_0_xxyyyyy_0_xxxyyz_0,  \
                             g_0_xxyyyyy_0_xxxyyz_1,  \
                             g_0_xxyyyyy_0_xxxyz_1,   \
                             g_0_xxyyyyy_0_xxxyzz_0,  \
                             g_0_xxyyyyy_0_xxxyzz_1,  \
                             g_0_xxyyyyy_0_xxxzz_1,   \
                             g_0_xxyyyyy_0_xxxzzz_0,  \
                             g_0_xxyyyyy_0_xxxzzz_1,  \
                             g_0_xxyyyyy_0_xxyyy_1,   \
                             g_0_xxyyyyy_0_xxyyyy_0,  \
                             g_0_xxyyyyy_0_xxyyyy_1,  \
                             g_0_xxyyyyy_0_xxyyyz_0,  \
                             g_0_xxyyyyy_0_xxyyyz_1,  \
                             g_0_xxyyyyy_0_xxyyz_1,   \
                             g_0_xxyyyyy_0_xxyyzz_0,  \
                             g_0_xxyyyyy_0_xxyyzz_1,  \
                             g_0_xxyyyyy_0_xxyzz_1,   \
                             g_0_xxyyyyy_0_xxyzzz_0,  \
                             g_0_xxyyyyy_0_xxyzzz_1,  \
                             g_0_xxyyyyy_0_xxzzz_1,   \
                             g_0_xxyyyyy_0_xxzzzz_0,  \
                             g_0_xxyyyyy_0_xxzzzz_1,  \
                             g_0_xxyyyyy_0_xyyyy_1,   \
                             g_0_xxyyyyy_0_xyyyyy_0,  \
                             g_0_xxyyyyy_0_xyyyyy_1,  \
                             g_0_xxyyyyy_0_xyyyyz_0,  \
                             g_0_xxyyyyy_0_xyyyyz_1,  \
                             g_0_xxyyyyy_0_xyyyz_1,   \
                             g_0_xxyyyyy_0_xyyyzz_0,  \
                             g_0_xxyyyyy_0_xyyyzz_1,  \
                             g_0_xxyyyyy_0_xyyzz_1,   \
                             g_0_xxyyyyy_0_xyyzzz_0,  \
                             g_0_xxyyyyy_0_xyyzzz_1,  \
                             g_0_xxyyyyy_0_xyzzz_1,   \
                             g_0_xxyyyyy_0_xyzzzz_0,  \
                             g_0_xxyyyyy_0_xyzzzz_1,  \
                             g_0_xxyyyyy_0_xzzzz_1,   \
                             g_0_xxyyyyy_0_xzzzzz_0,  \
                             g_0_xxyyyyy_0_xzzzzz_1,  \
                             g_0_xxyyyyy_0_yyyyy_1,   \
                             g_0_xxyyyyy_0_yyyyyy_0,  \
                             g_0_xxyyyyy_0_yyyyyy_1,  \
                             g_0_xxyyyyy_0_yyyyyz_0,  \
                             g_0_xxyyyyy_0_yyyyyz_1,  \
                             g_0_xxyyyyy_0_yyyyz_1,   \
                             g_0_xxyyyyy_0_yyyyzz_0,  \
                             g_0_xxyyyyy_0_yyyyzz_1,  \
                             g_0_xxyyyyy_0_yyyzz_1,   \
                             g_0_xxyyyyy_0_yyyzzz_0,  \
                             g_0_xxyyyyy_0_yyyzzz_1,  \
                             g_0_xxyyyyy_0_yyzzz_1,   \
                             g_0_xxyyyyy_0_yyzzzz_0,  \
                             g_0_xxyyyyy_0_yyzzzz_1,  \
                             g_0_xxyyyyy_0_yzzzz_1,   \
                             g_0_xxyyyyy_0_yzzzzz_0,  \
                             g_0_xxyyyyy_0_yzzzzz_1,  \
                             g_0_xxyyyyy_0_zzzzz_1,   \
                             g_0_xxyyyyy_0_zzzzzz_0,  \
                             g_0_xxyyyyy_0_zzzzzz_1,  \
                             g_0_xxyyyyyz_0_xxxxxx_0, \
                             g_0_xxyyyyyz_0_xxxxxy_0, \
                             g_0_xxyyyyyz_0_xxxxxz_0, \
                             g_0_xxyyyyyz_0_xxxxyy_0, \
                             g_0_xxyyyyyz_0_xxxxyz_0, \
                             g_0_xxyyyyyz_0_xxxxzz_0, \
                             g_0_xxyyyyyz_0_xxxyyy_0, \
                             g_0_xxyyyyyz_0_xxxyyz_0, \
                             g_0_xxyyyyyz_0_xxxyzz_0, \
                             g_0_xxyyyyyz_0_xxxzzz_0, \
                             g_0_xxyyyyyz_0_xxyyyy_0, \
                             g_0_xxyyyyyz_0_xxyyyz_0, \
                             g_0_xxyyyyyz_0_xxyyzz_0, \
                             g_0_xxyyyyyz_0_xxyzzz_0, \
                             g_0_xxyyyyyz_0_xxzzzz_0, \
                             g_0_xxyyyyyz_0_xyyyyy_0, \
                             g_0_xxyyyyyz_0_xyyyyz_0, \
                             g_0_xxyyyyyz_0_xyyyzz_0, \
                             g_0_xxyyyyyz_0_xyyzzz_0, \
                             g_0_xxyyyyyz_0_xyzzzz_0, \
                             g_0_xxyyyyyz_0_xzzzzz_0, \
                             g_0_xxyyyyyz_0_yyyyyy_0, \
                             g_0_xxyyyyyz_0_yyyyyz_0, \
                             g_0_xxyyyyyz_0_yyyyzz_0, \
                             g_0_xxyyyyyz_0_yyyzzz_0, \
                             g_0_xxyyyyyz_0_yyzzzz_0, \
                             g_0_xxyyyyyz_0_yzzzzz_0, \
                             g_0_xxyyyyyz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_xxxxxx_0[i] = g_0_xxyyyyy_0_xxxxxx_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxy_0[i] = g_0_xxyyyyy_0_xxxxxy_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxxz_0[i] = g_0_xxyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxxz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxyy_0[i] = g_0_xxyyyyy_0_xxxxyy_0[i] * pb_z + g_0_xxyyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxyz_0[i] = g_0_xxyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxyz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxzz_0[i] =
            2.0 * g_0_xxyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyyy_0[i] = g_0_xxyyyyy_0_xxxyyy_0[i] * pb_z + g_0_xxyyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyyz_0[i] = g_0_xxyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyyz_0[i] * pb_z + g_0_xxyyyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyzz_0[i] =
            2.0 * g_0_xxyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxzzz_0[i] =
            3.0 * g_0_xxyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyyy_0[i] = g_0_xxyyyyy_0_xxyyyy_0[i] * pb_z + g_0_xxyyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyyz_0[i] = g_0_xxyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyyz_0[i] * pb_z + g_0_xxyyyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyzz_0[i] =
            2.0 * g_0_xxyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyzz_0[i] * pb_z + g_0_xxyyyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyzzz_0[i] =
            3.0 * g_0_xxyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxzzzz_0[i] =
            4.0 * g_0_xxyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyyy_0[i] = g_0_xxyyyyy_0_xyyyyy_0[i] * pb_z + g_0_xxyyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyyz_0[i] = g_0_xxyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyyz_0[i] * pb_z + g_0_xxyyyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyzz_0[i] =
            2.0 * g_0_xxyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyzz_0[i] * pb_z + g_0_xxyyyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyzzz_0[i] =
            3.0 * g_0_xxyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyzzz_0[i] * pb_z + g_0_xxyyyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyzzzz_0[i] =
            4.0 * g_0_xxyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xzzzzz_0[i] =
            5.0 * g_0_xxyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyyy_0[i] = g_0_xxyyyyy_0_yyyyyy_0[i] * pb_z + g_0_xxyyyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyyz_0[i] = g_0_xxyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyyyz_0[i] * pb_z + g_0_xxyyyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyzz_0[i] =
            2.0 * g_0_xxyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyyzz_0[i] * pb_z + g_0_xxyyyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyzzz_0[i] =
            3.0 * g_0_xxyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyzzz_0[i] * pb_z + g_0_xxyyyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyzzzz_0[i] =
            4.0 * g_0_xxyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyzzzz_0[i] * pb_z + g_0_xxyyyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yzzzzz_0[i] =
            5.0 * g_0_xxyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_zzzzzz_0[i] =
            6.0 * g_0_xxyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_zzzzzz_0[i] * pb_z + g_0_xxyyyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 644-672 components of targeted buffer : SLSI

    auto g_0_xxyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 644);

    auto g_0_xxyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 645);

    auto g_0_xxyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 646);

    auto g_0_xxyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 647);

    auto g_0_xxyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 648);

    auto g_0_xxyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 649);

    auto g_0_xxyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 650);

    auto g_0_xxyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 651);

    auto g_0_xxyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 652);

    auto g_0_xxyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 653);

    auto g_0_xxyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 654);

    auto g_0_xxyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 655);

    auto g_0_xxyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 656);

    auto g_0_xxyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 657);

    auto g_0_xxyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 658);

    auto g_0_xxyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 659);

    auto g_0_xxyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 660);

    auto g_0_xxyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 661);

    auto g_0_xxyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 662);

    auto g_0_xxyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 663);

    auto g_0_xxyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 664);

    auto g_0_xxyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 665);

    auto g_0_xxyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 666);

    auto g_0_xxyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 667);

    auto g_0_xxyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 668);

    auto g_0_xxyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 669);

    auto g_0_xxyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 670);

    auto g_0_xxyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 671);

#pragma omp simd aligned(g_0_xxyyyy_0_xxxxxy_0,       \
                             g_0_xxyyyy_0_xxxxxy_1,   \
                             g_0_xxyyyy_0_xxxxyy_0,   \
                             g_0_xxyyyy_0_xxxxyy_1,   \
                             g_0_xxyyyy_0_xxxyyy_0,   \
                             g_0_xxyyyy_0_xxxyyy_1,   \
                             g_0_xxyyyy_0_xxyyyy_0,   \
                             g_0_xxyyyy_0_xxyyyy_1,   \
                             g_0_xxyyyy_0_xyyyyy_0,   \
                             g_0_xxyyyy_0_xyyyyy_1,   \
                             g_0_xxyyyyz_0_xxxxxy_0,  \
                             g_0_xxyyyyz_0_xxxxxy_1,  \
                             g_0_xxyyyyz_0_xxxxyy_0,  \
                             g_0_xxyyyyz_0_xxxxyy_1,  \
                             g_0_xxyyyyz_0_xxxyyy_0,  \
                             g_0_xxyyyyz_0_xxxyyy_1,  \
                             g_0_xxyyyyz_0_xxyyyy_0,  \
                             g_0_xxyyyyz_0_xxyyyy_1,  \
                             g_0_xxyyyyz_0_xyyyyy_0,  \
                             g_0_xxyyyyz_0_xyyyyy_1,  \
                             g_0_xxyyyyzz_0_xxxxxx_0, \
                             g_0_xxyyyyzz_0_xxxxxy_0, \
                             g_0_xxyyyyzz_0_xxxxxz_0, \
                             g_0_xxyyyyzz_0_xxxxyy_0, \
                             g_0_xxyyyyzz_0_xxxxyz_0, \
                             g_0_xxyyyyzz_0_xxxxzz_0, \
                             g_0_xxyyyyzz_0_xxxyyy_0, \
                             g_0_xxyyyyzz_0_xxxyyz_0, \
                             g_0_xxyyyyzz_0_xxxyzz_0, \
                             g_0_xxyyyyzz_0_xxxzzz_0, \
                             g_0_xxyyyyzz_0_xxyyyy_0, \
                             g_0_xxyyyyzz_0_xxyyyz_0, \
                             g_0_xxyyyyzz_0_xxyyzz_0, \
                             g_0_xxyyyyzz_0_xxyzzz_0, \
                             g_0_xxyyyyzz_0_xxzzzz_0, \
                             g_0_xxyyyyzz_0_xyyyyy_0, \
                             g_0_xxyyyyzz_0_xyyyyz_0, \
                             g_0_xxyyyyzz_0_xyyyzz_0, \
                             g_0_xxyyyyzz_0_xyyzzz_0, \
                             g_0_xxyyyyzz_0_xyzzzz_0, \
                             g_0_xxyyyyzz_0_xzzzzz_0, \
                             g_0_xxyyyyzz_0_yyyyyy_0, \
                             g_0_xxyyyyzz_0_yyyyyz_0, \
                             g_0_xxyyyyzz_0_yyyyzz_0, \
                             g_0_xxyyyyzz_0_yyyzzz_0, \
                             g_0_xxyyyyzz_0_yyzzzz_0, \
                             g_0_xxyyyyzz_0_yzzzzz_0, \
                             g_0_xxyyyyzz_0_zzzzzz_0, \
                             g_0_xxyyyzz_0_xxxxxx_0,  \
                             g_0_xxyyyzz_0_xxxxxx_1,  \
                             g_0_xxyyyzz_0_xxxxxz_0,  \
                             g_0_xxyyyzz_0_xxxxxz_1,  \
                             g_0_xxyyyzz_0_xxxxzz_0,  \
                             g_0_xxyyyzz_0_xxxxzz_1,  \
                             g_0_xxyyyzz_0_xxxzzz_0,  \
                             g_0_xxyyyzz_0_xxxzzz_1,  \
                             g_0_xxyyyzz_0_xxzzzz_0,  \
                             g_0_xxyyyzz_0_xxzzzz_1,  \
                             g_0_xxyyyzz_0_xzzzzz_0,  \
                             g_0_xxyyyzz_0_xzzzzz_1,  \
                             g_0_xxyyzz_0_xxxxxx_0,   \
                             g_0_xxyyzz_0_xxxxxx_1,   \
                             g_0_xxyyzz_0_xxxxxz_0,   \
                             g_0_xxyyzz_0_xxxxxz_1,   \
                             g_0_xxyyzz_0_xxxxzz_0,   \
                             g_0_xxyyzz_0_xxxxzz_1,   \
                             g_0_xxyyzz_0_xxxzzz_0,   \
                             g_0_xxyyzz_0_xxxzzz_1,   \
                             g_0_xxyyzz_0_xxzzzz_0,   \
                             g_0_xxyyzz_0_xxzzzz_1,   \
                             g_0_xxyyzz_0_xzzzzz_0,   \
                             g_0_xxyyzz_0_xzzzzz_1,   \
                             g_0_xyyyyzz_0_xxxxyz_0,  \
                             g_0_xyyyyzz_0_xxxxyz_1,  \
                             g_0_xyyyyzz_0_xxxyyz_0,  \
                             g_0_xyyyyzz_0_xxxyyz_1,  \
                             g_0_xyyyyzz_0_xxxyz_1,   \
                             g_0_xyyyyzz_0_xxxyzz_0,  \
                             g_0_xyyyyzz_0_xxxyzz_1,  \
                             g_0_xyyyyzz_0_xxyyyz_0,  \
                             g_0_xyyyyzz_0_xxyyyz_1,  \
                             g_0_xyyyyzz_0_xxyyz_1,   \
                             g_0_xyyyyzz_0_xxyyzz_0,  \
                             g_0_xyyyyzz_0_xxyyzz_1,  \
                             g_0_xyyyyzz_0_xxyzz_1,   \
                             g_0_xyyyyzz_0_xxyzzz_0,  \
                             g_0_xyyyyzz_0_xxyzzz_1,  \
                             g_0_xyyyyzz_0_xyyyyz_0,  \
                             g_0_xyyyyzz_0_xyyyyz_1,  \
                             g_0_xyyyyzz_0_xyyyz_1,   \
                             g_0_xyyyyzz_0_xyyyzz_0,  \
                             g_0_xyyyyzz_0_xyyyzz_1,  \
                             g_0_xyyyyzz_0_xyyzz_1,   \
                             g_0_xyyyyzz_0_xyyzzz_0,  \
                             g_0_xyyyyzz_0_xyyzzz_1,  \
                             g_0_xyyyyzz_0_xyzzz_1,   \
                             g_0_xyyyyzz_0_xyzzzz_0,  \
                             g_0_xyyyyzz_0_xyzzzz_1,  \
                             g_0_xyyyyzz_0_yyyyyy_0,  \
                             g_0_xyyyyzz_0_yyyyyy_1,  \
                             g_0_xyyyyzz_0_yyyyyz_0,  \
                             g_0_xyyyyzz_0_yyyyyz_1,  \
                             g_0_xyyyyzz_0_yyyyz_1,   \
                             g_0_xyyyyzz_0_yyyyzz_0,  \
                             g_0_xyyyyzz_0_yyyyzz_1,  \
                             g_0_xyyyyzz_0_yyyzz_1,   \
                             g_0_xyyyyzz_0_yyyzzz_0,  \
                             g_0_xyyyyzz_0_yyyzzz_1,  \
                             g_0_xyyyyzz_0_yyzzz_1,   \
                             g_0_xyyyyzz_0_yyzzzz_0,  \
                             g_0_xyyyyzz_0_yyzzzz_1,  \
                             g_0_xyyyyzz_0_yzzzz_1,   \
                             g_0_xyyyyzz_0_yzzzzz_0,  \
                             g_0_xyyyyzz_0_yzzzzz_1,  \
                             g_0_xyyyyzz_0_zzzzzz_0,  \
                             g_0_xyyyyzz_0_zzzzzz_1,  \
                             g_0_yyyyzz_0_xxxxyz_0,   \
                             g_0_yyyyzz_0_xxxxyz_1,   \
                             g_0_yyyyzz_0_xxxyyz_0,   \
                             g_0_yyyyzz_0_xxxyyz_1,   \
                             g_0_yyyyzz_0_xxxyzz_0,   \
                             g_0_yyyyzz_0_xxxyzz_1,   \
                             g_0_yyyyzz_0_xxyyyz_0,   \
                             g_0_yyyyzz_0_xxyyyz_1,   \
                             g_0_yyyyzz_0_xxyyzz_0,   \
                             g_0_yyyyzz_0_xxyyzz_1,   \
                             g_0_yyyyzz_0_xxyzzz_0,   \
                             g_0_yyyyzz_0_xxyzzz_1,   \
                             g_0_yyyyzz_0_xyyyyz_0,   \
                             g_0_yyyyzz_0_xyyyyz_1,   \
                             g_0_yyyyzz_0_xyyyzz_0,   \
                             g_0_yyyyzz_0_xyyyzz_1,   \
                             g_0_yyyyzz_0_xyyzzz_0,   \
                             g_0_yyyyzz_0_xyyzzz_1,   \
                             g_0_yyyyzz_0_xyzzzz_0,   \
                             g_0_yyyyzz_0_xyzzzz_1,   \
                             g_0_yyyyzz_0_yyyyyy_0,   \
                             g_0_yyyyzz_0_yyyyyy_1,   \
                             g_0_yyyyzz_0_yyyyyz_0,   \
                             g_0_yyyyzz_0_yyyyyz_1,   \
                             g_0_yyyyzz_0_yyyyzz_0,   \
                             g_0_yyyyzz_0_yyyyzz_1,   \
                             g_0_yyyyzz_0_yyyzzz_0,   \
                             g_0_yyyyzz_0_yyyzzz_1,   \
                             g_0_yyyyzz_0_yyzzzz_0,   \
                             g_0_yyyyzz_0_yyzzzz_1,   \
                             g_0_yyyyzz_0_yzzzzz_0,   \
                             g_0_yyyyzz_0_yzzzzz_1,   \
                             g_0_yyyyzz_0_zzzzzz_0,   \
                             g_0_yyyyzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_xxxxxx_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxxxx_0[i] * pb_y + g_0_xxyyyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxxxy_0[i] = g_0_xxyyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxxxy_0[i] * pb_z +
                                     g_0_xxyyyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxxxz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxxxz_0[i] * pb_y + g_0_xxyyyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxxyy_0[i] = g_0_xxyyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxxyy_0[i] * pb_z +
                                     g_0_xxyyyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxxyz_0[i] = g_0_yyyyzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxxzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxxzz_0[i] * pb_y + g_0_xxyyyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxyyy_0[i] = g_0_xxyyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxyyy_0[i] * pb_z +
                                     g_0_xxyyyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxyyz_0[i] = g_0_yyyyzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxyzz_0[i] = g_0_yyyyzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxzzz_0[i] * pb_y + g_0_xxyyyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxyyyy_0[i] = g_0_xxyyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxyyyy_0[i] * pb_z +
                                     g_0_xxyyyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxyyyz_0[i] = g_0_yyyyzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxyyzz_0[i] = g_0_yyyyzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxyzzz_0[i] = g_0_yyyyzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxzzzz_0[i] * pb_y + g_0_xxyyyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xyyyyy_0[i] = g_0_xxyyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xyyyyy_0[i] * pb_z +
                                     g_0_xxyyyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xyyyyz_0[i] = g_0_yyyyzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyzz_0_xyyyyz_0[i] * pb_x + g_0_xyyyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyyyzz_0[i] = g_0_yyyyzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyzz_0_xyyyzz_0[i] * pb_x + g_0_xyyyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyyzzz_0[i] = g_0_yyyyzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyzz_0_xyyzzz_0[i] * pb_x + g_0_xyyyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyzzzz_0[i] = g_0_yyyyzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyyzz_0_xyzzzz_0[i] * pb_x + g_0_xyyyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xzzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xzzzzz_0[i] * pb_y + g_0_xxyyyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_yyyyyy_0[i] = g_0_yyyyzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyyy_0[i] * pb_x +
                                     g_0_xyyyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyyyz_0[i] = g_0_yyyyzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyyz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyyzz_0[i] = g_0_yyyyzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyzzz_0[i] = g_0_yyyyzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyzzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyzzzz_0[i] = g_0_yyyyzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzzzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yzzzzz_0[i] = g_0_yyyyzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzzzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_zzzzzz_0[i] = g_0_yyyyzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_zzzzzz_0[i] * pb_x +
                                     g_0_xyyyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 672-700 components of targeted buffer : SLSI

    auto g_0_xxyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 672);

    auto g_0_xxyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 673);

    auto g_0_xxyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 674);

    auto g_0_xxyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 675);

    auto g_0_xxyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 676);

    auto g_0_xxyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 677);

    auto g_0_xxyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 678);

    auto g_0_xxyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 679);

    auto g_0_xxyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 680);

    auto g_0_xxyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 681);

    auto g_0_xxyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 682);

    auto g_0_xxyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 683);

    auto g_0_xxyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 684);

    auto g_0_xxyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 685);

    auto g_0_xxyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 686);

    auto g_0_xxyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 687);

    auto g_0_xxyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 688);

    auto g_0_xxyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 689);

    auto g_0_xxyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 690);

    auto g_0_xxyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 691);

    auto g_0_xxyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 692);

    auto g_0_xxyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 693);

    auto g_0_xxyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 694);

    auto g_0_xxyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 695);

    auto g_0_xxyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 696);

    auto g_0_xxyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 697);

    auto g_0_xxyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 698);

    auto g_0_xxyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 699);

#pragma omp simd aligned(g_0_xxyyyz_0_xxxxxy_0,       \
                             g_0_xxyyyz_0_xxxxxy_1,   \
                             g_0_xxyyyz_0_xxxxyy_0,   \
                             g_0_xxyyyz_0_xxxxyy_1,   \
                             g_0_xxyyyz_0_xxxyyy_0,   \
                             g_0_xxyyyz_0_xxxyyy_1,   \
                             g_0_xxyyyz_0_xxyyyy_0,   \
                             g_0_xxyyyz_0_xxyyyy_1,   \
                             g_0_xxyyyz_0_xyyyyy_0,   \
                             g_0_xxyyyz_0_xyyyyy_1,   \
                             g_0_xxyyyzz_0_xxxxxy_0,  \
                             g_0_xxyyyzz_0_xxxxxy_1,  \
                             g_0_xxyyyzz_0_xxxxyy_0,  \
                             g_0_xxyyyzz_0_xxxxyy_1,  \
                             g_0_xxyyyzz_0_xxxyyy_0,  \
                             g_0_xxyyyzz_0_xxxyyy_1,  \
                             g_0_xxyyyzz_0_xxyyyy_0,  \
                             g_0_xxyyyzz_0_xxyyyy_1,  \
                             g_0_xxyyyzz_0_xyyyyy_0,  \
                             g_0_xxyyyzz_0_xyyyyy_1,  \
                             g_0_xxyyyzzz_0_xxxxxx_0, \
                             g_0_xxyyyzzz_0_xxxxxy_0, \
                             g_0_xxyyyzzz_0_xxxxxz_0, \
                             g_0_xxyyyzzz_0_xxxxyy_0, \
                             g_0_xxyyyzzz_0_xxxxyz_0, \
                             g_0_xxyyyzzz_0_xxxxzz_0, \
                             g_0_xxyyyzzz_0_xxxyyy_0, \
                             g_0_xxyyyzzz_0_xxxyyz_0, \
                             g_0_xxyyyzzz_0_xxxyzz_0, \
                             g_0_xxyyyzzz_0_xxxzzz_0, \
                             g_0_xxyyyzzz_0_xxyyyy_0, \
                             g_0_xxyyyzzz_0_xxyyyz_0, \
                             g_0_xxyyyzzz_0_xxyyzz_0, \
                             g_0_xxyyyzzz_0_xxyzzz_0, \
                             g_0_xxyyyzzz_0_xxzzzz_0, \
                             g_0_xxyyyzzz_0_xyyyyy_0, \
                             g_0_xxyyyzzz_0_xyyyyz_0, \
                             g_0_xxyyyzzz_0_xyyyzz_0, \
                             g_0_xxyyyzzz_0_xyyzzz_0, \
                             g_0_xxyyyzzz_0_xyzzzz_0, \
                             g_0_xxyyyzzz_0_xzzzzz_0, \
                             g_0_xxyyyzzz_0_yyyyyy_0, \
                             g_0_xxyyyzzz_0_yyyyyz_0, \
                             g_0_xxyyyzzz_0_yyyyzz_0, \
                             g_0_xxyyyzzz_0_yyyzzz_0, \
                             g_0_xxyyyzzz_0_yyzzzz_0, \
                             g_0_xxyyyzzz_0_yzzzzz_0, \
                             g_0_xxyyyzzz_0_zzzzzz_0, \
                             g_0_xxyyzzz_0_xxxxxx_0,  \
                             g_0_xxyyzzz_0_xxxxxx_1,  \
                             g_0_xxyyzzz_0_xxxxxz_0,  \
                             g_0_xxyyzzz_0_xxxxxz_1,  \
                             g_0_xxyyzzz_0_xxxxzz_0,  \
                             g_0_xxyyzzz_0_xxxxzz_1,  \
                             g_0_xxyyzzz_0_xxxzzz_0,  \
                             g_0_xxyyzzz_0_xxxzzz_1,  \
                             g_0_xxyyzzz_0_xxzzzz_0,  \
                             g_0_xxyyzzz_0_xxzzzz_1,  \
                             g_0_xxyyzzz_0_xzzzzz_0,  \
                             g_0_xxyyzzz_0_xzzzzz_1,  \
                             g_0_xxyzzz_0_xxxxxx_0,   \
                             g_0_xxyzzz_0_xxxxxx_1,   \
                             g_0_xxyzzz_0_xxxxxz_0,   \
                             g_0_xxyzzz_0_xxxxxz_1,   \
                             g_0_xxyzzz_0_xxxxzz_0,   \
                             g_0_xxyzzz_0_xxxxzz_1,   \
                             g_0_xxyzzz_0_xxxzzz_0,   \
                             g_0_xxyzzz_0_xxxzzz_1,   \
                             g_0_xxyzzz_0_xxzzzz_0,   \
                             g_0_xxyzzz_0_xxzzzz_1,   \
                             g_0_xxyzzz_0_xzzzzz_0,   \
                             g_0_xxyzzz_0_xzzzzz_1,   \
                             g_0_xyyyzzz_0_xxxxyz_0,  \
                             g_0_xyyyzzz_0_xxxxyz_1,  \
                             g_0_xyyyzzz_0_xxxyyz_0,  \
                             g_0_xyyyzzz_0_xxxyyz_1,  \
                             g_0_xyyyzzz_0_xxxyz_1,   \
                             g_0_xyyyzzz_0_xxxyzz_0,  \
                             g_0_xyyyzzz_0_xxxyzz_1,  \
                             g_0_xyyyzzz_0_xxyyyz_0,  \
                             g_0_xyyyzzz_0_xxyyyz_1,  \
                             g_0_xyyyzzz_0_xxyyz_1,   \
                             g_0_xyyyzzz_0_xxyyzz_0,  \
                             g_0_xyyyzzz_0_xxyyzz_1,  \
                             g_0_xyyyzzz_0_xxyzz_1,   \
                             g_0_xyyyzzz_0_xxyzzz_0,  \
                             g_0_xyyyzzz_0_xxyzzz_1,  \
                             g_0_xyyyzzz_0_xyyyyz_0,  \
                             g_0_xyyyzzz_0_xyyyyz_1,  \
                             g_0_xyyyzzz_0_xyyyz_1,   \
                             g_0_xyyyzzz_0_xyyyzz_0,  \
                             g_0_xyyyzzz_0_xyyyzz_1,  \
                             g_0_xyyyzzz_0_xyyzz_1,   \
                             g_0_xyyyzzz_0_xyyzzz_0,  \
                             g_0_xyyyzzz_0_xyyzzz_1,  \
                             g_0_xyyyzzz_0_xyzzz_1,   \
                             g_0_xyyyzzz_0_xyzzzz_0,  \
                             g_0_xyyyzzz_0_xyzzzz_1,  \
                             g_0_xyyyzzz_0_yyyyyy_0,  \
                             g_0_xyyyzzz_0_yyyyyy_1,  \
                             g_0_xyyyzzz_0_yyyyyz_0,  \
                             g_0_xyyyzzz_0_yyyyyz_1,  \
                             g_0_xyyyzzz_0_yyyyz_1,   \
                             g_0_xyyyzzz_0_yyyyzz_0,  \
                             g_0_xyyyzzz_0_yyyyzz_1,  \
                             g_0_xyyyzzz_0_yyyzz_1,   \
                             g_0_xyyyzzz_0_yyyzzz_0,  \
                             g_0_xyyyzzz_0_yyyzzz_1,  \
                             g_0_xyyyzzz_0_yyzzz_1,   \
                             g_0_xyyyzzz_0_yyzzzz_0,  \
                             g_0_xyyyzzz_0_yyzzzz_1,  \
                             g_0_xyyyzzz_0_yzzzz_1,   \
                             g_0_xyyyzzz_0_yzzzzz_0,  \
                             g_0_xyyyzzz_0_yzzzzz_1,  \
                             g_0_xyyyzzz_0_zzzzzz_0,  \
                             g_0_xyyyzzz_0_zzzzzz_1,  \
                             g_0_yyyzzz_0_xxxxyz_0,   \
                             g_0_yyyzzz_0_xxxxyz_1,   \
                             g_0_yyyzzz_0_xxxyyz_0,   \
                             g_0_yyyzzz_0_xxxyyz_1,   \
                             g_0_yyyzzz_0_xxxyzz_0,   \
                             g_0_yyyzzz_0_xxxyzz_1,   \
                             g_0_yyyzzz_0_xxyyyz_0,   \
                             g_0_yyyzzz_0_xxyyyz_1,   \
                             g_0_yyyzzz_0_xxyyzz_0,   \
                             g_0_yyyzzz_0_xxyyzz_1,   \
                             g_0_yyyzzz_0_xxyzzz_0,   \
                             g_0_yyyzzz_0_xxyzzz_1,   \
                             g_0_yyyzzz_0_xyyyyz_0,   \
                             g_0_yyyzzz_0_xyyyyz_1,   \
                             g_0_yyyzzz_0_xyyyzz_0,   \
                             g_0_yyyzzz_0_xyyyzz_1,   \
                             g_0_yyyzzz_0_xyyzzz_0,   \
                             g_0_yyyzzz_0_xyyzzz_1,   \
                             g_0_yyyzzz_0_xyzzzz_0,   \
                             g_0_yyyzzz_0_xyzzzz_1,   \
                             g_0_yyyzzz_0_yyyyyy_0,   \
                             g_0_yyyzzz_0_yyyyyy_1,   \
                             g_0_yyyzzz_0_yyyyyz_0,   \
                             g_0_yyyzzz_0_yyyyyz_1,   \
                             g_0_yyyzzz_0_yyyyzz_0,   \
                             g_0_yyyzzz_0_yyyyzz_1,   \
                             g_0_yyyzzz_0_yyyzzz_0,   \
                             g_0_yyyzzz_0_yyyzzz_1,   \
                             g_0_yyyzzz_0_yyzzzz_0,   \
                             g_0_yyyzzz_0_yyzzzz_1,   \
                             g_0_yyyzzz_0_yzzzzz_0,   \
                             g_0_yyyzzz_0_yzzzzz_1,   \
                             g_0_yyyzzz_0_zzzzzz_0,   \
                             g_0_yyyzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_xxxxxx_0[i] = 2.0 * g_0_xxyzzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxxxx_0[i] * pb_y + g_0_xxyyzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxxxy_0[i] = 2.0 * g_0_xxyyyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxxxy_0[i] * pb_z + g_0_xxyyyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxxxz_0[i] = 2.0 * g_0_xxyzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxxxz_0[i] * pb_y + g_0_xxyyzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxxyy_0[i] = 2.0 * g_0_xxyyyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxxyy_0[i] * pb_z + g_0_xxyyyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxxyz_0[i] = g_0_yyyzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxxzz_0[i] = 2.0 * g_0_xxyzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxxzz_0[i] * pb_y + g_0_xxyyzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxyyy_0[i] = 2.0 * g_0_xxyyyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxxyyy_0[i] * pb_z + g_0_xxyyyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxyyz_0[i] = g_0_yyyzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxyzz_0[i] = g_0_yyyzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxzzz_0[i] = 2.0 * g_0_xxyzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxzzz_0[i] * pb_y + g_0_xxyyzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_xxyyyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xxyyyy_0[i] * pb_z + g_0_xxyyyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxyyyz_0[i] = g_0_yyyzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxyyzz_0[i] = g_0_yyyzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxyzzz_0[i] = g_0_yyyzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxzzzz_0[i] = 2.0 * g_0_xxyzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxzzzz_0[i] * pb_y + g_0_xxyyzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xyyyyy_0[i] = 2.0 * g_0_xxyyyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyyzz_0_xyyyyy_0[i] * pb_z + g_0_xxyyyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xyyyyz_0[i] = g_0_yyyzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                     g_0_xyyyzzz_0_xyyyyz_0[i] * pb_x + g_0_xyyyzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyyyzz_0[i] = g_0_yyyzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyzzz_0_xyyyzz_0[i] * pb_x + g_0_xyyyzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyyzzz_0[i] = g_0_yyyzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyzzz_0_xyyzzz_0[i] * pb_x + g_0_xyyyzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyzzzz_0[i] = g_0_yyyzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyyzzz_0_xyzzzz_0[i] * pb_x + g_0_xyyyzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xzzzzz_0[i] = 2.0 * g_0_xxyzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xzzzzz_0[i] * pb_y + g_0_xxyyzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_yyyyyy_0[i] = g_0_yyyzzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyyy_0[i] * pb_x +
                                     g_0_xyyyzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyyyz_0[i] = g_0_yyyzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyyz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyyzz_0[i] = g_0_yyyzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyzzz_0[i] = g_0_yyyzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyzzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyzzzz_0[i] = g_0_yyyzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzzzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yzzzzz_0[i] = g_0_yyyzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzzzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_zzzzzz_0[i] = g_0_yyyzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_zzzzzz_0[i] * pb_x +
                                     g_0_xyyyzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 700-728 components of targeted buffer : SLSI

    auto g_0_xxyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 700);

    auto g_0_xxyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 701);

    auto g_0_xxyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 702);

    auto g_0_xxyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 703);

    auto g_0_xxyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 704);

    auto g_0_xxyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 705);

    auto g_0_xxyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 706);

    auto g_0_xxyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 707);

    auto g_0_xxyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 708);

    auto g_0_xxyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 709);

    auto g_0_xxyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 710);

    auto g_0_xxyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 711);

    auto g_0_xxyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 712);

    auto g_0_xxyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 713);

    auto g_0_xxyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 714);

    auto g_0_xxyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 715);

    auto g_0_xxyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 716);

    auto g_0_xxyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 717);

    auto g_0_xxyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 718);

    auto g_0_xxyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 719);

    auto g_0_xxyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 720);

    auto g_0_xxyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 721);

    auto g_0_xxyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 722);

    auto g_0_xxyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 723);

    auto g_0_xxyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 724);

    auto g_0_xxyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 725);

    auto g_0_xxyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 726);

    auto g_0_xxyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 727);

#pragma omp simd aligned(g_0_xxyyzz_0_xxxxxy_0,       \
                             g_0_xxyyzz_0_xxxxxy_1,   \
                             g_0_xxyyzz_0_xxxxyy_0,   \
                             g_0_xxyyzz_0_xxxxyy_1,   \
                             g_0_xxyyzz_0_xxxyyy_0,   \
                             g_0_xxyyzz_0_xxxyyy_1,   \
                             g_0_xxyyzz_0_xxyyyy_0,   \
                             g_0_xxyyzz_0_xxyyyy_1,   \
                             g_0_xxyyzz_0_xyyyyy_0,   \
                             g_0_xxyyzz_0_xyyyyy_1,   \
                             g_0_xxyyzzz_0_xxxxxy_0,  \
                             g_0_xxyyzzz_0_xxxxxy_1,  \
                             g_0_xxyyzzz_0_xxxxyy_0,  \
                             g_0_xxyyzzz_0_xxxxyy_1,  \
                             g_0_xxyyzzz_0_xxxyyy_0,  \
                             g_0_xxyyzzz_0_xxxyyy_1,  \
                             g_0_xxyyzzz_0_xxyyyy_0,  \
                             g_0_xxyyzzz_0_xxyyyy_1,  \
                             g_0_xxyyzzz_0_xyyyyy_0,  \
                             g_0_xxyyzzz_0_xyyyyy_1,  \
                             g_0_xxyyzzzz_0_xxxxxx_0, \
                             g_0_xxyyzzzz_0_xxxxxy_0, \
                             g_0_xxyyzzzz_0_xxxxxz_0, \
                             g_0_xxyyzzzz_0_xxxxyy_0, \
                             g_0_xxyyzzzz_0_xxxxyz_0, \
                             g_0_xxyyzzzz_0_xxxxzz_0, \
                             g_0_xxyyzzzz_0_xxxyyy_0, \
                             g_0_xxyyzzzz_0_xxxyyz_0, \
                             g_0_xxyyzzzz_0_xxxyzz_0, \
                             g_0_xxyyzzzz_0_xxxzzz_0, \
                             g_0_xxyyzzzz_0_xxyyyy_0, \
                             g_0_xxyyzzzz_0_xxyyyz_0, \
                             g_0_xxyyzzzz_0_xxyyzz_0, \
                             g_0_xxyyzzzz_0_xxyzzz_0, \
                             g_0_xxyyzzzz_0_xxzzzz_0, \
                             g_0_xxyyzzzz_0_xyyyyy_0, \
                             g_0_xxyyzzzz_0_xyyyyz_0, \
                             g_0_xxyyzzzz_0_xyyyzz_0, \
                             g_0_xxyyzzzz_0_xyyzzz_0, \
                             g_0_xxyyzzzz_0_xyzzzz_0, \
                             g_0_xxyyzzzz_0_xzzzzz_0, \
                             g_0_xxyyzzzz_0_yyyyyy_0, \
                             g_0_xxyyzzzz_0_yyyyyz_0, \
                             g_0_xxyyzzzz_0_yyyyzz_0, \
                             g_0_xxyyzzzz_0_yyyzzz_0, \
                             g_0_xxyyzzzz_0_yyzzzz_0, \
                             g_0_xxyyzzzz_0_yzzzzz_0, \
                             g_0_xxyyzzzz_0_zzzzzz_0, \
                             g_0_xxyzzzz_0_xxxxxx_0,  \
                             g_0_xxyzzzz_0_xxxxxx_1,  \
                             g_0_xxyzzzz_0_xxxxxz_0,  \
                             g_0_xxyzzzz_0_xxxxxz_1,  \
                             g_0_xxyzzzz_0_xxxxzz_0,  \
                             g_0_xxyzzzz_0_xxxxzz_1,  \
                             g_0_xxyzzzz_0_xxxzzz_0,  \
                             g_0_xxyzzzz_0_xxxzzz_1,  \
                             g_0_xxyzzzz_0_xxzzzz_0,  \
                             g_0_xxyzzzz_0_xxzzzz_1,  \
                             g_0_xxyzzzz_0_xzzzzz_0,  \
                             g_0_xxyzzzz_0_xzzzzz_1,  \
                             g_0_xxzzzz_0_xxxxxx_0,   \
                             g_0_xxzzzz_0_xxxxxx_1,   \
                             g_0_xxzzzz_0_xxxxxz_0,   \
                             g_0_xxzzzz_0_xxxxxz_1,   \
                             g_0_xxzzzz_0_xxxxzz_0,   \
                             g_0_xxzzzz_0_xxxxzz_1,   \
                             g_0_xxzzzz_0_xxxzzz_0,   \
                             g_0_xxzzzz_0_xxxzzz_1,   \
                             g_0_xxzzzz_0_xxzzzz_0,   \
                             g_0_xxzzzz_0_xxzzzz_1,   \
                             g_0_xxzzzz_0_xzzzzz_0,   \
                             g_0_xxzzzz_0_xzzzzz_1,   \
                             g_0_xyyzzzz_0_xxxxyz_0,  \
                             g_0_xyyzzzz_0_xxxxyz_1,  \
                             g_0_xyyzzzz_0_xxxyyz_0,  \
                             g_0_xyyzzzz_0_xxxyyz_1,  \
                             g_0_xyyzzzz_0_xxxyz_1,   \
                             g_0_xyyzzzz_0_xxxyzz_0,  \
                             g_0_xyyzzzz_0_xxxyzz_1,  \
                             g_0_xyyzzzz_0_xxyyyz_0,  \
                             g_0_xyyzzzz_0_xxyyyz_1,  \
                             g_0_xyyzzzz_0_xxyyz_1,   \
                             g_0_xyyzzzz_0_xxyyzz_0,  \
                             g_0_xyyzzzz_0_xxyyzz_1,  \
                             g_0_xyyzzzz_0_xxyzz_1,   \
                             g_0_xyyzzzz_0_xxyzzz_0,  \
                             g_0_xyyzzzz_0_xxyzzz_1,  \
                             g_0_xyyzzzz_0_xyyyyz_0,  \
                             g_0_xyyzzzz_0_xyyyyz_1,  \
                             g_0_xyyzzzz_0_xyyyz_1,   \
                             g_0_xyyzzzz_0_xyyyzz_0,  \
                             g_0_xyyzzzz_0_xyyyzz_1,  \
                             g_0_xyyzzzz_0_xyyzz_1,   \
                             g_0_xyyzzzz_0_xyyzzz_0,  \
                             g_0_xyyzzzz_0_xyyzzz_1,  \
                             g_0_xyyzzzz_0_xyzzz_1,   \
                             g_0_xyyzzzz_0_xyzzzz_0,  \
                             g_0_xyyzzzz_0_xyzzzz_1,  \
                             g_0_xyyzzzz_0_yyyyyy_0,  \
                             g_0_xyyzzzz_0_yyyyyy_1,  \
                             g_0_xyyzzzz_0_yyyyyz_0,  \
                             g_0_xyyzzzz_0_yyyyyz_1,  \
                             g_0_xyyzzzz_0_yyyyz_1,   \
                             g_0_xyyzzzz_0_yyyyzz_0,  \
                             g_0_xyyzzzz_0_yyyyzz_1,  \
                             g_0_xyyzzzz_0_yyyzz_1,   \
                             g_0_xyyzzzz_0_yyyzzz_0,  \
                             g_0_xyyzzzz_0_yyyzzz_1,  \
                             g_0_xyyzzzz_0_yyzzz_1,   \
                             g_0_xyyzzzz_0_yyzzzz_0,  \
                             g_0_xyyzzzz_0_yyzzzz_1,  \
                             g_0_xyyzzzz_0_yzzzz_1,   \
                             g_0_xyyzzzz_0_yzzzzz_0,  \
                             g_0_xyyzzzz_0_yzzzzz_1,  \
                             g_0_xyyzzzz_0_zzzzzz_0,  \
                             g_0_xyyzzzz_0_zzzzzz_1,  \
                             g_0_yyzzzz_0_xxxxyz_0,   \
                             g_0_yyzzzz_0_xxxxyz_1,   \
                             g_0_yyzzzz_0_xxxyyz_0,   \
                             g_0_yyzzzz_0_xxxyyz_1,   \
                             g_0_yyzzzz_0_xxxyzz_0,   \
                             g_0_yyzzzz_0_xxxyzz_1,   \
                             g_0_yyzzzz_0_xxyyyz_0,   \
                             g_0_yyzzzz_0_xxyyyz_1,   \
                             g_0_yyzzzz_0_xxyyzz_0,   \
                             g_0_yyzzzz_0_xxyyzz_1,   \
                             g_0_yyzzzz_0_xxyzzz_0,   \
                             g_0_yyzzzz_0_xxyzzz_1,   \
                             g_0_yyzzzz_0_xyyyyz_0,   \
                             g_0_yyzzzz_0_xyyyyz_1,   \
                             g_0_yyzzzz_0_xyyyzz_0,   \
                             g_0_yyzzzz_0_xyyyzz_1,   \
                             g_0_yyzzzz_0_xyyzzz_0,   \
                             g_0_yyzzzz_0_xyyzzz_1,   \
                             g_0_yyzzzz_0_xyzzzz_0,   \
                             g_0_yyzzzz_0_xyzzzz_1,   \
                             g_0_yyzzzz_0_yyyyyy_0,   \
                             g_0_yyzzzz_0_yyyyyy_1,   \
                             g_0_yyzzzz_0_yyyyyz_0,   \
                             g_0_yyzzzz_0_yyyyyz_1,   \
                             g_0_yyzzzz_0_yyyyzz_0,   \
                             g_0_yyzzzz_0_yyyyzz_1,   \
                             g_0_yyzzzz_0_yyyzzz_0,   \
                             g_0_yyzzzz_0_yyyzzz_1,   \
                             g_0_yyzzzz_0_yyzzzz_0,   \
                             g_0_yyzzzz_0_yyzzzz_1,   \
                             g_0_yyzzzz_0_yzzzzz_0,   \
                             g_0_yyzzzz_0_yzzzzz_1,   \
                             g_0_yyzzzz_0_zzzzzz_0,   \
                             g_0_yyzzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_xxxxxx_0[i] = g_0_xxzzzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxxx_0[i] * pb_y +
                                     g_0_xxyzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxxxy_0[i] = 3.0 * g_0_xxyyzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxxxy_0[i] * pb_z + g_0_xxyyzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxxxz_0[i] = g_0_xxzzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxxz_0[i] * pb_y +
                                     g_0_xxyzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxxyy_0[i] = 3.0 * g_0_xxyyzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxxyy_0[i] * pb_z + g_0_xxyyzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxxyz_0[i] = g_0_yyzzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xyyzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxxzz_0[i] = g_0_xxzzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxzz_0[i] * pb_y +
                                     g_0_xxyzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_xxyyzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxxyyy_0[i] * pb_z + g_0_xxyyzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxyyz_0[i] = g_0_yyzzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxyzz_0[i] = g_0_yyzzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xyyzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxzzz_0[i] = g_0_xxzzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxzzz_0[i] * pb_y +
                                     g_0_xxyzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxyyyy_0[i] = 3.0 * g_0_xxyyzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xxyyyy_0[i] * pb_z + g_0_xxyyzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxyyyz_0[i] = g_0_yyzzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxyyzz_0[i] = g_0_yyzzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxyzzz_0[i] = g_0_yyzzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xyyzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxzzzz_0[i] = g_0_xxzzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxzzzz_0[i] * pb_y +
                                     g_0_xxyzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xyyyyy_0[i] = 3.0 * g_0_xxyyzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxyyzzz_0_xyyyyy_0[i] * pb_z + g_0_xxyyzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xyyyyz_0[i] = g_0_yyzzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                     g_0_xyyzzzz_0_xyyyyz_0[i] * pb_x + g_0_xyyzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyyyzz_0[i] = g_0_yyzzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzzzz_0_xyyyzz_0[i] * pb_x + g_0_xyyzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyyzzz_0[i] = g_0_yyzzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzzzz_0_xyyzzz_0[i] * pb_x + g_0_xyyzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyzzzz_0[i] = g_0_yyzzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                     g_0_xyyzzzz_0_xyzzzz_0[i] * pb_x + g_0_xyyzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xzzzzz_0[i] = g_0_xxzzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xzzzzz_0[i] * pb_y +
                                     g_0_xxyzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_yyyyyy_0[i] = g_0_yyzzzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyyy_0[i] * pb_x +
                                     g_0_xyyzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyyyz_0[i] = g_0_yyzzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyyz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyyzz_0[i] = g_0_yyzzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyzzz_0[i] = g_0_yyzzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyzzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyzzzz_0[i] = g_0_yyzzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzzzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yzzzzz_0[i] = g_0_yyzzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzzzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_zzzzzz_0[i] = g_0_yyzzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_zzzzzz_0[i] * pb_x +
                                     g_0_xyyzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 728-756 components of targeted buffer : SLSI

    auto g_0_xxyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 728);

    auto g_0_xxyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 729);

    auto g_0_xxyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 730);

    auto g_0_xxyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 731);

    auto g_0_xxyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 732);

    auto g_0_xxyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 733);

    auto g_0_xxyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 734);

    auto g_0_xxyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 735);

    auto g_0_xxyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 736);

    auto g_0_xxyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 737);

    auto g_0_xxyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 738);

    auto g_0_xxyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 739);

    auto g_0_xxyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 740);

    auto g_0_xxyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 741);

    auto g_0_xxyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 742);

    auto g_0_xxyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 743);

    auto g_0_xxyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 744);

    auto g_0_xxyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 745);

    auto g_0_xxyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 746);

    auto g_0_xxyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 747);

    auto g_0_xxyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 748);

    auto g_0_xxyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 749);

    auto g_0_xxyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 750);

    auto g_0_xxyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 751);

    auto g_0_xxyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 752);

    auto g_0_xxyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 753);

    auto g_0_xxyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 754);

    auto g_0_xxyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 755);

#pragma omp simd aligned(g_0_xxyzzzzz_0_xxxxxx_0,     \
                             g_0_xxyzzzzz_0_xxxxxy_0, \
                             g_0_xxyzzzzz_0_xxxxxz_0, \
                             g_0_xxyzzzzz_0_xxxxyy_0, \
                             g_0_xxyzzzzz_0_xxxxyz_0, \
                             g_0_xxyzzzzz_0_xxxxzz_0, \
                             g_0_xxyzzzzz_0_xxxyyy_0, \
                             g_0_xxyzzzzz_0_xxxyyz_0, \
                             g_0_xxyzzzzz_0_xxxyzz_0, \
                             g_0_xxyzzzzz_0_xxxzzz_0, \
                             g_0_xxyzzzzz_0_xxyyyy_0, \
                             g_0_xxyzzzzz_0_xxyyyz_0, \
                             g_0_xxyzzzzz_0_xxyyzz_0, \
                             g_0_xxyzzzzz_0_xxyzzz_0, \
                             g_0_xxyzzzzz_0_xxzzzz_0, \
                             g_0_xxyzzzzz_0_xyyyyy_0, \
                             g_0_xxyzzzzz_0_xyyyyz_0, \
                             g_0_xxyzzzzz_0_xyyyzz_0, \
                             g_0_xxyzzzzz_0_xyyzzz_0, \
                             g_0_xxyzzzzz_0_xyzzzz_0, \
                             g_0_xxyzzzzz_0_xzzzzz_0, \
                             g_0_xxyzzzzz_0_yyyyyy_0, \
                             g_0_xxyzzzzz_0_yyyyyz_0, \
                             g_0_xxyzzzzz_0_yyyyzz_0, \
                             g_0_xxyzzzzz_0_yyyzzz_0, \
                             g_0_xxyzzzzz_0_yyzzzz_0, \
                             g_0_xxyzzzzz_0_yzzzzz_0, \
                             g_0_xxyzzzzz_0_zzzzzz_0, \
                             g_0_xxzzzzz_0_xxxxx_1,   \
                             g_0_xxzzzzz_0_xxxxxx_0,  \
                             g_0_xxzzzzz_0_xxxxxx_1,  \
                             g_0_xxzzzzz_0_xxxxxy_0,  \
                             g_0_xxzzzzz_0_xxxxxy_1,  \
                             g_0_xxzzzzz_0_xxxxxz_0,  \
                             g_0_xxzzzzz_0_xxxxxz_1,  \
                             g_0_xxzzzzz_0_xxxxy_1,   \
                             g_0_xxzzzzz_0_xxxxyy_0,  \
                             g_0_xxzzzzz_0_xxxxyy_1,  \
                             g_0_xxzzzzz_0_xxxxyz_0,  \
                             g_0_xxzzzzz_0_xxxxyz_1,  \
                             g_0_xxzzzzz_0_xxxxz_1,   \
                             g_0_xxzzzzz_0_xxxxzz_0,  \
                             g_0_xxzzzzz_0_xxxxzz_1,  \
                             g_0_xxzzzzz_0_xxxyy_1,   \
                             g_0_xxzzzzz_0_xxxyyy_0,  \
                             g_0_xxzzzzz_0_xxxyyy_1,  \
                             g_0_xxzzzzz_0_xxxyyz_0,  \
                             g_0_xxzzzzz_0_xxxyyz_1,  \
                             g_0_xxzzzzz_0_xxxyz_1,   \
                             g_0_xxzzzzz_0_xxxyzz_0,  \
                             g_0_xxzzzzz_0_xxxyzz_1,  \
                             g_0_xxzzzzz_0_xxxzz_1,   \
                             g_0_xxzzzzz_0_xxxzzz_0,  \
                             g_0_xxzzzzz_0_xxxzzz_1,  \
                             g_0_xxzzzzz_0_xxyyy_1,   \
                             g_0_xxzzzzz_0_xxyyyy_0,  \
                             g_0_xxzzzzz_0_xxyyyy_1,  \
                             g_0_xxzzzzz_0_xxyyyz_0,  \
                             g_0_xxzzzzz_0_xxyyyz_1,  \
                             g_0_xxzzzzz_0_xxyyz_1,   \
                             g_0_xxzzzzz_0_xxyyzz_0,  \
                             g_0_xxzzzzz_0_xxyyzz_1,  \
                             g_0_xxzzzzz_0_xxyzz_1,   \
                             g_0_xxzzzzz_0_xxyzzz_0,  \
                             g_0_xxzzzzz_0_xxyzzz_1,  \
                             g_0_xxzzzzz_0_xxzzz_1,   \
                             g_0_xxzzzzz_0_xxzzzz_0,  \
                             g_0_xxzzzzz_0_xxzzzz_1,  \
                             g_0_xxzzzzz_0_xyyyy_1,   \
                             g_0_xxzzzzz_0_xyyyyy_0,  \
                             g_0_xxzzzzz_0_xyyyyy_1,  \
                             g_0_xxzzzzz_0_xyyyyz_0,  \
                             g_0_xxzzzzz_0_xyyyyz_1,  \
                             g_0_xxzzzzz_0_xyyyz_1,   \
                             g_0_xxzzzzz_0_xyyyzz_0,  \
                             g_0_xxzzzzz_0_xyyyzz_1,  \
                             g_0_xxzzzzz_0_xyyzz_1,   \
                             g_0_xxzzzzz_0_xyyzzz_0,  \
                             g_0_xxzzzzz_0_xyyzzz_1,  \
                             g_0_xxzzzzz_0_xyzzz_1,   \
                             g_0_xxzzzzz_0_xyzzzz_0,  \
                             g_0_xxzzzzz_0_xyzzzz_1,  \
                             g_0_xxzzzzz_0_xzzzz_1,   \
                             g_0_xxzzzzz_0_xzzzzz_0,  \
                             g_0_xxzzzzz_0_xzzzzz_1,  \
                             g_0_xxzzzzz_0_yyyyy_1,   \
                             g_0_xxzzzzz_0_yyyyyy_0,  \
                             g_0_xxzzzzz_0_yyyyyy_1,  \
                             g_0_xxzzzzz_0_yyyyyz_0,  \
                             g_0_xxzzzzz_0_yyyyyz_1,  \
                             g_0_xxzzzzz_0_yyyyz_1,   \
                             g_0_xxzzzzz_0_yyyyzz_0,  \
                             g_0_xxzzzzz_0_yyyyzz_1,  \
                             g_0_xxzzzzz_0_yyyzz_1,   \
                             g_0_xxzzzzz_0_yyyzzz_0,  \
                             g_0_xxzzzzz_0_yyyzzz_1,  \
                             g_0_xxzzzzz_0_yyzzz_1,   \
                             g_0_xxzzzzz_0_yyzzzz_0,  \
                             g_0_xxzzzzz_0_yyzzzz_1,  \
                             g_0_xxzzzzz_0_yzzzz_1,   \
                             g_0_xxzzzzz_0_yzzzzz_0,  \
                             g_0_xxzzzzz_0_yzzzzz_1,  \
                             g_0_xxzzzzz_0_zzzzz_1,   \
                             g_0_xxzzzzz_0_zzzzzz_0,  \
                             g_0_xxzzzzz_0_zzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_xxxxxx_0[i] = g_0_xxzzzzz_0_xxxxxx_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxy_0[i] = g_0_xxzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxxy_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxxz_0[i] = g_0_xxzzzzz_0_xxxxxz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxyy_0[i] =
            2.0 * g_0_xxzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyy_0[i] * pb_y + g_0_xxzzzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxyz_0[i] = g_0_xxzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxyz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxzz_0[i] = g_0_xxzzzzz_0_xxxxzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_xxzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyy_0[i] * pb_y + g_0_xxzzzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyyz_0[i] =
            2.0 * g_0_xxzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyyz_0[i] * pb_y + g_0_xxzzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyzz_0[i] = g_0_xxzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxzzz_0[i] = g_0_xxzzzzz_0_xxxzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyyy_0[i] =
            4.0 * g_0_xxzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyy_0[i] * pb_y + g_0_xxzzzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyyz_0[i] =
            3.0 * g_0_xxzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyyz_0[i] * pb_y + g_0_xxzzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_xxzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyzz_0[i] * pb_y + g_0_xxzzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyzzz_0[i] = g_0_xxzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxzzzz_0[i] = g_0_xxzzzzz_0_xxzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyyy_0[i] =
            5.0 * g_0_xxzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyy_0[i] * pb_y + g_0_xxzzzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyyz_0[i] =
            4.0 * g_0_xxzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyyz_0[i] * pb_y + g_0_xxzzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyzz_0[i] =
            3.0 * g_0_xxzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyzz_0[i] * pb_y + g_0_xxzzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyzzz_0[i] =
            2.0 * g_0_xxzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyzzz_0[i] * pb_y + g_0_xxzzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyzzzz_0[i] = g_0_xxzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xzzzzz_0[i] = g_0_xxzzzzz_0_xzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyyy_0[i] =
            6.0 * g_0_xxzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyyy_0[i] * pb_y + g_0_xxzzzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyyz_0[i] =
            5.0 * g_0_xxzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyyz_0[i] * pb_y + g_0_xxzzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyzz_0[i] =
            4.0 * g_0_xxzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyzz_0[i] * pb_y + g_0_xxzzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyzzz_0[i] =
            3.0 * g_0_xxzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyzzz_0[i] * pb_y + g_0_xxzzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyzzzz_0[i] =
            2.0 * g_0_xxzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyzzzz_0[i] * pb_y + g_0_xxzzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yzzzzz_0[i] = g_0_xxzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_zzzzzz_0[i] = g_0_xxzzzzz_0_zzzzzz_0[i] * pb_y + g_0_xxzzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 756-784 components of targeted buffer : SLSI

    auto g_0_xxzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 756);

    auto g_0_xxzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 757);

    auto g_0_xxzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 758);

    auto g_0_xxzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 759);

    auto g_0_xxzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 760);

    auto g_0_xxzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 761);

    auto g_0_xxzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 762);

    auto g_0_xxzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 763);

    auto g_0_xxzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 764);

    auto g_0_xxzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 765);

    auto g_0_xxzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 766);

    auto g_0_xxzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 767);

    auto g_0_xxzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 768);

    auto g_0_xxzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 769);

    auto g_0_xxzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 770);

    auto g_0_xxzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 771);

    auto g_0_xxzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 772);

    auto g_0_xxzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 773);

    auto g_0_xxzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 774);

    auto g_0_xxzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 775);

    auto g_0_xxzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 776);

    auto g_0_xxzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 777);

    auto g_0_xxzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 778);

    auto g_0_xxzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 779);

    auto g_0_xxzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 780);

    auto g_0_xxzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 781);

    auto g_0_xxzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 782);

    auto g_0_xxzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 783);

#pragma omp simd aligned(g_0_xxzzzz_0_xxxxxx_0,       \
                             g_0_xxzzzz_0_xxxxxx_1,   \
                             g_0_xxzzzz_0_xxxxxy_0,   \
                             g_0_xxzzzz_0_xxxxxy_1,   \
                             g_0_xxzzzz_0_xxxxyy_0,   \
                             g_0_xxzzzz_0_xxxxyy_1,   \
                             g_0_xxzzzz_0_xxxyyy_0,   \
                             g_0_xxzzzz_0_xxxyyy_1,   \
                             g_0_xxzzzz_0_xxyyyy_0,   \
                             g_0_xxzzzz_0_xxyyyy_1,   \
                             g_0_xxzzzz_0_xyyyyy_0,   \
                             g_0_xxzzzz_0_xyyyyy_1,   \
                             g_0_xxzzzzz_0_xxxxxx_0,  \
                             g_0_xxzzzzz_0_xxxxxx_1,  \
                             g_0_xxzzzzz_0_xxxxxy_0,  \
                             g_0_xxzzzzz_0_xxxxxy_1,  \
                             g_0_xxzzzzz_0_xxxxyy_0,  \
                             g_0_xxzzzzz_0_xxxxyy_1,  \
                             g_0_xxzzzzz_0_xxxyyy_0,  \
                             g_0_xxzzzzz_0_xxxyyy_1,  \
                             g_0_xxzzzzz_0_xxyyyy_0,  \
                             g_0_xxzzzzz_0_xxyyyy_1,  \
                             g_0_xxzzzzz_0_xyyyyy_0,  \
                             g_0_xxzzzzz_0_xyyyyy_1,  \
                             g_0_xxzzzzzz_0_xxxxxx_0, \
                             g_0_xxzzzzzz_0_xxxxxy_0, \
                             g_0_xxzzzzzz_0_xxxxxz_0, \
                             g_0_xxzzzzzz_0_xxxxyy_0, \
                             g_0_xxzzzzzz_0_xxxxyz_0, \
                             g_0_xxzzzzzz_0_xxxxzz_0, \
                             g_0_xxzzzzzz_0_xxxyyy_0, \
                             g_0_xxzzzzzz_0_xxxyyz_0, \
                             g_0_xxzzzzzz_0_xxxyzz_0, \
                             g_0_xxzzzzzz_0_xxxzzz_0, \
                             g_0_xxzzzzzz_0_xxyyyy_0, \
                             g_0_xxzzzzzz_0_xxyyyz_0, \
                             g_0_xxzzzzzz_0_xxyyzz_0, \
                             g_0_xxzzzzzz_0_xxyzzz_0, \
                             g_0_xxzzzzzz_0_xxzzzz_0, \
                             g_0_xxzzzzzz_0_xyyyyy_0, \
                             g_0_xxzzzzzz_0_xyyyyz_0, \
                             g_0_xxzzzzzz_0_xyyyzz_0, \
                             g_0_xxzzzzzz_0_xyyzzz_0, \
                             g_0_xxzzzzzz_0_xyzzzz_0, \
                             g_0_xxzzzzzz_0_xzzzzz_0, \
                             g_0_xxzzzzzz_0_yyyyyy_0, \
                             g_0_xxzzzzzz_0_yyyyyz_0, \
                             g_0_xxzzzzzz_0_yyyyzz_0, \
                             g_0_xxzzzzzz_0_yyyzzz_0, \
                             g_0_xxzzzzzz_0_yyzzzz_0, \
                             g_0_xxzzzzzz_0_yzzzzz_0, \
                             g_0_xxzzzzzz_0_zzzzzz_0, \
                             g_0_xzzzzzz_0_xxxxxz_0,  \
                             g_0_xzzzzzz_0_xxxxxz_1,  \
                             g_0_xzzzzzz_0_xxxxyz_0,  \
                             g_0_xzzzzzz_0_xxxxyz_1,  \
                             g_0_xzzzzzz_0_xxxxz_1,   \
                             g_0_xzzzzzz_0_xxxxzz_0,  \
                             g_0_xzzzzzz_0_xxxxzz_1,  \
                             g_0_xzzzzzz_0_xxxyyz_0,  \
                             g_0_xzzzzzz_0_xxxyyz_1,  \
                             g_0_xzzzzzz_0_xxxyz_1,   \
                             g_0_xzzzzzz_0_xxxyzz_0,  \
                             g_0_xzzzzzz_0_xxxyzz_1,  \
                             g_0_xzzzzzz_0_xxxzz_1,   \
                             g_0_xzzzzzz_0_xxxzzz_0,  \
                             g_0_xzzzzzz_0_xxxzzz_1,  \
                             g_0_xzzzzzz_0_xxyyyz_0,  \
                             g_0_xzzzzzz_0_xxyyyz_1,  \
                             g_0_xzzzzzz_0_xxyyz_1,   \
                             g_0_xzzzzzz_0_xxyyzz_0,  \
                             g_0_xzzzzzz_0_xxyyzz_1,  \
                             g_0_xzzzzzz_0_xxyzz_1,   \
                             g_0_xzzzzzz_0_xxyzzz_0,  \
                             g_0_xzzzzzz_0_xxyzzz_1,  \
                             g_0_xzzzzzz_0_xxzzz_1,   \
                             g_0_xzzzzzz_0_xxzzzz_0,  \
                             g_0_xzzzzzz_0_xxzzzz_1,  \
                             g_0_xzzzzzz_0_xyyyyz_0,  \
                             g_0_xzzzzzz_0_xyyyyz_1,  \
                             g_0_xzzzzzz_0_xyyyz_1,   \
                             g_0_xzzzzzz_0_xyyyzz_0,  \
                             g_0_xzzzzzz_0_xyyyzz_1,  \
                             g_0_xzzzzzz_0_xyyzz_1,   \
                             g_0_xzzzzzz_0_xyyzzz_0,  \
                             g_0_xzzzzzz_0_xyyzzz_1,  \
                             g_0_xzzzzzz_0_xyzzz_1,   \
                             g_0_xzzzzzz_0_xyzzzz_0,  \
                             g_0_xzzzzzz_0_xyzzzz_1,  \
                             g_0_xzzzzzz_0_xzzzz_1,   \
                             g_0_xzzzzzz_0_xzzzzz_0,  \
                             g_0_xzzzzzz_0_xzzzzz_1,  \
                             g_0_xzzzzzz_0_yyyyyy_0,  \
                             g_0_xzzzzzz_0_yyyyyy_1,  \
                             g_0_xzzzzzz_0_yyyyyz_0,  \
                             g_0_xzzzzzz_0_yyyyyz_1,  \
                             g_0_xzzzzzz_0_yyyyz_1,   \
                             g_0_xzzzzzz_0_yyyyzz_0,  \
                             g_0_xzzzzzz_0_yyyyzz_1,  \
                             g_0_xzzzzzz_0_yyyzz_1,   \
                             g_0_xzzzzzz_0_yyyzzz_0,  \
                             g_0_xzzzzzz_0_yyyzzz_1,  \
                             g_0_xzzzzzz_0_yyzzz_1,   \
                             g_0_xzzzzzz_0_yyzzzz_0,  \
                             g_0_xzzzzzz_0_yyzzzz_1,  \
                             g_0_xzzzzzz_0_yzzzz_1,   \
                             g_0_xzzzzzz_0_yzzzzz_0,  \
                             g_0_xzzzzzz_0_yzzzzz_1,  \
                             g_0_xzzzzzz_0_zzzzz_1,   \
                             g_0_xzzzzzz_0_zzzzzz_0,  \
                             g_0_xzzzzzz_0_zzzzzz_1,  \
                             g_0_zzzzzz_0_xxxxxz_0,   \
                             g_0_zzzzzz_0_xxxxxz_1,   \
                             g_0_zzzzzz_0_xxxxyz_0,   \
                             g_0_zzzzzz_0_xxxxyz_1,   \
                             g_0_zzzzzz_0_xxxxzz_0,   \
                             g_0_zzzzzz_0_xxxxzz_1,   \
                             g_0_zzzzzz_0_xxxyyz_0,   \
                             g_0_zzzzzz_0_xxxyyz_1,   \
                             g_0_zzzzzz_0_xxxyzz_0,   \
                             g_0_zzzzzz_0_xxxyzz_1,   \
                             g_0_zzzzzz_0_xxxzzz_0,   \
                             g_0_zzzzzz_0_xxxzzz_1,   \
                             g_0_zzzzzz_0_xxyyyz_0,   \
                             g_0_zzzzzz_0_xxyyyz_1,   \
                             g_0_zzzzzz_0_xxyyzz_0,   \
                             g_0_zzzzzz_0_xxyyzz_1,   \
                             g_0_zzzzzz_0_xxyzzz_0,   \
                             g_0_zzzzzz_0_xxyzzz_1,   \
                             g_0_zzzzzz_0_xxzzzz_0,   \
                             g_0_zzzzzz_0_xxzzzz_1,   \
                             g_0_zzzzzz_0_xyyyyz_0,   \
                             g_0_zzzzzz_0_xyyyyz_1,   \
                             g_0_zzzzzz_0_xyyyzz_0,   \
                             g_0_zzzzzz_0_xyyyzz_1,   \
                             g_0_zzzzzz_0_xyyzzz_0,   \
                             g_0_zzzzzz_0_xyyzzz_1,   \
                             g_0_zzzzzz_0_xyzzzz_0,   \
                             g_0_zzzzzz_0_xyzzzz_1,   \
                             g_0_zzzzzz_0_xzzzzz_0,   \
                             g_0_zzzzzz_0_xzzzzz_1,   \
                             g_0_zzzzzz_0_yyyyyy_0,   \
                             g_0_zzzzzz_0_yyyyyy_1,   \
                             g_0_zzzzzz_0_yyyyyz_0,   \
                             g_0_zzzzzz_0_yyyyyz_1,   \
                             g_0_zzzzzz_0_yyyyzz_0,   \
                             g_0_zzzzzz_0_yyyyzz_1,   \
                             g_0_zzzzzz_0_yyyzzz_0,   \
                             g_0_zzzzzz_0_yyyzzz_1,   \
                             g_0_zzzzzz_0_yyzzzz_0,   \
                             g_0_zzzzzz_0_yyzzzz_1,   \
                             g_0_zzzzzz_0_yzzzzz_0,   \
                             g_0_zzzzzz_0_yzzzzz_1,   \
                             g_0_zzzzzz_0_zzzzzz_0,   \
                             g_0_zzzzzz_0_zzzzzz_1,   \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_xxxxxx_0[i] = 5.0 * g_0_xxzzzz_0_xxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_xxxxxx_0[i] * pb_z + g_0_xxzzzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxxy_0[i] = 5.0 * g_0_xxzzzz_0_xxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_xxxxxy_0[i] * pb_z + g_0_xxzzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxxz_0[i] = g_0_zzzzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     5.0 * g_0_xzzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxxz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxyy_0[i] = 5.0 * g_0_xxzzzz_0_xxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_xxxxyy_0[i] * pb_z + g_0_xxzzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxyz_0[i] = g_0_zzzzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xzzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxyz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxxzz_0[i] = g_0_zzzzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_xzzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_xxxyyy_0[i] * pb_z + g_0_xxzzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxyyz_0[i] = g_0_zzzzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxyyz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxyzz_0[i] = g_0_zzzzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxyzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxzzz_0[i] = g_0_zzzzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_xzzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_xxyyyy_0[i] * pb_z + g_0_xxzzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxyyyz_0[i] = g_0_zzzzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyyyz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyyzz_0[i] = g_0_zzzzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyyzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyzzz_0[i] = g_0_zzzzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxzzzz_0[i] = g_0_zzzzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_xzzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxzzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_xxzzzzz_0_xyyyyy_0[i] * pb_z + g_0_xxzzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xyyyyz_0[i] = g_0_zzzzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyz_1[i] * fi_abcd_0 +
                                     g_0_xzzzzzz_0_xyyyyz_0[i] * pb_x + g_0_xzzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyyzz_0[i] = g_0_zzzzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzzzz_0_xyyyzz_0[i] * pb_x + g_0_xzzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyzzz_0[i] = g_0_zzzzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzzzz_0_xyyzzz_0[i] * pb_x + g_0_xzzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyzzzz_0[i] = g_0_zzzzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzzzz_0_xyzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xzzzzz_0[i] = g_0_zzzzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                     g_0_xzzzzzz_0_xzzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyyy_0[i] = g_0_zzzzzz_0_yyyyyy_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyyy_0[i] * pb_x +
                                     g_0_xzzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyyz_0[i] = g_0_zzzzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyyz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyzz_0[i] = g_0_zzzzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyzzz_0[i] = g_0_zzzzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyzzzz_0[i] = g_0_zzzzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yzzzzz_0[i] = g_0_zzzzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_zzzzzz_0[i] = g_0_zzzzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzzzz_0[i] * pb_x +
                                     g_0_xzzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 784-812 components of targeted buffer : SLSI

    auto g_0_xyyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 784);

    auto g_0_xyyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 785);

    auto g_0_xyyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 786);

    auto g_0_xyyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 787);

    auto g_0_xyyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 788);

    auto g_0_xyyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 789);

    auto g_0_xyyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 790);

    auto g_0_xyyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 791);

    auto g_0_xyyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 792);

    auto g_0_xyyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 793);

    auto g_0_xyyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 794);

    auto g_0_xyyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 795);

    auto g_0_xyyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 796);

    auto g_0_xyyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 797);

    auto g_0_xyyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 798);

    auto g_0_xyyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 799);

    auto g_0_xyyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 800);

    auto g_0_xyyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 801);

    auto g_0_xyyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 802);

    auto g_0_xyyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 803);

    auto g_0_xyyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 804);

    auto g_0_xyyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 805);

    auto g_0_xyyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 806);

    auto g_0_xyyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 807);

    auto g_0_xyyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 808);

    auto g_0_xyyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 809);

    auto g_0_xyyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 810);

    auto g_0_xyyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 811);

#pragma omp simd aligned(g_0_xyyyyyyy_0_xxxxxx_0,     \
                             g_0_xyyyyyyy_0_xxxxxy_0, \
                             g_0_xyyyyyyy_0_xxxxxz_0, \
                             g_0_xyyyyyyy_0_xxxxyy_0, \
                             g_0_xyyyyyyy_0_xxxxyz_0, \
                             g_0_xyyyyyyy_0_xxxxzz_0, \
                             g_0_xyyyyyyy_0_xxxyyy_0, \
                             g_0_xyyyyyyy_0_xxxyyz_0, \
                             g_0_xyyyyyyy_0_xxxyzz_0, \
                             g_0_xyyyyyyy_0_xxxzzz_0, \
                             g_0_xyyyyyyy_0_xxyyyy_0, \
                             g_0_xyyyyyyy_0_xxyyyz_0, \
                             g_0_xyyyyyyy_0_xxyyzz_0, \
                             g_0_xyyyyyyy_0_xxyzzz_0, \
                             g_0_xyyyyyyy_0_xxzzzz_0, \
                             g_0_xyyyyyyy_0_xyyyyy_0, \
                             g_0_xyyyyyyy_0_xyyyyz_0, \
                             g_0_xyyyyyyy_0_xyyyzz_0, \
                             g_0_xyyyyyyy_0_xyyzzz_0, \
                             g_0_xyyyyyyy_0_xyzzzz_0, \
                             g_0_xyyyyyyy_0_xzzzzz_0, \
                             g_0_xyyyyyyy_0_yyyyyy_0, \
                             g_0_xyyyyyyy_0_yyyyyz_0, \
                             g_0_xyyyyyyy_0_yyyyzz_0, \
                             g_0_xyyyyyyy_0_yyyzzz_0, \
                             g_0_xyyyyyyy_0_yyzzzz_0, \
                             g_0_xyyyyyyy_0_yzzzzz_0, \
                             g_0_xyyyyyyy_0_zzzzzz_0, \
                             g_0_yyyyyyy_0_xxxxx_1,   \
                             g_0_yyyyyyy_0_xxxxxx_0,  \
                             g_0_yyyyyyy_0_xxxxxx_1,  \
                             g_0_yyyyyyy_0_xxxxxy_0,  \
                             g_0_yyyyyyy_0_xxxxxy_1,  \
                             g_0_yyyyyyy_0_xxxxxz_0,  \
                             g_0_yyyyyyy_0_xxxxxz_1,  \
                             g_0_yyyyyyy_0_xxxxy_1,   \
                             g_0_yyyyyyy_0_xxxxyy_0,  \
                             g_0_yyyyyyy_0_xxxxyy_1,  \
                             g_0_yyyyyyy_0_xxxxyz_0,  \
                             g_0_yyyyyyy_0_xxxxyz_1,  \
                             g_0_yyyyyyy_0_xxxxz_1,   \
                             g_0_yyyyyyy_0_xxxxzz_0,  \
                             g_0_yyyyyyy_0_xxxxzz_1,  \
                             g_0_yyyyyyy_0_xxxyy_1,   \
                             g_0_yyyyyyy_0_xxxyyy_0,  \
                             g_0_yyyyyyy_0_xxxyyy_1,  \
                             g_0_yyyyyyy_0_xxxyyz_0,  \
                             g_0_yyyyyyy_0_xxxyyz_1,  \
                             g_0_yyyyyyy_0_xxxyz_1,   \
                             g_0_yyyyyyy_0_xxxyzz_0,  \
                             g_0_yyyyyyy_0_xxxyzz_1,  \
                             g_0_yyyyyyy_0_xxxzz_1,   \
                             g_0_yyyyyyy_0_xxxzzz_0,  \
                             g_0_yyyyyyy_0_xxxzzz_1,  \
                             g_0_yyyyyyy_0_xxyyy_1,   \
                             g_0_yyyyyyy_0_xxyyyy_0,  \
                             g_0_yyyyyyy_0_xxyyyy_1,  \
                             g_0_yyyyyyy_0_xxyyyz_0,  \
                             g_0_yyyyyyy_0_xxyyyz_1,  \
                             g_0_yyyyyyy_0_xxyyz_1,   \
                             g_0_yyyyyyy_0_xxyyzz_0,  \
                             g_0_yyyyyyy_0_xxyyzz_1,  \
                             g_0_yyyyyyy_0_xxyzz_1,   \
                             g_0_yyyyyyy_0_xxyzzz_0,  \
                             g_0_yyyyyyy_0_xxyzzz_1,  \
                             g_0_yyyyyyy_0_xxzzz_1,   \
                             g_0_yyyyyyy_0_xxzzzz_0,  \
                             g_0_yyyyyyy_0_xxzzzz_1,  \
                             g_0_yyyyyyy_0_xyyyy_1,   \
                             g_0_yyyyyyy_0_xyyyyy_0,  \
                             g_0_yyyyyyy_0_xyyyyy_1,  \
                             g_0_yyyyyyy_0_xyyyyz_0,  \
                             g_0_yyyyyyy_0_xyyyyz_1,  \
                             g_0_yyyyyyy_0_xyyyz_1,   \
                             g_0_yyyyyyy_0_xyyyzz_0,  \
                             g_0_yyyyyyy_0_xyyyzz_1,  \
                             g_0_yyyyyyy_0_xyyzz_1,   \
                             g_0_yyyyyyy_0_xyyzzz_0,  \
                             g_0_yyyyyyy_0_xyyzzz_1,  \
                             g_0_yyyyyyy_0_xyzzz_1,   \
                             g_0_yyyyyyy_0_xyzzzz_0,  \
                             g_0_yyyyyyy_0_xyzzzz_1,  \
                             g_0_yyyyyyy_0_xzzzz_1,   \
                             g_0_yyyyyyy_0_xzzzzz_0,  \
                             g_0_yyyyyyy_0_xzzzzz_1,  \
                             g_0_yyyyyyy_0_yyyyy_1,   \
                             g_0_yyyyyyy_0_yyyyyy_0,  \
                             g_0_yyyyyyy_0_yyyyyy_1,  \
                             g_0_yyyyyyy_0_yyyyyz_0,  \
                             g_0_yyyyyyy_0_yyyyyz_1,  \
                             g_0_yyyyyyy_0_yyyyz_1,   \
                             g_0_yyyyyyy_0_yyyyzz_0,  \
                             g_0_yyyyyyy_0_yyyyzz_1,  \
                             g_0_yyyyyyy_0_yyyzz_1,   \
                             g_0_yyyyyyy_0_yyyzzz_0,  \
                             g_0_yyyyyyy_0_yyyzzz_1,  \
                             g_0_yyyyyyy_0_yyzzz_1,   \
                             g_0_yyyyyyy_0_yyzzzz_0,  \
                             g_0_yyyyyyy_0_yyzzzz_1,  \
                             g_0_yyyyyyy_0_yzzzz_1,   \
                             g_0_yyyyyyy_0_yzzzzz_0,  \
                             g_0_yyyyyyy_0_yzzzzz_1,  \
                             g_0_yyyyyyy_0_zzzzz_1,   \
                             g_0_yyyyyyy_0_zzzzzz_0,  \
                             g_0_yyyyyyy_0_zzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_xxxxxx_0[i] =
            6.0 * g_0_yyyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxx_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxy_0[i] =
            5.0 * g_0_yyyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxy_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxxz_0[i] =
            5.0 * g_0_yyyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxyy_0[i] =
            4.0 * g_0_yyyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyy_0[i] * pb_x + g_0_yyyyyyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxyz_0[i] =
            4.0 * g_0_yyyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxzz_0[i] =
            4.0 * g_0_yyyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyyy_0[i] =
            3.0 * g_0_yyyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyy_0[i] * pb_x + g_0_yyyyyyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyyz_0[i] =
            3.0 * g_0_yyyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyz_0[i] * pb_x + g_0_yyyyyyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyzz_0[i] =
            3.0 * g_0_yyyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxzzz_0[i] =
            3.0 * g_0_yyyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyyy_0[i] =
            2.0 * g_0_yyyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyy_0[i] * pb_x + g_0_yyyyyyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyyz_0[i] =
            2.0 * g_0_yyyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyz_0[i] * pb_x + g_0_yyyyyyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyzz_0[i] * pb_x + g_0_yyyyyyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyzzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxzzzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyyy_0[i] = g_0_yyyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyy_0[i] * pb_x + g_0_yyyyyyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyyz_0[i] = g_0_yyyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyz_0[i] * pb_x + g_0_yyyyyyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyzz_0[i] = g_0_yyyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyzz_0[i] * pb_x + g_0_yyyyyyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyzzz_0[i] = g_0_yyyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzzz_0[i] * pb_x + g_0_yyyyyyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyzzzz_0[i] = g_0_yyyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xzzzzz_0[i] = g_0_yyyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyyy_0[i] = g_0_yyyyyyy_0_yyyyyy_0[i] * pb_x + g_0_yyyyyyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyyz_0[i] = g_0_yyyyyyy_0_yyyyyz_0[i] * pb_x + g_0_yyyyyyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyzz_0[i] = g_0_yyyyyyy_0_yyyyzz_0[i] * pb_x + g_0_yyyyyyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyzzz_0[i] = g_0_yyyyyyy_0_yyyzzz_0[i] * pb_x + g_0_yyyyyyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyzzzz_0[i] = g_0_yyyyyyy_0_yyzzzz_0[i] * pb_x + g_0_yyyyyyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yzzzzz_0[i] = g_0_yyyyyyy_0_yzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_zzzzzz_0[i] = g_0_yyyyyyy_0_zzzzzz_0[i] * pb_x + g_0_yyyyyyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 812-840 components of targeted buffer : SLSI

    auto g_0_xyyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 812);

    auto g_0_xyyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 813);

    auto g_0_xyyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 814);

    auto g_0_xyyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 815);

    auto g_0_xyyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 816);

    auto g_0_xyyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 817);

    auto g_0_xyyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 818);

    auto g_0_xyyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 819);

    auto g_0_xyyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 820);

    auto g_0_xyyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 821);

    auto g_0_xyyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 822);

    auto g_0_xyyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 823);

    auto g_0_xyyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 824);

    auto g_0_xyyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 825);

    auto g_0_xyyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 826);

    auto g_0_xyyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 827);

    auto g_0_xyyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 828);

    auto g_0_xyyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 829);

    auto g_0_xyyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 830);

    auto g_0_xyyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 831);

    auto g_0_xyyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 832);

    auto g_0_xyyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 833);

    auto g_0_xyyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 834);

    auto g_0_xyyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 835);

    auto g_0_xyyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 836);

    auto g_0_xyyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 837);

    auto g_0_xyyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 838);

    auto g_0_xyyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 839);

#pragma omp simd aligned(g_0_xyyyyyy_0_xxxxxx_0,      \
                             g_0_xyyyyyy_0_xxxxxx_1,  \
                             g_0_xyyyyyy_0_xxxxxy_0,  \
                             g_0_xyyyyyy_0_xxxxxy_1,  \
                             g_0_xyyyyyy_0_xxxxyy_0,  \
                             g_0_xyyyyyy_0_xxxxyy_1,  \
                             g_0_xyyyyyy_0_xxxyyy_0,  \
                             g_0_xyyyyyy_0_xxxyyy_1,  \
                             g_0_xyyyyyy_0_xxyyyy_0,  \
                             g_0_xyyyyyy_0_xxyyyy_1,  \
                             g_0_xyyyyyy_0_xyyyyy_0,  \
                             g_0_xyyyyyy_0_xyyyyy_1,  \
                             g_0_xyyyyyyz_0_xxxxxx_0, \
                             g_0_xyyyyyyz_0_xxxxxy_0, \
                             g_0_xyyyyyyz_0_xxxxxz_0, \
                             g_0_xyyyyyyz_0_xxxxyy_0, \
                             g_0_xyyyyyyz_0_xxxxyz_0, \
                             g_0_xyyyyyyz_0_xxxxzz_0, \
                             g_0_xyyyyyyz_0_xxxyyy_0, \
                             g_0_xyyyyyyz_0_xxxyyz_0, \
                             g_0_xyyyyyyz_0_xxxyzz_0, \
                             g_0_xyyyyyyz_0_xxxzzz_0, \
                             g_0_xyyyyyyz_0_xxyyyy_0, \
                             g_0_xyyyyyyz_0_xxyyyz_0, \
                             g_0_xyyyyyyz_0_xxyyzz_0, \
                             g_0_xyyyyyyz_0_xxyzzz_0, \
                             g_0_xyyyyyyz_0_xxzzzz_0, \
                             g_0_xyyyyyyz_0_xyyyyy_0, \
                             g_0_xyyyyyyz_0_xyyyyz_0, \
                             g_0_xyyyyyyz_0_xyyyzz_0, \
                             g_0_xyyyyyyz_0_xyyzzz_0, \
                             g_0_xyyyyyyz_0_xyzzzz_0, \
                             g_0_xyyyyyyz_0_xzzzzz_0, \
                             g_0_xyyyyyyz_0_yyyyyy_0, \
                             g_0_xyyyyyyz_0_yyyyyz_0, \
                             g_0_xyyyyyyz_0_yyyyzz_0, \
                             g_0_xyyyyyyz_0_yyyzzz_0, \
                             g_0_xyyyyyyz_0_yyzzzz_0, \
                             g_0_xyyyyyyz_0_yzzzzz_0, \
                             g_0_xyyyyyyz_0_zzzzzz_0, \
                             g_0_yyyyyyz_0_xxxxxz_0,  \
                             g_0_yyyyyyz_0_xxxxxz_1,  \
                             g_0_yyyyyyz_0_xxxxyz_0,  \
                             g_0_yyyyyyz_0_xxxxyz_1,  \
                             g_0_yyyyyyz_0_xxxxz_1,   \
                             g_0_yyyyyyz_0_xxxxzz_0,  \
                             g_0_yyyyyyz_0_xxxxzz_1,  \
                             g_0_yyyyyyz_0_xxxyyz_0,  \
                             g_0_yyyyyyz_0_xxxyyz_1,  \
                             g_0_yyyyyyz_0_xxxyz_1,   \
                             g_0_yyyyyyz_0_xxxyzz_0,  \
                             g_0_yyyyyyz_0_xxxyzz_1,  \
                             g_0_yyyyyyz_0_xxxzz_1,   \
                             g_0_yyyyyyz_0_xxxzzz_0,  \
                             g_0_yyyyyyz_0_xxxzzz_1,  \
                             g_0_yyyyyyz_0_xxyyyz_0,  \
                             g_0_yyyyyyz_0_xxyyyz_1,  \
                             g_0_yyyyyyz_0_xxyyz_1,   \
                             g_0_yyyyyyz_0_xxyyzz_0,  \
                             g_0_yyyyyyz_0_xxyyzz_1,  \
                             g_0_yyyyyyz_0_xxyzz_1,   \
                             g_0_yyyyyyz_0_xxyzzz_0,  \
                             g_0_yyyyyyz_0_xxyzzz_1,  \
                             g_0_yyyyyyz_0_xxzzz_1,   \
                             g_0_yyyyyyz_0_xxzzzz_0,  \
                             g_0_yyyyyyz_0_xxzzzz_1,  \
                             g_0_yyyyyyz_0_xyyyyz_0,  \
                             g_0_yyyyyyz_0_xyyyyz_1,  \
                             g_0_yyyyyyz_0_xyyyz_1,   \
                             g_0_yyyyyyz_0_xyyyzz_0,  \
                             g_0_yyyyyyz_0_xyyyzz_1,  \
                             g_0_yyyyyyz_0_xyyzz_1,   \
                             g_0_yyyyyyz_0_xyyzzz_0,  \
                             g_0_yyyyyyz_0_xyyzzz_1,  \
                             g_0_yyyyyyz_0_xyzzz_1,   \
                             g_0_yyyyyyz_0_xyzzzz_0,  \
                             g_0_yyyyyyz_0_xyzzzz_1,  \
                             g_0_yyyyyyz_0_xzzzz_1,   \
                             g_0_yyyyyyz_0_xzzzzz_0,  \
                             g_0_yyyyyyz_0_xzzzzz_1,  \
                             g_0_yyyyyyz_0_yyyyyy_0,  \
                             g_0_yyyyyyz_0_yyyyyy_1,  \
                             g_0_yyyyyyz_0_yyyyyz_0,  \
                             g_0_yyyyyyz_0_yyyyyz_1,  \
                             g_0_yyyyyyz_0_yyyyz_1,   \
                             g_0_yyyyyyz_0_yyyyzz_0,  \
                             g_0_yyyyyyz_0_yyyyzz_1,  \
                             g_0_yyyyyyz_0_yyyzz_1,   \
                             g_0_yyyyyyz_0_yyyzzz_0,  \
                             g_0_yyyyyyz_0_yyyzzz_1,  \
                             g_0_yyyyyyz_0_yyzzz_1,   \
                             g_0_yyyyyyz_0_yyzzzz_0,  \
                             g_0_yyyyyyz_0_yyzzzz_1,  \
                             g_0_yyyyyyz_0_yzzzz_1,   \
                             g_0_yyyyyyz_0_yzzzzz_0,  \
                             g_0_yyyyyyz_0_yzzzzz_1,  \
                             g_0_yyyyyyz_0_zzzzz_1,   \
                             g_0_yyyyyyz_0_zzzzzz_0,  \
                             g_0_yyyyyyz_0_zzzzzz_1,  \
                             wp_x,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyz_0_xxxxxx_0[i] = g_0_xyyyyyy_0_xxxxxx_0[i] * pb_z + g_0_xyyyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxxy_0[i] = g_0_xyyyyyy_0_xxxxxy_0[i] * pb_z + g_0_xyyyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxxz_0[i] =
            5.0 * g_0_yyyyyyz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxxz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxyy_0[i] = g_0_xyyyyyy_0_xxxxyy_0[i] * pb_z + g_0_xyyyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxyz_0[i] =
            4.0 * g_0_yyyyyyz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxyz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxxzz_0[i] =
            4.0 * g_0_yyyyyyz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxyyy_0[i] = g_0_xyyyyyy_0_xxxyyy_0[i] * pb_z + g_0_xyyyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxyyz_0[i] =
            3.0 * g_0_yyyyyyz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxyyz_0[i] * pb_x + g_0_yyyyyyz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxyzz_0[i] =
            3.0 * g_0_yyyyyyz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxyzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxzzz_0[i] =
            3.0 * g_0_yyyyyyz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyyyy_0[i] = g_0_xyyyyyy_0_xxyyyy_0[i] * pb_z + g_0_xyyyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxyyyz_0[i] =
            2.0 * g_0_yyyyyyz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyyyz_0[i] * pb_x + g_0_yyyyyyz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyyzz_0[i] =
            2.0 * g_0_yyyyyyz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyyzz_0[i] * pb_x + g_0_yyyyyyz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyzzz_0[i] =
            2.0 * g_0_yyyyyyz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxzzzz_0[i] =
            2.0 * g_0_yyyyyyz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyyyy_0[i] = g_0_xyyyyyy_0_xyyyyy_0[i] * pb_z + g_0_xyyyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xyyyyz_0[i] = g_0_yyyyyyz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyyyz_0[i] * pb_x + g_0_yyyyyyz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyyzz_0[i] = g_0_yyyyyyz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyyzz_0[i] * pb_x + g_0_yyyyyyz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyzzz_0[i] = g_0_yyyyyyz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyzzz_0[i] * pb_x + g_0_yyyyyyz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyzzzz_0[i] = g_0_yyyyyyz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xzzzzz_0[i] = g_0_yyyyyyz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyyy_0[i] = g_0_yyyyyyz_0_yyyyyy_0[i] * pb_x + g_0_yyyyyyz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyyz_0[i] = g_0_yyyyyyz_0_yyyyyz_0[i] * pb_x + g_0_yyyyyyz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyzz_0[i] = g_0_yyyyyyz_0_yyyyzz_0[i] * pb_x + g_0_yyyyyyz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyzzz_0[i] = g_0_yyyyyyz_0_yyyzzz_0[i] * pb_x + g_0_yyyyyyz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyzzzz_0[i] = g_0_yyyyyyz_0_yyzzzz_0[i] * pb_x + g_0_yyyyyyz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yzzzzz_0[i] = g_0_yyyyyyz_0_yzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_zzzzzz_0[i] = g_0_yyyyyyz_0_zzzzzz_0[i] * pb_x + g_0_yyyyyyz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 840-868 components of targeted buffer : SLSI

    auto g_0_xyyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 840);

    auto g_0_xyyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 841);

    auto g_0_xyyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 842);

    auto g_0_xyyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 843);

    auto g_0_xyyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 844);

    auto g_0_xyyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 845);

    auto g_0_xyyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 846);

    auto g_0_xyyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 847);

    auto g_0_xyyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 848);

    auto g_0_xyyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 849);

    auto g_0_xyyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 850);

    auto g_0_xyyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 851);

    auto g_0_xyyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 852);

    auto g_0_xyyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 853);

    auto g_0_xyyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 854);

    auto g_0_xyyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 855);

    auto g_0_xyyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 856);

    auto g_0_xyyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 857);

    auto g_0_xyyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 858);

    auto g_0_xyyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 859);

    auto g_0_xyyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 860);

    auto g_0_xyyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 861);

    auto g_0_xyyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 862);

    auto g_0_xyyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 863);

    auto g_0_xyyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 864);

    auto g_0_xyyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 865);

    auto g_0_xyyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 866);

    auto g_0_xyyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 867);

#pragma omp simd aligned(g_0_xyyyyyzz_0_xxxxxx_0,     \
                             g_0_xyyyyyzz_0_xxxxxy_0, \
                             g_0_xyyyyyzz_0_xxxxxz_0, \
                             g_0_xyyyyyzz_0_xxxxyy_0, \
                             g_0_xyyyyyzz_0_xxxxyz_0, \
                             g_0_xyyyyyzz_0_xxxxzz_0, \
                             g_0_xyyyyyzz_0_xxxyyy_0, \
                             g_0_xyyyyyzz_0_xxxyyz_0, \
                             g_0_xyyyyyzz_0_xxxyzz_0, \
                             g_0_xyyyyyzz_0_xxxzzz_0, \
                             g_0_xyyyyyzz_0_xxyyyy_0, \
                             g_0_xyyyyyzz_0_xxyyyz_0, \
                             g_0_xyyyyyzz_0_xxyyzz_0, \
                             g_0_xyyyyyzz_0_xxyzzz_0, \
                             g_0_xyyyyyzz_0_xxzzzz_0, \
                             g_0_xyyyyyzz_0_xyyyyy_0, \
                             g_0_xyyyyyzz_0_xyyyyz_0, \
                             g_0_xyyyyyzz_0_xyyyzz_0, \
                             g_0_xyyyyyzz_0_xyyzzz_0, \
                             g_0_xyyyyyzz_0_xyzzzz_0, \
                             g_0_xyyyyyzz_0_xzzzzz_0, \
                             g_0_xyyyyyzz_0_yyyyyy_0, \
                             g_0_xyyyyyzz_0_yyyyyz_0, \
                             g_0_xyyyyyzz_0_yyyyzz_0, \
                             g_0_xyyyyyzz_0_yyyzzz_0, \
                             g_0_xyyyyyzz_0_yyzzzz_0, \
                             g_0_xyyyyyzz_0_yzzzzz_0, \
                             g_0_xyyyyyzz_0_zzzzzz_0, \
                             g_0_yyyyyzz_0_xxxxx_1,   \
                             g_0_yyyyyzz_0_xxxxxx_0,  \
                             g_0_yyyyyzz_0_xxxxxx_1,  \
                             g_0_yyyyyzz_0_xxxxxy_0,  \
                             g_0_yyyyyzz_0_xxxxxy_1,  \
                             g_0_yyyyyzz_0_xxxxxz_0,  \
                             g_0_yyyyyzz_0_xxxxxz_1,  \
                             g_0_yyyyyzz_0_xxxxy_1,   \
                             g_0_yyyyyzz_0_xxxxyy_0,  \
                             g_0_yyyyyzz_0_xxxxyy_1,  \
                             g_0_yyyyyzz_0_xxxxyz_0,  \
                             g_0_yyyyyzz_0_xxxxyz_1,  \
                             g_0_yyyyyzz_0_xxxxz_1,   \
                             g_0_yyyyyzz_0_xxxxzz_0,  \
                             g_0_yyyyyzz_0_xxxxzz_1,  \
                             g_0_yyyyyzz_0_xxxyy_1,   \
                             g_0_yyyyyzz_0_xxxyyy_0,  \
                             g_0_yyyyyzz_0_xxxyyy_1,  \
                             g_0_yyyyyzz_0_xxxyyz_0,  \
                             g_0_yyyyyzz_0_xxxyyz_1,  \
                             g_0_yyyyyzz_0_xxxyz_1,   \
                             g_0_yyyyyzz_0_xxxyzz_0,  \
                             g_0_yyyyyzz_0_xxxyzz_1,  \
                             g_0_yyyyyzz_0_xxxzz_1,   \
                             g_0_yyyyyzz_0_xxxzzz_0,  \
                             g_0_yyyyyzz_0_xxxzzz_1,  \
                             g_0_yyyyyzz_0_xxyyy_1,   \
                             g_0_yyyyyzz_0_xxyyyy_0,  \
                             g_0_yyyyyzz_0_xxyyyy_1,  \
                             g_0_yyyyyzz_0_xxyyyz_0,  \
                             g_0_yyyyyzz_0_xxyyyz_1,  \
                             g_0_yyyyyzz_0_xxyyz_1,   \
                             g_0_yyyyyzz_0_xxyyzz_0,  \
                             g_0_yyyyyzz_0_xxyyzz_1,  \
                             g_0_yyyyyzz_0_xxyzz_1,   \
                             g_0_yyyyyzz_0_xxyzzz_0,  \
                             g_0_yyyyyzz_0_xxyzzz_1,  \
                             g_0_yyyyyzz_0_xxzzz_1,   \
                             g_0_yyyyyzz_0_xxzzzz_0,  \
                             g_0_yyyyyzz_0_xxzzzz_1,  \
                             g_0_yyyyyzz_0_xyyyy_1,   \
                             g_0_yyyyyzz_0_xyyyyy_0,  \
                             g_0_yyyyyzz_0_xyyyyy_1,  \
                             g_0_yyyyyzz_0_xyyyyz_0,  \
                             g_0_yyyyyzz_0_xyyyyz_1,  \
                             g_0_yyyyyzz_0_xyyyz_1,   \
                             g_0_yyyyyzz_0_xyyyzz_0,  \
                             g_0_yyyyyzz_0_xyyyzz_1,  \
                             g_0_yyyyyzz_0_xyyzz_1,   \
                             g_0_yyyyyzz_0_xyyzzz_0,  \
                             g_0_yyyyyzz_0_xyyzzz_1,  \
                             g_0_yyyyyzz_0_xyzzz_1,   \
                             g_0_yyyyyzz_0_xyzzzz_0,  \
                             g_0_yyyyyzz_0_xyzzzz_1,  \
                             g_0_yyyyyzz_0_xzzzz_1,   \
                             g_0_yyyyyzz_0_xzzzzz_0,  \
                             g_0_yyyyyzz_0_xzzzzz_1,  \
                             g_0_yyyyyzz_0_yyyyy_1,   \
                             g_0_yyyyyzz_0_yyyyyy_0,  \
                             g_0_yyyyyzz_0_yyyyyy_1,  \
                             g_0_yyyyyzz_0_yyyyyz_0,  \
                             g_0_yyyyyzz_0_yyyyyz_1,  \
                             g_0_yyyyyzz_0_yyyyz_1,   \
                             g_0_yyyyyzz_0_yyyyzz_0,  \
                             g_0_yyyyyzz_0_yyyyzz_1,  \
                             g_0_yyyyyzz_0_yyyzz_1,   \
                             g_0_yyyyyzz_0_yyyzzz_0,  \
                             g_0_yyyyyzz_0_yyyzzz_1,  \
                             g_0_yyyyyzz_0_yyzzz_1,   \
                             g_0_yyyyyzz_0_yyzzzz_0,  \
                             g_0_yyyyyzz_0_yyzzzz_1,  \
                             g_0_yyyyyzz_0_yzzzz_1,   \
                             g_0_yyyyyzz_0_yzzzzz_0,  \
                             g_0_yyyyyzz_0_yzzzzz_1,  \
                             g_0_yyyyyzz_0_zzzzz_1,   \
                             g_0_yyyyyzz_0_zzzzzz_0,  \
                             g_0_yyyyyzz_0_zzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_xxxxxx_0[i] =
            6.0 * g_0_yyyyyzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxx_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxy_0[i] =
            5.0 * g_0_yyyyyzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxy_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxxz_0[i] =
            5.0 * g_0_yyyyyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxxz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxyy_0[i] =
            4.0 * g_0_yyyyyzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyy_0[i] * pb_x + g_0_yyyyyzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxyz_0[i] =
            4.0 * g_0_yyyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxzz_0[i] =
            4.0 * g_0_yyyyyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyyy_0[i] =
            3.0 * g_0_yyyyyzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyy_0[i] * pb_x + g_0_yyyyyzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyyz_0[i] =
            3.0 * g_0_yyyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyz_0[i] * pb_x + g_0_yyyyyzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyzz_0[i] =
            3.0 * g_0_yyyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxzzz_0[i] =
            3.0 * g_0_yyyyyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyyy_0[i] =
            2.0 * g_0_yyyyyzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyy_0[i] * pb_x + g_0_yyyyyzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyyz_0[i] =
            2.0 * g_0_yyyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyz_0[i] * pb_x + g_0_yyyyyzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyzz_0[i] =
            2.0 * g_0_yyyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyzz_0[i] * pb_x + g_0_yyyyyzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyzzz_0[i] =
            2.0 * g_0_yyyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxzzzz_0[i] =
            2.0 * g_0_yyyyyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyyy_0[i] = g_0_yyyyyzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyy_0[i] * pb_x + g_0_yyyyyzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyyz_0[i] = g_0_yyyyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyz_0[i] * pb_x + g_0_yyyyyzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyzz_0[i] = g_0_yyyyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyzz_0[i] * pb_x + g_0_yyyyyzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyzzz_0[i] = g_0_yyyyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyzzz_0[i] * pb_x + g_0_yyyyyzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyzzzz_0[i] = g_0_yyyyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xzzzzz_0[i] = g_0_yyyyyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyyy_0[i] = g_0_yyyyyzz_0_yyyyyy_0[i] * pb_x + g_0_yyyyyzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyyz_0[i] = g_0_yyyyyzz_0_yyyyyz_0[i] * pb_x + g_0_yyyyyzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyzz_0[i] = g_0_yyyyyzz_0_yyyyzz_0[i] * pb_x + g_0_yyyyyzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyzzz_0[i] = g_0_yyyyyzz_0_yyyzzz_0[i] * pb_x + g_0_yyyyyzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyzzzz_0[i] = g_0_yyyyyzz_0_yyzzzz_0[i] * pb_x + g_0_yyyyyzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yzzzzz_0[i] = g_0_yyyyyzz_0_yzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_zzzzzz_0[i] = g_0_yyyyyzz_0_zzzzzz_0[i] * pb_x + g_0_yyyyyzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 868-896 components of targeted buffer : SLSI

    auto g_0_xyyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 868);

    auto g_0_xyyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 869);

    auto g_0_xyyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 870);

    auto g_0_xyyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 871);

    auto g_0_xyyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 872);

    auto g_0_xyyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 873);

    auto g_0_xyyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 874);

    auto g_0_xyyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 875);

    auto g_0_xyyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 876);

    auto g_0_xyyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 877);

    auto g_0_xyyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 878);

    auto g_0_xyyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 879);

    auto g_0_xyyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 880);

    auto g_0_xyyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 881);

    auto g_0_xyyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 882);

    auto g_0_xyyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 883);

    auto g_0_xyyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 884);

    auto g_0_xyyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 885);

    auto g_0_xyyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 886);

    auto g_0_xyyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 887);

    auto g_0_xyyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 888);

    auto g_0_xyyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 889);

    auto g_0_xyyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 890);

    auto g_0_xyyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 891);

    auto g_0_xyyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 892);

    auto g_0_xyyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 893);

    auto g_0_xyyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 894);

    auto g_0_xyyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 895);

#pragma omp simd aligned(g_0_xyyyyzzz_0_xxxxxx_0,     \
                             g_0_xyyyyzzz_0_xxxxxy_0, \
                             g_0_xyyyyzzz_0_xxxxxz_0, \
                             g_0_xyyyyzzz_0_xxxxyy_0, \
                             g_0_xyyyyzzz_0_xxxxyz_0, \
                             g_0_xyyyyzzz_0_xxxxzz_0, \
                             g_0_xyyyyzzz_0_xxxyyy_0, \
                             g_0_xyyyyzzz_0_xxxyyz_0, \
                             g_0_xyyyyzzz_0_xxxyzz_0, \
                             g_0_xyyyyzzz_0_xxxzzz_0, \
                             g_0_xyyyyzzz_0_xxyyyy_0, \
                             g_0_xyyyyzzz_0_xxyyyz_0, \
                             g_0_xyyyyzzz_0_xxyyzz_0, \
                             g_0_xyyyyzzz_0_xxyzzz_0, \
                             g_0_xyyyyzzz_0_xxzzzz_0, \
                             g_0_xyyyyzzz_0_xyyyyy_0, \
                             g_0_xyyyyzzz_0_xyyyyz_0, \
                             g_0_xyyyyzzz_0_xyyyzz_0, \
                             g_0_xyyyyzzz_0_xyyzzz_0, \
                             g_0_xyyyyzzz_0_xyzzzz_0, \
                             g_0_xyyyyzzz_0_xzzzzz_0, \
                             g_0_xyyyyzzz_0_yyyyyy_0, \
                             g_0_xyyyyzzz_0_yyyyyz_0, \
                             g_0_xyyyyzzz_0_yyyyzz_0, \
                             g_0_xyyyyzzz_0_yyyzzz_0, \
                             g_0_xyyyyzzz_0_yyzzzz_0, \
                             g_0_xyyyyzzz_0_yzzzzz_0, \
                             g_0_xyyyyzzz_0_zzzzzz_0, \
                             g_0_yyyyzzz_0_xxxxx_1,   \
                             g_0_yyyyzzz_0_xxxxxx_0,  \
                             g_0_yyyyzzz_0_xxxxxx_1,  \
                             g_0_yyyyzzz_0_xxxxxy_0,  \
                             g_0_yyyyzzz_0_xxxxxy_1,  \
                             g_0_yyyyzzz_0_xxxxxz_0,  \
                             g_0_yyyyzzz_0_xxxxxz_1,  \
                             g_0_yyyyzzz_0_xxxxy_1,   \
                             g_0_yyyyzzz_0_xxxxyy_0,  \
                             g_0_yyyyzzz_0_xxxxyy_1,  \
                             g_0_yyyyzzz_0_xxxxyz_0,  \
                             g_0_yyyyzzz_0_xxxxyz_1,  \
                             g_0_yyyyzzz_0_xxxxz_1,   \
                             g_0_yyyyzzz_0_xxxxzz_0,  \
                             g_0_yyyyzzz_0_xxxxzz_1,  \
                             g_0_yyyyzzz_0_xxxyy_1,   \
                             g_0_yyyyzzz_0_xxxyyy_0,  \
                             g_0_yyyyzzz_0_xxxyyy_1,  \
                             g_0_yyyyzzz_0_xxxyyz_0,  \
                             g_0_yyyyzzz_0_xxxyyz_1,  \
                             g_0_yyyyzzz_0_xxxyz_1,   \
                             g_0_yyyyzzz_0_xxxyzz_0,  \
                             g_0_yyyyzzz_0_xxxyzz_1,  \
                             g_0_yyyyzzz_0_xxxzz_1,   \
                             g_0_yyyyzzz_0_xxxzzz_0,  \
                             g_0_yyyyzzz_0_xxxzzz_1,  \
                             g_0_yyyyzzz_0_xxyyy_1,   \
                             g_0_yyyyzzz_0_xxyyyy_0,  \
                             g_0_yyyyzzz_0_xxyyyy_1,  \
                             g_0_yyyyzzz_0_xxyyyz_0,  \
                             g_0_yyyyzzz_0_xxyyyz_1,  \
                             g_0_yyyyzzz_0_xxyyz_1,   \
                             g_0_yyyyzzz_0_xxyyzz_0,  \
                             g_0_yyyyzzz_0_xxyyzz_1,  \
                             g_0_yyyyzzz_0_xxyzz_1,   \
                             g_0_yyyyzzz_0_xxyzzz_0,  \
                             g_0_yyyyzzz_0_xxyzzz_1,  \
                             g_0_yyyyzzz_0_xxzzz_1,   \
                             g_0_yyyyzzz_0_xxzzzz_0,  \
                             g_0_yyyyzzz_0_xxzzzz_1,  \
                             g_0_yyyyzzz_0_xyyyy_1,   \
                             g_0_yyyyzzz_0_xyyyyy_0,  \
                             g_0_yyyyzzz_0_xyyyyy_1,  \
                             g_0_yyyyzzz_0_xyyyyz_0,  \
                             g_0_yyyyzzz_0_xyyyyz_1,  \
                             g_0_yyyyzzz_0_xyyyz_1,   \
                             g_0_yyyyzzz_0_xyyyzz_0,  \
                             g_0_yyyyzzz_0_xyyyzz_1,  \
                             g_0_yyyyzzz_0_xyyzz_1,   \
                             g_0_yyyyzzz_0_xyyzzz_0,  \
                             g_0_yyyyzzz_0_xyyzzz_1,  \
                             g_0_yyyyzzz_0_xyzzz_1,   \
                             g_0_yyyyzzz_0_xyzzzz_0,  \
                             g_0_yyyyzzz_0_xyzzzz_1,  \
                             g_0_yyyyzzz_0_xzzzz_1,   \
                             g_0_yyyyzzz_0_xzzzzz_0,  \
                             g_0_yyyyzzz_0_xzzzzz_1,  \
                             g_0_yyyyzzz_0_yyyyy_1,   \
                             g_0_yyyyzzz_0_yyyyyy_0,  \
                             g_0_yyyyzzz_0_yyyyyy_1,  \
                             g_0_yyyyzzz_0_yyyyyz_0,  \
                             g_0_yyyyzzz_0_yyyyyz_1,  \
                             g_0_yyyyzzz_0_yyyyz_1,   \
                             g_0_yyyyzzz_0_yyyyzz_0,  \
                             g_0_yyyyzzz_0_yyyyzz_1,  \
                             g_0_yyyyzzz_0_yyyzz_1,   \
                             g_0_yyyyzzz_0_yyyzzz_0,  \
                             g_0_yyyyzzz_0_yyyzzz_1,  \
                             g_0_yyyyzzz_0_yyzzz_1,   \
                             g_0_yyyyzzz_0_yyzzzz_0,  \
                             g_0_yyyyzzz_0_yyzzzz_1,  \
                             g_0_yyyyzzz_0_yzzzz_1,   \
                             g_0_yyyyzzz_0_yzzzzz_0,  \
                             g_0_yyyyzzz_0_yzzzzz_1,  \
                             g_0_yyyyzzz_0_zzzzz_1,   \
                             g_0_yyyyzzz_0_zzzzzz_0,  \
                             g_0_yyyyzzz_0_zzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_xxxxxx_0[i] =
            6.0 * g_0_yyyyzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxx_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxy_0[i] =
            5.0 * g_0_yyyyzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxy_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxxz_0[i] =
            5.0 * g_0_yyyyzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxxz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxyy_0[i] =
            4.0 * g_0_yyyyzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyy_0[i] * pb_x + g_0_yyyyzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxyz_0[i] =
            4.0 * g_0_yyyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxzz_0[i] =
            4.0 * g_0_yyyyzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyyy_0[i] =
            3.0 * g_0_yyyyzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyy_0[i] * pb_x + g_0_yyyyzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyyz_0[i] =
            3.0 * g_0_yyyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyz_0[i] * pb_x + g_0_yyyyzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyzz_0[i] =
            3.0 * g_0_yyyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxzzz_0[i] =
            3.0 * g_0_yyyyzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyyy_0[i] =
            2.0 * g_0_yyyyzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyy_0[i] * pb_x + g_0_yyyyzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyyz_0[i] =
            2.0 * g_0_yyyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyz_0[i] * pb_x + g_0_yyyyzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyzz_0[i] =
            2.0 * g_0_yyyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyzz_0[i] * pb_x + g_0_yyyyzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyzzz_0[i] =
            2.0 * g_0_yyyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxzzzz_0[i] =
            2.0 * g_0_yyyyzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyyy_0[i] = g_0_yyyyzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyy_0[i] * pb_x + g_0_yyyyzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyyz_0[i] = g_0_yyyyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyz_0[i] * pb_x + g_0_yyyyzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyzz_0[i] = g_0_yyyyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyzz_0[i] * pb_x + g_0_yyyyzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyzzz_0[i] = g_0_yyyyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyzzz_0[i] * pb_x + g_0_yyyyzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyzzzz_0[i] = g_0_yyyyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xzzzzz_0[i] = g_0_yyyyzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyyy_0[i] = g_0_yyyyzzz_0_yyyyyy_0[i] * pb_x + g_0_yyyyzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyyz_0[i] = g_0_yyyyzzz_0_yyyyyz_0[i] * pb_x + g_0_yyyyzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyzz_0[i] = g_0_yyyyzzz_0_yyyyzz_0[i] * pb_x + g_0_yyyyzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyzzz_0[i] = g_0_yyyyzzz_0_yyyzzz_0[i] * pb_x + g_0_yyyyzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyzzzz_0[i] = g_0_yyyyzzz_0_yyzzzz_0[i] * pb_x + g_0_yyyyzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yzzzzz_0[i] = g_0_yyyyzzz_0_yzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_zzzzzz_0[i] = g_0_yyyyzzz_0_zzzzzz_0[i] * pb_x + g_0_yyyyzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 896-924 components of targeted buffer : SLSI

    auto g_0_xyyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 896);

    auto g_0_xyyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 897);

    auto g_0_xyyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 898);

    auto g_0_xyyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 899);

    auto g_0_xyyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 900);

    auto g_0_xyyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 901);

    auto g_0_xyyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 902);

    auto g_0_xyyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 903);

    auto g_0_xyyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 904);

    auto g_0_xyyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 905);

    auto g_0_xyyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 906);

    auto g_0_xyyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 907);

    auto g_0_xyyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 908);

    auto g_0_xyyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 909);

    auto g_0_xyyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 910);

    auto g_0_xyyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 911);

    auto g_0_xyyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 912);

    auto g_0_xyyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 913);

    auto g_0_xyyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 914);

    auto g_0_xyyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 915);

    auto g_0_xyyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 916);

    auto g_0_xyyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 917);

    auto g_0_xyyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 918);

    auto g_0_xyyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 919);

    auto g_0_xyyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 920);

    auto g_0_xyyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 921);

    auto g_0_xyyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 922);

    auto g_0_xyyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 923);

#pragma omp simd aligned(g_0_xyyyzzzz_0_xxxxxx_0,     \
                             g_0_xyyyzzzz_0_xxxxxy_0, \
                             g_0_xyyyzzzz_0_xxxxxz_0, \
                             g_0_xyyyzzzz_0_xxxxyy_0, \
                             g_0_xyyyzzzz_0_xxxxyz_0, \
                             g_0_xyyyzzzz_0_xxxxzz_0, \
                             g_0_xyyyzzzz_0_xxxyyy_0, \
                             g_0_xyyyzzzz_0_xxxyyz_0, \
                             g_0_xyyyzzzz_0_xxxyzz_0, \
                             g_0_xyyyzzzz_0_xxxzzz_0, \
                             g_0_xyyyzzzz_0_xxyyyy_0, \
                             g_0_xyyyzzzz_0_xxyyyz_0, \
                             g_0_xyyyzzzz_0_xxyyzz_0, \
                             g_0_xyyyzzzz_0_xxyzzz_0, \
                             g_0_xyyyzzzz_0_xxzzzz_0, \
                             g_0_xyyyzzzz_0_xyyyyy_0, \
                             g_0_xyyyzzzz_0_xyyyyz_0, \
                             g_0_xyyyzzzz_0_xyyyzz_0, \
                             g_0_xyyyzzzz_0_xyyzzz_0, \
                             g_0_xyyyzzzz_0_xyzzzz_0, \
                             g_0_xyyyzzzz_0_xzzzzz_0, \
                             g_0_xyyyzzzz_0_yyyyyy_0, \
                             g_0_xyyyzzzz_0_yyyyyz_0, \
                             g_0_xyyyzzzz_0_yyyyzz_0, \
                             g_0_xyyyzzzz_0_yyyzzz_0, \
                             g_0_xyyyzzzz_0_yyzzzz_0, \
                             g_0_xyyyzzzz_0_yzzzzz_0, \
                             g_0_xyyyzzzz_0_zzzzzz_0, \
                             g_0_yyyzzzz_0_xxxxx_1,   \
                             g_0_yyyzzzz_0_xxxxxx_0,  \
                             g_0_yyyzzzz_0_xxxxxx_1,  \
                             g_0_yyyzzzz_0_xxxxxy_0,  \
                             g_0_yyyzzzz_0_xxxxxy_1,  \
                             g_0_yyyzzzz_0_xxxxxz_0,  \
                             g_0_yyyzzzz_0_xxxxxz_1,  \
                             g_0_yyyzzzz_0_xxxxy_1,   \
                             g_0_yyyzzzz_0_xxxxyy_0,  \
                             g_0_yyyzzzz_0_xxxxyy_1,  \
                             g_0_yyyzzzz_0_xxxxyz_0,  \
                             g_0_yyyzzzz_0_xxxxyz_1,  \
                             g_0_yyyzzzz_0_xxxxz_1,   \
                             g_0_yyyzzzz_0_xxxxzz_0,  \
                             g_0_yyyzzzz_0_xxxxzz_1,  \
                             g_0_yyyzzzz_0_xxxyy_1,   \
                             g_0_yyyzzzz_0_xxxyyy_0,  \
                             g_0_yyyzzzz_0_xxxyyy_1,  \
                             g_0_yyyzzzz_0_xxxyyz_0,  \
                             g_0_yyyzzzz_0_xxxyyz_1,  \
                             g_0_yyyzzzz_0_xxxyz_1,   \
                             g_0_yyyzzzz_0_xxxyzz_0,  \
                             g_0_yyyzzzz_0_xxxyzz_1,  \
                             g_0_yyyzzzz_0_xxxzz_1,   \
                             g_0_yyyzzzz_0_xxxzzz_0,  \
                             g_0_yyyzzzz_0_xxxzzz_1,  \
                             g_0_yyyzzzz_0_xxyyy_1,   \
                             g_0_yyyzzzz_0_xxyyyy_0,  \
                             g_0_yyyzzzz_0_xxyyyy_1,  \
                             g_0_yyyzzzz_0_xxyyyz_0,  \
                             g_0_yyyzzzz_0_xxyyyz_1,  \
                             g_0_yyyzzzz_0_xxyyz_1,   \
                             g_0_yyyzzzz_0_xxyyzz_0,  \
                             g_0_yyyzzzz_0_xxyyzz_1,  \
                             g_0_yyyzzzz_0_xxyzz_1,   \
                             g_0_yyyzzzz_0_xxyzzz_0,  \
                             g_0_yyyzzzz_0_xxyzzz_1,  \
                             g_0_yyyzzzz_0_xxzzz_1,   \
                             g_0_yyyzzzz_0_xxzzzz_0,  \
                             g_0_yyyzzzz_0_xxzzzz_1,  \
                             g_0_yyyzzzz_0_xyyyy_1,   \
                             g_0_yyyzzzz_0_xyyyyy_0,  \
                             g_0_yyyzzzz_0_xyyyyy_1,  \
                             g_0_yyyzzzz_0_xyyyyz_0,  \
                             g_0_yyyzzzz_0_xyyyyz_1,  \
                             g_0_yyyzzzz_0_xyyyz_1,   \
                             g_0_yyyzzzz_0_xyyyzz_0,  \
                             g_0_yyyzzzz_0_xyyyzz_1,  \
                             g_0_yyyzzzz_0_xyyzz_1,   \
                             g_0_yyyzzzz_0_xyyzzz_0,  \
                             g_0_yyyzzzz_0_xyyzzz_1,  \
                             g_0_yyyzzzz_0_xyzzz_1,   \
                             g_0_yyyzzzz_0_xyzzzz_0,  \
                             g_0_yyyzzzz_0_xyzzzz_1,  \
                             g_0_yyyzzzz_0_xzzzz_1,   \
                             g_0_yyyzzzz_0_xzzzzz_0,  \
                             g_0_yyyzzzz_0_xzzzzz_1,  \
                             g_0_yyyzzzz_0_yyyyy_1,   \
                             g_0_yyyzzzz_0_yyyyyy_0,  \
                             g_0_yyyzzzz_0_yyyyyy_1,  \
                             g_0_yyyzzzz_0_yyyyyz_0,  \
                             g_0_yyyzzzz_0_yyyyyz_1,  \
                             g_0_yyyzzzz_0_yyyyz_1,   \
                             g_0_yyyzzzz_0_yyyyzz_0,  \
                             g_0_yyyzzzz_0_yyyyzz_1,  \
                             g_0_yyyzzzz_0_yyyzz_1,   \
                             g_0_yyyzzzz_0_yyyzzz_0,  \
                             g_0_yyyzzzz_0_yyyzzz_1,  \
                             g_0_yyyzzzz_0_yyzzz_1,   \
                             g_0_yyyzzzz_0_yyzzzz_0,  \
                             g_0_yyyzzzz_0_yyzzzz_1,  \
                             g_0_yyyzzzz_0_yzzzz_1,   \
                             g_0_yyyzzzz_0_yzzzzz_0,  \
                             g_0_yyyzzzz_0_yzzzzz_1,  \
                             g_0_yyyzzzz_0_zzzzz_1,   \
                             g_0_yyyzzzz_0_zzzzzz_0,  \
                             g_0_yyyzzzz_0_zzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_xxxxxx_0[i] =
            6.0 * g_0_yyyzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxx_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxy_0[i] =
            5.0 * g_0_yyyzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxy_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxxz_0[i] =
            5.0 * g_0_yyyzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxxz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxyy_0[i] =
            4.0 * g_0_yyyzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyy_0[i] * pb_x + g_0_yyyzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxyz_0[i] =
            4.0 * g_0_yyyzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxzz_0[i] =
            4.0 * g_0_yyyzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_yyyzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyy_0[i] * pb_x + g_0_yyyzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyyz_0[i] =
            3.0 * g_0_yyyzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyz_0[i] * pb_x + g_0_yyyzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyzz_0[i] =
            3.0 * g_0_yyyzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxzzz_0[i] =
            3.0 * g_0_yyyzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyyy_0[i] =
            2.0 * g_0_yyyzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyy_0[i] * pb_x + g_0_yyyzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyyz_0[i] =
            2.0 * g_0_yyyzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyz_0[i] * pb_x + g_0_yyyzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_yyyzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyzz_0[i] * pb_x + g_0_yyyzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyzzz_0[i] =
            2.0 * g_0_yyyzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxzzzz_0[i] =
            2.0 * g_0_yyyzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyyy_0[i] = g_0_yyyzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyy_0[i] * pb_x + g_0_yyyzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyyz_0[i] = g_0_yyyzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyz_0[i] * pb_x + g_0_yyyzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyzz_0[i] = g_0_yyyzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyzz_0[i] * pb_x + g_0_yyyzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyzzz_0[i] = g_0_yyyzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyzzz_0[i] * pb_x + g_0_yyyzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyzzzz_0[i] = g_0_yyyzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xzzzzz_0[i] = g_0_yyyzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyyy_0[i] = g_0_yyyzzzz_0_yyyyyy_0[i] * pb_x + g_0_yyyzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyyz_0[i] = g_0_yyyzzzz_0_yyyyyz_0[i] * pb_x + g_0_yyyzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyzz_0[i] = g_0_yyyzzzz_0_yyyyzz_0[i] * pb_x + g_0_yyyzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyzzz_0[i] = g_0_yyyzzzz_0_yyyzzz_0[i] * pb_x + g_0_yyyzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyzzzz_0[i] = g_0_yyyzzzz_0_yyzzzz_0[i] * pb_x + g_0_yyyzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yzzzzz_0[i] = g_0_yyyzzzz_0_yzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_zzzzzz_0[i] = g_0_yyyzzzz_0_zzzzzz_0[i] * pb_x + g_0_yyyzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 924-952 components of targeted buffer : SLSI

    auto g_0_xyyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 924);

    auto g_0_xyyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 925);

    auto g_0_xyyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 926);

    auto g_0_xyyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 927);

    auto g_0_xyyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 928);

    auto g_0_xyyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 929);

    auto g_0_xyyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 930);

    auto g_0_xyyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 931);

    auto g_0_xyyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 932);

    auto g_0_xyyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 933);

    auto g_0_xyyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 934);

    auto g_0_xyyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 935);

    auto g_0_xyyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 936);

    auto g_0_xyyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 937);

    auto g_0_xyyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 938);

    auto g_0_xyyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 939);

    auto g_0_xyyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 940);

    auto g_0_xyyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 941);

    auto g_0_xyyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 942);

    auto g_0_xyyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 943);

    auto g_0_xyyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 944);

    auto g_0_xyyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 945);

    auto g_0_xyyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 946);

    auto g_0_xyyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 947);

    auto g_0_xyyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 948);

    auto g_0_xyyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 949);

    auto g_0_xyyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 950);

    auto g_0_xyyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 951);

#pragma omp simd aligned(g_0_xyyzzzzz_0_xxxxxx_0,     \
                             g_0_xyyzzzzz_0_xxxxxy_0, \
                             g_0_xyyzzzzz_0_xxxxxz_0, \
                             g_0_xyyzzzzz_0_xxxxyy_0, \
                             g_0_xyyzzzzz_0_xxxxyz_0, \
                             g_0_xyyzzzzz_0_xxxxzz_0, \
                             g_0_xyyzzzzz_0_xxxyyy_0, \
                             g_0_xyyzzzzz_0_xxxyyz_0, \
                             g_0_xyyzzzzz_0_xxxyzz_0, \
                             g_0_xyyzzzzz_0_xxxzzz_0, \
                             g_0_xyyzzzzz_0_xxyyyy_0, \
                             g_0_xyyzzzzz_0_xxyyyz_0, \
                             g_0_xyyzzzzz_0_xxyyzz_0, \
                             g_0_xyyzzzzz_0_xxyzzz_0, \
                             g_0_xyyzzzzz_0_xxzzzz_0, \
                             g_0_xyyzzzzz_0_xyyyyy_0, \
                             g_0_xyyzzzzz_0_xyyyyz_0, \
                             g_0_xyyzzzzz_0_xyyyzz_0, \
                             g_0_xyyzzzzz_0_xyyzzz_0, \
                             g_0_xyyzzzzz_0_xyzzzz_0, \
                             g_0_xyyzzzzz_0_xzzzzz_0, \
                             g_0_xyyzzzzz_0_yyyyyy_0, \
                             g_0_xyyzzzzz_0_yyyyyz_0, \
                             g_0_xyyzzzzz_0_yyyyzz_0, \
                             g_0_xyyzzzzz_0_yyyzzz_0, \
                             g_0_xyyzzzzz_0_yyzzzz_0, \
                             g_0_xyyzzzzz_0_yzzzzz_0, \
                             g_0_xyyzzzzz_0_zzzzzz_0, \
                             g_0_yyzzzzz_0_xxxxx_1,   \
                             g_0_yyzzzzz_0_xxxxxx_0,  \
                             g_0_yyzzzzz_0_xxxxxx_1,  \
                             g_0_yyzzzzz_0_xxxxxy_0,  \
                             g_0_yyzzzzz_0_xxxxxy_1,  \
                             g_0_yyzzzzz_0_xxxxxz_0,  \
                             g_0_yyzzzzz_0_xxxxxz_1,  \
                             g_0_yyzzzzz_0_xxxxy_1,   \
                             g_0_yyzzzzz_0_xxxxyy_0,  \
                             g_0_yyzzzzz_0_xxxxyy_1,  \
                             g_0_yyzzzzz_0_xxxxyz_0,  \
                             g_0_yyzzzzz_0_xxxxyz_1,  \
                             g_0_yyzzzzz_0_xxxxz_1,   \
                             g_0_yyzzzzz_0_xxxxzz_0,  \
                             g_0_yyzzzzz_0_xxxxzz_1,  \
                             g_0_yyzzzzz_0_xxxyy_1,   \
                             g_0_yyzzzzz_0_xxxyyy_0,  \
                             g_0_yyzzzzz_0_xxxyyy_1,  \
                             g_0_yyzzzzz_0_xxxyyz_0,  \
                             g_0_yyzzzzz_0_xxxyyz_1,  \
                             g_0_yyzzzzz_0_xxxyz_1,   \
                             g_0_yyzzzzz_0_xxxyzz_0,  \
                             g_0_yyzzzzz_0_xxxyzz_1,  \
                             g_0_yyzzzzz_0_xxxzz_1,   \
                             g_0_yyzzzzz_0_xxxzzz_0,  \
                             g_0_yyzzzzz_0_xxxzzz_1,  \
                             g_0_yyzzzzz_0_xxyyy_1,   \
                             g_0_yyzzzzz_0_xxyyyy_0,  \
                             g_0_yyzzzzz_0_xxyyyy_1,  \
                             g_0_yyzzzzz_0_xxyyyz_0,  \
                             g_0_yyzzzzz_0_xxyyyz_1,  \
                             g_0_yyzzzzz_0_xxyyz_1,   \
                             g_0_yyzzzzz_0_xxyyzz_0,  \
                             g_0_yyzzzzz_0_xxyyzz_1,  \
                             g_0_yyzzzzz_0_xxyzz_1,   \
                             g_0_yyzzzzz_0_xxyzzz_0,  \
                             g_0_yyzzzzz_0_xxyzzz_1,  \
                             g_0_yyzzzzz_0_xxzzz_1,   \
                             g_0_yyzzzzz_0_xxzzzz_0,  \
                             g_0_yyzzzzz_0_xxzzzz_1,  \
                             g_0_yyzzzzz_0_xyyyy_1,   \
                             g_0_yyzzzzz_0_xyyyyy_0,  \
                             g_0_yyzzzzz_0_xyyyyy_1,  \
                             g_0_yyzzzzz_0_xyyyyz_0,  \
                             g_0_yyzzzzz_0_xyyyyz_1,  \
                             g_0_yyzzzzz_0_xyyyz_1,   \
                             g_0_yyzzzzz_0_xyyyzz_0,  \
                             g_0_yyzzzzz_0_xyyyzz_1,  \
                             g_0_yyzzzzz_0_xyyzz_1,   \
                             g_0_yyzzzzz_0_xyyzzz_0,  \
                             g_0_yyzzzzz_0_xyyzzz_1,  \
                             g_0_yyzzzzz_0_xyzzz_1,   \
                             g_0_yyzzzzz_0_xyzzzz_0,  \
                             g_0_yyzzzzz_0_xyzzzz_1,  \
                             g_0_yyzzzzz_0_xzzzz_1,   \
                             g_0_yyzzzzz_0_xzzzzz_0,  \
                             g_0_yyzzzzz_0_xzzzzz_1,  \
                             g_0_yyzzzzz_0_yyyyy_1,   \
                             g_0_yyzzzzz_0_yyyyyy_0,  \
                             g_0_yyzzzzz_0_yyyyyy_1,  \
                             g_0_yyzzzzz_0_yyyyyz_0,  \
                             g_0_yyzzzzz_0_yyyyyz_1,  \
                             g_0_yyzzzzz_0_yyyyz_1,   \
                             g_0_yyzzzzz_0_yyyyzz_0,  \
                             g_0_yyzzzzz_0_yyyyzz_1,  \
                             g_0_yyzzzzz_0_yyyzz_1,   \
                             g_0_yyzzzzz_0_yyyzzz_0,  \
                             g_0_yyzzzzz_0_yyyzzz_1,  \
                             g_0_yyzzzzz_0_yyzzz_1,   \
                             g_0_yyzzzzz_0_yyzzzz_0,  \
                             g_0_yyzzzzz_0_yyzzzz_1,  \
                             g_0_yyzzzzz_0_yzzzz_1,   \
                             g_0_yyzzzzz_0_yzzzzz_0,  \
                             g_0_yyzzzzz_0_yzzzzz_1,  \
                             g_0_yyzzzzz_0_zzzzz_1,   \
                             g_0_yyzzzzz_0_zzzzzz_0,  \
                             g_0_yyzzzzz_0_zzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_xxxxxx_0[i] =
            6.0 * g_0_yyzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxx_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxy_0[i] =
            5.0 * g_0_yyzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxy_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxxz_0[i] =
            5.0 * g_0_yyzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxxz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxyy_0[i] =
            4.0 * g_0_yyzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyy_0[i] * pb_x + g_0_yyzzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxyz_0[i] =
            4.0 * g_0_yyzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxzz_0[i] =
            4.0 * g_0_yyzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_yyzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyy_0[i] * pb_x + g_0_yyzzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyyz_0[i] =
            3.0 * g_0_yyzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyz_0[i] * pb_x + g_0_yyzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyzz_0[i] =
            3.0 * g_0_yyzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxzzz_0[i] =
            3.0 * g_0_yyzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyyy_0[i] =
            2.0 * g_0_yyzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyy_0[i] * pb_x + g_0_yyzzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyyz_0[i] =
            2.0 * g_0_yyzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyz_0[i] * pb_x + g_0_yyzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_yyzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyzz_0[i] * pb_x + g_0_yyzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyzzz_0[i] =
            2.0 * g_0_yyzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxzzzz_0[i] =
            2.0 * g_0_yyzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyyy_0[i] = g_0_yyzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyy_0[i] * pb_x + g_0_yyzzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyyz_0[i] = g_0_yyzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyz_0[i] * pb_x + g_0_yyzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyzz_0[i] = g_0_yyzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyzz_0[i] * pb_x + g_0_yyzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyzzz_0[i] = g_0_yyzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyzzz_0[i] * pb_x + g_0_yyzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyzzzz_0[i] = g_0_yyzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xzzzzz_0[i] = g_0_yyzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyyy_0[i] = g_0_yyzzzzz_0_yyyyyy_0[i] * pb_x + g_0_yyzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyyz_0[i] = g_0_yyzzzzz_0_yyyyyz_0[i] * pb_x + g_0_yyzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyzz_0[i] = g_0_yyzzzzz_0_yyyyzz_0[i] * pb_x + g_0_yyzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyzzz_0[i] = g_0_yyzzzzz_0_yyyzzz_0[i] * pb_x + g_0_yyzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyzzzz_0[i] = g_0_yyzzzzz_0_yyzzzz_0[i] * pb_x + g_0_yyzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yzzzzz_0[i] = g_0_yyzzzzz_0_yzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_zzzzzz_0[i] = g_0_yyzzzzz_0_zzzzzz_0[i] * pb_x + g_0_yyzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 952-980 components of targeted buffer : SLSI

    auto g_0_xyzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 952);

    auto g_0_xyzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 953);

    auto g_0_xyzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 954);

    auto g_0_xyzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 955);

    auto g_0_xyzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 956);

    auto g_0_xyzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 957);

    auto g_0_xyzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 958);

    auto g_0_xyzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 959);

    auto g_0_xyzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 960);

    auto g_0_xyzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 961);

    auto g_0_xyzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 962);

    auto g_0_xyzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 963);

    auto g_0_xyzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 964);

    auto g_0_xyzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 965);

    auto g_0_xyzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 966);

    auto g_0_xyzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 967);

    auto g_0_xyzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 968);

    auto g_0_xyzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 969);

    auto g_0_xyzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 970);

    auto g_0_xyzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 971);

    auto g_0_xyzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 972);

    auto g_0_xyzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 973);

    auto g_0_xyzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 974);

    auto g_0_xyzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 975);

    auto g_0_xyzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 976);

    auto g_0_xyzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 977);

    auto g_0_xyzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 978);

    auto g_0_xyzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 979);

#pragma omp simd aligned(g_0_xyzzzzzz_0_xxxxxx_0,     \
                             g_0_xyzzzzzz_0_xxxxxy_0, \
                             g_0_xyzzzzzz_0_xxxxxz_0, \
                             g_0_xyzzzzzz_0_xxxxyy_0, \
                             g_0_xyzzzzzz_0_xxxxyz_0, \
                             g_0_xyzzzzzz_0_xxxxzz_0, \
                             g_0_xyzzzzzz_0_xxxyyy_0, \
                             g_0_xyzzzzzz_0_xxxyyz_0, \
                             g_0_xyzzzzzz_0_xxxyzz_0, \
                             g_0_xyzzzzzz_0_xxxzzz_0, \
                             g_0_xyzzzzzz_0_xxyyyy_0, \
                             g_0_xyzzzzzz_0_xxyyyz_0, \
                             g_0_xyzzzzzz_0_xxyyzz_0, \
                             g_0_xyzzzzzz_0_xxyzzz_0, \
                             g_0_xyzzzzzz_0_xxzzzz_0, \
                             g_0_xyzzzzzz_0_xyyyyy_0, \
                             g_0_xyzzzzzz_0_xyyyyz_0, \
                             g_0_xyzzzzzz_0_xyyyzz_0, \
                             g_0_xyzzzzzz_0_xyyzzz_0, \
                             g_0_xyzzzzzz_0_xyzzzz_0, \
                             g_0_xyzzzzzz_0_xzzzzz_0, \
                             g_0_xyzzzzzz_0_yyyyyy_0, \
                             g_0_xyzzzzzz_0_yyyyyz_0, \
                             g_0_xyzzzzzz_0_yyyyzz_0, \
                             g_0_xyzzzzzz_0_yyyzzz_0, \
                             g_0_xyzzzzzz_0_yyzzzz_0, \
                             g_0_xyzzzzzz_0_yzzzzz_0, \
                             g_0_xyzzzzzz_0_zzzzzz_0, \
                             g_0_xzzzzzz_0_xxxxxx_0,  \
                             g_0_xzzzzzz_0_xxxxxx_1,  \
                             g_0_xzzzzzz_0_xxxxxz_0,  \
                             g_0_xzzzzzz_0_xxxxxz_1,  \
                             g_0_xzzzzzz_0_xxxxzz_0,  \
                             g_0_xzzzzzz_0_xxxxzz_1,  \
                             g_0_xzzzzzz_0_xxxzzz_0,  \
                             g_0_xzzzzzz_0_xxxzzz_1,  \
                             g_0_xzzzzzz_0_xxzzzz_0,  \
                             g_0_xzzzzzz_0_xxzzzz_1,  \
                             g_0_xzzzzzz_0_xzzzzz_0,  \
                             g_0_xzzzzzz_0_xzzzzz_1,  \
                             g_0_yzzzzzz_0_xxxxxy_0,  \
                             g_0_yzzzzzz_0_xxxxxy_1,  \
                             g_0_yzzzzzz_0_xxxxy_1,   \
                             g_0_yzzzzzz_0_xxxxyy_0,  \
                             g_0_yzzzzzz_0_xxxxyy_1,  \
                             g_0_yzzzzzz_0_xxxxyz_0,  \
                             g_0_yzzzzzz_0_xxxxyz_1,  \
                             g_0_yzzzzzz_0_xxxyy_1,   \
                             g_0_yzzzzzz_0_xxxyyy_0,  \
                             g_0_yzzzzzz_0_xxxyyy_1,  \
                             g_0_yzzzzzz_0_xxxyyz_0,  \
                             g_0_yzzzzzz_0_xxxyyz_1,  \
                             g_0_yzzzzzz_0_xxxyz_1,   \
                             g_0_yzzzzzz_0_xxxyzz_0,  \
                             g_0_yzzzzzz_0_xxxyzz_1,  \
                             g_0_yzzzzzz_0_xxyyy_1,   \
                             g_0_yzzzzzz_0_xxyyyy_0,  \
                             g_0_yzzzzzz_0_xxyyyy_1,  \
                             g_0_yzzzzzz_0_xxyyyz_0,  \
                             g_0_yzzzzzz_0_xxyyyz_1,  \
                             g_0_yzzzzzz_0_xxyyz_1,   \
                             g_0_yzzzzzz_0_xxyyzz_0,  \
                             g_0_yzzzzzz_0_xxyyzz_1,  \
                             g_0_yzzzzzz_0_xxyzz_1,   \
                             g_0_yzzzzzz_0_xxyzzz_0,  \
                             g_0_yzzzzzz_0_xxyzzz_1,  \
                             g_0_yzzzzzz_0_xyyyy_1,   \
                             g_0_yzzzzzz_0_xyyyyy_0,  \
                             g_0_yzzzzzz_0_xyyyyy_1,  \
                             g_0_yzzzzzz_0_xyyyyz_0,  \
                             g_0_yzzzzzz_0_xyyyyz_1,  \
                             g_0_yzzzzzz_0_xyyyz_1,   \
                             g_0_yzzzzzz_0_xyyyzz_0,  \
                             g_0_yzzzzzz_0_xyyyzz_1,  \
                             g_0_yzzzzzz_0_xyyzz_1,   \
                             g_0_yzzzzzz_0_xyyzzz_0,  \
                             g_0_yzzzzzz_0_xyyzzz_1,  \
                             g_0_yzzzzzz_0_xyzzz_1,   \
                             g_0_yzzzzzz_0_xyzzzz_0,  \
                             g_0_yzzzzzz_0_xyzzzz_1,  \
                             g_0_yzzzzzz_0_yyyyy_1,   \
                             g_0_yzzzzzz_0_yyyyyy_0,  \
                             g_0_yzzzzzz_0_yyyyyy_1,  \
                             g_0_yzzzzzz_0_yyyyyz_0,  \
                             g_0_yzzzzzz_0_yyyyyz_1,  \
                             g_0_yzzzzzz_0_yyyyz_1,   \
                             g_0_yzzzzzz_0_yyyyzz_0,  \
                             g_0_yzzzzzz_0_yyyyzz_1,  \
                             g_0_yzzzzzz_0_yyyzz_1,   \
                             g_0_yzzzzzz_0_yyyzzz_0,  \
                             g_0_yzzzzzz_0_yyyzzz_1,  \
                             g_0_yzzzzzz_0_yyzzz_1,   \
                             g_0_yzzzzzz_0_yyzzzz_0,  \
                             g_0_yzzzzzz_0_yyzzzz_1,  \
                             g_0_yzzzzzz_0_yzzzz_1,   \
                             g_0_yzzzzzz_0_yzzzzz_0,  \
                             g_0_yzzzzzz_0_yzzzzz_1,  \
                             g_0_yzzzzzz_0_zzzzzz_0,  \
                             g_0_yzzzzzz_0_zzzzzz_1,  \
                             wp_x,                    \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzzz_0_xxxxxx_0[i] = g_0_xzzzzzz_0_xxxxxx_0[i] * pb_y + g_0_xzzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxxxy_0[i] =
            5.0 * g_0_yzzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxxy_0[i] * pb_x + g_0_yzzzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxxz_0[i] = g_0_xzzzzzz_0_xxxxxz_0[i] * pb_y + g_0_xzzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxxyy_0[i] =
            4.0 * g_0_yzzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyy_0[i] * pb_x + g_0_yzzzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxyz_0[i] =
            4.0 * g_0_yzzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxyz_0[i] * pb_x + g_0_yzzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxzz_0[i] = g_0_xzzzzzz_0_xxxxzz_0[i] * pb_y + g_0_xzzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_yzzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyy_0[i] * pb_x + g_0_yzzzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxyyz_0[i] =
            3.0 * g_0_yzzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyz_0[i] * pb_x + g_0_yzzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxyzz_0[i] =
            3.0 * g_0_yzzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyzz_0[i] * pb_x + g_0_yzzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxzzz_0[i] = g_0_xzzzzzz_0_xxxzzz_0[i] * pb_y + g_0_xzzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxyyyy_0[i] =
            2.0 * g_0_yzzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyy_0[i] * pb_x + g_0_yzzzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyyyz_0[i] =
            2.0 * g_0_yzzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyz_0[i] * pb_x + g_0_yzzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_yzzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyzz_0[i] * pb_x + g_0_yzzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyzzz_0[i] =
            2.0 * g_0_yzzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyzzz_0[i] * pb_x + g_0_yzzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxzzzz_0[i] = g_0_xzzzzzz_0_xxzzzz_0[i] * pb_y + g_0_xzzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xyyyyy_0[i] = g_0_yzzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyy_0[i] * pb_x + g_0_yzzzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyyyz_0[i] = g_0_yzzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyz_0[i] * pb_x + g_0_yzzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyyzz_0[i] = g_0_yzzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyzz_0[i] * pb_x + g_0_yzzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyzzz_0[i] = g_0_yzzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyzzz_0[i] * pb_x + g_0_yzzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyzzzz_0[i] = g_0_yzzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzzzz_0[i] * pb_x + g_0_yzzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xzzzzz_0[i] = g_0_xzzzzzz_0_xzzzzz_0[i] * pb_y + g_0_xzzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_yyyyyy_0[i] = g_0_yzzzzzz_0_yyyyyy_0[i] * pb_x + g_0_yzzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyyyz_0[i] = g_0_yzzzzzz_0_yyyyyz_0[i] * pb_x + g_0_yzzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyyzz_0[i] = g_0_yzzzzzz_0_yyyyzz_0[i] * pb_x + g_0_yzzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyzzz_0[i] = g_0_yzzzzzz_0_yyyzzz_0[i] * pb_x + g_0_yzzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyzzzz_0[i] = g_0_yzzzzzz_0_yyzzzz_0[i] * pb_x + g_0_yzzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yzzzzz_0[i] = g_0_yzzzzzz_0_yzzzzz_0[i] * pb_x + g_0_yzzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_zzzzzz_0[i] = g_0_yzzzzzz_0_zzzzzz_0[i] * pb_x + g_0_yzzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 980-1008 components of targeted buffer : SLSI

    auto g_0_xzzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 980);

    auto g_0_xzzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 981);

    auto g_0_xzzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 982);

    auto g_0_xzzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 983);

    auto g_0_xzzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 984);

    auto g_0_xzzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 985);

    auto g_0_xzzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 986);

    auto g_0_xzzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 987);

    auto g_0_xzzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 988);

    auto g_0_xzzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 989);

    auto g_0_xzzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 990);

    auto g_0_xzzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 991);

    auto g_0_xzzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 992);

    auto g_0_xzzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 993);

    auto g_0_xzzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 994);

    auto g_0_xzzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 995);

    auto g_0_xzzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 996);

    auto g_0_xzzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 997);

    auto g_0_xzzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 998);

    auto g_0_xzzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 999);

    auto g_0_xzzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1000);

    auto g_0_xzzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1001);

    auto g_0_xzzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1002);

    auto g_0_xzzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1003);

    auto g_0_xzzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1004);

    auto g_0_xzzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1005);

    auto g_0_xzzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1006);

    auto g_0_xzzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1007);

#pragma omp simd aligned(g_0_xzzzzzzz_0_xxxxxx_0,     \
                             g_0_xzzzzzzz_0_xxxxxy_0, \
                             g_0_xzzzzzzz_0_xxxxxz_0, \
                             g_0_xzzzzzzz_0_xxxxyy_0, \
                             g_0_xzzzzzzz_0_xxxxyz_0, \
                             g_0_xzzzzzzz_0_xxxxzz_0, \
                             g_0_xzzzzzzz_0_xxxyyy_0, \
                             g_0_xzzzzzzz_0_xxxyyz_0, \
                             g_0_xzzzzzzz_0_xxxyzz_0, \
                             g_0_xzzzzzzz_0_xxxzzz_0, \
                             g_0_xzzzzzzz_0_xxyyyy_0, \
                             g_0_xzzzzzzz_0_xxyyyz_0, \
                             g_0_xzzzzzzz_0_xxyyzz_0, \
                             g_0_xzzzzzzz_0_xxyzzz_0, \
                             g_0_xzzzzzzz_0_xxzzzz_0, \
                             g_0_xzzzzzzz_0_xyyyyy_0, \
                             g_0_xzzzzzzz_0_xyyyyz_0, \
                             g_0_xzzzzzzz_0_xyyyzz_0, \
                             g_0_xzzzzzzz_0_xyyzzz_0, \
                             g_0_xzzzzzzz_0_xyzzzz_0, \
                             g_0_xzzzzzzz_0_xzzzzz_0, \
                             g_0_xzzzzzzz_0_yyyyyy_0, \
                             g_0_xzzzzzzz_0_yyyyyz_0, \
                             g_0_xzzzzzzz_0_yyyyzz_0, \
                             g_0_xzzzzzzz_0_yyyzzz_0, \
                             g_0_xzzzzzzz_0_yyzzzz_0, \
                             g_0_xzzzzzzz_0_yzzzzz_0, \
                             g_0_xzzzzzzz_0_zzzzzz_0, \
                             g_0_zzzzzzz_0_xxxxx_1,   \
                             g_0_zzzzzzz_0_xxxxxx_0,  \
                             g_0_zzzzzzz_0_xxxxxx_1,  \
                             g_0_zzzzzzz_0_xxxxxy_0,  \
                             g_0_zzzzzzz_0_xxxxxy_1,  \
                             g_0_zzzzzzz_0_xxxxxz_0,  \
                             g_0_zzzzzzz_0_xxxxxz_1,  \
                             g_0_zzzzzzz_0_xxxxy_1,   \
                             g_0_zzzzzzz_0_xxxxyy_0,  \
                             g_0_zzzzzzz_0_xxxxyy_1,  \
                             g_0_zzzzzzz_0_xxxxyz_0,  \
                             g_0_zzzzzzz_0_xxxxyz_1,  \
                             g_0_zzzzzzz_0_xxxxz_1,   \
                             g_0_zzzzzzz_0_xxxxzz_0,  \
                             g_0_zzzzzzz_0_xxxxzz_1,  \
                             g_0_zzzzzzz_0_xxxyy_1,   \
                             g_0_zzzzzzz_0_xxxyyy_0,  \
                             g_0_zzzzzzz_0_xxxyyy_1,  \
                             g_0_zzzzzzz_0_xxxyyz_0,  \
                             g_0_zzzzzzz_0_xxxyyz_1,  \
                             g_0_zzzzzzz_0_xxxyz_1,   \
                             g_0_zzzzzzz_0_xxxyzz_0,  \
                             g_0_zzzzzzz_0_xxxyzz_1,  \
                             g_0_zzzzzzz_0_xxxzz_1,   \
                             g_0_zzzzzzz_0_xxxzzz_0,  \
                             g_0_zzzzzzz_0_xxxzzz_1,  \
                             g_0_zzzzzzz_0_xxyyy_1,   \
                             g_0_zzzzzzz_0_xxyyyy_0,  \
                             g_0_zzzzzzz_0_xxyyyy_1,  \
                             g_0_zzzzzzz_0_xxyyyz_0,  \
                             g_0_zzzzzzz_0_xxyyyz_1,  \
                             g_0_zzzzzzz_0_xxyyz_1,   \
                             g_0_zzzzzzz_0_xxyyzz_0,  \
                             g_0_zzzzzzz_0_xxyyzz_1,  \
                             g_0_zzzzzzz_0_xxyzz_1,   \
                             g_0_zzzzzzz_0_xxyzzz_0,  \
                             g_0_zzzzzzz_0_xxyzzz_1,  \
                             g_0_zzzzzzz_0_xxzzz_1,   \
                             g_0_zzzzzzz_0_xxzzzz_0,  \
                             g_0_zzzzzzz_0_xxzzzz_1,  \
                             g_0_zzzzzzz_0_xyyyy_1,   \
                             g_0_zzzzzzz_0_xyyyyy_0,  \
                             g_0_zzzzzzz_0_xyyyyy_1,  \
                             g_0_zzzzzzz_0_xyyyyz_0,  \
                             g_0_zzzzzzz_0_xyyyyz_1,  \
                             g_0_zzzzzzz_0_xyyyz_1,   \
                             g_0_zzzzzzz_0_xyyyzz_0,  \
                             g_0_zzzzzzz_0_xyyyzz_1,  \
                             g_0_zzzzzzz_0_xyyzz_1,   \
                             g_0_zzzzzzz_0_xyyzzz_0,  \
                             g_0_zzzzzzz_0_xyyzzz_1,  \
                             g_0_zzzzzzz_0_xyzzz_1,   \
                             g_0_zzzzzzz_0_xyzzzz_0,  \
                             g_0_zzzzzzz_0_xyzzzz_1,  \
                             g_0_zzzzzzz_0_xzzzz_1,   \
                             g_0_zzzzzzz_0_xzzzzz_0,  \
                             g_0_zzzzzzz_0_xzzzzz_1,  \
                             g_0_zzzzzzz_0_yyyyy_1,   \
                             g_0_zzzzzzz_0_yyyyyy_0,  \
                             g_0_zzzzzzz_0_yyyyyy_1,  \
                             g_0_zzzzzzz_0_yyyyyz_0,  \
                             g_0_zzzzzzz_0_yyyyyz_1,  \
                             g_0_zzzzzzz_0_yyyyz_1,   \
                             g_0_zzzzzzz_0_yyyyzz_0,  \
                             g_0_zzzzzzz_0_yyyyzz_1,  \
                             g_0_zzzzzzz_0_yyyzz_1,   \
                             g_0_zzzzzzz_0_yyyzzz_0,  \
                             g_0_zzzzzzz_0_yyyzzz_1,  \
                             g_0_zzzzzzz_0_yyzzz_1,   \
                             g_0_zzzzzzz_0_yyzzzz_0,  \
                             g_0_zzzzzzz_0_yyzzzz_1,  \
                             g_0_zzzzzzz_0_yzzzz_1,   \
                             g_0_zzzzzzz_0_yzzzzz_0,  \
                             g_0_zzzzzzz_0_yzzzzz_1,  \
                             g_0_zzzzzzz_0_zzzzz_1,   \
                             g_0_zzzzzzz_0_zzzzzz_0,  \
                             g_0_zzzzzzz_0_zzzzzz_1,  \
                             wp_x,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_xxxxxx_0[i] =
            6.0 * g_0_zzzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxx_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxy_0[i] =
            5.0 * g_0_zzzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxy_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxxz_0[i] =
            5.0 * g_0_zzzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxyy_0[i] =
            4.0 * g_0_zzzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyy_0[i] * pb_x + g_0_zzzzzzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxyz_0[i] =
            4.0 * g_0_zzzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxzz_0[i] =
            4.0 * g_0_zzzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_zzzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyy_0[i] * pb_x + g_0_zzzzzzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyyz_0[i] =
            3.0 * g_0_zzzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyz_0[i] * pb_x + g_0_zzzzzzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyzz_0[i] =
            3.0 * g_0_zzzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxzzz_0[i] =
            3.0 * g_0_zzzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyyy_0[i] =
            2.0 * g_0_zzzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyy_0[i] * pb_x + g_0_zzzzzzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyyz_0[i] =
            2.0 * g_0_zzzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyz_0[i] * pb_x + g_0_zzzzzzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_zzzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyzz_0[i] * pb_x + g_0_zzzzzzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyzzz_0[i] =
            2.0 * g_0_zzzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxzzzz_0[i] =
            2.0 * g_0_zzzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyyy_0[i] = g_0_zzzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyy_0[i] * pb_x + g_0_zzzzzzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyyz_0[i] = g_0_zzzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyz_0[i] * pb_x + g_0_zzzzzzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyzz_0[i] = g_0_zzzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyzz_0[i] * pb_x + g_0_zzzzzzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyzzz_0[i] = g_0_zzzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzzz_0[i] * pb_x + g_0_zzzzzzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyzzzz_0[i] = g_0_zzzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xzzzzz_0[i] = g_0_zzzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyyy_0[i] = g_0_zzzzzzz_0_yyyyyy_0[i] * pb_x + g_0_zzzzzzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyyz_0[i] = g_0_zzzzzzz_0_yyyyyz_0[i] * pb_x + g_0_zzzzzzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyzz_0[i] = g_0_zzzzzzz_0_yyyyzz_0[i] * pb_x + g_0_zzzzzzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyzzz_0[i] = g_0_zzzzzzz_0_yyyzzz_0[i] * pb_x + g_0_zzzzzzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyzzzz_0[i] = g_0_zzzzzzz_0_yyzzzz_0[i] * pb_x + g_0_zzzzzzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yzzzzz_0[i] = g_0_zzzzzzz_0_yzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_zzzzzz_0[i] = g_0_zzzzzzz_0_zzzzzz_0[i] * pb_x + g_0_zzzzzzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 1008-1036 components of targeted buffer : SLSI

    auto g_0_yyyyyyyy_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1008);

    auto g_0_yyyyyyyy_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1009);

    auto g_0_yyyyyyyy_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1010);

    auto g_0_yyyyyyyy_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1011);

    auto g_0_yyyyyyyy_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1012);

    auto g_0_yyyyyyyy_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1013);

    auto g_0_yyyyyyyy_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1014);

    auto g_0_yyyyyyyy_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1015);

    auto g_0_yyyyyyyy_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1016);

    auto g_0_yyyyyyyy_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1017);

    auto g_0_yyyyyyyy_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1018);

    auto g_0_yyyyyyyy_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1019);

    auto g_0_yyyyyyyy_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1020);

    auto g_0_yyyyyyyy_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1021);

    auto g_0_yyyyyyyy_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1022);

    auto g_0_yyyyyyyy_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1023);

    auto g_0_yyyyyyyy_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1024);

    auto g_0_yyyyyyyy_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1025);

    auto g_0_yyyyyyyy_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1026);

    auto g_0_yyyyyyyy_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1027);

    auto g_0_yyyyyyyy_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1028);

    auto g_0_yyyyyyyy_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1029);

    auto g_0_yyyyyyyy_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1030);

    auto g_0_yyyyyyyy_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1031);

    auto g_0_yyyyyyyy_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1032);

    auto g_0_yyyyyyyy_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1033);

    auto g_0_yyyyyyyy_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1034);

    auto g_0_yyyyyyyy_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1035);

#pragma omp simd aligned(g_0_yyyyyy_0_xxxxxx_0,       \
                             g_0_yyyyyy_0_xxxxxx_1,   \
                             g_0_yyyyyy_0_xxxxxy_0,   \
                             g_0_yyyyyy_0_xxxxxy_1,   \
                             g_0_yyyyyy_0_xxxxxz_0,   \
                             g_0_yyyyyy_0_xxxxxz_1,   \
                             g_0_yyyyyy_0_xxxxyy_0,   \
                             g_0_yyyyyy_0_xxxxyy_1,   \
                             g_0_yyyyyy_0_xxxxyz_0,   \
                             g_0_yyyyyy_0_xxxxyz_1,   \
                             g_0_yyyyyy_0_xxxxzz_0,   \
                             g_0_yyyyyy_0_xxxxzz_1,   \
                             g_0_yyyyyy_0_xxxyyy_0,   \
                             g_0_yyyyyy_0_xxxyyy_1,   \
                             g_0_yyyyyy_0_xxxyyz_0,   \
                             g_0_yyyyyy_0_xxxyyz_1,   \
                             g_0_yyyyyy_0_xxxyzz_0,   \
                             g_0_yyyyyy_0_xxxyzz_1,   \
                             g_0_yyyyyy_0_xxxzzz_0,   \
                             g_0_yyyyyy_0_xxxzzz_1,   \
                             g_0_yyyyyy_0_xxyyyy_0,   \
                             g_0_yyyyyy_0_xxyyyy_1,   \
                             g_0_yyyyyy_0_xxyyyz_0,   \
                             g_0_yyyyyy_0_xxyyyz_1,   \
                             g_0_yyyyyy_0_xxyyzz_0,   \
                             g_0_yyyyyy_0_xxyyzz_1,   \
                             g_0_yyyyyy_0_xxyzzz_0,   \
                             g_0_yyyyyy_0_xxyzzz_1,   \
                             g_0_yyyyyy_0_xxzzzz_0,   \
                             g_0_yyyyyy_0_xxzzzz_1,   \
                             g_0_yyyyyy_0_xyyyyy_0,   \
                             g_0_yyyyyy_0_xyyyyy_1,   \
                             g_0_yyyyyy_0_xyyyyz_0,   \
                             g_0_yyyyyy_0_xyyyyz_1,   \
                             g_0_yyyyyy_0_xyyyzz_0,   \
                             g_0_yyyyyy_0_xyyyzz_1,   \
                             g_0_yyyyyy_0_xyyzzz_0,   \
                             g_0_yyyyyy_0_xyyzzz_1,   \
                             g_0_yyyyyy_0_xyzzzz_0,   \
                             g_0_yyyyyy_0_xyzzzz_1,   \
                             g_0_yyyyyy_0_xzzzzz_0,   \
                             g_0_yyyyyy_0_xzzzzz_1,   \
                             g_0_yyyyyy_0_yyyyyy_0,   \
                             g_0_yyyyyy_0_yyyyyy_1,   \
                             g_0_yyyyyy_0_yyyyyz_0,   \
                             g_0_yyyyyy_0_yyyyyz_1,   \
                             g_0_yyyyyy_0_yyyyzz_0,   \
                             g_0_yyyyyy_0_yyyyzz_1,   \
                             g_0_yyyyyy_0_yyyzzz_0,   \
                             g_0_yyyyyy_0_yyyzzz_1,   \
                             g_0_yyyyyy_0_yyzzzz_0,   \
                             g_0_yyyyyy_0_yyzzzz_1,   \
                             g_0_yyyyyy_0_yzzzzz_0,   \
                             g_0_yyyyyy_0_yzzzzz_1,   \
                             g_0_yyyyyy_0_zzzzzz_0,   \
                             g_0_yyyyyy_0_zzzzzz_1,   \
                             g_0_yyyyyyy_0_xxxxx_1,   \
                             g_0_yyyyyyy_0_xxxxxx_0,  \
                             g_0_yyyyyyy_0_xxxxxx_1,  \
                             g_0_yyyyyyy_0_xxxxxy_0,  \
                             g_0_yyyyyyy_0_xxxxxy_1,  \
                             g_0_yyyyyyy_0_xxxxxz_0,  \
                             g_0_yyyyyyy_0_xxxxxz_1,  \
                             g_0_yyyyyyy_0_xxxxy_1,   \
                             g_0_yyyyyyy_0_xxxxyy_0,  \
                             g_0_yyyyyyy_0_xxxxyy_1,  \
                             g_0_yyyyyyy_0_xxxxyz_0,  \
                             g_0_yyyyyyy_0_xxxxyz_1,  \
                             g_0_yyyyyyy_0_xxxxz_1,   \
                             g_0_yyyyyyy_0_xxxxzz_0,  \
                             g_0_yyyyyyy_0_xxxxzz_1,  \
                             g_0_yyyyyyy_0_xxxyy_1,   \
                             g_0_yyyyyyy_0_xxxyyy_0,  \
                             g_0_yyyyyyy_0_xxxyyy_1,  \
                             g_0_yyyyyyy_0_xxxyyz_0,  \
                             g_0_yyyyyyy_0_xxxyyz_1,  \
                             g_0_yyyyyyy_0_xxxyz_1,   \
                             g_0_yyyyyyy_0_xxxyzz_0,  \
                             g_0_yyyyyyy_0_xxxyzz_1,  \
                             g_0_yyyyyyy_0_xxxzz_1,   \
                             g_0_yyyyyyy_0_xxxzzz_0,  \
                             g_0_yyyyyyy_0_xxxzzz_1,  \
                             g_0_yyyyyyy_0_xxyyy_1,   \
                             g_0_yyyyyyy_0_xxyyyy_0,  \
                             g_0_yyyyyyy_0_xxyyyy_1,  \
                             g_0_yyyyyyy_0_xxyyyz_0,  \
                             g_0_yyyyyyy_0_xxyyyz_1,  \
                             g_0_yyyyyyy_0_xxyyz_1,   \
                             g_0_yyyyyyy_0_xxyyzz_0,  \
                             g_0_yyyyyyy_0_xxyyzz_1,  \
                             g_0_yyyyyyy_0_xxyzz_1,   \
                             g_0_yyyyyyy_0_xxyzzz_0,  \
                             g_0_yyyyyyy_0_xxyzzz_1,  \
                             g_0_yyyyyyy_0_xxzzz_1,   \
                             g_0_yyyyyyy_0_xxzzzz_0,  \
                             g_0_yyyyyyy_0_xxzzzz_1,  \
                             g_0_yyyyyyy_0_xyyyy_1,   \
                             g_0_yyyyyyy_0_xyyyyy_0,  \
                             g_0_yyyyyyy_0_xyyyyy_1,  \
                             g_0_yyyyyyy_0_xyyyyz_0,  \
                             g_0_yyyyyyy_0_xyyyyz_1,  \
                             g_0_yyyyyyy_0_xyyyz_1,   \
                             g_0_yyyyyyy_0_xyyyzz_0,  \
                             g_0_yyyyyyy_0_xyyyzz_1,  \
                             g_0_yyyyyyy_0_xyyzz_1,   \
                             g_0_yyyyyyy_0_xyyzzz_0,  \
                             g_0_yyyyyyy_0_xyyzzz_1,  \
                             g_0_yyyyyyy_0_xyzzz_1,   \
                             g_0_yyyyyyy_0_xyzzzz_0,  \
                             g_0_yyyyyyy_0_xyzzzz_1,  \
                             g_0_yyyyyyy_0_xzzzz_1,   \
                             g_0_yyyyyyy_0_xzzzzz_0,  \
                             g_0_yyyyyyy_0_xzzzzz_1,  \
                             g_0_yyyyyyy_0_yyyyy_1,   \
                             g_0_yyyyyyy_0_yyyyyy_0,  \
                             g_0_yyyyyyy_0_yyyyyy_1,  \
                             g_0_yyyyyyy_0_yyyyyz_0,  \
                             g_0_yyyyyyy_0_yyyyyz_1,  \
                             g_0_yyyyyyy_0_yyyyz_1,   \
                             g_0_yyyyyyy_0_yyyyzz_0,  \
                             g_0_yyyyyyy_0_yyyyzz_1,  \
                             g_0_yyyyyyy_0_yyyzz_1,   \
                             g_0_yyyyyyy_0_yyyzzz_0,  \
                             g_0_yyyyyyy_0_yyyzzz_1,  \
                             g_0_yyyyyyy_0_yyzzz_1,   \
                             g_0_yyyyyyy_0_yyzzzz_0,  \
                             g_0_yyyyyyy_0_yyzzzz_1,  \
                             g_0_yyyyyyy_0_yzzzz_1,   \
                             g_0_yyyyyyy_0_yzzzzz_0,  \
                             g_0_yyyyyyy_0_yzzzzz_1,  \
                             g_0_yyyyyyy_0_zzzzz_1,   \
                             g_0_yyyyyyy_0_zzzzzz_0,  \
                             g_0_yyyyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyyyy_0_xxxxxx_0, \
                             g_0_yyyyyyyy_0_xxxxxy_0, \
                             g_0_yyyyyyyy_0_xxxxxz_0, \
                             g_0_yyyyyyyy_0_xxxxyy_0, \
                             g_0_yyyyyyyy_0_xxxxyz_0, \
                             g_0_yyyyyyyy_0_xxxxzz_0, \
                             g_0_yyyyyyyy_0_xxxyyy_0, \
                             g_0_yyyyyyyy_0_xxxyyz_0, \
                             g_0_yyyyyyyy_0_xxxyzz_0, \
                             g_0_yyyyyyyy_0_xxxzzz_0, \
                             g_0_yyyyyyyy_0_xxyyyy_0, \
                             g_0_yyyyyyyy_0_xxyyyz_0, \
                             g_0_yyyyyyyy_0_xxyyzz_0, \
                             g_0_yyyyyyyy_0_xxyzzz_0, \
                             g_0_yyyyyyyy_0_xxzzzz_0, \
                             g_0_yyyyyyyy_0_xyyyyy_0, \
                             g_0_yyyyyyyy_0_xyyyyz_0, \
                             g_0_yyyyyyyy_0_xyyyzz_0, \
                             g_0_yyyyyyyy_0_xyyzzz_0, \
                             g_0_yyyyyyyy_0_xyzzzz_0, \
                             g_0_yyyyyyyy_0_xzzzzz_0, \
                             g_0_yyyyyyyy_0_yyyyyy_0, \
                             g_0_yyyyyyyy_0_yyyyyz_0, \
                             g_0_yyyyyyyy_0_yyyyzz_0, \
                             g_0_yyyyyyyy_0_yyyzzz_0, \
                             g_0_yyyyyyyy_0_yyzzzz_0, \
                             g_0_yyyyyyyy_0_yzzzzz_0, \
                             g_0_yyyyyyyy_0_zzzzzz_0, \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_xxxxxx_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxx_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxxxx_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxy_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxy_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxxz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxxz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxxxz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxyy_0[i] = 7.0 * g_0_yyyyyy_0_xxxxyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxyy_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyy_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xxxxyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxyz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxxzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyyy_0[i] = 7.0 * g_0_yyyyyy_0_xxxyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyyy_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyy_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xxxyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyyz_0[i] = 7.0 * g_0_yyyyyy_0_xxxyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxxzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyyy_0[i] = 7.0 * g_0_yyyyyy_0_xxyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyyy_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyy_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xxyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyyz_0[i] = 7.0 * g_0_yyyyyy_0_xxyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyzz_0[i] = 7.0 * g_0_yyyyyy_0_xxyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyzz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xxzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyyy_0[i] = 7.0 * g_0_yyyyyy_0_xyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyyy_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyy_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyyz_0[i] = 7.0 * g_0_yyyyyy_0_xyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyzz_0[i] = 7.0 * g_0_yyyyyy_0_xyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyzz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyzzz_0[i] = 7.0 * g_0_yyyyyy_0_xyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzzz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_xzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyyy_0[i] = 7.0 * g_0_yyyyyy_0_yyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyyy_1[i] * fti_ab_0 +
                                     6.0 * g_0_yyyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyy_0[i] * pb_y +
                                     g_0_yyyyyyy_0_yyyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyyz_0[i] = 7.0 * g_0_yyyyyy_0_yyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyzz_0[i] = 7.0 * g_0_yyyyyy_0_yyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyzz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyzzz_0[i] = 7.0 * g_0_yyyyyy_0_yyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyzzz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyzzzz_0[i] = 7.0 * g_0_yyyyyy_0_yyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzzzz_0[i] * pb_y +
                                     g_0_yyyyyyy_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_yzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_zzzzzz_0[i] = 7.0 * g_0_yyyyyy_0_zzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyyy_0_zzzzzz_0[i] * pb_y + g_0_yyyyyyy_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1036-1064 components of targeted buffer : SLSI

    auto g_0_yyyyyyyz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1036);

    auto g_0_yyyyyyyz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1037);

    auto g_0_yyyyyyyz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1038);

    auto g_0_yyyyyyyz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1039);

    auto g_0_yyyyyyyz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1040);

    auto g_0_yyyyyyyz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1041);

    auto g_0_yyyyyyyz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1042);

    auto g_0_yyyyyyyz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1043);

    auto g_0_yyyyyyyz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1044);

    auto g_0_yyyyyyyz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1045);

    auto g_0_yyyyyyyz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1046);

    auto g_0_yyyyyyyz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1047);

    auto g_0_yyyyyyyz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1048);

    auto g_0_yyyyyyyz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1049);

    auto g_0_yyyyyyyz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1050);

    auto g_0_yyyyyyyz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1051);

    auto g_0_yyyyyyyz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1052);

    auto g_0_yyyyyyyz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1053);

    auto g_0_yyyyyyyz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1054);

    auto g_0_yyyyyyyz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1055);

    auto g_0_yyyyyyyz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1056);

    auto g_0_yyyyyyyz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1057);

    auto g_0_yyyyyyyz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1058);

    auto g_0_yyyyyyyz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1059);

    auto g_0_yyyyyyyz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1060);

    auto g_0_yyyyyyyz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1061);

    auto g_0_yyyyyyyz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1062);

    auto g_0_yyyyyyyz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1063);

#pragma omp simd aligned(g_0_yyyyyyy_0_xxxxx_1,       \
                             g_0_yyyyyyy_0_xxxxxx_0,  \
                             g_0_yyyyyyy_0_xxxxxx_1,  \
                             g_0_yyyyyyy_0_xxxxxy_0,  \
                             g_0_yyyyyyy_0_xxxxxy_1,  \
                             g_0_yyyyyyy_0_xxxxxz_0,  \
                             g_0_yyyyyyy_0_xxxxxz_1,  \
                             g_0_yyyyyyy_0_xxxxy_1,   \
                             g_0_yyyyyyy_0_xxxxyy_0,  \
                             g_0_yyyyyyy_0_xxxxyy_1,  \
                             g_0_yyyyyyy_0_xxxxyz_0,  \
                             g_0_yyyyyyy_0_xxxxyz_1,  \
                             g_0_yyyyyyy_0_xxxxz_1,   \
                             g_0_yyyyyyy_0_xxxxzz_0,  \
                             g_0_yyyyyyy_0_xxxxzz_1,  \
                             g_0_yyyyyyy_0_xxxyy_1,   \
                             g_0_yyyyyyy_0_xxxyyy_0,  \
                             g_0_yyyyyyy_0_xxxyyy_1,  \
                             g_0_yyyyyyy_0_xxxyyz_0,  \
                             g_0_yyyyyyy_0_xxxyyz_1,  \
                             g_0_yyyyyyy_0_xxxyz_1,   \
                             g_0_yyyyyyy_0_xxxyzz_0,  \
                             g_0_yyyyyyy_0_xxxyzz_1,  \
                             g_0_yyyyyyy_0_xxxzz_1,   \
                             g_0_yyyyyyy_0_xxxzzz_0,  \
                             g_0_yyyyyyy_0_xxxzzz_1,  \
                             g_0_yyyyyyy_0_xxyyy_1,   \
                             g_0_yyyyyyy_0_xxyyyy_0,  \
                             g_0_yyyyyyy_0_xxyyyy_1,  \
                             g_0_yyyyyyy_0_xxyyyz_0,  \
                             g_0_yyyyyyy_0_xxyyyz_1,  \
                             g_0_yyyyyyy_0_xxyyz_1,   \
                             g_0_yyyyyyy_0_xxyyzz_0,  \
                             g_0_yyyyyyy_0_xxyyzz_1,  \
                             g_0_yyyyyyy_0_xxyzz_1,   \
                             g_0_yyyyyyy_0_xxyzzz_0,  \
                             g_0_yyyyyyy_0_xxyzzz_1,  \
                             g_0_yyyyyyy_0_xxzzz_1,   \
                             g_0_yyyyyyy_0_xxzzzz_0,  \
                             g_0_yyyyyyy_0_xxzzzz_1,  \
                             g_0_yyyyyyy_0_xyyyy_1,   \
                             g_0_yyyyyyy_0_xyyyyy_0,  \
                             g_0_yyyyyyy_0_xyyyyy_1,  \
                             g_0_yyyyyyy_0_xyyyyz_0,  \
                             g_0_yyyyyyy_0_xyyyyz_1,  \
                             g_0_yyyyyyy_0_xyyyz_1,   \
                             g_0_yyyyyyy_0_xyyyzz_0,  \
                             g_0_yyyyyyy_0_xyyyzz_1,  \
                             g_0_yyyyyyy_0_xyyzz_1,   \
                             g_0_yyyyyyy_0_xyyzzz_0,  \
                             g_0_yyyyyyy_0_xyyzzz_1,  \
                             g_0_yyyyyyy_0_xyzzz_1,   \
                             g_0_yyyyyyy_0_xyzzzz_0,  \
                             g_0_yyyyyyy_0_xyzzzz_1,  \
                             g_0_yyyyyyy_0_xzzzz_1,   \
                             g_0_yyyyyyy_0_xzzzzz_0,  \
                             g_0_yyyyyyy_0_xzzzzz_1,  \
                             g_0_yyyyyyy_0_yyyyy_1,   \
                             g_0_yyyyyyy_0_yyyyyy_0,  \
                             g_0_yyyyyyy_0_yyyyyy_1,  \
                             g_0_yyyyyyy_0_yyyyyz_0,  \
                             g_0_yyyyyyy_0_yyyyyz_1,  \
                             g_0_yyyyyyy_0_yyyyz_1,   \
                             g_0_yyyyyyy_0_yyyyzz_0,  \
                             g_0_yyyyyyy_0_yyyyzz_1,  \
                             g_0_yyyyyyy_0_yyyzz_1,   \
                             g_0_yyyyyyy_0_yyyzzz_0,  \
                             g_0_yyyyyyy_0_yyyzzz_1,  \
                             g_0_yyyyyyy_0_yyzzz_1,   \
                             g_0_yyyyyyy_0_yyzzzz_0,  \
                             g_0_yyyyyyy_0_yyzzzz_1,  \
                             g_0_yyyyyyy_0_yzzzz_1,   \
                             g_0_yyyyyyy_0_yzzzzz_0,  \
                             g_0_yyyyyyy_0_yzzzzz_1,  \
                             g_0_yyyyyyy_0_zzzzz_1,   \
                             g_0_yyyyyyy_0_zzzzzz_0,  \
                             g_0_yyyyyyy_0_zzzzzz_1,  \
                             g_0_yyyyyyyz_0_xxxxxx_0, \
                             g_0_yyyyyyyz_0_xxxxxy_0, \
                             g_0_yyyyyyyz_0_xxxxxz_0, \
                             g_0_yyyyyyyz_0_xxxxyy_0, \
                             g_0_yyyyyyyz_0_xxxxyz_0, \
                             g_0_yyyyyyyz_0_xxxxzz_0, \
                             g_0_yyyyyyyz_0_xxxyyy_0, \
                             g_0_yyyyyyyz_0_xxxyyz_0, \
                             g_0_yyyyyyyz_0_xxxyzz_0, \
                             g_0_yyyyyyyz_0_xxxzzz_0, \
                             g_0_yyyyyyyz_0_xxyyyy_0, \
                             g_0_yyyyyyyz_0_xxyyyz_0, \
                             g_0_yyyyyyyz_0_xxyyzz_0, \
                             g_0_yyyyyyyz_0_xxyzzz_0, \
                             g_0_yyyyyyyz_0_xxzzzz_0, \
                             g_0_yyyyyyyz_0_xyyyyy_0, \
                             g_0_yyyyyyyz_0_xyyyyz_0, \
                             g_0_yyyyyyyz_0_xyyyzz_0, \
                             g_0_yyyyyyyz_0_xyyzzz_0, \
                             g_0_yyyyyyyz_0_xyzzzz_0, \
                             g_0_yyyyyyyz_0_xzzzzz_0, \
                             g_0_yyyyyyyz_0_yyyyyy_0, \
                             g_0_yyyyyyyz_0_yyyyyz_0, \
                             g_0_yyyyyyyz_0_yyyyzz_0, \
                             g_0_yyyyyyyz_0_yyyzzz_0, \
                             g_0_yyyyyyyz_0_yyzzzz_0, \
                             g_0_yyyyyyyz_0_yzzzzz_0, \
                             g_0_yyyyyyyz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_xxxxxx_0[i] = g_0_yyyyyyy_0_xxxxxx_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxy_0[i] = g_0_yyyyyyy_0_xxxxxy_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxxz_0[i] = g_0_yyyyyyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxxz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxyy_0[i] = g_0_yyyyyyy_0_xxxxyy_0[i] * pb_z + g_0_yyyyyyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxyz_0[i] = g_0_yyyyyyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxyz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyyy_0[i] = g_0_yyyyyyy_0_xxxyyy_0[i] * pb_z + g_0_yyyyyyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyyz_0[i] = g_0_yyyyyyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyyz_0[i] * pb_z + g_0_yyyyyyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxzzz_0[i] =
            3.0 * g_0_yyyyyyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyyy_0[i] = g_0_yyyyyyy_0_xxyyyy_0[i] * pb_z + g_0_yyyyyyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyyz_0[i] = g_0_yyyyyyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyyz_0[i] * pb_z + g_0_yyyyyyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyzz_0[i] * pb_z + g_0_yyyyyyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyzzz_0[i] =
            3.0 * g_0_yyyyyyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxzzzz_0[i] =
            4.0 * g_0_yyyyyyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyyy_0[i] = g_0_yyyyyyy_0_xyyyyy_0[i] * pb_z + g_0_yyyyyyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyyz_0[i] = g_0_yyyyyyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyyz_0[i] * pb_z + g_0_yyyyyyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyzz_0[i] =
            2.0 * g_0_yyyyyyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyzz_0[i] * pb_z + g_0_yyyyyyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyzzz_0[i] =
            3.0 * g_0_yyyyyyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzzz_0[i] * pb_z + g_0_yyyyyyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyzzzz_0[i] =
            4.0 * g_0_yyyyyyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xzzzzz_0[i] =
            5.0 * g_0_yyyyyyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyyy_0[i] = g_0_yyyyyyy_0_yyyyyy_0[i] * pb_z + g_0_yyyyyyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyyz_0[i] = g_0_yyyyyyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyyz_0[i] * pb_z + g_0_yyyyyyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyzz_0[i] =
            2.0 * g_0_yyyyyyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyzz_0[i] * pb_z + g_0_yyyyyyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyzzz_0[i] =
            3.0 * g_0_yyyyyyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyzzz_0[i] * pb_z + g_0_yyyyyyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyzzzz_0[i] =
            4.0 * g_0_yyyyyyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzzzz_0[i] * pb_z + g_0_yyyyyyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yzzzzz_0[i] =
            5.0 * g_0_yyyyyyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_zzzzzz_0[i] =
            6.0 * g_0_yyyyyyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_zzzzzz_0[i] * pb_z + g_0_yyyyyyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 1064-1092 components of targeted buffer : SLSI

    auto g_0_yyyyyyzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1064);

    auto g_0_yyyyyyzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1065);

    auto g_0_yyyyyyzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1066);

    auto g_0_yyyyyyzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1067);

    auto g_0_yyyyyyzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1068);

    auto g_0_yyyyyyzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1069);

    auto g_0_yyyyyyzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1070);

    auto g_0_yyyyyyzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1071);

    auto g_0_yyyyyyzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1072);

    auto g_0_yyyyyyzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1073);

    auto g_0_yyyyyyzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1074);

    auto g_0_yyyyyyzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1075);

    auto g_0_yyyyyyzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1076);

    auto g_0_yyyyyyzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1077);

    auto g_0_yyyyyyzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1078);

    auto g_0_yyyyyyzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1079);

    auto g_0_yyyyyyzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1080);

    auto g_0_yyyyyyzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1081);

    auto g_0_yyyyyyzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1082);

    auto g_0_yyyyyyzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1083);

    auto g_0_yyyyyyzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1084);

    auto g_0_yyyyyyzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1085);

    auto g_0_yyyyyyzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1086);

    auto g_0_yyyyyyzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1087);

    auto g_0_yyyyyyzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1088);

    auto g_0_yyyyyyzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1089);

    auto g_0_yyyyyyzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1090);

    auto g_0_yyyyyyzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1091);

#pragma omp simd aligned(g_0_yyyyyy_0_xxxxxy_0,       \
                             g_0_yyyyyy_0_xxxxxy_1,   \
                             g_0_yyyyyy_0_xxxxyy_0,   \
                             g_0_yyyyyy_0_xxxxyy_1,   \
                             g_0_yyyyyy_0_xxxyyy_0,   \
                             g_0_yyyyyy_0_xxxyyy_1,   \
                             g_0_yyyyyy_0_xxyyyy_0,   \
                             g_0_yyyyyy_0_xxyyyy_1,   \
                             g_0_yyyyyy_0_xyyyyy_0,   \
                             g_0_yyyyyy_0_xyyyyy_1,   \
                             g_0_yyyyyy_0_yyyyyy_0,   \
                             g_0_yyyyyy_0_yyyyyy_1,   \
                             g_0_yyyyyyz_0_xxxxxy_0,  \
                             g_0_yyyyyyz_0_xxxxxy_1,  \
                             g_0_yyyyyyz_0_xxxxyy_0,  \
                             g_0_yyyyyyz_0_xxxxyy_1,  \
                             g_0_yyyyyyz_0_xxxyyy_0,  \
                             g_0_yyyyyyz_0_xxxyyy_1,  \
                             g_0_yyyyyyz_0_xxyyyy_0,  \
                             g_0_yyyyyyz_0_xxyyyy_1,  \
                             g_0_yyyyyyz_0_xyyyyy_0,  \
                             g_0_yyyyyyz_0_xyyyyy_1,  \
                             g_0_yyyyyyz_0_yyyyyy_0,  \
                             g_0_yyyyyyz_0_yyyyyy_1,  \
                             g_0_yyyyyyzz_0_xxxxxx_0, \
                             g_0_yyyyyyzz_0_xxxxxy_0, \
                             g_0_yyyyyyzz_0_xxxxxz_0, \
                             g_0_yyyyyyzz_0_xxxxyy_0, \
                             g_0_yyyyyyzz_0_xxxxyz_0, \
                             g_0_yyyyyyzz_0_xxxxzz_0, \
                             g_0_yyyyyyzz_0_xxxyyy_0, \
                             g_0_yyyyyyzz_0_xxxyyz_0, \
                             g_0_yyyyyyzz_0_xxxyzz_0, \
                             g_0_yyyyyyzz_0_xxxzzz_0, \
                             g_0_yyyyyyzz_0_xxyyyy_0, \
                             g_0_yyyyyyzz_0_xxyyyz_0, \
                             g_0_yyyyyyzz_0_xxyyzz_0, \
                             g_0_yyyyyyzz_0_xxyzzz_0, \
                             g_0_yyyyyyzz_0_xxzzzz_0, \
                             g_0_yyyyyyzz_0_xyyyyy_0, \
                             g_0_yyyyyyzz_0_xyyyyz_0, \
                             g_0_yyyyyyzz_0_xyyyzz_0, \
                             g_0_yyyyyyzz_0_xyyzzz_0, \
                             g_0_yyyyyyzz_0_xyzzzz_0, \
                             g_0_yyyyyyzz_0_xzzzzz_0, \
                             g_0_yyyyyyzz_0_yyyyyy_0, \
                             g_0_yyyyyyzz_0_yyyyyz_0, \
                             g_0_yyyyyyzz_0_yyyyzz_0, \
                             g_0_yyyyyyzz_0_yyyzzz_0, \
                             g_0_yyyyyyzz_0_yyzzzz_0, \
                             g_0_yyyyyyzz_0_yzzzzz_0, \
                             g_0_yyyyyyzz_0_zzzzzz_0, \
                             g_0_yyyyyzz_0_xxxxxx_0,  \
                             g_0_yyyyyzz_0_xxxxxx_1,  \
                             g_0_yyyyyzz_0_xxxxxz_0,  \
                             g_0_yyyyyzz_0_xxxxxz_1,  \
                             g_0_yyyyyzz_0_xxxxyz_0,  \
                             g_0_yyyyyzz_0_xxxxyz_1,  \
                             g_0_yyyyyzz_0_xxxxz_1,   \
                             g_0_yyyyyzz_0_xxxxzz_0,  \
                             g_0_yyyyyzz_0_xxxxzz_1,  \
                             g_0_yyyyyzz_0_xxxyyz_0,  \
                             g_0_yyyyyzz_0_xxxyyz_1,  \
                             g_0_yyyyyzz_0_xxxyz_1,   \
                             g_0_yyyyyzz_0_xxxyzz_0,  \
                             g_0_yyyyyzz_0_xxxyzz_1,  \
                             g_0_yyyyyzz_0_xxxzz_1,   \
                             g_0_yyyyyzz_0_xxxzzz_0,  \
                             g_0_yyyyyzz_0_xxxzzz_1,  \
                             g_0_yyyyyzz_0_xxyyyz_0,  \
                             g_0_yyyyyzz_0_xxyyyz_1,  \
                             g_0_yyyyyzz_0_xxyyz_1,   \
                             g_0_yyyyyzz_0_xxyyzz_0,  \
                             g_0_yyyyyzz_0_xxyyzz_1,  \
                             g_0_yyyyyzz_0_xxyzz_1,   \
                             g_0_yyyyyzz_0_xxyzzz_0,  \
                             g_0_yyyyyzz_0_xxyzzz_1,  \
                             g_0_yyyyyzz_0_xxzzz_1,   \
                             g_0_yyyyyzz_0_xxzzzz_0,  \
                             g_0_yyyyyzz_0_xxzzzz_1,  \
                             g_0_yyyyyzz_0_xyyyyz_0,  \
                             g_0_yyyyyzz_0_xyyyyz_1,  \
                             g_0_yyyyyzz_0_xyyyz_1,   \
                             g_0_yyyyyzz_0_xyyyzz_0,  \
                             g_0_yyyyyzz_0_xyyyzz_1,  \
                             g_0_yyyyyzz_0_xyyzz_1,   \
                             g_0_yyyyyzz_0_xyyzzz_0,  \
                             g_0_yyyyyzz_0_xyyzzz_1,  \
                             g_0_yyyyyzz_0_xyzzz_1,   \
                             g_0_yyyyyzz_0_xyzzzz_0,  \
                             g_0_yyyyyzz_0_xyzzzz_1,  \
                             g_0_yyyyyzz_0_xzzzz_1,   \
                             g_0_yyyyyzz_0_xzzzzz_0,  \
                             g_0_yyyyyzz_0_xzzzzz_1,  \
                             g_0_yyyyyzz_0_yyyyyz_0,  \
                             g_0_yyyyyzz_0_yyyyyz_1,  \
                             g_0_yyyyyzz_0_yyyyz_1,   \
                             g_0_yyyyyzz_0_yyyyzz_0,  \
                             g_0_yyyyyzz_0_yyyyzz_1,  \
                             g_0_yyyyyzz_0_yyyzz_1,   \
                             g_0_yyyyyzz_0_yyyzzz_0,  \
                             g_0_yyyyyzz_0_yyyzzz_1,  \
                             g_0_yyyyyzz_0_yyzzz_1,   \
                             g_0_yyyyyzz_0_yyzzzz_0,  \
                             g_0_yyyyyzz_0_yyzzzz_1,  \
                             g_0_yyyyyzz_0_yzzzz_1,   \
                             g_0_yyyyyzz_0_yzzzzz_0,  \
                             g_0_yyyyyzz_0_yzzzzz_1,  \
                             g_0_yyyyyzz_0_zzzzz_1,   \
                             g_0_yyyyyzz_0_zzzzzz_0,  \
                             g_0_yyyyyzz_0_zzzzzz_1,  \
                             g_0_yyyyzz_0_xxxxxx_0,   \
                             g_0_yyyyzz_0_xxxxxx_1,   \
                             g_0_yyyyzz_0_xxxxxz_0,   \
                             g_0_yyyyzz_0_xxxxxz_1,   \
                             g_0_yyyyzz_0_xxxxyz_0,   \
                             g_0_yyyyzz_0_xxxxyz_1,   \
                             g_0_yyyyzz_0_xxxxzz_0,   \
                             g_0_yyyyzz_0_xxxxzz_1,   \
                             g_0_yyyyzz_0_xxxyyz_0,   \
                             g_0_yyyyzz_0_xxxyyz_1,   \
                             g_0_yyyyzz_0_xxxyzz_0,   \
                             g_0_yyyyzz_0_xxxyzz_1,   \
                             g_0_yyyyzz_0_xxxzzz_0,   \
                             g_0_yyyyzz_0_xxxzzz_1,   \
                             g_0_yyyyzz_0_xxyyyz_0,   \
                             g_0_yyyyzz_0_xxyyyz_1,   \
                             g_0_yyyyzz_0_xxyyzz_0,   \
                             g_0_yyyyzz_0_xxyyzz_1,   \
                             g_0_yyyyzz_0_xxyzzz_0,   \
                             g_0_yyyyzz_0_xxyzzz_1,   \
                             g_0_yyyyzz_0_xxzzzz_0,   \
                             g_0_yyyyzz_0_xxzzzz_1,   \
                             g_0_yyyyzz_0_xyyyyz_0,   \
                             g_0_yyyyzz_0_xyyyyz_1,   \
                             g_0_yyyyzz_0_xyyyzz_0,   \
                             g_0_yyyyzz_0_xyyyzz_1,   \
                             g_0_yyyyzz_0_xyyzzz_0,   \
                             g_0_yyyyzz_0_xyyzzz_1,   \
                             g_0_yyyyzz_0_xyzzzz_0,   \
                             g_0_yyyyzz_0_xyzzzz_1,   \
                             g_0_yyyyzz_0_xzzzzz_0,   \
                             g_0_yyyyzz_0_xzzzzz_1,   \
                             g_0_yyyyzz_0_yyyyyz_0,   \
                             g_0_yyyyzz_0_yyyyyz_1,   \
                             g_0_yyyyzz_0_yyyyzz_0,   \
                             g_0_yyyyzz_0_yyyyzz_1,   \
                             g_0_yyyyzz_0_yyyzzz_0,   \
                             g_0_yyyyzz_0_yyyzzz_1,   \
                             g_0_yyyyzz_0_yyzzzz_0,   \
                             g_0_yyyyzz_0_yyzzzz_1,   \
                             g_0_yyyyzz_0_yzzzzz_0,   \
                             g_0_yyyyzz_0_yzzzzz_1,   \
                             g_0_yyyyzz_0_zzzzzz_0,   \
                             g_0_yyyyzz_0_zzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_xxxxxx_0[i] = 5.0 * g_0_yyyyzz_0_xxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxxxx_0[i] * pb_y + g_0_yyyyyzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxxy_0[i] = g_0_yyyyyy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxxxy_0[i] * pb_z +
                                     g_0_yyyyyyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxxxz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxxxz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxyy_0[i] = g_0_yyyyyy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxxyy_0[i] * pb_z +
                                     g_0_yyyyyyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxxyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxyz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxxzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxyyy_0[i] = g_0_yyyyyy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxyyy_0[i] * pb_z +
                                     g_0_yyyyyyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxyyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyyz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxyzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyyyy_0[i] = g_0_yyyyyy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxyyyy_0[i] * pb_z +
                                     g_0_yyyyyyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxyyyz_0[i] = 5.0 * g_0_yyyyzz_0_xxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyyz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyyzz_0[i] = 5.0 * g_0_yyyyzz_0_xxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyzz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyyyy_0[i] = g_0_yyyyyy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xyyyyy_0[i] * pb_z +
                                     g_0_yyyyyyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xyyyyz_0[i] = 5.0 * g_0_yyyyzz_0_xyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyyzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyyz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyyzz_0[i] = 5.0 * g_0_yyyyzz_0_xyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyzz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyzzz_0[i] = 5.0 * g_0_yyyyzz_0_xyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyzzz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyyyy_0[i] = g_0_yyyyyy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_yyyyyy_0[i] * pb_z +
                                     g_0_yyyyyyz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_yyyyyz_0[i] = 5.0 * g_0_yyyyzz_0_yyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyyzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyyyz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyyzz_0[i] = 5.0 * g_0_yyyyzz_0_yyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyyzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyyzz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyzzz_0[i] = 5.0 * g_0_yyyyzz_0_yyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyyzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyzzz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyzzzz_0[i] = 5.0 * g_0_yyyyzz_0_yyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyyzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyzzzz_0[i] * pb_y +
                                     g_0_yyyyyzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_yzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_zzzzzz_0[i] = 5.0 * g_0_yyyyzz_0_zzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_zzzzzz_0[i] * pb_y + g_0_yyyyyzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1092-1120 components of targeted buffer : SLSI

    auto g_0_yyyyyzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1092);

    auto g_0_yyyyyzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1093);

    auto g_0_yyyyyzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1094);

    auto g_0_yyyyyzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1095);

    auto g_0_yyyyyzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1096);

    auto g_0_yyyyyzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1097);

    auto g_0_yyyyyzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1098);

    auto g_0_yyyyyzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1099);

    auto g_0_yyyyyzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1100);

    auto g_0_yyyyyzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1101);

    auto g_0_yyyyyzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1102);

    auto g_0_yyyyyzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1103);

    auto g_0_yyyyyzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1104);

    auto g_0_yyyyyzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1105);

    auto g_0_yyyyyzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1106);

    auto g_0_yyyyyzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1107);

    auto g_0_yyyyyzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1108);

    auto g_0_yyyyyzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1109);

    auto g_0_yyyyyzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1110);

    auto g_0_yyyyyzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1111);

    auto g_0_yyyyyzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1112);

    auto g_0_yyyyyzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1113);

    auto g_0_yyyyyzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1114);

    auto g_0_yyyyyzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1115);

    auto g_0_yyyyyzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1116);

    auto g_0_yyyyyzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1117);

    auto g_0_yyyyyzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1118);

    auto g_0_yyyyyzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1119);

#pragma omp simd aligned(g_0_yyyyyz_0_xxxxxy_0,       \
                             g_0_yyyyyz_0_xxxxxy_1,   \
                             g_0_yyyyyz_0_xxxxyy_0,   \
                             g_0_yyyyyz_0_xxxxyy_1,   \
                             g_0_yyyyyz_0_xxxyyy_0,   \
                             g_0_yyyyyz_0_xxxyyy_1,   \
                             g_0_yyyyyz_0_xxyyyy_0,   \
                             g_0_yyyyyz_0_xxyyyy_1,   \
                             g_0_yyyyyz_0_xyyyyy_0,   \
                             g_0_yyyyyz_0_xyyyyy_1,   \
                             g_0_yyyyyz_0_yyyyyy_0,   \
                             g_0_yyyyyz_0_yyyyyy_1,   \
                             g_0_yyyyyzz_0_xxxxxy_0,  \
                             g_0_yyyyyzz_0_xxxxxy_1,  \
                             g_0_yyyyyzz_0_xxxxyy_0,  \
                             g_0_yyyyyzz_0_xxxxyy_1,  \
                             g_0_yyyyyzz_0_xxxyyy_0,  \
                             g_0_yyyyyzz_0_xxxyyy_1,  \
                             g_0_yyyyyzz_0_xxyyyy_0,  \
                             g_0_yyyyyzz_0_xxyyyy_1,  \
                             g_0_yyyyyzz_0_xyyyyy_0,  \
                             g_0_yyyyyzz_0_xyyyyy_1,  \
                             g_0_yyyyyzz_0_yyyyyy_0,  \
                             g_0_yyyyyzz_0_yyyyyy_1,  \
                             g_0_yyyyyzzz_0_xxxxxx_0, \
                             g_0_yyyyyzzz_0_xxxxxy_0, \
                             g_0_yyyyyzzz_0_xxxxxz_0, \
                             g_0_yyyyyzzz_0_xxxxyy_0, \
                             g_0_yyyyyzzz_0_xxxxyz_0, \
                             g_0_yyyyyzzz_0_xxxxzz_0, \
                             g_0_yyyyyzzz_0_xxxyyy_0, \
                             g_0_yyyyyzzz_0_xxxyyz_0, \
                             g_0_yyyyyzzz_0_xxxyzz_0, \
                             g_0_yyyyyzzz_0_xxxzzz_0, \
                             g_0_yyyyyzzz_0_xxyyyy_0, \
                             g_0_yyyyyzzz_0_xxyyyz_0, \
                             g_0_yyyyyzzz_0_xxyyzz_0, \
                             g_0_yyyyyzzz_0_xxyzzz_0, \
                             g_0_yyyyyzzz_0_xxzzzz_0, \
                             g_0_yyyyyzzz_0_xyyyyy_0, \
                             g_0_yyyyyzzz_0_xyyyyz_0, \
                             g_0_yyyyyzzz_0_xyyyzz_0, \
                             g_0_yyyyyzzz_0_xyyzzz_0, \
                             g_0_yyyyyzzz_0_xyzzzz_0, \
                             g_0_yyyyyzzz_0_xzzzzz_0, \
                             g_0_yyyyyzzz_0_yyyyyy_0, \
                             g_0_yyyyyzzz_0_yyyyyz_0, \
                             g_0_yyyyyzzz_0_yyyyzz_0, \
                             g_0_yyyyyzzz_0_yyyzzz_0, \
                             g_0_yyyyyzzz_0_yyzzzz_0, \
                             g_0_yyyyyzzz_0_yzzzzz_0, \
                             g_0_yyyyyzzz_0_zzzzzz_0, \
                             g_0_yyyyzzz_0_xxxxxx_0,  \
                             g_0_yyyyzzz_0_xxxxxx_1,  \
                             g_0_yyyyzzz_0_xxxxxz_0,  \
                             g_0_yyyyzzz_0_xxxxxz_1,  \
                             g_0_yyyyzzz_0_xxxxyz_0,  \
                             g_0_yyyyzzz_0_xxxxyz_1,  \
                             g_0_yyyyzzz_0_xxxxz_1,   \
                             g_0_yyyyzzz_0_xxxxzz_0,  \
                             g_0_yyyyzzz_0_xxxxzz_1,  \
                             g_0_yyyyzzz_0_xxxyyz_0,  \
                             g_0_yyyyzzz_0_xxxyyz_1,  \
                             g_0_yyyyzzz_0_xxxyz_1,   \
                             g_0_yyyyzzz_0_xxxyzz_0,  \
                             g_0_yyyyzzz_0_xxxyzz_1,  \
                             g_0_yyyyzzz_0_xxxzz_1,   \
                             g_0_yyyyzzz_0_xxxzzz_0,  \
                             g_0_yyyyzzz_0_xxxzzz_1,  \
                             g_0_yyyyzzz_0_xxyyyz_0,  \
                             g_0_yyyyzzz_0_xxyyyz_1,  \
                             g_0_yyyyzzz_0_xxyyz_1,   \
                             g_0_yyyyzzz_0_xxyyzz_0,  \
                             g_0_yyyyzzz_0_xxyyzz_1,  \
                             g_0_yyyyzzz_0_xxyzz_1,   \
                             g_0_yyyyzzz_0_xxyzzz_0,  \
                             g_0_yyyyzzz_0_xxyzzz_1,  \
                             g_0_yyyyzzz_0_xxzzz_1,   \
                             g_0_yyyyzzz_0_xxzzzz_0,  \
                             g_0_yyyyzzz_0_xxzzzz_1,  \
                             g_0_yyyyzzz_0_xyyyyz_0,  \
                             g_0_yyyyzzz_0_xyyyyz_1,  \
                             g_0_yyyyzzz_0_xyyyz_1,   \
                             g_0_yyyyzzz_0_xyyyzz_0,  \
                             g_0_yyyyzzz_0_xyyyzz_1,  \
                             g_0_yyyyzzz_0_xyyzz_1,   \
                             g_0_yyyyzzz_0_xyyzzz_0,  \
                             g_0_yyyyzzz_0_xyyzzz_1,  \
                             g_0_yyyyzzz_0_xyzzz_1,   \
                             g_0_yyyyzzz_0_xyzzzz_0,  \
                             g_0_yyyyzzz_0_xyzzzz_1,  \
                             g_0_yyyyzzz_0_xzzzz_1,   \
                             g_0_yyyyzzz_0_xzzzzz_0,  \
                             g_0_yyyyzzz_0_xzzzzz_1,  \
                             g_0_yyyyzzz_0_yyyyyz_0,  \
                             g_0_yyyyzzz_0_yyyyyz_1,  \
                             g_0_yyyyzzz_0_yyyyz_1,   \
                             g_0_yyyyzzz_0_yyyyzz_0,  \
                             g_0_yyyyzzz_0_yyyyzz_1,  \
                             g_0_yyyyzzz_0_yyyzz_1,   \
                             g_0_yyyyzzz_0_yyyzzz_0,  \
                             g_0_yyyyzzz_0_yyyzzz_1,  \
                             g_0_yyyyzzz_0_yyzzz_1,   \
                             g_0_yyyyzzz_0_yyzzzz_0,  \
                             g_0_yyyyzzz_0_yyzzzz_1,  \
                             g_0_yyyyzzz_0_yzzzz_1,   \
                             g_0_yyyyzzz_0_yzzzzz_0,  \
                             g_0_yyyyzzz_0_yzzzzz_1,  \
                             g_0_yyyyzzz_0_zzzzz_1,   \
                             g_0_yyyyzzz_0_zzzzzz_0,  \
                             g_0_yyyyzzz_0_zzzzzz_1,  \
                             g_0_yyyzzz_0_xxxxxx_0,   \
                             g_0_yyyzzz_0_xxxxxx_1,   \
                             g_0_yyyzzz_0_xxxxxz_0,   \
                             g_0_yyyzzz_0_xxxxxz_1,   \
                             g_0_yyyzzz_0_xxxxyz_0,   \
                             g_0_yyyzzz_0_xxxxyz_1,   \
                             g_0_yyyzzz_0_xxxxzz_0,   \
                             g_0_yyyzzz_0_xxxxzz_1,   \
                             g_0_yyyzzz_0_xxxyyz_0,   \
                             g_0_yyyzzz_0_xxxyyz_1,   \
                             g_0_yyyzzz_0_xxxyzz_0,   \
                             g_0_yyyzzz_0_xxxyzz_1,   \
                             g_0_yyyzzz_0_xxxzzz_0,   \
                             g_0_yyyzzz_0_xxxzzz_1,   \
                             g_0_yyyzzz_0_xxyyyz_0,   \
                             g_0_yyyzzz_0_xxyyyz_1,   \
                             g_0_yyyzzz_0_xxyyzz_0,   \
                             g_0_yyyzzz_0_xxyyzz_1,   \
                             g_0_yyyzzz_0_xxyzzz_0,   \
                             g_0_yyyzzz_0_xxyzzz_1,   \
                             g_0_yyyzzz_0_xxzzzz_0,   \
                             g_0_yyyzzz_0_xxzzzz_1,   \
                             g_0_yyyzzz_0_xyyyyz_0,   \
                             g_0_yyyzzz_0_xyyyyz_1,   \
                             g_0_yyyzzz_0_xyyyzz_0,   \
                             g_0_yyyzzz_0_xyyyzz_1,   \
                             g_0_yyyzzz_0_xyyzzz_0,   \
                             g_0_yyyzzz_0_xyyzzz_1,   \
                             g_0_yyyzzz_0_xyzzzz_0,   \
                             g_0_yyyzzz_0_xyzzzz_1,   \
                             g_0_yyyzzz_0_xzzzzz_0,   \
                             g_0_yyyzzz_0_xzzzzz_1,   \
                             g_0_yyyzzz_0_yyyyyz_0,   \
                             g_0_yyyzzz_0_yyyyyz_1,   \
                             g_0_yyyzzz_0_yyyyzz_0,   \
                             g_0_yyyzzz_0_yyyyzz_1,   \
                             g_0_yyyzzz_0_yyyzzz_0,   \
                             g_0_yyyzzz_0_yyyzzz_1,   \
                             g_0_yyyzzz_0_yyzzzz_0,   \
                             g_0_yyyzzz_0_yyzzzz_1,   \
                             g_0_yyyzzz_0_yzzzzz_0,   \
                             g_0_yyyzzz_0_yzzzzz_1,   \
                             g_0_yyyzzz_0_zzzzzz_0,   \
                             g_0_yyyzzz_0_zzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_xxxxxx_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxxxx_0[i] * pb_y + g_0_yyyyzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxxy_0[i] = 2.0 * g_0_yyyyyz_0_xxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxxxy_0[i] * pb_z + g_0_yyyyyzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxxxz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxxxz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxyy_0[i] = 2.0 * g_0_yyyyyz_0_xxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxxyy_0[i] * pb_z + g_0_yyyyyzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxxyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxyz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxxzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxyyy_0[i] = 2.0 * g_0_yyyyyz_0_xxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxxyyy_0[i] * pb_z + g_0_yyyyyzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyyz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyyyy_0[i] = 2.0 * g_0_yyyyyz_0_xxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xxyyyy_0[i] * pb_z + g_0_yyyyyzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxyyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyyz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyzz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyyyy_0[i] = 2.0 * g_0_yyyyyz_0_xyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_xyyyyy_0[i] * pb_z + g_0_yyyyyzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xyyyyz_0[i] = 4.0 * g_0_yyyzzz_0_xyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyyz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyyzz_0[i] = 4.0 * g_0_yyyzzz_0_xyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyzz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyzzz_0[i] = 4.0 * g_0_yyyzzz_0_xyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyzzz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyyyy_0[i] = 2.0 * g_0_yyyyyz_0_yyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyyzz_0_yyyyyy_0[i] * pb_z + g_0_yyyyyzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_yyyyyz_0[i] = 4.0 * g_0_yyyzzz_0_yyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyyzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyyyz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyyzz_0[i] = 4.0 * g_0_yyyzzz_0_yyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyyzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyyzz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyzzz_0[i] = 4.0 * g_0_yyyzzz_0_yyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyyzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyzzz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyzzzz_0[i] = 4.0 * g_0_yyyzzz_0_yyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyyzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyzzzz_0[i] * pb_y +
                                     g_0_yyyyzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_yzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_zzzzzz_0[i] = 4.0 * g_0_yyyzzz_0_zzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_zzzzzz_0[i] * pb_y + g_0_yyyyzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1120-1148 components of targeted buffer : SLSI

    auto g_0_yyyyzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1120);

    auto g_0_yyyyzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1121);

    auto g_0_yyyyzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1122);

    auto g_0_yyyyzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1123);

    auto g_0_yyyyzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1124);

    auto g_0_yyyyzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1125);

    auto g_0_yyyyzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1126);

    auto g_0_yyyyzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1127);

    auto g_0_yyyyzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1128);

    auto g_0_yyyyzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1129);

    auto g_0_yyyyzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1130);

    auto g_0_yyyyzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1131);

    auto g_0_yyyyzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1132);

    auto g_0_yyyyzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1133);

    auto g_0_yyyyzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1134);

    auto g_0_yyyyzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1135);

    auto g_0_yyyyzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1136);

    auto g_0_yyyyzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1137);

    auto g_0_yyyyzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1138);

    auto g_0_yyyyzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1139);

    auto g_0_yyyyzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1140);

    auto g_0_yyyyzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1141);

    auto g_0_yyyyzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1142);

    auto g_0_yyyyzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1143);

    auto g_0_yyyyzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1144);

    auto g_0_yyyyzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1145);

    auto g_0_yyyyzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1146);

    auto g_0_yyyyzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1147);

#pragma omp simd aligned(g_0_yyyyzz_0_xxxxxy_0,       \
                             g_0_yyyyzz_0_xxxxxy_1,   \
                             g_0_yyyyzz_0_xxxxyy_0,   \
                             g_0_yyyyzz_0_xxxxyy_1,   \
                             g_0_yyyyzz_0_xxxyyy_0,   \
                             g_0_yyyyzz_0_xxxyyy_1,   \
                             g_0_yyyyzz_0_xxyyyy_0,   \
                             g_0_yyyyzz_0_xxyyyy_1,   \
                             g_0_yyyyzz_0_xyyyyy_0,   \
                             g_0_yyyyzz_0_xyyyyy_1,   \
                             g_0_yyyyzz_0_yyyyyy_0,   \
                             g_0_yyyyzz_0_yyyyyy_1,   \
                             g_0_yyyyzzz_0_xxxxxy_0,  \
                             g_0_yyyyzzz_0_xxxxxy_1,  \
                             g_0_yyyyzzz_0_xxxxyy_0,  \
                             g_0_yyyyzzz_0_xxxxyy_1,  \
                             g_0_yyyyzzz_0_xxxyyy_0,  \
                             g_0_yyyyzzz_0_xxxyyy_1,  \
                             g_0_yyyyzzz_0_xxyyyy_0,  \
                             g_0_yyyyzzz_0_xxyyyy_1,  \
                             g_0_yyyyzzz_0_xyyyyy_0,  \
                             g_0_yyyyzzz_0_xyyyyy_1,  \
                             g_0_yyyyzzz_0_yyyyyy_0,  \
                             g_0_yyyyzzz_0_yyyyyy_1,  \
                             g_0_yyyyzzzz_0_xxxxxx_0, \
                             g_0_yyyyzzzz_0_xxxxxy_0, \
                             g_0_yyyyzzzz_0_xxxxxz_0, \
                             g_0_yyyyzzzz_0_xxxxyy_0, \
                             g_0_yyyyzzzz_0_xxxxyz_0, \
                             g_0_yyyyzzzz_0_xxxxzz_0, \
                             g_0_yyyyzzzz_0_xxxyyy_0, \
                             g_0_yyyyzzzz_0_xxxyyz_0, \
                             g_0_yyyyzzzz_0_xxxyzz_0, \
                             g_0_yyyyzzzz_0_xxxzzz_0, \
                             g_0_yyyyzzzz_0_xxyyyy_0, \
                             g_0_yyyyzzzz_0_xxyyyz_0, \
                             g_0_yyyyzzzz_0_xxyyzz_0, \
                             g_0_yyyyzzzz_0_xxyzzz_0, \
                             g_0_yyyyzzzz_0_xxzzzz_0, \
                             g_0_yyyyzzzz_0_xyyyyy_0, \
                             g_0_yyyyzzzz_0_xyyyyz_0, \
                             g_0_yyyyzzzz_0_xyyyzz_0, \
                             g_0_yyyyzzzz_0_xyyzzz_0, \
                             g_0_yyyyzzzz_0_xyzzzz_0, \
                             g_0_yyyyzzzz_0_xzzzzz_0, \
                             g_0_yyyyzzzz_0_yyyyyy_0, \
                             g_0_yyyyzzzz_0_yyyyyz_0, \
                             g_0_yyyyzzzz_0_yyyyzz_0, \
                             g_0_yyyyzzzz_0_yyyzzz_0, \
                             g_0_yyyyzzzz_0_yyzzzz_0, \
                             g_0_yyyyzzzz_0_yzzzzz_0, \
                             g_0_yyyyzzzz_0_zzzzzz_0, \
                             g_0_yyyzzzz_0_xxxxxx_0,  \
                             g_0_yyyzzzz_0_xxxxxx_1,  \
                             g_0_yyyzzzz_0_xxxxxz_0,  \
                             g_0_yyyzzzz_0_xxxxxz_1,  \
                             g_0_yyyzzzz_0_xxxxyz_0,  \
                             g_0_yyyzzzz_0_xxxxyz_1,  \
                             g_0_yyyzzzz_0_xxxxz_1,   \
                             g_0_yyyzzzz_0_xxxxzz_0,  \
                             g_0_yyyzzzz_0_xxxxzz_1,  \
                             g_0_yyyzzzz_0_xxxyyz_0,  \
                             g_0_yyyzzzz_0_xxxyyz_1,  \
                             g_0_yyyzzzz_0_xxxyz_1,   \
                             g_0_yyyzzzz_0_xxxyzz_0,  \
                             g_0_yyyzzzz_0_xxxyzz_1,  \
                             g_0_yyyzzzz_0_xxxzz_1,   \
                             g_0_yyyzzzz_0_xxxzzz_0,  \
                             g_0_yyyzzzz_0_xxxzzz_1,  \
                             g_0_yyyzzzz_0_xxyyyz_0,  \
                             g_0_yyyzzzz_0_xxyyyz_1,  \
                             g_0_yyyzzzz_0_xxyyz_1,   \
                             g_0_yyyzzzz_0_xxyyzz_0,  \
                             g_0_yyyzzzz_0_xxyyzz_1,  \
                             g_0_yyyzzzz_0_xxyzz_1,   \
                             g_0_yyyzzzz_0_xxyzzz_0,  \
                             g_0_yyyzzzz_0_xxyzzz_1,  \
                             g_0_yyyzzzz_0_xxzzz_1,   \
                             g_0_yyyzzzz_0_xxzzzz_0,  \
                             g_0_yyyzzzz_0_xxzzzz_1,  \
                             g_0_yyyzzzz_0_xyyyyz_0,  \
                             g_0_yyyzzzz_0_xyyyyz_1,  \
                             g_0_yyyzzzz_0_xyyyz_1,   \
                             g_0_yyyzzzz_0_xyyyzz_0,  \
                             g_0_yyyzzzz_0_xyyyzz_1,  \
                             g_0_yyyzzzz_0_xyyzz_1,   \
                             g_0_yyyzzzz_0_xyyzzz_0,  \
                             g_0_yyyzzzz_0_xyyzzz_1,  \
                             g_0_yyyzzzz_0_xyzzz_1,   \
                             g_0_yyyzzzz_0_xyzzzz_0,  \
                             g_0_yyyzzzz_0_xyzzzz_1,  \
                             g_0_yyyzzzz_0_xzzzz_1,   \
                             g_0_yyyzzzz_0_xzzzzz_0,  \
                             g_0_yyyzzzz_0_xzzzzz_1,  \
                             g_0_yyyzzzz_0_yyyyyz_0,  \
                             g_0_yyyzzzz_0_yyyyyz_1,  \
                             g_0_yyyzzzz_0_yyyyz_1,   \
                             g_0_yyyzzzz_0_yyyyzz_0,  \
                             g_0_yyyzzzz_0_yyyyzz_1,  \
                             g_0_yyyzzzz_0_yyyzz_1,   \
                             g_0_yyyzzzz_0_yyyzzz_0,  \
                             g_0_yyyzzzz_0_yyyzzz_1,  \
                             g_0_yyyzzzz_0_yyzzz_1,   \
                             g_0_yyyzzzz_0_yyzzzz_0,  \
                             g_0_yyyzzzz_0_yyzzzz_1,  \
                             g_0_yyyzzzz_0_yzzzz_1,   \
                             g_0_yyyzzzz_0_yzzzzz_0,  \
                             g_0_yyyzzzz_0_yzzzzz_1,  \
                             g_0_yyyzzzz_0_zzzzz_1,   \
                             g_0_yyyzzzz_0_zzzzzz_0,  \
                             g_0_yyyzzzz_0_zzzzzz_1,  \
                             g_0_yyzzzz_0_xxxxxx_0,   \
                             g_0_yyzzzz_0_xxxxxx_1,   \
                             g_0_yyzzzz_0_xxxxxz_0,   \
                             g_0_yyzzzz_0_xxxxxz_1,   \
                             g_0_yyzzzz_0_xxxxyz_0,   \
                             g_0_yyzzzz_0_xxxxyz_1,   \
                             g_0_yyzzzz_0_xxxxzz_0,   \
                             g_0_yyzzzz_0_xxxxzz_1,   \
                             g_0_yyzzzz_0_xxxyyz_0,   \
                             g_0_yyzzzz_0_xxxyyz_1,   \
                             g_0_yyzzzz_0_xxxyzz_0,   \
                             g_0_yyzzzz_0_xxxyzz_1,   \
                             g_0_yyzzzz_0_xxxzzz_0,   \
                             g_0_yyzzzz_0_xxxzzz_1,   \
                             g_0_yyzzzz_0_xxyyyz_0,   \
                             g_0_yyzzzz_0_xxyyyz_1,   \
                             g_0_yyzzzz_0_xxyyzz_0,   \
                             g_0_yyzzzz_0_xxyyzz_1,   \
                             g_0_yyzzzz_0_xxyzzz_0,   \
                             g_0_yyzzzz_0_xxyzzz_1,   \
                             g_0_yyzzzz_0_xxzzzz_0,   \
                             g_0_yyzzzz_0_xxzzzz_1,   \
                             g_0_yyzzzz_0_xyyyyz_0,   \
                             g_0_yyzzzz_0_xyyyyz_1,   \
                             g_0_yyzzzz_0_xyyyzz_0,   \
                             g_0_yyzzzz_0_xyyyzz_1,   \
                             g_0_yyzzzz_0_xyyzzz_0,   \
                             g_0_yyzzzz_0_xyyzzz_1,   \
                             g_0_yyzzzz_0_xyzzzz_0,   \
                             g_0_yyzzzz_0_xyzzzz_1,   \
                             g_0_yyzzzz_0_xzzzzz_0,   \
                             g_0_yyzzzz_0_xzzzzz_1,   \
                             g_0_yyzzzz_0_yyyyyz_0,   \
                             g_0_yyzzzz_0_yyyyyz_1,   \
                             g_0_yyzzzz_0_yyyyzz_0,   \
                             g_0_yyzzzz_0_yyyyzz_1,   \
                             g_0_yyzzzz_0_yyyzzz_0,   \
                             g_0_yyzzzz_0_yyyzzz_1,   \
                             g_0_yyzzzz_0_yyzzzz_0,   \
                             g_0_yyzzzz_0_yyzzzz_1,   \
                             g_0_yyzzzz_0_yzzzzz_0,   \
                             g_0_yyzzzz_0_yzzzzz_1,   \
                             g_0_yyzzzz_0_zzzzzz_0,   \
                             g_0_yyzzzz_0_zzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_xxxxxx_0[i] = 3.0 * g_0_yyzzzz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxxxx_0[i] * pb_y + g_0_yyyzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxxy_0[i] = 3.0 * g_0_yyyyzz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxxxy_0[i] * pb_z + g_0_yyyyzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxxxz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxxxz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxyy_0[i] = 3.0 * g_0_yyyyzz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxxyy_0[i] * pb_z + g_0_yyyyzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxxyz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxyz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxxzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxxyyy_0[i] * pb_z + g_0_yyyyzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyyz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xxyyyy_0[i] * pb_z + g_0_yyyyzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyyz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyzz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_xyyyyy_0[i] * pb_z + g_0_yyyyzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xyyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyyz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyzz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyzzz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyyyy_0[i] = 3.0 * g_0_yyyyzz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyyzzz_0_yyyyyy_0[i] * pb_z + g_0_yyyyzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_yyyyyz_0[i] = 3.0 * g_0_yyzzzz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyyzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyyyz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyyzz_0[i] = 3.0 * g_0_yyzzzz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyyzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyyzz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyzzz_0[i] = 3.0 * g_0_yyzzzz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyyzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyzzz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyzzzz_0[i] = 3.0 * g_0_yyzzzz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyyzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyzzzz_0[i] * pb_y +
                                     g_0_yyyzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_zzzzzz_0[i] = 3.0 * g_0_yyzzzz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_zzzzzz_0[i] * pb_y + g_0_yyyzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1148-1176 components of targeted buffer : SLSI

    auto g_0_yyyzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1148);

    auto g_0_yyyzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1149);

    auto g_0_yyyzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1150);

    auto g_0_yyyzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1151);

    auto g_0_yyyzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1152);

    auto g_0_yyyzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1153);

    auto g_0_yyyzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1154);

    auto g_0_yyyzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1155);

    auto g_0_yyyzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1156);

    auto g_0_yyyzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1157);

    auto g_0_yyyzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1158);

    auto g_0_yyyzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1159);

    auto g_0_yyyzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1160);

    auto g_0_yyyzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1161);

    auto g_0_yyyzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1162);

    auto g_0_yyyzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1163);

    auto g_0_yyyzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1164);

    auto g_0_yyyzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1165);

    auto g_0_yyyzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1166);

    auto g_0_yyyzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1167);

    auto g_0_yyyzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1168);

    auto g_0_yyyzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1169);

    auto g_0_yyyzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1170);

    auto g_0_yyyzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1171);

    auto g_0_yyyzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1172);

    auto g_0_yyyzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1173);

    auto g_0_yyyzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1174);

    auto g_0_yyyzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1175);

#pragma omp simd aligned(g_0_yyyzzz_0_xxxxxy_0,       \
                             g_0_yyyzzz_0_xxxxxy_1,   \
                             g_0_yyyzzz_0_xxxxyy_0,   \
                             g_0_yyyzzz_0_xxxxyy_1,   \
                             g_0_yyyzzz_0_xxxyyy_0,   \
                             g_0_yyyzzz_0_xxxyyy_1,   \
                             g_0_yyyzzz_0_xxyyyy_0,   \
                             g_0_yyyzzz_0_xxyyyy_1,   \
                             g_0_yyyzzz_0_xyyyyy_0,   \
                             g_0_yyyzzz_0_xyyyyy_1,   \
                             g_0_yyyzzz_0_yyyyyy_0,   \
                             g_0_yyyzzz_0_yyyyyy_1,   \
                             g_0_yyyzzzz_0_xxxxxy_0,  \
                             g_0_yyyzzzz_0_xxxxxy_1,  \
                             g_0_yyyzzzz_0_xxxxyy_0,  \
                             g_0_yyyzzzz_0_xxxxyy_1,  \
                             g_0_yyyzzzz_0_xxxyyy_0,  \
                             g_0_yyyzzzz_0_xxxyyy_1,  \
                             g_0_yyyzzzz_0_xxyyyy_0,  \
                             g_0_yyyzzzz_0_xxyyyy_1,  \
                             g_0_yyyzzzz_0_xyyyyy_0,  \
                             g_0_yyyzzzz_0_xyyyyy_1,  \
                             g_0_yyyzzzz_0_yyyyyy_0,  \
                             g_0_yyyzzzz_0_yyyyyy_1,  \
                             g_0_yyyzzzzz_0_xxxxxx_0, \
                             g_0_yyyzzzzz_0_xxxxxy_0, \
                             g_0_yyyzzzzz_0_xxxxxz_0, \
                             g_0_yyyzzzzz_0_xxxxyy_0, \
                             g_0_yyyzzzzz_0_xxxxyz_0, \
                             g_0_yyyzzzzz_0_xxxxzz_0, \
                             g_0_yyyzzzzz_0_xxxyyy_0, \
                             g_0_yyyzzzzz_0_xxxyyz_0, \
                             g_0_yyyzzzzz_0_xxxyzz_0, \
                             g_0_yyyzzzzz_0_xxxzzz_0, \
                             g_0_yyyzzzzz_0_xxyyyy_0, \
                             g_0_yyyzzzzz_0_xxyyyz_0, \
                             g_0_yyyzzzzz_0_xxyyzz_0, \
                             g_0_yyyzzzzz_0_xxyzzz_0, \
                             g_0_yyyzzzzz_0_xxzzzz_0, \
                             g_0_yyyzzzzz_0_xyyyyy_0, \
                             g_0_yyyzzzzz_0_xyyyyz_0, \
                             g_0_yyyzzzzz_0_xyyyzz_0, \
                             g_0_yyyzzzzz_0_xyyzzz_0, \
                             g_0_yyyzzzzz_0_xyzzzz_0, \
                             g_0_yyyzzzzz_0_xzzzzz_0, \
                             g_0_yyyzzzzz_0_yyyyyy_0, \
                             g_0_yyyzzzzz_0_yyyyyz_0, \
                             g_0_yyyzzzzz_0_yyyyzz_0, \
                             g_0_yyyzzzzz_0_yyyzzz_0, \
                             g_0_yyyzzzzz_0_yyzzzz_0, \
                             g_0_yyyzzzzz_0_yzzzzz_0, \
                             g_0_yyyzzzzz_0_zzzzzz_0, \
                             g_0_yyzzzzz_0_xxxxxx_0,  \
                             g_0_yyzzzzz_0_xxxxxx_1,  \
                             g_0_yyzzzzz_0_xxxxxz_0,  \
                             g_0_yyzzzzz_0_xxxxxz_1,  \
                             g_0_yyzzzzz_0_xxxxyz_0,  \
                             g_0_yyzzzzz_0_xxxxyz_1,  \
                             g_0_yyzzzzz_0_xxxxz_1,   \
                             g_0_yyzzzzz_0_xxxxzz_0,  \
                             g_0_yyzzzzz_0_xxxxzz_1,  \
                             g_0_yyzzzzz_0_xxxyyz_0,  \
                             g_0_yyzzzzz_0_xxxyyz_1,  \
                             g_0_yyzzzzz_0_xxxyz_1,   \
                             g_0_yyzzzzz_0_xxxyzz_0,  \
                             g_0_yyzzzzz_0_xxxyzz_1,  \
                             g_0_yyzzzzz_0_xxxzz_1,   \
                             g_0_yyzzzzz_0_xxxzzz_0,  \
                             g_0_yyzzzzz_0_xxxzzz_1,  \
                             g_0_yyzzzzz_0_xxyyyz_0,  \
                             g_0_yyzzzzz_0_xxyyyz_1,  \
                             g_0_yyzzzzz_0_xxyyz_1,   \
                             g_0_yyzzzzz_0_xxyyzz_0,  \
                             g_0_yyzzzzz_0_xxyyzz_1,  \
                             g_0_yyzzzzz_0_xxyzz_1,   \
                             g_0_yyzzzzz_0_xxyzzz_0,  \
                             g_0_yyzzzzz_0_xxyzzz_1,  \
                             g_0_yyzzzzz_0_xxzzz_1,   \
                             g_0_yyzzzzz_0_xxzzzz_0,  \
                             g_0_yyzzzzz_0_xxzzzz_1,  \
                             g_0_yyzzzzz_0_xyyyyz_0,  \
                             g_0_yyzzzzz_0_xyyyyz_1,  \
                             g_0_yyzzzzz_0_xyyyz_1,   \
                             g_0_yyzzzzz_0_xyyyzz_0,  \
                             g_0_yyzzzzz_0_xyyyzz_1,  \
                             g_0_yyzzzzz_0_xyyzz_1,   \
                             g_0_yyzzzzz_0_xyyzzz_0,  \
                             g_0_yyzzzzz_0_xyyzzz_1,  \
                             g_0_yyzzzzz_0_xyzzz_1,   \
                             g_0_yyzzzzz_0_xyzzzz_0,  \
                             g_0_yyzzzzz_0_xyzzzz_1,  \
                             g_0_yyzzzzz_0_xzzzz_1,   \
                             g_0_yyzzzzz_0_xzzzzz_0,  \
                             g_0_yyzzzzz_0_xzzzzz_1,  \
                             g_0_yyzzzzz_0_yyyyyz_0,  \
                             g_0_yyzzzzz_0_yyyyyz_1,  \
                             g_0_yyzzzzz_0_yyyyz_1,   \
                             g_0_yyzzzzz_0_yyyyzz_0,  \
                             g_0_yyzzzzz_0_yyyyzz_1,  \
                             g_0_yyzzzzz_0_yyyzz_1,   \
                             g_0_yyzzzzz_0_yyyzzz_0,  \
                             g_0_yyzzzzz_0_yyyzzz_1,  \
                             g_0_yyzzzzz_0_yyzzz_1,   \
                             g_0_yyzzzzz_0_yyzzzz_0,  \
                             g_0_yyzzzzz_0_yyzzzz_1,  \
                             g_0_yyzzzzz_0_yzzzz_1,   \
                             g_0_yyzzzzz_0_yzzzzz_0,  \
                             g_0_yyzzzzz_0_yzzzzz_1,  \
                             g_0_yyzzzzz_0_zzzzz_1,   \
                             g_0_yyzzzzz_0_zzzzzz_0,  \
                             g_0_yyzzzzz_0_zzzzzz_1,  \
                             g_0_yzzzzz_0_xxxxxx_0,   \
                             g_0_yzzzzz_0_xxxxxx_1,   \
                             g_0_yzzzzz_0_xxxxxz_0,   \
                             g_0_yzzzzz_0_xxxxxz_1,   \
                             g_0_yzzzzz_0_xxxxyz_0,   \
                             g_0_yzzzzz_0_xxxxyz_1,   \
                             g_0_yzzzzz_0_xxxxzz_0,   \
                             g_0_yzzzzz_0_xxxxzz_1,   \
                             g_0_yzzzzz_0_xxxyyz_0,   \
                             g_0_yzzzzz_0_xxxyyz_1,   \
                             g_0_yzzzzz_0_xxxyzz_0,   \
                             g_0_yzzzzz_0_xxxyzz_1,   \
                             g_0_yzzzzz_0_xxxzzz_0,   \
                             g_0_yzzzzz_0_xxxzzz_1,   \
                             g_0_yzzzzz_0_xxyyyz_0,   \
                             g_0_yzzzzz_0_xxyyyz_1,   \
                             g_0_yzzzzz_0_xxyyzz_0,   \
                             g_0_yzzzzz_0_xxyyzz_1,   \
                             g_0_yzzzzz_0_xxyzzz_0,   \
                             g_0_yzzzzz_0_xxyzzz_1,   \
                             g_0_yzzzzz_0_xxzzzz_0,   \
                             g_0_yzzzzz_0_xxzzzz_1,   \
                             g_0_yzzzzz_0_xyyyyz_0,   \
                             g_0_yzzzzz_0_xyyyyz_1,   \
                             g_0_yzzzzz_0_xyyyzz_0,   \
                             g_0_yzzzzz_0_xyyyzz_1,   \
                             g_0_yzzzzz_0_xyyzzz_0,   \
                             g_0_yzzzzz_0_xyyzzz_1,   \
                             g_0_yzzzzz_0_xyzzzz_0,   \
                             g_0_yzzzzz_0_xyzzzz_1,   \
                             g_0_yzzzzz_0_xzzzzz_0,   \
                             g_0_yzzzzz_0_xzzzzz_1,   \
                             g_0_yzzzzz_0_yyyyyz_0,   \
                             g_0_yzzzzz_0_yyyyyz_1,   \
                             g_0_yzzzzz_0_yyyyzz_0,   \
                             g_0_yzzzzz_0_yyyyzz_1,   \
                             g_0_yzzzzz_0_yyyzzz_0,   \
                             g_0_yzzzzz_0_yyyzzz_1,   \
                             g_0_yzzzzz_0_yyzzzz_0,   \
                             g_0_yzzzzz_0_yyzzzz_1,   \
                             g_0_yzzzzz_0_yzzzzz_0,   \
                             g_0_yzzzzz_0_yzzzzz_1,   \
                             g_0_yzzzzz_0_zzzzzz_0,   \
                             g_0_yzzzzz_0_zzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_xxxxxx_0[i] = 2.0 * g_0_yzzzzz_0_xxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxxxx_0[i] * pb_y + g_0_yyzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxxy_0[i] = 4.0 * g_0_yyyzzz_0_xxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxxxy_0[i] * pb_z + g_0_yyyzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxxxz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxxxz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxxyy_0[i] * pb_z + g_0_yyyzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxxyz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxyz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxxzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxxyyy_0[i] * pb_z + g_0_yyyzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxyyz_0[i] = 2.0 * g_0_yzzzzz_0_xxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyyz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxyzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xxyyyy_0[i] * pb_z + g_0_yyyzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyyz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyzz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyyyy_0[i] = 4.0 * g_0_yyyzzz_0_xyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_xyyyyy_0[i] * pb_z + g_0_yyyzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xyyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyyz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyzz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyzzz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyyyy_0[i] = 4.0 * g_0_yyyzzz_0_yyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyyzzzz_0_yyyyyy_0[i] * pb_z + g_0_yyyzzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_yyyyyz_0[i] = 2.0 * g_0_yzzzzz_0_yyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yyzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyyyz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyyzz_0[i] = 2.0 * g_0_yzzzzz_0_yyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yyzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyyzz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyzzz_0[i] = 2.0 * g_0_yzzzzz_0_yyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yyzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyzzz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyzzzz_0[i] = 2.0 * g_0_yzzzzz_0_yyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yyzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyzzzz_0[i] * pb_y +
                                     g_0_yyzzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_yzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_zzzzzz_0[i] = 2.0 * g_0_yzzzzz_0_zzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_zzzzzz_0[i] * pb_y + g_0_yyzzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1176-1204 components of targeted buffer : SLSI

    auto g_0_yyzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1176);

    auto g_0_yyzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1177);

    auto g_0_yyzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1178);

    auto g_0_yyzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1179);

    auto g_0_yyzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1180);

    auto g_0_yyzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1181);

    auto g_0_yyzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1182);

    auto g_0_yyzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1183);

    auto g_0_yyzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1184);

    auto g_0_yyzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1185);

    auto g_0_yyzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1186);

    auto g_0_yyzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1187);

    auto g_0_yyzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1188);

    auto g_0_yyzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1189);

    auto g_0_yyzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1190);

    auto g_0_yyzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1191);

    auto g_0_yyzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1192);

    auto g_0_yyzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1193);

    auto g_0_yyzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1194);

    auto g_0_yyzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1195);

    auto g_0_yyzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1196);

    auto g_0_yyzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1197);

    auto g_0_yyzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1198);

    auto g_0_yyzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1199);

    auto g_0_yyzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1200);

    auto g_0_yyzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1201);

    auto g_0_yyzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1202);

    auto g_0_yyzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1203);

#pragma omp simd aligned(g_0_yyzzzz_0_xxxxxy_0,       \
                             g_0_yyzzzz_0_xxxxxy_1,   \
                             g_0_yyzzzz_0_xxxxyy_0,   \
                             g_0_yyzzzz_0_xxxxyy_1,   \
                             g_0_yyzzzz_0_xxxyyy_0,   \
                             g_0_yyzzzz_0_xxxyyy_1,   \
                             g_0_yyzzzz_0_xxyyyy_0,   \
                             g_0_yyzzzz_0_xxyyyy_1,   \
                             g_0_yyzzzz_0_xyyyyy_0,   \
                             g_0_yyzzzz_0_xyyyyy_1,   \
                             g_0_yyzzzz_0_yyyyyy_0,   \
                             g_0_yyzzzz_0_yyyyyy_1,   \
                             g_0_yyzzzzz_0_xxxxxy_0,  \
                             g_0_yyzzzzz_0_xxxxxy_1,  \
                             g_0_yyzzzzz_0_xxxxyy_0,  \
                             g_0_yyzzzzz_0_xxxxyy_1,  \
                             g_0_yyzzzzz_0_xxxyyy_0,  \
                             g_0_yyzzzzz_0_xxxyyy_1,  \
                             g_0_yyzzzzz_0_xxyyyy_0,  \
                             g_0_yyzzzzz_0_xxyyyy_1,  \
                             g_0_yyzzzzz_0_xyyyyy_0,  \
                             g_0_yyzzzzz_0_xyyyyy_1,  \
                             g_0_yyzzzzz_0_yyyyyy_0,  \
                             g_0_yyzzzzz_0_yyyyyy_1,  \
                             g_0_yyzzzzzz_0_xxxxxx_0, \
                             g_0_yyzzzzzz_0_xxxxxy_0, \
                             g_0_yyzzzzzz_0_xxxxxz_0, \
                             g_0_yyzzzzzz_0_xxxxyy_0, \
                             g_0_yyzzzzzz_0_xxxxyz_0, \
                             g_0_yyzzzzzz_0_xxxxzz_0, \
                             g_0_yyzzzzzz_0_xxxyyy_0, \
                             g_0_yyzzzzzz_0_xxxyyz_0, \
                             g_0_yyzzzzzz_0_xxxyzz_0, \
                             g_0_yyzzzzzz_0_xxxzzz_0, \
                             g_0_yyzzzzzz_0_xxyyyy_0, \
                             g_0_yyzzzzzz_0_xxyyyz_0, \
                             g_0_yyzzzzzz_0_xxyyzz_0, \
                             g_0_yyzzzzzz_0_xxyzzz_0, \
                             g_0_yyzzzzzz_0_xxzzzz_0, \
                             g_0_yyzzzzzz_0_xyyyyy_0, \
                             g_0_yyzzzzzz_0_xyyyyz_0, \
                             g_0_yyzzzzzz_0_xyyyzz_0, \
                             g_0_yyzzzzzz_0_xyyzzz_0, \
                             g_0_yyzzzzzz_0_xyzzzz_0, \
                             g_0_yyzzzzzz_0_xzzzzz_0, \
                             g_0_yyzzzzzz_0_yyyyyy_0, \
                             g_0_yyzzzzzz_0_yyyyyz_0, \
                             g_0_yyzzzzzz_0_yyyyzz_0, \
                             g_0_yyzzzzzz_0_yyyzzz_0, \
                             g_0_yyzzzzzz_0_yyzzzz_0, \
                             g_0_yyzzzzzz_0_yzzzzz_0, \
                             g_0_yyzzzzzz_0_zzzzzz_0, \
                             g_0_yzzzzzz_0_xxxxxx_0,  \
                             g_0_yzzzzzz_0_xxxxxx_1,  \
                             g_0_yzzzzzz_0_xxxxxz_0,  \
                             g_0_yzzzzzz_0_xxxxxz_1,  \
                             g_0_yzzzzzz_0_xxxxyz_0,  \
                             g_0_yzzzzzz_0_xxxxyz_1,  \
                             g_0_yzzzzzz_0_xxxxz_1,   \
                             g_0_yzzzzzz_0_xxxxzz_0,  \
                             g_0_yzzzzzz_0_xxxxzz_1,  \
                             g_0_yzzzzzz_0_xxxyyz_0,  \
                             g_0_yzzzzzz_0_xxxyyz_1,  \
                             g_0_yzzzzzz_0_xxxyz_1,   \
                             g_0_yzzzzzz_0_xxxyzz_0,  \
                             g_0_yzzzzzz_0_xxxyzz_1,  \
                             g_0_yzzzzzz_0_xxxzz_1,   \
                             g_0_yzzzzzz_0_xxxzzz_0,  \
                             g_0_yzzzzzz_0_xxxzzz_1,  \
                             g_0_yzzzzzz_0_xxyyyz_0,  \
                             g_0_yzzzzzz_0_xxyyyz_1,  \
                             g_0_yzzzzzz_0_xxyyz_1,   \
                             g_0_yzzzzzz_0_xxyyzz_0,  \
                             g_0_yzzzzzz_0_xxyyzz_1,  \
                             g_0_yzzzzzz_0_xxyzz_1,   \
                             g_0_yzzzzzz_0_xxyzzz_0,  \
                             g_0_yzzzzzz_0_xxyzzz_1,  \
                             g_0_yzzzzzz_0_xxzzz_1,   \
                             g_0_yzzzzzz_0_xxzzzz_0,  \
                             g_0_yzzzzzz_0_xxzzzz_1,  \
                             g_0_yzzzzzz_0_xyyyyz_0,  \
                             g_0_yzzzzzz_0_xyyyyz_1,  \
                             g_0_yzzzzzz_0_xyyyz_1,   \
                             g_0_yzzzzzz_0_xyyyzz_0,  \
                             g_0_yzzzzzz_0_xyyyzz_1,  \
                             g_0_yzzzzzz_0_xyyzz_1,   \
                             g_0_yzzzzzz_0_xyyzzz_0,  \
                             g_0_yzzzzzz_0_xyyzzz_1,  \
                             g_0_yzzzzzz_0_xyzzz_1,   \
                             g_0_yzzzzzz_0_xyzzzz_0,  \
                             g_0_yzzzzzz_0_xyzzzz_1,  \
                             g_0_yzzzzzz_0_xzzzz_1,   \
                             g_0_yzzzzzz_0_xzzzzz_0,  \
                             g_0_yzzzzzz_0_xzzzzz_1,  \
                             g_0_yzzzzzz_0_yyyyyz_0,  \
                             g_0_yzzzzzz_0_yyyyyz_1,  \
                             g_0_yzzzzzz_0_yyyyz_1,   \
                             g_0_yzzzzzz_0_yyyyzz_0,  \
                             g_0_yzzzzzz_0_yyyyzz_1,  \
                             g_0_yzzzzzz_0_yyyzz_1,   \
                             g_0_yzzzzzz_0_yyyzzz_0,  \
                             g_0_yzzzzzz_0_yyyzzz_1,  \
                             g_0_yzzzzzz_0_yyzzz_1,   \
                             g_0_yzzzzzz_0_yyzzzz_0,  \
                             g_0_yzzzzzz_0_yyzzzz_1,  \
                             g_0_yzzzzzz_0_yzzzz_1,   \
                             g_0_yzzzzzz_0_yzzzzz_0,  \
                             g_0_yzzzzzz_0_yzzzzz_1,  \
                             g_0_yzzzzzz_0_zzzzz_1,   \
                             g_0_yzzzzzz_0_zzzzzz_0,  \
                             g_0_yzzzzzz_0_zzzzzz_1,  \
                             g_0_zzzzzz_0_xxxxxx_0,   \
                             g_0_zzzzzz_0_xxxxxx_1,   \
                             g_0_zzzzzz_0_xxxxxz_0,   \
                             g_0_zzzzzz_0_xxxxxz_1,   \
                             g_0_zzzzzz_0_xxxxyz_0,   \
                             g_0_zzzzzz_0_xxxxyz_1,   \
                             g_0_zzzzzz_0_xxxxzz_0,   \
                             g_0_zzzzzz_0_xxxxzz_1,   \
                             g_0_zzzzzz_0_xxxyyz_0,   \
                             g_0_zzzzzz_0_xxxyyz_1,   \
                             g_0_zzzzzz_0_xxxyzz_0,   \
                             g_0_zzzzzz_0_xxxyzz_1,   \
                             g_0_zzzzzz_0_xxxzzz_0,   \
                             g_0_zzzzzz_0_xxxzzz_1,   \
                             g_0_zzzzzz_0_xxyyyz_0,   \
                             g_0_zzzzzz_0_xxyyyz_1,   \
                             g_0_zzzzzz_0_xxyyzz_0,   \
                             g_0_zzzzzz_0_xxyyzz_1,   \
                             g_0_zzzzzz_0_xxyzzz_0,   \
                             g_0_zzzzzz_0_xxyzzz_1,   \
                             g_0_zzzzzz_0_xxzzzz_0,   \
                             g_0_zzzzzz_0_xxzzzz_1,   \
                             g_0_zzzzzz_0_xyyyyz_0,   \
                             g_0_zzzzzz_0_xyyyyz_1,   \
                             g_0_zzzzzz_0_xyyyzz_0,   \
                             g_0_zzzzzz_0_xyyyzz_1,   \
                             g_0_zzzzzz_0_xyyzzz_0,   \
                             g_0_zzzzzz_0_xyyzzz_1,   \
                             g_0_zzzzzz_0_xyzzzz_0,   \
                             g_0_zzzzzz_0_xyzzzz_1,   \
                             g_0_zzzzzz_0_xzzzzz_0,   \
                             g_0_zzzzzz_0_xzzzzz_1,   \
                             g_0_zzzzzz_0_yyyyyz_0,   \
                             g_0_zzzzzz_0_yyyyyz_1,   \
                             g_0_zzzzzz_0_yyyyzz_0,   \
                             g_0_zzzzzz_0_yyyyzz_1,   \
                             g_0_zzzzzz_0_yyyzzz_0,   \
                             g_0_zzzzzz_0_yyyzzz_1,   \
                             g_0_zzzzzz_0_yyzzzz_0,   \
                             g_0_zzzzzz_0_yyzzzz_1,   \
                             g_0_zzzzzz_0_yzzzzz_0,   \
                             g_0_zzzzzz_0_yzzzzz_1,   \
                             g_0_zzzzzz_0_zzzzzz_0,   \
                             g_0_zzzzzz_0_zzzzzz_1,   \
                             wp_y,                    \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_xxxxxx_0[i] = g_0_zzzzzz_0_xxxxxx_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxxx_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxxy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxxxy_0[i] * pb_z + g_0_yyzzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxxxz_0[i] = g_0_zzzzzz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxxz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxxyy_0[i] * pb_z + g_0_yyzzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxxyz_0[i] = g_0_zzzzzz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxyz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxz_1[i] * fi_abcd_0 +
                                     g_0_yzzzzzz_0_xxxxyz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxzz_0[i] = g_0_zzzzzz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxyyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxxyyy_0[i] * pb_z + g_0_yyzzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxyyz_0[i] = g_0_zzzzzz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyyz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxyzz_0[i] = g_0_zzzzzz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzzzz_0_xxxyzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxzzz_0[i] = g_0_zzzzzz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyyyy_0[i] = 5.0 * g_0_yyzzzz_0_xxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xxyyyy_0[i] * pb_z + g_0_yyzzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxyyyz_0[i] = g_0_zzzzzz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyyz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyyzz_0[i] = g_0_zzzzzz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyzzz_0[i] = g_0_zzzzzz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzzzz_0_xxyzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxzzzz_0[i] = g_0_zzzzzz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyyyy_0[i] = 5.0 * g_0_yyzzzz_0_xyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_xyyyyy_0[i] * pb_z + g_0_yyzzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xyyyyz_0[i] = g_0_zzzzzz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yzzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyyz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyyzz_0[i] = g_0_zzzzzz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyzzz_0[i] = g_0_zzzzzz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyzzzz_0[i] = g_0_zzzzzz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzzzz_0_xyzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xzzzzz_0[i] = g_0_zzzzzz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyyyy_0[i] = 5.0 * g_0_yyzzzz_0_yyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_yyzzzzz_0_yyyyyy_0[i] * pb_z + g_0_yyzzzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_yyyyyz_0[i] = g_0_zzzzzz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     5.0 * g_0_yzzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyyyz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyyzz_0[i] = g_0_zzzzzz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_yzzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyyzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyzzz_0[i] = g_0_zzzzzz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_yzzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyzzzz_0[i] = g_0_zzzzzz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_yzzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyzzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yzzzzz_0[i] = g_0_zzzzzz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzzz_1[i] * fi_abcd_0 +
                                     g_0_yzzzzzz_0_yzzzzz_0[i] * pb_y + g_0_yzzzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_zzzzzz_0[i] = g_0_zzzzzz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzzzz_0[i] * pb_y +
                                     g_0_yzzzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1204-1232 components of targeted buffer : SLSI

    auto g_0_yzzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1204);

    auto g_0_yzzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1205);

    auto g_0_yzzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1206);

    auto g_0_yzzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1207);

    auto g_0_yzzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1208);

    auto g_0_yzzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1209);

    auto g_0_yzzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1210);

    auto g_0_yzzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1211);

    auto g_0_yzzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1212);

    auto g_0_yzzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1213);

    auto g_0_yzzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1214);

    auto g_0_yzzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1215);

    auto g_0_yzzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1216);

    auto g_0_yzzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1217);

    auto g_0_yzzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1218);

    auto g_0_yzzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1219);

    auto g_0_yzzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1220);

    auto g_0_yzzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1221);

    auto g_0_yzzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1222);

    auto g_0_yzzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1223);

    auto g_0_yzzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1224);

    auto g_0_yzzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1225);

    auto g_0_yzzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1226);

    auto g_0_yzzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1227);

    auto g_0_yzzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1228);

    auto g_0_yzzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1229);

    auto g_0_yzzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1230);

    auto g_0_yzzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1231);

#pragma omp simd aligned(g_0_yzzzzzzz_0_xxxxxx_0,     \
                             g_0_yzzzzzzz_0_xxxxxy_0, \
                             g_0_yzzzzzzz_0_xxxxxz_0, \
                             g_0_yzzzzzzz_0_xxxxyy_0, \
                             g_0_yzzzzzzz_0_xxxxyz_0, \
                             g_0_yzzzzzzz_0_xxxxzz_0, \
                             g_0_yzzzzzzz_0_xxxyyy_0, \
                             g_0_yzzzzzzz_0_xxxyyz_0, \
                             g_0_yzzzzzzz_0_xxxyzz_0, \
                             g_0_yzzzzzzz_0_xxxzzz_0, \
                             g_0_yzzzzzzz_0_xxyyyy_0, \
                             g_0_yzzzzzzz_0_xxyyyz_0, \
                             g_0_yzzzzzzz_0_xxyyzz_0, \
                             g_0_yzzzzzzz_0_xxyzzz_0, \
                             g_0_yzzzzzzz_0_xxzzzz_0, \
                             g_0_yzzzzzzz_0_xyyyyy_0, \
                             g_0_yzzzzzzz_0_xyyyyz_0, \
                             g_0_yzzzzzzz_0_xyyyzz_0, \
                             g_0_yzzzzzzz_0_xyyzzz_0, \
                             g_0_yzzzzzzz_0_xyzzzz_0, \
                             g_0_yzzzzzzz_0_xzzzzz_0, \
                             g_0_yzzzzzzz_0_yyyyyy_0, \
                             g_0_yzzzzzzz_0_yyyyyz_0, \
                             g_0_yzzzzzzz_0_yyyyzz_0, \
                             g_0_yzzzzzzz_0_yyyzzz_0, \
                             g_0_yzzzzzzz_0_yyzzzz_0, \
                             g_0_yzzzzzzz_0_yzzzzz_0, \
                             g_0_yzzzzzzz_0_zzzzzz_0, \
                             g_0_zzzzzzz_0_xxxxx_1,   \
                             g_0_zzzzzzz_0_xxxxxx_0,  \
                             g_0_zzzzzzz_0_xxxxxx_1,  \
                             g_0_zzzzzzz_0_xxxxxy_0,  \
                             g_0_zzzzzzz_0_xxxxxy_1,  \
                             g_0_zzzzzzz_0_xxxxxz_0,  \
                             g_0_zzzzzzz_0_xxxxxz_1,  \
                             g_0_zzzzzzz_0_xxxxy_1,   \
                             g_0_zzzzzzz_0_xxxxyy_0,  \
                             g_0_zzzzzzz_0_xxxxyy_1,  \
                             g_0_zzzzzzz_0_xxxxyz_0,  \
                             g_0_zzzzzzz_0_xxxxyz_1,  \
                             g_0_zzzzzzz_0_xxxxz_1,   \
                             g_0_zzzzzzz_0_xxxxzz_0,  \
                             g_0_zzzzzzz_0_xxxxzz_1,  \
                             g_0_zzzzzzz_0_xxxyy_1,   \
                             g_0_zzzzzzz_0_xxxyyy_0,  \
                             g_0_zzzzzzz_0_xxxyyy_1,  \
                             g_0_zzzzzzz_0_xxxyyz_0,  \
                             g_0_zzzzzzz_0_xxxyyz_1,  \
                             g_0_zzzzzzz_0_xxxyz_1,   \
                             g_0_zzzzzzz_0_xxxyzz_0,  \
                             g_0_zzzzzzz_0_xxxyzz_1,  \
                             g_0_zzzzzzz_0_xxxzz_1,   \
                             g_0_zzzzzzz_0_xxxzzz_0,  \
                             g_0_zzzzzzz_0_xxxzzz_1,  \
                             g_0_zzzzzzz_0_xxyyy_1,   \
                             g_0_zzzzzzz_0_xxyyyy_0,  \
                             g_0_zzzzzzz_0_xxyyyy_1,  \
                             g_0_zzzzzzz_0_xxyyyz_0,  \
                             g_0_zzzzzzz_0_xxyyyz_1,  \
                             g_0_zzzzzzz_0_xxyyz_1,   \
                             g_0_zzzzzzz_0_xxyyzz_0,  \
                             g_0_zzzzzzz_0_xxyyzz_1,  \
                             g_0_zzzzzzz_0_xxyzz_1,   \
                             g_0_zzzzzzz_0_xxyzzz_0,  \
                             g_0_zzzzzzz_0_xxyzzz_1,  \
                             g_0_zzzzzzz_0_xxzzz_1,   \
                             g_0_zzzzzzz_0_xxzzzz_0,  \
                             g_0_zzzzzzz_0_xxzzzz_1,  \
                             g_0_zzzzzzz_0_xyyyy_1,   \
                             g_0_zzzzzzz_0_xyyyyy_0,  \
                             g_0_zzzzzzz_0_xyyyyy_1,  \
                             g_0_zzzzzzz_0_xyyyyz_0,  \
                             g_0_zzzzzzz_0_xyyyyz_1,  \
                             g_0_zzzzzzz_0_xyyyz_1,   \
                             g_0_zzzzzzz_0_xyyyzz_0,  \
                             g_0_zzzzzzz_0_xyyyzz_1,  \
                             g_0_zzzzzzz_0_xyyzz_1,   \
                             g_0_zzzzzzz_0_xyyzzz_0,  \
                             g_0_zzzzzzz_0_xyyzzz_1,  \
                             g_0_zzzzzzz_0_xyzzz_1,   \
                             g_0_zzzzzzz_0_xyzzzz_0,  \
                             g_0_zzzzzzz_0_xyzzzz_1,  \
                             g_0_zzzzzzz_0_xzzzz_1,   \
                             g_0_zzzzzzz_0_xzzzzz_0,  \
                             g_0_zzzzzzz_0_xzzzzz_1,  \
                             g_0_zzzzzzz_0_yyyyy_1,   \
                             g_0_zzzzzzz_0_yyyyyy_0,  \
                             g_0_zzzzzzz_0_yyyyyy_1,  \
                             g_0_zzzzzzz_0_yyyyyz_0,  \
                             g_0_zzzzzzz_0_yyyyyz_1,  \
                             g_0_zzzzzzz_0_yyyyz_1,   \
                             g_0_zzzzzzz_0_yyyyzz_0,  \
                             g_0_zzzzzzz_0_yyyyzz_1,  \
                             g_0_zzzzzzz_0_yyyzz_1,   \
                             g_0_zzzzzzz_0_yyyzzz_0,  \
                             g_0_zzzzzzz_0_yyyzzz_1,  \
                             g_0_zzzzzzz_0_yyzzz_1,   \
                             g_0_zzzzzzz_0_yyzzzz_0,  \
                             g_0_zzzzzzz_0_yyzzzz_1,  \
                             g_0_zzzzzzz_0_yzzzz_1,   \
                             g_0_zzzzzzz_0_yzzzzz_0,  \
                             g_0_zzzzzzz_0_yzzzzz_1,  \
                             g_0_zzzzzzz_0_zzzzz_1,   \
                             g_0_zzzzzzz_0_zzzzzz_0,  \
                             g_0_zzzzzzz_0_zzzzzz_1,  \
                             wp_y,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_xxxxxx_0[i] = g_0_zzzzzzz_0_xxxxxx_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxy_0[i] = g_0_zzzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxy_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxxz_0[i] = g_0_zzzzzzz_0_xxxxxz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxyy_0[i] =
            2.0 * g_0_zzzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyy_0[i] * pb_y + g_0_zzzzzzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxyz_0[i] = g_0_zzzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxzz_0[i] = g_0_zzzzzzz_0_xxxxzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyyy_0[i] =
            3.0 * g_0_zzzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyy_0[i] * pb_y + g_0_zzzzzzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyyz_0[i] =
            2.0 * g_0_zzzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyz_0[i] * pb_y + g_0_zzzzzzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyzz_0[i] = g_0_zzzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxzzz_0[i] = g_0_zzzzzzz_0_xxxzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyyy_0[i] =
            4.0 * g_0_zzzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyy_0[i] * pb_y + g_0_zzzzzzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyyz_0[i] =
            3.0 * g_0_zzzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyz_0[i] * pb_y + g_0_zzzzzzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyzz_0[i] =
            2.0 * g_0_zzzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyzz_0[i] * pb_y + g_0_zzzzzzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyzzz_0[i] = g_0_zzzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxzzzz_0[i] = g_0_zzzzzzz_0_xxzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyyy_0[i] =
            5.0 * g_0_zzzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyy_0[i] * pb_y + g_0_zzzzzzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyyz_0[i] =
            4.0 * g_0_zzzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyz_0[i] * pb_y + g_0_zzzzzzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyzz_0[i] =
            3.0 * g_0_zzzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyzz_0[i] * pb_y + g_0_zzzzzzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyzzz_0[i] =
            2.0 * g_0_zzzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzzz_0[i] * pb_y + g_0_zzzzzzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyzzzz_0[i] = g_0_zzzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xzzzzz_0[i] = g_0_zzzzzzz_0_xzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyyy_0[i] =
            6.0 * g_0_zzzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyy_0[i] * pb_y + g_0_zzzzzzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyyz_0[i] =
            5.0 * g_0_zzzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyz_0[i] * pb_y + g_0_zzzzzzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyzz_0[i] =
            4.0 * g_0_zzzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyzz_0[i] * pb_y + g_0_zzzzzzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyzzz_0[i] =
            3.0 * g_0_zzzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyzzz_0[i] * pb_y + g_0_zzzzzzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyzzzz_0[i] =
            2.0 * g_0_zzzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzzzz_0[i] * pb_y + g_0_zzzzzzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yzzzzz_0[i] = g_0_zzzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_zzzzzz_0[i] = g_0_zzzzzzz_0_zzzzzz_0[i] * pb_y + g_0_zzzzzzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 1232-1260 components of targeted buffer : SLSI

    auto g_0_zzzzzzzz_0_xxxxxx_0 = pbuffer.data(idx_eri_0_slsi + 1232);

    auto g_0_zzzzzzzz_0_xxxxxy_0 = pbuffer.data(idx_eri_0_slsi + 1233);

    auto g_0_zzzzzzzz_0_xxxxxz_0 = pbuffer.data(idx_eri_0_slsi + 1234);

    auto g_0_zzzzzzzz_0_xxxxyy_0 = pbuffer.data(idx_eri_0_slsi + 1235);

    auto g_0_zzzzzzzz_0_xxxxyz_0 = pbuffer.data(idx_eri_0_slsi + 1236);

    auto g_0_zzzzzzzz_0_xxxxzz_0 = pbuffer.data(idx_eri_0_slsi + 1237);

    auto g_0_zzzzzzzz_0_xxxyyy_0 = pbuffer.data(idx_eri_0_slsi + 1238);

    auto g_0_zzzzzzzz_0_xxxyyz_0 = pbuffer.data(idx_eri_0_slsi + 1239);

    auto g_0_zzzzzzzz_0_xxxyzz_0 = pbuffer.data(idx_eri_0_slsi + 1240);

    auto g_0_zzzzzzzz_0_xxxzzz_0 = pbuffer.data(idx_eri_0_slsi + 1241);

    auto g_0_zzzzzzzz_0_xxyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1242);

    auto g_0_zzzzzzzz_0_xxyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1243);

    auto g_0_zzzzzzzz_0_xxyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1244);

    auto g_0_zzzzzzzz_0_xxyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1245);

    auto g_0_zzzzzzzz_0_xxzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1246);

    auto g_0_zzzzzzzz_0_xyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1247);

    auto g_0_zzzzzzzz_0_xyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1248);

    auto g_0_zzzzzzzz_0_xyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1249);

    auto g_0_zzzzzzzz_0_xyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1250);

    auto g_0_zzzzzzzz_0_xyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1251);

    auto g_0_zzzzzzzz_0_xzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1252);

    auto g_0_zzzzzzzz_0_yyyyyy_0 = pbuffer.data(idx_eri_0_slsi + 1253);

    auto g_0_zzzzzzzz_0_yyyyyz_0 = pbuffer.data(idx_eri_0_slsi + 1254);

    auto g_0_zzzzzzzz_0_yyyyzz_0 = pbuffer.data(idx_eri_0_slsi + 1255);

    auto g_0_zzzzzzzz_0_yyyzzz_0 = pbuffer.data(idx_eri_0_slsi + 1256);

    auto g_0_zzzzzzzz_0_yyzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1257);

    auto g_0_zzzzzzzz_0_yzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1258);

    auto g_0_zzzzzzzz_0_zzzzzz_0 = pbuffer.data(idx_eri_0_slsi + 1259);

#pragma omp simd aligned(g_0_zzzzzz_0_xxxxxx_0,       \
                             g_0_zzzzzz_0_xxxxxx_1,   \
                             g_0_zzzzzz_0_xxxxxy_0,   \
                             g_0_zzzzzz_0_xxxxxy_1,   \
                             g_0_zzzzzz_0_xxxxxz_0,   \
                             g_0_zzzzzz_0_xxxxxz_1,   \
                             g_0_zzzzzz_0_xxxxyy_0,   \
                             g_0_zzzzzz_0_xxxxyy_1,   \
                             g_0_zzzzzz_0_xxxxyz_0,   \
                             g_0_zzzzzz_0_xxxxyz_1,   \
                             g_0_zzzzzz_0_xxxxzz_0,   \
                             g_0_zzzzzz_0_xxxxzz_1,   \
                             g_0_zzzzzz_0_xxxyyy_0,   \
                             g_0_zzzzzz_0_xxxyyy_1,   \
                             g_0_zzzzzz_0_xxxyyz_0,   \
                             g_0_zzzzzz_0_xxxyyz_1,   \
                             g_0_zzzzzz_0_xxxyzz_0,   \
                             g_0_zzzzzz_0_xxxyzz_1,   \
                             g_0_zzzzzz_0_xxxzzz_0,   \
                             g_0_zzzzzz_0_xxxzzz_1,   \
                             g_0_zzzzzz_0_xxyyyy_0,   \
                             g_0_zzzzzz_0_xxyyyy_1,   \
                             g_0_zzzzzz_0_xxyyyz_0,   \
                             g_0_zzzzzz_0_xxyyyz_1,   \
                             g_0_zzzzzz_0_xxyyzz_0,   \
                             g_0_zzzzzz_0_xxyyzz_1,   \
                             g_0_zzzzzz_0_xxyzzz_0,   \
                             g_0_zzzzzz_0_xxyzzz_1,   \
                             g_0_zzzzzz_0_xxzzzz_0,   \
                             g_0_zzzzzz_0_xxzzzz_1,   \
                             g_0_zzzzzz_0_xyyyyy_0,   \
                             g_0_zzzzzz_0_xyyyyy_1,   \
                             g_0_zzzzzz_0_xyyyyz_0,   \
                             g_0_zzzzzz_0_xyyyyz_1,   \
                             g_0_zzzzzz_0_xyyyzz_0,   \
                             g_0_zzzzzz_0_xyyyzz_1,   \
                             g_0_zzzzzz_0_xyyzzz_0,   \
                             g_0_zzzzzz_0_xyyzzz_1,   \
                             g_0_zzzzzz_0_xyzzzz_0,   \
                             g_0_zzzzzz_0_xyzzzz_1,   \
                             g_0_zzzzzz_0_xzzzzz_0,   \
                             g_0_zzzzzz_0_xzzzzz_1,   \
                             g_0_zzzzzz_0_yyyyyy_0,   \
                             g_0_zzzzzz_0_yyyyyy_1,   \
                             g_0_zzzzzz_0_yyyyyz_0,   \
                             g_0_zzzzzz_0_yyyyyz_1,   \
                             g_0_zzzzzz_0_yyyyzz_0,   \
                             g_0_zzzzzz_0_yyyyzz_1,   \
                             g_0_zzzzzz_0_yyyzzz_0,   \
                             g_0_zzzzzz_0_yyyzzz_1,   \
                             g_0_zzzzzz_0_yyzzzz_0,   \
                             g_0_zzzzzz_0_yyzzzz_1,   \
                             g_0_zzzzzz_0_yzzzzz_0,   \
                             g_0_zzzzzz_0_yzzzzz_1,   \
                             g_0_zzzzzz_0_zzzzzz_0,   \
                             g_0_zzzzzz_0_zzzzzz_1,   \
                             g_0_zzzzzzz_0_xxxxx_1,   \
                             g_0_zzzzzzz_0_xxxxxx_0,  \
                             g_0_zzzzzzz_0_xxxxxx_1,  \
                             g_0_zzzzzzz_0_xxxxxy_0,  \
                             g_0_zzzzzzz_0_xxxxxy_1,  \
                             g_0_zzzzzzz_0_xxxxxz_0,  \
                             g_0_zzzzzzz_0_xxxxxz_1,  \
                             g_0_zzzzzzz_0_xxxxy_1,   \
                             g_0_zzzzzzz_0_xxxxyy_0,  \
                             g_0_zzzzzzz_0_xxxxyy_1,  \
                             g_0_zzzzzzz_0_xxxxyz_0,  \
                             g_0_zzzzzzz_0_xxxxyz_1,  \
                             g_0_zzzzzzz_0_xxxxz_1,   \
                             g_0_zzzzzzz_0_xxxxzz_0,  \
                             g_0_zzzzzzz_0_xxxxzz_1,  \
                             g_0_zzzzzzz_0_xxxyy_1,   \
                             g_0_zzzzzzz_0_xxxyyy_0,  \
                             g_0_zzzzzzz_0_xxxyyy_1,  \
                             g_0_zzzzzzz_0_xxxyyz_0,  \
                             g_0_zzzzzzz_0_xxxyyz_1,  \
                             g_0_zzzzzzz_0_xxxyz_1,   \
                             g_0_zzzzzzz_0_xxxyzz_0,  \
                             g_0_zzzzzzz_0_xxxyzz_1,  \
                             g_0_zzzzzzz_0_xxxzz_1,   \
                             g_0_zzzzzzz_0_xxxzzz_0,  \
                             g_0_zzzzzzz_0_xxxzzz_1,  \
                             g_0_zzzzzzz_0_xxyyy_1,   \
                             g_0_zzzzzzz_0_xxyyyy_0,  \
                             g_0_zzzzzzz_0_xxyyyy_1,  \
                             g_0_zzzzzzz_0_xxyyyz_0,  \
                             g_0_zzzzzzz_0_xxyyyz_1,  \
                             g_0_zzzzzzz_0_xxyyz_1,   \
                             g_0_zzzzzzz_0_xxyyzz_0,  \
                             g_0_zzzzzzz_0_xxyyzz_1,  \
                             g_0_zzzzzzz_0_xxyzz_1,   \
                             g_0_zzzzzzz_0_xxyzzz_0,  \
                             g_0_zzzzzzz_0_xxyzzz_1,  \
                             g_0_zzzzzzz_0_xxzzz_1,   \
                             g_0_zzzzzzz_0_xxzzzz_0,  \
                             g_0_zzzzzzz_0_xxzzzz_1,  \
                             g_0_zzzzzzz_0_xyyyy_1,   \
                             g_0_zzzzzzz_0_xyyyyy_0,  \
                             g_0_zzzzzzz_0_xyyyyy_1,  \
                             g_0_zzzzzzz_0_xyyyyz_0,  \
                             g_0_zzzzzzz_0_xyyyyz_1,  \
                             g_0_zzzzzzz_0_xyyyz_1,   \
                             g_0_zzzzzzz_0_xyyyzz_0,  \
                             g_0_zzzzzzz_0_xyyyzz_1,  \
                             g_0_zzzzzzz_0_xyyzz_1,   \
                             g_0_zzzzzzz_0_xyyzzz_0,  \
                             g_0_zzzzzzz_0_xyyzzz_1,  \
                             g_0_zzzzzzz_0_xyzzz_1,   \
                             g_0_zzzzzzz_0_xyzzzz_0,  \
                             g_0_zzzzzzz_0_xyzzzz_1,  \
                             g_0_zzzzzzz_0_xzzzz_1,   \
                             g_0_zzzzzzz_0_xzzzzz_0,  \
                             g_0_zzzzzzz_0_xzzzzz_1,  \
                             g_0_zzzzzzz_0_yyyyy_1,   \
                             g_0_zzzzzzz_0_yyyyyy_0,  \
                             g_0_zzzzzzz_0_yyyyyy_1,  \
                             g_0_zzzzzzz_0_yyyyyz_0,  \
                             g_0_zzzzzzz_0_yyyyyz_1,  \
                             g_0_zzzzzzz_0_yyyyz_1,   \
                             g_0_zzzzzzz_0_yyyyzz_0,  \
                             g_0_zzzzzzz_0_yyyyzz_1,  \
                             g_0_zzzzzzz_0_yyyzz_1,   \
                             g_0_zzzzzzz_0_yyyzzz_0,  \
                             g_0_zzzzzzz_0_yyyzzz_1,  \
                             g_0_zzzzzzz_0_yyzzz_1,   \
                             g_0_zzzzzzz_0_yyzzzz_0,  \
                             g_0_zzzzzzz_0_yyzzzz_1,  \
                             g_0_zzzzzzz_0_yzzzz_1,   \
                             g_0_zzzzzzz_0_yzzzzz_0,  \
                             g_0_zzzzzzz_0_yzzzzz_1,  \
                             g_0_zzzzzzz_0_zzzzz_1,   \
                             g_0_zzzzzzz_0_zzzzzz_0,  \
                             g_0_zzzzzzz_0_zzzzzz_1,  \
                             g_0_zzzzzzzz_0_xxxxxx_0, \
                             g_0_zzzzzzzz_0_xxxxxy_0, \
                             g_0_zzzzzzzz_0_xxxxxz_0, \
                             g_0_zzzzzzzz_0_xxxxyy_0, \
                             g_0_zzzzzzzz_0_xxxxyz_0, \
                             g_0_zzzzzzzz_0_xxxxzz_0, \
                             g_0_zzzzzzzz_0_xxxyyy_0, \
                             g_0_zzzzzzzz_0_xxxyyz_0, \
                             g_0_zzzzzzzz_0_xxxyzz_0, \
                             g_0_zzzzzzzz_0_xxxzzz_0, \
                             g_0_zzzzzzzz_0_xxyyyy_0, \
                             g_0_zzzzzzzz_0_xxyyyz_0, \
                             g_0_zzzzzzzz_0_xxyyzz_0, \
                             g_0_zzzzzzzz_0_xxyzzz_0, \
                             g_0_zzzzzzzz_0_xxzzzz_0, \
                             g_0_zzzzzzzz_0_xyyyyy_0, \
                             g_0_zzzzzzzz_0_xyyyyz_0, \
                             g_0_zzzzzzzz_0_xyyyzz_0, \
                             g_0_zzzzzzzz_0_xyyzzz_0, \
                             g_0_zzzzzzzz_0_xyzzzz_0, \
                             g_0_zzzzzzzz_0_xzzzzz_0, \
                             g_0_zzzzzzzz_0_yyyyyy_0, \
                             g_0_zzzzzzzz_0_yyyyyz_0, \
                             g_0_zzzzzzzz_0_yyyyzz_0, \
                             g_0_zzzzzzzz_0_yyyzzz_0, \
                             g_0_zzzzzzzz_0_yyzzzz_0, \
                             g_0_zzzzzzzz_0_yzzzzz_0, \
                             g_0_zzzzzzzz_0_zzzzzz_0, \
                             wp_z,                    \
                             c_exps,                  \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_xxxxxx_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxx_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxx_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxxxx_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxy_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxy_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxxxy_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxxz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxxz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxxz_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxxz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxxz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxyy_0[i] = 7.0 * g_0_zzzzzz_0_xxxxyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxyy_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxxyy_0[i] * pb_z + g_0_zzzzzzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxyz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxyz_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxyz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzzzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xxxxzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyyy_0[i] = 7.0 * g_0_zzzzzz_0_xxxyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxyyy_0[i] * pb_z + g_0_zzzzzzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyyz_0[i] = 7.0 * g_0_zzzzzz_0_xxxyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyyz_0[i] * pb_z + g_0_zzzzzzz_0_xxxyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzzzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xxxyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzzzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xxxzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyyy_0[i] = 7.0 * g_0_zzzzzz_0_xxyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxyyyy_0[i] * pb_z + g_0_zzzzzzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyyz_0[i] = 7.0 * g_0_zzzzzz_0_xxyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyyz_0[i] * pb_z + g_0_zzzzzzz_0_xxyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyzz_0[i] = 7.0 * g_0_zzzzzz_0_xxyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzzzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xxyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzzzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xxyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzzzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xxzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyyy_0[i] = 7.0 * g_0_zzzzzz_0_xyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xyyyyy_0[i] * pb_z + g_0_zzzzzzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyyz_0[i] = 7.0 * g_0_zzzzzz_0_xyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyyz_0[i] * pb_z + g_0_zzzzzzz_0_xyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyzz_0[i] = 7.0 * g_0_zzzzzz_0_xyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzzzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyzzz_0[i] = 7.0 * g_0_zzzzzz_0_xyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzzzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzzzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xzzzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_zzzzzzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_xzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyyy_0[i] = 7.0 * g_0_zzzzzz_0_yyyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyyy_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_yyyyyy_0[i] * pb_z + g_0_zzzzzzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyyz_0[i] = 7.0 * g_0_zzzzzz_0_yyyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyyz_1[i] * fti_ab_0 +
                                     g_0_zzzzzzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyyz_0[i] * pb_z + g_0_zzzzzzz_0_yyyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyzz_0[i] = 7.0 * g_0_zzzzzz_0_yyyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyzz_1[i] * fti_ab_0 +
                                     2.0 * g_0_zzzzzzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_yyyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyzzz_0[i] = 7.0 * g_0_zzzzzz_0_yyyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyzzz_1[i] * fti_ab_0 +
                                     3.0 * g_0_zzzzzzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_yyyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyzzzz_0[i] = 7.0 * g_0_zzzzzz_0_yyzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyzzzz_1[i] * fti_ab_0 +
                                     4.0 * g_0_zzzzzzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_yyzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_yzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yzzzzz_1[i] * fti_ab_0 +
                                     5.0 * g_0_zzzzzzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_yzzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_zzzzzz_0[i] = 7.0 * g_0_zzzzzz_0_zzzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_zzzzzz_1[i] * fti_ab_0 +
                                     6.0 * g_0_zzzzzzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_zzzzzz_0[i] * pb_z +
                                     g_0_zzzzzzz_0_zzzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
